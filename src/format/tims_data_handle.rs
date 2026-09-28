// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Frame-level access to a Bruker timsTOF `.d` directory: the part of the
//! third-party opentims library (`TimsDataHandle`, `TimsFrame`,
//! `OpenSourceTof2MzConverter`, `OpenSourceScan2ImConverterFactory`) that the
//! source `BrukerTimsFile` (`FORMAT/BrukerTimsFile.h`) links against.
//!
//! OpenMS does not implement the TDF frame codec itself: it fetches opentims
//! v1.2.0b4 (C++) and calls it. There is no suitable Rust crate to replace it
//! (`timsrust` 0.4 pins `rusqlite` 0.31, whose `libsqlite3-sys` cannot link
//! beside this crate's 0.40, and uses its own calibration; see
//! `docs/THIRD_PARTY_CRATE_DECISIONS.md`), so the codec is ported here on top of
//! the crates it needs: `rusqlite` for `analysis.tdf` and `ruzstd` (pure Rust)
//! for the zstd-compressed frames of `analysis.tdf_bin`. The module is behind
//! the `bruker-tdf` feature.
//!
//! # The TDF frame format (opentims `TimsFrame::decompress`, `save_to_buffs`)
//!
//! `Frames.TimsId` is the byte offset of the frame's packet in
//! `analysis.tdf_bin`. A packet is a little-endian `u32` packet size (header
//! included), a `u32` scan count, then one zstd frame. Decompressed, it holds
//! `num_scans + 2 * num_peaks` 32-bit words stored as four byte planes (all
//! lowest bytes, then all second bytes, ...). Word `s + 1` (for `s` below the
//! last scan) is twice the peak count of scan `s`; the last scan takes the
//! remaining peaks. From word `num_scans` on, each peak is a TOF delta and an
//! intensity; the TOF index is the running sum of the deltas of its scan,
//! started at `u32::MAX` so that a first delta of 1 is TOF index 0.
//!
//! Intensities are then scaled by `100 / Frames.AccumulationTime` and rounded
//! half up into a `u32` (`intensity * correction + 0.5`, truncated), as the
//! vendored opentims does.
//!
//! See `docs/BRUKER_TIMS_SUPPORT.md` for the member mapping and the evidence.

use crate::format::rational_scan2im_converter::{
    LinearScan2ImConverter, RationalScan2ImConverter, Scan2InvIonMobilityConverter,
    try_create_rational_converter,
};
use crate::{Error, Result};
use rusqlite::types::ValueRef;
use rusqlite::{Connection, OpenFlags};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

/// Default ceiling on the peaks of one frame, [`TdfLimits::max_frame_peaks`].
///
/// A timsTOF frame has at most about 10^6 peaks (the diaPASEF test run: 927
/// scans, at most 1335 peaks per scan); 2^25 peaks is a 256 MiB decompression
/// buffer.
pub const DEFAULT_MAX_FRAME_PEAKS: u64 = 1 << 25;

/// Default ceiling on the frames of one file, [`TdfLimits::max_frames`].
///
/// A one-hour run has about 10^5 frames; the frame Id is a `u32`.
pub const DEFAULT_MAX_FRAMES: u64 = 1 << 24;

/// Resource ceilings of [`TimsDataHandle`]; all native, opentims has none.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TdfLimits {
    /// Largest `Frames.NumPeaks` accepted, checked before a frame's
    /// decompression buffer is allocated. Default [`DEFAULT_MAX_FRAME_PEAKS`].
    pub max_frame_peaks: u64,
    /// Largest number of `Frames` rows accepted. Default
    /// [`DEFAULT_MAX_FRAMES`].
    pub max_frames: u64,
}

impl Default for TdfLimits {
    fn default() -> Self {
        Self {
            max_frame_peaks: DEFAULT_MAX_FRAME_PEAKS,
            max_frames: DEFAULT_MAX_FRAMES,
        }
    }
}

/// One row of the `Frames` table, the metadata opentims keeps per frame
/// (`TimsFrame`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TdfFrame {
    /// `Frames.Id`, 1-based.
    pub id: u32,
    /// `Frames.NumScans`, the number of TIMS scans (ion-mobility steps).
    pub num_scans: u32,
    /// `Frames.NumPeaks`, the number of peaks over all scans.
    pub num_peaks: u32,
    /// `Frames.MsMsType`: 0 for MS1; 8 for DDA-PASEF and 9 for DIA-PASEF MS2.
    pub msms_type: u32,
    /// `100 / Frames.AccumulationTime`, the opentims `intensity_correction`.
    pub intensity_correction: f64,
    /// `Frames.Time`, the retention time in seconds.
    pub time: f64,
    /// `Frames.TimsId`, the byte offset of the frame in `analysis.tdf_bin`.
    pub tims_id: u64,
}

/// The raw peaks of one decoded frame, in storage order (by scan, then by
/// TOF within a scan). The three vectors are parallel.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FramePeaks {
    /// TIMS scan index of each peak, 0-based.
    pub scan_ids: Vec<u32>,
    /// TOF index of each peak.
    pub tofs: Vec<u32>,
    /// Corrected intensity of each peak (see the module documentation).
    pub intensities: Vec<u32>,
}

/// opentims' open-source TOF-to-m/z model, `OpenSourceTof2MzConverter`: linear
/// in `sqrt(m/z)`.
///
/// `sqrt(m/z) = intercept + slope * tof` with `intercept = sqrt(mz_min)` and
/// `slope = (sqrt(mz_max) - sqrt(mz_min)) / tof_max_index`. For data acquired
/// by `Bruker otofControl` the range is widened by 5 on each side first, as
/// opentims does.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tof2MzConverter {
    intercept: f64,
    slope: f64,
}

impl Tof2MzConverter {
    /// The opentims constructor `OpenSourceTof2MzConverter(mz_min, mz_max,
    /// tof_max_index, is_otof_control)`. Nothing is checked, as in opentims; the
    /// factory [`Tof2MzConverter::from_global_metadata`] checks.
    pub fn new(mz_min: f64, mz_max: f64, tof_max_index: u32, is_otof_control: bool) -> Self {
        let (mz_min, mz_max) = if is_otof_control {
            (mz_min - 5.0, mz_max + 5.0)
        } else {
            (mz_min, mz_max)
        };
        let intercept = mz_min.sqrt();
        Self {
            intercept,
            slope: (mz_max.sqrt() - intercept) / f64::from(tof_max_index),
        }
    }

    /// The opentims `OpenSourceTof2MzConverterFactory::produce`: the model from
    /// the `GlobalMetadata` keys `MzAcqRangeLower`, `MzAcqRangeUpper`,
    /// `DigitizerNumSamples` and `AcquisitionSoftware`, read as C `atof` /
    /// `atol` read them.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the lower bound is not positive, the upper
    /// bound is not above it, or the sample count is 0 (opentims:
    /// "invalid calibration metadata"); the SQLite error text as
    /// [`Error::Io`] when the table cannot be read.
    pub fn from_global_metadata(db: &Connection) -> Result<Self> {
        let mut mz_min = 0.0;
        let mut mz_max = 0.0;
        let mut tof_max = 0u32;
        let mut is_otof = false;
        for (key, value) in global_metadata(
            db,
            "SELECT Key, Value FROM GlobalMetadata WHERE Key IN \
             ('MzAcqRangeLower','MzAcqRangeUpper','DigitizerNumSamples','AcquisitionSoftware')",
        )? {
            match key.as_str() {
                "MzAcqRangeLower" => mz_min = c_atof(&value),
                "MzAcqRangeUpper" => mz_max = c_atof(&value),
                // `static_cast<uint32_t>(std::atol(val))`: the low 32 bits.
                "DigitizerNumSamples" => tof_max = c_atol(&value) as u32,
                "AcquisitionSoftware" => is_otof = value == "Bruker otofControl",
                _ => {}
            }
        }
        if mz_min <= 0.0 || mz_max <= mz_min || tof_max == 0 {
            return Err(Error::InvalidValue(
                "OpenSourceTof2MzConverterFactory: invalid calibration metadata".into(),
            ));
        }
        Ok(Self::new(mz_min, mz_max, tof_max, is_otof))
    }

    /// `sqrt(m/z)` at TOF index 0.
    pub fn intercept(&self) -> f64 {
        self.intercept
    }

    /// Change of `sqrt(m/z)` per TOF index.
    pub fn slope(&self) -> f64 {
        self.slope
    }

    /// m/z of each TOF index, `(intercept + slope * tof)^2`; the opentims
    /// `convert` for `uint32_t` input.
    pub fn convert(&self, tofs: &[u32]) -> Vec<f64> {
        tofs.iter()
            .map(|&tof| {
                let value = self.intercept + self.slope * f64::from(tof);
                value * value
            })
            .collect()
    }

    /// The nearest TOF index of each m/z, the opentims `inverse_convert`:
    /// `val > 0 ? (uint32_t)(val + 0.5) : 0` with `val = (sqrt(mz) -
    /// intercept) / slope`. NaN and non-positive values give 0; values past
    /// `u32::MAX` saturate, where the source's cast is undefined.
    pub fn inverse_convert(&self, mzs: &[f64]) -> Vec<u32> {
        mzs.iter()
            .map(|&mz| {
                let value = (mz.sqrt() - self.intercept) / self.slope;
                if value > 0.0 { (value + 0.5) as u32 } else { 0 }
            })
            .collect()
    }

    /// `OpenSourceTof2MzConverter (linear-in-sqrt)`, the opentims text.
    pub fn description(&self) -> String {
        "OpenSourceTof2MzConverter (linear-in-sqrt)".into()
    }
}

/// The scan-to-1/K0 model a handle uses.
#[derive(Clone, Debug, PartialEq)]
pub enum ScanConverter {
    /// opentims' `OpenSourceScan2ImConverter` from `GlobalMetadata`.
    Linear(LinearScan2ImConverter),
    /// The source's per-frame `RationalScan2ImConverter`.
    Rational(RationalScan2ImConverter),
}

impl ScanConverter {
    /// The model as the trait object the conversions go through.
    pub fn as_converter(&self) -> &dyn Scan2InvIonMobilityConverter {
        match self {
            Self::Linear(linear) => linear,
            Self::Rational(rational) => rational,
        }
    }
}

/// The opentims `OpenSourceScan2ImConverterFactory::produce`: the linear model
/// from the `GlobalMetadata` keys `OneOverK0AcqRangeLower` /
/// `OneOverK0AcqRangeUpper` and `MAX(Frames.NumScans)`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the lower bound is not positive, the upper
/// bound is not above it, or the scan maximum is 0 (opentims: "invalid
/// calibration metadata"); the SQLite error as [`Error::Io`].
pub fn linear_scan_converter(db: &Connection) -> Result<LinearScan2ImConverter> {
    let mut im_min = 0.0;
    let mut im_max = 0.0;
    for (key, value) in global_metadata(
        db,
        "SELECT Key, Value FROM GlobalMetadata WHERE Key IN \
         ('OneOverK0AcqRangeLower','OneOverK0AcqRangeUpper')",
    )? {
        match key.as_str() {
            "OneOverK0AcqRangeLower" => im_min = c_atof(&value),
            "OneOverK0AcqRangeUpper" => im_max = c_atof(&value),
            _ => {}
        }
    }
    let scan_max: u32 = db
        .query_row("SELECT MAX(NumScans) FROM Frames", [], |row| {
            row.get_ref(0).map(sql_text_u32)
        })
        .map_err(sql_error)?;
    if im_min <= 0.0 || im_max <= im_min || scan_max == 0 {
        return Err(Error::InvalidValue(
            "OpenSourceScan2ImConverterFactory: invalid calibration metadata".into(),
        ));
    }
    Ok(LinearScan2ImConverter::new(im_min, im_max, scan_max))
}

/// An open timsTOF `.d` directory: the frame table, the frame file and the two
/// calibrations. The opentims class `TimsDataHandle`.
///
/// Frames are decoded on demand by [`Self::read_frame`] from
/// `analysis.tdf_bin`; nothing but the frame table is held in memory.
#[derive(Debug)]
pub struct TimsDataHandle {
    dir: PathBuf,
    frames: BTreeMap<u32, TdfFrame>,
    bin: File,
    bin_len: u64,
    tof2mz: Tof2MzConverter,
    scan2im: ScanConverter,
    limits: TdfLimits,
}

impl TimsDataHandle {
    /// Open `dir` with the linear calibrations, as the source's
    /// `openTimsDataHandle` always does first (opentims constructor with the
    /// open-source factories).
    ///
    /// Reads `analysis.tdf` (`Frames`: `Id, NumScans, NumPeaks, MsMsType,
    /// AccumulationTime, Time, TimsId`; `GlobalMetadata`:
    /// `TimsCompressionType`, the m/z and 1/K0 ranges) and opens
    /// `analysis.tdf_bin`.
    ///
    /// # Errors
    ///
    /// - [`Error::Io`] when `analysis.tdf` or `analysis.tdf_bin` cannot be
    ///   opened, or SQLite fails (kind `NotFound` for a missing file);
    /// - [`Error::Unsupported`] when `TimsCompressionType` is not 2 (zstd),
    ///   opentims' "not (yet) supported";
    /// - [`Error::Parse`] for a NULL `Frames` column (opentims: "null value
    ///   in SQL column"), a repeated frame Id, or more rows than
    ///   [`TdfLimits::max_frames`];
    /// - the errors of [`Tof2MzConverter::from_global_metadata`] and
    ///   [`linear_scan_converter`].
    pub fn open(dir: impl AsRef<Path>, limits: TdfLimits) -> Result<Self> {
        let dir = dir.as_ref().to_path_buf();
        let tdf = dir.join("analysis.tdf");
        let bin_path = dir.join("analysis.tdf_bin");
        if !tdf.is_file() {
            return Err(Error::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("{} not found", tdf.display()),
            )));
        }
        let db = open_tdf(&tdf)?;
        let frames = read_frames(&db, &limits)?;
        for (_, value) in global_metadata(
            &db,
            "SELECT Key, Value FROM GlobalMetadata WHERE Key == 'TimsCompressionType'",
        )? {
            // `atoi(row[0]) != 2`.
            if c_atol(&value) as i32 != 2 {
                return Err(Error::Unsupported(format!(
                    "Compression algorithm used in your TDF dataset: {value} is not (yet) supported. \
                     Right now only algorithm 2 (zstd) is supported."
                )));
            }
        }
        let bin = File::open(&bin_path)?;
        let bin_len = bin.metadata()?.len();
        let tof2mz = Tof2MzConverter::from_global_metadata(&db)?;
        let scan2im = ScanConverter::Linear(linear_scan_converter(&db)?);
        Ok(Self {
            dir,
            frames,
            bin,
            bin_len,
            tof2mz,
            scan2im,
            limits,
        })
    }

    /// Replace the scan-to-1/K0 model with the per-frame rational one when the
    /// `TimsCalibration` table allows it ([`try_create_rational_converter`]);
    /// returns whether it did.
    pub fn use_rational_calibration(&mut self) -> bool {
        match try_create_rational_converter(&self.dir) {
            Some(rational) => {
                self.scan2im = ScanConverter::Rational(rational);
                true
            }
            None => false,
        }
    }

    /// The `.d` directory.
    pub fn dir(&self) -> &Path {
        &self.dir
    }

    /// Every frame by Id.
    pub fn frames(&self) -> &BTreeMap<u32, TdfFrame> {
        &self.frames
    }

    /// The frame with Id `id`, the opentims `has_frame` / `get_frame`.
    pub fn frame(&self, id: u32) -> Option<&TdfFrame> {
        self.frames.get(&id)
    }

    /// Smallest frame Id; `u32::MAX` for a file without frames, as opentims'
    /// `min_frame_id` initialisation.
    pub fn min_frame_id(&self) -> u32 {
        self.frames.keys().next().copied().unwrap_or(u32::MAX)
    }

    /// Largest frame Id; 0 for a file without frames, as opentims'
    /// `max_frame_id` initialisation.
    pub fn max_frame_id(&self) -> u32 {
        self.frames.keys().next_back().copied().unwrap_or(0)
    }

    /// The TOF-to-m/z model.
    pub fn tof2mz(&self) -> &Tof2MzConverter {
        &self.tof2mz
    }

    /// The scan-to-1/K0 model.
    pub fn scan2im(&self) -> &ScanConverter {
        &self.scan2im
    }

    /// Decode frame `id`: its scan indices, TOF indices and corrected
    /// intensities, the opentims `save_to_buffs`.
    ///
    /// A frame with no peaks is not read, as in opentims.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] for an unknown frame Id, or a frame over
    ///   [`TdfLimits::max_frame_peaks`] (checked before allocation);
    /// - [`Error::Parse`] for a packet that does not fit the file, whose scan
    ///   count disagrees with `Frames.NumScans`, that does not decompress, that
    ///   decompresses to other than `4 * (num_scans + 2 * num_peaks)` bytes, or
    ///   whose per-scan counts exceed `num_peaks`, and for a frame with peaks
    ///   but no scans. opentims checks none of these but the decompression
    ///   error: its scan-count check is an `assert` compiled out of Release
    ///   builds, a short decompression leaves the previous frame's bytes in
    ///   its buffer, and the other cases write past its buffers;
    /// - [`Error::Io`] when the file cannot be read.
    pub fn read_frame(&self, id: u32) -> Result<FramePeaks> {
        let frame = *self
            .frames
            .get(&id)
            .ok_or_else(|| Error::InvalidValue(format!("frame {id} does not exist")))?;
        let num_peaks = u64::from(frame.num_peaks);
        if num_peaks == 0 {
            return Ok(FramePeaks::default());
        }
        if num_peaks > self.limits.max_frame_peaks {
            return Err(Error::InvalidValue(format!(
                "frame {id} has {num_peaks} peaks, above the ceiling {}",
                self.limits.max_frame_peaks
            )));
        }
        if frame.num_scans == 0 {
            return Err(frame_error(id, "has peaks but no scans"));
        }
        let words = u64::from(frame.num_scans) + 2 * num_peaks;
        let bytes = usize::try_from(4 * words)
            .map_err(|_| frame_error(id, "is too large for this platform"))?;
        let packet = self.read_packet(&frame)?;
        let mut data = vec![0u8; bytes];
        let mut decoder = ruzstd::decoding::FrameDecoder::new();
        let written = decoder.decode_all(&packet, &mut data).map_err(|error| {
            frame_error(
                id,
                &format!(
                    "could not be uncompressed ({error}). File is either corrupted, or in a \
                         (yet) unsupported variant of the format."
                ),
            )
        })?;
        // libzstd, under opentims, verifies a frame's content checksum when the
        // frame carries one; ruzstd leaves the comparison to the caller.
        if decoder
            .get_checksum_from_data()
            .zip(decoder.get_calculated_checksum())
            .is_some_and(|(stored, computed)| stored != computed)
        {
            return Err(frame_error(id, "fails its zstd content checksum"));
        }
        if written != bytes {
            return Err(frame_error(
                id,
                &format!("uncompressed to {written} bytes, expected {bytes}"),
            ));
        }
        decode_words(&frame, &data)
    }

    /// The zstd payload of a frame's packet, after the 8-byte header.
    fn read_packet(&self, frame: &TdfFrame) -> Result<Vec<u8>> {
        let id = frame.id;
        let fits = |length: u64| {
            frame
                .tims_id
                .checked_add(length)
                .is_some_and(|end| end <= self.bin_len)
        };
        if !fits(8) {
            return Err(frame_error(id, "starts past the end of analysis.tdf_bin"));
        }
        let mut file = &self.bin;
        file.seek(SeekFrom::Start(frame.tims_id))?;
        let mut header = [0u8; 8];
        file.read_exact(&mut header)?;
        let packet_size = u32::from_le_bytes([header[0], header[1], header[2], header[3]]);
        let scans = u32::from_le_bytes([header[4], header[5], header[6], header[7]]);
        if scans != frame.num_scans {
            return Err(frame_error(
                id,
                &format!(
                    "packet holds {scans} scans, Frames.NumScans says {}",
                    frame.num_scans
                ),
            ));
        }
        if packet_size < 8 || !fits(u64::from(packet_size)) {
            return Err(frame_error(
                id,
                &format!("packet size {packet_size} does not fit analysis.tdf_bin"),
            ));
        }
        let mut payload = vec![0u8; packet_size as usize - 8];
        file.read_exact(&mut payload)?;
        Ok(payload)
    }
}

/// Split the decompressed byte planes into peaks, the opentims
/// `save_to_buffs` loop.
fn decode_words(frame: &TdfFrame, data: &[u8]) -> Result<FramePeaks> {
    let words = data.len() / 4;
    let planes = [
        &data[..words],
        &data[words..2 * words],
        &data[2 * words..3 * words],
        &data[3 * words..],
    ];
    let word = |index: usize| {
        u32::from_le_bytes([
            planes[0][index],
            planes[1][index],
            planes[2][index],
            planes[3][index],
        ])
    };
    let num_peaks = frame.num_peaks as usize;
    let num_scans = frame.num_scans as usize;
    let mut peaks = FramePeaks {
        scan_ids: Vec::with_capacity(num_peaks),
        tofs: Vec::with_capacity(num_peaks),
        intensities: Vec::with_capacity(num_peaks),
    };
    let mut read_offset = num_scans;
    let push_scan = |peaks: &mut FramePeaks, scan: u32, count: usize, read_offset: &mut usize| {
        let mut tof = u32::MAX;
        for _ in 0..count {
            tof = tof.wrapping_add(word(*read_offset));
            let intensity = word(*read_offset + 1);
            *read_offset += 2;
            peaks.scan_ids.push(scan);
            peaks.tofs.push(tof);
            // `static_cast<double>(i) * intensity_correction + 0.5` stored back
            // into a `uint32_t`: truncation. Out of range saturates here; the
            // source's conversion is undefined there.
            peaks
                .intensities
                .push((f64::from(intensity) * frame.intensity_correction + 0.5) as u32);
        }
    };
    for scan in 0..num_scans - 1 {
        let count = (word(scan + 1) / 2) as usize;
        if peaks.tofs.len() + count > num_peaks {
            return Err(frame_error(
                frame.id,
                "scan peak counts exceed Frames.NumPeaks",
            ));
        }
        push_scan(&mut peaks, scan as u32, count, &mut read_offset);
    }
    let remaining = num_peaks - peaks.tofs.len();
    push_scan(
        &mut peaks,
        (num_scans - 1) as u32,
        remaining,
        &mut read_offset,
    );
    Ok(peaks)
}

fn frame_error(id: u32, message: &str) -> Error {
    Error::Parse {
        line: 0,
        message: format!("TDF frame {id} {message}"),
    }
}

fn open_tdf(path: &Path) -> Result<Connection> {
    Connection::open_with_flags(
        path,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
    )
    .map_err(sql_error)
}

/// An SQLite failure as an I/O error carrying SQLite's message.
pub(crate) fn sql_error(error: rusqlite::Error) -> Error {
    Error::Io(std::io::Error::other(error.to_string()))
}

/// `(Key, Value)` rows of a `GlobalMetadata` query, both as text, skipping
/// rows with a NULL, as the opentims callbacks do.
fn global_metadata(db: &Connection, sql: &str) -> Result<Vec<(String, String)>> {
    let mut statement = db.prepare(sql).map_err(sql_error)?;
    let mut rows = statement.query([]).map_err(sql_error)?;
    let mut out = Vec::new();
    while let Some(row) = rows.next().map_err(sql_error)? {
        let key = row.get_ref(0).map_err(sql_error)?;
        let value = row.get_ref(1).map_err(sql_error)?;
        if let (Some(key), Some(value)) = (sql_text(key), sql_text(value)) {
            out.push((key, value));
        }
    }
    Ok(out)
}

/// A column as SQLite's text conversion gives it to a `sqlite3_exec` callback;
/// `None` for NULL.
fn sql_text(value: ValueRef<'_>) -> Option<String> {
    match value {
        ValueRef::Null => None,
        ValueRef::Integer(i) => Some(i.to_string()),
        ValueRef::Real(r) => Some(sqlite_real_text(r)),
        ValueRef::Text(t) | ValueRef::Blob(t) => Some(String::from_utf8_lossy(t).into_owned()),
    }
}

/// SQLite renders a REAL as `%!.15g`; the digits `atof` then reads are the
/// value rounded to 15 significant digits.
fn sqlite_real_text(value: f64) -> String {
    if !value.is_finite() {
        return if value.is_nan() {
            String::new()
        } else if value > 0.0 {
            "Inf".into()
        } else {
            "-Inf".into()
        };
    }
    let text = format!("{value:.14e}");
    let rounded: f64 = text.parse().unwrap_or(value);
    let plain = format!("{rounded}");
    if plain.contains('.') || plain.contains('e') || plain.contains("inf") {
        plain
    } else {
        format!("{plain}.0")
    }
}

/// `MAX(NumScans)` as `atol` of its text, low 32 bits; 0 for NULL.
fn sql_text_u32(value: ValueRef<'_>) -> u32 {
    sql_text(value).map_or(0, |text| c_atol(&text) as u32)
}

/// C `atof`: the longest leading decimal floating-point number after leading
/// whitespace, 0 when there is none; `inf` and `nan` as `strtod` reads them.
/// Hexadecimal floats, which `strtod` also reads, are read as 0 (the TDF
/// metadata never holds one).
pub(crate) fn c_atof(text: &str) -> f64 {
    let text = text.trim_start_matches([' ', '\t', '\n', '\r', '\x0b', '\x0c']);
    let bytes = text.as_bytes();
    let mut end = 0;
    if end < bytes.len() && (bytes[end] == b'+' || bytes[end] == b'-') {
        end += 1;
    }
    let lower = text[end..].to_ascii_lowercase();
    if lower.starts_with("inf") || lower.starts_with("nan") {
        let sign = if text.starts_with('-') { -1.0 } else { 1.0 };
        return if lower.starts_with("inf") {
            sign * f64::INFINITY
        } else {
            f64::NAN
        };
    }
    let digits_start = end;
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        end += 1;
    }
    if end < bytes.len() && bytes[end] == b'.' {
        end += 1;
        while end < bytes.len() && bytes[end].is_ascii_digit() {
            end += 1;
        }
    }
    let mantissa_digits = text[digits_start..end]
        .chars()
        .filter(char::is_ascii_digit)
        .count();
    if mantissa_digits == 0 {
        return 0.0;
    }
    if end < bytes.len() && (bytes[end] == b'e' || bytes[end] == b'E') {
        let mut exponent_end = end + 1;
        if exponent_end < bytes.len()
            && (bytes[exponent_end] == b'+' || bytes[exponent_end] == b'-')
        {
            exponent_end += 1;
        }
        let exponent_digits = exponent_end;
        while exponent_end < bytes.len() && bytes[exponent_end].is_ascii_digit() {
            exponent_end += 1;
        }
        if exponent_end > exponent_digits {
            end = exponent_end;
        }
    }
    text[..end].parse().unwrap_or(0.0)
}

/// C `atol` on a 64-bit `long`: optional sign and leading decimal digits after
/// whitespace, 0 when there are none; saturating where C is undefined.
pub(crate) fn c_atol(text: &str) -> i64 {
    let text = text.trim_start_matches([' ', '\t', '\n', '\r', '\x0b', '\x0c']);
    let bytes = text.as_bytes();
    let mut end = 0;
    if end < bytes.len() && (bytes[end] == b'+' || bytes[end] == b'-') {
        end += 1;
    }
    let digits = end;
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        end += 1;
    }
    if end == digits {
        return 0;
    }
    text[..end].parse().unwrap_or(if text.starts_with('-') {
        i64::MIN
    } else {
        i64::MAX
    })
}

/// The `Frames` table as opentims reads it: every column as text through
/// `atol` / `atof` / `strtoull`, a NULL refused.
fn read_frames(db: &Connection, limits: &TdfLimits) -> Result<BTreeMap<u32, TdfFrame>> {
    let mut statement = db
        .prepare(
            "SELECT Id, NumScans, NumPeaks, MsMsType, AccumulationTime, Time, TimsId from Frames;",
        )
        .map_err(sql_error)?;
    let mut rows = statement.query([]).map_err(sql_error)?;
    let mut frames = BTreeMap::new();
    while let Some(row) = rows.next().map_err(sql_error)? {
        let mut text = Vec::with_capacity(7);
        for column in 0..7 {
            let value = row.get_ref(column).map_err(sql_error)?;
            text.push(sql_text(value).ok_or_else(|| Error::Parse {
                line: 0,
                message: format!("TimsFrameFromSql: null value in SQL column {column}"),
            })?);
        }
        if frames.len() as u64 >= limits.max_frames {
            return Err(Error::Parse {
                line: 0,
                message: format!("more than {} frames", limits.max_frames),
            });
        }
        let frame = TdfFrame {
            id: c_atol(&text[0]) as u32,
            num_scans: c_atol(&text[1]) as u32,
            num_peaks: c_atol(&text[2]) as u32,
            msms_type: c_atol(&text[3]) as u32,
            intensity_correction: 100.0 / c_atof(&text[4]),
            time: c_atof(&text[5]),
            tims_id: c_strtoull(&text[6]),
        };
        if frames.insert(frame.id, frame).is_some() {
            return Err(Error::Parse {
                line: 0,
                message: format!("repeated frame Id {}", frame.id),
            });
        }
    }
    Ok(frames)
}

/// C `strtoull(text, nullptr, 10)`: leading whitespace, optional sign (a `-`
/// negates modulo 2^64), decimal digits; saturates at `u64::MAX` on overflow.
fn c_strtoull(text: &str) -> u64 {
    let text = text.trim_start_matches([' ', '\t', '\n', '\r', '\x0b', '\x0c']);
    let (negative, digits) = match text.as_bytes().first() {
        Some(b'-') => (true, &text[1..]),
        Some(b'+') => (false, &text[1..]),
        _ => (false, text),
    };
    let end = digits.bytes().take_while(u8::is_ascii_digit).count();
    if end == 0 {
        return 0;
    }
    let value = digits[..end].parse::<u64>().unwrap_or(u64::MAX);
    if negative {
        value.wrapping_neg()
    } else {
        value
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn c_number_parsing_follows_libc() {
        assert_eq!(c_atof("  99.998563"), 99.998563);
        assert_eq!(c_atof("1700.000000"), 1700.0);
        assert_eq!(c_atof("1e3x"), 1000.0);
        assert_eq!(c_atof("1e"), 1.0);
        assert_eq!(c_atof("abc"), 0.0);
        assert_eq!(c_atof("-.5"), -0.5);
        assert_eq!(c_atol("396860"), 396860);
        assert_eq!(c_atol(" -7x"), -7);
        assert_eq!(c_atol("x"), 0);
        assert_eq!(c_strtoull("123"), 123);
        assert_eq!(c_strtoull("-1"), u64::MAX);
        assert_eq!(sqlite_real_text(100.0), "100.0");
        assert_eq!(sqlite_real_text(0.1), "0.1");
    }

    #[test]
    fn tof_model_matches_opentims() {
        let model = Tof2MzConverter::new(100.0, 1600.0, 1000, false);
        assert_eq!(model.intercept(), 10.0);
        assert_eq!(model.slope(), 0.03);
        assert_eq!(model.convert(&[0]), vec![100.0]);
        assert_eq!(model.inverse_convert(&[100.0, 0.0]), vec![0, 0]);
        let otof = Tof2MzConverter::new(105.0, 1595.0, 1000, true);
        assert_eq!(otof, model);
    }
}
