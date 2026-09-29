// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! On-disk cache of spectra and chromatograms, ported from
//! `FORMAT/HANDLERS/CachedMzMLHandler.h` (`OpenMS::Internal`).
//!
//! The cache is a time-efficient binary format: reading the data items is
//! fast and possible in any order once the in-memory index is built. Its
//! layout, written and read here field for field as the source does:
//!
//! ```text
//! i32  file identifier, CACHED_MZML_FILE_IDENTIFIER (8094)
//! per spectrum:     u64 n, u64 arrays, i32 ms_level, f64 rt,
//!                   n x f64 m/z, n x f64 intensity, arrays x ARRAY
//! per chromatogram: u64 n, u64 arrays,
//!                   n x f64 rt, n x f64 intensity, arrays x ARRAY
//! ARRAY:            u64 len, u64 name_len, name_len bytes of name, len x f64
//! u64  number of spectra
//! u64  number of chromatograms
//! ```
//!
//! The source writes host byte order with the host's `Size`; every platform it
//! supports is little-endian with a 64-bit `Size`, which this port writes and
//! reads explicitly, so its files are those of the source on those platforms.
//! A spectrum's float arrays are written before its integer arrays, all as
//! doubles; reading returns them all as float arrays. String arrays and all
//! metadata other than MS level and retention time are not cached: the
//! metadata lives in the companion mzML file of [`crate::format::cached_mzml`].
//! See `docs/CACHED_MZML_SUPPORT.md`.

use crate::concept::progress_logger::{ProgressLogger, ProgressReporter};
use crate::kernel::{
    ChromatogramPeak, DataArray, MSChromatogram, MSExperiment, MSSpectrum, Peak1D,
};
use crate::{Error, Result};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

/// Magic number at the start of every cache file (source
/// `CACHED_MZML_FILE_IDENTIFIER`).
pub const CACHED_MZML_FILE_IDENTIFIER: i32 = 8094;

/// Longest array name the source reads; a longer name is skipped there, and
/// refused here (see [`CachedMzMLHandler::read_spectrum_fast`]).
pub const MAX_ARRAY_NAME_BYTES: u64 = 1023;

const SIZE: u64 = 8;
const DATUM: u64 = 8;
const INT: u64 = 4;

/// One array read from the cache: the values as doubles and, for the extra
/// arrays, their name (the source's `OpenSwath::BinaryDataArray`, whose
/// `data` and `description` these are).
///
/// Mirrors `OPENSWATHALGO/DATAACCESS/DataStructures.h`, which this lane does
/// not port; when that header is ported this type should become a re-export.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BinaryDataArray {
    /// The values.
    pub data: Vec<f64>,
    /// The array name; empty for the m/z (or RT) and intensity arrays.
    pub description: String,
}

/// Reads and writes the cache format (source `Internal::CachedMzMLHandler`).
///
/// The handler holds the two indices [`Self::create_memdump_index`] builds:
/// the byte offset of every spectrum and every chromatogram. The source class
/// also derives from `ProgressLogger`; its three progress-reporting members
/// are the `*_with_progress` variants ([`Self::write_memdump_with_progress`],
/// [`Self::read_memdump_with_progress`],
/// [`Self::create_memdump_index_with_progress`]), which take the logger the
/// source object would be.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CachedMzMLHandler {
    spectra_index: Vec<u64>,
    chrom_index: Vec<u64>,
}

fn read_exact_or<R: Read>(reader: &mut R, buffer: &mut [u8], message: &str) -> Result<()> {
    reader
        .read_exact(buffer)
        .map_err(|error| match error.kind() {
            std::io::ErrorKind::UnexpectedEof => parse(message),
            _ => Error::Io(error),
        })
}

fn read_u64<R: Read>(reader: &mut R, message: &str) -> Result<u64> {
    let mut bytes = [0u8; 8];
    read_exact_or(reader, &mut bytes, message)?;
    Ok(u64::from_le_bytes(bytes))
}

fn read_i32<R: Read>(reader: &mut R, message: &str) -> Result<i32> {
    let mut bytes = [0u8; 4];
    read_exact_or(reader, &mut bytes, message)?;
    Ok(i32::from_le_bytes(bytes))
}

fn read_f64<R: Read>(reader: &mut R, message: &str) -> Result<f64> {
    let mut bytes = [0u8; 8];
    read_exact_or(reader, &mut bytes, message)?;
    Ok(f64::from_le_bytes(bytes))
}

fn parse(message: &str) -> Error {
    Error::Parse {
        line: 0,
        message: message.to_owned(),
    }
}

const WRONG_MAGIC: &str =
    "File might not be a cached mzML file (wrong file magic number). Aborting!";
const BAD_SPECTRUM: &str =
    "filestream in: Read an invalid spectrum length, something is wrong here. Aborting.";
const BAD_CHROMATOGRAM: &str =
    "filestream in: Read an invalid chromatogram length, something is wrong here. Aborting.";

/// Bytes between the current position and the end of `reader`.
fn remaining<R: Seek>(reader: &mut R) -> Result<u64> {
    let position = reader.stream_position()?;
    let end = reader.seek(SeekFrom::End(0))?;
    reader.seek(SeekFrom::Start(position))?;
    Ok(end.saturating_sub(position))
}

/// Reads `count` doubles, after checking that the stream holds them.
fn read_doubles<R: Read + Seek>(reader: &mut R, count: u64, message: &str) -> Result<Vec<f64>> {
    let bytes = count.checked_mul(DATUM).ok_or_else(|| parse(message))?;
    if bytes > remaining(reader)? {
        return Err(parse(message));
    }
    let count = usize::try_from(count).map_err(|_| parse(message))?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| Error::InvalidValue("cached mzML allocation failed".into()))?;
    let mut buffer = [0u8; 8];
    for _ in 0..count {
        read_exact_or(reader, &mut buffer, message)?;
        values.push(f64::from_le_bytes(buffer));
    }
    Ok(values)
}

/// A file count as the source's signed progress value: its `Size` converted
/// implicitly to `SignedSize`, which wraps a count above `i64::MAX` around.
fn source_signed(count: u64) -> i64 {
    count as i64
}

/// The source's `static_cast<int>(size) < 0` test, on the low 32 bits.
fn negative_as_int(size: u64) -> bool {
    (size as u32 as i32) < 0
}

impl CachedMzMLHandler {
    /// A handler with empty indices (source default constructor).
    pub fn new() -> Self {
        Self::default()
    }

    /// Write complete spectra and chromatograms as a dump to `out` (source
    /// `writeMemdump`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when the file cannot be created or written; the
    /// source does not check its stream and writes nothing silently.
    pub fn write_memdump(&self, exp: &MSExperiment, out: impl AsRef<Path>) -> Result<()> {
        Self::write_memdump_reporting(exp, out.as_ref(), ProgressReporter::silent())
    }

    /// [`Self::write_memdump`], reporting progress to `logger` as the source's
    /// `writeMemdump` does (`CachedMzMLHandler.cpp:41-57`).
    ///
    /// The calls are the source's: after the file identifier,
    /// `startProgress(0, spectra + chromatograms, "storing binary data")`;
    /// `setProgress(i)` before spectrum `i`, then `setProgress(i)` before
    /// chromatogram `i`, the index starting again at 0 as in the source; and
    /// `endProgress()` once the counts are written and the file is closed. The
    /// written bytes and every error are those of [`Self::write_memdump`]. A
    /// file that cannot be created fails before the start; the source, whose
    /// stream is unchecked, still makes every call. A write failure after the
    /// start leaves the section open, as the source's exceptions do.
    ///
    /// # Errors
    ///
    /// As [`Self::write_memdump`], plus the errors of the progress calls.
    pub fn write_memdump_with_progress(
        &self,
        exp: &MSExperiment,
        out: impl AsRef<Path>,
        logger: &mut ProgressLogger,
    ) -> Result<()> {
        Self::write_memdump_reporting(exp, out.as_ref(), ProgressReporter::new(Some(logger)))
    }

    fn write_memdump_reporting(
        exp: &MSExperiment,
        out: &Path,
        mut progress: ProgressReporter<'_>,
    ) -> Result<()> {
        let mut writer = BufWriter::new(File::create(out)?);
        Self::write_dump(exp, &mut writer, &mut progress)?;
        writer.flush()?;
        drop(writer);
        progress.end()
    }

    /// Write the dump of `exp` to `writer`.
    ///
    /// A spectrum's float arrays follow its peaks, then its integer arrays,
    /// each converted to double. The source skips the arrays of a spectrum or
    /// chromatogram without peaks although it counts them in the header, so
    /// its reader then takes the next record for the arrays; this port writes
    /// them, which its reader and the source's both read correctly.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when writing fails.
    pub fn write_memdump_to(&self, exp: &MSExperiment, writer: impl Write) -> Result<()> {
        Self::write_dump(exp, writer, &mut ProgressReporter::silent())
    }

    /// The dump without the final `endProgress`, which follows the close.
    fn write_dump(
        exp: &MSExperiment,
        mut writer: impl Write,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<()> {
        writer.write_all(&CACHED_MZML_FILE_IDENTIFIER.to_le_bytes())?;
        progress.start_count(
            exp.spectra.len() + exp.chromatograms.len(),
            "storing binary data",
        )?;
        for (i, spectrum) in exp.spectra.iter().enumerate() {
            progress.set_count(i)?;
            Self::write_spectrum(spectrum, &mut writer)?;
        }
        for (i, chromatogram) in exp.chromatograms.iter().enumerate() {
            progress.set_count(i)?;
            Self::write_chromatogram(chromatogram, &mut writer)?;
        }
        writer.write_all(&(exp.spectra.len() as u64).to_le_bytes())?;
        writer.write_all(&(exp.chromatograms.len() as u64).to_le_bytes())?;
        Ok(())
    }

    /// Write only the metadata of `exp` to the mzML file `out_meta` (source
    /// `writeMetadata`): every spectrum and chromatogram keeps its settings
    /// and loses its peaks and data arrays.
    ///
    /// With `add_cache_meta_value`, each spectrum and chromatogram receives one
    /// shared data-processing entry with the action "File format conversion"
    /// and the meta value `cached_data` = `true`, which marks the file as the
    /// metadata half of a cache.
    ///
    /// # Errors
    ///
    /// The errors of the mzML writer, [`crate::format::mzml::store`].
    #[cfg(feature = "mzml")]
    pub fn write_metadata(
        &self,
        exp: &MSExperiment,
        out_meta: impl AsRef<Path>,
        add_cache_meta_value: bool,
    ) -> Result<()> {
        let mut meta = exp.clone();
        for spectrum in &mut meta.spectra {
            spectrum.clear(false);
        }
        for chromatogram in &mut meta.chromatograms {
            chromatogram.clear(false);
        }
        if add_cache_meta_value {
            let mut processing = crate::metadata::DataProcessing::default();
            processing
                .actions
                .insert(crate::metadata::ProcessingAction::FormatConversion);
            processing.metadata.insert(
                "cached_data".into(),
                crate::metadata::MetaValue::from("true"),
            );
            let processing = std::sync::Arc::new(processing);
            for spectrum in &mut meta.spectra {
                spectrum.data_processing.push(processing.clone());
            }
            for chromatogram in &mut meta.chromatograms {
                chromatogram.data_processing.push(processing.clone());
            }
        }
        crate::format::mzml::store(out_meta, &meta)
    }

    /// The source's second implementation of [`Self::write_metadata`]
    /// (`writeMetadata_x`), which copies the experimental settings and each
    /// spectrum and chromatogram instead of the whole experiment. The written
    /// file is the same.
    ///
    /// # Errors
    ///
    /// As [`Self::write_metadata`].
    #[cfg(feature = "mzml")]
    pub fn write_metadata_x(
        &self,
        exp: &MSExperiment,
        out_meta: impl AsRef<Path>,
        add_cache_meta_value: bool,
    ) -> Result<()> {
        self.write_metadata(exp, out_meta, add_cache_meta_value)
    }

    /// Read all spectra and chromatograms from the dump `filename` (source
    /// `readMemdump`). Each spectrum has its MS level, retention time, peaks
    /// and float arrays; each chromatogram its peaks and float arrays.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when the file cannot be opened (the source's
    /// `FileNotFound`, `FileNotReadable` or `IOException`), and
    /// [`Error::Parse`] when it does not start with the file identifier ("File
    /// might not be a cached mzML file (wrong file magic number). Aborting!")
    /// or a record is truncated or invalid (see [`Self::read_spectrum_fast`]).
    pub fn read_memdump(&self, filename: impl AsRef<Path>) -> Result<MSExperiment> {
        Self::read_memdump_reporting(filename.as_ref(), ProgressReporter::silent())
    }

    /// [`Self::read_memdump`], reporting progress to `logger` as the source's
    /// `readMemdump` does (`CachedMzMLHandler.cpp:97-116`).
    ///
    /// The calls are the source's: once the identifier and the two trailing
    /// counts are read, `startProgress(0, spectra + chromatograms, "reading
    /// binary data")`; `setProgress(i)` before spectrum `i`, then before
    /// chromatogram `i`, the index starting again at 0; and `endProgress()`
    /// after the last record. The counts are the file's own, added and
    /// converted to the signed progress type with the source's wrap-around, so
    /// a forged count reaches the logger as the source passes it. The result
    /// and every error are those of [`Self::read_memdump`]; a refused header
    /// makes no call, and a failure in a record leaves the section open, as in
    /// the source.
    ///
    /// # Errors
    ///
    /// As [`Self::read_memdump`], plus the errors of the progress calls.
    pub fn read_memdump_with_progress(
        &self,
        filename: impl AsRef<Path>,
        logger: &mut ProgressLogger,
    ) -> Result<MSExperiment> {
        Self::read_memdump_reporting(filename.as_ref(), ProgressReporter::new(Some(logger)))
    }

    fn read_memdump_reporting(
        filename: &Path,
        mut progress: ProgressReporter<'_>,
    ) -> Result<MSExperiment> {
        let mut reader = BufReader::new(File::open(filename)?);
        let (spectra, chromatograms) = Self::read_header(&mut reader)?;
        let mut exp = MSExperiment::new();
        progress.start(
            0,
            source_signed(spectra.wrapping_add(chromatograms)),
            "reading binary data",
        )?;
        for i in 0..spectra {
            progress.set(source_signed(i))?;
            let mut spectrum = MSSpectrum::default();
            Self::read_spectrum(&mut spectrum, &mut reader)?;
            exp.spectra.push(spectrum);
        }
        for i in 0..chromatograms {
            progress.set(source_signed(i))?;
            let mut chromatogram = MSChromatogram::default();
            Self::read_chromatogram(&mut chromatogram, &mut reader)?;
            exp.chromatograms.push(chromatogram);
        }
        progress.end()?;
        Ok(exp)
    }

    /// Check the identifier, read the two trailing counts and position the
    /// stream after the identifier.
    fn read_header<R: Read + Seek>(reader: &mut R) -> Result<(u64, u64)> {
        if read_i32(reader, WRONG_MAGIC)? != CACHED_MZML_FILE_IDENTIFIER {
            return Err(parse(WRONG_MAGIC));
        }
        let end = reader.seek(SeekFrom::End(0))?;
        if end < INT + 2 * SIZE {
            return Err(parse(WRONG_MAGIC));
        }
        reader.seek(SeekFrom::Start(end - 2 * SIZE))?;
        let spectra = read_u64(reader, WRONG_MAGIC)?;
        let chromatograms = read_u64(reader, WRONG_MAGIC)?;
        reader.seek(SeekFrom::Start(INT))?;
        Ok((spectra, chromatograms))
    }

    /// Create an index on the location of all the spectra and chromatograms
    /// of the dump `filename` (source `createMemdumpIndex`), replacing the
    /// previous indices.
    ///
    /// # Errors
    ///
    /// As [`Self::read_memdump`]; also [`Error::Parse`] when a record would
    /// extend past the end of the file, where the source records offsets
    /// beyond it. On an error the previous indices are kept.
    pub fn create_memdump_index(&mut self, filename: impl AsRef<Path>) -> Result<()> {
        self.create_memdump_index_reporting(filename.as_ref(), ProgressReporter::silent())
    }

    /// [`Self::create_memdump_index`], reporting progress to `logger` as the
    /// source's `createMemdumpIndex` does (`CachedMzMLHandler.cpp:175-221`).
    ///
    /// The calls are the source's: once the identifier and the two trailing
    /// counts are read, `startProgress(0, spectra + chromatograms, "Creating
    /// index for binary spectra")` (the counts converted as at
    /// [`Self::read_memdump_with_progress`]); `setProgress(i)` before the
    /// offset of spectrum `i`, then before that of chromatogram `i`, the index
    /// starting again at 0; and `endProgress()` after the last record. The
    /// indices and every error are those of [`Self::create_memdump_index`]; a
    /// refused header makes no call, and a failure in a record leaves the
    /// section open and the previous indices in place.
    ///
    /// # Errors
    ///
    /// As [`Self::create_memdump_index`], plus the errors of the progress
    /// calls.
    pub fn create_memdump_index_with_progress(
        &mut self,
        filename: impl AsRef<Path>,
        logger: &mut ProgressLogger,
    ) -> Result<()> {
        self.create_memdump_index_reporting(filename.as_ref(), ProgressReporter::new(Some(logger)))
    }

    fn create_memdump_index_reporting(
        &mut self,
        filename: &Path,
        mut progress: ProgressReporter<'_>,
    ) -> Result<()> {
        let mut reader = BufReader::new(File::open(filename)?);
        let (spectra, chromatograms) = Self::read_header(&mut reader)?;
        let end = reader.seek(SeekFrom::End(0))?;
        reader.seek(SeekFrom::Start(INT))?;
        let mut spectra_index = Vec::new();
        let mut chrom_index = Vec::new();
        progress.start(
            0,
            source_signed(spectra.wrapping_add(chromatograms)),
            "Creating index for binary spectra",
        )?;
        for (count, index, extra, message) in [
            (spectra, &mut spectra_index, DATUM + INT, BAD_SPECTRUM),
            (chromatograms, &mut chrom_index, 0, BAD_CHROMATOGRAM),
        ] {
            for i in 0..count {
                progress.set(source_signed(i))?;
                let start = reader.stream_position()?;
                if start >= end {
                    return Err(parse(message));
                }
                index.try_reserve(1).map_err(|_| {
                    Error::InvalidValue("cached mzML index allocation failed".into())
                })?;
                index.push(start);
                let size = read_u64(&mut reader, message)?;
                let arrays = read_u64(&mut reader, message)?;
                let skip = size
                    .checked_mul(2 * DATUM)
                    .and_then(|n| n.checked_add(extra))
                    .ok_or_else(|| parse(message))?;
                Self::skip(&mut reader, skip, end, message)?;
                for _ in 0..arrays {
                    let len = read_u64(&mut reader, message)?;
                    let name = read_u64(&mut reader, message)?;
                    let skip = len
                        .checked_mul(DATUM)
                        .and_then(|n| n.checked_add(name))
                        .ok_or_else(|| parse(message))?;
                    Self::skip(&mut reader, skip, end, message)?;
                }
            }
        }
        self.spectra_index = spectra_index;
        self.chrom_index = chrom_index;
        progress.end()
    }

    fn skip<R: Seek>(reader: &mut R, bytes: u64, end: u64, message: &str) -> Result<()> {
        let position = reader.stream_position()?;
        let target = position.checked_add(bytes).ok_or_else(|| parse(message))?;
        if target > end {
            return Err(parse(message));
        }
        reader.seek(SeekFrom::Start(target))?;
        Ok(())
    }

    /// The byte offset of every spectrum, from [`Self::create_memdump_index`]
    /// (source `getSpectraIndex`).
    pub fn spectra_index(&self) -> &[u64] {
        &self.spectra_index
    }

    /// The byte offset of every chromatogram, from
    /// [`Self::create_memdump_index`] (source `getChromatogramIndex`).
    pub fn chromatogram_index(&self) -> &[u64] {
        &self.chrom_index
    }

    /// Fast access to the spectrum at the current position of `reader`
    /// (source `readSpectrumFast`), returning its arrays, MS level and
    /// retention time.
    ///
    /// The first two arrays are m/z and intensity, with empty names; the
    /// following ones are the extra arrays with their names.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Parse`] with the source's message "filestream in: Read
    /// an invalid spectrum length, something is wrong here. Aborting." when the
    /// record header cannot be read (the source's stream then leaves the
    /// length at -1) or the length's low 32 bits are negative as an `int`, as
    /// the source tests; and with the same message, as native checks, when the
    /// arrays extend past the end of the stream or a name is longer than
    /// [`MAX_ARRAY_NAME_BYTES`] (the source skips such a name and reuses the
    /// previous one) or not UTF-8. [`Error::Io`] when reading fails otherwise.
    pub fn read_spectrum_fast<R: Read + Seek>(
        reader: &mut R,
    ) -> Result<(Vec<BinaryDataArray>, i32, f64)> {
        let size = read_u64(reader, BAD_SPECTRUM)?;
        let arrays = read_u64(reader, BAD_SPECTRUM)?;
        let ms_level = read_i32(reader, BAD_SPECTRUM)?;
        let rt = read_f64(reader, BAD_SPECTRUM)?;
        if negative_as_int(size) {
            return Err(parse(BAD_SPECTRUM));
        }
        let data = Self::read_data_fast(reader, size, arrays, BAD_SPECTRUM)?;
        Ok((data, ms_level, rt))
    }

    /// Fast access to the chromatogram at the current position of `reader`
    /// (source `readChromatogramFast`): the RT and intensity arrays, then the
    /// extra arrays with their names.
    ///
    /// # Errors
    ///
    /// As [`Self::read_spectrum_fast`], with the message "filestream in: Read
    /// an invalid chromatogram length, something is wrong here. Aborting."
    pub fn read_chromatogram_fast<R: Read + Seek>(reader: &mut R) -> Result<Vec<BinaryDataArray>> {
        let size = read_u64(reader, BAD_CHROMATOGRAM)?;
        let arrays = read_u64(reader, BAD_CHROMATOGRAM)?;
        if negative_as_int(size) {
            return Err(parse(BAD_CHROMATOGRAM));
        }
        Self::read_data_fast(reader, size, arrays, BAD_CHROMATOGRAM)
    }

    /// The source's `readDataFast_`.
    fn read_data_fast<R: Read + Seek>(
        reader: &mut R,
        size: u64,
        arrays: u64,
        message: &str,
    ) -> Result<Vec<BinaryDataArray>> {
        // Each extra array needs at least its two length fields.
        if arrays > remaining(reader)? / (2 * SIZE) {
            return Err(parse(message));
        }
        let mut data = vec![
            BinaryDataArray {
                data: read_doubles(reader, size, message)?,
                description: String::new(),
            },
            BinaryDataArray::default(),
        ];
        data[1].data = read_doubles(reader, size, message)?;
        for _ in 0..arrays {
            let len = read_u64(reader, message)?;
            let name_len = read_u64(reader, message)?;
            if name_len > MAX_ARRAY_NAME_BYTES {
                return Err(parse(message));
            }
            let mut name = vec![0u8; name_len as usize];
            read_exact_or(reader, &mut name, message)?;
            let description = String::from_utf8(name).map_err(|_| parse(message))?;
            data.push(BinaryDataArray {
                data: read_doubles(reader, len, message)?,
                description,
            });
        }
        Ok(data)
    }

    /// Read the spectrum at the current position of `reader` into `spectrum`
    /// (source `readSpectrum`): its MS level and retention time are set, the
    /// peaks are appended (intensities narrowed to `f32`), and one float array
    /// per extra array is appended with its name.
    ///
    /// # Errors
    ///
    /// As [`Self::read_spectrum_fast`]; `spectrum` is unchanged on an error.
    pub fn read_spectrum<R: Read + Seek>(spectrum: &mut MSSpectrum, reader: &mut R) -> Result<()> {
        let (data, ms_level, rt) = Self::read_spectrum_fast(reader)?;
        let level = u32::try_from(ms_level).map_err(|_| parse(BAD_SPECTRUM))?;
        let mut data = data.into_iter();
        let mz = data.next().unwrap_or_default().data;
        let intensity = data.next().unwrap_or_default().data;
        spectrum
            .peaks
            .try_reserve(mz.len())
            .map_err(|_| Error::InvalidValue("cached mzML allocation failed".into()))?;
        spectrum.ms_level = level;
        spectrum.rt = rt;
        spectrum.peaks.extend(
            mz.iter()
                .zip(&intensity)
                .map(|(&mz, &intensity)| Peak1D::new(mz, intensity as f32)),
        );
        for array in data {
            spectrum.float_data_arrays.push(DataArray::new(
                array.description,
                array.data.iter().map(|&value| value as f32).collect(),
            ));
        }
        Ok(())
    }

    /// Read the chromatogram at the current position of `reader` into
    /// `chromatogram` (source `readChromatogram`): the peaks are appended and
    /// the float arrays replaced by one array per extra array.
    ///
    /// The source replaces them with arrays that carry the names but no
    /// values: its copy loop iterates over the new, empty array instead of the
    /// data it read. This port fills in the values, which the file holds.
    ///
    /// # Errors
    ///
    /// As [`Self::read_chromatogram_fast`]; `chromatogram` is unchanged on an
    /// error.
    pub fn read_chromatogram<R: Read + Seek>(
        chromatogram: &mut MSChromatogram,
        reader: &mut R,
    ) -> Result<()> {
        let data = Self::read_chromatogram_fast(reader)?;
        let mut data = data.into_iter();
        let rt = data.next().unwrap_or_default().data;
        let intensity = data.next().unwrap_or_default().data;
        chromatogram
            .peaks
            .try_reserve(rt.len())
            .map_err(|_| Error::InvalidValue("cached mzML allocation failed".into()))?;
        chromatogram.peaks.extend(
            rt.iter()
                .zip(&intensity)
                .map(|(&rt, &intensity)| ChromatogramPeak::new(rt, intensity as f32)),
        );
        chromatogram.float_data_arrays = data
            .map(|array| {
                DataArray::new(
                    array.description,
                    array.data.iter().map(|&value| value as f32).collect(),
                )
            })
            .collect();
        Ok(())
    }

    fn write_arrays(
        floats: &[DataArray<f32>],
        integers: &[DataArray<i32>],
        writer: &mut impl Write,
    ) -> Result<()> {
        for array in floats {
            Self::write_array(
                &array.name,
                array.data.iter().map(|&v| f64::from(v)),
                array.data.len(),
                writer,
            )?;
        }
        for array in integers {
            Self::write_array(
                &array.name,
                array.data.iter().map(|&v| f64::from(v)),
                array.data.len(),
                writer,
            )?;
        }
        Ok(())
    }

    fn write_array(
        name: &str,
        values: impl Iterator<Item = f64>,
        len: usize,
        writer: &mut impl Write,
    ) -> Result<()> {
        writer.write_all(&(len as u64).to_le_bytes())?;
        writer.write_all(&(name.len() as u64).to_le_bytes())?;
        writer.write_all(name.as_bytes())?;
        for value in values {
            writer.write_all(&value.to_le_bytes())?;
        }
        Ok(())
    }

    /// The source's `writeSpectrum_`.
    fn write_spectrum(spectrum: &MSSpectrum, writer: &mut impl Write) -> Result<()> {
        writer.write_all(&(spectrum.peaks.len() as u64).to_le_bytes())?;
        let arrays = spectrum.float_data_arrays.len() + spectrum.integer_data_arrays.len();
        writer.write_all(&(arrays as u64).to_le_bytes())?;
        // The source stores the level in an `int`, truncating as C++ does.
        writer.write_all(&(spectrum.ms_level as i32).to_le_bytes())?;
        writer.write_all(&spectrum.rt.to_le_bytes())?;
        for peak in &spectrum.peaks {
            writer.write_all(&peak.mz.to_le_bytes())?;
        }
        for peak in &spectrum.peaks {
            writer.write_all(&f64::from(peak.intensity).to_le_bytes())?;
        }
        Self::write_arrays(
            &spectrum.float_data_arrays,
            &spectrum.integer_data_arrays,
            writer,
        )
    }

    /// The source's `writeChromatogram_`.
    fn write_chromatogram(chromatogram: &MSChromatogram, writer: &mut impl Write) -> Result<()> {
        writer.write_all(&(chromatogram.peaks.len() as u64).to_le_bytes())?;
        let arrays = chromatogram.float_data_arrays.len() + chromatogram.integer_data_arrays.len();
        writer.write_all(&(arrays as u64).to_le_bytes())?;
        for peak in &chromatogram.peaks {
            writer.write_all(&peak.rt.to_le_bytes())?;
        }
        for peak in &chromatogram.peaks {
            writer.write_all(&f64::from(peak.intensity).to_le_bytes())?;
        }
        Self::write_arrays(
            &chromatogram.float_data_arrays,
            &chromatogram.integer_data_arrays,
            writer,
        )
    }
}
