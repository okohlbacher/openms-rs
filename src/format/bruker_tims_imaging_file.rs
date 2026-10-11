// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reader for Bruker timsTOF MALDI imaging `.d` directories, the source
//! `BrukerTimsImagingFile` (`FORMAT/BrukerTimsImagingFile.h`), behind the
//! `bruker-tdf` and `mzml` features (the latter for the imaging experiment).
//!
//! Pixel coordinates come from the `MaldiFrameInfo` table of `analysis.tdf`
//! and are bound to the per-frame spectra that [`BrukerTimsFile`](crate::format::bruker_tims_file::BrukerTimsFile) produces in
//! `Frame` export mode. No proprietary MALDI SDK is involved: the positions
//! live entirely in the SQLite metadata.
//!
//! Coordinates are normalised to zero-based pixels by subtracting the minimum
//! observed `XIndexPos` and `YIndexPos`. The geometry is strictly 2D: datasets
//! with more than one distinct `ZIndexPos` are refused; load each section into
//! its own experiment.
//!
//! See `docs/BRUKER_TIMS_SUPPORT.md` for the member mapping and the evidence.

use crate::concept::progress_logger::ProgressLogType;
use crate::concept::string_conversion;
use crate::format::bruker_tims_file::{self, BrukerTimsFile, ExportMode, sql_int};
use crate::format::imzml_file::ImagingExperiment;
use crate::format::tims_data_handle::{sql_error, sql_text};
use crate::kernel::on_disc_imzml_experiment::ImagingGeometry;
use crate::{Error, Result};
use rusqlite::{Connection, OpenFlags, OptionalExtension};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// One row of `MaldiFrameInfo`, the source `MaldiPixel`. `ZIndexPos` is
/// intentionally not exposed.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MaldiPixel {
    /// Bruker frame Id, 1-based in `analysis.tdf`.
    pub frame_id: i32,
    /// Raw `XIndexPos`.
    pub x: i32,
    /// Raw `YIndexPos`.
    pub y: i32,
}

/// Processing and export configuration, the source
/// `BrukerTimsImagingFile::Config`.
#[derive(Clone, Debug, PartialEq)]
pub struct Config {
    /// Configuration of the inner [`BrukerTimsFile`]; its export mode is
    /// forced to `Frame` and `load_ms1` to true before delegation.
    pub inner_config: bruker_tims_file::Config,
    /// Refuse datasets whose `MaldiApplicationType` is not `Imaging`.
    /// Default true.
    pub strict_imaging_only: bool,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            inner_config: bruker_tims_file::Config::default(),
            strict_imaging_only: true,
        }
    }
}

/// Reader for timsTOF MALDI imaging data, the source class
/// `BrukerTimsImagingFile` (a `ProgressLogger`).
#[derive(Clone, Debug, Default)]
pub struct BrukerTimsImagingFile {
    log_type: ProgressLogType,
}

impl BrukerTimsImagingFile {
    /// A reader with no progress output, the source default constructor.
    pub fn new() -> Self {
        Self::default()
    }

    /// Progress reporting mode, from the source's `ProgressLogger` base; handed
    /// to the inner [`BrukerTimsFile`].
    pub fn log_type(&self) -> ProgressLogType {
        self.log_type
    }

    /// Set the progress reporting mode.
    pub fn set_log_type(&mut self, log_type: ProgressLogType) {
        self.log_type = log_type;
    }

    /// Load a MALDI imaging `.d` directory, the source `load(path, exp,
    /// config)`; the source overload without a config is
    /// `load(path, &Config::default())`. The result replaces, rather than
    /// extends, any previous experiment.
    ///
    /// Every frame is loaded in `Frame` export mode; each `MaldiFrameInfo` row
    /// binds pixel `(x - x_min, y - y_min)` to the spectrum of its frame. The
    /// grid is the bounding box of the observed coordinates, and the pixel size
    /// is `GlobalMetadata.MaldiSpotSize_um` (else `ImagingAreaSpotSize_um`,
    /// else 1) micrometres; a value that does not parse as a number is 1, as
    /// in the source.
    ///
    /// A pixel whose frame was not loaded (for example because
    /// `inner_config` restricts the frame range) is skipped and only counted
    /// into one warning, as the source's `@note` says: the image may be
    /// incomplete without `load` failing.
    ///
    /// # Errors
    ///
    /// - [`Error::Io`] of kind `NotFound` when the directory or its
    ///   `analysis.tdf` is missing (the source's `FileNotReadable`), checked
    ///   before the imaging-mode test so a missing file is never reported as
    ///   "not imaging";
    /// - [`Error::InvalidValue`] when `strict_imaging_only` is set and the
    ///   dataset is not MALDI imaging, or when `MaldiFrameInfo` holds more than
    ///   one distinct `ZIndexPos`;
    /// - the errors of [`Self::read_maldi_frame_info`] and of
    ///   [`BrukerTimsFile::load`];
    /// - [`Error::InvalidValue`] when the bounding box does not fit `u32`, is
    ///   above the geometry's pixel ceiling, or two rows share a coordinate
    ///   (the source's `addPixel` throws `InvalidValue` for the latter).
    pub fn load(&self, path: impl AsRef<Path>, config: &Config) -> Result<ImagingExperiment> {
        let path = path.as_ref();
        if !path.exists() {
            return Err(not_readable(path));
        }
        let tdf = resolve_tdf(path)?;
        if config.strict_imaging_only && !Self::is_imaging_dataset(path) {
            return Err(Error::InvalidValue(format!(
                "Not a MALDI imaging dataset (MaldiApplicationType != 'Imaging'): {}",
                path.display()
            )));
        }
        {
            let db = open(&tdf)?;
            if column_exists(&db, "MaldiFrameInfo", "ZIndexPos")? {
                let distinct: i32 = db
                    .query_row(
                        "SELECT COUNT(DISTINCT ZIndexPos) FROM MaldiFrameInfo;",
                        [],
                        |row| row.get_ref(0).map(sql_int),
                    )
                    .map_err(sql_error)?;
                if distinct > 1 {
                    return Err(Error::InvalidValue(format!(
                        "Multi-section MALDI imaging (>1 distinct ZIndexPos) is out of scope for the \
                         2D MSImagingGeometry; split into per-section .d folders or filter \
                         MaldiFrameInfo before loading: {}",
                        path.display()
                    )));
                }
            }
        }
        let rows = Self::read_maldi_frame_info(path)?;

        let mut inner_config = config.inner_config.clone();
        inner_config.export_mode = ExportMode::Frame;
        inner_config.load_ms1 = true;
        let mut inner = BrukerTimsFile::new();
        inner.set_log_type(self.log_type);
        let experiment = inner.load(path, &inner_config)?;

        let mut frame_to_index = BTreeMap::new();
        for (index, spectrum) in experiment.spectra.iter().enumerate() {
            if let Some(frame) = parse_frame_id(&spectrum.native_id) {
                frame_to_index.insert(frame, index);
            }
        }

        let (x_min, y_min, width, height) = bounding_box(&rows)?;
        let mut geometry = ImagingGeometry::new();
        geometry.set_dimensions(width, height)?;
        let mut spot = Self::read_global_metadata_value(path, "MaldiSpotSize_um")?;
        if spot.is_empty() {
            spot = Self::read_global_metadata_value(path, "ImagingAreaSpotSize_um")?;
        }
        let spot_size = if spot.is_empty() {
            1.0
        } else {
            to_double(&spot).unwrap_or(1.0)
        };
        geometry.set_pixel_size(spot_size, spot_size, "micrometer")?;

        let mut dropped = 0usize;
        for row in &rows {
            match frame_to_index.get(&row.frame_id) {
                Some(&index) => geometry.add_pixel(
                    (i64::from(row.x) - x_min) as u32,
                    (i64::from(row.y) - y_min) as u32,
                    index,
                )?,
                None => dropped += 1,
            }
        }
        if dropped > 0 {
            let _ = crate::openms_log_warn!(
                "BrukerTimsImagingFile: dropped {dropped} pixel(s) whose Frame ID was not present in the loaded spectra."
            );
        }
        let mut exp = ImagingExperiment::new();
        exp.set_ms_experiment(experiment);
        exp.set_geometry(geometry);
        exp.validate()?;
        Ok(exp)
    }

    /// Whether `GlobalMetadata.MaldiApplicationType` is `Imaging`, the source
    /// `isImagingDataset`: a non-failing probe, false on any I/O error. Use
    /// [`Self::read_global_metadata_value`] to tell "not imaging" from
    /// "unreadable".
    pub fn is_imaging_dataset(path: impl AsRef<Path>) -> bool {
        Self::read_global_metadata_value(path, "MaldiApplicationType")
            .is_ok_and(|value| value == "Imaging")
    }

    /// `MaldiFrameInfo` rows in frame order, raw coordinates, the source
    /// `readMaldiFrameInfo`. Values are read with SQLite's integer conversion
    /// (`sqlite3_column_int`). `ZIndexPos` is not read.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] (kind `NotFound` for a missing `analysis.tdf`) when the
    /// database cannot be opened or read; [`Error::Parse`] when
    /// `MaldiFrameInfo` is missing or empty.
    pub fn read_maldi_frame_info(d_folder: impl AsRef<Path>) -> Result<Vec<MaldiPixel>> {
        let tdf = resolve_tdf(d_folder.as_ref())?;
        let db = open(&tdf)?;
        if !table_exists(&db, "MaldiFrameInfo")? {
            return Err(parse_error(&tdf, "MaldiFrameInfo table missing"));
        }
        let mut statement = db
            .prepare("SELECT Frame, XIndexPos, YIndexPos FROM MaldiFrameInfo ORDER BY Frame;")
            .map_err(sql_error)?;
        let mut rows = statement.query([]).map_err(sql_error)?;
        let mut pixels = Vec::new();
        while let Some(row) = rows.next().map_err(sql_error)? {
            let column = |index: usize| row.get_ref(index).map(sql_int).map_err(sql_error);
            pixels.push(MaldiPixel {
                frame_id: column(0)?,
                x: column(1)?,
                y: column(2)?,
            });
        }
        if pixels.is_empty() {
            return Err(parse_error(&tdf, "MaldiFrameInfo is empty"));
        }
        Ok(pixels)
    }

    /// The value of one `GlobalMetadata` key as text, the source
    /// `readGlobalMetadataValue`; empty when the key or the table is absent
    /// (or the value is NULL).
    ///
    /// # Errors
    ///
    /// [`Error::Io`] (kind `NotFound` for a missing `analysis.tdf`) when the
    /// database cannot be opened or read.
    pub fn read_global_metadata_value(d_folder: impl AsRef<Path>, key: &str) -> Result<String> {
        let tdf = resolve_tdf(d_folder.as_ref())?;
        let db = open(&tdf)?;
        if !table_exists(&db, "GlobalMetadata")? {
            return Ok(String::new());
        }
        let value = db
            .query_row(
                "SELECT Value FROM GlobalMetadata WHERE Key = ?1 LIMIT 1;",
                [key],
                |row| row.get_ref(0).map(sql_text),
            )
            .optional()
            .map_err(sql_error)?;
        Ok(value.flatten().unwrap_or_default())
    }
}

fn not_readable(path: &Path) -> Error {
    Error::Io(std::io::Error::new(
        std::io::ErrorKind::NotFound,
        format!("{} is not readable (does not exist)", path.display()),
    ))
}

fn parse_error(tdf: &Path, message: &str) -> Error {
    Error::Parse {
        line: 0,
        message: format!("{}: {message}", tdf.display()),
    }
}

/// The source's `resolveTdfPath_`.
fn resolve_tdf(d_folder: &Path) -> Result<PathBuf> {
    let tdf = d_folder.join("analysis.tdf");
    if tdf.exists() {
        Ok(tdf)
    } else {
        Err(not_readable(&tdf))
    }
}

fn open(tdf: &Path) -> Result<Connection> {
    Connection::open_with_flags(
        tdf,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
    )
    .map_err(sql_error)
}

fn table_exists(db: &Connection, table: &str) -> Result<bool> {
    db.query_row(
        "SELECT COUNT(*) FROM sqlite_master WHERE type = 'table' AND name = ?1",
        [table],
        |row| row.get::<_, i64>(0),
    )
    .map(|count| count > 0)
    .map_err(sql_error)
}

fn column_exists(db: &Connection, table: &str, column: &str) -> Result<bool> {
    if !table_exists(db, table)? {
        return Ok(false);
    }
    let mut statement = db
        .prepare(&format!("PRAGMA table_info(\"{table}\")"))
        .map_err(sql_error)?;
    let names = statement
        .query_map([], |row| row.get::<_, String>(1))
        .map_err(sql_error)?;
    for name in names {
        if name.map_err(sql_error)? == column {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The source's `parseFrameId`: the integer after a leading `frame=`, up to
/// the first space, through `StringUtils::toInt32` (the shared
/// [`crate::concept::string_conversion::to_int32`]); `None` where the
/// source returns -1.
fn parse_frame_id(native_id: &str) -> Option<i32> {
    let rest = native_id.strip_prefix("frame=")?;
    let field = rest.split(' ').next().unwrap_or("");
    string_conversion::to_int32(field).ok()
}

/// `StringUtils::toDouble` through the shared
/// [`crate::concept::string_conversion::to_double`]; `None` where the
/// source throws (a finite literal that overflows or a nonzero one that
/// underflows to zero is out of range for the Release build's
/// `std::from_chars`).
fn to_double(text: &str) -> Option<f64> {
    string_conversion::to_double(text).ok()
}

/// The source's `computeBoundingBox_`, in 64-bit arithmetic: `(x_min, y_min,
/// width, height)`.
fn bounding_box(rows: &[MaldiPixel]) -> Result<(i64, i64, u32, u32)> {
    let x_min = rows.iter().map(|p| i64::from(p.x)).min().unwrap_or(0);
    let x_max = rows.iter().map(|p| i64::from(p.x)).max().unwrap_or(0);
    let y_min = rows.iter().map(|p| i64::from(p.y)).min().unwrap_or(0);
    let y_max = rows.iter().map(|p| i64::from(p.y)).max().unwrap_or(0);
    let extent = |span: i64| {
        u32::try_from(span + 1).map_err(|_| {
            Error::InvalidValue(format!("imaging extent {} does not fit u32", span + 1))
        })
    };
    Ok((x_min, y_min, extent(x_max - x_min)?, extent(y_max - y_min)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frame_ids_parse_as_the_source() {
        assert_eq!(parse_frame_id("frame=12"), Some(12));
        assert_eq!(parse_frame_id("frame=7 windowGroup=1 scan=0"), Some(7));
        assert_eq!(parse_frame_id("frame=+3"), Some(3));
        assert_eq!(parse_frame_id(" frame=3"), None);
        assert_eq!(parse_frame_id("frame=x"), None);
        assert_eq!(parse_frame_id("frame=++3"), None);
        assert_eq!(parse_frame_id("frame=+-3"), Some(-3));
        assert_eq!(parse_frame_id("frame=-4"), Some(-4));
        assert_eq!(parse_frame_id("frame=3\t"), Some(3));
        assert_eq!(parse_frame_id("frame=3\x0c"), None);
        assert_eq!(parse_frame_id("frame=2147483648"), None);
        assert_eq!(parse_frame_id("scan=3"), None);
    }

    #[test]
    fn doubles_parse_as_the_source() {
        assert_eq!(to_double(" 50 "), Some(50.0));
        assert_eq!(to_double("+2.5"), Some(2.5));
        assert_eq!(to_double("++2"), None);
        assert_eq!(to_double("20um"), None);
        assert_eq!(to_double("1e999"), None);
        assert_eq!(to_double("inf"), Some(f64::INFINITY));
        assert_eq!(to_double("1e-999"), None);
        assert_eq!(to_double("0e5"), Some(0.0));
        assert_eq!(to_double("\x0c5"), None);
        assert!(to_double("NaN").is_some_and(f64::is_nan));
    }
}
