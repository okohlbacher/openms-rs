// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reader for extracted ion peak-map Parquet files (`.xipm`).
//!
//! Ports `FORMAT/XIPMParquetFile.h` / `XIPMParquetFile.cpp`. A `.xipm` row is
//! one targeted m/z–RT–ion-mobility peak cloud with its precursor and
//! transition annotation and its extraction window, as `XIPMParquetConsumer`
//! writes it ([`XIPMSchema`](crate::format::arrow_schema_registry::XIPMSchema)).
//! Unlike its `.xic` and `.xim` siblings the reader has no filter
//! expressions: peak maps are selected by the metadata fields of
//! [`XipmQuery`] alone, and every file is checked against the XIPM schema in
//! subset mode (extra columns allowed, and missing ones too, but a present
//! column must have its registered type).
//!
//! ```no_run
//! use openms::format::xipm_parquet_file::{XipmParquetFile, XipmQuery};
//!
//! # fn main() -> openms::Result<()> {
//! let xipm = XipmParquetFile::new("run.xipm")?;
//! let query = XipmQuery {
//!     peakmap_type: "precursor".into(),
//!     ..XipmQuery::default()
//! };
//! for peak_map in xipm.peak_maps(&query)? {
//!     assert_eq!(peak_map.mz.len(), peak_map.intensity.len());
//! }
//! # Ok(())
//! # }
//! ```

use crate::format::arrow_schema_registry::XIPMSchema;
use crate::format::arrow_schema_registry::validation::{self, Mode};
use crate::format::parquet_reader_helpers::{self as prh, Codec, DecodeBudget};
pub use crate::format::parquet_reader_helpers::{ReaderLimits, RunInfo as XipmRunInfo};
use crate::{Error, Result};
use arrow_array::RecordBatch;
use std::collections::HashSet;
use std::path::{Path, PathBuf};

/// One peak-map row of a `.xipm` file.
///
/// Ports `XIPMParquetFile::XIPMPeakMap`. Each `has_x` / `x` pair of the
/// source is one `Option`. A null or absent string reads as `""`, and a null
/// or absent `RUN_ID`, `MS_LEVEL` or `TARGET_MZ` as zero, as in the source.
///
/// When the row has no `TARGET_RT` but has an extraction window whose end
/// lies after its start, [`target_rt`](Self::target_rt) is the window's
/// midpoint, as in the source.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct XipmPeakMap {
    /// `RUN_ID`.
    pub run_id: i64,
    /// `SOURCE_FILE`.
    pub source_file: String,
    /// `MS_LEVEL`: 1 for a precursor peak map, 2 for a fragment peak map.
    pub ms_level: i64,
    /// `PEAKMAP_TYPE`: `precursor` or `transition`.
    pub peakmap_type: String,
    /// `PRECURSOR_ID`.
    pub precursor_id: Option<i64>,
    /// `TRANSITION_ID`.
    pub transition_id: Option<i64>,
    /// `MODIFIED_SEQUENCE`.
    pub modified_sequence: String,
    /// `PRECURSOR_CHARGE`.
    pub precursor_charge: Option<i64>,
    /// `PRODUCT_CHARGE`.
    pub product_charge: Option<i64>,
    /// `DETECTING_TRANSITION`.
    pub detecting_transition: Option<i64>,
    /// `PRECURSOR_DECOY`.
    pub precursor_decoy: Option<i64>,
    /// `PRODUCT_DECOY`.
    pub product_decoy: Option<i64>,
    /// `TRANSITION_ORDINAL`.
    pub transition_ordinal: Option<i64>,
    /// `TRANSITION_TYPE`.
    pub transition_type: String,
    /// `ANNOTATION`: the transition annotation (`y7^1`) or the precursor's
    /// (`Precursor_i0`).
    pub annotation: String,
    /// `TARGET_MZ`, the extraction target m/z.
    pub target_mz: f64,
    /// `TARGET_RT`, or the extraction window's midpoint.
    pub target_rt: Option<f64>,
    /// `TARGET_ION_MOBILITY`.
    pub target_ion_mobility: Option<f64>,
    /// `RT_START`, the extraction window's start.
    pub rt_start: Option<f64>,
    /// `RT_END`, the extraction window's end.
    pub rt_end: Option<f64>,
    /// Decoded m/z values.
    pub mz: Vec<f64>,
    /// Decoded retention times.
    pub rt: Vec<f64>,
    /// Decoded ion mobilities.
    pub ion_mobility: Vec<f64>,
    /// Decoded intensities; all four arrays have one length.
    pub intensity: Vec<f64>,
}

/// Metadata selection for [`XipmParquetFile::peak_maps`].
///
/// The optional arguments of the source's `getPeakMaps`, evaluated row by
/// row. `None`, `""` and a negative value ignore a field, as the source's
/// defaults (`-1`, `""`) do.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct XipmQuery {
    /// `PRECURSOR_ID` equals this.
    pub precursor_id: Option<i64>,
    /// `TRANSITION_ID` equals this.
    pub transition_id: Option<i64>,
    /// `MODIFIED_SEQUENCE` equals this (empty to ignore).
    pub modified_sequence: String,
    /// `PRECURSOR_CHARGE` equals this.
    pub precursor_charge: Option<i64>,
    /// `PRODUCT_CHARGE` equals this.
    pub product_charge: Option<i64>,
    /// `MS_LEVEL` equals this.
    pub ms_level: Option<i64>,
    /// `RUN_ID` equals this.
    pub run_id: Option<i64>,
    /// `PEAKMAP_TYPE` equals this (empty to ignore).
    pub peakmap_type: String,
}

/// Reader for extracted ion peak-map Parquet files (`.xipm`).
///
/// Ports `OpenMS::XIPMParquetFile`. The files are opened on every call; the
/// struct only holds their names and the resource ceilings.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct XipmParquetFile {
    filename: PathBuf,
    filenames: Vec<PathBuf>,
    /// Resource ceilings for every read (native; see [`ReaderLimits`]).
    pub limits: ReaderLimits,
}

const CODEC: Codec = Codec {
    pic_zlib: false,
    label: "binary",
    unknown: "Unknown compression code",
};

/// `validateXIPMTable_`: the table against the XIPM schema, subset mode.
fn validate(table: &RecordBatch) -> Result<()> {
    let result = validation::validate(table, &XIPMSchema::schema(), Mode::Subset);
    if result.valid {
        Ok(())
    } else {
        Err(Error::InvalidValue(format!(
            "Invalid XIPM parquet schema: {result}"
        )))
    }
}

fn selected(value: Option<i64>) -> Option<i64> {
    value.filter(|value| *value >= 0)
}

impl XipmParquetFile {
    /// A reader for one `.xipm` file.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] of kind `NotFound` when the file does not exist (the
    /// source's `Exception::FileNotFound`).
    pub fn new(filename: impl AsRef<Path>) -> Result<Self> {
        Self::from_files(&[filename])
    }

    /// A reader for several `.xipm` files, read as one table.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] of kind `NotFound` for an empty list ("No .xipm files
    /// given", the source's `Exception::FileNotFound`) and naming the first
    /// missing file.
    pub fn from_files<P: AsRef<Path>>(filenames: &[P]) -> Result<Self> {
        let filenames = prh::owned_paths(filenames);
        let filename = filenames.first().cloned().ok_or_else(|| {
            Error::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "No .xipm files given",
            ))
        })?;
        prh::require_files(&filenames)?;
        Ok(Self {
            filename,
            filenames,
            limits: ReaderLimits::default(),
        })
    }

    /// The primary file: the first of several (`getFilename`).
    pub fn filename(&self) -> &Path {
        &self.filename
    }

    /// Every file of this reader, in the order given (`getFilenames`).
    pub fn filenames(&self) -> &[PathBuf] {
        &self.filenames
    }

    /// Every peak map of the files (`load`).
    ///
    /// # Errors
    ///
    /// As [`peak_maps`](Self::peak_maps).
    pub fn load(&self) -> Result<Vec<XipmPeakMap>> {
        self.peak_maps(&XipmQuery::default())
    }

    /// The peak maps `query` selects (`getPeakMaps`), in file order.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] "Invalid XIPM parquet schema" when a column
    /// has another type than the XIPM schema registers, for an unknown
    /// compression identifier ("Unknown compression code"), a corrupt array,
    /// arrays of different lengths ("Peak-map array size mismatch") and the
    /// [`limits`](Self::limits); [`Error::MissingInformation`] when a required
    /// column (`RUN_ID`, `MS_LEVEL`, `TARGET_MZ`, the four `*_DATA` and four
    /// `*_COMPRESSION` columns) is missing from a non-empty table; and the
    /// file errors of reading.
    pub fn peak_maps(&self, query: &XipmQuery) -> Result<Vec<XipmPeakMap>> {
        let table = prh::read_tables(&self.filenames, None, &self.limits)?;
        validate(&table)?;
        if table.num_rows() == 0 {
            return Ok(Vec::new());
        }
        let col = |name: &str, required: bool| prh::column(&table, name, required);
        let run_id = col("RUN_ID", true)?;
        let source_file = col("SOURCE_FILE", false)?;
        let ms_level = col("MS_LEVEL", true)?;
        let peakmap_type = col("PEAKMAP_TYPE", false)?;
        let precursor_id = col("PRECURSOR_ID", false)?;
        let transition_id = col("TRANSITION_ID", false)?;
        let modified_sequence = col("MODIFIED_SEQUENCE", false)?;
        let precursor_charge = col("PRECURSOR_CHARGE", false)?;
        let product_charge = col("PRODUCT_CHARGE", false)?;
        let detecting_transition = col("DETECTING_TRANSITION", false)?;
        let precursor_decoy = col("PRECURSOR_DECOY", false)?;
        let product_decoy = col("PRODUCT_DECOY", false)?;
        let transition_ordinal = col("TRANSITION_ORDINAL", false)?;
        let transition_type = col("TRANSITION_TYPE", false)?;
        let annotation = col("ANNOTATION", false)?;
        let target_mz = col("TARGET_MZ", true)?;
        let target_rt = col("TARGET_RT", false)?;
        let target_im = col("TARGET_ION_MOBILITY", false)?;
        let rt_start = col("RT_START", false)?;
        let rt_end = col("RT_END", false)?;
        let mz_data = col("MZ_DATA", true)?;
        let rt_data = col("RT_DATA", true)?;
        let mobility_data = col("MOBILITY_DATA", true)?;
        let intensity_data = col("INTENSITY_DATA", true)?;
        let mz_compression = col("MZ_COMPRESSION", true)?;
        let rt_compression = col("RT_COMPRESSION", true)?;
        let mobility_compression = col("MOBILITY_COMPRESSION", true)?;
        let intensity_compression = col("INTENSITY_COMPRESSION", true)?;
        let mut budget = DecodeBudget::new(self.limits);
        let mut output = Vec::new();
        for row in 0..table.num_rows() {
            let int =
                |array: &Option<arrow_array::ArrayRef>| prh::optional_int(array.as_ref(), row);
            let real =
                |array: &Option<arrow_array::ArrayRef>| prh::optional_double(array.as_ref(), row);
            let text =
                |array: &Option<arrow_array::ArrayRef>| prh::string_or_empty(array.as_ref(), row);
            let mut peak_map = XipmPeakMap {
                run_id: int(&run_id)?.unwrap_or(0),
                source_file: text(&source_file)?,
                peakmap_type: text(&peakmap_type)?,
                precursor_id: int(&precursor_id)?,
                transition_id: int(&transition_id)?,
                modified_sequence: text(&modified_sequence)?,
                precursor_charge: int(&precursor_charge)?,
                product_charge: int(&product_charge)?,
                detecting_transition: int(&detecting_transition)?,
                precursor_decoy: int(&precursor_decoy)?,
                product_decoy: int(&product_decoy)?,
                transition_ordinal: int(&transition_ordinal)?,
                transition_type: text(&transition_type)?,
                annotation: text(&annotation)?,
                ms_level: int(&ms_level)?.unwrap_or(0),
                target_mz: real(&target_mz)?.unwrap_or(0.0),
                target_rt: real(&target_rt)?,
                target_ion_mobility: real(&target_im)?,
                rt_start: real(&rt_start)?,
                rt_end: real(&rt_end)?,
                ..XipmPeakMap::default()
            };
            if peak_map.target_rt.is_none() {
                if let (Some(start), Some(end)) = (peak_map.rt_start, peak_map.rt_end) {
                    if end > start {
                        peak_map.target_rt = Some((start + end) / 2.0);
                    }
                }
            }
            if !matches_query(&peak_map, query) {
                continue;
            }
            let code = |array: &Option<arrow_array::ArrayRef>| -> Result<i64> {
                Ok(prh::optional_int(array.as_ref(), row)?.unwrap_or(0))
            };
            let mz_code = code(&mz_compression)?;
            let rt_code = code(&rt_compression)?;
            let mobility_code = code(&mobility_compression)?;
            let intensity_code = code(&intensity_compression)?;
            peak_map.mz = prh::decode(
                prh::binary(mz_data.as_ref(), row)?,
                mz_code,
                CODEC,
                &mut budget,
            )?;
            peak_map.rt = prh::decode(
                prh::binary(rt_data.as_ref(), row)?,
                rt_code,
                CODEC,
                &mut budget,
            )?;
            peak_map.ion_mobility = prh::decode(
                prh::binary(mobility_data.as_ref(), row)?,
                mobility_code,
                CODEC,
                &mut budget,
            )?;
            peak_map.intensity = prh::decode(
                prh::binary(intensity_data.as_ref(), row)?,
                intensity_code,
                CODEC,
                &mut budget,
            )?;
            let n = peak_map.mz.len();
            if peak_map.rt.len() != n
                || peak_map.ion_mobility.len() != n
                || peak_map.intensity.len() != n
            {
                return Err(Error::InvalidValue(format!(
                    "Peak-map array size mismatch: {row}"
                )));
            }
            output.push(peak_map);
        }
        Ok(output)
    }

    /// Distinct `(RUN_ID, SOURCE_FILE)` pairs in first-seen order
    /// (`getRuns`); no array is decoded.
    ///
    /// Unlike the `.xic` and `.xim` readers, a null `RUN_ID` reads as 0, as
    /// in the source.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when the files have neither column, or
    /// rows but no `RUN_ID`; [`Error::InvalidValue`] "Invalid XIPM parquet
    /// schema".
    pub fn runs(&self) -> Result<Vec<XipmRunInfo>> {
        let table = prh::read_tables(
            &self.filenames,
            Some(&["RUN_ID", "SOURCE_FILE"]),
            &self.limits,
        )?;
        validate(&table)?;
        if table.num_rows() == 0 {
            return Ok(Vec::new());
        }
        let run_id = prh::column(&table, "RUN_ID", true)?;
        let source_file = prh::column(&table, "SOURCE_FILE", false)?;
        let mut seen = HashSet::new();
        let mut output = Vec::new();
        for row in 0..table.num_rows() {
            let info = XipmRunInfo {
                run_id: prh::optional_int(run_id.as_ref(), row)?.unwrap_or(0),
                source_file: prh::string_or_empty(source_file.as_ref(), row)?,
            };
            if seen.insert(format!("{}\t{}", info.run_id, info.source_file)) {
                output.push(info);
            }
        }
        Ok(output)
    }

    /// The column names of the files' schema (`getColumns`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a schema cannot be read or the files'
    /// schemas differ.
    pub fn columns(&self) -> Result<Vec<String>> {
        prh::column_names(&self.filenames)
    }
}

/// The source's per-row selection.
fn matches_query(peak_map: &XipmPeakMap, query: &XipmQuery) -> bool {
    let optional = |wanted: Option<i64>, actual: Option<i64>| match selected(wanted) {
        Some(wanted) => actual == Some(wanted),
        None => true,
    };
    let plain = |wanted: Option<i64>, actual: i64| match selected(wanted) {
        Some(wanted) => actual == wanted,
        None => true,
    };
    optional(query.precursor_id, peak_map.precursor_id)
        && optional(query.transition_id, peak_map.transition_id)
        && (query.modified_sequence.is_empty()
            || peak_map.modified_sequence == query.modified_sequence)
        && optional(query.precursor_charge, peak_map.precursor_charge)
        && optional(query.product_charge, peak_map.product_charge)
        && plain(query.ms_level, peak_map.ms_level)
        && plain(query.run_id, peak_map.run_id)
        && (query.peakmap_type.is_empty() || peak_map.peakmap_type == query.peakmap_type)
}
