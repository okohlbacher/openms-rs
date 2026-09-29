// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reader for OpenSWATH chromatogram Parquet files (`.xic`).
//!
//! Ports `FORMAT/XICParquetFile.h` / `XICParquetFile.cpp`. One or more files
//! are read as one table; chromatograms can be selected by the common
//! metadata fields ([`XicQuery`](crate::format::xic_parquet_file::XicQuery)), by a textual filter or by a typed
//! [`ParquetFilter`](crate::format::parquet_filter::ParquetFilter), and the
//! selection happens before the compressed RT and intensity arrays are
//! decoded. The `.xic` schema is the one
//! [`MSChromatogramParquetConsumer`](crate::format::ms_chromatogram_parquet_consumer::MSChromatogramParquetConsumer)
//! writes ([`XICSchema`](crate::format::arrow_schema_registry::XICSchema)).
//!
//! # Filter syntax
//!
//! The textual filter is a left-to-right chain of comparisons over column
//! names, without precedence or parentheses:
//!
//! - comparison: `=`, `==`, `!=`, `<`, `<=`, `>`, `>=`;
//! - set membership: `in [v1, v2, ...]` (or `in (...)`);
//! - connectors: `AND` / `OR`, also `&&`, `||`, `&`, `|`.
//!
//! Values are integers or strings; a string may be unquoted when it holds no
//! space or comma (`annotation=y3^1`), otherwise it is quoted with `'` or
//! `"`. Column names are case-insensitive. Integer columns: `RUN_ID`,
//! `MS_LEVEL`, `PRECURSOR_ID`, `TRANSITION_ID`, `PRECURSOR_CHARGE`,
//! `PRODUCT_CHARGE`, `DETECTING_TRANSITION`, `PRECURSOR_DECOY`,
//! `PRODUCT_DECOY`, `TRANSITION_ORDINAL`, `RT_COMPRESSION`,
//! `INTENSITY_COMPRESSION`; string columns: `SOURCE_FILE`,
//! `MODIFIED_SEQUENCE`, `TRANSITION_TYPE`, `ANNOTATION`. A column the files
//! do not have is an error. The RT and intensity arrays are not filterable.
//! How the operands are compared is described in `docs/OPENSWATH_PARQUET_SUPPORT.md`.
//!
//! ```no_run
//! use openms::format::xic_parquet_file::{XicParquetFile, XicQuery};
//!
//! # fn main() -> openms::Result<()> {
//! let xic = XicParquetFile::new("run.xic")?;
//! let query = XicQuery {
//!     filter: "PRECURSOR_ID=2 OR TRANSITION_ID in [2,3]".into(),
//!     ..XicQuery::default()
//! };
//! for chromatogram in xic.chromatograms(&query)? {
//!     assert_eq!(chromatogram.rt.len(), chromatogram.intensity.len());
//! }
//! # Ok(())
//! # }
//! ```

use crate::format::parquet_filter::{
    ColumnType, FilterExpression, ParquetFilter, ParquetFilterBuilder,
};
use crate::format::parquet_reader_helpers::{
    self as prh, Codec, DecodeBudget, FilterPruningOptions,
};
pub use crate::format::parquet_reader_helpers::{
    Analyte as XicAnalyte, ReaderLimits, RunInfo as XicRunInfo,
};
use crate::{Error, Result};
use std::path::{Path, PathBuf};

/// One chromatogram row of a `.xic` file.
///
/// Ports `XICParquetFile::XICChromatogram`. Each `has_x` / `x` pair of the
/// source is one `Option`; a null or absent string column reads as `""` and
/// a null or absent `MS_LEVEL` as `0`, as in the source.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct XicChromatogram {
    /// `RUN_ID`.
    pub run_id: i64,
    /// `SOURCE_FILE`.
    pub source_file: String,
    /// `MS_LEVEL`: 1 for a precursor trace, 2 for a fragment trace.
    pub ms_level: i64,
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
    /// `ANNOTATION`.
    pub annotation: String,
    /// Decoded retention times.
    pub rt: Vec<f64>,
    /// Decoded intensities, as many as retention times.
    pub intensity: Vec<f64>,
}

/// Metadata selection for [`XicParquetFile::chromatograms`].
///
/// The optional arguments of the source's `getChromatograms`. The source
/// ignores a negative integer (`-1` by default) and an empty string; here
/// `None` and `""` ignore a field, and so does a negative value, as in the
/// source. Every selection is combined with AND, and with `filter`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct XicQuery {
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
    /// A textual filter (see the [module documentation](self)); empty for
    /// none.
    pub filter: String,
}

/// Reader for OpenSWATH chromatogram Parquet files (`.xic`).
///
/// Ports `OpenMS::XICParquetFile`. The files are opened on every call, as in
/// the source; the struct only holds their names and the resource ceilings.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct XicParquetFile {
    filename: PathBuf,
    filenames: Vec<PathBuf>,
    /// Resource ceilings for every read (native; see [`ReaderLimits`]).
    pub limits: ReaderLimits,
}

/// `filterPruningOptions_`: `RT` and `INTENSITY` are known but not
/// filterable (no `.xic` file has columns of those names, so in practice a
/// filter on them is an unknown-column error), dropped without error.
const PRUNING: FilterPruningOptions = FilterPruningOptions {
    unsupported_columns: &["RT", "INTENSITY"],
    reject_unsupported_after_pruning: false,
    reject_empty_after_pruning: false,
};

const COLUMN_TYPES: [(&str, ColumnType); 16] = [
    ("RUN_ID", ColumnType::Int),
    ("MS_LEVEL", ColumnType::Int),
    ("PRECURSOR_ID", ColumnType::Int),
    ("TRANSITION_ID", ColumnType::Int),
    ("PRECURSOR_CHARGE", ColumnType::Int),
    ("PRODUCT_CHARGE", ColumnType::Int),
    ("DETECTING_TRANSITION", ColumnType::Int),
    ("PRECURSOR_DECOY", ColumnType::Int),
    ("PRODUCT_DECOY", ColumnType::Int),
    ("TRANSITION_ORDINAL", ColumnType::Int),
    ("RT_COMPRESSION", ColumnType::Int),
    ("INTENSITY_COMPRESSION", ColumnType::Int),
    ("SOURCE_FILE", ColumnType::String),
    ("MODIFIED_SEQUENCE", ColumnType::String),
    ("TRANSITION_TYPE", ColumnType::String),
    ("ANNOTATION", ColumnType::String),
];

const CODEC: Codec = Codec {
    pic_zlib: true,
    label: "RT/intensity",
    unknown: "Unsupported compression type",
};

/// A selection value the source honours: a non-negative integer.
fn selected(value: Option<i64>) -> Option<i64> {
    value.filter(|value| *value >= 0)
}

/// `buildFilterFromArgs_`.
fn args_filter(query: &XicQuery) -> FilterExpression {
    let mut filter = ParquetFilter::new();
    if let Some(value) = selected(query.precursor_id) {
        filter.eq_int("PRECURSOR_ID", value);
    }
    if let Some(value) = selected(query.transition_id) {
        filter.eq_int("TRANSITION_ID", value);
    }
    if !query.modified_sequence.is_empty() {
        filter.eq_str("MODIFIED_SEQUENCE", &query.modified_sequence);
    }
    if let Some(value) = selected(query.precursor_charge) {
        filter.eq_int("PRECURSOR_CHARGE", value);
    }
    if let Some(value) = selected(query.product_charge) {
        filter.eq_int("PRODUCT_CHARGE", value);
    }
    if let Some(value) = selected(query.ms_level) {
        filter.eq_int("MS_LEVEL", value);
    }
    if let Some(value) = selected(query.run_id) {
        filter.eq_int("RUN_ID", value);
    }
    filter.expression().clone()
}

impl XicParquetFile {
    /// A reader for one `.xic` file.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] of kind `NotFound` when the file does not exist (the
    /// source's `Exception::FileNotFound`).
    pub fn new(filename: impl AsRef<Path>) -> Result<Self> {
        Self::from_files(&[filename])
    }

    /// A reader for several `.xic` files, read as one table.
    ///
    /// The files must share one schema; that is checked when they are read.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] "No parquet files provided" for an empty list;
    /// [`Error::Io`] of kind `NotFound` naming the first missing file.
    pub fn from_files<P: AsRef<Path>>(filenames: &[P]) -> Result<Self> {
        let filenames = prh::owned_paths(filenames);
        let filename = filenames
            .first()
            .cloned()
            .ok_or_else(|| Error::InvalidValue("No parquet files provided".into()))?;
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

    /// Every chromatogram of the files (`load`).
    ///
    /// # Errors
    ///
    /// As [`chromatograms`](Self::chromatograms).
    pub fn load(&self) -> Result<Vec<XicChromatogram>> {
        self.chromatograms(&XicQuery::default())
    }

    /// The chromatograms `query` selects (`getChromatograms` with its
    /// optional arguments), in file order.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a filter that does not parse, names a
    /// column the files lack or compares incompatible types; for a null
    /// `RUN_ID` ("NULL value in required RUN_ID column"), an unknown
    /// compression identifier ("Unsupported compression type"), a corrupt
    /// array, an RT/intensity length mismatch ("RT and intensity array size
    /// mismatch") and the [`limits`](Self::limits);
    /// [`Error::MissingInformation`] when a required column (`RUN_ID`,
    /// `PRECURSOR_ID`, `RT_DATA`, `INTENSITY_DATA`, `RT_COMPRESSION`,
    /// `INTENSITY_COMPRESSION`) is missing; and the file errors of reading.
    pub fn chromatograms(&self, query: &XicQuery) -> Result<Vec<XicChromatogram>> {
        self.chromatograms_impl(&FilterExpression::default(), query)
    }

    /// The chromatograms a typed filter selects (`getChromatograms(output,
    /// const ParquetFilter&)`).
    ///
    /// # Errors
    ///
    /// As [`chromatograms`](Self::chromatograms).
    pub fn chromatograms_with_filter(
        &self,
        filter: &ParquetFilter,
    ) -> Result<Vec<XicChromatogram>> {
        self.chromatograms_impl(filter.expression(), &XicQuery::default())
    }

    /// The chromatograms a filter builder selects (`getChromatograms(output,
    /// const ParquetFilterBuilder&)`).
    ///
    /// # Errors
    ///
    /// As [`chromatograms`](Self::chromatograms).
    pub fn chromatograms_with_builder(
        &self,
        filter: &ParquetFilterBuilder,
    ) -> Result<Vec<XicChromatogram>> {
        self.chromatograms_with_filter(filter.filter())
    }

    /// The source's `getChromatograms_`.
    fn chromatograms_impl(
        &self,
        extra: &FilterExpression,
        query: &XicQuery,
    ) -> Result<Vec<XicChromatogram>> {
        let parsed = prh::parse_filter(&query.filter, &COLUMN_TYPES)?;
        let combined = ParquetFilter::merge(&args_filter(query), &parsed, "AND");
        let combined = ParquetFilter::merge(&combined, extra, "AND");
        let context = if query.filter.is_empty() {
            "<typed/args filter>"
        } else {
            query.filter.as_str()
        };
        let table = prh::read_tables(&self.filenames, None, &self.limits)?;
        let table = prh::apply_filter(table, &combined, context, &PRUNING)?;
        let col = |name: &str, required: bool| prh::column(&table, name, required);
        let run_id = col("RUN_ID", true)?;
        let source_file = col("SOURCE_FILE", false)?;
        let ms_level = col("MS_LEVEL", false)?;
        let precursor_id = col("PRECURSOR_ID", true)?;
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
        let rt_data = col("RT_DATA", true)?;
        let intensity_data = col("INTENSITY_DATA", true)?;
        let rt_compression = col("RT_COMPRESSION", true)?;
        let intensity_compression = col("INTENSITY_COMPRESSION", true)?;
        let mut budget = DecodeBudget::new(self.limits);
        let mut output = Vec::new();
        for row in 0..table.num_rows() {
            let mut chromatogram = XicChromatogram {
                run_id: prh::optional_int(run_id.as_ref(), row)?.ok_or_else(|| {
                    Error::InvalidValue(format!("NULL value in required RUN_ID column: {row}"))
                })?,
                precursor_id: prh::optional_int(precursor_id.as_ref(), row)?,
                transition_id: prh::optional_int(transition_id.as_ref(), row)?,
                precursor_charge: prh::optional_int(precursor_charge.as_ref(), row)?,
                product_charge: prh::optional_int(product_charge.as_ref(), row)?,
                detecting_transition: prh::optional_int(detecting_transition.as_ref(), row)?,
                precursor_decoy: prh::optional_int(precursor_decoy.as_ref(), row)?,
                product_decoy: prh::optional_int(product_decoy.as_ref(), row)?,
                transition_ordinal: prh::optional_int(transition_ordinal.as_ref(), row)?,
                source_file: prh::string_or_empty(source_file.as_ref(), row)?,
                modified_sequence: prh::string_or_empty(modified_sequence.as_ref(), row)?,
                transition_type: prh::string_or_empty(transition_type.as_ref(), row)?,
                annotation: prh::string_or_empty(annotation.as_ref(), row)?,
                ms_level: prh::optional_int(ms_level.as_ref(), row)?.unwrap_or(0),
                rt: Vec::new(),
                intensity: Vec::new(),
            };
            if !matches_query(&chromatogram, query) {
                continue;
            }
            let rt_code = prh::optional_int(rt_compression.as_ref(), row)?.unwrap_or(0);
            let intensity_code =
                prh::optional_int(intensity_compression.as_ref(), row)?.unwrap_or(0);
            chromatogram.rt = prh::decode(
                prh::binary(rt_data.as_ref(), row)?,
                rt_code,
                CODEC,
                &mut budget,
            )?;
            chromatogram.intensity = prh::decode(
                prh::binary(intensity_data.as_ref(), row)?,
                intensity_code,
                CODEC,
                &mut budget,
            )?;
            if chromatogram.rt.len() != chromatogram.intensity.len() {
                return Err(Error::InvalidValue(format!(
                    "RT and intensity array size mismatch: {row}"
                )));
            }
            output.push(chromatogram);
        }
        Ok(output)
    }

    /// Distinct `(RUN_ID, SOURCE_FILE)` pairs in first-seen order
    /// (`getRuns`); no array is decoded.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] without a `RUN_ID` column;
    /// [`Error::InvalidValue`] "NULL value in required RUN_ID column".
    pub fn runs(&self) -> Result<Vec<XicRunInfo>> {
        prh::runs(&self.filenames, &self.limits)
    }

    /// Distinct analytes (`getAnalytes`); no array is decoded.
    ///
    /// `columns` selects analyte columns by case-insensitive name, empty for
    /// all of `PRECURSOR_ID`, `MODIFIED_SEQUENCE`, `PRECURSOR_CHARGE`,
    /// `PRECURSOR_DECOY`, `TRANSITION_ID`, `PRODUCT_CHARGE`,
    /// `TRANSITION_ORDINAL`, `DETECTING_TRANSITION`, `PRODUCT_DECOY`,
    /// `TRANSITION_TYPE` and `ANNOTATION` (the source's default argument).
    /// Without `nest_transitions`, each result is one distinct combination of
    /// the requested columns; with it (the source's default), each result is
    /// one distinct combination of the requested precursor columns, holding
    /// its distinct transitions in the vector members. See [`XicAnalyte`].
    ///
    /// A nested listing that requests no precursor column collapses every
    /// row into one analyte, as in the source.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] "Unsupported analyte column" and "Column not
    /// found in parquet schema"; the schema and read errors.
    pub fn analytes<S: AsRef<str>>(
        &self,
        columns: &[S],
        nest_transitions: bool,
    ) -> Result<Vec<XicAnalyte>> {
        prh::analytes(
            &self.filenames,
            columns,
            nest_transitions,
            false,
            &self.limits,
        )
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

/// The source's per-row re-check of the query arguments.
fn matches_query(chromatogram: &XicChromatogram, query: &XicQuery) -> bool {
    let optional = |wanted: Option<i64>, actual: Option<i64>| match selected(wanted) {
        Some(wanted) => actual == Some(wanted),
        None => true,
    };
    let plain = |wanted: Option<i64>, actual: i64| match selected(wanted) {
        Some(wanted) => actual == wanted,
        None => true,
    };
    optional(query.precursor_id, chromatogram.precursor_id)
        && optional(query.transition_id, chromatogram.transition_id)
        && (query.modified_sequence.is_empty()
            || chromatogram.modified_sequence == query.modified_sequence)
        && optional(query.precursor_charge, chromatogram.precursor_charge)
        && optional(query.product_charge, chromatogram.product_charge)
        && plain(query.ms_level, chromatogram.ms_level)
        && plain(query.run_id, chromatogram.run_id)
}
