// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reader for OpenSWATH mobilogram Parquet files (`.xim`).
//!
//! Ports `FORMAT/XIMParquetFile.h` / `XIMParquetFile.cpp`, the ion-mobility
//! sibling of [`xic_parquet_file`](crate::format::xic_parquet_file). One or
//! more files are read as one table; mobilograms can be selected by the
//! common metadata fields ([`XimQuery`](crate::format::xim_parquet_file::XimQuery)), by a textual filter or by a typed
//! [`ParquetFilter`](crate::format::parquet_filter::ParquetFilter), and the
//! selection happens before the compressed mobility and intensity arrays are
//! decoded. The `.xim` schema is the one
//! [`MobilogramParquetConsumer`](crate::format::mobilogram_parquet_consumer::MobilogramParquetConsumer)
//! writes ([`XIMSchema`](crate::format::arrow_schema_registry::XIMSchema)).
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
//! `PRODUCT_DECOY`, `TRANSITION_ORDINAL`, `FEATURE_ID`,
//! `MOBILITY_COMPRESSION`, `INTENSITY_COMPRESSION`; string columns:
//! `SOURCE_FILE`, `MOBILOGRAM_TYPE`, `MODIFIED_SEQUENCE`, `TRANSITION_TYPE`,
//! `ANNOTATION`. A column the files do not have is an error. The mobility and
//! intensity arrays are not filterable, and neither is `FEATURE_RT`: a filter
//! naming `MOBILITY`, `INTENSITY` or `FEATURE_RT` is refused (select a
//! feature RT with [`XimQuery::feature_rt`](crate::format::xim_parquet_file::XimQuery::feature_rt) instead).
//! How the operands are compared is described in `docs/OPENSWATH_PARQUET_SUPPORT.md`.
//!
//! ```no_run
//! use openms::format::xim_parquet_file::{XimParquetFile, XimQuery};
//!
//! # fn main() -> openms::Result<()> {
//! let xim = XimParquetFile::new("run.xim")?;
//! let query = XimQuery {
//!     mobilogram_type: "ms2".into(),
//!     filter: "PRECURSOR_ID in [1, 2]".into(),
//!     ..XimQuery::default()
//! };
//! for mobilogram in xim.mobilograms(&query)? {
//!     assert_eq!(mobilogram.mobility.len(), mobilogram.intensity.len());
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
    Analyte as XimAnalyte, ReaderLimits, RunInfo as XimRunInfo,
};
use crate::{Error, Result};
use std::path::{Path, PathBuf};

/// One mobilogram row of a `.xim` file.
///
/// Ports `XIMParquetFile::XIMMobilogram`. Each `has_x` / `x` pair of the
/// source is one `Option`; a null or absent string column reads as `""` and
/// a null or absent `MS_LEVEL` as `0`, as in the source.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct XimMobilogram {
    /// `RUN_ID`.
    pub run_id: i64,
    /// `SOURCE_FILE`.
    pub source_file: String,
    /// `MS_LEVEL`.
    pub ms_level: i64,
    /// `MOBILOGRAM_TYPE`, for example `ms1`, `ms2` or `ms1_contrast`.
    pub mobilogram_type: String,
    /// `PRECURSOR_ID`.
    pub precursor_id: Option<i64>,
    /// `TRANSITION_ID`.
    pub transition_id: Option<i64>,
    /// `FEATURE_ID`.
    pub feature_id: Option<i64>,
    /// `FEATURE_RT`, read from a `double` or `float` column.
    pub feature_rt: Option<f64>,
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
    /// Decoded ion mobilities.
    pub mobility: Vec<f64>,
    /// Decoded intensities, as many as mobilities.
    pub intensity: Vec<f64>,
}

/// Metadata selection for [`XimParquetFile::mobilograms`].
///
/// The optional arguments of the source's `getMobilograms`. The source
/// ignores a negative integer (`-1` by default), a negative `feature_rt` and
/// an empty string; here `None` and `""` ignore a field, and so does a
/// negative value, as in the source. Every selection is combined with AND,
/// and with `filter`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct XimQuery {
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
    /// `MOBILOGRAM_TYPE` equals this (empty to ignore).
    pub mobilogram_type: String,
    /// `FEATURE_ID` equals this.
    pub feature_id: Option<i64>,
    /// `FEATURE_RT` is within 1e-9 of this. Checked row by row only: the
    /// source never turns it into a filter condition.
    pub feature_rt: Option<f64>,
    /// A textual filter (see the [module documentation](self)); empty for
    /// none.
    pub filter: String,
}

/// Reader for OpenSWATH mobilogram Parquet files (`.xim`).
///
/// Ports `OpenMS::XIMParquetFile`. The files are opened on every call, as in
/// the source; the struct only holds their names and the resource ceilings.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct XimParquetFile {
    filename: PathBuf,
    filenames: Vec<PathBuf>,
    /// Resource ceilings for every read (native; see [`ReaderLimits`]).
    pub limits: ReaderLimits,
}

/// `filterPruningOptions_`: `MOBILITY`, `INTENSITY` and `FEATURE_RT` are
/// known but not filterable, and a filter naming one is refused (a filter on
/// `MOBILITY` or `INTENSITY` is an unknown-column error first, since no
/// `.xim` file has columns of those names).
const PRUNING: FilterPruningOptions = FilterPruningOptions {
    unsupported_columns: &["MOBILITY", "INTENSITY", "FEATURE_RT"],
    reject_unsupported_after_pruning: true,
    reject_empty_after_pruning: true,
};

const COLUMN_TYPES: [(&str, ColumnType); 18] = [
    ("RUN_ID", ColumnType::Int),
    ("MS_LEVEL", ColumnType::Int),
    ("PRECURSOR_ID", ColumnType::Int),
    ("TRANSITION_ID", ColumnType::Int),
    ("FEATURE_ID", ColumnType::Int),
    ("PRECURSOR_CHARGE", ColumnType::Int),
    ("PRODUCT_CHARGE", ColumnType::Int),
    ("DETECTING_TRANSITION", ColumnType::Int),
    ("PRECURSOR_DECOY", ColumnType::Int),
    ("PRODUCT_DECOY", ColumnType::Int),
    ("TRANSITION_ORDINAL", ColumnType::Int),
    ("MOBILITY_COMPRESSION", ColumnType::Int),
    ("INTENSITY_COMPRESSION", ColumnType::Int),
    ("SOURCE_FILE", ColumnType::String),
    ("MOBILOGRAM_TYPE", ColumnType::String),
    ("MODIFIED_SEQUENCE", ColumnType::String),
    ("TRANSITION_TYPE", ColumnType::String),
    ("ANNOTATION", ColumnType::String),
];

const CODEC: Codec = Codec {
    pic_zlib: true,
    label: "mobility/intensity",
    unknown: "Unsupported compression type",
};

/// A selection value the source honours: a non-negative integer.
fn selected(value: Option<i64>) -> Option<i64> {
    value.filter(|value| *value >= 0)
}

/// `buildFilterFromArgs_`.
fn args_filter(query: &XimQuery) -> FilterExpression {
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
    if !query.mobilogram_type.is_empty() {
        filter.eq_str("MOBILOGRAM_TYPE", &query.mobilogram_type);
    }
    if let Some(value) = selected(query.feature_id) {
        filter.eq_int("FEATURE_ID", value);
    }
    filter.expression().clone()
}

impl XimParquetFile {
    /// A reader for one `.xim` file.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] of kind `NotFound` when the file does not exist (the
    /// source's `Exception::FileNotFound`).
    pub fn new(filename: impl AsRef<Path>) -> Result<Self> {
        Self::from_files(&[filename])
    }

    /// A reader for several `.xim` files, read as one table.
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

    /// Every mobilogram of the files (`load`).
    ///
    /// # Errors
    ///
    /// As [`mobilograms`](Self::mobilograms).
    pub fn load(&self) -> Result<Vec<XimMobilogram>> {
        self.mobilograms(&XimQuery::default())
    }

    /// The mobilograms `query` selects (`getMobilograms` with its optional
    /// arguments), in file order.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a filter that does not parse, names a
    /// column the files lack or compares incompatible types; for a null
    /// `RUN_ID` ("NULL value in required RUN_ID column"), an unknown
    /// compression identifier ("Unsupported compression type"), a corrupt
    /// array, a mobility/intensity length mismatch ("Mobility and intensity
    /// array size mismatch") and the [`limits`](Self::limits); for a filter
    /// naming `FEATURE_RT` ("Unsupported filter columns after pruning");
    /// [`Error::MissingInformation`] when a required column (`RUN_ID`,
    /// `PRECURSOR_ID`, `MOBILITY_DATA`, `INTENSITY_DATA`,
    /// `MOBILITY_COMPRESSION`, `INTENSITY_COMPRESSION`) is missing; and the
    /// file errors of reading.
    pub fn mobilograms(&self, query: &XimQuery) -> Result<Vec<XimMobilogram>> {
        self.mobilograms_impl(&FilterExpression::default(), query)
    }

    /// The mobilograms a typed filter selects (`getMobilograms(output, const
    /// ParquetFilter&)`).
    ///
    /// # Errors
    ///
    /// As [`mobilograms`](Self::mobilograms).
    pub fn mobilograms_with_filter(&self, filter: &ParquetFilter) -> Result<Vec<XimMobilogram>> {
        self.mobilograms_impl(filter.expression(), &XimQuery::default())
    }

    /// The mobilograms a filter builder selects (`getMobilograms(output,
    /// const ParquetFilterBuilder&)`).
    ///
    /// # Errors
    ///
    /// As [`mobilograms`](Self::mobilograms).
    pub fn mobilograms_with_builder(
        &self,
        filter: &ParquetFilterBuilder,
    ) -> Result<Vec<XimMobilogram>> {
        self.mobilograms_with_filter(filter.filter())
    }

    /// The source's `getMobilograms_`.
    fn mobilograms_impl(
        &self,
        extra: &FilterExpression,
        query: &XimQuery,
    ) -> Result<Vec<XimMobilogram>> {
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
        let mobilogram_type = col("MOBILOGRAM_TYPE", false)?;
        let precursor_id = col("PRECURSOR_ID", true)?;
        let transition_id = col("TRANSITION_ID", false)?;
        let feature_id = col("FEATURE_ID", false)?;
        let feature_rt = col("FEATURE_RT", false)?;
        let modified_sequence = col("MODIFIED_SEQUENCE", false)?;
        let precursor_charge = col("PRECURSOR_CHARGE", false)?;
        let product_charge = col("PRODUCT_CHARGE", false)?;
        let detecting_transition = col("DETECTING_TRANSITION", false)?;
        let precursor_decoy = col("PRECURSOR_DECOY", false)?;
        let product_decoy = col("PRODUCT_DECOY", false)?;
        let transition_ordinal = col("TRANSITION_ORDINAL", false)?;
        let transition_type = col("TRANSITION_TYPE", false)?;
        let annotation = col("ANNOTATION", false)?;
        let mobility_data = col("MOBILITY_DATA", true)?;
        let intensity_data = col("INTENSITY_DATA", true)?;
        let mobility_compression = col("MOBILITY_COMPRESSION", true)?;
        let intensity_compression = col("INTENSITY_COMPRESSION", true)?;
        let mut budget = DecodeBudget::new(self.limits);
        let mut output = Vec::new();
        for row in 0..table.num_rows() {
            let mut mobilogram = XimMobilogram {
                run_id: prh::optional_int(run_id.as_ref(), row)?.ok_or_else(|| {
                    Error::InvalidValue(format!("NULL value in required RUN_ID column: {row}"))
                })?,
                precursor_id: prh::optional_int(precursor_id.as_ref(), row)?,
                transition_id: prh::optional_int(transition_id.as_ref(), row)?,
                feature_id: prh::optional_int(feature_id.as_ref(), row)?,
                feature_rt: prh::optional_double(feature_rt.as_ref(), row)?,
                precursor_charge: prh::optional_int(precursor_charge.as_ref(), row)?,
                product_charge: prh::optional_int(product_charge.as_ref(), row)?,
                detecting_transition: prh::optional_int(detecting_transition.as_ref(), row)?,
                precursor_decoy: prh::optional_int(precursor_decoy.as_ref(), row)?,
                product_decoy: prh::optional_int(product_decoy.as_ref(), row)?,
                transition_ordinal: prh::optional_int(transition_ordinal.as_ref(), row)?,
                source_file: prh::string_or_empty(source_file.as_ref(), row)?,
                mobilogram_type: prh::string_or_empty(mobilogram_type.as_ref(), row)?,
                modified_sequence: prh::string_or_empty(modified_sequence.as_ref(), row)?,
                transition_type: prh::string_or_empty(transition_type.as_ref(), row)?,
                annotation: prh::string_or_empty(annotation.as_ref(), row)?,
                ms_level: prh::optional_int(ms_level.as_ref(), row)?.unwrap_or(0),
                mobility: Vec::new(),
                intensity: Vec::new(),
            };
            if !matches_query(&mobilogram, query) {
                continue;
            }
            let mobility_code = prh::optional_int(mobility_compression.as_ref(), row)?.unwrap_or(0);
            let intensity_code =
                prh::optional_int(intensity_compression.as_ref(), row)?.unwrap_or(0);
            mobilogram.mobility = prh::decode(
                prh::binary(mobility_data.as_ref(), row)?,
                mobility_code,
                CODEC,
                &mut budget,
            )?;
            mobilogram.intensity = prh::decode(
                prh::binary(intensity_data.as_ref(), row)?,
                intensity_code,
                CODEC,
                &mut budget,
            )?;
            if mobilogram.mobility.len() != mobilogram.intensity.len() {
                return Err(Error::InvalidValue(format!(
                    "Mobility and intensity array size mismatch: {row}"
                )));
            }
            output.push(mobilogram);
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
    pub fn runs(&self) -> Result<Vec<XimRunInfo>> {
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
    /// its distinct transitions in the vector members. See [`XimAnalyte`].
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] "Unsupported analyte column", "Column not found
    /// in parquet schema", and, for a nested listing that requests none of
    /// `PRECURSOR_ID`, `MODIFIED_SEQUENCE`, `PRECURSOR_CHARGE` and
    /// `PRECURSOR_DECOY`, "nest_transitions=true requires at least one
    /// precursor discriminator column"; the schema and read errors.
    pub fn analytes<S: AsRef<str>>(
        &self,
        columns: &[S],
        nest_transitions: bool,
    ) -> Result<Vec<XimAnalyte>> {
        prh::analytes(
            &self.filenames,
            columns,
            nest_transitions,
            true,
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
fn matches_query(mobilogram: &XimMobilogram, query: &XimQuery) -> bool {
    let optional = |wanted: Option<i64>, actual: Option<i64>| match selected(wanted) {
        Some(wanted) => actual == Some(wanted),
        None => true,
    };
    let plain = |wanted: Option<i64>, actual: i64| match selected(wanted) {
        Some(wanted) => actual == wanted,
        None => true,
    };
    optional(query.precursor_id, mobilogram.precursor_id)
        && optional(query.transition_id, mobilogram.transition_id)
        && (query.modified_sequence.is_empty()
            || mobilogram.modified_sequence == query.modified_sequence)
        && optional(query.precursor_charge, mobilogram.precursor_charge)
        && optional(query.product_charge, mobilogram.product_charge)
        && plain(query.ms_level, mobilogram.ms_level)
        && plain(query.run_id, mobilogram.run_id)
        && (query.mobilogram_type.is_empty() || mobilogram.mobilogram_type == query.mobilogram_type)
        && optional(query.feature_id, mobilogram.feature_id)
        && match query.feature_rt.filter(|rt| *rt >= 0.0) {
            Some(wanted) => mobilogram
                .feature_rt
                .is_some_and(|actual| (actual - wanted).abs() <= 1e-9),
            None => true,
        }
}
