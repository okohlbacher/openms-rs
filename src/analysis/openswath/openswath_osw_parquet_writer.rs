// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Write OpenSwathWorkflow results as an OSW Parquet container (`.oswpq`).
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.h` and its
//! translation unit. The writer mirrors the OSW SQLite tables in Parquet: the
//! assay library (written with
//! [`TransitionParquetFile`](crate::analysis::openswath::transition_parquet_file::TransitionParquetFile)),
//! a run index, and per run the scored features, their MS1 isotope
//! intensities and their per-transition scores (identification transitions
//! included when UIS scoring is on):
//!
//! ```text
//! <output>.oswpq
//! ├── metadata.json
//! ├── library/
//! │   ├── precursors.parquet
//! │   ├── transitions.parquet
//! │   └── metadata.json
//! └── runs/
//!     ├── runs.parquet
//!     └── run_id=<id>/
//!         ├── features.parquet
//!         ├── feature_precursor.parquet
//!         └── feature_transition.parquet
//! ```
//!
//! The output is a directory, or a store-only ZIP archive of that layout
//! with a JSON size index named `<archive name>.idx.json`. Writing into a
//! container that already exists adds a new `run_id=<id>` partition,
//! extends `runs.parquet` and refreshes `metadata.json`; existing runs are
//! not modified, and the library must describe exactly the same canonical
//! precursor and transition ids.
//!
//! Feature values come from the meta values the OpenSWATH scoring attaches
//! to each feature and to its subordinates (one per transition, and one per
//! MS1 isotope). Every numeric meta value is read as the source reads it:
//! rendered as `DataValue::toString()` and parsed back with
//! `StringUtils::toDouble`, so a value survives with the 15 significant
//! digits that rendering keeps, and a text value that parses counts as a
//! number. A missing, empty, unparseable or non-finite value is written as
//! null.
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_PARQUET_SUPPORT.md).
//!
//! # Differences from the source
//!
//! - The three run tables are written one after the other; the source writes
//!   them concurrently (`std::async`). Either way a failure removes the run
//!   directory, and the files written are the same.
//! - A container that is an archive is assembled with its entries in sorted
//!   path order; the source follows `std::filesystem::recursive_directory_iterator`,
//!   whose order the file system decides.
//! - `runs.parquet` is read through the type-checking column accessors; the
//!   source casts its `run_id` and `filename` columns unchecked
//!   (`std::static_pointer_cast`), which is undefined for any other type.
//! - The warning the source writes to `OPENMS_LOG_WARN` when it appends to an
//!   existing archive is returned to the caller instead.
//! - `metadata.json`'s `build_time` is written as `""`, as in
//!   [`TransitionParquetFile`](crate::analysis::openswath::transition_parquet_file::TransitionParquetFile).
//! - `feature_precursor.parquet`'s `precursor_apex_intensity` is the MS1
//!   subordinate's own `peak_apex_int`. The source writes
//!   `appendOptionalFloat_(b, extractMetaDouble_(sub, "peak_apex_int", value), value, ..)`
//!   with `value` passed by value, and C++ leaves the order in which the
//!   arguments are evaluated unspecified: the executed Linux Release build
//!   (GCC 14.4, `-O2`) copies `value` before the extraction runs and writes
//!   the double extracted before it on the same feature (measured: `7.0`, a
//!   preceding transition's `points_across_half_height`, where the
//!   subordinate carries `3.25`). The null pattern agrees. Every other
//!   column written the same way agrees with the Release build.
//!
//! Reproduced source behaviour worth knowing: an error raised while the
//! features are converted (a missing `PeptideRef`, an unknown transition)
//! comes after the run directory was created, and the directory stays
//! behind, so the same `run_id` is then refused until it is removed; the
//! same holds for an error while `runs.parquet` is updated.

use crate::analysis::openswath::openswath_library_id_normalizer::{
    OpenSwathLibraryIDNormalizer, SourceIDMapping,
};
use crate::analysis::openswath::transition_experiment::LightTargetedExperiment;
use crate::analysis::openswath::transition_parquet_file::TransitionParquetFile;
use crate::format::arrow_io_helpers::data_value_text;
use crate::format::arrow_schema_registry::validation::{self, Mode};
use crate::format::arrow_schema_registry::{
    OSWFeaturePrecursorSchema, OSWFeatureSchema, OSWFeatureTransitionSchema, OSWPrecursorSchema,
    OSWRunSchema, OSWTransitionSchema,
};
use crate::format::parquet_file::ParquetFile;
use crate::format::zip_archive_file::ZipArchiveFile;
use crate::kernel::{Feature, FeatureMap};
use crate::metadata::{MetaInfo, MetaValueData};
use crate::param::value::format_float;
use crate::system::file::{self, CopyOptions, TempDir};
use crate::{Error, Result};
use arrow_array::{ArrayRef, Float64Array, Int32Array, Int64Array, RecordBatch, StringArray};
use arrow_schema::{DataType, SchemaRef};
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Writer of OSW Parquet containers.
///
/// Ports `OpenSwathOSWParquetWriter`. Its one setting,
/// [`set_preserve_existing`](Self::set_preserve_existing), decides whether
/// an existing archive (not a directory) is extended or replaced.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct OpenSwathOSWParquetWriter {
    preserve_existing: bool,
}

/// Per-run row counts of the metadata (the source's `RunCounts`).
#[derive(Clone, Copy, Default)]
struct RunCounts {
    features: i64,
    feature_precursor: i64,
    feature_transition: i64,
}

/// A table under construction: integer and nullable double columns by name,
/// assembled in the order of a registry schema.
struct TableBuilder {
    schema: SchemaRef,
    int64: HashMap<&'static str, Vec<i64>>,
    int32: HashMap<&'static str, Vec<i32>>,
    float64: HashMap<String, Vec<Option<f64>>>,
}

impl TableBuilder {
    fn new(schema: SchemaRef) -> Self {
        Self {
            schema,
            int64: HashMap::new(),
            int32: HashMap::new(),
            float64: HashMap::new(),
        }
    }

    fn i64(&mut self, column: &'static str, value: i64) {
        self.int64.entry(column).or_default().push(value);
    }

    fn i32(&mut self, column: &'static str, value: i32) {
        self.int32.entry(column).or_default().push(value);
    }

    fn f64(&mut self, column: &str, value: Option<f64>) {
        self.float64
            .entry(column.to_owned())
            .or_default()
            .push(value);
    }

    /// The table in schema order; a column never appended to is empty.
    fn finish(mut self, what: &str) -> Result<RecordBatch> {
        let mut columns: Vec<ArrayRef> = Vec::with_capacity(self.schema.fields().len());
        for field in self.schema.fields() {
            let name = field.name().as_str();
            let column: ArrayRef = match field.data_type() {
                DataType::Int64 => Arc::new(Int64Array::from(
                    self.int64.remove(name).unwrap_or_default(),
                )),
                DataType::Int32 => Arc::new(Int32Array::from(
                    self.int32.remove(name).unwrap_or_default(),
                )),
                DataType::Float64 => Arc::new(Float64Array::from(
                    self.float64.remove(name).unwrap_or_default(),
                )),
                other => {
                    return Err(Error::InvalidValue(format!(
                        "{what}: unexpected column type {other} for '{name}'"
                    )));
                }
            };
            columns.push(column);
        }
        if !self.int64.is_empty() || !self.int32.is_empty() || !self.float64.is_empty() {
            return Err(Error::InvalidValue(format!(
                "{what}: values for a column the schema does not have"
            )));
        }
        let table = RecordBatch::try_new(self.schema.clone(), columns)
            .map_err(|e| Error::InvalidValue(format!("{what}: {e}")))?;
        let result = validation::validate(&table, self.schema.as_ref(), Mode::Strict);
        if !result.valid {
            return Err(Error::InvalidValue(format!(
                "{what} schema validation failed: {result}"
            )));
        }
        Ok(table)
    }
}

/// `StringUtils::toDouble` (`StringUtils.cpp:239-276`): whether `text`
/// converts, and to what.
///
/// Surrounding space, tab, CR and LF are skipped; `nan` (any case, with an
/// optional parenthesised payload) is accepted first; then one `+` may
/// precede a `std::from_chars` general-format number: an optional `-`,
/// digits with an optional fraction and exponent, or `inf`/`infinity`/`nan`.
/// A finite literal that overflows fails, as `result_out_of_range` does.
fn to_double(text: &str) -> Option<f64> {
    let bytes = text.as_bytes();
    let first = skip_whitespace(bytes, 0);
    if first == bytes.len() {
        return None;
    }
    if let Some(length) = word_at(bytes, first, "nan") {
        let mut end = first + length;
        if bytes.get(end) == Some(&b'(') {
            match bytes[end..].iter().position(|b| *b == b')') {
                Some(close) => end += close + 1,
                None => end = first,
            }
        }
        if end != first && skip_whitespace(bytes, end) == bytes.len() {
            return Some(f64::NAN);
        }
    }
    let mut start = first;
    if bytes[start] == b'+' {
        start += 1;
    }
    let mut cursor = start;
    if bytes.get(cursor) == Some(&b'-') {
        cursor += 1;
    }
    let infinity = word_at(bytes, cursor, "infinity").or_else(|| word_at(bytes, cursor, "inf"));
    let nan = word_at(bytes, cursor, "nan");
    if let Some(length) = infinity.or(nan) {
        cursor += length;
    } else {
        let integer_start = cursor;
        while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
            cursor += 1;
        }
        let mut digits = cursor - integer_start;
        if bytes.get(cursor) == Some(&b'.') {
            let fraction_start = cursor + 1;
            let mut fraction_end = fraction_start;
            while bytes.get(fraction_end).is_some_and(u8::is_ascii_digit) {
                fraction_end += 1;
            }
            if digits > 0 || fraction_end > fraction_start {
                digits += fraction_end - fraction_start;
                cursor = fraction_end;
            }
        }
        if digits == 0 {
            return None;
        }
        if matches!(bytes.get(cursor), Some(b'e' | b'E')) {
            let mut exponent = cursor + 1;
            if matches!(bytes.get(exponent), Some(b'+' | b'-')) {
                exponent += 1;
            }
            let exponent_digits = exponent;
            while bytes.get(exponent).is_some_and(u8::is_ascii_digit) {
                exponent += 1;
            }
            if exponent > exponent_digits {
                cursor = exponent;
            }
        }
    }
    let literal = text.get(start..cursor)?;
    let value: f64 = if nan.is_some() && infinity.is_none() {
        f64::NAN
    } else {
        literal.parse().ok()?
    };
    if value.is_infinite() && infinity.is_none() {
        return None;
    }
    (skip_whitespace(bytes, cursor) == bytes.len()).then_some(value)
}

/// `StringUtils::toInt64` (`StringUtils.cpp:168-196`): surrounding
/// whitespace skipped, one `+` allowed, then a complete decimal `i64`.
pub(crate) fn to_int64(text: &str) -> Option<i64> {
    let bytes = text.as_bytes();
    let mut cursor = skip_whitespace(bytes, 0);
    if cursor == bytes.len() {
        return None;
    }
    if bytes[cursor] == b'+' {
        cursor += 1;
    }
    let number_start = cursor;
    if bytes.get(cursor) == Some(&b'-') {
        cursor += 1;
    }
    let digits_start = cursor;
    while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
        cursor += 1;
    }
    if cursor == digits_start {
        return None;
    }
    let value = text.get(number_start..cursor)?.parse::<i64>().ok()?;
    (skip_whitespace(bytes, cursor) == bytes.len()).then_some(value)
}

/// The source's four whitespace characters, skipped from `index` on.
fn skip_whitespace(bytes: &[u8], mut index: usize) -> usize {
    while bytes
        .get(index)
        .is_some_and(|b| matches!(*b, b' ' | b'\t' | b'\n' | b'\r'))
    {
        index += 1;
    }
    index
}

/// Length of a case-insensitive ASCII `word` at `bytes[index..]`, if present.
fn word_at(bytes: &[u8], index: usize, word: &str) -> Option<usize> {
    let candidate = bytes.get(index..index.checked_add(word.len())?)?;
    candidate
        .eq_ignore_ascii_case(word.as_bytes())
        .then_some(word.len())
}

/// `extractMetaDouble_`: a present, non-empty meta value that converts.
fn meta_double(meta: &MetaInfo, key: &str) -> Option<f64> {
    let value = meta.get(key)?;
    if value.is_empty() {
        return None;
    }
    to_double(&data_value_text(value))
}

/// `appendOptionalFloat_`: the value when present and finite, else null.
fn finite(value: Option<f64>) -> Option<f64> {
    value.filter(|v| v.is_finite())
}

/// `appendFeatureScore_`: the meta value, else the finite fallback, else null.
/// A meta value that converts to a non-finite number is null, not the
/// fallback.
fn feature_score(meta: &MetaInfo, key: &str, fallback: Option<f64>) -> Option<f64> {
    match meta_double(meta, key) {
        Some(value) => finite(Some(value)),
        None => finite(fallback),
    }
}

/// `getSeparateScore_`: a list meta value as strings (numbers rendered by
/// `StringUtils::toStr`), a scalar as its one string, missing or empty as
/// nothing.
fn separate_score(meta: &MetaInfo, key: &str) -> Vec<String> {
    let Some(value) = meta.get(key) else {
        return Vec::new();
    };
    match value.data() {
        MetaValueData::Empty => Vec::new(),
        MetaValueData::StringList(values) => values.clone(),
        MetaValueData::IntegerList(values) => values.iter().map(i64::to_string).collect(),
        MetaValueData::FloatList(values) => values.iter().map(|v| format_float(*v, true)).collect(),
        _ => vec![data_value_text(value)],
    }
}

/// `appendOptionalFloatFromList_`: element `index` when present and it
/// converts to a finite number, else null.
fn from_list(values: &[String], index: usize) -> Option<f64> {
    finite(values.get(index).and_then(|text| to_double(text)))
}

/// `getMetaValue(key) == "<level>"`: a string meta value without a unit.
fn has_level(meta: &MetaInfo, level: &str) -> bool {
    meta.get("FeatureLevel").is_some_and(|value| {
        value.unit().is_none() && matches!(value.data(), MetaValueData::String(s) if s == level)
    })
}

/// `StringUtils::toStr(getMetaValue(key))`: the text, empty when missing.
fn meta_text(meta: &MetaInfo, key: &str) -> String {
    meta.get(key).map(data_value_text).unwrap_or_default()
}

/// `SqliteHelper::clearSignBit`: the value with its top bit cleared.
fn clear_sign_bit(value: u64) -> i64 {
    // Cannot fail: the top bit is clear.
    i64::try_from(value & !(1u64 << 63)).unwrap_or(i64::MAX)
}

/// The feature-level double columns whose value is a meta value of the
/// feature, written only when the feature carries MS1 scores
/// (`var_ms1_ppm_diff`): (column, meta key).
const MS1_COLUMNS: [(&str, &str); 17] = [
    ("ms1_area_intensity", "ms1_area_intensity"),
    ("ms1_apex_intensity", "ms1_apex_intensity"),
    ("ms1_exp_im", "im_ms1_drift"),
    ("ms1_delta_im", "im_ms1_delta"),
    ("var_ms1_massdev_score", "var_ms1_ppm_diff"),
    ("var_ms1_im_ms1_delta_score", "var_im_ms1_delta_score"),
    ("var_ms1_mi_score", "var_ms1_mi_score"),
    ("var_ms1_mi_contrast_score", "var_ms1_mi_contrast_score"),
    ("var_ms1_mi_combined_score", "var_ms1_mi_combined_score"),
    (
        "var_ms1_isotope_correlation_score",
        "var_ms1_isotope_correlation",
    ),
    ("var_ms1_isotope_overlap_score", "var_ms1_isotope_overlap"),
    ("var_ms1_xcorr_coelution", "var_ms1_xcorr_coelution"),
    (
        "var_ms1_xcorr_coelution_contrast",
        "var_ms1_xcorr_coelution_contrast",
    ),
    (
        "var_ms1_xcorr_coelution_combined",
        "var_ms1_xcorr_coelution_combined",
    ),
    ("var_ms1_xcorr_shape", "var_ms1_xcorr_shape"),
    (
        "var_ms1_xcorr_shape_contrast",
        "var_ms1_xcorr_shape_contrast",
    ),
    (
        "var_ms1_xcorr_shape_combined",
        "var_ms1_xcorr_shape_combined",
    ),
];

/// The MS2 feature-level double columns read from meta values after
/// `ms2_area_intensity` (the feature intensity): (column, meta key).
const MS2_COLUMNS: [(&str, &str); 36] = [
    ("ms2_total_area_intensity", "total_xic"),
    ("ms2_apex_intensity", "peak_apices_sum"),
    ("ms2_exp_im", "im_drift"),
    ("ms2_exp_im_leftwidth", "im_drift_left"),
    ("ms2_exp_im_rightwidth", "im_drift_right"),
    ("ms2_delta_im", "im_delta"),
    ("ms2_total_mi", "total_mi"),
    ("var_ms2_bseries_score", "var_bseries_score"),
    ("var_ms2_dotprod_score", "var_dotprod_score"),
    ("var_ms2_intensity_score", "var_intensity_score"),
    (
        "var_ms2_isotope_correlation_score",
        "var_isotope_correlation_score",
    ),
    ("var_ms2_isotope_overlap_score", "var_isotope_overlap_score"),
    ("var_ms2_library_corr", "var_library_corr"),
    ("var_ms2_library_dotprod", "var_library_dotprod"),
    ("var_ms2_library_manhattan", "var_library_manhattan"),
    ("var_ms2_library_rmsd", "var_library_rmsd"),
    (
        "var_ms2_library_rootmeansquare",
        "var_library_rootmeansquare",
    ),
    ("var_ms2_library_sangle", "var_library_sangle"),
    ("var_ms2_log_sn_score", "var_log_sn_score"),
    ("var_ms2_manhattan_score", "var_manhatt_score"),
    ("var_ms2_massdev_score", "var_massdev_score"),
    (
        "var_ms2_massdev_score_weighted",
        "var_massdev_score_weighted",
    ),
    ("var_ms2_mi_score", "var_mi_score"),
    ("var_ms2_mi_weighted_score", "var_mi_weighted_score"),
    ("var_ms2_mi_ratio_score", "var_mi_ratio_score"),
    ("var_ms2_norm_rt_score", "var_norm_rt_score"),
    ("var_ms2_xcorr_coelution", "var_xcorr_coelution"),
    (
        "var_ms2_xcorr_coelution_weighted",
        "var_xcorr_coelution_weighted",
    ),
    ("var_ms2_xcorr_shape", "var_xcorr_shape"),
    ("var_ms2_xcorr_shape_weighted", "var_xcorr_shape_weighted"),
    ("var_ms2_yseries_score", "var_yseries_score"),
    (
        "var_ms2_elution_model_fit_score",
        "var_elution_model_fit_score",
    ),
    ("var_ms2_im_xcorr_shape", "var_im_xcorr_shape"),
    ("var_ms2_im_xcorr_coelution", "var_im_xcorr_coelution"),
    ("var_ms2_im_delta_score", "var_im_delta_score"),
    ("var_ms2_im_log_intensity", "im_log_intensity"),
];

/// Per-transition score columns an MS2 subordinate row leaves null.
const FT_NULL_SCORE_COLUMNS: [&str; 21] = [
    "var_intensity_score",
    "var_intensity_ratio_score",
    "var_log_intensity",
    "var_xcorr_coelution",
    "var_xcorr_shape",
    "var_log_sn_score",
    "var_massdev_score",
    "var_mi_score",
    "var_mi_ratio_score",
    "var_isotope_correlation_score",
    "var_isotope_overlap_score",
    "exp_im",
    "exp_im_leftwidth",
    "exp_im_rightwidth",
    "delta_im",
    "var_im_delta_score",
    "var_im_log_intensity",
    "var_im_xcorr_coelution_contrast",
    "var_im_xcorr_shape_contrast",
    "var_im_xcorr_coelution_combined",
    "var_im_xcorr_shape_combined",
];

/// The peak-shape columns: (column, subordinate meta key, UIS list suffix).
const PEAK_SHAPE_COLUMNS: [&str; 13] = [
    "start_position_at_5",
    "end_position_at_5",
    "start_position_at_10",
    "end_position_at_10",
    "start_position_at_50",
    "end_position_at_50",
    "total_width",
    "tailing_factor",
    "asymmetry_factor",
    "slope_of_baseline",
    "baseline_delta_2_height",
    "points_across_baseline",
    "points_across_half_height",
];

/// The per-transition columns of an identification (UIS) row, in the
/// source's order: (column, list meta key suffix after `id_<target|decoy>_`).
const UIS_COLUMNS: [(&str, &str); 28] = [
    ("area_intensity", "area_intensity"),
    ("total_area_intensity", "total_area_intensity"),
    ("apex_intensity", "apex_intensity"),
    ("apex_rt", "peak_apex_position"),
    ("rt_fwhm", "width_at_50"),
    ("masserror_ppm", "ind_massdev_score"),
    ("total_mi", "total_mi"),
    ("var_intensity_score", "intensity_score"),
    ("var_intensity_ratio_score", "intensity_ratio_score"),
    ("var_log_intensity", "ind_log_intensity"),
    ("var_xcorr_coelution", "ind_xcorr_coelution"),
    ("var_xcorr_shape", "ind_xcorr_shape"),
    ("var_log_sn_score", "ind_log_sn_score"),
    ("var_massdev_score", "ind_massdev_score"),
    ("var_mi_score", "ind_mi_score"),
    ("var_mi_ratio_score", "ind_mi_ratio_score"),
    ("var_isotope_correlation_score", "ind_isotope_correlation"),
    ("var_isotope_overlap_score", "ind_isotope_overlap"),
    ("exp_im", "ind_im_drift"),
    ("exp_im_leftwidth", "ind_im_drift_left"),
    ("exp_im_rightwidth", "ind_im_drift_right"),
    ("delta_im", "ind_im_delta"),
    ("var_im_delta_score", "ind_im_delta_score"),
    ("var_im_log_intensity", "ind_im_log_intensity"),
    (
        "var_im_xcorr_coelution_contrast",
        "ind_im_contrast_coelution",
    ),
    ("var_im_xcorr_shape_contrast", "ind_im_contrast_shape"),
    (
        "var_im_xcorr_coelution_combined",
        "ind_im_sum_contrast_coelution",
    ),
    ("var_im_xcorr_shape_combined", "ind_im_sum_contrast_shape"),
];

/// The canonical id maps of the library.
struct LibraryIds {
    compound_to_precursor: HashMap<String, i64>,
    transition_to_id: HashMap<String, i64>,
    valid_transition_ids: HashSet<i64>,
}

fn canonical_int(text: &str, what: &str) -> Result<i64> {
    to_int64(text).ok_or_else(|| Error::InvalidValue(format!("Non-canonical {what} id '{text}'")))
}

impl LibraryIds {
    fn new(exp: &LightTargetedExperiment) -> Result<Self> {
        let mut compound_to_precursor = HashMap::with_capacity(exp.compounds.len());
        for compound in &exp.compounds {
            let id = canonical_int(&compound.id, "precursor")?;
            compound_to_precursor
                .entry(compound.id.clone())
                .or_insert(id);
        }
        let mut transition_to_id = HashMap::with_capacity(exp.transitions.len());
        let mut valid_transition_ids = HashSet::with_capacity(exp.transitions.len());
        for transition in &exp.transitions {
            let id = canonical_int(&transition.transition_name, "transition")?;
            transition_to_id
                .entry(transition.transition_name.clone())
                .or_insert(id);
            valid_transition_ids.insert(id);
        }
        Ok(Self {
            compound_to_precursor,
            transition_to_id,
            valid_transition_ids,
        })
    }
}

/// The three run tables under construction.
struct RunTables {
    features: TableBuilder,
    feature_transition: TableBuilder,
    feature_precursor: TableBuilder,
}

impl RunTables {
    fn new() -> Self {
        Self {
            features: TableBuilder::new(OSWFeatureSchema::schema()),
            feature_transition: TableBuilder::new(OSWFeatureTransitionSchema::schema()),
            feature_precursor: TableBuilder::new(OSWFeaturePrecursorSchema::schema()),
        }
    }

    /// The feature row, its MS2 transition rows, its MS1 isotope rows and,
    /// with `uis`, its identification transition rows.
    fn add_feature(
        &mut self,
        feature: &Feature,
        run_id: i64,
        ids: &LibraryIds,
        uis: bool,
    ) -> Result<()> {
        let meta = &feature.metadata;
        let feature_id = clear_sign_bit(feature.unique_id);
        if !meta.contains_key("PeptideRef") {
            return Err(Error::MissingInformation(
                "Feature missing PeptideRef meta value".into(),
            ));
        }
        let peptide_ref = meta_text(meta, "PeptideRef");
        let Some(&precursor_id) = ids.compound_to_precursor.get(&peptide_ref) else {
            return Err(Error::InvalidValue(format!(
                "Feature references unknown peptide: {peptide_ref}"
            )));
        };

        let f = &mut self.features;
        f.i64("feature_id", feature_id);
        f.i64("run_id", run_id);
        f.i64("precursor_id", precursor_id);
        f.f64("exp_rt", finite(Some(feature.rt)));
        f.f64("exp_im", feature_score(meta, "im_drift", None));
        f.f64("norm_rt", feature_score(meta, "norm_RT", Some(-1.0)));
        f.f64("delta_rt", feature_score(meta, "delta_rt", Some(-1.0)));
        f.f64("left_width", feature_score(meta, "leftWidth", None));
        f.f64("right_width", feature_score(meta, "rightWidth", None));
        f.f64(
            "exp_im_leftwidth",
            feature_score(meta, "im_drift_left", None),
        );
        f.f64(
            "exp_im_rightwidth",
            feature_score(meta, "im_drift_right", None),
        );
        let has_ms1 = meta.contains_key("var_ms1_ppm_diff");
        for (column, key) in MS1_COLUMNS {
            let value = if has_ms1 {
                finite(meta_double(meta, key))
            } else {
                None
            };
            f.f64(column, value);
        }
        f.f64(
            "ms2_area_intensity",
            finite(Some(f64::from(feature.intensity))),
        );
        for (column, key) in MS2_COLUMNS {
            f.f64(column, finite(meta_double(meta, key)));
        }

        let masserror_ppm = separate_score(meta, "masserror_ppm");
        // Several per-MS2 lists hold entries for MS2 subordinates only.
        let mut ms2_index = 0usize;
        for sub in &feature.subordinates {
            let sub_meta = &sub.metadata;
            if has_level(sub_meta, "MS2") {
                if !sub_meta.contains_key("native_id") {
                    continue;
                }
                let native_id = meta_text(sub_meta, "native_id");
                let transition_id = match ids.transition_to_id.get(&native_id) {
                    Some(&id) => id,
                    None => match to_int64(&native_id) {
                        Some(id) if ids.valid_transition_ids.contains(&id) => id,
                        Some(_) => {
                            return Err(Error::InvalidValue(format!(
                                "Transition references unknown numeric id: {native_id}"
                            )));
                        }
                        None => {
                            return Err(Error::InvalidValue(format!(
                                "Transition references unknown id: {native_id}"
                            )));
                        }
                    },
                };
                let ft = &mut self.feature_transition;
                ft.i64("feature_id", feature_id);
                ft.i64("run_id", run_id);
                ft.i64("transition_id", transition_id);
                ft.f64("area_intensity", finite(Some(f64::from(sub.intensity))));
                ft.f64(
                    "total_area_intensity",
                    finite(meta_double(sub_meta, "total_xic")),
                );
                ft.f64(
                    "apex_intensity",
                    finite(meta_double(sub_meta, "peak_apex_int")),
                );
                ft.f64(
                    "apex_rt",
                    finite(meta_double(sub_meta, "peak_apex_position")),
                );
                ft.f64("rt_fwhm", finite(meta_double(sub_meta, "width_at_50")));
                ft.f64("masserror_ppm", from_list(&masserror_ppm, ms2_index));
                ms2_index += 1;
                ft.f64("total_mi", finite(meta_double(sub_meta, "total_mi")));
                for column in FT_NULL_SCORE_COLUMNS {
                    ft.f64(column, None);
                }
                let has_peak_shape = sub_meta.contains_key("start_position_at_5");
                for column in PEAK_SHAPE_COLUMNS {
                    let value = if has_peak_shape {
                        finite(meta_double(sub_meta, column))
                    } else {
                        None
                    };
                    ft.f64(column, value);
                }
            } else if has_level(sub_meta, "MS1") && sub.intensity > 0.0 {
                let native_id = meta_text(sub_meta, "native_id");
                // `StringUtils::split(s, "Precursor_i", parts)`; an empty
                // string gives no parts.
                let parts: Vec<&str> = if native_id.is_empty() {
                    Vec::new()
                } else {
                    native_id.split("Precursor_i").collect()
                };
                if parts.len() < 2 {
                    continue;
                }
                let Some(isotope) = to_int64(parts[1]) else {
                    continue;
                };
                let fp = &mut self.feature_precursor;
                fp.i64("feature_id", feature_id);
                fp.i64("run_id", run_id);
                // `static_cast<int32_t>`, which wraps.
                fp.i32("precursor_isotope", isotope as i32);
                fp.f64(
                    "precursor_area_intensity",
                    finite(Some(f64::from(sub.intensity))),
                );
                // The subordinate's own value; the Release build writes a stale
                // one here (argument evaluation order, see the module docs).
                fp.f64(
                    "precursor_apex_intensity",
                    finite(meta_double(sub_meta, "peak_apex_int")),
                );
            }
        }

        if uis {
            for kind in ["target", "decoy"] {
                self.add_identification_rows(meta, kind, feature_id, run_id, ids);
            }
        }
        Ok(())
    }

    /// The identification-transition rows of one class (`target` or `decoy`).
    fn add_identification_rows(
        &mut self,
        meta: &MetaInfo,
        kind: &str,
        feature_id: i64,
        run_id: i64,
        ids: &LibraryIds,
    ) {
        let list = |suffix: &str| separate_score(meta, &format!("id_{kind}_{suffix}"));
        let names = list("transition_names");
        let values: Vec<(&str, Vec<String>)> = UIS_COLUMNS
            .iter()
            .map(|(column, suffix)| (*column, list(suffix)))
            .collect();
        let shape: Vec<(&str, Vec<String>)> = PEAK_SHAPE_COLUMNS
            .iter()
            .map(|column| (*column, list(&format!("ind_{column}"))))
            .collect();
        let peak_shape = shape[0].1.first().is_some_and(|first| first != "0");
        for (i, name) in names.iter().enumerate() {
            let Some(&transition_id) = ids.transition_to_id.get(name) else {
                continue;
            };
            let ft = &mut self.feature_transition;
            ft.i64("feature_id", feature_id);
            ft.i64("run_id", run_id);
            ft.i64("transition_id", transition_id);
            for (column, list) in &values {
                ft.f64(column, from_list(list, i));
            }
            for (column, list) in &shape {
                let value = if peak_shape { from_list(list, i) } else { None };
                ft.f64(column, value);
            }
        }
    }
}

/// Every regular file below `base`, as (`/`-separated relative name, path),
/// sorted by name.
fn files_below(base: &Path) -> Result<Vec<(String, PathBuf)>> {
    let mut out = Vec::new();
    let mut pending = vec![(base.to_path_buf(), String::new())];
    while let Some((dir, prefix)) = pending.pop() {
        for entry in std::fs::read_dir(&dir)? {
            let entry = entry?;
            let name = entry.file_name().to_string_lossy().into_owned();
            let relative = if prefix.is_empty() {
                name
            } else {
                format!("{prefix}/{name}")
            };
            let path = entry.path();
            if path.is_dir() {
                pending.push((path, relative));
            } else {
                out.push((relative, path));
            }
        }
    }
    out.sort();
    Ok(out)
}

/// The runs of `runs.parquet` (the source's `readRuns_`), empty when the
/// file does not exist.
fn read_runs(runs_parquet: &Path) -> Result<Vec<(i64, String)>> {
    if !runs_parquet.exists() {
        return Ok(Vec::new());
    }
    let table = ParquetFile::read_table(runs_parquet)?;
    let run_ids = ParquetFile::get_column(&table, OSWRunSchema::RUN_ID)?;
    let filenames = ParquetFile::get_column(&table, OSWRunSchema::FILENAME)?;
    (0..table.num_rows())
        .map(|row| {
            Ok((
                ParquetFile::get_int64(Some(run_ids.as_ref()), row, 0, false)?,
                ParquetFile::get_string(Some(filenames.as_ref()), row)?,
            ))
        })
        .collect()
}

/// The container `metadata.json` (the source's `writeMetadata_`), byte for
/// byte except `build_time`.
fn container_metadata(runs: &[(i64, String)], counts: &[RunCounts], total: RunCounts) -> String {
    let version = ParquetFile::json_escape(crate::CORE_SDK_VERSION);
    let mut out = String::new();
    out.push_str("{\n  \"openms\": {\n    \"schema_version\": 1,\n");
    out.push_str("    \"generator\": \"OpenSwathOSWParquetWriter\",\n");
    out.push_str(&format!("    \"openms_version\": \"{version}\",\n"));
    out.push_str("    \"build_time\": \"\",\n");
    out.push_str(&format!(
        "    \"tool\": {{\"name\": \"OpenSwathWorkflow\", \"version\": \"{version}\"}},\n"
    ));
    out.push_str("    \"counts\": {\n");
    out.push_str(&format!("      \"runs\": {},\n", runs.len()));
    out.push_str(&format!("      \"features\": {},\n", total.features));
    out.push_str(&format!(
        "      \"feature_precursor\": {},\n",
        total.feature_precursor
    ));
    out.push_str(&format!(
        "      \"feature_transition\": {}\n",
        total.feature_transition
    ));
    out.push_str("    }\n  },\n  \"runs\": [\n");
    for (i, ((run_id, filename), c)) in runs.iter().zip(counts).enumerate() {
        out.push_str(&format!(
            "    {{\"id\": {run_id}, \"filename\": \"{}\", \"counts\": {{\"features\": {}, \"feature_precursor\": {}, \"feature_transition\": {}}}}}",
            ParquetFile::json_escape(filename),
            c.features,
            c.feature_precursor,
            c.feature_transition
        ));
        if i + 1 < runs.len() {
            out.push(',');
        }
        out.push('\n');
    }
    out.push_str("  ]\n}\n");
    out
}

impl OpenSwathOSWParquetWriter {
    /// A writer that replaces an existing archive (the source's default).
    pub fn new() -> Self {
        Self::default()
    }

    /// Whether an existing archive at the output path is extended (`true`)
    /// or replaced (`false`, the default). Ports `setPreserveExisting`.
    ///
    /// An existing *directory* is always extended.
    pub fn set_preserve_existing(&mut self, preserve: bool) {
        self.preserve_existing = preserve;
    }

    /// Whether an existing archive is extended; see
    /// [`set_preserve_existing`](Self::set_preserve_existing). Native getter.
    pub fn preserve_existing(&self) -> bool {
        self.preserve_existing
    }

    /// Write the library and one run's features.
    ///
    /// Ports the six-argument `write`; see
    /// [`write_with_source_ids`](Self::write_with_source_ids).
    ///
    /// # Errors
    ///
    /// As [`write_with_source_ids`](Self::write_with_source_ids).
    pub fn write(
        &self,
        output_path: impl AsRef<Path>,
        assay_library: &LightTargetedExperiment,
        feature_map: &FeatureMap,
        run_id: u64,
        input_filename: &str,
        enable_uis_scoring: bool,
    ) -> Result<Vec<String>> {
        self.write_with_source_ids(
            output_path,
            assay_library,
            feature_map,
            run_id,
            input_filename,
            enable_uis_scoring,
            None,
        )
    }

    /// Write the library and one run's features, with optional source-id
    /// provenance for the library's `traml_id` columns.
    ///
    /// Ports the seven-argument `write`. `run_id` has its top bit cleared
    /// (`SqliteHelper::clearSignBit`), as has every feature's unique id. When
    /// `output_path` is an existing directory the container is written (or
    /// extended) there; otherwise it is assembled in a temporary directory,
    /// or, when [`preserve_existing`](Self::preserve_existing) is set and an
    /// archive exists, in the unpacked archive, and added entry by entry to
    /// the archive at `output_path`, followed by its size index.
    ///
    /// The library is written only when the container has none; an existing
    /// one must hold exactly the precursor ids and the transition to
    /// precursor references of `assay_library`.
    ///
    /// Returns the warnings the source writes to `OPENMS_LOG_WARN` (one, when
    /// an existing archive is extended).
    ///
    /// # Errors
    ///
    /// - [`OpenSwathLibraryIDNormalizer::validate_canonical_ids`]'s error when
    ///   the library ids are not canonical.
    /// - [`Error::InvalidValue`] "Existing OSWPQ library contains duplicate
    ///   precursor_id" / "… transition_id", "Existing OSWPQ library uses a
    ///   different canonical precursor/transition ID domain. …", "Run output
    ///   directory already exists", "Feature references unknown peptide",
    ///   "Transition references unknown numeric id" / "… unknown id", "Run id
    ///   already present in runs.parquet" and "Failed to append runs table",
    ///   as the source.
    /// - [`Error::MissingInformation`] "Feature missing PeptideRef meta
    ///   value", as the source.
    /// - The errors of [`TransitionParquetFile`], [`ParquetFile`] and
    ///   [`ZipArchiveFile`], and [`Error::Io`] for a file that cannot be
    ///   written.
    ///
    /// [`TransitionParquetFile`]: crate::analysis::openswath::transition_parquet_file::TransitionParquetFile
    /// [`ParquetFile`]: crate::format::parquet_file::ParquetFile
    /// [`ZipArchiveFile`]: crate::format::zip_archive_file::ZipArchiveFile
    #[allow(clippy::too_many_arguments)]
    pub fn write_with_source_ids(
        &self,
        output_path: impl AsRef<Path>,
        assay_library: &LightTargetedExperiment,
        feature_map: &FeatureMap,
        run_id: u64,
        input_filename: &str,
        enable_uis_scoring: bool,
        source_ids: Option<&SourceIDMapping>,
    ) -> Result<Vec<String>> {
        OpenSwathLibraryIDNormalizer::validate_canonical_ids(assay_library)?;
        let output_path = output_path.as_ref();
        let mut warnings = Vec::new();
        let run_id_clean = clear_sign_bit(run_id);
        let output_is_dir = file::is_directory(output_path);
        let mut temp_dir: Option<TempDir> = None;
        let base_dir: PathBuf = if output_is_dir {
            output_path.to_path_buf()
        } else if file::exists(output_path) && self.preserve_existing {
            let unpacked = ZipArchiveFile::unzip_directory(output_path, &mut temp_dir)?;
            warnings.push(format!(
                "Existing archive '{}' was detected — its contents will be preserved and new runs will be appended.\nTo overwrite the archive remove it before running the writer.",
                file::absolute_path(output_path)?.display()
            ));
            unpacked
        } else {
            let dir = TempDir::new(false)?;
            let base = dir.path().join("oswpq_output");
            temp_dir = Some(dir);
            base
        };

        let library_dir = base_dir.join("library");
        let runs_dir = base_dir.join("runs");
        if !file::exists(&base_dir) {
            file::make_dir(&base_dir)?;
        }
        if !file::exists(&runs_dir) {
            file::make_dir(&runs_dir)?;
        }

        let precursors_path = library_dir.join("precursors.parquet");
        let transitions_path = library_dir.join("transitions.parquet");
        let library_ready = file::exists(&precursors_path) && file::exists(&transitions_path);
        if library_ready {
            check_library_domain(
                assay_library,
                &precursors_path,
                &transitions_path,
                &library_dir,
            )?;
        } else {
            let library_tmp_dir = base_dir.join("library_tmp");
            if file::exists(&library_tmp_dir) {
                file::remove_dir_recursively(&library_tmp_dir)?;
            }
            if file::exists(&library_dir) {
                file::remove_dir_recursively(&library_dir)?;
            }
            file::make_dir(&library_tmp_dir)?;
            file::make_dir(&library_dir)?;
            TransitionParquetFile::new()
                .convert_light_targeted_experiment_to_parquet_with_source_ids(
                    &library_tmp_dir,
                    assay_library,
                    source_ids,
                )?;
            file::copy_dir_recursively(
                library_tmp_dir.join("library"),
                &library_dir,
                CopyOptions::Overwrite,
            )?;
            file::remove_dir_recursively(&library_tmp_dir)?;
        }

        let run_path = runs_dir.join(format!("run_id={run_id_clean}"));
        if file::exists(&run_path) {
            return Err(Error::InvalidValue(format!(
                "Run output directory already exists: {}",
                run_path.display()
            )));
        }
        file::make_dir(&run_path)?;

        // The library has crossed the canonical id boundary: run-level
        // foreign keys use its ids exactly.
        let ids = LibraryIds::new(assay_library)?;
        let mut tables = RunTables::new();
        for feature in &feature_map.features {
            tables.add_feature(feature, run_id_clean, &ids, enable_uis_scoring)?;
        }

        let written = (|| -> Result<()> {
            let features = tables.features.finish("Features table")?;
            let feature_transition = tables
                .feature_transition
                .finish("Feature transition table")?;
            let feature_precursor = tables.feature_precursor.finish("Feature precursor table")?;
            ParquetFile::write_table(&features, run_path.join("features.parquet"))?;
            ParquetFile::write_table(
                &feature_transition,
                run_path.join("feature_transition.parquet"),
            )?;
            ParquetFile::write_table(
                &feature_precursor,
                run_path.join("feature_precursor.parquet"),
            )
        })();
        if let Err(error) = written {
            if file::exists(&run_path) {
                // The source removes the partial run so a retry is possible.
                let _ = file::remove_dir_recursively(&run_path);
            }
            return Err(error);
        }

        let runs_parquet = runs_dir.join("runs.parquet");
        append_run(&runs_parquet, run_id_clean, input_filename)?;

        let runs = read_runs(&runs_parquet)?;
        let mut run_counts = Vec::with_capacity(runs.len());
        let mut total = RunCounts::default();
        for (id, _) in &runs {
            let path = runs_dir.join(format!("run_id={id}"));
            let counts = RunCounts {
                features: ParquetFile::row_count(path.join("features.parquet"))?,
                feature_precursor: ParquetFile::row_count(path.join("feature_precursor.parquet"))?,
                feature_transition: ParquetFile::row_count(
                    path.join("feature_transition.parquet"),
                )?,
            };
            total.features += counts.features;
            total.feature_precursor += counts.feature_precursor;
            total.feature_transition += counts.feature_transition;
            run_counts.push(counts);
        }
        std::fs::write(
            base_dir.join("metadata.json"),
            container_metadata(&runs, &run_counts, total),
        )?;

        if !output_is_dir {
            let output_zip = file::absolute_path(output_path)?;
            if file::exists(&output_zip) && !self.preserve_existing {
                file::remove(&output_zip)?;
            }
            for (name, path) in files_below(&base_dir)? {
                ZipArchiveFile::add_or_replace_from_file(output_path, &name, &path)?;
            }
            ZipArchiveFile::write_sidecar_index(&output_zip)?;
        }
        // The source's `TempDir` destructor removes the staging tree and
        // ignores a failure to do so; dropping the guard does the same.
        drop(temp_dir);
        Ok(warnings)
    }
}

/// The existing library must describe the same canonical id domain.
fn check_library_domain(
    exp: &LightTargetedExperiment,
    precursors_path: &Path,
    transitions_path: &Path,
    library_dir: &Path,
) -> Result<()> {
    let mut expected_precursors = HashSet::with_capacity(exp.compounds.len());
    for compound in &exp.compounds {
        expected_precursors.insert(canonical_int(&compound.id, "precursor")?);
    }
    let mut expected_refs: HashMap<i64, i64> = HashMap::with_capacity(exp.transitions.len());
    for transition in &exp.transitions {
        let id = canonical_int(&transition.transition_name, "transition")?;
        let reference = canonical_int(&transition.peptide_ref, "precursor")?;
        expected_refs.entry(id).or_insert(reference);
    }

    let table = ParquetFile::read_table(precursors_path)?;
    let column = ParquetFile::get_column(&table, OSWPrecursorSchema::PRECURSOR_ID)?;
    let mut existing_precursors = HashSet::with_capacity(table.num_rows());
    for row in 0..table.num_rows() {
        let id = ParquetFile::get_int64(Some(column.as_ref()), row, -1, false)?;
        if !existing_precursors.insert(id) {
            return Err(Error::InvalidValue(format!(
                "Existing OSWPQ library contains duplicate precursor_id: {id}"
            )));
        }
    }
    let table = ParquetFile::read_table(transitions_path)?;
    let ids = ParquetFile::get_column(&table, OSWTransitionSchema::TRANSITION_ID)?;
    let refs = ParquetFile::get_column(&table, OSWTransitionSchema::PRECURSOR_ID)?;
    let mut existing_refs: HashMap<i64, i64> = HashMap::with_capacity(table.num_rows());
    for row in 0..table.num_rows() {
        let id = ParquetFile::get_int64(Some(ids.as_ref()), row, -1, false)?;
        let reference = ParquetFile::get_int64(Some(refs.as_ref()), row, -1, false)?;
        if existing_refs.insert(id, reference).is_some() {
            return Err(Error::InvalidValue(format!(
                "Existing OSWPQ library contains duplicate transition_id: {id}"
            )));
        }
    }
    if existing_precursors != expected_precursors || existing_refs != expected_refs {
        return Err(Error::InvalidValue(format!(
            "Existing OSWPQ library uses a different canonical precursor/transition ID domain. Please rebuild the archive or use a different output path.: {}",
            library_dir.display()
        )));
    }
    Ok(())
}

/// Add one run to `runs.parquet`, creating it when missing.
fn append_run(runs_parquet: &Path, run_id: i64, input_filename: &str) -> Result<()> {
    let mut run_ids = Vec::new();
    let mut filenames = Vec::new();
    if runs_parquet.exists() {
        let existing = ParquetFile::read_table(runs_parquet)?;
        let expected = OSWRunSchema::schema();
        // `arrow::ConcatenateTables` refuses tables whose schemas differ.
        let same_schema = existing.schema().fields().len() == expected.fields().len()
            && existing
                .schema()
                .fields()
                .iter()
                .zip(expected.fields().iter())
                .all(|(a, e)| a.name() == e.name() && a.data_type() == e.data_type());
        let column = ParquetFile::get_column(&existing, OSWRunSchema::RUN_ID)?;
        for row in 0..existing.num_rows() {
            if ParquetFile::get_int64(Some(column.as_ref()), row, 0, false)? == run_id {
                return Err(Error::InvalidValue(format!(
                    "Run id already present in runs.parquet: {run_id}"
                )));
            }
        }
        if !same_schema {
            return Err(Error::InvalidValue(
                "Failed to append runs table: schemas differ".into(),
            ));
        }
        for (id, name) in read_runs(runs_parquet)? {
            run_ids.push(id);
            filenames.push(name);
        }
    }
    run_ids.push(run_id);
    filenames.push(input_filename.to_owned());
    let table = RecordBatch::try_new(
        OSWRunSchema::schema(),
        vec![
            Arc::new(Int64Array::from(run_ids)),
            Arc::new(StringArray::from(filenames)),
        ],
    )
    .map_err(|e| Error::InvalidValue(format!("Runs table: {e}")))?;
    let result = validation::validate(&table, OSWRunSchema::schema().as_ref(), Mode::Strict);
    if !result.valid {
        return Err(Error::InvalidValue(format!(
            "Runs table schema validation failed: {result}"
        )));
    }
    ParquetFile::write_table(&table, runs_parquet)
}
