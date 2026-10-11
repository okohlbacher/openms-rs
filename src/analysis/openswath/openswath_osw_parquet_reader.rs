// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Read OpenSwathWorkflow results from an OSW Parquet container (`.oswpq`).
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathOSWParquetReader.h` and its
//! translation unit. The reader reads the layout
//! [`OpenSwathOSWParquetWriter`](crate::analysis::openswath::openswath_osw_parquet_writer::OpenSwathOSWParquetWriter)
//! produces (`library/precursors.parquet`, `library/transitions.parquet`,
//! `runs/runs.parquet` and the per-run tables), as a directory or a ZIP
//! archive, and returns flat, column-oriented tables for downstream scoring
//! (PyProphet) or single-table exports:
//!
//! - [`OpenSwathOSWParquetReader::load`](crate::analysis::openswath::openswath_osw_parquet_reader::OpenSwathOSWParquetReader::load) and [`OpenSwathOSWParquetReader::rows`](crate::analysis::openswath::openswath_osw_parquet_reader::OpenSwathOSWParquetReader::rows):
//!   one [`Row`](crate::analysis::openswath::openswath_osw_parquet_reader::Row) per feature, with its precursor's charge and decoy flag and
//!   the number of detecting library transitions;
//! - [`OpenSwathOSWParquetReader::fetch_peak_group_features`](crate::analysis::openswath::openswath_osw_parquet_reader::OpenSwathOSWParquetReader::fetch_peak_group_features): the feature
//!   scores (`var_ms2_*`, `ms2_*`, optionally `var_ms1_*`, `ms1_*`), sorted by
//!   run, precursor and retention time;
//! - [`OpenSwathOSWParquetReader::fetch_transition_features`](crate::analysis::openswath::openswath_osw_parquet_reader::OpenSwathOSWParquetReader::fetch_transition_features): the
//!   per-transition scores, sorted by run, precursor, retention time and
//!   transition;
//! - [`OpenSwathOSWParquetReader::fetch_unscored_data`](crate::analysis::openswath::openswath_osw_parquet_reader::OpenSwathOSWParquetReader::fetch_unscored_data): the PyProphet
//!   "unscored" table, in file order.
//!
//! Only the columns a function reads must exist; a score column missing from
//! one run reads as NaN there, a missing optional value as the documented
//! default.
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_PARQUET_SUPPORT.md).
//!
//! # Differences from the source
//!
//! - Both sorts reproduce the Release build's `std::sort` comparison by
//!   comparison ([`source_sort_permutation`](crate::math::source_sort::source_sort_permutation)),
//!   so rows that compare equal (the same run, precursor and retention time)
//!   keep the order the source leaves them in.
//! - Archive entries are read by extracting them to a temporary file, the
//!   source's fallback path; its first attempt (`ZipRandomAccessFile`) reads
//!   the same bytes and is not ported.
//! - [`OpenSwathOSWParquetReader::load`](crate::analysis::openswath::openswath_osw_parquet_reader::OpenSwathOSWParquetReader::load) leaves [`rows`](crate::analysis::openswath::openswath_osw_parquet_reader::OpenSwathOSWParquetReader::rows)
//!   empty when it fails; the source leaves the rows read before the error.
//! - The warning the source writes to `OPENMS_LOG_WARN` for a `main_score`
//!   that is not among the score columns is returned in
//!   [`PeakGroupFeatureScoresResult::warnings`](crate::analysis::openswath::openswath_osw_parquet_reader::PeakGroupFeatureScoresResult::warnings).
//! - The C++ result members `RT`, `assay_RT`, `delta_RT`, `Charge`,
//!   `Intensity`, `aggr_prec_Peak_Area`, `aggr_prec_Peak_Apex`, `leftWidth`,
//!   `rightWidth`, `EXP_IM`, `IM_leftWidth` and `IM_rightWidth` of
//!   `UnscoredResult` are snake case here; see [`UnscoredResult`](crate::analysis::openswath::openswath_osw_parquet_reader::UnscoredResult).

use crate::format::container_io::{self as file, TempDir};
use crate::format::parquet_file::ParquetFile;
use crate::format::zip_archive_file::ZipArchiveFile;
use crate::math::source_sort::source_sort_permutation;
use crate::{Error, Result};
use arrow_array::{Array, ArrayRef, RecordBatch};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

/// One feature combined with its precursor and run metadata (the source's
/// `Row`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Row {
    /// Feature id.
    pub feature_id: i64,
    /// Run id, as the features table stores it.
    pub run_id: i64,
    /// Precursor id.
    pub precursor_id: i64,
    /// Experimental retention time; 0 when null.
    pub exp_rt: f64,
    /// MS2 area intensity; 0 when null or absent.
    pub ms2_area_intensity: f64,
    /// MS2 total area intensity; 0 when null or absent.
    pub ms2_total_area_intensity: f64,
    /// MS2 apex intensity; 0 when null or absent.
    pub ms2_apex_intensity: f64,
    /// Precursor charge from the library.
    pub precursor_charge: i32,
    /// Precursor decoy flag from the library.
    pub decoy: bool,
    /// Number of detecting library transitions of the precursor.
    pub transition_count: i64,
    /// `<run_id>_<precursor_id>_<feature_id>`.
    pub group_id: String,
}

/// Feature-level scores across all runs, column-oriented (the source's
/// `PeakGroupFeatureScoresResult`). Every vector has one entry per feature.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PeakGroupFeatureScoresResult {
    /// Feature ids.
    pub feature_id: Vec<i64>,
    /// Run ids, from `runs.parquet`.
    pub run_id: Vec<i64>,
    /// Precursor ids.
    pub precursor_id: Vec<i64>,
    /// Experimental retention times; 0 when null.
    pub exp_rt: Vec<f64>,
    /// Precursor charges from the library.
    pub precursor_charge: Vec<i32>,
    /// Precursor decoy flags from the library.
    pub decoy: Vec<bool>,
    /// Detecting library transitions per precursor.
    pub transition_count: Vec<i64>,
    /// `<run_id>_<precursor_id>_<feature_id>`.
    pub group_id: Vec<String>,
    /// The `var_ms2_*` and `ms2_*` columns found in any run, sorted, the
    /// requested main score first.
    pub ms2_columns: Vec<String>,
    /// One vector per entry of [`Self::ms2_columns`]; NaN where null or absent.
    pub ms2_values: Vec<Vec<f64>>,
    /// The `var_ms1_*` and `ms1_*` columns found in any run, sorted; only for
    /// level `ms1ms2`.
    pub ms1_columns: Vec<String>,
    /// One vector per entry of [`Self::ms1_columns`].
    pub ms1_values: Vec<Vec<f64>>,
    /// Native: the lines the source writes to `OPENMS_LOG_WARN`.
    pub warnings: Vec<String>,
}

/// Transition-level scores across all runs, column-oriented (the source's
/// `TransitionFeaturesResult`). Every vector has one entry per
/// feature-transition row.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TransitionFeaturesResult {
    /// Feature ids.
    pub feature_id: Vec<i64>,
    /// Run ids, from `runs.parquet`.
    pub run_id: Vec<i64>,
    /// Precursor ids of the features.
    pub precursor_id: Vec<i64>,
    /// Retention times of the features; 0 when null.
    pub exp_rt: Vec<f64>,
    /// Precursor charges from the library.
    pub precursor_charge: Vec<i32>,
    /// Transition ids.
    pub transition_id: Vec<i64>,
    /// Product charges from the library.
    pub product_charge: Vec<i32>,
    /// Transition decoy flags from the library.
    pub decoy: Vec<bool>,
    /// Transition peak area; NaN when null or absent.
    pub area_intensity: Vec<f64>,
    /// Total transition area; NaN when null or absent.
    pub total_area_intensity: Vec<f64>,
    /// Apex intensity; NaN when null or absent.
    pub apex_intensity: Vec<f64>,
    /// Apex retention time; NaN when null or absent.
    pub apex_rt: Vec<f64>,
    /// Full width at half maximum; NaN when null or absent.
    pub rt_fwhm: Vec<f64>,
    /// Mass error in ppm; NaN when null or absent.
    pub masserror_ppm: Vec<f64>,
    /// Total mutual information; NaN when null or absent.
    pub total_mi: Vec<f64>,
    /// The `var_*` columns found in any run's `feature_transition.parquet`, sorted.
    pub transition_var_columns: Vec<String>,
    /// One vector per entry of [`Self::transition_var_columns`].
    pub transition_var_values: Vec<Vec<f64>>,
    /// `<run_id>_<feature_id>_<precursor_id>_<transition_id>`.
    pub group_id: Vec<String>,
}

/// The PyProphet "unscored" feature table, column-oriented (the source's
/// `UnscoredResult`). Every vector has one entry per feature, in file order;
/// a value the files do not provide is NaN, 0, `false` or empty.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct UnscoredResult {
    /// Run id (`id_run`).
    pub id_run: Vec<i64>,
    /// Peptide id from `library/precursor_peptide_mapping.parquet`; 0 when unknown (`id_peptide`).
    pub id_peptide: Vec<i64>,
    /// Precursor id (`transition_group_id`).
    pub transition_group_id: Vec<i64>,
    /// Precursor decoy flag (`decoy`).
    pub decoy: Vec<bool>,
    /// Run id (`run_id`).
    pub run_id: Vec<i64>,
    /// Run file name (`filename`).
    pub filename: Vec<String>,
    /// Feature retention time `FEATURE.EXP_RT` (source `RT`).
    pub rt: Vec<f64>,
    /// `EXP_RT - DELTA_RT` (source `assay_rt`).
    pub assay_rt: Vec<f64>,
    /// `FEATURE.DELTA_RT` (source `delta_rt`).
    pub delta_rt: Vec<f64>,
    /// Library retention time `PRECURSOR.LIBRARY_RT` (source `assay_RT`).
    pub assay_rt_library: Vec<f64>,
    /// `FEATURE.NORM_RT - PRECURSOR.LIBRARY_RT` (source `delta_RT`).
    pub delta_rt_library: Vec<f64>,
    /// Feature id `FEATURE.ID` (source `id`).
    pub id: Vec<i64>,
    /// Precursor charge (source `Charge`).
    pub charge: Vec<i32>,
    /// Precursor m/z (source `mz`).
    pub mz: Vec<f64>,
    /// MS2 area intensity (source `Intensity`).
    pub intensity: Vec<f64>,
    /// MS1 area intensity (source `aggr_prec_Peak_Area`).
    pub aggr_prec_peak_area: Vec<f64>,
    /// MS1 apex intensity (source `aggr_prec_Peak_Apex`).
    pub aggr_prec_peak_apex: Vec<f64>,
    /// Left peak boundary (source `leftWidth`).
    pub left_width: Vec<f64>,
    /// Right peak boundary (source `rightWidth`).
    pub right_width: Vec<f64>,
    /// Ion mobility (source `EXP_IM`).
    pub exp_im: Vec<f64>,
    /// Left ion-mobility boundary (source `IM_leftWidth`).
    pub im_left_width: Vec<f64>,
    /// Right ion-mobility boundary (source `IM_rightWidth`).
    pub im_right_width: Vec<f64>,
    /// The `var_ms2_*` and `ms2_*` columns found in any run, sorted.
    pub ms2_columns: Vec<String>,
    /// One vector per entry of [`Self::ms2_columns`].
    pub ms2_values: Vec<Vec<f64>>,
    /// The `var_ms1_*` and `ms1_*` columns found in any run, sorted.
    pub ms1_columns: Vec<String>,
    /// One vector per entry of [`Self::ms1_columns`].
    pub ms1_values: Vec<Vec<f64>>,
}

/// Reader of OSW Parquet containers.
///
/// Ports `OpenSwathOSWParquetReader`: [`load`](Self::load) fills
/// [`rows`](Self::rows) and remembers the path; the three `fetch_*`
/// functions read a container of their own argument and keep no state.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OpenSwathOSWParquetReader {
    rows: Vec<Row>,
    oswpq_dir: PathBuf,
}

/// Tables of one container, extracted on demand (the source's
/// `open_table_from_entry`).
struct Container<'a> {
    path: &'a Path,
    temp_dir: Option<TempDir>,
}

impl<'a> Container<'a> {
    fn new(path: &'a Path) -> Self {
        Self {
            path,
            temp_dir: None,
        }
    }

    fn table(&mut self, entry: &str) -> Result<RecordBatch> {
        let path =
            ZipArchiveFile::extract_entry_to_temp_file(self.path, entry, &mut self.temp_dir)?;
        ParquetFile::read_table(path)
    }
}

fn features_entry(run_id: i64) -> String {
    format!("runs/run_id={run_id}/features.parquet")
}

fn feature_transition_entry(run_id: i64) -> String {
    format!("runs/run_id={run_id}/feature_transition.parquet")
}

fn col(array: &Option<ArrayRef>) -> Option<&dyn Array> {
    array.as_deref()
}

fn req(array: &ArrayRef) -> Option<&dyn Array> {
    Some(array.as_ref())
}

/// The library precursors: id to (charge, decoy), the first row of an id
/// winning.
fn precursor_info(container: &mut Container<'_>) -> Result<HashMap<i64, (i32, bool)>> {
    let table = container.table("library/precursors.parquet")?;
    let ids = ParquetFile::get_column(&table, "precursor_id")?;
    let charges = ParquetFile::get_column(&table, "charge")?;
    let decoys = ParquetFile::get_optional_column(&table, "decoy");
    let mut info = HashMap::with_capacity(table.num_rows());
    for row in 0..table.num_rows() {
        let id = ParquetFile::get_int64(req(&ids), row, 0, false)?;
        // `static_cast<int>`, which wraps.
        let charge = ParquetFile::get_int64(req(&charges), row, 0, false)? as i32;
        let decoy = ParquetFile::get_bool(col(&decoys), row, false, true)?;
        info.entry(id).or_insert((charge, decoy));
    }
    Ok(info)
}

/// Detecting library transitions per precursor.
fn transition_counts(container: &mut Container<'_>) -> Result<HashMap<i64, i64>> {
    let table = container.table("library/transitions.parquet")?;
    let precursors = ParquetFile::get_column(&table, "precursor_id")?;
    let detecting = ParquetFile::get_column(&table, "detecting")?;
    let mut counts: HashMap<i64, i64> = HashMap::with_capacity(table.num_rows());
    for row in 0..table.num_rows() {
        let id = ParquetFile::get_int64(req(&precursors), row, 0, false)?;
        if ParquetFile::get_bool(req(&detecting), row, true, true)? {
            *counts.entry(id).or_insert(0) += 1;
        }
    }
    Ok(counts)
}

/// The run ids of `runs/runs.parquet`, in file order.
fn run_ids(runs: &RecordBatch) -> Result<Vec<i64>> {
    let ids = ParquetFile::get_column(runs, "run_id")?;
    (0..runs.num_rows())
        .map(|row| ParquetFile::get_int64(req(&ids), row, 0, false))
        .collect()
}

/// Add every column name of `table` accepted by `accept` that is not yet in
/// `names`, keeping first-seen order.
fn discover(table: &RecordBatch, names: &mut Vec<String>, accept: impl Fn(&str) -> bool) {
    for field in table.schema().fields() {
        let name = field.name();
        if accept(name) && !names.iter().any(|n| n == name) {
            names.push(name.clone());
        }
    }
}

fn is_ms2(name: &str) -> bool {
    name.starts_with("var_ms2_") || name.starts_with("ms2_")
}

fn is_ms1(name: &str) -> bool {
    name.starts_with("var_ms1_") || name.starts_with("ms1_")
}

/// A nullable double defaulting to NaN.
fn nan_or(array: Option<&dyn Array>, row: usize) -> Result<f64> {
    ParquetFile::get_double(array, row, f64::NAN, true)
}

/// Apply a permutation (position `k` takes element `order[k]`).
fn permute<T: Clone>(values: &mut Vec<T>, order: &[usize]) {
    let permuted: Vec<T> = order.iter().map(|&i| values[i].clone()).collect();
    *values = permuted;
}

impl OpenSwathOSWParquetReader {
    /// An empty reader (the source's default constructor).
    pub fn new() -> Self {
        Self::default()
    }

    /// A reader that has [`load`](Self::load)ed `oswpq_dir` (the source's
    /// convenience constructor).
    ///
    /// # Errors
    ///
    /// As [`load`](Self::load).
    pub fn from_path(oswpq_dir: impl AsRef<Path>) -> Result<Self> {
        let mut reader = Self::new();
        reader.load(oswpq_dir)?;
        Ok(reader)
    }

    /// Read one [`Row`] per feature of every run of the container at
    /// `oswpq_dir` (directory or archive), in run and file order, replacing
    /// the previous rows, and remember the path.
    ///
    /// Ports `load`. Each row carries its precursor's library charge and
    /// decoy flag (the first library row of a precursor id wins) and the
    /// number of its detecting library transitions; `exp_rt` and the three
    /// MS2 intensities default to 0.
    ///
    /// # Errors
    ///
    /// - [`Error::Io`] of kind `NotFound` for a missing table (the source's
    ///   `Exception::FileNotFound`); [`Error::InvalidValue`] for an archive
    ///   or table that cannot be read.
    /// - [`Error::MissingInformation`] for a missing required column or a
    ///   null required value, and "Unknown precursor_id=\<id\> for
    ///   run_id=\<id\>" for a feature whose precursor the library lacks, as
    ///   the source.
    pub fn load(&mut self, oswpq_dir: impl AsRef<Path>) -> Result<()> {
        let path = oswpq_dir.as_ref();
        self.rows.clear();
        self.oswpq_dir = path.to_path_buf();
        let mut container = Container::new(path);
        let precursors = precursor_info(&mut container)?;
        let counts = transition_counts(&mut container)?;
        let runs = container.table("runs/runs.parquet")?;
        let mut rows = Vec::new();
        for run_id in run_ids(&runs)? {
            let table = container.table(&features_entry(run_id))?;
            let feature_ids = ParquetFile::get_column(&table, "feature_id")?;
            let feature_runs = ParquetFile::get_column(&table, "run_id")?;
            let precursor_ids = ParquetFile::get_column(&table, "precursor_id")?;
            let exp_rt = ParquetFile::get_column(&table, "exp_rt")?;
            let ms2_area = ParquetFile::get_optional_column(&table, "ms2_area_intensity");
            let ms2_total = ParquetFile::get_optional_column(&table, "ms2_total_area_intensity");
            let ms2_apex = ParquetFile::get_optional_column(&table, "ms2_apex_intensity");
            rows.reserve(table.num_rows());
            for row in 0..table.num_rows() {
                let mut out = Row {
                    feature_id: ParquetFile::get_int64(req(&feature_ids), row, 0, false)?,
                    run_id: ParquetFile::get_int64(req(&feature_runs), row, 0, false)?,
                    precursor_id: ParquetFile::get_int64(req(&precursor_ids), row, 0, false)?,
                    exp_rt: ParquetFile::get_double(req(&exp_rt), row, 0.0, true)?,
                    ms2_area_intensity: ParquetFile::get_double(col(&ms2_area), row, 0.0, true)?,
                    ms2_total_area_intensity: ParquetFile::get_double(
                        col(&ms2_total),
                        row,
                        0.0,
                        true,
                    )?,
                    ms2_apex_intensity: ParquetFile::get_double(col(&ms2_apex), row, 0.0, true)?,
                    ..Row::default()
                };
                let Some(&(charge, decoy)) = precursors.get(&out.precursor_id) else {
                    return Err(Error::MissingInformation(format!(
                        "Unknown precursor_id={} for run_id={}",
                        out.precursor_id, out.run_id
                    )));
                };
                out.precursor_charge = charge;
                out.decoy = decoy;
                out.transition_count = counts.get(&out.precursor_id).copied().unwrap_or(0);
                out.group_id = format!("{}_{}_{}", out.run_id, out.precursor_id, out.feature_id);
                rows.push(out);
            }
        }
        self.rows = rows;
        Ok(())
    }

    /// The path the last [`load`](Self::load) read (empty before any). Ports
    /// `oswpqPath`.
    pub fn oswpq_path(&self) -> &Path {
        &self.oswpq_dir
    }

    /// The rows the last successful [`load`](Self::load) read. Ports `rows`.
    pub fn rows(&self) -> &[Row] {
        &self.rows
    }

    /// Feature-level scores of every run of the container at `oswpq_dir`.
    ///
    /// Ports `fetchPeakGroupFeatures(oswpq_dir, level, main_score)`, whose
    /// defaults are `level = "ms2"` and `main_score = ""`. The score columns
    /// are every `var_ms2_*` or `ms2_*` column of any run's
    /// `features.parquet` (and, with level `"ms1ms2"`, every `var_ms1_*` or
    /// `ms1_*` column), sorted by name; a run without a column reads NaN
    /// there. The rows are sorted by run id, precursor id and retention
    /// time. A `main_score` among the MS2 columns is moved to the front; one
    /// that is not is reported in [`PeakGroupFeatureScoresResult::warnings`].
    ///
    /// # Errors
    ///
    /// As [`load`](Self::load), the unknown-precursor message naming the run
    /// id of `runs.parquet`.
    pub fn fetch_peak_group_features(
        &self,
        oswpq_dir: impl AsRef<Path>,
        level: &str,
        main_score: &str,
    ) -> Result<PeakGroupFeatureScoresResult> {
        let mut container = Container::new(oswpq_dir.as_ref());
        let precursors = precursor_info(&mut container)?;
        let counts = transition_counts(&mut container)?;
        let runs = run_ids(&container.table("runs/runs.parquet")?)?;

        let mut ms2_columns = Vec::new();
        let mut ms1_columns = Vec::new();
        for &run_id in &runs {
            let table = container.table(&features_entry(run_id))?;
            for field in table.schema().fields() {
                let name = field.name();
                if is_ms2(name) {
                    if !ms2_columns.iter().any(|n| n == name) {
                        ms2_columns.push(name.clone());
                    }
                } else if level == "ms1ms2"
                    && is_ms1(name)
                    && !ms1_columns.iter().any(|n| n == name)
                {
                    ms1_columns.push(name.clone());
                }
            }
        }
        ms2_columns.sort();
        ms1_columns.sort();

        let mut r = PeakGroupFeatureScoresResult {
            ms2_values: vec![Vec::new(); ms2_columns.len()],
            ms1_values: vec![Vec::new(); ms1_columns.len()],
            ..PeakGroupFeatureScoresResult::default()
        };
        for &run_id in &runs {
            let table = container.table(&features_entry(run_id))?;
            let feature_ids = ParquetFile::get_column(&table, "feature_id")?;
            let precursor_ids = ParquetFile::get_column(&table, "precursor_id")?;
            let exp_rt = ParquetFile::get_column(&table, "exp_rt")?;
            let ms2: Vec<Option<ArrayRef>> = ms2_columns
                .iter()
                .map(|c| ParquetFile::get_optional_column(&table, c))
                .collect();
            let ms1: Vec<Option<ArrayRef>> = ms1_columns
                .iter()
                .map(|c| ParquetFile::get_optional_column(&table, c))
                .collect();
            for row in 0..table.num_rows() {
                let fid = ParquetFile::get_int64(req(&feature_ids), row, 0, false)?;
                let pid = ParquetFile::get_int64(req(&precursor_ids), row, 0, false)?;
                let rt = ParquetFile::get_double(req(&exp_rt), row, 0.0, true)?;
                let Some(&(charge, decoy)) = precursors.get(&pid) else {
                    return Err(Error::MissingInformation(format!(
                        "Unknown precursor_id={pid} for run_id={run_id}"
                    )));
                };
                r.feature_id.push(fid);
                r.run_id.push(run_id);
                r.precursor_id.push(pid);
                r.exp_rt.push(rt);
                r.precursor_charge.push(charge);
                r.decoy.push(decoy);
                r.transition_count
                    .push(counts.get(&pid).copied().unwrap_or(0));
                r.group_id.push(format!("{run_id}_{pid}_{fid}"));
                for (values, array) in r.ms2_values.iter_mut().zip(&ms2) {
                    values.push(nan_or(col(array), row)?);
                }
                for (values, array) in r.ms1_values.iter_mut().zip(&ms1) {
                    values.push(nan_or(col(array), row)?);
                }
            }
        }

        let order = source_sort_permutation(r.feature_id.len(), |a, b| {
            if r.run_id[a] != r.run_id[b] {
                return r.run_id[a] < r.run_id[b];
            }
            if r.precursor_id[a] != r.precursor_id[b] {
                return r.precursor_id[a] < r.precursor_id[b];
            }
            r.exp_rt[a] < r.exp_rt[b]
        })?;
        permute(&mut r.feature_id, &order);
        permute(&mut r.run_id, &order);
        permute(&mut r.precursor_id, &order);
        permute(&mut r.exp_rt, &order);
        permute(&mut r.precursor_charge, &order);
        permute(&mut r.decoy, &order);
        permute(&mut r.transition_count, &order);
        permute(&mut r.group_id, &order);
        for values in r.ms2_values.iter_mut().chain(r.ms1_values.iter_mut()) {
            permute(values, &order);
        }
        r.ms2_columns = ms2_columns;
        r.ms1_columns = ms1_columns;

        if !main_score.is_empty() && !r.ms2_columns.is_empty() {
            match r.ms2_columns.iter().position(|c| c == main_score) {
                Some(0) => {}
                Some(found) => {
                    let name = r.ms2_columns.remove(found);
                    let values = r.ms2_values.remove(found);
                    r.ms2_columns.insert(0, name);
                    r.ms2_values.insert(0, values);
                }
                None => r.warnings.push(format!(
                    "Requested main_score '{main_score}' not found among discovered MS2 columns"
                )),
            }
        }
        Ok(r)
    }

    /// Transition-level scores of every run of the container at `oswpq_dir`.
    ///
    /// Ports `fetchTransitionFeatures`. Every row of a run's
    /// `feature_transition.parquet` is joined with its feature (precursor
    /// and retention time from `features.parquet`), the library precursor
    /// (charge) and the library transition (product charge, decoy flag). The
    /// `var_*` columns of any run are collected, sorted by name. Missing
    /// values are NaN. The rows are sorted by run id, precursor id,
    /// retention time and transition id.
    ///
    /// # Errors
    ///
    /// As [`load`](Self::load); [`Error::MissingInformation`]
    /// "feature_transition.parquet references unknown feature_id=\<id\> for
    /// run_id=\<id\>", "Unknown precursor_id=…" and "Unknown
    /// transition_id=\<id\> for run_id=\<id\>", as the source.
    pub fn fetch_transition_features(
        &self,
        oswpq_dir: impl AsRef<Path>,
    ) -> Result<TransitionFeaturesResult> {
        let mut container = Container::new(oswpq_dir.as_ref());
        let precursor_charge = {
            let table = container.table("library/precursors.parquet")?;
            let ids = ParquetFile::get_column(&table, "precursor_id")?;
            let charges = ParquetFile::get_column(&table, "charge")?;
            let mut map = HashMap::with_capacity(table.num_rows());
            for row in 0..table.num_rows() {
                let id = ParquetFile::get_int64(req(&ids), row, 0, false)?;
                let charge = ParquetFile::get_int64(req(&charges), row, 0, false)? as i32;
                map.entry(id).or_insert(charge);
            }
            map
        };
        let transition_info = {
            let table = container.table("library/transitions.parquet")?;
            let ids = ParquetFile::get_column(&table, "transition_id")?;
            let charges = ParquetFile::get_column(&table, "charge")?;
            let decoys = ParquetFile::get_optional_column(&table, "decoy");
            let mut map = HashMap::with_capacity(table.num_rows());
            for row in 0..table.num_rows() {
                let id = ParquetFile::get_int64(req(&ids), row, 0, false)?;
                let charge = ParquetFile::get_int64(req(&charges), row, 0, false)? as i32;
                let decoy = ParquetFile::get_bool(col(&decoys), row, false, true)?;
                map.entry(id).or_insert((charge, decoy));
            }
            map
        };
        let runs = run_ids(&container.table("runs/runs.parquet")?)?;

        let mut var_columns = Vec::new();
        for &run_id in &runs {
            let table = container.table(&feature_transition_entry(run_id))?;
            discover(&table, &mut var_columns, |name| name.starts_with("var_"));
        }
        var_columns.sort();

        let mut r = TransitionFeaturesResult {
            transition_var_values: vec![Vec::new(); var_columns.len()],
            ..TransitionFeaturesResult::default()
        };
        for &run_id in &runs {
            let mut features: HashMap<i64, (i64, f64)> = HashMap::new();
            {
                let table = container.table(&features_entry(run_id))?;
                let ids = ParquetFile::get_column(&table, "feature_id")?;
                let precursors = ParquetFile::get_column(&table, "precursor_id")?;
                let rts = ParquetFile::get_column(&table, "exp_rt")?;
                for row in 0..table.num_rows() {
                    let fid = ParquetFile::get_int64(req(&ids), row, 0, false)?;
                    let pid = ParquetFile::get_int64(req(&precursors), row, 0, false)?;
                    let rt = ParquetFile::get_double(req(&rts), row, 0.0, true)?;
                    features.entry(fid).or_insert((pid, rt));
                }
            }
            let table = container.table(&feature_transition_entry(run_id))?;
            let feature_ids = ParquetFile::get_column(&table, "feature_id")?;
            let transition_ids = ParquetFile::get_column(&table, "transition_id")?;
            let optional = |name: &str| ParquetFile::get_optional_column(&table, name);
            let (area, total_area, apex, apex_rt, fwhm, masserror, total_mi) = (
                optional("area_intensity"),
                optional("total_area_intensity"),
                optional("apex_intensity"),
                optional("apex_rt"),
                optional("rt_fwhm"),
                optional("masserror_ppm"),
                optional("total_mi"),
            );
            let vars: Vec<Option<ArrayRef>> = var_columns.iter().map(|c| optional(c)).collect();
            for row in 0..table.num_rows() {
                let fid = ParquetFile::get_int64(req(&feature_ids), row, 0, false)?;
                let tid = ParquetFile::get_int64(req(&transition_ids), row, 0, false)?;
                let Some(&(pid, rt)) = features.get(&fid) else {
                    return Err(Error::MissingInformation(format!(
                        "feature_transition.parquet references unknown feature_id={fid} for run_id={run_id}"
                    )));
                };
                let Some(&charge) = precursor_charge.get(&pid) else {
                    return Err(Error::MissingInformation(format!(
                        "Unknown precursor_id={pid} for run_id={run_id}"
                    )));
                };
                let Some(&(product_charge, decoy)) = transition_info.get(&tid) else {
                    return Err(Error::MissingInformation(format!(
                        "Unknown transition_id={tid} for run_id={run_id}"
                    )));
                };
                r.feature_id.push(fid);
                r.run_id.push(run_id);
                r.precursor_id.push(pid);
                r.exp_rt.push(rt);
                r.precursor_charge.push(charge);
                r.transition_id.push(tid);
                r.product_charge.push(product_charge);
                r.decoy.push(decoy);
                r.area_intensity.push(nan_or(col(&area), row)?);
                r.total_area_intensity.push(nan_or(col(&total_area), row)?);
                r.apex_intensity.push(nan_or(col(&apex), row)?);
                r.apex_rt.push(nan_or(col(&apex_rt), row)?);
                r.rt_fwhm.push(nan_or(col(&fwhm), row)?);
                r.masserror_ppm.push(nan_or(col(&masserror), row)?);
                r.total_mi.push(nan_or(col(&total_mi), row)?);
                r.group_id.push(format!("{run_id}_{fid}_{pid}_{tid}"));
                for (values, array) in r.transition_var_values.iter_mut().zip(&vars) {
                    values.push(nan_or(col(array), row)?);
                }
            }
        }

        let order = source_sort_permutation(r.feature_id.len(), |a, b| {
            if r.run_id[a] != r.run_id[b] {
                return r.run_id[a] < r.run_id[b];
            }
            if r.precursor_id[a] != r.precursor_id[b] {
                return r.precursor_id[a] < r.precursor_id[b];
            }
            if r.exp_rt[a] != r.exp_rt[b] {
                return r.exp_rt[a] < r.exp_rt[b];
            }
            r.transition_id[a] < r.transition_id[b]
        })?;
        permute(&mut r.feature_id, &order);
        permute(&mut r.run_id, &order);
        permute(&mut r.precursor_id, &order);
        permute(&mut r.exp_rt, &order);
        permute(&mut r.precursor_charge, &order);
        permute(&mut r.transition_id, &order);
        permute(&mut r.product_charge, &order);
        permute(&mut r.decoy, &order);
        permute(&mut r.area_intensity, &order);
        permute(&mut r.total_area_intensity, &order);
        permute(&mut r.apex_intensity, &order);
        permute(&mut r.apex_rt, &order);
        permute(&mut r.rt_fwhm, &order);
        permute(&mut r.masserror_ppm, &order);
        permute(&mut r.total_mi, &order);
        permute(&mut r.group_id, &order);
        for values in &mut r.transition_var_values {
            permute(values, &order);
        }
        r.transition_var_columns = var_columns;
        Ok(r)
    }

    /// The PyProphet "unscored" table of every run of the container at
    /// `oswpq_dir`, in run and file order.
    ///
    /// Ports `fetchUnscoredData`. Library values come from
    /// `precursors.parquet` (for a repeated precursor id the last row wins),
    /// whose `charge`, `decoy`, `precursor_mz` and `library_rt` columns are
    /// all optional; `id_peptide` from `library/precursor_peptide_mapping.parquet`
    /// when it and `library/peptides.parquet` both exist. The `var_ms2_*`,
    /// `ms2_*`, `var_ms1_*` and `ms1_*` columns of any run are collected,
    /// sorted by name. Missing values are NaN, 0, `false` or empty.
    ///
    /// # Errors
    ///
    /// As [`load`](Self::load), for the tables this reads.
    pub fn fetch_unscored_data(&self, oswpq_dir: impl AsRef<Path>) -> Result<UnscoredResult> {
        let path = oswpq_dir.as_ref();
        let mut container = Container::new(path);
        let table = container.table("library/precursors.parquet")?;
        let ids = ParquetFile::get_column(&table, "precursor_id")?;
        let charge_col = ParquetFile::get_optional_column(&table, "charge");
        let decoy_col = ParquetFile::get_optional_column(&table, "decoy");
        let mz_col = ParquetFile::get_optional_column(&table, "precursor_mz");
        let rt_col = ParquetFile::get_optional_column(&table, "library_rt");
        let mut charge: HashMap<i64, i32> = HashMap::new();
        let mut decoy: HashMap<i64, bool> = HashMap::new();
        let mut mz: HashMap<i64, f64> = HashMap::new();
        let mut library_rt: HashMap<i64, f64> = HashMap::new();
        for row in 0..table.num_rows() {
            let pid = ParquetFile::get_int64(req(&ids), row, 0, false)?;
            if charge_col.is_some() {
                charge.insert(
                    pid,
                    ParquetFile::get_int64(col(&charge_col), row, 0, false)? as i32,
                );
            }
            if decoy_col.is_some() {
                decoy.insert(
                    pid,
                    ParquetFile::get_bool(col(&decoy_col), row, false, true)?,
                );
            }
            if mz_col.is_some() {
                mz.insert(pid, nan_or(col(&mz_col), row)?);
            }
            if rt_col.is_some() {
                library_rt.insert(pid, nan_or(col(&rt_col), row)?);
            }
        }

        const MAPPING: &str = "library/precursor_peptide_mapping.parquet";
        const PEPTIDES: &str = "library/peptides.parquet";
        let (has_mapping, has_peptides) = if file::is_directory(path) {
            (
                file::exists(path.join(MAPPING)),
                file::exists(path.join(PEPTIDES)),
            )
        } else {
            let entries = ZipArchiveFile::list_entries(path);
            (
                entries.iter().any(|e| e == MAPPING),
                entries.iter().any(|e| e == PEPTIDES),
            )
        };
        let mut peptide: HashMap<i64, i64> = HashMap::new();
        if has_mapping && has_peptides {
            let mapping = container.table(MAPPING)?;
            // Read to check that it is present and readable, as the source does.
            container.table(PEPTIDES)?;
            let precursors = ParquetFile::get_column(&mapping, "precursor_id")?;
            let peptides = ParquetFile::get_column(&mapping, "peptide_id")?;
            for row in 0..mapping.num_rows() {
                let pid = ParquetFile::get_int64(req(&precursors), row, 0, false)?;
                let pep = ParquetFile::get_int64(req(&peptides), row, 0, false)?;
                peptide.insert(pid, pep);
            }
        }

        let runs_table = container.table("runs/runs.parquet")?;
        let runs = run_ids(&runs_table)?;
        let filenames = ParquetFile::get_optional_column(&runs_table, "filename");

        let mut ms2_columns = Vec::new();
        let mut ms1_columns = Vec::new();
        for &run_id in &runs {
            let table = container.table(&features_entry(run_id))?;
            for field in table.schema().fields() {
                let name = field.name();
                if is_ms2(name) {
                    if !ms2_columns.iter().any(|n| n == name) {
                        ms2_columns.push(name.clone());
                    }
                } else if is_ms1(name) && !ms1_columns.iter().any(|n| n == name) {
                    ms1_columns.push(name.clone());
                }
            }
        }
        ms2_columns.sort();
        ms1_columns.sort();

        let mut r = UnscoredResult {
            ms2_values: vec![Vec::new(); ms2_columns.len()],
            ms1_values: vec![Vec::new(); ms1_columns.len()],
            ..UnscoredResult::default()
        };
        for (run_index, &run_id) in runs.iter().enumerate() {
            let filename = if filenames.is_some() {
                ParquetFile::get_string(col(&filenames), run_index)?
            } else {
                String::new()
            };
            let table = container.table(&features_entry(run_id))?;
            let optional = |name: &str| ParquetFile::get_optional_column(&table, name);
            let feature_ids = ParquetFile::get_column(&table, "feature_id")?;
            let precursor_ids = ParquetFile::get_column(&table, "precursor_id")?;
            let exp_rt = optional("exp_rt");
            let delta_rt = optional("delta_rt");
            let norm_rt = optional("norm_rt");
            let ms2_area = optional("ms2_area_intensity");
            let ms1_area = optional("ms1_area_intensity");
            let ms1_apex = optional("ms1_apex_intensity");
            let left = optional("left_width");
            let right = optional("right_width");
            let exp_im = optional("exp_im");
            let im_left = optional("exp_im_leftwidth");
            let im_right = optional("exp_im_rightwidth");
            let ms2: Vec<Option<ArrayRef>> = ms2_columns.iter().map(|c| optional(c)).collect();
            let ms1: Vec<Option<ArrayRef>> = ms1_columns.iter().map(|c| optional(c)).collect();
            for row in 0..table.num_rows() {
                let fid = ParquetFile::get_int64(req(&feature_ids), row, 0, false)?;
                let pid = ParquetFile::get_int64(req(&precursor_ids), row, 0, false)?;
                let rt = nan_or(col(&exp_rt), row)?;
                let drt = nan_or(col(&delta_rt), row)?;
                let nrt = nan_or(col(&norm_rt), row)?;
                let lib_rt = library_rt.get(&pid).copied().unwrap_or(f64::NAN);
                r.id_run.push(run_id);
                r.id_peptide.push(peptide.get(&pid).copied().unwrap_or(0));
                r.transition_group_id.push(pid);
                r.decoy.push(decoy.get(&pid).copied().unwrap_or(false));
                r.run_id.push(run_id);
                r.filename.push(filename.clone());
                r.rt.push(rt);
                r.assay_rt.push(rt - drt);
                r.delta_rt.push(drt);
                r.assay_rt_library.push(lib_rt);
                r.delta_rt_library.push(nrt - lib_rt);
                r.id.push(fid);
                r.charge.push(charge.get(&pid).copied().unwrap_or(0));
                r.mz.push(mz.get(&pid).copied().unwrap_or(f64::NAN));
                r.intensity.push(nan_or(col(&ms2_area), row)?);
                r.aggr_prec_peak_area.push(nan_or(col(&ms1_area), row)?);
                r.aggr_prec_peak_apex.push(nan_or(col(&ms1_apex), row)?);
                r.left_width.push(nan_or(col(&left), row)?);
                r.right_width.push(nan_or(col(&right), row)?);
                r.exp_im.push(nan_or(col(&exp_im), row)?);
                r.im_left_width.push(nan_or(col(&im_left), row)?);
                r.im_right_width.push(nan_or(col(&im_right), row)?);
                for (values, array) in r.ms2_values.iter_mut().zip(&ms2) {
                    values.push(nan_or(col(array), row)?);
                }
                for (values, array) in r.ms1_values.iter_mut().zip(&ms1) {
                    values.push(nan_or(col(array), row)?);
                }
            }
        }
        r.ms2_columns = ms2_columns;
        r.ms1_columns = ms1_columns;
        Ok(r)
    }
}
