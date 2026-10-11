// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `OpenSwathOSWParquetWriter` and `OpenSwathOSWParquetReader`.
//!
//! Expected values come from:
//!
//! - the class tests `OpenSwathOSWParquetWriter_test.cpp`,
//!   `OpenSwathOSWParquetReader_test.cpp` and
//!   `OpenSwathOSWParquetRoundTrip_test.cpp` at core `bc9cc12`, transcribed,
//!   with the upstream fixture `OpenSwathWorkflow_tworuns_1_17.output.oswpq`;
//! - the executed Release oracle `../oracle/O5c-oswpq` (ibminode06,
//!   `openms4-release-bc9cc12-c19e494-174b576`, Arrow C++ 23.0.1, run twice
//!   with identical records and files): `tests/data/openswath_osw_parquet/`
//!   holds its record stream and every file the Release writer produced. The
//!   Rust writer runs on the same inputs and every table is compared value by
//!   value; the Rust reader runs on the Release-written containers and every
//!   returned column is compared with the Release reader's.
#![cfg(all(feature = "arrow", feature = "file-compression"))]

use openms::Error;
use openms::analysis::openswath::openswath_library_id_normalizer::{
    OpenSwathLibraryIDNormalizer, SourceIDMapping,
};
use openms::analysis::openswath::openswath_osw_parquet_reader::OpenSwathOSWParquetReader;
use openms::analysis::openswath::openswath_osw_parquet_writer::OpenSwathOSWParquetWriter;
use openms::analysis::openswath::transition_experiment::{
    LightCompound, LightProtein, LightTargetedExperiment, LightTransition,
};
use openms::analysis::openswath::transition_parquet_file::TransitionParquetFile;
use openms::format::parquet_file::ParquetFile;
use openms::format::zip_archive_file::ZipArchiveFile;
use openms::kernel::{Feature, FeatureMap};
use openms::metadata::MetaValue;
use openms::system::file::TempDir;

use arrow_array::{
    Array, ArrayRef, BooleanArray, Float64Array, Int32Array, Int64Array, RecordBatch,
};
use arrow_schema::{Field, Schema};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;

const DATA: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/data/openswath_osw_parquet"
);
const TPF_DATA: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/data/transition_parquet_file"
);
const TWORUNS: &str = "OpenSwathWorkflow_tworuns_1_17.output.oswpq";

// ------------------------------------------------------------------ records

fn unescape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c == '\\' {
            match chars.next() {
                Some('t') => out.push('\t'),
                Some('n') => out.push('\n'),
                Some('\\') => out.push('\\'),
                Some(other) => {
                    out.push('\\');
                    out.push(other);
                }
                None => out.push('\\'),
            }
        } else {
            out.push(c);
        }
    }
    out
}

fn read_records(dir: &str) -> Vec<Vec<String>> {
    std::fs::read_to_string(Path::new(dir).join("records.tsv"))
        .unwrap()
        .lines()
        .map(|line| line.split('\t').map(unescape).collect())
        .collect()
}

fn records() -> Vec<Vec<String>> {
    read_records(DATA)
}

fn of<'a>(all: &'a [Vec<String>], kind: &str, case: &str) -> Vec<&'a Vec<String>> {
    all.iter()
        .filter(|r| r[0] == kind && r.get(1).map(String::as_str) == Some(case))
        .collect()
}

fn outcome<'a>(all: &'a [Vec<String>], case: &str) -> &'a Vec<String> {
    all.iter()
        .find(|r| (r[0] == "OK" || r[0] == "E") && r[1] == case)
        .unwrap_or_else(|| panic!("no outcome record for {case}"))
}

fn fact(all: &[Vec<String>], case: &str, key: &str) -> String {
    all.iter()
        .find(|r| r[0] == "X" && r[1] == case && r[2] == key)
        .unwrap_or_else(|| panic!("no fact {key} for {case}"))[3]
        .clone()
}

fn num(s: &str) -> f64 {
    match s {
        "nan" | "-nan" => f64::NAN,
        "inf" => f64::INFINITY,
        "-inf" => f64::NEG_INFINITY,
        _ => s.parse().unwrap_or_else(|_| panic!("bad number {s}")),
    }
}

fn same(a: f64, b: f64) -> bool {
    (a.is_nan() && b.is_nan()) || a.to_bits() == b.to_bits()
}

fn assert_same(ctx: &str, got: f64, want: &str) {
    let want = num(want);
    assert!(same(got, want), "{ctx}: got {got:?}, want {want:?}");
}

fn assert_values(ctx: &str, got: &[Vec<f64>], row: usize, want: &str) {
    let want: Vec<&str> = if want.is_empty() {
        Vec::new()
    } else {
        want.split('|').collect()
    };
    assert_eq!(got.len(), want.len(), "{ctx}: value column count");
    for (column, text) in got.iter().zip(want) {
        assert_same(ctx, column[row], text);
    }
}

/// The error a C++ exception record corresponds to.
fn assert_error_matches(case: &str, err: &Error, record: &[String]) {
    let kind = record[2].as_str();
    let what = record[3].as_str();
    let (value, message) = match what.strip_prefix("the value '") {
        Some(rest) => rest.split_once("' was used but is not valid; ").unwrap(),
        None => ("", what),
    };
    // The work directory differs; compare the value's last path component.
    let value = value.rsplit('/').next().unwrap_or(value);
    match (kind, err) {
        ("OpenMS::Exception::InvalidValue", Error::InvalidValue(m)) => {
            assert!(m.contains(message), "{case}: '{m}' lacks '{message}'");
            assert!(m.contains(value), "{case}: '{m}' lacks '{value}'");
        }
        ("OpenMS::Exception::MissingInformation", Error::MissingInformation(m)) => {
            assert_eq!(m, message, "{case}");
        }
        ("OpenMS::Exception::FileNotFound", Error::Io(e)) => {
            assert_eq!(e.kind(), std::io::ErrorKind::NotFound, "{case}: {e}");
        }
        _ => panic!("{case}: C++ threw {kind} ({what}), Rust returned {err:?}"),
    }
}

fn tmp() -> TempDir {
    TempDir::new_in(std::env::temp_dir(), false).unwrap()
}

/// A container directory assembled from the oracle's flattened files of `source`.
fn container(all: &[Vec<String>], tmp: &Path, source: &str, name: &str) -> PathBuf {
    let dir = tmp.join(name);
    let files = Path::new(DATA).join("files");
    let entries = of(all, "F", source);
    assert!(!entries.is_empty(), "no files for {source}");
    for r in entries {
        let target = dir.join(&r[2]);
        std::fs::create_dir_all(target.parent().unwrap()).unwrap();
        std::fs::copy(files.join(&r[3]), target).unwrap();
    }
    dir
}

// ------------------------------------------------------------------ file comparison

/// Documented divergence (writer module docs, "Differences from the
/// source"). The source writes every extracted score as
/// `appendOptionalFloat_(b, extractMetaDouble_(x, key, value), value, ..)`
/// with `value` (declared `double value = 0.0;` once per feature) passed by
/// value, and C++ leaves the order in which arguments are evaluated
/// unspecified. The executed Release build (GCC 14.4) copies `value` before
/// the extraction runs, so each such cell carries the value of the previous
/// successful extraction on the same feature (0.0 for the first), while the
/// null pattern (the extraction's success) is right. The port writes the
/// cell's own value.
///
/// `STALE` lists the successful extractions of the features that have any,
/// in the source's order: feature A of `rich_features()` (run 42) and the
/// feature of `second_run_features()` (runs 7 and 8), as
/// `(file, column, row, extracted value, first extraction of its feature)`.
/// Extractions that fail (a list, an empty or unparseable text) are not in
/// it: they leave `value` alone. The
/// Release cell is `finite(previous)` and the port's `finite(own)`; the cells
/// where they differ must be exactly the differences between the two files.
const STALE: &[(&str, &str, usize, f64, bool)] = &[
    // ms1_area_intensity
    (
        "runs/run_id=42/features.parquet",
        "ms1_area_intensity",
        0,
        123456.78901234567,
        true,
    ),
    // ms1_apex_intensity
    (
        "runs/run_id=42/features.parquet",
        "ms1_apex_intensity",
        0,
        1e-300,
        false,
    ),
    // im_ms1_drift
    (
        "runs/run_id=42/features.parquet",
        "ms1_exp_im",
        0,
        5e-324,
        false,
    ),
    // im_ms1_delta
    (
        "runs/run_id=42/features.parquet",
        "ms1_delta_im",
        0,
        1e308,
        false,
    ),
    // var_ms1_ppm_diff
    (
        "runs/run_id=42/features.parquet",
        "var_ms1_massdev_score",
        0,
        1.5,
        false,
    ),
    // var_im_ms1_delta_score
    (
        "runs/run_id=42/features.parquet",
        "var_ms1_im_ms1_delta_score",
        0,
        42.0,
        false,
    ),
    // var_ms1_mi_score "nan"
    (
        "runs/run_id=42/features.parquet",
        "var_ms1_mi_score",
        0,
        f64::NAN,
        false,
    ),
    // var_ms1_mi_contrast_score "inf"
    (
        "runs/run_id=42/features.parquet",
        "var_ms1_mi_contrast_score",
        0,
        f64::INFINITY,
        false,
    ),
    // var_ms1_mi_combined_score "+2.5e3"
    (
        "runs/run_id=42/features.parquet",
        "var_ms1_mi_combined_score",
        0,
        2500.0,
        false,
    ),
    // var_ms1_xcorr_coelution
    (
        "runs/run_id=42/features.parquet",
        "var_ms1_xcorr_coelution",
        0,
        0.3,
        false,
    ),
    // var_ms1_xcorr_shape
    (
        "runs/run_id=42/features.parquet",
        "var_ms1_xcorr_shape",
        0,
        2.675,
        false,
    ),
    // total_xic
    (
        "runs/run_id=42/features.parquet",
        "ms2_total_area_intensity",
        0,
        1234567890123.0,
        false,
    ),
    // peak_apices_sum
    (
        "runs/run_id=42/features.parquet",
        "ms2_apex_intensity",
        0,
        9007199254740992.0,
        false,
    ),
    // im_drift
    (
        "runs/run_id=42/features.parquet",
        "ms2_exp_im",
        0,
        0.1,
        false,
    ),
    // im_drift_left " 7 "
    (
        "runs/run_id=42/features.parquet",
        "ms2_exp_im_leftwidth",
        0,
        7.0,
        false,
    ),
    // var_bseries_score
    (
        "runs/run_id=42/features.parquet",
        "var_ms2_bseries_score",
        0,
        3.0,
        false,
    ),
    // var_library_corr
    (
        "runs/run_id=42/features.parquet",
        "var_ms2_library_corr",
        0,
        0.987654321987654,
        false,
    ),
    // var_manhatt_score
    (
        "runs/run_id=42/features.parquet",
        "var_ms2_manhattan_score",
        0,
        12.0,
        false,
    ),
    // var_xcorr_shape
    (
        "runs/run_id=42/features.parquet",
        "var_ms2_xcorr_shape",
        0,
        0.8,
        false,
    ),
    // var_xcorr_shape_weighted
    (
        "runs/run_id=42/features.parquet",
        "var_ms2_xcorr_shape_weighted",
        0,
        0.3,
        false,
    ),
    // im_log_intensity
    (
        "runs/run_id=42/features.parquet",
        "var_ms2_im_log_intensity",
        0,
        7.25,
        false,
    ),
    // s1 total_xic
    (
        "runs/run_id=42/feature_transition.parquet",
        "total_area_intensity",
        0,
        20.0,
        false,
    ),
    // s1 peak_apex_int
    (
        "runs/run_id=42/feature_transition.parquet",
        "apex_intensity",
        0,
        5.0,
        false,
    ),
    // s1 peak_apex_position
    (
        "runs/run_id=42/feature_transition.parquet",
        "apex_rt",
        0,
        100.3,
        false,
    ),
    // s1 width_at_50
    (
        "runs/run_id=42/feature_transition.parquet",
        "rt_fwhm",
        0,
        2.2,
        false,
    ),
    // s1 total_mi
    (
        "runs/run_id=42/feature_transition.parquet",
        "total_mi",
        0,
        0.7,
        false,
    ),
    // s1 start_position_at_5
    (
        "runs/run_id=42/feature_transition.parquet",
        "start_position_at_5",
        0,
        99.0,
        false,
    ),
    // s1 end_position_at_5
    (
        "runs/run_id=42/feature_transition.parquet",
        "end_position_at_5",
        0,
        101.0,
        false,
    ),
    // s1 tailing_factor "1.25"
    (
        "runs/run_id=42/feature_transition.parquet",
        "tailing_factor",
        0,
        1.25,
        false,
    ),
    // s1 points_across_half_height
    (
        "runs/run_id=42/feature_transition.parquet",
        "points_across_half_height",
        0,
        7.0,
        false,
    ),
    // s_iso peak_apex_int
    (
        "runs/run_id=42/feature_precursor.parquet",
        "precursor_apex_intensity",
        1,
        3.25,
        false,
    ),
    // second_run_features(): its feature's only extraction, written as run 7
    // (appended) and as run 8 (the replaced archive).
    (
        "runs/run_id=7/features.parquet",
        "var_ms2_xcorr_shape",
        0,
        0.5,
        true,
    ),
    (
        "runs/run_id=8/features.parquet",
        "var_ms2_xcorr_shape",
        0,
        0.5,
        true,
    ),
];

fn finite(v: f64) -> Option<f64> {
    v.is_finite().then_some(v)
}

/// The Release and port cells of `STALE[k]`.
fn stale_pair(k: usize) -> (Option<f64>, Option<f64>) {
    let previous = if STALE[k].4 { 0.0 } else { STALE[k - 1].3 };
    (finite(previous), finite(STALE[k].3))
}

fn same_cell(a: Option<f64>, b: Option<f64>) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => x.to_bits() == y.to_bits(),
        (None, None) => true,
        _ => false,
    }
}

/// The indices of `STALE` in the given runs whose Release and port cells
/// differ.
fn stale_differences(runs: &[i64]) -> BTreeSet<usize> {
    (0..STALE.len())
        .filter(|&k| {
            let (release, port) = stale_pair(k);
            runs.iter()
                .any(|r| STALE[k].0.starts_with(&format!("runs/run_id={r}/")))
                && !same_cell(release, port)
        })
        .collect()
}

fn float_at(a: &ArrayRef, row: usize) -> Option<f64> {
    let a = a.as_any().downcast_ref::<Float64Array>()?;
    if a.is_null(row) {
        None
    } else {
        Some(a.value(row))
    }
}

/// Compares every column of `got` with `want`; returns the [`STALE`]
/// entries it had to apply. Panics listing every other difference.
fn assert_tables_equal(
    ctx: &str,
    rel: &str,
    got: &Path,
    want: &Path,
    problems: &mut Vec<String>,
) -> BTreeSet<usize> {
    let g = ParquetFile::read_table(got).unwrap_or_else(|e| panic!("{ctx}: {e}"));
    let w = ParquetFile::read_table(want).unwrap();
    assert_eq!(g.num_rows(), w.num_rows(), "{ctx} rows");
    let fields = |s: &Schema| {
        s.fields()
            .iter()
            .map(|f| (f.name().clone(), f.data_type().clone(), f.is_nullable()))
            .collect::<Vec<_>>()
    };
    assert_eq!(fields(&g.schema()), fields(&w.schema()), "{ctx} schema");
    let mut used = BTreeSet::new();
    for (i, field) in w.schema().fields().iter().enumerate() {
        let (gc, wc) = (g.column(i), w.column(i));
        if gc.to_data() == wc.to_data() {
            continue;
        }
        if field.data_type() != &arrow_schema::DataType::Float64 {
            problems.push(format!("{ctx}: {} (non-float column)", field.name()));
            continue;
        }
        for row in 0..wc.len() {
            let (gv, wv) = (float_at(gc, row), float_at(wc, row));
            let equal = same_cell(gv, wv);
            if equal {
                continue;
            }
            let stale = (0..STALE.len()).find(|&k| {
                let (file, column, r, _, _) = STALE[k];
                let (release, port) = stale_pair(k);
                file == rel
                    && column == field.name()
                    && r == row
                    && same_cell(wv, release)
                    && same_cell(gv, port)
            });
            match stale {
                Some(k) => {
                    used.insert(k);
                }
                None => problems.push(format!(
                    "{ctx}: {} row {row}: got {gv:?}, want {wv:?}",
                    field.name()
                )),
            }
        }
    }
    used
}

fn mask_build_time(text: &str) -> String {
    text.lines()
        .map(|line| {
            if line.trim_start().starts_with("\"build_time\"") {
                "    \"build_time\": <masked>,".to_owned()
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn json_keys(text: &str) -> BTreeSet<String> {
    text.lines()
        .filter_map(|l| l.trim().strip_prefix('"'))
        .map(|l| l.split('"').next().unwrap().to_owned())
        .collect()
}

/// Every file the Release writer produced for `case` exists in `dir` and
/// agrees with it; `dir` holds no other file.
/// Returns the [`STALE`] entries applied.
fn assert_container_matches(all: &[Vec<String>], case: &str, dir: &Path) -> BTreeSet<usize> {
    let mut used = BTreeSet::new();
    let mut problems = Vec::new();
    let files = Path::new(DATA).join("files");
    let entries = of(all, "F", case);
    assert!(!entries.is_empty(), "no files for {case}");
    let mut want_names = BTreeSet::new();
    for r in entries {
        let rel = &r[2];
        want_names.insert(rel.clone());
        let ctx = format!("{case} {rel}");
        let got = dir.join(rel);
        let want = files.join(&r[3]);
        assert!(got.exists(), "{ctx}: missing");
        if rel.ends_with(".parquet") {
            used.extend(assert_tables_equal(&ctx, rel, &got, &want, &mut problems));
        } else if rel.ends_with(".idx.json") {
            let g = std::fs::read_to_string(&got).unwrap();
            let w = std::fs::read_to_string(&want).unwrap();
            assert_eq!(json_keys(&g), json_keys(&w), "{ctx}");
        } else {
            let g = std::fs::read_to_string(&got).unwrap();
            let w = std::fs::read_to_string(&want).unwrap();
            assert!(g.contains("\"build_time\": \"\""), "{ctx}");
            assert_eq!(mask_build_time(&g), mask_build_time(&w), "{ctx}");
        }
    }
    let mut got_names = BTreeSet::new();
    let mut pending = vec![dir.to_path_buf()];
    while let Some(d) = pending.pop() {
        for entry in std::fs::read_dir(&d).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                pending.push(path);
            } else {
                let rel = path
                    .strip_prefix(dir)
                    .unwrap()
                    .to_string_lossy()
                    .replace('\\', "/");
                got_names.insert(rel);
            }
        }
    }
    assert_eq!(got_names, want_names, "{case}: file set");
    assert!(problems.is_empty(), "{problems:#?}");
    used
}

/// Exactly the [`STALE`] differences of `runs` were applied.
fn assert_all_stale(case: &str, used: BTreeSet<usize>, runs: &[i64]) {
    assert_eq!(used, stale_differences(runs), "{case}: stale entries");
}

// ------------------------------------------------------------------ inputs

fn f(value: f64) -> MetaValue {
    MetaValue::try_from(value).unwrap()
}

fn fl(values: &[f64]) -> MetaValue {
    MetaValue::try_from(values.to_vec()).unwrap()
}

fn sl(values: &[&str]) -> MetaValue {
    MetaValue::from(values.iter().map(|s| s.to_string()).collect::<Vec<_>>())
}

fn il(values: &[i64]) -> MetaValue {
    MetaValue::from(values.to_vec())
}

fn set(feature: &mut Feature, key: &str, value: MetaValue) {
    feature.metadata.insert(key.to_owned(), value);
}

fn feature(unique_id: u64, rt: f64, intensity: f32) -> Feature {
    let mut feature = Feature::default();
    feature.unique_id = unique_id;
    feature.rt = rt;
    feature.intensity = intensity;
    feature
}

/// The driver's `sub(level, native_id, intensity)`.
fn sub(level: &str, native_id: Option<MetaValue>, intensity: f32) -> Feature {
    let mut s = feature(0, 0.0, intensity);
    set(&mut s, "FeatureLevel", MetaValue::from(level));
    if let Some(id) = native_id {
        set(&mut s, "native_id", id);
    }
    s
}

fn class_library() -> LightTargetedExperiment {
    let compound = |id: &str, charge: i32, rt: f64, sequence: &str| LightCompound {
        id: id.into(),
        charge,
        rt,
        sequence: sequence.into(),
        ..LightCompound::default()
    };
    let transition =
        |name: &str, reference: &str, precursor_mz: f64, product_mz: f64| LightTransition {
            transition_name: name.into(),
            peptide_ref: reference.into(),
            precursor_mz,
            product_mz,
            ..LightTransition::default()
        };
    LightTargetedExperiment {
        compounds: vec![
            compound("0", 2, 10.0, "PEPTIDEA"),
            compound("7", 1, 20.0, "PEPTIDEB"),
        ],
        transitions: vec![
            transition("3", "0", 500.0, 100.0),
            transition("100", "7", 600.0, 200.0),
        ],
        proteins: Vec::new(),
    }
}

fn class_features() -> FeatureMap {
    let mut fmap = FeatureMap::default();
    let mut f1 = feature(1, 10.0, 100.0);
    set(&mut f1, "PeptideRef", MetaValue::from("0"));
    let mut s1 = feature(0, 10.0, 10.0);
    set(&mut s1, "FeatureLevel", MetaValue::from("MS2"));
    set(&mut s1, "native_id", MetaValue::from("3"));
    f1.subordinates.push(s1);
    let mut f2 = feature(2, 20.0, 200.0);
    set(&mut f2, "PeptideRef", MetaValue::from("7"));
    let mut s2 = feature(0, 20.0, 20.0);
    set(&mut s2, "FeatureLevel", MetaValue::from("MS2"));
    set(&mut s2, "native_id", MetaValue::from("100"));
    f2.subordinates.push(s2);
    fmap.features.push(f1);
    fmap.features.push(f2);
    fmap
}

fn rich_library() -> LightTargetedExperiment {
    let compound =
        |id: &str, sequence: &str, charge: i32, rt: f64, proteins: &[&str]| LightCompound {
            id: id.into(),
            sequence: sequence.into(),
            charge,
            rt,
            protein_refs: proteins.iter().map(|p| p.to_string()).collect(),
            ..LightCompound::default()
        };
    let transition = |name: &str,
                      reference: &str,
                      pmz: f64,
                      mz: f64,
                      charge: i8,
                      decoy: bool,
                      detecting: bool| {
        let mut t = LightTransition {
            transition_name: name.into(),
            peptide_ref: reference.into(),
            precursor_mz: pmz,
            product_mz: mz,
            fragment_charge: charge,
            fragment_nr: 3,
            library_intensity: 10.0,
            ..LightTransition::default()
        };
        t.set_fragment_type("y");
        t.set_decoy(decoy);
        t.set_detecting_transition(detecting);
        t
    };
    let protein = |id: &str| LightProtein {
        id: id.into(),
        ..LightProtein::default()
    };
    LightTargetedExperiment {
        compounds: vec![
            compound("0", "PEPTIDEK", 2, 100.5, &["P1"]),
            compound("5", "ELVISLIVESK", 3, 250.25, &["P2"]),
            compound("9", "KEDITPEP", 2, 101.0, &["DECOY_P1"]),
        ],
        transitions: vec![
            transition("2", "0", 450.25, 300.1, 1, false, true),
            transition("4", "0", 450.25, 400.2, 2, false, true),
            transition("7", "0", 450.25, 500.3, 1, false, false),
            transition("11", "5", 600.5, 350.3, 1, false, true),
            transition("12", "5", 600.5, 360.3, 3, false, true),
            transition("20", "9", 450.5, 300.2, 1, true, true),
            transition("0", "9", 450.5, 410.2, 2, true, true),
        ],
        proteins: vec![protein("P1"), protein("P2"), protein("DECOY_P1")],
    }
}

fn rich_ids() -> SourceIDMapping {
    let mut ids = SourceIDMapping::default();
    ids.precursor_canonical_to_source
        .insert("5".into(), "ELVIS_3".into());
    ids.transition_canonical_to_source
        .insert("20".into(), "decoy_tr".into());
    ids
}

/// The driver's `richFeatures()`, value type for value type.
fn rich_features() -> FeatureMap {
    let mut fmap = FeatureMap::default();
    let s = |v: &str| MetaValue::from(v);
    let i = |v: i64| MetaValue::from(v);
    {
        let mut a = feature((1u64 << 63) | 5, 100.25, 1.1f32);
        set(&mut a, "PeptideRef", s("0"));
        set(&mut a, "im_drift", f(0.1));
        set(&mut a, "norm_RT", f(1.0 / 3.0));
        set(&mut a, "delta_rt", f(-0.0));
        set(&mut a, "leftWidth", f(99.5));
        set(&mut a, "rightWidth", s("101.25"));
        set(&mut a, "im_drift_left", s(" 7 "));
        set(&mut a, "im_drift_right", s("abc"));
        set(&mut a, "var_ms1_ppm_diff", f(1.5));
        set(&mut a, "ms1_area_intensity", f(123456.789012345678));
        set(&mut a, "ms1_apex_intensity", f(1e-300));
        set(&mut a, "im_ms1_drift", f(5e-324));
        set(&mut a, "im_ms1_delta", f(1e308));
        set(&mut a, "var_im_ms1_delta_score", i(42));
        set(&mut a, "var_ms1_mi_score", s("nan"));
        set(&mut a, "var_ms1_mi_contrast_score", s("inf"));
        set(&mut a, "var_ms1_mi_combined_score", s("+2.5e3"));
        set(&mut a, "var_ms1_isotope_correlation", fl(&[1.0, 2.0]));
        set(&mut a, "var_ms1_isotope_overlap", s(""));
        set(&mut a, "var_ms1_xcorr_coelution", f(0.30000000000000004));
        set(&mut a, "var_ms1_xcorr_shape", f(2.675));
        set(&mut a, "total_xic", f(1234567890123.0));
        set(&mut a, "peak_apices_sum", i(9007199254740993));
        set(&mut a, "var_bseries_score", i(3));
        set(&mut a, "var_library_corr", f(0.987654321987654321));
        set(&mut a, "var_manhatt_score", f(12.0));
        set(&mut a, "var_xcorr_shape", f(0.8));
        set(&mut a, "var_xcorr_shape_weighted", f(0.1 + 0.2));
        set(&mut a, "im_log_intensity", f(7.25));
        set(&mut a, "masserror_ppm", fl(&[1.5, -2.25]));
        set(&mut a, "id_target_transition_names", sl(&["2", "999", "4"]));
        set(&mut a, "id_target_area_intensity", fl(&[10.5, 11.5, 12.5]));
        set(
            &mut a,
            "id_target_total_area_intensity",
            sl(&["1", "x", "3"]),
        );
        set(&mut a, "id_target_apex_intensity", il(&[7, 8]));
        set(&mut a, "id_target_peak_apex_position", f(5.5));
        set(
            &mut a,
            "id_target_ind_massdev_score",
            fl(&[0.25, 0.5, 0.75]),
        );
        set(
            &mut a,
            "id_target_ind_log_intensity",
            fl(&[1.0 / 7.0, 2.0, 3.0]),
        );
        set(
            &mut a,
            "id_target_ind_start_position_at_5",
            sl(&["0", "1", "2"]),
        );
        set(
            &mut a,
            "id_target_ind_end_position_at_5",
            fl(&[9.0, 9.0, 9.0]),
        );
        set(&mut a, "id_decoy_transition_names", sl(&["20"]));
        set(&mut a, "id_decoy_area_intensity", fl(&[99.0]));
        set(&mut a, "id_decoy_ind_start_position_at_5", fl(&[3.5]));
        set(&mut a, "id_decoy_ind_end_position_at_5", fl(&[4.5]));
        set(&mut a, "id_decoy_ind_total_width", sl(&["1e-5"]));
        let mut s1 = sub("MS2", Some(s("2")), 10.0);
        set(&mut s1, "total_xic", f(20.0));
        set(&mut s1, "peak_apex_int", f(5.0));
        set(&mut s1, "peak_apex_position", f(100.3));
        set(&mut s1, "width_at_50", f(2.2));
        set(&mut s1, "total_mi", f(0.7));
        set(&mut s1, "start_position_at_5", f(99.0));
        set(&mut s1, "end_position_at_5", f(101.0));
        set(&mut s1, "tailing_factor", s("1.25"));
        set(&mut s1, "points_across_half_height", i(7));
        a.subordinates.push(s1);
        a.subordinates
            .push(sub("MS1", Some(s("0_Precursor_i0")), 50.0));
        a.subordinates.push(sub("MS2", Some(s("+4")), 11.0));
        a.subordinates.push(sub("MS2", None, 12.0));
        a.subordinates
            .push(sub("MS1", Some(s("0_Precursor_i1")), 0.0));
        let mut iso = sub("MS1", Some(s("0_Precursor_i2")), 30.5);
        set(&mut iso, "peak_apex_int", f(3.25));
        a.subordinates.push(iso);
        a.subordinates.push(sub("MS1", Some(s("noiso")), 30.0));
        a.subordinates
            .push(sub("MS1", Some(s("x_Precursor_iabc")), 30.0));
        a.subordinates
            .push(sub("MS1", Some(s("Precursor_i2_Precursor_i3")), 30.0));
        a.subordinates
            .push(sub("MS1", Some(s("Precursor_i 5 ")), 30.0));
        a.subordinates.push(sub("MS1", None, 30.0));
        a.subordinates.push(sub("ms2", Some(s("7")), 13.0));
        a.subordinates.push(sub("MS2", Some(i(7)), 14.0));
        fmap.features.push(a);
    }
    {
        let mut b = feature(77, 250.0, 0.0);
        set(&mut b, "PeptideRef", i(5));
        set(&mut b, "norm_RT", s("nan"));
        set(&mut b, "masserror_ppm", s("3.5"));
        set(&mut b, "var_norm_rt_score", il(&[1]));
        b.subordinates.push(sub("MS2", Some(s("11")), 1.0));
        b.subordinates.push(sub("MS2", Some(s("12")), 2.0));
        fmap.features.push(b);
    }
    {
        let mut c = feature(3, 101.0, 5.0);
        set(&mut c, "PeptideRef", s("9"));
        fmap.features.push(c);
    }
    fmap
}

fn second_run_features() -> FeatureMap {
    let mut fmap = FeatureMap::default();
    let mut f1 = feature(1000, 10.0, 7.5);
    set(&mut f1, "PeptideRef", MetaValue::from("5"));
    set(&mut f1, "var_xcorr_shape", f(0.5));
    f1.subordinates
        .push(sub("MS2", Some(MetaValue::from("11")), 3.0));
    fmap.features.push(f1);
    fmap
}

fn one_feature(peptide_ref: Option<&str>, native_id: Option<&str>) -> FeatureMap {
    let mut fmap = FeatureMap::default();
    let mut f1 = feature(1, 0.0, 0.0);
    if let Some(r) = peptide_ref {
        set(&mut f1, "PeptideRef", MetaValue::from(r));
    }
    if let Some(n) = native_id {
        f1.subordinates
            .push(sub("MS2", Some(MetaValue::from(n)), 1.0));
    }
    fmap.features.push(f1);
    fmap
}

// ------------------------------------------------------------------ class tests

/// `OpenSwathOSWParquetWriter_test.cpp`. The section's last part, a
/// cross-check against `MSChromatogramParquetConsumer`'s XIC output, needs
/// that class, which is not ported.
#[test]
fn class_test_writer() {
    let tmp = tmp();
    let base = tmp.path().join("oswpq");
    std::fs::create_dir_all(&base).unwrap();
    let lib = class_library();
    let fmap = class_features();
    let writer = OpenSwathOSWParquetWriter::new();
    writer
        .write(&base, &lib, &fmap, 90, "input.mzML", false)
        .unwrap();
    let lib_prec = base.join("library/precursors.parquet");
    let run_features = base.join("runs/run_id=90/features.parquet");
    assert!(lib_prec.exists());
    assert!(run_features.exists());

    let ids = |path: &Path, column: &str, default: i64| -> Vec<i64> {
        let table = ParquetFile::read_table(path).unwrap();
        let array = ParquetFile::get_column(&table, column).unwrap();
        (0..table.num_rows())
            .map(|r| ParquetFile::get_int64(Some(array.as_ref()), r, default, false).unwrap())
            .collect()
    };
    let lib_ids: HashSet<i64> = ids(&lib_prec, "precursor_id", 0).into_iter().collect();
    for pid in ids(&run_features, "precursor_id", 0) {
        assert!(lib_ids.contains(&pid));
    }
    assert_eq!(lib_ids.len(), 2);
    assert!(lib_ids.contains(&0) && lib_ids.contains(&7));

    let lib_trans = base.join("library/transitions.parquet");
    let map: BTreeMap<i64, i64> = ids(&lib_trans, "transition_id", 0)
        .into_iter()
        .zip(ids(&lib_trans, "precursor_id", 0))
        .collect();
    assert_eq!(map.len(), 2);
    assert_eq!(map[&3], 0);
    assert_eq!(map[&100], 7);

    let ft_path = base.join("runs/run_id=90/feature_transition.parquet");
    assert!(ft_path.exists());
    let ft_ids: HashSet<i64> = ids(&ft_path, "transition_id", -1).into_iter().collect();
    assert_eq!(ft_ids.len(), 2);
    assert!(ft_ids.contains(&3) && ft_ids.contains(&100));

    // An existing archive with a different canonical id domain is refused.
    let mut incompatible = lib.clone();
    incompatible.compounds[0].id = "1".into();
    incompatible.compounds[1].id = "8".into();
    incompatible.transitions[0].transition_name = "4".into();
    incompatible.transitions[0].peptide_ref = "1".into();
    incompatible.transitions[1].transition_name = "101".into();
    incompatible.transitions[1].peptide_ref = "8".into();
    OpenSwathLibraryIDNormalizer::validate_canonical_ids(&incompatible).unwrap();
    let incompatible_base = tmp.path().join("incompatible_oswpq");
    std::fs::create_dir_all(&incompatible_base).unwrap();
    TransitionParquetFile::new()
        .convert_light_targeted_experiment_to_parquet(&incompatible_base, &incompatible)
        .unwrap();
    let err = writer
        .write(&incompatible_base, &lib, &fmap, 91, "input.mzML", false)
        .unwrap_err();
    assert!(matches!(err, Error::InvalidValue(_)), "{err:?}");
}

fn write_cols(path: &Path, columns: Vec<(&str, ArrayRef)>) {
    let fields: Vec<Field> = columns
        .iter()
        .map(|(name, array)| Field::new(*name, array.data_type().clone(), true))
        .collect();
    let arrays: Vec<ArrayRef> = columns.into_iter().map(|(_, a)| a).collect();
    let table = RecordBatch::try_new(Arc::new(Schema::new(fields)), arrays).unwrap();
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    ParquetFile::write_table(&table, path).unwrap();
}

fn i64s(values: &[i64]) -> ArrayRef {
    Arc::new(Int64Array::from(values.to_vec()))
}

fn i32s(values: &[i32]) -> ArrayRef {
    Arc::new(Int32Array::from(values.to_vec()))
}

fn f64s(values: &[f64]) -> ArrayRef {
    Arc::new(Float64Array::from(values.to_vec()))
}

fn bools(values: &[bool]) -> ArrayRef {
    Arc::new(BooleanArray::from(values.to_vec()))
}

/// `OpenSwathOSWParquetReader_test.cpp`.
#[test]
fn class_test_reader() {
    let tmp = tmp();
    let base = tmp.path().join("test.oswpq");
    write_cols(
        &base.join("library/precursors.parquet"),
        vec![
            ("precursor_id", i64s(&[1])),
            ("precursor_mz", f64s(&[500.0])),
            ("charge", i32s(&[2])),
            ("decoy", bools(&[false])),
        ],
    );
    write_cols(
        &base.join("library/transitions.parquet"),
        vec![
            ("transition_id", i64s(&[1])),
            ("precursor_id", i64s(&[1])),
            ("charge", i32s(&[1])),
            ("decoy", bools(&[false])),
            ("detecting", bools(&[true])),
        ],
    );
    write_cols(
        &base.join("runs/runs.parquet"),
        vec![("run_id", i64s(&[1]))],
    );
    write_cols(
        &base.join("runs/run_id=1/features.parquet"),
        vec![
            ("feature_id", i64s(&[1, 2])),
            ("run_id", i64s(&[1, 1])),
            ("precursor_id", i64s(&[1, 1])),
            ("exp_rt", f64s(&[10.0, 11.0])),
        ],
    );
    write_cols(
        &base.join("runs/run_id=1/feature_transition.parquet"),
        vec![
            ("feature_id", i64s(&[1, 2])),
            ("transition_id", i64s(&[1, 1])),
            ("area_intensity", f64s(&[100.0, 110.0])),
        ],
    );
    let reader = OpenSwathOSWParquetReader::from_path(&base).unwrap();
    let rows = reader.rows();
    assert_eq!(rows.len(), 2);
    assert!(rows[0].feature_id == 1 || rows[1].feature_id == 1);
    let result = reader.fetch_peak_group_features(&base, "ms2", "").unwrap();
    assert_eq!(result.feature_id.len(), 2);
    assert_eq!(result.precursor_charge.len(), 2);
    assert_eq!(result.precursor_charge[0], 2);
    let t = reader.fetch_transition_features(&base).unwrap();
    assert_eq!(t.feature_id.len(), 2);
    assert_eq!(t.transition_id.len(), 2);
    assert_eq!(t.product_charge.len(), 2);
    assert_eq!(t.product_charge[0], 1);
    assert_eq!(t.precursor_charge[0], 2);
    assert_eq!(t.group_id[0], "1_1_1_1");
    let u = reader.fetch_unscored_data(&base).unwrap();
    assert_eq!(u.id.len(), 2);
    assert_eq!(u.rt.len(), 2);

    // The upstream OpenSwathWorkflow output.
    let real = Path::new(DATA).join(TWORUNS);
    let real_reader = OpenSwathOSWParquetReader::from_path(&real).unwrap();
    let pf = real_reader
        .fetch_peak_group_features(&real, "ms2", "")
        .unwrap();
    assert!(!pf.feature_id.is_empty());
    assert_eq!(pf.feature_id.len(), pf.exp_rt.len());
    for i in 0..pf.feature_id.len().min(5) {
        assert!(pf.feature_id[i] > 0);
        assert!(pf.exp_rt[i].is_finite());
    }
    let tf = real_reader.fetch_transition_features(&real).unwrap();
    assert!(!tf.feature_id.is_empty());
    assert_eq!(tf.feature_id.len(), tf.area_intensity.len());
    for i in 0..tf.feature_id.len().min(5) {
        assert!(tf.area_intensity[i].is_finite());
        assert!(tf.area_intensity[i] >= 0.0);
    }
    let uf = real_reader.fetch_unscored_data(&real).unwrap();
    assert!(!uf.id.is_empty());
    assert_eq!(uf.id.len(), uf.rt.len());
    if uf.left_width.len() == uf.id.len() {
        for i in 0..uf.id.len().min(5) {
            assert!(uf.left_width[i].is_finite());
            assert!(uf.right_width[i].is_finite());
        }
    }
    assert!(pf.feature_id.len() > 1 || tf.feature_id.len() > 1 || uf.id.len() > 1);
}

/// `OpenSwathOSWParquetRoundTrip_test.cpp`. The class test builds its
/// library with `convertTargetedExp`, which is not ported; the light
/// experiment the Release build produced from the same TraML file (record
/// `traml_input` of the TransitionParquetFile oracle) stands in. Its check
/// that `ZipRandomAccessFile::Open` succeeds on the archive has no
/// counterpart (that class is not ported; archive entries are extracted).
#[test]
fn class_test_roundtrip() {
    let all = read_records(TPF_DATA);
    let mut light = LightTargetedExperiment::default();
    for r in all.iter().filter(|r| r[0] == "C" && r[1] == "traml_input") {
        light.compounds.push(LightCompound {
            id: r[3].clone(),
            rt: num(&r[4]),
            drift_time: num(&r[5]),
            charge: r[6].parse().unwrap(),
            sequence: r[7].clone(),
            compound_name: r[9].clone(),
            sum_formula: r[15].clone(),
            ..LightCompound::default()
        });
    }
    for r in all.iter().filter(|r| r[0] == "T" && r[1] == "traml_input") {
        let mut t = LightTransition {
            transition_name: r[3].clone(),
            peptide_ref: r[4].clone(),
            library_intensity: num(&r[5]),
            product_mz: num(&r[6]),
            precursor_mz: num(&r[7]),
            precursor_im: num(&r[8]),
            fragment_charge: r[9].parse().unwrap(),
            fragment_nr: r[14].parse().unwrap(),
            ..LightTransition::default()
        };
        t.set_decoy(r[10] == "1");
        t.set_detecting_transition(r[11] == "1");
        t.set_identifying_transition(r[12] == "1");
        t.set_quantifying_transition(r[13] == "1");
        t.set_fragment_type(&r[15]);
        light.transitions.push(t);
    }
    assert!(!light.compounds.is_empty());

    let tmp = tmp();
    let out_archive = tmp.path().join("roundtrip.oswpq");
    let source_ids = OpenSwathLibraryIDNormalizer::normalize_source_ids(&mut light).unwrap();
    OpenSwathOSWParquetWriter::new()
        .write_with_source_ids(
            &out_archive,
            &light,
            &FeatureMap::default(),
            1,
            "test_input",
            false,
            Some(&source_ids),
        )
        .unwrap();
    assert!(out_archive.exists());
    let entries = ZipArchiveFile::list_entries(&out_archive);
    assert!(entries.iter().any(|e| e == "roundtrip.oswpq.idx.json"));

    let mut roundtrip = LightTargetedExperiment::default();
    let mut roundtrip_ids = SourceIDMapping::default();
    TransitionParquetFile::new()
        .convert_parquet_to_targeted_experiment_with_source_ids(
            &out_archive,
            &mut roundtrip,
            Some(&mut roundtrip_ids),
        )
        .unwrap();
    OpenSwathLibraryIDNormalizer::validate_canonical_ids(&roundtrip).unwrap();
    assert_eq!(
        roundtrip_ids.precursor_source_to_canonical.len(),
        source_ids.precursor_source_to_canonical.len()
    );
    assert_eq!(roundtrip.compounds.len(), light.compounds.len());
    assert_eq!(roundtrip.transitions.len(), light.transitions.len());
}

// ------------------------------------------------------------------ oracle: writer

#[test]
fn oracle_writer_class_and_rich_directories() {
    let all = records();
    let tmp = tmp();
    let writer = OpenSwathOSWParquetWriter::new();

    assert_eq!(outcome(&all, "w_class")[0], "OK");
    let dir = tmp.path().join("w_class.oswpq");
    std::fs::create_dir_all(&dir).unwrap();
    writer
        .write(
            &dir,
            &class_library(),
            &class_features(),
            90,
            "input.mzML",
            false,
        )
        .unwrap();
    assert!(assert_container_matches(&all, "w_class", &dir).is_empty());

    let rich = rich_library();
    let dir = tmp.path().join("w_rich_dir.oswpq");
    std::fs::create_dir_all(&dir).unwrap();
    assert_eq!(outcome(&all, "w_rich_dir")[0], "OK");
    let warnings = writer
        .write_with_source_ids(
            &dir,
            &rich,
            &rich_features(),
            (1u64 << 63) | 42,
            "run\t\"one\".mzML",
            true,
            Some(&rich_ids()),
        )
        .unwrap();
    assert!(warnings.is_empty());
    assert_all_stale(
        "w_rich_dir",
        assert_container_matches(&all, "w_rich_dir", &dir),
        &[42],
    );

    assert_eq!(outcome(&all, "w_rich_dir_append")[0], "OK");
    writer
        .write(&dir, &rich, &second_run_features(), 7, "second.mzML", false)
        .unwrap();
    assert_all_stale(
        "w_rich_dir_append",
        assert_container_matches(&all, "w_rich_dir_append", &dir),
        &[42, 7],
    );

    let err = writer
        .write(&dir, &rich, &second_run_features(), 7, "again.mzML", false)
        .unwrap_err();
    assert_error_matches("w_rich_dir_again", &err, outcome(&all, "w_rich_dir_again"));

    let dir = tmp.path().join("w_empty_fmap.oswpq");
    std::fs::create_dir_all(&dir).unwrap();
    writer
        .write(&dir, &rich, &FeatureMap::default(), 1, "", false)
        .unwrap();
    assert!(assert_container_matches(&all, "w_empty_fmap", &dir).is_empty());
}

fn unzip(tmp: &Path, zip: &Path, name: &str) -> PathBuf {
    let dir = tmp.join(name);
    for entry in ZipArchiveFile::list_entries(zip) {
        if entry.ends_with('/') {
            continue;
        }
        let mut guard = None;
        let path = ZipArchiveFile::extract_entry_to_temp_file(zip, &entry, &mut guard).unwrap();
        let target = dir.join(&entry);
        std::fs::create_dir_all(target.parent().unwrap()).unwrap();
        std::fs::copy(path, target).unwrap();
    }
    dir
}

fn zip_entries(all: &[Vec<String>], case: &str) -> BTreeSet<String> {
    of(all, "Z", case)
        .into_iter()
        .map(|r| r[3].clone())
        .collect()
}

/// Archives: new, extended with `preserve_existing`, and replaced. The
/// entry set is the Release build's; the order is sorted here (see the
/// support document), where the Release build's follows the file system.
#[test]
fn oracle_writer_archives() {
    let all = records();
    let tmp = tmp();
    let rich = rich_library();
    let zip = tmp.path().join("w_rich.oswpq");

    OpenSwathOSWParquetWriter::new()
        .write(&zip, &rich, &rich_features(), 42, "zip.mzML", true)
        .unwrap();
    let entries: BTreeSet<String> = ZipArchiveFile::list_entries(&zip).into_iter().collect();
    assert_eq!(entries, zip_entries(&all, "w_rich_zip"));
    let used = assert_container_matches(&all, "w_rich_zip", &unzip(tmp.path(), &zip, "first"));
    assert_all_stale("w_rich_zip", used, &[42]);

    let mut preserving = OpenSwathOSWParquetWriter::new();
    preserving.set_preserve_existing(true);
    assert!(preserving.preserve_existing());
    let warnings = preserving
        .write(&zip, &rich, &second_run_features(), 7, "second.mzML", false)
        .unwrap();
    assert_eq!(warnings.len(), 1);
    assert!(
        warnings[0].starts_with("Existing archive '"),
        "{warnings:?}"
    );
    let entries: BTreeSet<String> = ZipArchiveFile::list_entries(&zip).into_iter().collect();
    assert_eq!(entries, zip_entries(&all, "w_rich_zip_preserve"));
    let used = assert_container_matches(
        &all,
        "w_rich_zip_preserve",
        &unzip(tmp.path(), &zip, "preserved"),
    );
    assert_all_stale("w_rich_zip_preserve", used, &[42, 7]);

    OpenSwathOSWParquetWriter::new()
        .write(&zip, &rich, &second_run_features(), 8, "third.mzML", false)
        .unwrap();
    let entries: BTreeSet<String> = ZipArchiveFile::list_entries(&zip).into_iter().collect();
    assert_eq!(entries, zip_entries(&all, "w_rich_zip_replace"));
    let used = assert_container_matches(
        &all,
        "w_rich_zip_replace",
        &unzip(tmp.path(), &zip, "replaced"),
    );
    assert_all_stale("w_rich_zip_replace", used, &[8]);
}

#[test]
fn oracle_writer_errors() {
    let all = records();
    let tmp = tmp();
    let rich = rich_library();
    let writer = OpenSwathOSWParquetWriter::new();

    // Errors while the features are converted leave the run directory behind.
    for (case, fmap) in [
        ("w_no_peptideref", one_feature(None, None)),
        ("w_unknown_peptide", one_feature(Some("12"), None)),
        ("w_unknown_numeric", one_feature(Some("0"), Some("12345"))),
        ("w_unknown_id", one_feature(Some("0"), Some("abc"))),
    ] {
        let dir = tmp.path().join(format!("{case}.oswpq"));
        std::fs::create_dir_all(&dir).unwrap();
        let err = writer
            .write(&dir, &rich, &fmap, 3, "x.mzML", false)
            .unwrap_err();
        assert_error_matches(case, &err, outcome(&all, case));
        assert_eq!(
            dir.join("runs/run_id=3").exists(),
            fact(&all, case, "run_dir_exists") == "1",
            "{case}"
        );
        assert_eq!(
            dir.join("runs/runs.parquet").exists(),
            fact(&all, case, "runs_exists") == "1",
            "{case}"
        );
    }

    // A run id already in runs.parquet whose directory is gone.
    let dir = tmp.path().join("w_runs_dup.oswpq");
    std::fs::create_dir_all(&dir).unwrap();
    writer
        .write(&dir, &rich, &second_run_features(), 5, "a.mzML", false)
        .unwrap();
    std::fs::remove_dir_all(dir.join("runs/run_id=5")).unwrap();
    let err = writer
        .write(&dir, &rich, &second_run_features(), 5, "b.mzML", false)
        .unwrap_err();
    assert_error_matches("w_runs_dup", &err, outcome(&all, "w_runs_dup"));
    assert_eq!(
        dir.join("runs/run_id=5").exists(),
        fact(&all, "w_runs_dup", "run_dir_exists") == "1"
    );

    // An existing library with duplicated ids.
    let dir = tmp.path().join("w_existing_dup.oswpq");
    write_cols(
        &dir.join("library/precursors.parquet"),
        vec![("precursor_id", i64s(&[0, 0, 5, 9]))],
    );
    write_cols(
        &dir.join("library/transitions.parquet"),
        vec![("transition_id", i64s(&[2])), ("precursor_id", i64s(&[0]))],
    );
    let err = writer
        .write(&dir, &rich, &second_run_features(), 5, "a.mzML", false)
        .unwrap_err();
    assert_error_matches(
        "w_existing_dup_precursor",
        &err,
        outcome(&all, "w_existing_dup_precursor"),
    );
    let dir = tmp.path().join("w_existing_dup_tr.oswpq");
    write_cols(
        &dir.join("library/precursors.parquet"),
        vec![("precursor_id", i64s(&[0, 5, 9]))],
    );
    write_cols(
        &dir.join("library/transitions.parquet"),
        vec![
            ("transition_id", i64s(&[2, 2])),
            ("precursor_id", i64s(&[0, 0])),
        ],
    );
    let err = writer
        .write(&dir, &rich, &second_run_features(), 5, "a.mzML", false)
        .unwrap_err();
    assert_error_matches(
        "w_existing_dup_transition",
        &err,
        outcome(&all, "w_existing_dup_transition"),
    );

    // A different canonical id domain (the class test's incompatible library).
    let lib = class_library();
    let mut inc = lib.clone();
    inc.compounds[0].id = "1".into();
    inc.compounds[1].id = "8".into();
    inc.transitions[0].transition_name = "4".into();
    inc.transitions[0].peptide_ref = "1".into();
    inc.transitions[1].transition_name = "101".into();
    inc.transitions[1].peptide_ref = "8".into();
    let base = tmp.path().join("incompatible.oswpq");
    std::fs::create_dir_all(&base).unwrap();
    TransitionParquetFile::new()
        .convert_light_targeted_experiment_to_parquet(&base, &inc)
        .unwrap();
    let err = writer
        .write(&base, &lib, &class_features(), 91, "input.mzML", false)
        .unwrap_err();
    assert_error_matches(
        "w_class_incompatible",
        &err,
        outcome(&all, "w_class_incompatible"),
    );
    assert_eq!(
        base.join("runs/run_id=91").exists(),
        fact(&all, "w_class_incompatible", "run_dir_exists") == "1"
    );

    // Non-canonical library ids.
    let mut lib = rich.clone();
    lib.compounds[0].id = "PEP".into();
    let dir = tmp.path().join("w_noncanonical.oswpq");
    std::fs::create_dir_all(&dir).unwrap();
    let err = writer
        .write(&dir, &lib, &FeatureMap::default(), 1, "", false)
        .unwrap_err();
    assert_error_matches("w_noncanonical", &err, outcome(&all, "w_noncanonical"));
}

// ------------------------------------------------------------------ oracle: reader

fn check_reader(all: &[Vec<String>], case: &str, path: &Path) {
    let reader_case = |suffix: &str| format!("{case}/{suffix}");

    // load
    let c = reader_case("load");
    let mut reader = OpenSwathOSWParquetReader::new();
    let result = reader.load(path);
    match outcome(all, &c)[0].as_str() {
        "OK" => {
            result.unwrap_or_else(|e| panic!("{c}: {e:?}"));
            assert_eq!(reader.oswpq_path(), path);
            let want = of(all, "L", case);
            assert_eq!(reader.rows().len(), want.len(), "{c}");
            for (row, w) in reader.rows().iter().zip(want) {
                let ctx = format!("{c} row {}", w[2]);
                assert_eq!(row.feature_id.to_string(), w[3], "{ctx}");
                assert_eq!(row.run_id.to_string(), w[4], "{ctx}");
                assert_eq!(row.precursor_id.to_string(), w[5], "{ctx}");
                assert_same(&ctx, row.exp_rt, &w[6]);
                assert_same(&ctx, row.ms2_area_intensity, &w[7]);
                assert_same(&ctx, row.ms2_total_area_intensity, &w[8]);
                assert_same(&ctx, row.ms2_apex_intensity, &w[9]);
                assert_eq!(row.precursor_charge.to_string(), w[10], "{ctx}");
                assert_eq!(u8::from(row.decoy).to_string(), w[11], "{ctx}");
                assert_eq!(row.transition_count.to_string(), w[12], "{ctx}");
                assert_eq!(row.group_id, w[13], "{ctx}");
            }
        }
        _ => {
            assert_error_matches(&c, &result.unwrap_err(), outcome(all, &c));
            assert!(reader.rows().is_empty());
        }
    }

    // fetch_peak_group_features
    for (tag, level, main_score) in [
        ("ms2", "ms2", ""),
        ("ms1ms2", "ms1ms2", ""),
        ("main", "ms2", "var_ms2_xcorr_shape"),
        ("main_missing", "ms2", "no_such_score"),
    ] {
        let c = reader_case(&format!("pg_{tag}"));
        let result = reader.fetch_peak_group_features(path, level, main_score);
        if outcome(all, &c)[0] != "OK" {
            assert_error_matches(&c, &result.unwrap_err(), outcome(all, &c));
            continue;
        }
        let r = result.unwrap_or_else(|e| panic!("{c}: {e:?}"));
        let columns = |kind: &str| -> Vec<String> {
            of(all, "PGC", &c)
                .into_iter()
                .filter(|x| x[2] == kind)
                .map(|x| x[4].clone())
                .collect()
        };
        assert_eq!(r.ms2_columns, columns("ms2"), "{c}");
        assert_eq!(r.ms1_columns, columns("ms1"), "{c}");
        let want = of(all, "PG", &c);
        assert_eq!(r.feature_id.len(), want.len(), "{c}");
        for (i, w) in want.iter().enumerate() {
            let ctx = format!("{c} row {i}");
            assert_eq!(r.feature_id[i].to_string(), w[3], "{ctx}");
            assert_eq!(r.run_id[i].to_string(), w[4], "{ctx}");
            assert_eq!(r.precursor_id[i].to_string(), w[5], "{ctx}");
            assert_same(&ctx, r.exp_rt[i], &w[6]);
            assert_eq!(r.precursor_charge[i].to_string(), w[7], "{ctx}");
            assert_eq!(u8::from(r.decoy[i]).to_string(), w[8], "{ctx}");
            assert_eq!(r.transition_count[i].to_string(), w[9], "{ctx}");
            assert_eq!(r.group_id[i], w[10], "{ctx}");
            assert_values(&ctx, &r.ms2_values, i, &w[11]);
            assert_values(&ctx, &r.ms1_values, i, &w[12]);
        }
        // The source warns (OPENMS_LOG_WARN) when a requested main score is
        // not among the discovered MS2 columns, and only when there are any.
        let missing = !main_score.is_empty()
            && !r.ms2_columns.is_empty()
            && !r.ms2_columns.iter().any(|x| x == main_score);
        assert_eq!(r.warnings.len(), usize::from(missing), "{c}");
    }

    // fetch_transition_features
    let c = reader_case("tf");
    let result = reader.fetch_transition_features(path);
    if outcome(all, &c)[0] == "OK" {
        let r = result.unwrap_or_else(|e| panic!("{c}: {e:?}"));
        let columns: Vec<String> = of(all, "TFC", &c)
            .into_iter()
            .map(|x| x[3].clone())
            .collect();
        assert_eq!(r.transition_var_columns, columns, "{c}");
        let want = of(all, "TF", &c);
        assert_eq!(r.feature_id.len(), want.len(), "{c}");
        for (i, w) in want.iter().enumerate() {
            let ctx = format!("{c} row {i}");
            assert_eq!(r.feature_id[i].to_string(), w[3], "{ctx}");
            assert_eq!(r.run_id[i].to_string(), w[4], "{ctx}");
            assert_eq!(r.precursor_id[i].to_string(), w[5], "{ctx}");
            assert_same(&ctx, r.exp_rt[i], &w[6]);
            assert_eq!(r.precursor_charge[i].to_string(), w[7], "{ctx}");
            assert_eq!(r.transition_id[i].to_string(), w[8], "{ctx}");
            assert_eq!(r.product_charge[i].to_string(), w[9], "{ctx}");
            assert_eq!(u8::from(r.decoy[i]).to_string(), w[10], "{ctx}");
            assert_same(&ctx, r.area_intensity[i], &w[11]);
            assert_same(&ctx, r.total_area_intensity[i], &w[12]);
            assert_same(&ctx, r.apex_intensity[i], &w[13]);
            assert_same(&ctx, r.apex_rt[i], &w[14]);
            assert_same(&ctx, r.rt_fwhm[i], &w[15]);
            assert_same(&ctx, r.masserror_ppm[i], &w[16]);
            assert_same(&ctx, r.total_mi[i], &w[17]);
            assert_eq!(r.group_id[i], w[18], "{ctx}");
            assert_values(&ctx, &r.transition_var_values, i, &w[19]);
        }
    } else {
        assert_error_matches(&c, &result.unwrap_err(), outcome(all, &c));
    }

    // fetch_unscored_data
    let c = reader_case("u");
    let result = reader.fetch_unscored_data(path);
    if outcome(all, &c)[0] == "OK" {
        let r = result.unwrap_or_else(|e| panic!("{c}: {e:?}"));
        let columns = |kind: &str| -> Vec<String> {
            of(all, "UC", &c)
                .into_iter()
                .filter(|x| x[2] == kind)
                .map(|x| x[4].clone())
                .collect()
        };
        assert_eq!(r.ms2_columns, columns("ms2"), "{c}");
        assert_eq!(r.ms1_columns, columns("ms1"), "{c}");
        let want = of(all, "U", &c);
        assert_eq!(r.id.len(), want.len(), "{c}");
        for (i, w) in want.iter().enumerate() {
            let ctx = format!("{c} row {i}");
            assert_eq!(r.id_run[i].to_string(), w[3], "{ctx}");
            assert_eq!(r.id_peptide[i].to_string(), w[4], "{ctx}");
            assert_eq!(r.transition_group_id[i].to_string(), w[5], "{ctx}");
            assert_eq!(u8::from(r.decoy[i]).to_string(), w[6], "{ctx}");
            assert_eq!(r.run_id[i].to_string(), w[7], "{ctx}");
            assert_eq!(r.filename[i], w[8], "{ctx}");
            assert_same(&ctx, r.rt[i], &w[9]);
            assert_same(&ctx, r.assay_rt[i], &w[10]);
            assert_same(&ctx, r.delta_rt[i], &w[11]);
            assert_same(&ctx, r.assay_rt_library[i], &w[12]);
            assert_same(&ctx, r.delta_rt_library[i], &w[13]);
            assert_eq!(r.id[i].to_string(), w[14], "{ctx}");
            assert_eq!(r.charge[i].to_string(), w[15], "{ctx}");
            assert_same(&ctx, r.mz[i], &w[16]);
            assert_same(&ctx, r.intensity[i], &w[17]);
            assert_same(&ctx, r.aggr_prec_peak_area[i], &w[18]);
            assert_same(&ctx, r.aggr_prec_peak_apex[i], &w[19]);
            assert_same(&ctx, r.left_width[i], &w[20]);
            assert_same(&ctx, r.right_width[i], &w[21]);
            assert_same(&ctx, r.exp_im[i], &w[22]);
            assert_same(&ctx, r.im_left_width[i], &w[23]);
            assert_same(&ctx, r.im_right_width[i], &w[24]);
            assert_values(&ctx, &r.ms2_values, i, &w[25]);
            assert_values(&ctx, &r.ms1_values, i, &w[26]);
        }
    } else {
        assert_error_matches(&c, &result.unwrap_err(), outcome(all, &c));
    }
}

#[test]
fn oracle_reader_on_release_containers() {
    let all = records();
    let tmp = tmp();
    check_reader(&all, "r_tworuns", &Path::new(DATA).join(TWORUNS));
    let cases = [
        ("r_class", "w_class"),
        ("r_rich_dir", "w_rich_dir_append"),
        ("r_minimal", "r_minimal_files"),
        ("r_ties", "r_ties_files"),
        ("r_err_unknown_precursor", "r_err_unknown_precursor_files"),
        ("r_err_unknown_feature", "r_err_unknown_feature_files"),
        ("r_err_unknown_transition", "r_err_unknown_transition_files"),
        ("r_err_no_runs", "r_err_no_runs_files"),
    ];
    for (case, source) in cases {
        let dir = container(&all, tmp.path(), source, &format!("{case}.oswpq"));
        check_reader(&all, case, &dir);
    }
    // The same containers read through archives.
    for (case, source) in [
        ("r_ties_zip", "r_ties_files"),
        ("r_rich_zip_first", "w_rich_zip"),
        ("r_rich_zip_replaced", "w_rich_zip_replace"),
    ] {
        let dir = container(&all, tmp.path(), source, &format!("{case}_src"));
        let zip = tmp.path().join(format!("{case}.oswpq"));
        ZipArchiveFile::zip_directory(&dir, &zip).unwrap();
        check_reader(&all, case, &zip);
    }
}

/// Documented divergence: the Release reader cannot read the archive its own
/// writer extends with `preserve_existing` (measured: libzip rewrites every
/// replaced entry deflated, and `ZipRandomAccessFile` opens a deflated entry
/// and then fails to seek in it, "Failed to create parquet reader from
/// RandomAccessFile", instead of falling back to extraction). The port
/// extracts entries and reads the archive. The archive here starts as the
/// Release writer's first archive (`w_rich_zip`, run 42) and the port's
/// writer extends it with run 7; its runs are then those of the directory
/// case `r_rich_dir` (the Release writer's runs 42 and 7 of the same library
/// and features), whose Release load, peak-group and transition dumps it
/// matches. Starting from the Release run 42 keeps the Release writer's
/// stale values (see [`STALE`]) on both sides; run 7 is the port's, and its
/// one stale cell (`var_ms2_xcorr_shape`, its `STALE` entry) is checked
/// as modelled.
#[test]
fn oracle_reader_preserved_archive_divergence() {
    let all = records();
    for suffix in ["load", "pg_ms2", "tf", "u"] {
        let record = outcome(&all, &format!("r_rich_zip/{suffix}"));
        assert_eq!(record[0], "E");
        assert!(record[3].contains("Failed to create parquet reader from RandomAccessFile"));
    }
    let tmp = tmp();
    let rich = rich_library();
    let first = container(&all, tmp.path(), "w_rich_zip", "release_first");
    let zip = tmp.path().join("w_rich.oswpq");
    ZipArchiveFile::zip_directory(&first, &zip).unwrap();
    let mut preserving = OpenSwathOSWParquetWriter::new();
    preserving.set_preserve_existing(true);
    preserving
        .write(&zip, &rich, &second_run_features(), 7, "second.mzML", false)
        .unwrap();
    let reader = OpenSwathOSWParquetReader::from_path(&zip).unwrap();
    let want = of(&all, "L", "r_rich_dir");
    assert_eq!(reader.rows().len(), want.len());
    for (row, w) in reader.rows().iter().zip(want) {
        assert_eq!(row.group_id, w[13]);
        assert_same("preserved load", row.exp_rt, &w[6]);
    }
    let pg = reader
        .fetch_peak_group_features(&zip, "ms1ms2", "")
        .unwrap();
    let want = of(&all, "PG", "r_rich_dir/pg_ms1ms2");
    assert_eq!(pg.feature_id.len(), want.len());
    let stale_column = pg
        .ms2_columns
        .iter()
        .position(|c| c == "var_ms2_xcorr_shape")
        .unwrap();
    let run7 = STALE
        .iter()
        .position(|e| e.0 == "runs/run_id=7/features.parquet")
        .unwrap();
    let (release, port) = stale_pair(run7);
    let mut stale_rows = 0;
    for (i, w) in want.iter().enumerate() {
        assert_eq!(pg.group_id[i], w[10]);
        assert_values("preserved pg", &pg.ms1_values, i, &w[12]);
        let texts: Vec<&str> = w[11].split('|').collect();
        assert_eq!(pg.ms2_values.len(), texts.len());
        for (c, text) in texts.iter().enumerate() {
            if pg.run_id[i] == 7 && c == stale_column {
                stale_rows += 1;
                assert_same("preserved pg Release", release.unwrap(), text);
                assert_eq!(Some(pg.ms2_values[c][i]), port);
            } else {
                assert_same("preserved pg", pg.ms2_values[c][i], text);
            }
        }
    }
    assert_eq!(stale_rows, 1);
    let tf = reader.fetch_transition_features(&zip).unwrap();
    let want = of(&all, "TF", "r_rich_dir/tf");
    assert_eq!(tf.feature_id.len(), want.len());
    for (i, w) in want.iter().enumerate() {
        assert_eq!(tf.group_id[i], w[18]);
        assert_values("preserved tf", &tf.transition_var_values, i, &w[19]);
    }
}
