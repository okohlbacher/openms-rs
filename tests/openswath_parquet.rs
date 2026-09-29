// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

#![cfg(feature = "arrow")]

//! The OpenSWATH Parquet readers (`.xic`, `.xim`, `.xipm`) against the
//! upstream class tests.
//!
//! Every section of `XICParquetFile_test.cpp`, `XIMParquetFile_test.cpp` and
//! `XIPMParquetFile_test.cpp` at core `bc9cc12` is reproduced with its
//! literals (evidence tier 3). The `.xic` and `.xim` inputs are the retained
//! upstream fixtures `XICParquetFile_1_input.xic`, `XICParquetFile_2_input.xic`
//! and `XIMParquetFile_23_input.xim`, which the C++ writers produced. The
//! XIPM class test writes its input with `XIPMParquetConsumer`, which needs
//! the unported `PeakMapExtractor::ExtractedPeakMap`; its sections run here
//! on the retained `XIPMParquetFile_1_input.xipm`, which holds the same two
//! peak maps (checked field by field below with an independent reader's
//! listing, pyarrow 19) but predates the `TARGET_RT` column. Counts marked
//! "fixture listing" are that pyarrow listing, not Rust output. See
//! `tests/data/openswath_parquet_provenance.json`.

use arrow_array::builder::{
    BinaryBuilder, Float32Builder, Float64Builder, Int64Builder, LargeStringBuilder, StringBuilder,
    UInt64Builder,
};
use arrow_array::{ArrayRef, RecordBatch};
use arrow_schema::{DataType, Field, Schema};
use openms::Error;
use openms::format::parquet_file::ParquetFile;
use openms::format::parquet_filter::ParquetFilter;
use openms::format::xic_parquet_file::{XicParquetFile, XicQuery};
use openms::format::xim_parquet_file::{XimParquetFile, XimQuery};
use openms::format::xipm_parquet_file::{XipmParquetFile, XipmQuery};
use std::path::PathBuf;
use std::sync::Arc;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/openswath_parquet")
        .join(name)
}

fn scratch(test: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms-oswpq-{}-{test}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The class tests' `encodeDoubles_`: raw host-order (little-endian) doubles.
fn encode_doubles(values: &[f64]) -> Vec<u8> {
    values
        .iter()
        .flat_map(|value| value.to_le_bytes())
        .collect()
}

fn binary(values: &[f64]) -> ArrayRef {
    let mut builder = BinaryBuilder::new();
    builder.append_value(encode_doubles(values));
    Arc::new(builder.finish())
}

fn int64(value: i64) -> ArrayRef {
    let mut builder = Int64Builder::new();
    builder.append_value(value);
    Arc::new(builder.finish())
}

fn large_utf8(value: &str) -> ArrayRef {
    let mut builder = LargeStringBuilder::new();
    builder.append_value(value);
    Arc::new(builder.finish())
}

/// The class tests' `writeParquetTable_` (`parquet::arrow::WriteTable`
/// without properties), which `ParquetFile::write_table` ports.
fn write(path: &PathBuf, fields: Vec<(&str, DataType)>, columns: Vec<ArrayRef>) {
    let schema = Arc::new(Schema::new(
        fields
            .into_iter()
            .map(|(name, data_type)| Field::new(name, data_type, true))
            .collect::<Vec<_>>(),
    ));
    let table = RecordBatch::try_new(schema, columns).unwrap();
    ParquetFile::write_table(&table, path).unwrap();
}

fn is_invalid<T: std::fmt::Debug>(result: openms::Result<T>) -> bool {
    matches!(result, Err(Error::InvalidValue(_)))
}

fn is_not_found<T: std::fmt::Debug>(result: openms::Result<T>) -> bool {
    matches!(result, Err(Error::Io(ref e)) if e.kind() == std::io::ErrorKind::NotFound)
}

fn filter(text: &str) -> XicQuery {
    XicQuery {
        filter: text.into(),
        ..XicQuery::default()
    }
}

// ---------------------------------------------------------------- XIC ----

/// `START_SECTION(void load(std::vector<XICChromatogram>& output) const)`.
#[test]
fn xic_load() {
    let xic = XicParquetFile::new(data("XICParquetFile_1_input.xic")).unwrap();
    let chroms = xic.load().unwrap();
    assert_eq!(chroms.len(), 18);
    assert!(!chroms[0].rt.is_empty());
    assert!(!chroms[0].intensity.is_empty());
    for c in &chroms {
        assert_eq!(c.rt.len(), c.intensity.len());
    }
}

/// `START_SECTION(void getChromatograms(...) const)`.
#[test]
fn xic_get_chromatograms() {
    let xic = XicParquetFile::new(data("XICParquetFile_1_input.xic")).unwrap();
    let all = xic.chromatograms(&XicQuery::default()).unwrap();
    assert_eq!(all.len(), 18);
    assert!(
        xic.chromatograms(&filter("precursor_id=99999"))
            .unwrap()
            .is_empty()
    );
    let filtered = xic.chromatograms(&filter("precursor_id=2")).unwrap();
    assert!(!filtered.is_empty());
    assert!(filtered.len() <= all.len());
    // Fixture listing: precursor 2 owns transitions 1, 2 and 3.
    assert_eq!(filtered.len(), 3);
    let eqeq = xic.chromatograms(&filter("precursor_id==2")).unwrap();
    assert_eq!(eqeq.len(), filtered.len());
    let mut typed = ParquetFilter::new();
    typed.eq_int("PRECURSOR_ID", 2);
    let typed = xic.chromatograms_with_filter(&typed).unwrap();
    assert!(!typed.is_empty());
    assert!(typed.len() <= all.len());
    assert!(is_invalid(xic.chromatograms(&filter("precursor_id=="))));
    assert!(is_invalid(xic.chromatograms(&filter("precursor_id IN []"))));
    assert!(is_invalid(xic.chromatograms(&filter("precusor_id=2"))));
}

/// `START_SECTION(void getChromatograms_multi_file)`.
#[test]
fn xic_multi_file() {
    let xic = XicParquetFile::from_files(&[
        data("XICParquetFile_1_input.xic"),
        data("XICParquetFile_2_input.xic"),
    ])
    .unwrap();
    assert_eq!(xic.chromatograms(&XicQuery::default()).unwrap().len(), 36);
    assert_eq!(xic.filenames().len(), 2);
    assert_eq!(xic.filename(), data("XICParquetFile_1_input.xic"));
}

/// `START_SECTION(void getChromatograms_large_string_columns)`.
#[test]
fn xic_large_string_columns() {
    let path = scratch("xic-large").join("large.xic");
    write(
        &path,
        vec![
            ("RUN_ID", DataType::Int64),
            ("PRECURSOR_ID", DataType::Int64),
            ("SOURCE_FILE", DataType::LargeUtf8),
            ("MODIFIED_SEQUENCE", DataType::LargeUtf8),
            ("RT_DATA", DataType::Binary),
            ("INTENSITY_DATA", DataType::Binary),
            ("RT_COMPRESSION", DataType::Int64),
            ("INTENSITY_COMPRESSION", DataType::Int64),
        ],
        vec![
            int64(7),
            int64(42),
            large_utf8("large_source.raw"),
            large_utf8("PEPTIDE"),
            binary(&[100.0, 101.0]),
            binary(&[1000.0, 1001.0]),
            int64(0),
            int64(0),
        ],
    );
    let xic = XicParquetFile::new(&path).unwrap();
    let chroms = xic
        .chromatograms(&filter("MODIFIED_SEQUENCE=PEPTIDE"))
        .unwrap();
    assert_eq!(chroms.len(), 1);
    assert_eq!(chroms[0].source_file, "large_source.raw");
    assert_eq!(chroms[0].modified_sequence, "PEPTIDE");
    assert_eq!(chroms[0].rt.len(), 2);
    assert_eq!(chroms[0].intensity.len(), 2);
    assert_eq!(chroms[0].rt[0], 100.0);
    assert_eq!(chroms[0].intensity[1], 1001.0);
}

/// `START_SECTION(void load_uint64_overflow_throws)`.
#[test]
fn xic_load_uint64_overflow_throws() {
    let path = scratch("xic-u64").join("overflow.xic");
    let mut run_id = UInt64Builder::new();
    run_id.append_value(i64::MAX as u64 + 1);
    let mut source = StringBuilder::new();
    source.append_value("overflow.raw");
    write(
        &path,
        vec![
            ("RUN_ID", DataType::UInt64),
            ("PRECURSOR_ID", DataType::Int64),
            ("SOURCE_FILE", DataType::Utf8),
            ("RT_DATA", DataType::Binary),
            ("INTENSITY_DATA", DataType::Binary),
            ("RT_COMPRESSION", DataType::Int64),
            ("INTENSITY_COMPRESSION", DataType::Int64),
        ],
        vec![
            Arc::new(run_id.finish()),
            int64(42),
            Arc::new(source.finish()),
            binary(&[100.0]),
            binary(&[1000.0]),
            int64(0),
            int64(0),
        ],
    );
    let xic = XicParquetFile::new(&path).unwrap();
    assert!(is_invalid(xic.load()));
}

/// `START_SECTION(void getRuns(std::vector<XICRunInfo>& output) const)`.
#[test]
fn xic_get_runs() {
    let xic = XicParquetFile::new(data("XICParquetFile_1_input.xic")).unwrap();
    let runs = xic.runs().unwrap();
    assert_eq!(runs.len(), 1);
    assert_ne!(runs[0].run_id, 0);
}

/// `START_SECTION(void getAnalytes(std::vector<XICAnalyte>& output, bool) const)`.
#[test]
fn xic_get_analytes() {
    let xic = XicParquetFile::new(data("XICParquetFile_1_input.xic")).unwrap();
    let none: [&str; 0] = [];
    let exploded = xic.analytes(&none, false).unwrap();
    assert_eq!(exploded.len(), 18);
    let nested = xic.analytes(&none, true).unwrap();
    assert_eq!(nested.len(), 6);
    let transitions: usize = nested.iter().map(|a| a.transition_ids.len()).sum();
    assert_eq!(transitions, exploded.len());
}

/// `START_SECTION(void load_invalid_path)`.
#[test]
fn xic_load_invalid_path() {
    assert!(is_not_found(XicParquetFile::new("no_such_file.xic")));
}

// ---------------------------------------------------------------- XIM ----

fn xim_query(filter: &str) -> XimQuery {
    XimQuery {
        filter: filter.into(),
        ..XimQuery::default()
    }
}

/// `START_SECTION(void load(std::vector<XIMMobilogram>& output) const)`.
#[test]
fn xim_load() {
    let xim = XimParquetFile::new(data("XIMParquetFile_23_input.xim")).unwrap();
    let mobilograms = xim.load().unwrap();
    assert!(!mobilograms.is_empty());
    assert!(!mobilograms[0].mobility.is_empty());
    assert!(!mobilograms[0].intensity.is_empty());
    for m in &mobilograms {
        assert_eq!(m.mobility.len(), m.intensity.len());
    }
}

/// `START_SECTION(void getMobilograms(...) const)`.
#[test]
fn xim_get_mobilograms() {
    let xim = XimParquetFile::new(data("XIMParquetFile_23_input.xim")).unwrap();
    let all = xim.mobilograms(&XimQuery::default()).unwrap();
    assert!(!all.is_empty());
    assert!(
        xim.mobilograms(&xim_query("feature_id=99999"))
            .unwrap()
            .is_empty()
    );
    let by_run = XimQuery {
        run_id: Some(all[0].run_id),
        ..XimQuery::default()
    };
    let filtered_string = xim.mobilograms(&by_run).unwrap();
    assert!(!filtered_string.is_empty());
    assert!(filtered_string.len() <= all.len());
    let filtered_args = xim.mobilograms(&by_run).unwrap();
    assert!(!filtered_args.is_empty());
    assert!(filtered_args.len() <= all.len());
    let mut typed = ParquetFilter::new();
    typed.eq_int("RUN_ID", all[0].run_id);
    let typed = xim.mobilograms_with_filter(&typed).unwrap();
    assert!(!typed.is_empty());
    assert!(typed.len() <= all.len());
    assert!(is_invalid(xim.mobilograms(&xim_query("feature_id=="))));
}

/// `START_SECTION(void getMobilograms_multi_file)`.
#[test]
fn xim_multi_file() {
    let file = data("XIMParquetFile_23_input.xim");
    let single = XimParquetFile::new(&file).unwrap().load().unwrap();
    let multi = XimParquetFile::from_files(&[&file, &file])
        .unwrap()
        .load()
        .unwrap();
    assert_eq!(multi.len(), 2 * single.len());
}

/// `START_SECTION(void getMobilograms_large_string_and_double_feature_rt)`.
#[test]
fn xim_large_string_and_double_feature_rt() {
    let path = scratch("xim-large").join("large.xim");
    let mut feature_rt = Float64Builder::new();
    feature_rt.append_value(12.5);
    write(
        &path,
        vec![
            ("RUN_ID", DataType::Int64),
            ("PRECURSOR_ID", DataType::Int64),
            ("FEATURE_ID", DataType::Int64),
            ("FEATURE_RT", DataType::Float64),
            ("SOURCE_FILE", DataType::LargeUtf8),
            ("MOBILOGRAM_TYPE", DataType::LargeUtf8),
            ("MODIFIED_SEQUENCE", DataType::LargeUtf8),
            ("MOBILITY_DATA", DataType::Binary),
            ("INTENSITY_DATA", DataType::Binary),
            ("MOBILITY_COMPRESSION", DataType::Int64),
            ("INTENSITY_COMPRESSION", DataType::Int64),
        ],
        vec![
            int64(11),
            int64(101),
            int64(202),
            Arc::new(feature_rt.finish()),
            large_utf8("large_mobi.raw"),
            large_utf8("precursor"),
            large_utf8("PEPTIDE"),
            binary(&[0.8, 0.9]),
            binary(&[10.0, 20.0]),
            int64(0),
            int64(0),
        ],
    );
    let xim = XimParquetFile::new(&path).unwrap();
    let query = XimQuery {
        feature_rt: Some(12.5),
        filter: "MOBILOGRAM_TYPE=precursor".into(),
        ..XimQuery::default()
    };
    let mobilograms = xim.mobilograms(&query).unwrap();
    assert_eq!(mobilograms.len(), 1);
    assert_eq!(mobilograms[0].source_file, "large_mobi.raw");
    assert_eq!(mobilograms[0].mobilogram_type, "precursor");
    assert_eq!(mobilograms[0].feature_rt, Some(12.5));
    assert_eq!(mobilograms[0].mobility.len(), 2);
    assert_eq!(mobilograms[0].intensity.len(), 2);
}

/// `START_SECTION(void getMobilograms_float_feature_rt)`.
#[test]
fn xim_float_feature_rt() {
    let path = scratch("xim-float").join("float.xim");
    let mut feature_rt = Float32Builder::new();
    feature_rt.append_value(42.5);
    write(
        &path,
        vec![
            ("RUN_ID", DataType::Int64),
            ("PRECURSOR_ID", DataType::Int64),
            ("FEATURE_RT", DataType::Float32),
            ("MOBILITY_DATA", DataType::Binary),
            ("INTENSITY_DATA", DataType::Binary),
            ("MOBILITY_COMPRESSION", DataType::Int64),
            ("INTENSITY_COMPRESSION", DataType::Int64),
        ],
        vec![
            int64(12),
            int64(303),
            Arc::new(feature_rt.finish()),
            binary(&[1.1, 1.2]),
            binary(&[30.0, 40.0]),
            int64(0),
            int64(0),
        ],
    );
    let xim = XimParquetFile::new(&path).unwrap();
    let query = XimQuery {
        feature_rt: Some(42.5),
        ..XimQuery::default()
    };
    let mobilograms = xim.mobilograms(&query).unwrap();
    assert_eq!(mobilograms.len(), 1);
    assert_eq!(mobilograms[0].feature_rt, Some(42.5));
}

/// `START_SECTION(void getRuns(std::vector<XIMRunInfo>& output) const)`.
#[test]
fn xim_get_runs() {
    let xim = XimParquetFile::new(data("XIMParquetFile_23_input.xim")).unwrap();
    let runs = xim.runs().unwrap();
    assert!(!runs.is_empty());
    assert_ne!(runs[0].run_id, 0);
}

/// `START_SECTION(void load_invalid_path)`.
#[test]
fn xim_load_invalid_path() {
    assert!(is_not_found(XimParquetFile::new("no_such_file.xim")));
}

/// `START_SECTION(void getColumns(std::vector<std::string>&) const)`.
#[test]
fn xim_get_columns() {
    let xim = XimParquetFile::new(data("XIMParquetFile_23_input.xim")).unwrap();
    let columns = xim.columns().unwrap();
    assert!(!columns.is_empty());
    for wanted in ["RUN_ID", "MOBILITY_DATA", "INTENSITY_DATA"] {
        assert!(columns.iter().any(|c| c == wanted), "{wanted}");
    }
}

/// `START_SECTION(void getAnalytes(...) const)`.
#[test]
fn xim_get_analytes() {
    let xim = XimParquetFile::new(data("XIMParquetFile_23_input.xim")).unwrap();
    let none: [&str; 0] = [];
    assert!(!xim.analytes(&none, true).unwrap().is_empty());
    assert!(!xim.analytes(&none, false).unwrap().is_empty());
    let filtered = xim
        .analytes(
            &["PRECURSOR_ID", "MODIFIED_SEQUENCE", "PRECURSOR_CHARGE"],
            false,
        )
        .unwrap();
    assert!(!filtered.is_empty());
    assert!(is_invalid(xim.analytes(&["RUN_ID"], false)));
}

// --------------------------------------------------------------- XIPM ----

/// `START_SECTION(void load(std::vector<XIPMPeakMap>& output) const)`, on
/// the retained fixture (see the file documentation).
#[test]
fn xipm_load() {
    let xipm = XipmParquetFile::new(data("XIPMParquetFile_1_input.xipm")).unwrap();
    let peak_maps = xipm.load().unwrap();
    assert_eq!(peak_maps.len(), 2);
    let first = &peak_maps[0];
    assert_eq!(first.mz.len(), first.rt.len());
    assert_eq!(first.mz.len(), first.ion_mobility.len());
    assert_eq!(first.mz.len(), first.intensity.len());
    // The fixture has no TARGET_RT column: the value is the RT_START/RT_END
    // midpoint, (95 + 105) / 2, which is what the class test's own writer
    // stores as TARGET_RT.
    assert_eq!(first.target_rt, Some(100.0));
}

/// `START_SECTION(void getPeakMaps(...filters...) const)`.
#[test]
fn xipm_get_peak_maps() {
    let xipm = XipmParquetFile::new(data("XIPMParquetFile_1_input.xipm")).unwrap();
    let precursor = xipm
        .peak_maps(&XipmQuery {
            ms_level: Some(1),
            run_id: Some(7),
            peakmap_type: "precursor".into(),
            ..XipmQuery::default()
        })
        .unwrap();
    assert_eq!(precursor.len(), 1);
    assert_eq!(precursor[0].peakmap_type, "precursor");
    let transition = xipm
        .peak_maps(&XipmQuery {
            transition_id: Some(1),
            ms_level: Some(2),
            run_id: Some(7),
            peakmap_type: "transition".into(),
            ..XipmQuery::default()
        })
        .unwrap();
    assert_eq!(transition.len(), 1);
    assert_eq!(transition[0].annotation, "y7^1");
    assert_eq!(transition[0].target_rt, Some(100.0));
}

/// `START_SECTION(void getPeakMaps_multi_file)`.
#[test]
fn xipm_multi_file() {
    let file = data("XIPMParquetFile_1_input.xipm");
    let xipm = XipmParquetFile::from_files(&[&file, &file]).unwrap();
    assert_eq!(xipm.load().unwrap().len(), 4);
}

/// `START_SECTION(void getRuns(std::vector<XIPMRunInfo>& output) const)`.
#[test]
fn xipm_get_runs() {
    let xipm = XipmParquetFile::new(data("XIPMParquetFile_1_input.xipm")).unwrap();
    let runs = xipm.runs().unwrap();
    assert_eq!(runs.len(), 1);
    assert_eq!(runs[0].run_id, 7);
}

/// `START_SECTION(void getColumns(std::vector<std::string>& output) const)`,
/// without its `TARGET_RT` check: the retained fixture predates that column.
#[test]
fn xipm_get_columns() {
    let xipm = XipmParquetFile::new(data("XIPMParquetFile_1_input.xipm")).unwrap();
    let columns = xipm.columns().unwrap();
    for wanted in ["MZ_DATA", "RT_DATA", "MOBILITY_DATA"] {
        assert!(columns.iter().any(|c| c == wanted), "{wanted}");
    }
}

/// `START_SECTION(void load_invalid_path)`.
#[test]
fn xipm_load_invalid_path() {
    assert!(is_not_found(XipmParquetFile::new("no_such_file.xipm")));
    let empty: [&str; 0] = [];
    assert!(is_not_found(XipmParquetFile::from_files(&empty)));
}
