// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Class tests for the three feature-text adapters `KroenikFile`,
//! `MsInspectFile` and `SpecArrayFile`.
//!
//! Every expected value is transcribed from the upstream class tests at the
//! core pin (tier 2): `KroenikFile_test.cpp`, `MsInspectFile_test.cpp` and
//! `SpecArrayFile_test.cpp`, run against their unmodified fixtures (see
//! `tests/data/search_engine_files_provenance.json`). `TEST_EQUAL` becomes an
//! exact comparison and `TEST_REAL_SIMILAR` the ClassTest tolerance. The
//! error-path cases that the upstream tests do not cover are derived from the
//! source's control flow and say so.

use openms::format::kroenik_file::KroenikFile;
use openms::format::ms_inspect_file::MsInspectFile;
use openms::format::spec_array_file::SpecArrayFile;
use openms::kernel::features::FeatureMap;
use openms::metadata::{MetaValue, MetaValueData};
use openms::{Error, MSSpectrum};
use std::io::Cursor;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// ClassTest `isRealSimilar` (`ClassTest.cpp:364`, testframework): absolute
/// difference at most 1e-5, or ratio at most 1 + 1e-5, with the source's zero
/// cases.
fn real_similar(actual: f64, expected: f64) -> bool {
    if actual.is_nan() || expected.is_nan() {
        return false;
    }
    let small = (actual - expected).abs() <= 1e-5;
    if actual == 0.0 {
        return expected == 0.0 || small;
    }
    if expected == 0.0 {
        return small;
    }
    let mut ratio = actual / expected;
    if ratio < 0.0 {
        return small;
    }
    if ratio < 1.0 {
        ratio = 1.0 / ratio;
    }
    ratio <= 1.0 + 1e-5 || small
}

fn assert_similar(actual: f64, expected: f64) {
    assert!(
        real_similar(actual, expected),
        "{actual} is not similar to {expected}"
    );
}

fn meta<'a>(map: &'a FeatureMap, index: usize, name: &str) -> &'a MetaValue {
    map.features[index]
        .metadata
        .get(name)
        .unwrap_or_else(|| panic!("feature {index} lacks meta value {name}"))
}

/// A meta value that must be stored as a double, as the source's
/// `setMetaValue(name, toDouble(...))` stores it.
fn float_meta(map: &FeatureMap, index: usize, name: &str) -> f64 {
    match meta(map, index, name).data() {
        MetaValueData::Float(value) => *value,
        other => panic!("meta value {name} of feature {index} is not a double: {other:?}"),
    }
}

/// A meta value that must be stored as an integer.
fn int_meta(map: &FeatureMap, index: usize, name: &str) -> i64 {
    match meta(map, index, name).data() {
        MetaValueData::Integer(value) => *value,
        other => panic!("meta value {name} of feature {index} is not an integer: {other:?}"),
    }
}

/// A meta value that must be stored as a string.
fn string_meta<'a>(map: &'a FeatureMap, index: usize, name: &str) -> &'a str {
    match meta(map, index, name).data() {
        MetaValueData::String(value) => value,
        other => panic!("meta value {name} of feature {index} is not a string: {other:?}"),
    }
}

#[test]
fn kroenik_load_matches_the_class_test() {
    // KroenikFile_test.cpp, section load.
    let map = KroenikFile::new()
        .load(data("KroenikFile_test_1.krf"))
        .unwrap();
    assert_eq!(map.features.len(), 3);
    let first = &map.features[0];
    assert_eq!(first.rt, 63.2);
    assert_similar(first.mz, 1002.11);
    assert_eq!(first.intensity, 999_999.0);
    assert_eq!(first.charge, 1);
    assert_eq!(string_meta(&map, 0, "AveragineModifications"), "Carbamido");
    let second = &map.features[1];
    assert_eq!(second.rt, 62.2);
    assert_similar(second.mz, 252.057);
    assert_eq!(second.intensity, 9999.0);
    assert_eq!(second.charge, 2);
    assert_eq!(string_meta(&map, 1, "AveragineModifications"), "Carbamido2");
}

#[test]
fn kroenik_load_refuses_a_line_with_fifteen_columns() {
    // KroenikFile_test.cpp: TEST_EXCEPTION(ParseError) on KroenikFile_test_2.krf,
    // whose third line (the second data line) carries a fifteenth column.
    match KroenikFile::new().load(data("KroenikFile_test_2.krf")) {
        Err(Error::Parse { line, message }) => {
            assert_eq!(line, 3);
            assert!(
                message.starts_with(
                    "Failed parsing in line 3: missing 14 tab-separated entries (got 15)"
                ),
                "{message}"
            );
        }
        other => panic!("expected a parse error, got {other:?}"),
    }
}

#[test]
fn kroenik_load_of_a_missing_file_is_an_io_error() {
    // KroenikFile_test.cpp: TEST_EXCEPTION(FileNotFound).
    let error = KroenikFile::new()
        .load(data("KroenikFile_test_2_doesnotexist.edta"))
        .unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error:?}");
}

#[test]
fn kroenik_store_is_not_implemented() {
    // KroenikFile_test.cpp: TEST_EXCEPTION(NotImplemented, f.store("bla", spec)).
    let error = KroenikFile::new()
        .store("bla", &MSSpectrum::default())
        .unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error:?}");
}

#[test]
fn kroenik_feature_carries_the_source_hull_and_meta_values() {
    // Derived from KroenikFile.cpp: the hull spans First..Last RTime by
    // m/z..m/z + 3/charge, and the meta values keep the source's types
    // (FirstScan and NumOfScans double, LastScan integer).
    let map = KroenikFile::new()
        .load(data("KroenikFile_test_1.krf"))
        .unwrap();
    let second = &map.features[1];
    assert_eq!(second.quality, 0.8_f32);
    assert_eq!(second.convex_hulls.len(), 1);
    let bounds = second.convex_hulls[0].bounding_box().unwrap();
    assert_eq!(bounds.min().rt, 61.1);
    assert_eq!(bounds.max().rt, 64.2);
    assert_eq!(bounds.min().mz, second.mz);
    assert_eq!(bounds.max().mz, second.mz + 3.0 / 2.0);
    assert_eq!(float_meta(&map, 1, "Mass"), 502.1);
    assert_eq!(float_meta(&map, 1, "FirstScan"), 8.0);
    assert_eq!(int_meta(&map, 1, "LastScan"), 9);
    assert_eq!(float_meta(&map, 1, "NumOfScans"), 2.0);
}

#[test]
fn kroenik_empty_input_yields_an_empty_map() {
    // Derived from KroenikFile.cpp: `if (it == input.end()) return;`.
    let map = KroenikFile::new().load_reader(Cursor::new("")).unwrap();
    assert!(map.features.is_empty());
}

#[test]
fn kroenik_refuses_a_zero_charge() {
    // Native refusal: the source divides the mass by zero and stores an
    // infinite m/z.
    let text = "header\nf\t3\t5\t3\t0\t1001.1\t1000.0\t1\t2\t60.1\t66.2\t63.2\t0.98\tmod\n";
    assert!(matches!(
        KroenikFile::new().load_reader(Cursor::new(text)),
        Err(Error::Parse { line: 2, .. })
    ));
}

#[test]
fn ms_inspect_load_matches_the_class_test() {
    // MsInspectFile_test.cpp, section load.
    let map = MsInspectFile::new()
        .load(data("MSInspectFile_test_1.msi"))
        .unwrap();
    assert_eq!(map.features.len(), 2);
    let first = &map.features[0];
    assert_similar(first.rt, 12.92);
    assert_similar(first.mz, 501.51);
    assert_similar(f64::from(first.intensity), 45677.0);
    assert_similar(f64::from(first.quality), 0.98);
    assert_eq!(float_meta(&map, 0, "background"), 0.11);
    let second = &map.features[1];
    assert_similar(second.rt, 22.92);
    assert_similar(second.mz, 601.51);
    assert_similar(f64::from(second.intensity), 245_677.0);
    assert_similar(f64::from(second.quality), 0.99);
    assert_eq!(float_meta(&map, 1, "background"), 0.22);
}

#[test]
fn ms_inspect_meta_values_keep_the_source_types() {
    // Derived from MsInspectFile.h: accurateMZ and description are strings,
    // chargeStates and the scan columns integers, mass a double.
    let map = MsInspectFile::new()
        .load(data("MSInspectFile_test_1.msi"))
        .unwrap();
    assert_eq!(string_meta(&map, 0, "accurateMZ"), "501.2");
    assert_eq!(float_meta(&map, 0, "mass"), 1001.0);
    assert_eq!(int_meta(&map, 1, "chargeStates"), 6);
    assert_eq!(int_meta(&map, 1, "scanCount"), 5);
    assert_eq!(float_meta(&map, 1, "totalIntensity"), 111_111_111.0);
    assert_eq!(string_meta(&map, 1, "description"), "example_feature_2");
    assert_eq!(map.features[1].charge, 5);
}

#[test]
fn ms_inspect_reports_the_failing_column() {
    // Derived from MsInspectFile.h: the ParseError names column_to_convert + 1.
    let text = "scan\ttime\n3\t12.92\t501.51\t501.2\t1001\t45677\ttwo\t2\t0.98\t0.11\t211\t231\t1\t5\t5\t3\t2\td\n";
    match MsInspectFile::new().load_reader(Cursor::new(text)) {
        Err(Error::Parse { line, message }) => {
            assert_eq!(line, 2);
            assert_eq!(
                message,
                "Failed to convert value in column 7 into a number (line '2)"
            );
        }
        other => panic!("expected a parse error, got {other:?}"),
    }
    let short = "# c\nheader\n1\t2\t3\n";
    match MsInspectFile::new().load_reader(Cursor::new(short)) {
        Err(Error::Parse { line, message }) => {
            assert_eq!(line, 3);
            assert_eq!(
                message,
                "Failed to convert line 3. Not enough columns (expected 18 or more, got 3)"
            );
        }
        other => panic!("expected a parse error, got {other:?}"),
    }
}

#[test]
fn ms_inspect_store_is_not_implemented() {
    // MsInspectFile_test.cpp: TEST_EXCEPTION(NotImplemented, f.store("bla", spec)).
    let error = MsInspectFile::new()
        .store("bla", &MSSpectrum::default())
        .unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error:?}");
}

#[test]
fn spec_array_load_matches_the_class_test() {
    // SpecArrayFile_test.cpp, section load.
    let map = SpecArrayFile::new()
        .load(data("SpecArrayFile_test_1.peplist"))
        .unwrap();
    assert_eq!(map.features.len(), 2);
    let first = &map.features[0];
    assert_eq!(first.rt, 60.1 * 60.0);
    assert_similar(first.mz, 500.1);
    assert_eq!(first.intensity, 4_343_534.0);
    assert_eq!(first.charge, 5);
    assert_eq!(float_meta(&map, 0, "s/n"), 3.2);
    let second = &map.features[1];
    assert_eq!(second.rt, 40.1 * 60.0);
    assert_similar(second.mz, 700.1);
    assert_eq!(second.intensity, 222_432.0);
    assert_eq!(second.charge, 3);
    assert_eq!(float_meta(&map, 1, "s/n"), 2.2);
}

#[test]
fn spec_array_load_refuses_a_short_line() {
    // SpecArrayFile_test.cpp: TEST_EXCEPTION(ParseError) on
    // SpecArrayFile_test_2.peplist, whose second data line has four columns.
    match SpecArrayFile::new().load(data("SpecArrayFile_test_2.peplist")) {
        Err(Error::Parse { line, message }) => {
            assert_eq!(line, 3);
            assert_eq!(
                message,
                "Failed to convert line3not enough columns (expected 5 or more, got 4)"
            );
        }
        other => panic!("expected a parse error, got {other:?}"),
    }
}

#[test]
fn spec_array_load_of_a_missing_file_is_an_io_error() {
    // SpecArrayFile_test.cpp: TEST_EXCEPTION(FileNotFound).
    let error = SpecArrayFile::new()
        .load(data("SpecArrayFile_test_2_doesnotexist.peplist"))
        .unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error:?}");
}

#[test]
fn spec_array_store_is_not_implemented() {
    // SpecArrayFile_test.cpp: TEST_EXCEPTION(NotImplemented, f.store("bla", spec)).
    let error = SpecArrayFile::new()
        .store("bla", &MSSpectrum::default())
        .unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error:?}");
}

#[test]
fn spec_array_refuses_a_non_numeric_value() {
    // Derived from SpecArrayFile.h: any conversion failure becomes one
    // ParseError naming the line.
    let text = "m/z\trt(min)\tsnr\tcharge\tintensity\n500.1\t60.1\t3.2\t5.0\t1\n";
    match SpecArrayFile::new().load_reader(Cursor::new(text)) {
        Err(Error::Parse { line, message }) => {
            assert_eq!(line, 2);
            assert_eq!(message, "Failed to convert value into a number (line '2)");
        }
        other => panic!("expected a parse error, got {other:?}"),
    }
}
