// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `StringUtils::toInt32` / `toFloat` / `toDouble` and `ProteinIdentification::sort`
//! against the Linux x86_64 Release build (tier 1).
//!
//! `tests/data/string_utils_release/cpp_release_records.tsv` is the output of
//! `../oracle/integ-F-stringutils/driver.cpp`, built against
//! `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576` with the
//! conda-forge g++ 14.4.0 and libOpenMS's core flags on ibminode06, and run twice
//! with identical output. Each record is `kind <TAB> input <TAB> ... <TAB> result`,
//! where the input escapes tab, CR, LF and bytes outside printable ASCII as `\t`,
//! `\r`, `\n` and `\xNN`, and a result starting `ERR:` is the text of the
//! `ConversionError` the source threw. Provenance:
//! `tests/data/string_utils_release_provenance.json`.

use openms::data_structures::string_utils::{to_double, to_float, to_int32};
use openms::identification::{ProteinHit, ProteinIdentification};

const RECORDS: &str = include_str!("data/string_utils_release/cpp_release_records.tsv");

fn unescape(text: &str) -> String {
    let mut out = Vec::new();
    let bytes = text.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'\\' && i + 1 < bytes.len() {
            match bytes[i + 1] {
                b't' => {
                    out.push(b'\t');
                    i += 2;
                    continue;
                }
                b'n' => {
                    out.push(b'\n');
                    i += 2;
                    continue;
                }
                b'r' => {
                    out.push(b'\r');
                    i += 2;
                    continue;
                }
                b'x' => {
                    let hex = std::str::from_utf8(&bytes[i + 2..i + 4]).unwrap();
                    out.push(u8::from_str_radix(hex, 16).unwrap());
                    i += 4;
                    continue;
                }
                _ => {}
            }
        }
        out.push(bytes[i]);
        i += 1;
    }
    String::from_utf8(out).unwrap()
}

/// The source's `%.17g` / `%.9g` rendering of the driver, for the values the
/// records hold: `nan` / `-nan`, `inf` / `-inf`, or a number compared by value.
fn same_double(expected: &str, actual: f64) -> bool {
    match expected {
        "nan" => actual.is_nan() && actual.is_sign_positive(),
        "-nan" => actual.is_nan() && actual.is_sign_negative(),
        "inf" => actual == f64::INFINITY,
        "-inf" => actual == f64::NEG_INFINITY,
        "-0" => actual == 0.0 && actual.is_sign_negative(),
        "0" => actual == 0.0 && actual.is_sign_positive(),
        _ => expected.parse::<f64>().unwrap() == actual,
    }
}

fn same_float(expected: &str, actual: f32) -> bool {
    match expected {
        "nan" => actual.is_nan() && actual.is_sign_positive(),
        "-nan" => actual.is_nan() && actual.is_sign_negative(),
        "inf" => actual == f32::INFINITY,
        "-inf" => actual == f32::NEG_INFINITY,
        "-0" => actual == 0.0 && actual.is_sign_negative(),
        "0" => actual == 0.0 && actual.is_sign_positive(),
        // `%.9g` round-trips binary32.
        _ => expected.parse::<f32>().unwrap() == actual,
    }
}

fn records(kind: &str) -> Vec<(String, String)> {
    RECORDS
        .lines()
        .filter_map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            (fields[0] == kind).then(|| (unescape(fields[1]), (*fields.last().unwrap()).to_owned()))
        })
        .collect()
}

#[test]
fn to_int32_matches_the_release_build() {
    let cases = records("toInt32");
    assert_eq!(cases.len(), 82);
    for (input, expected) in cases {
        match (expected.strip_prefix("ERR:"), to_int32(&input)) {
            (Some(message), Err(error)) => assert_eq!(error, message, "{input:?}"),
            (None, Ok(value)) => assert_eq!(value.to_string(), expected, "{input:?}"),
            (want, got) => panic!("{input:?}: source {want:?} / {expected}, port {got:?}"),
        }
    }
}

#[test]
fn to_double_matches_the_release_build() {
    let cases = records("toDouble");
    assert_eq!(cases.len(), 82);
    for (input, expected) in cases {
        match (expected.strip_prefix("ERR:"), to_double(&input)) {
            (Some(message), Err(error)) => assert_eq!(error, message, "{input:?}"),
            (None, Ok(value)) => assert!(same_double(&expected, value), "{input:?}: {value:?}"),
            (want, got) => panic!("{input:?}: source {want:?} / {expected}, port {got:?}"),
        }
    }
}

#[test]
fn to_float_matches_the_release_build() {
    let cases = records("toFloat");
    assert_eq!(cases.len(), 82);
    for (input, expected) in cases {
        match (expected.strip_prefix("ERR:"), to_float(&input)) {
            (Some(message), Err(error)) => assert_eq!(error, message, "{input:?}"),
            (None, Ok(value)) => assert!(same_float(&expected, value), "{input:?}: {value:?}"),
            (want, got) => panic!("{input:?}: source {want:?} / {expected}, port {got:?}"),
        }
    }
}

#[test]
fn protein_identification_sort_matches_the_release_build() {
    let sorted = records("ProteinIdentification::sort");
    assert_eq!(sorted.len(), 2);
    for (direction, expected) in sorted {
        let higher = match direction.as_str() {
            "higher_score_better=true" => true,
            "higher_score_better=false" => false,
            other => panic!("unexpected record {other}"),
        };
        // The driver's insertion order.
        let hits = [
            (1.0, "P1"),
            (1.0, "P2"),
            (2.0, "P3"),
            (1.0, "P10"),
            (2.0, "A"),
            (0.5, "Z"),
            (1.0, "p1"),
        ];
        let mut identification = ProteinIdentification {
            higher_score_better: higher,
            hits: hits
                .iter()
                .map(|&(score, accession)| ProteinHit::new(score, 0, accession, "").unwrap())
                .collect(),
            ..ProteinIdentification::default()
        };
        identification.sort().unwrap();
        let got = identification
            .hits
            .iter()
            .map(|hit| format!("{}:{}", hit.accession, hit.score))
            .collect::<Vec<_>>()
            .join(",");
        assert_eq!(got, expected, "{direction}");
    }
}
