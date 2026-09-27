// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! `ParamJSONFile::load` on 115 documents against the C++ Release build: the
//! B8c lane's oracle (`../oracle/B8c-numbers`, mode `json`: the Linux x86_64
//! Release install at core `bc9cc12`, cli `c19e494`, topp `174b576` on
//! ibminode06, gcc 14.4.0, two identical runs; retained in
//! `tests/data/b8c_release/json`, provenance
//! `tests/data/b8c_release_provenance.json`).
//!
//! The documents (`make_inputs.py`) take every path of `nlohmann::json`
//! 3.11.2's scanner and parser that a load can meet: the byte order mark,
//! literals, strings with escapes, surrogates, control characters and each
//! class of UTF-8 lead byte, every state of the number automaton and a number
//! beyond `double`, objects, arrays and their separators, positions over
//! several lines, and the NUL byte the scanner takes for the end of input.
//! Each was loaded onto a tree with the parameters `T:1:d`, `T:1:i`, `T:1:s`,
//! `T:1:a:d` and `T:1:dl`; `release/<case>.kind` is the exception's class
//! (`ok` when the load succeeded) and `release/<case>.what` its `what()`, byte
//! for byte, or the values the load set.
//!
//! A refused document is compared on the whole text: this port's
//! `Error::Parse` message is the source's `what()` without the ` in: ` that
//! `ParseError` appends for its empty expression (`Exception.cpp:160-164`),
//! with the oracle's input directory mapped to this repository's. A byte that
//! is not UTF-8 in the source's `last read: '…'` is U+FFFD here, where the
//! message is text.
#![cfg(feature = "paramxml")]

use openms::Error;
use openms::format::param_json;
use openms::param::{Param, ParamValue};
use std::fs;
use std::path::{Path, PathBuf};

/// Where the oracle read its documents.
const ORACLE_INPUTS: &str = "/scratch/kohlbach/b8c-oracle/inputs/json/";

fn data(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/b8c_release/json")
        .join(relative)
}

/// The driver's tree (`driver.cpp`, `json`).
fn tree() -> Param {
    let mut param = Param::new();
    for (key, value) in [
        ("T:1:d", ParamValue::Float(0.0)),
        ("T:1:i", ParamValue::Integer(0)),
        ("T:1:s", ParamValue::String(String::new())),
        ("T:1:a:d", ParamValue::Float(0.0)),
        ("T:1:dl", ParamValue::FloatList(Vec::new())),
    ] {
        param.set_value(key, value, "", &[]).unwrap();
    }
    param
}

/// The values the driver printed after a successful load:
/// `d=<d> i=<i> s=<s> a:d=<d> dl=[<d>, ...]`, doubles at precision 17.
fn assert_values(case: &str, param: &Param, printed: &str) {
    let rest = printed.strip_prefix("d=").unwrap();
    let (d, rest) = rest.split_once(" i=").unwrap();
    let (i, rest) = rest.split_once(" s=").unwrap();
    let (s, rest) = rest.split_once(" a:d=").unwrap();
    let (nested, list) = rest.split_once(" dl=").unwrap();
    let float = |key: &str| param.value(key).unwrap().to_f64().unwrap();
    assert_eq!(float("T:1:d"), d.parse::<f64>().unwrap(), "{case}");
    assert_eq!(
        param.value("T:1:i").unwrap(),
        &ParamValue::Integer(i.parse().unwrap()),
        "{case}"
    );
    assert_eq!(
        param.value("T:1:s").unwrap(),
        &ParamValue::String(s.to_owned()),
        "{case}"
    );
    assert_eq!(float("T:1:a:d"), nested.parse::<f64>().unwrap(), "{case}");
    let items = list.strip_prefix('[').unwrap().strip_suffix(']').unwrap();
    let expected: Vec<f64> = if items.is_empty() {
        Vec::new()
    } else {
        items
            .split(", ")
            .map(|item| item.parse().unwrap())
            .collect()
    };
    assert_eq!(
        param.value("T:1:dl").unwrap(),
        &ParamValue::FloatList(expected),
        "{case}"
    );
}

#[test]
fn every_document_loads_or_fails_as_on_the_release_build() {
    let cases = fs::read_to_string(data("json_cases.tsv")).unwrap();
    let (mut refused, mut loaded) = (0, 0);
    for line in cases.lines() {
        let (case, file) = line.split_once('\t').unwrap();
        let kind = fs::read_to_string(data(&format!("release/{case}.kind"))).unwrap();
        let what = fs::read(data(&format!("release/{case}.what"))).unwrap();
        let input = data(&format!("inputs/{file}"));
        let mut param = tree();
        let result = param_json::load(&input, &mut param);
        match kind.trim_end() {
            "ok" => {
                assert!(result.is_ok(), "{case}: {result:?}");
                assert_values(case, &param, std::str::from_utf8(&what).unwrap());
                loaded += 1;
            }
            "Parse Error" => {
                let release = String::from_utf8_lossy(&what)
                    .replace(ORACLE_INPUTS, &format!("{}/", data("inputs").display()));
                let release = release
                    .strip_suffix(" in: ")
                    .unwrap_or_else(|| panic!("{case}: {release}"));
                match result {
                    Err(Error::Parse { line: 0, message }) => {
                        assert_eq!(message, release, "{case}");
                    }
                    other => panic!("{case}: {other:?}, the Release build: {release}"),
                }
                assert_eq!(param, tree(), "{case}: a failed load changes nothing");
                refused += 1;
            }
            other => panic!("{case}: {other}"),
        }
    }
    assert_eq!((refused, loaded), (106, 9));
}
