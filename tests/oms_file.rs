// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! OMS files (`FORMAT/OMSFile.h`).
//!
//! Tier 1: the retained JSONExporter outputs of the upstream TOPP suite
//! (test-data 0cb15f2 `topp/CMakeLists.txt:3664-3672`) are reproduced from the
//! retained `.oms` inputs through `OmsFile::export_to_json`, which is all the
//! tool does (topp 174b576 `src/JSONExporter.cpp`), and compared with the
//! upstream `${DIFF}` contract (FuzzyDiff, ratio 1.01, absdiff 0.01, the INI's
//! whitelist).

#![cfg(feature = "sqlite")]

#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;

use openms::concept::progress_logger::ProgressLogType;
use openms::format::oms_file::{OmsFile, SCHEMA_VERSION};
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/oms_file")
        .join(name)
}

fn assert_fuzzy_equal(produced: &str, retained: &str) {
    let settings = fuzzy::FuzzyDiffSettings::upstream().unwrap();
    assert_eq!(settings.ratio, 1.01);
    assert_eq!(settings.absdiff, 0.01);
    let expected = std::fs::read(data(retained)).unwrap();
    if let Err(log) = settings.compare_bytes(produced.as_bytes(), &expected) {
        panic!("{retained}: FuzzyDiff failed:\n{log}");
    }
}

#[test]
fn retained_json_exporter_1_protein() {
    // TOPP_JSONExporter_1: -in JSONExporter_protein.oms
    let text = OmsFile::default()
        .json_text(data("JSONExporter_protein.oms"))
        .unwrap();
    assert_fuzzy_equal(&text, "JSONExporter_1_output.json");
}

#[test]
fn retained_json_exporter_2_rna() {
    // TOPP_JSONExporter_2: -in JSONExporter_RNA.oms
    let text = OmsFile::default()
        .json_text(data("JSONExporter_RNA.oms"))
        .unwrap();
    assert_fuzzy_equal(&text, "JSONExporter_2_output.json");
}

#[test]
fn export_to_json_writes_the_rendered_text() {
    let dir = std::env::temp_dir().join(format!("openms-oms-json-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let out = dir.join("protein.json");
    let adapter = OmsFile::new(ProgressLogType::None);
    adapter
        .export_to_json(data("JSONExporter_protein.oms"), &out)
        .unwrap();
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        adapter.json_text(data("JSONExporter_protein.oms")).unwrap()
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn retained_inputs_declare_the_pinned_schema_version() {
    let adapter = OmsFile::default();
    assert_eq!(
        adapter
            .schema_version(data("JSONExporter_protein.oms"))
            .unwrap(),
        SCHEMA_VERSION
    );
    assert_eq!(
        adapter.schema_version(data("JSONExporter_RNA.oms")).unwrap(),
        SCHEMA_VERSION
    );
}

#[test]
fn a_file_that_is_not_an_oms_database_is_refused() {
    let adapter = OmsFile::default();
    // Not an SQLite database at all.
    assert!(adapter
        .json_text(data("JSONExporter_1_output.json"))
        .is_err());
    // A missing file is not created by the read-only open.
    let missing = data("does_not_exist.oms");
    assert!(adapter.json_text(&missing).is_err());
    assert!(!missing.exists());
}
