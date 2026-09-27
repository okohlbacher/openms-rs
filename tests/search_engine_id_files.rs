// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Class tests for the search-engine identification adapters.
//!
//! Every expected value is transcribed from the upstream class test named at
//! each case, at the core pin (tier 2), and run against the unmodified
//! upstream fixtures described in `tests/data/search_engine_files_provenance.json`.
//! Cases the upstream tests do not cover are derived from the source's control
//! flow and say so.

use openms::Error;
use openms::chemistry::AASequence;
use openms::format::omssa_csv_file::OMSSACSVFile;
use std::io::Cursor;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

#[test]
fn omssa_csv_load_matches_the_class_test() {
    // OMSSACSVFile_test.cpp, section load: one identification (the protein
    // identification the source leaves untouched is not part of the Rust API).
    let ids = OMSSACSVFile::new()
        .load(data("OMSSACSVFile_test_1.csv"))
        .unwrap();
    assert_eq!(ids.len(), 1);
    // Derived from OMSSACSVFile.cpp: the three lines of spectrum 0 are three
    // hits with score type OMSSA, the P-value as score and charge 1.
    let id = &ids[0];
    assert_eq!(id.score_type, "OMSSA");
    assert_eq!(id.hits.len(), 3);
    for hit in &id.hits {
        assert_eq!(hit.sequence, AASequence::parse("MSHHWGYGK").unwrap());
        assert_eq!(hit.score, 1.3086e-08);
        assert_eq!(hit.charge, 1);
    }
}

#[test]
fn omssa_csv_refuses_an_unbalanced_quote() {
    // OMSSACSVFile_test.cpp: TEST_EXCEPTION(ParseError) on the unbalanced-quote
    // regression fixture (15 fields, 5 quoted commas).
    match OMSSACSVFile::new().load(data("OMSSACSVFile_test_unbalanced_quote.csv")) {
        Err(Error::Parse { line, message }) => {
            assert_eq!(line, 2);
            assert_eq!(message, "number of columns should be 14 in line 1");
        }
        other => panic!("expected a parse error, got {other:?}"),
    }
}

#[test]
fn omssa_csv_reads_shifted_columns_after_a_quoted_defline() {
    // Derived from OMSSACSVFile.cpp: two quoted commas shift the charge and
    // P-value two fields right; a new spectrum number starts a new
    // identification.
    let text = "Spectrum number, ...\n\
                0,,PEPTIDE,1,2,0,3,1,7,plain,,2,3,0.5\n\
                1,,PEPTIDER,1,2,0,3,1,8,\"a, b, c\",,3,4,0.25\n";
    let ids = OMSSACSVFile::new().load_reader(Cursor::new(text)).unwrap();
    assert_eq!(ids.len(), 2);
    assert_eq!(ids[0].hits[0].charge, 2);
    assert_eq!(ids[0].hits[0].score, 0.5);
    assert_eq!(ids[1].hits[0].charge, 3);
    assert_eq!(ids[1].hits[0].score, 0.25);
}

#[test]
fn omssa_csv_requires_the_spectrum_header() {
    // Derived from OMSSACSVFile.cpp: hasPrefix(line, "Spectrum") on the first
    // line, an empty file included.
    for text in ["", "number, ...\n"] {
        match OMSSACSVFile::new().load_reader(Cursor::new(text)) {
            Err(Error::Parse { line, message }) => {
                assert_eq!(line, 1);
                assert_eq!(message, "first line does not contain a description");
            }
            other => panic!("expected a parse error, got {other:?}"),
        }
    }
}

#[test]
fn omssa_csv_missing_file_is_an_io_error() {
    // Derived from OMSSACSVFile.cpp: FileNotFound when the stream fails.
    let error = OMSSACSVFile::new()
        .load(data("OMSSACSVFile_test_doesnotexist.csv"))
        .unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error:?}");
}
