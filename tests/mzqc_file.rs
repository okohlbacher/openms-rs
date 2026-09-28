// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! mzQC files (`FORMAT/MzQCFile.h`) and `FileHandler::computeFileHash`.
//!
//! Tier 1: TOPP_QCCalculator_2 (test-data 0cb15f2
//! `topp/CMakeLists.txt:5387-5389`) is reproduced from its retained inputs
//! through the library, as the tool (topp 174b576 `src/QCCalculator.cpp`)
//! drives it - load the mzML, sort its spectra, load and RT-sort the
//! features, load the identifications, `FileHandler::storeQC` into
//! `MzQCFile::store` - and compared with the retained
//! `QCCalculator_2_output.mzQC` under the upstream `${DIFF}` contract
//! (`FuzzyDiff -ini FuzzyDiff.ini`: ratio 1.01, absdiff 0.01) with the test's
//! whitelist `creationDate`, `location`, `version`. There is no MzQCFile class
//! test upstream. Tier 4: the documented behaviour without inputs, and the
//! checksum of an unreadable file.

#![cfg(all(feature = "mzml", feature = "featurexml", feature = "idxml"))]

#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;

use openms::data_structures::DateTime;
use openms::format::FileHandler;
use openms::format::FileType;
use openms::format::idxml::IdXmlDocument;
use openms::format::mzqc_file::{MzQcFile, MzQcInput};
use openms::kernel::{FeatureMap, MSExperiment};
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/mzqc_file")
        .join(name)
}

fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms-mzqc-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The three inputs as QCCalculator loads them (topp 174b576
/// `src/QCCalculator.cpp:139-165`).
fn tool_inputs() -> (MSExperiment, FeatureMap, IdXmlDocument) {
    let mut exp =
        FileHandler::load_experiment(data("QCCalculator_input.mzML"), &[FileType::MzMl]).unwrap();
    exp.sort_spectra(true).unwrap();
    let mut features = FileHandler::load_feature_map(
        data("FeatureFinderMetaboIdent_1_output.featureXML"),
        &[FileType::FeatureXml],
    )
    .unwrap();
    features.sort_by_rt().unwrap();
    let ids = FileHandler::load_identifications(data("OpenPepXL_output.idXML"), &[FileType::IdXml])
        .unwrap();
    (exp, features, ids)
}

#[test]
fn qccalculator_2_reproduces_the_retained_mzqc() {
    let (exp, features, mut ids) = tool_inputs();
    let input_file = data("QCCalculator_input.mzML");
    let input_file = input_file.to_str().unwrap();
    let input = MzQcInput {
        input_file,
        experiment: &exp,
        contact_name: "name",
        contact_address: "address",
        description: "description",
        label: "label",
        feature_map: &features,
    };
    let out = temp_dir("qccalculator-2").join("QCCalculator_2.mzQC");
    let report = MzQcFile::new()
        .store(
            &out,
            &input,
            &ids.protein_identifications,
            &mut ids.peptide_identifications,
        )
        .unwrap();
    // Every metric's inputs are present, so nothing is withheld.
    assert!(
        report.log.iter().all(|line| !line.contains("cannot run")),
        "{:?}",
        report.log
    );
    let produced = std::fs::read(&out).unwrap();
    let expected = std::fs::read(data("QCCalculator_2_output.mzQC")).unwrap();
    let settings = fuzzy::FuzzyDiffSettings::upstream()
        .unwrap()
        .with_whitelist(&["creationDate", "location", "version"]);
    assert_eq!(settings.ratio, 1.01);
    assert_eq!(settings.absdiff, 0.01);
    if let Err(log) = settings.compare_bytes(&produced, &expected) {
        panic!(
            "FuzzyDiff failed:\n{log}\n--- produced ---\n{}",
            String::from_utf8_lossy(&produced)
        );
    }
}

/// The lines FuzzyDiff skips are checked here for what they must hold: the
/// given creation date, the absolute input location, mzQC version 1.0.0 and
/// the core SDK version; and the input checksum is the SHA-1 the retained
/// output records (and `shasum` computes) for the retained input.
#[test]
fn whitelisted_lines_hold_what_the_source_writes() {
    let (exp, features, mut ids) = tool_inputs();
    let input_file = data("QCCalculator_input.mzML");
    let input = MzQcInput {
        input_file: input_file.to_str().unwrap(),
        experiment: &exp,
        contact_name: "",
        contact_address: "",
        description: "",
        label: "label",
        feature_map: &features,
    };
    let date = DateTime::parse("2025-05-16T14:22:42").unwrap();
    let (text, _) = MzQcFile::new()
        .text(
            &input,
            &ids.protein_identifications,
            &mut ids.peptide_identifications,
            &date,
        )
        .unwrap();
    assert!(text.starts_with(
        "{\n  \"mzQC\": {\n    \"creationDate\": \"2025-05-16T14:22:42\",\n    \"version\": \"1.0.0\",\n    \"runQualities\": ["
    ));
    let location = format!("\"location\": \"{}\"", input_file.to_str().unwrap());
    assert!(text.contains(&location), "{text}");
    assert!(text.contains("\"value\": \"24aeada3affec9df36568c0054e4197ac3515642\""));
    assert!(text.contains(&format!("\"version\": \"{}\"", openms::CORE_SDK_VERSION)));
    assert!(!text.ends_with('\n'));
    assert_eq!(
        FileHandler::compute_file_hash(&input_file),
        "24aeada3affec9df36568c0054e4197ac3515642"
    );
}

/// Without an input file, features or identifications, only the metrics that
/// need none of them are written, and each withheld metric logs the source's
/// note once.
#[test]
fn metrics_without_their_inputs_are_withheld_with_a_note() {
    let (exp, _, _) = tool_inputs();
    let features = FeatureMap::default();
    let input = MzQcInput {
        input_file: "",
        experiment: &exp,
        contact_name: "",
        contact_address: "",
        description: "",
        label: "",
        feature_map: &features,
    };
    let (text, report) = MzQcFile::new()
        .text(
            &input,
            &[],
            &mut [],
            &DateTime::parse("2025-01-01T00:00:00").unwrap(),
        )
        .unwrap();
    for present in ["QC:4000135", "QC:4000053", "QC:4000138"] {
        assert!(text.contains(present), "{present}");
    }
    for absent in ["QC:4000059", "QC:4000067", "QC:4000257", "QC:4000186"] {
        assert!(!text.contains(absent), "{absent}");
    }
    assert_eq!(
        report.log,
        [
            "Note: Metric 'SpectrumCount' cannot run because input data 'raw.mzML' is missing!",
            "Note: Metric 'TIC' cannot run because input data 'raw.mzML' is missing!",
            "Note: Metric 'Summary of features from featureXML file' cannot run because input data 'preFDR.featureXML' is missing!",
            "Note: Metric 'Summary of detected Proteins and Peptides from idXML file' cannot run because input data 'id.idXML' is missing!",
        ]
    );
    // An empty input file name: an empty checksum, the source's unreadable-file answer.
    assert!(text.contains("\"name\": \"SHA-1\""));
    assert!(text.contains("\"value\": \"\""));
}

#[test]
fn an_unreadable_file_has_an_empty_checksum() {
    assert_eq!(
        FileHandler::compute_file_hash(data("no-such-file.mzML")),
        ""
    );
}

/// Native: the document is complete before anything is written, so a failed
/// computation leaves an existing file untouched (the source has already
/// truncated it).
#[test]
fn a_failed_store_leaves_the_output_unchanged() {
    let (exp, features, mut ids) = tool_inputs();
    for prot in &mut ids.protein_identifications {
        prot.search_parameters.digestion_enzyme = "unknown_enzyme".into();
    }
    let out = temp_dir("failed-store").join("kept.mzQC");
    std::fs::write(&out, "previous").unwrap();
    let input = MzQcInput {
        input_file: "",
        experiment: &exp,
        contact_name: "",
        contact_address: "",
        description: "",
        label: "",
        feature_map: &features,
    };
    let result = MzQcFile::new().store(
        &out,
        &input,
        &ids.protein_identifications,
        &mut ids.peptide_identifications,
    );
    assert!(result.is_err());
    assert_eq!(std::fs::read_to_string(&out).unwrap(), "previous");
}
