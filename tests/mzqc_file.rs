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

use openms::chemistry::{ModificationsDB, ResidueModification};
use openms::data_structures::DateTime;
use openms::format::FileHandler;
use openms::format::FileType;
use openms::format::idxml;
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
///
/// The experiment is read from `QCCalculator_input_strict_reader.mzML`, the
/// retained input with the eight `dataProcessingRef` attributes of its m/z,
/// intensity and time arrays deleted and nothing else changed: the strict
/// mzML reader refuses processing on a primary array ("primary-array
/// processing has no independent native owner", `docs/FILE_INFO_SUPPORT.md`),
/// while the source's `MzMLHandler` (core bc9cc12, lines 433-434 and 733-739)
/// copies only the meta values of those arrays and drops their processing, so
/// it builds the same experiment from either file. The mzQC input-file
/// record (location, name, SHA-1) is written for the retained file name, as
/// the tool writes it.
fn tool_inputs() -> (MSExperiment, FeatureMap, IdXmlDocument) {
    let mut exp = FileHandler::load_experiment(
        data("QCCalculator_input_strict_reader.mzML"),
        &[FileType::MzMl],
    )
    .unwrap();
    exp.sort_spectra(true).unwrap();
    let mut features = FileHandler::load_feature_map(
        data("FeatureFinderMetaboIdent_1_output.featureXML"),
        &[FileType::FeatureXml],
    )
    .unwrap();
    features.sort_by_rt().unwrap();
    let ids = read_ids("OpenPepXL_output.idXML");
    (exp, features, ids)
}

/// The idXML reader charges every named sequence annotation as if it scanned
/// the whole registry (`AASequence::parse_with_budget`), so with the global
/// registry its fixed 50 M work budget refuses `OpenPepXL_output.idXML` (19
/// annotations) with "peptide parent parsing resource limit exceeded" - the
/// reader defect `tests/mztab_export.rs` and `docs/MZTAB_SUPPORT.md` record.
/// Reading through the global registry's own records for the names the file
/// uses resolves every modification to the same record.
fn read_ids(input: &str) -> IdXmlDocument {
    let records: Vec<ResidueModification> = ModificationsDB::global()
        .entries()
        .iter()
        .filter(|record| {
            [
                "Oxidation",
                "Carbamidomethyl",
                "Xlink:DSS[155]",
                "Xlink:DSS[156]",
            ]
            .contains(&record.name())
        })
        .map(|record| ResidueModification::clone(record))
        .collect();
    let registry = ModificationsDB::from_records(records).unwrap();
    let file = std::io::BufReader::new(std::fs::File::open(data(input)).unwrap());
    idxml::read_with_registry(file, &idxml::ReadOptions::default(), &registry).unwrap()
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

/// Exact agreement with the Linux x86_64 Release build (oracle `F2-mzqc`,
/// ibminode06, run twice with identical results; `release/host.txt`), which
/// ran QCCalculator's body on the same inputs: every line is identical but the
/// creation date and the input's absolute location, which depend on the run.
/// This pins the one number the retained output differs in: the RT mean shift
/// 0.916678249835968 of the pinned build (the retained file, written by an
/// older build, has 0.9166780114173889, within FuzzyDiff's tolerance).
#[test]
fn text_is_the_release_builds_output_line_for_line() {
    let (exp, features, mut ids) = tool_inputs();
    let input_file = data("QCCalculator_input.mzML");
    let input = MzQcInput {
        input_file: input_file.to_str().unwrap(),
        experiment: &exp,
        contact_name: "name",
        contact_address: "address",
        description: "description",
        label: "label",
        feature_map: &features,
    };
    let (text, report) = MzQcFile::new()
        .text(
            &input,
            &ids.protein_identifications,
            &mut ids.peptide_identifications,
            &DateTime::parse("2026-09-28T15:17:44").unwrap(),
        )
        .unwrap();
    assert!(report.log.is_empty(), "{:?}", report.log);
    let release = std::fs::read_to_string(data("release/QCCalculator_2_release.mzQC")).unwrap();
    let run_dependent = |line: &&str| !line.contains("\"location\"");
    let produced: Vec<&str> = text.split('\n').filter(run_dependent).collect();
    let expected: Vec<&str> = release.split('\n').filter(run_dependent).collect();
    assert_eq!(produced, expected);
    assert!(text.contains("\"value\": 0.916678249835968\n"));
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
