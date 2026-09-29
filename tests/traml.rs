// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! TraMLFile / TraMLHandler: the upstream class test and tier-1 comparisons
//! against TraML documents the C++ writer produced (retained TOPP outputs at
//! test-data `0cb15f2`). See `tests/data/traml_provenance.json`.

#![cfg(feature = "traml")]

use openms::Error;
use openms::analysis::openswath::reaction_monitoring_transition::{
    DecoyTransitionType, ReactionMonitoringTransition,
};
use openms::analysis::openswath::targeted_experiment::TargetedExperiment;
use openms::analysis::openswath::targeted_experiment_helper::{Peptide, Protein, RTType, RTUnit};
use openms::format::traml::{self, ReadOptions};
use openms::format::traml_handler::{from_xsd_string, read_document, write_xml_escape};
use openms::metadata::MetaValueData;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/traml")
        .join(name)
}

fn temp(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms-traml-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

fn to_text(exp: &TargetedExperiment) -> String {
    let mut out = Vec::new();
    traml::write(&mut out, exp).unwrap();
    String::from_utf8(out).unwrap()
}

/// Lines of a document with each run of consecutive `userParam` lines
/// sorted: the source writes a record's userParams in `MetaInfoRegistry`
/// order, this port by name, and nothing else about the bytes may differ.
fn normalized_lines(text: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut run: Vec<String> = Vec::new();
    for line in text.lines() {
        if line.trim_start().starts_with("<userParam ") {
            run.push(line.to_owned());
            continue;
        }
        run.sort();
        out.append(&mut run);
        out.push(line.to_owned());
    }
    run.sort();
    out.append(&mut run);
    out
}

// --- TraMLFile_test.cpp -----------------------------------------------------

/// `store`: ToyExample1 survives a store/load round trip unchanged.
#[test]
fn class_test_store_round_trips_toy_example() {
    let original = traml::load(data("ToyExample1.traML")).unwrap();
    let path = temp("toy_roundtrip.traML");
    traml::store(&path, &original).unwrap();
    let reloaded = traml::load(&path).unwrap();
    assert!(reloaded == original);
}

/// `store`: a default experiment reloads equal once the two CVs the writer
/// adds are taken over, as the class test does.
#[test]
fn class_test_store_minimal_experiment() {
    let mut minimal = TargetedExperiment::new();
    let path = temp("minimal.traML");
    traml::store(&path, &minimal).unwrap();
    let reloaded = traml::load(&path).unwrap();
    assert_eq!(reloaded.cvs().len(), 2);
    minimal.set_cvs(reloaded.cvs().to_vec());
    assert!(reloaded == minimal);
}

/// `store`: one default protein, peptide and transition.
#[test]
fn class_test_store_minimal_protein_peptide_transition() {
    let mut minimal = TargetedExperiment::new();
    minimal.add_protein(Protein::default());
    minimal.add_peptide(Peptide::default());
    minimal.add_transition(ReactionMonitoringTransition::default());
    let path = temp("minimal_ppt.traML");
    traml::store(&path, &minimal).unwrap();
    let reloaded = traml::load(&path).unwrap();
    minimal.set_cvs(reloaded.cvs().to_vec());
    assert!(reloaded == minimal);
}

/// `equal`, `assign` and `add`.
#[test]
fn class_test_equal_assign_add() {
    let original = traml::load(data("ToyExample1.traML")).unwrap();
    let second = traml::load(data("ToyExample1.traML")).unwrap();
    assert!(second == original);

    let assigned = original.clone();
    assert_eq!(original.target_cv_terms().terms().len(), 1);
    assert_eq!(assigned.target_cv_terms().terms().len(), 1);
    assert!(assigned == original);

    let mut added = TargetedExperiment::new();
    added.merge(&original).unwrap();
    assert!(added == original);
}

/// `isSemanticallyValid`: a stored empty experiment has no errors and no
/// warnings.
#[cfg(feature = "semantic-validation")]
#[test]
fn class_test_semantically_valid_empty_document() {
    let path = temp("empty_semantic.traML");
    traml::store(&path, &TargetedExperiment::new()).unwrap();
    let report = traml::is_semantically_valid(&path).unwrap();
    assert!(report.is_valid(), "{report:?}");
    assert!(report.errors.is_empty());
    assert!(report.warnings.is_empty());
}

// --- Tier 1: documents the C++ writer produced --------------------------------

/// Every retained TraML output written by the C++ `TraMLHandler::writeTo` is
/// reproduced line for line by loading it and storing it again; the only
/// permitted difference is the order of consecutive `userParam` lines.
#[test]
fn rewriting_retained_cpp_output_reproduces_it() {
    for name in [
        "ConvertTSVToTraML_output.TraML",
        "ConvertTSVToTraML_2_output.TraML",
        "ConvertTSVToTraML_3_output.TraML",
        "ConvertTSVToTraML_5_output.TraML",
        "ConvertTSVToTraML_6_output.TraML",
        "TargetedFileConverter_3_output.TraML",
        "TargetedFileConverter_4_output.TraML",
        "TargetedFileConverter_13_output.TraML",
        "FileMerger_9_output.traML",
        "OpenSwathDecoyGenerator_output.TraML",
    ] {
        let expected = std::fs::read_to_string(data(name)).unwrap();
        let exp = traml::load(data(name)).unwrap_or_else(|e| panic!("{name}: {e}"));
        let written = to_text(&exp);
        let (want, got) = (normalized_lines(&expected), normalized_lines(&written));
        for (index, (a, b)) in want.iter().zip(&got).enumerate() {
            assert_eq!(a, b, "{name}: line {}", index + 1);
        }
        assert_eq!(want.len(), got.len(), "{name}: line count");
        // And the rewritten document decodes to the same experiment.
        let reread = traml::read_with_options(written.as_bytes(), &ReadOptions::default())
            .unwrap()
            .experiment;
        assert!(reread == exp, "{name}: round trip");
    }
}

/// FileMerger_9: the two inputs appended (`operator+=`) and stored give the
/// retained C++ output.
#[test]
fn file_merger_9_merge_matches_retained_output() {
    let mut merged = traml::load(data("FileMerger_9_input1.traML")).unwrap();
    let second = traml::load(data("FileMerger_9_input2.traML")).unwrap();
    merged.merge(&second).unwrap();
    // FileMerger.cpp:400-404 then appends its own Software record; the
    // version is the one the retained output was produced with (the upstream
    // FuzzyDiff whitelists that line through its `id=`).
    merged.add_software(openms::metadata::Software {
        name: "FileMerger".into(),
        version: "2.3.0".into(),
        ..Default::default()
    });
    let expected = traml::load(data("FileMerger_9_output.traML")).unwrap();
    assert!(merged == expected);
}

/// Values of ConvertTSVToTraML_output read through the handler's routing.
#[test]
fn convert_tsv_output_decodes_as_the_source_routes_it() {
    let exp = traml::load(data("ConvertTSVToTraML_output.TraML")).unwrap();
    assert_eq!(exp.proteins().len(), 2);
    let peptide = &exp.peptides()[0];
    assert_eq!(peptide.id, "tr_gr>1");
    assert_eq!(peptide.charge_state(), Some(2));
    assert_eq!(peptide.peptide_group_label(), "tr_gr1");
    assert_eq!(peptide.drift_time(), 4.5);
    let rt = &peptide.rts[0];
    assert_eq!(rt.rt(), Some(0.44));
    assert_eq!(rt.retention_time_type, RTType::Irt);
    // An iRT term has no unit; the reader records "unknown".
    assert_eq!(rt.retention_time_unit, RTUnit::Unknown);
    let modified = &exp.peptides()[2];
    // TraML location 0 is the N-terminus (-1), 20 the C-terminus.
    assert_eq!(modified.mods[0].location, -1);
    assert_eq!(modified.mods[0].unimod_id, 5);
    assert_eq!(modified.mods[1].location, 19);
    let transition = &exp.transitions()[0];
    assert_eq!(
        transition.decoy_transition_type(),
        DecoyTransitionType::Target
    );
    assert!(!exp.contains_invalid_references());
}

// --- Native behaviour ---------------------------------------------------------

const HEAD: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<TraML version="1.0.0" xmlns="http://psi.hupo.org/ms/traml">"#;

#[test]
fn unknown_elements_are_refused_by_default_and_reported_in_source_mode() {
    let text = format!("{HEAD}<Bogus/></TraML>");
    let error = read_document(text.as_bytes(), &ReadOptions::default(), None).unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error}");
    let load = read_document(text.as_bytes(), &ReadOptions::source(), None).unwrap();
    assert_eq!(
        load.warnings,
        vec![
            "TraMLHandler: unknown tag opening: 'Bogus'".to_owned(),
            "TraMLHandler: unknown tag closing: 'Bogus'".to_owned()
        ]
    );
}

#[test]
fn product_user_params_land_on_the_transition_as_in_the_source() {
    let text = format!(
        "{HEAD}<TransitionList><Transition id=\"t\"><Precursor>\
         <cvParam cvRef=\"MS\" accession=\"MS:1000827\" name=\"isolation window target m/z\" value=\"500\"/>\
         </Precursor><Product><userParam name=\"k\" type=\"xsd:double\" value=\"1.5\"/>\
         <cvParam cvRef=\"MS\" accession=\"MS:1000041\" name=\"charge state\" value=\"2.9\"/>\
         </Product></Transition></TransitionList></TraML>"
    );
    let exp = read_document(text.as_bytes(), &ReadOptions::default(), None)
        .unwrap()
        .experiment;
    let transition = &exp.transitions()[0];
    assert_eq!(transition.precursor_mz(), 500.0);
    assert!(matches!(
        transition.cv_terms.metadata["k"].data(),
        MetaValueData::Float(x) if *x == 1.5
    ));
    // The source converts the product charge through toDouble into an int.
    assert_eq!(transition.product().charge_state(), Some(2));
}

#[test]
fn limits_are_checked() {
    let text = format!("{HEAD}</TraML>");
    let options = ReadOptions {
        max_input_bytes: 10,
        ..ReadOptions::default()
    };
    assert!(read_document(text.as_bytes(), &options, None).is_err());
    let options = ReadOptions {
        max_depth: 1,
        ..ReadOptions::default()
    };
    let nested = format!("{HEAD}<cvList/></TraML>");
    assert!(read_document(nested.as_bytes(), &options, None).is_err());
    assert!(read_document("<NotTraML/>".as_bytes(), &ReadOptions::default(), None).is_err());
    let doctype = "<!DOCTYPE TraML [<!ENTITY x \"y\">]><TraML/>";
    assert!(read_document(doctype.as_bytes(), &ReadOptions::default(), None).is_err());
}

#[test]
fn xsd_values_and_escaping_follow_the_source_helpers() {
    assert!(matches!(
        from_xsd_string("xsd:int", "42").unwrap().data(),
        MetaValueData::Integer(42)
    ));
    assert!(matches!(
        from_xsd_string("xsd:long", "-7").unwrap().data(),
        MetaValueData::Integer(-7)
    ));
    assert!(matches!(
        from_xsd_string("xsd:boolean", "true").unwrap().data(),
        MetaValueData::String(s) if s == "true"
    ));
    assert!(from_xsd_string("xsd:int", "x").is_err());
    assert_eq!(
        write_xml_escape("a&b<c>\"d'"),
        "a&amp;b&lt;c&gt;&quot;d&apos;"
    );
}

#[test]
fn a_retention_time_without_a_value_is_refused_on_store() {
    use openms::analysis::openswath::targeted_experiment_helper::RetentionTime;
    let mut exp = TargetedExperiment::new();
    let mut peptide = Peptide::default();
    peptide.rts.push(RetentionTime::default());
    exp.add_peptide(peptide);
    let path = temp("unset_rt.traML");
    assert!(traml::store(&path, &exp).is_err());
    assert!(!path.exists());
}
