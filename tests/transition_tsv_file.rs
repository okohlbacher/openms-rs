// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! TransitionTSVFile: the upstream class test and tier-1 reproductions of the
//! `TargetedFileConverter` heavy-path TOPP tests (test-data `0cb15f2`,
//! `topp/CMakeLists.txt:1736-1864`) against their retained C++ outputs. See
//! `tests/data/transition_tsv_file_provenance.json`.

#![cfg(feature = "traml")]

use openms::analysis::openswath::targeted_experiment::TargetedExperiment;
use openms::analysis::openswath::transition_tsv_file::{
    HEADER_NAMES, RetentionTimeInterpretation, TransitionTSVFile,
};
use openms::format::FileType;
use openms::format::traml;
use openms::param::ParamValue;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/traml")
        .join(name)
}

fn temp(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms-tsv-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

/// Equal lines, or lines whose delimited tokens are equal or read as the same
/// `f64` (the retained outputs of older C++ releases spell some numbers
/// differently; the upstream comparison is FuzzyDiff, which accepts both).
fn same_line(a: &str, b: &str, separators: &[char]) -> bool {
    if a == b {
        return true;
    }
    let ta: Vec<&str> = a.split(separators).collect();
    let tb: Vec<&str> = b.split(separators).collect();
    ta.len() == tb.len()
        && ta.iter().zip(&tb).all(|(x, y)| {
            x == y || matches!((x.parse::<f64>(), y.parse::<f64>()), (Ok(p), Ok(q)) if p == q)
        })
}

/// Consecutive userParam lines sorted (registry order upstream, name order
/// here).
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

fn assert_same_document(name: &str, expected: &str, got: &str, separators: &[char]) {
    let (want, got) = (normalized_lines(expected), normalized_lines(got));
    for (index, (a, b)) in want.iter().zip(&got).enumerate() {
        assert!(
            same_line(a, b, separators),
            "{name}: line {}\n want: {a}\n  got: {b}",
            index + 1
        );
    }
    assert_eq!(want.len(), got.len(), "{name}: line count");
}

/// TargetedFileConverter heavy path, TSV/MRM to TraML
/// (`TargetedFileConverter.cpp` in openswath `e29e7ee`, lines 234-241 and
/// 274-276): read, validate, store.
fn tsv_to_traml(reader: &TransitionTSVFile, input: &str, file_type: FileType) -> String {
    let exp = reader
        .convert_tsv_to_targeted_experiment(data(input), file_type)
        .unwrap_or_else(|e| panic!("{input}: {e}"));
    reader.validate_targeted_experiment(&exp).unwrap();
    let mut out = Vec::new();
    traml::write(&mut out, &exp).unwrap();
    String::from_utf8(out).unwrap()
}

/// TraML to TSV (lines 257-266): load, convert.
fn traml_to_tsv(input: &str) -> String {
    let exp = traml::load(data(input)).unwrap();
    let mut out = Vec::new();
    TransitionTSVFile::new()
        .write_targeted_experiment(&mut out, &exp)
        .unwrap();
    String::from_utf8(out).unwrap()
}

// --- TransitionTSVFile_test.cpp ---------------------------------------------

#[test]
fn class_test_convert_targeted_experiment_to_tsv() {
    let exp = traml::load(data("MRMAssay_detectingTransistionCompound_input.TraML")).unwrap();
    let path = temp("class_out.tsv");
    TransitionTSVFile::new()
        .convert_targeted_experiment_to_tsv(&path, &exp)
        .unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.lines().count() > 1);
}

#[test]
fn class_test_convert_tsv_to_targeted_experiment() {
    let exp = traml::load(data("MRMAssay_detectingTransistionCompound_input.TraML")).unwrap();
    let path = temp("class_roundtrip.tsv");
    let handler = TransitionTSVFile::new();
    handler
        .convert_targeted_experiment_to_tsv(&path, &exp)
        .unwrap();
    let back = handler
        .convert_tsv_to_targeted_experiment(&path, FileType::Tsv)
        .unwrap();
    assert!(!back.transitions().is_empty());
    assert!(!back.compounds().is_empty() || !back.peptides().is_empty());
}

#[test]
fn defaults_match_the_source_parameters() {
    let defaults = TransitionTSVFile::defaults().unwrap();
    assert_eq!(
        defaults.value("retentionTimeInterpretation").unwrap(),
        &ParamValue::String("iRT".into())
    );
    assert_eq!(
        defaults.value("override_group_label_check").unwrap(),
        &ParamValue::String("false".into())
    );
    assert_eq!(
        defaults.value("force_invalid_mods").unwrap(),
        &ParamValue::String("false".into())
    );
    let mut param = defaults.clone();
    param
        .set_value(
            "retentionTimeInterpretation",
            ParamValue::String("minutes".into()),
            "",
            &[],
        )
        .unwrap();
    let handler = TransitionTSVFile::from_param(&param).unwrap();
    assert_eq!(
        handler.retention_time_interpretation,
        RetentionTimeInterpretation::Minutes
    );
    param
        .set_value(
            "retentionTimeInterpretation",
            ParamValue::String("hours".into()),
            "",
            &[],
        )
        .unwrap();
    assert!(TransitionTSVFile::from_param(&param).is_err());
}

// --- Tier 1: TargetedFileConverter TOPP tests --------------------------------

/// ConvertTSVToTraML tests 1, 5 and 6 (default parameters), 2 (minutes) and
/// 3 (seconds), 4 (SpectraST MRM) and TargetedFileConverter test 13.
#[test]
fn tsv_to_traml_reproduces_retained_outputs() {
    let default = TransitionTSVFile::new();
    let minutes = TransitionTSVFile {
        retention_time_interpretation: RetentionTimeInterpretation::Minutes,
        ..TransitionTSVFile::new()
    };
    let seconds = TransitionTSVFile {
        retention_time_interpretation: RetentionTimeInterpretation::Seconds,
        ..TransitionTSVFile::new()
    };
    for (reader, input, file_type, expected) in [
        (
            &default,
            "ConvertTSVToTraML_1_input.tsv",
            FileType::Tsv,
            "ConvertTSVToTraML_output.TraML",
        ),
        (
            &minutes,
            "ConvertTSVToTraML_2_input.tsv",
            FileType::Tsv,
            "ConvertTSVToTraML_2_output.TraML",
        ),
        (
            &seconds,
            "ConvertTSVToTraML_2_input.tsv",
            FileType::Tsv,
            "ConvertTSVToTraML_3_output.TraML",
        ),
        (
            &default,
            "ConvertTSVToTraML_4_input.mrm",
            FileType::Mrm,
            "ConvertTSVToTraML_4_output.TraML",
        ),
        (
            &default,
            "ConvertTSVToTraML_5_input.tsv",
            FileType::Tsv,
            "ConvertTSVToTraML_5_output.TraML",
        ),
        (
            &default,
            "ConvertTSVToTraML_6_input.tsv",
            FileType::Tsv,
            "ConvertTSVToTraML_6_output.TraML",
        ),
        (
            &default,
            "TargetedFileConverter_12_input.tsv",
            FileType::Tsv,
            "TargetedFileConverter_13_output.TraML",
        ),
    ] {
        let got = tsv_to_traml(reader, input, file_type);
        let want = std::fs::read_to_string(data(expected)).unwrap();
        assert_same_document(expected, &want, &got, &['"']);
    }
}

/// TargetedFileConverter tests 5 and 6: TraML to TSV.
#[test]
fn traml_to_tsv_reproduces_retained_outputs() {
    for (input, expected) in [
        (
            "ConvertTSVToTraML_output.TraML",
            "ConvertTraMLToTSV_output.tsv",
        ),
        (
            "ConvertTSVToTraML_5_output.TraML",
            "ConvertTraMLToTSV_output_2.tsv",
        ),
    ] {
        let got = traml_to_tsv(input);
        let want = std::fs::read_to_string(data(expected)).unwrap();
        assert_same_document(expected, &want, &got, &['\t']);
    }
}

/// TargetedFileConverter test 8: TSV -> TraML -> TSV with
/// `override_group_label_check` gives the input back.
#[test]
fn tsv_traml_tsv_round_trip_matches_test_8() {
    let reader = TransitionTSVFile {
        override_group_label_check: true,
        ..TransitionTSVFile::new()
    };
    let traml_text = tsv_to_traml(&reader, "TargetedFileConverter_8_input.tsv", FileType::Tsv);
    let exp: TargetedExperiment =
        traml::read_with_options(traml_text.as_bytes(), &traml::ReadOptions::default())
            .unwrap()
            .experiment;
    let mut out = Vec::new();
    reader.write_targeted_experiment(&mut out, &exp).unwrap();
    let got = String::from_utf8(out).unwrap();
    let want = std::fs::read_to_string(data("TargetedFileConverter_8_input.tsv")).unwrap();
    assert_same_document("TargetedFileConverter_8", &want, &got, &['\t']);
}

// --- Native behaviour ---------------------------------------------------------

#[test]
fn header_detection_follows_the_source() {
    let (delimiter, header) = TransitionTSVFile::tsv_header("a;b;c;d").unwrap();
    assert_eq!(delimiter, ';');
    assert_eq!(header["d"], 3);
    let (delimiter, header) = TransitionTSVFile::tsv_header("\"a\",'b',c,d,").unwrap();
    assert_eq!(delimiter, ',');
    // Quotes are removed and a trailing delimiter adds no column.
    assert_eq!(header.len(), 4);
    assert!(header.contains_key("a") && header.contains_key("b"));
    let error = TransitionTSVFile::tsv_header("a\tb")
        .unwrap_err()
        .to_string();
    assert!(error.contains("delimiter ','"), "{error}");
    assert_eq!(HEADER_NAMES.len(), 29);
}

#[test]
fn a_short_line_and_a_missing_column_are_refused() {
    let reader = TransitionTSVFile::new();
    let short = "PrecursorMz\tProductMz\tLibraryIntensity\tiRT\n500\t600\t1\n";
    assert!(
        reader
            .read_transitions(short.as_bytes(), FileType::Tsv)
            .is_err()
    );
    let no_rt = "PrecursorMz\tProductMz\tLibraryIntensity\tOther\n500\t600\t1\tx\n";
    assert!(
        reader
            .read_transitions(no_rt.as_bytes(), FileType::Tsv)
            .is_err()
    );
}

#[test]
fn a_mixed_label_group_is_split_unless_overridden() {
    let text = "PrecursorMz\tProductMz\tLibraryIntensity\tiRT\tPeptideSequence\tPeptideGroupLabel\tTransitionGroupId\n\
                500\t600\t1\t1\tPEPTIDE\tg\tA\n\
                501\t601\t1\t1\tPEPTIDEK\tg\tB\n";
    let (exp, warnings) = TransitionTSVFile::new()
        .read_targeted_experiment(text.as_bytes(), FileType::Tsv)
        .unwrap();
    assert_eq!(exp.peptides()[1].peptide_group_label(), "B");
    assert_eq!(warnings.len(), 1);
    let overriding = TransitionTSVFile {
        override_group_label_check: true,
        ..TransitionTSVFile::new()
    };
    let (exp, _) = overriding
        .read_targeted_experiment(text.as_bytes(), FileType::Tsv)
        .unwrap();
    assert_eq!(exp.peptides()[1].peptide_group_label(), "g");
}

#[test]
fn an_experiment_with_invalid_references_is_not_written() {
    use openms::analysis::openswath::reaction_monitoring_transition::ReactionMonitoringTransition;
    let mut exp = TargetedExperiment::new();
    let mut transition = ReactionMonitoringTransition::default();
    transition.set_peptide_ref("missing");
    exp.add_transition(transition);
    let handler = TransitionTSVFile::new();
    assert!(handler.write_targeted_experiment(Vec::new(), &exp).is_err());
    assert!(handler.validate_targeted_experiment(&exp).is_err());
}
