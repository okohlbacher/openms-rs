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

// ---------------------------------------------------------------------------
// The LightTargetedExperiment members (group O integration)
// ---------------------------------------------------------------------------
//
// `../oracle/integ-O-tsvlight` (ibminode06, Release install
// `openms4-release-bc9cc12-c19e494-174b576`, run twice, byte-identical) loaded
// ten cases through `convertTSVToTargetedExperiment(Light)` and wrote each
// experiment with `convertLightTargetedExperimentToTSV`; the dump of the light
// experiment (`%.17g`) and the TSV are the expected values.

fn light_dump(
    exp: &openms::analysis::openswath::transition_experiment::LightTargetedExperiment,
) -> String {
    use openms::format::file_info::text_format::ostream_g;
    use std::fmt::Write as _;
    let g = |v: f64| ostream_g(v, 17);
    let mut out = String::new();
    let _ = writeln!(out, "n_transitions\t{}", exp.transitions.len());
    let _ = writeln!(out, "n_compounds\t{}", exp.compounds.len());
    let _ = writeln!(out, "n_proteins\t{}", exp.proteins.len());
    for (i, t) in exp.transitions.iter().enumerate() {
        let _ = writeln!(
            out,
            "transition.{i}\tname={}\tref={}\tlib={}\tprod={}\tprec={}\tim={}\tfcharge={}\tdecoy={}\tdet={}\tquant={}\tident={}\tnr={}\ttype={}\tann={}\tforms={}",
            t.transition_name,
            t.peptide_ref,
            g(t.library_intensity),
            g(t.product_mz),
            g(t.precursor_mz),
            g(t.precursor_im),
            t.fragment_charge,
            u8::from(t.decoy()),
            u8::from(t.is_detecting_transition()),
            u8::from(t.is_quantifying_transition()),
            u8::from(t.is_identifying_transition()),
            t.fragment_nr,
            t.fragment_type_str(),
            t.annotation(),
            t.peptidoforms.join("|"),
        );
    }
    for (i, c) in exp.compounds.iter().enumerate() {
        let mods: String = c
            .modifications
            .iter()
            .map(|m| format!("{}:{};", m.location, m.unimod_id))
            .collect();
        let _ = writeln!(
            out,
            "compound.{i}\tid={}\tim={}\trt={}\trts={}\trte={}\tcharge={}\tseq={}\tprots={}\tlabel={}\tgene={}\tsum={}\tname={}\tltype={}\tsmiles={}\tadducts={}\tmods={mods}",
            c.id,
            g(c.drift_time),
            g(c.rt),
            g(c.rt_start),
            g(c.rt_end),
            c.charge,
            c.sequence,
            c.protein_refs.join(";"),
            c.peptide_group_label,
            c.gene_name,
            c.sum_formula,
            c.compound_name,
            c.label_type,
            c.smiles,
            c.adducts,
        );
    }
    for (i, p) in exp.proteins.iter().enumerate() {
        let _ = writeln!(
            out,
            "protein.{i}\tid={}\tseq={}\tuniprot={}",
            p.id, p.sequence, p.uniprot_id
        );
    }
    out
}

/// `convertTSVToTargetedExperiment(Light)` and `convertLightTargetedExperimentToTSV`
/// against the Release build, for the TSV and SpectraST MRM inputs of the
/// heavy-path tests, with the three parameters varied once.
#[test]
fn light_path_matches_the_release_build() {
    let expected = std::fs::read_to_string(
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/transition_tsv_light/oracle.stdout"),
    )
    .unwrap();
    let cases = std::fs::read_to_string(
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/transition_tsv_light/cases.tsv"),
    )
    .unwrap();
    let mut sections: Vec<(String, String)> = Vec::new();
    for line in cases.lines().filter(|l| !l.is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        let mut reader = TransitionTSVFile::new();
        for kv in &f[3..] {
            let (key, value) = kv.split_once('=').unwrap();
            match key {
                "retentionTimeInterpretation" => {
                    reader.retention_time_interpretation = match value {
                        "minutes" => RetentionTimeInterpretation::Minutes,
                        "seconds" => RetentionTimeInterpretation::Seconds,
                        _ => RetentionTimeInterpretation::Irt,
                    }
                }
                "override_group_label_check" => reader.override_group_label_check = value == "true",
                "force_invalid_mods" => reader.force_invalid_mods = value == "true",
                other => panic!("unknown parameter {other}"),
            }
        }
        let file = data(f[1].rsplit('/').next().unwrap());
        let file_type = if f[2] == "MRM" {
            FileType::Mrm
        } else {
            FileType::Tsv
        };
        let mut text = String::new();
        match reader.convert_tsv_to_light_targeted_experiment(&file, file_type) {
            Ok((exp, _warnings)) => {
                text.push_str(&light_dump(&exp));
                text.push_str("--- tsv\n");
                let mut bytes = Vec::new();
                reader
                    .write_light_targeted_experiment(&mut bytes, &exp)
                    .unwrap();
                text.push_str(&String::from_utf8(bytes).unwrap());
            }
            Err(error) => text.push_str(&format!("ERROR\t{error}\n")),
        }
        sections.push((f[0].to_owned(), text));
    }
    let want: Vec<(String, String)> = expected
        .split("### ")
        .filter(|s| !s.is_empty())
        .map(|s| {
            let (name, body) = s.split_once('\n').unwrap();
            (name.to_owned(), body.to_owned())
        })
        .collect();
    assert_eq!(sections.len(), want.len());
    assert!(sections.len() >= 10);
    for ((name, got), (want_name, want_body)) in sections.iter().zip(&want) {
        assert_eq!(name, want_name);
        // the numbers of the TSV are `toStr(double)` text in both
        assert_eq!(got, want_body, "case {name}");
    }
}

/// TransitionTSVFile_test.cpp "[EXTRA] Light path TSV roundtrip": TraML -> heavy
/// -> TSV -> light -> TSV loads data and writes a non-trivial file.
#[test]
fn class_test_light_path_tsv_roundtrip() {
    let heavy = traml::load(data("MRMAssay_detectingTransistionCompound_input.TraML")).unwrap();
    let reader = TransitionTSVFile::new();
    let tsv = temp("light_roundtrip.tsv");
    reader
        .convert_targeted_experiment_to_tsv(&tsv, &heavy)
        .unwrap();
    let (light, _) = reader
        .convert_tsv_to_light_targeted_experiment(&tsv, FileType::Tsv)
        .unwrap();
    assert!(!light.transitions.is_empty());
    assert!(!light.compounds.is_empty());
    let out = temp("light_roundtrip_out.tsv");
    reader
        .convert_light_targeted_experiment_to_tsv(&out, &light)
        .unwrap();
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.lines().count() > 1);
    // one line per transition after the header
    assert_eq!(text.lines().count(), light.transitions.len() + 1);
}
