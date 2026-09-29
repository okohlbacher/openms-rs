// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! TransitionPQPFile: tier-1 reproductions of the `TargetedFileConverter`
//! PQP tests (test-data `0cb15f2`, `topp/CMakeLists.txt:1762-1851`) and the
//! heavy-path sections of `TransitionPQPFile_test.cpp`. See
//! `tests/data/transition_pqp_file_provenance.json`.

#![cfg(all(feature = "sqlite", feature = "traml"))]

use openms::analysis::openswath::targeted_experiment::TargetedExperiment;
use openms::analysis::openswath::transition_pqp_file::TransitionPQPFile;
use openms::analysis::openswath::transition_tsv_file::TransitionTSVFile;
use openms::format::traml;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/traml")
        .join(name)
}

fn temp(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms-pqp-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

/// Lines agree when equal after trimming surrounding whitespace (FuzzyDiff
/// splits on whitespace; the hand-edited `TargetedFileConverter_11_input`
/// indents one `</Product>` by five spaces), or when every token is equal or
/// reads as the same `f64`.
fn same_line(a: &str, b: &str, separators: &[char]) -> bool {
    let (a, b) = (a.trim(), b.trim());
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

/// Compare as the upstream FuzzyDiff does, with numbers required to be the
/// identical double; lines containing a whitelisted term are skipped on both
/// sides, as FuzzyDiff's `-whitelist` does.
fn assert_same_document(
    name: &str,
    expected: &str,
    got: &str,
    separators: &[char],
    whitelist: &[&str],
) {
    let keep = |lines: Vec<String>| -> Vec<String> {
        lines
            .into_iter()
            .filter(|l| !whitelist.iter().any(|w| l.contains(w)))
            .collect()
    };
    let (want, got) = (
        keep(normalized_lines(expected)),
        keep(normalized_lines(got)),
    );
    for (index, (a, b)) in want.iter().zip(&got).enumerate() {
        assert!(
            same_line(a, b, separators),
            "{name}: line {}\n want: {a}\n  got: {b}",
            index + 1
        );
    }
    assert_eq!(want.len(), got.len(), "{name}: line count");
}

fn traml_text(exp: &TargetedExperiment) -> String {
    let mut out = Vec::new();
    traml::write(&mut out, exp).unwrap();
    String::from_utf8(out).unwrap()
}

/// `TargetedFileConverter` TraML -> PQP (openswath `e29e7ee`,
/// `TargetedFileConverter.cpp:257-273`), then PQP -> TraML (243-251, 274-276).
fn traml_pqp_traml(
    pqp: &TransitionPQPFile,
    input: &str,
    scratch: &str,
    legacy_traml_id: bool,
) -> String {
    let exp = traml::load(data(input)).unwrap();
    let path = temp(scratch);
    pqp.convert_targeted_experiment_to_pqp(&path, &exp).unwrap();
    let (back, _) = pqp
        .convert_pqp_to_targeted_experiment(&path, legacy_traml_id)
        .unwrap();
    pqp.validate_targeted_experiment(&back).unwrap();
    traml_text(&back)
}

/// TargetedFileConverter tests 3 (numeric ids) and 4 (`-legacy_traml_id`):
/// metabolomics TraML through PQP back to TraML.
#[test]
fn metabolomics_round_trip_reproduces_tests_3_and_4() {
    let pqp = TransitionPQPFile::new();
    let got = traml_pqp_traml(&pqp, "TargetedFileConverter_3_input.TraML", "t3.pqp", false);
    let want = std::fs::read_to_string(data("TargetedFileConverter_3_output.TraML")).unwrap();
    assert_same_document("TargetedFileConverter_3", &want, &got, &['"'], &[]);
    let got = traml_pqp_traml(&pqp, "TargetedFileConverter_3_input.TraML", "t4.pqp", true);
    let want = std::fs::read_to_string(data("TargetedFileConverter_4_output.TraML")).unwrap();
    assert_same_document("TargetedFileConverter_4", &want, &got, &['"'], &[]);
}

/// TargetedFileConverter tests 11 and 12: TraML through PQP gives the input
/// back, ignoring `id=` lines as the upstream whitelist does.
#[test]
fn proteomics_round_trip_reproduces_tests_11_and_12() {
    let pqp = TransitionPQPFile::new();
    let got = traml_pqp_traml(
        &pqp,
        "TargetedFileConverter_11_input.TraML",
        "t11.pqp",
        false,
    );
    let want = std::fs::read_to_string(data("TargetedFileConverter_11_input.TraML")).unwrap();
    assert_same_document("TargetedFileConverter_11", &want, &got, &['"'], &["id="]);

    let overriding = TransitionPQPFile {
        tsv: TransitionTSVFile {
            override_group_label_check: true,
            ..TransitionTSVFile::new()
        },
    };
    let got = traml_pqp_traml(
        &overriding,
        "TargetedFileConverter_12_input.TraML",
        "t12.pqp",
        false,
    );
    let want = std::fs::read_to_string(data("TargetedFileConverter_12_input.TraML")).unwrap();
    assert_same_document("TargetedFileConverter_12", &want, &got, &['"'], &["id="]);
}

/// TargetedFileConverter test 10: a PQP the C++ writer produced, read with
/// legacy ids into TraML and on to TSV, gives the retained TSV.
#[test]
fn cpp_written_pqp_reads_as_test_10() {
    let pqp = TransitionPQPFile::new();
    let (exp, _) = pqp
        .convert_pqp_to_targeted_experiment(data("TargetedFileConverter_10_input.pqp"), true)
        .unwrap();
    // The upstream test stores TraML and reads it back before writing TSV.
    let reread = traml::read_with_options(traml_text(&exp).as_bytes(), &Default::default())
        .unwrap()
        .experiment;
    let mut out = Vec::new();
    TransitionTSVFile::new()
        .write_targeted_experiment(&mut out, &reread)
        .unwrap();
    let got = String::from_utf8(out).unwrap();
    let want = std::fs::read_to_string(data("TargetedFileConverter_10_output.tsv")).unwrap();
    assert_same_document("TargetedFileConverter_10", &want, &got, &['\t'], &[]);
}

/// Native: the two id maps cover both tables of a written file, and a table
/// name other than PRECURSOR or TRANSITION is refused (as
/// `TransitionPQPFile_test.cpp:313-316` checks on its light-path file).
#[test]
fn id_maps_cover_precursors_and_transitions() {
    let exp = traml::load(data("TargetedFileConverter_3_input.TraML")).unwrap();
    let path = temp("idmap.pqp");
    TransitionPQPFile::new()
        .convert_targeted_experiment_to_pqp(&path, &exp)
        .unwrap();
    let transitions = TransitionPQPFile::pqp_id_to_traml_id_map(&path, "TRANSITION").unwrap();
    assert_eq!(transitions.len(), exp.transitions().len());
    for transition in exp.transitions() {
        assert!(transitions.contains_key(transition.native_id()));
    }
    let current = TransitionPQPFile::pqp_current_id_to_traml_id_map(&path, "PRECURSOR").unwrap();
    assert_eq!(current.len(), exp.compounds().len() + exp.peptides().len());
    assert!(TransitionPQPFile::pqp_id_to_traml_id_map(&path, "PEPTIDE").is_err());
    assert!(TransitionPQPFile::pqp_current_id_to_traml_id_map(&path, "GENE").is_err());
}

#[test]
fn an_experiment_with_invalid_references_is_not_written() {
    use openms::analysis::openswath::reaction_monitoring_transition::ReactionMonitoringTransition;
    let mut exp = TargetedExperiment::new();
    let mut transition = ReactionMonitoringTransition::default();
    transition.set_compound_ref("missing");
    exp.add_transition(transition);
    let path = temp("invalid.pqp");
    assert!(
        TransitionPQPFile::new()
            .convert_targeted_experiment_to_pqp(&path, &exp)
            .is_err()
    );
    assert!(!path.exists());
}

const SCHEMA: &str = "CREATE TABLE VERSION(ID INT NOT NULL);\
CREATE TABLE GENE(ID INT PRIMARY KEY NOT NULL, GENE_NAME TEXT NOT NULL, DECOY INT NOT NULL);\
CREATE TABLE PEPTIDE_GENE_MAPPING(PEPTIDE_ID INT NOT NULL, GENE_ID INT NOT NULL);\
CREATE TABLE PROTEIN(ID INT PRIMARY KEY NOT NULL, PROTEIN_ACCESSION TEXT NOT NULL, DECOY INT NOT NULL);\
CREATE TABLE PEPTIDE_PROTEIN_MAPPING(PEPTIDE_ID INT NOT NULL, PROTEIN_ID INT NOT NULL);\
CREATE TABLE PEPTIDE(ID INT PRIMARY KEY NOT NULL, UNMODIFIED_SEQUENCE TEXT NOT NULL, MODIFIED_SEQUENCE TEXT NOT NULL, DECOY INT NOT NULL);\
CREATE TABLE PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID INT NOT NULL, PEPTIDE_ID INT NOT NULL);\
CREATE TABLE COMPOUND(ID INT PRIMARY KEY NOT NULL, COMPOUND_NAME TEXT NOT NULL, SUM_FORMULA TEXT NOT NULL, SMILES TEXT NOT NULL, ADDUCTS TEXT NOT NULL, DECOY INT NOT NULL);\
CREATE TABLE PRECURSOR_COMPOUND_MAPPING(PRECURSOR_ID INT NOT NULL, COMPOUND_ID INT NOT NULL);\
CREATE TABLE PRECURSOR(ID INT PRIMARY KEY NOT NULL, TRAML_ID TEXT NULL, GROUP_LABEL TEXT NULL, PRECURSOR_MZ REAL NOT NULL, CHARGE INT NULL, LIBRARY_INTENSITY REAL NULL, LIBRARY_RT REAL NULL, LIBRARY_DRIFT_TIME REAL NULL, DECOY INT NOT NULL);\
CREATE TABLE TRANSITION_PRECURSOR_MAPPING(TRANSITION_ID INT NOT NULL, PRECURSOR_ID INT NOT NULL);\
CREATE TABLE TRANSITION_PEPTIDE_MAPPING(TRANSITION_ID INT NOT NULL, PEPTIDE_ID INT NOT NULL);\
CREATE TABLE TRANSITION(ID INT PRIMARY KEY NOT NULL, TRAML_ID TEXT NULL, PRODUCT_MZ REAL NOT NULL, CHARGE INT NULL, TYPE CHAR(1) NULL, ANNOTATION TEXT NULL, ORDINAL INT NULL, DETECTING INT NOT NULL, IDENTIFYING INT NOT NULL, QUANTIFYING INT NOT NULL, LIBRARY_INTENSITY REAL NULL, DECOY INT NOT NULL);";

/// The minimal database of `TransitionPQPFile_test.cpp:80-126`, plus the
/// statements given.
fn minimal_pqp(name: &str, extra: &[&str]) -> PathBuf {
    use openms::format::sqlite_connector::SqliteConnector;
    let path = temp(name);
    let _ = std::fs::remove_file(&path);
    let conn = SqliteConnector::new(&path).unwrap();
    conn.execute_statement(SCHEMA).unwrap();
    for statement in [
        "INSERT INTO VERSION (ID) VALUES (3);",
        "INSERT INTO PROTEIN (ID, PROTEIN_ACCESSION, DECOY) VALUES (0, 'TestProtein', 0);",
        "INSERT INTO PEPTIDE (ID, UNMODIFIED_SEQUENCE, MODIFIED_SEQUENCE, DECOY) VALUES (0, 'PEPTIDE', 'PEPTIDE', 0);",
        "INSERT INTO PEPTIDE_PROTEIN_MAPPING (PEPTIDE_ID, PROTEIN_ID) VALUES (0, 0);",
        "INSERT INTO PRECURSOR (ID, TRAML_ID, GROUP_LABEL, PRECURSOR_MZ, CHARGE, LIBRARY_RT, LIBRARY_DRIFT_TIME, DECOY) VALUES (0, 'precursor_1', 'group1', 500.0, 2, 100.0, -1, 0);",
        "INSERT INTO PRECURSOR_PEPTIDE_MAPPING (PRECURSOR_ID, PEPTIDE_ID) VALUES (0, 0);",
        "INSERT INTO TRANSITION (ID, TRAML_ID, PRODUCT_MZ, CHARGE, TYPE, ORDINAL, DETECTING, IDENTIFYING, QUANTIFYING, LIBRARY_INTENSITY, DECOY) VALUES (0, 'transition_1', 300.0, 1, 'y', 3, 1, 0, 1, 1000.0, 0);",
        "INSERT INTO TRANSITION_PRECURSOR_MAPPING (TRANSITION_ID, PRECURSOR_ID) VALUES (0, 0);",
    ]
    .iter()
    .chain(extra)
    {
        conn.execute_statement(statement).unwrap();
    }
    path
}

/// `TransitionPQPFile_test.cpp:71-137`: an empty GENE table (issue #8687)
/// must not empty the result through its inner join.
#[test]
fn class_test_empty_gene_table() {
    let path = minimal_pqp("empty_gene.pqp", &[]);
    let (exp, _) = TransitionPQPFile::new()
        .convert_pqp_to_targeted_experiment(&path, true)
        .unwrap();
    assert_eq!(exp.transitions().len(), 1);
    assert_eq!(exp.peptides().len(), 1);
    assert_eq!(exp.proteins().len(), 1);
    assert_eq!(exp.transitions()[0].native_id(), "transition_1");
    assert_eq!(exp.peptides()[0].id, "precursor_1");
}

/// `TransitionPQPFile_test.cpp:163-206`: two genes on one peptide do not
/// duplicate its transition.
#[test]
fn class_test_multiple_gene_mappings() {
    let path = minimal_pqp(
        "two_genes.pqp",
        &[
            "INSERT INTO GENE (ID, GENE_NAME, DECOY) VALUES (0, 'GeneA', 0);",
            "INSERT INTO GENE (ID, GENE_NAME, DECOY) VALUES (1, 'GeneB', 0);",
            "INSERT INTO PEPTIDE_GENE_MAPPING (PEPTIDE_ID, GENE_ID) VALUES (0, 0);",
            "INSERT INTO PEPTIDE_GENE_MAPPING (PEPTIDE_ID, GENE_ID) VALUES (0, 1);",
        ],
    );
    let (exp, _) = TransitionPQPFile::new()
        .convert_pqp_to_targeted_experiment(&path, true)
        .unwrap();
    assert_eq!(exp.transitions().len(), 1);
    assert_eq!(exp.peptides().len(), 1);
}
