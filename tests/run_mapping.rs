// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::identification::{
    IdentifierMSRunMapper, PeptideIdentification, ProteinHit, ProteinIdentification,
};
use openms::metadata::{MetaValue, MetaValueData};
use std::borrow::Cow;

fn run(identifier: &str, paths: &[&str]) -> ProteinIdentification {
    ProteinIdentification {
        identifier: identifier.into(),
        primary_ms_run_paths: paths.iter().map(|s| s.to_string()).collect(),
        ..Default::default()
    }
}
fn peptide(identifier: &str) -> PeptideIdentification {
    PeptideIdentification {
        identifier: identifier.into(),
        ..Default::default()
    }
}
fn source_runs() -> Vec<ProteinIdentification> {
    vec![
        run("ID_A", &["/data/runA.mzML"]),
        run("ID_B", &["/data/runB1.mzML", "/data/runB2.mzML"]),
    ]
}
fn strings(paths: &[&str]) -> Vec<String> {
    paths.iter().map(|s| s.to_string()).collect()
}

#[test]
fn source_lookup_literals_and_complete_replacement() {
    let mut m = IdentifierMSRunMapper::new();
    assert!(m.is_empty());
    assert_eq!(m.len(), 0);
    m.create(&source_runs()).unwrap();
    assert_eq!(m.len(), 2);
    assert_eq!(m.identifiers().collect::<Vec<_>>(), ["ID_A", "ID_B"]);
    assert!(m.has_identifier("ID_A"));
    assert!(!m.has_identifier("ID_X"));
    assert_eq!(m.ms_run_paths("ID_A"), ["/data/runA.mzML"]);
    assert_eq!(
        m.ms_run_paths("ID_B"),
        ["/data/runB1.mzML", "/data/runB2.mzML"]
    );
    assert!(m.ms_run_paths("ID_X").is_empty());
    let path = strings(&["/data/runA.mzML"]);
    assert!(m.has_run_path(&path));
    assert_eq!(m.identifier(&path).unwrap(), "ID_A");
    assert_eq!(m.try_identifier(&path), Some("ID_A"));
    assert_eq!(
        m.identifier(&strings(&["/data/runB1.mzML", "/data/runB2.mzML"]))
            .unwrap(),
        "ID_B"
    );
    assert!(m.identifier(&strings(&["/data/runB1.mzML"])).is_err());
    assert!(!m.has_run_path(&strings(&["/data/runB1.mzML"])));
    assert_eq!(m.try_identifier(&strings(&["/missing.mzML"])), None);
    let copy = m.clone();
    m.create(&[run("ONLY", &["only.mzML"])]).unwrap();
    assert_eq!(m.len(), 1);
    assert!(!m.has_identifier("ID_A"));
    assert_eq!(copy.len(), 2);
    m.create(&[]).unwrap();
    assert!(m.is_empty());
}

#[test]
fn source_default_index_explicit_index_and_legacy_fallback() {
    let m = IdentifierMSRunMapper::from_runs(&source_runs()).unwrap();
    let mut p = peptide("ID_B");
    assert_eq!(m.primary_ms_run_path(&p).unwrap(), "/data/runB1.mzML");
    assert!(matches!(
        m.primary_ms_run_path(&p).unwrap(),
        Cow::Borrowed(_)
    ));
    p.metadata.insert("id_merge_index".into(), 1i64.into());
    assert_eq!(m.primary_ms_run_path(&p).unwrap(), "/data/runB2.mzML");
    p.metadata.insert("id_merge_index".into(), 42i64.into());
    assert_eq!(m.primary_ms_run_path(&p).unwrap(), "");
    p.metadata
        .insert("base_name".into(), "legacy_run.mzML".into());
    assert_eq!(m.primary_ms_run_path(&p).unwrap(), "legacy_run.mzML");
    let mut both = peptide("ID_A");
    both.metadata
        .insert("base_name".into(), "legacy_run.mzML".into());
    assert_eq!(m.primary_ms_run_path(&both).unwrap(), "/data/runA.mzML");
    let mut missing = peptide("missing");
    assert_eq!(m.primary_ms_run_path(&missing).unwrap(), "");
    missing
        .metadata
        .insert("base_name".into(), "legacy_run.mzML".into());
    missing
        .metadata
        .insert("id_merge_index".into(), "unused wrong type".into());
    assert_eq!(m.primary_ms_run_path(&missing).unwrap(), "legacy_run.mzML");
}

#[test]
fn negative_and_wrong_type_indices_error_before_fallback_only_when_consumed() {
    let m = IdentifierMSRunMapper::from_runs(&[run("one", &["a"]), run("empty", &[])]).unwrap();
    for value in [
        (-1i64).into(),
        "0".into(),
        MetaValue::new(MetaValueData::Float(0.0)).unwrap(),
        MetaValue::default(),
    ] {
        let mut p = peptide("one");
        p.metadata.insert("id_merge_index".into(), value.clone());
        p.metadata.insert("base_name".into(), "fallback".into());
        assert!(m.primary_ms_run_path(&p).is_err());
        // The source validator deliberately exempts single-file runs.
        m.validate_merge_index(&p, 7).unwrap();
        p.identifier = "empty".into();
        assert_eq!(m.primary_ms_run_path(&p).unwrap(), "fallback");
        m.validate_merge_index(&p, 7).unwrap();
    }
}

#[test]
fn source_duplicate_error_publishes_all_forward_and_only_prefix_reverse_records() {
    let mut m = IdentifierMSRunMapper::from_runs(&[run("old", &["old"])]).unwrap();
    let runs = [
        run("A", &["a"]),
        run("B", &["shared"]),
        run("C", &["shared"]),
        run("D", &["d"]),
    ];
    assert!(m.create(&runs).is_err());
    assert_eq!(m.identifiers().collect::<Vec<_>>(), ["A", "B", "C", "D"]);
    assert_eq!(m.primary_ms_run_path(&peptide("D")).unwrap(), "d");
    assert_eq!(m.primary_ms_run_path(&peptide("C")).unwrap(), "shared");
    assert!(!m.has_identifier("old"));
    assert_eq!(m.identifier(&strings(&["shared"])).unwrap(), "B");
    assert!(m.has_run_path(&strings(&["a"])));
    assert!(!m.has_run_path(&strings(&["d"])));
    assert!(IdentifierMSRunMapper::from_runs(&runs).is_err());
    m.create(&[run("recovered", &["new"])]).unwrap();
    assert_eq!(m.identifiers().collect::<Vec<_>>(), ["recovered"]);
}

#[test]
fn repeated_identifiers_and_empty_path_lists_follow_two_distinct_source_maps() {
    let mut m =
        IdentifierMSRunMapper::from_runs(&[run("A", &["first"]), run("A", &["last"])]).unwrap();
    assert_eq!(m.len(), 1);
    assert_eq!(m.ms_run_paths("A"), ["last"]);
    assert_eq!(m.identifier(&strings(&["first"])).unwrap(), "A");
    assert_eq!(m.identifier(&strings(&["last"])).unwrap(), "A");
    assert!(
        m.create(&[run("A", &["same"]), run("A", &["same"])])
            .is_err()
    );
    assert_eq!(m.len(), 1);
    assert_eq!(m.identifier(&strings(&["same"])).unwrap(), "A");
    assert!(m.create(&[run("A", &[]), run("B", &[])]).is_err());
    assert_eq!(m.len(), 2);
    assert_eq!(m.identifier(&[]).unwrap(), "A");
    assert!(m.has_identifier("B"));
}

#[test]
fn merge_validator_rejects_only_ambiguous_multi_file_assignments() {
    let m = IdentifierMSRunMapper::from_runs(&[
        run("merged", &["a", "b"]),
        run("single", &["c"]),
        run("none", &[]),
    ])
    .unwrap();
    let mut p = peptide("merged");
    let error = m.validate_merge_index(&p, 123).unwrap_err().to_string();
    assert!(error.contains("123") && error.contains("missing id_merge_index"));
    for value in [
        (-1i64).into(),
        2i64.into(),
        i64::MAX.into(),
        "1".into(),
        MetaValue::new(MetaValueData::Float(1.0)).unwrap(),
        MetaValue::default(),
    ] {
        p.metadata.insert("id_merge_index".into(), value);
        assert!(m.validate_merge_index(&p, 1).is_err());
    }
    for value in [0i64, 1] {
        p.metadata.insert("id_merge_index".into(), value.into());
        m.validate_merge_index(&p, 1).unwrap();
    }
    p.metadata
        .insert("id_merge_index".into(), "not an integer".into());
    for id in ["single", "none", "unknown"] {
        p.identifier = id.into();
        m.validate_merge_index(&p, usize::MAX).unwrap();
    }
}

#[test]
fn lexical_paths_preserve_order_duplicates_spelling_and_unicode() {
    let m = IdentifierMSRunMapper::from_runs(&[
        run("α", &["./a", "a", "a"]),
        run("z", &["a", "./a", "a"]),
        run("", &[" C:\\λ "]),
    ])
    .unwrap();
    assert_eq!(m.identifiers().collect::<Vec<_>>(), ["", "z", "α"]);
    assert_eq!(m.identifier(&strings(&["./a", "a", "a"])).unwrap(), "α");
    assert_eq!(m.identifier(&strings(&["a", "./a", "a"])).unwrap(), "z");
    assert_eq!(m.primary_ms_run_path(&peptide("")).unwrap(), " C:\\λ ");
}

#[test]
fn lenient_legacy_fallback_has_source_scalar_and_list_formatting() {
    let m = IdentifierMSRunMapper::new();
    let mut p = peptide("legacy");
    for (value, expected) in [
        (42i64.into(), "42"),
        (MetaValue::new(MetaValueData::Float(-0.0)).unwrap(), "-0.0"),
        (
            MetaValue::new(MetaValueData::FloatList(vec![1e-5, 2.0])).unwrap(),
            "[1.0e-05, 2.0]",
        ),
        (
            MetaValue::new(MetaValueData::IntegerList(vec![-2, 3])).unwrap(),
            "[-2, 3]",
        ),
        (vec!["a,b".to_string(), "λ".into()].into(), "[a,b, λ]"),
        (MetaValue::default(), ""),
    ] {
        p.metadata.insert("base_name".into(), value);
        assert_eq!(m.primary_ms_run_path(&p).unwrap(), expected);
    }
}

#[test]
fn construction_ignores_unrelated_payload_and_owns_input_paths() {
    let mut input = source_runs();
    input[0].hits.push(ProteinHit {
        score: f64::NAN,
        ..Default::default()
    });
    input[0].raw_ms_run_paths.push("unused.raw".into());
    input[0]
        .metadata
        .insert("spectra_data".into(), "unused legacy metadata".into());
    let m = IdentifierMSRunMapper::from_runs(&input).unwrap();
    input[0].primary_ms_run_paths[0].push('x');
    assert_eq!(m.ms_run_paths("ID_A"), ["/data/runA.mzML"]);
}

#[test]
fn resource_failures_preserve_previous_mapping_before_clone() {
    let mut m = IdentifierMSRunMapper::from_runs(&source_runs()).unwrap();
    let before = m.clone();
    let mut oversized = run("huge", &[]);
    oversized.primary_ms_run_paths = vec![String::new(); IdentifierMSRunMapper::MAX_ITEMS];
    assert!(m.create(&[oversized]).is_err());
    assert_eq!(m, before);
    let mut oversized = run("huge", &[]);
    oversized
        .primary_ms_run_paths
        .push("x".repeat(IdentifierMSRunMapper::MAX_BYTES / 2));
    assert!(m.create(&[oversized]).is_err());
    assert_eq!(m, before);
}
