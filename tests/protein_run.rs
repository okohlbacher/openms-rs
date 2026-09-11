// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::comparison::Tolerance;
use openms::identification::{
    EnzymeTermSpecificity, ProteinGroup, ProteinHit, ProteinIdentification,
};
use openms::metadata::{MetaValue, MetaValueData, Unit};

fn run(engine: &str) -> ProteinIdentification {
    ProteinIdentification {
        search_engine: engine.into(),
        search_engine_version: "v1".into(),
        ..Default::default()
    }
}
fn hit(accession: &str, score: f64) -> ProteinHit {
    ProteinHit {
        accession: accession.into(),
        score,
        ..Default::default()
    }
}

#[test]
fn inference_engine_heuristics_and_explicit_overrides_are_source_exact() {
    for engine in [
        "Fido",
        "BayesianProteinInference",
        "Epifany",
        "ProteinInference",
    ] {
        let r = run(engine);
        assert!(r.has_inference_engine_as_search_engine());
        assert_eq!(r.inference_engine().unwrap(), engine);
        assert_eq!(r.inference_engine_version().unwrap(), "v1");
        assert!(r.has_inference_data().unwrap());
    }
    for engine in [
        "",
        "Comet",
        "Percolator",
        "fido",
        "ConsensusID",
        "FidoAdapter",
    ] {
        let r = run(engine);
        assert!(!r.has_inference_engine_as_search_engine());
        assert_eq!(r.inference_engine().unwrap(), "");
        assert_eq!(r.inference_engine_version().unwrap(), "");
    }
    let mut r = run("Percolator");
    r.indistinguishable_groups.push(ProteinGroup::default());
    assert!(r.has_inference_data().unwrap());
    r.set_inference_engine("");
    assert!(!r.has_inference_data().unwrap());
    assert_eq!(r.inference_engine_version().unwrap(), "");
    r.set_inference_engine("Custom");
    assert_eq!(r.inference_engine_version().unwrap(), "v1");
    r.set_inference_engine_version("v2");
    assert_eq!(r.inference_engine_version().unwrap(), "v2");
    r.set_inference_engine("");
    assert_eq!(r.inference_engine_version().unwrap(), "v2");
}

#[test]
fn inference_values_use_strict_types_but_explicit_version_needs_no_engine() {
    let mut r = run("Fido");
    r.search_parameters
        .metadata
        .insert("InferenceEngine".into(), 2i64.into());
    assert!(r.inference_engine().is_err());
    assert!(r.has_inference_data().is_err());
    assert!(r.inference_engine_version().is_err());
    r.set_inference_engine_version("independent");
    assert_eq!(r.inference_engine_version().unwrap(), "independent");
    r.search_parameters
        .metadata
        .insert("InferenceEngineVersion".into(), MetaValue::default());
    assert!(r.inference_engine_version().is_err());
}

#[test]
fn original_engine_uses_key_order_case_and_substrings() {
    let mut r = run("Comet");
    assert_eq!(r.original_search_engine_name(), "Comet");
    r.search_engine = "preConsensusIDsuffix".into();
    assert_eq!(r.original_search_engine_name(), "Unknown");
    for key in ["SE:percolator", "SE:Zeta", "SE:Alpha", "SE:ApercolatorB"] {
        r.search_parameters
            .metadata
            .insert(key.into(), 42i64.into());
    }
    assert_eq!(r.original_search_engine_name(), "Alpha");
    r.search_parameters.metadata.insert("SE:".into(), "".into());
    assert_eq!(r.original_search_engine_name(), "");
    r.search_parameters.metadata.clear();
    r.search_parameters
        .metadata
        .insert("SE:Percolator".into(), "".into());
    assert_eq!(r.original_search_engine_name(), "Percolator");
    r.search_engine = "consensusid".into();
    assert_eq!(r.original_search_engine_name(), "consensusid");
}

#[test]
fn standard_settings_keep_source_order_units_and_modification_lists() {
    let mut r = run("Comet");
    let p = &mut r.search_parameters;
    p.database = "/db/proteins.fasta".into();
    p.database_version = "v7".into();
    p.fragment_tolerance = Tolerance::Absolute(0.02);
    p.precursor_tolerance = Tolerance::Ppm(20.0);
    p.digestion_enzyme = "Trypsin".into();
    p.enzyme_specificity = EnzymeTermSpecificity::Full;
    p.charges = "2,3".into();
    p.missed_cleavages = 2;
    p.fixed_modifications = vec!["B".into(), "A".into(), "B".into()];
    let expected = [
        ("db", "/db/proteins.fasta"),
        ("db_version", "v7"),
        ("fragment_mass_tolerance", "0.02"),
        ("fragment_mass_tolerance_unit", "Da"),
        ("precursor_mass_tolerance", "20.0"),
        ("precursor_mass_tolerance_unit", "ppm"),
        ("enzyme", "Trypsin"),
        ("enzyme_term_specificity", "full"),
        ("charges", "2,3"),
        ("missed_cleavages", "2"),
        ("fixed_modifications", "B,A,B"),
        ("variable_modifications", ""),
    ]
    .map(|(k, v)| (k.to_string(), v.to_string()))
    .to_vec();
    assert_eq!(
        r.search_engine_settings_as_pairs("Comet").unwrap(),
        expected
    );
    assert_eq!(r.search_engine_settings_as_pairs("").unwrap(), expected);
    for engine in ["Percolator", "ConsensusID", "ConsensusIDPEP"] {
        r.search_engine = engine.into();
        assert_eq!(r.search_engine_settings_as_pairs("").unwrap(), expected);
        assert!(
            r.search_engine_settings_as_pairs(engine)
                .unwrap()
                .is_empty()
        );
    }
    for (value, text) in [
        (EnzymeTermSpecificity::Unknown, "unknown"),
        (EnzymeTermSpecificity::Semi, "semi"),
        (EnzymeTermSpecificity::None, "none"),
    ] {
        r.search_parameters.enzyme_specificity = value;
        assert_eq!(r.search_engine_settings_as_pairs("").unwrap()[7].1, text);
    }
}

#[test]
fn source_mztab_non_string_settings_regression_and_all_metadata_types() {
    // Literal regression from ProteinIdentification_test.cpp, with further
    // independent type/unit cases for the lenient DataValue conversion.
    let mut r = run("ConsensusID");
    let m = &mut r.search_parameters.metadata;
    m.insert("Comet:missed_cleavages".into(), 2i64.into());
    m.insert(
        "Comet:fragment_bin_tol".into(),
        MetaValue::new(MetaValueData::Float(0.02)).unwrap(),
    );
    m.insert("Comet:enzyme".into(), "Trypsin".into());
    m.insert("Comet:empty".into(), MetaValue::default());
    m.insert(
        "Comet:strings".into(),
        vec!["a,b".to_string(), "λ".into()].into(),
    );
    m.insert(
        "Comet:integers".into(),
        MetaValue::new(MetaValueData::IntegerList(vec![-2, 3])).unwrap(),
    );
    m.insert(
        "Comet:floats".into(),
        MetaValue::new(MetaValueData::FloatList(vec![-0.0, 1e-5]))
            .unwrap()
            .with_unit(Unit::new("UO:0000010", "second", "UO").unwrap())
            .unwrap(),
    );
    let values = r.search_engine_settings_as_pairs("Comet").unwrap();
    assert_eq!(
        values,
        [
            ("empty", ""),
            ("enzyme", "Trypsin"),
            ("floats", "[-0.0, 1.0e-05]"),
            ("fragment_bin_tol", "0.02"),
            ("integers", "[-2, 3]"),
            ("missed_cleavages", "2"),
            ("strings", "[a,b, λ]"),
        ]
        .map(|(k, v)| (k.into(), v.into()))
        .to_vec()
    );
}

#[test]
fn settings_literal_prefix_rule_clamps_and_checks_utf8_byte_boundaries() {
    let mut r = run("Merged");
    for key in ["Comet", "CometXfoo", "Comet:λ", "Other:Comet"] {
        r.search_parameters.metadata.insert(key.into(), "v".into());
    }
    assert_eq!(
        r.search_engine_settings_as_pairs("Comet").unwrap(),
        [
            ("".into(), "v".into()),
            ("λ".into(), "v".into()),
            ("foo".into(), "v".into()),
        ]
    );
    r.search_parameters
        .metadata
        .insert("Cometλ".into(), "v".into());
    assert!(r.search_engine_settings_as_pairs("Comet").is_err());
}

#[test]
fn singleton_groups_deduplicate_in_hit_order_and_preserve_existing_payload() {
    let mut r = run("Comet");
    r.indistinguishable_groups.push(ProteinGroup {
        probability: 0.75,
        accessions: vec!["B".into(), "A".into()],
        ..Default::default()
    });
    let existing = r.indistinguishable_groups[0].clone();
    r.hits = vec![
        hit("B", f64::NAN),
        hit("C", 7.0),
        hit("C", 99.0),
        hit("", -2.0),
        hit("D", 3.0),
    ];
    assert_eq!(
        r.fill_indistinguishable_groups_with_singletons().unwrap(),
        3
    );
    assert_eq!(r.indistinguishable_groups[0], existing);
    assert_eq!(
        r.indistinguishable_groups[1..]
            .iter()
            .map(|g| (g.accessions[0].as_str(), g.probability))
            .collect::<Vec<_>>(),
        [("C", 7.0), ("", -2.0), ("D", 3.0)]
    );
    assert_eq!(
        r.fill_indistinguishable_groups_with_singletons().unwrap(),
        0
    );
    r.hits.push(hit("bad", f64::NAN));
    let before = r.indistinguishable_groups.clone();
    assert!(r.fill_indistinguishable_groups_with_singletons().is_err());
    assert_eq!(r.indistinguishable_groups, before);
}

#[test]
fn mergeability_requires_engine_and_version_and_delegates_label_rules() {
    let a = run("Comet");
    let mut b = a.clone();
    assert!(a.peptide_ids_mergeable(&b, "label-free").unwrap());
    b.search_engine_version = "v2".into();
    assert!(!a.peptide_ids_mergeable(&b, "label-free").unwrap());
    b.search_engine_version = "v1".into();
    b.search_engine = "Other".into();
    assert!(!a.peptide_ids_mergeable(&b, "label-free").unwrap());
    b.search_engine = "Comet".into();
    b.search_parameters.variable_modifications.push("M".into());
    assert!(!a.peptide_ids_mergeable(&b, "label-free").unwrap());
    assert!(a.peptide_ids_mergeable(&b, "labeled_MS1").unwrap());
    b.search_parameters.charges = "invalid".into();
    assert!(a.peptide_ids_mergeable(&b, "labeled_MS1").is_err());
}

#[test]
fn metadata_copy_keeps_target_results_and_copies_paths_and_owned_payload() {
    let mut source = run("Fido");
    source.identifier = "source".into();
    source.date_time = Some("2026-01-02T03:04:05".into());
    source.score_type = "probability".into();
    source.higher_score_better = false;
    source.significance_threshold = -0.0;
    source.primary_ms_run_paths = vec!["run.mzML".into()];
    source.raw_ms_run_paths = vec!["run.raw".into()];
    source
        .metadata
        .insert("owned".into(), vec!["old".to_string()].into());
    source.set_inference_engine("Inference");
    source.search_parameters.fixed_modifications = vec!["mod".into()];
    source.hits = vec![hit("ignored", f64::NAN)];
    let mut target = run("Comet");
    target.hits.push(hit("retained", 3.0));
    target.protein_groups.push(ProteinGroup::default());
    target
        .indistinguishable_groups
        .push(ProteinGroup::default());
    let ptr = target.hits.as_ptr();
    let (groups, indist) = (
        target.protein_groups.clone(),
        target.indistinguishable_groups.clone(),
    );
    target.copy_metadata_only(&source).unwrap();
    assert_eq!(ptr, target.hits.as_ptr());
    assert_eq!(target.hits[0].accession, "retained");
    assert_eq!(target.protein_groups, groups);
    assert_eq!(target.indistinguishable_groups, indist);
    let mut expected = source.clone();
    expected.hits = target.hits.clone();
    expected.protein_groups = target.protein_groups.clone();
    expected.indistinguishable_groups = target.indistinguishable_groups.clone();
    assert_eq!(target, expected);
    assert_eq!(target.significance_threshold.to_bits(), (-0.0f64).to_bits());
    source.set_inference_engine("changed");
    source.primary_ms_run_paths[0].push('x');
    assert_eq!(target.inference_engine().unwrap(), "Inference");
    assert_eq!(target.primary_ms_run_paths, ["run.mzML"]);
}

#[test]
fn oversized_nested_metadata_fails_before_publication() {
    let mut source = run("Merged");
    source.search_parameters.metadata.insert(
        "Engine:list".into(),
        MetaValue::new(MetaValueData::IntegerList(vec![0; 1_000_001])).unwrap(),
    );
    assert!(source.search_engine_settings_as_pairs("Engine").is_err());
    let mut target = run("Keep");
    target.hits.push(hit("P", 9.0));
    let before = target.clone();
    assert!(target.copy_metadata_only(&source).is_err());
    assert_eq!(target, before);
    // Unused metadata is irrelevant to standard settings and group creation.
    assert_eq!(
        source.search_engine_settings_as_pairs("").unwrap().len(),
        12
    );
    assert_eq!(
        source
            .fill_indistinguishable_groups_with_singletons()
            .unwrap(),
        0
    );
}
