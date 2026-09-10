// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::analysis::scores::*;
use openms::chemistry::AASequence;
use openms::identification::{
    PeptideHit, PeptideIdentification, ProteinHit, ProteinIdentification,
};
use openms::kernel::{BaseFeature, ConsensusFeature, ConsensusMap, Feature, FeatureMap};
use openms::metadata::MetaValue;

fn id(score: f64, pep: f64) -> PeptideIdentification {
    let mut hit = PeptideHit::new(score, 3, 2, AASequence::parse("PEPTIDE").unwrap()).unwrap();
    hit.metadata
        .insert("pep".into(), MetaValue::try_from(pep).unwrap());
    PeptideIdentification {
        hits: vec![hit],
        score_type: "XTandem".into(),
        ..Default::default()
    }
}

#[test]
fn source_score_registry_categories_suffixes_and_case() {
    let names: Vec<_> = all_score_names().collect();
    assert_eq!(names.len(), 29);
    for kind in ScoreType::ALL {
        assert!(kind.names().windows(2).all(|pair| pair[0] < pair[1]));
        for name in kind.names() {
            assert_eq!(ScoreType::from_name(name), Some(kind));
            assert!(kind.matches(&format!("{name}_score")));
            assert_eq!(ScoreType::from_name(&format!("{name}_score")), None);
        }
    }
    assert!(ScoreType::Raw.matches("XTandem"));
    assert!(!ScoreType::Raw.matches("xtandem")); // Implementation is case-sensitive despite C++ header text.
    assert!(ScoreType::Raw.higher_is_better());
    assert!(ScoreType::PosteriorProbability.higher_is_better());
    assert!(!ScoreType::PosteriorErrorProbability.higher_is_better());
    for (text, expected) in [
        ("RAW", ScoreType::Raw),
        ("Raw E-Value", ScoreType::RawEValue),
        ("q-value_score", ScoreType::QValue),
        (
            "Posterior_Error_Probability",
            ScoreType::PosteriorErrorProbability,
        ),
        ("False discovery rate", ScoreType::Fdr),
        ("PP", ScoreType::PosteriorProbability),
    ] {
        assert_eq!(ScoreType::parse(text).unwrap(), expected);
    }
    assert!(ScoreType::parse("hyperscore").is_err()); // Score name and category are distinct APIs.
}

#[test]
fn score_detection_uses_source_main_and_first_hit_priority() {
    let mut value = id(10.0, 0.1);
    assert_eq!(
        find_peptide_score(&value, ScoreType::Raw).unwrap(),
        ScoreSearchResult {
            is_main_score: true,
            name: "XTandem".into()
        }
    );
    assert_eq!(
        find_peptide_score(&value, ScoreType::PosteriorErrorProbability)
            .unwrap()
            .name,
        "pep"
    );
    value.hits[0].metadata.insert("PEP".into(), 0_i64.into());
    assert_eq!(
        find_peptide_score(&value, ScoreType::PosteriorErrorProbability)
            .unwrap()
            .name,
        "PEP"
    );
    value.hits[0]
        .metadata
        .insert("MS:1001493_score".into(), 1_i64.into());
    assert_eq!(
        find_peptide_score(&value, ScoreType::PosteriorErrorProbability)
            .unwrap()
            .name,
        "MS:1001493_score"
    );
    assert!(find_peptide_score(&value, ScoreType::QValue).is_none());
    value.hits.clear();
    assert!(find_peptide_score(&value, ScoreType::PosteriorErrorProbability).is_none());
}

#[test]
fn switch_and_restore_preserve_score_metadata_hit_order_and_ranks() {
    let mut ids = vec![id(10.0, 0.1), id(20.0, 0.02)];
    let switcher = ScoreSwitcher::new("pep", false);
    assert_eq!(switcher.switch_peptides(&mut ids).unwrap(), 2);
    assert_eq!(ids[0].hits[0].score, 0.1);
    assert_eq!(ids[1].hits[0].score, 0.02);
    assert_eq!(ids[0].hits[0].metadata["XTandem"].as_f64().unwrap(), 10.0);
    assert_eq!(ids[0].hits[0].rank, 3);
    assert_eq!(ids[0].score_type, "pep");
    assert!(!ids[0].higher_score_better);
    ScoreSwitcher::new("XTandem", true)
        .switch_peptides(&mut ids)
        .unwrap();
    assert_eq!(ids[0].hits[0].score, 10.0);
    assert_eq!(ids[0].hits[0].metadata["pep"].as_f64().unwrap(), 0.1);
}

#[test]
fn backup_collisions_keep_original_scores_and_fail_atomically_when_occupied() {
    let mut ids = vec![id(10.0, 0.1)];
    ids[0].hits[0]
        .metadata
        .insert("XTandem".into(), 20_i64.into());
    ScoreSwitcher::new("pep", false)
        .switch_peptides(&mut ids)
        .unwrap();
    assert_eq!(ids[0].hits[0].metadata["XTandem"].as_i64().unwrap(), 20);
    assert_eq!(ids[0].hits[0].metadata["XTandem~"].as_f64().unwrap(), 10.0);
    let mut ids = vec![id(1.0, 0.1), id(10.0, 0.2)];
    ids[1].hits[0]
        .metadata
        .insert("XTandem".into(), 20_i64.into());
    ids[1].hits[0]
        .metadata
        .insert("XTandem~".into(), 30_i64.into());
    let before = ids.clone();
    assert!(
        ScoreSwitcher::new("pep", false)
            .switch_peptides(&mut ids)
            .is_err()
    );
    assert_eq!(ids, before);
    // A compatible occupied backup retains its exact numeric type and unit.
    let backup = MetaValue::from(10_i64)
        .with_unit(openms::metadata::Unit::new("UO:0000186", "dimensionless unit", "UO").unwrap())
        .unwrap();
    ids[1].hits[0]
        .metadata
        .insert("XTandem~".into(), backup.clone());
    ScoreSwitcher::new("pep", false)
        .switch_peptides(&mut ids)
        .unwrap();
    assert_eq!(ids[1].hits[0].metadata["XTandem~"], backup);
}

#[test]
fn score_backup_relative_tolerance_handles_zero_signs_and_large_finite_values() {
    for (main, existing, backup) in [
        (0.0, 0.0, false),
        (10.0, 10.000001, false),
        (10.0, 10.1, true),
        (-1.0, 1.0, true),
        (f64::MAX, f64::MAX / 2.0, true),
    ] {
        let mut ids = vec![id(main, 0.1)];
        ids[0].hits[0]
            .metadata
            .insert("XTandem".into(), MetaValue::try_from(existing).unwrap());
        ScoreSwitcher::new("pep", false)
            .switch_peptides(&mut ids)
            .unwrap();
        assert_eq!(ids[0].hits[0].metadata.contains_key("XTandem~"), backup);
    }
}

#[test]
fn malformed_or_missing_late_scores_never_partially_change_records() {
    let mut ids = vec![id(10.0, 0.1), id(20.0, 0.2)];
    ids[1].hits[0].metadata.insert("pep".into(), "0.2".into());
    let before = ids.clone();
    assert!(
        ScoreSwitcher::new("pep", false)
            .switch_peptides(&mut ids)
            .is_err()
    );
    assert_eq!(ids, before); // Typed numeric conversion does not silently parse strings.
    ids[1].hits[0].metadata.remove("pep");
    let before = ids.clone();
    assert!(
        ScoreSwitcher::new("pep", false)
            .switch_peptides(&mut ids)
            .is_err()
    );
    assert_eq!(ids, before);
    assert!(
        ScoreSwitcher::new("", false)
            .switch_peptides(&mut [])
            .is_err()
    );
}

#[test]
fn score_backups_cannot_create_invalid_target_decoy_annotations() {
    let mut ids = vec![id(10.0, 0.1), id(20.0, 0.2)];
    ids[1].score_type = "target_decoy".into();
    let before = ids.clone();
    assert!(
        ScoreSwitcher::new("pep", false)
            .switch_peptides(&mut ids)
            .is_err()
    );
    assert_eq!(ids, before);

    let mut hit = ProteinHit::new(9.0, 0, "P1", "PEPTIDE").unwrap();
    hit.metadata
        .insert("pep".into(), MetaValue::try_from(0.1).unwrap());
    let mut runs = vec![ProteinIdentification {
        hits: vec![hit],
        score_type: "Mascot".into(),
        ..Default::default()
    }];
    let before = runs.clone();
    let switcher = ScoreSwitcher {
        old_score: Some("target_decoy".into()),
        ..ScoreSwitcher::new("pep", false)
    };
    assert!(switcher.switch_proteins(&mut runs).is_err());
    assert_eq!(runs, before);
    runs[0].validate().unwrap();
}

#[test]
fn category_switch_checks_each_record_even_if_first_already_matches() {
    let mut first = id(0.01, 0.1);
    first.score_type = "PEP".into();
    first.higher_score_better = false;
    let mut second = id(10.0, 0.2);
    second.hits[0]
        .metadata
        .insert("PEP_score".into(), MetaValue::try_from(0.05).unwrap());
    let mut ids = vec![first.clone(), second];
    assert_eq!(
        switch_peptides_to_category(&mut ids, ScoreType::PosteriorErrorProbability).unwrap(),
        1
    );
    assert_eq!(ids[0], first);
    assert_eq!(ids[1].score_type, "PEP");
    assert_eq!(ids[1].hits[0].score, 0.05);
    assert!(!ids[1].higher_score_better);
}

#[test]
fn protein_and_map_adapters_preserve_nonselected_records_and_atomicity() {
    let mut protein = ProteinHit::new(9.0, 0, "P1", "PEPTIDE").unwrap();
    protein
        .metadata
        .insert("FDR".into(), MetaValue::try_from(0.02).unwrap());
    let mut runs = vec![ProteinIdentification {
        hits: vec![protein],
        score_type: "Mascot".into(),
        ..Default::default()
    }];
    assert_eq!(
        switch_proteins_to_category(&mut runs, ScoreType::Fdr).unwrap(),
        1
    );
    assert_eq!(runs[0].hits[0].score, 0.02);
    let feature = Feature {
        base: BaseFeature {
            peptide_identifications: vec![id(9.0, 0.03)],
            ..Default::default()
        },
        subordinates: vec![Feature {
            base: BaseFeature {
                peptide_identifications: vec![id(8.0, 0.04)],
                ..Default::default()
            },
            ..Default::default()
        }],
        ..Default::default()
    };
    let mut map = FeatureMap {
        features: vec![feature],
        unassigned_peptide_identifications: vec![id(7.0, 0.05)],
        ..Default::default()
    };
    assert_eq!(
        ScoreSwitcher::new("pep", false)
            .switch_feature_map(&mut map, false)
            .unwrap(),
        1
    );
    assert_eq!(
        map.features[0].peptide_identifications[0].hits[0].score,
        0.03
    );
    assert_eq!(
        map.features[0].subordinates[0].peptide_identifications[0].hits[0].score,
        8.0
    );
    assert_eq!(map.unassigned_peptide_identifications[0].hits[0].score, 7.0);
    let mut feature = ConsensusFeature::default();
    feature.peptide_identifications = vec![id(9.0, 0.03)];
    let mut map = ConsensusMap {
        features: vec![feature],
        unassigned_peptide_identifications: vec![id(8.0, 0.04)],
        ..Default::default()
    };
    map.unassigned_peptide_identifications[0].hits[0]
        .metadata
        .remove("pep");
    let before = map.clone();
    assert!(
        ScoreSwitcher::new("pep", false)
            .switch_consensus_map(&mut map, true)
            .is_err()
    );
    assert_eq!(map, before);
}
