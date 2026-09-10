// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::analysis::id_conflict_resolver::*;
use openms::chemistry::AASequence;
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::kernel::{BaseFeature, ConsensusFeature, ConsensusMap, Feature, FeatureMap};

fn id(reference: &str, higher: bool, hits: &[(&str, f64, i32)]) -> PeptideIdentification {
    let mut result = PeptideIdentification {
        score_type: "score".into(),
        higher_score_better: higher,
        hits: hits
            .iter()
            .map(|(seq, score, charge)| {
                PeptideHit::new(*score, 42, *charge, AASequence::parse(seq).unwrap()).unwrap()
            })
            .collect(),
        ..Default::default()
    };
    result.set_spectrum_reference(reference);
    result
}
fn feature(intensity: f32, charge: i32, ids: Vec<PeptideIdentification>) -> Feature {
    Feature::from(BaseFeature {
        intensity,
        charge,
        peptide_identifications: ids,
        ..Default::default()
    })
}
fn reference(ids: &[PeptideIdentification]) -> Vec<String> {
    ids.iter()
        .map(PeptideIdentification::spectrum_reference)
        .collect()
}

#[test]
fn source_same_spectrum_direction_chimeras_charge_and_order() {
    for higher in [false, true] {
        let mut ids = vec![
            id("first", higher, &[("AAAK", 0.5, 2)]),
            id("scan=2075", higher, &[("PEPTIDEK", 0.9, 2)]),
            id("scan=2075", higher, &[("PEPTIDEK", 0.1, 2)]),
            id("scan=2075", higher, &[("PEPTIDER", 0.2, 2)]),
            id("scan=2075", higher, &[("PEPTIDEK", 0.3, 3)]),
            id("last", higher, &[("CCCK", 0.5, 2)]),
        ];
        let report = reduce_to_one_per_spectrum(&mut ids).unwrap();
        assert_eq!(report.removed, 1);
        assert_eq!(report.multiply_identified_spectra, 1);
        assert_eq!(ids[1].hits[0].score, if higher { 0.9 } else { 0.1 });
        assert_eq!(
            reference(&ids),
            ["first", "scan=2075", "scan=2075", "scan=2075", "last"]
        );
        assert_eq!(report.example, "scan=2075 / PEPTIDEK / charge 2");
    }
}

#[test]
fn source_stored_top_hit_empty_missing_and_inconsistent_groups() {
    let mut ids = vec![
        id("", false, &[("PEPTIDEK", 0.1, 2)]),
        id("", false, &[("PEPTIDEK", 0.9, 2)]),
        id("scan", false, &[]),
        id("scan", false, &[("PEPTIDEK", 0.1, 2)]),
        id("scan", true, &[("PEPTIDEK", 0.9, 2)]),
    ];
    let before = ids.clone();
    let report = reduce_to_one_per_spectrum(&mut ids).unwrap();
    assert_eq!(ids, before);
    assert_eq!(report.without_spectrum_reference, 2);
    assert_eq!(report.inconsistent_score_direction, 1);
    assert_eq!(report.multiply_identified_spectra, 0);
    let mut ids = vec![
        id("scan", false, &[("AAAK", 0.5, 2), ("PEPTIDEK", 0.01, 2)]),
        id("scan", false, &[("PEPTIDEK", 0.1, 2)]),
    ];
    assert_eq!(reduce_to_one_per_spectrum(&mut ids).unwrap().removed, 0);
    assert_eq!(ids[0].hits[0].sequence.as_str(), "AAAK"); // no implicit sorting
    assert_eq!(
        reduce_to_one_per_spectrum(&mut vec![]).unwrap(),
        UnresolvedIdentifications::default()
    );
}

#[test]
fn reduction_preserves_modified_forms_ties_and_rejects_mixed_runs_atomically() {
    let first = id("scan", true, &[("ACMK", 1., 2)]);
    let mut duplicate = first.clone();
    duplicate.metadata.insert("second".into(), "yes".into());
    let modified = id("scan", true, &[("ACM(Oxidation)K", 100., 2)]);
    let mut ids = vec![first.clone(), duplicate, modified.clone()];
    let report = reduce_to_one_per_spectrum(&mut ids).unwrap();
    assert_eq!(report.removed, 1);
    assert_eq!(ids, [first, modified]);
    ids[1].identifier = "another run".into();
    let before = ids.clone();
    assert!(reduce_to_one_per_spectrum(&mut ids).is_err());
    assert_eq!(ids, before);
    ids[1].identifier.clear();
    ids[1].hits[0].sequence = ids[0].hits[0].sequence.clone();
    ids[1].score_type = "other scale".into();
    let before = ids.clone();
    let report = reduce_to_one_per_spectrum(&mut ids).unwrap();
    assert_eq!(ids, before);
    assert_eq!(report.inconsistent_score_type, 1);
    assert_eq!(report.multiply_identified_spectra, 0);
}

#[test]
fn best_score_resolution_retains_metadata_and_source_unassigned_order() {
    for higher in [false, true] {
        let mut ids = vec![
            id("a", higher, &[("AAAK", 2., 2), ("AAAR", 5., 2)]),
            id("b", higher, &[("CCCK", 3., 2)]),
            id("c", higher, &[("DDDK", 1., 3)]),
        ];
        let mut removed = vec![id("existing", false, &[])];
        resolve_identifications(&mut ids, &mut removed, 100, ResolutionMethod::BestScore).unwrap();
        assert_eq!(reference(&ids), if higher { vec!["a"] } else { vec!["c"] });
        assert_eq!(
            reference(&removed),
            if higher {
                vec!["existing", "b", "c"]
            } else {
                vec!["existing", "a", "b"]
            }
        );
        assert_eq!(ids[0].hits.len(), 1);
        assert_eq!(ids[0].hits[0].rank, 42);
        assert_eq!(ids[0].metadata["feature_id"].as_str().unwrap(), "100");
        assert!(!removed[0].metadata.contains_key("feature_id"));
        for old in &removed[1..] {
            assert_eq!(old.hits.len(), 1);
            assert_eq!(old.metadata["feature_id"].as_str().unwrap(), "100");
        }
    }
}

#[test]
fn keep_matching_preserves_source_winner_alternatives_and_rejected_hits() {
    let mut ids = vec![
        id("matching", true, &[("CCCK", 2., 2), ("AAAK", 1., 3)]),
        id("rejected", true, &[("DDDK", 3., 2), ("EEEK", 2., 2)]),
        id("winner", true, &[("CCCK", 4., 2), ("AAAK", 10., 2)]),
    ];
    let mut removed = vec![];
    resolve_identifications(&mut ids, &mut removed, 17, ResolutionMethod::KeepMatching).unwrap();
    assert_eq!(reference(&ids), ["winner", "matching"]);
    assert_eq!(ids[0].hits.len(), 2); // actual source overload keeps winner's alternatives
    assert_eq!(ids[0].hits[0].sequence.as_str(), "AAAK");
    assert_eq!(ids[1].hits.len(), 1);
    assert_eq!(ids[1].hits[0].charge, 3); // matching ignores charge
    assert!(!ids[0].metadata.contains_key("feature_id"));
    assert!(!ids[1].metadata.contains_key("feature_id"));
    assert_eq!(removed[0].hits.len(), 2);
    assert_eq!(removed[0].metadata["feature_id"].as_str().unwrap(), "17");
}

#[test]
fn source_rank_aggregation_beats_best_single_record_and_keeps_original_score() {
    // Source SEQB uses unsupported ambiguous B; SEQC preserves the independent
    // rank/count calculation while both represented native peptides are valid.
    let ids = vec![
        id("one", true, &[("SEQC", 0.99, 2), ("SEQA", 0.5, 2)]),
        id("two", true, &[("SEQA", 0.8, 2), ("SEQC", 0.1, 2)]),
        id("three", true, &[("SEQA", 0.7, 2), ("SEQC", 0.05, 2)]),
    ];
    let mut best = ids.clone();
    resolve_identifications(&mut best, &mut vec![], 1, ResolutionMethod::BestScore).unwrap();
    assert_eq!(best[0].hits[0].sequence.as_str(), "SEQC");
    let mut map = ConsensusMap {
        features: vec![ConsensusFeature::from(BaseFeature {
            unique_id: 1,
            peptide_identifications: ids,
            ..Default::default()
        })],
        ..Default::default()
    };
    resolve_consensus_map(&mut map, ResolutionMethod::RankAggregation).unwrap();
    let result = &map.features[0].peptide_identifications;
    assert_eq!(reference(result), ["two"]);
    assert_eq!(result[0].hits.len(), 1);
    assert_eq!(result[0].hits[0].sequence.as_str(), "SEQA");
    assert_eq!(result[0].hits[0].score, 0.8); // aggregate 5/6 is selection-only
    assert_eq!(
        reference(&map.unassigned_peptide_identifications),
        ["one", "three"]
    );
    assert_eq!(
        map.unassigned_peptide_identifications[0].hits[0]
            .sequence
            .as_str(),
        "SEQC"
    );
}

#[test]
fn rank_missing_penalties_duplicates_and_sequence_length_tie_order() {
    let mut ids = vec![
        id(
            "a",
            false,
            &[("VVV", 1., 2), ("AAAK", 2., 2), ("AAAK", 3., 3)],
        ),
        id("b", false, &[("AAAK", 1., 2)]),
        id("empty", false, &[]),
    ];
    let mut removed = vec![];
    resolve_identifications(&mut ids, &mut removed, 9, ResolutionMethod::RankAggregation).unwrap();
    // max_hits3,N3; AAAK ranks1+0+3=4; VVV ranks0+3+3=6.
    assert_eq!(ids[0].hits[0].sequence.as_str(), "AAAK");
    assert_eq!(ids[0].hits[0].score, 1.);
    assert_eq!(removed[1].hits.len(), 0);
    let mut ids = vec![
        id("long", true, &[("AAAA", 1., 2)]),
        id("short", true, &[("VVV", 1., 2)]),
    ];
    resolve_identifications(&mut ids, &mut vec![], 0, ResolutionMethod::RankAggregation).unwrap();
    assert_eq!(reference(&ids), ["short"]); // source AASequence compares length before letters
}

#[test]
fn empty_identifications_are_safe_in_all_methods() {
    for method in [
        ResolutionMethod::BestScore,
        ResolutionMethod::KeepMatching,
        ResolutionMethod::RankAggregation,
    ] {
        let mut ids = vec![
            id("empty", false, &[]),
            id("populated", false, &[("AAAK", 0.1, 2)]),
        ];
        let mut removed = vec![];
        resolve_identifications(&mut ids, &mut removed, 1, method).unwrap();
        assert_eq!(reference(&ids), ["populated"]);
        assert_eq!(reference(&removed), ["empty"]);
        let mut ids = vec![id("first", false, &[]), id("second", true, &[])];
        resolve_identifications(&mut ids, &mut removed, 1, method).unwrap();
        assert_eq!(reference(&ids), ["first"]);
        assert!(ids[0].hits.is_empty());
    }
}

#[test]
fn source_between_feature_intensity_charge_and_modification_case_is_fully_exercised() {
    let original = id("unmodified", true, &[("MORRISSEY", 23., 0)]);
    let modified = id("modified", true, &[("M(Oxidation)ORRISSEY", 23., 0)]);
    let mut map = FeatureMap {
        features: vec![
            feature(1000., 2, vec![original.clone()]),
            feature(10000., 2, vec![original.clone()]),
            feature(1000., 3, vec![original.clone()]),
            feature(1001., 2, vec![modified.clone()]),
            feature(10000., 2, vec![original]),
        ],
        ..Default::default()
    };
    resolve_between_features(&mut map).unwrap();
    assert_eq!(map.features.len(), 5);
    assert!(map.features[0].peptide_identifications.is_empty());
    assert_eq!(map.features[1].peptide_identifications.len(), 1);
    assert_eq!(map.features[2].peptide_identifications.len(), 1);
    assert_eq!(map.features[3].peptide_identifications, [modified]);
    assert!(map.features[4].peptide_identifications.is_empty()); // tie retains first
    assert_eq!(map.unassigned_peptide_identifications.len(), 2);
}

#[test]
fn map_annotations_subordinates_and_late_failures_are_transactional() {
    let mut f = feature(
        100.,
        2,
        vec![
            id("a", false, &[("AAAK", 0.1, 2)]),
            id("b", false, &[("CCCK", 0.2, 2)]),
        ],
    );
    f.base.unique_id = 99;
    f.subordinates
        .push(feature(1., 2, vec![id("child", true, &[("DDDK", 1., 2)])]));
    let child = f.subordinates.clone();
    let mut map = FeatureMap {
        features: vec![f],
        unassigned_peptide_identifications: vec![id("existing", true, &[])],
        ..Default::default()
    };
    resolve_feature_map(&mut map, ResolutionMethod::BestScore).unwrap();
    assert_eq!(map.features[0].subordinates, child);
    assert_eq!(map.features[0].metadata["feature_id"], "99");
    assert_eq!(
        map.unassigned_peptide_identifications[0].metadata["feature_id"]
            .as_str()
            .unwrap(),
        "not mapped"
    );
    assert_eq!(
        map.unassigned_peptide_identifications[1].metadata["feature_id"]
            .as_str()
            .unwrap(),
        "99"
    );
    let mut invalid = feature(
        1.,
        2,
        vec![
            id("late", true, &[("AAAK", 1., 2)]),
            id("bad", false, &[("AAAK", 1., 2)]),
        ],
    );
    invalid.base.unique_id = 100;
    map.features.push(invalid);
    let before = map.clone();
    assert!(resolve_feature_map(&mut map, ResolutionMethod::BestScore).is_err());
    assert_eq!(map, before);
    assert!(resolve_between_features(&mut map).is_err());
    assert_eq!(map, before);
}

#[test]
fn consensus_between_resolution_sorts_and_moves_complete_identifications() {
    let mut map = ConsensusMap {
        features: [10., 20.]
            .into_iter()
            .map(|intensity| {
                ConsensusFeature::from(BaseFeature {
                    intensity,
                    charge: 2,
                    peptide_identifications: vec![id(
                        "scan",
                        true,
                        &[("VVVK", 0.1, 2), ("AAAK", 2., 2)],
                    )],
                    ..Default::default()
                })
            })
            .collect(),
        ..Default::default()
    };
    resolve_between_consensus_features(&mut map).unwrap();
    assert!(map.features[0].peptide_identifications.is_empty());
    assert_eq!(map.features[1].peptide_identifications[0].hits.len(), 2);
    assert_eq!(
        map.unassigned_peptide_identifications[0].hits[0]
            .sequence
            .as_str(),
        "AAAK"
    );
    assert_eq!(map.unassigned_peptide_identifications[0].hits.len(), 2);
}
