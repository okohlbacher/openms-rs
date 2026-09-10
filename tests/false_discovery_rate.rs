// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::analysis::false_discovery_rate::*;
use openms::chemistry::AASequence;
use openms::identification::{
    PeptideHit, PeptideIdentification, ProteinGroup, ProteinHit, ProteinIdentification,
    TargetDecoyType,
};
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t, "{a} != {b} (tolerance {t})");
}
fn label(value: &str) -> TargetDecoyType {
    match value {
        "target" => TargetDecoyType::Target,
        "decoy" => TargetDecoyType::Decoy,
        "target+decoy" => TargetDecoyType::TargetAndDecoy,
        _ => panic!("bad fixture label"),
    }
}
fn peptide(score: f64, td: &str, sequence: &str, higher: bool) -> PeptideIdentification {
    let mut hit = PeptideHit::new(score, 0, 2, AASequence::parse(sequence).unwrap()).unwrap();
    hit.set_target_decoy_type(label(td));
    PeptideIdentification {
        score_type: "score".into(),
        higher_score_better: higher,
        hits: vec![hit],
        ..Default::default()
    }
}
fn protein(score: f64, td: &str, accession: &str) -> ProteinHit {
    let mut hit = ProteinHit::new(score, 0, accession, "").unwrap();
    hit.set_target_decoy_type(label(td)).unwrap();
    hit
}
fn table(text: &str) -> Vec<Vec<&str>> {
    text.lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').collect())
        .collect()
}

#[test]
fn upstream_omssa_concatenated_psm_qvalues() {
    let mut ids: Vec<_> = table(include_str!("data/fdr_omssa.tsv"))
        .iter()
        .map(|r| {
            let mut id = peptide(r[1].parse().unwrap(), r[3], "PEPTIDE", false);
            id.score_type = "OMSSA".into();
            id.hits[0].charge = r[2].parse().unwrap();
            id
        })
        .collect();
    let source_count = ids.len();
    FalseDiscoveryRate::default()
        .apply_peptides(&mut ids, false)
        .unwrap();
    assert_eq!(source_count, 1534);
    close(ids[0].hits[0].score, 0.0730478589420655, 1e-14);
    close(ids[5].hits[0].score, 0.409926470588235, 1e-14);
    assert!(ids[9].hits.is_empty());
    let (mut small, mut large) = (0, 0);
    for id in &ids {
        assert_eq!(id.score_type, "q-value");
        assert!(!id.higher_score_better);
        for hit in &id.hits {
            let old = hit.metadata["OMSSA_score"].as_f64().unwrap();
            if old <= 1e-3 {
                close(hit.score, 0., 1e-15);
                small += 1;
            }
            if old >= 1000. && hit.charge != 1 {
                assert!(hit.score > 0.1);
                large += 1;
            }
        }
    }
    assert!(small > 0 && large > 0);
}

#[test]
fn upstream_xtandem_separate_psm_and_concatenated_protein_goldens() {
    let mut forward: Vec<_> = table(include_str!("data/fdr_xtandem_fwd_peptides.tsv"))
        .iter()
        .map(|r| peptide(r[1].parse().unwrap(), "target", "PEPTIDE", true))
        .collect();
    let mut reverse: Vec<_> = table(include_str!("data/fdr_xtandem_rev_peptides.tsv"))
        .iter()
        .map(|r| peptide(r[1].parse().unwrap(), "decoy", "PEPTIDE", true))
        .collect();
    let before = reverse.clone();
    FalseDiscoveryRate::default()
        .apply_separate_peptides(&mut forward, &mut reverse)
        .unwrap();
    assert_eq!(reverse, before);
    let (mut high, mut boundary) = (0, 0);
    for id in &forward {
        let hit = &id.hits[0];
        let old = hit.metadata["score_score"].as_f64().unwrap();
        if old >= 39.4 {
            close(hit.score, 0., 1e-15);
            high += 1;
        }
        if (old - 37.9).abs() < 1e-4 {
            close(hit.score, 0.08, 1e-14);
            boundary += 1;
        }
    }
    assert!(high > 0 && boundary > 0);
    let mut ids = Vec::new();
    for (fixture, td) in [
        (include_str!("data/fdr_xtandem_fwd_proteins.tsv"), "target"),
        (include_str!("data/fdr_xtandem_rev_proteins.tsv"), "decoy"),
    ] {
        ids.push(ProteinIdentification {
            score_type: "XTandem".into(),
            higher_score_better: false,
            hits: table(fixture)
                .iter()
                .map(|r| protein(r[1].parse().unwrap(), td, r[0]))
                .collect(),
            ..Default::default()
        });
    }
    FalseDiscoveryRate::default()
        .apply_proteins(&mut ids)
        .unwrap();
    assert!(ids[1].hits.is_empty());
    let (mut high, mut boundary, mut low) = (0, 0, 0);
    for hit in &ids[0].hits {
        let old = hit.metadata["XTandem_score"].as_f64().unwrap();
        if old < -1.8 {
            close(hit.score, 0., 1e-15);
            high += 1;
        }
        if old == 0. {
            close(hit.score, 0.897384, 1e-6);
            boundary += 1;
        }
        if old > -1.2 {
            assert!(hit.score > 0.1);
            low += 1;
        }
    }
    assert!(high > 0 && boundary > 0 && low > 0);
}

#[test]
fn upstream_picked_protein_literal_values_and_decoy_affix_validation() {
    let mut id = ProteinIdentification {
        score_type: "OMSSA".into(),
        higher_score_better: false,
        hits: table(include_str!("data/fdr_picked.tsv"))
            .iter()
            .map(|r| protein(r[1].parse().unwrap(), r[2], r[0]))
            .collect(),
        ..Default::default()
    };
    FalseDiscoveryRate::default()
        .apply_picked_protein(&mut id, DecoyAffix::Prefix("decoy_"), true)
        .unwrap();
    assert_eq!(id.hits.len(), 6);
    for (hit, expected) in id.hits.iter().zip([0.25, 0.25, 0.25, 0.4, 0.4, 0.5]) {
        close(hit.score, expected, 1e-14);
    }
    let mut bad = ProteinIdentification {
        hits: vec![protein(2., "target", "A"), protein(1., "decoy", "A_decoy")],
        ..Default::default()
    };
    let before = bad.clone();
    assert!(
        FalseDiscoveryRate::default()
            .apply_picked_protein(&mut bad, DecoyAffix::Prefix("decoy_"), false)
            .is_err()
    );
    assert_eq!(bad, before);
    FalseDiscoveryRate::default()
        .apply_picked_protein(&mut bad, DecoyAffix::Suffix("_decoy"), false)
        .unwrap();
    assert_eq!(bad.hits.len(), 1);
}

#[test]
fn distinct_legacy_and_basic_formulas_ties_directions_and_fractional_labels() {
    let targets = [10., 8., 6.];
    let decoys = [9., 7.];
    let q = FalseDiscoveryRate::default()
        .calculate_legacy(&targets, &decoys, true)
        .unwrap();
    for (score, value) in [
        (10., 0.),
        (9., 0.5),
        (8., 0.5),
        (7., 2. / 3.),
        (6., 2. / 3.),
    ] {
        close(q.value(score).unwrap(), value, 1e-14);
    }
    let raw = FalseDiscoveryRate {
        output: FdrOutput::Fdr,
        ..Default::default()
    };
    let raw_curve = raw.calculate_legacy(&targets, &decoys, true).unwrap();
    // Legacy raw decoys inherit the first target's value under its source ordering.
    assert_eq!(raw_curve.value(9.).unwrap(), 0.);
    close(raw_curve.value(8.).unwrap(), 0.5, 1e-14);
    // A raw decoy/target tie overwrites the shared source lookup entry as in C++.
    let tied = raw.calculate_legacy(&targets, &[8.], true).unwrap();
    assert_eq!(tied.value(8.).unwrap(), 0.);
    let observations = [
        ScoreLabel::new(10., true),
        ScoreLabel::new(9., false),
        ScoreLabel::new(8., true),
        ScoreLabel::new(8., true),
        ScoreLabel::new(7., false),
        ScoreLabel::new(6., true),
    ];
    let basic = raw.calculate_basic(&observations, true).unwrap();
    for (score, value) in [(10., 0.5), (9., 1.), (8., 0.5), (7., 0.75), (6., 0.6)] {
        close(basic.value(score).unwrap(), value, 1e-14);
    }
    let q = FalseDiscoveryRate::default()
        .calculate_basic(&observations, true)
        .unwrap();
    for (score, value) in [(10., 0.5), (9., 0.5), (8., 0.5), (7., 0.6), (6., 0.6)] {
        close(q.value(score).unwrap(), value, 1e-14);
    }
    let mirrored: Vec<_> = observations
        .iter()
        .map(|p| ScoreLabel {
            score: -p.score,
            ..*p
        })
        .collect();
    let low = FalseDiscoveryRate::default()
        .calculate_basic(&mirrored, false)
        .unwrap();
    for p in &observations {
        close(
            low.value(-p.score).unwrap(),
            q.value(p.score).unwrap(),
            1e-14,
        );
    }
    let alternate = FalseDiscoveryRate {
        conservative: false,
        ..raw
    }
    .calculate_basic(&observations, true)
    .unwrap();
    close(alternate.value(8.).unwrap(), 2. / 5., 1e-14);
    let fractions = [
        ScoreLabel {
            score: 10.,
            target_fraction: 0.5,
        },
        ScoreLabel::new(9., true),
    ];
    let curve = FalseDiscoveryRate::default()
        .calculate_basic(&fractions, true)
        .unwrap();
    assert_eq!(curve.value(10.).unwrap(), 0.6);
    let all_decoys = raw
        .calculate_basic(
            &[ScoreLabel::new(1., false), ScoreLabel::new(2., false)],
            true,
        )
        .unwrap();
    assert_eq!(all_decoys.value(1.).unwrap(), 3.);
    let near = raw
        .calculate_basic(
            &[
                ScoreLabel::new(1., true),
                ScoreLabel::new(1. - 5e-13, false),
            ],
            true,
        )
        .unwrap();
    assert_eq!(near.len(), 1);
    assert_eq!(near.value(1. - 5e-13).unwrap(), 1.);
}

#[test]
fn upstream_peptide_level_best_representative_is_order_independent_and_prefers_target_ties() {
    fn build(reverse: bool, first_label: &str, first_score: f64) -> Vec<PeptideIdentification> {
        let mut ids = vec![
            peptide(0.9, "target", "AAAAAAAK", true),
            peptide(0.8, "target", "CCCCCCCK", true),
            peptide(0.7, "target", "DDDDDDDK", true),
            peptide(0.5, "decoy", "FFFFFFFK", true),
            peptide(0.3, "decoy", "HHHHHHHK", true),
        ];
        let mut duplicates = vec![
            peptide(first_score, first_label, "EEEEEEEK", true),
            peptide(0.95, "target", "EEEEEEEK", true),
        ];
        if reverse {
            duplicates.reverse();
        }
        ids.extend(duplicates);
        ids
    }
    for (first_label, first_score) in [("target", 0.1), ("decoy", 0.1), ("decoy", 0.95)] {
        let mut a = build(false, first_label, first_score);
        let mut b = build(true, first_label, first_score);
        FalseDiscoveryRate::default()
            .apply_basic_peptide_level(&mut a)
            .unwrap();
        FalseDiscoveryRate::default()
            .apply_basic_peptide_level(&mut b)
            .unwrap();
        let score = |ids: &[PeptideIdentification]| {
            ids.iter()
                .filter_map(|id| id.hits.first())
                .find(|h| h.sequence.as_str() == "EEEEEEEK")
                .unwrap()
                .score
        };
        close(score(&a), score(&b), 1e-14);
        close(score(&a), 0.2, 1e-14);
        for id in a.iter().filter(|id| !id.hits.is_empty()) {
            assert_eq!(id.score_type, "peptide q-value");
            assert!(!id.higher_score_better);
            assert!(id.hits[0].metadata.contains_key("score"));
        }
    }
}

#[test]
fn legacy_annotations_best_hit_selection_run_charge_pools_and_metadata() {
    let mut ids = vec![
        peptide(10., "target", "M(Oxidation)K", true),
        peptide(8., "target", "MK", true),
        peptide(9., "decoy", "KK", true),
    ];
    let decoy = ids[2].hits[0].clone();
    ids[0].hits.push(decoy);
    ids[0].hits.reverse();
    ids[0].rt = Some(4.);
    ids[0].metadata.insert("sample".into(), "x".into());
    FalseDiscoveryRate {
        add_decoy_peptides: true,
        ..Default::default()
    }
    .apply_peptides(&mut ids, true)
    .unwrap();
    assert_eq!(ids[0].hits.len(), 1);
    assert_eq!(ids[0].rt, Some(4.));
    assert_eq!(ids[0].metadata["sample"].as_str().unwrap(), "x");
    assert_eq!(
        ids[0].hits[0].metadata["peptide q-value"],
        ids[1].hits[0].metadata["peptide q-value"]
    );
    let mut pools = vec![
        peptide(1., "target", "PEPTIDE", true),
        peptide(2., "decoy", "PEPTIDE", true),
        peptide(1., "target+decoy", "PEPTIDE", false),
    ];
    pools[0].identifier = "a".into();
    pools[1].identifier = "a".into();
    pools[1].hits[0].charge = 3;
    pools[2].identifier = "b".into();
    let options = FalseDiscoveryRate {
        split_charge_variants: true,
        treat_runs_separately: true,
        add_decoy_peptides: true,
        ..Default::default()
    };
    options.apply_peptides(&mut pools, false).unwrap();
    assert_eq!(pools[0].hits[0].score, 0.);
    assert!(pools[1].hits.is_empty());
    assert_eq!(pools[2].hits[0].score, 0.);
    // The legacy all-decoy special case removes decoys even when retention is set.
}

#[test]
fn basic_psm_first_hit_counting_keeps_secondary_hits_and_checks_sorting() {
    let mut ids = vec![
        peptide(10., "target", "AK", true),
        peptide(9., "decoy", "CK", true),
        peptide(8., "target", "DK", true),
    ];
    ids[0]
        .hits
        .push(peptide(7., "target", "EK", true).hits.remove(0));
    FalseDiscoveryRate::default()
        .apply_basic_peptides(&mut ids)
        .unwrap();
    assert_eq!(ids[0].hits.len(), 2);
    assert!(ids[1].hits.is_empty());
    assert_eq!(
        ids[0].hits[0].metadata["score_score"].as_f64().unwrap(),
        10.
    );
    assert_eq!(ids[0].hits[1].metadata["score_score"].as_f64().unwrap(), 7.);
    let mut unsorted = vec![peptide(1., "target", "AK", true)];
    unsorted[0]
        .hits
        .push(peptide(2., "target", "CK", true).hits.remove(0));
    let before = unsorted.clone();
    assert!(
        FalseDiscoveryRate::default()
            .apply_basic_peptides(&mut unsorted)
            .is_err()
    );
    assert_eq!(unsorted, before);
}

#[test]
fn basic_protein_groups_labels_and_non_fdr_fields_are_preserved() {
    let mut id = ProteinIdentification {
        score_type: "score".into(),
        hits: vec![
            protein(10., "target", "A"),
            protein(9., "decoy", "B"),
            protein(8., "target", "C"),
        ],
        indistinguishable_groups: vec![
            ProteinGroup {
                probability: 10.,
                accessions: vec!["A".into(), "B".into()],
                ..Default::default()
            },
            ProteinGroup {
                probability: 9.,
                accessions: vec!["B".into()],
                ..Default::default()
            },
            ProteinGroup {
                probability: 8.,
                accessions: vec!["C".into()],
                ..Default::default()
            },
        ],
        ..Default::default()
    };
    id.protein_groups.push(ProteinGroup {
        probability: 42.,
        ..Default::default()
    });
    id.metadata.insert("dataset".into(), "kept".into());
    id.primary_ms_run_paths.push("run.mzML".into());
    FalseDiscoveryRate::default()
        .apply_basic_protein(&mut id, true)
        .unwrap();
    assert_eq!(id.hits.len(), 2);
    assert_eq!(id.indistinguishable_groups.len(), 3);
    assert_eq!(id.protein_groups[0].probability, 42.);
    assert_eq!(id.metadata["dataset"].as_str().unwrap(), "kept");
    assert_eq!(id.primary_ms_run_paths, ["run.mzML"]);
    for (group, value) in id
        .indistinguishable_groups
        .iter()
        .zip([0.5, 2. / 3., 2. / 3.])
    {
        close(group.probability, value, 1e-14);
    }
}

#[test]
fn estimated_probability_scores_ties_and_roc_cutoff_groups() {
    let estimator = FalseDiscoveryRate {
        add_decoy_proteins: true,
        ..Default::default()
    };
    let pep = [0.01, 0.02, 0.02, 0.2];
    let pp = pep.map(|v| 1. - v);
    let a = estimator.calculate_estimated(&pep, false).unwrap();
    let b = estimator.calculate_estimated(&pp, true).unwrap();
    for &p in &pep {
        close(a.value(p).unwrap(), b.value(1. - p).unwrap(), 1e-14);
    }
    close(a.value(0.02).unwrap(), 0.05 / 3., 1e-14);
    let mut id = ProteinIdentification {
        score_type: "Posterior Error Probability".into(),
        higher_score_better: false,
        hits: pep
            .iter()
            .enumerate()
            .map(|(i, &s)| ProteinHit::new(s, 0, format!("P{i}"), "").unwrap())
            .collect(),
        ..Default::default()
    };
    estimator.apply_estimated_protein(&mut id).unwrap();
    assert_eq!(id.score_type, "Estimated Q-Values");
    assert!(!id.higher_score_better);
    close(id.hits[3].score, 0.0625, 1e-14);
    let observations = [
        ScoreLabel::new(3., true),
        ScoreLabel::new(2., true),
        ScoreLabel::new(2., false),
        ScoreLabel::new(1., false),
    ];
    close(
        estimator.roc_n(&observations, true, None).unwrap(),
        0.875,
        1e-14,
    );
    close(
        estimator.roc_n(&observations, true, Some(1)).unwrap(),
        0.75,
        1e-14,
    );
    assert_eq!(
        estimator
            .roc_n(&[ScoreLabel::new(1., true)], true, None)
            .unwrap(),
        1.
    );
    assert_eq!(estimator.roc_n(&[], true, None).unwrap(), 0.);
}

#[test]
fn peptide_level_score_backup_preserves_reserved_labels_atomically() {
    for output in [FdrOutput::QValue, FdrOutput::Fdr] {
        for keep_decoys in [false, true] {
            let options = FalseDiscoveryRate {
                output,
                add_decoy_peptides: keep_decoys,
                ..Default::default()
            };
            let mut ids = vec![
                peptide(0.9, "decoy", "AK", true),
                peptide(0.8, "target", "CK", true),
            ];
            for id in &mut ids {
                id.score_type = "target_decoy".into();
                id.validate().unwrap();
            }
            let before = ids.clone();
            assert!(options.apply_basic_peptide_level(&mut ids).is_err());
            assert_eq!(ids, before);

            // Other methods append _score, so this score type is safe there.
            options.apply_basic_peptides(&mut ids).unwrap();
            for id in &ids {
                id.validate().unwrap();
                for hit in &id.hits {
                    assert!(hit.metadata.contains_key("target_decoy_score"));
                    assert_ne!(hit.target_decoy_type().unwrap(), TargetDecoyType::Unknown);
                }
            }
        }
    }
}

#[test]
fn errors_limits_and_late_group_failures_are_atomic() {
    let options = FalseDiscoveryRate::default();
    for score in [f64::NAN, f64::INFINITY] {
        assert!(options.calculate_legacy(&[score], &[], true).is_err());
        assert!(
            options
                .calculate_basic(&[ScoreLabel::new(score, true)], true)
                .is_err()
        );
    }
    assert!(
        options
            .calculate_basic(
                &[ScoreLabel {
                    score: 1.,
                    target_fraction: 1.1
                }],
                true
            )
            .is_err()
    );
    assert!(options.calculate_estimated(&[-0.1], false).is_err());
    assert!(
        options
            .roc_n(&[ScoreLabel::new(1., false)], true, None)
            .is_err()
    );
    assert!(
        FalseDiscoveryRate {
            max_hits: 1,
            ..options
        }
        .calculate_basic(
            &[ScoreLabel::new(1., true), ScoreLabel::new(2., false)],
            true
        )
        .is_err()
    );
    let mut missing = vec![
        peptide(1., "target", "PEPTIDE", true),
        PeptideIdentification {
            hits: vec![PeptideHit::default()],
            ..Default::default()
        },
    ];
    let before = missing.clone();
    assert!(options.apply_peptides(&mut missing, false).is_err());
    assert_eq!(missing, before);
    let mut mixed = vec![
        peptide(1., "target", "PEPTIDE", true),
        peptide(2., "decoy", "PEPTIDE", false),
    ];
    let before = mixed.clone();
    assert!(options.apply_peptides(&mut mixed, false).is_err());
    assert_eq!(mixed, before);
    let mut groups = ProteinIdentification {
        hits: vec![protein(1., "target", "A")],
        indistinguishable_groups: vec![ProteinGroup {
            probability: 1.,
            accessions: vec!["missing".into()],
            ..Default::default()
        }],
        ..Default::default()
    };
    let before = groups.clone();
    assert!(options.apply_basic_protein(&mut groups, true).is_err());
    assert_eq!(groups, before);
    let mut ids = vec![
        peptide(1., "target", "PEPTIDE", true),
        peptide(2., "decoy", "PEPTIDE", true),
    ];
    ids[1].identifier = "other".into();
    let before = ids.clone();
    assert!(
        FalseDiscoveryRate {
            max_groups: 1,
            treat_runs_separately: true,
            ..options
        }
        .apply_peptides(&mut ids, false)
        .is_err()
    );
    assert_eq!(ids, before);
}

#[test]
fn separate_protein_searches_do_not_need_labels_and_keep_reverse_records() {
    let targets: Vec<_> = table(include_str!("data/fdr_xtandem_fwd_proteins.tsv"))
        .iter()
        .map(|r| ProteinHit::new(r[1].parse().unwrap(), 0, r[0], "").unwrap())
        .collect();
    let decoys: Vec<_> = table(include_str!("data/fdr_xtandem_rev_proteins.tsv"))
        .iter()
        .map(|r| ProteinHit::new(r[1].parse().unwrap(), 0, r[0], "").unwrap())
        .collect();
    let mut forward = vec![ProteinIdentification {
        score_type: "XTandem".into(),
        higher_score_better: false,
        hits: targets,
        ..Default::default()
    }];
    let reverse = vec![ProteinIdentification {
        score_type: "XTandem".into(),
        higher_score_better: false,
        hits: decoys,
        ..Default::default()
    }];
    let before = reverse.clone();
    FalseDiscoveryRate::default()
        .apply_separate_proteins(&mut forward, &reverse)
        .unwrap();
    assert_eq!(reverse, before);
    let mut checked = 0;
    for hit in &forward[0].hits {
        if hit.metadata["XTandem_score"].as_f64().unwrap() == 0. {
            close(hit.score, 0.897384, 1e-6);
            checked += 1;
        }
    }
    assert_eq!(checked, 398);
    let mut peptides = vec![peptide(10., "target", "AK", true)];
    let mut decoys = vec![peptide(9., "decoy", "CK", true)];
    FalseDiscoveryRate {
        add_decoy_peptides: true,
        ..Default::default()
    }
    .apply_separate_peptides(&mut peptides, &mut decoys)
    .unwrap();
    assert_eq!(decoys[0].score_type, "q-value");
    assert!(!decoys[0].higher_score_better);
    assert!(decoys[0].hits[0].metadata.contains_key("score_score"));
}

#[test]
fn source_picked_mixed_group_order_contribution_and_checked_pool_limits() {
    fn record(decoy_first: bool) -> ProteinIdentification {
        let accessions = if decoy_first {
            vec!["decoy_B".into(), "A".into()]
        } else {
            vec!["A".into(), "decoy_B".into()]
        };
        ProteinIdentification {
            score_type: "score".into(),
            hits: vec![
                protein(10., "target", "A"),
                protein(9., "decoy", "decoy_A"),
                protein(11., "decoy", "decoy_B"),
            ],
            indistinguishable_groups: vec![
                ProteinGroup {
                    probability: 11.,
                    accessions: vec!["decoy_B".into()],
                    ..Default::default()
                },
                ProteinGroup {
                    probability: 10.,
                    accessions: vec!["A".into()],
                    ..Default::default()
                },
                ProteinGroup {
                    probability: 9.,
                    accessions,
                    ..Default::default()
                },
            ],
            ..Default::default()
        }
    }
    let mut a = record(true);
    let mut b = record(false);
    let options = FalseDiscoveryRate::default();
    options
        .apply_picked_protein(&mut a, DecoyAffix::Prefix("decoy_"), true)
        .unwrap();
    options
        .apply_picked_protein(&mut b, DecoyAffix::Prefix("decoy_"), true)
        .unwrap();
    assert_eq!(a.indistinguishable_groups[2].probability, 1.);
    close(b.indistinguishable_groups[2].probability, 2. / 3., 1e-14);
    let mut forward = vec![peptide(1., "target", "AK", true)];
    let mut reverse = vec![peptide(1., "decoy", "CK", true)];
    let before = forward.clone();
    assert!(
        FalseDiscoveryRate {
            max_records: 1,
            ..options
        }
        .apply_separate_peptides(&mut forward, &mut reverse)
        .is_err()
    );
    assert_eq!(forward, before);
    assert!(options.calculate_basic(&[], true).unwrap().is_empty());
    assert!(options.calculate_legacy(&[], &[], true).unwrap().is_empty());
    assert_eq!(
        options
            .calculate_legacy(&[], &[1.], true)
            .unwrap()
            .value(1.)
            .unwrap(),
        1.
    );
    assert_eq!(
        options
            .calculate_legacy(&[1.], &[], true)
            .unwrap()
            .value(1.)
            .unwrap(),
        0.
    );
}
