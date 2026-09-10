// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::analysis::protein_inference::{AggregationMethod, BasicProteinInference};
use openms::analysis::scores::ScoreType;
use openms::chemistry::AASequence;
use openms::identification::{
    PeptideEvidence, PeptideHit, PeptideIdentification, ProteinHit, ProteinIdentification,
    TargetDecoyType,
};
use openms::kernel::{ConsensusFeature, ConsensusMap};
use openms::metadata::MetaValue;

fn hit(sequence: &str, score: f64, charge: i32, accessions: &[&str]) -> PeptideHit {
    PeptideHit {
        sequence: AASequence::parse(sequence).unwrap(),
        score,
        charge,
        evidences: accessions
            .iter()
            .map(|accession| PeptideEvidence {
                protein_accession: (*accession).into(),
                ..Default::default()
            })
            .collect(),
        ..Default::default()
    }
}
fn id(hit: PeptideHit) -> PeptideIdentification {
    PeptideIdentification {
        identifier: "run".into(),
        score_type: "hyperscore".into(),
        hits: vec![hit],
        ..Default::default()
    }
}
fn proteins(accessions: &[&str]) -> Vec<ProteinIdentification> {
    vec![ProteinIdentification {
        identifier: "run".into(),
        hits: accessions
            .iter()
            .map(|a| {
                let mut hit = ProteinHit {
                    accession: (*a).into(),
                    ..Default::default()
                };
                hit.set_target_decoy_type(TargetDecoyType::Target).unwrap();
                hit
            })
            .collect(),
        ..Default::default()
    }]
}
fn fixture() -> (Vec<PeptideIdentification>, Vec<ProteinIdentification>) {
    let mut p = proteins(&[]);
    let mut ids = Vec::new();
    for line in include_str!("data/protein_inference_merger.tsv")
        .lines()
        .skip(2)
    {
        let row: Vec<_> = line.split('\t').collect();
        if row[0] == "protein" {
            let mut h = ProteinHit {
                accession: row[2].into(),
                score: row[3].parse().unwrap(),
                ..Default::default()
            };
            h.metadata.insert("target_decoy".into(), row[6].into());
            p[0].hits.push(h);
        } else {
            let index = row[1].parse::<usize>().unwrap();
            let accessions: Vec<_> = if row[5] == "-" {
                vec![]
            } else {
                row[5].split(',').collect()
            };
            let mut h = hit(
                row[2],
                row[3].parse().unwrap(),
                row[4].parse().unwrap(),
                &accessions,
            );
            if row[6] != "-" {
                h.metadata.insert("target_decoy".into(), row[6].into());
            }
            h.metadata
                .insert("protein_references".into(), row[7].into());
            h.metadata.insert(
                "XTandem_score".into(),
                MetaValue::try_from(row[8].parse::<f64>().unwrap()).unwrap(),
            );
            if index == ids.len() {
                let mut next = id(h);
                next.score_type = "Posterior Error Probability".into();
                next.higher_score_better = false;
                ids.push(next);
            } else {
                ids[index].hits.push(h);
            }
        }
    }
    (ids, p)
}
fn scores(proteins: &[ProteinIdentification]) -> Vec<f64> {
    proteins[0].hits.iter().map(|h| h.score).collect()
}
fn counts(proteins: &[ProteinIdentification]) -> Vec<i64> {
    proteins[0]
        .hits
        .iter()
        .map(|h| h.metadata["nr_found_peptides"].as_i64().unwrap())
        .collect()
}
fn close(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        assert!((a - b).abs() < 1e-12, "{a} != {b}");
    }
}

#[test]
fn source_defaults_and_merger_golden() {
    let algorithm = BasicProteinInference::default();
    assert_eq!(algorithm.aggregation, AggregationMethod::Best);
    assert_eq!(algorithm.min_peptides_per_protein, 1);
    assert!(algorithm.annotate_indistinguishable_groups);
    assert!(!algorithm.greedy_group_resolution);
    let (mut ids, mut p) = fixture();
    let report = algorithm.run(&mut ids, &mut p).unwrap();
    close(&scores(&p), &[0.6, 0.6, 0.8, 0.6, 0.9]);
    assert_eq!(counts(&p), [1, 1, 2, 1, 1]);
    close(
        &p[0]
            .indistinguishable_groups
            .iter()
            .map(|g| g.probability)
            .collect::<Vec<_>>(),
        &[0.9, 0.8, 0.6, 0.6],
    );
    assert_eq!(p[0].indistinguishable_groups[3].accessions, ["A", "D"]);
    assert_eq!(p[0].score_type, "Posterior Probability");
    assert!(p[0].higher_score_better);
    assert_eq!(
        ids.iter().map(|i| i.hits.len()).collect::<Vec<_>>(),
        [1, 1, 1, 1]
    );
    close(
        &ids.iter().map(|i| i.hits[0].score).collect::<Vec<_>>(),
        &[0.4, 0.2, 0.4, 0.1],
    );
    assert_eq!(report.proteins_removed, 1);
    assert_eq!(report.peptide_hits_removed, 4);
    assert_eq!(
        p[0].search_parameters.metadata["TOPPProteinInference:aggregation_method"]
            .as_str()
            .unwrap(),
        "best"
    );
}

#[test]
fn source_greedy_merger_golden_including_minimum_zero() {
    let (mut ids, mut p) = fixture();
    let config = BasicProteinInference {
        min_peptides_per_protein: 0,
        greedy_group_resolution: true,
        ..Default::default()
    };
    let report = config.run(&mut ids, &mut p).unwrap();
    close(&scores(&p), &[0.6, 0.6, 0.8, 0.9]);
    assert_eq!(counts(&p), [1, 1, 2, 1]);
    assert_eq!(
        p[0].hits
            .iter()
            .map(|p| p.accession.as_str())
            .collect::<Vec<_>>(),
        ["A", "D", "C", "BSA4"]
    );
    close(
        &p[0]
            .indistinguishable_groups
            .iter()
            .map(|g| g.probability)
            .collect::<Vec<_>>(),
        &[0.9, 0.8, 0.6],
    );
    assert_eq!(
        ids[2].hits[0]
            .protein_accessions()
            .into_iter()
            .collect::<Vec<_>>(),
        ["C"]
    );
    assert_eq!(
        ids[2].hits[0].metadata["protein_references"]
            .as_str()
            .unwrap(),
        "non-unique"
    ); // source keeps annotation
    assert_eq!(report.resolved_peptide_hits, 1);
}

#[test]
fn source_raw_score_golden_and_exact_restoration() {
    for greedy in [false, true] {
        let (mut ids, mut p) = fixture();
        let config = BasicProteinInference {
            score_type: Some(ScoreType::Raw),
            greedy_group_resolution: greedy,
            ..Default::default()
        };
        config.run(&mut ids, &mut p).unwrap();
        close(
            &scores(&p),
            if greedy {
                &[2.5, 2.5, 5.0, 10.0]
            } else {
                &[2.5, 2.5, 5.0, 2.5, 10.0]
            },
        );
        assert_eq!(p[0].score_type, "XTandem");
        for id in &ids {
            assert_eq!(id.score_type, "Posterior Error Probability");
            assert!(!id.higher_score_better);
        }
        close(
            &ids.iter().map(|i| i.hits[0].score).collect::<Vec<_>>(),
            &[0.4, 0.2, 0.4, 0.1],
        );
        assert_eq!(ids[0].hits[0].metadata["XTandem"].as_f64().unwrap(), 2.5);
        assert_eq!(
            ids[0].hits[0].metadata["Posterior Error Probability"]
                .as_f64()
                .unwrap(),
            0.4
        );
    }
}

#[test]
fn source_shared_peptide_exclusion_and_minimum_cleanup() {
    let (mut ids, mut p) = fixture();
    BasicProteinInference {
        use_shared_peptides: false,
        ..Default::default()
    }
    .run(&mut ids, &mut p)
    .unwrap();
    close(&scores(&p), &[0.8, 0.9]);
    assert_eq!(counts(&p), [1, 1]);
    assert!(ids[0].hits.is_empty());
    assert_eq!(
        ids[2].hits[0]
            .protein_accessions()
            .into_iter()
            .collect::<Vec<_>>(),
        ["C"]
    );
}

#[test]
fn source_mean_and_product_rules_including_zero_contributions() {
    for (aggregation, expected) in [
        (AggregationMethod::Best, 0.8),
        (AggregationMethod::Product, 0.48),
        (AggregationMethod::Mean, 1.4 / 3.0),
    ] {
        let mut ids = vec![
            id(hit("AA", 0.2, 2, &["A"])),
            id(hit("CC", 0.4, 2, &["A"])),
            id(hit("DD", 1.0, 2, &["A"])),
        ];
        for id in &mut ids {
            id.score_type = "PEP".into();
            id.higher_score_better = false;
        }
        let mut p = proteins(&["A"]);
        BasicProteinInference {
            aggregation,
            ..Default::default()
        }
        .run(&mut ids, &mut p)
        .unwrap();
        close(&scores(&p), &[expected]);
        assert_eq!(counts(&p), [3]);
    }
    assert_eq!(
        AggregationMethod::from_source_name("sum").unwrap(),
        AggregationMethod::Mean
    );
    assert_eq!(
        AggregationMethod::from_source_name("maximum").unwrap(),
        AggregationMethod::Best
    );
    assert!(AggregationMethod::from_source_name("median").is_err());
    let mut ids = vec![
        id(hit("AA", -3.0, 2, &["A"])),
        id(hit("CC", 0.0, 2, &["A"])),
        id(hit("DD", 2.0, 2, &["A"])),
    ];
    let mut p = proteins(&["A", "B"]);
    BasicProteinInference {
        aggregation: AggregationMethod::Product,
        min_peptides_per_protein: 0,
        ..Default::default()
    }
    .run(&mut ids, &mut p)
    .unwrap();
    assert_eq!(scores(&p), [2.0, 1.0]);
    assert_eq!(counts(&p), [3, 0]);
}

#[test]
fn representative_charge_and_modification_collapsing_and_first_ties() {
    for (charge, modification, expected_count) in [
        (true, true, 3),
        (false, true, 2),
        (true, false, 2),
        (false, false, 1),
    ] {
        let mut ids = vec![
            id(hit("ACM", 2.0, 2, &["A"])),
            id(hit("ACM", 3.0, 3, &["A"])),
            id(hit("ACM(Oxidation)", 4.0, 2, &["A"])),
            id(hit("ACM", 2.0, 2, &["A"])),
        ];
        ids[0].hits[0]
            .metadata
            .insert("tie_marker".into(), "first".into());
        let mut p = proteins(&["A"]);
        BasicProteinInference {
            treat_charge_variants_separately: charge,
            treat_modification_variants_separately: modification,
            ..Default::default()
        }
        .run(&mut ids, &mut p)
        .unwrap();
        assert_eq!(scores(&p), [4.0]);
        assert_eq!(counts(&p), [expected_count]);
        assert_eq!(ids.len(), 4);
        assert!(ids.iter().all(|id| id.hits.len() == 1));
    }
    // Equal representatives retain the first accession set; no union is fabricated.
    let mut ids = vec![id(hit("AA", 3.0, 2, &["A"])), id(hit("AA", 3.0, 2, &["B"]))];
    let mut p = proteins(&["A", "B"]);
    BasicProteinInference::default()
        .run(&mut ids, &mut p)
        .unwrap();
    assert_eq!(p[0].hits[0].accession, "A");
    assert_eq!(p[0].hits.len(), 1);
    assert!(ids[1].hits.is_empty());
}

#[test]
fn original_main_score_selects_candidate_before_switching() {
    let mut first = hit("AA", 10.0, 2, &["A"]);
    first
        .metadata
        .insert("PEP".into(), MetaValue::try_from(0.8).unwrap());
    let mut second = hit("CC", 9.0, 2, &["B"]);
    second
        .metadata
        .insert("PEP".into(), MetaValue::try_from(0.1).unwrap());
    let mut ids = vec![id(first)];
    ids[0].hits.push(second);
    let mut p = proteins(&["A", "B"]);
    BasicProteinInference {
        score_type: Some(ScoreType::PosteriorErrorProbability),
        ..Default::default()
    }
    .run(&mut ids, &mut p)
    .unwrap();
    assert_eq!(ids[0].hits[0].sequence.as_str(), "AA");
    assert_eq!(ids[0].hits[0].score, 10.0);
    close(&scores(&p), &[0.2]);
    assert_eq!(p[0].hits[0].accession, "A");
}

#[test]
fn score_restore_honors_backup_collisions_without_losing_original() {
    let mut h = hit("AA", 10.0, 2, &["A"]);
    h.metadata
        .insert("hyperscore".into(), MetaValue::try_from(9.0).unwrap());
    h.metadata
        .insert("PEP".into(), MetaValue::try_from(0.2).unwrap());
    let mut ids = vec![id(h)];
    let mut p = proteins(&["A"]);
    let config = BasicProteinInference {
        score_type: Some(ScoreType::PosteriorErrorProbability),
        ..Default::default()
    };
    config.run(&mut ids, &mut p).unwrap();
    assert_eq!(ids[0].score_type, "hyperscore");
    assert_eq!(ids[0].hits[0].score, 10.0);
    assert_eq!(ids[0].hits[0].metadata["hyperscore"].as_f64().unwrap(), 9.0);
    assert_eq!(
        ids[0].hits[0].metadata["hyperscore~"].as_f64().unwrap(),
        10.0
    );
}

#[test]
fn grouping_uses_psm_neighborhoods_and_source_negative_floor() {
    let mut ids = vec![
        id(hit("AA", -3.0, 2, &["A", "C"])),
        id(hit("CC", -2.0, 2, &["B"])),
    ];
    let mut p = proteins(&["C", "B", "A"]);
    BasicProteinInference::default()
        .run(&mut ids, &mut p)
        .unwrap();
    assert_eq!(p[0].indistinguishable_groups.len(), 2);
    assert!(
        p[0].indistinguishable_groups
            .iter()
            .all(|g| g.probability == -1.0)
    );
    assert_eq!(p[0].indistinguishable_groups[0].accessions, ["B"]);
    assert_eq!(p[0].indistinguishable_groups[1].accessions, ["A", "C"]);
}

#[test]
fn greedy_ties_prefer_target_then_current_psm_count_then_lexical_accession() {
    for (decoy, extra, expected) in [(true, false, "Z"), (false, true, "Z"), (false, false, "A")] {
        let mut ids = vec![
            id(hit("AA", 10.0, 2, &["A", "Z"])),
            id(hit("CC", 5.0, 2, &["A"])),
            id(hit("DD", 5.0, 2, &["Z"])),
        ];
        if extra {
            ids.push(id(hit("DD", 5.0, 2, &["Z"])));
        }
        let mut p = proteins(&["Z", "A"]);
        if decoy {
            p[0].hits[1]
                .set_target_decoy_type(TargetDecoyType::Decoy)
                .unwrap();
        }
        BasicProteinInference {
            greedy_group_resolution: true,
            ..Default::default()
        }
        .run(&mut ids, &mut p)
        .unwrap();
        assert_eq!(
            ids[0].hits[0]
                .protein_accessions()
                .into_iter()
                .collect::<Vec<_>>(),
            [expected]
        );
    }
}

#[test]
fn negative_greedy_scores_outrank_unscored_proteins() {
    let mut shared = hit("AA", -5.0, 2, &["A", "B"]);
    shared
        .metadata
        .insert("protein_references".into(), "non-unique".into());
    let mut unique = hit("CC", -2.0, 2, &["A"]);
    unique
        .metadata
        .insert("protein_references".into(), "unique".into());
    let mut ids = vec![id(shared), id(unique)];
    let mut p = proteins(&["A", "B"]);
    BasicProteinInference {
        min_peptides_per_protein: 0,
        use_shared_peptides: false,
        greedy_group_resolution: true,
        ..Default::default()
    }
    .run(&mut ids, &mut p)
    .unwrap();
    assert_eq!(scores(&p), [-2.0]);
    assert_eq!(p[0].hits[0].accession, "A");
    assert_eq!(
        ids[0].hits[0]
            .protein_accessions()
            .into_iter()
            .collect::<Vec<_>>(),
        ["A"]
    );
    assert_eq!(p[0].indistinguishable_groups[0].probability, -2.0);
}

#[test]
fn count_suppression_does_not_disable_minimum_filter() {
    let mut ids = vec![id(hit("AA", 3.0, 2, &["A"]))];
    let mut p = proteins(&["A", "B"]);
    p[0].hits[1]
        .metadata
        .insert("nr_found_peptides".into(), 100i64.into());
    BasicProteinInference {
        skip_count_annotation: true,
        ..Default::default()
    }
    .run(&mut ids, &mut p)
    .unwrap();
    assert_eq!(p[0].hits.len(), 1);
    assert!(!p[0].hits[0].metadata.contains_key("nr_found_peptides"));
}

#[test]
fn multiple_runs_are_independent_and_can_use_different_main_scores() {
    let mut ids = vec![id(hit("AA", 3.0, 2, &["A"]))];
    let mut p = proteins(&["A"]);
    let mut other = id(hit("CC", 0.2, 2, &["A"]));
    other.identifier = "other".into();
    other.score_type = "PEP".into();
    other.higher_score_better = false;
    ids.push(other);
    let mut other = p[0].clone();
    other.identifier = "other".into();
    p.push(other);
    BasicProteinInference::default()
        .run(&mut ids, &mut p)
        .unwrap();
    assert_eq!(p[0].hits[0].score, 3.0);
    close(&[p[1].hits[0].score], &[0.8]);
    assert_eq!(p[0].score_type, "hyperscore");
    assert_eq!(p[1].score_type, "Posterior Probability");
}

#[test]
fn nonfinite_undefined_and_invalid_paths_are_atomic() {
    let (base_ids, base_p) = fixture();
    for aggregation in [AggregationMethod::Best, AggregationMethod::Mean] {
        let mut ids = base_ids.clone();
        let mut p = base_p.clone();
        assert!(
            BasicProteinInference {
                aggregation,
                min_peptides_per_protein: 0,
                ..Default::default()
            }
            .run(&mut ids, &mut p)
            .is_err()
        );
        assert_eq!(ids, base_ids);
        assert_eq!(p, base_p);
    }
    for case in 0..9 {
        let mut ids = vec![id(hit("AA", 3.0, 2, &["A"]))];
        let mut p = proteins(&["A"]);
        let mut config = BasicProteinInference::default();
        match case {
            0 => {
                p.push(p[0].clone());
            }
            1 => {
                let duplicate = p[0].hits[0].clone();
                p[0].hits.push(duplicate);
            }
            2 => {
                ids[0].identifier = "unknown".into();
            }
            3 => {
                ids.push(id(hit("AA", 4.0, 3, &["B"])));
            }
            4 => {
                ids[0].score_type = "PEP".into();
                ids[0].higher_score_better = false;
            }
            5 => {
                config.max_input_hits = 1;
            }
            6 => {
                ids.push(id(hit("CC", f64::MAX, 2, &["A"])));
                ids[0].hits[0].score = f64::MAX;
                config.aggregation = AggregationMethod::Mean;
            }
            7 => {
                ids[0].hits[0]
                    .metadata
                    .insert("protein_references".into(), 10i64.into());
                config.use_shared_peptides = false;
            }
            _ => {
                config.greedy_group_resolution = true;
                p[0].hits[0].metadata.remove("target_decoy");
            }
        }
        let before_ids = ids.clone();
        let before_p = p.clone();
        assert!(config.run(&mut ids, &mut p).is_err(), "case {case}");
        assert_eq!(ids, before_ids);
        assert_eq!(p, before_p);
    }
}

#[test]
fn empty_inputs_lower_better_and_group_suppression() {
    assert_eq!(
        BasicProteinInference::default()
            .run(&mut [], &mut [])
            .unwrap()
            .protein_runs,
        0
    );
    let mut p = proteins(&["A"]);
    BasicProteinInference::default()
        .run(&mut [], &mut p)
        .unwrap();
    assert!(p[0].hits.is_empty());
    let mut ids = vec![id(hit("AA", 0.2, 2, &["A"])), id(hit("CC", 0.1, 2, &["A"]))];
    for id in &mut ids {
        id.score_type = "q-value".into();
        id.higher_score_better = false;
    }
    let mut p = proteins(&["A"]);
    BasicProteinInference {
        annotate_indistinguishable_groups: false,
        ..Default::default()
    }
    .run(&mut ids, &mut p)
    .unwrap();
    assert_eq!(scores(&p), [0.1]);
    assert!(!p[0].higher_score_better);
    assert!(p[0].indistinguishable_groups.is_empty());
    assert!(
        BasicProteinInference {
            greedy_group_resolution: true,
            ..Default::default()
        }
        .run(&mut ids, &mut p)
        .is_err()
    );
}

#[test]
fn single_run_switches_before_best_selection_and_retains_other_candidates() {
    let mut a = hit("AA", 10.0, 2, &["A"]);
    a.metadata
        .insert("PEP".into(), MetaValue::try_from(0.8).unwrap());
    let mut b = hit("CC", 9.0, 2, &["B"]);
    b.metadata
        .insert("PEP".into(), MetaValue::try_from(0.1).unwrap());
    let mut first = id(a);
    first.hits.push(b);
    let mut second = id(hit("DD", 8.0, 2, &["A"]));
    second.hits[0]
        .metadata
        .insert("PEP".into(), MetaValue::try_from(0.2).unwrap());
    let mut unrelated = id(hit("EE", 5.0, 2, &["OTHER"]));
    unrelated.identifier = "other".into();
    let before_unrelated = unrelated.clone();
    let mut ids = vec![first, second, unrelated];
    let mut p = proteins(&["A", "B"]);
    BasicProteinInference {
        score_type: Some(ScoreType::PosteriorErrorProbability),
        ..Default::default()
    }
    .run_single(&mut ids, &mut p[0])
    .unwrap();
    close(&scores(&p), &[0.8, 0.9]);
    assert_eq!(counts(&p), [1, 1]);
    assert_eq!(ids[0].hits.len(), 2);
    assert_eq!(ids[0].hits[0].sequence.as_str(), "CC");
    assert_eq!(ids[0].hits[0].score, 9.0);
    assert_eq!(ids[0].hits[1].sequence.as_str(), "AA");
    assert_eq!(ids[0].hits[1].score, 10.0);
    assert_eq!(ids[2], before_unrelated);
    assert_eq!(ids[0].score_type, "hyperscore");
}

#[test]
fn single_run_greedy_cleanup_counts_retained_lower_candidates() {
    let mut first = id(hit("AA", 0.8, 2, &["A"]));
    first.hits.push(hit("CC", 0.7, 2, &["B"]));
    let mut ids = vec![first];
    let mut p = proteins(&["A", "B"]);
    BasicProteinInference {
        aggregation: AggregationMethod::Product,
        min_peptides_per_protein: 0,
        greedy_group_resolution: true,
        ..Default::default()
    }
    .run_single(&mut ids, &mut p[0])
    .unwrap();
    assert_eq!(scores(&p), [0.8, 1.0]);
    assert_eq!(counts(&p), [1, 0]);
    assert_eq!(ids[0].hits.len(), 2);
    assert_eq!(p[0].indistinguishable_groups[0].accessions, ["B"]);
}

#[test]
fn consensus_union_ignores_run_ids_optionally_includes_unassigned_and_sorts_proteins() {
    for include_unassigned in [false, true] {
        let mut first = id(hit("AA", 3.0, 2, &["A"]));
        first.identifier = "different".into();
        let mut alternate = hit("CC", 1.0, 2, &["B"]);
        alternate
            .metadata
            .insert("extra".into(), "preserved".into());
        first.hits.push(alternate);
        let mut unassigned = id(hit("DD", 5.0, 2, &["B"]));
        unassigned.identifier = "third".into();
        unassigned.hits.push(hit("EE", 4.0, 2, &["B"]));
        let mut feature = ConsensusFeature::default();
        feature.rt = 42.0;
        feature.peptide_identifications.push(first);
        let mut map = ConsensusMap {
            features: vec![feature],
            unassigned_peptide_identifications: vec![unassigned.clone()],
            protein_identifications: proteins(&["OLD"]),
            ..Default::default()
        };
        let before_proteins = map.protein_identifications.clone();
        let mut p = proteins(&["A", "B"]);
        BasicProteinInference::default()
            .run_consensus_map(&mut map, &mut p[0], include_unassigned)
            .unwrap();
        assert_eq!(map.features[0].rt, 42.0);
        assert_eq!(
            map.features[0].peptide_identifications[0].identifier,
            "different"
        );
        assert_eq!(map.features[0].peptide_identifications[0].hits.len(), 1);
        assert_eq!(map.protein_identifications, before_proteins);
        if include_unassigned {
            assert_eq!(scores(&p), [5.0, 3.0]);
            assert_eq!(map.unassigned_peptide_identifications[0].hits.len(), 1);
            assert_eq!(
                map.unassigned_peptide_identifications[0].identifier,
                "third"
            );
        } else {
            assert_eq!(scores(&p), [3.0]);
            assert_eq!(map.unassigned_peptide_identifications[0], unassigned);
        }
    }
}

#[test]
fn unused_tilde_metadata_is_not_interpreted_as_score() {
    let mut h = hit("AA", 10.0, 2, &["A"]);
    h.metadata
        .insert("hyperscore~".into(), "unrelated annotation".into());
    h.metadata
        .insert("PEP".into(), MetaValue::try_from(0.2).unwrap());
    let mut ids = vec![id(h)];
    let mut p = proteins(&["A"]);
    BasicProteinInference {
        score_type: Some(ScoreType::PosteriorErrorProbability),
        ..Default::default()
    }
    .run(&mut ids, &mut p)
    .unwrap();
    assert_eq!(ids[0].hits[0].score, 10.0);
    assert_eq!(
        ids[0].hits[0].metadata["hyperscore~"].as_str().unwrap(),
        "unrelated annotation"
    );
}

#[test]
fn grouping_keeps_distinct_psm_neighbors_even_for_identical_sequences() {
    let mut ids = vec![id(hit("AA", 0.9, 2, &["A"])), id(hit("AA", 0.8, 2, &["B"]))];
    let mut p = proteins(&["A", "B"]);
    BasicProteinInference {
        aggregation: AggregationMethod::Product,
        min_peptides_per_protein: 0,
        ..Default::default()
    }
    .run(&mut ids, &mut p)
    .unwrap();
    assert_eq!(scores(&p), [0.9, 1.0]);
    assert_eq!(counts(&p), [1, 0]);
    assert_eq!(p[0].indistinguishable_groups.len(), 2);
    assert!(
        p[0].indistinguishable_groups
            .iter()
            .all(|group| group.accessions.len() == 1)
    );
}

#[test]
fn adapter_errors_leave_map_and_protein_records_unchanged() {
    let mut feature = ConsensusFeature::default();
    feature
        .peptide_identifications
        .push(id(hit("AA", 3.0, 2, &["A"])));
    let mut map = ConsensusMap {
        features: vec![feature],
        ..Default::default()
    };
    let mut p = proteins(&["A", "B"]);
    let before_map = map.clone();
    let before_p = p.clone();
    let config = BasicProteinInference {
        min_peptides_per_protein: 0,
        ..Default::default()
    };
    assert!(config.run_consensus_map(&mut map, &mut p[0], true).is_err());
    assert_eq!(map, before_map);
    assert_eq!(p, before_p);
    let mut ids = vec![id(hit("AA", 3.0, 2, &["A"]))];
    let before_ids = ids.clone();
    assert!(config.run_single(&mut ids, &mut p[0]).is_err());
    assert_eq!(ids, before_ids);
    assert_eq!(p, before_p);
}
