// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Executed differential (tier 1) for the `FalseDiscoveryRate.h` evaluation
//! members: `diffEstimatedEmpirical`, `rocN` on sorted pairs, on peptide
//! identification lists and on consensus maps (with and without a run
//! identifier and unassigned identifications), and the three
//! `applyEvaluateProteinIDs` overloads.
//!
//! `tests/data/fdr_evaluation_oracle.tsv` is the output of
//! `../oracle/I2-fdr/driver.cpp` linked against the Linux x86_64 Release
//! `libOpenMS.so` on ibminode06 (two runs, identical bytes); it holds every
//! input and result at `%.17g`. Contract: bit-identical results.

use openms::analysis::false_discovery_rate::FalseDiscoveryRate;
use openms::identification::{
    PeptideHit, PeptideIdentification, ProteinHit, ProteinIdentification,
};
use openms::kernel::{ConsensusFeature, ConsensusMap};
use openms::metadata::MetaValue;
use std::path::Path;

fn number(text: &str) -> f64 {
    text.parse().unwrap()
}

fn parse_pairs(fields: &[&str]) -> Vec<(f64, f64)> {
    fields[2..]
        .iter()
        .map(|f| {
            let (s, l) = f.split_once(':').unwrap();
            (number(s), number(l))
        })
        .collect()
}

fn fdr(conservative: bool, use_all_hits: bool) -> FalseDiscoveryRate {
    FalseDiscoveryRate {
        conservative,
        use_all_hits,
        ..FalseDiscoveryRate::default()
    }
}

fn assert_bits(actual: f64, expected: f64, what: &str) {
    assert!(
        actual.to_bits() == expected.to_bits(),
        "{what}: {actual:.17e} vs {expected:.17e}"
    );
}

#[test]
fn evaluation_members_match_the_release_build() {
    let text = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/fdr_evaluation_oracle.tsv"),
    )
    .unwrap();
    let mut lines = text.lines().peekable();
    let (mut checked, mut logged, mut context) = (0usize, 0usize, String::new());
    let mut calculator = FalseDiscoveryRate::default();
    let mut pairs: Vec<(f64, f64)> = Vec::new();
    let mut ids: Vec<PeptideIdentification> = Vec::new();
    let mut kind = String::new();
    while let Some(line) = lines.next() {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "case" => {
                context = line.to_owned();
                kind = f[1].to_owned();
                calculator = if kind.starts_with("roc_") {
                    fdr(true, f[2] == "1")
                } else {
                    fdr(f[2] == "1", false)
                };
            }
            "pairs" => pairs = parse_pairs(&f),
            "ids" => {
                let n: usize = f[1].parse().unwrap();
                ids.clear();
                for _ in 0..n {
                    let g: Vec<&str> = lines.next().unwrap().split('\t').collect();
                    let hits = g[3..]
                        .iter()
                        .map(|h| {
                            let (s, td) = h.split_once(':').unwrap();
                            let mut hit = PeptideHit {
                                score: number(s),
                                ..PeptideHit::default()
                            };
                            hit.metadata
                                .insert("target_decoy".into(), MetaValue::from(td));
                            hit
                        })
                        .collect();
                    ids.push(PeptideIdentification {
                        identifier: g[1].into(),
                        higher_score_better: g[2] == "1",
                        score_type: "score".into(),
                        hits,
                        ..Default::default()
                    });
                }
            }
            "diff" => {
                let got = calculator.diff_estimated_empirical(&pairs, number(f[1]));
                assert_bits(got, number(f[2]), &format!("{context} diff {}", f[1]));
                checked += 1;
            }
            "roc" if kind == "diff_roc" => {
                let got = calculator.roc_n_sorted(&pairs, f[1].parse().unwrap());
                assert_bits(got, number(f[2]), &format!("{context} roc {}", f[1]));
                checked += 1;
            }
            "eval" | "eval2" if kind == "eval_pairs" => {
                let mut copy = pairs.clone();
                let (cut, fp, w) = if f[0] == "eval" {
                    (1.0, 50, 0.2)
                } else {
                    (0.5, 10, 0.4)
                };
                let got = calculator
                    .apply_evaluate_score_labels(&mut copy, cut, fp, w)
                    .unwrap();
                assert_bits(got, number(f[1]), &format!("{context} {}", f[0]));
                checked += 1;
            }
            "eval" | "eval_vec" => {
                let mut protein = ProteinIdentification {
                    score_type: "Posterior Probability".into(),
                    higher_score_better: true,
                    ..Default::default()
                };
                for (i, &(score, label)) in pairs.iter().enumerate() {
                    let mut hit = ProteinHit {
                        accession: format!("P{i}"),
                        score,
                        ..ProteinHit::default()
                    };
                    hit.metadata.insert(
                        "target_decoy".into(),
                        MetaValue::from(if label > 0.0 { "target" } else { "decoy" }),
                    );
                    protein.hits.push(hit);
                }
                let got = if f[0] == "eval" {
                    calculator.apply_evaluate_protein_ids(&protein, 1.0, 50, 0.2)
                } else {
                    calculator.apply_evaluate_protein_ids_first(&[protein], 1.0, 50, 0.2)
                }
                .unwrap();
                assert_bits(got, number(f[1]), &format!("{context} {}", f[0]));
                checked += 1;
            }
            "roc" | "roc_runB" if kind == "roc_ids" => {
                let run = (f[0] == "roc_runB").then_some("runB");
                let got = calculator
                    .roc_n_peptides(&ids, f[1].parse().unwrap(), run)
                    .unwrap();
                assert_bits(got, number(f[2]), &format!("{context} {} {}", f[0], f[1]));
                checked += 1;
            }
            "roc" | "roc_runB" if kind == "roc_cmap" => {
                let mut map = ConsensusMap::default();
                for i in 0..20 {
                    let mut feature = ConsensusFeature::new();
                    feature.base.peptide_identifications = ids[i * 3..i * 3 + 3].to_vec();
                    map.features.push(feature);
                }
                map.unassigned_peptide_identifications = ids[60..].to_vec();
                let run = (f[0] == "roc_runB").then_some("runB");
                let got = calculator
                    .roc_n_consensus(&map, f[1].parse().unwrap(), run, f[2] == "1")
                    .unwrap();
                assert_bits(
                    got,
                    number(f[3]),
                    &format!("{context} {} {} {}", f[0], f[1], f[2]),
                );
                checked += 1;
            }
            // The source's OPENMS_LOG_INFO lines of the single-run and pairs
            // overloads, flushed at the end of the driver's output.
            other if other.starts_with("Evaluation of protein probabilities") => logged += 1,
            other => panic!("unexpected oracle record {other} in {context}"),
        }
    }
    assert_eq!(checked, 2 * 3 * 2 * (2 + 3 + 2 + 2) + 2 * 2 * (6 + 12));
    assert_eq!(logged, 2 * 3 * 2 * 3);
}

#[test]
fn evaluation_edge_cases() {
    // Rust-only, from the source text: empty inputs, missing annotations and
    // the posterior-probability requirement.
    let calculator = FalseDiscoveryRate::default();
    assert_eq!(calculator.diff_estimated_empirical(&[], 1.0), 1.0);
    assert_eq!(calculator.roc_n_sorted(&[], 50), 0.0);
    assert_eq!(calculator.roc_n_sorted(&[(0.9, 1.0), (0.8, 1.0)], 50), 1.0);
    let pep = ProteinIdentification {
        score_type: "Posterior Error Probability".into(),
        ..Default::default()
    };
    assert!(
        calculator
            .apply_evaluate_protein_ids(&pep, 1.0, 50, 0.2)
            .is_err()
    );
    assert!(
        calculator
            .apply_evaluate_protein_ids_first(&[], 1.0, 50, 0.2)
            .is_err()
    );
    assert!(calculator.roc_n_peptides(&[], 0, None).is_err());
    let unannotated = PeptideIdentification {
        hits: vec![PeptideHit::default()],
        ..Default::default()
    };
    assert!(calculator.roc_n_peptides(&[unannotated], 0, None).is_err());
}
