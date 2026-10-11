// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Executed differential (tier 1) for the source-exact `FalseDiscoveryRate.h`
//! members: `applyBasic(ids, higher, charge, identifier, only_best_per_pep)`,
//! `applyBasic(run_info, ids)`, `applyBasic(ConsensusMap&, bool)`, both
//! `applyBasicPeptideLevel` overloads, `DecoyStringHelper::findDecoyString`
//! and `applyPickedProteinFDR` with an empty decoy string.
//!
//! `tests/data/fdr_basic_oracle.tsv` is the output of
//! `../oracle/I2-fdr/driver2.cpp` linked against the Linux x86_64 Release
//! `libOpenMS.so` on ibminode06 (two runs, identical bytes). Every case prints
//! its input identifications and the resulting ones (score bits, charge,
//! sequence and every meta value), the source's exception, or `crash` where
//! the source's undefined behaviour (reading the first hit of an emptied
//! record) killed the child process that ran the case. Contract: identical
//! records and bit-identical scores; the port returns an error exactly where
//! the source threw or crashed.

use openms::analysis::false_discovery_rate::{FalseDiscoveryRate, FdrOutput, find_decoy_string};
use openms::chemistry::AASequence;
use openms::identification::{
    PeptideHit, PeptideIdentification, ProteinHit, ProteinIdentification,
};
use openms::kernel::{ConsensusFeature, ConsensusMap};
use openms::metadata::{MetaValue, MetaValueData};
use std::iter::Peekable;
use std::path::Path;
use std::str::Lines;

type Input<'a> = Peekable<Lines<'a>>;

fn read_ids(lines: &mut Input<'_>, tag: &str) -> Vec<PeptideIdentification> {
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    assert_eq!(head[0], tag);
    let n: usize = head[1].parse().unwrap();
    (0..n)
        .map(|_| {
            let p: Vec<&str> = lines.next().unwrap().split('\t').collect();
            assert_eq!(p[0], "pid");
            let hits = (0..p[4].parse::<usize>().unwrap())
                .map(|_| {
                    let h: Vec<&str> = lines.next().unwrap().split('\t').collect();
                    assert_eq!(h[0], "hit");
                    let mut hit = PeptideHit {
                        score: h[1].parse().unwrap(),
                        charge: h[2].parse().unwrap(),
                        sequence: AASequence::parse(h[3]).unwrap(),
                        ..PeptideHit::default()
                    };
                    for kv in &h[5..] {
                        let (key, typed) = kv.split_once('=').unwrap();
                        let (kind, value) = typed.split_once(':').unwrap();
                        let value = match kind {
                            "d" => MetaValue::try_from(value.parse::<f64>().unwrap()).unwrap(),
                            "i" => MetaValue::from(value.parse::<i64>().unwrap()),
                            _ => MetaValue::from(value),
                        };
                        hit.metadata.insert(key.into(), value);
                    }
                    hit
                })
                .collect();
            PeptideIdentification {
                identifier: p[1].into(),
                score_type: p[2].into(),
                higher_score_better: p[3] == "1",
                hits,
                ..Default::default()
            }
        })
        .collect()
}

fn same_value(a: &MetaValue, b: &MetaValue) -> bool {
    match (a.data(), b.data()) {
        (MetaValueData::Float(x), MetaValueData::Float(y)) => x.to_bits() == y.to_bits(),
        _ => a == b,
    }
}

fn assert_same(actual: &[PeptideIdentification], expected: &[PeptideIdentification], what: &str) {
    assert_eq!(actual.len(), expected.len(), "{what}: record count");
    for (i, (a, e)) in actual.iter().zip(expected).enumerate() {
        assert_eq!(a.identifier, e.identifier, "{what} record {i}");
        assert_eq!(a.score_type, e.score_type, "{what} record {i} score type");
        assert_eq!(
            a.higher_score_better, e.higher_score_better,
            "{what} record {i} direction"
        );
        assert_eq!(a.hits.len(), e.hits.len(), "{what} record {i} hit count");
        for (j, (ah, eh)) in a.hits.iter().zip(&e.hits).enumerate() {
            assert_eq!(
                ah.score.to_bits(),
                eh.score.to_bits(),
                "{what} record {i} hit {j}: {} vs {}",
                ah.score,
                eh.score
            );
            assert_eq!(ah.charge, eh.charge);
            assert_eq!(ah.sequence.as_str(), eh.sequence.as_str());
            let ak: Vec<_> = ah.metadata.keys().collect();
            let ek: Vec<_> = eh.metadata.keys().collect();
            assert_eq!(ak, ek, "{what} record {i} hit {j} meta keys");
            for (key, value) in &eh.metadata {
                assert!(
                    same_value(&ah.metadata[key], value),
                    "{what} record {i} hit {j} {key}: {:?} vs {value:?}",
                    ah.metadata[key]
                );
            }
        }
    }
}

fn make_cmap(ids: &[PeptideIdentification], features: usize) -> ConsensusMap {
    let mut map = ConsensusMap::default();
    for i in 0..features {
        let mut feature = ConsensusFeature::new();
        feature.base.peptide_identifications = ids[i * 2..i * 2 + 2].to_vec();
        map.features.push(feature);
    }
    map.unassigned_peptide_identifications = ids[features * 2..].to_vec();
    map
}

fn runs(second_charges: &str) -> Vec<ProteinIdentification> {
    let mut a = ProteinIdentification {
        identifier: "runA".into(),
        ..Default::default()
    };
    a.search_parameters.charges = "1:3".into();
    let mut b = ProteinIdentification {
        identifier: "runB".into(),
        ..Default::default()
    };
    b.search_parameters.charges = second_charges.into();
    vec![a, b]
}

/// The outcome record that follows a case: result records, or the source's
/// failure.
enum Outcome {
    Ids(Vec<PeptideIdentification>),
    Map(Vec<Vec<PeptideIdentification>>, Vec<PeptideIdentification>),
    Failed(String),
}

fn read_outcome(lines: &mut Input<'_>) -> Outcome {
    let next = *lines.peek().unwrap();
    let f: Vec<&str> = next.split('\t').collect();
    match f[0] {
        "result" => Outcome::Ids(read_ids(lines, "result")),
        "cmap" => {
            lines.next();
            let n: usize = f[1].parse().unwrap();
            let features = (0..n).map(|_| read_ids(lines, "feature")).collect();
            let unassigned = read_ids(lines, "unassigned");
            Outcome::Map(features, unassigned)
        }
        "error" | "crash" => {
            lines.next();
            Outcome::Failed(next.to_owned())
        }
        other => panic!("unexpected outcome record {other}"),
    }
}

fn fdr(
    all: bool,
    add: bool,
    conservative: bool,
    qvalues: bool,
    treat: bool,
    split: bool,
) -> FalseDiscoveryRate {
    FalseDiscoveryRate {
        use_all_hits: all,
        add_decoy_peptides: add,
        conservative,
        output: if qvalues {
            FdrOutput::QValue
        } else {
            FdrOutput::Fdr
        },
        treat_runs_separately: treat,
        split_charge_variants: split,
        ..FalseDiscoveryRate::default()
    }
}

#[test]
fn basic_members_match_the_release_build() {
    let text = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/fdr_basic_oracle.tsv"),
    )
    .unwrap();
    let mut lines = text.lines().peekable();
    let mut counts = std::collections::BTreeMap::<String, usize>::new();
    let (mut affixes, mut logged) = (Vec::<String>::new(), Vec::<String>::new());
    while let Some(line) = lines.next() {
        let f: Vec<&str> = line.split('\t').collect();
        let flag = |i: usize| f[i] == "1";
        match (f[0], f.get(1).copied()) {
            ("case", Some("basic_selected")) => {
                let calculator = fdr(flag(2), flag(3), flag(4), flag(5), false, false);
                let mut ids = read_ids(&mut lines, "input");
                let result = calculator.apply_basic_selected(
                    &mut ids,
                    flag(6),
                    f[7].parse().unwrap(),
                    f.get(8).copied().unwrap_or(""),
                    f.get(9) == Some(&"1"),
                );
                check(result, &ids, None, read_outcome(&mut lines), line);
            }
            ("case", Some("basic_runs")) => {
                let calculator = fdr(false, flag(4), true, true, flag(2), flag(3));
                let mut ids = read_ids(&mut lines, "input");
                let result = calculator.apply_basic_runs(&runs("2-4"), &mut ids);
                check(result, &ids, None, read_outcome(&mut lines), line);
            }
            ("case", Some(kind @ ("basic_cmap" | "peptide_cmap"))) => {
                let calculator = fdr(flag(5), flag(4), true, true, flag(2), flag(3));
                let ids = read_ids(&mut lines, "input");
                let mut map = make_cmap(&ids, 12);
                map.protein_identifications = runs("+2");
                let result = if kind == "basic_cmap" {
                    calculator.apply_basic_consensus(&mut map, flag(7))
                } else {
                    calculator.apply_basic_peptide_level_consensus(&mut map, flag(7))
                };
                check(result, &[], Some(&map), read_outcome(&mut lines), line);
            }
            ("case", Some("peptide_list")) => {
                let calculator = fdr(flag(5), flag(4), true, true, flag(2), flag(3));
                let mut ids = read_ids(&mut lines, "input");
                let result = calculator.apply_basic_peptide_level_source(&mut ids);
                check(result, &ids, None, read_outcome(&mut lines), line);
            }
            ("decoy", _) => {
                let protein = ProteinIdentification {
                    hits: f[5..]
                        .iter()
                        .map(|acc| ProteinHit {
                            accession: (*acc).into(),
                            ..ProteinHit::default()
                        })
                        .collect(),
                    ..Default::default()
                };
                let found = find_decoy_string(&protein);
                assert_eq!(found.success, f[2] == "1", "{line}");
                assert_eq!(found.name, f[3], "{line}");
                assert_eq!(found.is_prefix, f[4] == "1", "{line}");
            }
            ("case", Some("picked_auto")) => {
                let p: Vec<&str> = lines.next().unwrap().split('\t').collect();
                let mut protein = ProteinIdentification {
                    score_type: "Posterior Probability".into(),
                    higher_score_better: true,
                    hits: p[2..]
                        .iter()
                        .map(|h| {
                            let mut parts = h.splitn(3, ':');
                            let accession = parts.next().unwrap();
                            let score: f64 = parts.next().unwrap().parse().unwrap();
                            let mut hit = ProteinHit {
                                accession: accession.into(),
                                score,
                                ..ProteinHit::default()
                            };
                            hit.metadata.insert(
                                "target_decoy".into(),
                                MetaValue::from(parts.next().unwrap()),
                            );
                            hit
                        })
                        .collect(),
                    ..Default::default()
                };
                let expected: Vec<&str> = lines.next().unwrap().split('\t').collect();
                let used = FalseDiscoveryRate::default()
                    .apply_picked_protein_auto(&mut protein, false)
                    .unwrap();
                affixes.push(format!(
                    "Using {} decoy string '{}'",
                    if used.is_prefix { "prefix" } else { "suffix" },
                    used.name
                ));
                assert_eq!(protein.score_type, expected[1]);
                assert_eq!(protein.higher_score_better, expected[2] == "1");
                let got: Vec<(String, u64)> = protein
                    .hits
                    .iter()
                    .map(|h| (h.accession.clone(), h.score.to_bits()))
                    .collect();
                let want: Vec<(String, u64)> = expected[4..]
                    .iter()
                    .map(|h| {
                        let (a, s) = h.rsplit_once(':').unwrap();
                        (a.to_owned(), s.parse::<f64>().unwrap().to_bits())
                    })
                    .collect();
                assert_eq!(got, want, "{line}");
            }
            // The source's OPENMS_LOG_INFO lines, flushed at the end with
            // the log stream's repeat suppression.
            (other, _) if other.starts_with("Using ") || other.starts_with("<Using ") => {
                logged.push(other.to_owned());
                continue;
            }
            (other, _) => panic!("unexpected oracle record {other}"),
        }
        *counts.entry(f[0..2.min(f.len())].join(" ")).or_default() += 1;
    }
    let total: usize = counts.values().sum();
    assert_eq!(total, 160 + 16 + 64 + 64 + 32 + 11 + 3, "{counts:?}");
    // Three picked cases: the affix each one used, as the source logged it
    // (its second identical line reported as a repeat).
    assert_eq!(
        affixes,
        [
            "Using prefix decoy string 'DECOY_'",
            "Using suffix decoy string '_rev'",
            "Using prefix decoy string 'DECOY_'"
        ]
    );
    let mut logged_sorted = logged.clone();
    logged_sorted.sort();
    assert_eq!(
        logged_sorted,
        [
            "<Using prefix decoy string 'DECOY_'> occurred 2 times",
            "Using prefix decoy string 'DECOY_'",
            "Using suffix decoy string '_rev'"
        ]
    );
}

fn check(
    result: openms::Result<()>,
    ids: &[PeptideIdentification],
    map: Option<&ConsensusMap>,
    outcome: Outcome,
    context: &str,
) {
    match outcome {
        Outcome::Failed(reason) => assert!(
            result.is_err(),
            "{context}: the source failed ({reason}) but the port succeeded"
        ),
        Outcome::Ids(expected) => {
            result.unwrap_or_else(|e| panic!("{context}: {e}"));
            assert_same(ids, &expected, context);
        }
        Outcome::Map(features, unassigned) => {
            result.unwrap_or_else(|e| panic!("{context}: {e}"));
            let map = map.unwrap();
            assert_eq!(map.features.len(), features.len());
            for (i, (feature, expected)) in map.features.iter().zip(&features).enumerate() {
                assert_same(
                    &feature.base.peptide_identifications,
                    expected,
                    &format!("{context} feature {i}"),
                );
            }
            assert_same(
                &map.unassigned_peptide_identifications,
                &unassigned,
                &format!("{context} unassigned"),
            );
        }
    }
}
