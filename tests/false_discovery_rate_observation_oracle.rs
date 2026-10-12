// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Executed differential (tier 1) for
//! `FalseDiscoveryRate::applyToObservationMatches`.
//!
//! `tests/data/fdr_observation_oracle.tsv` is the output of
//! `../oracle/I2-fdr/driver3.cpp` linked against the Linux x86_64 Release
//! `libOpenMS.so` on ibminode06 (two runs, identical bytes): six synthetic
//! identification graphs (10 proteins of which 4 decoys, 16 peptides, 40
//! observations with up to three peptide matches and some compound matches,
//! a tenth of the matches unscored), each under all eight combinations of
//! `use_all_hits`, `add_decoy_peptides` and q-value/FDR output. Contract: the
//! same registered score type and bit-identical values on the same matches.
//!
//! Scores are distinct within one observation (ties across observations
//! stay): the source breaks a best-match tie by molecule heap address, and a
//! first fixture with such ties gave a different winner in the first of eight
//! copies of the same graph, in two of six graphs.

use openms::analysis::false_discovery_rate::{FalseDiscoveryRate, FdrOutput};
use openms::chemistry::AASequence;
use openms::identification::graph::{
    IdentificationData, IdentifiedCompound, IdentifiedMolecule, IdentifiedPeptide, InputFile,
    MoleculeType, Observation, ObservationMatch, ParentMatch, ParentSequence, ScoreType,
    ScoreTypeId,
};
use openms::metadata::CVTerm;
use std::collections::BTreeMap;
use std::path::Path;

type Key = (String, String, String);

/// Rebuild one graph from its `protein`/`peptide`/`match` records, in the
/// source's registration order.
fn build(records: &[Vec<&str>], higher: bool) -> (IdentificationData, ScoreTypeId) {
    let mut data = IdentificationData::new().unwrap();
    let file = data
        .register_input_file(InputFile::new("test.mzML"))
        .unwrap();
    let score = data
        .register_score_type(ScoreType {
            cv_term: CVTerm::new("", "score", ""),
            higher_better: higher,
            ..ScoreType::default()
        })
        .unwrap();
    let (mut proteins, mut peptides, mut observations) = (Vec::new(), Vec::new(), BTreeMap::new());
    // The compound is registered where the source registers it: after the
    // peptides, before the first observation.
    let mut compound = None;
    for g in records {
        match g[0] {
            "protein" => {
                let parent = ParentSequence {
                    molecule_type: MoleculeType::Protein,
                    is_decoy: g[2] == "1",
                    ..ParentSequence::new(format!("P{}", g[1]))
                };
                proteins.push(data.register_parent_sequence(parent).unwrap());
            }
            "peptide" => {
                let mut peptide = IdentifiedPeptide::new(AASequence::parse(g[2]).unwrap());
                for parent in [g[3], g[4]] {
                    let index: i64 = parent.parse().unwrap();
                    if index >= 0 {
                        peptide
                            .parent_matches
                            .entry(proteins[index as usize])
                            .or_default()
                            .insert(ParentMatch::default());
                    }
                }
                peptides.push(data.register_identified_peptide(peptide).unwrap());
            }
            "match" => {
                if compound.is_none() {
                    compound = Some(
                        data.register_identified_compound(IdentifiedCompound::new("glucose"))
                            .unwrap(),
                    );
                }
                let o: usize = g[1].parse().unwrap();
                let observation = match observations.get(&o) {
                    Some(&id) => id,
                    None => {
                        let id = data
                            .register_observation(Observation::new(format!("spectrum_{o}"), file))
                            .unwrap();
                        observations.insert(o, id);
                        id
                    }
                };
                let (molecule, charge): (IdentifiedMolecule, i32) = if g[2] == "p" {
                    (peptides[g[3].parse::<usize>().unwrap()].into(), 2)
                } else {
                    (compound.unwrap().into(), 1)
                };
                let mut m = ObservationMatch::new(molecule, observation);
                m.charge = charge;
                if g[4] != "-" {
                    m.result
                        .add_score(score, g[4].parse().unwrap(), None)
                        .unwrap();
                }
                data.register_observation_match(m).unwrap();
            }
            other => panic!("unexpected graph record {other}"),
        }
    }
    (data, score)
}

#[test]
fn observation_match_fdr_matches_the_release_build() {
    let text = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/fdr_observation_oracle.tsv"),
    )
    .unwrap();
    let mut lines = text.lines().peekable();
    let mut runs_checked = 0;
    while let Some(line) = lines.next() {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f[0], "graph", "{line}");
        let higher = f[2] == "1";
        let mut records = Vec::new();
        while let Some(next) = lines.peek() {
            let g: Vec<&str> = next.split('\t').collect();
            if !matches!(g[0], "protein" | "peptide" | "match") {
                break;
            }
            records.push(g);
            lines.next();
        }
        while let Some(next) = lines.peek() {
            let r: Vec<&str> = next.split('\t').collect();
            if r[0] != "run" {
                break;
            }
            lines.next();
            let calculator = FalseDiscoveryRate {
                use_all_hits: r[1] == "1",
                add_decoy_peptides: r[2] == "1",
                output: if r[3] == "1" {
                    FdrOutput::QValue
                } else {
                    FdrOutput::Fdr
                },
                ..FalseDiscoveryRate::default()
            };
            let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
            assert_eq!(head[0], "score_type");
            let (mut data, score) = build(&records, higher);
            let fdr_type = calculator
                .apply_to_observation_matches(&mut data, score)
                .unwrap();
            let registered = data.score_type(fdr_type).unwrap();
            assert_eq!(registered.cv_term.accession, head[1]);
            assert_eq!(registered.cv_term.name, head[2]);
            assert_eq!(registered.higher_better, head[3] == "1");
            let mut expected: BTreeMap<Key, Option<f64>> = BTreeMap::new();
            while let Some(next) = lines.peek() {
                let e: Vec<&str> = next.split('\t').collect();
                if e[0] != "fdr" {
                    break;
                }
                lines.next();
                let value = (e[4] != "-").then(|| e[4].parse::<f64>().unwrap());
                expected.insert((e[1].into(), e[2].into(), e[3].into()), value);
            }
            let mut got: BTreeMap<Key, Option<f64>> = BTreeMap::new();
            for (_, record) in data.observation_matches() {
                let data_id = data
                    .observation(record.observation)
                    .unwrap()
                    .data_id
                    .clone();
                let (kind, name) = match record.identified_molecule {
                    IdentifiedMolecule::Peptide(p) => (
                        "p".to_owned(),
                        data.peptide(p).unwrap().sequence.as_str().to_owned(),
                    ),
                    _ => ("c".to_owned(), "0".to_owned()),
                };
                got.insert((data_id, kind, name), record.result.score(fdr_type));
            }
            assert_eq!(got.len(), expected.len(), "{line} run {r:?}");
            for (key, want) in &expected {
                let have = got[key];
                assert_eq!(
                    have.map(f64::to_bits),
                    want.map(f64::to_bits),
                    "{line} run {r:?} match {key:?}: {have:?} vs {want:?}"
                );
            }
            runs_checked += 1;
        }
    }
    assert_eq!(runs_checked, 6 * 8);
}
