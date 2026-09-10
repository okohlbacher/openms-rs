// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Complete-value annotation identity across caller-owned registries. Source
//! ProteinIdentification collects (position, ResidueModification) values, while
//! native conflict resolution must not merge unequal annotated sequences.

use openms::analysis::id_conflict_resolver::{
    ResolutionMethod, reduce_to_one_per_spectrum, resolve_identifications,
};
use openms::analysis::id_filter::{DuplicatePolicy, remove_duplicate_peptide_hits};
use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationRecord, ModificationsDB, ResidueModification,
};
use openms::identification::{
    PeptideEvidence, PeptideHit, PeptideIdentification, ProteinHit, ProteinIdentification,
};
use std::collections::BTreeSet;

fn peptide(formula: &str, accession: &str) -> AASequence {
    let formula = EmpiricalFormula::parse(formula).unwrap();
    let record = ResidueModification::from_record(ModificationRecord {
        name: "Lab".into(),
        full_name: "Laboratory annotation".into(),
        origin: Some('M'),
        diff_mono_mass: formula.mono_mass(),
        diff_average_mass: formula.average_mass(),
        diff_formula: formula,
        obo_accession: Some(accession.into()),
        ..Default::default()
    })
    .unwrap();
    let database = ModificationsDB::from_records(vec![record]).unwrap();
    // Return the sequence after dropping the caller-owned database.
    AASequence::parse_with_registry("AM(Lab)A", &database).unwrap()
}
fn id(hits: Vec<(AASequence, f64)>) -> PeptideIdentification {
    let mut id = PeptideIdentification {
        identifier: "run".into(),
        score_type: "synthetic".into(),
        higher_score_better: true,
        hits: hits
            .into_iter()
            .map(|(peptide, score)| PeptideHit {
                sequence: peptide,
                score,
                charge: 2,
                evidences: vec![PeptideEvidence::new("P", 0..=2).unwrap()],
                ..Default::default()
            })
            .collect(),
        ..Default::default()
    };
    id.set_spectrum_reference("scan=1");
    id
}
fn proteins() -> ProteinIdentification {
    ProteinIdentification {
        hits: vec![ProteinHit::new(0.0, 0, "P", "AMA").unwrap()],
        ..Default::default()
    }
}

#[test]
fn owned_hit_identity_agrees_with_complete_sequence_and_charge_equality() {
    let first = peptide("O", "MOD:90001");
    let different_mass = peptide("O2", "MOD:90001");
    let different_vocabulary = peptide("O", "MOD:90002");
    assert_eq!(first.to_string(), different_mass.to_string());
    assert_eq!(first.to_string(), different_vocabulary.to_string());
    let mut hits = id(vec![
        (first.clone(), 4.0),
        (different_mass, 3.0),
        (different_vocabulary, 2.0),
        (first, 1.0),
    ])
    .hits;
    let mut changed_charge = hits[0].clone();
    changed_charge.charge = 3;
    hits.push(changed_charge);
    for first in &hits {
        for second in &hits {
            assert_eq!(
                first.identity_key() == second.identity_key(),
                first.same_sequence_and_charge(second)
            );
        }
    }
    let keys: BTreeSet<_> = hits.iter().map(PeptideHit::identity_key).collect();
    assert_eq!(keys.len(), 4);
    drop(hits);
    assert_eq!(keys.iter().filter(|(_, charge)| *charge == 2).count(), 3);
}

#[test]
fn sequence_duplicate_policy_keeps_distinct_chemistry_and_first_equal_value() {
    let first = peptide("O", "MOD:90001");
    let second = peptide("O2", "MOD:90001");
    let different_vocabulary = peptide("O", "MOD:90002");
    let mut ids = vec![id(vec![
        (first.clone(), 1.0),
        (second.clone(), 2.0),
        (different_vocabulary.clone(), 3.0),
        (peptide("O", "MOD:90001"), 4.0),
    ])];
    // Sequence-only filtering ignores charge as well as score, but it must not
    // flatten owned chemistry to the displayed name.
    ids[0].hits[3].charge = 3;
    remove_duplicate_peptide_hits(&mut ids, DuplicatePolicy::Sequence).unwrap();
    assert_eq!(ids[0].hits.len(), 3);
    for (hit, expected) in ids[0]
        .hits
        .iter()
        .zip([first, second, different_vocabulary])
    {
        assert_eq!(hit.sequence, expected);
    }
    assert_eq!(ids[0].hits[0].score, 1.0);
    assert_eq!(ids[0].hits[0].charge, 2);
}

#[test]
fn protein_observations_keep_distinct_chemistry_and_vocabularies_at_one_position() {
    let first = peptide("O", "MOD:90001");
    let different_mass = peptide("O2", "MOD:90001");
    let different_vocabulary = peptide("O", "MOD:90002");
    let identical_copy = peptide("O", "MOD:90001");
    assert_eq!(first, identical_copy);
    assert_ne!(first, different_mass);
    let mut hits = vec![
        (first.clone(), 4.0),
        (different_mass, 3.0),
        (different_vocabulary, 2.0),
        (identical_copy, 1.0),
    ];
    let mut forward = proteins();
    forward
        .compute_modifications(&[id(hits.clone())], &BTreeSet::new())
        .unwrap();
    let modifications = &forward.hits[0].modifications;
    assert_eq!(modifications.len(), 3);
    assert!(
        modifications
            .iter()
            .all(|m| m.position == 1 && m.modification.full_id() == "Lab (M)")
    );
    let observed: BTreeSet<_> = modifications
        .iter()
        .map(|m| {
            let record = m.modification.known().unwrap();
            (
                record.obo_accession().unwrap().to_owned(),
                record.diff_formula().to_string(),
            )
        })
        .collect();
    assert_eq!(
        observed,
        BTreeSet::from([
            (
                "MOD:90001".into(),
                EmpiricalFormula::parse("O").unwrap().to_string()
            ),
            (
                "MOD:90001".into(),
                EmpiricalFormula::parse("O2").unwrap().to_string()
            ),
            (
                "MOD:90002".into(),
                EmpiricalFormula::parse("O").unwrap().to_string()
            ),
        ])
    );
    hits.reverse();
    let mut reverse = proteins();
    reverse
        .compute_modifications(&[id(hits)], &BTreeSet::new())
        .unwrap();
    assert_eq!(
        reverse, forward,
        "observations must not depend on input order"
    );
    let mut skipped = proteins();
    skipped
        .compute_modifications(
            &[id(vec![(first, 1.0)])],
            &BTreeSet::from(["Lab (M)".into()]),
        )
        .unwrap();
    assert!(skipped.hits[0].modifications.is_empty());
}

#[test]
fn complete_observation_collection_remains_atomic_on_a_later_invalid_position() {
    let first = peptide("O", "MOD:90001");
    let second = peptide("O2", "MOD:90001");
    let mut proteins = proteins();
    proteins
        .compute_modifications(&[id(vec![(first.clone(), 1.0)])], &BTreeSet::new())
        .unwrap();
    let before = proteins.clone();
    let mut ids = vec![id(vec![(first, 2.0), (second, 1.0)])];
    ids[0].hits[1].evidences = vec![PeptideEvidence::new("P", 3..=5).unwrap()];
    assert!(
        proteins
            .compute_modifications(&ids, &BTreeSet::new())
            .is_err()
    );
    assert_eq!(proteins, before);
}

#[test]
fn spectrum_reduction_deduplicates_equal_values_without_merging_distinct_chemistry() {
    let first = peptide("O", "MOD:90001");
    let second = peptide("O2", "MOD:90001");
    let duplicate = peptide("O", "MOD:90001");
    let mut ids = vec![
        id(vec![(first, 1.0)]),
        id(vec![(second.clone(), 8.0)]),
        id(vec![(duplicate.clone(), 9.0)]),
    ];
    let report = reduce_to_one_per_spectrum(&mut ids).unwrap();
    assert_eq!(report.removed, 1);
    assert_eq!(ids.len(), 2);
    assert_eq!(ids[0].hits[0].sequence, second);
    assert_eq!(ids[1].hits[0].sequence, duplicate);
    assert_eq!(ids[1].hits[0].score, 9.0);
}

#[test]
fn rank_aggregation_counts_each_chemical_peptidoform_separately() {
    let first = peptide("O", "MOD:90001");
    let second = peptide("O2", "MOD:90001");
    let mut ids = vec![
        id(vec![(first, 100.0), (second.clone(), 99.0)]),
        id(vec![(second.clone(), 80.0)]),
        id(vec![(second.clone(), 70.0)]),
    ];
    let mut removed = Vec::new();
    resolve_identifications(
        &mut ids,
        &mut removed,
        17,
        ResolutionMethod::RankAggregation,
    )
    .unwrap();
    // With maximum rank-list length2 and3 lists, the second form has rank sum1
    // versus4 (including missing-hit penalties) for the first, so it wins.
    assert_eq!(ids.len(), 1);
    assert_eq!(ids[0].hits.len(), 1);
    assert_eq!(ids[0].hits[0].sequence, second);
    assert_eq!(ids[0].hits[0].score, 99.0);
    assert_eq!(removed.len(), 2);
}
