// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Identification operations preserve unresolved sequences and owned mass tags.

use openms::Error;
use openms::analysis::id_conflict_resolver::{
    ResolutionMethod, reduce_to_one_per_spectrum, resolve_identifications,
};
use openms::analysis::id_filter::{
    MatchAction, filter_peptides_by_length, filter_peptides_by_modifications,
    filter_peptides_by_mz_error, filter_peptides_by_sequences,
};
use openms::chemistry::AASequence;
use openms::comparison::Tolerance;
use openms::identification::{
    PeptideEvidence, PeptideHit, PeptideIdentification, ProteinHit, ProteinIdentification,
};
use std::collections::BTreeSet;

fn hit(sequence: &str, score: f64) -> PeptideHit {
    PeptideHit::new(score, 7, 2, AASequence::parse(sequence).unwrap()).unwrap()
}
fn id(hits: Vec<PeptideHit>) -> PeptideIdentification {
    let mut result = PeptideIdentification {
        identifier: "run".into(),
        score_type: "synthetic".into(),
        hits,
        ..Default::default()
    };
    result.set_spectrum_reference("scan=1");
    result
}
fn names(values: &[&str]) -> BTreeSet<String> {
    values.iter().map(|s| (*s).into()).collect()
}

#[test]
fn unresolved_residues_support_non_mass_filters_without_becoming_modifications() {
    let mut ids = vec![id(vec![hit("BZX", 3.0), hit("AK", 2.0)])];
    filter_peptides_by_length(&mut ids, 3, Some(3)).unwrap();
    assert_eq!(ids[0].hits[0].sequence.as_str(), "BZX");
    filter_peptides_by_sequences(&mut ids, &names(&["BZX"]), false, MatchAction::Keep).unwrap();
    assert_eq!(ids[0].hits.len(), 1);
    assert!(!ids[0].hits[0].sequence.is_modified());
    filter_peptides_by_modifications(&mut ids, &BTreeSet::new(), MatchAction::Keep).unwrap();
    assert!(ids[0].hits.is_empty());
}

#[test]
fn modification_filters_match_known_and_anonymous_full_ids_including_termini() {
    // Exact unknown full IDs follow ResidueModification::createUnknownFromMassString;
    // HasMatchingModification compares full IDs without requiring a registry lookup.
    let original = vec![id(vec![
        hit("C(Carbamidomethyl)K", 3.0),
        hit("X[123.456789]K", 2.0),
        hit(".[+12.3456789]AK", 1.0),
        hit("AK.[+34.5678912]", 0.5),
        hit("BZX", 0.1),
    ])];
    for (full_id, expected) in [
        ("Carbamidomethyl (C)", 0),
        ("X[123.456789]", 1),
        (".n[+12.3456789]", 2),
        (".c[+34.5678912]", 3),
    ] {
        let mut ids = original.clone();
        filter_peptides_by_modifications(&mut ids, &names(&[full_id]), MatchAction::Keep).unwrap();
        assert_eq!(ids[0].hits, [original[0].hits[expected].clone()]);
    }
    let mut ids = original.clone();
    filter_peptides_by_modifications(&mut ids, &BTreeSet::new(), MatchAction::Keep).unwrap();
    assert_eq!(ids[0].hits.len(), 4);
    let mut ids = original;
    filter_peptides_by_modifications(
        &mut ids,
        &names(&["Carbamidomethyl", "123.456789"]),
        MatchAction::Keep,
    )
    .unwrap();
    assert!(ids[0].hits.is_empty()); // short names are not full modification IDs
}

#[test]
fn precursor_filter_rejects_later_unresolved_mass_atomically() {
    let mut first = id(vec![hit("AA", 2.0)]);
    first.mz = Some(1000.0); // would be removed
    let mut second = id(vec![hit("BZX", 1.0)]);
    second.mz = Some(1000.0);
    let mut ids = vec![first, second];
    let before = ids.clone();
    assert!(matches!(
        filter_peptides_by_mz_error(&mut ids, Tolerance::Absolute(0.01)),
        Err(Error::Unsupported(_))
    ));
    assert_eq!(ids, before);
}

#[test]
fn absolute_mass_tag_filters_with_known_monoisotopic_mass_without_formula() {
    let sequence = AASequence::parse("AX[123.456789]K").unwrap();
    assert!(sequence.formula().is_err());
    assert!(sequence.average_mass().is_err());
    // (A_internal + 123.456789 + K_internal + H2O + 2 protons) / 2,
    // using source residue masses rounded to six decimal places.
    let observed = 171.306_992_3;
    assert!((sequence.mz(2).unwrap() - observed).abs() < 1e-6);
    let mut record = id(vec![PeptideHit::new(1.0, 7, 2, sequence.clone()).unwrap()]);
    record.mz = Some(observed);
    let mut ids = vec![record];
    filter_peptides_by_mz_error(&mut ids, Tolerance::Ppm(0.1)).unwrap();
    assert_eq!(ids[0].hits[0].sequence, sequence);
    ids[0].mz = Some(observed + 1.0);
    filter_peptides_by_mz_error(&mut ids, Tolerance::Absolute(0.01)).unwrap();
    assert!(ids[0].hits.is_empty());
}

#[test]
fn protein_modifications_own_mass_tags_preserve_known_records_and_honor_skip_sets() {
    let mut peptide = hit(
        ".[+12.3456789]AX[123.456789]C(Carbamidomethyl).[+34.5678912]",
        1.0,
    );
    peptide.evidences = vec![PeptideEvidence::new("P", 2..=4).unwrap(); 2];
    let mut proteins = ProteinIdentification {
        hits: vec![ProteinHit::new(0.0, 0, "P", "QQAXCQQ").unwrap()],
        ..Default::default()
    };
    let ids = vec![id(vec![peptide])];
    proteins
        .compute_modifications(&ids, &BTreeSet::new())
        .unwrap();
    let observed: Vec<_> = proteins.hits[0]
        .modifications
        .iter()
        .map(|m| (m.position, m.modification.full_id()))
        .collect();
    assert_eq!(
        observed,
        [
            (2, ".n[+12.3456789]"),
            (3, "X[123.456789]"),
            (4, ".c[+34.5678912]"),
            (4, "Carbamidomethyl (C)")
        ]
    );
    assert!(
        proteins.hits[0].modifications[1]
            .modification
            .known()
            .is_none()
    );
    assert_eq!(
        proteins.hits[0].modifications[1]
            .modification
            .mass_tag()
            .unwrap()
            .mass(),
        123.456789
    );
    assert!(
        proteins.hits[0].modifications[3]
            .modification
            .known()
            .is_some()
    );
    let cloned = proteins.clone();
    proteins
        .compute_modifications(
            &ids,
            &names(&[".n[+12.3456789]", "X[123.456789]", "Carbamidomethyl"]),
        )
        .unwrap();
    assert_eq!(proteins.hits[0].modifications.len(), 1);
    assert_eq!(
        proteins.hits[0].modifications[0].modification.full_id(),
        ".c[+34.5678912]"
    );
    drop(ids);
    drop(proteins);
    cloned.validate().unwrap();
    assert_eq!(cloned.hits[0].modifications.len(), 4); // no annotation borrowed from dropped PSMs
}

#[test]
fn invalid_protein_position_keeps_owned_modifications_unchanged() {
    let mut valid = hit("X[123.456789]", 1.0);
    valid.evidences = vec![PeptideEvidence::new("P", 0..=0).unwrap()];
    let mut proteins = ProteinIdentification {
        hits: vec![ProteinHit::new(0.0, 0, "P", "X").unwrap()],
        ..Default::default()
    };
    proteins
        .compute_modifications(&[id(vec![valid])], &BTreeSet::new())
        .unwrap();
    let before = proteins.clone();
    let mut invalid = hit("AX[124.567891]", 1.0);
    invalid.evidences = vec![PeptideEvidence::new("P", 0..=0).unwrap()];
    assert!(
        proteins
            .compute_modifications(&[id(vec![invalid])], &BTreeSet::new())
            .is_err()
    );
    assert_eq!(proteins, before);
    let mut malformed = before.clone();
    malformed.hits[0].modifications[0].position = 1;
    assert!(malformed.validate().is_err());
}

#[test]
fn spectrum_conflicts_distinguish_owned_mass_tags_and_preserve_unresolved_sequences() {
    let mut ids = vec![
        id(vec![hit("AX[123.456789]", 1.0)]),
        id(vec![hit("AX[124.567891]", 2.0)]),
        id(vec![hit("AX[123.456789]", 3.0)]),
        id(vec![hit("BZX", 4.0)]),
    ];
    let report = reduce_to_one_per_spectrum(&mut ids).unwrap();
    assert_eq!(report.removed, 1);
    assert_eq!(report.multiply_identified_spectra, 1);
    assert_eq!(ids.len(), 3);
    assert_eq!(ids[1].hits[0].score, 3.0);
    assert_eq!(
        ids[0].hits[0]
            .sequence
            .residue_modification(1)
            .unwrap()
            .unwrap()
            .full_id(),
        "X[124.567891]"
    );
    assert_eq!(ids[2].hits[0].sequence.as_str(), "BZX");
}

#[test]
fn rank_resolution_ties_use_distinct_deterministic_mass_tag_keys() {
    let first = "AX[123.456789]";
    let second = "AX[124.567891]";
    let mut ids = vec![
        id(vec![hit(first, 3.0), hit(second, 2.0)]),
        id(vec![hit(second, 3.0), hit(first, 2.0)]),
    ];
    let mut rejected = Vec::new();
    resolve_identifications(
        &mut ids,
        &mut rejected,
        42,
        ResolutionMethod::RankAggregation,
    )
    .unwrap();
    assert_eq!(ids.len(), 1);
    assert_eq!(ids[0].hits.len(), 1);
    assert_eq!(ids[0].hits[0].sequence, AASequence::parse(first).unwrap());
    assert!(!rejected.is_empty());
}
