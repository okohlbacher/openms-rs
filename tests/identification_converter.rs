// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

// Source literals and independent derivations are distinguished in
// tests/data/identification_converter_provenance.json.
use openms::chemistry::AASequence;
use openms::format::fasta::FASTAEntry;
use openms::identification::graph::*;
use openms::identification::{FlankingResidue, PeakAnnotation, PeptideEvidence, PeptideHit};
use std::collections::{BTreeMap, BTreeSet};

fn fasta(id: &str, description: &str, sequence: &str) -> FASTAEntry {
    FASTAEntry {
        identifier: id.into(),
        description: description.into(),
        sequence: sequence.into(),
    }
}
fn hit() -> PeptideHit {
    PeptideHit::new(7.0, 2, 3, AASequence::parse("PEPTIDE").unwrap()).unwrap()
}
fn parent(graph: &mut IdentificationData, accession: &str) -> ParentId {
    graph
        .register_parent_sequence(ParentSequence::new(accession))
        .unwrap()
}
fn match_map(id: ParentId, matches: impl IntoIterator<Item = ParentMatch>) -> ParentMatches {
    BTreeMap::from([(id, matches.into_iter().collect())])
}

#[test]
fn source_five_record_fasta_projection_preserves_all_literal_fields() {
    let entries: Vec<_> = include_str!("data/identification_converter_fasta.tsv")
        .lines()
        .skip(1)
        .map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            assert_eq!(fields.len(), 4);
            fasta(fields[1], fields[2], fields[3])
        })
        .collect();
    assert_eq!(entries.len(), 5);
    let mut graph = IdentificationData::new().unwrap();
    IdentificationDataConverter::import_sequences(&mut graph, &entries, MoleculeType::Protein, "")
        .unwrap();
    assert_eq!(graph.parent_count(), 5); // literal IdentificationDataConverter_test.cpp:135
    for entry in entries {
        let stored = graph
            .parents()
            .find(|(_, p)| p.accession == entry.identifier)
            .unwrap()
            .1;
        assert_eq!(stored.sequence, entry.sequence);
        assert_eq!(stored.description, entry.description);
        assert_eq!(stored.molecule_type, MoleculeType::Protein);
        assert_eq!(stored.coverage, 0.0);
        assert!(!stored.is_decoy);
    }
}

#[test]
fn source_import_does_not_parse_alphabet_and_decoy_match_is_case_sensitive_substring() {
    let entries = [
        fasta("prefixDECOYsuffix", "", "a[custom]U*"),
        fasta("decoy", "", ""),
        fasta("target", "", "???"),
    ];
    for kind in [
        MoleculeType::Protein,
        MoleculeType::RNA,
        MoleculeType::Compound,
    ] {
        let mut graph = IdentificationData::new().unwrap();
        IdentificationDataConverter::import_sequences(&mut graph, &entries, kind, "DECOY").unwrap();
        let records: Vec<_> = graph.parents().map(|(_, p)| p).collect();
        assert_eq!(records.iter().filter(|p| p.is_decoy).count(), 1);
        for entry in &entries {
            let p = records
                .iter()
                .find(|p| p.accession == entry.identifier)
                .unwrap();
            assert_eq!(p.sequence, entry.sequence);
            assert_eq!(p.molecule_type, kind);
        }
    }
}

#[test]
fn import_uses_existing_parent_merge_and_current_step_provenance() {
    let mut graph = IdentificationData::new().unwrap();
    let mut existing = ParentSequence::new("P1");
    existing.molecule_type = MoleculeType::RNA;
    existing.coverage = 0.25;
    let id = graph.register_parent_sequence(existing).unwrap();
    let sw = graph
        .register_processing_software(ProcessingSoftware::new("Import", "1"))
        .unwrap();
    let step = graph
        .register_processing_step(ProcessingStep::new(sw), None)
        .unwrap();
    graph.set_current_processing_step(step).unwrap();
    IdentificationDataConverter::import_sequences(
        &mut graph,
        &[fasta("P1", "description", "ACGU")],
        MoleculeType::Protein,
        "P",
    )
    .unwrap();
    let result = graph.parent(id).unwrap();
    assert_eq!(result.sequence, "ACGU");
    assert_eq!(result.description, "description");
    assert_eq!(result.molecule_type, MoleculeType::RNA);
    assert_eq!(result.coverage, 0.25);
    assert!(result.is_decoy);
    assert_eq!(result.result.steps_and_scores.len(), 1);
    assert_eq!(
        result.result.steps_and_scores[0].processing_step,
        Some(step)
    );
    IdentificationDataConverter::import_sequences(
        &mut graph,
        &[fasta("P1", "", "")],
        MoleculeType::Protein,
        "",
    )
    .unwrap();
    assert!(graph.parent(id).unwrap().is_decoy);
    assert_eq!(graph.parent(id).unwrap().sequence, "ACGU");
}

#[test]
fn late_import_conflicts_and_limits_roll_back_every_new_parent() {
    let mut graph = IdentificationData::new().unwrap();
    let mut existing = ParentSequence::new("P1");
    existing.sequence = "AAAA".into();
    let id = graph.register_parent_sequence(existing).unwrap();
    let before = graph.parent(id).unwrap().clone();
    let entries = [fasta("new", "", "GG"), fasta("P1", "changed", "BBBB")];
    assert!(
        IdentificationDataConverter::import_sequences(
            &mut graph,
            &entries,
            MoleculeType::Protein,
            "P"
        )
        .is_err()
    );
    assert_eq!(graph.parent_count(), 1);
    assert_eq!(graph.parent(id).unwrap(), &before);
    let mut limited = IdentificationData::with_limits(GraphLimits {
        max_records: 1,
        ..GraphLimits::default()
    })
    .unwrap();
    assert!(
        IdentificationDataConverter::import_sequences(
            &mut limited,
            &[fasta("a", "", "A"), fasta("b", "", "B")],
            MoleculeType::Protein,
            ""
        )
        .is_err()
    );
    assert!(limited.is_empty());
    // Source empty input never examines the otherwise unused pattern.
    let mut zero = IdentificationData::with_limits(GraphLimits {
        max_work: 0,
        max_bytes: 0,
        ..GraphLimits::default()
    })
    .unwrap();
    IdentificationDataConverter::import_sequences(&mut zero, &[], MoleculeType::Protein, "unused")
        .unwrap();
    assert!(zero.is_empty());
}

#[test]
fn export_appends_sorts_all_evidences_and_preserves_duplicates_and_first_flank_bytes() {
    let mut graph = IdentificationData::new().unwrap();
    let z = parent(&mut graph, "Z");
    let a = parent(&mut graph, "A");
    let mut context = ParentMatch::new(Some(2), Some(5));
    context.left_neighbor = "[ignored-tail".into();
    context.right_neighbor = "Rtail".into();
    let mut unknown = ParentMatch::default();
    unknown.left_neighbor.clear();
    unknown.right_neighbor.clear();
    let matches = BTreeMap::from([
        (z, BTreeSet::from([context.clone()])),
        (a, BTreeSet::from([context, unknown])),
    ]);
    let mut value = hit();
    let old = PeptideEvidence {
        protein_accession: "Z".into(),
        start: Some(2),
        end: Some(5),
        aa_before: FlankingResidue::NTerminus,
        aa_after: FlankingResidue::Residue('R'),
    };
    value.evidences.push(old.clone());
    IdentificationDataConverter::export_parent_matches(&graph, &matches, &mut value).unwrap();
    assert_eq!(value.evidences.len(), 4);
    assert_eq!(
        value
            .evidences
            .iter()
            .map(|e| e.protein_accession.as_str())
            .collect::<Vec<_>>(),
        vec!["A", "A", "Z", "Z"]
    );
    assert_eq!(value.evidences[0].start, None); // source unknown -1 sorts before known positions
    assert_eq!(value.evidences[0].aa_before, FlankingResidue::Unknown);
    assert_eq!(value.evidences[1].start, Some(2));
    assert_eq!(value.evidences[2], old);
    assert_eq!(value.evidences[3], old);
    assert_eq!(value.score, 7.0);
    assert_eq!(value.sequence.as_str(), "PEPTIDE");
}

#[test]
fn export_empty_input_still_sorts_existing_evidence_and_ignores_unconsumed_hit_payload() {
    let graph = IdentificationData::new().unwrap();
    let mut value = hit();
    value.evidences = vec![
        PeptideEvidence::new("Z", 2..=2).unwrap(),
        PeptideEvidence::new("A", 1..=1).unwrap(),
    ];
    value.score = f64::NAN;
    value.peak_annotations.push(PeakAnnotation {
        mz: f64::NAN,
        ..Default::default()
    });
    IdentificationDataConverter::export_parent_matches(&graph, &ParentMatches::new(), &mut value)
        .unwrap();
    assert_eq!(value.evidences[0].protein_accession, "A");
    assert!(value.score.is_nan());
    assert!(value.peak_annotations[0].mz.is_nan());
}

#[test]
fn export_checks_unrepresentable_flanks_intervals_and_legacy_position_overflow_atomically() {
    let mut graph = IdentificationData::new().unwrap();
    let id = parent(&mut graph, "P");
    let mut bad_flank = ParentMatch::new(Some(1), Some(2));
    bad_flank.left_neighbor = "m6A".into();
    let mut utf8 = ParentMatch::new(Some(1), Some(2));
    utf8.right_neighbor = "é".into();
    for invalid in [
        bad_flank,
        utf8,
        ParentMatch::new(Some(5), Some(2)),
        ParentMatch::new(Some(i32::MAX as usize + 1), None),
    ] {
        let matches = match_map(id, [ParentMatch::new(Some(0), Some(0)), invalid]);
        let mut value = hit();
        value
            .evidences
            .push(PeptideEvidence::new("existing", 7..=7).unwrap());
        let saved = value.clone();
        assert!(
            IdentificationDataConverter::export_parent_matches(&graph, &matches, &mut value)
                .is_err()
        );
        assert_eq!(value, saved);
    }
    let mut value = hit();
    IdentificationDataConverter::export_parent_matches(
        &graph,
        &match_map(
            id,
            [ParentMatch::new(
                Some(i32::MAX as usize),
                Some(i32::MAX as usize),
            )],
        ),
        &mut value,
    )
    .unwrap();
    assert_eq!(value.evidences[0].end, Some(i32::MAX as usize));
}

#[test]
fn export_rejects_foreign_empty_parent_refs_and_total_output_limit_without_mutation() {
    let mut graph = IdentificationData::with_limits(GraphLimits {
        max_edges: 1,
        ..GraphLimits::default()
    })
    .unwrap();
    let id = parent(&mut graph, "P");
    let mut foreign = IdentificationData::new().unwrap();
    let foreign = parent(&mut foreign, "P");
    let mut value = hit();
    value
        .evidences
        .push(PeptideEvidence::new("existing", 0..=0).unwrap());
    let saved = value.clone();
    assert!(
        IdentificationDataConverter::export_parent_matches(
            &graph,
            &match_map(foreign, []),
            &mut value
        )
        .is_err()
    );
    assert_eq!(value, saved);
    assert!(
        IdentificationDataConverter::export_parent_matches(
            &graph,
            &match_map(id, [ParentMatch::new(Some(0), Some(0))]),
            &mut value
        )
        .is_err()
    );
    assert_eq!(value, saved);
}
