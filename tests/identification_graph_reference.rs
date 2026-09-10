// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Source literals and independent semantic cases from the retained 6bfc SDK.
//! Provenance: tests/data/identification_graph_provenance.json. No C++ execution.

use openms::chemistry::{AASequence, NASequence, RNaseDigestion};
use openms::identification::graph::*;
use std::collections::{BTreeMap, BTreeSet};

fn peptide(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}
fn oligo(text: &str) -> NASequence {
    NASequence::parse(text).unwrap()
}
fn matched(start: usize, end: usize) -> ParentMatch {
    ParentMatch::new(Some(start), Some(end))
}
fn parent(accession: &str, sequence: &str, molecule_type: MoleculeType) -> ParentSequence {
    let mut value = ParentSequence::new(accession);
    value.sequence = sequence.into();
    value.molecule_type = molecule_type;
    value
}

#[test]
fn source_registration_literals_keep_distinct_semantic_keys() {
    // IdentificationData_test.cpp:70-183. Explicit None time replaces source now().
    let mut graph = IdentificationData::new().unwrap();
    let file = InputFile::new("test.mzML");
    let file_id = graph.register_input_file(file.clone()).unwrap();
    assert_eq!(graph.register_input_file(file).unwrap(), file_id);
    assert_eq!(graph.input_file_count(), 1);
    let software = ProcessingSoftware::new("Tool", "1.0");
    let software_id = graph
        .register_processing_software(software.clone())
        .unwrap();
    assert_eq!(
        graph.register_processing_software(software).unwrap(),
        software_id
    );
    assert_eq!(graph.processing_software_count(), 1);
    let parameters = DBSearchParam {
        database: "test-db.fasta".into(),
        precursor_mass_tolerance: 1.0,
        fragment_mass_tolerance: 2.0,
        ..DBSearchParam::default()
    };
    let search = graph.register_db_search_param(parameters.clone()).unwrap();
    assert_eq!(graph.register_db_search_param(parameters).unwrap(), search);
    assert_eq!(graph.db_search_param_count(), 1);
    let mut first = ProcessingStep::new(software_id);
    first.input_files.push(file_id);
    let first_id = graph.register_processing_step(first.clone(), None).unwrap();
    assert_eq!(
        graph.register_processing_step(first, None).unwrap(),
        first_id
    );
    let second = ProcessingStep::new(software_id);
    let second_id = graph
        .register_processing_step(second.clone(), Some(search))
        .unwrap();
    assert_eq!(
        graph
            .register_processing_step(second, Some(search))
            .unwrap(),
        second_id
    );
    assert_eq!(graph.processing_step_count(), 2);
    assert_eq!(graph.search_param_for_step(first_id).unwrap(), None);
    assert_eq!(
        graph.search_param_for_step(second_id).unwrap(),
        Some(search)
    );
    let score = ScoreType::new("test_score", true);
    assert!(score.cv_term.accession.is_empty());
    let score_id = graph.register_score_type(score.clone()).unwrap();
    assert_eq!(graph.register_score_type(score).unwrap(), score_id);
    assert_eq!(graph.score_type_count(), 1);
    assert!(
        graph
            .register_score_type(ScoreType::new("test_score", false))
            .is_err()
    );
    assert!(graph.score_type(score_id).unwrap().higher_better);
}

#[test]
fn source_parentless_nodes_and_literal_protein_coverages() {
    // IdentificationData_test.cpp:210-332,499-510; integer endpoints are inclusive.
    let mut graph = IdentificationData::new().unwrap();
    assert!(
        graph
            .register_parent_sequence(ParentSequence::new(""))
            .is_err()
    );
    let protein = graph
        .register_parent_sequence(parent("protein_1", "TESTPEPTIDEAAA", MoleculeType::Protein))
        .unwrap();
    let rna = parent("rna_1", "", MoleculeType::RNA);
    let rna_id = graph.register_parent_sequence(rna.clone()).unwrap();
    assert_eq!(graph.register_parent_sequence(rna).unwrap(), rna_id);
    assert_eq!(graph.parent_count(), 2);
    assert!(
        graph
            .register_identified_peptide(IdentifiedPeptide::new(peptide("")))
            .is_err()
    );
    assert!(
        graph
            .register_identified_oligo(IdentifiedOligo::new(oligo("")))
            .is_err()
    );
    graph
        .register_identified_peptide(IdentifiedPeptide::new(peptide("TEST")))
        .unwrap();
    graph
        .register_identified_oligo(IdentifiedOligo::new(oligo("ACGU")))
        .unwrap();
    let mut hit = IdentifiedPeptide::new(peptide("PEPTIDE"));
    hit.parent_matches
        .insert(protein, BTreeSet::from([matched(4, 10)]));
    let hit_id = graph.register_identified_peptide(hit.clone()).unwrap();
    assert_eq!(
        graph.register_identified_peptide(hit.clone()).unwrap(),
        hit_id
    );
    let mut rna_hit = IdentifiedOligo::new(oligo("UGCA"));
    rna_hit.parent_matches.insert(rna_id, BTreeSet::new());
    graph.register_identified_oligo(rna_hit.clone()).unwrap();
    assert_eq!((graph.peptide_count(), graph.oligo_count()), (2, 2));
    hit.parent_matches.insert(rna_id, BTreeSet::new());
    rna_hit.parent_matches.insert(protein, BTreeSet::new());
    assert!(graph.register_identified_peptide(hit).is_err());
    assert!(graph.register_identified_oligo(rna_hit).is_err());
    graph.calculate_coverages(false).unwrap();
    assert_eq!(graph.parent(protein).unwrap().coverage, 0.5);
    let mut overlapping = IdentifiedPeptide::new(peptide("TESTPEP"));
    overlapping
        .parent_matches
        .insert(protein, BTreeSet::from([matched(0, 6)]));
    graph.register_identified_peptide(overlapping).unwrap();
    graph.calculate_coverages(false).unwrap();
    assert_eq!(graph.parent(protein).unwrap().coverage, 11.0 / 14.0);
}

#[test]
fn source_parent_and_match_merges_keep_first_nonkey_payload() {
    // ParentSequence.h and IdentifiedSequence.h define these finite merge results.
    let mut graph = IdentificationData::new().unwrap();
    let mut original = parent("parent", "", MoleculeType::Protein);
    original.coverage = 0.25;
    original
        .result
        .metadata
        .insert("origin".into(), "first".into());
    let id = graph.register_parent_sequence(original).unwrap();
    let mut update = parent("parent", "PEPTIDE", MoleculeType::RNA);
    update.description = "filled".into();
    update.coverage = 0.75;
    update.is_decoy = true;
    update
        .result
        .metadata
        .insert("origin".into(), "second".into());
    assert_eq!(graph.register_parent_sequence(update).unwrap(), id);
    let merged = graph.parent(id).unwrap();
    assert_eq!(merged.molecule_type, MoleculeType::Protein);
    assert_eq!(merged.coverage, 0.25);
    assert!(merged.is_decoy);
    assert_eq!(merged.sequence, "PEPTIDE");
    assert_eq!(merged.description, "filled");
    assert_eq!(merged.result.metadata["origin"], "second".into());
    let before = merged.clone();
    let mut conflict = parent("parent", "EDIT", MoleculeType::Protein);
    conflict
        .result
        .metadata
        .insert("origin".into(), "must not commit".into());
    assert!(graph.register_parent_sequence(conflict).is_err());
    assert_eq!(graph.parent(id).unwrap(), &before);

    let mut first = matched(0, 6);
    first.left_neighbor = "[".into();
    first.right_neighbor = "]".into();
    first.metadata.insert("source".into(), "first".into());
    let mut second = first.clone();
    second.left_neighbor = "K".into();
    second.metadata.insert("source".into(), "second".into());
    let mut hit = IdentifiedPeptide::new(peptide("PEPTIDE"));
    hit.parent_matches
        .insert(id, BTreeSet::from([first.clone()]));
    let hit_id = graph.register_identified_peptide(hit.clone()).unwrap();
    hit.parent_matches.insert(id, BTreeSet::from([second]));
    assert_eq!(graph.register_identified_peptide(hit).unwrap(), hit_id);
    let kept = graph.peptide(hit_id).unwrap().parent_matches[&id]
        .first()
        .unwrap();
    assert_eq!(kept.left_neighbor, "[");
    assert_eq!(kept.metadata["source"], "first".into());
    assert!(
        graph
            .peptide(hit_id)
            .unwrap()
            .all_parents_are_decoys(|p| Ok(graph.parent(p)?.is_decoy))
            .unwrap()
    );
}

#[test]
fn source_score_history_updates_do_not_reorder_and_priorities_can_repeat() {
    // AppliedProcessingStep.h:58-96; ScoredProcessingResult.h:43-70,153-190.
    let mut graph = IdentificationData::new().unwrap();
    let a = graph
        .register_score_type(ScoreType::new("a", true))
        .unwrap();
    let b = graph
        .register_score_type(ScoreType::new("b", false))
        .unwrap();
    let mut software = ProcessingSoftware::new("Tool", "1.0");
    software.assigned_scores = vec![b, b, a];
    let software_id = graph
        .register_processing_software(software.clone())
        .unwrap();
    software.assigned_scores = vec![a];
    assert_eq!(
        graph.register_processing_software(software).unwrap(),
        software_id
    );
    assert_eq!(
        graph
            .processing_software(software_id)
            .unwrap()
            .assigned_scores,
        [b, b, a]
    );
    let first = graph
        .register_processing_step(ProcessingStep::new(software_id), None)
        .unwrap();
    let file = graph
        .register_input_file(InputFile::new("test.mzML"))
        .unwrap();
    let mut next = ProcessingStep::new(software_id);
    next.input_files.push(file);
    let second = graph.register_processing_step(next, None).unwrap();
    let mut result = ScoredProcessingResult::default();
    result.add_score(a, 1.0, Some(first)).unwrap();
    result.add_score(a, 2.0, Some(second)).unwrap();
    result.add_score(a, 3.0, Some(first)).unwrap();
    result.add_score(b, 4.0, Some(second)).unwrap();
    assert_eq!(result.score(a), Some(2.0));
    assert_eq!(result.score_at_step(a, Some(first)), Some(3.0));
    assert_eq!(result.number_of_scores(), 3);
    assert_eq!(
        result
            .steps_and_scores
            .iter()
            .map(|s| s.processing_step)
            .collect::<Vec<_>>(),
        [Some(first), Some(second)]
    );
    assert_eq!(
        result.steps_and_scores[1]
            .scores_in_order(&[b, b, a], false)
            .unwrap(),
        [(b, 4.0), (b, 4.0), (a, 2.0)]
    );
    assert_eq!(
        result
            .most_recent_score(|id| Ok(&graph
                .processing_software(graph.processing_step(id)?.software)?
                .assigned_scores))
            .unwrap(),
        Some((b, 4.0))
    );
    result.add_score(a, 9.0, None).unwrap();
    assert_eq!(result.score_and_step(a), Some((9.0, None)));
    assert_eq!(
        result.steps_and_scores[2]
            .scores_in_order(&[b], false)
            .unwrap(),
        [(a, 9.0)]
    );
    // Source current-step literal EDIT is attached even without a score.
    graph.set_current_processing_step(first).unwrap();
    let id = graph
        .register_identified_peptide(IdentifiedPeptide::new(peptide("EDIT")))
        .unwrap();
    assert_eq!(graph.peptide(id).unwrap().result.steps_and_scores.len(), 1);
    assert_eq!(
        graph.peptide(id).unwrap().result.steps_and_scores[0].processing_step,
        Some(first)
    );
}

#[test]
fn source_merge_and_copy_translate_references_with_distinct_current_step_rules() {
    // IdentificationData.cpp:1036-1090,1102-1141,1209-1220.
    let mut incoming = IdentificationData::new().unwrap();
    incoming
        .set_metadata(BTreeMap::from([("graph".into(), "incoming".into())]))
        .unwrap();
    let software = incoming
        .register_processing_software(ProcessingSoftware::new("search", "1"))
        .unwrap();
    let old = incoming
        .register_db_search_param(DBSearchParam {
            database: "old.fasta".into(),
            ..DBSearchParam::default()
        })
        .unwrap();
    let new = incoming
        .register_db_search_param(DBSearchParam {
            database: "new.fasta".into(),
            ..DBSearchParam::default()
        })
        .unwrap();
    let step = incoming
        .register_processing_step(ProcessingStep::new(software), Some(old))
        .unwrap();
    incoming
        .register_processing_step(ProcessingStep::new(software), Some(new))
        .unwrap();
    assert_eq!(incoming.search_param_for_step(step).unwrap(), Some(old));
    incoming.set_current_processing_step(step).unwrap();
    let source_parent = incoming
        .register_parent_sequence(parent("rna", "ACGU", MoleculeType::RNA))
        .unwrap();
    let mut hit = IdentifiedOligo::new(oligo("ACG"));
    hit.parent_matches
        .insert(source_parent, BTreeSet::from([matched(0, 2)]));
    let source_hit = incoming.register_identified_oligo(hit).unwrap();

    let (copy, translated) = incoming.try_clone_with_translation().unwrap();
    assert_eq!(copy.metadata()["graph"], "incoming".into());
    assert_eq!(
        copy.current_processing_step(),
        Some(translated.processing_step(step).unwrap())
    );
    let copied = copy.oligo(translated.oligo(source_hit).unwrap()).unwrap();
    assert_eq!(copied.result.steps_and_scores.len(), 1);
    assert_eq!(
        copied.parent_matches.keys().copied().collect::<Vec<_>>(),
        [translated.parent(source_parent).unwrap()]
    );
    assert!(copy.oligo(source_hit).is_err());
    assert!(
        incoming
            .oligo(translated.oligo(source_hit).unwrap())
            .is_err()
    );

    let mut destination = IdentificationData::new().unwrap();
    destination
        .set_metadata(BTreeMap::from([("graph".into(), "destination".into())]))
        .unwrap();
    let destination_software = destination
        .register_processing_software(ProcessingSoftware::new("search", "1"))
        .unwrap();
    let other_search = destination
        .register_db_search_param(DBSearchParam {
            database: "destination.fasta".into(),
            ..DBSearchParam::default()
        })
        .unwrap();
    let equivalent = destination
        .register_processing_step(
            ProcessingStep::new(destination_software),
            Some(other_search),
        )
        .unwrap();
    let audit_software = destination
        .register_processing_software(ProcessingSoftware::new("audit", "1"))
        .unwrap();
    let current = destination
        .register_processing_step(ProcessingStep::new(audit_software), None)
        .unwrap();
    destination.set_current_processing_step(current).unwrap();
    let translation = destination.merge_from(&incoming).unwrap();
    assert_eq!(destination.metadata()["graph"], "destination".into());
    assert_eq!(destination.current_processing_step(), Some(current));
    assert_eq!(translation.processing_step(step).unwrap(), equivalent);
    assert_eq!(
        destination.search_param_for_step(equivalent).unwrap(),
        Some(translation.db_search_param(old).unwrap())
    );
    let merged = destination
        .oligo(translation.oligo(source_hit).unwrap())
        .unwrap();
    assert_eq!(
        merged
            .result
            .steps_and_scores
            .iter()
            .map(|value| value.processing_step)
            .collect::<Vec<_>>(),
        [Some(equivalent), Some(current)]
    );
    let preserved = translation.oligo(source_hit).unwrap();
    // Repeated source merge neither duplicates a node nor moves an existing step.
    assert_eq!(
        destination
            .merge_from(&incoming)
            .unwrap()
            .oligo(source_hit)
            .unwrap(),
        preserved
    );
    assert_eq!(
        destination
            .oligo(preserved)
            .unwrap()
            .result
            .steps_and_scores
            .len(),
        2
    );
}

#[test]
fn source_coverage_filters_intervals_and_retains_empty_parent_break() {
    let mut graph = IdentificationData::new().unwrap();
    let empty = graph
        .register_parent_sequence(parent("empty", "", MoleculeType::Protein))
        .unwrap();
    let full = graph
        .register_parent_sequence(parent("full", "PEPTIDE", MoleculeType::Protein))
        .unwrap();
    let mut hit = IdentifiedPeptide::new(peptide("PEP"));
    hit.parent_matches.insert(empty, BTreeSet::new());
    hit.parent_matches
        .insert(full, BTreeSet::from([matched(0, 2)]));
    graph.register_identified_peptide(hit).unwrap();
    graph.calculate_coverages(false).unwrap();
    assert_eq!(graph.parent(full).unwrap().coverage, 0.0);

    let mut other = IdentifiedPeptide::new(peptide("EDIT"));
    // Wrong length is accepted unless requested; reversed/unknown/outside never count.
    other.parent_matches.insert(
        full,
        BTreeSet::from([
            matched(2, 4),
            matched(6, 1),
            matched(6, 8),
            ParentMatch::default(),
        ]),
    );
    graph.register_identified_peptide(other).unwrap();
    graph.calculate_coverages(false).unwrap();
    assert_eq!(graph.parent(full).unwrap().coverage, 3.0 / 7.0);
    graph.calculate_coverages(true).unwrap();
    assert_eq!(graph.parent(full).unwrap().coverage, 0.0);
    assert!(!ParentMatch::new(Some(usize::MAX), Some(usize::MAX)).has_valid_positions(0, 7));
}

#[test]
fn native_owned_chemical_sequences_do_not_collapse_equal_display_strings() {
    use openms::chemistry::{
        EmpiricalFormula, Ribonucleotide, RibonucleotideDB, RibonucleotideRecord,
    };
    // This is a value-identity native boundary: source pointers are not serialized IDs.
    let isotope = Ribonucleotide::from_record(RibonucleotideRecord {
        name: "labeled A".into(),
        code: "A".into(),
        origin: 'A',
        formula: EmpiricalFormula::parse("(13)C10H13N5O4").unwrap(),
        ..RibonucleotideRecord::default()
    })
    .unwrap();
    let registry = RibonucleotideDB::from_records(vec![isotope]).unwrap();
    let labeled = NASequence::parse_with_registry("AA", &registry).unwrap();
    let ordinary = oligo("AA");
    assert_eq!(labeled.to_string(), ordinary.to_string());
    assert_ne!(labeled, ordinary);
    let mut graph = IdentificationData::new().unwrap();
    let first = graph
        .register_identified_oligo(IdentifiedOligo::new(labeled.clone()))
        .unwrap();
    let second = graph
        .register_identified_oligo(IdentifiedOligo::new(ordinary))
        .unwrap();
    assert_ne!(first, second);
    assert_eq!(graph.oligo_count(), 2);
    drop(registry);
    assert_eq!(graph.oligo(first).unwrap().sequence, labeled);
    let (copy, trans) = graph.try_clone_with_translation().unwrap();
    assert_eq!(copy.oligo_count(), 2);
    assert_eq!(
        copy.oligo(trans.oligo(first).unwrap()).unwrap().sequence,
        labeled
    );
    graph.clear().unwrap();
    graph
        .register_identified_oligo(IdentifiedOligo::new(oligo("AA")))
        .unwrap();
    assert!(graph.oligo(first).is_err());
    assert!(graph.oligo(second).is_err());
}

#[test]
fn source_rnase_graph_literals_and_derived_complete_parent_matches() {
    // Literal source counts/starts: RNaseDigestion_test.cpp:163-198.
    // Product strings, inclusive ends and flanks are independently derived from
    // RNaseDigestion.cpp:190-229, and are distinguished in the manifest.
    let mut digestion = RNaseDigestion::new("RNase_T1").unwrap();
    let mut graph = IdentificationData::new().unwrap();
    let rna = graph
        .register_parent_sequence(parent("test", "pAUGUCGCAG", MoleculeType::RNA))
        .unwrap();
    graph
        .register_parent_sequence(parent(
            "skip",
            "unparsed protein text",
            MoleculeType::Protein,
        ))
        .unwrap();
    digestion.digest_identification_data(&mut graph).unwrap();
    assert_eq!(graph.oligo_count(), 3);
    let actual: BTreeMap<_, _> = graph
        .oligos()
        .map(|(_, hit)| {
            assert_eq!(hit.parent_matches.len(), 1);
            let matches = &hit.parent_matches[&rna];
            assert_eq!(matches.len(), 1);
            let m = matches.first().unwrap();
            (
                hit.sequence.to_string(),
                (
                    m.start_pos.unwrap(),
                    m.end_pos.unwrap(),
                    m.left_neighbor.clone(),
                    m.right_neighbor.clone(),
                ),
            )
        })
        .collect();
    assert_eq!(
        actual,
        BTreeMap::from([
            ("pAUGp".into(), (0, 2, "[".into(), "U".into())),
            ("UCGp".into(), (3, 5, "G".into(), "C".into())),
            ("CAG".into(), (6, 8, "G".into(), "]".into())),
        ])
    );
    digestion.digest_identification_data(&mut graph).unwrap();
    assert_eq!(graph.oligo_count(), 3);

    let mut repeated = IdentificationData::new().unwrap();
    let rna = repeated
        .register_parent_sequence(parent("test", "ACUGACUGG", MoleculeType::RNA))
        .unwrap();
    digestion.min_length = 2;
    digestion.digest_identification_data(&mut repeated).unwrap();
    assert_eq!(repeated.oligo_count(), 1);
    let (_, hit) = repeated.oligos().next().unwrap();
    assert_eq!(hit.sequence.to_string(), "ACUGp");
    assert_eq!(hit.parent_matches.len(), 1);
    let matches: Vec<_> = hit.parent_matches[&rna]
        .iter()
        .map(|m| {
            (
                m.start_pos.unwrap(),
                m.end_pos.unwrap(),
                m.left_neighbor.as_str(),
                m.right_neighbor.as_str(),
            )
        })
        .collect();
    assert_eq!(matches, [(0, 3, "[", "A"), (4, 7, "G", "G")]);
    repeated.calculate_coverages(true).unwrap();
    assert_eq!(repeated.parent(rna).unwrap().coverage, 8.0 / 9.0);
}

#[test]
fn native_late_conflicts_leave_merges_and_rna_digestion_unchanged() {
    // Atomicity is an explicit native guarantee, not a C++ exception gold standard.
    let mut destination = IdentificationData::new().unwrap();
    let original = destination
        .register_parent_sequence(parent("z", "ACG", MoleculeType::RNA))
        .unwrap();
    let mut incoming = IdentificationData::new().unwrap();
    incoming
        .register_parent_sequence(parent("a", "UG", MoleculeType::RNA))
        .unwrap();
    incoming
        .register_parent_sequence(parent("z", "UG", MoleculeType::RNA))
        .unwrap();
    let before = destination.parent(original).unwrap().clone();
    assert!(destination.merge_from(&incoming).is_err());
    assert_eq!(destination.parent_count(), 1);
    assert_eq!(destination.parent(original).unwrap(), &before);

    let mut graph = IdentificationData::new().unwrap();
    graph
        .register_parent_sequence(parent("a", "ACUG", MoleculeType::RNA))
        .unwrap();
    graph
        .register_parent_sequence(parent("z", "[missing-rna-code]", MoleculeType::RNA))
        .unwrap();
    let preserved = graph
        .register_identified_oligo(IdentifiedOligo::new(oligo("AA")))
        .unwrap();
    let digestion = RNaseDigestion::new("RNase_T1").unwrap();
    assert!(digestion.digest_identification_data(&mut graph).is_err());
    assert_eq!(graph.oligo_count(), 1);
    assert_eq!(graph.oligo(preserved).unwrap().sequence, oligo("AA"));
}
