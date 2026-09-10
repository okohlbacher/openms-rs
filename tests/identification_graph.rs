// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::chemistry::{
    AASequence, ModificationsDB, NASequence, Ribonucleotide, RibonucleotideDB, RibonucleotideRecord,
};
use openms::identification::graph::*;
use openms::metadata::MetaValue;
use std::collections::{BTreeMap, BTreeSet};

fn step(graph: &mut IdentificationData, name: &str) -> ProcessingStepId {
    let software = graph
        .register_processing_software(ProcessingSoftware::new(name, "1"))
        .unwrap();
    graph
        .register_processing_step(ProcessingStep::new(software), None)
        .unwrap()
}
fn rna(name: &str, sequence: &str) -> ParentSequence {
    let mut parent = ParentSequence::new(name);
    parent.molecule_type = MoleculeType::RNA;
    parent.sequence = sequence.into();
    parent
}

#[test]
fn files_scores_and_software_preserve_source_keys_and_first_payload() {
    let mut graph = IdentificationData::new().unwrap();
    let score = graph
        .register_score_type(ScoreType::new("test_score", true))
        .unwrap();
    assert!(
        graph
            .score_type(score)
            .unwrap()
            .cv_term
            .accession
            .is_empty()
    );
    assert!(
        graph
            .register_score_type(ScoreType::new("test_score", false))
            .is_err()
    );
    let mut file = InputFile::new("test.mzML");
    file.primary_files.insert("one.raw".into());
    let id = graph.register_input_file(file).unwrap();
    let mut file = InputFile::new("test.mzML");
    file.experimental_design_id = "sample1".into();
    file.primary_files.insert("two.raw".into());
    assert_eq!(graph.register_input_file(file.clone()).unwrap(), id);
    assert_eq!(graph.input_file(id).unwrap().primary_files.len(), 2);
    file.experimental_design_id = "sample2".into();
    let before = graph.input_file(id).unwrap().clone();
    assert!(graph.register_input_file(file).is_err());
    assert_eq!(graph.input_file(id).unwrap(), &before);
    let software = ProcessingSoftware::new("Tool", "1.0");
    let software_id = graph
        .register_processing_software(software.clone())
        .unwrap();
    let mut later = software;
    later.assigned_scores.push(score);
    assert_eq!(
        graph.register_processing_software(later).unwrap(),
        software_id
    );
    assert!(
        graph
            .processing_software(software_id)
            .unwrap()
            .assigned_scores
            .is_empty()
    );
}

#[test]
fn foreign_stale_and_failed_step_links_never_become_local_references() {
    let mut graph = IdentificationData::new().unwrap();
    let mut foreign = IdentificationData::new().unwrap();
    let other_score = foreign
        .register_score_type(ScoreType::new("foreign", true))
        .unwrap();
    let mut software = ProcessingSoftware::new("Tool", "1");
    software.assigned_scores.push(other_score);
    assert!(graph.register_processing_software(software).is_err());
    assert!(graph.is_empty());
    let local = graph
        .register_processing_software(ProcessingSoftware::new("Tool", "1"))
        .unwrap();
    let foreign_search = foreign
        .register_db_search_param(DBSearchParam::default())
        .unwrap();
    assert!(
        graph
            .register_processing_step(ProcessingStep::new(local), Some(foreign_search))
            .is_err()
    );
    assert_eq!(graph.processing_step_count(), 0);
    let id = graph
        .register_parent_sequence(ParentSequence::new("P"))
        .unwrap();
    assert!(foreign.parent(id).is_err());
    graph.clear().unwrap();
    let fresh = graph
        .register_parent_sequence(ParentSequence::new("P"))
        .unwrap();
    assert_ne!(fresh, id);
    assert!(graph.parent(id).is_err());
}

#[test]
fn parent_merges_preserve_original_type_coverage_and_atomic_conflicts() {
    let mut graph = IdentificationData::new().unwrap();
    let mut first = ParentSequence::new("P");
    first.coverage = 0.3;
    first
        .result
        .metadata
        .insert("key".into(), MetaValue::from("old"));
    let id = graph.register_parent_sequence(first).unwrap();
    let mut update = rna("P", "ACGU");
    update.coverage = 0.9;
    update.is_decoy = true;
    update.description = "description".into();
    update
        .result
        .metadata
        .insert("key".into(), MetaValue::from("new"));
    graph.register_parent_sequence(update.clone()).unwrap();
    let result = graph.parent(id).unwrap();
    assert_eq!(result.coverage, 0.3);
    assert_eq!(result.molecule_type, MoleculeType::Protein);
    assert!(result.is_decoy);
    assert_eq!(result.sequence, "ACGU");
    assert_eq!(result.result.metadata["key"].as_str().unwrap(), "new");
    let before = result.clone();
    update.description = "conflict".into();
    update
        .result
        .metadata
        .insert("late".into(), MetaValue::from(1));
    assert!(graph.register_parent_sequence(update).is_err());
    assert_eq!(graph.parent(id).unwrap(), &before);
}

#[test]
fn sequences_union_positions_and_check_even_empty_parent_edges() {
    let mut graph = IdentificationData::new().unwrap();
    let parent = graph.register_parent_sequence(rna("R", "ACGU")).unwrap();
    let protein = graph
        .register_parent_sequence(ParentSequence::new("P"))
        .unwrap();
    let mut oligo = IdentifiedOligo::new("AC".parse().unwrap());
    let mut first = ParentMatch::new(Some(0), Some(1));
    first.left_neighbor = "first".into();
    oligo.parent_matches.insert(parent, BTreeSet::from([first]));
    let id = graph.register_identified_oligo(oligo.clone()).unwrap();
    let mut later = ParentMatch::new(Some(0), Some(1));
    later.left_neighbor = "later".into();
    oligo.parent_matches.insert(
        parent,
        BTreeSet::from([
            later,
            ParentMatch::new(None, None),
            ParentMatch::new(Some(3), Some(2)),
        ]),
    );
    assert_eq!(graph.register_identified_oligo(oligo).unwrap(), id);
    let matches = &graph.oligo(id).unwrap().parent_matches[&parent];
    assert_eq!(matches.len(), 3);
    assert_eq!(matches.first().unwrap().left_neighbor, "first");
    let mut invalid = IdentifiedOligo::new("UG".parse().unwrap());
    invalid.parent_matches.insert(protein, BTreeSet::new());
    assert!(graph.register_identified_oligo(invalid).is_err());
    assert_eq!(graph.oligo_count(), 1);
    assert!(
        graph
            .register_identified_oligo(IdentifiedOligo::new(NASequence::new()))
            .is_err()
    );
}

#[test]
fn current_step_attaches_on_reregistration_without_reordering_history() {
    let mut graph = IdentificationData::new().unwrap();
    let one = step(&mut graph, "one");
    let two = step(&mut graph, "two");
    let score = graph
        .register_score_type(ScoreType::new("score", true))
        .unwrap();
    graph.set_current_processing_step(one).unwrap();
    let id = graph
        .register_identified_peptide(IdentifiedPeptide::new("EDIT".parse().unwrap()))
        .unwrap();
    graph.add_peptide_score(id, score, 1.0).unwrap();
    graph.set_current_processing_step(two).unwrap();
    graph
        .register_identified_peptide(IdentifiedPeptide::new("EDIT".parse().unwrap()))
        .unwrap();
    graph.add_peptide_score(id, score, 2.0).unwrap();
    graph.set_current_processing_step(one).unwrap();
    let mut updated = IdentifiedPeptide::new("EDIT".parse().unwrap());
    updated.result.add_score(score, 3.0, Some(one)).unwrap();
    graph.register_identified_peptide(updated).unwrap();
    let result = &graph.peptide(id).unwrap().result;
    assert_eq!(
        result
            .steps_and_scores
            .iter()
            .map(|s| s.processing_step)
            .collect::<Vec<_>>(),
        vec![Some(one), Some(two)]
    );
    assert_eq!(result.score(score), Some(2.0));
    assert_eq!(result.score_at_step(score, Some(one)), Some(3.0));
    let before = result.clone();
    assert!(graph.add_peptide_score(id, score, f64::NAN).is_err());
    assert_eq!(&graph.peptide(id).unwrap().result, &before);
}

#[test]
fn distinct_custom_chemistry_is_not_deduplicated_by_display() {
    let make = |formula: &str| {
        RibonucleotideDB::from_records(vec![
            Ribonucleotide::from_record(RibonucleotideRecord {
                code: "A".into(),
                origin: 'A',
                formula: formula.parse().unwrap(),
                ..Default::default()
            })
            .unwrap(),
        ])
        .unwrap()
    };
    let first = NASequence::parse_with_registry("A", &make("C1H1")).unwrap();
    let second = NASequence::parse_with_registry("A", &make("C2H1")).unwrap();
    assert_eq!(first.to_string(), second.to_string());
    assert_ne!(first, second);
    let mut graph = IdentificationData::new().unwrap();
    assert_ne!(
        graph
            .register_identified_oligo(IdentifiedOligo::new(first))
            .unwrap(),
        graph
            .register_identified_oligo(IdentifiedOligo::new(second))
            .unwrap()
    );
    let original = ModificationsDB::global()
        .get_modification("Oxidation (M)", None, None)
        .unwrap()
        .clone();
    let modified = original.clone().with_absolute_masses(111.0, 112.0).unwrap();
    let first = AASequence::parse_with_registry(
        "M(Oxidation)",
        &ModificationsDB::from_records(vec![original]).unwrap(),
    )
    .unwrap();
    let second = AASequence::parse_with_registry(
        "M(Oxidation)",
        &ModificationsDB::from_records(vec![modified]).unwrap(),
    )
    .unwrap();
    assert_eq!(first.to_string(), second.to_string());
    assert_ne!(
        graph
            .register_identified_peptide(IdentifiedPeptide::new(first))
            .unwrap(),
        graph
            .register_identified_peptide(IdentifiedPeptide::new(second))
            .unwrap()
    );
}

#[test]
fn coverage_counts_residues_and_unions_intervals_atomically() {
    let mut graph = IdentificationData::new().unwrap();
    let parent = graph
        .register_parent_sequence(rna("RNA", "A[m3U]GCA"))
        .unwrap();
    let mut oligo = IdentifiedOligo::new("AGC".parse().unwrap());
    oligo.parent_matches.insert(
        parent,
        BTreeSet::from([
            ParentMatch::new(Some(1), Some(3)),
            ParentMatch::new(Some(2), Some(4)),
        ]),
    );
    graph.register_identified_oligo(oligo).unwrap();
    graph.calculate_coverages(true).unwrap();
    assert_eq!(graph.parent(parent).unwrap().coverage, 0.8);
    let invalid = graph
        .register_parent_sequence(rna("late", "[missing]"))
        .unwrap();
    let mut oligo = IdentifiedOligo::new("UU".parse().unwrap());
    oligo.parent_matches.insert(invalid, BTreeSet::new());
    graph.register_identified_oligo(oligo).unwrap();
    assert!(graph.calculate_coverages(false).is_err());
    assert_eq!(graph.parent(parent).unwrap().coverage, 0.8);
}

#[test]
fn empty_parent_preserves_source_break_in_native_id_order() {
    let mut graph = IdentificationData::new().unwrap();
    let empty = graph.register_parent_sequence(rna("empty", "")).unwrap();
    let present = graph
        .register_parent_sequence(rna("present", "ACGU"))
        .unwrap();
    let mut oligo = IdentifiedOligo::new("AC".parse().unwrap());
    oligo.parent_matches.insert(empty, BTreeSet::new());
    oligo.parent_matches.insert(
        present,
        BTreeSet::from([ParentMatch::new(Some(0), Some(1))]),
    );
    graph.register_identified_oligo(oligo).unwrap();
    graph.calculate_coverages(false).unwrap();
    assert_eq!(graph.parent(present).unwrap().coverage, 0.0);
}

#[test]
fn copy_and_merge_translate_every_reference_and_keep_distinct_current_steps() {
    let mut source = IdentificationData::new().unwrap();
    let source_step = step(&mut source, "source");
    let score = source
        .register_score_type(ScoreType::new("score", true))
        .unwrap();
    source.set_current_processing_step(source_step).unwrap();
    source
        .set_metadata(BTreeMap::from([("root".into(), MetaValue::from("source"))]))
        .unwrap();
    let parent = source.register_parent_sequence(rna("R", "ACGU")).unwrap();
    let mut oligo = IdentifiedOligo::new("AC".parse().unwrap());
    oligo
        .parent_matches
        .insert(parent, BTreeSet::from([ParentMatch::new(Some(0), Some(1))]));
    oligo
        .result
        .add_score(score, 4.0, Some(source_step))
        .unwrap();
    let id = source.register_identified_oligo(oligo).unwrap();
    let (copy, translation) = source.try_clone_with_translation().unwrap();
    assert!(copy.oligo(id).is_err());
    let copy_id = translation.oligo(id).unwrap();
    assert!(
        copy.oligo(copy_id)
            .unwrap()
            .parent_matches
            .contains_key(&translation.parent(parent).unwrap())
    );
    assert_eq!(
        copy.current_processing_step(),
        Some(translation.processing_step(source_step).unwrap())
    );
    assert_eq!(copy.metadata(), source.metadata());
    let mut destination = IdentificationData::new().unwrap();
    let current = step(&mut destination, "destination");
    destination.set_current_processing_step(current).unwrap();
    let map = destination.merge_from(&source).unwrap();
    assert_eq!(destination.current_processing_step(), Some(current));
    assert!(destination.metadata().is_empty());
    let result = &destination.oligo(map.oligo(id).unwrap()).unwrap().result;
    assert_eq!(
        result.steps_and_scores.last().unwrap().processing_step,
        Some(current)
    );
    assert_eq!(result.score(map.score_type(score).unwrap()), Some(4.0));
}

#[test]
fn merge_conflict_and_resource_errors_leave_destination_records_unchanged() {
    let mut destination = IdentificationData::new().unwrap();
    let parent = destination
        .register_parent_sequence(rna("Z", "AC"))
        .unwrap();
    let mut source = IdentificationData::new().unwrap();
    source.register_parent_sequence(rna("A", "UG")).unwrap();
    source.register_parent_sequence(rna("Z", "UG")).unwrap();
    assert!(destination.merge_from(&source).is_err());
    assert_eq!(destination.parent_count(), 1);
    assert_eq!(destination.parent(parent).unwrap().sequence, "AC");
    let mut bounded = IdentificationData::with_limits(GraphLimits {
        max_records: 1,
        ..Default::default()
    })
    .unwrap();
    let id = bounded
        .register_parent_sequence(ParentSequence::new("first"))
        .unwrap();
    assert!(
        bounded
            .register_parent_sequence(ParentSequence::new("second"))
            .is_err()
    );
    assert_eq!(bounded.parent_count(), 1);
    assert_eq!(bounded.parent(id).unwrap().accession, "first");
}

#[test]
fn ordinary_step_association_is_first_wins_but_merge_overwrites() {
    let mut graph = IdentificationData::new().unwrap();
    let software = graph
        .register_processing_software(ProcessingSoftware::new("tool", "1"))
        .unwrap();
    let a = DBSearchParam {
        database: "a".into(),
        ..Default::default()
    };
    let b = DBSearchParam {
        database: "b".into(),
        ..Default::default()
    };
    let a_id = graph.register_db_search_param(a.clone()).unwrap();
    let b_id = graph.register_db_search_param(b.clone()).unwrap();
    let step = graph
        .register_processing_step(ProcessingStep::new(software), Some(a_id))
        .unwrap();
    graph
        .register_processing_step(ProcessingStep::new(software), Some(b_id))
        .unwrap();
    assert_eq!(graph.search_param_for_step(step).unwrap(), Some(a_id));
    let mut other = IdentificationData::new().unwrap();
    let software = other
        .register_processing_software(ProcessingSoftware::new("tool", "1"))
        .unwrap();
    let search = other.register_db_search_param(b).unwrap();
    other
        .register_processing_step(ProcessingStep::new(software), Some(search))
        .unwrap();
    graph.merge_from(&other).unwrap();
    assert_eq!(graph.search_param_for_step(step).unwrap(), Some(b_id));
}

#[test]
fn returned_views_follow_semantic_order_without_reassigning_slots() {
    let mut graph = IdentificationData::new().unwrap();
    let z = graph
        .register_parent_sequence(ParentSequence::new("Z"))
        .unwrap();
    let a = graph
        .register_parent_sequence(ParentSequence::new("A"))
        .unwrap();
    assert_eq!(
        graph.parents().map(|(id, _)| id).collect::<Vec<_>>(),
        vec![a, z]
    );
    assert_eq!(graph.parent(z).unwrap().accession, "Z");
    let peptide = graph
        .register_identified_peptide(IdentifiedPeptide::new("TEST".parse().unwrap()))
        .unwrap();
    assert!(graph.peptide_parents_are_decoys(peptide).is_err());
}
