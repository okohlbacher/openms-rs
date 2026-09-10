// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

// Source contracts: Core SDK 6bfc0e4711105f4eda2fea86812a83af7c7e791f,
// METADATA/ID/{InputFile,ScoreType,AppliedProcessingStep,ScoredProcessingResult,
// ParentSequence,ParentMatch,IdentifiedSequence,DBSearchParam}.h. The simple
// test_score/Tool/test.mzML values also appear in IdentificationData_test.cpp.
use openms::chemistry::{
    AASequence, DigestionEnzymeRNA, DigestionEnzymeRNARecord, DigestionSpecificity, NASequence,
    ProteaseDB, Ribonucleotide, RibonucleotideRecord,
};
use openms::identification::graph::*;
use openms::metadata::{CVTerm, MetaValue, MetaValueData, ProcessingAction};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

fn scores_and_steps() -> (
    IdentificationData,
    ScoreTypeId,
    ScoreTypeId,
    ProcessingStepId,
    ProcessingStepId,
) {
    let mut graph = IdentificationData::new().unwrap();
    let first = graph
        .register_score_type(ScoreType::new("test_score", true))
        .unwrap();
    let second = graph
        .register_score_type(ScoreType::new("another_score", false))
        .unwrap();
    let mut software = ProcessingSoftware::new("Tool", "1.0");
    software.assigned_scores = vec![second, first, second];
    let software = graph.register_processing_software(software).unwrap();
    let mut step = ProcessingStep::new(software);
    step.date_time = Some("2024-01-01 12:00:00".parse().unwrap());
    let a = graph.register_processing_step(step.clone(), None).unwrap();
    step.date_time = Some("2024-01-02 12:00:00".parse().unwrap());
    let b = graph.register_processing_step(step, None).unwrap();
    (graph, first, second, a, b)
}

#[test]
fn source_input_file_merge_fills_unions_and_rejects_conflicts_atomically() {
    let mut original = InputFile::new("test.mzML");
    original.primary_files.insert("sample.raw".into());
    let mut incoming = InputFile::new("test.mzML");
    incoming.experimental_design_id = "experiment-A".into();
    incoming.primary_files.insert("sample2.raw".into());
    original.merge(&incoming).unwrap();
    assert_eq!(original.experimental_design_id, "experiment-A");
    assert_eq!(original.primary_files.len(), 2);
    let saved = original.clone();
    incoming.experimental_design_id = "experiment-B".into();
    incoming.primary_files.insert("late.raw".into());
    assert!(original.merge(&incoming).is_err());
    assert_eq!(original, saved);
}

#[test]
fn score_name_only_cv_is_valid_without_weakening_generic_cv_validation() {
    let mut graph = IdentificationData::new().unwrap();
    let mut score = ScoreType::new("test_score", true);
    score.cv_term.value = "score description".into();
    assert!(score.cv_term.validate().is_err());
    let id = graph.register_score_type(score.clone()).unwrap();
    score.cv_term.value = "replacement ignored by source".into();
    score.metadata.insert("new".into(), 1.into());
    assert_eq!(graph.register_score_type(score).unwrap(), id);
    assert_eq!(
        graph
            .score_type(id)
            .unwrap()
            .cv_term
            .value
            .as_str()
            .unwrap(),
        "score description"
    );
    assert!(graph.score_type(id).unwrap().metadata.is_empty());
    assert!(
        graph
            .register_score_type(ScoreType::new("test_score", false))
            .is_err()
    );
    assert!(graph.register_score_type(ScoreType::default()).is_err());
    assert!(ScoreType::new("high", true).is_better_score(2.0, 1.0));
    assert!(ScoreType::new("low", false).is_better_score(1.0, 2.0));
    assert!(!ScoreType::new("high", true).is_better_score(1.0, 1.0));
}

#[test]
fn software_priority_duplicates_and_unlisted_scores_follow_source_order() {
    let (_graph, first, second, step, _) = scores_and_steps();
    let applied = AppliedProcessingStep {
        processing_step: Some(step),
        scores: BTreeMap::from([(first, 10.0), (second, 20.0)]),
    };
    assert_eq!(
        applied.scores_in_order(&[second, second], false).unwrap(),
        vec![(second, 20.0), (second, 20.0), (first, 10.0)]
    );
    assert_eq!(
        applied.scores_in_order(&[second, first], true).unwrap(),
        vec![(second, 20.0)]
    );
    let no_step = AppliedProcessingStep {
        processing_step: None,
        ..applied
    };
    assert_eq!(
        no_step.scores_in_order(&[second], false).unwrap(),
        vec![(first, 10.0), (second, 20.0)]
    );
}

#[test]
fn score_history_updates_do_not_move_existing_steps_and_metadata_overwrites() {
    let (graph, first, second, a, b) = scores_and_steps();
    let mut result = ScoredProcessingResult::default();
    result.add_score(first, 1.0, Some(a)).unwrap();
    result.add_score(first, 2.0, Some(b)).unwrap();
    result.add_score(first, 99.0, Some(a)).unwrap();
    assert_eq!(result.score(first), Some(2.0));
    assert_eq!(result.score_and_step(first), Some((2.0, Some(b))));
    assert_eq!(result.score_at_step(first, Some(a)), Some(99.0));
    assert_eq!(
        result
            .steps_and_scores
            .iter()
            .map(|s| s.processing_step)
            .collect::<Vec<_>>(),
        vec![Some(a), Some(b)]
    );
    result.add_score(second, 3.0, Some(b)).unwrap();
    result.metadata.insert("sample".into(), "old".into());
    let mut incoming = ScoredProcessingResult::default();
    incoming.metadata.insert("sample".into(), "new".into());
    incoming.add_score(second, 4.0, Some(a)).unwrap();
    result.merge(&incoming).unwrap();
    assert_eq!(result.metadata["sample"].as_str().unwrap(), "new");
    assert_eq!(result.number_of_scores(), 4);
    assert_eq!(
        result
            .most_recent_score(|id| Ok(graph
                .processing_software(graph.processing_step(id)?.software)?
                .assigned_scores
                .as_slice()))
            .unwrap(),
        Some((second, 3.0))
    );
    result.add_score(first, 5.0, None).unwrap();
    assert_eq!(
        result
            .steps_by_processing_step()
            .unwrap()
            .iter()
            .map(|s| s.processing_step)
            .collect::<Vec<_>>(),
        vec![None, Some(a), Some(b)]
    );
    assert_eq!(result.score_and_step(first), Some((5.0, None)));
    assert_eq!(
        result
            .most_recent_score(|id| Ok(graph
                .processing_software(graph.processing_step(id)?.software)?
                .assigned_scores
                .as_slice()))
            .unwrap(),
        Some((first, 5.0))
    );
}

#[test]
fn score_only_step_deduplicates_and_invalid_updates_preserve_all_state() {
    let (_graph, first, second, _, _) = scores_and_steps();
    let mut result = ScoredProcessingResult::default();
    result.add_score(first, -1.0, None).unwrap();
    result.add_score(second, 0.0, None).unwrap();
    assert_eq!(result.steps_and_scores.len(), 1);
    let saved = result.clone();
    assert!(result.add_score(first, f64::NAN, None).is_err());
    assert_eq!(result, saved);
    let mut bad = ScoredProcessingResult::default();
    bad.metadata.insert("late".into(), "must not appear".into());
    bad.steps_and_scores.push(AppliedProcessingStep {
        processing_step: None,
        scores: BTreeMap::from([(first, f64::INFINITY)]),
    });
    assert!(result.merge(&bad).is_err());
    assert_eq!(result, saved);
    bad.steps_and_scores = vec![
        AppliedProcessingStep::default(),
        AppliedProcessingStep::default(),
    ];
    assert!(result.merge(&bad).is_err());
    assert_eq!(result, saved);
}

#[test]
fn parent_merge_retains_coverage_type_and_key_but_adds_decoy_and_history() {
    let (_graph, score, _, step, _) = scores_and_steps();
    let mut original = ParentSequence::new("rna_1");
    original.molecule_type = MoleculeType::RNA;
    original.coverage = 0.25;
    let mut incoming = ParentSequence::new("rna_1");
    incoming.sequence = "ACGU".into();
    incoming.description = "RNA parent".into();
    incoming.coverage = 0.75;
    incoming.is_decoy = true;
    incoming.result.add_score(score, 3.0, Some(step)).unwrap();
    original.merge(&incoming).unwrap();
    assert_eq!(original.sequence, "ACGU");
    assert_eq!(original.coverage, 0.25);
    assert_eq!(original.molecule_type, MoleculeType::RNA);
    assert!(original.is_decoy);
    assert_eq!(original.result.score(score), Some(3.0));
    let saved = original.clone();
    incoming.description = "conflicting description".into();
    incoming
        .result
        .metadata
        .insert("late".into(), "discard".into());
    assert!(original.merge(&incoming).is_err());
    assert_eq!(original, saved);
    incoming.description.clone_from(&original.description);
    incoming.sequence = "GGGG".into();
    assert!(original.merge(&incoming).is_err());
    assert_eq!(original, saved);
}

#[test]
fn parent_matches_use_inclusive_positions_only_and_checked_extremes() {
    let mut first = ParentMatch::new(Some(4), Some(10));
    first.left_neighbor = "first".into();
    first.metadata.insert("origin".into(), "first".into());
    let mut equal = first.clone();
    equal.left_neighbor = "replacement".into();
    equal.metadata.insert("origin".into(), "replacement".into());
    assert_eq!(first, equal);
    let mut set = BTreeSet::new();
    assert!(set.insert(first.clone()));
    assert!(!set.insert(equal));
    assert_eq!(set.len(), 1);
    assert_eq!(set.first().unwrap().left_neighbor, "first");
    assert!(first.has_valid_positions(7, 11));
    assert!(!first.has_valid_positions(6, 11));
    assert!(!first.has_valid_positions(7, 10));
    assert!(first.has_valid_positions(0, 0));
    assert!(!ParentMatch::default().has_valid_positions(0, 0));
    assert!(!ParentMatch::new(Some(2), Some(1)).has_valid_positions(0, 0));
    assert!(!ParentMatch::new(Some(0), Some(usize::MAX)).has_valid_positions(0, 0));
    // Native Option identity keeps an explicitly supplied maximum separate from unknown.
    assert!(ParentMatch::new(Some(usize::MAX), Some(usize::MAX)).has_valid_positions(0, 0));
    assert!(
        !ParentMatch::new(Some(usize::MAX), Some(usize::MAX)).has_valid_positions(0, usize::MAX)
    );
    assert!(first < ParentMatch::default());
    assert!(ParentMatch::new(Some(usize::MAX), Some(usize::MAX)) < ParentMatch::default());
}

#[test]
fn identified_sequence_merges_keep_first_match_payload_and_preserve_custom_chemistry() {
    let mut graph = IdentificationData::new().unwrap();
    let mut parent = ParentSequence::new("rna");
    parent.molecule_type = MoleculeType::RNA;
    parent.is_decoy = true;
    let parent = graph.register_parent_sequence(parent).unwrap();
    let sequence = NASequence::parse("ACGU").unwrap();
    let mut first = IdentifiedOligo::new(sequence.clone());
    let mut a = ParentMatch::new(Some(0), Some(3));
    a.left_neighbor = "first".into();
    first.parent_matches.insert(parent, BTreeSet::from([a]));
    let mut second = IdentifiedOligo::new(sequence);
    let mut b = ParentMatch::new(Some(0), Some(3));
    b.left_neighbor = "second".into();
    second.parent_matches.insert(
        parent,
        BTreeSet::from([b, ParentMatch::new(Some(4), Some(7))]),
    );
    first.merge(&second).unwrap();
    assert_eq!(first.parent_matches[&parent].len(), 2);
    assert_eq!(
        first.parent_matches[&parent].first().unwrap().left_neighbor,
        "first"
    );
    assert!(
        first
            .all_parents_are_decoys(|id| Ok(graph.parent(id)?.is_decoy))
            .unwrap()
    );
    assert!(
        IdentifiedOligo::new(NASequence::parse("A").unwrap())
            .all_parents_are_decoys(|_| Ok(true))
            .is_err()
    );
    let make = |mass| {
        Arc::new(
            Ribonucleotide::from_record(RibonucleotideRecord {
                code: "same".into(),
                mono_mass: mass,
                ..RibonucleotideRecord::default()
            })
            .unwrap(),
        )
    };
    let left = IdentifiedOligo::new(NASequence::from_records(vec![make(1.0)]).unwrap());
    let right = IdentifiedOligo::new(NASequence::from_records(vec![make(2.0)]).unwrap());
    assert_eq!(left.sequence.to_string(), right.sequence.to_string());
    assert_ne!(left.sequence, right.sequence);
    assert_ne!(
        graph.register_identified_oligo(left).unwrap(),
        graph.register_identified_oligo(right).unwrap()
    );
}

#[test]
fn search_parameter_key_has_all_source_fields_but_no_metadata() {
    let mut graph = IdentificationData::new().unwrap();
    let mut parameter = DBSearchParam {
        database: "test-db.fasta".into(),
        precursor_mass_tolerance: 1.0,
        fragment_mass_tolerance: 2.0,
        ..DBSearchParam::default()
    };
    let base = graph.register_db_search_param(parameter.clone()).unwrap();
    parameter.metadata.insert("non-key".into(), 7.into());
    assert_eq!(
        graph.register_db_search_param(parameter.clone()).unwrap(),
        base
    );
    assert!(graph.db_search_param(base).unwrap().metadata.is_empty());
    let mut variations = Vec::new();
    macro_rules! vary {
        ($field:ident,$value:expr) => {{
            let mut p = parameter.clone();
            p.$field = $value;
            variations.push(p);
        }};
    }
    vary!(molecule_type, MoleculeType::RNA);
    vary!(mass_type, MassType::Average);
    vary!(database, "different".into());
    vary!(database_version, "v2".into());
    vary!(taxonomy, "human".into());
    vary!(charges, BTreeSet::from([-2, 2]));
    vary!(fixed_mods, BTreeSet::from(["x".into()]));
    vary!(variable_mods, BTreeSet::from(["x".into()]));
    vary!(precursor_mass_tolerance, 3.0);
    vary!(fragment_mass_tolerance, 4.0);
    vary!(precursor_tolerance_ppm, true);
    vary!(fragment_tolerance_ppm, true);
    vary!(
        digestion_enzyme,
        Some(GraphEnzyme::Protein(
            *ProteaseDB::global().get_enzyme("Trypsin").unwrap()
        ))
    );
    vary!(enzyme_term_specificity, Some(DigestionSpecificity::Full));
    vary!(missed_cleavages, 2);
    vary!(min_length, 5);
    vary!(max_length, 20);
    let mut ids = BTreeSet::from([base]);
    for value in variations {
        assert!(ids.insert(graph.register_db_search_param(value).unwrap()));
    }
    assert_eq!(ids.len(), 18);
    let zero = graph
        .register_db_search_param(DBSearchParam::default())
        .unwrap();
    let neg_zero = DBSearchParam {
        precursor_mass_tolerance: -0.0,
        fragment_mass_tolerance: -0.0,
        ..DBSearchParam::default()
    };
    assert_eq!(graph.register_db_search_param(neg_zero).unwrap(), zero);
    assert!(
        graph
            .register_db_search_param(DBSearchParam {
                precursor_mass_tolerance: -1.0,
                ..DBSearchParam::default()
            })
            .is_ok()
    );
    assert!(
        graph
            .register_db_search_param(DBSearchParam {
                fragment_mass_tolerance: f64::NAN,
                ..DBSearchParam::default()
            })
            .is_err()
    );
}

#[test]
fn owned_rna_enzyme_complete_value_is_part_of_search_key() {
    let mut graph = IdentificationData::new().unwrap();
    let enzyme = |after: &str| {
        Arc::new(
            DigestionEnzymeRNA::from_record(DigestionEnzymeRNARecord {
                name: "same enzyme".into(),
                cuts_after: after.into(),
                ..DigestionEnzymeRNARecord::default()
            })
            .unwrap(),
        )
    };
    let a = DBSearchParam {
        digestion_enzyme: Some(GraphEnzyme::RNA(enzyme("G"))),
        ..DBSearchParam::default()
    };
    let b = DBSearchParam {
        digestion_enzyme: Some(GraphEnzyme::RNA(enzyme("A"))),
        ..DBSearchParam::default()
    };
    let id = graph.register_db_search_param(a.clone()).unwrap();
    assert_eq!(graph.register_db_search_param(a).unwrap(), id);
    assert_ne!(graph.register_db_search_param(b).unwrap(), id);
}

#[test]
fn step_key_preserves_timestamp_input_order_and_actions_but_ignores_metadata() {
    let mut graph = IdentificationData::new().unwrap();
    let software = graph
        .register_processing_software(ProcessingSoftware::new("Tool", "1"))
        .unwrap();
    let a = graph.register_input_file(InputFile::new("a")).unwrap();
    let b = graph.register_input_file(InputFile::new("b")).unwrap();
    let mut step = ProcessingStep::new(software);
    step.input_files = vec![a, b];
    step.date_time = Some("2024-01-01 00:00:00".parse().unwrap());
    let first = graph.register_processing_step(step.clone(), None).unwrap();
    step.metadata.insert("new".into(), 1.into());
    assert_eq!(
        graph.register_processing_step(step.clone(), None).unwrap(),
        first
    );
    assert!(graph.processing_step(first).unwrap().metadata.is_empty());
    step.input_files.reverse();
    assert_ne!(
        graph.register_processing_step(step.clone(), None).unwrap(),
        first
    );
    step.input_files.reverse();
    step.actions.insert(ProcessingAction::Identification);
    assert_ne!(
        graph.register_processing_step(step.clone(), None).unwrap(),
        first
    );
    step.actions.clear();
    step.date_time = None;
    assert_ne!(graph.register_processing_step(step, None).unwrap(), first);
}

#[test]
fn draft_payloads_are_validated_and_shared_work_caps_precede_graph_mutation() {
    let mut graph = IdentificationData::new().unwrap();
    let mut parent = ParentSequence::new("valid");
    parent.coverage = f64::NAN;
    assert!(graph.register_parent_sequence(parent).is_err());
    let mut score = ScoreType::new("finite", true);
    assert!(MetaValue::new(MetaValueData::Float(f64::INFINITY)).is_err());
    score.cv_term.accession = "MS: bad".into();
    assert!(graph.register_score_type(score).is_err());
    let mut software = ProcessingSoftware::new("invalid cv", "1");
    assert!(software.software.cv_terms.add(CVTerm::default()).is_err());
    let mut limited = IdentificationData::with_limits(GraphLimits {
        max_work: 1024,
        ..GraphLimits::default()
    })
    .unwrap();
    let mut large = InputFile::new("oversized");
    large.primary_files.insert("x".repeat(4096));
    assert!(limited.register_input_file(large).is_err());
    assert!(limited.is_empty());
    assert!(
        graph
            .register_identified_peptide(IdentifiedPeptide::new(AASequence::parse("").unwrap()))
            .is_err()
    );
    assert!(
        graph
            .register_identified_oligo(IdentifiedOligo::new(NASequence::parse("").unwrap()))
            .is_err()
    );
}
