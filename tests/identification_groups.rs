// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Source literal cases and derived branch expectations are distinguished in
// tests/data/identification_groups_provenance.json. No C++ execution.
use openms::chemistry::{AASequence, AdductInfo};
use openms::identification::graph::*;
use std::collections::{BTreeMap, BTreeSet};

fn parent(graph: &mut IdentificationData, name: &str, kind: MoleculeType) -> ParentId {
    let mut value = ParentSequence::new(name);
    value.molecule_type = kind;
    graph.register_parent_sequence(value).unwrap()
}
fn step(graph: &mut IdentificationData, name: &str) -> ProcessingStepId {
    let software = graph
        .register_processing_software(ProcessingSoftware::new(name, "1"))
        .unwrap();
    graph
        .register_processing_step(ProcessingStep::new(software), None)
        .unwrap()
}
fn matches(graph: &mut IdentificationData) -> [ObservationMatchId; 4] {
    let file = graph.register_input_file(InputFile::new("file")).unwrap();
    let first = graph
        .register_observation(Observation::new("one", file))
        .unwrap();
    let second = graph
        .register_observation(Observation::new("two", file))
        .unwrap();
    let p = graph
        .register_identified_peptide(IdentifiedPeptide::new(
            AASequence::parse("PEPTIDE").unwrap(),
        ))
        .unwrap();
    let c = graph
        .register_identified_compound(IdentifiedCompound::new("compound"))
        .unwrap();
    let a = graph
        .register_observation_match(ObservationMatch::new(p, first))
        .unwrap();
    let b = graph
        .register_observation_match(ObservationMatch::new(p, second))
        .unwrap();
    let d = graph
        .register_observation_match(ObservationMatch::new(c, first))
        .unwrap();
    let adduct = graph
        .register_adduct(AdductInfo::parse("M+H;1+").unwrap())
        .unwrap();
    let mut ion = ObservationMatch::new(p, first);
    ion.adduct = Some(adduct);
    let e = graph.register_observation_match(ion).unwrap();
    [a, b, d, e]
}

#[test]
fn source_parent_group_literal_mixed_types_and_duplicate_labels_append() {
    let mut graph = IdentificationData::new().unwrap();
    let protein = parent(&mut graph, "protein_1", MoleculeType::Protein);
    let rna = parent(&mut graph, "rna_1", MoleculeType::RNA);
    let mut groups = ParentGroupSet::new("test_grouping");
    groups
        .groups
        .push(ParentGroup::new(BTreeSet::from([protein, rna])));
    let first = graph.register_parent_group_set(groups.clone()).unwrap();
    assert_eq!(graph.parent_group_set_count(), 1);
    assert_eq!(graph.parent_group_count(), 1);
    assert_eq!(graph.parent_group(first, 0).unwrap().parent_refs.len(), 2);
    let second = graph.register_parent_group_set(groups).unwrap();
    assert_ne!(first, second);
    assert_eq!(
        graph
            .parent_group_sets()
            .map(|(id, _)| id)
            .collect::<Vec<_>>(),
        [first, second]
    );
    assert_eq!(graph.parent_group_set_count(), 2);
    assert!(graph.parent_group(first, 1).is_err());
    let empty = graph
        .register_parent_group_set(ParentGroupSet::default())
        .unwrap();
    assert!(graph.parent_group_set(empty).unwrap().groups.is_empty());
    assert_eq!(graph.record_count(), 7); // two parents, three sets, two nested groups
}

#[test]
fn parent_group_normalization_keeps_first_score_map_in_reference_order() {
    let mut graph = IdentificationData::new().unwrap();
    let first = parent(&mut graph, "z", MoleculeType::Protein);
    let second = parent(&mut graph, "a", MoleculeType::Protein);
    let score = graph
        .register_score_type(ScoreType::new("score", true))
        .unwrap();
    let make = |refs, score_value| ParentGroup {
        parent_refs: refs,
        scores: BTreeMap::from([(score, score_value)]),
    };
    let value = ParentGroupSet {
        groups: vec![
            make(BTreeSet::from([second]), 5.0),
            make(BTreeSet::from([first]), 3.0),
            make(BTreeSet::from([second]), 99.0),
            ParentGroup::default(),
        ],
        ..ParentGroupSet::default()
    };
    let id = graph.register_parent_group_set(value).unwrap();
    let groups = &graph.parent_group_set(id).unwrap().groups;
    assert_eq!(groups.len(), 3);
    assert!(groups[0].parent_refs.is_empty());
    assert_eq!(groups[1].parent_refs, BTreeSet::from([first]));
    assert_eq!(groups[2].scores[&score], 5.0);
}

#[test]
fn group_history_applies_current_without_reordering_and_membership_deduplicates() {
    let mut graph = IdentificationData::new().unwrap();
    let [a, b, c, _] = matches(&mut graph);
    let score = graph
        .register_score_type(ScoreType::new("score", true))
        .unwrap();
    let first = step(&mut graph, "first");
    let second = step(&mut graph, "second");
    graph.set_current_processing_step(first).unwrap();
    let mut value = ObservationMatchGroup::new(BTreeSet::from([a, b, c]));
    value.result.add_score(score, 1.0, None).unwrap();
    value.result.metadata.insert("key".into(), "before".into());
    let id = graph
        .register_observation_match_group(value.clone())
        .unwrap();
    assert_eq!(graph.observation_match_group_count(), 1);
    graph.set_current_processing_step(second).unwrap();
    value.result.metadata.insert("key".into(), "after".into());
    value.result.add_score(score, 2.0, Some(first)).unwrap();
    assert_eq!(graph.register_observation_match_group(value).unwrap(), id);
    graph.set_current_processing_step(first).unwrap();
    graph
        .register_observation_match_group(ObservationMatchGroup::new(BTreeSet::from([c, a, b])))
        .unwrap();
    let got = graph.observation_match_group(id).unwrap();
    assert_eq!(got.result.metadata["key"], "after".into());
    assert_eq!(got.result.score(score), Some(2.0));
    assert_eq!(
        got.result
            .steps_and_scores
            .iter()
            .map(|s| s.processing_step)
            .collect::<Vec<_>>(),
        [None, Some(first), Some(second)]
    );
    let parent_set = ParentGroupSet {
        result: got.result.clone(),
        ..ParentGroupSet::default()
    };
    let set = graph.register_parent_group_set(parent_set).unwrap();
    assert_eq!(
        graph
            .parent_group_set(set)
            .unwrap()
            .result
            .steps_and_scores
            .len(),
        3
    );
}

#[test]
fn same_molecule_and_query_use_references_with_empty_singleton_source_rules() {
    let mut graph = IdentificationData::new().unwrap();
    let [a, b, c, d] = matches(&mut graph);
    for (refs, same_molecule, same_query) in [
        (vec![], true, true),
        (vec![a], true, true),
        (vec![a, b], true, false),
        (vec![a, c], false, true),
        (vec![a, d], true, true),
        (vec![a, b, c], false, false),
    ] {
        let id = graph
            .register_observation_match_group(ObservationMatchGroup::new(
                refs.into_iter().collect(),
            ))
            .unwrap();
        assert_eq!(
            graph.match_group_all_same_molecule(id).unwrap(),
            same_molecule
        );
        assert_eq!(graph.match_group_all_same_query(id).unwrap(), same_query);
    }
    let singleton = ObservationMatchGroup::new(BTreeSet::from([a]));
    assert!(
        singleton
            .all_same_molecule(|_| panic!("source singleton must not resolve"))
            .unwrap()
    );
    let empty = ObservationMatchGroup::default();
    assert!(
        empty
            .all_same_query(|_| panic!("source empty must not resolve"))
            .unwrap()
    );
}

#[test]
fn group_references_nonfinite_scores_and_discarded_entries_are_checked_atomically() {
    let mut graph = IdentificationData::new().unwrap();
    let [a, ..] = matches(&mut graph);
    let local = parent(&mut graph, "local", MoleculeType::Protein);
    let score = graph
        .register_score_type(ScoreType::new("score", true))
        .unwrap();
    let mut other = IdentificationData::new().unwrap();
    let [foreign_match, ..] = matches(&mut other);
    let foreign_parent = parent(&mut other, "foreign", MoleculeType::Protein);
    let foreign_score = other
        .register_score_type(ScoreType::new("foreign", true))
        .unwrap();
    let foreign_step = step(&mut other, "foreign");
    let before = graph.record_count();
    for group in [
        ParentGroup {
            parent_refs: BTreeSet::from([foreign_parent]),
            scores: BTreeMap::new(),
        },
        ParentGroup {
            parent_refs: BTreeSet::from([local]),
            scores: BTreeMap::from([(foreign_score, 1.0)]),
        },
        ParentGroup {
            parent_refs: BTreeSet::from([local]),
            scores: BTreeMap::from([(score, f64::NAN)]),
        },
    ] {
        let value = ParentGroupSet {
            groups: vec![ParentGroup::new(BTreeSet::from([local])), group],
            ..ParentGroupSet::default()
        };
        assert!(graph.register_parent_group_set(value).is_err());
        assert_eq!(graph.record_count(), before);
    }
    let id = graph
        .register_observation_match_group(ObservationMatchGroup::new(BTreeSet::from([a])))
        .unwrap();
    let saved = graph.observation_match_group(id).unwrap().clone();
    assert!(
        graph
            .register_observation_match_group(ObservationMatchGroup::new(BTreeSet::from([
                foreign_match
            ])))
            .is_err()
    );
    let mut invalid = saved.clone();
    invalid
        .result
        .steps_and_scores
        .push(AppliedProcessingStep::new(Some(foreign_step)));
    assert!(graph.register_observation_match_group(invalid).is_err());
    assert_eq!(graph.observation_match_group(id).unwrap(), &saved);
}

#[test]
fn graph_copy_and_merge_translate_both_group_families_and_preserve_current_steps() {
    let mut graph = IdentificationData::new().unwrap();
    let [a, b, ..] = matches(&mut graph);
    let p = parent(&mut graph, "p", MoleculeType::Protein);
    let score = graph
        .register_score_type(ScoreType::new("group score", false))
        .unwrap();
    let current = step(&mut graph, "source-step");
    graph.set_current_processing_step(current).unwrap();
    let mut set = ParentGroupSet::new("repeated-label");
    set.groups.push(ParentGroup {
        parent_refs: BTreeSet::from([p]),
        scores: BTreeMap::from([(score, -7.0)]),
    });
    set.result.add_score(score, 0.2, None).unwrap();
    let set_id = graph.register_parent_group_set(set.clone()).unwrap();
    let another = graph.register_parent_group_set(set).unwrap();
    let mut group = ObservationMatchGroup::new(BTreeSet::from([a, b]));
    group.result.add_score(score, 0.5, None).unwrap();
    group
        .result
        .metadata
        .insert("payload".into(), "retained".into());
    let group_id = graph.register_observation_match_group(group).unwrap();
    let (copy, translation) = graph.try_clone_with_translation().unwrap();
    assert_eq!(copy.parent_group_set_count(), 2);
    assert_eq!(copy.observation_match_group_count(), 1);
    let copied = copy
        .parent_group_set(translation.parent_group_set(set_id).unwrap())
        .unwrap();
    assert_eq!(
        copied.groups[0].parent_refs,
        BTreeSet::from([translation.parent(p).unwrap()])
    );
    assert_eq!(
        copied.groups[0].scores[&translation.score_type(score).unwrap()],
        -7.0
    );
    assert_eq!(copied.result.steps_and_scores.len(), 2);
    assert_ne!(
        translation.parent_group_set(set_id).unwrap(),
        translation.parent_group_set(another).unwrap()
    );
    let copied_group = copy
        .observation_match_group(translation.observation_match_group(group_id).unwrap())
        .unwrap();
    assert_eq!(
        copied_group.observation_match_refs,
        BTreeSet::from([
            translation.observation_match(a).unwrap(),
            translation.observation_match(b).unwrap()
        ])
    );
    assert_eq!(copied_group.result.metadata["payload"], "retained".into());
    assert_eq!(
        copy.current_processing_step(),
        Some(translation.processing_step(current).unwrap())
    );
    let mut destination = IdentificationData::new().unwrap();
    let destination_step = step(&mut destination, "destination-step");
    destination
        .set_current_processing_step(destination_step)
        .unwrap();
    let translated = destination.merge_from(&graph).unwrap();
    let result = &destination
        .observation_match_group(translated.observation_match_group(group_id).unwrap())
        .unwrap()
        .result;
    assert_eq!(result.steps_and_scores.len(), 3);
    assert_eq!(
        result.steps_and_scores.last().unwrap().processing_step,
        Some(destination_step)
    );
    graph.clear().unwrap();
    assert_eq!(graph.parent_group_count(), 0);
    assert_eq!(graph.observation_match_group_count(), 0);
    assert!(
        graph.parent_group_set(set_id).is_err() && graph.observation_match_group(group_id).is_err()
    );
    assert!(copy.parent_group_set(set_id).is_err());
}

#[test]
fn nested_groups_count_toward_record_and_edge_limits_and_failed_merge_rolls_back() {
    let mut graph = IdentificationData::with_limits(GraphLimits {
        max_records: 2,
        ..GraphLimits::default()
    })
    .unwrap();
    let value = ParentGroupSet {
        groups: vec![ParentGroup::default()],
        ..ParentGroupSet::default()
    };
    graph.register_parent_group_set(value).unwrap();
    assert_eq!(graph.record_count(), 2);
    assert!(
        graph
            .register_parent_group_set(ParentGroupSet::default())
            .is_err()
    );
    let mut graph = IdentificationData::with_limits(GraphLimits {
        max_edges: 0,
        ..GraphLimits::default()
    })
    .unwrap();
    let p = parent(&mut graph, "p", MoleculeType::Protein);
    let value = ParentGroupSet {
        groups: vec![ParentGroup::new(BTreeSet::from([p]))],
        ..ParentGroupSet::default()
    };
    assert!(graph.register_parent_group_set(value).is_err());
    assert_eq!(graph.record_count(), 1);
    let mut source = IdentificationData::new().unwrap();
    let [a, ..] = matches(&mut source);
    source
        .register_parent_group_set(ParentGroupSet {
            groups: vec![ParentGroup::default()],
            ..ParentGroupSet::default()
        })
        .unwrap();
    source
        .register_observation_match_group(ObservationMatchGroup::new(BTreeSet::from([a])))
        .unwrap();
    let mut destination = IdentificationData::with_limits(GraphLimits {
        max_records: source.record_count() - 1,
        ..GraphLimits::default()
    })
    .unwrap();
    assert!(destination.merge_from(&source).is_err());
    assert!(destination.is_empty());
}

#[test]
fn native_group_equality_includes_metadata_and_direct_merge_preserves_membership() {
    let mut graph = IdentificationData::new().unwrap();
    let [a, b, ..] = matches(&mut graph);
    let mut one = ObservationMatchGroup::new(BTreeSet::from([a]));
    let mut other = ObservationMatchGroup::new(BTreeSet::from([b]));
    other.result.metadata.insert("key".into(), "value".into());
    one.merge(&other).unwrap();
    assert_eq!(one.observation_match_refs, BTreeSet::from([a]));
    assert_eq!(one.result.metadata["key"], "value".into());
    let mut changed = one.clone();
    changed
        .result
        .metadata
        .insert("key".into(), "different".into());
    assert_ne!(one, changed);
}
