// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! `pickScoreType`, `applyToObservations` and `RefTranslator::allow_missing`
//! of `IdentificationData.h` (core bc9cc12).
//!
//! Tier 1: every expected value below is a record of the executed Linux x86_64
//! Release build (install `openms4-release-bc9cc12-c19e494-174b576`, g++ 14.4.0,
//! run twice on ibminode06 with identical output), oracle
//! `../oracle/A5b-iddata` (driver.cpp sha256 b22c98b0…0472563, records.tsv
//! sha256 d5501799…0670fc4). Each test names its oracle case. The two cases
//! where this port deliberately differs (`pick_all_count_b2`,
//! `pick_any_count_a2`, `apply_rename_collides`) assert the documented
//! contract and state the source's executed answer beside it. No class-test
//! section of `IdentificationData_test.cpp` covers these members.

use openms::chemistry::AASequence;
use openms::identification::graph::*;
use std::collections::BTreeMap;

struct Scores {
    graph: IdentificationData,
    a: ScoreTypeId,
    b: ScoreTypeId,
    c: ScoreTypeId,
    step: ProcessingStepId,
}

/// Oracle setup: score_a (higher better), score_b (lower better), score_c
/// registered in that order (the executed build placed their set nodes at
/// ascending addresses: `setup address_rank_0..2` = a, b, c); software
/// "engine" assigns [b, a]; one processing step.
fn scores() -> Scores {
    let mut graph = IdentificationData::new().unwrap();
    let file = graph
        .register_input_file(InputFile::new("file.mzML"))
        .unwrap();
    let a = graph
        .register_score_type(ScoreType::new("score_a", true))
        .unwrap();
    let b = graph
        .register_score_type(ScoreType::new("score_b", false))
        .unwrap();
    let c = graph
        .register_score_type(ScoreType::new("score_c", true))
        .unwrap();
    let mut software = ProcessingSoftware::new("engine", "1.0");
    software.assigned_scores = vec![b, a];
    let software = graph.register_processing_software(software).unwrap();
    let mut step = ProcessingStep::new(software);
    step.input_files.push(file);
    let step = graph.register_processing_step(step, None).unwrap();
    Scores {
        graph,
        a,
        b,
        c,
        step,
    }
}

/// One applied step: its optional processing step and its scores.
type Step<'a> = (Option<ProcessingStepId>, &'a [(ScoreTypeId, f64)]);

fn result(steps: &[Step]) -> ScoredProcessingResult {
    let mut value = ScoredProcessingResult::default();
    for (step, scores) in steps {
        value
            .add_processing_step(AppliedProcessingStep {
                processing_step: *step,
                scores: scores.iter().copied().collect::<BTreeMap<_, _>>(),
            })
            .unwrap();
    }
    value
}

/// p1..p5 of the oracle driver.
fn elements(s: &Scores) -> [ScoredProcessingResult; 5] {
    [
        result(&[(Some(s.step), &[(s.a, 1.0), (s.b, 2.0)])]),
        result(&[(None, &[(s.c, 5.0)])]),
        result(&[(None, &[(s.a, 3.0), (s.c, 1.0)])]),
        result(&[(Some(s.step), &[(s.a, 1.0)]), (None, &[(s.c, 2.0)])]),
        result(&[]),
    ]
}

fn pick(
    s: &Scores,
    container: &[&ScoredProcessingResult],
    all: bool,
    any: bool,
) -> Option<ScoreTypeId> {
    s.graph
        .pick_score_type(container.iter().copied(), all, any)
        .unwrap()
}

#[test]
fn pick_score_type_first_element_follows_the_most_recent_primary_score() {
    let s = scores();
    let [p1, p2, p3, p4, p5] = elements(&s);
    // pick_first_skips_unscored: score_b (the software ranks b before a).
    assert_eq!(pick(&s, &[&p5, &p1, &p2], false, false), Some(s.b));
    // pick_first_no_step: score_c.
    assert_eq!(pick(&s, &[&p2, &p1], false, false), Some(s.c));
    // pick_empty and pick_only_unscored: getScoreTypes().end().
    assert_eq!(pick(&s, &[], false, false), None);
    assert_eq!(pick(&s, &[&p5], false, false), None);
    // pick_first_map_order: score_a (no step, so the lowest key).
    assert_eq!(pick(&s, &[&p3], false, false), Some(s.a));
    // pick_first_later_step: score_c (the most recent step wins).
    assert_eq!(pick(&s, &[&p4], false, false), Some(s.c));
}

#[test]
fn pick_score_type_counting_agrees_with_the_source_where_the_count_and_the_key_agree() {
    let s = scores();
    let [p1, p2, p3, p4, p5] = elements(&s);
    // pick_all_count_c2: score_c (b 1, c 2, a 1).
    assert_eq!(pick(&s, &[&p1, &p2, &p4, &p3], true, false), Some(s.c));
    // pick_all_single_a: score_a.
    assert_eq!(pick(&s, &[&p3], true, false), Some(s.a));
    // pick_any_single_c: score_c.
    assert_eq!(pick(&s, &[&p2], false, true), Some(s.c));
    // pick_any_tie_ab: score_b (a 1, b 1; the later-registered type wins the tie).
    assert_eq!(pick(&s, &[&p1], false, true), Some(s.b));
    // pick_any_unscored: end().
    assert_eq!(pick(&s, &[&p5], true, true), None);
}

#[test]
fn pick_score_type_counting_follows_the_documented_contract_not_the_source_address_rule() {
    let s = scores();
    let [p1, p2, p3, _, _] = elements(&s);
    // pick_all_count_b2: b is the most recent score of two elements, c of one.
    // The documented contract picks b; the executed source returned score_c,
    // the counted type with the highest set-node address (max_element over
    // std::pair compares the key first).
    assert_eq!(pick(&s, &[&p1, &p1, &p2], true, false), Some(s.b));
    // pick_any_count_a2: a 2, b 1, c 1. Contract: a; executed source: score_c.
    assert_eq!(pick(&s, &[&p1, &p3], false, true), Some(s.a));
}

#[test]
fn pick_score_type_refuses_foreign_ids() {
    let s = scores();
    let mut other = IdentificationData::new().unwrap();
    let foreign = other
        .register_score_type(ScoreType::new("foreign", true))
        .unwrap();
    let element = result(&[(None, &[(foreign, 1.0)])]);
    assert!(s.graph.pick_score_type([&element], false, false).is_err());
    assert!(s.graph.pick_score_type([&element], true, true).is_err());
}

fn observations() -> (IdentificationData, InputFileId, ObservationMatchId) {
    let mut graph = IdentificationData::new().unwrap();
    let file = graph
        .register_input_file(InputFile::new("obs.mzML"))
        .unwrap();
    let mut first = None;
    for (id, rt, mz) in [
        ("spec=1", 10.0, 500.0),
        ("spec=2", 20.0, 600.0),
        ("spec=3", 30.0, 700.0),
    ] {
        let mut observation = Observation::new(id, file);
        observation.rt = Some(rt);
        observation.mz = Some(mz);
        let observation = graph.register_observation(observation).unwrap();
        first.get_or_insert(observation);
    }
    let peptide = graph
        .register_identified_peptide(IdentifiedPeptide::new(
            AASequence::parse("PEPTIDE").unwrap(),
        ))
        .unwrap();
    let mut matched = ObservationMatch::new(peptide, first.unwrap());
    matched.charge = 2;
    let matched = graph.register_observation_match(matched).unwrap();
    (graph, file, matched)
}

fn listing(graph: &IdentificationData) -> Vec<(String, f64)> {
    graph
        .observations()
        .map(|(_, o)| (o.data_id.clone(), o.rt.unwrap()))
        .collect()
}

fn expect(pairs: &[(&str, f64)]) -> Vec<(String, f64)> {
    pairs.iter().map(|(id, rt)| (id.to_string(), *rt)).collect()
}

#[test]
fn apply_to_observations_modifies_in_key_order_and_reorders_changed_keys() {
    let (mut graph, _, matched) = observations();
    // apply_setup.
    assert_eq!(
        listing(&graph),
        expect(&[("spec=1", 10.0), ("spec=2", 20.0), ("spec=3", 30.0)])
    );
    // apply_shift_rt: every retention time + 5, order unchanged.
    graph
        .apply_to_observations(|o| o.rt = o.rt.map(|rt| rt + 5.0))
        .unwrap();
    assert_eq!(
        listing(&graph),
        expect(&[("spec=1", 15.0), ("spec=2", 25.0), ("spec=3", 35.0)])
    );
    // apply_rename_reorders: spec=3 becomes spec=0 and moves to the front.
    let renamed = graph
        .observations()
        .find(|(_, o)| o.data_id == "spec=3")
        .unwrap()
        .0;
    graph
        .apply_to_observations(|o| {
            if o.data_id == "spec=3" {
                o.data_id = "spec=0".into();
            }
        })
        .unwrap();
    assert_eq!(
        listing(&graph),
        expect(&[("spec=0", 35.0), ("spec=1", 15.0), ("spec=2", 25.0)])
    );
    assert_eq!(graph.observation(renamed).unwrap().data_id, "spec=0");
    assert_eq!(graph.observation_count(), 3);
    assert_eq!(graph.observation_match_count(), 1);
    assert_eq!(
        graph
            .observation(graph.observation_match(matched).unwrap().observation)
            .unwrap()
            .data_id,
        "spec=1"
    );
}

#[test]
fn apply_to_observations_refuses_a_key_collision_the_source_resolves_by_erasing() {
    let (mut graph, _, _) = observations();
    let before = listing(&graph);
    // apply_rename_collides: the executed source erased the renamed element
    // (count 3 -> 2). Here the call is refused and nothing changes, not even
    // the observations the closure visited before the collision.
    let result = graph.apply_to_observations(|o| {
        o.rt = o.rt.map(|rt| rt + 1.0);
        if o.data_id == "spec=2" {
            o.data_id = "spec=1".into();
        }
    });
    assert!(result.is_err());
    assert_eq!(listing(&graph), before);
    assert_eq!(graph.observation_count(), 3);
}

#[test]
fn apply_to_observations_refuses_invalid_records_atomically() {
    let (mut graph, _, _) = observations();
    let before = listing(&graph);
    assert!(
        graph
            .apply_to_observations(|o| o.rt = Some(f64::NAN))
            .is_err()
    );
    assert!(graph.apply_to_observations(|o| o.data_id.clear()).is_err());
    let mut other = IdentificationData::new().unwrap();
    let foreign = other.register_input_file(InputFile::new("x")).unwrap();
    assert!(
        graph
            .apply_to_observations(|o| o.input_file = foreign)
            .is_err()
    );
    assert_eq!(listing(&graph), before);
}

#[test]
fn reference_translator_allow_missing_governs_molecules_and_matches_only() {
    let mut source = IdentificationData::new().unwrap();
    let file = source
        .register_input_file(InputFile::new("src.mzML"))
        .unwrap();
    let observation = source
        .register_observation(Observation::new("scan=7", file))
        .unwrap();
    let peptide = source
        .register_identified_peptide(IdentifiedPeptide::new(
            AASequence::parse("SAMPLER").unwrap(),
        ))
        .unwrap();
    let mut matched = ObservationMatch::new(peptide, observation);
    matched.charge = 3;
    let matched = source.register_observation_match(matched).unwrap();

    let mut foreign = IdentificationData::new().unwrap();
    let foreign_file = foreign
        .register_input_file(InputFile::new("foreign.mzML"))
        .unwrap();
    let foreign_observation = foreign
        .register_observation(Observation::new("scan=9", foreign_file))
        .unwrap();
    let foreign_peptide = foreign
        .register_identified_peptide(IdentifiedPeptide::new(
            AASequence::parse("FOREIGN").unwrap(),
        ))
        .unwrap();
    let mut foreign_match = ObservationMatch::new(foreign_peptide, foreign_observation);
    foreign_match.charge = 1;
    let foreign_match = foreign.register_observation_match(foreign_match).unwrap();

    let mut merged = IdentificationData::new().unwrap();
    let mut translator = merged.merge_from(&source).unwrap();
    assert!(!translator.allow_missing);

    // translate_molecule_present: the merged peptide, sequence SAMPLER.
    let IdentifiedMolecule::Peptide(new_peptide) = translator.molecule(peptide.into()).unwrap()
    else {
        panic!("a peptide translates to a peptide");
    };
    assert_eq!(
        merged.peptide(new_peptide).unwrap().sequence.to_string(),
        "SAMPLER"
    );
    assert_eq!(merged.peptides().next().unwrap().0, new_peptide);
    // translate_match_present: the merged match.
    let new_match = translator.observation_match(matched).unwrap();
    assert_eq!(merged.observation_matches().next().unwrap().0, new_match);

    // allow_missing_false: MissingInformation for both.
    assert!(translator.molecule(foreign_peptide.into()).is_err());
    assert!(translator.observation_match(foreign_match).is_err());
    // allow_missing_true: the old reference comes back unchanged.
    translator.allow_missing = true;
    assert_eq!(
        translator.molecule(foreign_peptide.into()).unwrap(),
        IdentifiedMolecule::Peptide(foreign_peptide)
    );
    assert_eq!(
        translator.observation_match(foreign_match).unwrap(),
        foreign_match
    );
    // The plain tables are std::map lookups in the source and stay errors.
    assert!(translator.peptide(foreign_peptide).is_err());
    assert!(translator.observation(foreign_observation).is_err());
}
