// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Source: Core SDK 6bfc0e4, METADATA/ID/{Observation,IdentifiedCompound,IdentificationData}.
// Literal spectrum/compound values come from IdentificationData_test.cpp; boundary,
// ordering and transaction expectations are independently derived from source branches.
use openms::chemistry::{
    AASequence, AdductInfo, EmpiricalFormula, NAFragmentType, NASequence, PeptideFragmentType,
};
use openms::identification::PeakAnnotation;
use openms::identification::graph::*;

fn input(graph: &mut IdentificationData, name: &str) -> InputFileId {
    graph.register_input_file(InputFile::new(name)).unwrap()
}
fn observation(graph: &mut IdentificationData, file: InputFileId, name: &str) -> ObservationId {
    graph
        .register_observation(Observation::new(name, file))
        .unwrap()
}
fn compound(graph: &mut IdentificationData, name: &str) -> CompoundId {
    graph
        .register_identified_compound(IdentifiedCompound::new(name))
        .unwrap()
}
fn step(graph: &mut IdentificationData, name: &str) -> ProcessingStepId {
    let software = graph
        .register_processing_software(ProcessingSoftware::new(name, "1"))
        .unwrap();
    graph
        .register_processing_step(ProcessingStep::new(software), None)
        .unwrap()
}
fn formula(text: &str) -> EmpiricalFormula {
    text.parse().unwrap()
}

#[test]
fn observation_identity_and_missing_coordinate_overwrite_are_source_faithful() {
    let mut graph = IdentificationData::new().unwrap();
    let file = input(&mut graph, "data.mzML");
    let first_step = step(&mut graph, "acquisition");
    graph.set_current_processing_step(first_step).unwrap();
    let mut value = Observation::new("spectrum_1", file);
    value.rt = Some(100.0);
    value.mz = Some(1000.0);
    value.metadata.insert("old".into(), 1.into());
    let id = graph.register_observation(value).unwrap();
    let mut update = Observation::new("spectrum_1", file);
    update.rt = Some(-100.0);
    update.metadata.insert("old".into(), 2.into());
    update.metadata.insert("new".into(), "retained".into());
    assert_eq!(graph.register_observation(update).unwrap(), id);
    let got = graph.observation(id).unwrap();
    assert_eq!((got.rt, got.mz), (Some(-100.0), None));
    assert_eq!(got.metadata["old"], 2.into());
    let second_file = input(&mut graph, "other.mzML");
    let distinct = observation(&mut graph, second_file, "spectrum_1");
    assert_ne!(id, distinct);
    assert_eq!(graph.observation_count(), 2);
    graph
        .set_observation_coordinates(id, None, Some(-1.0))
        .unwrap();
    graph
        .set_observation_meta_value(id, "old".into(), 3.into())
        .unwrap();
    assert_eq!(graph.observation(id).unwrap().mz, Some(-1.0));
    assert_eq!(
        graph.observation(id).unwrap().metadata["new"],
        "retained".into()
    );
    let saved = graph.observation(id).unwrap().clone();
    assert!(
        graph
            .set_observation_coordinates(id, Some(3.0), Some(f64::INFINITY))
            .is_err()
    );
    let mut invalid = saved.clone();
    invalid.rt = Some(f64::NAN);
    assert!(graph.register_observation(invalid).is_err());
    assert_eq!(graph.observation(id).unwrap(), &saved);
    let mut other = IdentificationData::new().unwrap();
    let foreign = input(&mut other, "data.mzML");
    assert!(
        graph
            .register_observation(Observation::new("foreign", foreign))
            .is_err()
    );
    assert!(
        graph
            .register_observation(Observation::new("", file))
            .is_err()
    );
}

#[test]
fn compound_first_chemical_payload_survives_but_results_and_current_step_merge() {
    let mut graph = IdentificationData::new().unwrap();
    let score = graph
        .register_score_type(ScoreType::new("score", true))
        .unwrap();
    let first = step(&mut graph, "first");
    let later = step(&mut graph, "later");
    graph.set_current_processing_step(first).unwrap();
    let mut value = IdentifiedCompound::new("compound_1");
    value.formula = formula("C2H5OH");
    value.name = "ethanol".into();
    value.smile = "CCO".into();
    value.inchi = "initial".into();
    value
        .result
        .metadata
        .insert("source".into(), "original".into());
    let id = graph.register_identified_compound(value.clone()).unwrap();
    graph.set_current_processing_step(later).unwrap();
    let mut update = IdentifiedCompound::new("compound_1");
    update.formula = formula("C6H12O6");
    update.name = "glucose".into();
    update.result.add_score(score, 7.0, Some(first)).unwrap();
    update
        .result
        .metadata
        .insert("source".into(), "updated".into());
    assert_eq!(graph.register_identified_compound(update).unwrap(), id);
    let got = graph.compound(id).unwrap();
    assert_eq!(
        (&got.formula, &got.name, &got.smile, &got.inchi),
        (&value.formula, &value.name, &value.smile, &value.inchi)
    );
    assert_eq!(got.result.metadata["source"], "updated".into());
    assert_eq!(
        got.result
            .steps_and_scores
            .iter()
            .map(|s| s.processing_step)
            .collect::<Vec<_>>(),
        [Some(first), Some(later)]
    );
    graph.add_compound_score(id, score, 9.0).unwrap();
    assert_eq!(graph.compound(id).unwrap().result.score(score), Some(9.0));
    let empty = compound(&mut graph, "empty-first");
    let mut filled = IdentifiedCompound::new("empty-first");
    filled.formula = formula("H2O");
    filled.name = "water".into();
    graph.register_identified_compound(filled).unwrap();
    assert!(graph.compound(empty).unwrap().formula.is_empty());
    assert_eq!(graph.compound(empty).unwrap().name, "");
    assert!(
        graph
            .register_identified_compound(IdentifiedCompound::new(""))
            .is_err()
    );
}

#[test]
fn adduct_registration_ignores_name_and_multiplier_but_preserves_formula_and_charge_keys() {
    let mut graph = IdentificationData::new().unwrap();
    let first = graph
        .register_adduct(AdductInfo::new("Na+", formula("Na"), 1, 1).unwrap())
        .unwrap();
    assert_eq!(
        graph
            .register_adduct(AdductInfo::new("another name", formula("Na"), 1, 2).unwrap())
            .unwrap(),
        first
    );
    assert_eq!(graph.adduct(first).unwrap().name(), "Na+");
    assert_eq!(graph.adduct(first).unwrap().mol_multiplier(), 1);
    let negative = graph
        .register_adduct(AdductInfo::new("Na+", formula("Na"), -1, 1).unwrap())
        .unwrap();
    let labeled = graph
        .register_adduct(AdductInfo::new("Na+", formula("(13)C"), 1, 1).unwrap())
        .unwrap();
    let carbon = graph
        .register_adduct(AdductInfo::new("Na+", formula("C"), 1, 1).unwrap())
        .unwrap();
    assert_ne!(labeled, carbon);
    assert_eq!(graph.adducts().next().unwrap().0, negative);
    assert_eq!(graph.adduct_count(), 4);
}

#[test]
fn source_best_match_rules_cover_mixed_types_ties_directions_and_unscored_groups() {
    let mut graph = IdentificationData::new().unwrap();
    let file = input(&mut graph, "file");
    let obs = observation(&mut graph, file, "spectrum_1");
    let one = observation(&mut graph, file, "spectrum_2");
    let two = observation(&mut graph, file, "spectrum_3");
    let empty = observation(&mut graph, file, "spectrum_4");
    let peptide = graph
        .register_identified_peptide(IdentifiedPeptide::new(
            AASequence::parse("PEPTIDE").unwrap(),
        ))
        .unwrap();
    let cmp = compound(&mut graph, "compound");
    let rna = graph
        .register_identified_oligo(IdentifiedOligo::new(NASequence::parse("ACGU").unwrap()))
        .unwrap();
    let higher = graph
        .register_score_type(ScoreType::new("higher", true))
        .unwrap();
    let lower = graph
        .register_score_type(ScoreType::new("lower", false))
        .unwrap();
    // Insert opposite to molecule arm order. Source key order, not insertion order, wins ties.
    let rna_match = graph
        .register_observation_match(ObservationMatch::new(rna, obs))
        .unwrap();
    let compound_match = graph
        .register_observation_match(ObservationMatch::new(cmp, obs))
        .unwrap();
    let peptide_match = graph
        .register_observation_match(ObservationMatch::new(peptide, obs))
        .unwrap();
    for id in [rna_match, compound_match, peptide_match] {
        graph
            .add_observation_match_score(id, higher, 100.0)
            .unwrap();
    }
    graph
        .add_observation_match_score(rna_match, lower, -10.0)
        .unwrap();
    graph
        .add_observation_match_score(peptide_match, lower, 0.0)
        .unwrap();
    let single = graph
        .register_observation_match(ObservationMatch::new(rna, one))
        .unwrap();
    graph
        .register_observation_match(ObservationMatch::new(rna, two))
        .unwrap();
    graph
        .register_observation_match(ObservationMatch::new(cmp, two))
        .unwrap();
    assert_eq!(
        graph
            .matches_for_observation(obs)
            .unwrap()
            .map(|(id, _)| id)
            .collect::<Vec<_>>(),
        [peptide_match, compound_match, rna_match]
    );
    assert_eq!(graph.matches_for_observation(empty).unwrap().count(), 0);
    assert_eq!(graph.best_matches(higher, true).unwrap(), [peptide_match]);
    assert_eq!(
        graph.best_matches(higher, false).unwrap(),
        [peptide_match, single]
    );
    assert_eq!(graph.best_matches(lower, true).unwrap(), [rna_match]);
    graph
        .add_observation_match_score(rna_match, higher, 200.0)
        .unwrap();
    assert_eq!(
        graph.best_matches(higher, false).unwrap(),
        [rna_match, single]
    );
    let foreign = IdentificationData::new().unwrap();
    assert!(foreign.best_matches(higher, false).is_err());
    assert!(foreign.matches_for_observation(obs).is_err());
}

#[test]
fn annotations_and_all_new_references_translate_and_old_ids_fail_after_clear() {
    let mut graph = IdentificationData::new().unwrap();
    let file = input(&mut graph, "file");
    let obs = observation(&mut graph, file, "scan");
    let cmp = compound(&mut graph, "compound");
    let score = graph
        .register_score_type(ScoreType::new("score", true))
        .unwrap();
    let current = step(&mut graph, "step");
    graph.set_current_processing_step(current).unwrap();
    let adduct = graph
        .register_adduct(AdductInfo::new("H+", formula("H"), 1, 1).unwrap())
        .unwrap();
    let mut value = ObservationMatch::new(cmp, obs);
    value.adduct = Some(adduct);
    value.charge = -2;
    value.result.add_score(score, 3.0, None).unwrap();
    let annotation = PeakAnnotation {
        annotation: "duplicate".into(),
        mz: -1.0,
        intensity: -3.0,
        charge: -2,
    };
    value
        .peak_annotations
        .insert(None, vec![annotation.clone(), annotation.clone()]);
    value.peak_annotations.insert(Some(current), Vec::new());
    let id = graph.register_observation_match(value).unwrap();
    graph
        .set_molecule_meta_value(cmp, "label".into(), "chemical".into())
        .unwrap();
    graph
        .set_observation_match_meta_value(id, "delete".into(), 1.into())
        .unwrap();
    graph
        .remove_observation_match_meta_value(id, "delete")
        .unwrap();
    let (copy, trans) = graph.try_clone_with_translation().unwrap();
    let got = copy
        .observation_match(trans.observation_match(id).unwrap())
        .unwrap();
    assert_eq!(got.identified_molecule, trans.molecule(cmp.into()).unwrap());
    assert_eq!(got.observation, trans.observation(obs).unwrap());
    assert_eq!(
        copy.observation(got.observation).unwrap().input_file,
        trans.input_file(file).unwrap()
    );
    assert_eq!(got.adduct, Some(trans.adduct(adduct).unwrap()));
    assert_eq!(
        got.peak_annotations[&None],
        [annotation.clone(), annotation]
    );
    assert!(got.peak_annotations[&Some(trans.processing_step(current).unwrap())].is_empty());
    assert_eq!(
        got.result.score(trans.score_type(score).unwrap()),
        Some(3.0)
    );
    assert_eq!(got.result.steps_and_scores.len(), 2);
    assert!(!got.result.metadata.contains_key("delete"));
    assert_eq!(
        copy.current_processing_step(),
        Some(trans.processing_step(current).unwrap())
    );
    assert_eq!(
        copy.compound(trans.compound(cmp).unwrap())
            .unwrap()
            .result
            .metadata["label"],
        "chemical".into()
    );
    graph.clear().unwrap();
    assert!(graph.is_empty());
    assert!(
        graph.observation(obs).is_err()
            && graph.compound(cmp).is_err()
            && graph.adduct(adduct).is_err()
            && graph.observation_match(id).is_err()
    );
    assert!(copy.observation_match(id).is_err());
}

#[test]
fn every_optional_match_reference_is_checked_before_registration() {
    let mut graph = IdentificationData::new().unwrap();
    let file = input(&mut graph, "file");
    let obs = observation(&mut graph, file, "scan");
    let cmp = compound(&mut graph, "c");
    let base = ObservationMatch::new(cmp, obs);
    let id = graph.register_observation_match(base.clone()).unwrap();
    let mut foreign = IdentificationData::new().unwrap();
    let other_step = step(&mut foreign, "foreign");
    let other_adduct = foreign
        .register_adduct(AdductInfo::new("H+", formula("H"), 1, 1).unwrap())
        .unwrap();
    let other_cmp = compound(&mut foreign, "c");
    for case in 0..3 {
        let mut value = base.clone();
        match case {
            0 => {
                value.peak_annotations.insert(Some(other_step), Vec::new());
            }
            1 => value.adduct = Some(other_adduct),
            _ => value.identified_molecule = other_cmp.into(),
        }
        assert!(graph.register_observation_match(value).is_err());
        assert_eq!(graph.observation_match(id).unwrap(), &base);
    }
    assert_eq!(graph.observation_match_count(), 1);
}

#[test]
fn late_match_charge_conflict_rolls_back_all_dependency_and_payload_changes() {
    fn setup(
        charge: i32,
    ) -> (
        IdentificationData,
        ObservationId,
        CompoundId,
        ObservationMatchId,
    ) {
        let mut graph = IdentificationData::new().unwrap();
        let file = input(&mut graph, "same-file");
        let obs = observation(&mut graph, file, "same-scan");
        let cmp = compound(&mut graph, "same-compound");
        let mut value = ObservationMatch::new(cmp, obs);
        value.charge = charge;
        let id = graph.register_observation_match(value).unwrap();
        (graph, obs, cmp, id)
    }
    let (mut destination, obs, cmp, id) = setup(2);
    destination
        .set_observation_coordinates(obs, Some(10.0), Some(100.0))
        .unwrap();
    destination
        .set_compound_meta_value(cmp, "old".into(), "preserved".into())
        .unwrap();
    let before = destination.record_count();
    let before_obs = destination.observation(obs).unwrap().clone();
    let before_cmp = destination.compound(cmp).unwrap().clone();
    let before_match = destination.observation_match(id).unwrap().clone();
    let (mut source, source_obs, source_cmp, _) = setup(3);
    source
        .set_observation_meta_value(source_obs, "new".into(), 1.into())
        .unwrap();
    source
        .set_compound_meta_value(source_cmp, "old".into(), "changed".into())
        .unwrap();
    input(&mut source, "new-file");
    compound(&mut source, "new-compound");
    source
        .register_adduct(AdductInfo::new("Na+", formula("Na"), 1, 1).unwrap())
        .unwrap();
    let error = destination.merge_from(&source).unwrap_err();
    assert!(error.to_string().contains("charges"));
    assert_eq!(destination.record_count(), before);
    assert_eq!(destination.observation(obs).unwrap(), &before_obs);
    assert_eq!(destination.compound(cmp).unwrap(), &before_cmp);
    assert_eq!(destination.observation_match(id).unwrap(), &before_match);
    assert_eq!(destination.adduct_count(), 0);
}

#[test]
fn molecule_dispatch_uses_separate_source_fragment_enums_and_stored_compound_formula() {
    let mut graph = IdentificationData::new().unwrap();
    let peptide = graph
        .register_identified_peptide(IdentifiedPeptide::new(
            AASequence::parse(".(Acetyl)AC").unwrap(),
        ))
        .unwrap();
    assert_eq!(graph.molecule_string(peptide).unwrap(), ".(Acetyl)AC");
    assert_eq!(
        graph
            .molecule_formula(
                peptide,
                MoleculeFormulaKind::Peptide(PeptideFragmentType::BIon),
                2
            )
            .unwrap(),
        formula("C8H12N2O3S").with_charge(2)
    );
    assert_eq!(
        graph
            .molecule_formula(
                peptide,
                MoleculeFormulaKind::Peptide(PeptideFragmentType::Internal),
                -1
            )
            .unwrap(),
        formula("C6H10N2O2S").with_charge(-1)
    );
    assert!(
        graph
            .molecule_formula(peptide, MoleculeFormulaKind::Rna(NAFragmentType::Full), 0)
            .is_err()
    );
    let rna = graph
        .register_identified_oligo(IdentifiedOligo::new(NASequence::parse("AC").unwrap()))
        .unwrap();
    assert_eq!(
        graph.molecule_full_formula(rna, -1).unwrap(),
        formula("C19H24N8O11P")
    );
    assert_eq!(graph.molecule_string(rna).unwrap(), "AC");
    let mut value = IdentifiedCompound::new("compound-id");
    value.formula = formula("H-1C2").with_charge(-3);
    value.name = "not-the-ID".into();
    let expected = value.formula.clone();
    let cmp = graph.register_identified_compound(value).unwrap();
    assert_eq!(graph.molecule_string(cmp).unwrap(), "compound-id");
    assert_eq!(
        graph
            .molecule_formula(
                cmp,
                MoleculeFormulaKind::Rna(NAFragmentType::AIon),
                i32::MAX
            )
            .unwrap(),
        expected
    );
    let unknown = graph
        .register_identified_peptide(IdentifiedPeptide::new(
            AASequence::parse("X[123.45]").unwrap(),
        ))
        .unwrap();
    assert!(graph.molecule_full_formula(unknown, 0).is_err());
}

#[test]
fn record_edge_and_discarded_payload_limits_fail_atomically() {
    let mut graph = IdentificationData::with_limits(GraphLimits {
        max_records: 3,
        max_edges: 1,
        ..GraphLimits::default()
    })
    .unwrap();
    let file = input(&mut graph, "f");
    let obs = observation(&mut graph, file, "o");
    let cmp = compound(&mut graph, "c");
    assert!(
        graph
            .register_observation_match(ObservationMatch::new(cmp, obs))
            .is_err()
    );
    assert_eq!(graph.record_count(), 3);
    assert_eq!(graph.observation_match_count(), 0);
    let mut graph = IdentificationData::with_limits(GraphLimits {
        max_work: 100_000,
        ..GraphLimits::default()
    })
    .unwrap();
    let cmp = compound(&mut graph, "retained");
    let mut discarded = IdentifiedCompound::new("retained");
    discarded.smile = "x".repeat(100_001);
    assert!(graph.register_identified_compound(discarded).is_err());
    assert_eq!(graph.compound(cmp).unwrap().smile, "");
    let file = input(&mut graph, "f");
    let obs = observation(&mut graph, file, "o");
    let mut first = ObservationMatch::new(cmp, obs);
    first.peak_annotations.insert(None, Vec::new());
    let id = graph.register_observation_match(first.clone()).unwrap();
    let mut discarded = first.clone();
    discarded.peak_annotations.insert(
        None,
        vec![PeakAnnotation {
            annotation: "x".repeat(100_001),
            ..PeakAnnotation::default()
        }],
    );
    assert!(graph.register_observation_match(discarded).is_err());
    assert_eq!(graph.observation_match(id).unwrap(), &first);
}
