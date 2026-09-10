// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

// Source: Core SDK 6bfc0e4711105f4eda2fea86812a83af7c7e791f,
// METADATA/ID/{IdentifiedMolecule,ObservationMatch}.h and IdentificationData_test.cpp.
// Literal fixture names/formulas are from IdentifiedMolecule_test.cpp:36–38;
// merge expectations below are independent deductions from ObservationMatch::merge.
use openms::chemistry::{AASequence, AdductInfo, EmpiricalFormula, NASequence};
use openms::identification::PeakAnnotation;
use openms::identification::graph::*;

struct Fixture {
    graph: IdentificationData,
    peptide: PeptideId,
    compound: CompoundId,
    oligo: OligoId,
    observation: ObservationId,
    adduct: AdductId,
    score: ScoreTypeId,
    step: ProcessingStepId,
    later_step: ProcessingStepId,
}
fn fixture() -> Fixture {
    let mut graph = IdentificationData::new().unwrap();
    let file = graph
        .register_input_file(InputFile::new("test.mzML"))
        .unwrap();
    let mut observation = Observation::new("spectrum_1", file);
    observation.rt = Some(100.0);
    observation.mz = Some(1000.0);
    let observation = graph.register_observation(observation).unwrap();
    let peptide = graph
        .register_identified_peptide(IdentifiedPeptide::new(
            AASequence::parse("PEPTIDE").unwrap(),
        ))
        .unwrap();
    let mut compound = IdentifiedCompound::new("comp_id");
    compound.formula = EmpiricalFormula::parse("C6H12O6").unwrap();
    compound.name = "glucose".into();
    let compound = graph.register_identified_compound(compound).unwrap();
    let oligo = graph
        .register_identified_oligo(IdentifiedOligo::new(NASequence::parse("ACGU").unwrap()))
        .unwrap();
    let adduct = graph
        .register_adduct(AdductInfo::new("Na+", "Na".parse().unwrap(), 1, 1).unwrap())
        .unwrap();
    let score = graph
        .register_score_type(ScoreType::new("test_score", true))
        .unwrap();
    let software = graph
        .register_processing_software(ProcessingSoftware::new("Tool", "1.0"))
        .unwrap();
    let step = graph
        .register_processing_step(ProcessingStep::new(software), None)
        .unwrap();
    let software = graph
        .register_processing_software(ProcessingSoftware::new("Tool", "2.0"))
        .unwrap();
    let later_step = graph
        .register_processing_step(ProcessingStep::new(software), None)
        .unwrap();
    Fixture {
        graph,
        peptide,
        compound,
        oligo,
        observation,
        adduct,
        score,
        step,
        later_step,
    }
}
fn annotation(text: &str) -> PeakAnnotation {
    PeakAnnotation {
        mz: 123.5,
        intensity: 17.0,
        charge: -2,
        annotation: text.into(),
    }
}

#[test]
fn source_molecule_variants_and_checked_getters_use_reference_identity() {
    let mut f = fixture();
    let peptide = IdentifiedMolecule::from(f.peptide);
    let compound = IdentifiedMolecule::from(f.compound);
    let oligo = IdentifiedMolecule::from(f.oligo);
    assert_eq!(peptide.molecule_type(), MoleculeType::Protein);
    assert_eq!(compound.molecule_type(), MoleculeType::Compound);
    assert_eq!(oligo.molecule_type(), MoleculeType::RNA);
    assert_eq!(peptide.peptide().unwrap(), f.peptide);
    assert_eq!(compound.compound().unwrap(), f.compound);
    assert_eq!(oligo.oligo().unwrap(), f.oligo);
    assert!(peptide.compound().is_err());
    assert!(peptide.oligo().is_err());
    assert!(compound.peptide().is_err());
    assert!(compound.oligo().is_err());
    assert!(oligo.peptide().is_err());
    assert!(oligo.compound().is_err());
    assert!(peptide < compound && compound < oligo);
    let another = f
        .graph
        .register_identified_peptide(IdentifiedPeptide::new(AASequence::parse("AAAA").unwrap()))
        .unwrap();
    // Native stable slot order deliberately replaces source allocation addresses.
    assert!(peptide < IdentifiedMolecule::from(another));
    assert!("AAAA" < "PEPTIDE");
    let mut other_graph = IdentificationData::new().unwrap();
    let same_text = other_graph
        .register_identified_peptide(IdentifiedPeptide::new(
            AASequence::parse("PEPTIDE").unwrap(),
        ))
        .unwrap();
    assert_ne!(peptide, IdentifiedMolecule::from(same_text));
}

#[test]
fn source_match_defaults_and_literal_charge_adduct_are_independent() {
    let mut f = fixture();
    let default = ObservationMatch::new(f.peptide, f.observation);
    assert_eq!(default.charge, 0);
    assert_eq!(default.adduct, None);
    assert!(default.peak_annotations.is_empty());
    assert!(default.result.steps_and_scores.is_empty());
    let mut oligo_match = ObservationMatch::new(f.oligo, f.observation);
    oligo_match.charge = 2;
    oligo_match.adduct = Some(f.adduct);
    let id = f
        .graph
        .register_observation_match(oligo_match.clone())
        .unwrap();
    assert_eq!(f.graph.observation_match(id).unwrap(), &oligo_match);
    assert_eq!(f.graph.adduct(f.adduct).unwrap().charge(), 1);
    assert_eq!(f.graph.adduct(f.adduct).unwrap().name(), "Na+");
}

#[test]
fn source_annotation_map_keeps_first_whole_vector_including_empty_and_duplicates() {
    let f = fixture();
    let mut left = ObservationMatch::new(f.peptide, f.observation);
    let repeated = annotation("b3");
    left.peak_annotations
        .insert(None, vec![repeated.clone(), repeated]);
    left.peak_annotations.insert(Some(f.step), Vec::new());
    let saved = left.peak_annotations.clone();
    let mut incoming = ObservationMatch::new(f.peptide, f.observation);
    incoming
        .peak_annotations
        .insert(None, vec![annotation("replacement")]);
    incoming
        .peak_annotations
        .insert(Some(f.step), vec![annotation("ignored")]);
    let later = vec![annotation("z9"), annotation("a1")];
    incoming
        .peak_annotations
        .insert(Some(f.later_step), later.clone());
    left.merge(&incoming).unwrap();
    assert_eq!(left.peak_annotations[&None], saved[&None]);
    assert!(left.peak_annotations[&Some(f.step)].is_empty());
    assert_eq!(left.peak_annotations[&Some(f.later_step)], later);
    assert_eq!(
        left.peak_annotations.keys().copied().collect::<Vec<_>>(),
        vec![None, Some(f.step), Some(f.later_step)]
    );
}

#[test]
fn source_match_merge_updates_scores_without_reordering_history() {
    let f = fixture();
    let mut left = ObservationMatch::new(f.peptide, f.observation);
    left.result.add_score(f.score, 100.0, Some(f.step)).unwrap();
    left.result
        .add_score(f.score, 200.0, Some(f.later_step))
        .unwrap();
    left.result.metadata.insert("origin".into(), "first".into());
    let mut incoming = ObservationMatch::new(f.peptide, f.observation);
    incoming
        .result
        .add_score(f.score, 300.0, Some(f.step))
        .unwrap();
    incoming
        .result
        .metadata
        .insert("origin".into(), "incoming".into());
    incoming.result.metadata.insert("new".into(), 1.into());
    left.merge(&incoming).unwrap();
    assert_eq!(
        left.result.score_at_step(f.score, Some(f.step)),
        Some(300.0)
    );
    assert_eq!(left.result.score(f.score), Some(200.0));
    assert_eq!(
        left.result
            .steps_and_scores
            .iter()
            .map(|s| s.processing_step)
            .collect::<Vec<_>>(),
        vec![Some(f.step), Some(f.later_step)]
    );
    assert_eq!(left.result.metadata["origin"].as_str().unwrap(), "incoming");
    assert_eq!(left.result.metadata.len(), 2);
}

#[test]
fn source_charge_merge_is_asymmetric_and_conflicts_are_atomic() {
    let f = fixture();
    for charge in [0, 2, -2, i32::MIN, i32::MAX] {
        let mut left = ObservationMatch::new(f.peptide, f.observation);
        let mut incoming = left.clone();
        incoming.charge = charge;
        left.merge(&incoming).unwrap();
        assert_eq!(left.charge, charge);
        left.merge(&incoming).unwrap();
        if charge != 0 {
            let saved = left.clone();
            incoming.charge = 0;
            incoming.result.add_score(f.score, 100.0, None).unwrap();
            incoming.result.metadata.insert("late".into(), 1.into());
            incoming
                .peak_annotations
                .insert(None, vec![annotation("late")]);
            assert!(left.merge(&incoming).is_err());
            assert_eq!(left, saved);
        }
    }
    let mut left = ObservationMatch::new(f.peptide, f.observation);
    left.charge = 2;
    let mut incoming = left.clone();
    incoming.charge = 3;
    let saved = left.clone();
    assert!(left.merge(&incoming).is_err());
    assert_eq!(left, saved);
}

#[test]
fn direct_record_adduct_merge_fills_only_absent_adduct_and_preserves_identity_fields() {
    let f = fixture();
    let mut left = ObservationMatch::new(f.peptide, f.observation);
    let mut incoming = ObservationMatch::new(f.oligo, f.observation);
    incoming.adduct = Some(f.adduct);
    left.merge(&incoming).unwrap();
    assert_eq!(
        left.identified_molecule,
        IdentifiedMolecule::from(f.peptide)
    );
    assert_eq!(left.observation, f.observation);
    assert_eq!(left.adduct, Some(f.adduct));
    let saved = left.clone();
    incoming.adduct = None;
    incoming
        .result
        .metadata
        .insert("late".into(), "discard".into());
    assert!(left.merge(&incoming).is_err());
    assert_eq!(left, saved);
}

#[test]
fn nonfinite_annotation_in_ignored_incoming_vector_is_checked_before_mutation() {
    let f = fixture();
    for (mz, intensity) in [(f64::NAN, 1.0), (1.0, f64::INFINITY)] {
        let mut left = ObservationMatch::new(f.peptide, f.observation);
        left.peak_annotations.insert(None, Vec::new());
        let saved = left.clone();
        let mut incoming = left.clone();
        incoming
            .result
            .metadata
            .insert("late".into(), "discard".into());
        incoming.peak_annotations.insert(
            None,
            vec![PeakAnnotation {
                mz,
                intensity,
                ..annotation("bad")
            }],
        );
        assert!(left.merge(&incoming).is_err());
        assert_eq!(left, saved);
    }
}

#[test]
fn full_record_equality_includes_payload_beyond_source_registration_keys() {
    let f = fixture();
    let original = ObservationMatch::new(f.peptide, f.observation);
    let mut changed = original.clone();
    changed.charge = 2;
    assert_ne!(original, changed);
    changed = original.clone();
    changed.peak_annotations.insert(None, Vec::new());
    assert_ne!(original, changed);
    changed = original.clone();
    changed.result.add_score(f.score, 1.0, None).unwrap();
    assert_ne!(original, changed);
}
