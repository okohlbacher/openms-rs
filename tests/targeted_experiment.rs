// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `TargetedExperimentHelper_test.cpp`, `ReactionMonitoringTransition_test.cpp`,
//! `IncludeExcludeTarget_test.cpp` and `TargetedExperiment_test.cpp` at core pin
//! `bc9cc12`, plus native checks. Literals are transcribed from those class tests
//! (evidence tier 3, source review) unless a test says otherwise. No C++ was run
//! for these checks.

use openms::Error;
use openms::analysis::openswath::include_exclude_target::IncludeExcludeTarget;
use openms::analysis::openswath::reaction_monitoring_transition::{
    DecoyTransitionType, Product, ReactionMonitoringTransition,
};
use openms::analysis::openswath::targeted_experiment::TargetedExperiment;
use openms::analysis::openswath::targeted_experiment_helper::{
    CV, Compound, Configuration, Contact, Instrument, Interpretation, Peptide, PeptideCompound,
    PeptideModification, Prediction, Protein, Publication, RTType, RTUnit, RetentionTime,
    TraMLProduct, UnresolvedModificationPolicy, get_aa_sequence, get_aa_sequence_with,
    set_modification,
};
use openms::chemistry::{AASequence, ModificationsDB, PeptideFragmentType};
use openms::kernel::MSChromatogram;
use openms::kernel::mrm::MRMTransitionGroup;
use openms::metadata::{CVTerm, CVTermList, MetaValue, Software, SourceFile};
use std::collections::hash_map::DefaultHasher;
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};

fn hash_of<T: Hash>(value: &T) -> u64 {
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

/// The shared `charge_cv` fixture of `ReactionMonitoringTransition_test.cpp`.
fn charge_cv(value: i32) -> CVTerm {
    let mut term = CVTerm::new("MS:1000041", "charge state", "MS");
    term.value = MetaValue::from(value);
    term
}

fn predicted_rt(rt: f64) -> RetentionTime {
    let mut value = RetentionTime::new();
    value.set_rt(rt);
    value.retention_time_unit = RTUnit::Second;
    value.retention_time_type = RTType::Predicted;
    value
}

// ---------------------------------------------------------------------------
// TargetedExperimentHelper_test.cpp
// ---------------------------------------------------------------------------

/// START_SECTION((TargetedExperiment::RetentionTime))
#[test]
fn helper_retention_time() {
    let mut rt = RetentionTime::new();
    assert!(!rt.is_rt_set());
    assert_eq!(rt.rt(), None);
    rt.set_rt(5.0);
    assert!(rt.is_rt_set());
    assert_eq!(rt.rt(), Some(5.0));
    // Header defaults: the SIZE_OF enumerators, an empty software reference.
    let fresh = RetentionTime::default();
    assert_eq!(fresh.retention_time_unit, RTUnit::SizeOfRTUnit);
    assert_eq!(fresh.retention_time_type, RTType::SizeOfRTType);
    assert_eq!(fresh.software_ref, "");
    assert_eq!(RTUnit::SizeOfRTUnit as i8, 3);
    assert_eq!(RTType::SizeOfRTType as i8, 6);
    assert_eq!(RTType::Irt as i8, 4);
}

/// START_SECTION((TargetedExperiment::Peptide))
#[test]
fn helper_peptide() {
    let mut p = Peptide::default();
    assert!(!p.has_retention_time());
    assert_eq!(p.rts.len(), 0);

    let rt = predicted_rt(5.0);
    p.rts.push(rt.clone());
    assert_eq!(p.rts.len(), 1);
    assert_eq!(p.rts[0], rt);
    assert_eq!(p.rts[0].retention_time_unit, RTUnit::Second);
    assert_eq!(p.rts[0].retention_time_type, RTType::Predicted);
    assert_eq!(p.rts[0].rt(), Some(5.0));

    assert!(p.has_retention_time());
    assert_eq!(p.retention_time().unwrap(), 5.0);
    assert_eq!(p.retention_time_unit().unwrap(), RTUnit::Second);
    assert_eq!(p.retention_time_type().unwrap(), RTType::Predicted);

    assert_eq!(p.peptide_group_label(), "");
    p.set_peptide_group_label("test1");
    assert_eq!(p.peptide_group_label(), "test1");

    assert!(!p.has_charge());
    p.set_charge_state(-1);
    assert_eq!(p.charge_state(), Some(-1));
    p.set_charge_state(2);
    assert_eq!(p.charge_state(), Some(2));
}

/// START_SECTION((TargetedExperiment::Compound))
#[test]
fn helper_compound() {
    let mut p = Compound::default();
    assert!(!p.has_retention_time());
    assert_eq!(p.rts.len(), 0);
    let rt = predicted_rt(5.0);
    p.rts.push(rt.clone());
    assert_eq!(p.rts.len(), 1);
    assert_eq!(p.rts[0], rt);
    assert!(p.has_retention_time());
    assert_eq!(p.retention_time().unwrap(), 5.0);
    assert_eq!(p.retention_time_unit().unwrap(), RTUnit::Second);
    assert_eq!(p.retention_time_type().unwrap(), RTType::Predicted);
    assert!(!p.has_charge());
    p.set_charge_state(-1);
    assert_eq!(p.charge_state(), Some(-1));
    p.set_charge_state(2);
    assert_eq!(p.charge_state(), Some(2));
    assert_eq!(p.theoretical_mass, 0.0);
}

/// Native: the throwing retention-time getters (source
/// `Exception::IllegalArgument`), and an entry without a value.
#[test]
fn helper_retention_time_getters_refuse_missing_value() {
    let mut p = PeptideCompound::new();
    assert!(matches!(p.retention_time(), Err(Error::InvalidValue(_))));
    assert!(matches!(
        p.retention_time_type(),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        p.retention_time_unit(),
        Err(Error::InvalidValue(_))
    ));
    // An entry that exists but has no value set does not count.
    p.rts.push(RetentionTime::new());
    assert!(!p.has_retention_time());
    assert!(p.retention_time().is_err());
    // Drift time defaults to -1 and is excluded from equality, as the source.
    assert_eq!(p.drift_time(), -1.0);
    let mut q = p.clone();
    q.set_drift_time(1.5);
    assert_eq!(q.drift_time(), 1.5);
    assert_eq!(p, q);
    assert_eq!(hash_of(&p), hash_of(&q));
}

fn assert_hash_contract<T: Hash + Eq + Clone + std::fmt::Debug>(a: T, b: T, c: T) {
    assert_eq!(a, b);
    assert_eq!(hash_of(&a), hash_of(&b));
    assert_ne!(a, c);
    assert_ne!(hash_of(&a), hash_of(&c));
    let set: HashSet<T> = [a, b, c].into_iter().collect();
    assert_eq!(set.len(), 2);
}

/// START_SECTION((std::hash<TargetedExperimentHelper::CV>)) and the eleven
/// sibling hash sections: equal values hash equally, the source's "different"
/// value hashes differently, and a set of the three holds two.
#[test]
fn helper_hashes() {
    assert_hash_contract(
        CV::new("id1", "fullname1", "version1", "URI1"),
        CV::new("id1", "fullname1", "version1", "URI1"),
        CV::new("id2", "fullname2", "version2", "URI2"),
    );

    let protein = |id: &str, seq: &str| Protein {
        id: id.into(),
        sequence: seq.into(),
        ..Protein::default()
    };
    assert_hash_contract(
        protein("protein1", "ACDEFGH"),
        protein("protein1", "ACDEFGH"),
        protein("protein2", "IKJLMNPQRST"),
    );

    let rt = |value: f64, software: &str, unit, kind| {
        let mut rt = RetentionTime::new();
        rt.set_rt(value);
        rt.software_ref = software.into();
        rt.retention_time_unit = unit;
        rt.retention_time_type = kind;
        rt
    };
    assert_hash_contract(
        rt(5.0, "software1", RTUnit::Second, RTType::Predicted),
        rt(5.0, "software1", RTUnit::Second, RTType::Predicted),
        rt(10.0, "software2", RTUnit::Minute, RTType::Local),
    );

    let pc = |id: &str, charge: i32| {
        let mut pc = PeptideCompound::new();
        pc.id = id.into();
        pc.set_charge_state(charge);
        pc
    };
    assert_hash_contract(pc("compound1", 2), pc("compound1", 2), pc("compound2", 3));

    let compound = |id: &str, formula: &str, smiles: &str, mass: f64| {
        let mut c = Compound::default();
        c.id = id.into();
        c.molecular_formula = formula.into();
        c.smiles_string = smiles.into();
        c.theoretical_mass = mass;
        c
    };
    assert_hash_contract(
        compound(
            "compound1",
            "C6H12O6",
            "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
            180.063,
        ),
        compound(
            "compound1",
            "C6H12O6",
            "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
            180.063,
        ),
        compound("compound2", "H2O", "O", 18.015),
    );

    let peptide = |id: &str, seq: &str, label: &str| {
        let mut p = Peptide::default();
        p.id = id.into();
        p.sequence = seq.into();
        p.set_peptide_group_label(label);
        p
    };
    assert_hash_contract(
        peptide("peptide1", "ACDEFGHIK", "group1"),
        peptide("peptide1", "ACDEFGHIK", "group1"),
        peptide("peptide2", "LMNPQRSTVWY", "group2"),
    );

    let contact = |id: &str| Contact {
        id: id.into(),
        ..Contact::default()
    };
    assert_hash_contract(
        contact("contact1"),
        contact("contact1"),
        contact("contact2"),
    );
    let publication = |id: &str| Publication {
        id: id.into(),
        ..Publication::default()
    };
    assert_hash_contract(
        publication("pub1"),
        publication("pub1"),
        publication("pub2"),
    );
    let instrument = |id: &str| Instrument {
        id: id.into(),
        ..Instrument::default()
    };
    assert_hash_contract(
        instrument("inst1"),
        instrument("inst1"),
        instrument("inst2"),
    );
    let prediction = |software: &str, contact: &str| Prediction {
        software_ref: software.into(),
        contact_ref: contact.into(),
        ..Prediction::default()
    };
    assert_hash_contract(
        prediction("software1", "contact1"),
        prediction("software1", "contact1"),
        prediction("software2", "contact2"),
    );

    let interpretation = |ordinal, rank, iontype| Interpretation {
        ordinal,
        rank,
        iontype,
        ..Interpretation::default()
    };
    assert_hash_contract(
        interpretation(1, 1, PeptideFragmentType::BIon),
        interpretation(1, 1, PeptideFragmentType::BIon),
        interpretation(2, 2, PeptideFragmentType::YIon),
    );

    let product = |mz: f64, charge: i32| {
        let mut p = TraMLProduct::new();
        p.set_mz(mz);
        p.set_charge_state(charge);
        p
    };
    assert_hash_contract(product(500.5, 2), product(500.5, 2), product(600.6, 3));
}

/// Native: the TraMLProduct accessors and defaults the header states.
#[test]
fn helper_traml_product_accessors() {
    let mut p = TraMLProduct::new();
    assert!(!p.has_charge());
    assert_eq!(p.charge_state(), None);
    assert_eq!(p.mz(), 0.0);
    assert!(p.configuration_list().is_empty());
    p.add_configuration(Configuration::default());
    p.add_interpretation(Interpretation::default());
    assert_eq!(p.configuration_list().len(), 1);
    assert_eq!(p.interpretation_list().len(), 1);
    assert_eq!(
        p.interpretation_list()[0].iontype,
        PeptideFragmentType::Unannotated
    );
    p.reset_interpretations();
    assert!(p.interpretation_list().is_empty());
    assert_eq!(p.configuration_list().len(), 1);
}

/// Native divergence, documented at `Configuration` and `PeptideModification`:
/// the source compares only their CV terms; this port compares every field.
#[test]
fn helper_configuration_and_modification_equality_cover_all_fields() {
    let a = Configuration::default();
    let b = Configuration {
        instrument_ref: "other".into(),
        ..Configuration::default()
    };
    assert_ne!(a, b);
    let m = PeptideModification::default();
    assert_eq!(m.location, -1);
    assert_eq!(m.unimod_id, -1);
    let n = PeptideModification {
        location: 3,
        ..PeptideModification::default()
    };
    assert_ne!(m, n);
    // Both signs of zero are equal and hash equally.
    let z1 = PeptideModification {
        mono_mass_delta: 0.0,
        ..PeptideModification::default()
    };
    let z2 = PeptideModification {
        mono_mass_delta: -0.0,
        ..PeptideModification::default()
    };
    assert_eq!(z1, z2);
    assert_eq!(hash_of(&z1), hash_of(&z2));
}

fn peptide_with_mod(sequence: &str, location: i32, avg: f64, mono: f64) -> Peptide {
    let mut peptide = Peptide::default();
    peptide.sequence = sequence.into();
    peptide.mods.push(PeptideModification {
        avg_mass_delta: avg,
        location,
        mono_mass_delta: mono,
        ..PeptideModification::default()
    });
    peptide
}

/// START_SECTION(OpenMS::AASequence getAASequence(const OpenMS::TargetedExperiment::Peptide &peptide))
/// (in `TargetedExperiment_test.cpp`). The source asserts only the unmodified
/// string; the modified strings it leaves commented out are asserted here as
/// well, as tier 4: each mass is the Unimod monoisotopic delta of the named
/// modification, so the closest match within 1 Da is independently that one.
#[test]
fn helper_get_aa_sequence() {
    let peptide = peptide_with_mod("TESTPEPTIDE", 2, 79.9799, 79.966331);
    let aas = get_aa_sequence(&peptide).unwrap();
    assert_eq!(aas.as_str(), peptide.sequence);
    assert_eq!(aas.to_string(), "TES(Phospho)TPEPTIDE");

    let peptide2 = peptide_with_mod("TESTPEPTIDER", 11, 9.9296, 10.008269);
    let aas2 = get_aa_sequence(&peptide2).unwrap();
    assert_eq!(aas2.as_str(), peptide2.sequence);
    assert_eq!(aas2.to_string(), "TESTPEPTIDER(Label:13C(6)15N(4))");

    let peptide3 = peptide_with_mod("TESTMPEPTIDE", 4, 15.9994, 15.994915);
    let aas3 = get_aa_sequence(&peptide3).unwrap();
    assert_eq!(aas3.as_str(), peptide3.sequence);
    assert_eq!(aas3.to_string(), "TESTM(Oxidation)PEPTIDE");
}

/// Native: the UniMod-id path, the terminal locations of `setModification`,
/// the refused out-of-range location, and both unresolved-modification
/// policies.
#[test]
fn helper_get_aa_sequence_paths_and_policies() {
    // UniMod:35 is Oxidation.
    let mut peptide = peptide_with_mod("PEPMTIDE", 3, 0.0, 0.0);
    peptide.mods[0].unimod_id = 35;
    let db = ModificationsDB::global();
    let converted =
        get_aa_sequence_with(&peptide, db, UnresolvedModificationPolicy::Reject).unwrap();
    assert_eq!(converted.sequence.to_string(), "PEPM(Oxidation)TIDE");
    assert!(converted.warnings.is_empty());

    // Mass lookup records the source's warning.
    let by_mass = peptide_with_mod("PEPMTIDE", 3, 15.9994, 15.994915);
    let converted =
        get_aa_sequence_with(&by_mass, db, UnresolvedModificationPolicy::Reject).unwrap();
    assert_eq!(converted.warnings.len(), 1);
    assert!(converted.warnings[0].contains("No UniMod id set"));

    // A mass nothing matches: rejected natively, skipped under SourceSkip.
    let unresolved = peptide_with_mod("PEPMTIDE", 3, 0.0, 5000.0);
    assert!(matches!(
        get_aa_sequence(&unresolved),
        Err(Error::InvalidValue(_))
    ));
    let skipped =
        get_aa_sequence_with(&unresolved, db, UnresolvedModificationPolicy::SourceSkip).unwrap();
    assert_eq!(skipped.sequence.to_string(), "PEPMTIDE");
    assert_eq!(skipped.warnings.len(), 3);
    assert_eq!(skipped.warnings[2], "Skipping this modification");

    // Terminal location without a UniMod id has no residue to match.
    let terminal = peptide_with_mod("PEPTIDE", -1, 42.0106, 42.010565);
    assert!(get_aa_sequence(&terminal).is_err());

    // setModification: N-terminus, C-terminus, residue, out of range.
    let mut aas = AASequence::parse("PEPTIDEK").unwrap();
    set_modification(-1, 8, "Acetyl", &mut aas).unwrap();
    set_modification(8, 8, "Amidated", &mut aas).unwrap();
    set_modification(7, 8, "Label:13C(6)15N(2)", &mut aas).unwrap();
    assert_eq!(aas.as_str(), "PEPTIDEK");
    assert!(aas.to_string().starts_with(".(Acetyl)"));
    let before = aas.clone();
    assert!(matches!(
        set_modification(9, 8, "Acetyl", &mut aas),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        set_modification(-2, 8, "Acetyl", &mut aas),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(aas, before);
}

// ---------------------------------------------------------------------------
// ReactionMonitoringTransition_test.cpp
// ---------------------------------------------------------------------------

/// Defaults: getName/getPeptideRef/getCompoundRef empty, getPrecursorMZ and
/// getProductMZ 0, getProduct == Product(), the three flag getters, and the
/// constructor's library intensity (-101) and decoy type (UNKNOWN) from the
/// `.cpp`.
#[test]
fn rmt_defaults() {
    let transition = ReactionMonitoringTransition::new();
    assert_eq!(transition.name(), "");
    assert_eq!(transition.native_id(), "");
    assert_eq!(transition.peptide_ref(), "");
    assert_eq!(transition.compound_ref(), "");
    assert_eq!(transition.precursor_mz(), 0.0);
    assert_eq!(transition.product_mz(), 0.0);
    assert_eq!(*transition.product(), Product::default());
    assert!(transition.is_detecting_transition());
    assert!(!transition.is_identifying_transition());
    assert!(transition.is_quantifying_transition());
    assert_eq!(transition.library_intensity(), -101.0);
    assert_eq!(
        transition.decoy_transition_type(),
        DecoyTransitionType::Unknown
    );
    assert!(!transition.has_precursor_cv_terms());
    assert!(!transition.has_prediction());
    assert!(transition.intermediate_products().is_empty());
    assert!(!transition.is_product_charge_state_set());
    assert_eq!(transition.product_charge_state(), None);
}

/// START_SECTION((ReactionMonitoringTransition(const ReactionMonitoringTransition &rhs)))
#[test]
fn rmt_copy() {
    let mut tr1 = ReactionMonitoringTransition::new();
    tr1.add_precursor_cv_term(charge_cv(3)).unwrap();
    tr1.set_precursor_mz(42.0);
    let tr2 = tr1.clone();
    assert_eq!(tr1, tr2);
    let p = Prediction {
        contact_ref: "dummy".into(),
        ..Prediction::default()
    };
    tr1.set_prediction(p);
    tr1.set_identifying_transition(false);
    tr1.set_detecting_transition(false);
    tr1.set_quantifying_transition(false);
    let tr3 = tr1.clone();
    assert_eq!(tr1, tr3);
    assert_ne!(tr1, tr2);
}

/// START_SECTION((ReactionMonitoringTransition(ReactionMonitoringTransition &&rhs)))
/// A Rust move leaves nothing behind to inspect; `std::mem::take` stands in
/// for the moved-from source object, whose members the test checks are gone.
#[test]
fn rmt_move() {
    let mut tr1 = ReactionMonitoringTransition::new();
    tr1.set_prediction(Prediction {
        contact_ref: "dummy".into(),
        ..Prediction::default()
    });
    tr1.add_precursor_cv_term(charge_cv(3)).unwrap();
    tr1.set_precursor_mz(42.0);
    tr1.set_compound_ref("test_ref");
    let mut orig = tr1.clone();
    let tr2 = std::mem::take(&mut tr1);
    assert_eq!(orig, tr2);
    assert!(tr2.has_precursor_cv_terms());
    assert!(tr2.has_prediction());
    assert_eq!(tr2.prediction().unwrap().contact_ref, "dummy");
    assert!(tr2.precursor_cv_term_list().unwrap().contains("MS:1000041"));
    assert_eq!(tr2.compound_ref(), "test_ref");
    assert!(!tr1.has_precursor_cv_terms());
    assert!(!tr1.has_prediction());
    assert_eq!(tr1.compound_ref(), "");

    orig.set_prediction(Prediction {
        contact_ref: "dummy".into(),
        ..Prediction::default()
    });
    orig.set_identifying_transition(false);
    orig.set_detecting_transition(false);
    orig.set_quantifying_transition(false);
    let tr3 = orig.clone();
    assert_eq!(orig, tr3);
    assert_ne!(orig, tr2);
}

/// START_SECTION((ReactionMonitoringTransition& operator=(const ReactionMonitoringTransition &rhs)))
#[test]
fn rmt_assignment() {
    let mut tr1 = ReactionMonitoringTransition::new();
    tr1.add_precursor_cv_term(charge_cv(3)).unwrap();
    tr1.set_precursor_mz(42.0);
    let tr2 = tr1.clone();
    assert_eq!(tr1, tr2);
    tr1.set_prediction(Prediction {
        contact_ref: "dummy".into(),
        ..Prediction::default()
    });
    let mut tr3 = tr1.clone();
    assert_eq!(tr1, tr3);
    assert_ne!(tr1, tr2);
    tr1.set_detecting_transition(false);
    assert_ne!(tr1, tr3);
    tr1.set_identifying_transition(true);
    tr1.set_quantifying_transition(false);
    assert_ne!(tr1, tr3);
    tr3 = tr1.clone();
    assert_eq!(tr1, tr3);
}

/// The setter/getter sections: setName, setPeptideRef, setCompoundRef,
/// setPrecursorMZ, setPrecursorCVTermList, hasPrecursorCVTerms,
/// addPrecursorCVTerm, getPrecursorCVTermList, setProductMZ, setProduct,
/// addProductCVTerm and the three flag setters.
#[test]
fn rmt_accessors() {
    let mut tr = ReactionMonitoringTransition::new();
    tr.set_name("test_tr");
    assert_eq!(tr.name(), "test_tr");
    tr.set_peptide_ref("test_ref");
    assert_eq!(tr.peptide_ref(), "test_ref");
    tr.set_compound_ref("test_ref");
    assert_eq!(tr.compound_ref(), "test_ref");
    tr.set_precursor_mz(42.0);
    assert_eq!(tr.precursor_mz(), 42.0);

    let mut list = CVTermList::new();
    list.add(charge_cv(3)).unwrap();
    let mut tr = ReactionMonitoringTransition::new();
    assert!(!tr.has_precursor_cv_terms());
    tr.set_precursor_cv_term_list(list.clone());
    assert!(tr.has_precursor_cv_terms());
    assert!(tr.precursor_cv_term_list().unwrap().contains("MS:1000041"));
    assert_eq!(tr.precursor_cv_term_list().unwrap(), &list);

    let mut tr = ReactionMonitoringTransition::new();
    assert!(!tr.has_precursor_cv_terms());
    tr.add_precursor_cv_term(charge_cv(3)).unwrap();
    assert!(tr.has_precursor_cv_terms());
    assert!(tr.precursor_cv_term_list().unwrap().contains("MS:1000041"));

    tr.set_product_mz(42.0);
    assert_eq!(tr.product_mz(), 42.0);
    let mut product = Product::default();
    product.set_mz(2.0);
    product.set_charge_state(2);
    tr.set_product(product.clone());
    assert_eq!(*tr.product(), product);
    assert!(tr.is_product_charge_state_set());
    assert_eq!(tr.product_charge_state(), Some(2));
    tr.add_product_cv_term(charge_cv(1)).unwrap();
    assert!(tr.product().cv_terms.contains("MS:1000041"));

    tr.set_detecting_transition(false);
    assert!(!tr.is_detecting_transition());
    tr.set_detecting_transition(true);
    assert!(tr.is_detecting_transition());
    tr.set_identifying_transition(true);
    assert!(tr.is_identifying_transition());
    tr.set_identifying_transition(false);
    assert!(!tr.is_identifying_transition());
    tr.set_quantifying_transition(false);
    assert!(!tr.is_quantifying_transition());
    tr.set_quantifying_transition(true);
    assert!(tr.is_quantifying_transition());

    // Native: remaining accessors, and an invalid CV term leaves no list.
    tr.set_library_intensity(7.5);
    assert_eq!(tr.library_intensity(), 7.5);
    tr.set_decoy_transition_type(DecoyTransitionType::Decoy);
    assert_eq!(tr.decoy_transition_type(), DecoyTransitionType::Decoy);
    tr.add_intermediate_product(TraMLProduct::new());
    assert_eq!(tr.intermediate_products().len(), 1);
    tr.set_intermediate_products(Vec::new());
    assert!(tr.intermediate_products().is_empty());
    tr.set_retention_time(predicted_rt(12.0));
    assert_eq!(tr.retention_time().rt(), Some(12.0));
    tr.set_native_id("native");
    assert_eq!(tr.name(), "native");
    let mut fresh = ReactionMonitoringTransition::new();
    assert!(
        fresh
            .add_precursor_cv_term(CVTerm::new("", "bad", "MS"))
            .is_err()
    );
    assert!(!fresh.has_precursor_cv_terms());
    assert!(
        fresh
            .add_prediction_term(CVTerm::new("", "bad", "MS"))
            .is_err()
    );
    assert!(!fresh.has_prediction());
    fresh.add_prediction_term(charge_cv(2)).unwrap();
    assert!(fresh.has_prediction());
}

/// START_SECTION((bool operator==(const ReactionMonitoringTransition &rhs) const ))
/// and operator!=.
#[test]
fn rmt_equality() {
    let mut tr1 = ReactionMonitoringTransition::new();
    let mut tr2 = ReactionMonitoringTransition::new();
    assert_eq!(tr1, tr2);
    tr1.add_precursor_cv_term(charge_cv(3)).unwrap();
    assert_ne!(tr1, tr2);
    tr2.add_precursor_cv_term(charge_cv(3)).unwrap();
    assert_eq!(tr1, tr2);
    tr1.set_detecting_transition(false);
    assert_ne!(tr1, tr2);
    tr2.set_detecting_transition(false);
    assert_eq!(tr1, tr2);
}

/// START_SECTION(([EXTRA] std::hash<ReactionMonitoringTransition>))
#[test]
fn rmt_hash() {
    let mut tr1 = ReactionMonitoringTransition::new();
    let mut tr2 = ReactionMonitoringTransition::new();
    assert_eq!(hash_of(&tr1), hash_of(&tr2));
    tr1.set_name("test_transition");
    tr2.set_name("test_transition");
    assert_eq!(hash_of(&tr1), hash_of(&tr2));
    tr1.set_peptide_ref("peptide_ref_1");
    tr2.set_peptide_ref("peptide_ref_1");
    assert_eq!(hash_of(&tr1), hash_of(&tr2));
    tr1.set_precursor_mz(500.0);
    tr2.set_precursor_mz(500.0);
    assert_eq!(hash_of(&tr1), hash_of(&tr2));
    tr1.set_product_mz(250.0);
    tr2.set_product_mz(250.0);
    assert_eq!(hash_of(&tr1), hash_of(&tr2));
    tr1.add_precursor_cv_term(charge_cv(2)).unwrap();
    tr2.add_precursor_cv_term(charge_cv(2)).unwrap();
    assert_eq!(hash_of(&tr1), hash_of(&tr2));
    let pred = Prediction {
        contact_ref: "contact_1".into(),
        software_ref: "software_1".into(),
        ..Prediction::default()
    };
    tr1.set_prediction(pred.clone());
    tr2.set_prediction(pred);
    assert_eq!(hash_of(&tr1), hash_of(&tr2));
    tr1.set_detecting_transition(false);
    tr2.set_detecting_transition(false);
    assert_eq!(hash_of(&tr1), hash_of(&tr2));
    tr1.set_identifying_transition(true);
    tr2.set_identifying_transition(true);
    assert_eq!(hash_of(&tr1), hash_of(&tr2));

    let mut tr3 = ReactionMonitoringTransition::new();
    tr3.set_name("different_name");
    assert_ne!(hash_of(&tr1), hash_of(&tr3));
    let mut tr4 = ReactionMonitoringTransition::new();
    tr4.set_precursor_mz(999.0);
    assert_ne!(hash_of(&tr1), hash_of(&tr4));

    let mut set = HashSet::new();
    set.insert(tr1.clone());
    set.insert(tr2.clone());
    set.insert(tr3.clone());
    set.insert(tr4.clone());
    assert_eq!(set.len(), 3);

    let mut map = HashMap::new();
    map.insert(tr1.clone(), 1);
    map.insert(tr2.clone(), 2);
    map.insert(tr3.clone(), 3);
    assert_eq!(map.len(), 2);
    assert_eq!(map[&tr1], 2);
    assert_eq!(map[&tr3], 3);

    assert!(set.contains(&tr1));
    assert!(set.contains(&tr2));
    assert!(set.contains(&tr3));
    let mut tr5 = ReactionMonitoringTransition::new();
    tr5.set_name("not_in_set");
    assert!(!set.contains(&tr5));
}

/// START_SECTION((const Prediction& getPrediction() const))
#[test]
fn rmt_get_prediction() {
    let mut tr = ReactionMonitoringTransition::new();
    assert!(!tr.has_prediction());
    assert!(matches!(tr.prediction(), Err(Error::MissingInformation(_))));
    tr.set_prediction(Prediction {
        contact_ref: "dummy".into(),
        ..Prediction::default()
    });
    assert!(tr.has_prediction());
    assert_eq!(tr.prediction().unwrap().contact_ref, "dummy");
}

/// `MRMTransitionGroup_test.cpp` instantiates the group with
/// `ReactionMonitoringTransition`; its `getLibraryIntensity` section
/// (3 and -2 in, 3 and 0 out) now runs with the real transition type.
#[test]
fn rmt_builds_an_mrm_transition_group() {
    let mut group: MRMTransitionGroup<MSChromatogram, ReactionMonitoringTransition> =
        MRMTransitionGroup::new();
    let mut a = ReactionMonitoringTransition::new();
    a.set_library_intensity(3.0);
    let mut b = ReactionMonitoringTransition::new();
    b.set_library_intensity(-2.0);
    group.add_transition(a, "dummy1").unwrap();
    group.add_transition(b, "dummy2").unwrap();
    let result = group.library_intensity().unwrap();
    assert_eq!(result, vec![3.0, 0.0]);
}

// ---------------------------------------------------------------------------
// IncludeExcludeTarget_test.cpp (its accessor sections are TODO upstream)
// ---------------------------------------------------------------------------

/// Defaults from `IncludeExcludeTarget.cpp` and the accessors.
#[test]
fn include_exclude_target_defaults_and_accessors() {
    let mut t = IncludeExcludeTarget::new();
    assert_eq!(t.precursor_mz(), f64::MAX);
    assert_eq!(t.product_mz(), f64::MAX);
    assert_eq!(t.name(), "");
    t.set_name("target");
    t.set_peptide_ref("pep");
    t.set_compound_ref("comp");
    t.set_precursor_mz(500.0);
    t.set_product_mz(250.0);
    t.add_precursor_cv_term(charge_cv(2)).unwrap();
    t.add_product_cv_term(charge_cv(1)).unwrap();
    t.add_prediction_term(charge_cv(1)).unwrap();
    t.add_interpretation(CVTermList::new());
    t.add_configuration(Configuration::default());
    t.set_retention_time(predicted_rt(3.0));
    assert_eq!(t.name(), "target");
    assert_eq!(t.peptide_ref(), "pep");
    assert_eq!(t.compound_ref(), "comp");
    assert_eq!(t.precursor_mz(), 500.0);
    assert_eq!(t.product_mz(), 250.0);
    assert!(t.precursor_cv_term_list().contains("MS:1000041"));
    assert!(t.product_cv_term_list().contains("MS:1000041"));
    assert!(t.prediction().contains("MS:1000041"));
    assert_eq!(t.interpretations().len(), 1);
    assert_eq!(t.configurations().len(), 1);
    assert_eq!(t.retention_time().rt(), Some(3.0));
    t.set_interpretations(Vec::new());
    t.set_configurations(Vec::new());
    t.set_prediction(CVTermList::new());
    t.set_precursor_cv_term_list(CVTermList::new());
    t.set_product_cv_term_list(CVTermList::new());
    assert!(t.interpretations().is_empty());
    assert!(t.configurations().is_empty());
    assert!(t.prediction().is_empty());
}

/// START_SECTION(([EXTRA] std::hash<IncludeExcludeTarget>))
#[test]
fn include_exclude_target_hash() {
    let target = |name: &str, pre: f64, pro: f64, pep: &str, comp: &str| {
        let mut t = IncludeExcludeTarget::new();
        t.set_name(name);
        t.set_precursor_mz(pre);
        t.set_product_mz(pro);
        t.set_peptide_ref(pep);
        t.set_compound_ref(comp);
        t
    };
    let t1 = target("target1", 500.0, 200.0, "peptide_ref_1", "compound_ref_1");
    let t2 = target("target1", 500.0, 200.0, "peptide_ref_1", "compound_ref_1");
    assert_eq!(t1, t2);
    assert_eq!(hash_of(&t1), hash_of(&t2));
    let t3 = target("target2", 500.0, 200.0, "peptide_ref_1", "compound_ref_1");
    assert_ne!(hash_of(&t1), hash_of(&t3));
    let t4 = target("target1", 600.0, 200.0, "peptide_ref_1", "compound_ref_1");
    assert_ne!(hash_of(&t1), hash_of(&t4));
    let mut set = HashSet::new();
    set.insert(t1.clone());
    assert_eq!(set.len(), 1);
    set.insert(t2.clone());
    assert_eq!(set.len(), 1);
    set.insert(t3.clone());
    assert_eq!(set.len(), 2);
    let mut map = HashMap::new();
    map.insert(t1.clone(), 42);
    assert_eq!(map[&t1], 42);
    assert_eq!(map[&t2], 42);
    map.insert(t3.clone(), 99);
    assert_eq!(map[&t3], 99);
    assert_eq!(map.len(), 2);
    let mut t5 = IncludeExcludeTarget::new();
    let mut t6 = IncludeExcludeTarget::new();
    let mut rt = RetentionTime::new();
    rt.set_rt(100.0);
    t5.set_retention_time(rt.clone());
    t6.set_retention_time(rt);
    assert_eq!(hash_of(&t5), hash_of(&t6));
    let mut t7 = IncludeExcludeTarget::new();
    let mut rt2 = RetentionTime::new();
    rt2.set_rt(200.0);
    t7.set_retention_time(rt2);
    assert_ne!(hash_of(&t5), hash_of(&t7));
}

// ---------------------------------------------------------------------------
// TargetedExperiment_test.cpp
// ---------------------------------------------------------------------------

/// START_SECTION((TargetedExperiment(const TargetedExperiment &rhs))),
/// operator==, operator= and the list setter/getter/adder sections.
#[test]
fn te_lists_copy_and_equality() {
    let mut te = TargetedExperiment::new();
    assert!(te.proteins().is_empty());
    assert!(te.cvs().is_empty());
    te.add_cv(CV::new("MS", "PSI-MS", "1.0", "url"));
    te.set_cvs(vec![CV::new("MS", "PSI-MS", "1.0", "url")]);
    assert_eq!(te.cvs().len(), 1);
    te.add_contact(Contact::default());
    te.set_contacts(vec![Contact::default()]);
    assert_eq!(te.contacts().len(), 1);
    te.add_publication(Publication::default());
    te.set_publications(vec![Publication::default()]);
    assert_eq!(te.publications().len(), 1);
    te.add_instrument(Instrument::default());
    te.set_instruments(vec![Instrument::default()]);
    assert_eq!(te.instruments().len(), 1);
    te.add_software(Software::default());
    te.set_software(vec![Software::default()]);
    assert_eq!(te.software().len(), 1);
    te.set_proteins(vec![Protein::default()]);
    assert_eq!(te.proteins().len(), 1);
    te.add_protein(Protein::default());
    assert_eq!(te.proteins().len(), 2);
    te.set_compounds(vec![Compound::default()]);
    assert_eq!(te.compounds().len(), 1);
    te.add_compound(Compound::default());
    assert_eq!(te.compounds().len(), 2);
    te.set_peptides(vec![Peptide::default()]);
    assert_eq!(te.peptides().len(), 1);
    te.add_peptide(Peptide::default());
    assert_eq!(te.peptides().len(), 2);
    te.set_transitions(vec![ReactionMonitoringTransition::new()]);
    assert_eq!(te.transitions().len(), 1);
    te.add_transition(ReactionMonitoringTransition::new());
    assert_eq!(te.transitions().len(), 2);
    te.add_include_target(IncludeExcludeTarget::new());
    te.add_exclude_target(IncludeExcludeTarget::new());
    te.set_include_targets(vec![IncludeExcludeTarget::new(); 2]);
    te.set_exclude_targets(vec![IncludeExcludeTarget::new(); 3]);
    assert_eq!(te.include_targets().len(), 2);
    assert_eq!(te.exclude_targets().len(), 3);
    te.add_source_file(SourceFile::default());
    te.set_source_files(vec![SourceFile::default(); 2]);
    assert_eq!(te.source_files().len(), 2);
    te.add_target_cv_term(charge_cv(2)).unwrap();
    te.set_target_meta_value("key", MetaValue::from("value"));
    assert!(te.target_cv_terms().contains("MS:1000041"));
    assert!(te.target_cv_terms().metadata.contains_key("key"));

    let copy = te.clone();
    assert_eq!(copy, te);
    let mut other = te.clone();
    other.add_transition(ReactionMonitoringTransition::new());
    assert_ne!(other, te);
    let mut replaced = TargetedExperiment::new();
    assert_ne!(replaced, te);
    replaced.set_target_cv_terms(te.target_cv_terms().clone());
    assert_eq!(replaced.target_cv_terms(), te.target_cv_terms());
}

/// START_SECTION((bool hasProtein/hasPeptide/hasCompound)) and
/// START_SECTION((const Protein& getProteinByRef(const std::string& ref) const)),
/// including its cache-invalidation regression block.
#[test]
fn te_lookup_by_reference_and_invalidation() {
    let mut t = TargetedExperiment::new();
    let mut p = Protein::default();
    p.id = "myProtein".into();
    t.add_protein(p);
    assert!(t.has_protein("myProtein"));
    assert_eq!(t.protein_by_ref("myProtein").unwrap().id, "myProtein");
    assert!(matches!(
        t.protein_by_ref("does_not_exist"),
        Err(Error::MissingInformation(_))
    ));

    let mut pep = Peptide::default();
    pep.id = "myPep".into();
    t.add_peptide(pep);
    assert!(t.has_peptide("myPep"));
    assert_eq!(t.peptide_by_ref("myPep").unwrap().id, "myPep");
    assert!(t.peptide_by_ref("does_not_exist").is_err());

    let mut comp = Compound::default();
    comp.id = "myComp".into();
    t.add_compound(comp);
    assert!(t.has_compound("myComp"));
    assert_eq!(t.compound_by_ref("myComp").unwrap().id, "myComp");
    assert!(t.compound_by_ref("does_not_exist").is_err());

    let mut te2 = TargetedExperiment::new();
    let compound = |id: &str| {
        let mut c = Compound::default();
        c.id = id.into();
        c
    };
    te2.add_compound(compound("c1"));
    assert_eq!(te2.compound_by_ref("c1").unwrap().id, "c1");
    te2.add_compound(compound("c2"));
    assert_eq!(te2.compound_by_ref("c2").unwrap().id, "c2");
    te2.set_compounds(vec![compound("c3")]);
    assert_eq!(te2.compound_by_ref("c3").unwrap().id, "c3");
    assert!(te2.compound_by_ref("c1").is_err());

    let protein = |id: &str| Protein {
        id: id.into(),
        ..Protein::default()
    };
    te2.add_protein(protein("p1"));
    assert_eq!(te2.protein_by_ref("p1").unwrap().id, "p1");
    te2.set_proteins(vec![protein("p2")]);
    assert_eq!(te2.protein_by_ref("p2").unwrap().id, "p2");
    assert!(te2.protein_by_ref("p1").is_err());

    // Native: a repeated id resolves to the last occurrence, as the source's
    // map-building loop leaves it.
    let mut dup = TargetedExperiment::new();
    dup.add_protein(Protein {
        id: "x".into(),
        sequence: "FIRST".into(),
        ..Protein::default()
    });
    dup.add_protein(Protein {
        id: "x".into(),
        sequence: "LAST".into(),
        ..Protein::default()
    });
    assert_eq!(dup.protein_by_ref("x").unwrap().sequence, "LAST");
}

/// START_SECTION( bool TargetedExperiment::containsInvalidReferences() )
#[test]
fn te_contains_invalid_references() {
    {
        let mut tr = TargetedExperiment::new();
        assert!(!tr.contains_invalid_references());
        let mut peptide = Peptide::default();
        peptide.id = "peptide1".into();
        peptide.protein_refs.push("protein1".into());
        tr.add_peptide(peptide);
        assert!(tr.contains_invalid_references());
        tr.add_protein(Protein {
            id: "protein1".into(),
            ..Protein::default()
        });
        assert!(!tr.contains_invalid_references());
        let mut t = ReactionMonitoringTransition::new();
        t.set_native_id("tr1");
        t.set_peptide_ref("peptide1");
        tr.add_transition(t);
        assert!(!tr.contains_invalid_references());
        let mut t2 = ReactionMonitoringTransition::new();
        t2.set_native_id("tr2");
        t2.set_compound_ref("compound1");
        tr.add_transition(t2);
        assert!(tr.contains_invalid_references());
        assert_eq!(
            tr.first_invalid_reference().unwrap(),
            "Compound compound1 is not present in the provided data structure."
        );
        let mut comp = Compound::default();
        comp.id = "compound1".into();
        tr.add_compound(comp);
        assert!(!tr.contains_invalid_references());
    }
    {
        let mut tr = TargetedExperiment::new();
        assert!(!tr.contains_invalid_references());
        let mut peptide = Peptide::default();
        peptide.id = "peptide1".into();
        tr.add_peptide(peptide.clone());
        assert!(!tr.contains_invalid_references());
        tr.add_peptide(peptide);
        assert!(tr.contains_invalid_references());
    }
    {
        let mut tr = TargetedExperiment::new();
        assert!(!tr.contains_invalid_references());
        let mut t = ReactionMonitoringTransition::new();
        t.set_native_id("tr1");
        tr.add_transition(t.clone());
        assert!(!tr.contains_invalid_references());
        tr.add_transition(t);
        assert!(tr.contains_invalid_references());
    }
    {
        let mut tr = TargetedExperiment::new();
        assert!(!tr.contains_invalid_references());
        let protein = Protein {
            id: "protein1".into(),
            ..Protein::default()
        };
        tr.add_protein(protein.clone());
        assert!(!tr.contains_invalid_references());
        tr.add_protein(protein);
        assert!(tr.contains_invalid_references());
        assert_eq!(
            tr.first_invalid_reference().unwrap(),
            "Found duplicate protein id (must be unique): protein1"
        );
    }
}

/// Native, independently derived from `TargetedExperiment.cpp`: merge
/// appends every list and the target CV terms but not their metadata; the
/// owned merge agrees with the borrowed one; clear(false) keeps everything
/// but the transitions.
#[test]
fn te_merge_and_clear() {
    let mut a = TargetedExperiment::new();
    a.add_protein(Protein {
        id: "p1".into(),
        ..Protein::default()
    });
    a.add_transition(ReactionMonitoringTransition::new());
    let mut b = TargetedExperiment::new();
    b.add_protein(Protein {
        id: "p2".into(),
        ..Protein::default()
    });
    b.add_target_cv_term(charge_cv(2)).unwrap();
    b.set_target_meta_value("only_in_b", MetaValue::from(1));
    let mut t = ReactionMonitoringTransition::new();
    t.set_native_id("b1");
    b.add_transition(t);

    // Build a's protein index first: the merge must invalidate it.
    assert!(a.protein_by_ref("p2").is_err());
    let merged = a.merged(&b).unwrap();
    assert_eq!(merged.proteins().len(), 2);
    assert_eq!(merged.transitions().len(), 2);
    assert!(merged.target_cv_terms().contains("MS:1000041"));
    assert!(!merged.target_cv_terms().metadata.contains_key("only_in_b"));
    assert_eq!(merged.protein_by_ref("p2").unwrap().id, "p2");

    let mut owned = a.clone();
    owned.merge_owned(b.clone()).unwrap();
    assert_eq!(owned, merged);

    let mut cleared = merged.clone();
    cleared.clear(false);
    assert!(cleared.transitions().is_empty());
    assert_eq!(cleared.proteins().len(), 2);
    assert!(cleared.target_cv_terms().contains("MS:1000041"));
    cleared.clear(true);
    assert_eq!(cleared, TargetedExperiment::new());
    assert!(!cleared.has_protein("p1"));
}

/// Native: stable sorting by product m/z and by name.
#[test]
fn te_sorting() {
    let mut te = TargetedExperiment::new();
    for (name, mz) in [("c", 300.0), ("a", 100.0), ("b", 300.0), ("d", 200.0)] {
        let mut t = ReactionMonitoringTransition::new();
        t.set_name(name);
        t.set_product_mz(mz);
        te.add_transition(t);
    }
    te.sort_transitions_by_product_mz();
    let names: Vec<_> = te
        .transitions()
        .iter()
        .map(|t| t.name().to_owned())
        .collect();
    assert_eq!(names, ["a", "d", "c", "b"]);
    te.sort_transitions_by_name();
    let names: Vec<_> = te
        .transitions()
        .iter()
        .map(|t| t.name().to_owned())
        .collect();
    assert_eq!(names, ["a", "b", "c", "d"]);
}

/// Native, derived from `TargetedExperiment.cpp` `operator<<`: 2 target and
/// 1 decoy transition give 66.667% and 33.333% with `toStr(double, false)`.
#[test]
fn te_summary_statistics() {
    let mut te = TargetedExperiment::new();
    for (i, kind) in [
        DecoyTransitionType::Target,
        DecoyTransitionType::Decoy,
        DecoyTransitionType::Target,
    ]
    .into_iter()
    .enumerate()
    {
        let mut t = ReactionMonitoringTransition::new();
        t.set_native_id(format!("t{i}"));
        t.set_decoy_transition_type(kind);
        te.add_transition(t);
    }
    let s = te.summary();
    assert_eq!(s.transition_count, 3);
    assert_eq!(s.decoy_counts[&DecoyTransitionType::Target], 2);
    assert_eq!(s.decoy_counts[&DecoyTransitionType::Decoy], 1);
    assert!(!s.decoy_counts.contains_key(&DecoyTransitionType::Unknown));
    assert!(!s.contains_invalid_references);
    assert_eq!(
        s.to_string(),
        "# Proteins: 0\n# Peptides: 0\n# Compounds: 0\n# Transitions: 3\nTransition Type: 66.667% (target), 33.333% (decoy)\nAll internal references valid: yes\n"
    );
    let empty = TargetedExperiment::new().summary();
    assert!(empty.to_string().contains("Transition Type: \n"));
}
