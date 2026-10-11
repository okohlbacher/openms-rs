// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `MRMDecoy_test.cpp` at core pin `bc9cc12`, plus native checks.
//! Literals are transcribed from that class test (evidence tier 3, source
//! review) unless a test says otherwise. The seeded shuffle literals exercise
//! the whole Boost `mt19937` / `uniform_int<>` chain, so they pin the random
//! stream bit for bit.
//!
//! The field-by-field set-up after `default()` mirrors the C++ tests, and
//! `Peptide` fields set through `Deref` cannot move into a struct literal, so
//! the lint that asks for one is allowed here.

#![allow(clippy::field_reassign_with_default)]

use openms::analysis::openswath::mrm_decoy::{DecoyMethod, DecoyOptions, MRM_DECOY_NAME, MRMDecoy};
use openms::analysis::openswath::mrm_ion_series::IonSeriesOptions;
use openms::analysis::openswath::reaction_monitoring_transition::{
    DecoyTransitionType, ReactionMonitoringTransition,
};
use openms::analysis::openswath::targeted_experiment::TargetedExperiment;
use openms::analysis::openswath::targeted_experiment_helper::{
    Peptide, PeptideModification, Protein, TraMLProduct,
};
use openms::metadata::CVTerm;

fn modification(location: i32, mono: f64, avg: f64) -> PeptideModification {
    let mut m = PeptideModification::default();
    m.location = location;
    m.mono_mass_delta = mono;
    m.avg_mass_delta = avg;
    m
}

fn peptide(sequence: &str) -> Peptide {
    let mut p = Peptide::default();
    p.sequence = sequence.into();
    p
}

fn shuffled(sequence: &str, threshold: f64, seed: u32, attempts: usize) -> String {
    MRMDecoy::new()
        .shuffle_peptide(&peptide(sequence), threshold, seed, attempts)
        .unwrap()
        .sequence
}

// MRMDecoy_test.cpp:86-107
#[test]
fn fixed_residues() {
    let decoy = MRMDecoy::new();
    assert_eq!(decoy.fixed_residues("TRESTPEPTIKDE"), vec![1, 5, 7, 10]);
    assert_eq!(
        decoy.fixed_and_terminal_residues("TRESTPEPTIKDE"),
        vec![0, 1, 5, 7, 10, 12]
    );
}

// MRMDecoy_test.cpp:110-248 (default max_attempts is 100 in the header)
#[test]
fn shuffle_peptide() {
    let decoy = MRMDecoy::new();
    let mut p = peptide("TESTPEPTIDE");
    p.mods.push(modification(2, 79.966331, 79.9799));
    let s = decoy.shuffle_peptide(&p, 0.7, 43, 100).unwrap();
    assert_eq!(s.sequence, "TIDEPEPSTTE");
    assert_eq!(s.mods.len(), 1);
    assert_eq!(s.mods[0].location, 7);

    assert_eq!(shuffled("TESTPEPTIDE", 1.1, 42, 100), "TESTPEPTIDE");
    assert_eq!(shuffled("TESTPEPTIDE", 0.7, 42, 100), "TTETPEPIDSE");
    assert_eq!(shuffled("TESTPEPTIDE", 0.9, 42, 100), "TTETPEPIDSE");
    assert_eq!(shuffled("TESTPEPTIDE", 0.2, 42, 10000), "TNGCADQQEAE");
    assert_eq!(shuffled("TESTPEPTIDE", 0.0, 42, 20), "TEIEPAPTQTE");
    assert_eq!(shuffled("TESTPEPTIDE", 0.0, 42, 2000), "TNDQIADNNEE");
    assert_eq!(shuffled("TESTPEPTIDEK", 0.7, 42, 20), "TTETPEPEDSIK");
    assert_eq!(shuffled("TESTPEPTIDER", 0.7, 42, 20), "TTETPEPEDSIR");

    for (location, expected_location) in [(13, 13), (14, 14)] {
        let mut p = peptide("EPAHLMSLFGGKPM");
        let mut m = PeptideModification::default();
        m.location = location;
        m.unimod_id = 35;
        p.mods.push(m);
        let s = decoy.shuffle_peptide(&p, 0.7, 42, 20).unwrap();
        assert_eq!(s.sequence, "EPSALMGGHLFKPM");
        assert_eq!(s.mods.len(), 1);
        assert_eq!(s.mods[0].location, expected_location);
    }
}

// MRMDecoy_test.cpp:252-285, regression for JIRA ABL-749
#[test]
fn shuffle_peptide_with_modifications_and_two_attempts() {
    let mut p = peptide("GPPSEDGPGVPPPSPR");
    p.mods.push(modification(3, 79.966331, 79.9799));
    p.mods.push(modification(13, 79.966331, 79.9799));
    let s = MRMDecoy::new().shuffle_peptide(&p, 0.7, 130, 100).unwrap();
    assert_eq!(s.sequence, "GPPGDSEPGSPPPVPR");
    assert_eq!(s.mods.len(), 2);
    assert_eq!(s.mods[1].location, 9);
    assert_eq!(s.mods[0].location, 5);
}

// MRMDecoy_test.cpp:289-321
#[test]
fn shuffle_peptide_with_terminal_modifications() {
    let mut p = peptide("TESTPEPTIDE");
    p.mods.push(modification(-1, 135.983029, 136.1265));
    p.mods.push(modification(11, 4.008491, 3.9995));
    let s = MRMDecoy::new().shuffle_peptide(&p, 0.7, 43, 100).unwrap();
    assert_eq!(s.sequence, "TIDEPEPSTTE");
    assert_eq!(s.mods.len(), 2);
    assert_eq!(s.mods[0].location, -1);
    assert_eq!(s.mods[1].location, 11);
}

// MRMDecoy_test.cpp:325-333
#[test]
fn shuffle_peptide_with_kpr() {
    assert_eq!(shuffled("KPRKPRPK", 0.7, 130, 17), "KNRKPRPK");
}

// MRMDecoy_test.cpp:337-346
#[test]
fn aa_sequence_identity() {
    let identity = MRMDecoy::aa_sequence_identity("TESTPEPTIDE", "EDITPEPTSET").unwrap();
    assert!((f64::from(identity) - f64::from(0.454545_f32)).abs() <= 1e-5);
    assert!(MRMDecoy::aa_sequence_identity("AB", "A").is_err());
    assert!(MRMDecoy::aa_sequence_identity("", "").unwrap().is_nan());
}

// MRMDecoy_test.cpp:350-406
#[test]
fn reverse_peptide() {
    let mut p = peptide("TESTPEPTIDE");
    p.mods.push(modification(2, 79.966331, 79.9799));
    let r = MRMDecoy::reverse_peptide(&p, false, true, "").unwrap();
    assert_eq!(r.sequence, "DITPEPTSETE");
    assert_eq!(r.mods.len(), 1);
    assert_eq!(r.mods[0].location, 7);

    p.mods.push(modification(0, 49.966331, 49.9799));
    let r = MRMDecoy::reverse_peptide(&p, true, true, "").unwrap();
    assert_eq!(r.sequence, "TDITPEPTSEE");
    assert_eq!(r.mods[0].location, 8);
    assert_eq!(r.mods[0].mono_mass_delta, 79.966331);
    assert_eq!(r.mods[1].location, 0);
    assert_eq!(r.mods[1].mono_mass_delta, 49.966331);

    let r = MRMDecoy::reverse_peptide(&p, true, true, "I").unwrap();
    assert_eq!(r.sequence, "TDTPEPTSIEE");
    assert_eq!(r.mods[0].location, 7);
    assert_eq!(r.mods[1].location, 0);
}

// MRMDecoy_test.cpp:409-467
#[test]
fn pseudoreverse_and_reverse() {
    let mut p = peptide("TESTPEPTIDE");
    p.mods.push(modification(2, 79.966331, 79.9799));
    let r = MRMDecoy::pseudoreverse_peptide(&p).unwrap();
    assert_eq!(r.sequence, "DITPEPTSETE");
    assert_eq!(r.mods[0].location, 7);
    let r = MRMDecoy::reverse_peptide_full(&p).unwrap();
    assert_eq!(r.sequence, "EDITPEPTSET");
    assert_eq!(r.mods[0].location, 8);
    assert_eq!(
        MRMDecoy::pseudoreverse_peptide(&peptide("TESTPEPTIDE"))
            .unwrap()
            .sequence,
        "DITPEPTSETE"
    );
    assert_eq!(
        MRMDecoy::reverse_peptide_full(&peptide("TESTPEPTIDE"))
            .unwrap()
            .sequence,
        "EDITPEPTSET"
    );
}

// MRMDecoy_test.cpp:890-922 (switchKRLight, which switchKR delegates to)
#[test]
fn switch_kr() {
    let mut p = peptide("TESTPEPTIDEK");
    MRMDecoy::switch_kr(&mut p);
    assert_eq!(p.sequence, "TESTPEPTIDER");
    let mut p = peptide("TESTPEPTIDER");
    MRMDecoy::switch_kr(&mut p);
    assert_eq!(p.sequence, "TESTPEPTIDEK");
    let mut p = peptide("TESTPEPTIDEX");
    MRMDecoy::switch_kr(&mut p);
    assert_eq!(p.sequence.len(), 12);
    assert_ne!(p.sequence.as_bytes()[11], b'X');
    let mut empty = peptide("");
    MRMDecoy::switch_kr(&mut empty);
    assert_eq!(empty.sequence, "");
}

// MRMDecoy_test.cpp:924-1000 (hasCNterminalModsLight_ cases on heavy peptides)
#[test]
fn has_cn_terminal_mods() {
    let with = |locations: &[i32]| {
        let mut p = peptide("TESTPEPTIDE");
        for &l in locations {
            p.mods.push(modification(l, 1.0, 1.0));
        }
        p
    };
    assert!(MRMDecoy::has_cn_terminal_mods(&with(&[-1]), false));
    assert!(MRMDecoy::has_cn_terminal_mods(&with(&[11]), false));
    assert!(MRMDecoy::has_cn_terminal_mods(&with(&[10]), true));
    assert!(!MRMDecoy::has_cn_terminal_mods(&with(&[10]), false));
    assert!(!MRMDecoy::has_cn_terminal_mods(&with(&[5]), false));
    assert!(!MRMDecoy::has_cn_terminal_mods(&with(&[]), false));
    assert!(MRMDecoy::has_cn_terminal_mods(&with(&[-1, 11]), false));
}

fn transition(
    id: &str,
    peptide_ref: &str,
    product_mz: f64,
    precursor_mz: f64,
) -> ReactionMonitoringTransition {
    let mut t = ReactionMonitoringTransition::new();
    t.set_native_id(id);
    t.set_peptide_ref(peptide_ref);
    t.set_product_mz(product_mz);
    t.set_precursor_mz(precursor_mz);
    t.set_detecting_transition(true);
    t.add_product_cv_term(CVTerm::new("MS:1001220", "frag: y ion", "MS"))
        .unwrap();
    let mut product = TraMLProduct::new();
    product.set_charge_state(1);
    t.set_product(product);
    t
}

// MRMDecoy_test.cpp:1143-1276, regression for modified-sequence duplicate
// detection on the heavy path.
#[test]
fn generate_decoys_detects_duplicates_by_modified_sequence() {
    let mut exp = TargetedExperiment::new();
    let mut protein = Protein::default();
    protein.id = "protein1".into();
    exp.add_protein(protein);
    for (id, location, unimod, mono, avg) in [
        ("peptide1", 2, 21, 79.966331, 79.9799),
        ("peptide2", 5, 35, 15.994915, 15.9994),
    ] {
        let mut p = peptide("TESTPEPTIDE");
        p.id = id.into();
        p.set_charge_state(2);
        p.protein_refs.push("protein1".into());
        let mut m = modification(location, mono, avg);
        m.unimod_id = unimod;
        p.mods.push(m);
        exp.add_peptide(p);
    }
    exp.add_transition(transition("transition1", "peptide1", 500.0, 600.0));
    exp.add_transition(transition("transition2", "peptide2", 501.0, 601.0));

    let mut options =
        DecoyOptions::new(IonSeriesOptions::new(["y", "b"], [1, 2]).with_losses(true, false));
    options.method = DecoyMethod::PseudoReverse;
    options.max_attempts = 10;
    options.product_mz_threshold = 500.0;
    let mut dec = TargetedExperiment::new();
    MRMDecoy::new()
        .generate_decoys(&exp, &mut dec, &options)
        .unwrap();

    assert_eq!(dec.peptides().len(), 2);
    let (p1, p2) = (&dec.peptides()[0], &dec.peptides()[1]);
    assert!(!p1.mods.is_empty() && !p2.mods.is_empty());
    assert_ne!(p1.mods[0].location, p2.mods[0].location);
    assert_eq!(dec.transitions().len(), 2);
    // Native checks of the documented bookkeeping.
    assert_eq!(p1.id, "DECOY_peptide1");
    assert_eq!(p1.protein_refs, vec!["DECOY_protein1".to_owned()]);
    assert_eq!(dec.proteins().len(), 1);
    assert_eq!(dec.proteins()[0].id, "DECOY_protein1");
    for t in dec.transitions() {
        assert_eq!(t.decoy_transition_type(), DecoyTransitionType::Decoy);
        assert!(t.native_id().starts_with("DECOY_transition"));
        assert!(t.peptide_ref().starts_with("DECOY_peptide"));
    }
}

// Native: the documented error paths leave the decoy experiment unchanged;
// "shift" keeps the sequence, so each decoy equals its target and is
// excluded, as in the source.
#[test]
fn generate_decoys_errors_and_shift() {
    let mut exp = TargetedExperiment::new();
    let mut p = peptide("PEPTIDEK");
    p.id = "p".into();
    p.set_charge_state(2);
    exp.add_peptide(p);
    exp.add_transition(transition("t", "p", 500.0, 600.0));
    let mut dec = TargetedExperiment::new();
    let mut marker = Protein::default();
    marker.id = "keep".into();
    dec.add_protein(marker);
    let before = dec.clone();

    let mut options = DecoyOptions::new(IonSeriesOptions::new(["y", "b"], [1]));
    options.aim_decoy_fraction = 1.5;
    assert!(
        MRMDecoy::new()
            .generate_decoys(&exp, &mut dec, &options)
            .is_err()
    );
    options.aim_decoy_fraction = f64::NAN;
    assert!(
        MRMDecoy::new()
            .generate_decoys(&exp, &mut dec, &options)
            .is_err()
    );
    assert_eq!(dec, before);

    options.aim_decoy_fraction = 1.0;
    options.method = DecoyMethod::Shift;
    MRMDecoy::new()
        .generate_decoys(&exp, &mut dec, &options)
        .unwrap();
    assert!(dec.peptides().is_empty());
    assert!(dec.transitions().is_empty());

    assert!(DecoyMethod::from_source_str("pseudo-reverse").is_ok());
    assert!(DecoyMethod::from_source_str("mirror").is_err());
}

// Native: parameters round-trip and reject values outside the valid strings.
#[test]
fn parameters() {
    let defaults = MRMDecoy::defaults().unwrap();
    let decoy = MRMDecoy::from_param(&defaults).unwrap();
    assert_eq!(decoy, MRMDecoy::new());
    assert_eq!(decoy.keep_const_pattern, "KRP");
    assert!(decoy.keep_n && decoy.keep_c);
    assert_eq!(MRM_DECOY_NAME, "MRMDecoy");
    let mut changed = MRMDecoy::new();
    changed.keep_c = false;
    changed.keep_const_pattern = "P".into();
    assert_eq!(
        MRMDecoy::from_param(&changed.to_param().unwrap()).unwrap(),
        changed
    );
}
