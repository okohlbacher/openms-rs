// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `MRMAssay_test.cpp` at core pin `bc9cc12`, plus native checks.
//! Literals are transcribed from that class test (evidence tier 3, source
//! review) unless a test says otherwise. The class test's TraML-based sections
//! (`reannotateTransitions`, `restrictTransitions`, `detectingTransitions`,
//! `filterMinMaxTransitionsCompound`, `uisTransitions`) need `TraMLFile`,
//! which is not ported; they are replaced here by native checks whose
//! expectations are derived by hand from the source algorithm (tier 4).
//!
//! The field-by-field set-up after `default()` mirrors the C++ tests, and
//! `Peptide` fields set through `Deref` cannot move into a struct literal, so
//! the lint that asks for one is allowed here.

#![allow(clippy::field_reassign_with_default)]

use openms::analysis::openswath::mrm_assay::{MRMAssay, UisOptions};
use openms::analysis::openswath::mrm_ion_series::{IonSeriesOptions, MRMIonSeries};
use openms::analysis::openswath::reaction_monitoring_transition::{
    DecoyTransitionType, ReactionMonitoringTransition,
};
use openms::analysis::openswath::targeted_experiment::TargetedExperiment;
use openms::analysis::openswath::targeted_experiment_helper::{
    Compound, Peptide, Protein, RetentionTime,
};
use openms::chemistry::{AASequence, PeptideFragmentType};
use openms::metadata::MetaValue;

fn swathes() -> Vec<(f64, f64)> {
    let mut swathes = vec![(400.0, 425.0)];
    let mut lower = 424.0;
    while lower < 1174.5 {
        swathes.push((lower, lower + 26.0));
        lower += 25.0;
    }
    swathes
}

// MRMAssay_test.cpp:88-112
#[test]
fn matching_peptidoforms() {
    let mut ions = vec![
        (100.00_f64, "PEPTIDEK".to_owned()),
        (100.01, "PEPTIDEK".to_owned()),
        (100.10, "PEPT(UniMod:21)IDEK".to_owned()),
        (100.12, "PEPTIDEK".to_owned()),
        (100.11, "PEPTIDEK".to_owned()),
    ];
    ions.sort_by(|a, b| a.0.total_cmp(&b.0));
    assert!(MRMAssay::matching_peptidoforms(100.06, &ions, 0.03).is_empty());
    assert_eq!(
        MRMAssay::matching_peptidoforms(100.06, &ions, 0.06),
        vec!["PEPT(UniMod:21)IDEK".to_owned(), "PEPTIDEK".to_owned()]
    );
}

// MRMAssay_test.cpp:116-202
#[test]
fn swath_index_and_is_in_swath() {
    let swathes = swathes();
    assert_eq!(swathes.len(), 32);
    assert_eq!(swathes[31], (1174.0, 1200.0));
    assert_eq!(MRMAssay::swath_index(&swathes, 427.229959), 1);
    assert_eq!(MRMAssay::swath_index(&swathes, 449.0), 2);
    assert_eq!(MRMAssay::swath_index(&swathes, 449.229959), 2);
    assert_eq!(MRMAssay::swath_index(&swathes, 685.8547721), 11);
    assert_eq!(MRMAssay::swath_index(&swathes, 1685.8547721), -1);
    assert_eq!(MRMAssay::swath_index(&swathes, -41.1), -1);
    assert!(!MRMAssay::is_in_swath(&swathes, 685.8547721, 427.229959));
    assert!(MRMAssay::is_in_swath(&swathes, 685.8547721, 689.0));
    assert_eq!(MRMAssay::swath_index(&[], 500.0), -1);
}

// MRMAssay_test.cpp:206-220
#[test]
fn random_sequence() {
    assert_eq!(MRMAssay::random_sequence_with_seed(10, 42), "CHLNHHQQNE");
}

// MRMAssay_test.cpp:222-296
#[test]
fn n_choose_k_combinations() {
    let n: Vec<usize> = (1..=16).collect();
    let combinations = MRMAssay::n_choose_k_combinations(&n, 5).unwrap();
    assert_eq!(combinations.len(), 4368);
    assert_eq!(combinations[0].len(), 5);

    let n: Vec<usize> = (1..=5).collect();
    let res = MRMAssay::n_choose_k_combinations(&n, 3).unwrap();
    let expected: Vec<Vec<usize>> = vec![
        vec![1, 2, 3],
        vec![1, 2, 4],
        vec![1, 2, 5],
        vec![1, 3, 4],
        vec![1, 3, 5],
        vec![1, 4, 5],
        vec![2, 3, 4],
        vec![2, 3, 5],
        vec![2, 4, 5],
        vec![3, 4, 5],
    ];
    assert_eq!(res, expected);
    // Native: the source's bitmask truncation for k > n, and the ceiling.
    assert_eq!(
        MRMAssay::n_choose_k_combinations(&[1, 2], 5).unwrap(),
        vec![vec![1, 2]]
    );
    assert_eq!(
        MRMAssay::n_choose_k_combinations(&[1, 2], 0).unwrap(),
        vec![Vec::<usize>::new()]
    );
    let many: Vec<usize> = (0..60).collect();
    assert!(MRMAssay::n_choose_k_combinations(&many, 10).is_err());
}

// MRMAssay_test.cpp:298-344
#[test]
fn add_modifications_sequences() {
    let sequences = vec![AASequence::parse("PEPTDIEK").unwrap()];
    let oxidation = MRMAssay::n_choose_k_combinations(&[1, 3, 5, 8, 9], 1).unwrap();
    let sequences =
        MRMAssay::add_modifications_sequences(&sequences, &oxidation, "Oxidation").unwrap();
    let phospho = MRMAssay::n_choose_k_combinations(&[4, 5, 8], 1).unwrap();
    let sequences = MRMAssay::add_modifications_sequences(&sequences, &phospho, "Phospho").unwrap();
    let strings: Vec<String> = sequences.iter().map(ToString::to_string).collect();
    assert_eq!(
        strings,
        vec![
            "P(Oxidation)EPT(Phospho)DIEK",
            "P(Oxidation)EPTD(Phospho)IEK",
            "P(Oxidation)EPTDIEK(Phospho)",
            "PEP(Oxidation)T(Phospho)DIEK",
            "PEP(Oxidation)TD(Phospho)IEK",
            "PEP(Oxidation)TDIEK(Phospho)",
            "PEPT(Phospho)D(Oxidation)IEK",
            "PEPTD(Oxidation)IEK(Phospho)",
            "PEPT(Phospho)DIEK(Oxidation)",
            "PEPTD(Phospho)IEK(Oxidation)",
        ]
    );
    assert!(!strings.iter().any(|s| s == "PEPT(Phospho)DIEK.(Oxidation)"));
}

fn parse_all(texts: &[&str]) -> Vec<AASequence> {
    texts
        .iter()
        .map(|t| AASequence::parse(t).unwrap())
        .collect()
}

// MRMAssay_test.cpp:346-364
#[test]
fn theoretical_peptidoforms() {
    let sequences = MRMAssay::theoretical_peptidoforms(
        &AASequence::parse(".(Acetyl)PEPT(Phospho)DIEK").unwrap(),
    )
    .unwrap();
    assert_eq!(
        sequences,
        parse_all(&[
            ".(Acetyl)PE(Phospho)PTDIEK",
            ".(Acetyl)PEPT(Phospho)DIEK",
            ".(Acetyl)PEPTD(Phospho)IEK",
            ".(Acetyl)PEPTDIE(Phospho)K",
            ".(Acetyl)PEPTDIEK(Phospho)",
            "PE(Phospho)PT(Acetyl)DIEK",
            "PEPT(Acetyl)D(Phospho)IEK",
            "PEPT(Acetyl)DIE(Phospho)K",
            "PEPT(Acetyl)DIEK(Phospho)",
            "PE(Phospho)PTDIEK(Acetyl)",
            "PEPT(Phospho)DIEK(Acetyl)",
            "PEPTD(Phospho)IEK(Acetyl)",
            "PEPTDIE(Phospho)K(Acetyl)",
        ])
    );
}

// MRMAssay_test.cpp:366-384
#[test]
fn theoretical_peptidoforms_decoy() {
    let sequences = MRMAssay::theoretical_peptidoforms_decoy(
        &AASequence::parse(".(Acetyl)PEPT(Phospho)DIEK").unwrap(),
        &AASequence::parse("PESTDIEK").unwrap(),
    )
    .unwrap();
    assert_eq!(
        sequences,
        parse_all(&[
            ".(Acetyl)PE(Phospho)STDIEK",
            ".(Acetyl)PEST(Phospho)DIEK",
            ".(Acetyl)PESTD(Phospho)IEK",
            ".(Acetyl)PESTDIE(Phospho)K",
            ".(Acetyl)PESTDIEK(Phospho)",
            "PE(Phospho)ST(Acetyl)DIEK",
            "PEST(Acetyl)D(Phospho)IEK",
            "PEST(Acetyl)DIE(Phospho)K",
            "PEST(Acetyl)DIEK(Phospho)",
            "PE(Phospho)STDIEK(Acetyl)",
            "PEST(Phospho)DIEK(Acetyl)",
            "PESTD(Phospho)IEK(Acetyl)",
            "PESTDIE(Phospho)K(Acetyl)",
        ])
    );
}

fn peptide(id: &str, sequence: &str, charge: i32, protein: &str) -> Peptide {
    let mut p = Peptide::default();
    p.id = id.into();
    p.sequence = sequence.into();
    p.set_charge_state(charge);
    p.protein_refs.push(protein.into());
    let mut rt = RetentionTime::new();
    rt.set_rt(12.5);
    p.rts.push(rt);
    p
}

fn protein(id: &str) -> Protein {
    let mut p = Protein::default();
    p.id = id.into();
    p
}

fn transition(
    id: &str,
    peptide_ref: &str,
    precursor: f64,
    product: f64,
    intensity: f64,
) -> ReactionMonitoringTransition {
    let mut t = ReactionMonitoringTransition::new();
    t.set_native_id(id);
    t.set_name(id);
    t.set_peptide_ref(peptide_ref);
    t.set_precursor_mz(precursor);
    t.set_product_mz(product);
    t.set_library_intensity(intensity);
    t
}

// Native (tier 4): the PEPTIDEK y5^3 and y1^1 values of the MRMIonSeries
// class test; a transition 0.1 Th off every ion is dropped, an off precursor
// likewise; kept transitions snap to theory and gain their CV terms.
#[test]
fn reannotate_transitions() {
    let mut exp = TargetedExperiment::new();
    exp.add_protein(protein("P1"));
    exp.add_peptide(peptide("pep", "PEPTIDEK", 3, "P1"));
    let precursor = AASequence::parse("PEPTIDEK").unwrap().mz(3).unwrap();
    exp.add_transition(transition("t1", "pep", precursor + 0.01, 202.44, 10.0));
    exp.add_transition(transition("t2", "pep", precursor, 147.1, 10.0));
    exp.add_transition(transition("t3", "pep", precursor, 300.0, 10.0));
    exp.add_transition(transition("t4", "pep", precursor + 1.0, 202.44, 10.0));
    let options = IonSeriesOptions::new(["b", "y"], [1, 2, 3]);
    MRMAssay::reannotate_transitions(&mut exp, 0.05, 0.05, &options).unwrap();
    let names: Vec<&str> = exp.transitions().iter().map(|t| t.native_id()).collect();
    assert_eq!(names, vec!["t1", "t2"]);
    let t1 = &exp.transitions()[0];
    assert_eq!(t1.product_mz(), 202.4429);
    assert_eq!(
        t1.cv_terms.metadata.get("annotation"),
        Some(&MetaValue::from("y5^3"))
    );
    assert_eq!(t1.product().charge_state(), Some(3));
    let t2 = &exp.transitions()[1];
    assert_eq!(t2.product_mz(), 147.1128);
    assert_eq!(
        t2.product().interpretation_list()[0].iontype,
        PeptideFragmentType::YIon
    );
    assert_eq!(exp.peptides().len(), 1);

    // A transition of a missing peptide is an error and leaves exp unchanged.
    exp.add_transition(transition("t5", "missing", precursor, 202.44, 1.0));
    let before = exp.clone();
    assert!(MRMAssay::reannotate_transitions(&mut exp, 0.05, 0.05, &options).is_err());
    assert_eq!(exp, before);
}

// Native (tier 4): the three restriction rules in order.
#[test]
fn restrict_transitions() {
    let mut exp = TargetedExperiment::new();
    exp.add_peptide(peptide("pep", "PEPTIDEK", 2, "P1"));
    let mut unannotated = transition("unannotated", "pep", 500.0, 700.0, 1.0);
    MRMIonSeries::annotate_transition_cv(&mut unannotated, "q3^1").unwrap();
    exp.add_transition(unannotated);
    exp.add_transition(transition("in_swath", "pep", 500.0, 510.0, 1.0));
    exp.add_transition(transition("low", "pep", 500.0, 350.0, 1.0));
    exp.add_transition(transition("high", "pep", 500.0, 2100.0, 1.0));
    exp.add_transition(transition("keep", "pep", 500.0, 700.0, 1.0));
    exp.add_transition(transition("no_swath", "pep", 3000.0, 700.0, 1.0));
    let options = IonSeriesOptions::new(["y"], [1]);
    let mut without = exp.clone();
    MRMAssay::restrict_transitions(&mut without, 400.0, 2000.0, &[], &options).unwrap();
    let names: Vec<&str> = without
        .transitions()
        .iter()
        .map(|t| t.native_id())
        .collect();
    assert_eq!(names, vec!["in_swath", "keep", "no_swath"]);
    MRMAssay::restrict_transitions(&mut exp, 400.0, 2000.0, &swathes(), &options).unwrap();
    let names: Vec<&str> = exp.transitions().iter().map(|t| t.native_id()).collect();
    assert_eq!(names, vec!["keep"]);
}

// Native (tier 4): top-N by intensity with the tie cap, decoy exclusion,
// minimum count, and peptide/protein pruning.
#[test]
fn detecting_transitions() {
    let mut exp = TargetedExperiment::new();
    exp.add_protein(protein("P1"));
    exp.add_protein(protein("P2"));
    exp.add_peptide(peptide("a", "PEPTIDEK", 2, "P1"));
    exp.add_peptide(peptide("b", "PEPTIDER", 2, "P2"));
    for (id, intensity) in [
        ("a1", 5.0),
        ("a2", 9.0),
        ("a3", 9.0),
        ("a4", 1.0),
        ("a5", 9.0),
    ] {
        exp.add_transition(transition(id, "a", 500.0, 600.0, intensity));
    }
    let mut decoy = transition("a6", "a", 500.0, 600.0, 20.0);
    decoy.set_decoy_transition_type(DecoyTransitionType::Decoy);
    exp.add_transition(decoy);
    exp.add_transition(transition("b1", "b", 500.0, 600.0, 3.0));
    MRMAssay::detecting_transitions(&mut exp, 2, 3);
    let names: Vec<&str> = exp.transitions().iter().map(|t| t.native_id()).collect();
    // Top three intensities {20, 9, 9} -> set {20, 9}; the decoy is skipped,
    // the three 9s fill the limit, and b has too few transitions.
    assert_eq!(names, vec!["a2", "a3", "a5"]);
    assert!(
        exp.transitions()
            .iter()
            .all(|t| t.is_detecting_transition())
    );
    assert_eq!(exp.peptides().len(), 1);
    assert_eq!(exp.proteins().len(), 1);
    assert_eq!(exp.proteins()[0].id, "P1");
}

fn compound(id: &str) -> Compound {
    let mut c = Compound::default();
    c.id = id.into();
    c
}

fn compound_transition(
    id: &str,
    compound_ref: &str,
    intensity: f64,
    decoy: bool,
) -> ReactionMonitoringTransition {
    let mut t = transition(id, "", 100.0, 50.0, intensity);
    t.set_compound_ref(compound_ref);
    if decoy {
        t.set_decoy_transition_type(DecoyTransitionType::Decoy);
    }
    t
}

// Native (tier 4): compound top-N, decoy groups below the minimum, and the
// unreferenced-decoy filter.
#[test]
fn compound_filters() {
    let mut exp = TargetedExperiment::new();
    for id in [
        "1_Cmp_[M+H]+_1",
        "1_Cmp_decoy_[M+H]+_1",
        "2_Other_decoy_[M+H]+_2",
    ] {
        exp.add_compound(compound(id));
    }
    exp.add_transition(compound_transition("t1", "1_Cmp_[M+H]+_1", 1.0, false));
    exp.add_transition(compound_transition("t2", "1_Cmp_[M+H]+_1", 3.0, false));
    exp.add_transition(compound_transition("t3", "1_Cmp_[M+H]+_1", 2.0, false));
    exp.add_transition(compound_transition("d1", "1_Cmp_decoy_[M+H]+_1", 1.0, true));
    exp.add_transition(compound_transition(
        "d2",
        "2_Other_decoy_[M+H]+_2",
        1.0,
        true,
    ));

    let mut filtered = exp.clone();
    MRMAssay::filter_min_max_transitions_compound(&mut filtered, 3, 2);
    let names: Vec<&str> = filtered
        .transitions()
        .iter()
        .map(|t| t.native_id())
        .collect();
    assert_eq!(names, vec!["t2", "t3", "d1", "d2"]);
    assert_eq!(filtered.compounds().len(), 3);

    MRMAssay::filter_unreferenced_decoys_compound(&mut exp).unwrap();
    let ids: Vec<&str> = exp.compounds().iter().map(|c| c.id.as_str()).collect();
    assert_eq!(ids, vec!["1_Cmp_[M+H]+_1", "1_Cmp_decoy_[M+H]+_1"]);
    let names: Vec<&str> = exp.transitions().iter().map(|t| t.native_id()).collect();
    assert_eq!(names, vec!["t1", "t2", "t3", "d1"]);

    let mut bad = TargetedExperiment::new();
    bad.add_compound(compound("mydecoy"));
    let before = bad.clone();
    assert!(MRMAssay::filter_unreferenced_decoys_compound(&mut bad).is_err());
    assert_eq!(bad, before);
}

// Native (tier 4): UIS structure for one phosphopeptide. PEPS(Phospho)TIDEK
// has two phospho sites (S and T) and so two peptidoforms; every appended
// transition is identifying, non-detecting, non-quantifying, named after
// its index, peptidoforms, m/z, RT and annotation, and the decoys are
// reproducible for a fixed seed.
#[test]
fn uis_transitions() {
    let mut exp = TargetedExperiment::new();
    let mut p = peptide("pep", "PEPSTIDEK", 2, "P1");
    let mut m =
        openms::analysis::openswath::targeted_experiment_helper::PeptideModification::default();
    m.location = 3;
    m.unimod_id = 21;
    m.mono_mass_delta = 79.966331;
    p.mods.push(m);
    exp.add_peptide(p);
    let existing = transition("existing", "pep", 500.0, 600.0, 1.0);
    exp.add_transition(existing.clone());

    let mut options = UisOptions::new(IonSeriesOptions::new(["y", "b"], [1]), 0.05);
    options.shuffle_seed = Some(42);
    let mut first = exp.clone();
    MRMAssay::uis_transitions(&mut first, &options).unwrap();
    let mut second = exp.clone();
    MRMAssay::uis_transitions(&mut second, &options).unwrap();
    assert_eq!(first, second);

    let transitions = first.transitions();
    assert_eq!(transitions[0], existing);
    let appended = &transitions[1..];
    assert!(!appended.is_empty());
    let targets: Vec<_> = appended
        .iter()
        .filter(|t| t.decoy_transition_type() != DecoyTransitionType::Decoy)
        .collect();
    let decoys: Vec<_> = appended
        .iter()
        .filter(|t| t.decoy_transition_type() == DecoyTransitionType::Decoy)
        .collect();
    assert!(!targets.is_empty() && !decoys.is_empty());
    for t in appended {
        assert!(t.is_identifying_transition());
        assert!(!t.is_detecting_transition());
        assert!(!t.is_quantifying_transition());
        assert_eq!(
            t.cv_terms.metadata.get("insilico_transition"),
            Some(&MetaValue::from("true"))
        );
        assert_eq!(t.name(), t.native_id());
        assert!(
            t.name()
                .ends_with(t.cv_terms.metadata["annotation"].as_str().unwrap())
        );
        assert!(t.name().contains("_12.5_"));
    }
    // Every peptidoform carries exactly one phospho group on the unchanged
    // backbone, and the listed peptidoforms are sorted and distinct.
    for t in &targets {
        assert!(t.name().contains("_UIS_{"));
        let forms: Vec<&str> = t.cv_terms.metadata["Peptidoforms"]
            .as_str()
            .unwrap()
            .split('|')
            .collect();
        assert!(forms.windows(2).all(|w| w[0] < w[1]));
        for form in forms {
            assert_eq!(form.matches("(Phospho)").count(), 1, "{form}");
            assert_eq!(AASequence::parse(form).unwrap().as_str(), "PEPSTIDEK");
        }
    }
    // The first target transition has index 0 and the rounded precursor m/z.
    let precursor = AASequence::parse("PEPS(Phospho)TIDEK")
        .unwrap()
        .mz(2)
        .unwrap();
    // Source Math::roundDecimal(x, -4): floor(0.5 + x / 10^-4) * 10^-4.
    let scale = 10.0_f64.powf(-4.0);
    let rounded = (0.5 + precursor / scale).floor() * scale;
    assert_eq!(targets[0].precursor_mz(), rounded);
    assert!(targets[0].name().starts_with("0_UIS_{"));
    assert!(decoys.iter().all(|t| t.name().contains("_UISDECOY_{")));
    assert!(decoys.iter().all(|t| t.precursor_mz() == rounded));

    // Disabling decoys appends only the targets.
    let mut no_decoys = exp.clone();
    options.disable_decoy_transitions = true;
    MRMAssay::uis_transitions(&mut no_decoys, &options).unwrap();
    assert_eq!(no_decoys.transitions().len(), 1 + targets.len());
}
