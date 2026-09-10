// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Literal mass expectations: AASequence_test.cpp:348–410, OpenMS4-core 6bfc0e4.
//! Fragment/terminal cases derive from AASequence.cpp:480–620 and Residue.h.

use openms::Error;
use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationRecord, ModificationsDB, PROTON_MASS_U,
    PeptideFragmentType as F, ResidueModification, TermSpecificity,
};

const TYPES: [F; 19] = [
    F::Full,
    F::Internal,
    F::NTerminal,
    F::CTerminal,
    F::AIon,
    F::BIon,
    F::CIon,
    F::XIon,
    F::YIon,
    F::ZIon,
    F::Zp1Ion,
    F::Zp2Ion,
    F::Precursor,
    F::BIonMinusH2O,
    F::YIonMinusH2O,
    F::BIonMinusNH3,
    F::YIonMinusNH3,
    F::NonIdentified,
    F::Unannotated,
];
fn sequence(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}
fn formula(text: &str) -> EmpiricalFormula {
    EmpiricalFormula::parse(text).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.14} vs {b:.14}");
}

#[test]
fn source_literal_average_mono_and_terminal_masses() {
    let peptide = sequence("DFPIANGER");
    near(
        peptide.average_mass_for(F::Full, 0).unwrap(),
        1018.08088,
        0.01,
    );
    near(peptide.average_mass_for(F::YIon, 1).unwrap(), 1019.09, 0.01);
    near(peptide.mono_mass_for(F::Full, 0).unwrap(), 1017.48796, 1e-5);
    // Source's short decimal is checked with its original relative tolerance.
    near(peptide.mono_mass_for(F::YIon, 1).unwrap(), 1018.4952, 0.001);
    near(
        sequence("(NIC)DFPIANGER")
            .mono_mass_for(F::Full, 0)
            .unwrap(),
        1122.51,
        0.01,
    );
    near(
        sequence("(dNIC)DFPIANGER")
            .mono_mass_for(F::Full, 0)
            .unwrap(),
        1017.48796 + 109.048119,
        1e-5,
    );
    near(
        sequence("(UniMod:51)CPEPTIDE")
            .mono_mass_for(F::Full, 0)
            .unwrap(),
        902.3691545801998 + 788.725777,
        1e-6,
    );
}

#[test]
fn alanine_source_ion_algebra_and_all_type_charge_combinations() {
    for (kind, neutral) in [
        (F::Internal, "C3H5NO"),
        (F::Full, "C3H7NO2"),
        (F::AIon, "C2H5N"),
        (F::BIon, "C3H5NO"),
        (F::YIon, "C3H7NO2"),
        (F::ZIon, "C3H4O2"),
    ] {
        near(
            sequence("A").mono_mass_for(kind, 1).unwrap(),
            formula(neutral).mono_mass() + PROTON_MASS_U,
            1e-12,
        );
    }
    let peptide = sequence("ACDEF");
    for kind in TYPES {
        for charge in [i32::MIN, -2, 0, 1, 3, i32::MAX] {
            let expected = peptide.formula_for(kind, charge).unwrap();
            near(
                peptide.mono_mass_for(kind, charge).unwrap(),
                expected.mono_mass(),
                1e-6,
            );
            assert_eq!(
                peptide.average_mass_for(kind, charge).unwrap().to_bits(),
                expected.average_mass().to_bits()
            );
        }
    }
}

#[test]
fn terminal_mass_tags_follow_type_selection_but_never_invent_average_mass() {
    let base = sequence("AG");
    let n = sequence("n[+12.3456789]AG");
    let c = sequence("AGc[+23.4567891]");
    for kind in TYPES {
        let keeps_n = matches!(kind, F::Full | F::NTerminal | F::AIon | F::BIon | F::CIon);
        let keeps_c = matches!(kind, F::Full | F::CTerminal | F::XIon | F::YIon | F::ZIon);
        for (modified, keeps, delta) in [(&n, keeps_n, 12.3456789), (&c, keeps_c, 23.4567891)] {
            near(
                modified.mono_mass_for(kind, 2).unwrap(),
                base.mono_mass_for(kind, 2).unwrap() + if keeps { delta } else { 0. },
                1e-12,
            );
            if keeps {
                assert!(matches!(
                    modified.average_mass_for(kind, 2),
                    Err(Error::Unsupported(_))
                ));
            } else {
                assert_eq!(
                    modified.average_mass_for(kind, 2).unwrap(),
                    base.average_mass_for(kind, 2).unwrap()
                );
            }
        }
    }
}

#[test]
fn mass_only_and_absolute_custom_residues_propagate_without_formula() {
    let db = ModificationsDB::from_records(vec![
        ResidueModification::from_record(ModificationRecord {
            name: "Delta".into(),
            origin: Some('M'),
            diff_mono_mass: 12.5,
            ..Default::default()
        })
        .unwrap(),
        ResidueModification::from_record(ModificationRecord {
            name: "Absolute".into(),
            origin: Some('X'),
            diff_mono_mass: 1.,
            mono_mass: 300.,
            ..Default::default()
        })
        .unwrap(),
        ResidueModification::from_record(ModificationRecord {
            name: "DeclaredTerminal".into(),
            term_specificity: TermSpecificity::NTerm,
            diff_mono_mass: 12.5,
            diff_formula: formula("O"),
            ..Default::default()
        })
        .unwrap(),
    ])
    .unwrap();
    let changed = AASequence::parse_with_registry("AM(Delta)K", &db).unwrap();
    let absolute = AASequence::parse_with_registry("AX(Absolute)K", &db).unwrap();
    for kind in TYPES {
        near(
            changed.mono_mass_for(kind, 0).unwrap(),
            sequence("AMK").mono_mass_for(kind, 0).unwrap() + 12.5,
            1e-10,
        );
        near(
            absolute.mono_mass_for(kind, 0).unwrap(),
            sequence("AK").mono_mass_for(kind, 0).unwrap() + 300. - formula("H2O").mono_mass(),
            1e-10,
        );
        assert!(matches!(
            changed.average_mass_for(kind, 0),
            Err(Error::Unsupported(_))
        ));
        assert!(matches!(
            absolute.average_mass_for(kind, 0),
            Err(Error::Unsupported(_))
        ));
    }
    let terminal = AASequence::parse_with_registry(".(DeclaredTerminal)AG", &db).unwrap();
    near(
        terminal.mono_mass_for(F::BIon, 1).unwrap(),
        sequence("AG").mono_mass_for(F::BIon, 1).unwrap() + 12.5,
        1e-12,
    );
    near(
        terminal.average_mass_for(F::BIon, 1).unwrap(),
        sequence("AG").average_mass_for(F::BIon, 1).unwrap() + formula("O").average_mass(),
        1e-12,
    );
    assert_eq!(
        terminal.mono_mass_for(F::YIon, 0).unwrap(),
        sequence("AG").mono_mass_for(F::YIon, 0).unwrap()
    );
    let tagged = sequence("AX[999]K");
    near(
        tagged.prefix(2).unwrap().mono_mass_for(F::BIon, 1).unwrap(),
        sequence("A").mono_mass_for(F::BIon, 1).unwrap() + 999.,
        1e-12,
    );
    near(
        tagged.suffix(2).unwrap().mono_mass_for(F::YIon, 1).unwrap(),
        sequence("K").mono_mass_for(F::YIon, 1).unwrap() + 999.,
        1e-12,
    );
}

#[test]
fn empty_and_unknown_residues_and_signed_finite_outputs() {
    for kind in TYPES {
        for charge in [i32::MIN, -1, 0, 1, i32::MAX] {
            assert_eq!(
                AASequence::default().mono_mass_for(kind, charge).unwrap(),
                0.
            );
            assert_eq!(
                AASequence::default()
                    .average_mass_for(kind, charge)
                    .unwrap(),
                0.
            );
        }
        for text in ["ABG", "AZG", "AXG"] {
            assert!(matches!(
                sequence(text).mono_mass_for(kind, 0),
                Err(Error::Unsupported(_))
            ));
        }
    }
    let balanced = sequence("n[+1000.123456789]Ac[-1000.123456789]");
    assert!(balanced.mono_mass_for(F::YIon, 0).unwrap() < 0.);
    assert!(sequence("A").mono_mass_for(F::Full, -100).unwrap() < 0.);
    assert!(sequence("A").average_mass_for(F::Full, -100).unwrap() < 0.);
    // Existing m/z API remains a checked physical mass-to-charge query.
    assert!(sequence("A").mz(-100).is_err());
}

#[test]
fn large_terminal_source_order_is_explicit_and_queries_do_not_mutate() {
    let large = format!("1{}", "0".repeat(200));
    let peptide = sequence(&format!("n[+{large}]AGc[-{large}]"));
    let original = peptide.clone();
    // This new API reproduces literal source scalar order: a small initial
    // proton contribution is rounded away before opposite terminals cancel.
    assert_eq!(
        peptide.mono_mass_for(F::Full, 1).unwrap(),
        peptide.mono_mass_for(F::Full, 0).unwrap()
    );
    assert_eq!(
        peptide.mono_mass().unwrap(),
        sequence("AG").mono_mass().unwrap()
    );
    assert!(matches!(
        peptide.average_mass_for(F::Full, 1),
        Err(Error::Unsupported(_))
    ));
    assert_eq!(peptide, original);
}
