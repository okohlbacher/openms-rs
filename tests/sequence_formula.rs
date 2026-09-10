// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Source literals: AASequence_test.cpp:341–346 at OpenMS4-core 6bfc0e4.
//! Other expectations derive directly from AASequence.cpp:383–478 and the
//! elemental corrections in Residue.h:63–143; no C++ execution is claimed.

use openms::Error;
use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationRecord, ModificationsDB, PeptideFragmentType as F,
    ResidueModification, TermSpecificity,
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
fn formula(text: &str) -> EmpiricalFormula {
    EmpiricalFormula::parse(text).unwrap()
}
fn sequence(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}

#[test]
fn source_literal_full_charged_and_b_formulas() {
    let peptide = sequence("ACDEF");
    assert_eq!(
        peptide.formula_for(F::Full, 0).unwrap(),
        formula("O10SH33N5C24")
    );
    assert_eq!(
        peptide.formula_for(F::Full, 1).unwrap(),
        formula("O10SH33N5C24+")
    );
    assert_eq!(
        peptide.formula_for(F::BIon, 0).unwrap(),
        formula("O9SH31N5C24")
    );
    assert_eq!(
        peptide.formula_for(F::Full, 0).unwrap(),
        peptide.formula().unwrap()
    );
}

#[test]
fn every_formula_type_has_source_correction_and_terminal_selection() {
    let peptide = sequence("(Acetyl)AG.(Amidated)");
    // AG internal = C5H8N2O2; N acetyl = C2H2O; C amide = HNO-1.
    let expected = [
        "C7H13N3O3",
        "C5H8N2O2",
        "C7H11N2O3",
        "C5H10N3O2",
        "C6H10N2O2",
        "C7H10N2O3",
        "C7H13N3O3",
        "C6H9N3O3",
        "C5H11N3O2",
        "C5H8N2O2",
        "C5H8N2O2",
        "C5H8N2O2",
        "C5H8N2O2",
        "C5H8N2O2",
        "C5H8N2O2",
        "C5H8N2O2",
        "C5H8N2O2",
        "C5H8N2O2",
        "C5H8N2O2",
    ];
    for (kind, expected) in TYPES.into_iter().zip(expected) {
        for charge in [-2, 0, 3] {
            assert_eq!(
                peptide.formula_for(kind, charge).unwrap(),
                formula(expected).with_charge(charge),
                "{kind:?}, {charge}"
            );
        }
    }
    assert_eq!(
        peptide.formula_for(F::Full, 0).unwrap(),
        peptide.formula().unwrap()
    );
}

#[test]
fn charge_is_metadata_including_extremes_and_does_not_add_hydrogen() {
    let peptide = sequence("PEPTIDE");
    for kind in TYPES {
        let neutral = peptide.formula_for(kind, 0).unwrap();
        for charge in [i32::MIN, -3, 1, i32::MAX] {
            let charged = peptide.formula_for(kind, charge).unwrap();
            assert_eq!(charged.charge(), charge);
            assert_eq!(charged.clone().with_charge(0), neutral);
            assert_eq!(charged.count("H").unwrap(), neutral.count("H").unwrap());
        }
    }
}

#[test]
fn fragment_kind_uses_all_supplied_residues_and_slicing_selects_the_length() {
    let whole = sequence("(Acetyl)AC(Carbamidomethyl)M(Oxidation)K.(Amidated)");
    let prefix = whole.prefix(3).unwrap();
    let suffix = whole.suffix(2).unwrap();
    assert_eq!(
        prefix.formula_for(F::BIon, 0).unwrap(),
        sequence("(Acetyl)AC(Carbamidomethyl)M(Oxidation)")
            .formula()
            .unwrap()
            .checked_sub(&formula("H2O"))
            .unwrap()
    );
    assert_eq!(
        suffix.formula_for(F::YIon, 0).unwrap(),
        sequence("M(Oxidation)K.(Amidated)").formula().unwrap()
    );
    assert_ne!(
        whole.formula_for(F::BIon, 0).unwrap(),
        prefix.formula_for(F::BIon, 0).unwrap()
    );
}

#[test]
fn discarded_unknown_terminal_composition_does_not_block_fragments() {
    let base = sequence("AG");
    let unknown_n = sequence("n[+12.3456789]AG");
    let unknown_c = sequence("AGc[+12.3456789]");
    for kind in TYPES {
        let keeps_n = matches!(kind, F::Full | F::NTerminal | F::AIon | F::BIon | F::CIon);
        let keeps_c = matches!(kind, F::Full | F::CTerminal | F::XIon | F::YIon | F::ZIon);
        for (peptide, retained) in [(&unknown_n, keeps_n), (&unknown_c, keeps_c)] {
            if retained {
                assert!(
                    matches!(peptide.formula_for(kind, 2), Err(Error::Unsupported(_))),
                    "{kind:?}"
                );
            } else {
                assert_eq!(
                    peptide.formula_for(kind, 2).unwrap(),
                    base.formula_for(kind, 2).unwrap()
                );
            }
        }
    }
    for input in ["ABG", "AZG", "AXG", "AX[999]G", "AG[+12.3456789]"] {
        for kind in TYPES {
            assert!(matches!(
                sequence(input).formula_for(kind, 0),
                Err(Error::Unsupported(_))
            ));
        }
    }
}

#[test]
fn custom_absolute_residue_formula_and_unknown_terminal_follow_existing_chemistry() {
    let db = ModificationsDB::from_records(vec![
        ResidueModification::from_record(ModificationRecord {
            name: "Restore".into(),
            origin: Some('X'),
            diff_mono_mass: 1.0,
            absolute_formula: Some(formula("C3H7NO2")),
            ..Default::default()
        })
        .unwrap(),
        ResidueModification::from_record(ModificationRecord {
            name: "MassOnly".into(),
            origin: Some('M'),
            diff_mono_mass: 12.5,
            ..Default::default()
        })
        .unwrap(),
        ResidueModification::from_record(ModificationRecord {
            name: "Terminal".into(),
            term_specificity: TermSpecificity::NTerm,
            diff_mono_mass: 12.5,
            absolute_formula: Some(formula("C3H7NO2")),
            ..Default::default()
        })
        .unwrap(),
    ])
    .unwrap();
    let restored = AASequence::parse_with_registry("GX(Restore)", &db).unwrap();
    for kind in TYPES {
        assert_eq!(
            restored.formula_for(kind, -1).unwrap(),
            sequence("GA").formula_for(kind, -1).unwrap()
        );
        let missing = AASequence::parse_with_registry("GM(MassOnly)", &db).unwrap();
        assert!(matches!(
            missing.formula_for(kind, 0),
            Err(Error::Unsupported(_))
        ));
    }
    let terminal = AASequence::parse_with_registry(".(Terminal)GA", &db).unwrap();
    assert!(matches!(
        terminal.formula_for(F::Full, 0),
        Err(Error::Unsupported(_))
    ));
    assert_eq!(
        terminal.formula_for(F::YIon, 0).unwrap(),
        sequence("GA").formula().unwrap()
    );
}

#[test]
fn source_order_checked_charge_overflow_and_signed_fragment_atom_counts() {
    let db = ModificationsDB::from_records(vec![
        ResidueModification::from_record(ModificationRecord {
            name: "Plus".into(),
            term_specificity: TermSpecificity::NTerm,
            diff_formula: EmpiricalFormula::default().with_charge(1),
            ..Default::default()
        })
        .unwrap(),
        ResidueModification::from_record(ModificationRecord {
            name: "Minus".into(),
            term_specificity: TermSpecificity::CTerm,
            diff_formula: EmpiricalFormula::default().with_charge(-1),
            ..Default::default()
        })
        .unwrap(),
        ResidueModification::from_record(ModificationRecord {
            name: "Carbon".into(),
            origin: Some('X'),
            diff_mono_mass: 1.,
            absolute_formula: Some(formula("CH2O")),
            ..Default::default()
        })
        .unwrap(),
    ])
    .unwrap();
    let peptide = AASequence::parse_with_registry(".(Plus)AG.(Minus)", &db).unwrap();
    let original = peptide.clone();
    assert_eq!(
        peptide.formula_for(F::Full, 0).unwrap(),
        sequence("AG").formula().unwrap()
    );
    // N charge is added first, so MAX + 1 fails before the C charge cancels it.
    assert!(matches!(
        peptide.formula_for(F::Full, i32::MAX),
        Err(Error::InvalidValue(_))
    ));
    assert!(peptide.formula_for(F::YIon, i32::MIN).is_err());
    assert_eq!(
        peptide.formula_for(F::Internal, i32::MAX).unwrap().charge(),
        i32::MAX
    );
    assert_eq!(peptide, original);
    let signed = AASequence::parse_with_registry("X(Carbon)", &db).unwrap();
    // Full residue CH2O leaves internal C. Source z correction is H-1N-1O;
    // a fragment formula is signed algebra, not a new validated sequence.
    assert_eq!(
        signed.formula_for(F::ZIon, -1).unwrap(),
        formula("CH-1N-1O").with_charge(-1)
    );
}

#[test]
fn empty_sequence_ignores_every_formula_type_and_charge() {
    for kind in TYPES {
        for charge in [i32::MIN, -1, 0, 1, i32::MAX] {
            assert_eq!(
                AASequence::default().formula_for(kind, charge).unwrap(),
                EmpiricalFormula::default()
            );
        }
    }
}
