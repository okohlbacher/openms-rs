// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::adduct_info::{AdductInfo, MAX_ADDUCT_TERMS, MAX_ADDUCT_TEXT_BYTES};
use openms::chemistry::{ELECTRON_MASS_U, EmpiricalFormula, PROTON_MASS_U};

fn formula(text: &str) -> EmpiricalFormula {
    EmpiricalFormula::parse(text).unwrap()
}

#[test]
fn programmatic_construction_preserves_names_and_source_four_field_equality() {
    let parsed = AdductInfo::parse("M+Na;1+").unwrap();
    let same = AdductInfo::new("M+Na;1+", formula("Na"), 1, 1).unwrap();
    assert_eq!(same, parsed);
    let other_name = AdductInfo::new("sodium", formula("Na"), 1, 1).unwrap();
    assert_ne!(same, other_name);
    assert_eq!(same.empirical_formula(), other_name.empirical_formula());
    assert_ne!(same, AdductInfo::parse("1M+Na;01+").unwrap());
    let unnamed = AdductInfo::new("", EmpiricalFormula::default(), -1, u32::MAX).unwrap();
    assert_eq!(unnamed.name(), "");
    assert_eq!(unnamed.mol_multiplier(), u32::MAX);
    assert!(AdductInfo::new("invalid charge", formula("H").with_charge(1), 1, 1).is_err());
    assert!(AdductInfo::new("zero", formula("H"), 0, 1).is_err());
    assert!(AdductInfo::new("zero multiplier", formula("H"), 1, 0).is_err());
    assert!(AdductInfo::new("minimum charge", formula("H"), i32::MIN, 1).is_err());
}

#[test]
fn mixed_terms_keep_neutral_atoms_and_distinct_charge_arithmetic() {
    let adduct: AdductInfo = "2M+2K-H;3+".parse().unwrap();
    assert_eq!(adduct.empirical_formula(), &formula("H-1K2"));
    assert_eq!(adduct.empirical_formula().charge(), 0);
    assert_eq!(adduct.charge(), 3);
    assert_eq!(adduct.mol_multiplier(), 2);
    let mass = 300.125;
    let expected = (mass * 2.0 + formula("H-1K2").mono_mass() - 3.0 * ELECTRON_MASS_U) / 3.0;
    assert_eq!(adduct.mz(mass).unwrap(), expected);
    assert!((adduct.neutral_mass(expected).unwrap() - mass).abs() < 1e-12);
    let bare = AdductInfo::parse("M;1+").unwrap();
    assert_eq!(bare.mz(mass).unwrap(), mass - ELECTRON_MASS_U);
    assert_ne!(bare.mz(mass).unwrap(), mass + PROTON_MASS_U);
}

#[test]
fn parser_keeps_exact_stripped_spelling_and_source_numeric_empty_fragments() {
    let text = " \t2 M + 02H - H\n; -02+\r";
    let adduct = AdductInfo::parse(text).unwrap();
    assert_eq!(adduct.name(), "2M+02H-H;-02+");
    assert_eq!(adduct.empirical_formula(), &formula("H"));
    assert_eq!(adduct.charge(), 2);
    for text in ["M+2;1+", "M+0;1+", "M+00C;1+", "M+H-H;1+"] {
        assert!(
            AdductInfo::parse(text)
                .unwrap()
                .empirical_formula()
                .is_empty(),
            "{text}"
        );
    }
    assert_eq!(
        AdductInfo::parse("M+H-2;1+").unwrap().empirical_formula(),
        &formula("H")
    );
    assert!(AdductInfo::parse("M+0Notanelement;1+").is_err());
}

#[test]
fn mass_conversions_are_inverses_for_signed_charges_multipliers_and_finite_signed_inputs() {
    for notation in ["M+H;1+", "3M-2H;2-", "2M+Cl;1-", "M+2Na;2+", "M;3-"] {
        let adduct = AdductInfo::parse(notation).unwrap();
        for neutral in [-123.25, 0.0, 180.0634, 1250.75] {
            let observed = adduct.mz(neutral).unwrap();
            assert!(
                (adduct.neutral_mass(observed).unwrap() - neutral).abs() < 1e-10,
                "{notation}, {neutral}"
            );
        }
    }
}

#[test]
fn nonfinite_arguments_and_intermediate_overflow_are_rejected() {
    let adduct = AdductInfo::parse("2M+H;2+").unwrap();
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(adduct.mz(value).is_err());
        assert!(adduct.neutral_mass(value).is_err());
    }
    assert!(adduct.mz(f64::MAX).is_err());
    assert!(adduct.neutral_mass(f64::MAX).is_err());
    let saved = adduct.clone();
    assert!(adduct.mz(-f64::MAX).is_err());
    assert_eq!(adduct, saved);
}

#[test]
fn compatibility_uses_signed_counts_and_minimum_count_cannot_match_any_native_candidate() {
    let large_loss = AdductInfo::new("loss", formula("H-2147483648"), -1, 1).unwrap();
    assert!(!large_loss.is_compatible(&formula("H2147483647")));
    let sodium = AdductInfo::parse("M+Na;1+").unwrap();
    assert!(sodium.is_compatible(&formula("Na-1")));
    assert!(!sodium.is_compatible(&formula("Na-2")));
    // The source check is against one formula, not the n-mer's atom supply.
    let dimer_loss = AdductInfo::parse("2M-2H;1-").unwrap();
    assert!(!dimer_loss.is_compatible(&formula("H")));
    assert!(dimer_loss.is_compatible(&formula("H2").with_charge(12)));
}

#[test]
fn malformed_and_unicode_boundary_inputs_return_errors_without_panics() {
    for text in [
        "",
        "M",
        ";1+",
        "M;",
        "M;+",
        "M;1",
        "M;0+",
        "0M;1+",
        "M+;1+",
        "M-;1+",
        "+M;1+",
        "M++H;1+",
        "M+-H;1+",
        "M-+H;1+",
        "M--H;1+",
        "M%H;1+",
        "[M+H]+",
        "M+H;1+;",
        "M;1++;",
        "M;++2+",
        "M;2147483648+",
        "M;-2147483648-",
        "2147483648M;1+",
        "M+2147483648H;1+",
        "M+H2147483647+H;1+",
        "M+H2147483647H;1+",
        "M+2H2147483647;1+",
        "M+\u{b}H;1+",
        "M+\u{c}H;1+",
        "M+H\u{a0};1+",
        "M;\u{b}1+",
        "µM+H;1+",
        "M+☃;1+",
        "M+\0;1+",
        "M+(13C;1+",
    ] {
        assert!(AdductInfo::parse(text).is_err(), "{text:?}");
    }
}

#[test]
fn text_term_and_repeated_formula_work_limits_are_checked_before_large_growth() {
    let name = "x".repeat(MAX_ADDUCT_TEXT_BYTES);
    assert!(AdductInfo::new(&name, formula("H"), 1, 1).is_ok());
    assert!(AdductInfo::new(format!("{name}x"), formula("H"), 1, 1).is_err());
    assert!(
        AdductInfo::parse(&format!("{name}x"))
            .unwrap_err()
            .to_string()
            .contains("text limit")
    );
    let max_terms = format!("M{};1+", "+0".repeat(MAX_ADDUCT_TERMS));
    assert!(AdductInfo::parse(&max_terms).is_ok());
    let excess_terms = format!("M{};1+", "+0".repeat(MAX_ADDUCT_TERMS + 1));
    assert!(
        AdductInfo::parse(&excess_terms)
            .unwrap_err()
            .to_string()
            .contains("term limit")
    );
    let repeated = format!("M{};1+", "+H".repeat(1000));
    assert!(
        AdductInfo::parse(&repeated)
            .unwrap_err()
            .to_string()
            .contains("work limit")
    );
}
