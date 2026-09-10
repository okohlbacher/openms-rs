// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Independent source literals and atomic-constant arithmetic. Provenance and
// distinctions from the source's loose mass checks: data/adduct_provenance.json.

use openms::chemistry::EmpiricalFormula;
use openms::chemistry::adduct_info::AdductInfo;

// Literal values from pinned Constants.h and ElementDB.cpp, independent of the
// Rust constants/formula mass APIs under test.
const PROTON: f64 = 1.007_276_466_771;
const ELECTRON: f64 = 1.0 / 1_822.888_502_047_7;
const HYDROGEN: f64 = 1.007_825_031_9;
const SODIUM: f64 = 22.989_769_280_9;
const OXYGEN: f64 = 15.994_915;

fn formula(text: &str) -> EmpiricalFormula {
    EmpiricalFormula::parse(text).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17} != {expected:.17}"
    );
}

#[test]
fn source_class_test_literals_cover_forward_inverse_shift_and_compatibility() {
    for (name, charge, multiplier, expected_mz) in [
        ("M+H;1+", 1, 1, 100.0 + PROTON),
        ("2M+H;1+", 1, 2, 200.0 + PROTON),
        ("M-H;1-", -1, 1, 100.0 - PROTON),
    ] {
        let adduct = AdductInfo::parse(name).unwrap();
        assert_eq!(adduct.name(), name);
        assert_eq!(adduct.charge(), charge);
        assert_eq!(adduct.mol_multiplier(), multiplier);
        assert_eq!(adduct.empirical_formula().charge(), 0);
        // Upstream AdductInfo_test.cpp explicitly sets absolute tolerance 1e-3.
        close(adduct.mz(100.0).unwrap(), expected_mz, 1e-3);
    }
    for name in ["M+H;1+", "M+Na;1+", "2M+H;1+", "M-H;1-"] {
        let adduct = AdductInfo::parse(name).unwrap();
        close(
            adduct.neutral_mass(adduct.mz(523.25).unwrap()).unwrap(),
            523.25,
            1e-3,
        );
    }
    let proton = AdductInfo::new("my_proton", formula("H"), 1, 1).unwrap();
    assert_eq!(proton.name(), "my_proton");
    close(proton.mass_shift(false).unwrap(), 0.0, 1e-3);
    let sodium = AdductInfo::parse("M+Na;1+").unwrap();
    close(sodium.mass_shift(false).unwrap(), SODIUM - HYDROGEN, 1e-3);
    assert!(sodium.mass_shift(true).unwrap() > 0.0);
    let loss = AdductInfo::parse("M-H;1-").unwrap();
    assert!(loss.is_compatible(&formula("C6H12O6")));
    assert!(!loss.is_compatible(&formula("O2")));
    assert!(sodium.is_compatible(&formula("O2")));
}

#[test]
fn exact_atom_electron_arithmetic_preserves_small_shift_and_average_conventions() {
    let proton = AdductInfo::parse("M+H;1+").unwrap();
    close(
        proton.mz(100.0).unwrap(),
        (100.0 + HYDROGEN) - ELECTRON,
        1e-13,
    );
    let shift = HYDROGEN - (PROTON + ELECTRON);
    assert_ne!(
        shift, 0.0,
        "source atomic and particle constants are independently tabulated"
    );
    close(proton.mass_shift(false).unwrap(), shift, 1e-15);
    let average_hydrogen = HYDROGEN * 0.999_885 + 2.014_101_78 * 0.000_115;
    close(
        proton.mass_shift(true).unwrap(),
        average_hydrogen - (PROTON + ELECTRON),
        1e-15,
    );
    let sodium = AdductInfo::parse("M+Na;1+").unwrap();
    // Sodium has one isotope with abundance 1: the source comment suggesting
    // different sodium average/mono results is not an assertion or a fact.
    assert_eq!(
        sodium.mass_shift(false).unwrap(),
        sodium.mass_shift(true).unwrap()
    );

    let dimer = AdductInfo::parse("2M+Na-H;2-").unwrap();
    let adduct_mass = SODIUM - HYDROGEN;
    let expected = ((123.5 * 2.0 + adduct_mass) + 2.0 * ELECTRON) / 2.0;
    close(dimer.mz(123.5).unwrap(), expected, 1e-12);
    close(dimer.neutral_mass(expected).unwrap(), 123.5, 1e-12);
    close(
        dimer.mass_shift(false).unwrap(),
        adduct_mass + 2.0 * (PROTON + ELECTRON),
        1e-12,
    );
    let water_loss = AdductInfo::parse("M-H2O;1+").unwrap();
    close(
        water_loss.mz(100.0).unwrap(),
        (100.0 - 2.0 * HYDROGEN - OXYGEN) - ELECTRON,
        1e-12,
    );
    let isotope = AdductInfo::parse("M+(13)C+D;1+").unwrap();
    close(
        isotope.mz(100.0).unwrap(),
        (100.0 + 13.003_355 + 2.014_101_78) - ELECTRON,
        1e-12,
    );
}

#[test]
fn source_parser_keeps_numeric_noops_signed_charge_and_exact_name_identity() {
    for (text, expected, charge, multiplier) in [
        ("M;1+", "", 1, 1),
        ("M+2;1+", "", 1, 1),
        ("M+0H+Na;1+", "Na", 1, 1),
        ("M+H-1;1+", "H", 1, 1),
        ("M+H2-3;1+", "H2", 1, 1),
        ("M+H-H;1+", "", 1, 1),
        ("M;001+", "", 1, 1),
        ("M;-2+", "", 2, 1),
        ("M;+2-", "", -2, 1),
        ("M;+-2+", "", 2, 1),
        ("M;+-2-", "", -2, 1),
        ("002M+03H;2+", "H3", 2, 2),
        ("M+(13)C+D;1+", "(13)C(2)H", 1, 1),
    ] {
        let adduct: AdductInfo = text.parse().unwrap();
        assert_eq!(adduct.name(), text);
        assert_eq!(adduct.empirical_formula(), &formula(expected), "{text}");
        assert_eq!(adduct.charge(), charge, "{text}");
        assert_eq!(adduct.mol_multiplier(), multiplier, "{text}");
    }
    let compact = AdductInfo::parse("2M+H;2+").unwrap();
    assert_eq!(AdductInfo::parse("\t2 M + H ; \n2 +\r").unwrap(), compact);
    let alternate = AdductInfo::parse("02M+H;02+").unwrap();
    assert_eq!(alternate.empirical_formula(), compact.empirical_formula());
    assert_eq!(alternate.charge(), compact.charge());
    assert_ne!(
        alternate, compact,
        "equality includes the retained textual name"
    );
}

#[test]
fn compatibility_compares_signed_counts_without_multimer_or_charge_adjustment() {
    let dimer_loss = AdductInfo::parse("2M-2H;1-").unwrap();
    assert!(
        !dimer_loss.is_compatible(&formula("H")),
        "source tests the monomer counts directly"
    );
    assert!(dimer_loss.is_compatible(&formula("H2")));
    assert!(dimer_loss.is_compatible(&formula("H2").with_charge(999)));
    let gain = AdductInfo::parse("M+2H;1+").unwrap();
    assert!(gain.is_compatible(&formula("H-2")));
    assert!(
        !gain.is_compatible(&formula("H-3")),
        "plus-only compatibility assumes nonnegative candidate counts"
    );
    let isotope_loss = AdductInfo::parse("M-(13)C;1-").unwrap();
    assert!(!isotope_loss.is_compatible(&formula("C100")));
    assert!(isotope_loss.is_compatible(&formula("(13)C")));
    let limit_loss = AdductInfo::new("large loss", formula("H-2147483648"), -1, 1).unwrap();
    assert!(!limit_loss.is_compatible(&formula("H2147483647")));
    let bare = AdductInfo::parse("M;1+").unwrap();
    assert!(bare.is_compatible(&formula("H-2147483648")));
}

#[test]
fn source_finite_signed_mass_arithmetic_is_kept_and_undefined_numeric_cases_error() {
    let bare = AdductInfo::parse("M;2+").unwrap();
    for mass in [-10.0, 0.0, 10.0] {
        close(bare.mz(mass).unwrap(), (mass - 2.0 * ELECTRON) / 2.0, 1e-14);
        close(
            bare.neutral_mass(mass).unwrap(),
            mass * 2.0 + 2.0 * ELECTRON,
            1e-14,
        );
    }
    let before = bare.clone();
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(bare.mz(value).is_err());
        assert!(bare.neutral_mass(value).is_err());
    }
    assert!(bare.neutral_mass(f64::MAX).is_err());
    assert_eq!(bare, before);
    assert!(AdductInfo::new("zero", formula("H"), 0, 1).is_err());
    assert!(AdductInfo::new("zero", formula("H"), 1, 0).is_err());
    assert!(AdductInfo::new("charged formula", formula("H+"), 1, 1).is_err());
    assert!(AdductInfo::new("undefined abs", formula(""), i32::MIN, 1).is_err());
    assert_eq!(
        AdductInfo::new("large multimer", formula(""), 1, u32::MAX)
            .unwrap()
            .mol_multiplier(),
        u32::MAX
    );
    for invalid in [
        "",
        "M",
        "M;1",
        "M;+",
        "M;0+",
        "0M;1+",
        "-2M;1+",
        "H;1+",
        "M;;1+",
        "M+;1+",
        "M++H;1+",
        "M+-H;1+",
        "M+%H;1+",
        "M+Notanelement;1+",
        "M;++2+",
        "M;1.0+",
        "M;\u{b}2+",
        "M;\u{c}2+",
        "2147483648M;1+",
        "M+2147483648;1+",
        "M;2147483648-",
        "M;-2147483648+",
        "M;-2147483648-",
        "M+2147483647H+H;1+",
    ] {
        assert!(AdductInfo::parse(invalid).is_err(), "{invalid:?}");
    }
}
