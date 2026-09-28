// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The `CHEMISTRY/EmpiricalFormula.h` members added by lane A5a.
//!
//! Expected values are the assertions of `EmpiricalFormula_test.cpp` at core
//! `bc9cc12` (tier 2; the section names are cited per test). The text forms and
//! the theoretical isotope count have no exact upstream value; they are derived
//! from the `EmpiricalFormula.cpp` definitions at the same pin (lines 79-103,
//! 258-267, 407-440) and the tabulated natural isotopes (tier 4).

use openms::chemistry::{EmpiricalFormula, element};
use openms::concept::constants::NEUTRON_MASS_U;
use std::collections::HashSet;

fn ef(text: &str) -> EmpiricalFormula {
    EmpiricalFormula::parse(text).unwrap()
}

/// `TEST_REAL_SIMILAR` with the class-test defaults, as `tests/binned_spectrum.rs`
/// states them: absolute difference within 1e-5, or ratio within 1 + 1e-5.
fn assert_real_similar(actual: f64, expected: f64) {
    let absdiff = (actual - expected).abs();
    let ratio = if actual.abs() > expected.abs() {
        actual / expected
    } else {
        expected / actual
    };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{actual} is not similar to {expected}"
    );
}

/// `EmpiricalFormula(SignedSize number, const Element* element, SignedSize charge=0)`.
#[test]
fn from_element_is_the_parsed_formula() {
    let c4 = EmpiricalFormula::from_element(4, element("C").unwrap(), 0);
    assert_eq!(c4, ef("C4"));
    assert_eq!(c4.charge(), 0);
    let charged = EmpiricalFormula::from_element(6, element("C").unwrap(), 1);
    assert_eq!(charged, ef("C6+"));
}

/// `getNumberOf(const Element*)`, `hasElement(const Element*)`.
#[test]
fn number_of_and_has_element() {
    let c4 = ef("C4");
    assert_eq!(c4.number_of(element("C").unwrap()), 4);
    assert!(c4.has_element(element("C").unwrap()));
    assert!(!c4.has_element(element("H").unwrap()));
    // A labelled isotope is a different ElementDB element in the source.
    let labelled = ef("(13)C2");
    assert!(!labelled.has_element(element("C").unwrap()));
    assert_eq!(labelled.count("(13)C").unwrap(), 2);
}

/// `setCharge`, `getCharge`, `isCharged`.
#[test]
fn set_charge_and_is_charged() {
    let mut c4 = ef("C4");
    c4.set_charge(1);
    assert_eq!(c4.charge(), 1);
    assert!(c4.is_charged());
    c4.set_charge(0);
    assert!(!c4.is_charged());
    assert_eq!(ef("C2+3").charge(), 3);
}

/// `operator<`.
#[test]
fn ordering_matches_the_source_cases() {
    assert!(ef("C5H2") < ef("C6H2"));
    assert!(ef("C5H2") < ef("C5H3"));
    assert!(ef("C5") < ef("C5H2"));
    assert!(ef("C5H2") >= ef("C4H2"));
    assert!(ef("C5") >= ef("C5"));
}

/// `addChargeAdduct(Int count, const std::string& adduct)`.
#[test]
fn add_charge_adduct_matches_the_source_cases() {
    let h = element("H").unwrap();
    let mut glucose = ef("C6H12O6");
    glucose.add_charge_adduct(2, "H").unwrap();
    assert_eq!(glucose.charge(), 0);
    assert_eq!(glucose.number_of(h), 14);
    assert_eq!(glucose, ef("C6H14O6"));
    assert_real_similar(glucose.mono_mass(), ef("C6H14O6").mono_mass());

    let mut na = ef("C6H12O6");
    na.add_charge_adduct(1, "Na").unwrap();
    assert_eq!(na.number_of(element("Na").unwrap()), 1);
    assert_eq!(na.charge(), 0);

    let mut deprotonated = ef("C6H12O6");
    deprotonated.add_charge_adduct(-1, "H").unwrap();
    assert_eq!(deprotonated.number_of(h), 11);

    let mut charged = ef("C6H12O6");
    charged.set_charge(2);
    let count = charged.charge();
    charged.add_charge_adduct(count, "H").unwrap();
    assert_eq!(charged.charge(), 0);
    assert_eq!(charged.number_of(h), 14);
}

/// A charge written in the adduct is ignored; an unparsable adduct leaves the
/// formula unchanged (source: `EmpiricalFormula add(adduct); add.setCharge(0)`).
#[test]
fn add_charge_adduct_ignores_adduct_charge_and_rejects_bad_text() {
    let mut formula = ef("C2+");
    formula.add_charge_adduct(1, "Na+").unwrap();
    assert_eq!(formula, ef("C2Na"));
    let before = formula.clone();
    assert!(formula.add_charge_adduct(1, "Xx").is_err());
    assert_eq!(formula, before);
}

/// `getLightestIsotopeWeight`.
#[test]
fn lightest_isotope_mass_matches_the_source_case() {
    let formula = ef("C636H1007N167O178S5PtH-2");
    assert_real_similar(
        formula.lightest_isotope_mass(),
        formula.mono_mass() - NEUTRON_MASS_U * 3.0,
    );
}

/// `toString`, `toMap`, `operator<<`.
#[test]
fn text_forms_and_symbol_map() {
    let formula = ef("C2H5");
    assert_eq!(formula.to_source_string(), "C2H5");
    let map = formula.to_map();
    assert_eq!(map["C"], 2);
    assert_eq!(map["H"], 5);
    assert!(formula.source_stream_text().contains("H5"));
    assert!(formula.source_stream_text().contains("C2"));

    // toString writes every count and never the charge; operator<< writes
    // counts above one only, then the charge.
    assert_eq!(ef("CH4+").to_source_string(), "C1H4");
    assert_eq!(ef("CH4+").source_stream_text(), "CH4+");
    assert_eq!(ef("C2+3").source_stream_text(), "C2+3");
    assert_eq!(ef("C2-").source_stream_text(), "C2-");
    assert_eq!(ef("H-2O").to_source_string(), "H-2O1");
    assert_eq!(ef("H-2O").source_stream_text(), "HO");
    // `os << "-" << charge_` with a negative charge (EmpiricalFormula.cpp:431).
    assert_eq!(ef("C2").with_charge(-2).source_stream_text(), "C2--2");
    // Symbols sort as strings: "(13)C" before "C".
    assert_eq!(ef("C2(13)C").to_source_string(), "(13)C1C2");
}

/// `operator==` with the `D` and `T` aliases.
#[test]
fn deuterium_and_tritium_aliases_are_equal() {
    assert_eq!(ef("C5(2)H5"), ef("C5D5"));
    assert_eq!(ef("C5(3)H5"), ef("C5T5"));
    assert_eq!(ef("C5D5").to_map()["(2)H"], 5);
}

/// `begin()` / `end()`.
#[test]
fn entries_yield_every_count() {
    let formula = ef("C6H12O6");
    let expected = [("C", 6), ("H", 12), ("O", 6)];
    let entries: Vec<_> = formula.entries().collect();
    assert_eq!(entries.len(), 3);
    for entry in entries {
        let (_, count) = expected
            .iter()
            .find(|(symbol, _)| *symbol == entry.symbol())
            .unwrap();
        assert_eq!(entry.count, *count);
        assert_eq!(entry.isotope, None);
    }
    let labelled: Vec<_> = ef("(13)C2").entries().collect();
    assert_eq!(labelled[0].symbol(), "(13)C");
    assert_eq!(labelled[0].isotope, Some(13));
}

/// `static hydrogen(int n_atoms = 1)` and `static water(int n_molecules = 1)`.
#[test]
fn hydrogen_and_water_factories() {
    assert_eq!(EmpiricalFormula::hydrogen(1), ef("H"));
    assert_eq!(EmpiricalFormula::water(1).unwrap(), ef("H2O"));
    assert_eq!(EmpiricalFormula::water(3).unwrap(), ef("H6O3"));
    assert!(EmpiricalFormula::water(i32::MAX).is_err());
}

/// `[EXTRA] std::hash<EmpiricalFormula>`.
#[test]
fn hash_matches_equality() {
    let glucose = ef("C6H12O6");
    assert_eq!(glucose.source_hash64(), ef("C6H12O6").source_hash64());
    assert_eq!(glucose.source_hash64(), ef("H12C6O6").source_hash64());
    assert_ne!(glucose.source_hash64(), ef("C6H12O7").source_hash64());
    assert_eq!(
        ef("C6H12O6+").source_hash64(),
        ef("C6H12O6+").source_hash64()
    );
    assert_ne!(
        ef("C6H12O6+").source_hash64(),
        ef("C6H12O6-").source_hash64()
    );
    assert_ne!(glucose.source_hash64(), ef("C6H12O6+").source_hash64());
    assert_eq!(
        EmpiricalFormula::default().source_hash64(),
        EmpiricalFormula::default().source_hash64()
    );

    let mut set = HashSet::new();
    set.insert(glucose.clone());
    set.insert(ef("C6H12O6"));
    set.insert(ef("C6H12O7"));
    assert_eq!(set.len(), 2);
    assert!(set.contains(&ef("H12C6O6")));
}

/// `calculateTheoreticalIsotopesNumber` (EmpiricalFormula.cpp:79-103): the
/// product of `C(count, k)` over entries, `k` the number of isotopes with a
/// nonzero natural abundance (C: 12C, 13C; H: 1H, 2H; O: 16O, 17O, 18O).
#[test]
fn theoretical_isotopes_number_is_the_binomial_product() {
    assert_eq!(
        EmpiricalFormula::default()
            .theoretical_isotopes_number()
            .unwrap(),
        1.0
    );
    // count == 1 takes the `count * k` branch.
    assert_eq!(ef("C").theoretical_isotopes_number().unwrap(), 2.0);
    // C(6,2) * C(12,2) * C(6,3) = 15 * 66 * 20.
    assert_eq!(
        ef("C6H12O6").theoretical_isotopes_number().unwrap(),
        19_800.0
    );
    // A labelled isotope has one isotope: count * 1.
    assert_eq!(ef("(13)C4").theoretical_isotopes_number().unwrap(), 4.0);
    // Boost raises a domain error for k > n; a negative count is refused.
    assert!(ef("O2").theoretical_isotopes_number().is_err());
    assert!(ef("C-2").theoretical_isotopes_number().is_err());
}

/// `getConditionalFragmentIsotopeDist(precursor, precursor_isotopes, solver)`.
#[test]
fn conditional_fragment_isotope_dist_matches_the_class_test() {
    use openms::chemistry::isotopes::{CoarseIsotopePatternGenerator as Generator, CoarseMassMode};

    let precursor = ef("C2");
    let fragment = ef("C");
    let solver = Generator::default();
    let check = |isotopes: &[usize], probabilities: &[f64], masses: Option<&[f64]>| {
        let iso = fragment
            .conditional_fragment_isotope_dist(&precursor, isotopes, &solver)
            .unwrap();
        assert_eq!(iso.len(), probabilities.len(), "{isotopes:?}");
        for (i, peak) in iso.peaks().iter().enumerate() {
            assert_real_similar(peak.probability, probabilities[i]);
            if let Some(masses) = masses {
                assert_real_similar(peak.mass, masses[i]);
            }
        }
    };
    check(&[0], &[1.0], Some(&[12.0]));
    check(&[1], &[0.5, 0.5], Some(&[12.0, 13.0033548378]));
    check(&[0, 1], &[0.98941, 0.01059], Some(&[12.0, 13.0033548378]));
    check(&[0, 1, 2], &[0.9893, 0.0107], None);
    // M+3 is unreachable for C2 and ignored.
    check(&[0, 1, 2, 3], &[0.9893, 0.0107], None);

    let rounded = Generator::new(None, CoarseMassMode::Nominal).unwrap();
    let iso = fragment
        .conditional_fragment_isotope_dist(&precursor, &[0, 1, 2, 3], &rounded)
        .unwrap();
    let masses: Vec<f64> = iso.peaks().iter().map(|p| p.mass).collect();
    assert_eq!(masses, [12.0, 13.0]);

    let precursor = ef("C10H10N10O10S2");
    let big = ef("C9H9N9O9S1")
        .conditional_fragment_isotope_dist(&precursor, &[1], &solver)
        .unwrap();
    let small = ef("C1H1N1O1S1")
        .conditional_fragment_isotope_dist(&precursor, &[1], &solver)
        .unwrap();
    assert!(big.peaks()[0].probability < 0.2);
    assert!(big.peaks()[1].probability > 0.8);
    assert!(small.peaks()[0].probability > 0.8);
    assert!(small.peaks()[1].probability < 0.2);
    for (b, s) in big.peaks().iter().zip(small.peaks().iter().rev()) {
        assert_real_similar(b.probability, s.probability);
    }

    assert!(
        fragment
            .conditional_fragment_isotope_dist(&precursor, &[], &solver)
            .is_err()
    );
}
