// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::fine_isotopes::{
    FineIsotopePatternGenerator as Fine, FineIsotopeStop as Stop, MAX_FINE_ATOMS, MAX_FINE_PEAKS,
};
use openms::chemistry::{EmpiricalFormula, IsotopeDistribution, element};

fn run(formula: &str, stop: Stop) -> IsotopeDistribution {
    Fine { stop }
        .run(&EmpiricalFormula::parse(formula).unwrap())
        .unwrap()
}

#[test]
fn empty_and_fixed_support_obey_inclusive_threshold_and_coverage_contracts() {
    for formula in ["", "Na", "(13)C2", "(3)H"] {
        let distribution = run(formula, Stop::AbsoluteThreshold(1.0));
        assert_eq!(distribution.len(), 1);
        assert_eq!(distribution.peaks()[0].probability, 1.0);
        assert!(run(formula, Stop::RelativeThreshold(1.0001)).is_empty());
        assert!(run(formula, Stop::UnexplainedProbability(1.0)).is_empty());
        assert_eq!(run(formula, Stop::AbsoluteThreshold(0.0)), distribution);
        assert_eq!(
            run(formula, Stop::UnexplainedProbability(0.0)),
            distribution
        );
    }
    assert_eq!(run("", Stop::RelativeThreshold(0.0)).peaks()[0].mass, 0.0);
    let tritium = run("(3)H", Stop::RelativeThreshold(0.0));
    assert_eq!(tritium.peaks()[0].mass, 3.01604927);
}

#[test]
fn natural_zero_abundances_are_skipped_but_explicit_labels_have_unit_weight() {
    let h = run("H2", Stop::AbsoluteThreshold(0.0));
    assert_eq!(h.len(), 3); // natural tritium is absent
    let mixed = run("H2(3)H", Stop::AbsoluteThreshold(0.0));
    assert_eq!(mixed.len(), 3);
    for (a, b) in h.peaks().iter().zip(mixed.peaks()) {
        assert!((b.mass - a.mass - 3.01604927).abs() < 1e-14);
        assert_eq!(a.probability, b.probability);
    }
    assert_eq!(
        run("(2)H3", Stop::RelativeThreshold(0.0)),
        run("D3", Stop::RelativeThreshold(0.0))
    );
}

#[test]
fn positive_charge_adds_natural_hydrogen_without_mz_division() {
    let neutral = EmpiricalFormula::parse("C2O2").unwrap();
    let charged = EmpiricalFormula::parse("C2O2+2").unwrap();
    let method = Fine {
        stop: Stop::RelativeThreshold(0.0),
    };
    assert_eq!(method.run(&charged).unwrap(), run("C2O2H2", method.stop));
    assert_ne!(method.run(&charged).unwrap(), method.run(&neutral).unwrap());
    // Labelled hydrogen stays fixed while the charge contributes a separate
    // natural-hydrogen isotope category.
    let labelled = EmpiricalFormula::parse("D2+").unwrap();
    assert_eq!(method.run(&labelled).unwrap(), run("D2H", method.stop));
    let negative = EmpiricalFormula::parse("C2O2-").unwrap();
    assert!(method.run(&negative).is_err());
}

fn compositions_two_atoms(isotopes: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    for first in 0..isotopes {
        for second in first..isotopes {
            let mut counts = vec![0; isotopes];
            counts[first] += 1;
            counts[second] += 1;
            result.push(counts);
        }
    }
    result
}

#[test]
fn full_small_multinomial_product_matches_independent_direct_combinatorics() {
    // All 3*3*6 configurations, using elementary two-atom probabilities rather
    // than logarithms, heap traversal or the coarse convolution implementation.
    let mut expected = vec![(0.0_f64, 1.0_f64)];
    for symbol in ["C", "H", "O"] {
        let isotopes: Vec<_> = element(symbol)
            .unwrap()
            .isotopes()
            .iter()
            .filter(|isotope| isotope.abundance > 0.0)
            .collect();
        let mut next = Vec::new();
        for counts in compositions_two_atoms(isotopes.len()) {
            let mass: f64 = counts
                .iter()
                .zip(&isotopes)
                .map(|(&count, isotope)| count as f64 * isotope.mass)
                .sum();
            let distinct = counts.iter().filter(|&&count| count > 0).count() == 2;
            let probability = counts.iter().zip(&isotopes).fold(
                if distinct { 2.0 } else { 1.0 },
                |value, (&count, isotope)| {
                    value * f64::from(isotope.abundance as f32).powi(count as i32)
                },
            );
            for &(previous_mass, previous_probability) in &expected {
                next.push((previous_mass + mass, previous_probability * probability));
            }
        }
        expected = next;
    }
    expected.sort_by(|a, b| a.0.total_cmp(&b.0));
    let actual = run("C2H2O2", Stop::RelativeThreshold(0.0));
    assert_eq!(actual.len(), 54);
    for (peak, &(mass, probability)) in actual.peaks().iter().zip(&expected) {
        assert!((peak.mass - mass).abs() < 3e-14);
        assert_eq!(peak.probability, f64::from(probability as f32));
    }
    assert!(
        actual
            .peaks()
            .windows(2)
            .all(|pair| pair[0].mass <= pair[1].mass)
    );
}

#[test]
fn abundance_rounding_and_output_rounding_do_not_renormalize_the_pattern() {
    let distribution = run("C2H2O2", Stop::RelativeThreshold(0.0));
    let sum: f64 = distribution.peaks().iter().map(|p| p.probability).sum();
    assert_ne!(sum, 1.0);
    assert!(
        distribution
            .peaks()
            .iter()
            .all(|p| p.probability == f64::from(p.probability as f32))
    );
    // Coverage one can exhaust a support whose stored weights fall below one.
    let oxygen = run("O", Stop::AbsoluteThreshold(0.0));
    assert_eq!(run("O", Stop::UnexplainedProbability(0.0)), oxygen);
}

#[test]
fn coverage_is_a_minimal_probability_prefix_and_runs_are_deterministic() {
    let all = run("C2H2O2", Stop::RelativeThreshold(0.0));
    let mut probabilities: Vec<_> = all.peaks().iter().map(|p| p.probability).collect();
    probabilities.sort_by(|a, b| b.total_cmp(a));
    for target in [0.5, 0.99, 0.9999] {
        let selected = run("C2H2O2", Stop::UnexplainedProbability(1.0 - target));
        let sum: f64 = probabilities[..selected.len()].iter().sum();
        assert!(sum >= target);
        if selected.len() > 1 {
            assert!(sum - probabilities[selected.len() - 1] < target);
        }
        assert_eq!(
            selected,
            run("C2H2O2", Stop::UnexplainedProbability(1.0 - target))
        );
    }
}

#[test]
fn invalid_settings_formulas_and_cardinality_fail_without_mutation() {
    let formula = EmpiricalFormula::parse("C6H12O6").unwrap();
    let before = formula.clone();
    for stop in [
        Stop::AbsoluteThreshold(-1.0),
        Stop::AbsoluteThreshold(f64::INFINITY),
        Stop::RelativeThreshold(f64::NAN),
        Stop::UnexplainedProbability(-0.01),
        Stop::UnexplainedProbability(1.01),
    ] {
        assert!(Fine { stop }.run(&formula).is_err());
    }
    assert_eq!(formula, before);
    for text in ["C-1", "C-1H2", &format!("C{}", MAX_FINE_ATOMS + 1)] {
        let invalid = EmpiricalFormula::parse(text).unwrap();
        assert!(Fine::default().run(&invalid).is_err());
    }
    let full = Fine {
        stop: Stop::RelativeThreshold(0.0),
    };
    assert!(
        full.run(&EmpiricalFormula::parse(&format!("C{MAX_FINE_PEAKS}")).unwrap())
            .is_err()
    );
    assert!(full.run(&EmpiricalFormula::parse("Sn20").unwrap()).is_err());
    let charged = EmpiricalFormula::parse(&format!("C+{MAX_FINE_ATOMS}")).unwrap();
    assert!(Fine::default().run(&charged).is_err());
}

#[test]
fn enormous_positive_threshold_search_stops_at_a_checked_resource_limit() {
    // Tin has ten natural isotopes: even 100 atoms have vastly more states
    // than can be materialized. A positive cutoff does not use the full-support
    // preflight, so this exercises actual bounded frontier growth.
    let formula = EmpiricalFormula::parse("Sn100").unwrap();
    let before = formula.clone();
    let error = Fine {
        stop: Stop::RelativeThreshold(1e-200),
    }
    .run(&formula)
    .unwrap_err();
    assert!(error.to_string().contains("limit"));
    assert_eq!(formula, before);
}
