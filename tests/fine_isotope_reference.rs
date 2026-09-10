// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Source literals and independent chemical references, without C++ execution.
//! See data/fine_isotope_provenance.json for provenance and precision policies.

use openms::chemistry::{
    EmpiricalFormula, FineIsotopePatternGenerator as Generator, FineIsotopeStop as Stop,
    IsotopeDistribution, IsotopePeak,
};

const COUNTS: &str = include_str!("data/fine_isotope_counts.tsv");
const FRUCTOSE: &str = include_str!("data/fine_isotope_fructose.tsv");
const BROMINE: &str = include_str!("data/fine_isotope_bromine.tsv");

fn rows(table: &str) -> impl Iterator<Item = Vec<&str>> {
    table
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .skip(1)
        .map(|line| line.split('\t').collect())
}
fn number(text: &str) -> f64 {
    text.parse().unwrap()
}
fn bits64(text: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(text, 16).unwrap())
}
fn bits32(text: &str) -> f64 {
    f64::from(f32::from_bits(u32::from_str_radix(text, 16).unwrap()))
}
fn run(formula: &str, stop: Stop) -> IsotopeDistribution {
    Generator { stop }
        .run(&formula.parse::<EmpiricalFormula>().unwrap())
        .unwrap()
}
fn source_relative(actual: f64, expected: f64) {
    // Exactly the materialized comparison rule from IsoSpec_test.cpp.
    const EPSILON: f64 = 0.000_000_1;
    assert!(
        actual * (1.0 - EPSILON) <= expected && expected <= actual * (1.0 + EPSILON),
        "{actual:.17e} != {expected:.17e} within source relative tolerance"
    );
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e} (absolute tolerance {tolerance})"
    );
}
fn fructose_reference() -> Vec<IsotopePeak> {
    rows(FRUCTOSE)
        .map(|row| {
            let mass = bits64(row[2]);
            let probability = bits32(row[3]);
            assert_eq!(mass, number(row[0]));
            assert_eq!(probability, f64::from(number(row[1]) as f32));
            IsotopePeak { mass, probability }
        })
        .collect()
}
fn compare(actual: &IsotopePeak, expected: &IsotopePeak) {
    source_relative(actual.mass, expected.mass);
    source_relative(actual.probability, expected.probability);
    // The native container stores f64, but source-compatible fine output has
    // already passed through the source Peak1D probability type, f32.
    assert_eq!(actual.probability, f64::from(actual.probability as f32));
}

#[test]
fn all_44_active_source_count_cases() {
    let mut count = 0;
    for row in rows(COUNTS) {
        let value = number(row[2]);
        let stop = match row[1] {
            "relative_threshold" => Stop::RelativeThreshold(value),
            "absolute_threshold" => Stop::AbsoluteThreshold(value),
            "trimmed_coverage" => Stop::UnexplainedProbability(1.0 - value),
            _ => panic!("unknown fixture mode"),
        };
        let formula: EmpiricalFormula = row[0].parse().unwrap();
        let distribution = Generator { stop }
            .run(&formula)
            .unwrap_or_else(|error| panic!("{} {} {}: {error}", row[0], row[1], row[2]));
        assert_eq!(
            distribution.len(),
            row[3].parse::<usize>().unwrap(),
            "{} {} {} from {}",
            row[0],
            row[1],
            row[2],
            row[4]
        );
        assert!(
            distribution
                .peaks()
                .windows(2)
                .all(|p| p[0].mass <= p[1].mass)
        );
        count += 1;
    }
    assert_eq!(count, 44);
}

#[test]
fn complete_fructose_literal_table_matches_both_threshold_modes() {
    let expected = fructose_reference();
    assert_eq!(expected.len(), 14);
    for stop in [Stop::RelativeThreshold(1e-5), Stop::AbsoluteThreshold(1e-5)] {
        let actual = run("C6H12O6", stop);
        assert_eq!(actual.len(), expected.len());
        for (actual, expected) in actual.peaks().iter().zip(&expected) {
            compare(actual, expected);
        }
    }
}

#[test]
fn source_coverage_adds_three_fructose_states_beyond_the_14_row_reference() {
    let distribution = run("C6H12O6", Stop::UnexplainedProbability(0.000_01));
    assert_eq!(distribution.len(), 17);
    let mut actual = distribution.peaks().to_vec();
    actual.sort_by(|a, b| b.probability.total_cmp(&a.probability));
    let mut expected = fructose_reference();
    expected.sort_by(|a, b| b.probability.total_cmp(&a.probability));
    // Upstream supplies exact masses/probabilities for the fourteen most likely
    // states only; no invented source goldens for the remaining three.
    for (actual, expected) in actual.iter().zip(&expected) {
        compare(actual, expected);
    }
    let default = Generator::default()
        .run(&"C6H12O6".parse().unwrap())
        .unwrap();
    assert_eq!(default.len(), 3);
    let mut default = default.peaks().to_vec();
    default.sort_by(|a, b| b.probability.total_cmp(&a.probability));
    for (actual, expected) in default.iter().zip(&expected[..3]) {
        compare(actual, expected);
    }
}

#[test]
fn source_insulin_coverage_boundary_uses_stored_binary32_probabilities() {
    let target = 0.999_99;
    let distribution = run(
        "C520H817N139O147S8",
        Stop::UnexplainedProbability(1.0 - target),
    );
    assert_eq!(distribution.len(), 19_615); // Exact upstream assertion.
    let mut probabilities: Vec<_> = distribution.peaks().iter().map(|p| p.probability).collect();
    probabilities.sort_by(|a, b| b.total_cmp(a));
    assert!(probabilities.iter().all(|&p| p == f64::from(p as f32)));
    let sum: f64 = probabilities.iter().sum();
    // Independent Python multinomial prototype: raw f64 probabilities sum to
    // 0.9999899999250583 at 19615, still below target. Source f32 trimming crosses
    // it here. This expected sum is derived evidence, not captured C++ output.
    near(sum, 0.999_990_000_538_514_1, 1e-12);
    assert!(sum >= target);
    assert!(sum - probabilities.last().unwrap() < target);
}

#[test]
fn gapped_bromine_matches_source_literals_and_closed_form_probabilities() {
    let distribution = run("CBr2", Stop::RelativeThreshold(1e-3));
    let expected: Vec<_> = rows(BROMINE).collect();
    assert_eq!(distribution.len(), 6); // Two carbon choices, three Br2 compositions.
    for (actual, row) in distribution.peaks().iter().zip(expected) {
        assert_eq!(actual.mass.round(), number(row[0]));
        near(actual.probability, number(row[1]), 1e-7);
        // These additional exact-value oracles are independently computed from
        // source isotope masses and the n=2 binomial, without heap enumeration.
        assert_eq!(number(row[2]), bits64(row[3]));
        near(actual.mass, bits64(row[3]), 1e-10);
        assert_eq!(actual.probability, bits32(row[4]));
    }
    near(distribution.peaks()[0].mass, 169.836_674_2, 1e-10);
    near(distribution.peaks()[1].mass, 170.840_029_2, 1e-10);
}

#[test]
fn relative_threshold_uses_the_carbon_mode_instead_of_the_monoisotopologue() {
    // C100's one-13C/zero-13C probability ratio is 100*p13/p12 > 1.
    // The reverse ratio is about 0.92458, so a 0.95-relative cutoff excludes
    // the lightest configuration. The two-13C state is also below the cutoff.
    let distribution = run("C100", Stop::RelativeThreshold(0.95));
    assert_eq!(distribution.len(), 1);
    near(distribution.peaks()[0].mass, 1_201.003_355, 1e-10);
}

#[test]
fn source_full_carbon_tail_is_selected_before_output_probability_underflow() {
    let distribution = run("C100", Stop::RelativeThreshold(1e-250));
    assert_eq!(distribution.len(), 101);
    assert_eq!(distribution.peaks()[0].mass, 1200.0);
    let last = distribution.peaks().last().unwrap();
    near(last.mass, 1_300.335_500_000_000_1, 1e-10);
    // The roughly 8.67e-198 backend probability passes the nonzero threshold,
    // then narrows to zero in source Peak1D. Its configuration is retained.
    assert_eq!(last.probability, 0.0);
    assert!(distribution.peaks().iter().any(|p| p.probability > 0.0));
}
