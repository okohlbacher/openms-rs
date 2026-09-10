// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Independent source and closed-form references; no C++ execution.
//! Provenance: data/fine_isotope_stream_provenance.json.

use openms::chemistry::{
    EmpiricalFormula, FineIsotopeConfiguration as Configuration, FineIsotopeIterator as Isotopes,
    FineIsotopePatternGenerator as Generator, FineIsotopeStop as Stop,
};

const FRUCTOSE: &str = include_str!("data/fine_isotope_fructose.tsv");

fn formula(text: &str) -> EmpiricalFormula {
    text.parse().unwrap()
}
fn source_relative(actual: f64, expected: f64) {
    // IsoSpec_test.cpp:32-34, with no absolute floor.
    const EPSILON: f64 = 0.000_000_1;
    assert!(
        actual * (1.0 - EPSILON) <= expected && expected <= actual * (1.0 + EPSILON),
        "{actual:.17e} != {expected:.17e} within source relative tolerance"
    );
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; absolute tolerance {tolerance}"
    );
}
fn check_accessors(configuration: &Configuration) {
    let peak = configuration.to_peak().unwrap();
    assert_eq!(peak.mz, configuration.mass);
    assert_eq!(peak.intensity, configuration.probability as f32);
    source_relative(
        configuration.probability,
        configuration.log_probability.exp(),
    );
}
fn fructose_reference() -> Vec<(f64, f64)> {
    let mut expected: Vec<_> = FRUCTOSE
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .skip(1)
        .map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            let mass = f64::from_bits(u64::from_str_radix(fields[2], 16).unwrap());
            let probability =
                f64::from(f32::from_bits(u32::from_str_radix(fields[3], 16).unwrap()));
            (mass, probability)
        })
        .collect();
    expected.sort_by(|a, b| b.1.total_cmp(&a.1));
    assert_eq!(expected.len(), 14);
    expected
}
fn check_fructose_prefix(actual: &[Configuration]) {
    for (actual, (mass, probability)) in actual.iter().zip(fructose_reference()) {
        source_relative(actual.mass, mass);
        source_relative(actual.probability, probability);
        check_accessors(actual);
    }
}

// Source ElementDB.cpp literals, deliberately independent of Rust element().
// Rows use C/H/N/O/S order. Decimal strings preserve the original literals.
const MASSES: [&[&str]; 5] = [
    &["12.0", "13.003355000000001"],
    &["1.0078250319", "2.01410178"],
    &["14.003074", "15.000109"],
    &["15.994915000000001", "16.999132", "17.999168999999998"],
    &[
        "31.972070729999999",
        "32.971457999999998",
        "33.967866999999998",
        "35.967081",
    ],
];
const ABUNDANCES: [&[&str]; 5] = [
    &["0.9893000000000001", "0.010700000000000001"],
    &["0.999885", "0.000115"],
    &["0.9963200000000001", "0.00368"],
    &[
        "0.9975700000000001",
        "0.00037999999999999997",
        "0.0020499999999999997",
    ],
    &["0.9493", "0.0076", "0.0429", "0.0002"],
];
fn natural_tables(counts: [u32; 5]) -> (Vec<u32>, Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let mut active_counts = Vec::new();
    let mut masses = Vec::new();
    let mut probabilities = Vec::new();
    for (index, count) in counts.into_iter().enumerate() {
        if count != 0 {
            active_counts.push(count);
            masses.push(MASSES[index].iter().map(|x| x.parse().unwrap()).collect());
            probabilities.push(
                ABUNDANCES[index]
                    .iter()
                    .map(|x| f64::from(x.parse::<f64>().unwrap() as f32))
                    .collect(),
            );
        }
    }
    (active_counts, masses, probabilities)
}
fn custom_natural(counts: [u32; 5]) -> Isotopes {
    let (counts, masses, probabilities) = natural_tables(counts);
    Isotopes::from_isotopes(&counts, &masses, &probabilities).unwrap()
}

#[test]
fn source_ordered_fructose_prefix_and_complete_support_from_both_inputs() {
    for mut stream in [
        Isotopes::from_formula(&formula("C6H12O6")).unwrap(),
        custom_natural([6, 12, 0, 6, 0]),
    ] {
        let mut prefix = Vec::new();
        let mut previous_log = f64::INFINITY;
        let mut previous_probability = f64::INFINITY;
        let mut count = 0;
        for item in &mut stream {
            let configuration = item.unwrap();
            check_accessors(&configuration);
            assert!(configuration.log_probability <= previous_log);
            assert!(configuration.probability <= previous_probability);
            previous_log = configuration.log_probability;
            previous_probability = configuration.probability;
            if count < 14 {
                prefix.push(configuration);
            }
            count += 1;
        }
        assert_eq!(count, 2_548); // Exact ordered-generator source assertion.
        assert_eq!(prefix.len(), 14);
        check_fructose_prefix(&prefix);
        assert!(stream.next().is_none());
        assert!(stream.next().is_none());
    }
}

#[test]
fn source_insulin_ten_thousand_configuration_prefix_retains_raw_accessors() {
    let mut stream = Isotopes::from_formula(&formula("C520H817N139O147S8")).unwrap();
    let mut previous_log = f64::INFINITY;
    let mut previous_probability = f64::INFINITY;
    for _ in 0..10_000 {
        let configuration = stream.next().unwrap().unwrap();
        check_accessors(&configuration);
        assert!(configuration.log_probability <= previous_log);
        assert!(configuration.probability <= previous_probability);
        previous_log = configuration.log_probability;
        previous_probability = configuration.probability;
    }
    // Upstream only provides accessor checks for this prefix, not literal masses.
    // Its helper's extra uninspected advance is not a scientific assertion.
}

#[test]
fn source_threshold_counts_hold_for_formula_and_custom_natural_inputs() {
    for (text, counts, absolute, threshold, expected) in [
        ("C6H12O6", [6, 12, 0, 6, 0], false, 1e-5, 14),
        ("C6H12O6", [6, 12, 0, 6, 0], true, 1e-5, 14),
        (
            "C520H817N139O147S8",
            [520, 817, 139, 147, 8],
            false,
            1e-5,
            5_513,
        ),
        (
            "C520H817N139O147S8",
            [520, 817, 139, 147, 8],
            false,
            0.01,
            267,
        ),
        (
            "C520H817N139O147S8",
            [520, 817, 139, 147, 8],
            true,
            1e-5,
            1_734,
        ),
        (
            "C520H817N139O147S8",
            [520, 817, 139, 147, 8],
            true,
            0.01,
            21,
        ),
    ] {
        for stream in [
            Isotopes::from_formula(&formula(text)).unwrap(),
            custom_natural(counts),
        ] {
            let selected = if absolute {
                stream.with_absolute_threshold(threshold).unwrap()
            } else {
                stream.with_relative_threshold(threshold).unwrap()
            };
            let actual: Vec<_> = selected.collect::<openms::Result<_>>().unwrap();
            assert_eq!(
                actual.len(),
                expected,
                "{text}, absolute={absolute}, {threshold}"
            );
            for configuration in &actual {
                check_accessors(configuration);
            }
            if text == "C6H12O6" {
                check_fructose_prefix(&actual);
            }
        }
    }
}

#[test]
fn source_custom_natural_materialization_keeps_source_output_rounding() {
    let (counts, masses, probabilities) = natural_tables([6, 12, 0, 6, 0]);
    for (stop, count) in [
        (Stop::RelativeThreshold(1e-5), 14),
        (Stop::AbsoluteThreshold(1e-5), 14),
        (Stop::UnexplainedProbability(0.000_01), 17),
    ] {
        let actual = Generator { stop }
            .run_with_isotopes(&counts, &masses, &probabilities)
            .unwrap();
        assert_eq!(actual.len(), count);
        let mut actual = actual.peaks().to_vec();
        actual.sort_by(|a, b| b.probability.total_cmp(&a.probability));
        for (actual, (mass, probability)) in actual.iter().zip(fructose_reference()) {
            source_relative(actual.mass, mass);
            source_relative(actual.probability, probability);
            assert_eq!(actual.probability, f64::from(actual.probability as f32));
        }
    }
}

#[test]
fn custom_binary64_weights_are_not_rounded_before_enumeration() {
    let probabilities = [vec![0.500_000_000_000_000_1, 0.499_999_999_999_999_9]];
    let configurations: Vec<_> = Isotopes::from_isotopes(&[1], &[vec![10.0, 20.0]], &probabilities)
        .unwrap()
        .collect::<openms::Result<_>>()
        .unwrap();
    assert_eq!(configurations.len(), 2);
    assert_eq!(configurations[0].mass, 10.0);
    assert_eq!(configurations[1].mass, 20.0);
    assert!(configurations[0].probability > 0.5);
    assert!(configurations[1].probability < 0.5);
    assert!(configurations[0].log_probability > configurations[1].log_probability);
    for configuration in &configurations {
        check_accessors(configuration);
        assert_eq!(configuration.to_peak().unwrap().intensity, 0.5);
    }
}

#[test]
fn independent_two_component_multinomial_retains_equal_mass_configurations() {
    let mut configurations: Vec<_> = Isotopes::from_isotopes(
        &[2, 1],
        &[vec![10.0, 11.0], vec![100.0, 102.0]],
        &[vec![0.25, 0.75], vec![0.6, 0.4]],
    )
    .unwrap()
    .collect::<openms::Result<_>>()
    .unwrap();
    assert_eq!(configurations.len(), 6);
    assert!(
        configurations
            .windows(2)
            .all(|p| p[0].probability >= p[1].probability)
    );
    // Hand product of the three n=2 binomial configurations with two n=1
    // configurations. Equal .225 probabilities do not impose a source tie order.
    configurations.sort_by(|a, b| {
        a.mass
            .total_cmp(&b.mass)
            .then(a.probability.total_cmp(&b.probability))
    });
    let expected = [
        (120.0, 0.0375),
        (121.0, 0.225),
        (122.0, 0.025),
        (122.0, 0.3375),
        (123.0, 0.15),
        (124.0, 0.225),
    ];
    for (actual, (mass, probability)) in configurations.iter().zip(expected) {
        near(actual.mass, mass, 1e-12);
        near(actual.probability, probability, 1e-15);
        near(actual.log_probability, probability.ln(), 1e-14);
        check_accessors(actual);
    }
}
