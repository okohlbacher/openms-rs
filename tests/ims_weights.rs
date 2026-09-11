// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::IMSWeights;
use openms::chemistry::ims_weights::{MAX_IMS_WEIGHTS, MAX_IMS_WEIGHTS_OUTPUT_BYTES};

const SOURCE_MASSES: [f64; 5] = [71.0456, 180.0312, 1.0186, 4284.36894, 255.0];

#[test]
fn source_all_constructor_copy_access_and_precision_literals() {
    let empty = IMSWeights::new();
    assert_eq!(empty.len(), 0);
    assert!(empty.is_empty());
    assert_eq!(empty.precision(), None);
    let mut weights = IMSWeights::from_masses(&SOURCE_MASSES, 0.01).unwrap();
    for (precision, expected) in [
        (0.01, [7105, 18003, 102, 428437, 25500]),
        (0.1, [710, 1800, 10, 42844, 2550]),
        (1.0, [71, 180, 1, 4284, 255]),
        (0.0001, [710456, 1800312, 10186, 42843689, 2550000]),
        (0.01, [7105, 18003, 102, 428437, 25500]),
    ] {
        weights.set_precision(precision).unwrap();
        assert_eq!(weights.precision(), Some(precision));
        assert_eq!(weights.weights(), expected);
        for i in 0..5 {
            assert_eq!(weights.weight(i).unwrap(), expected[i]);
            assert_eq!(weights.alphabet_mass(i).unwrap(), SOURCE_MASSES[i]);
        }
        assert_eq!(weights.back().unwrap(), expected[4]);
        let cloned = weights.clone();
        assert_eq!(cloned, weights);
        assert_ne!(cloned.masses().as_ptr(), weights.masses().as_ptr());
        assert_ne!(cloned.weights().as_ptr(), weights.weights().as_ptr());
    }
    weights.set_precision(0.00025).unwrap();
    assert_eq!(weights.precision(), Some(0.00025));
}

#[test]
fn source_parent_mass_swap_and_stream_literals() {
    let weights = IMSWeights::from_masses(&SOURCE_MASSES, 0.01).unwrap();
    for i in 0..5 {
        let mut counts = [0; 5];
        counts[i] = 1;
        assert_eq!(weights.parent_mass(&counts).unwrap(), SOURCE_MASSES[i]);
        counts[i] = 2;
        assert_eq!(
            weights.parent_mass(&counts).unwrap(),
            SOURCE_MASSES[i] * 2.0
        );
    }
    assert!(
        weights
            .parent_mass(&[0; 3])
            .unwrap_err()
            .to_string()
            .contains("Expected 5 but got 3.")
    );
    assert_eq!(
        weights.to_text().unwrap(),
        "7105\n18003\n102\n428437\n25500\n"
    );
    let mut swapped = weights.clone();
    swapped.swap(0, 1).unwrap();
    swapped.swap(1, 3).unwrap();
    assert_eq!(
        swapped.masses(),
        [
            SOURCE_MASSES[1],
            SOURCE_MASSES[3],
            SOURCE_MASSES[2],
            SOURCE_MASSES[0],
            SOURCE_MASSES[4]
        ]
    );
    assert_eq!(swapped.weights(), [18003, 428437, 102, 7105, 25500]);
    swapped.set_precision(1.0).unwrap();
    assert_eq!(swapped.weights(), [180, 4284, 1, 71, 255]);
}

#[test]
fn source_gcd_and_relative_error_literals() {
    let mut weights = IMSWeights::from_masses(&[3.0, 5.0, 8.0], 0.1).unwrap();
    assert!(weights.divide_by_gcd().unwrap());
    assert_eq!(weights.weights(), [3, 5, 8]);
    assert_eq!(weights.precision(), Some(1.0));
    assert!(!weights.divide_by_gcd().unwrap());
    let mut primes = IMSWeights::from_masses(&[1.13, 1.67, 2.41], 0.01).unwrap();
    assert!(!primes.divide_by_gcd().unwrap());
    let mut singleton = IMSWeights::from_masses(&[40.0], 0.01).unwrap();
    assert!(!singleton.divide_by_gcd().unwrap());
    let mut values = IMSWeights::from_masses(&SOURCE_MASSES, 0.01).unwrap();
    assert!((values.min_rounding_error().unwrap() - -6.6655113114361e-06).abs() < 1e-15);
    assert!((values.max_rounding_error().unwrap() - 0.00137443549970554).abs() < 1e-15);
    values.set_precision(10.0).unwrap();
    assert_eq!(values.min_rounding_error().unwrap(), -1.0);
    assert!((values.max_rounding_error().unwrap() - 0.0196078431372549).abs() < 1e-15);
}

#[test]
fn gcd_keeps_source_coprime_pair_and_zero_weight_distinctions() {
    for (masses, expected) in [
        (&[3.0, 5.0][..], true),
        (&[3.0, 5.0, 7.0][..], false),
        (&[0.0, 1.0][..], true),
        (&[0.0, 0.0, 1.0][..], false),
    ] {
        let mut weights = IMSWeights::from_masses(masses, 1.0).unwrap();
        let before = weights.clone();
        assert_eq!(weights.divide_by_gcd().unwrap(), expected);
        assert_eq!(weights, before);
    }
    let mut with_zero = IMSWeights::from_masses(&[0.0, 6.0, 9.0], 1.0).unwrap();
    assert!(with_zero.divide_by_gcd().unwrap());
    assert_eq!(with_zero.weights(), [0, 2, 3]);
    assert_eq!(with_zero.precision(), Some(3.0));
    for count in [2, 3] {
        let mut zeros = IMSWeights::from_masses(&vec![0.0; count], 1.0).unwrap();
        let before = zeros.clone();
        assert!(zeros.divide_by_gcd().is_err());
        assert_eq!(zeros, before);
    }
}

#[test]
fn empty_zero_and_signed_finite_source_operations_are_defined_explicitly() {
    let mut empty = IMSWeights::new();
    assert_eq!(empty.parent_mass(&[]).unwrap(), 0.0);
    assert_eq!(empty.min_rounding_error().unwrap(), 0.0);
    assert_eq!(empty.max_rounding_error().unwrap(), 0.0);
    assert_eq!(empty.to_text().unwrap(), "");
    assert!(!empty.divide_by_gcd().unwrap());
    assert!(empty.back().is_err());
    empty.set_precision(0.0).unwrap();
    assert_eq!(empty.precision(), Some(0.0));
    let zeros = IMSWeights::from_masses(&[0.0, -0.0], 2.0).unwrap();
    assert_eq!(zeros.weights(), [0, 0]);
    assert_eq!(zeros.masses()[1].to_bits(), (-0.0f64).to_bits());
    assert_eq!(zeros.min_rounding_error().unwrap(), 0.0); // Source NaN comparisons ignore zero/zero.
    assert_eq!(zeros.max_rounding_error().unwrap(), 0.0);
    let tiny_negative = IMSWeights::from_masses(&[-0.5], 1.0).unwrap();
    assert_eq!(tiny_negative.weights(), [0]);
    assert_eq!(tiny_negative.min_rounding_error().unwrap(), -1.0);
    assert_eq!(tiny_negative.parent_mass(&[2]).unwrap(), -1.0);
    let mut negative = IMSWeights::from_masses(&[-1.0, -2.0], -0.5).unwrap();
    assert_eq!(negative.weights(), [2, 4]);
    assert!(negative.divide_by_gcd().unwrap());
    assert_eq!(negative.weights(), [1, 2]);
    assert_eq!(negative.precision(), Some(-1.0));
    assert_eq!(negative.parent_mass(&[2, 1]).unwrap(), -4.0);
}

#[test]
fn literal_half_rounding_and_full_u64_boundary() {
    let half = 2.5f64;
    let masses = [
        f64::from_bits(half.to_bits() - 1),
        half,
        f64::from_bits(half.to_bits() + 1),
    ];
    assert_eq!(
        IMSWeights::from_masses(&masses, 1.0).unwrap().weights(),
        [2, 3, 3]
    );
    let end = 18_446_744_073_709_551_616.0f64;
    let below = f64::from_bits(end.to_bits() - 1);
    let weights = IMSWeights::from_masses(&[below], 1.0).unwrap();
    assert_eq!(weights.weight(0).unwrap(), 18_446_744_073_709_549_568);
    assert_eq!(weights.alphabet_mass(0).unwrap().to_bits(), below.to_bits());
    assert!(IMSWeights::from_masses(&[end], 1.0).is_err());
    assert!(IMSWeights::from_masses(&[f64::MAX], 1.0).is_err());
}

#[test]
fn failed_mutations_and_arithmetic_overflow_are_checked_atomically() {
    let mut weights = IMSWeights::from_masses(&[0.001, 100.0], 1.0).unwrap();
    let old = weights.clone();
    let pointer = weights.weights().as_ptr();
    for precision in [0.0, -1.0, f64::INFINITY, f64::NAN, f64::MIN_POSITIVE] {
        assert!(weights.set_precision(precision).is_err());
        assert_eq!(weights, old);
        assert_eq!(weights.weights().as_ptr(), pointer);
    }
    assert!(weights.swap(0, 2).is_err());
    assert_eq!(weights, old);
    assert!(weights.weight(2).is_err());
    assert!(weights.alphabet_mass(usize::MAX).is_err());
    for mass in [f64::INFINITY, f64::NAN, -1.0] {
        assert!(IMSWeights::from_masses(&[mass], 1.0).is_err());
    }
    let huge = IMSWeights::from_masses(&[f64::MAX], f64::MAX).unwrap();
    assert!(huge.parent_mass(&[2]).is_err());
    let precision = f64::MAX / 1.5;
    let mut huge = IMSWeights::from_masses(&[f64::MAX, f64::MAX], precision).unwrap();
    assert_eq!(huge.weights(), [2, 2]);
    let saved = huge.clone();
    assert!(huge.divide_by_gcd().is_err());
    assert_eq!(huge, saved);
    assert_eq!(huge.min_rounding_error().unwrap(), 0.0); // Source ignores the positive infinity here.
    assert!(huge.max_rounding_error().is_err());
}

#[test]
fn independent_integer_oracles_parent_order_and_resource_limits() {
    // Quarter masses and integer precisions: independently use exact integers.
    for precision in 1..=7u64 {
        let masses: Vec<_> = (0..400u64).map(|n| n as f64 / 4.0).collect();
        let expected: Vec<_> = (0..400u64)
            .map(|n| (n + 2 * precision) / (4 * precision))
            .collect();
        assert_eq!(
            IMSWeights::from_masses(&masses, precision as f64)
                .unwrap()
                .weights(),
            expected
        );
    }
    let mut weights = IMSWeights::from_masses(&[9_007_199_254_740_992.0, 1.0, 1.0], 1.0).unwrap();
    assert_eq!(
        weights.parent_mass(&[1, 1, 1]).unwrap(),
        9_007_199_254_740_992.0
    );
    weights.swap(0, 2).unwrap();
    assert_eq!(
        weights.parent_mass(&[1, 1, 1]).unwrap(),
        9_007_199_254_740_994.0
    );
    assert!(IMSWeights::from_masses(&vec![1.0; MAX_IMS_WEIGHTS + 1], 1.0).is_err());
    let count = MAX_IMS_WEIGHTS_OUTPUT_BYTES / 21 + 1;
    let mass = f64::from_bits(18_446_744_073_709_551_616.0f64.to_bits() - 1);
    let weights = IMSWeights::from_masses(&vec![mass; count], 1.0).unwrap();
    assert!(weights.to_text().is_err());
    assert_eq!(weights.len(), count);
}
