// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `LinearRegressionWithoutIntercept_test.cpp` at core pin `bc9cc12`
//! (`src/tests/class_tests/openms/source/`, evidence tier 3), section by
//! section, plus the refusal that replaces the source's unchecked read.

use openms::math::ml::linear_regression_without_intercept::LinearRegressionWithoutIntercept;

/// ClassTest `isRealSimilar`: absolute difference at most 1e-5, or ratio at
/// most 1 + 1e-5.
fn similar(actual: f64, expected: f64) {
    let ratio = if actual.abs() > expected.abs() {
        actual / expected
    } else {
        expected / actual
    };
    assert!(
        (actual - expected).abs() <= 1e-5 || (ratio - 1.0).abs() <= 1e-5,
        "{actual} !~ {expected}"
    );
}

#[test]
fn default_constructor_has_no_slope() {
    // "no data yet -> slope is NaN"
    assert!(LinearRegressionWithoutIntercept::new().slope().is_nan());
    assert!(LinearRegressionWithoutIntercept::default().slope().is_nan());
}

#[test]
fn add_data_single_observations() {
    let mut lr = LinearRegressionWithoutIntercept::new();
    // a single observation is not enough -> NaN
    lr.add_data(1.0, 2.0);
    assert!(lr.slope().is_nan());

    // perfectly collinear through the origin: y = 2x;
    // slope = sum(x*y) / sum(x*x) = 28 / 14 = 2
    lr.add_data(2.0, 4.0);
    lr.add_data(3.0, 6.0);
    similar(lr.slope(), 2.0);
}

#[test]
fn add_data_vectors() {
    let mut lr = LinearRegressionWithoutIntercept::new();
    let x = [1.0, 2.0];
    let y = [2.1, 3.9];
    lr.add_data_slices(&x, &y).unwrap();
    // slope = (1*2.1 + 2*3.9) / (1*1 + 2*2) = 9.9 / 5 = 1.98
    similar(lr.slope(), 1.98);
}

#[test]
fn get_slope() {
    // fewer than two observations -> NaN
    let lr0 = LinearRegressionWithoutIntercept::new();
    assert!(lr0.slope().is_nan());

    let mut lr1 = LinearRegressionWithoutIntercept::new();
    lr1.add_data(5.0, 10.0);
    assert!(lr1.slope().is_nan());

    // two observations, scaled line y = 3x
    let mut lr2 = LinearRegressionWithoutIntercept::new();
    lr2.add_data(1.0, 3.0);
    lr2.add_data(2.0, 6.0);
    similar(lr2.slope(), 3.0);
}

#[test]
fn shorter_dependent_vector_is_refused_and_adds_nothing() {
    // Tier 4: the source reads `y[i]` past the end; this port refuses and
    // leaves the accumulator unchanged.
    let mut lr = LinearRegressionWithoutIntercept::new();
    lr.add_data(1.0, 1.0);
    assert!(lr.add_data_slices(&[1.0, 2.0, 3.0], &[1.0, 2.0]).is_err());
    assert!(lr.slope().is_nan(), "still one observation");
    // A longer y is ignored beyond x, as the source's loop over x does.
    lr.add_data_slices(&[2.0], &[2.0, 99.0]).unwrap();
    similar(lr.slope(), 1.0);
}
