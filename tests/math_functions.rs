// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//
// Expected values marked "source literal" are transcribed from the pinned
// MathFunctions_test.cpp at revision bc9cc12514c768385ce121d6ca4bb710fe1983c4
// (tier 3, source review: no C++ was built or executed). Values marked
// "independently derived" are computed here from a closed form, an exact
// rational evaluation or an algebraic invariant, and are stronger evidence.

use openms::Error;
use openms::concept::math_functions::{
    Bin, MAX_BINOMIAL_TRIALS, MAX_BINS, MAX_DECIMAL_DIGITS, approximately_equal,
    binomial_cdf_complement, ceil_decimal, contains, create_bins, extend_range, extended_gcd, gcd,
    interval_transformation, is_odd, linear_to_log10, log_binomial_coef, log_sum_exp,
    log10_to_linear, percent_of, ppm, ppm_abs, ppm_to_mass, ppm_to_mass_abs, quantile, round,
    round_decimal, round_to, tolerance_window, zoom_in,
};

/// The class test's default `TOLERANCE_RELATIVE(1.0 + 1e-5)`, with a unit-sized
/// absolute floor so that comparisons against zero are meaningful.
const SOURCE_TOLERANCE: f64 = 1e-5;

#[track_caller]
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "{actual} != {expected} within {relative}"
    );
}

// START_SECTION((double log_binomial_coef(unsigned n, unsigned k)))
#[test]
fn log_binomial_coef_matches_source_values_and_is_exactly_symmetric() {
    // Source literals.
    close(log_binomial_coef(10, 5).unwrap(), 5.5294, SOURCE_TOLERANCE);
    close(
        log_binomial_coef(20, 10).unwrap(),
        12.12679,
        SOURCE_TOLERANCE,
    );
    assert_eq!(log_binomial_coef(5, 0).unwrap(), 0.0);
    assert_eq!(log_binomial_coef(5, 5).unwrap(), 0.0);
    // The source asserts approximate symmetry; the `k > n/2` swap makes the two
    // calls evaluate the identical expression, so equality is exact.
    assert_eq!(
        log_binomial_coef(10, 3).unwrap(),
        log_binomial_coef(10, 7).unwrap()
    );
    // Independently derived: ln(C(10,5)) = ln(252), ln(C(20,10)) = ln(184756).
    close(log_binomial_coef(10, 5).unwrap(), 252.0_f64.ln(), 1e-14);
    close(
        log_binomial_coef(20, 10).unwrap(),
        184_756.0_f64.ln(),
        1e-14,
    );
    // Source throws std::invalid_argument; this port returns InvalidValue.
    assert!(matches!(
        log_binomial_coef(5, 6),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(log_binomial_coef(0, 0).unwrap(), 0.0);
}

// START_SECTION((double log_sum_exp(double x, double y)))
#[test]
fn log_sum_exp_is_stable_and_treats_negative_infinity_as_an_identity() {
    // Source literals.
    close(log_sum_exp(1.0, 2.0), 2.31326169, SOURCE_TOLERANCE);
    close(log_sum_exp(10.0, 10.0), 10.6931472, SOURCE_TOLERANCE);
    close(log_sum_exp(100.0, 0.0), 100.0, SOURCE_TOLERANCE);
    close(log_sum_exp(0.0, 100.0), 100.0, SOURCE_TOLERANCE);
    assert_eq!(log_sum_exp(f64::NEG_INFINITY, 5.0), 5.0);
    assert_eq!(log_sum_exp(5.0, f64::NEG_INFINITY), 5.0);
    // Independently derived: ln(e^x + e^x) = x + ln 2, exactly representable
    // as an addition here, and ln(e^1 + e^2) = 2 + ln(1 + e^-1).
    assert_eq!(log_sum_exp(10.0, 10.0), 10.0 + 2.0_f64.ln());
    close(log_sum_exp(1.0, 2.0), 2.0 + (-1.0_f64).exp().ln_1p(), 1e-15);
    // Factoring out the maximum is what keeps large arguments finite; the naive
    // form would overflow here.
    assert!(log_sum_exp(1000.0, 1000.0).is_finite());
    assert_eq!(log_sum_exp(1000.0, 1000.0), 1000.0 + 2.0_f64.ln());
}

// START_SECTION((ceilDecimal))
#[test]
fn ceil_decimal_shifts_ceils_and_shifts_back() {
    // Source literals. Non-negative powers are exact: pow(10, 0..2) is exact and
    // the ceiling of an exactly scaled value is integral.
    assert_eq!(ceil_decimal(12345.67, 0).unwrap(), 12346.0);
    assert_eq!(ceil_decimal(12345.67, 1).unwrap(), 12350.0);
    assert_eq!(ceil_decimal(12345.67, 2).unwrap(), 12400.0);
    // Negative powers carry the rounding of pow(10, -n).
    close(ceil_decimal(12345.671, -2).unwrap(), 12345.68, 1e-12);
    close(ceil_decimal(12345.67, -1).unwrap(), 12345.7, 1e-12);
    // Independently derived: the header's own worked example.
    assert_eq!(ceil_decimal(123.0, 1).unwrap(), 130.0);
    assert_eq!(ceil_decimal(123.0, 2).unwrap(), 200.0);
    close(ceil_decimal(0.123, -2).unwrap(), 0.13, 1e-12);
    // Rounds up, never to nearest: it is not roundDecimal.
    assert_eq!(ceil_decimal(12345.01, 2).unwrap(), 12400.0);
    assert!(ceil_decimal(f64::NAN, 0).is_err());
    assert!(ceil_decimal(1.0, 400).is_err());
    assert!(ceil_decimal(1.0, -400).is_err());
}

// START_SECTION((roundDecimal))
#[test]
fn round_decimal_rounds_to_nearest_and_keeps_the_negative_zero_branch() {
    // Source literals.
    assert_eq!(round_decimal(12345.67, 0).unwrap(), 12346.0);
    assert_eq!(round_decimal(12345.67, 1).unwrap(), 12350.0);
    assert_eq!(round_decimal(12345.67, 2).unwrap(), 12300.0);
    close(round_decimal(12345.671, -2).unwrap(), 12345.67, 1e-12);
    close(round_decimal(12345.67, -1).unwrap(), 12345.7, 1e-12);
    // Independently derived: the header's worked example, and the mirrored
    // magnitude for negatives, which takes the second branch.
    assert_eq!(round_decimal(123.0, 1).unwrap(), 120.0);
    assert_eq!(round_decimal(123.0, 2).unwrap(), 100.0);
    assert_eq!(round_decimal(-12345.67, 2).unwrap(), -12300.0);
    assert_eq!(
        round_decimal(-12345.67, 2).unwrap(),
        -round_decimal(12345.67, 2).unwrap()
    );
    // Zero is not > 0, so it takes the negating branch and returns -0.0.
    let zero = round_decimal(0.0, 0).unwrap();
    assert_eq!(zero, 0.0);
    assert!(zero.is_sign_negative());
    assert_eq!(1.0_f64.copysign(zero), -1.0);
}

// START_SECTION((intervalTransformation))
#[test]
fn interval_transformation_maps_between_two_spans() {
    // Source literals; all five are exactly representable.
    assert_eq!(
        interval_transformation(0.5, 0.0, 1.0, 0.0, 600.0).unwrap(),
        300.0
    );
    assert_eq!(
        interval_transformation(0.5, 0.25, 1.0, 0.0, 600.0).unwrap(),
        200.0
    );
    assert_eq!(
        interval_transformation(0.5, 0.0, 0.75, 0.0, 600.0).unwrap(),
        400.0
    );
    assert_eq!(
        interval_transformation(0.5, 0.0, 1.0, 150.0, 600.0).unwrap(),
        375.0
    );
    assert_eq!(
        interval_transformation(0.5, 0.0, 1.0, 0.0, 450.0).unwrap(),
        225.0
    );
    // Independently derived: endpoints map to endpoints, for any spans.
    assert_eq!(
        interval_transformation(3.0, 3.0, 7.0, -2.0, 11.0).unwrap(),
        -2.0
    );
    assert_eq!(
        interval_transformation(7.0, 3.0, 7.0, -2.0, 11.0).unwrap(),
        11.0
    );
    // The source divides by right1 - left1 unguarded.
    assert!(interval_transformation(0.5, 1.0, 1.0, 0.0, 600.0).is_err());
    assert!(interval_transformation(f64::NAN, 0.0, 1.0, 0.0, 1.0).is_err());
}

// START_SECTION((linear2log))
#[test]
fn linear_to_log10_adds_one_before_taking_the_logarithm() {
    // Source literals, all exact powers of ten after the shift.
    assert_eq!(linear_to_log10(0.0).unwrap(), 0.0);
    assert_eq!(linear_to_log10(9.0).unwrap(), 1.0);
    assert_eq!(linear_to_log10(99.0).unwrap(), 2.0);
    assert_eq!(linear_to_log10(999.0).unwrap(), 3.0);
    // The source returns -inf at x == -1 and NaN below; this port refuses.
    assert!(linear_to_log10(-1.0).is_err());
    assert!(linear_to_log10(-2.0).is_err());
    assert!(linear_to_log10(f64::INFINITY).is_err());
}

// START_SECTION((log2linear))
#[test]
fn log10_to_linear_inverts_linear_to_log10() {
    // Source literals.
    assert_eq!(log10_to_linear(0.0).unwrap(), 0.0);
    assert_eq!(log10_to_linear(1.0).unwrap(), 9.0);
    assert_eq!(log10_to_linear(2.0).unwrap(), 99.0);
    assert_eq!(log10_to_linear(3.0).unwrap(), 999.0);
    // Independently derived: the pair round-trips on the source's own literals.
    for x in [0.0, 9.0, 99.0, 999.0, 1.5, 12345.0] {
        close(
            log10_to_linear(linear_to_log10(x).unwrap()).unwrap(),
            x,
            1e-12,
        );
    }
    assert!(log10_to_linear(400.0).is_err());
}

// START_SECTION((isOdd))
#[test]
fn is_odd_tests_the_low_bit() {
    // Source literals.
    assert!(!is_odd(0));
    assert!(is_odd(1));
    assert!(!is_odd(2));
    assert!(is_odd(3));
    // Independently derived: the predicate is the low bit for every width.
    assert!(is_odd(u32::MAX));
    assert!(!is_odd(u32::MAX - 1));
}

// START_SECTION((template <typename T> T round (T x)))
#[test]
fn round_breaks_ties_away_from_zero() {
    // Source literals. The first two are float in the class test; the float
    // instantiation maps to f32::round, which is asserted alongside.
    assert_eq!(round(f64::from(14.49_f32)), 14.0);
    assert_eq!(round(f64::from(14.50_f32)), 15.0);
    assert_eq!(14.49_f32.round(), 14.0);
    assert_eq!(14.50_f32.round(), 15.0);
    assert_eq!(round(-999.49), -999.0);
    assert_eq!(round(-675.77), -676.0);
    // Independently derived: halves go away from zero, not to even.
    assert_eq!(round(0.5), 1.0);
    assert_eq!(round(-0.5), -1.0);
    assert_eq!(round(1.5), 2.0);
    assert_eq!(round(2.5), 3.0);
    assert!(round(f64::NAN).is_nan());
}

// START_SECTION(template<typename T> T roundTo(const T value, int digits))
#[test]
#[allow(clippy::approx_constant)] // The source's own literal is 3.14159265.
fn round_to_uses_the_source_scaling_loop() {
    // Source literals.
    close(round_to(3.14159265, 2).unwrap(), 3.14, 1e-12);
    close(round_to(1234.9, -2).unwrap(), 1200.0, SOURCE_TOLERANCE);
    assert_eq!(round_to(1234.9, 0).unwrap(), 1235.0);
    close(round_to(1234.9, -1).unwrap(), 1230.0, SOURCE_TOLERANCE);
    close(round_to(1234.9, -3).unwrap(), 1000.0, SOURCE_TOLERANCE);
    // Independently derived: the factor is the repeated-division value, not
    // 10f64.powi(digits). Reproducing the loop here reproduces the result bit
    // for bit; the powi form is only equal to within the source tolerance.
    let mut factor = 1.0_f64;
    for _ in 0..3 {
        factor /= 10.0;
    }
    assert_eq!(
        round_to(1234.9, -3).unwrap(),
        (1234.9 * factor).round() / factor
    );
    // Zero digits is the identity on integers and rounds halves away from zero.
    assert_eq!(round_to(-2.5, 0).unwrap(), -3.0);
    assert!(round_to(f64::NAN, 2).is_err());
    assert!(round_to(1.0, 320).is_err());
    // The source's factor loop runs |digits| times, so i32::MIN would spin over
    // two billion iterations; MAX_DECIMAL_DIGITS refuses before the loop.
    assert!(round_to(1.0, MAX_DECIMAL_DIGITS as i32 + 1).is_err());
    assert!(round_to(1.0, i32::MIN).is_err());
}

// START_SECTION(template<typename T> double percentOf(T value, T total, int digits))
#[test]
fn percent_of_rounds_and_refuses_negative_arguments() {
    // Source literals.
    close(percent_of(1.0 / 3.0, 1.0, 2).unwrap(), 33.33, 1e-12);
    close(percent_of(1.0 / 3.0, 1.0, 3).unwrap(), 33.333, 1e-12);
    close(percent_of(1.0 / 3.0, 1.0, 4).unwrap(), 33.3333, 1e-12);
    close(percent_of(166.6666, 1000.0, 1).unwrap(), 16.7, 1e-12);
    // Source throws Exception::InvalidValue for either negative argument.
    assert!(matches!(
        percent_of(-1.0, 1000.0, 2),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        percent_of(1.0, -1000.0, 2),
        Err(Error::InvalidValue(_))
    ));
    // The documented zero-total shortcut: 0.0, not a division by zero.
    assert_eq!(percent_of(5.0, 0.0, 2).unwrap(), 0.0);
    assert_eq!(percent_of(0.0, 0.0, 2).unwrap(), 0.0);
    // Independently derived: a part equal to the whole is one hundred percent.
    assert_eq!(percent_of(7.0, 7.0, 6).unwrap(), 100.0);
    assert!(percent_of(f64::NAN, 1.0, 2).is_err());
}

// START_SECTION((bool approximatelyEqual(double a, double b, double tol)))
#[test]
fn approximately_equal_is_an_absolute_comparison() {
    // Source literals.
    assert!(approximately_equal(1.1, 1.1002, 0.1));
    assert!(approximately_equal(1.1, 1.1002, 0.01));
    assert!(approximately_equal(1.1, 1.1002, 0.001));
    assert!(!approximately_equal(1.1, 1.1002, 0.0001));
    // Independently derived: the tolerance is absolute, so the same relative
    // difference fails at a larger magnitude.
    assert!(!approximately_equal(1.1e6, 1.1002e6, 0.1));
    assert!(approximately_equal(5.0, 5.0, 0.0));
    assert!(!approximately_equal(f64::NAN, f64::NAN, 1e300));
}

// START_SECTION((template <typename T> T getPPM(T mz_obs, T mz_ref)))
#[test]
fn ppm_divides_by_the_reference_and_is_not_antisymmetric() {
    // Source literals: 1/1000*1e6 and -1/1000*1e6.
    assert_eq!(ppm(1001.0, 1000.0).unwrap(), 1000.0);
    assert_eq!(ppm(999.0, 1000.0).unwrap(), -1000.0);
    // Independently derived asymmetry: the divisor is always the *reference*,
    // so swapping the arguments is not a sign flip. -1/1001*1e6 is about
    // -999.001, which is 1 ppm away from the naive -1000.
    let swapped = ppm(1000.0, 1001.0).unwrap();
    close(swapped, -1.0 / 1001.0 * 1e6, 1e-12);
    assert_ne!(swapped, -ppm(1001.0, 1000.0).unwrap());
    assert!((swapped + 1000.0).abs() > 0.99);
    // The observed value is the numerator only; equal inputs give exactly zero.
    assert_eq!(ppm(1234.5, 1234.5).unwrap(), 0.0);
    assert!(ppm(1000.0, 0.0).is_err());
    assert!(ppm(f64::NAN, 1000.0).is_err());
}

// START_SECTION((template <typename T> T getPPMAbs(T mz_obs, T mz_ref)))
#[test]
fn ppm_abs_drops_the_sign_only() {
    // Source literals.
    assert_eq!(ppm_abs(1001.0, 1000.0).unwrap(), 1000.0);
    assert_eq!(ppm_abs(999.0, 1000.0).unwrap(), 1000.0);
    // Independently derived: it is exactly |ppm|, for both signs.
    for (obs, theo) in [(500.001, 500.0), (499.999, 500.0), (2000.0, 1999.0)] {
        assert_eq!(ppm_abs(obs, theo).unwrap(), ppm(obs, theo).unwrap().abs());
        assert!(ppm_abs(obs, theo).unwrap() >= 0.0);
    }
    assert!(ppm_abs(1.0, 0.0).is_err());
}

// START_SECTION((pair<double, double> getTolWindow(double val, double tol, bool ppm)))
#[test]
fn tolerance_window_is_asymmetric_in_ppm_mode_and_symmetric_in_dalton_mode() {
    // Source literals.
    let (left, right) = tolerance_window(1000.0, 10.0, true).unwrap();
    close(left, 999.99, 1e-12);
    close(right, 1000.0100001, SOURCE_TOLERANCE);
    let (left_da, right_da) = tolerance_window(1000.0, 10.0, false).unwrap();
    assert_eq!(left_da, 990.0);
    assert_eq!(right_da, 1010.0);
    let (left5, right5) = tolerance_window(500.0, 5.0, true).unwrap();
    close(left5, 499.9975, 1e-12);
    close(right5, 500.0025000125, SOURCE_TOLERANCE);

    // Independently derived: the ppm window is wider to the right, and the
    // right edge is exactly the largest x whose own ppm window still holds val.
    assert!(right - 1000.0 > 1000.0 - left);
    close(ppm(1000.0, right).unwrap(), -10.0, 1e-9);
    close(ppm(left, 1000.0).unwrap(), -10.0, 1e-9);
    // The Dalton window is symmetric by construction.
    assert_eq!(right_da - 1000.0, 1000.0 - left_da);
    // tol == 1e6 ppm makes the source's denominator zero.
    assert!(tolerance_window(1000.0, 1e6, true).is_err());
    assert!(tolerance_window(1000.0, f64::NAN, false).is_err());
}

// START_SECTION((double binomial_cdf_complement(unsigned N, unsigned n, double p)))
#[test]
fn binomial_cdf_complement_matches_source_values_and_an_exact_rational_oracle() {
    // Source literals.
    close(
        binomial_cdf_complement(10, 5, 0.5).unwrap(),
        0.623046875,
        SOURCE_TOLERANCE,
    );
    close(
        binomial_cdf_complement(20, 10, 0.4).unwrap(),
        0.24466,
        SOURCE_TOLERANCE,
    );
    assert_eq!(binomial_cdf_complement(10, 0, 0.3).unwrap(), 1.0);
    close(
        binomial_cdf_complement(10, 10, 0.7).unwrap(),
        0.0282475249,
        SOURCE_TOLERANCE,
    );
    assert_eq!(binomial_cdf_complement(10, 0, 0.0).unwrap(), 1.0);
    assert_eq!(binomial_cdf_complement(10, 1, 0.0).unwrap(), 0.0);
    assert_eq!(binomial_cdf_complement(10, 0, 1.0).unwrap(), 1.0);
    assert_eq!(binomial_cdf_complement(10, 10, 1.0).unwrap(), 1.0);
    // The AScore cases from the source's own comment.
    close(
        binomial_cdf_complement(1, 1, 0.1).unwrap(),
        0.1,
        SOURCE_TOLERANCE,
    );
    close(
        binomial_cdf_complement(3, 1, 0.1).unwrap(),
        0.271,
        SOURCE_TOLERANCE,
    );
    close(
        binomial_cdf_complement(100, 60, 0.5).unwrap(),
        0.02844,
        SOURCE_TOLERANCE,
    );
    // Source throws std::invalid_argument for all three.
    assert!(binomial_cdf_complement(10, 11, 0.5).is_err());
    assert!(binomial_cdf_complement(10, 5, -0.1).is_err());
    assert!(binomial_cdf_complement(10, 5, 1.1).is_err());
    // The source's p range test passes a NaN straight into Boost; this refuses.
    assert!(binomial_cdf_complement(10, 5, f64::NAN).is_err());
    assert!(binomial_cdf_complement(MAX_BINOMIAL_TRIALS + 1, 1, 0.5).is_err());

    // Independently derived: P(X >= n) for B(10, 1/2) is the dyadic rational
    // sum_{k=5}^{10} C(10,k) / 2^10 = (252+210+120+45+10+1)/1024 = 638/1024.
    // The log-domain sum reproduces it to about 3.4e-15 relative, four orders
    // of magnitude inside the source test's own tolerance; it is not bit-exact
    // because the terms pass through lgamma.
    close(
        binomial_cdf_complement(10, 5, 0.5).unwrap(),
        638.0 / 1024.0,
        1e-13,
    );
    // P(X >= 1) = 1 - (1-p)^N, and P(X >= N) = p^N.
    for (trials, p) in [(1_u32, 0.1), (3, 0.1), (10, 0.7), (7, 0.25)] {
        close(
            binomial_cdf_complement(trials, 1, p).unwrap(),
            1.0 - (1.0 - p).powi(trials as i32),
            1e-12,
        );
        close(
            binomial_cdf_complement(trials, trials, p).unwrap(),
            p.powi(trials as i32),
            1e-12,
        );
    }
    // A probability, and monotone decreasing in the success threshold.
    let mut previous = 1.0;
    for successes in 0..=20 {
        let value = binomial_cdf_complement(20, successes, 0.4).unwrap();
        assert!((0.0..=1.0).contains(&value));
        assert!(value <= previous);
        previous = value;
    }
}

// No START_SECTION covers extendRange, contains, zoomIn, createBins, gcd,
// quantile or ppmToMass upstream. The following tests are independently
// derived from the header's documented behaviour and its inline bodies.

#[test]
fn extend_range_and_contains_follow_the_source_comparison_order() {
    let (mut min, mut max) = (1.0, 2.0);
    assert!(extend_range(&mut min, &mut max, 0.5));
    assert_eq!((min, max), (0.5, 2.0));
    assert!(extend_range(&mut min, &mut max, 3.0));
    assert_eq!((min, max), (0.5, 3.0));
    assert!(!extend_range(&mut min, &mut max, 1.0));
    assert_eq!((min, max), (0.5, 3.0));
    // A NaN extends nothing, because both comparisons are false.
    assert!(!extend_range(&mut min, &mut max, f64::NAN));
    assert_eq!((min, max), (0.5, 3.0));
    // The source returns after the first satisfied branch, so one call on an
    // inverted interval lowers min only.
    let (mut lo, mut hi) = (10.0, 0.0);
    assert!(extend_range(&mut lo, &mut hi, 5.0));
    assert_eq!((lo, hi), (5.0, 0.0));

    assert!(contains(1.0, 1.0, 2.0));
    assert!(contains(2.0, 1.0, 2.0));
    assert!(contains(1.5, 1.0, 2.0));
    assert!(!contains(0.9, 1.0, 2.0));
    assert!(!contains(2.1, 1.0, 2.0));
    assert!(!contains(f64::NAN, 1.0, 2.0));
    assert!(!contains(1.5, 2.0, 1.0));
}

#[test]
fn zoom_in_round_trips_and_keeps_the_single_precision_subtraction() {
    // The header's own @code example: zooming in by a factor and back out by
    // its inverse restores the interval.
    let (a2, b2) = zoom_in(10.0, 20.0, 0.5, 0.5).unwrap();
    assert_eq!((a2, b2), (12.5, 17.5));
    assert_eq!(zoom_in(a2, b2, 2.0, 0.5).unwrap(), (10.0, 20.0));
    // align 0 keeps the left edge, align 1 keeps the right edge.
    assert_eq!(zoom_in(10.0, 20.0, 0.5, 0.0).unwrap(), (10.0, 15.0));
    assert_eq!(zoom_in(10.0, 20.0, 0.5, 1.0).unwrap(), (15.0, 20.0));
    // A factor above one zooms out.
    assert_eq!(zoom_in(10.0, 20.0, 2.0, 0.5).unwrap(), (5.0, 25.0));
    // (1.0f - factor) is a float subtraction in the source. At this scale the
    // difference from the f64 subtraction is over twenty units, so the choice
    // is observable rather than academic.
    let (left, _) = zoom_in(0.0, 1e9, 0.1, 1.0).unwrap();
    assert_eq!(left, f64::from(1.0_f32 - 0.1_f32) * 1e9);
    assert!((left - (1.0 - f64::from(0.1_f32)) * 1e9).abs() > 20.0);
    // The source states these as debug-only preconditions.
    assert!(zoom_in(0.0, 1.0, -0.5, 0.5).is_err());
    assert!(zoom_in(0.0, 1.0, 0.5, 1.5).is_err());
    assert!(zoom_in(0.0, 1.0, 0.5, -0.1).is_err());
    assert!(zoom_in(0.0, 1.0, f32::NAN, 0.5).is_err());
    assert!(zoom_in(f64::NAN, 1.0, 0.5, 0.5).is_err());
}

#[test]
fn create_bins_partitions_overlaps_and_never_extends_the_outer_borders() {
    let bins = create_bins(0.0, 10.0, 4, 0.0).unwrap();
    assert_eq!(bins.len(), 4);
    assert_eq!(bins[0], Bin { min: 0.0, max: 2.5 });
    assert_eq!(
        bins[3],
        Bin {
            min: 7.5,
            max: 10.0
        }
    );
    // Adjacent bins touch exactly, and the union is the original interval.
    for pair in bins.windows(2) {
        assert_eq!(pair[0].max, pair[1].min);
    }

    // A margin overlaps neighbours by 2 * margin but leaves the outer borders.
    let overlapping = create_bins(0.0, 10.0, 2, 1.0).unwrap();
    assert_eq!(overlapping[0], Bin { min: 0.0, max: 6.0 });
    assert_eq!(
        overlapping[1],
        Bin {
            min: 4.0,
            max: 10.0
        }
    );
    assert_eq!(overlapping[0].max - overlapping[1].min, 2.0);

    // A negative margin shrinks, and may leave an interior bin empty; that is
    // the source's documented "feature", not an error.
    let shrunk = create_bins(0.0, 3.0, 3, -0.6).unwrap();
    assert!(shrunk[1].is_empty());
    assert_eq!(shrunk[1], Bin { min: 1.6, max: 1.4 });
    assert!(!shrunk[0].is_empty());
    assert_eq!(shrunk[0].min, 0.0);
    assert_eq!(shrunk[2].max, 3.0);

    // A single bin is the whole interval.
    assert_eq!(
        create_bins(-1.0, 1.0, 1, 0.0).unwrap(),
        vec![Bin {
            min: -1.0,
            max: 1.0
        }]
    );
    // Bin::contains mirrors RangeBase::contains.
    let bin = Bin { min: 1.0, max: 2.0 };
    assert!(bin.contains(1.0) && bin.contains(2.0) && bin.contains(1.5));
    assert!(!bin.contains(2.5) && !bin.contains(f64::NAN));
    assert!(!bin.is_empty());

    // Preconditions the source only checks in debug builds.
    assert!(matches!(
        create_bins(1.0, 1.0, 2, 0.0),
        Err(Error::InvalidRange(_))
    ));
    assert!(matches!(
        create_bins(2.0, 1.0, 2, 0.0),
        Err(Error::InvalidRange(_))
    ));
    assert!(create_bins(0.0, 1.0, 0, 0.0).is_err());
    assert!(create_bins(0.0, 1.0, MAX_BINS + 1, 0.0).is_err());
    assert!(create_bins(0.0, f64::INFINITY, 2, 0.0).is_err());
}

#[test]
fn gcd_and_extended_gcd_follow_knuth_and_check_their_arithmetic() {
    assert_eq!(gcd(240, 46).unwrap(), 2);
    assert_eq!(gcd(46, 240).unwrap(), 2);
    assert_eq!(gcd(17, 5).unwrap(), 1);
    assert_eq!(gcd(12, 0).unwrap(), 12);
    assert_eq!(gcd(0, 12).unwrap(), 12);
    // The source does not normalise the sign; it is whatever `a` holds at exit.
    assert_eq!(gcd(-4, 2).unwrap(), 2);
    assert_eq!(gcd(-4, 0).unwrap(), -4);

    // Independently derived: the Bezout identity holds for every pair, and the
    // returned divisor agrees with the plain Euclidean algorithm.
    for (a, b) in [(240_i64, 46_i64), (17, 5), (1071, 462), (12, 0), (0, 12)] {
        let result = extended_gcd(a, b).unwrap();
        assert_eq!(result.gcd, gcd(a, b).unwrap());
        assert_eq!(a * result.u1 + b * result.u2, result.gcd);
    }
    let result = extended_gcd(240, 46).unwrap();
    assert_eq!((result.gcd, result.u1, result.u2), (2, -9, 47));

    // i64::MIN % -1 overflows; the C++ template has undefined behaviour here.
    assert!(gcd(i64::MIN, -1).is_err());
    assert!(extended_gcd(i64::MIN, -1).is_err());
}

#[test]
fn ppm_to_mass_carries_the_sign_and_inverts_the_ppm_computation() {
    assert_eq!(ppm_to_mass(1000.0, 1000.0).unwrap(), 1.0);
    assert_eq!(ppm_to_mass(-1000.0, 1000.0).unwrap(), -1.0);
    assert_eq!(ppm_to_mass_abs(-1000.0, 1000.0).unwrap(), 1.0);
    assert_eq!(ppm_to_mass_abs(1000.0, 1000.0).unwrap(), 1.0);
    // Independently derived: adding the converted mass to the reference
    // reproduces the ppm, because both use the reference as the anchor.
    for (tolerance, reference) in [(10.0, 500.0), (-10.0, 500.0), (5.0, 1200.5)] {
        let delta = ppm_to_mass(tolerance, reference).unwrap();
        close(ppm(reference + delta, reference).unwrap(), tolerance, 1e-9);
    }
    // No divisor, so a zero reference is a zero mass difference, not an error.
    assert_eq!(ppm_to_mass(10.0, 0.0).unwrap(), 0.0);
    assert!(ppm_to_mass(f64::INFINITY, 1.0).is_err());
    assert!(ppm_to_mass_abs(1.0, f64::NAN).is_err());
}

#[test]
fn quantile_uses_the_source_index_convention_and_checks_its_precondition() {
    let sample = [1.0, 2.0, 3.0, 4.0, 5.0];
    // The source index is max(0, n*q - 1), not (n-1)*q, so the "median" of five
    // ascending values is 2.5 rather than 3.0. Reproduced deliberately.
    assert_eq!(quantile(&sample, 0.5).unwrap(), 2.5);
    assert_eq!(quantile(&sample, 0.0).unwrap(), 1.0);
    assert_eq!(quantile(&sample, 1.0).unwrap(), 5.0);
    assert_eq!(quantile(&sample, 0.6).unwrap(), 3.0);
    // q outside [0, 1] is clamped, as in the source.
    assert_eq!(
        quantile(&sample, 5.0).unwrap(),
        quantile(&sample, 1.0).unwrap()
    );
    assert_eq!(
        quantile(&sample, -5.0).unwrap(),
        quantile(&sample, 0.0).unwrap()
    );
    // A single sample is its own every quantile.
    assert_eq!(quantile(&[7.0], 0.0).unwrap(), 7.0);
    assert_eq!(quantile(&[7.0], 1.0).unwrap(), 7.0);
    // Source throws Exception::InvalidParameter on an empty container.
    assert!(matches!(quantile(&[], 0.5), Err(Error::InvalidValue(_))));
    // The source states sortedness in its @brief and never checks it.
    assert!(matches!(
        quantile(&[3.0, 1.0, 2.0], 0.5),
        Err(Error::UnsortedData)
    ));
    assert!(quantile(&[1.0, f64::NAN], 0.5).is_err());
    assert!(quantile(&sample, f64::NAN).is_err());
}
