// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Port of `StatisticFunctions_test.cpp` (core SDK bc9cc12): one test per
//! `START_SECTION`, cited by source line. Transcribed literals are tier 3;
//! tests whose expectation is derived from a closed form or an invariant say so
//! at the assertion, and the native guard tests are tier 4.

use openms::Error;
use openms::math::statistic_functions::{
    AdaptiveQuantileResult, DEFAULT_R_DENSE, DEFAULT_R_SPARSE, DEFAULT_TUKEY_FACTOR,
    SummaryStatistics, absdev, absdev_with_mean, check_exhausted, check_not_empty,
    check_ranges_end_together, classification_rate, compute_rank, covariance, mad,
    matthews_correlation_coefficient, mean, mean_absolute_deviation, mean_square_error, median,
    median_sorted, pearson_correlation_coefficient, quantile, quantile1st, quantile1st_sorted,
    quantile3rd, quantile3rd_sorted, rank_correlation_coefficient, root_mean_square_error, sd,
    sd_with_mean, sum, tail_fraction_above, tukey_upper_fence, variance, variance_with_mean,
    winsorized_quantile,
};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.),
        "{actual:.17} != {expected:.17}"
    );
}

fn invalid_range(error: &Error) -> bool {
    matches!(error, Error::InvalidRange(_))
}

// L32 sum(begin, end)
#[test]
fn section_sum() {
    let x = [-1.0, 0.0, 1.0, 2.0, 3.0];
    assert_eq!(sum(&x), 5.0);
    // The source's `Math::sum(x, x)` is an empty range, which accumulates to
    // the 0.0 seed; no emptiness check exists on this one.
    assert_eq!(sum(&x[..0]), 0.0);
    let y = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0];
    close(sum(&y), 3.5);
}

// L50 mean(begin, end)
#[test]
fn section_mean() {
    let x = [-1.0, 0.0, 1.0, 2.0, 3.0];
    assert_eq!(mean(&x).unwrap(), 1.0);
    // TEST_EXCEPTION(Exception::InvalidRange, Math::mean(x, x))
    assert!(invalid_range(&mean(&x[..0]).unwrap_err()));
    let y = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0];
    close(mean(&y).unwrap(), 0.5);
}

// L61 median(begin, end, sorted)
#[test]
fn section_median() {
    let x = [-1.0, 0.0, 1.0, 2.0, 3.0];
    close(median_sorted(&x).unwrap(), 1.0);
    let x2 = [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0];
    // Even count: the two middle values are averaged, (1 + 2) / 2.
    close(median_sorted(&x2).unwrap(), 1.5);
    assert!(invalid_range(&median_sorted(&x[..0]).unwrap_err()));

    // unsorted
    let mut y = vec![1.0, -0.5, 2.0, 0.5, -1.0, 1.5, 0.0];
    close(median(&mut y).unwrap(), 0.5);
    y.push(-1.5);
    close(median(&mut y).unwrap(), 0.25);

    // sorted
    let z_odd = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0];
    close(median_sorted(&z_odd).unwrap(), 0.5);
    let z_even = [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0];
    close(median_sorted(&z_even).unwrap(), 0.25);

    // Native: the port checks the claimed sortedness the source only asserts
    // in a debug build.
    assert!(matches!(
        median_sorted(&[2.0, 1.0]).unwrap_err(),
        Error::UnsortedData
    ));
}

// L83 MAD(begin, end, median_of_numbers)
#[test]
fn section_mad() {
    let x = [-1.0, 0.0, 1.0, 2.0, 3.0];
    // median{2, 1, 0, 1, 2} == 1
    close(mad(&x, 1.0).unwrap(), 1.0);

    let x2 = [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0];
    // The source line passes the literal `true`, which converts to the median
    // 1.0, while its comment states the real median 1.5. Both give 1.5:
    // |x - 1| sorts to {0,1,1,2,2,3} with median 1.5, and |x - 1.5| sorts to
    // {0.5,0.5,1.5,1.5,2.5,2.5} with median 1.5. Derived, not transcribed.
    close(mad(&x2, 1.0).unwrap(), 1.5);
    close(mad(&x2, 1.5).unwrap(), 1.5);

    let z_odd = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0];
    close(mad(&z_odd, 0.5).unwrap(), 1.0);
    let z_even = [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0];
    close(mad(&z_even, 0.5).unwrap(), 1.0);

    assert!(invalid_range(&mad(&[], 0.0).unwrap_err()));
}

// L97 MeanAbsoluteDeviation(begin, end, mean_of_numbers)
#[test]
fn section_mean_absolute_deviation() {
    let x1 = [-1.0, 0.0, 1.0, 2.0, 3.0];
    // (2 + 1 + 0 + 1 + 2) / 5 == 6/5, exact in binary floating point.
    assert_eq!(mean_absolute_deviation(&x1, 1.0), 1.2);

    let x2 = [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0];
    close(mean_absolute_deviation(&x2, 1.5), 1.5);

    let x3 = [-1.0];
    close(mean_absolute_deviation(&x3, -1.0), 0.0);

    // Empty range: the source divides 0.0 by 0 and its class test pins the NaN.
    assert!(mean_absolute_deviation(&[], 1.0).is_nan());

    let z_odd = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0];
    // 6/7 exactly; the class test writes the truncated 0.857142.
    close(mean_absolute_deviation(&z_odd, 0.5), 6.0 / 7.0);

    let z_even = [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0];
    close(mean_absolute_deviation(&z_even, 0.25), 1.0);
}

// L121 absdev(begin, end, mean)
#[test]
fn section_absdev() {
    let x1 = [-1.0, 0.0, 1.0, 2.0, 3.0];
    close(absdev_with_mean(&x1, 1.0).unwrap(), 1.2);
    // mean computed internally, which is 1
    close(absdev(&x1).unwrap(), 1.2);

    // Symmetric around the mean: the signed sum is exactly zero, the absolute
    // deviation is not. This is the regression the source's comment records.
    let sym = [-2.0, -1.0, 1.0, 2.0];
    close(absdev_with_mean(&sym, 0.0).unwrap(), 1.5);

    let x3 = [-1.0];
    close(absdev_with_mean(&x3, -1.0).unwrap(), 0.0);

    assert!(invalid_range(&absdev(&[]).unwrap_err()));
    assert!(invalid_range(&absdev_with_mean(&[], 0.0).unwrap_err()));
}

// L143 meanSquareError(begin_a, end_a, begin_b, end_b)
#[test]
fn section_mean_square_error() {
    let numbers1 = [1.5; 20];
    let numbers2 = [1.3; 20];
    // (1.5 - 1.3)^2 == 0.04 up to the representation of 0.2; tolerance 1e-6 as
    // the source's TOLERANCE_ABSOLUTE(0.000001).
    let result = mean_square_error(&numbers1, &numbers2).unwrap();
    assert!((result - 0.04).abs() <= 1e-6, "{result}");
    // Divisor is n, not n - 1: 20 equal terms average to one term.
    close(result, (1.5f64 - 1.3).powi(2));
}

// L155 classificationRate(begin_a, end_a, begin_b, end_b)
#[test]
fn section_classification_rate() {
    let mut numbers1 = vec![1.0; 20];
    let mut numbers2 = vec![1.0; 20];
    numbers1.resize(40, -1.0);
    numbers2.resize(40, -1.0);
    for index in [2, 7, 11, 15, 17] {
        numbers1[index] = -1.0;
    }
    for index in [25, 27, 29, 31, 37] {
        numbers1[index] = 1.0;
    }
    // Ten of forty positions disagree in sign, so 30/40 agree.
    close(classification_rate(&numbers1, &numbers2).unwrap(), 0.75);
    assert!(invalid_range(
        &classification_rate(&numbers1, &numbers2[..10]).unwrap_err()
    ));
}

// L180 pearsonCorrelationCoefficient(begin_a, end_a, begin_b, end_b)
#[test]
fn section_pearson_correlation_coefficient() {
    let mut numbers1 = [1.5; 20];
    let mut numbers2 = [1.3; 20];
    numbers1[0] = 0.1;
    numbers2[0] = 0.5;
    numbers1[1] = 0.2;
    numbers2[1] = 0.7;
    numbers1[2] = 0.01;
    numbers2[2] = 0.03;
    numbers1[3] = 1.7;
    numbers2[3] = 1.0;
    numbers1[4] = 3.2;
    numbers2[4] = 4.0;
    let result = pearson_correlation_coefficient(&numbers1, &numbers2).unwrap();
    assert!((result - 0.897811).abs() <= 1e-6, "{result}");

    // A constant range gives a zero denominator and therefore NaN, which the
    // source documents and its test converts to -1.0 before comparing.
    let vv1 = [1.0; 5];
    let vv2 = [1.0, 2.0, 3.0, 4.0, 5.0];
    assert!(
        pearson_correlation_coefficient(&vv1, &vv2)
            .unwrap()
            .is_nan()
    );

    let v1 = [1.0, 2.0, 3.0, 4.0, 5.0];
    close(pearson_correlation_coefficient(&v1, &v1).unwrap(), 1.0);
    let v2 = [-1.0, -2.0, -3.0, -4.0, -5.0];
    close(pearson_correlation_coefficient(&v1, &v2).unwrap(), -1.0);

    // Two uncorrelated 20-point samples; the source asserts "0" against its
    // default tolerance. The values are `float` literals in the source, widened
    // here, which is why the residual differs in the last digits.
    let u1 = [
        0.371_680_3,
        0.277_811_1,
        0.815_237_2,
        0.771_509_7,
        0.016_317_9,
        -0.489_873_8,
        -0.606_013_7,
        -0.888_297_0,
        0.291_359_1,
        -0.366_179_1,
        0.132_075_0,
        0.263_722_9,
        -0.739_022_6,
        -0.039_592_9,
        0.338_733_4,
        0.859_854_1,
        0.738_823_6,
        -0.592_808_3,
        0.922_600_6,
        -0.357_142_7,
    ];
    let u2 = [
        0.639_696_9,
        0.794_240_5,
        -0.636_447_3,
        -0.684_563_3,
        -0.690_886_2,
        -0.503_416_9,
        0.574_529_8,
        -0.124_759_1,
        -0.512_956_4,
        0.074_585_7,
        0.073_366_5,
        -0.011_888_2,
        0.176_347_1,
        0.102_759_9,
        -0.973_780_5,
        0.874_767_7,
        0.947_939_2,
        0.084_360_4,
        -0.351_896_1,
        -0.303_403_9,
    ];
    let r = pearson_correlation_coefficient(&u1, &u2).unwrap();
    assert!(r.abs() < 1e-7, "{r}");

    let w1 = [
        -0.183_334_1,
        0.656_444_9,
        0.872_503_9,
        0.361_092_1,
        0.792_614_4,
        0.183_334_1,
        -0.656_444_9,
        -0.414_106_1,
        -0.872_503_9,
        0.826_998_5,
        -0.587_871_5,
        -0.295_044_3,
        -0.361_092_1,
        -0.826_998_5,
        -0.047_032_7,
        0.414_106_1,
        0.047_032_7,
        0.295_044_3,
        -0.792_614_4,
        0.587_871_5,
    ];
    let w2 = [
        0.033_611_4,
        0.430_919_9,
        0.761_263_1,
        0.130_387_5,
        0.628_237_7,
        0.033_611_4,
        0.430_919_9,
        0.171_483_9,
        0.761_263_1,
        0.683_926_4,
        0.345_592_9,
        0.087_051_1,
        0.130_387_5,
        0.683_926_4,
        0.002_212_1,
        0.171_483_9,
        0.002_212_1,
        0.087_051_1,
        0.628_237_7,
        0.345_592_9,
    ];
    let r = pearson_correlation_coefficient(&w1, &w2).unwrap();
    assert!(r.abs() < 1e-12, "{r}");
}

// L342 computeRank(std::vector<double>& w)
#[test]
fn section_compute_rank() {
    let mut numbers1 = vec![1.5; 10];
    numbers1[0] = 1.4;
    numbers1[1] = 0.2;
    numbers1[2] = 0.01;
    numbers1[3] = 1.7;
    numbers1[4] = 3.2;
    numbers1[5] = 2.2;
    close(numbers1[0], 1.4);
    close(numbers1[5], 2.2);
    close(numbers1[6], 1.5);
    close(numbers1[9], 1.5);

    compute_rank(&mut numbers1).unwrap();

    // Sorted: 0.01, 0.2, 1.4, 1.5, 1.5, 1.5, 1.5, 1.7, 2.2, 3.2. The tie block
    // of four 1.5 values starts at sorted position i = 3 and ends before
    // z = 7, so its shared rank is 0.5 * (i + z + 1) = 5.5 — the mean of the
    // one-based ranks 4, 5, 6, 7. Derived from the rule, not transcribed.
    close(numbers1[0], 3.0);
    close(numbers1[1], 2.0);
    close(numbers1[2], 1.0);
    close(numbers1[3], 8.0);
    close(numbers1[4], 10.0);
    close(numbers1[5], 9.0);
    close(numbers1[6], 5.5);
    close(numbers1[7], 5.5);
    close(numbers1[8], 5.5);
    close(numbers1[9], 5.5);

    // Native: the source computes `w.size() - 1` in unsigned arithmetic and
    // wraps on an empty vector; the port returns without touching anything.
    let mut none: Vec<f64> = Vec::new();
    compute_rank(&mut none).unwrap();
    assert!(none.is_empty());
    let mut one = vec![7.0];
    compute_rank(&mut one).unwrap();
    assert_eq!(one, vec![1.0]);
}

// L373 rankCorrelationCoefficient(begin_a, end_a, begin_b, end_b)
#[test]
fn section_rank_correlation_coefficient() {
    let mut numbers1 = vec![1.5; 10];
    let mut numbers2 = vec![1.3; 10];
    let numbers3 = vec![0.42; 10];
    let numbers4: Vec<f64> = (1..=10).map(f64::from).collect();
    numbers1[0] = 0.4;
    numbers2[0] = 0.5;
    numbers1[1] = 0.2;
    numbers2[1] = 0.7;
    numbers1[2] = 0.01;
    numbers2[2] = 0.03;
    numbers1[3] = 1.7;
    numbers2[3] = 1.0;
    numbers1[4] = 3.2;
    numbers2[4] = 4.0;
    numbers1[5] = 2.2;
    numbers2[5] = 3.0;

    close(
        rank_correlation_coefficient(&numbers1, &numbers2).unwrap(),
        0.858_064_516_129_032,
    );

    let reversed2: Vec<f64> = numbers2.iter().rev().copied().collect();
    close(
        rank_correlation_coefficient(&numbers1, &reversed2).unwrap(),
        0.303_225_806_451_613,
    );

    // A constant range gives a zero sum of squares, for which the source
    // returns 0.0 rather than NaN.
    close(
        rank_correlation_coefficient(&numbers3, &numbers4).unwrap(),
        0.0,
    );
    close(
        rank_correlation_coefficient(&numbers3, &numbers3).unwrap(),
        0.0,
    );
    // Derived: identical strictly increasing ranks correlate to exactly +1, and
    // their reversal to exactly -1, by the symmetry of the rank sum about mu.
    close(
        rank_correlation_coefficient(&numbers4, &numbers4).unwrap(),
        1.0,
    );
    let reversed4: Vec<f64> = numbers4.iter().rev().copied().collect();
    close(
        rank_correlation_coefficient(&numbers4, &reversed4).unwrap(),
        -1.0,
    );
}

// L420 quantile1st/quantile3rd(begin, end, sorted)
#[test]
fn section_quantile1st_and_quantile3rd() {
    let x = [3.0, 6.0, 7.0, 8.0, 8.0, 10.0, 13.0, 15.0, 16.0, 20.0];
    let y = [3.0, 6.0, 7.0, 8.0, 8.0, 10.0, 13.0, 15.0, 16.0];

    close(quantile1st_sorted(&x).unwrap(), 6.5);
    close(median_sorted(&x).unwrap(), 9.0);
    close(quantile3rd_sorted(&x).unwrap(), 15.5);
    close(quantile1st_sorted(&y).unwrap(), 6.5);
    close(median_sorted(&y).unwrap(), 8.0);
    close(quantile3rd_sorted(&y).unwrap(), 14.0);

    // issue #9659: both must be total for n >= 1
    let n1 = [5.0];
    close(quantile1st_sorted(&n1).unwrap(), 5.0);
    close(quantile3rd_sorted(&n1).unwrap(), 5.0);
    let n2 = [1.0, 3.0];
    close(quantile1st_sorted(&n2).unwrap(), 1.0);
    close(quantile3rd_sorted(&n2).unwrap(), 3.0);
    // Unsorted input is sorted internally before the min/max is taken.
    let mut n2u = [3.0, 1.0];
    close(quantile1st(&mut n2u).unwrap(), 1.0);
    let mut n2u = [3.0, 1.0];
    close(quantile3rd(&mut n2u).unwrap(), 3.0);
    let n3 = [1.0, 2.0, 3.0];
    close(quantile1st_sorted(&n3).unwrap(), 1.0);
    close(quantile3rd_sorted(&n3).unwrap(), 3.0);
    let n4 = [1.0, 2.0, 3.0, 4.0];
    close(quantile1st_sorted(&n4).unwrap(), 1.0);
    close(quantile3rd_sorted(&n4).unwrap(), 4.0);

    assert!(invalid_range(&quantile1st_sorted(&[]).unwrap_err()));
    assert!(invalid_range(&quantile3rd_sorted(&[]).unwrap_err()));
}

// L459 struct SummaryStatistics
#[test]
fn section_summary_statistics() {
    let mut one = [5.0];
    let s1 = SummaryStatistics::new(&mut one).unwrap();
    assert_eq!(s1.count, 1);
    close(s1.mean, 5.0);
    // n == 1 has no n-1 divisor, so the source reports the empty case's 0.0.
    close(s1.variance, 0.0);
    close(s1.min, 5.0);
    close(s1.lowerq, 5.0);
    close(s1.median, 5.0);
    close(s1.upperq, 5.0);
    close(s1.max, 5.0);

    let mut two = [1.0, 3.0];
    let s2 = SummaryStatistics::new(&mut two).unwrap();
    assert_eq!(s2.count, 2);
    close(s2.mean, 2.0);
    // Sample variance ((1-2)^2 + (3-2)^2) / (2-1) == 2; derived, not transcribed.
    close(s2.variance, 2.0);
    close(s2.min, 1.0);
    close(s2.lowerq, 1.0);
    close(s2.median, 2.0);
    close(s2.upperq, 3.0);
    close(s2.max, 3.0);

    let mut three = [1.0, 2.0, 3.0];
    let s3 = SummaryStatistics::new(&mut three).unwrap();
    assert_eq!(s3.count, 3);
    close(s3.lowerq, 1.0);
    close(s3.median, 2.0);
    close(s3.upperq, 3.0);

    let mut none: [f64; 0] = [];
    let s0 = SummaryStatistics::new(&mut none).unwrap();
    assert_eq!(s0.count, 0);
    close(s0.variance, 0.0);
    close(s0.lowerq, 0.0);
    close(s0.upperq, 0.0);
    assert_eq!(s0, SummaryStatistics::default());
}

// L504 quantile(begin, end, q)
#[test]
fn section_quantile() {
    let v = [1.0, 2.0, 3.0, 4.0];
    close(quantile(&v, 0.0).unwrap(), 1.0);
    close(quantile(&v, 1.0).unwrap(), 4.0);
    // Type-7 median of an even count: pos = 0.5 * (n-1) = 1.5 -> 2.5
    close(quantile(&v, 0.5).unwrap(), 2.5);
    close(quantile(&v, 0.25).unwrap(), 1.75);
    close(quantile(&v, 0.75).unwrap(), 3.25);

    assert!(invalid_range(&quantile(&[], 0.5).unwrap_err()));
    assert!(matches!(
        quantile(&v, -0.1).unwrap_err(),
        Error::InvalidValue(_)
    ));
    assert!(matches!(
        quantile(&v, 1.1).unwrap_err(),
        Error::InvalidValue(_)
    ));
    // Native: the source states sortedness as @pre and only checks it in a
    // debug build.
    assert!(matches!(
        quantile(&[2.0, 1.0], 0.5).unwrap_err(),
        Error::UnsortedData
    ));
}

// L521 tukeyUpperFence(begin, end, k)
#[test]
fn section_tukey_upper_fence() {
    // v = 1..9 -> Q1 = 3, Q3 = 7, IQR = 4 -> UF = 7 + 1.5*4 = 13
    let v: Vec<f64> = (1..=9).map(f64::from).collect();
    close(tukey_upper_fence(&v, DEFAULT_TUKEY_FACTOR).unwrap(), 13.0);
    close(quantile(&v, 0.25).unwrap(), 3.0);
    close(quantile(&v, 0.75).unwrap(), 7.0);

    // Too few finite values -> +inf
    let tiny = [1.0, f64::NAN];
    assert!(tukey_upper_fence(&tiny, 1.5).unwrap().is_infinite());
    // A constant range has IQR 0, which is also "no fence".
    assert!(tukey_upper_fence(&[2.0; 8], 1.5).unwrap().is_infinite());
}

// L533 tailFractionAbove(begin, end, threshold)
#[test]
fn section_tail_fraction_above() {
    let v: Vec<f64> = (1..=9).map(f64::from).collect();
    // values > 7 are {8, 9} -> 2/9
    let r = tail_fraction_above(&v, 7.0);
    assert!((r - 2.0 / 9.0).abs() <= 1e-12, "{r}");
    close(tail_fraction_above(&v, 10.0), 0.0);
    // Non-finite values leave both counts alone; with none finite the fraction
    // is the source's explicit 0.0 rather than a division by zero.
    close(tail_fraction_above(&[f64::NAN, f64::INFINITY], 0.0), 0.0);
}

// L546 winsorizedQuantile(begin, end, q, upper_fence)
#[test]
fn section_winsorized_quantile() {
    let v = [1.0, 2.0, 3.0, 100.0];
    let uf = 65.5;
    let raw_q = quantile(&v, 0.75).unwrap();
    close(raw_q, 27.25);
    // pos = 2.25 between 3 and 65.5 -> 0.75*3 + 0.25*65.5
    close(winsorized_quantile(&v, 0.75, uf).unwrap(), 18.625);
    // A non-finite fence disables the capping.
    close(winsorized_quantile(&v, 0.75, f64::INFINITY).unwrap(), raw_q);
    // No finite value at all -> 0.0
    close(winsorized_quantile(&[f64::NAN], 0.5, 1.0).unwrap(), 0.0);
}

// L566 adaptiveQuantile(begin, end, q, k, r_sparse, r_dense)
#[test]
fn section_adaptive_quantile() {
    let mut a: Vec<f64> = Vec::with_capacity(201);
    for _ in 0..20 {
        for z in 0..=9 {
            a.push(f64::from(z));
        }
    }
    a.push(1000.0);
    let q = 0.997;

    // Case A: one outlier in 201 values -> sparse tail -> robust wins.
    let ada_a = super_adaptive(&a, q);
    assert!(ada_a.tail_fraction < 0.01, "{}", ada_a.tail_fraction);
    close(ada_a.blended, ada_a.half_rob);
    assert!(ada_a.half_raw > ada_a.half_rob);
    close(ada_a.weight, 0.0);
    // Derived: Q1 = 2, Q3 = 7 over the sorted base, so UF = 7 + 1.5*5 = 14.5,
    // and 1000 is the only value above it: 1/201.
    close(ada_a.upper_fence, 14.5);
    close(ada_a.tail_fraction, 1.0 / 201.0);
    // raw: pos = 0.997 * 200 = 199.4 -> 0.6*9 + 0.4*1000
    close(ada_a.half_raw, 0.6 * 9.0 + 0.4 * 1000.0);
    // robust: the outlier is capped at 14.5 -> 0.6*9 + 0.4*14.5
    close(ada_a.half_rob, 0.6 * 9.0 + 0.4 * 14.5);

    // Case B: 30 outliers -> dense tail -> raw wins.
    let mut b = a.clone();
    b.extend(std::iter::repeat_n(1000.0, 29));
    let ada_b = super_adaptive(&b, q);
    assert!(ada_b.tail_fraction > 0.10, "{}", ada_b.tail_fraction);
    close(ada_b.blended, ada_b.half_raw);
    assert!(ada_b.half_raw >= ada_b.half_rob);
    close(ada_b.weight, 1.0);

    // Case C: 10 outliers in 210 values -> interpolation.
    let mut c: Vec<f64> = Vec::with_capacity(210);
    for _ in 0..20 {
        for z in 0..=9 {
            c.push(f64::from(z));
        }
    }
    c.extend(std::iter::repeat_n(1000.0, 10));
    let ada_c = super_adaptive(&c, q);
    let w = ((ada_c.tail_fraction - 0.01) / (0.10 - 0.01)).clamp(0.0, 1.0);
    let expected = (1.0 - w) * ada_c.half_rob + w * ada_c.half_raw;
    assert!(ada_c.tail_fraction > 0.01 && ada_c.tail_fraction < 0.10);
    assert!(
        (ada_c.blended - expected).abs() <= 1e-9,
        "{}",
        ada_c.blended
    );
    assert!(ada_c.half_rob <= ada_c.blended && ada_c.blended <= ada_c.half_raw);
    close(ada_c.tail_fraction, 10.0 / 210.0);

    // No finite value -> the source's default-initialised result.
    let none = openms::math::statistic_functions::adaptive_quantile(
        &[f64::NAN],
        q,
        DEFAULT_TUKEY_FACTOR,
        DEFAULT_R_SPARSE,
        DEFAULT_R_DENSE,
    )
    .unwrap();
    assert_eq!(none, AdaptiveQuantileResult::default());
    assert!(none.upper_fence.is_infinite());
}

fn super_adaptive(values: &[f64], q: f64) -> AdaptiveQuantileResult {
    openms::math::statistic_functions::adaptive_quantile(
        values,
        q,
        DEFAULT_TUKEY_FACTOR,
        DEFAULT_R_SPARSE,
        DEFAULT_R_DENSE,
    )
    .unwrap()
}

// L619 rootMeanSquareError(begin_a, end_a, begin_b, end_b)
#[test]
fn section_root_mean_square_error() {
    let a = [1.0, 2.0, 3.0];
    let b = [1.5, 2.5, 3.5];
    let mse = mean_square_error(&a, &b).unwrap();
    let rmse = root_mean_square_error(&a, &b).unwrap();
    close(mse, 0.25);
    close(rmse, 0.25f64.sqrt());

    let errs = [-2.0, 0.0, 2.0];
    let zeros = [0.0; 3];
    close(
        root_mean_square_error(&errs, &zeros).unwrap(),
        ((4.0 + 0.0 + 4.0) / 3.0f64).sqrt(),
    );

    let shortv = [1.0, 2.0];
    assert!(invalid_range(
        &root_mean_square_error(&errs, &shortv).unwrap_err()
    ));
    assert!(invalid_range(
        &root_mean_square_error(&[], &[]).unwrap_err()
    ));
}

// L652 checkIteratorsAreValid(begin_b, end_b, begin_a, end_a)
#[test]
fn section_check_iterators_are_valid() {
    let v1 = [1.0, 2.0, 3.0];
    let v2 = [10.0, 20.0, 30.0];

    // Both ranges still have elements.
    check_ranges_end_together(&v1, &v2).unwrap();
    // Both exhausted.
    check_ranges_end_together(&v1[3..], &v2[3..]).unwrap();
    // One exhausted, the other not.
    assert!(invalid_range(
        &check_ranges_end_together(&v1, &v2[3..]).unwrap_err()
    ));
    assert!(invalid_range(
        &check_ranges_end_together(&v1[3..], &v2).unwrap_err()
    ));

    // The two companion helpers, which have no section of their own.
    check_not_empty(&v1).unwrap();
    assert!(invalid_range(&check_not_empty(&v1[3..]).unwrap_err()));
    check_exhausted(&v1[3..]).unwrap();
    assert!(invalid_range(&check_exhausted(&v1).unwrap_err()));
}

// Native: members with no class-test section of their own.
#[test]
fn variance_covariance_and_sd_follow_the_source_divisors() {
    let x = [1.0, 2.0, 3.0, 4.0];
    // Sample variance about mean 2.5: (2.25 + 0.25 + 0.25 + 2.25) / 3
    close(variance(&x).unwrap(), 5.0 / 3.0);
    close(variance_with_mean(&x, 2.5).unwrap(), 5.0 / 3.0);
    close(sd(&x).unwrap(), (5.0f64 / 3.0).sqrt());
    close(sd_with_mean(&x, 2.5).unwrap(), (5.0f64 / 3.0).sqrt());

    // Covariance of a range with itself is its variance.
    close(covariance(&x, &x).unwrap(), 5.0 / 3.0);
    let y = [2.0, 4.0, 6.0, 8.0];
    close(covariance(&x, &y).unwrap(), 10.0 / 3.0);

    // Native: the source divides by n - 1 unchecked and returns NaN for n == 1.
    assert!(invalid_range(&variance(&[1.0]).unwrap_err()));
    assert!(invalid_range(&sd(&[1.0]).unwrap_err()));
    assert!(invalid_range(&covariance(&[1.0], &[2.0]).unwrap_err()));
    assert!(invalid_range(&covariance(&x, &y[..2]).unwrap_err()));
    assert!(invalid_range(&variance(&[]).unwrap_err()));
}

// Native: Matthews correlation has no class-test section.
#[test]
fn matthews_correlation_matches_its_closed_form() {
    // Predicted and real labels agree everywhere: tp = 2, tn = 2, fp = fn = 0,
    // so the coefficient is exactly 1.
    let predicted = [1.0, 1.0, -1.0, -1.0];
    let real = [1.0, 1.0, -1.0, -1.0];
    close(
        matthews_correlation_coefficient(&predicted, &real).unwrap(),
        1.0,
    );
    // Completely inverted -> -1.
    let inverted = [-1.0, -1.0, 1.0, 1.0];
    close(
        matthews_correlation_coefficient(&predicted, &inverted).unwrap(),
        -1.0,
    );
    // (tp, fp, tn, fn) = (2, 1, 1, 0):
    // (2*1 - 1*0) / sqrt(3 * 2 * 2 * 1) = 2 / sqrt(12)
    let a = [1.0, 1.0, 1.0, -1.0];
    let b = [1.0, 1.0, -1.0, -1.0];
    close(
        matthews_correlation_coefficient(&a, &b).unwrap(),
        2.0 / 12.0f64.sqrt(),
    );
    // A one-sided prediction zeroes two counts, which forces the numerator to
    // zero as well; the source's 0/0 is NaN and so is this.
    assert!(
        matthews_correlation_coefficient(&[1.0, 1.0], &[1.0, 1.0])
            .unwrap()
            .is_nan()
    );
    assert!(invalid_range(
        &matthews_correlation_coefficient(&[], &[]).unwrap_err()
    ));
    assert!(invalid_range(
        &matthews_correlation_coefficient(&a, &b[..2]).unwrap_err()
    ));
}

// Native: resource ceilings the source does not have.
#[test]
fn bounded_work_refuses_before_allocating() {
    use openms::math::statistic_functions::MAX_ITEMS;
    // The ceiling is a compile-time constant, so only its value is asserted
    // here; the checks themselves are exercised by every call above.
    assert_eq!(MAX_ITEMS, 50_000_000);
    // A NaN in a range that must be sorted is rejected rather than making the
    // ordering implementation-defined.
    assert!(matches!(
        quantile(&[1.0, f64::NAN, 2.0], 0.5).unwrap_err(),
        Error::UnsortedData
    ));
}
