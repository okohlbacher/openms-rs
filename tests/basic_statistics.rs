// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Port of `BasicStatistics_test.cpp` (core SDK bc9cc12): one test per
//! `START_SECTION`, cited by source line. The file opens 18 sections; a plain
//! `grep -c START_SECTION` reports 20, because two further occurrences at its
//! lines 341 and 356 sit inside comments that refer to a section rather than
//! opening one. Transcribed literals are tier 3; the closed-form checks say so
//! at the assertion.

use openms::Error;
use openms::math::basic_statistics::BasicStatistics;

/// The 195-value sample the class test declares at L35.
const DVECTOR_DATA: [f64; 195] = [
    82.70033, 18.53697, 130.43985, 71.42455, 50.63099, 20.31581, 30.19521, 36.79161, 135.08596,
    84.68491, 124.30681, 71.33620, 126.07538, 73.61598, 130.07241, 88.97545, 112.80919, 81.12736,
    170.80468, 74.20200, 29.40524, 44.20175, 124.63237, 84.51534, 165.35688, 79.33067, 68.44432,
    18.62523, 112.01351, 77.03597, 29.93905, 49.71414, 30.82335, 61.01894, 113.46661, 78.16001,
    162.25406, 89.78833, 158.70900, 74.51220, 73.57289, 124.63237, 84.51534, 165.35688, 79.33067,
    68.44432, 18.62523, 112.01351, 77.03597, 29.93905, 49.71414, 30.82335, 61.01894, 113.46661,
    78.16001, 162.25406, 89.78833, 158.70900, 74.51220, 73.57289, 17.14514, 130.14515, 83.68410,
    29.89634, 47.08373, 76.58917, 29.00928, 57.22767, 22.04459, 108.34564, 79.49656, 140.83229,
    67.81030, 28.82848, 78.72329, 31.32767, 62.28604, 29.48579, 76.01188, 142.99623, 71.69667,
    140.45532, 78.81924, 57.99051, 19.66125, 29.71268, 63.73135, 65.07940, 27.78494, 127.22279,
    67.27982, 29.50484, 142.99623, 71.69667, 140.45532, 78.81924, 57.99051, 19.66125, 29.71268,
    63.73135, 65.07940, 27.78494, 127.22279, 67.27982, 29.50484, 142.99623, 71.69667, 140.45532,
    78.81924, 57.99051, 19.66125, 29.71268, 63.73135, 65.07940, 27.78494, 127.22279, 67.27982,
    29.50484, 54.54108, 30.53517, 86.44319, 67.76178, 18.95834, 123.73745, 77.66034, 30.29570,
    60.94120, 142.92731, 82.77405, 141.99247, 76.17666, 157.02459, 78.28177, 96.25540, 19.82469,
    27.72561, 53.91157, 29.91151, 60.05424, 61.35466, 16.14011, 163.18400, 77.86948, 153.28102,
    91.43451, 29.32177, 83.93723, 111.66644, 80.25561, 129.31559, 90.71809, 107.97381, 75.83463,
    147.61897, 78.47707, 29.93856, 68.92398, 177.78189, 81.44311, 68.58626, 24.30645, 132.16980,
    79.22136, 28.12488, 78.71920, 151.88722, 83.39256, 29.69833, 71.72692, 52.76207, 15.71214,
    116.18279, 75.74875, 115.52147, 91.14405, 127.02429, 95.27849, 67.42286, 20.34733, 102.67339,
    93.84615, 128.95366, 69.28015, 138.62953, 94.72963, 129.24376, 66.28535, 27.90273, 58.98529,
    29.84631, 47.59564, 118.73823, 77.77458, 72.75859, 18.41622,
];

/// The `float` coordinate array the class test builds at L101 and L195,
/// `1000 - i` rounded to `f32` exactly as the source stores it.
fn coordinates() -> Vec<f64> {
    (0..195).map(|i| f64::from(1000.0f32 - i as f32)).collect()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17} != {expected:.17}"
    );
}

// L86 BasicStatistics()
#[test]
fn section_default_constructor() {
    let mut stats = BasicStatistics::new();
    stats.update(&DVECTOR_DATA).unwrap();

    assert_eq!(DVECTOR_DATA.len(), 195);
    // TOLERANCE_ABSOLUTE(0.1)
    near(stats.sum(), 15228.2, 0.1);
    near(stats.mean(), 96.4639, 0.1);
    near(stats.variance(), 3276.51, 0.1);

    let mut stats2 = BasicStatistics::default();
    stats2
        .update_with_coordinates(&DVECTOR_DATA, &coordinates())
        .unwrap();
    near(stats2.sum(), stats.sum(), 0.1);
    // Coordinates 1000 - i mirror the positions, so the weighted mean mirrors
    // too; derived from the reflection, not transcribed.
    near(stats2.mean(), 1000.0 - stats.mean(), 0.1);
    near(stats2.variance(), 3276.51, 0.1);

    assert_eq!(BasicStatistics::new(), BasicStatistics::default());
}

// L114 BasicStatistics(BasicStatistics const& arg)
#[test]
fn section_copy_constructor() {
    let mut stats = BasicStatistics::new();
    stats.update(&DVECTOR_DATA).unwrap();
    // C++ copy construction is a `Copy` here; the type holds three doubles.
    let stats_copy = stats;
    near(stats_copy.sum(), stats.sum(), 0.1);
    near(stats_copy.mean(), stats.mean(), 0.1);
    near(stats_copy.variance(), stats.variance(), 0.1);
    assert_eq!(stats_copy, stats);
}

// L132 operator=(BasicStatistics const& arg)
#[test]
fn section_assignment() {
    let mut stats = BasicStatistics::new();
    stats.update(&DVECTOR_DATA).unwrap();
    let mut stats_copy = BasicStatistics::new();
    assert_eq!(stats_copy, BasicStatistics::default());
    stats_copy = stats;
    near(stats_copy.sum(), stats.sum(), 0.1);
    near(stats_copy.mean(), stats.mean(), 0.1);
    near(stats_copy.variance(), stats.variance(), 0.1);
}

// L152 void clear()
#[test]
fn section_clear() {
    let mut stats = BasicStatistics::new();
    stats.update(&DVECTOR_DATA).unwrap();
    near(stats.sum(), 15228.2, 0.1);
    near(stats.mean(), 96.4639, 0.1);
    near(stats.variance(), 3276.51, 0.1);
    stats.clear();
    assert_eq!(stats.sum(), 0.0);
    assert_eq!(stats.mean(), 0.0);
    assert_eq!(stats.variance(), 0.0);
}

// L168 update(probability_begin, probability_end)
#[test]
fn section_update_from_probabilities() {
    let mut stats = BasicStatistics::new();
    assert_eq!(stats.sum(), 0.0);
    assert_eq!(stats.mean(), 0.0);
    assert_eq!(stats.variance(), 0.0);

    stats.update(&DVECTOR_DATA).unwrap();
    near(stats.sum(), 15228.2, 0.1);
    near(stats.mean(), 96.4639, 0.1);
    near(stats.variance(), 3276.51, 0.1);

    // An all-zero probability vector has zero mass; the source detects the
    // resulting non-finite mean and substitutes zeros.
    stats.update(&[0.0; 12]).unwrap();
    assert_eq!(stats.sum(), 0.0);
    assert_eq!(stats.mean(), 0.0);
    assert_eq!(stats.variance(), 0.0);

    // Native: an empty input takes the same path.
    stats.update(&DVECTOR_DATA).unwrap();
    stats.update(&[]).unwrap();
    assert_eq!(stats.sum(), 0.0);
    assert_eq!(stats.mean(), 0.0);
    assert_eq!(stats.variance(), 0.0);
}

// L187 update(probability_begin, probability_end, coordinate_begin)
#[test]
fn section_update_with_coordinates() {
    let mut stats = BasicStatistics::new();
    assert_eq!(stats.sum(), 0.0);
    assert_eq!(stats.mean(), 0.0);
    assert_eq!(stats.variance(), 0.0);

    stats
        .update_with_coordinates(&DVECTOR_DATA, &coordinates())
        .unwrap();
    near(stats.sum(), 15228.2, 0.1);
    near(stats.mean(), 1000.0 - 96.4639, 0.1);
    near(stats.variance(), 3276.51, 0.1);

    // Native: the source has no coordinate end iterator and reads past the end
    // of a short coordinate range; the port refuses.
    assert!(matches!(
        stats
            .update_with_coordinates(&DVECTOR_DATA, &coordinates()[..10])
            .unwrap_err(),
        Error::InvalidValue(_)
    ));
    // The refused call left the previous state intact.
    near(stats.sum(), 15228.2, 0.1);

    // This overload has no zero-mass guard in the source, so the moments stay
    // NaN where `update` would report zeros.
    let mut zero_mass = BasicStatistics::new();
    zero_mass
        .update_with_coordinates(&[0.0, 0.0], &[1.0, 2.0])
        .unwrap();
    assert_eq!(zero_mass.sum(), 0.0);
    assert!(zero_mass.mean().is_nan());
    assert!(zero_mass.variance().is_nan());
}

// L209 RealType mean() const
#[test]
fn section_mean_getter() {
    let bid = BasicStatistics::new();
    assert_eq!(bid.mean(), 0.0);
}

// L216 void setMean(RealType const& mean)
#[test]
fn section_set_mean() {
    let mut bid = BasicStatistics::new();
    assert_eq!(bid.mean(), 0.0);
    bid.set_mean(17.0);
    assert_eq!(bid.mean(), 17.0);
    // The setter touches nothing else.
    assert_eq!(bid.variance(), 0.0);
    assert_eq!(bid.sum(), 0.0);
}

// L224 RealType variance() const
#[test]
fn section_variance_getter() {
    let bid = BasicStatistics::new();
    assert_eq!(bid.variance(), 0.0);
}

// L231 void setVariance(RealType const& variance)
#[test]
fn section_set_variance() {
    let mut bid = BasicStatistics::new();
    assert_eq!(bid.variance(), 0.0);
    bid.set_variance(18.0);
    assert_eq!(bid.variance(), 18.0);
    assert_eq!(bid.mean(), 0.0);
    assert_eq!(bid.sum(), 0.0);
}

// L239 RealType sum() const
#[test]
fn section_sum_getter() {
    let bid = BasicStatistics::new();
    assert_eq!(bid.sum(), 0.0);
}

// L246 void setSum(RealType const& sum)
#[test]
fn section_set_sum() {
    let mut bid = BasicStatistics::new();
    assert_eq!(bid.sum(), 0.0);
    bid.set_sum(19.0);
    assert_eq!(bid.sum(), 19.0);
    assert_eq!(bid.mean(), 0.0);
    assert_eq!(bid.variance(), 0.0);
}

// L256 static RealType sqrt2pi()
#[test]
fn section_sqrt2pi() {
    // The source asserts exact equality with the literal
    // 2.50662827463100050240, which parses to this f64.
    assert_eq!(BasicStatistics::sqrt2pi(), 2.506_628_274_631_000_7);
    assert_eq!(BasicStatistics::SQRT_2PI, BasicStatistics::sqrt2pi());
    // Derived: that literal is one ulp above the correctly rounded
    // sqrt(2 pi) = 2.5066282746310002, and the port keeps the source's value.
    let exact = (2.0 * std::f64::consts::PI).sqrt();
    assert_ne!(BasicStatistics::SQRT_2PI, exact);
    assert_eq!(
        BasicStatistics::SQRT_2PI.to_bits(),
        exact.to_bits() + 1,
        "the source constant must stay exactly one ulp above sqrt(2 pi)"
    );
}

// L262 RealType normalDensity_sqrt2pi(RealType coordinate) const
#[test]
fn section_normal_density_sqrt2pi() {
    let mut bid = BasicStatistics::new();
    bid.set_mean(10.0);
    bid.set_variance(3.0);
    // TOLERANCE_ABSOLUTE(.0001)
    near(bid.normal_density_sqrt2pi(10.0), 1.0, 1e-4);
    near(
        bid.normal_density_sqrt2pi(7.0),
        0.223_130_160_148_429_82,
        1e-4,
    );
    near(bid.normal_density_sqrt2pi(9.0), 0.846_481_724_890_614, 1e-4);
    near(
        bid.normal_density_sqrt2pi(11.0),
        0.846_481_724_890_614,
        1e-4,
    );
    // Derived: exp(-(c - mu)^2 / 2 / sigma^2) is symmetric about the mean.
    assert_eq!(
        bid.normal_density_sqrt2pi(9.0),
        bid.normal_density_sqrt2pi(11.0)
    );

    // Native: a zero variance is the source's unchecked division, whose value
    // is 0 away from the mean and NaN at it.
    let mut degenerate = BasicStatistics::new();
    degenerate.set_mean(1.0);
    assert_eq!(degenerate.normal_density_sqrt2pi(2.0), 0.0);
    assert!(degenerate.normal_density_sqrt2pi(1.0).is_nan());
}

// L275 RealType normalDensity(RealType const coordinate) const
#[test]
fn section_normal_density() {
    let mut bid = BasicStatistics::new();
    bid.set_mean(10.0);
    bid.set_variance(3.0);
    let sqrt2pi = BasicStatistics::SQRT_2PI;
    near(bid.normal_density(10.0), 1.0 / sqrt2pi, 1e-4);
    near(
        bid.normal_density(7.0),
        0.223_130_160_148_429_82 / sqrt2pi,
        1e-4,
    );
    near(
        bid.normal_density(9.0),
        0.846_481_724_890_614 / sqrt2pi,
        1e-4,
    );
    near(
        bid.normal_density(11.0),
        0.846_481_724_890_614 / sqrt2pi,
        1e-4,
    );
}

// L289 normalApproximation(probability, size)
#[test]
fn section_normal_approximation_with_size() {
    let dvector2_data = [0.0, 1.0, 3.0, 2.0, 0.0];
    let mut stats = BasicStatistics::new();
    stats.update(&dvector2_data).unwrap();

    assert_eq!(dvector2_data.len(), 5);
    near(stats.sum(), 6.0, 0.1);
    // Derived: mean = (1*1 + 3*2 + 2*3) / 6 = 13/6 exactly.
    near(stats.mean(), 13.0 / 6.0, 0.1);
    near(stats.variance(), 17.0 / 36.0, 0.1);
    near(BasicStatistics::sqrt2pi(), 2.506_628_274_631_000_7, 1e-6);

    let probs = stats.normal_approximation(6).unwrap();
    let good_probs = [
        0.0241689,
        0.824253,
        3.38207,
        1.66963,
        0.0991695,
        0.000708684,
    ];
    assert_eq!(probs.len(), 6);
    for (actual, expected) in probs.iter().zip(good_probs.iter()) {
        near(*actual, *expected, 0.1);
    }
    // Derived: the entries carry the distribution's own mass, so they sum to
    // `sum()` up to rounding.
    let total: f64 = probs.iter().sum();
    near(total, stats.sum(), 1e-9);

    // The source's size-less overload is this one with the existing length.
    let probs2 = stats.normal_approximation(6).unwrap();
    assert_eq!(probs, probs2);
}

// L354 normalApproximation(probability)
#[test]
fn section_normal_approximation_in_place() {
    // The source marks this NOT_TESTABLE and refers to the sized overload; the
    // port has a single function, so the same call covers both.
    let mut stats = BasicStatistics::new();
    stats.update(&[0.0, 1.0, 3.0, 2.0, 0.0]).unwrap();
    let probs = stats.normal_approximation(6).unwrap();
    assert_eq!(probs.len(), 6);
    near(probs[2], 3.38207, 0.1);
}

// L361 normalApproximation(probability, coordinate)
#[test]
fn section_normal_approximation_at_coordinates() {
    let magic1 = 200usize;
    let magic2 = 100.0f64;

    let mut data = vec![0.0; magic1];
    data.extend_from_slice(&DVECTOR_DATA);

    let mut stats = BasicStatistics::new();
    stats.update(&data).unwrap();
    let fit = stats.normal_approximation(data.len() + magic1).unwrap();
    let mut stats2 = BasicStatistics::new();
    stats2.update(&fit).unwrap();

    near(stats.sum(), stats2.sum(), 0.1);
    near(stats.mean(), stats2.mean(), 0.1);
    near(stats.variance(), stats2.variance(), 0.1);

    // Sampling the same approximation 100x more densely scales the index-based
    // mean by 100 and the variance by 100^2.
    let mut pos2 = Vec::new();
    let mut i = 0.0f64;
    while i < fit.len() as f64 {
        pos2.push(i);
        i += 1.0 / magic2;
    }
    let fit2 = stats.normal_approximation_at(&pos2).unwrap();
    stats2.update(&fit2).unwrap();
    near(stats.sum(), stats2.sum(), 0.1);
    near(stats.mean(), stats2.mean() / magic2, 0.1);
    near(stats.variance(), stats2.variance() / magic2 / magic2, 0.1);
}

// Native: guards the source does not have.
#[test]
fn normal_approximation_refuses_a_degenerate_density() {
    let mut stats = BasicStatistics::new();
    stats.set_sum(1.0);
    stats.set_mean(-1.0);
    // A zero variance makes every density zero away from the mean, so the
    // normalising sum is zero and the source's `density / gaussSum` is 0/0.
    assert!(matches!(
        stats.normal_approximation(4).unwrap_err(),
        Error::InvalidValue(_)
    ));
    assert!(matches!(
        stats.normal_approximation_at(&[1.0, 2.0]).unwrap_err(),
        Error::InvalidValue(_)
    ));
}

// Derived from `normalApproximationHelper_` (BasicStatistics.h:236-252): for
// `size == 0` both `for (i = 0; i < size; ++i)` loops run zero times, so
// `gaussSum` stays 0 but nothing ever divides by it and `probability` is left
// empty. The source returns normally, and so must this — an empty request is
// not a degenerate density.
#[test]
fn normal_approximation_of_nothing_is_empty_not_an_error() {
    let mut stats = BasicStatistics::new();
    stats.update(&[0.0, 1.0, 3.0, 2.0, 0.0]).unwrap();
    assert!(stats.normal_approximation(0).unwrap().is_empty());
    assert!(stats.normal_approximation_at(&[]).unwrap().is_empty());

    // The same holds when the density itself is degenerate: the source never
    // reaches the division, so there is nothing to refuse.
    let mut degenerate = BasicStatistics::new();
    degenerate.set_sum(1.0);
    degenerate.set_mean(0.0);
    assert!(degenerate.normal_approximation(0).unwrap().is_empty());
    assert!(degenerate.normal_approximation_at(&[]).unwrap().is_empty());

    // And on a default-constructed instance, which is the source's own
    // `normalApproximation(probability)` on an empty container.
    let fresh = BasicStatistics::new();
    assert!(fresh.normal_approximation(0).unwrap().is_empty());
    assert!(fresh.normal_approximation_at(&[]).unwrap().is_empty());
}

// Native: the debugging stream operator.
#[test]
fn display_reports_all_three_parameters() {
    let mut stats = BasicStatistics::new();
    stats.set_mean(1.0);
    stats.set_variance(2.0);
    stats.set_sum(3.0);
    assert_eq!(
        stats.to_string(),
        "BasicStatistics:  mean=1  variance=2  sum=3"
    );
}
