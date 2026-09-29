// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `ML/REGRESSION` and `ML/RANSAC`, plus `Math::RandomShuffler`.
//!
//! Expected values are the class-test literals of core `bc9cc12`
//! (`LinearRegression_test.cpp`, `LinearRegressionWithoutIntercept_test.cpp`,
//! `QuadraticRegression_test.cpp`, `RANSACModelLinear_test.cpp`,
//! `RANSACModelQuadratic_test.cpp`) unless a test says otherwise; see
//! `docs/RANSAC_SUPPORT.md`.

// The class tests' literals are transcribed verbatim, including digits beyond
// f64 precision.
#![allow(clippy::excessive_precision)]

use openms::math::ml::linear_regression::LinearRegression;
use openms::math::ml::linear_regression_without_intercept::LinearRegressionWithoutIntercept;
use openms::math::ml::quadratic_regression::QuadraticRegression;
use openms::math::ml::ransac::{
    DPair, Ransac, RansacModel, RansacModelLinear, RansacModelQuadratic, RansacParam,
};
use openms::math::random_shuffler::RandomShuffler;

/// The class tests' `TEST_REAL_SIMILAR` with its default tolerances.
fn real_similar(a: f64, b: f64) {
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

// ---------------------------------------------------------------- regression

fn linear_axes() -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let x: Vec<f64> = (0..10).map(f64::from).collect();
    let y = x.iter().map(|i| 2.0 * i + 4.0).collect();
    let y0 = x.iter().map(|i| 2.0 * i).collect();
    let w = x.iter().map(|i| 1.0 + i).collect();
    (x, y, y0, w)
}

#[test]
fn linear_regression_class_test() {
    let (x, y, y0, w) = linear_axes();
    let fresh = LinearRegression::new();
    assert_eq!(fresh.slope(), 0.0);
    assert_eq!(fresh.r_squared(), 0.0);

    let mut lin_reg = LinearRegression::new();
    let mut lin_reg2 = LinearRegression::new();
    lin_reg.compute_regression(0.95, &x, &y, true).unwrap();
    real_similar(lin_reg.slope(), 2.0);
    real_similar(lin_reg.intercept(), 4.0);
    real_similar(lin_reg.chi_squared(), 0.0);
    lin_reg2.compute_regression(0.95, &x, &y0, true).unwrap();
    real_similar(lin_reg2.slope(), 2.0);
    real_similar(lin_reg2.intercept(), 0.0);
    real_similar(lin_reg2.chi_squared(), 0.0);

    lin_reg
        .compute_regression_weighted(0.95, &x, &y, &w, false)
        .unwrap();
    real_similar(lin_reg.slope(), 2.0);
    real_similar(lin_reg.intercept(), 4.0);
    lin_reg2
        .compute_regression_weighted(0.95, &x, &y0, &w, false)
        .unwrap();
    real_similar(lin_reg2.slope(), 2.0);
    real_similar(lin_reg2.intercept(), 0.0);

    lin_reg
        .compute_regression_weighted(0.95, &x, &y, &w, true)
        .unwrap();
    real_similar(lin_reg.chi_squared(), 0.0);
    real_similar(lin_reg.intercept(), 4.0);
    real_similar(lin_reg.lower(), -2.0);
    real_similar(lin_reg.upper(), -2.0);
    real_similar(lin_reg.slope(), 2.0);
    real_similar(lin_reg.stand_dev_res(), 0.0);
    real_similar(lin_reg.stand_err_slope(), 0.0);
    real_similar(lin_reg.r_squared(), 1.0);
    real_similar(lin_reg.t_value(), 2.306);
    real_similar(lin_reg.x_intercept(), -2.0);
    real_similar(lin_reg.rsd(), 0.0);
    real_similar(lin_reg.mean_res(), 0.0);
}

#[test]
fn linear_regression_t_quantile_matches_published_tables() {
    // Independently derived: two-sided 95 % Student's t critical values from
    // the standard table (e.g. NIST/SEMATECH e-Handbook 1.3.6.7.2), n - 2
    // degrees of freedom. Five-digit table precision.
    for (n, expected) in [
        (3, 12.706),
        (4, 4.303),
        (7, 2.571),
        (12, 2.228),
        (32, 2.042),
    ] {
        let x: Vec<f64> = (0..n).map(f64::from).collect();
        let y: Vec<f64> = x.iter().map(|v| 3.0 * v + 1.0 + (v * 7.0).sin()).collect();
        let mut regression = LinearRegression::new();
        regression.compute_regression(0.95, &x, &y, true).unwrap();
        assert!(
            (regression.t_value() - expected).abs() < 5e-4,
            "n = {n}: {} vs {expected}",
            regression.t_value()
        );
    }
}

#[test]
fn linear_regression_unable_to_fit_and_skipped_goodness() {
    let mut regression = LinearRegression::new();
    // One distinct x: the height-line fit reports failure (UnableToFit).
    let error = regression
        .compute_regression(0.95, &[2.0, 2.0, 2.0], &[1.0, 2.0, 3.0], true)
        .unwrap_err();
    assert!(error.to_string().contains("UnableToFit-LinearRegression"));
    assert!(regression.compute_regression(0.95, &[], &[], true).is_err());
    assert!(
        regression
            .compute_regression(0.95, &[1.0, 2.0], &[1.0], true)
            .is_err()
    );
    // Two points: goodness statistics are skipped and R-squared keeps 0.
    regression
        .compute_regression(0.95, &[1.0, 3.0], &[2.0, 6.0], true)
        .unwrap();
    assert_eq!(regression.slope(), 2.0);
    assert_eq!(regression.intercept(), 0.0);
    assert_eq!(regression.r_squared(), 0.0);
    assert_eq!(LinearRegression::compute_point_y(3.0, 2.0, 1.0), 7.0);
}

#[test]
fn linear_regression_without_intercept_class_test() {
    let empty = LinearRegressionWithoutIntercept::new();
    assert_eq!(empty.slope(), None);

    let mut lr = LinearRegressionWithoutIntercept::new();
    lr.add_data(1.0, 2.0);
    assert_eq!(lr.slope(), None);
    lr.add_data(2.0, 4.0);
    lr.add_data(3.0, 6.0);
    real_similar(lr.slope().unwrap(), 2.0);

    let mut lr = LinearRegressionWithoutIntercept::new();
    lr.add_data_slices(&[1.0, 2.0], &[2.1, 3.9]).unwrap();
    real_similar(lr.slope().unwrap(), 1.98);

    let mut lr1 = LinearRegressionWithoutIntercept::new();
    lr1.add_data(5.0, 10.0);
    assert_eq!(lr1.slope(), None);
    let mut lr2 = LinearRegressionWithoutIntercept::new();
    lr2.add_data(1.0, 3.0);
    lr2.add_data(2.0, 6.0);
    real_similar(lr2.slope().unwrap(), 3.0);

    // Native: mismatched slices are refused and change nothing.
    let before = lr2.clone();
    assert!(lr2.add_data_slices(&[1.0], &[]).is_err());
    assert_eq!(lr2, before);
    assert_eq!(lr2.len(), 2);
}

#[test]
fn quadratic_regression_class_test() {
    let x: Vec<f64> = (0..10).map(f64::from).collect();
    let y: Vec<f64> = x.iter().map(|i| 5.5 * i * i + 2.0 * i + 4.0).collect();
    let y0: Vec<f64> = x.iter().map(|i| 5.5 * i * i + 2.0 * i).collect();
    let w: Vec<f64> = x.iter().map(|i| 1.0 + i).collect();

    let fresh = QuadraticRegression::new();
    assert_eq!(
        (fresh.a(), fresh.b(), fresh.c(), fresh.chi_squared()),
        (0.0, 0.0, 0.0, 0.0)
    );

    let mut q_reg = QuadraticRegression::new();
    let mut q_reg2 = QuadraticRegression::new();
    q_reg.compute_regression(&x, &y).unwrap();
    real_similar(q_reg.a(), 4.0);
    real_similar(q_reg.b(), 2.0);
    real_similar(q_reg.c(), 5.5);
    real_similar(q_reg.chi_squared(), 0.0);
    q_reg2.compute_regression(&x, &y0).unwrap();
    real_similar(q_reg2.a(), 0.0);
    real_similar(q_reg2.b(), 2.0);
    real_similar(q_reg2.c(), 5.5);
    real_similar(q_reg2.chi_squared(), 0.0);

    q_reg.compute_regression_weighted(&x, &y, &w).unwrap();
    real_similar(q_reg.a(), 4.0);
    real_similar(q_reg.b(), 2.0);
    real_similar(q_reg.c(), 5.5);
    real_similar(q_reg.chi_squared(), 0.0);
    q_reg2.compute_regression_weighted(&x, &y0, &w).unwrap();
    real_similar(q_reg2.a(), 0.0);
    real_similar(q_reg2.b(), 2.0);
    real_similar(q_reg2.c(), 5.5);
    real_similar(q_reg2.chi_squared(), 0.0);

    let at = 100.0;
    real_similar(q_reg.eval(at), at * at * 5.5 + at * 2.0 + 4.0);
    real_similar(
        QuadraticRegression::eval_with(4.0, 2.0, 5.5, at),
        at * at * 5.5 + at * 2.0 + 4.0,
    );
}

#[test]
fn quadratic_regression_singular_system_is_unable_to_fit() {
    let mut regression = QuadraticRegression::new();
    regression
        .compute_regression(&[0.0, 1.0, 2.0], &[1.0, 2.0, 5.0])
        .unwrap();
    let before = regression.clone();
    // Two distinct x values cannot determine three coefficients.
    let error = regression
        .compute_regression(&[1.0, 1.0, 2.0, 2.0], &[1.0, 2.0, 3.0, 4.0])
        .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("UnableToFit-QuadraticRegression")
    );
    assert!(regression.compute_regression(&[], &[]).is_err());
    assert_eq!(regression, before);
}

// -------------------------------------------------------------------- RANSAC

fn linear_model_pairs() -> Vec<DPair> {
    vec![
        (7.66217066e+00, 3.32871078e+02),
        (1.88986378e+01, 8.41782838e+02),
        (1.43387751e+01, 6.48336013e+02),
        (1.04946477e+01, 5.30115032e+02),
        (2.40052860e+00, 1.36793947e+02),
        (2.65925164e+00, 1.38532208e+02),
        (7.00156815e+00, 3.03487855e+01),
        (1.76671412e+01, 7.67575677e+02),
        (1.02592601e+01, 5.32449429e+02),
        (1.29020672e+01, -1.74450591e+01),
        (2.66076055e-02, 1.78205080e+01),
        (1.87212750e+01, 8.59152499e+02),
        (1.81219758e+01, -5.79165989e-01),
        (5.27778174e+00, 1.88005119e+02),
        (4.56777946e+00, 1.61530045e+02),
        (2.82887267e+00, 1.64411907e+02),
        (5.77563248e+00, 2.69781852e+02),
        (1.08263921e+01, 4.65275655e+02),
        (9.61444550e+00, 3.82697907e+02),
        (5.34540857e+00, 2.56156813e+02),
    ]
}

fn linear_model_new_pairs() -> Vec<DPair> {
    vec![
        (1.20513989e+01, 5.42172984e+02),
        (1.68354224e+00, 1.23674095e+02),
        (4.64668635e+00, 2.61350113e+02),
        (8.13976269e+00, 3.24462812e+02),
        (1.04776397e+01, 4.04452477e+02),
        (1.56315091e+01, 6.95756737e+02),
        (1.27266524e+01, 6.53571377e+01),
        (1.33784812e+01, 3.03064682e+01),
        (9.73484306e+00, -1.55933991e+00),
        (1.29040386e+00, 4.19535249e+01),
        (1.36889336e+01, 5.37472495e+02),
        (3.37465643e+00, 1.52514434e+02),
        (2.86567552e+00, 5.62442618e+01),
        (1.63579656e+01, 8.41451166e+02),
        (2.01345432e+01, 8.57894838e+02),
        (1.62549940e+01, 7.15378774e+02),
        (5.79326803e+00, 2.69370208e+02),
        (2.04520306e+00, 8.66527618e+01),
        (1.16970916e+01, 6.05836392e+02),
        (8.68788731e+00, 9.52993526e+00),
        (2.79787727e+00, 1.08213952e+02),
        (1.95778572e+01, 1.39196902e+02),
        (1.69500204e-01, 3.09473207e+01),
        (1.17974170e+01, 2.51798532e+01),
        (4.67384259e+00, 2.30870376e+02),
        (1.41658478e+01, 5.86317425e+02),
        (5.00923637e+00, -1.86559595e+01),
        (9.87160022e+00, 4.61676941e+02),
        (1.14474730e+01, 4.83241860e+02),
        (3.79416666e+00, 1.64038065e+02),
    ]
}

#[test]
fn ransac_model_linear_class_test() {
    let model = RansacModelLinear;
    let coeff = model.rm_fit(&linear_model_pairs()).unwrap();
    real_similar(coeff[0], 46.03865245);
    real_similar(coeff[1], 31.20358812);
    let rss = model.rm_rss(&linear_model_pairs(), &coeff).unwrap();
    real_similar(rss, 864089.67832345);

    let inliers = model
        .rm_inliers(&linear_model_new_pairs(), &coeff, 7e3)
        .unwrap();
    let expected = [
        1.68354224e+00,
        4.64668635e+00,
        8.13976269e+00,
        1.04776397e+01,
        1.29040386e+00,
        1.36889336e+01,
        3.37465643e+00,
        2.86567552e+00,
        5.79326803e+00,
        2.04520306e+00,
        2.79787727e+00,
        1.69500204e-01,
        4.67384259e+00,
        1.14474730e+01,
        3.79416666e+00,
    ];
    assert_eq!(inliers.len(), 15);
    for (inlier, x) in inliers.iter().zip(expected) {
        real_similar(inlier.0, x);
    }
    // Too few coefficients are refused rather than read out of bounds.
    assert!(model.rm_rss(&linear_model_pairs(), &[1.0]).is_err());
    assert!(model.rm_inliers(&linear_model_pairs(), &[], 1.0).is_err());
}

const TX: [f64; 50] = [
    7.66217066,
    18.8986378,
    14.3387751,
    10.4946477,
    2.4005286,
    2.65925164,
    7.00156815,
    17.6671412,
    10.2592601,
    12.9020672,
    0.0266076055,
    18.721275,
    18.1219758,
    5.27778174,
    4.56777946,
    2.82887267,
    5.77563248,
    10.8263921,
    9.6144455,
    5.34540857,
    12.0513989,
    1.68354224,
    4.64668635,
    8.13976269,
    10.4776397,
    15.6315091,
    12.7266524,
    13.3784812,
    9.73484306,
    1.29040386,
    13.6889336,
    3.37465643,
    2.86567552,
    16.3579656,
    20.1345432,
    16.254994,
    5.79326803,
    2.04520306,
    11.6970916,
    8.68788731,
    2.79787727,
    19.5778572,
    0.169500204,
    11.797417,
    4.67384259,
    14.1658478,
    5.00923637,
    9.87160022,
    11.447473,
    3.79416666,
];
const TY: [f64; 50] = [
    332.871078,
    841.782838,
    648.336013,
    530.115032,
    136.793947,
    138.532208,
    30.3487855,
    767.575677,
    532.449429,
    17.4450591,
    17.820508,
    859.152499,
    0.579165989,
    188.005119,
    161.530045,
    164.411907,
    269.781852,
    465.275655,
    382.697907,
    256.156813,
    542.172984,
    123.674095,
    261.350113,
    324.462812,
    404.452477,
    695.756737,
    65.3571377,
    30.3064682,
    1.55933991,
    41.9535249,
    537.472495,
    152.514434,
    56.2442618,
    841.451166,
    857.894838,
    715.378774,
    269.370208,
    86.6527618,
    605.836392,
    9.52993526,
    108.213952,
    139.196902,
    30.9473207,
    25.1798532,
    230.870376,
    586.317425,
    18.6559595,
    461.676941,
    483.24186,
    164.038065,
];
const TY_OUT: [f64; 40] = [
    17.8205, 30.9473, 41.9535, 123.674, 86.6528, 136.794, 138.532, 108.214, 164.412, 56.2443,
    152.514, 164.038, 161.53, 261.35, 230.87, 188.005, 256.157, 269.782, 269.37, 332.871, 324.463,
    382.698, 461.677, 532.449, 404.452, 530.115, 465.276, 483.242, 605.836, 542.173, 537.472,
    586.317, 648.336, 695.757, 715.379, 841.451, 767.576, 859.152, 841.783, 857.895,
];

fn sorted(mut pairs: Vec<DPair>) -> Vec<DPair> {
    pairs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    pairs
}

#[test]
fn ransac_linear_class_test_is_seed_independent() {
    // The class test constructs `RANSAC<RansacModelLinear>` with its default
    // seed, `time(nullptr)`, and still expects exactly these 40 points: the
    // result must not depend on the seed. Several seeds, including extremes
    // and ones that differ only above bit 31, exercise that.
    let pairs = sorted(TX.iter().copied().zip(TY).collect());
    for seed in [0, 1, 42, 1_727_600_000, u64::MAX, 1 << 40] {
        let mut ransac: Ransac<RansacModelLinear> = Ransac::new(seed);
        let out = sorted(
            ransac
                .ransac(&pairs, 2, 1200, 100.0 * 100.0, 10, false)
                .unwrap(),
        );
        assert_eq!(out.len(), 40, "seed {seed}");
        for (point, y) in out.iter().zip(TY_OUT) {
            real_similar(point.1, y);
        }
    }
}

fn quadratic_pairs(outliers: bool) -> Vec<DPair> {
    let tx: Vec<f64> = (0..19).map(|i| 300.0 + 50.0 * f64::from(i)).collect();
    let ty = [
        327.899691196695,
        380.224568059509,
        432.60001240954,
        485.025329595228,
        537.500950780092,
        590.025780790517,
        642.599328692568,
        695.226226189721,
        747.901404527071,
        800.628123375242,
        853.398402241557,
        906.221690051125,
        959.098650822155,
        1012.02444287071,
        1064.99750095039,
        1118.02244136503,
        1171.10218429403,
        1224.22697636988,
        1277.3989501643,
    ];
    let tyo = [
        327.899129765285,
        380.22505352209,
        3107.92239745832,
        485.025787154647,
        4105.48012991713,
        590.022752890109,
        642.600208487572,
        695.225450449224,
        747.90063358956,
        800.627835435032,
        853.39930288619,
        906.221242296453,
        4676.76961970711,
        1012.02448829764,
        1765.7626054147,
        1118.02258525226,
        1171.09984462503,
        1224.22603819705,
        1277.40083805864,
    ];
    tx.into_iter()
        .zip(if outliers { tyo } else { ty })
        .collect()
}

#[test]
fn ransac_model_quadratic_class_test() {
    let model = RansacModelQuadratic;
    let pairs = quadratic_pairs(false);
    let coeff = model.rm_fit(&pairs).unwrap();
    real_similar(coeff[0], 15.0009);
    real_similar(coeff[1], 1.04);
    real_similar(coeff[2], 0.00001);
    let rss = model.rm_rss(&pairs, &coeff).unwrap();
    real_similar(rss, 5.2254915523925468e-005);

    assert_eq!(
        model.rm_inliers(&pairs, &coeff, 0.5).unwrap().len(),
        pairs.len()
    );
    let inliers = model
        .rm_inliers(&quadratic_pairs(true), &coeff, 0.5)
        .unwrap();
    assert_eq!(inliers.len(), 15);
    real_similar(inliers[2].0, 450.0);
    real_similar(inliers[3].0, 550.0);
    real_similar(inliers[10].0, 950.0);
    real_similar(inliers[11].0, 1050.0);
}

#[test]
fn ransac_quadratic_class_test_seed_zero() {
    let mut ransac: Ransac<RansacModelQuadratic> = Ransac::new(0);
    let out = sorted(
        ransac
            .ransac(&quadratic_pairs(true), 5, 50, 2.0, 3, false)
            .unwrap(),
    );
    assert_eq!(out.len(), 15);
    real_similar(out[2].0, 450.0);
    real_similar(out[3].0, 550.0);
    real_similar(out[10].0, 950.0);
    real_similar(out[11].0, 1050.0);
}

#[test]
fn ransac_preconditions() {
    let pairs = quadratic_pairs(true);
    let mut ransac: Ransac<RansacModelQuadratic> = Ransac::new(0);
    // n must be smaller than the number of points.
    assert!(ransac.ransac(&pairs, 19, 5, 2.0, 3, false).is_err());
    // A relative d of 100 % or more is a precondition violation.
    assert!(ransac.ransac(&pairs, 3, 5, 2.0, 100, true).is_err());
    assert!(RansacParam::new(3, 5, 2.0, 100, true).is_err());
    let param = RansacParam::new(3, 5, 2.0, 99, true).unwrap();
    assert_eq!(
        param.to_string(),
        "RANSAC param:\n  n: 3\n  k: 5 iterations\n  t: 2 threshold\n  d: 99 inliers\n\n"
    );
    assert_eq!(RansacParam::default().k, 0);
    // Non-finite points are refused.
    let mut bad = pairs.clone();
    bad[4].1 = f64::NAN;
    assert!(ransac.ransac(&bad, 3, 5, 2.0, 3, false).is_err());
    // The work ceiling is checked before anything runs.
    assert!(ransac.ransac(&pairs, 3, usize::MAX, 2.0, 3, false).is_err());
    // Zero iterations find nothing.
    assert!(
        ransac
            .ransac(&pairs, 3, 0, 2.0, 3, false)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn ransac_is_reproducible_and_param_alias_agrees() {
    let pairs = quadratic_pairs(true);
    let param = RansacParam::new(5, 50, 2.0, 3, false).unwrap();
    let mut a: Ransac<RansacModelQuadratic> = Ransac::new(7);
    let mut b: Ransac<RansacModelQuadratic> = Ransac::new(7);
    assert_eq!(
        a.ransac_with_param(&pairs, &param).unwrap(),
        b.ransac(&pairs, 5, 50, 2.0, 3, false).unwrap()
    );
    // The constructor keeps only the low 32 bits of the seed (the source's
    // `int` narrowing); set_seed keeps all 64.
    let mut narrowed: Ransac<RansacModelQuadratic> = Ransac::new((1 << 32) | 7);
    let mut reseeded: Ransac<RansacModelQuadratic> = Ransac::new(0);
    reseeded.set_seed(7);
    let reference = Ransac::<RansacModelQuadratic>::new(7).ransac(&pairs, 5, 50, 2.0, 3, false);
    assert_eq!(
        narrowed.ransac(&pairs, 5, 50, 2.0, 3, false).unwrap(),
        reference.as_ref().unwrap().clone()
    );
    assert_eq!(
        reseeded.ransac(&pairs, 5, 50, 2.0, 3, false).unwrap(),
        reference.unwrap()
    );
}

// ------------------------------------------------------------ RandomShuffler

#[test]
fn random_shuffler_seeding() {
    let mut a = RandomShuffler::new(4711);
    let mut b = RandomShuffler::default();
    assert_eq!(b, RandomShuffler::with_seed(RandomShuffler::DEFAULT_SEED));
    b.seed(4711);
    assert_eq!(a, b);
    let mut left: Vec<u8> = b"ABCDEFGHIJKLMNOPQRST".to_vec();
    let mut right = left.clone();
    a.portable_random_shuffle(&mut left);
    b.portable_random_shuffle(&mut right);
    assert_eq!(left, right);
    // Same permutation as the decoy generator's engine, whose unit tests pin
    // it against independent integer arithmetic for this seed.
    assert_eq!(left, b"IHEDJOFCGMSTAKRPNLBQ");
    let mut empty: [u8; 0] = [];
    let before = a.clone();
    a.portable_random_shuffle(&mut empty);
    a.portable_random_shuffle(&mut [1_u8]);
    assert_eq!(a, before);
}
