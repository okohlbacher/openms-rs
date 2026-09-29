// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `ANALYSIS/OPENSWATH/MRMRTNormalizer`.
//!
//! Expected values are the literals of `MRMRTNormalizer_test.cpp` at core
//! `bc9cc12` unless a test says otherwise. The class test does not cover
//! `removeOutliersRANSAC` or `computeBinnedCoverage`; their tests here are
//! independently derived and say so.

use openms::analysis::openswath::mrm_rt_normalizer::{
    OutlierMethod, chauvenet, chauvenet_probability, compute_binned_coverage,
    jackknife_outlier_candidate, remove_outliers_iterative, remove_outliers_ransac,
    residual_outlier_candidate,
};

fn real_similar(a: f64, b: f64) {
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

const X1: [f64; 6] = [1.1, 2.0, 3.3, 3.9, 4.9, 6.2];
const Y1: [f64; 6] = [0.9, 1.9, 3.0, 3.7, 5.2, 6.1];
const X2: [f64; 6] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
const Y_OUTLIER: [f64; 6] = [0.9, 1.9, 7.0, 3.7, 5.2, 6.1];

#[test]
fn jackknife_outlier_candidate_class_test() {
    assert_eq!(jackknife_outlier_candidate(&X1, &Y1).unwrap(), 4);
    assert_eq!(jackknife_outlier_candidate(&X2, &X2).unwrap(), 0);
}

#[test]
fn residual_outlier_candidate_class_test() {
    assert_eq!(residual_outlier_candidate(&X1, &Y1).unwrap(), 4);
    assert_eq!(residual_outlier_candidate(&X2, &X2).unwrap(), 0);
}

fn pairs(x: &[f64], y: &[f64]) -> Vec<(f64, f64)> {
    x.iter().copied().zip(y.iter().copied()).collect()
}

fn thirty() -> Vec<(f64, f64)> {
    let x3 = [
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
        17.0, 18.0, 19.0, 1.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0,
    ];
    let y3 = [
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
        17.0, 18.0, 19.0, 20.0, 1.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0,
    ];
    pairs(&x3, &y3)
}

#[test]
fn remove_outliers_iterative_class_test() {
    let input1 = pairs(&X1, &Y1);
    let output1 =
        remove_outliers_iterative(&input1, 0.9, 0.5, true, OutlierMethod::IterResidual).unwrap();
    assert_eq!(output1.len(), input1.len());

    let input2 = pairs(&X1, &Y_OUTLIER);
    let input3 = thirty();
    for (use_chauvenet, method) in [
        (true, OutlierMethod::IterResidual),
        (false, OutlierMethod::IterResidual),
        (false, OutlierMethod::IterJackknife),
    ] {
        let output2 = remove_outliers_iterative(&input2, 0.9, 0.5, use_chauvenet, method).unwrap();
        assert_eq!(output2.len(), input2.len() - 1);
        assert_eq!(output2[0].0, input2[0].0);
        assert_eq!(output2[1].1, input2[1].1);
        assert_eq!(output2[2].0, input2[3].0);
        assert_eq!(output2[3].1, input2[4].1);

        let output3 = remove_outliers_iterative(&input3, 0.9, 0.2, use_chauvenet, method).unwrap();
        assert_eq!(output3.len(), input3.len() - 2);
        assert_eq!(output3[18].0, input3[18].0);
        assert_eq!(output3[19].1, input3[21].1);
    }
}

#[test]
fn remove_outliers_iterative_errors() {
    assert!(
        remove_outliers_iterative(
            &pairs(&X1[..2], &Y1[..2]),
            0.9,
            0.5,
            true,
            OutlierMethod::IterResidual
        )
        .is_err()
    );
    // An R-squared limit above 1 can never be met: UnableToFit once coverage
    // stops the loop.
    let error = remove_outliers_iterative(
        &pairs(&X1, &Y_OUTLIER),
        1.1,
        0.9,
        false,
        OutlierMethod::IterResidual,
    )
    .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("UnableToFit-LinearRegression-RTNormalizer")
    );
    assert_eq!(
        "iter_jackknife".parse::<OutlierMethod>().unwrap(),
        OutlierMethod::IterJackknife
    );
    assert_eq!(
        "iter_residual".parse::<OutlierMethod>().unwrap().name(),
        "iter_residual"
    );
    assert!("iter_other".parse::<OutlierMethod>().is_err());
}

const DATA1: [f64; 15] = [
    1.0, 2.0, 3.0, 4.0, 2.0, 10.0, 11.0, 75.0, 5.0, 8.0, 3.0, 5.0, 6.0, 9.0, 130.0,
];

#[test]
fn chauvenet_probability_class_test() {
    let expected = [
        0.61831553, 0.6387955, 0.65955473, 0.68057951, 0.6387955, 0.81146293, 0.8339146,
        0.10161557, 0.70185552, 0.76703896, 0.65955473, 0.70185552, 0.72336784, 0.78916526,
        0.00126358,
    ];
    for (pos, value) in expected.iter().enumerate() {
        real_similar(chauvenet_probability(&DATA1, pos).unwrap(), *value);
    }
    assert!(chauvenet_probability(&DATA1, 15).is_err());
    assert!(chauvenet_probability(&[], 0).is_err());
}

#[test]
fn chauvenet_class_test() {
    for pos in 0..14 {
        assert!(!chauvenet(&DATA1, pos).unwrap(), "position {pos}");
    }
    assert!(chauvenet(&DATA1, 14).unwrap());
}

/// Thirty points on `y = 2x + 5` with three gross outliers, `x` spread over
/// `[0, 100)`.
fn rt_pairs() -> Vec<(f64, f64)> {
    let mut data: Vec<(f64, f64)> = (0..30)
        .map(|i| {
            let x = f64::from(i) * 3.3 + 0.7;
            // Deterministic small noise in [-0.5, 0.5).
            let noise = f64::from((i * 37) % 11) / 11.0 - 0.5;
            (x, 2.0 * x + 5.0 + noise)
        })
        .collect();
    data[4].1 += 80.0;
    data[17].1 -= 120.0;
    data[25].1 += 55.0;
    data
}

#[test]
fn remove_outliers_ransac_drops_gross_outliers() {
    // Independently derived: the three displaced points deviate from the line
    // by 55 or more, every other point by at most 0.5, so with a deviation
    // threshold of 5 RANSAC keeps exactly the 27 clean points whatever the
    // seed, and their R-squared exceeds 0.999.
    let data = rt_pairs();
    for seed in [0, 1, 1_727_600_000] {
        let kept = remove_outliers_ransac(&data, 0.99, 0.6, 1000, 5.0, 5, seed).unwrap();
        assert_eq!(kept.len(), 27, "seed {seed}");
        for outlier in [data[4], data[17], data[25]] {
            assert!(!kept.contains(&outlier));
        }
    }
}

#[test]
fn remove_outliers_ransac_limits() {
    let data = rt_pairs();
    // Fewer than 5 sampled peptides.
    let error = remove_outliers_ransac(&data, 0.9, 0.6, 100, 5.0, 4, 0).unwrap_err();
    assert!(error.to_string().contains("below limit of 5 peptides"));
    // Fewer than 30 input peptides.
    let error = remove_outliers_ransac(&data[..29], 0.9, 0.6, 100, 5.0, 5, 0).unwrap_err();
    assert!(error.to_string().contains("below limit of 30 peptides"));
    // Coverage demands more points than the clean ones (27 < floor(0.95*30)).
    // RANSAC then finds no set with more than d = 28 inliers besides the
    // sample, and the empty result cannot be fitted.
    assert!(remove_outliers_ransac(&data, 0.9, 0.95, 200, 5.0, 5, 0).is_err());
    // R-squared limit above 1 is never met.
    let error = remove_outliers_ransac(&data, 1.5, 0.6, 200, 5.0, 5, 0).unwrap_err();
    assert!(error.to_string().contains("rsq"));
    assert!(remove_outliers_ransac(&data, 0.9, f64::NAN, 200, 5.0, 5, 0).is_err());
    assert!(remove_outliers_ransac(&data, 0.9, -0.1, 200, 5.0, 5, 0).is_err());
}

#[test]
fn compute_binned_coverage_counts_bins() {
    // Independently derived from the definition. Range [0, 100), 4 bins of
    // width 25; theoretical RTs (second members) 1, 2, 30, 31, 32, 99 and
    // 100 (the range end, clamped into the last bin).
    let data: Vec<(f64, f64)> = [1.0, 2.0, 30.0, 31.0, 32.0, 99.0, 100.0]
        .iter()
        .map(|rt| (0.0, *rt))
        .collect();
    // Bin counts: [2, 3, 0, 2].
    assert!(compute_binned_coverage((0.0, 100.0), &data, 4, 2, 3).unwrap());
    assert!(!compute_binned_coverage((0.0, 100.0), &data, 4, 2, 4).unwrap());
    assert!(!compute_binned_coverage((0.0, 100.0), &data, 4, 3, 2).unwrap());
    assert!(compute_binned_coverage((0.0, 100.0), &data, 4, 3, 1).unwrap());
    // Slightly below the range still truncates into the first bin.
    assert!(compute_binned_coverage((0.0, 100.0), &[(0.0, -1.0)], 4, 1, 1).unwrap());
    // Far below the range, a zero-width range and no bins are refused.
    assert!(compute_binned_coverage((0.0, 100.0), &[(0.0, -30.0)], 4, 1, 1).is_err());
    assert!(compute_binned_coverage((5.0, 5.0), &[(0.0, 5.0)], 4, 1, 1).is_err());
    assert!(compute_binned_coverage((0.0, 100.0), &data, 0, 1, 1).is_err());
}
