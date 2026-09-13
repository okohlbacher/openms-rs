// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Port of `MultipleTesting_test.cpp` (core SDK bc9cc12): one test per
//! `START_SECTION`, cited by source line.
//!
//! Most literals are tier 3, transcribed from the class test, which itself
//! took them from PyProphet's and the R `qvalue` package's test suites. Three
//! sections compare against whole CSV fixtures the upstream test ships —
//! `test_qvalue_ref_data.csv` is R `qvalue` output and `test_lfdr_ref_data.csv`
//! is PyProphet output — so they check the port against a third-party
//! implementation rather than against the C++. That is still tier 3: the
//! fixture was retained by the upstream test, not produced by running the C++
//! here.
//!
//! Where a value can be derived, it is, and the comment says so.
//!
//! The smoothing spline that `pi0_est` needs is supplied here, through
//! `Pi0Smoother`, by the crate's port of `BSplineSmoothingSpline`. That is the
//! seam `src/math/multiple_testing.rs` describes: `math` may not depend on
//! `processing`, so the composition happens in this test and in whatever
//! caller assembles the two. Reproducing the class test's pi0 literals through
//! it is what shows the seam carries the source's behaviour.

use openms::Error;
use openms::math::multiple_testing::{
    DEFAULT_SMOOTH_DF, LfdrOptions, LfdrTransform, Pi0Method, Pi0Smoother, compute_model_fdr, lfdr,
    p_emp, p_norm, pi0_est, q_value,
};
use openms::processing::spline::smoothing::BSplineSmoothingSpline;
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// The class test's `TOLERANCE_ABSOLUTE(1e-4)`, loosened to match "the
/// reference comparison used in pyprophet (4 decimal places)".
fn similar(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-4,
        "{actual} is not within 1e-4 of {expected}"
    );
}

/// The source's `BSplineSmoothingSpline spl(xs, ys, -1.0, smooth_df)` followed
/// by `spl.eval(max_lambda)`, with `!spl.ok()` mapped to `None`.
struct BSplineSmoother;

impl Pi0Smoother for BSplineSmoother {
    fn smooth_eval(&self, x: &[f64], y: &[f64], smooth_df: i32, at: f64) -> Option<f64> {
        BSplineSmoothingSpline::with_smoothing(x, y, -1.0, smooth_df)
            .ok()
            .and_then(|spline| spline.eval(at).ok())
    }
}

fn read_column(path: &Path, column: usize) -> Vec<f64> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut out = Vec::new();
    for line in text.lines().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split(',').map(|p| p.trim_matches('"')).collect();
        if parts.len() <= column {
            continue;
        }
        out.push(parts[column].parse().unwrap());
    }
    out
}

/// The rows of a fixture, sorted ascending by the p-value in column zero, as
/// every section that reads one does before comparing.
fn read_sorted_rows(path: &Path, columns: usize) -> Vec<Vec<f64>> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut rows: Vec<Vec<f64>> = Vec::new();
    for line in text.lines().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split(',').map(|p| p.trim_matches('"')).collect();
        if parts.len() < columns {
            continue;
        }
        rows.push(
            parts[..columns]
                .iter()
                .map(|p| p.parse().unwrap())
                .collect(),
        );
    }
    rows.sort_by(|a, b| a[0].total_cmp(&b[0]));
    rows
}

// L34 template<class T> std::vector<double> computeModelFDR(const std::vector<T>&)
#[test]
fn section_compute_model_fdr() {
    let fdr = compute_model_fdr(&[0.1, 0.2, 0.3]).unwrap();
    assert_eq!(fdr.len(), 3);
    similar(fdr[0], 0.1);
    similar(fdr[1], 0.15);
    similar(fdr[2], 0.2);

    // Unsorted input: the smallest PEP gets the smallest q, in its own slot.
    let unsorted = compute_model_fdr(&[0.3, 0.1, 0.2]).unwrap();
    assert_eq!(unsorted.len(), 3);
    similar(unsorted[0], 0.2);
    similar(unsorted[1], 0.1);
    similar(unsorted[2], 0.15);

    // Ties share one value, the cumulative sum at the last of the tie over the
    // max rank: 0.5/3 = 1/6.
    let tied = compute_model_fdr(&[0.1, 0.2, 0.2, 0.3]).unwrap();
    assert_eq!(tied.len(), 4);
    similar(tied[0], 0.1);
    similar(tied[1], 1.0 / 6.0);
    similar(tied[2], 1.0 / 6.0);
    similar(tied[3], 0.2);

    let both = compute_model_fdr(&[0.3, 0.2, 0.2, 0.1]).unwrap();
    assert_eq!(both.len(), 4);
    similar(both[0], 0.2);
    similar(both[1], 1.0 / 6.0);
    similar(both[2], 1.0 / 6.0);
    similar(both[3], 0.1);

    // Any NaN invalidates the whole result.
    let with_nan = compute_model_fdr(&[0.1, f64::NAN, 0.2]).unwrap();
    assert_eq!(with_nan.len(), 3);
    for value in &with_nan {
        assert!(value.is_nan());
    }

    assert!(compute_model_fdr(&[]).unwrap().is_empty());
}

// L84 "pNorm: pyprophet reference vector"
#[test]
fn section_p_norm_pyprophet_vector() {
    let stat = [0.0, 1.0, 3.0, 2.0, 0.1, 0.5, 0.6, 0.3, 0.5, 0.6, 0.2, 0.5];
    let stat0 = [0.4, 0.2, 0.5, 1.0, 0.5, 0.7, 0.2, 0.4];
    let out = p_norm(&stat, &stat0).unwrap();
    let expected = [
        9.674763e-01,
        2.621760e-02,
        0.0,
        5.201675e-09,
        9.287418e-01,
        4.811347e-01,
        3.351438e-01,
        7.610205e-01,
        4.811347e-01,
        3.351438e-01,
        8.617105e-01,
        4.811347e-01,
    ];
    assert_eq!(out.len(), expected.len());
    for (got, want) in out.iter().zip(expected.iter()) {
        similar(*got, *want);
    }
    // Derived: equal statistics must receive equal tails, and the tail must be
    // strictly decreasing in the statistic.
    assert_eq!(out[5], out[8]);
    assert_eq!(out[6], out[9]);
    for pair in [(0usize, 4usize), (4, 10), (10, 7), (7, 5), (5, 6), (6, 1)] {
        assert!(out[pair.0] > out[pair.1], "{pair:?}");
    }
}

// L96 "lfdr CSV reference (pyprophet) regression"
#[test]
fn section_lfdr_pyprophet_csv() {
    let rows = read_sorted_rows(&data("test_lfdr_ref_data.csv"), 5);
    assert!(!rows.is_empty());
    let p: Vec<f64> = rows.iter().map(|r| r[0]).collect();
    let pi0 = 0.669926026474838;

    let default = lfdr(&p, pi0, &LfdrOptions::default()).unwrap();
    let no_monotone = lfdr(
        &p,
        pi0,
        &LfdrOptions {
            monotone: false,
            ..LfdrOptions::default()
        },
    )
    .unwrap();
    let logit = lfdr(
        &p,
        pi0,
        &LfdrOptions {
            transform: LfdrTransform::Logit,
            ..LfdrOptions::default()
        },
    )
    .unwrap();
    let big_eps = lfdr(
        &p,
        pi0,
        &LfdrOptions {
            eps: 1e-2,
            ..LfdrOptions::default()
        },
    )
    .unwrap();

    assert_eq!(default.len(), rows.len());
    // The class test's tolerance, "match python test decimal=2".
    let tol = 1e-2;
    for (index, row) in rows.iter().enumerate() {
        assert!((default[index] - row[1]).abs() <= tol, "default[{index}]");
        assert!(
            (no_monotone[index] - row[2]).abs() <= tol,
            "monotone=false[{index}]"
        );
        assert!((logit[index] - row[3]).abs() <= tol, "logit[{index}]");
        assert!((big_eps[index] - row[4]).abs() <= tol, "eps[{index}]");
    }
}

// L153 "pemp: pyprophet reference vector"
#[test]
fn section_p_emp_pyprophet_vector() {
    let stat = [0.0, 1.0, 3.0, 2.0, 0.1, 0.5, 0.6, 0.3, 0.5, 0.6, 0.2, 0.5];
    let stat0 = [0.4, 0.2, 0.5, 1.0, 0.5, 0.7, 0.2, 0.4];
    let out = p_emp(&stat, &stat0).unwrap();
    let expected = [
        1.0, 0.125, 0.125, 0.125, 1.0, 0.25, 0.25, 0.75, 0.25, 0.25, 0.75, 0.25,
    ];
    assert_eq!(out.len(), expected.len());
    for (got, want) in out.iter().zip(expected.iter()) {
        similar(*got, *want);
    }
    // Derived: no p-value may fall below 1/|stat0| = 0.125.
    for value in &out {
        assert!(*value >= 0.125 - 1e-12);
    }
}

// L164 "lfdr basic checks (probit & logit)"
#[test]
fn section_lfdr_basic_checks() {
    let p = [0.001, 0.01, 0.05, 0.2, 0.8, 0.95];
    let pi0 = 0.8;
    let options = LfdrOptions {
        gridsize: 256,
        ..LfdrOptions::default()
    };

    let probit = lfdr(&p, pi0, &options).unwrap();
    assert_eq!(probit.len(), p.len());
    for value in &probit {
        assert!(value.is_finite());
        assert!(*value >= 0.0);
        assert!(*value <= 1.0);
    }
    let mut order: Vec<usize> = (0..p.len()).collect();
    order.sort_by(|a, b| p[*a].total_cmp(&p[*b]));
    for window in order.windows(2) {
        assert!(probit[window[1]] >= probit[window[0]]);
    }

    let logit = lfdr(
        &p,
        pi0,
        &LfdrOptions {
            transform: LfdrTransform::Logit,
            ..options
        },
    )
    .unwrap();
    assert_eq!(logit.len(), p.len());
    for value in &logit {
        assert!(value.is_finite());
        assert!(*value >= 0.0);
        assert!(*value <= 1.0);
    }

    // Native: a non-finite p-value is dropped and its slot comes back NaN,
    // while the finite ones are unaffected.
    let padded = [0.001, f64::NAN, 0.01, 0.05, 0.2, 0.8, 0.95];
    let out = lfdr(&padded, pi0, &options).unwrap();
    assert!(out[1].is_nan());
    assert_eq!(out[0], probit[0]);
    assert_eq!(out[2], probit[1]);

    // The two `std::invalid_argument` throws.
    assert!(matches!(
        lfdr(&[1.5], 0.8, &options),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        lfdr(&[0.5], 1.5, &options),
        Err(Error::InvalidValue(_))
    ));
}

// L189 "bw_nrd0: pyprophet reference value"
#[test]
fn section_bw_nrd0_pyprophet_value() {
    let stat = [0.0, 1.0, 3.0, 2.0, 0.1, 0.5, 0.6, 0.3, 0.5, 0.6, 0.2, 0.5];
    let bw = openms::math::kernel_density::bw_nrd0(&stat).unwrap();
    similar(bw, 0.1736562);
}

// L197 "forrt/revrt roundtrip small"
#[test]
fn section_for_rt_rev_rt_roundtrip() {
    use openms::math::kernel_density::{for_rt, rev_rt};
    let x = [0.0, 1.0, 2.0, 3.0];
    let packed = for_rt(&x, x.len()).unwrap();
    let restored = rev_rt(&packed, x.len()).unwrap();
    assert_eq!(restored.len(), x.len());
    for (got, want) in restored.iter().zip(x.iter()) {
        similar(*got, *want);
    }
}

// L210 "pi0est: pyprophet reference checks"
#[test]
fn section_pi0_est_pyprophet_checks() {
    let mut p = read_column(&data("test_lfdr_ref_data.csv"), 0);
    p.sort_by(f64::total_cmp);

    // A single lambda needs no spline at all.
    let single = pi0_est(
        &p,
        &[0.4],
        Pi0Method::Smoother,
        DEFAULT_SMOOTH_DF,
        false,
        None,
    )
    .unwrap();
    similar(single.pi0, 0.697161);
    assert!(!single.pi0_smooth);

    // The default lambda grid with the smoother.
    let smoothed = pi0_est(
        &p,
        &[],
        Pi0Method::Smoother,
        DEFAULT_SMOOTH_DF,
        false,
        Some(&BSplineSmoother),
    )
    .unwrap();
    similar(smoothed.pi0, 0.6685638);
    assert!(smoothed.pi0_smooth);
    assert_eq!(smoothed.lambda.len(), 19);
    // Derived: the source builds the grid by repeated addition, so the third
    // entry carries the accumulated rounding.
    assert_eq!(smoothed.lambda[2], 0.15000000000000002);

    // lambda 0.4..0.95 step 0.05 with smooth_log_pi0.
    let mut lambda = Vec::new();
    let mut l = 0.4;
    while l < 1.0 - 1e-12 {
        lambda.push(l);
        l += 0.05;
    }
    let log_smoothed = pi0_est(
        &p,
        &lambda,
        Pi0Method::Smoother,
        DEFAULT_SMOOTH_DF,
        true,
        Some(&BSplineSmoother),
    )
    .unwrap();
    similar(log_smoothed.pi0, 0.6658949);
    assert!(log_smoothed.pi0_smooth);

    // Without a smoother the source's `!spl.ok()` fallback is taken, which is
    // the conservative min over the per-lambda estimates.
    let fallback = pi0_est(&p, &[], Pi0Method::Smoother, DEFAULT_SMOOTH_DF, false, None).unwrap();
    assert!(!fallback.pi0_smooth);
    let minimum = smoothed
        .pi0_lambda
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min)
        .min(1.0);
    assert_eq!(fallback.pi0, minimum);

    // Bootstrap needs no spline; it must land on one of the per-lambda values.
    let bootstrap = pi0_est(
        &p,
        &[],
        Pi0Method::Bootstrap,
        DEFAULT_SMOOTH_DF,
        false,
        None,
    )
    .unwrap();
    assert!(!bootstrap.pi0_smooth);
    assert!(
        bootstrap
            .pi0_lambda
            .iter()
            .any(|v| v.min(1.0) == bootstrap.pi0)
    );
}

// L248 double bwNrd0(const std::vector<double>&)
#[test]
fn section_bw_nrd0_five_point_ladder() {
    let x = [0.0, 1.0, 2.0, 3.0, 4.0];
    let bw = openms::math::kernel_density::bw_nrd0(&x).unwrap();
    similar(bw, 0.9735846228506357);
    // Derived: sd = sqrt(2.5), IQR = 2 so IQR/1.34 = 1.4925... < sd, and the
    // rule is 0.9 * (2/1.34) * 5^-0.2.
    similar(bw, 0.9 * (2.0 / 1.34) * 5.0_f64.powf(-0.2));
}

// L256 std::vector<double> linBin(...)
#[test]
fn section_lin_bin_histogram() {
    use openms::math::kernel_density::lin_bin;
    let x = [0.1, 0.4, 0.9];
    let bins = lin_bin(&x, 0.0, 1.0, 5, None).unwrap();
    assert_eq!(bins.len(), 5);
    similar(bins[0], 1.0);
    similar(bins[1], 0.0);
    similar(bins[2], 1.0);
    similar(bins[3], 0.0);
    similar(bins[4], 1.0);

    let weights = [2.0, 3.0, 5.0];
    let weighted = lin_bin(&x, 0.0, 1.0, 5, Some(&weights)).unwrap();
    similar(weighted[0], 2.0);
    similar(weighted[2], 3.0);
    similar(weighted[4], 5.0);
}

// L274 std::vector<double> pNorm(const std::vector<double>&, const std::vector<double>&)
#[test]
fn section_p_norm_small() {
    let stat0 = [0.0, 1.0];
    let stat = [0.5, 1.5, -0.5, f64::NAN];
    let out = p_norm(&stat, &stat0).unwrap();
    assert_eq!(out.len(), stat.len());

    // mu = 0.5, sample variance 0.5 with one degree of freedom removed.
    let mu = 0.5;
    let sigma = 0.5_f64.sqrt();
    let sqrt2 = 2.0_f64.sqrt();
    for index in 0..3 {
        let z = (stat[index] - mu) / sigma;
        let expected = 1.0 - 0.5 * (1.0 + libm::erf(z / sqrt2));
        similar(out[index], expected);
    }
    assert!(out[3].is_nan());
    // Derived: the statistic equal to the mean has an upper tail of exactly a
    // half, by symmetry.
    similar(out[0], 0.5);

    // The two `std::invalid_argument` throws.
    assert!(matches!(p_norm(&stat, &[]), Err(Error::InvalidValue(_))));
    assert!(matches!(
        p_norm(&stat, &[f64::NAN, f64::INFINITY]),
        Err(Error::InvalidValue(_))
    ));
    // A degenerate null is a point mass: 1 below the mean, 0 at or above it.
    let degenerate = p_norm(&[0.5, 1.0, 1.5], &[1.0, 1.0]).unwrap();
    assert_eq!(degenerate, vec![1.0, 0.0, 0.0]);
}

// L299 template<class T> std::vector<double> pEmp(const std::vector<T>&, const std::vector<T>&)
#[test]
fn section_p_emp_small() {
    let p = p_emp(&[0.9, 0.8], &[0.1, 0.2]).unwrap();
    assert_eq!(p.len(), 2);
    similar(p[0], 0.5);
    similar(p[1], 0.5);
    // Both raw p-values are zero and the floor of 1/|stat0| = 0.5 lifts them.
    assert!(matches!(p_emp(&[], &[0.1]), Err(Error::InvalidValue(_))));
    assert!(matches!(p_emp(&[0.1], &[]), Err(Error::InvalidValue(_))));
}

// L308 std::vector<double> qValue(const std::vector<double>&, double, bool)
#[test]
fn section_q_value_small() {
    let p = [0.01, 0.02, 0.03];
    let q = q_value(&p, 1.0, false).unwrap();
    assert_eq!(q.len(), p.len());
    similar(q[0], 0.03);
    similar(q[1], 0.03);
    similar(q[2], 0.03);
    // Derived: the raw values are 3*0.01/1, 3*0.02/2 and 3*0.03/3, all exactly
    // 0.03, so the monotone sweep changes nothing.
    for value in &q {
        assert!((*value - 0.03).abs() <= 1e-15);
    }

    // Non-finite entries are dropped and their slots come back NaN.
    let with_nan = q_value(&[0.01, f64::NAN, 0.02, 0.03], 1.0, false).unwrap();
    assert!(with_nan[1].is_nan());
    similar(with_nan[0], 0.03);
    assert!(q_value(&[], 1.0, false).unwrap().is_empty());

    // The two `std::invalid_argument` throws.
    assert!(matches!(
        q_value(&[1.5], 1.0, false),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        q_value(&[0.5], 1.5, false),
        Err(Error::InvalidValue(_))
    ));
}

// L318 Pi0Result pi0Est(const std::vector<double>&, const std::vector<double>&)
#[test]
fn section_pi0_est_clipped_to_one() {
    let p = [0.9, 0.8, 0.85, 0.95];
    let lambda = [0.5, 0.6];
    let result = pi0_est(
        &p,
        &lambda,
        Pi0Method::Smoother,
        DEFAULT_SMOOTH_DF,
        false,
        Some(&BSplineSmoother),
    )
    .unwrap();
    assert_eq!(result.pi0, 1.0);
    // Two lambdas is fewer than the four the smoother needs, so the fallback
    // runs and the estimate is min(pi0_lambda, 1). Every p exceeds both
    // lambdas, so the raw estimates are 1/(1-0.5) = 2 and 1/(1-0.6) = 2.5.
    assert!(!result.pi0_smooth);
    similar(result.pi0_lambda[0], 2.0);
    similar(result.pi0_lambda[1], 2.5);

    // No finite p-value at all is the third `std::invalid_argument`.
    assert!(matches!(
        pi0_est(
            &[f64::NAN],
            &[],
            Pi0Method::Smoother,
            DEFAULT_SMOOTH_DF,
            false,
            None
        ),
        Err(Error::InvalidValue(_))
    ));
    // A lambda outside [0, 1) likewise.
    assert!(matches!(
        pi0_est(
            &[0.5],
            &[1.0],
            Pi0Method::Smoother,
            DEFAULT_SMOOTH_DF,
            false,
            None
        ),
        Err(Error::InvalidValue(_))
    ));
}

// L327 "R reference qvalue CSV matches C++ implementation"
#[test]
fn section_q_value_r_reference_csv() {
    let rows = read_sorted_rows(&data("test_qvalue_ref_data.csv"), 3);
    assert!(!rows.is_empty());
    let p: Vec<f64> = rows.iter().map(|r| r[0]).collect();
    let pi0 = 0.669926026474838;

    let default = q_value(&p, pi0, false).unwrap();
    let positive = q_value(&p, pi0, true).unwrap();
    assert_eq!(default.len(), rows.len());
    assert_eq!(positive.len(), rows.len());
    for (index, row) in rows.iter().enumerate() {
        similar(default[index], row[1]);
        similar(positive[index], row[2]);
    }
    // Derived: q is non-decreasing in p and never above one.
    for window in default.windows(2) {
        assert!(window[1] >= window[0] - 1e-12);
    }
    assert!(default.iter().all(|v| *v <= 1.0));
}

// L385 "silvermanKernelFFT properties"
#[test]
fn section_silverman_kernel_properties() {
    use openms::math::kernel_density::silverman_kernel_fft;
    let m = 8usize;
    let kernel = silverman_kernel_fft(1.0, m, 10.0).unwrap();
    assert_eq!(kernel.len(), m);
    similar(kernel[0], 1.0);
    let half = m / 2;
    for k in 0..=half {
        assert!(kernel[k].is_finite());
        assert!(kernel[k] > 0.0);
        if k > 0 {
            assert!(kernel[k] <= kernel[k - 1]);
        }
    }
    for k in 1..half {
        similar(kernel[k], kernel[half + k]);
    }
}

// L408 "gridKdeFFT basic integration and sizes"
#[test]
fn section_grid_kde_basic_integration() {
    use openms::math::kernel_density::{bw_nrd0, grid_kde_fft};
    let x = [0.0, 1.0, 2.0];
    let bw = bw_nrd0(&x).unwrap();
    let result = grid_kde_fft(&x, bw, 64, 3.0).unwrap();
    assert_eq!(result.density.len(), result.grid.len());
    assert!(!result.density.is_empty());
    let delta = if result.grid.len() > 1 {
        result.grid[1] - result.grid[0]
    } else {
        1.0
    };
    let mut sum = 0.0;
    for value in &result.density {
        assert!(value.is_finite());
        assert!(*value >= 0.0);
        sum += *value;
    }
    similar(sum * delta, 1.0);
}

// L425 "kdeFFTEval matches spline-interpolated grid density at sample points"
#[test]
fn section_kde_eval_matches_grid_spline() {
    use openms::math::kernel_density::{bw_nrd0, grid_kde_fft, kde_fft_eval};
    use openms::processing::spline::cubic::CubicSpline2d;
    let x = [0.0, 1.0, 2.0];
    let bw = bw_nrd0(&x).unwrap();
    let grid = grid_kde_fft(&x, bw, 64, 3.0).unwrap();
    let spline = CubicSpline2d::new(&grid.grid, &grid.density).unwrap();
    let expected: Vec<f64> = x.iter().map(|v| spline.eval(*v).unwrap()).collect();
    let got = kde_fft_eval(&x, bw, 64, 3.0).unwrap();
    assert_eq!(got.len(), expected.len());
    for (got, want) in got.iter().zip(expected.iter()) {
        similar(*got, *want);
    }
}

/// Native: the enum conversions the source spells as free functions over
/// strings, including the two `std::invalid_argument` throws.
#[test]
fn enum_conversions_round_trip() {
    assert_eq!(Pi0Method::Smoother.as_str(), "smoother");
    assert_eq!(Pi0Method::Bootstrap.as_str(), "bootstrap");
    assert_eq!(Pi0Method::parse("SMOOTHER").unwrap(), Pi0Method::Smoother);
    assert_eq!(Pi0Method::parse("Bootstrap").unwrap(), Pi0Method::Bootstrap);
    assert!(matches!(
        Pi0Method::parse("spline"),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(Pi0Method::default(), Pi0Method::Smoother);

    assert_eq!(LfdrTransform::Probit.as_str(), "probit");
    assert_eq!(LfdrTransform::Logit.as_str(), "logit");
    assert_eq!(
        LfdrTransform::parse("PROBIT").unwrap(),
        LfdrTransform::Probit
    );
    assert_eq!(LfdrTransform::parse("logit").unwrap(), LfdrTransform::Logit);
    assert!(matches!(
        LfdrTransform::parse("cloglog"),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(LfdrTransform::default(), LfdrTransform::Probit);

    let defaults = LfdrOptions::default();
    assert!(defaults.truncate);
    assert!(defaults.monotone);
    assert_eq!(defaults.transform, LfdrTransform::Probit);
    assert_eq!(defaults.adj, 1.5);
    assert_eq!(defaults.eps, 1e-8);
    assert_eq!(defaults.gridsize, 512);
    assert_eq!(defaults.cut, 3.0);
}
