// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Port of `KernelDensityEstimation_test.cpp` (core SDK bc9cc12): one test per
//! `START_SECTION`, cited by source line.
//!
//! The literals are tier 3 — transcribed from the class test — except where a
//! comment says otherwise. Two sections carry a stronger oracle: the
//! `kde_reference_data.csv` fixture is `statsmodels` output the class test
//! itself compares against, and the scipy section's vector was produced by
//! `scipy.stats.gaussian_kde`. Both are third-party reference data retained by
//! the upstream test, not output of the C++ under test, so they are tier 3 as
//! well; what they add is independence from the C++ implementation, not from a
//! transcription step.
//!
//! Several expectations are tier 4, derived here: `bw_nrd0` of `0..=4` is
//! recomputed in closed form, the Munro packing is checked against the crate's
//! own FFT, and the cubic-spline interpolation inside `kde_fft_eval` is checked
//! against the ported `CubicSpline2d` in `crate::processing`, which is the
//! class the C++ uses and which `src/math/kernel_density.rs` may not depend on.

use openms::Error;
use openms::math::fft::{Complex, real_fft};
use openms::math::kernel_density::{
    DEFAULT_CUT, DEFAULT_GRIDSIZE, MAX_GRID, bw_nrd0, for_rt, grid_kde_fft, kde_fft_eval, lin_bin,
    rev_rt, silverman_kernel_fft,
};
use openms::processing::spline::cubic::CubicSpline2d;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// The class test's `TOLERANCE_ABSOLUTE(1e-6)`.
fn similar(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-6,
        "{actual} is not within 1e-6 of {expected}"
    );
}

// L34 double bwNrd0(const std::vector<double>&)
#[test]
fn section_bw_nrd0() {
    // 100 evenly spaced values: the rule must produce something positive.
    let uniform: Vec<f64> = (0..100).map(|i| i as f64).collect();
    assert!(bw_nrd0(&uniform).unwrap() > 0.0);

    // n = 9, mean 0, sum of squares 15, so the sample variance is 15/8 and the
    // rule is checkable in closed form. The quartiles are the exact grid points
    // -1 and 1, so IQR/1.34 = 1.4925... exceeds sd = sqrt(15/8) = 1.3693... and
    // the rule takes the standard deviation. The class test only brackets the
    // answer in (0.3, 1.0); this pins it.
    let normal = [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0];
    let bw = bw_nrd0(&normal).unwrap();
    assert!(bw > 0.3 && bw < 1.0);
    similar(bw, 0.9 * (15.0_f64 / 8.0).sqrt() * 9.0_f64.powf(-0.2));

    // Constant data: sd and IQR are both zero, so the fallback chain runs.
    assert!(bw_nrd0(&[5.0; 10]).unwrap() > 0.0);

    // L59 all-zero data reaches the last link of the chain, |x[0]| == 0, and
    // then the literal 1.0.
    similar(bw_nrd0(&[0.0; 10]).unwrap(), 0.9 * 10.0_f64.powf(-0.2));

    // Fewer than two finite values yields 0.0, not an error.
    similar(bw_nrd0(&[1.0]).unwrap(), 0.0);
    similar(bw_nrd0(&[]).unwrap(), 0.0);

    // NaN and infinity are filtered before anything is computed.
    let with_nan = [1.0, 2.0, f64::NAN, 3.0, f64::INFINITY, 4.0];
    assert!(bw_nrd0(&with_nan).unwrap() > 0.0);
    // Derived: the survivors are exactly [1, 2, 3, 4], so the answer must equal
    // the bandwidth of that vector.
    similar(
        bw_nrd0(&with_nan).unwrap(),
        bw_nrd0(&[1.0, 2.0, 3.0, 4.0]).unwrap(),
    );
}

// L86 std::vector<double> linBin(const std::vector<double>&, double, double,
//     std::size_t, const std::vector<double>*)
#[test]
fn section_lin_bin() {
    let values = [0.0, 0.5, 1.0];
    let bins = lin_bin(&values, 0.0, 1.0, 3, None).unwrap();
    assert_eq!(bins.len(), 3);
    assert!(bins[0] > 0.0 && bins[1] > 0.0 && bins[2] > 0.0);
    similar(bins.iter().sum::<f64>(), 3.0);
    // Derived: width is 1/3, so each value lands whole in its own bin and the
    // counts are exactly one apiece - which is what shows this is a histogram
    // and not the proportional allocation the header describes.
    assert_eq!(bins, vec![1.0, 1.0, 1.0]);

    let weights = [1.0, 2.0, 3.0];
    let bins = lin_bin(&values, 0.0, 1.0, 3, Some(&weights)).unwrap();
    assert_eq!(bins.len(), 3);
    similar(bins.iter().sum::<f64>(), 6.0);

    // Out-of-range values are ignored; only 0.5 survives.
    let out_of_range = [-1.0, 0.5, 2.0];
    let bins = lin_bin(&out_of_range, 0.0, 1.0, 5, None).unwrap();
    similar(bins.iter().sum::<f64>(), 1.0);

    let bins = lin_bin(&[], 0.0, 1.0, 5, None).unwrap();
    assert_eq!(bins.len(), 5);
    similar(bins.iter().sum::<f64>(), 0.0);

    // MultipleTesting_test.cpp L259 pins the histogram behaviour explicitly.
    let bins = lin_bin(&[0.1, 0.4, 0.9], 0.0, 1.0, 5, None).unwrap();
    assert_eq!(bins, vec![1.0, 0.0, 1.0, 0.0, 1.0]);
    let weights = [2.0, 3.0, 5.0];
    let bins = lin_bin(&[0.1, 0.4, 0.9], 0.0, 1.0, 5, Some(&weights)).unwrap();
    assert_eq!(bins, vec![2.0, 0.0, 3.0, 0.0, 5.0]);

    // A weight slice of the wrong length is ignored, not rejected - the
    // source's `weights->size() == x.size()` test.
    let bins = lin_bin(&[0.1, 0.4, 0.9], 0.0, 1.0, 5, Some(&[7.0])).unwrap();
    assert_eq!(bins, vec![1.0, 0.0, 1.0, 0.0, 1.0]);

    // The two `std::invalid_argument` throws.
    assert!(matches!(
        lin_bin(&[0.5], 0.0, 1.0, 0, None),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        lin_bin(&[0.5], 1.0, 0.0, 4, None),
        Err(Error::InvalidValue(_))
    ));
}

// L131 std::vector<double> forRt(const std::vector<double>&, std::size_t)
#[test]
fn section_for_rt() {
    let signal = [1.0, 2.0, 3.0, 4.0];
    let transformed = for_rt(&signal, 4).unwrap();
    assert_eq!(transformed.len(), 4);
    assert!(transformed[0] > 0.0);
    // Derived: the packing is [Re Y0, Re Y1, Re Y2, Im Y1] and the transform of
    // [1,2,3,4] is [10, -2+2i, -2, -2-2i], so the packed vector is exact.
    assert_eq!(transformed, vec![10.0, -2.0, -2.0, 2.0]);

    // Zero padding to a longer transform.
    let transformed = for_rt(&signal, 8).unwrap();
    assert_eq!(transformed.len(), 8);

    // A constant signal has only a DC component.
    let constant = [5.0; 8];
    let transformed = for_rt(&constant, 8).unwrap();
    assert!(transformed[0] > 0.0);
    similar(transformed[0], 40.0);
    for value in &transformed[1..=4] {
        similar(*value, 0.0);
    }

    // Derived: the packed layout must agree with the crate's own real FFT.
    let values: Vec<f64> = (0..16).map(|i| ((i * 7) % 11) as f64 - 5.0).collect();
    let packed = for_rt(&values, 16).unwrap();
    let spectrum: Vec<Complex> = real_fft(&values).unwrap();
    for k in 0..=8 {
        similar(packed[k], spectrum[k].re);
    }
    for k in 1..8 {
        similar(packed[8 + k], spectrum[k].im);
    }

    // A length of zero yields an empty vector; a non-power-of-two length is
    // refused where the source silently transforms a different size.
    assert!(for_rt(&[], 0).unwrap().is_empty());
    assert!(matches!(for_rt(&signal, 3), Err(Error::InvalidValue(_))));
}

// L157 std::vector<double> revRt(const std::vector<double>&, std::size_t)
#[test]
fn section_rev_rt() {
    let signal = [1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0, 0.0];
    let transformed = for_rt(&signal, 8).unwrap();
    let restored = rev_rt(&transformed, 8).unwrap();
    assert_eq!(restored.len(), signal.len());
    for (got, want) in restored.iter().zip(signal.iter()) {
        similar(*got, *want);
    }

    let small = [1.0, -1.0, 1.0, -1.0];
    let transformed = for_rt(&small, 4).unwrap();
    let restored = rev_rt(&transformed, 4).unwrap();
    for (got, want) in restored.iter().zip(small.iter()) {
        similar(*got, *want);
    }

    // MultipleTesting_test.cpp L206 round-trips [0, 1, 2, 3].
    let x = [0.0, 1.0, 2.0, 3.0];
    let restored = rev_rt(&for_rt(&x, x.len()).unwrap(), x.len()).unwrap();
    for (got, want) in restored.iter().zip(x.iter()) {
        similar(*got, *want);
    }

    // The `std::invalid_argument` for a packed vector shorter than M.
    assert!(matches!(
        rev_rt(&[1.0, 2.0], 4),
        Err(Error::InvalidValue(_))
    ));
    assert!(rev_rt(&[], 0).unwrap().is_empty());
}

// L186 std::vector<double> silvermanKernelFFT(double, std::size_t, double)
#[test]
fn section_silverman_kernel_fft() {
    let m = 8usize;
    let kernel = silverman_kernel_fft(1.0, m, 10.0).unwrap();
    assert_eq!(kernel.len(), m);
    // Bin zero is exp(0) / 1, exactly one - the kernel preserves total mass.
    similar(kernel[0], 1.0);
    assert_eq!(kernel[0], 1.0);
    for value in &kernel[1..=m / 2] {
        assert!(*value <= 1.0);
    }

    let kernel = silverman_kernel_fft(0.5, m, 10.0).unwrap();
    assert_eq!(kernel.len(), m);
    similar(kernel[0], 1.0);

    let kernel = silverman_kernel_fft(1.0, 64, 10.0).unwrap();
    assert_eq!(kernel.len(), 64);

    // MultipleTesting_test.cpp L400: positive, non-increasing over the first
    // half, and mirrored in the packed second half.
    let kernel = silverman_kernel_fft(1.0, 8, 10.0).unwrap();
    let half = 4usize;
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

    // Derived: bin j is exp(-2 (pi bw j / range)^2) / (1 - (j pi / M)^2 / 3).
    let bw = 1.0;
    let range = 10.0;
    for (j, value) in kernel.iter().enumerate().take(half + 1).skip(1) {
        let jf = j as f64;
        let numerator = (-2.0 * (std::f64::consts::PI * bw * jf / range).powi(2)).exp();
        let scaled = jf * std::f64::consts::PI / 8.0;
        let expected = numerator / (1.0 - scaled * scaled / 3.0);
        similar(*value, expected);
    }

    assert!(silverman_kernel_fft(1.0, 0, 10.0).unwrap().is_empty());
    // Native refusals where the source divides by zero.
    assert!(matches!(
        silverman_kernel_fft(1.0, 8, 0.0),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        silverman_kernel_fft(f64::NAN, 8, 10.0),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        silverman_kernel_fft(1.0, MAX_GRID * 2, 10.0),
        Err(Error::InvalidRange(_))
    ));
}

// L213 (gridKdeFFT)
#[test]
fn section_grid_kde_fft() {
    let values: Vec<f64> = (0..20).map(|i| i as f64).collect();
    let result = grid_kde_fft(&values, 1.0, 64, 3.0).unwrap();
    assert_eq!(result.density.len(), result.grid.len());
    // The requested 64 is a lower bound; the source raises it to 512.
    assert!(result.density.len() >= 64);
    assert_eq!(result.density.len(), 512);
    for pair in result.grid.windows(2) {
        assert!(pair[1] > pair[0]);
    }
    for value in &result.density {
        assert!(*value >= 0.0);
    }
    let spacing =
        (result.grid[result.grid.len() - 1] - result.grid[0]) / (result.density.len() - 1) as f64;
    let integral: f64 = result.density.iter().map(|v| v * spacing).sum();
    assert!(integral > 0.5 && integral < 2.0);
    // Derived: the source renormalises so sum * spacing is exactly one.
    similar(integral, 1.0);

    // A single point gives a Gaussian-shaped bump peaking near it.
    let result = grid_kde_fft(&[0.0], 1.0, 32, 3.0).unwrap();
    assert!(result.density.len() >= 32);
    assert!(result.grid.len() >= 32);
    let peak = result
        .density
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(index, _)| index)
        .unwrap();
    assert!(result.grid[peak].abs() < 1.0);

    let bimodal = [-5.0, -4.9, -5.1, 5.0, 4.9, 5.1];
    let result = grid_kde_fft(&bimodal, 0.5, 128, 3.0).unwrap();
    assert!(result.density.len() >= 128);

    // Empty data still produces a grid, centred on zero.
    let result = grid_kde_fft(&[], 1.0, 32, 3.0).unwrap();
    assert!(result.density.len() >= 32);
    similar(result.grid[0], -3.0);
    similar(result.grid[result.grid.len() - 1], 3.0);

    // MultipleTesting_test.cpp L425: the estimate integrates to one.
    let x = [0.0, 1.0, 2.0];
    let bw = bw_nrd0(&x).unwrap();
    let result = grid_kde_fft(&x, bw, 64, 3.0).unwrap();
    let delta = result.grid[1] - result.grid[0];
    let mut sum = 0.0;
    for value in &result.density {
        assert!(value.is_finite());
        assert!(*value >= 0.0);
        sum += *value;
    }
    similar(sum * delta, 1.0);

    // Native refusals.
    assert!(matches!(
        grid_kde_fft(&[0.0, f64::NAN], 1.0, 64, 3.0),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        grid_kde_fft(&[1.0, 1.0], 0.0, 64, 3.0),
        Err(Error::InvalidValue(_))
    ));
}

// L268 (kdeFFTEval)
#[test]
fn section_kde_fft_eval() {
    let values = [0.0, 1.0, 2.0, 3.0, 4.0];
    let densities = kde_fft_eval(&values, 1.0, 64, 3.0).unwrap();
    assert_eq!(densities.len(), values.len());
    for value in &densities {
        assert!(*value > 0.0);
        assert!(value.is_finite());
    }
    let mean_density = densities.iter().sum::<f64>() / densities.len() as f64;
    for value in &densities {
        assert!(*value > mean_density * 0.5 && *value < mean_density * 1.5);
    }

    let mut clustered: Vec<f64> = (0..50).map(|i| i as f64 * 0.01).collect();
    clustered.extend((0..50).map(|i| 10.0 + i as f64 * 0.01));
    let densities = kde_fft_eval(&clustered, 0.5, 128, 3.0).unwrap();
    assert_eq!(densities.len(), clustered.len());
    let first: f64 = densities[..50].iter().sum::<f64>() / 50.0;
    let second: f64 = densities[50..].iter().sum::<f64>() / 50.0;
    assert!(first > 0.0);
    assert!(second > 0.0);

    let densities = kde_fft_eval(&[5.0], 1.0, 64, 3.0).unwrap();
    assert_eq!(densities.len(), 1);
    assert!(densities[0] > 0.0);

    // Identical points all receive the same density.
    let identical = [3.0; 10];
    let densities = kde_fft_eval(&identical, 0.5, 64, 3.0).unwrap();
    assert_eq!(densities.len(), 10);
    for value in &densities {
        assert!(*value > 0.0);
        similar(*value, densities[0]);
    }

    assert!(kde_fft_eval(&[], 1.0, 64, 3.0).unwrap().is_empty());
}

// MultipleTesting_test.cpp L446 "kdeFFTEval matches spline-interpolated grid
// density at sample points". Tier 4 here: the expectation is computed with the
// crate's own port of `CubicSpline2d`, the very class the C++ uses, which
// `src/math/kernel_density.rs` cannot reach across the module graph and so
// reimplements privately. This is what proves the reimplementation agrees.
#[test]
fn section_kde_eval_matches_the_ported_cubic_spline() {
    let x = [0.0, 1.0, 2.0];
    let bw = bw_nrd0(&x).unwrap();
    let grid = grid_kde_fft(&x, bw, 64, 3.0).unwrap();
    let spline = CubicSpline2d::new(&grid.grid, &grid.density).unwrap();
    let expected: Vec<f64> = x.iter().map(|v| spline.eval(*v).unwrap()).collect();
    let got = kde_fft_eval(&x, bw, 64, 3.0).unwrap();
    assert_eq!(got.len(), expected.len());
    for (index, (a, b)) in got.iter().zip(expected.iter()).enumerate() {
        assert_eq!(a, b, "sample {index}: {a} != {b}");
    }

    // And on a larger, less regular sample, where a coefficient recurrence that
    // had drifted would show up.
    let sample: Vec<f64> = (0..64).map(|i| ((i * 37) % 101) as f64 * 0.1).collect();
    let bw = bw_nrd0(&sample).unwrap();
    let grid = grid_kde_fft(&sample, bw, DEFAULT_GRIDSIZE, DEFAULT_CUT).unwrap();
    let spline = CubicSpline2d::new(&grid.grid, &grid.density).unwrap();
    let got = kde_fft_eval(&sample, bw, DEFAULT_GRIDSIZE, DEFAULT_CUT).unwrap();
    for (index, value) in sample.iter().enumerate() {
        let want = spline.eval(*value).unwrap().max(0.0);
        assert_eq!(got[index], want, "sample {index}");
    }
}

// L339 (Integration test: Full KDE pipeline)
#[test]
fn section_integration_pipeline() {
    let mut values = vec![-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0];
    values.extend([4.0, 4.5, 5.0, 5.5, 6.0]);

    let bw = bw_nrd0(&values).unwrap();
    assert!(bw > 0.0);

    let grid = grid_kde_fft(&values, bw, 256, 3.0).unwrap();
    assert!(grid.density.len() >= 256);
    assert!(grid.grid.len() >= 256);

    let densities = kde_fft_eval(&values, bw, 256, 3.0).unwrap();
    assert_eq!(densities.len(), values.len());
    for value in &densities {
        assert!(*value > 0.0);
        assert!(value.is_finite());
    }
    let max = densities.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let min = densities.iter().copied().fold(f64::INFINITY, f64::min);
    assert!(max > min);
}

// L381 "kde reference (statsmodels) regression"
#[test]
fn section_statsmodels_reference() {
    let text = std::fs::read_to_string(data("kde_reference_data.csv")).unwrap();
    let mut datasets: BTreeMap<String, (Vec<f64>, Vec<f64>, f64)> = BTreeMap::new();
    for line in text.lines().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split(',').map(|p| p.trim_matches('"')).collect();
        if parts.len() < 5 {
            continue;
        }
        let entry = datasets.entry(parts[0].to_string()).or_default();
        entry.0.push(parts[2].parse().unwrap());
        entry.1.push(parts[4].parse().unwrap());
        entry.2 = parts[3].parse().unwrap();
    }
    assert!(!datasets.is_empty());

    // The class test's tolerances.
    let bw_rel_tol = 0.05;
    let dens_rel_tol = 0.05;
    let dens_abs_tol = 0.005;
    for (name, (points, reference, reference_bw)) in &datasets {
        let our_bw = bw_nrd0(points).unwrap();
        assert!(our_bw > 0.0, "{name}");
        let bw_error = if *reference_bw > 1e-10 {
            (our_bw - reference_bw).abs() / reference_bw
        } else {
            (our_bw - reference_bw).abs()
        };
        assert!(
            bw_error <= bw_rel_tol,
            "{name}: bandwidth {our_bw} vs {reference_bw}"
        );

        let ours = kde_fft_eval(points, our_bw, DEFAULT_GRIDSIZE, DEFAULT_CUT).unwrap();
        assert_eq!(ours.len(), reference.len(), "{name}");
        for (index, (got, want)) in ours.iter().zip(reference.iter()).enumerate() {
            let relative = if *want > 1e-10 {
                (got - want).abs() / want.abs()
            } else {
                (got - want).abs()
            };
            assert!(
                relative <= dens_rel_tol || (got - want).abs() <= dens_abs_tol,
                "{name}[{index}]: {got} vs {want}"
            );
        }
    }
}

// L466 "kde density vs scipy.stats.gaussian_kde reference"
#[test]
fn section_scipy_reference() {
    let values = [
        0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 2.5, 3.5, 3.0, 4.5,
    ];
    let bw = 1.2298580108;
    let scipy = [
        0.05184986, 0.08447517, 0.12205155, 0.14815310, 0.14174787, 0.11287841, 0.08684046,
        0.07360770, 0.06438698, 0.04732923, 0.13832075, 0.14919717, 0.14815310, 0.12843312,
    ];

    let our_bw = bw_nrd0(&values).unwrap();
    assert!((our_bw - bw).abs() / bw <= 0.05, "bandwidth {our_bw}");

    let ours = kde_fft_eval(&values, bw, DEFAULT_GRIDSIZE, DEFAULT_CUT).unwrap();
    assert_eq!(ours.len(), scipy.len());
    for (index, (got, want)) in ours.iter().zip(scipy.iter()).enumerate() {
        let relative = (got - want).abs() / want.abs();
        assert!(
            relative <= 0.10 || (got - want).abs() <= 0.01,
            "sample {index}: {got} vs {want}"
        );
    }
}
