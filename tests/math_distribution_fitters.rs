// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Class-test parity for the four `MATH/STATISTICS` distribution fitters.
//!
//! Every `START_SECTION` of `GaussFitter_test.cpp`,
//! `GammaDistributionFitter_test.cpp`, `GumbelDistributionFitter_test.cpp` and
//! `GumbelMaxLikelihoodFitter_test.cpp` at revision `bc9cc12` has a test here.
//! Literals transcribed from those files are tier-3 evidence (source review);
//! the tests that derive a value from a closed form or an invariant say so and
//! are tier 4. See `docs/DISTRIBUTION_FITTERS_SUPPORT.md` and
//! `tests/data/distribution_fitters_provenance.json`.

// The expected values are transcribed character for character from the class
// tests, several of which carry more decimal digits than `f64` can hold.
// Shortening them would break the correspondence a reviewer checks.
#![allow(clippy::excessive_precision)]

use openms::math::fitters::gamma::{GammaDistributionFitResult, GammaDistributionFitter};
use openms::math::fitters::gauss::{GaussFitResult, GaussFitter};
use openms::math::fitters::gumbel::{GumbelDistributionFitResult, GumbelDistributionFitter};
use openms::math::fitters::gumbel_max_likelihood::{
    GumbelDistributionFitResult as MleFitResult, GumbelMaxLikelihoodFitter,
};

fn close(got: f64, want: f64, tolerance: f64) {
    assert!(
        (got - want).abs() <= tolerance,
        "got {got:?}, want {want:?}, difference {:.3e} exceeds {tolerance:.3e}",
        (got - want).abs()
    );
}

fn close_relative(got: f64, want: f64, tolerance: f64) {
    let relative = (got - want).abs() / want.abs();
    assert!(
        relative <= tolerance,
        "got {got:?}, want {want:?}, relative difference {relative:.3e} exceeds {tolerance:.3e}"
    );
}

// ---------------------------------------------------------------------------
// GaussFitter_test.cpp
// ---------------------------------------------------------------------------

/// The m/z grid of the second `fit` case, `GaussFitter_test.cpp:41-48`.
const GAUSS_MZ: [f64; 7] = [
    240.1000470172,
    240.1002675493,
    240.1004880817,
    240.1007086145,
    240.1009291475,
    240.1011496808,
    240.1013702145,
];

/// The intensities of the second `fit` case, `GaussFitter_test.cpp:50-58`.
const GAUSS_INTENSITIES: [f64; 7] = [
    61134.39453125,
    111288.5390625,
    163761.46875,
    165861.4375,
    162133.46875,
    120060.5234375,
    71102.1328125,
];

/// The initial guess of the second `fit` case, `GaussFitter_test.cpp:61-65`.
fn gauss_peak_guess() -> GaussFitResult {
    GaussFitResult::new(168324.0, 240.10051, 0.000375375)
}

/// `START_SECTION(GaussFitter())`
#[test]
fn gauss_fitter_default_construction() {
    let fitter = GaussFitter::new();
    // The C++ only asserts the pointer is non-null; the observable state of a
    // default-constructed fitter is its initial guess, `GaussFitter.cpp:23`.
    assert_eq!(
        fitter.initial_parameters(),
        GaussFitResult::new(0.06, 3.0, 0.5)
    );
    assert_eq!(fitter, GaussFitter::default());
}

/// `START_SECTION((virtual ~GaussFitter()))`
#[test]
fn gauss_fitter_owns_no_resources_to_destroy() {
    // The C++ section is NOT_TESTABLE: it only deletes the pointer. The Rust
    // counterpart of "the destructor is trivial" is that the type is a plain
    // `Copy` value of exactly its three parameters, with no `Drop`.
    assert_eq!(
        std::mem::size_of::<GaussFitter>(),
        3 * std::mem::size_of::<f64>()
    );
    let fitter = GaussFitter::new();
    let copy = fitter;
    assert_eq!(copy.initial_parameters().sigma, 0.5);
}

/// `START_SECTION((GaussFitResult fit(std::vector< DPosition< 2 > >& points) const))`
#[test]
fn gauss_fit_reproduces_both_published_cases() {
    // Case one: the default initial guess (0.06, 3.0, 0.5) and six points,
    // `GaussFitter_test.cpp:70-97`.
    let points = [
        (0.0, 0.01),
        (0.05, 0.2),
        (0.16, 0.63),
        (0.28, 0.99),
        (0.66, 0.03),
        (0.50, 0.36),
    ];
    let result = GaussFitter::new().fit(&points).unwrap();
    close_relative(result.a, 1.01898275662372, 1e-9);
    close_relative(result.x0, 0.300612870901173, 1e-9);
    close_relative(result.sigma, 0.136316330927453, 1e-9);

    // Case two: the source's own comment says this one settles on a negative
    // sigma internally and needs the `fabs`, `GaussFitter_test.cpp:99-119`.
    let peak: Vec<(f64, f64)> = GAUSS_MZ
        .iter()
        .copied()
        .zip(GAUSS_INTENSITIES.iter().copied())
        .collect();
    let mut fitter = GaussFitter::new();
    fitter.set_initial_parameters(gauss_peak_guess());
    let fitted = fitter.fit(&peak).unwrap();
    close_relative(fitted.a, 175011.893006749, 1e-11);
    close_relative(fitted.x0, 240.1007246725147, 1e-11);
    close_relative(fitted.sigma, 0.00046642320683761701, 1e-11);
    // The reported width is the absolute value, so it is usable as a width.
    assert!(fitted.sigma > 0.0);
}

/// `START_SECTION((void setInitialParameters(const GaussFitResult& result)))`
#[test]
fn gauss_set_initial_parameters_is_read_back_and_used() {
    // The C++ section is NOT_TESTABLE and says the setter is implicitly tested
    // by `fit`. Here it is tested directly and through `fit`: the same points
    // fitted from two different guesses reach two different results.
    let mut fitter = GaussFitter::new();
    fitter.set_initial_parameters(GaussFitResult::new(-1.0, -1.0, -1.0));
    assert_eq!(
        fitter.initial_parameters(),
        GaussFitResult::new(-1.0, -1.0, -1.0)
    );

    let peak: Vec<(f64, f64)> = GAUSS_MZ
        .iter()
        .copied()
        .zip(GAUSS_INTENSITIES.iter().copied())
        .collect();
    let mut guided = GaussFitter::new();
    guided.set_initial_parameters(gauss_peak_guess());
    let from_guess = guided.fit(&peak).unwrap();
    let from_default = GaussFitter::new().fit(&peak).unwrap();
    close_relative(from_guess.x0, 240.1007246725147, 1e-11);
    assert!(
        (from_default.x0 - from_guess.x0).abs() > 1.0,
        "the default guess must not land on the guided result: {from_default:?}"
    );
}

/// `START_SECTION((static std::vector<double> eval(const std::vector<double>& evaluation_points, const GaussFitResult& model)))`
#[test]
fn gauss_static_eval_reproduces_the_published_intensities() {
    // `GaussFitter_test.cpp:134-150`, evaluated at the m/z grid with the
    // initial guess as the model.
    let expected = [
        78670.515322697669,
        136633.77791868619,
        168037.29915800504,
        146337.00743127937,
        90240.802825824489,
        39405.008909696895,
        12184.248044493703,
    ];
    let got = GaussFitter::eval(&GAUSS_MZ, &gauss_peak_guess()).unwrap();
    assert_eq!(got.len(), expected.len());
    for (&got, &want) in got.iter().zip(expected.iter()) {
        close_relative(got, want, 1e-14);
    }
    // The member form must agree with the static one.
    let model = gauss_peak_guess();
    for (&point, &want) in GAUSS_MZ.iter().zip(expected.iter()) {
        close_relative(model.eval(point).unwrap(), want, 1e-14);
    }
}

/// Derived rather than transcribed: the log density is the log of the density
/// divided by the amplitude, which is what "no normalize" means.
#[test]
fn gauss_log_eval_is_the_log_of_the_unit_amplitude_density() {
    let model = GaussFitResult::new(3.0, 1.5, 0.4);
    let unit = GaussFitResult::new(1.0, 1.5, 0.4);
    for &x in &[0.5, 1.5, 2.7] {
        let density = unit.eval(x).unwrap() / (0.4 * (2.0 * std::f64::consts::PI).sqrt());
        close_relative(model.log_eval_no_normalize(x).unwrap(), density.ln(), 1e-13);
    }
}

// ---------------------------------------------------------------------------
// GammaDistributionFitter_test.cpp
// ---------------------------------------------------------------------------

/// The 40 points of the `fit` section, `GammaDistributionFitter_test.cpp:45-84`.
const GAMMA_POINTS: [(f64, f64); 40] = [
    (0.0001, 0.1),
    (0.0251, 0.3),
    (0.0501, 0.0),
    (0.0751, 0.7),
    (0.1001, 0.0),
    (0.1251, 1.6),
    (0.1501, 0.0),
    (0.1751, 2.1),
    (0.2001, 0.0),
    (0.2251, 3.7),
    (0.2501, 0.0),
    (0.2751, 4.0),
    (0.3001, 0.0),
    (0.3251, 3.0),
    (0.3501, 0.0),
    (0.3751, 2.6),
    (0.4001, 0.0),
    (0.4251, 3.0),
    (0.4501, 0.0),
    (0.4751, 3.0),
    (0.5001, 0.0),
    (0.5251, 2.5),
    (0.5501, 0.0),
    (0.5751, 1.7),
    (0.6001, 0.0),
    (0.6251, 1.0),
    (0.6501, 0.0),
    (0.6751, 0.5),
    (0.7001, 0.0),
    (0.7251, 0.3),
    (0.7501, 0.0),
    (0.7751, 0.4),
    (0.8001, 0.0),
    (0.8251, 0.0),
    (0.8501, 0.0),
    (0.8751, 0.1),
    (0.9001, 0.0),
    (0.9251, 0.1),
    (0.9501, 0.0),
    (0.9751, 0.2),
];

/// `START_SECTION(GammaDistributionFitter())`
#[test]
fn gamma_fitter_default_construction() {
    let fitter = GammaDistributionFitter::new();
    // `GammaDistributionFitter.cpp:27`.
    assert_eq!(
        fitter.initial_parameters(),
        GammaDistributionFitResult::new(1.0, 5.0)
    );
    assert_eq!(fitter, GammaDistributionFitter::default());
}

/// `START_SECTION(virtual ~GammaDistributionFitter())`
#[test]
fn gamma_fitter_owns_no_resources_to_destroy() {
    // The C++ section only deletes the pointer.
    assert_eq!(
        std::mem::size_of::<GammaDistributionFitter>(),
        2 * std::mem::size_of::<f64>()
    );
    let fitter = GammaDistributionFitter::new();
    let copy = fitter;
    assert_eq!(copy.initial_parameters().p, 5.0);
}

/// `START_SECTION((GammaDistributionFitResult fit(std::vector< DPosition< 2 > > & points)))`
#[test]
fn gamma_fit_reproduces_the_published_parameters() {
    // `GammaDistributionFitter_test.cpp:86-94`: initial guess (1.0, 3.0),
    // expected b = 7.25 and p = 3.11 at TOLERANCE_ABSOLUTE(0.01).
    let mut fitter = GammaDistributionFitter::new();
    fitter.set_initial_parameters(GammaDistributionFitResult::new(1.0, 3.0));
    let result = fitter.fit(&GAMMA_POINTS).unwrap();
    close(result.b, 7.25, 0.01);
    close(result.p, 3.11, 0.01);
    // Both parameters stay in the domain where the customized density is the
    // real Gamma density, so the fit is a distribution and not the zero branch.
    assert!(result.b > 0.0 && result.p > 0.0);
}

/// `START_SECTION((void setInitialParameters(const GammaDistributionFitResult & result)))`
#[test]
fn gamma_set_initial_parameters_is_read_back_and_used() {
    let mut fitter = GammaDistributionFitter::new();
    fitter.set_initial_parameters(GammaDistributionFitResult::new(1.0, 5.0));
    assert_eq!(
        fitter.initial_parameters(),
        GammaDistributionFitResult::new(1.0, 5.0)
    );
    // The C++ calls this implicitly tested by `fit`: the published result is
    // reached from the published guess, and the default guess is a different
    // starting point on the same surface.
    let from_default = GammaDistributionFitter::new().fit(&GAMMA_POINTS).unwrap();
    let mut guided = GammaDistributionFitter::new();
    guided.set_initial_parameters(GammaDistributionFitResult::new(1.0, 3.0));
    let from_guess = guided.fit(&GAMMA_POINTS).unwrap();
    close(from_guess.b, 7.25, 0.01);
    assert!(from_default.b.is_finite());
}

// ---------------------------------------------------------------------------
// GumbelDistributionFitter_test.cpp
// ---------------------------------------------------------------------------

/// The 30 points of the first `fit` case, `GumbelDistributionFitter_test.cpp:42-71`.
const GUMBEL_POINTS: [(f64, f64); 30] = [
    (-2.7, 0.017),
    (-2.5, 0.025),
    (-2.0, 0.052),
    (-1.0, 0.127),
    (-0.7, 0.147),
    (-0.01, 0.178),
    (0.0, 0.178),
    (0.2, 0.182),
    (0.5, 0.184),
    (1.0, 0.179),
    (1.3, 0.171),
    (1.9, 0.151),
    (2.5, 0.127),
    (2.6, 0.123),
    (2.7, 0.119),
    (2.8, 0.115),
    (2.9, 0.111),
    (3.0, 0.108),
    (3.5, 0.089),
    (3.9, 0.076),
    (4.01, 0.073),
    (4.22, 0.067),
    (4.7, 0.054),
    (4.9, 0.05),
    (5.0, 0.047),
    (6.0, 0.03),
    (7.0, 0.017),
    (7.5, 0.015),
    (7.9, 0.012),
    (8.03, 0.011),
];

/// The 10 points of the second `fit` case, `GumbelDistributionFitter_test.cpp:87-96`.
const GUMBEL_POINTS_2: [(f64, f64); 10] = [
    (0.0, 0.18),
    (0.2, 0.24),
    (0.5, 0.32),
    (1.0, 0.37),
    (1.3, 0.35),
    (1.9, 0.27),
    (2.5, 0.18),
    (2.6, 0.16),
    (3.0, 0.12),
    (5.0, 0.02),
];

/// `START_SECTION(GumbelDistributionFitter())`
#[test]
fn gumbel_fitter_default_construction() {
    // `GumbelDistributionFitter.cpp:32`.
    assert_eq!(
        GumbelDistributionFitter::new().initial_parameters(),
        GumbelDistributionFitResult::new(0.25, 0.1)
    );
}

/// `START_SECTION((virtual ~GumbelDistributionFitter()))`
#[test]
fn gumbel_fitter_owns_no_resources_to_destroy() {
    assert_eq!(
        std::mem::size_of::<GumbelDistributionFitter>(),
        2 * std::mem::size_of::<f64>()
    );
    let fitter = GumbelDistributionFitter::new();
    let copy = fitter;
    assert_eq!(copy.initial_parameters().b, 0.1);
}

/// `START_SECTION((GumbelDistributionFitResult fit(std::vector<DPosition<2> >& points)))`
#[test]
fn gumbel_fit_reproduces_both_published_cases() {
    // `GumbelDistributionFitter_test.cpp:76-84`: guess (1.0, 3.0), the data
    // were generated from a = 0.5, b = 2.0, TOLERANCE_ABSOLUTE(0.1).
    let mut fitter = GumbelDistributionFitter::new();
    fitter.set_initial_parameters(GumbelDistributionFitResult::new(1.0, 3.0));
    let result = fitter.fit(&GUMBEL_POINTS).unwrap();
    close(result.a, 0.5, 0.1);
    close(result.b, 2.0, 0.1);

    // `GumbelDistributionFitter_test.cpp:99-106`: guess (3.0, 3.0), generating
    // parameters a = 1.0, b = 1.0.
    let mut second = GumbelDistributionFitter::new();
    second.set_initial_parameters(GumbelDistributionFitResult::new(3.0, 3.0));
    let result = second.fit(&GUMBEL_POINTS_2).unwrap();
    close(result.a, 1.0, 0.1);
    close(result.b, 1.0, 0.1);
}

/// `START_SECTION((void setInitialParameters(const GumbelDistributionFitResult& result)))`
#[test]
fn gumbel_set_initial_parameters_is_read_back_and_used() {
    let mut fitter = GumbelDistributionFitter::new();
    // The C++ section passes a default-constructed result, which is (1.0, 2.0).
    fitter.set_initial_parameters(GumbelDistributionFitResult::default());
    assert_eq!(
        fitter.initial_parameters(),
        GumbelDistributionFitResult::new(1.0, 2.0)
    );
    let from_default_result = fitter.fit(&GUMBEL_POINTS).unwrap();
    close(from_default_result.a, 0.5, 0.1);
}

/// `START_SECTION((GumbelDistributionFitter(const GumbelDistributionFitter& rhs)))`
#[test]
fn gumbel_fitter_copy_carries_the_initial_guess() {
    // The C++ declares the copy constructor private and never defines it, so
    // the section is NOT_TESTABLE there. Rust has no reason to forbid copying a
    // two-field value, so the port derives `Copy` and this asserts the copy is
    // a real one.
    let mut fitter = GumbelDistributionFitter::new();
    fitter.set_initial_parameters(GumbelDistributionFitResult::new(5.0, 4.0));
    let copy = fitter;
    assert_eq!(
        copy.initial_parameters(),
        GumbelDistributionFitResult::new(5.0, 4.0)
    );
}

/// `START_SECTION((GumbelDistributionFitter& operator = (const GumbelDistributionFitter& rhs)))`
#[test]
fn gumbel_fitter_assignment_replaces_the_initial_guess() {
    // Also NOT_TESTABLE in the C++, for the same reason.
    let mut target = GumbelDistributionFitter::new();
    assert_eq!(
        target.initial_parameters(),
        GumbelDistributionFitResult::new(0.25, 0.1)
    );
    let mut source = GumbelDistributionFitter::new();
    source.set_initial_parameters(GumbelDistributionFitResult::new(3.0, 2.2));
    target = source;
    assert_eq!(
        target.initial_parameters(),
        GumbelDistributionFitResult::new(3.0, 2.2)
    );
}

/// `START_SECTION((GumbelDistributionFitter::GumbelDistributionFitResult()))`
#[test]
fn gumbel_result_default_construction() {
    // `GumbelDistributionFitter_test.cpp:130-134` asserts a = 1.0, b = 2.0.
    let result = GumbelDistributionFitResult::default();
    close(result.a, 1.0, 1e-12);
    close(result.b, 2.0, 1e-12);
}

/// `START_SECTION((GumbelDistributionFitter::GumbelDistributionFitResult(const GumbelDistributionFitter::GumbelDistributionFitResult& rhs)))`
#[test]
fn gumbel_result_copy_construction() {
    // `GumbelDistributionFitter_test.cpp:137-142`.
    let original = GumbelDistributionFitResult::new(5.0, 4.0);
    let copy = original;
    close(copy.a, 5.0, 1e-12);
    close(copy.b, 4.0, 1e-12);
}

/// `START_SECTION((GumbelDistributionFitter::GumbelDistributionFitResult& operator = (const GumbelDistributionFitter::GumbelDistributionFitResult& rhs)))`
#[test]
fn gumbel_result_assignment() {
    // `GumbelDistributionFitter_test.cpp:145-151`.
    let source = GumbelDistributionFitResult::new(3.0, 2.2);
    let mut target = GumbelDistributionFitResult::default();
    close(target.a, 1.0, 1e-12);
    target = source;
    close(target.a, 3.0, 1e-12);
    close(target.b, 2.2, 1e-12);
}

/// `START_SECTION(MLE)`
#[test]
fn gumbel_maximum_likelihood_on_the_published_sample() {
    // `GumbelDistributionFitter_test.cpp:154-178`: 1200 samples from
    // `Gumbel_1D.csv`, unit weights, guess (4.0, 2.0). The section inherits
    // TOLERANCE_ABSOLUTE(0.1) from the `fit` section above it.
    let samples = read_gumbel_samples();
    assert_eq!(samples.len(), 1200);
    let weights = vec![1.0; samples.len()];
    let mut fitter =
        GumbelMaxLikelihoodFitter::with_initial_parameters(MleFitResult::new(4.0, 2.0));
    let result = fitter.fit_weighted(&samples, &weights).unwrap();
    close(result.a, 2.0, 0.1);
    close(result.b, 0.6, 0.1);
    // The fit is written back into the fitter, as the source does.
    assert_eq!(fitter.initial_parameters(), result);
}

fn read_gumbel_samples() -> Vec<f64> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/gumbel_1d.csv");
    let text = std::fs::read_to_string(path).expect("gumbel_1d.csv");
    text.split([',', '\n', '\r'])
        .map(str::trim)
        .filter(|token| !token.is_empty())
        .map(|token| token.parse::<f64>().expect("sample"))
        .collect()
}

// ---------------------------------------------------------------------------
// GumbelMaxLikelihoodFitter_test.cpp
// ---------------------------------------------------------------------------

/// `START_SECTION((GumbelMaxLikelihoodFitter()))`
#[test]
fn mle_fitter_default_construction() {
    // `GumbelMaxLikelihoodFitter_test.cpp:33-41`: with no data points,
    // `fitWeighted` returns the default initial parameters.
    let mut fitter = GumbelMaxLikelihoodFitter::new();
    let result = fitter.fit_weighted(&[], &[]).unwrap();
    close(result.a, 0.25, 1e-12);
    close(result.b, 0.1, 1e-12);
}

/// `START_SECTION((~GumbelMaxLikelihoodFitter()))`
#[test]
fn mle_fitter_owns_no_resources_to_destroy() {
    assert_eq!(
        std::mem::size_of::<GumbelMaxLikelihoodFitter>(),
        2 * std::mem::size_of::<f64>()
    );
    let fitter = GumbelMaxLikelihoodFitter::new();
    let copy = fitter;
    assert_eq!(copy.initial_parameters().a, 0.25);
}

/// `START_SECTION((GumbelMaxLikelihoodFitter(GumbelDistributionFitResult init)))`
#[test]
fn mle_fitter_construction_from_initial_parameters() {
    // `GumbelMaxLikelihoodFitter_test.cpp:50-57`.
    let mut fitter =
        GumbelMaxLikelihoodFitter::with_initial_parameters(MleFitResult::new(3.0, 0.5));
    let result = fitter.fit_weighted(&[], &[]).unwrap();
    close(result.a, 3.0, 1e-12);
    close(result.b, 0.5, 1e-12);
}

/// `START_SECTION((void setInitialParameters(const GumbelDistributionFitResult & result)))`
#[test]
fn mle_set_initial_parameters() {
    // `GumbelMaxLikelihoodFitter_test.cpp:60-67`.
    let mut fitter = GumbelMaxLikelihoodFitter::new();
    fitter.set_initial_parameters(MleFitResult::new(4.0, 0.7));
    let result = fitter.fit_weighted(&[], &[]).unwrap();
    close(result.a, 4.0, 1e-12);
    close(result.b, 0.7, 1e-12);
}

/// `START_SECTION((GumbelDistributionFitResult fitWeighted(const std::vector<double> & x, const std::vector<double> & w)))`
#[test]
fn mle_fit_weighted_recovers_the_generating_parameters() {
    // `GumbelMaxLikelihoodFitter_test.cpp:70-91`: a Gumbel(2.0, 0.8) density
    // sampled on a 0.1 grid over [-2, 8] and used as its own weights, guess
    // (1.0, 1.0), TOLERANCE_ABSOLUTE(0.05).
    const A: f64 = 2.0;
    const B: f64 = 0.8;
    let mut x = Vec::new();
    let mut w = Vec::new();
    let mut xi = -2.0f64;
    while xi <= 8.0 {
        let z = (xi - A) / B;
        w.push((1.0 / B) * (-(z + (-z).exp())).exp());
        x.push(xi);
        xi += 0.1;
    }
    let mut fitter =
        GumbelMaxLikelihoodFitter::with_initial_parameters(MleFitResult::new(1.0, 1.0));
    let result = fitter.fit_weighted(&x, &w).unwrap();
    close(result.a, A, 0.05);
    close(result.b, B, 0.05);
    assert!(result.b > 0.0);
}

/// `START_SECTION(([EXTRA] GumbelDistributionFitResult(double a, double b) + log_eval_no_normalize))`
#[test]
fn mle_result_stores_its_parameters_and_evaluates_the_log_density() {
    // `GumbelMaxLikelihoodFitter_test.cpp:94-117`. The expected values are
    // derived from the closed form rather than transcribed.
    let log_gumbel = |x: f64, a: f64, b: f64| {
        let diff = (x - a) / b;
        -b.ln() - diff - (-diff).exp()
    };
    let result = MleFitResult::new(2.0, 1.0);
    close(result.a, 2.0, 1e-12);
    close(result.b, 1.0, 1e-12);
    for &x in &[2.0, 3.0, 0.0] {
        close(
            result.log_eval_no_normalize(x).unwrap(),
            log_gumbel(x, 2.0, 1.0),
            1e-13,
        );
    }
    // The log density at the mode is exactly -ln(b) - 1.
    close(result.log_eval_no_normalize(2.0).unwrap(), -1.0, 1e-13);
    assert!(
        result.log_eval_no_normalize(2.0).unwrap() > result.log_eval_no_normalize(3.0).unwrap()
    );
    assert!(
        result.log_eval_no_normalize(2.0).unwrap() > result.log_eval_no_normalize(1.0).unwrap()
    );
    let wider = MleFitResult::new(2.0, 2.0);
    close(
        wider.log_eval_no_normalize(2.0).unwrap(),
        log_gumbel(2.0, 2.0, 2.0),
        1e-13,
    );
}

// ---------------------------------------------------------------------------
// Native boundaries, beyond the class tests
// ---------------------------------------------------------------------------

/// Every fitter refuses input it cannot use rather than returning NaN.
#[test]
fn fitters_refuse_degenerate_input() {
    assert!(GaussFitter::new().fit(&[(0.0, 1.0), (1.0, 2.0)]).is_err());
    assert!(GammaDistributionFitter::new().fit(&[(1.0, 1.0)]).is_err());
    assert!(GumbelDistributionFitter::new().fit(&[(1.0, 1.0)]).is_err());
    let mut mle = GumbelMaxLikelihoodFitter::new();
    assert!(mle.fit_weighted(&[1.0, 2.0], &[1.0]).is_err());
    assert!(
        GaussFitter::new()
            .fit(&[(0.0, 0.0), (1.0, f64::INFINITY), (2.0, 0.0)])
            .is_err()
    );
    assert!(
        GumbelDistributionFitter::new()
            .fit(&[(0.0, 0.0), (f64::NAN, 1.0)])
            .is_err()
    );
    assert!(
        GammaDistributionFitter::new()
            .fit(&[(-1.0, 1.0), (1.0, 1.0)])
            .is_err()
    );
}

/// Independent of every transcribed literal: the parameters the
/// maximum-likelihood fitter reports are a local minimum of the weighted
/// negative log-likelihood the source defines, not merely the point at which
/// its solver happened to stop.
///
/// This matters because the source hands Levenberg-Marquardt a residual vector
/// whose only non-zero entry is that scalar, so the solver minimizes its
/// square. Where the objective stays strictly positive - as it does on all the
/// class-test data - the square has its stationary points exactly where the
/// objective does, and the fit is the maximum-likelihood estimate.
#[test]
fn mle_result_is_a_local_minimum_of_the_weighted_log_likelihood() {
    let samples = read_gumbel_samples();
    let weights = vec![1.0; samples.len()];
    let negative_log_likelihood = |a: f64, b: f64| {
        let sigma = b.abs();
        let log_sigma = sigma.ln();
        let mut sum = 0.0;
        for (&x, &w) in samples.iter().zip(weights.iter()) {
            let diff = (x - a) / sigma;
            sum += w * (-log_sigma - diff - (-diff).exp());
        }
        -sum
    };
    let mut fitter =
        GumbelMaxLikelihoodFitter::with_initial_parameters(MleFitResult::new(4.0, 2.0));
    let result = fitter.fit_weighted(&samples, &weights).unwrap();
    let at_optimum = negative_log_likelihood(result.a, result.b);
    // The objective is strictly positive here, which is the condition under
    // which minimizing its square is minimizing it.
    assert!(at_optimum > 0.0, "{at_optimum}");
    for &step in &[1e-3, 1e-2] {
        for &(da, db) in &[(step, 0.0), (-step, 0.0), (0.0, step), (0.0, -step)] {
            let neighbour = negative_log_likelihood(result.a + da, result.b + db);
            assert!(
                neighbour > at_optimum,
                "moving by ({da}, {db}) lowered the objective: {neighbour} < {at_optimum}"
            );
        }
    }
}

/// A repeated maximum-likelihood fit continues from the previous result,
/// because the source writes the result back into its start parameters.
#[test]
fn mle_fit_is_stateful_across_calls() {
    let samples = read_gumbel_samples();
    let weights = vec![1.0; samples.len()];
    let mut fitter =
        GumbelMaxLikelihoodFitter::with_initial_parameters(MleFitResult::new(4.0, 2.0));
    let first = fitter.fit_weighted(&samples, &weights).unwrap();
    assert_eq!(fitter.initial_parameters(), first);
    let second = fitter.fit_weighted(&samples, &weights).unwrap();
    // Restarting at the optimum keeps it there.
    close(second.a, first.a, 1e-6);
    close(second.b, first.b, 1e-6);
}
