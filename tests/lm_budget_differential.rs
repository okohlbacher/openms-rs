// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Levenberg-Marquardt evaluation-budget differential (package B3-LM).
//!
//! `TraceFitter::optimize_` (`TraceFitter.cpp:101-135` at core `bc9cc12`) hands
//! `fit:max_iterations` to `Eigen::LevenbergMarquardt` as `maxfev`, so a fit that
//! hits the budget stops at a point that depends on exactly how Eigen counts
//! evaluations and where it tests the count. This file pins that accounting for
//! the backend the crate ships, `openms::math::fitters::levenberg_marquardt::minimize`,
//! against the C2 class-level oracle (`../oracle/featurefinder-picked`, manifest
//! sha256 `7f6adefb...`), budget by budget.
//!
//! * **Tier 1 (executed differential).** `tests/data/lm_budget_differential/c2_trace_fits.txt`
//!   holds, for the eight `GaussTraceFitter_test`/`EGHTraceFitter_test` fits
//!   (theoretical intensities 0.8/0.2 and 0.4/0.6, weighted and unweighted), the
//!   50 Gauss and EGH fits of the 25 `FeatureFinderCentroided_1` seeds, and four
//!   degenerate inputs, the Eigen outcome (status, nfev, njev, x) at every
//!   `max_fev` from 1 to 500. `x` is library output; status, nfev and njev come
//!   from C2's `optimize_` replica, which equals the library fit bit for bit at
//!   every recorded budget. The transcription must reproduce status, nfev and
//!   njev exactly at all 29,004 budgets, and `x` within the declared tolerance.
//! * **Tier 4 (derived from Eigen's source).** For the Gauss, Gamma, Gumbel and
//!   Gumbel maximum-likelihood fitter problems, which have no C++ budget sweep,
//!   the budget outcomes must follow Eigen's counting rule (see
//!   `distribution_fitter_budgets_follow_eigen_accounting`).
//! * **Measurement, ignored by default.** The `levenberg-marquardt =0.14.0`
//!   adapter with exact `maxfev` emulation that B3-LM evaluated, and the gate
//!   report that rejected it. Reproduce with
//!   `cargo test --test lm_budget_differential -- --ignored --nocapture`.
//!
//! The residuals and Jacobians are written out here from the class-test
//! functors (`GaussTraceFitter.cpp:145-204`, `EGHTraceFitter.cpp:35-145`) and
//! from `src/math/fitters/{gauss,gamma,gumbel,gumbel_max_likelihood}.rs`; nothing
//! is imported from the trace-fitter modules.
//!
//! **Tolerance for `x`.** Relative `1e-9` with an absolute floor of `1e-12`,
//! per the per-group rule of `docs/DIFFERENTIAL_VALIDATION.md`: the oracle ran on
//! macOS arm64 (Apple `libm`, NEON `stableNorm`), the gate runs on Linux x86_64.
//! Measured on the gate node: every component is within `6.4e-10` relative
//! except the EGH class-test `tau`, whose true value is zero (the data are a
//! symmetric Gaussian); it is `~4e-15 s` and differs by at most `1.01e-15 s`.
//! Status, nfev and njev are exact. See `docs/DISTRIBUTION_FITTERS_SUPPORT.md`.

// The class-test point lists carry more digits than `f64` holds; they are kept
// character for character, as in `tests/math_distribution_fitters.rs`.
#![allow(clippy::excessive_precision)]

use std::cell::Cell;

use openms::math::fitters::levenberg_marquardt::{
    DenseMatrix, LmParameters, LmStatus, minimize, numerical_jacobian, stable_norm,
};

/// Relative tolerance for fitted parameters compared with the C2 oracle.
const RELATIVE: f64 = 1e-9;
/// Absolute floor for fitted parameters whose value is zero up to rounding.
const ABSOLUTE_FLOOR: f64 = 1e-12;
/// The largest budget C2 swept.
const N_MAX: usize = 500;

// ---------------------------------------------------------------------------
// C2 fixture
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq)]
enum Shape {
    Gauss,
    Egh,
}

#[derive(Debug)]
struct Trace {
    theoretical: f64,
    rt: Vec<f64>,
    intensity: Vec<f32>,
}

#[derive(Clone, Debug)]
struct Outcome {
    budget: usize,
    status: i32,
    nfev: usize,
    njev: usize,
    x: Vec<f64>,
}

#[derive(Debug)]
struct Problem {
    name: String,
    shape: Shape,
    weighted: bool,
    baseline: f64,
    traces: Vec<Trace>,
    x_init: Vec<f64>,
    residuals: Option<Vec<f64>>,
    jacobian: Option<Vec<f64>>,
    fnorm: Option<f64>,
    outcomes: Vec<Outcome>,
    natural: Option<(usize, usize)>,
}

impl Problem {
    fn values(&self) -> usize {
        self.traces.iter().map(|t| t.rt.len()).sum()
    }

    /// The oracle outcome at `budget`, expanding the `natural K 500` line: C2
    /// verified that every budget from `K` to 500 ends as budget `K` does.
    fn oracle(&self, budget: usize) -> Option<&Outcome> {
        if let Some(found) = self.outcomes.iter().find(|o| o.budget == budget) {
            return Some(found);
        }
        let (k, last) = self.natural?;
        if (k..=last).contains(&budget) {
            self.outcomes.iter().find(|o| o.budget == k)
        } else {
            None
        }
    }

    fn budgets(&self) -> impl Iterator<Item = (usize, &Outcome)> {
        (1..=N_MAX).filter_map(move |b| self.oracle(b).map(|o| (b, o)))
    }
}

fn bits64(token: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(token, 16).expect("f64 bit pattern"))
}

fn bits32(token: &str) -> f32 {
    f32::from_bits(u32::from_str_radix(token, 16).expect("f32 bit pattern"))
}

fn load_fixture() -> Vec<Problem> {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/lm_budget_differential/c2_trace_fits.txt"
    );
    let text = std::fs::read_to_string(path).expect("C2 fixture");
    let mut problems = Vec::new();
    let mut current: Option<Problem> = None;
    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let tokens: Vec<&str> = line.split(' ').collect();
        let (tag, rest) = (tokens[0], &tokens[1..]);
        if tag == "problem" {
            assert!(current.is_none(), "unterminated problem before {line}");
            current = Some(Problem {
                name: rest[0].to_string(),
                shape: match rest[1] {
                    "gauss" => Shape::Gauss,
                    "egh" => Shape::Egh,
                    other => panic!("unknown shape {other}"),
                },
                weighted: rest[2] == "1",
                baseline: 0.0,
                traces: Vec::new(),
                x_init: Vec::new(),
                residuals: None,
                jacobian: None,
                fnorm: None,
                outcomes: Vec::new(),
                natural: None,
            });
            continue;
        }
        let problem = current.as_mut().expect("line outside a problem");
        let floats = || rest.iter().map(|t| bits64(t)).collect::<Vec<f64>>();
        match tag {
            "baseline" => problem.baseline = bits64(rest[0]),
            "trace" => {
                let mut trace = Trace {
                    theoretical: bits64(rest[0]),
                    rt: Vec::new(),
                    intensity: Vec::new(),
                };
                assert_eq!(rest[1..].len() % 2, 0, "unpaired peak in {}", problem.name);
                for pair in rest[1..].chunks(2) {
                    trace.rt.push(bits64(pair[0]));
                    trace.intensity.push(bits32(pair[1]));
                }
                problem.traces.push(trace);
            }
            "x_init" => problem.x_init = floats(),
            "residuals" => problem.residuals = Some(floats()),
            "jacobian" => problem.jacobian = Some(floats()),
            "fnorm" => problem.fnorm = Some(bits64(rest[0])),
            "outcome" => problem.outcomes.push(Outcome {
                budget: rest[0].parse().expect("budget"),
                status: rest[1].parse().expect("status"),
                nfev: rest[2].parse().expect("nfev"),
                njev: rest[3].parse().expect("njev"),
                x: rest[4..].iter().map(|t| bits64(t)).collect(),
            }),
            "natural" => {
                problem.natural = Some((
                    rest[0].parse().expect("natural budget"),
                    rest[1].parse().expect("last budget"),
                ))
            }
            "end" => problems.push(current.take().expect("problem")),
            other => panic!("unknown fixture tag {other}"),
        }
    }
    assert!(current.is_none(), "fixture ends inside a problem");
    problems
}

// ---------------------------------------------------------------------------
// Trace functors, written out from GaussTraceFitter.cpp:145-204 and
// EGHTraceFitter.cpp:35-145 (core bc9cc12) in the source's operation order.
// ---------------------------------------------------------------------------

fn trace_residuals(p: &Problem, x: &[f64], fvec: &mut [f64]) {
    match p.shape {
        Shape::Gauss => gauss_trace_residuals(p, x, fvec),
        Shape::Egh => egh_trace_residuals(p, x, fvec),
    }
}

fn trace_jacobian(p: &Problem, x: &[f64], jac: &mut DenseMatrix) {
    match p.shape {
        Shape::Gauss => gauss_trace_jacobian(p, x, jac),
        Shape::Egh => egh_trace_jacobian(p, x, jac),
    }
}

/// `(baseline + theo * height * exp(c_fac * pow2(rt - x0)) - intensity) * weight`
/// with `c_fac = -0.5 / pow2(sig)`.
fn gauss_trace_residuals(p: &Problem, x: &[f64], fvec: &mut [f64]) {
    let (height, x0, sig) = (x[0], x[1], x[2]);
    let c_fac = -0.5 / (sig * sig);
    let mut count = 0;
    for trace in &p.traces {
        let weight = if p.weighted { trace.theoretical } else { 1.0 };
        for (rt, intensity) in trace.rt.iter().zip(&trace.intensity) {
            let d = rt - x0;
            fvec[count] = (p.baseline + trace.theoretical * height * (c_fac * (d * d)).exp()
                - f64::from(*intensity))
                * weight;
            count += 1;
        }
    }
}

/// The source's Jacobian, including the `0.125` factor on the sigma column.
fn gauss_trace_jacobian(p: &Problem, x: &[f64], jac: &mut DenseMatrix) {
    let (height, x0, sig) = (x[0], x[1], x[2]);
    let sig_sq = sig * sig;
    let inv_siq2 = 1.0 / sig_sq;
    let sig_3 = sig * sig_sq;
    let inv_sig3 = 1.0 / sig_3;
    let c_fac = -0.5 / sig_sq;
    let mut count = 0;
    for trace in &p.traces {
        let weight = if p.weighted { trace.theoretical } else { 1.0 };
        for &rt in &trace.rt {
            let e = (c_fac * ((rt - x0) * (rt - x0))).exp();
            jac.set(count, 0, trace.theoretical * e * weight);
            jac.set(
                count,
                1,
                trace.theoretical * height * e * (rt - x0) * inv_siq2 * weight,
            );
            jac.set(
                count,
                2,
                0.125
                    * trace.theoretical
                    * height
                    * e
                    * ((rt - x0) * (rt - x0))
                    * inv_sig3
                    * weight,
            );
            count += 1;
        }
    }
}

/// EGH residual: zero model where `2 sigma^2 + tau (t - tR) <= 0`; signed sigma.
fn egh_trace_residuals(p: &Problem, x: &[f64], fvec: &mut [f64]) {
    let (h, tr, sigma, tau) = (x[0], x[1], x[2], x[3]);
    let mut count = 0;
    for trace in &p.traces {
        let weight = if p.weighted { trace.theoretical } else { 1.0 };
        for (rt, intensity) in trace.rt.iter().zip(&trace.intensity) {
            let t_diff = rt - tr;
            let t_diff2 = t_diff * t_diff;
            let denominator = 2.0 * sigma * sigma + tau * t_diff;
            let fegh = if denominator > 0.0 {
                p.baseline + trace.theoretical * h * (-t_diff2 / denominator).exp()
            } else {
                0.0
            };
            fvec[count] = (fegh - f64::from(*intensity)) * weight;
            count += 1;
        }
    }
}

/// EGH Jacobian: uses `|sigma|` where the residual uses signed sigma.
fn egh_trace_jacobian(p: &Problem, x: &[f64], jac: &mut DenseMatrix) {
    let (h, tr, sigma, tau) = (x[0], x[1], x[2].abs(), x[3]);
    let mut count = 0;
    for trace in &p.traces {
        let weight = if p.weighted { trace.theoretical } else { 1.0 };
        for &rt in &trace.rt {
            let t_diff = rt - tr;
            let t_diff2 = t_diff * t_diff;
            let denominator = 2.0 * sigma * sigma + tau * t_diff;
            let (d_h, d_tr, d_sigma, d_tau) = if denominator > 0.0 {
                let exp1 = (-t_diff2 / denominator).exp();
                (
                    trace.theoretical * exp1,
                    trace.theoretical * h * exp1 * ((4.0 * sigma * sigma + tau * t_diff) * t_diff)
                        / (denominator * denominator),
                    trace.theoretical * h * exp1 * 4.0 * sigma * t_diff2
                        / (denominator * denominator),
                    trace.theoretical * h * exp1 * t_diff * t_diff2 / (denominator * denominator),
                )
            } else {
                (0.0, 0.0, 0.0, 0.0)
            };
            jac.set(count, 0, d_h * weight);
            jac.set(count, 1, d_tr * weight);
            jac.set(count, 2, d_sigma * weight);
            jac.set(count, 3, d_tau * weight);
            count += 1;
        }
    }
}

// ---------------------------------------------------------------------------
// Counting runner
// ---------------------------------------------------------------------------

/// One fit: the Eigen status code, nfev and njev as Eigen counts them, and `x`.
#[derive(Clone, Debug, PartialEq)]
struct Run {
    status: i32,
    nfev: usize,
    njev: usize,
    x: Vec<f64>,
}

/// A backend with `minimize`'s signature, taking type-erased closures.
type Backend = fn(
    &mut [f64],
    usize,
    &mut dyn FnMut(&[f64], &mut [f64]),
    &mut dyn FnMut(&[f64], &mut DenseMatrix) -> usize,
    &LmParameters,
) -> LmStatus;

/// The shipped backend, the Eigen/MINPACK transcription.
fn transcription(
    x: &mut [f64],
    values: usize,
    residuals: &mut dyn FnMut(&[f64], &mut [f64]),
    jacobian: &mut dyn FnMut(&[f64], &mut DenseMatrix) -> usize,
    parameters: &LmParameters,
) -> LmStatus {
    minimize(x, values, residuals, jacobian, parameters)
}

fn budget(max_fev: usize) -> LmParameters {
    LmParameters {
        max_fev,
        ..LmParameters::default()
    }
}

/// Fit a trace problem, counting residual and Jacobian calls. With an analytic
/// Jacobian Eigen's nfev is the number of residual evaluations and njev the
/// number of Jacobian evaluations (`LevenbergMarquardt.h` lines 182, 207, 262).
fn run_trace(problem: &Problem, max_fev: usize, backend: Backend) -> Run {
    let fev = Cell::new(0usize);
    let jev = Cell::new(0usize);
    let mut x = problem.x_init.clone();
    let mut residuals = |v: &[f64], f: &mut [f64]| {
        fev.set(fev.get() + 1);
        trace_residuals(problem, v, f);
    };
    let mut jacobian = |v: &[f64], jac: &mut DenseMatrix| {
        jev.set(jev.get() + 1);
        trace_jacobian(problem, v, jac);
        0
    };
    let status = backend(
        &mut x,
        problem.values(),
        &mut residuals,
        &mut jacobian,
        &budget(max_fev),
    );
    Run {
        status: status.code(),
        nfev: fev.get(),
        njev: jev.get(),
        x,
    }
}

/// Whether `got` matches `want` under the declared tolerance.
fn close(got: f64, want: f64) -> bool {
    if got == want || got.to_bits() == want.to_bits() || (got.is_nan() && want.is_nan()) {
        return true;
    }
    let difference = (got - want).abs();
    difference <= RELATIVE * got.abs().max(want.abs()) || difference <= ABSOLUTE_FLOOR
}

fn assert_budget_sweeps(filter: impl Fn(&Problem) -> bool, expected_problems: usize) {
    let problems = load_fixture();
    let mut checked_problems = 0usize;
    let mut checked_budgets = 0usize;
    for problem in problems.iter().filter(|p| filter(p)) {
        checked_problems += 1;
        for (b, oracle) in problem.budgets() {
            checked_budgets += 1;
            let run = run_trace(problem, b, transcription);
            assert_eq!(
                (run.status, run.nfev, run.njev),
                (oracle.status, oracle.nfev, oracle.njev),
                "{} at max_fev {b}: (status, nfev, njev) differs from C2",
                problem.name
            );
            for (i, (got, want)) in run.x.iter().zip(&oracle.x).enumerate() {
                assert!(
                    close(*got, *want),
                    "{} at max_fev {b}: x[{i}] = {got:e}, C2 {want:e}",
                    problem.name
                );
            }
        }
    }
    assert_eq!(checked_problems, expected_problems);
    assert!(checked_budgets >= expected_problems);
}

// ---------------------------------------------------------------------------
// Tier 1: the transcription against C2
// ---------------------------------------------------------------------------

/// The fixture holds what the generator promised: 8 class-test and 50
/// `FeatureFinderCentroided_1` sweeps over budgets 1..500, and 4 degenerate
/// fits at 500 - 29,004 budgets in all.
#[test]
fn the_c2_fixture_covers_every_budget_of_every_problem() {
    let problems = load_fixture();
    let count = |prefix: &str| {
        problems
            .iter()
            .filter(|p| p.name.starts_with(prefix))
            .count()
    };
    assert_eq!(count("classtest/"), 8);
    assert_eq!(count("ffc1/"), 50);
    assert_eq!(count("degenerate/"), 4);
    let mut budgets = 0usize;
    for problem in &problems {
        assert_eq!(
            problem.x_init.len(),
            if problem.shape == Shape::Gauss { 3 } else { 4 }
        );
        let swept = problem.natural.is_some();
        assert_eq!(
            swept,
            !problem.name.starts_with("degenerate/"),
            "{}",
            problem.name
        );
        let covered = problem.budgets().count();
        assert_eq!(covered, if swept { N_MAX } else { 1 }, "{}", problem.name);
        budgets += covered;
        for outcome in &problem.outcomes {
            assert_eq!(outcome.x.len(), problem.x_init.len());
        }
    }
    assert_eq!(budgets, 29_004);
}

/// The written-out functors evaluate to C2's recorded residuals and Jacobians at
/// the class-test start vectors, and to C2's residual norm at the start of each
/// `FeatureFinderCentroided_1` fit. Measured on the gate node: the class-test
/// residuals and Jacobians are bit-identical; the norms differ by at most
/// `5.4e-16` relative, which is Eigen's vectorized `stableNorm` accumulation.
#[test]
fn the_written_out_functors_reproduce_the_c2_start_evaluations() {
    let problems = load_fixture();
    let mut compared = 0usize;
    for problem in &problems {
        let (m, n) = (problem.values(), problem.x_init.len());
        let mut fvec = vec![0.0; m];
        trace_residuals(problem, &problem.x_init, &mut fvec);
        if let Some(expected) = &problem.residuals {
            assert_eq!(expected.len(), m);
            for (i, (got, want)) in fvec.iter().zip(expected).enumerate() {
                assert!(
                    close(*got, *want),
                    "{} residual {i}: {got:e} vs {want:e}",
                    problem.name
                );
            }
            compared += 1;
        }
        if let Some(expected) = &problem.jacobian {
            assert_eq!(expected.len(), m * n);
            let mut jac = DenseMatrix::zeros(m, n);
            trace_jacobian(problem, &problem.x_init, &mut jac);
            for j in 0..n {
                for i in 0..m {
                    let (got, want) = (jac.at(i, j), expected[j * m + i]);
                    assert!(
                        close(got, want),
                        "{} J({i},{j}): {got:e} vs {want:e}",
                        problem.name
                    );
                }
            }
        }
        if let Some(expected) = problem.fnorm {
            assert!(
                close(stable_norm(&fvec), expected),
                "{} start norm",
                problem.name
            );
            compared += 1;
        }
    }
    assert_eq!(compared, 58);
}

/// `GaussTraceFitter_test` and `EGHTraceFitter_test` fits, `max_fev` 1..500.
#[test]
fn transcription_reproduces_the_class_test_budget_sweeps() {
    assert_budget_sweeps(|p| p.name.starts_with("classtest/"), 8);
}

/// The 25 Gauss fits of `FeatureFinderCentroided_1`, `max_fev` 1..500. Seed 24
/// needs 131 evaluations, the largest natural budget in the fixture.
#[test]
fn transcription_reproduces_the_ffc1_gauss_budget_sweeps() {
    assert_budget_sweeps(
        |p| p.name.starts_with("ffc1/") && p.shape == Shape::Gauss,
        25,
    );
}

/// The 25 EGH fits of `FeatureFinderCentroided_1`, `max_fev` 1..500.
#[test]
fn transcription_reproduces_the_ffc1_egh_budget_sweeps() {
    assert_budget_sweeps(|p| p.name.starts_with("ffc1/") && p.shape == Shape::Egh, 25);
}

/// Degenerate inputs at `max_fev` 500: flat and three-point traces. The EGH
/// starts carry a NaN sigma and an infinite tau, for which Eigen stops at once
/// with `CosinusTooSmall`; the flat Gauss fit drives sigma to about `1.1e5` and
/// stops on the rank-deficient Jacobian with `RelativeErrorTooSmall`.
#[test]
fn transcription_reproduces_the_degenerate_fits() {
    assert_budget_sweeps(|p| p.name.starts_with("degenerate/"), 4);
}

// ---------------------------------------------------------------------------
// Tier 4: Eigen's accounting on the distribution fitter problems
// ---------------------------------------------------------------------------

const GAUSS_MZ: [f64; 7] = [
    240.1000470172,
    240.1002675493,
    240.1004880817,
    240.1007086145,
    240.1009291475,
    240.1011496808,
    240.1013702145,
];

const GAUSS_INTENSITIES: [f64; 7] = [
    61134.39453125,
    111288.5390625,
    163761.46875,
    165861.4375,
    162133.46875,
    120060.5234375,
    71102.1328125,
];

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

/// The first `GaussFitter` class-test case, `GaussFitter_test.cpp:70-97`.
const GAUSS_POINTS: [(f64, f64); 6] = [
    (0.0, 0.01),
    (0.05, 0.2),
    (0.16, 0.63),
    (0.28, 0.99),
    (0.66, 0.03),
    (0.50, 0.36),
];

/// Copy of the private series in `src/math/fitters/gamma.rs`.
fn digamma(x: f64) -> f64 {
    let mut value = x;
    let mut result = 0.0f64;
    while value < 10.0 {
        result -= 1.0 / value;
        value += 1.0;
    }
    let f = 1.0 / (value * value);
    let series = f
        * (1.0 / 12.0
            - f * (1.0 / 120.0
                - f * (1.0 / 252.0
                    - f * (1.0 / 240.0
                        - f * (1.0 / 132.0 - f * (691.0 / 32760.0 - f * (1.0 / 12.0)))))));
    result + value.ln() - 0.5 / value - series
}

fn gumbel_samples() -> Vec<f64> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/gumbel_1d.csv");
    let text = std::fs::read_to_string(path).expect("gumbel_1d.csv");
    text.split([',', '\n', '\r'])
        .map(str::trim)
        .filter(|token| !token.is_empty())
        .map(|token| token.parse::<f64>().expect("sample"))
        .collect()
}

enum Model {
    Gauss(Vec<(f64, f64)>),
    Gamma(Vec<(f64, f64)>),
    Gumbel(Vec<(f64, f64)>),
    MaximumLikelihood(Vec<f64>, Vec<f64>),
}

struct DistributionProblem {
    name: &'static str,
    model: Model,
    x_init: Vec<f64>,
}

/// Every fit the four distribution-fitter class tests and
/// `tests/math_distribution_fitters.rs` run, with their start vectors.
fn distribution_problems() -> Vec<DistributionProblem> {
    let peak: Vec<(f64, f64)> = GAUSS_MZ
        .iter()
        .copied()
        .zip(GAUSS_INTENSITIES.iter().copied())
        .collect();
    let samples = gumbel_samples();
    let unit = vec![1.0; samples.len()];
    let (mut grid_x, mut grid_w) = (Vec::new(), Vec::new());
    let mut xi = -2.0f64;
    while xi <= 8.0 {
        let z = (xi - 2.0) / 0.8;
        grid_w.push((1.0 / 0.8) * (-(z + (-z).exp())).exp());
        grid_x.push(xi);
        xi += 0.1;
    }
    let gauss = |name, x_init: [f64; 3]| DistributionProblem {
        name,
        model: Model::Gauss(GAUSS_POINTS.to_vec()),
        x_init: x_init.to_vec(),
    };
    let gumbel = |name, points: &[(f64, f64)], x_init: [f64; 2]| DistributionProblem {
        name,
        model: Model::Gumbel(points.to_vec()),
        x_init: x_init.to_vec(),
    };
    let gamma = |name, x_init: [f64; 2]| DistributionProblem {
        name,
        model: Model::Gamma(GAMMA_POINTS.to_vec()),
        x_init: x_init.to_vec(),
    };
    vec![
        gauss("gauss/published", [0.06, 3.0, 0.5]),
        DistributionProblem {
            name: "gauss/peak",
            model: Model::Gauss(peak),
            x_init: vec![168324.0, 240.10051, 0.000375375],
        },
        gauss("gauss/far", [0.06, 10.0, 0.5]),
        gauss("gauss/left_edge", [0.06, 0.0, 0.5]),
        gauss("gauss/basin_1", [1.0, 0.3, 0.2]),
        gauss("gauss/basin_2", [1.0, 1.0, 1.0]),
        gauss("gauss/basin_3", [0.5, -1.0, 2.0]),
        gauss("gauss/negative", [-1.0, -1.0, -1.0]),
        gamma("gamma/published", [1.0, 3.0]),
        gamma("gamma/default", [1.0, 5.0]),
        gumbel("gumbel/first", &GUMBEL_POINTS, [1.0, 3.0]),
        gumbel("gumbel/second", &GUMBEL_POINTS_2, [3.0, 3.0]),
        gumbel("gumbel/default_result", &GUMBEL_POINTS, [1.0, 2.0]),
        gumbel("gumbel/default", &GUMBEL_POINTS, [0.25, 0.1]),
        DistributionProblem {
            name: "gumbel_mle/published",
            model: Model::MaximumLikelihood(samples, unit),
            x_init: vec![4.0, 2.0],
        },
        DistributionProblem {
            name: "gumbel_mle/grid",
            model: Model::MaximumLikelihood(grid_x, grid_w),
            x_init: vec![1.0, 1.0],
        },
    ]
}

/// Fit a distribution problem with the closures of its fitter. nfev counts the
/// residual calls plus the evaluations a numerical Jacobian reports consuming,
/// as Eigen's `nfev += df_ret` does.
fn run_distribution(problem: &DistributionProblem, max_fev: usize, backend: Backend) -> Run {
    let fev = Cell::new(0usize);
    let consumed = Cell::new(0usize);
    let jev = Cell::new(0usize);
    let mut x = problem.x_init.clone();
    let parameters = budget(max_fev);
    let status = match &problem.model {
        Model::Gauss(points) => {
            let mut residuals = |p: &[f64], fvec: &mut [f64]| {
                fev.set(fev.get() + 1);
                let (amplitude, center, sigma) = (p[0], p[1], p[2]);
                let sig2 = 2.0 * sigma * sigma;
                for (slot, &(px, py)) in fvec.iter_mut().zip(points) {
                    *slot = amplitude * (-(px - center) * (px - center) / sig2).exp() - py;
                }
            };
            let mut jacobian = |p: &[f64], jac: &mut DenseMatrix| {
                jev.set(jev.get() + 1);
                let (amplitude, center, sigma) = (p[0], p[1], p[2]);
                let sig2 = 2.0 * sigma * sigma;
                let sig3 = 2.0 * sig2 * sigma;
                for (row, &(px, _)) in points.iter().enumerate() {
                    let xd = px - center;
                    let xd2 = xd * xd;
                    let j0 = (-xd2 / sig2).exp();
                    jac.set(row, 0, j0);
                    jac.set(
                        row,
                        1,
                        amplitude * j0 * (-(-2.0 * px + 2.0 * center) / sig2),
                    );
                    jac.set(row, 2, amplitude * j0 * (xd2 / sig3));
                }
                0
            };
            backend(
                &mut x,
                points.len(),
                &mut residuals,
                &mut jacobian,
                &parameters,
            )
        }
        Model::Gamma(points) => {
            let mut residuals = |v: &[f64], fvec: &mut [f64]| {
                fev.set(fev.get() + 1);
                let (b, p) = (v[0], v[1]);
                if b > 0.0 && p > 0.0 {
                    for (slot, &(px, py)) in fvec.iter_mut().zip(points) {
                        *slot =
                            b.powf(p) / libm::tgamma(p) * px.powf(p - 1.0) * (-b * px).exp() - py;
                    }
                } else {
                    for (slot, &(_, py)) in fvec.iter_mut().zip(points) {
                        *slot = -py;
                    }
                }
            };
            let mut jacobian = |v: &[f64], jac: &mut DenseMatrix| {
                jev.set(jev.get() + 1);
                let (b, p) = (v[0], v[1]);
                if b > 0.0 && p > 0.0 {
                    for (row, &(px, _)) in points.iter().enumerate() {
                        let part_dev_b = px.powf(p - 1.0) * (-px * b).exp() / libm::tgamma(p)
                            * (p * b.powf(p - 1.0) - px * b.powf(p));
                        jac.set(row, 0, part_dev_b);
                        let factor = (-b * px).exp() * px.powf(p - 1.0) * b.powf(p)
                            / (libm::tgamma(p) * libm::tgamma(p));
                        let argument =
                            (b.ln() + px.ln()) * libm::tgamma(p) - libm::tgamma(p) * digamma(p);
                        jac.set(row, 1, factor * argument);
                    }
                } else {
                    for row in 0..points.len() {
                        jac.set(row, 0, 0.0);
                        jac.set(row, 1, 0.0);
                    }
                }
                0
            };
            backend(
                &mut x,
                points.len(),
                &mut residuals,
                &mut jacobian,
                &parameters,
            )
        }
        Model::Gumbel(points) => {
            let mut residuals = |v: &[f64], fvec: &mut [f64]| {
                fev.set(fev.get() + 1);
                let (a, b) = (v[0], v[1]);
                for (slot, &(px, py)) in fvec.iter_mut().zip(points) {
                    let z = ((a - px) / b).exp();
                    *slot = (z * (-z).exp()) / b - py;
                }
            };
            let mut jacobian = |v: &[f64], jac: &mut DenseMatrix| {
                jev.set(jev.get() + 1);
                let (a, b) = (v[0], v[1]);
                for (row, &(px, _)) in points.iter().enumerate() {
                    let z = ((a - px) / b).exp();
                    let f = z * (-z).exp();
                    jac.set(row, 0, (f - z * z * (-z).exp()) / (b * b));
                    let dev_z = (px - a) / (b * b);
                    let cum = f * dev_z;
                    jac.set(row, 1, ((cum - z * cum) * b - f) / (b * b));
                }
                0
            };
            backend(
                &mut x,
                points.len(),
                &mut residuals,
                &mut jacobian,
                &parameters,
            )
        }
        Model::MaximumLikelihood(samples, weights) => {
            let objective = |v: &[f64], fvec: &mut [f64]| {
                let sigma = v[1].abs();
                let logsigma = sigma.ln();
                let mut sum = 0.0f64;
                for (&sample, &weight) in samples.iter().zip(weights) {
                    let diff = (sample - v[0]) / sigma;
                    sum += weight * (-logsigma - diff - (-diff).exp());
                }
                fvec[0] = -sum;
                fvec[1] = 0.0;
            };
            let mut residuals = |v: &[f64], fvec: &mut [f64]| {
                fev.set(fev.get() + 1);
                objective(v, fvec);
            };
            let mut jacobian = |v: &[f64], jac: &mut DenseMatrix| {
                jev.set(jev.get() + 1);
                let spent = numerical_jacobian(v, 2, jac, objective);
                consumed.set(consumed.get() + spent);
                spent
            };
            backend(&mut x, 2, &mut residuals, &mut jacobian, &parameters)
        }
    };
    Run {
        status: status.code(),
        nfev: fev.get() + consumed.get(),
        njev: jev.get(),
        x,
    }
}

/// Eigen's budget rule on the four distribution fitters, which have no C++ sweep.
///
/// From `LevenbergMarquardt.h` (the lines C2 records): nfev is 1 after the
/// start, grows by the evaluations a numerical Jacobian reports and by 1 per
/// trial step, and `nfev >= maxfev` is tested once after each trial, after the
/// `ftol`/`xtol` tests. Hence, for every budget `b` from 1 to 500:
///
/// * if the unbounded fit terminates on its own after `N` evaluations, every
///   `b >= N` reproduces it exactly (status, nfev, njev and every bit of `x`);
/// * a `TooManyFunctionEvaluation` stop happens at the first post-trial count
///   `>= b`: exactly `max(b, 2)` with an analytic Jacobian, and within
///   `[b, max(b + n + 1, n + 3)]` with the `n + 1`-evaluation numerical one,
///   whose consecutive post-trial counts differ by 1 or `n + 2`;
/// * any other stop spends no more than that upper bound.
///
/// Fifteen of the 16 problems terminate on their own within 500 evaluations;
/// `gauss/left_edge` spends every budget, which `GaussFitter` turns into
/// `UnableToFit`. The test asserts that both kinds of outcome occur, so the
/// sweep cannot pass vacuously.
#[test]
fn distribution_fitter_budgets_follow_eigen_accounting() {
    let mut natural_problems = 0usize;
    let mut budget_stops = 0usize;
    for problem in distribution_problems() {
        let n = problem.x_init.len();
        let numerical = matches!(problem.model, Model::MaximumLikelihood(..));
        let free = run_distribution(&problem, N_MAX, transcription);
        let natural = free.status != LmStatus::TooManyFunctionEvaluation.code();
        if natural {
            natural_problems += 1;
        }
        for b in 1..=N_MAX {
            let run = run_distribution(&problem, b, transcription);
            if natural && b >= free.nfev {
                assert_eq!(run.status, free.status, "{} at {b}", problem.name);
                assert_eq!(
                    (run.nfev, run.njev),
                    (free.nfev, free.njev),
                    "{} at {b}",
                    problem.name
                );
                let same = run
                    .x
                    .iter()
                    .zip(&free.x)
                    .all(|(a, f)| a.to_bits() == f.to_bits());
                assert!(same, "{} at {b}: {:?} vs {:?}", problem.name, run.x, free.x);
            }
            let upper = if numerical {
                (b + n + 1).max(n + 3)
            } else {
                b.max(2)
            };
            if run.status == LmStatus::TooManyFunctionEvaluation.code() {
                budget_stops += 1;
                assert!(
                    run.nfev >= b && run.nfev <= upper,
                    "{} at {b}: nfev {}",
                    problem.name,
                    run.nfev
                );
                if !numerical {
                    assert_eq!(run.nfev, b.max(2), "{} at {b}", problem.name);
                }
            } else {
                assert!(
                    run.nfev <= upper,
                    "{} at {b}: nfev {}",
                    problem.name,
                    run.nfev
                );
            }
        }
    }
    assert!(natural_problems > 0 && budget_stops > 0);
}

// ---------------------------------------------------------------------------
// Measurement: the levenberg-marquardt crate candidate B3-LM rejected
// ---------------------------------------------------------------------------

/// The `levenberg-marquardt =0.14.0` adapter behind `minimize`'s signature that
/// B3-LM measured, kept so the gate can be re-run when the decision is reopened.
///
/// Configuration: `ftol = xtol = sqrt(f64::EPSILON)` through `with_ftol` and
/// `with_xtol`, `gtol` 0, `stepbound` = Eigen's `factor`, diagonal scaling on,
/// `patience = ceil(max_fev / (n + 1))`. Budget emulation: the problem counts
/// Eigen's nfev (1 at the start, `+ consumed` per Jacobian, `+1` per trial); the
/// first call after a trial applies Eigen's post-trial test `nfev >= max_fev`,
/// returning `None` and restoring the last accepted `x` (snapshotted at every
/// Jacobian call, which the crate makes only at accepted points).
mod candidate {
    use std::cell::{Cell, RefCell};

    use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt, TerminationReason};
    use nalgebra::{DMatrix, DVector, Dyn, storage::Owned};
    use openms::math::fitters::levenberg_marquardt::{DenseMatrix, LmParameters, LmStatus};

    #[derive(Clone, Copy, Debug)]
    struct Budget {
        max_fev: usize,
        nfev: usize,
        started: bool,
        pending: bool,
        stopped: bool,
    }

    impl Budget {
        /// Apply the test Eigen makes after the most recent trial, once.
        fn exhausted(&mut self) -> bool {
            if self.pending {
                self.pending = false;
                if self.nfev >= self.max_fev {
                    self.stopped = true;
                    return true;
                }
            }
            false
        }
    }

    struct Adapter<R, J> {
        params: DVector<f64>,
        values: usize,
        residuals: RefCell<R>,
        jacobian: RefCell<J>,
        fvec: RefCell<Vec<f64>>,
        fjac: RefCell<DenseMatrix>,
        budget: Cell<Budget>,
        snapshot: RefCell<DVector<f64>>,
    }

    impl<R, J> LeastSquaresProblem<f64, Dyn, Dyn> for Adapter<R, J>
    where
        R: FnMut(&[f64], &mut [f64]),
        J: FnMut(&[f64], &mut DenseMatrix) -> usize,
    {
        type ResidualStorage = Owned<f64, Dyn>;
        type JacobianStorage = Owned<f64, Dyn, Dyn>;
        type ParameterStorage = Owned<f64, Dyn>;

        fn set_params(&mut self, x: &DVector<f64>) {
            self.params.copy_from(x);
        }

        fn params(&self) -> DVector<f64> {
            self.params.clone()
        }

        fn residuals(&self) -> Option<DVector<f64>> {
            let mut budget = self.budget.get();
            if budget.started && budget.exhausted() {
                self.budget.set(budget);
                return None;
            }
            let mut fvec = self.fvec.try_borrow_mut().ok()?;
            let mut residuals = self.residuals.try_borrow_mut().ok()?;
            (residuals)(self.params.as_slice(), &mut fvec);
            if budget.started {
                budget.nfev = budget.nfev.saturating_add(1);
                budget.pending = true;
            } else {
                budget.started = true;
                budget.nfev = 1;
            }
            self.budget.set(budget);
            Some(DVector::from_column_slice(&fvec))
        }

        fn jacobian(&self) -> Option<DMatrix<f64>> {
            let mut budget = self.budget.get();
            if let Ok(mut snapshot) = self.snapshot.try_borrow_mut() {
                snapshot.copy_from(&self.params);
            }
            if budget.exhausted() {
                self.budget.set(budget);
                return None;
            }
            let mut fjac = self.fjac.try_borrow_mut().ok()?;
            let mut jacobian = self.jacobian.try_borrow_mut().ok()?;
            let consumed = (jacobian)(self.params.as_slice(), &mut fjac);
            if consumed > 0 {
                budget.nfev = budget.nfev.saturating_add(consumed);
            }
            self.budget.set(budget);
            Some(DMatrix::from_fn(self.values, self.params.len(), |i, j| {
                fjac.at(i, j)
            }))
        }
    }

    /// `minimize` on the crate, with Eigen's status codes.
    pub fn minimize(
        x: &mut [f64],
        values: usize,
        residuals: &mut dyn FnMut(&[f64], &mut [f64]),
        jacobian: &mut dyn FnMut(&[f64], &mut DenseMatrix) -> usize,
        parameters: &LmParameters,
    ) -> LmStatus {
        let n = x.len();
        if n == 0
            || values < n
            || parameters.ftol < 0.0
            || parameters.xtol < 0.0
            || parameters.gtol < 0.0
            || parameters.max_fev == 0
            || parameters.factor <= 0.0
        {
            return LmStatus::ImproperInputParameters;
        }
        let width = n + 1;
        let patience = parameters.max_fev.div_ceil(width).min(usize::MAX / width);
        // `abs` clears the sign of -0.0 and of a negative NaN, which the
        // crate's sign-based assertions would otherwise reject.
        let solver = LevenbergMarquardt::new()
            .with_ftol(parameters.ftol.abs())
            .with_xtol(parameters.xtol.abs())
            .with_gtol(parameters.gtol.abs())
            .with_stepbound(parameters.factor.abs())
            .with_patience(patience);
        let problem = Adapter {
            params: DVector::from_column_slice(x),
            values,
            residuals: RefCell::new(residuals),
            jacobian: RefCell::new(jacobian),
            fvec: RefCell::new(vec![0.0; values]),
            fjac: RefCell::new(DenseMatrix::zeros(values, n)),
            budget: Cell::new(Budget {
                max_fev: parameters.max_fev,
                nfev: 0,
                started: false,
                pending: false,
                stopped: false,
            }),
            snapshot: RefCell::new(DVector::from_column_slice(x)),
        };
        let (problem, report) = solver.minimize(problem);
        let budget = problem.budget.get();
        let over_budget = budget.pending && budget.nfev >= budget.max_fev;
        let status = if budget.stopped {
            LmStatus::TooManyFunctionEvaluation
        } else {
            match report.termination {
                TerminationReason::Converged {
                    ftol: true,
                    xtol: true,
                } => LmStatus::RelativeErrorAndReductionTooSmall,
                TerminationReason::Converged { ftol: true, .. } => {
                    LmStatus::RelativeReductionTooSmall
                }
                TerminationReason::Converged { .. } => LmStatus::RelativeErrorTooSmall,
                TerminationReason::LostPatience => LmStatus::TooManyFunctionEvaluation,
                TerminationReason::NoImprovementPossible(_) if over_budget => {
                    LmStatus::TooManyFunctionEvaluation
                }
                TerminationReason::NoImprovementPossible("ftol") => LmStatus::FtolTooSmall,
                TerminationReason::NoImprovementPossible("xtol") => LmStatus::XtolTooSmall,
                TerminationReason::NoImprovementPossible(_) => LmStatus::GtolTooSmall,
                TerminationReason::Orthogonal | TerminationReason::ResidualsZero => {
                    LmStatus::CosinusTooSmall
                }
                TerminationReason::Numerical(_) => LmStatus::TooManyFunctionEvaluation,
                _ => LmStatus::ImproperInputParameters,
            }
        };
        if budget.stopped {
            if let Ok(snapshot) = problem.snapshot.try_borrow() {
                x.copy_from_slice(snapshot.as_slice());
            }
        } else {
            x.copy_from_slice(problem.params.as_slice());
        }
        status
    }
}

fn max_relative(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| {
            if x == y || x.to_bits() == y.to_bits() || (x.is_nan() && y.is_nan()) {
                0.0
            } else {
                (x - y).abs() / x.abs().max(y.abs())
            }
        })
        .fold(0.0, f64::max)
}

/// The B3-LM gate for the crate candidate (acceptance 5 of the package).
///
/// Status and nfev must equal the reference at every budget; `x` must agree
/// with the transcription within `1e-12` relative, or be no further from C2
/// than the transcription. Prints the report; it asserts nothing, because the
/// candidate failed and is not the shipped backend. Measured on the gate node
/// on 2026-09-14, the report reads: status, nfev and njev equal C2 at 29,003 of
/// 29,004 trace budgets (not `degenerate/flat3_gauss`) and the transcription at
/// all 8,000 distribution budgets; the `x` clause fails at 3,811 trace budgets,
/// 9 of them at budget 500, and 2,533 distribution budgets lie beyond `1e-12`
/// of the transcription, where no C++ sweep can decide the other clause.
///
/// The report ends with the wall time of the 62 trace fits at budget 500 on
/// each backend; that line means something only in a `--release` run.
#[test]
#[ignore = "B3-LM gate report for the rejected levenberg-marquardt candidate; run with --ignored --nocapture"]
fn levenberg_marquardt_crate_candidate_gate_report() {
    let problems = load_fixture();
    let (mut budgets, mut accounting, mut within, mut no_further) =
        (0usize, 0usize, 0usize, 0usize);
    for problem in &problems {
        for (b, oracle) in problem.budgets() {
            budgets += 1;
            let t = run_trace(problem, b, transcription);
            let a = run_trace(problem, b, candidate::minimize);
            if (a.status, a.nfev, a.njev) == (oracle.status, oracle.nfev, oracle.njev) {
                accounting += 1;
            } else {
                println!(
                    "accounting differs: {} at {b}: {a:?} vs C2 {oracle:?}",
                    problem.name
                );
            }
            let versus_transcription = max_relative(&a.x, &t.x);
            if versus_transcription <= 1e-12 {
                within += 1;
            } else if max_relative(&a.x, &oracle.x) <= max_relative(&t.x, &oracle.x) {
                no_further += 1;
            } else if b == N_MAX {
                println!(
                    "x clause fails at the natural end of {}: candidate {:?}, transcription {:?}, C2 {:?}",
                    problem.name, a.x, t.x, oracle.x
                );
            }
        }
    }
    println!(
        "trace budgets {budgets}: accounting equal to C2 {accounting}; x within 1e-12 of the transcription {within}; else no further from C2 {no_further}; x clause failing {}",
        budgets - within - no_further
    );
    let (mut dist_budgets, mut dist_accounting, mut dist_within) = (0usize, 0usize, 0usize);
    for problem in distribution_problems() {
        for b in 1..=N_MAX {
            dist_budgets += 1;
            let t = run_distribution(&problem, b, transcription);
            let a = run_distribution(&problem, b, candidate::minimize);
            if (a.status, a.nfev, a.njev) == (t.status, t.nfev, t.njev) {
                dist_accounting += 1;
            }
            if max_relative(&a.x, &t.x) <= 1e-12 {
                dist_within += 1;
            }
        }
    }
    println!(
        "distribution budgets {dist_budgets}: accounting equal to the transcription {dist_accounting}; x within 1e-12 {dist_within} (no C++ sweep exists for the other clause)"
    );
    const REPEATS: usize = 200;
    for (label, backend) in [
        ("transcription", transcription as Backend),
        ("candidate", candidate::minimize as Backend),
    ] {
        let start = std::time::Instant::now();
        let mut evaluations = 0usize;
        for _ in 0..REPEATS {
            for problem in &problems {
                evaluations += run_trace(problem, N_MAX, backend).nfev;
            }
        }
        println!(
            "timing {label}: {REPEATS} x 62 fits at budget 500, {evaluations} residual evaluations, {:?}",
            start.elapsed()
        );
    }
}
