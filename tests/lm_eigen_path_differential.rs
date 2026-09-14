// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Levenberg-Marquardt evaluation-path differential (package B3b-LM-FIDELITY).
//!
//! `tests/lm_budget_differential.rs` pins how many evaluations Eigen spends and
//! where it stops; this file pins *which* points it evaluates. For 141 trace
//! fits - the eight `GaussTraceFitter_test`/`EGHTraceFitter_test` fits, the 50
//! `FeatureFinderCentroided_1` seed fits and four degenerate fits of the C2
//! oracle, and the 79 inputs of the B4-GAUSS review (edge cases, short budgets
//! and seeded random traces) - the fixture holds, per oracle build of Eigen
//! 5.0.1 configured as `TraceFitter::optimize_` configures it, the status,
//! `nfev`, `njev`, the final parameters and a digest of every residual-evaluation
//! argument in order.
//!
//! * **Tier 1 (executed differential).** The oracle is
//!   `../oracle/lm-eigen-path` (`run.sh`, `manifest.json`): stock
//!   `Eigen::LevenbergMarquardt` around the libraries' own `GaussTraceFunctor`
//!   and `EGHTraceFunctor`, byte-identical to the same driver around the
//!   written-out functors below, whose final `x` equals
//!   `GaussTraceFitter::fit`/`EGHTraceFitter::fit` of the library it was built
//!   against in 141 of 141 fits on both platforms:
//!   - `linux-x86_64-release`: the Linux x86_64 Release build the benchmark uses
//!     (gcc 14.4 `-O3 -mssse3 -ffp-contract=off`, SSE2 packets, no FMA);
//!   - `macos-arm64-sdk`: the product SDK (Apple clang `-O0 -ffp-contract=off`,
//!     NEON packets with `EIGEN_VECTORIZE_FMA`);
//!   - `macos-arm64-nofma`: an experiment, the SDK configuration with
//!     `__ARM_FEATURE_FMA` hidden from Eigen.
//! * **What is asserted.** On Linux x86_64 with glibc, every fit must reproduce
//!   `linux-x86_64-release` exactly: status, `nfev`, `njev`, the bits of `x` and
//!   of every evaluation argument. On macOS arm64 every fit must reproduce
//!   `macos-arm64-nofma` exactly. The residuals depend on the platform's `exp`,
//!   so each expectation is checked only where the oracle's libm is the one
//!   linked; the start residual norm is compared first so that a different libm
//!   is reported as such.
//! * **What is not asserted.** Against `macos-arm64-sdk` the solver agrees in
//!   21 of 141 paths: Eigen fuses its packet multiply-adds with FMA on arm64 and
//!   this port does not (a pending decision, recorded in
//!   `docs/DISTRIBUTION_FITTERS_SUPPORT.md` section 1). The ignored
//!   `macos_arm64_sdk_gap_report` measures it.
//!
//! NaN compares equal to NaN in every comparison: neither Eigen nor rustc
//! preserves NaN sign or payload through arithmetic, and one degenerate fit
//! (`review/inf_rt_first`) differs from the C++ only there.

use std::cell::{Cell, RefCell};
use std::collections::BTreeMap;

use openms::math::fitters::levenberg_marquardt::{
    DenseMatrix, LmParameters, minimize, stable_norm,
};

const PROBLEMS: usize = 141;
const CANONICAL_NAN: u64 = 0x7ff8_0000_0000_0000;

// ---------------------------------------------------------------------------
// Fixture
// ---------------------------------------------------------------------------

#[derive(Debug)]
struct Trace {
    theoretical: f64,
    rt: Vec<f64>,
    intensity: Vec<f32>,
}

#[derive(Debug)]
struct Problem {
    name: String,
    egh: bool,
    weighted: bool,
    budget: usize,
    baseline: f64,
    traces: Vec<Trace>,
    x_init: Vec<f64>,
}

impl Problem {
    fn values(&self) -> usize {
        self.traces.iter().map(|t| t.rt.len()).sum()
    }
}

#[derive(Clone, Debug, PartialEq)]
struct Outcome {
    budget: usize,
    status: i32,
    nfev: usize,
    njev: usize,
    x: Vec<u64>,
    evaluations: usize,
    path: u64,
    start_fnorm: u64,
}

fn bits64(token: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(token, 16).expect("f64 bit pattern"))
}

fn bits32(token: &str) -> f32 {
    f32::from_bits(u32::from_str_radix(token, 16).expect("f32 bit pattern"))
}

fn data(name: &str) -> String {
    let path = format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"))
}

/// Problems in the fixture format of `lm_budget_differential/c2_trace_fits.txt`
/// (`problem`, `budget`, `baseline`, `trace`, `x_init`, `end`; other tags are
/// the C2 sweep and are skipped). `budget` defaults to 500, the C2 natural end.
fn parse_problems(text: &str, prefix: &str) -> Vec<Problem> {
    let mut problems = Vec::new();
    let mut current: Option<Problem> = None;
    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let tokens: Vec<&str> = line.split(' ').collect();
        match tokens[0] {
            "problem" => {
                assert!(current.is_none(), "unterminated problem before {line}");
                let name = tokens[1].strip_prefix("review/").unwrap_or(tokens[1]);
                current = Some(Problem {
                    name: format!("{prefix}{name}"),
                    egh: tokens[2] == "egh",
                    weighted: tokens[3] == "1",
                    budget: 500,
                    baseline: 0.0,
                    traces: Vec::new(),
                    x_init: Vec::new(),
                });
            }
            "end" => problems.push(current.take().expect("end without problem")),
            tag => {
                let problem = current.as_mut().expect("line outside a problem");
                match tag {
                    "budget" => problem.budget = tokens[1].parse().expect("budget"),
                    "baseline" => problem.baseline = bits64(tokens[1]),
                    "trace" => {
                        let mut trace = Trace {
                            theoretical: bits64(tokens[1]),
                            rt: Vec::new(),
                            intensity: Vec::new(),
                        };
                        assert_eq!(tokens[2..].len() % 2, 0, "unpaired peak");
                        for pair in tokens[2..].chunks(2) {
                            trace.rt.push(bits64(pair[0]));
                            trace.intensity.push(bits32(pair[1]));
                        }
                        problem.traces.push(trace);
                    }
                    "x_init" => problem.x_init = tokens[1..].iter().map(|t| bits64(t)).collect(),
                    _ => {}
                }
            }
        }
    }
    assert!(current.is_none(), "fixture ends inside a problem");
    problems
}

fn load_problems() -> Vec<Problem> {
    let mut problems = parse_problems(&data("lm_budget_differential/c2_trace_fits.txt"), "c2/");
    problems.extend(parse_problems(
        &data("lm_eigen_path_differential/review_problems.txt"),
        "review/",
    ));
    problems
}

struct Fixture {
    outcomes: BTreeMap<(String, String), Outcome>,
    library: BTreeMap<(String, String), (String, Vec<u64>)>,
}

fn hex_list(token: &str) -> Vec<u64> {
    token
        .split(',')
        .map(|t| u64::from_str_radix(t, 16).expect("bit pattern"))
        .collect()
}

fn load_fixture() -> Fixture {
    let mut outcomes = BTreeMap::new();
    let mut library = BTreeMap::new();
    for line in data("lm_eigen_path_differential/outcomes.txt").lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let t: Vec<&str> = line.split(' ').collect();
        match t[0] {
            "outcome" => {
                let outcome = Outcome {
                    budget: t[3].parse().expect("budget"),
                    status: t[4].parse().expect("status"),
                    nfev: t[5].parse().expect("nfev"),
                    njev: t[6].parse().expect("njev"),
                    x: hex_list(t[7]),
                    evaluations: t[8].parse().expect("evaluations"),
                    path: u64::from_str_radix(t[9], 16).expect("digest"),
                    start_fnorm: u64::from_str_radix(t[10], 16).expect("fnorm"),
                };
                let key = (t[1].to_string(), t[2].to_string());
                assert!(outcomes.insert(key, outcome).is_none(), "duplicate {line}");
            }
            "libfit" => {
                let key = (t[1].to_string(), t[2].to_string());
                let value = (t[3].to_string(), hex_list(t[4]));
                assert!(library.insert(key, value).is_none(), "duplicate {line}");
            }
            other => panic!("unknown fixture tag {other}"),
        }
    }
    Fixture { outcomes, library }
}

fn canonical(bits: u64) -> u64 {
    if f64::from_bits(bits).is_nan() {
        CANONICAL_NAN
    } else {
        bits
    }
}

/// FNV-1a 64 over the little-endian bytes of each value, NaN canonicalised,
/// as `../oracle/lm-eigen-path/make_fixture.py` digests the C++ paths.
fn digest(values: &[f64]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for value in values {
        for byte in canonical(value.to_bits()).to_le_bytes() {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    hash
}

// ---------------------------------------------------------------------------
// Functors: GaussTraceFitter.cpp:145-204 and EGHTraceFitter.cpp:35-145 (core
// bc9cc12), operation for operation, as in tests/lm_budget_differential.rs.
// ---------------------------------------------------------------------------

fn residuals(p: &Problem, x: &[f64], fvec: &mut [f64]) {
    let mut count = 0;
    if !p.egh {
        let (height, x0, sig) = (x[0], x[1], x[2]);
        let c_fac = -0.5 / (sig * sig);
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
        return;
    }
    let (h, tr, sigma, tau) = (x[0], x[1], x[2], x[3]);
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

fn jacobian(p: &Problem, x: &[f64], jac: &mut DenseMatrix) {
    let mut count = 0;
    if !p.egh {
        let (height, x0, sig) = (x[0], x[1], x[2]);
        let sig_sq = sig * sig;
        let inv_siq2 = 1.0 / sig_sq;
        let sig_3 = sig * sig_sq;
        let inv_sig3 = 1.0 / sig_3;
        let c_fac = -0.5 / sig_sq;
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
        return;
    }
    let (h, tr, sigma, tau) = (x[0], x[1], x[2].abs(), x[3]);
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

/// Fit one problem as `TraceFitter::optimize_` does and summarise it in the
/// fixture's terms.
fn run(problem: &Problem) -> Outcome {
    let path = RefCell::new(Vec::new());
    let fev = Cell::new(0usize);
    let jev = Cell::new(0usize);
    let mut x = problem.x_init.clone();
    let parameters = LmParameters {
        max_fev: problem.budget,
        ..LmParameters::default()
    };
    let status = minimize(
        &mut x,
        problem.values(),
        |v: &[f64], f: &mut [f64]| {
            fev.set(fev.get() + 1);
            path.borrow_mut().extend_from_slice(v);
            residuals(problem, v, f);
        },
        |v: &[f64], jac: &mut DenseMatrix| {
            jev.set(jev.get() + 1);
            jacobian(problem, v, jac);
            0
        },
        &parameters,
    );
    let path = path.into_inner();
    let mut start = vec![0.0; problem.values()];
    residuals(problem, &problem.x_init, &mut start);
    Outcome {
        budget: problem.budget,
        status: status.code(),
        nfev: fev.get(),
        njev: jev.get(),
        x: x.iter().map(|v| canonical(v.to_bits())).collect(),
        evaluations: path.len() / problem.x_init.len(),
        path: digest(&path),
        start_fnorm: canonical(stable_norm(&start).to_bits()),
    }
}

fn expected(fixture: &Fixture, oracle: &str, name: &str) -> Outcome {
    let mut outcome = fixture
        .outcomes
        .get(&(oracle.to_string(), name.to_string()))
        .unwrap_or_else(|| panic!("{oracle} has no outcome for {name}"))
        .clone();
    outcome.x = outcome.x.iter().map(|&b| canonical(b)).collect();
    outcome.start_fnorm = canonical(outcome.start_fnorm);
    outcome
}

/// Reproduce one oracle exactly, after checking that this platform's libm gives
/// the oracle's start residual norms.
fn assert_reproduces(oracle: &str) {
    let problems = load_problems();
    let fixture = load_fixture();
    let runs: Vec<Outcome> = problems.iter().map(run).collect();
    let libm: Vec<&str> = problems
        .iter()
        .zip(&runs)
        .filter(|(p, r)| r.start_fnorm != expected(&fixture, oracle, &p.name).start_fnorm)
        .map(|(p, _)| p.name.as_str())
        .collect();
    assert!(
        libm.is_empty(),
        "the start residual norms differ from {oracle} for {} problems ({:?}): this platform's \
         exp is not the oracle's, so its evaluation paths cannot be compared",
        libm.len(),
        &libm[..libm.len().min(5)]
    );
    let mut mismatches = Vec::new();
    for (problem, got) in problems.iter().zip(&runs) {
        let want = expected(&fixture, oracle, &problem.name);
        if *got != want {
            mismatches.push(format!("{}: got {got:?}, {oracle} {want:?}", problem.name));
        }
    }
    assert!(
        mismatches.is_empty(),
        "{} of {PROBLEMS} fits leave the {oracle} evaluation path:\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
}

// ---------------------------------------------------------------------------
// The fixture itself
// ---------------------------------------------------------------------------

/// Every oracle covers every problem at the problem's budget, and every library
/// oracle records a fit for each.
#[test]
fn the_fixture_covers_every_problem_for_every_oracle() {
    let problems = load_problems();
    let fixture = load_fixture();
    assert_eq!(problems.len(), PROBLEMS);
    assert_eq!(
        problems
            .iter()
            .filter(|p| p.name.starts_with("c2/"))
            .count(),
        62
    );
    for oracle in [
        "linux-x86_64-release",
        "macos-arm64-sdk",
        "macos-arm64-nofma",
    ] {
        for problem in &problems {
            let outcome = expected(&fixture, oracle, &problem.name);
            assert_eq!(outcome.budget, problem.budget, "{oracle} {}", problem.name);
            assert_eq!(
                outcome.x.len(),
                problem.x_init.len(),
                "{oracle} {}",
                problem.name
            );
            // Every Jacobian here is analytic, so each evaluation is one nfev.
            assert_eq!(
                outcome.evaluations, outcome.nfev,
                "{oracle} {}",
                problem.name
            );
        }
    }
    assert_eq!(fixture.outcomes.len(), 3 * PROBLEMS);
    assert_eq!(fixture.library.len(), 2 * PROBLEMS);
}

/// The oracle's replica is the library: for both library builds, the replica's
/// final `x` is the parameter set `GaussTraceFitter::fit` (with `|sigma|`) or
/// `EGHTraceFitter::fit` leaves behind, in every fit.
#[test]
fn the_oracle_replicas_equal_the_library_fits() {
    let problems = load_problems();
    let fixture = load_fixture();
    for oracle in ["linux-x86_64-release", "macos-arm64-sdk"] {
        for problem in &problems {
            let replica = expected(&fixture, oracle, &problem.name);
            let (_, params) = &fixture.library[&(oracle.to_string(), problem.name.clone())];
            let mut x: Vec<f64> = replica.x.iter().map(|&b| f64::from_bits(b)).collect();
            if !problem.egh {
                x[2] = x[2].abs();
            }
            for (i, (&got, &want)) in x.iter().zip(params).enumerate() {
                let want = f64::from_bits(want);
                assert!(
                    got.to_bits() == want.to_bits() || (got.is_nan() && want.is_nan()),
                    "{oracle} {} parameter {i}: replica {got:e}, library {want:e}",
                    problem.name
                );
            }
        }
    }
}

/// On one machine the SDK build and its FMA-hidden twin share `exp` and
/// `stableNorm` (which has no multiply-add), so their start norms agree bit for
/// bit; only the solver's fused lanes separate them. Checks the oracle, not Rust.
#[test]
fn fma_changes_the_macos_paths_but_not_their_start() {
    let problems = load_problems();
    let fixture = load_fixture();
    let mut paths_differ = 0;
    for problem in &problems {
        let sdk = expected(&fixture, "macos-arm64-sdk", &problem.name);
        let nofma = expected(&fixture, "macos-arm64-nofma", &problem.name);
        assert_eq!(sdk.start_fnorm, nofma.start_fnorm, "{}", problem.name);
        paths_differ += usize::from(sdk.path != nofma.path);
    }
    assert!(
        paths_differ > 0,
        "the FMA axis must matter for this fixture"
    );
}

// ---------------------------------------------------------------------------
// Tier 1: the solver against the executed Eigen
// ---------------------------------------------------------------------------

/// The Linux x86_64 Release build of OpenMS, the benchmark reference: every one
/// of the 141 fits follows its evaluation path bit for bit.
#[cfg(all(target_os = "linux", target_arch = "x86_64", target_env = "gnu"))]
#[test]
fn minimize_reproduces_the_linux_x86_64_release_paths() {
    assert_reproduces("linux-x86_64-release");
}

/// macOS arm64 with Apple libm: every fit follows Eigen's path when Eigen does
/// not fuse its packet multiply-adds.
#[cfg(all(target_os = "macos", target_arch = "aarch64"))]
#[test]
fn minimize_reproduces_eigen_without_fma_on_macos_arm64() {
    assert_reproduces("macos-arm64-nofma");
}

/// Measurement only: agreement with the product SDK itself, whose Eigen fuses
/// the packet multiply-adds. Run with
/// `cargo test --test lm_eigen_path_differential -- --ignored --nocapture`.
#[cfg(all(target_os = "macos", target_arch = "aarch64"))]
#[test]
#[ignore = "measurement of the pending FMA decision; asserts nothing"]
fn macos_arm64_sdk_gap_report() {
    let problems = load_problems();
    let fixture = load_fixture();
    let (mut paths, mut xs, mut statuses) = (0, 0, 0);
    for problem in &problems {
        let got = run(problem);
        let want = expected(&fixture, "macos-arm64-sdk", &problem.name);
        paths += usize::from(got == want);
        xs += usize::from(got.x == want.x);
        statuses += usize::from(got.status == want.status);
    }
    eprintln!(
        "macos-arm64-sdk: identical paths {paths}/{PROBLEMS}, identical final x {xs}/{PROBLEMS}, \
         identical status {statuses}/{PROBLEMS}"
    );
}
