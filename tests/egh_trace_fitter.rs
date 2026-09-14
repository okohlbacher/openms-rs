// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FEATUREFINDER/EGHTraceFitter` against its class test and two executed
//! product-SDK oracles.
//!
//! Every `START_SECTION` of `EGHTraceFitter_test.cpp` (core bc9cc12) is
//! transcribed with its literals and comparison macros. `TEST_EQUAL` is exact
//! equality; `TEST_REAL_SIMILAR` is ClassTest's `isRealSimilar` with its default
//! tolerances, absolute 1e-5 and ratio 1 + 1e-5. The class test fits once at
//! file scope and mutates the traces inside the `fit` section; each test here
//! rebuilds the state that section left.
//!
//! The oracle replays read two fixtures under `tests/data/egh_trace_fitter/`:
//!
//! - `egh_trace_fitter_c2.tsv`, the EGH records of the C2 class-level oracle
//!   (`../oracle/featurefinder-picked`), extracted by
//!   `../oracle/egh-trace-fitter/extract_c2.py`: the class-test traces at
//!   theoretical intensities 0.8/0.2 and 0.4/0.6, weighted and unweighted, and
//!   three degenerate inputs;
//! - `egh_trace_fitter_oracle.tsv`, the output of
//!   `../oracle/egh-trace-fitter/driver.cpp`: asymmetric traces, edge profiles,
//!   non-positive budgets and hand-picked parameter vectors.
//!
//! Both record the inputs bit for bit, so nothing is regenerated here. Results
//! are compared in three classes, each printing how many comparisons were
//! bit-identical (visible with `--nocapture`):
//!
//! - Exact: booleans, integers, error messages, gnuplot formulas of the
//!   oracle's own parameters, and the budget-boundary pattern (the Rust fit at
//!   `max_iteration = n` equals the Rust fit at 500 bit for bit exactly where the
//!   C++ fits do).
//! - A single evaluation on recorded inputs (functor rows, start points, the
//!   queries of a model set to the oracle's parameters) within
//!   [`DIRECT_RELATIVE`] = 1e-14, residuals on the scale of their terms.
//! - A Levenberg-Marquardt result against the C++ fit within [`FIT_RELATIVE`] =
//!   1e-9, `tau` on the scale of `sigma` only where `tau` is rounding noise.
//!
//! NaN matches NaN whatever its sign or payload. The Rust results are the same
//! on every platform up to the sign and payload of a NaN (x86-64 produces
//! negative default NaNs, AArch64 positive ones): the `libm` crate is pure Rust
//! and Rust never contracts floating-point operations. The support document
//! records the measured differences and why bit identity with the C++ is not
//! reached.
//!
//! The fits run through the shared driver `trace_fitter::optimize`, whose
//! solver is not yet bit-faithful to Eigen beyond recorded fixtures such as
//! these (`docs/TRACE_FITTER_SUPPORT.md`, "Known gap: solver fidelity beyond
//! the fixtures").

// The class-test literals are transcribed verbatim, including digits beyond
// the precision of `f32`, so the values match the C++ literals exactly.
#![allow(clippy::excessive_precision)]

use std::collections::BTreeMap;

use openms::Error;
use openms::analysis::feature_finder_picked::egh_trace_fitter::{EGHTraceFitter, EGHTraceFunctor};
use openms::analysis::feature_finder_picked::helper_structs::{MassTrace, MassTraces, TracePeak};
use openms::analysis::feature_finder_picked::trace_fitter::{TraceFitter, TraceFitterParams};

/// The C2 EGH records, extracted by `../oracle/egh-trace-fitter/extract_c2.py`.
const C2_ORACLE: &str = include_str!("data/egh_trace_fitter/egh_trace_fitter_c2.tsv");

/// Output of `../oracle/egh-trace-fitter/driver.cpp`.
const B5_ORACLE: &str = include_str!("data/egh_trace_fitter/egh_trace_fitter_oracle.tsv");

// ---------------------------------------------------------------------------
// Class test
// ---------------------------------------------------------------------------

/// Retention times of both class-test traces.
const CT_RT: [f64; 21] = [
    677.1, 677.4, 677.7, 678.0, 678.3, 678.6, 678.9, 679.2, 679.5, 679.8, 680.1, 680.4, 680.7,
    681.0, 681.3, 681.6, 681.9, 682.2, 682.5, 682.8, 683.1,
];

/// `p1_1` to `p1_21`, m/z 1000.
const CT_MT1: [f32; 21] = [
    1.08268226589,
    1.58318959267,
    2.22429840363,
    3.00248879081,
    3.89401804768,
    4.8522452777,
    5.80919229659,
    6.68216169129,
    7.38493077109,
    7.84158938645,
    8.0,
    7.84158938645,
    7.38493077109,
    6.68216169129,
    5.80919229659,
    4.8522452777,
    3.89401804768,
    3.00248879081,
    2.22429840363,
    1.58318959267,
    1.08268226589,
];

/// `p2_1` to `p2_21`, m/z 1001.
const CT_MT2: [f32; 21] = [
    0.270670566473,
    0.395797398167,
    0.556074600906,
    0.750622197703,
    0.97350451192,
    1.21306131943,
    1.45229807415,
    1.67054042282,
    1.84623269277,
    1.96039734661,
    2.0,
    1.96039734661,
    1.84623269277,
    1.67054042282,
    1.45229807415,
    1.21306131943,
    0.97350451192,
    0.750622197703,
    0.556074600906,
    0.395797398167,
    0.270670566473,
];

const EXPECTED_SIGMA: f64 = 1.5;
const EXPECTED_H: f64 = 10.0;
const EXPECTED_X0: f64 = 680.1;
const EXPECTED_TAU: f64 = 0.0;

/// ClassTest `isRealSimilar` with the default `absdiff_max_allowed` 1e-5 and
/// `ratio_max_allowed` 1 + 1e-5 (ClassTest.cpp:35-38 and 364-470). The source
/// evaluates in `long double`; the transcribed cases are far from both limits.
fn is_real_similar(number_1: f64, number_2: f64) -> bool {
    const ABSDIFF_MAX_ALLOWED: f64 = 1e-5;
    const RATIO_MAX_ALLOWED: f64 = 1.0 + 1e-5;
    if number_1.is_nan() || number_2.is_nan() {
        return false;
    }
    let is_absdiff_small = (number_1 - number_2).abs() <= ABSDIFF_MAX_ALLOWED;
    if number_1 == 0.0 || number_2 == 0.0 {
        return (number_1 == 0.0 && number_2 == 0.0) || is_absdiff_small;
    }
    let mut ratio = number_1 / number_2;
    if ratio < 0.0 {
        return is_absdiff_small;
    }
    if ratio < 1.0 {
        ratio = 1.0 / ratio;
    }
    ratio <= RATIO_MAX_ALLOWED || is_absdiff_small
}

fn assert_real_similar(actual: f64, expected: f64, what: &str) {
    assert!(
        is_real_similar(actual, expected),
        "{what}: {actual} is not similar to {expected}"
    );
}

/// One class-test trace with `update_maximum` applied, as the file-scope setup
/// does before pushing it.
fn class_trace(theoretical_int: f64, mz: f64, intensities: &[f32; 21], trace: usize) -> MassTrace {
    let mut mt = MassTrace {
        theoretical_int,
        ..MassTrace::default()
    };
    for (k, (&rt, &intensity)) in CT_RT.iter().zip(intensities).enumerate() {
        mt.peaks.push(TracePeak::new(trace, k, rt, mz, intensity));
    }
    mt.update_maximum();
    mt
}

/// `mts` of the file-scope setup: `mt1` (0.8) and `mt2` (0.2), baseline 0,
/// `max_trace` 0.
fn class_mts() -> MassTraces {
    let mut mts = MassTraces::new();
    mts.push(class_trace(0.8, 1000.0, &CT_MT1, 1));
    mts.push(class_trace(0.2, 1001.0, &CT_MT2, 2));
    mts.baseline = 0.0;
    mts.max_trace = 0;
    mts
}

/// `egh_trace_fitter` of the file-scope setup: `max_iteration` 500 set on an
/// otherwise default `Param`, then fitted to `mts`.
fn class_fitter() -> EGHTraceFitter {
    let mut fitter = EGHTraceFitter::new();
    fitter.set_parameters(TraceFitterParams {
        max_iteration: 500,
        weighted: false,
    });
    fitter
        .fit(&class_mts())
        .expect("the class-test fit succeeds");
    fitter
}

#[test]
fn section_default_constructor() {
    // TEST_NOT_EQUAL(ptr, nullPointer): construction succeeds; the defaults are
    // the source's.
    let fitter = EGHTraceFitter::new();
    assert_eq!(
        fitter.parameters(),
        &TraceFitterParams {
            max_iteration: 500,
            weighted: false,
        }
    );
    assert_eq!(EGHTraceFitter::default(), fitter);
}

#[test]
fn section_destructor() {
    let fitter = Box::new(EGHTraceFitter::new());
    drop(fitter);
}

#[test]
fn section_copy_constructor() {
    let egh_trace_fitter = class_fitter();
    let egh1 = egh_trace_fitter.clone();
    assert_eq!(egh1.center(), egh_trace_fitter.center());
    assert_eq!(egh1.height(), egh_trace_fitter.height());
    assert_eq!(egh1.lower_rt_bound(), egh_trace_fitter.lower_rt_bound());
    assert_eq!(egh1.upper_rt_bound(), egh_trace_fitter.upper_rt_bound());
}

#[test]
fn section_assignment_operator() {
    let egh_trace_fitter = class_fitter();
    let mut egh1 = EGHTraceFitter::new();
    assert_eq!(egh1.height(), 0.0);
    egh1.clone_from(&egh_trace_fitter);
    assert_eq!(egh1.center(), egh_trace_fitter.center());
    assert_eq!(egh1.height(), egh_trace_fitter.height());
    assert_eq!(egh1.lower_rt_bound(), egh_trace_fitter.lower_rt_bound());
    assert_eq!(egh1.upper_rt_bound(), egh_trace_fitter.upper_rt_bound());
}

#[test]
fn section_fit() {
    // fit was already done before
    let egh_trace_fitter = class_fitter();
    assert_real_similar(egh_trace_fitter.center(), EXPECTED_X0, "center");
    assert_real_similar(egh_trace_fitter.height(), EXPECTED_H, "height");

    let mut mts = class_mts();
    let mut weighted_fitter = EGHTraceFitter::new();
    let mut params = *weighted_fitter.parameters();
    params.weighted = true;
    weighted_fitter.set_parameters(params);
    weighted_fitter.fit(&mts).expect("weighted fit");
    assert_real_similar(weighted_fitter.center(), EXPECTED_X0, "weighted center");
    assert_real_similar(weighted_fitter.height(), EXPECTED_H, "weighted height");

    mts[0].theoretical_int = 0.4;
    mts[1].theoretical_int = 0.6;
    weighted_fitter.fit(&mts).expect("weighted fit 0.4/0.6");
    assert_real_similar(
        weighted_fitter.center(),
        EXPECTED_X0,
        "weighted 0.4/0.6 center",
    );
    assert_real_similar(weighted_fitter.height(), 6.0825, "weighted 0.4/0.6 height");
}

#[test]
fn section_get_lower_rt_bound() {
    assert_real_similar(
        class_fitter().lower_rt_bound(),
        EXPECTED_X0 - 2.5 * EXPECTED_SIGMA,
        "lower bound",
    );
}

#[test]
fn section_get_upper_rt_bound() {
    assert_real_similar(
        class_fitter().upper_rt_bound(),
        EXPECTED_X0 + 2.5 * EXPECTED_SIGMA,
        "upper bound",
    );
}

#[test]
fn section_get_height() {
    assert_real_similar(class_fitter().height(), EXPECTED_H, "height");
}

#[test]
fn section_get_center() {
    assert_real_similar(class_fitter().center(), EXPECTED_X0, "center");
}

#[test]
fn section_get_tau() {
    assert_real_similar(class_fitter().tau(), EXPECTED_TAU, "tau");
}

#[test]
fn section_get_sigma() {
    assert_real_similar(class_fitter().sigma(), EXPECTED_SIGMA, "sigma");
}

#[test]
fn section_get_value() {
    assert_real_similar(class_fitter().value(EXPECTED_X0), EXPECTED_H, "value");
}

#[test]
fn section_compute_theoretical() {
    let egh_trace_fitter = class_fitter();
    let mut mt = MassTrace {
        theoretical_int: 0.8,
        ..MassTrace::default()
    };
    mt.peaks.push(TracePeak::new(0, 0, EXPECTED_X0, 0.0, 8.0));
    // theoretical should be expected_H * theoretical_int at position expected_x0
    assert_real_similar(
        egh_trace_fitter
            .compute_theoretical(&mt, 0)
            .expect("index 0 is in range"),
        mt.theoretical_int * EXPECTED_H,
        "computeTheoretical",
    );
}

#[test]
fn section_check_maximal_rt_span() {
    let egh_trace_fitter = class_fitter();
    // Maximum RT span in relation to extended area that the model is allowed to have
    // 5.0 * sigma_ > max_rt_span * region_rt_span_
    let mt1 = class_trace(0.8, 1000.0, &CT_MT1, 1);
    let region_rt_span = mt1.peaks[mt1.peaks.len() - 1].rt - mt1.peaks[0].rt;
    let mut max_rt_span = 5.0 * EXPECTED_SIGMA / region_rt_span;
    assert!(!egh_trace_fitter.check_maximal_rt_span(max_rt_span));
    max_rt_span -= 0.1; // accept only smaller regions
    assert!(egh_trace_fitter.check_maximal_rt_span(max_rt_span));
}

#[test]
fn section_check_minimal_rt_span() {
    let egh_trace_fitter = class_fitter();
    // (rt_bounds.second-rt_bounds.first) < min_rt_span * 5.0 * sigma_;
    let rt_bounds = (0.0, 4.0);
    let mut min_rt_span = 0.5;
    assert!(!egh_trace_fitter.check_minimal_rt_span(rt_bounds, min_rt_span));
    min_rt_span += 0.5;
    assert!(egh_trace_fitter.check_minimal_rt_span(rt_bounds, min_rt_span));
}

#[test]
fn section_get_area() {
    assert_real_similar(
        class_fitter().area(),
        (2.0 * std::f64::consts::PI).sqrt() * EXPECTED_SIGMA * EXPECTED_H,
        "area",
    );
}

#[test]
fn section_get_gnuplot_formula() {
    let egh_trace_fitter = class_fitter();
    // The fit section has already set mts[0].theoretical_int to 0.4.
    let mut mts = class_mts();
    mts[0].theoretical_int = 0.4;
    mts[1].theoretical_int = 0.6;
    let formula = egh_trace_fitter.gnuplot_formula(&mts[0], 'f', 0.0, 0.0);
    // should look like -- f(x)= 0 + (((4.5 + 3.93096e-15 * (x - 680.1 )) > 0) ? 8 * exp(-1 * (x - 680.1)**2 / ( 4.5 + 3.93096e-15 * (x - 680.1 ))) : 0) --
    assert!(formula.starts_with("f(x)= 0 + ((("), "{formula}");
    assert!(formula.contains(" )) > 0) ? "), "{formula}");
    assert!(formula.contains(" * exp(-1 * ("), "{formula}");
    assert!(formula.contains(")**2 / ( "), "{formula}");
    assert!(formula.ends_with(" ))) : 0)"), "{formula}");
}

#[test]
fn section_get_fwhm() {
    assert_real_similar(class_fitter().fwhm(), 3.53223007592464, "FWHM");
}

// ---------------------------------------------------------------------------
// Oracle fixtures
// ---------------------------------------------------------------------------

/// Rows of one oracle fixture, keyed by case and quantity.
struct Oracle {
    rows: BTreeMap<(String, String), String>,
    cases: Vec<String>,
}

impl Oracle {
    fn parse(text: &str) -> Self {
        let mut rows = BTreeMap::new();
        let mut cases: Vec<String> = Vec::new();
        for line in text.lines() {
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            let mut parts = line.splitn(3, '\t');
            let (Some(case), Some(quantity), Some(value)) =
                (parts.next(), parts.next(), parts.next())
            else {
                panic!("malformed oracle row: {line}");
            };
            if cases.last().map(String::as_str) != Some(case) && !cases.iter().any(|c| c == case) {
                cases.push(case.to_owned());
            }
            let previous = rows.insert((case.to_owned(), quantity.to_owned()), value.to_owned());
            assert!(previous.is_none(), "duplicate oracle row {case} {quantity}");
        }
        Self { rows, cases }
    }

    fn has(&self, case: &str, quantity: &str) -> bool {
        self.rows
            .contains_key(&(case.to_owned(), quantity.to_owned()))
    }

    fn raw(&self, case: &str, quantity: &str) -> &str {
        self.rows
            .get(&(case.to_owned(), quantity.to_owned()))
            .unwrap_or_else(|| panic!("missing oracle row {case} {quantity}"))
    }

    fn tagged(&self, case: &str, quantity: &str, tag: &str) -> &str {
        let raw = self.raw(case, quantity);
        raw.strip_prefix(tag)
            .unwrap_or_else(|| panic!("oracle row {case} {quantity} is not {tag}: {raw}"))
    }

    fn f64(&self, case: &str, quantity: &str) -> f64 {
        bits64(self.tagged(case, quantity, "f:"))
    }

    fn f64s(&self, case: &str, quantity: &str) -> Vec<f64> {
        let list = self.tagged(case, quantity, "F:");
        if list.is_empty() {
            return Vec::new();
        }
        list.split(',').map(bits64).collect()
    }

    fn f32s(&self, case: &str, quantity: &str) -> Vec<f32> {
        let list = self.tagged(case, quantity, "G:");
        if list.is_empty() {
            return Vec::new();
        }
        list.split(',')
            .map(|hex| f32::from_bits(u32::from_str_radix(hex, 16).expect("f32 bit pattern")))
            .collect()
    }

    fn i64(&self, case: &str, quantity: &str) -> i64 {
        self.tagged(case, quantity, "i:").parse().expect("integer")
    }

    fn bool(&self, case: &str, quantity: &str) -> bool {
        match self.tagged(case, quantity, "b:") {
            "true" => true,
            "false" => false,
            other => panic!("oracle row {case} {quantity} is not a bool: {other}"),
        }
    }

    fn text(&self, case: &str, quantity: &str) -> &str {
        self.tagged(case, quantity, "s:")
    }

    /// `(name, what)` of a recorded OpenMS exception.
    fn exception(&self, case: &str, quantity: &str) -> (&str, &str) {
        self.tagged(case, quantity, "e:")
            .split_once('|')
            .expect("exception name and message")
    }

    /// Cases that record mass traces.
    fn trace_cases(&self) -> Vec<&str> {
        self.cases
            .iter()
            .map(String::as_str)
            .filter(|case| self.has(case, "traces.count"))
            .collect()
    }

    /// The recorded mass traces of `case`, with `update_maximum` applied as the
    /// drivers do.
    fn traces(&self, case: &str) -> MassTraces {
        let count = usize::try_from(self.i64(case, "traces.count")).expect("trace count");
        let mut traces = MassTraces::new();
        for t in 0..count {
            let prefix = format!("trace.{t}.");
            let rt = self.f64s(case, &format!("{prefix}rt"));
            let mz = self.f64s(case, &format!("{prefix}mz"));
            let intensity = self.f32s(case, &format!("{prefix}intensity"));
            assert_eq!(rt.len(), mz.len());
            assert_eq!(rt.len(), intensity.len());
            let mut trace = MassTrace {
                theoretical_int: self.f64(case, &format!("{prefix}theoretical_int")),
                ..MassTrace::default()
            };
            for k in 0..rt.len() {
                trace
                    .peaks
                    .push(TracePeak::new(t, k, rt[k], mz[k], intensity[k]));
            }
            trace.update_maximum();
            traces.push(trace);
        }
        traces.baseline = self.f64(case, "traces.baseline");
        traces.max_trace = usize::try_from(self.i64(case, "traces.max_trace")).expect("max_trace");
        traces
    }
}

fn bits64(hex: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(hex, 16).expect("f64 bit pattern"))
}

fn both() -> [(&'static str, Oracle); 2] {
    [
        ("C2", Oracle::parse(C2_ORACLE)),
        ("B5", Oracle::parse(B5_ORACLE)),
    ]
}

/// Relative tolerance for a single evaluation on inputs recorded bit for bit:
/// functor rows, start points and the queries of a model set to the oracle's
/// own parameters.
///
/// These evaluate the same expressions in the same order as the C++, so they
/// differ only where a `libm` crate function returns a neighbouring double of
/// the platform `exp`, `log` or `atan` the product SDK called (the `libm` crate's
/// `exp` is not correctly rounded; Apple's was on 3594 of the 3600 functor
/// inputs). Such a difference is one or two units in the last place of the
/// affected term, a relative 5e-16 at most in the measurements.
const DIRECT_RELATIVE: f64 = 1e-14;

/// Relative tolerance for a Levenberg-Marquardt result against the C++ fit:
/// the bundle plan's 1e-9 for fits compared off the oracle's own platform.
///
/// Unit-in-the-last-place differences in the residuals steer each iteration
/// slightly differently. The largest measured difference is 2.0e-11, `tau` on
/// the `apex_at_last_scan` budget sweep at budget 240; after all 500
/// evaluations that case differs by at most 1.4e-11 in its parameters (`tau`)
/// and 1.5e-11 in a derived quantity (the area).
const FIT_RELATIVE: f64 = 1e-9;

/// Comparison of results against the oracle, collecting every mismatch before
/// failing and counting how many results agree bit for bit.
#[derive(Default)]
struct Checker {
    compared: usize,
    bit_identical: usize,
    mismatches: Vec<String>,
}

impl Checker {
    /// Records `actual` against `expected`: equal bits, NaN against NaN
    /// (whatever the sign or payload), equal infinities, or finite values with
    /// `|actual - expected| <= relative * max(|actual|, |expected|, scale)`.
    fn close(&mut self, what: &str, actual: f64, expected: f64, relative: f64, scale: f64) {
        self.compared += 1;
        if actual.to_bits() == expected.to_bits() || (actual.is_nan() && expected.is_nan()) {
            self.bit_identical += 1;
            return;
        }
        let magnitude = actual.abs().max(expected.abs()).max(scale.abs());
        let difference = (actual - expected).abs();
        if actual.is_finite() && expected.is_finite() && difference <= relative * magnitude {
            return;
        }
        self.mismatches.push(format!(
            "{what}: got {actual:e} ({:016x}), expected {expected:e} ({:016x}), \
             difference {difference:e}, allowed {:e}",
            actual.to_bits(),
            expected.to_bits(),
            relative * magnitude
        ));
    }

    /// A single evaluation on recorded inputs; see [`DIRECT_RELATIVE`].
    fn direct(&mut self, what: &str, actual: f64, expected: f64) {
        self.close(what, actual, expected, DIRECT_RELATIVE, 0.0);
    }

    fn direct_all(&mut self, what: &str, actual: &[f64], expected: &[f64]) {
        if self.same_length(what, actual.len(), expected.len()) {
            for (i, (&a, &e)) in actual.iter().zip(expected).enumerate() {
                self.direct(&format!("{what}[{i}]"), a, e);
            }
        }
    }

    /// Fitted parameters `[H, t_R, sigma, tau]` against the C++ fit; see
    /// [`FIT_RELATIVE`].
    ///
    /// `tau` is compared on its own scale, except where both values are
    /// rounding noise, below [`FIT_RELATIVE`] times `sigma`: there it is
    /// compared on the scale of `sigma`. The model sees `tau` only through
    /// `tau * t` next to `2 sigma^2`, and a symmetric peak, such as the class
    /// test's, leaves `tau` at rounding noise of about `4e-15` whose relative
    /// value carries no information. A `tau` above that floor gets the plain
    /// relative bound.
    fn fit_params(&mut self, what: &str, actual: &[f64], expected: &[f64]) {
        if !self.same_length(what, actual.len(), expected.len()) || actual.len() != 4 {
            return;
        }
        for (i, name) in ["height", "center", "sigma", "tau"].into_iter().enumerate() {
            let scale = if i == 3 {
                let sigma = actual[2].abs().max(expected[2].abs());
                let tau = actual[3].abs().max(expected[3].abs());
                if tau < FIT_RELATIVE * sigma {
                    sigma
                } else {
                    0.0
                }
            } else {
                0.0
            };
            self.close(
                &format!("{what}.{name}"),
                actual[i],
                expected[i],
                FIT_RELATIVE,
                scale,
            );
        }
    }

    /// A quantity derived from a fitted model against the C++ fit's.
    fn fitted(&mut self, what: &str, actual: f64, expected: f64) {
        self.close(what, actual, expected, FIT_RELATIVE, 0.0);
    }

    fn same_length(&mut self, what: &str, actual: usize, expected: usize) -> bool {
        if actual == expected {
            return true;
        }
        self.compared += 1;
        self.mismatches
            .push(format!("{what}: length {actual} != {expected}"));
        false
    }

    fn eq<T: PartialEq + std::fmt::Debug>(&mut self, what: &str, actual: T, expected: T) {
        self.compared += 1;
        if actual == expected {
            self.bit_identical += 1;
        } else {
            self.mismatches
                .push(format!("{what}: got {actual:?}, expected {expected:?}"));
        }
    }

    fn finish(self, label: &str) {
        assert!(self.compared > 0, "{label}: nothing compared");
        // Shown with --nocapture; the support document records these counts.
        eprintln!(
            "{label}: {} comparisons, {} bit-identical or exactly equal",
            self.compared, self.bit_identical
        );
        assert!(
            self.mismatches.is_empty(),
            "{label}: {} of {} comparisons differ:\n{}",
            self.mismatches.len(),
            self.compared,
            self.mismatches.join("\n")
        );
    }
}

fn fitter_for(max_iteration: i64, weighted: bool) -> EGHTraceFitter {
    EGHTraceFitter::with_parameters(TraceFitterParams {
        max_iteration,
        weighted,
    })
}

fn fitted_params(fitter: &EGHTraceFitter) -> [f64; 4] {
    [
        fitter.height(),
        fitter.center(),
        fitter.sigma(),
        fitter.tau(),
    ]
}

#[test]
fn oracle_defaults_match_the_source() {
    let oracle = Oracle::parse(B5_ORACLE);
    let fitter = EGHTraceFitter::new();
    assert_eq!(
        fitter.parameters().max_iteration,
        oracle.i64("defaults", "max_iteration")
    );
    assert_eq!(
        if fitter.parameters().weighted {
            "true"
        } else {
            "false"
        },
        oracle.text("defaults", "weighted")
    );
}

#[test]
fn oracle_initial_parameters() {
    let mut check = Checker::default();
    for (label, oracle) in both() {
        for case in oracle.trace_cases() {
            let traces = oracle.traces(case);
            let init = EGHTraceFitter::initial_parameters(&traces).expect("initial parameters");
            let what = format!("{label} {case} init");
            for (name, actual) in [
                ("height", init.height),
                ("apex_rt", init.apex_rt),
                ("sigma", init.sigma),
                ("tau", init.tau),
                ("region_rt_span", init.region_rt_span),
            ] {
                check.direct(
                    &format!("{what}.{name}"),
                    actual,
                    oracle.f64(case, &format!("init.{name}")),
                );
            }
        }
    }
    check.finish("initial parameters");
}

#[test]
fn oracle_functor_residuals_and_jacobians() {
    let mut check = Checker::default();
    for (label, oracle) in both() {
        for case in oracle.trace_cases() {
            let traces = oracle.traces(case);
            let weighted = oracle.bool(case, "weighted");
            let functor = EGHTraceFunctor::new(&traces, weighted);
            // |I * w| per row: a residual is the difference of the model term and
            // this, so its rounding is measured on the terms' scale.
            let observed: Vec<f64> = traces
                .iter()
                .flat_map(|trace| {
                    let weight = if weighted { trace.theoretical_int } else { 1.0 };
                    trace
                        .peaks
                        .iter()
                        .map(move |peak| (f64::from(peak.intensity) * weight).abs())
                })
                .collect();
            let prefixes: Vec<String> = oracle
                .rows
                .keys()
                .filter(|(c, q)| c == case && q.starts_with("functor.") && q.ends_with(".x"))
                .filter_map(|(_, q)| q.strip_suffix('x').map(str::to_owned))
                .collect();
            for prefix in prefixes {
                let what = format!("{label} {case} {prefix}");
                let x = oracle.f64s(case, &format!("{prefix}x"));
                check.eq(
                    &format!("{what}values"),
                    i64::try_from(functor.values()).expect("values"),
                    oracle.i64(case, &format!("{prefix}values")),
                );
                check.eq(
                    &format!("{what}inputs"),
                    i64::try_from(functor.inputs()).expect("inputs"),
                    oracle.i64(case, &format!("{prefix}inputs")),
                );
                check.eq(
                    &format!("{what}residual_return"),
                    0,
                    oracle.i64(case, &format!("{prefix}residual_return")),
                );
                check.eq(
                    &format!("{what}jacobian_return"),
                    0,
                    oracle.i64(case, &format!("{prefix}jacobian_return")),
                );
                let mut fvec = vec![0.0; functor.values()];
                functor.residuals(&x, &mut fvec).expect("residuals");
                let expected = oracle.f64s(case, &format!("{prefix}residuals"));
                if check.same_length(&format!("{what}residuals"), fvec.len(), expected.len()) {
                    for (row, ((&a, &e), &scale)) in
                        fvec.iter().zip(&expected).zip(&observed).enumerate()
                    {
                        check.close(
                            &format!("{what}residuals[{row}]"),
                            a,
                            e,
                            DIRECT_RELATIVE,
                            scale,
                        );
                    }
                }
                let mut jacobian = vec![0.0; functor.values() * functor.inputs()];
                functor.jacobian(&x, &mut jacobian).expect("jacobian");
                check.direct_all(
                    &format!("{what}jacobian"),
                    &jacobian,
                    &oracle.f64s(case, &format!("{prefix}jacobian_column_major")),
                );
            }
        }
    }
    check.finish("functor");
}

/// The queries recorded under `prefix`, evaluated on `model`, which holds the
/// oracle's own parameters, so every comparison is a single evaluation on
/// recorded inputs. Formulas are compared as text.
fn check_model_queries(
    check: &mut Checker,
    oracle: &Oracle,
    case: &str,
    prefix: &str,
    what: &str,
    model: &EGHTraceFitter,
    traces: Option<&MassTraces>,
) {
    let q = |name: &str| format!("{prefix}{name}");
    let has = |name: &str| oracle.has(case, &q(name));
    for (name, actual) in [
        ("lower_rt_bound", model.lower_rt_bound()),
        ("upper_rt_bound", model.upper_rt_bound()),
        ("sigma_5_bound_first", model.lower_rt_bound()),
        ("sigma_5_bound_second", model.upper_rt_bound()),
        ("fwhm", model.fwhm()),
        ("area", model.area()),
        ("alpha_0_5_first", model.alpha_boundaries(0.5).0),
        ("alpha_0_5_second", model.alpha_boundaries(0.5).1),
    ] {
        if has(name) {
            check.direct(
                &format!("{what}.{name}"),
                actual,
                oracle.f64(case, &q(name)),
            );
        }
    }
    if has("value_rts") {
        let rts = oracle.f64s(case, &q("value_rts"));
        let values: Vec<f64> = rts.iter().map(|&rt| model.value(rt)).collect();
        check.direct_all(
            &format!("{what}.values"),
            &values,
            &oracle.f64s(case, &q("values")),
        );
    }
    let mut i = 0;
    while has(&format!("check_min.{i}.result")) {
        let p = format!("check_min.{i}.");
        let bounds = (
            oracle.f64(case, &q(&format!("{p}lower"))),
            oracle.f64(case, &q(&format!("{p}upper"))),
        );
        let span = oracle.f64(case, &q(&format!("{p}min_rt_span")));
        check.eq(
            &format!("{what}.{p}result"),
            model.check_minimal_rt_span(bounds, span),
            oracle.bool(case, &q(&format!("{p}result"))),
        );
        i += 1;
    }
    let mut i = 0;
    while has(&format!("check_max.{i}.result")) {
        let p = format!("check_max.{i}.");
        let span = oracle.f64(case, &q(&format!("{p}max_rt_span")));
        check.eq(
            &format!("{what}.{p}result"),
            model.check_maximal_rt_span(span),
            oracle.bool(case, &q(&format!("{p}result"))),
        );
        i += 1;
    }
    let Some(traces) = traces else {
        return;
    };
    for (t, trace) in traces.iter().enumerate() {
        let p = format!("trace.{t}.");
        if !has(&format!("{p}compute_theoretical")) {
            continue;
        }
        let theoretical: Vec<f64> = (0..trace.peaks.len())
            .map(|k| model.compute_theoretical(trace, k).expect("peak index"))
            .collect();
        check.direct_all(
            &format!("{what}.{p}compute_theoretical"),
            &theoretical,
            &oracle.f64s(case, &q(&format!("{p}compute_theoretical"))),
        );
        let name = oracle
            .text(case, &q(&format!("{p}gnuplot_name")))
            .chars()
            .next()
            .expect("name");
        let formula = model.gnuplot_formula(
            trace,
            name,
            oracle.f64(case, &q(&format!("{p}gnuplot_baseline"))),
            oracle.f64(case, &q(&format!("{p}gnuplot_rt_shift"))),
        );
        check.eq(
            &format!("{what}.{p}gnuplot_formula"),
            formula.as_str(),
            oracle.text(case, &q(&format!("{p}gnuplot_formula"))),
        );
    }
    if has("gnuplot_f_0_0") {
        let formula = model.gnuplot_formula(&traces[0], 'f', 0.0, 0.0);
        check.eq(
            &format!("{what}.gnuplot_f_0_0"),
            formula.as_str(),
            oracle.text(case, &q("gnuplot_f_0_0")),
        );
    }
}

/// `UnableToFit-FinalSet: <what>`, the message of a source `UnableToFit`.
fn expect_unable_to_fit(
    check: &mut Checker,
    what: &str,
    result: openms::Result<()>,
    name: &str,
    message: &str,
) {
    match result {
        Err(Error::InvalidValue(text)) => {
            check.eq(what, text, format!("{name}: {message}"));
        }
        other => check.eq(
            what,
            format!("{other:?}"),
            format!("Err(InvalidValue({name}: {message}))"),
        ),
    }
}

/// The C++ fit's `[H, t_R, sigma, tau]` of `case`.
fn oracle_fit_params(oracle: &Oracle, case: &str) -> Vec<f64> {
    if oracle.has(case, "fit.params") {
        return oracle.f64s(case, "fit.params");
    }
    ["height", "center", "sigma", "tau"]
        .iter()
        .map(|name| oracle.f64(case, &format!("fit.{name}")))
        .collect()
}

#[test]
fn oracle_fit_and_queries() {
    let mut check = Checker::default();
    for (label, oracle) in both() {
        for case in oracle.trace_cases() {
            let traces = oracle.traces(case);
            let what = format!("{label} {case} fit");
            let mut fitter = fitter_for(
                oracle.i64(case, "fit.max_iteration"),
                oracle.bool(case, "weighted"),
            );
            let before = fitter.clone();
            let result = fitter.fit(&traces);
            if oracle.bool(case, "fit.threw") {
                let (name, message) = oracle.exception(case, "fit.exception");
                expect_unable_to_fit(&mut check, &what, result, name, message);
                check.eq(
                    &format!("{what} leaves the fitter unchanged"),
                    &fitter,
                    &before,
                );
                continue;
            }
            if let Err(error) = result {
                check.eq(&what, format!("{error:?}"), "Ok".to_owned());
                continue;
            }
            let expected = oracle_fit_params(&oracle, case);
            check.fit_params(&what, &fitted_params(&fitter), &expected);
            if oracle.has(case, "fit.member_max_iterations") {
                check.eq(
                    &format!("{what}.member_max_iterations"),
                    fitter.parameters().max_iteration,
                    oracle.i64(case, "fit.member_max_iterations"),
                );
                check.eq(
                    &format!("{what}.member_weighted"),
                    fitter.parameters().weighted,
                    oracle.bool(case, "fit.member_weighted"),
                );
            }
            // The fitted model's own derived quantities, against the C++ fit's.
            for (name, actual) in [
                ("lower_rt_bound", fitter.lower_rt_bound()),
                ("upper_rt_bound", fitter.upper_rt_bound()),
                ("fwhm", fitter.fwhm()),
                ("area", fitter.area()),
            ] {
                if oracle.has(case, &format!("fit.{name}")) {
                    check.fitted(
                        &format!("{what}.{name} of the Rust fit"),
                        actual,
                        oracle.f64(case, &format!("fit.{name}")),
                    );
                }
            }
            for (t, trace) in traces.iter().enumerate() {
                let key = format!("fit.trace.{t}.compute_theoretical");
                if !oracle.has(case, &key) {
                    continue;
                }
                for (k, &e) in oracle.f64s(case, &key).iter().enumerate() {
                    let actual = fitter.compute_theoretical(trace, k).expect("peak index");
                    check.fitted(
                        &format!("{what}.trace.{t}.compute_theoretical[{k}] of the Rust fit"),
                        actual,
                        e,
                    );
                }
            }
            // Every query on the C++ fit's parameters, with the region span of
            // the Rust fit, which is a plain difference of recorded times.
            let mut model = fitter.clone();
            let oracle_x: [f64; 4] = expected.try_into().expect("four parameters");
            model.set_optimized_parameters(oracle_x);
            check_model_queries(
                &mut check,
                &oracle,
                case,
                "fit.",
                &what,
                &model,
                Some(&traces),
            );
        }
    }
    check.finish("fit");
}

#[test]
fn oracle_budget_boundaries() {
    let mut check = Checker::default();
    for (label, oracle) in both() {
        for case in oracle.trace_cases() {
            if !oracle.has(case, "sweep.n_max") {
                continue;
            }
            let traces = oracle.traces(case);
            let weighted = oracle.bool(case, "weighted");
            let n_max = oracle.i64(case, "sweep.n_max");
            let smallest = oracle.i64(case, "sweep.smallest_identical");
            let fit = |n: i64| {
                let mut fitter = fitter_for(n, weighted);
                fitter.fit(&traces).map(|()| fitted_params(&fitter))
            };
            let reference = fit(n_max).expect("fit at n_max");
            for n in 1..=n_max {
                let what = format!("{label} {case} max_iteration {n}");
                let result = fit(n);
                let p = format!("sweep.{n}.");
                if oracle.has(case, &format!("{p}params")) {
                    match &result {
                        Ok(params) => check.fit_params(
                            &what,
                            params,
                            &oracle.f64s(case, &format!("{p}params")),
                        ),
                        Err(error) => check.eq(&what, format!("{error:?}"), "Ok".to_owned()),
                    }
                }
                // The Rust result at n must equal the Rust result at n_max bit
                // for bit exactly where the C++ result at n equals the C++ one
                // at n_max, and differ where it differs.
                let identical = match &result {
                    Ok(params) => params
                        .iter()
                        .zip(&reference)
                        .all(|(a, b)| a.to_bits() == b.to_bits()),
                    Err(_) => false,
                };
                let expected = if n >= smallest {
                    true
                } else {
                    oracle.bool(case, &format!("{p}equals_n_max"))
                };
                check.eq(
                    &format!("{what} equals max_iteration {n_max}"),
                    identical,
                    expected,
                );
            }
        }
    }
    check.finish("budget boundaries");
}

#[test]
fn oracle_parameter_vectors() {
    let mut check = Checker::default();
    for (label, oracle) in both() {
        let cases: Vec<&str> = oracle
            .cases
            .iter()
            .map(String::as_str)
            .filter(|case| oracle.has(case, "optimized.input"))
            .collect();
        for case in cases {
            let what = format!("{label} {case} optimized");
            let input: [f64; 4] = oracle
                .f64s(case, "optimized.input")
                .try_into()
                .expect("four parameters");
            let mut model = EGHTraceFitter::new();
            model.set_optimized_parameters(input);
            // The getters return the vector unchanged.
            for (name, actual) in [
                ("height", model.height()),
                ("center", model.center()),
                ("apex_rt", model.center()),
                ("sigma", model.sigma()),
                ("tau", model.tau()),
            ] {
                if oracle.has(case, &format!("optimized.{name}")) {
                    let expected = oracle.f64(case, &format!("optimized.{name}"));
                    check.eq(
                        &format!("{what}.{name} bits"),
                        actual.to_bits(),
                        expected.to_bits(),
                    );
                }
            }
            check_model_queries(&mut check, &oracle, case, "optimized.", &what, &model, None);
            if oracle.has(case, "optimized.synthetic_trace.rt") {
                let mut trace = MassTrace {
                    theoretical_int: oracle.f64(case, "optimized.synthetic_trace.theoretical_int"),
                    ..MassTrace::default()
                };
                for (k, rt) in oracle
                    .f64s(case, "optimized.synthetic_trace.rt")
                    .into_iter()
                    .enumerate()
                {
                    trace.peaks.push(TracePeak::new(0, k, rt, 500.0, 1.0));
                }
                let theoretical: Vec<f64> = (0..trace.peaks.len())
                    .map(|k| model.compute_theoretical(&trace, k).expect("peak index"))
                    .collect();
                check.direct_all(
                    &format!("{what}.synthetic_trace.compute_theoretical"),
                    &theoretical,
                    &oracle.f64s(case, "optimized.synthetic_trace.compute_theoretical"),
                );
                check.eq(
                    &format!("{what}.gnuplot_h_2.5_17.25"),
                    model.gnuplot_formula(&trace, 'h', 2.5, 17.25).as_str(),
                    oracle.text(case, "optimized.gnuplot_h_2.5_17.25"),
                );
                check.eq(
                    &format!("{what}.gnuplot_f_0_0"),
                    model.gnuplot_formula(&trace, 'f', 0.0, 0.0).as_str(),
                    oracle.text(case, "optimized.gnuplot_f_0_0"),
                );
            }
            if oracle.has(case, "alpha_boundaries.alpha") {
                let alphas = oracle.f64s(case, "alpha_boundaries.alpha");
                let (first, second): (Vec<f64>, Vec<f64>) = alphas
                    .iter()
                    .map(|&alpha| model.alpha_boundaries(alpha))
                    .unzip();
                check.direct_all(
                    &format!("{what}.alpha_boundaries.first"),
                    &first,
                    &oracle.f64s(case, "alpha_boundaries.first"),
                );
                check.direct_all(
                    &format!("{what}.alpha_boundaries.second"),
                    &second,
                    &oracle.f64s(case, "alpha_boundaries.second"),
                );
            }
        }
    }
    check.finish("parameter vectors");
}

/// The package's acceptance values from the C2 oracle, written out in decimal:
/// the class-test fit (unweighted, 0.8/0.2) and the weighted 0.4/0.6 height,
/// within [`FIT_RELATIVE`].
#[test]
fn c2_acceptance_values() {
    let within = |actual: f64, expected: f64, what: &str| {
        assert!(
            (actual - expected).abs() <= FIT_RELATIVE * expected.abs(),
            "{what}: {actual:e} is not within 1e-9 of {expected:e}"
        );
    };
    let fitter = class_fitter();
    within(fitter.height(), 9.9999999742433339, "H");
    within(fitter.sigma(), 1.5000000035584744, "sigma");
    within(fitter.fwhm(), 3.5322300759260088, "FWHM");
    within(fitter.area(), 37.599425992354142, "area");

    let mut mts = class_mts();
    mts[0].theoretical_int = 0.4;
    mts[1].theoretical_int = 0.6;
    let mut weighted = fitter_for(500, true);
    weighted.fit(&mts).expect("weighted 0.4/0.6 fit");
    within(weighted.height(), 6.0824742044607856, "weighted 0.4/0.6 H");
}

// ---------------------------------------------------------------------------
// Native contracts
// ---------------------------------------------------------------------------

fn unable_to_fit_message(result: openms::Result<()>) -> String {
    match result {
        Err(Error::InvalidValue(message)) => message,
        other => panic!("expected Error::InvalidValue, got {other:?}"),
    }
}

#[test]
fn traces_without_peaks_are_unable_to_fit_and_leave_the_fitter_unchanged() {
    let mut fitter = class_fitter();
    let before = fitter.clone();

    let empty = MassTraces::new();
    assert_eq!(
        unable_to_fit_message(fitter.fit(&empty)),
        "UnableToFit-FinalSet: Skipping feature, we always expect N>=p"
    );
    assert_eq!(fitter, before);

    let mut no_peaks = MassTraces::new();
    no_peaks.push(MassTrace::default());
    no_peaks.push(MassTrace::default());
    assert_eq!(
        unable_to_fit_message(fitter.fit(&no_peaks)),
        "UnableToFit-FinalSet: Skipping feature, we always expect N>=p"
    );
    assert_eq!(fitter, before);

    assert!(matches!(
        EGHTraceFitter::initial_parameters(&no_peaks),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn a_failed_fit_keeps_the_previous_model() {
    let mut fitter = class_fitter();
    let before = fitter.clone();
    let mut params = *fitter.parameters();
    params.max_iteration = 0;
    fitter.set_parameters(params);
    assert_eq!(
        unable_to_fit_message(fitter.fit(&class_mts())),
        "UnableToFit-FinalSet: Could not fit the gaussian to the data: Error 0"
    );
    assert_eq!(fitted_params(&fitter), fitted_params(&before));
    assert_eq!(fitter.lower_rt_bound(), before.lower_rt_bound());
    assert_eq!(fitter.upper_rt_bound(), before.upper_rt_bound());
    assert_eq!(fitter.parameters().max_iteration, 0);

    params.max_iteration = i64::MIN;
    fitter.set_parameters(params);
    assert!(fitter.fit(&class_mts()).is_err());
    params.max_iteration = i64::MAX;
    fitter.set_parameters(params);
    fitter.fit(&class_mts()).expect("an unbounded budget fits");
    assert_eq!(fitted_params(&fitter), fitted_params(&before));
}

#[test]
fn a_nan_retention_time_the_profile_cannot_merge_is_an_error() {
    let mut traces = class_mts();
    traces[1].peaks[3].rt = f64::NAN;
    let mut fitter = class_fitter();
    let before = fitter.clone();
    assert!(matches!(fitter.fit(&traces), Err(Error::InvalidValue(_))));
    assert_eq!(fitter, before);
}

#[test]
fn compute_theoretical_refuses_an_out_of_range_peak() {
    let fitter = class_fitter();
    let traces = class_mts();
    assert!(fitter.compute_theoretical(&traces[0], 20).is_ok());
    assert!(matches!(
        fitter.compute_theoretical(&traces[0], 21),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn functor_refuses_mis_sized_buffers() {
    let traces = class_mts();
    let functor = EGHTraceFunctor::new(&traces, false);
    assert_eq!(functor.values(), 42);
    assert_eq!(functor.inputs(), 4);
    let x = [10.0, 680.1, 1.5, 0.0];
    let mut fvec = vec![7.0; 41];
    assert!(functor.residuals(&x, &mut fvec).is_err());
    assert!(fvec.iter().all(|&v| v == 7.0));
    let mut fvec = vec![0.0; 42];
    assert!(functor.residuals(&x[..3], &mut fvec).is_err());
    assert!(functor.residuals(&x, &mut fvec).is_ok());
    let mut jacobian = vec![7.0; 42 * 4 - 1];
    assert!(functor.jacobian(&x, &mut jacobian).is_err());
    assert!(jacobian.iter().all(|&v| v == 7.0));
    let mut jacobian = vec![0.0; 42 * 4];
    assert!(functor.jacobian(&[1.0; 5], &mut jacobian).is_err());
    assert!(functor.jacobian(&x, &mut jacobian).is_ok());
}

#[test]
fn the_fitter_works_through_the_trait_object() {
    let mut fitter: Box<dyn TraceFitter> = Box::new(EGHTraceFitter::new());
    fitter.set_parameters(TraceFitterParams {
        max_iteration: 500,
        weighted: true,
    });
    assert!(fitter.parameters().weighted);
    fitter
        .fit(&class_mts())
        .expect("fit through dyn TraceFitter");
    assert_real_similar(fitter.center(), EXPECTED_X0, "center");
    assert_real_similar(fitter.height(), EXPECTED_H, "height");
}

#[test]
fn a_stored_model_rebuilds_the_fitted_one() {
    let fitted = class_fitter();
    let mut rebuilt = EGHTraceFitter::new();
    rebuilt.set_optimized_parameters(fitted_params(&fitted));
    assert_eq!(rebuilt.lower_rt_bound(), fitted.lower_rt_bound());
    assert_eq!(rebuilt.upper_rt_bound(), fitted.upper_rt_bound());
    assert_eq!(rebuilt.fwhm(), fitted.fwhm());
    assert_eq!(rebuilt.area(), fitted.area());
    // The region span is not part of the parameter vector.
    assert!(!rebuilt.check_maximal_rt_span(f64::INFINITY));
    assert!(rebuilt.check_maximal_rt_span(0.0));
}
