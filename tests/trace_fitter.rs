// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FEATUREFINDER/TraceFitter`: the trait, its parameters and the shared
//! Levenberg-Marquardt driver.
//!
//! `TraceFitter_test.cpp` (core bc9cc12) has 16 `START_SECTION`s. It derives a
//! subclass whose every override throws `Exception::NotImplemented` and checks
//! that each call throws; the copy constructor and assignment are
//! `NOT_TESTABLE`. Here the subclass is `DerivedTraceFitter`, whose fallible
//! methods return `Error::Unsupported` and whose plain queries return NaN and
//! count the call, so each section checks that the trait dispatches to the
//! implementor and that the base supplies nothing. That the trait cannot be
//! implemented partially is checked at compile time by the `compile_fail`
//! doctest on `TraceFitter`.
//!
//! The parameter tests replay `tests/data/gauss_trace_fitter/gauss_extra.tsv`,
//! the executed product-SDK `getDefaults` and `setParameters` results (tier 1).
//! The driver tests compare [`optimize`] with a direct call of the solver the
//! source configures, and check its refusals and the state they leave.

use std::cell::Cell;
use std::collections::BTreeMap;

use openms::analysis::feature_finder_picked::helper_structs::{MassTrace, MassTraces, TracePeak};
use openms::analysis::feature_finder_picked::trace_fitter::{
    FEWER_RESIDUALS_THAN_PARAMETERS, MAX_RESIDUAL_WORK, ProfileSmoothing, TraceFitter,
    TraceFitterParams, UNABLE_TO_FIT_FINAL_SET, compute_theoretical, initial_shape, optimize,
    optimize_with_status, stream_number, unable_to_fit,
};
use openms::math::fitters::levenberg_marquardt::{
    DenseMatrix, LmParameters, LmStatus, MAX_POINTS, minimize,
};
use openms::param::{Param, ParamValue};
use openms::{Error, Result};

const EXTRA: &str = include_str!("data/gauss_trace_fitter/gauss_extra.tsv");

// ---------------------------------------------------------------------------
// START_SECTIONs of TraceFitter_test.cpp
// ---------------------------------------------------------------------------

/// The class test's `DerivedTraceFitter`: nothing is implemented.
#[derive(Default)]
struct DerivedTraceFitter {
    parameters: TraceFitterParams,
    queries: Cell<usize>,
}

impl DerivedTraceFitter {
    fn not_implemented(&self) -> f64 {
        self.queries.set(self.queries.get() + 1);
        f64::NAN
    }
}

impl TraceFitter for DerivedTraceFitter {
    fn parameters(&self) -> &TraceFitterParams {
        &self.parameters
    }
    fn set_parameters(&mut self, parameters: TraceFitterParams) {
        self.parameters = parameters;
    }
    fn fit(&mut self, _traces: &MassTraces) -> Result<()> {
        Err(Error::Unsupported("NotImplemented".into()))
    }
    fn lower_rt_bound(&self) -> f64 {
        self.not_implemented()
    }
    fn upper_rt_bound(&self) -> f64 {
        self.not_implemented()
    }
    fn height(&self) -> f64 {
        self.not_implemented()
    }
    fn center(&self) -> f64 {
        self.not_implemented()
    }
    fn fwhm(&self) -> f64 {
        self.not_implemented()
    }
    fn value(&self, _rt: f64) -> f64 {
        self.not_implemented()
    }
    fn area(&self) -> f64 {
        self.not_implemented()
    }
    fn check_minimal_rt_span(&self, _rt_bounds: (f64, f64), _min_rt_span: f64) -> bool {
        self.not_implemented();
        false
    }
    fn check_maximal_rt_span(&self, _max_rt_span: f64) -> bool {
        self.not_implemented();
        false
    }
    fn gnuplot_formula(&self, _: &MassTrace, _: char, _: f64, _: f64) -> String {
        self.not_implemented();
        String::new()
    }
    fn compute_theoretical(&self, trace: &MassTrace, k: usize) -> Result<f64> {
        compute_theoretical(self, trace, k)
    }
}

#[track_caller]
fn assert_not_implemented(fitter: &DerivedTraceFitter, value: f64, calls: usize) {
    assert!(value.is_nan());
    assert_eq!(fitter.queries.get(), calls);
}

/// `START_SECTION(TraceFitter())` and `START_SECTION(~TraceFitter())`: an
/// implementor is constructed, used through the trait object and dropped.
#[test]
fn section_constructor_and_destructor() {
    let ptr: Box<dyn TraceFitter> = Box::new(DerivedTraceFitter::default());
    assert_eq!(*ptr.parameters(), TraceFitterParams::default());
    drop(ptr);
}

/// `START_SECTION((TraceFitter(const TraceFitter& source)))` and
/// `START_SECTION((virtual TraceFitter& operator=(const TraceFitter& source)))`
/// are `NOT_TESTABLE` in the source ("has no public members to check if copy has
/// same proberties"). The trait has no copy; each fitter derives `Clone`.
#[test]
fn section_copy_and_assignment_are_not_testable() {}

/// `START_SECTION((virtual void fit(...)=0))`.
#[test]
fn section_fit() {
    let mut fitter = DerivedTraceFitter::default();
    let m = MassTraces::new();
    assert!(matches!(fitter.fit(&m), Err(Error::Unsupported(_))));
}

/// `getLowerRTBound`, `getUpperRTBound`, `getHeight`, `getCenter`, `getValue`,
/// `getArea` and `getFWHM`: one section each in the source.
#[test]
fn sections_queries() {
    let fitter = DerivedTraceFitter::default();
    assert_not_implemented(&fitter, fitter.lower_rt_bound(), 1);
    assert_not_implemented(&fitter, fitter.upper_rt_bound(), 2);
    assert_not_implemented(&fitter, fitter.height(), 3);
    assert_not_implemented(&fitter, fitter.center(), 4);
    assert_not_implemented(&fitter, fitter.value(0.0), 5);
    assert_not_implemented(&fitter, fitter.area(), 6);
    assert_not_implemented(&fitter, fitter.fwhm(), 7);
}

/// `START_SECTION((double computeTheoretical(...)))`: the non-virtual member
/// calls `getValue`, which the source's subclass makes throw.
#[test]
fn section_compute_theoretical() {
    let fitter = DerivedTraceFitter::default();
    let mut mt = MassTrace::default();
    mt.peaks.push(TracePeak::new(0, 0, 1.0, 0.0, 0.0));
    assert_not_implemented(&fitter, fitter.compute_theoretical(&mt, 0).unwrap(), 1);
}

/// `START_SECTION((virtual bool checkMinimalRTSpan(...)=0))` and
/// `START_SECTION((virtual bool checkMaximalRTSpan(...)=0))`.
#[test]
fn sections_span_checks() {
    let fitter = DerivedTraceFitter::default();
    let p = (0.0, 0.0);
    let x = 0.0;
    assert!(!fitter.check_minimal_rt_span(p, x));
    assert_eq!(fitter.queries.get(), 1);
    assert!(!fitter.check_maximal_rt_span(x));
    assert_eq!(fitter.queries.get(), 2);
}

/// `START_SECTION((virtual std::string getGnuplotFormula(...)=0))`.
#[test]
fn section_gnuplot_formula() {
    let fitter = DerivedTraceFitter::default();
    let mt = MassTrace::default();
    assert!(fitter.gnuplot_formula(&mt, 'f', 0.0, 0.0).is_empty());
    assert_eq!(fitter.queries.get(), 1);
}

// ---------------------------------------------------------------------------
// Parameters against the executed getDefaults and setParameters
// ---------------------------------------------------------------------------

fn fixture() -> BTreeMap<(String, String), (String, String)> {
    let mut rows = BTreeMap::new();
    for line in EXTRA.lines() {
        let fields: Vec<&str> = line.splitn(4, '\t').collect();
        assert_eq!(fields.len(), 4, "{line:?}");
        rows.insert(
            (fields[0].to_owned(), fields[1].to_owned()),
            (fields[2].to_owned(), fields[3].to_owned()),
        );
    }
    rows
}

#[track_caller]
fn cell<'a>(
    rows: &'a BTreeMap<(String, String), (String, String)>,
    case: &str,
    quantity: &str,
) -> &'a str {
    &rows
        .get(&(case.to_owned(), quantity.to_owned()))
        .unwrap_or_else(|| panic!("no fixture row {case} {quantity}"))
        .1
}

#[test]
fn the_default_record_is_the_source_default() {
    let defaults = TraceFitterParams::default();
    assert_eq!(defaults.max_iteration, 500);
    assert_eq!(
        defaults.max_iteration,
        TraceFitterParams::DEFAULT_MAX_ITERATION
    );
    assert!(!defaults.weighted);
    assert_eq!(TraceFitterParams::HANDLER_NAME, "TraceFitter");
}

#[test]
fn defaults_match_the_executed_get_defaults() {
    let rows = fixture();
    let defaults = TraceFitterParams::defaults().unwrap();
    let size: usize = cell(&rows, "defaults", "size").parse().unwrap();
    assert_eq!(defaults.size(), size);
    let items: Vec<_> = defaults.iter().unwrap().collect();
    assert_eq!(items.len(), size);
    for (index, item) in items.iter().enumerate() {
        let q = |name: &str| cell(&rows, "defaults", &format!("entry[{index}].{name}"));
        assert_eq!(item.key, q("name"));
        let entry = item.entry;
        match (&entry.value, q("type")) {
            (ParamValue::Integer(value), "int") => assert_eq!(value.to_string(), q("value")),
            (ParamValue::String(value), "string") => assert_eq!(value, q("value")),
            (value, kind) => panic!("{value:?} is not a {kind}"),
        }
        assert_eq!(entry.description, q("description"));
        let tags: Vec<&str> = entry.tags.iter().map(String::as_str).collect();
        assert_eq!(tags.join(","), q("tags"));
        assert_eq!(entry.valid_strings.join(","), q("valid_strings"));
        assert_eq!(entry.min_int.to_string(), q("min_int"));
        assert_eq!(entry.max_int.to_string(), q("max_int"));
    }
    // A fresh source fitter holds the defaults in its members.
    let (members, warnings) = TraceFitterParams::from_param(&Param::new()).unwrap();
    assert!(warnings.is_empty());
    assert_eq!(
        members.max_iteration.to_string(),
        cell(&rows, "defaults", "member.max_iterations")
    );
    assert_eq!(
        members.weighted.to_string(),
        cell(&rows, "defaults", "member.weighted")
    );
    assert_eq!(members, TraceFitterParams::default());
}

#[test]
fn set_parameters_matches_the_executed_source() {
    let rows = fixture();
    let cases: [(&str, &str, ParamValue); 8] = [
        (
            "only_max_iteration_7",
            "max_iteration",
            ParamValue::Integer(7),
        ),
        (
            "weighted_true",
            "weighted",
            ParamValue::String("true".into()),
        ),
        (
            "max_iteration_negative",
            "max_iteration",
            ParamValue::Integer(-3),
        ),
        ("unknown_key", "foo", ParamValue::Integer(1)),
        (
            "weighted_invalid",
            "weighted",
            ParamValue::String("maybe".into()),
        ),
        ("weighted_int", "weighted", ParamValue::Integer(1)),
        (
            "max_iteration_float",
            "max_iteration",
            ParamValue::Float(500.0),
        ),
        (
            "max_iteration_string",
            "max_iteration",
            ParamValue::String("500".into()),
        ),
    ];
    for (name, key, value) in cases {
        let case = format!("set_parameters.{name}");
        let mut param = Param::new();
        param.set_value(key, value, "", &[]).unwrap();
        let result = TraceFitterParams::from_param(&param);
        if cell(&rows, &case, "threw") == "true" {
            // The source throws InvalidParameter; the message wording is the
            // Param module's, not the source's.
            assert!(
                matches!(result, Err(Error::InvalidValue(_))),
                "{case}: expected an error, got {result:?}"
            );
            assert_eq!(cell(&rows, &case, "exception.name"), "InvalidParameter");
            continue;
        }
        let (members, warnings) = result.unwrap_or_else(|e| panic!("{case}: {e}"));
        assert_eq!(
            members.max_iteration.to_string(),
            cell(&rows, &case, "member.max_iterations"),
            "{case}"
        );
        assert_eq!(
            members.weighted.to_string(),
            cell(&rows, &case, "member.weighted"),
            "{case}"
        );
        assert_eq!(
            members.max_iteration.to_string(),
            cell(&rows, &case, "parameters.max_iteration"),
            "{case}"
        );
        assert_eq!(
            if members.weighted { "true" } else { "false" },
            cell(&rows, &case, "parameters.weighted"),
            "{case}"
        );
        if name == "unknown_key" {
            // The source warns "TraceFitter received the unknown parameter 'foo'!".
            assert_eq!(warnings.len(), 1, "{warnings:?}");
            assert!(
                warnings[0].contains("TraceFitter") && warnings[0].contains("'foo'"),
                "{warnings:?}"
            );
        } else {
            assert!(warnings.is_empty(), "{case}: {warnings:?}");
        }
    }
}

#[test]
fn to_param_round_trips_and_carries_the_restrictions() {
    for record in [
        TraceFitterParams::default(),
        TraceFitterParams {
            max_iteration: -7,
            weighted: true,
        },
        TraceFitterParams {
            max_iteration: i64::from(i32::MAX),
            weighted: false,
        },
    ] {
        let param = record.to_param().unwrap();
        assert_eq!(param.size(), 2);
        assert_eq!(
            *param.value("max_iteration").unwrap(),
            ParamValue::Integer(record.max_iteration)
        );
        assert_eq!(
            param.valid_strings("weighted").unwrap(),
            ["true".to_owned(), "false".to_owned()]
        );
        assert!(param.has_tag("weighted", "advanced").unwrap());
        assert!(param.has_tag("max_iteration", "advanced").unwrap());
        let (back, warnings) = TraceFitterParams::from_param(&param).unwrap();
        assert!(warnings.is_empty());
        assert_eq!(back, record);
    }
}

/// The current boundary of the `Param` round trip, a native difference
/// (`docs/TRACE_FITTER_SUPPORT.md`, native difference 8).
///
/// `to_param` writes every `i64`. `from_param` refuses values outside `i32`,
/// because `crate::param`'s restriction check converts integer entries to
/// `i32`. The executed product SDK
/// (`../oracle/gauss-trace-fitter/param-range/results/out.tsv`) accepts each
/// value below and stores it unchanged in `param_` and in `max_iterations_`.
/// When `crate::param` narrows as the source does, this test must change to
/// the source's outcome.
#[test]
fn to_param_and_from_param_disagree_beyond_i32() {
    for max_iteration in [i64::from(i32::MAX), i64::from(i32::MIN)] {
        let record = TraceFitterParams {
            max_iteration,
            weighted: false,
        };
        let (back, _) = TraceFitterParams::from_param(&record.to_param().unwrap()).unwrap();
        assert_eq!(back, record);
    }
    for max_iteration in [
        i64::from(i32::MAX) + 1,
        3_000_000_000,
        i64::from(i32::MIN) - 1,
        i64::MAX,
    ] {
        let record = TraceFitterParams {
            max_iteration,
            weighted: true,
        };
        let param = record.to_param().unwrap();
        assert_eq!(
            *param.value("max_iteration").unwrap(),
            ParamValue::Integer(max_iteration)
        );
        match TraceFitterParams::from_param(&param) {
            Err(Error::InvalidValue(message)) => {
                assert_eq!(message, "parameter value cannot be converted to i32");
            }
            other => panic!("{max_iteration}: {other:?}"),
        }
    }
}

// ---------------------------------------------------------------------------
// The driver
// ---------------------------------------------------------------------------

const TIMES: [f64; 6] = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5];
const DATA: [f64; 6] = [3.9, 3.1, 2.3, 1.9, 1.4, 1.2];

fn residual(x: &[f64], f: &mut [f64]) {
    for ((slot, t), y) in f.iter_mut().zip(TIMES).zip(DATA) {
        *slot = x[0] * libm::exp(-x[1] * t) + x[2] - y;
    }
}

fn jacobian(x: &[f64], jac: &mut DenseMatrix) {
    for (row, t) in TIMES.into_iter().enumerate() {
        let e = libm::exp(-x[1] * t);
        jac.set(row, 0, e);
        jac.set(row, 1, -x[0] * t * e);
        jac.set(row, 2, 1.0);
    }
}

const START: [f64; 3] = [1.0, 0.1, 0.0];

/// `optimize` is `minimize` with Eigen's defaults, `max_fev = max_iteration`
/// and an analytic Jacobian that consumes no budget, at every budget until
/// natural termination.
#[test]
fn optimize_is_the_configured_solver() {
    for max_iteration in 1..=60 {
        for weighted in [false, true] {
            let mut ours = START;
            let status = optimize_with_status(
                &mut ours,
                TIMES.len(),
                residual,
                jacobian,
                &TraceFitterParams {
                    max_iteration,
                    weighted,
                },
            )
            .unwrap();
            let mut direct = START;
            let expected = minimize(
                &mut direct,
                TIMES.len(),
                residual,
                |x: &[f64], jac: &mut DenseMatrix| {
                    jacobian(x, jac);
                    0
                },
                &LmParameters {
                    max_fev: usize::try_from(max_iteration).unwrap(),
                    ..LmParameters::default()
                },
            );
            assert_eq!(status, expected, "max_iteration {max_iteration}");
            assert_eq!(
                ours.map(f64::to_bits),
                direct.map(f64::to_bits),
                "max_iteration {max_iteration}"
            );

            let mut plain = START;
            optimize(
                &mut plain,
                TIMES.len(),
                residual,
                jacobian,
                &TraceFitterParams {
                    max_iteration,
                    weighted,
                },
            )
            .unwrap();
            assert_eq!(plain.map(f64::to_bits), ours.map(f64::to_bits));
        }
    }
    let defaults = LmParameters::default();
    assert_eq!(defaults.factor, 100.0);
    assert_eq!(defaults.ftol, f64::EPSILON.sqrt());
    assert_eq!(defaults.xtol, f64::EPSILON.sqrt());
    assert_eq!(defaults.gtol, 0.0);
}

/// An exhausted budget is accepted with the parameters reached, which
/// `optimize_is_the_configured_solver` shows are the solver's.
#[test]
fn an_exhausted_budget_is_accepted() {
    for max_iteration in [1, 2, 5] {
        let mut x = START;
        let status = optimize_with_status(
            &mut x,
            TIMES.len(),
            residual,
            jacobian,
            &TraceFitterParams {
                max_iteration,
                weighted: false,
            },
        )
        .unwrap();
        assert_eq!(
            status,
            LmStatus::TooManyFunctionEvaluation,
            "{max_iteration}"
        );
    }
}

fn counted(x: &mut [f64], values: usize, max_iteration: i64) -> (Result<LmStatus>, usize) {
    let calls = Cell::new(0usize);
    let result = optimize_with_status(
        x,
        values,
        |point: &[f64], f: &mut [f64]| {
            calls.set(calls.get() + 1);
            residual(point, f);
        },
        |point: &[f64], jac: &mut DenseMatrix| {
            calls.set(calls.get() + 1);
            jacobian(point, jac);
        },
        &TraceFitterParams {
            max_iteration,
            weighted: false,
        },
    );
    (result, calls.get())
}

#[test]
fn fewer_residuals_than_parameters_is_unable_to_fit() {
    let mut x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
    let (result, calls) = counted(&mut x, 6, 500);
    match result {
        Err(Error::InvalidValue(message)) => {
            assert_eq!(
                message,
                "UnableToFit-FinalSet: Skipping feature, we always expect N>=p"
            );
        }
        other => panic!("{other:?}"),
    }
    assert_eq!(calls, 0);
    assert_eq!(x, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]);
}

/// Eigen refuses `maxfev <= 0` and an empty parameter vector as
/// `ImproperInputParameters` before evaluating anything, and the source throws
/// "Error 0" (the executed product SDK agrees for 0 and -1; see
/// `tests/gauss_trace_fitter.rs`).
#[test]
fn improper_input_is_unable_to_fit_with_status_zero() {
    for max_iteration in [0, -1, i64::MIN] {
        let mut x = START;
        let (result, calls) = counted(&mut x, TIMES.len(), max_iteration);
        match result {
            Err(Error::InvalidValue(message)) => assert_eq!(
                message,
                "UnableToFit-FinalSet: Could not fit the gaussian to the data: Error 0"
            ),
            other => panic!("{max_iteration}: {other:?}"),
        }
        assert_eq!(calls, 0);
        assert_eq!(x, START);
    }
    let mut empty: [f64; 0] = [];
    let (result, calls) = counted(&mut empty, 6, 500);
    assert!(matches!(result, Err(Error::InvalidValue(ref m)) if m.ends_with("Error 0")));
    assert_eq!(calls, 0);
}

/// The solver's point ceiling is checked before anything is evaluated.
#[test]
fn oversized_problems_are_refused_before_evaluating() {
    let mut x = START;
    let (result, calls) = counted(&mut x, MAX_POINTS + 1, 500);
    match result {
        Err(Error::InvalidValue(message)) => {
            assert!(!message.starts_with(UNABLE_TO_FIT_FINAL_SET), "{message}")
        }
        other => panic!("{other:?}"),
    }
    assert_eq!(calls, 0);
    assert_eq!(x, START);
    assert_eq!(MAX_RESIDUAL_WORK, 1 << 30);
}

#[test]
fn a_huge_budget_behaves_like_any_budget_beyond_natural_termination() {
    let mut reference = START;
    let natural = optimize_with_status(
        &mut reference,
        TIMES.len(),
        residual,
        jacobian,
        &TraceFitterParams::default(),
    )
    .unwrap();
    let mut huge = START;
    let status = optimize_with_status(
        &mut huge,
        TIMES.len(),
        residual,
        jacobian,
        &TraceFitterParams {
            max_iteration: i64::MAX,
            weighted: false,
        },
    )
    .unwrap();
    assert_eq!(status, natural);
    assert_eq!(huge.map(f64::to_bits), reference.map(f64::to_bits));
}

#[test]
fn unable_to_fit_names_the_source_exception() {
    match unable_to_fit(FEWER_RESIDUALS_THAN_PARAMETERS) {
        Error::InvalidValue(message) => {
            assert_eq!(
                message,
                "UnableToFit-FinalSet: Skipping feature, we always expect N>=p"
            );
        }
        other => panic!("{other:?}"),
    }
}

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

#[test]
fn compute_theoretical_is_the_source_formula() {
    struct Constant(TraceFitterParams);
    impl TraceFitter for Constant {
        fn parameters(&self) -> &TraceFitterParams {
            &self.0
        }
        fn set_parameters(&mut self, parameters: TraceFitterParams) {
            self.0 = parameters;
        }
        fn fit(&mut self, _: &MassTraces) -> Result<()> {
            Ok(())
        }
        fn lower_rt_bound(&self) -> f64 {
            0.0
        }
        fn upper_rt_bound(&self) -> f64 {
            0.0
        }
        fn height(&self) -> f64 {
            0.0
        }
        fn center(&self) -> f64 {
            0.0
        }
        fn fwhm(&self) -> f64 {
            0.0
        }
        fn value(&self, rt: f64) -> f64 {
            rt * 3.0
        }
        fn area(&self) -> f64 {
            0.0
        }
        fn check_minimal_rt_span(&self, _: (f64, f64), _: f64) -> bool {
            false
        }
        fn check_maximal_rt_span(&self, _: f64) -> bool {
            false
        }
        fn gnuplot_formula(&self, _: &MassTrace, _: char, _: f64, _: f64) -> String {
            String::new()
        }
        fn compute_theoretical(&self, trace: &MassTrace, k: usize) -> Result<f64> {
            compute_theoretical(self, trace, k)
        }
    }
    let fitter = Constant(TraceFitterParams::default());
    let trace = MassTrace {
        theoretical_int: 0.25,
        peaks: vec![
            TracePeak::new(0, 0, 2.0, 0.0, 1.0),
            TracePeak::new(0, 1, 4.0, 0.0, 1.0),
        ],
        ..MassTrace::default()
    };
    assert_eq!(fitter.compute_theoretical(&trace, 0).unwrap(), 0.25 * 6.0);
    assert_eq!(fitter.compute_theoretical(&trace, 1).unwrap(), 0.25 * 12.0);
    assert!(matches!(
        fitter.compute_theoretical(&trace, 2),
        Err(Error::InvalidValue(_))
    ));
    let dynamic: &dyn TraceFitter = &fitter;
    assert_eq!(compute_theoretical(dynamic, &trace, 1).unwrap(), 3.0);
}

#[test]
fn stream_number_is_the_default_ostream() {
    assert_eq!(stream_number(0.0), "0");
    assert_eq!(stream_number(-0.0), "-0");
    assert_eq!(stream_number(680.1), "680.1");
    assert_eq!(stream_number(7.999665827), "7.99967");
    assert_eq!(stream_number(1234567.0), "1.23457e+06");
    assert_eq!(stream_number(0.000123456789), "0.000123457");
    assert_eq!(stream_number(f64::NAN), "nan");
}

fn single_trace(rts: &[f64], intensities: &[f32], baseline: f64) -> MassTraces {
    let mut trace = MassTrace {
        theoretical_int: 1.0,
        ..MassTrace::default()
    };
    for (k, (&rt, &intensity)) in rts.iter().zip(intensities).enumerate() {
        trace.peaks.push(TracePeak::new(0, k, rt, 500.0, intensity));
    }
    let mut traces = MassTraces::new();
    traces.push(trace);
    traces.baseline = baseline;
    traces
}

/// Worked by hand from the source loops (tier 4); the executed start values of
/// the Gaussian, which use the same helper, are in `tests/gauss_trace_fitter.rs`.
#[test]
fn initial_shape_follows_the_source_loops() {
    // Three points: the Gaussian takes the sums unsmoothed.
    let traces = single_trace(&[10.0, 11.0, 12.0], &[2.0, 5.0, 5.0], 0.0);
    let short = initial_shape(&traces, ProfileSmoothing::SkipShortProfiles).unwrap();
    assert_eq!(short.max_index, 1, "the first of equal maxima");
    assert_eq!(short.height, 5.0);
    assert_eq!(short.apex_rt, 11.0);
    assert_eq!(short.region_rt_span, 2.0);
    assert_eq!(
        (short.left_index, short.left_height, short.left_rt),
        (0, 2.0, 10.0)
    );
    assert_eq!(
        (short.right_index, short.right_height, short.right_rt),
        (2, 5.0, 12.0)
    );

    // EGH smooths every profile: totals [0, 0, 2, 5, 5, 0, 0], the running sum
    // starts at 0 + 2 + 5 and every smoothed value is 12 / 5.
    let always = initial_shape(&traces, ProfileSmoothing::Always).unwrap();
    assert_eq!(always.max_index, 0);
    assert_eq!(always.height, 12.0 / 5.0);
    assert_eq!(always.apex_rt, 10.0);
    assert_eq!((always.left_index, always.right_index), (0, 2));
    assert_eq!((always.right_height, always.right_rt), (12.0 / 5.0, 12.0));

    // Four points are smoothed by both; the baseline is subtracted.
    let traces = single_trace(&[10.0, 11.0, 12.0, 13.0], &[1.0, 3.0, 3.0, 1.0], 0.5);
    for smoothing in [
        ProfileSmoothing::SkipShortProfiles,
        ProfileSmoothing::Always,
    ] {
        let shape = initial_shape(&traces, smoothing).unwrap();
        assert_eq!(shape.max_index, 1);
        assert_eq!(shape.height, 8.0 / 5.0 - 0.5);
        assert_eq!((shape.left_index, shape.right_index), (0, 3));
    }

    assert!(matches!(
        initial_shape(&MassTraces::new(), ProfileSmoothing::Always),
        Err(Error::InvalidValue(_))
    ));
}
