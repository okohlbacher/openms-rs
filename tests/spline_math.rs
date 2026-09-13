// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Differential test of the OpenMS `MATH/MISC` spline family against an
//! executed probe of the pinned C++ sources.
//!
//! `tests/data/spline_math_cpp_probe.tsv` is the stdout of a driver that
//! compiles `CubicSpline2d.cpp`, `BSpline2d.cpp`, `BSplineSmoothingSpline.cpp`,
//! `SplineBisection.h` and the vendored `eol-bspline` library unmodified from
//! openms4-core bc9cc12, and prints every observation at `%.17g`. The driver
//! lives outside this repository under `../oracle/msc_splines/`; its sha256 and
//! build flags are recorded in `tests/data/spline_math_provenance.json`.
//!
//! Every comparison here is bit-for-bit. The probe is built with
//! `-ffp-contract=off` so no multiply-add is fused, and two strict builds at
//! `-O0` and `-O2` agree byte for byte. The same sources built with clang's
//! default contraction on arm64 differ in 1 234 of the 1 828 recorded values,
//! by a median of `8.9e-15` relative and by order unity on the dozen or so that
//! are numerically zero; that is a portability property of the C++ and is
//! measured in `docs/BSPLINE2D_SUPPORT.md` rather than papered over with a
//! tolerance.

use openms::processing::spline::{
    BSpline2d, BSplineSmoothingSpline, BoundaryCondition, CubicSpline2d, SplineFunction,
    spline_bisection,
};
use std::collections::BTreeMap;

const PROBE: &str = include_str!("data/spline_math_cpp_probe.tsv");
const SINUS: &str = include_str!("data/BSpline2d_test_sinus.txt");

/// Probe rows grouped by case, preserving the order they were printed in.
struct Probe(BTreeMap<String, Vec<(String, String)>>);

impl Probe {
    fn load() -> Self {
        let mut cases: BTreeMap<String, Vec<(String, String)>> = BTreeMap::new();
        for line in PROBE.lines() {
            let mut parts = line.split('\t');
            let (Some(case), Some(key), Some(value)) = (parts.next(), parts.next(), parts.next())
            else {
                continue;
            };
            cases
                .entry(case.to_string())
                .or_default()
                .push((key.to_string(), value.to_string()));
        }
        Self(cases)
    }
    fn rows(&self, case: &str) -> &[(String, String)] {
        self.0
            .get(case)
            .unwrap_or_else(|| panic!("probe has no case {case}"))
    }
    fn raw(&self, case: &str, key: &str) -> &str {
        self.rows(case)
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
            .unwrap_or_else(|| panic!("probe case {case} has no key {key}"))
    }
    fn number(&self, case: &str, key: &str) -> f64 {
        parse(self.raw(case, key))
    }
    fn integer(&self, case: &str, key: &str) -> i64 {
        self.raw(case, key).parse().unwrap()
    }
}

fn parse(text: &str) -> f64 {
    text.parse()
        .unwrap_or_else(|_| panic!("probe value {text} is not a number"))
}

fn sinus() -> (Vec<f64>, Vec<f64>) {
    let mut x = Vec::new();
    let mut y = Vec::new();
    for line in SINUS.lines() {
        let mut fields = line.split_whitespace();
        if let (Some(a), Some(b)) = (fields.next(), fields.next()) {
            x.push(parse(a));
            y.push(parse(b));
        }
    }
    (x, y)
}

fn upstream_peak() -> (Vec<f64>, Vec<f64>) {
    (
        vec![
            486.784, 486.787, 486.790, 486.793, 486.795, 486.797, 486.800, 486.802, 486.805,
            486.808, 486.811,
        ],
        vec![
            0.0, 154683.17, 620386.5, 1701390.12, 2848879.25, 3564045.5, 2744585.7, 1605583.0,
            1518984.0, 1591352.21, 1691345.1,
        ],
    )
}

fn squares(n: usize) -> (Vec<f64>, Vec<f64>) {
    let x: Vec<f64> = (0..n).map(|i| i as f64).collect();
    let y: Vec<f64> = x.iter().map(|v| v * v).collect();
    (x, y)
}

#[test]
fn probe_fixture_matches_the_sinus_data_file() {
    let probe = Probe::load();
    let (x, y) = sinus();
    assert_eq!(x.len() as i64, probe.integer("fixture", "rows"));
    assert_eq!(x.len(), y.len());
    assert_eq!(x.len(), 202);
}

// --------------------------------------------------------------- CubicSpline2d

fn check_cubic(probe: &Probe, case: &str, spline: &CubicSpline2d, derivatives: bool) {
    let mut checked = 0usize;
    for (key, value) in probe.rows(case) {
        let expected = parse(value);
        if let Some(position) = key.strip_prefix("eval@") {
            assert_eq!(
                spline.eval(parse(position)).unwrap(),
                expected,
                "{case} {key}"
            );
        } else if let Some(position) = key.strip_prefix("derivative_alias@") {
            assert_eq!(
                spline.derivative(parse(position), 1).unwrap(),
                expected,
                "{case} {key}"
            );
        } else if derivatives {
            for (prefix, order) in [("d1@", 1u8), ("d2@", 2), ("d3@", 3)] {
                if let Some(position) = key.strip_prefix(prefix) {
                    assert_eq!(
                        spline.derivative(parse(position), order).unwrap(),
                        expected,
                        "{case} {key}"
                    );
                }
            }
        }
        checked += 1;
    }
    assert!(checked > 0, "{case} had no comparable rows");
}

#[test]
fn cubic_spline_reproduces_every_probe_row() {
    let probe = Probe::load();
    let (mz, intensity) = upstream_peak();
    let vector = CubicSpline2d::new(&mz, &intensity).unwrap();
    check_cubic(&probe, "cubic_upstream", &vector, true);

    let pairs: Vec<(f64, f64)> = mz.iter().copied().zip(intensity.iter().copied()).collect();
    let mapped = CubicSpline2d::from_pairs(&pairs).unwrap();
    check_cubic(&probe, "cubic_upstream_map", &mapped, true);

    // The sine grid of the class test, built exactly as the test builds it.
    let (x_min, x_max) = (-0.5f64, 1.5f64);
    let n = 10usize;
    let x: Vec<f64> = (0..=n)
        .map(|i| x_min + (i as f64) / 10.0 * (x_max - x_min))
        .collect();
    let y: Vec<f64> = x.iter().map(|v| v.sin()).collect();
    let sine = CubicSpline2d::new(&x, &y).unwrap();
    check_cubic(&probe, "cubic_sine", &sine, true);
    assert_eq!(
        sine.derivative(x[0], 2).unwrap(),
        probe.number("cubic_sine", "d2_first")
    );
    assert_eq!(
        sine.derivative(x[n], 2).unwrap(),
        probe.number("cubic_sine", "d2_last")
    );
}

#[test]
fn cubic_spline_rejects_what_the_probe_records_as_throwing() {
    let probe = Probe::load();
    let (mz, intensity) = upstream_peak();
    let spline = CubicSpline2d::new(&mz, &intensity).unwrap();
    let thrown = |name: &str| probe.integer("cubic_errors", name) == 1;

    assert!(thrown("eval_below") && spline.eval(486.783).is_err());
    assert!(thrown("eval_above") && spline.eval(486.812).is_err());
    assert!(thrown("order_zero") && spline.derivative(486.79, 0).is_err());
    assert!(thrown("order_four") && spline.derivative(486.79, 4).is_err());
    assert!(thrown("unsorted") && CubicSpline2d::new(&[1.0, 0.0, 2.0], &[1.0, 2.0, 3.0]).is_err());
    assert!(thrown("too_short") && CubicSpline2d::new(&[1.0], &[1.0]).is_err());
    assert!(thrown("size_mismatch") && CubicSpline2d::new(&[1.0, 2.0], &[1.0]).is_err());

    // The one case where this port is stricter: the C++ check accepts a
    // repeated abscissa and then evaluates to NaN.
    assert_eq!(probe.integer("cubic_duplicate_x", "threw"), 0);
    assert_eq!(probe.integer("cubic_duplicate_x", "eval_is_nan"), 1);
    assert!(CubicSpline2d::new(&[0.0, 1.0, 1.0, 2.0], &[0.0, 1.0, 2.0, 3.0]).is_err());
}

// ------------------------------------------------------------------ BSpline2d

fn check_bspline(
    probe: &Probe,
    case: &str,
    x: &[f64],
    y: &[f64],
    wavelength: f64,
    boundary_condition: BoundaryCondition,
    num_nodes: usize,
) -> BSpline2d {
    let spline = BSpline2d::with_options(x, y, wavelength, boundary_condition, num_nodes)
        .unwrap_or_else(|error| panic!("{case} failed to fit: {error}"));
    assert_eq!(probe.integer(case, "ok"), 1, "{case} ok");
    assert_eq!(
        spline.node_count() as i64,
        probe.integer(case, "nNodes"),
        "{case} nNodes"
    );
    assert_eq!(x.len() as i64, probe.integer(case, "nX"), "{case} nX");
    assert_eq!(spline.domain().0, probe.number(case, "Xmin"), "{case} Xmin");
    assert_eq!(spline.domain().1, probe.number(case, "Xmax"), "{case} Xmax");
    assert_eq!(spline.alpha(), probe.number(case, "alpha"), "{case} alpha");
    for i in 0..spline.node_count() {
        assert_eq!(
            spline.coefficient(i),
            probe.number(case, &format!("coeff{i}")),
            "{case} coeff{i}"
        );
    }
    for (key, value) in probe.rows(case) {
        let expected = parse(value);
        if let Some(position) = key.strip_prefix("eval@") {
            assert_eq!(
                spline.eval(parse(position)).unwrap(),
                expected,
                "{case} {key}"
            );
        } else if let Some(position) = key.strip_prefix("slope@") {
            assert_eq!(
                spline.derivative(parse(position)).unwrap(),
                expected,
                "{case} {key}"
            );
        }
    }
    spline
}

#[test]
fn bspline_reproduces_the_upstream_sinus_fixture() {
    let probe = Probe::load();
    let (x, y) = sinus();
    check_bspline(
        &probe,
        "bspline_sinus_default",
        &x,
        &y,
        0.0,
        BoundaryCondition::ZeroSecond,
        0,
    );
    check_bspline(
        &probe,
        "bspline_sinus_wl2",
        &x,
        &y,
        2.0,
        BoundaryCondition::ZeroSecond,
        0,
    );
    check_bspline(
        &probe,
        "bspline_sinus_bc0_wl10",
        &x,
        &y,
        10.0,
        BoundaryCondition::ZeroEndpoints,
        0,
    );
    check_bspline(
        &probe,
        "bspline_sinus_bc1_wl1",
        &x,
        &y,
        1.0,
        BoundaryCondition::ZeroFirst,
        0,
    );
    // The class test also builds `BSpline2d(x, y, 100, BC_ZERO_SECOND)` and
    // discards it without calling ok(). It cannot succeed: 100 exceeds the
    // 11.57 span of the abscissae, so setDomain returns false and every
    // evaluation yields zero. The probe records that, and the port errors.
    assert_eq!(probe.integer("bspline_sinus_bc2_wl100", "ok"), 0);
    assert_eq!(probe.number("bspline_sinus_bc2_wl100", "eval@0"), 0.0);
    assert!(BSpline2d::with_options(&x, &y, 100.0, BoundaryCondition::ZeroSecond, 0).is_err());
}

#[test]
fn bspline_smoothing_beats_the_noise_by_the_margin_the_class_test_requires() {
    let probe = Probe::load();
    let (x, y) = sinus();

    let mut noisy = 0.0;
    for i in 0..x.len() {
        let error = y[i] - 10.0 * x[i].sin();
        noisy += error * error;
    }
    noisy /= x.len() as f64;
    assert_eq!(noisy, probe.number("bspline_mse", "noisy"));

    for (case, wavelength) in [("smoothed_default", 0.0), ("smoothed_wl2", 2.0)] {
        let spline =
            BSpline2d::with_options(&x, &y, wavelength, BoundaryCondition::ZeroSecond, 0).unwrap();
        let mut smoothed = 0.0;
        for &xi in &x {
            let error = spline.eval(xi).unwrap() - 10.0 * xi.sin();
            smoothed += error * error;
        }
        smoothed /= x.len() as f64;
        assert_eq!(smoothed, probe.number("bspline_mse", case));
        // The class test's assertion, reproduced.
        assert!(smoothed < 0.5 * noisy);
    }

    let spline = BSpline2d::with_options(&x, &y, 0.0, BoundaryCondition::ZeroSecond, 0).unwrap();
    let mut derivative_error = 0.0;
    for &xi in &x {
        derivative_error += (spline.derivative(xi).unwrap() - 10.0 * xi.cos()).abs();
    }
    derivative_error /= x.len() as f64;
    assert_eq!(
        derivative_error,
        probe.number("bspline_mse", "mean_abs_derivative_error")
    );
    assert!(derivative_error < 10.0 * 0.2);
}

#[test]
fn bspline_solve_refits_the_same_domain() {
    let probe = Probe::load();
    let (x, y) = sinus();
    let mut spline = BSpline2d::new(&x, &y).unwrap();
    let clean: Vec<f64> = x.iter().map(|v| 10.0 * v.sin()).collect();
    spline.solve(&clean).unwrap();
    assert_eq!(probe.integer("bspline_solve", "returned"), 1);
    for (key, value) in probe.rows("bspline_solve") {
        if let Some(position) = key.strip_prefix("eval@") {
            assert_eq!(
                spline.eval(parse(position)).unwrap(),
                parse(value),
                "bspline_solve {key}"
            );
        }
    }
}

#[test]
fn bspline_is_insensitive_to_the_order_of_the_abscissae() {
    let probe = Probe::load();
    let ascending: Vec<f64> = (0..10).map(|i| i as f64).collect();
    check_bspline(
        &probe,
        "bspline_ramp_asc",
        &ascending,
        &ascending,
        0.0,
        BoundaryCondition::ZeroSecond,
        100,
    );
    let descending: Vec<f64> = (1..=10).rev().map(|i| i as f64).collect();
    check_bspline(
        &probe,
        "bspline_ramp_desc",
        &descending,
        &descending,
        0.0,
        BoundaryCondition::ZeroSecond,
        100,
    );
}

#[test]
fn bspline_node_counts_wavelengths_and_boundary_conditions_match_the_probe() {
    let probe = Probe::load();
    let (x, y) = squares(8);
    check_bspline(
        &probe,
        "bspline_small_auto",
        &x,
        &y,
        0.0,
        BoundaryCondition::ZeroSecond,
        0,
    );
    check_bspline(
        &probe,
        "bspline_small_n4",
        &x,
        &y,
        0.0,
        BoundaryCondition::ZeroSecond,
        4,
    );
    check_bspline(
        &probe,
        "bspline_small_n6_bc0",
        &x,
        &y,
        0.0,
        BoundaryCondition::ZeroEndpoints,
        6,
    );
    check_bspline(
        &probe,
        "bspline_small_n6_bc1",
        &x,
        &y,
        0.0,
        BoundaryCondition::ZeroFirst,
        6,
    );

    let (wx, wy) = squares(21);
    check_bspline(
        &probe,
        "bspline_wave5",
        &wx,
        &wy,
        5.0,
        BoundaryCondition::ZeroSecond,
        0,
    );
    check_bspline(
        &probe,
        "bspline_wave12",
        &wx,
        &wy,
        12.0,
        BoundaryCondition::ZeroFirst,
        0,
    );
}

/// `num_nodes == 2` makes `BSplineBase::calculateQ`'s `Q.setup(M + 1, 3)` fail —
/// the dimension is below the bandwidth — and the C++ ignores the return value,
/// so the penalty matrix keeps the 1x1 shape of its default constructor, the
/// second coefficient is never solved for, and the object still reports `ok()`.
/// The port reproduces that rather than refusing it; these are the probe's own
/// numbers for the release (`NDEBUG`) build the oracle records. With assertions
/// enabled the same C++ aborts in `Beta`, because this path reads outside
/// `BoundaryConditions`; every one of those reads is consumed by a write the
/// banded storage discards, which is why the values below are reproducible at
/// all.
#[test]
fn bspline_two_node_grids_reproduce_the_degenerate_source_fit() {
    let probe = Probe::load();
    let (x, y) = squares(8);
    let two = check_bspline(
        &probe,
        "bspline_small_n2",
        &x,
        &y,
        0.0,
        BoundaryCondition::ZeroSecond,
        2,
    );
    assert!(two.ok());
    assert_eq!(two.node_count(), 2);
    assert_eq!(two.node_spacing(), 7.0);
    check_bspline(
        &probe,
        "bspline_small_n2_bc0",
        &x,
        &y,
        0.0,
        BoundaryCondition::ZeroEndpoints,
        2,
    );
    check_bspline(
        &probe,
        "bspline_small_n2_bc1",
        &x,
        &y,
        0.0,
        BoundaryCondition::ZeroFirst,
        2,
    );
    // Three nodes is the smallest grid whose matrix the C++ does shape, and it
    // is a fit rather than an artefact: at 3.5 the two-node curve reports 23.7
    // where the data is 12.25, the three-node one 11.68.
    let three = check_bspline(
        &probe,
        "bspline_small_n3",
        &x,
        &y,
        0.0,
        BoundaryCondition::ZeroSecond,
        3,
    );
    assert!(three.eval(3.5).unwrap() < two.eval(3.5).unwrap());

    // The unsolved second coefficient is the raw right-hand side, so negating
    // the ordinates negates it exactly like the solved one. Derived from the
    // linearity of the accumulation, not transcribed.
    let mut flipped = two;
    let negated: Vec<f64> = y.iter().map(|v| -v).collect();
    flipped.solve(&negated).unwrap();
    assert_eq!(flipped.coefficient(0), 8.011010956384775);
    assert_eq!(flipped.coefficient(1), -15.00000000000001);
}

#[test]
fn bspline_setup_failures_become_errors() {
    let probe = Probe::load();
    assert_eq!(probe.integer("bspline_wl_too_long", "ok"), 0);
    assert_eq!(probe.number("bspline_wl_too_long", "eval@1.5"), 0.0);
    let x = [0.0, 1.0, 2.0, 3.0];
    let y = [0.0, 1.0, 4.0, 9.0];
    assert!(BSpline2d::with_options(&x, &y, 100.0, BoundaryCondition::ZeroSecond, 0).is_err());

    assert_eq!(probe.integer("bspline_small_wl3", "ok"), 0);
    let (sx, sy) = squares(8);
    assert!(BSpline2d::with_options(&sx, &sy, 3.0, BoundaryCondition::ZeroSecond, 0).is_err());
}

// ------------------------------------------------- BSplineSmoothingSpline

fn check_smoothing(probe: &Probe, case: &str, x: &[f64], y: &[f64], s: f64, degree: i32) {
    let expected_ok = probe.integer(case, "ok") == 1;
    let fitted = BSplineSmoothingSpline::with_smoothing(x, y, s, degree);
    assert_eq!(fitted.is_ok(), expected_ok, "{case} ok");
    let Ok(spline) = fitted else { return };
    assert_eq!(
        i64::from(spline.num_interior_knots()),
        probe.integer(case, "num_interior_knots"),
        "{case} num_interior_knots"
    );
    assert_eq!(spline.rss(), probe.number(case, "rss"), "{case} rss");
    assert_eq!(
        spline.smoothing_param(),
        probe.number(case, "smoothing_param"),
        "{case} smoothing_param"
    );
    for (key, value) in probe.rows(case) {
        if let Some(position) = key.strip_prefix("eval@") {
            assert_eq!(
                spline.eval(parse(position)).unwrap(),
                parse(value),
                "{case} {key}"
            );
        }
    }
}

#[test]
fn smoothing_spline_reproduces_every_probe_case() {
    let probe = Probe::load();
    let ramp5: Vec<f64> = (0..5).map(|i| i as f64).collect();
    let squares5 = vec![0.0, 1.0, 4.0, 9.0, 16.0];
    let ramp6: Vec<f64> = (0..6).map(|i| i as f64).collect();
    let shifted5: Vec<f64> = (1..6).map(|i| i as f64).collect();
    let shifted6: Vec<f64> = (1..7).map(|i| i as f64).collect();

    check_smoothing(&probe, "smooth_linear5_auto", &ramp5, &ramp5, -1.0, 3);
    check_smoothing(&probe, "smooth_linear5_zero", &ramp5, &ramp5, 0.0, 3);
    check_smoothing(&probe, "smooth_linear5_s2", &ramp5, &ramp5, 2.0, 3);
    check_smoothing(&probe, "smooth_square5_auto", &ramp5, &squares5, -1.0, 3);
    check_smoothing(&probe, "smooth_ramp6_zero", &ramp6, &shifted6, 0.0, 3);
    check_smoothing(&probe, "smooth_ramp6_s10", &ramp6, &shifted6, 10.0, 3);
    check_smoothing(&probe, "smooth_ramp6_auto", &ramp6, &shifted6, -1.0, 3);
    check_smoothing(&probe, "smooth_ramp6_s5", &ramp6, &shifted6, 5.0, 3);
    check_smoothing(&probe, "smooth_ramp5_auto", &ramp5, &shifted5, -1.0, 3);
    check_smoothing(&probe, "smooth_ramp5_s2", &ramp5, &shifted5, 2.0, 3);

    let noisy_x: Vec<f64> = (0..20).map(|i| f64::from(i) * 0.3).collect();
    let noisy_y: Vec<f64> = (0..20)
        .map(|i| (f64::from(i) * 0.3).sin() + 0.2 * f64::from(i % 3 - 1))
        .collect();
    check_smoothing(&probe, "smooth_noisy20_s20", &noisy_x, &noisy_y, 20.0, 3);
    check_smoothing(&probe, "smooth_noisy20_auto", &noisy_x, &noisy_y, -1.0, 3);

    let wiggle = vec![1.0, 3.0, 2.0, 4.0, 3.5];
    check_smoothing(&probe, "smooth_wiggle5_zero", &ramp5, &wiggle, 0.0, 3);
    check_smoothing(&probe, "smooth_wiggle5_s2", &ramp5, &wiggle, 2.0, 3);

    let line: Vec<f64> = ramp6.iter().map(|v| 2.0 * v + 1.0).collect();
    check_smoothing(&probe, "smooth_line6_auto", &ramp6, &line, -1.0, 3);

    let quad_x: Vec<f64> = (-2..4).map(f64::from).collect();
    let quad_y: Vec<f64> = quad_x.iter().map(|v| v * v).collect();
    check_smoothing(&probe, "smooth_quad6_s1", &quad_x, &quad_y, 1.0, 3);

    let degrees = vec![1.0, 2.0, 1.5, 3.0, 2.5, 4.0];
    check_smoothing(&probe, "smooth_degrees_k1", &ramp6, &degrees, -1.0, 1);
    check_smoothing(&probe, "smooth_degrees_k2", &ramp6, &degrees, -1.0, 2);
    check_smoothing(&probe, "smooth_degrees_k3", &ramp6, &degrees, -1.0, 3);

    let alt_x: Vec<f64> = (0..12).map(f64::from).collect();
    let alt_y: Vec<f64> = (0..12).map(|i| f64::from(i % 2)).collect();
    check_smoothing(&probe, "smooth_alternating12_s05", &alt_x, &alt_y, 0.5, 3);
    check_smoothing(&probe, "smooth_alternating12_s3", &alt_x, &alt_y, 3.0, 3);
    check_smoothing(&probe, "smooth_alternating12_auto", &alt_x, &alt_y, -1.0, 3);

    check_smoothing(&probe, "smooth_n2", &[0.0, 1.0], &[0.0, 1.0], -1.0, 3);
    check_smoothing(
        &probe,
        "smooth_n3",
        &[0.0, 1.0, 2.0],
        &[1.0, 2.0, 3.0],
        -1.0,
        3,
    );

    check_smoothing(
        &probe,
        "smooth_duplicate_x",
        &[1.0, 1.0, 2.0, 3.0],
        &[1.0, 2.0, 3.0, 4.0],
        -1.0,
        3,
    );
    check_smoothing(
        &probe,
        "smooth_unsorted",
        &[0.0, 2.0, 1.0, 3.0],
        &[1.0, 2.0, 3.0, 4.0],
        -1.0,
        3,
    );
    check_smoothing(&probe, "smooth_size_mismatch", &ramp5, &[1.0, 2.0], -1.0, 3);
    check_smoothing(&probe, "smooth_single_point", &[0.0], &[1.0], -1.0, 3);
}

/// The candidate `std::sort`'s comparator is not a strict weak ordering, so
/// which candidate survives at position zero is not settled by the comparator
/// alone. These four cases pin it against the executed C++ instead of against an
/// argument about what a standard library does with a short range — the oracle
/// was built with Apple clang and therefore libc++, whose `std::sort` sends
/// three-, four- and five-element ranges through the `__sort3`/`__sort4`/
/// `__sort5` networks rather than through an insertion sort, and libstdc++ would
/// not be obliged to agree.
///
/// * `smooth_tie_break` — eight samples of `exp(-x^2)`, budget `6e-4`. The
///   six- and eight-node grids are both inside the budget and 3.9e-4 apart, so
///   the "prefer fewer knots" branch takes the six-node fit even though the
///   eight-node fit is nearer the budget.
/// * `smooth_tie_break_off` — the same data at `5.5e-4`, which puts the
///   eight-node fit over budget so the branch cannot fire and it wins instead.
/// * `smooth_cycle` — nine samples of `6*sin(x)` on abscissae around 1000, so
///   the polynomial branch's normal equations are useless and the node search
///   always runs. At `4.24e-3` the relation over the three leading candidates is
///   a cycle: six beats eight and eight beats nine on knot count, while nine
///   beats six on closeness.
/// * `smooth_cycle_off` — the same data at `4.20e-3`, where the nine-node fit
///   is over budget and the cycle disappears.
#[test]
fn smoothing_spline_candidate_selection_is_pinned_where_the_comparator_is_not() {
    let probe = Probe::load();

    let gauss_x: Vec<f64> = (0..8).map(|i| f64::from(i) * 0.25).collect();
    let gauss_y: Vec<f64> = gauss_x.iter().map(|v| (-v * v).exp()).collect();
    check_smoothing(&probe, "smooth_tie_break", &gauss_x, &gauss_y, 0.0006, 3);
    check_smoothing(
        &probe,
        "smooth_tie_break_off",
        &gauss_x,
        &gauss_y,
        0.00055,
        3,
    );
    // The two budgets differ by 5e-5 and select different fits, which is the
    // whole point: the knot-count branch, not closeness, decided the first.
    let tight = BSplineSmoothingSpline::with_smoothing(&gauss_x, &gauss_y, 0.0006, 3).unwrap();
    let loose = BSplineSmoothingSpline::with_smoothing(&gauss_x, &gauss_y, 0.00055, 3).unwrap();
    assert_eq!(tight.num_interior_knots(), 4);
    assert_eq!(loose.num_interior_knots(), 6);
    assert!(tight.rss() < loose.rss());

    let shift_x: Vec<f64> = (0..9).map(|i| 1000.0 + f64::from(i) * 0.25).collect();
    let shift_y: Vec<f64> = (0..9).map(|i| 6.0 * (f64::from(i) * 0.25).sin()).collect();
    check_smoothing(&probe, "smooth_cycle", &shift_x, &shift_y, 0.00424, 3);
    check_smoothing(&probe, "smooth_cycle_off", &shift_x, &shift_y, 0.0042, 3);
    let cyclic = BSplineSmoothingSpline::with_smoothing(&shift_x, &shift_y, 0.00424, 3).unwrap();
    assert_eq!(cyclic.num_interior_knots(), 4);
    let acyclic = BSplineSmoothingSpline::with_smoothing(&shift_x, &shift_y, 0.0042, 3).unwrap();
    assert_eq!(acyclic.num_interior_knots(), 7);
}

// ------------------------------------------------------------- SplineBisection

struct Parabola {
    peak: f64,
    height: f64,
    a: f64,
}

impl SplineFunction for Parabola {
    fn eval(&self, x: f64) -> openms::Result<f64> {
        Ok(-self.a * (x - self.peak) * (x - self.peak) + self.height)
    }
    fn first_derivative(&self, x: f64) -> openms::Result<f64> {
        Ok(-2.0 * self.a * (x - self.peak))
    }
}

#[test]
fn spline_bisection_reproduces_every_probe_bracket() {
    let probe = Probe::load();
    let cases: [(&str, f64, f64, f64, f64, f64, f64); 8] = [
        ("centered", 500.0, 1000.0, 1.0, 499.0, 501.0, 1e-6),
        ("offcenter", 500.3, 750.0, 2.0, 499.0, 501.0, 1e-6),
        ("tight", 500.123456, 100.0, 5.0, 499.0, 501.0, 1e-9),
        ("outside_right", 505.0, 50.0, 1.0, 499.0, 501.0, 1e-6),
        ("outside_left", 495.0, 50.0, 1.0, 499.0, 501.0, 1e-6),
        ("apex_on_left_edge", 499.0, 10.0, 1.0, 499.0, 501.0, 1e-6),
        ("reversed_bracket", 500.0, 1000.0, 1.0, 501.0, 499.0, 1e-6),
        ("coarse_threshold", 500.3, 750.0, 2.0, 499.0, 501.0, 0.1),
    ];
    for (name, peak, height, a, left, right, threshold) in cases {
        let parabola = Parabola { peak, height, a };
        let (position, value) = spline_bisection(&parabola, left, right, threshold).unwrap();
        assert_eq!(
            position,
            probe.number("bisection", &format!("{name}_mz")),
            "{name} position"
        );
        assert_eq!(
            value,
            probe.number("bisection", &format!("{name}_int")),
            "{name} value"
        );
    }
}

#[test]
fn spline_bisection_over_a_cubic_spline_finds_the_probe_apex() {
    let probe = Probe::load();
    let (mz, intensity) = upstream_peak();
    let spline = CubicSpline2d::new(&mz, &intensity).unwrap();
    let (position, value) = spline_bisection(&spline, 486.795, 486.800, 1e-6).unwrap();
    assert_eq!(position, probe.number("bisection_cubic", "mz"));
    assert_eq!(value, probe.number("bisection_cubic", "int"));
    let (position, value) = spline_bisection(&spline, 486.795, 486.800, 1e-9).unwrap();
    assert_eq!(position, probe.number("bisection_cubic", "mz_tight"));
    assert_eq!(value, probe.number("bisection_cubic", "int_tight"));

    // A bracket outside the knot range propagates the spline's own error rather
    // than evaluating undefined coefficients.
    assert!(spline_bisection(&spline, 486.0, 487.0, 1e-6).is_err());
}

#[test]
fn spline_bisection_also_drives_a_bspline() {
    let (x, y) = squares(21);
    let spline = BSpline2d::with_options(&x, &y, 5.0, BoundaryCondition::ZeroSecond, 0).unwrap();
    // A B-spline never fails on range, so the bisection can run outside the
    // fitted domain; the derivative there is zero, which stops it immediately.
    let (position, _) = spline_bisection(&spline, -60.0, -40.0, 1e-6).unwrap();
    assert_eq!(position, -50.0);
}
