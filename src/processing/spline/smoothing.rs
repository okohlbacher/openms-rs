// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Smoothing spline with a residual-sum-of-squares budget.
//!
//! Ports `src/openms/include/OpenMS/MATH/MISC/BSplineSmoothingSpline.h` and
//! `src/openms/source/MATH/MISC/BSplineSmoothingSpline.cpp`. See
//! `docs/BSPLINE_SMOOTHING_SPLINE_SUPPORT.md`.
//!
//! The class exists to approximate `scipy.interpolate.UnivariateSpline` for
//! PyProphet compatibility. It does not solve the penalised problem
//! `min RSS(f) + lambda * integral f''^2` that its own documentation states;
//! what it actually does is try a cubic *polynomial* first and, only if that
//! misses the residual budget, search a short fixed list of node counts for the
//! [`BSpline2d`](crate::processing::spline::BSpline2d) fit whose residual sum of
//! squares comes closest to the budget. That search is ported as written,
//! because its output is what callers see.

use crate::processing::spline::b_spline::{BSpline2d, BoundaryCondition};
use crate::{Error, Result};

/// Degree of the polynomial the first fitting attempt uses.
const POLYNOMIAL_DEGREE: usize = 3;

#[derive(Clone, Debug)]
enum Fit {
    /// Coefficients `[a0, a1, a2, a3]` of `a0 + a1 x + a2 x^2 + a3 x^3`.
    Polynomial(Vec<f64>),
    Spline(Box<BSpline2d>),
}

/// Smoothing spline balancing closeness to the data against model size.
///
/// The smoothing parameter `s` is a budget on the residual sum of squares, not
/// a curvature weight:
///
/// * `s = 0` asks for interpolation — in practice the densest
///   [`BSpline2d`] grid, which passes
///   close to but not exactly through the data;
/// * `s > 0` allows the curve to deviate from the data;
/// * a negative `s` requests the scipy default `m - sqrt(2 m)` for `m` points;
/// * larger `s` approaches a plain cubic polynomial fit.
///
/// For fewer than four points the constructor falls back to the interpolating
/// branch, as the source's `@note` records.
///
/// # Differences from the source
///
/// * Construction returns a `Result` rather than leaving a value whose `ok()` is
///   false; every failure path that sets `ok_ = false` in C++ is an `Err` here,
///   and [`BSplineSmoothingSpline::eval`] therefore never returns the source's
///   `NaN`.
/// * The source sorts its candidate fits with `std::sort` under a comparator
///   that is not a strict weak ordering, so which candidate ends up first is
///   left to the standard library and genuinely differs between
///   implementations. The comparator is reproduced verbatim and the sort is a
///   stable insertion sort; the agreement with the C++ is established by the
///   executed probe, which was built against libc++, rather than argued from
///   the size of the candidate list. `docs/BSPLINE_SMOOTHING_SPLINE_SUPPORT.md`
///   has the detail and `tests/spline_math.rs` the cases that pin it.
#[derive(Clone, Debug)]
pub struct BSplineSmoothingSpline {
    fit: Fit,
    num_interior_knots: i32,
    rss: f64,
    s: f64,
    degree: i32,
}

impl BSplineSmoothingSpline {
    /// Maximum number of data points, inherited from
    /// [`BSpline2d::MAX_POINTS`](crate::processing::spline::BSpline2d::MAX_POINTS).
    ///
    /// Native addition; the source bounds nothing.
    pub const MAX_POINTS: usize = BSpline2d::MAX_POINTS;

    /// Fit a smoothing spline with the scipy default smoothing parameter and a
    /// cubic degree, equivalent to the C++ default arguments `s = -1.0`,
    /// `k = 3`.
    ///
    /// # Errors
    ///
    /// As [`BSplineSmoothingSpline::with_smoothing`].
    pub fn new(x: &[f64], y: &[f64]) -> Result<Self> {
        Self::with_smoothing(x, y, -1.0, 3)
    }

    /// Fit a smoothing spline with automatic knot selection.
    ///
    /// # Arguments
    ///
    /// * `x` — the data positions; must be strictly increasing. The source's
    ///   documentation says "sorted", but its check rejects equal neighbours
    ///   too, and its class test relies on that.
    /// * `y` — the values at those positions, one per `x`.
    /// * `s` — the residual-sum-of-squares budget. Negative selects the scipy
    ///   default `m - sqrt(2 m)`.
    /// * `degree` — the source's `k`. Retained for signature parity and
    ///   reported by [`BSplineSmoothingSpline::degree`], but it selects nothing:
    ///   the C++ member is declared `[[maybe_unused]]` and every fit is cubic.
    ///   Its class test checks `k` in 1..=3 and gets identical curves, which
    ///   this port reproduces.
    ///
    /// # How the fit is chosen
    ///
    /// With fewer than four points, or once `s` resolves to zero or less, the
    /// constructor fits a single [`BSpline2d`]
    /// with an automatic node count and reports `n - 2` interior knots without
    /// checking how well it fits.
    ///
    /// Otherwise it first fits a cubic polynomial by normal equations with
    /// partial pivoting and accepts it when its residual sum of squares is
    /// within ten percent of the budget — the source's margin — reporting zero
    /// interior knots. Failing that it builds one B-spline per entry of the node
    /// list `4, 6, 8, max(4, n/2), max(4, 3n/4), n` — only the two middle terms
    /// are clamped, the divisions truncate, and the list is then sorted and
    /// deduplicated — keeps those that fit, and picks the one whose residual sum
    /// of squares is closest to the budget, preferring fewer interior knots when
    /// two are within `0.001` of each other and both under budget. That last
    /// rule is not a tiebreak on an otherwise total order; see the note on the
    /// type.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when fewer than two points are supplied,
    /// when `x` and `y` differ in length, when `x` is not strictly increasing,
    /// when a value or `s` is not finite, when more than
    /// [`BSplineSmoothingSpline::MAX_POINTS`] points are supplied, or when no
    /// candidate fit can be built at all — the last being the source's "No valid
    /// spline configurations found".
    pub fn with_smoothing(x: &[f64], y: &[f64], s: f64, degree: i32) -> Result<Self> {
        let n = x.len();
        if n < 2 || y.len() != n {
            return Err(Error::InvalidValue(
                "smoothing spline needs at least two matching points".into(),
            ));
        }
        if n > Self::MAX_POINTS {
            return Err(Error::InvalidValue(
                "smoothing spline point count exceeds the maximum".into(),
            ));
        }
        if x.iter().chain(y).any(|v| !v.is_finite()) || !s.is_finite() {
            return Err(Error::InvalidValue(
                "smoothing spline data and smoothing parameter must be finite".into(),
            ));
        }
        for i in 1..n {
            if x[i] <= x[i - 1] {
                return Err(Error::UnsortedData);
            }
        }

        let s_ = if s < 0.0 {
            let m = n as f64;
            m - (2.0 * m).sqrt()
        } else {
            s
        };

        if n < 4 || s_ <= 0.0 {
            let spline = BSpline2d::with_options(x, y, 0.0, BoundaryCondition::ZeroSecond, 0)?;
            let rss = spline_rss(&spline, x, y)?;
            return Ok(Self {
                fit: Fit::Spline(Box::new(spline)),
                num_interior_knots: (n as i32) - 2,
                rss,
                s: s_,
                degree,
            });
        }

        Self::fit_smoothing_spline(x, y, s_, degree)
    }

    fn fit_smoothing_spline(x: &[f64], y: &[f64], s_target: f64, degree: i32) -> Result<Self> {
        if let Some(coefficients) = try_polynomial_fit(x, y, s_target, POLYNOMIAL_DEGREE) {
            let rss = polynomial_rss(&coefficients, x, y);
            return Ok(Self {
                fit: Fit::Polynomial(coefficients),
                num_interior_knots: 0,
                rss,
                s: s_target,
                degree,
            });
        }

        let n = x.len() as i64;
        let mut node_counts = vec![4i64, 6, 8, 4.max(n / 2), 4.max(3 * n / 4), n];
        node_counts.sort_unstable();
        node_counts.dedup();

        struct Candidate {
            num_interior_knots: i32,
            rss: f64,
            spline: BSpline2d,
        }
        let mut candidates: Vec<Candidate> = Vec::new();
        for num_nodes in node_counts {
            let Ok(count) = usize::try_from(num_nodes) else {
                continue;
            };
            let Ok(spline) =
                BSpline2d::with_options(x, y, 0.0, BoundaryCondition::ZeroSecond, count)
            else {
                continue;
            };
            let Ok(rss) = spline_rss(&spline, x, y) else {
                continue;
            };
            candidates.push(Candidate {
                num_interior_knots: (num_nodes - 2) as i32,
                rss,
                spline,
            });
        }
        if candidates.is_empty() {
            return Err(Error::InvalidValue(
                "smoothing spline found no valid spline configuration".into(),
            ));
        }

        let keys: Vec<(f64, i32)> = candidates
            .iter()
            .map(|c| (c.rss, c.num_interior_knots))
            .collect();
        let best = candidates.swap_remove(best_candidate(&keys, s_target));
        Ok(Self {
            fit: Fit::Spline(Box::new(best.spline)),
            num_interior_knots: best.num_interior_knots,
            rss: best.rss,
            s: s_target,
            degree,
        })
    }

    /// Evaluate the smoothing spline at `x`.
    ///
    /// Defined everywhere: a polynomial fit extrapolates as a cubic, and a
    /// B-spline fit decays to the mean of the ordinates more than two node
    /// intervals outside the data, as
    /// [`BSpline2d::eval`](crate::processing::spline::BSpline2d::eval)
    /// describes. The source returns `NaN` when the fit failed; here a failed
    /// fit is an `Err` from the constructor instead, so no value of this type
    /// can be in that state.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `x` is not finite or the evaluation
    /// overflows.
    pub fn eval(&self, x: f64) -> Result<f64> {
        match &self.fit {
            Fit::Polynomial(coefficients) => {
                if !x.is_finite() {
                    return Err(Error::InvalidValue(
                        "smoothing spline evaluation position must be finite".into(),
                    ));
                }
                let value = eval_polynomial(coefficients, x);
                if value.is_finite() {
                    Ok(value)
                } else {
                    Err(Error::InvalidValue(
                        "smoothing spline evaluation is not finite".into(),
                    ))
                }
            }
            Fit::Spline(spline) => spline.eval(x),
        }
    }

    /// Number of interior knots the fit selected.
    ///
    /// Zero for the polynomial branch, `n - 2` for the interpolating branch and
    /// `num_nodes - 2` for a chosen B-spline candidate — the source's own
    /// bookkeeping, which counts nodes of a uniform grid rather than knots of a
    /// B-spline basis.
    pub fn num_interior_knots(&self) -> i32 {
        self.num_interior_knots
    }

    /// Residual sum of squares the chosen fit achieves on the input data.
    pub fn rss(&self) -> f64 {
        self.rss
    }

    /// Smoothing parameter in force, after a negative request has been replaced
    /// by the scipy default `m - sqrt(2 m)`.
    pub fn smoothing_param(&self) -> f64 {
        self.s
    }

    /// The `degree` the fit was asked for.
    ///
    /// Native accessor for the source's `[[maybe_unused]] int k_`. It reports
    /// what the caller passed; it does not describe the curve, which is always
    /// cubic.
    pub fn degree(&self) -> i32 {
        self.degree
    }

    /// Whether the fit is the polynomial branch rather than a B-spline.
    ///
    /// Native accessor for the source's private `fit_type_`.
    pub fn is_polynomial(&self) -> bool {
        matches!(self.fit, Fit::Polynomial(_))
    }
}

/// The source's candidate comparator, transcribed. `a` and `b` are
/// `(rss, num_interior_knots)`.
///
/// It is not a strict weak ordering: when two candidates are both inside the
/// budget and their distances to it differ by less than `0.001`, the knot count
/// decides instead of the distance, and that rule can contradict the distance
/// rule applied to a third candidate. See [`best_candidate`].
fn candidate_precedes(a: (f64, i32), b: (f64, i32), s_target: f64) -> bool {
    let err_a = (a.0 - s_target).abs();
    let err_b = (b.0 - s_target).abs();
    if a.0 <= s_target && b.0 <= s_target && (err_a - err_b).abs() < 0.001 {
        return a.1 < b.1;
    }
    err_a < err_b
}

/// Index of the candidate the source's `std::sort` leaves at position zero,
/// which is the only position it then reads.
///
/// The source sorts the candidates with `std::sort` under
/// [`candidate_precedes`], which is not a strict weak ordering, so the standard
/// library is free to produce any permutation and different implementations do.
/// The candidate list holds one to six entries; in libc++'s `__algorithm/sort.h`
/// that range is dispatched by a `switch (__len)` to the `__sort3`, `__sort4`
/// and `__sort5` comparison networks, and only a length of six or more reaches
/// `__insertion_sort` — and libstdc++'s insertion sort is a different algorithm
/// again, so it is not obliged to agree with either. Which library ran is
/// therefore part of the answer, and no argument about "small ranges use
/// insertion sort" settles it.
///
/// This port does a stable insertion sort — element `j` walks left while it
/// precedes its neighbour — and the agreement with the C++ is established by
/// measurement, not by that argument: the executed probe was built against
/// libc++ and `tests/spline_math.rs` reproduces its selection for every probed
/// dataset, including `smooth_tie_break`, where the knot-count branch overrides
/// the distance rule, and `smooth_cycle`, where the comparator is cyclic over
/// the three leading candidates so the winner depends on the sort algorithm
/// itself. Against another standard library the C++ may well select a different
/// candidate there; this port's choice is the one pinned by those tests.
fn best_candidate(keys: &[(f64, i32)], s_target: f64) -> usize {
    let mut order: Vec<usize> = (0..keys.len()).collect();
    for i in 1..order.len() {
        let mut j = i;
        while j > 0 && candidate_precedes(keys[order[j]], keys[order[j - 1]], s_target) {
            order.swap(j, j - 1);
            j -= 1;
        }
    }
    order.first().copied().unwrap_or(0)
}

/// `BSplineSmoothingSpline::eval_polynomial`, in its original accumulation
/// order: ascending powers, with the power carried in a running product.
fn eval_polynomial(coefficients: &[f64], x: f64) -> f64 {
    let mut result = 0.0;
    let mut xpow = 1.0;
    for &c in coefficients {
        result += c * xpow;
        xpow *= x;
    }
    result
}

/// `BSplineSmoothingSpline::compute_polynomial_rss`.
fn polynomial_rss(coefficients: &[f64], x: &[f64], y: &[f64]) -> f64 {
    let mut rss = 0.0;
    for i in 0..x.len() {
        let fitted = eval_polynomial(coefficients, x[i]);
        let residual = y[i] - fitted;
        rss += residual * residual;
    }
    rss
}

/// `BSplineSmoothingSpline::compute_rss`.
fn spline_rss(spline: &BSpline2d, x: &[f64], y: &[f64]) -> Result<f64> {
    let mut rss = 0.0;
    for i in 0..x.len() {
        let fitted = spline.eval(x[i])?;
        let residual = y[i] - fitted;
        rss += residual * residual;
    }
    Ok(rss)
}

/// `BSplineSmoothingSpline::try_polynomial_fit`: least squares through the
/// normal equations `X^T X c = X^T y`, solved by Gaussian elimination with
/// partial pivoting.
///
/// Returns the coefficients only when the residual sum of squares is within the
/// source's ten-percent margin of `s_target`; `None` otherwise, including when
/// there are fewer points than coefficients or the pivot falls below `1e-10`.
fn try_polynomial_fit(x: &[f64], y: &[f64], s_target: f64, degree: usize) -> Option<Vec<f64>> {
    let n = x.len();
    let m = degree + 1;
    if n < m {
        return None;
    }

    let mut xtx = vec![0.0; m * m];
    let mut xty = vec![0.0; m];
    let mut powers = vec![0.0; m];
    for i in 0..n {
        let xi = x[i];
        let yi = y[i];
        powers[0] = 1.0;
        for p in 1..m {
            powers[p] = powers[p - 1] * xi;
        }
        for r in 0..m {
            xty[r] += powers[r] * yi;
            for c in 0..m {
                xtx[r * m + c] += powers[r] * powers[c];
            }
        }
    }

    for i in 0..m {
        let mut pivot = i;
        let mut maxval = xtx[i * m + i].abs();
        for j in i + 1..m {
            let value = xtx[j * m + i].abs();
            if value > maxval {
                maxval = value;
                pivot = j;
            }
        }
        if maxval < 1e-10 {
            return None;
        }
        if pivot != i {
            for c in 0..m {
                xtx.swap(i * m + c, pivot * m + c);
            }
            xty.swap(i, pivot);
        }
        let diag = xtx[i * m + i];
        for j in i + 1..m {
            let factor = xtx[j * m + i] / diag;
            for c in i..m {
                xtx[j * m + c] -= factor * xtx[i * m + c];
            }
            xty[j] -= factor * xty[i];
        }
    }

    let mut coefficients = vec![0.0; m];
    for i in (0..m).rev() {
        let mut sum = xty[i];
        for j in i + 1..m {
            sum -= xtx[i * m + j] * coefficients[j];
        }
        coefficients[i] = sum / xtx[i * m + i];
    }
    if coefficients.iter().any(|v| !v.is_finite()) {
        return None;
    }

    let rss = polynomial_rss(&coefficients, x, y);
    if rss <= s_target * 1.1 {
        Some(coefficients)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Expected values are from the executed probe of the pinned C++ sources;
    // see tests/data/spline_math_cpp_probe.tsv.
    #[test]
    fn the_scipy_default_budget_selects_a_cubic_polynomial() {
        let x = [0.0, 1.0, 2.0, 3.0, 4.0];
        let y = [0.0, 1.0, 2.0, 3.0, 4.0];
        let s = BSplineSmoothingSpline::new(&x, &y).unwrap();
        assert!(s.is_polynomial());
        assert_eq!(s.num_interior_knots(), 0);
        assert_eq!(s.smoothing_param(), 1.8377223398316205);
        assert_eq!(s.rss(), 0.0);
        assert_eq!(s.eval(-0.5).unwrap(), -0.5);
        assert_eq!(s.eval(1.5).unwrap(), 1.5);
        assert_eq!(s.eval(4.5).unwrap(), 4.5);
    }

    #[test]
    fn a_zero_budget_takes_the_interpolating_branch_which_does_not_interpolate() {
        let x = [0.0, 1.0, 2.0, 3.0, 4.0];
        let y = [0.0, 1.0, 2.0, 3.0, 4.0];
        let s = BSplineSmoothingSpline::with_smoothing(&x, &y, 0.0, 3).unwrap();
        assert!(!s.is_polynomial());
        assert_eq!(s.num_interior_knots(), 3);
        assert_eq!(s.smoothing_param(), 0.0);
        assert_eq!(s.rss(), 1.6094558499620675e-18);
        assert_eq!(s.eval(-0.5).unwrap(), 0.39062498970116977);
        assert_eq!(s.eval(0.0).unwrap(), -7.763834020124705e-11);
        assert_eq!(s.eval(4.5).unwrap(), 4.588541666767149);
    }

    #[test]
    fn a_polynomial_within_the_ten_percent_margin_wins_over_a_closer_spline() {
        let x = [0.0, 1.0, 2.0, 3.0, 4.0];
        let y = [1.0, 3.0, 2.0, 4.0, 3.5];
        let s = BSplineSmoothingSpline::with_smoothing(&x, &y, 2.0, 3).unwrap();
        assert!(s.is_polynomial());
        assert_eq!(s.rss(), 1.8892857142857145);
        assert_eq!(s.eval(0.0).unwrap(), 1.164285714285722);
        assert_eq!(s.eval(4.0).unwrap(), 3.6642857142857106);

        let interp = BSplineSmoothingSpline::with_smoothing(&x, &y, 0.0, 3).unwrap();
        assert_eq!(interp.num_interior_knots(), 3);
        assert_eq!(interp.rss(), 0.003749027381497478);
        assert_eq!(interp.eval(1.0).unwrap(), 2.968598081229554);
    }

    #[test]
    fn a_budget_no_polynomial_can_meet_runs_the_node_search() {
        let x: Vec<f64> = (0..12).map(|i| i as f64).collect();
        let y: Vec<f64> = (0..12).map(|i| f64::from(i % 2)).collect();

        // The cubic polynomial misses 0.5 by a wide margin, so the candidate
        // list is built and the densest grid, num_nodes = n, wins.
        let searched = BSplineSmoothingSpline::with_smoothing(&x, &y, 0.5, 3).unwrap();
        assert!(!searched.is_polynomial());
        assert_eq!(searched.num_interior_knots(), 10);
        assert_eq!(searched.rss(), 0.0018445489048118143);
        assert_eq!(searched.eval(0.0).unwrap(), 0.0027631640190408646);
        assert_eq!(searched.eval(2.5).unwrap(), 0.5175076479989152);
        assert_eq!(searched.eval(11.0).unwrap(), 0.9972368359809591);
        // Far outside the node domain the B-spline decays to the fitted mean.
        assert_eq!(searched.eval(13.0).unwrap(), 0.5);

        // With a budget of 3 the same polynomial is inside the margin.
        let relaxed = BSplineSmoothingSpline::with_smoothing(&x, &y, 3.0, 3).unwrap();
        assert!(relaxed.is_polynomial());
        assert_eq!(relaxed.num_interior_knots(), 0);
        assert_eq!(relaxed.rss(), 2.784770784770785);
        assert_eq!(relaxed.eval(-1.0).unwrap(), -0.13131313131316444);
        assert_eq!(relaxed.eval(13.0).unwrap(), 1.609168609168652);
    }

    #[test]
    fn small_datasets_fall_back_to_the_interpolating_branch() {
        let two = BSplineSmoothingSpline::new(&[0.0, 1.0], &[0.0, 1.0]).unwrap();
        assert_eq!(two.num_interior_knots(), 0);
        // n = 2 makes the scipy default exactly zero.
        assert_eq!(two.smoothing_param(), 0.0);
        assert_eq!(two.rss(), 0.0033376576827463334);
        assert_eq!(two.eval(0.0).unwrap(), 0.04085130158725875);
        assert_eq!(two.eval(0.5).unwrap(), 0.5);

        let three = BSplineSmoothingSpline::new(&[0.0, 1.0, 2.0], &[1.0, 2.0, 3.0]).unwrap();
        assert_eq!(three.num_interior_knots(), 1);
        assert_eq!(three.smoothing_param(), 0.5505102572168221);
        assert_eq!(three.eval(0.0).unwrap(), 1.0000000001835452);
        assert_eq!(three.eval(2.0).unwrap(), 2.9999999998164566);
    }

    #[test]
    fn the_degree_argument_changes_nothing() {
        let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [1.0, 2.0, 1.5, 3.0, 2.5, 4.0];
        let mut values = Vec::new();
        for degree in 1..=3 {
            let s = BSplineSmoothingSpline::with_smoothing(&x, &y, -1.0, degree).unwrap();
            assert_eq!(s.degree(), degree);
            assert_eq!(s.rss(), 0.9623015873015867);
            values.push(s.eval(2.0).unwrap());
        }
        assert_eq!(values, vec![2.0793650793650973; 3]);
    }

    #[test]
    fn rejects_the_inputs_the_source_reports_as_not_ok() {
        let x = [0.0, 1.0, 2.0, 3.0, 4.0];
        assert!(BSplineSmoothingSpline::new(&x, &[1.0, 2.0]).is_err());
        assert!(BSplineSmoothingSpline::new(&[0.0], &[1.0]).is_err());
        assert!(BSplineSmoothingSpline::new(&[0.0, 2.0, 1.0, 3.0], &[1.0, 2.0, 3.0, 4.0]).is_err());
        assert!(BSplineSmoothingSpline::new(&[1.0, 1.0, 2.0, 3.0], &[1.0, 2.0, 3.0, 4.0]).is_err());
        assert!(BSplineSmoothingSpline::with_smoothing(&x, &x, f64::NAN, 3).is_err());
    }

    #[test]
    fn the_candidate_comparator_is_not_a_strict_weak_ordering() {
        // Three candidates, budget 1.0, listed the way the search builds them:
        // ascending node count, so ascending interior-knot count. The finer
        // grids fit worse here, which is what makes the relation cyclic and is
        // exactly the shape the `smooth_cycle` probe case has.
        let a = (0.9975, 2); // furthest from the budget, fewest knots
        let b = (0.9983, 4);
        let c = (0.9991, 6); // closest to the budget, most knots
        let s = 1.0;
        // Neighbours are within 0.001 of each other, so the knot count decides
        // and the coarser candidate wins each pair.
        assert!(candidate_precedes(a, b, s));
        assert!(candidate_precedes(b, c, s));
        // The outer pair is 0.0016 apart, so closeness decides instead, and it
        // reverses the relation the other two imply. That is the cycle.
        assert!(!candidate_precedes(a, c, s));
        assert!(candidate_precedes(c, a, s));
        // Position zero is therefore a property of the sort, not of the
        // comparator. This port's insertion sort leaves the first element that
        // nothing displaces.
        assert_eq!(best_candidate(&[a, b, c], s), 0);
        // And it depends on the order the candidates arrive in, which is why
        // the node-count list is sorted before the splines are built.
        assert_eq!(best_candidate(&[c, b, a], s), 1);
    }

    #[test]
    fn the_knot_count_branch_overrides_closeness_only_inside_the_budget() {
        // Both inside the budget and 0.0004 apart: the coarser candidate wins
        // even though the finer one is nearer the budget.
        let coarse = (0.9990, 2);
        let fine = (0.9994, 6);
        assert_eq!(best_candidate(&[coarse, fine], 1.0), 0);
        // Move the budget below the finer candidate's residual and the branch
        // can no longer fire, so plain closeness takes over and it wins.
        assert_eq!(best_candidate(&[coarse, fine], 0.9993), 1);
        // Inside the budget but more than 0.001 apart: closeness again.
        assert_eq!(best_candidate(&[(0.9980, 2), (0.9995, 6)], 1.0), 1);
    }

    #[test]
    fn repeated_construction_is_deterministic() {
        let x = [0.0, 1.0, 2.0, 3.0, 4.0];
        let y = [1.0, 2.0, 3.0, 4.0, 5.0];
        let a = BSplineSmoothingSpline::with_smoothing(&x, &y, 2.0, 3).unwrap();
        let b = BSplineSmoothingSpline::with_smoothing(&x, &y, 2.0, 3).unwrap();
        assert_eq!(a.rss(), b.rss());
        assert_eq!(a.num_interior_knots(), b.num_interior_knots());
        for &xi in &x {
            assert_eq!(a.eval(xi).unwrap(), b.eval(xi).unwrap());
        }
    }
}
