// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//
// Translated from Percolator (https://github.com/percolator/percolator,
// commit eb157f7, as vendored by OpenMS), Copyright 2006-2012 Lukas Käll and
// the Percolator authors, licensed under the Apache License, Version 2.0; see
// docs/PERCOLATOR_SUPPORT.md for the notice.

//! Posterior error probabilities from q-values: the vendored `InferPEP`
//! (`IsotonicPEP.cpp`) with its two monotone regressors, `PAVARegressor`
//! (pool-adjacent-violators) and `ISplineTRRRegressor` (a non-negative ridge
//! fit of I-spline coefficients, solved by a trust-region reflective method).
//!
//! The source solves the I-spline problem with Eigen dense products
//! (`X^T X` is a blocked GEMM whose block sizes Eigen derives from the host's
//! cache sizes); this port uses plain in-order sums, so the I-spline PEPs
//! agree with the C++ to rounding, not to the bit. The PAVA path is exact.

use super::engine::{std_max, std_min};
use crate::{Error, Result};

/// `InferPEP`: q-values to PEPs by monotone regression.
pub(super) struct InferPep {
    use_ispline: bool,
}

impl InferPep {
    /// `InferPEP(use_ispline)`: the I-spline regressor, or PAVA.
    pub(super) fn new(use_ispline: bool) -> Self {
        Self { use_ispline }
    }

    /// `q_to_pep`: raw PEPs are the increments of `q_i * (i + 1)`, plus a
    /// pseudocount of half a false discovery spread over all entries, fitted
    /// monotonically with a lower clip of `1e-10`.
    pub(super) fn q_to_pep(&self, q_values: &[f64]) -> Result<Vec<f64>> {
        let n = q_values.len();
        let mut raw = vec![0.0; n];
        let mut previous = 0.0;
        for i in 0..n {
            let qn = q_values[i] * (i + 1) as f64;
            raw[i] = qn - previous;
            previous = qn;
        }
        if n > 0 {
            let pseudo = 0.5 / n as f64;
            for value in &mut raw {
                *value += pseudo;
            }
        }
        const EPS: f64 = 1e-10;
        if self.use_ispline {
            ispline_fit_y(&raw, EPS, 1.0)
        } else {
            Ok(pava(&raw, EPS, 1.0))
        }
    }

    /// `InferPEP::interpolate`: linear in q between two known PEPs, clamped
    /// to `[0, 1]`; the right PEP when the q-values coincide.
    pub(super) fn interpolate(q: f64, q1: f64, q2: f64, pep1: f64, pep2: f64) -> f64 {
        if (q2 - q1).abs() <= 1e-12 {
            return std_max(0.0, std_min(1.0, pep2));
        }
        let value = pep1 + (q - q1) * (pep2 - pep1) / (q2 - q1);
        std_max(0.0, std_min(1.0, value))
    }
}

/// `PAVARegressor::pava_impl`: non-decreasing isotonic fit by merging
/// adjacent violators, then clipped to `[min, max]`.
pub(super) fn pava(values: &[f64], min_value: f64, max_value: f64) -> Vec<f64> {
    struct Block {
        sum: f64,
        count: i64,
        avg: f64,
    }
    let mut stack: Vec<Block> = Vec::with_capacity(values.len());
    for &value in values {
        stack.push(Block {
            sum: value,
            count: 1,
            avg: value,
        });
        while stack.len() > 1 {
            let top = &stack[stack.len() - 1];
            let second = &stack[stack.len() - 2];
            let sum = second.sum + top.sum;
            let count = second.count + top.count;
            let avg = sum / count as f64;
            if second.avg > top.avg {
                stack.pop();
                stack.pop();
                stack.push(Block { sum, count, avg });
            } else {
                break;
            }
        }
    }
    let mut result = Vec::with_capacity(values.len());
    for block in &stack {
        let value = std_max(std_min(block.avg, max_value), min_value);
        for _ in 0..block.count {
            result.push(value);
        }
    }
    result
}

/// `MonotoneParams` defaults: ridge `1e-4`, smoothing `1e-3`, cubic
/// I-splines with an intercept in column 0, y decreasing in x.
const RIDGE_LAMBDA: f64 = 1e-4;
const SMOOTH_LAMBDA: f64 = 1e-3;
const DEGREE: usize = 3;

/// `ISplineTRRRegressor::fit_y`: x is the reversed rank `N, ..., 1`.
fn ispline_fit_y(y: &[f64], clip_lo: f64, clip_hi: f64) -> Result<Vec<f64>> {
    let n = y.len();
    let x: Vec<f64> = (0..n).map(|i| (n - i) as f64).collect();
    ispline_fit_xy(&x, y, clip_lo, clip_hi)
}

/// `ispline_detail::normalize_to_unit_interval`.
fn normalize_to_unit_interval(x: &[f64]) -> Vec<f64> {
    if x.is_empty() {
        return Vec::new();
    }
    let (mut lo, mut hi) = (x[0], x[0]);
    for &v in &x[1..] {
        if v < lo {
            lo = v;
        }
        if !(v < hi) {
            hi = v;
        }
    }
    let span = hi - lo;
    if !(span > 0.0) {
        return vec![0.5; x.len()];
    }
    x.iter().map(|v| (v - lo) / span).collect()
}

/// `make_default_knots`: a clamped cubic knot vector with at most 200
/// internal knots at data quantiles.
fn default_knots(scores: &[f64]) -> Vec<f64> {
    let n = scores.len();
    if n == 0 {
        return Vec::new();
    }
    let mut sorted = scores.to_vec();
    sorted.sort_by(f64::total_cmp);
    let lo = sorted[0];
    let mut hi = sorted[n - 1];
    if hi <= lo {
        hi = lo + 1e-6;
    }
    let num_internal = 200.min((n as f64).sqrt() as i64) as usize;
    let mut internal = Vec::with_capacity(num_internal);
    let mut previous = lo;
    for k in 1..=num_internal {
        let q = k as f64 / (num_internal + 1) as f64;
        let index = ((q * (n - 1) as f64) as usize).min(n - 1);
        let value = sorted[index];
        if value > lo + 1e-12 && value < hi - 1e-12 && value > previous + 1e-12 {
            internal.push(value);
            previous = value;
        }
    }
    let order = DEGREE + 1;
    let mut knots = Vec::with_capacity(2 * order + internal.len());
    knots.extend(std::iter::repeat_n(lo, order));
    knots.extend(internal);
    knots.extend(std::iter::repeat_n(hi, order));
    knots
}

/// One row of `build_ispline_design`: the intercept, then the I-spline
/// columns `I_j(x) = sum_{m >= j} B_m(x)` from the de Boor recursion.
fn design_row(x: f64, knots: &[f64], row: &mut [f64]) {
    let n_knots = knots.len();
    let d = DEGREE;
    let n_bspline = n_knots as i64 - d as i64 - 1;
    if n_bspline <= 1 {
        row[0] = 1.0;
        return;
    }
    let n_bspline = n_bspline as usize;
    row[0] = 1.0;
    let t = knots;
    let k = if x >= t[n_knots - d - 1] {
        n_knots - d - 2
    } else if x <= t[d] {
        d
    } else {
        let (mut lo, mut hi) = (d, n_knots - d - 1);
        while lo < hi {
            let mid = (lo + hi) / 2;
            if t[mid + 1] <= x {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    };
    let mut b = [0.0f64; DEGREE + 1];
    b[0] = 1.0;
    for j in 1..=d {
        let mut saved = 0.0;
        for r in 0..j {
            let left = x - t[k + 1 + r - j];
            let right = t[k + 1 + r] - x;
            let denom = left + right;
            let temp = if denom > 0.0 { b[r] / denom } else { 0.0 };
            b[r] = saved + right * temp;
            saved = left * temp;
        }
        b[j] = saved;
    }
    let mut cum = [0.0f64; DEGREE + 1];
    cum[d] = b[d];
    for r in (0..d).rev() {
        cum[r] = cum[r + 1] + b[r];
    }
    for j in 1..n_bspline {
        let value = if j + d <= k {
            1.0
        } else if j > k {
            0.0
        } else {
            cum[j + d - k]
        };
        row[j] = value;
    }
}

/// Number of design columns `build_ispline_design` produces.
fn design_columns(knots: &[f64]) -> usize {
    let n_bspline = knots.len() as i64 - DEGREE as i64 - 1;
    if n_bspline <= 1 {
        1
    } else {
        n_bspline as usize
    }
}

/// `ISplineTRRRegressor::fit_xy` with the default parameters.
fn ispline_fit_xy(x: &[f64], y: &[f64], clip_lo: f64, clip_hi: f64) -> Result<Vec<f64>> {
    if x.len() != y.len() {
        return Err(Error::InvalidValue(
            "I-spline regression needs one x per y".into(),
        ));
    }
    let n = x.len();
    // y decreases in x: fit against -x.
    let negated: Vec<f64> = x.iter().map(|v| -v).collect();
    let x_norm = normalize_to_unit_interval(&negated);
    let knots = default_knots(&x_norm);
    let p = design_columns(&knots);
    let inv_n = 1.0 / n.max(1) as f64;
    let mut g = vec![vec![0.0; p]; p];
    let mut h = vec![0.0; p];
    let mut c = 0.0;
    let mut row = vec![0.0; p];
    for i in 0..n {
        row.iter_mut().for_each(|v| *v = 0.0);
        design_row(x_norm[i], &knots, &mut row);
        for a in 0..p {
            let xa = row[a];
            if xa == 0.0 {
                continue;
            }
            for b in 0..p {
                g[a][b] += xa * row[b];
            }
            h[a] += xa * y[i];
        }
        c += y[i] * y[i];
    }
    for a in 0..p {
        for b in 0..p {
            g[a][b] *= inv_n;
        }
        h[a] *= inv_n;
    }
    c *= inv_n;
    for (a, g_row) in g.iter_mut().enumerate() {
        g_row[a] += RIDGE_LAMBDA;
    }
    // Second-difference smoothing on the I-spline coefficients.
    let n_ispline = p - 1;
    if n_ispline > 2 {
        let nd = n_ispline - 2;
        for di in 0..nd {
            let cols = [(1 + di, 1.0), (2 + di, -2.0), (3 + di, 1.0)];
            for &(ca, va) in &cols {
                for &(cb, vb) in &cols {
                    g[ca][cb] += SMOOTH_LAMBDA * (va * vb);
                }
            }
        }
    }
    let mut lb = vec![0.0; p];
    lb[0] = f64::NEG_INFINITY;
    let ub = vec![f64::INFINITY; p];
    let solution = trr_boxed_qp(&g, &h, c, &lb, &ub)?;
    let mut fitted = Vec::with_capacity(n);
    for i in 0..n {
        row.iter_mut().for_each(|v| *v = 0.0);
        design_row(x_norm[i], &knots, &mut row);
        let mut value = 0.0;
        for a in 0..p {
            value += row[a] * solution[a];
        }
        if value < clip_lo {
            value = clip_lo;
        }
        if value > clip_hi {
            value = clip_hi;
        }
        fitted.push(value);
    }
    Ok(fitted)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn mat_vec(m: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    m.iter().map(|row| dot(row, v)).collect()
}

/// `trr_boxed_qp`: minimise `x^T G x - 2 h^T x + c` subject to
/// `lb <= x <= ub` from `x = 0`, with Steihaug-Toint CG steps on the free
/// variables; at most 1000 iterations, projected-gradient tolerance `1e-8`.
fn trr_boxed_qp(g: &[Vec<f64>], h: &[f64], c: f64, lb: &[f64], ub: &[f64]) -> Result<Vec<f64>> {
    const DELTA_MAX: f64 = 1e6;
    const ETA: f64 = 1e-3;
    const GTOL: f64 = 1e-8;
    const MAX_ITERS: usize = 1000;
    let p = h.len();
    let mut x = vec![0.0; p];
    for i in 0..p {
        if x[i] < lb[i] {
            x[i] = lb[i];
        } else if x[i] > ub[i] {
            x[i] = ub[i];
        }
    }
    let mut delta = 1.0;
    let hess: Vec<Vec<f64>> = g
        .iter()
        .map(|row| row.iter().map(|v| 2.0 * v).collect())
        .collect();
    let cost_grad = |x: &[f64]| {
        let gx = mat_vec(g, x);
        let f = dot(x, &gx) - 2.0 * dot(h, x) + c;
        let grad: Vec<f64> = gx.iter().zip(h).map(|(a, b)| 2.0 * (a - b)).collect();
        (f, grad)
    };
    for _ in 0..MAX_ITERS {
        let (f_old, grad) = cost_grad(&x);
        let mut norm_inf: f64 = 0.0;
        for i in 0..p {
            let mut gi = grad[i];
            if x[i] <= lb[i] + 1e-12 {
                gi = std_min(0.0, gi);
            } else if x[i] >= ub[i] - 1e-12 {
                gi = std_max(0.0, gi);
            }
            norm_inf = std_max(norm_inf, gi.abs());
        }
        if norm_inf < GTOL {
            return Ok(x);
        }
        let free: Vec<usize> = (0..p)
            .filter(|&i| {
                let at_l = (x[i] - lb[i]).abs() <= 1e-12;
                let at_u = (x[i] - ub[i]).abs() <= 1e-12;
                (!at_l && !at_u) || (at_l && grad[i] < 0.0) || (at_u && grad[i] > 0.0)
            })
            .collect();
        if free.is_empty() {
            return Ok(x);
        }
        let (step, hit_boundary) = steihaug_cg(&hess, &grad, &free, delta, p);
        if dot(&step, &step).sqrt() == 0.0 {
            return Ok(x);
        }
        let mut alpha = 1.0;
        for i in 0..p {
            if step[i] > 0.0 {
                alpha = std_min(alpha, (ub[i] - x[i]) / step[i]);
            } else if step[i] < 0.0 {
                alpha = std_min(alpha, (lb[i] - x[i]) / step[i]);
            }
        }
        let alpha_box = std_max(0.0, alpha);
        let s: Vec<f64> = step.iter().map(|v| alpha_box * v).collect();
        let trial: Vec<f64> = x.iter().zip(&s).map(|(a, b)| a + b).collect();
        let (f_new, _) = cost_grad(&trial);
        let ared = f_old - f_new;
        let hs = mat_vec(&hess, &s);
        let pred = -(dot(&grad, &s) + 0.5 * dot(&s, &hs));
        let rho = if pred > 0.0 {
            ared / pred
        } else {
            f64::NEG_INFINITY
        };
        if rho < 0.25 {
            delta *= 0.25;
        } else if rho > 0.75 && (alpha_box - 1.0).abs() < 1e-12 && hit_boundary {
            delta = std_min(2.0 * delta, DELTA_MAX);
        }
        if rho >= ETA {
            x = trial;
        }
    }
    Ok(x)
}

/// `steihaugCG`: truncated CG on the free coordinates within `||p|| <= delta`;
/// the full-space step and whether it reached the boundary.
fn steihaug_cg(
    hess: &[Vec<f64>],
    grad: &[f64],
    free: &[usize],
    delta: f64,
    p: usize,
) -> (Vec<f64>, bool) {
    let hmul = |v_free: &[f64]| -> Vec<f64> {
        let mut v = vec![0.0; p];
        for (k, &i) in free.iter().enumerate() {
            v[i] = v_free[k];
        }
        let hv = mat_vec(hess, &v);
        free.iter().map(|&i| hv[i]).collect()
    };
    let nf = free.len();
    let g_free: Vec<f64> = free.iter().map(|&i| grad[i]).collect();
    let mut p_free = vec![0.0; nf];
    let mut r: Vec<f64> = g_free.iter().map(|v| -v).collect();
    let mut d = r.clone();
    let mut rr = dot(&r, &r);
    let mut hit_boundary = false;
    let boundary = |p_free: &[f64], d: &[f64]| {
        let a = dot(d, d);
        let bq = 2.0 * dot(p_free, d);
        let c = dot(p_free, p_free) - delta * delta;
        (-bq + std_max(0.0, bq * bq - 4.0 * a * c).sqrt()) / (2.0 * a)
    };
    let mut iteration = 0;
    while iteration < nf && rr > 1e-12 {
        iteration += 1;
        let hd = hmul(&d);
        let d_hd = dot(&d, &hd);
        if d_hd <= 0.0 {
            let tau = boundary(&p_free, &d);
            for k in 0..nf {
                p_free[k] += tau * d[k];
            }
            hit_boundary = true;
            break;
        }
        let alpha = rr / d_hd;
        let trial: Vec<f64> = p_free.iter().zip(&d).map(|(a, b)| a + alpha * b).collect();
        if dot(&trial, &trial).sqrt() >= delta {
            let tau = boundary(&p_free, &d);
            for k in 0..nf {
                p_free[k] += tau * d[k];
            }
            hit_boundary = true;
            break;
        }
        p_free = trial;
        let r_new: Vec<f64> = r.iter().zip(&hd).map(|(a, b)| a - alpha * b).collect();
        let rr_new = dot(&r_new, &r_new);
        let beta = rr_new / rr;
        d = r_new.iter().zip(&d).map(|(a, b)| a + beta * b).collect();
        r = r_new;
        rr = rr_new;
    }
    let mut step = vec![0.0; p];
    for (k, &i) in free.iter().enumerate() {
        step[i] = p_free[k];
    }
    (step, hit_boundary)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pava_pools_violators() {
        assert_eq!(pava(&[1.0, 3.0, 2.0, 4.0], 0.0, 10.0), [1.0, 2.5, 2.5, 4.0]);
        assert_eq!(pava(&[0.5, 0.2], 0.0, 1.0), [0.35, 0.35]);
        assert!(pava(&[], 0.0, 1.0).is_empty());
    }

    #[test]
    fn ispline_design_rows_are_monotone_and_bounded() {
        let knots = default_knots(&(0..100).map(|i| i as f64 / 99.0).collect::<Vec<_>>());
        let p = design_columns(&knots);
        let mut row = vec![0.0; p];
        design_row(0.37, &knots, &mut row);
        assert_eq!(row[0], 1.0);
        for j in 2..p {
            assert!(row[j] <= row[j - 1] + 1e-15);
            assert!((0.0..=1.0 + 1e-15).contains(&row[j]));
        }
    }

    #[test]
    fn ispline_fit_is_monotone_in_rank() {
        let q: Vec<f64> = (0..50).map(|i| f64::from(i) / 100.0).collect();
        let pep = InferPep::new(true).q_to_pep(&q).unwrap();
        assert_eq!(pep.len(), 50);
        for w in pep.windows(2) {
            assert!(w[1] >= w[0] - 1e-12);
        }
    }
}
