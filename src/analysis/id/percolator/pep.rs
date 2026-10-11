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

use crate::math::fitters::levenberg_marquardt::{
    eigen_dot, eigen_gemv_row, eigen_squared_norm, eigen_sum,
};
use crate::math::libstdcxx::{std_max, std_min};
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
        if v.partial_cmp(&hi) != Some(std::cmp::Ordering::Less) {
            hi = v;
        }
    }
    let span = hi - lo;
    if span.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
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

/// The L1 data-cache size Eigen's GEMM blocking is derived from: 32 KiB, as
/// `Eigen::l1CacheSize()` reports on the reference host (AMD EPYC 7763,
/// ibminode06). The source queries it at run time with CPUID, so its
/// I-spline PEPs depend on the host once more than 504 targets are fitted.
const EIGEN_L1: i64 = 32 * 1024;
const EIGEN_L2: i64 = 512 * 1024;
/// `gebp_traits<double, double>` without FMA on x86-64: `mr = 4`, `nr = 4`.
const GEBP_MR: i64 = 4;
const GEBP_NR: i64 = 4;

/// `evaluateProductBlockingSizesHeuristic<double, double, 1>` with one
/// thread: the `(kc, mc, nc)` blocking of a `m x k` by `k x n` product.
fn gemm_blocking(k: i64, m: i64, n: i64) -> (i64, i64, i64) {
    let (mut k, mut m, mut n) = (k, m, n);
    if k.max(m.max(n)) < 48 {
        return (k, m, n);
    }
    let k_div = GEBP_MR * 8 + GEBP_NR * 8;
    let k_sub = GEBP_MR * GEBP_NR * 8;
    let max_kc = (((EIGEN_L1 - k_sub) / k_div) & !7).max(1);
    let old_k = k;
    if k > max_kc {
        k = if k % max_kc == 0 {
            max_kc
        } else {
            max_kc - 8 * ((max_kc - 1 - (k % max_kc)) / (8 * (k / max_kc + 1)))
        };
    }
    let actual_l2: i64 = 1_572_864;
    let lhs_bytes = m * k * 8;
    let remaining_l1 = EIGEN_L1 - k_sub - lhs_bytes;
    let max_nc = if remaining_l1 >= GEBP_NR * 8 * k {
        remaining_l1 / (k * 8)
    } else {
        (3 * actual_l2) / (2 * 2 * max_kc * 8)
    };
    let nc = (actual_l2 / (2 * k * 8)).min(max_nc) & !(GEBP_NR - 1);
    if n > nc {
        n = if n % nc == 0 {
            nc
        } else {
            nc - GEBP_NR * ((nc - (n % nc)) / (GEBP_NR * (n / nc + 1)))
        };
    } else if old_k == k {
        let problem_size = k * n * 8;
        let mut actual_lm = actual_l2;
        let mut max_mc = m;
        if problem_size <= 1024 {
            actual_lm = EIGEN_L1;
        } else if problem_size <= 32768 {
            actual_lm = EIGEN_L2;
            max_mc = max_mc.min(576);
        }
        let mut mc = (actual_lm / (3 * k * 8)).min(max_mc);
        if mc > GEBP_MR {
            mc -= mc % GEBP_MR;
        } else if mc == 0 {
            return (k, m, n);
        }
        m = if m % mc == 0 {
            mc
        } else {
            mc - GEBP_MR * ((mc - (m % mc)) / (GEBP_MR * (m / mc + 1)))
        };
    }
    (k, m, n)
}

/// `inv_n * (X^T X)` as the serial Eigen 5 GEMM computes it for `double`
/// on x86-64 (SSE2, no FMA): per `kc` depth panel, each entry accumulates
/// its products from zero in depth order, except that the rows the
/// one-packet micro kernel handles against a four-column block split the
/// peeled part of the panel into even and odd depths (`C` and `D`
/// accumulators, added after the peeled loop); the panel sum is then
/// scaled by `alpha` and added to the entry. Below Eigen's coefficient-based
/// threshold (`n + 2p < 20`) the product is `alpha` times the vectorised
/// sum of the products instead. Zero products are skipped, which cannot
/// change a sum of non-negative terms.
fn eigen_xtx(rows: &[Vec<f64>], p: usize, alpha: f64) -> Vec<Vec<f64>> {
    let n = rows.len();
    let mut g = vec![vec![0.0; p]; p];
    if n + 2 * p < 20 {
        if n == 0 {
            return g;
        }
        for (a, g_row) in g.iter_mut().enumerate() {
            for (b, entry) in g_row.iter_mut().enumerate() {
                *entry = alpha * eigen_sum(n, |k| rows[k][a] * rows[k][b]);
            }
        }
        return g;
    }
    let (kc, mc, nc) = gemm_blocking(n as i64, p as i64, p as i64);
    let (kc, mc, nc) = (kc.max(1) as usize, mc.max(1) as usize, nc.max(1) as usize);
    // Which entries the 1pX4 kernel computes with split accumulators.
    let mut split = vec![vec![false; p]; p];
    for i2 in (0..p).step_by(mc) {
        let r_len = mc.min(p - i2);
        let peeled_mc2 = (r_len / 4) * 4;
        let peeled_mc1 = peeled_mc2 + ((r_len - peeled_mc2) / 2) * 2;
        for j2 in (0..p).step_by(nc) {
            let c_len = nc.min(p - j2);
            let packet_cols4 = (c_len / 4) * 4;
            for r in peeled_mc2..peeled_mc1 {
                for c in 0..packet_cols4 {
                    split[i2 + r][j2 + c] = true;
                }
            }
        }
    }
    let mut acc_c = vec![vec![0.0; p]; p];
    let mut acc_d = vec![vec![0.0; p]; p];
    let mut nonzero = Vec::with_capacity(p);
    for k2 in (0..n).step_by(kc) {
        let depth = kc.min(n - k2);
        let peeled_kc = depth & !7;
        acc_c
            .iter_mut()
            .for_each(|r| r.iter_mut().for_each(|v| *v = 0.0));
        acc_d
            .iter_mut()
            .for_each(|r| r.iter_mut().for_each(|v| *v = 0.0));
        let merge = |acc_c: &mut Vec<Vec<f64>>, acc_d: &Vec<Vec<f64>>| {
            for a in 0..p {
                for b in 0..p {
                    if split[a][b] {
                        acc_c[a][b] += acc_d[a][b];
                    }
                }
            }
        };
        for local in 0..depth {
            if local == peeled_kc {
                merge(&mut acc_c, &acc_d);
            }
            let row = &rows[k2 + local];
            nonzero.clear();
            nonzero.extend((0..p).filter(|&a| row[a] != 0.0));
            let odd_peeled = local < peeled_kc && local % 2 == 1;
            for &a in &nonzero {
                for &b in &nonzero {
                    let product = row[a] * row[b];
                    if odd_peeled && split[a][b] {
                        acc_d[a][b] += product;
                    } else {
                        acc_c[a][b] += product;
                    }
                }
            }
        }
        if peeled_kc == depth {
            merge(&mut acc_c, &acc_d);
        }
        for a in 0..p {
            for b in 0..p {
                g[a][b] += acc_c[a][b] * alpha;
            }
        }
    }
    g
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
    let rows: Vec<Vec<f64>> = x_norm
        .iter()
        .map(|&v| {
            let mut row = vec![0.0; p];
            design_row(v, &knots, &mut row);
            row
        })
        .collect();
    let inv_n = 1.0 / n.max(1) as f64;
    let mut g = eigen_xtx(&rows, p, inv_n);
    // X^T y: the row-major GEMV of the transposed design, alpha = 1/n.
    let h: Vec<f64> = (0..p)
        .map(|a| 0.0 + inv_n * eigen_gemv_row(n, |k| rows[k][a], |k| y[k]))
        .collect();
    let c = inv_n * eigen_squared_norm(n, |k| y[k]);
    for (a, g_row) in g.iter_mut().enumerate() {
        g_row[a] += RIDGE_LAMBDA;
    }
    // Second-difference smoothing on the I-spline coefficients: the integer
    // matrix D^T D, scaled once and added (the product's alpha).
    let n_ispline = p - 1;
    if n_ispline > 2 {
        let nd = n_ispline - 2;
        let mut dtd = vec![vec![0.0f64; p]; p];
        for di in 0..nd {
            let cols = [(1 + di, 1.0), (2 + di, -2.0), (3 + di, 1.0)];
            for &(ca, va) in &cols {
                for &(cb, vb) in &cols {
                    dtd[ca][cb] += va * vb;
                }
            }
        }
        for a in 0..p {
            for b in 0..p {
                g[a][b] += SMOOTH_LAMBDA * dtd[a][b];
            }
        }
    }
    let mut lb = vec![0.0; p];
    lb[0] = f64::NEG_INFINITY;
    let ub = vec![f64::INFINITY; p];
    let solution = trr_boxed_qp(&g, &h, c, &lb, &ub)?;
    let mut fitted = Vec::with_capacity(n);
    for row in &rows {
        let mut value = row[0] * solution[0];
        for a in 1..p {
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

/// Eigen's vectorised `a.dot(b)`.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    eigen_dot(a.len(), |i| a[i], |i| b[i])
}

/// Eigen's `squaredNorm()`.
fn squared_norm(a: &[f64]) -> f64 {
    eigen_squared_norm(a.len(), |i| a[i])
}

/// Eigen's column-major GEMV `m * v` (`GeneralMatrixVector.h`): the
/// columns in one block below 128 columns, else in blocks of 16; per block
/// each row sums its products from zero in column order and the block sum
/// is added to the zero-initialised result.
fn mat_vec(m: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    let cols = v.len();
    let block = if cols < 128 { cols.max(1) } else { 16 };
    m.iter()
        .map(|row| {
            let mut res = 0.0;
            for j2 in (0..cols).step_by(block) {
                let mut acc = 0.0;
                for j in j2..(j2 + block).min(cols) {
                    acc += row[j] * v[j];
                }
                res += acc;
            }
            res
        })
        .collect()
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
        if squared_norm(&step).sqrt() == 0.0 {
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
        let trial: Vec<f64> = x
            .iter()
            .zip(&step)
            .map(|(a, b)| a + alpha_box * b)
            .collect();
        let (f_new, _) = cost_grad(&trial);
        let ared = f_old - f_new;
        let s: Vec<f64> = step.iter().map(|v| alpha_box * v).collect();
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
    let mut rr = squared_norm(&r);
    let mut hit_boundary = false;
    let boundary = |p_free: &[f64], d: &[f64]| {
        let a = squared_norm(d);
        let bq = 2.0 * dot(p_free, d);
        let c = squared_norm(p_free) - delta * delta;
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
        if squared_norm(&trial).sqrt() >= delta {
            let tau = boundary(&p_free, &d);
            for k in 0..nf {
                p_free[k] += tau * d[k];
            }
            hit_boundary = true;
            break;
        }
        p_free = trial;
        let r_new: Vec<f64> = r.iter().zip(&hd).map(|(a, b)| a - alpha * b).collect();
        let rr_new = squared_norm(&r_new);
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

    #[test]
    fn gemm_blocking_matches_the_probed_eigen_values() {
        // (k, m, n) -> (kc, mc, nc) as Eigen 5 reported them on ibminode06
        // (tests/data/percolator/ispline, ../oracle/I2-ispline-probe).
        assert_eq!(gemm_blocking(300, 21, 21), (300, 12, 21));
        assert_eq!(gemm_blocking(600, 28, 28), (304, 28, 28));
        assert_eq!(gemm_blocking(400, 24, 24), (400, 24, 24));
        assert_eq!(gemm_blocking(100, 14, 14), (100, 8, 14));
        assert_eq!(gemm_blocking(1500, 42, 42), (504, 42, 42));
        assert_eq!(gemm_blocking(20000, 145, 145), (504, 145, 145));
        assert_eq!(gemm_blocking(17, 8, 8), (17, 8, 8));
    }

    /// Executed probe (tier 2): the vendored `ISplineTRRRegressor`, compiled
    /// with the Release build's `OpenMS_Percolator` flags and run serially on
    /// ibminode06, against this port, bit for bit.
    #[test]
    fn ispline_pep_matches_the_serial_release_probe() {
        let dir =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/percolator/ispline");
        let read = |path: std::path::PathBuf| -> Vec<f64> {
            std::fs::read_to_string(path)
                .unwrap()
                .lines()
                .filter(|l| !l.is_empty())
                .map(|l| l.parse().unwrap())
                .collect()
        };
        for name in [
            "rescore_default",
            "rescore_nested",
            "rescore_none",
            "rescore_subsample",
            "syn_1500",
            "syn_20000",
            "syn_one",
            "syn_small17",
            "syn_tiny5",
            "tron",
        ] {
            let q = read(dir.join(format!("{name}.q")));
            let expected = read(dir.join(format!("{name}.pep")));
            let got = InferPep::new(true).q_to_pep(&q).unwrap();
            assert_eq!(got.len(), expected.len(), "{name}");
            for (i, (a, b)) in got.iter().zip(&expected).enumerate() {
                assert_eq!(a.to_bits(), b.to_bits(), "{name} row {i}: {a:e} vs {b:e}");
            }
        }
    }
}
