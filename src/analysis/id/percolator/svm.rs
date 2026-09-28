// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//
// Translated from the TRON solver of LIBLINEAR v2.11 (Copyright (c)
// 2007-2015 The LIBLINEAR Project, BSD 3-Clause License) as integrated into
// Percolator by J. Halloran and vendored by OpenMS (ssl.cpp, tron.cpp); see
// docs/PERCOLATOR_SUPPORT.md for the notice.

//! The L2-loss linear SVM the vendored Percolator trains
//! (`L2_SVM_MFN`, which the OpenMS copy forwards to TRON).
//!
//! The objective is liblinear's L2-loss SVC with per-example costs
//! `C_i = c_i / (2 lambda)` (`lambda = 1`), which has the same minimiser as
//! Percolator's original modified-finite-Newton objective. The four BLAS-1
//! helpers of the source are Eigen maps; they are reproduced with the Eigen
//! reduction kernels of [`crate::math::fitters::levenberg_marquardt`], so the
//! arithmetic follows the x86_64 Release build's order.

use crate::math::fitters::levenberg_marquardt::{eigen_dot, eigen_gemv_row, eigen_squared_norm};
use crate::{Error, Result};

/// `EPSILON`, the outer TRON tolerance of the Percolator configuration.
const EPSILON: f64 = 1e-7;
/// `MFNITERMAX`, TRON's outer iteration cap.
const MAX_ITER: usize = 50;
/// `kTronEpsCg`, the inner CG tolerance.
const EPS_CG: f64 = 0.1;

fn nrm2(x: &[f64]) -> f64 {
    eigen_squared_norm(x.len(), |i| x[i]).sqrt()
}

fn dot(x: &[f64], y: &[f64]) -> f64 {
    eigen_dot(x.len(), |i| x[i], |i| y[i])
}

fn axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    for (yi, xi) in y.iter_mut().zip(x) {
        *yi += alpha * xi;
    }
}

/// `TronFunction`: `f(w) = 0.5 ||w||^2 + sum_{margin_i < 1} C_i (1 -
/// margin_i)^2` over a row-major design with a trailing bias column.
struct TronFunction<'a> {
    x: &'a [f64],
    t: &'a [f64],
    c: &'a [f64],
    m: usize,
    n: usize,
    margin: Vec<f64>,
    active: Vec<usize>,
}

impl TronFunction<'_> {
    fn row(&self, i: usize) -> &[f64] {
        &self.x[i * self.n..(i + 1) * self.n]
    }

    fn fun(&mut self, w: &[f64]) -> f64 {
        for i in 0..self.m {
            let row = &self.x[i * self.n..(i + 1) * self.n];
            self.margin[i] = eigen_gemv_row(self.n, |j| row[j], |j| w[j]) * self.t[i];
        }
        let mut f = 0.5 * eigen_squared_norm(w.len(), |i| w[i]);
        for i in 0..self.m {
            if self.margin[i] < 1.0 {
                let d = 1.0 - self.margin[i];
                f += self.c[i] * d * d;
            }
        }
        f
    }

    fn grad(&mut self, w: &[f64], g: &mut [f64]) {
        self.active.clear();
        g.iter_mut().for_each(|v| *v = 0.0);
        for i in 0..self.m {
            if self.margin[i] < 1.0 {
                self.active.push(i);
                let z = 2.0 * self.c[i] * self.t[i] * (self.margin[i] - 1.0);
                let row = self.row(i);
                for j in 0..self.n {
                    g[j] += z * row[j];
                }
            }
        }
        for j in 0..self.n {
            g[j] += w[j];
        }
    }

    fn hv(&self, s: &[f64], hs: &mut [f64]) {
        hs.iter_mut().for_each(|v| *v = 0.0);
        for &i in &self.active {
            let row = self.row(i);
            let mut xts = 0.0;
            for j in 0..self.n {
                xts += row[j] * s[j];
            }
            let scale = 2.0 * self.c[i] * xts;
            for j in 0..self.n {
                hs[j] += scale * row[j];
            }
        }
        for j in 0..self.n {
            hs[j] += s[j];
        }
    }
}

/// liblinear's own `min` and `max` templates (`(x<y)?x:y`, `(x>y)?x:y`).
fn tmin(x: f64, y: f64) -> f64 {
    if x < y { x } else { y }
}
fn tmax(x: f64, y: f64) -> f64 {
    if x > y { x } else { y }
}

/// Train the L2-loss SVM: `rows` are the training feature rows (without the
/// bias), `labels` are +1/-1. Returns `n_features + 1` weights, bias last,
/// from a zero start.
///
/// # Errors
///
/// Native guard: [`Error::InvalidValue`] when the inner CG or the outer
/// trust-region loop fails to terminate (at most `20 n + 1000` CG steps and
/// 100,000 outer steps); the source loops on forever there, which only a
/// NaN in the data can cause.
pub(super) fn l2_svm(
    rows: &[&[f64]],
    labels: &[f64],
    n_features: usize,
    cpos: f64,
    cneg: f64,
) -> Result<Vec<f64>> {
    let m = rows.len();
    let n = n_features + 1;
    let mut x = Vec::with_capacity(m * n);
    for row in rows {
        x.extend_from_slice(&row[..n_features]);
        x.push(1.0);
    }
    let inv_2lambda = 1.0 / (2.0 * 1.0);
    let c: Vec<f64> = labels
        .iter()
        .map(|&y| (if y > 0.0 { cpos } else { cneg }) * inv_2lambda)
        .collect();
    let mut function = TronFunction {
        x: &x,
        t: labels,
        c: &c,
        m,
        n,
        margin: vec![0.0; m],
        active: Vec::with_capacity(m),
    };
    let mut w = vec![0.0; n];
    tron(&mut function, &mut w)?;
    Ok(w)
}

/// `TRON::tron`.
fn tron(function: &mut TronFunction<'_>, w: &mut [f64]) -> Result<()> {
    let (eta0, eta1, eta2) = (1e-4, 0.25, 0.75);
    let (sigma1, sigma2, sigma3) = (0.25, 0.5, 4.0);
    let n = function.n;
    let mut s = vec![0.0; n];
    let mut r = vec![0.0; n];
    let mut g = vec![0.0; n];
    let w0 = vec![0.0; n];
    function.fun(&w0);
    function.grad(&w0, &mut g);
    let gnorm0 = nrm2(&g);
    let mut f = function.fun(w);
    function.grad(w, &mut g);
    let mut delta = nrm2(&g);
    let search = !(delta <= EPSILON * gnorm0);
    let mut iter = 1usize;
    let mut w_new = vec![0.0; n];
    let mut steps = 0usize;
    while iter <= MAX_ITER && search {
        steps += 1;
        if steps > 100_000 {
            return Err(Error::InvalidValue(
                "TRON trust-region loop does not terminate".into(),
            ));
        }
        let reach_boundary = trcg(function, delta, &g, &mut s, &mut r)?;
        w_new.copy_from_slice(w);
        axpy(1.0, &s, &mut w_new);
        let gs = dot(&g, &s);
        let prered = -0.5 * (gs - dot(&s, &r));
        let fnew = function.fun(&w_new);
        let actred = f - fnew;
        let snorm = nrm2(&s);
        if iter == 1 {
            delta = tmin(delta, snorm);
        }
        let alpha = if fnew - f - gs <= 0.0 {
            sigma3
        } else {
            tmax(sigma1, -0.5 * (gs / (fnew - f - gs)))
        };
        if actred < eta0 * prered {
            delta = tmin(tmax(alpha, sigma1) * snorm, sigma2 * delta);
        } else if actred < eta1 * prered {
            delta = tmax(sigma1 * delta, tmin(alpha * snorm, sigma2 * delta));
        } else if actred < eta2 * prered {
            delta = tmax(sigma1 * delta, tmin(alpha * snorm, sigma3 * delta));
        } else if reach_boundary {
            delta *= sigma3;
        } else {
            delta = tmax(delta, tmin(alpha * snorm, sigma3 * delta));
        }
        if actred > eta0 * prered {
            iter += 1;
            w.copy_from_slice(&w_new);
            f = fnew;
            function.grad(w, &mut g);
            if nrm2(&g) <= EPSILON * gnorm0 {
                break;
            }
        }
        if f < -1.0e32 {
            break;
        }
        if prered <= 0.0 {
            break;
        }
        if actred.abs() <= 1.0e-12 * f.abs() && prered.abs() <= 1.0e-12 * f.abs() {
            break;
        }
    }
    Ok(())
}

/// `TRON::trcg`: CG on the Newton system inside the trust region; whether
/// the step reached the boundary.
fn trcg(
    function: &TronFunction<'_>,
    delta: f64,
    g: &[f64],
    s: &mut [f64],
    r: &mut [f64],
) -> Result<bool> {
    let n = function.n;
    let mut d = vec![0.0; n];
    let mut hd = vec![0.0; n];
    for i in 0..n {
        s[i] = 0.0;
        r[i] = -g[i];
        d[i] = r[i];
    }
    let cgtol = EPS_CG * nrm2(g);
    let mut rtr = dot(r, r);
    let cap = 20 * n + 1000;
    let mut iterations = 0usize;
    loop {
        if nrm2(r) <= cgtol {
            break;
        }
        iterations += 1;
        if iterations > cap {
            return Err(Error::InvalidValue(
                "TRON conjugate gradient does not terminate".into(),
            ));
        }
        function.hv(&d, &mut hd);
        let mut alpha = rtr / dot(&d, &hd);
        axpy(alpha, &d, s);
        if nrm2(s) > delta {
            alpha = -alpha;
            axpy(alpha, &d, s);
            let std = dot(s, &d);
            let sts = dot(s, s);
            let dtd = dot(&d, &d);
            let dsq = delta * delta;
            let rad = (std * std + dtd * (dsq - sts)).sqrt();
            alpha = if std >= 0.0 {
                (dsq - sts) / (std + rad)
            } else {
                (rad - std) / dtd
            };
            axpy(alpha, &d, s);
            alpha = -alpha;
            axpy(alpha, &hd, r);
            return Ok(true);
        }
        alpha = -alpha;
        axpy(alpha, &hd, r);
        let rnew = dot(r, r);
        let beta = rnew / rtr;
        for di in d.iter_mut() {
            *di *= beta;
        }
        axpy(1.0, r, &mut d);
        rtr = rnew;
    }
    Ok(false)
}
