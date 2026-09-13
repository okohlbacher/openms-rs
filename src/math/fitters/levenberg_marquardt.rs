// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Levenberg-Marquardt least-squares solver reproducing the one the C++
//! distribution fitters call.
//!
//! `GaussFitter.cpp`, `GammaDistributionFitter.cpp`,
//! `GumbelDistributionFitter.cpp` and `GumbelMaxLikelihoodFitter.cpp` each
//! construct an `Eigen::LevenbergMarquardt<Functor>` and call `minimize`.
//! Eigen's solver is a transcription of MINPACK `lmder`: column-pivoted
//! Householder QR of the Jacobian, the `lmpar` trust-region parameter search
//! with a Givens-rotation `qrsolv`, and MINPACK's termination tests. Because
//! the crate may not take a dependency on a non-linear optimizer, that
//! algorithm is reproduced here, step for step and in the same arithmetic
//! order, so that the fitters converge to the published parameters rather than
//! to some other point of a non-convex surface.
//!
//! Defaults match `Eigen::LevenbergMarquardt::Parameters`:
//! `factor = 100`, `maxfev = 400`, `ftol = xtol = sqrt(f64::EPSILON)`,
//! `gtol = 0`, and no external scaling. See
//! [`LmParameters`](crate::math::fitters::levenberg_marquardt::LmParameters).
//!
//! The solver itself performs no allocation proportional to anything but the
//! problem it is handed; the callers preflight their point counts through
//! [`preflight_points`](crate::math::fitters::levenberg_marquardt::preflight_points).
//!
//! See `docs/DISTRIBUTION_FITTERS_SUPPORT.md`.

// The loops below are a line-by-line transcription of MINPACK's `lmder`,
// `lmpar` and `qrsolv` and of Eigen's column-pivoted QR. Their index ranges are
// correlated across several arrays and a matrix at once - `for i in (k+1)..n`
// walks `sdiag`, `s(i, k)` and `wa` together - so rewriting them as iterator
// chains would hide, not clarify, the correspondence with the reference
// algorithm that this module exists to preserve.
#![allow(clippy::needless_range_loop)]

use crate::{Error, Result};

/// Largest number of (x, y) observations a fitter accepts.
///
/// The C++ fitters have no ceiling; the residual and Jacobian buffers scale
/// with the observation count, so one is imposed here and checked before any
/// allocation.
pub const MAX_POINTS: usize = 1_000_000;

/// Largest total dense-Jacobian allocation a fit may request, in bytes.
pub const MAX_BYTES: usize = 64 * 1024 * 1024;

/// Reject an observation count or Jacobian size beyond the fitters' ceilings.
///
/// `values` is the number of residuals (one per observation) and `inputs` the
/// number of fitted parameters. Nothing is allocated before this returns.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when more than [`MAX_POINTS`] observations
/// are supplied or when the dense Jacobian would exceed [`MAX_BYTES`].
pub fn preflight_points(values: usize, inputs: usize) -> Result<()> {
    if values > MAX_POINTS {
        return Err(limit());
    }
    let bytes = values
        .checked_mul(inputs)
        .and_then(|n| n.checked_mul(std::mem::size_of::<f64>()))
        .ok_or_else(limit)?;
    if bytes > MAX_BYTES {
        return Err(limit());
    }
    Ok(())
}

fn limit() -> Error {
    Error::InvalidValue("fitter input exceeds its point or byte limit".into())
}

/// Tuning constants of the iteration, named as in `Eigen::LevenbergMarquardt`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LmParameters {
    /// Initial trust-region scale: `delta = factor * ||diag * x||`.
    pub factor: f64,
    /// Maximum number of residual evaluations before the iteration gives up.
    pub max_fev: usize,
    /// Relative reduction of the residual norm below which the fit is accepted.
    pub ftol: f64,
    /// Relative parameter change below which the fit is accepted.
    pub xtol: f64,
    /// Scaled gradient norm below which the fit is accepted.
    pub gtol: f64,
}

impl Default for LmParameters {
    /// `factor = 100`, `max_fev = 400`, `ftol = xtol = sqrt(f64::EPSILON)`,
    /// `gtol = 0`, exactly as Eigen's defaults, which every OpenMS fitter uses
    /// because none of them touches `lmSolver.parameters`.
    fn default() -> Self {
        Self {
            factor: 100.0,
            max_fev: 400,
            ftol: f64::EPSILON.sqrt(),
            xtol: f64::EPSILON.sqrt(),
            gtol: 0.0,
        }
    }
}

/// Why the iteration stopped, mirroring `Eigen::LevenbergMarquardtSpace::Status`.
///
/// The C++ fitters branch on the numeric value of this enum, so
/// [`LmStatus::code`](crate::math::fitters::levenberg_marquardt::LmStatus::code)
/// exposes it. `NotStarted` (-2), `Running` (-1) and `UserAsked` (9) cannot be
/// observed here: `minimize` loops until it is no longer running and the
/// residual closures cannot abort.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LmStatus {
    /// Dimensions or tolerances are unusable, in particular fewer residuals
    /// than parameters. Eigen's value 0.
    ImproperInputParameters,
    /// Both the actual and the predicted reduction fell below `ftol`. Value 1.
    RelativeReductionTooSmall,
    /// The trust region shrank below `xtol * ||diag * x||`. Value 2.
    RelativeErrorTooSmall,
    /// Both of the two preceding tests passed at once. Value 3.
    RelativeErrorAndReductionTooSmall,
    /// The scaled gradient norm fell to `gtol`. Value 4.
    CosinusTooSmall,
    /// `max_fev` residual evaluations were spent. Value 5.
    TooManyFunctionEvaluation,
    /// `ftol` is too small to be meaningful at this precision. Value 6.
    FtolTooSmall,
    /// `xtol` is too small to be meaningful at this precision. Value 7.
    XtolTooSmall,
    /// `gtol` is too small to be meaningful at this precision. Value 8.
    GtolTooSmall,
}

impl LmStatus {
    /// The numeric value Eigen gives this status.
    ///
    /// `GaussFitter` rejects 0 and 5; the other three fitters reject anything
    /// `<= 0`, which after a completed `minimize` can only be 0.
    pub fn code(self) -> i32 {
        match self {
            Self::ImproperInputParameters => 0,
            Self::RelativeReductionTooSmall => 1,
            Self::RelativeErrorTooSmall => 2,
            Self::RelativeErrorAndReductionTooSmall => 3,
            Self::CosinusTooSmall => 4,
            Self::TooManyFunctionEvaluation => 5,
            Self::FtolTooSmall => 6,
            Self::XtolTooSmall => 7,
            Self::GtolTooSmall => 8,
        }
    }
}

/// Dense row-major `f64` matrix, the solver's Jacobian and QR storage.
///
/// Only the operations the algorithm needs are provided. Element access is
/// bounds-checked; every index the solver forms is derived from the row and
/// column counts it allocated, so a failure would be an internal defect rather
/// than a reachable input condition.
///
/// Named `DenseMatrix` rather than `Matrix` deliberately: this is Eigen's
/// `MatrixXd` for one solver, not a port of `DATASTRUCTURES/Matrix.h`, and the
/// coverage ledger maps candidate Rust types to headers by name.
#[derive(Clone, Debug, PartialEq)]
pub struct DenseMatrix {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}

impl DenseMatrix {
    /// A `rows`-by-`cols` matrix of zeros.
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![0.0; rows.saturating_mul(cols)],
        }
    }

    /// Number of rows.
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Number of columns.
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Element `(row, col)`.
    pub fn at(&self, row: usize, col: usize) -> f64 {
        self.data[row * self.cols + col]
    }

    /// Overwrite element `(row, col)`.
    pub fn set(&mut self, row: usize, col: usize, value: f64) {
        self.data[row * self.cols + col] = value;
    }
}

/// Euclidean norm as Eigen's `stableNorm()` computes it.
///
/// A single scaling pass by the largest magnitude, then the sum of squares of
/// the scaled entries: `scale * sqrt(sum((v / scale)^2))`. A one-element vector
/// short-circuits to its magnitude and an all-zero vector to zero, both as in
/// Eigen. This is not the same expression as [`blue_norm`] and the two are used
/// exactly where Eigen uses each.
pub fn stable_norm(v: &[f64]) -> f64 {
    match v.len() {
        0 => 0.0,
        1 => v[0].abs(),
        _ => {
            let mut scale = 0.0f64;
            for &x in v {
                let ax = x.abs();
                if ax > scale {
                    scale = ax;
                }
            }
            if scale == 0.0 {
                return 0.0;
            }
            let inv = 1.0 / scale;
            let mut ssq = 0.0f64;
            for &x in v {
                let t = x * inv;
                ssq += t * t;
            }
            scale * ssq.sqrt()
        }
    }
}

// Blue's algorithm constants for IEEE binary64, derived exactly as Eigen's
// `blueNorm_impl` derives them from radix 2, 53 digits and exponents -1021 and
// 1024, including the C++ truncating integer divisions.
const BLUE_B1: f64 = 1.4916681462400413e-154; // 2^-511, lower edge of the medium range
const BLUE_B2: f64 = 1.997919072202235e146; // 2^486, upper edge before division by n
const BLUE_S1M: f64 = 6.703903964971299e153; // 2^511, scaling for the small bin
const BLUE_S2M: f64 = 1.1113793747425387e-162; // 2^-538, scaling for the large bin

/// Euclidean norm as Eigen's `blueNorm()` computes it.
///
/// Blue's algorithm: magnitudes are binned into small, medium and large ranges
/// with separate scalings so that neither the squares of the large entries
/// overflow nor those of the small entries underflow. For data whose
/// magnitudes all fall in the medium range - which is every case these fitters
/// meet - the result is `sqrt(sum(v^2))` accumulated in order. A NaN anywhere
/// in the medium bin propagates, as in Eigen.
pub fn blue_norm(v: &[f64]) -> f64 {
    let relerr = f64::EPSILON.sqrt();
    let ab2 = if v.is_empty() {
        BLUE_B2
    } else {
        BLUE_B2 / v.len() as f64
    };
    let mut asml = 0.0f64;
    let mut amed = 0.0f64;
    let mut abig = 0.0f64;
    for &x in v {
        let ax = x.abs();
        if ax > ab2 {
            let t = ax * BLUE_S2M;
            abig += t * t;
        } else if ax < BLUE_B1 {
            let t = ax * BLUE_S1M;
            asml += t * t;
        } else {
            amed += ax * ax;
        }
    }
    if amed.is_nan() {
        return amed;
    }
    if abig > 0.0 {
        abig = abig.sqrt();
        if !abig.is_finite() {
            return abig;
        }
        if amed > 0.0 {
            abig /= BLUE_S2M;
            amed = amed.sqrt();
        } else {
            return abig / BLUE_S2M;
        }
    } else if asml > 0.0 {
        if amed > 0.0 {
            abig = amed.sqrt();
            amed = asml.sqrt() / BLUE_S1M;
        } else {
            return asml.sqrt() / BLUE_S1M;
        }
    } else {
        return amed.sqrt();
    }
    let small = abig.min(amed);
    let large = abig.max(amed);
    if small <= large * relerr {
        large
    } else {
        let ratio = small / large;
        large * (1.0 + ratio * ratio).sqrt()
    }
}

/// Column-pivoted Householder QR, reproducing `Eigen::ColPivHouseholderQR`.
struct ColPivQr {
    m: usize,
    n: usize,
    /// `R` in the upper triangle, the essential Householder vectors below it.
    qr: DenseMatrix,
    /// Householder scaling factors, one per pivot.
    tau: Vec<f64>,
    /// `ind[j]` is the original column now in position `j` (MINPACK's `ipvt`).
    ind: Vec<usize>,
    nonzero_pivots: usize,
    maxpivot: f64,
}

/// The reflector `H = I - tau v v^T` with `v = [1, essential]` that maps
/// `column` to `[beta, 0, ...]`, exactly as `MatrixBase::makeHouseholder`.
fn make_householder(column: &[f64]) -> (f64, f64, Vec<f64>) {
    let c0 = column[0];
    let tail = &column[1..];
    let mut tail_sq = 0.0f64;
    for &x in tail {
        tail_sq += x * x;
    }
    if tail_sq <= f64::MIN_POSITIVE {
        return (0.0, c0, vec![0.0; tail.len()]);
    }
    let mut beta = (c0 * c0 + tail_sq).sqrt();
    if c0 >= 0.0 {
        beta = -beta;
    }
    let denominator = c0 - beta;
    let essential = tail.iter().map(|x| x / denominator).collect();
    ((beta - c0) / beta, beta, essential)
}

/// Apply `H` to the sub-block of `mat` whose first row is `row0` and whose
/// first column is `col0`, as `applyHouseholderOnTheLeft`.
fn apply_householder_left(
    mat: &mut DenseMatrix,
    row0: usize,
    col0: usize,
    essential: &[f64],
    tau: f64,
) {
    let rows = mat.rows() - row0;
    let cols = mat.cols().saturating_sub(col0);
    if cols == 0 {
        return;
    }
    if rows == 1 {
        for col in col0..mat.cols() {
            let value = mat.at(row0, col) * (1.0 - tau);
            mat.set(row0, col, value);
        }
        return;
    }
    if tau == 0.0 {
        return;
    }
    let mut tmp = vec![0.0f64; cols];
    for (offset, slot) in tmp.iter_mut().enumerate() {
        let col = col0 + offset;
        let mut sum = 0.0f64;
        for (below, &e) in essential.iter().enumerate() {
            sum += e * mat.at(row0 + 1 + below, col);
        }
        *slot = sum + mat.at(row0, col);
    }
    for (offset, &t) in tmp.iter().enumerate() {
        let col = col0 + offset;
        let value = mat.at(row0, col) - tau * t;
        mat.set(row0, col, value);
    }
    for (below, &e) in essential.iter().enumerate() {
        let row = row0 + 1 + below;
        for (offset, &t) in tmp.iter().enumerate() {
            let col = col0 + offset;
            let value = mat.at(row, col) - tau * e * t;
            mat.set(row, col, value);
        }
    }
}

/// Apply `H` to the tail of a vector starting at `row0`.
fn apply_householder_vector(w: &mut [f64], row0: usize, essential: &[f64], tau: f64) {
    let rows = w.len() - row0;
    if rows == 1 {
        w[row0] *= 1.0 - tau;
        return;
    }
    if tau == 0.0 {
        return;
    }
    let mut sum = 0.0f64;
    for (below, &e) in essential.iter().enumerate() {
        sum += e * w[row0 + 1 + below];
    }
    sum += w[row0];
    w[row0] -= tau * sum;
    for (below, &e) in essential.iter().enumerate() {
        w[row0 + 1 + below] -= tau * e * sum;
    }
}

fn plain_norm(v: &[f64]) -> f64 {
    let mut sum = 0.0f64;
    for &x in v {
        sum += x * x;
    }
    sum.sqrt()
}

impl ColPivQr {
    /// Factorize `a` in place, choosing at each step the remaining column of
    /// largest updated norm, with LAPACK's norm-downdating rule and its
    /// recomputation threshold of `sqrt(f64::EPSILON)`.
    fn new(a: &DenseMatrix) -> Self {
        let m = a.rows();
        let n = a.cols();
        let size = m.min(n);
        let mut qr = a.clone();
        let mut tau = vec![0.0f64; size];
        let mut updated = Vec::with_capacity(n);
        for col in 0..n {
            let column: Vec<f64> = (0..m).map(|row| qr.at(row, col)).collect();
            updated.push(plain_norm(&column));
        }
        let mut direct = updated.clone();
        let mut biggest = 0.0f64;
        for &value in &updated {
            if value > biggest {
                biggest = value;
            }
        }
        let helper = {
            let scaled = biggest * f64::EPSILON;
            scaled * scaled / m as f64
        };
        let downdate = f64::EPSILON.sqrt();
        let mut nonzero_pivots = size;
        let mut maxpivot = 0.0f64;
        let mut ind: Vec<usize> = (0..n).collect();
        for k in 0..size {
            let mut best = k;
            for j in (k + 1)..n {
                if updated[j] > updated[best] {
                    best = j;
                }
            }
            let biggest_sq = updated[best] * updated[best];
            if nonzero_pivots == size && biggest_sq < helper * (m - k) as f64 {
                nonzero_pivots = k;
            }
            if k != best {
                for row in 0..m {
                    let left = qr.at(row, k);
                    let right = qr.at(row, best);
                    qr.set(row, k, right);
                    qr.set(row, best, left);
                }
                updated.swap(k, best);
                direct.swap(k, best);
                ind.swap(k, best);
            }
            let column: Vec<f64> = (k..m).map(|row| qr.at(row, k)).collect();
            let (t, beta, essential) = make_householder(&column);
            tau[k] = t;
            qr.set(k, k, beta);
            for (below, &e) in essential.iter().enumerate() {
                qr.set(k + 1 + below, k, e);
            }
            if beta.abs() > maxpivot {
                maxpivot = beta.abs();
            }
            apply_householder_left(&mut qr, k, k + 1, &essential, t);
            for j in (k + 1)..n {
                if updated[j] == 0.0 {
                    continue;
                }
                let mut ratio = qr.at(k, j).abs() / updated[j];
                ratio = (1.0 + ratio) * (1.0 - ratio);
                if ratio < 0.0 {
                    ratio = 0.0;
                }
                let scaled = updated[j] / direct[j];
                let guard = ratio * scaled * scaled;
                if guard <= downdate {
                    let rest: Vec<f64> = ((k + 1)..m).map(|row| qr.at(row, j)).collect();
                    direct[j] = plain_norm(&rest);
                    updated[j] = direct[j];
                } else {
                    updated[j] *= ratio.sqrt();
                }
            }
        }
        Self {
            m,
            n,
            qr,
            tau,
            ind,
            nonzero_pivots,
            maxpivot,
        }
    }

    /// Numerical rank at Eigen's default threshold of
    /// `|maxpivot| * f64::EPSILON * min(rows, cols)`.
    fn rank(&self) -> usize {
        let threshold = self.maxpivot.abs() * f64::EPSILON * self.m.min(self.n) as f64;
        (0..self.nonzero_pivots)
            .filter(|&i| self.qr.at(i, i).abs() > threshold)
            .count()
    }

    /// `Q^T w`, by applying the reflectors in factorization order.
    fn transpose_apply(&self, w: &[f64]) -> Vec<f64> {
        let mut out = w.to_vec();
        for k in 0..self.tau.len() {
            let essential: Vec<f64> = ((k + 1)..self.m).map(|row| self.qr.at(row, k)).collect();
            apply_householder_vector(&mut out, k, &essential, self.tau[k]);
        }
        out
    }
}

/// The Givens rotation of `Eigen::JacobiRotation::makeGivens` for reals,
/// returning `(c, s)`.
fn make_givens(p: f64, q: f64) -> (f64, f64) {
    if q == 0.0 {
        return (if p < 0.0 { -1.0 } else { 1.0 }, 0.0);
    }
    if p == 0.0 {
        return (0.0, if q < 0.0 { 1.0 } else { -1.0 });
    }
    if p.abs() > q.abs() {
        let t = q / p;
        let mut u = (1.0 + t * t).sqrt();
        if p < 0.0 {
            u = -u;
        }
        let c = 1.0 / u;
        (c, -t * c)
    } else {
        let t = p / q;
        let mut u = (1.0 + t * t).sqrt();
        if q < 0.0 {
            u = -u;
        }
        let s = -1.0 / u;
        (-t * s, s)
    }
}

/// MINPACK `qrsolv`: solve `(R^T R + D^2) x = R^T Q^T b` by eliminating the
/// diagonal `D` with Givens rotations.
///
/// `s` is the `n`-by-`n` leading block of the QR factor; the routine mirrors
/// its upper triangle into the strict lower triangle, transforms the lower
/// triangle and diagonal, then restores the original diagonal and reports the
/// transformed one in `sdiag`, exactly as Eigen does.
fn qrsolv(
    s: &mut DenseMatrix,
    ind: &[usize],
    diag: &[f64],
    qtb: &[f64],
    n: usize,
) -> (Vec<f64>, Vec<f64>) {
    let saved: Vec<f64> = (0..n).map(|j| s.at(j, j)).collect();
    let mut wa = qtb[..n].to_vec();
    for i in 0..n {
        for j in 0..i {
            let mirrored = s.at(j, i);
            s.set(i, j, mirrored);
        }
    }
    // Eigen leaves this workspace uninitialized; every element is written
    // before the `nsing` scan reads it unless a zero scaling factor breaks the
    // loop early, which the driver's positive `diag` makes unreachable.
    let mut sdiag = vec![0.0f64; n];
    for j in 0..n {
        let l = ind[j];
        if diag[l] == 0.0 {
            break;
        }
        for slot in sdiag.iter_mut().skip(j) {
            *slot = 0.0;
        }
        sdiag[j] = diag[l];
        let mut qtbpj = 0.0f64;
        for k in j..n {
            let (c, sn) = make_givens(-s.at(k, k), sdiag[k]);
            let updated = c * s.at(k, k) + sn * sdiag[k];
            s.set(k, k, updated);
            let temp = c * wa[k] + sn * qtbpj;
            qtbpj = -sn * wa[k] + c * qtbpj;
            wa[k] = temp;
            for i in (k + 1)..n {
                let temp = c * s.at(i, k) + sn * sdiag[i];
                sdiag[i] = -sn * s.at(i, k) + c * sdiag[i];
                s.set(i, k, temp);
            }
        }
    }
    let mut nsing = 0usize;
    while nsing < n && sdiag[nsing] != 0.0 {
        nsing += 1;
    }
    for slot in wa.iter_mut().skip(nsing) {
        *slot = 0.0;
    }
    for k in (0..nsing).rev() {
        let mut sum = 0.0f64;
        for i in (k + 1)..nsing {
            sum += s.at(i, k) * wa[i];
        }
        wa[k] = (wa[k] - sum) / s.at(k, k);
    }
    let transformed: Vec<f64> = (0..n).map(|j| s.at(j, j)).collect();
    for (j, &value) in saved.iter().enumerate() {
        s.set(j, j, value);
    }
    let mut x = vec![0.0f64; n];
    for j in 0..n {
        x[ind[j]] = wa[j];
    }
    (x, transformed)
}

/// MINPACK `lmpar`, in Eigen's `lmpar2` form: find the Levenberg parameter for
/// which the scaled step length matches the trust-region radius `delta`,
/// capped at ten bisection-style iterations.
fn lmpar(qr: &ColPivQr, diag: &[f64], qtb: &[f64], delta: f64, par_in: f64) -> (Vec<f64>, f64) {
    let n = qr.n;
    let rank = qr.rank();
    let mut wa1 = qtb[..n].to_vec();
    for slot in wa1.iter_mut().skip(rank) {
        *slot = 0.0;
    }
    for k in (0..rank).rev() {
        let mut sum = 0.0f64;
        for i in (k + 1)..rank {
            sum += qr.qr.at(k, i) * wa1[i];
        }
        wa1[k] = (wa1[k] - sum) / qr.qr.at(k, k);
    }
    let mut x = vec![0.0f64; n];
    for j in 0..n {
        x[qr.ind[j]] = wa1[j];
    }
    let mut iter = 0usize;
    let mut wa2: Vec<f64> = (0..n).map(|i| diag[i] * x[i]).collect();
    let mut dxnorm = blue_norm(&wa2);
    let mut fp = dxnorm - delta;
    if fp <= 0.1 * delta {
        return (x, 0.0);
    }
    let mut parl = 0.0f64;
    if rank == n {
        let mut work: Vec<f64> = (0..n)
            .map(|j| diag[qr.ind[j]] * (wa2[qr.ind[j]] / dxnorm))
            .collect();
        for j in 0..n {
            let mut sum = 0.0f64;
            for i in 0..j {
                sum += qr.qr.at(i, j) * work[i];
            }
            work[j] = (work[j] - sum) / qr.qr.at(j, j);
        }
        let temp = blue_norm(&work);
        parl = fp / delta / temp / temp;
    }
    let mut upper = vec![0.0f64; n];
    for (j, slot) in upper.iter_mut().enumerate() {
        let mut sum = 0.0f64;
        for i in 0..=j {
            sum += qr.qr.at(i, j) * qtb[i];
        }
        *slot = sum / diag[qr.ind[j]];
    }
    let gnorm = stable_norm(&upper);
    let mut paru = gnorm / delta;
    if paru == 0.0 {
        paru = f64::MIN_POSITIVE / delta.min(0.1);
    }
    let mut par = par_in.max(parl).min(paru);
    if par == 0.0 {
        par = gnorm / dxnorm;
    }
    let mut s = DenseMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            s.set(i, j, qr.qr.at(i, j));
        }
    }
    loop {
        iter += 1;
        if par == 0.0 {
            par = f64::MIN_POSITIVE.max(0.001 * paru);
        }
        let scaled: Vec<f64> = diag.iter().map(|d| par.sqrt() * d).collect();
        let (next, sdiag) = qrsolv(&mut s, &qr.ind, &scaled, qtb, n);
        x = next;
        wa2 = (0..n).map(|i| diag[i] * x[i]).collect();
        dxnorm = blue_norm(&wa2);
        let previous = fp;
        fp = dxnorm - delta;
        if fp.abs() <= 0.1 * delta
            || (parl == 0.0 && fp <= previous && previous < 0.0)
            || iter == 10
        {
            break;
        }
        let mut work: Vec<f64> = (0..n)
            .map(|j| diag[qr.ind[j]] * (wa2[qr.ind[j]] / dxnorm))
            .collect();
        for j in 0..n {
            work[j] /= sdiag[j];
            let temp = work[j];
            for i in (j + 1)..n {
                work[i] -= s.at(i, j) * temp;
            }
        }
        let temp = blue_norm(&work);
        let parc = fp / delta / temp / temp;
        if fp > 0.0 {
            parl = parl.max(par);
        }
        if fp < 0.0 {
            paru = paru.min(par);
        }
        par = parl.max(par + parc);
    }
    if iter == 0 {
        par = 0.0;
    }
    (x, par)
}

/// Minimize `sum(residuals^2)` over `x`, reproducing
/// `Eigen::LevenbergMarquardt::minimize`.
///
/// `values` is the number of residuals. `residuals` fills its second argument
/// from the parameter vector. `jacobian` fills the `values`-by-`x.len()`
/// matrix and returns how many residual evaluations it consumed - zero for an
/// analytic Jacobian, which makes the call a Jacobian evaluation instead, and
/// `n + 1` for the forward-difference approximation Eigen's `NumericalDiff`
/// performs. That count is what `max_fev` is compared against, as in Eigen.
///
/// `x` holds the initial guess on entry and the fitted parameters on return,
/// including when the returned status reports a failure - the C++ fitters read
/// the vector back the same way and only then decide whether to throw.
///
/// The iteration is serial. The source is serial here too: no `#pragma omp`
/// appears in Eigen's non-linear optimization module or in the four fitters.
pub fn minimize<R, J>(
    x: &mut [f64],
    values: usize,
    mut residuals: R,
    mut jacobian: J,
    parameters: &LmParameters,
) -> LmStatus
where
    R: FnMut(&[f64], &mut [f64]),
    J: FnMut(&[f64], &mut DenseMatrix) -> usize,
{
    let n = x.len();
    let m = values;
    if n == 0
        || m < n
        || parameters.ftol < 0.0
        || parameters.xtol < 0.0
        || parameters.gtol < 0.0
        || parameters.max_fev == 0
        || parameters.factor <= 0.0
    {
        return LmStatus::ImproperInputParameters;
    }

    let mut nfev = 1usize;
    let mut fvec = vec![0.0f64; m];
    residuals(x, &mut fvec);
    let mut fnorm = stable_norm(&fvec);
    let mut par = 0.0f64;
    let mut iter = 1usize;
    let mut diag = vec![0.0f64; n];
    let mut delta = 0.0f64;
    let mut xnorm = 0.0f64;
    let mut fjac = DenseMatrix::zeros(m, n);
    let mut trial = vec![0.0f64; m];

    loop {
        let consumed = jacobian(x, &mut fjac);
        if consumed > 0 {
            nfev = nfev.saturating_add(consumed);
        }
        let column_norms: Vec<f64> = (0..n)
            .map(|j| {
                let column: Vec<f64> = (0..m).map(|i| fjac.at(i, j)).collect();
                blue_norm(&column)
            })
            .collect();
        let qr = ColPivQr::new(&fjac);
        if iter == 1 {
            for (j, slot) in diag.iter_mut().enumerate() {
                *slot = if column_norms[j] == 0.0 {
                    1.0
                } else {
                    column_norms[j]
                };
            }
            let scaled: Vec<f64> = (0..n).map(|j| diag[j] * x[j]).collect();
            xnorm = stable_norm(&scaled);
            delta = parameters.factor * xnorm;
            if delta == 0.0 {
                delta = parameters.factor;
            }
        }
        let projected = qr.transpose_apply(&fvec);
        let qtf = &projected[..n];
        let mut gnorm = 0.0f64;
        if fnorm != 0.0 {
            for j in 0..n {
                if column_norms[qr.ind[j]] == 0.0 {
                    continue;
                }
                let mut sum = 0.0f64;
                for i in 0..=j {
                    sum += qr.qr.at(i, j) * (qtf[i] / fnorm);
                }
                let candidate = (sum / column_norms[qr.ind[j]]).abs();
                if candidate > gnorm {
                    gnorm = candidate;
                }
            }
        }
        if gnorm <= parameters.gtol {
            return LmStatus::CosinusTooSmall;
        }
        for (j, slot) in diag.iter_mut().enumerate() {
            if column_norms[j] > *slot {
                *slot = column_norms[j];
            }
        }

        loop {
            let (step, next_par) = lmpar(&qr, &diag, qtf, delta, par);
            par = next_par;
            let step: Vec<f64> = step.iter().map(|v| -v).collect();
            let candidate: Vec<f64> = (0..n).map(|j| x[j] + step[j]).collect();
            let scaled_step: Vec<f64> = (0..n).map(|j| diag[j] * step[j]).collect();
            let pnorm = stable_norm(&scaled_step);
            if iter == 1 {
                delta = delta.min(pnorm);
            }
            residuals(&candidate, &mut trial);
            nfev = nfev.saturating_add(1);
            let fnorm1 = stable_norm(&trial);
            let mut actred = -1.0f64;
            if 0.1 * fnorm1 < fnorm {
                let ratio = fnorm1 / fnorm;
                actred = 1.0 - ratio * ratio;
            }
            let permuted: Vec<f64> = (0..n).map(|j| step[qr.ind[j]]).collect();
            let mut projected_step = vec![0.0f64; n];
            for (i, slot) in projected_step.iter_mut().enumerate() {
                let mut sum = 0.0f64;
                for j in i..n {
                    sum += qr.qr.at(i, j) * permuted[j];
                }
                *slot = sum;
            }
            let temp1 = {
                let t = stable_norm(&projected_step) / fnorm;
                t * t
            };
            let temp2 = {
                let t = par.sqrt() * pnorm / fnorm;
                t * t
            };
            let prered = temp1 + temp2 / 0.5;
            let dirder = -(temp1 + temp2);
            let mut ratio = 0.0f64;
            if prered != 0.0 {
                ratio = actred / prered;
            }
            if ratio <= 0.25 {
                let mut temp = 0.5f64;
                if actred < 0.0 {
                    temp = 0.5 * dirder / (dirder + 0.5 * actred);
                }
                if 0.1 * fnorm1 >= fnorm || temp < 0.1 {
                    temp = 0.1;
                }
                delta = temp * delta.min(pnorm / 0.1);
                par /= temp;
            } else if !(par != 0.0 && ratio < 0.75) {
                delta = pnorm / 0.5;
                par *= 0.5;
            }
            if ratio >= 1e-4 {
                x.copy_from_slice(&candidate);
                let scaled: Vec<f64> = (0..n).map(|j| diag[j] * x[j]).collect();
                xnorm = stable_norm(&scaled);
                fvec.copy_from_slice(&trial);
                fnorm = fnorm1;
                iter += 1;
            }
            let reduction_small =
                actred.abs() <= parameters.ftol && prered <= parameters.ftol && 0.5 * ratio <= 1.0;
            if reduction_small && delta <= parameters.xtol * xnorm {
                return LmStatus::RelativeErrorAndReductionTooSmall;
            }
            if reduction_small {
                return LmStatus::RelativeReductionTooSmall;
            }
            if delta <= parameters.xtol * xnorm {
                return LmStatus::RelativeErrorTooSmall;
            }
            if nfev >= parameters.max_fev {
                return LmStatus::TooManyFunctionEvaluation;
            }
            if actred.abs() <= f64::EPSILON && prered <= f64::EPSILON && 0.5 * ratio <= 1.0 {
                return LmStatus::FtolTooSmall;
            }
            if delta <= f64::EPSILON * xnorm {
                return LmStatus::XtolTooSmall;
            }
            if gnorm <= f64::EPSILON {
                return LmStatus::GtolTooSmall;
            }
            if ratio >= 1e-4 {
                break;
            }
        }
    }
}

/// Forward-difference Jacobian, reproducing `Eigen::NumericalDiff` in its
/// default `Forward` mode with `epsfcn = 0`.
///
/// The step is `sqrt(f64::EPSILON) * |x[j]|`, or `sqrt(f64::EPSILON)` when that
/// is zero. Returns `n + 1`, the number of residual evaluations it spent, which
/// the driver adds to its `max_fev` budget exactly as Eigen does.
pub fn numerical_jacobian<R>(
    x: &[f64],
    values: usize,
    jac: &mut DenseMatrix,
    mut residuals: R,
) -> usize
where
    R: FnMut(&[f64], &mut [f64]),
{
    let n = x.len();
    let eps = f64::EPSILON.sqrt();
    let mut base = vec![0.0f64; values];
    residuals(x, &mut base);
    let mut shifted = vec![0.0f64; values];
    let mut probe = x.to_vec();
    for j in 0..n {
        let mut h = eps * x[j].abs();
        if h == 0.0 {
            h = eps;
        }
        probe[j] = x[j] + h;
        residuals(&probe, &mut shifted);
        probe[j] = x[j];
        for i in 0..values {
            jac.set(i, j, (shifted[i] - base[i]) / h);
        }
    }
    n + 1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blue_norm_constants_are_the_powers_of_two_eigen_derives() {
        assert_eq!(BLUE_B1, 2.0f64.powi(-511));
        assert_eq!(BLUE_B2, 2.0f64.powi(486));
        assert_eq!(BLUE_S1M, 2.0f64.powi(511));
        assert_eq!(BLUE_S2M, 2.0f64.powi(-538));
    }

    #[test]
    fn stable_norm_matches_the_plain_norm_for_ordinary_data() {
        let v = [3.0, 4.0];
        assert!((stable_norm(&v) - 5.0).abs() < 1e-15);
        assert!((blue_norm(&v) - 5.0).abs() < 1e-15);
        assert_eq!(stable_norm(&[]), 0.0);
        assert_eq!(stable_norm(&[-7.0]), 7.0);
        assert_eq!(stable_norm(&[0.0, 0.0]), 0.0);
    }

    #[test]
    fn stable_norm_survives_magnitudes_that_overflow_a_naive_sum() {
        let v = [1e200, 1e200];
        let expected = 1e200 * 2.0f64.sqrt();
        assert!((stable_norm(&v) / expected - 1.0).abs() < 1e-15);
        assert!((blue_norm(&v) / expected - 1.0).abs() < 1e-15);
        let tiny = [1e-200, 1e-200];
        let expected = 1e-200 * 2.0f64.sqrt();
        assert!((blue_norm(&tiny) / expected - 1.0).abs() < 1e-15);
    }

    #[test]
    fn fewer_residuals_than_parameters_is_improper() {
        let mut x = [1.0, 1.0, 1.0];
        let status = minimize(
            &mut x,
            2,
            |_, f: &mut [f64]| f.fill(0.0),
            |_, _: &mut DenseMatrix| 0,
            &LmParameters::default(),
        );
        assert_eq!(status, LmStatus::ImproperInputParameters);
        assert_eq!(status.code(), 0);
    }

    #[test]
    fn a_linear_least_squares_problem_is_solved_exactly() {
        // Fit y = a + b x through (0,1), (1,3), (2,5): the exact solution is
        // a = 1, b = 2 and the residuals vanish, so this is independent of any
        // transcribed literal.
        let points = [(0.0, 1.0), (1.0, 3.0), (2.0, 5.0)];
        let mut x = [0.0, 0.0];
        let status = minimize(
            &mut x,
            3,
            |p: &[f64], f: &mut [f64]| {
                for (i, &(px, py)) in points.iter().enumerate() {
                    f[i] = p[0] + p[1] * px - py;
                }
            },
            |_: &[f64], jac: &mut DenseMatrix| {
                for (i, &(px, _)) in points.iter().enumerate() {
                    jac.set(i, 0, 1.0);
                    jac.set(i, 1, px);
                }
                0
            },
            &LmParameters::default(),
        );
        assert_ne!(status, LmStatus::ImproperInputParameters);
        assert!((x[0] - 1.0).abs() < 1e-10, "{x:?}");
        assert!((x[1] - 2.0).abs() < 1e-10, "{x:?}");
    }

    #[test]
    fn point_preflight_rejects_oversized_input() {
        assert!(preflight_points(10, 3).is_ok());
        // 1e6 points times three parameters is 24 MiB, inside MAX_BYTES.
        assert!(preflight_points(MAX_POINTS, 3).is_ok());
        assert!(preflight_points(MAX_POINTS + 1, 3).is_err());
        // A wider model reaches the byte ceiling before the point ceiling.
        assert!(preflight_points(MAX_POINTS, 9).is_err());
        assert!(preflight_points(usize::MAX, 3).is_err());
    }
}
