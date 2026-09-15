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
//! with a Givens-rotation `qrsolv`, and MINPACK's termination tests. That
//! algorithm is reproduced here, step for step and in the same arithmetic
//! order, so that the fitters converge to the published parameters rather than
//! to some other point of a non-convex surface. `TraceFitter::optimize_`
//! (`FEATUREFINDER/TraceFitter.cpp`) calls the same Eigen class, with
//! `maxfev` set from `max_iteration`, and uses this module too.
//!
//! **Why not the `levenberg-marquardt` crate.** Package B3-LM measured
//! `levenberg-marquardt =0.14.0` behind this signature, with Eigen's `maxfev`
//! emulated exactly, against the C2 class-level oracle and did not adopt it.
//! Its evaluation accounting matched Eigen at 29,003 of 29,004 trace-fit
//! budgets and this transcription at all 8,000 distribution-fitter budgets, but
//! its fitted parameters were further from the executed C++ than this
//! transcription's at 3,811 of 29,004 trace-fit budgets (up to `2.2e-9`
//! relative where the transcription is at `6.4e-10`), and on a degenerate
//! flat-trace fit it took a step Eigen does not take, ending with a different
//! status after a different number of evaluations.
//! This transcription reproduces Eigen's status, `nfev` and `njev` at all
//! 29,004 budgets. The measurements are in `docs/DISTRIBUTION_FITTERS_SUPPORT.md`
//! and `tests/lm_budget_differential.rs`.
//!
//! **Arithmetic order.** Step for step is not enough for the same bits: Eigen
//! accumulates every `squaredNorm()`, `dot()`, matrix-vector product and
//! triangular solve in SIMD lanes. Package B3b-LM-FIDELITY traced both paths on
//! 141 trace fits and found the first divergence in the residual norm and the
//! QR column norms (lane order), then in the Householder projections (Eigen's
//! row-major matrix-vector kernel), the Gauss-Newton back substitution (Eigen
//! subtracts columns), and the predicted-reduction norm (Eigen's `wa3` has `m`
//! rows, not `n`). The kernels here reproduce Eigen 5.0.1 with two-lane
//! packets, which is what both OpenMS reference builds use. With them every
//! evaluation path is bit-identical to the Linux x86_64 Release build; on arm64
//! Eigen additionally fuses the lanes with FMA, which is not modelled (see
//! [`minimize`](crate::math::fitters::levenberg_marquardt::minimize)).
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
    ///
    /// Counted as Eigen counts `nfev`: 1 for the start, 1 per trial step, plus
    /// whatever a Jacobian evaluation reports consuming (0 for an analytic
    /// Jacobian, `n + 1` for [`numerical_jacobian`]). The count is compared
    /// once after each trial step, after the `ftol` and `xtol` tests, so a fit
    /// that converges on the step that reaches the budget still reports
    /// convergence, and a budget of 1 still takes one trial step.
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

// ---------------------------------------------------------------------------
// Eigen 5.0.1 reduction kernels
//
// Eigen does not accumulate its sums and inner products left to right. Every
// `squaredNorm()`, `dot()`, row-major matrix-vector product and triangular
// solve on this path runs through a SIMD kernel whose summation order is fixed
// by the packet width. Both reference builds of OpenMS use two-lane `f64`
// packets: NEON `Packet2d` on arm64, and SSE `Packet2d` on x86_64, where
// `cmake/compiler_flags.cmake` passes `-mssse3` and deliberately no AVX
// ("AVX's 256-bit reductions change Eigen's floating-point evaluation order").
// The kernels below reproduce those two-lane orders. B3b-LM-FIDELITY measured
// this against the executed Eigen on both platforms; see
// `docs/DISTRIBUTION_FITTERS_SUPPORT.md` section 1 and
// `tests/lm_eigen_path_differential.rs`.
// ---------------------------------------------------------------------------

/// One SIMD lane of Eigen's `pmadd(a, b, c) = a * b + c`, the accumulation
/// step of the inner-product and row-major matrix-vector kernels.
///
/// Not fused. That is what Eigen does wherever `EIGEN_VECTORIZE_FMA` is
/// undefined, which includes the x86_64 OpenMS builds (`-mssse3`). On arm64
/// Eigen defines it from `__ARM_FEATURE_FMA` and `pmadd` becomes
/// `vfmaq_f64`, a fused multiply-add, whatever `-ffp-contract` says; this
/// helper does not model that (see the platform note in
/// `docs/DISTRIBUTION_FITTERS_SUPPORT.md` section 1). The scalar tails of the
/// same kernels are unfused on both platforms, because
/// `EIGEN_SCALAR_MADD_USE_FMA` is fixed before the FMA detection runs.
#[inline]
fn lane_madd(a: f64, b: f64, accumulator: f64) -> f64 {
    a * b + accumulator
}

/// `DenseBase::sum()` of `size >= 1` terms: `redux_impl` with
/// `LinearVectorizedTraversal` (`Redux.h:275-322`), two-lane packets and an
/// aligned start of 0, which it always is here because every summand is an
/// expression without direct access.
///
/// Two packet accumulators take four terms per step; they are added lane-wise,
/// a trailing packet joins lane-wise, the two lanes are added, and any odd term
/// is added last.
fn eigen_sum(size: usize, term: impl Fn(usize) -> f64) -> f64 {
    let aligned = size / 2 * 2;
    if aligned == 0 {
        return term(0);
    }
    let aligned4 = size / 4 * 4;
    let mut lane0 = term(0);
    let mut lane1 = term(1);
    if aligned > 2 {
        let mut lane2 = term(2);
        let mut lane3 = term(3);
        let mut index = 4;
        while index < aligned4 {
            lane0 += term(index);
            lane1 += term(index + 1);
            lane2 += term(index + 2);
            lane3 += term(index + 3);
            index += 4;
        }
        lane0 += lane2;
        lane1 += lane3;
        if aligned > aligned4 {
            lane0 += term(aligned4);
            lane1 += term(aligned4 + 1);
        }
    }
    let mut result = lane0 + lane1;
    for index in aligned..size {
        result += term(index);
    }
    result
}

/// `squaredNorm()` of `size` coefficients (`Dot.h:21-27`): the reduction of
/// the squares. Zero for an empty vector, as Eigen's `sum()`.
fn eigen_squared_norm(size: usize, coeff: impl Fn(usize) -> f64) -> f64 {
    if size == 0 {
        return 0.0;
    }
    eigen_sum(size, |i| {
        let v = coeff(i);
        v * v
    })
}

/// `a.dot(b)` over `size` coefficients: `inner_product_impl`
/// (`InnerProduct.h:117-172`), which Eigen 5 uses instead of `redux`.
///
/// Four packet accumulators of two lanes each; the first four packets seed
/// them, later packets are multiply-added in, and the accumulators are folded
/// as `a2 += a3; a1 += a2; a0 += a1` before the two lanes of `a0` are added.
/// The odd tail is multiply-added one coefficient at a time.
fn eigen_dot(size: usize, a: impl Fn(usize) -> f64, b: impl Fn(usize) -> f64) -> f64 {
    if size == 0 {
        return 0.0;
    }
    if size < 2 {
        return a(0) * b(0);
    }
    let packet_end = size / 2 * 2;
    let quad_end = size / 8 * 8;
    let packets = size / 2;
    let remaining = (packet_end - quad_end) / 2;
    let seed = |i: usize| [a(i) * b(i), a(i + 1) * b(i + 1)];
    let madd = |acc: [f64; 2], i: usize| {
        [
            lane_madd(a(i), b(i), acc[0]),
            lane_madd(a(i + 1), b(i + 1), acc[1]),
        ]
    };
    let add = |x: [f64; 2], y: [f64; 2]| [x[0] + y[0], x[1] + y[1]];
    let mut acc0 = seed(0);
    let mut acc1 = [0.0; 2];
    let mut acc2 = [0.0; 2];
    if packets >= 2 {
        acc1 = seed(2);
    }
    if packets >= 3 {
        acc2 = seed(4);
    }
    if packets >= 4 {
        let mut acc3 = seed(6);
        let mut k = 8;
        while k < quad_end {
            acc0 = madd(acc0, k);
            acc1 = madd(acc1, k + 2);
            acc2 = madd(acc2, k + 4);
            acc3 = madd(acc3, k + 6);
            k += 8;
        }
        if remaining >= 1 {
            acc0 = madd(acc0, quad_end);
        }
        if remaining >= 2 {
            acc1 = madd(acc1, quad_end + 2);
        }
        if remaining == 3 {
            acc2 = madd(acc2, quad_end + 4);
        }
        acc2 = add(acc2, acc3);
    }
    if packets >= 3 {
        acc1 = add(acc1, acc2);
    }
    if packets >= 2 {
        acc0 = add(acc0, acc1);
    }
    let mut result = acc0[0] + acc0[1];
    for k in packet_end..size {
        result += a(k) * b(k);
    }
    result
}

/// One result coefficient of the row-major `general_matrix_vector_product`
/// (`GeneralMatrixVector.h:298-462`): two lanes multiply-added over the
/// packet-aligned prefix starting from zero, the lanes added, then a scalar
/// tail. The `1 * cc` scaling and the addition to a zeroed result are exact.
fn eigen_gemv_row(size: usize, a: impl Fn(usize) -> f64, b: impl Fn(usize) -> f64) -> f64 {
    let full = size / 2 * 2;
    let (mut lane0, mut lane1) = (0.0f64, 0.0f64);
    let mut j = 0;
    while j < full {
        lane0 = lane_madd(a(j), b(j), lane0);
        lane1 = lane_madd(a(j + 1), b(j + 1), lane1);
        j += 2;
    }
    let mut result = lane0 + lane1;
    for j in full..size {
        result += a(j) * b(j);
    }
    0.0 + result
}

/// Euclidean norm as Eigen's `stableNorm()` computes it.
///
/// A single scaling pass by the largest magnitude, then the sum of squares of
/// the scaled entries: `scale * sqrt(sum((v / scale)^2))`, in blocks of 4096
/// coefficients as `stable_norm_impl_inner_step` walks them, with the sum
/// accumulated in Eigen's two-lane reduction order. A one-element vector
/// short-circuits to its magnitude and an all-zero vector to zero, both as in
/// Eigen, and a NaN in the first coefficient of a block becomes the scale and
/// makes the norm NaN, as `maxCoeff` does.
///
/// Eigen's path through this algorithm uses three different norm expressions,
/// not two, and each is reproduced where Eigen uses it: `stableNorm` for the
/// residual, step and scaled-`x` norms and for `lmpar`'s `gnorm`; [`blue_norm`]
/// for the Jacobian column norms the driver turns into `diag` and for the two
/// scaled-step norms inside `lmpar`; and `sqrt(squaredNorm())` of
/// `MatrixBase::norm()`, the private `plain_norm`, for the pivot column norms
/// inside the column-pivoted QR.
pub fn stable_norm(v: &[f64]) -> f64 {
    stable_norm_by(v.len(), |i| v[i])
}

/// [`stable_norm`] over `len` coefficients read through `coeff`.
fn stable_norm_by(len: usize, coeff: impl Fn(usize) -> f64) -> f64 {
    const BLOCK: usize = 4096;
    if len == 0 {
        return 0.0;
    }
    if len == 1 {
        return coeff(0).abs();
    }
    let mut scale = 0.0f64;
    let mut inv_scale = 1.0f64;
    let mut ssq = 0.0f64;
    let mut start = 0;
    while start < len {
        let size = (len - start).min(BLOCK);
        // `bl.cwiseAbs().maxCoeff()`: the first coefficient seeds the maximum
        // and only a strictly larger one replaces it, so NaN wins only there.
        let mut max_coeff = coeff(start).abs();
        for i in 1..size {
            let candidate = coeff(start + i).abs();
            if candidate > max_coeff {
                max_coeff = candidate;
            }
        }
        // `stable_norm_kernel` does not take the plain reciprocal: it guards
        // the two ends of the range first. A subnormal largest coefficient
        // makes `1 / maxCoeff` overflow, and an infinite one makes it zero; in
        // either case Eigen substitutes a usable pair. Unreachable from the
        // fitters, whose inputs are validated finite, but transcribed rather
        // than assumed away.
        if max_coeff > scale {
            let ratio = scale / max_coeff;
            ssq *= ratio * ratio;
            let reciprocal = 1.0 / max_coeff;
            if reciprocal > f64::MAX {
                inv_scale = f64::MAX;
                scale = 1.0 / inv_scale;
            } else if max_coeff > f64::MAX {
                inv_scale = 1.0;
                scale = max_coeff;
            } else {
                scale = max_coeff;
                inv_scale = reciprocal;
            }
        } else if max_coeff.is_nan() {
            scale = max_coeff;
        }
        if scale > 0.0 {
            ssq += eigen_squared_norm(size, |i| coeff(start + i) * inv_scale);
        }
        start += size;
    }
    scale * ssq.sqrt()
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
///
/// This is *not* the norm the column-pivoted QR uses for pivoting; see
/// [`stable_norm`] for which of the three goes where.
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

/// `(std::max)(a, b)` as Eigen's solver calls it: `a < b ? b : a`. Unlike
/// `f64::max`, a NaN in `a` is returned and a NaN in `b` is ignored.
#[inline]
fn std_max(a: f64, b: f64) -> f64 {
    if a < b { b } else { a }
}

/// `(std::min)(a, b)`: `b < a ? b : a`, with `std_max`'s NaN asymmetry.
#[inline]
fn std_min(a: f64, b: f64) -> f64 {
    if b < a { b } else { a }
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
    let tail_sq = eigen_squared_norm(tail.len(), |i| tail[i]);
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
///
/// `tmp = essential^T * bottom` is Eigen's `GemvProduct`: the row-major
/// matrix-vector kernel for a block of two or more columns, and the runtime
/// fallback to `dot()` for a single column (`ProductEvaluators.h:380-384`).
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
    let len = essential.len();
    for (offset, slot) in tmp.iter_mut().enumerate() {
        let col = col0 + offset;
        let below = |b: usize| mat.at(row0 + 1 + b, col);
        let product = if cols == 1 {
            0.0 + eigen_dot(len, |b| essential[b], below)
        } else {
            eigen_gemv_row(len, below, |b| essential[b])
        };
        *slot = product + mat.at(row0, col);
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

/// Apply `H` to the tail of a vector starting at `row0`, the essential vector
/// being column `k` of `qr` below its diagonal.
///
/// For a vector Eigen's `essential^T * bottom` is an `InnerProduct`, so the
/// projection is [`eigen_dot`].
fn apply_householder_vector(w: &mut [f64], row0: usize, qr: &DenseMatrix, k: usize, tau: f64) {
    let rows = w.len() - row0;
    if rows == 1 {
        w[row0] *= 1.0 - tau;
        return;
    }
    if tau == 0.0 {
        return;
    }
    let len = rows - 1;
    let mut sum = eigen_dot(len, |b| qr.at(k + 1 + b, k), |b| w[row0 + 1 + b]);
    sum += w[row0];
    w[row0] -= tau * sum;
    for below in 0..len {
        let e = qr.at(k + 1 + below, k);
        w[row0 + 1 + below] -= tau * e * sum;
    }
}

/// Euclidean norm as Eigen's `MatrixBase::norm()` computes it:
/// `sqrt(squaredNorm())`, with no scaling pass, over `len` coefficients read
/// through `coeff`.
///
/// The third of the three norms on this path. `ColPivHouseholderQR` uses it,
/// and only it, for the initial column norms and for the direct recomputation
/// the LAPACK downdating rule falls back to - never `stableNorm` or
/// `blueNorm`.
fn plain_norm(len: usize, coeff: impl Fn(usize) -> f64) -> f64 {
    eigen_squared_norm(len, coeff).sqrt()
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
            updated.push(plain_norm(m, |row| qr.at(row, col)));
        }
        let mut direct = updated.clone();
        // `m_colNormsUpdated.maxCoeff()`: seeded with the first norm.
        let mut biggest = updated.first().copied().unwrap_or(0.0);
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
                    direct[j] = plain_norm(m - k - 1, |row| qr.at(k + 1 + row, j));
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
            apply_householder_vector(&mut out, k, &self.qr, k, self.tau[k]);
        }
        out
    }
}

/// `triangular_solve_vector<OnTheLeft, Upper, ColMajor>` for `size <= 16`
/// (one `EIGEN_TUNE_TRIANGULAR_PANEL_WIDTH` panel): back substitution that
/// divides the pivot and then subtracts `w[i] * a(j, i)` from every `j < i`,
/// leaving a zero right-hand side untouched. `a(i, j)` reads the triangle.
fn solve_col_major_upper(size: usize, a: impl Fn(usize, usize) -> f64, w: &mut [f64]) {
    for k in 0..size {
        let i = size - k - 1;
        if w[i] != 0.0 {
            w[i] /= a(i, i);
            let pivot = w[i];
            for j in 0..i {
                w[j] -= pivot * a(j, i);
            }
        }
    }
}

/// `triangular_solve_vector<OnTheLeft, Lower, RowMajor>` for `size <= 16`:
/// forward substitution whose inner sum `sum_j a(i, j) w[j]` is an
/// [`eigen_sum`] reduction, then division unless the right-hand side is zero.
fn solve_row_major_lower(size: usize, a: impl Fn(usize, usize) -> f64, w: &mut [f64]) {
    for i in 0..size {
        if i > 0 {
            let sum = eigen_sum(i, |j| a(i, j) * w[j]);
            w[i] -= sum;
        }
        if w[i] != 0.0 {
            w[i] /= a(i, i);
        }
    }
}

/// `triangular_solve_vector<OnTheLeft, Upper, RowMajor>` for `size <= 16`:
/// back substitution with an [`eigen_sum`] reduction over the `k` solved
/// coefficients to the right of row `i`.
fn solve_row_major_upper(size: usize, a: impl Fn(usize, usize) -> f64, w: &mut [f64]) {
    for k in 0..size {
        let i = size - k - 1;
        if k > 0 {
            let start = i + 1;
            let sum = eigen_sum(k, |j| a(i, start + j) * w[start + j]);
            w[i] -= sum;
        }
        if w[i] != 0.0 {
            w[i] /= a(i, i);
        }
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
    // `s.topLeftCorner(nsing, nsing).transpose().triangularView<Upper>()
    // .solveInPlace(wa.head(nsing))`: the transpose is row-major, so this is
    // `triangular_solve_vector<OnTheLeft, Upper, RowMajor>`
    // (`TriangularSolverVector.h:30-71`) - each row's inner sum is a reduction,
    // and a zero right-hand side skips its division.
    solve_row_major_upper(nsing, |i, j| s.at(j, i), &mut wa);
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
    // `triangularView<Upper>().solveInPlace` on the column-major R:
    // `triangular_solve_vector<OnTheLeft, Upper, ColMajor>`
    // (`TriangularSolverVector.h:74-118`) divides each pivot and then subtracts
    // its column from the rows above, skipping a zero right-hand side.
    solve_col_major_upper(rank, |i, j| qr.qr.at(i, j), &mut wa1);
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
        // `lmpar.h:198` spells this branch `P^-1 * diag.cwiseProduct(wa2) /
        // dxnorm`, so the product is formed first and the quotient taken
        // second. The Newton correction further down spells the same quantity
        // `P^-1 * diag.cwiseProduct(wa2 / dxnorm)` and therefore divides
        // first. The two associations differ in the last bit; each site keeps
        // the one its line of Eigen has. See the solver-deviation section of
        // `docs/DISTRIBUTION_FITTERS_SUPPORT.md`.
        let mut work: Vec<f64> = (0..n)
            .map(|j| (diag[qr.ind[j]] * wa2[qr.ind[j]]) / dxnorm)
            .collect();
        // `topLeftCorner(n, n).transpose().triangularView<Lower>()`: row-major
        // lower, `triangular_solve_vector<OnTheLeft, Lower, RowMajor>`.
        solve_row_major_lower(n, |i, j| qr.qr.at(j, i), &mut work);
        let temp = blue_norm(&work);
        parl = fp / delta / temp / temp;
    }
    let mut upper = vec![0.0f64; n];
    for (j, slot) in upper.iter_mut().enumerate() {
        let sum = eigen_dot(j + 1, |i| qr.qr.at(i, j), |i| qtb[i]);
        *slot = sum / diag[qr.ind[j]];
    }
    let gnorm = stable_norm(&upper);
    let mut paru = gnorm / delta;
    if paru == 0.0 {
        paru = f64::MIN_POSITIVE / std_min(delta, 0.1);
    }
    let mut par = std_min(std_max(par_in, parl), paru);
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
            par = std_max(f64::MIN_POSITIVE, 0.001 * paru);
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
        // `lmpar.h:241`: `P^-1 * diag.cwiseProduct(wa2 / dxnorm)` - the
        // quotient first here, unlike the `parl` branch above.
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
            parl = std_max(parl, par);
        }
        if fp < 0.0 {
            paru = std_min(paru, par);
        }
        par = std_max(parl, par + parc);
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
/// the vector back the same way and only then decide whether to throw. After a
/// [`LmStatus::TooManyFunctionEvaluation`] stop that is the last accepted
/// point, not the rejected trial.
///
/// The evaluation accounting is checked budget by budget against the executed
/// C++: for the eight `GaussTraceFitter`/`EGHTraceFitter` class-test fits and
/// the 50 trace fits of `FeatureFinderCentroided_1`, the status and the
/// residual and Jacobian evaluation counts equal Eigen's at every `max_fev`
/// from 1 to 500, and for four degenerate fits at 500
/// (`tests/lm_budget_differential.rs`); the statuses reached there are 1 to 5.
///
/// The arithmetic is checked point by point: for those 62 fits and 79 more,
/// every residual-evaluation argument, the final parameters and the counts are
/// bit-identical to Eigen 5.0.1 as the Linux x86_64 Release build of OpenMS
/// compiles it, and to the same Eigen on macOS arm64 with its FMA lanes
/// disabled (`tests/lm_eigen_path_differential.rs`). Against the product SDK on
/// arm64, whose Eigen fuses the packet multiply-adds, 21 of the 141 paths are
/// identical; see `docs/DISTRIBUTION_FITTERS_SUPPORT.md` section 1. Two limits
/// bound the claim, and no OpenMS caller reaches the first: the triangular
/// solves are Eigen's single-panel form, exact for at most 16 parameters; and
/// where a NaN enters a norm, Eigen's vectorized `maxCoeff` is itself
/// platform-dependent, so only the NaN is reproduced, not its sign or payload.
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
                let sum = eigen_dot(j + 1, |i| qr.qr.at(i, j), |i| qtf[i] / fnorm);
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
                delta = std_min(delta, pnorm);
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
            // `wa3.noalias() = fjac.triangularView<Upper>() * (...)` resizes
            // `wa3` to the Jacobian's `m` rows and zero-fills the `m - n` below
            // the triangle, so `wa3.stableNorm()` reduces `m` coefficients.
            let temp1 = {
                let padded = |i: usize| if i < n { projected_step[i] } else { 0.0 };
                let t = stable_norm_by(m, padded) / fnorm;
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
                delta = temp * std_min(delta, pnorm / 0.1);
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
            // `do { ... } while (ratio < Scalar(1e-4))`: only a ratio that
            // compares below the threshold retries inside this step, so a NaN
            // ratio leaves the loop and the Jacobian is evaluated again, as in
            // Eigen. Spelled without a negated comparison for `clippy`.
            if ratio >= 1e-4 || ratio.is_nan() {
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

    /// The two guards `stable_norm_kernel` puts around `1 / maxCoeff`, which
    /// the fitters themselves never reach because they validate their input
    /// finite: a subnormal largest coefficient, where the reciprocal
    /// overflows, and an infinite one, where it underflows to zero.
    #[test]
    fn stable_norm_guards_the_reciprocal_at_both_ends_of_the_range() {
        let subnormal = [f64::from_bits(1), f64::from_bits(1)];
        let got = stable_norm(&subnormal);
        assert!(got.is_finite(), "{got:?}");
        assert!(got >= 0.0);
        let infinite = [f64::INFINITY, 1.0];
        assert_eq!(stable_norm(&infinite), f64::INFINITY);
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

    /// `lmpar` scales `diag * x` by `1 / dxnorm` twice, and Eigen associates
    /// the two occurrences differently: `lmpar.h:198` computes
    /// `(diag * wa2) / dxnorm` for the `parl` lower bound, `lmpar.h:241`
    /// computes `diag * (wa2 / dxnorm)` for the Newton correction. The
    /// associations are not interchangeable in binary64, and this pins a
    /// triple where they differ, so that collapsing the two sites onto one
    /// spelling fails here rather than silently moving the last bits of every
    /// fitted parameter.
    #[test]
    fn the_two_lmpar_scalings_are_not_the_same_association() {
        let (diag, wa2, dxnorm) = (0.1_f64, 1.1_f64, 7.0_f64);
        assert_ne!((diag * wa2) / dxnorm, diag * (wa2 / dxnorm));
    }

    /// Left to right, `1e16 + 1` loses the one. Eigen's two-lane reduction
    /// (`Redux.h:275-322`) pairs term 0 with term 2 and term 1 with term 3
    /// before the lanes meet, so the cancellation happens first and both ones
    /// survive; the odd fifth term comes last.
    #[test]
    fn eigen_sum_pairs_the_terms_as_the_two_lane_reduction() {
        let terms = [1e16, 1.0, -1e16, 1.0, 1.0];
        let sequential = terms.iter().fold(0.0, |acc, t| acc + t);
        assert_eq!(sequential, 2.0);
        assert_eq!(eigen_sum(terms.len(), |i| terms[i]), 3.0);
        assert_eq!(eigen_sum(1, |_| 7.0), 7.0);
        assert_eq!(eigen_squared_norm(0, |_| unreachable!()), 0.0);
    }

    /// `inner_product_impl` (`InnerProduct.h:117-172`) seeds four two-lane
    /// accumulators with the first eight products, multiply-adds the ninth and
    /// tenth into the first, and folds the accumulators from the back. Here
    /// the `-1e16` meets the `1e16` before any one is added to either, so all
    /// eight ones survive; left to right only one does.
    #[test]
    fn eigen_dot_folds_four_accumulators_from_the_back() {
        let a = [1e16, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1e16, 1.0];
        let b = [1.0; 10];
        let sequential = a.iter().zip(&b).fold(0.0, |acc, (x, y)| acc + x * y);
        assert_eq!(sequential, 1.0);
        assert_eq!(eigen_dot(a.len(), |i| a[i], |i| b[i]), 8.0);
        assert_eq!(eigen_dot(1, |_| 3.0, |_| 2.0), 6.0);
    }

    /// The row-major matrix-vector kernel (`GeneralMatrixVector.h:298-462`)
    /// runs one accumulator per lane over the paired prefix and adds the odd
    /// tail after the lanes are combined.
    #[test]
    fn eigen_gemv_row_accumulates_one_lane_per_packet_slot() {
        let a = [1e16, 1.0, -1e16, 1.0, 1.0];
        let b = [1.0; 5];
        let sequential = a.iter().zip(&b).fold(0.0, |acc, (x, y)| acc + x * y);
        assert_eq!(sequential, 2.0);
        assert_eq!(eigen_gemv_row(a.len(), |i| a[i], |i| b[i]), 3.0);
    }

    /// Eigen's triangular vector solvers leave a zero right-hand side alone
    /// instead of dividing it, so a zero pivot above a zero entry does not
    /// turn into `0 / 0`.
    #[test]
    fn triangular_solves_skip_a_zero_right_hand_side() {
        // Upper triangle [[2, 1], [0, 0]], right-hand side [4, 0].
        let upper = [[2.0, 1.0], [0.0, 0.0]];
        let mut w = [4.0, 0.0];
        solve_col_major_upper(2, |i, j| upper[i][j], &mut w);
        assert_eq!(w, [2.0, 0.0]);
        let mut w = [0.0, 4.0];
        solve_row_major_lower(2, |i, j| upper[j][i], &mut w);
        assert!(w[1].is_infinite(), "{w:?}");
        assert_eq!(w[0], 0.0);
        let mut w = [4.0, 0.0];
        solve_row_major_upper(2, |i, j| upper[i][j], &mut w);
        assert_eq!(w, [2.0, 0.0]);
    }

    /// `std::max(a, b)` is `a < b ? b : a`: a NaN on the left is kept and one
    /// on the right is ignored, unlike `f64::max`, which drops either.
    #[test]
    fn std_min_and_max_keep_eigens_nan_asymmetry() {
        assert!(std_max(f64::NAN, 1.0).is_nan());
        assert_eq!(std_max(1.0, f64::NAN), 1.0);
        assert!(std_min(f64::NAN, 1.0).is_nan());
        assert_eq!(std_min(1.0, f64::NAN), 1.0);
        assert!(!f64::NAN.max(1.0).is_nan());
    }

    /// `stable_norm_kernel` seeds `maxCoeff` with the first magnitude: a NaN
    /// there becomes the scale and the norm; a NaN later is not a maximum but
    /// poisons the sum of squares unless every other entry is zero.
    #[test]
    fn stable_norm_propagates_nan_as_eigen_does() {
        assert!(stable_norm(&[f64::NAN, 1.0]).is_nan());
        assert!(stable_norm(&[1.0, f64::NAN]).is_nan());
        assert_eq!(stable_norm(&[0.0, f64::NAN]), 0.0);
        assert!(stable_norm(&[f64::NAN, 0.0]).is_nan());
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
