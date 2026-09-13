// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Cubic B-spline smoothing over a uniform node grid.
//!
//! Ports `src/openms/include/OpenMS/MATH/MISC/BSpline2d.h` and
//! `src/openms/source/MATH/MISC/BSpline2d.cpp`. In C++ that pair is a PIMPL
//! wrapper around the vendored `eol-bspline` template library
//! (`src/openms/extern/eol-bspline/BSpline/`, UCAR, BSD-3-Clause); this crate
//! takes no third-party dependency, so the library's algorithm is reimplemented
//! here from that pinned source. See `docs/BSPLINE2D_SUPPORT.md` for what that
//! means for the evidence, including the executed probe the numbers are checked
//! against.
//!
//! The formulation is Ooyama's cubic B-spline (Monthly Weather Review 115,
//! October 1987): the curve is a sum of basis functions centred on `M + 1`
//! equally spaced nodes, and the coefficients solve a diagonally banded system
//! `(P + Q) a = b` where `P` is the least-squares normal matrix over the data
//! and `Q` a second-derivative penalty scaled by the cutoff wavelength.
//!
//! Use [`crate::processing::spline::CubicSpline2d`] instead when the curve must
//! pass through the knots; this type fits, it does not interpolate.

use crate::processing::spline::bisection::SplineFunction;
use crate::{Error, Result};

/// Truncated value of pi used by the source's `BSplineBase<T>::PI`.
///
/// The vendored library hard-codes `3.1415927` rather than using
/// `std::numbers::pi`, and that value feeds the wavelength-to-penalty
/// conversion, so reproducing it is required for agreement past the seventh
/// significant digit. Do not replace it with [`std::f64::consts::PI`].
// The whole point of this constant is that it is not pi.
#[allow(clippy::approx_constant)]
const EOL_PI: f64 = 3.1415927;

/// Beta coefficients of the boundary constraint, indexed by boundary condition
/// and then by the node position 0, 1, M-1, M.
const BOUNDARY_CONDITIONS: [[f64; 4]; 3] = [
    [-4.0, -1.0, -1.0, -4.0],
    [0.0, 1.0, 1.0, 0.0],
    [2.0, -1.0, -1.0, 2.0],
];

/// Integrals of the products of the second derivatives of two normalised basis
/// functions `m` nodes apart, per unit domain. This is the `K == 2` plane of the
/// source's `qparts[3][4][4]`; `K` is fixed at 2 by the `BSplineBase`
/// constructor, so the other two planes are unreachable.
const QPARTS_K2: [[f64; 4]; 4] = [
    [0.75, 2.25, 2.25, 0.75],
    [0.0, -1.125, -1.125, -1.125],
    [0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.375],
];

/// Boundary condition applied at the two ends of the node domain.
///
/// The source warns not to change the constants because they are passed through
/// to the B-spline implementation; that coupling is internal to this port, but
/// the discriminants are kept identical so a caller reading a stored integer
/// still maps to the same case.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum BoundaryCondition {
    /// Set the endpoints of the spline to zero (`BC_ZERO_ENDPOINTS`, 0).
    ZeroEndpoints = 0,
    /// Set the first derivative of the spline to zero at the endpoints
    /// (`BC_ZERO_FIRST`, 1).
    ZeroFirst = 1,
    /// Set the second derivative to zero (`BC_ZERO_SECOND`, 2). The default,
    /// as in the C++ constructor signature.
    #[default]
    ZeroSecond = 2,
}

impl BoundaryCondition {
    /// Row of the beta table this condition selects, equal to the C++ enum value.
    pub fn index(self) -> usize {
        self as usize
    }
}

/// Diagonally banded square matrix holding three bands either side of the
/// diagonal, the storage the source calls `BandedMatrix`.
///
/// Reads and writes outside the band are silently ignored and read back as
/// zero. That is not a convenience: the C++ `BandedMatrix::element` returns a
/// shared `out_of_bounds` scratch reference for such coordinates, so a write
/// lands nowhere and a read yields whatever was written there last. Every
/// access the ported algorithms make is provably inside the band, so the two
/// behaviours coincide; see `docs/BSPLINE2D_SUPPORT.md`.
#[derive(Clone, Debug)]
struct Band {
    n: i64,
    data: Vec<f64>,
}

impl Band {
    fn new(n: i64) -> Self {
        Self {
            n,
            data: vec![0.0; (n.max(0) as usize) * 7],
        }
    }
    fn slot(&self, i: i64, j: i64) -> Option<usize> {
        let d = j - i;
        if i < 0 || j < 0 || i >= self.n || j >= self.n || !(-3..=3).contains(&d) {
            None
        } else {
            Some((i as usize) * 7 + (d + 3) as usize)
        }
    }
    fn get(&self, i: i64, j: i64) -> f64 {
        self.slot(i, j).map_or(0.0, |k| self.data[k])
    }
    fn set(&mut self, i: i64, j: i64, value: f64) {
        if let Some(k) = self.slot(i, j) {
            self.data[k] = value;
        }
    }
    fn add(&mut self, i: i64, j: i64, value: f64) {
        if let Some(k) = self.slot(i, j) {
            self.data[k] += value;
        }
    }
    /// One-based access, matching the `A(i, j)` spelling of the source's LU code.
    fn at(&self, i: i64, j: i64) -> f64 {
        self.get(i - 1, j - 1)
    }
    fn put(&mut self, i: i64, j: i64, value: f64) {
        self.set(i - 1, j - 1, value);
    }
}

/// Everything the basis functions need: the node grid and the penalty weight.
#[derive(Clone, Debug)]
struct Domain {
    xmin: f64,
    xmax: f64,
    m: i64,
    dx: f64,
    alpha: f64,
    wave_length: f64,
    bc: usize,
}

impl Domain {
    /// `BSplineBase<T>::Beta`.
    fn beta(&self, m: i64) -> f64 {
        if m > 1 && m < self.m - 1 {
            return 0.0;
        }
        let mut k = m;
        if k >= self.m - 1 {
            k -= self.m - 3;
        }
        if !(0..4).contains(&k) {
            // Unreachable for 0 <= m <= M, which is the only range the ported
            // callers use. The C++ asserts here and reads out of bounds when
            // NDEBUG is set; see docs/BSPLINE2D_SUPPORT.md.
            return 0.0;
        }
        BOUNDARY_CONDITIONS[self.bc][k as usize]
    }

    /// `BSplineBase<T>::Basis`: the closed cubic basis function at node `m`.
    fn basis(&self, m: i64, x: f64) -> f64 {
        let mut y = 0.0;
        let xm = self.xmin + (m as f64) * self.dx;
        let mut z = ((x - xm) / self.dx).abs();
        if z < 2.0 {
            z = 2.0 - z;
            y = 0.25 * (z * z * z);
            z -= 1.0;
            if z > 0.0 {
                y -= z * z * z;
            }
        }
        if m == 0 || m == 1 {
            y += self.beta(m) * self.basis(-1, x);
        } else if m == self.m - 1 || m == self.m {
            y += self.beta(m) * self.basis(self.m + 1, x);
        }
        y
    }

    /// `BSplineBase<T>::DBasis`: the derivative of [`Domain::basis`].
    fn dbasis(&self, m: i64, x: f64) -> f64 {
        let mut dy = 0.0;
        let xm = self.xmin + (m as f64) * self.dx;
        let delta = (x - xm) / self.dx;
        let mut z = delta.abs();
        if z < 2.0 {
            z = 2.0 - z;
            dy = 0.25 * z * z;
            z -= 1.0;
            if z > 0.0 {
                dy -= z * z;
            }
            dy *= (if delta > 0.0 { -1.0 } else { 1.0 }) * 3.0 / self.dx;
        }
        if m == 0 || m == 1 {
            dy += self.beta(m) * self.dbasis(-1, x);
        } else if m == self.m - 1 || m == self.m {
            dy += self.beta(m) * self.dbasis(self.m + 1, x);
        }
        dy
    }

    /// `BSplineBase<T>::qDelta`.
    fn q_delta(&self, m1: i64, m2: i64) -> f64 {
        let (m1, m2) = if m1 > m2 { (m2, m1) } else { (m1, m2) };
        if m2 - m1 > 3 {
            return 0.0;
        }
        let row = match usize::try_from(m2 - m1).ok().and_then(|r| QPARTS_K2.get(r)) {
            Some(row) => row,
            None => return 0.0,
        };
        let mut q = 0.0;
        let mut m = (m1 - 2).max(0);
        let end = (m1 + 2).min(self.m);
        while m < end {
            if let Some(value) = usize::try_from(m - m1 + 2).ok().and_then(|c| row.get(c)) {
                q += value;
            }
            m += 1;
        }
        q * self.alpha
    }

    /// `BSplineBase<T>::calculateQ`: the derivative-constraint matrix.
    fn calculate_q(&self) -> Band {
        let mut q = Band::new(self.m + 1);
        if self.alpha == 0.0 {
            return q;
        }
        for i in 0..=self.m {
            q.set(i, i, self.q_delta(i, i));
            let mut j = 1;
            while j < 4 && i + j <= self.m {
                let value = self.q_delta(i, i + j);
                q.set(i, i + j, value);
                q.set(i + j, i, value);
                j += 1;
            }
        }
        // Upper-left boundary block. `b1`, `b2` and the accumulator are `float`
        // in the source, and the accumulator is rounded to `float` after every
        // `+=`; that rounding is reproduced because it is observable.
        for i in 0..=1i64 {
            let b1 = self.beta(i) as f32;
            for j in i..i + 4 {
                if j > self.m {
                    // The C++ evaluates Beta(j) past the end of the beta table
                    // here and writes the result to a banded coordinate that
                    // discards it. Skipping is numerically identical.
                    continue;
                }
                let b2 = self.beta(j) as f32;
                let mut acc: f32 = 0.0;
                if i + 1 < 4 {
                    acc = (f64::from(acc) + f64::from(b2) * self.q_delta(-1, i)) as f32;
                }
                if j + 1 < 4 {
                    acc = (f64::from(acc) + f64::from(b1) * self.q_delta(-1, j)) as f32;
                }
                acc = (f64::from(acc) + f64::from(b1 * b2) * self.q_delta(-1, -1)) as f32;
                let value = q.get(i, j) + f64::from(acc);
                q.set(i, j, value);
                q.set(j, i, value);
            }
        }
        // Lower-right boundary block.
        for i in self.m - 1..=self.m {
            let b1 = self.beta(i) as f32;
            for j in i - 3..=i {
                if j < 0 {
                    // Same discarded-write argument as above; here the C++ reads
                    // BoundaryConditions[bc][-1].
                    continue;
                }
                let b2 = self.beta(j) as f32;
                let mut acc: f32 = 0.0;
                if self.m + 1 - i < 4 {
                    acc = (f64::from(acc) + f64::from(b2) * self.q_delta(i, self.m + 1)) as f32;
                }
                if self.m + 1 - j < 4 {
                    acc = (f64::from(acc) + f64::from(b1) * self.q_delta(j, self.m + 1)) as f32;
                }
                acc = (f64::from(acc) + f64::from(b1 * b2) * self.q_delta(self.m + 1, self.m + 1))
                    as f32;
                let value = q.get(i, j) + f64::from(acc);
                q.set(i, j, value);
                q.set(j, i, value);
            }
        }
        q
    }

    /// `BSplineBase<T>::addP`: add the least-squares normal matrix in place.
    ///
    /// The source holds each basis value and each product in a `float` before
    /// widening it back to `double` for the accumulation; that truncation is
    /// worth roughly seven significant digits and is reproduced.
    fn add_p(&self, q: &mut Band, x: &[f64]) {
        for &xi in x {
            let mx = ((xi - self.xmin) / self.dx) as i64;
            let hi = mx.saturating_add(2).min(self.m);
            let mut m = mx.saturating_sub(1).max(0);
            while m <= hi {
                let pm = self.basis(m, xi) as f32;
                let sum = pm * pm;
                q.add(m, m, f64::from(sum));
                let mut n = m + 1;
                while n <= hi {
                    let pn = self.basis(n, xi) as f32;
                    let sum = pm * pn;
                    q.add(m, n, f64::from(sum));
                    q.add(n, m, f64::from(sum));
                    n += 1;
                }
                m += 1;
            }
        }
    }
}

/// `LU_factor_banded` from the source's `BandedMatrix.h`, Crout's algorithm
/// restricted to the three bands either side of the diagonal.
fn lu_factor_banded(a: &mut Band) -> bool {
    let n = a.n;
    for j in 1..=n {
        if a.at(j, j) == 0.0 {
            return false;
        }
        let start = if j > 3 { j - 3 } else { 1 };
        for i in start..=j {
            let mut sum = 0.0;
            let mut k = start;
            while k < i {
                sum += a.at(i, k) * a.at(k, j);
                k += 1;
            }
            let value = a.at(i, j) - sum;
            a.put(i, j, value);
        }
        let mut i = j + 1;
        while i <= n && i <= j + 3 {
            let mut sum = 0.0;
            let mut k = if i > 3 { i - 3 } else { 1 };
            while k < j {
                sum += a.at(i, k) * a.at(k, j);
                k += 1;
            }
            let value = (a.at(i, j) - sum) / a.at(j, j);
            a.put(i, j, value);
            i += 1;
        }
    }
    true
}

/// `LU_solve_banded`: forward then backward substitution, in place on `b`.
fn lu_solve_banded(a: &Band, b: &mut [f64]) -> bool {
    let n = a.n;
    if n <= 0 || b.len() != n as usize {
        return false;
    }
    for i in 2..=n {
        let mut sum = b[(i - 1) as usize];
        let mut j = if i > 3 { i - 3 } else { 1 };
        while j < i {
            sum -= a.at(i, j) * b[(j - 1) as usize];
            j += 1;
        }
        b[(i - 1) as usize] = sum;
    }
    // The source divides by A(M, M) here without checking it, unlike the loop
    // below; a zero pivot yields an infinity that the finiteness check in
    // BSpline2d::solve turns into an error.
    b[(n - 1) as usize] /= a.at(n, n);
    let mut i = n - 1;
    while i >= 1 {
        if a.at(i, i) == 0.0 {
            return false;
        }
        let mut sum = b[(i - 1) as usize];
        let mut j = i + 1;
        while j <= n && j <= i + 3 {
            sum -= a.at(i, j) * b[(j - 1) as usize];
            j += 1;
        }
        b[(i - 1) as usize] = sum / a.at(i, i);
        i -= 1;
    }
    true
}

/// Cubic B-spline fitted to scattered data by least squares, with an optional
/// second-derivative penalty acting as a low-pass filter.
///
/// The curve is *not* an interpolant: it is the least-squares cubic B-spline
/// over a uniform grid of nodes spanning the data, so it does not generally pass
/// through the input points. The cutoff wavelength adds Ooyama's derivative
/// constraint, which suppresses structure shorter than that wavelength.
///
/// # Evaluation outside the fitted domain
///
/// Evaluation is defined everywhere and never fails on range. Each basis
/// function is zero more than two node intervals from its node, so past that
/// distance the sum is empty and [`BSpline2d::eval`] returns exactly the mean of
/// the ordinates it was fitted to, while [`BSpline2d::derivative`] returns zero.
/// Between the last node and that point the curve decays towards the mean.
/// This is the source's behaviour, reproduced deliberately: contrast
/// [`CubicSpline2d`](crate::processing::spline::CubicSpline2d), which rejects
/// any query outside its knot range.
///
/// # Differences from the source
///
/// * The C++ constructor always yields an object and reports failure through
///   `ok()`, with `eval` then returning zero. Here a domain that cannot be set
///   up or factored is an `Err` from the constructor, so a value of this type
///   always started life usable. [`BSpline2d::ok`] remains, because
///   [`BSpline2d::solve`] can still fail later, and it keeps the source's
///   "evaluate to zero once not ok" contract.
/// * Non-finite inputs and non-finite solutions are rejected rather than
///   propagated.
#[derive(Clone, Debug)]
pub struct BSpline2d {
    domain: Domain,
    x: Vec<f64>,
    lu: Band,
    coefficients: Vec<f64>,
    mean: f64,
    ok: bool,
}

impl BSpline2d {
    /// Maximum number of data points accepted by a constructor.
    ///
    /// Native addition: the source bounds nothing. Chosen together with
    /// [`BSpline2d::MAX_NODES`] so that the automatic node count of a
    /// wavelength-free fit, `2n + 1`, still fits.
    pub const MAX_POINTS: usize = 250_000;

    /// Maximum number of nodes, `M + 1`, accepted or derived.
    ///
    /// Native addition. The banded matrix holds seven `f64` per node, so this
    /// caps that allocation at 28 MiB.
    pub const MAX_NODES: usize = 500_001;

    /// Fit a spline with the source's defaults: no derivative constraint, a
    /// zero second derivative at the endpoints, and an automatic node count.
    ///
    /// Equivalent to `BSpline2d(x, y)` in C++.
    ///
    /// # Errors
    ///
    /// As [`BSpline2d::with_options`].
    pub fn new(x: &[f64], y: &[f64]) -> Result<Self> {
        Self::with_options(x, y, 0.0, BoundaryCondition::ZeroSecond, 0)
    }

    /// Fit a spline, setting up the node domain and solving for `y` in one step.
    ///
    /// # Arguments
    ///
    /// * `x` — the abscissae of the domain. They need not be sorted; only their
    ///   minimum and maximum set the node grid, and repeats are allowed.
    /// * `y` — the ordinates, one per `x`, in the same order.
    /// * `wavelength` — the cutoff wavelength in the units of `x`. Zero is
    ///   documented as disabling the derivative constraint; see the note below,
    ///   because it does not.
    /// * `boundary_condition` — the constraint at the two ends of the node grid.
    /// * `num_nodes` — the number of nodes for the cubic B-spline. Below two, a
    ///   node count is derived from the data and the cutoff wavelength.
    ///
    /// # Notes
    ///
    /// A `wavelength` of zero does not disable the derivative constraint, even
    /// though both the OpenMS and the eol-bspline documentation say so: the
    /// setup rewrites a zero wavelength to `1.0` in the units of `x` before the
    /// penalty weight is computed, so the weight is `(1 / (2 pi dx))^4`, which
    /// for a dense grid is large rather than zero. The behaviour is reproduced;
    /// the claim is not.
    ///
    /// With an automatic node count and a zero wavelength the grid gets `2n`
    /// intervals for `n` points, which is where the fit comes closest to
    /// interpolating. With an explicit `num_nodes` the wavelength is likewise
    /// rewritten to `1.0`, so two fits that differ only in `num_nodes` carry
    /// different penalty weights.
    ///
    /// When a cutoff wavelength is given and no node count is, the source
    /// requires the wavelength to be shorter than the span of `x` and searches
    /// for a grid with at least two — preferably four — intervals per
    /// wavelength while keeping at least one, preferably two, points per
    /// interval. That search is ported unchanged, including its use of the
    /// wavelength-to-interval ratio from the *last* trial when deciding whether
    /// to continue.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `x` and `y` differ in length or are
    /// empty, when a value is not finite, when `wavelength` is negative or not
    /// finite, when all abscissae are equal so the node spacing would be zero,
    /// when the wavelength exceeds the span of `x` or the node search cannot
    /// keep one point per interval (both `setDomain` failures in C++), when
    /// fewer than three nodes result — the banded matrix refuses to be shaped
    /// and the C++ leaves `ok()` false — when the matrix cannot be factored, or
    /// when [`BSpline2d::MAX_POINTS`] or [`BSpline2d::MAX_NODES`] is exceeded.
    pub fn with_options(
        x: &[f64],
        y: &[f64],
        wavelength: f64,
        boundary_condition: BoundaryCondition,
        num_nodes: usize,
    ) -> Result<Self> {
        if x.len() != y.len() {
            return Err(Error::InvalidValue(
                "B-spline x and y must have the same length".into(),
            ));
        }
        if x.is_empty() {
            return Err(Error::InvalidValue(
                "B-spline needs at least one point".into(),
            ));
        }
        if x.len() > Self::MAX_POINTS {
            return Err(Error::InvalidValue(
                "B-spline point count exceeds the maximum".into(),
            ));
        }
        if num_nodes > Self::MAX_NODES {
            return Err(Error::InvalidValue(
                "B-spline node count exceeds the maximum".into(),
            ));
        }
        if x.iter().chain(y).any(|v| !v.is_finite()) {
            return Err(Error::InvalidValue("B-spline data must be finite".into()));
        }
        if !wavelength.is_finite() || wavelength < 0.0 {
            return Err(Error::InvalidValue(
                "B-spline cutoff wavelength must be finite and non-negative".into(),
            ));
        }

        let domain = setup(x, wavelength, num_nodes, boundary_condition)?;
        let mut lu = domain.calculate_q();
        domain.add_p(&mut lu, x);
        if !lu_factor_banded(&mut lu) {
            return Err(Error::InvalidValue(
                "B-spline normal matrix could not be factored".into(),
            ));
        }

        let mut spline = Self {
            coefficients: vec![0.0; (domain.m + 1) as usize],
            domain,
            x: x.to_vec(),
            lu,
            mean: 0.0,
            ok: true,
        };
        spline.solve(y)?;
        Ok(spline)
    }

    /// Solve the spline curve again for a new set of ordinates over the same
    /// domain.
    ///
    /// `y` must have one value per abscissa the spline was constructed with;
    /// the C++ states that as a precondition which release builds do not check
    /// and which then reads past the end of the vector.
    ///
    /// The mean of `y` is subtracted before the fit and added back on
    /// evaluation, exactly as in the source, so the reported value far outside
    /// the domain is that mean.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `y` has the wrong length, holds a
    /// non-finite value, or when the banded back-substitution fails or produces
    /// a non-finite coefficient. On any error the spline is left not
    /// [`ok`](BSpline2d::ok) with zeroed coefficients, which is the state the
    /// C++ leaves behind, and further evaluation returns zero rather than a
    /// half-updated curve.
    pub fn solve(&mut self, y: &[f64]) -> Result<()> {
        if y.len() != self.x.len() {
            return Err(Error::InvalidValue(
                "B-spline solve needs one ordinate per abscissa".into(),
            ));
        }
        if y.iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidValue(
                "B-spline ordinates must be finite".into(),
            ));
        }
        let mut a = vec![0.0; (self.domain.m + 1) as usize];
        let mut mean = 0.0;
        for &value in y {
            mean += value;
        }
        mean /= y.len() as f64;

        for (j, &xj) in self.x.iter().enumerate() {
            let yj = y[j] - mean;
            let mx = ((xj - self.domain.xmin) / self.domain.dx) as i64;
            let hi = mx.saturating_add(2).min(self.domain.m);
            let mut m = mx.saturating_sub(1).max(0);
            while m <= hi {
                a[m as usize] += yj * self.domain.basis(m, xj);
                m += 1;
            }
        }

        let solved = lu_solve_banded(&self.lu, &mut a);
        if !solved || a.iter().any(|v| !v.is_finite()) {
            self.ok = false;
            self.coefficients = vec![0.0; (self.domain.m + 1) as usize];
            self.mean = 0.0;
            return Err(Error::InvalidValue(
                "B-spline coefficients could not be solved".into(),
            ));
        }
        self.coefficients = a;
        self.mean = mean;
        self.ok = true;
        Ok(())
    }

    /// Evaluate the smoothed curve at `x`.
    ///
    /// Returns zero when the spline is not [`ok`](BSpline2d::ok), as the source
    /// does. See the type-level note for what happens outside the fitted domain.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `x` is not finite — the C++ casts
    /// it to `int` to pick a node, which is undefined behaviour for `NaN` — or
    /// when the sum overflows.
    pub fn eval(&self, x: f64) -> Result<f64> {
        if !x.is_finite() {
            return Err(Error::InvalidValue(
                "B-spline evaluation position must be finite".into(),
            ));
        }
        if !self.ok {
            return Ok(0.0);
        }
        let mut y = 0.0;
        for i in self.window(x) {
            y += self.coefficients[i as usize] * self.domain.basis(i, x);
        }
        y += self.mean;
        finite(y, "B-spline evaluation is not finite")
    }

    /// First derivative of the spline curve at `x`.
    ///
    /// Returns zero when the spline is not [`ok`](BSpline2d::ok), and also more
    /// than two node intervals outside the fitted domain, where every basis
    /// derivative vanishes. The fitted mean is not added here, because the
    /// derivative of a constant is zero.
    ///
    /// # Errors
    ///
    /// As [`BSpline2d::eval`].
    pub fn derivative(&self, x: f64) -> Result<f64> {
        if !x.is_finite() {
            return Err(Error::InvalidValue(
                "B-spline evaluation position must be finite".into(),
            ));
        }
        if !self.ok {
            return Ok(0.0);
        }
        let mut dy = 0.0;
        for i in self.window(x) {
            dy += self.coefficients[i as usize] * self.domain.dbasis(i, x);
        }
        finite(dy, "B-spline derivative is not finite")
    }

    /// Whether the last fit succeeded.
    ///
    /// A constructed value starts `true`, because a failed setup is an `Err`
    /// rather than a not-ok object here; a failed [`BSpline2d::solve`] sets it
    /// `false` and makes evaluation return zero, as in the source.
    pub fn ok(&self) -> bool {
        self.ok
    }

    /// Number of nodes, one more than the number of node intervals.
    ///
    /// Native accessor for `BSplineBase::nNodes`, which the OpenMS wrapper
    /// hides behind its PIMPL pointer.
    pub fn node_count(&self) -> usize {
        (self.domain.m + 1) as usize
    }

    /// Spacing between neighbouring nodes, `BSplineBase::DX`.
    ///
    /// Native accessor.
    pub fn node_spacing(&self) -> f64 {
        self.domain.dx
    }

    /// Node domain as `(Xmin, Xmax)`, where `Xmax` is reconstructed as
    /// `Xmin + M * DX` exactly as `BSplineBase::Xmax` does.
    ///
    /// Native accessor. The reconstructed upper end can differ from the largest
    /// abscissa in the last bits, which is why it is reported rather than the
    /// stored maximum.
    pub fn domain(&self) -> (f64, f64) {
        (
            self.domain.xmin,
            self.domain.xmin + (self.domain.m as f64) * self.domain.dx,
        )
    }

    /// Smallest and largest abscissa actually supplied, as the source's
    /// protected `xmin` and `xmax` members hold them.
    ///
    /// Native accessor. The upper end differs from the one
    /// [`BSpline2d::domain`] reports whenever `Xmin + M * DX` does not round
    /// back to it exactly.
    pub fn data_range(&self) -> (f64, f64) {
        (self.domain.xmin, self.domain.xmax)
    }

    /// Weight of the derivative constraint, `BSplineBase::Alpha`.
    ///
    /// Native accessor. Zero means the penalty matrix is skipped entirely; see
    /// the note on [`BSpline2d::with_options`] about why a zero cutoff
    /// wavelength does not produce a zero weight.
    pub fn alpha(&self) -> f64 {
        self.domain.alpha
    }

    /// Cutoff wavelength actually in force, after the setup's rewrite of zero
    /// to `1.0`.
    ///
    /// Native accessor.
    pub fn wavelength(&self) -> f64 {
        self.domain.wave_length
    }

    /// Mean of the ordinates of the last fit, added back by
    /// [`BSpline2d::eval`].
    ///
    /// Native accessor for the source's protected `mean` member.
    pub fn fitted_mean(&self) -> f64 {
        self.mean
    }

    /// Basis coefficient `n`, from zero to `node_count() - 1`.
    ///
    /// Returns zero outside that range and when the spline is not
    /// [`ok`](BSpline2d::ok), as `BSpline<T>::coefficient` does.
    pub fn coefficient(&self, n: usize) -> f64 {
        if self.ok {
            self.coefficients.get(n).copied().unwrap_or(0.0)
        } else {
            0.0
        }
    }

    /// Node indices whose basis function can be non-zero at `x`.
    fn window(&self, x: f64) -> std::ops::RangeInclusive<i64> {
        let n = ((x - self.domain.xmin) / self.domain.dx) as i64;
        let lo = n.saturating_sub(1).max(0);
        let hi = n.saturating_add(2).min(self.domain.m);
        lo..=hi
    }
}

impl SplineFunction for BSpline2d {
    fn eval(&self, x: f64) -> Result<f64> {
        BSpline2d::eval(self, x)
    }
    fn first_derivative(&self, x: f64) -> Result<f64> {
        self.derivative(x)
    }
}

fn finite(value: f64, message: &str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::InvalidValue(message.into()))
    }
}

/// `BSplineBase<T>::Ratiod`: points per node interval, and, through its second
/// return, nodes per cutoff wavelength for the trial interval count.
fn ratios(ni: i64, xmin: f64, xmax: f64, wave_length: f64, nx: usize) -> (f64, f64) {
    let deltax = (xmax - xmin) / ni as f64;
    let ratiof = wave_length / deltax;
    let ratiod = nx as f64 / (ni as f64 + 1.0);
    (ratiod, ratiof)
}

/// `BSplineBase<T>::Setup` followed by `Alpha`: choose the node grid.
fn setup(
    x: &[f64],
    wavelength: f64,
    num_nodes: usize,
    boundary_condition: BoundaryCondition,
) -> Result<Domain> {
    let nx = x.len();
    let mut xmin = x[0];
    let mut xmax = x[0];
    for &value in &x[1..] {
        if value < xmin {
            xmin = value;
        } else if value > xmax {
            xmax = value;
        }
    }
    if xmax <= xmin {
        return Err(Error::InvalidValue(
            "B-spline needs at least two distinct abscissae".into(),
        ));
    }

    let mut wave_length = wavelength;
    let mut ni: i64 = 9;
    if num_nodes >= 2 {
        ni = num_nodes as i64 - 1;
        if wave_length == 0.0 {
            wave_length = 1.0;
        }
    } else if wave_length == 0.0 {
        ni = (nx as i64).saturating_mul(2);
        wave_length = 1.0;
    } else if wave_length > xmax - xmin {
        return Err(Error::InvalidValue(
            "B-spline cutoff wavelength exceeds the span of the abscissae".into(),
        ));
    } else {
        const FMIN: f64 = 2.0;
        let ceiling = BSpline2d::MAX_NODES as i64;
        // Raise the interval count until there are at least two intervals per
        // cutoff wavelength, while at least one point per interval remains.
        loop {
            ni += 1;
            if ni > ceiling {
                return Err(Error::InvalidValue(
                    "B-spline node search exceeded the node maximum".into(),
                ));
            }
            let (ratiod, ratiof) = ratios(ni, xmin, xmax, wave_length, nx);
            if ratiod < 1.0 {
                return Err(Error::InvalidValue(
                    "B-spline cannot keep one point per node interval at this cutoff wavelength"
                        .into(),
                ));
            }
            // `while (ratiof < fmin)` in the source, negated to a break; a
            // non-comparable ratio ends the loop there too.
            if ratiof >= FMIN || ratiof.is_nan() {
                break;
            }
        }
        // Keep raising it towards four intervals per wavelength while at least
        // two points per interval remain and the grid is not already far finer
        // than the wavelength needs.
        loop {
            ni += 1;
            if ni > ceiling {
                return Err(Error::InvalidValue(
                    "B-spline node search exceeded the node maximum".into(),
                ));
            }
            let (ratiod, ratiof) = ratios(ni, xmin, xmax, wave_length, nx);
            if ratiod < 1.0 || ratiof > 15.0 {
                ni -= 1;
                break;
            }
            if !(ratiof < 4.0 || ratiod > 2.0) {
                break;
            }
        }
    }

    if ni < 2 {
        // M + 1 < 3 leaves the banded matrix unshaped in C++, after which
        // factoring fails and ok() stays false.
        return Err(Error::InvalidValue(
            "B-spline needs at least three nodes".into(),
        ));
    }
    if ni + 1 > BSpline2d::MAX_NODES as i64 {
        return Err(Error::InvalidValue(
            "B-spline node count exceeds the maximum".into(),
        ));
    }

    let m = ni;
    let dx = (xmax - xmin) / m as f64;
    if !dx.is_finite() || dx <= 0.0 {
        return Err(Error::InvalidValue(
            "B-spline node spacing is not positive and finite".into(),
        ));
    }
    // Alpha(waveLength) with the derivative-constraint degree K fixed at 2.
    let a = wave_length / ((2.0 * EOL_PI) * dx);
    let a = a * a;
    let alpha = a * a;
    if !alpha.is_finite() {
        return Err(Error::InvalidValue(
            "B-spline derivative constraint weight overflows".into(),
        ));
    }
    Ok(Domain {
        xmin,
        xmax,
        m,
        dx,
        alpha,
        wave_length,
        bc: boundary_condition.index(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn squares() -> (Vec<f64>, Vec<f64>) {
        let x: Vec<f64> = (0..8).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|v| v * v).collect();
        (x, y)
    }

    // Every expected value in this module comes from the executed probe of the
    // pinned C++ sources; see tests/data/spline_math_cpp_probe.tsv.
    #[test]
    fn small_automatic_grid_matches_the_probe() {
        let (x, y) = squares();
        let s = BSpline2d::new(&x, &y).unwrap();
        assert_eq!(s.node_count(), 17);
        assert_eq!(s.domain(), (0.0, 7.0));
        assert_eq!(s.alpha(), 0.017513311466234586);
        assert_eq!(s.fitted_mean(), 17.5);
        assert_eq!(s.coefficient(13), 9.83198038112768);
        assert_eq!(s.coefficient(16), 20.997538675730194);
        assert_eq!(s.eval(-0.5).unwrap(), 4.158566027203834);
        assert_eq!(s.eval(0.0).unwrap(), -0.0036919460307487384);
        assert_eq!(s.eval(1.25).unwrap(), 1.5417880342313275);
        assert_eq!(s.eval(3.5).unwrap(), 12.247107986791324);
        assert_eq!(s.eval(7.0).unwrap(), 48.99630801359529);
        assert_eq!(s.derivative(-0.5).unwrap(), -25.053447218911195);
        assert_eq!(s.derivative(3.5).unwrap(), 6.999999496223682);
        assert_eq!(s.derivative(7.0).unwrap(), 13.411085010792291);
    }

    #[test]
    fn beyond_two_node_intervals_the_curve_is_the_fitted_mean() {
        let (x, y) = squares();
        let s = BSpline2d::new(&x, &y).unwrap();
        // Independently derived: mean of 0,1,4,...,49 is 140/8.
        assert_eq!(s.fitted_mean(), 17.5);
        assert_eq!(s.eval(-2.0).unwrap(), 17.5);
        assert_eq!(s.eval(9.0).unwrap(), 17.5);
        assert_eq!(s.derivative(-2.0).unwrap(), 0.0);
        assert_eq!(s.derivative(9.0).unwrap(), 0.0);
    }

    #[test]
    fn explicit_node_counts_and_boundary_conditions_match_the_probe() {
        let (x, y) = squares();
        let n4 = BSpline2d::with_options(&x, &y, 0.0, BoundaryCondition::ZeroSecond, 4).unwrap();
        assert_eq!(n4.node_count(), 4);
        assert_eq!(n4.alpha(), 2.1645785961380018e-05);
        assert_eq!(n4.coefficient(0), -11.833571244441133);
        assert_eq!(n4.coefficient(3), 20.83361959531383);
        assert_eq!(n4.eval(1.25).unwrap(), 1.9513238463546756);
        assert_eq!(n4.eval(9.0).unwrap(), 64.35243573004475);
        assert_eq!(n4.derivative(3.5).unwrap(), 6.999548900296078);

        let bc0 =
            BSpline2d::with_options(&x, &y, 0.0, BoundaryCondition::ZeroEndpoints, 6).unwrap();
        assert_eq!(bc0.alpha(), 0.00016701995340571017);
        assert_eq!(bc0.coefficient(0), -11.757589211611233);
        assert_eq!(bc0.coefficient(5), 20.460038997247644);
        // BC_ZERO_ENDPOINTS pulls the curve back to the fitted mean at the ends.
        assert_eq!(bc0.eval(0.0).unwrap(), 17.5);
        assert_eq!(bc0.eval(7.0).unwrap(), 17.499999999999936);
        assert_eq!(bc0.eval(1.25).unwrap(), 1.2255450471611553);
        assert_eq!(bc0.derivative(-0.5).unwrap(), -50.753565622482505);
    }

    #[test]
    fn a_cutoff_wavelength_selects_its_own_node_count() {
        let x: Vec<f64> = (0..=20).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|v| v * v).collect();

        let s = BSpline2d::with_options(&x, &y, 5.0, BoundaryCondition::ZeroSecond, 0).unwrap();
        assert_eq!(s.node_count(), 17);
        assert_eq!(s.wavelength(), 5.0);
        assert_eq!(s.alpha(), 0.16425570636886425);
        assert_eq!(s.coefficient(0), -91.40721658156676);
        assert_eq!(s.coefficient(16), 175.25945222994875);
        assert_eq!(s.eval(3.25).unwrap(), 10.586981824450163);
        assert_eq!(s.eval(17.5).unwrap(), 306.351466543918);
        assert_eq!(s.derivative(10.0).unwrap(), 20.000000641486807);

        // A longer wavelength keeps fewer nodes, and BC_ZERO_FIRST flattens the
        // curve at both ends: an invariant of the boundary condition, not a
        // transcribed number.
        let flat = BSpline2d::with_options(&x, &y, 12.0, BoundaryCondition::ZeroFirst, 0).unwrap();
        assert_eq!(flat.node_count(), 12);
        assert_eq!(flat.derivative(0.0).unwrap(), 0.0);
        assert!(flat.derivative(20.0).unwrap().abs() < 1e-9);
        assert_eq!(flat.eval(10.0).unwrap(), 99.50513840705506);
    }

    #[test]
    fn setup_failures_the_source_reports_through_ok() {
        let x = [0.0, 1.0, 2.0, 3.0];
        let y = [0.0, 1.0, 4.0, 9.0];
        // Cutoff wavelength longer than the span of the abscissae.
        assert!(BSpline2d::with_options(&x, &y, 100.0, BoundaryCondition::ZeroSecond, 0).is_err());
        // The node search starts at ten intervals, so an automatic node count
        // with a positive wavelength needs at least eleven points to keep one
        // point per interval. Eight points fail, as the probe records.
        let (sx, sy) = squares();
        assert!(BSpline2d::with_options(&sx, &sy, 3.0, BoundaryCondition::ZeroSecond, 0).is_err());
    }

    #[test]
    fn solve_replaces_the_curve_and_keeps_the_domain() {
        let (x, y) = squares();
        let mut s = BSpline2d::new(&x, &y).unwrap();
        let nodes = s.node_count();
        let flipped: Vec<f64> = y.iter().map(|v| -v).collect();
        s.solve(&flipped).unwrap();
        assert!(s.ok());
        assert_eq!(s.node_count(), nodes);
        assert_eq!(s.fitted_mean(), -17.5);
        // By linearity of the normal equations, negating every ordinate negates
        // every coefficient and the fitted mean: derived, not transcribed.
        assert_eq!(s.eval(-2.0).unwrap(), -17.5);
        assert_eq!(s.eval(3.5).unwrap(), -12.247107986791324);
        // A rejected solve leaves the previous curve intact.
        assert!(s.solve(&y[..3]).is_err());
        assert!(s.ok());
        assert_eq!(s.eval(3.5).unwrap(), -12.247107986791324);
    }

    #[test]
    fn rejects_degenerate_input() {
        assert!(BSpline2d::new(&[], &[]).is_err());
        assert!(BSpline2d::new(&[1.0, 2.0], &[1.0]).is_err());
        assert!(BSpline2d::new(&[1.0, 1.0, 1.0], &[1.0, 2.0, 3.0]).is_err());
        assert!(BSpline2d::new(&[1.0, f64::NAN], &[1.0, 2.0]).is_err());
        let (x, y) = squares();
        assert!(BSpline2d::with_options(&x, &y, -1.0, BoundaryCondition::ZeroSecond, 0).is_err());
        assert!(
            BSpline2d::with_options(&x, &y, f64::NAN, BoundaryCondition::ZeroSecond, 0).is_err()
        );
        // Two nodes leave the C++ banded matrix unshaped and ok() false.
        assert!(BSpline2d::with_options(&x, &y, 0.0, BoundaryCondition::ZeroSecond, 2).is_err());
        let s = BSpline2d::new(&x, &y).unwrap();
        assert!(s.eval(f64::NAN).is_err());
        assert!(s.derivative(f64::INFINITY).is_err());
        assert_eq!(s.coefficient(10_000), 0.0);
    }
}
