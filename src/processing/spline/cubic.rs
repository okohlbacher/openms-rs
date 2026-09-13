// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Natural cubic-spline interpolation of a 2D data set.
//!
//! Ports `src/openms/include/OpenMS/MATH/MISC/CubicSpline2d.h` and
//! `src/openms/source/MATH/MISC/CubicSpline2d.cpp`. See
//! `docs/CUBIC_SPLINE2D_SUPPORT.md` for the API mapping and the evidence behind
//! every documented quirk.
//!
//! The type is re-exported as
//! [`crate::processing::peak_picking::CubicSpline2d`] because the peak picker
//! and the retention-time transformations were its first consumers.

use crate::processing::spline::bisection::SplineFunction;
use crate::{Error, Result};

fn checked(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::InvalidValue(
            "cubic spline arithmetic is not finite".into(),
        ))
    }
}

/// Natural cubic-spline interpolation of a 2D data set.
///
/// Fits a piecewise cubic polynomial through the supplied `(x, y)` knots so the
/// resulting spline is twice continuously differentiable. The construction
/// follows R. L. Burden and J. D. Faires, *Numerical Analysis*, 4th ed.,
/// PWS-Kent 1989, ISBN 0-53491-585-X, pp. 126-131, and this port reproduces the
/// C++ recurrence term by term and in its original evaluation order, because
/// floating-point addition is not associative.
///
/// Construction is the expensive step; the source notes that subsequent
/// evaluations are roughly fifty times cheaper. After construction the spline
/// can be sampled anywhere in the closed interval `[x_first, x_last]` with
/// [`CubicSpline2d::eval`], and its first three derivatives are available from
/// [`CubicSpline2d::derivative`].
///
/// # The boundary condition
///
/// The spline is *natural*: the second derivative vanishes at both ends. The
/// source sets the trailing quadratic coefficient `c_.back()` to zero and its
/// recurrence leaves `c_[0]` at zero as well.
///
/// The two ends are not equally exact, and the difference is observable.
/// `f''(x_first)` is `2*c_[0] + 6*d_[0]*0`, which is exactly `0.0` for every
/// input. `f''(x_last)` is evaluated on the *last segment* at its right end, as
/// `2*c_[n-1] + 6*d_[n-1]*h`, and although `d_[n-1]` is
/// `(c_[n] - c_[n-1]) / (3*h)` with `c_[n]` exactly zero, three roundings stand
/// between that and `2*c_[n]`. It is therefore zero only up to rounding: the
/// probe records exactly `0` for the class test's uniform sine grid but
/// `-3.814697265625e-06` for the upstream peak, whose interior second
/// derivatives are of order `1e11`. The C++ has the same behaviour, and
/// `CubicSpline2d_test.cpp` asserts the condition with a tolerance.
///
/// # Differences from the source
///
/// * The C++ accepts *non-decreasing* abscissae, so two equal `x` values pass
///   its check and then divide by a zero interval width, producing a spline of
///   `NaN` coefficients that no later call reports. This port requires strictly
///   increasing, finite abscissae and rejects the rest up front.
/// * Every intermediate of the recurrence is checked for finiteness, so an
///   overflow is an error rather than a silently poisoned spline.
#[derive(Clone, Debug)]
pub struct CubicSpline2d {
    x: Vec<f64>,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
}

impl CubicSpline2d {
    /// Default ceiling on the number of knots, used by
    /// [`CubicSpline2d::new`] and [`CubicSpline2d::from_pairs`].
    ///
    /// Native addition: the source allocates five vectors proportional to the
    /// knot count with no bound at all.
    pub const MAX_POINTS: usize = 1_000_000;

    /// Build the spline from parallel `x` / `y` slices, with at most
    /// [`CubicSpline2d::MAX_POINTS`] knots.
    ///
    /// `x` holds the knot abscissae and must contain at least two entries;
    /// `y` holds the knot ordinates and must have the same length.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the slices differ in length, when
    /// fewer than two knots are supplied, when a value is not finite, or when
    /// `x` is not strictly increasing — the three conditions the C++ constructor
    /// reports as `Exception::IllegalArgument`, plus the two the C++ does not
    /// check at all (see the type-level note on non-decreasing abscissae).
    pub fn new(x: &[f64], y: &[f64]) -> Result<Self> {
        Self::with_max_points(x, y, Self::MAX_POINTS)
    }

    /// Build the spline from parallel `x` / `y` slices with an explicit knot
    /// ceiling.
    ///
    /// `max_points` must be at least two. Native addition; see
    /// [`CubicSpline2d::MAX_POINTS`].
    ///
    /// # Errors
    ///
    /// As [`CubicSpline2d::new`], and additionally when `x` holds more than
    /// `max_points` knots.
    pub fn with_max_points(x: &[f64], y: &[f64], max_points: usize) -> Result<Self> {
        if x.len() != y.len() || x.len() < 2 || max_points < 2 || x.len() > max_points {
            return Err(Error::InvalidValue(
                "spline needs matching arrays of 2..=max_points knots".into(),
            ));
        }
        if x.iter().chain(y).any(|v| !v.is_finite()) {
            return Err(Error::InvalidValue("spline knots must be finite".into()));
        }
        if x.windows(2).any(|p| p[0] >= p[1]) {
            return Err(Error::InvalidValue(
                "spline coordinates must be strictly increasing".into(),
            ));
        }
        let n = x.len() - 1;
        let h: Vec<_> = x.windows(2).map(|p| p[1] - p[0]).collect();
        if h.iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidValue("spline knot spacing overflows".into()));
        }
        let mut mu = vec![0.0; n];
        let mut z = vec![0.0; n];
        for i in 1..n {
            let span = checked(x[i + 1] - x[i - 1])?;
            let l = checked(2.0 * span - h[i - 1] * mu[i - 1])?;
            mu[i] = checked(h[i] / l)?;
            z[i] = checked(
                (3.0 * (y[i + 1] * h[i - 1] - y[i] * span + y[i - 1] * h[i]) / (h[i - 1] * h[i])
                    - h[i - 1] * z[i - 1])
                    / l,
            )?;
        }
        let mut b = vec![0.0; n];
        let mut c = vec![0.0; n + 1];
        let mut d = vec![0.0; n];
        for j in (0..n).rev() {
            c[j] = checked(z[j] - mu[j] * c[j + 1])?;
            b[j] = checked((y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0)?;
            d[j] = checked((c[j + 1] - c[j]) / (3.0 * h[j]))?;
        }
        Ok(Self {
            x: x.to_vec(),
            a: y[..n].to_vec(),
            b,
            c,
            d,
        })
    }

    /// Build the spline from unordered `(x, y)` knots, the port of the C++
    /// `CubicSpline2d(const std::map<double, double>&)` constructor.
    ///
    /// The C++ takes a `std::map` keyed by `x`, so the container has already
    /// sorted the knots by abscissa and collapsed repeats by the time the
    /// constructor sees them: `std::map::insert` keeps the *first* value stored
    /// for a key. This function reproduces both effects — it sorts by `x` and,
    /// for equal abscissae, keeps the ordinate that appears first in `points` —
    /// so that a caller with a sequence of pairs gets the spline the C++ would
    /// have built from the equivalent map.
    ///
    /// As a consequence, and unlike [`CubicSpline2d::new`], repeated abscissae
    /// are not an error here: they are a single knot.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when fewer than two *distinct* abscissae
    /// remain, matching the C++ check on the already-deduplicated map size, when
    /// a coordinate is not finite, or when more than
    /// [`CubicSpline2d::MAX_POINTS`] pairs are supplied.
    pub fn from_pairs(points: &[(f64, f64)]) -> Result<Self> {
        if points.len() > Self::MAX_POINTS {
            return Err(Error::InvalidValue(
                "spline knot count exceeds the maximum".into(),
            ));
        }
        if points.iter().any(|p| !p.0.is_finite() || !p.1.is_finite()) {
            return Err(Error::InvalidValue("spline knots must be finite".into()));
        }
        let mut order: Vec<usize> = (0..points.len()).collect();
        // Stable by construction: equal abscissae keep their input order, so the
        // first one wins the deduplication below, as std::map::insert does.
        order.sort_by(|&l, &r| points[l].0.total_cmp(&points[r].0));
        let mut x = Vec::with_capacity(points.len());
        let mut y = Vec::with_capacity(points.len());
        for index in order {
            let (px, py) = points[index];
            if x.last().is_some_and(|&last: &f64| last == px) {
                continue;
            }
            x.push(px);
            y.push(py);
        }
        if x.len() < 2 {
            return Err(Error::InvalidValue(
                "spline needs two or more distinct abscissae".into(),
            ));
        }
        Self::with_max_points(&x, &y, Self::MAX_POINTS)
    }

    /// Closed interval the spline is defined on, as `(first_knot, last_knot)`.
    ///
    /// Native accessor; the C++ keeps its knot vector private.
    pub fn domain(&self) -> (f64, f64) {
        (self.x[0], self.x[self.x.len() - 1])
    }

    /// Number of cubic segments, one fewer than the number of knots.
    ///
    /// Native accessor.
    pub fn segment_count(&self) -> usize {
        self.a.len()
    }

    fn interval(&self, x: f64) -> Result<usize> {
        let (lo, hi) = self.domain();
        if !x.is_finite() || x < lo || x > hi {
            return Err(Error::InvalidValue("query outside spline domain".into()));
        }
        // At a knot use the interval starting at that knot; at the last knot use the preceding interval.
        Ok(self
            .x
            .partition_point(|&v| v <= x)
            .saturating_sub(1)
            .min(self.a.len() - 1))
    }

    /// Evaluate the spline at `x`.
    ///
    /// `x` must lie in the closed interval `[first_knot, last_knot]`. Inside a
    /// segment the value is Horner's form of the stored coefficients,
    /// `((d*t + c)*t + b)*t + a` with `t = x - x_i`, exactly as the C++ writes
    /// it. At an interior knot the segment starting there is used, so the value
    /// is that knot's ordinate; at the last knot the preceding segment is used
    /// and evaluated at its right end.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `x` lies outside the knot range —
    /// the C++ `Exception::IllegalArgument` — when `x` is not finite, which the
    /// C++ does not check, or when the evaluation overflows.
    pub fn eval(&self, x: f64) -> Result<f64> {
        let i = self.interval(x)?;
        let dx = x - self.x[i];
        checked(((self.d[i] * dx + self.c[i]) * dx + self.b[i]) * dx + self.a[i])
    }

    /// Derivative of the spline at `x` of order one, two or three.
    ///
    /// This single method covers both C++ overloads: `derivative(x)` is
    /// `derivative(x, 1)`, and `derivatives(x, order)` is this call. A cubic
    /// spline has meaningful first, second and third derivatives; higher orders
    /// are zero everywhere and are rejected, as in the source.
    ///
    /// The third derivative is piecewise constant and therefore discontinuous at
    /// the knots; the first and second are continuous there, which is what makes
    /// the fit a cubic spline rather than merely piecewise cubic.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `x` lies outside the knot range or
    /// is not finite, and when `order` is not 1, 2 or 3. The range is checked
    /// before the order, as in the C++.
    pub fn derivative(&self, x: f64, order: u8) -> Result<f64> {
        let i = self.interval(x)?;
        let dx = x - self.x[i];
        checked(match order {
            1 => self.b[i] + 2.0 * self.c[i] * dx + 3.0 * self.d[i] * dx * dx,
            2 => 2.0 * self.c[i] + 6.0 * self.d[i] * dx,
            3 => 6.0 * self.d[i],
            _ => {
                return Err(Error::InvalidValue(
                    "spline derivative order must be 1, 2, or 3".into(),
                ));
            }
        })
    }

    /// OpenMS derivative bisection for a peak bracket, assuming a positive
    /// derivative on its left side. This is not a global spline maximizer.
    ///
    /// Native predecessor of
    /// [`spline_bisection`](crate::processing::spline::spline_bisection), kept
    /// because the peak picker depends on its extra guarantees: it validates
    /// that both bracket ends are inside the domain, stops when the midpoint
    /// stops moving, and fails rather than looping when 128 halvings are not
    /// enough. Prefer `spline_bisection` for a faithful port of
    /// `SplineBisection.h`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a bracket end lies outside the
    /// domain, when `left > right`, when `tolerance` is not finite and positive,
    /// or when the search does not converge in 128 iterations.
    pub fn peak_maximum(
        &self,
        mut left: f64,
        mut right: f64,
        tolerance: f64,
    ) -> Result<(f64, f64)> {
        self.interval(left)?;
        self.interval(right)?;
        if left > right || !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(Error::InvalidValue(
                "invalid spline maximum bracket or tolerance".into(),
            ));
        }
        for _ in 0..128 {
            let mid = left / 2.0 + right / 2.0;
            let derivative = self.derivative(mid, 1)?;
            if derivative.abs() <= f64::EPSILON || mid == left || mid == right {
                return Ok((mid, self.eval(mid)?));
            }
            if derivative < 0.0 {
                right = mid;
            } else {
                left = mid;
            }
            if right - left <= tolerance {
                let mid = left / 2.0 + right / 2.0;
                return Ok((mid, self.eval(mid)?));
            }
        }
        Err(Error::InvalidValue(
            "spline maximum did not converge within 128 iterations".into(),
        ))
    }
}

impl SplineFunction for CubicSpline2d {
    fn eval(&self, x: f64) -> Result<f64> {
        CubicSpline2d::eval(self, x)
    }
    fn first_derivative(&self, x: f64) -> Result<f64> {
        self.derivative(x, 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn upstream() -> (Vec<f64>, Vec<f64>) {
        (
            vec![
                486.784, 486.787, 486.79, 486.793, 486.795, 486.797, 486.8, 486.802, 486.805,
                486.808, 486.811,
            ],
            vec![
                0.0, 154683.17, 620386.5, 1701390.12, 2848879.25, 3564045.5, 2744585.7, 1605583.0,
                1518984.0, 1591352.21, 1691345.1,
            ],
        )
    }

    // Full-precision values from the executed probe of the pinned C++ source,
    // case `cubic_upstream`; the class test prints the same numbers truncated.
    #[test]
    fn matches_the_cpp_probe_bit_for_bit() {
        let (x, y) = upstream();
        let s = CubicSpline2d::new(&x, &y).unwrap();
        assert_eq!(s.eval(486.785).unwrap(), 35173.18417789844);
        assert_eq!(s.eval(486.794).unwrap(), 2271426.9331624084);
        assert_eq!(s.derivative(486.785, 1).unwrap(), 39270152.29962468);
        assert_eq!(s.derivative(486.785, 2).unwrap(), 12290904368.273579);
        assert_eq!(s.derivative(486.794, 1).unwrap(), 594825947.1542642);
        assert_eq!(s.derivative(486.794, 2).unwrap(), 7415503644.895798);
        assert_eq!(s.derivative(486.784, 3).unwrap(), 12290904367865.562);
        assert_eq!(s.eval(486.7855).unwrap(), 56600.683880853714);
        assert_eq!(s.derivative(486.7855, 2).unwrap(), 18436356552.06104);
    }

    #[test]
    fn the_map_constructor_sorts_and_keeps_the_first_ordinate_per_abscissa() {
        let (x, y) = upstream();
        let mut pairs: Vec<(f64, f64)> = x.iter().copied().zip(y.iter().copied()).collect();
        pairs.reverse();
        // A repeat of an abscissa already present is dropped, as std::map does.
        pairs.push((486.79, -1.0));
        let s = CubicSpline2d::from_pairs(&pairs).unwrap();
        assert_eq!(s.eval(486.785).unwrap(), 35173.18417789844);
        assert_eq!(s.eval(486.79).unwrap(), 620386.5);
        assert_eq!(s.segment_count(), 10);
        assert!(CubicSpline2d::from_pairs(&[(1.0, 1.0), (1.0, 2.0)]).is_err());
        assert!(CubicSpline2d::from_pairs(&[(1.0, 1.0)]).is_err());
    }

    #[test]
    fn the_natural_boundary_condition_holds_at_both_ends() {
        let n = 10usize;
        let (x_min, x_max) = (-0.5f64, 1.5f64);
        let x: Vec<f64> = (0..=n)
            .map(|i| x_min + (i as f64) / 10.0 * (x_max - x_min))
            .collect();
        let y: Vec<f64> = x.iter().map(|v| v.sin()).collect();
        let s = CubicSpline2d::new(&x, &y).unwrap();
        // The first knot is exact for every input: `mu[0]` and `z[0]` are never
        // written, so `c[0]` is exactly zero and the reported value is `2*c[0]`.
        assert_eq!(s.derivative(x[0], 2).unwrap(), 0.0);
        // The last knot is evaluated on the last segment at its right end, so it
        // is only zero up to rounding. On this uniform grid it lands on exact
        // zero, which the probe also records.
        assert_eq!(s.derivative(x[n], 2).unwrap(), 0.0);

        // On the non-uniform upstream peak it does not, and the value the probe
        // records is asserted here rather than hidden behind a tolerance: it is
        // 3e-17 of the interior second derivatives, but it is not zero.
        let (ux, uy) = upstream();
        let peak = CubicSpline2d::new(&ux, &uy).unwrap();
        assert_eq!(peak.derivative(ux[0], 2).unwrap(), 0.0);
        assert_eq!(
            peak.derivative(ux[ux.len() - 1], 2).unwrap(),
            -3.814697265625e-06
        );
    }

    #[test]
    fn first_and_second_derivatives_are_continuous_across_every_knot() {
        let (x, y) = upstream();
        let s = CubicSpline2d::new(&x, &y).unwrap();
        for knot in &x[1..x.len() - 1] {
            let step = 1e-9;
            for order in [1u8, 2] {
                let left = s.derivative(knot - step, order).unwrap();
                let right = s.derivative(knot + step, order).unwrap();
                let scale = left.abs().max(right.abs()).max(1.0);
                assert!(
                    (left - right).abs() <= 1e-4 * scale,
                    "order {order} jumps at {knot}: {left} vs {right}"
                );
            }
        }
    }

    #[test]
    fn rejects_what_the_source_rejects_and_what_it_silently_accepts() {
        let (x, y) = upstream();
        let s = CubicSpline2d::new(&x, &y).unwrap();
        assert!(s.eval(486.783).is_err());
        assert!(s.eval(486.812).is_err());
        assert!(s.eval(f64::NAN).is_err());
        assert!(s.derivative(486.79, 0).is_err());
        assert!(s.derivative(486.79, 4).is_err());
        assert!(CubicSpline2d::new(&[1.0, 0.0, 2.0], &[1.0, 2.0, 3.0]).is_err());
        assert!(CubicSpline2d::new(&[1.0], &[1.0]).is_err());
        assert!(CubicSpline2d::new(&[1.0, 2.0], &[1.0]).is_err());
        // The C++ accepts this (its check is "non-decreasing") and then divides
        // by a zero interval width, producing NaN coefficients.
        assert!(CubicSpline2d::new(&[0.0, 1.0, 1.0, 2.0], &[0.0, 1.0, 2.0, 3.0]).is_err());
    }
}
