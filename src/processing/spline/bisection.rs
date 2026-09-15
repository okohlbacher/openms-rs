// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bisection search for the maximum of a spline.
//!
//! Ports `src/openms/include/OpenMS/MATH/MISC/SplineBisection.h`, a header-only
//! function template. See `docs/SPLINE_BISECTION_SUPPORT.md`.
//!
//! The C++ template accepts any type exposing `eval(double)` and
//! `derivative(double)`; the Rust equivalent is the
//! [`SplineFunction`](crate::processing::spline::bisection::SplineFunction)
//! trait, implemented by
//! [`CubicSpline2d`](crate::processing::spline::CubicSpline2d) and
//! [`BSpline2d`](crate::processing::spline::BSpline2d) — the two types the
//! source documents it as working with.

use crate::{Error, Result};

/// A spline that can be evaluated and differentiated once at a position.
///
/// The C++ `Math::spline_bisection` is a template whose only requirement on its
/// argument is the pair of member functions `eval(double) const` and
/// `derivative(double) const`. This trait is that requirement, made explicit.
/// `first_derivative` is spelled out rather than called `derivative` because
/// [`CubicSpline2d::derivative`](crate::processing::spline::CubicSpline2d::derivative)
/// already takes a derivative order, which the C++ `derivatives` overload
/// carries instead.
pub trait SplineFunction {
    /// Value of the spline at `x`.
    ///
    /// # Errors
    ///
    /// Implementations that are only defined on a bounded domain return
    /// [`Error::InvalidValue`] outside it, matching their C++ counterparts.
    fn eval(&self, x: f64) -> Result<f64>;

    /// First derivative of the spline at `x`.
    ///
    /// # Errors
    ///
    /// As [`SplineFunction::eval`].
    fn first_derivative(&self, x: f64) -> Result<f64>;
}

/// A shared reference to a spline is a spline.
///
/// Forwarding only; it adds no behaviour. It exists because
/// [`spline_bisection`] is generic, and a generic parameter does not get the
/// deref coercion that turns `&&T` into `&T` at an ordinary call. A caller
/// holding a borrowed spline — the peak picker, once it fits through
/// [`CubicSpline2dFitter`](crate::processing::spline::CubicSpline2dFitter),
/// holds `&CubicSpline2d` rather than an owned one — can therefore keep writing
/// `spline_bisection(&spline, ..)` unchanged.
impl<T: SplineFunction + ?Sized> SplineFunction for &T {
    fn eval(&self, x: f64) -> Result<f64> {
        (**self).eval(x)
    }
    fn first_derivative(&self, x: f64) -> Result<f64> {
        (**self).first_derivative(x)
    }
}

/// Iteration ceiling for [`spline_bisection`].
///
/// Native addition. The C++ `do`/`while` has no iteration cap: it halves the
/// bracket until `righthand - lefthand > threshold` fails. With a `threshold`
/// below the spacing of the two `f64` bracket ends the midpoint stops moving
/// and the loop never terminates, so this port bounds it. A bisection over the
/// whole `f64` range needs at most about 2 100 halvings, so the ceiling is
/// never reached by a search that the C++ would finish.
pub const MAX_BISECTION_STEPS: usize = 4096;

/// Locate the maximum of a spline inside a bracket by bisecting its first
/// derivative, returning `(position, value)`.
///
/// Port of `OpenMS::Math::spline_bisection`. The C++ writes its two results
/// through the `max_peak_mz` and `max_peak_int` out-parameters; this returns
/// them instead.
///
/// # Arguments
///
/// * `spline` — the curve to search; `left_neighbor` and `right_neighbor` must
///   lie inside its domain if it has one.
/// * `left_neighbor`, `right_neighbor` — the bracket ends, named after the
///   neighbouring m/z values the source passes in.
/// * `threshold` — bracket width at which to stop, the source's default is
///   `1e-6`. Must be finite and strictly positive.
///
/// # The termination condition, exactly as the source writes it
///
/// The loop is a `do`/`while`, so the body always runs at least once — even
/// when `right_neighbor <= left_neighbor`, where the width test is false from
/// the start. Each pass takes the midpoint, evaluates the first derivative
/// there and stops early when `!(|f'(mid)| > f64::EPSILON)`; that spelling also
/// stops on a `NaN` derivative, which a plain `<=` would not. Otherwise the
/// sign of the derivative selects the half to keep. The source seeds
/// `lefthand_sign = true` and never updates it, so the exclusive-or reduces to
/// "move the right end down when the derivative is negative, the left end up
/// otherwise"; a derivative of exactly `-0.0` counts as non-negative because
/// the source tests `< 0.0`. The reported position is the midpoint of the final
/// bracket, which after an early break is the bracket *before* the midpoint that
/// triggered it, not that midpoint.
///
/// # Notes
///
/// The search assumes a single interior maximum. When the true apex lies outside
/// the bracket the derivative never changes sign and the result converges onto
/// the bracket end nearest the apex; the class test pins that behaviour.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a bracket end or `threshold` is not
/// finite, when `threshold` is not strictly positive, or when
/// [`MAX_BISECTION_STEPS`] is exhausted. Errors raised by `spline` — most often
/// a query outside its domain — are propagated unchanged.
pub fn spline_bisection<T: SplineFunction + ?Sized>(
    spline: &T,
    left_neighbor: f64,
    right_neighbor: f64,
    threshold: f64,
) -> Result<(f64, f64)> {
    if !left_neighbor.is_finite() || !right_neighbor.is_finite() {
        return Err(Error::InvalidValue(
            "spline bisection bracket must be finite".into(),
        ));
    }
    if !threshold.is_finite() || threshold <= 0.0 {
        return Err(Error::InvalidValue(
            "spline bisection threshold must be finite and positive".into(),
        ));
    }

    let mut lefthand = left_neighbor;
    let mut righthand = right_neighbor;
    let eps = f64::EPSILON;

    // do { ... } while (righthand - lefthand > threshold)
    let mut steps = 0usize;
    loop {
        let mid = (lefthand + righthand) / 2.0;
        let midpoint_deriv_val = spline.first_derivative(mid)?;
        // The source writes `!(fabs(d) > eps)`, which is also true for a NaN
        // derivative; both halves of that are spelled out here.
        let magnitude = midpoint_deriv_val.abs();
        if magnitude <= eps || magnitude.is_nan() {
            break;
        }
        // `lefthand_sign` is fixed at true in the source, so this is
        // `true ^ midpoint_sign`, i.e. the negation of the midpoint sign.
        let midpoint_sign = midpoint_deriv_val >= 0.0;
        if !midpoint_sign {
            righthand = mid;
        } else {
            lefthand = mid;
        }
        steps += 1;
        if steps >= MAX_BISECTION_STEPS {
            return Err(Error::InvalidValue(
                "spline bisection did not reach its threshold within the iteration ceiling".into(),
            ));
        }
        // `while (righthand - lefthand > threshold)`, negated to a break.
        let width = righthand - lefthand;
        if width <= threshold || width.is_nan() {
            break;
        }
    }

    let max_peak_mz = (lefthand + righthand) / 2.0;
    let max_peak_int = spline.eval(max_peak_mz)?;
    Ok((max_peak_mz, max_peak_int))
}

/// The source's default bracket width for [`spline_bisection`].
pub const DEFAULT_BISECTION_THRESHOLD: f64 = 1e-6;

#[cfg(test)]
mod tests {
    use super::*;

    struct Parabola {
        peak: f64,
        height: f64,
        a: f64,
    }
    impl SplineFunction for Parabola {
        fn eval(&self, x: f64) -> Result<f64> {
            Ok(-self.a * (x - self.peak) * (x - self.peak) + self.height)
        }
        fn first_derivative(&self, x: f64) -> Result<f64> {
            Ok(-2.0 * self.a * (x - self.peak))
        }
    }

    // Values transcribed from the executed probe of the pinned C++ source
    // (tests/data/spline_math_cpp_probe.tsv, case `bisection`).
    #[test]
    fn reproduces_the_probe_brackets() {
        let cases: [(Parabola, f64, f64, f64, f64, f64); 8] = [
            (
                Parabola {
                    peak: 500.0,
                    height: 1000.0,
                    a: 1.0,
                },
                499.0,
                501.0,
                DEFAULT_BISECTION_THRESHOLD,
                500.0,
                1000.0,
            ),
            (
                Parabola {
                    peak: 500.3,
                    height: 750.0,
                    a: 2.0,
                },
                499.0,
                501.0,
                DEFAULT_BISECTION_THRESHOLD,
                500.2999997138977,
                749.9999999999999,
            ),
            (
                Parabola {
                    peak: 500.123456,
                    height: 100.0,
                    a: 5.0,
                },
                499.0,
                501.0,
                1e-9,
                500.12345599988475,
                100.0,
            ),
            (
                Parabola {
                    peak: 505.0,
                    height: 50.0,
                    a: 1.0,
                },
                499.0,
                501.0,
                DEFAULT_BISECTION_THRESHOLD,
                500.99999952316284,
                33.99999618530251,
            ),
            (
                Parabola {
                    peak: 495.0,
                    height: 50.0,
                    a: 1.0,
                },
                499.0,
                501.0,
                DEFAULT_BISECTION_THRESHOLD,
                499.00000047683716,
                33.99999618530251,
            ),
            (
                Parabola {
                    peak: 499.0,
                    height: 10.0,
                    a: 1.0,
                },
                499.0,
                501.0,
                DEFAULT_BISECTION_THRESHOLD,
                499.00000047683716,
                9.999999999999773,
            ),
            (
                Parabola {
                    peak: 500.0,
                    height: 1000.0,
                    a: 1.0,
                },
                501.0,
                499.0,
                DEFAULT_BISECTION_THRESHOLD,
                500.0,
                1000.0,
            ),
            (
                Parabola {
                    peak: 500.3,
                    height: 750.0,
                    a: 2.0,
                },
                499.0,
                501.0,
                0.1,
                500.28125,
                749.999296875,
            ),
        ];
        for (spline, left, right, threshold, mz, intensity) in cases {
            let (a, b) = spline_bisection(&spline, left, right, threshold).unwrap();
            assert_eq!(a, mz);
            assert_eq!(b, intensity);
        }
    }

    #[test]
    fn rejects_degenerate_arguments() {
        let s = Parabola {
            peak: 500.0,
            height: 1.0,
            a: 1.0,
        };
        assert!(spline_bisection(&s, f64::NAN, 501.0, 1e-6).is_err());
        assert!(spline_bisection(&s, 499.0, f64::INFINITY, 1e-6).is_err());
        assert!(spline_bisection(&s, 499.0, 501.0, 0.0).is_err());
        assert!(spline_bisection(&s, 499.0, 501.0, -1.0).is_err());
        assert!(spline_bisection(&s, 499.0, 501.0, f64::NAN).is_err());
    }

    #[test]
    fn a_threshold_below_the_float_spacing_hits_the_ceiling_instead_of_hanging() {
        // The bracket ends are already adjacent f64 values, and their exact
        // midpoint 1 + 2^-53 is a tie that rounds to the even mantissa, i.e.
        // back to the left end. An apex far to the right keeps the derivative
        // positive, so the source assigns `lefthand = mid` forever without the
        // bracket ever narrowing: the C++ loop does not terminate here.
        let left = 1.0f64;
        let right = f64::from_bits(left.to_bits() + 1);
        let s = Parabola {
            peak: 1.0e6,
            height: 1.0,
            a: 1.0,
        };
        assert!(right - left > 1e-30);
        assert_eq!((left + right) / 2.0, left);
        assert!(spline_bisection(&s, left, right, 1e-30).is_err());
    }
}
