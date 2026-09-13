// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Andreas Bertsch, Chris Bielow, OpenMS Rust contributors $

//! Least-squares fit of a Gaussian to a set of (x, y) points.
//!
//! Ported from `MATH/STATISTICS/GaussFitter.h` and its `.cpp` at revision
//! `bc9cc12`. The fitted model is `A * exp(-(x - x0)^2 / (2 sigma^2))`, with
//! `A` an amplitude rather than a normalizing constant: see
//! [`GaussFitResult::eval`](crate::math::fitters::gauss::GaussFitResult::eval).
//!
//! The residual, the analytic Jacobian, the initial guess and the termination
//! rule are the source's; the optimizer is
//! [`levenberg_marquardt`](crate::math::fitters::levenberg_marquardt), which
//! reproduces the `Eigen::LevenbergMarquardt` the C++ calls. See
//! `docs/DISTRIBUTION_FITTERS_SUPPORT.md`.

use super::levenberg_marquardt::{DenseMatrix, LmParameters, LmStatus, minimize, preflight_points};
use crate::{Error, Result};

/// `0.5 * ln(2 * pi)`, the source's per-instance `halflogtwopi` member.
///
/// In C++ this is a non-static data member of the result struct initialized to
/// `0.5 * log(2.0 * Constants::PI)`, so every instance carries its own copy of
/// a constant. Here it is a constant.
const HALF_LOG_TWO_PI: f64 = 0.918_938_533_204_672_7;

/// `sqrt(2 * pi)`, the divisor Boost's normal density actually forms.
///
/// `boost/math/distributions/normal.hpp` writes the last step of `pdf` as
/// `result /= sd * sqrt(2 * constants::pi<RealType>())`: it evaluates the
/// square root of the rounded `2 * pi` at run time and does **not** use
/// `constants::root_two_pi`. That distinction is worth a constant of its own,
/// because Boost's `root_two_pi` decimal literal rounds to
/// `2.506_628_274_631_000_7`, one unit in the last place *above* this value,
/// and using it costs `GaussFitter::eval` up to one unit in the last place on
/// every point. The unit test below pins this constant to
/// `(2.0 * PI).sqrt()`, which is the expression the source evaluates.
const SQRT_TWO_PI: f64 = 2.506_628_274_631_000_2;

/// Number of fitted parameters: amplitude, center and width.
const PARAMETERS: usize = 3;

/// Amplitude, center and width of a Gaussian.
///
/// Doubles as the initial guess and as the fit result, exactly like the
/// source's `GaussFitter::GaussFitResult`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GaussFitResult {
    /// Parameter `A`: the height of the model at its center, not an integral.
    pub a: f64,
    /// Parameter `x0`: the center position.
    pub x0: f64,
    /// Parameter `sigma`: the width.
    ///
    /// [`GaussFitter::fit`](crate::math::fitters::gauss::GaussFitter::fit)
    /// always reports a non-negative width; the optimizer is unconstrained and
    /// a negative width describes the same curve.
    pub sigma: f64,
}

impl Default for GaussFitResult {
    /// `(-1, -1, -1)`, the source's default constructor. These are deliberately
    /// invalid values: `sigma = -1` is not a usable width, so a default-built
    /// result is a marker, not a model.
    fn default() -> Self {
        Self {
            a: -1.0,
            x0: -1.0,
            sigma: -1.0,
        }
    }
}

impl GaussFitResult {
    /// A result with the given amplitude, center and width.
    ///
    /// Nothing is validated here, matching the source's three-argument
    /// constructor; the consumers - `eval`, `log_eval_no_normalize` and `fit` -
    /// check what they need.
    pub fn new(a: f64, x0: f64, sigma: f64) -> Self {
        Self { a, x0, sigma }
    }

    /// Evaluate the density model at `x`.
    ///
    /// Returns the intensity, that is the normal probability density scaled so
    /// that its maximum at `x0` equals `a`. The source spells this out as
    /// `pdf(x) * (A / pdf(x0))`, with a comment warning that "simply
    /// multiplying the CDF with A is wrong"; the same two-step form is used
    /// here, in the same order, so the result is bit-comparable.
    ///
    /// This may be called with any parameters - the initial guess, to see a
    /// before-fit state, or the fitted ones.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `sigma` is not finite and strictly
    /// positive, or when `a`, `x0` or `x` is not finite. The source hands these
    /// to `boost::math::normal_distribution`, whose default policy throws
    /// `std::domain_error` on a non-positive or non-finite scale; a non-finite
    /// `x` is not rejected there but would yield NaN, which this port refuses
    /// to return silently.
    pub fn eval(&self, x: f64) -> Result<f64> {
        self.check_density()?;
        if !x.is_finite() {
            return Err(bad("Gaussian evaluation point must be finite"));
        }
        Ok(
            normal_pdf(x, self.x0, self.sigma)
                * (self.a / normal_pdf(self.x0, self.x0, self.sigma)),
        )
    }

    /// Evaluate the log density of the unit-amplitude Gaussian at `x`.
    ///
    /// Computes `-ln(sigma) - 0.5 ln(2 pi) - 0.5 ((x - x0) / sigma)^2`. The
    /// amplitude `a` deliberately takes no part: this is the log of a proper
    /// normal density, which is what a likelihood needs, and is why the source
    /// calls it "no normalize".
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `sigma` is not finite and strictly
    /// positive or when `x0` or `x` is not finite. The source takes `log` of a
    /// negative or zero `sigma` and returns NaN or negative infinity without
    /// comment.
    pub fn log_eval_no_normalize(&self, x: f64) -> Result<f64> {
        self.check_density()?;
        if !x.is_finite() {
            return Err(bad("Gaussian evaluation point must be finite"));
        }
        let scaled = (x - self.x0) / self.sigma;
        // The source writes `pow(scaled, 2.0)`; for a correctly rounded `pow`
        // that is the same value as `scaled * scaled`.
        Ok(-self.sigma.ln() - HALF_LOG_TWO_PI - 0.5 * (scaled * scaled))
    }

    fn check_density(&self) -> Result<()> {
        if !self.sigma.is_finite() || self.sigma <= 0.0 {
            return Err(bad("Gaussian sigma must be finite and positive"));
        }
        if !self.x0.is_finite() || !self.a.is_finite() {
            return Err(bad("Gaussian parameters must be finite"));
        }
        Ok(())
    }
}

/// Boost's `normal_distribution` density, in its arithmetic order.
///
/// Statement for statement `boost/math/distributions/normal.hpp`: form the
/// deviation, negate-and-square it in place, divide by `2 * sd * sd`,
/// exponentiate, then divide by `sd * sqrt(2 * pi)`. The last divisor is
/// [`SQRT_TWO_PI`], the value Boost computes there - not its `root_two_pi`
/// literal, which is a different `f64`.
fn normal_pdf(x: f64, mean: f64, sd: f64) -> f64 {
    let mut exponent = x - mean;
    exponent *= -exponent;
    exponent /= 2.0 * sd * sd;
    let mut result = exponent.exp();
    result /= sd * SQRT_TWO_PI;
    result
}

/// Fits a Gaussian to a set of points by non-linear least squares.
///
/// The initial guess matters: the residual surface is not convex and the
/// source's default guess, `(A, x0, sigma) = (0.06, 3.0, 0.5)`, is carried over
/// verbatim so that a fit started without
/// [`set_initial_parameters`](crate::math::fitters::gauss::GaussFitter::set_initial_parameters)
/// lands where the C++ lands.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GaussFitter {
    initial: GaussFitResult,
}

impl Default for GaussFitter {
    /// The source's constructor guess `(0.06, 3.0, 0.5)`.
    fn default() -> Self {
        Self {
            initial: GaussFitResult::new(0.06, 3.0, 0.5),
        }
    }
}

impl GaussFitter {
    /// A fitter carrying the source's default initial guess.
    pub fn new() -> Self {
        Self::default()
    }

    /// Replace the initial guess the next `fit` starts from.
    pub fn set_initial_parameters(&mut self, result: GaussFitResult) {
        self.initial = result;
    }

    /// The initial guess the next `fit` will start from.
    ///
    /// The source has no such accessor; its `init_param_` is protected.
    pub fn initial_parameters(&self) -> GaussFitResult {
        self.initial
    }

    /// Fit the Gaussian to `points`, each an `(x, y)` observation.
    ///
    /// The residual of observation `i` is
    /// `A exp(-(x_i - x0)^2 / (2 sigma^2)) - y_i` and the Jacobian is the
    /// analytic one the source supplies. The reported width is `|sigma|`: the
    /// optimizer is unconstrained and can settle on a negative width, which
    /// describes the same curve, and the source applies the same `fabs` with
    /// the comment that the absolute value is the correct solution. One of the
    /// two class-test cases exercises exactly that.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when
    ///
    /// - fewer than three points are supplied, or more than
    ///   [`MAX_POINTS`](crate::math::fitters::levenberg_marquardt::MAX_POINTS);
    /// - any coordinate is not finite;
    /// - the initial guess is not finite or its `sigma` is zero, which would
    ///   make the residual divide by zero;
    /// - the optimizer reports `ImproperInputParameters` or
    ///   `TooManyFunctionEvaluation`, the two states the source treats as
    ///   failure, throwing `Exception::UnableToFit`. Every other termination
    ///   state is accepted, as the source's comment records;
    /// - the fitted parameters are not finite. The source returns them as they
    ///   are.
    ///
    /// The source has no bound on the number of points; the ceiling is checked
    /// before anything is allocated, so a rejected call touches nothing.
    pub fn fit(&self, points: &[(f64, f64)]) -> Result<GaussFitResult> {
        preflight_points(points.len(), PARAMETERS)?;
        if points.len() < PARAMETERS {
            return Err(bad(
                "Gaussian fitting needs at least three points, one per parameter",
            ));
        }
        for &(x, y) in points {
            if !x.is_finite() || !y.is_finite() {
                return Err(bad("Gaussian fitting needs finite point coordinates"));
            }
        }
        let start = self.initial;
        if !start.a.is_finite() || !start.x0.is_finite() || !start.sigma.is_finite() {
            return Err(bad("Gaussian initial parameters must be finite"));
        }
        if start.sigma == 0.0 {
            return Err(bad("Gaussian initial sigma must not be zero"));
        }

        let mut x = [start.a, start.x0, start.sigma];
        let status = minimize(
            &mut x,
            points.len(),
            |p: &[f64], fvec: &mut [f64]| {
                let (amplitude, center, sigma) = (p[0], p[1], p[2]);
                let sig2 = 2.0 * sigma * sigma;
                for (slot, &(px, py)) in fvec.iter_mut().zip(points) {
                    *slot = amplitude * (-(px - center) * (px - center) / sig2).exp() - py;
                }
            },
            |p: &[f64], jac: &mut DenseMatrix| {
                let (amplitude, center, sigma) = (p[0], p[1], p[2]);
                let sig2 = 2.0 * sigma * sigma;
                let sig3 = 2.0 * sig2 * sigma;
                for (row, &(px, _)) in points.iter().enumerate() {
                    let xd = px - center;
                    let xd2 = xd * xd;
                    let j0 = (-xd2 / sig2).exp();
                    jac.set(row, 0, j0);
                    jac.set(
                        row,
                        1,
                        amplitude * j0 * (-(-2.0 * px + 2.0 * center) / sig2),
                    );
                    jac.set(row, 2, amplitude * j0 * (xd2 / sig3));
                }
                0
            },
            &LmParameters::default(),
        );
        if status == LmStatus::ImproperInputParameters
            || status == LmStatus::TooManyFunctionEvaluation
        {
            return Err(bad(&format!(
                "UnableToFit-GaussFitter: could not fit the Gaussian to the data: Error {}",
                status.code()
            )));
        }
        let result = GaussFitResult::new(x[0], x[1], x[2].abs());
        if !result.a.is_finite() || !result.x0.is_finite() || !result.sigma.is_finite() {
            return Err(bad(
                "UnableToFit-GaussFitter: the fit produced non-finite parameters",
            ));
        }
        Ok(result)
    }

    /// Evaluate `model` at every point of `evaluation_points`.
    ///
    /// Equivalent to calling
    /// [`GaussFitResult::eval`](crate::math::fitters::gauss::GaussFitResult::eval)
    /// on each, and the source implements it that way too, hoisting the
    /// normalization factor out of the loop. Any parameters may be used: the
    /// initial guess to see a before-fit state, or the fitted ones.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] under the same conditions as
    /// [`GaussFitResult::eval`](crate::math::fitters::gauss::GaussFitResult::eval),
    /// and when more than
    /// [`MAX_POINTS`](crate::math::fitters::levenberg_marquardt::MAX_POINTS)
    /// points are supplied. The source reserves the output vector and never
    /// fails.
    pub fn eval(evaluation_points: &[f64], model: &GaussFitResult) -> Result<Vec<f64>> {
        preflight_points(evaluation_points.len(), 1)?;
        model.check_density()?;
        let normalization = model.a / normal_pdf(model.x0, model.x0, model.sigma);
        let mut out = Vec::with_capacity(evaluation_points.len());
        for &point in evaluation_points {
            if !point.is_finite() {
                return Err(bad("Gaussian evaluation point must be finite"));
            }
            out.push(normal_pdf(point, model.x0, model.sigma) * normalization);
        }
        Ok(out)
    }
}

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn half_log_two_pi_matches_the_source_expression() {
        assert_eq!(HALF_LOG_TWO_PI, 0.5 * (2.0 * std::f64::consts::PI).ln());
    }

    /// Boost's `pdf` divides by `sd * sqrt(2 * constants::pi<RealType>())`, so
    /// the divisor is the square root of the rounded `2 * pi`, evaluated at run
    /// time. It is *not* `constants::root_two_pi`, whose decimal literal rounds
    /// to a different `f64` one unit in the last place higher. Both halves are
    /// asserted so that a future edit cannot quietly swap one for the other.
    #[test]
    fn the_density_divisor_is_the_square_root_boost_evaluates() {
        assert_eq!(SQRT_TWO_PI, (2.0 * std::f64::consts::PI).sqrt());
        let root_two_pi_literal = 2.506_628_274_631_000_7_f64;
        assert_ne!(SQRT_TWO_PI, root_two_pi_literal);
        assert_eq!(
            f64::from_bits(SQRT_TWO_PI.to_bits() + 1),
            root_two_pi_literal
        );
    }

    #[test]
    fn the_default_result_is_the_sources_invalid_marker() {
        let result = GaussFitResult::default();
        assert_eq!(result, GaussFitResult::new(-1.0, -1.0, -1.0));
        assert!(result.eval(0.0).is_err());
    }

    #[test]
    fn the_default_initial_guess_is_the_sources() {
        assert_eq!(
            GaussFitter::new().initial_parameters(),
            GaussFitResult::new(0.06, 3.0, 0.5)
        );
    }

    #[test]
    fn eval_reaches_the_amplitude_at_the_center() {
        // Independent of any transcribed literal: the model is A at x0 and
        // A * exp(-1/2) one standard deviation away.
        let model = GaussFitResult::new(7.0, 2.0, 0.5);
        assert!((model.eval(2.0).unwrap() - 7.0).abs() < 1e-12);
        let one_sigma = model.eval(2.5).unwrap();
        assert!((one_sigma - 7.0 * (-0.5f64).exp()).abs() < 1e-12);
        assert!((model.eval(1.5).unwrap() - one_sigma).abs() < 1e-15);
    }

    #[test]
    fn a_degenerate_width_is_refused_everywhere() {
        let model = GaussFitResult::new(1.0, 0.0, 0.0);
        assert!(model.eval(0.0).is_err());
        assert!(model.log_eval_no_normalize(0.0).is_err());
        assert!(GaussFitter::eval(&[0.0], &model).is_err());
        let negative = GaussFitResult::new(1.0, 0.0, -1.0);
        assert!(negative.log_eval_no_normalize(0.0).is_err());
    }

    #[test]
    fn fitting_refuses_input_it_cannot_use() {
        let fitter = GaussFitter::new();
        assert!(fitter.fit(&[]).is_err());
        assert!(fitter.fit(&[(0.0, 1.0), (1.0, 2.0)]).is_err());
        assert!(
            fitter
                .fit(&[(0.0, 1.0), (1.0, f64::NAN), (2.0, 3.0)])
                .is_err()
        );
        let mut zero_sigma = GaussFitter::new();
        zero_sigma.set_initial_parameters(GaussFitResult::new(1.0, 0.0, 0.0));
        assert!(
            zero_sigma
                .fit(&[(0.0, 1.0), (1.0, 2.0), (2.0, 3.0)])
                .is_err()
        );
    }
}
