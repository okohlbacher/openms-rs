// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Andreas Bertsch, OpenMS Rust contributors $

//! Least-squares fit of a Gamma distribution to a set of (x, y) points.
//!
//! Ported from `MATH/STATISTICS/GammaDistributionFitter.h` and its `.cpp` at
//! revision `bc9cc12`. The model is the two-parameter density
//! `b^p / Gamma(p) * x^(p - 1) * exp(-b x)`.
//!
//! The source carries a `@note` that matters: the fitted function is a
//! *customized* Gamma density which is defined to be zero whenever `b` or `p`
//! drops to zero or below, so that an unconstrained optimizer can be used on a
//! distribution that is only defined for positive parameters. That branch, and
//! the matching all-zero Jacobian, are reproduced here - see
//! [`GammaDistributionFitter::fit`](crate::math::fitters::gamma::GammaDistributionFitter::fit).
//!
//! See `docs/DISTRIBUTION_FITTERS_SUPPORT.md`.

use super::levenberg_marquardt::{DenseMatrix, LmParameters, LmStatus, minimize, preflight_points};
use crate::{Error, Result};

/// Number of fitted parameters: rate and shape.
const PARAMETERS: usize = 2;

/// Rate and shape of a Gamma distribution.
///
/// Doubles as the initial guess and as the fit result, exactly like the
/// source's `GammaDistributionFitter::GammaDistributionFitResult`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GammaDistributionFitResult {
    /// Parameter `b`, the rate. The density is zero for `b <= 0`.
    pub b: f64,
    /// Parameter `p`, the shape. The density is zero for `p <= 0`.
    pub p: f64,
}

impl GammaDistributionFitResult {
    /// A result with the given rate and shape.
    ///
    /// Nothing is validated, matching the source's only constructor. The
    /// source has no default constructor for this struct, so neither does this
    /// type.
    pub fn new(b: f64, p: f64) -> Self {
        Self { b, p }
    }

    /// Evaluate the customized Gamma density at `x`.
    ///
    /// Returns `b^p / Gamma(p) * x^(p - 1) * exp(-b x)`, or zero when `b <= 0`
    /// or `p <= 0`, in the source's arithmetic order. The source has no such
    /// public method - the expression lives inside the optimizer functor - but
    /// exposing it lets a caller see the model a fit produced without
    /// rebuilding the formula.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `b`, `p` or `x` is not finite, or
    /// when the expression is not finite - for example when `x` is zero and
    /// `p < 1`, where `x^(p - 1)` is infinite.
    pub fn eval(&self, x: f64) -> Result<f64> {
        if !self.b.is_finite() || !self.p.is_finite() || !x.is_finite() {
            return Err(bad("Gamma parameters and evaluation point must be finite"));
        }
        let value = gamma_density(self.b, self.p, x);
        if !value.is_finite() {
            return Err(bad("Gamma density is not finite at this point"));
        }
        Ok(value)
    }
}

/// The source's residual model: zero outside the positive quadrant.
fn gamma_density(b: f64, p: f64, x: f64) -> f64 {
    if b > 0.0 && p > 0.0 {
        b.powf(p) / libm::tgamma(p) * x.powf(p - 1.0) * (-b * x).exp()
    } else {
        0.0
    }
}

/// The digamma function `psi(x) = d/dx ln Gamma(x)`, for `x > 0`.
///
/// The source calls `boost::math::digamma`. Boost uses rational minimax
/// approximations; this is the recurrence `psi(x) = psi(x + 1) - 1/x` up to
/// `x >= 10` followed by the standard asymptotic series through the Bernoulli
/// number `B14`. Both are accurate to a few units in the last place.
///
/// That substitution is not free of consequences and is not claimed to be. The
/// value enters the shape column of the Jacobian, and in Levenberg-Marquardt
/// the Jacobian sets the `diag` scaling, the trust-region radius, the gradient
/// test and every termination test, so a perturbed Jacobian can stop the
/// iteration at a slightly different point - the same argument the solver
/// module makes for keeping Eigen's stopping rule. What makes the substitution
/// safe here is the size of the slack, not an absence of effect: the agreement
/// with the closed forms is better than `1e-14` absolute, while the Gamma class
/// test asserts only the parameters its data were generated from, `b = 7.25`
/// and `p = 3.11`, at `0.01` absolute, and the port lands `2.7e-3` and `6.9e-3`
/// away. Nothing in this group publishes a C++-produced Gamma parameter, so no
/// tighter statement is available and none is made.
fn digamma(x: f64) -> f64 {
    let mut value = x;
    let mut result = 0.0f64;
    while value < 10.0 {
        result -= 1.0 / value;
        value += 1.0;
    }
    let f = 1.0 / (value * value);
    let series = f
        * (1.0 / 12.0
            - f * (1.0 / 120.0
                - f * (1.0 / 252.0
                    - f * (1.0 / 240.0
                        - f * (1.0 / 132.0 - f * (691.0 / 32760.0 - f * (1.0 / 12.0)))))));
    result + value.ln() - 0.5 / value - series
}

/// Fits a Gamma distribution to a set of points by non-linear least squares.
///
/// The source's default initial guess is `(b, p) = (1.0, 5.0)`, carried over
/// verbatim because the residual surface is not convex.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GammaDistributionFitter {
    initial: GammaDistributionFitResult,
}

impl Default for GammaDistributionFitter {
    /// The source's constructor guess `(1.0, 5.0)`.
    fn default() -> Self {
        Self {
            initial: GammaDistributionFitResult::new(1.0, 5.0),
        }
    }
}

impl GammaDistributionFitter {
    /// A fitter carrying the source's default initial guess.
    pub fn new() -> Self {
        Self::default()
    }

    /// Replace the start parameters `b` and `p` the next `fit` uses.
    pub fn set_initial_parameters(&mut self, result: GammaDistributionFitResult) {
        self.initial = result;
    }

    /// The initial guess the next `fit` will start from.
    ///
    /// The source has no such accessor; its `init_param_` is protected.
    pub fn initial_parameters(&self) -> GammaDistributionFitResult {
        self.initial
    }

    /// Fit the Gamma distribution to `points`, each an `(x, y)` observation.
    ///
    /// The residual of observation `i` is
    /// `b^p / Gamma(p) * x_i^(p - 1) * exp(-b x_i) - y_i` while both parameters
    /// are positive, and `-y_i` otherwise - the source's way of keeping the
    /// search unconstrained without ever evaluating an undefined density. The
    /// Jacobian is the source's analytic one, and is set to zero in the same
    /// non-positive branch, which makes the optimizer take no step away from a
    /// non-positive parameter except through the trust region.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when
    ///
    /// - fewer than two points are supplied, or more than
    ///   [`MAX_POINTS`](crate::math::fitters::levenberg_marquardt::MAX_POINTS);
    /// - any coordinate is not finite, or any `x` is negative - the source
    ///   raises `x` to a real power, which is NaN for a negative base;
    /// - the initial guess is not finite;
    /// - the optimizer reports `ImproperInputParameters`. The source throws
    ///   `Exception::UnableToFit` for any status less than or equal to that
    ///   one; the smaller values, `NotStarted` and `Running`, cannot be
    ///   returned by a completed `minimize`;
    /// - the fitted parameters are not finite.
    pub fn fit(&self, points: &[(f64, f64)]) -> Result<GammaDistributionFitResult> {
        preflight_points(points.len(), PARAMETERS)?;
        if points.len() < PARAMETERS {
            return Err(bad(
                "Gamma fitting needs at least two points, one per parameter",
            ));
        }
        for &(x, y) in points {
            if !x.is_finite() || !y.is_finite() {
                return Err(bad("Gamma fitting needs finite point coordinates"));
            }
            if x < 0.0 {
                return Err(bad(
                    "Gamma fitting needs non-negative abscissae; the density is undefined below zero",
                ));
            }
        }
        let start = self.initial;
        if !start.b.is_finite() || !start.p.is_finite() {
            return Err(bad("Gamma initial parameters must be finite"));
        }

        let mut x = [start.b, start.p];
        let status = minimize(
            &mut x,
            points.len(),
            |v: &[f64], fvec: &mut [f64]| {
                let (b, p) = (v[0], v[1]);
                if b > 0.0 && p > 0.0 {
                    for (slot, &(px, py)) in fvec.iter_mut().zip(points) {
                        *slot =
                            b.powf(p) / libm::tgamma(p) * px.powf(p - 1.0) * (-b * px).exp() - py;
                    }
                } else {
                    for (slot, &(_, py)) in fvec.iter_mut().zip(points) {
                        *slot = -py;
                    }
                }
            },
            |v: &[f64], jac: &mut DenseMatrix| {
                let (b, p) = (v[0], v[1]);
                if b > 0.0 && p > 0.0 {
                    for (row, &(px, _)) in points.iter().enumerate() {
                        let part_dev_b = px.powf(p - 1.0) * (-px * b).exp() / libm::tgamma(p)
                            * (p * b.powf(p - 1.0) - px * b.powf(p));
                        jac.set(row, 0, part_dev_b);
                        let factor = (-b * px).exp() * px.powf(p - 1.0) * b.powf(p)
                            / (libm::tgamma(p) * libm::tgamma(p));
                        let argument =
                            (b.ln() + px.ln()) * libm::tgamma(p) - libm::tgamma(p) * digamma(p);
                        jac.set(row, 1, factor * argument);
                    }
                } else {
                    for row in 0..points.len() {
                        jac.set(row, 0, 0.0);
                        jac.set(row, 1, 0.0);
                    }
                }
                0
            },
            &LmParameters::default(),
        );
        if status == LmStatus::ImproperInputParameters {
            return Err(bad(
                "UnableToFit-GammaDistributionFitter: could not fit the gamma distribution to the data",
            ));
        }
        let result = GammaDistributionFitResult::new(x[0], x[1]);
        if !result.b.is_finite() || !result.p.is_finite() {
            return Err(bad(
                "UnableToFit-GammaDistributionFitter: the fit produced non-finite parameters",
            ));
        }
        Ok(result)
    }
}

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn digamma_reproduces_its_closed_form_values() {
        // psi(1) = -gamma, psi(1/2) = -gamma - 2 ln 2, psi(n+1) = psi(n) + 1/n.
        let euler = 0.577_215_664_901_532_9_f64;
        assert!((digamma(1.0) + euler).abs() < 1e-14);
        assert!((digamma(0.5) + euler + 2.0 * 2.0f64.ln()).abs() < 1e-14);
        for n in 1..12 {
            let x = n as f64;
            assert!(
                (digamma(x + 1.0) - digamma(x) - 1.0 / x).abs() < 1e-13,
                "recurrence at {x}"
            );
        }
    }

    #[test]
    fn the_default_initial_guess_is_the_sources() {
        assert_eq!(
            GammaDistributionFitter::new().initial_parameters(),
            GammaDistributionFitResult::new(1.0, 5.0)
        );
    }

    #[test]
    fn the_customized_density_is_zero_outside_the_positive_quadrant() {
        assert_eq!(
            GammaDistributionFitResult::new(-1.0, 3.0)
                .eval(1.0)
                .unwrap(),
            0.0
        );
        assert_eq!(
            GammaDistributionFitResult::new(1.0, 0.0).eval(1.0).unwrap(),
            0.0
        );
        // Exponential special case p = 1: density is b exp(-b x).
        let model = GammaDistributionFitResult::new(2.0, 1.0);
        assert!((model.eval(0.0).unwrap() - 2.0).abs() < 1e-14);
        assert!((model.eval(1.0).unwrap() - 2.0 * (-2.0f64).exp()).abs() < 1e-14);
    }

    #[test]
    fn fitting_refuses_input_it_cannot_use() {
        let fitter = GammaDistributionFitter::new();
        assert!(fitter.fit(&[]).is_err());
        assert!(fitter.fit(&[(1.0, 1.0)]).is_err());
        assert!(fitter.fit(&[(-1.0, 1.0), (1.0, 1.0)]).is_err());
        assert!(fitter.fit(&[(1.0, f64::INFINITY), (2.0, 1.0)]).is_err());
    }
}
