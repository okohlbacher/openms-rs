// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: David Wojnar, OpenMS Rust contributors $

//! Least-squares fit of a Gumbel distribution to a set of (x, y) points.
//!
//! Ported from `MATH/STATISTICS/GumbelDistributionFitter.h` and its `.cpp` at
//! revision `bc9cc12`. The model is
//! `(1 / b) exp((a - x) / b) exp(-exp((a - x) / b))`, with `a` the location and
//! `b` the scale.
//!
//! For the maximum-likelihood fitter of the same distribution, which is a
//! separate class with a separate result type, see
//! [`gumbel_max_likelihood`](crate::math::fitters::gumbel_max_likelihood).
//!
//! Two members declared by the C++ header have no definition anywhere in the
//! SDK and are therefore not ported: `GumbelDistributionFitResult::eval` and
//! `GumbelDistributionFitter::fitWeighted`. Calling either from C++ fails to
//! link. `docs/DISTRIBUTION_FITTERS_SUPPORT.md` records both.

use super::levenberg_marquardt::{DenseMatrix, LmParameters, LmStatus, minimize, preflight_points};
use crate::{Error, Result};

/// Number of fitted parameters: location and scale.
const PARAMETERS: usize = 2;

/// Location and scale of a Gumbel distribution.
///
/// Doubles as the initial guess and as the fit result, exactly like the
/// source's `GumbelDistributionFitter::GumbelDistributionFitResult`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GumbelDistributionFitResult {
    /// Location parameter `a`, the mode of the density.
    pub a: f64,
    /// Scale parameter `b`.
    pub b: f64,
}

impl Default for GumbelDistributionFitResult {
    /// `(a, b) = (1.0, 2.0)`, the source's default arguments. The class test
    /// asserts exactly these two values for a default-constructed result.
    fn default() -> Self {
        Self { a: 1.0, b: 2.0 }
    }
}

impl GumbelDistributionFitResult {
    /// A result with the given location and scale.
    ///
    /// Nothing is validated, matching the source's constructor.
    pub fn new(a: f64, b: f64) -> Self {
        Self { a, b }
    }

    /// Evaluate the log density at `x`.
    ///
    /// Computes `-ln(b) - (x - a)/b - exp(-(x - a)/b)`, the source's
    /// expression and order. The "no normalize" in the name refers to the
    /// absent amplitude: this is a proper log density, which is what a
    /// likelihood needs.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `b` is not finite and strictly
    /// positive, or when `a` or `x` is not finite. The source takes `log` of a
    /// non-positive `b` and returns NaN without comment.
    pub fn log_eval_no_normalize(&self, x: f64) -> Result<f64> {
        if !self.b.is_finite() || self.b <= 0.0 {
            return Err(bad("Gumbel scale must be finite and positive"));
        }
        if !self.a.is_finite() || !x.is_finite() {
            return Err(bad("Gumbel location and evaluation point must be finite"));
        }
        let diff = (x - self.a) / self.b;
        Ok(-self.b.ln() - diff - (-diff).exp())
    }

    /// Evaluate the density at `x`.
    ///
    /// Computes `exp((a - x)/b) exp(-exp((a - x)/b)) / b`, the expression the
    /// source's optimizer functor uses for its residual. The source *declares*
    /// a method called `eval` on this struct but never defines it, so no C++
    /// caller can use one; this is the residual model made reachable, not a
    /// port of that declaration.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `b` is not finite and strictly
    /// positive, when `a` or `x` is not finite, or when the result is not
    /// finite.
    pub fn eval(&self, x: f64) -> Result<f64> {
        if !self.b.is_finite() || self.b <= 0.0 {
            return Err(bad("Gumbel scale must be finite and positive"));
        }
        if !self.a.is_finite() || !x.is_finite() {
            return Err(bad("Gumbel location and evaluation point must be finite"));
        }
        let value = gumbel_density(self.a, self.b, x);
        if !value.is_finite() {
            return Err(bad("Gumbel density is not finite at this point"));
        }
        Ok(value)
    }
}

/// The source functor's residual model, in its arithmetic order.
fn gumbel_density(a: f64, b: f64, x: f64) -> f64 {
    let z = ((a - x) / b).exp();
    (z * (-z).exp()) / b
}

/// Fits a Gumbel distribution to a set of points by non-linear least squares.
///
/// The source's default initial guess is `(a, b) = (0.25, 0.1)`, carried over
/// verbatim because the residual surface is not convex. Note that the fitter's
/// default guess is not the result type's default `(1.0, 2.0)`: the C++
/// constructor assigns `GumbelDistributionFitResult(0.25, 0.1)` over the
/// default-constructed member.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GumbelDistributionFitter {
    initial: GumbelDistributionFitResult,
}

impl Default for GumbelDistributionFitter {
    /// The source's constructor guess `(0.25, 0.1)`.
    fn default() -> Self {
        Self {
            initial: GumbelDistributionFitResult::new(0.25, 0.1),
        }
    }
}

impl GumbelDistributionFitter {
    /// A fitter carrying the source's default initial guess.
    pub fn new() -> Self {
        Self::default()
    }

    /// Replace the start parameters `a` and `b` the next `fit` uses.
    pub fn set_initial_parameters(&mut self, result: GumbelDistributionFitResult) {
        self.initial = result;
    }

    /// The initial guess the next `fit` will start from.
    ///
    /// The source has no such accessor; its `init_param_` is protected.
    pub fn initial_parameters(&self) -> GumbelDistributionFitResult {
        self.initial
    }

    /// Fit the Gumbel distribution to `points`, each an `(x, y)` observation.
    ///
    /// The residual of observation `i` is
    /// `exp((a - x_i)/b) exp(-exp((a - x_i)/b)) / b - y_i` and the Jacobian is
    /// the source's analytic one. Unlike
    /// [`GaussFitter::fit`](crate::math::fitters::gauss::GaussFitter::fit) the
    /// source applies no `fabs` here, so a fit that settles on a negative scale
    /// reports it; that parameter set is not a density, and the caller sees it
    /// rather than a silently corrected one.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when
    ///
    /// - fewer than two points are supplied, or more than
    ///   [`MAX_POINTS`](crate::math::fitters::levenberg_marquardt::MAX_POINTS);
    /// - any coordinate is not finite;
    /// - the initial guess is not finite or its scale is zero, which would make
    ///   the residual divide by zero;
    /// - the optimizer reports `ImproperInputParameters`, which is the only
    ///   reachable member of the status range the source rejects;
    /// - the fitted parameters are not finite.
    pub fn fit(&self, points: &[(f64, f64)]) -> Result<GumbelDistributionFitResult> {
        preflight_points(points.len(), PARAMETERS)?;
        if points.len() < PARAMETERS {
            return Err(bad(
                "Gumbel fitting needs at least two points, one per parameter",
            ));
        }
        for &(x, y) in points {
            if !x.is_finite() || !y.is_finite() {
                return Err(bad("Gumbel fitting needs finite point coordinates"));
            }
        }
        let start = self.initial;
        if !start.a.is_finite() || !start.b.is_finite() {
            return Err(bad("Gumbel initial parameters must be finite"));
        }
        if start.b == 0.0 {
            return Err(bad("Gumbel initial scale must not be zero"));
        }

        let mut x = [start.a, start.b];
        let status = minimize(
            &mut x,
            points.len(),
            |v: &[f64], fvec: &mut [f64]| {
                let (a, b) = (v[0], v[1]);
                for (slot, &(px, py)) in fvec.iter_mut().zip(points) {
                    let z = ((a - px) / b).exp();
                    *slot = (z * (-z).exp()) / b - py;
                }
            },
            |v: &[f64], jac: &mut DenseMatrix| {
                let (a, b) = (v[0], v[1]);
                for (row, &(px, _)) in points.iter().enumerate() {
                    let z = ((a - px) / b).exp();
                    let f = z * (-z).exp();
                    let part_dev_a = (f - z * z * (-z).exp()) / (b * b);
                    jac.set(row, 0, part_dev_a);
                    let dev_z = (px - a) / (b * b);
                    let cum = f * dev_z;
                    let part_dev_b = ((cum - z * cum) * b - f) / (b * b);
                    jac.set(row, 1, part_dev_b);
                }
                0
            },
            &LmParameters::default(),
        );
        if status == LmStatus::ImproperInputParameters {
            return Err(bad(
                "UnableToFit-GumbelDistributionFitter: could not fit the gumbel distribution to the data",
            ));
        }
        let result = GumbelDistributionFitResult::new(x[0], x[1]);
        if !result.a.is_finite() || !result.b.is_finite() {
            return Err(bad(
                "UnableToFit-GumbelDistributionFitter: the fit produced non-finite parameters",
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
    fn the_default_result_is_the_sources_default_arguments() {
        assert_eq!(
            GumbelDistributionFitResult::default(),
            GumbelDistributionFitResult::new(1.0, 2.0)
        );
    }

    #[test]
    fn the_fitters_default_guess_differs_from_the_results_default() {
        assert_eq!(
            GumbelDistributionFitter::new().initial_parameters(),
            GumbelDistributionFitResult::new(0.25, 0.1)
        );
        assert_ne!(
            GumbelDistributionFitter::new().initial_parameters(),
            GumbelDistributionFitResult::default()
        );
    }

    #[test]
    fn the_density_peaks_at_the_location_parameter() {
        // Independent invariants: the mode of a Gumbel is a, the density there
        // is 1/(b e), and the log density is the log of the density.
        let model = GumbelDistributionFitResult::new(2.0, 0.5);
        let peak = model.eval(2.0).unwrap();
        assert!((peak - 1.0 / (0.5 * std::f64::consts::E)).abs() < 1e-14);
        assert!(model.eval(1.5).unwrap() < peak);
        assert!(model.eval(2.5).unwrap() < peak);
        assert!(
            (model.log_eval_no_normalize(2.3).unwrap() - model.eval(2.3).unwrap().ln()).abs()
                < 1e-13
        );
    }

    #[test]
    fn a_degenerate_scale_is_refused() {
        assert!(
            GumbelDistributionFitResult::new(0.0, 0.0)
                .eval(0.0)
                .is_err()
        );
        assert!(
            GumbelDistributionFitResult::new(0.0, -1.0)
                .log_eval_no_normalize(0.0)
                .is_err()
        );
    }

    #[test]
    fn fitting_refuses_input_it_cannot_use() {
        let fitter = GumbelDistributionFitter::new();
        assert!(fitter.fit(&[]).is_err());
        assert!(fitter.fit(&[(1.0, 1.0)]).is_err());
        assert!(fitter.fit(&[(1.0, 1.0), (2.0, f64::NAN)]).is_err());
        let mut zero_scale = GumbelDistributionFitter::new();
        zero_scale.set_initial_parameters(GumbelDistributionFitResult::new(1.0, 0.0));
        assert!(zero_scale.fit(&[(1.0, 1.0), (2.0, 1.0)]).is_err());
    }
}
