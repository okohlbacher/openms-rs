// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Julianus Pfeuffer, OpenMS Rust contributors $

//! Weighted maximum-likelihood fit of a Gumbel distribution to raw samples.
//!
//! Ported from `MATH/STATISTICS/GumbelMaxLikelihoodFitter.h` and its `.cpp` at
//! revision `bc9cc12`. Unlike
//! [`gumbel`](crate::math::fitters::gumbel), which fits a density curve to
//! (x, y) points by least squares, this class takes raw observations with
//! weights and drives the weighted negative log-likelihood
//! `-sum_i w_i (-ln b - z_i - exp(-z_i))`, with `z_i = (x_i - a) / |b|`,
//! towards zero with the same Levenberg-Marquardt solver.
//!
//! Towards zero, not towards its minimum: the source hands the solver a
//! two-element residual vector whose first element is that scalar and whose
//! second is a constant zero, and Levenberg-Marquardt minimizes the sum of
//! squares of what it is given. This is a quirk of the source, reproduced
//! here, and is why the class test's expectations are loose. See
//! [`GumbelMaxLikelihoodFitter::fit_weighted`](crate::math::fitters::gumbel_max_likelihood::GumbelMaxLikelihoodFitter::fit_weighted)
//! and `docs/DISTRIBUTION_FITTERS_SUPPORT.md`.

use super::levenberg_marquardt::{
    DenseMatrix, LmParameters, LmStatus, minimize, numerical_jacobian, preflight_points,
};
use crate::{Error, Result};

/// Number of residuals the source's functor declares.
///
/// The likelihood is one scalar, but Eigen's solver rejects a problem with
/// fewer residuals than parameters, and there are two of those - location and
/// scale - so the source declares two residuals and leaves the second
/// permanently zero.
const RESIDUALS: usize = 2;

/// Location and scale of a Gumbel distribution, as this fitter reports them.
///
/// This is a distinct type from
/// [`gumbel::GumbelDistributionFitResult`](crate::math::fitters::gumbel::GumbelDistributionFitResult):
/// the C++ declares two same-named structs nested in two different classes,
/// and this one has no default constructor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GumbelDistributionFitResult {
    /// Location parameter `a`, the mode of the density.
    pub a: f64,
    /// Scale parameter `b`.
    pub b: f64,
}

impl GumbelDistributionFitResult {
    /// A result with the given location and scale.
    ///
    /// Nothing is validated, matching the source's only constructor.
    pub fn new(a: f64, b: f64) -> Self {
        Self { a, b }
    }

    /// Evaluate the log density at `x`.
    ///
    /// Computes `-ln(b) - (x - a)/b - exp(-(x - a)/b)`, character for character
    /// the same body as
    /// [`gumbel::GumbelDistributionFitResult::log_eval_no_normalize`](crate::math::fitters::gumbel::GumbelDistributionFitResult::log_eval_no_normalize);
    /// the source duplicates it across the two classes.
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
}

/// Fits a Gumbel distribution to weighted samples by maximum likelihood.
///
/// The source's default initial guess is `(a, b) = (0.25, 0.1)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GumbelMaxLikelihoodFitter {
    initial: GumbelDistributionFitResult,
}

impl Default for GumbelMaxLikelihoodFitter {
    /// The source's default constructor guess `(0.25, 0.1)`.
    fn default() -> Self {
        Self {
            initial: GumbelDistributionFitResult::new(0.25, 0.1),
        }
    }
}

impl GumbelMaxLikelihoodFitter {
    /// A fitter carrying the source's default initial guess `(0.25, 0.1)`.
    pub fn new() -> Self {
        Self::default()
    }

    /// A fitter starting from `initial`, the source's second constructor.
    pub fn with_initial_parameters(initial: GumbelDistributionFitResult) -> Self {
        Self { initial }
    }

    /// Replace the start parameters `a` and `b` the next fit uses.
    pub fn set_initial_parameters(&mut self, result: GumbelDistributionFitResult) {
        self.initial = result;
    }

    /// The parameters the next fit will start from.
    ///
    /// After a successful
    /// [`fit_weighted`](crate::math::fitters::gumbel_max_likelihood::GumbelMaxLikelihoodFitter::fit_weighted)
    /// these are the fitted parameters: the source overwrites `init_param_`
    /// with the result, so a second fit continues from where the first stopped.
    /// The source has no such accessor; its member is protected.
    pub fn initial_parameters(&self) -> GumbelDistributionFitResult {
        self.initial
    }

    /// Fit the Gumbel distribution to the weighted samples `x` with weights
    /// `w`, by driving the weighted negative log-likelihood towards zero.
    ///
    /// The objective is the source's, in the source's accumulation order: a
    /// running sum over the samples in index order of
    /// `w_i * (-ln|b| - z_i - exp(-z_i))` with `z_i = (x_i - a) / |b|`,
    /// negated once at the end. The scale enters through its absolute value, so
    /// the search is unconstrained, and the reported scale is `|b|`, as the
    /// source reports it.
    ///
    /// The Jacobian is the forward-difference approximation of
    /// `Eigen::NumericalDiff` at its default step, not an analytic derivative,
    /// because that is what the source uses and the step size changes where the
    /// iteration stops.
    ///
    /// Empty input is not an error: with no samples the objective is
    /// identically zero, its gradient vanishes, the solver stops immediately
    /// with `CosinusTooSmall`, and the initial parameters come back unchanged.
    /// Three of the class-test sections rely on exactly that to read back the
    /// initial guess.
    ///
    /// On success the fitter's own start parameters are replaced by the
    /// result, which is why this takes `&mut self`; the source does the same in
    /// `fitWeighted`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when
    ///
    /// - `x` and `w` have different lengths. The source iterates the weights
    ///   with a second iterator advanced in lockstep and reads past the end of
    ///   the shorter one;
    /// - more than
    ///   [`MAX_POINTS`](crate::math::fitters::levenberg_marquardt::MAX_POINTS)
    ///   samples are supplied;
    /// - any sample or weight is not finite;
    /// - the initial guess is not finite or its scale is zero, which would make
    ///   the objective divide by zero;
    /// - the optimizer reports `ImproperInputParameters`, the only reachable
    ///   member of the status range the source rejects with
    ///   `Exception::UnableToFit`;
    /// - the fitted parameters are not finite.
    ///
    /// A rejected call leaves the fitter's start parameters untouched.
    pub fn fit_weighted(&mut self, x: &[f64], w: &[f64]) -> Result<GumbelDistributionFitResult> {
        if x.len() != w.len() {
            return Err(bad(
                "Gumbel maximum-likelihood fitting needs one weight per sample",
            ));
        }
        preflight_points(x.len(), 1)?;
        for (&sample, &weight) in x.iter().zip(w) {
            if !sample.is_finite() || !weight.is_finite() {
                return Err(bad(
                    "Gumbel maximum-likelihood fitting needs finite samples and weights",
                ));
            }
        }
        let start = self.initial;
        if !start.a.is_finite() || !start.b.is_finite() {
            return Err(bad("Gumbel initial parameters must be finite"));
        }
        if start.b == 0.0 {
            return Err(bad("Gumbel initial scale must not be zero"));
        }

        let objective = |v: &[f64], fvec: &mut [f64]| {
            let sigma = v[1].abs();
            let logsigma = sigma.ln();
            let mut sum = 0.0f64;
            for (&sample, &weight) in x.iter().zip(w) {
                let diff = (sample - v[0]) / sigma;
                sum += weight * (-logsigma - diff - (-diff).exp());
            }
            fvec[0] = -sum;
            fvec[1] = 0.0;
        };

        let mut parameters = [start.a, start.b];
        let status = minimize(
            &mut parameters,
            RESIDUALS,
            &objective,
            |v: &[f64], jac: &mut DenseMatrix| numerical_jacobian(v, RESIDUALS, jac, &objective),
            &LmParameters::default(),
        );
        if status == LmStatus::ImproperInputParameters {
            return Err(bad(
                "UnableToFit-GumbelMaxLikelihoodFitter: could not fit the gumbel distribution to the data",
            ));
        }
        let result = GumbelDistributionFitResult::new(parameters[0], parameters[1].abs());
        if !result.a.is_finite() || !result.b.is_finite() {
            return Err(bad(
                "UnableToFit-GumbelMaxLikelihoodFitter: the fit produced non-finite parameters",
            ));
        }
        self.initial = result;
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
    fn the_default_initial_guess_is_the_sources() {
        assert_eq!(
            GumbelMaxLikelihoodFitter::new().initial_parameters(),
            GumbelDistributionFitResult::new(0.25, 0.1)
        );
    }

    #[test]
    fn the_log_density_is_maximal_at_the_location_parameter() {
        // The mode of a Gumbel is a and the log density there is -ln b - 1.
        let model = GumbelDistributionFitResult::new(2.0, 1.0);
        assert!((model.log_eval_no_normalize(2.0).unwrap() + 1.0).abs() < 1e-14);
        assert!(
            model.log_eval_no_normalize(2.0).unwrap() > model.log_eval_no_normalize(3.0).unwrap()
        );
        assert!(
            model.log_eval_no_normalize(2.0).unwrap() > model.log_eval_no_normalize(1.0).unwrap()
        );
        let wider = GumbelDistributionFitResult::new(2.0, 2.0);
        assert!((wider.log_eval_no_normalize(2.0).unwrap() + 2.0f64.ln() + 1.0).abs() < 1e-14);
    }

    #[test]
    fn mismatched_weights_are_refused() {
        let mut fitter = GumbelMaxLikelihoodFitter::new();
        assert!(fitter.fit_weighted(&[1.0, 2.0], &[1.0]).is_err());
        assert!(fitter.fit_weighted(&[f64::NAN], &[1.0]).is_err());
        // A rejected call leaves the start parameters alone.
        assert_eq!(
            fitter.initial_parameters(),
            GumbelDistributionFitResult::new(0.25, 0.1)
        );
    }

    #[test]
    fn empty_input_returns_the_initial_parameters() {
        let mut fitter = GumbelMaxLikelihoodFitter::with_initial_parameters(
            GumbelDistributionFitResult::new(3.0, 0.5),
        );
        let result = fitter.fit_weighted(&[], &[]).unwrap();
        assert_eq!(result, GumbelDistributionFitResult::new(3.0, 0.5));
    }
}
