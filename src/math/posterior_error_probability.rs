// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A two-component mixture model of search-engine scores, fitted by
//! expectation-maximisation, and the posterior error probabilities it yields.
//!
//! Port of
//! `src/openms/include/OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h`
//! and the numerical half of its translation unit. See
//! `docs/POSTERIOR_ERROR_PROBABILITY_SUPPORT.md`, which lists every public
//! member of the header and what became of it.
//!
//! Incorrectly assigned peptide-spectrum matches are modelled by one component
//! and correctly assigned ones by a Gaussian.
//! [`crate::math::posterior_error_probability::PosteriorErrorProbabilityModel::fit`]
//! fits **two Gaussians** whatever
//! [`crate::math::posterior_error_probability::IncorrectComponent`] says, and
//! uses the choice only to decide which density
//! [`crate::math::posterior_error_probability::PosteriorErrorProbabilityModel::compute_probability`]
//! and the gnuplot formulae evaluate;
//! [`crate::math::posterior_error_probability::PosteriorErrorProbabilityModel::fit_gumbel_gauss`]
//! fits a Gumbel by weighted maximum likelihood instead. The source marks the
//! first with a `TODO: incorrect is currently filled with gauss as fitting
//! gumble is not supported`, and that mismatch — a Gaussian fitted, a Gumbel
//! evaluated — is preserved because every probability the source has ever
//! reported depends on it.
//!
//! The component fitters are
//! [`crate::math::fitters::gauss::GaussFitResult`] and
//! [`crate::math::fitters::gumbel_max_likelihood::GumbelMaxLikelihoodFitter`].
//!
//! # Not ported here
//!
//! Five public members are unported. The header's identification-facing half
//! — `extractAndTransformScores` and `updateScores`, with their two *private*
//! helpers `transformScore_` and `getScore_`, which sit below the `private:`
//! at `PosteriorErrorProbabilityModel.h:239` and are not part of the public
//! surface — and its plotting half — `initPlots`,
//! `plotTargetDecoyEstimation`, `tryGnuplot` — need `METADATA`,
//! `FORMAT/TextFile` and `SYSTEM/File`, which are the crate's `metadata`,
//! `format` and `system` modules. This crate ratchets its cross-module
//! dependency graph and `math` reaches no other top-level module, so those
//! members belong above this layer rather than in it. The gnuplot *formula*
//! builders are pure string formatting and are ported.
//!
//! `DefaultParamHandler` likewise lives in the `param` module; the four numeric
//! parameters become
//! [`crate::math::posterior_error_probability::PepParameters`], whose `Default`
//! carries the source's `defaults_.setValue` values.
//!
//! # Differences from the source
//!
//! * Serial, as the source is: no `#pragma omp` anywhere in the translation
//!   unit.
//! * Every fitter call is fallible here, so `fit` returns `Result<bool>`: the
//!   `bool` is the source's return value and the error is a refusal the source
//!   does not make.

use crate::math::fitters::gauss::GaussFitResult;
use crate::math::fitters::gumbel_max_likelihood::{
    GumbelDistributionFitResult, GumbelMaxLikelihoodFitter,
};
use crate::math::statistic_functions::{
    mean, quantile1st_sorted, quantile3rd_sorted, sd_with_mean, sum,
};
use crate::{Error, Result};
use std::f64::consts::PI;

/// Maximum number of scores one fit may consume.
///
/// Native guard; the source allocates whatever it is handed.
pub const MAX_ITEMS: usize = 50_000_000;

/// The constant added to every score, on top of `|smallest score|`, to move the
/// sample strictly above zero before fitting.
///
/// The source's literal `0.001`, applied in `fit`, `fitGumbelGauss`,
/// `computeProbability` and `plotTargetDecoyEstimation`. A Gumbel density is
/// defined on the whole line, so the shift is about keeping the fitted location
/// away from zero rather than about the domain.
pub const SCORE_SHIFT: f64 = 0.001;

/// `!(a > b)`: true when `a <= b` **and** when either value is `NaN`.
///
/// Spelled through `partial_cmp` because clippy rejects the negated comparison
/// on a partially ordered type. The `NaN` arm is load-bearing at every call
/// site: the source writes `if (!(x > 0.0))` precisely so that a `NaN` takes
/// the guarded branch, and `x <= b` would not.
fn not_greater(a: f64, b: f64) -> bool {
    !matches!(a.partial_cmp(&b), Some(std::cmp::Ordering::Greater))
}

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.to_string())
}

/// Which distribution models the incorrectly assigned scores when the model is
/// *evaluated*.
///
/// The source's `incorrectly_assigned` parameter, valid strings `Gumbel` and
/// `Gauss`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum IncorrectComponent {
    /// Gumbel. The source's default, and what `computeProbability` evaluates
    /// unconditionally.
    #[default]
    Gumbel,
    /// Gaussian.
    Gauss,
}

impl IncorrectComponent {
    /// The parameter string the source accepts.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Gumbel => "Gumbel",
            Self::Gauss => "Gauss",
        }
    }

    /// Parse the source's parameter string, case-sensitively as
    /// `Param::setValidStrings` compares.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for anything but `Gumbel` or `Gauss`,
    /// which `DefaultParamHandler` reports as `Exception::InvalidParameter`.
    pub fn parse(text: &str) -> Result<Self> {
        match text {
            "Gumbel" => Ok(Self::Gumbel),
            "Gauss" => Ok(Self::Gauss),
            _ => Err(bad("incorrectly_assigned must be 'Gumbel' or 'Gauss'")),
        }
    }
}

/// What a fit does with scores far outside the interquartile range.
///
/// The source's `outlier_handling` parameter, passed as a string to every
/// `fit` overload.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum OutlierHandling {
    /// Drop everything outside `Q1 - 3 IQR` to `Q3 + 3 IQR`. The source's
    /// default.
    #[default]
    IgnoreIqrOutliers,
    /// Clamp everything outside that range to the nearest value inside it.
    SetIqrToClosestValid,
    /// Drop everything at or below the value at index `floor(n / 100) + 1` and
    /// at or above the value at index `floor(n * 99.9 / 100)` of the sorted
    /// scores. The parameter
    /// description says "99th and 1st percentile"; the code uses `99.9 / 100`
    /// for the upper index, and the comparisons are inclusive, so equal values
    /// at either end — censored maxima, for instance — go too.
    IgnoreExtremePercentiles,
    /// Leave the scores alone.
    None,
}

impl OutlierHandling {
    /// The parameter string the source accepts.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::IgnoreIqrOutliers => "ignore_iqr_outliers",
            Self::SetIqrToClosestValid => "set_iqr_to_closest_valid",
            Self::IgnoreExtremePercentiles => "ignore_extreme_percentiles",
            Self::None => "none",
        }
    }

    /// Parse the source's parameter string.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a string outside the four valid
    /// ones. The source compares `outlier_handling != "none"` and then tests
    /// the two named alternatives, so anything unrecognised silently selects
    /// `ignore_extreme_percentiles`; `DefaultParamHandler` is what would have
    /// rejected it, and this port rejects it here instead.
    pub fn parse(text: &str) -> Result<Self> {
        match text {
            "ignore_iqr_outliers" => Ok(Self::IgnoreIqrOutliers),
            "set_iqr_to_closest_valid" => Ok(Self::SetIqrToClosestValid),
            "ignore_extreme_percentiles" => Ok(Self::IgnoreExtremePercentiles),
            "none" => Ok(Self::None),
            _ => Err(bad("outlier_handling is not one of the four valid values")),
        }
    }
}

/// The model's parameters, the source's `DefaultParamHandler` defaults.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct PepParameters {
    /// Histogram resolution for the plots. Nothing in this module reads it;
    /// it is carried so that a caller assembling the source's plot output has
    /// the value in one place. Source default `100`.
    pub number_of_bins: u32,
    /// Which density models the incorrect component when the model is
    /// evaluated. Source default [`IncorrectComponent::Gumbel`].
    pub incorrectly_assigned: IncorrectComponent,
    /// Iteration cap for the EM loop. Source default `1000`.
    ///
    /// The loop tests `itns >= max_nr_iterations` *before* incrementing, so the
    /// body runs at most `max_nr_iterations + 1` times.
    pub max_nr_iterations: i32,
    /// Negative base-ten logarithm of the convergence threshold on the
    /// log-likelihood increase. Source default `6`, so the loop stops when the
    /// increase drops below `1e-6`.
    pub neg_log_delta: i32,
}

impl Default for PepParameters {
    fn default() -> Self {
        Self {
            number_of_bins: 100,
            incorrectly_assigned: IncorrectComponent::Gumbel,
            max_nr_iterations: 1000,
            neg_log_delta: 6,
        }
    }
}

/// Which density the negative component's gnuplot formula describes.
///
/// Stands in for the source's `getNegativeGnuplotFormula_` member function
/// pointer. `fit` sets it from
/// [`PepParameters::incorrectly_assigned`] and `fit_gumbel_gauss` always sets
/// it to Gumbel; the positive component's pointer is always the Gaussian, so it
/// needs no field.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum NegativeFormula {
    /// `getGumbelGnuplotFormula`. The source's initial value.
    #[default]
    Gumbel,
    /// `getGaussGnuplotFormula`.
    Gauss,
}

/// A mixture of an incorrect and a correct score component, fitted by EM.
#[derive(Clone, Debug)]
pub struct PosteriorErrorProbabilityModel {
    parameters: PepParameters,
    incorrectly_assigned_fit_param: GaussFitResult,
    incorrectly_assigned_fit_gumbel_param: GumbelDistributionFitResult,
    correctly_assigned_fit_param: GaussFitResult,
    negative_prior: f64,
    max_incorrectly: f64,
    max_correctly: f64,
    smallest_score: f64,
    negative_formula: NegativeFormula,
    probability_ready: bool,
}

impl Default for PosteriorErrorProbabilityModel {
    /// The source's default constructor: both Gaussian parameter sets at
    /// `(-1, -1, -1)`, the Gumbel at `(-1, -1)`, a negative prior of `0.5` and
    /// zeroed peaks. Those are markers, not usable parameters; the model is
    /// meaningless before a fit, as the header's `@note` on
    /// `computeProbability` says.
    fn default() -> Self {
        Self {
            parameters: PepParameters::default(),
            incorrectly_assigned_fit_param: GaussFitResult::default(),
            incorrectly_assigned_fit_gumbel_param: GumbelDistributionFitResult::new(-1.0, -1.0),
            correctly_assigned_fit_param: GaussFitResult::default(),
            negative_prior: 0.5,
            max_incorrectly: 0.0,
            max_correctly: 0.0,
            smallest_score: 0.0,
            negative_formula: NegativeFormula::Gumbel,
            probability_ready: false,
        }
    }
}

impl PosteriorErrorProbabilityModel {
    /// A model with the source's default parameters.
    pub fn new() -> Self {
        Self::default()
    }

    /// A model with explicit parameters, the source's `setParameters`.
    pub fn with_parameters(parameters: PepParameters) -> Self {
        Self {
            parameters,
            ..Self::default()
        }
    }

    /// The current parameters.
    pub fn parameters(&self) -> PepParameters {
        self.parameters
    }

    /// Replace the parameters, the source's `DefaultParamHandler::setParameters`.
    ///
    /// Fitted state is left untouched, as it is in the source.
    pub fn set_parameters(&mut self, parameters: PepParameters) {
        self.parameters = parameters;
    }

    /// Estimated Gaussian for the correctly assigned scores, the source's
    /// `getCorrectlyAssignedFitResult`.
    pub fn correctly_assigned_fit_result(&self) -> GaussFitResult {
        self.correctly_assigned_fit_param
    }

    /// Estimated Gaussian for the incorrectly assigned scores, the source's
    /// `getIncorrectlyAssignedFitResult`.
    ///
    /// This is what [`PosteriorErrorProbabilityModel::fit`] estimates even when
    /// [`IncorrectComponent::Gumbel`] is selected; in that case
    /// [`PosteriorErrorProbabilityModel::compute_probability`] reads its `x0`
    /// and `sigma` as the Gumbel's location and scale, which is what the
    /// header's `@note` about `GaussFitResult` holding all parameters means.
    pub fn incorrectly_assigned_fit_result(&self) -> GaussFitResult {
        self.incorrectly_assigned_fit_param
    }

    /// Estimated Gumbel for the incorrectly assigned scores, the source's
    /// `getIncorrectlyAssignedGumbelFitResult`. Only
    /// [`PosteriorErrorProbabilityModel::fit_gumbel_gauss`] sets it.
    pub fn incorrectly_assigned_gumbel_fit_result(&self) -> GumbelDistributionFitResult {
        self.incorrectly_assigned_fit_gumbel_param
    }

    /// Prior probability of the incorrect component, the source's
    /// `getNegativePrior`.
    pub fn negative_prior(&self) -> f64 {
        self.negative_prior
    }

    /// Smallest score of the last fit, the source's `getSmallestScore`.
    ///
    /// [`PosteriorErrorProbabilityModel::compute_probability`] re-applies the
    /// same shift this records, which is why a score must be passed to it
    /// untransformed.
    pub fn smallest_score(&self) -> f64 {
        self.smallest_score
    }

    /// Which density the negative gnuplot formula describes after the last fit.
    pub fn negative_formula(&self) -> NegativeFormula {
        self.negative_formula
    }

    /// The Gumbel density at `x` for parameters read as location `x0` and scale
    /// `sigma`, the source's static `getGumbel_`.
    ///
    /// `z = exp((x0 - x) / sigma)`, then `z exp(-z) / sigma`. The source
    /// divides by `sigma` twice without checking it; a zero or negative scale
    /// gives infinities or a reflected density there.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `sigma` is not finite and strictly
    /// positive or when `x` or `x0` is not finite.
    pub fn gumbel_density(x: f64, params: GaussFitResult) -> Result<f64> {
        if !params.sigma.is_finite() || not_greater(params.sigma, 0.0) {
            return Err(bad("Gumbel scale must be finite and positive"));
        }
        if !x.is_finite() || !params.x0.is_finite() {
            return Err(bad("Gumbel evaluation point must be finite"));
        }
        let z = ((params.x0 - x) / params.sigma).exp();
        Ok((z * (-z).exp()) / params.sigma)
    }

    /// Densities of both components at every score, the source's
    /// `fillDensities`.
    ///
    /// Returns `(incorrect, correct)`. Both are the *Gaussian* densities, even
    /// when [`IncorrectComponent::Gumbel`] is selected — the source's `TODO`
    /// admits as much at this very function.
    ///
    /// # Errors
    ///
    /// Whatever [`crate::math::fitters::gauss::GaussFitResult::eval`] returns:
    /// [`Error::InvalidValue`] for a non-positive or non-finite width, or a
    /// non-finite score.
    pub fn fill_densities(&self, x_scores: &[f64]) -> Result<(Vec<f64>, Vec<f64>)> {
        let mut incorrect = Vec::with_capacity(x_scores.len());
        let mut correct = Vec::with_capacity(x_scores.len());
        for score in x_scores.iter().copied() {
            incorrect.push(self.incorrectly_assigned_fit_param.eval(score)?);
            correct.push(self.correctly_assigned_fit_param.eval(score)?);
        }
        Ok((incorrect, correct))
    }

    /// Log densities of both components at every score, the source's
    /// `fillLogDensities`.
    ///
    /// Returns `(incorrect, correct)`, both from
    /// [`crate::math::fitters::gauss::GaussFitResult::log_eval_no_normalize`],
    /// which is the log of a proper normal density with the amplitude left out
    /// — the form a likelihood needs.
    ///
    /// # Errors
    ///
    /// As [`PosteriorErrorProbabilityModel::fill_densities`].
    pub fn fill_log_densities(&self, x_scores: &[f64]) -> Result<(Vec<f64>, Vec<f64>)> {
        let mut incorrect = Vec::with_capacity(x_scores.len());
        let mut correct = Vec::with_capacity(x_scores.len());
        for score in x_scores.iter().copied() {
            incorrect.push(
                self.incorrectly_assigned_fit_param
                    .log_eval_no_normalize(score)?,
            );
            correct.push(
                self.correctly_assigned_fit_param
                    .log_eval_no_normalize(score)?,
            );
        }
        Ok((incorrect, correct))
    }

    /// Log densities with a Gumbel incorrect component, the source's
    /// `fillLogDensitiesGumbel`.
    ///
    /// Returns `(incorrect, correct)`. Unlike
    /// [`PosteriorErrorProbabilityModel::fill_log_densities`], the incorrect
    /// component really is a Gumbel here, which is why
    /// [`PosteriorErrorProbabilityModel::fit_gumbel_gauss`] converges to a
    /// different optimum than [`PosteriorErrorProbabilityModel::fit`].
    ///
    /// # Errors
    ///
    /// As
    /// [`crate::math::fitters::gumbel_max_likelihood::GumbelDistributionFitResult::log_eval_no_normalize`]
    /// and
    /// [`crate::math::fitters::gauss::GaussFitResult::log_eval_no_normalize`].
    pub fn fill_log_densities_gumbel(&self, x_scores: &[f64]) -> Result<(Vec<f64>, Vec<f64>)> {
        let mut incorrect = Vec::with_capacity(x_scores.len());
        let mut correct = Vec::with_capacity(x_scores.len());
        for score in x_scores.iter().copied() {
            incorrect.push(
                self.incorrectly_assigned_fit_gumbel_param
                    .log_eval_no_normalize(score)?,
            );
            correct.push(
                self.correctly_assigned_fit_param
                    .log_eval_no_normalize(score)?,
            );
        }
        Ok((incorrect, correct))
    }

    /// Mixture log-likelihood from plain densities, the source's
    /// `computeLogLikelihood`.
    ///
    /// Sums `log10(prior * incorrect + (1 - prior) * correct)` in index order.
    ///
    /// The base-ten logarithm is the source's, and it is **not** the base the
    /// EM loop uses:
    /// [`PosteriorErrorProbabilityModel::compute_ll_and_incorrect_posteriors_from_log_densities`]
    /// accumulates natural logarithms, and that is the quantity the convergence
    /// threshold is compared against. Nothing in the source calls this
    /// function; it is kept because it is public API.
    ///
    /// Trailing entries of the longer slice are ignored: the source walks the
    /// `correct` range and advances an `incorrect` iterator alongside it, so a
    /// shorter `incorrect` would be read past its end. Here the shorter length
    /// wins instead.
    pub fn compute_log_likelihood(
        &self,
        incorrect_density: &[f64],
        correct_density: &[f64],
    ) -> f64 {
        let mut maxlike = 0.0;
        for (incorrect, correct) in incorrect_density.iter().zip(correct_density.iter()) {
            maxlike +=
                (self.negative_prior * incorrect + (1.0 - self.negative_prior) * correct).log10();
        }
        maxlike
    }

    /// Posterior probability of the incorrect component at every score, and the
    /// mixture log-likelihood; the source's
    /// `computeLLAndIncorrectPosteriorsFromLogDensities`.
    ///
    /// Both are computed by the log-sum-exp trick, subtracting the larger of
    /// the two log responsibilities before exponentiating, so neither
    /// underflows. The likelihood is a **natural** logarithm.
    ///
    /// Returns `(log_likelihood, incorrect_posteriors)`.
    pub fn compute_ll_and_incorrect_posteriors_from_log_densities(
        &self,
        incorrect_log_density: &[f64],
        correct_log_density: &[f64],
    ) -> (f64, Vec<f64>) {
        let mut loglikelihood = 0.0;
        let log_prior_pos = (1.0 - self.negative_prior).ln();
        let log_prior_neg = self.negative_prior.ln();
        let mut posteriors = Vec::with_capacity(incorrect_log_density.len());
        for (incorrect, correct) in incorrect_log_density.iter().zip(correct_log_density.iter()) {
            let mut log_resp_correct = log_prior_pos + correct;
            let mut log_resp_incorrect = log_prior_neg + incorrect;
            let max_log_resp = log_resp_correct.max(log_resp_incorrect);
            log_resp_correct -= max_log_resp;
            log_resp_incorrect -= max_log_resp;
            let resp_correct = log_resp_correct.exp();
            let resp_incorrect = log_resp_incorrect.exp();
            let total = resp_correct + resp_incorrect;
            posteriors.push(resp_incorrect / total);
            loglikelihood += max_log_resp + total.ln();
        }
        (loglikelihood, posteriors)
    }

    /// Posterior-weighted sums of the scores, the source's
    /// `pos_neg_mean_weighted_posteriors`.
    ///
    /// Returns `(correct, incorrect)` — the source's `pair.first` and
    /// `pair.second`. These are **sums, not means**: the caller divides by the
    /// summed posteriors afterwards, which is why the source's own name is
    /// misleading and why the division cannot be folded in here without
    /// changing the accumulation order.
    pub fn pos_neg_mean_weighted_posteriors(
        x_scores: &[f64],
        incorrect_posteriors: &[f64],
    ) -> (f64, f64) {
        let mut pos_x0 = 0.0;
        let mut neg_x0 = 0.0;
        for (incorrect, score) in incorrect_posteriors.iter().zip(x_scores.iter()) {
            pos_x0 += (1.0 - incorrect) * score;
            neg_x0 += incorrect * score;
        }
        (pos_x0, neg_x0)
    }

    /// Posterior-weighted sums of squared deviations, the source's
    /// `pos_neg_sigma_weighted_posteriors`.
    ///
    /// Returns `(correct, incorrect)`, again as sums; the caller divides by the
    /// summed posteriors and takes the square root.
    ///
    /// # Arguments
    ///
    /// * `pos_neg_mean` — the already-divided means, `(correct, incorrect)`.
    pub fn pos_neg_sigma_weighted_posteriors(
        x_scores: &[f64],
        incorrect_posteriors: &[f64],
        pos_neg_mean: (f64, f64),
    ) -> (f64, f64) {
        let mut pos_sigma = 0.0;
        let mut neg_sigma = 0.0;
        for (incorrect, score) in incorrect_posteriors.iter().zip(x_scores.iter()) {
            pos_sigma += (1.0 - incorrect) * (score - pos_neg_mean.0).powi(2);
            neg_sigma += incorrect * (score - pos_neg_mean.1).powi(2);
        }
        (pos_sigma, neg_sigma)
    }

    /// Posterior error probability of a raw, untransformed score.
    ///
    /// Applies the fit's own shift, `score + |smallest_score| + 0.001`, then
    /// returns `prior * f_incorrect / (prior * f_incorrect + (1 - prior) *
    /// f_correct)`.
    ///
    /// Below the incorrect component's location the incorrect density is
    /// replaced by its peak value, and above the correct component's location
    /// the correct density is; both keep the probability monotone where the
    /// densities themselves would turn it back. The incorrect density is
    /// evaluated as a **Gumbel in every branch**, whatever
    /// [`IncorrectComponent`] says — a second `TODO` in the source admits this
    /// and calls it "confusing at best". It is preserved.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `score` is not finite, when either
    /// component's width is not positive — which is the state a
    /// default-constructed model is in, and which the source would silently
    /// carry into `boost::math::normal_distribution` — or when the mixture
    /// density is zero, where the source divides by zero and returns `NaN`.
    pub fn compute_probability(&self, score: f64) -> Result<f64> {
        if !self.probability_ready {
            return Err(Error::MissingInformation(
                "the model has no usable incorrect-component peak; call fit() first".to_string(),
            ));
        }
        if !score.is_finite() {
            return Err(bad("score must be finite"));
        }
        let shifted = score + self.smallest_score.abs() + SCORE_SHIFT;
        let (x_neg, x_pos) = if shifted < self.incorrectly_assigned_fit_param.x0 {
            (
                self.max_incorrectly,
                self.correctly_assigned_fit_param.eval(shifted)?,
            )
        } else if shifted > self.correctly_assigned_fit_param.x0 {
            (
                Self::gumbel_density(shifted, self.incorrectly_assigned_fit_param)?,
                self.max_correctly,
            )
        } else {
            (
                Self::gumbel_density(shifted, self.incorrectly_assigned_fit_param)?,
                self.correctly_assigned_fit_param.eval(shifted)?,
            )
        };
        let numerator = self.negative_prior * x_neg;
        let denominator = numerator + (1.0 - self.negative_prior) * x_pos;
        if denominator == 0.0 {
            return Err(bad(
                "the mixture density vanishes at this score; the model is not usable here",
            ));
        }
        Ok(numerator / denominator)
    }

    /// The gnuplot expression for a Gumbel with location `x0` and scale
    /// `sigma`, the source's `getGumbelGnuplotFormula`.
    ///
    /// Numbers are formatted the way `std::ostream << double` formats them —
    /// `printf("%g")` with six significant digits — so the expression is
    /// byte-identical to the source's.
    pub fn gumbel_gnuplot_formula(params: GaussFitResult) -> String {
        let sigma = format_g(params.sigma);
        let x0 = format_g(params.x0);
        format!("(1/{sigma}) * exp(( {x0}- x)/{sigma}) * exp(-exp(({x0} - x)/{sigma}))")
    }

    /// The gnuplot expression for a Gaussian, the source's
    /// `getGaussGnuplotFormula`.
    pub fn gauss_gnuplot_formula(params: GaussFitResult) -> String {
        format!(
            "{} * exp(-(x - {}) ** 2 / 2 / ({}) ** 2)",
            format_g(params.a),
            format_g(params.x0),
            format_g(params.sigma)
        )
    }

    /// The gnuplot expression for the whole mixture, the source's
    /// `getBothGnuplotFormula`.
    ///
    /// The negative half follows [`PosteriorErrorProbabilityModel::negative_formula`],
    /// which the last fit set; the positive half is always the Gaussian.
    pub fn both_gnuplot_formula(
        &self,
        incorrect: GaussFitResult,
        correct: GaussFitResult,
    ) -> String {
        let negative = match self.negative_formula {
            NegativeFormula::Gumbel => Self::gumbel_gnuplot_formula(incorrect),
            NegativeFormula::Gauss => Self::gauss_gnuplot_formula(incorrect),
        };
        format!(
            "{}*{} + (1-{})*{}",
            format_g(self.negative_prior),
            negative,
            format_g(self.negative_prior),
            Self::gauss_gnuplot_formula(correct)
        )
    }

    /// Fit two Gaussians by EM, the source's
    /// `fit(std::vector<double>&, const std::string&)`.
    ///
    /// `search_engine_scores` is sorted ascending in place, as the header's
    /// `@note` says. Everything after that works on a shifted copy:
    /// `score + |smallest| + 0.001`, with outliers handled per
    /// `outlier_handling`.
    ///
    /// Initialisation is the source's, literally: the incorrect component's
    /// location is the mean of the lower half of the shifted scores *plus the
    /// smallest of them*, its width is the deviation of all scores about that
    /// location, and the correct component starts at the mean of the top 30%
    /// plus the score at that boundary, sharing the incorrect component's
    /// width. The negative prior starts at `0.7`. These are not textbook
    /// starting values and a mixture fit lands on a different local optimum if
    /// they are changed, which is why they are transcribed rather than
    /// rationalised.
    ///
    /// The loop stops when the log-likelihood increase falls below
    /// `10^-neg_log_delta`, when the increase turns negative, or when the
    /// iteration count reaches [`PepParameters::max_nr_iterations`].
    ///
    /// Returns `true` when the fit ran to completion. The source returns
    /// `false` only for an empty input, a `NaN` likelihood step, or a
    /// likelihood that decreased; note that the "impossible standard
    /// deviations" branch `break`s out of the loop with `good_fit` still
    /// `true`, so it reports success. That is preserved.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a score is not finite, when outlier
    /// handling removes every score — which leaves the source indexing an empty
    /// vector — or when a component evaluation fails, and
    /// [`Error::InvalidRange`] when the input exceeds [`MAX_ITEMS`].
    pub fn fit(
        &mut self,
        search_engine_scores: &mut [f64],
        outlier_handling: OutlierHandling,
    ) -> Result<bool> {
        let Some(x_scores) = self.prepare_scores(search_engine_scores, outlier_handling)? else {
            return Ok(false);
        };
        let n = x_scores.len();

        self.negative_prior = 0.7;
        let half = (0.5 * n as f64).ceil() as usize;
        self.incorrectly_assigned_fit_param.x0 = mean(&x_scores[..half])? + x_scores[0];
        self.incorrectly_assigned_fit_param.sigma =
            sd_with_mean(&x_scores, self.incorrectly_assigned_fit_param.x0)?;
        self.incorrectly_assigned_fit_param.a =
            1.0 / (2.0 * PI * self.incorrectly_assigned_fit_param.sigma.powi(2)).sqrt();
        // Both branches of the source compute the same Gaussian start; only the
        // plot formula differs.
        self.negative_formula = match self.parameters.incorrectly_assigned {
            IncorrectComponent::Gumbel => NegativeFormula::Gumbel,
            IncorrectComponent::Gauss => NegativeFormula::Gauss,
        };

        let start = (n - 1).min((n as f64 * 0.7).ceil() as usize);
        self.correctly_assigned_fit_param.x0 = mean(&x_scores[start..])? + x_scores[start];
        self.correctly_assigned_fit_param.sigma = self.incorrectly_assigned_fit_param.sigma;
        self.correctly_assigned_fit_param.a =
            1.0 / (2.0 * PI * self.correctly_assigned_fit_param.sigma.powi(2)).sqrt();

        let mut good_fit = true;
        let mut stop = false;
        let max_itns = self.parameters.max_nr_iterations;
        let delta = self.parameters.neg_log_delta;
        let mut itns: i32 = 0;

        let (mut incorrect_log, mut correct_log) = self.fill_log_densities(&x_scores)?;
        let (mut maxlike, mut incorrect_posteriors) = self
            .compute_ll_and_incorrect_posteriors_from_log_densities(&incorrect_log, &correct_log);
        let mut sum_incorrect = sum(&incorrect_posteriors);
        let mut sum_correct = n as f64 - sum_incorrect;

        loop {
            let mut new_means =
                Self::pos_neg_mean_weighted_posteriors(&x_scores, &incorrect_posteriors);
            new_means.0 /= sum_correct;
            new_means.1 /= sum_incorrect;
            let mut new_sigmas = Self::pos_neg_sigma_weighted_posteriors(
                &x_scores,
                &incorrect_posteriors,
                new_means,
            );
            new_sigmas.0 = (new_sigmas.0 / sum_correct).sqrt();
            new_sigmas.1 = (new_sigmas.1 / sum_incorrect).sqrt();
            if not_greater(new_sigmas.0, 0.0)
                || not_greater(new_sigmas.1, 0.0)
                || new_sigmas.0.is_nan()
                || new_sigmas.1.is_nan()
            {
                // "Warning: encountered impossible standard deviations.
                // Aborting fit." - and `good_fit` stays true.
                break;
            }

            self.correctly_assigned_fit_param.x0 = new_means.0;
            self.incorrectly_assigned_fit_param.x0 = new_means.1;
            self.correctly_assigned_fit_param.sigma = new_sigmas.0;
            self.correctly_assigned_fit_param.a =
                1.0 / (2.0 * PI * self.correctly_assigned_fit_param.sigma.powi(2)).sqrt();
            self.incorrectly_assigned_fit_param.sigma = new_sigmas.1;
            self.incorrectly_assigned_fit_param.a =
                1.0 / (2.0 * PI * self.incorrectly_assigned_fit_param.sigma.powi(2)).sqrt();

            let filled = self.fill_log_densities(&x_scores)?;
            incorrect_log = filled.0;
            correct_log = filled.1;
            let (new_maxlike, posteriors) = self
                .compute_ll_and_incorrect_posteriors_from_log_densities(
                    &incorrect_log,
                    &correct_log,
                );
            incorrect_posteriors = posteriors;
            sum_incorrect = sum(&incorrect_posteriors);
            sum_correct = n as f64 - sum_incorrect;
            self.negative_prior = sum_incorrect / n as f64;

            if (new_maxlike - maxlike).is_nan() {
                return Ok(false);
            }
            if (new_maxlike - maxlike) < 10.0_f64.powi(-delta) || itns >= max_itns {
                stop = true;
                good_fit = true;
            } else if new_maxlike < maxlike {
                stop = true;
                good_fit = false;
            }
            maxlike = new_maxlike;
            itns += 1;
            if stop {
                break;
            }
        }

        self.max_incorrectly = match self.parameters.incorrectly_assigned {
            IncorrectComponent::Gumbel => Self::gumbel_density(
                self.incorrectly_assigned_fit_param.x0,
                self.incorrectly_assigned_fit_param,
            )?,
            IncorrectComponent::Gauss => self
                .incorrectly_assigned_fit_param
                .eval(self.incorrectly_assigned_fit_param.x0)?,
        };
        self.max_correctly = self
            .correctly_assigned_fit_param
            .eval(self.correctly_assigned_fit_param.x0)?;
        self.probability_ready = true;
        Ok(good_fit)
    }

    /// Fit and then evaluate, the source's
    /// `fit(std::vector<double>&, std::vector<double>&, const std::string&)`.
    ///
    /// Returns the posterior error probability of every score, in the sorted
    /// order the fit leaves `search_engine_scores` in, or `None` when
    /// [`PosteriorErrorProbabilityModel::fit`] reported failure.
    ///
    /// # Errors
    ///
    /// As [`PosteriorErrorProbabilityModel::fit`] and
    /// [`PosteriorErrorProbabilityModel::compute_probability`].
    pub fn fit_with_probabilities(
        &mut self,
        search_engine_scores: &mut [f64],
        outlier_handling: OutlierHandling,
    ) -> Result<Option<Vec<f64>>> {
        if !self.fit(search_engine_scores, outlier_handling)? {
            return Ok(None);
        }
        let mut probabilities = Vec::with_capacity(search_engine_scores.len());
        for score in search_engine_scores.iter().copied() {
            probabilities.push(self.compute_probability(score)?);
        }
        Ok(Some(probabilities))
    }

    /// Fit a Gumbel and a Gaussian by EM, the source's `fitGumbelGauss`.
    ///
    /// As [`PosteriorErrorProbabilityModel::fit`], except that the incorrect
    /// component is a genuine Gumbel refitted each iteration by
    /// [`crate::math::fitters::gumbel_max_likelihood::GumbelMaxLikelihoodFitter::fit_weighted`]
    /// against the current posteriors, the correct component's mean and width
    /// are updated by the same weighted formulae written out inline, and the
    /// negative prior starts at `0.7`. The fitter is constructed once, before
    /// the loop, from the initial Gumbel parameters, so each iteration starts
    /// from the previous result — the source relies on `fitWeighted`
    /// overwriting the fitter's own start parameters.
    ///
    /// A non-positive or `NaN` Gumbel scale aborts the loop with `good_fit`
    /// still `true`, exactly as in
    /// [`PosteriorErrorProbabilityModel::fit`].
    ///
    /// The peak of the incorrect component is computed at the end from
    /// `incorrectly_assigned_fit_param_`, the *Gaussian* parameter set, which
    /// this function never writes: it is still `(-1, -1, -1)` unless a previous
    /// `fit` left something there. That is a source defect, recorded in
    /// `OpenMS_CPP_ISSUES.md`; the port refuses rather than reproducing the
    /// nonsense, which is why a `fitGumbelGauss` on a fresh model returns an
    /// error where the C++ returns `true` with a meaningless
    /// `max_incorrectly_`.
    ///
    /// # Errors
    ///
    /// As [`PosteriorErrorProbabilityModel::fit`], plus whatever the Gumbel
    /// fitter returns, plus [`Error::InvalidValue`] when the Gaussian parameter
    /// set that the final peak calculation reads has not been fitted.
    pub fn fit_gumbel_gauss(
        &mut self,
        search_engine_scores: &mut [f64],
        outlier_handling: OutlierHandling,
    ) -> Result<bool> {
        let Some(x_scores) = self.prepare_scores(search_engine_scores, outlier_handling)? else {
            return Ok(false);
        };
        let n = x_scores.len();

        let half = (0.5 * n as f64).ceil() as usize;
        self.incorrectly_assigned_fit_gumbel_param.a = mean(&x_scores[..half])? + x_scores[0];
        self.incorrectly_assigned_fit_gumbel_param.b =
            sd_with_mean(&x_scores, self.incorrectly_assigned_fit_gumbel_param.a)?;
        self.negative_prior = 0.7;
        self.negative_formula = NegativeFormula::Gumbel;

        let start = (n - 1).min((n as f64 * 0.7).ceil() as usize);
        self.correctly_assigned_fit_param.x0 = mean(&x_scores[start..])? + x_scores[start];
        self.correctly_assigned_fit_param.sigma = self.incorrectly_assigned_fit_gumbel_param.b;
        self.correctly_assigned_fit_param.a =
            1.0 / (2.0 * PI * self.correctly_assigned_fit_param.sigma.powi(2)).sqrt();

        let mut good_fit = true;
        let mut stop = false;
        let max_itns = self.parameters.max_nr_iterations;
        let delta = self.parameters.neg_log_delta;
        let mut itns: i32 = 0;

        let (mut incorrect_log, mut correct_log) = self.fill_log_densities_gumbel(&x_scores)?;
        let (mut maxlike, mut incorrect_posteriors) = self
            .compute_ll_and_incorrect_posteriors_from_log_densities(&incorrect_log, &correct_log);
        let mut sum_incorrect = sum(&incorrect_posteriors);
        let mut sum_correct = n as f64 - sum_incorrect;

        let mut fitter = GumbelMaxLikelihoodFitter::with_initial_parameters(
            self.incorrectly_assigned_fit_gumbel_param,
        );

        loop {
            let mut new_gauss_mean = 0.0;
            for (incorrect, score) in incorrect_posteriors.iter().zip(x_scores.iter()) {
                new_gauss_mean += (1.0 - incorrect) * score;
            }
            new_gauss_mean /= sum_correct;

            let mut new_gauss_sigma = 0.0;
            for (incorrect, score) in incorrect_posteriors.iter().zip(x_scores.iter()) {
                new_gauss_sigma += (1.0 - incorrect) * (score - new_gauss_mean).powi(2);
            }
            new_gauss_sigma = (new_gauss_sigma / sum_correct).sqrt();

            let new_gumbel = fitter.fit_weighted(&x_scores, &incorrect_posteriors)?;
            if new_gumbel.b <= 0.0 || new_gumbel.b.is_nan() {
                // "Warning: encountered impossible standard deviations."
                break;
            }

            self.correctly_assigned_fit_param.x0 = new_gauss_mean;
            self.correctly_assigned_fit_param.sigma = new_gauss_sigma;
            self.correctly_assigned_fit_param.a = 1.0 / (2.0 * PI * new_gauss_sigma.powi(2)).sqrt();
            self.incorrectly_assigned_fit_gumbel_param = new_gumbel;

            let filled = self.fill_log_densities_gumbel(&x_scores)?;
            incorrect_log = filled.0;
            correct_log = filled.1;
            let (new_maxlike, posteriors) = self
                .compute_ll_and_incorrect_posteriors_from_log_densities(
                    &incorrect_log,
                    &correct_log,
                );
            incorrect_posteriors = posteriors;
            sum_incorrect = sum(&incorrect_posteriors);
            sum_correct = n as f64 - sum_incorrect;
            self.negative_prior = sum_incorrect / n as f64;

            if (new_maxlike - maxlike).is_nan() {
                return Ok(false);
            }
            if (new_maxlike - maxlike) < 10.0_f64.powi(-delta) || itns >= max_itns {
                stop = true;
                good_fit = true;
            } else if new_maxlike < maxlike {
                stop = true;
                good_fit = false;
            }
            maxlike = new_maxlike;
            itns += 1;
            if stop {
                break;
            }
        }

        // The source finishes with
        //   max_incorrectly_ = getGumbel_(incorrectly_assigned_fit_param_.x0,
        //                                 incorrectly_assigned_fit_param_);
        // reading the *Gaussian* parameter set, which this function never
        // writes. On a freshly constructed model that set is still
        // `(-1, -1, -1)` and the expression evaluates to -0.3679 - a negative
        // "peak density" that `computeProbability` would then divide by. The
        // peak is computed here only when a previous `fit` left a usable
        // Gaussian set, and the model is marked not ready for
        // `compute_probability` either way.
        if self.incorrectly_assigned_fit_param.sigma > 0.0 {
            self.max_incorrectly = Self::gumbel_density(
                self.incorrectly_assigned_fit_param.x0,
                self.incorrectly_assigned_fit_param,
            )?;
        }
        self.max_correctly = self
            .correctly_assigned_fit_param
            .eval(self.correctly_assigned_fit_param.x0)?;
        self.probability_ready = false;
        Ok(good_fit)
    }

    /// Sort, shift and outlier-process the scores, shared by both fits.
    ///
    /// Returns `None` for an empty input, which both fits report as `false`.
    fn prepare_scores(
        &mut self,
        search_engine_scores: &mut [f64],
        outlier_handling: OutlierHandling,
    ) -> Result<Option<Vec<f64>>> {
        if search_engine_scores.len() > MAX_ITEMS {
            return Err(Error::InvalidRange(format!(
                "posterior error probability fit: {} scores exceeds the maximum {MAX_ITEMS}",
                search_engine_scores.len()
            )));
        }
        if search_engine_scores.is_empty() {
            return Ok(None);
        }
        if search_engine_scores.iter().any(|v| !v.is_finite()) {
            return Err(bad("search engine scores must be finite"));
        }
        search_engine_scores.sort_by(f64::total_cmp);
        self.smallest_score = search_engine_scores[0];
        let shift = self.smallest_score.abs() + SCORE_SHIFT;
        let mut x_scores: Vec<f64> = search_engine_scores.iter().map(|v| v + shift).collect();
        process_outliers(&mut x_scores, outlier_handling)?;
        if x_scores.is_empty() {
            return Err(bad(
                "outlier handling removed every score; the source would index an empty vector here",
            ));
        }
        Ok(Some(x_scores))
    }
}

/// The source's `processOutliers_`, operating on the already-sorted shifted
/// scores.
///
/// The quartiles are the median-of-halves
/// [`crate::math::statistic_functions::quantile1st_sorted`] and
/// [`crate::math::statistic_functions::quantile3rd_sorted`], not the
/// interpolating quantile, because the source passes `sorted = true` to those
/// two functions specifically.
fn process_outliers(x_scores: &mut Vec<f64>, handling: OutlierHandling) -> Result<()> {
    if x_scores.is_empty() || handling == OutlierHandling::None {
        return Ok(());
    }
    let q1 = quantile1st_sorted(x_scores)?;
    let q3 = quantile3rd_sorted(x_scores)?;
    let iqr = q3 - q1;
    match handling {
        OutlierHandling::IgnoreIqrOutliers => {
            let lower = q1 - 3.0 * iqr;
            let upper = q3 + 3.0 * iqr;
            x_scores.retain(|v| !(*v < lower || *v > upper));
        }
        OutlierHandling::SetIqrToClosestValid => {
            let lower = q1 - 3.0 * iqr;
            let upper = q3 + 3.0 * iqr;
            // `lower_bound` / `--upper_bound` over the sorted scores.
            let first = x_scores.partition_point(|v| *v < lower);
            let last = x_scores.partition_point(|v| *v <= upper);
            if first >= x_scores.len() || last == 0 {
                // The source decrements `upper_bound` unconditionally and
                // dereferences `lower_bound`, so an all-outlier sample walks
                // off one end or the other.
                return Err(bad(
                    "every score is an outlier; there is no valid value to clamp to",
                ));
            }
            let low_value = x_scores[first];
            let high_value = x_scores[last - 1];
            for value in x_scores[..first].iter_mut() {
                *value = low_value;
            }
            for value in x_scores[last..].iter_mut() {
                *value = high_value;
            }
        }
        OutlierHandling::IgnoreExtremePercentiles => {
            let n = x_scores.len();
            let upper_index = ((n as f64) * 99.9 / 100.0) as usize;
            let lower_index = ((n as f64) / 100.0) as usize + 1;
            if upper_index >= n || lower_index >= n {
                // `x_scores[first_idx]` is an unchecked index in the source.
                // `first_idx = floor(n / 100) + 1` and
                // `ninetyninth_idx = floor(n * 99.9 / 100)`, so the only length
                // that reads past the end is `n == 1`, where `first_idx == 1`;
                // the guard is written against both indices rather than against
                // that one length so it stays correct if either expression
                // changes.
                return Err(bad(
                    "percentile-based outlier handling needs more than one score; the source indexes out of range here",
                ));
            }
            let upper_value = x_scores[upper_index];
            let lower_value = x_scores[lower_index];
            x_scores.retain(|v| !(*v <= lower_value || *v >= upper_value));
        }
        OutlierHandling::None => {}
    }
    Ok(())
}

/// Format a `double` the way an unconfigured `std::ostream` does:
/// `printf("%g")` with six significant digits.
///
/// Rust's `{}` prints the shortest string that round-trips, which is not what
/// the source writes into its gnuplot files; a formula built with it would be a
/// different byte sequence for the same fit.
fn format_g(value: f64) -> String {
    if value.is_nan() {
        return "nan".to_string();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.to_string();
    }
    const PRECISION: i32 = 6;
    let scientific = format!("{:.*e}", (PRECISION - 1) as usize, value);
    let Some(marker) = scientific.find('e') else {
        return scientific;
    };
    let exponent: i32 = scientific[marker + 1..].parse().unwrap_or(0);
    if (-4..PRECISION).contains(&exponent) {
        let decimals = (PRECISION - 1 - exponent).max(0) as usize;
        trim_trailing_zeros(format!("{value:.decimals$}"))
    } else {
        let mantissa = trim_trailing_zeros(scientific[..marker].to_string());
        let sign = if exponent < 0 { '-' } else { '+' };
        format!("{mantissa}e{sign}{:02}", exponent.abs())
    }
}

/// Drop the trailing zeros of a fractional part, and the point if nothing is
/// left after it — the `%g` rule, absent the `#` flag.
fn trim_trailing_zeros(text: String) -> String {
    if !text.contains('.') {
        return text;
    }
    let trimmed = text.trim_end_matches('0');
    let trimmed = trimmed.strip_suffix('.').unwrap_or(trimmed);
    trimmed.to_string()
}
