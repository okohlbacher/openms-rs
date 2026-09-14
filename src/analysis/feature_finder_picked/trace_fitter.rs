// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The retention-time shape model shared by the picked feature finder's trace
//! fitters (`FEATUREFINDER/TraceFitter.h`).
//!
//! Registered ahead of its early-TOPP-bundle work package. The integrator
//! transcribed the public interface of the pinned header (core `bc9cc12`) as the
//! `TraceFitter` trait and the `TraceFitterParams` record below, so that package
//! B4-GAUSS (`GaussTraceFitter`) and package B5-EGH (`EGHTraceFitter`) can build
//! against one contract at the same time. Only declarations live here so far:
//! B4-GAUSS owns this file and fills it with the Levenberg-Marquardt driver
//! (source `optimize_`), the parameter defaults and their `Param` mapping, and
//! any helper the two fitters share. B5-EGH reads this file and does not edit
//! it.
//!
//! Changing a signature of the trait or the fields of the parameter record needs
//! the integrator, because both fitters and the later feature-extension package
//! depend on them.
//!
//! # Source members that are not trait methods
//!
//! - `GenericFunctor`, the raw-pointer residual and Jacobian functor, becomes
//!   the pair of closures passed to the crate's Levenberg-Marquardt `minimize`.
//! - `ModelData`, the protected bundle of the traces and the weighting flag, is
//!   internal to each fitter.
//! - `getOptimizedParameters_`, the protected hook that copies the optimised
//!   vector back, is a private method of each fitter.
//! - `optimize_` becomes a free function of this module with the planned
//!   signature `optimize(x: &mut [f64], values: usize, residual, jacobian,
//!   parameters: &TraceFitterParams) -> Result<()>`. It refuses fewer residuals
//!   than parameters, and every Levenberg-Marquardt status up to and including
//!   `ImproperInputParameters`, with the source's `UnableToFit-FinalSet`
//!   condition; it accepts every later status, including the exhausted
//!   evaluation budget.
//! - The copy constructor and assignment operator become `Clone` on each
//!   fitter; `updateMembers_` becomes the typed parameter setter.
//!
//! The source is serial, and so are the fitters.

use crate::Result;
use crate::analysis::feature_finder_picked::helper_structs::{MassTrace, MassTraces};

/// The parameters every trace fitter reads: source `TraceFitter`'s
/// `DefaultParamHandler` defaults and the members `max_iterations_` and
/// `weighted_` that `updateMembers_` refreshes from them.
///
/// The source registers `max_iteration` with the default `500` and the
/// `advanced` tag, and `weighted` as the string `"false"` with the valid
/// strings `"true"` and `"false"`, also `advanced`. Package B4-GAUSS adds the
/// default value and the mapping to and from `Param`; the record deliberately
/// derives no `Default`, because a derived default would give `0` and `false`
/// instead of the source's `500` and `false`.
///
/// `FeatureFinderAlgorithmPicked` sets only `max_iteration`, from its own
/// `fit:max_iterations` (at least 1), so the tool path always fits unweighted.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TraceFitterParams {
    /// Evaluation budget of the Levenberg-Marquardt fit: source parameter
    /// `max_iteration` (singular), default 500.
    ///
    /// Despite its name the source passes it to Eigen as `maxfev`, the number
    /// of residual evaluations, not iterations. It is signed, as the source's
    /// `SignedSize` member: the source parameter carries no minimum, and a value
    /// of zero or below makes every fit fail, because Eigen refuses
    /// `maxfev <= 0` as improper input.
    pub max_iteration: i64,
    /// Whether each trace's residuals are weighted by its theoretical
    /// intensity during the fit: source parameter `weighted`, the string
    /// `"true"` or `"false"`, default `"false"`.
    pub weighted: bool,
}

/// A retention-time shape model fitted to the mass traces of one feature
/// candidate: source abstract class `TraceFitter`.
///
/// Implementors supply a concrete profile, a Gaussian (`GaussTraceFitter`) or an
/// exponential-Gaussian hybrid (`EGHTraceFitter`), fit it with [`Self::fit`],
/// and answer the geometric queries, the span checks and the plot helper from
/// the fitted parameters. `FeatureFinderAlgorithmPicked` chooses the
/// implementation from `feature:rt_shape`, sets its parameters, fits the traces
/// and then calls only the methods of this trait, so the trait is object safe:
/// no method is generic and none returns `Self`.
///
/// The queries are meaningful only after a successful [`Self::fit`]; before it
/// the source reads uninitialised members. Model-specific accessors, such as the
/// EGH `getTau` and `getSigma` that the algorithm reaches through a dynamic
/// cast, are inherent methods of the implementing type.
///
/// Signatures here are part of the wave-2 integration contract. A change needs
/// the integrator, because package B4-GAUSS, package B5-EGH and the later
/// feature-extension package build against them in parallel.
pub trait TraceFitter {
    /// The parameters the next [`Self::fit`] uses: the typed counterpart of the
    /// source's inherited `DefaultParamHandler::getParameters`.
    fn parameters(&self) -> &TraceFitterParams;

    /// Replace the parameters: the typed counterpart of the source's
    /// `DefaultParamHandler::setParameters` followed by `updateMembers_`.
    ///
    /// Validation of textual parameter values, which the source performs
    /// inside `setParameters`, happens where a `Param` is converted into a
    /// [`TraceFitterParams`], so this setter cannot fail.
    fn set_parameters(&mut self, parameters: TraceFitterParams);

    /// Fit the model to `traces` with Levenberg-Marquardt: source pure virtual
    /// `fit`.
    ///
    /// The implementation derives initial parameters from the traces, runs the
    /// shared driver with the evaluation budget
    /// [`TraceFitterParams::max_iteration`], and stores the optimised
    /// parameters, which the query methods then read. The source takes the
    /// traces by non-const reference and says subclasses may set bookkeeping on
    /// them; neither ported subclass writes to them, so they are borrowed
    /// immutably here.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] whose message starts with
    /// `UnableToFit-FinalSet` where the source throws `Exception::UnableToFit`:
    /// when the traces hold fewer peaks than the model has parameters, or when
    /// the Levenberg-Marquardt status is `ImproperInputParameters` or earlier.
    /// Errors from the helper structures, such as their resource ceilings,
    /// propagate unchanged. Each implementation documents what its fitted
    /// parameters hold after an error.
    fn fit(&mut self, traces: &MassTraces) -> Result<()>;

    /// Lower retention-time bound of the fitted model, in seconds: source
    /// `getLowerRTBound`.
    fn lower_rt_bound(&self) -> f64;

    /// Upper retention-time bound of the fitted model, in seconds: source
    /// `getUpperRTBound`.
    fn upper_rt_bound(&self) -> f64;

    /// Height of the fitted model at its centre: source `getHeight`.
    fn height(&self) -> f64;

    /// Centre retention time of the fitted model, in seconds: source
    /// `getCenter`.
    fn center(&self) -> f64;

    /// Full width at half maximum of the fitted model, in seconds: source
    /// `getFWHM`.
    fn fwhm(&self) -> f64;

    /// Model intensity at retention time `rt`, in seconds: source `getValue`.
    fn value(&self, rt: f64) -> f64;

    /// Integrated intensity of the fitted model: source `getArea`.
    ///
    /// The source declares it non-const, but neither subclass changes state, so
    /// it takes `&self`.
    fn area(&self) -> f64;

    /// The minimal retention-time span check: source `checkMinimalRTSpan`.
    ///
    /// `rt_bounds` is a `(lower, upper)` retention-time range in seconds, and
    /// `min_rt_span` a fraction of the fitted model's width.
    /// `FeatureFinderAlgorithmPicked` passes the retention-time bounds of the
    /// cropped feature traces. The source header documents `true` as the model
    /// spanning enough of the range, but both subclasses return `true` when the
    /// range is narrower than `min_rt_span` times the model's width (five sigma
    /// for the Gaussian; for EGH the distance between its bounds at relative
    /// height 0.043937, where a Gaussian is 2.5 sigma from its centre), and the
    /// algorithm rejects the fit on `true` ("Less than 'min_rt_span' left after
    /// fit"). Implementations keep the subclasses' comparison. The source
    /// declares it non-const; it reads state only.
    fn check_minimal_rt_span(&self, rt_bounds: (f64, f64), min_rt_span: f64) -> bool;

    /// The maximal retention-time span check: source `checkMaximalRTSpan`.
    ///
    /// `max_rt_span` is a fraction of the retention-time span of the region the
    /// fit was initialised from. The source header documents `true` as the model
    /// staying within that fraction, but both subclasses return `true` when the
    /// model's width exceeds it, and the algorithm rejects the fit on `true`
    /// ("Fitted model is bigger than 'max_rt_span'"). Implementations keep the
    /// subclasses' comparison. The source declares it non-const; it reads state
    /// only.
    fn check_maximal_rt_span(&self, max_rt_span: f64) -> bool;

    /// A gnuplot expression of the fitted model for `trace`: source
    /// `getGnuplotFormula`.
    ///
    /// `function_name` names the function (`'f'` gives `f(x)= ...`), `baseline`
    /// is the intensity added to the model, and `rt_shift` is added to the
    /// model's centre so several traces can be plotted side by side (typically
    /// `0` for the first trace). The trace's theoretical intensity scales the
    /// height. Numbers are written as the source's default `std::ostream`
    /// does, precision 6 in `%g` style. The source declares it non-const; it
    /// reads state only.
    fn gnuplot_formula(
        &self,
        trace: &MassTrace,
        function_name: char,
        baseline: f64,
        rt_shift: f64,
    ) -> String;

    /// The model value at the retention time of peak `k` of `trace`, scaled by
    /// the trace's theoretical intensity: source `computeTheoretical`, that is
    /// `trace.theoretical_int * value(trace.peaks[k].rt)`.
    ///
    /// This is the source's one non-virtual member; implementors keep the
    /// source formula exactly.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] when `k` is not an index into
    /// `trace.peaks`. The source indexes the vector without a check, which is
    /// undefined behaviour; every call in `FeatureFinderAlgorithmPicked` stays
    /// in range.
    fn compute_theoretical(&self, trace: &MassTrace, k: usize) -> Result<f64>;
}
