// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The retention-time shape model shared by the picked feature finder's trace
//! fitters (`FEATUREFINDER/TraceFitter.h`).
//!
//! The integrator transcribed the public interface of the pinned header (core
//! `bc9cc12`) as the `TraceFitter` trait and the `TraceFitterParams` record, so
//! that package B4-GAUSS (`GaussTraceFitter`) and package B5-EGH
//! (`EGHTraceFitter`) could build against one contract at the same time.
//! B4-GAUSS filled the rest of the file: the parameter defaults and their
//! `Param` mapping, the Levenberg-Marquardt driver [`optimize`] (source
//! `optimize_`), and the helpers the two fitters share. Changing a signature of
//! the trait or the fields of the parameter record needs the integrator,
//! because both fitters and the later feature-extension package depend on them.
//!
//! The API mapping, the preserved source conventions, the native differences
//! and the evidence are in `docs/TRACE_FITTER_SUPPORT.md`.
//!
//! # Source members that are not trait methods
//!
//! - `GenericFunctor`, the raw-pointer residual and Jacobian functor, becomes
//!   the pair of closures passed to [`optimize`], which hands them to
//!   [`crate::math::fitters::levenberg_marquardt::minimize`]. The functor's
//!   `inputs()` is the length of the parameter vector and `values()` the
//!   `values` argument.
//! - `ModelData`, the protected bundle of the traces and the weighting flag, is
//!   internal to each fitter's functor.
//! - `getOptimizedParameters_`, the protected hook that copies the optimised
//!   vector back, is a method of each fitter, which calls it after [`optimize`]
//!   succeeds.
//! - `optimize_` becomes the free function [`optimize`] with the signature
//!   fixed for the wave-2 contract: `optimize(x: &mut [f64], values: usize,
//!   residual, jacobian, parameters: &TraceFitterParams) -> Result<()>`, where
//!   `residual` is `FnMut(&[f64], &mut [f64])` and `jacobian` is
//!   `FnMut(&[f64], &mut DenseMatrix)`. It refuses fewer residuals than
//!   parameters, and every Levenberg-Marquardt status up to and including
//!   `ImproperInputParameters`, with the source's `UnableToFit-FinalSet`
//!   condition; it accepts every later status, including the exhausted
//!   evaluation budget. [`optimize_with_status`] is the same driver returning
//!   the accepted status.
//! - `computeTheoretical`, the source's one non-virtual member, is the required
//!   trait method [`TraceFitter::compute_theoretical`]; [`compute_theoretical`]
//!   is the shared implementation both fitters delegate to.
//! - The copy constructor and assignment operator become `Clone` on each
//!   fitter; `updateMembers_` becomes the typed parameter setter.
//!
//! # Shared helpers
//!
//! - [`TraceFitterParams::defaults`], [`TraceFitterParams::from_param`] and
//!   [`TraceFitterParams::to_param`] are the `DefaultParamHandler` surface:
//!   `getDefaults`, `setParameters` with `updateMembers_`, and `getParameters`.
//! - [`initial_shape`] is the part of `setInitialParameters_` that the Gaussian
//!   and EGH fitters share: the intensity profile, the five-point running mean,
//!   the first strict maximum and the half-maximum walk.
//! - [`unable_to_fit`] builds the error that stands in for
//!   `Exception::UnableToFit` with the name `UnableToFit-FinalSet`, and
//!   [`stream_number`] writes a number as the default `std::ostream` does, for
//!   the gnuplot formulas.
//!
//! The source is serial, and so are the fitters.
//!
//! [`optimize`]: crate::analysis::feature_finder_picked::trace_fitter::optimize
//! [`optimize_with_status`]: crate::analysis::feature_finder_picked::trace_fitter::optimize_with_status
//! [`compute_theoretical`]: crate::analysis::feature_finder_picked::trace_fitter::compute_theoretical
//! [`initial_shape`]: crate::analysis::feature_finder_picked::trace_fitter::initial_shape
//! [`unable_to_fit`]: crate::analysis::feature_finder_picked::trace_fitter::unable_to_fit
//! [`stream_number`]: crate::analysis::feature_finder_picked::trace_fitter::stream_number
//! [`TraceFitter::compute_theoretical`]: crate::analysis::feature_finder_picked::trace_fitter::TraceFitter::compute_theoretical
//! [`TraceFitterParams::defaults`]: crate::analysis::feature_finder_picked::trace_fitter::TraceFitterParams::defaults
//! [`TraceFitterParams::from_param`]: crate::analysis::feature_finder_picked::trace_fitter::TraceFitterParams::from_param
//! [`TraceFitterParams::to_param`]: crate::analysis::feature_finder_picked::trace_fitter::TraceFitterParams::to_param

use crate::analysis::feature_finder_picked::helper_structs::{MassTrace, MassTraces};
use crate::format::file_info::text_format::ostream_g;
use crate::math::fitters::levenberg_marquardt::{
    DenseMatrix, LmParameters, LmStatus, minimize, preflight_points,
};
use crate::param::{DefaultParamHandler, Param, ParamValue};
use crate::{Error, Result};

/// The parameters every trace fitter reads: source `TraceFitter`'s
/// `DefaultParamHandler` defaults and the members `max_iterations_` and
/// `weighted_` that `updateMembers_` refreshes from them.
///
/// The source registers `max_iteration` with the default `500` and the
/// `advanced` tag, and `weighted` as the string `"false"` with the valid
/// strings `"true"` and `"false"`, also `advanced`. [`Default`] gives those
/// values, [`Self::defaults`] the source's `Param` with descriptions, tags and
/// restrictions, and [`Self::from_param`] and [`Self::to_param`] convert
/// between the record and a `Param`. The record derives no `Default`, because a
/// derived default would give `0` and `false` instead of the source's `500` and
/// `false`; the implementation is written out.
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
    ///
    /// The record and [`Self::to_param`] take the whole `i64` range, but
    /// [`Self::from_param`] reads back only values within `i32`; see there.
    pub max_iteration: i64,
    /// Whether each trace's residuals are weighted by its theoretical
    /// intensity during the fit: source parameter `weighted`, the string
    /// `"true"` or `"false"`, default `"false"`.
    pub weighted: bool,
}

impl Default for TraceFitterParams {
    /// `max_iteration` 500 and `weighted` false: the source defaults, which a
    /// freshly constructed source fitter holds, because its constructor ends in
    /// `defaultsToParam_` and hence `updateMembers_`.
    fn default() -> Self {
        Self {
            max_iteration: Self::DEFAULT_MAX_ITERATION,
            weighted: false,
        }
    }
}

impl TraceFitterParams {
    /// The source default of `max_iteration`.
    pub const DEFAULT_MAX_ITERATION: i64 = 500;

    /// The `DefaultParamHandler` name of the source class, used in parameter
    /// warnings and errors.
    pub const HANDLER_NAME: &'static str = "TraceFitter";

    /// Description of `max_iteration`, verbatim from `TraceFitter.cpp`.
    pub const MAX_ITERATION_DESCRIPTION: &'static str =
        "Maximum number of iterations used by the Levenberg-Marquardt algorithm.";

    /// Description of `weighted`, verbatim from `TraceFitter.cpp`.
    pub const WEIGHTED_DESCRIPTION: &'static str =
        "Weight mass traces according to their theoretical intensities.";

    /// The source defaults as a `Param`: source `getDefaults`.
    ///
    /// `max_iteration` is the integer `500` and `weighted` the string
    /// `"false"` restricted to `"true"` and `"false"`; both carry the source's
    /// description and the `advanced` tag. The source sets no numeric minimum
    /// or maximum on `max_iteration`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] only if the `Param` resource ceilings
    /// were smaller than these two entries, which they are not.
    pub fn defaults() -> Result<Param> {
        let advanced = [String::from("advanced")];
        let mut param = Param::new();
        param.set_value(
            "max_iteration",
            ParamValue::Integer(Self::DEFAULT_MAX_ITERATION),
            Self::MAX_ITERATION_DESCRIPTION,
            &advanced,
        )?;
        param.set_value(
            "weighted",
            ParamValue::String("false".into()),
            Self::WEIGHTED_DESCRIPTION,
            &advanced,
        )?;
        param.set_valid_strings("weighted", &[String::from("true"), String::from("false")])?;
        Ok(param)
    }

    /// The record for `parameters`: source `setParameters` followed by
    /// `updateMembers_`.
    ///
    /// Missing entries take their defaults, as `Param::setDefaults` fills
    /// them. The staged tree is then checked against [`Self::defaults`] as
    /// `Param::checkDefaults` checks it: an unknown key is a warning, returned
    /// in the vector, and a value of the wrong type or outside its
    /// restrictions is an error. Finally `max_iteration` is read as an integer
    /// and `weighted` is true exactly when it is the string `"true"`, as the
    /// source's `param_.getValue("weighted") == "true"`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] where the source throws
    /// `Exception::InvalidParameter`: when `max_iteration` is not an integer,
    /// or `weighted` is not one of the strings `"true"` and `"false"`. Also
    /// returns it when a `Param` resource ceiling is reached.
    ///
    /// Also returns it, with the message `"parameter value cannot be converted
    /// to i32"`, when `max_iteration` lies outside the `i32` range, where the
    /// source accepts it. This is a native difference inherited from
    /// `crate::param`, whose restriction check converts every integer entry to
    /// `i32` and fails on overflow. The source's `ParamEntry::isValid` narrows
    /// the value to `int` without a check, and `updateMembers_` then reads the
    /// full 64-bit value. The executed product SDK accepts and stores 2^31,
    /// 3,000,000,000, -2^31 - 1 and `i64::MAX`.
    ///
    /// `FeatureFinderAlgorithmPicked` does not reach the difference. It passes
    /// `fit:max_iterations` as a `UInt`, and the source checks that parameter's
    /// minimum of 1 on the `int`-narrowed value, so the value it hands the
    /// fitter stays below 2^31.
    pub fn from_param(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(Self::HANDLER_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_checked_param)
    }

    /// The record as a `Param` with the source's descriptions, tags and
    /// restrictions: source `getParameters` of a fitter holding these values.
    ///
    /// Every `max_iteration` is written, but [`Self::from_param`] reads the
    /// result back only when `max_iteration` lies within `i32`. Outside that
    /// range the round trip fails, where the source's round trip succeeds.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] only if a `Param` resource ceiling is
    /// reached, which two entries cannot do.
    pub fn to_param(&self) -> Result<Param> {
        let mut param = Self::defaults()?;
        let advanced = [String::from("advanced")];
        param.set_value(
            "max_iteration",
            ParamValue::Integer(self.max_iteration),
            Self::MAX_ITERATION_DESCRIPTION,
            &advanced,
        )?;
        param.set_value(
            "weighted",
            ParamValue::String(if self.weighted { "true" } else { "false" }.into()),
            Self::WEIGHTED_DESCRIPTION,
            &advanced,
        )?;
        param.set_valid_strings("weighted", &[String::from("true"), String::from("false")])?;
        Ok(param)
    }

    /// `updateMembers_` on a tree that already passed the default check.
    fn from_checked_param(parameters: &Param) -> Result<Self> {
        let max_iteration = parameters.value("max_iteration")?.to_i64()?;
        let weighted = *parameters.value("weighted")? == ParamValue::String("true".into());
        Ok(Self {
            max_iteration,
            weighted,
        })
    }
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
///
/// The trait provides no method, so an implementor that leaves one out does not
/// compile. That is the Rust form of the source's pure virtual members, which
/// `TraceFitter_test` probes with a subclass whose every override throws
/// `Exception::NotImplemented`:
///
/// ```compile_fail,E0046
/// use openms::analysis::feature_finder_picked::helper_structs::MassTraces;
/// use openms::analysis::feature_finder_picked::trace_fitter::{TraceFitter, TraceFitterParams};
///
/// struct Incomplete(TraceFitterParams);
///
/// impl TraceFitter for Incomplete {
///     fn parameters(&self) -> &TraceFitterParams {
///         &self.0
///     }
///     fn set_parameters(&mut self, parameters: TraceFitterParams) {
///         self.0 = parameters;
///     }
///     fn fit(&mut self, _traces: &MassTraces) -> openms::Result<()> {
///         Ok(())
///     }
/// }
/// ```
pub trait TraceFitter {
    /// The parameters the next [`Self::fit`] uses: the typed counterpart of the
    /// source's inherited `DefaultParamHandler::getParameters`.
    ///
    /// [`TraceFitterParams::to_param`] gives the `Param` form.
    fn parameters(&self) -> &TraceFitterParams;

    /// Replace the parameters: the typed counterpart of the source's
    /// `DefaultParamHandler::setParameters` followed by `updateMembers_`.
    ///
    /// Validation of textual parameter values, which the source performs
    /// inside `setParameters`, happens where a `Param` is converted into a
    /// [`TraceFitterParams`], in [`TraceFitterParams::from_param`], so this
    /// setter cannot fail.
    fn set_parameters(&mut self, parameters: TraceFitterParams);

    /// Fit the model to `traces` with Levenberg-Marquardt: source pure virtual
    /// `fit`.
    ///
    /// The implementation derives initial parameters from the traces, runs the
    /// shared driver [`optimize`] with the evaluation budget
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
    /// propagate unchanged, as do [`optimize`]'s resource ceilings. Each
    /// implementation documents what its fitted parameters hold after an error.
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
    /// does, precision 6 in `%g` style ([`stream_number`]). The source declares
    /// it non-const; it reads state only.
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
    /// source formula exactly, most simply by delegating to the shared
    /// [`compute_theoretical`].
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] when `k` is not an index into
    /// `trace.peaks`. The source indexes the vector without a check, which is
    /// undefined behaviour; every call in `FeatureFinderAlgorithmPicked` stays
    /// in range.
    fn compute_theoretical(&self, trace: &MassTrace, k: usize) -> Result<f64>;
}

/// The name of the source's `Exception::UnableToFit` thrown by `optimize_`,
/// which starts the message of every [`unable_to_fit`] error.
pub const UNABLE_TO_FIT_FINAL_SET: &str = "UnableToFit-FinalSet";

/// Most residual evaluations times residuals that one [`optimize`] call may
/// spend, 2^30.
///
/// The source bounds a fit only by `max_iteration`, which has no maximum, so a
/// fit that does not terminate on its own could run for as many evaluations as
/// the parameter allows. [`optimize`] therefore passes the smaller of
/// `max_iteration` and `MAX_RESIDUAL_WORK / values` (at least 1) to the solver
/// and refuses, rather than accepts, a fit that exhausts this ceiling before
/// the configured budget.
///
/// A fit that ends before the ceiling is unaffected: its evaluations, status
/// and parameters are those of the uncapped budget. So is a fit that ends at
/// exactly the ceiling's evaluation with status 1, 2 or 3, because the solver
/// tests those before the budget. A fit whose uncapped run would end at exactly
/// that evaluation with status 4 (`CosinusTooSmall`), 6, 7 or 8 is refused. The
/// solver tests the budget first: before `FtolTooSmall`, `XtolTooSmall` and
/// `GtolTooSmall` in the same step, and before the next iteration's
/// `CosinusTooSmall` test, as Eigen does. This needs at least
/// `MAX_RESIDUAL_WORK / values` evaluations, 1,073 or more within the solver's
/// point ceiling. At the feature finder's sizes (tens to a few thousand
/// residuals) the ceiling is hundreds of thousands of evaluations or more,
/// against natural termination after tens.
pub const MAX_RESIDUAL_WORK: usize = 1 << 30;

/// The error standing in for the source's `Exception::UnableToFit` named
/// `UnableToFit-FinalSet`: [`Error::InvalidValue`] with the message
/// `"UnableToFit-FinalSet: <message>"`.
///
/// `optimize_` throws it with the messages `"Skipping feature, we always
/// expect N>=p"` and `"Could not fit the gaussian to the data: Error <status>"`,
/// the latter also for the EGH model.
pub fn unable_to_fit(message: &str) -> Error {
    Error::InvalidValue(format!("{UNABLE_TO_FIT_FINAL_SET}: {message}"))
}

/// The message `optimize_` throws with when there are fewer residuals than
/// parameters, verbatim from `TraceFitter.cpp`. The fitters also use it for
/// traces without any peak, which the source cannot reach without undefined
/// behaviour.
pub const FEWER_RESIDUALS_THAN_PARAMETERS: &str = "Skipping feature, we always expect N>=p";

/// Run Levenberg-Marquardt from `x` on analytic residuals and Jacobian: source
/// `TraceFitter::optimize_`.
///
/// `x` holds the initial parameters on entry and the optimised parameters after
/// success. `values` is the number of residuals, the source functor's
/// `values()`, and `x.len()` its `inputs()`. `residual(x, fvec)` fills the
/// `values` residuals at `x`; `jacobian(x, jac)` fills the `values`-by-`x.len()`
/// Jacobian, row `i` and column `j` being the derivative of residual `i` with
/// respect to parameter `j`. The Jacobian is analytic, so each call counts as a
/// Jacobian evaluation and none of the budget, as Eigen counts a functor whose
/// `df` returns 0.
///
/// The solver is [`crate::math::fitters::levenberg_marquardt::minimize`] with
/// Eigen's defaults (`factor` 100, `ftol` and `xtol` `sqrt(f64::EPSILON)`,
/// `gtol` 0) and `max_fev` set to [`TraceFitterParams::max_iteration`], which
/// is what the source sets `lmSolver.parameters.maxfev` to. The weighting flag
/// of `parameters` is not read here: it belongs to the residual functor.
///
/// The configuration is the source's, but the solver is not yet bit-faithful
/// to the executed Eigen 5.0.1. On most executed inputs the transcription's
/// path departs from Eigen's in the last bits at the first trial step, so
/// fitted parameters can differ well beyond 1e-9 and the status and the number
/// of evaluations can differ too. The agreement recorded for the class-test
/// and FeatureFinderCentroided_1 fits holds for those fixtures only. The root
/// cause is in `src/math/fitters/levenberg_marquardt.rs`, under investigation
/// in lane B3b; `docs/TRACE_FITTER_SUPPORT.md` ("Known gap") has the
/// measurements.
///
/// The source believes, after reading Eigen, that every status except
/// `NotStarted`, `Running` and `ImproperInputParameters` is a good termination;
/// so an exhausted budget (`TooManyFunctionEvaluation`) is accepted with the
/// parameters reached, exactly as a converged fit.
///
/// # Errors
///
/// Every error leaves `x` unchanged, as the source copies the solver's vector
/// back only after its checks.
///
/// - [`unable_to_fit`] with `"Skipping feature, we always expect N>=p"` when
///   `values < x.len()`, before anything is evaluated.
/// - [`unable_to_fit`] with `"Could not fit the gaussian to the data: Error 0"`
///   when the solver reports `ImproperInputParameters`: when `x` is empty or
///   `max_iteration <= 0`. The source's message says "gaussian" whichever
///   model is fitted.
/// - [`Error::InvalidValue`] from
///   [`crate::math::fitters::levenberg_marquardt::preflight_points`] when
///   `values` or the dense Jacobian exceed the solver's ceilings, checked
///   before anything is evaluated. The source has no ceiling.
/// - [`Error::InvalidValue`] when the solver stops on the work ceiling derived
///   from [`MAX_RESIDUAL_WORK`] while the configured budget is larger. That
///   includes a fit whose natural end with status 4, 6, 7 or 8 falls on exactly
///   the ceiling's evaluation, because the budget test comes first; see that
///   constant. The source would continue or accept.
pub fn optimize<R, J>(
    x: &mut [f64],
    values: usize,
    residual: R,
    jacobian: J,
    parameters: &TraceFitterParams,
) -> Result<()>
where
    R: FnMut(&[f64], &mut [f64]),
    J: FnMut(&[f64], &mut DenseMatrix),
{
    optimize_with_status(x, values, residual, jacobian, parameters).map(|_| ())
}

/// [`optimize`], returning the accepted Levenberg-Marquardt status.
///
/// The source discards the status once it has decided to accept it; this
/// returns it, so that a caller or a test can tell a converged fit from one
/// that spent its budget. The status is never `ImproperInputParameters`, which
/// is an error.
///
/// # Errors
///
/// As [`optimize`].
pub fn optimize_with_status<R, J>(
    x: &mut [f64],
    values: usize,
    residual: R,
    jacobian: J,
    parameters: &TraceFitterParams,
) -> Result<LmStatus>
where
    R: FnMut(&[f64], &mut [f64]),
    J: FnMut(&[f64], &mut DenseMatrix),
{
    optimize_bounded(x, values, residual, jacobian, parameters, MAX_RESIDUAL_WORK)
}

/// The driver with the work ceiling `residual_work` in place of
/// [`MAX_RESIDUAL_WORK`], so the ceiling can be tested at a small size.
fn optimize_bounded<R, J>(
    x: &mut [f64],
    values: usize,
    residual: R,
    mut jacobian: J,
    parameters: &TraceFitterParams,
    residual_work: usize,
) -> Result<LmStatus>
where
    R: FnMut(&[f64], &mut [f64]),
    J: FnMut(&[f64], &mut DenseMatrix),
{
    if values < x.len() {
        return Err(unable_to_fit(FEWER_RESIDUALS_THAN_PARAMETERS));
    }
    preflight_points(values, x.len())?;
    if x.is_empty() || parameters.max_iteration <= 0 {
        // Eigen's `minimizeInit` refuses `n <= 0` and `maxfev <= 0` before it
        // evaluates anything; the source then throws with that status.
        return Err(improper_status(LmStatus::ImproperInputParameters));
    }
    let budget = usize::try_from(parameters.max_iteration).unwrap_or(usize::MAX);
    let ceiling = (residual_work / values.max(1)).max(1);
    let max_fev = budget.min(ceiling);
    let mut working = x.to_vec();
    let status = minimize(
        &mut working,
        values,
        residual,
        |point: &[f64], jac: &mut DenseMatrix| {
            jacobian(point, jac);
            0
        },
        &LmParameters {
            max_fev,
            ..LmParameters::default()
        },
    );
    if status.code() <= LmStatus::ImproperInputParameters.code() {
        return Err(improper_status(status));
    }
    if status == LmStatus::TooManyFunctionEvaluation && max_fev < budget {
        return Err(Error::InvalidValue(format!(
            "the trace fit spent its work ceiling of {max_fev} residual evaluations over \
             {values} residuals before its configured budget of {budget}"
        )));
    }
    x.copy_from_slice(&working);
    Ok(status)
}

/// The `UnableToFit-FinalSet` error for a refused solver status.
fn improper_status(status: LmStatus) -> Error {
    unable_to_fit(&format!(
        "Could not fit the gaussian to the data: Error {}",
        status.code()
    ))
}

/// The shared source formula of `TraceFitter::computeTheoretical`:
/// `trace.theoretical_int * fitter.value(trace.peaks[k].rt)`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `k` is not an index into
/// `trace.peaks`, where the source indexes without a check.
pub fn compute_theoretical<F: TraceFitter + ?Sized>(
    fitter: &F,
    trace: &MassTrace,
    k: usize,
) -> Result<f64> {
    let peak = trace.peaks.get(k).ok_or_else(|| {
        Error::InvalidValue(format!(
            "peak index {k} is out of range for a mass trace of {} peaks",
            trace.peaks.len()
        ))
    })?;
    Ok(trace.theoretical_int * fitter.value(peak.rt))
}

/// A number as the source's default `std::ostream` writes a `double`:
/// precision 6 in `%g` style, as `getGnuplotFormula` streams it.
///
/// This is `crate::format::file_info::text_format::ostream_g` at precision 6,
/// which documents the spelling of non-finite values and the one class of
/// exact decimal ties where Apple libc differs from the C standard.
pub fn stream_number(value: f64) -> String {
    ostream_g(value, 6)
}

/// Whether [`initial_shape`] smooths the intensity profile.
///
/// The source's Gaussian and EGH `setInitialParameters_` differ only here.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ProfileSmoothing {
    /// Smooth only profiles of more than three retention times, and take a
    /// shorter profile's summed intensities unchanged: `GaussTraceFitter`,
    /// whose source guards the moving average with `N <= LEN + 1`.
    SkipShortProfiles,
    /// Smooth every profile: `EGHTraceFitter`, whose source has no guard.
    Always,
}

/// The start-value estimate that `setInitialParameters_` derives from the
/// traces before each model turns it into parameters.
///
/// Every field is computed with the source's arithmetic, in its order; see
/// [`initial_shape`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InitialShape {
    /// Index into the profile of the first strict maximum of the smoothed
    /// intensities: source `max_index`.
    pub max_index: usize,
    /// Smoothed maximum minus the traces' baseline: source `height_`.
    pub height: f64,
    /// Retention time of the profile entry at [`Self::max_index`]: source
    /// `x0_` (Gaussian) or `apex_rt_` (EGH).
    pub apex_rt: f64,
    /// Last minus first retention time of the profile: source
    /// `region_rt_span_`.
    pub region_rt_span: f64,
    /// Where the leftward half-maximum walk stopped: source `left_index`
    /// (Gaussian) or `index` (EGH).
    pub left_index: usize,
    /// Smoothed intensity at [`Self::left_index`]: source `left_height`.
    pub left_height: f64,
    /// Retention time at [`Self::left_index`]: source `left_rt`.
    pub left_rt: f64,
    /// Where the rightward half-maximum walk stopped.
    pub right_index: usize,
    /// Smoothed intensity at [`Self::right_index`]: source `right_height`.
    pub right_height: f64,
    /// Retention time at [`Self::right_index`]: source `right_rt`.
    pub right_rt: f64,
}

/// Half-width of the moving-average window: source `LEN`, for a window of
/// `2 * LEN + 1` points.
const SMOOTHING_LEN: usize = 2;

/// The shared start-value estimate of `GaussTraceFitter::setInitialParameters_`
/// and `EGHTraceFitter::setInitialParameters_`.
///
/// The steps, as both sources write them:
///
/// 1. The traces' intensity profile, [`MassTraces::intensity_profile`]: one
///    summed intensity per retention time, `N` entries.
/// 2. Smoothing with a centred five-point moving average over `totals`, the
///    profile intensities padded with two zeros at each end, kept as a running
///    sum: the sum starts as `0.0 + totals[2] + totals[3]`, and for each `i`
///    adds `totals[i + 4]`, stores the sum divided by 5 and then subtracts
///    `totals[i]`. With [`ProfileSmoothing::SkipShortProfiles`] a profile of
///    `N <= 3` entries is taken unsmoothed instead. The padding is read through
///    an index function rather than allocated; the values and the order of the
///    additions are the source's.
/// 3. The first strict maximum of the smoothed values (`>`), so the earliest
///    of equal maxima wins and a NaN never replaces it.
/// 4. `height` is that maximum minus [`MassTraces::baseline`], `apex_rt` its
///    retention time and `region_rt_span` the last minus the first profile
///    retention time.
/// 5. From the maximum, walk left while the index is above 0 and the smoothed
///    value is greater than `height * 0.5`, and likewise right while the index
///    is below `N - 1`; the stops give the half-maximum heights and retention
///    times.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the profile is empty, which happens
/// exactly when the traces hold no peak. The source reads the first element of
/// the empty profile, which is undefined behaviour; both fitters refuse such
/// traces before calling this. Errors of [`MassTraces::intensity_profile`]
/// propagate unchanged.
pub fn initial_shape(traces: &MassTraces, smoothing: ProfileSmoothing) -> Result<InitialShape> {
    let profile = traces.intensity_profile()?;
    let n = profile.len();
    let (Some(first), Some(last)) = (profile.first(), profile.last()) else {
        return Err(Error::InvalidValue(
            "an empty intensity profile has no initial shape".into(),
        ));
    };
    // `totals[i]` of the source: the profile intensities padded with LEN zeros.
    let totals = |index: usize| -> f64 {
        if index < SMOOTHING_LEN || index >= n + SMOOTHING_LEN {
            0.0
        } else {
            profile[index - SMOOTHING_LEN].1
        }
    };
    let mut smoothed = Vec::with_capacity(n);
    let mut max_index = 0usize;
    if smoothing == ProfileSmoothing::SkipShortProfiles && n <= SMOOTHING_LEN + 1 {
        for i in 0..n {
            smoothed.push(totals(i + SMOOTHING_LEN));
            if smoothed[i] > smoothed[max_index] {
                max_index = i;
            }
        }
    } else {
        let window = (2 * SMOOTHING_LEN + 1) as f64;
        let mut sum = 0.0;
        for index in SMOOTHING_LEN..2 * SMOOTHING_LEN {
            sum += totals(index);
        }
        for i in 0..n {
            sum += totals(i + 2 * SMOOTHING_LEN);
            smoothed.push(sum / window);
            sum -= totals(i);
            if smoothed[i] > smoothed[max_index] {
                max_index = i;
            }
        }
    }
    let height = smoothed[max_index] - traces.baseline;
    let apex_rt = profile[max_index].0;
    let region_rt_span = last.0 - first.0;
    let half = height * 0.5;
    let mut left_index = max_index;
    while left_index > 0 && smoothed[left_index] > half {
        left_index -= 1;
    }
    let mut right_index = max_index;
    while right_index < n - 1 && smoothed[right_index] > half {
        right_index += 1;
    }
    Ok(InitialShape {
        max_index,
        height,
        apex_rt,
        region_rt_span,
        left_index,
        left_height: smoothed[left_index],
        left_rt: profile[left_index].0,
        right_index,
        right_height: smoothed[right_index],
        right_rt: profile[right_index].0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const TIMES: [f64; 5] = [0.0, 1.0, 2.0, 3.0, 4.0];
    const DATA: [f64; 5] = [1.1, 2.6, 7.7, 19.5, 55.3];
    const START: [f64; 2] = [0.2, 0.1];

    /// `a * exp(b * t) - y` over five points that no exponential passes
    /// through, so the fit ends with a non-zero residual after some evaluations.
    fn exponential(x: &[f64], f: &mut [f64]) {
        for ((slot, t), y) in f.iter_mut().zip(TIMES).zip(DATA) {
            *slot = x[0] * libm::exp(x[1] * t) - y;
        }
    }

    fn exponential_jacobian(x: &[f64], jac: &mut DenseMatrix) {
        for (row, t) in TIMES.into_iter().enumerate() {
            let e = libm::exp(x[1] * t);
            jac.set(row, 0, e);
            jac.set(row, 1, x[0] * t * e);
        }
    }

    fn run(max_iteration: i64, residual_work: usize) -> (Result<LmStatus>, [f64; 2], usize) {
        let mut x = START;
        let mut evaluations = 0usize;
        let result = optimize_bounded(
            &mut x,
            5,
            |point: &[f64], f: &mut [f64]| {
                evaluations += 1;
                exponential(point, f);
            },
            exponential_jacobian,
            &TraceFitterParams {
                max_iteration,
                weighted: false,
            },
            residual_work,
        );
        (result, x, evaluations)
    }

    /// A fit that ends within the work ceiling is the uncapped fit; one that
    /// would need more is refused and leaves the parameters unchanged, unless
    /// the configured budget is the ceiling itself, which the source accepts.
    #[test]
    fn the_work_ceiling_refuses_only_fits_that_outlast_it() {
        let (natural, fitted, evaluations) = run(500, MAX_RESIDUAL_WORK);
        let natural = natural.unwrap();
        // Termination inside the trial loop, where the ftol and xtol tests
        // precede the budget test.
        assert!((1..=3).contains(&natural.code()), "{natural:?}");
        assert!(evaluations > 3, "{evaluations}");

        // Ceiling exactly at natural termination: the same fit.
        let (status, x, _) = run(500, 5 * evaluations);
        assert_eq!(status.unwrap(), natural);
        assert_eq!(x.map(f64::to_bits), fitted.map(f64::to_bits));

        // One evaluation short with a larger budget: refused, x unchanged.
        let (status, x, spent) = run(500, 5 * (evaluations - 1));
        assert!(matches!(status, Err(Error::InvalidValue(ref m)) if m.contains("work ceiling")));
        assert_eq!(x, START);
        assert_eq!(spent, evaluations - 1);

        // The same cap as the configured budget: accepted, as the source does.
        let budget = i64::try_from(evaluations - 1).unwrap();
        let (status, capped, _) = run(budget, 5 * (evaluations - 1));
        assert_eq!(status.unwrap(), LmStatus::TooManyFunctionEvaluation);
        let (uncapped_status, uncapped, _) = run(budget, MAX_RESIDUAL_WORK);
        assert_eq!(
            uncapped_status.unwrap(),
            LmStatus::TooManyFunctionEvaluation
        );
        assert_eq!(capped.map(f64::to_bits), uncapped.map(f64::to_bits));
    }

    /// `x - 1` and a constant zero: one Gauss-Newton step reaches a zero
    /// residual, and the next iteration stops with `CosinusTooSmall` after two
    /// evaluations.
    fn linear(x: &[f64], f: &mut [f64]) {
        f[0] = x[0] - 1.0;
        f[1] = 0.0;
    }

    fn linear_jacobian(_x: &[f64], jac: &mut DenseMatrix) {
        jac.set(0, 0, 1.0);
        jac.set(1, 0, 0.0);
    }

    fn run_linear(max_iteration: i64, residual_work: usize) -> (Result<LmStatus>, [f64; 1], usize) {
        let mut x = [0.0];
        let mut evaluations = 0usize;
        let result = optimize_bounded(
            &mut x,
            2,
            |point: &[f64], f: &mut [f64]| {
                evaluations += 1;
                linear(point, f);
            },
            linear_jacobian,
            &TraceFitterParams {
                max_iteration,
                weighted: false,
            },
            residual_work,
        );
        (result, x, evaluations)
    }

    /// The budget test precedes `CosinusTooSmall` (and `FtolTooSmall`,
    /// `XtolTooSmall`, `GtolTooSmall`), so a fit whose natural end with one of
    /// those statuses falls on exactly the ceiling's evaluation is refused, as
    /// `MAX_RESIDUAL_WORK` documents; one evaluation more of ceiling gives the
    /// uncapped fit.
    #[test]
    fn a_late_status_at_exactly_the_ceiling_is_refused() {
        let (natural, fitted, evaluations) = run_linear(500, MAX_RESIDUAL_WORK);
        assert_eq!(natural.unwrap(), LmStatus::CosinusTooSmall);
        assert_eq!(evaluations, 2);
        assert_eq!(fitted, [1.0]);

        let (status, x, _) = run_linear(500, 2 * 3);
        assert_eq!(status.unwrap(), LmStatus::CosinusTooSmall);
        assert_eq!(x.map(f64::to_bits), fitted.map(f64::to_bits));

        let (status, x, spent) = run_linear(500, 2 * 2);
        assert!(matches!(status, Err(Error::InvalidValue(ref m)) if m.contains("work ceiling")));
        assert_eq!(x, [0.0]);
        assert_eq!(spent, 2);

        // The same stop from the configured budget is accepted, as the source
        // accepts `TooManyFunctionEvaluation`.
        let (status, x, _) = run_linear(2, MAX_RESIDUAL_WORK);
        assert_eq!(status.unwrap(), LmStatus::TooManyFunctionEvaluation);
        assert_eq!(x.map(f64::to_bits), fitted.map(f64::to_bits));
    }

    /// A ceiling below one evaluation per residual still allows the initial
    /// evaluation and one trial, and then refuses the fit.
    #[test]
    fn a_ceiling_below_one_evaluation_still_allows_one() {
        let (status, x, spent) = run(500, 0);
        assert!(matches!(status, Err(Error::InvalidValue(_))));
        assert_eq!(x, START);
        assert_eq!(spent, 2);
    }
}
