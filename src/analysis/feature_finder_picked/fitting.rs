// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Trace fitting, cropping, quality checks and feature creation of the picked
//! feature finder (`FEATUREFINDER/FeatureFinderAlgorithmPicked.h`, steps 3.3.2
//! to 3.3.5 of `run_`).
//!
//! After the seed's mass traces are extended
//! ([`crate::analysis::feature_finder_picked::extension`]) the candidate goes
//! through four stages, all here:
//!
//! - **fit** — [`FittedModel`](crate::analysis::feature_finder_picked::fitting::FittedModel) selects the retention-time model from
//!   `feature:rt_shape` (source `chooseTraceFitter_`) and fits it to the
//!   traces;
//! - **crop** — [`crop_feature`](crate::analysis::feature_finder_picked::fitting::crop_feature) keeps the peaks inside the model's
//!   retention-time bounds and drops traces that the model describes badly
//!   (source `cropFeature_`);
//! - **check** — [`check_feature_quality`](crate::analysis::feature_finder_picked::fitting::check_feature_quality) applies the five acceptance rules in
//!   source order (source `checkFeatureQuality_`);
//! - **create** — [`build_feature`](crate::analysis::feature_finder_picked::fitting::build_feature) turns the surviving traces into a
//!   [`Feature`](crate::kernel::Feature) (source step 3.3.5).
//!
//! The reasons a candidate is rejected are the source's verbatim strings; they
//! are the keys the source counts in `aborts_` and this port returns in
//! [`crate::analysis::feature_finder_picked::algorithm::RunOutput::aborts`].
//!
//! The module is serial: the source parallelises one level above it, over
//! seeds. See `docs/FEATURE_FINDER_PICKED_SUPPORT.md`.

use crate::analysis::feature_finder_picked::algorithm::{ReportedMz, RtShape, Settings};
use crate::analysis::feature_finder_picked::debug::{LogSink, NoLog, g, put_all};
use crate::analysis::feature_finder_picked::egh_trace_fitter::EGHTraceFitter;
use crate::analysis::feature_finder_picked::gauss_trace_fitter::GaussTraceFitter;
use crate::analysis::feature_finder_picked::helper_structs::{
    IsotopePattern, MassTrace, MassTraces,
};
use crate::analysis::feature_finder_picked::seeds::IsotopeWindows;
use crate::analysis::feature_finder_picked::trace_fitter::{TraceFitter, TraceFitterParams};
use crate::concept::constants::PROTON_MASS_U;
use crate::concept::constants::user_param::NUM_OF_DATAPOINTS;
use crate::kernel::Feature;
use crate::math::statistic_functions::pearson_correlation_coefficient;
use crate::metadata::MetaValue;
use crate::{Error, Result};

/// Abort reason of a seed whose best isotope fit stays below
/// `feature:min_isotope_fit`, verbatim from the source.
pub const ABORT_NO_ISOTOPE_PATTERN: &str =
    "Could not find good enough isotope pattern containing the seed";

/// Abort reason of a seed whose extended traces do not describe it any more,
/// verbatim from the source.
pub const ABORT_COULD_NOT_EXTEND: &str = "Could not extend seed";

/// Abort reason of a fit wider than `feature:max_rt_span`, verbatim from the
/// source.
pub const ABORT_MAX_RT_SPAN: &str = "Invalid fit: Fitted model is bigger than 'max_rt_span'";

/// Abort reason of a cropped candidate with too few traces or peaks, verbatim
/// from the source.
pub const ABORT_TOO_FEW_TRACES: &str = "Invalid feature after fit - too few traces or peaks left";

/// Abort reason of a fit whose centre lies outside the cropped traces, verbatim
/// from the source.
pub const ABORT_CENTER_OUTSIDE: &str = "Invalid fit: Center outside of feature bounds";

/// Abort reason of a candidate narrower than `feature:min_rt_span`, verbatim
/// from the source.
pub const ABORT_MIN_RT_SPAN: &str = "Invalid fit: Less than 'min_rt_span' left after fit";

/// Abort reason of a candidate below `feature:min_score`, verbatim from the
/// source.
pub const ABORT_QUALITY_TOO_LOW: &str = "Feature quality too low after fit";

/// `std::max(0.0, value)` as the source spells it: `(0.0 < value) ? value :
/// 0.0`, so a NaN gives `0.0` and `-0.0` gives `0.0`.
fn max0(value: f64) -> f64 {
    if 0.0 < value { value } else { 0.0 }
}

/// The retention-time model of one feature candidate: source
/// `chooseTraceFitter_`.
///
/// `feature:rt_shape = symmetric` selects the Gaussian and `asymmetric` the
/// exponential-Gaussian hybrid. The source returns a `TraceFitter` pointer and
/// signals the asymmetric choice by setting its `tau` output to `-1.0`, which
/// it later tests with `egh_tau != 0.0` before casting the pointer back to
/// `EGHTraceFitter` for the three `EGH_*` meta values. This enum carries the
/// choice in the type instead, so no cast is needed; [`Self::egh`] is the cast.
#[derive(Clone, Debug, PartialEq)]
pub enum FittedModel {
    /// `symmetric`: the Gaussian model, the default.
    Gauss(GaussTraceFitter),
    /// `asymmetric`: the exponential-Gaussian hybrid.
    Egh(EGHTraceFitter),
}

impl FittedModel {
    /// The model `shape` selects, holding `parameters`: source
    /// `chooseTraceFitter_` followed by `fitter->setParameters(...)`.
    pub fn new(shape: RtShape, parameters: TraceFitterParams) -> Self {
        match shape {
            RtShape::Symmetric => Self::Gauss(GaussTraceFitter::with_parameters(parameters)),
            RtShape::Asymmetric => Self::Egh(EGHTraceFitter::with_parameters(parameters)),
        }
    }

    /// The model as the trait the algorithm calls it through.
    pub fn as_fitter(&self) -> &dyn TraceFitter {
        match self {
            Self::Gauss(fitter) => fitter,
            Self::Egh(fitter) => fitter,
        }
    }

    /// The exponential-Gaussian hybrid, or `None` for the Gaussian: source
    /// `std::dynamic_pointer_cast<EGHTraceFitter>(fitter)`, which the source
    /// performs only when its `tau` flag is non-zero.
    pub fn egh(&self) -> Option<&EGHTraceFitter> {
        match self {
            Self::Gauss(_) => None,
            Self::Egh(fitter) => Some(fitter),
        }
    }

    /// Fit the model to `traces`: source `fitter->fit(traces)`.
    ///
    /// # `Exception::UnableToFit` is unreachable from the seed loop
    ///
    /// `TraceFitter::optimize_` throws in two places, and the source's seed
    /// loop (`FeatureFinderAlgorithmPicked.cpp:595-670`) catches neither; the
    /// exception would leave the OpenMP region and terminate the process. No
    /// input reaches either:
    ///
    /// - `TraceFitter.cpp:111`, fewer residuals than parameters. The residuals
    ///   are the peaks of the traces the loop fits, and the loop only fits
    ///   traces that pass `MassTraces::isValid` (`:638`), so there are at
    ///   least two. `extendMassTraces_` (`:1381-1479`) appends a trace with
    ///   fewer than three peaks (`MassTrace::isValid`) only at pattern index 0
    ///   while the maximum trace, which has at least three, is not yet
    ///   appended (`MassTraces::max_trace` is still 0 then and `p ==
    ///   max_trace`); every later short trace clears the list or ends it. So at
    ///   most one trace is short, and it holds at least its start peak: at
    ///   least `1 + 3 = 4` residuals, and the Gaussian has 3 parameters, the
    ///   EGH model 4.
    /// - `TraceFitter.cpp:129`, a solver status up to
    ///   `ImproperInputParameters`. Eigen's `LevenbergMarquardt::minimize`
    ///   returns that status only from `minimizeInit`, for `n <= 0`, `m < n`,
    ///   a negative tolerance, `maxfev <= 0` or `factor <= 0` (Eigen
    ///   `NonLinearOptimization/LevenbergMarquardt.h`, the install's
    ///   `deps/include/eigen3`); `minimizeOneStep` returns only `Running` or a
    ///   positive status. Here `n` is 3 or 4, `m >= n` by the first check, the
    ///   tolerances and `factor` are Eigen's defaults, and `maxfev` is
    ///   `fit:max_iterations`, which `DefaultParamHandler::setParameters`
    ///   restricts to at least 1 (`:89`): `ParamValue::operator int` and
    ///   `operator unsigned int` keep the same low 32 bits, so a value that
    ///   passes the restriction reaches the solver unchanged and positive.
    ///
    /// The residual count the source checks is `int`
    /// (`static_cast<int>(getPeakCount())`, `GaussTraceFitter.cpp:140`, and
    /// `EGHTraceFitter.cpp:29` through the `int` constructor): traces with more
    /// than `INT_MAX` peaks would also throw at `:111`, but the solver's
    /// `MAX_POINTS` ceiling refuses such traces first.
    ///
    /// # Errors
    ///
    /// As [`TraceFitter::fit`] of the selected model. From the seed loop, that
    /// is one of the port's resource ceilings (the solver's point, byte and
    /// work ceilings), which the source does not have, or the start point's
    /// refusal of a NaN retention time in the intensity profile
    /// ([`MassTraces::intensity_profile`]), where the source loops forever; the
    /// algorithm returns such an error instead of recording an abort reason.
    pub fn fit(&mut self, traces: &MassTraces) -> Result<()> {
        match self {
            Self::Gauss(fitter) => fitter.fit(traces),
            Self::Egh(fitter) => fitter.fit(traces),
        }
    }
}

/// Keep the peaks inside the fitted model's retention-time bounds and drop
/// badly described traces: source `cropFeature_`.
///
/// A trace keeps the peaks whose retention time lies within
/// `[lower_rt_bound, upper_rt_bound]`, inclusive. Its score is
/// `sqrt(correlation * max(0, 1 - mean relative deviation))`, where the
/// deviation of one peak is `|real - theoretical| / theoretical` with the
/// theoretical intensity `baseline + theoretical_int * model(rt)` and the
/// correlation is the Pearson correlation of those two series, floored at zero.
/// A trace with fewer than three remaining peaks or a score below
/// `feature:min_trace_score` is *bad*, and the source then decides by position:
///
/// - before [`MassTraces::max_trace`]: every trace kept so far is discarded and
///   the bad trace is skipped;
/// - at `max_trace`: everything is discarded and the cropping stops, which
///   makes the candidate fail the next check;
/// - after `max_trace`: the cropping stops, keeping what was collected.
///
/// The result's [`MassTraces::max_trace`] is the position of the surviving
/// maximum trace, and its [`MassTraces::baseline`] is copied from `traces`, as
/// the source copies it last.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] from
/// [`TraceFitter::compute_theoretical`] and from
/// [`pearson_correlation_coefficient`], and the ceilings of
/// [`MassTraces::reserve`].
pub fn crop_feature(
    fitter: &dyn TraceFitter,
    traces: &MassTraces,
    min_trace_score: f64,
) -> Result<MassTraces> {
    crop_feature_logged(fitter, traces, min_trace_score, &mut NoLog)
}

/// [`crop_feature`] writing the source's debug lines to `log`.
///
/// The per-trace line prints the correlation where it says `final score`, as
/// the source does.
pub(crate) fn crop_feature_logged<L: LogSink>(
    fitter: &dyn TraceFitter,
    traces: &MassTraces,
    min_trace_score: f64,
    log: &mut L,
) -> Result<MassTraces> {
    let low_bound = fitter.lower_rt_bound();
    let high_bound = fitter.upper_rt_bound();
    if log.enabled() {
        put_all(
            log,
            &[
                "    => RT bounds: ",
                &g(low_bound),
                " - ",
                &g(high_bound),
                "\n",
            ],
        );
    }
    let mut new_traces = MassTraces::new();
    let mut theoretical: Vec<f64> = Vec::new();
    let mut real: Vec<f64> = Vec::new();
    for t in 0..traces.len() {
        let trace = &traces[t];
        if log.enabled() {
            put_all(
                log,
                &[
                    "   - Trace ",
                    &t.to_string(),
                    ": (",
                    &g(trace.theoretical_int),
                    ")\n",
                ],
            );
        }
        let mut new_trace = MassTrace::default();
        let mut deviation = 0.0;
        theoretical.clear();
        real.clear();
        for k in 0..trace.peaks.len() {
            let peak = trace.peaks[k];
            if peak.rt >= low_bound && peak.rt <= high_bound {
                new_trace.peaks.push(peak);
                let theo = traces.baseline + fitter.compute_theoretical(trace, k)?;
                theoretical.push(theo);
                let measured = f64::from(peak.intensity);
                real.push(measured);
                deviation += (measured - theo).abs() / theo;
            }
        }
        let mut fit_score = 0.0;
        let mut correlation = 0.0;
        let mut final_score = 0.0;
        if !new_trace.peaks.is_empty() {
            fit_score = deviation / new_trace.peaks.len() as f64;
            correlation = max0(pearson_correlation_coefficient(&theoretical, &real)?);
            final_score = (correlation * max0(1.0 - fit_score)).sqrt();
        }
        if log.enabled() {
            put_all(
                log,
                &[
                    "     - peaks: ",
                    &new_trace.peaks.len().to_string(),
                    " / ",
                    &trace.peaks.len().to_string(),
                    " - relative deviation: ",
                    &g(fit_score),
                    " - correlation: ",
                    &g(correlation),
                    " - final score: ",
                    &g(correlation),
                    "\n",
                ],
            );
        }
        if !new_trace.is_valid() || final_score < min_trace_score {
            if t < traces.max_trace {
                new_traces = MassTraces::new();
                put_all(
                    log,
                    &["     - removed this and previous traces due to bad fit\n"],
                );
                continue;
            } else if t == traces.max_trace {
                new_traces = MassTraces::new();
                put_all(log, &["     - aborting (max trace was removed)\n"]);
                break;
            }
            put_all(
                log,
                &["     - removed due to bad fit => omitting the rest\n"],
            );
            break;
        }
        new_trace.theoretical_int = trace.theoretical_int;
        new_traces.push(new_trace);
        if t == traces.max_trace {
            new_traces.max_trace = new_traces.len() - 1;
        }
    }
    new_traces.baseline = traces.baseline;
    Ok(new_traces)
}

/// The quality of an accepted candidate: source `checkFeatureQuality_`'s three
/// output parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FeatureQuality {
    /// `max(0, 1 - mean relative deviation)` over every peak of the cropped
    /// traces; written as the meta value `score_fit`.
    pub fit_score: f64,
    /// Pearson correlation of the theoretical and measured intensities, floored
    /// at zero; written as the meta value `score_correlation`.
    pub correlation: f64,
    /// `sqrt(correlation * fit_score)`, the feature's overall quality.
    pub final_score: f64,
}

/// Whether a candidate passed [`check_feature_quality`].
#[derive(Clone, Debug, PartialEq)]
pub enum QualityOutcome {
    /// The candidate becomes a feature.
    Accepted(FeatureQuality),
    /// The candidate is dropped with this source abort reason, one of
    /// [`ABORT_MAX_RT_SPAN`], [`ABORT_TOO_FEW_TRACES`],
    /// [`ABORT_CENTER_OUTSIDE`], [`ABORT_MIN_RT_SPAN`] and
    /// [`ABORT_QUALITY_TOO_LOW`].
    Rejected(&'static str),
}

/// The five acceptance rules of a fitted candidate, in source order: source
/// `checkFeatureQuality_`.
///
/// 1. the model must not be wider than `feature:max_rt_span` of the region it
///    was initialised from ([`TraceFitter::check_maximal_rt_span`]);
/// 2. the cropped traces must still describe the seed
///    ([`MassTraces::is_valid`] against `mass_trace:mz_tolerance`);
/// 3. the model's centre must lie within the cropped traces' retention-time
///    bounds;
/// 4. those bounds must cover at least `feature:min_rt_span` of the model's
///    width ([`TraceFitter::check_minimal_rt_span`]);
/// 5. `sqrt(correlation * fit_score)` over every remaining peak must reach
///    `feature:min_score`.
///
/// Rule 5's deviation is `|real - theoretical| / theoretical` per peak, summed
/// and divided by the total peak count, then subtracted from one and floored at
/// zero; the correlation is floored at zero the same way. A NaN from either
/// floor becomes zero, as `std::max` gives.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] from [`MassTraces::rt_bounds`],
/// [`TraceFitter::compute_theoretical`] and
/// [`pearson_correlation_coefficient`].
pub fn check_feature_quality(
    fitter: &dyn TraceFitter,
    traces: &MassTraces,
    seed_mz: f64,
    settings: &Settings,
) -> Result<QualityOutcome> {
    check_feature_quality_logged(fitter, traces, seed_mz, settings, &mut NoLog)
}

/// [`check_feature_quality`] writing the source's `Quality estimation:` block
/// to `log`, which the source writes once the first four checks pass.
pub(crate) fn check_feature_quality_logged<L: LogSink>(
    fitter: &dyn TraceFitter,
    traces: &MassTraces,
    seed_mz: f64,
    settings: &Settings,
    log: &mut L,
) -> Result<QualityOutcome> {
    let rejected = |reason| Ok(QualityOutcome::Rejected(reason));
    if fitter.check_maximal_rt_span(settings.max_rt_span) {
        return rejected(ABORT_MAX_RT_SPAN);
    }
    if !traces.is_valid(seed_mz, settings.trace_tolerance) {
        return rejected(ABORT_TOO_FEW_TRACES);
    }
    let rt_bounds = traces.rt_bounds()?;
    if fitter.center() < rt_bounds.0 || fitter.center() > rt_bounds.1 {
        return rejected(ABORT_CENTER_OUTSIDE);
    }
    // The source reads the bounds a second time for the next check; they cannot
    // have changed.
    if fitter.check_minimal_rt_span(rt_bounds, settings.min_rt_span) {
        return rejected(ABORT_MIN_RT_SPAN);
    }
    let mut theoretical: Vec<f64> = Vec::new();
    let mut real: Vec<f64> = Vec::new();
    let mut deviation = 0.0;
    for trace in traces.iter() {
        for k in 0..trace.peaks.len() {
            let theo = traces.baseline + fitter.compute_theoretical(trace, k)?;
            theoretical.push(theo);
            let measured = f64::from(trace.peaks[k].intensity);
            real.push(measured);
            deviation += (measured - theo).abs() / theo;
        }
    }
    let fit_score = max0(1.0 - (deviation / traces.peak_count() as f64));
    let correlation = max0(pearson_correlation_coefficient(&theoretical, &real)?);
    let final_score = (correlation * fit_score).sqrt();
    if log.enabled() {
        put_all(log, &["Quality estimation:\n"]);
        put_all(log, &[" - relative deviation: ", &g(fit_score), "\n"]);
        put_all(log, &[" - correlation: ", &g(correlation), "\n"]);
        put_all(log, &[" => final score: ", &g(final_score), "\n"]);
    }
    if final_score < settings.min_feature_score {
        return rejected(ABORT_QUALITY_TOO_LOW);
    }
    Ok(QualityOutcome::Accepted(FeatureQuality {
        fit_score,
        correlation,
        final_score,
    }))
}

/// Everything source step 3.3.5 reads while it builds one feature.
#[derive(Clone, Copy, Debug)]
pub struct FeatureInput<'a> {
    /// The fitted retention-time model.
    pub model: &'a FittedModel,
    /// The cropped mass traces, which the source has assigned back to `traces`
    /// before this step.
    pub traces: &'a MassTraces,
    /// The best isotope pattern of the seed; only its theoretical pattern's
    /// `trimmed_left` is read.
    pub pattern: &'a IsotopePattern,
    /// The precalculated isotope windows, read once by m/z for the intensity.
    pub windows: &'a IsotopeWindows,
    /// The typed parameters.
    pub settings: &'a Settings,
    /// The charge hypothesis of this candidate.
    pub charge: i32,
    /// The provisional label: the source's `plot_nr`, overwritten with the
    /// feature number when the candidate survives the containment pass.
    pub plot_nr: i64,
    /// The quality [`check_feature_quality`] returned.
    pub quality: FeatureQuality,
}

/// Build the feature of one accepted candidate: source step 3.3.5.
///
/// Fields, in the source's order: the provisional label (meta key `label`, the
/// `MetaInfoRegistry` index 3), the charge, the overall quality (narrowed to
/// `f32` by `setOverallQuality`), the meta values `score_fit` and
/// `score_correlation`, the retention time (the model's centre), the width (the
/// model's FWHM, narrowed to `f32`, which also writes the meta value `FWHM`),
/// the meta value `num_of_datapoints`, the three `EGH_*` meta values for an
/// asymmetric fit, the m/z, the intensity and one convex hull per trace.
///
/// The m/z follows `feature:reported_mz`:
///
/// - `maximum` — the average m/z of the trace with the highest theoretical
///   intensity;
/// - `average` — the intensity-weighted average m/z of every remaining peak;
/// - `monoisotopic` — the `maximum` value minus `PROTON_MASS_U / charge` times
///   the sum of that trace's *index* and the pattern's `trimmed_left`.
///
/// The intensity is the model's area divided by the maximum of the isotope
/// window **of the feature's m/z**, not of its mass (B7 candidate 4 in
/// `docs/FEATURE_FINDER_PICKED_SUPPORT.md`). The monoisotopic spacing likewise
/// uses `PROTON_MASS_U` and a trace index rather than the isotope mass
/// difference (B7 candidate 3).
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the traces are empty
/// ([`MassTraces::theoretical_max_position`]); when the isotope window of the
/// feature's m/z was not precalculated, which a NaN m/z causes (an infinite
/// intensity in the reported traces makes the average m/z `inf / inf`) and
/// where the source's exception escapes its parallel region and terminates the
/// process; when the point count of a hull exceeds
/// its ceiling, and when the model's FWHM is not finite or is negative, which
/// [`crate::kernel::BaseFeature::set_width`] refuses and the source stores. A
/// non-finite FWHM needs a non-finite fit, which the source's quality checks let
/// through because every comparison against a NaN is false.
pub fn build_feature(input: FeatureInput<'_>) -> Result<Feature> {
    let FeatureInput {
        model,
        traces,
        pattern,
        windows,
        settings,
        charge,
        plot_nr,
        quality,
    } = input;
    let fitter = model.as_fitter();
    let mut feature = Feature::new(0.0, 0.0, 0.0);
    feature
        .metadata
        .insert("label".into(), MetaValue::from(plot_nr));
    feature.charge = charge;
    feature.quality = quality.final_score as f32;
    feature
        .metadata
        .insert("score_fit".into(), MetaValue::try_from(quality.fit_score)?);
    feature.metadata.insert(
        "score_correlation".into(),
        MetaValue::try_from(quality.correlation)?,
    );
    feature.rt = fitter.center();
    feature.base.set_width(fitter.fwhm() as f32)?;
    let datapoints: usize = traces.iter().map(|trace| trace.peaks.len()).sum();
    feature
        .metadata
        .insert(NUM_OF_DATAPOINTS.into(), MetaValue::from(datapoints as i64));
    if let Some(egh) = model.egh() {
        feature
            .metadata
            .insert("EGH_tau".into(), MetaValue::try_from(egh.tau())?);
        feature
            .metadata
            .insert("EGH_height".into(), MetaValue::try_from(egh.height())?);
        feature
            .metadata
            .insert("EGH_sigma".into(), MetaValue::try_from(egh.sigma())?);
    }
    let c = f64::from(charge);
    feature.mz = match settings.reported_mz {
        ReportedMz::Maximum => traces[traces.theoretical_max_position()?].avg_mz(),
        ReportedMz::Average => {
            let mut total_intensity = 0.0;
            let mut average_mz = 0.0;
            for trace in traces.iter() {
                for peak in &trace.peaks {
                    average_mz += peak.mz * f64::from(peak.intensity);
                    total_intensity += f64::from(peak.intensity);
                }
            }
            average_mz / total_intensity
        }
        ReportedMz::Monoisotopic => {
            let position = traces.theoretical_max_position()?;
            let mono = traces[position].avg_mz();
            mono - (PROTON_MASS_U / c)
                * (position + pattern.theoretical_pattern.trimmed_left) as f64
        }
    };
    // Source `getIsotopeDistribution_(f.getMZ())` inside the seed loop's
    // OpenMP region (`FeatureFinderAlgorithmPicked.cpp:790`): its
    // `Exception::InvalidValue` is not caught there, so the source terminates.
    let window = windows.get(feature.mz).map_err(|error| {
        Error::InvalidValue(format!(
            "FeatureFinderAlgorithmPicked step 3.3.5: the feature m/z {} has no isotope window \
             ({error}); the source throws this inside its OpenMP region, where std::terminate \
             ends the process",
            feature.mz
        ))
    })?;
    feature.intensity = (fitter.area() / window.max) as f32;
    let mut hulls = Vec::new();
    hulls
        .try_reserve_exact(traces.len())
        .map_err(|_| Error::InvalidValue("cannot allocate the feature's convex hulls".into()))?;
    for trace in traces.iter() {
        hulls.push(trace.convex_hull()?);
    }
    feature.convex_hulls = hulls;
    Ok(feature)
}
