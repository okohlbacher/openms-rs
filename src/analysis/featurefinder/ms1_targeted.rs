// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The steps `FeatureFinderIdentificationAlgorithm` and
//! `FeatureFinderAlgorithmMetaboIdent` have in common.
//!
//! Both classes build an assay library for MS1 traces, extract the chromatograms
//! from the MS1 spectra, hand them to `MRMFeatureFinderScoring`, annotate the
//! picked features and merge FAIMS groups. The source repeats these steps in
//! each class; this module holds one copy of each so that the two ports cannot
//! drift apart. Nothing here is public API: the two algorithm modules re-state
//! what they document.

use crate::analysis::featurefinder::elution_model_fitter::ElutionModelFitter;
use crate::analysis::openswath::chromatogram_extractor::{
    ChromatogramExtractor, prepare_coordinates, return_chromatogram,
};
use crate::analysis::openswath::chromatogram_extractor_algorithm::ExtractionCoordinates;
use crate::analysis::openswath::data_structures::Chromatogram;
use crate::analysis::openswath::spectrum_access::SpectrumAccessPtr;
use crate::analysis::openswath::transition_experiment::LightTargetedExperiment;
use crate::kernel::geometry::{ConvexHull2D, Point2D};
use crate::kernel::{Feature, FeatureMap, MSChromatogram, MSExperiment, MSSpectrum};
use crate::param::{Param, ParamValue};
use crate::processing::feature_overlap_filter::FeatureOverlapFilter;
use crate::{Error, Result};

/// Largest retention time difference (seconds) of two FAIMS features that are
/// merged (`FeatureOverlapFilter::mergeFAIMSFeatures` argument of both classes).
pub(crate) const FAIMS_MERGE_RT_TOLERANCE: f64 = 5.0;

/// Largest m/z difference (Da) of two FAIMS features that are merged.
pub(crate) const FAIMS_MERGE_MZ_TOLERANCE: f64 = 0.05;

/// Keep only the MS1 spectra (source `setMSData`, both overloads).
pub(crate) fn keep_ms1(experiment: &mut MSExperiment) {
    experiment.spectra.retain(|spectrum| spectrum.ms_level == 1);
}

/// A feature's RT range as the `leftWidth` and `rightWidth` meta values the
/// OpenSWATH picker writes.
pub(crate) fn feature_bounds(feature: &Feature) -> Result<(f64, f64)> {
    let read = |key: &str| -> Result<f64> {
        feature
            .metadata
            .get(key)
            .ok_or_else(|| {
                Error::MissingInformation(format!("a picked feature has no '{key}' meta value"))
            })?
            .as_f64()
    };
    Ok((read("leftWidth")?, read("rightWidth")?))
}

/// Give a feature a hull per mass trace when it has none (source
/// `ensureConvexHulls_`, identical in both classes): a rectangle spanning the
/// feature's RT bounds and the m/z extraction window around the trace.
///
/// # Errors
///
/// [`Error::MissingInformation`] without the `leftWidth` or `rightWidth` meta
/// values, and [`Error::InvalidValue`] for a non-finite hull corner.
pub(crate) fn ensure_convex_hulls(
    feature: &mut Feature,
    mz_window: f64,
    mz_window_ppm: bool,
) -> Result<()> {
    if !feature.convex_hulls.is_empty() {
        return Ok(());
    }
    // The source reads the bounds even when the feature has no subordinate.
    let (rt_min, rt_max) = feature_bounds(feature)?;
    let mut hulls = Vec::with_capacity(feature.subordinates.len());
    for sub in &feature.subordinates {
        let mut abs_mz_tol = mz_window / 2.0;
        if mz_window_ppm {
            abs_mz_tol = sub.mz * abs_mz_tol * 1.0e-6;
        }
        let mut hull = ConvexHull2D::new();
        hull.add_point(Point2D::new(rt_min, sub.mz - abs_mz_tol))?;
        hull.add_point(Point2D::new(rt_min, sub.mz + abs_mz_tol))?;
        hull.add_point(Point2D::new(rt_max, sub.mz - abs_mz_tol))?;
        hull.add_point(Point2D::new(rt_max, sub.mz + abs_mz_tol))?;
        hulls.push(hull);
    }
    feature.convex_hulls = hulls;
    Ok(())
}

/// The finder parameters both classes set before picking (the lines that are
/// identical in the two `run` bodies).
///
/// The MS1 level makes the fragment-ion scores meaningless, and the RT of an
/// assay is not reliable, so those scores are switched off. `EMGScoring`
/// settings are propagated from the caller's parameters.
pub(crate) fn apply_common_finder_parameters(
    finder: &mut Param,
    emg_max_iteration: &ParamValue,
    emg_init_mom: &ParamValue,
    gauss_width: f64,
    min_peak_width: f64,
    signal_to_noise: f64,
) -> Result<()> {
    let mut set = |key: &str, value: ParamValue| finder.set_value(key, value, "", &[]);
    // return all features
    set("stop_report_after_feature", ParamValue::Integer(-1))?;
    // propagate the settings to the sub algorithms
    set("EMGScoring:max_iteration", emg_max_iteration.clone())?;
    set("EMGScoring:init_mom", emg_init_mom.clone())?;
    // RT may not be reliable
    set("Scores:use_rt_score", ParamValue::String("false".into()))?;
    // only MS1 spectra are used, so these scores are useless
    set(
        "Scores:use_ionseries_scores",
        ParamValue::String("false".into()),
    )?;
    set(
        "Scores:use_ms2_isotope_scores",
        ParamValue::String("false".into()),
    )?;
    // redundant to the "MS2" correlation and, because the precursor transition
    // is the first product transition, additionally biased
    set(
        "Scores:use_ms1_correlation",
        ParamValue::String("false".into()),
    )?;
    // the same, on the MS1 level only the "MS1 fullscan" scores matter
    set("Scores:use_ms1_mi", ParamValue::String("false".into()))?;
    set(
        "TransitionGroupPicker:PeakPickerChromatogram:gauss_width",
        ParamValue::Float(gauss_width),
    )?;
    set(
        "TransitionGroupPicker:min_peak_width",
        ParamValue::Float(min_peak_width),
    )?;
    // a zero signal-to-noise threshold breaks the OpenSWATH feature detection
    set(
        "TransitionGroupPicker:PeakPickerChromatogram:signal_to_noise",
        ParamValue::Float(signal_to_noise),
    )?;
    set(
        "TransitionGroupPicker:PeakPickerChromatogram:write_sn_log_messages",
        ParamValue::String("false".into()),
    )?;
    set(
        "TransitionGroupPicker:recalculate_peaks",
        ParamValue::String("true".into()),
    )?;
    set(
        "TransitionGroupPicker:PeakPickerChromatogram:peak_width",
        ParamValue::Float(-1.0),
    )?;
    set(
        "TransitionGroupPicker:PeakPickerChromatogram:method",
        ParamValue::String("corrected".into()),
    )?;
    Ok(())
}

/// Extract the MS1 traces of every coordinate of `light_library` (the common
/// extraction block of both `run` bodies): prepare the coordinates with an
/// unlimited RT window, extract with a top-hat filter, with the ion mobility
/// window when `im_window` is given, and convert the traces to annotated
/// chromatograms using the first spectrum for the settings.
///
/// The chromatograms are appended to `output`; nothing is appended on error.
///
/// # Errors
///
/// The errors of `prepare_coordinates`, the extraction and `return_chromatogram`.
pub(crate) fn extract_ms1_chromatograms(
    spectra: &SpectrumAccessPtr,
    settings: &MSSpectrum,
    light_library: &LightTargetedExperiment,
    mz_window: f64,
    mz_window_ppm: bool,
    im_window: Option<f64>,
    output: &mut Vec<MSChromatogram>,
) -> Result<()> {
    let extractor = ChromatogramExtractor::new();
    let mut chrom_temp: Vec<Chromatogram> = Vec::new();
    let mut coords: Vec<ExtractionCoordinates> = Vec::new();
    // The NaN window relies on rt_start/rt_end carried on the light compound,
    // populated by convertTargetedExp from the RTs of the library.
    prepare_coordinates(
        &mut chrom_temp,
        &mut coords,
        light_library,
        f64::NAN,
        false,
        0,
    )?;
    match im_window {
        Some(window) => extractor.extract_chromatograms_im(
            spectra.as_ref(),
            &mut chrom_temp,
            &coords,
            mz_window,
            mz_window_ppm,
            window,
            "tophat",
        )?,
        None => extractor.extract_chromatograms(
            spectra.as_ref(),
            &mut chrom_temp,
            &coords,
            mz_window,
            mz_window_ppm,
            "tophat",
        )?,
    }
    return_chromatogram(
        &chrom_temp,
        &coords,
        light_library,
        settings,
        output,
        false,
        0.0,
    )
}

/// The elution model fitter configured as both classes configure it: the
/// `model:` section of the algorithm parameters without `type`, with
/// `asymmetric` following the model type.
///
/// # Errors
///
/// The errors of the parameter tree and of the fitter's parameter validation.
pub(crate) fn elution_model_fitter(param: &Param, model_type: &str) -> Result<ElutionModelFitter> {
    let mut emf_params = param.copy("model:", true)?;
    emf_params.remove("type")?;
    emf_params.set_value(
        "asymmetric",
        ParamValue::String(
            if model_type == "asymmetric" {
                "true"
            } else {
                "false"
            }
            .into(),
        ),
        "",
        &[],
    )?;
    ElutionModelFitter::from_param(&emf_params)
}

/// Merge the features of one analyte seen at several FAIMS compensation voltages
/// when requested (the block both `run` bodies share).
///
/// Returns how many features were merged away.
///
/// # Errors
///
/// The errors of `FeatureOverlapFilter::merge_faims_features`; the map is
/// unchanged then.
pub(crate) fn merge_faims_features(features: &mut FeatureMap, merge: bool) -> Result<usize> {
    if !merge {
        return Ok(0);
    }
    let before = features.features.len();
    FeatureOverlapFilter::merge_faims_features(
        features,
        FAIMS_MERGE_RT_TOLERANCE,
        FAIMS_MERGE_MZ_TOLERANCE,
    )?;
    Ok(before.saturating_sub(features.features.len()))
}
