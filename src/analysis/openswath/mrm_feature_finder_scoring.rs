// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Finding and scoring of MRM features among co-eluting chromatograms.
//!
//! Ports `ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h` and its `.cpp`.
//! [`MRMFeatureFinderScoring`] maps an input of chromatograms and a transition
//! list onto transition groups, picks the peak groups of every group with
//! [`MRMTransitionGroupPicker`] and scores each with [`OpenSwathScoring`],
//! adding the scores and the identification of the assay to the features of a
//! [`FeatureMap`].
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/MRM_FEATURE_FINDER_SUPPORT.md).
//!
//! # Example
//!
//! ```
//! use openms::analysis::openswath::mrm_feature_finder_scoring::MRMFeatureFinderScoring;
//!
//! let finder = MRMFeatureFinderScoring::new().unwrap();
//! assert_eq!(finder.rt_normalization_factor, 1.0);
//! assert_eq!(finder.scoring_model, "default");
//! ```
//!
//! # Parallelism
//!
//! The source scores the peak groups of one transition group in an OpenMP
//! parallel loop and collects them under a critical section, so the order in
//! which the scored features reach its list depends on the thread schedule
//! before the final sort by quality. This port scores them in order on one
//! thread; its result is one the source can produce with a single thread.
//!
//! # Not ported
//!
//! The `MobilogramParquetConsumer` of the source is a
//! [`MobilogramSink`]; the Parquet implementation is the backlog item of
//! `FORMAT/DATAACCESS/MobilogramParquetConsumer.h`.
//!
//! [`MRMFeatureFinderScoring`]: crate::analysis::openswath::mrm_feature_finder_scoring::MRMFeatureFinderScoring
//! [`MRMTransitionGroupPicker`]: crate::analysis::openswath::mrm_transition_group_picker::MRMTransitionGroupPicker
//! [`OpenSwathScoring`]: crate::analysis::openswath::openswath_scoring::OpenSwathScoring
//! [`FeatureMap`]: crate::kernel::FeatureMap
//! [`MobilogramSink`]: crate::analysis::openswath::ion_mobility_scoring::MobilogramSink

use crate::analysis::feature_finder_picked::glibc_libm;
use crate::analysis::featurefinder::emg_scoring::EmgScoring;
use crate::analysis::openswath::data_access_helper::{
    convert_targeted_exp, convert_to_openms_chromatogram, convert_to_openms_chromatogram_filter,
};
use crate::analysis::openswath::dia_scoring::DiaScoring;
use crate::analysis::openswath::ion_mobility_scoring::MobilogramSink;
use crate::analysis::openswath::itransition::SignalToNoisePtr;
use crate::analysis::openswath::mrm_feature_access::{MRMFeatureOpenMS, SignalToNoiseOpenMS};
use crate::analysis::openswath::mrm_feature_support::meta_double;
use crate::analysis::openswath::mrm_scoring::MRMScoring;
use crate::analysis::openswath::mrm_transition_group_picker::MRMTransitionGroupPicker;
use crate::analysis::openswath::openswath_scores::OpenSwathScoresUsage;
use crate::analysis::openswath::openswath_scoring::OpenSwathScoring;
use crate::analysis::openswath::spectrum_access::SpectrumAccessPtr;
use crate::analysis::openswath::spectrum_access_openms::SimpleOpenMSSpectraFactory;
use crate::analysis::openswath::swath_map::SwathMap;
use crate::analysis::openswath::targeted_experiment::TargetedExperiment;
use crate::analysis::openswath::transition_experiment::{
    LightCompound, LightTargetedExperiment, LightTransition,
};
use crate::analysis::transformations::TransformationDescription;
use crate::chemistry::AASequence;
use crate::chemistry::digestion::{ProteaseDB, ProteaseDigestion};
use crate::concept::unique_id::UniqueIdGenerator;
use crate::identification::{
    PeptideEvidence, PeptideHit, PeptideIdentification, ProteinHit, ProteinIdentification,
};
use crate::kernel::mrm::{MRMFeature, MRMTransitionGroup, OpenSwathIndScores, OpenSwathScores};
use crate::kernel::ranges::RangeBase;
use crate::kernel::{Feature, FeatureMap, MSChromatogram, MSExperiment, Precursor};
use crate::math::source_sort::{source_sort_by, source_sort_reversed_by};
use crate::metadata::MetaValue;
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

/// The source `DefaultParamHandler` name.
pub const MRM_FEATURE_FINDER_SCORING_NAME: &str = "MRMFeatureFinderScoring";

/// The identifier of the single protein identification run the finder
/// writes (source macro `run_identifier`).
pub const RUN_IDENTIFIER: &str = "unique_run_identifier";

/// A transition group holds the chromatographic data and peaks across the
/// chromatograms of one compound (source `MRMTransitionGroupType`).
pub type MRMTransitionGroupType = MRMTransitionGroup<MSChromatogram, LightTransition>;

/// Transition groups by transition group ID (source `TransitionGroupMapType`).
pub type TransitionGroupMapType = BTreeMap<String, MRMTransitionGroupType>;

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

fn log(value: f64) -> f64 {
    glibc_libm::log(value)
}

fn entry(name: &str, value: ParamValue, description: &str, advanced: bool) -> ParamEntry {
    let mut e = ParamEntry {
        name: name.into(),
        description: description.into(),
        value,
        ..ParamEntry::default()
    };
    if advanced {
        e.tags.insert("advanced".into());
    }
    e
}

fn flag(name: &str, value: bool, description: &str, advanced: bool) -> ParamEntry {
    let mut e = entry(
        name,
        ParamValue::String(if value { "true" } else { "false" }.into()),
        description,
        advanced,
    );
    e.valid_strings = vec!["true".into(), "false".into()];
    e
}

fn choice(
    name: &str,
    value: &str,
    description: &str,
    advanced: bool,
    valid: &[&str],
) -> ParamEntry {
    let mut e = entry(
        name,
        ParamValue::String(value.into()),
        description,
        advanced,
    );
    e.valid_strings = valid.iter().map(|s| (*s).to_owned()).collect();
    e
}

fn float(
    name: &str,
    value: f64,
    description: &str,
    advanced: bool,
    min: Option<f64>,
) -> ParamEntry {
    let mut e = entry(name, ParamValue::Float(value), description, advanced);
    if let Some(min) = min {
        e.min_float = min;
    }
    e
}

fn int(name: &str, value: i64, description: &str, advanced: bool, min: Option<i32>) -> ParamEntry {
    let mut e = entry(name, ParamValue::Integer(value), description, advanced);
    if let Some(min) = min {
        e.min_int = min;
    }
    e
}

/// A shorter reborrow of an optional sink, so that one sink can serve many
/// calls.
fn reborrow<'s>(
    sink: &'s mut Option<&mut dyn MobilogramSink>,
) -> Option<&'s mut dyn MobilogramSink> {
    match sink {
        Some(s) => Some(&mut **s),
        None => None,
    }
}

/// `SqliteHelper::clearSignBit`: the unique ID with its top bit cleared, as a
/// non-negative `Int64`.
fn clear_sign_bit(id: u64) -> i64 {
    (id & !(1_u64 << 63)) as i64
}

/// Add a score to a feature (source `MRMFeature::addScore`, which is
/// `setMetaValue`): the value is stored as it is, non-finite ones included, as
/// the source's `DataValue` stores them.
fn add_score(feature: &mut MRMFeature, name: &str, value: f64) {
    feature
        .feature
        .base
        .metadata
        .insert(name.to_owned(), MetaValue::source_float(value));
}

fn set_meta(feature: &mut Feature, name: &str, value: MetaValue) {
    feature.base.metadata.insert(name.to_owned(), value);
}

fn float_meta(feature: &Feature, name: &str) -> Result<f64> {
    meta_double(&feature.base.metadata, name)
}

/// Group-invariant values cached once per transition group (source private
/// `TransitionGroupCache`).
struct TransitionGroupCache {
    normalized_library_intensity: Vec<f64>,
    transition_native_ids: Vec<String>,
    precursor_ids: Vec<String>,
}

/// Finds and scores peak groups (source `MRMFeatureFinderScoring`).
///
/// The parameters live in a [`Param`] tree with the sections
/// `TransitionGroupPicker:`, `DIAScoring:`, `EMGScoring:` and `Scores:`; the
/// public fields are the source's `updateMembers_` copies of them. Set them
/// through [`set_parameters`](Self::set_parameters) so that tree and members
/// stay in step.
#[derive(Clone)]
pub struct MRMFeatureFinderScoring {
    param: Param,
    /// Only extract RT around the expected value: `+/-` this many seconds
    /// (`rt_extraction_window`; `-1` extracts the whole range).
    pub rt_extraction_window: f64,
    /// Cutoff in m/z below which peaks are not used for quantification
    /// (`quantification_cutoff`).
    pub quantification_cutoff: f64,
    /// Stop reporting after this many features, ordered by quality (`-1`
    /// reports all) (`stop_report_after_feature`).
    pub stop_report_after_feature: i32,
    /// Write the points of every feature (`write_convex_hull`).
    pub write_convex_hull: bool,
    /// Error on, rather than skip, a transition without a chromatogram
    /// (`strict`).
    pub strict: bool,
    /// Extract ion mobility in MS1 (`use_ms1_ion_mobility`).
    pub use_ms1_ion_mobility: bool,
    /// Peak-pick the extracted ion mobilograms (`apply_im_peak_picking`).
    pub apply_im_peak_picking: bool,
    /// `default` or `single_transition` (`scoring_model`).
    pub scoring_model: String,
    /// The enzyme used for counting missed cleavages (`enzyme`).
    pub enzyme: String,
    /// Range of the normalised RT (`rt_normalization_factor`).
    pub rt_normalization_factor: f64,
    /// Spectra added up around the apex (`add_up_spectra`).
    pub add_up_spectra: i32,
    /// `simple` or `resample` (`spectrum_addition_method`).
    pub spectrum_addition_method: String,
    /// `fixed` or `dynamic` (`spectrum_merge_method_type`).
    pub spectrum_merge_method_type: String,
    /// Spacing for resampling added spectra (`spacing_for_spectra_resampling`).
    pub spacing_for_spectra_resampling: f64,
    /// Fraction of the peak width to merge (`merge_spectra_by_peak_width_fraction`).
    pub merge_spectra_by_peak_width_fraction: f64,
    /// S/N threshold of identification transitions (`uis_threshold_sn`).
    pub uis_threshold_sn: f64,
    /// Peak area threshold of identification transitions
    /// (`uis_threshold_peak_area`).
    pub uis_threshold_peak_area: f64,
    /// Window length of the S/N estimate, from
    /// `TransitionGroupPicker:PeakPickerChromatogram:sn_win_len`.
    pub sn_win_len: f64,
    /// Bin count of the S/N estimate, from `...:sn_bin_count`.
    pub sn_bin_count: u32,
    /// Whether the S/N estimator logs, from `...:write_sn_log_messages`.
    pub write_log_messages: bool,
    /// Extra drift time extracted for ion mobility scoring (`im_extra_drift`).
    pub im_extra_drift: f64,
    peptide_ref_map: BTreeMap<String, LightCompound>,
    su: OpenSwathScoresUsage,
    diascoring: DiaScoring,
    emgscoring: EmgScoring,
    ms1_map: Option<SpectrumAccessPtr>,
    ids: Arc<Mutex<UniqueIdGenerator>>,
}

impl std::fmt::Debug for MRMFeatureFinderScoring {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MRMFeatureFinderScoring")
            .field("scoring_model", &self.scoring_model)
            .field("enzyme", &self.enzyme)
            .field("rt_normalization_factor", &self.rt_normalization_factor)
            .finish_non_exhaustive()
    }
}

const SCORES_FLAGS: [(&str, bool, &str); 20] = [
    (
        "use_shape_score",
        true,
        "Use the shape score (this score measures the similarity in shape of the transitions using a cross-correlation)",
    ),
    (
        "use_coelution_score",
        true,
        "Use the coelution score (this score measures the similarity in coelution of the transitions using a cross-correlation)",
    ),
    (
        "use_rt_score",
        true,
        "Use the retention time score (this score measure the difference in retention time)",
    ),
    ("use_library_score", true, "Use the library score"),
    (
        "use_elution_model_score",
        true,
        "Use the elution model (EMG) score (this score fits a gaussian model to the peak and checks the fit)",
    ),
    ("use_intensity_score", true, "Use the intensity score"),
    ("use_nr_peaks_score", true, "Use the number of peaks score"),
    ("use_total_xic_score", true, "Use the total XIC score"),
    ("use_total_mi_score", false, "Use the total MI score"),
    ("use_sn_score", true, "Use the SN (signal to noise) score"),
    (
        "use_mi_score",
        false,
        "Use the MI (mutual information) score",
    ),
    (
        "use_dia_scores",
        true,
        "Use the DIA (SWATH) scores. If turned off, will not use fragment ion spectra for scoring.",
    ),
    (
        "use_ms1_correlation",
        false,
        "Use the correlation scores with the MS1 elution profiles",
    ),
    (
        "use_ion_mobility_scores",
        false,
        "Use the scores for Ion Mobility scans",
    ),
    (
        "use_ms1_fullscan",
        false,
        "Use the full MS1 scan at the peak apex for scoring (ppm accuracy of precursor and isotopic pattern)",
    ),
    ("use_ms1_mi", false, "Use the MS1 MI score"),
    (
        "use_uis_scores",
        false,
        "Use UIS scores for peptidoform identification",
    ),
    (
        "use_peak_shape_metrics",
        false,
        "Use peak shape metrics for scoring",
    ),
    (
        "use_ionseries_scores",
        true,
        "Use MS2-level b/y ion-series scores for peptidoform identification",
    ),
    (
        "use_ms2_isotope_scores",
        true,
        "Use MS2-level isotope scores (pearson & manhattan) across product transitions (based on ID if annotated or averagine)",
    ),
];

impl MRMFeatureFinderScoring {
    /// A finder with the source defaults.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure, which the constant defaults
    /// cannot trigger.
    pub fn new() -> Result<Self> {
        Self::from_param(&Param::new())
    }

    /// The source defaults as a [`Param`] (source `getDefaults`): the scalar
    /// parameters and the sections `TransitionGroupPicker:`
    /// ([`MRMTransitionGroupPicker::defaults`]), `DIAScoring:`
    /// ([`DiaScoring::defaults`]), `EMGScoring:` ([`EmgScoring::defaults`]) and
    /// `Scores:` (one switch per score).
    ///
    /// # Errors
    ///
    /// A [`Param`] construction failure from the nested defaults.
    pub fn defaults() -> Result<Param> {
        let enzyme_names: Vec<String> = ProteaseDB::global()
            .names()
            .into_iter()
            .map(str::to_owned)
            .collect();
        let mut enzyme = entry(
            "enzyme",
            ParamValue::String("Trypsin".into()),
            "Enzyme used for counting missed cleavages in peptide sequences",
            false,
        );
        enzyme.valid_strings = enzyme_names;
        let scores_entries: Vec<ParamEntry> = SCORES_FLAGS
            .iter()
            .map(|(name, value, description)| {
                let mut e = flag(name, *value, description, true);
                if *name == "use_peak_shape_metrics" {
                    // the source gives this one no valid strings
                    e.valid_strings.clear();
                }
                e
            })
            .collect();
        let mut top = Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries: vec![
                int(
                    "stop_report_after_feature",
                    -1,
                    "Stop reporting after feature (ordered by quality; -1 means do not stop).",
                    false,
                    None,
                ),
                float(
                    "rt_extraction_window",
                    -1.0,
                    "Only extract RT around this value (-1 means extract over the whole range, a value of 500 means to extract around +/- 500 s of the expected elution). For this to work, the TraML input file needs to contain normalized RT values.",
                    false,
                    None,
                ),
                float(
                    "rt_normalization_factor",
                    1.0,
                    "The normalized RT is expected to be between 0 and 1. If your normalized RT has a different range, pass this here (e.g. it goes from 0 to 100, set this value to 100)",
                    false,
                    None,
                ),
                float(
                    "quantification_cutoff",
                    0.0,
                    "Cutoff in m/z below which peaks should not be used for quantification any more",
                    true,
                    Some(0.0),
                ),
                flag(
                    "write_convex_hull",
                    false,
                    "Whether to write out all points of all features into the featureXML",
                    true,
                ),
                choice(
                    "spectrum_addition_method",
                    "simple",
                    "For spectrum addition, either use simple concatenation or use peak resampling",
                    true,
                    &["simple", "resample"],
                ),
                choice(
                    "spectrum_merge_method_type",
                    "fixed",
                    "For spectrum addition, either use a fixed number of spectra or dynamically select the number of spectra to add around the peak apex based on the merge_spectra_by_peak_width_fraction.",
                    true,
                    &["fixed", "dynamic"],
                ),
                int(
                    "add_up_spectra",
                    1,
                    "Add up spectra on the left and right around the retention time peak apex.",
                    true,
                    Some(1),
                ),
                float(
                    "spacing_for_spectra_resampling",
                    0.005,
                    "If spectra are to be added, use this spacing to add them up",
                    true,
                    Some(0.0),
                ),
                {
                    let mut e = float(
                        "merge_spectra_by_peak_width_fraction",
                        0.15,
                        "If spectra are to be added based on the peak width of peak, construct number of spectra to be added based on N percent of number of points of peak width.",
                        true,
                        Some(0.0001),
                    );
                    e.max_float = 1.0;
                    e
                },
                int(
                    "uis_threshold_sn",
                    -1,
                    "S/N threshold to consider identification transition (set to -1 to consider all)",
                    false,
                    None,
                ),
                int(
                    "uis_threshold_peak_area",
                    0,
                    "Peak area threshold to consider identification transition (set to -1 to consider all)",
                    false,
                    None,
                ),
                choice(
                    "scoring_model",
                    "default",
                    "Scoring model to use",
                    true,
                    &["default", "single_transition"],
                ),
                float(
                    "im_extra_drift",
                    0.0,
                    "Extra drift time to extract for IM scoring (as a fraction, e.g. 0.25 means 25% extra on each side)",
                    true,
                    Some(0.0),
                ),
                flag(
                    "strict",
                    true,
                    "Whether to error (true) or skip (false) if a transition in a transition group does not have a corresponding chromatogram.",
                    true,
                ),
                flag(
                    "use_ms1_ion_mobility",
                    true,
                    "Performs ion mobility extraction in MS1. Set to false if MS1 spectra do not contain ion mobility",
                    true,
                ),
                flag(
                    "apply_im_peak_picking",
                    false,
                    "Perform peak picking on the extracted ion mobilograms. This is useful for reducing intefering signals from co-eluting analytes in the ion mobility dimension. The peak picking will take the highest peak and discard the remaining peaks for ion mobility scoring. ",
                    true,
                ),
                enzyme,
            ],
            nodes: Vec::new(),
        })?;
        top.insert(
            "TransitionGroupPicker:",
            &MRMTransitionGroupPicker::defaults()?,
        )?;
        top.insert("DIAScoring:", &DiaScoring::defaults()?)?;
        top.insert("EMGScoring:", &EmgScoring::new().defaults()?)?;
        let scores = Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries: scores_entries,
            nodes: Vec::new(),
        })?;
        top.insert("Scores:", &scores)?;
        let tail = Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries: vec![
                float(
                    "mz_extraction_window",
                    -1.0,
                    "m/z extraction window (full width) to be used for matching/extraction. If -1, defaults in workflow will be used.",
                    true,
                    None,
                ),
                choice(
                    "mz_extraction_window_unit",
                    "Th",
                    "Unit for mz_extraction_window: 'Th' or 'ppm'",
                    true,
                    &["Th", "ppm"],
                ),
                float(
                    "mz_extraction_window_ms1",
                    -1.0,
                    "m/z extraction window for MS1 (full width). If -1, defaults in workflow will be used.",
                    true,
                    None,
                ),
                choice(
                    "mz_extraction_window_ms1_unit",
                    "Th",
                    "Unit for mz_extraction_window_ms1: 'Th' or 'ppm'",
                    true,
                    &["Th", "ppm"],
                ),
                float(
                    "irt_mz_extraction_window",
                    -1.0,
                    "m/z extraction window (full width) to be used specifically for iRT matching/extraction. If -1, defaults in workflow will be used.",
                    true,
                    None,
                ),
                choice(
                    "irt_mz_extraction_window_unit",
                    "Th",
                    "Unit for irt_mz_extraction_window: 'Th' or 'ppm'",
                    true,
                    &["Th", "ppm"],
                ),
            ],
            nodes: Vec::new(),
        })?;
        top.merge(&tail)?;
        Ok(top)
    }

    /// The current parameters (source `getParameters`).
    pub fn parameters(&self) -> &Param {
        &self.param
    }

    /// Build a finder from parameters, as source `setParameters` followed by
    /// `updateMembers_`; missing ones take their defaults.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a parameter has the wrong type or violates
    /// its restriction, and the errors of building the nested sections.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`from_param`](Self::from_param), also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`from_param`](Self::from_param).
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(MRM_FEATURE_FINDER_SCORING_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    /// Replace the parameters (source `setParameters`); the members follow.
    ///
    /// # Errors
    ///
    /// As [`from_param`](Self::from_param); the finder is unchanged on error.
    /// Its transition lookup and MS1 map are kept.
    pub fn set_parameters(&mut self, parameters: &Param) -> Result<Vec<String>> {
        let (mut fresh, warnings) = Self::from_param_with_warnings(parameters)?;
        fresh.peptide_ref_map = std::mem::take(&mut self.peptide_ref_map);
        fresh.ms1_map = self.ms1_map.take();
        fresh.ids = Arc::clone(&self.ids);
        *self = fresh;
        Ok(warnings)
    }

    fn from_complete_param(param: &Param) -> Result<Self> {
        let text = |key: &str| -> Result<String> { Ok(param.value(key)?.as_str()?.to_owned()) };
        let flag = |key: &str| -> Result<bool> { param.value(key)?.to_bool() };
        let int = |key: &str| -> Result<i32> {
            i32::try_from(param.value(key)?.to_i64()?)
                .map_err(|_| bad(format!("parameter '{key}' does not fit an int")))
        };
        let double = |key: &str| -> Result<f64> { param.value(key)?.to_f64() };
        let su = OpenSwathScoresUsage {
            use_coelution_score: flag("Scores:use_coelution_score")?,
            use_shape_score: flag("Scores:use_shape_score")?,
            use_rt_score: flag("Scores:use_rt_score")?,
            use_library_score: flag("Scores:use_library_score")?,
            use_elution_model_score: flag("Scores:use_elution_model_score")?,
            use_intensity_score: flag("Scores:use_intensity_score")?,
            use_total_xic_score: flag("Scores:use_total_xic_score")?,
            use_total_mi_score: flag("Scores:use_total_mi_score")?,
            use_nr_peaks_score: flag("Scores:use_nr_peaks_score")?,
            use_sn_score: flag("Scores:use_sn_score")?,
            use_mi_score: flag("Scores:use_mi_score")?,
            use_dia_scores: flag("Scores:use_dia_scores")?,
            use_im_scores: flag("Scores:use_ion_mobility_scores")?,
            use_ms1_correlation: flag("Scores:use_ms1_correlation")?,
            use_ms1_fullscan: flag("Scores:use_ms1_fullscan")?,
            use_ms1_mi: flag("Scores:use_ms1_mi")?,
            use_uis_scores: flag("Scores:use_uis_scores")?,
            use_ionseries_scores: flag("Scores:use_ionseries_scores")?,
            use_ms2_isotope_scores: flag("Scores:use_ms2_isotope_scores")?,
            use_peak_shape_metrics: flag("Scores:use_peak_shape_metrics")?,
        };
        let sn_bin_count = param
            .value("TransitionGroupPicker:PeakPickerChromatogram:sn_bin_count")?
            .to_i64()?;
        let mut emgscoring = EmgScoring::new();
        emgscoring.set_fitter_param(param.copy("EMGScoring:", true)?);
        Ok(Self {
            param: param.clone(),
            rt_extraction_window: double("rt_extraction_window")?,
            quantification_cutoff: double("quantification_cutoff")?,
            stop_report_after_feature: int("stop_report_after_feature")?,
            write_convex_hull: flag("write_convex_hull")?,
            strict: flag("strict")?,
            use_ms1_ion_mobility: flag("use_ms1_ion_mobility")?,
            apply_im_peak_picking: flag("apply_im_peak_picking")?,
            scoring_model: text("scoring_model")?,
            enzyme: text("enzyme")?,
            rt_normalization_factor: double("rt_normalization_factor")?,
            add_up_spectra: int("add_up_spectra")?,
            spectrum_addition_method: text("spectrum_addition_method")?,
            spectrum_merge_method_type: text("spectrum_merge_method_type")?,
            spacing_for_spectra_resampling: double("spacing_for_spectra_resampling")?,
            merge_spectra_by_peak_width_fraction: double("merge_spectra_by_peak_width_fraction")?,
            uis_threshold_sn: double("uis_threshold_sn")?,
            uis_threshold_peak_area: double("uis_threshold_peak_area")?,
            sn_win_len: double("TransitionGroupPicker:PeakPickerChromatogram:sn_win_len")?,
            // `(unsigned int)param_.getValue(...)`: a negative value wraps.
            sn_bin_count: sn_bin_count as u32,
            write_log_messages: flag(
                "TransitionGroupPicker:PeakPickerChromatogram:write_sn_log_messages",
            )?,
            im_extra_drift: double("im_extra_drift")?,
            peptide_ref_map: BTreeMap::new(),
            su,
            diascoring: DiaScoring::from_param(&param.copy("DIAScoring:", true)?)?,
            emgscoring,
            ms1_map: None,
            ids: Arc::new(Mutex::new(UniqueIdGenerator::new())),
        })
    }

    /// Set the strict flag (source `setStrictFlag`): error on, or skip, a
    /// transition without a chromatogram.
    pub fn set_strict_flag(&mut self, strict: bool) {
        self.strict = strict;
    }

    /// Provide an MS1 map for DIA precursor scores (source `setMS1Map`).
    pub fn set_ms1_map(&mut self, ms1_map: Option<SpectrumAccessPtr>) {
        self.ms1_map = ms1_map;
    }

    /// Replace the generator of unique feature IDs (the source draws them
    /// from the global `UniqueIdGenerator`).
    pub fn set_unique_id_generator(&mut self, generator: UniqueIdGenerator) {
        self.ids = Arc::new(Mutex::new(generator));
    }

    /// Pick and score the features of one experiment given as OpenMS
    /// structures (source `pickExperiment(const PeakMap&, FeatureMap&,
    /// const TargetedExperiment&, ...)`, the overload "for wrapping in
    /// Python"): the chromatograms, the transition list, an optional
    /// transformation of the retention time to the normalised space and an
    /// optional SWATH (DIA) map. The transition group map is not returned.
    ///
    /// # Errors
    ///
    /// The errors of [`convert_targeted_exp`] and of
    /// [`pick_experiment_light`](Self::pick_experiment_light).
    pub fn pick_experiment(
        &mut self,
        chromatograms: &MSExperiment,
        output: &mut FeatureMap,
        transition_exp: &TargetedExperiment,
        trafo: &TransformationDescription,
        swath_map: &MSExperiment,
    ) -> Result<()> {
        let transitions = convert_targeted_exp(transition_exp)?;
        let chromatogram_ptr = SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(Arc::new(
            chromatograms.clone(),
        ))?;
        let swath_ptr =
            SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(Arc::new(swath_map.clone()))?;
        let swath = SwathMap {
            sptr: Some(swath_ptr),
            ..SwathMap::default()
        };
        let mut transition_group_map = TransitionGroupMapType::new();
        self.pick_experiment_light(
            &chromatogram_ptr,
            output,
            &transitions,
            trafo,
            &[swath],
            &mut transition_group_map,
        )
    }

    /// Pick and score the features of one experiment (source
    /// `pickExperiment(const SpectrumAccessPtr&, FeatureMap&,
    /// const LightTargetedExperiment&, ...)`).
    ///
    /// Stores the proteins of the transition list as one protein
    /// identification run in `output`, maps the chromatograms onto the
    /// transition groups ([`map_experiment_to_transition_list`](Self::map_experiment_to_transition_list)),
    /// then for every group with chromatograms and transitions picks the peak
    /// groups with the group picker and scores them
    /// ([`score_peakgroups`](Self::score_peakgroups)). `transition_group_map`
    /// receives the groups with their peak groups. With `Scores:use_total_mi_score`
    /// the picker is told to compute the total mutual information.
    ///
    /// # Errors
    ///
    /// The errors of mapping, picking and scoring.
    pub fn pick_experiment_light(
        &mut self,
        input: &SpectrumAccessPtr,
        output: &mut FeatureMap,
        transition_exp: &LightTargetedExperiment,
        trafo: &TransformationDescription,
        swath_maps: &[SwathMap],
        transition_group_map: &mut TransitionGroupMapType,
    ) -> Result<()> {
        // Step 1: store the peptide retention times in an intermediate map
        self.prepare_protein_peptide_maps(transition_exp);
        // store the proteins of the input in the output feature map
        let prot_hits: Vec<ProteinHit> = transition_exp
            .proteins()
            .iter()
            .map(|prot| ProteinHit {
                sequence: prot.sequence.clone(),
                accession: prot.id.clone(),
                ..ProteinHit::default()
            })
            .collect();
        let prot_id = ProteinIdentification {
            hits: prot_hits,
            identifier: RUN_IDENTIFIER.to_owned(),
            ..ProteinIdentification::default()
        };
        output.protein_identifications.push(prot_id);

        // Step 2: create all transition groups from the individual transitions
        self.map_experiment_to_transition_list(
            input,
            transition_exp,
            transition_group_map,
            trafo.clone(),
            self.rt_extraction_window,
        )?;

        // Step 3: create the consensus features, then score them
        let mut picker_param = self.param.copy("TransitionGroupPicker:", true)?;
        // If use_total_mi_score is defined, the picker has to compute it
        if self.su.use_total_mi_score {
            picker_param.set_value(
                "compute_total_mi",
                ParamValue::String("true".into()),
                "",
                &[],
            )?;
        }
        let picker = MRMTransitionGroupPicker::from_param(&picker_param)?;
        let keys: Vec<String> = transition_group_map.keys().cloned().collect();
        for key in keys {
            let Some(transition_group) = transition_group_map.get_mut(&key) else {
                continue;
            };
            if transition_group.chromatograms().is_empty()
                || transition_group.transitions().is_empty()
            {
                continue;
            }
            {
                let mut ids = self
                    .ids
                    .lock()
                    .map_err(|_| bad("the unique ID generator is poisoned"))?;
                picker.pick_transition_group(transition_group, &mut ids)?;
            }
            self.score_peakgroups(transition_group, trafo, swath_maps, output, false, None)?;
        }
        Ok(())
    }

    /// Prepare the internal mapping of compound ID to compound (source
    /// `prepareProteinPeptideMaps_`). Required before
    /// [`score_peakgroups`](Self::score_peakgroups); a repeated ID keeps its
    /// last compound.
    pub fn prepare_protein_peptide_maps(&mut self, transition_exp: &LightTargetedExperiment) {
        for compound in transition_exp.compounds() {
            self.peptide_ref_map
                .insert(compound.id.clone(), compound.clone());
        }
    }

    /// Map the chromatograms of `input` onto the transitions of
    /// `transition_exp` by identifier (source `mapExperimentToTransitionList`):
    /// the native ID of a chromatogram is the native ID of its transition.
    ///
    /// With a positive `rt_extraction_window` each chromatogram is cut to the
    /// de-normalised expected RT of its compound `+/-` the window (requires
    /// [`prepare_protein_peptide_maps`](Self::prepare_protein_peptide_maps)).
    /// The chromatogram and its transition are added to the group of the
    /// transition's peptide reference, which is created on first use.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] with `strict` when a transition has no
    /// chromatogram or a chromatogram has no points in the window; always when
    /// a group's transitions and chromatograms differ in number afterwards
    /// ("Could not match all transition to all chromatograms"), and when the
    /// compound of a windowed transition is unknown.
    pub fn map_experiment_to_transition_list(
        &self,
        input: &SpectrumAccessPtr,
        transition_exp: &LightTargetedExperiment,
        transition_group_map: &mut TransitionGroupMapType,
        mut trafo: TransformationDescription,
        rt_extraction_window: f64,
    ) -> Result<()> {
        trafo.invert()?;
        let nr_chromatograms = input.nr_chromatograms();
        let mut chromatogram_map: BTreeMap<String, usize> = BTreeMap::new();
        for i in 0..nr_chromatograms {
            chromatogram_map.insert(input.chromatogram_native_id(i)?, i);
        }
        // Iterate through all transitions and store the transition with the
        // corresponding chromatogram in the corresponding transition group
        for transition in transition_exp.transitions() {
            let Some(&index) = chromatogram_map.get(&transition.transition_name) else {
                if self.strict {
                    return Err(bad(format!(
                        "Error: Transition {} from group {} does not have a corresponding chromatogram",
                        transition.transition_name, transition.peptide_ref
                    )));
                }
                continue;
            };
            // Retrieve the chromatogram and filter it by the desired RT
            let cptr = input.chromatogram_by_id(index)?;
            let mut chromatogram = MSChromatogram::default();
            // The expected RT is transformed back to the real RT space (the
            // transformation was inverted above) and the window taken around it.
            if rt_extraction_window > 0.0 {
                let compound = self
                    .peptide_ref_map
                    .get(&transition.peptide_ref)
                    .ok_or_else(|| {
                        bad(format!(
                            "the compound '{}' of transition '{}' is unknown; call prepare_protein_peptide_maps first",
                            transition.peptide_ref, transition.transition_name
                        ))
                    })?;
                let de_normalized_experimental_rt = trafo.apply(compound.rt)?;
                let rt_max = de_normalized_experimental_rt + rt_extraction_window;
                let rt_min = de_normalized_experimental_rt - rt_extraction_window;
                convert_to_openms_chromatogram_filter(&mut chromatogram, &cptr, rt_min, rt_max)?;
            } else {
                convert_to_openms_chromatogram(&cptr, &mut chromatogram)?;
            }
            // Check for empty chromatograms (e.g. the RT transformation is off)
            if chromatogram.peaks.is_empty() && self.strict {
                return Err(bad(format!(
                    "Error: Could not find any points for chromatogram {}. Maybe your retention time transformation is off?",
                    transition.transition_name
                )));
            }
            chromatogram.metadata.insert(
                "product_mz".to_owned(),
                MetaValue::source_float(transition.product_mz),
            );
            chromatogram.metadata.insert(
                "precursor_mz".to_owned(),
                MetaValue::source_float(transition.precursor_mz),
            );
            chromatogram.precursor = Precursor {
                mz: transition.precursor_mz,
                ..Precursor::default()
            };
            chromatogram.product.mz = transition.product_mz;
            chromatogram.native_id = transition.transition_name.clone();

            // Create a new transition group if there is none for this peptide
            let transition_group = transition_group_map
                .entry(transition.peptide_ref.clone())
                .or_insert_with(|| {
                    let mut group = MRMTransitionGroupType::new();
                    group.set_transition_group_id(transition.peptide_ref.clone());
                    group
                });
            // Now add the transition and the chromatogram to the group
            transition_group.add_transition(transition.clone(), &transition.transition_name)?;
            let id = chromatogram.native_id.clone();
            transition_group.add_chromatogram(chromatogram, &id)?;
        }
        // The assumption is that each transition has exactly one chromatogram
        for group in transition_group_map.values() {
            if !group.is_internally_consistent() {
                return Err(bad(format!(
                    "Error: Could not match all transition to all chromatograms:\nFor chromatogram {} I found {} chromatograms but {} transitions.",
                    group.transition_group_id(),
                    group.chromatograms().len(),
                    group.transitions().len()
                )));
            }
        }
        Ok(())
    }

    /// Split the transitions of a group into the detecting ones (source
    /// `splitTransitionGroupsDetection_`): a copy of the whole group when all
    /// transitions detect, otherwise the subset of the detecting transitions.
    fn split_detection(&self, group: &MRMTransitionGroupType) -> Result<MRMTransitionGroupType> {
        let detecting: Vec<String> = group
            .transitions()
            .iter()
            .filter(|t| t.is_detecting_transition())
            .map(|t| t.native_id().to_owned())
            .collect();
        if detecting.len() == group.transitions().len() {
            Ok(group.clone())
        } else {
            group.subset(&detecting)
        }
    }

    /// Split the identifying transitions of a group into target and decoy
    /// groups (source `splitTransitionGroupsIdentification_`).
    fn split_identification(
        &self,
        group: &MRMTransitionGroupType,
    ) -> Result<(MRMTransitionGroupType, MRMTransitionGroupType)> {
        let mut identifying = Vec::new();
        let mut identifying_decoy = Vec::new();
        for t in group.transitions() {
            if t.is_identifying_transition() {
                if t.decoy() {
                    identifying_decoy.push(t.native_id().to_owned());
                } else {
                    identifying.push(t.native_id().to_owned());
                }
            }
        }
        Ok((
            group.subset_dependent(&identifying)?,
            group.subset_dependent(&identifying_decoy)?,
        ))
    }

    #[allow(clippy::too_many_arguments)]
    fn score_identification(
        &self,
        trgr_ident: &MRMTransitionGroupType,
        trgr_detect: &MRMTransitionGroupType,
        scorer: &OpenSwathScoring,
        feature_idx: usize,
        feature_id: i64,
        native_ids_detection: &[String],
        signal_noise_estimators_identification: &[SignalToNoisePtr],
        det_intensity_ratio_score: f64,
        det_mi_ratio_score: f64,
        swath_maps: &[SwathMap],
        drift_target: f64,
        im_range: &RangeBase,
        mut sink: Option<&mut dyn MobilogramSink>,
    ) -> Result<OpenSwathIndScores> {
        let idmrmfeature: &MRMFeature =
            trgr_ident.features().get(feature_idx).ok_or_else(|| {
                bad("the identification group holds fewer features than the detection group")
            })?;
        let transitions_identification = trgr_ident.transitions();
        if signal_noise_estimators_identification.len() != transitions_identification.len() {
            return Err(bad(
                "Error: Identification signal-to-noise estimators do not match the transition group.",
            ));
        }
        let mut native_ids_identification: Vec<String> = Vec::new();
        let mut selected: Vec<SignalToNoisePtr> = Vec::new();
        for (i, transition) in transitions_identification.iter().enumerate() {
            let native_id = transition.native_id();
            let snptr = &signal_noise_estimators_identification[i];
            let intensity = f64::from(idmrmfeature.feature(native_id)?.base.intensity);
            if snptr.value_at_rt(idmrmfeature.feature.base.rt) > self.uis_threshold_sn
                && intensity > self.uis_threshold_peak_area
            {
                selected.push(Arc::clone(snptr));
                native_ids_identification.push(native_id.to_owned());
            }
        }
        let mut idscores = OpenSwathIndScores::default();
        if native_ids_identification.is_empty() {
            return Ok(idscores);
        }
        {
            let idimrmfeature = MRMFeatureOpenMS::new(idmrmfeature);
            // chromatographic ID scores
            scorer.calculate_chromatographic_id_scores(
                &idimrmfeature,
                &native_ids_identification,
                native_ids_detection,
                &selected,
                &mut idscores,
            )?;
            let ind_mi_score = idscores.ind_mi_score.clone();
            for (i, native_id) in native_ids_identification.iter().enumerate() {
                idscores.ind_transition_names.push(native_id.clone());
                let f = idmrmfeature.feature(native_id)?;
                if f.base.intensity > 0.0 {
                    let intensity_score = f64::from(f.base.intensity) / float_meta(f, "total_xic")?;
                    let mut intensity_ratio = 0.0;
                    if det_intensity_ratio_score > 0.0 {
                        intensity_ratio = intensity_score / det_intensity_ratio_score;
                    }
                    if intensity_ratio > 1.0 {
                        intensity_ratio = 1.0 / intensity_ratio;
                    }
                    let mut total_mi = 0.0;
                    if self.su.use_total_mi_score {
                        total_mi = float_meta(f, "total_mi")?;
                    }
                    let mut mi_ratio = 0.0;
                    if self.su.use_mi_score && self.su.use_total_mi_score {
                        if det_mi_ratio_score > 0.0 && total_mi > 0.0 {
                            let mi = ind_mi_score.get(i).copied().unwrap_or(0.0);
                            mi_ratio = (mi / total_mi) / det_mi_ratio_score;
                        }
                        if mi_ratio > 1.0 {
                            mi_ratio = 1.0 / mi_ratio;
                        }
                    }
                    idscores
                        .ind_area_intensity
                        .push(f64::from(f.base.intensity));
                    idscores
                        .ind_total_area_intensity
                        .push(float_meta(f, "total_xic")?);
                    idscores.ind_intensity_score.push(intensity_score);
                    idscores
                        .ind_apex_intensity
                        .push(float_meta(f, "peak_apex_int")?);
                    idscores
                        .ind_apex_position
                        .push(float_meta(f, "peak_apex_position")?);
                    idscores.ind_fwhm.push(float_meta(f, "width_at_50")?);
                    idscores.ind_total_mi.push(total_mi);
                    // `std::log(float)` is the C library's `logf`
                    idscores.ind_log_intensity.push(f64::from(
                        crate::analysis::featurefinder::glibc_logf::logf(f.base.intensity),
                    ));
                    idscores.ind_intensity_ratio.push(intensity_ratio);
                    idscores.ind_mi_ratio.push(mi_ratio);
                    if self.su.use_peak_shape_metrics {
                        idscores
                            .ind_start_position_at_5
                            .push(float_meta(f, "start_position_at_5")?);
                        idscores
                            .ind_end_position_at_5
                            .push(float_meta(f, "end_position_at_5")?);
                        idscores
                            .ind_start_position_at_10
                            .push(float_meta(f, "start_position_at_10")?);
                        idscores
                            .ind_end_position_at_10
                            .push(float_meta(f, "end_position_at_10")?);
                        idscores
                            .ind_start_position_at_50
                            .push(float_meta(f, "start_position_at_50")?);
                        idscores
                            .ind_end_position_at_50
                            .push(float_meta(f, "end_position_at_50")?);
                        idscores.ind_total_width.push(float_meta(f, "total_width")?);
                        idscores
                            .ind_tailing_factor
                            .push(float_meta(f, "tailing_factor")?);
                        idscores
                            .ind_asymmetry_factor
                            .push(float_meta(f, "asymmetry_factor")?);
                        idscores
                            .ind_slope_of_baseline
                            .push(float_meta(f, "slope_of_baseline")?);
                        idscores
                            .ind_baseline_delta_2_height
                            .push(float_meta(f, "baseline_delta_2_height")?);
                        idscores
                            .ind_points_across_baseline
                            .push(float_meta(f, "points_across_baseline")?);
                        idscores
                            .ind_points_across_half_height
                            .push(float_meta(f, "points_across_half_height")?);
                    }
                } else {
                    idscores.ind_area_intensity.push(0.0);
                    idscores.ind_total_area_intensity.push(0.0);
                    idscores.ind_intensity_score.push(0.0);
                    idscores.ind_apex_intensity.push(0.0);
                    idscores.ind_apex_position.push(0.0);
                    idscores.ind_fwhm.push(0.0);
                    idscores.ind_total_mi.push(0.0);
                    idscores.ind_log_intensity.push(0.0);
                    idscores.ind_intensity_ratio.push(0.0);
                    idscores.ind_mi_ratio.push(0.0);
                    if self.su.use_peak_shape_metrics {
                        idscores.ind_start_position_at_5.push(0.0);
                        idscores.ind_end_position_at_5.push(0.0);
                        idscores.ind_start_position_at_10.push(0.0);
                        idscores.ind_end_position_at_10.push(0.0);
                        idscores.ind_start_position_at_50.push(0.0);
                        idscores.ind_end_position_at_50.push(0.0);
                        idscores.ind_total_width.push(0.0);
                        idscores.ind_tailing_factor.push(0.0);
                        idscores.ind_asymmetry_factor.push(0.0);
                        idscores.ind_slope_of_baseline.push(0.0);
                        idscores.ind_baseline_delta_2_height.push(0.0);
                        idscores.ind_points_across_baseline.push(0.0);
                        idscores.ind_points_across_half_height.push(0.0);
                    }
                }
            }
            idscores.ind_num_transitions = i32::try_from(native_ids_identification.len())
                .map_err(|_| bad("identification transition count does not fit an int"))?;

            // DIA scores only on the identification transitions
            let swath_present = swath_maps
                .first()
                .and_then(|m| m.sptr.as_ref())
                .is_some_and(|s| s.nr_spectra() > 0);
            if swath_present && self.su.use_dia_scores {
                for native_id in &native_ids_identification {
                    let mut tmp_scores = OpenSwathScores::default();
                    let transition = trgr_ident.transition(native_id)?;
                    scorer.calculate_dia_id_scores(
                        &idimrmfeature,
                        transition,
                        trgr_detect.transitions(),
                        swath_maps,
                        im_range,
                        &self.diascoring,
                        &mut tmp_scores,
                        drift_target,
                        reborrow(&mut sink),
                        feature_id,
                    )?;
                    idscores
                        .ind_isotope_correlation
                        .push(tmp_scores.isotope_correlation);
                    idscores
                        .ind_isotope_overlap
                        .push(tmp_scores.isotope_overlap);
                    idscores.ind_massdev_score.push(tmp_scores.massdev_score);
                    // ion mobility scores
                    idscores.ind_im_drift.push(tmp_scores.im_drift);
                    idscores.ind_im_drift_left.push(tmp_scores.im_drift_left);
                    idscores.ind_im_drift_right.push(tmp_scores.im_drift_right);
                    idscores.ind_im_delta.push(tmp_scores.im_delta);
                    idscores.ind_im_delta_score.push(tmp_scores.im_delta_score);
                    idscores
                        .ind_im_log_intensity
                        .push(tmp_scores.im_log_intensity);
                    idscores
                        .ind_im_contrast_coelution
                        .push(tmp_scores.im_ind_contrast_coelution);
                    idscores
                        .ind_im_contrast_shape
                        .push(tmp_scores.im_ind_contrast_shape);
                    idscores
                        .ind_im_sum_contrast_coelution
                        .push(tmp_scores.im_ind_sum_contrast_coelution);
                    idscores
                        .ind_im_sum_contrast_shape
                        .push(tmp_scores.im_ind_sum_contrast_shape);
                }
            }
        }
        Ok(idscores)
    }

    /// Score all peak groups of a transition group (source `scorePeakgroups`).
    ///
    /// Every peak group found along the chromatograms of the group is scored:
    /// chromatographic and library scores, the DIA scores with a non-empty
    /// SWATH map, the identification (UIS) scores, the elution model fit
    /// score, the LDA prescores that become the overall quality, then the
    /// score metavalues, the peptide identification and the subordinate
    /// features are written. The scored features are appended to `output`
    /// ordered by quality, best first (with `stop_report_after_feature` only
    /// the best ones). Afterwards `transition_group` is replaced by its
    /// detection part holding the scored peak groups.
    ///
    /// With `ms1only` only the precursor chromatograms are scored.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when
    /// [`prepare_protein_peptide_maps`](Self::prepare_protein_peptide_maps)
    /// was not called or the group has no ID ("Error: Peptide reference map is
    /// empty", "Error: Transition group id is empty"), when the group's compound
    /// is unknown, when a group has no chromatograms and `ms1only` is false, and
    /// the errors of the scores.
    pub fn score_peakgroups(
        &self,
        transition_group: &mut MRMTransitionGroupType,
        trafo: &TransformationDescription,
        swath_maps: &[SwathMap],
        output: &mut FeatureMap,
        ms1only: bool,
        mut mobilogram_consumer: Option<&mut dyn MobilogramSink>,
    ) -> Result<()> {
        if self.peptide_ref_map.is_empty() {
            return Err(bad(
                "Error: Peptide reference map is empty, please call prepareProteinPeptideMaps_ first.",
            ));
        }
        if transition_group.transition_group_id().is_empty() {
            return Err(bad("Error: Transition group id is empty, please set it."));
        }
        let mut group_detection = self.split_detection(transition_group)?;
        let (group_identification, group_identification_decoy) = if self.su.use_uis_scores {
            self.split_identification(transition_group)?
        } else {
            (MRMTransitionGroupType::new(), MRMTransitionGroupType::new())
        };

        let has_feature_report_limit = self.stop_report_after_feature >= 0;

        // the drift time window of the precursor of the first chromatogram
        let mut im_range = RangeBase::new();
        let mut drift_target = 0.0;
        let drift_range = |prec: &Precursor| -> Result<(f64, RangeBase)> {
            let drift_time = prec.drift_time.unwrap_or(-1.0);
            if drift_time > 0.0 {
                let lower = drift_time - prec.drift_window_lower_offset;
                let upper = drift_time + prec.drift_window_upper_offset;
                Ok((drift_time, RangeBase::from_min_max(lower, upper)?))
            } else {
                Ok((drift_time, RangeBase::new()))
            }
        };
        if let Some(chromatogram) = group_detection.chromatograms().first() {
            let (target, range) = drift_range(&chromatogram.precursor)?;
            drift_target = target;
            im_range = range;
        } else if let Some(chromatogram) = group_detection.precursor_chromatograms().first() {
            let (target, range) = drift_range(&chromatogram.precursor)?;
            drift_target = target;
            im_range = range;
        }

        let sn_estimator = |chromatogram: &MSChromatogram| -> Result<SignalToNoisePtr> {
            Ok(Arc::new(SignalToNoiseOpenMS::from_chromatogram(
                chromatogram,
                self.sn_win_len,
                self.sn_bin_count as usize,
                self.write_log_messages,
            )?))
        };
        let mut signal_noise_estimators: Vec<SignalToNoisePtr> = Vec::new();
        for chromatogram in group_detection.chromatograms() {
            signal_noise_estimators.push(sn_estimator(chromatogram)?);
        }
        // skip the MS1 noise estimator when scoring fragment ions
        let mut ms1_signal_noise_estimators: Vec<SignalToNoisePtr> = Vec::new();
        if ms1only {
            for chromatogram in group_detection.precursor_chromatograms() {
                ms1_signal_noise_estimators.push(sn_estimator(chromatogram)?);
            }
        }

        // the expected RT of this compound, in the experimental RT space
        let pep = self
            .peptide_ref_map
            .get(group_detection.transition_group_id())
            .ok_or_else(|| {
                bad(format!(
                    "the compound '{}' is not in the peptide reference map",
                    group_detection.transition_group_id()
                ))
            })?
            .clone();
        let mut newtr = trafo.clone();
        newtr.invert()?;
        let expected_rt = newtr.apply(pep.rt)?;

        let mut scorer = OpenSwathScoring::new();
        scorer.initialize(
            self.rt_normalization_factor,
            self.add_up_spectra,
            self.spacing_for_spectra_resampling,
            self.merge_spectra_by_peak_width_fraction,
            self.im_extra_drift,
            &self.su,
            &self.spectrum_addition_method,
            &self.spectrum_merge_method_type,
            self.use_ms1_ion_mobility,
            self.apply_im_peak_picking,
        )?;
        let mut pd = ProteaseDigestion::default();
        pd.set_enzyme(&self.enzyme)?;

        let n_features = group_detection.features().len();
        let defer_feature_output = has_feature_report_limit
            && usize::try_from(self.stop_report_after_feature)
                .is_ok_and(|limit| limit < n_features);

        // group-invariant values
        let swath_present = swath_maps
            .first()
            .and_then(|m| m.sptr.as_ref())
            .is_some_and(|s| s.nr_spectra() > 0);
        let detection_precursor_mz = group_detection
            .transitions()
            .first()
            .map_or(-1.0, LightTransition::precursor_mz);
        let group_size = group_detection.size();
        if group_size == 0 && !ms1only {
            return Err(bad(format!(
                "Error: Transition group {} has no chromatograms.",
                group_detection.transition_group_id()
            )));
        }

        // cache
        let cache = {
            let transitions = group_detection.transitions();
            let mut normalized: Vec<f64> = transitions
                .iter()
                .map(|t| 0.0_f64.max(t.library_intensity()))
                .collect();
            if !normalized.is_empty() {
                crate::analysis::openswath::scoring::normalize_sum(&mut normalized);
            }
            TransitionGroupCache {
                normalized_library_intensity: normalized,
                transition_native_ids: transitions
                    .iter()
                    .map(|t| t.native_id().to_owned())
                    .collect(),
                precursor_ids: group_detection
                    .precursor_chromatograms()
                    .iter()
                    .map(|c| c.native_id.clone())
                    .collect(),
            }
        };

        // identification signal-to-noise estimators, once per group
        let ident_estimators = |group: &MRMTransitionGroupType| -> Result<Vec<SignalToNoisePtr>> {
            group
                .transitions()
                .iter()
                .map(|t| sn_estimator(group.chromatogram(t.native_id())?))
                .collect()
        };
        let (sn_identification, sn_identification_decoy) = if self.su.use_uis_scores {
            (
                ident_estimators(&group_identification)?,
                ident_estimators(&group_identification_decoy)?,
            )
        } else {
            (Vec::new(), Vec::new())
        };

        // The features are taken out of the group while they are scored, so
        // the group can be read (elution model, identification) meanwhile.
        let mut mrmfeatures: Vec<MRMFeature> = std::mem::take(group_detection.features_mut());
        let mut feature_list: Vec<MRMFeature> = Vec::new();
        let mut scored_feature_indices: Vec<usize> = Vec::new();

        let scoring_result = (|| -> Result<()> {
            for (feature_idx, current) in mrmfeatures.iter_mut().enumerate() {
                if current.feature.base.unique_id == 0 {
                    current.feature.base.unique_id = self.ids_next();
                }
                let feature_id = clear_sign_bit(current.feature.base.unique_id);
                let precursor_mz: f64;
                let mut scores = *current.scores();
                let xx_lda_prescore;
                let mut overall_quality;
                let mut masserror_ppm: Option<Vec<f64>> = None;
                let mut id_target: Option<OpenSwathIndScores> = None;
                let mut id_decoy: Option<OpenSwathIndScores> = None;
                if ms1only {
                    let imrmfeature = MRMFeatureOpenMS::new(&*current);
                    precursor_mz = current.feature.base.mz;
                    // S/N scores
                    scores.sn_ratio =
                        MRMScoring::calc_sn_score(&imrmfeature, &ms1_signal_noise_estimators);
                    // everything below S/N 1 is set to zero (the log is safe)
                    scores.log_sn_score = if scores.sn_ratio < 1.0 {
                        0.0
                    } else {
                        log(scores.sn_ratio)
                    };
                    // RT scores: the score is the delta iRT
                    let normalized_experimental_rt = trafo.apply(current.feature.base.rt)?;
                    let rt_score = MRMScoring::calc_rt_score(&pep, normalized_experimental_rt);
                    scores.normalized_experimental_rt = normalized_experimental_rt;
                    scores.raw_rt_score = rt_score;
                    scores.norm_rt_score = rt_score / self.rt_normalization_factor;
                    // full spectra scores
                    if let Some(ms1) = self.ms1_map.as_ref().filter(|m| m.nr_spectra() > 0) {
                        if current.feature.base.mz > 0.0 {
                            scorer.calculate_precursor_dia_scores(
                                Some(ms1),
                                &self.diascoring,
                                precursor_mz,
                                current.feature.base.rt,
                                &pep,
                                im_range,
                                &mut scores,
                            )?;
                        }
                    }
                    xx_lda_prescore = self.lda_prescore(&scores);
                    overall_quality = xx_lda_prescore as f32;
                } else {
                    let imrmfeature = MRMFeatureOpenMS::with_ids(
                        &*current,
                        &cache.transition_native_ids,
                        &cache.precursor_ids,
                    )?;
                    // library and chromatographic scores
                    scorer.calculate_chromatographic_scores(
                        &imrmfeature,
                        &cache.transition_native_ids,
                        &cache.precursor_ids,
                        &cache.normalized_library_intensity,
                        &signal_noise_estimators,
                        &mut scores,
                    )?;
                    let normalized_experimental_rt = trafo.apply(current.feature.base.rt)?;
                    scorer.calculate_library_scores(
                        &imrmfeature,
                        group_detection.transitions(),
                        &pep,
                        normalized_experimental_rt,
                        &mut scores,
                    )?;
                    // DIA scores
                    if swath_present && self.su.use_dia_scores {
                        masserror_ppm = Some(scorer.calculate_dia_scores(
                            &imrmfeature,
                            group_detection.transitions(),
                            &cache.normalized_library_intensity,
                            swath_maps,
                            self.ms1_map.as_ref(),
                            &self.diascoring,
                            &pep,
                            &mut scores,
                            drift_target,
                            &im_range,
                            reborrow(&mut mobilogram_consumer),
                            feature_id,
                        )?);
                    }
                    let mut det_intensity_ratio_score = 0.0;
                    let total_xic = float_meta(&current.feature, "total_xic")?;
                    if total_xic > 0.0 {
                        det_intensity_ratio_score =
                            f64::from(current.feature.base.intensity) / total_xic;
                    }
                    // mutual information scores
                    let mut det_mi_ratio_score = 0.0;
                    if self.su.use_mi_score && self.su.use_total_mi_score {
                        let total_mi = float_meta(&current.feature, "total_mi")?;
                        if total_mi > 0.0 {
                            det_mi_ratio_score = scores.mi_score / total_mi;
                        }
                    }
                    // unique ion signature (UIS) scores
                    if self.su.use_uis_scores {
                        if !group_identification.transitions().is_empty() {
                            id_target = Some(self.score_identification(
                                &group_identification,
                                &group_detection,
                                &scorer,
                                feature_idx,
                                feature_id,
                                &cache.transition_native_ids,
                                &sn_identification,
                                det_intensity_ratio_score,
                                det_mi_ratio_score,
                                swath_maps,
                                drift_target,
                                &im_range,
                                reborrow(&mut mobilogram_consumer),
                            )?);
                        }
                        if !group_identification_decoy.transitions().is_empty() {
                            id_decoy = Some(self.score_identification(
                                &group_identification_decoy,
                                &group_detection,
                                &scorer,
                                feature_idx,
                                feature_id,
                                &cache.transition_native_ids,
                                &sn_identification_decoy,
                                det_intensity_ratio_score,
                                det_mi_ratio_score,
                                swath_maps,
                                drift_target,
                                &im_range,
                                reborrow(&mut mobilogram_consumer),
                            )?);
                        }
                    }
                    if self.su.use_elution_model_score {
                        scores.elution_model_fit_score = self
                            .emgscoring
                            .calc_elution_fit_score(&*current, &group_detection)?;
                    }
                    xx_lda_prescore = self.lda_prescore(&scores);
                    overall_quality = xx_lda_prescore as f32;
                    // add the DIA scores and ion mobility scores
                    if swath_present && self.su.use_dia_scores {
                        let xx_swath_prescore = -scores.calculate_swath_lda_prescore();
                        overall_quality = xx_swath_prescore as f32;
                    }
                    precursor_mz = detection_precursor_mz;
                }
                let mrmfeature = current;
                mrmfeature.set_scores(scores);
                mrmfeature.feature.base.quality = overall_quality;
                if let Some(ppm) = masserror_ppm {
                    set_meta(
                        &mut mrmfeature.feature,
                        "masserror_ppm",
                        MetaValue::try_from(ppm)?,
                    );
                }
                if let Some(ids) = &id_target {
                    mrmfeature.id_scores_as_meta_value(false, ids)?;
                }
                if let Some(ids) = &id_decoy {
                    mrmfeature.id_scores_as_meta_value(true, ids)?;
                }
                if defer_feature_output {
                    scored_feature_indices.push(feature_idx);
                } else {
                    // add the peptide hit information to the feature
                    self.add_score_meta_values(
                        mrmfeature,
                        &group_detection,
                        &signal_noise_estimators,
                        &ms1_signal_noise_estimators,
                        expected_rt,
                        swath_present,
                        ms1only,
                    )?;
                    self.prepare_scored_feature_output(
                        mrmfeature,
                        &pep,
                        &pd,
                        swath_present,
                        precursor_mz,
                        ms1only,
                    )?;
                    feature_list.push(mrmfeature.clone());
                }
            }
            Ok(())
        })();

        if let Err(error) = scoring_result {
            *group_detection.features_mut() = mrmfeatures;
            return Err(error);
        }

        if defer_feature_output {
            // Order by quality before the expensive output preparation; only
            // the best candidates are reported.
            let report_count = scored_feature_indices
                .len()
                .min(usize::try_from(self.stop_report_after_feature).unwrap_or(0));
            source_sort_by(&mut scored_feature_indices, |&lhs, &rhs| {
                mrmfeatures[rhs].feature.base.quality < mrmfeatures[lhs].feature.base.quality
            })?;
            for &idx in scored_feature_indices.iter().take(report_count) {
                let selected_precursor_mz = if ms1only {
                    mrmfeatures[idx].feature.base.mz
                } else {
                    detection_precursor_mz
                };
                let selected = &mut mrmfeatures[idx];
                self.add_score_meta_values(
                    selected,
                    &group_detection,
                    &signal_noise_estimators,
                    &ms1_signal_noise_estimators,
                    expected_rt,
                    swath_present,
                    ms1only,
                )?;
                self.prepare_scored_feature_output(
                    selected,
                    &pep,
                    &pd,
                    swath_present,
                    selected_precursor_mz,
                    ms1only,
                )?;
                output.features.push(selected.feature.clone());
            }
        } else {
            // order by quality (high to low, via the reversed range)
            source_sort_reversed_by(&mut feature_list, |a, b| {
                a.feature.base.quality < b.feature.base.quality
            })?;
            for f in &feature_list {
                output.features.push(f.feature.clone());
            }
        }
        // store all data manipulation performed on the group's features
        *group_detection.features_mut() = mrmfeatures;
        *transition_group = group_detection;
        Ok(())
    }

    /// The next unique ID of the shared generator (the source's global
    /// `UniqueIdGenerator::getUniqueId`).
    fn ids_next(&self) -> u64 {
        match self.ids.lock() {
            Ok(mut generator) => generator.get_unique_id(),
            Err(poisoned) => poisoned.into_inner().get_unique_id(),
        }
    }

    /// `-calculate_lda_prescore`, or the single transition model.
    fn lda_prescore(&self, scores: &OpenSwathScores) -> f64 {
        if self.scoring_model == "single_transition" {
            -scores.calculate_lda_single_transition()
        } else {
            -scores.calculate_lda_prescore()
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn add_score_meta_values(
        &self,
        mrmfeature: &mut MRMFeature,
        group: &MRMTransitionGroupType,
        signal_noise_estimators: &[SignalToNoisePtr],
        ms1_signal_noise_estimators: &[SignalToNoisePtr],
        expected_rt: f64,
        swath_present: bool,
        ms1only: bool,
    ) -> Result<()> {
        let scores = *mrmfeature.scores();
        let su = &self.su;
        let rt = mrmfeature.feature.base.rt;
        let xx_lda_prescore = self.lda_prescore(&scores);
        if ms1only {
            if su.use_sn_score {
                add_score(mrmfeature, "sn_ratio", scores.sn_ratio);
                add_score(mrmfeature, "var_log_sn_score", scores.log_sn_score);
                // compute the subfeature log-SN values
                for (k, chromatogram) in group.precursor_chromatograms().iter().enumerate() {
                    let estimator = ms1_signal_noise_estimators.get(k).ok_or_else(|| {
                        bad("fewer MS1 noise estimators than precursor chromatograms")
                    })?;
                    let mut sn_value = estimator.value_at_rt(rt);
                    if sn_value < 1.0 {
                        sn_value = 1.0;
                    }
                    let f = mrmfeature.precursor_feature_mut(&chromatogram.native_id)?;
                    set_meta(f, "logSN", MetaValue::source_float(log(sn_value)));
                }
            }
            if su.use_rt_score {
                add_score(mrmfeature, "delta_rt", rt - expected_rt);
                add_score(mrmfeature, "assay_rt", expected_rt);
                add_score(mrmfeature, "norm_RT", scores.normalized_experimental_rt);
                add_score(mrmfeature, "rt_score", scores.raw_rt_score);
                add_score(mrmfeature, "var_norm_rt_score", scores.norm_rt_score);
            }
            if su.use_ms1_fullscan {
                add_score(mrmfeature, "var_ms1_ppm_diff", scores.ms1_ppm_score);
                add_score(
                    mrmfeature,
                    "var_ms1_isotope_correlation",
                    scores.ms1_isotope_correlation,
                );
                add_score(
                    mrmfeature,
                    "var_ms1_isotope_overlap",
                    scores.ms1_isotope_overlap,
                );
            }
            add_score(mrmfeature, "main_var_xx_lda_prelim_score", xx_lda_prescore);
            add_score(mrmfeature, "xx_lda_prelim_score", xx_lda_prescore);
            return Ok(());
        }
        if su.use_coelution_score {
            add_score(
                mrmfeature,
                "var_xcorr_coelution",
                scores.xcorr_coelution_score,
            );
            add_score(
                mrmfeature,
                "var_xcorr_coelution_weighted",
                scores.weighted_coelution_score,
            );
        }
        if su.use_shape_score {
            add_score(mrmfeature, "var_xcorr_shape", scores.xcorr_shape_score);
            add_score(
                mrmfeature,
                "var_xcorr_shape_weighted",
                scores.weighted_xcorr_shape,
            );
        }
        if su.use_library_score {
            add_score(mrmfeature, "var_library_corr", scores.library_corr);
            add_score(
                mrmfeature,
                "var_library_rmsd",
                scores.library_norm_manhattan,
            );
            add_score(mrmfeature, "var_library_sangle", scores.library_sangle);
            add_score(
                mrmfeature,
                "var_library_rootmeansquare",
                scores.library_rootmeansquare,
            );
            add_score(
                mrmfeature,
                "var_library_manhattan",
                scores.library_manhattan,
            );
            add_score(mrmfeature, "var_library_dotprod", scores.library_dotprod);
        }
        if su.use_rt_score {
            add_score(mrmfeature, "delta_rt", rt - expected_rt);
            add_score(mrmfeature, "assay_rt", expected_rt);
            add_score(mrmfeature, "norm_RT", scores.normalized_experimental_rt);
            add_score(mrmfeature, "rt_score", scores.raw_rt_score);
            add_score(mrmfeature, "var_norm_rt_score", scores.norm_rt_score);
        }
        // TODO of the source: do we really want these intensity scores?
        if su.use_intensity_score {
            let total_xic = float_meta(&mrmfeature.feature, "total_xic")?;
            if total_xic > 0.0 {
                let intensity = f64::from(mrmfeature.feature.base.intensity);
                add_score(mrmfeature, "var_intensity_score", intensity / total_xic);
            } else {
                add_score(mrmfeature, "var_intensity_score", 0.0);
            }
        }
        if su.use_total_xic_score {
            let total_xic = float_meta(&mrmfeature.feature, "total_xic")?;
            add_score(mrmfeature, "total_xic", total_xic);
        }
        if su.use_total_mi_score {
            let total_mi = float_meta(&mrmfeature.feature, "total_mi")?;
            add_score(mrmfeature, "total_mi", total_mi);
        }
        if su.use_nr_peaks_score {
            add_score(mrmfeature, "nr_peaks", scores.nr_peaks);
        }
        if su.use_sn_score {
            add_score(mrmfeature, "sn_ratio", scores.sn_ratio);
            add_score(mrmfeature, "var_log_sn_score", scores.log_sn_score);
            // compute the subfeature log-SN values
            for (k, chromatogram) in group.chromatograms().iter().enumerate() {
                let estimator = signal_noise_estimators
                    .get(k)
                    .ok_or_else(|| bad("fewer noise estimators than chromatograms"))?;
                let mut sn_value = estimator.value_at_rt(rt);
                if sn_value < 1.0 {
                    sn_value = 1.0;
                }
                let f = mrmfeature.feature_mut(&chromatogram.native_id)?;
                set_meta(f, "logSN", MetaValue::source_float(log(sn_value)));
            }
        }
        if su.use_mi_score {
            add_score(mrmfeature, "var_mi_score", scores.mi_score);
            add_score(
                mrmfeature,
                "var_mi_weighted_score",
                scores.weighted_mi_score,
            );
            if su.use_total_mi_score {
                let total_mi = float_meta(&mrmfeature.feature, "total_mi")?;
                if total_mi > 0.0 {
                    add_score(mrmfeature, "var_mi_ratio_score", scores.mi_score / total_mi);
                } else {
                    add_score(mrmfeature, "var_mi_ratio_score", 0.0);
                }
            }
        }
        if su.use_elution_model_score {
            add_score(
                mrmfeature,
                "var_elution_model_fit_score",
                scores.elution_model_fit_score,
            );
        }
        if !swath_present {
            add_score(mrmfeature, "main_var_xx_lda_prelim_score", xx_lda_prescore);
        }
        add_score(mrmfeature, "xx_lda_prelim_score", xx_lda_prescore);
        if swath_present && su.use_dia_scores {
            if su.use_ms2_isotope_scores {
                add_score(
                    mrmfeature,
                    "var_isotope_correlation_score",
                    scores.isotope_correlation,
                );
                add_score(
                    mrmfeature,
                    "var_isotope_overlap_score",
                    scores.isotope_overlap,
                );
            }
            add_score(mrmfeature, "var_massdev_score", scores.massdev_score);
            add_score(
                mrmfeature,
                "var_massdev_score_weighted",
                scores.weighted_massdev_score,
            );
            if su.use_ionseries_scores {
                add_score(mrmfeature, "var_bseries_score", scores.bseries_score);
                add_score(mrmfeature, "var_yseries_score", scores.yseries_score);
            }
            if su.use_ms2_isotope_scores {
                add_score(mrmfeature, "var_dotprod_score", scores.dotprod_score_dia);
                add_score(mrmfeature, "var_manhatt_score", scores.manhatt_score_dia);
            }
            if su.use_ms1_correlation {
                if scores.ms1_xcorr_shape_score > -1.0 {
                    add_score(
                        mrmfeature,
                        "var_ms1_xcorr_shape",
                        scores.ms1_xcorr_shape_score,
                    );
                }
                if scores.ms1_xcorr_coelution_score > -1.0 {
                    add_score(
                        mrmfeature,
                        "var_ms1_xcorr_coelution",
                        scores.ms1_xcorr_coelution_score,
                    );
                }
                add_score(
                    mrmfeature,
                    "var_ms1_xcorr_shape_contrast",
                    scores.ms1_xcorr_shape_contrast_score,
                );
                add_score(
                    mrmfeature,
                    "var_ms1_xcorr_shape_combined",
                    scores.ms1_xcorr_shape_combined_score,
                );
                add_score(
                    mrmfeature,
                    "var_ms1_xcorr_coelution_contrast",
                    scores.ms1_xcorr_coelution_contrast_score,
                );
                add_score(
                    mrmfeature,
                    "var_ms1_xcorr_coelution_combined",
                    scores.ms1_xcorr_coelution_combined_score,
                );
            }
            if su.use_ms1_mi {
                if scores.ms1_mi_score > -1.0 {
                    add_score(mrmfeature, "var_ms1_mi_score", scores.ms1_mi_score);
                }
                add_score(
                    mrmfeature,
                    "var_ms1_mi_contrast_score",
                    scores.ms1_mi_contrast_score,
                );
                add_score(
                    mrmfeature,
                    "var_ms1_mi_combined_score",
                    scores.ms1_mi_combined_score,
                );
            }
            if su.use_ms1_fullscan {
                add_score(mrmfeature, "var_ms1_ppm_diff", scores.ms1_ppm_score);
                add_score(
                    mrmfeature,
                    "var_ms1_isotope_correlation",
                    scores.ms1_isotope_correlation,
                );
                add_score(
                    mrmfeature,
                    "var_ms1_isotope_overlap",
                    scores.ms1_isotope_overlap,
                );
            }
            let xx_swath_prescore = -scores.calculate_swath_lda_prescore();
            add_score(
                mrmfeature,
                "main_var_xx_swath_prelim_score",
                xx_swath_prescore,
            );
        }
        if swath_present && su.use_im_scores {
            add_score(
                mrmfeature,
                "var_im_xcorr_shape",
                scores.im_xcorr_shape_score,
            );
            add_score(
                mrmfeature,
                "var_im_xcorr_coelution",
                scores.im_xcorr_coelution_score,
            );
            add_score(mrmfeature, "var_im_delta_score", scores.im_delta_score);
            add_score(
                mrmfeature,
                "var_im_ms1_delta_score",
                scores.im_ms1_delta_score,
            );
            add_score(mrmfeature, "im_drift", scores.im_drift); // MS2 level
            add_score(mrmfeature, "im_drift_left", scores.im_drift_left); // MS2 level
            add_score(mrmfeature, "im_drift_right", scores.im_drift_right); // MS2 level
            add_score(mrmfeature, "im_drift_weighted", scores.im_drift_weighted); // MS2 level
            add_score(mrmfeature, "im_log_intensity", scores.im_log_intensity); // MS2 level
            add_score(mrmfeature, "im_ms1_drift", scores.im_ms1_drift); // MS1 level
            add_score(mrmfeature, "im_ms1_delta", scores.im_ms1_delta); // MS1 level
            add_score(mrmfeature, "im_delta", scores.im_delta); // MS2 level
        }
        Ok(())
    }

    fn prepare_scored_feature_output(
        &self,
        mrmfeature: &mut MRMFeature,
        pep: &LightCompound,
        pd: &ProteaseDigestion,
        swath_present: bool,
        precursor_mz: f64,
        ms1only: bool,
    ) -> Result<()> {
        let mut pep_id = PeptideIdentification::default();
        let mut pep_hit = PeptideHit::default();
        if pep.charge_state() != 0 {
            pep_hit.charge = pep.charge_state();
        }
        let xx_lda_prescore = if mrmfeature
            .feature
            .base
            .metadata
            .contains_key("xx_lda_prelim_score")
        {
            float_meta(&mrmfeature.feature, "xx_lda_prelim_score")?
        } else {
            f64::from(mrmfeature.feature.base.quality)
        };
        pep_hit.score = xx_lda_prescore;
        if swath_present
            && mrmfeature
                .feature
                .base
                .metadata
                .contains_key("xx_swath_prelim_score")
        {
            pep_hit.score = float_meta(&mrmfeature.feature, "xx_swath_prelim_score")?;
        }
        if pep.is_peptide() && !pep.sequence.is_empty() {
            let sequence = AASequence::parse(&pep.sequence)?;
            let count = pd.peptide_count(&sequence)? as i64;
            pep_hit.sequence = sequence;
            set_meta(
                &mut mrmfeature.feature,
                "missedCleavages",
                MetaValue::from(count - 1),
            );
        }
        // set the protein accession numbers
        for reference in &pep.protein_refs {
            pep_hit.evidences.push(PeptideEvidence {
                protein_accession: reference.clone(),
                ..PeptideEvidence::default()
            });
        }
        pep_id.hits.push(pep_hit);
        pep_id.identifier = RUN_IDENTIFIER.to_owned();
        mrmfeature.feature.base.peptide_identifications.push(pep_id);
        if mrmfeature.feature.base.unique_id == 0 {
            mrmfeature.feature.base.unique_id = self.ids_next();
        }
        set_meta(
            &mut mrmfeature.feature,
            "PrecursorMZ",
            MetaValue::source_float(precursor_mz),
        );
        self.prepare_feature_output(mrmfeature, ms1only, pep.charge_state())?;
        set_meta(
            &mut mrmfeature.feature,
            "xx_swath_prelim_score",
            MetaValue::source_float(0.0),
        );
        Ok(())
    }

    fn process_feature_for_output(
        &self,
        curr_feature: &mut Feature,
        total_intensity: &mut f64,
        total_peak_apices: &mut f64,
        ms_level: &str,
    ) -> Result<()> {
        // Save some space when writing out the featureXML
        if !self.write_convex_hull {
            curr_feature.convex_hulls.clear();
        }
        // Ensure a unique id is present
        if curr_feature.base.unique_id == 0 {
            curr_feature.base.unique_id = self.ids_next();
        }
        // Sum up the intensities of the features
        if curr_feature.base.mz > self.quantification_cutoff {
            *total_intensity += f64::from(curr_feature.base.intensity);
            *total_peak_apices += float_meta(curr_feature, "peak_apex_int")?;
        }
        set_meta(curr_feature, "FeatureLevel", MetaValue::from(ms_level));
        Ok(())
    }

    fn prepare_feature_output(
        &self,
        mrmfeature: &mut MRMFeature,
        ms1only: bool,
        charge: i32,
    ) -> Result<()> {
        // process all current features, then append all precursor subfeatures
        let mut all_features: Vec<Feature> = mrmfeature.features().to_vec();
        let (mut total_intensity, mut total_peak_apices) = (0.0, 0.0);
        let (mut ms1_total_intensity, mut ms1_total_peak_apices) = (0.0, 0.0);
        for f in &mut all_features {
            self.process_feature_for_output(
                f,
                &mut total_intensity,
                &mut total_peak_apices,
                "MS2",
            )?;
        }
        // also append the data of the MS1 precursors
        let precursor_ids: Vec<String> = mrmfeature
            .precursor_feature_ids()
            .map(str::to_owned)
            .collect();
        for id in &precursor_ids {
            let mut curr_feature = mrmfeature.precursor_feature(id)?.clone();
            if charge != 0 {
                curr_feature.base.charge = charge;
            }
            self.process_feature_for_output(
                &mut curr_feature,
                &mut ms1_total_intensity,
                &mut ms1_total_peak_apices,
                "MS1",
            )?;
            if ms1only {
                total_intensity += f64::from(curr_feature.base.intensity);
                total_peak_apices += float_meta(&curr_feature, "peak_apex_int")?;
            }
            all_features.push(curr_feature);
        }
        // add all the subfeatures as subordinates
        mrmfeature.feature.subordinates = all_features;
        // overwrite the reported intensities with those above the m/z cutoff
        mrmfeature.feature.base.intensity = total_intensity as f32;
        set_meta(
            &mut mrmfeature.feature,
            "peak_apices_sum",
            MetaValue::source_float(total_peak_apices),
        );
        set_meta(
            &mut mrmfeature.feature,
            "ms1_area_intensity",
            MetaValue::source_float(ms1_total_intensity),
        );
        set_meta(
            &mut mrmfeature.feature,
            "ms1_apex_intensity",
            MetaValue::source_float(ms1_total_peak_apices),
        );
        Ok(())
    }
}
