// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Selected/multiple reaction monitoring peak groups.
//!
//! Ports `KERNEL/MRMFeature.h` (with `KERNEL/MRMFeature.cpp`) and the
//! header-only template `KERNEL/MRMTransitionGroup.h`. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/MRM_SUPPORT.md)
//! for the full member mapping and evidence.
//!
//! A [`MRMFeature`](crate::kernel::mrm::MRMFeature) is a peak group: one
//! elution peak observed across several transition chromatograms at once. It
//! extends [`Feature`](crate::kernel::Feature) with the per-transition
//! features, the precursor features and the OpenSWATH peak-group score record.
//!
//! A [`MRMTransitionGroup`](crate::kernel::mrm::MRMTransitionGroup) collects
//! the transitions (assay metadata), their measured chromatograms, the
//! precursor chromatograms extracted from the MS1 map, and the peak groups
//! found across all of them. It keeps three parallel key/index maps, one per
//! list; the pair
//! [`is_internally_consistent`](crate::kernel::mrm::MRMTransitionGroup::is_internally_consistent)
//! and
//! [`chromatogram_ids_match`](crate::kernel::mrm::MRMTransitionGroup::chromatogram_ids_match)
//! is what checks that they agree.
//!
//! The group is generic over its transition type because the source's
//! `ReactionMonitoringTransition` (`ANALYSIS/MRM`) is not ported. The
//! [`Transition`](crate::kernel::mrm::Transition) trait carries only what the
//! group and its immediate consumers need;
//! [`SimpleTransition`](crate::kernel::mrm::SimpleTransition) is a concrete
//! implementation for tests and for callers with no assay library type of
//! their own. The chromatogram parameter is likewise generic over
//! [`NativeIdentified`](crate::kernel::mrm::NativeIdentified), because the
//! source header states that the group must also accept `MSSpectrum` as raw
//! data storage.
//!
//! The source contains no `#pragma omp`; both source and port are serial.

use crate::kernel::features::{BaseFeature, Feature};
use crate::kernel::{MSChromatogram, MSSpectrum};
use crate::metadata::{MetaValue, MetaValueData};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::ops::{Deref, DerefMut};

fn limit(what: &str) -> Error {
    Error::InvalidValue(format!("MRM {what} exceeds the configured ceiling"))
}

fn missing(what: &str, key: &str) -> Error {
    Error::MissingInformation(format!("no {what} registered under key '{key}'"))
}

/// What to do when a key is already registered in a lookup map.
///
/// The two source classes disagree: `MRMTransitionGroup::addTransition` and
/// both `add*Chromatogram` overloads throw `Exception::InvalidValue` on a
/// repeated key, while `MRMFeature::addFeature` and
/// `MRMFeature::addPrecursorFeature` overwrite the map entry and leave the
/// previously stored feature in the list but unreachable by any key. This enum
/// selects between the two for the feature lists; the transition-group lists
/// always reject, as their source does.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DuplicateKeyPolicy {
    /// Refuse the insertion with [`Error::InvalidValue`]; nothing is stored.
    /// This is the native default because the source alternative silently
    /// strands data.
    #[default]
    Reject,
    /// Reproduce `MRMFeature::addFeature`: append the new feature and point the
    /// key at it, leaving the previously keyed feature in the list without a
    /// key of its own.
    SourceOverwrite,
}

/// A value carrying a native identifier, as the source `getNativeID()`.
///
/// `MRMTransitionGroup` looks chromatograms up by a caller-supplied key and
/// separately compares that key against the stored native ID, so the raw data
/// type must expose one. Implemented here for
/// [`MSChromatogram`] and [`MSSpectrum`], the two types the source header
/// names.
pub trait NativeIdentified {
    /// The stored native identifier, empty when unset.
    fn native_id(&self) -> &str;
}

impl NativeIdentified for MSChromatogram {
    fn native_id(&self) -> &str {
        &self.native_id
    }
}

impl NativeIdentified for MSSpectrum {
    fn native_id(&self) -> &str {
        &self.native_id
    }
}

/// The assay metadata a transition group needs from one transition.
///
/// The source class is `ReactionMonitoringTransition` in `ANALYSIS/MRM`, which
/// is not ported; this trait is the minimal replacement. Only
/// [`Transition::native_id`] and [`Transition::library_intensity`] are read by
/// `MRMTransitionGroup` itself (in `subset`, `subsetDependent` and
/// `getLibraryIntensity`). The precursor and product m/z and the three
/// selection flags are carried because every OpenSWATH consumer of a group
/// reads them straight off the transitions it holds, and a trait without them
/// would force those callers back to a concrete type.
pub trait Transition {
    /// Unique identifier of the transition; the key `addTransition` expects.
    fn native_id(&self) -> &str;
    /// Precursor (Q1) mass-to-charge ratio.
    fn precursor_mz(&self) -> f64;
    /// Product (Q3) mass-to-charge ratio.
    fn product_mz(&self) -> f64;
    /// Expected relative intensity from the assay library. The source allows
    /// negative values here; `getLibraryIntensity` clamps them, this trait does
    /// not.
    fn library_intensity(&self) -> f64;
    /// Whether the transition takes part in peak-group detection.
    fn is_detecting(&self) -> bool;
    /// Whether the transition is used for peptidoform identification (IPF).
    fn is_identifying(&self) -> bool;
    /// Whether the transition contributes to quantification.
    fn is_quantifying(&self) -> bool;
}

/// A minimal concrete [`Transition`].
///
/// Not a port of `ReactionMonitoringTransition`: it holds the seven trait
/// values and nothing else — no peptide or compound reference, no decoy type,
/// no CV terms, no interpretation list. It exists so that
/// [`MRMTransitionGroup`] can be constructed and tested without the unported
/// assay-library types, and so that callers who only need the group's
/// bookkeeping need not define a type of their own.
///
/// The source class defaults `isDetectingTransition` and
/// `isQuantifyingTransition` to `true` and `isIdentifyingTransition` to
/// `false`; [`Default`] reproduces that. Library intensity defaults to `-1`,
/// as the source `library_intensity_` does.
#[derive(Clone, Debug, PartialEq)]
pub struct SimpleTransition {
    /// Unique identifier; the key under which the group registers it.
    pub native_id: String,
    /// Precursor (Q1) mass-to-charge ratio.
    pub precursor_mz: f64,
    /// Product (Q3) mass-to-charge ratio.
    pub product_mz: f64,
    /// Expected relative intensity from the assay library; `-1` means unset.
    pub library_intensity: f64,
    /// Whether the transition takes part in peak-group detection.
    pub detecting: bool,
    /// Whether the transition is used for peptidoform identification.
    pub identifying: bool,
    /// Whether the transition contributes to quantification.
    pub quantifying: bool,
}

impl Default for SimpleTransition {
    fn default() -> Self {
        Self {
            native_id: String::new(),
            precursor_mz: 0.0,
            product_mz: 0.0,
            library_intensity: -1.0,
            detecting: true,
            identifying: false,
            quantifying: true,
        }
    }
}

impl SimpleTransition {
    /// A transition with the given native ID and library intensity, otherwise
    /// default.
    pub fn new(native_id: impl Into<String>, library_intensity: f64) -> Self {
        Self {
            native_id: native_id.into(),
            library_intensity,
            ..Self::default()
        }
    }
}

impl Transition for SimpleTransition {
    fn native_id(&self) -> &str {
        &self.native_id
    }
    fn precursor_mz(&self) -> f64 {
        self.precursor_mz
    }
    fn product_mz(&self) -> f64 {
        self.product_mz
    }
    fn library_intensity(&self) -> f64 {
        self.library_intensity
    }
    fn is_detecting(&self) -> bool {
        self.detecting
    }
    fn is_identifying(&self) -> bool {
        self.identifying
    }
    fn is_quantifying(&self) -> bool {
        self.quantifying
    }
}

/// The OpenSWATH sub-scores computed for one peak group.
///
/// Ports the data members of `struct OpenSwath_Scores`
/// (`ANALYSIS/OPENSWATH/OpenSwathScores.h`), which `MRMFeature` stores by
/// value. Field names and the per-field defaults are the source's, including
/// the five that start at `-1` rather than `0` to mark "not computed".
///
/// The struct's four LDA scoring methods (`get_quick_lda_score`,
/// `calculate_lda_prescore`, `calculate_lda_single_transition`,
/// `calculate_swath_lda_prescore`) live in
/// `ANALYSIS/OPENSWATH/OpenSwathScores.cpp` and belong to the OpenSWATH
/// analysis port, not to this kernel work package; they are not ported here.
/// `OpenSwath_Scores_Usage` is likewise an analysis-side switchboard and is not
/// ported.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OpenSwathScores {
    pub elution_model_fit_score: f64,
    pub library_corr: f64,
    pub library_norm_manhattan: f64,
    pub library_rootmeansquare: f64,
    pub library_sangle: f64,
    pub norm_rt_score: f64,

    pub isotope_correlation: f64,
    pub isotope_overlap: f64,
    pub massdev_score: f64,
    pub xcorr_coelution_score: f64,
    pub xcorr_shape_score: f64,

    pub yseries_score: f64,
    pub bseries_score: f64,
    pub log_sn_score: f64,

    pub weighted_coelution_score: f64,
    pub weighted_xcorr_shape: f64,
    pub weighted_massdev_score: f64,

    /// Source default `-1`.
    pub ms1_xcorr_coelution_score: f64,
    pub ms1_xcorr_coelution_contrast_score: f64,
    pub ms1_xcorr_coelution_combined_score: f64,
    /// Source default `-1`.
    pub ms1_xcorr_shape_score: f64,
    pub ms1_xcorr_shape_contrast_score: f64,
    pub ms1_xcorr_shape_combined_score: f64,
    pub ms1_ppm_score: f64,
    pub ms1_isotope_correlation: f64,
    pub ms1_isotope_overlap: f64,
    /// Source default `-1`.
    pub ms1_mi_score: f64,
    pub ms1_mi_contrast_score: f64,
    pub ms1_mi_combined_score: f64,

    pub im_xcorr_coelution_score: f64,
    pub im_xcorr_shape_score: f64,
    pub im_delta_score: f64,
    pub im_ms1_delta_score: f64,
    pub im_drift: f64,
    pub im_drift_left: f64,
    pub im_drift_right: f64,
    pub im_drift_weighted: f64,
    /// Source default `-1`.
    pub im_delta: f64,
    pub im_log_intensity: f64,
    pub im_ms1_contrast_coelution: f64,
    pub im_ms1_contrast_shape: f64,
    pub im_ms1_sum_contrast_coelution: f64,
    pub im_ms1_sum_contrast_shape: f64,
    pub im_ms1_drift: f64,
    /// Source default `-1`.
    pub im_ms1_delta: f64,
    pub im_ind_contrast_coelution: f64,
    pub im_ind_contrast_shape: f64,
    pub im_ind_sum_contrast_coelution: f64,
    pub im_ind_sum_contrast_shape: f64,

    pub library_manhattan: f64,
    pub library_dotprod: f64,
    pub intensity: f64,
    pub total_xic: f64,
    pub nr_peaks: f64,
    pub sn_ratio: f64,
    pub mi_score: f64,
    pub weighted_mi_score: f64,

    pub rt_difference: f64,
    pub normalized_experimental_rt: f64,
    pub raw_rt_score: f64,

    pub dotprod_score_dia: f64,
    pub manhatt_score_dia: f64,
}

impl Default for OpenSwathScores {
    fn default() -> Self {
        Self {
            elution_model_fit_score: 0.0,
            library_corr: 0.0,
            library_norm_manhattan: 0.0,
            library_rootmeansquare: 0.0,
            library_sangle: 0.0,
            norm_rt_score: 0.0,
            isotope_correlation: 0.0,
            isotope_overlap: 0.0,
            massdev_score: 0.0,
            xcorr_coelution_score: 0.0,
            xcorr_shape_score: 0.0,
            yseries_score: 0.0,
            bseries_score: 0.0,
            log_sn_score: 0.0,
            weighted_coelution_score: 0.0,
            weighted_xcorr_shape: 0.0,
            weighted_massdev_score: 0.0,
            ms1_xcorr_coelution_score: -1.0,
            ms1_xcorr_coelution_contrast_score: 0.0,
            ms1_xcorr_coelution_combined_score: 0.0,
            ms1_xcorr_shape_score: -1.0,
            ms1_xcorr_shape_contrast_score: 0.0,
            ms1_xcorr_shape_combined_score: 0.0,
            ms1_ppm_score: 0.0,
            ms1_isotope_correlation: 0.0,
            ms1_isotope_overlap: 0.0,
            ms1_mi_score: -1.0,
            ms1_mi_contrast_score: 0.0,
            ms1_mi_combined_score: 0.0,
            im_xcorr_coelution_score: 0.0,
            im_xcorr_shape_score: 0.0,
            im_delta_score: 0.0,
            im_ms1_delta_score: 0.0,
            im_drift: 0.0,
            im_drift_left: 0.0,
            im_drift_right: 0.0,
            im_drift_weighted: 0.0,
            im_delta: -1.0,
            im_log_intensity: 0.0,
            im_ms1_contrast_coelution: 0.0,
            im_ms1_contrast_shape: 0.0,
            im_ms1_sum_contrast_coelution: 0.0,
            im_ms1_sum_contrast_shape: 0.0,
            im_ms1_drift: 0.0,
            im_ms1_delta: -1.0,
            im_ind_contrast_coelution: 0.0,
            im_ind_contrast_shape: 0.0,
            im_ind_sum_contrast_coelution: 0.0,
            im_ind_sum_contrast_shape: 0.0,
            library_manhattan: 0.0,
            library_dotprod: 0.0,
            intensity: 0.0,
            total_xic: 0.0,
            nr_peaks: 0.0,
            sn_ratio: 0.0,
            mi_score: 0.0,
            weighted_mi_score: 0.0,
            rt_difference: 0.0,
            normalized_experimental_rt: 0.0,
            raw_rt_score: 0.0,
            dotprod_score_dia: 0.0,
            manhatt_score_dia: 0.0,
        }
    }
}

/// Per-transition scores for the unique-ion-signature (IPF) workflow.
///
/// Ports the data members of `struct OpenSwath_Ind_Scores`
/// (`ANALYSIS/OPENSWATH/OpenSwathScores.h`). Every field except
/// `ind_num_transitions` and `ind_transition_names` is one value per
/// identifying transition; the source does not require the lists to agree in
/// length and neither does this port, because
/// [`MRMFeature::id_scores_as_meta_value`] writes each list independently.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OpenSwathIndScores {
    pub ind_num_transitions: i32,
    pub ind_transition_names: Vec<String>,
    pub ind_isotope_correlation: Vec<f64>,
    pub ind_isotope_overlap: Vec<f64>,
    pub ind_massdev_score: Vec<f64>,
    pub ind_xcorr_coelution_score: Vec<f64>,
    pub ind_xcorr_shape_score: Vec<f64>,
    pub ind_log_sn_score: Vec<f64>,
    pub ind_area_intensity: Vec<f64>,
    pub ind_total_area_intensity: Vec<f64>,
    pub ind_intensity_score: Vec<f64>,
    pub ind_apex_intensity: Vec<f64>,
    pub ind_apex_position: Vec<f64>,
    pub ind_fwhm: Vec<f64>,
    pub ind_total_mi: Vec<f64>,
    pub ind_log_intensity: Vec<f64>,
    pub ind_intensity_ratio: Vec<f64>,
    pub ind_mi_ratio: Vec<f64>,
    pub ind_mi_score: Vec<f64>,

    pub ind_im_drift: Vec<f64>,
    pub ind_im_drift_left: Vec<f64>,
    pub ind_im_drift_right: Vec<f64>,
    pub ind_im_delta: Vec<f64>,
    pub ind_im_delta_score: Vec<f64>,
    pub ind_im_log_intensity: Vec<f64>,
    pub ind_im_contrast_coelution: Vec<f64>,
    pub ind_im_contrast_shape: Vec<f64>,
    pub ind_im_sum_contrast_coelution: Vec<f64>,
    pub ind_im_sum_contrast_shape: Vec<f64>,

    pub ind_start_position_at_5: Vec<f64>,
    pub ind_end_position_at_5: Vec<f64>,
    pub ind_start_position_at_10: Vec<f64>,
    pub ind_end_position_at_10: Vec<f64>,
    pub ind_start_position_at_50: Vec<f64>,
    pub ind_end_position_at_50: Vec<f64>,
    pub ind_total_width: Vec<f64>,
    pub ind_tailing_factor: Vec<f64>,
    pub ind_asymmetry_factor: Vec<f64>,
    pub ind_slope_of_baseline: Vec<f64>,
    pub ind_baseline_delta_2_height: Vec<f64>,
    pub ind_points_across_baseline: Vec<f64>,
    pub ind_points_across_half_height: Vec<f64>,
}

impl OpenSwathIndScores {
    /// The 42 metadata keys written by [`MRMFeature::id_scores_as_meta_value`],
    /// without the `id_target_` / `id_decoy_` prefix, paired with the field
    /// they come from. Exposed so callers can find and clear the whole block.
    ///
    /// The source writes 43 `setMetaValue` calls but only 42 distinct keys:
    /// `transition_names` is written twice from the same field. See
    /// `OpenMS_CPP_ISSUES.md`.
    pub const KEY_SUFFIXES: [&'static str; 42] = [
        "transition_names",
        "num_transitions",
        "area_intensity",
        "total_area_intensity",
        "intensity_score",
        "intensity_ratio_score",
        "apex_intensity",
        "peak_apex_position",
        "width_at_50",
        "total_mi",
        "ind_log_intensity",
        "ind_xcorr_coelution",
        "ind_xcorr_shape",
        "ind_log_sn_score",
        "ind_isotope_correlation",
        "ind_isotope_overlap",
        "ind_massdev_score",
        "ind_mi_score",
        "ind_mi_ratio_score",
        "ind_im_drift",
        "ind_im_drift_left",
        "ind_im_drift_right",
        "ind_im_delta",
        "ind_im_delta_score",
        "ind_im_log_intensity",
        "ind_im_contrast_coelution",
        "ind_im_contrast_shape",
        "ind_im_sum_contrast_coelution",
        "ind_im_sum_contrast_shape",
        "ind_start_position_at_5",
        "ind_end_position_at_5",
        "ind_start_position_at_10",
        "ind_end_position_at_10",
        "ind_start_position_at_50",
        "ind_end_position_at_50",
        "ind_total_width",
        "ind_tailing_factor",
        "ind_asymmetry_factor",
        "ind_slope_of_baseline",
        "ind_baseline_delta_2_height",
        "ind_points_across_baseline",
        "ind_points_across_half_height",
    ];

    /// Total number of stored values across every list, for preflighting.
    fn value_count(&self) -> Option<usize> {
        let lists = self.float_lists();
        let mut total = self.ind_transition_names.len();
        for list in lists {
            total = total.checked_add(list.len())?;
        }
        Some(total)
    }

    /// Every `Vec<f64>` member in the order `IDScoresAsMetaValue` writes them.
    fn float_lists(&self) -> [&Vec<f64>; 40] {
        [
            &self.ind_area_intensity,
            &self.ind_total_area_intensity,
            &self.ind_intensity_score,
            &self.ind_intensity_ratio,
            &self.ind_apex_intensity,
            &self.ind_apex_position,
            &self.ind_fwhm,
            &self.ind_total_mi,
            &self.ind_log_intensity,
            &self.ind_xcorr_coelution_score,
            &self.ind_xcorr_shape_score,
            &self.ind_log_sn_score,
            &self.ind_isotope_correlation,
            &self.ind_isotope_overlap,
            &self.ind_massdev_score,
            &self.ind_mi_score,
            &self.ind_mi_ratio,
            &self.ind_im_drift,
            &self.ind_im_drift_left,
            &self.ind_im_drift_right,
            &self.ind_im_delta,
            &self.ind_im_delta_score,
            &self.ind_im_log_intensity,
            &self.ind_im_contrast_coelution,
            &self.ind_im_contrast_shape,
            &self.ind_im_sum_contrast_coelution,
            &self.ind_im_sum_contrast_shape,
            &self.ind_start_position_at_5,
            &self.ind_end_position_at_5,
            &self.ind_start_position_at_10,
            &self.ind_end_position_at_10,
            &self.ind_start_position_at_50,
            &self.ind_end_position_at_50,
            &self.ind_total_width,
            &self.ind_tailing_factor,
            &self.ind_asymmetry_factor,
            &self.ind_slope_of_baseline,
            &self.ind_baseline_delta_2_height,
            &self.ind_points_across_baseline,
            &self.ind_points_across_half_height,
        ]
    }
}

/// A peak group: one elution peak seen across several transition chromatograms.
///
/// Ports `MRMFeature`. The source derives from `Feature`; this port holds the
/// feature in the public [`MRMFeature::feature`] field and dereferences to it,
/// the same shape `Feature` itself uses for `BaseFeature`, so `mrm.rt`,
/// `mrm.intensity` and `mrm.metadata` all reach through.
///
/// The individual features found in each chromatogram are ordinary
/// [`Feature`] values addressed by the key they were
/// added under — in practice the transition's native ID. Precursor features,
/// extracted from the MS1 map, form a second, independent list with its own
/// keys. The peak-group score record is separate from both.
///
/// Both key maps and both lists are private because the source's are: the
/// index a key maps to must stay inside the list it indexes.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MRMFeature {
    /// The feature this peak group is, in the base kernel type. Source
    /// `MRMFeature` inherits from `Feature`; all of `Feature`'s API is reached
    /// through this field or through `Deref`.
    pub feature: Feature,
    features: Vec<Feature>,
    precursor_features: Vec<Feature>,
    scores: OpenSwathScores,
    feature_map: BTreeMap<String, usize>,
    precursor_feature_map: BTreeMap<String, usize>,
}

impl Deref for MRMFeature {
    type Target = Feature;
    fn deref(&self) -> &Self::Target {
        &self.feature
    }
}

impl DerefMut for MRMFeature {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.feature
    }
}

impl From<Feature> for MRMFeature {
    fn from(feature: Feature) -> Self {
        Self {
            feature,
            ..Self::default()
        }
    }
}

fn insert_keyed(
    list: &mut Vec<Feature>,
    map: &mut BTreeMap<String, usize>,
    feature: Feature,
    key: &str,
    policy: DuplicateKeyPolicy,
    what: &str,
) -> Result<()> {
    if list.len() >= MRMFeature::MAX_FEATURES {
        return Err(limit("feature count"));
    }
    let occupied = map.contains_key(key);
    if occupied && policy == DuplicateKeyPolicy::Reject {
        return Err(Error::InvalidValue(format!(
            "{what} with key '{key}' was already present"
        )));
    }
    list.push(feature);
    map.insert(key.to_owned(), list.len() - 1);
    Ok(())
}

impl MRMFeature {
    /// Maximum number of sub-features, checked before each insertion. Applies
    /// independently to the transition features and to the precursor features.
    pub const MAX_FEATURES: usize = 1_000_000;

    /// Maximum number of individual score values accepted by
    /// [`MRMFeature::id_scores_as_meta_value`] in one call.
    pub const MAX_SCORE_VALUES: usize = 4_000_000;

    /// Conservative cumulative owned payload accepted by
    /// [`MRMFeature::id_scores_as_meta_value`] in one call.
    pub const MAX_BYTES: usize = 64 * 1024 * 1024;

    /// An empty peak group, as the source default constructor.
    pub fn new() -> Self {
        Self::default()
    }

    /// The peak-group score record (source `getScores`).
    pub fn scores(&self) -> &OpenSwathScores {
        &self.scores
    }

    /// Mutable peak-group score record (source's non-const `getScores`).
    pub fn scores_mut(&mut self) -> &mut OpenSwathScores {
        &mut self.scores
    }

    /// Replace the peak-group score record (source `setScores`).
    ///
    /// The source copy constructor and assignment operator call this in
    /// addition to copying the member directly, which is redundant but
    /// harmless; the derived [`Clone`] here copies it once.
    pub fn set_scores(&mut self, scores: OpenSwathScores) {
        self.scores = scores;
    }

    /// Store a single named score.
    ///
    /// Source `addScore` is `setMetaValue(score_name, score)`: despite the
    /// name, scores added this way go into the feature's metadata, not into
    /// the [`OpenSwathScores`] record, and a repeated name overwrites.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `score` is not finite. The source
    /// stores a non-finite `DataValue` without complaint; this port refuses,
    /// because [`MetaValue`] guarantees finite floating members.
    pub fn add_score(&mut self, score_name: impl Into<String>, score: f64) -> Result<()> {
        let value = MetaValue::try_from(score)?;
        self.feature.base.metadata.insert(score_name.into(), value);
        Ok(())
    }

    /// Read back a score stored by [`MRMFeature::add_score`].
    ///
    /// Native: the source class once had `double getScore(const std::string&)`
    /// — the class test still carries a `NOT_TESTABLE` section for it — but the
    /// member is gone from the pinned header and callers reach the value
    /// through `getMetaValue`. `None` means no such key; a key holding a
    /// non-numeric value also yields `None`.
    pub fn score(&self, score_name: &str) -> Option<f64> {
        self.feature
            .base
            .metadata
            .get(score_name)
            .and_then(|value| value.as_f64().ok())
    }

    /// Add a feature found in one transition chromatogram, under `key`.
    ///
    /// Source `addFeature`. `key` is the identifier the feature is retrieved
    /// by; OpenSWATH passes the transition's native ID.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `key` is already in use, and when
    /// the list already holds [`MRMFeature::MAX_FEATURES`] features. The source
    /// accepts a repeated key, appends anyway and re-points the key, which
    /// leaves the previously keyed feature in the list and reachable only by
    /// iterating it; use
    /// [`add_feature_with`](MRMFeature::add_feature_with) with
    /// [`DuplicateKeyPolicy::SourceOverwrite`] for that behaviour.
    pub fn add_feature(&mut self, feature: Feature, key: impl AsRef<str>) -> Result<()> {
        self.add_feature_with(feature, key, DuplicateKeyPolicy::Reject)
    }

    /// [`add_feature`](MRMFeature::add_feature) with an explicit duplicate-key
    /// policy. On any error nothing is stored.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] on a duplicate key under
    /// [`DuplicateKeyPolicy::Reject`], or when the list is already at
    /// [`MRMFeature::MAX_FEATURES`].
    pub fn add_feature_with(
        &mut self,
        feature: Feature,
        key: impl AsRef<str>,
        policy: DuplicateKeyPolicy,
    ) -> Result<()> {
        insert_keyed(
            &mut self.features,
            &mut self.feature_map,
            feature,
            key.as_ref(),
            policy,
            "feature",
        )
    }

    /// The feature registered under `key` (source `getFeature`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when `key` is unknown. The
    /// source's const overload throws `std::out_of_range`; its non-const
    /// overload uses `operator[]` on the map, which inserts `key -> 0` and then
    /// returns the *first* feature, or throws `std::out_of_range` when the list
    /// is empty. Both are recorded in `OpenMS_CPP_ISSUES.md`; neither is
    /// reproduced.
    pub fn feature(&self, key: &str) -> Result<&Feature> {
        let index = *self
            .feature_map
            .get(key)
            .ok_or_else(|| missing("feature", key))?;
        self.features
            .get(index)
            .ok_or_else(|| missing("feature", key))
    }

    /// Mutable access to the feature registered under `key`.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when `key` is unknown.
    pub fn feature_mut(&mut self, key: &str) -> Result<&mut Feature> {
        let index = *self
            .feature_map
            .get(key)
            .ok_or_else(|| missing("feature", key))?;
        self.features
            .get_mut(index)
            .ok_or_else(|| missing("feature", key))
    }

    /// Whether a feature is registered under `key`. Native; the source offers
    /// no such test and callers reach for `getFeature` directly.
    pub fn has_feature(&self, key: &str) -> bool {
        self.feature_map.contains_key(key)
    }

    /// All transition features, in insertion order (source `getFeatures`).
    pub fn features(&self) -> &[Feature] {
        &self.features
    }

    /// Keys of the registered transition features, in lexicographic order.
    ///
    /// Source `getFeatureIDs(std::vector<std::string>& result)` *appends* to the
    /// caller's vector rather than clearing it, and walks a `std::map`, so the
    /// order is lexicographic by key, not insertion order. This returns a
    /// borrowed iterator instead of an out-parameter; a caller who wants the
    /// source's append can `extend` their own vector from it.
    pub fn feature_ids(&self) -> impl ExactSizeIterator<Item = &str> {
        self.feature_map.keys().map(String::as_str)
    }

    /// Add a precursor feature extracted from the MS1 map, under `key`.
    ///
    /// Source `addPrecursorFeature`. The precursor list and its keys are
    /// entirely independent of the transition-feature list: the same key may
    /// appear in both.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] on a duplicate key or at
    /// [`MRMFeature::MAX_FEATURES`]; see
    /// [`add_feature`](MRMFeature::add_feature) for the source's differing
    /// duplicate handling.
    pub fn add_precursor_feature(&mut self, feature: Feature, key: impl AsRef<str>) -> Result<()> {
        self.add_precursor_feature_with(feature, key, DuplicateKeyPolicy::Reject)
    }

    /// [`add_precursor_feature`](MRMFeature::add_precursor_feature) with an
    /// explicit duplicate-key policy. On any error nothing is stored.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] on a duplicate key under
    /// [`DuplicateKeyPolicy::Reject`], or at [`MRMFeature::MAX_FEATURES`].
    pub fn add_precursor_feature_with(
        &mut self,
        feature: Feature,
        key: impl AsRef<str>,
        policy: DuplicateKeyPolicy,
    ) -> Result<()> {
        insert_keyed(
            &mut self.precursor_features,
            &mut self.precursor_feature_map,
            feature,
            key.as_ref(),
            policy,
            "precursor feature",
        )
    }

    /// The precursor feature registered under `key`
    /// (source `getPrecursorFeature`).
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when `key` is unknown; see
    /// [`feature`](MRMFeature::feature) for what the source does instead.
    pub fn precursor_feature(&self, key: &str) -> Result<&Feature> {
        let index = *self
            .precursor_feature_map
            .get(key)
            .ok_or_else(|| missing("precursor feature", key))?;
        self.precursor_features
            .get(index)
            .ok_or_else(|| missing("precursor feature", key))
    }

    /// Mutable access to the precursor feature registered under `key`.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when `key` is unknown.
    pub fn precursor_feature_mut(&mut self, key: &str) -> Result<&mut Feature> {
        let index = *self
            .precursor_feature_map
            .get(key)
            .ok_or_else(|| missing("precursor feature", key))?;
        self.precursor_features
            .get_mut(index)
            .ok_or_else(|| missing("precursor feature", key))
    }

    /// Whether a precursor feature is registered under `key`. Native, as
    /// [`has_feature`](MRMFeature::has_feature).
    pub fn has_precursor_feature(&self, key: &str) -> bool {
        self.precursor_feature_map.contains_key(key)
    }

    /// All precursor features, in insertion order.
    ///
    /// Native: the source exposes the precursor list only through
    /// `getPrecursorFeatureIDs` and `getPrecursorFeature`, with no counterpart
    /// to `getFeatures`.
    pub fn precursor_features(&self) -> &[Feature] {
        &self.precursor_features
    }

    /// Keys of the registered precursor features, in lexicographic order.
    ///
    /// Source `getPrecursorFeatureIDs`; see
    /// [`feature_ids`](MRMFeature::feature_ids) for the out-parameter and
    /// ordering notes, which apply identically.
    pub fn precursor_feature_ids(&self) -> impl ExactSizeIterator<Item = &str> {
        self.precursor_feature_map.keys().map(String::as_str)
    }

    /// Write the IPF per-transition scores into the feature's metadata.
    ///
    /// Source `IDScoresAsMetaValue`. Keys are prefixed `id_decoy_` when `decoy`
    /// is true and `id_target_` otherwise, and the 42 suffixes are
    /// [`OpenSwathIndScores::KEY_SUFFIXES`]. Existing keys are overwritten;
    /// other metadata is untouched.
    ///
    /// Every value is written even when its list is empty, matching the source,
    /// so a reader can rely on the whole block being present.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when any score is not finite — the
    /// source stores non-finite `DataValue`s silently — and when the scores
    /// hold more than [`MRMFeature::MAX_SCORE_VALUES`] values or would exceed
    /// [`MRMFeature::MAX_BYTES`]. The check runs before anything is written and
    /// the block is committed in one step, so a rejected call leaves the
    /// feature's metadata exactly as it was.
    pub fn id_scores_as_meta_value(
        &mut self,
        decoy: bool,
        idscores: &OpenSwathIndScores,
    ) -> Result<()> {
        let prefix = if decoy { "id_decoy_" } else { "id_target_" };
        let values = idscores.value_count().ok_or_else(|| limit("score count"))?;
        if values > Self::MAX_SCORE_VALUES {
            return Err(limit("score count"));
        }
        let names: usize = idscores
            .ind_transition_names
            .iter()
            .try_fold(0usize, |total, name| total.checked_add(name.len()))
            .ok_or_else(|| limit("score payload"))?;
        let bytes = values
            .checked_mul(std::mem::size_of::<f64>() + std::mem::size_of::<MetaValue>())
            .and_then(|n| n.checked_add(names))
            .ok_or_else(|| limit("score payload"))?;
        if bytes > Self::MAX_BYTES {
            return Err(limit("score payload"));
        }

        let mut block: Vec<(String, MetaValue)> = Vec::with_capacity(42);
        let mut push_floats = |suffix: &str, list: &[f64]| -> Result<()> {
            let value = MetaValue::new(MetaValueData::FloatList(list.to_vec()))?;
            block.push((format!("{prefix}{suffix}"), value));
            Ok(())
        };

        // Source order, with the duplicated `transition_names` write folded to
        // one. The three non-float entries are handled separately below.
        push_floats("area_intensity", &idscores.ind_area_intensity)?;
        push_floats("total_area_intensity", &idscores.ind_total_area_intensity)?;
        push_floats("intensity_score", &idscores.ind_intensity_score)?;
        push_floats("intensity_ratio_score", &idscores.ind_intensity_ratio)?;
        push_floats("apex_intensity", &idscores.ind_apex_intensity)?;
        push_floats("peak_apex_position", &idscores.ind_apex_position)?;
        push_floats("width_at_50", &idscores.ind_fwhm)?;
        push_floats("total_mi", &idscores.ind_total_mi)?;
        push_floats("ind_log_intensity", &idscores.ind_log_intensity)?;
        push_floats("ind_xcorr_coelution", &idscores.ind_xcorr_coelution_score)?;
        push_floats("ind_xcorr_shape", &idscores.ind_xcorr_shape_score)?;
        push_floats("ind_log_sn_score", &idscores.ind_log_sn_score)?;
        push_floats("ind_isotope_correlation", &idscores.ind_isotope_correlation)?;
        push_floats("ind_isotope_overlap", &idscores.ind_isotope_overlap)?;
        push_floats("ind_massdev_score", &idscores.ind_massdev_score)?;
        push_floats("ind_mi_score", &idscores.ind_mi_score)?;
        push_floats("ind_mi_ratio_score", &idscores.ind_mi_ratio)?;
        push_floats("ind_im_drift", &idscores.ind_im_drift)?;
        push_floats("ind_im_drift_left", &idscores.ind_im_drift_left)?;
        push_floats("ind_im_drift_right", &idscores.ind_im_drift_right)?;
        push_floats("ind_im_delta", &idscores.ind_im_delta)?;
        push_floats("ind_im_delta_score", &idscores.ind_im_delta_score)?;
        push_floats("ind_im_log_intensity", &idscores.ind_im_log_intensity)?;
        push_floats(
            "ind_im_contrast_coelution",
            &idscores.ind_im_contrast_coelution,
        )?;
        push_floats("ind_im_contrast_shape", &idscores.ind_im_contrast_shape)?;
        push_floats(
            "ind_im_sum_contrast_coelution",
            &idscores.ind_im_sum_contrast_coelution,
        )?;
        push_floats(
            "ind_im_sum_contrast_shape",
            &idscores.ind_im_sum_contrast_shape,
        )?;
        push_floats("ind_start_position_at_5", &idscores.ind_start_position_at_5)?;
        push_floats("ind_end_position_at_5", &idscores.ind_end_position_at_5)?;
        push_floats(
            "ind_start_position_at_10",
            &idscores.ind_start_position_at_10,
        )?;
        push_floats("ind_end_position_at_10", &idscores.ind_end_position_at_10)?;
        push_floats(
            "ind_start_position_at_50",
            &idscores.ind_start_position_at_50,
        )?;
        push_floats("ind_end_position_at_50", &idscores.ind_end_position_at_50)?;
        push_floats("ind_total_width", &idscores.ind_total_width)?;
        push_floats("ind_tailing_factor", &idscores.ind_tailing_factor)?;
        push_floats("ind_asymmetry_factor", &idscores.ind_asymmetry_factor)?;
        push_floats("ind_slope_of_baseline", &idscores.ind_slope_of_baseline)?;
        push_floats(
            "ind_baseline_delta_2_height",
            &idscores.ind_baseline_delta_2_height,
        )?;
        push_floats(
            "ind_points_across_baseline",
            &idscores.ind_points_across_baseline,
        )?;
        push_floats(
            "ind_points_across_half_height",
            &idscores.ind_points_across_half_height,
        )?;

        block.push((
            format!("{prefix}transition_names"),
            MetaValue::from(idscores.ind_transition_names.clone()),
        ));
        block.push((
            format!("{prefix}num_transitions"),
            MetaValue::from(i64::from(idscores.ind_num_transitions)),
        ));

        self.feature.base.metadata.extend(block);
        Ok(())
    }
}

/// A group of transitions, their chromatograms and the peak groups on them.
///
/// Ports the header-only template `MRMTransitionGroup<ChromatogramType,
/// TransitionType>`. `C` is the raw-data type — the source header notes that
/// not every OpenMS function accepts `MSChromatogram`, so `MSSpectrum` must
/// work too — and `T` is the assay transition, generic because the source's
/// `ReactionMonitoringTransition` is not ported.
///
/// The group stores three keyed lists: the transitions, the fragment-ion
/// chromatograms measured for them, and the precursor chromatograms extracted
/// from the MS1 map. The source header states the invariant plainly: *for the
/// data structure to be consistent it needs the same identifiers for the
/// chromatograms as for the transitions.* Nothing in the source enforces it at
/// run time in a release build; see
/// [`is_internally_consistent`](MRMTransitionGroup::is_internally_consistent).
///
/// ```
/// use openms::kernel::MSChromatogram;
/// use openms::kernel::mrm::{MRMTransitionGroup, SimpleTransition};
///
/// let mut group: MRMTransitionGroup<MSChromatogram, SimpleTransition> =
///     MRMTransitionGroup::default();
/// group.set_transition_group_id("PEPTIDEK_2");
///
/// let transition = SimpleTransition::new("y7", 12_000.0);
/// group.add_transition(transition, "y7")?;
///
/// let chromatogram = MSChromatogram {
///     native_id: "y7".into(),
///     ..MSChromatogram::default()
/// };
/// group.add_chromatogram(chromatogram, "y7")?;
///
/// assert_eq!(group.size(), 1);
/// assert!(group.is_internally_consistent());
/// assert!(group.chromatogram_ids_match());
/// assert_eq!(group.library_intensity()?, vec![12_000.0]);
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct MRMTransitionGroup<C, T> {
    transition_group_id: String,
    transitions: Vec<T>,
    chromatograms: Vec<C>,
    precursor_chromatograms: Vec<C>,
    mrm_features: Vec<MRMFeature>,
    transition_map: BTreeMap<String, usize>,
    chromatogram_map: BTreeMap<String, usize>,
    precursor_chromatogram_map: BTreeMap<String, usize>,
}

impl<C, T> Default for MRMTransitionGroup<C, T> {
    fn default() -> Self {
        Self {
            transition_group_id: String::new(),
            transitions: Vec::new(),
            chromatograms: Vec::new(),
            precursor_chromatograms: Vec::new(),
            mrm_features: Vec::new(),
            transition_map: BTreeMap::new(),
            chromatogram_map: BTreeMap::new(),
            precursor_chromatogram_map: BTreeMap::new(),
        }
    }
}

fn insert_chromatogram<C>(
    list: &mut Vec<C>,
    map: &mut BTreeMap<String, usize>,
    value: C,
    key: &str,
    what: &str,
) -> Result<()> {
    if list.len() >= MAX_ITEMS {
        return Err(limit("list length"));
    }
    if map.contains_key(key) {
        return Err(Error::InvalidValue(format!(
            "Internal error: {what} with nativeID was already present: '{key}'"
        )));
    }
    list.push(value);
    map.insert(key.to_owned(), list.len() - 1);
    Ok(())
}

/// Shared ceiling on every list a transition group holds.
const MAX_ITEMS: usize = 1_000_000;

impl<C, T> MRMTransitionGroup<C, T> {
    /// Maximum entries in any one of the group's four lists, and maximum number
    /// of transition identifiers accepted by
    /// [`subset`](MRMTransitionGroup::subset) and
    /// [`subset_dependent`](MRMTransitionGroup::subset_dependent).
    pub const MAX_ITEMS: usize = MAX_ITEMS;

    /// Conservative cumulative owned payload a subset operation may build.
    pub const MAX_BYTES: usize = 256 * 1024 * 1024;

    /// An empty transition group, as the source default constructor.
    pub fn new() -> Self {
        Self::default()
    }

    /// The number of fragment-ion chromatograms.
    ///
    /// Source `size()` counts the chromatograms only: the transitions, the
    /// precursor chromatograms and the features are not part of it. When the
    /// group is internally consistent this equals the transition count.
    pub fn size(&self) -> usize {
        self.chromatograms.len()
    }

    /// Whether the group holds no fragment-ion chromatogram. Native, the
    /// companion predicate to [`size`](MRMTransitionGroup::size).
    pub fn is_empty(&self) -> bool {
        self.chromatograms.is_empty()
    }

    /// The peak-group identifier (source `getTransitionGroupID`).
    pub fn transition_group_id(&self) -> &str {
        &self.transition_group_id
    }

    /// Set the peak-group identifier (source `setTransitionGroupID`).
    pub fn set_transition_group_id(&mut self, id: impl Into<String>) {
        self.transition_group_id = id.into();
    }

    /// All transitions, in insertion order (source `getTransitions`).
    pub fn transitions(&self) -> &[T] {
        &self.transitions
    }

    /// Mutable transitions, in insertion order.
    ///
    /// Source `getTransitionsMuteable` hands out the `std::vector` itself, so a
    /// caller can push, erase or sort it and silently invalidate every stored
    /// index. This returns a slice: fields of a transition may be edited, the
    /// length cannot change. Reordering the slice or editing a native ID still
    /// breaks the mapping;
    /// [`is_internally_consistent`](MRMTransitionGroup::is_internally_consistent)
    /// and [`chromatogram_ids_match`](MRMTransitionGroup::chromatogram_ids_match)
    /// detect that afterwards.
    pub fn transitions_mut(&mut self) -> &mut [T] {
        &mut self.transitions
    }

    /// Register a transition under `key`.
    ///
    /// The source documents that transitions are mapped by their native ID,
    /// i.e. `TransitionType::getNativeID`, and that the same key must be used
    /// when querying. The key is not checked against the transition's own
    /// native ID here — the source does not check either — but a mismatch makes
    /// [`subset`](MRMTransitionGroup::subset) skip the transition, because that
    /// method looks itself up by native ID.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `key` is already in use, matching
    /// the source's `Exception::InvalidValue` ("Internal error: Transition with
    /// nativeID was already present!"), and when the list already holds
    /// [`MRMTransitionGroup::MAX_ITEMS`] transitions. Nothing is stored on
    /// error.
    pub fn add_transition(&mut self, transition: T, key: impl AsRef<str>) -> Result<()> {
        let key = key.as_ref();
        if self.transitions.len() >= MAX_ITEMS {
            return Err(limit("list length"));
        }
        if self.transition_map.contains_key(key) {
            return Err(Error::InvalidValue(format!(
                "Internal error: Transition with nativeID was already present: '{key}'"
            )));
        }
        self.transitions.push(transition);
        self.transition_map
            .insert(key.to_owned(), self.transitions.len() - 1);
        Ok(())
    }

    /// Whether a transition is registered under `key` (source `hasTransition`).
    pub fn has_transition(&self, key: &str) -> bool {
        self.transition_map.contains_key(key)
    }

    /// The transition registered under `key` (source `getTransition`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when `key` is unknown. The source
    /// guards both the key and the index with `OPENMS_PRECONDITION`, which is
    /// compiled out unless `OPENMS_ASSERTIONS` is set, so a release build
    /// indexes the vector with whatever `transition_map_[key]` default-inserts
    /// — zero — and returns the first transition or reads out of bounds on an
    /// empty group.
    pub fn transition(&self, key: &str) -> Result<&T> {
        let index = *self
            .transition_map
            .get(key)
            .ok_or_else(|| missing("transition", key))?;
        self.transitions
            .get(index)
            .ok_or_else(|| missing("transition", key))
    }

    /// All fragment-ion chromatograms, in insertion order
    /// (source `getChromatograms`, const overload).
    pub fn chromatograms(&self) -> &[C] {
        &self.chromatograms
    }

    /// Mutable fragment-ion chromatograms, in insertion order.
    ///
    /// The source's non-const `getChromatograms` returns the `std::vector`
    /// itself; this returns a slice, for the reason given on
    /// [`transitions_mut`](MRMTransitionGroup::transitions_mut).
    pub fn chromatograms_mut(&mut self) -> &mut [C] {
        &mut self.chromatograms
    }

    /// Register a fragment-ion chromatogram under `key`.
    ///
    /// The source documents `ChromatogramType::getNativeID` as a good choice of
    /// key and requires the same key when querying;
    /// [`chromatogram_ids_match`](MRMTransitionGroup::chromatogram_ids_match)
    /// is the check that the two actually agree.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `key` is already in use (source
    /// `Exception::InvalidValue`, "Internal error: Chromatogram with nativeID
    /// was already present!") or at [`MRMTransitionGroup::MAX_ITEMS`]. Nothing
    /// is stored on error.
    pub fn add_chromatogram(&mut self, chromatogram: C, key: impl AsRef<str>) -> Result<()> {
        insert_chromatogram(
            &mut self.chromatograms,
            &mut self.chromatogram_map,
            chromatogram,
            key.as_ref(),
            "Chromatogram",
        )
    }

    /// Whether a fragment-ion chromatogram is registered under `key`
    /// (source `hasChromatogram`).
    pub fn has_chromatogram(&self, key: &str) -> bool {
        self.chromatogram_map.contains_key(key)
    }

    /// The fragment-ion chromatogram registered under `key`
    /// (source `getChromatogram`).
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when `key` is unknown; the source's
    /// preconditions are compiled out in release builds, as described on
    /// [`transition`](MRMTransitionGroup::transition).
    pub fn chromatogram(&self, key: &str) -> Result<&C> {
        let index = *self
            .chromatogram_map
            .get(key)
            .ok_or_else(|| missing("chromatogram", key))?;
        self.chromatograms
            .get(index)
            .ok_or_else(|| missing("chromatogram", key))
    }

    /// Mutable access to the fragment-ion chromatogram registered under `key`
    /// (source's non-const `getChromatogram`).
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when `key` is unknown.
    pub fn chromatogram_mut(&mut self, key: &str) -> Result<&mut C> {
        let index = *self
            .chromatogram_map
            .get(key)
            .ok_or_else(|| missing("chromatogram", key))?;
        self.chromatograms
            .get_mut(index)
            .ok_or_else(|| missing("chromatogram", key))
    }

    /// All precursor chromatograms, in insertion order
    /// (source `getPrecursorChromatograms`, const overload).
    pub fn precursor_chromatograms(&self) -> &[C] {
        &self.precursor_chromatograms
    }

    /// Mutable precursor chromatograms, in insertion order; a slice, as
    /// [`chromatograms_mut`](MRMTransitionGroup::chromatograms_mut).
    pub fn precursor_chromatograms_mut(&mut self) -> &mut [C] {
        &mut self.precursor_chromatograms
    }

    /// Register a precursor chromatogram, extracted from an MS1 map, under
    /// `key`.
    ///
    /// # Arguments
    ///
    /// * `chromatogram` — chromatographic traces from the MS1 map to be added.
    /// * `key` — unique identifier of the chromatogram, e.g. its native ID.
    ///
    /// The precursor list has its own key map: a key used for a fragment-ion
    /// chromatogram may be reused here.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `key` is already used by another precursor
    /// chromatogram, or at [`MRMTransitionGroup::MAX_ITEMS`]. Nothing is stored
    /// on error.
    pub fn add_precursor_chromatogram(
        &mut self,
        chromatogram: C,
        key: impl AsRef<str>,
    ) -> Result<()> {
        insert_chromatogram(
            &mut self.precursor_chromatograms,
            &mut self.precursor_chromatogram_map,
            chromatogram,
            key.as_ref(),
            "Chromatogram",
        )
    }

    /// Whether a precursor chromatogram is registered under `key`
    /// (source `hasPrecursorChromatogram`).
    pub fn has_precursor_chromatogram(&self, key: &str) -> bool {
        self.precursor_chromatogram_map.contains_key(key)
    }

    /// The precursor chromatogram registered under `key`
    /// (source `getPrecursorChromatogram`).
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when `key` is unknown; the source's
    /// preconditions are compiled out in release builds.
    pub fn precursor_chromatogram(&self, key: &str) -> Result<&C> {
        let index = *self
            .precursor_chromatogram_map
            .get(key)
            .ok_or_else(|| missing("precursor chromatogram", key))?;
        self.precursor_chromatograms
            .get(index)
            .ok_or_else(|| missing("precursor chromatogram", key))
    }

    /// Mutable access to the precursor chromatogram registered under `key`.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when `key` is unknown.
    pub fn precursor_chromatogram_mut(&mut self, key: &str) -> Result<&mut C> {
        let index = *self
            .precursor_chromatogram_map
            .get(key)
            .ok_or_else(|| missing("precursor chromatogram", key))?;
        self.precursor_chromatograms
            .get_mut(index)
            .ok_or_else(|| missing("precursor chromatogram", key))
    }

    /// The peak groups found across this group's chromatograms
    /// (source `getFeatures`).
    pub fn features(&self) -> &[MRMFeature] {
        &self.mrm_features
    }

    /// Mutable peak groups (source `getFeaturesMuteable`).
    ///
    /// Unlike the transition and chromatogram lists, the feature list has no
    /// key map, so handing out the vector itself cannot desynchronise anything
    /// and the source's full mutability is preserved.
    pub fn features_mut(&mut self) -> &mut Vec<MRMFeature> {
        &mut self.mrm_features
    }

    /// Append a peak group (source `addFeature`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the list already holds
    /// [`MRMTransitionGroup::MAX_ITEMS`] features. The source has no ceiling.
    pub fn add_feature(&mut self, feature: MRMFeature) -> Result<()> {
        if self.mrm_features.len() >= MAX_ITEMS {
            return Err(limit("list length"));
        }
        self.mrm_features.push(feature);
        Ok(())
    }

    /// Whether the transition, chromatogram and key-map counts agree and every
    /// chromatogram key also names a transition.
    ///
    /// Source `isInternallyConsistent` states the three conditions as
    /// `OPENMS_PRECONDITION`s and then returns `true`. Those macros expand to
    /// nothing unless `OPENMS_ASSERTIONS` is defined, so in an ordinary release
    /// build the function is `return true;` and cannot report an inconsistent
    /// group; in an assertions build a violation throws
    /// `Exception::Precondition` rather than returning `false`. This port
    /// evaluates the conditions and returns the answer, which is what the name
    /// and return type promise. The cost is one pass over the chromatogram key
    /// map. See `OpenMS_CPP_ISSUES.md`.
    pub fn is_internally_consistent(&self) -> bool {
        self.transitions.len() == self.chromatograms.len()
            && self.transition_map.len() == self.chromatogram_map.len()
            && self.is_mapping_consistent()
    }

    /// Source `isMappingConsistent_`: equal map sizes, and every chromatogram
    /// key is also a transition key.
    fn is_mapping_consistent(&self) -> bool {
        self.transition_map.len() == self.chromatogram_map.len()
            && self
                .chromatogram_map
                .keys()
                .all(|key| self.transition_map.contains_key(key))
    }
}

impl<C: NativeIdentified, T> MRMTransitionGroup<C, T> {
    /// Whether every chromatogram's own native ID equals the key it is stored
    /// under, for both the fragment-ion and the precursor list.
    ///
    /// Source `chromatogramIdsMatch`. The fragment-ion map is walked first and
    /// the precursor map second; the first mismatch short-circuits. Storing the
    /// same chromatogram twice under two keys therefore fails, because at most
    /// one of the keys can equal its single native ID.
    pub fn chromatogram_ids_match(&self) -> bool {
        let fragment = self.chromatogram_map.iter().all(|(key, &index)| {
            self.chromatograms
                .get(index)
                .is_some_and(|chromatogram| chromatogram.native_id() == key)
        });
        if !fragment {
            return false;
        }
        self.precursor_chromatogram_map.iter().all(|(key, &index)| {
            self.precursor_chromatograms
                .get(index)
                .is_some_and(|chromatogram| chromatogram.native_id() == key)
        })
    }
}

impl<C, T: Transition> MRMTransitionGroup<C, T> {
    /// The library intensities of the transitions, in insertion order, with
    /// negatives replaced by zero.
    ///
    /// Source `getLibraryIntensity(std::vector<double>& result)` *appends* to
    /// the caller's vector and then clamps the whole vector, so pre-existing
    /// negative entries belonging to the caller are clamped as well. This
    /// returns a fresh vector and clamps only its own entries. The comment in
    /// the source explains the clamp: the library intensity should never be
    /// below zero.
    ///
    /// A non-finite library intensity is passed through unchanged, as in the
    /// source; `NaN < 0.0` is false, so a NaN is not clamped.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the group holds more than
    /// [`MRMTransitionGroup::MAX_ITEMS`] transitions, which can only happen if
    /// the list was filled by something other than
    /// [`add_transition`](MRMTransitionGroup::add_transition).
    pub fn library_intensity(&self) -> Result<Vec<f64>> {
        if self.transitions.len() > MAX_ITEMS {
            return Err(limit("list length"));
        }
        Ok(self
            .transitions
            .iter()
            .map(|transition| {
                let intensity = transition.library_intensity();
                if intensity < 0.0 { 0.0 } else { intensity }
            })
            .collect())
    }
}

impl<C: Clone, T: Transition + Clone> MRMTransitionGroup<C, T> {
    /// A new group holding only the named transitions, their chromatograms and
    /// the whole feature list restricted to those transitions.
    ///
    /// Source `subsetDependent`. It differs from
    /// [`subset`](MRMTransitionGroup::subset) in two ways, both of which the
    /// name refers to: the peak groups are copied *whole*, keeping the
    /// per-transition features of transitions that are not in `tr_ids`, and no
    /// precursor chromatogram is carried over. The transition group ID is
    /// copied.
    ///
    /// Transitions are selected by their own native ID, not by the key they
    /// were registered under, and each selected transition must have a
    /// chromatogram registered under that native ID — the source indexes
    /// `chromatogram_map_.at()` with no guard, unlike `subset`.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when a selected transition has no
    /// chromatogram under its native ID (the source throws
    /// `std::out_of_range`), [`Error::InvalidValue`] when two selected
    /// transitions share a native ID, and [`Error::InvalidValue`] when the
    /// input exceeds [`MRMTransitionGroup::MAX_ITEMS`] or
    /// [`MRMTransitionGroup::MAX_BYTES`]. The result is built separately and
    /// `self` is never touched.
    pub fn subset_dependent(&self, tr_ids: &[String]) -> Result<Self> {
        let wanted = self.preflight_subset(tr_ids, false)?;
        let mut subset = Self::default();
        subset.set_transition_group_id(self.transition_group_id.clone());
        for transition in &self.transitions {
            let native_id = transition.native_id();
            if !wanted.contains(native_id) {
                continue;
            }
            let index = *self
                .chromatogram_map
                .get(native_id)
                .ok_or_else(|| missing("chromatogram", native_id))?;
            let chromatogram = self
                .chromatograms
                .get(index)
                .ok_or_else(|| missing("chromatogram", native_id))?
                .clone();
            let key = native_id.to_owned();
            subset.add_transition(transition.clone(), &key)?;
            subset.add_chromatogram(chromatogram, &key)?;
        }
        for feature in &self.mrm_features {
            subset.add_feature(feature.clone())?;
        }
        Ok(subset)
    }

    /// Preflight shared by both subset operations; returns the selection set.
    fn preflight_subset<'a>(
        &self,
        tr_ids: &'a [String],
        with_precursors: bool,
    ) -> Result<BTreeSet<&'a str>> {
        if tr_ids.len() > MAX_ITEMS {
            return Err(limit("identifier count"));
        }
        let mut items = self.transitions.len();
        items = items
            .checked_add(self.chromatograms.len())
            .and_then(|n| n.checked_add(self.mrm_features.len()))
            .ok_or_else(|| limit("list length"))?;
        if with_precursors {
            items = items
                .checked_add(self.precursor_chromatograms.len())
                .ok_or_else(|| limit("list length"))?;
        }
        if items > MAX_ITEMS {
            return Err(limit("list length"));
        }
        let per_item = std::mem::size_of::<C>()
            .max(std::mem::size_of::<T>())
            .max(std::mem::size_of::<MRMFeature>())
            .checked_add(2 * std::mem::size_of::<(String, usize)>() + 128)
            .ok_or_else(|| limit("payload"))?;
        let bytes = items
            .checked_mul(per_item)
            .ok_or_else(|| limit("payload"))?;
        if bytes > Self::MAX_BYTES {
            return Err(limit("payload"));
        }
        Ok(tr_ids.iter().map(String::as_str).collect())
    }
}

impl<C: Clone + NativeIdentified, T: Transition + Clone> MRMTransitionGroup<C, T> {
    /// A new group holding only the named transitions and their chromatograms,
    /// with every precursor chromatogram and a reduced copy of every peak
    /// group.
    ///
    /// Source `subset`. The transition group ID is copied. Transitions are
    /// selected by their own native ID, and each selected transition and
    /// chromatogram is re-registered under that native ID rather than under the
    /// key it originally had. Both transfers are guarded independently, so a
    /// transition registered under a key other than its native ID is dropped
    /// while a chromatogram registered under the native ID is still carried
    /// over.
    ///
    /// Every precursor chromatogram is carried over regardless of `tr_ids`,
    /// re-keyed by its *own* native ID rather than by the key it was stored
    /// under.
    ///
    /// Each peak group is rebuilt rather than copied: only intensity, retention
    /// time and metadata are taken from the original, so quality, m/z, charge,
    /// width, unique ID, convex hulls, subordinates, peptide identifications
    /// and the [`OpenSwathScores`] record are left at their defaults. Its
    /// per-transition features are copied for the selected transitions only,
    /// and all of its precursor features are copied.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when a peak group has no feature for a
    /// selected transition (the source throws `std::out_of_range` from
    /// `std::map::at`), [`Error::InvalidValue`] when two selected transitions
    /// or two precursor chromatograms share a native ID (the source throws
    /// `Exception::InvalidValue`), and [`Error::InvalidValue`] when the input
    /// exceeds [`MRMTransitionGroup::MAX_ITEMS`] or
    /// [`MRMTransitionGroup::MAX_BYTES`]. The result is built separately and
    /// `self` is never touched.
    pub fn subset(&self, tr_ids: &[String]) -> Result<Self> {
        let wanted = self.preflight_subset(tr_ids, true)?;
        let mut subset = Self::default();
        subset.set_transition_group_id(self.transition_group_id.clone());

        for transition in &self.transitions {
            let native_id = transition.native_id();
            if !wanted.contains(native_id) {
                continue;
            }
            let key = native_id.to_owned();
            if self.transition_map.contains_key(&key) {
                subset.add_transition(transition.clone(), &key)?;
            }
            if let Some(&index) = self.chromatogram_map.get(&key) {
                let chromatogram = self
                    .chromatograms
                    .get(index)
                    .ok_or_else(|| missing("chromatogram", &key))?
                    .clone();
                subset.add_chromatogram(chromatogram, &key)?;
            }
        }

        for precursor in &self.precursor_chromatograms {
            let key = precursor.native_id().to_owned();
            subset.add_precursor_chromatogram(precursor.clone(), &key)?;
        }

        for source in &self.mrm_features {
            let mut feature = MRMFeature::from(Feature {
                base: BaseFeature {
                    rt: source.feature.base.rt,
                    intensity: source.feature.base.intensity,
                    metadata: source.feature.base.metadata.clone(),
                    ..BaseFeature::default()
                },
                ..Feature::default()
            });
            for transition in &self.transitions {
                let native_id = transition.native_id();
                if !wanted.contains(native_id) {
                    continue;
                }
                feature.add_feature(source.feature(native_id)?.clone(), native_id)?;
            }
            let precursor_ids: Vec<String> =
                source.precursor_feature_ids().map(str::to_owned).collect();
            for id in &precursor_ids {
                feature.add_precursor_feature(source.precursor_feature(id)?.clone(), id)?;
            }
            subset.add_feature(feature)?;
        }

        Ok(subset)
    }
}

impl<C, T> MRMTransitionGroup<C, T> {
    /// The peak group with the highest overall quality.
    ///
    /// Source `getBestFeature`, which compares `Feature::getOverallQuality()`
    /// — the `quality` field — with a strict `>`, so the first of several equal
    /// maxima wins. A non-finite quality never compares greater and so never
    /// wins; the source behaves the same way.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when the group holds no feature.
    /// The source states the requirement as an `OPENMS_PRECONDITION`, which is
    /// compiled out in release builds, where it then reads `getFeatures()[0]`
    /// on an empty vector.
    pub fn best_feature(&self) -> Result<&MRMFeature> {
        let mut best = self.mrm_features.first().ok_or_else(|| {
            Error::MissingInformation("cannot get best feature for empty transition group".into())
        })?;
        for feature in &self.mrm_features {
            if feature.feature.base.quality > best.feature.base.quality {
                best = feature;
            }
        }
        Ok(best)
    }
}
