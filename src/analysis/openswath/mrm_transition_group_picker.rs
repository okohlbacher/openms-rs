// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Peak-group picking across the chromatograms of one precursor.
//!
//! Ported from `ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h` and `.cpp` at
//! OpenMS4-core `bc9cc12`. [`MRMTransitionGroupPicker::pick_transition_group`]
//! takes an [`MRMTransitionGroup`] filled with chromatograms and
//!
//! 1. finds peaks in the individual chromatograms (smoothing and
//!    [`PeakPickerChromatogram`], which seeds with the `PeakPickerHiRes`), and
//! 2. merges them into consensus features that span the chromatograms: the
//!    largest (or widest) remaining peak defines the boundaries, every
//!    chromatogram is resampled onto the sampling of the seed chromatogram and
//!    integrated with [`PeakIntegrator`] between them, and the process repeats
//!    until no peak is left.
//!
//! The arithmetic follows the Linux x86-64 Release build: resampling
//! redistributes intensities with the source's own `dist_right / (dist_left +
//! dist_right)` weights and `f32` accumulation, overlapping and quality
//! computations keep the source's evaluation order, and three source defects
//! that change results are reproduced and documented at the items that carry
//! them. See `docs/MRM_TRANSITION_GROUP_PICKER_SUPPORT.md`.
//!
//! [`MRMTransitionGroupPicker::pick_transition_group`]: crate::analysis::openswath::mrm_transition_group_picker::MRMTransitionGroupPicker::pick_transition_group
//! [`MRMTransitionGroup`]: crate::kernel::mrm::MRMTransitionGroup
//! [`PeakPickerChromatogram`]: crate::processing::chromatogram::PeakPickerChromatogram
//! [`PeakIntegrator`]: crate::analysis::peak_integrator::PeakIntegrator

use crate::analysis::openswath::scoring::{
    compute_and_append_rank, normalized_cross_correlation, ranked_mutual_information,
    xcorr_array_get_max_peak,
};
use crate::analysis::openswath::stats_helpers::MeanAndStddev;
use crate::analysis::peak_integrator::{PeakArea, PeakIntegrator};
use crate::concept::unique_id::{HasUniqueId, UniqueIdGenerator};
use crate::kernel::mrm::{DuplicateKeyPolicy, MRMFeature, MRMTransitionGroup, Transition};
use crate::kernel::{ChromatogramPeak, ConvexHull2D, Feature, MSChromatogram, Point2D};
use crate::metadata::MetaValue;
use crate::param::{DefaultParamHandler, Param, ParamNode};
use crate::processing::chromatogram::PeakPickerChromatogram;
use crate::processing::peak_picking::{
    flag, flag_entry, float, float_entry, int_entry, string_entry,
};
use crate::{Error, Result};

/// The source class name, `DefaultParamHandler("MRMTransitionGroupPicker")`,
/// under which unknown-parameter warnings are reported.
pub const MRM_TRANSITION_GROUP_PICKER_NAME: &str = "MRMTransitionGroupPicker";

/// The metadata keys the picker writes (source namespace
/// `MRMTransitionGroupPickerMeta`).
///
/// The source registers each name once with `MetaInfo::registry()` and hands
/// out the numeric index; the port's metadata is keyed by name, so these are
/// the names themselves.
pub mod meta {
    /// `area_background_level`: background area subtracted from a transition.
    pub const AREA_BACKGROUND_LEVEL: &str = "area_background_level";
    /// `asymmetry_factor` of a transition peak.
    pub const ASYMMETRY_FACTOR: &str = "asymmetry_factor";
    /// `baseline_delta_2_height` of a transition peak.
    pub const BASELINE_DELTA_2_HEIGHT: &str = "baseline_delta_2_height";
    /// `end_position_at_5` of a transition peak.
    pub const END_POSITION_AT_5: &str = "end_position_at_5";
    /// `end_position_at_10` of a transition peak.
    pub const END_POSITION_AT_10: &str = "end_position_at_10";
    /// `end_position_at_50` of a transition peak.
    pub const END_POSITION_AT_50: &str = "end_position_at_50";
    /// `initialPeakQuality` of a peak group.
    pub const INITIAL_PEAK_QUALITY: &str = "initialPeakQuality";
    /// `leftWidth`: left boundary of a peak group.
    pub const LEFT_WIDTH: &str = "leftWidth";
    /// `native_id` of the chromatogram a transition feature comes from.
    pub const NATIVE_ID: &str = "native_id";
    /// `noise_background_level` subtracted from a transition's apex.
    pub const NOISE_BACKGROUND_LEVEL: &str = "noise_background_level";
    /// `peak_apex_int`: apex intensity of a transition peak.
    pub const PEAK_APEX_INT: &str = "peak_apex_int";
    /// `peak_apex_position` of a transition peak.
    pub const PEAK_APEX_POSITION: &str = "peak_apex_position";
    /// `peak_apices_sum` over the quantifying transitions of a peak group.
    pub const PEAK_APICES_SUM: &str = "peak_apices_sum";
    /// `PeptideRef`: the transition group ID of a peak group.
    pub const PEPTIDE_REF: &str = "PeptideRef";
    /// `points_across_baseline` of a transition peak.
    pub const POINTS_ACROSS_BASELINE: &str = "points_across_baseline";
    /// `points_across_half_height` of a transition peak.
    pub const POINTS_ACROSS_HALF_HEIGHT: &str = "points_across_half_height";
    /// `potentialOutlier`: native ID of a suspect transition, or `none`.
    pub const POTENTIAL_OUTLIER: &str = "potentialOutlier";
    /// `rightWidth`: right boundary of a peak group.
    pub const RIGHT_WIDTH: &str = "rightWidth";
    /// `slope_of_baseline` of a transition peak.
    pub const SLOPE_OF_BASELINE: &str = "slope_of_baseline";
    /// `start_position_at_5` of a transition peak.
    pub const START_POSITION_AT_5: &str = "start_position_at_5";
    /// `start_position_at_10` of a transition peak.
    pub const START_POSITION_AT_10: &str = "start_position_at_10";
    /// `start_position_at_50` of a transition peak.
    pub const START_POSITION_AT_50: &str = "start_position_at_50";
    /// `tailing_factor` of a transition peak.
    pub const TAILING_FACTOR: &str = "tailing_factor";
    /// `total_mi`: total mutual information of a transition or peak group.
    pub const TOTAL_MI: &str = "total_mi";
    /// `total_width` of a transition peak.
    pub const TOTAL_WIDTH: &str = "total_width";
    /// `total_xic`: summed chromatogram intensity of a transition or of the
    /// detecting transitions of a peak group.
    pub const TOTAL_XIC: &str = "total_xic";
    /// `width_at_5` of a transition peak.
    pub const WIDTH_AT_5: &str = "width_at_5";
    /// `width_at_10` of a transition peak.
    pub const WIDTH_AT_10: &str = "width_at_10";
    /// `width_at_50` of a transition peak.
    pub const WIDTH_AT_50: &str = "width_at_50";
}

/// Which chromatogram peak areas and heights are integrated on (source
/// parameter `peak_integration`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PeakIntegrationSignal {
    /// The raw chromatogram (`original`, the default).
    #[default]
    Original,
    /// The smoothed chromatogram the peak picker produced (`smoothed`).
    Smoothed,
}

impl PeakIntegrationSignal {
    /// The source parameter value.
    pub const fn name(self) -> &'static str {
        match self {
            Self::Original => "original",
            Self::Smoothed => "smoothed",
        }
    }
}

/// How background is removed from each transition's peak (source parameter
/// `background_subtraction`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum BackgroundSubtraction {
    /// No background removal (`none`, the default).
    #[default]
    None,
    /// The historical estimate (`original`): the mean of the raw
    /// chromatogram's two boundary intensities, times the number of raw
    /// samples between the boundaries.
    Original,
    /// [`PeakIntegrator`]'s background estimate on the integrated signal
    /// (`exact`).
    Exact,
}

impl BackgroundSubtraction {
    /// The source parameter value.
    pub const fn name(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Original => "original",
            Self::Exact => "exact",
        }
    }
}

/// How the seed peak of each peak group is chosen (source parameter
/// `boundary_selection_method`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum BoundarySelection {
    /// The most intense remaining picked peak (`largest`, the default).
    #[default]
    Largest,
    /// The remaining picked peak with the widest boundaries (`widest`).
    Widest,
}

impl BoundarySelection {
    /// The source parameter value.
    pub const fn name(self) -> &'static str {
        match self {
            Self::Largest => "largest",
            Self::Widest => "widest",
        }
    }
}

/// Sums accumulated over the chromatograms of one peak group: the source's
/// `total_intensity`, `total_xic`, `total_mi` and `total_peak_apices`
/// out-parameters of `pickFragmentChromatograms`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PeakGroupTotals {
    /// Summed integrated intensity of the quantifying transitions.
    pub total_intensity: f64,
    /// Summed chromatogram intensity of the detecting transitions.
    pub total_xic: f64,
    /// Mutual information summed over the detecting transitions.
    pub total_mi: f64,
    /// Summed apex intensity of the quantifying transitions.
    pub total_peak_apices: f64,
}

/// Per-group data that every peak group of one transition group reuses: the
/// precomputed arguments of the source's optimised `createMRMFeature`
/// overload.
#[derive(Clone, Debug)]
pub struct TransitionGroupInputs<'a> {
    /// Total chromatogram intensity per transition, in transition order.
    pub transition_total_xics: Vec<f64>,
    /// Sum of [`transition_total_xics`](Self::transition_total_xics) over the
    /// detecting transitions.
    pub detecting_total_xic: f64,
    /// The raw chromatogram of each transition, in transition order.
    pub transition_chromatograms: Vec<&'a MSChromatogram>,
    /// The raw chromatogram behind each picked chromatogram, in picked order.
    pub picked_input_chromatograms: Vec<&'a MSChromatogram>,
}

impl<'a> TransitionGroupInputs<'a> {
    /// Precompute the inputs for `group` and its `picked` chromatograms, as the
    /// source's compatibility `createMRMFeature` overload does: the raw
    /// chromatogram behind each picked one is found by its native ID.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a transition or picked chromatogram names a
    /// chromatogram the group does not hold (source `selectChromHelper_`'s
    /// `Exception::IllegalArgument`).
    pub fn new<T: Transition>(
        group: &'a MRMTransitionGroup<MSChromatogram, T>,
        picked: &[MSChromatogram],
    ) -> Result<Self> {
        let mut inputs = Self::for_transitions(group)?;
        for chromatogram in picked {
            inputs
                .picked_input_chromatograms
                .push(select_chrom_helper(group, &chromatogram.native_id)?);
        }
        Ok(inputs)
    }

    fn for_transitions<T: Transition>(
        group: &'a MRMTransitionGroup<MSChromatogram, T>,
    ) -> Result<Self> {
        let mut inputs = Self {
            transition_total_xics: Vec::with_capacity(group.transitions().len()),
            detecting_total_xic: 0.0,
            transition_chromatograms: Vec::with_capacity(group.transitions().len()),
            picked_input_chromatograms: Vec::new(),
        };
        for transition in group.transitions() {
            let chromatogram = select_chrom_helper(group, transition.native_id())?;
            inputs.transition_chromatograms.push(chromatogram);
            let total_xic = chromatogram
                .peaks
                .iter()
                .fold(0.0, |total, peak| total + f64::from(peak.intensity));
            inputs.transition_total_xics.push(total_xic);
            if transition.is_detecting() {
                inputs.detecting_total_xic += total_xic;
            }
        }
        Ok(inputs)
    }
}

/// Finds peak groups in the chromatograms that belong to the same precursor
/// (source `MRMTransitionGroupPicker`).
///
/// Step 1 finds features (peaks) in the individual chromatograms by smoothing
/// them and applying the `PeakPickerHiRes`, through
/// [`picker`](Self::picker). Step 2 merges them into consensus features that
/// span several chromatograms, by finding the largest peak overall and
/// propagating its boundaries through all chromatograms.
///
/// The source class is a `DefaultParamHandler`; its parameters, with the
/// `PeakPickerChromatogram:` and `PeakIntegrator:` subsections, are
/// [`MRMTransitionGroupPicker::defaults`], [`MRMTransitionGroupPicker::to_param`]
/// and [`MRMTransitionGroupPicker::from_param`]. The source's protected copy
/// assignment operator copies nothing ("don't copy parameters"); [`Clone`]
/// copies every field. The source's thread-local allocation pools
/// (`PickerPoolShared`) only avoid reallocation and have no counterpart.
#[derive(Clone, Debug)]
pub struct MRMTransitionGroupPicker {
    /// Source `peak_integration` (advanced): which signal is integrated.
    pub peak_integration: PeakIntegrationSignal,
    /// Source `background_subtraction` (advanced). The `original` method is
    /// only provided for historical purposes; the source advises `exact`, with
    /// the settings of [`integrator`](Self::integrator). The chromatogram
    /// [`peak_integration`](Self::peak_integration) selects is used for the
    /// estimate.
    pub background_subtraction: BackgroundSubtraction,
    /// Source `recalculate_peaks` (advanced, default `false`): look at the
    /// consistency of all picked peaks and use the consensus (median) border
    /// when the proposed one deviates too much. Consensus mode only.
    pub recalculate_peaks: bool,
    /// Source `use_precursors` (advanced, default `false`): also pick the
    /// precursor chromatograms, so that precursor signal may drive the
    /// picking.
    pub use_precursors: bool,
    /// Source `use_consensus` (advanced, default `true`): use the seed's
    /// boundaries for every transition; otherwise each transition keeps the
    /// boundaries of its own peak closest to the seed apex.
    pub use_consensus: bool,
    /// Source `compute_peak_quality` (advanced, default `false`): score each
    /// peak group's co-elution and shape consistency and flag an outlier
    /// transition. The score centres on zero; above zero is generally good,
    /// below -1 or -2 usually bad.
    pub compute_peak_quality: bool,
    /// Source `compute_peak_shape_metrics` (advanced, default `false`): store
    /// the shape metrics of each transition peak for QC. `width_at_50` is
    /// stored regardless, as the source does for its TOPP tests.
    pub compute_peak_shape_metrics: bool,
    /// Source `compute_total_mi` (advanced, default `false`): store the
    /// mutual information of each transition with the detecting transitions,
    /// for OpenSWATH and IPF scoring.
    pub compute_total_mi: bool,
    /// Source `minimal_quality` (advanced, default -10000): with
    /// [`compute_peak_quality`](Self::compute_peak_quality), peak groups below
    /// this quality are dropped.
    pub minimal_quality: f64,
    /// Source `stop_after_feature` (default -1): stop once more than this many
    /// peak groups were found; zero or less never stops. The source's
    /// `cnt > stop_after_feature` test collects one more peak group than the
    /// number given, and so does this port.
    pub stop_after_feature: i32,
    /// Source `stop_after_intensity_ratio` (default 0.0001): stop after a peak
    /// group whose intensity, relative to the detecting transitions' total
    /// chromatogram intensity, is below this ratio.
    pub stop_after_intensity_ratio: f64,
    /// Source `min_peak_width` (advanced, default 0.001 s): in consensus mode,
    /// peak groups narrower than this are dropped; zero or less disables it.
    pub min_peak_width: f64,
    /// Source `recalculate_peaks_max_z` (advanced, default 1.0): the z-score of
    /// a boundary above which [`recalculate_peaks`](Self::recalculate_peaks)
    /// falls back to the median.
    pub recalculate_peaks_max_z: f64,
    /// Source `resample_boundary` (advanced, default 15 s): how far beyond the
    /// peak the chromatograms are sampled for the peak quality.
    pub resample_boundary: f64,
    /// Source `boundary_selection_method` (advanced).
    pub boundary_selection_method: BoundarySelection,
    /// The chromatogram peak picker, source subsection
    /// `PeakPickerChromatogram:`.
    pub picker: PeakPickerChromatogram,
    /// The integrator, source subsection `PeakIntegrator:`.
    pub integrator: PeakIntegrator,
    /// Native ceiling on the number of fragment and precursor chromatograms
    /// in a group, checked before any work. Every other bound is enforced by
    /// [`picker`](Self::picker) and [`integrator`](Self::integrator).
    pub max_chromatograms: usize,
}

impl Default for MRMTransitionGroupPicker {
    fn default() -> Self {
        Self {
            peak_integration: PeakIntegrationSignal::Original,
            background_subtraction: BackgroundSubtraction::None,
            recalculate_peaks: false,
            use_precursors: false,
            use_consensus: true,
            compute_peak_quality: false,
            compute_peak_shape_metrics: false,
            compute_total_mi: false,
            minimal_quality: -10000.0,
            stop_after_feature: -1,
            stop_after_intensity_ratio: 0.0001,
            min_peak_width: 0.001,
            recalculate_peaks_max_z: 1.0,
            resample_boundary: 15.0,
            boundary_selection_method: BoundarySelection::Largest,
            picker: PeakPickerChromatogram::default(),
            integrator: PeakIntegrator::default(),
            max_chromatograms: 100_000,
        }
    }
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

fn set_float(feature: &mut Feature, key: &str, value: f64) -> Result<()> {
    let value = MetaValue::try_from(value)
        .map_err(|_| invalid(&format!("meta value '{key}' is not finite")))?;
    feature.base.metadata.insert(key.to_owned(), value);
    Ok(())
}

fn set_value(feature: &mut Feature, key: &str, value: MetaValue) {
    feature.base.metadata.insert(key.to_owned(), value);
}

fn float_meta(feature: &Feature, key: &str) -> Result<f64> {
    feature
        .base
        .metadata
        .get(key)
        .ok_or_else(|| Error::MissingInformation(format!("meta value '{key}' is missing")))?
        .as_f64()
}

fn count(value: usize) -> MetaValue {
    MetaValue::from(i64::try_from(value).unwrap_or(i64::MAX))
}

/// One float data array entry of a picked chromatogram (source
/// `getFloatDataArrays()[index][peak]`).
fn array_value(chromatogram: &MSChromatogram, index: usize, peak: usize) -> Result<f64> {
    chromatogram
        .float_data_arrays
        .get(index)
        .and_then(|array| array.data.get(peak))
        .map(|&value| f64::from(value))
        .ok_or_else(|| {
            invalid(&format!(
                "picked chromatogram '{}' has no float data array {index} value for peak {peak}",
                chromatogram.native_id
            ))
        })
}

/// Source `selectChromHelper_`: the fragment chromatogram with `native_id`,
/// else the precursor chromatogram with it.
fn select_chrom_helper<'a, T>(
    group: &'a MRMTransitionGroup<MSChromatogram, T>,
    native_id: &str,
) -> Result<&'a MSChromatogram> {
    if group.has_chromatogram(native_id) {
        group.chromatogram(native_id)
    } else if group.has_precursor_chromatogram(native_id) {
        group.precursor_chromatogram(native_id)
    } else {
        Err(invalid(&format!(
            "Did not find chromatogram for id '{native_id}'."
        )))
    }
}

/// libstdc++ `std::lower_bound` by position: well defined, if meaningless, on
/// a container that is not sorted by position.
fn lower_bound(peaks: &[ChromatogramPeak], value: f64) -> usize {
    let (mut first, mut len) = (0, peaks.len());
    while len > 0 {
        let half = len >> 1;
        let middle = first + half;
        if peaks[middle].rt < value {
            first = middle + 1;
            len = len - half - 1;
        } else {
            len = half;
        }
    }
    first
}

/// libstdc++ `std::upper_bound` by position; see [`lower_bound`].
fn upper_bound(peaks: &[ChromatogramPeak], value: f64) -> usize {
    let (mut first, mut len) = (0, peaks.len());
    while len > 0 {
        let half = len >> 1;
        let middle = first + half;
        if value < peaks[middle].rt {
            len = half;
        } else {
            first = middle + 1;
            len = len - half - 1;
        }
    }
    first
}

/// Source `MSChromatogram::sortByIntensity()`: a stable ascending sort of the
/// peaks that permutes every parallel data array with them.
fn sort_by_intensity(chromatogram: &mut MSChromatogram) -> Result<()> {
    if chromatogram
        .peaks
        .windows(2)
        .all(|pair| pair[0].intensity <= pair[1].intensity)
    {
        return Ok(());
    }
    chromatogram.validate()?;
    let mut order: Vec<usize> = (0..chromatogram.peaks.len()).collect();
    let peaks = &chromatogram.peaks;
    order.sort_by(|&a, &b| {
        if peaks[a].intensity < peaks[b].intensity {
            std::cmp::Ordering::Less
        } else if peaks[b].intensity < peaks[a].intensity {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    });
    fn permute<T: Clone>(values: &mut Vec<T>, order: &[usize]) {
        if values.len() == order.len() {
            *values = order.iter().map(|&i| values[i].clone()).collect();
        }
    }
    permute(&mut chromatogram.peaks, &order);
    for array in &mut chromatogram.float_data_arrays {
        permute(&mut array.data, &order);
    }
    for array in &mut chromatogram.integer_data_arrays {
        permute(&mut array.data, &order);
    }
    for array in &mut chromatogram.string_data_arrays {
        permute(&mut array.data, &order);
    }
    Ok(())
}

/// The source's inclusive-plus-one-neighbour sample range around
/// `[left, right]`, shared by `prepareMasterContainer_` and
/// `resampleChromatogram_`.
fn padded_range(peaks: &[ChromatogramPeak], left: f64, right: f64) -> (usize, usize) {
    let mut begin = 0;
    while begin < peaks.len() && peaks[begin].rt < left {
        begin += 1;
    }
    begin = begin.saturating_sub(1);
    let mut end = begin;
    while end < peaks.len() && peaks[end].rt < right {
        end += 1;
    }
    if end != peaks.len() {
        end += 1;
    }
    (begin, end)
}

/// Source `LinearResamplerAlign::raster(raw_begin, raw_end, out_begin,
/// out_end)`, reproduced step for step: intensities before the first grid
/// point are added to it in `f32`, each interior sample is split between the
/// two grid points around it with weights `dist_right / (dist_left +
/// dist_right)` and `dist_left / (...)` evaluated in `f64` and stored as
/// `f32`, and the rest is added to the last grid point in `f32`.
///
/// `crate::processing::resample_to_grid` computes the same split as
/// `1 - fraction` and `fraction`, which rounds differently; the peak-group
/// areas are compared against the Release build, so the source form is kept
/// here.
fn raster(raw: &[ChromatogramPeak], out: &mut [ChromatogramPeak]) {
    if out.is_empty() {
        return;
    }
    let mut r = 0;
    while r < raw.len() && raw[r].rt < out[0].rt {
        out[0].intensity += raw[r].intensity;
        r += 1;
    }
    let mut b = 0;
    while r < raw.len() {
        while b < out.len() && out[b].rt < raw[r].rt {
            b += 1;
        }
        b = b.saturating_sub(1);
        if b + 1 == out.len() {
            break;
        }
        let dist_left = (raw[r].rt - out[b].rt).abs();
        let dist_right = (raw[r].rt - out[b + 1].rt).abs();
        let intensity = f64::from(raw[r].intensity);
        out[b].intensity = (f64::from(out[b].intensity)
            + intensity * dist_right / (dist_left + dist_right)) as f32;
        out[b + 1].intensity = (f64::from(out[b + 1].intensity)
            + intensity * dist_left / (dist_left + dist_right))
            as f32;
        r += 1;
    }
    while r < raw.len() {
        out[b].intensity += raw[r].intensity;
        r += 1;
    }
}

impl MRMTransitionGroupPicker {
    /// The source defaults (`getDefaults`): fifteen parameters in declaration
    /// order, all but the two `stop_after_*` ones tagged advanced, followed by
    /// the `PeakPickerChromatogram:` and `PeakIntegrator:` subsections at
    /// those classes' defaults.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters of this picker in the layout of
    /// [`MRMTransitionGroupPicker::defaults`], carrying the current values
    /// (source `getParameters`), including those of
    /// [`picker`](Self::picker) and [`integrator`](Self::integrator).
    ///
    /// # Errors
    ///
    /// As [`PeakPickerChromatogram::to_param`].
    pub fn to_param(&self) -> Result<Param> {
        let flag_adv =
            |name: &str, value: bool, description: &str| flag_entry(name, value, description, true);
        let entries = vec![
            int_entry(
                "stop_after_feature",
                i64::from(self.stop_after_feature),
                "Stop finding after feature (ordered by intensity; -1 means do not stop).",
                false,
                None,
                None,
            ),
            float_entry(
                "stop_after_intensity_ratio",
                self.stop_after_intensity_ratio,
                "Stop after reaching intensity ratio",
                false,
                None,
                None,
            ),
            float_entry(
                "min_peak_width",
                self.min_peak_width,
                "Minimal peak width (s), discard all peaks below this value (-1 means no action).",
                true,
                None,
                None,
            ),
            string_entry(
                "peak_integration",
                self.peak_integration.name(),
                "Calculate the peak area and height either the smoothed or the raw chromatogram data.",
                true,
                &["original", "smoothed"],
            ),
            string_entry(
                "background_subtraction",
                self.background_subtraction.name(),
                "Remove background from peak signal using estimated noise levels. The 'original' method is only provided for historical purposes, please use the 'exact' method and set parameters using the PeakIntegrator: settings. The same original or smoothed chromatogram specified by peak_integration will be used for background estimation.",
                true,
                &["none", "original", "exact"],
            ),
            flag_adv(
                "recalculate_peaks",
                self.recalculate_peaks,
                "Tries to get better peak picking by looking at peak consistency of all picked peaks. Tries to use the consensus (median) peak border if the variation within the picked peaks is too large.",
            ),
            flag_adv(
                "use_precursors",
                self.use_precursors,
                "Use precursor chromatogram for peak picking (note that this may lead to precursor signal driving the peak picking)",
            ),
            flag_adv(
                "use_consensus",
                self.use_consensus,
                "Use consensus peak boundaries when computing transition group picking (if false, compute independent peak boundaries for each transition)",
            ),
            float_entry(
                "recalculate_peaks_max_z",
                self.recalculate_peaks_max_z,
                "Determines the maximal Z-Score (difference measured in standard deviations) that is considered too large for peak boundaries. If the Z-Score is above this value, the median is used for peak boundaries (default value 1.0).",
                true,
                None,
                None,
            ),
            float_entry(
                "minimal_quality",
                self.minimal_quality,
                "Only if compute_peak_quality is set, this parameter will not consider peaks below this quality threshold",
                true,
                None,
                None,
            ),
            float_entry(
                "resample_boundary",
                self.resample_boundary,
                "For computing peak quality, how many extra seconds should be sample left and right of the actual peak",
                true,
                None,
                None,
            ),
            flag_adv(
                "compute_peak_quality",
                self.compute_peak_quality,
                "Tries to compute a quality value for each peakgroup and detect outlier transitions. The resulting score is centered around zero and values above 0 are generally good and below -1 or -2 are usually bad.",
            ),
            flag_adv(
                "compute_peak_shape_metrics",
                self.compute_peak_shape_metrics,
                "Calculates various peak shape metrics (e.g., tailing) that can be used for downstream QC/QA.",
            ),
            flag_adv(
                "compute_total_mi",
                self.compute_total_mi,
                "Compute mutual information metrics for individual transitions that can be used for OpenSWATH/IPF scoring.",
            ),
            string_entry(
                "boundary_selection_method",
                self.boundary_selection_method.name(),
                "Method to use when selecting the best boundaries for peaks.",
                true,
                &["largest", "widest"],
            ),
        ];
        let mut param = Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries,
            nodes: Vec::new(),
        })?;
        param.insert("PeakPickerChromatogram:", &self.picker.to_param()?)?;
        param.insert("PeakIntegrator:", &self.integrator.to_param()?)?;
        Ok(param)
    }

    /// Build a picker from parameters, as source `setParameters` followed by
    /// `updateMembers_`, which hands the two subsections to
    /// [`PeakPickerChromatogram::from_param`] and
    /// [`PeakIntegrator::from_param`].
    ///
    /// Missing parameters take their defaults; unknown ones are ignored, as the
    /// source only warns (see
    /// [`MRMTransitionGroupPicker::from_param_with_warnings`]).
    /// [`max_chromatograms`](Self::max_chromatograms) keeps its default.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong value type
    /// or violates its restriction (the source `Exception::InvalidParameter`),
    /// or `stop_after_feature` does not fit 32 bits, and the errors of the two
    /// nested `from_param` calls.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`MRMTransitionGroupPicker::from_param`], also returning the
    /// unknown-parameter warnings the source logs. As in the source, an
    /// unknown key inside a subsection is reported twice: once by this class
    /// and once by the subsection's own `setParameters`.
    ///
    /// # Errors
    ///
    /// As [`MRMTransitionGroupPicker::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(MRM_TRANSITION_GROUP_PICKER_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        let ((picker, nested), mut warnings) =
            handler.set_parameters_with(parameters, Self::from_complete_param)?;
        warnings.extend(nested);
        Ok((picker, warnings))
    }

    fn from_complete_param(param: &Param) -> Result<(Self, Vec<String>)> {
        let (picker, mut warnings) = PeakPickerChromatogram::from_param_with_warnings(
            &param.copy("PeakPickerChromatogram:", true)?,
        )?;
        let (integrator, integrator_warnings) =
            PeakIntegrator::from_param_with_warnings(&param.copy("PeakIntegrator:", true)?)?;
        warnings.extend(integrator_warnings);
        let text = |key: &str| -> Result<String> { Ok(param.value(key)?.as_str()?.to_owned()) };
        let peak_integration = match text("peak_integration")?.as_str() {
            "smoothed" => PeakIntegrationSignal::Smoothed,
            "original" => PeakIntegrationSignal::Original,
            other => {
                return Err(invalid(&format!(
                    "Peak integration chromatogram {other} is not a valid method for MRMTransitionGroupPicker"
                )));
            }
        };
        let background_subtraction = match text("background_subtraction")?.as_str() {
            "none" => BackgroundSubtraction::None,
            "original" => BackgroundSubtraction::Original,
            "exact" => BackgroundSubtraction::Exact,
            other => {
                return Err(invalid(&format!(
                    "unknown background_subtraction '{other}'"
                )));
            }
        };
        let boundary_selection_method = match text("boundary_selection_method")?.as_str() {
            "largest" => BoundarySelection::Largest,
            "widest" => BoundarySelection::Widest,
            other => {
                return Err(invalid(&format!(
                    "unknown boundary_selection_method '{other}'"
                )));
            }
        };
        Ok((
            Self {
                peak_integration,
                background_subtraction,
                recalculate_peaks: flag(param, "recalculate_peaks")?,
                use_precursors: flag(param, "use_precursors")?,
                use_consensus: flag(param, "use_consensus")?,
                compute_peak_quality: flag(param, "compute_peak_quality")?,
                compute_peak_shape_metrics: flag(param, "compute_peak_shape_metrics")?,
                compute_total_mi: flag(param, "compute_total_mi")?,
                minimal_quality: float(param, "minimal_quality")?,
                stop_after_feature: param.value("stop_after_feature")?.to_i32()?,
                stop_after_intensity_ratio: float(param, "stop_after_intensity_ratio")?,
                min_peak_width: float(param, "min_peak_width")?,
                recalculate_peaks_max_z: float(param, "recalculate_peaks_max_z")?,
                resample_boundary: float(param, "resample_boundary")?,
                boundary_selection_method,
                picker,
                integrator,
                ..Self::default()
            },
            warnings,
        ))
    }

    /// Pick a group of chromatograms belonging to the same peptide (source
    /// `pickTransitionGroup`).
    ///
    /// Identifies peaks in a set of chromatograms that belong to the same
    /// peptide. Only chromatograms of detecting transitions (and, with
    /// [`use_precursors`](Self::use_precursors), the precursor chromatograms)
    /// are used for peak picking; identifying transitions are processed
    /// alongside but do not contribute to the metadata, e.g. `total_xic` or
    /// `peak_apices_sum`. A chromatogram with no transition registered under
    /// its native ID is picked as detecting.
    ///
    /// The resulting peak groups are added to the group, most intense first
    /// (or widest first), skipping any whose `[leftWidth, rightWidth]` lies
    /// inside an earlier one's. Each carries `PeptideRef`, `leftWidth`,
    /// `rightWidth`, `total_xic` (the detecting fragment traces' summed
    /// intensity) and `peak_apices_sum`, plus `total_mi` and the quality
    /// values when enabled, and a unique ID drawn from `ids` (source
    /// `ensureUniqueId`, which draws from a process-wide generator).
    ///
    /// The source states `isInternallyConsistent` and `chromatogramIdsMatch`
    /// as preconditions that a Release build does not check; neither is
    /// checked here either. A transition whose chromatogram is missing is an
    /// error, as in the source.
    ///
    /// # Errors
    ///
    /// On any error the group is unchanged. [`Error::InvalidValue`] when the
    /// group holds more than [`max_chromatograms`](Self::max_chromatograms)
    /// chromatograms; when a transition names no chromatogram of the group;
    /// with [`use_consensus`](Self::use_consensus) off, when a transition is
    /// not detecting (the source's `Exception::IllegalArgument`, "all
    /// transitions need to be detecting transitions"); when smoothed
    /// integration lacks a smoothed chromatogram (source
    /// `Exception::IllegalArgument`); and when a computed metadata value is
    /// not finite. Also the errors of [`picker`](Self::picker) and
    /// [`integrator`](Self::integrator) — including the index and range cases
    /// where the source reads out of bounds — and those of
    /// [`MRMTransitionGroupPicker::create_mrm_feature_with_inputs`].
    pub fn pick_transition_group<T: Transition>(
        &self,
        group: &mut MRMTransitionGroup<MSChromatogram, T>,
        ids: &mut UniqueIdGenerator,
    ) -> Result<()> {
        let features = self.find_peak_groups(group, ids)?;
        let mut accepted = Vec::new();
        for (i, feature) in features.iter().enumerate() {
            let left = float_meta(feature, meta::LEFT_WIDTH)?;
            let right = float_meta(feature, meta::RIGHT_WIDTH)?;
            let mut skip = false;
            for earlier in &features[..i] {
                if left >= float_meta(earlier, meta::LEFT_WIDTH)?
                    && right <= float_meta(earlier, meta::RIGHT_WIDTH)?
                {
                    skip = true;
                }
            }
            if feature.base.intensity > 0.0 && !skip {
                accepted.push(i);
            }
        }
        if group.features().len() + accepted.len()
            > MRMTransitionGroup::<MSChromatogram, T>::MAX_ITEMS
        {
            return Err(invalid("transition group feature limit exceeded"));
        }
        let mut features: Vec<Option<MRMFeature>> = features.into_iter().map(Some).collect();
        for i in accepted {
            if let Some(feature) = features[i].take() {
                group.add_feature(feature)?;
            }
        }
        Ok(())
    }

    /// Everything `pickTransitionGroup` does before it adds the peak groups.
    fn find_peak_groups<T: Transition>(
        &self,
        group: &MRMTransitionGroup<MSChromatogram, T>,
        ids: &mut UniqueIdGenerator,
    ) -> Result<Vec<MRMFeature>> {
        let chromatogram_count = group
            .chromatograms()
            .len()
            .saturating_add(group.precursor_chromatograms().len());
        if chromatogram_count > self.max_chromatograms {
            return Err(invalid("transition group exceeds the chromatogram limit"));
        }
        let mut picked = Vec::new();
        let mut smoothed = Vec::new();
        let mut picked_inputs = Vec::new();
        for chromatogram in group.chromatograms() {
            let id = chromatogram.native_id.as_str();
            if !group.transitions().is_empty()
                && group.has_transition(id)
                && !group.transition(id)?.is_detecting()
            {
                continue;
            }
            let (p, s) = self.pick_one(chromatogram, true)?;
            picked.push(p);
            smoothed.push(s);
            picked_inputs.push(chromatogram);
        }
        if self.use_precursors {
            for chromatogram in group.precursor_chromatograms() {
                let (p, s) = self.pick_one(chromatogram, false)?;
                picked.push(p);
                smoothed.push(s);
                picked_inputs.push(chromatogram);
            }
        }
        let mut inputs = TransitionGroupInputs::for_transitions(group)?;
        inputs.picked_input_chromatograms = picked_inputs;

        let mut largest: Vec<Option<usize>> = picked
            .iter()
            .map(|chromatogram| chromatogram.peaks.len().checked_sub(1))
            .collect();
        let mut features = Vec::new();
        let mut found: i64 = 0;
        loop {
            let seed = match self.boundary_selection_method {
                BoundarySelection::Largest => {
                    Self::find_largest_peak_cached(&picked, &mut largest)?
                }
                BoundarySelection::Widest => Self::find_widest_peak_indices(&picked)?,
            };
            let Some((chr_idx, peak_idx)) = seed else {
                break;
            };
            let feature = self.create_mrm_feature_with_inputs(
                group,
                &mut picked,
                &smoothed,
                &inputs,
                chr_idx,
                peak_idx,
                ids,
            )?;
            let intensity = f64::from(feature.base.intensity);
            let mut total_xic = 0.0;
            if intensity > 0.0 {
                total_xic = float_meta(&feature, meta::TOTAL_XIC)?;
                features.push(feature);
                found += 1;
            }
            if self.stop_after_feature > 0 && found > i64::from(self.stop_after_feature) {
                break;
            }
            if intensity > 0.0 && intensity / total_xic < self.stop_after_intensity_ratio {
                break;
            }
        }
        Ok(features)
    }

    /// Pick one chromatogram as the source loop does: `pickChromatogram`, then
    /// `sortByIntensity` on the picked result.
    fn pick_one(
        &self,
        chromatogram: &MSChromatogram,
        fragment: bool,
    ) -> Result<(MSChromatogram, MSChromatogram)> {
        if chromatogram.peaks.is_empty() {
            // `pickChromatogram` returns before touching either output; only a
            // fragment's smoothed trace had its native ID set beforehand.
            let smoothed = MSChromatogram {
                native_id: if fragment {
                    chromatogram.native_id.clone()
                } else {
                    String::new()
                },
                ..MSChromatogram::default()
            };
            return Ok((MSChromatogram::default(), smoothed));
        }
        let result = self.picker.pick_chromatogram(chromatogram)?;
        let mut picked = result.picked.chromatogram;
        sort_by_intensity(&mut picked)?;
        Ok((picked, result.smoothed))
    }

    /// Create a peak group from picked chromatograms and a seed peak (source
    /// `createMRMFeature`, compatibility overload for callers that do not
    /// precompute the per-group inputs).
    ///
    /// `picked` are picked chromatograms (as
    /// [`PeakPickerChromatogram::pick_chromatogram`] produces them, with the
    /// `leftWidth`, `rightWidth` and `IntegratedIntensity` arrays at the
    /// source indices); peaks the new group covers are marked as consumed by
    /// setting their intensity to zero. `smoothed` are the smoothed
    /// chromatograms used for [`PeakIntegrationSignal::Smoothed`]. The seed is
    /// peak `peak_idx` of `picked[chr_idx]`. The raw chromatogram behind each
    /// picked one is found by its native ID.
    ///
    /// Returns the peak group, or one with zero intensity when a filter
    /// (minimal width, minimal quality) rejects it.
    ///
    /// # Errors
    ///
    /// As [`TransitionGroupInputs::new`] and
    /// [`MRMTransitionGroupPicker::create_mrm_feature_with_inputs`].
    pub fn create_mrm_feature<T: Transition>(
        &self,
        group: &MRMTransitionGroup<MSChromatogram, T>,
        picked: &mut [MSChromatogram],
        smoothed: &[MSChromatogram],
        chr_idx: usize,
        peak_idx: usize,
        ids: &mut UniqueIdGenerator,
    ) -> Result<MRMFeature> {
        let inputs = TransitionGroupInputs::new(group, picked)?;
        self.create_mrm_feature_with_inputs(
            group, picked, smoothed, &inputs, chr_idx, peak_idx, ids,
        )
    }

    /// Create a peak group from picked chromatograms and precomputed
    /// per-group data (source `createMRMFeature`, the overload the optimised
    /// picking path uses so that repeated peak groups reuse the same data).
    ///
    /// The seed's boundaries (`leftWidth`, `rightWidth`) become the group's.
    /// In consensus mode, with [`recalculate_peaks`](Self::recalculate_peaks),
    /// they may be replaced by the median of the other chromatograms' borders
    /// (and the apex moved to their centre if it falls outside); then every
    /// picked peak inside or overlapping them is consumed. Otherwise each
    /// chromatogram keeps the boundaries of its own peak closest to the seed
    /// apex ([`MRMTransitionGroupPicker::pick_apex`]). The seed is always
    /// consumed. The raw chromatogram of the seed defines the sampling that
    /// every chromatogram is resampled onto before integration.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `chr_idx`, `peak_idx` or a picked
    /// chromatogram's arrays do not hold the seed (the source reads out of
    /// bounds), when the seed's raw chromatogram is empty (the source then
    /// resamples into an empty range), when a metadata value is not finite,
    /// and the errors of [`MRMTransitionGroupPicker::pick_apex`],
    /// [`MRMTransitionGroupPicker::pick_fragment_chromatograms`] and
    /// [`MRMTransitionGroupPicker::pick_precursor_chromatograms`].
    #[allow(clippy::too_many_arguments)]
    pub fn create_mrm_feature_with_inputs<T: Transition>(
        &self,
        group: &MRMTransitionGroup<MSChromatogram, T>,
        picked: &mut [MSChromatogram],
        smoothed: &[MSChromatogram],
        inputs: &TransitionGroupInputs<'_>,
        chr_idx: usize,
        peak_idx: usize,
        ids: &mut UniqueIdGenerator,
    ) -> Result<MRMFeature> {
        let seed = picked
            .get(chr_idx)
            .ok_or_else(|| invalid("seed chromatogram index out of range"))?;
        let mut best_left = array_value(seed, PeakPickerChromatogram::IDX_LEFTBORDER, peak_idx)?;
        let mut best_right = array_value(seed, PeakPickerChromatogram::IDX_RIGHTBORDER, peak_idx)?;
        let mut peak_apex = seed
            .peaks
            .get(peak_idx)
            .ok_or_else(|| invalid("seed peak index out of range"))?
            .rt;
        let mut feature = MRMFeature::new();
        feature.base.intensity = 0.0;

        if self.use_consensus && self.recalculate_peaks {
            Self::recalculate_peak_borders(
                picked,
                &mut best_left,
                &mut best_right,
                self.recalculate_peaks_max_z,
            )?;
            if peak_apex < best_left || peak_apex > best_right {
                peak_apex = (best_left + best_right) / 2.0;
            }
        }

        let mut left_edges = Vec::new();
        let mut right_edges = Vec::new();
        let mut min_left = best_left;
        let mut max_right = best_right;
        if self.use_consensus {
            Self::remove_overlapping_features(picked, best_left, best_right);
        } else {
            self.pick_apex(
                picked,
                best_left,
                best_right,
                peak_apex,
                &mut min_left,
                &mut max_right,
                &mut left_edges,
                &mut right_edges,
            )?;
        }
        picked[chr_idx].peaks[peak_idx].intensity = 0.0;

        if self.use_consensus {
            if self.min_peak_width > 0.0 && (best_right - best_left).abs() < self.min_peak_width {
                return Ok(feature);
            }
            if self.compute_peak_quality {
                let (quality, outlier) =
                    self.compute_quality(picked, inputs, chr_idx, best_left, best_right)?;
                if quality < self.minimal_quality {
                    return Ok(feature);
                }
                set_value(
                    &mut feature,
                    meta::POTENTIAL_OUTLIER,
                    MetaValue::from(outlier),
                );
                set_float(&mut feature, meta::INITIAL_PEAK_QUALITY, quality)?;
                feature.base.quality = quality as f32;
            }
        }

        let reference = inputs
            .picked_input_chromatograms
            .get(chr_idx)
            .ok_or_else(|| invalid("no raw chromatogram for the seed chromatogram"))?;
        let master = Self::prepare_master_container(reference, min_left, max_right)?;

        let mut totals = PeakGroupTotals::default();
        self.pick_fragment_chromatograms(
            group,
            picked,
            &mut feature,
            smoothed,
            best_left,
            best_right,
            self.use_consensus,
            &mut totals,
            &master,
            &left_edges,
            &right_edges,
            inputs,
            chr_idx,
            peak_idx,
        )?;
        self.pick_precursor_chromatograms(
            group,
            picked,
            &mut feature,
            smoothed,
            best_left,
            best_right,
            self.use_consensus,
            &mut totals.total_intensity,
            &master,
            &left_edges,
            &right_edges,
            chr_idx,
            peak_idx,
        )?;

        feature.base.rt = peak_apex;
        feature.base.intensity = totals.total_intensity as f32;
        set_value(
            &mut feature,
            meta::PEPTIDE_REF,
            MetaValue::from(group.transition_group_id()),
        );
        set_float(&mut feature, meta::LEFT_WIDTH, best_left)?;
        set_float(&mut feature, meta::RIGHT_WIDTH, best_right)?;
        set_float(&mut feature, meta::TOTAL_XIC, totals.total_xic)?;
        if self.compute_total_mi {
            set_float(&mut feature, meta::TOTAL_MI, totals.total_mi)?;
        }
        set_float(
            &mut feature,
            meta::PEAK_APICES_SUM,
            totals.total_peak_apices,
        )?;
        feature.feature.base.unique_id.ensure_unique_id(ids);
        Ok(feature)
    }

    /// Apex-based peak selection for non-consensus mode (source `pickApex`).
    ///
    /// For each picked chromatogram, the peak whose apex lies closest to
    /// `peak_apex` is chosen and consumed (only one per chromatogram, so a
    /// second peak is found later), and its boundaries are appended to
    /// `left_edges` and `right_edges`; a chromatogram without a peak takes the
    /// seed's `best_left` and `best_right`. `min_left` and `max_right` widen to
    /// cover every chosen boundary.
    ///
    /// Each apex is the [`PeakIntegrator`] apex of the picked (centroided)
    /// chromatogram between the peak's own boundaries, and only apices above
    /// `1e-11` count. The picked chromatograms are sorted by intensity, not by
    /// position, and the source still bounds this integration with the
    /// binary searches `PosBegin` and `PosEnd`; they are reproduced exactly
    /// (libstdc++ `lower_bound` and `upper_bound`), so a peak is selected
    /// exactly as the Release build selects it.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a picked chromatogram lacks its boundary
    /// arrays, or when those searches cross (the source then iterates past the
    /// end of the container); [`Error::Unsupported`] when the integrator fits
    /// an EMG model, which the source would fit to the centroided peaks and
    /// this port does not reproduce.
    #[allow(clippy::too_many_arguments)]
    pub fn pick_apex(
        &self,
        picked: &mut [MSChromatogram],
        best_left: f64,
        best_right: f64,
        peak_apex: f64,
        min_left: &mut f64,
        max_right: &mut f64,
        left_edges: &mut Vec<f64>,
        right_edges: &mut Vec<f64>,
    ) -> Result<()> {
        if self.integrator.emg.is_some() {
            return Err(Error::Unsupported(
                "non-consensus apex selection with EMG fitting of centroided peaks".into(),
            ));
        }
        for chromatogram in picked.iter_mut() {
            let mut distance_min = f64::MAX;
            let mut closest: Option<usize> = None;
            for i in 0..chromatogram.peaks.len() {
                let left = array_value(chromatogram, PeakPickerChromatogram::IDX_LEFTBORDER, i)?;
                let right = array_value(chromatogram, PeakPickerChromatogram::IDX_RIGHTBORDER, i)?;
                let begin = lower_bound(&chromatogram.peaks, left);
                let end = upper_bound(&chromatogram.peaks, right);
                if begin > end {
                    return Err(invalid(
                        "picked-peak boundaries select a crossed range in an intensity-sorted chromatogram",
                    ));
                }
                let mut apex = PeakArea {
                    apex_pos: (left + right) / 2.0,
                    ..PeakArea::default()
                };
                for peak in &chromatogram.peaks[begin..end] {
                    if apex.height < f64::from(peak.intensity) {
                        apex.height = f64::from(peak.intensity);
                        apex.apex_pos = peak.rt;
                    }
                }
                if apex.apex_pos > 1e-11 && (apex.apex_pos - peak_apex).abs() < distance_min {
                    distance_min = (apex.apex_pos - peak_apex).abs();
                    closest = Some(i);
                }
            }
            let (mut l, mut r) = (best_left, best_right);
            if let Some(i) = closest {
                l = array_value(chromatogram, PeakPickerChromatogram::IDX_LEFTBORDER, i)?;
                r = array_value(chromatogram, PeakPickerChromatogram::IDX_RIGHTBORDER, i)?;
                chromatogram.peaks[i].intensity = 0.0;
            }
            left_edges.push(l);
            right_edges.push(r);
            if l < *min_left {
                *min_left = l;
            }
            if r > *max_right {
                *max_right = r;
            }
        }
        Ok(())
    }

    /// Resample, integrate and annotate every transition's chromatogram for
    /// one peak group (source `pickFragmentChromatograms`).
    ///
    /// For each transition, in transition order: its raw chromatogram (or,
    /// with [`PeakIntegrationSignal::Smoothed`], `smoothed[k]`) is resampled
    /// onto `master` around the boundaries and integrated between them; the
    /// background is removed as [`background_subtraction`](Self::background_subtraction)
    /// says, clamping the area and apex at zero; and a sub-feature keyed by
    /// the chromatogram's native ID is added to `feature` with the seed's
    /// retention time, the product m/z (or the legacy `product_mz` metadata,
    /// also copied to `MZ`), `native_id`, `peak_apex_int`,
    /// `peak_apex_position`, `total_xic`, `width_at_50`, and the optional
    /// metrics. `feature` takes the precursor m/z. Quantifying transitions add
    /// to `totals.total_intensity` and `totals.total_peak_apices`; the
    /// detecting transitions' summed intensity is added to `totals.total_xic`.
    ///
    /// Outside consensus mode each transition uses its own entry of
    /// `left_edges` and `right_edges`, and its hull is the whole resampled
    /// chromatogram so that every transition yields equally long intensity
    /// vectors for cross-correlation scoring (source issue 9138).
    ///
    /// Two source behaviours are kept. `smoothed` is indexed by transition,
    /// though it holds only the picked (detecting) chromatograms in
    /// chromatogram order, so smoothed integration misaligns when a
    /// non-detecting transition precedes a detecting one. And with
    /// [`compute_total_mi`](Self::compute_total_mi), a detecting transition's
    /// mutual information, already averaged over the detecting transitions, is
    /// divided by their number again before it is added to
    /// `totals.total_mi`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] outside consensus mode for a non-detecting
    /// transition (source `Exception::IllegalArgument`) or a missing edge; for
    /// smoothed integration without `smoothed[k]` (source
    /// `Exception::IllegalArgument`); when the `original` background finds no
    /// raw sample between the boundaries (the source dereferences an end
    /// iterator); when `product_mz` is not numeric; and when a value is not
    /// finite. Also the errors of the [`integrator`](Self::integrator) and of
    /// [`ranked_mutual_information`] (chromatograms of different lengths).
    #[allow(clippy::too_many_arguments)]
    pub fn pick_fragment_chromatograms<T: Transition>(
        &self,
        group: &MRMTransitionGroup<MSChromatogram, T>,
        picked: &[MSChromatogram],
        feature: &mut MRMFeature,
        smoothed: &[MSChromatogram],
        best_left: f64,
        best_right: f64,
        use_consensus: bool,
        totals: &mut PeakGroupTotals,
        master: &MSChromatogram,
        left_edges: &[f64],
        right_edges: &[f64],
        inputs: &TransitionGroupInputs<'_>,
        chr_idx: usize,
        peak_idx: usize,
    ) -> Result<()> {
        let seed_rt = seed_position(picked, chr_idx, peak_idx)?;
        totals.total_xic += inputs.detecting_total_xic;
        let transitions = group.transitions();
        for (k, transition) in transitions.iter().enumerate() {
            let (mut local_left, mut local_right) = (best_left, best_right);
            if !use_consensus {
                if !transition.is_detecting() {
                    return Err(invalid(
                        "When using non-consensus peak picker, all transitions need to be detecting transitions.",
                    ));
                }
                local_left = *left_edges
                    .get(k)
                    .ok_or_else(|| invalid("missing left edge for a transition"))?;
                local_right = *right_edges
                    .get(k)
                    .ok_or_else(|| invalid("missing right edge for a transition"))?;
            }
            let chromatogram = *inputs
                .transition_chromatograms
                .get(k)
                .ok_or_else(|| invalid("missing chromatogram for a transition"))?;
            let transition_total_xic = *inputs
                .transition_total_xics
                .get(k)
                .ok_or_else(|| invalid("missing total XIC for a transition"))?;

            let mut transition_total_mi = 0.0;
            if self.compute_total_mi {
                let intensities: Vec<f64> = chromatogram
                    .peaks
                    .iter()
                    .map(|p| f64::from(p.intensity))
                    .collect();
                let mut det_ranked = Vec::new();
                let max_rank_det = compute_and_append_rank(&intensities, &mut det_ranked)?;
                let mut norm = 0_i32;
                let mut id_ranked = Vec::new();
                for (m, other) in transitions.iter().enumerate() {
                    if other.is_detecting() {
                        let detecting: Vec<f64> = inputs
                            .transition_chromatograms
                            .get(m)
                            .ok_or_else(|| invalid("missing chromatogram for a transition"))?
                            .peaks
                            .iter()
                            .map(|p| f64::from(p.intensity))
                            .collect();
                        let max_rank_id = compute_and_append_rank(&detecting, &mut id_ranked)?;
                        transition_total_mi += ranked_mutual_information(
                            &id_ranked,
                            &det_ranked,
                            max_rank_id,
                            max_rank_det,
                        )?;
                        norm += 1;
                    }
                }
                if norm > 0 {
                    transition_total_mi /= f64::from(norm);
                }
                if transition.is_detecting() {
                    totals.total_mi += transition_total_mi / f64::from(norm);
                }
            }

            let source = match self.peak_integration {
                PeakIntegrationSignal::Original => chromatogram,
                PeakIntegrationSignal::Smoothed => smoothed.get(k).ok_or_else(|| {
                    invalid(
                        "Tried to calculate peak area and height without any smoothed chromatograms",
                    )
                })?,
            };
            let used = Self::resample_chromatogram(source, master, local_left, local_right);

            // Source `setQuality(0, 0)` and `setOverallQuality(0)`: the defaults.
            let mut sub = Feature::default();
            let mut area =
                self.integrator
                    .integrate_chromatogram(&used, local_left, local_right)?;
            let mut peak_integral = area.area;
            let mut peak_apex_int = area.height;
            set_float(&mut sub, meta::PEAK_APEX_POSITION, area.apex_pos)?;
            if !use_consensus {
                area.hull_points = used
                    .peaks
                    .iter()
                    .map(|p| [p.rt, f64::from(p.intensity)])
                    .collect();
            }
            if self.background_subtraction != BackgroundSubtraction::None {
                let (background, noise) =
                    self.background(chromatogram, &used, local_left, local_right, area.apex_pos)?;
                peak_integral -= background;
                peak_apex_int -= noise;
                if peak_integral < 0.0 {
                    peak_integral = 0.0;
                }
                if peak_apex_int < 0.0 {
                    peak_apex_int = 0.0;
                }
                set_float(&mut sub, meta::AREA_BACKGROUND_LEVEL, background)?;
                set_float(&mut sub, meta::NOISE_BACKGROUND_LEVEL, noise)?;
            }
            sub.base.rt = seed_rt;
            sub.base.intensity = peak_integral as f32;
            sub.convex_hulls.push(hull(&area.hull_points)?);
            sub.base.mz = chromatogram.product.mz;
            feature.base.mz = chromatogram.precursor.mz;
            if let Some(product_mz) = chromatogram.metadata.get("product_mz") {
                sub.base.mz = product_mz.as_f64()?;
                set_value(&mut sub, "MZ", product_mz.clone());
            }
            set_value(
                &mut sub,
                meta::NATIVE_ID,
                MetaValue::from(chromatogram.native_id.as_str()),
            );
            set_float(&mut sub, meta::PEAK_APEX_INT, peak_apex_int)?;
            set_float(&mut sub, meta::TOTAL_XIC, transition_total_xic)?;
            if self.compute_total_mi {
                set_float(&mut sub, meta::TOTAL_MI, transition_total_mi)?;
            }
            if transition.is_quantifying() {
                totals.total_intensity += peak_integral;
                totals.total_peak_apices += peak_apex_int;
            }
            let shape = self.integrator.calculate_shape_metrics_chromatogram(
                &used,
                local_left,
                local_right,
                peak_apex_int,
                area.apex_pos,
            )?;
            set_float(&mut sub, meta::WIDTH_AT_50, shape.width_at_50)?;
            if self.compute_peak_shape_metrics {
                for (key, value) in [
                    (meta::WIDTH_AT_5, shape.width_at_5),
                    (meta::WIDTH_AT_10, shape.width_at_10),
                    (meta::START_POSITION_AT_5, shape.start_position_at_5),
                    (meta::START_POSITION_AT_10, shape.start_position_at_10),
                    (meta::START_POSITION_AT_50, shape.start_position_at_50),
                    (meta::END_POSITION_AT_5, shape.end_position_at_5),
                    (meta::END_POSITION_AT_10, shape.end_position_at_10),
                    (meta::END_POSITION_AT_50, shape.end_position_at_50),
                    (meta::TOTAL_WIDTH, shape.total_width),
                    (meta::TAILING_FACTOR, shape.tailing_factor),
                    (meta::ASYMMETRY_FACTOR, shape.asymmetry_factor),
                    (meta::SLOPE_OF_BASELINE, shape.slope_of_baseline),
                    (meta::BASELINE_DELTA_2_HEIGHT, shape.baseline_delta_2_height),
                ] {
                    set_float(&mut sub, key, value)?;
                }
                set_value(
                    &mut sub,
                    meta::POINTS_ACROSS_BASELINE,
                    count(shape.points_across_baseline),
                );
                set_value(
                    &mut sub,
                    meta::POINTS_ACROSS_HALF_HEIGHT,
                    count(shape.points_across_half_height),
                );
            }
            feature.add_feature_with(
                sub,
                &chromatogram.native_id,
                DuplicateKeyPolicy::SourceOverwrite,
            )?;
        }
        Ok(())
    }

    /// Resample, integrate and annotate the precursor chromatograms for one
    /// peak group (source `pickPrecursorChromatograms`).
    ///
    /// As [`MRMTransitionGroupPicker::pick_fragment_chromatograms`], for each
    /// precursor chromatogram, without shape metrics, `total_xic` or
    /// `total_mi`. The sub-feature takes the precursor m/z (or the legacy
    /// `precursor_mz` metadata), which the peak group also takes from the
    /// first (monoisotopic) precursor. A precursor's integral adds to
    /// `total_intensity` only when [`use_precursors`](Self::use_precursors) is
    /// set and the group has no transitions.
    ///
    /// The precursor's smoothed chromatogram and edges are looked up at index
    /// "number of fragment chromatograms + k", which is only right when every
    /// fragment chromatogram was picked (all transitions detecting) and
    /// [`use_precursors`](Self::use_precursors) is set; otherwise smoothed
    /// integration fails as in the source, and outside consensus mode the
    /// seed's boundaries are used when the edge is missing.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for smoothed integration without that smoothed
    /// chromatogram (source `Exception::IllegalArgument`), for an `original`
    /// background with no raw sample between the boundaries, when
    /// `precursor_mz` is not numeric, and when a value is not finite; the
    /// errors of the [`integrator`](Self::integrator).
    #[allow(clippy::too_many_arguments)]
    pub fn pick_precursor_chromatograms<T: Transition>(
        &self,
        group: &MRMTransitionGroup<MSChromatogram, T>,
        picked: &[MSChromatogram],
        feature: &mut MRMFeature,
        smoothed: &[MSChromatogram],
        best_left: f64,
        best_right: f64,
        use_consensus: bool,
        total_intensity: &mut f64,
        master: &MSChromatogram,
        left_edges: &[f64],
        right_edges: &[f64],
        chr_idx: usize,
        peak_idx: usize,
    ) -> Result<()> {
        let seed_rt = seed_position(picked, chr_idx, peak_idx)?;
        for (k, chromatogram) in group.precursor_chromatograms().iter().enumerate() {
            let prec_idx = group.chromatograms().len() + k;
            let (mut local_left, mut local_right) = (best_left, best_right);
            if !use_consensus && right_edges.len() > prec_idx && left_edges.len() > prec_idx {
                local_left = left_edges[prec_idx];
                local_right = right_edges[prec_idx];
            }
            let source = match self.peak_integration {
                PeakIntegrationSignal::Original => chromatogram,
                PeakIntegrationSignal::Smoothed => smoothed.get(prec_idx).ok_or_else(|| {
                    invalid(
                        "Tried to calculate peak area and height without any smoothed chromatograms for precursors",
                    )
                })?,
            };
            let used = Self::resample_chromatogram(source, master, local_left, local_right);
            // Source `setQuality(0, 0)` and `setOverallQuality(0)`: the defaults.
            let mut sub = Feature::default();
            let mut area =
                self.integrator
                    .integrate_chromatogram(&used, local_left, local_right)?;
            let mut peak_integral = area.area;
            let mut peak_apex_int = area.height;
            if !use_consensus {
                area.hull_points = used
                    .peaks
                    .iter()
                    .map(|p| [p.rt, f64::from(p.intensity)])
                    .collect();
            }
            if self.background_subtraction != BackgroundSubtraction::None {
                let (background, noise) =
                    self.background(chromatogram, &used, local_left, local_right, area.apex_pos)?;
                peak_integral -= background;
                peak_apex_int -= noise;
                if peak_integral < 0.0 {
                    peak_integral = 0.0;
                }
                if peak_apex_int < 0.0 {
                    peak_apex_int = 0.0;
                }
                set_float(&mut sub, meta::AREA_BACKGROUND_LEVEL, background)?;
                set_float(&mut sub, meta::NOISE_BACKGROUND_LEVEL, noise)?;
            }
            sub.base.mz = chromatogram.precursor.mz;
            if k == 0 {
                feature.base.mz = chromatogram.precursor.mz;
            }
            if let Some(precursor_mz) = chromatogram.metadata.get("precursor_mz") {
                let mz = precursor_mz.as_f64()?;
                sub.base.mz = mz;
                if k == 0 {
                    feature.base.mz = mz;
                }
            }
            sub.base.rt = seed_rt;
            sub.base.intensity = peak_integral as f32;
            sub.convex_hulls.push(hull(&area.hull_points)?);
            set_value(
                &mut sub,
                meta::NATIVE_ID,
                MetaValue::from(chromatogram.native_id.as_str()),
            );
            set_float(&mut sub, meta::PEAK_APEX_INT, peak_apex_int)?;
            if self.use_precursors && group.transitions().is_empty() {
                *total_intensity += peak_integral;
            }
            feature.add_precursor_feature_with(
                sub,
                &chromatogram.native_id,
                DuplicateKeyPolicy::SourceOverwrite,
            )?;
        }
        Ok(())
    }

    /// The background area and apex level of one chromatogram, as
    /// [`background_subtraction`](Self::background_subtraction) selects.
    fn background(
        &self,
        raw: &MSChromatogram,
        used: &MSChromatogram,
        left: f64,
        right: f64,
        apex: f64,
    ) -> Result<(f64, f64)> {
        match self.background_subtraction {
            BackgroundSubtraction::None => Ok((0.0, 0.0)),
            BackgroundSubtraction::Original => {
                // `PosBegin(left)` and `PosEnd(right) - 1` on the raw trace.
                let begin = raw.peaks.partition_point(|p| p.rt < left);
                let end = raw.peaks.partition_point(|p| p.rt <= right);
                if begin >= raw.peaks.len() || end == 0 || end < begin {
                    return Err(invalid(
                        "original background subtraction found no raw sample between the boundaries",
                    ));
                }
                let intensity_left = f64::from(raw.peaks[begin].intensity);
                let intensity_right = f64::from(raw.peaks[end - 1].intensity);
                let points = (end - begin) as f64;
                let noise = (intensity_right + intensity_left) / 2.0;
                Ok((noise * points, noise))
            }
            BackgroundSubtraction::Exact => {
                let background = self
                    .integrator
                    .estimate_background_chromatogram(used, left, right, apex)?;
                Ok((background.area, background.height))
            }
        }
    }

    /// Consume every picked peak inside or overlapping the current seed's
    /// boundaries (source `remove_overlapping_features`).
    ///
    /// A peak is inside when its position lies in `[best_left, best_right]`;
    /// it overlaps when either of its own borders lies strictly inside
    /// `(best_left, best_right)`. Directly adjacent peaks, which share a
    /// border, are kept. Consumed peaks get intensity zero; peaks already at or
    /// below `1e-11` are not checked for overlap.
    ///
    /// Native: a peak whose boundary arrays are missing is only checked for
    /// being inside, where the source would read out of bounds.
    pub fn remove_overlapping_features(
        picked: &mut [MSChromatogram],
        best_left: f64,
        best_right: f64,
    ) {
        for chromatogram in picked.iter_mut() {
            for peak in &mut chromatogram.peaks {
                if peak.rt >= best_left && peak.rt <= best_right {
                    peak.intensity = 0.0;
                }
            }
        }
        for chromatogram in picked.iter_mut() {
            for i in 0..chromatogram.peaks.len() {
                if f64::from(chromatogram.peaks[i].intensity) <= 1e-11 {
                    continue;
                }
                let (Ok(left), Ok(right)) = (
                    array_value(chromatogram, PeakPickerChromatogram::IDX_LEFTBORDER, i),
                    array_value(chromatogram, PeakPickerChromatogram::IDX_RIGHTBORDER, i),
                ) else {
                    continue;
                };
                if (left > best_left && left < best_right)
                    || (right > best_left && right < best_right)
                {
                    chromatogram.peaks[i].intensity = 0.0;
                }
            }
        }
    }

    /// The most intense peak over all picked chromatograms, as
    /// `(chromatogram, peak)` (source `findLargestPeak`); `None` when no peak
    /// has a positive intensity. Ties keep the first chromatogram and the
    /// first peak.
    pub fn find_largest_peak(picked: &[MSChromatogram]) -> Option<(usize, usize)> {
        let mut largest = 0.0_f32;
        let mut best = None;
        for (k, chromatogram) in picked.iter().enumerate() {
            for (i, peak) in chromatogram.peaks.iter().enumerate() {
                if peak.intensity > largest {
                    largest = peak.intensity;
                    best = Some((k, i));
                }
            }
        }
        best
    }

    /// [`MRMTransitionGroupPicker::find_largest_peak`] for chromatograms
    /// sorted by ascending intensity, with a per-chromatogram cursor (source
    /// `findLargestPeakCached_`).
    ///
    /// `cursors[k]` starts at the last peak of `picked[k]` and walks down past
    /// peaks whose intensity is zero or less; the positive candidates are then
    /// compared, the first chromatogram winning ties. Within a chromatogram
    /// the highest index wins, which differs from
    /// [`MRMTransitionGroupPicker::find_largest_peak`] only for equal
    /// intensities.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `cursors` and `picked` differ in length
    /// (a source precondition) or a cursor is out of range.
    pub fn find_largest_peak_cached(
        picked: &[MSChromatogram],
        cursors: &mut [Option<usize>],
    ) -> Result<Option<(usize, usize)>> {
        if cursors.len() != picked.len() {
            return Err(invalid(
                "Cached largest peak indices must match picked chromatogram count.",
            ));
        }
        let mut largest = 0.0_f64;
        let mut best = None;
        for (k, chromatogram) in picked.iter().enumerate() {
            let mut candidate = cursors[k];
            while let Some(index) = candidate {
                let peak = chromatogram
                    .peaks
                    .get(index)
                    .ok_or_else(|| invalid("cached peak index out of range"))?;
                if peak.intensity > 0.0 {
                    break;
                }
                candidate = index.checked_sub(1);
            }
            cursors[k] = candidate;
            let Some(index) = candidate else {
                continue;
            };
            let intensity = f64::from(chromatogram.peaks[index].intensity);
            if intensity > largest {
                largest = intensity;
                best = Some((k, index));
            }
        }
        Ok(best)
    }

    /// The picked peak with the widest boundaries over all chromatograms, as
    /// `(chromatogram, peak)` (source `findWidestPeakIndices`), skipping
    /// consumed peaks (intensity exactly zero). `None` when no remaining peak
    /// has a positive width. Ties keep the first.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a remaining peak has no boundary values
    /// (the source reads out of bounds).
    pub fn find_widest_peak_indices(picked: &[MSChromatogram]) -> Result<Option<(usize, usize)>> {
        let mut max_width = 0.0;
        let mut best = None;
        for (i, chromatogram) in picked.iter().enumerate() {
            for k in 0..chromatogram.peaks.len() {
                if chromatogram.peaks[k].intensity == 0.0 {
                    continue;
                }
                let left = array_value(chromatogram, PeakPickerChromatogram::IDX_LEFTBORDER, k)?;
                let right = array_value(chromatogram, PeakPickerChromatogram::IDX_RIGHTBORDER, k)?;
                let width = right - left;
                if width > max_width {
                    max_width = width;
                    best = Some((i, k));
                }
            }
        }
        Ok(best)
    }

    /// Transition-group quality, higher is better (source `computeQuality_`).
    ///
    /// Only the co-elution and internal consistency of the chromatograms are
    /// used, no library information: every picked chromatogram's raw trace is
    /// resampled onto the seed's sampling, [`resample_boundary`](Self::resample_boundary)
    /// seconds beyond the boundaries on either side; traces without signal are
    /// dropped; and each remaining trace is cross-correlated with every other.
    /// The score is the mean shape score, minus `(mean co-elution - 1) / 2`,
    /// minus the fraction of chromatograms without a picked peak inside the
    /// boundaries. When one trace has both the worst shape and the worst
    /// co-elution, its native ID is the potential outlier, else `none`.
    ///
    /// Three source behaviours are kept. The cross-correlation standardises
    /// the traces in place, so later correlations see already standardised
    /// data. The co-elution mean accumulates into the accumulator that already
    /// holds the shape scores, so it averages shapes and co-elutions together.
    /// And the outlier index counts only traces with signal but indexes all
    /// picked chromatograms. (The source's per-chromatogram border collection
    /// stores a comparison, `max_int = intensity > max_int`, where it means an
    /// intensity; the borders it selects are never used, so only the peak
    /// counts are computed here.)
    fn compute_quality(
        &self,
        picked: &[MSChromatogram],
        inputs: &TransitionGroupInputs<'_>,
        chr_idx: usize,
        best_left: f64,
        best_right: f64,
    ) -> Result<(f64, String)> {
        let boundary = self.resample_boundary;
        let reference = inputs
            .picked_input_chromatograms
            .get(chr_idx)
            .ok_or_else(|| invalid("no raw chromatogram for the seed chromatogram"))?;
        let master =
            Self::prepare_master_container(reference, best_left - boundary, best_right + boundary)?;
        let mut all_ints: Vec<Vec<f64>> = Vec::new();
        for k in 0..picked.len() {
            let raw = inputs
                .picked_input_chromatograms
                .get(k)
                .ok_or_else(|| invalid("no raw chromatogram for a picked chromatogram"))?;
            let used = Self::resample_chromatogram(
                raw,
                &master,
                best_left - boundary,
                best_right + boundary,
            );
            let intensities: Vec<f64> = used.peaks.iter().map(|p| f64::from(p.intensity)).collect();
            let tic = intensities.iter().fold(0.0, |sum, value| sum + value);
            if tic > 1e-11 {
                all_ints.push(intensities);
            }
        }

        let mut mean_shapes = Vec::with_capacity(all_ints.len());
        let mut mean_coel = Vec::with_capacity(all_ints.len());
        for k in 0..all_ints.len() {
            let mut shapes = Vec::new();
            let mut coel = Vec::new();
            for i in 0..all_ints.len() {
                if i == k {
                    continue;
                }
                let maxdelay = i32::try_from(all_ints[i].len())
                    .map_err(|_| invalid("trace too long for cross-correlation"))?;
                let (first, second) = if k < i {
                    let (head, tail) = all_ints.split_at_mut(i);
                    (&mut head[k], &mut tail[0])
                } else {
                    let (head, tail) = all_ints.split_at_mut(k);
                    (&mut tail[0], &mut head[i])
                };
                let result = normalized_cross_correlation(first, second, maxdelay, 1)?;
                let best = result.data[xcorr_array_get_max_peak(&result)?];
                coel.push(f64::from(best.0).abs());
                shapes.push(best.1.abs());
            }
            let mut accumulator = MeanAndStddev::new();
            accumulator.extend(shapes.iter().copied());
            let shapes_mean = accumulator.mean();
            accumulator.extend(coel.iter().copied());
            let coel_mean = accumulator.mean();
            mean_shapes.push(shapes_mean);
            mean_coel.push(coel_mean);
        }

        let mut min_index_shape = 0;
        for (i, value) in mean_shapes.iter().enumerate() {
            if *value < mean_shapes[min_index_shape] {
                min_index_shape = i;
            }
        }
        let mut max_index_coel = 0;
        for (i, value) in mean_coel.iter().enumerate() {
            if mean_coel[max_index_coel] < *value {
                max_index_coel = i;
            }
        }

        let mut missing_peaks = 0_u32;
        for chromatogram in picked {
            let found = chromatogram
                .peaks
                .iter()
                .filter(|p| p.rt >= best_left && p.rt <= best_right)
                .count();
            if found == 0 {
                missing_peaks += 1;
            }
        }

        let outlier = if min_index_shape == max_index_coel {
            picked
                .get(min_index_shape)
                .map(|c| c.native_id.clone())
                .unwrap_or_default()
        } else {
            "none".to_owned()
        };
        let shape_score =
            mean_shapes.iter().fold(0.0, |sum, value| sum + value) / mean_shapes.len() as f64;
        let mut coel_score =
            mean_coel.iter().fold(0.0, |sum, value| sum + value) / mean_coel.len() as f64;
        coel_score = (coel_score - 1.0) / 2.0;
        let score = shape_score - coel_score - 1.0 * f64::from(missing_peaks) / picked.len() as f64;
        if !score.is_finite() {
            return Err(invalid(
                "peak quality is not finite (no chromatogram has signal around the peak)",
            ));
        }
        Ok((score, outlier))
    }

    /// Recalculate the peak borders from the other chromatograms' peaks
    /// (source `recalculatePeakBorders_`).
    ///
    /// The borders of the most abundant (`IntegratedIntensity`) picked peak
    /// inside `[best_left, best_right]` of each chromatogram are collected;
    /// when a border lies more than `max_z` population standard deviations
    /// from their mean, it is replaced by their upper median ("pseudo
    /// median"). The standard deviation is `sqrt(E[x^2] - mean^2)` as the
    /// source computes it, so equal borders give NaN or infinite z-scores,
    /// which compare as the Release build compares them.
    fn recalculate_peak_borders(
        picked: &[MSChromatogram],
        best_left: &mut f64,
        best_right: &mut f64,
        max_z: f64,
    ) -> Result<()> {
        let mut left_borders = Vec::new();
        let mut right_borders = Vec::new();
        for chromatogram in picked {
            let mut max_int = -1.0;
            let (mut left, mut right) = (-1.0, -1.0);
            for (i, peak) in chromatogram.peaks.iter().enumerate() {
                if peak.rt >= *best_left && peak.rt <= *best_right {
                    let abundance =
                        array_value(chromatogram, PeakPickerChromatogram::IDX_ABUNDANCE, i)?;
                    if abundance > max_int {
                        max_int = abundance;
                        left =
                            array_value(chromatogram, PeakPickerChromatogram::IDX_LEFTBORDER, i)?;
                        right =
                            array_value(chromatogram, PeakPickerChromatogram::IDX_RIGHTBORDER, i)?;
                    }
                }
            }
            if max_int > -1.0 {
                left_borders.push(left);
                right_borders.push(right);
            }
        }
        if right_borders.is_empty() {
            return Ok(());
        }
        let z_replace = |borders: &mut Vec<f64>, best: &mut f64| {
            let n = borders.len() as f64;
            let mean = borders.iter().fold(0.0, |sum, value| sum + value) / n;
            let squares = borders.iter().fold(0.0, |sum, value| sum + value * value);
            let stdev = (squares / n - mean * mean).sqrt();
            borders.sort_by(f64::total_cmp);
            if (*best - mean).abs() / stdev > max_z {
                *best = borders[borders.len() / 2];
            }
        };
        z_replace(&mut right_borders, best_right);
        z_replace(&mut left_borders, best_left);
        Ok(())
    }

    /// An empty container with the reference's positions between the
    /// boundaries, plus one sample beyond each (source
    /// `prepareMasterContainer_`); all intensities are zero.
    fn prepare_master_container(
        reference: &MSChromatogram,
        left: f64,
        right: f64,
    ) -> Result<MSChromatogram> {
        let (begin, end) = padded_range(&reference.peaks, left, right);
        if begin >= end {
            return Err(invalid(
                "the seed's raw chromatogram has no samples to resample onto",
            ));
        }
        Ok(MSChromatogram {
            peaks: reference.peaks[begin..end]
                .iter()
                .map(|p| ChromatogramPeak::new(p.rt, 0.0))
                .collect(),
            ..MSChromatogram::default()
        })
    }

    /// `chromatogram` resampled at the positions of `master`, using the
    /// samples between the boundaries plus one beyond each (source
    /// `resampleChromatogram_`); see [`raster`].
    fn resample_chromatogram(
        chromatogram: &MSChromatogram,
        master: &MSChromatogram,
        left: f64,
        right: f64,
    ) -> MSChromatogram {
        let (begin, end) = padded_range(&chromatogram.peaks, left, right);
        let mut resampled = master.clone();
        raster(&chromatogram.peaks[begin..end], &mut resampled.peaks);
        resampled
    }
}

fn seed_position(picked: &[MSChromatogram], chr_idx: usize, peak_idx: usize) -> Result<f64> {
    picked
        .get(chr_idx)
        .and_then(|chromatogram| chromatogram.peaks.get(peak_idx))
        .map(|peak| peak.rt)
        .ok_or_else(|| invalid("seed peak index out of range"))
}

fn hull(points: &[[f64; 2]]) -> Result<ConvexHull2D> {
    let points: Vec<Point2D> = points.iter().map(|&[x, y]| Point2D::new(x, y)).collect();
    let mut hull = ConvexHull2D::new();
    hull.set_hull_points(&points)?;
    Ok(hull)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn trace(points: &[(f64, f32)]) -> Vec<ChromatogramPeak> {
        points
            .iter()
            .map(|&(rt, intensity)| ChromatogramPeak::new(rt, intensity))
            .collect()
    }

    #[test]
    fn source_raster_conserves_intensity_and_splits_by_distance() {
        let raw = trace(&[(0.5, 4.0), (1.0, 2.0), (1.25, 8.0), (3.5, 1.0)]);
        let mut out = trace(&[(1.0, 0.0), (2.0, 0.0), (3.0, 0.0)]);
        raster(&raw, &mut out);
        // 0.5 lies before the grid, 1.0 on its first point, 1.25 a quarter of
        // the way to 2.0, and 3.5 beyond the last point.
        assert_eq!(out[0].intensity, 4.0 + 2.0 + 6.0);
        assert_eq!(out[1].intensity, 2.0);
        assert_eq!(out[2].intensity, 1.0);
    }

    #[test]
    fn libstdcxx_bounds_on_unsorted_positions() {
        // Positions 30, 10, 20: lower_bound(15) probes the middle (10 < 15),
        // moves right, probes 20 and stops there.
        let peaks = trace(&[(30.0, 1.0), (10.0, 2.0), (20.0, 3.0)]);
        assert_eq!(lower_bound(&peaks, 15.0), 2);
        assert_eq!(upper_bound(&peaks, 15.0), 2);
        assert_eq!(lower_bound(&peaks, 35.0), 3);
        assert_eq!(upper_bound(&peaks, 5.0), 0);
    }

    #[test]
    fn padded_range_adds_one_neighbour_on_each_side() {
        let peaks = trace(&[(1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0), (5.0, 0.0)]);
        assert_eq!(padded_range(&peaks, 2.5, 3.5), (1, 4));
        assert_eq!(padded_range(&peaks, 0.0, 10.0), (0, 5));
        assert_eq!(padded_range(&peaks, 3.0, 3.0), (1, 3));
    }
}
