// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The class that calls the OpenSWATH scoring routines, from
//! `ANALYSIS/OPENSWATH/OpenSwathScoring.h` (class `OpenSwathScoring`).
//!
//! [`OpenSwathScoring`](crate::analysis::openswath::openswath_scoring::OpenSwathScoring) invokes the chromatographic
//! ([`MRMScoring`](crate::analysis::openswath::mrm_scoring::MRMScoring)), library, full-spectrum DIA
//! ([`DiaScoring`](crate::analysis::openswath::dia_scoring::DiaScoring)) and ion mobility
//! ([`IonMobilityScoring`](crate::analysis::openswath::ion_mobility_scoring::IonMobilityScoring)) scores for one
//! peak group, as selected by an [`OpenSwathScoresUsage`](crate::analysis::openswath::openswath_scores::OpenSwathScoresUsage), and fetches the
//! DIA spectra at the peak apex ([`OpenSwathScoring::fetch_spectrum_swath`](crate::analysis::openswath::openswath_scoring::OpenSwathScoring::fetch_spectrum_swath)).
//!
//! Every scoring call writes into its [`OpenSwathScores`](crate::kernel::mrm::OpenSwathScores) (or
//! [`OpenSwathIndScores`](crate::kernel::mrm::OpenSwathIndScores)) only on success; on an error the scores are
//! unchanged.
//!
//! The source keeps thread-local pools and caches (spectra, formulas,
//! sequences, a persistent `MRMScoring`); the port recomputes, which yields the
//! same values. The source logs debug messages; the port does not.
//!
//! # Gaps
//!
//! Ion mobility peak picking (`apply_im_peak_picking`) and the Parquet
//! mobilogram export are gaps of
//! [`ion_mobility_scoring`](crate::analysis::openswath::ion_mobility_scoring),
//! which this module calls.
//!
//! See `docs/OPENSWATH_SCORING_SUPPORT.md`.

use super::data_access_helper::convert_peptide_to_aa_sequence;
use super::data_structures::SpectrumPtr;
use super::dia_helper::range_bounds;
use super::dia_scoring::DiaScoring;
use super::ion_mobility_scoring::{
    DriftExtraction, IonMobilityScoring, MobilogramExport, MobilogramSink,
};
use super::itransition::{MRMFeatureAccess, SignalToNoisePtr};
use super::mrm_scoring::MRMScoring;
use super::openswath_scores::OpenSwathScoresUsage;
use super::spectrum_access::{
    SpectrumAccess, SpectrumAccessPtr, SpectrumSequence, filter_by_drift,
};
use super::spectrum_addition::{add_up_spectra, concatenate_spectra};
use super::swath_map::SwathMap;
use super::transition_experiment::{LightCompound, LightTransition};
use crate::chemistry::{EmpiricalFormula, PeptideFragmentType};
use crate::kernel::mrm::{OpenSwathIndScores, OpenSwathScores};
use crate::kernel::ranges::RangeBase;
use crate::{Error, Result};
use std::sync::Arc;

/// How the spectra around the apex are combined (source private enum
/// `SpectrumAdditionMethod`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum SpectrumAdditionMethod {
    /// Keep the spectra as a sequence; the scores integrate across all of them
    /// (`"simple"`).
    #[default]
    Addition,
    /// Merge them into one resampled (or concatenated, with ion mobility)
    /// spectrum (`"resample"`).
    Resample,
}

/// How many spectra are merged (source private enum `SpectrumMergeMethodType`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum SpectrumMergeMethodType {
    /// A fixed number, `add_up_spectra` (`"fixed"`).
    #[default]
    Fixed,
    /// A fraction of the chromatographic peak width (`"dynamic"`).
    Dynamic,
}

/// The OpenSWATH scoring driver (source `OpenSwathScoring`).
///
/// The source default constructor leaves `merge_spectra_by_peak_width_fraction_`,
/// `use_ms1_ion_mobility_` and `apply_im_peak_picking_` uninitialised; the port
/// defaults them to 0, `false` and `false`. Call
/// [`OpenSwathScoring::initialize`] to set everything, as the source's callers
/// do.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OpenSwathScoring {
    /// Range of the normalised retention-time space; the RT score is divided
    /// by it.
    pub rt_normalization_factor: f64,
    /// Spacing for spectrum resampling.
    pub spacing_for_spectra_resampling: f64,
    /// Fraction of the peak width that sets the number of spectra to merge in
    /// [`SpectrumMergeMethodType::Dynamic`] mode.
    pub merge_spectra_by_peak_width_fraction: f64,
    /// How many spectra to add up in [`SpectrumMergeMethodType::Fixed`] mode.
    pub add_up_spectra: i32,
    /// Spectrum addition method.
    pub spectra_addition_method: SpectrumAdditionMethod,
    /// Spectrum merge method.
    pub spectra_merge_method_type: SpectrumMergeMethodType,
    /// Extension of the ion-mobility window on each side, as a fraction.
    pub im_drift_extra_pcnt: f64,
    /// Which scores to compute.
    pub su: OpenSwathScoresUsage,
    /// Whether to extract across ion mobility in MS1.
    pub use_ms1_ion_mobility: bool,
    /// Whether to peak-pick the ion mobilograms.
    pub apply_im_peak_picking: bool,
}

impl Default for OpenSwathScoring {
    fn default() -> Self {
        Self {
            rt_normalization_factor: 1.0,
            spacing_for_spectra_resampling: 0.005,
            merge_spectra_by_peak_width_fraction: 0.0,
            add_up_spectra: 1,
            spectra_addition_method: SpectrumAdditionMethod::Addition,
            spectra_merge_method_type: SpectrumMergeMethodType::Fixed,
            im_drift_extra_pcnt: 0.0,
            su: OpenSwathScoresUsage::default(),
            use_ms1_ion_mobility: false,
            apply_im_peak_picking: false,
        }
    }
}

fn map_access(map: &SwathMap) -> Result<&SpectrumAccessPtr> {
    map.sptr
        .as_ref()
        .ok_or_else(|| Error::MissingInformation("SWATH map without spectrum access".into()))
}

fn fetch_count(n: i32) -> usize {
    usize::try_from(n).unwrap_or(0)
}

fn contains_drift_time_data(spectra: &[SpectrumPtr]) -> bool {
    spectra.iter().any(|s| s.drift_time_array().is_some())
}

/// Source `concatenateDriftFilteredSpectra_`.
fn concatenate_drift_filtered(
    spectra: &[SpectrumPtr],
    im_range: &RangeBase,
) -> Result<SpectrumPtr> {
    let (lo, hi) = range_bounds(im_range);
    let filtered = spectra
        .iter()
        .map(|s| Ok(Arc::new(filter_by_drift(s, lo, hi)?)))
        .collect::<Result<Vec<_>>>()?;
    concatenate_spectra(&filtered)
}

/// Merge spectra for scoring (source `mergeSpectraForScoring_`): with an IM
/// range, concatenate the drift-filtered spectra, which keeps the drift arrays
/// for IM-aware scoring; with drift data but no range, concatenate (a single
/// spectrum is returned as is); otherwise add up by resampling, dropping zero
/// intensities.
fn merge_spectra_for_scoring(
    spectra: &[SpectrumPtr],
    im_range: &RangeBase,
    sampling_rate: f64,
) -> Result<SpectrumPtr> {
    if !im_range.is_empty() {
        return concatenate_drift_filtered(spectra, im_range);
    }
    if contains_drift_time_data(spectra) {
        return if spectra.len() == 1 {
            Ok(Arc::clone(&spectra[0]))
        } else {
            concatenate_spectra(spectra)
        };
    }
    add_up_spectra(spectra, sampling_rate, true)
}

impl OpenSwathScoring {
    /// A driver with the source constructor defaults (see the struct notes).
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the scoring parameters (source `initialize`).
    ///
    /// `spectrum_addition_method` is `"simple"` or `"resample"`;
    /// `spectrum_merge_method_type` is `"fixed"` or `"dynamic"`. `drift_extra`
    /// widens the ion-mobility window beyond its bounds, in percent (as a
    /// fraction).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for any other method string, as the source
    /// throws `Exception::InvalidValue`. The driver is unchanged on error.
    #[allow(clippy::too_many_arguments)]
    pub fn initialize(
        &mut self,
        rt_normalization_factor: f64,
        add_up_spectra: i32,
        spacing_for_spectra_resampling: f64,
        merge_spectra_by_peak_width_fraction: f64,
        drift_extra: f64,
        su: &OpenSwathScoresUsage,
        spectrum_addition_method: &str,
        spectrum_merge_method_type: &str,
        use_ms1_ion_mobility: bool,
        apply_im_peak_picking: bool,
    ) -> Result<()> {
        let addition = match spectrum_addition_method {
            "simple" => SpectrumAdditionMethod::Addition,
            "resample" => SpectrumAdditionMethod::Resample,
            other => {
                return Err(Error::InvalidValue(format!(
                    "spectrum_addition_method must be simple or resample, got '{other}'"
                )));
            }
        };
        let merge = match spectrum_merge_method_type {
            "fixed" => SpectrumMergeMethodType::Fixed,
            "dynamic" => SpectrumMergeMethodType::Dynamic,
            other => {
                return Err(Error::InvalidValue(format!(
                    "spectrum_merge_method_type must be fixed or dynamic, got '{other}'"
                )));
            }
        };
        *self = Self {
            rt_normalization_factor,
            spacing_for_spectra_resampling,
            merge_spectra_by_peak_width_fraction,
            add_up_spectra,
            spectra_addition_method: addition,
            spectra_merge_method_type: merge,
            im_drift_extra_pcnt: drift_extra,
            su: *su,
            use_ms1_ion_mobility,
            apply_im_peak_picking,
        };
        Ok(())
    }

    /// Score a peak group using only chromatographic properties (coelution,
    /// signal to noise, mutual information), computing only the scores
    /// selected in [`OpenSwathScoring::su`] (source
    /// `calculateChromatographicScores`).
    ///
    /// `native_ids` give the canonical transition order;
    /// `normalized_library_intensity` weights the weighted scores. With
    /// precursor features, the MS1 cross-correlation and MI scores are added;
    /// the pure precursor scores need at least two precursor isotopes.
    ///
    /// # Errors
    ///
    /// The errors of [`MRMScoring`] (unknown IDs, trace lengths, weights not
    /// matching the transitions), and [`Error::InvalidValue`] when the feature
    /// count does not fit an `int` (the source's `boost::numeric_cast`).
    pub fn calculate_chromatographic_scores(
        &self,
        feature: &dyn MRMFeatureAccess,
        native_ids: &[String],
        precursor_ids: &[String],
        normalized_library_intensity: &[f64],
        signal_noise_estimators: &[SignalToNoisePtr],
        scores: &mut OpenSwathScores,
    ) -> Result<()> {
        let su = &self.su;
        let mut out = *scores;
        let mut mrm = MRMScoring::new();
        let has_precursors = !precursor_ids.is_empty();
        if su.use_coelution_score
            || su.use_shape_score
            || (has_precursors && su.use_ms1_correlation)
        {
            mrm.initialize_xcorr_matrix(feature, native_ids)?;
        }
        if su.use_coelution_score {
            out.xcorr_coelution_score = mrm.calc_xcorr_coelution_score();
            out.weighted_coelution_score =
                mrm.calc_xcorr_coelution_weighted_score(normalized_library_intensity)?;
        }
        if su.use_shape_score {
            out.xcorr_shape_score = mrm.calc_xcorr_shape_score();
            out.weighted_xcorr_shape =
                mrm.calc_xcorr_shape_weighted_score(normalized_library_intensity)?;
        }
        if has_precursors && su.use_ms1_correlation {
            if precursor_ids.len() > 1 {
                mrm.initialize_xcorr_precursor_matrix(feature, precursor_ids)?;
                out.ms1_xcorr_coelution_score = mrm.calc_xcorr_precursor_coelution_score();
                out.ms1_xcorr_shape_score = mrm.calc_xcorr_precursor_shape_score();
            }
            mrm.initialize_xcorr_precursor_contrast_matrix(feature, precursor_ids, native_ids)?;
            out.ms1_xcorr_coelution_contrast_score =
                mrm.calc_xcorr_precursor_contrast_coelution_score();
            out.ms1_xcorr_shape_contrast_score = mrm.calc_xcorr_precursor_contrast_shape_score();
            mrm.initialize_xcorr_precursor_combined_matrix(feature, precursor_ids, native_ids)?;
            out.ms1_xcorr_coelution_combined_score =
                mrm.calc_xcorr_precursor_combined_coelution_score();
            out.ms1_xcorr_shape_combined_score = mrm.calc_xcorr_precursor_combined_shape_score();
        }
        if su.use_nr_peaks_score {
            let n = i32::try_from(feature.len())
                .map_err(|_| Error::InvalidValue("feature count does not fit an int".into()))?;
            out.nr_peaks = f64::from(n);
        }
        if su.use_sn_score {
            out.sn_ratio = MRMScoring::calc_sn_score(feature, signal_noise_estimators);
            out.log_sn_score = if out.sn_ratio < 1.0 {
                0.0
            } else {
                out.sn_ratio.ln()
            };
        }
        if su.use_mi_score {
            mrm.initialize_mi_matrix(feature, native_ids)?;
            out.mi_score = mrm.calc_mi_score();
            out.weighted_mi_score = mrm.calc_mi_weighted_score(normalized_library_intensity)?;
        }
        if has_precursors && su.use_ms1_mi {
            if precursor_ids.len() > 1 {
                mrm.initialize_mi_precursor_matrix(feature, precursor_ids)?;
                out.ms1_mi_score = mrm.calc_mi_precursor_score();
            }
            mrm.initialize_mi_precursor_contrast_matrix(feature, precursor_ids, native_ids)?;
            out.ms1_mi_contrast_score = mrm.calc_mi_precursor_contrast_score();
            mrm.initialize_mi_precursor_combined_matrix(feature, precursor_ids, native_ids)?;
            out.ms1_mi_combined_score = mrm.calc_mi_precursor_combined_score();
        }
        *scores = out;
        Ok(())
    }

    /// Score identification transitions against the detection transitions of
    /// a peak group using only chromatographic properties (source
    /// `calculateChromatographicIdScores`).
    ///
    /// Per identification transition, the contrast coelution and shape scores
    /// against the detection transitions, the log S/N and the contrast mutual
    /// information are set, as selected in [`OpenSwathScoring::su`]. The
    /// contrast cross-correlation matrix is always built, as in the source.
    ///
    /// # Errors
    ///
    /// The errors of [`MRMScoring`].
    pub fn calculate_chromatographic_id_scores(
        &self,
        feature: &dyn MRMFeatureAccess,
        native_ids_identification: &[String],
        native_ids_detection: &[String],
        signal_noise_estimators: &[SignalToNoisePtr],
        idscores: &mut OpenSwathIndScores,
    ) -> Result<()> {
        let mut mrm = MRMScoring::new();
        mrm.initialize_xcorr_contrast_matrix(
            feature,
            native_ids_identification,
            native_ids_detection,
        )?;
        let coelution = self
            .su
            .use_coelution_score
            .then(|| mrm.calc_separate_xcorr_contrast_coelution_score());
        let shape = self
            .su
            .use_shape_score
            .then(|| mrm.calc_separate_xcorr_contrast_shape_score());
        let sn = self
            .su
            .use_sn_score
            .then(|| MRMScoring::calc_separate_sn_score(feature, signal_noise_estimators));
        let mi = if self.su.use_mi_score {
            mrm.initialize_mi_contrast_matrix(
                feature,
                native_ids_identification,
                native_ids_detection,
            )?;
            Some(mrm.calc_separate_mi_contrast_score())
        } else {
            None
        };
        if let Some(v) = coelution {
            idscores.ind_xcorr_coelution_score = v;
        }
        if let Some(v) = shape {
            idscores.ind_xcorr_shape_score = v;
        }
        if let Some(v) = sn {
            idscores.ind_log_sn_score = v;
        }
        if let Some(v) = mi {
            idscores.ind_mi_score = v;
        }
        Ok(())
    }

    /// Score a chromatographic feature against a spectral library: the library
    /// scores of [`MRMScoring::calc_library_score`] and the retention-time
    /// score, the latter divided by
    /// [`OpenSwathScoring::rt_normalization_factor`] into `norm_rt_score`
    /// (source `calculateLibraryScores`).
    ///
    /// # Errors
    ///
    /// The errors of [`MRMScoring::calc_library_score`].
    pub fn calculate_library_scores(
        &self,
        feature: &dyn MRMFeatureAccess,
        transitions: &[LightTransition],
        compound: &LightCompound,
        normalized_feature_rt: f64,
        scores: &mut OpenSwathScores,
    ) -> Result<()> {
        let mut out = *scores;
        if self.su.use_library_score {
            let lib = MRMScoring::calc_library_score(feature, transitions)?;
            out.library_corr = lib.correlation;
            out.library_norm_manhattan = lib.norm_manhattan;
            out.library_manhattan = lib.manhattan;
            out.library_dotprod = lib.dotprod;
            out.library_sangle = lib.spectral_angle;
            out.library_rootmeansquare = lib.rmsd;
        }
        if self.su.use_rt_score {
            let rt_score = MRMScoring::calc_rt_score(compound, normalized_feature_rt);
            out.normalized_experimental_rt = normalized_feature_rt;
            out.raw_rt_score = rt_score;
            out.norm_rt_score = rt_score / self.rt_normalization_factor;
        }
        *scores = out;
        Ok(())
    }

    /// The number of spectra to merge: `add_up_spectra`, or in dynamic mode
    /// `max(1, ceil((rightWidth - leftWidth) * fraction))` from the feature's
    /// meta values.
    fn n_merge_spectra(&self, feature: &dyn MRMFeatureAccess) -> Result<i32> {
        match self.spectra_merge_method_type {
            SpectrumMergeMethodType::Fixed => Ok(self.add_up_spectra),
            SpectrumMergeMethodType::Dynamic => {
                let left = feature.meta_value("leftWidth")?;
                let right = feature.meta_value("rightWidth")?;
                let n = ((right - left) * self.merge_spectra_by_peak_width_fraction).ceil();
                if n.is_nan() || n > f64::from(i32::MAX) {
                    return Err(Error::InvalidValue(
                        "the number of spectra to merge does not fit an int".into(),
                    ));
                }
                // n is at most i32::MAX; a value below 1 becomes 1, as the source's max.
                Ok(if n >= 1.0 { n as i32 } else { 1 })
            }
        }
    }

    fn extraction(&self, diascoring: &DiaScoring) -> DriftExtraction {
        DriftExtraction {
            dia_extraction_window: diascoring.dia_extraction_window(),
            dia_extraction_ppm: diascoring.is_dia_extraction_ppm(),
            drift_extra: self.im_drift_extra_pcnt,
        }
    }

    /// Score a chromatographic feature with the DIA (SWATH) full-spectrum
    /// scores (source `calculateDIAScores`). Returns the per-transition signed
    /// mass errors in ppm (-1 where no signal was found), the source's
    /// `masserror_ppm` out-parameter.
    ///
    /// The spectra around the feature apex are fetched from `swath_maps`
    /// ([`OpenSwathScoring::fetch_spectrum_swath_maps_into`]), with the number
    /// of spectra given by the merge mode. In order: ion mobility scores
    /// (`use_im_scores`), the mass-difference score (always), the isotope and
    /// theoretical-spectrum scores (`use_ms2_isotope_scores`), the b/y ion
    /// scores for a peptide with a known sequence (`use_ionseries_scores`, at
    /// fragment charge 1), then with a non-empty `ms1_map` the precursor scores
    /// ([`OpenSwathScoring::calculate_precursor_dia_scores`]) and the MS1 ion
    /// mobility scores. MS1 extraction uses `im_range` only with
    /// `use_ms1_ion_mobility`.
    ///
    /// `sink` receives the extracted mobilograms, tagged with the feature RT
    /// and `feature_id` (source `mobilogram_consumer`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for no transitions or no SWATH maps (debug-only
    /// preconditions in the source), and the errors of the spectrum access,
    /// merging and the individual scores.
    #[allow(clippy::too_many_arguments)]
    pub fn calculate_dia_scores(
        &self,
        feature: &dyn MRMFeatureAccess,
        transitions: &[LightTransition],
        normalized_library_intensity: &[f64],
        swath_maps: &[SwathMap],
        ms1_map: Option<&SpectrumAccessPtr>,
        diascoring: &DiaScoring,
        compound: &LightCompound,
        scores: &mut OpenSwathScores,
        drift_target: f64,
        im_range: &RangeBase,
        mut sink: Option<&mut dyn MobilogramSink>,
        feature_id: i64,
    ) -> Result<Vec<f64>> {
        let first = transitions.first().ok_or_else(|| {
            Error::InvalidValue("There needs to be at least one transition.".into())
        })?;
        if swath_maps.is_empty() {
            return Err(Error::InvalidValue(
                "There needs to be at least one swath map.".into(),
            ));
        }
        let mut out = *scores;
        let n_merge = self.n_merge_spectra(feature)?;
        let rt = feature.rt();
        let spectra = self.fetch_spectrum_swath_maps_into(swath_maps, rt, n_merge, im_range)?;
        let extraction = self.extraction(diascoring);
        if self.su.use_im_scores {
            IonMobilityScoring::drift_scoring(
                &spectra,
                transitions,
                &mut out,
                drift_target,
                *im_range,
                &extraction,
                self.apply_im_peak_picking,
                export(&mut sink, rt, feature_id),
            )?;
        }
        let massdiff = diascoring.dia_massdiff_score(
            transitions,
            &spectra,
            normalized_library_intensity,
            im_range,
        )?;
        out.massdev_score = massdiff.ppm_score;
        out.weighted_massdev_score = massdiff.ppm_score_weighted;
        if self.su.use_ms2_isotope_scores {
            let pre = diascoring.score_with_isotopes(&spectra, transitions, im_range)?;
            out.dotprod_score_dia = pre.dotprod;
            out.manhatt_score_dia = pre.manhattan;
            let iso = diascoring.dia_isotope_scores(transitions, &spectra, feature, im_range)?;
            out.isotope_correlation = iso.isotope_corr;
            out.isotope_overlap = iso.isotope_overlap;
        }
        if compound.is_peptide() && !compound.sequence.is_empty() && self.su.use_ionseries_scores {
            let sequence = convert_peptide_to_aa_sequence(compound)?;
            let by = diascoring.dia_by_ion_score(&spectra, &sequence, 1, im_range)?;
            out.bseries_score = by.bseries_score;
            out.yseries_score = by.yseries_score;
        }
        let im_range_ms1 = if self.use_ms1_ion_mobility {
            *im_range
        } else {
            RangeBase::new()
        };
        let ms1 = ms1_map.filter(|m| m.nr_spectra() > 0);
        if let Some(ms1) = ms1 {
            let precursor_spectra = self.precursor_scores(
                ms1,
                diascoring,
                first.precursor_mz,
                rt,
                compound,
                im_range_ms1,
                &mut out,
            )?;
            if self.su.use_im_scores {
                let ms1_spectra = if n_merge == self.add_up_spectra {
                    precursor_spectra
                } else {
                    self.fetch_spectrum_swath_into(ms1.as_ref(), rt, n_merge, &im_range_ms1)?
                };
                IonMobilityScoring::drift_scoring_ms1(
                    &ms1_spectra,
                    transitions,
                    &mut out,
                    drift_target,
                    im_range_ms1,
                    &extraction,
                    export(&mut sink, rt, feature_id),
                )?;
                IonMobilityScoring::drift_scoring_ms1_contrast(
                    &spectra,
                    &ms1_spectra,
                    transitions,
                    &mut out,
                    im_range_ms1,
                    &extraction,
                    export(&mut sink, rt, feature_id),
                )?;
            }
        }
        *scores = out;
        Ok(massdiff.diff_ppm)
    }

    /// Score a feature with the precursor map (source
    /// `calculatePrecursorDIAScores`).
    ///
    /// Without `use_ms1_ion_mobility`, `im_range` is cleared first. With a
    /// non-empty `ms1_map`, `add_up_spectra` spectra around `rt` are fetched
    /// and the MS1 mass error (`ms1_ppm_score`) and isotope scores
    /// (`ms1_isotope_correlation`, `ms1_isotope_overlap`) computed. The
    /// isotope pattern comes from the peptide's formula at the precursor
    /// charge (1 when unset), from the metabolite's sum formula (its
    /// `sequence`) with that charge, or from an averagine model when no
    /// sequence is known.
    ///
    /// # Errors
    ///
    /// The errors of the spectrum access and the scores, and of parsing the
    /// peptide sequence or the sum formula.
    #[allow(clippy::too_many_arguments)]
    pub fn calculate_precursor_dia_scores(
        &self,
        ms1_map: Option<&SpectrumAccessPtr>,
        diascoring: &DiaScoring,
        precursor_mz: f64,
        rt: f64,
        compound: &LightCompound,
        im_range: RangeBase,
        scores: &mut OpenSwathScores,
    ) -> Result<()> {
        let im_range = if self.use_ms1_ion_mobility {
            im_range
        } else {
            RangeBase::new()
        };
        if let Some(ms1) = ms1_map.filter(|m| m.nr_spectra() > 0) {
            let mut out = *scores;
            self.precursor_scores(
                ms1,
                diascoring,
                precursor_mz,
                rt,
                compound,
                im_range,
                &mut out,
            )?;
            *scores = out;
        }
        Ok(())
    }

    /// The body of [`OpenSwathScoring::calculate_precursor_dia_scores`] for a
    /// non-empty map and an already adjusted range; returns the fetched MS1
    /// spectra for reuse by the MS1 ion mobility scores.
    #[allow(clippy::too_many_arguments)]
    fn precursor_scores(
        &self,
        ms1: &SpectrumAccessPtr,
        diascoring: &DiaScoring,
        precursor_mz: f64,
        rt: f64,
        compound: &LightCompound,
        im_range: RangeBase,
        out: &mut OpenSwathScores,
    ) -> Result<SpectrumSequence> {
        let spectra =
            self.fetch_spectrum_swath_into(ms1.as_ref(), rt, self.add_up_spectra, &im_range)?;
        out.ms1_ppm_score = diascoring
            .dia_ms1_massdiff_score(precursor_mz, &spectra, &im_range)?
            .ppm_score;
        let charge = if compound.charge != 0 {
            compound.charge
        } else {
            1
        };
        let iso = if compound.sequence.is_empty() {
            diascoring.dia_ms1_isotope_scores_averagine(
                precursor_mz,
                &spectra,
                charge,
                &im_range,
            )?
        } else if compound.is_peptide() {
            let formula = convert_peptide_to_aa_sequence(compound)?
                .formula_for(PeptideFragmentType::Full, charge)?;
            diascoring.dia_ms1_isotope_scores(precursor_mz, &spectra, &im_range, &formula)?
        } else {
            let formula = EmpiricalFormula::parse(&compound.sequence)?.with_charge(charge);
            diascoring.dia_ms1_isotope_scores(precursor_mz, &spectra, &im_range, &formula)?
        };
        out.ms1_isotope_correlation = iso.isotope_corr;
        out.ms1_isotope_overlap = iso.isotope_overlap;
        Ok(spectra)
    }

    /// Score an identification transition with the DIA scores (source
    /// `calculateDIAIdScores`).
    ///
    /// The spectra around the feature apex are fetched as in
    /// [`OpenSwathScoring::calculate_dia_scores`]. The transition's product
    /// m/z is scored like a precursor: averagine isotope scores at its product
    /// charge (1 when unset) into `isotope_correlation`/`isotope_overlap`, and
    /// the MS1-style mass error into `massdev_score`. With `use_im_scores`,
    /// [`IonMobilityScoring::drift_id_scoring`] scores it against the detection
    /// transitions.
    ///
    /// The source takes the detection transition group and only reads its
    /// transitions; this takes the transitions directly.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for no SWATH maps (a debug-only precondition in
    /// the source), and the errors of the spectrum access and the scores.
    #[allow(clippy::too_many_arguments)]
    pub fn calculate_dia_id_scores(
        &self,
        feature: &dyn MRMFeatureAccess,
        transition: &LightTransition,
        detection_transitions: &[LightTransition],
        swath_maps: &[SwathMap],
        im_range: &RangeBase,
        diascoring: &DiaScoring,
        scores: &mut OpenSwathScores,
        drift_target: f64,
        mut sink: Option<&mut dyn MobilogramSink>,
        feature_id: i64,
    ) -> Result<()> {
        if swath_maps.is_empty() {
            return Err(Error::InvalidValue(
                "There needs to be at least one swath map.".into(),
            ));
        }
        let mut out = *scores;
        let n_merge = self.n_merge_spectra(feature)?;
        let rt = feature.rt();
        let spectra = self.fetch_spectrum_swath_maps_into(swath_maps, rt, n_merge, im_range)?;
        let charge = match transition.product_charge_state() {
            0 => 1,
            c => c,
        };
        let iso = diascoring.dia_ms1_isotope_scores_averagine(
            transition.product_mz,
            &spectra,
            charge,
            im_range,
        )?;
        out.isotope_correlation = iso.isotope_corr;
        out.isotope_overlap = iso.isotope_overlap;
        out.massdev_score = diascoring
            .dia_ms1_massdiff_score(transition.product_mz, &spectra, im_range)?
            .ppm_score;
        if self.su.use_im_scores {
            IonMobilityScoring::drift_id_scoring(
                &spectra,
                std::slice::from_ref(transition),
                detection_transitions,
                &mut out,
                drift_target,
                *im_range,
                &self.extraction(diascoring),
                self.apply_im_peak_picking,
                export(&mut sink, rt, feature_id),
            )?;
        }
        *scores = out;
        Ok(())
    }

    /// Prepare the spectra for DIA analysis from one map (source
    /// `fetchSpectrumSwath(SpectrumAccessPtr, ...)` returning a sequence).
    ///
    /// Fetches up to `nr_spectra_to_add` spectra around `rt` (one more for an
    /// even count, see
    /// [`SpectrumAccess::multiple_spectra`]). With
    /// [`SpectrumAdditionMethod::Addition`] they are returned as they are and
    /// merged downstream by the scores; with
    /// [`SpectrumAdditionMethod::Resample`] one merged spectrum is returned,
    /// even when nothing was fetched: the merge of no spectra is an empty
    /// spectrum here, where the source's `addUpSpectra` indexes the empty
    /// vector (undefined behaviour). `im_range` is used only for resampling, to
    /// drift-filter.
    ///
    /// # Errors
    ///
    /// The errors of the spectrum access and of the merge.
    pub fn fetch_spectrum_swath(
        &self,
        swath_map: &dyn SpectrumAccess,
        rt: f64,
        nr_spectra_to_add: i32,
        im_range: &RangeBase,
    ) -> Result<SpectrumSequence> {
        let all = swath_map.multiple_spectra(rt, fetch_count(nr_spectra_to_add))?;
        match self.spectra_addition_method {
            SpectrumAdditionMethod::Addition => Ok(all),
            SpectrumAdditionMethod::Resample => Ok(vec![merge_spectra_for_scoring(
                &all,
                im_range,
                self.spacing_for_spectra_resampling,
            )?]),
        }
    }

    /// As [`OpenSwathScoring::fetch_spectrum_swath`], in the form the scoring
    /// calls use (source `fetchSpectrumSwath(SpectrumAccessPtr, ...,
    /// SpectrumSequence& out)`): when nothing is fetched the result is empty
    /// in both addition modes.
    ///
    /// # Errors
    ///
    /// As [`OpenSwathScoring::fetch_spectrum_swath`].
    pub fn fetch_spectrum_swath_into(
        &self,
        swath_map: &dyn SpectrumAccess,
        rt: f64,
        nr_spectra_to_add: i32,
        im_range: &RangeBase,
    ) -> Result<SpectrumSequence> {
        let all = swath_map.multiple_spectra(rt, fetch_count(nr_spectra_to_add))?;
        if all.is_empty() {
            return Ok(all);
        }
        self.fetch_spectrum_swath_merge(all, im_range)
    }

    fn fetch_spectrum_swath_merge(
        &self,
        all: SpectrumSequence,
        im_range: &RangeBase,
    ) -> Result<SpectrumSequence> {
        match self.spectra_addition_method {
            SpectrumAdditionMethod::Addition => Ok(all),
            SpectrumAdditionMethod::Resample => Ok(vec![merge_spectra_for_scoring(
                &all,
                im_range,
                self.spacing_for_spectra_resampling,
            )?]),
        }
    }

    /// Per SWATH map, the fetched spectra concatenated (skipping maps with
    /// none), drift-filtered to `im_range` when it is not empty.
    fn per_map(
        &self,
        swath_maps: &[SwathMap],
        rt: f64,
        n: usize,
        im_range: &RangeBase,
    ) -> Result<Vec<SpectrumSequence>> {
        let mut out = Vec::with_capacity(swath_maps.len());
        for map in swath_maps {
            let access = map_access(map)?;
            let sequence = if im_range.is_empty() {
                access.multiple_spectra(rt, n)?
            } else {
                let (lo, hi) = range_bounds(im_range);
                access.multiple_spectra_in_drift_range(rt, n, lo, hi)?
            };
            if !sequence.is_empty() {
                out.push(sequence);
            }
        }
        Ok(out)
    }

    /// Prepare the spectra for DIA analysis from one or more maps (source
    /// `fetchSpectrumSwath(const std::vector<SwathMap>&, ...)` returning a
    /// sequence).
    ///
    /// One map behaves as [`OpenSwathScoring::fetch_spectrum_swath`]. Several
    /// maps (SONAR data) always give at most one spectrum: with an IM range,
    /// each map's drift-filtered spectra are concatenated, and the results
    /// concatenated; without, each map's spectra are concatenated (addition)
    /// or merged (resampling), and the results merged by resampling.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for no maps (a debug-only precondition in the
    /// source), [`Error::MissingInformation`] for a map without spectrum
    /// access, and the errors of the access and the merge.
    pub fn fetch_spectrum_swath_maps(
        &self,
        swath_maps: &[SwathMap],
        rt: f64,
        nr_spectra_to_add: i32,
        im_range: &RangeBase,
    ) -> Result<SpectrumSequence> {
        match swath_maps {
            [] => Err(Error::InvalidValue(
                "swath_maps vector cannot be empty".into(),
            )),
            [single] => self.fetch_spectrum_swath(
                map_access(single)?.as_ref(),
                rt,
                nr_spectra_to_add,
                im_range,
            ),
            _ => self.fetch_spectrum_swath_maps_multi(swath_maps, rt, nr_spectra_to_add, im_range),
        }
    }

    /// As [`OpenSwathScoring::fetch_spectrum_swath_maps`], in the form the
    /// scoring calls use (source `fetchSpectrumSwath(const
    /// std::vector<SwathMap>&, ..., SpectrumSequence& out)`): a single map
    /// behaves as [`OpenSwathScoring::fetch_spectrum_swath_into`]; several
    /// maps as in the returning form.
    ///
    /// # Errors
    ///
    /// As [`OpenSwathScoring::fetch_spectrum_swath_maps`].
    pub fn fetch_spectrum_swath_maps_into(
        &self,
        swath_maps: &[SwathMap],
        rt: f64,
        nr_spectra_to_add: i32,
        im_range: &RangeBase,
    ) -> Result<SpectrumSequence> {
        match swath_maps {
            [] => Err(Error::InvalidValue(
                "swath_maps vector cannot be empty".into(),
            )),
            [single] => self.fetch_spectrum_swath_into(
                map_access(single)?.as_ref(),
                rt,
                nr_spectra_to_add,
                im_range,
            ),
            _ => self.fetch_spectrum_swath_maps_multi(swath_maps, rt, nr_spectra_to_add, im_range),
        }
    }

    fn fetch_spectrum_swath_maps_multi(
        &self,
        swath_maps: &[SwathMap],
        rt: f64,
        nr_spectra_to_add: i32,
        im_range: &RangeBase,
    ) -> Result<SpectrumSequence> {
        let n = fetch_count(nr_spectra_to_add);
        let empty = RangeBase::new();
        let per_map = self.per_map(swath_maps, rt, n, im_range)?;
        let mut merged = Vec::with_capacity(per_map.len());
        for sequence in &per_map {
            let spectrum = if !im_range.is_empty()
                || self.spectra_addition_method == SpectrumAdditionMethod::Addition
            {
                concatenate_spectra(sequence)?
            } else {
                merge_spectra_for_scoring(sequence, &empty, self.spacing_for_spectra_resampling)?
            };
            merged.push(spectrum);
        }
        if merged.is_empty() {
            return Ok(Vec::new());
        }
        let result = if im_range.is_empty() {
            merge_spectra_for_scoring(&merged, &empty, self.spacing_for_spectra_resampling)?
        } else {
            concatenate_spectra(&merged)?
        };
        Ok(vec![result])
    }
}

fn export<'a>(
    sink: &'a mut Option<&mut dyn MobilogramSink>,
    feature_rt: f64,
    feature_id: i64,
) -> Option<MobilogramExport<'a>> {
    sink.as_deref_mut().map(|sink| MobilogramExport {
        sink,
        feature_rt,
        feature_id,
    })
}
