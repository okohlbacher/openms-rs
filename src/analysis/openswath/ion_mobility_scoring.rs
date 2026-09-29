// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The OpenSWATH ion mobility scoring routines, from
//! `ANALYSIS/OPENSWATH/IonMobilityScoring.h` (class `IonMobilityScoring` and
//! the free function `sumAlignedMobilograms`).
//!
//! These scores use the ion mobilograms of individual peptides in one (or
//! more) frames:
//!
//! - [`IonMobilityScoring::drift_scoring`](crate::analysis::openswath::ion_mobility_scoring::IonMobilityScoring::drift_scoring) scores fragment ion mobilograms
//!   extracted from a DIA frame;
//! - [`IonMobilityScoring::drift_scoring_ms1`](crate::analysis::openswath::ion_mobility_scoring::IonMobilityScoring::drift_scoring_ms1) scores the precursor ion
//!   mobility in an MS1 frame;
//! - [`IonMobilityScoring::drift_scoring_ms1_contrast`](crate::analysis::openswath::ion_mobility_scoring::IonMobilityScoring::drift_scoring_ms1_contrast) cross-correlates
//!   precursor against fragment ion mobilograms;
//! - [`IonMobilityScoring::drift_id_scoring`](crate::analysis::openswath::ion_mobility_scoring::IonMobilityScoring::drift_id_scoring) scores an identification
//!   transition against the detection transitions.
//!
//! Scores are written into an [`OpenSwathScores`](crate::kernel::mrm::OpenSwathScores). On an error the scores are
//! left unchanged.
//!
//! # Gaps
//!
//! - Ion-mobilogram peak picking (`apply_im_peak_picking = true`) needs
//!   `ANALYSIS/OPENSWATH/PeakPickerMobilogram.h`, which is not ported yet; the
//!   request fails with [`Error::Unsupported`](crate::Error::Unsupported) naming that dependency.
//! - The optional `MobilogramParquetConsumer` export (from
//!   `FORMAT/DATAACCESS/MobilogramParquetConsumer.h`) is represented by the
//!   [`MobilogramSink`](crate::analysis::openswath::ion_mobility_scoring::MobilogramSink) trait; the Parquet implementation of it is a backlog
//!   item of the format port.
//!
//! See `docs/OPENSWATH_SCORING_SUPPORT.md`.

use std::collections::BTreeMap;

use super::data_structures::SpectrumPtr;
use super::dia_helper::{create_mz_range_ppm, integrate_window_sequence};
use super::mrm_scoring::MRMScoring;
use super::transition_experiment::LightTransition;
use crate::kernel::mrm::OpenSwathScores;
use crate::kernel::ranges::RangeBase;
use crate::kernel::{MobilityPeak1D, Mobilogram};
use crate::{Error, Result};

/// A consumer of the aligned ion mobilograms the scoring routines extract
/// (source `MobilogramParquetConsumer::consumeMobilogram`, whose argument list
/// this mirrors).
pub trait MobilogramSink {
    /// Consume one mobilogram.
    ///
    /// `mobilogram_type` is for example `"ms1"`, `"ms2"`, `"ms1_contrast"` or
    /// `"identification"`; `ms_level` is 1 or 2; `transition_id` is -1 when
    /// unset (the scoring routines always pass -1); `transition_native_id` is
    /// empty when unset; `feature_rt` is NaN and `feature_id` -1 when unset.
    ///
    /// # Errors
    ///
    /// Implementation-specific; an error aborts the scoring call.
    #[allow(clippy::too_many_arguments)]
    fn consume_mobilogram(
        &mut self,
        mobilogram: &Mobilogram,
        mobilogram_type: &str,
        ms_level: i64,
        transition_id: i64,
        transition_native_id: &str,
        feature_rt: f64,
        feature_id: i64,
    ) -> Result<()>;
}

/// Where and how to export mobilograms: the source's three trailing optional
/// arguments `mobilogram_consumer`, `feature_rt` and `feature_id`.
pub struct MobilogramExport<'a> {
    /// The consumer.
    pub sink: &'a mut dyn MobilogramSink,
    /// Retention-time apex the mobilograms belong to; NaN when unset.
    pub feature_rt: f64,
    /// Feature id of that apex; -1 when unset.
    pub feature_id: i64,
}

/// The shared extraction settings of the scoring routines.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DriftExtraction {
    /// m/z extraction width, in Th or ppm.
    pub dia_extraction_window: f64,
    /// Whether [`Self::dia_extraction_window`] is in ppm.
    pub dia_extraction_ppm: bool,
    /// Fraction by which the ion-mobility range is widened on **each** side
    /// (the range is scaled by `2 * drift_extra + 1`).
    pub drift_extra: f64,
}

/// An extracted ion mobilogram with its weighted mean ion mobility and total
/// intensity (the outputs of source `computeIonMobilogram`).
#[derive(Clone, Debug, PartialEq)]
pub struct IonMobilogram {
    /// Intensity-weighted mean ion mobility; -1 without signal.
    pub im: f64,
    /// Total intensity; 0 without signal.
    pub intensity: f64,
    /// The mobilogram, one peak per ion-mobility bin in ascending order.
    pub mobilogram: Mobilogram,
}

/// The OpenSWATH ion mobility scoring routines (source `IonMobilityScoring`,
/// whose members are all static).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct IonMobilityScoring;

fn has_drift(spectra: &[SpectrumPtr]) -> bool {
    spectra.iter().all(|s| s.drift_time_array().is_some())
}

/// Widen the range by `drift_extra` on each side and derive the grid epsilon,
/// `max(1e-8, span * 1e-5)`, which is scale-invariant across IM units (about
/// 1e-5 for VSSC ranges of 0.8-1.5, about 0.002 for CCS ranges of 300-500).
fn widen(mut im_range: RangeBase, drift_extra: f64) -> Result<(RangeBase, f64)> {
    im_range.scale_by(drift_extra * 2.0 + 1.0)?;
    let span = im_range.span().unwrap_or(f64::NEG_INFINITY);
    Ok((im_range, f64::max(1e-8, span * 1e-5)))
}

fn peak_picking_unsupported() -> Error {
    Error::Unsupported(
        "ion mobilogram peak picking needs ANALYSIS/OPENSWATH/PeakPickerMobilogram.h, which is not ported yet"
            .into(),
    )
}

fn intensities(mobilogram: &Mobilogram) -> Vec<f64> {
    mobilogram
        .peaks
        .iter()
        .map(|p| f64::from(p.intensity))
        .collect()
}

/// Sum the aligned traces point by point.
fn sum_traces(reference_len: usize, traces: &[Vec<f64>]) -> Vec<f64> {
    let mut out = vec![0.0; reference_len];
    for (k, value) in out.iter_mut().enumerate() {
        for trace in traces {
            *value += trace.get(k).copied().unwrap_or(0.0);
        }
    }
    out
}

impl IonMobilityScoring {
    /// A scorer; the source class has only static members.
    pub fn new() -> Self {
        Self
    }

    /// Compute the ion mobilogram of `spectra` in `mz_range` and `im_range`,
    /// and integrate its intensity (source `computeIonMobilogram`).
    ///
    /// Peaks with `mz_range.min <= mz < mz_range.max` whose ion mobility lies in
    /// `im_range` contribute; they are binned by `trunc(im / eps)` across the
    /// whole sequence, and the bins, in ascending order, become the mobilogram
    /// peaks at mobility `bin * eps`. Intensities are summed in `f64` per bin
    /// and stored as the mobilogram's `f32`, as the source does. Without signal
    /// `im` is -1 and `intensity` 0.
    ///
    /// # Errors
    ///
    /// - [`Error::MissingInformation`] for a spectrum without ion-mobility
    ///   array (a debug-only precondition in the source, which otherwise
    ///   dereferences a null pointer).
    /// - [`Error::InvalidValue`] when the intensity or ion-mobility array is
    ///   shorter than the m/z array (debug-only in the source), for a
    ///   non-positive `eps`, or when a bin index does not fit the source's
    ///   `int` (undefined behaviour in the source).
    pub fn compute_ion_mobilogram(
        spectra: &[SpectrumPtr],
        mz_range: &RangeBase,
        im_range: &RangeBase,
        eps: f64,
    ) -> Result<IonMobilogram> {
        if !(eps.is_finite() && eps > 0.0) {
            return Err(Error::InvalidValue(
                "eps must be positive and finite".into(),
            ));
        }
        let multiplier = 1.0 / eps;
        let mut bins: BTreeMap<i32, f64> = BTreeMap::new();
        let (mut im, mut intensity) = (0.0, 0.0);
        let (lo, hi) = match (mz_range.min(), mz_range.max()) {
            (Ok(lo), Ok(hi)) => (lo, hi),
            _ => (f64::MAX, -f64::MAX),
        };
        for spectrum in spectra {
            let drift = spectrum.drift_time_array().ok_or_else(|| {
                Error::MissingInformation(
                    "Cannot filter by drift time if no drift time is available.".into(),
                )
            })?;
            let mz = &spectrum.mz_array().data;
            let ints = &spectrum.intensity_array().data;
            if ints.len() < mz.len() || drift.data.len() < mz.len() {
                return Err(Error::InvalidValue(
                    "MZ, Intensity and Drift Time arrays need to have the same length.".into(),
                ));
            }
            let start = mz.partition_point(|&v| v < lo);
            for i in start..mz.len() {
                if mz[i] >= hi {
                    break;
                }
                let drift_i = drift.data[i];
                if im_range.contains(drift_i) {
                    intensity += ints[i];
                    im += ints[i] * drift_i;
                    let scaled = (drift_i * multiplier).trunc();
                    if !(scaled >= f64::from(i32::MIN) && scaled <= f64::from(i32::MAX)) {
                        return Err(Error::InvalidValue(
                            "ion mobility bin does not fit a 32-bit integer".into(),
                        ));
                    }
                    *bins.entry(scaled as i32).or_insert(0.0) += ints[i];
                }
            }
        }
        if intensity > 0.0 {
            im /= intensity;
        } else {
            im = -1.0;
            intensity = 0.0;
        }
        let peaks = bins
            .into_iter()
            .map(|(bin, value)| MobilityPeak1D {
                mobility: f64::from(bin) / multiplier,
                intensity: value as f32,
            })
            .collect();
        Ok(IonMobilogram {
            im,
            intensity,
            mobilogram: Mobilogram::from_peaks(peaks),
        })
    }

    /// Discretise the ion mobilities of all `mobilograms` into a grid: sorted,
    /// keeping a value only when it exceeds its predecessor by more than `eps`
    /// (source private `computeGrid_`).
    pub fn compute_grid(mobilograms: &[Mobilogram], eps: f64) -> Vec<f64> {
        let mut values: Vec<f64> = mobilograms
            .iter()
            .flat_map(|m| m.peaks.iter().map(|p| p.mobility))
            .collect();
        values.sort_by(f64::total_cmp);
        let mut grid = Vec::new();
        if let Some(&first) = values.first() {
            grid.push(first);
            for pair in values.windows(2) {
                if (pair[1] - pair[0]).abs() > eps {
                    grid.push(pair[1]);
                }
            }
        }
        grid
    }

    /// Project `profile` onto `im_grid` (source private `alignToGrid_`).
    ///
    /// Walking the grid, a profile point within `10 * eps` of the grid value is
    /// taken over and the profile advances; otherwise a zero-intensity point at
    /// the grid value is inserted. Returns the aligned mobilogram and the grid
    /// index of the maximum, which the source computes from the *next* profile
    /// point after each step (kept as is).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when the profile falls behind the grid by more
    /// than `10 * eps`, where the source logs "This should never happen" and
    /// throws `Exception::OutOfRange`.
    pub fn align_to_grid(
        profile: &Mobilogram,
        im_grid: &[f64],
        eps: f64,
    ) -> Result<(Mobilogram, usize)> {
        let mut it = profile.peaks.iter().peekable();
        let mut max_peak_idx = 0;
        let mut max_int = 0.0_f32;
        let mut aligned = Vec::with_capacity(im_grid.len());
        for (k, &grid_value) in im_grid.iter().enumerate() {
            let peak = match it.peek() {
                Some(p) if (p.mobility - grid_value).abs() < eps * 10.0 => {
                    let p = **p;
                    it.next();
                    p
                }
                _ => MobilityPeak1D {
                    mobility: grid_value,
                    intensity: 0.0,
                },
            };
            if let Some(p) = it.peek() {
                if grid_value - p.mobility > eps * 10.0 {
                    return Err(Error::InvalidRange(format!(
                        "the profile has advanced past the grid: {grid_value} / {}",
                        p.mobility
                    )));
                }
                if p.intensity > max_int {
                    max_int = p.intensity;
                    max_peak_idx = k;
                }
            }
            aligned.push(peak);
        }
        Ok((Mobilogram::from_peaks(aligned), max_peak_idx))
    }

    /// The intensities of each mobilogram (source private
    /// `extractIntensities(const std::vector<Mobilogram>&, ...)`).
    pub fn extract_intensities(mobilograms: &[Mobilogram]) -> Vec<Vec<f64>> {
        mobilograms.iter().map(intensities).collect()
    }

    /// The intensities of one mobilogram (source private
    /// `extractIntensities(const Mobilogram&)`).
    pub fn extract_intensities_single(mobilogram: &Mobilogram) -> Vec<f64> {
        intensities(mobilogram)
    }

    fn aligned(mobilograms: &[Mobilogram], grid: &[f64], eps: f64) -> Result<Vec<Mobilogram>> {
        mobilograms
            .iter()
            .map(|m| Ok(Self::align_to_grid(m, grid, eps)?.0))
            .collect()
    }

    /// Score the ion-mobility dimension in MS2 (source `driftScoring`).
    ///
    /// The range `im_range` is widened by `drift_extra` on each side. For every
    /// transition the fragment mobilogram is extracted; transitions with signal
    /// contribute to the mean measured drift `im_drift`, the
    /// intensity-weighted drift `im_drift_weighted`, the mean signed and
    /// absolute deltas to `drift_target` (`im_delta`, `im_delta_score`) and
    /// `im_log_intensity = ln(total + 1)`; without any, those are -1 (and the
    /// log intensity 0). The mobilograms are aligned to a common grid and, with
    /// at least two, their cross-correlation coelution and shape scores
    /// (`im_xcorr_coelution_score`, `im_xcorr_shape_score`, which can be NaN)
    /// are set; with fewer, 0 and NaN.
    ///
    /// If any spectrum lacks an ion-mobility array, nothing is scored and
    /// `Ok(())` is returned, as the source returns early (logging at debug
    /// level).
    ///
    /// # Errors
    ///
    /// [`Error::Unsupported`] when `apply_im_peak_picking` is set (see the
    /// module gaps), the errors of [`Self::compute_ion_mobilogram`],
    /// [`Self::align_to_grid`], the cross-correlation and the sink.
    #[allow(clippy::too_many_arguments)]
    pub fn drift_scoring(
        spectra: &[SpectrumPtr],
        transitions: &[LightTransition],
        scores: &mut OpenSwathScores,
        drift_target: f64,
        im_range: RangeBase,
        extraction: &DriftExtraction,
        apply_im_peak_picking: bool,
        export: Option<MobilogramExport<'_>>,
    ) -> Result<()> {
        if !has_drift(spectra) {
            return Ok(());
        }
        let (im_range, eps) = widen(im_range, extraction.drift_extra)?;
        let mut out = *scores;
        let (mut delta, mut delta_abs, mut computed, mut weighted, mut sum_int) =
            (0.0, 0.0, 0.0, 0.0, 0.0);
        let mut used = 0u32;
        let mut mobilograms = Vec::with_capacity(transitions.len());
        for transition in transitions {
            let mz_range = create_mz_range_ppm(
                transition.product_mz,
                extraction.dia_extraction_window,
                extraction.dia_extraction_ppm,
            )?;
            let r = Self::compute_ion_mobilogram(spectra, &mz_range, &im_range, eps)?;
            mobilograms.push(r.mobilogram);
            if r.intensity <= 0.0 {
                continue;
            }
            used += 1;
            delta_abs += (drift_target - r.im).abs();
            delta += drift_target - r.im;
            computed += r.im;
            weighted += r.im * r.intensity;
            sum_int += r.intensity;
        }
        if used != 0 {
            let n = f64::from(used);
            delta /= n;
            delta_abs /= n;
            computed /= n;
            weighted /= sum_int;
        } else {
            (delta, delta_abs, computed, weighted) = (-1.0, -1.0, -1.0, -1.0);
        }
        out.im_delta_score = delta_abs;
        out.im_delta = delta;
        out.im_drift = computed;
        out.im_drift_weighted = weighted;
        out.im_log_intensity = (sum_int + 1.0).ln();

        let grid = Self::compute_grid(&mobilograms, eps);
        let aligned = Self::aligned(&mobilograms, &grid, eps)?;
        if let Some(export) = export {
            for (m, t) in aligned.iter().zip(transitions) {
                export.sink.consume_mobilogram(
                    m,
                    "ms2",
                    2,
                    -1,
                    &t.transition_name,
                    export.feature_rt,
                    export.feature_id,
                )?;
            }
        }
        if apply_im_peak_picking {
            return Err(peak_picking_unsupported());
        }
        if aligned.len() < 2 {
            out.im_xcorr_coelution_score = 0.0;
            out.im_xcorr_shape_score = f64::NAN;
        } else {
            let mut mrm = MRMScoring::new();
            mrm.initialize_xcorr_matrix_from_data(&Self::extract_intensities(&aligned))?;
            out.im_xcorr_coelution_score = mrm.calc_xcorr_coelution_score();
            out.im_xcorr_shape_score = mrm.calc_xcorr_shape_score();
        }
        *scores = out;
        Ok(())
    }

    /// Score the precursor ion mobility in MS1 (source `driftScoringMS1`).
    ///
    /// Integrates the widened `im_range` window at the first transition's
    /// precursor m/z and sets `im_ms1_drift` to the measured ion mobility (-1
    /// without signal), `im_ms1_delta = drift_target - im` and
    /// `im_ms1_delta_score = |im_ms1_delta|`. With `export`, the MS1
    /// mobilogram is also extracted and handed to the sink, as the source does
    /// only when a consumer is given.
    ///
    /// Returns `Ok(())` without scoring when a spectrum lacks an ion-mobility
    /// array, as the source does.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for no transitions (a debug-only precondition in
    /// the source, which otherwise reads out of bounds), and the integration,
    /// extraction and sink errors.
    pub fn drift_scoring_ms1(
        spectra: &[SpectrumPtr],
        transitions: &[LightTransition],
        scores: &mut OpenSwathScores,
        drift_target: f64,
        im_range: RangeBase,
        extraction: &DriftExtraction,
        export: Option<MobilogramExport<'_>>,
    ) -> Result<()> {
        if !has_drift(spectra) {
            return Ok(());
        }
        let first = transitions
            .first()
            .ok_or_else(|| Error::InvalidValue("Need at least one transition".into()))?;
        let (im_range, eps) = widen(im_range, extraction.drift_extra)?;
        let mz_range = create_mz_range_ppm(
            first.precursor_mz,
            extraction.dia_extraction_window,
            extraction.dia_extraction_ppm,
        )?;
        let im =
            integrate_window_sequence(spectra, &mz_range, &im_range, false)?.map_or(-1.0, |w| w.im);
        if let Some(export) = export {
            let r = Self::compute_ion_mobilogram(spectra, &mz_range, &im_range, eps)?;
            export.sink.consume_mobilogram(
                &r.mobilogram,
                "ms1",
                1,
                -1,
                "",
                export.feature_rt,
                export.feature_id,
            )?;
        }
        scores.im_ms1_drift = im;
        scores.im_ms1_delta_score = (drift_target - im).abs();
        scores.im_ms1_delta = drift_target - im;
        Ok(())
    }

    /// Cross-correlate the precursor ion mobilogram with the fragment ion
    /// mobilograms (source `driftScoringMS1Contrast`).
    ///
    /// The fragment mobilograms of every transition and the precursor
    /// mobilogram at the first transition's precursor m/z are aligned to one
    /// grid built from all of them. `im_ms1_contrast_coelution` and
    /// `im_ms1_contrast_shape` score the precursor against each fragment;
    /// `im_ms1_sum_contrast_coelution` and `im_ms1_sum_contrast_shape` score it
    /// against the summed fragments.
    ///
    /// Returns `Ok(())` without scoring when a spectrum of either sequence
    /// lacks an ion-mobility array, as the source does.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for no transitions (debug-only in the source),
    /// and the extraction, alignment, cross-correlation and sink errors.
    pub fn drift_scoring_ms1_contrast(
        spectra: &[SpectrumPtr],
        ms1_spectra: &[SpectrumPtr],
        transitions: &[LightTransition],
        scores: &mut OpenSwathScores,
        im_range: RangeBase,
        extraction: &DriftExtraction,
        export: Option<MobilogramExport<'_>>,
    ) -> Result<()> {
        if !has_drift(spectra) || !has_drift(ms1_spectra) {
            return Ok(());
        }
        let first = transitions
            .first()
            .ok_or_else(|| Error::InvalidValue("Need at least one transition".into()))?;
        let (im_range, eps) = widen(im_range, extraction.drift_extra)?;
        let mut mobilograms = Vec::with_capacity(transitions.len() + 1);
        for transition in transitions {
            let mz_range = create_mz_range_ppm(
                transition.product_mz,
                extraction.dia_extraction_window,
                extraction.dia_extraction_ppm,
            )?;
            mobilograms
                .push(Self::compute_ion_mobilogram(spectra, &mz_range, &im_range, eps)?.mobilogram);
        }
        let mz_range = create_mz_range_ppm(
            first.precursor_mz,
            extraction.dia_extraction_window,
            extraction.dia_extraction_ppm,
        )?;
        let ms1_profile =
            Self::compute_ion_mobilogram(ms1_spectra, &mz_range, &im_range, eps)?.mobilogram;
        mobilograms.push(ms1_profile);
        let grid = Self::compute_grid(&mobilograms, eps);
        let ms1_profile = mobilograms.pop().unwrap_or_default();
        let aligned_ms2 = Self::aligned(&mobilograms, &grid, eps)?;
        let (aligned_ms1, _) = Self::align_to_grid(&ms1_profile, &grid, eps)?;
        if let Some(export) = export {
            for (m, t) in aligned_ms2.iter().zip(transitions) {
                export.sink.consume_mobilogram(
                    m,
                    "ms1_contrast",
                    2,
                    -1,
                    &t.transition_name,
                    export.feature_rt,
                    export.feature_id,
                )?;
            }
            export.sink.consume_mobilogram(
                &aligned_ms1,
                "ms1",
                1,
                -1,
                "",
                export.feature_rt,
                export.feature_id,
            )?;
        }
        let ms1_int = intensities(&aligned_ms1);
        let ms2_int = Self::extract_intensities(&aligned_ms2);
        let (coelution, shape, sum_coelution, sum_shape) = contrast_scores(&ms1_int, &ms2_int)?;
        scores.im_ms1_contrast_coelution = coelution;
        scores.im_ms1_contrast_shape = shape;
        scores.im_ms1_sum_contrast_coelution = sum_coelution;
        scores.im_ms1_sum_contrast_shape = sum_shape;
        Ok(())
    }

    /// Score an identification transition against the detection transitions
    /// in the ion-mobility dimension (source `driftIdScoring`).
    ///
    /// The first of `transitions` is the identification transition: its
    /// measured drift, deltas to `drift_target` and `ln(1 + intensity)` set
    /// `im_drift`, `im_delta`, `im_delta_score` and `im_log_intensity`. When it
    /// has signal, its mobilogram and those of `detection_transitions` are
    /// aligned to one grid and the contrast scores `im_ind_contrast_*` (against
    /// each detection transition) and `im_ind_sum_contrast_*` (against their
    /// sum) are set; without signal those four are 0.
    ///
    /// The source takes the detection transition group
    /// (`MRMTransitionGroup<MSChromatogram, LightTransition>`) and only reads
    /// its transitions, so this takes the transitions directly.
    ///
    /// Returns `Ok(())` without scoring when a spectrum lacks an ion-mobility
    /// array, as the source does.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for no transitions (the source reads out of
    /// bounds), [`Error::Unsupported`] when `apply_im_peak_picking` is set and
    /// the identification transition has signal, and the extraction,
    /// alignment, cross-correlation and sink errors.
    #[allow(clippy::too_many_arguments)]
    pub fn drift_id_scoring(
        spectra: &[SpectrumPtr],
        transitions: &[LightTransition],
        detection_transitions: &[LightTransition],
        scores: &mut OpenSwathScores,
        drift_target: f64,
        im_range: RangeBase,
        extraction: &DriftExtraction,
        apply_im_peak_picking: bool,
        export: Option<MobilogramExport<'_>>,
    ) -> Result<()> {
        if !has_drift(spectra) {
            return Ok(());
        }
        let first = transitions
            .first()
            .ok_or_else(|| Error::InvalidValue("Need at least one transition".into()))?;
        let (im_range, eps) = widen(im_range, extraction.drift_extra)?;
        let mut out = *scores;
        let id_range = create_mz_range_ppm(
            first.product_mz,
            extraction.dia_extraction_window,
            extraction.dia_extraction_ppm,
        )?;
        let identification = Self::compute_ion_mobilogram(spectra, &id_range, &im_range, eps)?;
        out.im_drift = identification.im;
        out.im_delta_score = (drift_target - identification.im).abs();
        out.im_delta = drift_target - identification.im;
        out.im_log_intensity = identification.intensity.ln_1p();

        let mut mobilograms = Vec::with_capacity(detection_transitions.len() + 1);
        for detection in detection_transitions {
            let mz_range = create_mz_range_ppm(
                detection.product_mz,
                extraction.dia_extraction_window,
                extraction.dia_extraction_ppm,
            )?;
            mobilograms
                .push(Self::compute_ion_mobilogram(spectra, &mz_range, &im_range, eps)?.mobilogram);
        }
        if identification.im != -1.0 {
            mobilograms.push(identification.mobilogram.clone());
            let grid = Self::compute_grid(&mobilograms, eps);
            mobilograms.pop();
            let aligned = Self::aligned(&mobilograms, &grid, eps)?;
            let (aligned_id, _) = Self::align_to_grid(&identification.mobilogram, &grid, eps)?;
            if let Some(export) = export {
                for (m, t) in aligned.iter().zip(detection_transitions) {
                    export.sink.consume_mobilogram(
                        m,
                        "ms2",
                        2,
                        -1,
                        &t.transition_name,
                        export.feature_rt,
                        export.feature_id,
                    )?;
                }
                export.sink.consume_mobilogram(
                    &aligned_id,
                    "identification",
                    2,
                    -1,
                    &first.transition_name,
                    export.feature_rt,
                    export.feature_id,
                )?;
            }
            if apply_im_peak_picking {
                return Err(peak_picking_unsupported());
            }
            let (coelution, shape, sum_coelution, sum_shape) = contrast_scores(
                &intensities(&aligned_id),
                &Self::extract_intensities(&aligned),
            )?;
            out.im_ind_contrast_coelution = coelution;
            out.im_ind_contrast_shape = shape;
            out.im_ind_sum_contrast_coelution = sum_coelution;
            out.im_ind_sum_contrast_shape = sum_shape;
        } else {
            out.im_ind_contrast_coelution = 0.0;
            out.im_ind_contrast_shape = 0.0;
            out.im_ind_sum_contrast_coelution = 0.0;
            out.im_ind_sum_contrast_shape = 0.0;
        }
        *scores = out;
        Ok(())
    }
}

/// Contrast coelution and shape of `reference` against each trace and against
/// their point-wise sum: `(coelution, shape, sum_coelution, sum_shape)`.
///
/// The source comments the summed form as "horribly broken" because it passes
/// a 1x1 matrix to a score whose debug assertion wants 1x2, and calls the
/// `SumFrag` variants, which compute the same reductions, to avoid the
/// assertion; the port does the same.
fn contrast_scores(reference: &[f64], traces: &[Vec<f64>]) -> Result<(f64, f64, f64, f64)> {
    let mut mrm = MRMScoring::new();
    mrm.initialize_xcorr_precursor_contrast_matrix_from_data(&[reference.to_vec()], traces)?;
    let coelution = mrm.calc_xcorr_precursor_contrast_coelution_score();
    let shape = mrm.calc_xcorr_precursor_contrast_shape_score();
    let summed = sum_traces(reference.len(), traces);
    let mut sum = MRMScoring::new();
    sum.initialize_xcorr_precursor_contrast_matrix_from_data(&[reference.to_vec()], &[summed])?;
    Ok((
        coelution,
        shape,
        sum.calc_xcorr_precursor_contrast_sum_frag_coelution_score(),
        sum.calc_xcorr_precursor_contrast_sum_frag_shape_score(),
    ))
}

/// Sum aligned mobilograms point by point into one mobilogram with the first
/// mobilogram's mobilities (source free function `sumAlignedMobilograms`).
///
/// Intensities are added in `f32`, as the source adds `float` peak
/// intensities. An empty input gives an empty mobilogram.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the mobilograms differ in length; the source
/// checks this only in debug builds and otherwise skips missing points.
pub fn sum_aligned_mobilograms(aligned: &[Mobilogram]) -> Result<Mobilogram> {
    let Some(first) = aligned.first() else {
        return Ok(Mobilogram::default());
    };
    if aligned.iter().any(|m| m.peaks.len() != first.peaks.len()) {
        return Err(Error::InvalidValue(
            "All Mobilograms in aligned_mobilograms must have the same size.".into(),
        ));
    }
    let peaks = (0..first.peaks.len())
        .map(|j| MobilityPeak1D {
            mobility: first.peaks[j].mobility,
            intensity: aligned
                .iter()
                .fold(0.0_f32, |acc, m| acc + m.peaks[j].intensity),
        })
        .collect();
    Ok(Mobilogram::from_peaks(peaks))
}
