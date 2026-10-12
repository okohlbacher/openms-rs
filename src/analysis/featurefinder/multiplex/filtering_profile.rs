// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Filtering centroided and spline-interpolated profile data for multiplex
//! peak patterns.
//!
//! Native equivalent of `FEATUREFINDER/MultiplexFilteringProfile.h`. See
//! `docs/FEATUREFINDER_MULTIPLEX_SUPPORT.md`.

use std::ops::{Deref, DerefMut};

use super::filtered_ms_experiment::MultiplexFilteredMSExperiment;
use super::filtered_peak::{MultiplexFilteredPeak, SatelliteMultimap};
use super::filtering::{MultiplexFiltering, rt_begin, rt_end};
use super::isotopic_peak_pattern::MultiplexIsotopicPeakPattern;
use super::satellite::MultiplexSatelliteProfile;
use crate::concept::progress_logger::ProgressReporter;
use crate::kernel::MSExperiment;
use crate::processing::peak_picking::PeakBoundary;
use crate::processing::spline::spline_interpolated_peaks::{
    DEFAULT_NAVIGATOR_SCALING, SplineInterpolatedPeaks,
};
use crate::{Error, Result};

/// Filters centroided data and the spline-interpolated profile data around it
/// for peak patterns (source `MultiplexFilteringProfile`).
///
/// Searches the centroided data for the patterns of multiple peptides as
/// [`MultiplexFilteringCentroided`](super::filtering_centroided::MultiplexFilteringCentroided)
/// does, then scans the spline-interpolated profile of every candidate peak
/// from boundary to boundary and keeps the profile data points whose
/// satellites pass the averagine and peptide correlation filters. The shared
/// state and filters are those of [`MultiplexFiltering`], reachable through
/// [`Deref`].
#[derive(Clone, Debug)]
pub struct MultiplexFilteringProfile {
    base: MultiplexFiltering,
    /// Spline-interpolated profile spectra (source `exp_spline_profile_`).
    exp_spline_profile: Vec<SplineInterpolatedPeaks>,
    /// Peak boundaries of the peaks above the intensity cutoff (source
    /// `boundaries_`).
    boundaries: Vec<Vec<PeakBoundary>>,
}

impl Deref for MultiplexFilteringProfile {
    type Target = MultiplexFiltering;
    fn deref(&self) -> &MultiplexFiltering {
        &self.base
    }
}

impl DerefMut for MultiplexFilteringProfile {
    fn deref_mut(&mut self) -> &mut MultiplexFiltering {
        &mut self.base
    }
}

impl MultiplexFilteringProfile {
    /// Prepare filtering (source constructor).
    ///
    /// `exp_profile` is the profile data, `exp_centroided` the centroided data
    /// and `boundaries` the peak boundaries of `exp_centroided`, spectrum by
    /// spectrum and peak by peak, as `PeakPickerHiRes` reports them. The other
    /// arguments are those of [`MultiplexFiltering::new`].
    ///
    /// The boundaries of the peaks at or below the intensity cutoff are dropped,
    /// to keep the one-to-one mapping between peaks and boundaries, and each
    /// profile spectrum is spline-interpolated.
    ///
    /// The boundaries stay in the input order of `exp_centroided` and the
    /// splines in that of `exp_profile`, while the centroided data is sorted.
    /// The three line up only for input sorted by RT and m/z, which is what
    /// `FeatureFinderMultiplexAlgorithm` passes; the source has the same
    /// requirement and does not state it.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] when profile and centroided data do not hold
    ///   the same number of spectra, or the centroided data and its boundaries
    ///   (source `Exception::IllegalArgument`). The source checks the latter
    ///   only after it has indexed the boundaries; the port checks first.
    /// - [`Error::InvalidRange`] when a peak above the cutoff lies outside its
    ///   boundary (source `Exception::InvalidRange`), and
    ///   [`Error::InvalidValue`] when a spectrum has fewer boundaries than
    ///   peaks, which the source reads past.
    /// - The errors of [`MultiplexFiltering::new`] and of
    ///   [`SplineInterpolatedPeaks::from_spectrum`], which refuses a profile
    ///   spectrum of fewer than three points as the source does.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        exp_profile: &MSExperiment,
        exp_centroided: &MSExperiment,
        boundaries: &[Vec<PeakBoundary>],
        patterns: Vec<MultiplexIsotopicPeakPattern>,
        isotopes_per_peptide_min: usize,
        isotopes_per_peptide_max: usize,
        intensity_cutoff: f64,
        rt_band: f64,
        mz_tolerance: f64,
        mz_tolerance_unit_in_ppm: bool,
        peptide_similarity: f64,
        averagine_similarity: f64,
        averagine_similarity_scaling: f64,
        averagine_type: &str,
    ) -> Result<Self> {
        let base = MultiplexFiltering::new(
            exp_centroided,
            patterns,
            isotopes_per_peptide_min,
            isotopes_per_peptide_max,
            intensity_cutoff,
            rt_band,
            mz_tolerance,
            mz_tolerance_unit_in_ppm,
            peptide_similarity,
            averagine_similarity,
            averagine_similarity_scaling,
            averagine_type,
        )?;
        if exp_profile.len() != exp_centroided.len() {
            return Err(Error::InvalidValue(format!(
                "Profile and centroided data do not contain same number of spectra. ({}!={})",
                exp_profile.len(),
                exp_centroided.len()
            )));
        }
        if exp_centroided.len() != boundaries.len() {
            return Err(Error::InvalidValue(format!(
                "Centroided data and the corresponding list of peak boundaries do not contain same number of spectra. ({}!={})",
                exp_centroided.len(),
                boundaries.len()
            )));
        }
        let mut kept = Vec::with_capacity(boundaries.len());
        for (spectrum, spectrum_boundaries) in exp_centroided.spectra.iter().zip(boundaries) {
            let mut boundaries_temp = Vec::new();
            for (idx_mz, peak) in spectrum.peaks.iter().enumerate() {
                if f64::from(peak.intensity) > intensity_cutoff {
                    let boundary = *spectrum_boundaries.get(idx_mz).ok_or_else(|| {
                        Error::InvalidValue("a centroided peak has no peak boundary".into())
                    })?;
                    if boundary.min > peak.mz || peak.mz > boundary.max {
                        return Err(Error::InvalidRange(
                            "a centroided peak lies outside its peak boundary".into(),
                        ));
                    }
                    boundaries_temp.push(boundary);
                }
            }
            kept.push(boundaries_temp);
        }
        let exp_spline_profile = exp_profile
            .spectra
            .iter()
            .map(SplineInterpolatedPeaks::from_spectrum)
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            base,
            exp_spline_profile,
            boundaries: kept,
        })
    }

    /// The intensity-filtered peak boundaries (source `getPeakBoundaries`).
    pub fn peak_boundaries(&self) -> &[Vec<PeakBoundary>] {
        &self.boundaries
    }

    /// The spline-interpolated profile spectra (source `exp_spline_profile_`).
    pub fn spline_profile(&self) -> &[SplineInterpolatedPeaks] {
        &self.exp_spline_profile
    }

    /// Filter for the patterns: one result per pattern, in pattern order
    /// (source `filter`).
    ///
    /// # Errors
    ///
    /// As [`Self::filter_with_progress`].
    pub fn filter(&mut self) -> Result<Vec<MultiplexFilteredMSExperiment>> {
        self.filter_with_progress(&mut ProgressReporter::silent())
    }

    /// [`Self::filter`], reporting progress with the label "filtering LC-MS
    /// data" over `patterns * spectra` steps.
    ///
    /// A navigator with the default step scaling
    /// ([`DEFAULT_NAVIGATOR_SCALING`]) is made for every spline spectrum; the
    /// navigators live for the whole call, as in the source. For each pattern
    /// the white experiment is rebuilt, and every peak of every white spectrum
    /// that is not empty, has boundaries and a spline, is tried as a primary
    /// peak with [`MultiplexFiltering::filter_peak_positions`]. A candidate's
    /// profile is scanned from its lower to its upper boundary in the
    /// navigator's steps: at each position, every satellite is shifted by the
    /// same offset from its centroid and evaluated on its own spectrum's
    /// spline (RT and intensity narrowed to `f32`), and when these satellite
    /// data points pass
    /// [`MultiplexFiltering::filter_averagine_model_profile`] and
    /// [`MultiplexFiltering::filter_peptide_correlation_profile`] they are
    /// added to the peak. A peak with at least one such data point joins the
    /// result and its mass traces are blacklisted at once.
    ///
    /// As for the centroided filter, the source's `omp parallel for` over the
    /// peaks makes its blacklist racy; this is the serial (one-thread) result.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] "packages_ must not be empty" when a profile
    ///   spectrum has no spline package (source `Exception::InvalidSize` from
    ///   `getNavigator`).
    /// - [`Error::InvalidValue`] when the profile scan stops advancing, where
    ///   the source would loop forever.
    /// - The errors of the filters, the blacklisting and `progress`.
    pub fn filter_with_progress(
        &mut self,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<Vec<MultiplexFilteredMSExperiment>> {
        let Self {
            base,
            exp_spline_profile,
            boundaries,
        } = self;
        let steps = base.patterns.len().saturating_mul(exp_spline_profile.len());
        progress.start_count(steps, "filtering LC-MS data")?;
        let mut navigators = exp_spline_profile
            .iter()
            .map(|spline| spline.navigator(DEFAULT_NAVIGATOR_SCALING))
            .collect::<Result<Vec<_>>>()?;
        let out_of_range =
            || Error::InvalidValue("MultiplexFilteringProfile: index outside the data".into());

        let mut done = 0usize;
        let mut filter_results = Vec::with_capacity(base.patterns.len());
        for pattern_idx in 0..base.patterns.len() {
            let pattern = base.patterns[pattern_idx].clone();
            let mut result = MultiplexFilteredMSExperiment::new();
            base.update_white_ms_experiment();
            for idx_rt in 0..base.exp_centroided_white.spectra.len() {
                let white = &base.exp_centroided_white.spectra;
                let spectrum_boundaries = boundaries.get(idx_rt).ok_or_else(out_of_range)?;
                let spline_size = exp_spline_profile
                    .get(idx_rt)
                    .ok_or_else(out_of_range)?
                    .size();
                if white[idx_rt].peaks.is_empty()
                    || spectrum_boundaries.is_empty()
                    || spline_size == 0
                {
                    continue;
                }
                done += 1;
                progress.set_count(done)?;
                let rt = white[idx_rt].rt;
                let band = rt_begin(white, rt - base.rt_band / 2.0)
                    ..rt_end(white, rt + base.rt_band / 2.0);
                for s in 0..white[idx_rt].peaks.len() {
                    let mz = base.exp_centroided_white.spectra[idx_rt].peaks[s].mz;
                    let mz_idx = base.exp_centroided_mapping[idx_rt][s];
                    let mut peak = MultiplexFilteredPeak::new(mz, rt as f32, mz_idx, idx_rt);
                    if !base.filter_peak_positions(mz, band.clone(), &pattern, &mut peak)? {
                        continue;
                    }
                    let boundary = *boundaries[idx_rt].get(mz_idx).ok_or_else(out_of_range)?;
                    scan_profile(base, &mut navigators, &pattern, &mut peak, idx_rt, boundary)?;
                    if peak.satellite_profile_count() > 0 {
                        base.blacklist_peak(&peak, pattern_idx)?;
                        result.add_peak(peak);
                    }
                }
            }
            filter_results.push(result);
        }
        progress.end()?;
        Ok(filter_results)
    }
}

/// The profile scan of one candidate peak between its boundaries; see
/// [`MultiplexFilteringProfile::filter_with_progress`].
fn scan_profile(
    base: &MultiplexFiltering,
    navigators: &mut [crate::processing::spline::spline_interpolated_peaks::Navigator<'_>],
    pattern: &MultiplexIsotopicPeakPattern,
    peak: &mut MultiplexFilteredPeak,
    idx_rt: usize,
    boundary: PeakBoundary,
) -> Result<()> {
    let out_of_range =
        || Error::InvalidValue("MultiplexFilteringProfile: index outside the data".into());
    let peak_mz = peak.mz();
    let satellites: Vec<(usize, usize, usize)> = peak
        .satellite_entries()
        .map(|(key, s)| (key, s.rt_idx, s.mz_idx))
        .collect();
    let mut mz_profile = boundary.min;
    while mz_profile < boundary.max {
        let mz_shift = mz_profile - peak_mz;
        let mut satellites_profile: SatelliteMultimap<MultiplexSatelliteProfile> =
            SatelliteMultimap::new();
        for &(key, rt_idx, mz_idx) in &satellites {
            let spectrum = base
                .exp_centroided
                .spectra
                .get(rt_idx)
                .ok_or_else(out_of_range)?;
            let mz_satellite = spectrum.peaks.get(mz_idx).ok_or_else(out_of_range)?.mz;
            let mz = mz_satellite + mz_shift;
            let intensity = navigators
                .get_mut(rt_idx)
                .ok_or_else(out_of_range)?
                .eval(mz)?;
            satellites_profile
                .entry(key)
                .or_default()
                .push(MultiplexSatelliteProfile::new(
                    spectrum.rt as f32,
                    mz,
                    intensity as f32,
                ));
        }
        if base.filter_averagine_model_profile(pattern, peak, &satellites_profile)?
            && base.filter_peptide_correlation_profile(pattern, &satellites_profile)?
        {
            for (key, points) in satellites_profile {
                for point in points {
                    peak.add_satellite_profile_from(point, key);
                }
            }
        }
        let next = navigators
            .get_mut(idx_rt)
            .ok_or_else(out_of_range)?
            .next_pos(mz_profile)?;
        if next <= mz_profile && next < boundary.max {
            return Err(Error::InvalidValue(
                "MultiplexFilteringProfile: the profile scan does not advance (the source loops forever)".into(),
            ));
        }
        mz_profile = next;
    }
    Ok(())
}
