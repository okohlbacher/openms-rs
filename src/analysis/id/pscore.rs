// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `ANALYSIS/ID/PScore.h`: the Andromeda-style PScore, windowed
//! intensity ranks and peak-level spectra.
//!
//! The PScore counts how many theoretical fragment ions are matched in an
//! experimental spectrum that has been filtered to its top `level + 1` peaks
//! per m/z window, and turns the count into `-10 log10` of the cumulative
//! binomial probability of matching at least that many ions by chance
//! ([`AScore::compute_cumulative_score`]). The best score over all peak levels
//! (and, for [`PScore::compute_pscore_best_of`], over all theoretical spectra)
//! is returned.
//!
//! The source class has only static members; [`PScore`] is a unit struct that
//! carries them as associated functions. The mapping, the preserved source
//! conventions and the source defects are recorded in
//! `docs/PSCORE_ASCORE_SUPPORT.md`.

use crate::analysis::id::ascore::AScore;
use crate::comparison::{Tolerance, matched_alignment};
use crate::kernel::MSExperiment;
use crate::{Error, MSSpectrum, Result};
use std::collections::BTreeMap;

/// Largest peak count [`PScore::calculate_intensity_rank_in_mz_window`] ranks.
///
/// The rank vector is allocated up front and the neighbour scan is quadratic
/// in the number of peaks inside one window, so a ceiling bounds both. The
/// source has none. Checked before anything is allocated.
pub const MAX_RANK_PEAKS: usize = 10_000_000;

/// Default m/z window, in Th, of the rank and PScore functions (the source's
/// default argument `mz_window = 100`).
pub const DEFAULT_MZ_WINDOW: f64 = 100.0;

/// Default lowest peak level of [`PScore::calculate_peak_level_spectra`]
/// (source default `min_level = 1`).
pub const DEFAULT_MIN_LEVEL: usize = 1;

/// Default highest peak level of [`PScore::calculate_peak_level_spectra`]
/// (source default `max_level = 9`).
pub const DEFAULT_MAX_LEVEL: usize = 9;

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

/// The PScore of Andromeda (Cox et al., J. Proteome Res. 2011) and the
/// windowed peak ranks it is computed from.
///
/// Every source member is static; the source's implicit constructor and
/// destructor are the `Default`, `Clone` and drop of this unit struct.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PScore;

impl PScore {
    /// Calculate local (windowed) peak ranks.
    ///
    /// The peak rank is the number of neighbouring peaks within
    /// `+/- mz_window / 2` Th that have a strictly higher intensity. The result
    /// can be used to filter a spectrum efficiently for its top `1..n` peaks in
    /// mass windows. Ranks are zero based: the highest-intensity peak in its
    /// window has rank 0.
    ///
    /// # Arguments
    ///
    /// * `mz` - m/z positions of the peaks, ascending
    /// * `intensities` - intensities of the peaks, one per position
    /// * `mz_window` - window in Th, centred at each peak
    ///
    /// The neighbour scan stops at the first neighbour outside the window on
    /// each side, exactly as the source does, so it relies on `mz` being
    /// ascending. The source does not check that and neither does this port:
    /// unsorted positions give the same (meaningless) ranks the source gives.
    /// An empty input returns an empty vector.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the two slices differ in length
    /// (the source reads past the end of `intensities` when it is shorter), when
    /// any position, intensity or the window is not finite (the source's
    /// comparisons silently misbehave on NaN), or when there are more than
    /// [`MAX_RANK_PEAKS`] peaks.
    pub fn calculate_intensity_rank_in_mz_window(
        mz: &[f64],
        intensities: &[f64],
        mz_window: f64,
    ) -> Result<Vec<usize>> {
        if mz.len() != intensities.len() {
            return Err(bad(
                "PScore ranks need exactly one intensity per m/z position",
            ));
        }
        if mz.len() > MAX_RANK_PEAKS {
            return Err(bad("PScore rank input exceeds MAX_RANK_PEAKS"));
        }
        if !mz_window.is_finite() {
            return Err(bad("PScore m/z window must be finite"));
        }
        if mz.iter().chain(intensities).any(|value| !value.is_finite()) {
            return Err(bad("PScore rank input must be finite"));
        }
        let half_window = mz_window / 2.0;
        let mut ranks = Vec::with_capacity(mz.len());
        for (p, (&m, &i)) in mz.iter().zip(intensities).enumerate() {
            let mut rank = 0;
            // Neighbours to the left with higher intensity.
            for j in (0..p).rev() {
                if mz[j] < m - half_window {
                    break;
                }
                if intensities[j] > i {
                    rank += 1;
                }
            }
            // Neighbours to the right with higher intensity.
            for j in p + 1..mz.len() {
                if mz[j] > m + half_window {
                    break;
                }
                if intensities[j] > i {
                    rank += 1;
                }
            }
            ranks.push(rank);
        }
        Ok(ranks)
    }

    /// Precalculated windowed peak ranks for a whole experiment.
    ///
    /// For every spectrum of `peak_map`, in storage order, the ranks of
    /// [`PScore::calculate_intensity_rank_in_mz_window`] over its peaks'
    /// m/z and intensity. Ranks are zero based (the top peak has rank 0).
    /// `peak_map` is typically an experiment from which all MS1 spectra were
    /// removed; the source default for `mz_window` is [`DEFAULT_MZ_WINDOW`].
    ///
    /// # Errors
    ///
    /// As [`PScore::calculate_intensity_rank_in_mz_window`], for the first
    /// spectrum that fails.
    pub fn calculate_rank_map(peak_map: &MSExperiment, mz_window: f64) -> Result<Vec<Vec<usize>>> {
        let mut rank_map = Vec::with_capacity(peak_map.spectra.len());
        for spectrum in &peak_map.spectra {
            let mz: Vec<f64> = spectrum.peaks.iter().map(|peak| peak.mz).collect();
            let intensities: Vec<f64> = spectrum
                .peaks
                .iter()
                .map(|peak| f64::from(peak.intensity))
                .collect();
            rank_map.push(Self::calculate_intensity_rank_in_mz_window(
                &mz,
                &intensities,
                mz_window,
            )?);
        }
        Ok(rank_map)
    }

    /// Spectra for the peak levels `min_level..=max_level`, keyed by level.
    ///
    /// The spectrum of level `n` keeps the peaks whose windowed rank (from
    /// [`PScore::calculate_intensity_rank_in_mz_window`]) is at most `n`, that
    /// is the `n + 1` most intense peaks in the window centred at each peak.
    /// Levels are zero based: level 0 keeps only the top peak of each window,
    /// level 1 the top two. The source notes that its default range,
    /// [`DEFAULT_MIN_LEVEL`]`..=`[`DEFAULT_MAX_LEVEL`], comes from the
    /// Andromeda publication and is similar to the AScore one.
    ///
    /// A level only appears in the map once a peak is added to it, as the
    /// source's `operator[]` creates it on the first push: an empty `spec`, or
    /// `min_level > max_level`, gives an empty map. Each level spectrum carries
    /// the kept peaks in the input order and nothing else of `spec` (the source
    /// `push_back` copies the peak only).
    ///
    /// `ranks` is walked, not `spec`: a shorter `ranks` processes only the
    /// leading peaks, as in the source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `ranks` is longer than `spec` (the
    /// source reads past the end of the spectrum), or when a level or a rank
    /// exceeds `i32::MAX` (the source compares them after an `int` cast).
    pub fn calculate_peak_level_spectra(
        spec: &MSSpectrum,
        ranks: &[usize],
        min_level: usize,
        max_level: usize,
    ) -> Result<BTreeMap<usize, MSSpectrum>> {
        let mut peak_level_spectra: BTreeMap<usize, MSSpectrum> = BTreeMap::new();
        if spec.is_empty() {
            return Ok(peak_level_spectra);
        }
        if ranks.len() > spec.len() {
            return Err(bad(
                "PScore peak levels need a spectrum peak for every rank",
            ));
        }
        let limit = i32::MAX as usize;
        if min_level > limit || max_level > limit || ranks.iter().any(|&rank| rank > limit) {
            return Err(bad("PScore peak levels and ranks must fit the source int"));
        }
        for (peak, &rank) in spec.peaks.iter().zip(ranks) {
            // Start at the highest (least restrictive) level; a rank too high
            // for one level is too high for every lower level.
            let mut level = max_level;
            while level >= min_level {
                if rank > level {
                    break;
                }
                peak_level_spectra
                    .entry(level)
                    .or_default()
                    .peaks
                    .push(*peak);
                if level == 0 {
                    break;
                }
                level -= 1;
            }
        }
        Ok(peak_level_spectra)
    }

    /// The best PScore over several theoretical spectra (source
    /// `computePScore` overload taking a `std::vector<PeakSpectrum>`).
    ///
    /// As in Andromeda, the theoretical spectra can hold e.g. loss spectra or
    /// higher-charge spectra depending on the sequence; the best score over
    /// all of them and over all peak levels is returned, starting from `0.0`.
    ///
    /// # Arguments
    ///
    /// * `fragment_mass_tolerance` - tolerance for matching peaks
    /// * `fragment_mass_tolerance_unit_ppm` - whether the tolerance is in ppm
    ///   of the theoretical m/z (otherwise Th)
    /// * `peak_level_spectra` - spectra filtered by maximum rank, as from
    ///   [`PScore::calculate_peak_level_spectra`]
    /// * `theo_spectra` - theoretical spectra, e.g. from
    ///   [`TheoreticalSpectrumGenerator`](crate::chemistry::theoretical::TheoreticalSpectrumGenerator)
    /// * `mz_window` - window in Th centred at each peak; source default
    ///   [`DEFAULT_MZ_WINDOW`]
    ///
    /// This overload keeps the source's own matching, which differs from
    /// [`PScore::compute_pscore`] in two ways, both preserved:
    ///
    /// * each theoretical peak counts as matched when the *nearest*
    ///   experimental peak lies strictly closer than the tolerance
    ///   (`findNearest` and a `<` test in double precision), so two
    ///   theoretical peaks can share one experimental peak;
    /// * the match probability is `level / mz_window`, not
    ///   `(level + 1) / mz_window`. A level-0 spectrum therefore has `p = 0`,
    ///   and any match at that level scores `+inf`; no match scores `0`.
    ///
    /// The second is a source defect candidate (the single-spectrum overload
    /// and the AScore both use the `+ 1` form); it is kept for parity.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] when a peak-level spectrum is not sorted
    /// by m/z, and [`Error::InvalidValue`] when one is empty (the source's
    /// `findNearest` throws `Exception::Precondition`), when the tolerance is
    /// not finite, when a match probability exceeds one (the source's binomial
    /// throws `std::invalid_argument`), or when a count exceeds `u32::MAX`.
    pub fn compute_pscore_best_of(
        fragment_mass_tolerance: f64,
        fragment_mass_tolerance_unit_ppm: bool,
        peak_level_spectra: &BTreeMap<usize, MSSpectrum>,
        theo_spectra: &[MSSpectrum],
        mz_window: f64,
    ) -> Result<f64> {
        if !fragment_mass_tolerance.is_finite() {
            return Err(bad("PScore fragment tolerance must be finite"));
        }
        let mut best_pscore = 0.0;
        for theo_spectrum in theo_spectra {
            let n_trials = theo_spectrum.len();
            for (&level, exp_spectrum) in peak_level_spectra {
                let mut matched_peaks = 0;
                for theo_peak in &theo_spectrum.peaks {
                    let theo_mz = theo_peak.mz;
                    let max_dist_dalton = if fragment_mass_tolerance_unit_ppm {
                        theo_mz * fragment_mass_tolerance * 1e-6
                    } else {
                        fragment_mass_tolerance
                    };
                    let index = exp_spectrum.find_nearest(theo_mz)?.ok_or_else(|| {
                        bad("PScore peak-level spectrum must not be empty (source Precondition)")
                    })?;
                    let exp_mz = exp_spectrum.peaks[index].mz;
                    if (theo_mz - exp_mz).abs() < max_dist_dalton {
                        matched_peaks += 1;
                    }
                }
                let p = level as f64 / mz_window;
                let pscore =
                    -10.0 * AScore::compute_cumulative_score(n_trials, matched_peaks, p)?.log10();
                if pscore > best_pscore {
                    best_pscore = pscore;
                }
            }
        }
        Ok(best_pscore)
    }

    /// The PScore of one theoretical spectrum (source `computePScore`
    /// overload taking a single `PeakSpectrum`).
    ///
    /// For every peak level, ascending, the theoretical peaks matched in the
    /// level spectrum are counted with the source's `MatchedIterator` walk
    /// (`PpmTrait` or `DaTrait`: each theoretical peak is paired with the
    /// closest experimental peak the forward merge reaches, within an
    /// inclusive tolerance computed in single precision), and scored as
    /// `-10 log10 P(X >= matched)` for `X ~ B(N, (level + 1) / mz_window)`
    /// with `N` the theoretical peak count. The best score over the levels is
    /// returned, starting from `0.0`.
    ///
    /// # Arguments
    ///
    /// * `fragment_mass_tolerance` - tolerance for matching peaks
    /// * `fragment_mass_tolerance_unit_ppm` - whether the tolerance is in ppm
    ///   of the theoretical m/z (otherwise Th)
    /// * `peak_level_spectra` - spectra filtered by maximum rank
    /// * `theo_spectrum` - a theoretical spectrum, sorted by m/z
    /// * `mz_window` - window in Th centred at each peak; source default
    ///   [`DEFAULT_MZ_WINDOW`]
    ///
    /// A cumulative probability of zero (every ion matched at a vanishing
    /// probability) scores `+inf`, as in the source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] when the theoretical spectrum or a level
    /// spectrum is not sorted by m/z; the source's `MatchedIterator` documents
    /// sorted input as a requirement and does not check it. Returns
    /// [`Error::InvalidValue`] when the tolerance is negative, not finite or
    /// outside the `f32` range the source's iterator stores it in, when a match
    /// probability exceeds one (the source's binomial throws
    /// `std::invalid_argument`), or when a count exceeds `u32::MAX`.
    pub fn compute_pscore(
        fragment_mass_tolerance: f64,
        fragment_mass_tolerance_unit_ppm: bool,
        peak_level_spectra: &BTreeMap<usize, MSSpectrum>,
        theo_spectrum: &MSSpectrum,
        mz_window: f64,
    ) -> Result<f64> {
        let tolerance = if fragment_mass_tolerance_unit_ppm {
            Tolerance::Ppm(fragment_mass_tolerance)
        } else {
            Tolerance::Absolute(fragment_mass_tolerance)
        };
        if !theo_spectrum.is_sorted() {
            return Err(Error::UnsortedData);
        }
        let n_trials = theo_spectrum.len();
        let mut best_pscore = 0.0;
        for (&level, exp_spectrum) in peak_level_spectra {
            if !exp_spectrum.is_sorted() {
                return Err(Error::UnsortedData);
            }
            let matched_peaks = matched_alignment(theo_spectrum, exp_spectrum, tolerance)?.len();
            let p = (level as f64 + 1.0) / mz_window;
            let pscore =
                -10.0 * AScore::compute_cumulative_score(n_trials, matched_peaks, p)?.log10();
            if pscore > best_pscore {
                best_pscore = pscore;
            }
        }
        Ok(best_pscore)
    }

    /// Andromeda's additive mass correction, `0.024 * (mass - 600)`.
    ///
    /// One of the correction terms Andromeda adds to the PScore
    /// (`pscore + massC + cleaveC + modC - 100`; see the Andromeda source).
    /// The source notes that the constants might be instrument dependent.
    ///
    /// The source declares this as a static member but defines a free
    /// function of the same name in `PScore.cpp`, so `PScore::massCorrectionTerm`
    /// has no definition and a call fails to link. The port provides the body
    /// the source wrote.
    pub fn mass_correction_term(mass: f64) -> f64 {
        0.024 * (mass - 600.0)
    }

    /// Andromeda's additive correction for the type of cleavage.
    ///
    /// `53.2` for no cleavage, `42.1` for one consecutive cleavage and `31.1`
    /// for one non-consecutive cleavage, `17.0` for two, `0.0` beyond. The
    /// source notes that the constants might be instrument dependent, and has
    /// the same missing-definition defect as [`PScore::mass_correction_term`].
    pub fn cleavage_correction_term(cleavages: usize, consecutive_cleavage: bool) -> f64 {
        match cleavages {
            0 => 53.2,
            1 if consecutive_cleavage => 42.1,
            1 => 31.1,
            2 => 17.0,
            _ => 0.0,
        }
    }

    /// Andromeda's additive correction for the number of modifications.
    ///
    /// `42, 28, 22, 16, 9, 5, 2` for zero to six modifications and `0.0`
    /// beyond. The source notes that the constants might be instrument
    /// dependent, and has the same missing-definition defect as
    /// [`PScore::mass_correction_term`].
    pub fn modification_correction_term(modifications: usize) -> f64 {
        match modifications {
            0 => 42.0,
            1 => 28.0,
            2 => 22.0,
            3 => 16.0,
            4 => 9.0,
            5 => 5.0,
            6 => 2.0,
            _ => 0.0,
        }
    }
}
