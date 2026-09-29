// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Code common to filtering centroided and profile data for multiplex peak
//! patterns.
//!
//! Native equivalent of `FEATUREFINDER/MultiplexFiltering.h`. See
//! `docs/FEATUREFINDER_MULTIPLEX_SUPPORT.md`.

use std::collections::BTreeMap;
use std::ops::Range;

use super::filtered_peak::{MultiplexFilteredPeak, SatelliteMultimap};
use super::isotopic_peak_pattern::MultiplexIsotopicPeakPattern;
use super::satellite::MultiplexSatelliteProfile;
use crate::chemistry::isotopes::ProbabilityPrecision;
use crate::chemistry::{CoarseIsotopePatternGenerator, CoarseMassMode, IsotopeDistribution};
use crate::kernel::{MSExperiment, MSSpectrum, Peak1D};
use crate::math::source_sort::source_sort_by;
use crate::math::statistic_functions::{
    pearson_correlation_coefficient, rank_correlation_coefficient,
};
use crate::{Error, Result};

/// Ceiling on the peaks above the intensity cutoff, checked before the
/// filtered experiment, the blacklist and the white experiment (three copies
/// of every retained peak) are allocated. The source has no bound.
pub const MAX_FILTERED_PEAKS: usize = 200_000_000;

/// Index mapping from a white spectrum's peaks to the original spectrum's
/// (source `MultiplexFiltering::White2Original`, a
/// `std::vector<std::map<int, int>>`).
///
/// A white experiment is the centroided experiment without its blacklisted
/// peaks, so its peak indices are shifted; entry `i` of spectrum `s` is the
/// original index of white peak `i`. The source's map is dense (its keys are
/// `0..count`), so a vector is the same mapping.
pub type White2Original = Vec<Vec<usize>>;

/// The averagine model used by the averagine filter (source parameter
/// `averagine_type`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AveragineType {
    /// `peptide`: `CoarseIsotopePatternGenerator::estimateFromPeptideWeight`.
    Peptide,
    /// `RNA`: `estimateFromRNAWeight`.
    Rna,
    /// `DNA`: `estimateFromDNAWeight`.
    Dna,
}

impl AveragineType {
    /// The type named `name` (`peptide`, `RNA` or `DNA`, case-sensitive), or
    /// `None`.
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "peptide" => Some(Self::Peptide),
            "RNA" => Some(Self::Rna),
            "DNA" => Some(Self::Dna),
            _ => None,
        }
    }
}

/// Source `MSSpectrum::findNearest(mz)` on sorted, non-empty `peaks`: the
/// lower bound of `mz`, or the peak before it when that one is not farther
/// away (ties choose the lower m/z) (`MSSpectrum.cpp:329-358`).
pub(crate) fn find_nearest(peaks: &[Peak1D], mz: f64) -> Option<usize> {
    if peaks.is_empty() {
        return None;
    }
    let it = peaks.partition_point(|p| p.mz < mz);
    if it == 0 {
        return Some(0);
    }
    if it == peaks.len() {
        return Some(peaks.len() - 1);
    }
    if (peaks[it].mz - mz).abs() < (peaks[it - 1].mz - mz).abs() {
        Some(it)
    } else {
        Some(it - 1)
    }
}

/// Source `MSSpectrum::findNearest(mz, tolerance)` (`MSSpectrum.cpp:311-327`):
/// the nearest peak when it lies within `[mz - tolerance, mz + tolerance]`,
/// otherwise `None` (the source's `-1`). Unlike
/// [`MSSpectrum::find_nearest_with_tolerance`] it does not try the neighbour
/// on the other side, and it does not re-check sortedness on every call.
pub(crate) fn find_nearest_with_tolerance(
    peaks: &[Peak1D],
    mz: f64,
    tolerance: f64,
) -> Option<usize> {
    let i = find_nearest(peaks, mz)?;
    let found = peaks[i].mz;
    if found >= mz - tolerance && found <= mz + tolerance {
        Some(i)
    } else {
        None
    }
}

/// Source `MSExperiment::RTBegin(rt)` on RT-sorted spectra: the first spectrum
/// whose RT is not less than `rt`.
pub(crate) fn rt_begin(spectra: &[MSSpectrum], rt: f64) -> usize {
    spectra.partition_point(|s| s.rt < rt)
}

/// Source `MSExperiment::RTEnd(rt)` on RT-sorted spectra: the first spectrum
/// whose RT is greater than `rt`.
pub(crate) fn rt_end(spectra: &[MSSpectrum], rt: f64) -> usize {
    spectra.partition_point(|s| !(rt < s.rt))
}

/// The pattern's m/z shift at `i`, or an error where the source reads past the
/// shift list.
fn mz_shift(pattern: &MultiplexIsotopicPeakPattern, i: usize) -> Result<f64> {
    pattern.mz_shift_at(i).ok_or_else(|| {
        Error::InvalidValue(format!(
            "MultiplexFiltering: m/z shift {i} is outside the pattern's {} shifts (the source reads past the end); \
             patterns need isotopes_per_peptide_max peaks per peptide",
            pattern.mz_shift_count()
        ))
    })
}

/// Code common to filtering centroided and profile data for peak patterns
/// (source `MultiplexFiltering`).
///
/// The algorithm searches the data for patterns of multiple peptides, which
/// appear as characteristic patterns of isotopic peaks in MS1 spectra. It first
/// searches the centroided data, and optionally, in a second step, the
/// spline-interpolated profile data; each peak pattern yields one filter
/// result. The two filters
/// [`MultiplexFilteringCentroided`](super::filtering_centroided::MultiplexFilteringCentroided)
/// and [`MultiplexFilteringProfile`](super::filtering_profile::MultiplexFilteringProfile)
/// wrap this type and contain the mode-specific `filter` method.
///
/// The source derives from `ProgressLogger`; the filters take a
/// [`ProgressReporter`](crate::concept::progress_logger::ProgressReporter)
/// instead.
///
/// The source declares a single-argument `blacklistPeak_(const
/// MultiplexFilteredPeak&)` that it never defines; it cannot be called and is
/// not ported.
#[derive(Clone, Debug)]
pub struct MultiplexFiltering {
    /// Centroided data above the intensity cutoff, sorted (source
    /// `exp_centroided_`).
    pub(crate) exp_centroided: MSExperiment,
    /// Per spectrum and peak of `exp_centroided`: `-1` for white, otherwise the
    /// index of the mass trace in the pattern that blacklisted the peak (source
    /// `blacklist_`).
    pub(crate) blacklist: Vec<Vec<i32>>,
    /// `exp_centroided` without its blacklisted peaks (source
    /// `exp_centroided_white_`).
    pub(crate) exp_centroided_white: MSExperiment,
    /// White to original peak indices (source `exp_centroided_mapping_`).
    pub(crate) exp_centroided_mapping: White2Original,
    /// The peak patterns (source `patterns_`).
    pub(crate) patterns: Vec<MultiplexIsotopicPeakPattern>,
    /// Minimum number of isotopic peaks per peptide (source
    /// `isotopes_per_peptide_min_`).
    pub(crate) isotopes_per_peptide_min: usize,
    /// Maximum number of isotopic peaks per peptide (source
    /// `isotopes_per_peptide_max_`).
    pub(crate) isotopes_per_peptide_max: usize,
    /// Intensity cutoff (source `intensity_cutoff_`).
    pub(crate) intensity_cutoff: f64,
    /// RT range used for filtering (source `rt_band_`).
    pub(crate) rt_band: f64,
    /// m/z shift tolerance (source `mz_tolerance_`).
    pub(crate) mz_tolerance: f64,
    /// Unit of the tolerance: ppm when true, Th when false (source
    /// `mz_tolerance_unit_in_ppm_`).
    pub(crate) mz_tolerance_unit_in_ppm: bool,
    /// Peptide similarity (source `peptide_similarity_`).
    pub(crate) peptide_similarity: f64,
    /// Averagine similarity (source `averagine_similarity_`).
    pub(crate) averagine_similarity: f64,
    /// Averagine similarity scaling (source `averagine_similarity_scaling_`).
    pub(crate) averagine_similarity_scaling: f64,
    /// Type of averagine to use (source `averagine_type_`).
    pub(crate) averagine_type: String,
}

impl MultiplexFiltering {
    /// Prepare filtering `exp_centroided` for `patterns` (source constructor).
    ///
    /// Peaks at or below `intensity_cutoff` cannot be relevant and are dropped
    /// right away, which reduces memory and run time. Each retained spectrum
    /// keeps only its RT and its peaks (all other spectrum metadata is dropped,
    /// and its MS level is the default 1, as the source's fresh `MSSpectrum`).
    /// The spectra are then sorted by RT with the Release build's `std::sort`
    /// ([`source_sort_by`]) and each spectrum's peaks by m/z with a stable sort,
    /// as `MSExperiment::sortSpectra(true)` does, and every peak starts white.
    ///
    /// The arguments: `isotopes_per_peptide_min` and
    /// `isotopes_per_peptide_max` bound the isotopic peaks of a peptide;
    /// `rt_band` is the RT range used for filtering; `mz_tolerance` is the error
    /// margin for matching expected patterns to the data, in ppm when
    /// `mz_tolerance_unit_in_ppm` is true and in Th otherwise;
    /// `peptide_similarity` is the similarity required between two peptides of
    /// one multiplet; `averagine_similarity` the similarity required between a
    /// peptide's isotope pattern and the averagine model; and
    /// `averagine_similarity_scaling` the factor `x` that makes the averagine
    /// similarity `p` stricter for singlets, `p' = p + x(1 - p)`.
    /// `averagine_type` names the averagine model, `peptide`, `RNA` or `DNA`;
    /// as in the source it is checked only when the averagine filter first runs.
    ///
    /// The source takes the isotope counts as `int` and stores them as
    /// `size_t`, so a negative count wraps; the port takes `usize`.
    ///
    /// Every pattern must hold `isotopes_per_peptide_max` peaks per peptide:
    /// the filters address the shift of isotope `j` of peptide `i` as
    /// `i * isotopes_per_peptide_max + j`. The source does not check this and
    /// reads past the shift list otherwise; the port returns an error where the
    /// read would happen.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a non-finite retention time or peak
    /// m/z in `exp_centroided`, where the source's sorts and binary searches
    /// have no defined result, and [`Error::InvalidRange`] when more than
    /// [`MAX_FILTERED_PEAKS`] peaks lie above the cutoff.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        exp_centroided: &MSExperiment,
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
        let mut retained = 0usize;
        for spectrum in &exp_centroided.spectra {
            if !spectrum.rt.is_finite() {
                return Err(Error::InvalidValue(
                    "MultiplexFiltering: spectrum retention time is not finite".into(),
                ));
            }
            for peak in &spectrum.peaks {
                if f64::from(peak.intensity) > intensity_cutoff {
                    if !peak.mz.is_finite() {
                        return Err(Error::InvalidValue(
                            "MultiplexFiltering: peak m/z is not finite".into(),
                        ));
                    }
                    retained += 1;
                }
            }
        }
        if retained > MAX_FILTERED_PEAKS {
            return Err(Error::InvalidRange(format!(
                "MultiplexFiltering: more than {MAX_FILTERED_PEAKS} peaks above the intensity cutoff"
            )));
        }
        let mut spectra = Vec::with_capacity(exp_centroided.spectra.len());
        for spectrum in &exp_centroided.spectra {
            let mut filtered = MSSpectrum {
                rt: spectrum.rt,
                ..MSSpectrum::default()
            };
            filtered.peaks = spectrum
                .peaks
                .iter()
                .filter(|p| f64::from(p.intensity) > intensity_cutoff)
                .copied()
                .collect();
            spectra.push(filtered);
        }
        source_sort_by(&mut spectra, |a, b| a.rt < b.rt)?;
        for spectrum in &mut spectra {
            // `sortByPosition` is a `std::stable_sort` by m/z; every m/z is finite.
            spectrum.peaks.sort_by(|a, b| a.mz.total_cmp(&b.mz));
        }
        let blacklist = spectra.iter().map(|s| vec![-1; s.peaks.len()]).collect();
        let exp = MSExperiment {
            spectra,
            ..MSExperiment::default()
        };
        Ok(Self {
            exp_centroided: exp,
            blacklist,
            exp_centroided_white: MSExperiment::default(),
            exp_centroided_mapping: Vec::new(),
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
            averagine_type: averagine_type.to_owned(),
        })
    }

    /// The intensity-filtered, sorted centroided data (source
    /// `getCentroidedExperiment`).
    pub fn centroided_experiment(&self) -> &MSExperiment {
        &self.exp_centroided
    }

    /// The blacklisted peaks: every spectrum of the centroided data with only
    /// its peaks that are not white (source `getBlacklist`).
    pub fn blacklist_experiment(&self) -> MSExperiment {
        let spectra = self
            .exp_centroided
            .spectra
            .iter()
            .zip(&self.blacklist)
            .map(|(spectrum, marks)| MSSpectrum {
                rt: spectrum.rt,
                peaks: spectrum
                    .peaks
                    .iter()
                    .zip(marks)
                    .filter(|&(_, &mark)| mark != -1)
                    .map(|(peak, _)| *peak)
                    .collect(),
                ..MSSpectrum::default()
            })
            .collect();
        MSExperiment {
            spectra,
            ..MSExperiment::default()
        }
    }

    /// The intensity cutoff (source `intensity_cutoff_`).
    pub fn intensity_cutoff(&self) -> f64 {
        self.intensity_cutoff
    }

    /// The blacklist marks: per spectrum and peak, `-1` for white or the mass
    /// trace index that blacklisted it (source `blacklist_`).
    pub fn blacklist(&self) -> &[Vec<i32>] {
        &self.blacklist
    }

    /// The current white experiment (source `exp_centroided_white_`).
    pub fn white_experiment(&self) -> &MSExperiment {
        &self.exp_centroided_white
    }

    /// The current white-to-original index mapping (source
    /// `exp_centroided_mapping_`).
    pub fn white_mapping(&self) -> &White2Original {
        &self.exp_centroided_mapping
    }

    /// The peak patterns (source `patterns_`).
    pub fn patterns(&self) -> &[MultiplexIsotopicPeakPattern] {
        &self.patterns
    }

    /// Rebuild the white experiment from the peaks that are not blacklisted,
    /// and the mapping of white peak positions to their original positions
    /// (source protected `updateWhiteMSExperiment_`).
    pub fn update_white_ms_experiment(&mut self) {
        let mut white = Vec::with_capacity(self.exp_centroided.spectra.len());
        let mut mapping = Vec::with_capacity(self.exp_centroided.spectra.len());
        for (spectrum, marks) in self.exp_centroided.spectra.iter().zip(&self.blacklist) {
            let mut white_spectrum = MSSpectrum {
                rt: spectrum.rt,
                ..MSSpectrum::default()
            };
            let mut map = Vec::new();
            for (i, (peak, &mark)) in spectrum.peaks.iter().zip(marks).enumerate() {
                if mark == -1 {
                    white_spectrum.peaks.push(*peak);
                    map.push(i);
                }
            }
            white.push(white_spectrum);
            mapping.push(map);
        }
        self.exp_centroided_white = MSExperiment {
            spectra: white,
            ..MSExperiment::default()
        };
        self.exp_centroided_mapping = mapping;
    }

    /// The absolute m/z tolerance in Th at `mz`: `mz * tolerance * 1e-6` in
    /// ppm mode, the tolerance itself otherwise. The source computes it once
    /// per pattern, as it varies minimally within one.
    pub(crate) fn absolute_tolerance(&self, mz: f64) -> f64 {
        if self.mz_tolerance_unit_in_ppm {
            mz * self.mz_tolerance * 1e-6
        } else {
            self.mz_tolerance
        }
    }

    /// Whether a significant peak lies at `mz` in white spectrum `rt_idx`
    /// (source protected `checkForSignificantPeak_`).
    ///
    /// Returns the white peak index of the nearest peak within `mz_tolerance`
    /// when its intensity exceeds `0.3 * intensity_first_peak`, otherwise
    /// `None` (the source's `-1`).
    pub fn check_for_significant_peak(
        &self,
        mz: f64,
        mz_tolerance: f64,
        rt_idx: usize,
        intensity_first_peak: f64,
    ) -> Option<usize> {
        let peaks = &self.exp_centroided_white.spectra.get(rt_idx)?.peaks;
        let mz_idx = find_nearest_with_tolerance(peaks, mz, mz_tolerance)?;
        let threshold = 0.3;
        if f64::from(peaks[mz_idx].intensity) > threshold * intensity_first_peak {
            Some(mz_idx)
        } else {
            None
        }
    }

    /// Whether enough peaks in the RT band form `pattern` around the primary
    /// peak at `mz` (source protected `filterPeakPositions_`).
    ///
    /// `band` is the range of white spectra in the RT band. The primary peak's
    /// position must not be blacklisted, except as the mono-isotopic peak of the
    /// lightest (or only) peptide (mark `0`). Then, isotope by isotope, every
    /// peptide must show a peak at `mz` plus the pattern's shift in some
    /// spectrum of the band, a peak that is white or already seen as part of the
    /// same mass trace and not yet a satellite; each such peak is added to
    /// `peak` as a satellite. At least `isotopes_per_peptide_min` *consecutive*
    /// isotopes starting with the mono-isotopic one must be present in all
    /// peptides, since skipping one would make, say, a 2+ peptide look like a
    /// 4+ one.
    ///
    /// Finally, with twice the m/z tolerance, the primary peak's spectrum must
    /// hold no significant peak (above 30% of the peptide's first peak) left of
    /// the mono-isotopic peak that is not itself a satellite, and no peak that
    /// makes the pattern a higher charge state: a 2+ pattern must not be 4+ or
    /// 6+, a 3+ not 6+, a 1+ not 2+ to 6+ (the source stops at 6+ for run time
    /// and comments "up to 7+").
    ///
    /// Satellites record original peak indices, but the left-peak test
    /// compares the satellites against the **white** index of the left peak, as
    /// the source does (`MultiplexFiltering.cpp`, `checkSatellite(peak.getRTidx(),
    /// mz_idx)` with `mz_idx` from `checkForSignificantPeak_` on the white
    /// spectrum). Once peaks of a spectrum are blacklisted the two index spaces
    /// differ; the port reproduces the comparison.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] where the source reads past the
    /// pattern's shift list or past the experiment.
    pub fn filter_peak_positions(
        &self,
        mz: f64,
        band: Range<usize>,
        pattern: &MultiplexIsotopicPeakPattern,
        peak: &mut MultiplexFilteredPeak,
    ) -> Result<bool> {
        let out_of_range =
            || Error::InvalidValue("MultiplexFiltering: index outside the experiment".into());
        let primary_mark = *self
            .blacklist
            .get(peak.rt_idx())
            .and_then(|marks| marks.get(peak.mz_idx()))
            .ok_or_else(out_of_range)?;
        if primary_mark > 0 {
            return Ok(false);
        }
        let tolerance = self.absolute_tolerance(mz);
        let ipp_max = self.isotopes_per_peptide_max;
        let white = &self.exp_centroided_white.spectra;

        let mut length = 0usize;
        let mut interrupted = false;
        for isotope in 0..ipp_max {
            let mut found_in_all_peptides = true;
            for peptide in 0..pattern.mass_shift_count() {
                let mz_shift_idx = peptide * ipp_max + isotope;
                let shift = mz_shift(pattern, mz_shift_idx)?;
                let mut found = false;
                for rt_idx in band.clone() {
                    let spectrum = white.get(rt_idx).ok_or_else(out_of_range)?;
                    if let Some(i) =
                        find_nearest_with_tolerance(&spectrum.peaks, mz + shift, tolerance)
                    {
                        let mz_idx = *self
                            .exp_centroided_mapping
                            .get(rt_idx)
                            .and_then(|m| m.get(i))
                            .ok_or_else(out_of_range)?;
                        let mark = self.blacklist[rt_idx][mz_idx];
                        if (mark == -1 || i64::from(mark) == mz_shift_idx as i64)
                            && !peak.check_satellite(rt_idx, mz_idx)
                        {
                            found = true;
                            peak.add_satellite(rt_idx, mz_idx, mz_shift_idx);
                        }
                    }
                }
                if !found {
                    found_in_all_peptides = false;
                }
            }
            if found_in_all_peptides && !interrupted {
                length += 1;
            } else {
                interrupted = true;
                if length < self.isotopes_per_peptide_min {
                    return Ok(false);
                }
            }
        }

        let rt_idx = peak.rt_idx();
        let spectrum = white.get(rt_idx).ok_or_else(out_of_range)?;
        let peak_mz = peak.mz();
        for peptide in 0..pattern.mass_shift_count() {
            let shift_first = mz_shift(pattern, peptide * ipp_max)?;
            let mz_first_peak = peak_mz + shift_first;
            let Some(idx_first) =
                find_nearest_with_tolerance(&spectrum.peaks, mz_first_peak, tolerance)
            else {
                continue;
            };
            let intensity_first_peak = f64::from(spectrum.peaks[idx_first].intensity);
            let shift_second = mz_shift(pattern, peptide * ipp_max + 1)?;

            // Is there a significant zeroth peak?
            let left = peak_mz + 2.0 * shift_first - shift_second;
            if let Some(mz_idx) =
                self.check_for_significant_peak(left, 2.0 * tolerance, rt_idx, intensity_first_peak)
            {
                if !peak.check_satellite(rt_idx, mz_idx) {
                    return Ok(false);
                }
            }
            // Mistaken charge state identities: patterns are searched by most
            // common occurrence, not by decreasing charge.
            let mut probes = Vec::new();
            match pattern.charge() {
                2 => {
                    probes.push(peak_mz + shift_first / 2.0 + shift_second / 2.0);
                    probes.push(peak_mz + shift_first * 2.0 / 3.0 + shift_second / 3.0);
                }
                3 => probes.push(peak_mz + shift_first / 2.0 + shift_second / 2.0),
                1 => {
                    for c in 2..7 {
                        let c = f64::from(c);
                        probes.push(peak_mz + shift_first * (c - 1.0) / c + shift_second / c);
                    }
                }
                _ => {}
            }
            for probe in probes {
                if self
                    .check_for_significant_peak(
                        probe,
                        2.0 * tolerance,
                        rt_idx,
                        intensity_first_peak,
                    )
                    .is_some()
                {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    /// Blacklist the mass traces of `peak` found for pattern `pattern_idx`
    /// (source protected `blacklistPeak_(const MultiplexFilteredPeak&,
    /// unsigned)`).
    ///
    /// Each satellite belongs to a mass trace. For each trace, from the
    /// earliest to the latest spectrum holding one of its satellites, extended
    /// by `2 * rt_band` on both sides, every original peak nearest to the
    /// trace's m/z within the tolerance is marked with the trace index, whether
    /// or not it is a satellite. The later bound is `RTBegin(rt + 2 * rt_band)`
    /// plus one spectrum unless at the end, as in the source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a pattern or satellite index outside
    /// the patterns or the experiment, where the source reads out of bounds,
    /// and for a mass trace index beyond `i32`.
    pub fn blacklist_peak(
        &mut self,
        peak: &MultiplexFilteredPeak,
        pattern_idx: usize,
    ) -> Result<()> {
        let tolerance = self.absolute_tolerance(peak.mz());
        let pattern = self.patterns.get(pattern_idx).ok_or_else(|| {
            Error::InvalidValue("MultiplexFiltering: pattern index out of range".into())
        })?;
        let mut rt_boundaries: BTreeMap<usize, (usize, usize)> = BTreeMap::new();
        for (trace, satellite) in peak.satellite_entries() {
            let entry = rt_boundaries
                .entry(trace)
                .or_insert((satellite.rt_idx, satellite.rt_idx));
            entry.0 = entry.0.min(satellite.rt_idx);
            entry.1 = entry.1.max(satellite.rt_idx);
        }
        let spectra = &self.exp_centroided.spectra;
        let mut marks = Vec::new();
        for (&trace, &(first, last)) in &rt_boundaries {
            let mz = peak.mz() + mz_shift(pattern, trace)?;
            let out_of_range = || {
                Error::InvalidValue(
                    "MultiplexFiltering: satellite index outside the experiment".into(),
                )
            };
            let begin_rt = spectra.get(first).ok_or_else(out_of_range)?.rt;
            let end_rt = spectra.get(last).ok_or_else(out_of_range)?.rt;
            let begin = rt_begin(spectra, begin_rt - 2.0 * self.rt_band);
            let mut end = rt_begin(spectra, end_rt + 2.0 * self.rt_band);
            if end != spectra.len() {
                end += 1;
            }
            let mark = i32::try_from(trace).map_err(|_| {
                Error::InvalidValue("MultiplexFiltering: mass trace index beyond int".into())
            })?;
            for (rt_idx, spectrum) in spectra.iter().enumerate().take(end).skip(begin) {
                if let Some(idx_mz) = find_nearest_with_tolerance(&spectrum.peaks, mz, tolerance) {
                    marks.push((rt_idx, idx_mz, mark));
                }
            }
        }
        for (rt_idx, mz_idx, mark) in marks {
            self.blacklist[rt_idx][mz_idx] = mark;
        }
        Ok(())
    }

    /// The averagine distribution for a peptide of `mass`, retaining
    /// `isotopes_per_peptide_max` peaks (all when zero) in the source's binary32
    /// arithmetic.
    fn averagine(&self, mass: f64) -> Result<IsotopeDistribution> {
        let kind = AveragineType::from_name(&self.averagine_type)
            .ok_or_else(|| Error::InvalidValue("Invalid averagine type.".into()))?;
        let max = (self.isotopes_per_peptide_max > 0).then_some(self.isotopes_per_peptide_max);
        let solver = CoarseIsotopePatternGenerator::new(max, CoarseMassMode::default())?
            .with_precision(ProbabilityPrecision::SourceSingle);
        match kind {
            AveragineType::Peptide => solver.estimate_from_peptide_weight(mass),
            AveragineType::Rna => solver.estimate_from_rna_weight(mass),
            AveragineType::Dna => solver.estimate_from_dna_weight(mass),
        }
    }

    /// The similarity the averagine filter requires: `p` for multiplets and
    /// `p + x(1 - p)` for singlets.
    fn averagine_threshold(&self, pattern: &MultiplexIsotopicPeakPattern) -> f64 {
        if pattern.mass_shift_count() == 1 {
            self.averagine_similarity
                + self.averagine_similarity_scaling * (1.0 - self.averagine_similarity)
        } else {
            self.averagine_similarity
        }
    }

    /// The averagine test on per-isotope mean intensities: for each peptide,
    /// the model intensity and the mean of `intensities(idx)` of every isotope
    /// with data, both correlations at least the threshold.
    fn averagine_test(
        &self,
        pattern: &MultiplexIsotopicPeakPattern,
        mass: f64,
        mut intensities: impl FnMut(usize) -> Vec<f64>,
    ) -> Result<bool> {
        let distribution = self.averagine(mass)?;
        let similarity = self.averagine_threshold(pattern);
        let ipp_max = self.isotopes_per_peptide_max;
        for peptide in 0..pattern.mass_shift_count() {
            let mut intensities_model = Vec::new();
            let mut intensities_data = Vec::new();
            for isotope in 0..ipp_max {
                let values = intensities(peptide * ipp_max + isotope);
                if !values.is_empty() {
                    let model = distribution.peaks().get(isotope).ok_or_else(|| {
                        Error::InvalidValue(
                            "MultiplexFiltering: averagine distribution shorter than the isotope count \
                             (the source reads past its end)"
                                .into(),
                        )
                    })?;
                    let mut sum = 0.0;
                    for v in &values {
                        sum += v;
                    }
                    // `distribution[isotope].getIntensity()` is a binary32
                    // intensity widened to double.
                    intensities_model.push(f64::from(model.probability as f32));
                    intensities_data.push(sum / values.len() as f64);
                }
            }
            if intensities_model.len() < self.isotopes_per_peptide_min
                || intensities_data.len() < self.isotopes_per_peptide_min
            {
                return Err(Error::InvalidRange(
                    "intensity vectors too small for correlation calculation".into(),
                ));
            }
            let pearson = pearson_correlation_coefficient(&intensities_model, &intensities_data)?;
            let spearman = rank_correlation_coefficient(&intensities_model, &intensities_data)?;
            if pearson < similarity || spearman < similarity {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// Whether the satellite intensities of `peak` follow the averagine model
    /// (source protected `filterAveragineModel_`).
    ///
    /// The model is the averagine distribution of mass `mz * charge` (the
    /// lightest peptide's; the peptides of a multiplet are close in mass). For
    /// every peptide, each isotope with satellites contributes the model
    /// intensity and the mean satellite intensity, and both the Pearson and
    /// the Spearman rank correlation of the two must reach the averagine
    /// similarity, made stricter for singlets as
    /// [`MultiplexFiltering::new`] describes. A NaN correlation (a constant
    /// intensity vector) passes, as `NaN < p` is false in the source too.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an unknown averagine type ("Invalid
    /// averagine type.", source `Exception::InvalidParameter`) and when the
    /// distribution has fewer peaks than an isotope index needs (undefined in
    /// the source); [`Error::InvalidRange`] "intensity vectors too small for
    /// correlation calculation" when a peptide has fewer isotopes with data than
    /// `isotopes_per_peptide_min` (source `Exception::InvalidSize`); and the
    /// errors of the correlation functions and the isotope generator.
    pub fn filter_averagine_model(
        &self,
        pattern: &MultiplexIsotopicPeakPattern,
        peak: &MultiplexFilteredPeak,
    ) -> Result<bool> {
        let mass = peak.mz() * f64::from(pattern.charge());
        let spectra = &self.exp_centroided.spectra;
        let mut failure = None;
        let result = self.averagine_test(pattern, mass, |idx| {
            let mut values = Vec::new();
            for satellite in peak.satellites().get(&idx).into_iter().flatten() {
                match spectra
                    .get(satellite.rt_idx)
                    .and_then(|s| s.peaks.get(satellite.mz_idx))
                {
                    Some(p) => values.push(f64::from(p.intensity)),
                    None => failure = Some(()),
                }
            }
            values
        });
        if failure.is_some() {
            return Err(Error::InvalidValue(
                "MultiplexFiltering: satellite index outside the experiment".into(),
            ));
        }
        result
    }

    /// The averagine filter on spline-interpolated satellite data points (the
    /// source's private `MultiplexFilteringProfile::filterAveragineModel_`),
    /// the same test as [`Self::filter_averagine_model`] on the data points'
    /// intensities.
    ///
    /// # Errors
    ///
    /// As [`Self::filter_averagine_model`].
    pub fn filter_averagine_model_profile(
        &self,
        pattern: &MultiplexIsotopicPeakPattern,
        peak: &MultiplexFilteredPeak,
        satellites_profile: &SatelliteMultimap<MultiplexSatelliteProfile>,
    ) -> Result<bool> {
        let mass = peak.mz() * f64::from(pattern.charge());
        self.averagine_test(pattern, mass, |idx| {
            satellites_profile
                .get(&idx)
                .into_iter()
                .flatten()
                .map(|s| f64::from(s.intensity))
                .collect()
        })
    }

    /// The correlation test between every pair of peptides: for each isotope,
    /// every pair of satellites of the two peptides' mass traces that
    /// `same_spectrum` pairs contributes their intensities.
    fn peptide_correlation_test<T>(
        &self,
        pattern: &MultiplexIsotopicPeakPattern,
        satellites: &SatelliteMultimap<T>,
        mut pair: impl FnMut(&T, &T) -> Result<Option<(f64, f64)>>,
    ) -> Result<bool> {
        if pattern.mass_shift_count() < 2 {
            return Ok(true);
        }
        let ipp_max = self.isotopes_per_peptide_max;
        let count = pattern.mass_shift_count();
        for peptide_1 in 0..count - 1 {
            for peptide_2 in peptide_1 + 1..count {
                let mut intensities_1 = Vec::new();
                let mut intensities_2 = Vec::new();
                for isotope in 0..ipp_max {
                    let idx_1 = peptide_1 * ipp_max + isotope;
                    let idx_2 = peptide_2 * ipp_max + isotope;
                    for s1 in satellites.get(&idx_1).into_iter().flatten() {
                        for s2 in satellites.get(&idx_2).into_iter().flatten() {
                            if let Some((i1, i2)) = pair(s1, s2)? {
                                intensities_1.push(i1);
                                intensities_2.push(i2);
                            }
                        }
                    }
                }
                if intensities_1.is_empty() || intensities_2.is_empty() {
                    return Ok(false);
                }
                let pearson = pearson_correlation_coefficient(&intensities_1, &intensities_2)?;
                let spearman = rank_correlation_coefficient(&intensities_1, &intensities_2)?;
                if pearson < self.peptide_similarity || spearman < self.peptide_similarity {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    /// Whether corresponding satellites of the different peptides of `peak`
    /// correlate (source protected `filterPeptideCorrelation_`).
    ///
    /// Peptides of one multiplet share their sequence and so their isotope
    /// distribution. For every pair of peptides, the intensities of satellites
    /// of the same isotope in the same spectrum are paired, and both the
    /// Pearson and the Spearman rank correlation must reach the peptide
    /// similarity. No pair at all fails the filter. Singlets pass.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a satellite outside the experiment,
    /// and the errors of the correlation functions.
    pub fn filter_peptide_correlation(
        &self,
        pattern: &MultiplexIsotopicPeakPattern,
        peak: &MultiplexFilteredPeak,
    ) -> Result<bool> {
        let spectra = &self.exp_centroided.spectra;
        let intensity = |rt_idx: usize, mz_idx: usize| -> Result<f64> {
            spectra
                .get(rt_idx)
                .and_then(|s| s.peaks.get(mz_idx))
                .map(|p| f64::from(p.intensity))
                .ok_or_else(|| {
                    Error::InvalidValue(
                        "MultiplexFiltering: satellite index outside the experiment".into(),
                    )
                })
        };
        self.peptide_correlation_test(pattern, peak.satellites(), |s1, s2| {
            if s1.rt_idx == s2.rt_idx {
                Ok(Some((
                    intensity(s1.rt_idx, s1.mz_idx)?,
                    intensity(s2.rt_idx, s2.mz_idx)?,
                )))
            } else {
                Ok(None)
            }
        })
    }

    /// The peptide correlation filter on spline-interpolated satellite data
    /// points (the source's private
    /// `MultiplexFilteringProfile::filterPeptideCorrelation_`): data points
    /// pair when their `f32` retention times are equal.
    ///
    /// # Errors
    ///
    /// The errors of the correlation functions.
    pub fn filter_peptide_correlation_profile(
        &self,
        pattern: &MultiplexIsotopicPeakPattern,
        satellites_profile: &SatelliteMultimap<MultiplexSatelliteProfile>,
    ) -> Result<bool> {
        self.peptide_correlation_test(pattern, satellites_profile, |s1, s2| {
            // The source widens both `float` RTs to double before `==`.
            if f64::from(s1.rt) == f64::from(s2.rt) {
                Ok(Some((f64::from(s1.intensity), f64::from(s2.intensity))))
            } else {
                Ok(None)
            }
        })
    }
}
