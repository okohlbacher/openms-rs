// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Intensity, trace and isotope-pattern scores of the picked feature finder
//! (`FEATUREFINDER/FeatureFinderAlgorithmPicked.h`).
//!
//! Each peak of the input receives three scores between 0 and 1 before seeds are
//! selected, as in the source's steps 1, 2 and 3.1:
//!
//! - the **intensity score** says how significant the peak's intensity is in its
//!   local environment, interpolated from the 20-quantiles of intensity bins
//!   ([`IntensityThresholds`](crate::analysis::feature_finder_picked::scoring::IntensityThresholds));
//! - the **trace score** says how well the peak's m/z recurs in the neighbouring
//!   scans, and a flag records whether it is the local maximum of that trace;
//! - the **pattern score**, one per charge, says how well an averagine isotope
//!   pattern containing the peak fits the data
//!   ([`find_isotope`](crate::analysis::feature_finder_picked::scoring::find_isotope),
//!   [`isotope_score`](crate::analysis::feature_finder_picked::scoring::isotope_score)).
//!
//! The source stores them as `float` data arrays of each spectrum, named
//! `trace_score`, `intensity_score`, `local_max`, `pattern_score_<charge>` and
//! `overall_score_<charge>`. They are only ever read by the algorithm itself,
//! and written out only by the debug mode
//! ([`debug_experiment`](crate::analysis::feature_finder_picked::debug::debug_experiment)
//! rebuilds those arrays for it). They live in
//! [`ScoreArrays`](crate::analysis::feature_finder_picked::scoring::ScoreArrays)
//! instead, one flat `f32` array per score, so the input spectra
//! keep their own data arrays and no per-spectrum allocation is needed.
//!
//! Arithmetic follows the source operation by operation, including the `float`
//! narrowing of every stored score. See `docs/FEATURE_FINDER_PICKED_SUPPORT.md`.

use crate::analysis::feature_finder_picked::debug::{LogSink, NoLog, g, number, put_all};
use crate::analysis::feature_finder_picked::helper_structs::{
    IsotopePattern, PatternPeak, TheoreticalIsotopePattern,
};
use crate::kernel::{MSExperiment, MSSpectrum, NumericRange, Peak1D, nearest};
use crate::math::statistic_functions::pearson_correlation_coefficient;
use crate::{Error, Result};

/// Number of stored quantiles per intensity bin: the 0th to the 20th
/// 20-quantile.
pub const QUANTILE_COUNT: usize = 21;

/// A score between 0 and 1 for the m/z deviation of two peaks: source
/// `positionScore_`.
///
/// With `d = |pos1 - pos2|` and `a = allowed_deviation`: `0.1 * (0.5a - d) /
/// (0.5a) + 0.9` when `d <= 0.5a`, `0.9 * (a - d) / (0.5a)` when `d <= a`, and 0
/// otherwise. The operations are evaluated in the source's order. A deviation of
/// zero with a tolerance of zero divides zero by zero and gives NaN, as in the
/// source.
pub fn position_score(pos1: f64, pos2: f64, allowed_deviation: f64) -> f64 {
    let diff = (pos1 - pos2).abs();
    if diff <= 0.5 * allowed_deviation {
        0.1 * (0.5 * allowed_deviation - diff) / (0.5 * allowed_deviation) + 0.9
    } else if diff <= allowed_deviation {
        0.9 * (allowed_deviation - diff) / (0.5 * allowed_deviation)
    } else {
        0.0
    }
}

/// The index of the peak nearest to `pos`, searching linearly upwards from
/// `start`: source `nearest_`.
///
/// Moves to the next peak while it is strictly closer to `pos` and returns the
/// last index moved to. The walk therefore stops at the first local minimum of
/// the distance, which is the nearest peak when `start` lies at or below it in a
/// spectrum sorted by m/z. Ties keep the lower index. The second value is the
/// number of steps taken.
///
/// Returns `None` when `start` is not a peak index; the source reads past the
/// end of the spectrum there.
pub fn nearest_from(peaks: &[Peak1D], pos: f64, start: usize) -> Option<(usize, usize)> {
    let mut distance = (pos - peaks.get(start)?.mz).abs();
    let mut index = start + 1;
    while let Some(peak) = peaks.get(index) {
        let new_distance = (pos - peak.mz).abs();
        if new_distance < distance {
            distance = new_distance;
            index += 1;
        } else {
            break;
        }
    }
    Some((index - 1, index - 1 - start))
}

/// Precalculated intensity 20-quantiles of a regular RT by m/z grid (source
/// members `intensity_rt_step_`, `intensity_mz_step_` and
/// `intensity_thresholds_`, filled in step 1 of `run_`).
#[derive(Clone, Debug, PartialEq)]
pub struct IntensityThresholds {
    bins: usize,
    rt_start: f64,
    mz_start: f64,
    rt_step: f64,
    mz_step: f64,
    quantiles: Vec<[f64; QUANTILE_COUNT]>,
}

impl IntensityThresholds {
    /// Bin the intensities of `experiment` into `bins` by `bins` cells and store
    /// 21 quantiles per cell: step 1 of source `run_`.
    ///
    /// The grid spans the MS1 retention-time and m/z ranges: bin `i` of a
    /// dimension covers `[start + i * step, start + (i + 1) * step]` with
    /// `step = (max - start) / bins`, and both borders are inclusive, as for
    /// `MSExperiment::areaBeginConst`, so a peak on a border belongs to both
    /// cells. Each cell's intensities are promoted to `f64` and sorted, and
    /// quantile `i` is element `floor(0.05 * i * (n - 1))`. An empty cell keeps
    /// 21 zeros.
    ///
    /// The source walks each cell with an area iterator; this walks the same
    /// scans and peak ranges with two binary searches per cell and scan, so the
    /// whole experiment is not revalidated once per cell.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `bins` is zero, when a spectrum is
    /// not MS1 or holds a non-finite value, when there is no MS1 peak, or when
    /// the retention-time or m/z range has zero width, and
    /// [`Error::UnsortedData`] when the spectra are not sorted by RT and m/z. A
    /// zero-width range makes the source divide by a zero step and convert the
    /// resulting NaN to an unsigned integer, which is undefined behaviour.
    pub fn compute(experiment: &MSExperiment, bins: usize) -> Result<Self> {
        let mut work = Work::unlimited();
        Self::compute_with_work(experiment, bins, &mut work)
    }

    pub(crate) fn compute_with_work(
        experiment: &MSExperiment,
        bins: usize,
        work: &mut Work,
    ) -> Result<Self> {
        if bins == 0 {
            return Err(Error::InvalidValue(
                "intensity:bins must be positive".into(),
            ));
        }
        let (rt, mz) = ms1_ranges(experiment)?;
        let bins_f = bins as f64;
        let rt_step = (rt.max - rt.min) / bins_f;
        let mz_step = (mz.max - mz.min) / bins_f;
        if rt_step == 0.0 || mz_step == 0.0 {
            return Err(Error::InvalidValue(format!(
                "FeatureFinderAlgorithmPicked needs a retention-time and an m/z range of \
                 positive width (RT {} to {}, m/z {} to {}); the source divides by a zero bin \
                 width here",
                rt.min, rt.max, mz.min, mz.max
            )));
        }
        let cells = bins
            .checked_mul(bins)
            .ok_or_else(|| Error::InvalidValue("intensity bin count overflow".into()))?;
        let mut quantiles = Vec::new();
        quantiles
            .try_reserve_exact(cells)
            .map_err(|_| Error::InvalidValue("cannot allocate the intensity bins".into()))?;
        let spectra = &experiment.spectra;
        let mut values: Vec<f64> = Vec::new();
        for rt_bin in 0..bins {
            let min_rt = rt.min + rt_bin as f64 * rt_step;
            let max_rt = rt.min + (rt_bin + 1) as f64 * rt_step;
            let begin = spectra.partition_point(|spectrum| spectrum.rt < min_rt);
            let end = spectra.partition_point(|spectrum| spectrum.rt <= max_rt);
            for mz_bin in 0..bins {
                let min_mz = mz.min + mz_bin as f64 * mz_step;
                let max_mz = mz.min + (mz_bin + 1) as f64 * mz_step;
                values.clear();
                work.consume(end.saturating_sub(begin) as u64 + 1)?;
                for spectrum in spectra.get(begin..end).unwrap_or(&[]) {
                    let peaks = &spectrum.peaks;
                    let low = peaks.partition_point(|peak| peak.mz < min_mz);
                    let high = peaks.partition_point(|peak| peak.mz <= max_mz);
                    if let Some(selected) = peaks.get(low..high) {
                        values.extend(selected.iter().map(|peak| f64::from(peak.intensity)));
                    }
                }
                work.consume(values.len() as u64)?;
                let mut cell = [0.0; QUANTILE_COUNT];
                if !values.is_empty() {
                    values.sort_unstable_by(f64::total_cmp);
                    let last = (values.len() - 1) as f64;
                    for (i, quantile) in cell.iter_mut().enumerate() {
                        let index = (0.05 * i as f64 * last).floor() as usize;
                        *quantile = values[index.min(values.len() - 1)];
                    }
                }
                quantiles.push(cell);
            }
        }
        Ok(Self {
            bins,
            rt_start: rt.min,
            mz_start: mz.min,
            rt_step,
            mz_step,
            quantiles,
        })
    }

    /// Bins per dimension (source `intensity_bins_`).
    pub fn bins(&self) -> usize {
        self.bins
    }

    /// The smallest MS1 retention time, where the first RT bin starts.
    pub fn rt_start(&self) -> f64 {
        self.rt_start
    }

    /// The smallest MS1 m/z, where the first m/z bin starts.
    pub fn mz_start(&self) -> f64 {
        self.mz_start
    }

    /// RT bin width (source `intensity_rt_step_`).
    pub fn rt_step(&self) -> f64 {
        self.rt_step
    }

    /// m/z bin width (source `intensity_mz_step_`).
    pub fn mz_step(&self) -> f64 {
        self.mz_step
    }

    /// The 21 ascending quantiles of a cell, or `None` outside the grid.
    pub fn quantiles(&self, rt_bin: usize, mz_bin: usize) -> Option<&[f64; QUANTILE_COUNT]> {
        if rt_bin >= self.bins || mz_bin >= self.bins {
            return None;
        }
        self.quantiles.get(rt_bin * self.bins + mz_bin)
    }

    /// The intensity score of `intensity` in one cell: source
    /// `intensityScore_(rt_bin, mz_bin, intensity)`.
    ///
    /// Finds the first quantile `q[k]` not below `intensity`. Above the largest
    /// quantile the score is 1. At `k = 0` the bin score is
    /// `0.05 * intensity / q[0]`, otherwise `0.05 * (intensity - q[k-1]) / (q[k] -
    /// q[k-1])`; the result is that bin score plus `0.05 * (k - 1)`, clamped to
    /// `[0, 1]`. A NaN from `0 / 0` (a non-positive intensity against a zero first
    /// quantile) passes the clamp unchanged, as in the source.
    ///
    /// Returns `None` for a cell outside the grid.
    pub fn bin_score(&self, rt_bin: usize, mz_bin: usize, intensity: f64) -> Option<f64> {
        let quantiles = self.quantiles(rt_bin, mz_bin)?;
        let position = quantiles.partition_point(|&quantile| quantile < intensity);
        let Some(&upper) = quantiles.get(position) else {
            return Some(1.0);
        };
        let bin_score = if position == 0 {
            0.05 * intensity / upper
        } else {
            let lower = quantiles[position - 1];
            0.05 * (intensity - lower) / (upper - lower)
        };
        // `clamp` keeps NaN and -0.0, as the source's two comparisons do.
        Some((bin_score + 0.05 * (position as f64 - 1.0)).clamp(0.0, 1.0))
    }

    /// The intensity score of a peak: source `intensityScore_(spectrum, peak)`.
    ///
    /// The peak's position on a half-bin grid, `floor((x - start) / step * 2)`
    /// capped at `2 * bins - 1`, selects the two nearest bins per dimension (one
    /// at the outer half-bins). The four cell scores of [`Self::bin_score`] are
    /// weighted by `d = sqrt((1 - d_rt)^2 + (1 - d_mz)^2)`, where `d_rt` and
    /// `d_mz` are the distances of the peak to each bin centre in bin widths,
    /// and each weight is divided by the sum of the four. The squares are
    /// products where the source calls `std::pow(x, 2)`; the executed intensity
    /// scores with 7 and 10 bins agree with them bit for bit.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a half-bin position is negative or
    /// not finite, which happens only for a position outside the binned range;
    /// the source converts it to an unsigned integer, which is undefined.
    pub fn score(&self, rt: f64, mz: f64, intensity: f64) -> Result<f64> {
        let rt_bin = self.half_bin(rt, self.rt_start, self.rt_step)?;
        let mz_bin = self.half_bin(mz, self.mz_start, self.mz_step)?;
        let last = 2 * self.bins - 1;
        let neighbours = |bin: usize| -> (usize, usize) {
            if bin == 0 || bin == last {
                (bin / 2, bin / 2)
            } else if bin % 2 == 1 {
                (bin / 2, bin / 2 + 1)
            } else {
                (bin / 2 - 1, bin / 2)
            }
        };
        let (ml, mh) = neighbours(mz_bin);
        let (rl, rh) = neighbours(rt_bin);
        let drl = (self.rt_start + (0.5 + rl as f64) * self.rt_step - rt).abs() / self.rt_step;
        let drh = (self.rt_start + (0.5 + rh as f64) * self.rt_step - rt).abs() / self.rt_step;
        let dml = (self.mz_start + (0.5 + ml as f64) * self.mz_step - mz).abs() / self.mz_step;
        let dmh = (self.mz_start + (0.5 + mh as f64) * self.mz_step - mz).abs() / self.mz_step;
        let square = |x: f64| x * x;
        let d1 = (square(1.0 - drl) + square(1.0 - dml)).sqrt();
        let d2 = (square(1.0 - drh) + square(1.0 - dml)).sqrt();
        let d3 = (square(1.0 - drl) + square(1.0 - dmh)).sqrt();
        let d4 = (square(1.0 - drh) + square(1.0 - dmh)).sqrt();
        let d_sum = d1 + d2 + d3 + d4;
        let cell = |r: usize, m: usize| {
            self.bin_score(r, m, intensity)
                .ok_or_else(|| Error::InvalidValue("intensity bin outside the grid".into()))
        };
        Ok(cell(rl, ml)? * (d1 / d_sum)
            + cell(rh, ml)? * (d2 / d_sum)
            + cell(rl, mh)? * (d3 / d_sum)
            + cell(rh, mh)? * (d4 / d_sum))
    }

    fn half_bin(&self, x: f64, start: f64, step: f64) -> Result<usize> {
        let position = ((x - start) / step * 2.0).floor();
        if position.is_nan() || position.is_infinite() || position < 0.0 {
            return Err(Error::InvalidValue(format!(
                "position {x} lies outside the intensity bins starting at {start}"
            )));
        }
        let last = 2 * self.bins - 1;
        Ok(if position >= last as f64 {
            last
        } else {
            position as usize
        })
    }
}

/// The per-peak scores of one run: the source's float data arrays.
///
/// Scores of spectrum `s` are the slices returned for `s`, aligned with its
/// peaks. Charges are addressed by their index `charge - charge_low`.
/// Trace scores, local-maximum flags and overall scores stay zero for the first
/// and last `min_spectra` scans, as in the source.
#[derive(Clone, Debug, PartialEq)]
pub struct ScoreArrays {
    offsets: Vec<usize>,
    charge_low: i32,
    trace: Vec<f32>,
    intensity: Vec<f32>,
    local_max: Vec<f32>,
    pattern: Vec<Vec<f32>>,
    overall: Vec<Vec<f32>>,
}

impl ScoreArrays {
    /// Zero-initialised arrays for the peaks of `experiment` and `charge_count`
    /// charges starting at `charge_low`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the arrays would exceed
    /// `max_bytes` or cannot be allocated.
    pub fn new(
        experiment: &MSExperiment,
        charge_low: i32,
        charge_count: usize,
        max_bytes: usize,
    ) -> Result<Self> {
        let mut offsets = Vec::with_capacity(experiment.spectra.len() + 1);
        let mut total = 0usize;
        offsets.push(0);
        for spectrum in &experiment.spectra {
            total = total
                .checked_add(spectrum.peaks.len())
                .ok_or_else(|| Error::InvalidValue("peak count overflow".into()))?;
            offsets.push(total);
        }
        let arrays = charge_count
            .checked_mul(2)
            .and_then(|n| n.checked_add(3))
            .ok_or_else(|| Error::InvalidValue("score array count overflow".into()))?;
        let bytes = total
            .checked_mul(arrays)
            .and_then(|n| n.checked_mul(std::mem::size_of::<f32>()))
            .ok_or_else(|| Error::InvalidValue("score array size overflow".into()))?;
        if bytes > max_bytes {
            return Err(Error::InvalidValue(format!(
                "score arrays of {bytes} bytes exceed the limit of {max_bytes}"
            )));
        }
        let zeros = |length: usize| -> Result<Vec<f32>> {
            let mut array = Vec::new();
            array
                .try_reserve_exact(length)
                .map_err(|_| Error::InvalidValue("cannot allocate the score arrays".into()))?;
            array.resize(length, 0.0);
            Ok(array)
        };
        let mut pattern = Vec::with_capacity(charge_count);
        let mut overall = Vec::with_capacity(charge_count);
        for _ in 0..charge_count {
            pattern.push(zeros(total)?);
            overall.push(zeros(total)?);
        }
        Ok(Self {
            offsets,
            charge_low,
            trace: zeros(total)?,
            intensity: zeros(total)?,
            local_max: zeros(total)?,
            pattern,
            overall,
        })
    }

    /// Number of spectra.
    pub fn spectrum_count(&self) -> usize {
        self.offsets.len() - 1
    }

    /// Number of charges.
    pub fn charge_count(&self) -> usize {
        self.pattern.len()
    }

    /// The lowest charge; charge index 0.
    pub fn charge_low(&self) -> i32 {
        self.charge_low
    }

    /// The source names of the arrays, in the source's order: `trace_score`,
    /// `intensity_score`, `local_max`, then `pattern_score_<c>` and
    /// `overall_score_<c>` for every charge.
    pub fn array_names(&self) -> Vec<String> {
        let charges = (0..self.charge_count()).map(|i| i64::from(self.charge_low) + i as i64);
        ["trace_score", "intensity_score", "local_max"]
            .iter()
            .map(|name| (*name).to_string())
            .chain(charges.clone().map(|c| format!("pattern_score_{c}")))
            .chain(charges.map(|c| format!("overall_score_{c}")))
            .collect()
    }

    fn range(&self, spectrum: usize) -> Option<std::ops::Range<usize>> {
        Some(*self.offsets.get(spectrum)?..*self.offsets.get(spectrum + 1)?)
    }

    /// Trace scores of a spectrum's peaks (source array 0, `trace_score`).
    pub fn trace(&self, spectrum: usize) -> Option<&[f32]> {
        self.trace.get(self.range(spectrum)?)
    }

    /// Intensity scores of a spectrum's peaks (source array 1, `intensity_score`).
    pub fn intensity(&self, spectrum: usize) -> Option<&[f32]> {
        self.intensity.get(self.range(spectrum)?)
    }

    /// Local-maximum flags, 1 or 0, of a spectrum's peaks (source array 2,
    /// `local_max`).
    pub fn local_max(&self, spectrum: usize) -> Option<&[f32]> {
        self.local_max.get(self.range(spectrum)?)
    }

    /// Pattern scores for charge index `charge` (source `pattern_score_<c>`).
    pub fn pattern(&self, charge: usize, spectrum: usize) -> Option<&[f32]> {
        self.pattern.get(charge)?.get(self.range(spectrum)?)
    }

    /// Overall scores for charge index `charge` (source `overall_score_<c>`).
    pub fn overall(&self, charge: usize, spectrum: usize) -> Option<&[f32]> {
        self.overall.get(charge)?.get(self.range(spectrum)?)
    }

    pub(crate) fn offset(&self, spectrum: usize) -> usize {
        self.offsets[spectrum]
    }

    pub(crate) fn intensity_mut(&mut self) -> &mut [f32] {
        &mut self.intensity
    }

    pub(crate) fn trace_and_local_max_mut(&mut self) -> (&mut [f32], &mut [f32]) {
        (&mut self.trace, &mut self.local_max)
    }

    pub(crate) fn pattern_mut(&mut self, charge: usize) -> &mut [f32] {
        &mut self.pattern[charge]
    }

    /// The flat trace, intensity, local-maximum and pattern arrays of one
    /// charge, and its overall array mutably.
    #[allow(clippy::type_complexity)]
    pub(crate) fn seed_inputs(
        &mut self,
        charge: usize,
    ) -> (&[f32], &[f32], &[f32], &[f32], &mut [f32]) {
        (
            &self.trace,
            &self.intensity,
            &self.local_max,
            &self.pattern[charge],
            &mut self.overall[charge],
        )
    }
}

/// A running work budget; [`Error::InvalidValue`] once it is exhausted.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Work {
    remaining: u64,
}

impl Work {
    pub(crate) fn new(limit: u64) -> Self {
        Self { remaining: limit }
    }

    pub(crate) fn unlimited() -> Self {
        Self::new(u64::MAX)
    }

    pub(crate) fn consume(&mut self, units: u64) -> Result<()> {
        self.remaining = self.remaining.checked_sub(units).ok_or_else(|| {
            Error::InvalidValue(
                "FeatureFinderAlgorithmPicked scoring exceeded its work limit".into(),
            )
        })?;
        Ok(())
    }
}

/// The MS1 retention-time and m/z ranges of a validated experiment: source
/// `spectrumRanges().byMSLevel(1)`, computed on demand.
///
/// Checks that every spectrum is MS1, finite and sorted, then takes the RT of
/// the first and last spectrum (empty spectra included, as the source extends
/// the RT range with every spectrum) and the extreme peak m/z values.
///
/// # Errors
///
/// As [`IntensityThresholds::compute`], without the zero-width checks.
pub(crate) fn ms1_ranges(experiment: &MSExperiment) -> Result<(NumericRange, NumericRange)> {
    check_ms1_sorted(experiment)?;
    let (Some(first), Some(last)) = (experiment.spectra.first(), experiment.spectra.last()) else {
        return Err(no_peak());
    };
    let rt = NumericRange {
        min: first.rt,
        max: last.rt,
    };
    let mut mz: Option<NumericRange> = None;
    for spectrum in &experiment.spectra {
        let (Some(low), Some(high)) = (spectrum.peaks.first(), spectrum.peaks.last()) else {
            continue;
        };
        mz = Some(match mz {
            None => NumericRange {
                min: low.mz,
                max: high.mz,
            },
            Some(range) => NumericRange {
                min: if low.mz < range.min {
                    low.mz
                } else {
                    range.min
                },
                max: if high.mz > range.max {
                    high.mz
                } else {
                    range.max
                },
            },
        });
    }
    Ok((rt, mz.ok_or_else(no_peak)?))
}

fn no_peak() -> Error {
    Error::InvalidValue(
        "FeatureFinderAlgorithmPicked needs at least one MS1 peak; the source reads an empty \
         m/z range here"
            .into(),
    )
}

fn check_ms1_sorted(experiment: &MSExperiment) -> Result<()> {
    let mut previous_rt = f64::NEG_INFINITY;
    for spectrum in &experiment.spectra {
        if spectrum.ms_level != 1 {
            return Err(Error::InvalidValue(
                "FeatureFinderAlgorithmPicked scores MS1 spectra only".into(),
            ));
        }
        if !spectrum.rt.is_finite() {
            return Err(Error::InvalidValue("non-finite retention time".into()));
        }
        if spectrum.rt < previous_rt {
            return Err(Error::UnsortedData);
        }
        previous_rt = spectrum.rt;
        let mut previous_mz = f64::NEG_INFINITY;
        for peak in &spectrum.peaks {
            if !peak.mz.is_finite() || !peak.intensity.is_finite() {
                return Err(Error::InvalidValue(
                    "non-finite peak m/z or intensity".into(),
                ));
            }
            if peak.mz < previous_mz {
                return Err(Error::UnsortedData);
            }
            previous_mz = peak.mz;
        }
    }
    Ok(())
}

/// Store the intensity score of every peak: the second half of step 1.
pub(crate) fn fill_intensity_scores(
    experiment: &MSExperiment,
    thresholds: &IntensityThresholds,
    scores: &mut ScoreArrays,
    work: &mut Work,
) -> Result<()> {
    let target = scores.intensity_mut();
    let mut index = 0;
    for spectrum in &experiment.spectra {
        work.consume(spectrum.peaks.len() as u64)?;
        for peak in &spectrum.peaks {
            target[index] =
                thresholds.score(spectrum.rt, peak.mz, f64::from(peak.intensity))? as f32;
            index += 1;
        }
    }
    Ok(())
}

/// Trace scores and local-maximum flags: step 2 of source `run_`.
///
/// For each peak of the scans `min_spectra..len - min_spectra`, the nearest
/// peak of each of the `min_spectra` following and then preceding non-empty
/// scans contributes its [`position_score`] against `trace_tolerance`; the sum
/// is divided by `2 * min_spectra`. The peak is a local maximum unless a
/// contributing neighbour with a positive position score is strictly more
/// intense (compared as `f32`).
pub(crate) fn fill_trace_scores(
    experiment: &MSExperiment,
    min_spectra: usize,
    trace_tolerance: f64,
    scores: &mut ScoreArrays,
    work: &mut Work,
) -> Result<()> {
    let spectra = &experiment.spectra;
    let end = spectra.len() - min_spectra.min(spectra.len());
    let divisor = (2 * min_spectra) as f64;
    let offsets: Vec<usize> = (0..spectra.len()).map(|s| scores.offset(s)).collect();
    let (trace, local_max) = scores.trace_and_local_max_mut();
    for s in min_spectra..end {
        let spectrum = &spectra[s];
        work.consume(
            (spectrum.peaks.len() as u64).saturating_mul((min_spectra as u64).saturating_mul(2)),
        )?;
        for (p, peak) in spectrum.peaks.iter().enumerate() {
            let pos = peak.mz;
            let intensity = peak.intensity;
            let mut trace_score = 0.0;
            let mut is_max_peak = true;
            let neighbours = (1..=min_spectra)
                .map(|i| &spectra[s + i])
                .chain((1..=min_spectra).map(|i| &spectra[s - i]));
            for next in neighbours {
                let Some(index) = nearest(&next.peaks, pos, |peak| peak.mz) else {
                    continue;
                };
                let found = next.peaks[index];
                let score = position_score(pos, found.mz, trace_tolerance);
                if score > 0.0 && found.intensity > intensity {
                    is_max_peak = false;
                }
                trace_score += score;
            }
            trace_score /= divisor;
            let slot = offsets[s] + p;
            trace[slot] = trace_score as f32;
            local_max[slot] = if is_max_peak { 1.0 } else { 0.0 };
        }
    }
    Ok(())
}

/// Search one isotope peak in a scan and its two neighbours: source
/// `findIsotope_`.
///
/// In scan `spectrum_index` the peak nearest to `pos` is found by
/// [`nearest_from`] starting at `peak_index`, which is updated to it. That peak
/// and the nearest peaks of the preceding and following non-empty scans each
/// match when their [`position_score`] against `pattern_tolerance` is not zero
/// (a NaN score matches). The pattern stores `pos` as the theoretical m/z, the
/// matched peak of the central scan or else of the first neighbour that
/// matched, the mean intensity and the mean position score of the matches, or
/// [`PatternPeak::NotFound`] and zeros when nothing matched. The intensities are
/// summed in `f64` after promoting each `f32`.
///
/// Returns the number of work units spent: the linear search steps plus the
/// binary searches.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `spectrum_index`, `pattern_index` or
/// `peak_index` is out of range; the source does not check them.
pub fn find_isotope(
    spectra: &[MSSpectrum],
    pos: f64,
    spectrum_index: usize,
    pattern: &mut IsotopePattern,
    pattern_index: usize,
    peak_index: &mut usize,
    pattern_tolerance: f64,
) -> Result<u64> {
    find_isotope_logged(
        spectra,
        pos,
        spectrum_index,
        pattern,
        pattern_index,
        peak_index,
        pattern_tolerance,
        &mut NoLog,
    )
}

/// [`find_isotope`] writing the source's debug lines to `log`: `   - Isotope
/// <i>: `, each match's intensity with one decimal (suffixed `b` and `a` for
/// the preceding and following scan), and ` missing` or `=> <mean>`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn find_isotope_logged<L: LogSink>(
    spectra: &[MSSpectrum],
    pos: f64,
    spectrum_index: usize,
    pattern: &mut IsotopePattern,
    pattern_index: usize,
    peak_index: &mut usize,
    pattern_tolerance: f64,
    log: &mut L,
) -> Result<u64> {
    let out_of_range = || Error::InvalidValue("findIsotope_ index out of range".into());
    let spectrum = spectra.get(spectrum_index).ok_or_else(out_of_range)?;
    if pattern_index >= pattern.peak.len()
        || pattern_index >= pattern.spectrum.len()
        || pattern_index >= pattern.intensity.len()
        || pattern_index >= pattern.mz_score.len()
        || pattern_index >= pattern.theoretical_mz.len()
    {
        return Err(out_of_range());
    }
    put_all(log, &["   - Isotope ", &pattern_index.to_string(), ": "]);
    let (nearest_index, steps) =
        nearest_from(&spectrum.peaks, pos, *peak_index).ok_or_else(out_of_range)?;
    *peak_index = nearest_index;
    let mut work = steps as u64 + 1;
    let mut intensity = 0.0;
    let mut pos_score = 0.0;
    let mut matches: u32 = 0;
    let this_mz_score = position_score(pos, spectrum.peaks[nearest_index].mz, pattern_tolerance);
    pattern.theoretical_mz[pattern_index] = pos;
    if this_mz_score != 0.0 {
        if log.enabled() {
            let intensity = f64::from(spectrum.peaks[nearest_index].intensity);
            put_all(log, &[&number(intensity, 1), " "]);
        }
        pattern.peak[pattern_index] = PatternPeak::Found(nearest_index);
        pattern.spectrum[pattern_index] = spectrum_index;
        intensity += f64::from(spectrum.peaks[nearest_index].intensity);
        pos_score += this_mz_score;
        matches += 1;
    }
    let neighbours = [
        (spectrum_index.checked_sub(1), "b "),
        (
            spectrum_index
                .checked_add(1)
                .filter(|&index| index < spectra.len()),
            "a ",
        ),
    ];
    for (neighbour_index, suffix) in neighbours {
        let Some(neighbour_index) = neighbour_index else {
            continue;
        };
        let neighbour = &spectra[neighbour_index];
        let Some(index) = nearest(&neighbour.peaks, pos, |peak| peak.mz) else {
            continue;
        };
        work += 1;
        let mz_score = position_score(pos, neighbour.peaks[index].mz, pattern_tolerance);
        if mz_score != 0.0 {
            if log.enabled() {
                let found = f64::from(neighbour.peaks[index].intensity);
                put_all(log, &[&number(found, 1), suffix]);
            }
            intensity += f64::from(neighbour.peaks[index].intensity);
            pos_score += mz_score;
            matches += 1;
            if pattern.peak[pattern_index] == PatternPeak::NotFound {
                pattern.peak[pattern_index] = PatternPeak::Found(index);
                pattern.spectrum[pattern_index] = neighbour_index;
            }
        }
    }
    if matches == 0 {
        put_all(log, &[" missing\n"]);
        pattern.peak[pattern_index] = PatternPeak::NotFound;
        pattern.mz_score[pattern_index] = 0.0;
        pattern.intensity[pattern_index] = 0.0;
    } else {
        if log.enabled() {
            put_all(log, &["=> ", &g(intensity / f64::from(matches)), "\n"]);
        }
        pattern.mz_score[pattern_index] = pos_score / f64::from(matches);
        pattern.intensity[pattern_index] = intensity / f64::from(matches);
    }
    Ok(work)
}

/// A score between 0 and 1 for the correlation of a theoretical and a found
/// isotope pattern: source `isotopeScore_`.
///
/// 1. If a required isotope, one between `optional_begin` and `len -
///    optional_end`, is [`PatternPeak::NotFound`], the score is 0.
/// 2. The fit may leave out optional isotopes at either end, but no gap: the
///    search starts behind the last missing optional isotope at each end.
/// 3. For every candidate `b` leading and `e` trailing isotopes left out, with
///    more than two isotopes left (or exactly two for the starting candidate),
///    the Pearson correlation of the theoretical and found intensities is
///    computed; NaN counts as 0 and a two-isotope fit is capped at
///    `min_isotope_fit`. A candidate replaces the best when its score divided by
///    the best is at least `1 + optional_fit_improvement`; the best starts at
///    0.01.
/// 4. If the best fit leaves no isotope, the score is 0. Otherwise the left-out
///    isotopes become [`PatternPeak::Removed`] with zero intensity and m/z
///    score, and with `consider_mz_distances` the score is multiplied by the
///    mean m/z score of the kept isotopes.
///
/// The source reads the start of the inner candidate loop from the current best
/// trailing count on every outer iteration, so a new best fit narrows the
/// candidates of the following outer iterations; this does the same.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the pattern's vectors do not all have
/// the length of `isotopes`, or the optional counts exceed that length. The
/// source assumes both.
pub fn isotope_score(
    isotopes: &TheoreticalIsotopePattern,
    pattern: &mut IsotopePattern,
    consider_mz_distances: bool,
    min_isotope_fit: f64,
    optional_fit_improvement: f64,
) -> Result<f64> {
    Ok(isotope_score_with_work(
        isotopes,
        pattern,
        consider_mz_distances,
        min_isotope_fit,
        optional_fit_improvement,
    )?
    .0)
}

pub(crate) fn isotope_score_with_work(
    isotopes: &TheoreticalIsotopePattern,
    pattern: &mut IsotopePattern,
    consider_mz_distances: bool,
    min_isotope_fit: f64,
    optional_fit_improvement: f64,
) -> Result<(f64, u64)> {
    isotope_score_logged(
        isotopes,
        pattern,
        consider_mz_distances,
        min_isotope_fit,
        optional_fit_improvement,
        &mut NoLog,
    )
}

/// [`isotope_score`] writing the source's debug lines to `log`: the number of
/// peaks, a missing core peak, the starting `best_begin/end`, and every
/// candidate fit with ` - new best fit ` when it wins.
pub(crate) fn isotope_score_logged<L: LogSink>(
    isotopes: &TheoreticalIsotopePattern,
    pattern: &mut IsotopePattern,
    consider_mz_distances: bool,
    min_isotope_fit: f64,
    optional_fit_improvement: f64,
    log: &mut L,
) -> Result<(f64, u64)> {
    let size = isotopes.len();
    if pattern.peak.len() != size
        || pattern.intensity.len() != size
        || pattern.mz_score.len() != size
        || isotopes
            .optional_begin
            .checked_add(isotopes.optional_end)
            .is_none_or(|optional| optional > size)
    {
        return Err(Error::InvalidValue(
            "isotope pattern and theoretical pattern differ in length".into(),
        ));
    }
    let mut work = 0u64;
    if log.enabled() {
        put_all(
            log,
            &[
                "   - fitting ",
                &pattern.intensity.len().to_string(),
                " peaks\n",
            ],
        );
    }
    for iso in isotopes.optional_begin..size - isotopes.optional_end {
        if pattern.peak[iso] == PatternPeak::NotFound {
            put_all(log, &["   - aborting: core peak is missing\n"]);
            return Ok((0.0, work));
        }
    }
    let mut best_int_score = 0.01;
    let mut best_begin = 0;
    for i in (1..=isotopes.optional_begin).rev() {
        if pattern.peak[i - 1] == PatternPeak::NotFound {
            best_begin = i;
            break;
        }
    }
    let mut best_end = 0;
    for i in (1..=isotopes.optional_end).rev() {
        if pattern.peak[size - i] == PatternPeak::NotFound {
            best_end = i;
            break;
        }
    }
    if log.enabled() {
        put_all(
            log,
            &[
                "   - best_begin/end: ",
                &best_begin.to_string(),
                "/",
                &best_end.to_string(),
                "\n",
            ],
        );
    }
    let first_begin = best_begin;
    for b in first_begin..=isotopes.optional_begin {
        let mut e = best_end;
        while e <= isotopes.optional_end {
            let kept = size - b - e;
            if kept > 2 || (b == best_begin && e == best_end && kept > 1) {
                work += kept as u64;
                let mut int_score = pearson_correlation_coefficient(
                    &isotopes.intensity[b..size - e],
                    &pattern.intensity[b..size - e],
                )?;
                if int_score.is_nan() {
                    int_score = 0.0;
                }
                if kept == 2 && int_score > min_isotope_fit {
                    int_score = min_isotope_fit;
                }
                if log.enabled() {
                    put_all(
                        log,
                        &[
                            "   - fit (",
                            &b.to_string(),
                            "/",
                            &e.to_string(),
                            "): ",
                            &g(int_score),
                        ],
                    );
                }
                if int_score / best_int_score >= 1.0 + optional_fit_improvement {
                    put_all(log, &[" - new best fit "]);
                    best_int_score = int_score;
                    best_begin = b;
                    best_end = e;
                }
                put_all(log, &["\n"]);
            }
            e += 1;
        }
    }
    if size - best_begin - best_end == 0 {
        return Ok((0.0, work));
    }
    for i in 0..best_begin {
        pattern.peak[i] = PatternPeak::Removed;
        pattern.intensity[i] = 0.0;
        pattern.mz_score[i] = 0.0;
    }
    for i in 0..best_end {
        pattern.peak[size - 1 - i] = PatternPeak::Removed;
        pattern.intensity[size - 1 - i] = 0.0;
        pattern.mz_score[size - 1 - i] = 0.0;
    }
    if consider_mz_distances {
        let kept = &pattern.mz_score[best_begin..size - best_end];
        let sum = kept.iter().fold(0.0, |total, value| total + value);
        best_int_score *= sum / kept.len() as f64;
    }
    Ok((best_int_score, work))
}

/// Reset `pattern` to `size` isotopes with nothing matched, reusing its
/// allocations: the state of a fresh source `IsotopePattern(size)`.
pub(crate) fn reset_pattern(pattern: &mut IsotopePattern, size: usize) {
    pattern.peak.clear();
    pattern.peak.resize(size, PatternPeak::NotFound);
    pattern.spectrum.clear();
    pattern.spectrum.resize(size, 0);
    pattern.intensity.clear();
    pattern.intensity.resize(size, 0.0);
    pattern.mz_score.clear();
    pattern.mz_score.resize(size, 0.0);
    pattern.theoretical_mz.clear();
    pattern.theoretical_mz.resize(size, 0.0);
}
