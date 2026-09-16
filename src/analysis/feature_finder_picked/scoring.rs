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
//! and written out only by the debug mode, which the port refuses. They live in
//! [`ScoreArrays`](crate::analysis::feature_finder_picked::scoring::ScoreArrays)
//! instead, one flat `f32` array per score, so the input spectra
//! keep their own data arrays and no per-spectrum allocation is needed.
//!
//! Arithmetic follows the source operation by operation, including the `float`
//! narrowing of every stored score. See `docs/FEATURE_FINDER_PICKED_SUPPORT.md`.

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
    /// A step can be zero (every MS1 spectrum at one retention time, every MS1
    /// peak at one m/z, or an extent so small that the division underflows)
    /// or, for the retention time, infinite (an extent that overflows). The
    /// source computes the bins anyway (`FeatureFinderAlgorithmPicked.cpp:244-277`),
    /// and so does this: with a zero step every cell spans the whole extent,
    /// and with an infinite step the first cell's bounds are `NaN` and `inf`,
    /// which the inclusive walk, like the source's area iterator, reads as
    /// every spectrum, while the other cells start at `inf` and are empty. The
    /// executed Linux x86_64 Release build computes the same quantiles for
    /// all these grids. Scoring a peak on such a grid is where the source
    /// becomes undefined; see [`Self::score`] and
    /// [`DegenerateBinStep`](crate::analysis::feature_finder_picked::algorithm::DegenerateBinStep).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `bins` is zero, when a spectrum is
    /// not MS1 or holds a non-finite value, or when there is no MS1 peak, and
    /// [`Error::UnsortedData`] when the spectra are not sorted by RT and m/z.
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
        let (rt_step, mz_step) = bin_steps(&rt, &mz, bins);
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
    /// `[0, 1]`. A NaN from `0 / 0` (a zero intensity against a zero first
    /// quantile) or from a NaN intensity passes the clamp unchanged, as in the
    /// source, with the sign and payload the Linux x86_64 Release build gives
    /// it: the default NaN for `0 / 0`, the intensity's own NaN otherwise.
    ///
    /// Returns `None` for a cell outside the grid.
    pub fn bin_score(&self, rt_bin: usize, mz_bin: usize, intensity: f64) -> Option<f64> {
        let quantiles = self.quantiles(rt_bin, mz_bin)?;
        let position = quantiles.partition_point(|&quantile| quantile < intensity);
        let Some(&upper) = quantiles.get(position) else {
            return Some(1.0);
        };
        // The operations follow the source; a NaN result is the one the Linux
        // x86_64 Release build returns (see `x86_64`): a NaN intensity
        // propagates, and `0 / 0` is the default NaN.
        let bin_score = if position == 0 {
            x86_64::div(x86_64::mul(0.05, intensity), upper)
        } else {
            let lower = quantiles[position - 1];
            x86_64::div(
                x86_64::mul(0.05, x86_64::sub(intensity, lower)),
                upper - lower,
            )
        };
        // `clamp` keeps NaN and -0.0, as the source's two comparisons do.
        Some(x86_64::add(bin_score, 0.05 * (position as f64 - 1.0)).clamp(0.0, 1.0))
    }

    /// The intensity score of a peak: source `intensityScore_(spectrum, peak)`.
    ///
    /// The peak's position on a half-bin grid, `floor((x - start) / step * 2)`
    /// converted to `UInt` and capped at `2 * bins - 1`, selects the two nearest
    /// bins per dimension (one at the outer half-bins). The four cell scores of
    /// [`Self::bin_score`] are weighted by `d = sqrt((1 - d_rt)^2 + (1 -
    /// d_mz)^2)`, where `d_rt` and `d_mz` are the distances of the peak to each
    /// bin centre in bin widths, and each weight is divided by the sum of the
    /// four. The squares are products where the source calls `std::pow(x, 2)`,
    /// which GCC compiles as a product; the executed intensity scores with 1, 2,
    /// 7 and 10 bins agree with them bit for bit.
    ///
    /// # Undefined behaviour of the source, reproduced as the Linux x86_64 Release build computes it
    ///
    /// The conversion `(UInt) std::floor(...)` (`FeatureFinderAlgorithmPicked.cpp:1837-1838`)
    /// is undefined in C++ when the floored half-bin position is NaN,
    /// infinite, negative, or `2^32` and above. That happens for a peak
    /// outside the binned range, which the algorithm never scores, and for
    /// every peak when a bin step is zero or infinite
    /// ([`DegenerateBinStep`](crate::analysis::feature_finder_picked::algorithm::DegenerateBinStep)).
    /// GCC 14 at `-O3` compiles the conversion in `libOpenMS.so` as
    /// `cvttsd2si %xmm2,%rdi` (64-bit truncation, `0x8000000000000000` for NaN
    /// and for values outside the signed 64-bit range) followed by the low 32
    /// bits of the register, and the cap as an unsigned `cmovbe`. This
    /// function reproduces exactly that (the crate-private
    /// `x86_64::truncate_to_u32`): a negative position wraps to a large
    /// unsigned value and is capped at the last half-bin, and a NaN or
    /// infinite position selects half-bin 0. The
    /// cells read are then always inside the grid. With a zero or infinite step
    /// the distances are `0 / 0` or `inf / inf`, so the score is NaN whatever
    /// half-bin was selected.
    ///
    /// NaN results carry the sign and payload that build returns: each
    /// arithmetic step follows SSE2's NaN rule (the first NaN operand of the
    /// emitted instruction, quieted; the default NaN `0xfff8000000000000` for an
    /// invalid operation) in the operand order of the emitted instructions,
    /// which swaps some commutative operands of the source text. The executed
    /// probe `iscore_probe` (57 positions, including NaN and infinite ones, on
    /// four grids) and every peak of the degenerate captures agree with this
    /// bit for bit, NaN bits included.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] only when `2 * bins - 1` does not fit the
    /// source's `UInt`, which [`Self::compute`] cannot produce in practice
    /// because the quantile table of such a grid cannot be allocated.
    pub fn score(&self, rt: f64, mz: f64, intensity: f64) -> Result<f64> {
        use x86_64::{abs, add, div, mul, sqrt, sub};
        let last = u32::try_from(self.bins)
            .ok()
            .and_then(|bins| bins.checked_mul(2))
            .map(|twice| twice - 1)
            .ok_or_else(|| {
                Error::InvalidValue("intensity:bins exceeds the source's UInt".into())
            })?;
        // `floor((x - start) / step * 2.0)`; GCC emits the doubling as `q + q`.
        let half_bin = |x: f64, start: f64, step: f64| -> u32 {
            let quotient = div(sub(x, start), step);
            let position = add(quotient, quotient);
            last.min(x86_64::truncate_to_u32(position.floor()))
        };
        let rt_bin = half_bin(rt, self.rt_start, self.rt_step);
        let mz_bin = half_bin(mz, self.mz_start, self.mz_step);
        let neighbours = |bin: u32| -> (u32, u32) {
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
        // `|start + (0.5 + b) * step - x| / step`, in the emitted order.
        let distance = |b: u32, start: f64, step: f64, x: f64| -> f64 {
            div(
                abs(sub(add(mul(add(f64::from(b), 0.5), step), start), x)),
                step,
            )
        };
        let drl = distance(rl, self.rt_start, self.rt_step, rt);
        let drh = distance(rh, self.rt_start, self.rt_step, rt);
        let dml = distance(ml, self.mz_start, self.mz_step, mz);
        let dmh = distance(mh, self.mz_start, self.mz_step, mz);
        let square = |x: f64| mul(x, x);
        let a = square(sub(1.0, drl));
        let b = square(sub(1.0, drh));
        let c = square(sub(1.0, dml));
        let d = square(sub(1.0, dmh));
        // Source `d1 = sqrt(a + c)`, `d2 = sqrt(b + c)`, `d3 = sqrt(a + d)`,
        // `d4 = sqrt(b + d)`; GCC emits `c + a`, `c + b`, `a + d` and `d + b`.
        let d1 = sqrt(add(c, a));
        let d2 = sqrt(add(c, b));
        let d3 = sqrt(add(a, d));
        let d4 = sqrt(add(d, b));
        let d_sum = add(add(add(d1, d2), d3), d4);
        let cell = |r: u32, m: u32| {
            self.bin_score(r as usize, m as usize, intensity)
                .ok_or_else(|| Error::InvalidValue("intensity bin outside the grid".into()))
        };
        // Source `c1 * (d1 / d_sum) + c2 * (d2 / d_sum) + ...`, left to right;
        // GCC emits each product with the weight first and the first two sums
        // with the new term first.
        let t1 = mul(div(d1, d_sum), cell(rl, ml)?);
        let t2 = mul(div(d2, d_sum), cell(rh, ml)?);
        let t3 = mul(div(d3, d_sum), cell(rl, mh)?);
        let t4 = mul(div(d4, d_sum), cell(rh, mh)?);
        Ok(add(add(t3, add(t2, t1)), t4))
    }
}

/// Emulation of the instructions the Linux x86_64 Release build of
/// `libOpenMS.so` (`openms4-release-bc9cc12-c19e494-174b576`, GCC 14.4, `-O3
/// -mssse3 -ffp-contract=off`) emits for `intensityScore_`, where IEEE 754
/// alone does not fix the result.
///
/// The arithmetic helpers return the IEEE result whenever it is not NaN, so
/// they never change a number. For NaN they follow the SSE2 rule (Intel SDM
/// vol. 1, "Rules for handling NaNs"): the first NaN operand of the
/// instruction, quieted, or the default NaN `0xfff8000000000000` when the
/// operation itself is invalid. Rust's own arithmetic leaves those bits to the
/// host, and an arm64 host produces the positive default NaN.
pub(crate) mod x86_64 {
    /// SSE2's default ("real indefinite") NaN.
    pub(crate) const DEFAULT_NAN: f64 = f64::from_bits(0xfff8_0000_0000_0000);

    fn quiet(x: f64) -> f64 {
        f64::from_bits(x.to_bits() | 0x0008_0000_0000_0000)
    }

    /// The NaN rule of a two-operand SSE2 instruction whose destination
    /// register holds `first`.
    fn nan_rule(first: f64, second: f64, result: f64) -> f64 {
        if !result.is_nan() {
            result
        } else if first.is_nan() {
            quiet(first)
        } else if second.is_nan() {
            quiet(second)
        } else {
            DEFAULT_NAN
        }
    }

    /// `addsd`: `first + second`.
    pub(crate) fn add(first: f64, second: f64) -> f64 {
        nan_rule(first, second, first + second)
    }

    /// `subsd`: `first - second`.
    pub(crate) fn sub(first: f64, second: f64) -> f64 {
        nan_rule(first, second, first - second)
    }

    /// `mulsd`: `first * second`.
    pub(crate) fn mul(first: f64, second: f64) -> f64 {
        nan_rule(first, second, first * second)
    }

    /// `divsd`: `first / second`.
    pub(crate) fn div(first: f64, second: f64) -> f64 {
        nan_rule(first, second, first / second)
    }

    /// `sqrtsd`.
    pub(crate) fn sqrt(x: f64) -> f64 {
        nan_rule(x, x, x.sqrt())
    }

    /// `andpd` with the absolute-value mask: clears the sign bit, of a NaN too.
    pub(crate) fn abs(x: f64) -> f64 {
        f64::from_bits(x.to_bits() & !(1 << 63))
    }

    /// `cvttsd2si %xmm, %r64` followed by the low 32 bits of the register: how
    /// the Release build converts a `double` to `UInt`.
    ///
    /// The instruction truncates towards zero and returns
    /// `0x8000000000000000` for NaN and for every value outside the signed
    /// 64-bit range; the low 32 bits of that are 0. An in-range value keeps its
    /// low 32 bits, so `-1.0` becomes `0xffffffff` and `2^32 + 2` becomes 2.
    pub(crate) fn truncate_to_u32(x: f64) -> u32 {
        const TWO_TO_63: f64 = 9_223_372_036_854_775_808.0;
        let wide = if x.is_nan() || !(-TWO_TO_63..TWO_TO_63).contains(&x) {
            i64::MIN
        } else {
            x as i64
        };
        wide as u32
    }

    /// `cvtsd2ss`: `x` narrowed to `f32`. A NaN keeps its sign and the top 22
    /// bits of its payload and is quieted, as the instruction does; Rust's `as`
    /// leaves NaN bits to the host.
    pub(crate) fn narrow(x: f64) -> f32 {
        if x.is_nan() {
            let bits = x.to_bits();
            let sign = ((bits >> 32) as u32) & 0x8000_0000;
            let payload = ((bits >> 29) as u32) & 0x003f_ffff;
            f32::from_bits(sign | 0x7fc0_0000 | payload)
        } else {
            x as f32
        }
    }
}

/// The RT and m/z bin steps of step 1: the extents divided by `bins`
/// (`FeatureFinderAlgorithmPicked.cpp:244-245`).
pub(crate) fn bin_steps(rt: &NumericRange, mz: &NumericRange, bins: usize) -> (f64, f64) {
    let bins_f = bins as f64;
    ((rt.max - rt.min) / bins_f, (mz.max - mz.min) / bins_f)
}

/// Whether a bin step makes every intensity score undefined in the source:
/// zero, or not finite.
pub(crate) fn degenerate_step(step: f64) -> bool {
    !(step > 0.0 && step.is_finite())
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
#[derive(Clone, Copy, Debug)]
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
/// As [`IntensityThresholds::compute`].
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
            target[index] = x86_64::narrow(thresholds.score(
                spectrum.rt,
                peak.mz,
                f64::from(peak.intensity),
            )?);
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
        pattern.peak[pattern_index] = PatternPeak::Found(nearest_index);
        pattern.spectrum[pattern_index] = spectrum_index;
        intensity += f64::from(spectrum.peaks[nearest_index].intensity);
        pos_score += this_mz_score;
        matches += 1;
    }
    let neighbours = [
        spectrum_index.checked_sub(1),
        spectrum_index
            .checked_add(1)
            .filter(|&index| index < spectra.len()),
    ];
    for neighbour_index in neighbours.into_iter().flatten() {
        let neighbour = &spectra[neighbour_index];
        let Some(index) = nearest(&neighbour.peaks, pos, |peak| peak.mz) else {
            continue;
        };
        work += 1;
        let mz_score = position_score(pos, neighbour.peaks[index].mz, pattern_tolerance);
        if mz_score != 0.0 {
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
        pattern.peak[pattern_index] = PatternPeak::NotFound;
        pattern.mz_score[pattern_index] = 0.0;
        pattern.intensity[pattern_index] = 0.0;
    } else {
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
    for iso in isotopes.optional_begin..size - isotopes.optional_end {
        if pattern.peak[iso] == PatternPeak::NotFound {
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
                if int_score / best_int_score >= 1.0 + optional_fit_improvement {
                    best_int_score = int_score;
                    best_begin = b;
                    best_end = e;
                }
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
