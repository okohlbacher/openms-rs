// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Chris Bielow, OpenMS Rust contributors $

//! Iteratively clipped histogram-mean signal-to-noise estimation.
//!
//! Ports the header-only template
//! `src/openms/include/OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h`.
//! See docs/MEAN_NOISE_SUPPORT.md for the source's fixed denominator and
//! historical percentile conventions. The source's `DefaultParamHandler`
//! parameters map through [`DefaultParamHandling`], implemented on
//! [`SignalToNoiseEstimatorMeanIterative`].

use crate::kernel::{MSChromatogram, MSSpectrum};
use crate::param::handler::int_parameter;
use crate::param::{DefaultParamHandling, Param, ParamValue};
use crate::processing::peak_picking::NoiseRangeParameters;
use crate::{Error, Result};

/// How the histogram's upper intensity bound is found (source enum
/// `IntensityThresholdCalculation`, chosen by the parameter `auto_mode`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MeanNoiseHistogramRange {
    /// Source `AUTOMAXBYSTDEV` (`auto_mode` 0): mean plus `factor` standard
    /// deviations of the whole trace (`auto_max_stdev_factor`).
    StandardDeviation {
        /// `auto_max_stdev_factor`, in `0..=999`.
        factor: f64,
    },
    /// Source `MANUAL` (`auto_mode` -1): the fixed bound `max_intensity`,
    /// which must be positive. Intensities at or above it are not counted.
    Manual {
        /// `max_intensity`; the source parameter is an integer.
        max_intensity: f64,
    },
    /// Literal historical source heuristic, including its minimum-as-maximum
    /// and bounded histogram scan. Unsafe source indices return errors.
    /// This mode retains the source's f32 intensity arithmetic.
    ///
    /// Source `AUTOMAXBYPERCENT` (`auto_mode` 1).
    LegacyPercentile {
        /// `auto_max_percentile`, in `0..=100`.
        percentile: u8,
    },
}
impl Default for MeanNoiseHistogramRange {
    fn default() -> Self {
        Self::StandardDeviation { factor: 3.0 }
    }
}

/// Three histogram clipping passes with the original window count as denominator.
///
/// Every field but the three ceilings is a source parameter; see
/// [`DefaultParamHandling`] for the names, defaults and restrictions.
#[derive(Clone, Debug, PartialEq)]
pub struct SignalToNoiseEstimatorMeanIterative {
    /// The histogram's upper bound (source `auto_mode` with
    /// `auto_max_stdev_factor`, `max_intensity` or `auto_max_percentile`).
    pub histogram_range: MeanNoiseHistogramRange,
    /// The values of the three range parameters as configured. The mode in
    /// [`histogram_range`](Self::histogram_range) takes its own value from
    /// there; the other two are kept only so that
    /// [`DefaultParamHandling::to_param`] returns the parameters the estimator
    /// was built from.
    pub range_parameters: NoiseRangeParameters,
    /// Full coordinate width, at least one; left inclusive and right exclusive.
    /// Source `win_len`, default 200.
    pub window_length: f64,
    /// Histogram bins, at least three (source `bin_count`, default 30).
    pub bin_count: usize,
    /// Standard deviations above the mean at which each clipping pass cuts the
    /// histogram, in `0.01..=999` (source `stdev_mp`, default 3).
    pub stdev_multiplier: f64,
    /// Points a window needs before it is estimated rather than treated as
    /// sparse, at least one (source `min_required_elements`, default 10).
    pub min_required_elements: usize,
    /// The noise assigned to a sparse window, positive and finite (source
    /// `noise_for_empty_window`, default 1e20).
    pub noise_for_empty_window: f64,
    /// Native ceiling on the points of one trace; the source has none.
    pub max_points: usize,
    /// Native ceiling on [`bin_count`](Self::bin_count); the source has none.
    pub max_bins: usize,
    /// Output points, global intensity scans, histogram slots/updates and bin visits.
    pub max_work: usize,
}
impl Default for SignalToNoiseEstimatorMeanIterative {
    fn default() -> Self {
        Self {
            histogram_range: Default::default(),
            range_parameters: NoiseRangeParameters::default(),
            window_length: 200.0,
            bin_count: 30,
            stdev_multiplier: 3.0,
            min_required_elements: 10,
            noise_for_empty_window: 1e20,
            max_points: 1_000_000,
            max_bins: 1_000_000,
            max_work: 50_000_000,
        }
    }
}

/// The estimates of one trace, in input order (source `stn_estimates_` with
/// the diagnostics its `computeSTN_` logs).
#[derive(Clone, Debug, PartialEq)]
pub struct MeanNoiseEstimates {
    /// Signal-to-noise ratio of every point (source `getSignalToNoise`).
    pub signal_to_noise: Vec<f64>,
    /// The noise estimate of every point's window.
    pub noise: Vec<f64>,
    /// The histogram's upper bound the configured range produced.
    pub max_intensity: f64,
    /// Percentage of points whose window was sparse.
    pub sparse_window_percent: f64,
}

/// The source parameters of `SignalToNoiseEstimatorMeanIterative`, declared by
/// its constructor (`SignalToNoiseEstimatorMeanIterative.h:65-106`) and read by
/// `updateMembers_` (`:366-377`). The base class is constructed as
/// `DefaultParamHandler("SignalToNoiseEstimator")` and renamed by
/// `setName("SignalToNoiseEstimatorMeanIterative")`.
impl DefaultParamHandling for SignalToNoiseEstimatorMeanIterative {
    const NAME: &'static str = "SignalToNoiseEstimatorMeanIterative";

    /// `max_intensity` (-1, min -1), `auto_max_stdev_factor` (3.0, 0..999),
    /// `auto_max_percentile` (95, 0..100), `auto_mode` (0, -1..1), `win_len`
    /// (200.0, min 1), `bin_count` (30, min 3), `stdev_mp` (3.0, 0.01..999),
    /// `min_required_elements` (10, min 1) and `noise_for_empty_window` (1e20);
    /// all but `win_len`, `bin_count` and `min_required_elements` advanced.
    fn defaults() -> Result<Param> {
        let advanced = ["advanced".to_owned()];
        let mut defaults = Param::new();
        defaults.set_value(
            "max_intensity",
            ParamValue::Integer(-1),
            "maximal intensity considered for histogram construction. By default, it will be calculated automatically (see auto_mode). Only provide this parameter if you know what you are doing (and change 'auto_mode' to '-1')! All intensities EQUAL/ABOVE 'max_intensity' will not be added to the histogram. If you choose 'max_intensity' too small, the noise estimate might be too small as well. If chosen too big, the bins become quite large (which you could counter by increasing 'bin_count', which increases runtime).",
            &advanced,
        )?;
        defaults.set_min_int("max_intensity", -1)?;
        defaults.set_value(
            "auto_max_stdev_factor",
            ParamValue::Float(3.0),
            "parameter for 'max_intensity' estimation (if 'auto_mode' == 0): mean + 'auto_max_stdev_factor' * stdev",
            &advanced,
        )?;
        defaults.set_min_float("auto_max_stdev_factor", 0.0)?;
        defaults.set_max_float("auto_max_stdev_factor", 999.0)?;
        defaults.set_value(
            "auto_max_percentile",
            ParamValue::Integer(95),
            "parameter for 'max_intensity' estimation (if 'auto_mode' == 1): auto_max_percentile th percentile",
            &advanced,
        )?;
        defaults.set_min_int("auto_max_percentile", 0)?;
        defaults.set_max_int("auto_max_percentile", 100)?;
        defaults.set_value(
            "auto_mode",
            ParamValue::Integer(0),
            "method to use to determine maximal intensity: -1 --> use 'max_intensity'; 0 --> 'auto_max_stdev_factor' method (default); 1 --> 'auto_max_percentile' method",
            &advanced,
        )?;
        defaults.set_min_int("auto_mode", -1)?;
        defaults.set_max_int("auto_mode", 1)?;
        defaults.set_value(
            "win_len",
            ParamValue::Float(200.0),
            "window length in Thomson",
            &[],
        )?;
        defaults.set_min_float("win_len", 1.0)?;
        defaults.set_value(
            "bin_count",
            ParamValue::Integer(30),
            "number of bins for intensity values",
            &[],
        )?;
        defaults.set_min_int("bin_count", 3)?;
        defaults.set_value(
            "stdev_mp",
            ParamValue::Float(3.0),
            "multiplier for stdev",
            &advanced,
        )?;
        defaults.set_min_float("stdev_mp", 0.01)?;
        defaults.set_max_float("stdev_mp", 999.0)?;
        defaults.set_value(
            "min_required_elements",
            ParamValue::Integer(10),
            "minimum number of elements required in a window (otherwise it is considered sparse)",
            &[],
        )?;
        defaults.set_min_int("min_required_elements", 1)?;
        defaults.set_value(
            "noise_for_empty_window",
            ParamValue::Float(1e20),
            "noise value used for sparse windows",
            &advanced,
        )?;
        Ok(defaults)
    }

    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a manual `max_intensity` is not an integer
    /// in the 32-bit range, or a count exceeds the source's 32-bit `int`.
    fn to_param(&self) -> Result<Param> {
        let mut range = self.range_parameters;
        let mode = match self.histogram_range {
            MeanNoiseHistogramRange::StandardDeviation { factor } => {
                range.auto_max_stdev_factor = factor;
                0
            }
            MeanNoiseHistogramRange::LegacyPercentile { percentile } => {
                range.auto_max_percentile = i32::from(percentile);
                1
            }
            MeanNoiseHistogramRange::Manual { max_intensity } => {
                let integral = max_intensity.fract() == 0.0
                    && max_intensity >= f64::from(i32::MIN)
                    && max_intensity <= f64::from(i32::MAX);
                if !integral {
                    return Err(bad(&format!(
                        "'max_intensity' = {max_intensity} is not representable as the source integer parameter"
                    )));
                }
                range.max_intensity = max_intensity as i32;
                -1
            }
        };
        let mut param = Self::defaults()?;
        param.replace_value(
            "max_intensity",
            ParamValue::Integer(i64::from(range.max_intensity)),
        )?;
        param.replace_value(
            "auto_max_stdev_factor",
            ParamValue::Float(range.auto_max_stdev_factor),
        )?;
        param.replace_value(
            "auto_max_percentile",
            ParamValue::Integer(i64::from(range.auto_max_percentile)),
        )?;
        param.replace_value("auto_mode", ParamValue::Integer(mode))?;
        param.replace_value("win_len", ParamValue::Float(self.window_length))?;
        param.replace_value("bin_count", int_parameter(self.bin_count, "bin_count")?)?;
        param.replace_value("stdev_mp", ParamValue::Float(self.stdev_multiplier))?;
        param.replace_value(
            "min_required_elements",
            int_parameter(self.min_required_elements, "min_required_elements")?,
        )?;
        param.replace_value(
            "noise_for_empty_window",
            ParamValue::Float(self.noise_for_empty_window),
        )?;
        Ok(param)
    }

    /// Source `updateMembers_`. `auto_mode` 0 selects
    /// [`MeanNoiseHistogramRange::StandardDeviation`], 1
    /// [`MeanNoiseHistogramRange::LegacyPercentile`] and -1 (the only other
    /// value the restriction admits) [`MeanNoiseHistogramRange::Manual`], as
    /// `computeSTN_` branches on it. A manual bound that is not positive is
    /// accepted here, as by the source, and refused by estimation (the source
    /// throws there too).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a count below zero, which the restrictions
    /// already refuse.
    fn update_members(param: &Param) -> Result<Self> {
        let range = NoiseRangeParameters {
            max_intensity: param.value("max_intensity")?.to_i32()?,
            auto_max_stdev_factor: param.value("auto_max_stdev_factor")?.to_f64()?,
            auto_max_percentile: param.value("auto_max_percentile")?.to_i32()?,
        };
        let histogram_range = match param.value("auto_mode")?.to_i64()? {
            0 => MeanNoiseHistogramRange::StandardDeviation {
                factor: range.auto_max_stdev_factor,
            },
            1 => MeanNoiseHistogramRange::LegacyPercentile {
                percentile: u8::try_from(range.auto_max_percentile)
                    .map_err(|_| bad("'auto_max_percentile' must be in 0..=100"))?,
            },
            _ => MeanNoiseHistogramRange::Manual {
                max_intensity: f64::from(range.max_intensity),
            },
        };
        let count = |key: &str| -> Result<usize> {
            usize::try_from(param.value(key)?.to_i64()?)
                .map_err(|_| bad(&format!("parameter '{key}' must not be negative")))
        };
        Ok(Self {
            histogram_range,
            range_parameters: range,
            window_length: param.value("win_len")?.to_f64()?,
            bin_count: count("bin_count")?,
            stdev_multiplier: param.value("stdev_mp")?.to_f64()?,
            min_required_elements: count("min_required_elements")?,
            noise_for_empty_window: param.value("noise_for_empty_window")?.to_f64()?,
            ..Self::default()
        })
    }

    /// Keeps [`max_points`](Self::max_points), [`max_bins`](Self::max_bins)
    /// and [`max_work`](Self::max_work).
    fn keep_native_members(&self, updated: Self) -> Self {
        Self {
            max_points: self.max_points,
            max_bins: self.max_bins,
            max_work: self.max_work,
            ..updated
        }
    }
}

impl SignalToNoiseEstimatorMeanIterative {
    /// Check the options before any estimation.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a window below one or not finite, fewer than
    /// three bins or more than [`max_bins`](Self::max_bins), a zero ceiling,
    /// `min_required_elements` of zero, a `stdev_multiplier` outside
    /// `0.01..=999`, a sparse-window noise that is not positive and finite, a
    /// standard-deviation factor outside `0..=999`, a manual bound that is not
    /// positive and finite (the source's `computeSTN_` throws
    /// `Exception::InvalidValue` for it) and a percentile above 100.
    pub fn validate(&self) -> Result<()> {
        if !self.window_length.is_finite()
            || self.window_length < 1.0
            || self.bin_count < 3
            || self.bin_count > self.max_bins
            || self.bin_count > i32::MAX as usize
            || self.max_points == 0
            || self.max_work == 0
            || self.min_required_elements == 0
            || !self.stdev_multiplier.is_finite()
            || !(0.01..=999.0).contains(&self.stdev_multiplier)
            || !self.noise_for_empty_window.is_finite()
            || self.noise_for_empty_window <= 0.0
        {
            return Err(bad("invalid iterative mean-noise options or limits"));
        }
        match self.histogram_range {
            MeanNoiseHistogramRange::StandardDeviation { factor }
                if !factor.is_finite() || !(0.0..=999.0).contains(&factor) =>
            {
                Err(bad(
                    "mean-noise standard-deviation factor must be in 0..=999",
                ))
            }
            MeanNoiseHistogramRange::Manual { max_intensity }
                if !max_intensity.is_finite() || max_intensity <= 0.0 =>
            {
                Err(bad(
                    "manual mean-noise histogram maximum must be positive and finite",
                ))
            }
            MeanNoiseHistogramRange::LegacyPercentile { percentile } if percentile > 100 => {
                Err(bad("legacy noise percentile must be in 0..=100"))
            }
            _ => Ok(()),
        }
    }

    /// Estimate an ordered trace without modifying it. Coincident positions and
    /// finite signed intensities are supported; only bin assignment clamps negatives.
    pub fn estimate(&self, positions: &[f64], intensities: &[f64]) -> Result<MeanNoiseEstimates> {
        self.preflight(positions.len())?;
        if positions.len() != intensities.len()
            || positions.iter().chain(intensities).any(|v| !v.is_finite())
        {
            return Err(bad(
                "mean-noise coordinates and intensities must be finite and aligned",
            ));
        }
        if positions.windows(2).any(|v| v[0] > v[1]) {
            return Err(Error::UnsortedData);
        }
        let mut budget = Budget {
            remaining: self.max_work,
        };
        let max_intensity = self.maximum(intensities, &mut budget)?;
        if !max_intensity.is_finite() || max_intensity < 0.0 {
            return Err(bad("mean-noise histogram maximum is negative or nonfinite"));
        }
        let n = positions.len();
        budget.charge(n)?;
        if n != 0 {
            budget.charge(self.bin_count)?;
        }
        let mut result = MeanNoiseEstimates {
            signal_to_noise: Vec::with_capacity(n),
            noise: Vec::with_capacity(n),
            max_intensity,
            sparse_window_percent: 0.0,
        };
        if n == 0 {
            return Ok(result);
        }
        let bin_size = (max_intensity / self.bin_count as f64).max(1.0);
        let bin_value: Vec<_> = (0..self.bin_count)
            .map(|i| finite((i as f64 + 0.5) * bin_size))
            .collect::<Result<_>>()?;
        let mut histogram = vec![0usize; self.bin_count];
        let to_bin = |y: f64| {
            let b = y.max(0.0) / bin_size;
            (b < self.bin_count as f64).then_some(b as usize)
        };
        let (mut left, mut right, mut count, mut sparse) = (0, 0, 0usize, 0usize);
        for &position in positions {
            let low = finite(position - self.window_length / 2.0)?;
            let high = finite(position + self.window_length / 2.0)?;
            if high <= position {
                return Err(bad("mean-noise window cannot advance the upper coordinate"));
            }
            while left < right && positions[left] < low {
                budget.charge(1)?;
                if let Some(bin) = to_bin(intensities[left]) {
                    histogram[bin] = histogram[bin]
                        .checked_sub(1)
                        .ok_or_else(|| bad("mean-noise histogram underflow"))?;
                    count -= 1;
                }
                left += 1;
            }
            while right < n && positions[right] < high {
                budget.charge(1)?;
                if let Some(bin) = to_bin(intensities[right]) {
                    histogram[bin] += 1;
                    count += 1;
                }
                right += 1;
            }
            let noise = if count < self.min_required_elements {
                sparse += 1;
                self.noise_for_empty_window
            } else {
                let mut rightmost = self.bin_count;
                let mut mean = 0.0;
                for _ in 0..3 {
                    budget.charge(
                        rightmost
                            .checked_mul(2)
                            .ok_or_else(|| bad("mean-noise work count overflow"))?,
                    )?;
                    mean = 0.0;
                    for bin in 0..rightmost {
                        mean =
                            finite(mean + histogram[bin] as f64 / count as f64 * bin_value[bin])?;
                    }
                    let mut variance = 0.0;
                    for bin in 0..rightmost {
                        let difference = bin_value[bin] - mean;
                        variance = finite(
                            variance
                                + histogram[bin] as f64 / count as f64 * difference * difference,
                        )?;
                    }
                    // The denominator deliberately remains the original count.
                    let estimate = finite(
                        (mean + variance.sqrt() * self.stdev_multiplier - 1.0) / bin_size + 1.0,
                    )?;
                    if estimate < 0.0 || estimate >= i32::MAX as f64 + 1.0 {
                        return Err(bad(
                            "mean-noise clipping threshold exceeds source integer range",
                        ));
                    }
                    rightmost = (estimate as usize).min(self.bin_count);
                }
                mean.max(1.0)
            };
            let index = result.noise.len();
            result
                .signal_to_noise
                .push(finite(intensities[index] / noise)?);
            result.noise.push(noise);
        }
        result.sparse_window_percent = sparse as f64 * 100.0 / n as f64;
        Ok(result)
    }

    /// [`estimate`](Self::estimate) over a spectrum's m/z positions and
    /// intensities (source `init(spectrum)`).
    ///
    /// # Errors
    ///
    /// As [`estimate`](Self::estimate), and for an invalid spectrum.
    pub fn estimate_spectrum(&self, input: &MSSpectrum) -> Result<MeanNoiseEstimates> {
        self.preflight(input.len())?;
        input.validate()?;
        self.estimate(
            &input.peaks.iter().map(|p| p.mz).collect::<Vec<_>>(),
            &input
                .peaks
                .iter()
                .map(|p| f64::from(p.intensity))
                .collect::<Vec<_>>(),
        )
    }

    /// [`estimate`](Self::estimate) over a chromatogram's retention times and
    /// intensities.
    ///
    /// Native extension: the source template cannot be instantiated for
    /// `MSChromatogram`, because its `computeSTN_` calls `getMZ` on the peaks
    /// (`SignalToNoiseEstimatorMeanIterative.h:262`, `:279`), which
    /// `ChromatogramPeak` does not have.
    ///
    /// # Errors
    ///
    /// As [`estimate`](Self::estimate), and for an invalid chromatogram.
    pub fn estimate_chromatogram(&self, input: &MSChromatogram) -> Result<MeanNoiseEstimates> {
        self.preflight(input.len())?;
        input.validate()?;
        self.estimate(
            &input.peaks.iter().map(|p| p.rt).collect::<Vec<_>>(),
            &input
                .peaks
                .iter()
                .map(|p| f64::from(p.intensity))
                .collect::<Vec<_>>(),
        )
    }

    fn preflight(&self, n: usize) -> Result<()> {
        self.validate()?;
        if n > self.max_points {
            return Err(bad("mean-noise trace exceeds point limit"));
        }
        if n > self.max_work {
            return Err(bad("mean-noise trace exceeds minimum per-point work"));
        }
        Ok(())
    }

    fn maximum(&self, y: &[f64], budget: &mut Budget) -> Result<f64> {
        match self.histogram_range {
            MeanNoiseHistogramRange::Manual { max_intensity } => Ok(max_intensity),
            MeanNoiseHistogramRange::StandardDeviation { factor } => {
                if y.is_empty() {
                    return Ok(0.0);
                }
                budget.charge(
                    y.len()
                        .checked_mul(2)
                        .ok_or_else(|| bad("mean-noise work count overflow"))?,
                )?;
                let mut mean = 0.0;
                for &v in y {
                    mean = finite(mean + v)?;
                }
                mean /= y.len() as f64;
                let mut variance = 0.0;
                for &v in y {
                    let d = mean - v;
                    variance = finite(variance + d * d)?;
                }
                finite(mean + (variance / y.len() as f64).sqrt() * factor)
            }
            MeanNoiseHistogramRange::LegacyPercentile { percentile } => {
                if y.is_empty() {
                    return Err(bad("legacy noise percentile requires input"));
                }
                budget.charge(
                    y.len()
                        .checked_mul(2)
                        .and_then(|n| n.checked_add(100))
                        .ok_or_else(|| bad("mean-noise work count overflow"))?,
                )?;
                // Source reversed max_element comparator actually selects the minimum.
                let minimum = y.iter().copied().fold(f64::INFINITY, f64::min) as f32;
                let size = f64::from(minimum / 100.0);
                if !size.is_finite() || size <= 0.0 {
                    return Err(bad("legacy noise percentile has invalid histogram spacing"));
                }
                let mut histogram = [0usize; 100];
                for &value in y {
                    let value = value as f32;
                    let bin = (f64::from(value - 1.0) / size).trunc();
                    if !bin.is_finite() || !(0.0..100.0).contains(&bin) {
                        return Err(bad(
                            "legacy noise percentile would index outside its histogram",
                        ));
                    }
                    histogram[bin as usize] += 1;
                }
                let required = (f64::from(percentile) * y.len() as f64 / 100.0) as usize;
                let (mut seen, mut cursor) = (0, 0);
                // The source also stops after visiting as many bins as input points.
                while cursor < y.len() && seen < required {
                    budget.charge(1)?;
                    let count = histogram
                        .get(cursor)
                        .ok_or_else(|| bad("legacy noise percentile exhausted its histogram"))?;
                    seen += count;
                    cursor += 1;
                }
                finite((cursor as f64 - 0.5) * size)
            }
        }
    }
}

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("mean-noise numerical expression is nonfinite"))
    }
}
struct Budget {
    remaining: usize,
}
impl Budget {
    fn charge(&mut self, count: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(count)
            .ok_or_else(|| bad("mean-noise estimation exceeds work limit"))?;
        Ok(())
    }
}
