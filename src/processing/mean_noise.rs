// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Chris Bielow, OpenMS Rust contributors $

//! Iteratively clipped histogram-mean signal-to-noise estimation.
//! See docs/MEAN_NOISE_SUPPORT.md for the source's fixed denominator and
//! historical percentile conventions.

use crate::kernel::{MSChromatogram, MSSpectrum};
use crate::{Error, Result};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MeanNoiseHistogramRange {
    StandardDeviation {
        factor: f64,
    },
    Manual {
        max_intensity: f64,
    },
    /// Literal historical source heuristic, including its minimum-as-maximum
    /// and bounded histogram scan. Unsafe source indices return errors.
    /// This mode retains the source's f32 intensity arithmetic.
    LegacyPercentile {
        percentile: u8,
    },
}
impl Default for MeanNoiseHistogramRange {
    fn default() -> Self {
        Self::StandardDeviation { factor: 3.0 }
    }
}

/// Three histogram clipping passes with the original window count as denominator.
#[derive(Clone, Debug, PartialEq)]
pub struct SignalToNoiseEstimatorMeanIterative {
    pub histogram_range: MeanNoiseHistogramRange,
    /// Full coordinate width, at least one; left inclusive and right exclusive.
    pub window_length: f64,
    pub bin_count: usize,
    pub stdev_multiplier: f64,
    pub min_required_elements: usize,
    pub noise_for_empty_window: f64,
    pub max_points: usize,
    pub max_bins: usize,
    /// Output points, global intensity scans, histogram slots/updates and bin visits.
    pub max_work: usize,
}
impl Default for SignalToNoiseEstimatorMeanIterative {
    fn default() -> Self {
        Self {
            histogram_range: Default::default(),
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

#[derive(Clone, Debug, PartialEq)]
pub struct MeanNoiseEstimates {
    pub signal_to_noise: Vec<f64>,
    pub noise: Vec<f64>,
    pub max_intensity: f64,
    pub sparse_window_percent: f64,
}

impl SignalToNoiseEstimatorMeanIterative {
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
