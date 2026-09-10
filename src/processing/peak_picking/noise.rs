// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::{bad, validate_signal};
use crate::Result;

/// Histogram upper range. Values above the range go into its final bin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NoiseHistogramRange {
    /// Mean plus factor times population standard deviation.
    StandardDeviation {
        factor: f64,
    },
    Manual {
        max_intensity: f64,
    },
}
impl Default for NoiseHistogramRange {
    fn default() -> Self {
        Self::StandardDeviation { factor: 3.0 }
    }
}

/// Sliding-window histogram median, with within-bin median interpolation.
#[derive(Clone, Debug)]
pub struct SignalToNoiseEstimatorMedian {
    pub histogram_range: NoiseHistogramRange,
    /// Full window width in coordinate units (Th or seconds), at least one.
    pub window_length: f64,
    pub bin_count: usize,
    pub min_required_elements: usize,
    pub noise_for_empty_window: f64,
    pub max_points: usize,
    pub max_bins: usize,
    /// Histogram insertions, removals, and median-bin visits.
    pub max_work: usize,
}
impl Default for SignalToNoiseEstimatorMedian {
    fn default() -> Self {
        Self {
            histogram_range: Default::default(),
            window_length: 200.0,
            bin_count: 30,
            min_required_elements: 10,
            noise_for_empty_window: 1e20,
            max_points: 1_000_000,
            max_bins: 1_000_000,
            max_work: 50_000_000,
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct NoiseEstimates {
    pub signal_to_noise: Vec<f64>,
    pub noise: Vec<f64>,
    pub max_intensity: f64,
    pub sparse_window_percent: f64,
    pub histogram_rightmost_percent: f64,
}
impl SignalToNoiseEstimatorMedian {
    pub(super) fn validate(&self) -> Result<()> {
        if !self.window_length.is_finite()
            || self.window_length < 1.0
            || self.bin_count < 3
            || self.bin_count > self.max_bins
            || self.min_required_elements == 0
            || !self.noise_for_empty_window.is_finite()
            || self.noise_for_empty_window <= 0.0
            || self.max_points == 0
            || self.max_work == 0
        {
            return Err(bad("invalid median-noise options or resource limits"));
        }
        match self.histogram_range {
            NoiseHistogramRange::StandardDeviation { factor }
                if !factor.is_finite() || !(0.0..=999.0).contains(&factor) =>
            {
                Err(bad("noise standard-deviation factor must be in 0..=999"))
            }
            NoiseHistogramRange::Manual { max_intensity }
                if !max_intensity.is_finite() || max_intensity <= 0.0 =>
            {
                Err(bad(
                    "manual noise histogram maximum must be positive and finite",
                ))
            }
            _ => Ok(()),
        }
    }
    /// Input coordinates must be strictly increasing and intensities nonnegative.
    pub fn estimate(&self, positions: &[f64], intensities: &[f64]) -> Result<NoiseEstimates> {
        self.validate()?;
        validate_signal(positions, intensities, self.max_points)?;
        let n = positions.len();
        let max_intensity = match self.histogram_range {
            NoiseHistogramRange::Manual { max_intensity } => max_intensity,
            NoiseHistogramRange::StandardDeviation { factor } => {
                if n == 0 {
                    0.0
                } else {
                    let mean = intensities.iter().sum::<f64>() / n as f64;
                    let variance = intensities
                        .iter()
                        .map(|&y| (mean - y) * (mean - y))
                        .sum::<f64>()
                        / n as f64;
                    mean + variance.sqrt() * factor
                }
            }
        };
        if !max_intensity.is_finite() {
            return Err(bad("noise histogram maximum overflows"));
        }
        let width = (max_intensity / self.bin_count as f64).max(1.0);
        let to_bin = |y: f64| (y / width).min((self.bin_count - 1) as f64) as usize;
        let mut histogram = vec![0usize; self.bin_count];
        let mut result = NoiseEstimates {
            signal_to_noise: Vec::with_capacity(n),
            noise: Vec::with_capacity(n),
            max_intensity,
            sparse_window_percent: 0.0,
            histogram_rightmost_percent: 0.0,
        };
        let (mut left, mut right, mut work) = (0, 0, 0usize);
        let mut charge = || -> Result<()> {
            if work == self.max_work {
                Err(bad("noise estimation exceeds configured work limit"))
            } else {
                work += 1;
                Ok(())
            }
        };
        for i in 0..n {
            let low = positions[i] - self.window_length / 2.0;
            let high = positions[i] + self.window_length / 2.0;
            if !low.is_finite() || !high.is_finite() {
                return Err(bad("noise window coordinates overflow"));
            }
            while left < right && positions[left] < low {
                charge()?;
                histogram[to_bin(intensities[left])] -= 1;
                left += 1;
            }
            // A disjoint next window can skip samples never inserted in the previous window.
            while right < n && positions[right] < low {
                charge()?;
                right += 1;
                left = right;
            }
            while right < n && positions[right] <= high {
                charge()?;
                histogram[to_bin(intensities[right])] += 1;
                right += 1;
            }
            let count = right - left;
            let noise = if count < self.min_required_elements {
                result.sparse_window_percent += 1.0;
                self.noise_for_empty_window
            } else {
                let rank = count.div_ceil(2);
                let mut before = 0;
                let mut median = 0;
                while median < self.bin_count {
                    charge()?;
                    if before + histogram[median] >= rank {
                        break;
                    }
                    before += histogram[median];
                    median += 1;
                }
                if median == self.bin_count - 1 {
                    result.histogram_rightmost_percent += 1.0;
                }
                (median as f64 * width + (rank - before) as f64 / histogram[median] as f64 * width)
                    .max(1.0)
            };
            let ratio = intensities[i] / noise;
            if !noise.is_finite() || !ratio.is_finite() {
                return Err(bad("noise estimate is not finite"));
            }
            result.noise.push(noise);
            result.signal_to_noise.push(ratio);
        }
        if n != 0 {
            result.sparse_window_percent *= 100.0 / n as f64;
            result.histogram_rightmost_percent *= 100.0 / n as f64;
        }
        Ok(result)
    }
}
