// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Source sliding/jumping window intensity selection, preserving aligned arrays.
//! See `docs/WINDOW_MOWER_SUPPORT.md` for duplicate membership, final-window
//! quotas and the explicit computational-work accounting convention.

use super::SpectrumFilter;
use crate::{Error, MSExperiment, MSSpectrum, Result};
use std::cmp::Ordering;
use std::str::FromStr;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum WindowMowerMethod {
    /// Slide by one raw peak; stop after the first window reaching the end.
    #[default]
    Sliding,
    /// Start each nonoverlapping window at the next observed peak.
    Jumping,
}

impl FromStr for WindowMowerMethod {
    type Err = Error;
    fn from_str(value: &str) -> Result<Self> {
        match value {
            "slide" => Ok(Self::Sliding),
            "jump" => Ok(Self::Jumping),
            _ => Err(bad("window movement must be slide or jump")),
        }
    }
}

/// Retain the strongest observed peaks in source-compatible m/z windows.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WindowMower {
    /// Strictly positive finite width in Th. The right endpoint is excluded.
    pub window_size: f64,
    /// Number of highest-intensity peaks selected per full window; zero is valid.
    pub peak_count: usize,
    pub method: WindowMowerMethod,
    /// Maximum input peaks across all processed spectra in one call.
    pub max_points: usize,
    /// Work units across all processed spectra, checked before sorting/copying
    /// windows. Includes scans, linear selection and logarithmic sorting/search
    /// allowances; this is not a wall-clock or exact CPU-instruction limit.
    pub max_work: usize,
}

impl Default for WindowMower {
    fn default() -> Self {
        Self {
            window_size: 50.0,
            peak_count: 2,
            method: WindowMowerMethod::Sliding,
            max_points: 1_000_000,
            max_work: 50_000_000,
        }
    }
}

impl WindowMower {
    fn validate(&self, count: usize) -> Result<()> {
        if !self.window_size.is_finite() || self.window_size <= 0.0 {
            return Err(bad("window width must be finite and positive"));
        }
        if self.max_points == 0 || self.max_work == 0 || count > self.max_points {
            return Err(bad(
                "invalid window mower limits or input point limit exceeded",
            ));
        }
        Ok(())
    }

    /// Original input indices in output order. Sliding retains original order;
    /// jumping returns stable ascending m/z order. Coincident positions can
    /// retain more peaks than the requested quota; see the method documentation.
    pub fn retained_indices(&self, input: &MSSpectrum) -> Result<Vec<usize>> {
        self.validate(input.len())?;
        self.indices(input, &mut Work::new(self.max_work))
    }

    /// Produce an owned filtered spectrum without modifying the input.
    pub fn filtered_spectrum(&self, input: &MSSpectrum) -> Result<MSSpectrum> {
        let indices = self.retained_indices(input)?;
        let mut output = input.clone();
        output.select(&indices)?;
        Ok(output)
    }

    /// Apply sliding selection regardless of the configured movement method.
    pub fn filter_sliding(&self, input: &mut MSSpectrum) -> Result<()> {
        Self {
            method: WindowMowerMethod::Sliding,
            ..*self
        }
        .filter_spectrum(input)
    }

    /// Apply jumping selection regardless of the configured movement method.
    pub fn filter_jumping(&self, input: &mut MSSpectrum) -> Result<()> {
        Self {
            method: WindowMowerMethod::Jumping,
            ..*self
        }
        .filter_spectrum(input)
    }

    fn indices(&self, input: &MSSpectrum, work: &mut Work) -> Result<Vec<usize>> {
        let n = input.len();
        work.charge(n.checked_add(1).ok_or_else(|| bad("work size overflow"))?)?;
        input.validate()?;
        if n == 0 || self.peak_count == 0 {
            return Ok(Vec::new());
        }
        let minimum = input
            .peaks
            .iter()
            .map(|p| p.mz)
            .fold(f64::INFINITY, f64::min);
        let maximum = input
            .peaks
            .iter()
            .map(|p| p.mz)
            .fold(f64::NEG_INFINITY, f64::max);
        if !(maximum - minimum).is_finite() {
            return Err(bad("window mower coordinate span overflows"));
        }
        work.sort(n)?;
        let mut order: Vec<usize> = (0..n).collect();
        // Original index makes the ordering deterministic, including signed zeros.
        order.sort_unstable_by(|&a, &b| {
            input.peaks[a]
                .mz
                .partial_cmp(&input.peaks[b].mz)
                .unwrap()
                .then(a.cmp(&b))
        });
        let mut retained = vec![false; n];
        let mut window = Vec::new();
        let mut start = 0;
        let mut end = 0;
        while start < n {
            while end < n {
                work.charge(1)?;
                if input.peaks[order[end]].mz - input.peaks[order[start]].mz >= self.window_size {
                    break;
                }
                end += 1;
            }
            let count = if self.method == WindowMowerMethod::Jumping && end == n {
                let span = input.peaks[order[end - 1]].mz - input.peaks[order[start]].mz;
                let quota = (span / self.window_size * self.peak_count as f64).round();
                if !quota.is_finite() {
                    return Err(bad("final window peak quota is not finite"));
                }
                // Clamp to available points before converting, including huge quotas.
                quota.min((end - start) as f64) as usize
            } else {
                self.peak_count.min(end - start)
            };
            if count > 0 {
                let length = end - start;
                // Linear-time order-statistic selection replaces sorting every window.
                work.charge(
                    length
                        .checked_mul(4)
                        .ok_or_else(|| bad("window work overflow"))?,
                )?;
                window.clear();
                window.extend(start..end);
                if count < length {
                    window.select_nth_unstable_by(count, |&a, &b| {
                        input.peaks[order[b]]
                            .intensity
                            .partial_cmp(&input.peaks[order[a]].intensity)
                            .unwrap()
                            .then(a.cmp(&b))
                    });
                }
                for &slot in &window[..count] {
                    retained[slot] = true;
                }
            }
            // Source sliding stops immediately, not after processing shorter tail windows.
            if end == n {
                break;
            }
            start = match self.method {
                WindowMowerMethod::Sliding => start + 1,
                WindowMowerMethod::Jumping => end,
            };
        }
        match self.method {
            WindowMowerMethod::Sliding => {
                work.charge(
                    n.checked_mul(2)
                        .ok_or_else(|| bad("membership work overflow"))?,
                )?;
                let mut keep_original = vec![false; n];
                let mut begin = 0;
                while begin < n {
                    let mut end = begin + 1;
                    while end < n && input.peaks[order[begin]].mz == input.peaks[order[end]].mz {
                        end += 1;
                    }
                    if retained[begin..end].iter().any(|&v| v) {
                        for &original in &order[begin..end] {
                            keep_original[original] = true;
                        }
                    }
                    begin = end;
                }
                Ok((0..n).filter(|&i| keep_original[i]).collect())
            }
            WindowMowerMethod::Jumping => {
                let count = retained.iter().filter(|&&v| v).count();
                work.sort(count)?;
                work.search(n, count)?;
                let mut keys: Vec<(f64, f32)> = order
                    .iter()
                    .zip(retained)
                    .filter_map(|(&i, keep)| {
                        keep.then_some((input.peaks[i].mz, input.peaks[i].intensity))
                    })
                    .collect();
                keys.sort_unstable_by(compare_key);
                keys.dedup();
                // Peak1D equality includes both m/z and intensity, unlike sliding membership.
                Ok(order
                    .into_iter()
                    .filter(|&i| {
                        keys.binary_search_by(|key| {
                            compare_key(key, &(input.peaks[i].mz, input.peaks[i].intensity))
                        })
                        .is_ok()
                    })
                    .collect())
            }
        }
    }
}

impl SpectrumFilter for WindowMower {
    fn filter_spectrum(&self, input: &mut MSSpectrum) -> Result<()> {
        let indices = self.retained_indices(input)?;
        input.select(&indices)
    }

    fn filter_experiment(&self, input: &mut MSExperiment) -> Result<()> {
        let total = input.spectra.iter().try_fold(0usize, |count, spectrum| {
            count
                .checked_add(spectrum.len())
                .ok_or_else(|| bad("experiment point count overflows"))
        })?;
        self.validate(total)?;
        let mut work = Work::new(self.max_work);
        // Plan all selections before cloning records or committing any spectrum.
        let plans = input
            .spectra
            .iter()
            .map(|spectrum| self.indices(spectrum, &mut work))
            .collect::<Result<Vec<_>>>()?;
        let mut spectra = input.spectra.clone();
        for (spectrum, indices) in spectra.iter_mut().zip(plans) {
            spectrum.select(&indices)?;
        }
        input.spectra = spectra;
        Ok(())
    }
}

fn compare_key(a: &(f64, f32), b: &(f64, f32)) -> Ordering {
    a.0.partial_cmp(&b.0)
        .unwrap()
        .then_with(|| a.1.partial_cmp(&b.1).unwrap())
}

struct Work {
    remaining: usize,
}
impl Work {
    fn new(remaining: usize) -> Self {
        Self { remaining }
    }
    fn charge(&mut self, count: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(count)
            .ok_or_else(|| bad("window mower work limit exceeded"))?;
        Ok(())
    }
    fn sort(&mut self, length: usize) -> Result<()> {
        if length < 2 {
            return Ok(());
        }
        let levels = usize::BITS as usize - (length - 1).leading_zeros() as usize;
        self.charge(
            length
                .checked_mul(levels)
                .and_then(|v| v.checked_mul(4))
                .ok_or_else(|| bad("sort work overflow"))?,
        )
    }
    fn search(&mut self, length: usize, keys: usize) -> Result<()> {
        let levels = if keys == 0 {
            1
        } else {
            usize::BITS as usize - keys.leading_zeros() as usize + 1
        };
        self.charge(
            length
                .checked_mul(levels)
                .ok_or_else(|| bad("search work overflow"))?,
        )
    }
}
fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
