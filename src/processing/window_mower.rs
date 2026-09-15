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

/// How a window's start advances, as source `movetype` names it.
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
///
/// Both resource ceilings apply **per spectrum**, and the work ceiling is
/// derived from that spectrum's point count; see [`WindowMower::max_points`]
/// and [`WindowMower::work_per_point`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WindowMower {
    /// Strictly positive finite width in Th. The right endpoint is excluded.
    pub window_size: f64,
    /// Number of highest-intensity peaks selected per full window; zero is valid.
    pub peak_count: usize,
    /// How each window's start advances; see [`WindowMowerMethod`].
    pub method: WindowMowerMethod,
    /// Maximum input peaks in **one** spectrum.
    ///
    /// Source `WindowMower` has no ceiling at all. This one is a per-record
    /// bound, not a per-run one: a run-wide bound shrinks as the run grows and
    /// so rejects real data, while every spectrum of a real run sits far inside
    /// a per-record bound. The benchmark's 1.2 GB LTQ Orbitrap Velos run holds
    /// 88,434,492 peaks in 43,745 spectra, 88 times a run-wide million, while
    /// its largest single spectrum holds 16,766 — 60 times inside the default.
    /// See `docs/WINDOW_MOWER_SUPPORT.md`.
    pub max_points: usize,
    /// Work units available for a spectrum before its point count is credited,
    /// checked before sorting/copying windows. Includes scans, linear selection
    /// and logarithmic sorting/search allowances; this is not a wall-clock or
    /// exact CPU-instruction limit.
    pub max_work: usize,
    /// Work units credited for each input peak of the spectrum being filtered,
    /// so the ceiling for `n` points is `max_work + work_per_point * n`.
    ///
    /// This is the size-derived allowance pattern of the mzML reader
    /// (`src/format/mzml_scaling.rs`), with the spectrum's point count in place
    /// of consumed input bytes. Sliding cost is `Θ(n · w)` for mean window
    /// occupancy `w`, and `w` is set by the data, not by `n`, so a ceiling
    /// linear in `n` accepts every spectrum at a realistic occupancy however
    /// many points it has, and stops one that is quadratic in its own size
    /// after work linear in that size rather than letting it run to completion.
    /// The default is eight times the largest ratio measured over the benchmark
    /// inputs (3,493 work units per point, on the 600-spectrum Q Exactive
    /// profile slice), rounded up to a power of two. See
    /// `docs/WINDOW_MOWER_SUPPORT.md`.
    pub work_per_point: usize,
}

impl Default for WindowMower {
    fn default() -> Self {
        Self {
            window_size: 50.0,
            peak_count: 2,
            method: WindowMowerMethod::Sliding,
            max_points: 1_000_000,
            max_work: 50_000_000,
            work_per_point: 32_768,
        }
    }
}

impl WindowMower {
    /// The work ceiling for a spectrum of `count` points, saturating at
    /// `usize::MAX`.
    fn work_budget(&self, count: usize) -> usize {
        self.max_work
            .saturating_add(self.work_per_point.saturating_mul(count))
    }

    fn validate(&self, count: usize) -> Result<()> {
        if !self.window_size.is_finite() || self.window_size <= 0.0 {
            return Err(bad("window width must be finite and positive"));
        }
        if self.max_points == 0 || self.work_budget(count) == 0 || count > self.max_points {
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
        self.indices(input, &mut Work::new(self.work_budget(input.len())))
    }

    /// Produce an owned filtered spectrum without modifying the input.
    pub fn filtered_spectrum(&self, input: &MSSpectrum) -> Result<MSSpectrum> {
        let indices = self.retained_indices(input)?;
        super::AcquisitionCopies::default().spectrum(input)?;
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

    /// Filter every spectrum in order (source `filterPeakMap`).
    ///
    /// Every ceiling is metered **per spectrum**, not once for the whole
    /// experiment: a run-wide ledger shrinks as the run grows, so it rejects
    /// real data while each of its records stays far inside the same allowance.
    /// The benchmark's 1.2 GB Velos run holds 88,434,492 peaks, 88 times a
    /// run-wide million, while its largest single spectrum needs 16,766 points
    /// and 45,766,084 work units. The metadata-copy ledger has the same shape
    /// and had to move with it, as [`super::baseline::MorphologicalFilter`]
    /// already did: a build that differs only in running that ledger once over
    /// the whole experiment still fails on the Velos run and on the
    /// 40,856-spectrum UK222 run, both with "data array description resource
    /// limit exceeded", while every one of their spectra is far inside the same
    /// allowance on its own. The source has no
    /// ceiling here at all, and what the port bounds is work on data this
    /// experiment already holds in memory, so the record is the level the
    /// bound belongs at.
    ///
    /// # Errors
    ///
    /// Returns the first spectrum's error from [`Self::retained_indices`], or
    /// [`Error::InvalidValue`] when one spectrum's metadata exceeds the
    /// processing copy budget. The experiment is unchanged on error.
    fn filter_experiment(&self, input: &mut MSExperiment) -> Result<()> {
        // Plan all selections before cloning records or committing any spectrum.
        let plans = input
            .spectra
            .iter()
            .map(|spectrum| self.retained_indices(spectrum))
            .collect::<Result<Vec<_>>>()?;
        for spectrum in &input.spectra {
            super::AcquisitionCopies::default().spectrum(spectrum)?;
        }
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
