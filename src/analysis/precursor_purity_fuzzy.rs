// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Fuzzy single-scan and RT-interpolated purity from PrecursorPurity.cpp:22–244
//! at OpenMS4-core revision 7c029e8, with checked physical iterator boundaries.

use super::{
    MAX_PURITY_ISOTOPES, MAX_PURITY_PRECURSORS, PrecursorPurity, WorkBudget, finite, invalid,
    validate_precursor, validate_spectrum,
};
use crate::Result;
use crate::kernel::{MSExperiment, MSSpectrum, Peak1D, nearest};
use std::ops::Range;

const NEUTRON_MASS_U: f64 = 1.008_664_915_66;

struct Window {
    seed: usize,
    strict_lower: f64,
    strict_upper: f64,
    fuzzy_lower: f64,
    fuzzy_upper: f64,
    spacing: f64,
    left: Range<usize>,
    right: Range<usize>,
    steps: usize,
}

impl PrecursorPurity {
    /// Source fuzzy-window purity for every precursor of the selected spectrum.
    ///
    /// The initial peak is the globally nearest parent peak, even outside the
    /// isolation window. Isotopes use neutron mass divided by charge, strict ppm
    /// comparisons and half intensity at/outside the strict isolation borders.
    /// Sums and the final division use f32, then the result is widened to f64.
    /// An empty parent returns ones. A zero-width window stops processing and
    /// leaves that precursor and all following precursors at one.
    ///
    /// Source lookup compares the first peak at/above an expected mass with its
    /// physical successor, not its predecessor. This behavior is retained even
    /// outside the logical search range. At physical end, a lone remaining peak
    /// is used; an absent candidate means no match. A matched peak outside the
    /// fuzzy window can therefore yield a source ratio above one. Negative
    /// charge, invalid coordinates/intensities, zero total intensity, nonfinite arithmetic or a
    /// stalled isotope traversal return errors. Charge zero means one. Only the
    /// selected spectra are validated; their MS levels need not have a parent/
    /// child relationship. Work is conservatively preflighted using the shared
    /// purity limits before allocating output or traversing isotope ladders.
    pub fn compute_single_scan(
        experiment: &MSExperiment,
        ms2_index: usize,
        parent_index: usize,
        max_deviation_ppm: f64,
    ) -> Result<Vec<f64>> {
        let mut work = WorkBudget::new();
        Self::compute_single_scan_with_budget(
            experiment,
            ms2_index,
            parent_index,
            max_deviation_ppm,
            &mut work,
        )
    }

    /// Interpolate fuzzy purity by absolute RT distances, retaining source
    /// extrapolation when the selected RT is outside the parent RT interval.
    /// None, an out-of-range/non-MS1 next parent, or a zero/nonfinite parent RT
    /// difference returns the earlier scan result. Other numerical failures are
    /// errors. Both scan calculations share one work budget; results are not
    /// clamped to [0, 1].
    pub fn compute_interpolated(
        experiment: &MSExperiment,
        ms2_index: usize,
        parent_index: usize,
        next_parent: Option<usize>,
        max_deviation_ppm: f64,
    ) -> Result<Vec<f64>> {
        let mut work = WorkBudget::new();
        let early = Self::compute_single_scan_with_budget(
            experiment,
            ms2_index,
            parent_index,
            max_deviation_ppm,
            &mut work,
        )?;
        let Some(next_index) = next_parent.filter(|&i| i < experiment.spectra.len()) else {
            return Ok(early);
        };
        let next = &experiment.spectra[next_index];
        if next.ms_level != 1 {
            return Ok(early);
        }
        let denominator = (next.rt - experiment.spectra[parent_index].rt).abs();
        if !denominator.is_finite() || denominator <= 0.0 {
            return Ok(early);
        }
        let late = Self::compute_single_scan_with_budget(
            experiment,
            ms2_index,
            next_index,
            max_deviation_ppm,
            &mut work,
        )?;
        work.consume(early.len())?;
        let distance = finite(
            (experiment.spectra[ms2_index].rt - experiment.spectra[parent_index].rt).abs(),
            "purity interpolation RT distance",
        )?;
        early
            .into_iter()
            .zip(late)
            .map(|(early, late)| {
                finite(
                    distance * ((late - early) / denominator) + early,
                    "interpolated precursor purity",
                )
            })
            .collect()
    }

    pub(super) fn compute_single_scan_with_budget(
        experiment: &MSExperiment,
        ms2_index: usize,
        parent_index: usize,
        max_deviation_ppm: f64,
        work: &mut WorkBudget,
    ) -> Result<Vec<f64>> {
        let ms2 = experiment
            .spectra
            .get(ms2_index)
            .ok_or_else(|| invalid("purity spectrum index is out of range"))?;
        let parent = experiment
            .spectra
            .get(parent_index)
            .ok_or_else(|| invalid("purity parent index is out of range"))?;
        finite(max_deviation_ppm, "precursor isotope deviation")?;
        if max_deviation_ppm < 0.0 {
            return Err(invalid("precursor isotope deviation must be nonnegative"));
        }
        if ms2.precursors.len() > MAX_PURITY_PRECURSORS {
            return Err(invalid("fuzzy purity exceeds precursor limit"));
        }
        validate_spectrum(ms2, work)?;
        validate_spectrum(parent, work)?;
        work.consume(ms2.precursors.len())?;
        for precursor in &ms2.precursors {
            validate_precursor(precursor)?;
            if precursor.charge < 0 {
                return Err(invalid(
                    "fuzzy precursor purity requires nonnegative charge",
                ));
            }
        }
        if parent.is_empty() {
            return Ok(vec![1.0; ms2.precursors.len()]);
        }
        let windows = plan_windows(ms2, parent, max_deviation_ppm, work)?;
        let mut purities = vec![1.0; ms2.precursors.len()];
        for (output, window) in purities.iter_mut().zip(windows) {
            *output = compute_window(&parent.peaks, &window, max_deviation_ppm)?;
        }
        Ok(purities)
    }
}

fn plan_windows(
    ms2: &MSSpectrum,
    parent: &MSSpectrum,
    deviation: f64,
    work: &mut WorkBudget,
) -> Result<Vec<Window>> {
    let peaks = &parent.peaks;
    let search_work = (usize::BITS - peaks.len().leading_zeros()) as usize + 4;
    let mut windows = Vec::new();
    let mut isotope_steps = 0_usize;
    for precursor in &ms2.precursors {
        let strict_lower = precursor.mz - precursor.isolation_window_lower_offset;
        let strict_upper = finite(
            precursor.mz + precursor.isolation_window_upper_offset,
            "upper precursor isolation boundary",
        )?;
        if strict_lower < 0.0 {
            return Err(invalid(
                "fuzzy purity requires nonnegative isolation boundaries",
            ));
        }
        if strict_lower == strict_upper {
            break;
        }
        let fraction = deviation / 1e6;
        let fuzzy_lower = finite(strict_lower * (1.0 - fraction), "lower fuzzy boundary")?;
        let fuzzy_upper = finite(strict_upper * (1.0 + fraction), "upper fuzzy boundary")?;
        if fuzzy_lower < 0.0 {
            return Err(invalid(
                "fuzzy purity requires a nonnegative lower fuzzy boundary",
            ));
        }
        work.consume(search_work * 5)?;
        let seed = nearest(peaks, precursor.mz, |peak| peak.mz).expect("nonempty validated parent");
        let left = peaks.partition_point(|peak| peak.mz < fuzzy_lower)
            ..peaks.partition_point(|peak| peak.mz <= precursor.mz);
        let right = peaks.partition_point(|peak| peak.mz < precursor.mz)
            ..peaks.partition_point(|peak| peak.mz <= fuzzy_upper);
        let spacing = NEUTRON_MASS_U / f64::from(precursor.charge.max(1));
        let left_active = peaks[seed].mz - spacing > fuzzy_lower;
        let right_active =
            finite(peaks[seed].mz + spacing, "right isotope position")? < fuzzy_upper;
        let steps = if left_active || right_active {
            let span_steps = ((peaks[seed].mz - fuzzy_lower).max(0.0) / spacing).ceil()
                + ((fuzzy_upper - peaks[seed].mz).max(0.0) / spacing).ceil();
            // Every unmatched iteration advances by spacing. Successful strict
            // advances visit distinct peaks; allow two physical successors past
            // the logical range plus margin for the two initial positions.
            let bound = span_steps + (right.end - left.start) as f64 + 4.0;
            if !bound.is_finite() || bound > MAX_PURITY_ISOTOPES as f64 {
                return Err(invalid("fuzzy purity exceeds isotope traversal limit"));
            }
            bound as usize
        } else {
            0
        };
        isotope_steps = isotope_steps
            .checked_add(steps)
            .filter(|&steps| steps <= MAX_PURITY_ISOTOPES)
            .ok_or_else(|| invalid("fuzzy purity exceeds total isotope traversal limit"))?;
        let estimate = steps
            .checked_mul(search_work + 8)
            .and_then(|steps| steps.checked_add(peaks.len()))
            .ok_or_else(|| invalid("fuzzy purity work estimate overflow"))?;
        work.consume(estimate)?;
        windows.push(Window {
            seed,
            strict_lower,
            strict_upper,
            fuzzy_lower,
            fuzzy_upper,
            spacing,
            left,
            right,
            steps,
        });
    }
    Ok(windows)
}

fn compute_window(peaks: &[Peak1D], window: &Window, deviation: f64) -> Result<f64> {
    let seed = peaks[window.seed];
    let mut target = seed.intensity;
    let mut remaining = window.steps;
    for left in [true, false] {
        let range = if left { &window.left } else { &window.right };
        let step = if left {
            -window.spacing
        } else {
            window.spacing
        };
        let mut expected = finite(seed.mz + step, "expected isotope position")?;
        while if left {
            expected > window.fuzzy_lower
        } else {
            expected < window.fuzzy_upper
        } {
            remaining = remaining
                .checked_sub(1)
                .ok_or_else(|| invalid("fuzzy isotope traversal exhausted its bound"))?;
            let candidate = source_candidate(peaks, range, expected);
            let matched = if let Some(peak) = candidate {
                let ppm = finite(
                    (peak.mz - expected).abs() * 1_000_000.0 / expected,
                    "isotope ppm deviation",
                )?;
                (ppm < deviation).then_some(peak)
            } else {
                None
            };
            let next = if let Some(peak) = matched {
                let half = if left {
                    peak.mz <= window.strict_lower
                } else {
                    peak.mz >= window.strict_upper
                };
                target = add_intensity(target, peak.intensity, half)?;
                finite(peak.mz + step, "matched isotope progression")?
            } else {
                finite(expected + step, "isotope progression")?
            };
            if (left && next >= expected) || (!left && next <= expected) {
                return Err(invalid(
                    "fuzzy isotope search cannot make directional progress",
                ));
            }
            expected = next;
        }
    }
    let mut total = seed.intensity;
    for peak in peaks[..window.seed].iter().rev() {
        if peak.mz <= window.fuzzy_lower {
            break;
        }
        total = add_intensity(total, peak.intensity, peak.mz <= window.strict_lower)?;
    }
    for peak in &peaks[window.seed + 1..] {
        if peak.mz >= window.fuzzy_upper {
            break;
        }
        total = add_intensity(total, peak.intensity, peak.mz >= window.strict_upper)?;
    }
    if total <= 0.0 {
        return Err(invalid("fuzzy purity requires positive total intensity"));
    }
    finite(f64::from(target / total), "fuzzy precursor purity")
}

fn source_candidate<'a>(
    peaks: &'a [Peak1D],
    range: &Range<usize>,
    expected: f64,
) -> Option<&'a Peak1D> {
    let index = range.start + peaks[range.clone()].partition_point(|peak| peak.mz < expected);
    let first = peaks.get(index)?;
    let Some(second) = peaks.get(index + 1) else {
        return Some(first);
    };
    Some(
        if (first.mz - expected).abs() < (second.mz - expected).abs() {
            first
        } else {
            second
        },
    )
}

fn add_intensity(total: f32, intensity: f32, half: bool) -> Result<f32> {
    let next = if half {
        (f64::from(total) + 0.5 * f64::from(intensity)) as f32
    } else {
        total + intensity
    };
    finite(f64::from(next), "fuzzy intensity sum")?;
    Ok(next)
}
