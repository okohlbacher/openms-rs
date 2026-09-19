// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! OpenSWATH chromatogram picking with smoothed spline seeds, source boundary
//! extension and inclusive raw intensity sums. See CHROMATOGRAM_PICKING_SUPPORT.md.
//!
//! [`PeakPickerChromatogram::compatibility`] selects the source behaviours the
//! native default refuses; the internal noise estimate always reproduces the
//! source, because the signal it reads is this picker's own smoothed trace or a
//! chromatogram the picker has already validated.

use super::checked_intensity;
use super::peak_picking::{
    FwhmUnit, PeakPickerHiRes, PickedChromatogram, PickingCompatibility,
    SignalToNoiseEstimatorMedian,
};
use super::smoothing::{GaussFilter, GaussianWidth, SavitzkyGolayFilter};
use crate::kernel::{ChromatogramPeak, DataArray, MSChromatogram};
use crate::{Error, Result};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ChromatogramPickingMethod {
    /// Find boundaries on the original, unsmoothed signal.
    Legacy,
    /// Find boundaries on the smoothed signal; integrate original intensities.
    #[default]
    Corrected,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ChromatogramSmoothing {
    /// Width in seconds, equal to eight Gaussian standard deviations.
    Gaussian { width: f64 },
    /// Fits sample indices; even frame lengths are incremented, as in OpenMS.
    SavitzkyGolay {
        frame_length: usize,
        polynomial_order: usize,
    },
}
impl Default for ChromatogramSmoothing {
    fn default() -> Self {
        Self::Gaussian { width: 50.0 }
    }
}

/// Indices into the original input. Both boundary indices are inclusive.
/// Access original peaks to recover their full f64 RTs; the source metadata
/// arrays leftWidth/rightWidth store rounded f32 RT values.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ChromatogramPeakRegion {
    pub apex_index: usize,
    pub left_index: usize,
    pub right_index: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ChromatogramPickingResult {
    /// Spline seeds and five source arrays. `boundaries` are the HiRes seed
    /// support, while `regions` below holds the final integration boundaries.
    /// `omitted_arrays` explicitly lists unaggregated input profile arrays.
    pub picked: PickedChromatogram,
    /// Full original sampling and annotations, with smoothed intensities.
    pub smoothed: MSChromatogram,
    pub regions: Vec<ChromatogramPeakRegion>,
}

/// Native legacy/corrected PeakPickerChromatogram methods. Crawdad requires an
/// external backend and is not represented by either native method.
#[derive(Clone, Debug)]
pub struct PeakPickerChromatogram {
    pub method: ChromatogramPickingMethod,
    pub smoothing: ChromatogramSmoothing,
    /// Extend through nondecreasing intensities inside this RT distance from
    /// the spline apex. The S/N check still applies. None disables this rule.
    pub peak_width: Option<f64>,
    /// Boundary-extension S/N threshold. Zero disables the boundary S/N gate.
    pub signal_to_noise: f64,
    /// Boundary and apex-report estimator; source window defaults to 1000 s.
    pub noise_estimator: SignalToNoiseEstimatorMedian,
    /// Seed picking is independent: C++ configures HiRes once at S/N 1, with
    /// its 200 s default noise window, even after boundary settings change.
    pub seed_signal_to_noise: f64,
    pub seed_noise_estimator: SignalToNoiseEstimatorMedian,
    /// Report apex S/N even when the boundary S/N gate is disabled.
    pub report_sn: bool,
    pub remove_overlapping_peaks: bool,
    /// Which source behaviours this picker adopts where the native default
    /// refuses; see [`PickingCompatibility`].
    ///
    /// Two flags change what this picker accepts.
    /// [`allow_duplicate_positions`](PickingCompatibility::allow_duplicate_positions)
    /// accepts equal retention times, which `MSChromatogram::isSorted` lets
    /// through to `snt_.init`, and
    /// [`allow_negative_intensities`](PickingCompatibility::allow_negative_intensities)
    /// accepts negative sample intensities. The flag is also handed to the
    /// seed [`PeakPickerHiRes`], as the source's `pp_` member reproduces the
    /// source unconditionally.
    ///
    /// [`allow_unsorted_positions`](PickingCompatibility::allow_unsorted_positions)
    /// has no effect here: `pickChromatogram`
    /// (`ANALYSIS/OPENSWATH/PeakPickerChromatogram.cpp:68-72`) throws
    /// `Exception::IllegalArgument` for a chromatogram that is not sorted by
    /// position, so decreasing retention times stay refused in both profiles.
    ///
    /// The internal noise estimate does not consult this field; see
    /// [`PeakPickerChromatogram::pick_chromatogram`].
    pub compatibility: PickingCompatibility,
    pub max_points: usize,
    /// Per-stage work limit: conservative smoothing work, HiRes seed work,
    /// each noise estimate, and boundary/overlap/integration sample visits.
    pub max_work: usize,
}
impl Default for PeakPickerChromatogram {
    fn default() -> Self {
        Self {
            method: Default::default(),
            smoothing: Default::default(),
            peak_width: None,
            signal_to_noise: 1.0,
            noise_estimator: SignalToNoiseEstimatorMedian {
                window_length: 1000.0,
                ..Default::default()
            },
            seed_signal_to_noise: 1.0,
            seed_noise_estimator: Default::default(),
            report_sn: false,
            remove_overlapping_peaks: false,
            compatibility: PickingCompatibility::default(),
            max_points: 1_000_000,
            max_work: 50_000_000,
        }
    }
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn spend(remaining: &mut usize, amount: usize) -> Result<()> {
    *remaining = remaining
        .checked_sub(amount)
        .ok_or_else(|| invalid("chromatogram picking work limit exceeded"))?;
    Ok(())
}

impl PeakPickerChromatogram {
    /// Pick without changing input. Five arrays are returned in source order:
    /// FWHM, IntegratedIntensity, leftWidth, rightWidth, SN. IntegratedIntensity
    /// is the inclusive sum of raw f32 samples promoted to f64, not trapezoidal area.
    ///
    /// The boundary signal is the input under
    /// [`ChromatogramPickingMethod::Legacy`] and the smoothed trace under
    /// [`ChromatogramPickingMethod::Corrected`], matching source
    /// `pickChromatogram_(chromatogram | smoothed_chrom, ...)`. Its noise is
    /// estimated with [`PickingCompatibility::source`] whatever
    /// [`compatibility`](Self::compatibility) says, because `snt_` is a plain
    /// source estimator over a signal that is either picker-generated or
    /// already validated here. A `win_len` of NaN or `+inf`, which the source's
    /// `setMinFloat("win_len", 1.0)` restriction lets through, therefore picks
    /// rather than failing.
    ///
    /// # Errors
    ///
    /// * [`Error::UnsortedData`] when retention times decrease, as source
    ///   `pickChromatogram` throws `Exception::IllegalArgument` for a
    ///   chromatogram that is not sorted by position. No compatibility flag
    ///   lifts this.
    /// * [`Error::InvalidValue`] for invalid options or resource limits,
    ///   non-finite coordinates or intensities, and — unless
    ///   [`compatibility`](Self::compatibility) allows them — negative
    ///   intensities or duplicate retention times. Also for a `win_len` below
    ///   the source's minimum of one, a `bin_count` below three, an integration
    ///   bound outside the input, an intensity that leaves `f32` range, and an
    ///   exhausted work or point budget.
    /// * The errors of the configured smoother and of the seed
    ///   [`PeakPickerHiRes`].
    pub fn pick_chromatogram(&self, input: &MSChromatogram) -> Result<ChromatogramPickingResult> {
        if self.max_points == 0
            || self.max_work == 0
            || input.len() > self.max_points
            || [self.signal_to_noise, self.seed_signal_to_noise]
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0)
            || self.peak_width.is_some_and(|v| !v.is_finite() || v <= 0.0)
        {
            return Err(invalid(
                "invalid chromatogram picker options or point limit",
            ));
        }
        input.validate()?;
        if !self.compatibility.allow_negative_intensities
            && input.peaks.iter().any(|p| p.intensity < 0.0)
        {
            return Err(invalid(
                "chromatogram picking requires nonnegative intensities; negative intensities need PickingCompatibility::allow_negative_intensities",
            ));
        }
        // The source throws `Exception::IllegalArgument` here
        // (`PeakPickerChromatogram.cpp:68-72`), so this refusal is not one
        // `allow_unsorted_positions` may lift.
        if input.peaks.windows(2).any(|p| p[0].rt > p[1].rt) {
            return Err(Error::UnsortedData);
        }
        // `isSorted` accepts equal positions, so the source reaches both the
        // seed picker and `snt_.init` with duplicate retention times.
        if !self.compatibility.allow_duplicate_positions
            && input.peaks.windows(2).any(|p| p[0].rt == p[1].rt)
        {
            return Err(invalid(
                "chromatogram picking requires distinct RT samples; duplicates need PickingCompatibility::allow_duplicate_positions",
            ));
        }
        let mut copies = super::AcquisitionCopies::default();
        copies.chromatogram(input)?;
        let mut smoothed = input.clone();
        self.smooth(&mut smoothed)?;
        let mut seed_noise = self.seed_noise_estimator.clone();
        seed_noise.max_points = seed_noise.max_points.min(self.max_points);
        seed_noise.max_work = seed_noise.max_work.min(self.max_work);
        let mut picked = PeakPickerHiRes {
            signal_to_noise: self.seed_signal_to_noise,
            noise_estimator: seed_noise,
            spacing_difference: 0.0,
            spacing_difference_gap: 0.0,
            report_fwhm: Some(FwhmUnit::Absolute),
            compatibility: self.compatibility,
            max_points: self.max_points,
            max_work: self.max_work,
            ..Default::default()
        }
        .pick_chromatogram_with_acquisition(&smoothed, false, &mut copies)?;
        let boundary_signal = match self.method {
            ChromatogramPickingMethod::Legacy => input,
            ChromatogramPickingMethod::Corrected => &smoothed,
        };
        let noise = if self.signal_to_noise > 0.0 || self.report_sn {
            let mut estimator = self.noise_estimator.clone();
            estimator.max_points = estimator.max_points.min(self.max_points);
            estimator.max_work = estimator.max_work.min(self.max_work);
            // Source `snt_.init(chromatogram)`
            // (`ANALYSIS/OPENSWATH/PeakPickerChromatogram.cpp:171`), where
            // `chromatogram` is this same boundary signal. `snt_` is a plain
            // source estimator with no native variant, and the signal it reads
            // is not caller data: under `corrected` this picker produced it by
            // smoothing, and under `legacy` the caller's chromatogram already
            // passed the checks above. So the estimate reproduces the source
            // unconditionally rather than following `self.compatibility`.
            //
            // This matters on ordinary data: a Savitzky-Golay frame has
            // negative coefficients, so an entirely nonnegative chromatogram
            // smooths to a signal containing negative samples, which the strict
            // profile refused. Passing the peaks also widens each `f32` with
            // `cvtss2sd` semantics instead of `f64::from`.
            Some(
                estimator
                    .estimate_peaks(
                        &boundary_signal.peaks,
                        &PickingCompatibility::source(),
                        None,
                    )?
                    .signal_to_noise,
            )
        } else {
            None
        };
        let mut remaining = self.max_work;
        let (mut regions, apex_sn) = self.extend_boundaries(
            &boundary_signal.peaks,
            &picked.chromatogram.peaks,
            noise.as_deref(),
            &mut remaining,
        )?;
        if self.remove_overlapping_peaks {
            adjust_overlaps(&boundary_signal.peaks, &mut regions, &mut remaining)?;
        }
        let (mut areas, mut lefts, mut rights) = (Vec::new(), Vec::new(), Vec::new());
        for region in &regions {
            if region.left_index > region.right_index || region.right_index >= input.len() {
                return Err(invalid(
                    "chromatogram picker produced invalid integration bounds",
                ));
            }
            let mut area = 0.0;
            for sample in &input.peaks[region.left_index..=region.right_index] {
                spend(&mut remaining, 1)?;
                area += f64::from(sample.intensity);
            }
            areas.push(checked_intensity(area)?);
            lefts.push(checked_intensity(input.peaks[region.left_index].rt)?);
            rights.push(checked_intensity(input.peaks[region.right_index].rt)?);
        }
        // HiRes has already reported every discarded profile array. These
        // output arrays describe peaks and cannot reuse profile sample values.
        picked.chromatogram.float_data_arrays.extend([
            DataArray::new("IntegratedIntensity", areas),
            DataArray::new("leftWidth", lefts),
            DataArray::new("rightWidth", rights),
            DataArray::new("SN", apex_sn),
        ]);
        picked.chromatogram.validate()?;
        Ok(ChromatogramPickingResult {
            picked,
            smoothed,
            regions,
        })
    }

    /// Replace the input only after the complete result succeeds. Use the
    /// returning method above to retain smoothing, bounds and omitted-array reports.
    pub fn filter_chromatogram(&self, input: &mut MSChromatogram) -> Result<()> {
        *input = self.pick_chromatogram(input)?.picked.chromatogram;
        Ok(())
    }

    fn smooth(&self, input: &mut MSChromatogram) -> Result<()> {
        let mut remaining = self.max_work;
        spend(&mut remaining, input.len())?;
        match self.smoothing {
            ChromatogramSmoothing::Gaussian { width } => {
                let smoother = GaussFilter::new(GaussianWidth::Absolute(width))?;
                let coefficients = (4.0 * (width / 8.0) / 0.01).ceil();
                if !coefficients.is_finite()
                    || coefficients >= smoother.algorithm.max_coefficients as f64
                {
                    return Err(invalid(
                        "Gaussian chromatogram kernel exceeds coefficient limit",
                    ));
                }
                spend(&mut remaining, coefficients as usize + 1)?;
                // The fixed 0.01 s table extends by less than 0.02 s beyond
                // width/2. Count a conservative support bound before convolution.
                let support = width / 2.0 + 0.02;
                for peak in &input.peaks {
                    let low = peak.rt - support;
                    let high = peak.rt + support;
                    if !low.is_finite() || !high.is_finite() {
                        return Err(invalid("Gaussian chromatogram support overflows"));
                    }
                    let begin = input.peaks.partition_point(|p| p.rt < low);
                    let end = input.peaks.partition_point(|p| p.rt <= high);
                    spend(&mut remaining, end - begin)?;
                }
                smoother.filter_chromatogram(input)
            }
            ChromatogramSmoothing::SavitzkyGolay {
                frame_length,
                polynomial_order,
            } => {
                let frame = frame_length
                    .checked_add(usize::from(frame_length % 2 == 0))
                    .ok_or_else(|| invalid("Savitzky-Golay frame overflows"))?;
                let setup = frame
                    .checked_mul(frame)
                    .and_then(|n| n.checked_mul(polynomial_order.checked_add(1)?))
                    .ok_or_else(|| invalid("Savitzky-Golay setup work overflows"))?;
                spend(&mut remaining, setup)?;
                if frame <= input.len() {
                    spend(
                        &mut remaining,
                        input
                            .len()
                            .checked_mul(frame)
                            .ok_or_else(|| invalid("Savitzky-Golay convolution work overflows"))?,
                    )?;
                }
                SavitzkyGolayFilter::new(frame_length, polynomial_order)?.filter_chromatogram(input)
            }
        }
    }

    fn extend_boundaries(
        &self,
        samples: &[ChromatogramPeak],
        seeds: &[ChromatogramPeak],
        sn: Option<&[f64]>,
        remaining: &mut usize,
    ) -> Result<(Vec<ChromatogramPeakRegion>, Vec<f32>)> {
        let mut regions = Vec::new();
        let mut apex_sn = Vec::new();
        let mut current = 0;
        for seed in seeds {
            current = closest_peak(samples, seed.rt, current, remaining)?;
            if current == 0 || current >= samples.len().saturating_sub(1) {
                return Err(invalid(
                    "picked apex has no two neighboring chromatogram samples",
                ));
            }
            apex_sn.push(checked_intensity(
                sn.map_or(-1.0, |values| values[current]),
            )?);
            // Source always includes the immediate neighbors. The monotonic,
            // forced-width and S/N checks begin at the second sample away.
            let (mut left, mut right) = (current - 1, current + 1);
            while left > 0 {
                spend(remaining, 1)?;
                let candidate = left - 1;
                let forced = self
                    .peak_width
                    .is_some_and(|width| (samples[candidate].rt - seed.rt).abs() < width);
                if !(samples[candidate].intensity < samples[left].intensity || forced)
                    || (self.signal_to_noise > 0.0
                        && sn.expect("positive threshold initialized noise")[candidate]
                            < self.signal_to_noise)
                {
                    break;
                }
                left = candidate;
            }
            while right + 1 < samples.len() {
                spend(remaining, 1)?;
                let candidate = right + 1;
                let forced = self
                    .peak_width
                    .is_some_and(|width| (samples[candidate].rt - seed.rt).abs() < width);
                if !(samples[candidate].intensity < samples[right].intensity || forced)
                    || (self.signal_to_noise > 0.0
                        && sn.expect("positive threshold initialized noise")[candidate]
                            < self.signal_to_noise)
                {
                    break;
                }
                right = candidate;
            }
            regions.push(ChromatogramPeakRegion {
                apex_index: current,
                left_index: left,
                right_index: right,
            });
        }
        Ok((regions, apex_sn))
    }
}

/// Source lookup intentionally selects the right neighbor on ties and returns
/// len for a target at/beyond the final sample. Callers check this end sentinel.
fn closest_peak(
    samples: &[ChromatogramPeak],
    target: f64,
    mut current: usize,
    remaining: &mut usize,
) -> Result<usize> {
    while current < samples.len() {
        spend(remaining, 1)?;
        if target < samples[current].rt {
            if current > 0
                && (target - samples[current - 1].rt).abs() < (target - samples[current].rt).abs()
            {
                current -= 1;
            }
            return Ok(current);
        }
        current += 1;
    }
    Ok(current)
}

fn adjust_overlaps(
    samples: &[ChromatogramPeak],
    regions: &mut [ChromatogramPeakRegion],
    remaining: &mut usize,
) -> Result<()> {
    for index in 0..regions.len().saturating_sub(1) {
        spend(remaining, 1)?;
        if regions[index].right_index <= regions[index + 1].left_index {
            continue;
        }
        let mut right = regions[index].apex_index;
        while right + 1 < samples.len() {
            spend(remaining, 1)?;
            if samples[right + 1].intensity >= samples[right].intensity {
                break;
            }
            right += 1;
        }
        let mut left = regions[index + 1].apex_index;
        while left > 0 {
            spend(remaining, 1)?;
            if samples[left - 1].intensity >= samples[left].intensity {
                break;
            }
            left -= 1;
        }
        if left < right {
            // These assignments are sequential in C++. The second uses the
            // updated left value and can leave overlap; do not collapse both
            // borders to a common midpoint or silently claim disjoint regions.
            left += (right - left) / 2;
            right = left + (right - left) / 2;
        }
        regions[index].right_index = right;
        regions[index + 1].left_index = left;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_closest_lookup_ties_go_right_and_final_sample_returns_end() {
        let samples = [
            ChromatogramPeak::new(1.0, 1.0),
            ChromatogramPeak::new(3.0, 1.0),
        ];
        let mut work = 100;
        assert_eq!(closest_peak(&samples, 2.0, 0, &mut work).unwrap(), 1);
        assert_eq!(closest_peak(&samples, 1.0, 0, &mut work).unwrap(), 0);
        assert_eq!(closest_peak(&samples, 0.0, 0, &mut work).unwrap(), 0);
        assert_eq!(closest_peak(&samples, 3.0, 0, &mut work).unwrap(), 2);
        assert_eq!(closest_peak(&samples, 4.0, 0, &mut work).unwrap(), 2);
        assert!(closest_peak(&samples, 2.0, 0, &mut 0).is_err());
    }

    #[test]
    fn overlap_midpoints_retain_source_sequential_assignment() {
        // Deliberately supplied seed indices model the legacy mismatch between
        // smoothed apices and a monotonically descending raw signal.
        let samples: Vec<_> = [9., 8., 7., 6., 5., 4., 3., 2., 1.]
            .iter()
            .enumerate()
            .map(|(i, &intensity)| ChromatogramPeak::new(i as f64, intensity))
            .collect();
        let mut regions = [
            ChromatogramPeakRegion {
                apex_index: 2,
                left_index: 0,
                right_index: 8,
            },
            ChromatogramPeakRegion {
                apex_index: 4,
                left_index: 3,
                right_index: 8,
            },
        ];
        adjust_overlaps(&samples, &mut regions, &mut 100).unwrap();
        assert_eq!(regions[1].left_index, 6); // (4+8)/2
        assert_eq!(regions[0].right_index, 7); // (new left6+8)/2
    }
}
