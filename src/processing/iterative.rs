// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Hannes Roest, OpenMS Rust contributors $

//! Iterative integration and recentering of high-resolution profile peaks.
//!
//! Ported from `PROCESSING/CENTROIDING/PeakPickerIterative.h` at revision
//! `7c029e8`. The source's asymmetric rightward recentering and intermediate
//! f32 centroid rounding are retained. See `docs/ITERATIVE_PICKING_SUPPORT.md`.

use super::SpectrumFilter;
use super::peak_picking::{
    PeakBoundary, PeakPickerHiRes, PickedSpectrum, SignalToNoiseEstimatorMedian,
};
use crate::kernel::{DataArray, MSExperiment, MSSpectrum, Peak1D, SpectrumType};
use crate::{Error, Result};

/// Indices into the original profile, aligned with output centroid order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IterativePeakRegion {
    /// Ordinal in the original m/z-ordered HiRes seed list.
    pub seed_index: usize,
    /// First raw point strictly greater than the corresponding seed m/z.
    pub initial_center_index: usize,
    /// Final raw center chosen by the source's restricted nearest-point search.
    pub center_index: usize,
    /// Inclusive integration boundaries. Original peaks retain full f64 m/z.
    pub left_index: usize,
    pub right_index: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct IterativePickingResult {
    /// Full-f64 final integration boundaries and omitted input array names.
    pub picked: PickedSpectrum,
    pub regions: Vec<IterativePeakRegion>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct IterativeExperimentResult {
    pub experiment: MSExperiment,
    /// None means copied without picking because ms1_only excluded the spectrum.
    pub spectrum_regions: Vec<Option<Vec<IterativePeakRegion>>>,
    pub omitted_spectrum_arrays: Vec<Vec<String>>,
}

/// Native PeakPickerIterative options. Widths and windows use input m/z units.
#[derive(Clone, Debug)]
pub struct PeakPickerIterative {
    pub signal_to_noise: f64,
    /// Expected half width; zero disables extension through rising intensities.
    pub peak_width: f64,
    pub spacing_difference: f64,
    /// Refinement noise defaults to a 20-unit window and 30 histogram bins.
    /// HiRes seeds retain their own source-default 200-unit noise window.
    pub noise_estimator: SignalToNoiseEstimatorMedian,
    pub iterations: usize,
    pub check_width_internally: bool,
    /// Applies only to experiment picking; pick_spectrum always processes input.
    pub ms1_only: bool,
    /// Clear the three generated float arrays only in experiment picking.
    /// Exact regions remain available in the result.
    pub clear_meta_data: bool,
    /// Per-spectrum whole-input and seed point bound, checked before allocation.
    pub max_points: usize,
    /// Per-stage work bound for HiRes seed visits, each noise estimate, and
    /// combined seed association/refinement/sorting/overlap-suppression visits.
    pub max_work: usize,
}

impl Default for PeakPickerIterative {
    fn default() -> Self {
        Self {
            signal_to_noise: 1.0,
            peak_width: 0.0,
            spacing_difference: 1.5,
            noise_estimator: SignalToNoiseEstimatorMedian {
                window_length: 20.0,
                ..Default::default()
            },
            iterations: 5,
            check_width_internally: false,
            ms1_only: false,
            clear_meta_data: false,
            max_points: 1_000_000,
            max_work: 50_000_000,
        }
    }
}

#[derive(Clone, Debug)]
struct Candidate {
    region: IterativePeakRegion,
    priority: f32,
    intensity: f64,
    mz: f32,
    valid: bool,
}

impl PeakPickerIterative {
    fn validate(&self) -> Result<()> {
        if [
            self.signal_to_noise,
            self.peak_width,
            self.spacing_difference,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
            || self.iterations == 0
            || self.max_points == 0
            || self.max_work == 0
        {
            return Err(bad("invalid iterative picker options or resource limits"));
        }
        Ok(())
    }

    fn validate_input(&self, input: &MSSpectrum) -> Result<()> {
        self.validate()?;
        if input.len() > self.max_points {
            return Err(bad("iterative picker input exceeds point limit"));
        }
        input.validate()?;
        if input.peaks.iter().any(|p| p.mz < 0.0 || p.intensity < 0.0) {
            return Err(bad(
                "iterative picking requires nonnegative m/z and intensities",
            ));
        }
        if input.peaks.windows(2).any(|w| w[0].mz > w[1].mz) {
            return Err(Error::UnsortedData);
        }
        if input.peaks.windows(2).any(|w| w[0].mz == w[1].mz) {
            return Err(bad(
                "iterative picking requires distinct profile coordinates",
            ));
        }
        Ok(())
    }

    /// Pick a spectrum without mutating input. Output intensity is the inclusive
    /// raw sample sum, not a trapezoidal area. Equal-priority seeds retain order.
    pub fn pick_spectrum(&self, input: &MSSpectrum) -> Result<IterativePickingResult> {
        self.validate_input(input)?;
        if input.len() < 3 {
            return self.output(input, Vec::new());
        }
        let mut seed_noise = SignalToNoiseEstimatorMedian::default();
        if self.signal_to_noise > 0.0 {
            self.limit_noise(&mut seed_noise)?;
        }
        let seed_picker = PeakPickerHiRes {
            signal_to_noise: self.signal_to_noise,
            spacing_difference: self.spacing_difference,
            noise_estimator: seed_noise,
            max_points: self.max_points,
            max_work: self.max_work,
            ..Default::default()
        };
        // Only m/z and intensity seed this algorithm; all profile annotations
        // are reported as omitted, including ion mobility annotations.
        let seed_input = MSSpectrum {
            peaks: input.peaks.clone(),
            ..Default::default()
        };
        let seeds = seed_picker.pick_spectrum(&seed_input)?;
        self.refine_with_seeds(input, &seeds.spectrum)
    }

    /// Pick selected spectra and preserve chromatograms and experiment metadata.
    /// All selected-spectrum errors leave the original experiment unchanged.
    pub fn pick_experiment(&self, input: &MSExperiment) -> Result<IterativeExperimentResult> {
        self.validate()?;
        // Preflight selected records before validating or cloning the experiment.
        for spectrum in &input.spectra {
            if (!self.ms1_only || spectrum.ms_level == 1) && spectrum.len() > self.max_points {
                return Err(bad("iterative picker spectrum exceeds point limit"));
            }
        }
        input.validate()?;
        let mut result = IterativeExperimentResult {
            experiment: input.clone(),
            spectrum_regions: Vec::with_capacity(input.spectra.len()),
            omitted_spectrum_arrays: Vec::with_capacity(input.spectra.len()),
        };
        for (index, spectrum) in input.spectra.iter().enumerate() {
            if self.ms1_only && spectrum.ms_level != 1 {
                result.spectrum_regions.push(None);
                result.omitted_spectrum_arrays.push(Vec::new());
                continue;
            }
            let mut picked = self.pick_spectrum(spectrum)?;
            if self.clear_meta_data {
                picked.picked.spectrum.float_data_arrays.clear();
            }
            result.experiment.spectra[index] = picked.picked.spectrum;
            result.spectrum_regions.push(Some(picked.regions));
            result
                .omitted_spectrum_arrays
                .push(picked.picked.omitted_arrays);
        }
        Ok(result)
    }

    fn limit_noise(&self, noise: &mut SignalToNoiseEstimatorMedian) -> Result<()> {
        // Histogram allocation/initialization also consumes work, even when
        // there are no data points. The estimator bounds subsequent bin visits.
        noise.max_work = noise.max_work.min(
            self.max_work
                .checked_sub(noise.bin_count)
                .filter(|&remaining| remaining > 0)
                .ok_or_else(|| bad("iterative noise histogram exceeds work limit"))?,
        );
        noise.max_points = noise.max_points.min(self.max_points);
        Ok(())
    }

    fn refine_with_seeds(
        &self,
        input: &MSSpectrum,
        seeds: &MSSpectrum,
    ) -> Result<IterativePickingResult> {
        self.validate_input(input)?;
        if seeds.len() > self.max_points {
            return Err(bad("iterative seed count exceeds point limit"));
        }
        seeds.validate()?;
        if seeds.peaks.iter().any(|p| p.mz < 0.0 || p.intensity < 0.0) {
            return Err(bad("iterative seeds must be nonnegative"));
        }
        if seeds.peaks.windows(2).any(|p| p[0].mz > p[1].mz) {
            return Err(Error::UnsortedData);
        }
        if input.len() < 3 {
            return self.output(input, Vec::new());
        }
        let mut work = self.max_work;
        let mut candidates = Vec::new();
        let mut seed_index = 0;
        // Deliberately advance at most one seed per raw sample. This is neither
        // nearest-neighbor lookup nor an independent upper_bound for each seed.
        for (index, peak) in input.peaks.iter().enumerate() {
            if seed_index == seeds.len() {
                break;
            }
            spend(&mut work, 1)?;
            if peak.mz > seeds.peaks[seed_index].mz {
                if index == 0 || index + 1 >= input.len() {
                    return Err(bad(
                        "iterative seed association lacks two neighboring raw points",
                    ));
                }
                candidates.push(Candidate {
                    region: IterativePeakRegion {
                        seed_index,
                        initial_center_index: index,
                        center_index: index,
                        left_index: index - 1,
                        right_index: index + 1,
                    },
                    priority: seeds.peaks[seed_index].intensity,
                    intensity: 0.0,
                    mz: 0.0,
                    valid: true,
                });
                seed_index += 1;
            }
        }
        charge_sort(&mut work, candidates.len())?;
        candidates.sort_by(|a, b| b.priority.total_cmp(&a.priority));
        let noise = if self.signal_to_noise > 0.0 {
            let mut estimator = self.noise_estimator.clone();
            self.limit_noise(&mut estimator)?;
            let xs: Vec<_> = input.peaks.iter().map(|p| p.mz).collect();
            let ys: Vec<_> = input.peaks.iter().map(|p| f64::from(p.intensity)).collect();
            Some(estimator.estimate(&xs, &ys)?.signal_to_noise)
        } else {
            None
        };
        // Prevent a huge iteration count from doing unbounded empty work.
        if !candidates.is_empty() {
            let minimum = self
                .iterations
                .checked_mul(candidates.len())
                .ok_or_else(|| bad("iterative refinement count overflows"))?;
            if minimum > work {
                return Err(bad("iterative refinements exceed work limit"));
            }
            for _ in 0..self.iterations {
                for candidate in &mut candidates {
                    spend(&mut work, 1)?;
                    self.recenter(input, candidate, noise.as_deref(), &mut work)?;
                }
            }
        }
        for i in 0..candidates.len() {
            spend(&mut work, 1)?;
            if !candidates[i].valid {
                continue;
            }
            let left = input.peaks[candidates[i].region.left_index].mz;
            let right = input.peaks[candidates[i].region.right_index].mz;
            for candidate in &mut candidates[i + 1..] {
                spend(&mut work, 1)?;
                let mz = f64::from(candidate.mz);
                if left <= mz && mz <= right {
                    candidate.valid = false;
                }
            }
        }
        candidates.retain(|c| c.valid);
        charge_sort(&mut work, candidates.len())?;
        candidates.sort_by(|a, b| a.mz.total_cmp(&b.mz));
        self.output(input, candidates)
    }

    fn recenter(
        &self,
        input: &MSSpectrum,
        candidate: &mut Candidate,
        noise: Option<&[f64]>,
        work: &mut usize,
    ) -> Result<()> {
        let points = &input.peaks;
        let i = candidate.region.center_index;
        if i == 0 || i + 1 >= points.len() {
            return Err(bad("iterative center lacks two neighbors"));
        }
        spend(work, 3)?;
        let left_spacing = points[i].mz - points[i - 1].mz;
        let right_spacing = points[i + 1].mz - points[i].mz;
        if self.check_width_internally
            && (left_spacing > self.peak_width || right_spacing > self.peak_width)
        {
            candidate.valid = false;
            return Ok(());
        }
        let spacing = finite(self.spacing_difference * left_spacing.min(right_spacing))?;
        let mut left = i - 1;
        let mut right = i + 1;
        // Central sample and immediate neighbors are unconditional. All further
        // extensions use strict spacing/falling/width checks, with S/N equality passing.
        while left > 0 {
            spend(work, 1)?;
            let j = left - 1;
            if !(points[left].mz - points[j].mz < spacing
                && (points[j].intensity < points[left].intensity
                    || points[i].mz - points[j].mz < self.peak_width))
                || noise.is_some_and(|sn| sn[j] < self.signal_to_noise)
            {
                break;
            }
            left = j;
        }
        while right + 1 < points.len() {
            spend(work, 1)?;
            let j = right + 1;
            if !(points[j].mz - points[right].mz < spacing
                && (points[j].intensity < points[right].intensity
                    || points[j].mz - points[i].mz < self.peak_width))
                || noise.is_some_and(|sn| sn[j] < self.signal_to_noise)
            {
                break;
            }
            right = j;
        }
        let mut weighted = 0.0;
        let mut integrated = 0.0;
        for point in &points[left..=right] {
            spend(work, 1)?;
            weighted = finite(weighted + point.mz * f64::from(point.intensity))?;
            integrated = finite(integrated + f64::from(point.intensity))?;
        }
        if integrated <= 0.0 {
            return Err(bad("iterative candidate has zero integrated intensity"));
        }
        let weighted_mz = finite(weighted / integrated)?;
        candidate.intensity = integrated;
        candidate.mz = checked_f32(weighted_mz)?;
        candidate.region.left_index = left;
        candidate.region.right_index = right;
        let mut min_diff = (weighted_mz - points[i].mz).abs();
        let mut center = i;
        let mut m = 1;
        while m < i && points[left].mz < points[i - m].mz {
            spend(work, 1)?;
            let difference = (weighted_mz - points[i - m].mz).abs();
            if difference < min_diff {
                min_diff = difference;
                center = i - m;
            }
            m += 1;
        }
        m = 1;
        // The i-m > 0 condition is intentionally retained for the right search.
        // The additional upper guard rules out the source's unchecked read.
        while m < i && m < points.len() - i && points[right].mz > points[i + m].mz {
            spend(work, 1)?;
            let difference = (weighted_mz - points[i + m].mz).abs();
            if difference < min_diff {
                min_diff = difference;
                center = i + m;
            }
            m += 1;
        }
        candidate.region.center_index = center;
        Ok(())
    }

    fn output(
        &self,
        input: &MSSpectrum,
        candidates: Vec<Candidate>,
    ) -> Result<IterativePickingResult> {
        let mut spectrum = input.clone();
        spectrum.peaks.clear();
        let omitted_arrays = input
            .float_data_arrays
            .iter()
            .map(|a| a.name.clone())
            .chain(input.integer_data_arrays.iter().map(|a| a.name.clone()))
            .chain(input.string_data_arrays.iter().map(|a| a.name.clone()))
            .collect();
        spectrum.float_data_arrays.clear();
        spectrum.integer_data_arrays.clear();
        spectrum.string_data_arrays.clear();
        spectrum.spectrum_type = SpectrumType::Centroid;
        let mut intensity = DataArray::new("IntegratedIntensity", Vec::new());
        let mut left = DataArray::new("leftWidth", Vec::new());
        let mut right = DataArray::new("rightWidth", Vec::new());
        let mut boundaries = Vec::with_capacity(candidates.len());
        let mut regions = Vec::with_capacity(candidates.len());
        for candidate in candidates {
            let value = checked_f32(candidate.intensity)?;
            let lo = input.peaks[candidate.region.left_index].mz;
            let hi = input.peaks[candidate.region.right_index].mz;
            spectrum
                .peaks
                .push(Peak1D::new(f64::from(candidate.mz), value));
            intensity.data.push(value);
            left.data.push(checked_f32(lo)?);
            right.data.push(checked_f32(hi)?);
            boundaries.push(PeakBoundary { min: lo, max: hi });
            regions.push(candidate.region);
        }
        spectrum.float_data_arrays = vec![intensity, left, right];
        Ok(IterativePickingResult {
            picked: PickedSpectrum {
                spectrum,
                boundaries,
                omitted_arrays,
            },
            regions,
        })
    }
}

impl SpectrumFilter for PeakPickerIterative {
    fn filter_spectrum(&self, input: &mut MSSpectrum) -> Result<()> {
        *input = self.pick_spectrum(input)?.picked.spectrum;
        Ok(())
    }
    fn filter_experiment(&self, input: &mut MSExperiment) -> Result<()> {
        *input = self.pick_experiment(input)?.experiment;
        Ok(())
    }
}

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("iterative picker numerical expression overflows"))
    }
}
fn checked_f32(value: f64) -> Result<f32> {
    let value = value as f32;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("iterative picker value exceeds f32 range"))
    }
}
fn spend(work: &mut usize, amount: usize) -> Result<()> {
    *work = work
        .checked_sub(amount)
        .ok_or_else(|| bad("iterative picker work limit exceeded"))?;
    Ok(())
}
fn charge_sort(work: &mut usize, n: usize) -> Result<()> {
    if n < 2 {
        return Ok(());
    }
    let levels = usize::BITS as usize - (n - 1).leading_zeros() as usize;
    spend(
        work,
        n.checked_mul(levels)
            .and_then(|n| n.checked_mul(2))
            .ok_or_else(|| bad("iterative sort work overflows"))?,
    )
}

#[cfg(test)]
mod review_tests;
