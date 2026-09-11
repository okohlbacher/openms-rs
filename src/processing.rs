// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Peak filters ported from OpenMS `PROCESSING`, with typed parameters.
//!
//! Finite signed intensities are accepted. Operations validate before mutation.
//! All-zero normalization is a no-op; other zero divisors are errors rather than
//! producing NaN as the C++ implementation does. Sorting and filtering preserve
//! the association between peaks and their auxiliary data arrays.

use crate::kernel::{ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D};
use crate::{Error, Result};

pub mod baseline;
pub mod chromatogram;
pub mod deisotoping;
pub mod iterative;
pub mod mean_noise;
pub mod peak_picking;
pub mod smoothing;
pub mod window_mower;

/// A spectrum transformation with an atomic experiment convenience method.
pub trait SpectrumFilter {
    /// Transform a spectrum; errors leave it unchanged.
    fn filter_spectrum(&self, spectrum: &mut MSSpectrum) -> Result<()>;

    /// Transform every spectrum, retaining chromatograms and experiment metadata.
    ///
    /// Uses a temporary copy so an error leaves the whole experiment unchanged.
    fn filter_experiment(&self, experiment: &mut MSExperiment) -> Result<()> {
        AcquisitionCopies::default().spectra(&experiment.spectra)?;
        let mut spectra = experiment.spectra.clone();
        for spectrum in &mut spectra {
            self.filter_spectrum(spectrum)?;
        }
        experiment.spectra = spectra;
        Ok(())
    }
}

// Fixed acquisition-copy ceiling, separate from existing numerical algorithm
// budgets. Carry one ledger through nested bundled transformations and batches.
// Processing records behind Arc remain shared; their payload is not copied.
pub(super) struct AcquisitionCopies {
    work: usize,
    bytes: usize,
}
impl Default for AcquisitionCopies {
    fn default() -> Self {
        Self {
            work: 50_000_000,
            bytes: 256 * 1024 * 1024,
        }
    }
}
impl AcquisitionCopies {
    fn visit(&mut self, count: usize) -> Result<()> {
        self.work = self.work.checked_sub(count).ok_or_else(|| {
            Error::InvalidValue("processing acquisition copy work limit exceeded".into())
        })?;
        Ok(())
    }
    pub(super) fn spectrum(&mut self, input: &MSSpectrum) -> Result<()> {
        self.visit(1)?;
        input.acquisition_with_budget(&mut self.work, &mut self.bytes)
    }
    pub(super) fn chromatogram(&mut self, input: &MSChromatogram) -> Result<()> {
        self.visit(1)?;
        input.acquisition_with_budget(&mut self.work, &mut self.bytes)
    }
    pub(super) fn spectra(&mut self, spectra: &[MSSpectrum]) -> Result<()> {
        self.visit(spectra.len())?;
        for spectrum in spectra {
            self.spectrum(spectrum)?;
        }
        Ok(())
    }
    pub(super) fn experiment(&mut self, input: &MSExperiment) -> Result<()> {
        input
            .settings
            .with_budget(&mut self.work, &mut self.bytes)?;
        self.spectra(&input.spectra)?;
        self.visit(input.chromatograms.len())?;
        for chromatogram in &input.chromatograms {
            self.chromatogram(chromatogram)?;
        }
        Ok(())
    }
}

/// Divisor used by [`Normalizer`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum NormalizationMethod {
    /// Scale the maximum intensity to one (`to_one` in OpenMS).
    #[default]
    ToOne,
    /// Scale the sum of intensities to one (`to_TIC` in OpenMS).
    ToTic,
}

/// Spectrum-wise maximum or total-ion-current normalization.
#[derive(Clone, Copy, Debug, Default)]
pub struct Normalizer {
    /// Select maximum or TIC normalization.
    pub method: NormalizationMethod,
}

impl SpectrumFilter for Normalizer {
    fn filter_spectrum(&self, spectrum: &mut MSSpectrum) -> Result<()> {
        spectrum.validate()?;
        if spectrum.is_empty() {
            return Ok(());
        }
        let divisor = match self.method {
            NormalizationMethod::ToOne => spectrum
                .peaks
                .iter()
                .map(|p| f64::from(p.intensity))
                .fold(f64::NEG_INFINITY, f64::max),
            NormalizationMethod::ToTic => {
                spectrum.peaks.iter().map(|p| f64::from(p.intensity)).sum()
            }
        };
        if divisor == 0.0 {
            return if spectrum.peaks.iter().all(|p| p.intensity == 0.0) {
                Ok(())
            } else {
                Err(Error::InvalidValue("normalization divisor is zero".into()))
            };
        }
        let intensities: Result<Vec<f32>> = spectrum
            .peaks
            .iter()
            .map(|p| checked_intensity(f64::from(p.intensity) / divisor))
            .collect();
        for (peak, intensity) in spectrum.peaks.iter_mut().zip(intensities?) {
            peak.intensity = intensity;
        }
        Ok(())
    }
}

/// Retain peaks with intensity greater than or equal to a threshold.
#[derive(Clone, Copy, Debug)]
pub struct ThresholdMower {
    /// Inclusive intensity threshold; defaults to the source's 0.05.
    pub threshold: f64,
}

impl Default for ThresholdMower {
    fn default() -> Self {
        Self { threshold: 0.05 }
    }
}

impl SpectrumFilter for ThresholdMower {
    fn filter_spectrum(&self, spectrum: &mut MSSpectrum) -> Result<()> {
        if !self.threshold.is_finite() {
            return Err(Error::InvalidValue("threshold must be finite".into()));
        }
        spectrum.validate()?;
        spectrum.retain_peaks(|p| f64::from(p.intensity) >= self.threshold)
    }
}

/// Keep the N most intense peaks, sorted by decreasing intensity.
///
/// Like OpenMS, leaves the original order unchanged when `len() <= n`.
/// Equal intensities keep their input order in this port.
#[derive(Clone, Copy, Debug)]
pub struct NLargest {
    /// Maximum number of peaks to retain.
    pub n: usize,
}

impl Default for NLargest {
    fn default() -> Self {
        Self { n: 200 }
    }
}

impl SpectrumFilter for NLargest {
    fn filter_spectrum(&self, spectrum: &mut MSSpectrum) -> Result<()> {
        spectrum.validate()?;
        if spectrum.len() <= self.n {
            return Ok(());
        }
        let mut indices: Vec<usize> = (0..spectrum.len()).collect();
        indices.sort_by(|&a, &b| {
            spectrum.peaks[b]
                .intensity
                .partial_cmp(&spectrum.peaks[a].intensity)
                .unwrap()
        });
        indices.truncate(self.n);
        spectrum.select(&indices)
    }
}

/// Square-root intensities, clamping negative values to zero as OpenMS does.
#[derive(Clone, Copy, Debug, Default)]
pub struct SqrtScaler;

impl SpectrumFilter for SqrtScaler {
    fn filter_spectrum(&self, spectrum: &mut MSSpectrum) -> Result<()> {
        spectrum.validate()?;
        for peak in &mut spectrum.peaks {
            peak.intensity = f64::from(peak.intensity).max(0.0).sqrt() as f32;
        }
        Ok(())
    }
}

/// OpenMS rank scaling (equal intensities share a rank).
///
/// Sorts by increasing original intensity. Preserves OpenMS's unusual rank
/// offset: with N peaks the maximum is assigned N, except an all-zero maximum
/// is assigned N+1. Each distinct lower intensity decreases the rank by one.
#[derive(Clone, Copy, Debug, Default)]
pub struct RankScaler;

impl SpectrumFilter for RankScaler {
    fn filter_spectrum(&self, spectrum: &mut MSSpectrum) -> Result<()> {
        spectrum.sort_by_intensity(false)?;
        let mut count = spectrum.len() + 1;
        let mut last_intensity = 0.0;
        for peak in spectrum.peaks.iter_mut().rev() {
            if peak.intensity != last_intensity {
                count -= 1;
            }
            last_intensity = peak.intensity;
            peak.intensity = count as f32;
        }
        Ok(())
    }
}

/// Intensity-conserving redistribution onto a uniformly spaced grid.
///
/// Port of the absolute-spacing paths of OpenMS `LinearResamplerAlign`.
/// Per-peak auxiliary arrays cannot be meaningfully resampled and are rejected
/// when nonempty. Metadata and precursor information are preserved.
#[derive(Clone, Copy, Debug)]
pub struct LinearResamplerAlign {
    /// Absolute spacing in Th for spectra or seconds for chromatograms.
    pub spacing: f64,
    /// Bound on generated grid size (default ten million points).
    pub max_points: usize,
}

impl Default for LinearResamplerAlign {
    fn default() -> Self {
        Self {
            spacing: 0.05,
            max_points: 10_000_000,
        }
    }
}

impl LinearResamplerAlign {
    /// Create a resampler with finite, strictly positive absolute spacing.
    pub fn new(spacing: f64) -> Result<Self> {
        let value = Self {
            spacing,
            ..Self::default()
        };
        value.validate_options()?;
        Ok(value)
    }

    fn validate_options(&self) -> Result<()> {
        if !self.spacing.is_finite() || self.spacing <= 0.0 || self.max_points == 0 {
            return Err(Error::InvalidValue(
                "spacing and maximum grid size must be positive and finite".into(),
            ));
        }
        Ok(())
    }

    fn grid(&self, start: f64, end: f64) -> Result<Vec<f64>> {
        self.validate_options()?;
        if !start.is_finite() || !end.is_finite() || end < start {
            return Err(Error::InvalidValue(
                "resampling bounds must be finite and ordered".into(),
            ));
        }
        let count = ((end - start) / self.spacing).ceil() + 1.0;
        if !count.is_finite() || count > self.max_points as f64 {
            return Err(Error::InvalidValue(
                "resampling grid exceeds maximum point count".into(),
            ));
        }
        let grid: Vec<f64> = (0..count as usize)
            .map(|i| start + i as f64 * self.spacing)
            .collect();
        validate_grid(&grid)?;
        Ok(grid)
    }

    /// Redistribute peaks from the first through the last input m/z.
    pub fn raster(&self, spectrum: &mut MSSpectrum) -> Result<()> {
        self.validate_options()?;
        spectrum.validate()?;
        let Some(first) = spectrum.peaks.first() else {
            return Ok(());
        };
        self.raster_align(spectrum, first.mz, spectrum.peaks.last().unwrap().mz)
    }

    /// Redistribute peaks within the inclusive bounds, omitting outside peaks.
    /// The last generated grid point may lie just beyond `end`, as in OpenMS.
    /// Empty spectra remain empty, matching the upstream early return.
    pub fn raster_align(&self, spectrum: &mut MSSpectrum, start: f64, end: f64) -> Result<()> {
        self.validate_options()?;
        spectrum.validate()?;
        ensure_no_arrays(spectrum)?;
        if spectrum.is_empty() {
            return Ok(());
        }
        let grid = self.grid(start, end)?;
        validate_peaks(&spectrum.peaks)?;
        let selected: Vec<Peak1D> = spectrum
            .peaks
            .iter()
            .copied()
            .filter(|p| p.mz >= start && p.mz <= end)
            .collect();
        spectrum.peaks = resample_to_grid(&selected, &grid)?;
        Ok(())
    }

    /// Resample a chromatogram using spacing measured in seconds.
    pub fn raster_chromatogram(&self, chromatogram: &mut MSChromatogram) -> Result<()> {
        self.validate_options()?;
        chromatogram.validate()?;
        if chromatogram
            .float_data_arrays
            .iter()
            .any(|a| !a.data.is_empty())
            || chromatogram
                .integer_data_arrays
                .iter()
                .any(|a| !a.data.is_empty())
            || chromatogram
                .string_data_arrays
                .iter()
                .any(|a| !a.data.is_empty())
        {
            return Err(Error::Unsupported(
                "resampling nonempty auxiliary arrays".into(),
            ));
        }
        let Some(first) = chromatogram.peaks.first() else {
            return Ok(());
        };
        let grid = self.grid(first.rt, chromatogram.peaks.last().unwrap().rt)?;
        let input: Vec<Peak1D> = chromatogram
            .peaks
            .iter()
            .map(|p| Peak1D::new(p.rt, p.intensity))
            .collect();
        chromatogram.peaks = resample_to_grid(&input, &grid)?
            .into_iter()
            .map(|p| ChromatogramPeak::new(p.mz, p.intensity))
            .collect();
        Ok(())
    }
}

/// Redistribute to a strictly increasing grid, conserving total intensity.
///
/// Intensities outside the grid are added to the nearest boundary. Coordinates
/// must be finite; input peaks must be sorted (duplicate positions are allowed).
pub fn resample_to_grid(input: &[Peak1D], grid: &[f64]) -> Result<Vec<Peak1D>> {
    validate_grid(grid)?;
    validate_peaks(input)?;
    let mut output: Vec<Peak1D> = grid.iter().map(|&mz| Peak1D::new(mz, 0.0)).collect();
    let mut right = 0;
    for peak in input {
        while right < grid.len() && grid[right] < peak.mz {
            right += 1;
        }
        if right == 0 {
            output[0].intensity =
                checked_intensity(f64::from(output[0].intensity) + f64::from(peak.intensity))?;
        } else if right == grid.len() {
            let last = output.last_mut().unwrap();
            last.intensity =
                checked_intensity(f64::from(last.intensity) + f64::from(peak.intensity))?;
        } else {
            let fraction = (peak.mz - grid[right - 1]) / (grid[right] - grid[right - 1]);
            output[right - 1].intensity = checked_intensity(
                f64::from(output[right - 1].intensity)
                    + f64::from(peak.intensity) * (1.0 - fraction),
            )?;
            output[right].intensity = checked_intensity(
                f64::from(output[right].intensity) + f64::from(peak.intensity) * fraction,
            )?;
        }
    }
    Ok(output)
}

fn ensure_no_arrays(spectrum: &MSSpectrum) -> Result<()> {
    if spectrum
        .float_data_arrays
        .iter()
        .any(|a| !a.data.is_empty())
        || spectrum
            .integer_data_arrays
            .iter()
            .any(|a| !a.data.is_empty())
        || spectrum
            .string_data_arrays
            .iter()
            .any(|a| !a.data.is_empty())
    {
        return Err(Error::Unsupported(
            "resampling nonempty auxiliary arrays".into(),
        ));
    }
    Ok(())
}

fn validate_grid(grid: &[f64]) -> Result<()> {
    if grid.is_empty()
        || grid.iter().any(|x| !x.is_finite())
        || grid
            .windows(2)
            .any(|w| w[0] >= w[1] || !(w[1] - w[0]).is_finite())
    {
        return Err(Error::InvalidValue(
            "grid must be nonempty, finite and strictly increasing with finite intervals".into(),
        ));
    }
    Ok(())
}

fn validate_peaks(peaks: &[Peak1D]) -> Result<()> {
    if peaks
        .iter()
        .any(|p| !p.mz.is_finite() || !p.intensity.is_finite())
    {
        return Err(Error::InvalidValue("peaks must be finite".into()));
    }
    if peaks.windows(2).any(|p| p[0].mz > p[1].mz) {
        return Err(Error::UnsortedData);
    }
    Ok(())
}

fn checked_intensity(value: f64) -> Result<f32> {
    let result = value as f32;
    if result.is_finite() {
        Ok(result)
    } else {
        Err(Error::InvalidValue("intensity overflow".into()))
    }
}

#[cfg(test)]
mod acquisition_copy_tests {
    use super::*;
    use crate::metadata::DataProcessing;
    use std::sync::Arc;

    #[test]
    fn acquisition_copy_work_and_bytes_are_cumulative_before_cloning() {
        let mut input = MSSpectrum::default();
        input.source_file.name = "x".repeat(100);
        let mut copies = AcquisitionCopies {
            work: 150,
            bytes: 1_000_000,
        };
        copies.spectrum(&input).unwrap();
        assert!(copies.spectrum(&input).is_err());
        let mut copies = AcquisitionCopies {
            work: 150,
            bytes: 1_000_000,
        };
        assert!(copies.spectra(&[input.clone(), input.clone()]).is_err());
        let mut copies = AcquisitionCopies {
            work: 1_000_000,
            bytes: 1,
        };
        assert!(copies.spectrum(&input).is_err());
        assert_eq!(input.source_file.name.len(), 100);
    }

    #[test]
    fn nested_picking_uses_the_same_acquisition_copy_ledger() {
        let mut spectrum = MSSpectrum::default();
        spectrum.source_file.name = "x".repeat(100);
        let input = MSExperiment {
            spectra: vec![spectrum],
            ..Default::default()
        };
        let mut copies = AcquisitionCopies {
            work: 150,
            bytes: 1_000_000,
        };
        copies.experiment(&input).unwrap();
        let mut staged = input.clone();
        let result = peak_picking::PeakPickerHiRes::default().pick_spectrum_with_acquisition(
            &input.spectra[0],
            true,
            &mut copies,
        );
        assert!(result.is_err());
        // The rejected picked value cannot replace the already staged record.
        assert_eq!(staged, input);
        let picked = peak_picking::PeakPickerHiRes::default()
            .pick_spectrum(&input.spectra[0])
            .unwrap();
        staged.spectra[0] = picked.spectrum;
        assert_eq!(staged.spectra[0].source_file, input.spectra[0].source_file);
    }

    #[test]
    fn acquisition_copy_counts_processing_handles_without_copying_shared_payload() {
        let mut record = DataProcessing::default();
        record.software.name = "shared".repeat(10_000);
        let record = Arc::new(record);
        let input = MSSpectrum {
            data_processing: vec![Arc::clone(&record)],
            ..Default::default()
        };
        let mut copies = AcquisitionCopies {
            work: 100,
            bytes: 4096,
        };
        copies.spectrum(&input).unwrap();
        assert_eq!(Arc::strong_count(&record), 2);
        let output = input.clone();
        assert!(Arc::ptr_eq(
            &input.data_processing[0],
            &output.data_processing[0]
        ));
        assert_eq!(Arc::strong_count(&record), 3);
    }
}
