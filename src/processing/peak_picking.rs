// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! High-resolution centroiding with histogram noise estimation and natural cubic splines.
//! See `docs/PEAK_PICKING_SUPPORT.md` for source behavior and deliberate limits.

mod noise;
mod spline;
pub use noise::{NoiseEstimates, NoiseHistogramRange, SignalToNoiseEstimatorMedian};
pub use spline::CubicSpline2d;

use super::{SpectrumFilter, checked_intensity};
use crate::kernel::{
    ChromatogramPeak, DataArray, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, SpectrumType,
};
use crate::{Error, Result};

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn validate_signal(x: &[f64], y: &[f64], max_points: usize) -> Result<()> {
    if x.len() != y.len() || max_points == 0 || x.len() > max_points {
        return Err(bad("signal arrays differ in length or exceed point limit"));
    }
    if x.iter().chain(y).any(|v| !v.is_finite()) || y.iter().any(|&v| v < 0.0) {
        return Err(bad(
            "peak picking requires finite coordinates and nonnegative finite intensities",
        ));
    }
    if x.windows(2).any(|v| v[0] > v[1]) {
        return Err(Error::UnsortedData);
    }
    if x.windows(2).any(|v| v[0] == v[1]) {
        return Err(bad("peak picking requires distinct coordinates"));
    }
    Ok(())
}

/// Source extension boundary, including a final rejected missing sample.
/// Values use Th for spectra and seconds for chromatograms.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PeakBoundary {
    pub min: f64,
    pub max: f64,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FwhmUnit {
    Absolute,
    Ppm,
}
#[derive(Clone, Debug, PartialEq)]
pub struct PickedSpectrum {
    pub spectrum: MSSpectrum,
    pub boundaries: Vec<PeakBoundary>,
    /// Profile annotations that have no defined centroid aggregation rule.
    pub omitted_arrays: Vec<String>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct PickedChromatogram {
    pub chromatogram: MSChromatogram,
    pub boundaries: Vec<PeakBoundary>,
    pub omitted_arrays: Vec<String>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct PickedExperiment {
    pub experiment: MSExperiment,
    /// One entry per input spectrum; None means copied without picking.
    pub spectrum_boundaries: Vec<Option<Vec<PeakBoundary>>>,
    pub chromatogram_boundaries: Vec<Vec<PeakBoundary>>,
    pub omitted_spectrum_arrays: Vec<Vec<String>>,
    pub omitted_chromatogram_arrays: Vec<Vec<String>>,
}

/// OpenMS PeakPickerHiRes options. Zero signal-to-noise disables estimation.
#[derive(Clone, Debug)]
pub struct PeakPickerHiRes {
    pub signal_to_noise: f64,
    pub noise_estimator: SignalToNoiseEstimatorMedian,
    /// Zero disables this spacing constraint.
    pub spacing_difference: f64,
    /// Zero disables this spacing constraint.
    pub spacing_difference_gap: f64,
    pub missing: usize,
    pub allow_missing_flank: bool,
    pub report_fwhm: Option<FwhmUnit>,
    /// Empty means automatic type selection; otherwise pick these MS levels.
    pub ms_levels: Vec<u32>,
    pub check_spectrum_type: bool,
    /// Explicit array name for intensity-weighted ion mobility. None recognizes
    /// the source's "Ion Mobility" and "raw inverse reduced ion mobility array" names.
    pub ion_mobility_array: Option<String>,
    pub max_points: usize,
    /// Bound on apex candidates and extension samples for one record.
    pub max_work: usize,
}
impl Default for PeakPickerHiRes {
    fn default() -> Self {
        Self {
            signal_to_noise: 0.0,
            noise_estimator: Default::default(),
            spacing_difference: 1.5,
            spacing_difference_gap: 4.0,
            missing: 1,
            allow_missing_flank: false,
            report_fwhm: None,
            ms_levels: Vec::new(),
            check_spectrum_type: true,
            ion_mobility_array: None,
            max_points: 1_000_000,
            max_work: 10_000_000,
        }
    }
}
struct PickedSignal {
    positions: Vec<f64>,
    intensities: Vec<f32>,
    boundaries: Vec<PeakBoundary>,
    fwhm: Vec<f32>,
    mobility: Vec<f32>,
}
impl PeakPickerHiRes {
    fn validate(&self) -> Result<()> {
        if [
            self.signal_to_noise,
            self.spacing_difference,
            self.spacing_difference_gap,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
            || self.max_points == 0
            || self.max_work == 0
            || self.ms_levels.contains(&0)
        {
            return Err(bad("invalid peak picker parameters or resource limits"));
        }
        self.noise_estimator.validate()
    }
    fn pick_signal(
        &self,
        x: &[f64],
        y: &[f64],
        mut check_spacings: bool,
        mobility: Option<&[f32]>,
    ) -> Result<PickedSignal> {
        self.validate()?;
        validate_signal(x, y, self.max_points)?;
        if let Some(values) = mobility {
            if values.len() != x.len() || values.iter().any(|v| !v.is_finite()) {
                return Err(bad(
                    "ion mobility must be finite and aligned with input peaks",
                ));
            }
        }
        let mut out = PickedSignal {
            positions: Vec::new(),
            intensities: Vec::new(),
            boundaries: Vec::new(),
            fwhm: Vec::new(),
            mobility: Vec::new(),
        };
        if x.len() < 5 {
            return Ok(out);
        }
        let spacing = if self.spacing_difference == 0.0 {
            f64::INFINITY
        } else {
            self.spacing_difference
        };
        let gap = if self.spacing_difference_gap == 0.0 {
            f64::INFINITY
        } else {
            self.spacing_difference_gap
        };
        if spacing.is_infinite() && gap.is_infinite() {
            check_spacings = false;
        }
        let estimates = if self.signal_to_noise > 0.0 {
            Some(self.noise_estimator.estimate(x, y)?)
        } else {
            None
        };
        let passes_sn = |i: usize| {
            estimates
                .as_ref()
                .is_none_or(|s| s.signal_to_noise[i] >= self.signal_to_noise)
        };
        let mut used = 0usize;
        let mut charge = || -> Result<()> {
            if used == self.max_work {
                Err(bad("peak picking exceeds configured work limit"))
            } else {
                used += 1;
                Ok(())
            }
        };
        let mut i = 2;
        while i < x.len() - 2 {
            charge()?;
            if y[i - 1].abs() < f64::EPSILON || y[i + 1].abs() < f64::EPSILON {
                i += 1;
                continue;
            }
            let min_spacing = (x[i] - x[i - 1]).min(x[i + 1] - x[i]);
            if !min_spacing.is_finite() {
                return Err(bad("peak spacing overflows"));
            }
            let left_neighbor = !check_spacings || x[i] - x[i - 1] < spacing * min_spacing;
            let right_neighbor = !check_spacings || x[i + 1] - x[i] < spacing * min_spacing;
            let spacing_ok = if self.allow_missing_flank {
                left_neighbor || right_neighbor
            } else {
                left_neighbor && right_neighbor
            };
            if !(y[i] > y[i - 1]
                && y[i] > y[i + 1]
                && passes_sn(i)
                && passes_sn(i - 1)
                && passes_sn(i + 1)
                && spacing_ok)
            {
                i += 1;
                continue;
            }
            let left2 = x[i - 1] - x[i - 2] < spacing * min_spacing;
            let right2 = x[i + 2] - x[i + 1] < spacing * min_spacing;
            let outer_spacing = !check_spacings
                || if self.allow_missing_flank {
                    left2 || right2
                } else {
                    left2 && right2
                };
            if y[i - 1] < y[i - 2]
                && y[i + 1] < y[i + 2]
                && passes_sn(i - 2)
                && passes_sn(i + 2)
                && outer_spacing
            {
                i += 2;
                continue;
            }
            // Keep input indices so the spline and weighted mobility share exactly the same support.
            let mut support = vec![i];
            if left_neighbor {
                support.insert(0, i - 1);
            }
            if right_neighbor {
                support.push(i + 1);
            }
            let mut left_boundary = i - 1;
            let mut k = 2;
            let mut missing = 0;
            let mut previous_zero = false;
            let mut left_support = Vec::new();
            let mut first = support[0];
            while k <= i
                && !previous_zero
                && missing <= self.missing
                && y[i - k] <= y[first]
                && (!check_spacings || x[first] - x[i - k] < gap * min_spacing)
            {
                charge()?;
                let j = i - k;
                let good =
                    passes_sn(j) && (!check_spacings || x[first] - x[j] < spacing * min_spacing);
                if !good {
                    missing = missing
                        .checked_add(1)
                        .ok_or_else(|| bad("missing-point count overflow"))?;
                }
                if good || missing <= self.missing {
                    left_support.push(j);
                    first = j;
                }
                previous_zero = y[j] == 0.0;
                left_boundary = j;
                k += 1;
            }
            left_support.reverse();
            left_support.extend(support);
            support = left_support;
            let mut right_boundary = i + 1;
            k = 2;
            missing = 0;
            previous_zero = false;
            let mut last = *support.last().expect("peak core exists");
            while k < x.len() - i
                && !previous_zero
                && missing <= self.missing
                && y[i + k] <= y[last]
                && (!check_spacings || x[i + k] - x[last] < gap * min_spacing)
            {
                charge()?;
                let j = i + k;
                let good =
                    passes_sn(j) && (!check_spacings || x[j] - x[last] < spacing * min_spacing);
                if !good {
                    missing = missing
                        .checked_add(1)
                        .ok_or_else(|| bad("missing-point count overflow"))?;
                }
                if good || missing <= self.missing {
                    support.push(j);
                    last = j;
                }
                previous_zero = y[j] == 0.0;
                right_boundary = j;
                k += 1;
            }
            if support.len() < 3 {
                i += 1;
                continue;
            }
            let sx: Vec<_> = support.iter().map(|&j| x[j]).collect();
            let sy: Vec<_> = support.iter().map(|&j| y[j]).collect();
            let spline = CubicSpline2d::with_max_points(&sx, &sy, self.max_points)?;
            let left = if left_neighbor { x[i - 1] } else { x[i] };
            let right = if right_neighbor { x[i + 1] } else { x[i] };
            let (position, intensity) = spline.peak_maximum(left, right, 1e-6)?;
            if intensity <= 0.0 {
                return Err(bad("peak spline maximum is not positive"));
            }
            if let Some(unit) = self.report_fwhm {
                let (lo, hi) = spline.domain();
                let left = half_height(&spline, lo, position, intensity / 2.0)?;
                let right = half_height(&spline, hi, position, intensity / 2.0)?;
                let width = right - left;
                let width = match unit {
                    FwhmUnit::Absolute => width,
                    FwhmUnit::Ppm => {
                        if position <= 0.0 {
                            return Err(bad("ppm FWHM needs positive centroid position"));
                        }
                        width / position * 1e6
                    }
                };
                out.fwhm.push(checked_intensity(width)?);
            }
            if let Some(values) = mobility {
                let total = sy.iter().sum::<f64>();
                let weighted = support
                    .iter()
                    .map(|&j| f64::from(values[j]) * y[j])
                    .sum::<f64>();
                out.mobility.push(checked_intensity(weighted / total)?);
            }
            out.positions.push(position);
            out.intensities.push(checked_intensity(intensity)?);
            out.boundaries.push(PeakBoundary {
                min: x[left_boundary],
                max: x[right_boundary],
            });
            i += k;
        }
        Ok(out)
    }
    /// Pick a profile spectrum with source spacing checks enabled.
    pub fn pick_spectrum(&self, input: &MSSpectrum) -> Result<PickedSpectrum> {
        self.pick_spectrum_with_spacing(input, true)
    }
    pub fn pick_spectrum_with_spacing(
        &self,
        input: &MSSpectrum,
        check_spacings: bool,
    ) -> Result<PickedSpectrum> {
        self.pick_spectrum_with_acquisition(
            input,
            check_spacings,
            &mut super::AcquisitionCopies::default(),
        )
    }
    pub(super) fn pick_spectrum_with_acquisition(
        &self,
        input: &MSSpectrum,
        check_spacings: bool,
        copies: &mut super::AcquisitionCopies,
    ) -> Result<PickedSpectrum> {
        self.validate()?;
        if input.len() > self.max_points {
            return Err(bad("spectrum exceeds peak picker point limit"));
        }
        input.validate()?;
        let mut mobility = None;
        for (index, array) in input.float_data_arrays.iter().enumerate() {
            let selected = match &self.ion_mobility_array {
                Some(name) => array.name == *name,
                None => matches!(
                    array.name.as_str(),
                    "Ion Mobility" | "raw inverse reduced ion mobility array"
                ),
            };
            if selected && mobility.replace(index).is_some() {
                return Err(bad("multiple matching ion mobility arrays"));
            }
        }
        if self.ion_mobility_array.is_some() && mobility.is_none() {
            return Err(bad("requested ion mobility array is missing"));
        }
        let x: Vec<_> = input.peaks.iter().map(|p| p.mz).collect();
        let y: Vec<_> = input.peaks.iter().map(|p| f64::from(p.intensity)).collect();
        let picked = self.pick_signal(
            &x,
            &y,
            check_spacings,
            mobility.map(|j| input.float_data_arrays[j].data.as_slice()),
        )?;
        let omitted_arrays = input
            .float_data_arrays
            .iter()
            .enumerate()
            .filter(|(j, _)| Some(*j) != mobility)
            .map(|(_, a)| a.name.clone())
            .chain(input.integer_data_arrays.iter().map(|a| a.name.clone()))
            .chain(input.string_data_arrays.iter().map(|a| a.name.clone()))
            .collect();
        copies.spectrum(input)?;
        let mut output = input.clone();
        output.peaks = picked
            .positions
            .iter()
            .zip(&picked.intensities)
            .map(|(&x, &y)| Peak1D::new(x, y))
            .collect();
        output.spectrum_type = SpectrumType::Centroid;
        output.float_data_arrays.clear();
        output.integer_data_arrays.clear();
        output.string_data_arrays.clear();
        if let Some(j) = mobility {
            let mut array =
                DataArray::new(input.float_data_arrays[j].name.clone(), picked.mobility);
            input.float_data_arrays[j].copy_description_to(&mut array);
            output.float_data_arrays.push(array);
        }
        if let Some(unit) = self.report_fwhm {
            output
                .float_data_arrays
                .push(DataArray::new(fwhm_name(unit), picked.fwhm));
        }
        Ok(PickedSpectrum {
            spectrum: output,
            boundaries: picked.boundaries,
            omitted_arrays,
        })
    }
    /// Chromatogram picking disables spacing checks, as in OpenMS.
    pub fn pick_chromatogram(&self, input: &MSChromatogram) -> Result<PickedChromatogram> {
        self.pick_chromatogram_with_spacing(input, false)
    }
    pub fn pick_chromatogram_with_spacing(
        &self,
        input: &MSChromatogram,
        check_spacings: bool,
    ) -> Result<PickedChromatogram> {
        self.pick_chromatogram_with_acquisition(
            input,
            check_spacings,
            &mut super::AcquisitionCopies::default(),
        )
    }
    pub(super) fn pick_chromatogram_with_acquisition(
        &self,
        input: &MSChromatogram,
        check_spacings: bool,
        copies: &mut super::AcquisitionCopies,
    ) -> Result<PickedChromatogram> {
        self.validate()?;
        if input.len() > self.max_points {
            return Err(bad("chromatogram exceeds peak picker point limit"));
        }
        input.validate()?;
        let x: Vec<_> = input.peaks.iter().map(|p| p.rt).collect();
        let y: Vec<_> = input.peaks.iter().map(|p| f64::from(p.intensity)).collect();
        let picked = self.pick_signal(&x, &y, check_spacings, None)?;
        let omitted_arrays = input
            .float_data_arrays
            .iter()
            .map(|a| a.name.clone())
            .chain(input.integer_data_arrays.iter().map(|a| a.name.clone()))
            .chain(input.string_data_arrays.iter().map(|a| a.name.clone()))
            .collect();
        copies.chromatogram(input)?;
        let mut output = input.clone();
        output.peaks = picked
            .positions
            .iter()
            .zip(&picked.intensities)
            .map(|(&x, &y)| ChromatogramPeak::new(x, y))
            .collect();
        output.float_data_arrays.clear();
        output.integer_data_arrays.clear();
        output.string_data_arrays.clear();
        if let Some(unit) = self.report_fwhm {
            output
                .float_data_arrays
                .push(DataArray::new(fwhm_name(unit), picked.fwhm));
        }
        Ok(PickedChromatogram {
            chromatogram: output,
            boundaries: picked.boundaries,
            omitted_arrays,
        })
    }
    /// Pick all chromatograms and selected spectra without mutating input.
    /// Unknown spectrum types are inferred using the source's five-apex heuristic.
    pub fn pick_experiment(&self, input: &MSExperiment) -> Result<PickedExperiment> {
        self.validate()?;
        input.validate()?;
        let mut copies = super::AcquisitionCopies::default();
        copies.experiment(input)?;
        let mut result = PickedExperiment {
            experiment: input.clone(),
            spectrum_boundaries: Vec::new(),
            chromatogram_boundaries: Vec::new(),
            omitted_spectrum_arrays: Vec::new(),
            omitted_chromatogram_arrays: Vec::new(),
        };
        for (i, spectrum) in input.spectra.iter().enumerate() {
            let selected = self.ms_levels.is_empty() || self.ms_levels.contains(&spectrum.ms_level);
            if !selected {
                result.spectrum_boundaries.push(None);
                result.omitted_spectrum_arrays.push(Vec::new());
                continue;
            }
            let kind = if spectrum.spectrum_type == SpectrumType::Unknown {
                estimate_spectrum_type_with_limit(spectrum, self.max_points)?
            } else {
                spectrum.spectrum_type
            };
            if self.ms_levels.is_empty() && kind == SpectrumType::Centroid {
                result.spectrum_boundaries.push(None);
                result.omitted_spectrum_arrays.push(Vec::new());
                continue;
            }
            if !self.ms_levels.is_empty()
                && self.check_spectrum_type
                && kind == SpectrumType::Centroid
            {
                return Err(bad("centroid spectrum selected but profile input required"));
            }
            let picked = self.pick_spectrum_with_acquisition(spectrum, true, &mut copies)?;
            result.experiment.spectra[i] = picked.spectrum;
            result.spectrum_boundaries.push(Some(picked.boundaries));
            result.omitted_spectrum_arrays.push(picked.omitted_arrays);
        }
        for (i, chromatogram) in input.chromatograms.iter().enumerate() {
            let picked =
                self.pick_chromatogram_with_acquisition(chromatogram, false, &mut copies)?;
            result.experiment.chromatograms[i] = picked.chromatogram;
            result.chromatogram_boundaries.push(picked.boundaries);
            result
                .omitted_chromatogram_arrays
                .push(picked.omitted_arrays);
        }
        Ok(result)
    }
    pub fn filter_chromatogram(&self, input: &mut MSChromatogram) -> Result<()> {
        *input = self.pick_chromatogram(input)?.chromatogram;
        Ok(())
    }
}
impl SpectrumFilter for PeakPickerHiRes {
    fn filter_spectrum(&self, input: &mut MSSpectrum) -> Result<()> {
        *input = self.pick_spectrum(input)?.spectrum;
        Ok(())
    }
    fn filter_experiment(&self, input: &mut MSExperiment) -> Result<()> {
        *input = self.pick_experiment(input)?.experiment;
        Ok(())
    }
}
fn fwhm_name(unit: FwhmUnit) -> &'static str {
    match unit {
        FwhmUnit::Absolute => "FWHM",
        FwhmUnit::Ppm => "FWHM_ppm",
    }
}
fn half_height(spline: &CubicSpline2d, mut edge: f64, mut center: f64, height: f64) -> Result<f64> {
    if spline.eval(edge)? > height {
        return Ok(edge);
    }
    for _ in 0..128 {
        let mid = edge / 2.0 + center / 2.0;
        let value = spline.eval(mid)?;
        if (value - height).abs() <= 0.01 * height {
            return Ok(mid);
        }
        if mid == edge || mid == center {
            return Err(bad("FWHM cannot converge at coordinate precision"));
        }
        if value < height {
            edge = mid;
        } else {
            center = mid;
        }
    }
    Err(bad("FWHM did not converge within 128 iterations"))
}

/// Infer profile/centroid from up to five high-intensity peak shoulders.
/// Fewer than five samples are Unknown; no positive evidence is Centroid.
pub fn estimate_spectrum_type(input: &MSSpectrum) -> Result<SpectrumType> {
    estimate_spectrum_type_with_limit(input, 1_000_000)
}
fn estimate_spectrum_type_with_limit(
    input: &MSSpectrum,
    max_points: usize,
) -> Result<SpectrumType> {
    if input.len() > max_points {
        return Err(bad("spectrum type estimation exceeds point limit"));
    }
    input.validate()?;
    let x: Vec<_> = input.peaks.iter().map(|p| p.mz).collect();
    let mut y: Vec<_> = input.peaks.iter().map(|p| f64::from(p.intensity)).collect();
    validate_signal(&x, &y, max_points)?;
    Ok(crate::kernel::spectrum_type::estimate(&x, &mut y))
}
