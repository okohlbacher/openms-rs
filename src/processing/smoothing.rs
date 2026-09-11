// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Gaussian convolution and Savitzky–Golay least-squares smoothing.
//! Coordinates and aligned annotations remain unchanged. Gaussian convolution
//! integrates intervals; Savitzky–Golay fits equally spaced *sample indices*.

use super::{SpectrumFilter, checked_intensity};
use crate::kernel::{MSChromatogram, MSExperiment, MSSpectrum, SpectrumType};
use crate::{Error, Result};

/// Gaussian width is eight standard deviations, as in OpenMS's implementation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GaussianWidth {
    /// Full Gaussian width in coordinate units (Th or seconds).
    Absolute(f64),
    /// Full width in ppm, recalculated at each positive m/z position.
    Ppm(f64),
}

/// Numerical result before conversion to the kernel's f32 intensities.
#[derive(Clone, Debug, PartialEq)]
pub struct GaussianOutput {
    pub intensities: Vec<f64>,
    /// True if at least one integrated intensity is positive.
    pub found_signal: bool,
}

/// Tabulated Gaussian kernel and normalized trapezoidal convolution.
#[derive(Clone, Copy, Debug)]
pub struct GaussFilterAlgorithm {
    pub width: GaussianWidth,
    /// Spacing of the kernel lookup table, independently of input spacing.
    pub kernel_spacing: f64,
    /// Maximum one-sided kernel coefficients (default one million).
    pub max_coefficients: usize,
}
impl Default for GaussFilterAlgorithm {
    fn default() -> Self {
        Self {
            width: GaussianWidth::Absolute(0.8),
            kernel_spacing: 0.01,
            max_coefficients: 1_000_000,
        }
    }
}
impl GaussFilterAlgorithm {
    pub fn new(width: GaussianWidth, kernel_spacing: f64) -> Result<Self> {
        let algorithm = Self {
            width,
            kernel_spacing,
            ..Self::default()
        };
        algorithm.validate_options()?;
        Ok(algorithm)
    }

    fn validate_options(&self) -> Result<()> {
        let width = match self.width {
            GaussianWidth::Absolute(w) | GaussianWidth::Ppm(w) => w,
        };
        if !width.is_finite()
            || width <= 0.0
            || !self.kernel_spacing.is_finite()
            || self.kernel_spacing <= 0.0
            || self.max_coefficients == 0
        {
            return Err(Error::InvalidValue(
                "Gaussian width, kernel spacing, and coefficient limit must be finite and positive"
                    .into(),
            ));
        }
        Ok(())
    }

    fn coefficients(&self, width: f64) -> Result<Vec<f64>> {
        let sigma = width / 8.0;
        let count = (4.0 * sigma / self.kernel_spacing).ceil();
        if !sigma.is_finite()
            || sigma <= 0.0
            || !count.is_finite()
            || count < 0.0
            || count >= self.max_coefficients as f64
        {
            return Err(Error::InvalidValue(
                "Gaussian kernel exceeds coefficient limit or has an unrepresentable width".into(),
            ));
        }
        let variance_denominator = 2.0 * sigma * sigma;
        if !variance_denominator.is_finite() || variance_denominator <= 0.0 {
            return Err(Error::InvalidValue(
                "Gaussian variance is not representable".into(),
            ));
        }
        let scale = 1.0 / (sigma * (2.0 * std::f64::consts::PI).sqrt());
        if !scale.is_finite() {
            return Err(Error::InvalidValue(
                "Gaussian kernel amplitude overflows".into(),
            ));
        }
        let count = count as usize + 1;
        Ok((0..count)
            .map(|i| {
                let x = i as f64 * self.kernel_spacing;
                scale * (-(x * x) / variance_denominator).exp()
            })
            .collect())
    }

    /// Smooth parallel f64 arrays without changing their coordinates.
    ///
    /// OpenMS excludes intervals whose outside endpoint is exactly at the
    /// clipped integration boundary. This behavior is intentionally preserved.
    /// A zero/negative integrated numerator produces zero.
    pub fn filter(&self, positions: &[f64], intensities: &[f64]) -> Result<GaussianOutput> {
        self.validate_options()?;
        validate_signal(positions, intensities)?;
        if matches!(self.width, GaussianWidth::Ppm(_)) && positions.iter().any(|&x| x <= 0.0) {
            return Err(Error::InvalidValue(
                "ppm Gaussian smoothing needs positive m/z positions".into(),
            ));
        }
        let absolute = match self.width {
            GaussianWidth::Absolute(width) => Some(self.coefficients(width)?),
            _ => None,
        };
        let mut output = Vec::with_capacity(positions.len());
        for index in 0..positions.len() {
            let ppm_coefficients;
            let coefficients = if let Some(coefficients) = &absolute {
                coefficients.as_slice()
            } else {
                let GaussianWidth::Ppm(ppm) = self.width else {
                    unreachable!()
                };
                ppm_coefficients = self.coefficients((ppm / 1e6) * positions[index])?;
                &ppm_coefficients
            };
            let support = coefficients.len() as f64 * self.kernel_spacing;
            if !support.is_finite() {
                return Err(Error::InvalidValue("Gaussian support overflows".into()));
            }
            let start = (positions[index] - support).max(positions[0]);
            let end = (positions[index] + support).min(positions[positions.len() - 1]);
            let mut numerator = 0.0;
            let mut norm = 0.0;
            let coefficient = |distance: f64| -> Result<f64> {
                let left = (distance / self.kernel_spacing).floor() as usize;
                if left >= coefficients.len() {
                    return Err(Error::InvalidValue(
                        "Gaussian interpolation outside kernel".into(),
                    ));
                }
                let fraction =
                    (left as f64 * self.kernel_spacing - distance).abs() / self.kernel_spacing;
                Ok(if left + 1 < coefficients.len() {
                    (1.0 - fraction) * coefficients[left] + fraction * coefficients[left + 1]
                } else {
                    coefficients[left]
                })
            };
            let mut interval = |left: usize, right: usize| -> Result<()> {
                let a = coefficient((positions[index] - positions[left]).abs())?;
                let b = coefficient((positions[index] - positions[right]).abs())?;
                let half_span = (positions[right] - positions[left]).abs() / 2.0;
                norm += half_span * (a + b);
                numerator += half_span * (intensities[left] * a + intensities[right] * b);
                Ok(())
            };
            let mut cursor = index;
            while cursor > 0 && positions[cursor - 1] > start {
                interval(cursor - 1, cursor)?;
                cursor -= 1;
            }
            cursor = index;
            while cursor + 1 < positions.len() && positions[cursor + 1] < end {
                interval(cursor, cursor + 1)?;
                cursor += 1;
            }
            let value = if numerator > 0.0 {
                numerator / norm
            } else {
                0.0
            };
            if !numerator.is_finite() || !norm.is_finite() || !value.is_finite() {
                return Err(Error::InvalidValue("Gaussian integration overflows".into()));
            }
            output.push(value);
        }
        Ok(GaussianOutput {
            found_signal: output.iter().any(|&x| x > 0.0),
            intensities: output,
        })
    }
}

/// OpenMS Gaussian wrapper, including its all-zero-result preservation rule.
#[derive(Clone, Copy, Debug)]
pub struct GaussFilter {
    pub algorithm: GaussFilterAlgorithm,
}
impl Default for GaussFilter {
    fn default() -> Self {
        Self {
            algorithm: GaussFilterAlgorithm {
                width: GaussianWidth::Absolute(0.2),
                ..Default::default()
            },
        }
    }
}
impl GaussFilter {
    pub fn new(width: GaussianWidth) -> Result<Self> {
        Ok(Self {
            algorithm: GaussFilterAlgorithm::new(width, 0.01)?,
        })
    }

    /// Chromatograms support absolute widths in seconds; ppm is an error.
    pub fn filter_chromatogram(&self, chromatogram: &mut MSChromatogram) -> Result<()> {
        chromatogram.validate()?;
        if matches!(self.algorithm.width, GaussianWidth::Ppm(_)) {
            return Err(Error::InvalidValue(
                "ppm Gaussian smoothing is not defined for chromatograms".into(),
            ));
        }
        let positions: Vec<_> = chromatogram.peaks.iter().map(|p| p.rt).collect();
        let intensities: Vec<_> = chromatogram
            .peaks
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect();
        let output = self.algorithm.filter(&positions, &intensities)?;
        if output.found_signal || chromatogram.len() < 3 {
            let output = convert(output.intensities)?;
            for (peak, value) in chromatogram.peaks.iter_mut().zip(output) {
                peak.intensity = value;
            }
        }
        Ok(())
    }
}
impl SpectrumFilter for GaussFilter {
    fn filter_spectrum(&self, spectrum: &mut MSSpectrum) -> Result<()> {
        spectrum.validate()?;
        let positions: Vec<_> = spectrum.peaks.iter().map(|p| p.mz).collect();
        let intensities: Vec<_> = spectrum
            .peaks
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect();
        let output = self.algorithm.filter(&positions, &intensities)?;
        let output = if output.found_signal || spectrum.len() < 3 {
            Some(convert(output.intensities)?)
        } else {
            None
        };
        if let Some(output) = output {
            for (peak, value) in spectrum.peaks.iter_mut().zip(output) {
                peak.intensity = value;
            }
        }
        spectrum.spectrum_type = SpectrumType::Profile;
        Ok(())
    }

    /// Like OpenMS, smooth spectra AND chromatograms; any error is atomic.
    fn filter_experiment(&self, experiment: &mut MSExperiment) -> Result<()> {
        super::AcquisitionCopies::default().experiment(experiment)?;
        let mut output = experiment.clone();
        for spectrum in &mut output.spectra {
            self.filter_spectrum(spectrum)?;
        }
        for chromatogram in &mut output.chromatograms {
            self.filter_chromatogram(chromatogram)?;
        }
        *experiment = output;
        Ok(())
    }
}

fn validate_signal(positions: &[f64], intensities: &[f64]) -> Result<()> {
    if positions.len() != intensities.len() {
        return Err(Error::InvalidValue(
            "position and intensity lengths differ".into(),
        ));
    }
    if positions.iter().chain(intensities).any(|x| !x.is_finite()) {
        return Err(Error::InvalidValue("signal values must be finite".into()));
    }
    if positions.windows(2).any(|p| p[0] > p[1]) {
        return Err(Error::UnsortedData);
    }
    Ok(())
}
fn convert(values: Vec<f64>) -> Result<Vec<f32>> {
    values.into_iter().map(checked_intensity).collect()
}

/// Savitzky–Golay polynomial least-squares projection with asymmetric edges.
///
/// Coefficients are computed once with reorthogonalized QR on scaled indices.
/// This solves the same least-squares problem as the C++ SVD without squaring
/// its condition number through normal equations. No external solver is used.
#[derive(Clone, Debug)]
pub struct SavitzkyGolayFilter {
    frame_length: usize,
    polynomial_order: usize,
    coefficients: Vec<Vec<f64>>,
}
impl Default for SavitzkyGolayFilter {
    fn default() -> Self {
        Self::new(11, 4).expect("fixed default Savitzky–Golay parameters are valid")
    }
}
impl SavitzkyGolayFilter {
    /// Even frame lengths are incremented; degree must be below the odd length.
    /// Frames above 1023 or degrees above 32 are rejected to bound coefficient
    /// generation and avoid pretending ill-conditioned fits are reliable.
    pub fn new(frame_length: usize, polynomial_order: usize) -> Result<Self> {
        let frame_length = if frame_length % 2 == 0 {
            frame_length
                .checked_add(1)
                .ok_or_else(|| Error::InvalidValue("frame length overflows".into()))?
        } else {
            frame_length
        };
        if frame_length > 1023 || polynomial_order > 32 || polynomial_order >= frame_length {
            return Err(Error::InvalidValue(
                "Savitzky–Golay requires degree < odd frame length <= 1023 and degree <= 32".into(),
            ));
        }
        let columns = polynomial_order + 1;
        let scale = (frame_length / 2).max(1) as f64;
        let coordinates: Vec<_> = (0..frame_length)
            .map(|i| (i as f64 - (frame_length / 2) as f64) / scale)
            .collect();
        let mut q: Vec<Vec<f64>> = Vec::with_capacity(columns);
        for degree in 0..columns {
            let mut column: Vec<_> = coordinates.iter().map(|x| x.powi(degree as i32)).collect();
            // Reorthogonalization counters loss of orthogonality in Vandermonde bases.
            for _ in 0..2 {
                for basis in &q {
                    let projection: f64 = column.iter().zip(basis).map(|(a, b)| a * b).sum();
                    for (value, basis_value) in column.iter_mut().zip(basis) {
                        *value -= projection * basis_value;
                    }
                }
            }
            let norm = column.iter().map(|x| x * x).sum::<f64>().sqrt();
            if !norm.is_finite() || norm <= f64::EPSILON * 64.0 * (frame_length as f64).sqrt() {
                return Err(Error::InvalidValue(
                    "Savitzky–Golay fit is numerically rank deficient; reduce polynomial order"
                        .into(),
                ));
            }
            for value in &mut column {
                *value /= norm;
            }
            q.push(column);
        }
        let coefficients = (0..=frame_length / 2)
            .map(|evaluation| {
                (0..frame_length)
                    .map(|sample| {
                        q.iter()
                            .map(|basis| basis[evaluation] * basis[sample])
                            .sum()
                    })
                    .collect()
            })
            .collect();
        Ok(Self {
            frame_length,
            polynomial_order,
            coefficients,
        })
    }

    pub fn frame_length(&self) -> usize {
        self.frame_length
    }
    pub fn polynomial_order(&self) -> usize {
        self.polynomial_order
    }

    /// Smooth sample-index data; positions should describe uniform profile data.
    /// Fewer points than the full window are left unchanged; fitted negatives
    /// are clamped to zero, including at the asymmetric edge windows.
    pub fn filter(&self, positions: &[f64], intensities: &[f64]) -> Result<Vec<f64>> {
        validate_signal(positions, intensities)?;
        if self.frame_length > intensities.len() {
            return Ok(intensities.to_vec());
        }
        let middle = self.frame_length / 2;
        let mut output = Vec::with_capacity(intensities.len());
        for index in 0..intensities.len() {
            let start = index
                .saturating_sub(middle)
                .min(intensities.len() - self.frame_length);
            let evaluation = index - start;
            let reflected = evaluation > middle;
            let row = if reflected {
                self.frame_length - evaluation - 1
            } else {
                evaluation
            };
            let mut value = 0.0;
            for (j, intensity) in intensities[start..start + self.frame_length]
                .iter()
                .enumerate()
            {
                let coefficient_index = if reflected {
                    self.frame_length - j - 1
                } else {
                    j
                };
                value += intensity * self.coefficients[row][coefficient_index];
            }
            if !value.is_finite() {
                return Err(Error::InvalidValue(
                    "Savitzky–Golay convolution overflows".into(),
                ));
            }
            output.push(value.max(0.0));
        }
        Ok(output)
    }

    pub fn filter_chromatogram(&self, chromatogram: &mut MSChromatogram) -> Result<()> {
        chromatogram.validate()?;
        let positions: Vec<_> = chromatogram.peaks.iter().map(|p| p.rt).collect();
        let intensities: Vec<_> = chromatogram
            .peaks
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect();
        let output = convert(self.filter(&positions, &intensities)?)?;
        for (peak, intensity) in chromatogram.peaks.iter_mut().zip(output) {
            peak.intensity = intensity;
        }
        Ok(())
    }
}
impl SpectrumFilter for SavitzkyGolayFilter {
    fn filter_spectrum(&self, spectrum: &mut MSSpectrum) -> Result<()> {
        spectrum.validate()?;
        let positions: Vec<_> = spectrum.peaks.iter().map(|p| p.mz).collect();
        let intensities: Vec<_> = spectrum
            .peaks
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect();
        let output = convert(self.filter(&positions, &intensities)?)?;
        for (peak, intensity) in spectrum.peaks.iter_mut().zip(output) {
            peak.intensity = intensity;
        }
        Ok(())
    }

    /// Like OpenMS, smooth spectra AND chromatograms; any error is atomic.
    fn filter_experiment(&self, experiment: &mut MSExperiment) -> Result<()> {
        super::AcquisitionCopies::default().experiment(experiment)?;
        let mut output = experiment.clone();
        for spectrum in &mut output.spectra {
            self.filter_spectrum(spectrum)?;
        }
        for chromatogram in &mut output.chromatograms {
            self.filter_chromatogram(chromatogram)?;
        }
        *experiment = output;
        Ok(())
    }
}
