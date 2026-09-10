// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Flat-structuring-element morphology with clipped boundary windows.
//! Min/max operations use a linear-time monotonic deque; composition and signed
//! subtraction follow OpenMS, including `bothat = input - closing`.

use super::{SpectrumFilter, checked_intensity};
use crate::kernel::{MSChromatogram, MSSpectrum, SpectrumType};
use crate::{Error, Result};
use std::collections::VecDeque;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MorphologicalMethod {
    Identity,
    Erosion,
    Dilation,
    Opening,
    Closing,
    Gradient,
    #[default]
    TopHat,
    /// OpenMS's signed bottom-hat: input minus closing, generally nonpositive.
    BottomHat,
    /// Direct-window reference implementation with erosion semantics.
    ErosionSimple,
    /// Direct-window reference implementation with dilation semantics.
    DilationSimple,
}

/// Width is rounded up to an odd number of samples, matching the C++ wrapper.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum StructuringElement {
    DataPoints(usize),
    /// Spectrum width in Th, converted using global average sample spacing.
    Thomson(f64),
    /// Chromatogram extension: seconds, using global average sample spacing.
    Seconds(f64),
}
#[derive(Clone, Copy, Debug)]
pub struct MorphologicalFilter {
    pub method: MorphologicalMethod,
    pub structuring_element: StructuringElement,
}
impl Default for MorphologicalFilter {
    fn default() -> Self {
        Self {
            method: MorphologicalMethod::TopHat,
            structuring_element: StructuringElement::Thomson(3.0),
        }
    }
}
impl MorphologicalFilter {
    pub fn new(
        method: MorphologicalMethod,
        structuring_element: StructuringElement,
    ) -> Result<Self> {
        let filter = Self {
            method,
            structuring_element,
        };
        filter.validate_options()?;
        Ok(filter)
    }
    fn validate_options(&self) -> Result<()> {
        match self.structuring_element {
            StructuringElement::DataPoints(0) => Err(Error::InvalidValue(
                "structuring element must contain at least one point".into(),
            )),
            StructuringElement::Thomson(width) | StructuringElement::Seconds(width)
                if !width.is_finite() || width <= 0.0 =>
            {
                Err(Error::InvalidValue(
                    "structuring element width must be finite and positive".into(),
                ))
            }
            _ => Ok(()),
        }
    }

    /// Return the rounded odd window length. Width conversion is
    /// ceil(width * (N - 1) / (last_position - first_position)), then odd rounding.
    pub fn effective_window(&self, positions: &[f64]) -> Result<usize> {
        self.validate_options()?;
        validate_positions(positions)?;
        let points = match self.structuring_element {
            StructuringElement::DataPoints(points) => points,
            StructuringElement::Thomson(width) | StructuringElement::Seconds(width) => {
                if positions.len() < 2 {
                    return Ok(1);
                }
                let span = positions[positions.len() - 1] - positions[0];
                if !span.is_finite() || span <= 0.0 {
                    return Err(Error::InvalidValue(
                        "coordinate-based morphology needs a positive finite signal span".into(),
                    ));
                }
                let points = (width * (positions.len() - 1) as f64 / span).ceil();
                if !points.is_finite() || points >= usize::MAX as f64 {
                    return Err(Error::InvalidValue(
                        "structuring element length overflows".into(),
                    ));
                }
                points as usize
            }
        };
        if points % 2 == 0 {
            points
                .checked_add(1)
                .ok_or_else(|| Error::InvalidValue("structuring element length overflows".into()))
        } else {
            Ok(points)
        }
    }

    /// Morphology directly on intensities; requires DataPoints units.
    /// Unlike the spectrum convenience wrapper, a singleton TopHat is zero.
    pub fn filter_range(&self, intensities: &[f32]) -> Result<Vec<f32>> {
        self.validate_options()?;
        let StructuringElement::DataPoints(points) = self.structuring_element else {
            return Err(Error::InvalidValue(
                "intensity-only morphology requires DataPoints units".into(),
            ));
        };
        let points = if points % 2 == 0 { points + 1 } else { points };
        self.apply(intensities, points)
    }

    fn apply(&self, input: &[f32], points: usize) -> Result<Vec<f32>> {
        if input.iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidValue(
                "morphology intensities must be finite".into(),
            ));
        }
        use MorphologicalMethod::*;
        let erosion = |values: &[f32]| extrema(values, points, false);
        let dilation = |values: &[f32]| extrema(values, points, true);
        match self.method {
            Identity => Ok(input.to_vec()),
            Erosion => Ok(erosion(input)),
            Dilation => Ok(dilation(input)),
            ErosionSimple => Ok(extrema_simple(input, points, false)),
            DilationSimple => Ok(extrema_simple(input, points, true)),
            Opening => Ok(dilation(&erosion(input))),
            Closing => Ok(erosion(&dilation(input))),
            Gradient => subtract(&dilation(input), &erosion(input)),
            TopHat => subtract(input, &dilation(&erosion(input))),
            BottomHat => subtract(input, &erosion(&dilation(input))),
        }
    }

    /// Apply the same morphology to chromatograms. Seconds replaces Thomson;
    /// this convenience API extends the C++ class's spectrum-only interface.
    pub fn filter_chromatogram(&self, chromatogram: &mut MSChromatogram) -> Result<()> {
        self.validate_options()?;
        chromatogram.validate()?;
        if matches!(self.structuring_element, StructuringElement::Thomson(_)) {
            return Err(Error::InvalidValue(
                "chromatogram morphology uses Seconds or DataPoints".into(),
            ));
        }
        let positions: Vec<_> = chromatogram.peaks.iter().map(|p| p.rt).collect();
        let points = self.effective_window(&positions)?;
        if chromatogram.len() <= 1 {
            return Ok(());
        }
        let values: Vec<_> = chromatogram.peaks.iter().map(|p| p.intensity).collect();
        let output = self.apply(&values, points)?;
        for (peak, value) in chromatogram.peaks.iter_mut().zip(output) {
            peak.intensity = value;
        }
        Ok(())
    }
}
impl SpectrumFilter for MorphologicalFilter {
    fn filter_spectrum(&self, spectrum: &mut MSSpectrum) -> Result<()> {
        self.validate_options()?;
        spectrum.validate()?;
        if matches!(self.structuring_element, StructuringElement::Seconds(_)) {
            return Err(Error::InvalidValue(
                "spectrum morphology uses Thomson or DataPoints".into(),
            ));
        }
        let positions: Vec<_> = spectrum.peaks.iter().map(|p| p.mz).collect();
        let points = self.effective_window(&positions)?;
        if spectrum.len() > 1 {
            let values: Vec<_> = spectrum.peaks.iter().map(|p| p.intensity).collect();
            let output = self.apply(&values, points)?;
            for (peak, value) in spectrum.peaks.iter_mut().zip(output) {
                peak.intensity = value;
            }
        }
        spectrum.spectrum_type = SpectrumType::Profile;
        Ok(())
    }
}
fn validate_positions(positions: &[f64]) -> Result<()> {
    if positions.iter().any(|x| !x.is_finite()) {
        return Err(Error::InvalidValue("coordinates must be finite".into()));
    }
    if positions.windows(2).any(|p| p[0] > p[1]) {
        return Err(Error::UnsortedData);
    }
    Ok(())
}
fn subtract(left: &[f32], right: &[f32]) -> Result<Vec<f32>> {
    left.iter()
        .zip(right)
        .map(|(&a, &b)| checked_intensity(f64::from(a) - f64::from(b)))
        .collect()
}
fn extrema(values: &[f32], points: usize, maximum: bool) -> Vec<f32> {
    let radius = points / 2;
    let mut queue: VecDeque<usize> = VecDeque::new();
    let mut next = 0;
    let mut output = Vec::with_capacity(values.len());
    for center in 0..values.len() {
        let start = center.saturating_sub(radius);
        let end = center
            .saturating_add(radius)
            .saturating_add(1)
            .min(values.len());
        while next < end {
            while queue.back().is_some_and(|&i| {
                if maximum {
                    values[i] < values[next]
                } else {
                    values[i] > values[next]
                }
            }) {
                queue.pop_back();
            }
            queue.push_back(next);
            next += 1;
        }
        while queue.front().is_some_and(|&i| i < start) {
            queue.pop_front();
        }
        output.push(
            values[*queue
                .front()
                .expect("each nonempty window contains its center")],
        );
    }
    output
}
fn extrema_simple(values: &[f32], points: usize, maximum: bool) -> Vec<f32> {
    let radius = points / 2;
    (0..values.len())
        .map(|center| {
            let start = center.saturating_sub(radius);
            let end = center
                .saturating_add(radius)
                .saturating_add(1)
                .min(values.len());
            values[start + 1..end].iter().fold(values[start], |a, &b| {
                if maximum { a.max(b) } else { a.min(b) }
            })
        })
        .collect()
}
