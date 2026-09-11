// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Sampled peak integration, baseline estimates and shape metrics.
//!
//! Ported from `PeakIntegrator.h` at OpenMS4-core revision
//! `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. Bounds include observed points
//! without interpolation. Even-count Simpson integration can use a neighboring
//! sample outside the bounds. Optional EMG fitting integrates the full fitted
//! span, including extrapolated points, as in source EMGPreProcess_.

use crate::analysis::emg::EmgGradientDescent;
use crate::error::{Error, Result};
use crate::kernel::{MSChromatogram, MSSpectrum};
use std::fmt;
use std::ops::Range;
use std::str::FromStr;

/// Algorithm used for both peak integration and background area.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum IntegrationMethod {
    /// Sum observed intensities, without weighting by coordinate spacing.
    #[default]
    IntensitySum,
    /// Piecewise linear integral; adjacent intensities are added in `f32`.
    Trapezoid,
    /// Nonuniform Simpson integral, with the source's even-count averaging.
    Simpson,
}

impl IntegrationMethod {
    pub const fn name(self) -> &'static str {
        match self {
            Self::IntensitySum => "intensity_sum",
            Self::Trapezoid => "trapezoid",
            Self::Simpson => "simpson",
        }
    }
}

impl FromStr for IntegrationMethod {
    type Err = Error;
    fn from_str(value: &str) -> Result<Self> {
        match value {
            "intensity_sum" => Ok(Self::IntensitySum),
            "trapezoid" => Ok(Self::Trapezoid),
            "simpson" => Ok(Self::Simpson),
            _ => Err(invalid("unknown peak integration method")),
        }
    }
}

impl fmt::Display for IntegrationMethod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// Baseline determined from the first and last selected intensities.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum BaselineType {
    /// Line connecting the sampled endpoints.
    #[default]
    BaseToBase,
    /// Constant baseline at the smaller endpoint intensity.
    VerticalDivisionMin,
    /// Constant baseline at the larger endpoint intensity.
    VerticalDivisionMax,
}

impl BaselineType {
    pub const fn name(self) -> &'static str {
        match self {
            Self::BaseToBase => "base_to_base",
            Self::VerticalDivisionMin => "vertical_division_min",
            Self::VerticalDivisionMax => "vertical_division_max",
        }
    }
}

impl FromStr for BaselineType {
    type Err = Error;
    fn from_str(value: &str) -> Result<Self> {
        match value {
            "base_to_base" => Ok(Self::BaseToBase),
            "vertical_division" | "vertical_division_min" => Ok(Self::VerticalDivisionMin),
            "vertical_division_max" => Ok(Self::VerticalDivisionMax),
            _ => Err(invalid("unknown peak baseline type")),
        }
    }
}

impl fmt::Display for BaselineType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PeakArea {
    pub area: f64,
    /// First strictly positive maximum; zero if all selected intensities are nonpositive.
    pub height: f64,
    /// Position of that maximum, or `(left + right) / 2` if none is positive.
    pub apex_pos: f64,
    /// Selected `[position, intensity]` samples, in input order.
    pub hull_points: Vec<[f64; 2]>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PeakBackground {
    pub area: f64,
    pub height: f64,
}

/// Sampled shape measures; threshold positions are not interpolated.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PeakShapeMetrics {
    pub width_at_5: f64,
    pub width_at_10: f64,
    pub width_at_50: f64,
    pub start_position_at_5: f64,
    pub start_position_at_10: f64,
    pub start_position_at_50: f64,
    pub end_position_at_5: f64,
    pub end_position_at_10: f64,
    pub end_position_at_50: f64,
    pub total_width: f64,
    /// Zero when the apex equals the sampled 5% start, as in the source.
    pub tailing_factor: f64,
    /// Zero when the apex equals the sampled 10% start, as in the source.
    pub asymmetry_factor: f64,
    /// Right minus left intensity, computed in `f32`; not divided by width.
    pub slope_of_baseline: f64,
    /// Baseline intensity difference divided by supplied height; zero for zero height.
    pub baseline_delta_2_height: f64,
    pub points_across_baseline: usize,
    pub points_across_half_height: usize,
}

/// Read-only integration of a sorted spectrum or chromatogram.
///
/// All methods validate finite input, bounds, ordering and parallel arrays. The
/// source's finite negative areas are retained. Nonfinite calculated results and
/// divisions with no defined source result return errors. No input is changed.
#[derive(Clone, Debug, PartialEq)]
pub struct PeakIntegrator {
    pub integration_method: IntegrationMethod,
    pub baseline_type: BaselineType,
    /// Optional EMG preprocessing. Fit the requested inclusive bounds once per
    /// operation, then use the full fitted span. Supplied shape height and apex
    /// remain unchanged. `None` uses the observed samples directly.
    pub emg: Option<EmgGradientDescent>,
    /// Maximum whole-input point count, checked before validation or allocation.
    /// Must be positive. Work is linear in this count; hull storage is linear in
    /// selected points when fitting is disabled. Also caps fitted output points;
    /// EMG computation additionally uses the fitter's evaluation/iteration limits.
    pub max_points: usize,
}

impl Default for PeakIntegrator {
    fn default() -> Self {
        Self {
            integration_method: IntegrationMethod::default(),
            baseline_type: BaselineType::default(),
            emg: None,
            max_points: 1_000_000,
        }
    }
}

impl PeakIntegrator {
    /// Integrate sampled m/z points in inclusive `[left, right]`.
    pub fn integrate_spectrum(
        &self,
        input: &MSSpectrum,
        left: f64,
        right: f64,
    ) -> Result<PeakArea> {
        self.integrate(Trace::Spectrum(input), left, right)
    }

    /// Integrate sampled retention times in inclusive `[left, right]` (seconds).
    pub fn integrate_chromatogram(
        &self,
        input: &MSChromatogram,
        left: f64,
        right: f64,
    ) -> Result<PeakArea> {
        self.integrate(Trace::Chromatogram(input), left, right)
    }

    pub fn estimate_background_spectrum(
        &self,
        input: &MSSpectrum,
        left: f64,
        right: f64,
        apex: f64,
    ) -> Result<PeakBackground> {
        self.background(Trace::Spectrum(input), left, right, apex)
    }

    pub fn estimate_background_chromatogram(
        &self,
        input: &MSChromatogram,
        left: f64,
        right: f64,
        apex: f64,
    ) -> Result<PeakBackground> {
        self.background(Trace::Chromatogram(input), left, right, apex)
    }

    pub fn calculate_shape_metrics_spectrum(
        &self,
        input: &MSSpectrum,
        left: f64,
        right: f64,
        height: f64,
        apex: f64,
    ) -> Result<PeakShapeMetrics> {
        self.shape(Trace::Spectrum(input), left, right, height, apex)
    }

    pub fn calculate_shape_metrics_chromatogram(
        &self,
        input: &MSChromatogram,
        left: f64,
        right: f64,
        height: f64,
        apex: f64,
    ) -> Result<PeakShapeMetrics> {
        self.shape(Trace::Chromatogram(input), left, right, height, apex)
    }

    fn range(&self, trace: Trace<'_>, left: f64, right: f64) -> Result<Range<usize>> {
        if self.max_points == 0 || trace.len() > self.max_points {
            return Err(invalid("peak integration point limit exceeded or zero"));
        }
        if !left.is_finite() || !right.is_finite() || left > right {
            return Err(invalid("peak bounds must be finite and ordered"));
        }
        trace.validate()?;
        for i in 1..trace.len() {
            if trace.point(i - 1).0 > trace.point(i).0 {
                return Err(Error::UnsortedData);
            }
        }
        Ok(trace.bound(left, false)..trace.bound(right, true))
    }

    fn preprocess<'a>(&self, trace: Trace<'a>, left: f64, right: f64) -> Result<PreparedTrace<'a>> {
        // Bound the whole original trace before fitting allocates its working vectors.
        self.range(trace, left, right)?;
        match (&self.emg, trace) {
            (None, trace) => Ok(PreparedTrace::Borrowed(trace)),
            (Some(config), trace) => {
                let fitter = EmgGradientDescent {
                    max_points: config.max_points.min(self.max_points),
                    ..config.clone()
                };
                // Some(0.0) is a literal boundary, unlike the C++ zero sentinel.
                // The fitter returns typed parameter diagnostics, never an unaligned array.
                match trace {
                    Trace::Spectrum(input) => Ok(PreparedTrace::Spectrum(Box::new(
                        fitter
                            .fit_spectrum(input, Some(left), Some(right))?
                            .spectrum,
                    ))),
                    Trace::Chromatogram(input) => Ok(PreparedTrace::Chromatogram(Box::new(
                        fitter
                            .fit_chromatogram(input, Some(left), Some(right))?
                            .chromatogram,
                    ))),
                }
            }
        }
    }

    fn integrate(&self, trace: Trace<'_>, left: f64, right: f64) -> Result<PeakArea> {
        let prepared = self.preprocess(trace, left, right)?;
        let (trace, left, right) = prepared.view(left, right)?;
        let range = self.range(trace, left, right)?;
        let mut output = PeakArea {
            apex_pos: (left + right) / 2.0,
            hull_points: Vec::with_capacity(range.len()),
            ..Default::default()
        };
        for i in range.clone() {
            let (position, intensity) = trace.point(i);
            let intensity = f64::from(intensity);
            output.hull_points.push([position, intensity]);
            if output.height < intensity {
                output.height = intensity;
                output.apex_pos = position;
            }
        }
        output.area = match self.integration_method {
            IntegrationMethod::IntensitySum => {
                range.clone().map(|i| f64::from(trace.point(i).1)).sum()
            }
            IntegrationMethod::Trapezoid => trapezoid(trace, range),
            IntegrationMethod::Simpson if range.len() <= 2 => trapezoid(trace, range),
            IntegrationMethod::Simpson if range.len() % 2 == 1 => simpson(trace, range)?,
            IntegrationMethod::Simpson => {
                let mut areas = [None; 4];
                areas[0] = Some(simpson(trace, range.start..range.end - 1)?);
                areas[1] = Some(simpson(trace, range.start + 1..range.end)?);
                if range.start > 0 {
                    areas[2] = Some(simpson(trace, range.start - 1..range.end)?);
                }
                if range.end < trace.len() {
                    areas[3] = Some(simpson(trace, range.start..range.end + 1)?);
                }
                let mut sum = 0.0;
                let mut count = 0;
                for area in areas.into_iter().flatten() {
                    // Preserve the source's -1 sentinel collision for computed areas too.
                    if area != -1.0 {
                        sum += area;
                        count += 1;
                    }
                }
                if count == 0 {
                    return Err(invalid(
                        "all even-count Simpson subareas equal the source -1 sentinel",
                    ));
                }
                sum / f64::from(count)
            }
        };
        finite(output.area, "integrated area")?;
        finite(output.apex_pos, "peak apex")?;
        Ok(output)
    }

    fn background(
        &self,
        trace: Trace<'_>,
        left: f64,
        right: f64,
        apex: f64,
    ) -> Result<PeakBackground> {
        validate_apex(left, right, apex)?;
        let prepared = self.preprocess(trace, left, right)?;
        let (trace, left, right) = prepared.view(left, right)?;
        let range = self.range(trace, left, right)?;
        if range.is_empty() {
            return Err(invalid("background estimation requires selected samples"));
        }
        let (left_pos, left_int) = trace.point(range.start);
        let (right_pos, right_int) = trace.point(range.end - 1);
        let (left_int, right_int) = (f64::from(left_int), f64::from(right_int));
        let width = finite(right_pos - left_pos, "baseline width")?;
        let minimum = left_int.min(right_int);
        let count = range.len() as f64;
        let summed = self.integration_method == IntegrationMethod::IntensitySum;
        let output = match self.baseline_type {
            BaselineType::BaseToBase => {
                if width <= 0.0 {
                    return Err(invalid(
                        "base-to-base baseline requires distinct endpoint positions",
                    ));
                }
                let delta = right_int - left_int;
                let min_pos = if right_int <= left_int {
                    right_pos
                } else {
                    left_pos
                };
                let height = minimum + delta.abs() * (min_pos - apex).abs() / width;
                let area = if summed {
                    let pos_sum: f64 = range.map(|i| trace.point(i).0).sum();
                    (pos_sum - count * left_pos) * (delta / width) + count * left_int
                } else {
                    width * (minimum + 0.5 * delta.abs())
                };
                PeakBackground { area, height }
            }
            BaselineType::VerticalDivisionMin | BaselineType::VerticalDivisionMax => {
                let height = if self.baseline_type == BaselineType::VerticalDivisionMin {
                    minimum
                } else {
                    left_int.max(right_int)
                };
                PeakBackground {
                    area: height * if summed { count } else { width },
                    height,
                }
            }
        };
        finite(output.area, "background area")?;
        finite(output.height, "background height")?;
        Ok(output)
    }

    fn shape(
        &self,
        trace: Trace<'_>,
        left: f64,
        right: f64,
        height: f64,
        apex: f64,
    ) -> Result<PeakShapeMetrics> {
        self.range(trace, left, right)?;
        validate_apex(left, right, apex)?;
        if !height.is_finite() || height < 0.0 {
            return Err(invalid("shape height must be finite and nonnegative"));
        }
        if trace.len() == 0 {
            // Source returns empty shape metrics before invoking EMG preprocessing.
            return Ok(PeakShapeMetrics::default());
        }
        let prepared = self.preprocess(trace, left, right)?;
        let (trace, left, right) = prepared.view(left, right)?;
        let range = self.range(trace, left, right)?;
        let apex_index = trace.bound(apex, false);
        if range.is_empty() || !range.contains(&apex_index) {
            return Err(invalid(
                "shape apex must identify a selected sample at or above its position",
            ));
        }
        let starts = [0.05, 0.1, 0.5].map(|fraction| {
            threshold_position(trace, range.start..apex_index, height * fraction, true)
        });
        let ends = [0.05, 0.1, 0.5].map(|fraction| {
            threshold_position(trace, apex_index..range.end, height * fraction, false)
        });
        // Like the source, baseline slope is a float intensity subtraction, not dy/dx.
        let slope = f64::from(trace.point(range.end - 1).1 - trace.point(range.start).1);
        let mut output = PeakShapeMetrics {
            width_at_5: ends[0] - starts[0],
            width_at_10: ends[1] - starts[1],
            width_at_50: ends[2] - starts[2],
            start_position_at_5: starts[0],
            start_position_at_10: starts[1],
            start_position_at_50: starts[2],
            end_position_at_5: ends[0],
            end_position_at_10: ends[1],
            end_position_at_50: ends[2],
            total_width: trace.point(range.end - 1).0 - trace.point(range.start).0,
            slope_of_baseline: slope,
            points_across_baseline: range.len(),
            points_across_half_height: range
                .filter(|&i| f64::from(trace.point(i).1) >= 0.5 * height)
                .count(),
            ..Default::default()
        };
        if height != 0.0 {
            output.baseline_delta_2_height = slope / height;
        }
        if starts[0] != apex {
            let denominator = finite(2.0 * (apex - starts[0]), "tailing factor denominator")?;
            output.tailing_factor = output.width_at_5 / denominator;
        }
        if starts[1] != apex {
            output.asymmetry_factor = (ends[1] - apex) / (apex - starts[1]);
        }
        for value in [
            output.width_at_5,
            output.width_at_10,
            output.width_at_50,
            output.total_width,
            output.slope_of_baseline,
            output.baseline_delta_2_height,
            output.tailing_factor,
            output.asymmetry_factor,
        ] {
            finite(value, "peak shape metric")?;
        }
        Ok(output)
    }
}

// The fitted container is owned only for this operation; no-fit paths borrow input.
enum PreparedTrace<'a> {
    Borrowed(Trace<'a>),
    Spectrum(Box<MSSpectrum>),
    Chromatogram(Box<MSChromatogram>),
}

impl PreparedTrace<'_> {
    fn view(&self, left: f64, right: f64) -> Result<(Trace<'_>, f64, f64)> {
        let trace = match self {
            Self::Borrowed(trace) => return Ok((*trace, left, right)),
            Self::Spectrum(input) => Trace::Spectrum(input),
            Self::Chromatogram(input) => Trace::Chromatogram(input),
        };
        if trace.len() == 0 {
            return Err(invalid("EMG preprocessing returned no fitted points"));
        }
        Ok((trace, trace.point(0).0, trace.point(trace.len() - 1).0))
    }
}

// Borrow either native container without copying peaks or inventing a public trace framework.
#[derive(Clone, Copy)]
enum Trace<'a> {
    Spectrum(&'a MSSpectrum),
    Chromatogram(&'a MSChromatogram),
}

impl Trace<'_> {
    fn len(self) -> usize {
        match self {
            Self::Spectrum(input) => input.peaks.len(),
            Self::Chromatogram(input) => input.peaks.len(),
        }
    }
    fn point(self, index: usize) -> (f64, f32) {
        match self {
            Self::Spectrum(input) => (input.peaks[index].mz, input.peaks[index].intensity),
            Self::Chromatogram(input) => (input.peaks[index].rt, input.peaks[index].intensity),
        }
    }
    fn validate(self) -> Result<()> {
        match self {
            Self::Spectrum(input) => input.validate(),
            Self::Chromatogram(input) => input.validate(),
        }
    }
    fn bound(self, position: f64, upper: bool) -> usize {
        // Select the same inclusive sample range as source PosBegin/PosEnd.
        match self {
            Self::Spectrum(input) => input
                .peaks
                .partition_point(|p| p.mz < position || (upper && p.mz == position)),
            Self::Chromatogram(input) => input
                .peaks
                .partition_point(|p| p.rt < position || (upper && p.rt == position)),
        }
    }
}

fn trapezoid(trace: Trace<'_>, range: Range<usize>) -> f64 {
    let mut area = 0.0;
    for i in range.start..range.end.saturating_sub(1) {
        let (x, y) = trace.point(i);
        let (next_x, next_y) = trace.point(i + 1);
        area += (next_x - x) * (f64::from(y + next_y) / 2.0);
    }
    area
}

fn simpson(trace: Trace<'_>, range: Range<usize>) -> Result<f64> {
    let mut area = 0.0;
    for mid in (range.start + 1..range.end - 1).step_by(2) {
        let (left_x, left_y) = trace.point(mid - 1);
        let (x, y) = trace.point(mid);
        let (right_x, right_y) = trace.point(mid + 1);
        let h = x - left_x;
        let k = right_x - x;
        if h <= 0.0 || k <= 0.0 {
            return Err(invalid(
                "Simpson integration requires strictly increasing used sample positions",
            ));
        }
        area += (1.0 / 6.0)
            * (h + k)
            * ((2.0 - k / h) * f64::from(left_y)
                + ((h + k).powi(2) / (h * k)) * f64::from(y)
                + (2.0 - h / k) * f64::from(right_y));
    }
    finite(area, "Simpson subarea")
}

fn threshold_position(
    trace: Trace<'_>,
    range: Range<usize>,
    threshold: f64,
    left_half: bool,
) -> f64 {
    if range.is_empty() {
        return trace.point(range.start).0;
    }
    let mut closest = if left_half {
        range.start
    } else {
        range.end - 1
    };
    if left_half {
        for i in range {
            if f64::from(trace.point(i).1) > threshold {
                break;
            }
            closest = i;
        }
    } else {
        for i in range.rev() {
            if f64::from(trace.point(i).1) > threshold {
                break;
            }
            closest = i;
        }
    }
    trace.point(closest).0
}

fn validate_apex(left: f64, right: f64, apex: f64) -> Result<()> {
    if !apex.is_finite() || apex < left || apex > right {
        return Err(invalid(
            "peak apex must be finite and within the requested bounds",
        ));
    }
    Ok(())
}

fn finite(value: f64, label: &str) -> Result<f64> {
    if !value.is_finite() {
        return Err(invalid(&format!("{label} is not finite")));
    }
    Ok(value)
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
