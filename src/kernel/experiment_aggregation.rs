// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Window reductions from MSExperiment.h. Data remain borrowed and outputs are
//! staged locally; no acquisition or identification metadata are cloned.

use super::{
    ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, NumericRange, Peak1D, finite,
};
use crate::{Error, Result};
use std::{mem::size_of, ops::Range, str::FromStr};

/// Inclusive m/z and retention-time window, with RT measured in seconds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MzRtRegion {
    pub mz: NumericRange,
    pub rt: NumericRange,
}

impl MzRtRegion {
    pub fn new(min_mz: f64, max_mz: f64, min_rt: f64, max_rt: f64) -> Result<Self> {
        let region = Self {
            mz: NumericRange {
                min: min_mz,
                max: max_mz,
            },
            rt: NumericRange {
                min: min_rt,
                max: max_rt,
            },
        };
        region.validate()?;
        Ok(region)
    }

    fn validate(&self) -> Result<()> {
        for (range, label) in [(self.mz, "m/z region"), (self.rt, "RT region")] {
            finite(range.min, label)?;
            finite(range.max, label)?;
            if range.min > range.max {
                return Err(invalid("region minimum exceeds maximum"));
            }
        }
        Ok(())
    }
}

/// Matrix-wrapper reducers. Sum and mean accumulate widened intensities in f64,
/// unlike the f32 default of [`MSExperiment::aggregate`]. Empty slices yield 0.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MzAggregation {
    Sum,
    Max,
    Min,
    Mean,
}

impl FromStr for MzAggregation {
    type Err = Error;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "sum" => Ok(Self::Sum),
            "max" => Ok(Self::Max),
            "min" => Ok(Self::Min),
            "mean" => Ok(Self::Mean),
            _ => Err(invalid(
                "unknown m/z aggregation (expected sum, max, min or mean)",
            )),
        }
    }
}

impl MzAggregation {
    /// Reduce one slice with the source matrix-wrapper arithmetic.
    /// Standalone calls use the default work bound; callbacks supplied to a
    /// limited experiment operation remain responsible for their own work.
    pub fn reduce(self, peaks: &[Peak1D]) -> Result<f64> {
        let mut work = Work::new(AggregationLimits::default());
        work.consume(peaks.len())?;
        work.consume(peaks.len())?;
        for peak in peaks {
            finite(f64::from(peak.intensity), "aggregation intensity")?;
        }
        self.reduce_validated(peaks)
    }

    fn reduce_validated(self, peaks: &[Peak1D]) -> Result<f64> {
        let Some(first) = peaks.first() else {
            return Ok(0.0);
        };
        let value = match self {
            Self::Sum | Self::Mean => {
                let sum = peaks
                    .iter()
                    .fold(0.0, |sum, peak| sum + f64::from(peak.intensity));
                if self == Self::Mean {
                    sum / peaks.len() as f64
                } else {
                    sum
                }
            }
            Self::Max => f64::from(peaks.iter().fold(first.intensity, |best, peak| {
                if best < peak.intensity {
                    peak.intensity
                } else {
                    best
                }
            })),
            Self::Min => f64::from(peaks.iter().fold(first.intensity, |best, peak| {
                if peak.intensity < best {
                    peak.intensity
                } else {
                    best
                }
            })),
        };
        finite(value, "aggregation result")?;
        Ok(value)
    }
}

/// Per-operation limits shared by all windows, output rows and reductions.
/// Byte accounting includes conservative allocation overhead for every vector.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AggregationLimits {
    pub max_work: usize,
    pub max_bytes: usize,
    pub max_output_points: usize,
}

impl Default for AggregationLimits {
    fn default() -> Self {
        Self {
            max_work: 50_000_000,
            max_bytes: 256 * 1024 * 1024,
            max_output_points: 10_000_000,
        }
    }
}

impl MSExperiment {
    /// Sum intensities in each inclusive region, one f64 result per matching
    /// scan, using source f32 accumulation before widening. MS level is an
    /// exact filter: zero does not mean all levels. No regions or no matching
    /// level yields an empty outer vector; a window without RT scans has an
    /// empty row. Selected scans must be RT-sorted and visited peaks m/z-sorted.
    pub fn aggregate(&self, regions: &[MzRtRegion], ms_level: u32) -> Result<Vec<Vec<f64>>> {
        self.aggregate_with(regions, ms_level, sum_f32)
    }

    /// Apply a reducer to every selected peak slice, including empty slices.
    /// Calls occur in region order, then scan order. The source parallel API
    /// does not promise callback order. Callback side effects cannot be rolled
    /// back or bounded; the library bounds supplied slices and its own work.
    pub fn aggregate_with<F>(
        &self,
        regions: &[MzRtRegion],
        ms_level: u32,
        reducer: F,
    ) -> Result<Vec<Vec<f64>>>
    where
        F: FnMut(&[Peak1D]) -> Result<f64>,
    {
        self.aggregate_with_limits(regions, ms_level, reducer, &AggregationLimits::default())
    }

    /// Custom reduction with cumulative caller-selected resource limits.
    pub fn aggregate_with_limits<F>(
        &self,
        regions: &[MzRtRegion],
        ms_level: u32,
        reducer: F,
        limits: &AggregationLimits,
    ) -> Result<Vec<Vec<f64>>>
    where
        F: FnMut(&[Peak1D]) -> Result<f64>,
    {
        let mut work = Work::new(*limits);
        reduce(
            self,
            regions,
            ms_level,
            reducer,
            |_, value| Ok(value),
            &mut work,
        )
    }

    /// Matrix rows are `[min_mz, max_mz, min_rt, max_rt]`. Four-column shape is
    /// enforced by the Rust type. Sum/mean use f64 source-wrapper arithmetic.
    pub fn aggregate_from_matrix(
        &self,
        rows: &[[f64; 4]],
        ms_level: u32,
        aggregation: MzAggregation,
    ) -> Result<Vec<Vec<f64>>> {
        let mut work = Work::new(AggregationLimits::default());
        let regions = matrix_regions(rows, &mut work)?;
        reduce(
            self,
            &regions,
            ms_level,
            |peaks| aggregation.reduce_validated(peaks),
            |_, value| Ok(value),
            &mut work,
        )
    }

    /// Extract chromatograms with source f32 intensity sums and full-f64 RTs.
    /// Each product m/z is `(region.mz.min + region.mz.max) / 2`, including
    /// chromatograms whose RT window is empty. Other metadata have defaults.
    pub fn extract_xics(
        &self,
        regions: &[MzRtRegion],
        ms_level: u32,
    ) -> Result<Vec<MSChromatogram>> {
        self.extract_xics_with(regions, ms_level, sum_f32)
    }

    /// Custom XIC reduction. Finite f64 results narrow once to finite f32.
    /// Callback ordering and side-effect rules match [`Self::aggregate_with`].
    pub fn extract_xics_with<F>(
        &self,
        regions: &[MzRtRegion],
        ms_level: u32,
        reducer: F,
    ) -> Result<Vec<MSChromatogram>>
    where
        F: FnMut(&[Peak1D]) -> Result<f64>,
    {
        self.extract_xics_with_limits(regions, ms_level, reducer, &AggregationLimits::default())
    }

    /// Custom XIC reduction with limits shared across the entire operation.
    pub fn extract_xics_with_limits<F>(
        &self,
        regions: &[MzRtRegion],
        ms_level: u32,
        reducer: F,
        limits: &AggregationLimits,
    ) -> Result<Vec<MSChromatogram>>
    where
        F: FnMut(&[Peak1D]) -> Result<f64>,
    {
        xics(self, regions, ms_level, reducer, &mut Work::new(*limits))
    }

    /// XIC matrix wrapper, with f64 sum/mean followed by f32 peak storage.
    pub fn extract_xics_from_matrix(
        &self,
        rows: &[[f64; 4]],
        ms_level: u32,
        aggregation: MzAggregation,
    ) -> Result<Vec<MSChromatogram>> {
        let mut work = Work::new(AggregationLimits::default());
        let regions = matrix_regions(rows, &mut work)?;
        xics(
            self,
            &regions,
            ms_level,
            |peaks| aggregation.reduce_validated(peaks),
            &mut work,
        )
    }
}

fn matrix_regions(rows: &[[f64; 4]], work: &mut Work) -> Result<Vec<MzRtRegion>> {
    work.consume(rows.len())?;
    let mut result = work.vector(rows.len())?;
    for &[min_mz, max_mz, min_rt, max_rt] in rows {
        result.push(MzRtRegion::new(min_mz, max_mz, min_rt, max_rt)?);
    }
    Ok(result)
}

struct Work {
    remaining: usize,
    bytes: usize,
    points: usize,
}

impl Work {
    fn new(limits: AggregationLimits) -> Self {
        Self {
            remaining: limits.max_work,
            bytes: limits.max_bytes,
            points: limits.max_output_points,
        }
    }

    fn consume(&mut self, count: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(count)
            .ok_or_else(|| invalid("aggregation work limit exceeded"))?;
        Ok(())
    }

    fn vector<T>(&mut self, count: usize) -> Result<Vec<T>> {
        if count == 0 {
            return Ok(Vec::new());
        }
        let bytes = count
            .checked_mul(size_of::<T>())
            .and_then(|n| n.checked_add(64))
            .ok_or_else(|| invalid("aggregation allocation size overflow"))?;
        self.bytes = self
            .bytes
            .checked_sub(bytes)
            .ok_or_else(|| invalid("aggregation byte limit exceeded"))?;
        let mut result = Vec::new();
        result
            .try_reserve_exact(count)
            .map_err(|_| invalid("aggregation allocation failed"))?;
        Ok(result)
    }

    fn search(&mut self, count: usize) -> Result<()> {
        self.consume(if count == 0 {
            0
        } else {
            (usize::BITS - count.leading_zeros()) as usize + 1
        })
    }
}

/// A scan view and interval plan avoid both deep copies and a dense mapping of
/// all scans to all regions. Every output cell is counted before reducing.
struct Plan<'a> {
    spectra: Vec<&'a MSSpectrum>,
    intervals: Vec<Range<usize>>,
    validated: Vec<bool>,
}

fn plan<'a>(
    experiment: &'a MSExperiment,
    regions: &[MzRtRegion],
    ms_level: u32,
    work: &mut Work,
) -> Result<Plan<'a>> {
    let mut result = Plan {
        spectra: Vec::new(),
        intervals: Vec::new(),
        validated: Vec::new(),
    };
    if regions.is_empty() {
        return Ok(result);
    }
    work.consume(experiment.spectra.len())?;
    let count = experiment
        .spectra
        .iter()
        .filter(|s| s.ms_level == ms_level)
        .count();
    if count == 0 {
        return Ok(result);
    }
    work.consume(experiment.spectra.len())?;
    result.spectra = work.vector(count)?;
    let mut previous = None;
    for spectrum in &experiment.spectra {
        if spectrum.ms_level != ms_level {
            continue;
        }
        finite(spectrum.rt, "aggregation retention time")?;
        if previous.is_some_and(|rt| rt > spectrum.rt) {
            return Err(Error::UnsortedData);
        }
        previous = Some(spectrum.rt);
        result.spectra.push(spectrum);
    }
    work.consume(regions.len())?;
    result.intervals = work.vector(regions.len())?;
    for region in regions {
        region.validate()?;
        work.search(count)?;
        work.search(count)?;
        let first = result.spectra.partition_point(|s| s.rt < region.rt.min);
        let end = result.spectra.partition_point(|s| s.rt <= region.rt.max);
        work.points = work
            .points
            .checked_sub(end - first)
            .ok_or_else(|| invalid("aggregation output point limit exceeded"))?;
        result.intervals.push(first..end);
    }
    work.consume(count)?;
    result.validated = work.vector(count)?;
    result.validated.resize(count, false);
    Ok(result)
}

fn reduce<T, F, G>(
    experiment: &MSExperiment,
    regions: &[MzRtRegion],
    ms_level: u32,
    reducer: F,
    convert: G,
    work: &mut Work,
) -> Result<Vec<Vec<T>>>
where
    F: FnMut(&[Peak1D]) -> Result<f64>,
    G: FnMut(&MSSpectrum, f64) -> Result<T>,
{
    let plan = plan(experiment, regions, ms_level, work)?;
    reduce_plan(plan, regions, reducer, convert, work)
}

fn reduce_plan<T, F, G>(
    mut plan: Plan<'_>,
    regions: &[MzRtRegion],
    mut reducer: F,
    mut convert: G,
    work: &mut Work,
) -> Result<Vec<Vec<T>>>
where
    F: FnMut(&[Peak1D]) -> Result<f64>,
    G: FnMut(&MSSpectrum, f64) -> Result<T>,
{
    if plan.spectra.is_empty() {
        return Ok(Vec::new());
    }
    let mut result = work.vector(regions.len())?;
    for (region, interval) in regions.iter().zip(&plan.intervals) {
        work.consume(1)?;
        work.consume(interval.len())?;
        let mut row = work.vector(interval.len())?;
        for index in interval.clone() {
            let spectrum = plan.spectra[index];
            if !plan.validated[index] {
                work.consume(spectrum.peaks.len())?;
                let mut previous = None;
                for peak in &spectrum.peaks {
                    finite(peak.mz, "aggregation m/z")?;
                    if previous.is_some_and(|mz| mz > peak.mz) {
                        return Err(Error::UnsortedData);
                    }
                    previous = Some(peak.mz);
                }
                plan.validated[index] = true;
            }
            work.search(spectrum.len())?;
            work.search(spectrum.len())?;
            let start = spectrum
                .peaks
                .partition_point(|peak| peak.mz < region.mz.min);
            let end = spectrum
                .peaks
                .partition_point(|peak| peak.mz <= region.mz.max);
            let peaks = &spectrum.peaks[start..end];
            // Charge validation and one reduction traversal before either runs.
            work.consume(peaks.len())?;
            work.consume(peaks.len())?;
            for peak in peaks {
                finite(f64::from(peak.intensity), "aggregation intensity")?;
            }
            let value = reducer(peaks)?;
            finite(value, "aggregation result")?;
            row.push(convert(spectrum, value)?);
        }
        result.push(row);
    }
    Ok(result)
}

fn xics<F>(
    experiment: &MSExperiment,
    regions: &[MzRtRegion],
    ms_level: u32,
    reducer: F,
    work: &mut Work,
) -> Result<Vec<MSChromatogram>>
where
    F: FnMut(&[Peak1D]) -> Result<f64>,
{
    let plan = plan(experiment, regions, ms_level, work)?;
    if plan.spectra.is_empty() {
        return Ok(Vec::new());
    }
    work.consume(regions.len())?;
    let mut result = work.vector(regions.len())?;
    for region in regions {
        let mz = (region.mz.min + region.mz.max) / 2.0;
        finite(mz, "XIC product m/z")?;
        let mut chromatogram = MSChromatogram::new();
        chromatogram.product.mz = mz;
        result.push(chromatogram);
    }
    let rows = reduce_plan(
        plan,
        regions,
        reducer,
        |spectrum, value| {
            let intensity = value as f32;
            finite(f64::from(intensity), "XIC intensity")?;
            Ok(ChromatogramPeak::new(spectrum.rt, intensity))
        },
        work,
    )?;
    work.consume(rows.len())?;
    for (chromatogram, peaks) in result.iter_mut().zip(rows) {
        chromatogram.peaks = peaks;
    }
    Ok(result)
}

fn sum_f32(peaks: &[Peak1D]) -> Result<f64> {
    let sum = peaks.iter().fold(0.0_f32, |sum, peak| sum + peak.intensity);
    finite(f64::from(sum), "aggregation f32 sum")?;
    Ok(f64::from(sum))
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
