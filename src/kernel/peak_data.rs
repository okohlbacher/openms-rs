// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Source-compatible filtered f32 bulk exports. Row grouping follows the source
//! f64-versus-f32 comparison per peak, not spectrum identity.

use super::{AreaBounds, AreaIter, AreaOptions, MSExperiment};
use crate::{Error, Result};
use std::mem::size_of;

/// Bulk peak export as three parallel `f32` columns, one entry per peak.
///
/// Ports `MSExperiment::get2DPeakData`. Every peak of every selected spectrum
/// contributes one entry to each column, so the three vectors always share a
/// length and index together. Coordinates narrow to `f32` as in source.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FlatPeakData {
    /// Retention time of the spectrum each peak came from.
    pub rt: Vec<f32>,
    /// Mass-to-charge of each peak.
    pub mz: Vec<f32>,
    /// Intensity of each peak.
    pub intensity: Vec<f32>,
}

/// Bulk peak export grouped into one row per retention time.
///
/// Ports `MSExperiment::get2DPeakDataPerSpectrum`. The three vectors index
/// together, one entry per row.
///
/// Rows are **not** spectra. Source groups by comparing the `f64` retention
/// time against its `f32` narrowing per peak, so two spectra whose retention
/// times narrow to the same `f32` merge into one row, and a spectrum whose
/// retention time is not exactly representable in `f32` can split across rows.
/// This port preserves that grouping rather than spectrum identity.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpectrumPeakData {
    /// Retention time of each row.
    pub rt: Vec<f32>,
    /// Mass-to-charge values of each row's peaks.
    pub mz: Vec<Vec<f32>>,
    /// Intensities of each row's peaks.
    pub intensity: Vec<Vec<f32>>,
}

/// Cumulative ceilings for one export call.
///
/// Native bounds with no source counterpart. They cover area validation, any
/// output already present when appending, the newly produced output and the
/// scratch used to build it, and are checked before anything is written.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PeakDataLimits {
    /// Maximum spectra visited.
    pub max_spectra: usize,
    /// Maximum peaks visited.
    pub max_peaks: usize,
    /// Maximum points written, counting output already present when appending.
    pub max_output_points: usize,
    /// Maximum rows written, counting output already present when appending.
    pub max_output_rows: usize,
    /// Maximum weighted visits and comparisons.
    pub max_work: usize,
    /// Conservative ceiling on logical payload; an estimate, not measured.
    pub max_bytes: usize,
}
impl Default for PeakDataLimits {
    fn default() -> Self {
        Self {
            max_spectra: 1_000_000,
            max_peaks: 10_000_000,
            max_output_points: 10_000_000,
            max_output_rows: 1_000_000,
            max_work: 50_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}

impl MSExperiment {
    /// Every peak inside `bounds` at `ms_level`, as three parallel columns.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for invalid bounds or a nonfinite
    /// coordinate, or when [`PeakDataLimits`] is exceeded.
    pub fn get_2d_peak_data(&self, bounds: AreaBounds, ms_level: usize) -> Result<FlatPeakData> {
        self.get_2d_peak_data_with_limits(bounds, ms_level, PeakDataLimits::default())
    }
    /// As [`Self::get_2d_peak_data`], with explicit resource ceilings.
    pub fn get_2d_peak_data_with_limits(
        &self,
        bounds: AreaBounds,
        ms_level: usize,
        limits: PeakDataLimits,
    ) -> Result<FlatPeakData> {
        let mut result = FlatPeakData::default();
        self.append_2d_peak_data_with_limits(bounds, ms_level, &mut result, limits)?;
        Ok(result)
    }
    /// Append the selected peaks to existing columns, keeping their contents.
    ///
    /// # Errors
    ///
    /// As [`Self::get_2d_peak_data`]; the ceilings count what is already there.
    pub fn append_2d_peak_data(
        &self,
        bounds: AreaBounds,
        ms_level: usize,
        output: &mut FlatPeakData,
    ) -> Result<()> {
        self.append_2d_peak_data_with_limits(bounds, ms_level, output, PeakDataLimits::default())
    }
    /// As [`Self::append_2d_peak_data`], with explicit resource ceilings.
    pub fn append_2d_peak_data_with_limits(
        &self,
        bounds: AreaBounds,
        ms_level: usize,
        output: &mut FlatPeakData,
        limits: PeakDataLimits,
    ) -> Result<()> {
        let mut work = Work::new(limits);
        let points = selected(self, bounds, ms_level, limits, &mut work)?;
        if points.len() == 0 {
            return Ok(());
        }
        if output.rt.len() != output.mz.len() || output.rt.len() != output.intensity.len() {
            return Err(invalid("bulk peak output vectors are not aligned"));
        }
        let count = add(output.rt.len(), points.len())?;
        cap(count, limits.max_output_points)?;
        work.consume(mul(count, 3)?)?;
        let mut staged = FlatPeakData {
            rt: copied(&output.rt, count, &mut work)?,
            mz: copied(&output.mz, count, &mut work)?,
            intensity: copied(&output.intensity, count, &mut work)?,
        };
        for point in points {
            staged.rt.push(narrow(point.spectrum.rt, "RT")?);
            staged.mz.push(narrow(point.peak.mz, "m/z")?);
            finite_intensity(point.peak.intensity)?;
            staged.intensity.push(point.peak.intensity);
        }
        *output = staged;
        Ok(())
    }

    /// Every peak inside `bounds` at `ms_level`, grouped into rows by
    /// `f32`-narrowed retention time; see [`SpectrumPeakData`] on grouping.
    ///
    /// # Errors
    ///
    /// As [`MSExperiment::get_2d_peak_data`].
    pub fn get_2d_peak_data_per_spectrum(
        &self,
        bounds: AreaBounds,
        ms_level: usize,
    ) -> Result<SpectrumPeakData> {
        self.get_2d_peak_data_per_spectrum_with_limits(bounds, ms_level, PeakDataLimits::default())
    }
    /// As [`Self::get_2d_peak_data_per_spectrum`], with explicit ceilings.
    pub fn get_2d_peak_data_per_spectrum_with_limits(
        &self,
        bounds: AreaBounds,
        ms_level: usize,
        limits: PeakDataLimits,
    ) -> Result<SpectrumPeakData> {
        let mut result = SpectrumPeakData::default();
        self.append_2d_peak_data_per_spectrum_with_limits(bounds, ms_level, &mut result, limits)?;
        Ok(result)
    }
    /// Append the selected rows to existing row vectors, keeping their contents.
    ///
    /// # Errors
    ///
    /// As [`Self::get_2d_peak_data_per_spectrum`].
    pub fn append_2d_peak_data_per_spectrum(
        &self,
        bounds: AreaBounds,
        ms_level: usize,
        output: &mut SpectrumPeakData,
    ) -> Result<()> {
        self.append_2d_peak_data_per_spectrum_with_limits(
            bounds,
            ms_level,
            output,
            PeakDataLimits::default(),
        )
    }
    /// As [`Self::append_2d_peak_data_per_spectrum`], with explicit ceilings.
    pub fn append_2d_peak_data_per_spectrum_with_limits(
        &self,
        bounds: AreaBounds,
        ms_level: usize,
        output: &mut SpectrumPeakData,
        limits: PeakDataLimits,
    ) -> Result<()> {
        let mut work = Work::new(limits);
        let points = selected(self, bounds, ms_level, limits, &mut work)?;
        if points.len() == 0 {
            return Ok(());
        }
        let old_rows = output.rt.len();
        if old_rows != output.mz.len() || old_rows != output.intensity.len() {
            return Err(invalid("bulk peak row vectors are not aligned"));
        }
        cap(old_rows, limits.max_output_rows)?;
        work.consume(old_rows)?;
        let mut old_points = 0;
        for (mz, intensity) in output.mz.iter().zip(&output.intensity) {
            if mz.len() != intensity.len() {
                return Err(invalid("bulk peak row values are not aligned"));
            }
            old_points = add(old_points, mz.len())?;
            cap(old_points, limits.max_output_points)?;
        }
        cap(add(old_points, points.len())?, limits.max_output_points)?;
        // First pass validates conversion and determines exact source row sizes.
        // The shared area budget reserves one traversal; this extra pass and
        // copied/written scalar values are charged separately.
        work.consume(mul(points.len(), 4)?)?;
        // Count prior row descriptors and their eventual destruction as well
        // as the copied scalar buffers.
        work.consume(add(mul(old_points, 2)?, mul(old_rows, 5)?)?)?;
        let mut sizes: Vec<usize> = Vec::new();
        let mut old_last_extra = 0usize;
        let mut previous_rt = -1.0f32;
        for point in points.clone() {
            let rt = narrow(point.spectrum.rt, "RT")?;
            narrow(point.peak.mz, "m/z")?;
            finite_intensity(point.peak.intensity)?;
            if point.spectrum.rt != f64::from(previous_rt) {
                previous_rt = rt;
                cap(add(add(old_rows, sizes.len())?, 1)?, limits.max_output_rows)?;
                push_size(&mut sizes, &mut work)?;
            }
            if let Some(last) = sizes.last_mut() {
                *last = add(*last, 1)?;
            } else if old_rows != 0 {
                old_last_extra = add(old_last_extra, 1)?;
            } else {
                return Err(invalid(
                    "first selected RT is -1 but no existing output row exists",
                ));
            }
        }
        let rows = add(old_rows, sizes.len())?;
        work.consume(mul(sizes.len(), 2)?)?; // eventual nested-vector destruction
        // Old and new storage are charged together; only numeric payloads are
        // staged. Existing output is unchanged even on a late allocation error.
        let mut staged = SpectrumPeakData {
            rt: copied(&output.rt, rows, &mut work)?,
            mz: reserved(rows, &mut work)?,
            intensity: reserved(rows, &mut work)?,
        };
        for (i, (mz, intensity)) in output.mz.iter().zip(&output.intensity).enumerate() {
            let extra = if i + 1 == old_rows { old_last_extra } else { 0 };
            let capacity = add(mz.len(), extra)?;
            staged.mz.push(copied(mz, capacity, &mut work)?);
            staged
                .intensity
                .push(copied(intensity, capacity, &mut work)?);
        }
        previous_rt = -1.;
        let mut next_row = 0usize;
        for point in points {
            if point.spectrum.rt != f64::from(previous_rt) {
                previous_rt = point.spectrum.rt as f32;
                staged.rt.push(previous_rt);
                let size = *sizes
                    .get(next_row)
                    .ok_or_else(|| invalid("bulk peak row plan is inconsistent"))?;
                staged.mz.push(reserved(size, &mut work)?);
                staged.intensity.push(reserved(size, &mut work)?);
                next_row += 1;
            }
            // The counted plan guarantees a last row, including the source -1
            // case which appends directly to the caller's old final row.
            staged
                .mz
                .last_mut()
                .ok_or_else(|| invalid("bulk peak output has no current row"))?
                .push(point.peak.mz as f32);
            staged
                .intensity
                .last_mut()
                .ok_or_else(|| invalid("bulk peak output has no current row"))?
                .push(point.peak.intensity);
        }
        *output = staged;
        Ok(())
    }
}

fn selected<'a>(
    experiment: &'a MSExperiment,
    bounds: AreaBounds,
    ms_level: usize,
    limits: PeakDataLimits,
    work: &mut Work,
) -> Result<AreaIter<'a>> {
    // Source Size -> UInt -> uint8_t -> int8_t -> UInt, with explicit modular
    // narrowing on all native platforms. Zero is not a wildcard.
    let options = AreaOptions::source_compatible(bounds, ms_level as u32);
    experiment.area_iter_with_budget(
        options,
        limits.max_spectra,
        limits.max_peaks,
        &mut work.remaining,
        &mut work.bytes,
    )
}
fn narrow(value: f64, label: &str) -> Result<f32> {
    let result = value as f32;
    if !value.is_finite() || !result.is_finite() {
        return Err(invalid(&format!(
            "bulk peak {label} cannot be represented as finite f32"
        )));
    }
    Ok(result)
}
fn finite_intensity(value: f32) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(invalid("bulk peak intensity is nonfinite"))
    }
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn limit() -> Error {
    invalid("bulk peak export resource limit exceeded")
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or_else(limit)
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b).ok_or_else(limit)
}
fn cap(count: usize, max: usize) -> Result<()> {
    if count > max { Err(limit()) } else { Ok(()) }
}
struct Work {
    remaining: usize,
    bytes: usize,
}
impl Work {
    fn new(limits: PeakDataLimits) -> Self {
        Self {
            remaining: limits.max_work,
            bytes: limits.max_bytes,
        }
    }
    fn consume(&mut self, count: usize) -> Result<()> {
        self.remaining = self.remaining.checked_sub(count).ok_or_else(limit)?;
        Ok(())
    }
    fn allocate<T>(&mut self, count: usize) -> Result<()> {
        if count != 0 {
            self.bytes = self
                .bytes
                .checked_sub(add(mul(count, size_of::<T>())?, 64)?)
                .ok_or_else(limit)?;
        }
        Ok(())
    }
}
fn reserved<T>(count: usize, work: &mut Work) -> Result<Vec<T>> {
    work.allocate::<T>(count)?;
    let mut result = Vec::new();
    result.try_reserve_exact(count).map_err(|_| limit())?;
    Ok(result)
}
fn copied<T: Copy>(input: &[T], capacity: usize, work: &mut Work) -> Result<Vec<T>> {
    let mut result = reserved(capacity, work)?;
    result.extend_from_slice(input);
    Ok(result)
}
fn push_size(sizes: &mut Vec<usize>, work: &mut Work) -> Result<()> {
    if sizes.len() == sizes.capacity() {
        let capacity = add(sizes.capacity(), 1)?
            .checked_next_power_of_two()
            .ok_or_else(limit)?;
        work.allocate::<usize>(capacity)?;
        work.consume(sizes.len())?;
        sizes
            .try_reserve_exact(capacity - sizes.len())
            .map_err(|_| limit())?;
    }
    sizes.push(0);
    Ok(())
}
