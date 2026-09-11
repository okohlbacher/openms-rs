// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::{DataArray, MSExperiment, MSSpectrum, Peak1D, Peak2D, RichPeak2D};
use crate::metadata::{MetaInfo, MetaValueData};
use crate::{Error, Result};
use std::mem::size_of;

/// Limits for concrete, unfiltered MS1 point conversion. Returned prior import
/// ownership and ordinary caller clones/drops are outside these per-call bounds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Data2DLimits {
    pub max_spectra: usize,
    pub max_points: usize,
    pub max_arrays: usize,
    pub max_name_bytes: usize,
    pub max_work: usize,
    pub max_bytes: usize,
}
impl Default for Data2DLimits {
    fn default() -> Self {
        Self {
            max_spectra: 1_000_000,
            max_points: 10_000_000,
            max_arrays: 1_000_000,
            max_name_bytes: 1024 * 1024,
            max_work: 50_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}
impl MSExperiment {
    pub fn get_2d_data(&self) -> Result<Vec<Peak2D>> {
        self.get_2d_data_with_limits(Data2DLimits::default())
    }
    pub fn get_2d_data_with_limits(&self, limits: Data2DLimits) -> Result<Vec<Peak2D>> {
        let mut result = Vec::new();
        self.append_2d_data_with_limits(&mut result, limits)?;
        Ok(result)
    }
    pub fn append_2d_data(&self, output: &mut Vec<Peak2D>) -> Result<()> {
        self.append_2d_data_with_limits(output, Data2DLimits::default())
    }
    pub fn append_2d_data_with_limits(
        &self,
        output: &mut Vec<Peak2D>,
        limits: Data2DLimits,
    ) -> Result<()> {
        let mut work = Work::new(limits);
        cap(self.spectra.len(), limits.max_spectra)?;
        work.consume(self.spectra.len())?;
        let mut added = 0usize;
        for s in &self.spectra {
            if s.ms_level == 1 {
                added = add(added, s.peaks.len())?;
                cap(added, limits.max_points)?;
            }
        }
        if added == 0 {
            return Ok(());
        }
        let count = add(output.len(), added)?;
        cap(count, limits.max_points)?;
        work.consume(add(self.spectra.len(), count)?)?;
        let mut staged = reserved(count, &mut work)?;
        staged.extend_from_slice(output);
        for spectrum in &self.spectra {
            if spectrum.ms_level != 1 {
                continue;
            }
            for peak in &spectrum.peaks {
                let point = Peak2D::new(spectrum.rt, peak.mz, peak.intensity);
                validate_point(point)?;
                staged.push(point);
            }
        }
        *output = staged;
        Ok(())
    }
    /// Atomically install new MS1 spectra and clear current chromatograms and
    /// experiment metadata. The returned map owns all replaced data unchanged.
    pub fn set_2d_data(&mut self, input: &[Peak2D]) -> Result<MSExperiment> {
        self.set_2d_data_with_limits(input, Data2DLimits::default())
    }
    pub fn set_2d_data_with_limits(
        &mut self,
        input: &[Peak2D],
        limits: Data2DLimits,
    ) -> Result<MSExperiment> {
        let staged = import(input, &[], limits)?;
        Ok(std::mem::replace(self, staged))
    }
    /// Import selected numeric metadata as parallel f32 arrays in the requested
    /// order. Missing keys become quiet NaN; metadata and IDs are not copied.
    /// Returns all prior experiment ownership, like set_2d_data.
    pub fn set_2d_data_rich(
        &mut self,
        input: &[RichPeak2D],
        names: &[String],
    ) -> Result<MSExperiment> {
        self.set_2d_data_rich_with_limits(input, names, Data2DLimits::default())
    }
    pub fn set_2d_data_rich_with_limits(
        &mut self,
        input: &[RichPeak2D],
        names: &[String],
        limits: Data2DLimits,
    ) -> Result<MSExperiment> {
        let staged = import(input, names, limits)?;
        Ok(std::mem::replace(self, staged))
    }
}
trait PointSource {
    fn point(&self) -> Peak2D;
    fn meta(&self) -> Option<&MetaInfo> {
        None
    }
}
impl PointSource for Peak2D {
    fn point(&self) -> Peak2D {
        *self
    }
}
impl PointSource for RichPeak2D {
    fn point(&self) -> Peak2D {
        self.peak
    }
    fn meta(&self) -> Option<&MetaInfo> {
        Some(&self.metadata)
    }
}
fn import<T: PointSource>(
    input: &[T],
    names: &[String],
    limits: Data2DLimits,
) -> Result<MSExperiment> {
    let mut work = Work::new(limits);
    cap(input.len(), limits.max_points)?;
    work.consume(input.len())?;
    let mut sizes = Vec::<usize>::new();
    let mut previous = None;
    for item in input {
        let point = item.point();
        validate_point(point)?;
        if previous.is_some_and(|rt| rt > point.rt()) {
            return Err(Error::UnsortedData);
        }
        if previous != Some(point.rt()) {
            cap(add(sizes.len(), 1)?, limits.max_spectra)?;
            grow_push(&mut sizes, 1, &mut work)?;
        } else {
            let last = sizes
                .last_mut()
                .ok_or_else(|| invalid("missing 2D import group"))?;
            *last = add(*last, 1)?;
        }
        previous = Some(point.rt());
    }
    let mut staged = MSExperiment::default();
    if input.is_empty() {
        return Ok(staged);
    }
    let arrays = mul(sizes.len(), names.len())?;
    cap(arrays, limits.max_arrays)?;
    work.consume(names.len())?;
    let mut name_bytes = 0usize;
    for name in names {
        name_bytes = add(name_bytes, name.len())?;
        cap(name_bytes, limits.max_name_bytes)?;
    }
    // Charges scalar writes, requested metadata queries, row construction and
    // newly created descriptors. Unrequested metadata payload is never visited.
    // Include disposal of the newly built scalar/empty-description containers
    // on an error. Prior experiment ownership is returned, never disposed here.
    work.consume(add(
        input.len(),
        add(mul(sizes.len(), 2)?, mul(arrays, 3)?)?,
    )?)?;
    work.consume(mul(name_bytes, sizes.len())?)?;
    staged.spectra = reserved(sizes.len(), &mut work)?;
    let mut start = 0usize;
    for count in sizes {
        let end = add(start, count)?;
        let mut spectrum = MSSpectrum {
            rt: input[start].point().rt(),
            ms_level: 1,
            ..Default::default()
        };
        spectrum.peaks = reserved(count, &mut work)?;
        spectrum.float_data_arrays = reserved(names.len(), &mut work)?;
        for name in names {
            work.allocate::<u8>(name.len())?;
            let mut owned = String::new();
            owned.try_reserve_exact(name.len()).map_err(|_| limit())?;
            owned.push_str(name);
            spectrum
                .float_data_arrays
                .push(DataArray::new(owned, reserved(count, &mut work)?));
        }
        for item in &input[start..end] {
            let point = item.point();
            spectrum
                .peaks
                .push(Peak1D::new(point.mz(), point.intensity));
            for (name, array) in names.iter().zip(&mut spectrum.float_data_arrays) {
                // A conservative full-map key-comparison bound precharges each
                // BTree lookup without traversing unrelated values/units/lists.
                let meta = item
                    .meta()
                    .ok_or_else(|| invalid("metadata unavailable for 2D import"))?;
                work.consume(mul(add(meta.len(), 1)?, add(name.len(), 1)?)?)?;
                let value = match meta.get(name) {
                    None => f32::NAN,
                    Some(value) => {
                        let value = match value.data() {
                            // Source converts stored integer directly to float.
                            MetaValueData::Integer(value) => *value as f32,
                            MetaValueData::Float(value) => *value as f32,
                            _ => return Err(invalid("2D metadata value is not numeric")),
                        };
                        if !value.is_finite() {
                            return Err(invalid("2D metadata does not fit finite f32"));
                        }
                        value
                    }
                };
                array.data.push(value);
            }
        }
        staged.spectra.push(spectrum);
        start = end;
    }
    Ok(staged)
}
fn validate_point(point: Peak2D) -> Result<()> {
    if point.position.iter().all(|v| v.is_finite()) && point.intensity.is_finite() {
        Ok(())
    } else {
        Err(invalid("2D point coordinates and intensity must be finite"))
    }
}
fn invalid(s: &str) -> Error {
    Error::InvalidValue(s.into())
}
fn limit() -> Error {
    invalid("2D point conversion resource limit exceeded")
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or_else(limit)
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b).ok_or_else(limit)
}
fn cap(n: usize, max: usize) -> Result<()> {
    if n > max { Err(limit()) } else { Ok(()) }
}
struct Work {
    remaining: usize,
    bytes: usize,
}
impl Work {
    fn new(l: Data2DLimits) -> Self {
        Self {
            remaining: l.max_work,
            bytes: l.max_bytes,
        }
    }
    fn consume(&mut self, n: usize) -> Result<()> {
        self.remaining = self.remaining.checked_sub(n).ok_or_else(limit)?;
        Ok(())
    }
    fn allocate<T>(&mut self, n: usize) -> Result<()> {
        if n > 0 {
            self.bytes = self
                .bytes
                .checked_sub(add(mul(n, size_of::<T>())?, 64)?)
                .ok_or_else(limit)?;
        }
        Ok(())
    }
}
fn reserved<T>(n: usize, work: &mut Work) -> Result<Vec<T>> {
    work.allocate::<T>(n)?;
    let mut result = Vec::new();
    result.try_reserve_exact(n).map_err(|_| limit())?;
    Ok(result)
}
fn grow_push(values: &mut Vec<usize>, value: usize, work: &mut Work) -> Result<()> {
    if values.len() == values.capacity() {
        let capacity = add(values.len(), 1)?
            .checked_next_power_of_two()
            .ok_or_else(limit)?;
        work.allocate::<usize>(capacity)?;
        work.consume(values.len())?;
        values
            .try_reserve_exact(capacity - values.len())
            .map_err(|_| limit())?;
    }
    values.push(value);
    Ok(())
}
