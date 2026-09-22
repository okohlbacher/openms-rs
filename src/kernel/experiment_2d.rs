// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::{BaseFeature, DataArray, MSExperiment, MSSpectrum, Peak1D, Peak2D, RichPeak2D};
use crate::concept::constants::C13C12_MASSDIFF_U;
use crate::concept::constants::user_param::NUM_OF_MASSTRACES;
use crate::metadata::{MetaInfo, MetaValueData};
use crate::{Error, Result};
use std::fmt::Write as _;
use std::mem::size_of;
use std::ops::Deref;

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
    /// Every peak of the run as two-dimensional retention time / m/z points.
    ///
    /// Ports `MSExperiment::get2DData`. Unlike the filtered bulk exports in
    /// `peak_data`, this converts the whole run.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] on a nonfinite coordinate or an exceeded
    /// ceiling.
    pub fn get_2d_data(&self) -> Result<Vec<Peak2D>> {
        self.get_2d_data_with_limits(Data2DLimits::default())
    }
    /// As [`Self::get_2d_data`], with explicit resource ceilings.
    pub fn get_2d_data_with_limits(&self, limits: Data2DLimits) -> Result<Vec<Peak2D>> {
        let mut result = Vec::new();
        self.append_2d_data_with_limits(&mut result, limits)?;
        Ok(result)
    }
    /// Append every peak of the run to an existing vector, keeping its contents.
    ///
    /// # Errors
    ///
    /// As [`Self::get_2d_data`]; the ceilings count what is already there.
    pub fn append_2d_data(&self, output: &mut Vec<Peak2D>) -> Result<()> {
        self.append_2d_data_with_limits(output, Data2DLimits::default())
    }
    /// As [`Self::append_2d_data`], with explicit resource ceilings.
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
    /// As [`Self::set_2d_data`], with explicit resource ceilings.
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
    /// As [`Self::set_2d_data_rich`], with explicit resource ceilings.
    pub fn set_2d_data_rich_with_limits(
        &mut self,
        input: &[RichPeak2D],
        names: &[String],
        limits: Data2DLimits,
    ) -> Result<MSExperiment> {
        let staged = import(input, names, limits)?;
        Ok(std::mem::replace(self, staged))
    }
    /// Atomically install one MS1 spectrum per retention time of a feature
    /// container, expanding each feature's recorded mass traces into peaks, and
    /// clear current chromatograms and experiment metadata. The returned map
    /// owns all replaced data unchanged.
    ///
    /// Ports `set2DData<true>` (`MSExperiment.h:265-289`) and its
    /// `ContainerAdd_<T, true>::addData_` (`MSExperiment.h:1368-1393`), which
    /// `FileConverter` uses to turn a feature map into peaks. A feature with a
    /// `num_of_masstraces` meta value `n` becomes `n` peaks, the `i`-th at
    /// `mz + C13C12_MASSDIFF_U / charge * i` with the intensity of its
    /// `masstrace_intensity_<i>` meta value (charge `0` counts as `1`, and a
    /// negative charge spaces the peaks downward); the feature's own point is
    /// then not added. Any other feature becomes one peak at its m/z and
    /// intensity. As in the source, a spectrum opens for each new retention
    /// time before its features are added, so a retention time whose features
    /// all record zero mass traces yields an empty spectrum, and no metadata
    /// arrays are created.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] when the retention times decrease, and
    /// [`Error::InvalidValue`] when a `num_of_masstraces` is not a nonnegative
    /// integer (the source's `DataValue` conversion to `Size`), when an
    /// expected `masstrace_intensity_<i>` is missing (the source message
    /// `Meta value '…' expected but not found in container.`) or is neither an
    /// integer nor a float, on a nonfinite coordinate or intensity, or on an
    /// exceeded ceiling. Features are checked in order, so the first failing
    /// feature decides the error, as in the source. Nothing is replaced on an
    /// error, where the source has already cleared the experiment.
    ///
    /// ```
    /// use openms::kernel::{BaseFeature, Feature, MSExperiment};
    /// use openms::metadata::MetaValue;
    ///
    /// let mut feature = Feature::from(BaseFeature::new(18.5, 19.5, 17.5));
    /// feature.charge = 2;
    /// let meta = &mut feature.metadata;
    /// meta.insert("num_of_masstraces".into(), 2i64.into());
    /// meta.insert("masstrace_intensity_0".into(), MetaValue::try_from(11.0).unwrap());
    /// meta.insert("masstrace_intensity_1".into(), MetaValue::try_from(12.0).unwrap());
    /// let mut experiment = MSExperiment::new();
    /// experiment.set_2d_data_mass_traces(&[feature]).unwrap();
    /// let peaks = &experiment.spectra[0].peaks;
    /// assert_eq!((peaks[0].mz, peaks[0].intensity), (19.5, 11.0));
    /// assert_eq!((peaks[1].mz, peaks[1].intensity), (19.5 + 1.0033548378 / 2.0, 12.0));
    /// ```
    pub fn set_2d_data_mass_traces<F: Deref<Target = BaseFeature>>(
        &mut self,
        input: &[F],
    ) -> Result<MSExperiment> {
        self.set_2d_data_mass_traces_with_limits(input, Data2DLimits::default())
    }
    /// As [`Self::set_2d_data_mass_traces`], with explicit resource ceilings.
    /// `max_points` bounds the expanded peaks and `max_spectra` the retention
    /// times.
    pub fn set_2d_data_mass_traces_with_limits<F: Deref<Target = BaseFeature>>(
        &mut self,
        input: &[F],
        limits: Data2DLimits,
    ) -> Result<MSExperiment> {
        let staged = import_mass_traces(input, limits)?;
        Ok(std::mem::replace(self, staged))
    }
}
/// The source `set2DData<true>` loop: sort check and spectrum opening per new
/// retention time, then the feature's points, one feature at a time.
fn import_mass_traces<F: Deref<Target = BaseFeature>>(
    input: &[F],
    limits: Data2DLimits,
) -> Result<MSExperiment> {
    let mut work = Work::new(limits);
    work.consume(input.len())?;
    // (retention time, peak count) per spectrum, and the peaks in order.
    let mut groups = Vec::<(f64, usize)>::new();
    let mut peaks = Vec::<Peak1D>::new();
    let mut name = String::new();
    for feature in input {
        let rt = feature.rt;
        if !rt.is_finite() {
            return Err(invalid("2D point coordinates and intensity must be finite"));
        }
        let previous = groups.last().map(|&(previous, _)| previous);
        if previous != Some(rt) {
            if previous.is_some_and(|previous| previous > rt) {
                return Err(Error::UnsortedData);
            }
            cap(add(groups.len(), 1)?, limits.max_spectra)?;
            grow_push(&mut groups, (rt, 0), &mut work)?;
        }
        let before = peaks.len();
        match mass_trace_count(&feature.metadata)? {
            None => push_peak(&mut peaks, feature.mz, feature.intensity, limits, &mut work)?,
            Some(count) => {
                let charge = if feature.charge == 0 {
                    1
                } else {
                    feature.charge
                };
                let spacing = C13C12_MASSDIFF_U / f64::from(charge);
                // Bounded without trusting `count`: every iteration pushes a
                // peak under `max_points`, and at most `metadata.len() + 1`
                // distinct names can be looked up before one is missing.
                for i in 0..count {
                    name.clear();
                    write!(name, "masstrace_intensity_{i}").map_err(|_| limit())?;
                    let meta = &feature.metadata;
                    work.consume(mul(add(meta.len(), 1)?, add(name.len(), 1)?)?)?;
                    let value = meta.get(&name).ok_or_else(|| {
                        invalid(&format!(
                            "Meta value '{name}' expected but not found in container."
                        ))
                    })?;
                    let intensity = match value.data() {
                        // The source's DataValue conversion to float.
                        MetaValueData::Integer(value) => *value as f32,
                        MetaValueData::Float(value) => *value as f32,
                        _ => {
                            return Err(invalid(&format!(
                                "mass trace intensity '{name}' is not numeric"
                            )));
                        }
                    };
                    let mz = feature.mz + spacing * i as f64;
                    push_peak(&mut peaks, mz, intensity, limits, &mut work)?;
                }
            }
        }
        let added = peaks.len() - before;
        let last = groups
            .last_mut()
            .ok_or_else(|| invalid("missing 2D import group"))?;
        last.1 = add(last.1, added)?;
    }
    let mut staged = MSExperiment::default();
    work.consume(add(groups.len(), peaks.len())?)?;
    staged.spectra = reserved(groups.len(), &mut work)?;
    let mut start = 0usize;
    for (rt, count) in groups {
        let end = add(start, count)?;
        let mut spectrum = MSSpectrum {
            rt,
            ms_level: 1,
            ..Default::default()
        };
        spectrum.peaks = reserved(count, &mut work)?;
        spectrum.peaks.extend_from_slice(&peaks[start..end]);
        staged.spectra.push(spectrum);
        start = end;
    }
    Ok(staged)
}
/// The feature's `num_of_masstraces`, converted as the source's `DataValue`
/// conversion to `Size` does: only a nonnegative integer converts.
fn mass_trace_count(meta: &MetaInfo) -> Result<Option<usize>> {
    let Some(value) = meta.get(NUM_OF_MASSTRACES) else {
        return Ok(None);
    };
    match value.data() {
        MetaValueData::Integer(count) => usize::try_from(*count)
            .map(Some)
            .map_err(|_| invalid("num_of_masstraces is a negative integer")),
        _ => Err(invalid("num_of_masstraces is not an integer")),
    }
}
fn push_peak(
    peaks: &mut Vec<Peak1D>,
    mz: f64,
    intensity: f32,
    limits: Data2DLimits,
    work: &mut Work,
) -> Result<()> {
    if !mz.is_finite() || !intensity.is_finite() {
        return Err(invalid("2D point coordinates and intensity must be finite"));
    }
    cap(add(peaks.len(), 1)?, limits.max_points)?;
    grow_push(peaks, Peak1D::new(mz, intensity), work)
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
fn grow_push<T>(values: &mut Vec<T>, value: T, work: &mut Work) -> Result<()> {
    if values.len() == values.capacity() {
        let capacity = add(values.len(), 1)?
            .checked_next_power_of_two()
            .ok_or_else(limit)?;
        work.allocate::<T>(capacity)?;
        work.consume(values.len())?;
        values
            .try_reserve_exact(capacity - values.len())
            .map_err(|_| limit())?;
    }
    values.push(value);
    Ok(())
}
