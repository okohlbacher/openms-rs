// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Borrowed peak-area traversal, with inclusive scalar RT/m/z boundaries.
//! Spectrum-level mobility filtering requires a separate scan mobility model.

use super::{MSExperiment, MSSpectrum, MzRtRegion, NumericRange, Peak1D, finite};
use crate::{Error, Result};
use std::{
    iter::{Enumerate, FusedIterator},
    mem::size_of,
    slice,
    sync::Arc,
};

/// Inclusive area dimensions. None means unrestricted, including finite extrema.
/// This does not represent source RangeManager mobility or intensity dimensions.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AreaBounds {
    pub rt: Option<NumericRange>,
    pub mz: Option<NumericRange>,
}
impl AreaBounds {
    /// The argument order follows source areaBegin: RT first, then m/z.
    pub fn new(min_rt: f64, max_rt: f64, min_mz: f64, max_mz: f64) -> Result<Self> {
        let result = Self {
            rt: Some(NumericRange {
                min: min_rt,
                max: max_rt,
            }),
            mz: Some(NumericRange {
                min: min_mz,
                max: max_mz,
            }),
        };
        result.validate()?;
        Ok(result)
    }
    fn validate(self) -> Result<()> {
        for (range, label) in [
            (self.rt, "area RT boundary"),
            (self.mz, "area m/z boundary"),
        ] {
            if let Some(range) = range {
                finite(range.min, label)?;
                finite(range.max, label)?;
                if range.min > range.max {
                    return Err(invalid("area minimum exceeds maximum"));
                }
            }
        }
        Ok(())
    }
}
impl From<MzRtRegion> for AreaBounds {
    fn from(value: MzRtRegion) -> Self {
        Self {
            rt: Some(value.rt),
            mz: Some(value.mz),
        }
    }
}

/// Owned area settings. Exact level matching is used; zero is never a wildcard.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AreaOptions {
    pub bounds: AreaBounds,
    pub ms_level: u32,
}
impl Default for AreaOptions {
    fn default() -> Self {
        Self::new(AreaBounds::default(), 1)
    }
}
impl AreaOptions {
    /// Idiomatic native constructor: match the complete supplied u32 MS level.
    pub const fn new(bounds: AreaBounds, ms_level: u32) -> Self {
        Self { bounds, ms_level }
    }
    /// Reproduce source areaBegin's UInt -> uint8_t -> int8_t -> UInt conversion.
    /// For example 256 becomes 0, while 255 becomes u32::MAX. Ordinary new()
    /// intentionally does not narrow native requests.
    pub const fn source_compatible(bounds: AreaBounds, requested_ms_level: u32) -> Self {
        Self::new(bounds, (requested_ms_level as u8 as i8) as u32)
    }
}

/// Limits shared by validation, interval construction and reserved traversal.
/// No peaks or nested metadata are cloned. The plan has at most one record per scan.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AreaLimits {
    pub max_spectra: usize,
    /// Total raw peaks inspected, including nonmatching/out-of-area scans.
    pub max_peaks: usize,
    pub max_work: usize,
    pub max_bytes: usize,
}
impl Default for AreaLimits {
    fn default() -> Self {
        Self {
            max_spectra: 1_000_000,
            max_peaks: 10_000_000,
            max_work: 50_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}

/// One immutable area point; indices address the original experiment.
#[derive(Clone, Copy, Debug)]
pub struct AreaPeak<'a> {
    pub spectrum_index: usize,
    pub peak_index: usize,
    pub spectrum: &'a MSSpectrum,
    pub peak: &'a Peak1D,
}

/// One exclusive area point. RT/MS-level snapshots identify its scan without
/// aliasing an immutable whole-spectrum reference with its mutable peak.
#[derive(Debug)]
pub struct AreaPeakMut<'a> {
    pub spectrum_index: usize,
    pub peak_index: usize,
    pub rt: f64,
    pub ms_level: u32,
    pub peak: &'a mut Peak1D,
}

#[derive(Clone, Copy, Debug)]
struct ScanWindow {
    spectrum: usize,
    begin: usize,
    end: usize,
}
#[derive(Debug)]
struct Plan {
    windows: Vec<ScanWindow>,
    count: usize,
}

/// Forward, fused, exactly sized area iterator. Clone shares only its immutable
/// interval plan and keeps an independent cursor. Equality uses current peak
/// address (all exhausted iterators compare equal), matching the source.
#[derive(Clone, Debug, Default)]
pub struct AreaIter<'a> {
    spectra: &'a [MSSpectrum],
    windows: Option<Arc<Vec<ScanWindow>>>,
    window: usize,
    offset: usize,
    remaining: usize,
}
impl<'a> AreaIter<'a> {
    /// Source dereference/getSpectrum/getPeakIndex equivalents, checked at end.
    pub fn peek(&self) -> Option<AreaPeak<'a>> {
        if self.remaining == 0 {
            return None;
        }
        let window = self.windows.as_ref()?.get(self.window)?;
        let peak_index = window.begin + self.offset;
        let spectrum = &self.spectra[window.spectrum];
        Some(AreaPeak {
            spectrum_index: window.spectrum,
            peak_index,
            spectrum,
            peak: &spectrum.peaks[peak_index],
        })
    }
}
impl<'a> Iterator for AreaIter<'a> {
    type Item = AreaPeak<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        let result = self.peek()?;
        self.remaining -= 1;
        self.offset += 1;
        if let Some(window) = self.windows.as_ref().and_then(|w| w.get(self.window)) {
            if window.begin + self.offset == window.end {
                self.window += 1;
                self.offset = 0;
            }
        }
        Some(result)
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}
impl ExactSizeIterator for AreaIter<'_> {}
impl FusedIterator for AreaIter<'_> {}
impl PartialEq<AreaIter<'_>> for AreaIter<'_> {
    fn eq(&self, other: &AreaIter<'_>) -> bool {
        match (self.peek(), other.peek()) {
            (None, None) => true,
            (Some(a), Some(b)) => std::ptr::eq(a.peak, b.peak),
            _ => false,
        }
    }
}
impl Eq for AreaIter<'_> {}

#[derive(Debug)]
struct MutableScan<'a> {
    spectrum_index: usize,
    rt: f64,
    ms_level: u32,
    next_peak: usize,
    peaks: slice::IterMut<'a, Peak1D>,
}
/// Exclusive forward iterator. Each peak is yielded once; yielded references
/// may coexist. Cloning mutable iterator aliases is deliberately not exposed.
/// Selection endpoints are fixed before the first mutation, as in source scans.
#[derive(Debug)]
pub struct AreaIterMut<'a> {
    spectra: Enumerate<slice::IterMut<'a, MSSpectrum>>,
    windows: std::vec::IntoIter<ScanWindow>,
    current: Option<MutableScan<'a>>,
    remaining: usize,
}
impl<'a> Iterator for AreaIterMut<'a> {
    type Item = AreaPeakMut<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            return None;
        }
        loop {
            if let Some(scan) = &mut self.current {
                if let Some(peak) = scan.peaks.next() {
                    let peak_index = scan.next_peak;
                    scan.next_peak += 1;
                    self.remaining -= 1;
                    return Some(AreaPeakMut {
                        spectrum_index: scan.spectrum_index,
                        peak_index,
                        rt: scan.rt,
                        ms_level: scan.ms_level,
                        peak,
                    });
                }
            }
            self.current = None;
            let window = self.windows.next()?;
            for (index, spectrum) in self.spectra.by_ref() {
                if index == window.spectrum {
                    self.current = Some(MutableScan {
                        spectrum_index: index,
                        rt: spectrum.rt,
                        ms_level: spectrum.ms_level,
                        next_peak: window.begin,
                        peaks: spectrum.peaks[window.begin..window.end].iter_mut(),
                    });
                    break;
                }
            }
        }
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}
impl ExactSizeIterator for AreaIterMut<'_> {}
impl FusedIterator for AreaIterMut<'_> {}

impl MSExperiment {
    pub fn area_iter(&self, options: AreaOptions) -> Result<AreaIter<'_>> {
        self.area_iter_with_limits(options, AreaLimits::default())
    }
    pub fn area_iter_with_limits(
        &self,
        options: AreaOptions,
        limits: AreaLimits,
    ) -> Result<AreaIter<'_>> {
        let mut remaining_work = limits.max_work;
        let mut remaining_bytes = limits.max_bytes;
        self.area_iter_with_budget(
            options,
            limits.max_spectra,
            limits.max_peaks,
            &mut remaining_work,
            &mut remaining_bytes,
        )
    }
    /// Reuse area selection inside a larger operation without resetting budgets.
    pub(crate) fn area_iter_with_budget(
        &self,
        options: AreaOptions,
        max_spectra: usize,
        max_peaks: usize,
        remaining_work: &mut usize,
        remaining_bytes: &mut usize,
    ) -> Result<AreaIter<'_>> {
        let mut work = Work {
            remaining: *remaining_work,
            bytes: *remaining_bytes,
        };
        let result = plan_with_work(self, options, max_spectra, max_peaks, &mut work);
        *remaining_work = work.remaining;
        *remaining_bytes = work.bytes;
        let plan = result?;
        Ok(AreaIter {
            spectra: &self.spectra,
            windows: if plan.windows.is_empty() {
                None
            } else {
                Some(Arc::new(plan.windows))
            },
            window: 0,
            offset: 0,
            remaining: plan.count,
        })
    }
    pub fn area_iter_mut(&mut self, options: AreaOptions) -> Result<AreaIterMut<'_>> {
        self.area_iter_mut_with_limits(options, AreaLimits::default())
    }
    pub fn area_iter_mut_with_limits(
        &mut self,
        options: AreaOptions,
        limits: AreaLimits,
    ) -> Result<AreaIterMut<'_>> {
        let plan = plan(self, options, limits)?;
        Ok(AreaIterMut {
            spectra: self.spectra.iter_mut().enumerate(),
            windows: plan.windows.into_iter(),
            current: None,
            remaining: plan.count,
        })
    }
    /// Source-compatible scalar areaBeginConst wrapper, including MS-level narrowing.
    pub fn area_begin(
        &self,
        min_rt: f64,
        max_rt: f64,
        min_mz: f64,
        max_mz: f64,
        ms_level: u32,
    ) -> Result<AreaIter<'_>> {
        self.area_iter(AreaOptions::source_compatible(
            AreaBounds::new(min_rt, max_rt, min_mz, max_mz)?,
            ms_level,
        ))
    }
    /// Source-compatible scalar areaBegin wrapper with exclusive peak references.
    pub fn area_begin_mut(
        &mut self,
        min_rt: f64,
        max_rt: f64,
        min_mz: f64,
        max_mz: f64,
        ms_level: u32,
    ) -> Result<AreaIterMut<'_>> {
        self.area_iter_mut(AreaOptions::source_compatible(
            AreaBounds::new(min_rt, max_rt, min_mz, max_mz)?,
            ms_level,
        ))
    }
}

fn plan(experiment: &MSExperiment, options: AreaOptions, limits: AreaLimits) -> Result<Plan> {
    let mut work = Work {
        remaining: limits.max_work,
        bytes: limits.max_bytes,
    };
    plan_with_work(
        experiment,
        options,
        limits.max_spectra,
        limits.max_peaks,
        &mut work,
    )
}
fn plan_with_work(
    experiment: &MSExperiment,
    options: AreaOptions,
    max_spectra: usize,
    max_peaks: usize,
    work: &mut Work,
) -> Result<Plan> {
    options.bounds.validate()?;
    if experiment.spectra.len() > max_spectra {
        return Err(limit());
    }
    work.consume(experiment.spectra.len())?;
    let mut peak_count = 0usize;
    let mut previous_rt = None;
    // Source's isSorted(true) checks the whole experiment, not just selected
    // scans. Neither intensity, annotation alignment, nor metadata are consumed.
    for spectrum in &experiment.spectra {
        finite(spectrum.rt, "area spectrum RT")?;
        if previous_rt.is_some_and(|previous| previous > spectrum.rt) {
            return Err(Error::UnsortedData);
        }
        previous_rt = Some(spectrum.rt);
        peak_count = peak_count
            .checked_add(spectrum.peaks.len())
            .filter(|&count| count <= max_peaks)
            .ok_or_else(limit)?;
        work.consume(spectrum.peaks.len())?;
        let mut previous_mz = None;
        for peak in &spectrum.peaks {
            finite(peak.mz, "area peak m/z")?;
            if previous_mz.is_some_and(|previous| previous > peak.mz) {
                return Err(Error::UnsortedData);
            }
            previous_mz = Some(peak.mz);
        }
    }
    let (begin, end) = if let Some(rt) = options.bounds.rt {
        (
            bound(&experiment.spectra, rt.min, false, |s| s.rt, work)?,
            bound(&experiment.spectra, rt.max, true, |s| s.rt, work)?,
        )
    } else {
        (0, experiment.spectra.len())
    };
    let mut windows = Vec::new();
    let candidate_count = end - begin;
    if candidate_count != 0 {
        // One contiguous plan vector; fixed Arc header is charged even though
        // mutable iteration consumes the vector directly without an Arc.
        work.allocate(
            candidate_count
                .checked_mul(size_of::<ScanWindow>())
                .and_then(|n| n.checked_add(64))
                .ok_or_else(limit)?,
        )?;
        windows
            .try_reserve_exact(candidate_count)
            .map_err(|_| limit())?;
    }
    work.consume(candidate_count)?;
    let mut count = 0usize;
    for (index, spectrum) in experiment.spectra[begin..end].iter().enumerate() {
        if spectrum.ms_level != options.ms_level {
            continue;
        }
        let (first, last) = if let Some(mz) = options.bounds.mz {
            (
                bound(&spectrum.peaks, mz.min, false, |p| p.mz, work)?,
                bound(&spectrum.peaks, mz.max, true, |p| p.mz, work)?,
            )
        } else {
            (0, spectrum.peaks.len())
        };
        if first != last {
            count = count.checked_add(last - first).ok_or_else(limit)?;
            windows.push(ScanWindow {
                spectrum: begin + index,
                begin: first,
                end: last,
            });
        }
    }
    // Reserve all future scan and per-peak traversal before exposing mutable
    // references. Construction errors leave the experiment completely unchanged.
    work.consume(end)?;
    work.consume(count)?;
    Ok(Plan { windows, count })
}
fn bound<T>(
    values: &[T],
    query: f64,
    upper: bool,
    coordinate: impl Fn(&T) -> f64,
    work: &mut Work,
) -> Result<usize> {
    let (mut first, mut len) = (0, values.len());
    while len != 0 {
        work.consume(1)?;
        let half = len / 2;
        let middle = first + half;
        let value = coordinate(&values[middle]);
        if if upper { value <= query } else { value < query } {
            first = middle + 1;
            len -= half + 1;
        } else {
            len = half;
        }
    }
    Ok(first)
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn limit() -> Error {
    invalid("peak-area iteration resource limit exceeded")
}
struct Work {
    remaining: usize,
    bytes: usize,
}
impl Work {
    fn consume(&mut self, count: usize) -> Result<()> {
        self.remaining = self.remaining.checked_sub(count).ok_or_else(limit)?;
        Ok(())
    }
    fn allocate(&mut self, bytes: usize) -> Result<()> {
        self.bytes = self.bytes.checked_sub(bytes).ok_or_else(limit)?;
        Ok(())
    }
}
