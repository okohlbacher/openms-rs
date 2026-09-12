// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Borrowed peak-area traversal, with inclusive RT, m/z and scan-mobility
//! boundaries.
//!
//! Ports `KERNEL/AreaIterator.h` and the `areaBegin`/`areaBeginConst`/`areaEnd`
//! family of `KERNEL/MSExperiment.h`. The source's iterator pair becomes one
//! Rust iterator; the source's `AreaIterator::Param` named-parameter builder
//! becomes [`AreaOptions`](crate::kernel::AreaOptions). Scan mobility is
//! filtered on the spectrum's *scalar* drift time, exactly as the source's
//! `nextScan_` does, so a per-peak ion-mobility array never takes part in the
//! selection. `docs/AREA_ITERATION_SUPPORT.md` lists every source member and
//! its counterpart.

use super::ranges::{MSDim, RangeManager};
use super::{MSExperiment, MSSpectrum, MzRtRegion, NumericRange, Peak1D, finite};
use crate::{Error, Result};
use std::{
    iter::{Enumerate, FusedIterator},
    mem::size_of,
    slice,
    sync::Arc,
};

/// Inclusive RT and m/z area dimensions. `None` means unrestricted, including
/// the finite extrema.
///
/// Scan mobility is not one of these: it filters whole scans rather than peaks
/// and lives in [`AreaOptions::mobility`]. The source's `RangeManager`
/// intensity dimension has no counterpart at all, because the source's area
/// iterator never filters on intensity.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AreaBounds {
    /// Inclusive retention-time boundary of the area, in seconds.
    pub rt: Option<NumericRange>,
    /// Inclusive m/z boundary applied inside every selected scan, in Th.
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
    /// Inclusive RT and m/z boundaries.
    pub bounds: AreaBounds,
    /// Only scans of exactly this MS level are visited.
    pub ms_level: u32,
    /// Inclusive scan-mobility boundary on the spectrum's scalar drift time,
    /// or `None` for no mobility restriction.
    ///
    /// Source `AreaIterator::Param::lowIM`/`highIM` (`AreaIterator.h:88-100`),
    /// which the iterator turns into `RangeMobility{low_im_, high_im_}` and
    /// tests per scan with `containsMobility(getDriftTime())`
    /// (`AreaIterator.h:277-282`). The source defaults the pair to
    /// `lowest()`/`max()`, a range that contains every finite drift time, so
    /// `None` here selects the same scans. A restricted range makes the RT
    /// window span several ion-mobility frames and keep only the frames whose
    /// drift time falls inside it.
    pub mobility: Option<NumericRange>,
}
impl Default for AreaOptions {
    fn default() -> Self {
        Self::new(AreaBounds::default(), 1)
    }
}
impl AreaOptions {
    /// Idiomatic native constructor: match the complete supplied u32 MS level.
    pub const fn new(bounds: AreaBounds, ms_level: u32) -> Self {
        Self {
            bounds,
            ms_level,
            mobility: None,
        }
    }
    /// Reproduce source areaBegin's UInt -> uint8_t -> int8_t -> UInt conversion.
    /// For example 256 becomes 0, while 255 becomes u32::MAX. Ordinary new()
    /// intentionally does not narrow native requests.
    pub const fn source_compatible(bounds: AreaBounds, requested_ms_level: u32) -> Self {
        Self::new(bounds, source_level(requested_ms_level))
    }
    /// The same options restricted to an inclusive scan-mobility range.
    ///
    /// Chains the source's `Param::lowIM(min_im).highIM(max_im)` pair. Only the
    /// spectrum's scalar drift time is compared; a spectrum whose peaks carry an
    /// ion-mobility float data array but whose scalar drift time is unset still
    /// presents the source sentinel `-1`, and any range above `-1` therefore
    /// excludes it.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when either bound is not finite or
    /// `min_im` exceeds `max_im`. The source stores whatever it is given and
    /// lets `RangeBase::contains` decide, under which a reversed pair silently
    /// selects nothing and a NaN bound excludes every scan.
    pub fn with_mobility(mut self, min_im: f64, max_im: f64) -> Result<Self> {
        self.mobility = Some(NumericRange {
            min: min_im,
            max: max_im,
        });
        self.validate()?;
        Ok(self)
    }
    /// Area settings taken from a [`RangeManager`], as source
    /// `areaBegin(const RangeManagerType&, UInt)` does.
    ///
    /// The RT, m/z and mobility dimensions each restrict the area; an **empty**
    /// dimension does not, because the source reads it with
    /// `getNonEmptyRange()`, which answers `(lowest, max)` for an empty range
    /// (`RangeManager.h:278-284`). A dimension the manager does not carry at all
    /// is treated the same way — the source's manager is a fixed
    /// `RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility>` and
    /// cannot express a missing dimension, and "missing" and "empty" have the
    /// same effect. An intensity dimension is ignored: the source's area
    /// iterator has no intensity filter.
    ///
    /// `ms_level` is matched exactly. Use
    /// [`MSExperiment::area_begin_from_ranges`] for the source's byte-narrowing
    /// wrapper.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a non-empty dimension holds a
    /// non-finite bound.
    pub fn from_range_manager(range: &RangeManager, ms_level: u32) -> Result<Self> {
        let result = Self {
            bounds: AreaBounds {
                rt: dimension(range, MSDim::Rt)?,
                mz: dimension(range, MSDim::Mz)?,
            },
            ms_level,
            mobility: dimension(range, MSDim::Mobility)?,
        };
        result.validate()?;
        Ok(result)
    }
    fn validate(&self) -> Result<()> {
        self.bounds.validate()?;
        if let Some(range) = self.mobility {
            finite(range.min, "area mobility boundary")?;
            finite(range.max, "area mobility boundary")?;
            if range.min > range.max {
                return Err(invalid("area minimum exceeds maximum"));
            }
        }
        Ok(())
    }
}

/// Source narrowing of a requested MS level: `UInt` into the `uint8_t`
/// parameter of `AreaIterator::Param`, stored in its `int8_t ms_level_` field
/// and compared against the spectrum's unsigned level
/// (`AreaIterator.h:57`, `AreaIterator.h:124`, `AreaIterator.h:281`).
const fn source_level(requested: u32) -> u32 {
    (requested as u8 as i8) as u32
}

/// One dimension of a [`RangeManager`] as an optional inclusive area boundary.
/// A dimension that is absent or empty does not restrict the area.
fn dimension(range: &RangeManager, dim: MSDim) -> Result<Option<NumericRange>> {
    if !range.has_dim(dim) {
        return Ok(None);
    }
    let base = range.range_for_dim(dim)?;
    if base.is_empty() {
        return Ok(None);
    }
    Ok(Some(NumericRange {
        min: base.min()?,
        max: base.max()?,
    }))
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
    /// Index of the scan this peak belongs to, in the original experiment.
    pub spectrum_index: usize,
    /// Index of the peak inside that scan, in the original experiment.
    pub peak_index: usize,
    /// The scan itself (source `AreaIterator::getSpectrum`).
    pub spectrum: &'a MSSpectrum,
    /// The peak itself (source `AreaIterator::operator*`/`operator->`).
    pub peak: &'a Peak1D,
}
impl AreaPeak<'_> {
    /// Scan-wide ion mobility drift time, with the source sentinel `-1` when it
    /// is unset (source `AreaIterator::getDriftTime`, `AreaIterator.h:248-251`).
    ///
    /// Use [`MSSpectrum::drift_time_if_set`](crate::kernel::MSSpectrum::drift_time_if_set)
    /// through [`Self::spectrum`] to get `None` instead of the sentinel.
    pub fn drift_time(&self) -> f64 {
        self.spectrum.drift_time
    }
}

/// One exclusive area point. RT/MS-level snapshots identify its scan without
/// aliasing an immutable whole-spectrum reference with its mutable peak.
#[derive(Debug)]
pub struct AreaPeakMut<'a> {
    /// Index of the scan this peak belongs to, in the original experiment.
    pub spectrum_index: usize,
    /// Index of the peak inside that scan, in the original experiment.
    pub peak_index: usize,
    /// Retention time of that scan (source `AreaIterator::getRT`).
    pub rt: f64,
    /// MS level of that scan.
    pub ms_level: u32,
    /// Scan-wide ion mobility drift time of that scan, with the source
    /// sentinel `-1` when unset (source `AreaIterator::getDriftTime`).
    pub drift_time: f64,
    /// The peak itself, exclusively borrowed.
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
    drift_time: f64,
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
                        drift_time: scan.drift_time,
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
                        drift_time: spectrum.drift_time,
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
    /// Borrowed traversal of every peak inside an area.
    ///
    /// Ports `MSExperiment::areaBegin`/`areaEnd`, whose iterator pair becomes one
    /// Rust iterator.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for invalid bounds or an exceeded ceiling.
    pub fn area_iter(&self, options: AreaOptions) -> Result<AreaIter<'_>> {
        self.area_iter_with_limits(options, AreaLimits::default())
    }
    /// As [`Self::area_iter`], with explicit resource ceilings.
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
    /// Mutable traversal of every peak inside an area.
    ///
    /// # Errors
    ///
    /// As [`Self::area_iter`].
    pub fn area_iter_mut(&mut self, options: AreaOptions) -> Result<AreaIterMut<'_>> {
        self.area_iter_mut_with_limits(options, AreaLimits::default())
    }
    /// As [`Self::area_iter_mut`], with explicit resource ceilings.
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
    /// Borrowed traversal of an area given by a [`RangeManager`].
    ///
    /// Ports `areaBeginConst(const RangeManagerType& range, UInt ms_level)`
    /// (`MSExperiment.cpp:573-583`), including the source's MS-level byte
    /// narrowing. Empty and absent dimensions do not restrict the area; see
    /// [`AreaOptions::from_range_manager`].
    ///
    /// # Errors
    ///
    /// As [`Self::area_iter`], plus [`Error::InvalidValue`] for a non-finite
    /// bound in a non-empty dimension.
    pub fn area_begin_from_ranges(
        &self,
        range: &RangeManager,
        ms_level: u32,
    ) -> Result<AreaIter<'_>> {
        self.area_iter(AreaOptions::from_range_manager(
            range,
            source_level(ms_level),
        )?)
    }
    /// Exclusive traversal of an area given by a [`RangeManager`].
    ///
    /// Ports `areaBegin(const RangeManagerType& range, UInt ms_level)`
    /// (`MSExperiment.cpp:544-554`).
    ///
    /// # Errors
    ///
    /// As [`Self::area_begin_from_ranges`].
    pub fn area_begin_mut_from_ranges(
        &mut self,
        range: &RangeManager,
        ms_level: u32,
    ) -> Result<AreaIterMut<'_>> {
        self.area_iter_mut(AreaOptions::from_range_manager(
            range,
            source_level(ms_level),
        )?)
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
    options.validate()?;
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
    if options.mobility.is_some() {
        // One drift-time test per candidate scan, charged separately from the
        // level test above so a mobility-filtered call cannot exceed its budget.
        work.consume(candidate_count)?;
    }
    let mut count = 0usize;
    for (index, spectrum) in experiment.spectra[begin..end].iter().enumerate() {
        if spectrum.ms_level != options.ms_level {
            continue;
        }
        if let Some(mobility) = options.mobility {
            // Source nextScan_ skips a scan whose scalar drift time is outside
            // the mobility range (AreaIterator.h:281). A non-finite drift time
            // is refused rather than silently excluded, because RangeBase's
            // `min <= v & v <= max` answers false for NaN and the caller would
            // never learn that the scan was dropped.
            finite(spectrum.drift_time, "area spectrum drift time")?;
            if spectrum.drift_time < mobility.min || spectrum.drift_time > mobility.max {
                continue;
            }
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
