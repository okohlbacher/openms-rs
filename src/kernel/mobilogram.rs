// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Owned one-dimensional mobility peaks and checked mobilogram operations.

use super::{DataArray, NumericRange, finite};
use crate::{Error, Result, metadata::DriftTimeUnit};
use std::{
    fmt,
    hash::{Hash, Hasher},
    mem::size_of,
    ops::Range,
};

/// A mobility coordinate in the container's unit and an f32 intensity.
/// Scalar field comparisons implement the source's comparator overloads.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MobilityPeak1D {
    pub mobility: f64,
    pub intensity: f32,
}
impl MobilityPeak1D {
    pub const DIMENSION: usize = 1;
    pub const fn new(mobility: f64, intensity: f32) -> Self {
        Self {
            mobility,
            intensity,
        }
    }
}
impl fmt::Display for MobilityPeak1D {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match f.precision() {
            Some(p) => write!(
                f,
                "POS: {:.*} INT: {:.*}",
                p, self.mobility, p, self.intensity
            ),
            None => write!(f, "POS: {} INT: {}", self.mobility, self.intensity),
        }
    }
}
impl Hash for MobilityPeak1D {
    fn hash<H: Hasher>(&self, state: &mut H) {
        (if self.mobility == 0.0 {
            0
        } else {
            self.mobility.to_bits()
        })
        .hash(state);
        (if self.intensity == 0.0 {
            0
        } else {
            self.intensity.to_bits()
        })
        .hash(state);
    }
}

/// Current peak bounds. Public edits cannot leave a stale cached range.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MobilogramRanges {
    pub mobility: Option<NumericRange>,
    pub intensity: Option<NumericRange>,
}

/// One operation's input counts, actual visits/comparisons, and temporary bytes.
/// Borrowed annotations are moved, never recursively cloned by sorting/selection.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MobilogramLimits {
    pub max_peaks: usize,
    pub max_arrays: usize,
    pub max_work: usize,
    pub max_bytes: usize,
}
impl Default for MobilogramLimits {
    fn default() -> Self {
        Self {
            max_peaks: 10_000_000,
            max_arrays: 100_000,
            max_work: 50_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}

/// Source mobilogram fields. Ordinary Vec operations replace exported C++
/// container methods. Full native equality includes parallel array descriptions.
#[derive(Clone, Debug, PartialEq)]
pub struct Mobilogram {
    pub peaks: Vec<MobilityPeak1D>,
    pub rt: f64,
    pub drift_time_unit: DriftTimeUnit,
    pub float_data_arrays: Vec<DataArray<f32>>,
    pub integer_data_arrays: Vec<DataArray<i32>>,
    pub string_data_arrays: Vec<DataArray<String>>,
}
impl Default for Mobilogram {
    fn default() -> Self {
        Self {
            peaks: Vec::new(),
            rt: -1.0,
            drift_time_unit: DriftTimeUnit::None,
            float_data_arrays: Vec::new(),
            integer_data_arrays: Vec::new(),
            string_data_arrays: Vec::new(),
        }
    }
}
impl From<Vec<MobilityPeak1D>> for Mobilogram {
    fn from(peaks: Vec<MobilityPeak1D>) -> Self {
        Self::from_peaks(peaks)
    }
}
impl Mobilogram {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn from_peaks(peaks: Vec<MobilityPeak1D>) -> Self {
        Self {
            peaks,
            ..Self::default()
        }
    }
    pub fn len(&self) -> usize {
        self.peaks.len()
    }
    pub fn is_empty(&self) -> bool {
        self.peaks.is_empty()
    }
    pub fn drift_time_unit_as_str(&self) -> &'static str {
        self.drift_time_unit.name()
    }

    /// Source equality deliberately ignores all arrays and cached ranges.
    pub fn source_equal(&self, other: &Self) -> bool {
        self.peaks == other.peaks
            && self.rt == other.rt
            && self.drift_time_unit == other.drift_time_unit
    }
    /// Exact source partial swap. Arrays stay with their original owner; this
    /// can leave different-length arrays invalid. Use std::mem::swap for a full swap.
    pub fn swap_peak_data(&mut self, other: &mut Self) {
        std::mem::swap(&mut self.peaks, &mut other.peaks);
        std::mem::swap(&mut self.rt, &mut other.rt);
        std::mem::swap(&mut self.drift_time_unit, &mut other.drift_time_unit);
    }
    /// Clear peaks and all parallel arrays; retain RT and drift-time unit.
    pub fn clear(&mut self) {
        self.peaks.clear();
        self.float_data_arrays.clear();
        self.integer_data_arrays.clear();
        self.string_data_arrays.clear();
    }

    /// Validate numeric waveform values and array lengths; unrelated annotation
    /// metadata are preserved without recursive validation.
    pub fn validate(&self) -> Result<()> {
        self.validate_with_limits(MobilogramLimits::default())
    }
    pub fn validate_with_limits(&self, limits: MobilogramLimits) -> Result<()> {
        let mut work = Work::new(self, limits)?;
        finite(self.rt, "mobilogram RT")?;
        work.waveform(self, true, true)?;
        work.arrays(self)?;
        Ok(())
    }
    pub fn ranges(&self) -> Result<MobilogramRanges> {
        self.ranges_with_limits(MobilogramLimits::default())
    }
    pub fn ranges_with_limits(&self, limits: MobilogramLimits) -> Result<MobilogramRanges> {
        let mut work = Work::new(self, limits)?;
        work.waveform(self, true, true)?;
        let mut result = MobilogramRanges::default();
        work.consume(self.len())?;
        for peak in &self.peaks {
            extend(&mut result.mobility, peak.mobility);
            extend(&mut result.intensity, f64::from(peak.intensity));
        }
        Ok(result)
    }
    pub fn is_sorted(&self) -> Result<bool> {
        self.is_sorted_with_limits(MobilogramLimits::default())
    }
    pub fn is_sorted_with_limits(&self, limits: MobilogramLimits) -> Result<bool> {
        let mut work = Work::new(self, limits)?;
        work.waveform(self, true, false)?;
        self.sorted_inner(
            &mut |this, a, b| this.peaks[a].mobility < this.peaks[b].mobility,
            &mut work,
        )
    }
    /// The predicate receives the unmodified container and two original indices.
    /// The caller supplies a strict weak ordering; its internal work is external.
    pub fn is_sorted_by(&self, less: impl FnMut(&Self, usize, usize) -> bool) -> Result<bool> {
        self.is_sorted_by_with_limits(less, MobilogramLimits::default())
    }
    pub fn is_sorted_by_with_limits(
        &self,
        mut less: impl FnMut(&Self, usize, usize) -> bool,
        limits: MobilogramLimits,
    ) -> Result<bool> {
        let mut work = Work::new(self, limits)?;
        // Source isSorted(predicate) does not validate arrays, unlike sort(predicate).
        self.sorted_inner(&mut less, &mut work)
    }
    fn sorted_inner(
        &self,
        less: &mut impl FnMut(&Self, usize, usize) -> bool,
        work: &mut Work,
    ) -> Result<bool> {
        for i in 1..self.len() {
            work.consume(1)?;
            if less(self, i, i - 1) {
                return Ok(false);
            }
        }
        Ok(true)
    }
    pub fn sort_by_position(&mut self) -> Result<()> {
        self.sort_by_position_with_limits(MobilogramLimits::default())
    }
    pub fn sort_by_position_with_limits(&mut self, limits: MobilogramLimits) -> Result<()> {
        let mut work = Work::new(self, limits)?;
        work.waveform(self, true, false)?;
        let mut less =
            |this: &Self, a: usize, b: usize| this.peaks[a].mobility < this.peaks[b].mobility;
        // Preserve the source's already-sorted no-op, including unused malformed arrays.
        if self.sorted_inner(&mut less, &mut work)? {
            return Ok(());
        }
        self.sort_inner(less, &mut work)
    }
    pub fn sort_by_intensity(&mut self, reverse: bool) -> Result<()> {
        self.sort_by_intensity_with_limits(reverse, MobilogramLimits::default())
    }
    pub fn sort_by_intensity_with_limits(
        &mut self,
        reverse: bool,
        limits: MobilogramLimits,
    ) -> Result<()> {
        let mut work = Work::new(self, limits)?;
        work.waveform(self, false, true)?;
        let mut less = |this: &Self, a: usize, b: usize| {
            if reverse {
                this.peaks[b].intensity < this.peaks[a].intensity
            } else {
                this.peaks[a].intensity < this.peaks[b].intensity
            }
        };
        if self.sorted_inner(&mut less, &mut work)? {
            return Ok(());
        }
        self.sort_inner(less, &mut work)
    }
    pub fn sort_by(&mut self, less: impl FnMut(&Self, usize, usize) -> bool) -> Result<()> {
        self.sort_by_with_limits(less, MobilogramLimits::default())
    }
    pub fn sort_by_with_limits(
        &mut self,
        less: impl FnMut(&Self, usize, usize) -> bool,
        limits: MobilogramLimits,
    ) -> Result<()> {
        self.sort_inner(less, &mut Work::new(self, limits)?)
    }
    fn sort_inner(
        &mut self,
        mut less: impl FnMut(&Self, usize, usize) -> bool,
        work: &mut Work,
    ) -> Result<()> {
        work.arrays(self)?; // Before the first possibly array-indexing predicate.
        let n = self.len();
        let mut order = work.vector::<usize>(n)?;
        let mut scratch = work.vector::<usize>(n)?;
        work.consume(mul(n, 2)?)?;
        order.extend(0..n);
        scratch.resize(n, 0);
        let mut width = 1;
        while width < n {
            let mut start = 0;
            while start < n {
                let middle = start.saturating_add(width).min(n);
                let end = middle.saturating_add(width).min(n);
                let (mut left, mut right) = (start, middle);
                for slot in &mut scratch[start..end] {
                    work.consume(2)?;
                    if right == end || (left < middle && !less(self, order[right], order[left])) {
                        *slot = order[left];
                        left += 1;
                    } else {
                        *slot = order[right];
                        right += 1;
                    }
                }
                start = end;
            }
            std::mem::swap(&mut order, &mut scratch);
            width = width.saturating_mul(2);
        }
        self.select_inner(&order, work)
    }
    /// Keep unique indices in caller order, moving parallel values and retaining
    /// every annotation description. Invalid indices or limits leave self intact.
    pub fn select(&mut self, indices: &[usize]) -> Result<()> {
        self.select_with_limits(indices, MobilogramLimits::default())
    }
    pub fn select_with_limits(
        &mut self,
        indices: &[usize],
        limits: MobilogramLimits,
    ) -> Result<()> {
        self.select_inner(indices, &mut Work::new(self, limits)?)
    }
    fn select_inner(&mut self, indices: &[usize], work: &mut Work) -> Result<()> {
        work.consume(indices.len())?;
        if indices.len() > self.len() {
            return Err(invalid("too many selected mobility indices"));
        }
        for &index in indices {
            if index >= self.len() {
                return Err(invalid("mobility index out of bounds"));
            }
        }
        let arrays = work.arrays(self)?;
        let mut destination = work.vector::<usize>(self.len())?;
        work.consume(mul(self.len(), 3)?)?;
        destination.resize(self.len(), usize::MAX);
        for (to, &from) in indices.iter().enumerate() {
            if destination[from] != usize::MAX {
                return Err(invalid("duplicate mobility index"));
            }
            destination[from] = to;
        }
        let mut omitted = indices.len();
        for target in &mut destination {
            if *target == usize::MAX {
                *target = omitted;
                omitted += 1;
            }
        }
        // All remaining work (including dropped String destructors) is reserved
        // before the first swap. Metadata and shared processing handles stay put.
        work.consume(mul(self.len(), arrays.checked_add(3).ok_or_else(limit)?)?)?;
        for i in 0..destination.len() {
            while destination[i] != i {
                let other = destination[i];
                self.peaks.swap(i, other);
                swap_arrays(&mut self.float_data_arrays, i, other);
                swap_arrays(&mut self.integer_data_arrays, i, other);
                swap_arrays(&mut self.string_data_arrays, i, other);
                destination.swap(i, other);
            }
        }
        self.peaks.truncate(indices.len());
        truncate(&mut self.float_data_arrays, indices.len());
        truncate(&mut self.integer_data_arrays, indices.len());
        truncate(&mut self.string_data_arrays, indices.len());
        Ok(())
    }

    pub fn mobility_begin(&self, mobility: f64) -> Result<usize> {
        self.mobility_begin_in(mobility, 0..self.len())
    }
    pub fn mobility_end(&self, mobility: f64) -> Result<usize> {
        self.mobility_end_in(mobility, 0..self.len())
    }
    pub fn mobility_begin_in(&self, mobility: f64, range: Range<usize>) -> Result<usize> {
        self.bound_with_limits(mobility, range, false, MobilogramLimits::default())
    }
    pub fn mobility_end_in(&self, mobility: f64, range: Range<usize>) -> Result<usize> {
        self.bound_with_limits(mobility, range, true, MobilogramLimits::default())
    }
    /// lower_bound, or upper_bound when upper=true, within the supplied subrange.
    /// Returned indices address this mobilogram, not the subrange.
    pub fn bound_with_limits(
        &self,
        mobility: f64,
        range: Range<usize>,
        upper: bool,
        limits: MobilogramLimits,
    ) -> Result<usize> {
        finite(mobility, "mobility query")?;
        let mut work = Work::new(self, limits)?;
        self.check_sorted_range(&range, &mut work)?;
        work.consume(usize::BITS as usize)?;
        Ok(range.start
            + self.peaks[range].partition_point(|p| {
                if upper {
                    p.mobility <= mobility
                } else {
                    p.mobility < mobility
                }
            }))
    }
    pub fn find_nearest(&self, mobility: f64) -> Result<Option<usize>> {
        self.find_nearest_with_limits(mobility, MobilogramLimits::default())
    }
    /// Empty input returns None, replacing the source no-tolerance exception.
    pub fn find_nearest_with_limits(
        &self,
        mobility: f64,
        limits: MobilogramLimits,
    ) -> Result<Option<usize>> {
        finite(mobility, "mobility query")?;
        self.nearest_inner(mobility, &mut Work::new(self, limits)?)
    }
    fn nearest_inner(&self, mobility: f64, work: &mut Work) -> Result<Option<usize>> {
        self.check_sorted_range(&(0..self.len()), work)?;
        work.consume(usize::BITS as usize)?;
        Ok(super::nearest(&self.peaks, mobility, |p| p.mobility))
    }
    pub fn find_nearest_with_tolerance(
        &self,
        mobility: f64,
        tolerance: f64,
    ) -> Result<Option<usize>> {
        self.find_nearest_with_tolerance_and_limits(
            mobility,
            tolerance,
            MobilogramLimits::default(),
        )
    }
    pub fn find_nearest_with_tolerance_and_limits(
        &self,
        mobility: f64,
        tolerance: f64,
        limits: MobilogramLimits,
    ) -> Result<Option<usize>> {
        finite(mobility, "mobility query")?;
        finite(tolerance, "mobility tolerance")?;
        let index = self.nearest_inner(mobility, &mut Work::new(self, limits)?)?;
        Ok(index.filter(|&i| {
            self.peaks[i].mobility >= mobility - tolerance
                && self.peaks[i].mobility <= mobility + tolerance
        }))
    }
    pub fn find_nearest_in_window(
        &self,
        mobility: f64,
        left: f64,
        right: f64,
    ) -> Result<Option<usize>> {
        self.find_nearest_in_window_with_limits(mobility, left, right, MobilogramLimits::default())
    }
    /// Source one-sided checks are retained even for finite negative tolerances.
    pub fn find_nearest_in_window_with_limits(
        &self,
        mobility: f64,
        left: f64,
        right: f64,
        limits: MobilogramLimits,
    ) -> Result<Option<usize>> {
        window(mobility, left, right)?;
        let Some(mut i) = self.nearest_inner(mobility, &mut Work::new(self, limits)?)? else {
            return Ok(None);
        };
        if self.peaks[i].mobility < mobility {
            if self.peaks[i].mobility >= mobility - left {
                return Ok(Some(i));
            }
            i += 1;
            Ok((i < self.len() && self.peaks[i].mobility <= mobility + right).then_some(i))
        } else {
            if self.peaks[i].mobility <= mobility + right {
                return Ok(Some(i));
            }
            Ok(i.checked_sub(1)
                .filter(|&j| self.peaks[j].mobility >= mobility - left))
        }
    }
    pub fn find_highest_in_window(
        &self,
        mobility: f64,
        left: f64,
        right: f64,
    ) -> Result<Option<usize>> {
        self.find_highest_in_window_with_limits(mobility, left, right, MobilogramLimits::default())
    }
    pub fn find_highest_in_window_with_limits(
        &self,
        mobility: f64,
        left: f64,
        right: f64,
        limits: MobilogramLimits,
    ) -> Result<Option<usize>> {
        window(mobility, left, right)?;
        let mut work = Work::new(self, limits)?;
        self.check_sorted_range(&(0..self.len()), &mut work)?;
        if self.is_empty() {
            return Ok(None);
        }
        let (low, high) = (mobility - left, mobility + right);
        if low > high {
            return Err(invalid("inverted mobility window"));
        }
        work.consume(2 * usize::BITS as usize)?;
        let begin = self.peaks.partition_point(|p| p.mobility < low);
        let end = self.peaks.partition_point(|p| p.mobility <= high);
        self.base_inner(begin..end, &mut work)
    }
    pub fn base_peak_index(&self) -> Result<Option<usize>> {
        self.base_peak_index_with_limits(MobilogramLimits::default())
    }
    pub fn base_peak_index_with_limits(&self, limits: MobilogramLimits) -> Result<Option<usize>> {
        self.base_inner(0..self.len(), &mut Work::new(self, limits)?)
    }
    pub fn base_peak(&self) -> Result<Option<&MobilityPeak1D>> {
        Ok(self.base_peak_index()?.map(|i| &self.peaks[i]))
    }
    pub fn base_peak_mut(&mut self) -> Result<Option<&mut MobilityPeak1D>> {
        Ok(self.base_peak_index()?.map(|i| &mut self.peaks[i]))
    }
    fn base_inner(&self, range: Range<usize>, work: &mut Work) -> Result<Option<usize>> {
        work.consume(range.len())?;
        let mut best = None;
        for i in range {
            finite(f64::from(self.peaks[i].intensity), "mobility intensity")?;
            if best.is_none_or(|j: usize| self.peaks[j].intensity < self.peaks[i].intensity) {
                best = Some(i);
            }
        }
        Ok(best)
    }
    pub fn calculate_tic(&self) -> Result<f32> {
        self.calculate_tic_with_limits(MobilogramLimits::default())
    }
    pub fn calculate_tic_with_limits(&self, limits: MobilogramLimits) -> Result<f32> {
        let mut work = Work::new(self, limits)?;
        work.consume(self.len())?;
        let mut sum = 0.0_f32;
        for peak in &self.peaks {
            finite(f64::from(peak.intensity), "mobility intensity")?;
            sum += peak.intensity;
            finite(f64::from(sum), "mobility TIC")?;
        }
        Ok(sum)
    }
    fn check_sorted_range(&self, range: &Range<usize>, work: &mut Work) -> Result<()> {
        if range.start > range.end || range.end > self.len() {
            return Err(invalid("invalid mobility subrange"));
        }
        work.consume(range.len())?;
        let mut previous = None;
        for peak in &self.peaks[range.clone()] {
            finite(peak.mobility, "mobility coordinate")?;
            if previous.is_some_and(|p| p > peak.mobility) {
                return Err(Error::UnsortedData);
            }
            previous = Some(peak.mobility);
        }
        Ok(())
    }
}
impl fmt::Display for Mobilogram {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "-- MOBILOGRAM BEGIN --")?;
        for peak in &self.peaks {
            match f.precision() {
                Some(p) => writeln!(f, "{peak:.p$}")?,
                None => writeln!(f, "{peak}")?,
            }
        }
        writeln!(f, "-- MOBILOGRAM END --")
    }
}
fn swap_arrays<T>(arrays: &mut [DataArray<T>], a: usize, b: usize) {
    for array in arrays {
        if !array.data.is_empty() {
            array.data.swap(a, b);
        }
    }
}
fn truncate<T>(arrays: &mut [DataArray<T>], len: usize) {
    for array in arrays {
        array.data.truncate(len);
    }
}
fn window(center: f64, left: f64, right: f64) -> Result<()> {
    finite(center, "mobility query")?;
    finite(left, "left mobility tolerance")?;
    finite(right, "right mobility tolerance")
}
fn extend(range: &mut Option<NumericRange>, value: f64) {
    match range {
        Some(r) => {
            if value < r.min {
                r.min = value;
            }
            if value > r.max {
                r.max = value;
            }
        }
        None => {
            *range = Some(NumericRange {
                min: value,
                max: value,
            })
        }
    }
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn limit() -> Error {
    invalid("mobilogram resource limit exceeded")
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b).ok_or_else(limit)
}
struct Work {
    remaining: usize,
    bytes: usize,
    arrays: usize,
}
impl Work {
    fn new(value: &Mobilogram, limits: MobilogramLimits) -> Result<Self> {
        if value.len() > limits.max_peaks {
            return Err(limit());
        }
        Ok(Self {
            remaining: limits.max_work,
            bytes: limits.max_bytes,
            arrays: limits.max_arrays,
        })
    }
    fn consume(&mut self, amount: usize) -> Result<()> {
        self.remaining = self.remaining.checked_sub(amount).ok_or_else(limit)?;
        Ok(())
    }
    fn vector<T>(&mut self, size: usize) -> Result<Vec<T>> {
        self.bytes = self
            .bytes
            .checked_sub(mul(size, size_of::<T>())?)
            .ok_or_else(limit)?;
        let mut result = Vec::new();
        result.try_reserve_exact(size).map_err(|_| limit())?;
        Ok(result)
    }
    fn waveform(&mut self, value: &Mobilogram, mobility: bool, intensity: bool) -> Result<()> {
        self.consume(value.len())?;
        for peak in &value.peaks {
            if mobility {
                finite(peak.mobility, "mobility coordinate")?;
            }
            if intensity {
                finite(f64::from(peak.intensity), "mobility intensity")?;
            }
        }
        Ok(())
    }
    fn arrays(&mut self, value: &Mobilogram) -> Result<usize> {
        let count = value
            .float_data_arrays
            .len()
            .checked_add(value.integer_data_arrays.len())
            .and_then(|n| n.checked_add(value.string_data_arrays.len()))
            .ok_or_else(limit)?;
        if count > self.arrays {
            return Err(limit());
        }
        self.consume(count)?;
        for len in value
            .float_data_arrays
            .iter()
            .map(|a| a.data.len())
            .chain(value.integer_data_arrays.iter().map(|a| a.data.len()))
            .chain(value.string_data_arrays.iter().map(|a| a.data.len()))
        {
            if len != 0 && len != value.len() {
                return Err(invalid("unaligned mobility data array"));
            }
        }
        Ok(count)
    }
}
