// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Retention-time union of two chromatograms, and the remaining members of
//! `KERNEL/MSChromatogram.h` that had no native counterpart.
//!
//! Source, all at Core SDK `bc9cc12`: `MSChromatogram::mergePeaks` with its
//! file-local helper `setSumSimilarUnion`
//! (`src/openms/source/KERNEL/MSChromatogram.cpp:474-565`), the subrange search
//! overloads `RTBegin`/`RTEnd`/`PosBegin`/`PosEnd`
//! (`MSChromatogram.cpp:228-333`), the predicate sort
//! `template<class Predicate> void sort` (`MSChromatogram.h:240-252`),
//! `MSChromatogram::MZLess` (`MSChromatogram.cpp:37-40`) and
//! `MSChromatogram::operator==` (`MSChromatogram.cpp:61-69`). Every other public
//! member of that header already maps onto
//! [`MSChromatogram`](crate::kernel::MSChromatogram) in `src/kernel.rs`,
//! `src/kernel/acquisition_fields.rs` and `src/kernel/ranges.rs`; the complete
//! member table is in `docs/CHROMATOGRAM_MERGE_SUPPORT.md`.
//!
//! Two points merge when their retention times agree after scaling to
//! milliseconds and rounding — see
//! [`merge_rt_key`](crate::kernel::chromatogram_merge::merge_rt_key). The merged
//! point keeps the destination's retention time and carries the sum of the two
//! intensities, exactly as the source helper does.
//!
//! The merge is bounded by
//! [`ChromatogramMergeLimits`](crate::kernel::chromatogram_merge::ChromatogramMergeLimits)
//! and builds its result in a temporary, so any rejection leaves both
//! chromatograms untouched. The source parallelises nothing here and neither
//! does this module.

use super::{ChromatogramPeak, MSChromatogram};
use crate::concept::constants::user_param::MERGED_CHROMATOGRAM_MZS;
use crate::metadata::{MetaValue, MetaValueData};
use crate::{Error, Result};
use std::cmp::Ordering;
use std::mem::size_of;
use std::ops::Range;

/// Retention-time bucket two points must share to be summed into one.
///
/// Source `setSumSimilarUnion` compares `round(rt * 1000.0)`
/// (`MSChromatogram.cpp:495`), which the comment above it describes as "within
/// 1/1000 seconds". The description is looser than the code: the comparison is a
/// bucketing, not a distance. Retention times `0.00149` and `0.00151` are
/// 0.00002 s apart yet fall in buckets 1 and 2 and stay separate, while
/// `0.00050` and `0.00149` are 0.00099 s apart, share bucket 1 and are summed.
/// The bucketing is reproduced rather than replaced by a distance test, because
/// it decides which points are summed and any change would alter the merged
/// intensities.
///
/// Rust's [`f64::round`] breaks halves away from zero, as C's `round` does, so
/// the bucket boundaries agree.
///
/// A retention time large enough that `rt * 1000.0` overflows to infinity yields
/// an infinite key; two such points compare equal and are summed. The source
/// does the same multiplication and shares the behaviour.
pub fn merge_rt_key(rt: f64) -> f64 {
    (rt * 1000.0).round()
}

/// Ascending order by product m/z (source `MSChromatogram::MZLess`).
///
/// Source `MZLess::operator()` returns `a.getMZ() < b.getMZ()`, and `getMZ()`
/// returns `getProduct().getMZ()`, so this is a plain comparison of the public
/// [`MSChromatogram::product`](crate::kernel::MSChromatogram) m/z. It is a free
/// function rather than a unit struct for the same reason
/// [`mass_trace_mz_less`](crate::analysis::feature_hypothesis::mass_trace_mz_less)
/// is: Rust sorts take a closure, not a comparator type.
///
/// A NaN product m/z makes every comparison false, which is not a strict weak
/// ordering. The source comparator has the identical defect; sort a slice with
/// [`slice::sort_by`] only after checking the m/z values are finite.
pub fn chromatogram_mz_less(a: &MSChromatogram, b: &MSChromatogram) -> bool {
    a.product.mz < b.product.mz
}

/// What a merge does with the destination's per-point annotation arrays.
///
/// The source updates the peaks and nothing else, so the float, integer and
/// string data arrays keep their pre-merge length while the point count changes.
/// The header states this as `@note` "Peak level metadata stored in float_array
/// string_array and int_array of the destination MSChromatogram is not
/// guaranteed to be correct after merging" (`MSChromatogram.h:472`). The
/// consequence is stronger than the wording: the destination then fails its own
/// `checkDataArraySizes_`, so a later `sortByPosition()` or `select()` throws.
/// This port therefore makes the choice explicit instead of silent.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MergedDataArrays {
    /// Refuse the merge when either chromatogram carries a non-empty annotation
    /// array. This is the default: per-point annotations have no defined value
    /// for a summed point, so the operation declines rather than inventing one.
    /// Arrays with no entries are declared placeholders everywhere else in the
    /// kernel and do not block a merge; they survive it unchanged.
    #[default]
    Reject,
    /// Drop every annotation array from the destination. The result is
    /// consistent, and the loss is requested rather than silent.
    Drop,
    /// Leave the destination's annotation arrays exactly as they are, as the
    /// source does. The result can be inconsistent; see
    /// [`MSChromatogram::validate`](crate::kernel::MSChromatogram::validate).
    Source,
}

/// Whole-input ceilings for one [`MSChromatogram::merge_peaks`] call.
///
/// Native bounds with no source counterpart. They are checked before anything is
/// allocated or mutated, so exceeding one leaves both chromatograms unchanged.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ChromatogramMergeLimits {
    /// Combined point count of the two inputs, which bounds the result.
    pub max_peaks: usize,
    /// Entries the `merged_chromatogram_mzs` list may reach, including the one
    /// this call appends.
    pub max_merged_mzs: usize,
    /// Cumulative temporary storage for the merged points and the metadata list.
    pub max_bytes: usize,
}

impl Default for ChromatogramMergeLimits {
    fn default() -> Self {
        Self {
            max_peaks: 10_000_000,
            max_merged_mzs: 1_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}

/// Everything [`MSChromatogram::merge_peaks_with_options`] can vary.
///
/// [`Default`] is the native default: record nothing, refuse annotated inputs.
/// [`ChromatogramMergeOptions::source`] selects the source behaviour instead.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ChromatogramMergeOptions {
    /// Append the other chromatogram's product m/z to the destination's
    /// `merged_chromatogram_mzs` metadata list (source `add_meta`).
    pub add_meta: bool,
    /// How the destination's annotation arrays are treated.
    pub data_arrays: MergedDataArrays,
    /// Resource ceilings for this call.
    pub limits: ChromatogramMergeLimits,
}

/// Whole-input ceilings for one [`MSChromatogram::sort_by`] call.
///
/// Native bounds with no source counterpart, checked before the permutation is
/// allocated and before the predicate is first invoked. Evaluating the predicate
/// is the caller's cost and is not counted here, as in
/// `docs/MOBILOGRAM_SUPPORT.md`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ChromatogramSortLimits {
    /// Points the chromatogram may hold.
    pub max_peaks: usize,
    /// Cumulative temporary storage for the permutation.
    pub max_bytes: usize,
}

impl Default for ChromatogramSortLimits {
    fn default() -> Self {
        Self {
            max_peaks: 10_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}

impl ChromatogramMergeOptions {
    /// The source's own behaviour: annotation arrays are left untouched and may
    /// end up misaligned with the points. `add_meta` still defaults to false,
    /// matching the source's default argument.
    pub fn source() -> Self {
        Self {
            add_meta: false,
            data_arrays: MergedDataArrays::Source,
            limits: ChromatogramMergeLimits::default(),
        }
    }

    /// The same options with `add_meta` set, for the source's second argument.
    pub fn with_add_meta(mut self, add_meta: bool) -> Self {
        self.add_meta = add_meta;
        self
    }
}

impl MSChromatogram {
    /// Source equality: every field except the name.
    ///
    /// `MSChromatogram::operator==` carries the comment "name_ can differ => it
    /// is not checked; also ranges are not checked" (`MSChromatogram.cpp:63`).
    /// It compares the point vector, the inherited `ChromatogramSettings`
    /// (native id, comment, instrument settings, acquisition info, source file,
    /// precursor, product, deeply-compared processing handles, chromatogram type
    /// and the meta-info map) and the three annotation-array lists.
    ///
    /// The derived [`PartialEq`] compares the name as well, because a native
    /// value type that ignores one of its own fields is a trap. Use this method
    /// where the source predicate is what you need; the ranges are not part of
    /// either comparison, since this port computes them on demand rather than
    /// caching them (see
    /// [`MSChromatogram::range_manager`](crate::kernel::MSChromatogram::range_manager)).
    ///
    /// Shared processing handles compare by pointed-to value, not by address,
    /// matching the source's `boost::make_indirect_iterator` comparison.
    pub fn source_equal(&self, other: &Self) -> bool {
        // Every field of the struct except `name`. Spelled out rather than
        // derived so that adding a field is a compile error here, not a silent
        // hole in the predicate.
        let Self {
            instrument_settings,
            acquisition_info,
            source_file,
            data_processing,
            chromatogram_type,
            peaks,
            native_id,
            name: _,
            precursor,
            product,
            metadata,
            float_data_arrays,
            integer_data_arrays,
            string_data_arrays,
        } = self;
        *instrument_settings == other.instrument_settings
            && *acquisition_info == other.acquisition_info
            && *source_file == other.source_file
            && *data_processing == other.data_processing
            && *chromatogram_type == other.chromatogram_type
            && *peaks == other.peaks
            && *native_id == other.native_id
            && *precursor == other.precursor
            && *product == other.product
            && *metadata == other.metadata
            && *float_data_arrays == other.float_data_arrays
            && *integer_data_arrays == other.integer_data_arrays
            && *string_data_arrays == other.string_data_arrays
    }

    /// Index of the first point at or above `rt` within `range`, searching only
    /// that subrange (source `RTBegin(begin, rt, end)` and its `PosBegin` alias).
    ///
    /// The returned index addresses the chromatogram, not the subrange, so it
    /// can be used directly with [`MSChromatogram::peaks`]. An empty `range`
    /// returns `range.start`, matching the source's `lower_bound(begin, …,
    /// begin)`, which returns `begin`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] when the points inside `range` are not in
    /// non-decreasing retention-time order, and [`Error::InvalidValue`] when
    /// `rt` or a point's retention time is not finite, or when `range` is
    /// inverted or reaches past the end. The source `@note` documents unsorted
    /// input as undefined and checks nothing; an out-of-range iterator pair is
    /// undefined behaviour there.
    pub fn rt_begin_in(&self, rt: f64, range: Range<usize>) -> Result<usize> {
        self.rt_bound_in(rt, range, false)
    }

    /// Index of the first point strictly above `rt` within `range`, searching
    /// only that subrange (source `RTEnd(begin, rt, end)` and its `PosEnd`
    /// alias).
    ///
    /// The returned index addresses the chromatogram, not the subrange.
    ///
    /// # Errors
    ///
    /// As [`Self::rt_begin_in`].
    pub fn rt_end_in(&self, rt: f64, range: Range<usize>) -> Result<usize> {
        self.rt_bound_in(rt, range, true)
    }

    fn rt_bound_in(&self, rt: f64, range: Range<usize>, upper: bool) -> Result<usize> {
        if !rt.is_finite() {
            return Err(Error::InvalidValue(
                "query retention time must be finite".into(),
            ));
        }
        if range.start > range.end || range.end > self.peaks.len() {
            return Err(Error::InvalidValue(
                "retention time subrange is inverted or past the end".into(),
            ));
        }
        let points = &self.peaks[range.clone()];
        for peak in points {
            if !peak.rt.is_finite() {
                return Err(Error::InvalidValue(
                    "chromatogram retention time must be finite".into(),
                ));
            }
        }
        sorted(points)?;
        Ok(range.start
            + points.partition_point(|peak| if upper { peak.rt <= rt } else { peak.rt < rt }))
    }

    /// Stable sort by a caller-supplied comparison over point indices, moving
    /// every parallel annotation value with its point.
    ///
    /// This is source `template<class Predicate> void sort(const Predicate&
    /// lambda)`. The predicate receives the chromatogram and two indices into it
    /// — peaks or data arrays — and returns whether the first orders before the
    /// second; it must express a strict weak ordering. Equal keys keep their
    /// relative order, as the source's `std::stable_sort` does.
    ///
    /// The source `@note` "All data arrays are reordered alongside the peaks" is
    /// preserved. Its second `@note`, that cached ranges survive a permutation,
    /// has no counterpart: this port has no range cache and computes ranges on
    /// demand ([`MSChromatogram::range_manager`](crate::kernel::MSChromatogram::range_manager)).
    ///
    /// The predicate can be invoked twice for one comparison, because the
    /// standard-library sort is driven by a three-way ordering while the source
    /// predicate answers only "less". A predicate that is not a strict weak
    /// ordering panics here; the source's `std::stable_sort` is undefined in
    /// that case.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a non-empty annotation array's
    /// length differs from the point count, or when a ceiling in
    /// [`ChromatogramSortLimits`] is exceeded. As the source `@exception`
    /// requires, the array check happens before the predicate is ever invoked,
    /// so a mis-sized array leaves the chromatogram unchanged and the predicate
    /// unevaluated.
    pub fn sort_by(&mut self, less: impl FnMut(&Self, usize, usize) -> bool) -> Result<()> {
        self.sort_by_with_limits(less, ChromatogramSortLimits::default())
    }

    /// As [`Self::sort_by`], with explicit resource ceilings.
    ///
    /// # Errors
    ///
    /// As [`Self::sort_by`].
    pub fn sort_by_with_limits(
        &mut self,
        mut less: impl FnMut(&Self, usize, usize) -> bool,
        limits: ChromatogramSortLimits,
    ) -> Result<()> {
        let count = self.peaks.len();
        if count > limits.max_peaks {
            return Err(limit());
        }
        if count.checked_mul(size_of::<usize>()).ok_or_else(limit)? > limits.max_bytes {
            return Err(limit());
        }
        // Before the first, possibly array-indexing, predicate call.
        self.validate_data_arrays()?;
        let mut order: Vec<usize> = Vec::new();
        order.try_reserve_exact(count).map_err(|_| limit())?;
        order.extend(0..count);
        let this: &Self = self;
        order.sort_by(|&a, &b| {
            if less(this, a, b) {
                Ordering::Less
            } else if less(this, b, a) {
                Ordering::Greater
            } else {
                Ordering::Equal
            }
        });
        self.select(&order)
    }

    /// The recorded product m/z values of chromatograms merged into this one.
    ///
    /// Reads the `merged_chromatogram_mzs` meta value that
    /// [`Self::merge_peaks`] writes when `add_meta` is set. An absent key is not
    /// an error and yields an empty slice, because "nothing has been merged in"
    /// and "an empty merge list" are the same state.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the key exists but does not hold a
    /// float list. The source reaches the same condition through
    /// `DataValue::toDoubleList`, which throws `Exception::ConversionError`.
    pub fn merged_chromatogram_mzs(&self) -> Result<&[f64]> {
        match self.metadata.get(MERGED_CHROMATOGRAM_MZS) {
            None => Ok(&[]),
            Some(value) => value.as_float_list(),
        }
    }

    /// Add every point of `other`, summing points that share a retention-time
    /// bucket, and record the merge in this chromatogram's metadata.
    ///
    /// This is source `mergePeaks(MSChromatogram& other, bool add_meta = false)`
    /// with the native defaults: annotated inputs are refused rather than left
    /// misaligned (see [`MergedDataArrays`]). Call
    /// [`Self::merge_peaks_with_options`] with
    /// [`ChromatogramMergeOptions::source`] for the source behaviour.
    ///
    /// # Arguments
    ///
    /// * `other` — the chromatogram whose points are taken. The source signature
    ///   is a mutable reference, but the body only reads it, so this takes a
    ///   shared reference. That also makes `a.merge_peaks(&a, …)` a compile
    ///   error; the source permits the aliased call.
    /// * `add_meta` — when true, append `other`'s product m/z to this
    ///   chromatogram's `merged_chromatogram_mzs` list, creating the list if it
    ///   is absent. Read it back with [`Self::merged_chromatogram_mzs`].
    ///
    /// This chromatogram's own product m/z is unchanged, as the source states.
    ///
    /// # Errors
    ///
    /// As [`Self::merge_peaks_with_options`].
    pub fn merge_peaks(&mut self, other: &MSChromatogram, add_meta: bool) -> Result<()> {
        self.merge_peaks_with_options(
            other,
            ChromatogramMergeOptions {
                add_meta,
                ..Default::default()
            },
        )
    }

    /// [`Self::merge_peaks`] with an explicit annotation-array policy and
    /// explicit ceilings.
    ///
    /// Both inputs are validated and checked for retention-time order before
    /// anything is built, the merged points are assembled in a temporary, and
    /// the metadata value is constructed before either is committed. An error
    /// therefore leaves this chromatogram exactly as it was, and `other` is
    /// never written to at all.
    ///
    /// The source `@note` "Make sure BOTH chromatograms are sorted with respect
    /// to RT. Otherwise the result is undefined" is checked here instead of
    /// being left undefined: the cost is one pass over data the merge reads
    /// anyway, and an unnoticed unsorted input silently produces an unsorted,
    /// wrongly-summed result.
    ///
    /// # Errors
    ///
    /// * [`Error::UnsortedData`] when either chromatogram's points are not in
    ///   non-decreasing retention-time order.
    /// * [`Error::InvalidValue`] when either chromatogram fails
    ///   [`MSChromatogram::validate`], when a summed intensity is not finite
    ///   (the source lets it become infinite), when the existing
    ///   `merged_chromatogram_mzs` value is not a float list, when `other`'s
    ///   product m/z is not finite, or when a ceiling in
    ///   [`ChromatogramMergeLimits`] is exceeded.
    /// * [`Error::Unsupported`] when [`MergedDataArrays::Reject`] is in force
    ///   and either chromatogram carries a non-empty annotation array.
    pub fn merge_peaks_with_options(
        &mut self,
        other: &MSChromatogram,
        options: ChromatogramMergeOptions,
    ) -> Result<()> {
        let limits = options.limits;
        let total = self
            .peaks
            .len()
            .checked_add(other.peaks.len())
            .ok_or_else(limit)?;
        if total > limits.max_peaks {
            return Err(limit());
        }
        let mut bytes = total
            .checked_mul(size_of::<ChromatogramPeak>())
            .ok_or_else(limit)?;
        if options.data_arrays == MergedDataArrays::Reject
            && (has_values(&self.float_data_arrays)
                || has_values(&self.integer_data_arrays)
                || has_values(&self.string_data_arrays)
                || has_values(&other.float_data_arrays)
                || has_values(&other.integer_data_arrays)
                || has_values(&other.string_data_arrays))
        {
            return Err(Error::Unsupported(
                "merging chromatograms with per-point annotation arrays has no defined result; \
                 choose MergedDataArrays::Drop or ::Source"
                    .into(),
            ));
        }
        self.validate()?;
        other.validate()?;
        sorted(&self.peaks)?;
        sorted(&other.peaks)?;

        // The metadata list is prepared first: it is the only other fallible
        // step, and building it here keeps the commit below infallible.
        let merged_mzs = if options.add_meta {
            if !other.product.mz.is_finite() {
                return Err(Error::InvalidValue(
                    "merged chromatogram product m/z must be finite".into(),
                ));
            }
            let existing = self.merged_chromatogram_mzs()?;
            let length = existing.len().checked_add(1).ok_or_else(limit)?;
            if length > limits.max_merged_mzs {
                return Err(limit());
            }
            bytes = bytes
                .checked_add(length.checked_mul(size_of::<f64>()).ok_or_else(limit)?)
                .ok_or_else(limit)?;
            if bytes > limits.max_bytes {
                return Err(limit());
            }
            let mut list = Vec::new();
            list.try_reserve_exact(length).map_err(|_| limit())?;
            list.extend_from_slice(existing);
            list.push(other.product.mz);
            Some(MetaValue::new(MetaValueData::FloatList(list))?)
        } else {
            if bytes > limits.max_bytes {
                return Err(limit());
            }
            None
        };

        let mut merged: Vec<ChromatogramPeak> = Vec::new();
        merged.try_reserve_exact(total).map_err(|_| limit())?;
        let (mut left, mut right) = (0usize, 0usize);
        while left < self.peaks.len() || right < other.peaks.len() {
            // Source `setSumSimilarUnion` drains whichever side is exhausted
            // first, checking `first1 == last1` before `first2 == last2`.
            if left == self.peaks.len() {
                merged.push(other.peaks[right]);
                right += 1;
                continue;
            }
            if right == other.peaks.len() {
                merged.push(self.peaks[left]);
                left += 1;
                continue;
            }
            let here = merge_rt_key(self.peaks[left].rt);
            let there = merge_rt_key(other.peaks[right].rt);
            if here < there {
                merged.push(self.peaks[left]);
                left += 1;
            } else if there < here {
                merged.push(other.peaks[right]);
                right += 1;
            } else {
                // Approximately equal: the destination's retention time wins and
                // the intensities are added (`MSChromatogram.cpp:508-512`).
                let mut peak = self.peaks[left];
                peak.intensity += other.peaks[right].intensity;
                if !peak.intensity.is_finite() {
                    return Err(Error::InvalidValue(
                        "merged chromatogram intensity is not finite".into(),
                    ));
                }
                merged.push(peak);
                left += 1;
                right += 1;
            }
        }

        // Nothing below can fail.
        self.peaks = merged;
        if options.data_arrays == MergedDataArrays::Drop {
            self.float_data_arrays.clear();
            self.integer_data_arrays.clear();
            self.string_data_arrays.clear();
        }
        if let Some(value) = merged_mzs {
            self.metadata
                .insert(MERGED_CHROMATOGRAM_MZS.to_owned(), value);
        }
        Ok(())
    }
}

/// Source `operator<<(std::ostream&, const MSChromatogram&)`
/// (`MSChromatogram.cpp:19-35`): a banner, the chromatogram settings, one line
/// per point, a closing banner.
///
/// The settings block is reproduced exactly, which means two delimiters and
/// nothing between them: `operator<<(std::ostream&, const ChromatogramSettings&)`
/// (`ChromatogramSettings.cpp:149-154`) takes its argument unnamed and writes
/// only its own `BEGIN`/`END` lines, so no setting has ever appeared in this
/// dump. The banner is kept rather than dropped because callers and the class
/// test match on the surrounding text.
///
/// A formatter precision is forwarded to each point, as for
/// [`Mobilogram`](crate::kernel::Mobilogram). C++ stream locale and precision
/// state is not otherwise emulated.
impl std::fmt::Display for MSChromatogram {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(formatter, "-- MSCHROMATOGRAM BEGIN --")?;
        writeln!(formatter, "-- CHROMATOGRAMSETTINGS BEGIN --")?;
        writeln!(formatter, "-- CHROMATOGRAMSETTINGS END --")?;
        for point in &self.peaks {
            match formatter.precision() {
                Some(precision) => writeln!(formatter, "{point:.precision$}")?,
                None => writeln!(formatter, "{point}")?,
            }
        }
        writeln!(formatter, "-- MSCHROMATOGRAM END --")
    }
}

/// True when any array in the list holds at least one value. An array with no
/// entries is a declared placeholder, as everywhere else in the kernel, and does
/// not block a merge.
fn has_values<T>(arrays: &[super::DataArray<T>]) -> bool {
    arrays.iter().any(|array| !array.data.is_empty())
}

fn sorted(peaks: &[ChromatogramPeak]) -> Result<()> {
    if peaks.windows(2).any(|pair| pair[0].rt > pair[1].rt) {
        return Err(Error::UnsortedData);
    }
    Ok(())
}

fn limit() -> Error {
    Error::InvalidValue("chromatogram merge resource limit exceeded".into())
}
