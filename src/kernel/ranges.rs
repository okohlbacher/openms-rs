// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Multidimensional RT / m/z / intensity / mobility ranges.
//!
//! Ports `OpenMS/KERNEL/RangeManager.h`, `OpenMS/KERNEL/SpectrumRangeManager.h`
//! and `OpenMS/KERNEL/ChromatogramRangeManager.h` as pure value types. The
//! source keeps a mutable range cache inside every peak container and asks the
//! caller to refresh it with `updateRanges()`; this port never caches. Each
//! container computes its manager on demand ([`MSSpectrum::range_manager`],
//! [`MSChromatogram::range_manager`], [`Mobilogram::range_manager`],
//! [`MSExperiment::spectrum_range_manager`],
//! [`MSExperiment::chromatogram_range_manager`] and
//! [`MSExperiment::combined_range_manager`]), so a stale range cannot exist and
//! the source's `clearRanges()`-then-`updateRanges()` cycle has no counterpart.
//!
//! The source's `RangeManager<RangeBases...>` is a variadic template whose
//! dimension set is fixed at compile time. Rust has no variadic generics, so
//! [`RangeManager`](crate::kernel::ranges::RangeManager) carries its dimension
//! set at run time: an ordered set of [`MSDim`](crate::kernel::ranges::MSDim)
//! values, each with its own [`RangeBase`](crate::kernel::ranges::RangeBase).
//! Operations between two
//! managers act on the dimensions they have in common, exactly as the source's
//! fold expressions do, and "no dimension in common" is the same
//! `InvalidRange` error the source throws. See `docs/RANGES_SUPPORT.md`.

use super::{MSChromatogram, MSExperiment, MSSpectrum, Mobilogram, finite};
use crate::error::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

/// Dimensions of data acquisition for MS data, as the source `MSDim`.
///
/// The source `DIM_UNIT` (used by `RangeManager::clear`) is finer: its three
/// ion-mobility units `IM_MS`, `IM_VSSC` and `FAIMS_CV` all clear the one
/// mobility dimension. This port names the dimension directly.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MSDim {
    /// Retention time in seconds (source `MSDim::RT`, `RangeRT`).
    Rt,
    /// Mass-to-charge ratio (source `MSDim::MZ`, `RangeMZ`).
    Mz,
    /// Intensity (source `MSDim::INT`, `RangeIntensity`).
    Intensity,
    /// Ion mobility in the container's drift-time unit (source `MSDim::IM`,
    /// `RangeMobility`).
    Mobility,
}

impl MSDim {
    /// Every dimension, in the source declaration order.
    pub const ALL: [MSDim; 4] = [MSDim::Rt, MSDim::Mz, MSDim::Intensity, MSDim::Mobility];

    /// The label the source `operator<<` prints before a dimension's range.
    pub const fn label(self) -> &'static str {
        match self {
            MSDim::Rt => "rt",
            MSDim::Mz => "mz",
            MSDim::Intensity => "intensity",
            MSDim::Mobility => "mobility",
        }
    }

    const fn index(self) -> usize {
        match self {
            MSDim::Rt => 0,
            MSDim::Mz => 1,
            MSDim::Intensity => 2,
            MSDim::Mobility => 3,
        }
    }
}

/// State of the dimensions of a [`RangeManager`], as the source `HasRangeType`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum HasRangeType {
    /// All dimensions are filled.
    All,
    /// Some dimensions are empty, some are filled.
    Some,
    /// All dimensions are empty (cleared).
    None,
}

/// A simple closed range `[min, max]`, as the source `RangeBase`.
///
/// An empty range is represented, as in the source, by `min > max`; the
/// default is `min = f64::MAX`, `max = f64::MIN`. Shrinking operations can
/// produce other inverted pairs, and equality compares the stored endpoints
/// verbatim, so two empty ranges are equal only when their endpoints agree —
/// use [`Self::is_empty`] to test emptiness, as the source does.
///
/// Every stored endpoint is finite. The source stores any `double`: a NaN
/// endpoint silently breaks its `min <= max` invariant because every comparison
/// is false, and infinite endpoints turn later arithmetic into NaN. This port
/// rejects non-finite input on every entry point and rejects an arithmetic
/// result that leaves the finite domain, leaving the value unchanged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RangeBase {
    min: f64,
    max: f64,
}

impl Default for RangeBase {
    /// The empty range (source default constructor).
    fn default() -> Self {
        Self {
            min: f64::MAX,
            max: f64::MIN,
        }
    }
}

impl RangeBase {
    /// The empty range; [`Self::is_empty`] is `true`.
    pub fn new() -> Self {
        Self::default()
    }

    /// A singular range `[value, value]` (source `RangeBase(const double single)`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `value` is not finite; the source
    /// accepts any `double`.
    pub fn singular(value: f64) -> Result<Self> {
        finite(value, "range value")?;
        Ok(Self {
            min: value,
            max: value,
        })
    }

    /// The range `[min, max]` (source `RangeBase(const double min, const double max)`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when `min > max`, as the source throws
    /// `Exception::InvalidRange`, and [`Error::InvalidValue`] when either
    /// endpoint is not finite (a native check).
    pub fn from_min_max(min: f64, max: f64) -> Result<Self> {
        finite(min, "range minimum")?;
        finite(max, "range maximum")?;
        if min > max {
            return Err(Error::InvalidRange(
                "invalid initialization of range: minimum exceeds maximum".into(),
            ));
        }
        Ok(Self { min, max })
    }

    /// Make the range empty, so that [`Self::is_empty`] is `true`.
    pub fn clear(&mut self) {
        *self = Self::default();
    }

    /// Is the range empty, i.e. `min > max`?
    pub fn is_empty(&self) -> bool {
        self.min > self.max
    }

    /// Is `value` within `[min, max]`? Always `false` for an empty range.
    pub fn contains(&self, value: f64) -> bool {
        self.min <= value && value <= self.max
    }

    /// Is `inner` within `[min, max]`, i.e. are both of its endpoints contained?
    ///
    /// The source compares the raw endpoints, so an empty `inner` (`min > max`)
    /// is contained only when both of its inverted endpoints happen to lie
    /// inside this range. This port keeps that arithmetic; [`RangeManager::contains_all`]
    /// adds the source's "empty is contained" rule at the manager level.
    pub fn contains_range(&self, inner: &RangeBase) -> bool {
        self.contains(inner.min) && self.contains(inner.max)
    }

    /// Set the minimum, raising the maximum to it when the maximum is smaller
    /// (the source's "and the maximum, if uninitialized").
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `min` is not finite.
    pub fn set_min(&mut self, min: f64) -> Result<()> {
        finite(min, "range minimum")?;
        self.min = min;
        if self.max < min {
            self.max = min;
        }
        Ok(())
    }

    /// Set the maximum, lowering the minimum to it when the minimum is larger
    /// (the source's "and the minimum, if uninitialized").
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `max` is not finite.
    pub fn set_max(&mut self, max: f64) -> Result<()> {
        finite(max, "range maximum")?;
        self.max = max;
        if self.min > max {
            self.min = max;
        }
        Ok(())
    }

    /// The minimum value of the range.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when the range is empty, as the source
    /// throws `Exception::InvalidRange` ("Did you forget to call
    /// updateRanges()?"). Ranges are computed on demand here, so an empty range
    /// means the underlying data has no value in this dimension.
    pub fn min(&self) -> Result<f64> {
        if self.is_empty() {
            return Err(empty_range());
        }
        Ok(self.min)
    }

    /// The maximum value of the range.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when the range is empty, as [`Self::min`].
    pub fn max(&self) -> Result<f64> {
        if self.is_empty() {
            return Err(empty_range());
        }
        Ok(self.max)
    }

    /// Ensure the range includes the range of `other`.
    ///
    /// An empty `other` changes nothing. On equal endpoints the existing value
    /// is kept, as the source `std::min`/`std::max` keep their first argument
    /// (this preserves the sign of a zero endpoint).
    pub fn extend(&mut self, other: &RangeBase) {
        if other.min < self.min {
            self.min = other.min;
        }
        if other.max > self.max {
            self.max = other.max;
        }
    }

    /// Extend the range such that it includes `value`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `value` is not finite. The source
    /// silently ignores a NaN (both comparisons are false) and accepts
    /// infinities; this port refuses both so that a range can never hide a
    /// non-finite measurement.
    pub fn extend_value(&mut self, value: f64) -> Result<()> {
        finite(value, "range value")?;
        if value < self.min {
            self.min = value;
        }
        if value > self.max {
            self.max = value;
        }
        Ok(())
    }

    /// Extend the range by `by` units to the left and to the right.
    ///
    /// Negative values shrink the range; it may become empty. Calling this on
    /// an empty range has no effect.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `by` is not finite and
    /// [`Error::InvalidRange`] when an endpoint would leave the finite domain;
    /// the range is unchanged on error.
    pub fn extend_left_right(&mut self, by: f64) -> Result<()> {
        finite(by, "range extension")?;
        if self.is_empty() {
            return Ok(());
        }
        self.commit(self.min - by, self.max + by)
    }

    /// If the range is a single point (`min == max`), extend it by
    /// `min_span / 2` on either side, so that [`Self::span`] returns `min_span`.
    ///
    /// # Errors
    ///
    /// As [`Self::extend_left_right`].
    pub fn min_span_if_singular(&mut self, min_span: f64) -> Result<()> {
        finite(min_span, "minimum span")?;
        if self.min == self.max {
            self.extend_left_right(min_span / 2.0)?;
        }
        Ok(())
    }

    /// Ensure this range does not exceed the range of `other`.
    ///
    /// If `other` already contains this range, nothing changes. If this range
    /// lies entirely outside `other`, the result is empty. An empty range is
    /// not modified.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when `other` is empty, as the source
    /// throws `Exception::InvalidRange`. The source tests its own emptiness
    /// first, so an empty range clamped to an empty `other` is `Ok`; this port
    /// keeps that order.
    pub fn clamp_to(&mut self, other: &RangeBase) -> Result<()> {
        if self.is_empty() {
            return Ok(());
        }
        if other.is_empty() {
            return Err(Error::InvalidRange("cannot clamp to an empty range".into()));
        }
        self.min = self.min.max(other.min);
        self.max = self.max.min(other.max);
        Ok(())
    }

    /// Move this range into `sandbox` without changing its span, if possible.
    ///
    /// If the span exceeds the sandbox's span, the range is first cut to the
    /// sandbox's span (keeping its minimum); it is then shifted right when its
    /// minimum is below the sandbox or left when its maximum is above. An empty
    /// range is not modified.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when `sandbox` is empty, as the source
    /// throws `Exception::InvalidRange`, or when the shifted endpoints would
    /// leave the finite domain; the range is unchanged on error.
    pub fn push_into(&mut self, sandbox: &RangeBase) -> Result<()> {
        if self.is_empty() {
            return Ok(());
        }
        if sandbox.is_empty() {
            return Err(Error::InvalidRange(
                "cannot push into an empty range".into(),
            ));
        }
        if sandbox.contains_range(self) {
            return Ok(());
        }
        let mut min = self.min;
        let mut max = self.max;
        let sandbox_span = sandbox.max - sandbox.min;
        if max - min > sandbox_span {
            max = min + sandbox_span;
        }
        if min < sandbox.min {
            let distance = sandbox.min - min;
            min += distance;
            max += distance;
        } else if max > sandbox.max {
            let distance = sandbox.max - max;
            min += distance;
            max += distance;
        }
        self.commit(min, max)
    }

    /// Scale the range by `factor`; `> 1` widens it, `< 1` narrows it.
    ///
    /// With `d = max - min`, the new minimum is `min - d * (factor - 1) / 2`
    /// and the new maximum `max + d * (factor - 1) / 2`, so `scale_by(1.5)`
    /// extends the range by 25% on each side. Scaling an empty range has no
    /// effect; a singular range (`d == 0`) is unchanged.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `factor` is not finite and
    /// [`Error::InvalidRange`] when an endpoint would leave the finite domain;
    /// the range is unchanged on error.
    pub fn scale_by(&mut self, factor: f64) -> Result<()> {
        finite(factor, "scale factor")?;
        if self.is_empty() {
            return Ok(());
        }
        let extension = (self.max - self.min) * (factor - 1.0) / 2.0;
        self.commit(self.min - extension, self.max + extension)
    }

    /// Move the range by `distance`; negative values shift left. Shifting an
    /// empty range has no effect.
    ///
    /// # Errors
    ///
    /// As [`Self::scale_by`].
    pub fn shift(&mut self, distance: f64) -> Result<()> {
        finite(distance, "shift distance")?;
        if self.is_empty() {
            return Ok(());
        }
        self.commit(self.min + distance, self.max + distance)
    }

    /// The center point of the range, or `None` when empty.
    ///
    /// The source returns NaN for an empty range.
    pub fn center(&self) -> Option<f64> {
        if self.is_empty() {
            return None;
        }
        Some(self.min + (self.max - self.min) / 2.0)
    }

    /// The width `max - min` of the range, or `None` when empty.
    ///
    /// The source returns NaN for an empty range. The width of a range spanning
    /// nearly the whole finite domain can be `+inf`.
    pub fn span(&self) -> Option<f64> {
        if self.is_empty() {
            return None;
        }
        Some(self.max - self.min)
    }

    /// The current `(min, max)`, or the full finite domain
    /// `(f64::MIN, f64::MAX)` when empty, so that `min <= max` always holds.
    pub fn non_empty_range(&self) -> (f64, f64) {
        if self.is_empty() {
            (f64::MIN, f64::MAX)
        } else {
            (self.min, self.max)
        }
    }

    fn commit(&mut self, min: f64, max: f64) -> Result<()> {
        if !min.is_finite() || !max.is_finite() {
            return Err(Error::InvalidRange(
                "range arithmetic left the finite domain".into(),
            ));
        }
        self.min = min;
        self.max = max;
        Ok(())
    }
}

/// Formats as `[min, max]`, or `[, ]` for an empty range, as the source
/// `operator<<`. Numbers use Rust's shortest round-trip formatting, whereas the
/// source uses the stream's default six significant digits.
impl fmt::Display for RangeBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            write!(f, "[, ]")
        } else {
            write!(f, "[{}, {}]", self.min, self.max)
        }
    }
}

fn empty_range() -> Error {
    Error::InvalidRange(
        "empty or uninitialized range object: the data has no value in this dimension".into(),
    )
}

fn absent(dim: MSDim) -> Error {
    Error::InvalidValue(format!("range manager has no {} dimension", dim.label()))
}

fn no_overlap() -> Error {
    Error::InvalidRange("no dimensions in common".into())
}

fn limit() -> Error {
    Error::InvalidValue("range computation exceeds its item or byte limit".into())
}

/// A manager for a run-time set of ranges, one per [`MSDim`].
///
/// Ports the source `RangeManager<RangeBases...>`, whose dimension set is a
/// compile-time parameter pack: `RangeManager<RangeRT, RangeMZ>` for a
/// spectrum, and so on. Here the set is chosen at construction
/// ([`Self::new`] or one of the presets) and every typed accessor such as
/// [`Self::min_rt`] fails with [`Error::InvalidValue`] on a dimension the
/// manager does not carry — the source rejects that call at compile time.
///
/// Operations between two managers (`assign`, `extend`, `push_into`,
/// `clamp_to`, `contains_all`) act on the dimensions the two have in common
/// and leave the others untouched, as the source's `for_each_base_` folds do.
/// The `*_unsafe` variants report whether any dimension overlapped; the checked
/// variants return [`Error::InvalidRange`] when none did, as the source throws.
///
/// The source `RangeManagerContainer` adds `updateRanges()` and `getRange()`
/// to every peak container; those map to the on-demand `range_manager()`
/// accessors in this module.
#[derive(Clone, Copy, Debug)]
pub struct RangeManager {
    order: [Option<MSDim>; 4],
    ranges: [Option<RangeBase>; 4],
}

impl PartialEq for RangeManager {
    /// Equal when both carry the same dimensions with equal ranges. The
    /// declaration order affects only [`fmt::Display`].
    fn eq(&self, other: &Self) -> bool {
        self.ranges == other.ranges
    }
}

impl RangeManager {
    /// Items (peaks, mobility values, spectra and chromatograms) one on-demand
    /// range computation may visit. A native ceiling checked before any work.
    pub const MAX_ITEMS: usize = 100_000_000;
    /// Bytes of per-MS-level managers one [`MSExperiment::spectrum_range_manager`]
    /// call may allocate. A native ceiling checked before any allocation.
    pub const MAX_BYTES: usize = 64 * 1024 * 1024;

    /// A manager over `dims`, in that order, with every dimension empty.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `dims` is empty or names a
    /// dimension twice; both are compile errors in the source.
    pub fn new(dims: &[MSDim]) -> Result<Self> {
        if dims.is_empty() {
            return Err(Error::InvalidValue(
                "a range manager needs at least one dimension".into(),
            ));
        }
        let mut manager = Self {
            order: [None; 4],
            ranges: [None; 4],
        };
        for (slot, &dim) in dims.iter().enumerate() {
            if manager.ranges[dim.index()].is_some() {
                return Err(Error::InvalidValue(format!(
                    "duplicate {} dimension in range manager",
                    dim.label()
                )));
            }
            manager.ranges[dim.index()] = Some(RangeBase::default());
            manager.order[slot] = Some(dim);
        }
        Ok(manager)
    }

    fn with_dims(dims: &[MSDim]) -> Self {
        // Presets name at most four distinct dimensions, so `new` cannot fail.
        Self::new(dims).unwrap_or(Self {
            order: [None; 4],
            ranges: [None; 4],
        })
    }

    /// RT, m/z, intensity and mobility: the source `MSExperiment::RangeManagerType`
    /// and the class-test `RangeManagerContainer<RangeRT, RangeMZ, RangeIntensity, RangeMobility>`.
    pub fn experiment() -> Self {
        Self::with_dims(&MSDim::ALL)
    }

    /// m/z, intensity and mobility: the source `MSSpectrum::RangeManagerType`.
    pub fn spectrum() -> Self {
        Self::with_dims(&[MSDim::Mz, MSDim::Intensity, MSDim::Mobility])
    }

    /// RT and intensity: the source `MSChromatogram::RangeManagerType`.
    pub fn chromatogram() -> Self {
        Self::with_dims(&[MSDim::Rt, MSDim::Intensity])
    }

    /// Mobility and intensity: the source `Mobilogram::RangeManagerType`.
    pub fn mobilogram() -> Self {
        Self::with_dims(&[MSDim::Mobility, MSDim::Intensity])
    }

    /// RT, intensity and m/z: the source `ChromatogramRangeManager`
    /// (`RangeManager<RangeRT, RangeIntensity, RangeMZ>`), which adds nothing
    /// to its base beyond a `BaseType` alias.
    pub fn chromatogram_manager() -> Self {
        Self::with_dims(&[MSDim::Rt, MSDim::Intensity, MSDim::Mz])
    }

    /// m/z, intensity, mobility and RT: the source `SpectrumRangeManager::BaseType`.
    pub fn spectrum_manager() -> Self {
        Self::with_dims(&[MSDim::Mz, MSDim::Intensity, MSDim::Mobility, MSDim::Rt])
    }

    /// The dimensions this manager carries, in declaration order.
    pub fn dims(&self) -> impl Iterator<Item = MSDim> + '_ {
        self.order.iter().flatten().copied()
    }

    /// Does this manager carry `dim`?
    pub fn has_dim(&self, dim: MSDim) -> bool {
        self.ranges[dim.index()].is_some()
    }

    /// The range of `dim` (source `getRangeForDim`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the manager does not carry `dim`;
    /// the source only `assert`s.
    pub fn range_for_dim(&self, dim: MSDim) -> Result<&RangeBase> {
        self.ranges[dim.index()].as_ref().ok_or_else(|| absent(dim))
    }

    /// The mutable range of `dim` (source mutable `getRangeForDim`).
    ///
    /// # Errors
    ///
    /// As [`Self::range_for_dim`].
    pub fn range_for_dim_mut(&mut self, dim: MSDim) -> Result<&mut RangeBase> {
        self.ranges[dim.index()].as_mut().ok_or_else(|| absent(dim))
    }

    /// The minimum of `dim` (source `getMin*`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `dim` is absent, [`Error::InvalidRange`]
    /// when it is empty.
    pub fn min(&self, dim: MSDim) -> Result<f64> {
        self.range_for_dim(dim)?.min()
    }

    /// The maximum of `dim` (source `getMax*`).
    ///
    /// # Errors
    ///
    /// As [`Self::min`].
    pub fn max(&self, dim: MSDim) -> Result<f64> {
        self.range_for_dim(dim)?.max()
    }

    /// Set the minimum of `dim` (source `setMin*`); see [`RangeBase::set_min`].
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `dim` is absent or `min` is not finite.
    pub fn set_min(&mut self, dim: MSDim, min: f64) -> Result<()> {
        self.range_for_dim_mut(dim)?.set_min(min)
    }

    /// Set the maximum of `dim` (source `setMax*`); see [`RangeBase::set_max`].
    ///
    /// # Errors
    ///
    /// As [`Self::set_min`].
    pub fn set_max(&mut self, dim: MSDim, max: f64) -> Result<()> {
        self.range_for_dim_mut(dim)?.set_max(max)
    }

    /// Extend `dim` to include `value` (source `extend*`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `dim` is absent or `value` is not finite.
    pub fn extend_value(&mut self, dim: MSDim, value: f64) -> Result<()> {
        self.range_for_dim_mut(dim)?.extend_value(value)
    }

    /// Extend `dim` to include `range` (the source `extend(RangeMZ{100, 1500})`
    /// form, resolved through the inherited `RangeBase::extend`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `dim` is absent.
    pub fn extend_range(&mut self, dim: MSDim, range: &RangeBase) -> Result<()> {
        self.range_for_dim_mut(dim)?.extend(range);
        Ok(())
    }

    /// Is `value` within the range of `dim` (source `contains*`)?
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `dim` is absent.
    pub fn contains_value(&self, dim: MSDim, value: f64) -> Result<bool> {
        Ok(self.range_for_dim(dim)?.contains(value))
    }

    /// Is `inner` within the range of `dim` (source `contains*(const RangeBase&)`)?
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `dim` is absent.
    pub fn contains_range(&self, dim: MSDim, inner: &RangeBase) -> Result<bool> {
        Ok(self.range_for_dim(dim)?.contains_range(inner))
    }

    /// Is the range of `dim` empty (source `RangeRT::isEmpty()` and siblings)?
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `dim` is absent.
    pub fn is_dim_empty(&self, dim: MSDim) -> Result<bool> {
        Ok(self.range_for_dim(dim)?.is_empty())
    }

    /// Minimum RT (source `getMinRT`); see [`Self::min`].
    pub fn min_rt(&self) -> Result<f64> {
        self.min(MSDim::Rt)
    }
    /// Maximum RT (source `getMaxRT`); see [`Self::max`].
    pub fn max_rt(&self) -> Result<f64> {
        self.max(MSDim::Rt)
    }
    /// Set the minimum RT (source `setMinRT`); see [`Self::set_min`].
    pub fn set_min_rt(&mut self, min: f64) -> Result<()> {
        self.set_min(MSDim::Rt, min)
    }
    /// Set the maximum RT (source `setMaxRT`); see [`Self::set_max`].
    pub fn set_max_rt(&mut self, max: f64) -> Result<()> {
        self.set_max(MSDim::Rt, max)
    }
    /// Extend the RT range to include `rt` (source `extendRT`); see [`Self::extend_value`].
    pub fn extend_rt(&mut self, rt: f64) -> Result<()> {
        self.extend_value(MSDim::Rt, rt)
    }
    /// Is `rt` within the RT range (source `containsRT`)? See [`Self::contains_value`].
    pub fn contains_rt(&self, rt: f64) -> Result<bool> {
        self.contains_value(MSDim::Rt, rt)
    }
    /// Is `inner` within the RT range (source `containsRT(const RangeBase&)`)?
    pub fn contains_rt_range(&self, inner: &RangeBase) -> Result<bool> {
        self.contains_range(MSDim::Rt, inner)
    }

    /// Minimum m/z (source `getMinMZ`); see [`Self::min`].
    pub fn min_mz(&self) -> Result<f64> {
        self.min(MSDim::Mz)
    }
    /// Maximum m/z (source `getMaxMZ`); see [`Self::max`].
    pub fn max_mz(&self) -> Result<f64> {
        self.max(MSDim::Mz)
    }
    /// Set the minimum m/z (source `setMinMZ`); see [`Self::set_min`].
    pub fn set_min_mz(&mut self, min: f64) -> Result<()> {
        self.set_min(MSDim::Mz, min)
    }
    /// Set the maximum m/z (source `setMaxMZ`); see [`Self::set_max`].
    pub fn set_max_mz(&mut self, max: f64) -> Result<()> {
        self.set_max(MSDim::Mz, max)
    }
    /// Extend the m/z range to include `mz` (source `extendMZ`); see [`Self::extend_value`].
    pub fn extend_mz(&mut self, mz: f64) -> Result<()> {
        self.extend_value(MSDim::Mz, mz)
    }
    /// Is `mz` within the m/z range (source `containsMZ`)? See [`Self::contains_value`].
    pub fn contains_mz(&self, mz: f64) -> Result<bool> {
        self.contains_value(MSDim::Mz, mz)
    }
    /// Is `inner` within the m/z range (source `containsMZ(const RangeBase&)`)?
    pub fn contains_mz_range(&self, inner: &RangeBase) -> Result<bool> {
        self.contains_range(MSDim::Mz, inner)
    }

    /// Minimum intensity (source `getMinIntensity`); see [`Self::min`].
    pub fn min_intensity(&self) -> Result<f64> {
        self.min(MSDim::Intensity)
    }
    /// Maximum intensity (source `getMaxIntensity`); see [`Self::max`].
    pub fn max_intensity(&self) -> Result<f64> {
        self.max(MSDim::Intensity)
    }
    /// Set the minimum intensity (source `setMinIntensity`); see [`Self::set_min`].
    pub fn set_min_intensity(&mut self, min: f64) -> Result<()> {
        self.set_min(MSDim::Intensity, min)
    }
    /// Set the maximum intensity (source `setMaxIntensity`); see [`Self::set_max`].
    pub fn set_max_intensity(&mut self, max: f64) -> Result<()> {
        self.set_max(MSDim::Intensity, max)
    }
    /// Extend the intensity range to include `intensity` (source
    /// `extendIntensity`); see [`Self::extend_value`].
    pub fn extend_intensity(&mut self, intensity: f64) -> Result<()> {
        self.extend_value(MSDim::Intensity, intensity)
    }
    /// Is `intensity` within the intensity range (source `containsIntensity`)?
    pub fn contains_intensity(&self, intensity: f64) -> Result<bool> {
        self.contains_value(MSDim::Intensity, intensity)
    }
    /// Is `inner` within the intensity range (source
    /// `containsIntensity(const RangeBase&)`)?
    pub fn contains_intensity_range(&self, inner: &RangeBase) -> Result<bool> {
        self.contains_range(MSDim::Intensity, inner)
    }

    /// Minimum mobility (source `getMinMobility`); see [`Self::min`].
    pub fn min_mobility(&self) -> Result<f64> {
        self.min(MSDim::Mobility)
    }
    /// Maximum mobility (source `getMaxMobility`); see [`Self::max`].
    pub fn max_mobility(&self) -> Result<f64> {
        self.max(MSDim::Mobility)
    }
    /// Set the minimum mobility (source `setMinMobility`); see [`Self::set_min`].
    pub fn set_min_mobility(&mut self, min: f64) -> Result<()> {
        self.set_min(MSDim::Mobility, min)
    }
    /// Set the maximum mobility (source `setMaxMobility`); see [`Self::set_max`].
    pub fn set_max_mobility(&mut self, max: f64) -> Result<()> {
        self.set_max(MSDim::Mobility, max)
    }
    /// Extend the mobility range to include `mobility` (source
    /// `extendMobility`); see [`Self::extend_value`].
    pub fn extend_mobility(&mut self, mobility: f64) -> Result<()> {
        self.extend_value(MSDim::Mobility, mobility)
    }
    /// Is `mobility` within the mobility range (source `containsMobility`)?
    pub fn contains_mobility(&self, mobility: f64) -> Result<bool> {
        self.contains_value(MSDim::Mobility, mobility)
    }
    /// Is `inner` within the mobility range (source
    /// `containsMobility(const RangeBase&)`)?
    pub fn contains_mobility_range(&self, inner: &RangeBase) -> Result<bool> {
        self.contains_range(MSDim::Mobility, inner)
    }

    /// Dimensions carried by both managers, in this manager's order.
    fn common(&self, rhs: &RangeManager) -> [Option<MSDim>; 4] {
        let mut out = [None; 4];
        for (slot, dim) in self.dims().enumerate() {
            if rhs.has_dim(dim) {
                out[slot] = Some(dim);
            }
        }
        out
    }

    /// Copy every dimension in common from `rhs`; other dimensions are left
    /// untouched. Returns whether one or more dimensions overlapped.
    pub fn assign_unsafe(&mut self, rhs: &RangeManager) -> bool {
        let mut found = false;
        for dim in self.common(rhs).into_iter().flatten() {
            self.ranges[dim.index()] = rhs.ranges[dim.index()];
            found = true;
        }
        found
    }

    /// Copy every dimension in common from `rhs`; other dimensions are left
    /// untouched.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when no dimension overlapped, as the
    /// source throws.
    pub fn assign(&mut self, rhs: &RangeManager) -> Result<()> {
        if self.assign_unsafe(rhs) {
            Ok(())
        } else {
            Err(no_overlap())
        }
    }

    /// Extend every dimension in common to contain the range of `rhs`; other
    /// dimensions are left untouched. Returns whether any dimension overlapped
    /// (an empty `rhs` dimension counts as overlapping but changes nothing).
    pub fn extend_unsafe(&mut self, rhs: &RangeManager) -> bool {
        let mut found = false;
        for dim in self.common(rhs).into_iter().flatten() {
            if let (Some(mine), Some(theirs)) = (
                self.ranges[dim.index()].as_mut(),
                rhs.ranges[dim.index()].as_ref(),
            ) {
                mine.extend(theirs);
            }
            found = true;
        }
        found
    }

    /// Extend every dimension in common to contain the range of `rhs`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when no dimension overlapped, as the
    /// source throws.
    pub fn extend(&mut self, rhs: &RangeManager) -> Result<()> {
        if self.extend_unsafe(rhs) {
            Ok(())
        } else {
            Err(no_overlap())
        }
    }

    /// Apply [`RangeBase::scale_by`] to every dimension; empty and singular
    /// dimensions are unchanged.
    ///
    /// # Errors
    ///
    /// As [`RangeBase::scale_by`]; no dimension changes on error.
    pub fn scale_by(&mut self, factor: f64) -> Result<()> {
        let mut scaled = *self;
        for range in scaled.ranges.iter_mut().flatten() {
            range.scale_by(factor)?;
        }
        *self = scaled;
        Ok(())
    }

    /// Apply [`RangeBase::min_span_if_singular`] to every dimension; empty
    /// dimensions remain unchanged.
    ///
    /// # Errors
    ///
    /// As [`RangeBase::min_span_if_singular`]; no dimension changes on error.
    pub fn min_span_if_singular(&mut self, min_span: f64) -> Result<()> {
        let mut widened = *self;
        for range in widened.ranges.iter_mut().flatten() {
            range.min_span_if_singular(min_span)?;
        }
        *self = widened;
        Ok(())
    }

    /// Move every dimension in common into the corresponding dimension of
    /// `sandbox` without changing its span, if possible (see
    /// [`RangeBase::push_into`]). Dimensions absent from `sandbox` or empty in
    /// `sandbox` are left untouched. Returns whether any dimension overlapped.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when a shifted endpoint would leave the
    /// finite domain; no dimension changes on error.
    pub fn push_into_unsafe(&mut self, sandbox: &RangeManager) -> Result<bool> {
        let mut moved = *self;
        let mut found = false;
        for dim in self.common(sandbox).into_iter().flatten() {
            if let (Some(mine), Some(theirs)) = (
                moved.ranges[dim.index()].as_mut(),
                sandbox.ranges[dim.index()].as_ref(),
            ) {
                if !theirs.is_empty() {
                    mine.push_into(theirs)?;
                }
            }
            found = true;
        }
        *self = moved;
        Ok(found)
    }

    /// Move every dimension in common into `sandbox`, as [`Self::push_into_unsafe`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when no dimension overlapped, as the
    /// source throws, or as [`Self::push_into_unsafe`].
    pub fn push_into(&mut self, sandbox: &RangeManager) -> Result<()> {
        if self.push_into_unsafe(sandbox)? {
            Ok(())
        } else {
            Err(no_overlap())
        }
    }

    /// Clamp every dimension in common to the corresponding dimension of `rhs`
    /// (see [`RangeBase::clamp_to`]); this may tighten a dimension to a single
    /// point or empty it. Dimensions absent from `rhs` or empty in `rhs` are
    /// left untouched. Returns whether any dimension overlapped.
    pub fn clamp_to_unsafe(&mut self, rhs: &RangeManager) -> bool {
        let mut found = false;
        for dim in self.common(rhs).into_iter().flatten() {
            if let (Some(mine), Some(theirs)) = (
                self.ranges[dim.index()].as_mut(),
                rhs.ranges[dim.index()].as_ref(),
            ) {
                if !theirs.is_empty() {
                    // Cannot fail: `theirs` is non-empty and both are finite.
                    let _ = mine.clamp_to(theirs);
                }
            }
            found = true;
        }
        found
    }

    /// Clamp every dimension in common to `rhs`, as [`Self::clamp_to_unsafe`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when no dimension overlapped, as the
    /// source throws.
    pub fn clamp_to(&mut self, rhs: &RangeManager) -> Result<()> {
        if self.clamp_to_unsafe(rhs) {
            Ok(())
        } else {
            Err(no_overlap())
        }
    }

    /// Whether none, some or all of the dimensions are populated.
    pub fn has_range(&self) -> HasRangeType {
        let total = self.ranges.iter().flatten().count();
        let filled = self
            .ranges
            .iter()
            .flatten()
            .filter(|r| !r.is_empty())
            .count();
        if filled == 0 {
            HasRangeType::None
        } else if filled == total {
            HasRangeType::All
        } else {
            HasRangeType::Some
        }
    }

    /// Are all dimensions of `rhs` that overlap with this manager contained in
    /// this manager's ranges?
    ///
    /// An empty `rhs` dimension is considered contained (even when this
    /// manager's dimension is empty too); a non-empty `rhs` dimension is not
    /// contained in an empty dimension of this manager. If every overlapping
    /// dimension is empty in `rhs`, the result is `true`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when no dimension overlaps, as the
    /// source throws.
    pub fn contains_all(&self, rhs: &RangeManager) -> Result<bool> {
        let mut has_overlap = false;
        let mut contained = true;
        for dim in self.common(rhs).into_iter().flatten() {
            has_overlap = true;
            if let (Some(mine), Some(theirs)) =
                (&self.ranges[dim.index()], &rhs.ranges[dim.index()])
            {
                if theirs.is_empty() || mine.contains_range(theirs) {
                    continue;
                }
                contained = false;
            }
        }
        if !has_overlap {
            return Err(no_overlap());
        }
        Ok(contained)
    }

    /// Reset every dimension to empty (source `clearRanges`).
    pub fn clear_ranges(&mut self) {
        for range in self.ranges.iter_mut().flatten() {
            range.clear();
        }
    }

    /// Reset one dimension (source `clear(DIM_UNIT)`). If the manager does not
    /// carry `dim`, nothing happens, as in the source. The source's three
    /// ion-mobility units all name [`MSDim::Mobility`] here.
    pub fn clear_dim(&mut self, dim: MSDim) -> &mut Self {
        if let Some(range) = self.ranges[dim.index()].as_mut() {
            range.clear();
        }
        self
    }
}

/// Prints every dimension on its own line in declaration order, as the source
/// `printRange` / `operator<<`: `rt: [1, 1]\n`, `mz: [2, 2]\n`, and so on.
impl fmt::Display for RangeManager {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for dim in self.dims() {
            if let Some(range) = &self.ranges[dim.index()] {
                writeln!(f, "{}: {range}", dim.label())?;
            }
        }
        Ok(())
    }
}

/// Range manager for MS spectra with separate ranges for each MS level, as the
/// source `SpectrumRangeManager`.
///
/// A global manager over m/z, intensity, mobility and RT (the source base
/// class, [`RangeManager::spectrum_manager`]) covers every level, and a map from
/// MS level to a manager of the same shape covers each level separately. Level
/// `0` always addresses the global manager in the `extend*` operations; the
/// per-level map never holds a `0` key, so [`Self::by_ms_level`]`(0)` is
/// `None` — the source throws `Exception::InvalidValue` there — and the global
/// ranges are read through [`Self::global`].
///
/// The source header declares nothing beyond the members below: there is no
/// `clear(dim)` other than the inherited [`RangeManager::clear_dim`] (which the
/// source applies to the global ranges only), no `insert` and no `getRange`.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectrumRangeManager {
    global: RangeManager,
    ms_level_ranges: BTreeMap<u32, RangeManager>,
}

impl Default for SpectrumRangeManager {
    fn default() -> Self {
        Self::new()
    }
}

impl SpectrumRangeManager {
    /// An empty manager with no registered MS level.
    pub fn new() -> Self {
        Self {
            global: RangeManager::spectrum_manager(),
            ms_level_ranges: BTreeMap::new(),
        }
    }

    /// The global ranges over all MS levels (the source base-class subobject).
    pub fn global(&self) -> &RangeManager {
        &self.global
    }

    /// The mutable global ranges; use it for the inherited operations such as
    /// [`RangeManager::clear_dim`], which the source applies to the global
    /// ranges only.
    pub fn global_mut(&mut self) -> &mut RangeManager {
        &mut self.global
    }

    /// Clear the global and every MS-level-specific range (source `clearRanges`).
    pub fn clear_ranges(&mut self) {
        self.global.clear_ranges();
        self.ms_level_ranges.clear();
    }

    fn level_mut(&mut self, ms_level: u32) -> &mut RangeManager {
        if ms_level == 0 {
            &mut self.global
        } else {
            self.ms_level_ranges
                .entry(ms_level)
                .or_insert_with(RangeManager::spectrum_manager)
        }
    }

    /// Extend the ranges of `ms_level` (`0` for the global ranges) with `other`.
    /// A level not seen before is registered.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when `other` has no dimension in
    /// common, as the base `extend` throws. The source registers the level
    /// before throwing; this port registers it only on success.
    pub fn extend(&mut self, other: &RangeManager, ms_level: u32) -> Result<()> {
        let mut target = if ms_level == 0 {
            self.global
        } else {
            self.ms_level_ranges
                .get(&ms_level)
                .copied()
                .unwrap_or_else(RangeManager::spectrum_manager)
        };
        target.extend(other)?;
        *self.level_mut(ms_level) = target;
        Ok(())
    }

    /// The ranges of `ms_level`, or `None` when no ranges were recorded for it
    /// (source `byMSLevel`, which throws `Exception::InvalidValue`). Level `0`
    /// is never present; read [`Self::global`] instead.
    pub fn by_ms_level(&self, ms_level: u32) -> Option<&RangeManager> {
        self.ms_level_ranges.get(&ms_level)
    }

    /// Every MS level for which specific ranges exist, ascending (source
    /// `getMSLevels`). Global extends (level `0`) register no level.
    pub fn ms_levels(&self) -> BTreeSet<u32> {
        self.ms_level_ranges.keys().copied().collect()
    }

    /// Extend the RT range of `ms_level` (`0` for the global range) with `rt`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `rt` is not finite; nothing is
    /// registered on error.
    pub fn extend_rt(&mut self, rt: f64, ms_level: u32) -> Result<()> {
        finite(rt, "spectrum retention time")?;
        self.level_mut(ms_level).extend_rt(rt)
    }

    /// Extend the m/z range of `ms_level` (`0` for the global range) with `mz`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `mz` is not finite; nothing is
    /// registered on error.
    pub fn extend_mz(&mut self, mz: f64, ms_level: u32) -> Result<()> {
        finite(mz, "spectrum m/z")?;
        self.level_mut(ms_level).extend_mz(mz)
    }

    /// Extend the ranges of `ms_level` (`0` for the global ranges) with the
    /// current ranges of `spectrum` (source `extendUnsafe(const MSSpectrum&, UInt)`).
    ///
    /// The source reads the spectrum's cached range; this port computes it with
    /// [`MSSpectrum::range_manager`]. The spectrum's own RT is not part of its
    /// range, exactly as in the source, so callers add it with [`Self::extend_rt`].
    ///
    /// # Errors
    ///
    /// As [`MSSpectrum::range_manager`]; nothing is registered on error.
    pub fn extend_spectrum(&mut self, spectrum: &MSSpectrum, ms_level: u32) -> Result<()> {
        let ranges = spectrum.range_manager()?;
        self.level_mut(ms_level).extend_unsafe(&ranges);
        Ok(())
    }
}

/// Index of the first float data array that carries ion mobility, following
/// the source `IMDataArrayUtils::getIMUnit`: an exact PSI-MS child of
/// `MS:1002893 ! ion mobility array` at the pinned CV, or one of the
/// UserParam fallback prefixes (`"Ion Mobility"`, `"inverse reduced ion
/// mobility"`, `"mean inverse reduced ion mobility array"`).
fn ion_mobility_array_index(spectrum: &MSSpectrum) -> Option<usize> {
    const CV_CHILDREN: [&str; 9] = [
        "mean ion mobility drift time array",
        "mean ion mobility array",
        "mean inverse reduced ion mobility array",
        "raw ion mobility array",
        "raw inverse reduced ion mobility array",
        "raw ion mobility drift time array",
        "deconvoluted ion mobility array",
        "deconvoluted inverse reduced ion mobility array",
        "deconvoluted ion mobility drift time array",
    ];
    use crate::constants::user_param;
    spectrum.float_data_arrays.iter().position(|array| {
        let name = array.name.as_str();
        CV_CHILDREN.contains(&name)
            || name.starts_with(user_param::MEAN_INVERSE_REDUCED_ION_MOBILITY_ARRAY)
            || name.starts_with(user_param::INVERSE_REDUCED_ION_MOBILITY)
            || name.starts_with(user_param::ION_MOBILITY)
    })
}

/// The source `IMTypes::DRIFTTIME_NOT_SET`.
const DRIFTTIME_NOT_SET: f64 = -1.0;

fn checked_add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or_else(limit)
}

fn spectrum_items(spectrum: &MSSpectrum) -> Result<usize> {
    let mobility = ion_mobility_array_index(spectrum)
        .map_or(0, |index| spectrum.float_data_arrays[index].data.len());
    let items = checked_add(checked_add(spectrum.peaks.len(), mobility)?, 1)?;
    if items > RangeManager::MAX_ITEMS {
        return Err(limit());
    }
    Ok(items)
}

impl MSSpectrum {
    /// The current m/z, intensity and mobility ranges, computed on demand
    /// (source `updateRanges()` followed by `getRange()`).
    ///
    /// m/z and intensity come from the peaks. Mobility comes from the first
    /// ion-mobility float data array when the spectrum represents an IM frame
    /// (source `containsIMData()`), otherwise from the scalar
    /// [`MSSpectrum::drift_time`] when it is not the unset sentinel `-1`
    /// (`MSSpectrum.cpp:592-604`). The spectrum's own RT is not part of its
    /// range; an experiment adds it (`MSExperiment.cpp:696`).
    ///
    /// # Errors
    ///
    /// As [`MSSpectrum::validate`]; also [`Error::InvalidValue`] when the
    /// drift time or an ion-mobility array value is not finite, or when the
    /// peak and mobility count exceeds [`RangeManager::MAX_ITEMS`].
    pub fn range_manager(&self) -> Result<RangeManager> {
        spectrum_items(self)?;
        self.validate()?;
        let mut ranges = RangeManager::spectrum();
        for peak in &self.peaks {
            ranges.extend_mz(peak.mz)?;
            ranges.extend_intensity(f64::from(peak.intensity))?;
        }
        if let Some(index) = ion_mobility_array_index(self) {
            for &mobility in &self.float_data_arrays[index].data {
                ranges.extend_mobility(f64::from(mobility))?;
            }
        } else if self.drift_time != DRIFTTIME_NOT_SET {
            ranges.extend_mobility(self.drift_time)?;
        }
        Ok(ranges)
    }
}

impl MSChromatogram {
    /// The current RT and intensity ranges of the points, computed on demand
    /// (source `updateRanges()` followed by `getRange()`). The product m/z is
    /// not part of a chromatogram's range; an experiment adds it
    /// (`MSExperiment.cpp:713`).
    ///
    /// # Errors
    ///
    /// As [`MSChromatogram::validate`]; also [`Error::InvalidValue`] when the
    /// point count exceeds [`RangeManager::MAX_ITEMS`].
    pub fn range_manager(&self) -> Result<RangeManager> {
        if self.peaks.len() > RangeManager::MAX_ITEMS {
            return Err(limit());
        }
        self.validate()?;
        let mut ranges = RangeManager::chromatogram();
        for peak in &self.peaks {
            ranges.extend_rt(peak.rt)?;
            ranges.extend_intensity(f64::from(peak.intensity))?;
        }
        Ok(ranges)
    }
}

impl Mobilogram {
    /// The current mobility and intensity ranges of the peaks, computed on
    /// demand (source `updateRanges()` followed by `getRange()`,
    /// `Mobilogram.cpp:43-48`). The scalar RT is not part of the range.
    ///
    /// # Errors
    ///
    /// As [`Mobilogram::validate`]; also [`Error::InvalidValue`] when the peak
    /// count exceeds [`RangeManager::MAX_ITEMS`].
    pub fn range_manager(&self) -> Result<RangeManager> {
        if self.peaks.len() > RangeManager::MAX_ITEMS {
            return Err(limit());
        }
        self.validate()?;
        let mut ranges = RangeManager::mobilogram();
        for peak in &self.peaks {
            ranges.extend_mobility(peak.mobility)?;
            ranges.extend_intensity(f64::from(peak.intensity))?;
        }
        Ok(ranges)
    }
}

impl MSExperiment {
    /// Whole-input item and byte ceilings for one experiment-level range
    /// computation, checked before any allocation.
    fn preflight_ranges(&self, spectra: bool, chromatograms: bool) -> Result<()> {
        let mut items = 0usize;
        if spectra {
            let per_level = std::mem::size_of::<(u32, RangeManager)>()
                .checked_mul(2)
                .and_then(|n| n.checked_add(64))
                .ok_or_else(limit)?;
            let mut levels = BTreeSet::new();
            items = checked_add(items, self.spectra.len())?;
            for spectrum in &self.spectra {
                items = checked_add(items, spectrum_items(spectrum)?)?;
                if items > RangeManager::MAX_ITEMS {
                    return Err(limit());
                }
                if spectrum.ms_level != 0 && levels.insert(spectrum.ms_level) {
                    let bytes = levels.len().checked_mul(per_level).ok_or_else(limit)?;
                    if bytes > RangeManager::MAX_BYTES {
                        return Err(limit());
                    }
                }
            }
        }
        if chromatograms {
            items = checked_add(items, self.chromatograms.len())?;
            for chromatogram in &self.chromatograms {
                items = checked_add(items, chromatogram.peaks.len())?;
                if items > RangeManager::MAX_ITEMS {
                    return Err(limit());
                }
            }
        }
        Ok(())
    }

    /// The current spectrum ranges, global and per MS level (source
    /// `updateRanges()` followed by `spectrumRanges()`).
    ///
    /// Every spectrum extends the global ranges and the ranges of its own MS
    /// level with its m/z, intensity and mobility ([`MSSpectrum::range_manager`])
    /// and with its RT, so a spectrum without peaks still contributes its RT
    /// (`MSExperiment.cpp:696`). A spectrum with MS level `0` extends the
    /// global ranges twice and registers no level, because level `0` addresses
    /// the global ranges in the source. The source's `updateRanges()` also
    /// refreshes each spectrum's own cache as a side effect (`MSExperiment.cpp:691`);
    /// on-demand computation makes that implicit.
    ///
    /// # Errors
    ///
    /// As [`MSSpectrum::range_manager`] for any spectrum, and
    /// [`Error::InvalidValue`] when the visited items exceed
    /// [`RangeManager::MAX_ITEMS`] or the per-level managers would exceed
    /// [`RangeManager::MAX_BYTES`]. Nothing is allocated before the preflight
    /// passes.
    pub fn spectrum_range_manager(&self) -> Result<SpectrumRangeManager> {
        self.preflight_ranges(true, false)?;
        let mut manager = SpectrumRangeManager::new();
        for spectrum in &self.spectra {
            let ranges = spectrum.range_manager()?;
            manager.global_mut().extend_unsafe(&ranges);
            manager.extend_rt(spectrum.rt, 0)?;
            manager.level_mut(spectrum.ms_level).extend_unsafe(&ranges);
            manager.extend_rt(spectrum.rt, spectrum.ms_level)?;
        }
        Ok(manager)
    }

    /// The current chromatogram ranges over RT, intensity and m/z (source
    /// `updateRanges()` followed by `chromatogramRanges()`, a
    /// `ChromatogramRangeManager`).
    ///
    /// RT and intensity come from every chromatogram's points
    /// ([`MSChromatogram::range_manager`]); m/z is each chromatogram's
    /// `getMZ()`, which is its **product** m/z (`MSChromatogram.cpp:81-84`),
    /// added even when the chromatogram has no points (`MSExperiment.cpp:713`).
    /// A chromatogram without a product therefore contributes m/z `0`, as in
    /// the source.
    ///
    /// # Errors
    ///
    /// As [`MSChromatogram::range_manager`] for any chromatogram, and
    /// [`Error::InvalidValue`] when the visited items exceed
    /// [`RangeManager::MAX_ITEMS`].
    pub fn chromatogram_range_manager(&self) -> Result<RangeManager> {
        self.preflight_ranges(false, true)?;
        self.chromatogram_ranges_unchecked()
    }

    fn chromatogram_ranges_unchecked(&self) -> Result<RangeManager> {
        let mut manager = RangeManager::chromatogram_manager();
        for chromatogram in &self.chromatograms {
            manager.extend(&chromatogram.range_manager()?)?;
            manager.extend_mz(chromatogram.product.mz)?;
        }
        Ok(manager)
    }

    /// The combined ranges of all spectra and chromatograms over RT, m/z,
    /// intensity and mobility (source `updateRanges()` followed by
    /// `combinedRanges()`).
    ///
    /// The global spectrum ranges are merged first, then the chromatogram
    /// ranges (`MSExperiment.cpp:718-719`), so the combined RT and intensity
    /// include chromatogram points and the combined m/z includes every
    /// chromatogram's product m/z. On equal endpoints the spectrum value is
    /// kept, which preserves the sign of a zero.
    ///
    /// # Errors
    ///
    /// As [`Self::spectrum_range_manager`] and [`Self::chromatogram_range_manager`].
    pub fn combined_range_manager(&self) -> Result<RangeManager> {
        self.preflight_ranges(true, true)?;
        let spectra = self.spectrum_range_manager()?;
        let chromatograms = self.chromatogram_ranges_unchecked()?;
        let mut combined = RangeManager::experiment();
        combined.extend_unsafe(spectra.global());
        combined.extend_unsafe(&chromatograms);
        Ok(combined)
    }
}
