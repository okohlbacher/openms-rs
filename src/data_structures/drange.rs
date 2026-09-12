// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Half-open `D`-dimensional ranges.
//!
//! Native equivalent of `DATASTRUCTURES/DRange.h` (header-only template;
//! `DRange.cpp` only instantiates default objects). The source class derives
//! from `Internal::DIntervalBase`; this port composes a [`DIntervalBase`] and
//! delegates its public members explicitly, so every inherited operation is
//! visible on [`DRange`] and convertible in both directions. Work is over two
//! fixed-size arrays, so no bounded-work ceiling applies. See
//! `docs/DPOSITION_SUPPORT.md`.

use super::{DIntervalBase, DPosition};
use crate::{Error, Result};
use std::{
    fmt,
    hash::{Hash, Hasher},
    ops::{Add, AddAssign, Sub, SubAssign},
};

/// Kinds of intersection between two ranges (source `DRange::DRangeIntersection`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum DRangeIntersection {
    /// No intersection.
    Disjoint,
    /// Intersection.
    Intersects,
    /// One contains the other.
    Inside,
}

/// A `D`-dimensional half-open interval.
///
/// This class describes a range in `D`-dimensional space delimited by two
/// points (i.e. a `D`-dimensional hyper-rectangle). The two points define the
/// lower left and the upper right corner in 2D and analogous points in higher
/// dimensions.
///
/// A range is a pair of positions in `D`-space represented by [`DPosition`].
/// The two limiting points are accessed as [`min_position`](Self::min_position)
/// and [`max_position`](Self::max_position).
///
/// A range denotes a semi-open interval: the lower coordinate of each
/// dimension is part of the range, the higher coordinate is not.
///
/// The invariant `min_position()[i] <= max_position()[i]` and the
/// [`empty`](Self::empty) sentinel are those of [`DIntervalBase`].
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct DRange<const D: usize> {
    base: DIntervalBase<D>,
}

/// One-dimensional half-open range.
pub type DRange1 = DRange<1>;
/// Two-dimensional half-open range.
pub type DRange2 = DRange<2>;

impl<const D: usize> DRange<D> {
    /// Number of dimensions (source `DIMENSION` enumerator).
    pub const DIMENSION: usize = D;

    /// Range from two points (source `DRange(lower, upper)`); the corners are
    /// normalized per dimension like [`DIntervalBase::new`].
    pub fn new(lower: DPosition<D>, upper: DPosition<D>) -> Self {
        Self {
            base: DIntervalBase::new(lower, upper),
        }
    }

    /// The empty range (source inherited static `empty`). This is also what
    /// [`Default`] produces: the source default constructor's comment says
    /// "all coordinates zero", but its body calls the base default
    /// constructor, which yields the empty sentinel.
    pub const fn empty() -> Self {
        Self {
            base: DIntervalBase::empty(),
        }
    }

    /// The range with all positions zero (source inherited static `zero`).
    pub const fn zero() -> Self {
        Self {
            base: DIntervalBase::zero(),
        }
    }

    /// The closed-interval view of the same corners (source base-class
    /// subobject).
    pub const fn base(&self) -> &DIntervalBase<D> {
        &self.base
    }

    /// Mutable closed-interval view of the same corners.
    pub fn base_mut(&mut self) -> &mut DIntervalBase<D> {
        &mut self.base
    }

    /// Accessor to the minimum position (source `minPosition()`).
    pub const fn min_position(&self) -> &DPosition<D> {
        self.base.min_position()
    }

    /// Accessor to the maximum position (source `maxPosition()`).
    pub const fn max_position(&self) -> &DPosition<D> {
        self.base.max_position()
    }

    /// Sets the minimum position, raising the maximum where necessary; see
    /// [`DIntervalBase::set_min`].
    pub fn set_min(&mut self, position: DPosition<D>) {
        self.base.set_min(position);
    }

    /// Sets the maximum position, lowering the minimum where necessary; see
    /// [`DIntervalBase::set_max`].
    pub fn set_max(&mut self, position: DPosition<D>) {
        self.base.set_max(position);
    }

    /// Sets both corners with normalization; see [`DIntervalBase::set_min_max`].
    pub fn set_min_max(&mut self, min: DPosition<D>, max: DPosition<D>) {
        self.base.set_min_max(min, max);
    }

    /// Copies dimensions `0..min(D, D2)` from another interval; see
    /// [`DIntervalBase::assign`].
    pub fn assign<const D2: usize>(&mut self, rhs: &DIntervalBase<D2>) {
        self.base.assign(rhs);
    }

    /// Makes the range empty; see [`DIntervalBase::clear`].
    pub fn clear(&mut self) {
        self.base.clear();
    }

    /// Whether the range is the empty sentinel; see [`DIntervalBase::is_empty`].
    pub fn is_empty(&self) -> bool {
        self.base.is_empty()
    }

    /// Whether dimension `dim` is the empty sentinel pair; see
    /// [`DIntervalBase::is_empty_dim`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `dim >= D`.
    pub fn is_empty_dim(&self, dim: usize) -> Result<bool> {
        self.base.is_empty_dim(dim)
    }

    /// Sets a single dimension from a one-dimensional interval; see
    /// [`DIntervalBase::set_dim_min_max`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `dim >= D`.
    pub fn set_dim_min_max(&mut self, dim: usize, min_max: &DIntervalBase<1>) -> Result<()> {
        self.base.set_dim_min_max(dim, min_max)
    }

    /// The center `(min + max) / 2`; see [`DIntervalBase::center`].
    pub fn center(&self) -> DPosition<D> {
        self.base.center()
    }

    /// The diagonal `max - min`; see [`DIntervalBase::diagonal`].
    pub fn diagonal(&self) -> DPosition<D> {
        self.base.diagonal()
    }

    /// Minimum of dimension zero; see [`DIntervalBase::min_x`].
    pub fn min_x(&self) -> f64 {
        self.base.min_x()
    }

    /// Minimum of dimension one; see [`DIntervalBase::min_y`].
    pub fn min_y(&self) -> f64 {
        self.base.min_y()
    }

    /// Maximum of dimension zero; see [`DIntervalBase::max_x`].
    pub fn max_x(&self) -> f64 {
        self.base.max_x()
    }

    /// Maximum of dimension one; see [`DIntervalBase::max_y`].
    pub fn max_y(&self) -> f64 {
        self.base.max_y()
    }

    /// Sets the minimum of dimension zero; see [`DIntervalBase::set_min_x`].
    pub fn set_min_x(&mut self, c: f64) {
        self.base.set_min_x(c);
    }

    /// Sets the minimum of dimension one; see [`DIntervalBase::set_min_y`].
    pub fn set_min_y(&mut self, c: f64) {
        self.base.set_min_y(c);
    }

    /// Sets the maximum of dimension zero; see [`DIntervalBase::set_max_x`].
    pub fn set_max_x(&mut self, c: f64) {
        self.base.set_max_x(c);
    }

    /// Sets the maximum of dimension one; see [`DIntervalBase::set_max_y`].
    pub fn set_max_y(&mut self, c: f64) {
        self.base.set_max_y(c);
    }

    /// Width, the difference of dimension zero; see [`DIntervalBase::width`].
    pub fn width(&self) -> f64 {
        self.base.width()
    }

    /// Height, the difference of dimension one; see [`DIntervalBase::height`].
    pub fn height(&self) -> f64 {
        self.base.height()
    }

    /// Checks whether this range (half-open interval!) contains `position`
    /// (source `encloses(const PositionType&)`): every coordinate must be
    /// `>= min` and `< max`. A NaN coordinate fails neither comparison and is
    /// therefore enclosed, as in the source.
    pub fn encloses(&self, position: &DPosition<D>) -> bool {
        let (min, max) = (self.base.min_position(), self.base.max_position());
        for i in 0..D {
            if position[i] < min[i] {
                return false;
            }
            if position[i] >= max[i] {
                return false;
            }
        }
        true
    }

    /// The smallest range containing this range and `other` (source
    /// `united`): the per-dimension minimum of the minima and maximum of the
    /// maxima, passed through the normalizing [`set_min_max`](Self::set_min_max).
    ///
    /// Uniting a non-empty range with the [`empty`](Self::empty) sentinel
    /// returns that range. Uniting two empty ranges does **not** stay empty:
    /// the sentinel's inverted corners `(f64::MAX, f64::MIN)` are swapped by
    /// the normalization into the all-encompassing range `[f64::MIN,
    /// f64::MAX]`. This transcribes the source; it is recorded as a C++ issue
    /// candidate rather than corrected, because callers may depend on it.
    pub fn united(&self, other: &Self) -> Self {
        let (min, max) = (self.base.min_position(), self.base.max_position());
        let (other_min, other_max) = (other.base.min_position(), other.base.max_position());
        let mut united_min = DPosition::<D>::zero();
        let mut united_max = DPosition::<D>::zero();
        for i in 0..D {
            united_min[i] = if min[i] < other_min[i] {
                min[i]
            } else {
                other_min[i]
            };
            united_max[i] = if max[i] > other_max[i] {
                max[i]
            } else {
                other_max[i]
            };
        }
        let mut united = Self::empty();
        united.set_min_max(united_min, united_max);
        united
    }

    /// Checks how this range intersects with another `range` (source
    /// `intersects`).
    ///
    /// The decision follows the source exactly: if `range`'s minimum is
    /// [enclosed](Self::encloses) the result is [`Inside`](DRangeIntersection::Inside)
    /// unless some coordinate of `range`'s maximum exceeds this maximum
    /// ([`Intersects`](DRangeIntersection::Intersects)); otherwise it is
    /// [`Disjoint`](DRangeIntersection::Disjoint) when any `range.min[i] >=
    /// max[i]` or any `range.max[i] <= min[i]`, and `Intersects` else.
    pub fn intersects(&self, range: &Self) -> DRangeIntersection {
        let (min, max) = (self.base.min_position(), self.base.max_position());
        let (range_min, range_max) = (range.base.min_position(), range.base.max_position());
        if self.encloses(range_min) {
            if (0..D).any(|i| range_max[i] > max[i]) {
                return DRangeIntersection::Intersects;
            }
            return DRangeIntersection::Inside;
        }
        if (0..D).any(|i| range_min[i] >= max[i]) || (0..D).any(|i| range_max[i] <= min[i]) {
            return DRangeIntersection::Disjoint;
        }
        DRangeIntersection::Intersects
    }

    /// Whether the areas intersect, i.e. they intersect or one contains the
    /// other (source `isIntersected`). Equivalent to
    /// `intersects(range) != Disjoint`.
    pub fn is_intersected(&self, range: &Self) -> bool {
        self.intersects(range) != DRangeIntersection::Disjoint
    }

    /// Extends the range in all dimensions by a multiplier while keeping the
    /// center (source `extend(double factor)`).
    ///
    /// Examples for `D = 1`: `factor = 1.01` extends the range by 1% in
    /// total, i.e. 0.5% left and right; `factor = 2.0` doubles the total
    /// range, e.g. from `[0, 100]` to `[-50, 150]`. Each side moves by
    /// `(max - min) / 2 * (factor - 1)`; a factor below 1 shrinks the range,
    /// and `0` collapses it onto its center.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `factor` is negative (source
    /// `Exception::InvalidParameter`; the allowed domain is `[0, inf)`) or
    /// not finite. The source does not reject NaN or infinity and would
    /// poison the corners; the check happens before any mutation, so an
    /// error leaves the range unchanged.
    pub fn extend_by_factor(&mut self, factor: f64) -> Result<&mut Self> {
        if !factor.is_finite() || factor < 0. {
            return Err(Error::InvalidValue(format!(
                "DRange::extend(): factor must be finite and not negative, got {factor}"
            )));
        }
        let (min, max) = self.base.corners_mut();
        for i in 0..D {
            let extra = (max[i] - min[i]) / 2. * (factor - 1.);
            min[i] -= extra;
            max[i] += extra;
        }
        Ok(self)
    }

    /// Extends the range in all dimensions by an additive amount while
    /// keeping the center (source `extend(PositionType addition)`).
    ///
    /// Half of `addition[i]` is subtracted from the minimum and added to the
    /// maximum of dimension `i`, so `addition = 0.5` widens a 1D range by 1
    /// in total. `addition` may be negative: the range then shrinks and may
    /// end with `min == max`, but never `min > max` — a dimension that would
    /// invert is reduced to its single center point (`(min + max) / 2`).
    /// The source `@param` remark that "resulting invalid min/max are not
    /// fixed automatically" predates that collapse and is superseded by the
    /// implementation, which this port follows.
    pub fn extend_by(&mut self, mut addition: DPosition<D>) -> &mut Self {
        addition /= 2.;
        let (min, max) = self.base.corners_mut();
        *min -= addition;
        *max += addition;
        for i in 0..D {
            if min[i] > max[i] {
                let center = (min[i] + max[i]) / 2.;
                min[i] = center;
                max[i] = center;
            }
        }
        self
    }

    /// Ensures every dimension spans at least `min_span[i]` (source
    /// `ensureMinSpan`): a dimension narrower than that is widened
    /// symmetrically around its center by the missing amount, wider
    /// dimensions are unchanged.
    pub fn ensure_min_span(&mut self, min_span: DPosition<D>) -> &mut Self {
        let mut extend_by = DPosition::<D>::zero();
        {
            let (min, max) = (self.base.min_position(), self.base.max_position());
            for i in 0..D {
                if max[i] - min[i] < min_span[i] {
                    extend_by[i] = min_span[i] - (max[i] - min[i]);
                }
            }
        }
        self.extend_by(extend_by)
    }

    /// Pulls `point` into the current area (source `pullIn`, whose in/out
    /// parameter becomes the return value here): each coordinate is clamped
    /// to `[min[i], max[i]]`, the maximum included.
    ///
    /// The clamp is transcribed from the source's `std::max(min,
    /// std::min(point, max))`, so a NaN coordinate yields `min[i]`, unlike
    /// `f64::clamp`.
    pub fn pull_in(&self, mut point: DPosition<D>) -> DPosition<D> {
        let (min, max) = (self.base.min_position(), self.base.max_position());
        for i in 0..D {
            let upper = if max[i] < point[i] { max[i] } else { point[i] };
            point[i] = if min[i] < upper { upper } else { min[i] };
        }
        point
    }
}

impl DRange<2> {
    /// Convenient 2D constructor from four coordinates (source
    /// `DRange(minx, miny, maxx, maxy)`, restricted to `D == 2` by
    /// `static_assert`). Corners are normalized per dimension, so
    /// `DRange::xy(2., 3., -2., -3.)` has minimum `(-2, -3)` and maximum
    /// `(2, 3)`.
    pub fn xy(min_x: f64, min_y: f64, max_x: f64, max_y: f64) -> Self {
        Self::new(DPosition::xy(min_x, min_y), DPosition::xy(max_x, max_y))
    }

    /// 2D version of [`encloses`](Self::encloses) for convenience (source
    /// `encloses(x, y)`).
    pub fn encloses_xy(&self, x: f64, y: f64) -> bool {
        self.encloses(&DPosition::xy(x, y))
    }

    /// Swaps the dimensions of 2D data, i.e. x and y coordinates (source
    /// `swapDimensions`, restricted to `D == 2` by `static_assert`).
    pub fn swap_dimensions(&mut self) -> &mut Self {
        let (min, max) = self.base.corners_mut();
        min.coordinates.swap(0, 1);
        max.coordinates.swap(0, 1);
        self
    }
}

impl<const D: usize> From<DIntervalBase<D>> for DRange<D> {
    /// Source "copy constructor for the base class" and `operator=(const
    /// Base&)`: the corners are taken as they are.
    fn from(base: DIntervalBase<D>) -> Self {
        Self { base }
    }
}

impl<const D: usize> From<DRange<D>> for DIntervalBase<D> {
    /// Source derived-to-base conversion.
    fn from(range: DRange<D>) -> Self {
        range.base
    }
}

impl<const D: usize> PartialEq<DIntervalBase<D>> for DRange<D> {
    /// Source `operator==(const Base&)`: corner-wise equality.
    fn eq(&self, other: &DIntervalBase<D>) -> bool {
        self.base == *other
    }
}

impl<const D: usize> PartialEq<DRange<D>> for DIntervalBase<D> {
    /// Symmetric counterpart of the range/base comparison.
    fn eq(&self, other: &DRange<D>) -> bool {
        *self == other.base
    }
}

impl<const D: usize> Add<DPosition<D>> for DRange<D> {
    type Output = Self;

    /// Translates both corners by `point`. The source's inherited `operator+`
    /// returns a `DIntervalBase`; this keeps the range type.
    fn add(mut self, point: DPosition<D>) -> Self {
        self.base += point;
        self
    }
}

impl<const D: usize> AddAssign<DPosition<D>> for DRange<D> {
    /// Translates both corners by `point` in place.
    fn add_assign(&mut self, point: DPosition<D>) {
        self.base += point;
    }
}

impl<const D: usize> Sub<DPosition<D>> for DRange<D> {
    type Output = Self;

    /// Translates both corners by `-point`, keeping the range type.
    fn sub(mut self, point: DPosition<D>) -> Self {
        self.base -= point;
        self
    }
}

impl<const D: usize> SubAssign<DPosition<D>> for DRange<D> {
    /// Translates both corners by `-point` in place.
    fn sub_assign(&mut self, point: DPosition<D>) {
        self.base -= point;
    }
}

impl<const D: usize> Hash for DRange<D> {
    /// Hashes all minimum then all maximum coordinates with signed-zero
    /// normalization (source `std::hash<DRange>` via `hash_float`). No
    /// digest compatibility with the source FNV-1a combination is promised.
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.base.hash(state);
    }
}

impl<const D: usize> fmt::Display for DRange<D> {
    /// Source `operator<<`: four newline-terminated lines, `--DRANGE BEGIN--`,
    /// `MIN --> <min>`, `MAX --> <max>`, `--DRANGE END--`. A precision option
    /// is forwarded to the corners.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "--DRANGE BEGIN--")?;
        match f.precision() {
            Some(p) => {
                writeln!(f, "MIN --> {:.*}", p, self.min_position())?;
                writeln!(f, "MAX --> {:.*}", p, self.max_position())?;
            }
            None => {
                writeln!(f, "MIN --> {}", self.min_position())?;
                writeln!(f, "MAX --> {}", self.max_position())?;
            }
        }
        writeln!(f, "--DRANGE END--")
    }
}
