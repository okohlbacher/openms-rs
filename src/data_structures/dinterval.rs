// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Closed `D`-dimensional intervals delimited by two [`DPosition`]s.
//!
//! Native equivalent of `DATASTRUCTURES/DIntervalBase.h` (header-only
//! template in namespace `Internal`; `DIntervalBase.cpp` is empty). The
//! half-open [`DRange`](crate::data_structures::DRange) builds on this type. Work is over two
//! fixed-size arrays, so no bounded-work ceiling applies. See
//! `docs/DPOSITION_SUPPORT.md`.

use super::DPosition;
use crate::{Error, Result};
use std::{
    fmt,
    hash::{Hash, Hasher},
    ops::{Add, AddAssign, Sub, SubAssign},
};

/// A base class for `D`-dimensional intervals.
///
/// See `DIntervalBase` for a closed interval and [`DRange`](crate::data_structures::DRange)
/// for a half-open interval class.
///
/// **Invariant** (source `@invariant`): all methods maintain
/// `min_position()[x] <= max_position()[x]` for every dimension. The one
/// exception is the [`empty`](Self::empty) sentinel, whose minimum is
/// `f64::MAX` and whose maximum is `f64::MIN`; the source constructs it
/// through a protected non-normalizing constructor and every mutator that
/// takes a whole position restores the invariant from it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DIntervalBase<const D: usize> {
    /// Lower-left point (source `min_`).
    min: DPosition<D>,
    /// Upper-right point (source `max_`).
    max: DPosition<D>,
}

/// One-dimensional closed interval.
pub type DIntervalBase1 = DIntervalBase<1>;
/// Two-dimensional closed interval.
pub type DIntervalBase2 = DIntervalBase<2>;

impl<const D: usize> Default for DIntervalBase<D> {
    /// The empty interval (source default constructor, documented as "corners
    /// at infinity" although the corners are the finite extrema
    /// `f64::MAX`/`f64::MIN`).
    fn default() -> Self {
        Self::empty()
    }
}

impl<const D: usize> DIntervalBase<D> {
    /// Number of dimensions (source `DIMENSION` enumerator).
    pub const DIMENSION: usize = D;

    /// The empty instance (source static `empty`): minimum
    /// [`DPosition::max_positive`], maximum [`DPosition::min_negative`].
    pub const fn empty() -> Self {
        Self {
            min: DPosition::max_positive(),
            max: DPosition::min_negative(),
        }
    }

    /// The instance with all positions zero (source static `zero`).
    pub const fn zero() -> Self {
        Self {
            min: DPosition::zero(),
            max: DPosition::zero(),
        }
    }

    /// Interval from two corners (source `DIntervalBase(minimum, maximum)`).
    ///
    /// The corners are normalized per dimension: wherever `minimum[i] >
    /// maximum[i]` the two coordinates are swapped, so the invariant holds
    /// for any input order. NaN coordinates are never swapped.
    pub fn new(minimum: DPosition<D>, maximum: DPosition<D>) -> Self {
        let mut interval = Self {
            min: minimum,
            max: maximum,
        };
        interval.normalize();
        interval
    }

    /// Accessor to the minimum position (source `minPosition()`).
    pub const fn min_position(&self) -> &DPosition<D> {
        &self.min
    }

    /// Accessor to the maximum position (source `maxPosition()`).
    pub const fn max_position(&self) -> &DPosition<D> {
        &self.max
    }

    /// Sets the minimum position (source `setMin`).
    ///
    /// Source `@note`: the minimum given here is what `min_position()`
    /// returns afterwards; where necessary `max_position()` is raised to it,
    /// dimension by dimension.
    pub fn set_min(&mut self, position: DPosition<D>) {
        self.min = position;
        for i in 0..D {
            if self.min[i] > self.max[i] {
                self.max[i] = self.min[i];
            }
        }
    }

    /// Sets the maximum position (source `setMax`).
    ///
    /// Source `@note`: the maximum given here is what `max_position()`
    /// returns afterwards; where necessary `min_position()` is lowered to it,
    /// dimension by dimension.
    pub fn set_max(&mut self, position: DPosition<D>) {
        self.max = position;
        for i in 0..D {
            if self.min[i] > self.max[i] {
                self.min[i] = self.max[i];
            }
        }
    }

    /// Sets both corners and normalizes them per dimension (source `setMinMax`).
    pub fn set_min_max(&mut self, min: DPosition<D>, max: DPosition<D>) {
        self.min = min;
        self.max = max;
        self.normalize();
    }

    /// Assignment from an interval of a different dimension (source
    /// `assign<D2>`): only dimensions `0..min(D, D2)` are copied, the others
    /// keep their values. No normalization is applied, as in the source.
    pub fn assign<const D2: usize>(&mut self, rhs: &DIntervalBase<D2>) {
        for i in 0..D.min(D2) {
            self.min[i] = rhs.min[i];
            self.max[i] = rhs.max[i];
        }
    }

    /// Makes the interval empty (source `clear()`), i.e. equal to [`empty`](Self::empty).
    pub fn clear(&mut self) {
        *self = Self::empty();
    }

    /// Whether the interval is completely empty, i.e. cleared or default
    /// constructed (source `isEmpty()`). If `min == max` the interval is
    /// **not** empty. Compares against the sentinel, so a NaN corner is never
    /// empty.
    pub fn is_empty(&self) -> bool {
        *self == Self::empty()
    }

    /// Whether dimension `dim` is empty (source `isEmpty(UInt dim)`), i.e.
    /// its coordinates are the sentinel pair `f64::MAX`/`f64::MIN`. If
    /// `min == max` in that dimension the interval is **not** empty.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `dim >= D`; the source indexes
    /// its arrays without a check.
    pub fn is_empty_dim(&self, dim: usize) -> Result<bool> {
        Self::check_dim(dim)?;
        Ok(self.min[dim] == f64::MAX && self.max[dim] == f64::MIN)
    }

    /// Sets the interval of a single dimension from a one-dimensional
    /// interval (source `setDimMinMax`). The coordinates are copied as given,
    /// without normalization against the other dimensions.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `dim >= D`; the source indexes
    /// its arrays without a check.
    pub fn set_dim_min_max(&mut self, dim: usize, min_max: &DIntervalBase<1>) -> Result<()> {
        Self::check_dim(dim)?;
        self.min[dim] = min_max.min[0];
        self.max[dim] = min_max.max[0];
        Ok(())
    }

    /// The center of the interval, `(min + max) / 2` (source `center()`).
    /// For the empty sentinel this is zero, because `f64::MAX + f64::MIN`
    /// is `0`.
    pub fn center(&self) -> DPosition<D> {
        let mut center = self.min;
        center += self.max;
        center /= 2.;
        center
    }

    /// The diagonal of the area, i.e. `max - min` (source `diagonal()`).
    /// For the empty sentinel every coordinate is `-inf`.
    pub fn diagonal(&self) -> DPosition<D> {
        self.max - self.min
    }

    /// Minimum of dimension zero (source `minX()`).
    ///
    /// The 2D convenience accessors are defined for every `D` in the source
    /// and index out of bounds for smaller dimensions; here instantiating
    /// them with too few dimensions is a compile-time error.
    pub fn min_x(&self) -> f64 {
        const { assert!(D >= 1, "minX() needs at least one dimension") }
        self.min[0]
    }

    /// Minimum of dimension one (source `minY()`); compile-time error for `D < 2`.
    pub fn min_y(&self) -> f64 {
        const { assert!(D >= 2, "minY() needs at least two dimensions") }
        self.min[1]
    }

    /// Maximum of dimension zero (source `maxX()`); compile-time error for `D < 1`.
    pub fn max_x(&self) -> f64 {
        const { assert!(D >= 1, "maxX() needs at least one dimension") }
        self.max[0]
    }

    /// Maximum of dimension one (source `maxY()`); compile-time error for `D < 2`.
    pub fn max_y(&self) -> f64 {
        const { assert!(D >= 2, "maxY() needs at least two dimensions") }
        self.max[1]
    }

    /// Sets the minimum of dimension zero, raising the maximum to it if
    /// necessary (source `setMinX`); compile-time error for `D < 1`.
    pub fn set_min_x(&mut self, c: f64) {
        const { assert!(D >= 1, "setMinX() needs at least one dimension") }
        self.min[0] = c;
        if self.min[0] > self.max[0] {
            self.max[0] = self.min[0];
        }
    }

    /// Sets the minimum of dimension one, raising the maximum to it if
    /// necessary (source `setMinY`); compile-time error for `D < 2`.
    pub fn set_min_y(&mut self, c: f64) {
        const { assert!(D >= 2, "setMinY() needs at least two dimensions") }
        self.min[1] = c;
        if self.min[1] > self.max[1] {
            self.max[1] = self.min[1];
        }
    }

    /// Sets the maximum of dimension zero, lowering the minimum to it if
    /// necessary (source `setMaxX`); compile-time error for `D < 1`.
    pub fn set_max_x(&mut self, c: f64) {
        const { assert!(D >= 1, "setMaxX() needs at least one dimension") }
        self.max[0] = c;
        if self.min[0] > self.max[0] {
            self.min[0] = self.max[0];
        }
    }

    /// Sets the maximum of dimension one, lowering the minimum to it if
    /// necessary (source `setMaxY`); compile-time error for `D < 2`.
    pub fn set_max_y(&mut self, c: f64) {
        const { assert!(D >= 2, "setMaxY() needs at least two dimensions") }
        self.max[1] = c;
        if self.min[1] > self.max[1] {
            self.min[1] = self.max[1];
        }
    }

    /// Width of the area, i.e. the difference of dimension zero (source
    /// `width()`); compile-time error for `D < 1`.
    pub fn width(&self) -> f64 {
        const { assert!(D >= 1, "width() needs at least one dimension") }
        self.max[0] - self.min[0]
    }

    /// Height of the area, i.e. the difference of dimension one (source
    /// `height()`); compile-time error for `D < 2`.
    pub fn height(&self) -> f64 {
        const { assert!(D >= 2, "height() needs at least two dimensions") }
        self.max[1] - self.min[1]
    }

    /// Mutable corners for the half-open range type (source `min_`/`max_`
    /// are protected and `using`-imported by `DRange`).
    pub(super) fn corners_mut(&mut self) -> (&mut DPosition<D>, &mut DPosition<D>) {
        (&mut self.min, &mut self.max)
    }

    /// Source `normalize_()`: swap the coordinates of every dimension whose
    /// minimum exceeds its maximum.
    fn normalize(&mut self) {
        for i in 0..D {
            if self.min[i] > self.max[i] {
                let (min, max) = (self.min[i], self.max[i]);
                self.min[i] = max;
                self.max[i] = min;
            }
        }
    }

    fn check_dim(dim: usize) -> Result<()> {
        if dim < D {
            Ok(())
        } else {
            Err(Error::InvalidValue(format!(
                "dimension {dim} is out of range for a {D}-dimensional interval"
            )))
        }
    }
}

impl<const D: usize> Add<DPosition<D>> for DIntervalBase<D> {
    type Output = Self;

    /// Translates both corners by `point` (source `operator+`).
    fn add(mut self, point: DPosition<D>) -> Self {
        self += point;
        self
    }
}

impl<const D: usize> AddAssign<DPosition<D>> for DIntervalBase<D> {
    /// Translates both corners by `point` in place (source `operator+=`).
    fn add_assign(&mut self, point: DPosition<D>) {
        self.min += point;
        self.max += point;
    }
}

impl<const D: usize> Sub<DPosition<D>> for DIntervalBase<D> {
    type Output = Self;

    /// Translates both corners by `-point` (source `operator-`).
    fn sub(mut self, point: DPosition<D>) -> Self {
        self -= point;
        self
    }
}

impl<const D: usize> SubAssign<DPosition<D>> for DIntervalBase<D> {
    /// Translates both corners by `-point` in place (source `operator-=`).
    fn sub_assign(&mut self, point: DPosition<D>) {
        self.min -= point;
        self.max -= point;
    }
}

impl<const D: usize> Hash for DIntervalBase<D> {
    /// Hashes the minimum then the maximum corner with [`DPosition`]'s
    /// signed-zero normalization. The source declares no `std::hash` for
    /// `DIntervalBase`; this is the same ordering its `DRange` hash uses.
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.min.hash(state);
        self.max.hash(state);
    }
}

impl<const D: usize> fmt::Display for DIntervalBase<D> {
    /// Source `operator<<`: four lines, each terminated by a newline,
    /// `--DIntervalBase BEGIN--`, `MIN --> <min>`, `MAX --> <max>`,
    /// `--DIntervalBase END--`. A precision option is forwarded to the
    /// corners.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "--DIntervalBase BEGIN--")?;
        match f.precision() {
            Some(p) => {
                writeln!(f, "MIN --> {:.*}", p, self.min)?;
                writeln!(f, "MAX --> {:.*}", p, self.max)?;
            }
            None => {
                writeln!(f, "MIN --> {}", self.min)?;
                writeln!(f, "MAX --> {}", self.max)?;
            }
        }
        writeln!(f, "--DIntervalBase END--")
    }
}
