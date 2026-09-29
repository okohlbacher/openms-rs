// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A `D`-dimensional bounding box: `DATASTRUCTURES/DBoundingBox.h`.
//!
//! A [`DBoundingBox`] denotes a closed interval: upper and lower margins are
//! both contained. It is the source's class template over
//! [`DIntervalBase`], which it derives from; here the base is a private field
//! reached through `Deref`/`DerefMut`, so every base accessor and mutator
//! (`min_position`, `set_min_x`, `center`, ...) is available on the box, as it
//! is on the derived C++ class.
//!
//! `kernel::geometry::BoundingBox2D` is a separate, validated two-dimensional
//! box that the feature hulls use; it rejects inverted and non-finite corners
//! and represents emptiness as `Option`. This type is the source template
//! itself, for any dimension, with the source's sentinel emptiness.
//!
//! # API mapping
//!
//! | Source (`DBoundingBox.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `template <UInt D> class DBoundingBox` | [`DBoundingBox<D>`]; [`DBoundingBox1`], [`DBoundingBox2`] |
//! | `enum {DIMENSION = D}` | [`DBoundingBox::DIMENSION`] |
//! | `typedef Base`, `PositionType`, `CoordinateType` | [`DIntervalBase<D>`], [`DPosition<D>`], `f64` |
//! | `DBoundingBox()` | `Default`, the empty sentinel box |
//! | `DBoundingBox(const DBoundingBox&)`, `operator=(const DBoundingBox&)`, `~DBoundingBox()` | `Clone`, `Copy`, assignment, drop |
//! | `operator=(const Base&)` | `From<DIntervalBase<D>>` |
//! | `DBoundingBox(const PositionType& minimum, const PositionType& maximum)` | [`DBoundingBox::new`] |
//! | `enlarge(const PositionType&)` | [`DBoundingBox::enlarge`] |
//! | `enlarge(CoordinateType x, CoordinateType y)` | [`DBoundingBox::enlarge_xy`] (two dimensions) |
//! | `operator==(const DBoundingBox&)`, `operator==(const Base&)` | `PartialEq`, `PartialEq<DIntervalBase<D>>` |
//! | `encloses(const PositionType&)` | [`DBoundingBox::encloses`] |
//! | `encloses(CoordinateType x, CoordinateType y)` | [`DBoundingBox::encloses_xy`] (two dimensions) |
//! | `intersects(const DBoundingBox&)` | [`DBoundingBox::intersects`] |
//! | `isEmpty()` | [`DBoundingBox::is_empty`] |
//! | `operator<<(std::ostream&, const DBoundingBox<D>&)` | `Display` |
//! | `std::hash<DBoundingBox<D>>` | `Hash`, and [`DBoundingBox::source_hash64`] for the source's value |
//! | `Base::min_`, `Base::max_` (protected, `using`-imported) | the base's accessors through `Deref` |
//!
//! See `docs/DBOUNDING_BOX_SUPPORT.md`.

use std::fmt;
use std::hash::{Hash, Hasher};
use std::ops::{Deref, DerefMut};

use super::{DIntervalBase, DPosition};
use crate::concept::hash_utils::{hash_combine, hash_float};
use crate::concept::precision_wrapper::precision_wrapper;

/// A `D`-dimensional bounding box, a closed interval whose margins are
/// contained: source `DBoundingBox<D>`.
///
/// # Examples
///
/// ```
/// use openms::data_structures::{DBoundingBox2, DPosition2};
///
/// let mut bounds = DBoundingBox2::default();
/// assert!(bounds.is_empty());
/// bounds.enlarge(DPosition2::xy(11.0, 13.0));
/// bounds.enlarge_xy(9.0, 0.0);
/// assert!(bounds.encloses_xy(10.0, 1.0));
/// assert_eq!(bounds.min_position().x(), 9.0);
/// assert!(!bounds.is_empty());
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct DBoundingBox<const D: usize> {
    /// The interval the source class derives from.
    interval: DIntervalBase<D>,
}

/// One-dimensional bounding box, source `DBoundingBox<1>`.
pub type DBoundingBox1 = DBoundingBox<1>;
/// Two-dimensional bounding box, source `DBoundingBox<2>`.
pub type DBoundingBox2 = DBoundingBox<2>;

impl<const D: usize> DBoundingBox<D> {
    /// Number of dimensions: source enumerator `DIMENSION`.
    pub const DIMENSION: usize = D;

    /// The box between two corners: source `DBoundingBox(minimum, maximum)`.
    ///
    /// As the base constructor, the corners are normalised per dimension:
    /// wherever `minimum[i] > maximum[i]` the two coordinates are swapped.
    pub fn new(minimum: DPosition<D>, maximum: DPosition<D>) -> Self {
        Self {
            interval: DIntervalBase::new(minimum, maximum),
        }
    }

    /// The interval this box is, source `Base`.
    pub const fn as_interval(&self) -> &DIntervalBase<D> {
        &self.interval
    }

    /// Enlarges the box so that it contains `position`: source
    /// `enlarge(const PositionType&)`.
    ///
    /// Each dimension's minimum is lowered to the coordinate if it is smaller
    /// and its maximum raised if it is larger, so enlarging the default (empty
    /// sentinel) box by one point gives the degenerate box at that point. A NaN
    /// coordinate compares false both ways and changes nothing, as in the
    /// source.
    pub fn enlarge(&mut self, position: DPosition<D>) {
        let mut min = *self.interval.min_position();
        let mut max = *self.interval.max_position();
        for i in 0..D {
            if position[i] < min[i] {
                min[i] = position[i];
            }
            if position[i] > max[i] {
                max[i] = position[i];
            }
        }
        // min <= max holds per dimension after the loop unless the box was the
        // empty sentinel and `position` had a NaN coordinate; the corners are
        // written as they are, without normalising, as the source writes
        // `min_` and `max_` directly.
        let (min_corner, max_corner) = self.interval.corners_mut();
        *min_corner = min;
        *max_corner = max;
    }

    /// Whether `position` lies inside the box, margins included: source
    /// `encloses(const PositionType&)`.
    ///
    /// The source tests `position[i] < min_[i] || position[i] > max_[i]` for
    /// a rejection, so a NaN coordinate is never rejected and a position with
    /// one is enclosed; that is kept.
    pub fn encloses(&self, position: &DPosition<D>) -> bool {
        let min = self.interval.min_position();
        let max = self.interval.max_position();
        (0..D).all(|i| !(position[i] < min[i] || position[i] > max[i]))
    }

    /// Whether this box intersects `other`, touching margins included: source
    /// `intersects(const DBoundingBox&)`.
    ///
    /// As in the source, a NaN corner never rejects, and the empty sentinel
    /// box intersects nothing but boxes reaching the sentinel values.
    pub fn intersects(&self, other: &Self) -> bool {
        let min = self.interval.min_position();
        let max = self.interval.max_position();
        let other_min = other.interval.min_position();
        let other_max = other.interval.max_position();
        (0..D).all(|i| !(other_min[i] > max[i] || other_max[i] < min[i]))
    }

    /// Whether the box is empty: source `isEmpty()`, true when any dimension
    /// has `max <= min`.
    ///
    /// This differs from the base's [`DIntervalBase::is_empty`], which tests
    /// for the empty sentinel only, and hides it as the source's does: a
    /// zero-width box (`min == max` in some dimension) is empty here but not
    /// for the base. A zero-dimensional box is never empty. A NaN coordinate
    /// compares false, so it does not make a dimension empty.
    pub fn is_empty(&self) -> bool {
        let min = self.interval.min_position();
        let max = self.interval.max_position();
        (0..D).any(|i| max[i] <= min[i])
    }

    /// The source's `std::hash<DBoundingBox<D>>` value.
    ///
    /// Starting from 0, `hash_combine` folds in `hash_float` of every minimum
    /// coordinate and then of every maximum coordinate
    /// ([`crate::concept::hash_utils`]), so `-0.0` and `+0.0` hash alike. It
    /// is the value the source computes on a 64-bit little-endian host.
    pub fn source_hash64(&self) -> u64 {
        let mut seed = 0;
        for coordinate in self.interval.min_position().iter() {
            hash_combine(&mut seed, hash_float(*coordinate));
        }
        for coordinate in self.interval.max_position().iter() {
            hash_combine(&mut seed, hash_float(*coordinate));
        }
        seed
    }
}

impl DBoundingBox<2> {
    /// Enlarges the box to contain `(x, y)`: source
    /// `enlarge(CoordinateType x, CoordinateType y)`.
    pub fn enlarge_xy(&mut self, x: f64, y: f64) {
        self.enlarge(DPosition::xy(x, y));
    }

    /// Whether `(x, y)` lies inside the box: source
    /// `encloses(CoordinateType x, CoordinateType y)`, "for convenience only".
    pub fn encloses_xy(&self, x: f64, y: f64) -> bool {
        self.encloses(&DPosition::xy(x, y))
    }
}

impl<const D: usize> Deref for DBoundingBox<D> {
    type Target = DIntervalBase<D>;

    fn deref(&self) -> &DIntervalBase<D> {
        &self.interval
    }
}

impl<const D: usize> DerefMut for DBoundingBox<D> {
    fn deref_mut(&mut self) -> &mut DIntervalBase<D> {
        &mut self.interval
    }
}

impl<const D: usize> From<DIntervalBase<D>> for DBoundingBox<D> {
    /// The box that is `interval`: source `operator=(const Base&)`.
    fn from(interval: DIntervalBase<D>) -> Self {
        Self { interval }
    }
}

impl<const D: usize> From<DBoundingBox<D>> for DIntervalBase<D> {
    /// The interval a box is, the source's derived-to-base conversion.
    fn from(bounds: DBoundingBox<D>) -> Self {
        bounds.interval
    }
}

impl<const D: usize> PartialEq<DIntervalBase<D>> for DBoundingBox<D> {
    /// Source `operator==(const Base&)`: equal corners.
    fn eq(&self, other: &DIntervalBase<D>) -> bool {
        self.interval == *other
    }
}

impl<const D: usize> Hash for DBoundingBox<D> {
    /// Feeds [`DBoundingBox::source_hash64`] to the hasher, so equal boxes,
    /// including ones differing only in the sign of a zero, hash alike.
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u64(self.source_hash64());
    }
}

impl<const D: usize> fmt::Display for DBoundingBox<D> {
    /// Source `operator<<`: four lines, each ending in a newline,
    /// `--DBOUNDINGBOX BEGIN--`, `MIN --> <min>`, `MAX --> <max>` and
    /// `--DBOUNDINGBOX END--`.
    ///
    /// A corner is written as the source's `DPosition` output operator writes
    /// it: every coordinate through `precisionWrapper`, separated by single
    /// spaces, so `2` prints as `2.0`.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(formatter, "--DBOUNDINGBOX BEGIN--")?;
        formatter.write_str("MIN --> ")?;
        write_position(formatter, self.interval.min_position())?;
        formatter.write_str("\nMAX --> ")?;
        write_position(formatter, self.interval.max_position())?;
        formatter.write_str("\n--DBOUNDINGBOX END--\n")
    }
}

/// Source `operator<<(std::ostream&, const DPosition&)`: each coordinate
/// through `precisionWrapper`, separated by single spaces.
fn write_position<const D: usize>(
    formatter: &mut fmt::Formatter<'_>,
    position: &DPosition<D>,
) -> fmt::Result {
    for (index, coordinate) in position.iter().enumerate() {
        if index > 0 {
            formatter.write_str(" ")?;
        }
        write!(formatter, "{}", precision_wrapper(*coordinate))?;
    }
    Ok(())
}
