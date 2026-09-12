// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Coordinates in `D`-dimensional space.
//!
//! Native equivalent of `DATASTRUCTURES/DPosition.h` (header-only template;
//! `DPosition.cpp` only instantiates default objects). The source is a template
//! over the dimension and the coordinate type; this port fixes the coordinate
//! type to `f64`, the only instantiation the SDK's public API uses, and keeps
//! the dimension as a const generic. All work is over a fixed-size array, so no
//! bounded-work ceiling applies. See `docs/DPOSITION_SUPPORT.md`.

use std::{
    fmt,
    hash::{Hash, Hasher},
    ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign},
};

/// Representation of a coordinate in `D`-dimensional space.
///
/// Coordinates are `f64` (source `CoordinateType = double`). The array is
/// public: the source exposes each coordinate through a mutable `operator[]`,
/// so direct access hides nothing.
///
/// Equality is ordinary float equality; ordering ([`PartialOrd`]) is
/// lexicographic from dimension 0 to `D - 1`, exactly as the source's
/// `operator<`/`<=`/`>`/`>=` over `std::array`. Any NaN coordinate makes all
/// four ordered comparisons false, as it does in C++. Geometric comparison is
/// [`spatially_less_equal`](Self::spatially_less_equal) and
/// [`spatially_greater_equal`](Self::spatially_greater_equal).
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct DPosition<const D: usize> {
    /// The coordinates, dimension 0 first.
    pub coordinates: [f64; D],
}

/// One-dimensional position.
pub type DPosition1 = DPosition<1>;
/// Two-dimensional position.
pub type DPosition2 = DPosition<2>;

impl<const D: usize> Default for DPosition<D> {
    /// All coordinates zero (source default constructor); a manual impl
    /// because `[f64; D]` only derives `Default` for fixed small `D`.
    fn default() -> Self {
        Self::new()
    }
}

impl<const D: usize> DPosition<D> {
    /// Number of dimensions (source `DIMENSION` enumerator).
    pub const DIMENSION: usize = D;

    /// Position with all coordinates zero (source default constructor).
    pub const fn new() -> Self {
        Self {
            coordinates: [0.; D],
        }
    }

    /// Position with every dimension set to `x` (source `DPosition(CoordinateType x)`).
    pub const fn filled(x: f64) -> Self {
        Self {
            coordinates: [x; D],
        }
    }

    /// All zero (source `zero()`).
    pub const fn zero() -> Self {
        Self::filled(0.)
    }

    /// Smallest positive normal value in every dimension (source
    /// `minPositive()`, `std::numeric_limits<double>::min()`), i.e.
    /// [`f64::MIN_POSITIVE`].
    pub const fn min_positive() -> Self {
        Self::filled(f64::MIN_POSITIVE)
    }

    /// Most negative finite value in every dimension (source `minNegative()`,
    /// `std::numeric_limits<double>::lowest()`), i.e. [`f64::MIN`] = `-f64::MAX`.
    /// This is finite, not negative infinity.
    pub const fn min_negative() -> Self {
        Self::filled(f64::MIN)
    }

    /// Largest finite value in every dimension (source `maxPositive()`,
    /// `std::numeric_limits<double>::max()`), i.e. [`f64::MAX`]. Finite, not
    /// positive infinity.
    pub const fn max_positive() -> Self {
        Self::filled(f64::MAX)
    }

    /// Number of dimensions (source static `size()`).
    pub const fn size() -> usize {
        D
    }

    /// Coordinate `index`, or `None` when `index >= D`.
    ///
    /// The source's `operator[]` only checks the index in debug builds and is
    /// undefined otherwise; [`Index`] on this type panics, and this accessor
    /// never does.
    pub fn get(&self, index: usize) -> Option<f64> {
        self.coordinates.get(index).copied()
    }

    /// Mutable coordinate `index`, or `None` when `index >= D`.
    pub fn get_mut(&mut self, index: usize) -> Option<&mut f64> {
        self.coordinates.get_mut(index)
    }

    /// Makes every coordinate non-negative (source `abs()`).
    ///
    /// The source mutates in place and returns a reference; this consumes the
    /// `Copy` value and returns the result. `f64::abs` clears the sign of
    /// negative zero and keeps NaN as NaN.
    pub fn abs(mut self) -> Self {
        for value in &mut self.coordinates {
            *value = value.abs();
        }
        self
    }

    /// Sets all dimensions to zero (source `clear()`).
    pub fn clear(&mut self) {
        self.coordinates = [0.; D];
    }

    /// Spatially (geometrically) less or equal: every coordinate of `self` is
    /// `<=` the corresponding coordinate of `other` (source
    /// `spatiallyLessEqual`). Implemented as the source does, by returning
    /// `false` at the first coordinate that is `>`; a NaN pair therefore does
    /// not fail the test.
    pub fn spatially_less_equal(&self, other: &Self) -> bool {
        !self
            .coordinates
            .iter()
            .zip(&other.coordinates)
            .any(|(a, b)| a > b)
    }

    /// Spatially (geometrically) greater or equal: every coordinate of `self`
    /// is `>=` the corresponding coordinate of `other` (source
    /// `spatiallyGreaterEqual`). Returns `false` at the first coordinate that
    /// is `<`, as the source does.
    pub fn spatially_greater_equal(&self, other: &Self) -> bool {
        !self
            .coordinates
            .iter()
            .zip(&other.coordinates)
            .any(|(a, b)| a < b)
    }

    /// Inner product, summed from dimension 0 upward (source `operator*(const DPosition&)`).
    pub fn dot(&self, other: &Self) -> f64 {
        self.coordinates
            .iter()
            .zip(&other.coordinates)
            .fold(0., |acc, (a, b)| acc + b * a)
    }

    /// Borrowing iterator over the coordinates (source const `begin()`/`end()`).
    pub fn iter(&self) -> std::slice::Iter<'_, f64> {
        self.coordinates.iter()
    }

    /// Mutable iterator over the coordinates (source mutable `begin()`/`end()`).
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, f64> {
        self.coordinates.iter_mut()
    }

    /// The coordinates as a slice.
    pub fn as_slice(&self) -> &[f64] {
        &self.coordinates
    }
}

impl DPosition<2> {
    /// Two-dimensional position from `x` and `y` (source
    /// `DPosition(CoordinateType x, CoordinateType y)`, whose `static_assert`
    /// restricts it to `D == 2`; here the restriction is the `impl` block).
    pub const fn xy(x: f64, y: f64) -> Self {
        Self {
            coordinates: [x, y],
        }
    }

    /// First dimension (source `getX()`, "for visualization").
    pub const fn x(&self) -> f64 {
        self.coordinates[0]
    }

    /// Second dimension (source `getY()`).
    pub const fn y(&self) -> f64 {
        self.coordinates[1]
    }

    /// Sets the first dimension (source `setX()`).
    pub fn set_x(&mut self, c: f64) {
        self.coordinates[0] = c;
    }

    /// Sets the second dimension (source `setY()`).
    pub fn set_y(&mut self, c: f64) {
        self.coordinates[1] = c;
    }
}

impl DPosition<3> {
    /// Three-dimensional position from `x`, `y` and `z` (source
    /// `DPosition(CoordinateType x, CoordinateType y, CoordinateType z)`).
    pub const fn xyz(x: f64, y: f64, z: f64) -> Self {
        Self {
            coordinates: [x, y, z],
        }
    }
}

impl<const D: usize> From<[f64; D]> for DPosition<D> {
    fn from(coordinates: [f64; D]) -> Self {
        Self { coordinates }
    }
}

impl<const D: usize> From<DPosition<D>> for [f64; D] {
    fn from(position: DPosition<D>) -> Self {
        position.coordinates
    }
}

impl<const D: usize> Index<usize> for DPosition<D> {
    type Output = f64;

    /// Coordinate `index`.
    ///
    /// # Panics
    ///
    /// Panics when `index >= D`. The source's `OPENMS_PRECONDITION` reports
    /// this only in debug builds; use [`DPosition::get`] to avoid the panic.
    fn index(&self, index: usize) -> &f64 {
        &self.coordinates[index]
    }
}

impl<const D: usize> IndexMut<usize> for DPosition<D> {
    /// Mutable coordinate `index`; panics when `index >= D` like [`Index`].
    fn index_mut(&mut self, index: usize) -> &mut f64 {
        &mut self.coordinates[index]
    }
}

impl<'a, const D: usize> IntoIterator for &'a DPosition<D> {
    type Item = &'a f64;
    type IntoIter = std::slice::Iter<'a, f64>;

    fn into_iter(self) -> Self::IntoIter {
        self.coordinates.iter()
    }
}

impl<'a, const D: usize> IntoIterator for &'a mut DPosition<D> {
    type Item = &'a mut f64;
    type IntoIter = std::slice::IterMut<'a, f64>;

    fn into_iter(self) -> Self::IntoIter {
        self.coordinates.iter_mut()
    }
}

impl<const D: usize> IntoIterator for DPosition<D> {
    type Item = f64;
    type IntoIter = std::array::IntoIter<f64, D>;

    fn into_iter(self) -> Self::IntoIter {
        self.coordinates.into_iter()
    }
}

impl<const D: usize> Add for DPosition<D> {
    type Output = Self;

    /// Component-wise addition (source `operator+`).
    fn add(mut self, other: Self) -> Self {
        self += other;
        self
    }
}

impl<const D: usize> AddAssign for DPosition<D> {
    /// Component-wise addition in place (source `operator+=`).
    fn add_assign(&mut self, other: Self) {
        for (a, b) in self.coordinates.iter_mut().zip(other.coordinates) {
            *a += b;
        }
    }
}

impl<const D: usize> Sub for DPosition<D> {
    type Output = Self;

    /// Component-wise subtraction (source `operator-(const DPosition&)`).
    fn sub(mut self, other: Self) -> Self {
        self -= other;
        self
    }
}

impl<const D: usize> SubAssign for DPosition<D> {
    /// Component-wise subtraction in place (source `operator-=`).
    fn sub_assign(&mut self, other: Self) {
        for (a, b) in self.coordinates.iter_mut().zip(other.coordinates) {
            *a -= b;
        }
    }
}

impl<const D: usize> Neg for DPosition<D> {
    type Output = Self;

    /// Component-wise negation (source unary `operator-`). Zero becomes
    /// negative zero, which still compares equal.
    fn neg(mut self) -> Self {
        for value in &mut self.coordinates {
            *value = -*value;
        }
        self
    }
}

impl<const D: usize> Mul for DPosition<D> {
    type Output = f64;

    /// Inner product (source member `operator*(const DPosition&)`); see [`DPosition::dot`].
    fn mul(self, other: Self) -> f64 {
        self.dot(&other)
    }
}

impl<const D: usize> Mul<f64> for DPosition<D> {
    type Output = Self;

    /// Scalar multiplication (source free `operator*(DPosition, scalar)`).
    fn mul(mut self, scalar: f64) -> Self {
        self *= scalar;
        self
    }
}

impl<const D: usize> Mul<DPosition<D>> for f64 {
    type Output = DPosition<D>;

    /// Scalar multiplication (source free `operator*(scalar, DPosition)`).
    fn mul(self, position: DPosition<D>) -> DPosition<D> {
        position * self
    }
}

impl<const D: usize> MulAssign<f64> for DPosition<D> {
    /// Scalar multiplication in place (source `operator*=`).
    fn mul_assign(&mut self, scalar: f64) {
        for value in &mut self.coordinates {
            *value *= scalar;
        }
    }
}

impl<const D: usize> Div<f64> for DPosition<D> {
    type Output = Self;

    /// Scalar division (source free `operator/`). Division by zero follows
    /// IEEE rules, as in the source.
    fn div(mut self, scalar: f64) -> Self {
        self /= scalar;
        self
    }
}

impl<const D: usize> DivAssign<f64> for DPosition<D> {
    /// Scalar division in place (source `operator/=`).
    fn div_assign(&mut self, scalar: f64) {
        for value in &mut self.coordinates {
            *value /= scalar;
        }
    }
}

impl<const D: usize> Hash for DPosition<D> {
    /// Hashes every coordinate in order, normalizing `-0.0` to `+0.0` because
    /// the two compare equal (source `std::hash<DPosition>` via `hash_float`).
    /// No digest compatibility with the source FNV-1a combination is promised.
    fn hash<H: Hasher>(&self, state: &mut H) {
        for value in self.coordinates {
            (if value == 0. { 0 } else { value.to_bits() }).hash(state);
        }
    }
}

impl<const D: usize> fmt::Display for DPosition<D> {
    /// Coordinates separated by single spaces, no newline (source `operator<<`).
    ///
    /// The source writes each coordinate through `precisionWrapper`, i.e.
    /// with 15 significant digits. Without a precision option this uses Rust's
    /// round-trippable float display; `{:.N}` applies fixed-decimal
    /// formatting to every coordinate.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, value) in self.coordinates.iter().enumerate() {
            if i > 0 {
                f.write_str(" ")?;
            }
            match f.precision() {
                Some(p) => write!(f, "{value:.p$}")?,
                None => write!(f, "{value}")?,
            }
        }
        Ok(())
    }
}
