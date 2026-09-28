// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Sets of enum flags stored as bits of a `u64` (source
//! `DATASTRUCTURES/FlagSet.h`).
//!
//! The source class template stores combinations of enum values: each value
//! `r` is the bit `2^r`, and sets combine with bitwise operators that accept
//! either another set or a single enum value, never a bare number. This port
//! keeps that shape: [`FlagEnum`] maps an enum value to its bit position, and
//! [`FlagSet`] implements the source's operators through the standard
//! `std::ops` traits. The API mapping is in `docs/QC_SUPPORT.md`.

use std::fmt;
use std::marker::PhantomData;
use std::ops::{Add, AddAssign, BitAnd, BitAndAssign, BitOr, BitOrAssign, Sub, SubAssign};

/// An enum whose values can be members of a [`FlagSet`].
///
/// Source `FlagSet` converts an enum value `r` to the bit `pow(2, r)` and
/// states the precondition `0 <= r <= 63`, checked only in debug builds; a
/// larger value is undefined behaviour in a release build. Here a bit
/// position of 64 or more contributes no bit at all, so it can never alias
/// another flag.
pub trait FlagEnum: Copy {
    /// The bit position of this value, `0..=63`.
    fn bit(self) -> u32;
}

/// A set of [`FlagEnum`] values, one bit per value.
///
/// Source `FlagSet<ENUM>`. The source's copy constructor and assignment are
/// [`Clone`]/[`Copy`]; its deleted assignment from a bare enum value has no
/// Rust counterpart because assignment never converts implicitly. `|` and `+`
/// both add flags, `&` intersects, and `-` removes flags; each accepts either
/// another set or one enum value, as the source's overloads do.
pub struct FlagSet<E> {
    value: u64,
    marker: PhantomData<fn() -> E>,
}

impl<E> Clone for FlagSet<E> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<E> Copy for FlagSet<E> {}
impl<E> Default for FlagSet<E> {
    /// The empty set, source `FlagSet()`.
    fn default() -> Self {
        Self::new()
    }
}
/// Equality compares the stored bits, source `operator==`.
impl<E> PartialEq for FlagSet<E> {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}
impl<E> Eq for FlagSet<E> {}
impl<E> std::hash::Hash for FlagSet<E> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.value.hash(state);
    }
}
impl<E> fmt::Debug for FlagSet<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FlagSet({:#b})", self.value)
    }
}

/// The bit of `flag`, source `getPow_`; see [`FlagEnum`] for positions of 64
/// and above.
fn bit_of<E: FlagEnum>(flag: E) -> u64 {
    1_u64.checked_shl(flag.bit()).unwrap_or(0)
}

impl<E> FlagSet<E> {
    /// The empty set, source `FlagSet()`.
    pub const fn new() -> Self {
        Self {
            value: 0,
            marker: PhantomData,
        }
    }
    /// Whether no bit is set, source `empty`.
    pub fn is_empty(&self) -> bool {
        self.value == 0
    }
    /// The stored bits, source `value` ("mostly for illustrative purposes").
    pub fn value(&self) -> u64 {
        self.value
    }
    /// Whether this set holds every flag of `required`, source
    /// `isSuperSetOf(const FlagSet&)`. Every set is a superset of the empty
    /// set.
    pub fn is_super_set_of(&self, required: &Self) -> bool {
        (self.value | required.value) == self.value
    }
}

impl<E: FlagEnum> FlagSet<E> {
    /// The set holding `flag` alone, source `explicit FlagSet(const ENUM&)`.
    pub fn from_flag(flag: E) -> Self {
        Self {
            value: bit_of(flag),
            marker: PhantomData,
        }
    }
    /// Whether this set holds `required`, source `isSuperSetOf(const ENUM&)`.
    pub fn is_super_set_of_flag(&self, required: E) -> bool {
        let bit = bit_of(required);
        (self.value | bit) == self.value
    }
}

impl<E: FlagEnum> From<E> for FlagSet<E> {
    /// Source `explicit FlagSet(const ENUM&)`; the conversion is explicit in
    /// Rust too, through `from` or `into`.
    fn from(flag: E) -> Self {
        Self::from_flag(flag)
    }
}

impl<E: FlagEnum> BitOr<E> for FlagSet<E> {
    type Output = Self;
    /// Source `operator|(const ENUM&)`.
    fn bitor(mut self, flag: E) -> Self {
        self |= flag;
        self
    }
}
impl<E> BitOr for FlagSet<E> {
    type Output = Self;
    /// Source `operator|(const FlagSet&)`.
    fn bitor(mut self, rhs: Self) -> Self {
        self |= rhs;
        self
    }
}
impl<E: FlagEnum> BitOrAssign<E> for FlagSet<E> {
    /// Source `operator|=(const ENUM&)`.
    fn bitor_assign(&mut self, flag: E) {
        self.value |= bit_of(flag);
    }
}
impl<E> BitOrAssign for FlagSet<E> {
    /// Source `operator|=(const FlagSet&)`.
    fn bitor_assign(&mut self, rhs: Self) {
        self.value |= rhs.value;
    }
}
impl<E: FlagEnum> BitAnd<E> for FlagSet<E> {
    type Output = Self;
    /// Source `operator&(const ENUM&)`.
    fn bitand(mut self, flag: E) -> Self {
        self &= flag;
        self
    }
}
impl<E> BitAnd for FlagSet<E> {
    type Output = Self;
    /// Source `operator&(const FlagSet&)`.
    fn bitand(mut self, rhs: Self) -> Self {
        self &= rhs;
        self
    }
}
impl<E: FlagEnum> BitAndAssign<E> for FlagSet<E> {
    /// Source `operator&=(const ENUM&)`.
    fn bitand_assign(&mut self, flag: E) {
        self.value &= bit_of(flag);
    }
}
impl<E> BitAndAssign for FlagSet<E> {
    /// Source `operator&=(const FlagSet&)`.
    fn bitand_assign(&mut self, rhs: Self) {
        self.value &= rhs.value;
    }
}
impl<E: FlagEnum> Add<E> for FlagSet<E> {
    type Output = Self;
    /// Source `operator+(const ENUM&)`, the same as `|`.
    // The source defines `+` as `|`.
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn add(self, flag: E) -> Self {
        self | flag
    }
}
impl<E> Add for FlagSet<E> {
    type Output = Self;
    /// Source `operator+(const FlagSet&)`, the same as `|`.
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn add(self, rhs: Self) -> Self {
        self | rhs
    }
}
impl<E: FlagEnum> AddAssign<E> for FlagSet<E> {
    /// Source `operator+=(const ENUM&)`, the same as `|=`.
    // The source defines `+=` as `|=`.
    #[allow(clippy::suspicious_op_assign_impl)]
    fn add_assign(&mut self, flag: E) {
        *self |= flag;
    }
}
impl<E> AddAssign for FlagSet<E> {
    /// Source `operator+=(const FlagSet&)`, the same as `|=`.
    #[allow(clippy::suspicious_op_assign_impl)]
    fn add_assign(&mut self, rhs: Self) {
        *self |= rhs;
    }
}
impl<E: FlagEnum> Sub<E> for FlagSet<E> {
    type Output = Self;
    /// Source `operator-(const ENUM&)`: this set without `flag`.
    fn sub(mut self, flag: E) -> Self {
        self -= flag;
        self
    }
}
impl<E> Sub for FlagSet<E> {
    type Output = Self;
    /// Source `operator-(const FlagSet&)`: this set without the flags of
    /// `rhs`.
    fn sub(mut self, rhs: Self) -> Self {
        self -= rhs;
        self
    }
}
impl<E: FlagEnum> SubAssign<E> for FlagSet<E> {
    /// Source `operator-=(const ENUM&)`, which clears the overlapping bit with
    /// an exclusive or.
    fn sub_assign(&mut self, flag: E) {
        let overlap = self.value & bit_of(flag);
        self.value ^= overlap;
    }
}
impl<E> SubAssign for FlagSet<E> {
    /// Source `operator-=(const FlagSet&)`, which clears the overlapping bits
    /// with an exclusive or.
    fn sub_assign(&mut self, rhs: Self) {
        let overlap = self.value & rhs.value;
        self.value ^= overlap;
    }
}
