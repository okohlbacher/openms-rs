// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Enums shared across the library: `CONCEPT/CommonEnums.h`.
//!
//! The source collects these small enums in one header "to avoid big includes
//! of large classes and break circular dependencies"; this module keeps that
//! role, and `concept` depends on no other module of the crate.
//!
//! # API mapping
//!
//! | Source (`CommonEnums.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `enum class DIM_UNIT` (`RT`..`IM_CCS`) | [`DimUnit`](crate::concept::common_enums::DimUnit) |
//! | `DIM_UNIT::SIZE_OF_DIM_UNITS` | [`DimUnit::SIZE_OF_DIM_UNITS`](crate::concept::common_enums::DimUnit::SIZE_OF_DIM_UNITS), the length of [`DimUnit::ALL`](crate::concept::common_enums::DimUnit::ALL) |
//! | `DIM_NAMES[]` | [`DIM_NAMES`](crate::concept::common_enums::DIM_NAMES) and [`DimUnit::name`](crate::concept::common_enums::DimUnit::name) |
//! | `DIM_NAMES_SHORT[]` | [`DIM_NAMES_SHORT`](crate::concept::common_enums::DIM_NAMES_SHORT) and [`DimUnit::short_name`](crate::concept::common_enums::DimUnit::short_name) |
//! | `enum class MZ_UNITS` (`DA`, `PPM`) | [`MzUnits`](crate::concept::common_enums::MzUnits) |
//! | `MZ_UNITS::SIZE_OF_MZ_UNITS` | [`MzUnits::SIZE_OF_MZ_UNITS`](crate::concept::common_enums::MzUnits::SIZE_OF_MZ_UNITS), the length of [`MzUnits::ALL`](crate::concept::common_enums::MzUnits::ALL) |
//! | `MZ_UNIT_NAMES[]` | [`MZ_UNIT_NAMES`](crate::concept::common_enums::MZ_UNIT_NAMES) and [`MzUnits::name`](crate::concept::common_enums::MzUnits::name) |
//! | `enum class HydrophobicityScaleMethod` | `chemistry::hydrophobicity::HydrophobicityScale` |
//!
//! `HydrophobicityScaleMethod` lives with the table it indexes, in
//! `chemistry::hydrophobicity`, because `concept` must not depend on
//! `chemistry`: all seven values are there in source order with the source
//! discriminants 0 to 6, and the source's instruction to append a new scale
//! "with an enum value equal to the current maximum enum value plus one" holds
//! for that enum.
//!
//! The `SIZE_OF_*` enumerators are counts, not units. The source can pass one
//! wherever a `DIM_UNIT` is expected and then indexes past the end of the name
//! arrays; here the counts are associated constants, so every [`DimUnit`](crate::concept::common_enums::DimUnit) and
//! [`MzUnits`](crate::concept::common_enums::MzUnits) value is a real unit and [`DimUnit::name`](crate::concept::common_enums::DimUnit::name) cannot index out of
//! range. `(int)DIM_UNIT::X` is `DimUnit::X as usize` (or [`DimUnit::index`](crate::concept::common_enums::DimUnit::index)).

use std::fmt;
use std::str::FromStr;

use crate::Error;

/// A unit that can be shown on a plotting axis: source `enum class DIM_UNIT`.
///
/// The source says "the order is arbitrary", but the discriminants are the
/// source's (`RT = 0` counting up) because [`DIM_NAMES`] and
/// [`DIM_NAMES_SHORT`] are indexed by them.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum DimUnit {
    /// `RT`: retention time in seconds.
    Rt = 0,
    /// `MZ`: mass-to-charge ratio.
    Mz = 1,
    /// `INT`: intensity.
    Int = 2,
    /// `IM_MS`: ion mobility in milliseconds.
    ImMs = 3,
    /// `IM_VSSC`: ion mobility in volt-seconds per square centimetre, i.e. `1/K_0`.
    ImVssc = 4,
    /// `FAIMS_CV`: FAIMS compensation voltage.
    FaimsCv = 5,
    /// `IM_CCS`: collisional cross section in square angstrom.
    ImCcs = 6,
}

impl DimUnit {
    /// Every unit in source order, the order [`DIM_NAMES`] follows.
    pub const ALL: [Self; 7] = [
        Self::Rt,
        Self::Mz,
        Self::Int,
        Self::ImMs,
        Self::ImVssc,
        Self::FaimsCv,
        Self::ImCcs,
    ];

    /// The number of units: source `DIM_UNIT::SIZE_OF_DIM_UNITS`, 7.
    pub const SIZE_OF_DIM_UNITS: usize = Self::ALL.len();

    /// The source discriminant, `(int)unit`, used to index the name arrays.
    pub const fn index(self) -> usize {
        self as usize
    }

    /// The unit for a source discriminant, or `None` at or beyond
    /// [`Self::SIZE_OF_DIM_UNITS`].
    pub fn from_index(index: usize) -> Option<Self> {
        Self::ALL.get(index).copied()
    }

    /// The axis label, source `DIM_NAMES[(int)unit]`, for example `"RT [s]"`.
    pub const fn name(self) -> &'static str {
        DIM_NAMES[self as usize]
    }

    /// The short axis label, source `DIM_NAMES_SHORT[(int)unit]`, for example
    /// `"RT"`.
    ///
    /// Two units share a short name: [`DimUnit::ImMs`] and [`DimUnit::ImVssc`]
    /// are both `"IM"`, as in the source.
    pub const fn short_name(self) -> &'static str {
        DIM_NAMES_SHORT[self as usize]
    }
}

impl fmt::Display for DimUnit {
    /// Writes [`DimUnit::name`].
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}

/// Axis labels of the [`DimUnit`] values, indexed by their discriminant: source
/// `DIM_NAMES`.
pub const DIM_NAMES: [&str; DimUnit::SIZE_OF_DIM_UNITS] = [
    "RT [s]",
    "m/z [Th]",
    "intensity",
    "IM [milliseconds]",
    "IM [vs / cm2]",
    "FAIMS CV",
    "CCS [angstrom^2]",
];

/// Short axis labels of the [`DimUnit`] values, indexed by their discriminant:
/// source `DIM_NAMES_SHORT`.
pub const DIM_NAMES_SHORT: [&str; DimUnit::SIZE_OF_DIM_UNITS] =
    ["RT", "m/z", "int", "IM", "IM", "FCV", "CCS"];

/// The unit of an m/z tolerance: source `enum class MZ_UNITS`.
///
/// Tolerances themselves are `comparison::Tolerance::{Absolute, Ppm}`; this
/// enum names the unit, and [`MzUnits::name`] gives the text the source builds
/// from [`MZ_UNIT_NAMES`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum MzUnits {
    /// `DA`: Dalton.
    Da = 0,
    /// `PPM`: parts per million.
    Ppm = 1,
}

impl MzUnits {
    /// Both units in source order.
    pub const ALL: [Self; 2] = [Self::Da, Self::Ppm];

    /// The number of units: source `MZ_UNITS::SIZE_OF_MZ_UNITS`, 2.
    pub const SIZE_OF_MZ_UNITS: usize = Self::ALL.len();

    /// The source discriminant, `(int)unit`.
    pub const fn index(self) -> usize {
        self as usize
    }

    /// The unit for a source discriminant, or `None` at or beyond
    /// [`Self::SIZE_OF_MZ_UNITS`].
    pub fn from_index(index: usize) -> Option<Self> {
        Self::ALL.get(index).copied()
    }

    /// The unit name, source `MZ_UNIT_NAMES[(int)unit]`: `"Da"` or `"ppm"`.
    pub const fn name(self) -> &'static str {
        MZ_UNIT_NAMES[self as usize]
    }
}

impl fmt::Display for MzUnits {
    /// Writes [`MzUnits::name`].
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}

impl FromStr for MzUnits {
    type Err = Error;

    /// The unit whose [`MZ_UNIT_NAMES`] entry equals `text` exactly.
    ///
    /// The source has no parser of its own; callers look the text up in
    /// `MZ_UNIT_NAMES`, which is what this does. The comparison is
    /// case-sensitive, as that lookup is.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for any other text.
    fn from_str(text: &str) -> Result<Self, Self::Err> {
        Self::ALL
            .into_iter()
            .find(|unit| unit.name() == text)
            .ok_or_else(|| {
                Error::InvalidValue(format!(
                    "unknown m/z unit '{text}'; expected one of {MZ_UNIT_NAMES:?}"
                ))
            })
    }
}

/// Names of the [`MzUnits`] values, indexed by their discriminant: source
/// `MZ_UNIT_NAMES`.
pub const MZ_UNIT_NAMES: [&str; MzUnits::SIZE_OF_MZ_UNITS] = ["Da", "ppm"];
