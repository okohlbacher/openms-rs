// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The names `CHEMISTRY/EnzymaticDigestion.h` gives its `Specificity`
//! enumeration, its name table and its two enzyme-name constants.
//!
//! The one implementation is [`crate::chemistry::digestion`]:
//! [`DigestionSpecificity`] is the source `EnzymaticDigestion::Specificity`
//! with all six values and their numeric values (`SPEC_NONE` 0 ...
//! `SPEC_NONTERM` 9), `getSpecificityByName` is
//! [`DigestionSpecificity::by_name`], `filterByMissedCleavages` is
//! [`ProteaseDigestion::filter_by_missed_cleavages`] and the digestion
//! arithmetic itself is the `_source` family of [`ProteaseDigestion`]. This
//! module only re-exports them under the source's spelling, so that a caller
//! looking for `EnzymaticDigestion.h` finds them.
//!
//! # API mapping
//!
//! | Source (`EnzymaticDigestion.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `enum Specificity` | [`Specificity`] (`DigestionSpecificity`), [`DigestionSpecificity::value`] |
//! | `SIZE_OF_SPECIFICITY` | [`SIZE_OF_SPECIFICITY`] |
//! | `static const std::string NamesOfSpecificity[]` | [`NAMES_OF_SPECIFICITY`] |
//! | `static const std::string NoCleavage`, `UnspecificCleavage` | [`NO_CLEAVAGE`], [`UNSPECIFIC_CLEAVAGE`] |
//! | `static getSpecificityByName(const std::string&)` | [`Specificity::by_name`] |
//! | `getSpecificity`, `setSpecificity` | the `specificity` field of [`ProteaseDigestion`] |
//! | `filterByMissedCleavages` | [`ProteaseDigestion::filter_by_missed_cleavages`] |
//!
//! See `docs/DIGESTION_SUPPORT.md` for the rest of the class.

pub use super::digestion::{DigestionSpecificity as Specificity, NO_CLEAVAGE, UNSPECIFIC_CLEAVAGE};
#[allow(unused_imports)]
use super::digestion::{DigestionSpecificity, ProteaseDigestion};

/// Source `EnzymaticDigestion::SIZE_OF_SPECIFICITY`: the length of
/// [`NAMES_OF_SPECIFICITY`].
pub const SIZE_OF_SPECIFICITY: usize = 10;

/// Source `EnzymaticDigestion::NamesOfSpecificity`, indexed by
/// [`DigestionSpecificity::value`]; the unused values 3 to 7 are all
/// `unknown`.
pub const NAMES_OF_SPECIFICITY: [&str; SIZE_OF_SPECIFICITY] =
    DigestionSpecificity::NAMES_OF_SPECIFICITY;
