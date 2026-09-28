// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The remainder of `CHEMISTRY/EnzymaticDigestion.h`: the full `Specificity`
//! enum with its numeric values and name table, the two enzyme-name constants
//! and `filterByMissedCleavages`.
//!
//! The digestion itself is [`ProteaseDigestion`] in
//! [`crate::chemistry::digestion`], whose `specificity` field is a
//! [`DigestionSpecificity`]: the three values the source's digestion supports.
//!
//! # API mapping
//!
//! | Source (`EnzymaticDigestion.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `enum Specificity` (`SPEC_NONE` 0 ... `SPEC_NONTERM` 9) | [`Specificity`], [`Specificity::value`] |
//! | `SIZE_OF_SPECIFICITY` | [`SIZE_OF_SPECIFICITY`] |
//! | `static const std::string NamesOfSpecificity[]` | [`NAMES_OF_SPECIFICITY`] |
//! | `static const std::string NoCleavage`, `UnspecificCleavage` | [`NO_CLEAVAGE`], [`UNSPECIFIC_CLEAVAGE`] |
//! | `static getSpecificityByName(const std::string&)` | [`Specificity::by_name`]; the supported three also `DigestionSpecificity::from_str` |
//! | `EnzymaticDigestion()`, copy constructor, `operator=` | `ProteaseDigestion::default()` (trypsin, 0 missed cleavages, full), `Clone` |
//! | `getMissedCleavages`, `setMissedCleavages` | the `missed_cleavages` field |
//! | `getEnzymeName`, `setEnzyme(const DigestionEnzyme*)` | `enzyme.name()`, [`ProteaseDigestion::set_enzyme`] or the `enzyme` field |
//! | `getSpecificity`, `setSpecificity` | the `specificity` field; [`ProteaseDigestion::source_specificity`], [`ProteaseDigestion::set_source_specificity`] |
//! | both `digestUnmodified` overloads | [`ProteaseDigestion::digest_unmodified`], [`ProteaseDigestion::digest_ranges`] |
//! | `isValidProduct` | `ProteaseDigestion::is_valid_product_unmodified` |
//! | `countInternalCleavageSites` | [`ProteaseDigestion::count_internal_cleavage_sites`] |
//! | `filterByMissedCleavages` | [`ProteaseDigestion::filter_by_missed_cleavages`] |
//! | protected `isValidProduct_`, `tokenize_`, `semiSpecificDigestion_`, `digestAfterTokenize_`, `countMissedCleavages_` and the four data members | the private digestion engine of [`crate::chemistry::digestion`] |
//!
//! # Differences from the source
//!
//! - The source stores any `Specificity` and its digestion treats every value
//!   other than `SPEC_NONE` and `SPEC_SEMI` as fully specific; it documents
//!   `SPEC_NOCTERM` and `SPEC_NONTERM` as "currently not supported in the
//!   class", and `SPEC_UNKNOWN` is what a failed name lookup returns.
//!   [`ProteaseDigestion::set_source_specificity`] refuses those three values
//!   instead of digesting as if they were `SPEC_FULL`.
//! - `digestUnmodified` returns the number of length-discarded products; the
//!   Rust functions return the products (see [`crate::chemistry::digestion`]).

use super::digestion::{DigestionSpecificity, ProteaseDigestion};
use crate::{Error, Result};
use std::fmt;

/// Source `EnzymaticDigestion::SIZE_OF_SPECIFICITY`: the length of
/// [`NAMES_OF_SPECIFICITY`].
pub const SIZE_OF_SPECIFICITY: usize = 10;

/// Source `EnzymaticDigestion::NamesOfSpecificity`, indexed by
/// [`Specificity::value`]; the unused values 3 to 7 are all `unknown`.
pub const NAMES_OF_SPECIFICITY: [&str; SIZE_OF_SPECIFICITY] = [
    "none", "semi", "full", "unknown", "unknown", "unknown", "unknown", "unknown", "no-cterm",
    "no-nterm",
];

/// Source `EnzymaticDigestion::NoCleavage`: the name of the enzyme that does
/// not cleave.
pub const NO_CLEAVAGE: &str = "no cleavage";

/// Source `EnzymaticDigestion::UnspecificCleavage`: the name of the enzyme
/// that cleaves everywhere.
pub const UNSPECIFIC_CLEAVAGE: &str = "unspecific cleavage";

/// Source `EnzymaticDigestion::Specificity`, with its numeric values.
///
/// The first three values are the number of termini that must follow the
/// enzyme's rule, which is how some search engines report specificity.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Specificity {
    /// `SPEC_NONE` (0): no requirement on either terminus.
    None,
    /// `SPEC_SEMI` (1): one of the two termini must follow the rule.
    Semi,
    /// `SPEC_FULL` (2): both termini follow the rule or are protein termini.
    #[default]
    Full,
    /// `SPEC_UNKNOWN` (3): the result of an unknown name.
    Unknown,
    /// `SPEC_NOCTERM` (8): no requirement on the C-terminus; not supported by
    /// the source's digestion.
    NoCTerm,
    /// `SPEC_NONTERM` (9): no requirement on the N-terminus; not supported by
    /// the source's digestion.
    NoNTerm,
}

impl Specificity {
    /// Every value, in source numeric order.
    pub const ALL: [Self; 6] = [
        Self::None,
        Self::Semi,
        Self::Full,
        Self::Unknown,
        Self::NoCTerm,
        Self::NoNTerm,
    ];

    /// The source's numeric value.
    pub fn value(self) -> usize {
        match self {
            Self::None => 0,
            Self::Semi => 1,
            Self::Full => 2,
            Self::Unknown => 3,
            Self::NoCTerm => 8,
            Self::NoNTerm => 9,
        }
    }

    /// The source name, `NamesOfSpecificity[value]`.
    pub fn name(self) -> &'static str {
        NAMES_OF_SPECIFICITY[self.value()]
    }

    /// Source `getSpecificityByName`: the first entry of
    /// [`NAMES_OF_SPECIFICITY`] equal to `name`, and [`Self::Unknown`] when
    /// none is. The comparison is case-sensitive.
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::chemistry::enzymatic_digestion::Specificity;
    ///
    /// assert_eq!(Specificity::by_name("semi"), Specificity::Semi);
    /// assert_eq!(Specificity::by_name("no-cterm"), Specificity::NoCTerm);
    /// assert_eq!(Specificity::by_name("Full"), Specificity::Unknown);
    /// ```
    pub fn by_name(name: &str) -> Self {
        Self::ALL
            .into_iter()
            .find(|spec| spec.name() == name)
            .unwrap_or(Self::Unknown)
    }
}

impl fmt::Display for Specificity {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}

impl From<DigestionSpecificity> for Specificity {
    fn from(spec: DigestionSpecificity) -> Self {
        match spec {
            DigestionSpecificity::None => Self::None,
            DigestionSpecificity::Semi => Self::Semi,
            DigestionSpecificity::Full => Self::Full,
        }
    }
}

impl TryFrom<Specificity> for DigestionSpecificity {
    type Error = Error;

    /// The digestion's specificity for a source value.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Unsupported`] for [`Specificity::Unknown`],
    /// [`Specificity::NoCTerm`] and [`Specificity::NoNTerm`].
    fn try_from(spec: Specificity) -> Result<Self> {
        match spec {
            Specificity::None => Ok(Self::None),
            Specificity::Semi => Ok(Self::Semi),
            Specificity::Full => Ok(Self::Full),
            other => Err(Error::Unsupported(format!(
                "digestion specificity {:?} ({}) is not supported",
                other.name(),
                other.value()
            ))),
        }
    }
}

impl ProteaseDigestion {
    /// The specificity as the source enum: source `getSpecificity()`.
    pub fn source_specificity(&self) -> Specificity {
        self.specificity.into()
    }

    /// Sets the specificity from the source enum: source `setSpecificity`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Unsupported`] for [`Specificity::Unknown`],
    /// [`Specificity::NoCTerm`] and [`Specificity::NoNTerm`] (see the module
    /// notes); the specificity is then unchanged.
    pub fn set_source_specificity(&mut self, spec: Specificity) -> Result<()> {
        self.specificity = DigestionSpecificity::try_from(spec)?;
        Ok(())
    }

    /// Whether `filter` rejects `sequence` by its number of internal cleavage
    /// sites: source `filterByMissedCleavages(sequence, filter)`, which passes
    /// [`Self::count_internal_cleavage_sites`] to `filter` and returns its
    /// answer (`true` meaning "filter this sequence out").
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Self::count_internal_cleavage_sites`], and
    /// [`Error::InvalidValue`] when the count exceeds `i32` (the source's
    /// `Int` argument).
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::chemistry::digestion::ProteaseDigestion;
    ///
    /// let trypsin = ProteaseDigestion::default();
    /// // Two internal sites: after the first K and after the first R.
    /// assert!(trypsin.filter_by_missed_cleavages("ACKDERFGK", |mc| mc > 1).unwrap());
    /// assert!(!trypsin.filter_by_missed_cleavages("ACKDEFGK", |mc| mc > 1).unwrap());
    /// ```
    pub fn filter_by_missed_cleavages<F>(&self, sequence: &str, filter: F) -> Result<bool>
    where
        F: FnOnce(i32) -> bool,
    {
        let count = self.count_internal_cleavage_sites(sequence)?;
        let count = i32::try_from(count)
            .map_err(|_| Error::InvalidValue("missed cleavage count exceeds i32".into()))?;
        Ok(filter(count))
    }
}
