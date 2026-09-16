// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FeatureFinderDefs`: the definitions struct that
//! `FEATUREFINDER/FeatureFinderAlgorithmPicked.h` declares (lines 24-55) "to
//! provide definitions of classes and typedefs which are used throughout all
//! FeatureFinder classes".
//!
//! The source declares the same struct a second time, token for token, in
//! `FEATUREFINDER/FeatureFinderDefs.h`, a header nothing includes. The two
//! definitions cannot meet in one translation unit: including both headers is a
//! redefinition error (executed with the Linux x86_64 Release install's GCC 14,
//! `../oracle/ffap-sem-completion/drivers/defs_both_headers.cpp`). This module
//! is the port of the declaration in `FeatureFinderAlgorithmPicked.h`; the
//! picked feature finder itself uses none of it.
//!
//! The typedefs name `DATASTRUCTURES/IsotopeCluster.h` types, which are not
//! ported on their own; the aliases here spell them out.
//!
//! See `docs/FEATURE_FINDER_PICKED_SUPPORT.md`.

use std::collections::BTreeSet;
use std::fmt;
use std::ops::{Deref, DerefMut};
use std::panic::Location;

use crate::Error;

/// Index to a peak: the scan index and the peak index within that scan.
///
/// Source `FeatureFinderDefs::IndexPair`, a typedef of
/// `IsotopeCluster::IndexPair`, `std::pair<Size, Size>`. The source comment
/// calls the two members `UInt`s; they are `Size`.
pub type IndexPair = (usize, usize);

/// A set of peak indices: source `FeatureFinderDefs::IndexSet`
/// (`IsotopeCluster::IndexSet`, `std::set<IndexPair>`).
///
/// Ordered as the source's `std::set` orders its pairs: by scan index, then by
/// peak index, without duplicates.
pub type IndexSet = BTreeSet<IndexPair>;

/// A set of peak indices with a charge estimate: source
/// `FeatureFinderDefs::ChargedIndexSet` (`IsotopeCluster::ChargedIndexSet`).
///
/// The source derives it from [`IndexSet`]; here the set is the field
/// [`Self::indices`], and `Deref` gives the set's methods to the struct as the
/// base class does.
///
/// The source struct declares no comparison of its own, so `a == b` and `a <
/// b` reach `std::set`'s operators through the base class and compare the
/// index sets only, lexicographically; the charge takes no part. The
/// comparisons here do the same (executed against the Release build:
/// `../oracle/ffap-complete-fix1/drivers/defs_eq_probe.cpp`, where two sets
/// with equal indices and charges 1 and 2 compare equal).
#[derive(Clone, Debug, Default)]
pub struct ChargedIndexSet {
    /// The peak indices (the source's base class).
    pub indices: IndexSet,
    /// Charge estimate; zero, the default, means "no charge estimate".
    pub charge: i32,
}

impl PartialEq for ChargedIndexSet {
    /// `std::operator==` on the base `std::set`: the indices only.
    fn eq(&self, other: &Self) -> bool {
        self.indices == other.indices
    }
}

impl Eq for ChargedIndexSet {}

impl PartialOrd for ChargedIndexSet {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ChargedIndexSet {
    /// `std::operator<` on the base `std::set`: the indices only, compared
    /// lexicographically.
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.indices.cmp(&other.indices)
    }
}

impl Deref for ChargedIndexSet {
    type Target = IndexSet;

    fn deref(&self) -> &IndexSet {
        &self.indices
    }
}

impl DerefMut for ChargedIndexSet {
    fn deref_mut(&mut self) -> &mut IndexSet {
        &mut self.indices
    }
}

/// Whether a peak is already used in a feature: source
/// `FeatureFinderDefs::Flag`.
///
/// The discriminants are the source's enumerator values, and the
/// representation is the 4-byte `int` the Release build gives the enum.
#[repr(i32)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Flag {
    /// `UNUSED`, 0.
    Unused = 0,
    /// `USED`, 1.
    Used = 1,
}

/// The error for an index that has no successor or predecessor: source
/// exception `FeatureFinderDefs::NoSuccessor`, "thrown if a method an invalid
/// IndexPair is given".
///
/// The source constructor takes the throw location (`__FILE__`, `__LINE__`,
/// `OPENMS_PRETTY_FUNCTION`) and the index; [`Self::new`] records the caller's
/// file and line through `#[track_caller]` and has no function name to record.
/// The source constructor also stores its message in the process-wide
/// `GlobalExceptionHandler`; the port keeps no global exception state, so the
/// message is only [`Self::message`]. The source destructor is the default one
/// and has no counterpart. Nothing in the source throws this exception.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NoSuccessor {
    index: IndexPair,
    location: &'static Location<'static>,
}

impl NoSuccessor {
    /// The exception name, source `getName()`.
    pub const NAME: &'static str = "NoSuccessor";

    /// The error for `index`, recording the caller's location.
    #[track_caller]
    pub fn new(index: IndexPair) -> Self {
        Self {
            index,
            location: Location::caller(),
        }
    }

    /// The index without successor or predecessor: the source's protected
    /// member `index_`.
    pub fn index(&self) -> IndexPair {
        self.index
    }

    /// The exception name, [`Self::NAME`].
    pub fn name(&self) -> &'static str {
        Self::NAME
    }

    /// The message, source `what()`: `there is no successor/predecessor for the
    /// given Index: <scan>/<peak>`, both indices in decimal.
    pub fn message(&self) -> String {
        format!(
            "there is no successor/predecessor for the given Index: {}/{}",
            self.index.0, self.index.1
        )
    }

    /// The file of the call to [`Self::new`], source `getFile()`.
    pub fn file(&self) -> &'static str {
        self.location.file()
    }

    /// The line of the call to [`Self::new`], source `getLine()`.
    pub fn line(&self) -> u32 {
        self.location.line()
    }
}

impl fmt::Display for NoSuccessor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message())
    }
}

impl std::error::Error for NoSuccessor {}

impl From<NoSuccessor> for Error {
    /// [`Error::InvalidValue`] with the message `NoSuccessor: <message>`: the
    /// crate has no variant of its own for this exception, and the source
    /// throws it for an invalid argument.
    fn from(error: NoSuccessor) -> Self {
        Error::InvalidValue(format!("{}: {}", NoSuccessor::NAME, error.message()))
    }
}
