// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! An isotopic cluster of peaks: `DATASTRUCTURES/IsotopeCluster.h`.
//!
//! The source struct "stores information about an isotopic cluster (i.e.
//! potential peptide charge variants)": a group of related peaks that likely
//! originate from the same peptide but with different isotopic compositions,
//! the scans they appear in, and a charge estimate when one is available.
//!
//! # API mapping
//!
//! | Source (`IsotopeCluster.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `struct IsotopeCluster` | [`IsotopeCluster`] |
//! | `IsotopeCluster()` | `IsotopeCluster::default()` |
//! | member `peaks` (`ChargedIndexSet`) | [`IsotopeCluster::peaks`] |
//! | member `scans` (`std::vector<Size>`) | [`IsotopeCluster::scans`] |
//! | `typedef std::pair<Size, Size> IndexPair` | [`IndexPair`] |
//! | `typedef std::set<IndexPair> IndexSet` | [`IndexSet`] |
//! | `struct ChargedIndexSet : public IndexSet` | [`ChargedIndexSet`], the base as the field [`ChargedIndexSet::indices`] plus `Deref`/`DerefMut` |
//! | `ChargedIndexSet()` (charge 0) | `ChargedIndexSet::default()` |
//! | member `charge` (`Int`) | [`ChargedIndexSet::charge`] |
//!
//! The source declares the three types nested in `IsotopeCluster`; Rust has no
//! nested type declarations, so they are module items. `FeatureFinderDefs`
//! names the same three types, and
//! `analysis::feature_finder_picked::defs` re-exports these.
//!
//! The source struct declares no comparison, copy or output operators of its
//! own. `IsotopeCluster` derives `Clone` and `Debug`, which C++ supplies
//! implicitly (the copy constructor) or not at all (`Debug`); it has no
//! equality, as the source has none.
//!
//! See `docs/ISOTOPE_CLUSTER_SUPPORT.md`.

use std::collections::BTreeSet;
use std::ops::{Deref, DerefMut};

/// An index pair, typically `(scan_index, peak_index)` in an experiment:
/// source `IsotopeCluster::IndexPair`, `std::pair<Size, Size>`.
///
/// The first value usually refers to the scan (spectrum) index, the second to
/// the peak index within that scan.
pub type IndexPair = (usize, usize);

/// A set of index pairs, usually referring to peaks in an experiment: source
/// `IsotopeCluster::IndexSet`, `std::set<IndexPair>`.
///
/// It stores unique pairs pointing to specific peaks in specific scans, ordered
/// as the source's `std::set` orders them: by scan index, then by peak index.
pub type IndexSet = BTreeSet<IndexPair>;

/// An index set with an associated charge estimate: source
/// `IsotopeCluster::ChargedIndexSet`.
///
/// It extends the basic [`IndexSet`] with the charge state of the peaks, so
/// that peaks belonging to the same isotopic pattern and the charge they
/// represent are tracked together. The source derives it from [`IndexSet`];
/// here the set is the field [`Self::indices`], and `Deref` gives the set's
/// methods to the struct as the base class does.
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

/// An isotopic cluster: its peaks with their charge estimate and the scans it
/// appears in. Source `struct IsotopeCluster`.
///
/// The default value, source `IsotopeCluster()`, is an empty cluster with no
/// peaks, no scans and charge 0 ("no charge estimate").
///
/// # Examples
///
/// ```
/// use openms::data_structures::isotope_cluster::IsotopeCluster;
///
/// let mut cluster = IsotopeCluster::default();
/// cluster.peaks.insert((0, 3));
/// cluster.peaks.insert((1, 2));
/// cluster.peaks.charge = 2;
/// cluster.scans.extend([0, 1]);
/// assert_eq!(cluster.peaks.len(), 2);
/// assert_eq!(cluster.peaks.iter().next(), Some(&(0, 3)));
/// ```
#[derive(Clone, Debug, Default)]
pub struct IsotopeCluster {
    /// Peaks in this cluster, with their charge state information.
    pub peaks: ChargedIndexSet,
    /// The scan indices where this cluster appears.
    pub scans: Vec<usize>,
}
