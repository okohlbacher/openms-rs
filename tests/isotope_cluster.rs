// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `DATASTRUCTURES/IsotopeCluster.h`.
//!
//! There is no `IsotopeCluster_test.cpp` at core `bc9cc12`. The expected values
//! are the header's documented defaults (tier 3: an empty cluster, charge 0
//! meaning "no charge estimate") and the `std::set` / `std::pair` ordering the
//! typedefs name (tier 4). The `ChargedIndexSet` comparison and ordering rows
//! of the executed `FeatureFinderDefs` probe are in
//! `tests/feature_finder_picked_seeds.rs`.

use openms::analysis::feature_finder_picked::defs;
use openms::data_structures::IsotopeCluster as ReExported;
use openms::data_structures::isotope_cluster::{
    ChargedIndexSet, IndexPair, IndexSet, IsotopeCluster,
};

#[test]
fn default_cluster_is_empty_without_a_charge_estimate() {
    let cluster = IsotopeCluster::default();
    assert!(cluster.peaks.is_empty());
    assert!(cluster.peaks.indices.is_empty());
    assert_eq!(cluster.peaks.charge, 0);
    assert!(cluster.scans.is_empty());
    let _: ReExported = cluster;
}

#[test]
fn peaks_order_by_scan_then_peak_without_duplicates() {
    let mut cluster = IsotopeCluster::default();
    for pair in [(2, 1), (1, 5), (1, 2), (1, 2), (0, 9)] {
        cluster.peaks.insert(pair);
    }
    cluster.peaks.charge = 3;
    cluster.scans = vec![0, 1, 2];
    let order: Vec<IndexPair> = cluster.peaks.iter().copied().collect();
    assert_eq!(order, [(0, 9), (1, 2), (1, 5), (2, 1)]);
    assert_eq!(cluster.scans, [0, 1, 2]);

    let copy = cluster.clone();
    assert_eq!(copy.peaks, cluster.peaks);
    assert_eq!(copy.peaks.charge, 3);
    assert_eq!(copy.scans, cluster.scans);
}

#[test]
fn charged_index_set_compares_indices_only() {
    let indices: IndexSet = [(1, 2), (3, 4)].into_iter().collect();
    let singly = ChargedIndexSet {
        indices: indices.clone(),
        charge: 1,
    };
    let doubly = ChargedIndexSet { indices, charge: 2 };
    assert_eq!(singly, doubly);
    assert_eq!(singly.cmp(&doubly), std::cmp::Ordering::Equal);
    let mut smaller = singly.clone();
    smaller.remove(&(3, 4));
    assert!(smaller < singly);
}

#[test]
fn feature_finder_defs_name_the_same_types() {
    let mut from_defs = defs::ChargedIndexSet::default();
    from_defs.insert((0, 0));
    let cluster = IsotopeCluster {
        peaks: from_defs,
        scans: vec![0],
    };
    let pair: defs::IndexPair = (0, 0);
    let set: &defs::IndexSet = &cluster.peaks;
    assert!(set.contains(&pair));
}
