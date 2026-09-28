// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Section-by-section ports of `GridFeature_test.cpp` and `QTCluster_test.cpp`
//! at source revision bc9cc12514c768385ce121d6ca4bb710fe1983c4
//! (`src/tests/class_tests/openms/source/`).
//!
//! Evidence tier 3: every expected value is a transcribed class-test literal,
//! or, where marked, follows from the source's documented order. No expected
//! value was derived from this port's output.

use openms::analysis::mapmatching::grid_feature::{GridFeature, SourceOrderedSequence};
use openms::analysis::mapmatching::qt_cluster::{BulkData, Element, QTCluster};
use openms::chemistry::AASequence;
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::kernel::features::BaseFeature;

fn seq(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}

/// `bf` of both class tests: RT 1.1, m/z 2.2, charge 3 and two peptide
/// identifications whose first hits are AAA and CCC.
fn annotated_feature() -> BaseFeature {
    let mut bf = BaseFeature::new(1.1, 2.2, 0.0);
    bf.charge = 3;
    for sequence in ["AAA", "CCC"] {
        let mut identification = PeptideIdentification::default();
        identification.hits.push(PeptideHit {
            sequence: seq(sequence),
            ..PeptideHit::default()
        });
        bf.peptide_identifications.push(identification);
    }
    bf
}

// ---------------------------------------------------------------- GridFeature

#[test]
fn grid_feature_class_test() {
    // GridFeature(const BaseFeature& feature, Size map_index, Size feature_index)
    let bf = BaseFeature::default();
    let _gf = GridFeature::new(&bf, 0, 0);

    // const BaseFeature& getFeature() const
    let bf = BaseFeature {
        rt: 1.1,
        mz: 2.2,
        charge: 3,
        ..BaseFeature::default()
    };
    let gf = GridFeature::new(&bf, 0, 0);
    assert!(*gf.feature() == bf);

    let bf = BaseFeature::default();
    assert_eq!(GridFeature::new(&bf, 123, 0).map_index(), 123);
    assert_eq!(GridFeature::new(&bf, 0, 123).feature_index(), 123);
    assert_eq!(GridFeature::new(&bf, 0, 123).id(), 123);

    // const std::set<AASequence>& getAnnotations() const
    let mut bf = BaseFeature::default();
    assert_eq!(GridFeature::new(&bf, 0, 0).annotations().len(), 0);
    bf.peptide_identifications
        .resize(2, PeptideIdentification::default());
    bf.peptide_identifications[0].hits.push(PeptideHit {
        sequence: seq("AAA"),
        ..PeptideHit::default()
    });
    bf.peptide_identifications[1].hits.push(PeptideHit {
        sequence: seq("CCC"),
        ..PeptideHit::default()
    });
    let gf2 = GridFeature::new(&bf, 0, 0);
    assert_eq!(gf2.annotations().len(), 2);
    assert_eq!(gf2.annotations().iter().next().unwrap().0, seq("AAA"));
    assert_eq!(gf2.annotations().iter().next_back().unwrap().0, seq("CCC"));

    // double getRT() const / double getMZ() const
    let bf = BaseFeature {
        rt: 4.56,
        ..BaseFeature::default()
    };
    assert_eq!(GridFeature::new(&bf, 0, 123).rt(), 4.56);
    let bf = BaseFeature {
        mz: 4.56,
        ..BaseFeature::default()
    };
    assert_eq!(GridFeature::new(&bf, 0, 123).mz(), 4.56);
}

#[test]
fn source_sequence_order_is_length_first() {
    // Source `AASequence::operator<` (AASequence.cpp:328-334 at the pin)
    // compares residue counts before residues, so "CC" sorts before "AAA".
    let short = SourceOrderedSequence(seq("CC"));
    let long = SourceOrderedSequence(seq("AAA"));
    assert!(short < long);
    assert!(SourceOrderedSequence(seq("AAA")) < SourceOrderedSequence(seq("CCC")));
    assert!(SourceOrderedSequence(AASequence::default()) < short);
}

// ---------------------------------------------------------------- QTCluster

#[test]
fn qt_cluster_class_test() {
    let bf = annotated_feature();
    let gf = GridFeature::new(&bf, 123, 456);

    // BulkData(...) and QTCluster(BulkData*, bool)
    let _cluster = QTCluster::new(BulkData::new(&gf, 2, 11.1, 0, 0, 0), false);

    let mut cluster = QTCluster::new(BulkData::new(&gf, 2, 11.1, 7, 9, 1), true);
    assert_eq!(cluster.center_rt(), 1.1);
    assert_eq!(cluster.center_mz(), 2.2);
    assert_eq!(cluster.x_coord(), 7);
    assert_eq!(cluster.y_coord(), 9);
    assert_eq!(cluster.id(), 1);
    assert_eq!(cluster.size(), 1);

    let gf2 = GridFeature::new(&bf, 789, 1011);

    // void add(const GridFeature* const element, double distance)
    cluster.initialize_cluster();
    cluster.add(&gf2, 3.3);
    cluster.finalize_cluster();
    assert_eq!(cluster.size(), 2);

    // bool operator<(QTCluster& cluster)
    let cluster2 = QTCluster::new(BulkData::new(&gf, 2, 11.1, 0, 0, 2), false);
    assert!(cluster2 < cluster);

    // Elements getElements() const: either order.
    let elements = cluster.elements();
    assert_eq!(elements.len(), 2);
    let pair = (elements[0].feature, elements[1].feature);
    assert!(
        (std::ptr::eq(pair.0, &gf2) && std::ptr::eq(pair.1, &gf))
            || (std::ptr::eq(pair.0, &gf) && std::ptr::eq(pair.1, &gf2))
    );

    // Elements getAllNeighbors() const
    {
        let gf3 = GridFeature::new(&bf, 789, 1012);
        let gf4 = GridFeature::new(&bf, 222, 1011);
        let mut cluster2 = QTCluster::new(BulkData::new(&gf, 2, 11.1, 0, 0, 2), false);
        assert_eq!(cluster2.all_neighbors().len(), 0);
        cluster2.initialize_cluster();
        cluster2.add(&gf2, 3.3);
        cluster2.finalize_cluster();
        assert_eq!(cluster2.all_neighbors().len(), 1);
        assert!(std::ptr::eq(cluster2.all_neighbors()[0].feature, &gf2));

        // A better feature from the same map does not add a neighbour.
        cluster2.initialize_cluster();
        cluster2.add(&gf3, 3.0);
        cluster2.finalize_cluster();
        assert_eq!(cluster2.all_neighbors().len(), 1);
        assert!(std::ptr::eq(cluster2.all_neighbors()[0].feature, &gf3));

        // Features from a new map do.
        cluster2.initialize_cluster();
        cluster2.add(&gf4, 3.9);
        cluster2.add(&gf4, 3.2);
        cluster2.add(&gf4, 3.1);
        cluster2.add(&gf4, 3.8);
        cluster2.finalize_cluster();
        let neighbors = cluster2.all_neighbors();
        assert_eq!(neighbors.len(), 2);
        let pair = (neighbors[0].feature, neighbors[1].feature);
        assert!(
            (std::ptr::eq(pair.0, &gf4) && std::ptr::eq(pair.1, &gf3))
                || (std::ptr::eq(pair.0, &gf3) && std::ptr::eq(pair.1, &gf4))
        );
    }

    // bool update(const QTCluster::Elements& removed)
    let mut removed = vec![Element {
        map_index: 789,
        feature: &gf2,
    }];
    assert!(cluster.update(&removed));
    assert_eq!(cluster.size(), 1);
    removed.push(Element {
        map_index: 123,
        feature: &gf,
    });
    // Removing the centre invalidates the cluster.
    assert!(!cluster.update(&removed));
    assert!(cluster.is_invalid());

    // double getQuality()
    assert!(cluster.is_invalid());
    let mut cluster = QTCluster::new(BulkData::new(&gf, 2, 11.1, 7, 9, 3), true);
    cluster.initialize_cluster();
    cluster.add(&gf2, 3.3);
    cluster.finalize_cluster();
    assert_eq!(cluster.quality(), (11.1 - 3.3) / 11.1);
    assert!(!cluster.is_invalid());

    // double getCurrentQuality() const
    let current = cluster.current_quality();
    assert_eq!(current, cluster.quality());

    // const set<AASequence>& getAnnotations()
    assert_eq!(cluster.annotations().len(), 1);
    assert_eq!(cluster.annotations().iter().next().unwrap().0, seq("AAA"));
    let mut cluster2 = QTCluster::new(BulkData::new(&gf, 2, 11.1, 0, 0, 2), false);
    assert!(cluster2.annotations().is_empty());

    // inline bool isInvalid() const / void setInvalid()
    assert!(!cluster.is_invalid());
    cluster.set_invalid();
    assert!(cluster.is_invalid());
}

#[test]
fn qt_cluster_id_filter_and_missing_maps() {
    // Documented behaviour of `add` and `computeQuality_`, not class-test
    // literals: with IDs, a neighbour annotated with a disjoint sequence is
    // ignored, and a missing map counts as `max_distance`.
    let bf = annotated_feature();
    let mut other = BaseFeature::default();
    let mut identification = PeptideIdentification::default();
    identification.hits.push(PeptideHit {
        sequence: seq("DDD"),
        ..PeptideHit::default()
    });
    other.peptide_identifications.push(identification);
    let center = GridFeature::new(&bf, 0, 0);
    let disjoint = GridFeature::new(&other, 1, 0);
    let mut cluster = QTCluster::new(BulkData::new(&center, 3, 10.0, 0, 0, 0), true);
    cluster.initialize_cluster();
    cluster.add(&disjoint, 1.0);
    cluster.finalize_cluster();
    assert_eq!(cluster.size(), 1);

    let plain = BaseFeature::default();
    let unannotated = GridFeature::new(&plain, 1, 0);
    let same_map = GridFeature::new(&plain, 0, 1);
    let mut cluster = QTCluster::new(BulkData::new(&center, 3, 10.0, 0, 0, 0), false);
    cluster.initialize_cluster();
    cluster.add(&unannotated, 4.0);
    // An element of the centre's own map is never a neighbour.
    cluster.add(&same_map, 1.0);
    cluster.finalize_cluster();
    assert_eq!(cluster.size(), 2);
    // ((4 + 10) / 2) distance over two other maps; quality (10 - 7) / 10.
    assert_eq!(cluster.quality(), (10.0 - (4.0 + 10.0) / 2.0) / 10.0);
}
