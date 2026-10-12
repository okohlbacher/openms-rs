// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `FeatureGroupingAlgorithmKD_test.cpp` and
//! `MapAlignmentAlgorithmKD_test.cpp` at source revision
//! bc9cc12514c768385ce121d6ca4bb710fe1983c4
//! (`src/tests/class_tests/openms/source/`).
//!
//! Evidence tier 3 for the transcribed class-test sections; the sections
//! marked "tier 4" derive their expected values by hand from the source, as
//! shown at the assertion. Both `group` overloads and the aligner's
//! `addRTFitData`, `fitLOWESS` and `transform` are `NOT_TESTABLE` upstream
//! ("tested in the tool"); their tier-1 differential against the Release
//! build is `tests/mapmatching_kd_grouping_oracle.rs`, and against the
//! retained TOPP outputs `tests/mapmatching_kd_grouping_topp.rs`. No expected
//! value was derived from this port's output.

use openms::analysis::mapmatching::alignment_kd::MapAlignmentAlgorithmKD;
use openms::analysis::mapmatching::feature_grouping_algorithm::FeatureGroupingAlgorithm;
use openms::analysis::mapmatching::feature_grouping_algorithm_kd::{
    ClusterProxyKD, FeatureGroupingAlgorithmKD,
};
use openms::analysis::mapmatching::kd_tree_feature_maps::KDTreeFeatureMaps;
use openms::identification::{PeptideIdentification, ProteinIdentification};
use openms::kernel::features::{ConsensusMap, Feature, FeatureMap};
use openms::param::Param;

#[test]
fn feature_grouping_algorithm_kd_empty_maps_keep_identifications() {
    // Regression test upstream: two input maps with no features yield an
    // empty mass range; the identifications are still transferred.
    let mut algo = FeatureGroupingAlgorithmKD::default();
    let mut maps = vec![FeatureMap::new(), FeatureMap::new()];
    let prot = ProteinIdentification {
        identifier: "run0".into(),
        ..ProteinIdentification::default()
    };
    maps[0].protein_identifications.push(prot);
    let upep = PeptideIdentification {
        identifier: "run0".into(),
        ..PeptideIdentification::default()
    };
    maps[0].unassigned_peptide_identifications.push(upep);
    let mut out = ConsensusMap::new();
    algo.group(&maps, &mut out).unwrap();
    assert_eq!(out.len(), 0);
    assert_eq!(out.protein_identifications.len(), 1);
    assert_eq!(out.unassigned_peptide_identifications.len(), 1);
}

#[test]
fn feature_grouping_algorithm_kd_needs_two_maps() {
    // Tier 4: the source's `Exception::IllegalArgument`.
    let mut algo = FeatureGroupingAlgorithmKD::default();
    let mut out = ConsensusMap::new();
    assert!(algo.group(&[FeatureMap::new()], &mut out).is_err());
}

#[test]
fn feature_grouping_algorithm_kd_links_two_maps() {
    // Tier 4, derived from the source with warping off: two features 1 s and
    // 0.001 Da apart in two maps form one consensus feature of size 2; a
    // third feature 500 s away stays alone. Sizes sort first, so the pair
    // comes first; the consensus position is the mean.
    let mut algo = FeatureGroupingAlgorithmKD {
        warp_enabled: false,
        mz_ppm: false,
        link_mz_tol: 0.01,
        ..FeatureGroupingAlgorithmKD::default()
    };
    let f = |rt: f64, mz: f64, intensity: f32| {
        let mut feature = Feature::new(rt, mz, intensity);
        feature.charge = 2;
        feature
    };
    let mut a = FeatureMap::from_features(vec![f(100.0, 500.0, 10.0), f(600.0, 700.0, 5.0)]);
    let mut b = FeatureMap::from_features(vec![f(101.0, 500.001, 20.0)]);
    a.features[0].unique_id = 1;
    a.features[1].unique_id = 2;
    b.features[0].unique_id = 3;
    let mut out = ConsensusMap::new();
    algo.group(&[a, b], &mut out).unwrap();
    assert_eq!(out.len(), 2);
    let keys: Vec<_> = out.features[0].handles().iter().map(|h| h.key()).collect();
    assert_eq!(keys, vec![(0, 1), (1, 3)]);
    assert!((out.features[0].rt - 100.5).abs() < 1e-9);
    assert!((f64::from(out.features[0].intensity) - 15.0).abs() < 1e-6);
    let keys: Vec<_> = out.features[1].handles().iter().map(|h| h.key()).collect();
    assert_eq!(keys, vec![(0, 2)]);
}

#[test]
fn cluster_proxy_kd_order() {
    // Tier 4: the source's `operator<`: larger size first, then smaller
    // average distance, then larger centre index; the default is invalid.
    let a = ClusterProxyKD::new(3, 0.5, 1);
    let b = ClusterProxyKD::new(2, 0.1, 9);
    let c = ClusterProxyKD::new(3, 0.2, 0);
    let d = ClusterProxyKD::new(3, 0.2, 7);
    assert!(a < b);
    assert!(c < a);
    assert!(d < c);
    assert_eq!(ClusterProxyKD::new(3, 0.2, 7), d);
    assert_ne!(a, b);
    assert!(!ClusterProxyKD::default().is_valid());
    assert!(a.is_valid());
    assert_eq!((a.size(), a.avg_distance(), a.center_index()), (3, 0.5, 1));
}

#[test]
fn map_alignment_algorithm_kd_construction() {
    // `MapAlignmentAlgorithmKD(42, Param())`; without parameters the
    // aligner has no tolerances (uninitialised in the source), so using it
    // fails here.
    let mut aligner = MapAlignmentAlgorithmKD::new(42, &Param::new()).unwrap();
    assert_eq!(aligner.fit_data().len(), 42);
    let maps: Vec<FeatureMap> = Vec::new();
    let kd = KDTreeFeatureMaps::from_maps(&maps).unwrap();
    assert!(aligner.add_rt_fit_data(&kd).is_err());

    // With the KD defaults the aligner is configured.
    let defaults = FeatureGroupingAlgorithmKD::defaults().unwrap();
    let mut aligner = MapAlignmentAlgorithmKD::new(2, &defaults).unwrap();
    aligner.add_rt_fit_data(&kd).unwrap();
    // Fewer than 50 points per map: the identity, with one warning each.
    let warnings = aligner.fit_lowess().unwrap();
    assert_eq!(warnings.len(), 2);
}
