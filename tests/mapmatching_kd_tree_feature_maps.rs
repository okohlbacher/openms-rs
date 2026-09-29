// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `KDTreeFeatureMaps_test.cpp`, `KDTreeFeatureNode_test.cpp` and
//! `FeatureMapping_test.cpp` at source revision
//! bc9cc12514c768385ce121d6ca4bb710fe1983c4
//! (`src/tests/class_tests/openms/source/`).
//!
//! Evidence tier 3: every expected value is a transcribed class-test literal,
//! except the sections marked "tier 4", whose expected values are derived by
//! hand from the source as shown at the assertion. The source's
//! `optimizeTree`, `getNeighborhood`, `queryRegion` and
//! `applyTransformations` sections are `NOT_TESTABLE`; the query order and
//! the map-1 exclusion of `FeatureMapping` are checked against the Release
//! build in `tests/mapmatching_kd_tree_feature_maps_oracle.rs`. No expected
//! value was derived from this port's output.

use openms::analysis::mapmatching::feature_mapping::{
    FeatureMappingInfo, FeatureToMs2Indices, assign_ms2_index_to_feature,
};
use openms::analysis::mapmatching::kd_tree_feature_maps::{KDTreeFeatureMaps, KDTreeFeatureNode};
use openms::kernel::features::{Feature, FeatureMap};
use openms::kernel::{MSExperiment, MSSpectrum, Precursor};

fn similar(actual: f64, expected: f64) {
    let tolerance = 1e-5 * expected.abs().max(actual.abs()).max(1e-5);
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} is not similar to {expected}"
    );
}

fn feature(rt: f64, mz: f64, intensity: f32, charge: i32) -> Feature {
    let mut f = Feature::new(rt, mz, intensity);
    f.charge = charge;
    f
}

/// The class tests' map: `f1` (RT 1000, m/z 400, intensity 100, charge 2)
/// and `f2` (RT 2000, m/z 500, intensity 1000, charge 3).
fn fmaps() -> Vec<FeatureMap> {
    vec![FeatureMap::from_features(vec![
        feature(1000.0, 400.0, 100.0, 2),
        feature(2000.0, 500.0, 1000.0, 3),
    ])]
}

// ---------------------------------------------------------------------------
// KDTreeFeatureMaps_test.cpp
// ---------------------------------------------------------------------------

#[test]
fn kd_tree_feature_maps_class_test() {
    let empty = KDTreeFeatureMaps::new();
    assert_eq!(empty.num_maps(), 0);

    let maps = fmaps();
    let kd_data_1 = KDTreeFeatureMaps::from_maps(&maps).unwrap();

    // Copy constructor and assignment.
    let copy = kd_data_1.clone();
    assert_eq!(copy.size(), kd_data_1.size());
    assert_eq!(copy.size(), 2);
    assert_eq!(copy.mz(0).unwrap(), kd_data_1.mz(0).unwrap());
    assert_eq!(copy.mz(1).unwrap(), kd_data_1.mz(1).unwrap());

    // addMaps, addFeature.
    let mut kd_data_3 = KDTreeFeatureMaps::new();
    kd_data_3.add_maps(&maps).unwrap();
    assert_eq!(kd_data_3.size(), 2);
    let f3 = feature(500.0, 300.0, 0.0, 0);
    kd_data_3.add_feature(2, &f3.base).unwrap();
    assert_eq!(kd_data_3.size(), 3);

    // feature: the stored reference is the map's feature itself.
    assert!(std::ptr::eq(
        kd_data_1.feature(0).unwrap(),
        &maps[0].features[0].base
    ));
    assert!(std::ptr::eq(
        kd_data_1.feature(1).unwrap(),
        &maps[0].features[1].base
    ));

    similar(kd_data_1.rt(0).unwrap(), 1000.0);
    similar(kd_data_1.mz(0).unwrap(), 400.0);
    similar(f64::from(kd_data_1.intensity(0).unwrap()), 100.0);
    assert_eq!(kd_data_1.charge(0).unwrap(), 2);
    assert_eq!(kd_data_1.map_index(0).unwrap(), 0);
    assert_eq!(kd_data_1.size(), 2);
    assert_eq!(kd_data_3.size(), 3);
    assert_eq!(kd_data_1.tree_size(), 2);
    assert_eq!(kd_data_3.tree_size(), 3);
    assert_eq!(kd_data_1.num_maps(), 1);

    // clear, and the RT cache reset regression.
    kd_data_3.clear();
    assert_eq!(kd_data_3.size(), 0);
    assert_eq!(kd_data_3.tree_size(), 0);
    let g = feature(42.0, 200.0, 0.0, 0);
    kd_data_3.add_feature(0, &g.base).unwrap();
    assert_eq!(kd_data_3.size(), 1);
    similar(kd_data_3.rt(0).unwrap(), 42.0);
}

#[test]
fn kd_tree_feature_maps_queries() {
    // Tier 4, derived from the source: the region is inclusive on both
    // bounds; `ignored_map_index` drops that map's features; the neighborhood
    // of a feature excludes its own map by default; the log fold change
    // filter keeps |log10(I2/I1)| <= limit.
    let maps = vec![
        FeatureMap::from_features(vec![
            feature(100.0, 500.0, 100.0, 1),
            feature(110.0, 500.5, 1000.0, 1),
        ]),
        FeatureMap::from_features(vec![
            feature(105.0, 500.2, 300.0, 1),
            feature(300.0, 500.0, 100.0, 1),
        ]),
    ];
    let kd = KDTreeFeatureMaps::from_maps(&maps).unwrap();
    let mut found = vec![99];
    kd.query_region(100.0, 110.0, 500.0, 500.5, &mut found, None);
    found.sort();
    assert_eq!(found, vec![0, 1, 2]);
    kd.query_region(100.0, 110.0, 500.0, 500.5, &mut found, Some(0));
    assert_eq!(found, vec![2]);
    kd.query_region(100.1, 109.9, 500.0, 500.5, &mut found, None);
    assert_eq!(found, vec![2]);

    let mut neighbors = vec![7];
    kd.get_neighborhood(0, &mut neighbors, 10.0, 0.5, false, false, -1.0)
        .unwrap();
    // Appended; only map 1's feature 2 is within RT 90..110, m/z 499.5..500.5.
    assert_eq!(neighbors, vec![7, 2]);
    let mut neighbors = Vec::new();
    kd.get_neighborhood(0, &mut neighbors, 10.0, 0.5, false, true, -1.0)
        .unwrap();
    neighbors.sort();
    assert_eq!(neighbors, vec![0, 1, 2]);
    // log10(1000/100) = 1 > 0.5, log10(300/100) = 0.477 <= 0.5, log10(1) = 0.
    let mut neighbors = Vec::new();
    kd.get_neighborhood(0, &mut neighbors, 10.0, 0.5, false, true, 0.5)
        .unwrap();
    neighbors.sort();
    assert_eq!(neighbors, vec![0, 2]);

    // applyTransformations evaluates each map's transformation at the
    // feature's original RT; the features themselves are unchanged.
    let mut kd = kd;
    kd.apply_transformations(|map, rt| Ok(rt + 1000.0 * map as f64))
        .unwrap();
    similar(kd.rt(2).unwrap(), 1105.0);
    similar(kd.rt(0).unwrap(), 100.0);
    similar(maps[1].features[0].rt, 105.0);
}

// ---------------------------------------------------------------------------
// KDTreeFeatureNode_test.cpp
// ---------------------------------------------------------------------------

#[test]
fn kd_tree_feature_node_class_test() {
    let maps = fmaps();
    let kd_data = KDTreeFeatureMaps::from_maps(&maps).unwrap();
    let node_0 = KDTreeFeatureNode::new(&kd_data, 0);
    assert_eq!(node_0.index(), 0);
    let node_1 = KDTreeFeatureNode::new(&kd_data, 1);
    let copy = node_1;
    assert_eq!(copy.index(), node_1.index());
    similar(copy.get(0).unwrap(), node_1.get(0).unwrap());
    similar(copy.get(1).unwrap(), node_1.get(1).unwrap());
    assert_eq!(node_1.index(), 1);
    similar(node_1.get(0).unwrap(), 2000.0);
    similar(node_1.get(1).unwrap(), 500.0);
    // Tier 4: the source throws `ElementNotFound` for other indices.
    assert!(node_1.get(2).is_err());
}

// ---------------------------------------------------------------------------
// FeatureMapping_test.cpp
// ---------------------------------------------------------------------------

fn ms2(rt: f64, prec_mz: f64, has_prec: bool) -> MSSpectrum {
    let precursors = if has_prec {
        vec![Precursor {
            mz: prec_mz,
            ..Precursor::default()
        }]
    } else {
        Vec::new()
    };
    MSSpectrum {
        ms_level: 2,
        rt,
        precursors,
        ..MSSpectrum::default()
    }
}

fn experiment(spectra: Vec<MSSpectrum>) -> MSExperiment {
    MSExperiment {
        spectra,
        ..MSExperiment::default()
    }
}

#[test]
fn feature_mapping_default_state() {
    let info = FeatureMappingInfo::default();
    assert!(info.feature_maps.is_empty());
    let res = FeatureToMs2Indices::default();
    assert!(res.assigned_ms2.is_empty());
    assert!(res.unassigned_ms2.is_empty());
}

#[test]
fn feature_mapping_assign_ms2_index_to_feature() {
    // Two features at known (RT, m/z): F0 and F1.
    let maps = vec![FeatureMap::from_features(vec![
        feature(100.0, 500.0, 0.0, 0),
        feature(200.0, 600.0, 0.0, 0),
    ])];
    let info = FeatureMappingInfo::new(&maps).unwrap();
    let (f0, f1) = (0usize, 1usize);

    // Basic assignment, MS1 skipped, precursor-less MS2 dropped.
    let ms1 = MSSpectrum {
        ms_level: 1,
        rt: 150.0,
        ..MSSpectrum::default()
    };
    let exp = experiment(vec![
        ms2(100.0, 500.0, true),
        ms2(200.0, 600.0, true),
        ms2(300.0, 999.0, true),
        ms1,
        ms2(100.0, 0.0, false),
    ]);
    let r = assign_ms2_index_to_feature(&exp, &info, 0.5, 5.0, false).unwrap();
    assert_eq!(r.assigned_ms2.len(), 2);
    assert_eq!(r.assigned_ms2[&f0], vec![0]);
    assert_eq!(r.assigned_ms2[&f1], vec![1]);
    assert_eq!(r.unassigned_ms2, vec![2]);

    // RT outside tolerance -> unassigned.
    let exp = experiment(vec![ms2(120.0, 500.0, true)]);
    let r = assign_ms2_index_to_feature(&exp, &info, 0.5, 5.0, false).unwrap();
    assert!(r.assigned_ms2.is_empty());
    assert_eq!(r.unassigned_ms2, vec![0]);

    // ppm tolerance: in-window vs out-of-window.
    let in_win = experiment(vec![ms2(200.0, 600.05, true)]);
    let r = assign_ms2_index_to_feature(&in_win, &info, 200.0, 5.0, true).unwrap();
    assert_eq!(r.assigned_ms2.len(), 1);
    assert_eq!(r.assigned_ms2[&f1].len(), 1);
    let out_win = experiment(vec![ms2(200.0, 600.5, true)]);
    let r2 = assign_ms2_index_to_feature(&out_win, &info, 200.0, 5.0, true).unwrap();
    assert!(r2.assigned_ms2.is_empty());
    assert_eq!(r2.unassigned_ms2.len(), 1);
}

#[test]
fn feature_mapping_closest_mz_wins() {
    let maps = vec![FeatureMap::from_features(vec![
        feature(100.0, 500.0, 0.0, 0),
        feature(100.0, 500.3, 0.0, 0),
    ])];
    let info = FeatureMappingInfo::new(&maps).unwrap();
    let exp = experiment(vec![ms2(100.0, 500.1, true)]);
    let r = assign_ms2_index_to_feature(&exp, &info, 0.5, 5.0, false).unwrap();
    assert_eq!(r.assigned_ms2.len(), 1);
    assert!(r.assigned_ms2.contains_key(&0));
    assert!(!r.assigned_ms2.contains_key(&1));
    assert!(r.unassigned_ms2.is_empty());
}
