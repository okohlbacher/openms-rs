// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::PROTON_MASS_U;
use openms::kernel::NumericRange;
use openms::kernel::features::{
    BaseFeature, ColumnHeader, ConsensusFeature, ConsensusMap, Feature, FeatureHandle, FeatureMap,
};
use openms::kernel::geometry::{BoundingBox2D, ConvexHull2D, Point2D};

fn hull(points: &[(f64, f64)]) -> ConvexHull2D {
    ConvexHull2D::from_points(
        &points
            .iter()
            .map(|&(rt, mz)| Point2D::new(rt, mz))
            .collect::<Vec<_>>(),
    )
    .unwrap()
}
fn feature(id: u64, rt: f64, mz: f64, intensity: f32, charge: i32) -> Feature {
    let mut result = Feature::new(rt, mz, intensity);
    result.unique_id = id;
    result.charge = charge;
    result
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-10, "{a} != {b}");
}

#[test]
fn source_hull_containment_and_outline_order() {
    // ConvexHull2D_test.cpp: vec + vec2, including boundary/interpolation cases.
    let h = hull(&[
        (1.0, 2.0),
        (3.0, 4.0),
        (5.0, 0.0),
        (1.0, 1.0),
        (3.0, 1.0),
        (1.0, 3.0),
    ]);
    for &(rt, mz, expected) in &[
        (3.0, 3.0, true),
        (0.0, 0.0, false),
        (6.0, 0.0, false),
        (0.0, 6.0, false),
        (1.5, 1.5, true),
        (1.0, 1.0, true),
        (1.1, 1.0, true),
        (1.2, 2.5, true),
        (1.2, 3.21, false),
        (1.4, 0.99, false),
        (2.5, 1.2, true),
        (1.0, 1.1, true),
        (3.0, 1.0, true),
        (5.0, 0.0, true),
    ] {
        assert_eq!(
            h.encloses(Point2D::new(rt, mz)).unwrap(),
            expected,
            "{rt},{mz}"
        );
    }
    assert_eq!(
        h.hull_points(),
        vec![
            Point2D::new(1.0, 1.0),
            Point2D::new(3.0, 1.0),
            Point2D::new(5.0, 0.0),
            Point2D::new(3.0, 4.0),
            Point2D::new(1.0, 3.0)
        ]
    );
    let bounds = h.bounding_box().unwrap();
    assert_eq!(bounds.min(), Point2D::new(1.0, 0.0));
    assert_eq!(bounds.max(), Point2D::new(5.0, 4.0));
}

#[test]
fn hull_points_duplicates_single_scans_and_empty() {
    let mut h = ConvexHull2D::new();
    assert_eq!(h.bounding_box(), None);
    assert!(!h.encloses(Point2D::new(0.0, 0.0)).unwrap());
    assert_eq!(h.compress(), 0);
    h.expand_to_bounding_box();
    assert!(h.is_empty());
    assert!(h.add_point(Point2D::new(1.0, 2.0)).unwrap());
    assert!(!h.add_point(Point2D::new(1.0, 2.0)).unwrap());
    assert_eq!(h.hull_points(), vec![Point2D::new(1.0, 2.0)]);
    assert!(h.encloses(Point2D::new(1.0, 2.0)).unwrap());
    assert!(!h.encloses(Point2D::new(1.0, 2.001)).unwrap());
    assert!(h.add_point(Point2D::new(1.0, 3.0)).unwrap());
    assert!(!h.add_point(Point2D::new(1.0, 2.5)).unwrap());
    assert_eq!(h.hull_points().len(), 2);
    assert!(h.encloses(Point2D::new(1.0, 2.5)).unwrap());
    assert!(!h.encloses(Point2D::new(1.001, 2.5)).unwrap());
    assert_eq!(h, h.clone());
}

#[test]
fn source_compression_and_expansion() {
    let mut h = hull(&[
        (1.0, 1.0),
        (1.0, 10.0),
        (2.0, 1.0),
        (2.0, 10.0),
        (3.0, 1.0),
        (3.0, 10.0),
    ]);
    let bounds = h.bounding_box();
    let _ = h.hull_points(); // Reading must not create a stale cache after compress.
    assert_eq!(h.compress(), 1);
    assert_eq!(h.compress(), 0);
    assert_eq!(h.hull_points().len(), 4);
    assert_eq!(h.bounding_box(), bounds);
    h.add_points(&[
        Point2D::new(4.0, 1.0),
        Point2D::new(4.0, 10.0),
        Point2D::new(5.0, 2.0),
        Point2D::new(5.0, 10.0),
        Point2D::new(6.0, 1.0),
        Point2D::new(6.0, 10.0),
    ])
    .unwrap();
    assert_eq!(h.compress(), 1);
    assert_eq!(h.compress(), 0);
    for rt in [1.1, 2.1, 3.1, 4.1, 5.1, 5.9] {
        assert!(h.encloses(Point2D::new(rt, 5.0)).unwrap());
    }
    assert!(!h.encloses(Point2D::new(5.1, 1.0)).unwrap());
    let bounds = h.bounding_box();
    h.expand_to_bounding_box();
    assert_eq!(h.hull_points().len(), 4);
    assert_eq!(h.bounding_box(), bounds);
    assert!(h.encloses(Point2D::new(5.1, 1.0)).unwrap());
}

#[test]
fn source_exact_scan_neighbor_fallback_is_preserved() {
    let h = hull(&[(0.0, 0.0), (0.0, 10.0), (1.0, 5.0), (2.0, 0.0), (2.0, 10.0)]);
    // C++ checks strict neighboring scans after an unsuccessful exact-RT test.
    assert!(h.encloses(Point2D::new(1.0, 1.0)).unwrap());
    assert!(!h.encloses(Point2D::new(1.001, 1.0)).unwrap());
}

#[test]
fn outline_mode_and_geometry_validation_are_transactional() {
    let mut h = ConvexHull2D::new();
    h.set_hull_points(&[Point2D::new(1.0, 2.0), Point2D::new(3.0, 4.0)])
        .unwrap();
    let before = h.clone();
    assert!(h.encloses(Point2D::new(2.0, 3.0)).is_err());
    assert!(h.add_point(Point2D::new(2.0, 3.0)).is_err());
    assert_eq!(h, before);
    assert!(h.set_hull_points(&[Point2D::new(f64::NAN, 1.0)]).is_err());
    assert_eq!(h, before);
    h.expand_to_bounding_box();
    assert!(h.encloses(Point2D::new(2.0, 3.0)).unwrap());
    let before = h.clone();
    assert!(
        h.add_points(&[Point2D::new(0.0, 0.0), Point2D::new(5.0, f64::INFINITY)])
            .is_err()
    );
    assert_eq!(h, before);
    assert!(h.encloses(Point2D::new(f64::NAN, 0.0)).is_err());
    let huge = hull(&[(-f64::MAX, 1.0), (f64::MAX, 2.0)]);
    assert!(huge.encloses(Point2D::new(0.0, 1.5)).is_err());
    assert!(BoundingBox2D::new(Point2D::new(2.0, 0.0), Point2D::new(1.0, 0.0)).is_err());
}

#[test]
fn source_feature_hull_union_and_individual_trace_enclosure() {
    let mut f = Feature::default();
    assert!(f.convex_hull().is_empty());
    f.convex_hulls = vec![
        hull(&[(1.0, 2.0), (3.0, 4.0)]),
        hull(&[(0.5, 0.0), (1.0, 1.0)]),
    ];
    assert_eq!(
        f.convex_hull().hull_points(),
        vec![
            Point2D::new(0.5, 0.0),
            Point2D::new(3.0, 0.0),
            Point2D::new(3.0, 4.0),
            Point2D::new(0.5, 4.0)
        ]
    );
    f.convex_hulls[0].add_point(Point2D::new(3.0, 2.0)).unwrap();
    f.convex_hulls[1].add_point(Point2D::new(2.0, 1.0)).unwrap();
    for &(rt, mz, expected) in &[
        (0.0, 0.0, false),
        (1.0, 1.0, true),
        (2.0, 0.5, false),
        (2.0, 3.001, false),
        (2.0, 2.999, true),
        (2.0, 3.5, false),
        (4.0, 3.0, false),
        (1.5, 1.5, false),
        (2.0, 1.0, true),
        (0.5, 0.0, true),
        (3.0, 3.2, true),
    ] {
        assert_eq!(f.encloses(rt, mz).unwrap(), expected);
    }
    assert!(f.convex_hull().encloses(Point2D::new(1.5, 1.5)).unwrap());
    f.convex_hulls.truncate(1);
    assert_eq!(f.convex_hull(), f.convex_hulls[0]);
}

#[test]
fn feature_ranges_include_hulls_but_exclude_subordinates() {
    let mut f = feature(1, 2.0, 3.0, 10.0, 2);
    f.convex_hulls.push(hull(&[(-1.0, 1.0), (10.5, 3.123)]));
    f.subordinates.push(feature(1, 1_000.0, 5_000.0, 99.0, 2));
    let map = FeatureMap::from_features(vec![f, feature(2, 0.0, 0.0, 0.01, 0)]);
    let ranges = map.ranges().unwrap();
    assert_eq!(
        ranges.rt,
        Some(NumericRange {
            min: -1.0,
            max: 10.5
        })
    );
    assert_eq!(
        ranges.mz,
        Some(NumericRange {
            min: 0.0,
            max: 3.123
        })
    );
    assert_eq!(
        ranges.intensity,
        Some(NumericRange {
            min: f64::from(0.01f32),
            max: 10.0
        })
    );
    assert_eq!(FeatureMap::new().ranges().unwrap().rt, None);
}

#[test]
fn feature_maps_sort_stably_and_select_attached_subordinates() {
    let mut a = feature(1, 2.0, 2.0, 1.0, 0);
    a.subordinates.push(feature(99, 1.0, 1.0, 1.0, 0));
    a.metadata.insert("sample".into(), "a".into());
    let b = feature(2, 1.0, 3.0, 1.0, 0);
    let c = feature(3, 1.0, 2.0, 2.0, 0);
    let mut map = FeatureMap::from_features(vec![a, b, c]);
    map.sort_by_rt().unwrap();
    assert_eq!(
        map.features.iter().map(|f| f.unique_id).collect::<Vec<_>>(),
        vec![2, 3, 1]
    );
    map.sort_by_intensity(true).unwrap();
    assert_eq!(
        map.features.iter().map(|f| f.unique_id).collect::<Vec<_>>(),
        vec![3, 2, 1]
    );
    map.sort_by_position().unwrap();
    assert_eq!(
        map.features.iter().map(|f| f.unique_id).collect::<Vec<_>>(),
        vec![3, 2, 1]
    );
    let before = map.clone();
    assert!(map.select(&[0, 0]).is_err());
    assert_eq!(map, before);
    assert!(map.select(&[3]).is_err());
    assert_eq!(map, before);
    map.select(&[2, 0]).unwrap();
    assert_eq!(map.features[0].subordinates[0].unique_id, 99);
    assert_eq!(map.features[0].metadata["sample"].as_str().unwrap(), "a");
    map.select(&[]).unwrap();
    assert!(map.is_empty());
}

#[test]
fn feature_validation_ids_and_depth_are_checked_before_mutation() {
    let mut map = FeatureMap::from_features(vec![
        feature(1, 1.0, 2.0, 1.0, 1),
        feature(1, 0.0, 1.0, 1.0, 1),
    ]);
    assert!(map.validate().is_err());
    assert!(map.sort_by_rt().is_err());
    assert_eq!(map.features[0].rt, 1.0);
    assert!(map.unique_id_to_index(1).is_err());
    map.features[1].unique_id = 0;
    assert_eq!(map.unique_id_to_index(1).unwrap(), Some(0));
    assert_eq!(map.unique_id_to_index(9).unwrap(), None);
    assert!(map.unique_id_to_index(0).is_err());
    map.features[0].width = -1.0;
    assert!(map.select(&[1]).is_err());
    map.features[0].width = 0.0;
    map.features[1].quality_rt = f32::INFINITY;
    assert!(map.sort_by_quality(false).is_err());
    let mut deep = Feature::default();
    for _ in 0..=Feature::MAX_SUBORDINATE_DEPTH {
        deep = Feature {
            subordinates: vec![deep],
            ..Feature::default()
        };
    }
    assert!(deep.validate().is_err());
}

#[test]
fn handles_are_sorted_unique_and_merge_retains_existing_values() {
    let a = feature(5, 1.0, 2.0, 200.0, 2);
    let b = feature(3, 2.0, 3.0, 300.0, 3);
    let mut cf = ConsensusFeature::from_feature(4, &a).unwrap();
    assert_eq!(cf.base, a.base);
    cf.insert(FeatureHandle::new(2, &b)).unwrap();
    assert_eq!(
        cf.handles()
            .iter()
            .map(FeatureHandle::key)
            .collect::<Vec<_>>(),
        vec![(2, 3), (4, 5)]
    );
    let before = cf.clone();
    let mut changed = cf.handles()[0];
    changed.rt = 99.0;
    assert!(cf.insert(changed).is_err());
    assert_eq!(cf, before);
    assert!(cf.set_handles(vec![changed, changed]).is_err());
    assert_eq!(cf, before);
    let mut other = ConsensusFeature::new();
    other.insert(changed).unwrap();
    other.insert(FeatureHandle::new(1, &a)).unwrap();
    cf.merge(&other).unwrap();
    assert_eq!(cf.len(), 3);
    assert_eq!(cf.handles()[1].rt, 2.0);
    assert_eq!(cf.rt, 1.0); // Insertion never implicitly recomputes the summary.
    cf.clear();
    assert!(cf.is_empty());
    assert_eq!(cf.rt, 1.0);
}

#[test]
fn source_consensus_and_monoisotopic_means() {
    let mut cf = ConsensusFeature::new();
    for i in 1..=3 {
        cf.insert(FeatureHandle::new(
            i,
            &BaseFeature::new(i as f64, i as f64 + 1.0, (i + 1) as f32 * 100.0),
        ))
        .unwrap();
        cf.compute_consensus().unwrap();
        close(cf.rt, (i as f64 + 1.0) / 2.0);
        close(cf.mz, (i as f64 + 3.0) / 2.0);
        assert_eq!(cf.intensity, 150.0 + i as f32 * 50.0);
        cf.compute_monoisotopic_consensus().unwrap();
        close(cf.rt, (i as f64 + 1.0) / 2.0);
        close(cf.mz, 2.0);
    }
}

#[test]
fn consensus_charge_running_ties_and_signed_ranges() {
    let mut cf = ConsensusFeature::new();
    for (id, charge) in [2, -2, -2, 2].into_iter().enumerate() {
        cf.insert(FeatureHandle::new(
            id as u64,
            &feature(0, -(id as f64), -10.0, -(id as f32), charge),
        ))
        .unwrap();
    }
    cf.compute_consensus().unwrap();
    assert_eq!(cf.charge, -2); // The running winner, not a final ordered-count tie.
    assert_eq!(cf.intensity, -1.5);
    assert_eq!(
        cf.handle_ranges().rt,
        Some(NumericRange {
            min: -3.0,
            max: 0.0
        })
    );
    cf.insert(FeatureHandle::new(5, &feature(0, 1.0, 1.0, 1.0, i32::MIN)))
        .unwrap();
    cf.compute_monoisotopic_consensus().unwrap();
    assert_eq!(cf.charge, -2);
    cf.set_handles(vec![
        FeatureHandle::new(0, &feature(1, 0.0, 1.0, 1.0, -3)),
        FeatureHandle::new(1, &feature(1, 0.0, 1.0, 1.0, 2)),
    ])
    .unwrap();
    cf.compute_consensus().unwrap();
    assert_eq!(cf.charge, 2);
}

#[test]
fn consensus_empty_and_overflow_leave_summary_unchanged() {
    let mut cf = ConsensusFeature::from(BaseFeature::new(7.0, 8.0, 9.0));
    let before = cf.clone();
    assert!(cf.compute_consensus().is_err());
    assert!(cf.compute_monoisotopic_consensus().is_err());
    assert_eq!(cf, before);
    cf.set_handles(vec![
        FeatureHandle::new(0, &BaseFeature::new(f64::MAX, 1.0, 1.0)),
        FeatureHandle::new(1, &BaseFeature::new(f64::MAX, 1.0, 1.0)),
    ])
    .unwrap();
    let before = cf.clone();
    assert!(cf.compute_consensus().is_err());
    assert_eq!(cf, before);
    assert!(
        cf.insert(FeatureHandle::new(2, &BaseFeature::new(0.0, f64::NAN, 1.0)))
            .is_err()
    );
    assert_eq!(cf, before);
}

#[test]
fn source_decharge_consensus_with_adducts_and_weighting() {
    // ConsensusFeature_test.cpp's three ion states and neutral-mass offsets.
    let sodium = openms::chemistry::element("Na").unwrap().mono_mass();
    let neutral = [1000.5, 1001.0, 999.5];
    let charges = [3, 3, 5];
    let adduct = [
        3.0 * PROTON_MASS_U,
        PROTON_MASS_U + 2.0 * sodium,
        4.0 * PROTON_MASS_U + sodium,
    ];
    let rt = [100.0, 102.0, 101.0];
    let mut cf = ConsensusFeature::new();
    let mut map = FeatureMap::new();
    for i in 0..3 {
        let mut f = feature(
            i as u64 + 1,
            rt[i],
            (neutral[i] + adduct[i]) / f64::from(charges[i]),
            (i + 1) as f32 * 200.0,
            charges[i],
        );
        if i > 0 {
            f.metadata
                .insert("dc_charge_adduct_mass".into(), adduct[i].to_string().into());
        }
        cf.insert(FeatureHandle::new(if i == 0 { 2 } else { 4 }, &f))
            .unwrap();
        map.features.push(f);
    }
    cf.compute_decharge_consensus(&map, true).unwrap();
    close(cf.rt, 100.0 / 6.0 + 102.0 / 3.0 + 101.0 / 2.0);
    close(
        cf.mz,
        neutral[0] / 6.0 + neutral[1] / 3.0 + neutral[2] / 2.0,
    );
    assert_eq!(cf.intensity, 1200.0);
    assert_eq!(cf.charge, 0);
    cf.compute_decharge_consensus(&map, false).unwrap();
    close(cf.rt, 101.0);
    close(cf.mz, neutral.iter().sum::<f64>() / 3.0);
}

#[test]
fn decharge_missing_id_unknown_charge_zero_weight_and_overflow_are_atomic() {
    let mut f = feature(1, 1.0, 100.0, 0.0, 1);
    let mut cf = ConsensusFeature::from_feature(0, &f).unwrap();
    let before = cf.clone();
    assert!(
        cf.compute_decharge_consensus(&FeatureMap::new(), false)
            .is_err()
    );
    assert_eq!(cf, before);
    let mut map = FeatureMap::from_features(vec![f.clone()]);
    assert!(cf.compute_decharge_consensus(&map, true).is_err());
    assert_eq!(cf, before);
    map.features[0]
        .metadata
        .insert("dc_charge_adduct_mass".into(), "NaN".into());
    assert!(cf.compute_decharge_consensus(&map, false).is_err());
    assert_eq!(cf, before);
    f.charge = 0;
    cf.set_handles(vec![FeatureHandle::new(0, &f)]).unwrap();
    assert!(cf.compute_decharge_consensus(&map, false).is_err());
    f.charge = 2;
    f.intensity = f32::MAX;
    map.features[0].metadata.clear();
    map.features.push(feature(2, 2.0, 100.0, f32::MAX, 2));
    cf.set_handles(vec![
        FeatureHandle::new(0, &f),
        FeatureHandle::new(0, &map.features[1]),
    ])
    .unwrap();
    let before = cf.clone();
    assert!(cf.compute_decharge_consensus(&map, false).is_err());
    assert_eq!(cf, before);
}

#[test]
fn negative_charge_decharging_uses_absolute_charge_and_signed_protons() {
    let f = feature(42, -1.0, (100.0 - 2.0 * PROTON_MASS_U) / 2.0, 2.0, -2);
    let map = FeatureMap::from_features(vec![f.clone()]);
    let mut cf = ConsensusFeature::from_feature(99, &f).unwrap();
    cf.compute_decharge_consensus(&map, false).unwrap();
    close(cf.mz, 100.0);
}

#[test]
fn consensus_map_ranges_consistency_and_stable_sorting() {
    let f = feature(1, 1.0, 500.0, 20.0, 1);
    let mut a = ConsensusFeature::from_feature(4, &f).unwrap();
    a.insert(FeatureHandle::new(2, &feature(2, 1000.0, 1500.0, 200.0, 2)))
        .unwrap();
    a.rt = 10.0;
    a.mz = 600.0;
    a.intensity = 100.0;
    let mut b = ConsensusFeature::from_feature(2, &feature(2, -5.0, -7.0, -10.0, 1)).unwrap();
    b.rt = 10.0;
    let mut map = ConsensusMap::from_features(vec![a, b]);
    let ranges = map.ranges().unwrap();
    assert_eq!(
        ranges.rt,
        Some(NumericRange {
            min: -5.0,
            max: 1000.0
        })
    );
    assert_eq!(
        ranges.mz,
        Some(NumericRange {
            min: -7.0,
            max: 1500.0
        })
    );
    assert_eq!(
        ranges.intensity,
        Some(NumericRange {
            min: -10.0,
            max: 200.0
        })
    );
    assert!(map.validate_consistency().is_err());
    map.column_headers.insert(
        2,
        ColumnHeader {
            filename: "a.mzML".into(),
            ..ColumnHeader::default()
        },
    );
    map.column_headers.insert(
        4,
        ColumnHeader {
            filename: "b.mzML".into(),
            ..ColumnHeader::default()
        },
    );
    map.validate_consistency().unwrap(); // size 0 does not limit a unique ID.
    map.sort_by_rt().unwrap();
    assert_eq!(map.features[0].unique_id, 1);
    map.sort_by_size().unwrap();
    assert_eq!(map.features[0].len(), 2);
    map.sort_by_maps().unwrap();
    assert_eq!(map.features[0].len(), 1); // Prefix identity list orders first.
    let before = map.clone();
    assert!(map.select(&[99]).is_err());
    assert_eq!(map, before);
    map.select(&[1]).unwrap();
    assert_eq!(map.features[0].len(), 2);
    assert_eq!(map.column_headers.len(), 2);
    map.column_headers.get_mut(&4).unwrap().filename = "a.mzML".into();
    assert!(map.validate_consistency().is_err());
}

#[test]
fn clear_preserves_metadata_unless_requested() {
    let mut map = FeatureMap::from_features(vec![Feature::default()]);
    map.metadata.insert("a".into(), "b".into());
    map.unique_id = 99;
    map.clear(false);
    assert_eq!(map.unique_id, 99);
    assert_eq!(map.metadata["a"].as_str().unwrap(), "b");
    map.clear(true);
    assert_eq!(map, FeatureMap::new());
    let mut map = ConsensusMap::new();
    assert_eq!(map.experiment_type, "label-free");
    map.experiment_type = "invalid".into();
    assert!(map.validate().is_err());
    map.experiment_type = "labeled_MS2".into();
    map.validate().unwrap();
    map.column_headers.insert(1, ColumnHeader::default());
    map.clear(false);
    assert_eq!(map.column_headers.len(), 1);
    assert_eq!(map.experiment_type, "labeled_MS2");
    map.clear(true);
    assert_eq!(map, ConsensusMap::new());
}

#[test]
fn typed_map_metadata_processing_and_loaded_identity_follow_clear() {
    use openms::format::FileType;
    use openms::metadata::{DataProcessing, MetaValue, ProcessingAction};
    let mut feature = Feature::default();
    feature.set_width(1.25).unwrap();
    assert_eq!(feature.metadata["FWHM"].as_f64().unwrap(), 1.25);
    let before = feature.clone();
    assert!(feature.set_width(f32::NAN).is_err());
    assert_eq!(feature, before);
    let mut map = FeatureMap::from_features(vec![feature]);
    map.metadata
        .insert("samples".into(), MetaValue::from(vec![1_i64, 2]));
    let mut processing = DataProcessing::default();
    processing.actions.insert(ProcessingAction::FeatureGrouping);
    map.data_processing.push(processing.clone());
    map.loaded_file_path = "input.featureXML".into();
    map.loaded_file_type = FileType::FeatureXml;
    map.clear(false);
    assert!(map.is_empty());
    assert_eq!(map.metadata["samples"].as_integer_list().unwrap(), &[1, 2]);
    assert_eq!(map.data_processing, vec![processing]);
    assert_eq!(map.loaded_file_type, FileType::FeatureXml);
    map.clear(true);
    assert_eq!(map, FeatureMap::default());
    let mut consensus = ConsensusMap::default();
    consensus
        .column_headers
        .entry(0)
        .or_default()
        .metadata
        .insert("replicates".into(), 3_i64.into());
    consensus.loaded_file_path = "input.consensusXML".into();
    consensus.loaded_file_type = FileType::ConsensusXml;
    consensus.clear(false);
    assert_eq!(
        consensus.column_headers[&0].metadata["replicates"]
            .as_i64()
            .unwrap(),
        3
    );
    consensus.clear(true);
    assert_eq!(consensus, ConsensusMap::default());
}
