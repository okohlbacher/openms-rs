// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `PROCESSING/FEATURE/FeatureOverlapFilter` and `extern/Quadtree`.
//!
//! - Tier 3: the 14 START_SECTIONs of `FeatureOverlapFilter_test.cpp` at core
//!   `bc9cc12`, transcribed with their literals and checks unchanged.
//! - Tier 1: the C2 FAIMS facts (`../oracle/featurefinder-picked`,
//!   `faims_facts.jsonl`) and the B9 oracle (`../oracle/feature-overlap-filter`):
//!   78 cases executed against the product-sdk libOpenMS, replayed here and
//!   compared bit for bit, including every callback invocation in order. The
//!   four cases the Debug library aborts on are compared with the Release
//!   replica of the pinned source (tier 2), and so are the quadtree probes.
//! - Tier 4: atomicity, refusals and the relation between the public entry
//!   points.
//!
//! See `docs/FEATURE_OVERLAP_FILTER_SUPPORT.md` and
//! `tests/data/feature_overlap_filter_provenance.json`.

use openms::Error;
use openms::concept::constants::user_param::FAIMS_CV;
use openms::concept::unique_id::{HasUniqueId, UniqueIdGenerator};
use openms::identification::ProteinIdentification;
use openms::kernel::{ConvexHull2D, Feature, FeatureMap, Point2D};
use openms::metadata::{MetaValue, MetaValueData};
use openms::processing::feature_overlap_filter::quadtree::{QuadBox, Quadtree, Vector2};
use openms::processing::feature_overlap_filter::{
    CentroidTolerances, FAIMS_MERGE_COUNT, FaimsMergeCallback, FaimsMergeFidelity,
    FeatureOverlapFilter, FeatureOverlapMode, MERGED_CENTROID_IMS, MERGED_CENTROID_MZS,
    MERGED_CENTROID_RTS, MergeIntensityMode,
};
use std::collections::BTreeMap;

const ORACLE: &str = include_str!("data/feature_overlap_filter/feature_overlap_filter_oracle.tsv");
const QUADTREE_GRID: &str = include_str!("data/feature_overlap_filter/quadtree_grid.tsv");
const QUADTREE_NAN: &str = include_str!("data/feature_overlap_filter/quadtree_nan.tsv");

// ---------------------------------------------------------------------------
// Helpers

/// `TEST_REAL_SIMILAR` default: absolute difference within 1e-5, or ratio
/// within 1 + 1e-5.
fn real_similar(a: f64, b: f64) {
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

fn hull(points: &[(f64, f64)]) -> ConvexHull2D {
    let points: Vec<Point2D> = points
        .iter()
        .map(|&(rt, mz)| Point2D::new(rt, mz))
        .collect();
    ConvexHull2D::from_points(&points).unwrap()
}

/// The class test's `createTestFeature` (FeatureOverlapFilter_test.cpp:24-39).
fn create_test_feature(rt: f64, mz: f64, intensity: f64, charge: i32) -> Feature {
    let mut f = Feature::new(rt, mz, intensity as f32);
    f.charge = charge;
    f.convex_hulls = vec![hull(&[
        (rt - 1.0, mz - 0.01),
        (rt + 1.0, mz - 0.01),
        (rt + 1.0, mz + 0.01),
        (rt - 1.0, mz + 0.01),
    ])];
    f
}

fn set_cv(f: &mut Feature, cv: f64) {
    f.metadata
        .insert(FAIMS_CV.into(), MetaValue::try_from(cv).unwrap());
}

fn map_of(features: Vec<Feature>, generator: &mut UniqueIdGenerator) -> FeatureMap {
    let mut map = FeatureMap::from_features(features);
    for f in &mut map.features {
        f.unique_id.ensure_unique_id(generator);
    }
    map
}

fn list_len(f: &Feature, key: &str) -> usize {
    f.metadata[key].as_float_list().unwrap().len()
}

fn merge_overlapping(
    map: &mut FeatureMap,
    same_charge: bool,
    same_im: bool,
    mode: MergeIntensityMode,
    meta: bool,
) {
    FeatureOverlapFilter::merge_overlapping_features(
        map,
        5.0,
        0.05,
        same_charge,
        same_im,
        mode,
        meta,
    )
    .unwrap();
}

// ---------------------------------------------------------------------------
// Tier 3: FeatureOverlapFilter_test.cpp, every START_SECTION

#[test]
fn section_filter_feature_map() {
    let mut generator = UniqueIdGenerator::from_seed(1);
    let with = |quality: f32, rt: f64, mz: f64, points: &[(f64, f64)]| {
        let mut f = Feature::new(rt, mz, 0.5);
        f.quality = quality;
        f.convex_hulls = vec![hull(points)];
        f
    };
    let feature1 = with(8.0, 5.25, 1.5, &[(-1.0, 2.0), (4.0, 1.2), (5.0, 3.123)]);
    let feature2 = with(10.0, 5.25, 1.5, &[(-1.0, 2.0), (4.0, 1.2), (5.5, 3.123)]);
    let feature3 = with(7.0, 5.25, 1.5, &[(4.5, 2.0), (10.0, 1.2), (10.0, 3.123)]);
    let feature4 = with(7.0, 20.0, 10.0, &[(20.0, 5.0), (22.0, 10.0), (22.0, 14.0)]);
    let feature5 = with(0.0, 20.0, 11.0, &[(20.0, 12.0), (21.0, 16.0), (21.0, 18.0)]);
    let mut fmap = map_of(
        vec![feature1, feature2, feature3, feature4, feature5],
        &mut generator,
    );
    FeatureOverlapFilter::filter(
        &mut fmap,
        |left: &Feature, right: &Feature| left.quality > right.quality,
        |_: &mut Feature, _: &mut Feature| true,
        false,
    )
    .unwrap();
    assert_eq!(fmap.features[0].quality, 10.0);
    assert_eq!(fmap.features[1].quality, 7.0);
}

#[test]
fn section_merge_overlapping_features_basic_sum() {
    let mut generator = UniqueIdGenerator::from_seed(2);
    let mut fmap = map_of(
        vec![
            create_test_feature(100.0, 500.0, 1000.0, 2),
            create_test_feature(102.0, 500.02, 500.0, 2),
            create_test_feature(200.0, 600.0, 800.0, 2),
        ],
        &mut generator,
    );
    merge_overlapping(&mut fmap, true, false, MergeIntensityMode::Sum, true);
    assert_eq!(fmap.len(), 2);
    let mut found_merged = false;
    let mut found_separate = false;
    for f in &fmap.features {
        if (f64::from(f.intensity) - 1500.0).abs() < 0.01 {
            found_merged = true;
            assert!(f.metadata.contains_key(MERGED_CENTROID_RTS));
            assert!(f.metadata.contains_key(MERGED_CENTROID_MZS));
            assert_eq!(list_len(f, MERGED_CENTROID_RTS), 2);
        }
        if (f64::from(f.intensity) - 800.0).abs() < 0.01 {
            found_separate = true;
            assert!(!f.metadata.contains_key(MERGED_CENTROID_RTS));
        }
    }
    assert!(found_merged);
    assert!(found_separate);
}

#[test]
fn section_merge_overlapping_features_max_intensity() {
    let mut generator = UniqueIdGenerator::from_seed(3);
    let mut fmap = map_of(
        vec![
            create_test_feature(100.0, 500.0, 1000.0, 2),
            create_test_feature(102.0, 500.02, 500.0, 2),
        ],
        &mut generator,
    );
    merge_overlapping(&mut fmap, true, false, MergeIntensityMode::Max, true);
    assert_eq!(fmap.len(), 1);
    real_similar(f64::from(fmap.features[0].intensity), 1000.0);
}

#[test]
fn section_merge_overlapping_features_require_same_charge() {
    let mut generator = UniqueIdGenerator::from_seed(4);
    let pair = || {
        vec![
            create_test_feature(100.0, 500.0, 1000.0, 2),
            create_test_feature(101.0, 500.01, 500.0, 3),
        ]
    };
    let mut fmap = map_of(pair(), &mut generator);
    merge_overlapping(&mut fmap, true, false, MergeIntensityMode::Sum, true);
    assert_eq!(fmap.len(), 2);

    let mut fmap = map_of(pair(), &mut generator);
    merge_overlapping(&mut fmap, false, false, MergeIntensityMode::Sum, true);
    assert_eq!(fmap.len(), 1);
    real_similar(f64::from(fmap.features[0].intensity), 1500.0);
}

fn cv_pair(first: Option<f64>, second: Option<f64>) -> Vec<Feature> {
    let mut f1 = create_test_feature(100.0, 500.0, 1000.0, 2);
    if let Some(cv) = first {
        set_cv(&mut f1, cv);
    }
    let mut f2 = create_test_feature(101.0, 500.01, 500.0, 2);
    if let Some(cv) = second {
        set_cv(&mut f2, cv);
    }
    vec![f1, f2]
}

#[test]
fn section_merge_overlapping_features_require_same_im_with_faims_cv() {
    let mut generator = UniqueIdGenerator::from_seed(5);
    let mut fmap = map_of(cv_pair(Some(-45.0), Some(-60.0)), &mut generator);
    merge_overlapping(&mut fmap, true, true, MergeIntensityMode::Sum, true);
    assert_eq!(fmap.len(), 2);

    let mut fmap = map_of(cv_pair(Some(-45.0), Some(-60.0)), &mut generator);
    merge_overlapping(&mut fmap, true, false, MergeIntensityMode::Sum, true);
    assert_eq!(fmap.len(), 1);
    real_similar(f64::from(fmap.features[0].intensity), 1500.0);
    assert!(fmap.features[0].metadata.contains_key(MERGED_CENTROID_IMS));
    assert_eq!(list_len(&fmap.features[0], MERGED_CENTROID_IMS), 2);
    assert_eq!(
        fmap.features[0].metadata[FAIMS_MERGE_COUNT]
            .as_i64()
            .unwrap(),
        2
    );
}

#[test]
fn section_merge_overlapping_features_require_same_im_with_same_faims_cv() {
    let mut generator = UniqueIdGenerator::from_seed(6);
    let mut fmap = map_of(cv_pair(Some(-45.0), Some(-45.0)), &mut generator);
    merge_overlapping(&mut fmap, true, true, MergeIntensityMode::Sum, true);
    assert_eq!(fmap.len(), 1);
    real_similar(f64::from(fmap.features[0].intensity), 1500.0);
}

#[test]
fn section_merge_overlapping_features_without_faims_cv() {
    let mut generator = UniqueIdGenerator::from_seed(7);
    let mut fmap = map_of(cv_pair(None, None), &mut generator);
    merge_overlapping(&mut fmap, true, true, MergeIntensityMode::Sum, true);
    assert_eq!(fmap.len(), 1);
    real_similar(f64::from(fmap.features[0].intensity), 1500.0);
    assert!(!fmap.features[0].metadata.contains_key(MERGED_CENTROID_IMS));
    assert!(!fmap.features[0].metadata.contains_key(FAIMS_MERGE_COUNT));
}

#[test]
fn section_merge_overlapping_features_mixed_faims_cv_presence() {
    let mut generator = UniqueIdGenerator::from_seed(8);
    let mut fmap = map_of(cv_pair(Some(-45.0), None), &mut generator);
    merge_overlapping(&mut fmap, true, true, MergeIntensityMode::Sum, true);
    assert_eq!(fmap.len(), 2);
}

#[test]
fn section_merge_overlapping_features_write_meta_values_false() {
    let mut generator = UniqueIdGenerator::from_seed(9);
    let mut fmap = map_of(cv_pair(Some(-45.0), Some(-60.0)), &mut generator);
    merge_overlapping(&mut fmap, true, false, MergeIntensityMode::Sum, false);
    assert_eq!(fmap.len(), 1);
    real_similar(f64::from(fmap.features[0].intensity), 1500.0);
    for key in [
        MERGED_CENTROID_RTS,
        MERGED_CENTROID_MZS,
        MERGED_CENTROID_IMS,
        FAIMS_MERGE_COUNT,
    ] {
        assert!(!fmap.features[0].metadata.contains_key(key));
    }
}

#[test]
fn section_merge_overlapping_features_no_merge_outside_tolerance() {
    let mut generator = UniqueIdGenerator::from_seed(10);
    let mut fmap = map_of(
        vec![
            create_test_feature(100.0, 500.0, 1000.0, 2),
            create_test_feature(110.0, 500.0, 500.0, 2),
        ],
        &mut generator,
    );
    merge_overlapping(&mut fmap, true, false, MergeIntensityMode::Sum, true);
    assert_eq!(fmap.len(), 2);
}

#[test]
fn section_merge_overlapping_features_multiple_features() {
    let mut generator = UniqueIdGenerator::from_seed(11);
    let mut fmap = map_of(
        vec![
            create_test_feature(100.0, 500.0, 1000.0, 2),
            create_test_feature(101.0, 500.01, 500.0, 2),
            create_test_feature(102.0, 500.02, 300.0, 2),
        ],
        &mut generator,
    );
    merge_overlapping(&mut fmap, true, false, MergeIntensityMode::Sum, true);
    assert_eq!(fmap.len(), 1);
    real_similar(f64::from(fmap.features[0].intensity), 1800.0);
    assert_eq!(list_len(&fmap.features[0], MERGED_CENTROID_RTS), 3);
}

#[test]
fn section_merge_faims_features_only_merges_different_faims_cv() {
    let mut generator = UniqueIdGenerator::from_seed(12);
    let mut features = cv_pair(Some(-45.0), Some(-60.0));
    features.push(create_test_feature(100.5, 500.005, 800.0, 2));
    features.push(create_test_feature(101.5, 500.015, 400.0, 2));
    let mut fmap = map_of(features, &mut generator);
    FeatureOverlapFilter::merge_faims_features(&mut fmap, 5.0, 0.05).unwrap();
    assert_eq!(fmap.len(), 3);
    let mut faims_count = 0;
    let mut non_faims_count = 0;
    for f in &fmap.features {
        if f.metadata.contains_key(MERGED_CENTROID_IMS) || f.metadata.contains_key(FAIMS_CV) {
            faims_count += 1;
            if f.metadata.contains_key(MERGED_CENTROID_IMS) {
                real_similar(f64::from(f.intensity), 1500.0);
            }
        } else {
            non_faims_count += 1;
        }
    }
    assert_eq!(faims_count, 1);
    assert_eq!(non_faims_count, 2);
}

#[test]
fn section_merge_faims_features_does_not_merge_same_faims_cv() {
    let mut generator = UniqueIdGenerator::from_seed(13);
    let mut fmap = map_of(cv_pair(Some(-45.0), Some(-45.0)), &mut generator);
    FeatureOverlapFilter::merge_faims_features(&mut fmap, 5.0, 0.05).unwrap();
    assert_eq!(fmap.len(), 2);
    real_similar(f64::from(fmap.features[0].intensity), 1000.0);
    real_similar(f64::from(fmap.features[1].intensity), 500.0);
    assert!(!fmap.features[0].metadata.contains_key(MERGED_CENTROID_IMS));
    assert!(!fmap.features[1].metadata.contains_key(MERGED_CENTROID_IMS));
}

#[test]
fn section_merge_faims_features_no_op_on_non_faims_data() {
    let mut generator = UniqueIdGenerator::from_seed(14);
    let mut fmap = map_of(cv_pair(None, None), &mut generator);
    FeatureOverlapFilter::merge_faims_features(&mut fmap, 5.0, 0.05).unwrap();
    assert_eq!(fmap.len(), 2);
    real_similar(f64::from(fmap.features[0].intensity), 1000.0);
    real_similar(f64::from(fmap.features[1].intensity), 500.0);
}

// ---------------------------------------------------------------------------
// Tier 1: the C2 FAIMS facts (faims_facts.jsonl, the literals of faims_probe.cpp)

fn c2_feature(rt: f64, mz: f64, intensity: f32, cv: f64) -> Feature {
    let mut f = Feature::new(rt, mz, intensity);
    f.charge = 2;
    set_cv(&mut f, cv);
    f
}

fn float_list(f: &Feature, key: &str) -> Vec<u64> {
    f.metadata[key]
        .as_float_list()
        .unwrap()
        .iter()
        .map(|v| v.to_bits())
        .collect()
}

#[test]
fn c2_unique_id_zero_merge_removes_every_faims_feature() {
    let mut fmap = FeatureMap::from_features(vec![
        c2_feature(100.0, 500.0, 1000.0, -45.0),
        c2_feature(101.0, 500.01, 500.0, -60.0),
        c2_feature(900.0, 700.0, 300.0, -45.0),
    ]);
    FeatureOverlapFilter::merge_faims_features(&mut fmap, 5.0, 0.05).unwrap();
    // "before":3,"after":0,"features":[]
    assert!(fmap.is_empty());
}

#[test]
fn c2_valid_ids_pair_plus_far_control() {
    let mut fmap = FeatureMap::from_features(vec![
        c2_feature(100.0, 500.0, 1000.0, -45.0),
        c2_feature(101.0, 500.01, 500.0, -60.0),
        c2_feature(900.0, 700.0, 300.0, -45.0),
    ]);
    for (i, f) in fmap.features.iter_mut().enumerate() {
        f.unique_id = i as u64 + 1;
    }
    FeatureOverlapFilter::merge_faims_features(&mut fmap, 5.0, 0.05).unwrap();
    assert_eq!(fmap.len(), 2);
    let merged = &fmap.features[0];
    // "intensity_f32":"0x1.77p+10"
    assert_eq!(merged.intensity.to_bits(), 1500.0f32.to_bits());
    assert_eq!(
        float_list(merged, MERGED_CENTROID_RTS),
        [100.0f64.to_bits(), 101.0f64.to_bits()]
    );
    assert_eq!(
        float_list(merged, MERGED_CENTROID_MZS),
        [500.0f64.to_bits(), 0x407f_4028_f5c2_8f5c]
    );
    assert_eq!(
        float_list(merged, MERGED_CENTROID_IMS),
        [(-45.0f64).to_bits(), (-60.0f64).to_bits()]
    );
    assert_eq!(merged.metadata[FAIMS_MERGE_COUNT].as_i64().unwrap(), 2);
    assert!(!merged.metadata.contains_key(FAIMS_CV));
    let far = &fmap.features[1];
    assert_eq!(far.intensity, 300.0);
    assert_eq!(far.metadata[FAIMS_CV].as_f64().unwrap(), -45.0);
    assert_eq!(far.metadata.len(), 1);
}

#[test]
fn c2_three_voltages_merge_twice_into_1900_and_1700() {
    let mut fmap = FeatureMap::from_features(vec![
        c2_feature(100.0, 500.0, 1000.0, -45.0),
        c2_feature(100.5, 500.01, 900.0, -60.0),
        c2_feature(101.0, 500.02, 800.0, -75.0),
    ]);
    for (i, f) in fmap.features.iter_mut().enumerate() {
        f.unique_id = i as u64 + 1;
    }
    FeatureOverlapFilter::merge_faims_features(&mut fmap, 5.0, 0.05).unwrap();
    // "intensities_f32":["0x1.dbp+10","0x1.a9p+10"]
    let intensities: Vec<f32> = fmap.features.iter().map(|f| f.intensity).collect();
    assert_eq!(intensities, [1900.0, 1700.0]);
    let first = &fmap.features[0];
    assert_eq!(
        float_list(first, MERGED_CENTROID_RTS),
        [100.0f64.to_bits(), 100.5f64.to_bits()]
    );
    assert_eq!(
        float_list(first, MERGED_CENTROID_IMS),
        [(-45.0f64).to_bits(), (-60.0f64).to_bits()]
    );
    let second = &fmap.features[1];
    assert_eq!(second.rt, 101.0);
    assert_eq!(
        float_list(second, MERGED_CENTROID_RTS),
        [101.0f64.to_bits(), 100.5f64.to_bits()]
    );
    assert_eq!(
        float_list(second, MERGED_CENTROID_MZS),
        [0x407f_4051_eb85_1eb8, 0x407f_4028_f5c2_8f5c]
    );
    assert_eq!(
        float_list(second, MERGED_CENTROID_IMS),
        [(-75.0f64).to_bits(), (-60.0f64).to_bits()]
    );
    for f in &fmap.features {
        assert_eq!(f.metadata[FAIMS_MERGE_COUNT].as_i64().unwrap(), 2);
        assert!(!f.metadata.contains_key(FAIMS_CV));
    }
}

// ---------------------------------------------------------------------------
// Tier 1: the B9 oracle records

fn bits64(text: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(text, 16).unwrap())
}

fn bits32(text: &str) -> f32 {
    f32::from_bits(u32::from_str_radix(text, 16).unwrap())
}

fn unhex(text: &str) -> String {
    let bytes: Vec<u8> = (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).unwrap())
        .collect();
    String::from_utf8(bytes).unwrap()
}

fn hex(text: &str) -> String {
    use std::fmt::Write;
    text.bytes().fold(String::new(), |mut out, b| {
        let _ = write!(out, "{b:02x}");
        out
    })
}

fn fields(text: &str) -> BTreeMap<&str, &str> {
    text.split(';')
        .map(|field| field.split_once('=').unwrap())
        .collect()
}

fn counted(text: &str) -> (usize, &str) {
    let (count, rest) = text.split_once('#').unwrap();
    (count.parse().unwrap(), rest)
}

fn parse_hulls(text: &str) -> Vec<ConvexHull2D> {
    let (count, rest) = counted(text);
    if count == 0 {
        return Vec::new();
    }
    let hulls: Vec<ConvexHull2D> = rest
        .split('/')
        .map(|hull| {
            if hull.is_empty() {
                return ConvexHull2D::new();
            }
            let points: Vec<Point2D> = hull
                .split(',')
                .map(|point| {
                    let (rt, mz) = point.split_once(':').unwrap();
                    Point2D::new(bits64(rt), bits64(mz))
                })
                .collect();
            ConvexHull2D::from_points(&points).unwrap()
        })
        .collect();
    assert_eq!(hulls.len(), count);
    hulls
}

fn parse_meta(text: &str) -> openms::metadata::MetaInfo {
    let (count, rest) = counted(text);
    let mut meta = openms::metadata::MetaInfo::new();
    if count == 0 {
        return meta;
    }
    for entry in rest.split(',') {
        let mut parts = entry.splitn(3, '~');
        let key = parts.next().unwrap();
        let kind = parts.next().unwrap();
        let value = parts.next().unwrap();
        let value = match kind {
            "e" => MetaValue::default(),
            "d" => MetaValue::try_from(bits64(value)).unwrap(),
            "i" => MetaValue::from(value.parse::<i64>().unwrap()),
            "dl" => MetaValue::try_from(value.split(':').map(bits64).collect::<Vec<_>>()).unwrap(),
            "s" => MetaValue::from(unhex(value)),
            other => panic!("unknown meta type {other}"),
        };
        meta.insert(key.to_owned(), value);
    }
    assert_eq!(meta.len(), count);
    meta
}

fn parse_feature(text: &str) -> Feature {
    let f = fields(text);
    let mut feature = Feature::new(bits64(f["rt"]), bits64(f["mz"]), bits32(f["int"]));
    feature.quality = bits32(f["q"]);
    feature.charge = f["z"].parse().unwrap();
    feature.unique_id = f["uid"].parse().unwrap();
    feature.metadata = parse_meta(f["meta"]);
    feature.convex_hulls = parse_hulls(f["hulls"]);
    let (count, rest) = counted(f["subs"]);
    if count > 0 {
        for sub in rest.split('|') {
            feature.subordinates.push(Feature {
                convex_hulls: parse_hulls(sub),
                ..Feature::default()
            });
        }
    }
    assert_eq!(feature.subordinates.len(), count);
    feature
}

fn encode_hulls(hulls: &[ConvexHull2D]) -> String {
    let parts: Vec<String> = hulls
        .iter()
        .map(|hull| {
            hull.hull_points()
                .iter()
                .map(|p| format!("{:016x}:{:016x}", p.rt.to_bits(), p.mz.to_bits()))
                .collect::<Vec<_>>()
                .join(",")
        })
        .collect();
    format!("{}#{}", hulls.len(), parts.join("/"))
}

fn encode_meta(meta: &openms::metadata::MetaInfo) -> String {
    let parts: Vec<String> = meta
        .iter()
        .map(|(key, value)| {
            let encoded = match value.data() {
                MetaValueData::Empty => "e~".to_owned(),
                MetaValueData::Float(v) => format!("d~{:016x}", v.to_bits()),
                MetaValueData::Integer(v) => format!("i~{v}"),
                MetaValueData::FloatList(values) => format!(
                    "dl~{}",
                    values
                        .iter()
                        .map(|v| format!("{:016x}", v.to_bits()))
                        .collect::<Vec<_>>()
                        .join(":")
                ),
                MetaValueData::String(v) => format!("s~{}", hex(v)),
                other => panic!("unexpected meta value {other:?}"),
            };
            format!("{key}~{encoded}")
        })
        .collect();
    format!("{}#{}", meta.len(), parts.join(","))
}

fn encode_feature(f: &Feature) -> String {
    let subs: Vec<String> = f
        .subordinates
        .iter()
        .map(|sub| encode_hulls(&sub.convex_hulls))
        .collect();
    format!(
        "rt={:016x};mz={:016x};int={:08x};q={:08x};z={};uid={};meta={};hulls={};subs={}#{}",
        f.rt.to_bits(),
        f.mz.to_bits(),
        f.intensity.to_bits(),
        f.quality.to_bits(),
        f.charge,
        f.unique_id,
        encode_meta(&f.metadata),
        encode_hulls(&f.convex_hulls),
        f.subordinates.len(),
        subs.join("|")
    )
}

fn encode_map(map: &FeatureMap) -> String {
    format!(
        "id={};uid={};meta={};proteins={}",
        hex(&map.identifier),
        map.unique_id,
        encode_meta(&map.metadata),
        map.protein_identifications.len()
    )
}

fn apply_map_fields(map: &mut FeatureMap, text: &str) {
    let f = fields(text);
    map.identifier = unhex(f["id"]);
    map.unique_id = f["uid"].parse().unwrap();
    map.metadata = parse_meta(f["meta"]);
    let proteins: usize = f["proteins"].parse().unwrap();
    map.protein_identifications = vec![ProteinIdentification::default(); proteins];
}

fn tag_of(f: &Feature) -> i64 {
    f.metadata.get("tag").map_or(-1, |v| v.as_i64().unwrap())
}

#[derive(Default)]
struct OracleCase {
    name: String,
    source: String,
    library_exit: i32,
    inputs: Vec<String>,
    map_before: Option<String>,
    map_after: Option<String>,
    op: String,
    callbacks: Vec<(i64, i64, bool)>,
    exception: Option<String>,
    size: Option<usize>,
    outputs: Vec<String>,
}

fn oracle_cases() -> Vec<OracleCase> {
    let mut cases: Vec<OracleCase> = Vec::new();
    let mut current: Option<OracleCase> = None;
    for line in ORACLE.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        match parts[0] {
            "case" => {
                current = Some(OracleCase {
                    name: parts[1].to_owned(),
                    source: parts[2].trim_start_matches("source=").to_owned(),
                    library_exit: parts[3]
                        .trim_start_matches("library_exit=")
                        .parse()
                        .unwrap(),
                    ..OracleCase::default()
                });
            }
            "end" => cases.push(current.take().unwrap()),
            _ => {
                let case = current.as_mut().unwrap();
                match parts[0] {
                    "in" => case.inputs.push(parts[2].to_owned()),
                    "inref" => {
                        let earlier = cases.iter().find(|c| c.name == parts[1]).unwrap();
                        case.inputs = earlier.inputs.clone();
                    }
                    "map" if case.op.is_empty() => case.map_before = Some(parts[1].to_owned()),
                    "map" => case.map_after = Some(parts[1].to_owned()),
                    "op" => case.op = parts[1].to_owned(),
                    "cb" => case.callbacks.push((
                        parts[2].parse().unwrap(),
                        parts[3].parse().unwrap(),
                        parts[4] == "1",
                    )),
                    "exc" => case.exception = Some(parts[1].to_owned()),
                    "ok" => {}
                    "size" => case.size = Some(parts[1].parse().unwrap()),
                    "out" => case.outputs.push(parts[2].to_owned()),
                    other => panic!("unknown record {other}"),
                }
            }
        }
    }
    cases
}

/// Expand an output whose hulls are given as `@tag` from the input feature
/// with the same tag.
fn expand_output(output: &str, inputs: &[String]) -> String {
    let Some(head) = output.strip_suffix(";hulls=@tag") else {
        return output.to_owned();
    };
    let tag = tag_of(&parse_feature(&format!("{head};hulls=0#;subs=0#")));
    let input = inputs
        .iter()
        .find(|input| tag_of(&parse_feature(input)) == tag)
        .unwrap();
    let tail = &input[input.find(";hulls=").unwrap()..];
    format!("{head}{tail}")
}

fn op_fields(op: &str) -> (&str, BTreeMap<&str, &str>) {
    let mut words = op.split(' ');
    let name = words.next().unwrap();
    let mut map = BTreeMap::new();
    for word in words {
        for part in word.split(';') {
            if let Some((key, value)) = part.split_once('=') {
                map.insert(key, value);
            }
        }
    }
    (name, map)
}

fn run_oracle_op(
    op: &str,
    map: &mut FeatureMap,
    log: &mut Vec<(i64, i64, bool)>,
) -> openms::Result<()> {
    let (name, args) = op_fields(op);
    match name {
        "filter" => {
            let mode = match args["mode"] {
                "CONVEX_HULL" => FeatureOverlapMode::ConvexHull,
                "TRACE_LEVEL" => FeatureOverlapMode::TraceLevel,
                "CENTROID_BASED" => FeatureOverlapMode::CentroidBased,
                other => panic!("unknown mode {other}"),
            };
            let by_quality = args["comparator"] == "quality";
            let rule: u32 = args["rule"].parse().unwrap();
            let tolerances = CentroidTolerances {
                rt_tolerance: bits64(args["rt"]),
                mz_tolerance: bits64(args["mz"]),
                require_same_charge: args["charge"] == "1",
                require_same_im: args["im"] == "1",
            };
            FeatureOverlapFilter::filter_with_mode(
                map,
                |left: &Feature, right: &Feature| {
                    if by_quality {
                        left.quality > right.quality
                    } else {
                        left.intensity > right.intensity
                    }
                },
                |best: &mut Feature, other: &mut Feature| {
                    let (b, o) = (tag_of(best), tag_of(other));
                    let ret = match rule {
                        0 => false,
                        2 => true,
                        _ => (b * 31 + o) % 3 != 0,
                    };
                    if rule == 3 && ret {
                        best.rt += 0.5;
                    }
                    log.push((b, o, ret));
                    ret
                },
                mode,
                &tolerances,
            )
        }
        "filter(bool=false)" => FeatureOverlapFilter::filter(
            map,
            FeatureOverlapFilter::higher_overall_quality,
            FeatureOverlapFilter::always_overlapping,
            false,
        ),
        "filter(fmap)" => FeatureOverlapFilter::filter(
            map,
            FeatureOverlapFilter::higher_overall_quality,
            FeatureOverlapFilter::always_overlapping,
            true,
        ),
        "mergeOverlappingFeatures" => FeatureOverlapFilter::merge_overlapping_features(
            map,
            bits64(args["rt"]),
            bits64(args["mz"]),
            args["charge"] == "1",
            args["im"] == "1",
            if args["intensity"] == "SUM" {
                MergeIntensityMode::Sum
            } else {
                MergeIntensityMode::Max
            },
            args["meta"] == "1",
        ),
        "mergeFAIMSFeatures" => {
            FeatureOverlapFilter::merge_faims_features(map, bits64(args["rt"]), bits64(args["mz"]))
        }
        other => panic!("unknown operation {other}"),
    }
}

#[test]
fn oracle_cases_replay_bit_for_bit() {
    let cases = oracle_cases();
    assert_eq!(cases.len(), 78);
    let mut compared = 0;
    let mut refused = 0;
    for case in &cases {
        let mut map =
            FeatureMap::from_features(case.inputs.iter().map(|f| parse_feature(f)).collect());
        if let Some(fields) = &case.map_before {
            apply_map_fields(&mut map, fields);
        }
        // The inputs survive the round trip through the port's model.
        for (input, feature) in case.inputs.iter().zip(&map.features) {
            assert_eq!(&encode_feature(feature), input, "{}: input", case.name);
        }
        let before = map.clone();
        let mut log = Vec::new();
        let result = run_oracle_op(&case.op, &mut map, &mut log);

        if case.name == "edge_hull_less_convex_hull" {
            // Debug libOpenMS aborts on the quadtree assertion; the Release replica
            // ignores the feature without a hull. The port refuses it.
            assert_eq!(case.source, "replica");
            assert_eq!(case.library_exit, 134);
            assert_eq!(case.size, Some(2));
            assert!(case.callbacks.is_empty());
            assert!(
                matches!(result, Err(Error::MissingInformation(_))),
                "{result:?}"
            );
            assert_eq!(map, before);
            refused += 1;
            continue;
        }
        if case.source == "replica" {
            assert_eq!(case.library_exit, 134, "{}", case.name);
        } else {
            assert_eq!(case.library_exit, 0, "{}", case.name);
        }
        assert_eq!(log, case.callbacks, "{}: callbacks", case.name);
        match &case.exception {
            Some(exception) => {
                match (exception.as_str(), &result) {
                    ("InvalidRange", Err(Error::InvalidRange(_)))
                    | ("MissingInformation", Err(Error::MissingInformation(_)))
                    | ("ConversionError", Err(Error::InvalidValue(_))) => {}
                    _ => panic!("{}: C++ threw {exception}, port gave {result:?}", case.name),
                }
                // The source leaves the map sorted or stripped; the port leaves it unchanged.
                assert_eq!(map, before, "{}: map after an error", case.name);
                refused += 1;
            }
            None => {
                result.unwrap_or_else(|e| panic!("{}: {e}", case.name));
                assert_eq!(Some(map.len()), case.size, "{}: size", case.name);
                assert_eq!(case.outputs.len(), map.len(), "{}: records", case.name);
                for (k, (feature, output)) in map.features.iter().zip(&case.outputs).enumerate() {
                    assert_eq!(
                        encode_feature(feature),
                        expand_output(output, &case.inputs),
                        "{}: feature {k}",
                        case.name
                    );
                }
                if let Some(after) = &case.map_after {
                    assert_eq!(&encode_map(&map), after, "{}: map fields", case.name);
                }
                compared += 1;
            }
        }
    }
    assert_eq!((compared, refused), (71, 7));
}

#[test]
fn oracle_callback_counts_cover_the_quadtree_order() {
    let cases = oracle_cases();
    let callbacks: usize = cases.iter().map(|c| c.callbacks.len()).sum();
    assert_eq!(callbacks, 16917);
    let deep = cases
        .iter()
        .find(|c| c.name == "order_centroid_loose_rule0")
        .unwrap();
    assert_eq!(deep.inputs.len(), 240);
    assert!(deep.callbacks.len() > 1000);
}

#[test]
fn oracle_trace_with_rt_min_after_rt_max_is_skipped() {
    // getFeatureBounds skips a trace whose rt_min exceeds its rt_max
    // (FeatureOverlapFilter.cpp:84-87). In the first case the second feature's
    // only traces that would overlap the first feature's are such traces, so it is
    // kept without a callback; the control has a regular trace there and is
    // removed. Both are replayed bit for bit by oracle_cases_replay_bit_for_bit.
    let cases = oracle_cases();
    let skipped = cases
        .iter()
        .find(|c| c.name == "edge_trace_inverted_bounds_skipped")
        .unwrap();
    assert_eq!(skipped.source, "library");
    assert!(skipped.callbacks.is_empty());
    assert_eq!(skipped.size, Some(2));
    let control = cases
        .iter()
        .find(|c| c.name == "edge_trace_inverted_bounds_control")
        .unwrap();
    assert_eq!(control.source, "library");
    assert_eq!(control.callbacks, vec![(0, 1, true)]);
    assert_eq!(control.size, Some(1));

    // The first feature's trace hull ends at its last scan: [9.5, 12].
    let inputs: Vec<Feature> = skipped.inputs.iter().map(|f| parse_feature(f)).collect();
    let outline = inputs[0].subordinates[0].convex_hulls[0].hull_points();
    assert_eq!(outline.first().map(|p| p.rt), Some(9.5));
    assert_eq!(outline.last().map(|p| (p.rt, p.mz)), Some((12.0, 501.0)));
    // The skipped traces: no m/z above zero (rt_min 11 from the last point,
    // rt_max 10 from the first), and m/z 0 at the first scan's lower edge.
    let traces = &inputs[1].subordinates;
    let zero = traces[0].convex_hulls[0].hull_points();
    assert!(zero.iter().all(|p| p.mz <= 0.0));
    assert_eq!((zero[0].rt, zero[zero.len() - 1].rt), (10.0, 11.0));
    let lower_zero = traces[1].convex_hulls[0].hull_points();
    assert_eq!((lower_zero[0].rt, lower_zero[0].mz), (10.0, 0.0));

    for case in [skipped, control] {
        let mut map =
            FeatureMap::from_features(case.inputs.iter().map(|f| parse_feature(f)).collect());
        let mut log = Vec::new();
        run_oracle_op(&case.op, &mut map, &mut log).unwrap();
        assert_eq!(log, case.callbacks, "{}", case.name);
        assert_eq!(Some(map.len()), case.size, "{}", case.name);
    }
}

// ---------------------------------------------------------------------------
// Tier 2: the pinned quadtree header, executed

fn parse_box(parts: &[&str]) -> QuadBox {
    QuadBox::new(
        bits32(parts[0]),
        bits32(parts[1]),
        bits32(parts[2]),
        bits32(parts[3]),
    )
}

fn replay_quadtree(text: &str) -> (usize, usize, usize) {
    let mut boxes: Vec<QuadBox> = Vec::new();
    let mut tree: Option<Quadtree<usize>> = None;
    let (mut queries, mut removals, mut pair_lists) = (0, 0, 0);
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        match parts[0] {
            "scenario" | "end" => {}
            "root" => tree = Some(Quadtree::new(parse_box(&parts[1..5]))),
            "box" => {
                assert_eq!(parts[1].parse::<usize>().unwrap(), boxes.len());
                boxes.push(parse_box(&parts[2..6]));
            }
            "add" => {
                let value: usize = parts[1].parse().unwrap();
                let get_box = |&index: &usize| boxes[index];
                tree.as_mut().unwrap().add(value, &get_box).unwrap();
            }
            "query" => {
                let query_box = parse_box(&parts[2..6]);
                let expected: Vec<usize> = match parts.get(6) {
                    Some(list) if !list.is_empty() => {
                        list.split(',').map(|v| v.parse().unwrap()).collect()
                    }
                    _ => Vec::new(),
                };
                let get_box = |&index: &usize| boxes[index];
                assert_eq!(
                    tree.as_ref().unwrap().query(query_box, &get_box),
                    expected,
                    "query {}",
                    parts[1]
                );
                queries += 1;
            }
            "pairs" => {
                let expected: Vec<(usize, usize)> = match parts.get(1) {
                    Some(list) if !list.is_empty() => list
                        .split(',')
                        .map(|pair| {
                            let (a, b) = pair.split_once(':').unwrap();
                            (a.parse().unwrap(), b.parse().unwrap())
                        })
                        .collect(),
                    _ => Vec::new(),
                };
                let get_box = |&index: &usize| boxes[index];
                assert_eq!(
                    tree.as_ref()
                        .unwrap()
                        .find_all_intersections(&get_box)
                        .unwrap(),
                    expected
                );
                pair_lists += 1;
            }
            "remove" => {
                let value: usize = parts[1].parse().unwrap();
                let get_box = |&index: &usize| boxes[index];
                tree.as_mut().unwrap().remove(&value, &get_box).unwrap();
                removals += 1;
            }
            other => panic!("unknown record {other}"),
        }
    }
    (queries, removals, pair_lists)
}

#[test]
fn quadtree_grid_replays_the_pinned_header() {
    // 359 boxes; 29 queries after the adds, after removals 50, 100 and 150, and
    // after all 186 removals; the intersection list before and after.
    let (queries, removals, pairs) = replay_quadtree(QUADTREE_GRID);
    assert_eq!((queries, removals, pairs), (145, 186, 2));
}

#[test]
fn quadtree_nan_and_infinite_boxes_follow_the_release_header() {
    let (queries, removals, pairs) = replay_quadtree(QUADTREE_NAN);
    assert_eq!((queries, removals, pairs), (33, 0, 1));
}

// ---------------------------------------------------------------------------
// Tier 4: native contracts

#[test]
fn quadtree_constants_and_box_predicates_match_the_source() {
    assert_eq!(Quadtree::<usize>::THRESHOLD, 16);
    assert_eq!(Quadtree::<usize>::MAX_DEPTH, 8);
    let a = QuadBox::new(0.0, 0.0, 2.0, 2.0);
    let touching = QuadBox::new(2.0, 0.0, 2.0, 2.0);
    assert!(!a.intersects(touching));
    assert!(a.contains(QuadBox::new(0.0, 0.0, 2.0, 2.0)));
    let point = QuadBox::new(1.0, 1.0, 0.0, 0.0);
    assert!(!point.intersects(point));
    assert_eq!(a.center(), Vector2::new(1.0, 1.0));
    assert_eq!(a.size() / 2.0 + a.top_left(), Vector2::new(1.0, 1.0));
}

#[test]
fn quadtree_remove_of_an_absent_value_is_an_error_and_changes_nothing() {
    let boxes: Vec<QuadBox> = (0..40)
        .map(|i| QuadBox::new(i as f32 * 2.0, i as f32, 1.0, 1.0))
        .collect();
    let get_box = |&index: &usize| boxes[index];
    let mut tree = Quadtree::new(QuadBox::new(0.0, 0.0, 100.0, 100.0));
    for i in 0..30 {
        tree.add(i, &get_box).unwrap();
    }
    let everything = QuadBox::new(-1.0, -1.0, 102.0, 102.0);
    let before = tree.query(everything, &get_box);
    assert!(matches!(
        tree.remove(&35, &get_box),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(tree.query(everything, &get_box), before);
    assert_eq!(tree.len(), 30);
}

fn trace_feature(rt: f64, mz: f64, quality: f32, uid: u64, with_trace: bool) -> Feature {
    let mut f = create_test_feature(rt, mz, 100.0, 2);
    f.quality = quality;
    f.unique_id = uid;
    if with_trace {
        f.subordinates = vec![create_test_feature(rt, mz, 1.0, 2)];
    }
    f
}

#[test]
fn trace_mode_candidate_without_bounds_restores_the_map_after_a_callback() {
    let mut fmap = FeatureMap::from_features(vec![
        trace_feature(100.0, 500.0, 9.0, 1, true),
        trace_feature(100.2, 500.001, 1.0, 2, false),
        // Same retention time as the first: getFeatureBounds reduces both traces
        // to rt_min == rt_max == 99, so they overlap.
        trace_feature(100.0, 500.0, 5.0, 3, true),
    ]);
    let before = fmap.clone();
    let mut calls = 0;
    let result = FeatureOverlapFilter::filter(
        &mut fmap,
        FeatureOverlapFilter::higher_overall_quality,
        |best: &mut Feature, _: &mut Feature| {
            calls += 1;
            best.quality += 100.0;
            true
        },
        true,
    );
    assert!(matches!(result, Err(Error::InvalidValue(_))), "{result:?}");
    assert_eq!(calls, 1);
    assert_eq!(fmap, before);
}

#[test]
fn faims_merge_error_after_a_merge_restores_the_map() {
    let mut a = c2_feature(100.0, 500.0, 1000.0, -45.0);
    a.unique_id = 1;
    let mut b = c2_feature(100.5, 500.01, 900.0, -60.0);
    b.unique_id = 2;
    let mut c = c2_feature(300.0, 700.0, 800.0, 0.0);
    c.metadata.insert(FAIMS_CV.into(), MetaValue::default());
    c.unique_id = 3;
    let mut d = c2_feature(300.5, 700.01, 700.0, -45.0);
    d.unique_id = 4;
    let mut e = Feature::new(50.0, 400.0, 0.0);
    e.unique_id = 5;
    let mut fmap = FeatureMap::from_features(vec![e, d, c, b, a]);
    fmap.identifier = "kept".into();
    let before = fmap.clone();
    let result = FeatureOverlapFilter::merge_faims_features(&mut fmap, 5.0, 0.05);
    assert!(matches!(result, Err(Error::InvalidValue(_))), "{result:?}");
    assert_eq!(fmap, before);
}

#[test]
fn merge_error_after_repeated_changes_of_the_same_fields_restores_the_map() {
    // The survivor absorbs three features, changing its intensity, both centroid
    // lists, FAIMS_CV, merged_centroid_IMs and FAIMS_merge_count, before the
    // fourth sum leaves the f32 range. The journal keeps only each field's value
    // before the run, and that is what the map returns to.
    let mut generator = UniqueIdGenerator::from_seed(18);
    let mut features: Vec<Feature> = (0..6)
        .map(|k| create_test_feature(100.0 + 0.1 * f64::from(k), 500.0, 0.0, 2))
        .collect();
    for f in &mut features {
        f.intensity = f32::MAX / 4.0;
    }
    set_cv(&mut features[0], -45.0);
    features[0].metadata.insert(
        MERGED_CENTROID_RTS.into(),
        MetaValue::try_from(vec![7.0]).unwrap(),
    );
    let mut fmap = map_of(features, &mut generator);
    let before = fmap.clone();
    let result = FeatureOverlapFilter::merge_overlapping_features(
        &mut fmap,
        5.0,
        0.05,
        true,
        false,
        MergeIntensityMode::Sum,
        true,
    );
    assert!(matches!(result, Err(Error::InvalidValue(_))), "{result:?}");
    assert_eq!(fmap, before);
    for (after, original) in fmap.features.iter().zip(&before.features) {
        assert_eq!(after.intensity.to_bits(), original.intensity.to_bits());
    }
}

#[test]
fn fallible_callback_error_restores_the_map() {
    let mut generator = UniqueIdGenerator::from_seed(15);
    let mut fmap = map_of(
        vec![
            create_test_feature(100.0, 500.0, 1000.0, 2),
            create_test_feature(101.0, 500.01, 500.0, 2),
            create_test_feature(102.0, 500.02, 300.0, 2),
        ],
        &mut generator,
    );
    let before = fmap.clone();
    let mut calls = 0;
    let result = FeatureOverlapFilter::filter_with_fallible_callback(
        &mut fmap,
        |l: &Feature, r: &Feature| l.intensity > r.intensity,
        |best: &mut Feature, _: &mut Feature| {
            calls += 1;
            best.intensity += 1.0;
            if calls == 2 {
                Err(Error::InvalidValue("stop".into()))
            } else {
                Ok(true)
            }
        },
        FeatureOverlapMode::CentroidBased,
        &CentroidTolerances::default(),
    );
    assert!(matches!(result, Err(Error::InvalidValue(_))));
    assert_eq!(fmap, before);
}

#[test]
fn merged_intensity_outside_f32_is_refused() {
    let mut generator = UniqueIdGenerator::from_seed(16);
    let mut fmap = map_of(
        vec![
            create_test_feature(100.0, 500.0, 0.0, 2),
            create_test_feature(101.0, 500.01, 0.0, 2),
        ],
        &mut generator,
    );
    for f in &mut fmap.features {
        f.intensity = f32::MAX * 0.75;
    }
    let before = fmap.clone();
    let result = FeatureOverlapFilter::merge_overlapping_features(
        &mut fmap,
        5.0,
        0.05,
        true,
        false,
        MergeIntensityMode::Sum,
        true,
    );
    assert!(matches!(result, Err(Error::InvalidValue(_))));
    assert_eq!(fmap, before);
    FeatureOverlapFilter::merge_overlapping_features(
        &mut fmap,
        5.0,
        0.05,
        true,
        false,
        MergeIntensityMode::Max,
        true,
    )
    .unwrap();
    assert_eq!(fmap.len(), 1);
    assert_eq!(fmap.features[0].intensity, f32::MAX * 0.75);
}

fn refused(map: &FeatureMap, mode: FeatureOverlapMode, tolerances: CentroidTolerances) -> Error {
    let mut copy = map.clone();
    let error = FeatureOverlapFilter::filter_with_mode(
        &mut copy,
        FeatureOverlapFilter::higher_overall_quality,
        FeatureOverlapFilter::always_overlapping,
        mode,
        &tolerances,
    )
    .unwrap_err();
    assert_eq!(&copy, map);
    error
}

#[test]
fn invalid_input_is_refused_before_any_change() {
    let mut generator = UniqueIdGenerator::from_seed(17);
    let good = map_of(
        vec![
            create_test_feature(100.0, 500.0, 1000.0, 2),
            create_test_feature(101.0, 500.01, 500.0, 2),
        ],
        &mut generator,
    );
    let centroid = FeatureOverlapMode::CentroidBased;
    for (rt, mz) in [(-1.0, 0.05), (5.0, f64::NAN), (f64::INFINITY, 0.05)] {
        let tolerances = CentroidTolerances {
            rt_tolerance: rt,
            mz_tolerance: mz,
            ..CentroidTolerances::default()
        };
        assert!(matches!(
            refused(&good, centroid, tolerances),
            Error::InvalidValue(_)
        ));
    }

    let mut duplicate = good.clone();
    duplicate.features[1].unique_id = duplicate.features[0].unique_id;
    assert!(matches!(
        refused(&duplicate, centroid, CentroidTolerances::default()),
        Error::InvalidValue(_)
    ));

    let mut not_finite = good.clone();
    not_finite.features[0].rt = f64::INFINITY;
    assert!(matches!(
        refused(&not_finite, centroid, CentroidTolerances::default()),
        Error::InvalidValue(_)
    ));

    let mut too_large = good.clone();
    too_large.features[0].rt = 1e300;
    assert!(matches!(
        refused(&too_large, centroid, CentroidTolerances::default()),
        Error::InvalidValue(_)
    ));

    let mut hull_less = good.clone();
    hull_less.features[1].convex_hulls.clear();
    assert!(matches!(
        refused(
            &hull_less,
            FeatureOverlapMode::ConvexHull,
            CentroidTolerances::default()
        ),
        Error::MissingInformation(_)
    ));
    let mut empty_hull = good.clone();
    empty_hull.features[1]
        .convex_hulls
        .push(ConvexHull2D::new());
    assert!(matches!(
        refused(
            &empty_hull,
            FeatureOverlapMode::TraceLevel,
            CentroidTolerances::default()
        ),
        Error::MissingInformation(_)
    ));

    let mut more_subordinates = good.clone();
    more_subordinates.features[0].subordinates = vec![
        create_test_feature(100.0, 500.0, 1.0, 2),
        create_test_feature(100.0, 500.0, 1.0, 2),
    ];
    assert!(matches!(
        refused(
            &more_subordinates,
            FeatureOverlapMode::TraceLevel,
            CentroidTolerances::default()
        ),
        Error::InvalidValue(_)
    ));

    assert!(matches!(
        refused(&FeatureMap::new(), centroid, CentroidTolerances::default()),
        Error::InvalidRange(_)
    ));
}

#[test]
fn merge_faims_features_does_not_touch_a_map_without_faims_features() {
    let mut f = create_test_feature(100.0, 500.0, 1000.0, 2);
    f.rt = f64::NAN;
    let mut fmap = FeatureMap::from_features(vec![f.clone(), f]);
    let before_bits: Vec<u64> = fmap.features.iter().map(|f| f.rt.to_bits()).collect();
    FeatureOverlapFilter::merge_faims_features(&mut fmap, 5.0, 0.05).unwrap();
    let after_bits: Vec<u64> = fmap.features.iter().map(|f| f.rt.to_bits()).collect();
    assert_eq!(before_bits, after_bits);
    assert_eq!(fmap.len(), 2);
}

#[test]
fn faims_merge_callback_matches_merge_overlapping_features() {
    let cases = oracle_cases();
    let case = cases.iter().find(|c| c.name == "order_merge_sum").unwrap();
    let input = FeatureMap::from_features(case.inputs.iter().map(|f| parse_feature(f)).collect());
    let tolerances = CentroidTolerances::default();

    let mut by_merge = input.clone();
    FeatureOverlapFilter::merge_overlapping_features(
        &mut by_merge,
        tolerances.rt_tolerance,
        tolerances.mz_tolerance,
        true,
        false,
        MergeIntensityMode::Sum,
        true,
    )
    .unwrap();
    let mut by_callback = input.clone();
    FeatureOverlapFilter::filter_with_fallible_callback(
        &mut by_callback,
        |l: &Feature, r: &Feature| l.intensity > r.intensity,
        FeatureOverlapFilter::create_faims_merge_callback(MergeIntensityMode::Sum, true),
        FeatureOverlapMode::CentroidBased,
        &tolerances,
    )
    .unwrap();
    assert_eq!(by_merge, by_callback);
    assert!(by_merge.len() < input.len());
}

#[test]
fn faims_merge_callback_keeps_best_unchanged_on_a_conversion_error() {
    let mut best = create_test_feature(100.0, 500.0, 1000.0, 2);
    set_cv(&mut best, -45.0);
    best.metadata
        .insert(MERGED_CENTROID_MZS.into(), MetaValue::from("text"));
    let other = create_test_feature(101.0, 500.01, 500.0, 2);
    let before = best.clone();
    let callback = FaimsMergeCallback::new(MergeIntensityMode::Sum, true);
    assert!(matches!(
        callback.merge(&mut best, &other),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(best, before);

    let silent = FaimsMergeCallback::new(MergeIntensityMode::Max, false);
    assert!(silent.merge(&mut best, &other).unwrap());
    assert_eq!(best.intensity, 1000.0);
    assert_eq!(best.metadata, before.metadata);
}

// ---------------------------------------------------------------------------
// Tier 4: the corrected FAIMS merge (`FaimsMergeFidelity::Corrected`, CPP-283)
//
// The corrected merge has no C++ oracle: the C++ path is broken, so no C++
// build produces the merged features of a FAIMS run. It is pinned against the
// specification derived in `FaimsMergeFidelity` — a cluster of features at
// pairwise different voltages is one analyte and collapses to the member of
// highest intensity, whose intensity is the sum of the cluster with each
// member counted once — and against the executed source behaviour it departs
// from, which the `c2_*` cases above hold.
// ---------------------------------------------------------------------------

/// The three features of `c2_three_voltages_merge_twice_into_1900_and_1700`,
/// with the unique ids the merge keys on.
fn three_voltage_cluster() -> FeatureMap {
    let mut fmap = FeatureMap::from_features(vec![
        c2_feature(100.0, 500.0, 1000.0, -45.0),
        c2_feature(100.5, 500.01, 900.0, -60.0),
        c2_feature(101.0, 500.02, 800.0, -75.0),
    ]);
    for (i, f) in fmap.features.iter_mut().enumerate() {
        f.unique_id = i as u64 + 1;
    }
    fmap
}

/// The executed source keeps 1900 and 1700 for this input — 3600 units where
/// the input held 2700, because the survivor stops absorbing after its first
/// merge and the feature it removed is offered to the next survivor
/// (`CPP-283`). The corrected merge collapses the cluster to one feature of
/// 2700.
///
/// The three features are the only ones in the map, so the quadtree's root
/// node holds all of them (below its threshold of 16) and returns them in
/// insertion order, which after the intensity sort is 1000, 900, 800: the
/// survivor absorbs 900 and then 800, and the lists follow that order.
#[test]
fn the_corrected_merge_collapses_a_three_voltage_cluster_into_one_feature() {
    let mut fmap = three_voltage_cluster();
    FeatureOverlapFilter::merge_faims_features_with_fidelity(
        &mut fmap,
        5.0,
        0.05,
        FaimsMergeFidelity::Corrected,
    )
    .unwrap();
    assert_eq!(fmap.len(), 1);
    let merged = &fmap.features[0];
    assert_eq!(merged.intensity.to_bits(), 2700.0f32.to_bits());
    assert_eq!(merged.rt, 100.0);
    assert_eq!(merged.mz, 500.0);
    assert_eq!(merged.unique_id, 1);
    assert_eq!(
        float_list(merged, MERGED_CENTROID_RTS),
        [100.0f64.to_bits(), 100.5f64.to_bits(), 101.0f64.to_bits()]
    );
    assert_eq!(
        float_list(merged, MERGED_CENTROID_MZS),
        [500.0f64.to_bits(), 500.01f64.to_bits(), 500.02f64.to_bits()]
    );
    assert_eq!(
        float_list(merged, MERGED_CENTROID_IMS),
        [
            (-45.0f64).to_bits(),
            (-60.0f64).to_bits(),
            (-75.0f64).to_bits()
        ]
    );
    assert_eq!(merged.metadata[FAIMS_MERGE_COUNT].as_i64().unwrap(), 3);
    assert!(!merged.metadata.contains_key(FAIMS_CV));
}

/// `FaimsMergeFidelity::Source` is exactly `merge_faims_features`, so the
/// executed 1900/1700 answer is still available.
#[test]
fn the_source_fidelity_is_the_executed_merge() {
    let mut corrected = three_voltage_cluster();
    let mut source = three_voltage_cluster();
    let mut plain = three_voltage_cluster();
    FeatureOverlapFilter::merge_faims_features_with_fidelity(
        &mut corrected,
        5.0,
        0.05,
        FaimsMergeFidelity::Corrected,
    )
    .unwrap();
    FeatureOverlapFilter::merge_faims_features_with_fidelity(
        &mut source,
        5.0,
        0.05,
        FaimsMergeFidelity::Source,
    )
    .unwrap();
    FeatureOverlapFilter::merge_faims_features(&mut plain, 5.0, 0.05).unwrap();
    assert_eq!(source, plain);
    let intensities: Vec<f32> = source.features.iter().map(|f| f.intensity).collect();
    assert_eq!(intensities, [1900.0, 1700.0]);
    assert_ne!(corrected, source);
    assert_eq!(FaimsMergeFidelity::default(), FaimsMergeFidelity::Corrected);
}

/// Two voltages still give the source's answer: the defect needs a third
/// feature in the cluster, so every two-voltage cluster merges the same way.
#[test]
fn the_corrected_merge_equals_the_source_on_two_voltage_clusters() {
    let features = vec![
        c2_feature(100.0, 500.0, 1000.0, -45.0),
        c2_feature(101.0, 500.01, 500.0, -60.0),
        c2_feature(900.0, 700.0, 300.0, -45.0),
    ];
    let with_ids = |features: Vec<Feature>| {
        let mut fmap = FeatureMap::from_features(features);
        for (i, f) in fmap.features.iter_mut().enumerate() {
            f.unique_id = i as u64 + 1;
        }
        fmap
    };
    let mut corrected = with_ids(features.clone());
    let mut source = with_ids(features);
    FeatureOverlapFilter::merge_faims_features_with_fidelity(
        &mut corrected,
        5.0,
        0.05,
        FaimsMergeFidelity::Corrected,
    )
    .unwrap();
    FeatureOverlapFilter::merge_faims_features(&mut source, 5.0, 0.05).unwrap();
    assert_eq!(corrected, source);
    assert_eq!(corrected.len(), 2);
    assert_eq!(
        corrected.features[0].intensity.to_bits(),
        1500.0f32.to_bits()
    );
}

/// The corrected merge still keys removal on unique ids, so it refuses the
/// repeated id 0 that `FeatureFinderAlgorithmPicked` leaves, instead of
/// erasing every feature as the source does (`CPP-282`,
/// `c2_unique_id_zero_merge_removes_every_faims_feature`). The map is
/// unchanged.
#[test]
fn the_corrected_merge_refuses_features_that_share_a_unique_id() {
    let mut fmap = FeatureMap::from_features(vec![
        c2_feature(100.0, 500.0, 1000.0, -45.0),
        c2_feature(101.0, 500.01, 500.0, -60.0),
        c2_feature(900.0, 700.0, 300.0, -45.0),
    ]);
    let before = fmap.clone();
    let error = FeatureOverlapFilter::merge_faims_features_with_fidelity(
        &mut fmap,
        5.0,
        0.05,
        FaimsMergeFidelity::Corrected,
    )
    .unwrap_err();
    match &error {
        Error::InvalidValue(message) => assert!(message.contains("CPP-282"), "{message}"),
        other => panic!("{other:?}"),
    }
    assert_eq!(fmap, before);
}

/// One FAIMS feature cannot collide with itself, so the id check does not fire
/// and the merge is the no-op both fidelities make of it.
#[test]
fn the_corrected_merge_accepts_a_single_unassigned_faims_feature() {
    let mut fmap = FeatureMap::from_features(vec![
        c2_feature(100.0, 500.0, 1000.0, -45.0),
        create_test_feature(101.0, 500.01, 500.0, 2),
    ]);
    let before = fmap.clone();
    FeatureOverlapFilter::merge_faims_features_with_fidelity(
        &mut fmap,
        5.0,
        0.05,
        FaimsMergeFidelity::Corrected,
    )
    .unwrap();
    assert_eq!(fmap, before);
}

/// A survivor does not absorb another survivor: it carries no `FAIMS_CV` any
/// more, and it is the one of higher intensity, which the specification keeps.
/// Four features within the tolerances of each other, two per voltage, are
/// therefore two analytes of two voltages each and not one analyte of four.
#[test]
fn the_corrected_merge_does_not_absorb_a_survivor() {
    let mut fmap = FeatureMap::from_features(vec![
        c2_feature(100.0, 500.0, 1000.0, -45.0),
        c2_feature(100.2, 500.005, 900.0, -60.0),
        c2_feature(100.4, 500.01, 200.0, -45.0),
        c2_feature(100.6, 500.015, 100.0, -60.0),
    ]);
    for (i, f) in fmap.features.iter_mut().enumerate() {
        f.unique_id = i as u64 + 1;
    }
    FeatureOverlapFilter::merge_faims_features_with_fidelity(
        &mut fmap,
        5.0,
        0.05,
        FaimsMergeFidelity::Corrected,
    )
    .unwrap();
    // 1000 takes 900, the first -60 it meets, and then refuses both 200 (its
    // -45 is already in the survivor's list) and 100 (so is its -60). 200 then
    // queries, refuses the 1900 survivor because that carries no `FAIMS_CV`,
    // and takes 100, which no one had removed.
    let intensities: Vec<f32> = fmap.features.iter().map(|f| f.intensity).collect();
    assert_eq!(intensities, [1900.0, 300.0]);
    for merged in &fmap.features {
        assert!(!merged.metadata.contains_key(FAIMS_CV));
        assert_eq!(
            float_list(merged, MERGED_CENTROID_IMS),
            [(-45.0f64).to_bits(), (-60.0f64).to_bits()]
        );
        assert_eq!(merged.metadata[FAIMS_MERGE_COUNT].as_i64().unwrap(), 2);
    }
}
