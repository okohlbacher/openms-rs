// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `FeatureDistance_test.cpp`, `BaseGroupFinder_test.cpp`,
//! `StablePairFinder_test.cpp`, `FeatureGroupingAlgorithm_test.cpp` and
//! `FeatureGroupingAlgorithmUnlabeled_test.cpp` at source revision
//! bc9cc12514c768385ce121d6ca4bb710fe1983c4
//! (`src/tests/class_tests/openms/source/`).
//!
//! Evidence tier 3: every expected value is a transcribed class-test literal,
//! except the sections marked "tier 4", whose expected values are derived by
//! hand from the source as shown at the assertion. The unlabeled algorithm's
//! `group` is `NOT_TESTABLE` upstream ("tested extensively in TEST/TOPP"); its
//! tier-1 differential against the retained TOPP outputs is
//! `tests/mapmatching_feature_grouping_topp.rs`. No expected value was derived
//! from this port's output.

use openms::analysis::mapmatching::base_group_finder::{BaseGroupFinder, check_ids};
use openms::analysis::mapmatching::feature_distance::FeatureDistance;
use openms::analysis::mapmatching::feature_grouping_algorithm::{
    FeatureGroupingAlgorithm, transfer_subelements,
};
use openms::analysis::mapmatching::feature_grouping_algorithm_unlabeled::FeatureGroupingAlgorithmUnlabeled;
use openms::analysis::mapmatching::stable_pair_finder::StablePairFinder;
use openms::chemistry::AASequence;
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::kernel::features::{
    BaseFeature, ColumnHeader, ConsensusFeature, ConsensusMap, Feature, FeatureHandle, FeatureMap,
};
use openms::metadata::MetaValue;
use openms::param::{Param, ParamValue};

/// `TEST_REAL_SIMILAR` with the ClassTest default relative tolerance of 1e-5.
fn similar(actual: f64, expected: f64) {
    let tolerance = 1e-5 * expected.abs().max(actual.abs()).max(1e-5);
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} is not similar to {expected}"
    );
}

fn set(param: &mut Param, key: &str, value: ParamValue) {
    let description = param.description(key).unwrap().to_string();
    let tags: Vec<String> = param.tags(key).unwrap().iter().cloned().collect();
    param.set_value(key, value, &description, &tags).unwrap();
}

fn float(param: &mut Param, key: &str, value: f64) {
    set(param, key, ParamValue::Float(value));
}

fn string(param: &mut Param, key: &str, value: &str) {
    set(param, key, ParamValue::String(value.into()));
}

fn base(rt: f64, mz: f64, intensity: f32) -> BaseFeature {
    BaseFeature::new(rt, mz, intensity)
}

// ---------------------------------------------------------------------------
// FeatureDistance_test.cpp
// ---------------------------------------------------------------------------

#[test]
fn feature_distance_operator() {
    let mut param = FeatureDistance::defaults().unwrap();
    float(&mut param, "distance_RT:max_difference", 100.0);
    float(&mut param, "distance_MZ:max_difference", 1.0);
    float(&mut param, "distance_MZ:exponent", 1.0);
    float(&mut param, "distance_intensity:weight", 1.0);
    let dist = FeatureDistance::with_param(1000.0, false, &param)
        .unwrap()
        .0;

    let mut left = base(100.0, 100.0, 100.0);
    // All distance components vary by 10% of the maximum:
    let mut right = base(110.0, 100.1, 200.0);
    let result = dist.distance(&left, &right).unwrap();
    assert!(result.0);
    similar(result.1, 0.1);
    // No differences:
    let result = dist.distance(&left, &left).unwrap();
    assert!(result.0);
    similar(result.1, 0.0);
    // Differences at maximum:
    right.rt = 200.0;
    right.mz = 101.0;
    right.intensity = 1000.0;
    left.intensity = 0.0;
    let result = dist.distance(&left, &right).unwrap();
    assert!(result.0);
    similar(result.1, 1.0);
    // Differences beyond maximum:
    right.rt = 300.0;
    let result = dist.distance(&left, &right).unwrap();
    assert!(!result.0);
    similar(result.1, 1.33333333);
    let dist2 = FeatureDistance::new(1000.0, true);
    let result = dist2.distance(&left, &right).unwrap();
    assert!(!result.0);
    assert_eq!(result.1, FeatureDistance::INFINITY);

    // ppm for m/z
    float(&mut param, "distance_intensity:weight", 0.0);
    float(&mut param, "distance_RT:weight", 0.0);
    float(&mut param, "distance_MZ:weight", 1.0);
    float(&mut param, "distance_MZ:max_difference", 10.0);
    string(&mut param, "distance_MZ:unit", "ppm");
    let dist = FeatureDistance::with_param(1000.0, false, &param)
        .unwrap()
        .0;
    left.rt = 100.0;
    left.intensity = 100.0;
    left.mz = 100.0;
    right.rt = 110.0;
    right.intensity = 200.0;
    right.mz = 100.0 + 100.0 / 1e6 * 5.0; // 5 ppm off
    let result = dist.distance(&left, &right).unwrap();
    assert!(result.0);
    similar(result.1, 0.5);

    right.mz = 100.0 + 100.0 / 1e6 * 20.0; // 20 ppm off
    let result = dist.distance(&left, &right).unwrap();
    assert!(!result.0);
    similar(result.1, 2.0);

    // charge
    string(&mut param, "ignore_charge", "false");
    let dist = FeatureDistance::with_param(1000.0, false, &param)
        .unwrap()
        .0;
    right.mz = 100.0 + 100.0 / 1e6 * 5.0;
    right.charge = 1;
    left.charge = 2;
    let result = dist.distance(&left, &right).unwrap();
    assert!(!result.0);
    assert_eq!(result.1, FeatureDistance::INFINITY);
    right.charge = 1;
    left.charge = 0;
    let result = dist.distance(&left, &right).unwrap();
    assert!(result.0);
    similar(result.1, 0.5);
    string(&mut param, "ignore_charge", "true");
    let dist = FeatureDistance::with_param(1000.0, false, &param)
        .unwrap()
        .0;
    right.charge = 1;
    left.charge = 2;
    let result = dist.distance(&left, &right).unwrap();
    assert!(result.0);
    similar(result.1, 0.5);

    // adduct
    string(&mut param, "ignore_adduct", "false");
    let dist = FeatureDistance::with_param(1000.0, false, &param)
        .unwrap()
        .0;
    right.mz = 100.0 + 100.0 / 1e6 * 5.0;
    right
        .metadata
        .insert("dc_charge_adducts".into(), MetaValue::from("NH4"));
    let result = dist.distance(&left, &right).unwrap();
    assert!(result.0);
    similar(result.1, 0.5);
    left.metadata
        .insert("dc_charge_adducts".into(), MetaValue::from("H4N"));
    let result = dist.distance(&left, &right).unwrap();
    assert!(result.0);
    similar(result.1, 0.5);
    left.metadata
        .insert("dc_charge_adducts".into(), MetaValue::from("H5N"));
    let result = dist.distance(&left, &right).unwrap();
    assert!(!result.0);
    assert_eq!(result.1, FeatureDistance::INFINITY);
    string(&mut param, "ignore_adduct", "true");
    let dist = FeatureDistance::with_param(1000.0, false, &param)
        .unwrap()
        .0;
    let result = dist.distance(&left, &right).unwrap();
    assert!(result.0);
    similar(result.1, 0.5);
}

#[test]
fn feature_distance_assignment_keeps_parameters() {
    let mut param = FeatureDistance::defaults().unwrap();
    float(&mut param, "distance_RT:max_difference", 100.0);
    float(&mut param, "distance_MZ:max_difference", 1.0);
    float(&mut param, "distance_MZ:exponent", 1.0);
    float(&mut param, "distance_intensity:weight", 1.0);
    let dist = FeatureDistance::with_param(1000.0, true, &param).unwrap().0;
    let dist2 = dist.clone();
    assert!(
        dist.to_param()
            .unwrap()
            .source_equal(&dist2.to_param().unwrap())
            .unwrap()
    );
    assert_eq!(dist, dist2);
}

/// Tier 4: the defaults the source constructor registers.
#[test]
fn feature_distance_defaults() {
    let p = FeatureDistance::defaults().unwrap();
    assert_eq!(
        p.value("distance_RT:max_difference").unwrap(),
        &ParamValue::Float(100.0)
    );
    assert_eq!(
        p.value("distance_MZ:max_difference").unwrap(),
        &ParamValue::Float(0.3)
    );
    assert_eq!(
        p.value("distance_MZ:unit").unwrap(),
        &ParamValue::String("Da".into())
    );
    assert_eq!(
        p.value("distance_MZ:exponent").unwrap(),
        &ParamValue::Float(2.0)
    );
    assert_eq!(
        p.value("distance_intensity:weight").unwrap(),
        &ParamValue::Float(0.0)
    );
    assert_eq!(
        p.value("ignore_charge").unwrap(),
        &ParamValue::String("false".into())
    );
    assert_eq!(
        p.value("ignore_adduct").unwrap(),
        &ParamValue::String("true".into())
    );
    // A restriction violation is refused (source Exception::InvalidParameter).
    let mut bad = p.clone();
    float(&mut bad, "distance_RT:weight", -1.0);
    assert!(FeatureDistance::from_param(&bad).is_err());
}

/// Tier 4: log transformation compares `log10(x + 1)` and normalises by
/// `log10(max_intensity + 1)`: with max 999, |log10(100) - log10(10)| / 3 =
/// 1/3, weighted alone.
#[test]
fn feature_distance_log_transform() {
    let mut param = FeatureDistance::defaults().unwrap();
    float(&mut param, "distance_RT:weight", 0.0);
    float(&mut param, "distance_MZ:weight", 0.0);
    float(&mut param, "distance_intensity:weight", 1.0);
    string(&mut param, "distance_intensity:log_transform", "enabled");
    let dist = FeatureDistance::with_param(999.0, false, &param).unwrap().0;
    let result = dist
        .distance(&base(1.0, 1.0, 99.0), &base(1.0, 1.0, 9.0))
        .unwrap();
    assert!(result.0);
    similar(result.1, 1.0 / 3.0);
}

// ---------------------------------------------------------------------------
// BaseGroupFinder_test.cpp: the constructor and run are NOT_TESTABLE; the
// protected checkIds_ is exercised here (tier 4, from its body).
// ---------------------------------------------------------------------------

struct TestPairFinder;

impl BaseGroupFinder for TestPairFinder {
    fn run(&self, _input: &[ConsensusMap], _result: &mut ConsensusMap) -> openms::Result<()> {
        Ok(())
    }
}

#[test]
fn base_group_finder_check_ids() {
    let mut a = ConsensusMap::new();
    a.column_headers.insert(0, ColumnHeader::default());
    let mut b = ConsensusMap::new();
    b.column_headers.insert(1, ColumnHeader::default());
    assert!(check_ids(&[a.clone(), b.clone()]).is_ok());
    b.column_headers.insert(0, ColumnHeader::default());
    assert!(check_ids(&[a, b]).is_err());
    let mut result = ConsensusMap::new();
    TestPairFinder.run(&[], &mut result).unwrap();
}

// ---------------------------------------------------------------------------
// StablePairFinder_test.cpp
// ---------------------------------------------------------------------------

fn feature(rt: f64, mz: f64, intensity: f32, uid: u64) -> BaseFeature {
    let mut f = base(rt, mz, intensity);
    f.unique_id = uid;
    f
}

fn cons(map_index: u64, f: &BaseFeature) -> ConsensusFeature {
    ConsensusFeature::from_feature(map_index, f).unwrap()
}

fn same_handle(h: &FeatureHandle, expected: &FeatureHandle) {
    assert_eq!(h.key(), expected.key());
    assert_eq!(h.rt, expected.rt);
    assert_eq!(h.mz, expected.mz);
    assert_eq!(h.intensity, expected.intensity);
}

#[test]
fn stable_pair_finder_run() {
    let feat1 = feature(0.0, 0.0, 100.0, 0);
    let feat2 = feature(200.0, 300.0, 300.0, 1);
    let feat3 = feature(400.0, 500.0, 400.0, 2);
    let feat4 = feature(4.0, 0.04, 100.0, 0);
    let feat5 = feature(204.0, 300.04, 300.0, 1);
    let feat6 = feature(404.0, 500.04, 400.0, 2);
    let mut input = vec![ConsensusMap::new(), ConsensusMap::new()];
    for f in [&feat1, &feat2, &feat3] {
        input[0].features.push(cons(0, f));
    }
    for f in [&feat4, &feat5, &feat6] {
        input[1].features.push(cons(1, f));
    }
    let spf = StablePairFinder::from_param(&StablePairFinder::defaults().unwrap()).unwrap();
    let mut result = ConsensusMap::new();
    spf.run(&input, &mut result).unwrap();
    assert_eq!(result.len(), 3);
    let groups: Vec<&[FeatureHandle]> = result.features.iter().map(|c| c.handles()).collect();
    same_handle(&groups[0][0], &FeatureHandle::new(0, &feat1));
    same_handle(&groups[0][1], &FeatureHandle::new(1, &feat4));
    same_handle(&groups[1][0], &FeatureHandle::new(0, &feat2));
    same_handle(&groups[1][1], &FeatureHandle::new(1, &feat5));
    same_handle(&groups[2][0], &FeatureHandle::new(0, &feat3));
    same_handle(&groups[2][1], &FeatureHandle::new(1, &feat6));
}

#[test]
fn stable_pair_finder_quality() {
    let feat1 = feature(100.0, 100.0, 100.0, 0);
    let feat2 = feature(200.0, 200.0, 200.0, 1);
    let feat3 = feature(300.0, 300.0, 300.0, 2);
    let mut param = StablePairFinder::defaults().unwrap();
    float(&mut param, "distance_RT:max_difference", 1000.0);
    float(&mut param, "distance_MZ:max_difference", 1000.0);
    float(&mut param, "second_nearest_gap", 2.0);
    let spf = StablePairFinder::from_param(&param).unwrap();
    let mut result = ConsensusMap::new();

    // best case:
    let mut input = vec![ConsensusMap::new(), ConsensusMap::new()];
    input[0].features.push(cons(0, &feat1));
    input[1].features.push(cons(1, &feat1));
    spf.run(&input, &mut result).unwrap();
    assert_eq!(result.len(), 1);
    assert_eq!(result.features[0].len(), 2);
    assert_eq!(result.features[0].quality, 1.0);
    input[0] = result.clone();
    input[1].features[0] = cons(2, &feat1);
    spf.run(&input, &mut result).unwrap();
    assert_eq!(result.len(), 1);
    assert_eq!(result.features[0].len(), 3);
    assert_eq!(result.features[0].quality, 1.0);

    // worst case:
    input[0].clear(true);
    spf.run(&input, &mut result).unwrap();
    assert_eq!(result.len(), 1);
    assert_eq!(result.features[0].len(), 1);
    assert_eq!(result.features[0].quality, 0.0);

    // intermediate cases; basis: feat1 < feat2 < feat3
    input[1].clear(true);
    input[0].features.push(cons(0, &feat1));
    input[1].features.push(cons(1, &feat2));
    spf.run(&input, &mut result).unwrap();
    let cons1 = result.features[0].clone();
    assert_eq!(cons1.len(), 2);
    input[0] = result.clone();
    input[1].features[0] = cons(2, &feat3);
    spf.run(&input, &mut result).unwrap();
    let cons2 = result.features[0].clone();
    assert_eq!(cons2.len(), 3);
    assert!(cons1.quality > 0.0);
    assert!(cons2.quality > 0.0);
    assert!(cons1.quality < 1.0);
    assert!(cons2.quality < 1.0);
    // quality(feat1, feat2) > quality((feat1, feat2), feat3):
    assert!(cons1.quality > cons2.quality);
    input[0].clear(true);
    input[0].features.push(cons(1, &feat2));
    spf.run(&input, &mut result).unwrap();
    let cons3 = result.features[0].clone();
    // quality(feat2, feat3) > quality((feat1, feat2), feat3):
    assert!(cons3.quality > cons2.quality);
    input[0].features[0] = cons(0, &feat1);
    spf.run(&input, &mut result).unwrap();
    let cons4 = result.features[0].clone();
    // quality(feat1, feat3) < quality((feat1, feat2), feat3):
    assert!(cons4.quality < cons2.quality);
}

/// Tier 4: exactly two inputs are required, and shared column identifiers
/// are refused (source `Exception::IllegalArgument`).
#[test]
fn stable_pair_finder_input_checks() {
    let spf = StablePairFinder::default();
    let mut result = ConsensusMap::new();
    assert!(spf.run(&[ConsensusMap::new()], &mut result).is_err());
    let mut a = ConsensusMap::new();
    a.column_headers.insert(0, ColumnHeader::default());
    assert!(spf.run(&[a.clone(), a], &mut result).is_err());
}

/// Tier 4: with `use_identifications`, features whose best hits differ never
/// pair; a feature without identifications pairs with anything.
#[test]
fn stable_pair_finder_identifications() {
    let identified = |uid: u64, map: u64, sequence: &str| {
        let mut f = feature(100.0, 100.0, 100.0, uid);
        let mut identification = PeptideIdentification::default();
        identification.hits.push(PeptideHit {
            sequence: AASequence::parse(sequence).unwrap(),
            ..PeptideHit::default()
        });
        f.peptide_identifications.push(identification);
        cons(map, &f)
    };
    let mut param = StablePairFinder::defaults().unwrap();
    string(&mut param, "use_identifications", "true");
    let spf = StablePairFinder::from_param(&param).unwrap();
    let mut result = ConsensusMap::new();
    let mut input = vec![ConsensusMap::new(), ConsensusMap::new()];
    input[0].features.push(identified(1, 0, "PEPTIDE"));
    input[1].features.push(identified(2, 1, "PEPTIDER"));
    spf.run(&input, &mut result).unwrap();
    assert_eq!(result.len(), 2);
    input[1].features[0] = identified(2, 1, "PEPTIDE");
    spf.run(&input, &mut result).unwrap();
    assert_eq!(result.len(), 1);
    // The merged feature carries both features' identifications.
    assert_eq!(result.features[0].peptide_identifications.len(), 2);
    input[1].features[0] = cons(1, &feature(100.0, 100.0, 100.0, 2));
    spf.run(&input, &mut result).unwrap();
    assert_eq!(result.len(), 1);
}

// ---------------------------------------------------------------------------
// FeatureGroupingAlgorithm_test.cpp
// ---------------------------------------------------------------------------

struct Fga;

impl FeatureGroupingAlgorithm for Fga {
    fn group(&mut self, _maps: &[FeatureMap], map: &mut ConsensusMap) -> openms::Result<()> {
        let header = map.column_headers.entry(0).or_default();
        header.filename = "bla".into();
        header.size = 5;
        Ok(())
    }
}

#[test]
fn feature_grouping_algorithm_group() {
    let mut map = ConsensusMap::new();
    Fga.group(&[], &mut map).unwrap();
    assert_eq!(map.column_headers[&0].filename, "bla");
}

#[test]
fn feature_grouping_algorithm_transfer_subelements() {
    let mut maps = vec![ConsensusMap::new(), ConsensusMap::new()];
    for (m, names) in [(0usize, ["file1", "file2"]), (1, ["file3", "file4"])] {
        for (k, name) in names.iter().enumerate() {
            maps[m].column_headers.insert(
                k as u64,
                ColumnHeader {
                    filename: (*name).into(),
                    size: 1,
                    ..ColumnHeader::default()
                },
            );
        }
    }
    let (feat1, feat2, feat3, feat4) = (
        BaseFeature::default(),
        BaseFeature::default(),
        BaseFeature::default(),
        BaseFeature::default(),
    );
    let handle1 = FeatureHandle::new(0, &feat1);
    let handle2 = FeatureHandle::new(1, &feat2);
    let mut handle3 = FeatureHandle::new(0, &feat3);
    let mut handle4 = FeatureHandle::new(1, &feat4);
    let mut c0 = ConsensusFeature::new();
    c0.insert(handle1).unwrap();
    c0.insert(handle2).unwrap();
    c0.unique_id = 1;
    maps[0].features.push(c0);
    let mut c1 = ConsensusFeature::new();
    c1.insert(handle3).unwrap();
    c1.insert(handle4).unwrap();
    c1.unique_id = 2;
    maps[1].features.push(c1);

    let mut out = ConsensusMap::new();
    let handle5 = FeatureHandle::new(0, &maps[0].features[0].base);
    let handle6 = FeatureHandle::new(1, &maps[1].features[0].base);
    let mut o = ConsensusFeature::new();
    o.insert(handle5).unwrap();
    o.insert(handle6).unwrap();
    out.features.push(o);

    transfer_subelements(&maps, &mut out).unwrap();
    assert_eq!(out.column_headers.len(), 4);
    assert_eq!(out.column_headers[&0].filename, "file1");
    assert_eq!(out.column_headers[&3].filename, "file4");
    assert_eq!(out.len(), 1);
    assert_eq!(out.features[0].len(), 4);
    handle3.map_index = 2;
    handle4.map_index = 3;
    let group = out.features[0].handles();
    assert_eq!(group[0], handle1);
    assert_eq!(group[1], handle2);
    assert_eq!(group[2], handle3);
    assert_eq!(group[3], handle4);
}

/// Tier 4: `old_map_index` renumbering, and removal of `map_index` without it.
#[test]
fn feature_grouping_algorithm_transfer_reindexes_identifications() {
    let mut maps = vec![ConsensusMap::new(), ConsensusMap::new()];
    maps[0].column_headers.insert(0, ColumnHeader::default());
    maps[1].column_headers.insert(0, ColumnHeader::default());
    maps[1].column_headers.insert(5, ColumnHeader::default());
    let mut out = ConsensusMap::new();
    let mut with_old = PeptideIdentification::default();
    with_old
        .metadata
        .insert("map_index".into(), MetaValue::try_from(1u64).unwrap());
    with_old
        .metadata
        .insert("old_map_index".into(), MetaValue::try_from(5u64).unwrap());
    let mut without_old = PeptideIdentification::default();
    without_old
        .metadata
        .insert("map_index".into(), MetaValue::try_from(1u64).unwrap());
    out.unassigned_peptide_identifications = vec![with_old, without_old];
    transfer_subelements(&maps, &mut out).unwrap();
    // Columns: (0,0) -> 0, (1,0) -> 1, (1,5) -> 2.
    let ids = &out.unassigned_peptide_identifications;
    assert_eq!(
        ids[0].metadata.get("map_index").unwrap().as_i64().unwrap(),
        2
    );
    assert!(ids[0].metadata.get("old_map_index").is_none());
    assert!(ids[1].metadata.get("map_index").is_none());
}

// ---------------------------------------------------------------------------
// FeatureGroupingAlgorithmUnlabeled_test.cpp (group NOT_TESTABLE upstream;
// tier 4 here, tier 1 in the TOPP differential)
// ---------------------------------------------------------------------------

fn feature_map(features: &[(f64, f64, f32, u64)]) -> FeatureMap {
    FeatureMap::from_features(
        features
            .iter()
            .map(|&(rt, mz, intensity, uid)| {
                let mut f = Feature::new(rt, mz, intensity);
                f.unique_id = uid;
                f
            })
            .collect(),
    )
}

/// Map 1 has the most features and is the reference; map 0 then map 2 are
/// paired into it. The feature at (100, 500) occurs in all three maps, the
/// one at (300, 700) only in map 1, and (200, 600) in maps 0 and 1.
#[test]
fn unlabeled_group_links_three_maps() {
    let maps = vec![
        feature_map(&[(100.0, 500.0, 10.0, 1), (200.0, 600.0, 20.0, 2)]),
        feature_map(&[
            (101.0, 500.01, 11.0, 11),
            (201.0, 600.01, 21.0, 12),
            (300.0, 700.0, 30.0, 13),
        ]),
        feature_map(&[(99.0, 499.99, 9.0, 21)]),
    ];
    let mut algorithm = FeatureGroupingAlgorithmUnlabeled::default();
    let mut out = ConsensusMap::new();
    out.column_headers.insert(
        7,
        ColumnHeader {
            filename: "kept".into(),
            ..ColumnHeader::default()
        },
    );
    algorithm.group(&maps, &mut out).unwrap();
    assert_eq!(out.column_headers[&7].filename, "kept");
    assert_eq!(out.len(), 3);
    // Sorted by size descending first.
    let keys: Vec<Vec<(u64, u64)>> = out
        .features
        .iter()
        .map(|c| c.handles().iter().map(FeatureHandle::key).collect())
        .collect();
    assert_eq!(keys[0], vec![(0, 1), (1, 11), (2, 21)]);
    assert_eq!(keys[1], vec![(0, 2), (1, 12)]);
    assert_eq!(keys[2], vec![(1, 13)]);
    assert_eq!(out.features[2].quality, 0.0);
    assert!(algorithm.group(&maps[..1], &mut out).is_err());
}

/// The incremental path gives the same groups, unsorted and without
/// post-processing.
#[test]
fn unlabeled_incremental_path() {
    let reference = feature_map(&[(100.0, 500.0, 10.0, 1), (200.0, 600.0, 20.0, 2)]);
    let other = feature_map(&[(101.0, 500.01, 11.0, 11)]);
    let mut algorithm = FeatureGroupingAlgorithmUnlabeled::default();
    algorithm.set_reference(0, &reference).unwrap();
    assert!(algorithm.add_to_group(0, &other).is_err());
    algorithm.add_to_group(1, &other).unwrap();
    let result = algorithm.result_map();
    assert_eq!(result.len(), 2);
    // Sorted by m/z only: the pair (mean m/z 500.005) comes first.
    assert_eq!(result.features[0].len(), 2);
    assert_eq!(result.features[1].handles()[0].key(), (0, 2));
}
