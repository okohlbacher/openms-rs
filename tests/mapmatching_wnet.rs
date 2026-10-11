// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `WNetMatcher_test.cpp` and `FeatureGroupingAlgorithmWNet_test.cpp`
//! at source revision bc9cc12514c768385ce121d6ca4bb710fe1983c4
//! (`src/tests/class_tests/openms/source/`), plus checks of the parameter
//! tree and of the refusals this port adds.
//!
//! Evidence tier 3: every expected value is a transcribed class-test literal.
//! Sections marked "tier 4" are native checks derived from the model's
//! definition, not from this port's output. The differential against the
//! executed C++ is `tests/mapmatching_wnet_oracle.rs`; the retained TOPP
//! output is `tests/mapmatching_wnet_topp.rs`.

use std::collections::BTreeSet;

use openms::analysis::mapmatching::feature_grouping_algorithm::FeatureGroupingAlgorithm;
use openms::analysis::mapmatching::feature_grouping_algorithm_wnet::FeatureGroupingAlgorithmWNet;
use openms::analysis::mapmatching::wnet_aligner::{Spectrum, WNetAligner};
use openms::analysis::mapmatching::wnet_matcher::{DistanceMetric, WNetMatcher};
use openms::kernel::features::{ConsensusFeature, ConsensusMap, Feature, FeatureMap};
use openms::param::{Param, ParamValue};

fn set(param: &mut Param, key: &str, value: ParamValue) {
    let description = param.description(key).unwrap().to_string();
    let tags: Vec<String> = param.tags(key).unwrap().iter().cloned().collect();
    param.set_value(key, value, &description, &tags).unwrap();
}

fn pair_set(pairs: &[(usize, usize)]) -> BTreeSet<(usize, usize)> {
    pairs.iter().copied().collect()
}

// ---- WNetMatcher_test.cpp -------------------------------------------------

#[test]
fn metric_from_string() {
    assert_eq!(WNetMatcher::metric_from_string("L1"), DistanceMetric::L1);
    assert_eq!(WNetMatcher::metric_from_string("L2"), DistanceMetric::L2);
    assert_eq!(
        WNetMatcher::metric_from_string("LINF"),
        DistanceMetric::Linf
    );
    // Unknown strings fall back to LINF and emit a warning.
    let (metric, warning) = WNetMatcher::metric_from_string_with_warning("unknown");
    assert_eq!(metric, DistanceMetric::Linf);
    assert_eq!(
        warning.as_deref(),
        Some("WNetMatcher::metricFromString: unrecognized metric 'unknown', falling back to LINF.")
    );
    assert_eq!(WNetMatcher::metric_from_string(""), DistanceMetric::Linf);
    assert!(
        WNetMatcher::metric_from_string_with_warning("L1")
            .1
            .is_none()
    );
}

#[test]
fn match_empty_inputs() {
    let pos_b = [[1.0, 2.0]];
    let int_b = [1.0];
    let result = WNetMatcher::match_with_defaults(&[], &[], &pos_b, &int_b).unwrap();
    assert!(result.matched_pairs.is_empty());
    assert_eq!(result.cost, 0.0);
    let result = WNetMatcher::match_with_defaults(&pos_b, &int_b, &[], &[]).unwrap();
    assert!(result.matched_pairs.is_empty());
    assert_eq!(result.cost, 0.0);
}

#[test]
fn match_identical_point_sets() {
    let positions = [[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]];
    let intensities = [100.0, 200.0, 150.0];
    let result =
        WNetMatcher::match_with_defaults(&positions, &intensities, &positions, &intensities)
            .unwrap();
    assert_eq!(result.matched_pairs.len(), 3);
    assert!(result.cost.abs() < 1e-12, "cost {}", result.cost);
    // Each point must match its own index.
    let pairs = pair_set(&result.matched_pairs);
    assert!(pairs.contains(&(0, 0)));
    assert!(pairs.contains(&(1, 1)));
    assert!(pairs.contains(&(2, 2)));
}

#[test]
fn match_close_points_match_distant_points_do_not() {
    let pos_a = [[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]];
    let int_a = [100.0, 200.0, 150.0];
    let pos_b = [[1.1, 10.1], [2.1, 20.1], [50.0, 500.0]];
    let int_b = [110.0, 190.0, 50.0];
    let result = WNetMatcher::match_points(
        &pos_a,
        &int_a,
        &pos_b,
        &int_b,
        DistanceMetric::Linf,
        5.0,
        5.0,
    )
    .unwrap();
    // The first two pairs are close (LINF distance 0.1), the third is far.
    assert_eq!(result.matched_pairs.len(), 2);
    assert!(result.cost > 0.0);
    let pairs = pair_set(&result.matched_pairs);
    assert!(pairs.contains(&(0, 0)));
    assert!(pairs.contains(&(1, 1)));
}

#[test]
fn match_refuses_mismatched_lengths_and_bad_costs() {
    let pos = [[1.0, 1.0], [2.0, 2.0]];
    assert!(WNetMatcher::match_with_defaults(&pos, &[1.0], &pos, &[1.0, 1.0]).is_err());
    assert!(WNetMatcher::match_with_defaults(&pos, &[1.0, 1.0], &pos, &[1.0]).is_err());
    // A negative trash cost, and intensities that do not sum above zero.
    assert!(
        WNetMatcher::match_points(
            &pos,
            &[1.0, 1.0],
            &pos,
            &[1.0, 1.0],
            DistanceMetric::Linf,
            10.0,
            -1.0
        )
        .is_err()
    );
    assert!(WNetMatcher::match_with_defaults(&pos, &[0.0, 0.0], &pos, &[0.0, 0.0]).is_err());
    // A position that scaling makes non-finite never reaches the solver.
    assert!(
        WNetMatcher::match_with_defaults(&[[f64::NAN, 0.0]], &[1.0], &pos[..1], &[1.0]).is_err()
    );
}

// ---- Tier 4: the model's definition ---------------------------------------

/// One point on each side at distance d <= max_distance with equal weight:
/// they match, and the cost is `d` once the scaling is undone (the single
/// unit of flow travels `d`, instead of paying the trash cost).
#[test]
fn single_pair_cost_is_the_distance() {
    let a = Spectrum::new(vec![[0.0, 0.0]], vec![2.0]).unwrap();
    let b = Spectrum::new(vec![[3.0, 4.0]], vec![2.0]).unwrap();
    // L2 distance 5, L1 7, LINF 4.
    for (metric, distance) in [
        (DistanceMetric::L2, 5.0),
        (DistanceMetric::L1, 7.0),
        (DistanceMetric::Linf, 4.0),
    ] {
        let mut aligner = WNetAligner::new(&a, &[&b], metric, 10.0, 20.0).unwrap();
        aligner.set_point(&[1.0]).unwrap();
        let (e, t) = aligner.consensus_for_target(0).unwrap();
        assert_eq!((e, t), (vec![0], vec![0]), "{metric:?}");
        let cost = aligner.total_cost().unwrap();
        // Two units of weight travel `distance`; the truncations of the
        // scaled flow and distance change that by far less than 1e-6.
        assert!((cost - 2.0 * distance).abs() < 1e-6, "{metric:?}: {cost}");
    }
}

/// A pair farther apart than `max_distance` has no arc: both go unmatched,
/// the single component's flow passes through the trash arc and costs
/// nothing in a network without arcs.
#[test]
fn pairs_beyond_max_distance_are_not_connected() {
    let a = Spectrum::new(vec![[0.0, 0.0]], vec![1.0]).unwrap();
    let b = Spectrum::new(vec![[100.0, 0.0]], vec![1.0]).unwrap();
    let mut aligner = WNetAligner::new(&a, &[&b], DistanceMetric::Linf, 10.0, 10.0).unwrap();
    assert_eq!(aligner.no_subgraphs(), 0);
    aligner.set_point(&[1.0]).unwrap();
    let (e, t) = aligner.consensus_for_target(0).unwrap();
    assert!(e.is_empty() && t.is_empty());
    assert_eq!(aligner.total_cost().unwrap(), 0.0);
}

/// Two empirical points compete for one theoretical point: the nearer wins
/// the flow, and the consensus (one pair per peak) keeps only it.
#[test]
fn nearer_point_wins_the_competition() {
    let a = Spectrum::new(vec![[0.0, 0.0], [4.0, 0.0]], vec![1.0, 1.0]).unwrap();
    let b = Spectrum::new(vec![[1.0, 0.0]], vec![1.0]).unwrap();
    let mut aligner = WNetAligner::new(&a, &[&b], DistanceMetric::Linf, 10.0, 10.0).unwrap();
    aligner.set_point(&[1.0]).unwrap();
    let (e, t, f) = aligner.flows_for_target(0).unwrap();
    // Total empirical intensity 2 against theoretical 1: only the nearer
    // empirical point's unit travels, the other goes to the trash.
    assert_eq!(e, vec![0]);
    assert_eq!(t, vec![0]);
    assert_eq!(f.len(), 1);
    let (ce, ct) = aligner.consensus_for_target(0).unwrap();
    assert_eq!((ce, ct), (vec![0], vec![0]));
}

#[test]
fn aligner_refusals() {
    let a = Spectrum::new(vec![[0.0, 0.0]], vec![1.0]).unwrap();
    let b = Spectrum::new(vec![[0.0, 0.0]], vec![1.0]).unwrap();
    let empty = Spectrum::new(vec![], vec![]).unwrap();
    assert!(Spectrum::new(vec![[0.0, 0.0]], vec![]).is_err());
    assert!(WNetAligner::new(&a, &[], DistanceMetric::L1, 1.0, 1.0).is_err());
    assert!(WNetAligner::new(&empty, &[&b], DistanceMetric::L1, 1.0, 1.0).is_err());
    assert!(WNetAligner::new(&a, &[&empty], DistanceMetric::L1, 1.0, 1.0).is_err());
    assert!(WNetAligner::new(&a, &[&b], DistanceMetric::L1, 1.0, -1.0).is_err());
    // trash_cost 0: the scale factor is infinite.
    assert!(WNetAligner::new(&a, &[&b], DistanceMetric::L1, 1.0, 0.0).is_err());
    let mut aligner = WNetAligner::new(&a, &[&b], DistanceMetric::L1, 1.0, 1.0).unwrap();
    assert!(aligner.total_cost().is_err());
    assert!(aligner.flows_for_target(0).is_err());
    assert!(aligner.set_point(&[1.0, 2.0]).is_err());
    assert!(aligner.set_point(&[f64::NAN]).is_err());
    assert!(aligner.set_point(&[1.0]).is_ok());
    assert!(
        Spectrum::new(vec![[0.0, 0.0]], vec![0.0])
            .unwrap()
            .normalized()
            .is_err()
    );
}

// ---- FeatureGroupingAlgorithmWNet_test.cpp --------------------------------

fn feature(unique_id: u64, mz: f64, rt: f64, intensity: f32) -> Feature {
    let mut f = Feature::new(rt, mz, intensity);
    f.unique_id = unique_id;
    f
}

fn sizes(result: &ConsensusMap) -> (usize, usize) {
    let pairs = result.features.iter().filter(|f| f.len() == 2).count();
    let singletons = result.features.iter().filter(|f| f.len() == 1).count();
    (pairs, singletons)
}

#[test]
fn group_feature_maps_requires_two_maps() {
    let mut algorithm = FeatureGroupingAlgorithmWNet::default();
    let mut out = ConsensusMap::new();
    assert!(algorithm.group(&[], &mut out).is_err());
    assert!(algorithm.group(&[FeatureMap::new()], &mut out).is_err());
    assert!(algorithm.group_consensus(&[], &mut out).is_err());
}

fn class_test_maps() -> Vec<FeatureMap> {
    // Map 1: features at (mz=500, rt=100), (mz=600, rt=200), (mz=700, rt=300)
    // Map 2: features at (mz=500.1, rt=101), (mz=600.1, rt=201), (mz=900, rt=900)
    vec![
        FeatureMap::from_features(vec![
            feature(1, 500.0, 100.0, 1000.0),
            feature(2, 600.0, 200.0, 2000.0),
            feature(3, 700.0, 300.0, 1500.0),
        ]),
        FeatureMap::from_features(vec![
            feature(4, 500.1, 101.0, 1100.0),
            feature(5, 600.1, 201.0, 1900.0),
            feature(6, 900.0, 900.0, 500.0),
        ]),
    ]
}

#[test]
fn group_feature_maps_da_mode() {
    let mut param = FeatureGroupingAlgorithmWNet::defaults().unwrap();
    set(&mut param, "mz_unit", ParamValue::String("Da".into()));
    let mut algorithm = FeatureGroupingAlgorithmWNet::from_param(&param).unwrap();
    let maps = class_test_maps();
    let mut result = ConsensusMap::new();
    algorithm.group(&maps, &mut result).unwrap();
    // Two matched pairs and two singletons.
    assert_eq!(result.len(), 4);
    assert_eq!(sizes(&result), (2, 2));
    // The pairs link unique ids (1, 4) and (2, 5).
    let linked: BTreeSet<Vec<(u64, u64)>> = result
        .features
        .iter()
        .filter(|f| f.len() == 2)
        .map(|f| f.handles().iter().map(|h| h.key()).collect())
        .collect();
    assert!(linked.contains(&vec![(0, 1), (1, 4)]));
    assert!(linked.contains(&vec![(0, 2), (1, 5)]));
}

#[test]
fn group_feature_maps_ppm_mode() {
    let mut param = FeatureGroupingAlgorithmWNet::defaults().unwrap();
    set(&mut param, "mz_unit", ParamValue::String("ppm".into()));
    // 500 ppm at 500 Da is about 0.25 Da.
    set(&mut param, "max_mz_shift_ppm", ParamValue::Float(500.0));
    let mut algorithm = FeatureGroupingAlgorithmWNet::from_param(&param).unwrap();
    let mut result = ConsensusMap::new();
    algorithm.group(&class_test_maps(), &mut result).unwrap();
    assert_eq!(result.len(), 4);
    assert_eq!(sizes(&result), (2, 2));
}

#[test]
fn group_feature_maps_decharge_mode() {
    // z = 2 at 250.5 projects to [M+H]+ of about 499.993, which is within
    // 0.3 Da of the z = 1 feature at 500.0.
    let mut a = feature(11, 250.5, 100.0, 1000.0);
    a.charge = 2;
    let mut b = feature(12, 500.0, 100.0, 900.0);
    b.charge = 1;
    let maps = vec![
        FeatureMap::from_features(vec![a]),
        FeatureMap::from_features(vec![b]),
    ];
    let mut param = FeatureGroupingAlgorithmWNet::defaults().unwrap();
    set(&mut param, "mz_unit", ParamValue::String("Da".into()));
    set(&mut param, "decharge_mz", ParamValue::String("true".into()));
    let mut algorithm = FeatureGroupingAlgorithmWNet::from_param(&param).unwrap();
    let mut result = ConsensusMap::new();
    algorithm.group(&maps, &mut result).unwrap();
    assert_eq!(sizes(&result).0, 1);

    // Without decharging the two are 249 Da apart: no pair.
    let mut param = FeatureGroupingAlgorithmWNet::defaults().unwrap();
    set(&mut param, "mz_unit", ParamValue::String("Da".into()));
    let mut algorithm = FeatureGroupingAlgorithmWNet::from_param(&param).unwrap();
    let mut result = ConsensusMap::new();
    algorithm.group(&maps, &mut result).unwrap();
    assert_eq!(sizes(&result).0, 0);
}

#[test]
fn group_consensus_maps() {
    let mut param = FeatureGroupingAlgorithmWNet::defaults().unwrap();
    set(&mut param, "mz_unit", ParamValue::String("Da".into()));
    let mut algorithm = FeatureGroupingAlgorithmWNet::from_param(&param).unwrap();
    let consensus = |id: u64, mz: f64, rt: f64, intensity: f32| {
        let mut f = ConsensusFeature::new();
        f.unique_id = id;
        f.mz = mz;
        f.rt = rt;
        f.intensity = intensity;
        f
    };
    let mut maps = vec![ConsensusMap::new(), ConsensusMap::new()];
    maps[0].features.push(consensus(101, 400.0, 50.0, 500.0));
    maps[0].features.push(consensus(102, 500.0, 150.0, 800.0));
    maps[1].features.push(consensus(103, 400.1, 51.0, 520.0));
    maps[1].features.push(consensus(104, 900.0, 900.0, 100.0));
    let mut result = ConsensusMap::new();
    algorithm.group_consensus(&maps, &mut result).unwrap();
    // cf1 and cf3 match; cf2 and cf4 stay singletons.
    assert_eq!(result.len(), 3);
    assert_eq!(sizes(&result), (1, 2));
}

#[test]
fn group_leaves_empty_maps_alone() {
    let mut algorithm = FeatureGroupingAlgorithmWNet::default();
    let maps = vec![
        FeatureMap::new(),
        FeatureMap::from_features(vec![feature(1, 500.0, 10.0, 5.0)]),
        FeatureMap::from_features(vec![feature(2, 500.0, 10.0, 7.0)]),
    ];
    let mut result = ConsensusMap::new();
    algorithm.group(&maps, &mut result).unwrap();
    // The map-1 and map-2 features pair; the empty map contributes nothing.
    assert_eq!(result.len(), 1);
    assert_eq!(result.features[0].len(), 2);
}

#[test]
fn zero_tolerances_are_refused_not_undefined() {
    // An m/z tolerance of 0 makes the source's m/z scale infinite and every
    // scaled distance NaN; this port refuses.
    let mut param = FeatureGroupingAlgorithmWNet::defaults().unwrap();
    set(&mut param, "max_mz_shift_ppm", ParamValue::Float(0.0));
    let mut algorithm = FeatureGroupingAlgorithmWNet::from_param(&param).unwrap();
    let mut result = ConsensusMap::new();
    assert!(algorithm.group(&class_test_maps(), &mut result).is_err());
    assert!(result.is_empty());
    // max_rt_shift 0 with the derived trash cost 0 divides by zero.
    let mut param = FeatureGroupingAlgorithmWNet::defaults().unwrap();
    set(&mut param, "max_rt_shift", ParamValue::Float(0.0));
    let mut algorithm = FeatureGroupingAlgorithmWNet::from_param(&param).unwrap();
    assert!(algorithm.group(&class_test_maps(), &mut result).is_err());
}

// ---- Parameters -----------------------------------------------------------

#[test]
fn defaults_tree() {
    let p = FeatureGroupingAlgorithmWNet::defaults().unwrap();
    let get = |key: &str| p.value(key).unwrap().clone();
    assert_eq!(get("distance_metric"), ParamValue::String("LINF".into()));
    assert_eq!(get("max_rt_shift"), ParamValue::Float(100.0));
    assert_eq!(get("mz_unit"), ParamValue::String("ppm".into()));
    assert_eq!(get("max_mz_shift_da"), ParamValue::Float(0.3));
    assert_eq!(get("max_mz_shift_ppm"), ParamValue::Float(10.0));
    assert_eq!(get("trash_cost"), ParamValue::Float(0.0));
    assert_eq!(
        get("normalize_intensities"),
        ParamValue::String("true".into())
    );
    assert_eq!(get("decharge_mz"), ParamValue::String("false".into()));
    // advanced tags: distance_metric and trash_cost only.
    for key in [
        "distance_metric",
        "max_rt_shift",
        "mz_unit",
        "max_mz_shift_da",
        "max_mz_shift_ppm",
        "trash_cost",
        "normalize_intensities",
        "decharge_mz",
    ] {
        let advanced = p.tags(key).unwrap().iter().any(|tag| tag == "advanced");
        assert_eq!(
            advanced,
            key == "distance_metric" || key == "trash_cost",
            "{key}"
        );
    }
    // Restrictions: the metric and the unit are enumerations, the shifts are
    // non-negative, and unknown values are refused.
    let mut bad = p.clone();
    let refused = bad
        .set_value("distance_metric", ParamValue::String("L3".into()), "", &[])
        .is_err()
        || FeatureGroupingAlgorithmWNet::from_param(&bad).is_err();
    assert!(refused);
    let mut bad = p.clone();
    let refused = bad
        .set_value("max_rt_shift", ParamValue::Float(-1.0), "", &[])
        .is_err()
        || FeatureGroupingAlgorithmWNet::from_param(&bad).is_err();
    assert!(refused);
}

#[test]
fn parameters_round_trip() {
    let algorithm = FeatureGroupingAlgorithmWNet {
        distance_metric: DistanceMetric::L2,
        max_rt_shift: 42.0,
        mz_ppm: false,
        max_mz_shift_da: 0.05,
        max_mz_shift_ppm: 3.0,
        trash_cost: 12.5,
        normalize_intensities: false,
        decharge_mz: true,
    };
    let param = algorithm.to_param().unwrap();
    assert_eq!(
        FeatureGroupingAlgorithmWNet::from_param(&param).unwrap(),
        algorithm
    );
    // Unknown keys are reported, not fatal.
    let mut with_unknown = param.clone();
    with_unknown
        .set_value("not_a_parameter", ParamValue::Integer(1), "", &[])
        .unwrap();
    let (parsed, warnings) =
        FeatureGroupingAlgorithmWNet::from_param_with_warnings(&with_unknown).unwrap();
    assert_eq!(parsed, algorithm);
    assert_eq!(warnings.len(), 1);
}
