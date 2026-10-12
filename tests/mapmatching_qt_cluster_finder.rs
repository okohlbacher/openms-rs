// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `QTClusterFinder_test.cpp` and `FeatureGroupingAlgorithmQT_test.cpp`
//! at source revision bc9cc12514c768385ce121d6ca4bb710fe1983c4
//! (`src/tests/class_tests/openms/source/`).
//!
//! Evidence tier 3: every expected value is a transcribed class-test literal,
//! except sections marked "tier 4". The grouping algorithm's `group` sections
//! are `NOT_TESTABLE` upstream; its tier-1 differential against the retained
//! TOPP outputs is `tests/mapmatching_qt_cluster_finder_topp.rs`. No expected
//! value was derived from this port's output.

use openms::analysis::mapmatching::feature_grouping_algorithm::FeatureGroupingAlgorithm;
use openms::analysis::mapmatching::feature_grouping_algorithm_qt::FeatureGroupingAlgorithmQT;
use openms::analysis::mapmatching::qt_cluster_finder::QTClusterFinder;
use openms::chemistry::AASequence;
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::kernel::features::{ConsensusMap, Feature, FeatureHandle, FeatureMap};
use openms::param::{Param, ParamValue};

fn set(param: &mut Param, key: &str, value: ParamValue) {
    let description = param.description(key).unwrap().to_string();
    let tags: Vec<String> = param.tags(key).unwrap().iter().cloned().collect();
    param.set_value(key, value, &description, &tags).unwrap();
}

fn feature(rt: f64, mz: f64, uid: u64, sequence: Option<&str>) -> Feature {
    let mut f = Feature::new(rt, mz, 0.0);
    f.unique_id = uid;
    if let Some(sequence) = sequence {
        let mut identification = PeptideIdentification::default();
        identification.hits.push(PeptideHit {
            sequence: AASequence::parse(sequence).unwrap(),
            ..PeptideHit::default()
        });
        f.peptide_identifications.push(identification);
    }
    f
}

fn keys(result: &ConsensusMap, index: usize) -> Vec<(u64, u64)> {
    result.features[index]
        .handles()
        .iter()
        .map(FeatureHandle::key)
        .collect()
}

#[test]
fn qt_cluster_finder_run_feature_maps() {
    let feat1 = feature(0.0, 0.0, 0, Some("AAA"));
    let feat2 = feature(100.0, 200.0, 1, Some("CCC"));
    let feat3 = feature(4.0, 0.04, 0, Some("DDD"));
    let feat4 = feature(5.0, 0.05, 1, Some("AAA"));
    let feat5 = feature(104.0, 200.04, 2, None);
    let mut input = vec![
        FeatureMap::from_features(vec![feat1.clone(), feat2.clone()]),
        FeatureMap::from_features(vec![feat3.clone(), feat4.clone(), feat5.clone()]),
    ];
    let (ind1, ind2) = ((0, 0), (0, 1));
    let (ind3, ind4, ind5) = ((1, 0), (1, 1), (1, 2));

    let mut param = QTClusterFinder::defaults().unwrap();
    set(
        &mut param,
        "distance_RT:max_difference",
        ParamValue::Float(5.1),
    );
    set(
        &mut param,
        "distance_MZ:max_difference",
        ParamValue::Float(0.1),
    );
    set(&mut param, "nr_partitions", ParamValue::Integer(1));
    let mut finder = QTClusterFinder::from_param(&param).unwrap();
    let mut result = ConsensusMap::new();
    finder.run_feature_maps(&input, &mut result).unwrap();
    assert_eq!(result.len(), 3);
    assert_eq!(keys(&result, 0), vec![ind2, ind5]);
    assert_eq!(keys(&result, 1), vec![ind1, ind3]);
    assert_eq!(keys(&result, 2), vec![ind4]);
    // The handles carry the features' positions.
    let h = result.features[0].handles()[0];
    assert_eq!((h.rt, h.mz), (feat2.rt, feat2.mz));

    // Annotation-specific matching (simple case):
    set(
        &mut param,
        "use_identifications",
        ParamValue::String("true".into()),
    );
    let mut finder = QTClusterFinder::from_param(&param).unwrap();
    finder.run_feature_maps(&input, &mut result).unwrap();
    assert_eq!(result.len(), 3);
    assert_eq!(keys(&result, 0), vec![ind2, ind5]);
    assert_eq!(keys(&result, 1), vec![ind1, ind4]);
    assert_eq!(keys(&result, 2), vec![ind3]);

    // Annotation-specific matching (complex case):
    let feat6 = feature(104.0, 200.04, 0, Some("EEE"));
    let feat7 = feature(108.0, 200.08, 1, Some("CCC"));
    input.push(FeatureMap::from_features(vec![feat6, feat7]));
    finder.run_feature_maps(&input, &mut result).unwrap();
    assert_eq!(result.len(), 4);
    // "ind6" is closer, but its annotation doesn't match.
    assert_eq!(keys(&result, 0), vec![ind2, ind5, (2, 1)]);
}

#[test]
fn qt_cluster_finder_only_empty_maps() {
    let input = vec![FeatureMap::new(), FeatureMap::new()];
    let mut finder = QTClusterFinder::default();
    let mut result = ConsensusMap::new();
    finder.run_feature_maps(&input, &mut result).unwrap();
    assert_eq!(result.len(), 0);
    let mut p = finder.to_param().unwrap();
    set(&mut p, "nr_partitions", ParamValue::Integer(1));
    let mut finder = QTClusterFinder::from_param(&p).unwrap();
    let mut result_single = ConsensusMap::new();
    finder.run_feature_maps(&input, &mut result_single).unwrap();
    assert_eq!(result_single.len(), 0);
}

/// Tier 4: fewer than two maps is an error (source
/// `Exception::IllegalArgument`), the defaults are the source's, and the
/// grouping algorithm post-processes (sorted by size, then maps).
#[test]
fn qt_defaults_and_grouping() {
    let p = QTClusterFinder::defaults().unwrap();
    assert_eq!(p.value("nr_partitions").unwrap(), &ParamValue::Integer(100));
    assert_eq!(
        p.value("min_nr_diffs_per_bin").unwrap(),
        &ParamValue::Integer(50)
    );
    assert_eq!(
        p.value("min_IDscore_forTolCalc").unwrap(),
        &ParamValue::Float(1.0)
    );
    assert_eq!(p.value("noID_penalty").unwrap(), &ParamValue::Float(0.0));
    assert_eq!(
        p.value("distance_MZ:max_difference").unwrap(),
        &ParamValue::Float(0.3)
    );

    let one = vec![FeatureMap::from_features(vec![feature(
        1.0, 100.0, 1, None,
    )])];
    let mut finder = QTClusterFinder::default();
    let mut result = ConsensusMap::new();
    assert!(finder.run_feature_maps(&one, &mut result).is_err());

    let maps = vec![
        FeatureMap::from_features(vec![
            feature(10.0, 100.0, 1, None),
            feature(50.0, 300.0, 2, None),
        ]),
        FeatureMap::from_features(vec![feature(11.0, 100.01, 3, None)]),
    ];
    let mut algorithm = FeatureGroupingAlgorithmQT::default();
    let mut out = ConsensusMap::new();
    algorithm.group(&maps, &mut out).unwrap();
    assert_eq!(out.len(), 2);
    assert_eq!(keys(&out, 0), vec![(0, 1), (1, 3)]);
    assert_eq!(keys(&out, 1), vec![(0, 2)]);
    assert!(algorithm.group(&maps[..1], &mut out).is_err());
}
