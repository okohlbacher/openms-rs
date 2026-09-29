// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Evidence tier 1 for `FeatureGroupingAlgorithmQT` and `QTClusterFinder`:
//! the algorithm path of the upstream `TOPP_FeatureLinkerUnlabeledQT_1`,
//! `_3`, `_4`, `_5` and `_6` tests, compared with the retained C++ outputs.
//!
//! Upstream definition: `topp/CMakeLists.txt:2907-2925` at test-data revision
//! 0cb15f2 runs `FeatureLinkerUnlabeledQT -test` with the case's parameters
//! (case 4 adds `-algorithm:use_identifications`, case 6 runs on defaults
//! with `-algorithm:use_identifications -algorithm:distance_RT:max_difference
//! 200`) and compares with `FeatureLinkerUnlabeledQT_<n>_output.consensusXML`
//! through `${DIFF} -whitelist "id=" "href="`. Case 2 links consensusXML input
//! and is not reproduced here. The tool (`FeatureLinkerBase.cpp` at topp
//! revision 174b576) loads every featureXML without convex hulls and
//! subordinates, clears each feature's meta values except
//! `dc_charge_adducts` and `Group`, calls `group(maps, out_map)`, then assigns
//! fresh unique IDs and sorts identifications, which change neither the order
//! nor the values compared here.
//!
//! Compared as in `tests/mapmatching_feature_grouping_topp.rs`: order, handle
//! identities exactly, values within the six significant digits stored.

#![cfg(all(feature = "featurexml", feature = "consensusxml", feature = "paramxml"))]

use openms::analysis::mapmatching::feature_grouping_algorithm::FeatureGroupingAlgorithm;
use openms::analysis::mapmatching::feature_grouping_algorithm_qt::FeatureGroupingAlgorithmQT;
use openms::format::{consensusxml, featurexml, paramxml};
use openms::kernel::features::{ConsensusMap, FeatureMap};
use openms::param::{Param, ParamValue};

const DIR: &str = "tests/data";

fn close(actual: f64, expected: f64, relative: f64, what: &str) {
    let bound = relative * expected.abs() + 1e-6;
    assert!(
        (actual - expected).abs() <= bound,
        "{what}: produced {actual}, retained {expected}"
    );
}

fn set(param: &mut Param, key: &str, value: ParamValue) {
    let description = param.description(key).unwrap().to_string();
    let tags: Vec<String> = param.tags(key).unwrap().iter().cloned().collect();
    param.set_value(key, value, &description, &tags).unwrap();
}

/// Load an input as the tool does.
fn load(path: &str) -> FeatureMap {
    let mut map = featurexml::load(path).unwrap();
    for feature in &mut map.features {
        feature.subordinates.clear();
        feature.convex_hulls.clear();
        let adduct = feature.metadata.get("dc_charge_adducts").cloned();
        let group = feature.metadata.get("Group").cloned();
        feature.metadata.clear();
        if let Some(adduct) = adduct {
            feature.metadata.insert("dc_charge_adducts".into(), adduct);
        }
        if let Some(group) = group {
            feature.metadata.insert("Group".into(), group);
        }
    }
    map
}

fn ini_algorithm(ini: &str) -> Param {
    paramxml::load(format!("{DIR}/{ini}"))
        .unwrap()
        .copy("FeatureLinkerUnlabeledQT:1:algorithm:", true)
        .unwrap()
}

fn compare(out: &ConsensusMap, expected_file: &str, case: &str) {
    let expected = consensusxml::load(format!("{DIR}/{expected_file}")).unwrap();
    assert_eq!(out.len(), expected.len(), "{case}: consensus feature count");
    let float_rel = 5.0e-6 + f64::from(f32::EPSILON);
    for (index, (produced, retained)) in out.features.iter().zip(&expected.features).enumerate() {
        let at = format!("{case}, consensus feature {index}");
        let keys: Vec<_> = produced.handles().iter().map(|h| h.key()).collect();
        let retained_keys: Vec<_> = retained.handles().iter().map(|h| h.key()).collect();
        assert_eq!(keys, retained_keys, "{at}: handles");
        for (a, b) in produced.handles().iter().zip(retained.handles()) {
            close(a.rt, b.rt, 5.0e-6, &format!("{at} handle RT"));
            close(a.mz, b.mz, 5.0e-6, &format!("{at} handle m/z"));
            close(
                f64::from(a.intensity),
                f64::from(b.intensity),
                float_rel,
                &format!("{at} handle intensity"),
            );
        }
        close(produced.rt, retained.rt, 5.0e-6, &format!("{at} RT"));
        close(produced.mz, retained.mz, 5.0e-6, &format!("{at} m/z"));
        close(
            f64::from(produced.intensity),
            f64::from(retained.intensity),
            float_rel,
            &format!("{at} intensity"),
        );
        close(
            f64::from(produced.quality),
            f64::from(retained.quality),
            float_rel,
            &format!("{at} quality"),
        );
        assert_eq!(produced.charge, retained.charge, "{at}: charge");
    }
}

fn run(param: &Param, inputs: &[&str], expected: &str, case: &str) {
    let (mut algorithm, _warnings) =
        FeatureGroupingAlgorithmQT::from_param_with_warnings(param).unwrap();
    let maps: Vec<FeatureMap> = inputs.iter().map(|f| load(&format!("{DIR}/{f}"))).collect();
    let mut out = ConsensusMap::new();
    algorithm.group(&maps, &mut out).unwrap();
    compare(&out, expected, case);
}

const UNLABELED_1: [&str; 3] = [
    "feature_linker_unlabeled_1_input1.featureXML",
    "feature_linker_unlabeled_1_input2.featureXML",
    "feature_linker_unlabeled_1_input3.featureXML",
];
const QT_5: [&str; 3] = [
    "feature_linker_unlabeled_qt_5_input1.featureXML",
    "feature_linker_unlabeled_qt_5_input2.featureXML",
    "feature_linker_unlabeled_qt_5_input3.featureXML",
];

#[test]
fn topp_feature_linker_unlabeled_qt_1() {
    let param = ini_algorithm("feature_linker_unlabeled_qt_1_parameters.ini");
    run(
        &param,
        &UNLABELED_1,
        "feature_linker_unlabeled_qt_1_output.consensusXML",
        "QT_1",
    );
}

#[test]
fn topp_feature_linker_unlabeled_qt_3() {
    let param = ini_algorithm("feature_linker_unlabeled_qt_3_parameters.ini");
    run(
        &param,
        &[
            "feature_linker_unlabeled_qt_3_input1.featureXML",
            "feature_linker_unlabeled_qt_3_input2.featureXML",
        ],
        "feature_linker_unlabeled_qt_3_output.consensusXML",
        "QT_3",
    );
}

#[test]
fn topp_feature_linker_unlabeled_qt_4() {
    let mut param = ini_algorithm("feature_linker_unlabeled_qt_1_parameters.ini");
    set(
        &mut param,
        "use_identifications",
        ParamValue::String("true".into()),
    );
    run(
        &param,
        &UNLABELED_1,
        "feature_linker_unlabeled_qt_4_output.consensusXML",
        "QT_4",
    );
}

#[test]
fn topp_feature_linker_unlabeled_qt_5() {
    let param = ini_algorithm("feature_linker_unlabeled_qt_1_parameters.ini");
    run(
        &param,
        &QT_5,
        "feature_linker_unlabeled_qt_5_output.consensusXML",
        "QT_5",
    );
}

#[test]
fn topp_feature_linker_unlabeled_qt_6() {
    let mut param = FeatureGroupingAlgorithmQT::defaults().unwrap();
    set(
        &mut param,
        "use_identifications",
        ParamValue::String("true".into()),
    );
    set(
        &mut param,
        "distance_RT:max_difference",
        ParamValue::Float(200.0),
    );
    run(
        &param,
        &QT_5,
        "feature_linker_unlabeled_qt_6_output.consensusXML",
        "QT_6",
    );
}
