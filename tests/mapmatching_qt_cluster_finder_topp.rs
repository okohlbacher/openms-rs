// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Evidence tier 1 for `FeatureGroupingAlgorithmQT` and `QTClusterFinder`:
//! the algorithm path of the upstream `TOPP_FeatureLinkerUnlabeledQT_1`,
//! `_3`, `_4` and `_6` tests, compared with the retained C++ outputs. `_5` is
//! compared with the executed Release build in
//! `mapmatching_qt_tie_order_oracle.rs`; its retained output differs in one tie
//! (see below).
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
//! identities exactly, values within the six significant digits stored, and
//! each consensus feature's peptide identifications (map index and first-hit
//! sequence, in the order the tool sorts them).

#![cfg(all(feature = "featurexml", feature = "consensusxml", feature = "paramxml"))]

use openms::analysis::mapmatching::feature_grouping_algorithm::FeatureGroupingAlgorithm;
use openms::analysis::mapmatching::feature_grouping_algorithm_qt::FeatureGroupingAlgorithmQT;
use openms::format::{consensusxml, featurexml, paramxml};
use openms::kernel::features::{ConsensusMap, FeatureMap};
use openms::param::{Param, ParamValue};

const DIR: &str = "tests/data";
/// The source's null date attribute, which this port's readers refuse.
const NULL_DATE: &str = "date=\"0000-00-00T00:00:00\"";
/// The stand-in used only when a reader refuses [`NULL_DATE`].
const VALID_DATE: &str = "date=\"1970-01-01T00:00:00\"";

/// The peptide identifications of a consensus feature as `(map_index, first
/// hit sequence)`, ordered as the tool's `sortPeptideIdentificationsByMapIndex`
/// leaves them: stably by `map_index`, identifications without one last.
fn identifications(
    feature: &openms::kernel::features::ConsensusFeature,
) -> Vec<(Option<i64>, String)> {
    let mut ids: Vec<(Option<i64>, String)> = feature
        .peptide_identifications
        .iter()
        .map(|id| {
            (
                id.metadata.get("map_index").and_then(|v| v.as_i64().ok()),
                id.hits
                    .first()
                    .map(|h| h.sequence.to_string())
                    .unwrap_or_default(),
            )
        })
        .collect();
    ids.sort_by_key(|(index, _)| (index.is_none(), *index));
    ids
}

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
///
/// The `QT_5` inputs carry `IdentificationRun date="0000-00-00T00:00:00"`.
/// The source's `DateTime::fromString` turns that into a null date
/// (`DATASTRUCTURES/DateTime.cpp:565` at bc9cc12: invalid fields return an
/// invalid `DateTime`, no exception), but this port's featureXML reader
/// refuses it ("completion timestamp is outside calendar bounds", a reader
/// defect reported by this lane). Only in exactly that case is the attribute
/// replaced by a valid date before parsing; the date reaches the protein
/// identification's metadata and never the grouping. Once the reader accepts
/// the null date, the direct load succeeds and the replacement goes unused.
fn load(path: &str) -> FeatureMap {
    let mut map = match featurexml::load(path) {
        Ok(map) => map,
        Err(error) => {
            assert!(
                error
                    .to_string()
                    .contains("completion timestamp is outside calendar bounds"),
                "{path}: {error}"
            );
            let text = std::fs::read_to_string(path).unwrap();
            assert!(text.contains(NULL_DATE), "{path}: {error}");
            featurexml::read(text.replace(NULL_DATE, VALID_DATE).as_bytes()).unwrap()
        }
    };
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
    // The retained `QT_5`/`QT_6` outputs carry the same null date; see `load`.
    let path = format!("{DIR}/{expected_file}");
    let expected = match consensusxml::load(&path) {
        Ok(map) => map,
        Err(error) => {
            assert!(
                error
                    .to_string()
                    .contains("completion timestamp is outside calendar bounds"),
                "{path}: {error}"
            );
            let text = std::fs::read_to_string(&path).unwrap();
            assert!(text.contains(NULL_DATE), "{path}: {error}");
            consensusxml::read(text.replace(NULL_DATE, VALID_DATE).as_bytes()).unwrap()
        }
    };
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
        assert_eq!(
            identifications(produced),
            identifications(retained),
            "{at}: peptide identifications"
        );
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

// `TOPP_FeatureLinkerUnlabeledQT_5`: three clusters of exactly equal quality
// (0.974444) and identical handle keys (every input feature has unique ID 0).
// The extraction order of the three is the source's `boost::heap::
// fibonacci_heap` order over Boost's hash-ordered grid, which this port now
// reproduces; against the *executed* Linux x86_64 Release build it is exact
// (`mapmatching_qt_tie_order_oracle.rs`, case `topp_qt5`).
//
// The retained `FeatureLinkerUnlabeledQT_5_output.consensusXML` is NOT what
// that build produces. In the middle cluster (RT 300 / 31x / 320) the build
// takes the map-1 feature at RT 310.1 where the retained file has the one at
// RT 310.0: two equally near candidates, decided by container order or by
// another difference of the build that generated the file (what differs there
// is not established; the retained outputs of this family were not produced by
// the pinned Linux build). The test below
// states exactly that difference, with the retained values from the
// C++-generated file and the executed values from the oracle (the cluster's
// centroid RT 310.0333.. = (300 + 310.1 + 320) / 3 in `qt_oracle.tsv`): the
// two outer clusters agree in every handle, the middle cluster differs in the
// map-1 handle's RT alone.
#[test]
fn topp_feature_linker_unlabeled_qt_5_agrees_with_retained_except_one_tie() {
    let mut param = ini_algorithm("feature_linker_unlabeled_qt_1_parameters.ini");
    set(
        &mut param,
        "distance_RT:max_difference",
        ParamValue::Float(200.0),
    );
    let (mut algorithm, _warnings) =
        FeatureGroupingAlgorithmQT::from_param_with_warnings(&param).unwrap();
    let maps: Vec<FeatureMap> = QT_5.iter().map(|f| load(&format!("{DIR}/{f}"))).collect();
    let mut out = ConsensusMap::new();
    algorithm.group(&maps, &mut out).unwrap();

    let path = format!("{DIR}/feature_linker_unlabeled_qt_5_output.consensusXML");
    let text = std::fs::read_to_string(&path).unwrap();
    let retained = consensusxml::read(text.replace(NULL_DATE, VALID_DATE).as_bytes()).unwrap();
    assert_eq!(out.len(), 3);
    assert_eq!(retained.len(), 3);
    for (index, (produced, expected)) in out.features.iter().zip(&retained.features).enumerate() {
        let produced_rts: Vec<f64> = produced.handles().iter().map(|h| h.rt).collect();
        let expected_rts: Vec<f64> = expected.handles().iter().map(|h| h.rt).collect();
        if index == 1 {
            assert_eq!(expected_rts, vec![300.0, 310.0, 320.0], "retained");
            assert_eq!(produced_rts, vec![300.0, 310.1, 320.0], "executed build");
        } else {
            assert_eq!(produced_rts, expected_rts, "consensus feature {index}");
        }
        assert!(
            (f64::from(produced.quality) - f64::from(expected.quality)).abs() < 5e-6,
            "quality {index}"
        );
    }
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
