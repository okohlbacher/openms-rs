// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Evidence tier 1 for `FeatureGroupingAlgorithmKD` (and, in cases 4 to 8,
//! `MapAlignmentAlgorithmKD`): the algorithm path of the upstream
//! `TOPP_FeatureLinkerUnlabeledKD_1` to `_8` tests, compared with the
//! retained C++ outputs.
//!
//! Upstream definition: `topp/CMakeLists.txt:2927-2950` at test-data revision
//! 0cb15f2 runs `FeatureLinkerUnlabeledKD -test` with the case's `.ini` (cases
//! 4 to 7 share `_4_parameters.ini` with command-line merging overrides; case
//! 8 runs on defaults with `-algorithm:link:charge_merging Any`) and compares
//! with `FeatureLinkerUnlabeledKD_<n>_output.consensusXML` through `${DIFF}
//! -whitelist "id=" "href="`. The tool (`FeatureLinkerBase.cpp` at topp
//! revision 174b576) loads each featureXML without subordinates and convex
//! hulls, clears the meta values except `dc_charge_adducts` and `Group`, and
//! calls `group`; case 2 links two copies of the retained `_1` output as
//! consensus maps. It then assigns fresh unique IDs and sorts
//! identifications, which change neither the order nor the values compared.
//! The parameters are the tool defaults updated with the `.ini` section
//! (`Param::update`, which drops the `None` merging values of
//! `_4_parameters.ini` as invalid), then the command-line values.
//!
//! Compared: order, handle identities exactly, values within the six
//! significant digits stored (relative 5e-6, plus one `f32` unit for `float`
//! values, plus an absolute 1e-6), charge, and each consensus feature's
//! peptide identifications (map index and first-hit sequence, in the order
//! the tool sorts them). The executed Release build agrees with every
//! retained value within this tolerance (oracle `M2-labeled`, compared at
//! full precision in `tests/mapmatching_kd_grouping_oracle.rs`).

#![cfg(all(feature = "featurexml", feature = "consensusxml", feature = "paramxml"))]

use openms::analysis::mapmatching::feature_grouping_algorithm::{
    FeatureGroupingAlgorithm, transfer_subelements,
};
use openms::analysis::mapmatching::feature_grouping_algorithm_kd::FeatureGroupingAlgorithmKD;
use openms::format::{consensusxml, featurexml, paramxml};
use openms::kernel::features::{ConsensusFeature, ConsensusMap, FeatureMap};
use openms::param::{Param, ParamValue};

const DIR: &str = "tests/data";

fn close(actual: f64, expected: f64, relative: f64, what: &str) {
    let bound = relative * expected.abs() + 1e-6;
    assert!(
        (actual - expected).abs() <= bound,
        "{what}: produced {actual}, retained {expected}"
    );
}

/// The peptide identifications of a consensus feature as `(map_index, first
/// hit sequence)`, ordered as the tool's `sortPeptideIdentificationsByMapIndex`
/// leaves them: stably by `map_index`, identifications without one last.
fn identifications(feature: &ConsensusFeature) -> Vec<(Option<i64>, String)> {
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

fn set(param: &mut Param, key: &str, value: ParamValue) {
    let description = param.description(key).unwrap().to_string();
    let tags: Vec<String> = param.tags(key).unwrap().iter().cloned().collect();
    param.set_value(key, value, &description, &tags).unwrap();
}

fn load(file: &str) -> FeatureMap {
    let mut map = featurexml::load(format!("{DIR}/{file}")).unwrap();
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

fn algorithm(ini: Option<&str>, overrides: &[(&str, &str)]) -> FeatureGroupingAlgorithmKD {
    let mut p = FeatureGroupingAlgorithmKD::defaults().unwrap();
    if let Some(ini) = ini {
        let full = paramxml::load(format!("{DIR}/{ini}")).unwrap();
        p.update(
            &full
                .copy("FeatureLinkerUnlabeledKD:1:algorithm:", true)
                .unwrap(),
            false,
        )
        .unwrap();
    }
    for (key, value) in overrides {
        set(&mut p, key, ParamValue::String((*value).into()));
    }
    FeatureGroupingAlgorithmKD::from_param(&p).unwrap()
}

fn compare(out: &ConsensusMap, case: usize) {
    let expected = consensusxml::load(format!(
        "{DIR}/feature_linker_unlabeled_kd_{case}_output.consensusXML"
    ))
    .unwrap();
    assert_eq!(
        out.len(),
        expected.len(),
        "KD_{case}: consensus feature count"
    );
    let float_rel = 5.0e-6 + f64::from(f32::EPSILON);
    for (index, (produced, retained)) in out.features.iter().zip(&expected.features).enumerate() {
        let at = format!("KD_{case}, consensus feature {index}");
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

fn features(case: usize, mut algorithm: FeatureGroupingAlgorithmKD, inputs: &[&str]) {
    let maps: Vec<FeatureMap> = inputs.iter().map(|f| load(f)).collect();
    let mut out = ConsensusMap::new();
    algorithm.group(&maps, &mut out).unwrap();
    compare(&out, case);
}

const DC: [&str; 5] = [
    "feature_linker_unlabeled_kd_dc_input1.featureXML",
    "feature_linker_unlabeled_kd_dc_input2.featureXML",
    "feature_linker_unlabeled_kd_dc_input3.featureXML",
    "feature_linker_unlabeled_kd_dc_input1.featureXML",
    "feature_linker_unlabeled_kd_dc_input2.featureXML",
];

#[test]
fn topp_feature_linker_unlabeled_kd_1() {
    features(
        1,
        algorithm(Some("feature_linker_unlabeled_kd_1_parameters.ini"), &[]),
        &[
            "feature_linker_unlabeled_1_input1.featureXML",
            "feature_linker_unlabeled_1_input2.featureXML",
            "feature_linker_unlabeled_1_input3.featureXML",
        ],
    );
}

/// Case 2 sets `keep_subelements`: the tool copies every peptide
/// identification's `map_index` to `old_map_index` (assigned and unassigned,
/// `applyFunctionOnPeptideIDs(..., true)`), groups, and then replaces the
/// handles by the grouped consensus features' own handles with
/// `transferSubelements`.
#[test]
fn topp_feature_linker_unlabeled_kd_2() {
    let mut a = algorithm(Some("feature_linker_unlabeled_kd_2_parameters.ini"), &[]);
    let mut input = consensusxml::load(format!(
        "{DIR}/feature_linker_unlabeled_kd_1_output.consensusXML"
    ))
    .unwrap();
    let save_old_map_index = |id: &mut openms::identification::PeptideIdentification| {
        if let Some(value) = id.metadata.get("map_index").cloned() {
            id.metadata.insert("old_map_index".into(), value);
        }
    };
    for feature in &mut input.features {
        feature
            .peptide_identifications
            .iter_mut()
            .for_each(save_old_map_index);
    }
    input
        .unassigned_peptide_identifications
        .iter_mut()
        .for_each(save_old_map_index);
    let maps = vec![input.clone(), input];
    let mut out = ConsensusMap::new();
    a.group_consensus(&maps, &mut out).unwrap();
    transfer_subelements(&maps, &mut out).unwrap();
    compare(&out, 2);
}

#[test]
fn topp_feature_linker_unlabeled_kd_3() {
    features(
        3,
        algorithm(Some("feature_linker_unlabeled_kd_3_parameters.ini"), &[]),
        &[
            "feature_linker_unlabeled_qt_3_input1.featureXML",
            "feature_linker_unlabeled_qt_3_input2.featureXML",
        ],
    );
}

fn dc(case: usize, charge: &str, adduct: &str) {
    features(
        case,
        algorithm(
            Some("feature_linker_unlabeled_kd_4_parameters.ini"),
            &[
                ("link:charge_merging", charge),
                ("link:adduct_merging", adduct),
            ],
        ),
        &DC,
    );
}

#[test]
fn topp_feature_linker_unlabeled_kd_4() {
    dc(4, "Identical", "Any");
}

#[test]
fn topp_feature_linker_unlabeled_kd_5() {
    dc(5, "With_charge_zero", "Any");
}

#[test]
fn topp_feature_linker_unlabeled_kd_6() {
    dc(6, "Any", "With_unknown_adducts");
}

#[test]
fn topp_feature_linker_unlabeled_kd_7() {
    dc(7, "Any", "Identical");
}

#[test]
fn topp_feature_linker_unlabeled_kd_8() {
    features(
        8,
        algorithm(None, &[("link:charge_merging", "Any")]),
        &[
            "feature_linker_unlabeled_kd_charge_any_input1.featureXML",
            "feature_linker_unlabeled_kd_charge_any_input2.featureXML",
        ],
    );
}
