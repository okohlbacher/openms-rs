// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Evidence tier 1 for `FeatureGroupingAlgorithmWNet`: a differential against
//! executed C++.
//!
//! Oracle `M2-wnet` (`oracle/M2-wnet/driver.cpp`; host, build and library
//! versions as in `tests/mapmatching_wnet_oracle.rs`) ran the pinned
//! `FeatureGroupingAlgorithmWNet::group` twice with identical output.
//! `wnet_oracle_groups.tsv` retains, per case, `size` (consensus features,
//! protein and unassigned peptide identifications) and per consensus feature
//! `cf RT m/z intensity quality charge size handles...` and `pep` (count and
//! `sequence@map_index` of each identification), printed with `%.17g`:
//!
//! - `g_*`: inputs loaded as `FeatureLinkerBase.cpp` loads them (no
//!   subordinates or convex hulls; meta values cleared except
//!   `dc_charge_adducts` and `Group`), then `group` with the case's
//!   parameter overrides. `g_wnet1` is `TOPP_FeatureLinkerWNet_1`'s algorithm
//!   path (the FeatureLinkerUnlabeled_1 inputs, default parameters);
//!   `g_cons` links the three retained FeatureLinkerUnlabeled outputs as
//!   consensus maps; `g_dc*` use the adduct/charge test inputs of the KD
//!   tests, with and without `decharge_mz`.
//! - `class_*`: the sections of `FeatureGroupingAlgorithmWNet_test.cpp`.
//!
//! Compared: the counts, every handle's map index and unique ID exactly (in
//! order), RT and m/z to a relative 1e-12, intensity and quality (`float` in
//! the source) to one `f32` unit, the charge, and the identifications.

#![cfg(all(feature = "featurexml", feature = "consensusxml", feature = "paramxml"))]

use std::collections::BTreeMap;

use openms::analysis::mapmatching::feature_grouping_algorithm::FeatureGroupingAlgorithm;
use openms::analysis::mapmatching::feature_grouping_algorithm_wnet::FeatureGroupingAlgorithmWNet;
use openms::format::{consensusxml, featurexml, paramxml};
use openms::kernel::features::{ConsensusFeature, ConsensusMap, Feature, FeatureMap};
use openms::param::{Param, ParamValue};

const DIR: &str = "tests/data";

fn set(param: &mut Param, key: &str, value: ParamValue) {
    let description = param.description(key).unwrap().to_string();
    let tags: Vec<String> = param.tags(key).unwrap().iter().cloned().collect();
    param.set_value(key, value, &description, &tags).unwrap();
}

/// Load an input as `FeatureLinkerBase.cpp` does.
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

fn rel(actual: f64, expected: f64, relative: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "{what}: produced {actual:e}, oracle {expected:e}"
    );
}

struct Case {
    size: Vec<usize>,
    features: Vec<Vec<String>>,
    peptides: Vec<Vec<String>>,
}

fn oracle() -> BTreeMap<String, Case> {
    let text = std::fs::read_to_string(format!("{DIR}/wnet_oracle_groups.tsv")).unwrap();
    let mut out: BTreeMap<String, Case> = BTreeMap::new();
    let mut name = String::new();
    for line in text.lines() {
        let f: Vec<String> = line.split('\t').map(str::to_string).collect();
        match f[0].as_str() {
            "case" => {
                name = f[1].clone();
                out.insert(
                    name.clone(),
                    Case {
                        size: Vec::new(),
                        features: Vec::new(),
                        peptides: Vec::new(),
                    },
                );
            }
            "size" => {
                out.get_mut(&name).unwrap().size =
                    f[1..].iter().map(|x| x.parse().unwrap()).collect();
            }
            "cf" => out.get_mut(&name).unwrap().features.push(f),
            "pep" => out.get_mut(&name).unwrap().peptides.push(f),
            _ => {}
        }
    }
    out
}

fn compare(case: &str, out: &ConsensusMap) {
    let all = oracle();
    let expected = &all[case];
    assert_eq!(
        vec![
            out.len(),
            out.protein_identifications.len(),
            out.unassigned_peptide_identifications.len()
        ],
        expected.size,
        "{case}: sizes"
    );
    assert_eq!(out.features.len(), expected.features.len(), "{case}: count");
    let float_rel = f64::from(f32::EPSILON);
    for (index, produced) in out.features.iter().enumerate() {
        let at = format!("{case}, consensus feature {index}");
        let f = &expected.features[index];
        let n: usize = f[6].parse().unwrap();
        let keys: Vec<String> = produced
            .handles()
            .iter()
            .map(|h| format!("{}:{}", h.map_index, h.unique_id))
            .collect();
        assert_eq!(keys, f[7..7 + n].to_vec(), "{at}: handles");
        rel(
            produced.rt,
            f[1].parse().unwrap(),
            1e-12,
            &format!("{at} RT"),
        );
        rel(
            produced.mz,
            f[2].parse().unwrap(),
            1e-12,
            &format!("{at} m/z"),
        );
        rel(
            f64::from(produced.intensity),
            f[3].parse().unwrap(),
            float_rel,
            &format!("{at} intensity"),
        );
        assert_eq!(
            produced.quality,
            f[4].parse::<f64>().unwrap() as f32,
            "{at}: quality"
        );
        assert_eq!(
            produced.charge,
            f[5].parse::<i32>().unwrap(),
            "{at}: charge"
        );
        let peptides: Vec<String> = produced
            .peptide_identifications
            .iter()
            .map(|id| {
                let sequence = id
                    .hits
                    .first()
                    .map_or_else(|| "-".to_string(), |h| h.sequence.to_string());
                let map_index = id
                    .metadata
                    .get("map_index")
                    .and_then(|v| v.as_i64().ok())
                    .map_or_else(|| "-".to_string(), |v| v.to_string());
                format!("{sequence}@{map_index}")
            })
            .collect();
        let line = &expected.peptides[index];
        assert_eq!(
            line[1].parse::<usize>().unwrap(),
            peptides.len(),
            "{at}: peptide count"
        );
        assert_eq!(peptides, line[2..].to_vec(), "{at}: peptides");
    }
}

fn algorithm(overrides: &[(&str, ParamValue)]) -> FeatureGroupingAlgorithmWNet {
    let mut p = FeatureGroupingAlgorithmWNet::defaults().unwrap();
    for (key, value) in overrides {
        set(&mut p, key, value.clone());
    }
    FeatureGroupingAlgorithmWNet::from_param(&p).unwrap()
}

fn string(value: &str) -> ParamValue {
    ParamValue::String(value.into())
}

fn run_features(case: &str, overrides: &[(&str, ParamValue)], inputs: &[&str]) {
    let mut a = algorithm(overrides);
    let maps: Vec<FeatureMap> = inputs.iter().map(|f| load(f)).collect();
    let mut out = ConsensusMap::new();
    a.group(&maps, &mut out).unwrap();
    compare(case, &out);
}

const UNLABELED_1: [&str; 3] = [
    "feature_linker_unlabeled_1_input1.featureXML",
    "feature_linker_unlabeled_1_input2.featureXML",
    "feature_linker_unlabeled_1_input3.featureXML",
];
const UNLABELED_2: [&str; 3] = [
    "feature_linker_unlabeled_2_input1.featureXML",
    "feature_linker_unlabeled_2_input2.featureXML",
    "feature_linker_unlabeled_2_input3.featureXML",
];
const UNLABELED_3: [&str; 2] = [
    "feature_linker_unlabeled_3_input1.featureXML",
    "feature_linker_unlabeled_3_input2.featureXML",
];
const DC: [&str; 3] = [
    "feature_linker_unlabeled_kd_dc_input1.featureXML",
    "feature_linker_unlabeled_kd_dc_input2.featureXML",
    "feature_linker_unlabeled_kd_dc_input3.featureXML",
];
const ANY: [&str; 2] = [
    "feature_linker_unlabeled_kd_charge_any_input1.featureXML",
    "feature_linker_unlabeled_kd_charge_any_input2.featureXML",
];

#[test]
fn oracle_defaults_tree() {
    // `FeatureGroupingAlgorithmWNet().getParameters()` stored with
    // `ParamXMLFile`: every value, description, tag and restriction.
    let expected = paramxml::load(format!("{DIR}/wnet_oracle_defaults.ini")).unwrap();
    let produced = FeatureGroupingAlgorithmWNet::defaults().unwrap();
    assert_eq!(produced, expected);
}

#[test]
fn oracle_g_wnet1() {
    run_features("g_wnet1", &[], &UNLABELED_1);
}

#[test]
fn oracle_g_da() {
    run_features("g_da", &[("mz_unit", string("Da"))], &UNLABELED_1);
}

#[test]
fn oracle_g_l2() {
    run_features(
        "g_l2",
        &[
            ("distance_metric", string("L2")),
            ("max_rt_shift", ParamValue::Float(50.0)),
            ("max_mz_shift_ppm", ParamValue::Float(20.0)),
            ("trash_cost", ParamValue::Float(30.0)),
        ],
        &UNLABELED_1,
    );
}

#[test]
fn oracle_g_l1_nonorm() {
    run_features(
        "g_l1_nonorm",
        &[
            ("distance_metric", string("L1")),
            ("normalize_intensities", string("false")),
            ("max_rt_shift", ParamValue::Float(200.0)),
        ],
        &UNLABELED_1,
    );
}

#[test]
fn oracle_g_trash() {
    run_features(
        "g_trash",
        &[("trash_cost", ParamValue::Float(500.0))],
        &UNLABELED_1,
    );
}

#[test]
fn oracle_g_tight() {
    run_features(
        "g_tight",
        &[
            ("max_rt_shift", ParamValue::Float(5.0)),
            ("max_mz_shift_ppm", ParamValue::Float(2.0)),
        ],
        &UNLABELED_1,
    );
}

#[test]
fn oracle_g_u2() {
    run_features("g_u2", &[], &UNLABELED_2);
}

#[test]
fn oracle_g_u3() {
    run_features("g_u3", &[], &UNLABELED_3);
}

#[test]
fn oracle_g_dc() {
    run_features("g_dc", &[("decharge_mz", string("true"))], &DC);
}

#[test]
fn oracle_g_dc_da() {
    run_features(
        "g_dc_da",
        &[("decharge_mz", string("true")), ("mz_unit", string("Da"))],
        &DC,
    );
}

#[test]
fn oracle_g_nodc() {
    run_features("g_nodc", &[], &DC);
}

#[test]
fn oracle_g_any() {
    run_features("g_any", &[], &ANY);
}

#[test]
fn oracle_g_cons() {
    let mut a = algorithm(&[]);
    let maps: Vec<ConsensusMap> = [
        "feature_linker_unlabeled_1_output.consensusXML",
        "feature_linker_unlabeled_2_output.consensusXML",
        "feature_linker_unlabeled_3_output.consensusXML",
    ]
    .iter()
    .map(|f| consensusxml::load(format!("{DIR}/{f}")).unwrap())
    .collect();
    let mut out = ConsensusMap::new();
    a.group_consensus(&maps, &mut out).unwrap();
    compare("g_cons", &out);
}

fn feature(unique_id: u64, mz: f64, rt: f64, intensity: f32, charge: i32) -> Feature {
    let mut f = Feature::new(rt, mz, intensity);
    f.unique_id = unique_id;
    f.charge = charge;
    f
}

fn class_maps() -> Vec<FeatureMap> {
    vec![
        FeatureMap::from_features(vec![
            feature(1, 500.0, 100.0, 1000.0, 0),
            feature(2, 600.0, 200.0, 2000.0, 0),
            feature(3, 700.0, 300.0, 1500.0, 0),
        ]),
        FeatureMap::from_features(vec![
            feature(4, 500.1, 101.0, 1100.0, 0),
            feature(5, 600.1, 201.0, 1900.0, 0),
            feature(6, 900.0, 900.0, 500.0, 0),
        ]),
    ]
}

#[test]
fn oracle_class_da() {
    let mut a = algorithm(&[("mz_unit", string("Da"))]);
    let mut out = ConsensusMap::new();
    a.group(&class_maps(), &mut out).unwrap();
    compare("class_da", &out);
}

#[test]
fn oracle_class_ppm() {
    let mut a = algorithm(&[
        ("mz_unit", string("ppm")),
        ("max_mz_shift_ppm", ParamValue::Float(500.0)),
    ]);
    let mut out = ConsensusMap::new();
    a.group(&class_maps(), &mut out).unwrap();
    compare("class_ppm", &out);
}

#[test]
fn oracle_class_dc() {
    let maps = vec![
        FeatureMap::from_features(vec![feature(11, 250.5, 100.0, 1000.0, 2)]),
        FeatureMap::from_features(vec![feature(12, 500.0, 100.0, 900.0, 1)]),
    ];
    let mut a = algorithm(&[("mz_unit", string("Da")), ("decharge_mz", string("true"))]);
    let mut out = ConsensusMap::new();
    a.group(&maps, &mut out).unwrap();
    compare("class_dc", &out);
}

#[test]
fn oracle_class_consensus() {
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
    let mut a = algorithm(&[("mz_unit", string("Da"))]);
    let mut out = ConsensusMap::new();
    a.group_consensus(&maps, &mut out).unwrap();
    compare("class_consensus", &out);
}
