// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Evidence tier 1 for `FeatureGroupingAlgorithmKD`, `ClusterProxyKD` and
//! `MapAlignmentAlgorithmKD`: a differential against the executed C++
//! Release build.
//!
//! Oracle `M2-labeled` (`oracle/M2-labeled/driver.cpp`, cases `kddefaults`
//! and `kd1` to `kd8`) ran on ibminode06 (Linux x86_64) against the Release
//! install `openms4-release-bc9cc12-c19e494-174b576` (conda-forge g++
//! 14.4.0), twice, with identical output, retained verbatim:
//!
//! - `kd_grouping_oracle_defaults.ini`: `FeatureGroupingAlgorithmKD()
//!   .getParameters()` stored with `ParamXMLFile`;
//! - `kd_grouping_oracle_kd<n>.tsv`: the algorithm path of
//!   `TOPP_FeatureLinkerUnlabeledKD_<n>` (`topp/CMakeLists.txt:2927-2950` at
//!   test-data 0cb15f2): the defaults updated with the case's `.ini`
//!   algorithm section (`Param::update`, which drops the `None` merging
//!   values of `_4_parameters.ini` as invalid), the case's command-line
//!   overrides, inputs loaded as `FeatureLinkerBase.cpp` loads them (no
//!   subordinates or convex hulls; meta values cleared except
//!   `dc_charge_adducts` and `Group`), then `group`. Cases 4 to 8 warp with
//!   LOWESS; case 2 links two consensus maps; cases 4 to 7 carry adducts.
//!
//! Every consensus feature is printed as `kdcf RT m/z intensity quality
//! charge size handles... best-ion linked-groups peptide-count` with `%.17g`.
//! Compared: the counts, every handle's map index and unique ID exactly, the
//! RT and m/z to a relative 1e-12, intensity and quality (`float` in the
//! source) to one `f32` unit, the charge, the `best ion` and `LinkedGroups`
//! annotations exactly, and the number of peptide identifications.

#![cfg(all(feature = "featurexml", feature = "consensusxml", feature = "paramxml"))]

use openms::analysis::mapmatching::feature_grouping_algorithm::FeatureGroupingAlgorithm;
use openms::analysis::mapmatching::feature_grouping_algorithm_kd::FeatureGroupingAlgorithmKD;
use openms::format::{consensusxml, featurexml, paramxml};
use openms::kernel::features::{ConsensusMap, FeatureMap};
use openms::metadata::MetaValue;
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

/// The parameters of case `n`, as the oracle builds them.
fn parameters(ini: Option<&str>, overrides: &[(&str, &str)]) -> FeatureGroupingAlgorithmKD {
    let mut p = FeatureGroupingAlgorithmKD::defaults().unwrap();
    if let Some(ini) = ini {
        let full = paramxml::load(format!("{DIR}/{ini}")).unwrap();
        let section = full
            .copy("FeatureLinkerUnlabeledKD:1:algorithm:", true)
            .unwrap();
        p.update(&section, false).unwrap();
    }
    for (key, value) in overrides {
        set(&mut p, key, ParamValue::String((*value).into()));
    }
    FeatureGroupingAlgorithmKD::from_param(&p).unwrap()
}

fn rel(actual: f64, expected: f64, relative: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "{what}: produced {actual:e}, oracle {expected:e}"
    );
}

fn text(value: Option<&MetaValue>) -> String {
    match value {
        None => "-".into(),
        Some(v) => match v.as_string_list() {
            Ok(list) => list.join(","),
            Err(_) => v.as_str().unwrap().to_string(),
        },
    }
}

fn compare(case: &str, out: &ConsensusMap) {
    let expected = std::fs::read_to_string(format!("{DIR}/kd_grouping_oracle_{case}.tsv")).unwrap();
    let mut lines = expected.lines();
    let size: Vec<usize> = lines
        .next()
        .unwrap()
        .split('\t')
        .skip(1)
        .map(|x| x.parse().unwrap())
        .collect();
    assert_eq!(
        vec![
            out.len(),
            out.protein_identifications.len(),
            out.unassigned_peptide_identifications.len()
        ],
        size,
        "{case}: sizes"
    );
    let float_rel = f64::from(f32::EPSILON);
    for (index, (produced, line)) in out.features.iter().zip(lines).enumerate() {
        let at = format!("{case}, consensus feature {index}");
        let f: Vec<&str> = line.split('\t').collect();
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
        rel(
            f64::from(produced.quality),
            f[4].parse().unwrap(),
            float_rel,
            &format!("{at} quality"),
        );
        assert_eq!(
            produced.charge,
            f[5].parse::<i32>().unwrap(),
            "{at}: charge"
        );
        assert_eq!(
            text(produced.metadata.get("best ion")),
            f[7 + n],
            "{at}: best ion"
        );
        assert_eq!(
            text(produced.metadata.get("LinkedGroups")),
            f[8 + n],
            "{at}: linked groups"
        );
        assert_eq!(
            produced.peptide_identifications.len(),
            f[9 + n].parse::<usize>().unwrap(),
            "{at}: peptide identifications"
        );
    }
}

fn run_features(case: &str, algorithm: &mut FeatureGroupingAlgorithmKD, inputs: &[&str]) {
    let maps: Vec<FeatureMap> = inputs.iter().map(|f| load(f)).collect();
    let mut out = ConsensusMap::new();
    algorithm.group(&maps, &mut out).unwrap();
    compare(case, &out);
}

const UNLABELED_1: [&str; 3] = [
    "feature_linker_unlabeled_1_input1.featureXML",
    "feature_linker_unlabeled_1_input2.featureXML",
    "feature_linker_unlabeled_1_input3.featureXML",
];
const DC: [&str; 5] = [
    "feature_linker_unlabeled_kd_dc_input1.featureXML",
    "feature_linker_unlabeled_kd_dc_input2.featureXML",
    "feature_linker_unlabeled_kd_dc_input3.featureXML",
    "feature_linker_unlabeled_kd_dc_input1.featureXML",
    "feature_linker_unlabeled_kd_dc_input2.featureXML",
];

/// `ParamXMLFile` stores a double with 15 significant digits
/// (`LOWESS:span`, 2/3, is stored as `0.666666666666667`); every other
/// default is exactly representable in that form. The produced tree is
/// compared after the same rounding of that one value, and the unrounded
/// value is checked against the source literal `2/3.0`.
#[test]
fn oracle_kd_defaults() {
    let expected = paramxml::load(format!("{DIR}/kd_grouping_oracle_defaults.ini")).unwrap();
    let mut produced = FeatureGroupingAlgorithmKD::defaults().unwrap();
    let span = produced.value("LOWESS:span").unwrap().to_f64().unwrap();
    assert_eq!(span, 2.0 / 3.0);
    let stored: f64 = format!("{span:.14e}").parse().unwrap();
    set(&mut produced, "LOWESS:span", ParamValue::Float(stored));
    assert_eq!(produced, expected);
}

#[test]
fn oracle_kd_1() {
    let mut a = parameters(Some("feature_linker_unlabeled_kd_1_parameters.ini"), &[]);
    run_features("kd1", &mut a, &UNLABELED_1);
}

#[test]
fn oracle_kd_2_consensus_input() {
    let mut a = parameters(Some("feature_linker_unlabeled_kd_2_parameters.ini"), &[]);
    let input = consensusxml::load(format!(
        "{DIR}/feature_linker_unlabeled_kd_1_output.consensusXML"
    ))
    .unwrap();
    let maps = vec![input.clone(), input];
    let mut out = ConsensusMap::new();
    a.group_consensus(&maps, &mut out).unwrap();
    compare("kd2", &out);
}

#[test]
fn oracle_kd_3() {
    let mut a = parameters(Some("feature_linker_unlabeled_kd_3_parameters.ini"), &[]);
    run_features(
        "kd3",
        &mut a,
        &[
            "feature_linker_unlabeled_qt_3_input1.featureXML",
            "feature_linker_unlabeled_qt_3_input2.featureXML",
        ],
    );
}

fn dc_case(case: &str, charge: &str, adduct: &str) {
    let mut a = parameters(
        Some("feature_linker_unlabeled_kd_4_parameters.ini"),
        &[
            ("link:charge_merging", charge),
            ("link:adduct_merging", adduct),
        ],
    );
    run_features(case, &mut a, &DC);
}

#[test]
fn oracle_kd_4_warp_identical_charges() {
    dc_case("kd4", "Identical", "Any");
}

#[test]
fn oracle_kd_5_warp_charge_zero() {
    dc_case("kd5", "With_charge_zero", "Any");
}

#[test]
fn oracle_kd_6_warp_unknown_adducts() {
    dc_case("kd6", "Any", "With_unknown_adducts");
}

#[test]
fn oracle_kd_7_warp_identical_adducts() {
    dc_case("kd7", "Any", "Identical");
}

#[test]
fn oracle_kd_8_defaults_any_charge() {
    let mut a = parameters(None, &[("link:charge_merging", "Any")]);
    run_features(
        "kd8",
        &mut a,
        &[
            "feature_linker_unlabeled_kd_charge_any_input1.featureXML",
            "feature_linker_unlabeled_kd_charge_any_input2.featureXML",
        ],
    );
}
