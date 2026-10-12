// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Evidence tier 1 for the extraction order of `QTClusterFinder` and
//! `FeatureGroupingAlgorithmQT`, in particular among clusters of exactly equal
//! quality: a differential against the executed C++ Release build.
//!
//! Oracle `M2-qt` (`oracle/M2-qt/driver.cpp`) ran on ibminode06 (Linux
//! x86_64) against the Release install
//! `openms4-release-bc9cc12-c19e494-174b576` (conda-forge g++ 14.4.0,
//! libOpenMS's core flags), twice, with identical output.
//! `qt_instances.txt` holds 17 generated instances (`oracle/M2-qt/
//! gen_qt_instances.py`, seed 20261012): realistic noisy copies of a base map
//! under Da and ppm tolerances with several partitions, dense competition for
//! few cells, integer lattices with equal intensities (massive quality ties),
//! lattices with sub-tolerance jitter, exact duplicates, mirror pairs and
//! equidistant candidates. `qt_oracle.tsv` has, per instance, the consensus
//! features in the order `QTClusterFinder::run` extracts them (`raw`) and
//! after `FeatureGroupingAlgorithmQT::group`'s canonical sorts (`grp`), and
//! for the five `TOPP_FeatureLinkerUnlabeledQT_*` algorithm paths on the
//! retained inputs the same two views.
//!
//! Compared: the counts, every consensus feature's handles (map and unique
//! ID) in order, the quality (a `float` in the source) exactly, and RT, m/z
//! and intensity to 1e-12 relative and one `f32` unit.

#![cfg(all(feature = "featurexml", feature = "paramxml"))]

use std::collections::BTreeMap;

use openms::analysis::mapmatching::feature_grouping_algorithm::FeatureGroupingAlgorithm;
use openms::analysis::mapmatching::feature_grouping_algorithm_qt::FeatureGroupingAlgorithmQT;
use openms::analysis::mapmatching::qt_cluster_finder::QTClusterFinder;
use openms::format::{featurexml, paramxml};
use openms::kernel::features::{ConsensusMap, Feature, FeatureMap};
use openms::param::{Param, ParamValue};

const DIR: &str = "tests/data";

fn set(param: &mut Param, key: &str, value: ParamValue) {
    let description = param.description(key).unwrap().to_string();
    let tags: Vec<String> = param.tags(key).unwrap().iter().cloned().collect();
    param.set_value(key, value, &description, &tags).unwrap();
}

struct Instance {
    rt_tol: f64,
    mz_tol: f64,
    unit: String,
    parts: i64,
    maps: Vec<FeatureMap>,
}

fn instances() -> BTreeMap<String, Instance> {
    let text = std::fs::read_to_string(format!("{DIR}/qt_instances.txt")).unwrap();
    let mut out = BTreeMap::new();
    let mut name = String::new();
    let mut current: Option<Instance> = None;
    for line in text.lines() {
        let f: Vec<&str> = line.split_whitespace().collect();
        match f[0] {
            "instance" => {
                if let Some(done) = current.take() {
                    out.insert(name.clone(), done);
                }
                name = f[1].to_string();
            }
            "params" => {
                current = Some(Instance {
                    rt_tol: f[1].parse().unwrap(),
                    mz_tol: f[2].parse().unwrap(),
                    unit: f[3].to_string(),
                    parts: f[4].parse().unwrap(),
                    maps: Vec::new(),
                });
            }
            "F" => {
                let instance = current.as_mut().unwrap();
                let map: usize = f[1].parse().unwrap();
                while instance.maps.len() <= map {
                    instance.maps.push(FeatureMap::new());
                }
                let mut feature = Feature::new(
                    f[2].parse().unwrap(),
                    f[3].parse().unwrap(),
                    f[4].parse::<f64>().unwrap() as f32,
                );
                feature.unique_id = f[5].parse().unwrap();
                feature.charge = f[6].parse().unwrap();
                instance.maps[map].features.push(feature);
            }
            other => panic!("unexpected tag {other}"),
        }
    }
    if let Some(done) = current.take() {
        out.insert(name, done);
    }
    out
}

#[derive(Default)]
struct View {
    features: Vec<Vec<String>>,
}

#[derive(Default)]
struct Oracle {
    raw: View,
    grp: View,
}

fn oracle() -> BTreeMap<String, Oracle> {
    let text = std::fs::read_to_string(format!("{DIR}/qt_oracle.tsv")).unwrap();
    let mut out: BTreeMap<String, Oracle> = BTreeMap::new();
    let mut name = String::new();
    for line in text.lines() {
        let f: Vec<String> = line.split('\t').map(str::to_string).collect();
        match f[0].as_str() {
            "case" => {
                name = f[1].clone();
                out.insert(name.clone(), Oracle::default());
            }
            "raw" => out.get_mut(&name).unwrap().raw.features.push(f),
            "grp" => out.get_mut(&name).unwrap().grp.features.push(f),
            _ => {}
        }
    }
    out
}

fn rel(actual: f64, expected: f64, relative: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "{what}: produced {actual:e}, oracle {expected:e}"
    );
}

fn compare(case: &str, view: &str, out: &ConsensusMap, expected: &View) {
    assert_eq!(
        out.len(),
        expected.features.len(),
        "{case} {view}: consensus feature count"
    );
    for (index, (produced, f)) in out.features.iter().zip(&expected.features).enumerate() {
        let at = format!("{case} {view}, consensus feature {index}");
        let n: usize = f[6].parse().unwrap();
        let keys: Vec<String> = produced
            .handles()
            .iter()
            .map(|h| format!("{}:{}", h.map_index, h.unique_id))
            .collect();
        assert_eq!(keys, f[7..7 + n].to_vec(), "{at}: handles");
        assert_eq!(
            produced.quality,
            f[1].parse::<f64>().unwrap() as f32,
            "{at}: quality"
        );
        rel(
            produced.rt,
            f[2].parse().unwrap(),
            1e-12,
            &format!("{at} RT"),
        );
        rel(
            produced.mz,
            f[3].parse().unwrap(),
            1e-12,
            &format!("{at} m/z"),
        );
        rel(
            f64::from(produced.intensity),
            f[4].parse().unwrap(),
            f64::from(f32::EPSILON),
            &format!("{at} intensity"),
        );
        assert_eq!(
            produced.charge,
            f[5].parse::<i32>().unwrap(),
            "{at}: charge"
        );
    }
}

/// Which view of the oracle a test compares.
#[derive(Clone, Copy)]
enum Mode {
    /// `QTClusterFinder::run`: the order of extraction.
    Raw,
    /// `FeatureGroupingAlgorithmQT::group`: after the canonical sorts.
    Grouped,
}

fn check(case: &str, param: &Param, maps: &[FeatureMap], mode: Mode) {
    let all = oracle();
    let expected = &all[case];
    match mode {
        Mode::Raw => {
            let mut finder = QTClusterFinder::from_param(param).unwrap();
            let mut raw = ConsensusMap::new();
            finder.run_feature_maps(maps, &mut raw).unwrap();
            compare(case, "raw", &raw, &expected.raw);
        }
        Mode::Grouped => {
            let mut algorithm = FeatureGroupingAlgorithmQT::from_param(param).unwrap();
            let mut grp = ConsensusMap::new();
            algorithm.group(maps, &mut grp).unwrap();
            compare(case, "grp", &grp, &expected.grp);
        }
    }
}

fn instance_case(name: &str, mode: Mode) {
    let all = instances();
    let instance = &all[name];
    let mut p = QTClusterFinder::defaults().unwrap();
    set(
        &mut p,
        "distance_RT:max_difference",
        ParamValue::Float(instance.rt_tol),
    );
    set(
        &mut p,
        "distance_MZ:max_difference",
        ParamValue::Float(instance.mz_tol),
    );
    set(
        &mut p,
        "distance_MZ:unit",
        ParamValue::String(instance.unit.clone()),
    );
    set(&mut p, "nr_partitions", ParamValue::Integer(instance.parts));
    check(&format!("inst_{name}"), &p, &instance.maps, mode);
}

#[test]
fn oracle_covers_every_instance() {
    let inst = instances();
    let orc = oracle();
    assert_eq!(inst.len(), 17);
    let wanted: Vec<String> = inst.keys().map(|k| format!("inst_{k}")).collect();
    for name in &wanted {
        assert!(orc.contains_key(name), "{name}");
    }
}

macro_rules! instance_test {
    ($grouped:ident, $raw:ident, $name:literal) => {
        #[test]
        fn $grouped() {
            instance_case($name, Mode::Grouped);
        }
        #[test]
        fn $raw() {
            instance_case($name, Mode::Raw);
        }
    };
}

instance_test!(oracle_real3_da_grouped, oracle_real3_da_raw, "real3_da");
instance_test!(oracle_real3_ppm_grouped, oracle_real3_ppm_raw, "real3_ppm");
instance_test!(oracle_real5_da_grouped, oracle_real5_da_raw, "real5_da");
instance_test!(oracle_real5_ppm_grouped, oracle_real5_ppm_raw, "real5_ppm");
instance_test!(oracle_real8_da_grouped, oracle_real8_da_raw, "real8_da");
instance_test!(oracle_real8_ppm_grouped, oracle_real8_ppm_raw, "real8_ppm");
instance_test!(oracle_dense3_grouped, oracle_dense3_raw, "dense3");
instance_test!(oracle_dense4_grouped, oracle_dense4_raw, "dense4");
instance_test!(oracle_lat3_grouped, oracle_lat3_raw, "lat3");
instance_test!(oracle_lat4_grouped, oracle_lat4_raw, "lat4");
instance_test!(oracle_lat5_grouped, oracle_lat5_raw, "lat5");
instance_test!(oracle_latj4_grouped, oracle_latj4_raw, "latj4");
instance_test!(oracle_dups3_grouped, oracle_dups3_raw, "dups3");
instance_test!(oracle_mirror3_grouped, oracle_mirror3_raw, "mirror3");
instance_test!(oracle_equidist_grouped, oracle_equidist_raw, "equidist");
instance_test!(oracle_big4_ppm_grouped, oracle_big4_ppm_raw, "big4_ppm");
instance_test!(oracle_big4_da_grouped, oracle_big4_da_raw, "big4_da");

/// Load an input as the tool does (no subordinates or hulls; meta values
/// cleared except `dc_charge_adducts` and `Group`). The QT_5 inputs carry the
/// null date the readers refuse; see `mapmatching_qt_cluster_finder_topp.rs`.
fn load(file: &str) -> FeatureMap {
    const NULL_DATE: &str = "date=\"0000-00-00T00:00:00\"";
    const VALID_DATE: &str = "date=\"1970-01-01T00:00:00\"";
    let path = format!("{DIR}/{file}");
    let mut map = match featurexml::load(&path) {
        Ok(map) => map,
        Err(error) => {
            assert!(
                error
                    .to_string()
                    .contains("completion timestamp is outside calendar bounds"),
                "{path}: {error}"
            );
            let text = std::fs::read_to_string(&path).unwrap();
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

fn topp(n: u32, mode: Mode) {
    let (inputs, ini, use_ids, rt200): (Vec<String>, Option<&str>, bool, bool) = match n {
        1 | 4 => (
            (1..=3)
                .map(|i| format!("feature_linker_unlabeled_1_input{i}.featureXML"))
                .collect(),
            Some("feature_linker_unlabeled_qt_1_parameters.ini"),
            n == 4,
            false,
        ),
        3 => (
            (1..=2)
                .map(|i| format!("feature_linker_unlabeled_qt_3_input{i}.featureXML"))
                .collect(),
            Some("feature_linker_unlabeled_qt_3_parameters.ini"),
            false,
            false,
        ),
        _ => (
            (1..=3)
                .map(|i| format!("feature_linker_unlabeled_qt_5_input{i}.featureXML"))
                .collect(),
            if n == 5 {
                Some("feature_linker_unlabeled_qt_1_parameters.ini")
            } else {
                None
            },
            n == 6,
            true,
        ),
    };
    let mut p = QTClusterFinder::defaults().unwrap();
    if let Some(ini) = ini {
        let full = paramxml::load(format!("{DIR}/{ini}")).unwrap();
        let section = full
            .copy("FeatureLinkerUnlabeledQT:1:algorithm:", true)
            .unwrap();
        p.update(&section, false).unwrap();
    }
    if use_ids {
        set(
            &mut p,
            "use_identifications",
            ParamValue::String("true".into()),
        );
    }
    if rt200 {
        set(
            &mut p,
            "distance_RT:max_difference",
            ParamValue::Float(200.0),
        );
    }
    let maps: Vec<FeatureMap> = inputs.iter().map(|f| load(f)).collect();
    check(&format!("topp_qt{n}"), &p, &maps, mode);
}

#[test]
fn oracle_topp_qt1_grouped() {
    topp(1, Mode::Grouped);
}

#[test]
fn oracle_topp_qt1_raw() {
    topp(1, Mode::Raw);
}

#[test]
fn oracle_topp_qt3_grouped() {
    topp(3, Mode::Grouped);
}

#[test]
fn oracle_topp_qt3_raw() {
    topp(3, Mode::Raw);
}

#[test]
fn oracle_topp_qt4_grouped() {
    topp(4, Mode::Grouped);
}

#[test]
fn oracle_topp_qt4_raw() {
    topp(4, Mode::Raw);
}

#[test]
fn oracle_topp_qt5_grouped() {
    topp(5, Mode::Grouped);
}

#[test]
fn oracle_topp_qt5_raw() {
    topp(5, Mode::Raw);
}

#[test]
fn oracle_topp_qt6_grouped() {
    topp(6, Mode::Grouped);
}

#[test]
fn oracle_topp_qt6_raw() {
    topp(6, Mode::Raw);
}
