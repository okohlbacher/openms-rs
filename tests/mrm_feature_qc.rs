// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! MRMFeatureQC, MRMFeatureQCFile, MRMFeaturePicker and MRMFeaturePickerFile.
//!
//! Expected values are transcribed from the upstream class tests at the core
//! pin `bc9cc12` (evidence tier 3): `MRMFeatureQCFile_test.cpp` and
//! `MRMFeaturePickerFile_test.cpp`, with their data files. Tests marked
//! "native" check a behaviour the upstream tests do not reach (tier 4).

use openms::Error;
use openms::analysis::openswath::mrm_feature_picker::ComponentParams;
use openms::analysis::openswath::mrm_feature_qc::{ComponentGroupQCs, ComponentQCs, MRMFeatureQC};
use openms::format::mrm_feature_picker_file::MRMFeaturePickerFile;
use openms::format::mrm_feature_qc_file::{HeaderMap, MRMFeatureQCFile};
use openms::param::ParamValue;
use std::collections::BTreeMap;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn temp(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms_mrm_feature_qc_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

fn real_similar(a: f64, b: f64) {
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

fn strings(values: &[&str]) -> Vec<String> {
    values.iter().map(|s| (*s).to_owned()).collect()
}

fn headers(names: &[&str]) -> HeaderMap {
    names
        .iter()
        .enumerate()
        .map(|(i, n)| ((*n).to_owned(), i))
        .collect()
}

#[test]
fn defaults_are_the_member_initialisers() {
    // MRMFeatureQC.h: every lower bound 0.0 / 0, every upper bound 1e12 / 100.
    let c = ComponentQCs::default();
    assert_eq!(c.retention_time_l, 0.0);
    assert_eq!(c.retention_time_u, 1e12);
    assert_eq!(c.intensity_u, 1e12);
    assert_eq!(c.overall_quality_u, 1e12);
    let g = ComponentGroupQCs::default();
    assert_eq!((g.n_heavy_l, g.n_heavy_u), (0, 100));
    assert_eq!((g.n_transitions_l, g.n_transitions_u), (0, 100));
    assert_eq!(g.ion_ratio_u, 1e12);
    assert!(g.ion_ratio_pair_name_1.is_empty());
    assert!(MRMFeatureQC::default().component_group_pair_qcs.is_empty());
}

#[test]
fn qc_file_load() {
    // MRMFeatureQCFile_test.cpp, section load().
    let file = MRMFeatureQCFile;
    let mut qc = MRMFeatureQC::default();
    file.load(data("MRMFeatureQCFile_1.csv"), &mut qc, false)
        .unwrap();
    file.load(data("MRMFeatureQCFile_2.csv"), &mut qc, true)
        .unwrap();
    let c = &qc.component_qcs;
    let cg = &qc.component_group_qcs;
    assert_eq!(c.len(), 4);
    assert_eq!(c[0].component_name, "component1");
    assert_eq!(c[1].component_name, "component2");
    assert_eq!(c[2].component_name, "component3");
    // the row between component3 and component4 has an empty component_name
    assert_eq!(c[3].component_name, "component4");
    for (i, v) in [1.0, 3.0, 5.0, 0.0].into_iter().enumerate() {
        real_similar(c[i].retention_time_l, v);
    }
    for (i, v) in [2.0, 4.0, 6.0, 1e12].into_iter().enumerate() {
        real_similar(c[i].retention_time_u, v);
    }
    for (i, v) in [1000.0, 2000.0, 3000.0, 0.0].into_iter().enumerate() {
        real_similar(c[i].intensity_l, v);
    }
    for (i, v) in [1_000_000.0, 2_000_000.0, 3_000_000.0, 1e12]
        .into_iter()
        .enumerate()
    {
        real_similar(c[i].intensity_u, v);
    }
    for (i, v) in [2.0, 3.0, 4.0, 0.0].into_iter().enumerate() {
        real_similar(c[i].overall_quality_l, v);
    }
    for (i, v) in [5.0, 6.0, 7.0, 1e12].into_iter().enumerate() {
        real_similar(c[i].overall_quality_u, v);
    }
    for (i, v) in [1000.0, 2000.0, 3000.0, 0.0].into_iter().enumerate() {
        real_similar(c[i].meta_value_qc["peak_apex_int"].0, v);
    }
    for (i, v) in [1_000_000.0, 2_000_000.0, 3_000_000.0, 1e12]
        .into_iter()
        .enumerate()
    {
        real_similar(c[i].meta_value_qc["peak_apex_int"].1, v);
    }
    for (i, v) in [2.0, 5.0, 10.0, 0.0].into_iter().enumerate() {
        real_similar(c[i].meta_value_qc["sn_score"].0, v);
    }
    for (i, v) in [10.0, 20.0, 50.0, 1e12].into_iter().enumerate() {
        real_similar(c[i].meta_value_qc["sn_score"].1, v);
    }

    assert_eq!(cg.len(), 4);
    assert_eq!(cg[0].component_group_name, "componentGroup1");
    assert_eq!(cg[1].component_group_name, "componentGroup2");
    assert_eq!(cg[2].component_group_name, "componentGroup3");
    assert_eq!(cg[3].component_group_name, "componentGroup5");
    let counts = |q: &ComponentGroupQCs| {
        [
            q.n_heavy_l,
            q.n_heavy_u,
            q.n_light_l,
            q.n_light_u,
            q.n_detecting_l,
            q.n_detecting_u,
            q.n_quantifying_l,
            q.n_quantifying_u,
            q.n_identifying_l,
            q.n_identifying_u,
            q.n_transitions_l,
            q.n_transitions_u,
        ]
    };
    assert_eq!(counts(&cg[0]), [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);
    assert_eq!(counts(&cg[2]), [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]);
    assert_eq!(cg[0].ion_ratio_pair_name_1, "component1");
    assert_eq!(cg[2].ion_ratio_pair_name_1, "component5");
    assert_eq!(cg[0].ion_ratio_pair_name_2, "component2");
    assert_eq!(cg[2].ion_ratio_pair_name_2, "component6");
    real_similar(cg[0].ion_ratio_l, 0.5);
    real_similar(cg[2].ion_ratio_l, 2.5);
    real_similar(cg[0].ion_ratio_u, 0.6);
    real_similar(cg[2].ion_ratio_u, 2.6);
    assert_eq!(cg[0].ion_ratio_feature_name, "feature1");
    assert_eq!(cg[2].ion_ratio_feature_name, "feature3");
    for (i, v) in [1.0, 2.0, 3.0].into_iter().enumerate() {
        real_similar(cg[i].retention_time_l, v);
    }
    for (i, v) in [2.0, 3.0, 4.0].into_iter().enumerate() {
        real_similar(cg[i].retention_time_u, v);
    }
    for (i, v) in [1000.0, 1001.0, 1002.0].into_iter().enumerate() {
        real_similar(cg[i].intensity_l, v);
    }
    for (i, v) in [1_000_000.0, 1_000_001.0, 1_000_002.0]
        .into_iter()
        .enumerate()
    {
        real_similar(cg[i].intensity_u, v);
    }
    for (i, v) in [2.0, 3.0, 4.0].into_iter().enumerate() {
        real_similar(cg[i].overall_quality_l, v);
    }
    for (i, v) in [5.0, 6.0, 7.0].into_iter().enumerate() {
        real_similar(cg[i].overall_quality_u, v);
    }
    real_similar(cg[0].meta_value_qc["peak_apex_int"].0, 1000.0);
    real_similar(cg[2].meta_value_qc["peak_apex_int"].0, 1002.0);
    real_similar(cg[0].meta_value_qc["peak_apex_int"].1, 1_000_000.0);
    real_similar(cg[2].meta_value_qc["peak_apex_int"].1, 1_000_002.0);
    real_similar(cg[0].meta_value_qc["sn_score"].0, 2.0);
    real_similar(cg[2].meta_value_qc["sn_score"].0, 10.0);
    real_similar(cg[0].meta_value_qc["sn_score"].1, 10.0);
    real_similar(cg[2].meta_value_qc["sn_score"].1, 50.0);
    // componentGroup5 is a row of empty cells: every default
    let d = &cg[3];
    assert_eq!(counts(d), [0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100]);
    assert_eq!(d.ion_ratio_pair_name_1, "");
    assert_eq!(d.ion_ratio_pair_name_2, "");
    real_similar(d.ion_ratio_l, 0.0);
    real_similar(d.ion_ratio_u, 1e12);
    assert_eq!(d.ion_ratio_feature_name, "");
    real_similar(d.retention_time_l, 0.0);
    real_similar(d.retention_time_u, 1e12);
    real_similar(d.intensity_l, 0.0);
    real_similar(d.intensity_u, 1e12);
    real_similar(d.overall_quality_l, 0.0);
    real_similar(d.overall_quality_u, 1e12);
    real_similar(d.meta_value_qc["peak_apex_int"].0, 0.0);
    real_similar(d.meta_value_qc["peak_apex_int"].1, 1e12);
    real_similar(d.meta_value_qc["sn_score"].0, 0.0);
    real_similar(d.meta_value_qc["sn_score"].1, 1e12);
}

#[test]
fn qc_file_store_round_trips() {
    // MRMFeatureQCFile_test.cpp, section store().
    let file = MRMFeatureQCFile;
    let empty = MRMFeatureQC::default();
    let comp_path = temp("store_components.csv");
    let group_path = temp("store_groups.csv");
    let _ = std::fs::remove_file(&comp_path);
    let _ = std::fs::remove_file(&group_path);
    // an empty list writes no file, as in the source
    file.store(&comp_path, &empty, false).unwrap();
    file.store(&group_path, &empty, true).unwrap();
    assert!(!comp_path.exists() && !group_path.exists());

    let mut qc = MRMFeatureQC::default();
    file.load(data("MRMFeatureQCFile_1.csv"), &mut qc, false)
        .unwrap();
    file.load(data("MRMFeatureQCFile_2.csv"), &mut qc, true)
        .unwrap();
    file.store(&comp_path, &qc, false).unwrap();
    file.store(&group_path, &qc, true).unwrap();
    let mut again = MRMFeatureQC::default();
    file.load(&comp_path, &mut again, false).unwrap();
    file.load(&group_path, &mut again, true).unwrap();
    assert_eq!(qc.component_qcs.len(), again.component_qcs.len());
    assert_eq!(
        qc.component_group_qcs.len(),
        again.component_group_qcs.len()
    );
    for (a, b) in qc.component_qcs.iter().zip(&again.component_qcs) {
        assert_eq!(a, b);
    }
    for (a, b) in qc
        .component_group_qcs
        .iter()
        .zip(&again.component_group_qcs)
    {
        assert_eq!(a, b);
    }
}

#[test]
fn qc_file_store_matches_the_retained_component_file() {
    // MRMFeatureQCFile_1_test.csv ships with the upstream class-test data
    // (the class test does not read it). Its values are those of
    // MRMFeatureQCFile_1.csv spelled as StringUtils::toStr writes a double
    // (`1.0e06`); it is the reference for the stored text.
    let file = MRMFeatureQCFile;
    let mut qc = MRMFeatureQC::default();
    file.load(data("MRMFeatureQCFile_1.csv"), &mut qc, false)
        .unwrap();
    let path = temp("components_text.csv");
    file.store(&path, &qc, false).unwrap();
    let written = std::fs::read_to_string(&path).unwrap();
    let expected = std::fs::read_to_string(data("MRMFeatureQCFile_1_test.csv")).unwrap();
    let norm = |t: &str| {
        t.lines()
            .map(|l| l.trim_end().to_owned())
            .collect::<Vec<_>>()
    };
    assert_eq!(norm(&written), norm(&expected));
}

#[test]
fn qc_file_push_values_from_line_component() {
    // MRMFeatureQCFile_test.cpp, pushValuesFromLine_ (ComponentQCs).
    let h = headers(&[
        "component_name",
        "retention_time_l",
        "retention_time_u",
        "intensity_l",
        "intensity_u",
        "overall_quality_l",
        "overall_quality_u",
        "metaValue_peak_apex_int_l",
        "metaValue_peak_apex_int_u",
        "metaValue_sn_score_l",
        "metaValue_sn_score_u",
    ]);
    let sl1 = strings(&[
        "component1",
        "0.1",
        "0.2",
        "0.3",
        "0.4",
        "0.5",
        "0.6",
        "0.7",
        "0.8",
        "0.9",
        "1.0",
    ]);
    let sl2 = strings(&[
        "", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0",
    ]);
    let sl3 = strings(&["component3", "", "", "", "", "", "", "", "", "", ""]);
    let file = MRMFeatureQCFile;
    let mut c = Vec::new();
    file.push_values_from_line_component(&sl1, &h, &mut c)
        .unwrap();
    assert_eq!(c.len(), 1);
    assert_eq!(c[0].component_name, "component1");
    real_similar(c[0].retention_time_l, 0.1);
    real_similar(c[0].retention_time_u, 0.2);
    real_similar(c[0].intensity_l, 0.3);
    real_similar(c[0].intensity_u, 0.4);
    real_similar(c[0].overall_quality_l, 0.5);
    real_similar(c[0].overall_quality_u, 0.6);
    real_similar(c[0].meta_value_qc["peak_apex_int"].0, 0.7);
    real_similar(c[0].meta_value_qc["peak_apex_int"].1, 0.8);
    real_similar(c[0].meta_value_qc["sn_score"].0, 0.9);
    real_similar(c[0].meta_value_qc["sn_score"].1, 1.0);
    file.push_values_from_line_component(&sl2, &h, &mut c)
        .unwrap();
    assert_eq!(c.len(), 1);
    file.push_values_from_line_component(&sl3, &h, &mut c)
        .unwrap();
    assert_eq!(c.len(), 2);
    assert_eq!(c[1].component_name, "component3");
    real_similar(c[1].retention_time_l, 0.0);
    real_similar(c[1].retention_time_u, 1e12);
    real_similar(c[1].intensity_l, 0.0);
    real_similar(c[1].intensity_u, 1e12);
    real_similar(c[1].overall_quality_l, 0.0);
    real_similar(c[1].overall_quality_u, 1e12);
    real_similar(c[1].meta_value_qc["peak_apex_int"].0, 0.0);
    real_similar(c[1].meta_value_qc["peak_apex_int"].1, 1e12);
    real_similar(c[1].meta_value_qc["sn_score"].0, 0.0);
    real_similar(c[1].meta_value_qc["sn_score"].1, 1e12);
}

#[test]
fn qc_file_push_values_from_line_component_group() {
    // MRMFeatureQCFile_test.cpp, pushValuesFromLine_ (ComponentGroupQCs).
    let h = headers(&[
        "component_group_name",
        "n_heavy_l",
        "n_heavy_u",
        "n_light_l",
        "n_light_u",
        "n_detecting_l",
        "n_detecting_u",
        "n_quantifying_l",
        "n_quantifying_u",
        "n_identifying_l",
        "n_identifying_u",
        "n_transitions_l",
        "n_transitions_u",
        "ion_ratio_pair_name_1",
        "ion_ratio_pair_name_2",
        "ion_ratio_l",
        "ion_ratio_u",
        "ion_ratio_feature_name",
        "retention_time_l",
        "retention_time_u",
        "intensity_l",
        "intensity_u",
        "overall_quality_l",
        "overall_quality_u",
        "metaValue_peak_apex_int_l",
        "metaValue_peak_apex_int_u",
        "metaValue_sn_score_l",
        "metaValue_sn_score_u",
    ]);
    let sl1 = strings(&[
        "component_group_1",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "11",
        "12",
        "ionRatioPairName1",
        "ionRatioPairName2",
        "1.1",
        "1.2",
        "ionRatioFeatureName",
        "0.1",
        "0.2",
        "0.3",
        "0.4",
        "0.5",
        "0.6",
        "0.7",
        "0.8",
        "0.9",
        "1.0",
    ]);
    let mut sl2 = sl1.clone();
    sl2[0] = String::new();
    let mut sl3 = vec![String::new(); 28];
    sl3[0] = "component_group_3".into();
    let file = MRMFeatureQCFile;
    let mut cg = Vec::new();
    file.push_values_from_line_component_group(&sl1, &h, &mut cg)
        .unwrap();
    assert_eq!(cg.len(), 1);
    let q = &cg[0];
    assert_eq!(q.component_group_name, "component_group_1");
    assert_eq!(
        [
            q.n_heavy_l,
            q.n_heavy_u,
            q.n_light_l,
            q.n_light_u,
            q.n_detecting_l,
            q.n_detecting_u,
            q.n_quantifying_l,
            q.n_quantifying_u,
            q.n_identifying_l,
            q.n_identifying_u,
            q.n_transitions_l,
            q.n_transitions_u
        ],
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    );
    assert_eq!(q.ion_ratio_pair_name_1, "ionRatioPairName1");
    assert_eq!(q.ion_ratio_pair_name_2, "ionRatioPairName2");
    real_similar(q.ion_ratio_l, 1.1);
    real_similar(q.ion_ratio_u, 1.2);
    assert_eq!(q.ion_ratio_feature_name, "ionRatioFeatureName");
    real_similar(q.retention_time_l, 0.1);
    real_similar(q.retention_time_u, 0.2);
    real_similar(q.intensity_l, 0.3);
    real_similar(q.intensity_u, 0.4);
    real_similar(q.overall_quality_l, 0.5);
    real_similar(q.overall_quality_u, 0.6);
    real_similar(q.meta_value_qc["peak_apex_int"].0, 0.7);
    real_similar(q.meta_value_qc["peak_apex_int"].1, 0.8);
    real_similar(q.meta_value_qc["sn_score"].0, 0.9);
    real_similar(q.meta_value_qc["sn_score"].1, 1.0);
    file.push_values_from_line_component_group(&sl2, &h, &mut cg)
        .unwrap();
    assert_eq!(cg.len(), 1);
    file.push_values_from_line_component_group(&sl3, &h, &mut cg)
        .unwrap();
    assert_eq!(cg.len(), 2);
    let d = &cg[1];
    assert_eq!(d.component_group_name, "component_group_3");
    assert_eq!(
        [
            d.n_heavy_l,
            d.n_heavy_u,
            d.n_transitions_l,
            d.n_transitions_u
        ],
        [0, 100, 0, 100]
    );
    real_similar(d.ion_ratio_u, 1e12);
    real_similar(d.retention_time_u, 1e12);
    real_similar(d.meta_value_qc["sn_score"].1, 1e12);
}

#[test]
fn qc_file_set_pair_value() {
    // MRMFeatureQCFile_test.cpp, setPairValue_.
    let file = MRMFeatureQCFile;
    let mut m: BTreeMap<String, (f64, f64)> = BTreeMap::new();
    file.set_pair_value("meta1", "0.123", "u", &mut m).unwrap();
    assert_eq!(m.len(), 1);
    real_similar(m["meta1"].0, 0.0);
    real_similar(m["meta1"].1, 0.123);
    file.set_pair_value("meta1", "0.456", "l", &mut m).unwrap();
    assert_eq!(m.len(), 1);
    real_similar(m["meta1"].0, 0.456);
    real_similar(m["meta1"].1, 0.123);
    file.set_pair_value("meta1", "0.789", "u", &mut m).unwrap();
    real_similar(m["meta1"].0, 0.456);
    real_similar(m["meta1"].1, 0.789);
    file.set_pair_value("meta2", "0.111", "l", &mut m).unwrap();
    assert_eq!(m.len(), 2);
    real_similar(m["meta2"].0, 0.111);
    real_similar(m["meta2"].1, 1e12);
    file.set_pair_value("meta3", "0.222", "u", &mut m).unwrap();
    assert_eq!(m.len(), 3);
    real_similar(m["meta3"].0, 0.0);
    real_similar(m["meta3"].1, 0.222);
}

#[test]
fn qc_file_get_cast_value() {
    // MRMFeatureQCFile_test.cpp, the three getCastValue_ overloads.
    let file = MRMFeatureQCFile;
    let h = headers(&["component_group_name", "n_heavy_l", "n_heavy_u"]);
    let line1 = strings(&["componentgroup1", "2", "3"]);
    let line2 = strings(&["componentgroup2", "", "3"]);
    assert_eq!(
        file.get_cast_value_i32(&h, &line1, "n_heavy_l", 3).unwrap(),
        2
    );
    assert_eq!(
        file.get_cast_value_i32(&h, &line1, "n_light_l", 4).unwrap(),
        4
    );
    assert_eq!(
        file.get_cast_value_i32(&h, &line2, "n_heavy_l", 5).unwrap(),
        5
    );

    let h = headers(&["component_name", "retention_time_l", "retention_time_u"]);
    let line1 = strings(&["component1", "1.2", "1.3"]);
    let line2 = strings(&["component2", "", "1.3"]);
    assert_eq!(
        file.get_cast_value_f64(&h, &line1, "retention_time_l", 3.1)
            .unwrap(),
        1.2
    );
    assert_eq!(
        file.get_cast_value_f64(&h, &line1, "intensity_l", 4.1)
            .unwrap(),
        4.1
    );
    assert_eq!(
        file.get_cast_value_f64(&h, &line2, "retention_time_l", 5.1)
            .unwrap(),
        5.1
    );

    let h = headers(&["component_name", "ion_ratio_feature_name"]);
    let line1 = strings(&["component1", "name1"]);
    let line2 = strings(&["component2", ""]);
    assert_eq!(
        file.get_cast_value_string(&h, &line1, "ion_ratio_feature_name", "name30")
            .unwrap(),
        "name1"
    );
    assert_eq!(
        file.get_cast_value_string(&h, &line1, "intensity_l", "name30")
            .unwrap(),
        "name30"
    );
    assert_eq!(
        file.get_cast_value_string(&h, &line2, "ion_ratio_feature_name", "name30")
            .unwrap(),
        "name30"
    );
}

#[test]
fn qc_file_conversions_follow_std_stod_and_std_stoi() {
    // native: the source converts with std::stod / std::stoi, which read a
    // numeric prefix and throw on no number or on an out-of-range value.
    let file = MRMFeatureQCFile;
    let h = headers(&["v"]);
    let get = |text: &str| file.get_cast_value_f64(&h, &strings(&[text]), "v", -1.0);
    assert_eq!(get("1.5abc").unwrap(), 1.5);
    assert_eq!(get("  -2e2").unwrap(), -200.0);
    assert_eq!(get(".5").unwrap(), 0.5);
    assert!(get("abc").is_err());
    assert!(get("1e999").is_err());
    // strtod reads hexadecimal floats, and so does the Release build's
    // `ArrowIOHelpers::readMetaValues` (../oracle/F1a-arrow, d12: `0x1p3` is 8);
    // the one `std::stod` port reads them too.
    assert_eq!(get("0x10").unwrap(), 16.0);
    assert_eq!(get("0x1p3").unwrap(), 8.0);
    let geti = |text: &str| file.get_cast_value_i32(&h, &strings(&[text]), "v", -1);
    assert_eq!(geti("42xyz").unwrap(), 42);
    assert_eq!(geti(" -7").unwrap(), -7);
    assert!(geti("x").is_err());
    assert!(geti("3000000000").is_err());
}

#[test]
fn qc_file_refuses_a_short_row_and_misaligned_metavalues() {
    // native: the source reads past the end of a short row and writes past the
    // end of a row with more metavalues than the header; both are refused.
    let file = MRMFeatureQCFile;
    let h = headers(&["component_name", "retention_time_l"]);
    let mut c = Vec::new();
    let err = file
        .push_values_from_line_component(&strings(&["a"]), &h, &mut c)
        .unwrap_err();
    assert!(matches!(err, Error::InvalidValue(_)));

    let mut first = ComponentQCs {
        component_name: "a".into(),
        ..ComponentQCs::default()
    };
    first.meta_value_qc.insert("x".into(), (0.0, 1.0));
    let mut second = first.clone();
    second.component_name = "b".into();
    second.meta_value_qc.insert("y".into(), (0.0, 1.0));
    let qc = MRMFeatureQC {
        component_qcs: vec![first, second],
        ..MRMFeatureQC::default()
    };
    let err = file.store(temp("refused.csv"), &qc, false).unwrap_err();
    assert!(matches!(err, Error::InvalidValue(_)));
}

#[test]
fn qc_file_load_is_atomic_on_error() {
    // native: the source clears the list before parsing; a failed load here
    // leaves the destination untouched.
    let path = temp("bad_number.csv");
    std::fs::write(
        &path,
        "component_name,retention_time_l\ngood,1.0\nbad,not_a_number\n",
    )
    .unwrap();
    let mut qc = MRMFeatureQC::default();
    qc.component_qcs.push(ComponentQCs {
        component_name: "kept".into(),
        ..ComponentQCs::default()
    });
    assert!(MRMFeatureQCFile.load(&path, &mut qc, false).is_err());
    assert_eq!(qc.component_qcs.len(), 1);
    assert_eq!(qc.component_qcs[0].component_name, "kept");
}

#[test]
fn qc_file_missing_file_is_an_error() {
    let mut qc = MRMFeatureQC::default();
    let err = MRMFeatureQCFile
        .load(data("does_not_exist.csv"), &mut qc, false)
        .unwrap_err();
    assert!(matches!(err, Error::Io(_)), "{err}");
}

fn int(params: &ComponentParams, key: &str) -> i64 {
    match params.params.value(key).unwrap() {
        ParamValue::Integer(v) => *v,
        other => panic!("{key} is {other:?}"),
    }
}

#[test]
fn picker_file_load() {
    // MRMFeaturePickerFile_test.cpp, section load().
    let (cp, cgp) = MRMFeaturePickerFile
        .load(data("MRMFeaturePickerFile.csv"))
        .unwrap();
    assert_eq!(cp.len(), 11);
    assert_eq!(cgp.len(), 5);
    let float = |p: &openms::param::Param, k: &str| match p.value(k).unwrap() {
        ParamValue::Float(v) => *v,
        other => panic!("{k} is {other:?}"),
    };
    let text = |p: &openms::param::Param, k: &str| match p.value(k).unwrap() {
        ParamValue::String(v) => v.clone(),
        other => panic!("{k} is {other:?}"),
    };

    assert_eq!(cp[1].component_name, "arg-L.arg-L_1.Light");
    assert_eq!(cp[1].component_group_name, "arg-L");
    assert_eq!(int(&cp[1], "sgolay_frame_length"), 152);
    assert_eq!(int(&cp[1], "sgolay_polynomial_order"), 32);
    real_similar(float(&cp[1].params, "gauss_width"), 0.152);
    assert_eq!(text(&cp[1].params, "use_gauss"), "false");
    real_similar(float(&cp[1].params, "peak_width"), 0.12);
    real_similar(float(&cp[1].params, "signal_to_noise"), 0.012);
    real_similar(float(&cp[1].params, "sn_win_len"), 10002.0);
    assert_eq!(int(&cp[1], "sn_bin_count"), 302);
    assert_eq!(text(&cp[1].params, "write_sn_log_messages"), "false");
    assert_eq!(text(&cp[1].params, "remove_overlapping_peaks"), "true");
    assert_eq!(text(&cp[1].params, "method"), "corrected2");

    assert_eq!(cp[9].component_name, "ser-L.ser-L_2.Light");
    assert_eq!(cp[9].component_group_name, "ser-L");
    assert_eq!(int(&cp[9], "sgolay_frame_length"), 160);
    assert_eq!(int(&cp[9], "sgolay_polynomial_order"), 40);
    real_similar(float(&cp[9].params, "gauss_width"), 0.16);
    assert_eq!(text(&cp[9].params, "use_gauss"), "false");
    real_similar(float(&cp[9].params, "peak_width"), 0.2);
    real_similar(float(&cp[9].params, "signal_to_noise"), 0.02);
    real_similar(float(&cp[9].params, "sn_win_len"), 10010.0);
    assert_eq!(int(&cp[9], "sn_bin_count"), 310);
    assert_eq!(text(&cp[9].params, "write_sn_log_messages"), "false");
    assert_eq!(text(&cp[9].params, "remove_overlapping_peaks"), "true");
    assert_eq!(text(&cp[9].params, "method"), "corrected10");

    assert_eq!(cp[10].component_name, "component2");
    assert_eq!(cp[10].component_group_name, "group2");
    assert_eq!(int(&cp[10], "sgolay_polynomial_order"), 43);
    real_similar(float(&cp[10].params, "gauss_width"), 0.163);
    assert_eq!(text(&cp[10].params, "use_gauss"), "true");
    real_similar(float(&cp[10].params, "peak_width"), 0.23);
    real_similar(float(&cp[10].params, "signal_to_noise"), 0.023);
    real_similar(float(&cp[10].params, "sn_win_len"), 10013.0);
    assert_eq!(int(&cp[10], "sn_bin_count"), 313);
    assert_eq!(text(&cp[10].params, "write_sn_log_messages"), "true");
    assert_eq!(text(&cp[10].params, "remove_overlapping_peaks"), "false");
    assert_eq!(text(&cp[10].params, "method"), "corrected13");
    assert!(!cp[10].params.exists("sgolay_frame_length").unwrap());

    let g = |i: usize, k: &str| match cgp[i].params.value(k).unwrap() {
        ParamValue::Integer(v) => *v as f64,
        ParamValue::Float(v) => *v,
        other => panic!("{k} is {other:?}"),
    };
    assert_eq!(cgp[1].component_group_name, "orn");
    assert_eq!(g(1, "stop_after_feature"), 6.0);
    real_similar(g(1, "stop_after_intensity_ratio"), 0.0006);
    real_similar(g(1, "min_peak_width"), -6.0);
    assert_eq!(text(&cgp[1].params, "peak_integration"), "smoothed3");
    assert_eq!(text(&cgp[1].params, "background_subtraction"), "none3");
    assert_eq!(text(&cgp[1].params, "recalculate_peaks"), "false");
    assert_eq!(text(&cgp[1].params, "use_precursors"), "true");
    real_similar(g(1, "recalculate_peaks_max_z"), 3.0);
    real_similar(g(1, "minimal_quality"), -10003.0);
    real_similar(g(1, "resample_boundary"), 0.03);
    assert_eq!(text(&cgp[1].params, "compute_peak_quality"), "true");
    assert_eq!(text(&cgp[1].params, "compute_peak_shape_metrics"), "false");

    assert_eq!(cgp[3].component_group_name, "ser-L");
    assert_eq!(g(3, "stop_after_feature"), 16.0);
    real_similar(g(3, "stop_after_intensity_ratio"), 0.0016);
    real_similar(g(3, "min_peak_width"), -16.0);
    assert_eq!(text(&cgp[3].params, "peak_integration"), "smoothed8");
    assert_eq!(text(&cgp[3].params, "background_subtraction"), "none8");
    assert_eq!(text(&cgp[3].params, "recalculate_peaks"), "true");
    assert_eq!(text(&cgp[3].params, "use_precursors"), "false");
    real_similar(g(3, "recalculate_peaks_max_z"), 8.0);
    real_similar(g(3, "minimal_quality"), -10008.0);
    real_similar(g(3, "resample_boundary"), 0.08);
    assert_eq!(text(&cgp[3].params, "compute_peak_quality"), "false");
    assert_eq!(text(&cgp[3].params, "compute_peak_shape_metrics"), "true");

    assert_eq!(cgp[4].component_group_name, "group2");
    real_similar(g(4, "stop_after_intensity_ratio"), 0.0026);
    real_similar(g(4, "min_peak_width"), -26.0);
    assert_eq!(text(&cgp[4].params, "peak_integration"), "smoothed13");
    assert_eq!(text(&cgp[4].params, "background_subtraction"), "none13");
    assert_eq!(text(&cgp[4].params, "recalculate_peaks"), "false");
    assert_eq!(text(&cgp[4].params, "use_precursors"), "true");
    real_similar(g(4, "recalculate_peaks_max_z"), 13.0);
    real_similar(g(4, "minimal_quality"), -10013.0);
    real_similar(g(4, "resample_boundary"), 0.13);
    assert_eq!(text(&cgp[4].params, "compute_peak_quality"), "true");
    assert_eq!(text(&cgp[4].params, "compute_peak_shape_metrics"), "false");
    assert!(!cgp[4].params.exists("stop_after_feature").unwrap());
}

#[test]
fn picker_file_requires_the_name_columns() {
    // MRMFeaturePickerFile.h: MissingInformation when a name column is absent,
    // but only when the file has a header and at least one more line.
    let path = temp("picker_no_names.csv");
    std::fs::write(&path, "a,b\n1,2\n").unwrap();
    let err = MRMFeaturePickerFile.load(&path).unwrap_err();
    assert!(matches!(err, Error::MissingInformation(_)), "{err}");
    std::fs::write(&path, "a,b\n").unwrap();
    let (cp, cgp) = MRMFeaturePickerFile.load(&path).unwrap();
    assert!(cp.is_empty() && cgp.is_empty());
}

#[test]
fn picker_file_value_types() {
    // MRMFeaturePickerFile.h setCastValue_: the bool rule is exact-match, not
    // case-insensitive, and the unsigned cast truncates a double.
    let file = MRMFeaturePickerFile;
    let mut p = openms::param::Param::new();
    file.set_cast_value("use_gauss", "True", &mut p).unwrap();
    assert_eq!(
        p.value("use_gauss").unwrap(),
        &ParamValue::String("false".into())
    );
    file.set_cast_value("sn_bin_count", "301.9", &mut p)
        .unwrap();
    assert_eq!(p.value("sn_bin_count").unwrap(), &ParamValue::Integer(301));
    file.set_cast_value("stop_after_feature", "-3", &mut p)
        .unwrap();
    assert_eq!(
        p.value("stop_after_feature").unwrap(),
        &ParamValue::Integer(-3)
    );
    file.set_cast_value("peak_width", "", &mut p).unwrap();
    assert!(!p.exists("peak_width").unwrap());
    assert!(file.set_cast_value("peak_width", "x", &mut p).is_err());
    assert!(
        file.set_cast_value("stop_after_feature", "1.5", &mut p)
            .is_err()
    );
}
