// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Evidence tier 1 for `LabeledPairFinder` and
//! `FeatureGroupingAlgorithmLabeled`: a differential against the executed C++
//! Release build.
//!
//! Oracle `M2-labeled` (`oracle/M2-labeled/driver.cpp`) ran on ibminode06
//! (Linux x86_64) against the Release install
//! `openms4-release-bc9cc12-c19e494-174b576` built with conda-forge g++
//! 14.4.0, twice, with identical output. Its standard output is retained
//! verbatim:
//!
//! - `labeled_pair_finder_oracle_labeled1.tsv` and `..._labeled2.tsv`: the
//!   `TOPP_FeatureLinkerLabeled_1`/`_2` algorithm path (inputs from test-data
//!   0cb15f2, parameters transcribed from the case `.ini` files, feature meta
//!   values cleared as the tool does, headers 0 `light` and 1 `heavy`);
//! - `..._lpf250.tsv`: the estimating section of `LabeledPairFinder_test.cpp`
//!   (input `LabeledPairFinder.featureXML` from core bc9cc12, `mz_pair_dists`
//!   `[4]`, `mz_dev` 0.2, map index 0).
//!
//! For estimating cases the driver prints the fitted Gaussian (`fitted A x0
//! sigma`); every case prints each consensus feature as `cf RT m/z intensity
//! quality charge map:uid...` with `%.17g`.
//!
//! Compared: the consensus feature count, the protein and unassigned
//! identification counts, every handle's map index and unique ID exactly, the
//! RT and m/z to a relative 1e-12, intensity and quality (both `float` in the
//! source) to one `f32` unit, and the fitted RT distance and deviation to a
//! relative 1e-12.

#![cfg(feature = "featurexml")]

use openms::analysis::mapmatching::feature_grouping_algorithm_labeled::FeatureGroupingAlgorithmLabeled;
use openms::analysis::mapmatching::labeled_pair_finder::LabeledPairFinder;
use openms::concept::UniqueIdGenerator;
use openms::format::featurexml;
use openms::kernel::features::{ColumnHeader, ConsensusMap, FeatureMap};

fn rel(actual: f64, expected: f64, relative: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "{what}: produced {actual:e}, oracle {expected:e}"
    );
}

struct Oracle {
    fitted: Option<(f64, f64)>,
    size: (usize, usize, usize),
    features: Vec<Vec<String>>,
}

fn oracle(name: &str) -> Oracle {
    let text = std::fs::read_to_string(format!("tests/data/labeled_pair_finder_oracle_{name}.tsv"))
        .unwrap();
    let mut result = Oracle {
        fitted: None,
        size: (0, 0, 0),
        features: Vec::new(),
    };
    for line in text.lines() {
        let fields: Vec<String> = line.split('\t').map(str::to_string).collect();
        match fields[0].as_str() {
            "fitted" => {
                let x0: f64 = fields[2].parse().unwrap();
                let sigma: f64 = fields[3].parse().unwrap();
                result.fitted = Some((x0, sigma.abs() * 3.0));
            }
            "size" => {
                result.size = (
                    fields[1].parse().unwrap(),
                    fields[2].parse().unwrap(),
                    fields[3].parse().unwrap(),
                );
            }
            "cf" => result.features.push(fields),
            _ => {}
        }
    }
    result
}

fn light_heavy() -> ConsensusMap {
    let mut out = ConsensusMap::new();
    let header = ColumnHeader {
        filename: "x".into(),
        ..ColumnHeader::default()
    };
    out.column_headers.insert(0, header.clone());
    out.column_headers.insert(1, header);
    out.column_headers.get_mut(&0).unwrap().label = "light".into();
    out.column_headers.get_mut(&1).unwrap().label = "heavy".into();
    out
}

fn compare(name: &str, map: FeatureMap, finder: LabeledPairFinder) {
    let expected = oracle(name);
    let mut map = map;
    for feature in &mut map.features {
        feature.metadata.clear();
    }
    let algorithm = FeatureGroupingAlgorithmLabeled {
        pair_finder: finder,
    };
    let mut out = light_heavy();
    let mut generator = UniqueIdGenerator::from_seed(19_991_231_235_959);
    let estimate = algorithm
        .group_with_generator(std::slice::from_ref(&map), &mut out, &mut generator)
        .unwrap();
    match expected.fitted {
        Some((x0, deviation)) => {
            assert!(estimate.estimated, "{name}: estimated");
            rel(
                estimate.rt_pair_dist,
                x0,
                1e-12,
                &format!("{name} RT distance"),
            );
            rel(
                estimate.rt_dev_low,
                deviation,
                1e-12,
                &format!("{name} RT deviation"),
            );
            rel(
                estimate.rt_dev_high,
                deviation,
                1e-12,
                &format!("{name} RT deviation"),
            );
        }
        None => assert!(!estimate.estimated, "{name}: estimated"),
    }
    assert_eq!(
        (
            out.len(),
            out.protein_identifications.len(),
            out.unassigned_peptide_identifications.len()
        ),
        expected.size,
        "{name}: sizes"
    );
    let float_rel = f64::from(f32::EPSILON);
    for (index, (produced, fields)) in out.features.iter().zip(&expected.features).enumerate() {
        let at = format!("{name}, consensus feature {index}");
        rel(
            produced.rt,
            fields[1].parse().unwrap(),
            1e-12,
            &format!("{at} RT"),
        );
        rel(
            produced.mz,
            fields[2].parse().unwrap(),
            1e-12,
            &format!("{at} m/z"),
        );
        rel(
            f64::from(produced.intensity),
            fields[3].parse().unwrap(),
            float_rel,
            &format!("{at} intensity"),
        );
        rel(
            f64::from(produced.quality),
            fields[4].parse().unwrap(),
            float_rel,
            &format!("{at} quality"),
        );
        assert_eq!(
            produced.charge,
            fields[5].parse::<i32>().unwrap(),
            "{at} charge"
        );
        let keys: Vec<String> = produced
            .handles()
            .iter()
            .map(|h| format!("{}:{}", h.map_index, h.unique_id))
            .collect();
        assert_eq!(keys, fields[6..].to_vec(), "{at} handles");
    }
}

#[test]
fn oracle_feature_linker_labeled_1() {
    let finder = LabeledPairFinder {
        rt_estimate: false,
        rt_pair_dist: -1.0,
        rt_dev_low: 0.4,
        rt_dev_high: 0.4,
        mz_pair_dists: vec![8.0],
        mz_dev: 0.1,
        mrm: false,
    };
    let map = featurexml::load("tests/data/feature_linker_labeled_1_input.featureXML").unwrap();
    compare("labeled1", map, finder);
}

#[test]
fn oracle_feature_linker_labeled_2() {
    let finder = LabeledPairFinder {
        rt_estimate: true,
        rt_pair_dist: -121.0382,
        rt_dev_low: 123.3965,
        rt_dev_high: 123.3965,
        mz_pair_dists: vec![4.0],
        mz_dev: 0.1,
        mrm: false,
    };
    let map = featurexml::load("tests/data/feature_linker_labeled_2_input.featureXML").unwrap();
    compare("labeled2", map, finder);
}

#[test]
fn oracle_labeled_pair_finder_estimation() {
    let finder = LabeledPairFinder {
        mz_pair_dists: vec![4.0],
        mz_dev: 0.2,
        ..LabeledPairFinder::default()
    };
    let map = featurexml::load("tests/data/labeled_pair_finder.featureXML").unwrap();
    compare("lpf250", map, finder);
}
