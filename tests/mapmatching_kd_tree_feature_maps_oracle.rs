// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Evidence tier 1 for `KDTreeFeatureMaps` and `FeatureMapping`: a
//! differential against the executed C++ Release build.
//!
//! Oracle `M2-labeled` (`oracle/M2-labeled/driver.cpp`, cases `kdorder` and
//! `featuremapping`) ran on ibminode06 (Linux x86_64) against the Release
//! install `openms4-release-bc9cc12-c19e494-174b576` (conda-forge g++
//! 14.4.0), twice, with identical output, retained verbatim:
//!
//! - `kd_tree_feature_maps_oracle_kdorder.tsv`: a tree over two maps (map 0
//!   `LabeledPairFinder.featureXML` from core bc9cc12, map 1
//!   `FeatureLinkerLabeled_2_input.featureXML` from test-data 0cb15f2); for
//!   every 37th feature `i`, `getNeighborhood(i, 30, 0.5, Da, include same
//!   map)`, `getNeighborhood(i, 60, 20, ppm, exclude own map, log fold change
//!   0.5)` and `queryRegion(rt(i) ± 200, mz(i) ± 5, ignoring map 0)`, each
//!   printed as `tag i count indices...` **in the order returned**, which is
//!   the pre-order of libkdtree++'s tree after `optimise()`.
//! - `kd_tree_feature_maps_oracle_featuremapping.tsv`: three maps and four
//!   MS2 spectra through `assignMS2IndexToFeature` (m/z tolerance 0.2 Da, RT
//!   5 s), printed by global kd-tree index. It shows the source's map-1
//!   exclusion: spectrum 0 at m/z 500.03 goes to map 0's feature at 500.0
//!   (index 0) although map 1's feature at 500.02 (index 2) is closer, and
//!   spectrum 1, whose only candidate is in map 1, is unassigned. Spectrum 2
//!   ties between index 1 and index 5 and goes to the first in query order.
//!
//! Compared exactly: every index list, in order.

#![cfg(feature = "featurexml")]

use openms::analysis::mapmatching::feature_mapping::{
    FeatureMappingInfo, assign_ms2_index_to_feature,
};
use openms::analysis::mapmatching::kd_tree_feature_maps::KDTreeFeatureMaps;
use openms::format::featurexml;
use openms::kernel::features::{Feature, FeatureMap};
use openms::kernel::{MSExperiment, MSSpectrum, Precursor};

fn read(name: &str) -> Vec<Vec<String>> {
    std::fs::read_to_string(format!("tests/data/kd_tree_feature_maps_oracle_{name}.tsv"))
        .unwrap()
        .lines()
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect()
}

fn line(tag: &str, i: usize, v: &[usize]) -> Vec<String> {
    let mut out = vec![tag.to_string(), i.to_string(), v.len().to_string()];
    out.extend(v.iter().map(usize::to_string));
    out
}

#[test]
fn oracle_kd_tree_query_order() {
    let expected = read("kdorder");
    let maps = vec![
        featurexml::load("tests/data/labeled_pair_finder.featureXML").unwrap(),
        featurexml::load("tests/data/feature_linker_labeled_2_input.featureXML").unwrap(),
    ];
    let kd = KDTreeFeatureMaps::from_maps(&maps).unwrap();
    let mut produced = vec![vec![
        "kd".to_string(),
        kd.size().to_string(),
        kd.tree_size().to_string(),
        kd.num_maps().to_string(),
    ]];
    for i in (0..kd.size()).step_by(37) {
        let (mut a, mut b, mut c) = (Vec::new(), Vec::new(), Vec::new());
        kd.get_neighborhood(i, &mut a, 30.0, 0.5, false, true, -1.0)
            .unwrap();
        kd.get_neighborhood(i, &mut b, 60.0, 20.0, true, false, 0.5)
            .unwrap();
        let (rt, mz) = (kd.rt(i).unwrap(), kd.mz(i).unwrap());
        kd.query_region(rt - 200.0, rt + 200.0, mz - 5.0, mz + 5.0, &mut c, Some(0));
        produced.push(line("nbA", i, &a));
        produced.push(line("nbB", i, &b));
        produced.push(line("qr", i, &c));
    }
    assert_eq!(produced.len(), expected.len());
    for (p, e) in produced.iter().zip(&expected) {
        assert_eq!(p, e);
    }
}

fn ms2(rt: f64, mz: f64) -> MSSpectrum {
    let mut s = MSSpectrum::default();
    s.ms_level = 2;
    s.rt = rt;
    s.precursors = vec![Precursor {
        mz,
        ..Precursor::default()
    }];
    s
}

#[test]
fn oracle_feature_mapping_map_one_exclusion() {
    let expected = read("featuremapping");
    let f = |rt: f64, mz: f64| Feature::new(rt, mz, 0.0);
    let maps = vec![
        FeatureMap::from_features(vec![f(100.0, 500.0), f(300.0, 700.0)]),
        FeatureMap::from_features(vec![f(100.0, 500.02), f(200.0, 600.0)]),
        FeatureMap::from_features(vec![f(100.0, 500.1), f(300.0, 700.0)]),
    ];
    let info = FeatureMappingInfo::new(&maps).unwrap();
    let exp = MSExperiment {
        spectra: vec![
            ms2(100.0, 500.03),
            ms2(200.0, 600.0),
            ms2(300.0, 700.0),
            ms2(100.0, 500.09),
        ],
        ..MSExperiment::default()
    };
    let r = assign_ms2_index_to_feature(&exp, &info, 0.2, 5.0, false).unwrap();
    let mut produced = Vec::new();
    for (feature, spectra) in &r.assigned_ms2 {
        let mut row = vec!["assigned".to_string(), feature.to_string()];
        row.extend(spectra.iter().map(usize::to_string));
        produced.push(row);
    }
    let mut row = vec!["unassigned".to_string()];
    row.extend(r.unassigned_ms2.iter().map(usize::to_string));
    produced.push(row);
    assert_eq!(produced, expected);
}
