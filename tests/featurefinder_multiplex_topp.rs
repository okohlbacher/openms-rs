// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Differential test of `FeatureFinderMultiplexAlgorithm` against the retained
//! outputs of the eleven upstream TOPP tests `TOPP_FeatureFinderMultiplex_1`
//! to `_11` (`topp/CMakeLists.txt:463-526` of the test-data package at
//! 0cb15f2).
//!
//! Each case runs the algorithm as the TOPP tool does — MS1 spectra of the
//! input, the tool's `.ini` section `FeatureFinderMultiplex:1:` without its
//! file and logging options — and compares the feature map and the consensus
//! map with the retained `_output.featureXML` and `_output.consensusXML`.
//!
//! Evidence tier 1 (retained C++ output). The comparison is the upstream
//! contract: FuzzyDiff with `topp/FuzzyDiff.ini` (a number matches when its
//! ratio to the expected one is at most 1.01 or their difference at most
//! 0.01), with the whitelisted `id=`, `href=` and `spectra_data` excluded.
//! Counts, charges, map indices and column labels are compared exactly.
//! None of the inputs carries FAIMS data, so the tool's per-compensation-
//! voltage split is a single group and the algorithm's maps are the tool's.

#![cfg(all(
    feature = "mzml",
    feature = "paramxml",
    feature = "featurexml",
    feature = "consensusxml"
))]

use std::path::PathBuf;

use openms::analysis::featurefinder::multiplex::FeatureFinderMultiplexAlgorithm;
use openms::concept::unique_id::UniqueIdGenerator;
use openms::format::{consensusxml, featurexml, mzml, paramxml};
use openms::kernel::{ConsensusMap, FeatureMap};

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/feature_finder_multiplex")
        .join(name)
}

/// FuzzyDiff's number test with `ratio` 1.01 and `absdiff` 0.01.
fn fuzzy(actual: f64, expected: f64) -> bool {
    if actual == expected || (actual - expected).abs() <= 0.01 {
        return true;
    }
    if actual == 0.0 || expected == 0.0 || actual.signum() != expected.signum() {
        return false;
    }
    let (a, b) = (actual.abs(), expected.abs());
    a.max(b) / a.min(b) <= 1.01
}

fn check(what: &str, actual: f64, expected: f64) {
    assert!(
        fuzzy(actual, expected),
        "{what}: {actual} does not match {expected}"
    );
}

fn run(case: u32) -> (FeatureMap, ConsensusMap) {
    let mut exp = mzml::load(data(&format!("FeatureFinderMultiplex_{case}_input.mzML"))).unwrap();
    exp.spectra.retain(|s| s.ms_level == 1);
    let ini = paramxml::load(data(&format!(
        "FeatureFinderMultiplex_{case}_parameters.ini"
    )))
    .unwrap();
    let mut param = ini.copy("FeatureFinderMultiplex:1:", true).unwrap();
    let tool_options = [
        "in",
        "out",
        "out_multiplets",
        "out_blacklist",
        "log",
        "debug",
        "threads",
        "no_progress",
        "force",
        "test",
    ];
    for key in tool_options {
        if param.exists(key).unwrap() {
            param.remove(key).unwrap();
        }
    }
    let mut algorithm = FeatureFinderMultiplexAlgorithm::from_param(&param).unwrap();
    algorithm.unique_ids = UniqueIdGenerator::from_seed(0);
    algorithm.run(&mut exp, true).unwrap();
    (
        algorithm.feature_map().clone(),
        algorithm.consensus_map().clone(),
    )
}

fn compare_features(case: u32, actual: &FeatureMap) {
    let expected = featurexml::load(data(&format!(
        "FeatureFinderMultiplex_{case}_output.featureXML"
    )))
    .unwrap();
    assert_eq!(
        actual.features.len(),
        expected.features.len(),
        "case {case}: feature count"
    );
    for (i, (a, e)) in actual.features.iter().zip(&expected.features).enumerate() {
        let at = |q: &str| format!("case {case} feature {i} {q}");
        check(&at("RT"), a.rt, e.rt);
        check(&at("m/z"), a.mz, e.mz);
        check(
            &at("intensity"),
            f64::from(a.intensity),
            f64::from(e.intensity),
        );
        check(&at("quality"), f64::from(a.quality), f64::from(e.quality));
        assert_eq!(a.charge, e.charge, "{}", at("charge"));
        assert_eq!(
            a.convex_hulls.len(),
            e.convex_hulls.len(),
            "{}",
            at("hull count")
        );
        for (h, (ha, he)) in a.convex_hulls.iter().zip(&e.convex_hulls).enumerate() {
            let (ba, be) = (ha.bounding_box().unwrap(), he.bounding_box().unwrap());
            let hull = |q: &str| at(&format!("hull {h} {q}"));
            check(&hull("min RT"), ba.min().rt, be.min().rt);
            check(&hull("max RT"), ba.max().rt, be.max().rt);
            check(&hull("min m/z"), ba.min().mz, be.min().mz);
            check(&hull("max m/z"), ba.max().mz, be.max().mz);
        }
    }
}

fn compare_consensus(case: u32, actual: &ConsensusMap) {
    let expected = consensusxml::load(data(&format!(
        "FeatureFinderMultiplex_{case}_output.consensusXML"
    )))
    .unwrap();
    assert_eq!(
        actual.experiment_type, expected.experiment_type,
        "case {case}: experiment type"
    );
    assert_eq!(
        actual.column_headers.keys().collect::<Vec<_>>(),
        expected.column_headers.keys().collect::<Vec<_>>(),
        "case {case}: map indices"
    );
    for (k, header) in &actual.column_headers {
        let other = &expected.column_headers[k];
        assert_eq!(header.label, other.label, "case {case} map {k} label");
        assert_eq!(header.size, other.size, "case {case} map {k} size");
        assert_eq!(
            header.metadata.get("channel_id"),
            other.metadata.get("channel_id"),
            "case {case} map {k} channel_id"
        );
    }
    assert_eq!(
        actual.features.len(),
        expected.features.len(),
        "case {case}: consensus count"
    );
    for (i, (a, e)) in actual.features.iter().zip(&expected.features).enumerate() {
        let at = |q: &str| format!("case {case} consensus {i} {q}");
        check(&at("RT"), a.rt, e.rt);
        check(&at("m/z"), a.mz, e.mz);
        check(
            &at("intensity"),
            f64::from(a.intensity),
            f64::from(e.intensity),
        );
        check(&at("quality"), f64::from(a.quality), f64::from(e.quality));
        assert_eq!(a.charge, e.charge, "{}", at("charge"));
        assert_eq!(
            a.handles().len(),
            e.handles().len(),
            "{}",
            at("handle count")
        );
        for (h, (ha, he)) in a.handles().iter().zip(e.handles()).enumerate() {
            let handle = |q: &str| at(&format!("handle {h} {q}"));
            assert_eq!(ha.map_index, he.map_index, "{}", handle("map"));
            check(&handle("RT"), ha.rt, he.rt);
            check(&handle("m/z"), ha.mz, he.mz);
            check(
                &handle("intensity"),
                f64::from(ha.intensity),
                f64::from(he.intensity),
            );
            assert_eq!(ha.charge, he.charge, "{}", handle("charge"));
        }
    }
}

fn case(n: u32, with_consensus: bool) {
    let (features, consensus) = run(n);
    compare_features(n, &features);
    if with_consensus {
        compare_consensus(n, &consensus);
    }
}

#[test]
fn topp_feature_finder_multiplex_1() {
    case(1, true);
}

#[test]
fn topp_feature_finder_multiplex_2() {
    case(2, true);
}

#[test]
fn topp_feature_finder_multiplex_3() {
    case(3, true);
}

#[test]
fn topp_feature_finder_multiplex_4() {
    case(4, true);
}

#[test]
fn topp_feature_finder_multiplex_5() {
    case(5, true);
}

#[test]
fn topp_feature_finder_multiplex_6() {
    case(6, true);
}

#[test]
fn topp_feature_finder_multiplex_7() {
    case(7, true);
}

#[test]
fn topp_feature_finder_multiplex_8() {
    case(8, true);
}

#[test]
fn topp_feature_finder_multiplex_9() {
    case(9, true);
}

#[test]
fn topp_feature_finder_multiplex_10() {
    case(10, true);
}

/// The upstream test writes no consensus map for case 11.
#[test]
fn topp_feature_finder_multiplex_11() {
    case(11, false);
}
