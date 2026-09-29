// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Evidence tier 1 for `FeatureGroupingAlgorithmUnlabeled`,
//! `StablePairFinder` and `FeatureDistance`: the algorithm path of the
//! upstream `TOPP_FeatureLinkerUnlabeled_1` to `_3` tests, compared with the
//! retained C++ outputs.
//!
//! Upstream definition: `topp/CMakeLists.txt:2879-2887` at test-data revision
//! 0cb15f2 runs `FeatureLinkerUnlabeled -test -ini <case>_parameters.ini` on the
//! case's featureXML inputs and compares with `<case>_output.consensusXML`
//! through `${DIFF} -whitelist "id=" "href="` (FuzzyDiff, `ratio` 1.01,
//! `absdiff` 0.01). The tool (`FeatureLinkerBase.cpp` at topp revision 174b576)
//! loads each input, sets the column headers, passes the `algorithm:` section
//! of its parameters to `FeatureGroupingAlgorithmUnlabeled` and calls
//! `group(maps, out_map)`; it then assigns fresh unique IDs and sorts
//! identifications, neither of which changes the order or the values compared
//! here. The tool is not ported yet; this test runs the same library calls.
//!
//! Compared: the number and order of consensus features, each one's handles
//! (map index and unique ID, exactly, plus RT, m/z and intensity) and each
//! consensus feature's RT, m/z, intensity, charge and quality. The upstream
//! whitelist skips every line containing `id=`, which includes every handle
//! line, so this is stricter than the upstream test. The retained file stores
//! numbers with six significant digits: a stored value differs from the value
//! printed by at most a relative 5e-6, which is the tolerance here, plus
//! one `f32` rounding unit where the value is a `float`. Qualities near zero
//! are compared with the upstream `absdiff` of 0.01 scaled down to 1e-6.
//!
//! Case 2 sets `distance_intensity:weight` to 0.5, so its result depends on
//! the normalising maximum intensity, which the source takes from stale
//! stored ranges after the first pairing.

#![cfg(all(feature = "featurexml", feature = "consensusxml", feature = "paramxml"))]

use openms::analysis::mapmatching::feature_grouping_algorithm::FeatureGroupingAlgorithm;
use openms::analysis::mapmatching::feature_grouping_algorithm_unlabeled::FeatureGroupingAlgorithmUnlabeled;
use openms::format::{consensusxml, featurexml, paramxml};
use openms::kernel::features::{ColumnHeader, ConsensusMap};

fn close(actual: f64, expected: f64, relative: f64, what: &str) {
    let bound = relative * expected.abs() + 1e-6;
    assert!(
        (actual - expected).abs() <= bound,
        "{what}: produced {actual}, retained {expected}"
    );
}

fn run_case(case: usize, inputs: usize) {
    let dir = "tests/data";
    let ini = paramxml::load(format!(
        "{dir}/feature_linker_unlabeled_{case}_parameters.ini"
    ))
    .unwrap();
    let algorithm_param = ini
        .copy("FeatureLinkerUnlabeled:1:algorithm:", true)
        .unwrap();
    let (mut algorithm, _warnings) =
        FeatureGroupingAlgorithmUnlabeled::from_param_with_warnings(&algorithm_param).unwrap();

    let mut out = ConsensusMap::new();
    let mut maps = Vec::new();
    for i in 0..inputs {
        let mut map = featurexml::load(format!(
            "{dir}/feature_linker_unlabeled_{case}_input{}.featureXML",
            i + 1
        ))
        .unwrap();
        out.column_headers.insert(
            i as u64,
            ColumnHeader {
                size: map.len(),
                unique_id: map.unique_id,
                ..ColumnHeader::default()
            },
        );
        for feature in &mut map.features {
            feature.subordinates.clear();
            feature.convex_hulls.clear();
        }
        maps.push(map);
    }
    algorithm.group(&maps, &mut out).unwrap();

    let expected = consensusxml::load(format!(
        "{dir}/feature_linker_unlabeled_{case}_output.consensusXML"
    ))
    .unwrap();
    assert_eq!(
        out.len(),
        expected.len(),
        "case {case}: consensus feature count"
    );
    let float_rel = 5.0e-6 + f64::from(f32::EPSILON);
    for (index, (produced, retained)) in out.features.iter().zip(&expected.features).enumerate() {
        let at = format!("case {case}, consensus feature {index}");
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

#[test]
fn topp_feature_linker_unlabeled_1() {
    run_case(1, 3);
}

#[test]
fn topp_feature_linker_unlabeled_2() {
    run_case(2, 3);
}

#[test]
fn topp_feature_linker_unlabeled_3() {
    run_case(3, 2);
}
