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
//! of its parameters to `FeatureGroupingAlgorithmUnlabeled` and, for
//! featureXML input (`FeatureLinkerUnlabeled.cpp`), links incrementally:
//! `setReference` with the input holding the most features, `addToGroup` for
//! every other input in order, `getResultMap`, then `sortByMZ`. It then
//! assigns fresh unique IDs and sorts identifications, neither of which
//! changes the order or the values compared here. The tool is not ported yet;
//! this test runs the same library calls. `group` itself is exercised by the
//! class-test port.
//!
//! Compared: the number and order of consensus features, each one's handles
//! (map index and unique ID, exactly, plus RT, m/z and intensity) and each
//! consensus feature's RT, m/z, intensity, charge and quality, and its peptide
//! identifications (map index and first-hit sequence, in the order the tool
//! sorts them by map index). The upstream
//! whitelist skips every line containing `id=`, which includes every handle
//! line, so this is stricter than the upstream test. The retained file stores
//! numbers with six significant digits: a stored value differs from the value
//! printed by at most a relative 5e-6, which is the tolerance here, plus
//! one `f32` rounding unit where the value is a `float`, plus an absolute
//! 1e-6 so that a value printed as `0` compares.
//!
//! Case 2 sets `distance_intensity:weight` to 0.5, so its result depends on
//! the normalising maximum intensity, which the source takes from stale
//! stored ranges after the first pairing.

#![cfg(all(feature = "featurexml", feature = "consensusxml", feature = "paramxml"))]

use openms::analysis::mapmatching::feature_grouping_algorithm_unlabeled::FeatureGroupingAlgorithmUnlabeled;
use openms::format::{consensusxml, featurexml, paramxml};

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

    // The tool links incrementally: the map with the most features is the
    // reference (the first of equal ones), every other map is added in input
    // order, and the raw result is sorted by m/z (stably).
    let maps: Vec<_> = (0..inputs)
        .map(|i| {
            featurexml::load(format!(
                "{dir}/feature_linker_unlabeled_{case}_input{}.featureXML",
                i + 1
            ))
            .unwrap()
        })
        .collect();
    let mut reference = 0;
    let mut max_count = 0;
    for (i, map) in maps.iter().enumerate() {
        if map.len() > max_count {
            max_count = map.len();
            reference = i;
        }
    }
    algorithm
        .set_reference(reference as u64, &maps[reference])
        .unwrap();
    for (i, map) in maps.iter().enumerate() {
        if i != reference {
            algorithm.add_to_group(i as u64, map).unwrap();
        }
    }
    let mut out = algorithm.result_map().clone();
    out.features
        .sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(std::cmp::Ordering::Equal));

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
        assert_eq!(
            identifications(produced),
            identifications(retained),
            "{at}: peptide identifications"
        );
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
