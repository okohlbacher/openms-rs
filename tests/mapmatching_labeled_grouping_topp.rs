// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Evidence tier 1 for `FeatureGroupingAlgorithmLabeled` and
//! `LabeledPairFinder`: the algorithm path of the upstream
//! `TOPP_FeatureLinkerLabeled_1` and `_2` tests, compared with the retained
//! C++ outputs.
//!
//! Upstream definition: `topp/CMakeLists.txt:2872-2877` at test-data revision
//! 0cb15f2 runs `FeatureLinkerLabeled -test -ini <case>_parameters.ini -in
//! <case>_input.featureXML` and compares with `<case>_output.consensusXML`
//! through `${DIFF} -whitelist "id=" "href="`. The tool
//! (`FeatureLinkerBase.cpp` with `labeled = true`, topp revision 174b576)
//! loads the input, clears the features' meta values (keeping adduct
//! annotations, which these inputs lack), sets column header 0 from the
//! input, copies it to header 1, labels them `light` and `heavy`, and calls
//! `FeatureGroupingAlgorithmLabeled::group` with the `algorithm:` section. It
//! then assigns fresh unique IDs and sorts each consensus feature's
//! identifications stably by `map_index`, neither of which changes the order
//! or the values compared here. The tool is not ported yet; this test runs
//! the same library calls.
//!
//! Case 1 uses the manual RT window; case 2 estimates it with the Gaussian
//! fit (`rt_estimate` true).
//!
//! Case 2's retained qualities are not those of the pinned C++: the Release
//! build at bc9cc12 (oracle `M2-labeled`, ibminode06) reproduces the retained
//! pairs, positions and intensities, but its qualities differ from the
//! retained ones by up to a relative 3.7e-4 (first pair: 0.633799374 against
//! the retained 0.633831564), which the upstream FuzzyDiff (`ratio` 1.01)
//! accepts. The retained file predates the pin's `GaussFitter`. Case 2's
//! qualities are therefore compared with that Release oracle, bit for bit
//! up to one `f32` unit, in `tests/mapmatching_labeled_grouping_oracle.rs`,
//! and not here.
//!
//! Compared: the number and order of consensus features, each one's handles
//! (map index and unique ID exactly, plus RT, m/z and intensity), its RT,
//! m/z, intensity, charge and quality, and its peptide identifications (map
//! index and first-hit sequence, in order); and the protein and unassigned
//! peptide identification counts. The upstream whitelist skips every line
//! containing `id=`, including every handle line, so this is stricter than
//! the upstream test. The retained file stores numbers with six significant
//! digits (tolerance: relative 5e-6, plus one `f32` rounding unit for `float`
//! values, plus an absolute 1e-6 so that a value printed as `0` compares).

#![cfg(all(feature = "featurexml", feature = "consensusxml", feature = "paramxml"))]

use openms::analysis::mapmatching::feature_grouping_algorithm_labeled::FeatureGroupingAlgorithmLabeled;
use openms::concept::UniqueIdGenerator;
use openms::format::{consensusxml, featurexml, paramxml};
use openms::kernel::features::{ColumnHeader, ConsensusMap};

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

fn run_case(case: usize, compare_quality: bool) {
    let dir = "tests/data";
    let ini = paramxml::load(format!(
        "{dir}/feature_linker_labeled_{case}_parameters.ini"
    ))
    .unwrap();
    let algorithm_param = ini.copy("FeatureLinkerLabeled:1:algorithm:", true).unwrap();
    let (algorithm, warnings) =
        FeatureGroupingAlgorithmLabeled::from_param_with_warnings(&algorithm_param).unwrap();
    assert!(warnings.is_empty(), "{warnings:?}");

    let mut map = featurexml::load(format!(
        "{dir}/feature_linker_labeled_{case}_input.featureXML"
    ))
    .unwrap();
    for feature in &mut map.features {
        feature.metadata.clear();
    }
    let mut out = ConsensusMap::new();
    let header = ColumnHeader {
        filename: "input.featureXML".into(),
        size: map.len(),
        unique_id: map.unique_id,
        ..ColumnHeader::default()
    };
    out.column_headers.insert(1, header.clone());
    out.column_headers.insert(0, header);
    out.column_headers.get_mut(&0).unwrap().label = "light".into();
    out.column_headers.get_mut(&1).unwrap().label = "heavy".into();

    // The tool's `-test` seed; the IDs themselves are not compared.
    let mut generator = UniqueIdGenerator::from_seed(19_991_231_235_959);
    algorithm
        .group_with_generator(std::slice::from_ref(&map), &mut out, &mut generator)
        .unwrap();
    let expected = consensusxml::load(format!(
        "{dir}/feature_linker_labeled_{case}_output.consensusXML"
    ))
    .unwrap();
    assert_eq!(
        out.len(),
        expected.len(),
        "case {case}: consensus feature count"
    );
    assert_eq!(
        out.protein_identifications.len(),
        expected.protein_identifications.len(),
        "case {case}: protein identifications"
    );
    assert_eq!(
        out.unassigned_peptide_identifications.len(),
        expected.unassigned_peptide_identifications.len(),
        "case {case}: unassigned peptide identifications"
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
        if compare_quality {
            close(
                f64::from(produced.quality),
                f64::from(retained.quality),
                float_rel,
                &format!("{at} quality"),
            );
        }
        assert_eq!(produced.charge, retained.charge, "{at}: charge");
        assert_eq!(
            identifications(produced),
            identifications(retained),
            "{at}: peptide identifications"
        );
    }
}

#[test]
fn topp_feature_linker_labeled_1() {
    run_case(1, true);
}

#[test]
fn topp_feature_linker_labeled_2() {
    // Qualities: see the module documentation.
    run_case(2, false);
}
