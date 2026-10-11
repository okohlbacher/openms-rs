// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Evidence tier 1 for `ConsensusMapNormalizerAlgorithmThreshold`: the
//! algorithm path of the upstream `TOPP_ConsensusMapNormalizer_1` test,
//! compared with the retained C++ output.
//!
//! Upstream definition: `topp/CMakeLists.txt:337-339` at test-data revision
//! 0cb15f2 runs `ConsensusMapNormalizer -test` with default parameters on
//! `ConsensusMapNormalizer_input.consensusXML` and compares with
//! `ConsensusMapNormalizer_output.consensusXML` through `${DIFF}` (FuzzyDiff,
//! `ratio` 1.01, `absdiff` 0.01). The tool's defaults, at topp revision
//! 174b576 (`src/ConsensusMapNormalizer.cpp`), are `algorithm_type`
//! `robust_regression`, `ratio_threshold` 0.67 and empty filters, which run
//! `map.sortBySize()`, `computeCorrelation(map, 0.67, "", "")` and
//! `normalizeMaps(map, ratios)`. The tool itself is not ported yet; this test
//! runs the same three library calls.
//!
//! The retained file stores intensities with six significant digits, so a
//! stored value differs from the double it was printed from by at most half a
//! unit in the sixth digit, a relative 5e-6. That is the bound asserted here,
//! plus one `f32` rounding unit for the narrowing both sides perform. It is
//! far inside the FuzzyDiff ratio the upstream test uses.

#![cfg(feature = "consensusxml")]

use openms::analysis::mapmatching::consensus_map_normalizer_median::IdentificationFilter;
use openms::analysis::mapmatching::consensus_map_normalizer_threshold::ConsensusMapNormalizerAlgorithmThreshold as Threshold;
use openms::format::consensusxml;

const INPUT: &str = "tests/data/consensus_map_normalizer_input.consensusXML";
const EXPECTED: &str = "tests/data/consensus_map_normalizer_output.consensusXML";

#[test]
fn topp_consensus_map_normalizer_1_robust_regression() {
    let mut map = consensusxml::load(INPUT).unwrap();
    let expected = consensusxml::load(EXPECTED).unwrap();

    map.sort_by_size().unwrap();
    let filter = IdentificationFilter::new("", "").unwrap();
    let factors = Threshold::compute_correlation(&map, 0.67, &filter).unwrap();
    assert!(factors.enough_features);
    assert_eq!(factors.passing_features, map.len());
    Threshold::normalize_maps(&mut map, &factors.ratios).unwrap();

    assert_eq!(map.len(), expected.len());
    let mut compared = 0usize;
    let mut changed = 0usize;
    for (produced, retained) in map.features.iter().zip(&expected.features) {
        assert_eq!(produced.len(), retained.len());
        for (a, b) in produced.handles().iter().zip(retained.handles()) {
            assert_eq!(a.key(), b.key(), "consensus order differs after sortBySize");
            let (x, y) = (f64::from(a.intensity), f64::from(b.intensity));
            let bound = 5.0e-6 * y.abs() + f64::from(f32::EPSILON) * y.abs();
            assert!(
                (x - y).abs() <= bound,
                "handle {:?}: produced {x}, retained {y}",
                a.key()
            );
            compared += 1;
        }
    }
    // The fixture is only evidence if normalisation moved the intensities.
    let input: std::collections::BTreeMap<(u64, u64), f32> = consensusxml::load(INPUT)
        .unwrap()
        .features
        .iter()
        .flat_map(|cf| cf.handles().iter().map(|h| (h.key(), h.intensity)))
        .collect();
    for retained in expected.features.iter().flat_map(|cf| cf.handles()) {
        if input.get(&retained.key()) != Some(&retained.intensity) {
            changed += 1;
        }
    }
    assert!(compared > 0);
    assert!(changed > 0);
}
