// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Alignment precision and recall against the upstream class tests at core
//! pin `bc9cc12` (evidence tier 3): `MapAlignmentEvaluationAlgorithm_test.cpp`,
//! `MapAlignmentEvaluationAlgorithmPrecision_test.cpp` and
//! `MapAlignmentEvaluationAlgorithmRecall_test.cpp`, with their two
//! consensusXML inputs copied byte for byte (hashes in
//! `tests/data/map_alignment_provenance.json`).

#![cfg(feature = "consensusxml")]

use openms::analysis::mapmatching::alignment_evaluation::{
    MapAlignmentEvaluationAlgorithm, MapAlignmentEvaluationPrecision, MapAlignmentEvaluationRecall,
    is_same_handle,
};
use openms::format::consensusxml;
use openms::kernel::{BaseFeature, ConsensusMap, FeatureHandle};
use openms::{Error, Result};

fn maps() -> (ConsensusMap, ConsensusMap) {
    let input = consensusxml::read(
        include_bytes!("data/map_alignment_evaluation_in.consensusXML").as_slice(),
    )
    .unwrap();
    let truth = consensusxml::read(
        include_bytes!("data/map_alignment_evaluation_gt.consensusXML").as_slice(),
    )
    .unwrap();
    (input, truth)
}

/// ClassTest `isRealSimilar`: absolute difference at most 1e-5, or ratio at
/// most 1 + 1e-5 (`src/testframework/source/CONCEPT/ClassTest.cpp:364` at the
/// pin).
fn similar(actual: f64, expected: f64) {
    let ratio = if actual.abs() > expected.abs() {
        actual / expected
    } else {
        expected / actual
    };
    assert!(
        (actual - expected).abs() <= 1e-5 || (ratio - 1.0).abs() <= 1e-5,
        "{actual} !~ {expected}"
    );
}

#[test]
fn precision_of_the_upstream_maps() {
    let (input, truth) = maps();
    let out = MapAlignmentEvaluationPrecision
        .evaluate(&input, &truth, 0.1, 0.1, 100.0, true)
        .unwrap();
    similar(out, 0.757143);
}

#[test]
fn recall_of_the_upstream_maps_and_its_degenerate_inputs() {
    let (input, truth) = maps();
    let recall = MapAlignmentEvaluationRecall;
    similar(
        recall
            .evaluate(&input, &truth, 0.1, 0.1, 100.0, true)
            .unwrap(),
        0.5,
    );
    // An empty tool map recalls nothing; an empty ground truth is refused.
    let empty = ConsensusMap::new();
    similar(
        recall
            .evaluate(&empty, &truth, 0.1, 0.1, 100.0, true)
            .unwrap(),
        0.0,
    );
    assert!(matches!(
        recall.evaluate(&input, &empty, 0.1, 0.1, 100.0, true),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn precision_of_an_unusable_ground_truth_is_the_source_nan() {
    // Not in the class test: the source computes 1.0 / 0 * 0 without a guard.
    let (input, _) = maps();
    let out = MapAlignmentEvaluationPrecision
        .evaluate(&input, &ConsensusMap::new(), 0.1, 0.1, 100.0, true)
        .unwrap();
    assert!(out.is_nan());
}

#[test]
fn same_handle_by_position_intensity_and_charge() {
    // MapAlignmentEvaluationAlgorithm_test.cpp, section `isSameHandle`.
    let mut first = BaseFeature::new(100.0, 555.0, 200.0);
    first.charge = 3;
    first.unique_id = 1;
    let mut second = BaseFeature::new(101.0, 556.0, 1199.0);
    second.charge = 4;
    second.unique_id = 2;
    let a = FeatureHandle::new(0, &first);
    let b = FeatureHandle::new(0, &second);
    assert!(is_same_handle(&a, &b, 2.0, 1.5, 1000.0, false));
    assert!(!is_same_handle(&a, &b, 2.0, 1.5, 1000.0, true));
    second.charge = 3;
    let b2 = FeatureHandle::new(0, &second);
    assert!(is_same_handle(&a, &b2, 2.0, 1.5, 1000.0, false));
    assert!(is_same_handle(&a, &b2, 2.0, 1.5, 1000.0, true));
    // The trait's provided method is the same test.
    assert!(MapAlignmentEvaluationPrecision.is_same_handle(&a, &b2, 2.0, 1.5, 1000.0, true));
}

/// The class test's own subclass, whose `evaluate` returns 1.5 (section
/// `evaluate`): the interface is implementable outside the crate.
struct Constant;
impl MapAlignmentEvaluationAlgorithm for Constant {
    fn evaluate(
        &self,
        _: &ConsensusMap,
        _: &ConsensusMap,
        _: f64,
        _: f64,
        _: f32,
        _: bool,
    ) -> Result<f64> {
        Ok(1.5)
    }
}

#[test]
fn the_interface_is_implementable() {
    let empty = ConsensusMap::new();
    assert_eq!(
        Constant
            .evaluate(&empty, &empty, 0.0, 0.0, 0.0, false)
            .unwrap(),
        1.5
    );
}
