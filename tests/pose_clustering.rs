// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The pose-clustering superimposers against the upstream class tests at core
//! pin `bc9cc12` (evidence tier 3): `PoseClusteringShiftSuperimposer_test.cpp`
//! and `PoseClusteringAffineSuperimposer_test.cpp`, plus independently derived
//! checks of the parameter surface and the refusals.

use openms::Error;
use openms::analysis::mapmatching::pose_clustering_affine::PoseClusteringAffineSuperimposer;
use openms::analysis::mapmatching::pose_clustering_shift::PoseClusteringShiftSuperimposer;
use openms::analysis::mapmatching::superimposer::BaseSuperimposer;
use openms::analysis::transformations::TransformationDescription;
use openms::concept::progress_logger::ProgressReporter;
use openms::kernel::{BaseFeature, ConsensusFeature, ConsensusMap, Peak2D};
use openms::param::{Param, ParamValue};

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

/// `ConsensusFeature(feature)` of features at `(rt, mz)` with intensity 100.
fn consensus(points: &[(f64, f64)]) -> ConsensusMap {
    ConsensusMap::from_features(
        points
            .iter()
            .map(|&(rt, mz)| ConsensusFeature::from(BaseFeature::new(rt, mz, 100.0)))
            .collect(),
    )
}

fn check_linear(transformation: &TransformationDescription, slope: f64, intercept: f64) {
    assert_eq!(transformation.model_type(), "linear");
    let params = transformation.model_parameters();
    assert_eq!(params.size(), 2);
    similar(params.value("slope").unwrap().to_f64().unwrap(), slope);
    similar(
        params.value("intercept").unwrap().to_f64().unwrap(),
        intercept,
    );
}

fn param(values: &[(&str, ParamValue)]) -> Param {
    let mut p = Param::new();
    for (key, value) in values {
        p.set_value(key, value.clone(), "", &[]).unwrap();
    }
    p
}

#[test]
fn shift_superimposer_finds_the_upstream_shift() {
    // PoseClusteringShiftSuperimposer_test.cpp, section `run`.
    let model = consensus(&[(1.0, 1.0), (5.0, 5.0)]);
    let scene = consensus(&[(21.4, 1.02), (25.4, 5.02)]);
    let transformation = PoseClusteringShiftSuperimposer::default()
        .run(&model, &scene)
        .unwrap();
    check_linear(&transformation, 1.0, -20.4);
}

fn affine_parameters(extra: &[(&str, ParamValue)]) -> PoseClusteringAffineSuperimposer {
    let mut values = vec![
        ("scaling_bucket_size", ParamValue::Float(0.01)),
        ("shift_bucket_size", ParamValue::Float(0.1)),
    ];
    values.extend(extra.iter().cloned());
    PoseClusteringAffineSuperimposer::from_param(&param(&values)).unwrap()
}

fn peaks(rt: &[f64], mz: &[f64], intensity: &[f32]) -> Vec<Peak2D> {
    rt.iter()
        .zip(mz)
        .zip(intensity)
        .map(|((&rt, &mz), &intensity)| Peak2D::new(rt, mz, intensity))
        .collect()
}

#[test]
fn affine_superimposer_on_consensus_maps_and_points() {
    // PoseClusteringAffineSuperimposer_test.cpp, sections `run` (both
    // overloads) and the first `[EXTRA]` section.
    let superimposer = affine_parameters(&[]);
    let model = consensus(&[(1.0, 1.0), (5.0, 5.0)]);
    let scene = consensus(&[(1.4, 1.02), (5.4, 5.02)]);
    check_linear(&superimposer.run(&model, &scene).unwrap(), 1.0, -0.4);
    let model = peaks(&[1.0, 5.0], &[1.0, 5.0], &[100.0, 100.0]);
    let scene = peaks(&[1.4, 5.4], &[1.02, 5.02], &[100.0, 100.0]);
    let result = superimposer
        .superimpose_points(&model, &scene, &mut ProgressReporter::silent())
        .unwrap();
    check_linear(&result.transformation, 1.0, -0.4);
}

#[test]
fn affine_superimposer_selects_the_most_intense_points() {
    // PoseClusteringAffineSuperimposer_test.cpp, the second `[EXTRA]`
    // section: five points per map, reversed so that they are not sorted.
    let mut model = peaks(
        &[1.0, 5.0, 1.3, 2.2, 5.2],
        &[1.0, 5.0, 800.0, 900.0, 5.0],
        &[100.0, 100.0, 41.0, 20.0, 50.0],
    );
    let mut scene = peaks(
        &[1.4, 5.4, 4.4, 4.4, 5.8],
        &[1.02, 5.02, 800.0, 900.0, 5.02],
        &[100.0, 100.0, 40.0, 20.0, 50.0],
    );
    model.reverse();
    scene.reverse();
    let silent = || ProgressReporter::silent();
    let two = affine_parameters(&[("num_used_points", ParamValue::Integer(2))]);
    check_linear(
        &two.superimpose_points(&model, &scene, &mut silent())
            .unwrap()
            .transformation,
        1.0,
        -0.4,
    );
    let three = affine_parameters(&[("num_used_points", ParamValue::Integer(3))]);
    check_linear(
        &three
            .superimpose_points(&model, &scene, &mut silent())
            .unwrap()
            .transformation,
        0.977273,
        -0.368182,
    );
    // "What happens if we set the wrong parameters?": the source gets the
    // wrong results, and warns.
    let wrong = affine_parameters(&[
        ("num_used_points", ParamValue::Integer(3)),
        ("max_shift", ParamValue::Float(0.2)),
        ("max_scaling", ParamValue::Float(1.001)),
    ]);
    let result = wrong
        .superimpose_points(&model, &scene, &mut silent())
        .unwrap();
    check_linear(&result.transformation, 1.0, -0.4);
    assert_eq!(result.warnings.len(), 4, "{:?}", result.warnings);
    assert!(result.warnings[0].starts_with("WARNING: your map likely has a scaling around "));
    assert!(result.warnings[0].ends_with("maximal scaling of 1.001"));
    assert_eq!(
        result.warnings[3],
        "It is strongly advised to adjust your max_shift factor"
    );
}

#[test]
fn superimposer_parameters_round_trip_with_source_defaults() {
    let shift = PoseClusteringShiftSuperimposer::defaults().unwrap();
    assert_eq!(
        shift.value("mz_pair_max_distance").unwrap(),
        &ParamValue::Float(0.5)
    );
    assert_eq!(
        shift.value("num_used_points").unwrap(),
        &ParamValue::Integer(2000)
    );
    assert_eq!(shift.entry("num_used_points").unwrap().min_int, -1);
    assert!(shift.has_tag("max_shift", "advanced").unwrap());
    assert_eq!(
        PoseClusteringShiftSuperimposer::from_param(&shift).unwrap(),
        PoseClusteringShiftSuperimposer::default()
    );
    let affine = PoseClusteringAffineSuperimposer::defaults().unwrap();
    assert_eq!(
        affine.value("scaling_bucket_size").unwrap(),
        &ParamValue::Float(0.005)
    );
    assert_eq!(affine.entry("max_scaling").unwrap().min_float, 1.0);
    assert_eq!(
        affine.entry("rt_pair_distance_fraction").unwrap().max_float,
        1.0
    );
    assert_eq!(
        PoseClusteringAffineSuperimposer::from_param(&affine).unwrap(),
        PoseClusteringAffineSuperimposer::default()
    );
    // Unknown names are warnings, restriction violations errors.
    let (_, warnings) = PoseClusteringAffineSuperimposer::from_param_with_warnings(&param(&[(
        "no_such_parameter",
        ParamValue::Integer(1),
    )]))
    .unwrap();
    assert_eq!(warnings.len(), 1, "{warnings:?}");
    assert!(
        PoseClusteringAffineSuperimposer::from_param(&param(&[(
            "max_scaling",
            ParamValue::Float(0.5)
        )]))
        .is_err()
    );
}

#[test]
fn superimposers_refuse_empty_maps_debug_dumps_and_zero_buckets() {
    let model = consensus(&[(1.0, 1.0), (5.0, 5.0)]);
    let empty = ConsensusMap::new();
    assert!(matches!(
        PoseClusteringAffineSuperimposer::default().run(&model, &empty),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        PoseClusteringShiftSuperimposer::default().run(&empty, &model),
        Err(Error::InvalidValue(_))
    ));
    let dumping = PoseClusteringShiftSuperimposer {
        dump_buckets: "buckets".into(),
        ..Default::default()
    };
    assert!(matches!(
        dumping.run(&model, &model),
        Err(Error::Unsupported(_))
    ));
    let zero = PoseClusteringAffineSuperimposer {
        shift_bucket_size: 0.0,
        ..Default::default()
    };
    assert!(zero.run(&model, &model).is_err());
    // Two zero intensities make the pair similarity 0/0.
    let silent = ConsensusMap::from_features(vec![
        ConsensusFeature::from(BaseFeature::new(1.0, 1.0, 0.0)),
        ConsensusFeature::from(BaseFeature::new(5.0, 5.0, 0.0)),
    ]);
    assert!(
        PoseClusteringShiftSuperimposer::default()
            .run(&silent, &silent)
            .is_err()
    );
}
