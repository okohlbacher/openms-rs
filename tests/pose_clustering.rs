// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The pose-clustering superimposers against the upstream class tests at core
//! pin `bc9cc12` (evidence tier 3): `PoseClusteringShiftSuperimposer_test.cpp`
//! and `PoseClusteringAffineSuperimposer_test.cpp`, plus independently derived
//! checks of the parameter surface and the refusals.

use openms::Error;
use openms::analysis::mapmatching::alignment_pose_clustering::MapAlignmentAlgorithmPoseClustering;
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

/// The one file in `directory` whose name starts with `prefix` followed by a
/// serial number.
fn dump(directory: &std::path::Path, prefix: &str) -> String {
    let names: Vec<_> = std::fs::read_dir(directory)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .filter(|name| {
            name.strip_prefix(prefix)
                .is_some_and(|serial| serial.parse::<u32>().is_ok())
        })
        .collect();
    assert_eq!(names.len(), 1, "{prefix}: {names:?}");
    std::fs::read_to_string(directory.join(&names[0])).unwrap()
}

#[test]
fn shift_superimposer_writes_the_source_debug_dumps() {
    // Derived from the source's dump statements for the class-test maps:
    // one pair (i = 0, k = 0) of similarity 0.9 * 0.9, a 21-bucket histogram
    // offset to -30 in steps of 3, its only votes in buckets 3 and 4, and
    // three identical cut-off rounds at index mean 3.2 (key -20.4) and
    // standard deviation 0.4 (1.2 in key units) over buckets 2..5.
    let directory = std::env::temp_dir().join(format!("m1-shift-dumps-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let superimposer = PoseClusteringShiftSuperimposer {
        dump_buckets: directory.join("buckets").to_str().unwrap().to_owned(),
        dump_pairs: directory.join("pairs").to_str().unwrap().to_owned(),
        ..Default::default()
    };
    let model = consensus(&[(1.0, 1.0), (5.0, 5.0)]);
    let scene = consensus(&[(21.4, 1.02), (25.4, 5.02)]);
    let transformation = superimposer.run(&model, &scene).unwrap();
    check_linear(&transformation, 1.0, -20.4);
    let pairs = dump(&directory, "pairs");
    assert_eq!(pairs, "# i k\n0 1 1 0 21.4 1.02 0.81 \n");
    let buckets = dump(&directory, "buckets_");
    let lines: Vec<_> = buckets.lines().collect();
    assert!(lines[0].starts_with("# shift hash table buckets dump ( scale, height ) : "));
    assert_eq!(lines[1], "# unfiltered hash data");
    assert_eq!(lines[2], "0\t0\t-30\t0");
    assert_eq!(lines[5], "0\t3\t-21\t0.648");
    assert_eq!(lines[23], "");
    assert_eq!(lines[24], "# tophat filtered hash data");
    assert_eq!(lines[47], "# after freq_cutoff, which is: 0");
    assert_eq!(
        lines[70],
        "# loop: 0  mean: -20.4  stdev: 1.2  (mean-stdev): -21.6  (mean+stdev): -19.2  data_range_begin: 2  data_range_end: 5"
    );
    assert_eq!(lines[73], "# EOF");
    assert_eq!(lines.len(), 74);
    std::fs::remove_dir_all(&directory).unwrap();
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
fn superimposers_refuse_empty_maps_unwritable_dumps_and_zero_buckets() {
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
    // A dump into a directory that does not exist is an I/O error, where the
    // source's stream fails silently.
    let dumping = PoseClusteringShiftSuperimposer {
        dump_buckets: "/nonexistent-directory-for-openms-tests/buckets".into(),
        ..Default::default()
    };
    assert!(matches!(dumping.run(&model, &model), Err(Error::Io(_))));
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

/// A consensus map of single-handle features in map `index`, as
/// `MapConversion::convert` produces them.
fn converted(index: u64, points: &[(f64, f64)]) -> ConsensusMap {
    ConsensusMap::from_features(
        points
            .iter()
            .map(|&(rt, mz)| {
                ConsensusFeature::from_feature(index, &BaseFeature::new(rt, mz, 100.0)).unwrap()
            })
            .collect(),
    )
}

#[test]
fn pose_clustering_alignment_inverts_the_superimposer_for_the_pairs() {
    // Independently derived from the source's align(ConsensusMap): with the
    // class test's two-point maps the superimposer finds x - 0.4; the scene
    // handed to the pair finder is moved onto the reference, features and
    // first handles alike; each pair's scene RT is taken back through the
    // inverse (x + 0.4) before the final linear fit, which is therefore
    // x - 0.4 again. The pair finder here pairs features by position, in
    // place of StablePairFinder.
    let mut aligner = MapAlignmentAlgorithmPoseClustering::default();
    let mut params = MapAlignmentAlgorithmPoseClustering::defaults().unwrap();
    params
        .set_value(
            "superimposer:scaling_bucket_size",
            ParamValue::Float(0.01),
            "",
            &[],
        )
        .unwrap();
    params
        .set_value(
            "superimposer:shift_bucket_size",
            ParamValue::Float(0.1),
            "",
            &[],
        )
        .unwrap();
    assert!(aligner.set_parameters(&params).unwrap().is_empty());
    aligner.set_reference_consensus(&converted(0, &[(1.0, 1.0), (5.0, 5.0)]));
    let scene = converted(1, &[(1.4, 1.02), (5.4, 5.02)]);
    let mut seen = Vec::new();
    let trafo = aligner
        .align_consensus(
            &scene,
            |maps| {
                let mut pairs = Vec::new();
                for (reference, moved) in maps[0].features.iter().zip(&maps[1].features) {
                    seen.push((moved.rt, moved.handles()[0].rt));
                    let mut pair = reference.clone();
                    pair.merge(moved)?;
                    pairs.push(pair);
                }
                Ok(ConsensusMap::from_features(pairs))
            },
            &mut ProgressReporter::silent(),
        )
        .unwrap();
    for ((feature, handle), expected) in seen.iter().zip([1.0, 5.0]) {
        similar(*feature, expected);
        similar(*handle, expected);
    }
    assert_eq!(trafo.model_type(), "linear");
    assert_eq!(trafo.data_points().len(), 2);
    for point in trafo.data_points() {
        similar(point.x - point.y, 0.4);
    }
    similar(trafo.apply(10.0).unwrap(), 9.6);
    // No pairs: the linear fit has no data, as in the source.
    let unpaired = aligner.align_consensus(
        &scene,
        |_| Ok(ConsensusMap::new()),
        &mut ProgressReporter::silent(),
    );
    assert!(matches!(unpaired, Err(Error::InvalidValue(_))));
    // Parameter layout: the superimposer section and the peak count.
    let defaults = MapAlignmentAlgorithmPoseClustering::defaults().unwrap();
    assert_eq!(
        defaults.value("max_num_peaks_considered").unwrap(),
        &ParamValue::Integer(1000)
    );
    assert_eq!(
        defaults.value("superimposer:mz_pair_max_distance").unwrap(),
        &ParamValue::Float(0.5)
    );
}
