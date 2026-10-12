// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `MapAlignmentAlgorithmPoseClustering_test.cpp` at core pin `bc9cc12`
//! (`src/tests/class_tests/openms/source/`), run end to end with the
//! `StablePairFinder` the algorithm composes (evidence tier 3: the literals
//! are the class test's, 307 data points, slope 1.01164 and intercept
//! -32.0912).
//!
//! The input pair is the class test's data (`MapAlignmentAlgorithmPoseClustering_in1.mzML.gz`
//! and `_in2.mzML.gz`, copied byte for byte into `tests/data/pose_clustering/`;
//! hashes in `tests/data/map_alignment_provenance.json`). The class test's
//! other two sections (`align(FeatureMap)` and `align(ConsensusMap)`) are
//! `NOT_TESTABLE` upstream ("tested in TEST/TOPP").
#![cfg(feature = "mzml")]

use openms::analysis::mapmatching::alignment_pose_clustering::MapAlignmentAlgorithmPoseClustering;
use openms::concept::UniqueIdGenerator;
use openms::concept::progress_logger::ProgressReporter;
use openms::format::file_handler::FileHandler;
use openms::format::file_types::FileType;
use openms::format::peak_options::PeakFileOptions;

/// ClassTest `isRealSimilar`: absolute difference at most 1e-5, or ratio at
/// most 1 + 1e-5.
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

fn load(name: &str) -> openms::kernel::MSExperiment {
    FileHandler::load_experiment_with_options(
        format!("tests/data/pose_clustering/{name}"),
        &[FileType::MzMl],
        &PeakFileOptions::default(),
    )
    .unwrap()
}

#[test]
fn align_peak_map_class_test() {
    let reference = load("MapAlignmentAlgorithmPoseClustering_in1.mzML.gz");
    let scene = load("MapAlignmentAlgorithmPoseClustering_in2.mzML.gz");
    let mut generator = UniqueIdGenerator::default();

    let mut aligner = MapAlignmentAlgorithmPoseClustering::default();
    aligner
        .set_reference_peaks(&reference, &mut generator)
        .unwrap();
    let trafo = aligner
        .align_peaks(&scene, &mut generator, &mut ProgressReporter::silent())
        .unwrap();

    assert_eq!(trafo.model_type(), "linear");
    assert_eq!(trafo.data_points().len(), 307);
    let parameters = trafo.model_parameters();
    similar(
        parameters.value("slope").unwrap().to_f64().unwrap(),
        1.01164,
    );
    similar(
        parameters.value("intercept").unwrap().to_f64().unwrap(),
        -32.0912,
    );
}
