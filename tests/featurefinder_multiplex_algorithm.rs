// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `FeatureFinderMultiplexAlgorithm_test.cpp` at source revision
//! bc9cc12514c768385ce121d6ca4bb710fe1983c4
//! (`src/tests/class_tests/openms/source/`), on its inputs
//! `FeatureFinderMultiplex_1_input.mzML` and
//! `FeatureFinderMultiplex_1_parameters.ini` (class-test data at the same
//! revision).
//!
//! Evidence tier 3: the consensus size and the heavy-to-light ratio with its
//! absolute tolerance 0.2 are transcribed class-test literals.

#![cfg(all(feature = "mzml", feature = "paramxml"))]

use std::path::PathBuf;

use openms::analysis::featurefinder::multiplex::FeatureFinderMultiplexAlgorithm;
use openms::analysis::featurefinder::multiplex::algorithm::FeatureFinderMultiplexParameters;
use openms::concept::unique_id::UniqueIdGenerator;
use openms::format::{mzml, paramxml};

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join(name)
}

#[test]
fn feature_finder_multiplex_algorithm_class_test() {
    // FeatureFinderMultiplexAlgorithm()
    let defaults = FeatureFinderMultiplexAlgorithm::new();
    assert_eq!(
        defaults.parameters,
        FeatureFinderMultiplexParameters::default()
    );

    // run(): MzMLFile with addMSLevel(1)
    let mut exp = mzml::load(data("FeatureFinderMultiplex_1_input.mzML")).unwrap();
    exp.spectra.retain(|s| s.ms_level == 1);
    let param = paramxml::load(data("FeatureFinderMultiplex_1_parameters.ini")).unwrap();
    let mut param = param.copy("FeatureFinderMultiplex:1:", true).unwrap();
    let removed = [
        "in",
        "out",
        "out_multiplets",
        "log",
        "debug",
        "threads",
        "no_progress",
        "force",
        "test",
    ];
    for key in removed {
        param.remove(key).unwrap();
    }
    let mut algorithm = FeatureFinderMultiplexAlgorithm::from_param(&param).unwrap();
    algorithm.unique_ids = UniqueIdGenerator::from_seed(1);
    algorithm.run(&mut exp, true).unwrap();
    // The input was moved into the algorithm.
    assert!(exp.spectra.is_empty());

    let result = algorithm.consensus_map();
    assert_eq!(result.len(), 2);
    let handles = result.features[0].handles();
    let light = f64::from(handles[0].intensity);
    let heavy = f64::from(handles[1].intensity);
    // HEAVY:LIGHT is close to the expected 3:1 (TOLERANCE_ABSOLUTE(0.2)).
    assert!(
        (heavy / light - 3.0).abs() <= 0.2,
        "H/L = {}",
        heavy / light
    );
}

/// Tier 4: the parameter tree round-trips, and its sections carry the
/// source's descriptions and restrictions.
#[test]
fn feature_finder_multiplex_parameters_round_trip() {
    let defaults = FeatureFinderMultiplexAlgorithm::defaults().unwrap();
    assert_eq!(
        defaults
            .value("algorithm:labels")
            .unwrap()
            .as_str()
            .unwrap(),
        "[][Lys8,Arg10]"
    );
    assert_eq!(
        defaults
            .value("algorithm:rt_typical")
            .unwrap()
            .to_f64()
            .unwrap(),
        40.0
    );
    assert_eq!(
        defaults.value("labels:Lys8").unwrap().to_f64().unwrap(),
        8.0141988132
    );
    assert!(defaults.has_tag("labels:Lys8", "advanced").unwrap());
    assert_eq!(
        defaults.section_description("algorithm").unwrap(),
        "algorithmic parameters"
    );
    assert_eq!(
        defaults.section_description("labels").unwrap(),
        "mass shifts for all possible labels"
    );
    let parsed = FeatureFinderMultiplexAlgorithm::from_param(&defaults).unwrap();
    assert_eq!(
        parsed.parameters,
        FeatureFinderMultiplexParameters::default()
    );
    assert!(parsed.to_param().unwrap().source_equal(&defaults).unwrap());
}

/// Tier 4: an empty input is refused ("Error: No MS1 spectra in input
/// file.", source `Exception::FileEmpty`).
#[test]
fn feature_finder_multiplex_refuses_empty_input() {
    let mut algorithm = FeatureFinderMultiplexAlgorithm::new();
    let mut exp = openms::kernel::MSExperiment::default();
    let error = algorithm.run(&mut exp, false).unwrap_err();
    assert!(error.to_string().contains("No MS1 spectra"), "{error}");
}
