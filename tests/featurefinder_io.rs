// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The feature finder class tests that read upstream fixture files:
//! `PeakWidthEstimator_test.cpp` (`PeakPickerHiRes_orbitrap.mzML`, picked with
//! `ms_levels = [1]` and `signal_to_noise = 0`) and the consensus-map section of
//! `SeedListGenerator_test.cpp` (`ConsensusXMLFile_1.consensusXML`), both at
//! core `bc9cc12`. The two fixtures are the class-test data files, byte for
//! byte (sha256 `ca098bbd...df22` and `771a84b7...cd7d`), already carried for
//! the peak picker and consensusXML tests.
#![cfg(all(feature = "mzml", feature = "consensusxml"))]

use openms::analysis::featurefinder::peak_width_estimator::PeakWidthEstimator;
use openms::analysis::featurefinder::seed_list_generator::SeedListGenerator;
use openms::data_structures::DPosition2;
use openms::format::{consensusxml, mzml};
use openms::kernel::MSExperiment;
use openms::processing::peak_picking::{PeakBoundary, PeakPickerHiRes};

/// ClassTest `isRealSimilar` with its default tolerances (absolute 1e-5,
/// ratio 1 + 1e-5).
fn similar(actual: f64, expected: f64) {
    let abs_small = (actual - expected).abs() <= 1e-5;
    let mut ratio = actual / expected;
    if ratio < 1.0 {
        ratio = 1.0 / ratio;
    }
    assert!(
        abs_small || (ratio > 0.0 && ratio <= 1.0 + 1e-5),
        "{actual} is not similar to {expected}"
    );
}

fn orbitrap() -> MSExperiment {
    mzml::read_with_load_options(
        include_bytes!("data/peak_picking/PeakPickerHiRes_orbitrap.mzML").as_slice(),
        &mzml::LoadOptions::default(),
        &mzml::ReadOptions::default(),
    )
    .unwrap()
}

#[test]
fn peak_width_estimator() {
    let picker = PeakPickerHiRes {
        signal_to_noise: 0.0,
        ms_levels: vec![1],
        ..PeakPickerHiRes::default()
    };
    let picked = picker.pick_experiment(&orbitrap()).unwrap();
    // The source's pickExperiment appends boundaries for picked spectra only.
    let boundaries: Vec<Vec<PeakBoundary>> = picked
        .spectrum_boundaries
        .iter()
        .flatten()
        .cloned()
        .collect();
    let estimator = PeakWidthEstimator::new(&picked.experiment, &boundaries).unwrap();
    similar(estimator.peak_width(365.3).unwrap(), 0.00886469661896705);
    let estimator2 = PeakWidthEstimator::new(&picked.experiment, &boundaries).unwrap();
    similar(estimator2.peak_width(365.3).unwrap(), 0.00886469661896705);
    // outside m/z range
    similar(estimator2.peak_width(305.1).unwrap(), 0.00886699447290451);
    similar(estimator2.peak_width(405.1).unwrap(), 0.01184458329884600);
}

#[test]
fn peak_width_estimator_refuses_no_peaks() {
    assert!(PeakWidthEstimator::new(&MSExperiment::default(), &[]).is_err());
}

#[test]
fn seed_lists_from_consensus_map() {
    let consensus = consensusxml::read(
        include_bytes!("data/consensusxml/ConsensusXMLFile_1.consensusXML").as_slice(),
    )
    .unwrap();
    let seed_lists = SeedListGenerator::new()
        .generate_seed_lists(&consensus)
        .unwrap();
    assert_eq!(seed_lists.len(), 2);
    assert_eq!(seed_lists[&0].len(), 0);
    assert_eq!(seed_lists[&1].len(), 2);
    let pos = |rt: f64, mz: f64| DPosition2 {
        coordinates: [rt, mz],
    };
    assert_eq!(seed_lists[&1][0], pos(1273.27, 904.47));
    assert_eq!(seed_lists[&1][1], pos(1184.46, 953.368));
}
