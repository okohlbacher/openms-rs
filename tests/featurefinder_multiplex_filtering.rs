// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of the multiplex filtering class tests at source revision
//! bc9cc12514c768385ce121d6ca4bb710fe1983c4
//! (`src/tests/class_tests/openms/source/`): `MultiplexFiltering_test.cpp`
//! `MultiplexFilteringCentroided_test.cpp`, `MultiplexFilteringProfile_test.cpp`,
//! `MultiplexClustering_test.cpp` and `PeakWidthEstimator_test.cpp`, on their
//! inputs `MultiplexFiltering.mzML`, `MultiplexClustering.mzML` (byte-identical
//! to `FeatureFinderMultiplex_1_input.mzML`, sha256 `af1d7122...`, so only the
//! latter is kept) and `PeakPickerHiRes_orbitrap.mzML` (class-test data at the
//! same revision).
//!
//! Evidence tier 3: the expected result sizes are transcribed class-test
//! literals, which the C++ build asserts on the same input.

#![cfg(feature = "mzml")]

use std::path::PathBuf;

use openms::analysis::featurefinder::multiplex::{
    DeltaMass, LabelSet, MultiplexClustering, MultiplexDeltaMasses, MultiplexDistance,
    MultiplexFiltering, MultiplexFilteringCentroided, MultiplexFilteringProfile,
    MultiplexIsotopicPeakPattern,
};
use openms::analysis::featurefinder::peak_width_estimator::PeakWidthEstimator;
use openms::data_structures::dposition::DPosition;
use openms::format::mzml;
use openms::kernel::MSExperiment;
use openms::processing::peak_picking::PeakBoundary;
use openms::processing::peak_picking::{PeakPickerHiRes, PickingCompatibility};

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join(name)
}

/// The class tests' picking: `PeakPickerHiRes` on MS level 1 with
/// `signal_to_noise` 0. Returns the profile data, the picked data and the
/// spectrum boundaries, which the source appends for picked spectra only.
fn load_and_pick(path: PathBuf) -> (MSExperiment, MSExperiment, Vec<Vec<PeakBoundary>>) {
    let exp = mzml::load(path).unwrap();
    let picker = PeakPickerHiRes {
        ms_levels: vec![1],
        signal_to_noise: 0.0,
        // the C++ picker's behaviour on any input
        compatibility: PickingCompatibility::source(),
        ..PeakPickerHiRes::default()
    };
    let picked = picker.pick_experiment(&exp).unwrap();
    let boundaries = picked.spectrum_boundaries.into_iter().flatten().collect();
    (exp, picked.experiment, boundaries)
}

fn picked(name: &str) -> MSExperiment {
    load_and_pick(data(name)).1
}

/// `TEST_REAL_SIMILAR` at the ClassTest default relative tolerance 1e-5.
fn similar(actual: f64, expected: f64) {
    let tolerance = 1e-5 * expected.abs().max(actual.abs());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} is not similar to {expected}"
    );
}

/// The class tests' patterns: an Arg8 pair and an Arg8 double shift, for each
/// charge from `charge_max` down to `charge_min`.
fn patterns(
    charge_min: i32,
    charge_max: i32,
    isotopes_per_peptide_max: i32,
) -> Vec<MultiplexIsotopicPeakPattern> {
    let mut shifts1 = MultiplexDeltaMasses::new();
    shifts1
        .delta_masses
        .push(DeltaMass::with_label(0.0, "no_label"));
    shifts1
        .delta_masses
        .push(DeltaMass::with_label(8.0443702794, "Arg8"));
    let mut shifts2 = MultiplexDeltaMasses::new();
    shifts2
        .delta_masses
        .push(DeltaMass::with_label(0.0, "no_label"));
    let label_set: LabelSet = ["Arg8", "Arg8"].into_iter().collect();
    shifts2
        .delta_masses
        .push(DeltaMass::new(2.0 * 8.0443702794, label_set));
    let mut patterns = Vec::new();
    for c in (charge_min..=charge_max).rev() {
        patterns.push(
            MultiplexIsotopicPeakPattern::new(c, isotopes_per_peptide_max, shifts1.clone(), 0)
                .unwrap(),
        );
        patterns.push(
            MultiplexIsotopicPeakPattern::new(c, isotopes_per_peptide_max, shifts2.clone(), 1)
                .unwrap(),
        );
    }
    patterns
}

#[test]
fn multiplex_filtering_class_test() {
    let exp_picked = picked("MultiplexFiltering.mzML");
    let filtering = MultiplexFiltering::new(
        &exp_picked,
        patterns(1, 6, 6),
        3,
        6,
        10.0,
        2.0,
        40.0,
        true,
        0.8,
        0.75,
        0.75,
        "peptide",
    )
    .unwrap();
    // Constructor-only section; the intensity filter keeps every spectrum.
    assert_eq!(filtering.centroided_experiment().len(), exp_picked.len());
    assert_eq!(filtering.patterns().len(), 12);
    // Every peak starts white.
    assert!(filtering.blacklist().iter().flatten().all(|&m| m == -1));
    assert!(
        filtering
            .blacklist_experiment()
            .spectra
            .iter()
            .all(|s| s.peaks.is_empty())
    );
}

#[test]
fn multiplex_filtering_centroided_class_test() {
    let exp_picked = picked("MultiplexFiltering.mzML");
    let mut filtering = MultiplexFilteringCentroided::new(
        &exp_picked,
        patterns(1, 4, 6),
        3,
        6,
        10.0,
        3.0,
        40.0,
        true,
        0.8,
        0.75,
        0.75,
        "peptide",
    )
    .unwrap();
    let results = filtering.filter().unwrap();
    let sizes: Vec<usize> = results.iter().map(|r| r.len()).collect();
    assert_eq!(sizes, vec![0, 0, 0, 0, 4, 4, 4, 0]);
}

/// Tier 4: an unknown averagine type fails once the averagine filter first
/// runs ("Invalid averagine type."), as in the source.
#[test]
fn unknown_averagine_type_fails_when_used() {
    let exp_picked = picked("MultiplexFiltering.mzML");
    let mut filtering = MultiplexFilteringCentroided::new(
        &exp_picked,
        patterns(1, 4, 6),
        3,
        6,
        10.0,
        3.0,
        40.0,
        true,
        0.8,
        0.75,
        0.75,
        "protein",
    )
    .unwrap();
    let error = filtering.filter().unwrap_err();
    assert!(
        error.to_string().contains("Invalid averagine type."),
        "{error}"
    );
}

#[test]
fn multiplex_filtering_profile_class_test() {
    let (exp, exp_picked, boundaries) = load_and_pick(data("MultiplexFiltering.mzML"));
    let mut filtering = MultiplexFilteringProfile::new(
        &exp,
        &exp_picked,
        &boundaries,
        patterns(1, 4, 6),
        3,
        6,
        10.0,
        3.0,
        40.0,
        true,
        0.8,
        0.75,
        0.75,
        "peptide",
    )
    .unwrap();
    let results = filtering.filter().unwrap();
    let sizes: Vec<usize> = results.iter().map(|r| r.len()).collect();
    assert_eq!(sizes, vec![0, 0, 0, 0, 4, 5, 4, 0]);
}

#[test]
fn multiplex_clustering_class_test() {
    let (exp, exp_picked, boundaries) = load_and_pick(data("FeatureFinderMultiplex_1_input.mzML"));
    let mut filtering = MultiplexFilteringProfile::new(
        &exp,
        &exp_picked,
        &boundaries,
        patterns(1, 4, 6),
        3,
        6,
        10.0,
        3.0,
        40.0,
        true,
        0.8,
        0.75,
        0.75,
        "peptide",
    )
    .unwrap();
    let filter_results = filtering.filter().unwrap();
    let clustering =
        MultiplexClustering::from_profile(&exp, &exp_picked, &boundaries, 90.0).unwrap();
    let cluster_results = clustering.cluster(&filter_results).unwrap();
    let sizes: Vec<usize> = cluster_results.iter().map(|r| r.len()).collect();
    assert_eq!(sizes, vec![0, 0, 0, 0, 2, 0, 0, 0]);
}

/// Tier 4: the scaled Euclidean distance of the header,
/// `sqrt(dx^2 + s^2 dy^2)`, and its default scaling 1.
#[test]
fn multiplex_distance_is_scaled_euclidean() {
    let p1 = DPosition::xy(1.0, 10.0);
    let p2 = DPosition::xy(4.0, 12.0);
    assert_eq!(
        MultiplexDistance::default().distance(&p1, &p2),
        13.0_f64.sqrt()
    );
    assert_eq!(
        MultiplexDistance::new(2.0).distance(&p1, &p2),
        25.0_f64.sqrt()
    );
}

#[test]
fn peak_width_estimator_class_test() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/peak_picking/PeakPickerHiRes_orbitrap.mzML");
    let (_, exp_picked, boundaries) = load_and_pick(path);
    let estimator = PeakWidthEstimator::new(&exp_picked, &boundaries).unwrap();
    similar(estimator.peak_width(365.3).unwrap(), 0.00886469661896705);
    similar(estimator.peak_width(305.1).unwrap(), 0.00886699447290451); // outside m/z range
    similar(estimator.peak_width(405.1).unwrap(), 0.01184458329884600); // outside m/z range
}
