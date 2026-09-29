// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of the multiplex filtering class tests at source revision
//! bc9cc12514c768385ce121d6ca4bb710fe1983c4
//! (`src/tests/class_tests/openms/source/`): `MultiplexFiltering_test.cpp`
//! and `MultiplexFilteringCentroided_test.cpp`, on their input
//! `MultiplexFiltering.mzML` (class-test data at the same revision).
//!
//! Evidence tier 3: the expected result sizes are transcribed class-test
//! literals, which the C++ build asserts on the same input.

#![cfg(feature = "mzml")]

use std::path::PathBuf;

use openms::analysis::featurefinder::multiplex::{
    DeltaMass, LabelSet, MultiplexDeltaMasses, MultiplexFiltering, MultiplexFilteringCentroided,
    MultiplexIsotopicPeakPattern,
};
use openms::format::mzml;
use openms::kernel::MSExperiment;
use openms::processing::peak_picking::PeakPickerHiRes;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join(name)
}

/// The class tests' picking: `PeakPickerHiRes` on MS level 1 with
/// `signal_to_noise` 0.
fn picked(name: &str) -> MSExperiment {
    let exp = mzml::load(data(name)).unwrap();
    let picker = PeakPickerHiRes {
        ms_levels: vec![1],
        signal_to_noise: 0.0,
        ..PeakPickerHiRes::default()
    };
    picker.pick_experiment(&exp).unwrap().experiment
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
