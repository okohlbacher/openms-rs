// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Every `START_SECTION` of `PeakTypeEstimator_test.cpp` at Core SDK `bc9cc12`,
//! the executed product-sdk oracle for `PeakTypeEstimator::estimateType`, and
//! the native boundaries of `src/format/peak_type_estimator.rs`.
//!
//! Evidence labels:
//!
//! - Expectations citing a `PeakTypeEstimator_test.cpp` line are transcribed
//!   class-test literals (tier 3, source review).
//! - `tests/data/peak_type_estimator/oracle_estimates.tsv` is oracle-generated
//!   (tier 1 executed differential): the unmodified product-sdk `libOpenMS`
//!   (Debug, core 4fdec46, identical to `bc9cc12` for every file involved) run
//!   by `../oracle/pte-faims-helper/driver.cpp`. The manifest
//!   `tests/data/peak_type_estimator_provenance.json` records the hashes.
//! - The finiteness refusal and the ceilings are native (tier 4); the source
//!   has no analogue.

use openms::format::peak_type_estimator::PeakTypeEstimator;
use openms::kernel::{Peak1D, SpectrumType, SpectrumTypeQueryLimits};
use openms::processing::peak_picking::estimate_spectrum_type;
use openms::{Error, MSSpectrum};
use std::io::Cursor;

const ORACLE: &str = include_str!("data/peak_type_estimator/oracle_estimates.tsv");
const RAW: &str = include_str!("data/spectrum_type/PeakTypeEstimator_raw.dta");
const RAW_TOF: &str = include_str!("data/spectrum_type/PeakTypeEstimator_rawTOF.dta");
const PEAK: &str = include_str!("data/spectrum_type/PeakTypeEstimator_peak.dta");

/// The oracle records of one kind, split into tab-separated fields.
fn records(kind: &str) -> Vec<Vec<&'static str>> {
    ORACLE
        .lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|fields| fields[0] == kind)
        .collect()
}

/// Parses the C `printf("%a")` spelling the oracle writes (`0x1.9p+6`,
/// `-0x0p+0`, `inf`, `nan`). Every mantissa has at most 53 significant bits,
/// so the conversion is exact.
fn hex_f64(text: &str) -> f64 {
    let (negative, body) = match text.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, text),
    };
    let magnitude = match body {
        "inf" => f64::INFINITY,
        "nan" => f64::NAN,
        _ => {
            let body = body.strip_prefix("0x").expect("hexadecimal float");
            let (mantissa, exponent) = body.split_once('p').expect("binary exponent");
            let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
            let mut bits: u64 = 0;
            for digit in whole.chars().chain(fraction.chars()) {
                bits = bits * 16 + u64::from(digit.to_digit(16).expect("hexadecimal digit"));
            }
            let shift = exponent.parse::<i32>().expect("exponent") - 4 * fraction.len() as i32;
            let mut value = bits as f64;
            for _ in 0..shift.unsigned_abs() {
                value = if shift > 0 { value * 2.0 } else { value / 2.0 };
            }
            value
        }
    };
    if negative { -magnitude } else { magnitude }
}

fn spectrum_type(name: &str) -> SpectrumType {
    match name {
        "Unknown" => SpectrumType::Unknown,
        "Centroid" => SpectrumType::Centroid,
        "Profile" => SpectrumType::Profile,
        other => panic!("unexpected spectrum type {other}"),
    }
}

/// The peaks of one `pte_synthetic` oracle record, exactly as the driver built
/// them (the intensities are `float` values printed through `double`).
fn synthetic(name: &str) -> Vec<Peak1D> {
    let row = records("pte_synthetic")
        .into_iter()
        .find(|row| row[1] == name)
        .expect("synthetic oracle case");
    let mz = row[4].split(',').map(hex_f64);
    let intensity = row[5].split(',').map(|text| hex_f64(text) as f32);
    mz.zip(intensity)
        .map(|(mz, intensity)| Peak1D::new(mz, intensity))
        .collect()
}

#[test]
fn extra_constructor_and_destructor_sections() {
    // PeakTypeEstimator_test.cpp:33-40 allocates and deletes an instance of the
    // stateless class; the Rust type is a zero-sized unit struct.
    let constructed = PeakTypeEstimator;
    assert_eq!(constructed, PeakTypeEstimator);
    assert_eq!(std::mem::size_of::<PeakTypeEstimator>(), 0);
    assert!(!std::mem::needs_drop::<PeakTypeEstimator>());
    assert_eq!(PeakTypeEstimator::MIN_PEAKS, 5);
}

#[test]
fn estimate_type_section_classifies_the_three_upstream_spectra() {
    // PeakTypeEstimator_test.cpp:42-57: raw data with zeros, TOF raw data
    // without zeros, peak data, and too few data points after spec.resize(4).
    for (text, expected) in [
        (RAW, SpectrumType::Profile),
        (RAW_TOF, SpectrumType::Profile),
        (PEAK, SpectrumType::Centroid),
    ] {
        let mut spectrum = openms::format::dta::read(Cursor::new(text)).unwrap();
        assert_eq!(
            PeakTypeEstimator::estimate_type(&spectrum.peaks).unwrap(),
            expected
        );
        spectrum.peaks.truncate(4);
        assert_eq!(
            PeakTypeEstimator::estimate_type(&spectrum.peaks).unwrap(),
            SpectrumType::Unknown
        );
    }
}

#[test]
fn the_upstream_spectra_match_the_executed_oracle_record_for_record() {
    let full = records("pte_dta");
    let first4 = records("pte_dta_first4");
    assert_eq!(full.len(), 3);
    assert_eq!(first4.len(), 3);
    for ((text, row), short) in [RAW, RAW_TOF, PEAK].into_iter().zip(full).zip(first4) {
        let spectrum = openms::format::dta::read(Cursor::new(text)).unwrap();
        assert_eq!(spectrum.peaks.len().to_string(), row[2], "{}", row[1]);
        assert_eq!(
            PeakTypeEstimator::estimate_type(&spectrum.peaks).unwrap(),
            spectrum_type(row[3]),
            "{}",
            row[1]
        );
        assert_eq!(short[1], row[1]);
        assert_eq!(short[2], "4");
        assert_eq!(
            PeakTypeEstimator::estimate_type(&spectrum.peaks[..4]).unwrap(),
            spectrum_type(short[3])
        );
    }
}

#[cfg(feature = "mzml")]
fn load_mzml(label: &str) -> openms::MSExperiment {
    let relative = if label == "IM_FAIMS_test.mzML" {
        format!("tests/data/faims_helper/{label}")
    } else {
        format!("tests/data/{label}")
    };
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(relative);
    let file = std::fs::File::open(&path).unwrap();
    openms::format::mzml::read(std::io::BufReader::new(file)).unwrap()
}

#[cfg(feature = "mzml")]
fn mzml_labels() -> Vec<&'static str> {
    let mut labels: Vec<&str> = records("pte_mzml").iter().map(|row| row[1]).collect();
    labels.dedup();
    labels
}

#[cfg(feature = "mzml")]
#[test]
fn every_mzml_spectrum_matches_the_executed_oracle_and_msspectrum_get_type() {
    let rows = records("pte_mzml");
    let labels = mzml_labels();
    assert_eq!(labels.len(), 8);
    let mut compared = 0;
    for label in labels {
        let experiment = load_mzml(label);
        let expected: Vec<_> = rows.iter().filter(|row| row[1] == label).collect();
        assert_eq!(experiment.spectra.len(), expected.len(), "{label}");
        for (index, (spectrum, row)) in experiment.spectra.iter().zip(expected).enumerate() {
            assert_eq!(row[2], index.to_string());
            assert_eq!(spectrum.native_id, row[3], "{label} {index}");
            assert_eq!(spectrum.ms_level.to_string(), row[4], "{label} {index}");
            assert_eq!(spectrum.peaks.len().to_string(), row[5], "{label} {index}");
            let estimate = PeakTypeEstimator::estimate_type(&spectrum.peaks).unwrap();
            assert_eq!(estimate, spectrum_type(row[7]), "{label} {index}");
            // With nothing stored and no picking history, MSSpectrum::getType(true)
            // is exactly the estimator (MSSpectrum.cpp:142-165).
            let mut bare = spectrum.clone();
            bare.spectrum_type = SpectrumType::Unknown;
            bare.data_processing.clear();
            assert_eq!(bare.get_type(true).unwrap(), estimate, "{label} {index}");
            compared += 1;
        }
    }
    assert_eq!(compared, 65);
}

#[cfg(feature = "mzml")]
#[test]
fn stored_and_queried_types_match_the_executed_oracle() {
    // The oracle's getType(false) and getType(true) columns are reader
    // evidence rather than estimator evidence: they also depend on how each
    // reader stores the spectrum representation. None of these fixtures places
    // MS:1000525 inside a spectrum after MS:1000128, so the C++ reset to UNKNOWN
    // (MzMLHandler.cpp:1642-1645) is not exercised here.
    let rows = records("pte_mzml");
    for label in mzml_labels() {
        let experiment = load_mzml(label);
        let expected: Vec<_> = rows.iter().filter(|row| row[1] == label).collect();
        for (index, (spectrum, row)) in experiment.spectra.iter().zip(expected).enumerate() {
            assert_eq!(
                spectrum.spectrum_type,
                spectrum_type(row[6]),
                "{label} {index}"
            );
            assert_eq!(
                spectrum.get_type(true).unwrap(),
                spectrum_type(row[8]),
                "{label} {index}"
            );
        }
    }
}

#[test]
fn synthetic_shapes_match_the_executed_oracle_or_are_refused_when_not_finite() {
    let mut refused = Vec::new();
    let rows = records("pte_synthetic");
    assert_eq!(rows.len(), 11);
    for row in &rows {
        let peaks = synthetic(row[1]);
        assert_eq!(peaks.len().to_string(), row[2], "{}", row[1]);
        let finite = peaks
            .iter()
            .all(|peak| peak.mz.is_finite() && peak.intensity.is_finite());
        let result = PeakTypeEstimator::estimate_type(&peaks);
        if finite {
            assert_eq!(result.unwrap(), spectrum_type(row[3]), "{}", row[1]);
        } else {
            assert!(matches!(result, Err(Error::InvalidValue(_))), "{}", row[1]);
            refused.push((row[1], row[3]));
        }
    }
    // The source classifies both as CENTROID: a NaN is never a maximum and
    // poisons the total, and an infinite maximum makes every shoulder ratio NaN.
    assert_eq!(
        refused,
        [("nan_intensity", "Centroid"), ("inf_intensity", "Centroid")]
    );
}

#[test]
fn unsorted_duplicate_and_negative_input_is_classified_where_the_picker_refuses() {
    for name in [
        "descending_mz_profile",
        "duplicate_mz_profile",
        "negative_edges_profile",
    ] {
        let peaks = synthetic(name);
        assert_eq!(
            PeakTypeEstimator::estimate_type(&peaks).unwrap(),
            SpectrumType::Profile,
            "{name}"
        );
        assert!(
            estimate_spectrum_type(&MSSpectrum::from_peaks(peaks)).is_err(),
            "{name}"
        );
    }
}

#[test]
fn ceilings_and_their_order_are_those_of_msspectrum_get_type_with_limits() {
    let peaks = synthetic("ascending_profile");
    assert_eq!(peaks.len(), 7);
    let spectrum = MSSpectrum::from_peaks(peaks.clone());
    let generous = SpectrumTypeQueryLimits::default();
    let cases = [
        (
            SpectrumTypeQueryLimits {
                max_points: 6,
                ..generous
            },
            false,
        ),
        (
            SpectrumTypeQueryLimits {
                max_points: 7,
                ..generous
            },
            true,
        ),
        (
            SpectrumTypeQueryLimits {
                max_work: 7 * 32 - 1,
                ..generous
            },
            false,
        ),
        (
            SpectrumTypeQueryLimits {
                max_work: 7 * 32,
                ..generous
            },
            true,
        ),
        (
            SpectrumTypeQueryLimits {
                max_bytes: 7 * 16 - 1,
                ..generous
            },
            false,
        ),
        (
            SpectrumTypeQueryLimits {
                max_bytes: 7 * 16,
                ..generous
            },
            true,
        ),
    ];
    for (limits, accepted) in cases {
        let direct = PeakTypeEstimator::estimate_type_with_limits(&peaks, limits);
        let queried = spectrum.get_type_with_limits(true, limits);
        if accepted {
            assert_eq!(direct.unwrap(), SpectrumType::Profile, "{limits:?}");
            assert_eq!(queried.unwrap(), SpectrumType::Profile, "{limits:?}");
        } else {
            assert!(matches!(direct, Err(Error::InvalidValue(_))), "{limits:?}");
            assert!(matches!(queried, Err(Error::InvalidValue(_))), "{limits:?}");
        }
    }
    // Short input is Unknown before any ceiling is consulted, on both paths.
    let nothing = SpectrumTypeQueryLimits {
        max_points: 0,
        max_work: 0,
        max_bytes: 0,
    };
    assert_eq!(
        PeakTypeEstimator::estimate_type_with_limits(&peaks[..4], nothing).unwrap(),
        SpectrumType::Unknown
    );
    assert_eq!(
        MSSpectrum::from_peaks(peaks[..4].to_vec())
            .get_type_with_limits(true, nothing)
            .unwrap(),
        SpectrumType::Unknown
    );
}

#[test]
fn non_finite_values_are_refused_only_once_five_points_are_present() {
    let mut peaks = synthetic("five_flat_distinct");
    peaks[2].intensity = f32::NAN;
    // The source's short-input gate comes first, as in MSSpectrum::getType.
    assert_eq!(
        PeakTypeEstimator::estimate_type(&peaks[..4]).unwrap(),
        SpectrumType::Unknown
    );
    assert!(matches!(
        PeakTypeEstimator::estimate_type(&peaks),
        Err(Error::InvalidValue(_))
    ));
    assert!(MSSpectrum::from_peaks(peaks).get_type(true).is_err());

    let mut peaks = synthetic("five_flat_distinct");
    peaks[4].mz = f64::INFINITY;
    assert!(matches!(
        PeakTypeEstimator::estimate_type(&peaks),
        Err(Error::InvalidValue(_))
    ));
    peaks[4].mz = 500.0;
    peaks[0].intensity = f32::NEG_INFINITY;
    assert!(matches!(
        PeakTypeEstimator::estimate_type(&peaks),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn no_positive_intensity_and_the_one_thomson_window_follow_the_oracle() {
    // all_zero_five: no maximum is ever found, the evidence ratio is 0/0 = NaN,
    // and NaN > 0.75 is false. huge_mz_profile_shape: at 1e17, m/z + 1 == m/z,
    // so no shoulder point is within the window. ratio_exact_tenth: 10/100 is
    // not strictly more than 0.1.
    for (name, expected) in [
        ("all_zero_five", SpectrumType::Centroid),
        ("huge_mz_profile_shape", SpectrumType::Centroid),
        ("ratio_exact_tenth", SpectrumType::Centroid),
        ("five_flat_distinct", SpectrumType::Centroid),
        ("four_points", SpectrumType::Unknown),
    ] {
        assert_eq!(
            PeakTypeEstimator::estimate_type(&synthetic(name)).unwrap(),
            expected,
            "{name}"
        );
    }
}
