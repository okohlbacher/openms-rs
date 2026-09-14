// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The `START_SECTION`s of `IMDataConverter_test.cpp` that exercise
//! `splitByFAIMSCV` at Core SDK `bc9cc12`, the executed product-sdk oracle for
//! it, the C2 records of the FeatureFinderCentroided FAIMS path, and the native
//! boundaries of `src/kernel/im_data_converter.rs`.
//!
//! Evidence labels:
//!
//! - Expectations citing an `IMDataConverter_test.cpp` line are transcribed
//!   class-test literals (tier 3, source review).
//! - `tests/data/im_data_converter/oracle_cases.tsv` is oracle-generated (tier 1
//!   executed differential): the unmodified product-sdk `libOpenMS` (Debug, core
//!   4fdec46, identical to `bc9cc12` for every file involved) run by
//!   `../oracle/im-data-converter/driver.cpp` on six mzML files and 20 synthetic
//!   experiments.
//! - `tests/data/im_data_converter/c2_split_records.tsv` holds the split records
//!   of the C2 class-level oracle (`../oracle/featurefinder-picked/`, tier 1
//!   executed differential), extracted from its `faims_corrected_*.jsonl` and
//!   `faims_facts.jsonl`.
//! - The NaN refusal, the settings ceiling, the atomicity on error and the
//!   returned skipped spectra and chromatograms are native (tier 4). Where the
//!   source gives a result for input this port refuses, the test asserts the
//!   refusal *and* the recorded source result, so the divergence stays visible.

use openms::kernel::faims_helper::{CompensationVoltage, FaimsHelper};
use openms::kernel::im_data_converter::{
    FaimsGroupKey, FaimsSplit, FaimsSplitLogLevel, FaimsSplitMessage, ImDataConverter,
};
use openms::metadata::{DriftTimeUnit, ImTypes, MetaValue, Sample, to_drift_time_unit};
use openms::{
    ChromatogramPeak, Error, MSChromatogram, MSExperiment, MSSpectrum, Peak1D,
    data_structures::DateTime,
};
use std::collections::BTreeMap;

const ORACLE: &str = include_str!("data/im_data_converter/oracle_cases.tsv");
const C2_RECORDS: &str = include_str!("data/im_data_converter/c2_split_records.tsv");

const F: DriftTimeUnit = DriftTimeUnit::FaimsCompensationVoltage;
const N: DriftTimeUnit = DriftTimeUnit::None;
const MS: DriftTimeUnit = DriftTimeUnit::Millisecond;
const U: f64 = ImTypes::DRIFTTIME_NOT_SET;

// ---------------------------------------------------------------- oracle records

/// The records of one kind, split into tab-separated fields.
fn records<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|fields| fields[0] == kind)
        .collect()
}

/// The oracle records of one kind for one case (field 1 is the case name).
fn case_records(kind: &str, case: &str) -> Vec<Vec<&'static str>> {
    records(ORACLE, kind)
        .into_iter()
        .filter(|fields| fields[1] == case)
        .collect()
}

/// Parses the C `printf("%a")` spelling the oracle writes (`-0x1.68p+5`,
/// `-0x0p+0`, `inf`, `nan`). Every mantissa has at most 53 significant bits, so
/// the conversion is exact.
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

/// The oracle's spelling of a group key: `nan` for the unsplit group.
fn key_matches(key: FaimsGroupKey, recorded: &str) -> bool {
    match key {
        FaimsGroupKey::NotFaims => recorded == "nan",
        FaimsGroupKey::Voltage(voltage) => {
            recorded != "nan" && voltage.volts().to_bits() == hex_f64(recorded).to_bits()
        }
    }
}

/// The lines a C++ `LogStreamBuf` prints for `lines` written in order and a
/// final `clearCache()`, which is how the oracle driver captured each case.
///
/// `LogStream.cpp:299-360` (`syncLF_`): a line already in the repetition cache
/// only increments its counter; a new line first evicts the oldest entry when
/// the cache holds more than one (printing `<line> occurred N times` if it was
/// repeated, `addToCache_`, lines 207-234), then is printed. `clearCache`
/// (lines 236-253) prints the repeated entries in `std::map<std::string>` key
/// order.
fn emulate_log_stream(lines: &[String]) -> Vec<String> {
    struct Entry {
        counter: usize,
        timestamp: usize,
    }
    let mut cache: BTreeMap<String, Entry> = BTreeMap::new();
    let mut clock = 0usize;
    let mut printed = Vec::new();
    for line in lines {
        clock += 1;
        if let Some(entry) = cache.get_mut(line) {
            entry.counter += 1;
            entry.timestamp = clock;
            continue;
        }
        if cache.len() > 1 {
            let oldest = cache
                .iter()
                .min_by_key(|(_, entry)| entry.timestamp)
                .map(|(text, _)| text.clone())
                .unwrap();
            let entry = cache.remove(&oldest).unwrap();
            if entry.counter != 0 {
                printed.push(format!("<{oldest}> occurred {} times", entry.counter + 1));
            }
        }
        clock += 1;
        cache.insert(
            line.clone(),
            Entry {
                counter: 0,
                timestamp: clock,
            },
        );
        printed.push(line.clone());
    }
    for (text, entry) in &cache {
        if entry.counter != 0 {
            printed.push(format!("<{text}> occurred {} times", entry.counter + 1));
        }
    }
    printed
}

// ---------------------------------------------------------------- experiments

/// One synthetic spectrum per `(MS level, unit, drift time)`, exactly as the
/// oracle driver's `build()`: native ID `s<index>`, RT `10 * index`, one peak at
/// m/z `500 + index` with intensity 100.
fn build(spectra: &[(u32, DriftTimeUnit, f64)]) -> MSExperiment {
    let mut experiment = MSExperiment::new();
    for (index, &(ms_level, unit, drift_time)) in spectra.iter().enumerate() {
        experiment.spectra.push(MSSpectrum {
            native_id: format!("s{index}"),
            ms_level,
            rt: 10.0 * index as f64,
            drift_time,
            drift_time_unit: unit,
            peaks: vec![Peak1D::new(500.0 + index as f64, 100.0)],
            ..MSSpectrum::default()
        });
    }
    experiment
}

/// The settings-and-chromatograms experiment of the oracle driver.
fn with_settings_and_chromatogram(mut experiment: MSExperiment) -> MSExperiment {
    experiment.settings.comment = "settings comment".into();
    experiment.settings.date_time = DateTime::parse("2019-09-07T09:40:04").unwrap();
    experiment.settings.document.identifier = "document_identifier".into();
    experiment.settings.fraction_identifier = "fraction_1".into();
    experiment.sql_run_id = 42;
    experiment
        .settings
        .metadata
        .insert("custom_key".into(), MetaValue::from("custom_value"));
    experiment.chromatograms.push(MSChromatogram {
        native_id: "chromatogram_1".into(),
        peaks: vec![ChromatogramPeak::new(1.0, 2.0)],
        ..MSChromatogram::default()
    });
    experiment
}

/// The synthetic oracle cases other than the two class-test experiments and
/// the NaN cases, in driver order.
fn synthetic_cases() -> Vec<(&'static str, MSExperiment)> {
    let inf = f64::INFINITY;
    vec![
        ("empty", MSExperiment::new()),
        (
            "ms2_before_any_faims",
            build(&[(2, N, U), (1, F, -45.0), (2, N, U)]),
        ),
        (
            "cv_less_ms1_keeps_context",
            build(&[
                (1, F, -45.0),
                (1, N, U),
                (2, N, U),
                (1, MS, 12.5),
                (2, N, U),
            ]),
        ),
        (
            "sentinel_faims_spectrum",
            build(&[
                (1, F, -45.0),
                (2, N, U),
                (1, F, U),
                (2, N, U),
                (1, F, -60.0),
                (2, N, U),
            ]),
        ),
        ("only_sentinel", build(&[(1, F, U), (2, N, U)])),
        (
            "signed_zero_negative_first",
            build(&[(1, F, -0.0), (2, N, U), (1, F, 0.0), (2, N, U)]),
        ),
        (
            "signed_zero_positive_first",
            build(&[(1, F, 0.0), (2, N, U), (1, F, -0.0), (2, N, U)]),
        ),
        (
            "infinities",
            build(&[(1, F, inf), (1, F, -inf), (1, F, -45.0), (2, N, U)]),
        ),
        (
            "ms2_with_own_cv",
            build(&[
                (1, F, -45.0),
                (2, F, -60.0),
                (2, N, U),
                (1, F, -45.0),
                (2, N, U),
            ]),
        ),
        (
            "ms_level_zero_and_three",
            build(&[(1, F, -45.0), (0, N, U), (2, N, U), (3, N, U)]),
        ),
        (
            "interleaved_descending_cvs",
            build(&[
                (1, F, -30.0),
                (2, N, U),
                (1, F, -50.0),
                (2, N, U),
                (1, F, -70.0),
                (2, N, U),
                (1, F, -30.0),
            ]),
        ),
        (
            "unit_none_with_value",
            build(&[(1, N, -50.0), (1, F, -40.0), (2, N, -50.0)]),
        ),
        (
            "adjacent_doubles",
            build(&[(1, F, -50.0), (1, F, -49.99999999999999), (2, N, U)]),
        ),
        (
            "settings_and_chromatograms_faims",
            with_settings_and_chromatogram(build(&[
                (1, F, -45.0),
                (2, N, U),
                (1, F, -60.0),
                (1, N, U),
            ])),
        ),
        (
            "settings_and_chromatograms_non_faims",
            with_settings_and_chromatogram(build(&[(1, N, U), (2, N, U)])),
        ),
    ]
}

/// The class-test experiment of `IMDataConverter_test.cpp:205-231`: spectra
/// without peaks, native IDs or retention times.
fn class_ms2_last_seen_cv() -> MSExperiment {
    let spectrum = |ms_level: u32, unit: DriftTimeUnit, drift_time: f64| MSSpectrum {
        ms_level,
        drift_time_unit: unit,
        drift_time,
        ..MSSpectrum::default()
    };
    let mut experiment = MSExperiment::new();
    // ms1a, ms2a, ms2b, ms1b, ms2c; an MS2 spectrum's drift time stays unset.
    experiment.spectra = vec![
        spectrum(1, F, -55.0),
        spectrum(2, N, U),
        spectrum(2, N, U),
        spectrum(1, F, -45.0),
        spectrum(2, N, U),
    ];
    experiment
}

/// The class-test experiment of `IMDataConverter_test.cpp:259-263`.
fn class_non_faims() -> MSExperiment {
    let mut experiment = MSExperiment::new();
    experiment.spectra.push(MSSpectrum {
        ms_level: 1,
        drift_time_unit: MS,
        drift_time: 10.0,
        ..MSSpectrum::default()
    });
    experiment
}

#[cfg(feature = "mzml")]
fn fixture_bytes(path: &str) -> Vec<u8> {
    std::fs::read(std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(path)).unwrap()
}

#[cfg(feature = "mzml")]
fn read_mzml(bytes: &[u8]) -> MSExperiment {
    openms::format::mzml::read(bytes).unwrap()
}

/// FNV-1a, 64 bit: a cheap binding of a derived fixture to the bytes the C2
/// oracle executed on (whose sha256 the provenance manifest records).
#[cfg(feature = "mzml")]
fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, &byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

/// The C2 synthetic FAIMS copy of the FeatureFinderCentroided_1 input
/// (`../oracle/featurefinder-picked/results/inputs/synthetic/faims_{one,two}_cv.mzML`):
/// after every `<scan>` line an unindented scan-level `MS:1001581` cvParam,
/// cycling through `voltages`. The rule reproduces both C2 files byte for byte
/// (checked against their sha256 when the fixture was prepared; bound here by
/// length and FNV-1a).
#[cfg(feature = "mzml")]
fn derive_c2_fixture(voltages: &[&str]) -> Vec<u8> {
    let base = fixture_bytes("tests/data/mzml_mobility/FeatureFinderCentroided_1_input.mzML");
    let mut out = Vec::with_capacity(base.len() + 112 * 140);
    let mut inserted = 0usize;
    let lines: Vec<&[u8]> = base.split(|&byte| byte == b'\n').collect();
    for (index, line) in lines.iter().enumerate() {
        out.extend_from_slice(line);
        if index + 1 < lines.len() {
            out.push(b'\n');
        }
        if line.trim_ascii() == b"<scan>" {
            let volts = voltages[inserted % voltages.len()];
            inserted += 1;
            out.extend_from_slice(
                format!(
                    "<cvParam cvRef=\"MS\" accession=\"MS:1001581\" name=\"FAIMS compensation voltage\" value=\"{volts}\" unitAccession=\"UO:0000218\" unitName=\"volt\" unitCvRef=\"UO\"/>\n"
                )
                .as_bytes(),
            );
        }
    }
    assert_eq!(inserted, 112);
    out
}

// ---------------------------------------------------------------- comparison

/// The input records of a file case: native ID, MS level, drift time bits and
/// unit per spectrum, as the C++ reader produced them.
#[cfg(feature = "mzml")]
fn assert_oracle_input(case: &str, experiment: &MSExperiment) {
    let rows = case_records("input", case);
    assert_eq!(rows.len(), experiment.spectra.len(), "{case}");
    for (index, (row, spectrum)) in rows.iter().zip(&experiment.spectra).enumerate() {
        assert_eq!(row[2], index.to_string(), "{case}");
        assert_eq!(spectrum.native_id, row[3], "{case} {index}");
        assert_eq!(spectrum.ms_level.to_string(), row[4], "{case} {index}");
        assert_eq!(
            spectrum.drift_time.to_bits(),
            hex_f64(row[5]).to_bits(),
            "{case} {index}"
        );
        assert_eq!(
            spectrum.drift_time_unit,
            to_drift_time_unit(row[6]).unwrap(),
            "{case} {index}"
        );
    }
}

/// Splits `experiment` and asserts every oracle record of `case`: the case
/// record, each group record, each spectrum record in group order and the
/// printed log lines. Returns the split for further native assertions.
fn split_and_compare(case: &str, mut experiment: MSExperiment) -> FaimsSplit {
    let before = experiment.clone();
    let split = ImDataConverter::split_by_faims_cv(&mut experiment).unwrap();

    let case_rows = case_records("case", case);
    assert_eq!(case_rows.len(), 1, "{case}");
    let row = &case_rows[0];
    assert_eq!(row[2], before.spectra.len().to_string(), "{case}");
    assert_eq!(row[3], before.chromatograms.len().to_string(), "{case}");
    assert_eq!(row[4], split.groups.len().to_string(), "{case}");
    // The source empties the input in every case, chromatograms included.
    assert_eq!(row[5], "true", "{case}");
    assert_eq!(row[6], "0", "{case}");
    assert!(experiment.is_empty() && experiment.chromatograms.is_empty());
    // The source resets the settings with clear(true) when it splits; for the
    // unsplit case they are in an unspecified moved-from state (native
    // difference: the port always leaves the default).
    assert_eq!(experiment, MSExperiment::default(), "{case}");
    if split.has_faims() {
        assert_eq!(row[7], "true", "{case}");
    }

    let group_rows = case_records("group", case);
    assert_eq!(group_rows.len(), split.groups.len(), "{case}");
    let spectrum_rows = case_records("spectrum", case);
    let mut next_spectrum_row = 0usize;
    for (index, (row, group)) in group_rows.iter().zip(&split.groups).enumerate() {
        assert_eq!(row[2], index.to_string(), "{case}");
        assert!(
            key_matches(group.key, row[3]),
            "{case} group {index}: {:?} vs {}",
            group.key,
            row[3]
        );
        assert_eq!(
            row[4],
            group.experiment.spectra.len().to_string(),
            "{case} group {index}"
        );
        assert_eq!(
            row[5],
            group.experiment.chromatograms.len().to_string(),
            "{case} group {index}"
        );
        assert_eq!(row[6], "true", "{case} group {index}");
        assert_eq!(
            group.experiment.settings, before.settings,
            "{case} group {index}"
        );
        assert_eq!(
            group.experiment.settings.date_time.iso_string(),
            row[7],
            "{case} group {index}"
        );
        assert_eq!(
            group.experiment.settings.comment, row[8],
            "{case} group {index}"
        );
        assert_eq!(
            group.experiment.sql_run_id.to_string(),
            row[9],
            "{case} group {index}"
        );
        for (position, spectrum) in group.experiment.spectra.iter().enumerate() {
            let recorded = &spectrum_rows[next_spectrum_row];
            next_spectrum_row += 1;
            assert_eq!(recorded[2], index.to_string(), "{case}");
            assert_eq!(recorded[3], position.to_string(), "{case}");
            assert_eq!(
                spectrum.native_id, recorded[4],
                "{case} group {index} {position}"
            );
            assert_eq!(spectrum.ms_level.to_string(), recorded[5], "{case}");
            assert_eq!(
                spectrum.drift_time.to_bits(),
                hex_f64(recorded[6]).to_bits(),
                "{case}"
            );
            assert_eq!(
                spectrum.drift_time_unit,
                to_drift_time_unit(recorded[7]).unwrap(),
                "{case}"
            );
        }
        // Native difference: a source voltage group has no ranges
        // (`byMSLevel(1)` throws, the recorded `true`); here they are computed
        // on demand, so a group holding an MS1 peak always has MS1 ranges.
        if group.key != FaimsGroupKey::NotFaims {
            assert_eq!(row[10], "true", "{case} group {index}");
        }
        // (The range manager refuses a non-finite drift time, so the infinite
        // voltage groups are left out.)
        let spectra = &group.experiment.spectra;
        if spectra
            .iter()
            .any(|spectrum| spectrum.ms_level == 1 && !spectrum.peaks.is_empty())
            && spectra
                .iter()
                .all(|spectrum| spectrum.drift_time.is_finite())
        {
            let ranges = group.experiment.spectrum_range_manager().unwrap();
            assert!(ranges.by_ms_level(1).is_some(), "{case} group {index}");
        }
    }
    assert_eq!(next_spectrum_row, spectrum_rows.len(), "{case}");

    // Skipped spectra and dropped chromatograms: what the source destroys.
    let skipped_indices: Vec<usize> = split
        .messages
        .iter()
        .filter_map(FaimsSplitMessage::spectrum_index)
        .collect();
    assert_eq!(skipped_indices.len(), split.skipped_spectra.len(), "{case}");
    for (index, spectrum) in skipped_indices.iter().zip(&split.skipped_spectra) {
        assert_eq!(&before.spectra[*index], spectrum, "{case}");
    }
    let grouped: usize = split
        .groups
        .iter()
        .map(|group| group.experiment.spectra.len())
        .sum();
    assert_eq!(
        grouped + split.skipped_spectra.len(),
        before.spectra.len(),
        "{case}"
    );
    if split.has_faims() {
        assert_eq!(split.dropped_chromatograms, before.chromatograms, "{case}");
    } else {
        assert!(split.dropped_chromatograms.is_empty() && split.skipped_spectra.is_empty());
        assert_eq!(split.groups[0].experiment, before, "{case}");
    }

    // Log lines, per channel, through the C++ log stream's repetition filter.
    for (level, channel) in [
        (FaimsSplitLogLevel::Warning, "warn"),
        (FaimsSplitLogLevel::Info, "info"),
    ] {
        let lines: Vec<String> = split
            .messages
            .iter()
            .filter(|message| message.level() == level)
            .map(FaimsSplitMessage::text)
            .collect();
        let recorded: Vec<String> = case_records("log", case)
            .iter()
            .filter(|fields| fields[2] == channel)
            .map(|fields| fields[3].to_owned())
            .collect();
        assert_eq!(emulate_log_stream(&lines), recorded, "{case} {channel}");
    }
    split
}

// ---------------------------------------------------------------- class test

#[test]
fn constructor_and_destructor_sections() {
    // IMDataConverter_test.cpp:31-38 allocates and deletes an instance of the
    // stateless class; the Rust type is a zero-sized unit struct.
    let constructed = <ImDataConverter as Default>::default();
    assert_eq!(constructed, ImDataConverter);
    assert_eq!(std::mem::size_of::<ImDataConverter>(), 0);
    assert!(!std::mem::needs_drop::<ImDataConverter>());
}

#[cfg(feature = "mzml")]
#[test]
fn split_by_faims_cv_section() {
    // IMDataConverter_test.cpp:41-76.
    let mut exp = read_mzml(&fixture_bytes("tests/data/faims_helper/IM_FAIMS_test.mzML"));
    assert_eq!(exp.spectra.len(), 19);

    let split_peak_map = ImDataConverter::split_by_faims_cv(&mut exp).unwrap();
    assert!(exp.is_empty()); // moved out
    let groups = &split_peak_map.groups;
    assert_eq!(groups.len(), 3);

    // expect keys -65, -55, -45 in ascending order
    assert_eq!(groups[0].key.volts(), -65.0);
    assert_eq!(groups[1].key.volts(), -55.0);
    assert_eq!(groups[2].key.volts(), -45.0);

    assert_eq!(groups[0].experiment.len(), 4);
    assert_eq!(groups[1].experiment.len(), 9);
    assert_eq!(groups[2].experiment.len(), 6);

    for (group, volts) in groups.iter().zip([-65.0, -55.0, -45.0]) {
        for spectrum in &group.experiment.spectra {
            assert_eq!(spectrum.drift_time, volts);
        }
    }

    assert_eq!(
        groups[1].experiment.settings.date_time.iso_string(),
        "2019-09-07T09:40:04"
    );
}

#[test]
fn split_by_faims_cv_assigns_ms2_without_explicit_cv_to_last_seen_faims_cv_section() {
    // IMDataConverter_test.cpp:201-254.
    let mut exp_synth = class_ms2_last_seen_cv();
    let bins = ImDataConverter::split_by_faims_cv(&mut exp_synth)
        .unwrap()
        .groups;
    assert_eq!(bins.len(), 2);

    // bins ordered by ascending CV: -55 first, -45 second
    assert_eq!(bins[0].key.volts(), -55.0);
    assert_eq!(bins[1].key.volts(), -45.0);

    let bin_minus55 = &bins[0].experiment;
    let bin_minus45 = &bins[1].experiment;

    assert_eq!(bin_minus55.len(), 3); // ms1a + ms2a + ms2b
    assert_eq!(bin_minus45.len(), 2); // ms1b + ms2c

    let ms2_bin0 = bin_minus55
        .spectra
        .iter()
        .filter(|s| s.ms_level > 1)
        .count();
    assert_eq!(ms2_bin0, 2);
    let ms2_bin1 = bin_minus45
        .spectra
        .iter()
        .filter(|s| s.ms_level > 1)
        .count();
    assert_eq!(ms2_bin1, 1);
}

#[test]
fn split_by_faims_cv_returns_single_element_group_for_non_faims_dataset_section() {
    // IMDataConverter_test.cpp:256-269.
    let mut exp_nonfaims = class_non_faims();
    let bins = ImDataConverter::split_by_faims_cv(&mut exp_nonfaims)
        .unwrap()
        .groups;
    assert_eq!(bins.len(), 1);
    assert_eq!(bins[0].experiment.len(), 1);
}

// ---------------------------------------------------------------- executed oracle

#[test]
fn class_test_experiments_match_the_oracle() {
    split_and_compare("class_ms2_last_seen_cv", class_ms2_last_seen_cv());
    let split = split_and_compare("class_non_faims", class_non_faims());
    assert_eq!(split.messages, [FaimsSplitMessage::NoCompensationVoltages]);
}

#[test]
fn synthetic_experiments_match_the_oracle() {
    let cases = synthetic_cases();
    // Every synthetic oracle case is covered: the file cases, the two class
    // experiments above and the three NaN cases below account for the rest.
    let recorded = records(ORACLE, "case").len();
    assert_eq!(recorded, 6 + 2 + cases.len() + 3);
    for (case, experiment) in cases {
        split_and_compare(case, experiment);
    }
}

#[test]
fn the_adjacent_doubles_case_uses_the_driver_values() {
    // The driver writes std::nextafter(-50.0, 0.0); the literal above is that
    // double.
    assert_eq!((-49.99999999999999f64).to_bits(), (-50.0f64).to_bits() - 1);
}

#[test]
fn sentinel_case_records_in_detail() {
    let split = split_and_compare(
        "sentinel_faims_spectrum",
        build(&[
            (1, F, -45.0),
            (2, N, U),
            (1, F, U),
            (2, N, U),
            (1, F, -60.0),
            (2, N, U),
        ]),
    );
    assert_eq!(
        split.messages,
        [
            FaimsSplitMessage::CompensationVoltageWarning(
                FaimsHelper::MISSING_VOLTAGE_WARNING.into()
            ),
            FaimsSplitMessage::UnexpectedCompensationVoltage {
                spectrum_index: 2,
                volts: -1.0
            },
            FaimsSplitMessage::SpectrumWithoutCompensationVoltage { spectrum_index: 3 },
        ]
    );
    assert_eq!(
        split.messages[1].text(),
        "Encountered spectrum with unexpected FAIMS CV (not in detected set): -1"
    );
    let native: Vec<&str> = split
        .skipped_spectra
        .iter()
        .map(|s| s.native_id.as_str())
        .collect();
    assert_eq!(native, ["s2", "s3"]);
}

#[test]
fn nan_voltages_are_refused_where_the_source_breaks_its_ordering() {
    let cases = [
        (
            "nan_first",
            build(&[(1, F, f64::NAN), (1, F, -50.0), (1, F, -40.0)]),
        ),
        (
            "nan_middle",
            build(&[(1, F, -50.0), (1, F, f64::NAN), (2, N, U), (1, F, -40.0)]),
        ),
        (
            "nan_last",
            build(&[(1, F, -50.0), (1, F, -40.0), (1, F, f64::NAN), (2, N, U)]),
        ),
    ];
    for (case, experiment) in cases {
        let mut input = experiment.clone();
        let error = ImDataConverter::split_by_faims_cv(&mut input).unwrap_err();
        assert!(matches!(error, Error::InvalidValue(_)), "{case}");
        // Debug text, because NaN != NaN defeats PartialEq on the spectra.
        assert_eq!(
            format!("{input:?}"),
            format!("{experiment:?}"),
            "{case}: the input is unchanged"
        );
    }
    // The executed source results the refusal replaces: a NaN first leaves the
    // input unsplit under the NaN key with a spurious missing-voltage warning;
    // a later NaN spectrum joins the -50 group and the MS2 after it is skipped.
    let groups = |case: &str| -> Vec<(String, String)> {
        case_records("group", case)
            .iter()
            .map(|fields| (fields[3].to_owned(), fields[4].to_owned()))
            .collect()
    };
    assert_eq!(groups("nan_first"), [("nan".to_owned(), "3".to_owned())]);
    for case in ["nan_middle", "nan_last"] {
        assert_eq!(
            groups(case),
            [
                ("-0x1.9p+5".to_owned(), "2".to_owned()),
                ("-0x1.4p+5".to_owned(), "1".to_owned())
            ]
        );
        let members: Vec<&str> = case_records("spectrum", case)
            .iter()
            .filter(|fields| fields[2] == "0")
            .map(|fields| fields[6])
            .collect();
        assert!(members.contains(&"nan"), "{case}");
    }
}

#[cfg(feature = "mzml")]
#[test]
fn mzml_files_match_the_oracle() {
    let files = [
        (
            "im_faims_test",
            "tests/data/faims_helper/IM_FAIMS_test.mzML",
        ),
        (
            "ffc1_input",
            "tests/data/mzml_mobility/FeatureFinderCentroided_1_input.mzML",
        ),
        (
            "faims_test_data",
            "tests/data/mzml_mobility/FAIMS_test_data.mzML",
        ),
        (
            "faims_interleaved",
            "tests/data/mzml_mobility/FAIMS_CV-60C_V-45_Interleaved.mzML",
        ),
    ];
    for (case, path) in files {
        let experiment = read_mzml(&fixture_bytes(path));
        assert_oracle_input(case, &experiment);
        let split = split_and_compare(case, experiment);
        assert!(split.skipped_spectra.is_empty(), "{case}");
    }
}

#[cfg(feature = "mzml")]
#[test]
fn c2_fixtures_match_the_oracle_and_the_c2_records() {
    let one = derive_c2_fixture(&["-45"]);
    let two = derive_c2_fixture(&["-45", "-60"]);
    assert_eq!((one.len(), fnv1a64(&one)), (229_778, 0xd603_0001_77f0_1520));
    assert_eq!((two.len(), fnv1a64(&two)), (229_778, 0xf472_73d3_ff22_2316));

    let c2_inputs = records(C2_RECORDS, "c2_input");
    let c2_groups = records(C2_RECORDS, "c2_group");
    for (case, label, bytes) in [
        ("c2_faims_one_cv", "one_cv", one),
        ("c2_faims_two_cv", "two_cv", two),
        (
            "ffc1_input",
            "ffc1_non_faims",
            fixture_bytes("tests/data/mzml_mobility/FeatureFinderCentroided_1_input.mzML"),
        ),
    ] {
        let experiment = read_mzml(&bytes);
        assert_oracle_input(case, &experiment);

        let input = c2_inputs.iter().find(|fields| fields[1] == label).unwrap();
        assert_eq!(input[3], experiment.spectra.len().to_string(), "{label}");
        let voltages: Vec<String> = FaimsHelper::get_compensation_voltages(&experiment)
            .unwrap()
            .values()
            .map(|volts| format!("{:x}", volts.to_bits()))
            .collect();
        let recorded: Vec<String> = if input[4] == "-" {
            Vec::new()
        } else {
            input[4]
                .split(',')
                .map(|text| format!("{:x}", hex_f64(text).to_bits()))
                .collect()
        };
        assert_eq!(voltages, recorded, "{label}");

        let split = split_and_compare(case, experiment);
        let groups: Vec<&Vec<&str>> = c2_groups
            .iter()
            .filter(|fields| fields[1] == label)
            .collect();
        assert_eq!(groups.len(), split.groups.len(), "{label}");
        for (index, (fields, group)) in groups.iter().zip(&split.groups).enumerate() {
            assert_eq!(fields[2], "source_sequence");
            assert_eq!(fields[3], index.to_string());
            assert!(key_matches(group.key, fields[4]), "{label} {index}");
            assert_eq!(
                fields[5],
                group.experiment.len().to_string(),
                "{label} {index}"
            );
            // C2: every voltage group throws at byMSLevel(1) before
            // updateRanges; the unsplit loaded experiment does not.
            let throws = group.key != FaimsGroupKey::NotFaims;
            assert_eq!(fields[6], throws.to_string(), "{label} {index}");
            assert!(
                group
                    .experiment
                    .spectrum_range_manager()
                    .unwrap()
                    .by_ms_level(1)
                    .is_some()
            );
        }
    }
}

#[test]
fn c2_synthetic_fact_splits_as_recorded() {
    // faims_facts.cpp fact 1a: four MS1 spectra alternating -45 and -60, RT
    // 10 * i, one peak (500, 100).
    let mut experiment = MSExperiment::new();
    for index in 0..4 {
        experiment.spectra.push(MSSpectrum {
            ms_level: 1,
            rt: 10.0 * f64::from(index),
            drift_time: if index % 2 == 1 { -60.0 } else { -45.0 },
            drift_time_unit: F,
            peaks: vec![Peak1D::new(500.0, 100.0)],
            ..MSSpectrum::default()
        });
    }
    let split = ImDataConverter::split_by_faims_cv(&mut experiment).unwrap();
    let recorded: Vec<Vec<&str>> = records(C2_RECORDS, "c2_fact")
        .into_iter()
        .filter(|fields| fields[1] == "split_groups_have_no_ms_level_ranges_synthetic")
        .collect();
    assert_eq!(recorded.len(), split.groups.len());
    for (fields, group) in recorded.iter().zip(&split.groups) {
        assert!(key_matches(group.key, fields[3]));
        assert_eq!(fields[4], group.experiment.len().to_string());
        assert_eq!(fields[5], "true");
        let ranges = group.experiment.spectrum_range_manager().unwrap();
        assert!(ranges.by_ms_level(1).is_some());
    }
    // Fact 1b records the same two groups on the two-CV FFC_1 copy.
    let facts_1b = records(C2_RECORDS, "c2_fact")
        .into_iter()
        .filter(|fields| fields[1] == "split_groups_have_no_ms_level_ranges_two_cv_ffc1")
        .map(|fields| (fields[3], fields[4]))
        .collect::<Vec<_>>();
    assert_eq!(facts_1b, [("-0x1.ep+5", "56"), ("-0x1.68p+5", "56")]);
}

// ---------------------------------------------------------------- native boundaries

#[test]
fn group_keys() {
    assert!(FaimsGroupKey::NotFaims.volts().is_nan());
    assert_eq!(FaimsGroupKey::NotFaims.voltage(), None);
    let voltage = CompensationVoltage::new(-0.0).unwrap();
    let key = FaimsGroupKey::Voltage(voltage);
    assert_eq!(key.voltage(), Some(voltage));
    assert_eq!(key.volts().to_bits(), (-0.0f64).to_bits());
    assert!(FaimsGroupKey::NotFaims < key);
    assert_eq!(
        key,
        FaimsGroupKey::Voltage(CompensationVoltage::new(0.0).unwrap())
    );
}

#[test]
fn message_texts_and_levels() {
    let helper =
        FaimsSplitMessage::CompensationVoltageWarning(FaimsHelper::MISSING_VOLTAGE_WARNING.into());
    assert_eq!(helper.level(), FaimsSplitLogLevel::Warning);
    assert_eq!(helper.text(), FaimsHelper::MISSING_VOLTAGE_WARNING);
    assert_eq!(helper.spectrum_index(), None);

    let info = FaimsSplitMessage::NoCompensationVoltages;
    assert_eq!(info.level(), FaimsSplitLogLevel::Info);
    assert_eq!(
        info.to_string(),
        "Not FAIMS compensation voltages found in the data. Returning PeakMap as CV NaN."
    );
    assert_eq!(info.spectrum_index(), None);

    let unexpected = FaimsSplitMessage::UnexpectedCompensationVoltage {
        spectrum_index: 7,
        volts: -45.5,
    };
    assert_eq!(unexpected.level(), FaimsSplitLogLevel::Warning);
    assert_eq!(unexpected.spectrum_index(), Some(7));
    assert_eq!(
        unexpected.text(),
        "Encountered spectrum with unexpected FAIMS CV (not in detected set): -45.5"
    );

    let without = FaimsSplitMessage::SpectrumWithoutCompensationVoltage { spectrum_index: 3 };
    assert_eq!(without.level(), FaimsSplitLogLevel::Warning);
    assert_eq!(without.spectrum_index(), Some(3));
    assert_eq!(
        without.to_string(),
        "Skipping spectrum without FAIMS CV (no prior FAIMS CV context or unexpected layout)."
    );
}

#[test]
fn voltage_groups_copy_settings_and_run_id_and_return_what_the_source_drops() {
    let input = with_settings_and_chromatogram(build(&[
        (1, F, -45.0),
        (2, N, U),
        (1, F, -60.0),
        (1, N, U),
    ]));
    let mut experiment = input.clone();
    let split = ImDataConverter::split_by_faims_cv(&mut experiment).unwrap();
    assert_eq!(experiment, MSExperiment::default());
    assert!(split.has_faims());
    for group in &split.groups {
        assert_eq!(group.experiment.settings, input.settings);
        assert_eq!(group.experiment.sql_run_id, 42);
        assert!(group.experiment.chromatograms.is_empty());
    }
    assert_eq!(split.dropped_chromatograms, input.chromatograms);
    assert_eq!(split.skipped_spectra, [input.spectra[3].clone()]);

    let input = with_settings_and_chromatogram(build(&[(1, N, U), (2, N, U)]));
    let mut experiment = input.clone();
    let split = ImDataConverter::split_by_faims_cv(&mut experiment).unwrap();
    assert!(!split.has_faims());
    assert_eq!(split.groups.len(), 1);
    assert_eq!(split.groups[0].experiment, input);
    assert!(split.dropped_chromatograms.is_empty() && split.skipped_spectra.is_empty());
}

#[test]
fn an_empty_experiment_is_one_empty_unsplit_group() {
    let mut experiment = MSExperiment::new();
    let split = ImDataConverter::split_by_faims_cv(&mut experiment).unwrap();
    assert_eq!(split.groups.len(), 1);
    assert_eq!(split.groups[0].key, FaimsGroupKey::NotFaims);
    assert!(split.groups[0].experiment.is_empty());
    assert_eq!(split.messages, [FaimsSplitMessage::NoCompensationVoltages]);
    assert!(!split.has_faims());
}

#[test]
fn settings_beyond_the_resource_limits_are_refused_before_anything_moves() {
    let mut deep = Sample::default();
    for _ in 0..70 {
        deep = Sample {
            subsamples: vec![deep],
            ..Sample::default()
        };
    }
    let mut input = build(&[(1, F, -45.0), (2, N, U)]);
    input.settings.sample = deep;
    let mut experiment = input.clone();
    assert!(ImDataConverter::split_by_faims_cv(&mut experiment).is_err());
    assert_eq!(experiment, input);

    // Without voltages nothing is copied, so the same settings split.
    let mut plain = build(&[(1, N, U)]);
    plain.settings.sample = input.settings.sample.clone();
    let split = ImDataConverter::split_by_faims_cv(&mut plain).unwrap();
    assert_eq!(split.groups.len(), 1);
}

#[test]
fn ceiling_and_constants() {
    assert_eq!(ImDataConverter::MAX_SPECTRA, FaimsHelper::MAX_SPECTRA);
    assert_eq!(ImDataConverter::MAX_SPECTRA, 100_000_000);
    assert_eq!(
        ImDataConverter::UNEXPECTED_COMPENSATION_VOLTAGE_WARNING,
        "Encountered spectrum with unexpected FAIMS CV (not in detected set): "
    );
}

#[test]
fn the_log_stream_emulation_reproduces_its_repetition_rule() {
    // Checked against the oracle above (cv_less_ms1_keeps_context prints the
    // line once and then "occurred 2 times"); these pin the eviction rule.
    let lines = |texts: &[&str]| {
        texts
            .iter()
            .map(|text| (*text).to_owned())
            .collect::<Vec<_>>()
    };
    assert_eq!(
        emulate_log_stream(&lines(&["a", "a"])),
        lines(&["a", "<a> occurred 2 times"])
    );
    assert_eq!(
        emulate_log_stream(&lines(&["a", "a", "b", "c"])),
        lines(&["a", "b", "<a> occurred 2 times", "c"])
    );
    assert_eq!(
        emulate_log_stream(&lines(&["b", "a", "b"])),
        lines(&["b", "a", "<b> occurred 2 times"])
    );
}
