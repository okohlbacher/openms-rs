// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `src/format/chromeleon.rs`, the port of `FORMAT/ChromeleonFile.h`.
//!
//! Every `START_SECTION` of the upstream `ChromeleonFile_test.cpp` (core
//! `bc9cc12`) is reproduced with its literals. The fixtures are the upstream
//! class-test data files at the same revision, recorded with their sha256 in
//! `tests/data/chromeleon_file_provenance.json`. The mzML round trips of the
//! first two sections need the `mzml` feature and are gated on it.

use openms::Error;
use openms::format::chromeleon::{ChromeleonFile, ReadLimits};
use openms::kernel::MSExperiment;
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/chromeleon_file")
        .join(name)
}

/// The source's `TEST_REAL_SIMILAR` with its default tolerances: an absolute
/// difference of at most 1e-5 or a ratio within 1 +- 1e-5.
fn similar(a: f64, b: f64) -> bool {
    if (a - b).abs() <= 1e-5 {
        return true;
    }
    if a == 0.0 || b == 0.0 {
        return false;
    }
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    (1.0 / (1.0 + 1e-5)..=1.0 + 1e-5).contains(&ratio)
}

fn meta(experiment: &MSExperiment, key: &str) -> String {
    experiment.settings.metadata[key].to_string()
}

fn assert_points(experiment: &MSExperiment, expected: &[(usize, f64, f64)]) {
    let c = &experiment.chromatograms[0];
    for &(index, rt, intensity) in expected {
        assert!(
            similar(c.peaks[index].rt, rt),
            "rt[{index}] = {}",
            c.peaks[index].rt
        );
        assert!(
            similar(f64::from(c.peaks[index].intensity), intensity),
            "intensity[{index}] = {}",
            c.peaks[index].intensity
        );
    }
}

#[cfg(feature = "mzml")]
fn assert_mzml_round_trip(experiment: &MSExperiment, indices: &[usize]) {
    let dir = openms::system::file::TempDir::new(false).unwrap();
    let path = dir.path().join("chromeleon.mzML");
    // The source's mzML writer drops an instrument name that is no PSI-MS
    // instrument-model term and writes the generic `MS:1000031` instead
    // (`MzMLHandler.cpp:3961-3969`); this crate's writer refuses such a name
    // rather than lose it. The source section asserts only the chromatogram,
    // so the round trip stores the experiment without the name the source
    // would have dropped.
    let mut stored = experiment.clone();
    assert!(openms::format::mzml::store(&path, &stored).is_err());
    stored.settings.instrument.name.clear();
    openms::format::mzml::store(&path, &stored).unwrap();
    let read_exp = openms::format::mzml::load(&path).unwrap();
    assert_eq!(read_exp.chromatograms.len(), 1);
    let c1 = &experiment.chromatograms[0];
    let c2 = &read_exp.chromatograms[0];
    assert_eq!(c1.peaks.len(), c2.peaks.len());
    for &i in indices {
        assert!(similar(c1.peaks[i].rt, c2.peaks[i].rt));
        assert!(similar(
            f64::from(c1.peaks[i].intensity),
            f64::from(c2.peaks[i].intensity)
        ));
    }
}

// START_SECTION(ChromeleonFile()) / START_SECTION(~ChromeleonFile())
#[test]
fn construction_and_destruction() {
    let file = ChromeleonFile::new();
    assert_eq!(file, ChromeleonFile);
    let _ = file;
}

// START_SECTION(void load(const std::string& filename, MSExperiment& experiment) const)
#[test]
fn load_matches_the_source_section() {
    let experiment = ChromeleonFile::new()
        .load(data("20171013_HMP_C61_ISO_P1_GA1_UV_VIS_2.txt"))
        .unwrap();
    assert_eq!(meta(&experiment, "acq_method_name"), "UV_VIS_2");
    assert_eq!(meta(&experiment, "mzml_id"), "20171013_C61_ISO_P1_GA1");
    assert_eq!(
        experiment.settings.instrument.name,
        "HM_metode_ZorBax_0,02%_Acetic_acid_ver6"
    );
    assert_eq!(
        experiment.settings.instrument.software.name,
        "New ProcMethod"
    );
    assert_eq!(meta(&experiment, "injection_date"), "10/13/2017");
    assert_eq!(meta(&experiment, "injection_time"), "6:28:26 PM");
    assert_eq!(meta(&experiment, "detector"), "UV");
    assert_eq!(meta(&experiment, "signal_quantity"), "Absorbance");
    assert_eq!(meta(&experiment, "signal_unit"), "mAU");
    assert_eq!(meta(&experiment, "signal_info"), "WVL:280 nm");
    assert_eq!(experiment.chromatograms.len(), 1);
    assert_eq!(experiment.chromatograms[0].peaks.len(), 3301);
    assert_points(
        &experiment,
        &[
            (0, 0.0, 0.0),
            (660, 2.2, -0.812998),
            (1320, 4.4, -0.791189),
            (1980, 6.6, -0.285533),
            (2640, 8.8, -0.485941),
            (3300, 11.0, -0.130904),
        ],
    );
    #[cfg(feature = "mzml")]
    assert_mzml_round_trip(&experiment, &[0, 660, 1320, 1980, 2640, 3300]);
}

// START_SECTION(load_with_new_raw_data_header)
#[test]
fn load_with_new_raw_data_header() {
    let experiment = ChromeleonFile::new()
        .load(data("ChromeleonFile_new_header.txt"))
        .unwrap();
    assert_eq!(meta(&experiment, "acq_method_name"), "RID_Signal");
    assert_eq!(meta(&experiment, "mzml_id"), "S1");
    assert_eq!(experiment.settings.instrument.name, "SUGARS_MP.M");
    assert_eq!(experiment.settings.instrument.software.name, "SUGARS_CAL");
    assert_eq!(meta(&experiment, "injection_date"), "13/06/2019");
    assert_eq!(meta(&experiment, "injection_time"), "12:11:41 AM");
    assert_eq!(meta(&experiment, "detector"), "LCSystem");
    assert_eq!(meta(&experiment, "signal_quantity"), "");
    assert_eq!(meta(&experiment, "signal_unit"), "nRIU");
    assert_eq!(meta(&experiment, "signal_info"), "");
    assert_eq!(experiment.chromatograms.len(), 1);
    assert_eq!(experiment.chromatograms[0].peaks.len(), 10);
    assert_points(
        &experiment,
        &[
            (0, 0.0, 5.060000),
            (2, 0.014430, 5.450000),
            (4, 0.028860, 5.580000),
            (6, 0.043290, 5.380000),
            (9, 0.064935, 4.930000),
        ],
    );
    #[cfg(feature = "mzml")]
    assert_mzml_round_trip(&experiment, &[0, 2, 4, 6, 9]);
}

// START_SECTION(load_file_with_comma_thousands_separator)
#[test]
fn load_file_with_comma_thousands_separator() {
    let experiment = ChromeleonFile::new()
        .load(data("ChromeleonFile_commas.txt"))
        .unwrap();
    assert_eq!(meta(&experiment, "acq_method_name"), "RID_Signal");
    assert_eq!(meta(&experiment, "mzml_id"), "S2");
    assert_eq!(experiment.settings.instrument.name, "SUGARS_MP.M");
    assert_eq!(experiment.settings.instrument.software.name, "SUGARS_CAL");
    assert_eq!(meta(&experiment, "injection_date"), "12/06/2019");
    assert_eq!(meta(&experiment, "injection_time"), "11:49:36 PM");
    assert_eq!(meta(&experiment, "detector"), "LCSystem");
    assert_eq!(meta(&experiment, "signal_quantity"), "");
    assert_eq!(meta(&experiment, "signal_unit"), "nRIU");
    assert_eq!(meta(&experiment, "signal_info"), "");
    assert_eq!(experiment.chromatograms.len(), 1);
    assert_eq!(experiment.chromatograms[0].peaks.len(), 8);
    assert_points(
        &experiment,
        &[
            (0, 0.0, 1.4),
            (1, 8.300000, 1.6),
            (2, 8.513709, -18.980000),
            (3, 8.520924, -1234567.890000),
            (4, 9.855700, 1946.610000),
            (5, 9.884560, 2067.450000),
            (6, 9.898991, 2345678.900000),
            (7, 9.920635, 2028.580000),
        ],
    );
}

// ----------------------------------------------------------------------------
// Native regressions derived from ChromeleonFile.cpp
// ----------------------------------------------------------------------------

fn read(text: &str) -> openms::Result<MSExperiment> {
    ChromeleonFile::new().read(text.as_bytes())
}

#[test]
fn remove_commas_and_parse_double_edits_its_argument() {
    let file = ChromeleonFile::new();
    let mut number = String::from("-1,234,567.89");
    assert_eq!(
        file.remove_commas_and_parse_double(&mut number).unwrap(),
        -1_234_567.89
    );
    assert_eq!(number, "-1234567.89");
    let mut bad = String::from("1,2x");
    assert!(matches!(
        file.remove_commas_and_parse_double(&mut bad),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(bad, "12x");
}

#[test]
fn no_data_marker_gives_one_empty_chromatogram() {
    let experiment = read("Channel\tA\nDetector\tUV\n").unwrap();
    assert_eq!(experiment.chromatograms.len(), 1);
    assert!(experiment.chromatograms[0].peaks.is_empty());
    assert_eq!(meta(&experiment, "acq_method_name"), "A");
}

#[test]
fn chromatogram_data_marker_and_caption_line() {
    // The caption after the marker is skipped even when it looks like data.
    let experiment = read("Chromatogram Data:\n1\t2\t3\n4\t5\t6\n\n").unwrap();
    assert_eq!(experiment.chromatograms[0].peaks.len(), 1);
    assert_eq!(experiment.chromatograms[0].peaks[0].rt, 4.0);
    assert_eq!(experiment.chromatograms[0].peaks[0].intensity, 6.0);
}

#[test]
fn markers_and_keys_must_match_the_whole_line() {
    // `Raw Data: x` is not the marker; `Channel ` with a space is not a key.
    let experiment = read("Channel A\nRaw Data: x\n").unwrap();
    assert!(experiment.settings.metadata.is_empty());
    assert!(experiment.chromatograms[0].peaks.is_empty());
    // A later key overwrites an earlier one.
    let experiment = read("Channel\tA\nChannel\tB\n").unwrap();
    assert_eq!(meta(&experiment, "acq_method_name"), "B");
    // A form feed is a line separator that `.` does not match under no_mod_s.
    let experiment = read("Detector\tU\x0cV\n").unwrap();
    assert!(experiment.settings.metadata.is_empty());
}

#[test]
fn malformed_data_lines_are_refused() {
    match read("Raw Data:\ncaption\n1\t2\n") {
        Err(Error::Parse { line, message }) => {
            assert_eq!(line, 3);
            assert!(
                message.contains("Couldn't parse the raw data."),
                "{message}"
            );
        }
        other => panic!("{other:?}"),
    }
    assert!(matches!(
        read("Raw Data:\ncaption\n1\t2\tx\n"),
        Err(Error::Parse { line: 3, .. })
    ));
    assert!(matches!(
        read("Raw Data:\ncaption\n1\t2\tinf\n"),
        Err(Error::Parse { line: 3, .. })
    ));
}

#[test]
fn point_limit_is_enforced() {
    let limits = ReadLimits {
        max_points: 1,
        ..ReadLimits::default()
    };
    assert!(matches!(
        ChromeleonFile::new()
            .read_with_limits("Raw Data:\nc\n1\t0\t1\n2\t0\t2\n".as_bytes(), limits),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn missing_file_is_not_found() {
    match ChromeleonFile::new().load(data("does_not_exist.txt")) {
        Err(Error::Io(error)) => assert_eq!(error.kind(), std::io::ErrorKind::NotFound),
        other => panic!("{other:?}"),
    }
}
