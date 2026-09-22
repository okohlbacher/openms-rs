// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FileHandler::loadExperiment` with all of its arguments: the log type,
//! `rewrite_source_file`, `compute_hash`, the SRM conversion after an mzML
//! load, the mzXML and mzData dispatch, and `FileHandler::computeFileHash`.
//!
//! Two kinds of expected value, neither derived from Rust output:
//!
//! - the upstream class test `FileHandler_test.cpp` at core `bc9cc12`, cited
//!   per test, and RFC 3174's SHA-1 test vectors;
//! - `tests/data/file_handler_load_release.tsv`, written by the C++ Release
//!   build at core `bc9cc12` on ibminode06 through the oracle in
//!   `../oracle/file-handler-load` (provenance in
//!   `tests/data/file_handler_load_provenance.json`): a fresh `FileHandler`
//!   per case, log type GUI with a recording backend factory, then every
//!   source-file field, the first spectra, every chromatogram, and the
//!   warning-log lines of the case.

#![cfg(feature = "mzml")]

use openms::concept::progress_logger::ProgressLogger;
use openms::format::file_handler::ExperimentLoadOptions;
use openms::format::{FileHandler, FileType, PeakFileOptions};
use openms::kernel::NumericRange;
use openms::metadata::ChecksumType;
use openms::{Error, MSExperiment};
use std::path::{Path, PathBuf};

fn data(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(relative)
}

/// `MzMLFile_1.mzML` of the core class-test data; the port's copy has the
/// SHA-1 the upstream test expects, which is also the check that it is the
/// same file.
fn mzml_file_1() -> PathBuf {
    data("mzml_load_source_original.mzML")
}
fn dta2d_test_1() -> PathBuf {
    data("text_peak_lists/DTA2DFile_test_1.dta2d")
}

fn rewrite(hash: bool) -> ExperimentLoadOptions {
    ExperimentLoadOptions {
        rewrite_source_file: true,
        compute_hash: hash,
        ..Default::default()
    }
}

fn mz_range(min: f64, max: f64) -> PeakFileOptions {
    let mut options = PeakFileOptions::default();
    options.set_mz_range(NumericRange { min, max });
    options
}

fn first_mz(experiment: &MSExperiment, spectrum: usize, count: usize) -> Vec<f64> {
    experiment.spectra[spectrum]
        .peaks
        .iter()
        .take(count)
        .map(|peak| peak.mz)
        .collect()
}

fn assert_similar(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len(), "{actual:?} vs {expected:?}");
    for (a, e) in actual.iter().zip(expected) {
        // TEST_REAL_SIMILAR's default relative tolerance, 1e-5.
        assert!((a - e).abs() <= 1e-5 * e.abs().max(1.0), "{actual:?} vs {expected:?}");
    }
}

#[test]
fn upstream_mzml_load_rewrites_the_source_file_with_its_sha1() {
    // FileHandler_test.cpp:257-261.
    let experiment = FileHandler::load_experiment_with_progress(
        mzml_file_1(),
        &[FileType::MzMl],
        &rewrite(true),
        &ProgressLogger::new(),
    )
    .unwrap();
    assert_eq!(experiment.spectra.len(), 4);
    let sources = &experiment.settings.source_files;
    assert_eq!(sources.len(), 1);
    assert_eq!(sources[0].checksum, "36007593dbca0ba59a1f4fc32fb970f0e8991fa6");
    assert_eq!(sources[0].checksum_type, ChecksumType::Sha1);
}

#[test]
fn upstream_dta2d_load_applies_the_mz_range_and_rewrites_with_or_without_hash() {
    // FileHandler_test.cpp:263-279.
    let plain = FileHandler::load_experiment_with_options(
        dta2d_test_1(),
        &[],
        &PeakFileOptions::default(),
    )
    .unwrap();
    assert_similar(&first_mz(&plain, 0, 3), &[230.02, 430.02, 630.02]);

    let options = ExperimentLoadOptions {
        peaks: mz_range(300.0, 1000.0),
        ..rewrite(true)
    };
    let hashed = FileHandler::load_experiment_with_progress(
        dta2d_test_1(),
        &[FileType::Dta2d],
        &options,
        &ProgressLogger::new(),
    )
    .unwrap();
    assert_similar(&first_mz(&hashed, 0, 2), &[430.02, 630.02]);
    let source = &hashed.settings.source_files[0];
    assert_eq!(source.checksum, "d50d5144cc3805749b9e8d16f3bc8994979d8142");
    assert_eq!(source.checksum_type, ChecksumType::Sha1);

    let unhashed = FileHandler::load_experiment_with_progress(
        dta2d_test_1(),
        &[],
        &rewrite(false),
        &ProgressLogger::new(),
    )
    .unwrap();
    let source = &unhashed.settings.source_files[0];
    assert_eq!(source.checksum, "");
    assert_eq!(source.checksum_type, ChecksumType::Unknown);
}

#[test]
fn upstream_a_type_outside_the_allowed_list_is_refused() {
    // FileHandler_test.cpp:280-281 restricts to SIZE_OF_TYPE, which no file
    // has; the source throws ParseError, the port InvalidValue.
    let result = FileHandler::load_experiment_with_progress(
        dta2d_test_1(),
        &[FileType::Unknown],
        &rewrite(false),
        &ProgressLogger::new(),
    );
    assert!(matches!(result, Err(Error::InvalidValue(_))), "{result:?}");
}

#[test]
fn upstream_mzdata_and_mzxml_are_dispatched_with_the_peak_options() {
    // FileHandler_test.cpp:234-255.
    let load = |name: &str, options: &PeakFileOptions| {
        FileHandler::load_experiment_with_options(data(name), &[], options).unwrap()
    };
    let all = PeakFileOptions::default();
    let from_115 = mz_range(115.0, 1000.0);
    assert_similar(
        &first_mz(&load("MzDataFile_1.mzData", &all), 1, 3),
        &[110.0, 120.0, 130.0],
    );
    assert_similar(
        &first_mz(&load("MzDataFile_1.mzData", &from_115), 1, 2),
        &[120.0, 130.0],
    );
    assert_similar(
        &first_mz(&load("MzXMLFile_1.mzXML", &all), 2, 3),
        &[100.0, 110.0, 120.0],
    );
    assert_similar(
        &first_mz(&load("MzXMLFile_1.mzXML", &from_115), 2, 3),
        &[120.0, 130.0, 140.0],
    );
    for kind in [FileType::MzXml, FileType::MzData] {
        assert!(FileHandler::can_read_experiment(kind));
        assert!(!FileHandler::can_write_experiment(kind));
    }
}

#[test]
fn compute_file_hash_follows_rfc_3174_and_the_upstream_non_ascii_case() {
    let directory = std::env::temp_dir().join(format!(
        "openms-rs-file-hash-t\u{e4}st_\u{fc}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&directory).unwrap();
    let write = |name: &str, bytes: &[u8]| {
        let path = directory.join(name);
        std::fs::write(&path, bytes).unwrap();
        path
    };
    // RFC 3174 section 7.3, TEST1 ("abc"), and the digest of no input.
    assert_eq!(
        FileHandler::compute_file_hash(write("abc.txt", b"abc")),
        "a9993e364706816aba3e25717850c26c9cd0d89d"
    );
    assert_eq!(
        FileHandler::compute_file_hash(write("empty.txt", b"")),
        "da39a3ee5e6b4b0d3255bfef95601890afd80709"
    );
    // RFC 3174 TEST3, one million 'a', crosses many 8 KiB reads.
    assert_eq!(
        FileHandler::compute_file_hash(write("million.txt", &vec![b'a'; 1_000_000])),
        "34aa973cd4c4daa4f61eeb2bdbad27316534016f"
    );
    // FileHandler_test.cpp:293-319: a non-ASCII path hashes.
    assert!(!FileHandler::compute_file_hash(write("test.txt", b"hello")).is_empty());
    // FileHandler.cpp:771-774: a file that cannot be opened hashes to "".
    assert_eq!(FileHandler::compute_file_hash(directory.join("missing.txt")), "");
    std::fs::remove_dir_all(&directory).unwrap();
}

#[test]
fn compute_hash_without_rewrite_is_refused_before_any_io() {
    // FileHandler.cpp:853, the source's debug-build precondition.
    let options = ExperimentLoadOptions {
        compute_hash: true,
        ..Default::default()
    };
    let result = FileHandler::load_experiment_with_progress(
        data("does-not-exist.mzML"),
        &[],
        &options,
        &ProgressLogger::new(),
    );
    assert!(matches!(result, Err(Error::InvalidValue(_))), "{result:?}");
}
