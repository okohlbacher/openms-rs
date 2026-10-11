// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `src/format/ms_data_cached_consumer.rs`, the port of
//! `FORMAT/DATAACCESS/MSDataCachedConsumer.h`.
//!
//! Every `START_SECTION` of the upstream `MSDataCachedConsumer_test.cpp` (core
//! `bc9cc12`) is reproduced on the class-test file `MzMLFile_1.mzML`, which this
//! repository holds byte-identical as `tests/data/mzml_validator/MzMLFile_1.mzML`
//! (tier 3). `release_build_bytes` compares the consumer's file with the cache
//! the Linux x86_64 Release build wrote for the same input
//! (`tests/data/f4_progress/MzMLFile_1.cached` and `empty.cached`, from
//! `../oracle/F4-progress` on ibminode06, install
//! `openms4-release-bc9cc12-c19e494-174b576`): `CachedMzMLHandler::writeMemdump`
//! (`CachedMzMLHandler.cpp:33-58`) writes the identifier, every spectrum with
//! `writeSpectrum_`, every chromatogram with `writeChromatogram_` and the two
//! counts, which is what the consumer's constructor, `consumeSpectrum`,
//! `consumeChromatogram` and destructor write in turn
//! (`MSDataCachedConsumer.cpp:18-75`), so the bytes must agree (tier 1 bytes,
//! tier 4 equivalence). The `release_consumer_*` tests below compare with the
//! Release build's own `MSDataCachedConsumer` directly (tier 1,
//! `../oracle/integ-F-cachedconsumer`).
#![cfg(feature = "mzml")]

use openms::Error;
use openms::format::cached_mzml_handler::CachedMzMLHandler;
use openms::format::ms_data_cached_consumer::MSDataCachedConsumer;
use openms::interfaces::MSDataConsumer;
use openms::kernel::{DataArray, MSExperiment};
use openms::metadata::ExperimentalSettings;
use openms::system::file::TempDir;
use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};
use std::ops::ControlFlow;
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn load_exp() -> MSExperiment {
    openms::format::mzml::load(data("mzml_validator/MzMLFile_1.mzML")).unwrap()
}

/// Every consumed record answers `Continue`: this consumer never asks its
/// producer to stop.
fn go(result: openms::Result<ControlFlow<()>>) {
    assert_eq!(result.unwrap(), ControlFlow::Continue(()));
}

// START_SECTION((MSDataCachedConsumer()))
// START_SECTION((~MSDataCachedConsumer()))
#[test]
fn construction_and_destruction() {
    let dir = TempDir::new(false).unwrap();
    let tmp_filename = dir.path().join("constructed.cached");
    let cached_consumer = MSDataCachedConsumer::new(&tmp_filename).unwrap();
    assert!(cached_consumer.clear_data());
    assert_eq!(cached_consumer.filename(), tmp_filename.as_path());
    drop(cached_consumer);
    // The destructor wrote the footer: the file is a valid, empty cache.
    let mut cache = CachedMzMLHandler::new();
    cache.create_memdump_index(&tmp_filename).unwrap();
    assert!(cache.spectra_index().is_empty());
}

// START_SECTION((void consumeSpectrum(SpectrumType & s)))
#[test]
fn consume_spectrum() {
    let dir = TempDir::new(false).unwrap();
    let tmp_filename = dir.path().join("spectra.cached");
    let mut cached_consumer = MSDataCachedConsumer::with_clear_data(&tmp_filename, false).unwrap();

    let mut exp = load_exp();
    assert!(!exp.spectra.is_empty());

    cached_consumer.set_expected_size(2, 0).unwrap();
    assert_eq!(
        cached_consumer
            .consume_spectrum(&mut exp.spectra[0])
            .unwrap(),
        ControlFlow::Continue(())
    );
    go(cached_consumer.consume_spectrum(&mut exp.spectra[1]));
    assert_eq!(cached_consumer.spectra_written(), 2);
    cached_consumer.finish().unwrap();

    // Check whether it was written to disk correctly...
    let mut cache = CachedMzMLHandler::new();
    cache.create_memdump_index(&tmp_filename).unwrap();
    let spectra_index = cache.spectra_index().to_vec();
    let mut ifs = BufReader::new(File::open(&tmp_filename).unwrap());

    // retrieve the spectrum
    ifs.seek(SeekFrom::Start(spectra_index[0])).unwrap();
    let (data, _ms_level, _rt) = CachedMzMLHandler::read_spectrum_fast(&mut ifs).unwrap();
    assert_eq!(data[0].data.len(), exp.spectra[0].peaks.len());
    assert_eq!(data[1].data.len(), exp.spectra[0].peaks.len());

    // retrieve the spectrum
    ifs.seek(SeekFrom::Start(spectra_index[1])).unwrap();
    let (data, _ms_level, _rt) = CachedMzMLHandler::read_spectrum_fast(&mut ifs).unwrap();
    assert_eq!(data[0].data.len(), exp.spectra[1].peaks.len());
    assert_eq!(data[1].data.len(), exp.spectra[1].peaks.len());
}

// START_SECTION((void consumeChromatogram(ChromatogramType & c)))
#[test]
fn consume_chromatogram() {
    let dir = TempDir::new(false).unwrap();
    let tmp_filename = dir.path().join("chromatograms.cached");
    let mut cached_consumer = MSDataCachedConsumer::with_clear_data(&tmp_filename, false).unwrap();

    let mut exp = load_exp();
    assert!(!exp.chromatograms.is_empty());

    cached_consumer.set_expected_size(0, 1).unwrap();
    go(cached_consumer.consume_chromatogram(&mut exp.chromatograms[0]));
    assert_eq!(cached_consumer.chromatograms_written(), 1);
    drop(cached_consumer);

    // Check whether it was written to disk correctly...
    let mut cache = CachedMzMLHandler::new();
    cache.create_memdump_index(&tmp_filename).unwrap();
    let chrom_index = cache.chromatogram_index().to_vec();
    let mut ifs = BufReader::new(File::open(&tmp_filename).unwrap());

    // retrieve the chromatogram
    ifs.seek(SeekFrom::Start(chrom_index[0])).unwrap();
    let data = CachedMzMLHandler::read_chromatogram_fast(&mut ifs).unwrap();
    assert_eq!(data[0].data.len(), exp.chromatograms[0].peaks.len());
    assert_eq!(data[1].data.len(), exp.chromatograms[0].peaks.len());
}

// START_SECTION((MSDataCachedConsumer(std::string filename, bool clearData=true)))
#[test]
fn constructor_clear_data() {
    {
        let dir = TempDir::new(false).unwrap();
        let tmp_filename = dir.path().join("cleared.cached");
        let mut cached_consumer =
            MSDataCachedConsumer::with_clear_data(&tmp_filename, true).unwrap();

        let mut exp = load_exp();
        assert!(!exp.spectra.is_empty());
        let first_spectrum = exp.spectra[0].clone();

        cached_consumer.set_expected_size(2, 0).unwrap();

        assert!(!exp.spectra[0].peaks.is_empty());

        go(cached_consumer.consume_spectrum(&mut exp.spectra[0]));

        assert_eq!(exp.spectra[0].peaks.len(), 0);
        assert!(exp.spectra[0] != first_spectrum);
    }
    {
        let dir = TempDir::new(false).unwrap();
        let tmp_filename = dir.path().join("kept.cached");
        let mut cached_consumer =
            MSDataCachedConsumer::with_clear_data(&tmp_filename, false).unwrap();

        let mut exp = load_exp();
        assert!(!exp.spectra.is_empty());
        let first_spectrum = exp.spectra[0].clone();

        cached_consumer.set_expected_size(2, 0).unwrap();

        assert!(!exp.spectra[0].peaks.is_empty());

        go(cached_consumer.consume_spectrum(&mut exp.spectra[0]));

        assert!(!exp.spectra[0].peaks.is_empty());
        assert!(exp.spectra[0] == first_spectrum);
    }
}

// START_SECTION((void setExpectedSize(Size expectedSpectra, Size expectedChromatograms)))
// NOT_TESTABLE upstream ("tested above"); it is called in every section above.

// START_SECTION([EXTRA] test empty file)
#[test]
fn empty_file() {
    // try an empty file
    let dir = TempDir::new(false).unwrap();
    let tmp_filename = dir.path().join("empty.cached");
    let cached_consumer = MSDataCachedConsumer::with_clear_data(&tmp_filename, false).unwrap();
    drop(cached_consumer);

    // Check whether it was written to disk correctly...
    let mut cache = CachedMzMLHandler::new();
    cache.create_memdump_index(&tmp_filename).unwrap();
    assert_eq!(cache.spectra_index().len(), 0);
    assert_eq!(cache.chromatogram_index().len(), 0);
}

// START_SECTION((void setExperimentalSettings(const ExperimentalSettings&)))
#[test]
fn set_experimental_settings() {
    let dir = TempDir::new(false).unwrap();
    let tmp_filename = dir.path().join("settings.cached");
    let mut cached_consumer = MSDataCachedConsumer::with_clear_data(&tmp_filename, true).unwrap();

    cached_consumer.set_expected_size(2, 0).unwrap();
    let s = ExperimentalSettings::default();
    cached_consumer.set_experimental_settings(&s).unwrap();
    cached_consumer.finish().unwrap();
}

/// The consumer writes the bytes the Release build's `writeMemdump` wrote for
/// the same experiment, and a consumer given nothing writes the Release
/// build's empty cache.
#[test]
fn release_build_bytes() {
    let dir = TempDir::new(false).unwrap();
    let path = dir.path().join("MzMLFile_1.cached");
    let mut exp = load_exp();
    let mut cached_consumer = MSDataCachedConsumer::new(&path).unwrap();
    for spectrum in &mut exp.spectra {
        go(cached_consumer.consume_spectrum(spectrum));
    }
    for chromatogram in &mut exp.chromatograms {
        go(cached_consumer.consume_chromatogram(chromatogram));
    }
    cached_consumer.finish().unwrap();
    let expected = std::fs::read(data("f4_progress/MzMLFile_1.cached")).unwrap();
    assert!(std::fs::read(&path).unwrap() == expected);
    // clearData (the default) kept the metadata and removed the data.
    assert!(exp.spectra.iter().all(|s| s.peaks.is_empty()));
    assert!(exp.chromatograms.iter().all(|c| c.peaks.is_empty()));
    assert_eq!(exp.spectra[0].native_id, load_exp().spectra[0].native_id);

    let empty = dir.path().join("empty.cached");
    MSDataCachedConsumer::new(&empty).unwrap().finish().unwrap();
    let expected = std::fs::read(data("f4_progress/empty.cached")).unwrap();
    assert!(std::fs::read(&empty).unwrap() == expected);
}

/// `consumeSpectrum` after a chromatogram throws `IllegalArgument` with this
/// message (`MSDataCachedConsumer.cpp:42-46`) before writing anything.
#[test]
fn spectrum_after_chromatogram_is_refused() {
    let dir = TempDir::new(false).unwrap();
    let path = dir.path().join("order.cached");
    let mut exp = load_exp();
    let mut cached_consumer = MSDataCachedConsumer::new(&path).unwrap();
    go(cached_consumer.consume_chromatogram(&mut exp.chromatograms[0]));
    let before = exp.spectra[0].clone();
    match cached_consumer.consume_spectrum(&mut exp.spectra[0]) {
        Err(Error::InvalidValue(message)) => {
            assert_eq!(message, "Cannot write spectra after writing chromatograms.")
        }
        other => panic!("expected the source's refusal, got {other:?}"),
    }
    assert!(exp.spectra[0] == before);
    assert_eq!(cached_consumer.spectra_written(), 0);
    cached_consumer.finish().unwrap();
    let mut cache = CachedMzMLHandler::new();
    cache.create_memdump_index(&path).unwrap();
    assert_eq!(cache.spectra_index().len(), 0);
    assert_eq!(cache.chromatogram_index().len(), 1);
}

/// Float and integer data arrays are cached (as doubles, read back as named
/// arrays) and, with clearData, cleared from the caller's record together
/// with the peaks, as the header's note says.
#[test]
fn data_arrays_are_written_and_cleared() {
    let dir = TempDir::new(false).unwrap();
    let path = dir.path().join("arrays.cached");
    let mut exp = load_exp();
    let n = exp.spectra[0].peaks.len();
    let floats = DataArray::<f32> {
        name: "FWHM".into(),
        data: vec![0.5; n],
        ..DataArray::default()
    };
    let integers = DataArray::<i32> {
        name: "charge".into(),
        data: vec![2; n],
        ..DataArray::default()
    };
    exp.spectra[0].float_data_arrays = vec![floats];
    exp.spectra[0].integer_data_arrays = vec![integers];
    let mut cached_consumer = MSDataCachedConsumer::new(&path).unwrap();
    go(cached_consumer.consume_spectrum(&mut exp.spectra[0]));
    assert!(exp.spectra[0].float_data_arrays.is_empty());
    assert!(exp.spectra[0].integer_data_arrays.is_empty());
    cached_consumer.finish().unwrap();

    let mut cache = CachedMzMLHandler::new();
    cache.create_memdump_index(&path).unwrap();
    let mut ifs = BufReader::new(File::open(&path).unwrap());
    ifs.seek(SeekFrom::Start(cache.spectra_index()[0])).unwrap();
    let (data, _, _) = CachedMzMLHandler::read_spectrum_fast(&mut ifs).unwrap();
    assert_eq!(data.len(), 4);
    assert_eq!(data[2].description, "FWHM");
    assert_eq!(data[2].data, vec![0.5; n]);
    assert_eq!(data[3].description, "charge");
    assert_eq!(data[3].data, vec![2.0; n]);
}

/// A file that cannot be created is an error here; the source's `ofstream`
/// fails silently and every later write is lost.
#[test]
fn unwritable_path_is_an_error() {
    let dir = TempDir::new(false).unwrap();
    let missing = dir.path().join("no-such-directory").join("x.cached");
    assert!(matches!(
        MSDataCachedConsumer::new(&missing),
        Err(Error::Io(_))
    ));
}

// ---------------------------------------------------------------------------
// Executed differential against the Linux x86_64 Release build
// ---------------------------------------------------------------------------
//
// `tests/data/ms_data_cached_consumer/` holds the files and records the
// oracle driver `../oracle/integ-F-cachedconsumer` (outside the repository)
// wrote with the Release build's own `MSDataCachedConsumer`, linking the
// Linux x86_64 Release install `openms4-release-bc9cc12-c19e494-174b576` on
// ibminode06, run twice with identical files and records.

fn release(name: &str) -> Vec<u8> {
    std::fs::read(data(&format!("ms_data_cached_consumer/{name}"))).unwrap()
}

fn release_record(key: &str) -> String {
    let text = String::from_utf8(release("cpp_release_records.tsv")).unwrap();
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("{key}\t")).map(str::to_owned))
        .unwrap_or_else(|| panic!("no record {key}"))
}

/// Every record streamed with `clearData` (the default), spectra then
/// chromatograms: the same bytes, and the same cleared records (peaks and
/// all data arrays gone, native IDs kept).
#[test]
fn release_consumer_all_records_cleared() {
    let dir = TempDir::new(false).unwrap();
    let path = dir.path().join("all_clear.cached");
    let mut exp = load_exp();
    assert_eq!(release_record("spectra"), exp.spectra.len().to_string());
    assert_eq!(
        release_record("chromatograms"),
        exp.chromatograms.len().to_string()
    );
    let mut cached_consumer = MSDataCachedConsumer::new(&path).unwrap();
    cached_consumer
        .set_expected_size(exp.spectra.len(), exp.chromatograms.len())
        .unwrap();
    for spectrum in &mut exp.spectra {
        go(cached_consumer.consume_spectrum(spectrum));
    }
    for chromatogram in &mut exp.chromatograms {
        go(cached_consumer.consume_chromatogram(chromatogram));
    }
    drop(cached_consumer);
    assert!(std::fs::read(&path).unwrap() == release("all_clear.cached"));
    for (i, s) in exp.spectra.iter().enumerate() {
        assert_eq!(
            format!(
                "{}|{}|{}|{}|{}",
                s.peaks.len(),
                s.float_data_arrays.len(),
                s.integer_data_arrays.len(),
                s.string_data_arrays.len(),
                s.native_id
            ),
            release_record(&format!("all_clear.spectrum{i}"))
        );
    }
    for (i, c) in exp.chromatograms.iter().enumerate() {
        assert_eq!(
            format!("{}|{}", c.peaks.len(), c.native_id),
            release_record(&format!("all_clear.chromatogram{i}"))
        );
    }
}

/// Spectra 0 and 1 without `clearData`: the same bytes, the records
/// unchanged; and a consumer given nothing.
#[test]
fn release_consumer_kept_records_and_empty_file() {
    let dir = TempDir::new(false).unwrap();
    let path = dir.path().join("spectra_keep.cached");
    let original = load_exp();
    let mut exp = original.clone();
    let mut cached_consumer = MSDataCachedConsumer::with_clear_data(&path, false).unwrap();
    go(cached_consumer.consume_spectrum(&mut exp.spectra[0]));
    go(cached_consumer.consume_spectrum(&mut exp.spectra[1]));
    cached_consumer.finish().unwrap();
    assert!(std::fs::read(&path).unwrap() == release("spectra_keep.cached"));
    let unchanged = format!(
        "{}{}",
        u8::from(exp.spectra[0] == original.spectra[0]),
        u8::from(exp.spectra[1] == original.spectra[1])
    );
    assert_eq!(unchanged, release_record("spectra_keep.unchanged"));

    let empty = dir.path().join("empty.cached");
    MSDataCachedConsumer::with_clear_data(&empty, false)
        .unwrap()
        .finish()
        .unwrap();
    assert!(std::fs::read(&empty).unwrap() == release("empty.cached"));
}

/// A spectrum after a chromatogram: the Release build throws
/// `IllegalArgument` with this message, leaves the spectrum unchanged and
/// writes a file holding the one chromatogram.
#[test]
fn release_consumer_refuses_a_spectrum_after_a_chromatogram() {
    let dir = TempDir::new(false).unwrap();
    let path = dir.path().join("order.cached");
    let original = load_exp();
    let mut exp = original.clone();
    let mut cached_consumer = MSDataCachedConsumer::new(&path).unwrap();
    go(cached_consumer.consume_chromatogram(&mut exp.chromatograms[0]));
    let message = match cached_consumer.consume_spectrum(&mut exp.spectra[0]) {
        Err(Error::InvalidValue(message)) => message,
        other => panic!("expected the refusal, got {other:?}"),
    };
    assert_eq!(
        format!("IllegalArgument: {message}"),
        release_record("order.exception")
    );
    assert_eq!(
        u8::from(exp.spectra[0] == original.spectra[0]).to_string(),
        release_record("order.spectrum_unchanged")
    );
    drop(cached_consumer);
    assert!(std::fs::read(&path).unwrap() == release("order.cached"));
}
