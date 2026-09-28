// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FORMAT/IndexedMzMLFileLoader`.
//!
//! Every section of the upstream class test
//! (`src/tests/class_tests/openms/source/IndexedMzMLFileLoader_test.cpp` at
//! core `bc9cc12`) is transcribed below with its line numbers: source review,
//! tier 3. The fixtures are the unmodified upstream `IndexedmzMLFile_1.mzML`
//! (`tests/data/indexed_mzml/`) and `MzMLFile_1.mzML`
//! (`tests/data/mzml_validator/`), both already recorded in
//! `tests/data/on_disc_experiment_provenance.json`. The upstream test's
//! `VALIDATE_TMP_FILES` (line 184) checks the written files against the mzML schema; that
//! step runs here only with the `mzml-schema` feature, which carries the XSD
//! validator.
//!
//! The `native_*` tests are independently derived (tier 4): they pin the
//! refusals this port adds where the source would dereference a null pointer
//! or silently drop a requested encoding.

#![cfg(feature = "mzml")]

use openms::Error;
use openms::format::indexed_mzml_file_loader::IndexedMzMLFileLoader;
use openms::format::mzml;
use openms::format::peak_options::{NumpressCompression, NumpressConfig, PeakFileOptions};
use openms::kernel::MSExperiment;
use openms::kernel::on_disc_experiment::OnDiscPeakMap;
use openms::system::file::TempDir;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// The unmodified upstream `IndexedmzMLFile_1.mzML`: 2 spectra, 1 chromatogram.
fn indexed() -> PathBuf {
    data("indexed_mzml/IndexedmzMLFile_1.mzML")
}

/// The unmodified upstream `MzMLFile_1.mzML`, which carries no index.
fn not_indexed() -> PathBuf {
    data("mzml_validator/MzMLFile_1.mzML")
}

/// `MzMLFile().load(...)` of the upstream sections.
fn mzml_file_load(path: &Path) -> MSExperiment {
    mzml::read(BufReader::new(File::open(path).unwrap())).unwrap()
}

/// The comparison the three load/store sections share
/// (`IndexedMzMLFileLoader_test.cpp:80-95`, `:128-143`, `:162-177`).
fn assert_matches_whole_file_load(exp: &mut OnDiscPeakMap, exp2: &MSExperiment) {
    assert_eq!(exp.spectrum_count(), exp2.spectra.len());
    assert_eq!(exp.chromatogram_count(), exp2.chromatograms.len());
    assert_eq!(exp.spectrum_count(), 2);
    assert_eq!(exp.chromatogram_count(), 1);
    assert_eq!(exp.spectrum(0).unwrap(), exp2.spectra[0]);
    for i in 0..exp.spectrum_count() {
        assert_eq!(exp.spectrum(i).unwrap(), exp2.spectra[i], "spectrum {i}");
    }
    for i in 0..exp.chromatogram_count() {
        assert_eq!(
            exp.chromatogram(i).unwrap(),
            exp2.chromatograms[i],
            "chromatogram {i}"
        );
    }
    assert_eq!(exp.experimental_settings().unwrap(), &exp2.settings);
}

/// The upstream `VALIDATE_TMP_FILES` step for one written file.
#[cfg(feature = "mzml-schema")]
fn validate_tmp_file(path: &Path) {
    let report = mzml::validate_schema(path).unwrap();
    assert!(report.is_valid(), "{:?}", report.diagnostics);
}

#[cfg(not(feature = "mzml-schema"))]
fn validate_tmp_file(_path: &Path) {}

#[test]
fn constructs() {
    // IndexedMzMLFileLoader_test.cpp:37-44: construction and destruction.
    let loader = IndexedMzMLFileLoader::new();
    drop(loader);
}

#[test]
fn get_options_const() {
    // IndexedMzMLFileLoader_test.cpp:46-50.
    let file = IndexedMzMLFileLoader::new();
    let options = file.options();
    assert!(!options.has_ms_levels());
}

#[test]
fn get_options_mutable() {
    // IndexedMzMLFileLoader_test.cpp:52-56.
    let mut file = IndexedMzMLFileLoader::new();
    file.options_mut().add_ms_level(1).unwrap();
    assert!(file.options().has_ms_levels());
}

#[test]
fn set_options() {
    // IndexedMzMLFileLoader_test.cpp:58-67.
    let mut file = IndexedMzMLFileLoader::new();
    let options = file.options();
    assert!(!options.has_ms_levels());
    assert!(!file.options().has_ms_levels());
    let mut new_options: PeakFileOptions = file.options().clone();
    new_options.add_ms_level(1).unwrap();
    file.set_options(new_options);
    assert!(file.options().has_ms_levels());
}

#[test]
fn load() {
    // IndexedMzMLFileLoader_test.cpp:71-97.
    let file = IndexedMzMLFileLoader::new();
    let mut exp = OnDiscPeakMap::new();
    file.load(indexed(), &mut exp).unwrap();
    let exp2 = mzml_file_load(&indexed());
    assert_matches_whole_file_load(&mut exp, &exp2);
}

#[test]
fn check_parsing() {
    // IndexedMzMLFileLoader_test.cpp:99-110 ([EXTRA]CheckParsing).
    let file = IndexedMzMLFileLoader::new();
    let mut exp = OnDiscPeakMap::new();
    let success = file.load(not_indexed(), &mut exp).unwrap();
    assert!(!success);
    let success = file.load(indexed(), &mut exp).unwrap();
    assert!(success);
}

#[test]
fn store_on_disc() {
    // IndexedMzMLFileLoader_test.cpp:112-145.
    let mut file = IndexedMzMLFileLoader::new();
    let mut exp = OnDiscPeakMap::new();
    let mut exp_ = OnDiscPeakMap::new();
    file.load(indexed(), &mut exp_).unwrap();
    let dir = TempDir::new(false).unwrap();
    let tmp_filename = dir.path().join("stored_on_disc.mzML");
    file.store_on_disc(&tmp_filename, &mut exp_).unwrap();

    let success = file.load(&tmp_filename, &mut exp).unwrap();
    assert!(success);

    let exp2 = mzml_file_load(&indexed());
    assert_matches_whole_file_load(&mut exp, &exp2);
    validate_tmp_file(&tmp_filename);
}

#[test]
fn store_in_memory() {
    // IndexedMzMLFileLoader_test.cpp:147-179.
    let mut file = IndexedMzMLFileLoader::new();
    let mut exp = OnDiscPeakMap::new();
    let exp2 = mzml_file_load(&indexed());
    let dir = TempDir::new(false).unwrap();
    let tmp_filename = dir.path().join("stored_in_memory.mzML");
    file.store(&tmp_filename, &exp2).unwrap();

    let success = file.load(&tmp_filename, &mut exp).unwrap();
    assert!(success);
    assert_matches_whole_file_load(&mut exp, &exp2);
    validate_tmp_file(&tmp_filename);
}

#[test]
fn native_both_stores_set_write_index_persistently() {
    // IndexedMzMLFileLoader.cpp:48 and :65, `options_.setWriteIndex(true)`,
    // change the loader's own options, not a copy.
    let mut file = IndexedMzMLFileLoader::new();
    file.options_mut().write_index = false;
    let dir = TempDir::new(false).unwrap();
    let exp2 = mzml_file_load(&indexed());
    file.store(dir.path().join("a.mzML"), &exp2).unwrap();
    assert!(file.options().write_index);

    file.options_mut().write_index = false;
    let mut exp = OnDiscPeakMap::new();
    assert!(file.load(indexed(), &mut exp).unwrap());
    file.store_on_disc(dir.path().join("b.mzML"), &mut exp)
        .unwrap();
    assert!(file.options().write_index);
    // Both files are indexed, whatever the options said before the call.
    assert!(mzml::has_index(dir.path().join("a.mzML")).unwrap());
    assert!(mzml::has_index(dir.path().join("b.mzML")).unwrap());
}

#[test]
fn native_on_disc_store_without_metadata_is_refused() {
    // The source dereferences the null settings pointer
    // (IndexedMzMLFileLoader.cpp:47); here the store is refused before a
    // byte is written.
    let mut file = IndexedMzMLFileLoader::new();
    let mut exp = OnDiscPeakMap::new();
    assert!(exp.open_file(indexed(), true).unwrap());
    let dir = TempDir::new(false).unwrap();
    let target = dir.path().join("refused.mzML");
    let error = file.store_on_disc(&target, &mut exp).unwrap_err();
    assert!(matches!(error, Error::MissingInformation(_)), "{error:?}");
    assert!(!target.exists());

    // A never-opened experiment has no metadata either.
    let mut empty = OnDiscPeakMap::new();
    let error = file.store_on_disc(&target, &mut empty).unwrap_err();
    assert!(matches!(error, Error::MissingInformation(_)), "{error:?}");
    assert!(!target.exists());
}

fn numpress(compression: NumpressCompression) -> NumpressConfig {
    NumpressConfig {
        compression,
        ..NumpressConfig::default()
    }
}

#[test]
fn native_both_stores_write_the_same_bytes_under_every_encoding() {
    // The source's two stores reach the same `MzMLHandler` writer with the
    // same `options_` (`IndexedMzMLFileLoader.cpp:49`, `:66`), one record at
    // a time or the whole map at once. Here the streamed on-disc store and the
    // whole-document in-memory store of the same records must agree byte for
    // byte, under every writer flag the loader forwards.
    let variants: [(&str, fn(&mut PeakFileOptions)); 9] = [
        ("defaults", |_| {}),
        ("zlib", |o| o.zlib_compression = true),
        ("32-bit m/z", |o| o.mz_32_bit = true),
        ("64-bit intensity", |o| o.intensity_32_bit = false),
        ("TPP compatibility", |o| o.force_tpp_compatibility = true),
        ("Numpress linear m/z", |o| {
            let _ = o.set_numpress_configuration_mass_time(numpress(NumpressCompression::Linear));
        }),
        ("Numpress slof intensity", |o| {
            o.set_numpress_configuration_intensity(numpress(NumpressCompression::Slof))
        }),
        ("Numpress pic intensity with zlib", |o| {
            o.set_numpress_configuration_intensity(numpress(NumpressCompression::Pic));
            o.zlib_compression = true;
        }),
        ("Numpress on every array kind", |o| {
            let _ = o.set_numpress_configuration_mass_time(numpress(NumpressCompression::Linear));
            o.set_numpress_configuration_intensity(numpress(NumpressCompression::Slof));
            o.set_numpress_configuration_float_data_array(numpress(NumpressCompression::Slof));
        }),
    ];
    let dir = TempDir::new(false).unwrap();
    for (label, change) in variants {
        let mut file = IndexedMzMLFileLoader::new();
        change(file.options_mut());
        let mut exp = OnDiscPeakMap::new();
        assert!(file.load(indexed(), &mut exp).unwrap());
        let streamed = dir.path().join("streamed.mzML");
        file.store_on_disc(&streamed, &mut exp).unwrap();
        let whole = dir.path().join("whole.mzML");
        let records = exp.load_experiment().unwrap();
        file.store(&whole, &records).unwrap();
        assert_eq!(
            std::fs::read(&streamed).unwrap(),
            std::fs::read(&whole).unwrap(),
            "{label}"
        );
        // And the file is an indexed mzML the loader reopens.
        let mut reopened = OnDiscPeakMap::new();
        assert!(file.load(&streamed, &mut reopened).unwrap(), "{label}");
        assert_eq!(reopened.spectrum_count(), 2, "{label}");
        assert_eq!(reopened.chromatogram_count(), 1, "{label}");
    }
}

#[test]
fn native_lossless_encodings_round_trip_through_the_on_disc_store() {
    // Precision the file already has and zlib change no decoded value.
    let variants: [fn(&mut PeakFileOptions); 2] = [
        |o| o.zlib_compression = true,
        |o| o.intensity_32_bit = false,
    ];
    let dir = TempDir::new(false).unwrap();
    for change in variants {
        let mut file = IndexedMzMLFileLoader::new();
        change(file.options_mut());
        let mut exp_ = OnDiscPeakMap::new();
        assert!(file.load(indexed(), &mut exp_).unwrap());
        let target = dir.path().join("lossless.mzML");
        file.store_on_disc(&target, &mut exp_).unwrap();
        let mut exp = OnDiscPeakMap::new();
        assert!(file.load(&target, &mut exp).unwrap());
        assert_matches_whole_file_load(&mut exp, &mzml_file_load(&indexed()));
    }
}

#[test]
fn native_in_memory_store_of_no_records_is_refused() {
    // An indexed mzML needs a record to index; the source writes a dummy `-1`
    // offset instead (CPP-050).
    let mut file = IndexedMzMLFileLoader::new();
    let dir = TempDir::new(false).unwrap();
    let target = dir.path().join("empty.mzML");
    let error = file.store(&target, &MSExperiment::default()).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    assert!(!target.exists());
}
