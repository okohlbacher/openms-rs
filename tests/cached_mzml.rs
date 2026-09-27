// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `src/format/cached_mzml_handler.rs` and `src/format/cached_mzml.rs`,
//! the ports of `FORMAT/HANDLERS/CachedMzMLHandler.h` and `FORMAT/CachedMzML.h`.
//!
//! Every `START_SECTION` of the upstream `CachedMzMLHandler_test.cpp` and
//! `CachedMzML_test.cpp` (core `bc9cc12`) is reproduced. Both use the
//! class-test file `MzMLFile_1.mzML`, which this repository already holds
//! byte-identical as `tests/data/mzml_validator/MzMLFile_1.mzML` (sha256
//! `076fd42e...0bc2`). The mzML metadata half needs the `mzml` feature.
#![cfg(feature = "mzml")]

use openms::Error;
use openms::format::cached_mzml::CachedmzML;
use openms::format::cached_mzml_handler::{CACHED_MZML_FILE_IDENTIFIER, CachedMzMLHandler};
use openms::kernel::{DataArray, MSChromatogram, MSExperiment, MSSpectrum};
use openms::system::file::TempDir;
use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};
use std::path::{Path, PathBuf};

fn mzml_file_1() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/mzml_validator/MzMLFile_1.mzML")
}

fn load_exp() -> MSExperiment {
    let exp = openms::format::mzml::load(mzml_file_1()).unwrap();
    assert!(!exp.spectra.is_empty());
    assert!(!exp.chromatograms.is_empty());
    exp
}

/// The class test's `cacheFile`: dump `MzMLFile_1.mzML` into a temporary file.
fn cache_file(dir: &TempDir, cache: &CachedMzMLHandler) -> (PathBuf, MSExperiment) {
    let exp = load_exp();
    let path = dir.path().join("cache.cached");
    cache.write_memdump(&exp, &path).unwrap();
    (path, exp)
}

/// The source's `TEST_REAL_SIMILAR` with its default tolerances.
fn similar(a: f64, b: f64) -> bool {
    if (a - b).abs() <= 1e-5 {
        return true;
    }
    if a == 0.0 || b == 0.0 {
        return false;
    }
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    ratio <= 1.0 + 1e-5 && ratio >= 1.0 / (1.0 + 1e-5)
}

fn message(error: Error) -> String {
    match error {
        Error::Parse { message, .. } => message,
        other => panic!("expected a parse error, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// CachedMzMLHandler_test.cpp
// ---------------------------------------------------------------------------

// START_SECTION(CachedMzMLHandler()) / START_SECTION(~CachedMzMLHandler())
#[test]
fn handler_construction() {
    let handler = CachedMzMLHandler::new();
    assert!(handler.spectra_index().is_empty());
    assert!(handler.chromatogram_index().is_empty());
}

// START_SECTION(( [EXTRA] testCaching))
#[test]
fn handler_test_caching() {
    let dir = TempDir::new(false).unwrap();
    let mut cache = CachedMzMLHandler::new();
    let (tmp_filename, exp) = cache_file(&dir, &cache);

    // Check whether spectra were written to disk correctly...
    cache.create_memdump_index(&tmp_filename).unwrap();
    let spectra_index = cache.spectra_index().to_vec();
    assert_eq!(spectra_index.len(), 4);
    let mut ifs = BufReader::new(File::open(&tmp_filename).unwrap());
    for i in 0..4 {
        ifs.seek(SeekFrom::Start(spectra_index[i])).unwrap();
        let (darray, _ms_level, _rt) = CachedMzMLHandler::read_spectrum_fast(&mut ifs).unwrap();
        assert!(darray.len() >= 2);
        assert_eq!(darray[0].data.len(), exp.spectra[i].peaks.len());
        assert_eq!(darray[1].data.len(), exp.spectra[i].peaks.len());
    }
    // test spec 1
    ifs.seek(SeekFrom::Start(spectra_index[1])).unwrap();
    let (darray, _, _) = CachedMzMLHandler::read_spectrum_fast(&mut ifs).unwrap();
    assert_eq!(darray.len(), 4);
    assert_eq!(darray[0].description, ""); // mz
    assert_eq!(darray[1].description, ""); // intensity
    assert_eq!(darray[2].description, "signal to noise array");
    assert_eq!(darray[3].description, "user-defined name");

    // Check whether chromatograms were written to disk correctly...
    cache.create_memdump_index(&tmp_filename).unwrap();
    let chrom_index = cache.chromatogram_index().to_vec();
    assert_eq!(chrom_index.len(), 2);
    for i in 0..2 {
        ifs.seek(SeekFrom::Start(chrom_index[i])).unwrap();
        let darray = CachedMzMLHandler::read_chromatogram_fast(&mut ifs).unwrap();
        assert!(darray.len() >= 2);
        assert_eq!(darray[0].data.len(), exp.chromatograms[i].peaks.len());
        assert_eq!(darray[1].data.len(), exp.chromatograms[i].peaks.len());
    }
}

// START_SECTION(( void writeMemdump(MapType& exp, std::string out )))
// NOT_TESTABLE upstream (see testCaching); the file must at least start with
// the identifier and end with the two counts.
#[test]
fn handler_write_memdump() {
    let dir = TempDir::new(false).unwrap();
    let cache = CachedMzMLHandler::new();
    let (path, exp) = cache_file(&dir, &cache);
    let bytes = std::fs::read(path).unwrap();
    assert_eq!(&bytes[..4], &CACHED_MZML_FILE_IDENTIFIER.to_le_bytes());
    let n = bytes.len();
    assert_eq!(
        &bytes[n - 16..n - 8],
        &(exp.spectra.len() as u64).to_le_bytes()
    );
    assert_eq!(
        &bytes[n - 8..],
        &(exp.chromatograms.len() as u64).to_le_bytes()
    );
}

fn cleared_spectrum(spectrum: &MSSpectrum) -> MSSpectrum {
    let mut spectrum = spectrum.clone();
    spectrum.peaks.clear();
    spectrum.float_data_arrays.clear();
    spectrum.integer_data_arrays.clear();
    spectrum.string_data_arrays.clear();
    spectrum
}

fn cleared_chromatogram(chromatogram: &MSChromatogram) -> MSChromatogram {
    let mut chromatogram = chromatogram.clone();
    chromatogram.peaks.clear();
    chromatogram.float_data_arrays.clear();
    chromatogram.integer_data_arrays.clear();
    chromatogram.string_data_arrays.clear();
    chromatogram
}

// START_SECTION(( void writeMetadata(MapType exp, std::string out_meta, bool addCacheMetaValue=false) ))
#[test]
fn handler_write_metadata() {
    let dir = TempDir::new(false).unwrap();
    let tmp_filename = dir.path().join("meta.mzML");
    let exp = load_exp();
    let cache = CachedMzMLHandler::new();

    // without adding the cache value, the meta data of the two experiments should be equal
    cache.write_metadata(&exp, &tmp_filename, false).unwrap();
    let meta_exp = openms::format::mzml::load(&tmp_filename).unwrap();
    assert_eq!(meta_exp.settings, exp.settings);
    assert_eq!(meta_exp.spectra[0], cleared_spectrum(&exp.spectra[0]));
    assert_eq!(
        meta_exp.chromatograms[0],
        cleared_chromatogram(&exp.chromatograms[0])
    );

    // with the cache value, the experimental settings are still equal
    cache.write_metadata(&exp, &tmp_filename, true).unwrap();
    let meta_exp = openms::format::mzml::load(&tmp_filename).unwrap();
    assert_eq!(meta_exp.settings, exp.settings);
    let processing = meta_exp.spectra[0].data_processing.last().unwrap();
    assert_eq!(processing.metadata["cached_data"].to_string(), "true");
}

// START_SECTION(( void readMemdump(MapType& exp_reading, std::string filename) const ))
#[test]
fn handler_read_memdump() {
    let dir = TempDir::new(false).unwrap();
    let cache = CachedMzMLHandler::new();
    let (tmp_filename, exp) = cache_file(&dir, &cache);
    let exp_new = cache.read_memdump(&tmp_filename).unwrap();
    assert_eq!(exp_new.spectra.len(), exp.spectra.len());
    assert_eq!(exp_new.chromatograms.len(), exp.chromatograms.len());

    match cache.read_memdump(dir.path().join("unused.cached")) {
        Err(Error::Io(error)) => assert_eq!(error.kind(), std::io::ErrorKind::NotFound),
        other => panic!("{other:?}"),
    }
    let error = cache.read_memdump(mzml_file_1()).unwrap_err();
    assert_eq!(
        message(error),
        "File might not be a cached mzML file (wrong file magic number). Aborting!"
    );
}

// START_SECTION(( void createMemdumpIndex(std::string filename) ))
#[test]
fn handler_create_memdump_index() {
    let dir = TempDir::new(false).unwrap();
    let mut cache = CachedMzMLHandler::new();
    let (tmp_filename, _) = cache_file(&dir, &cache);
    cache.create_memdump_index(&tmp_filename).unwrap();
    assert_eq!(cache.spectra_index().len(), 4);
    assert_eq!(cache.chromatogram_index().len(), 2);

    match cache.create_memdump_index(dir.path().join("unused.cached")) {
        Err(Error::Io(error)) => assert_eq!(error.kind(), std::io::ErrorKind::NotFound),
        other => panic!("{other:?}"),
    }
    assert!(matches!(
        cache.create_memdump_index(mzml_file_1()),
        Err(Error::Parse { .. })
    ));
    // A failed call keeps the previous indices.
    assert_eq!(cache.spectra_index().len(), 4);
}

// START_SECTION(( const std::vector<std::streampos>& getSpectraIndex() const ))
// START_SECTION(( const std::vector<std::streampos>& getChromatogramIndex() const ))
#[test]
fn handler_index_accessors() {
    let dir = TempDir::new(false).unwrap();
    let mut cache = CachedMzMLHandler::new();
    let (tmp_filename, _) = cache_file(&dir, &cache);
    cache.create_memdump_index(&tmp_filename).unwrap();
    assert_eq!(cache.spectra_index().len(), 4);
    assert_eq!(cache.chromatogram_index().len(), 2);
}

// START_SECTION(static inline void readSpectrumFast(...))
#[test]
fn handler_read_spectrum_fast() {
    let dir = TempDir::new(false).unwrap();
    let mut cache = CachedMzMLHandler::new();
    let (tmp_filename, exp) = cache_file(&dir, &cache);
    cache.create_memdump_index(&tmp_filename).unwrap();
    let spectra_index = cache.spectra_index().to_vec();
    assert_eq!(spectra_index.len(), 4);
    let mut ifs = BufReader::new(File::open(&tmp_filename).unwrap());

    ifs.seek(SeekFrom::Start(spectra_index[0])).unwrap();
    let (data, ms_level, rt) = CachedMzMLHandler::read_spectrum_fast(&mut ifs).unwrap();
    let (mz_array, intensity_array) = (&data[0], &data[1]);
    assert!(!mz_array.data.is_empty());
    assert_eq!(mz_array.data.len(), exp.spectra[0].peaks.len());
    assert_eq!(intensity_array.data.len(), exp.spectra[0].peaks.len());
    assert_eq!(ms_level, 1);
    assert!(similar(rt, 5.1));
    for (i, peak) in exp.spectra[0].peaks.iter().enumerate() {
        assert!(similar(mz_array.data[i], peak.mz));
        assert!(similar(intensity_array.data[i], f64::from(peak.intensity)));
    }

    // Check error conditions. The source's `seekg(-1)` puts its stream into
    // the fail state, so the header read fails; a Rust stream cannot seek
    // before its start, and the failed read is reproduced at the end instead.
    let end = ifs.seek(SeekFrom::End(0)).unwrap();
    let expected =
        "filestream in: Read an invalid spectrum length, something is wrong here. Aborting.";
    assert_eq!(
        message(CachedMzMLHandler::read_spectrum_fast(&mut ifs).unwrap_err()),
        expected
    );
    // should not read after the file ends
    ifs.seek(SeekFrom::Start(spectra_index.last().unwrap() * 20))
        .unwrap();
    assert!(spectra_index.last().unwrap() * 20 > end);
    assert_eq!(
        message(CachedMzMLHandler::read_spectrum_fast(&mut ifs).unwrap_err()),
        expected
    );
}

// START_SECTION(static inline void readChromatogramFast(...))
#[test]
fn handler_read_chromatogram_fast() {
    let dir = TempDir::new(false).unwrap();
    let mut cache = CachedMzMLHandler::new();
    let (tmp_filename, exp) = cache_file(&dir, &cache);
    cache.create_memdump_index(&tmp_filename).unwrap();
    let chrom_index = cache.chromatogram_index().to_vec();
    assert_eq!(chrom_index.len(), 2);
    let mut ifs = BufReader::new(File::open(&tmp_filename).unwrap());

    ifs.seek(SeekFrom::Start(chrom_index[0])).unwrap();
    let data = CachedMzMLHandler::read_chromatogram_fast(&mut ifs).unwrap();
    let (time_array, intensity_array) = (&data[0], &data[1]);
    assert!(!time_array.data.is_empty());
    assert_eq!(time_array.data.len(), exp.chromatograms[0].peaks.len());
    assert_eq!(intensity_array.data.len(), exp.chromatograms[0].peaks.len());
    for (i, peak) in exp.chromatograms[0].peaks.iter().enumerate() {
        assert!(similar(time_array.data[i], peak.rt));
        assert!(similar(intensity_array.data[i], f64::from(peak.intensity)));
    }

    let expected =
        "filestream in: Read an invalid chromatogram length, something is wrong here. Aborting.";
    ifs.seek(SeekFrom::End(0)).unwrap();
    assert_eq!(
        message(CachedMzMLHandler::read_chromatogram_fast(&mut ifs).unwrap_err()),
        expected
    );
    ifs.seek(SeekFrom::Start(chrom_index.last().unwrap() * 20))
        .unwrap();
    assert_eq!(
        message(CachedMzMLHandler::read_chromatogram_fast(&mut ifs).unwrap_err()),
        expected
    );
}

// ---------------------------------------------------------------------------
// CachedMzML_test.cpp
// ---------------------------------------------------------------------------

// START_SECTION(CachedmzML()) / START_SECTION(~CachedmzML())
#[test]
fn cached_construction() {
    let cache = CachedmzML::new();
    assert_eq!(cache.nr_spectra(), 0);
    assert_eq!(cache.nr_chromatograms(), 0);
}

fn stored_example(dir: &TempDir) -> (PathBuf, MSExperiment) {
    let exp = load_exp();
    let tmp = dir.path().join("example.mzML");
    CachedmzML::store(&tmp, &exp).unwrap();
    (tmp, exp)
}

// START_SECTION(( [EXTRA] testCaching))
#[test]
fn cached_test_caching() {
    let dir = TempDir::new(false).unwrap();
    let (tmp_filename, exp) = stored_example(&dir);

    let mut cache = CachedmzML::new();
    CachedmzML::load(&tmp_filename, &mut cache).unwrap();
    assert_eq!(cache.nr_spectra(), 4);
    for i in 0..4 {
        assert_eq!(
            cache.spectrum(i).unwrap().peaks.len(),
            exp.spectra[i].peaks.len()
        );
        // identical except DataProcessing (and extra data arrays)
        let mut tmp1 = cache.spectrum(i).unwrap();
        let mut tmp2 = exp.spectra[i].clone();
        tmp1.data_processing.clear();
        tmp2.data_processing.clear();
        tmp1.float_data_arrays.clear();
        tmp2.float_data_arrays.clear();
        assert_eq!(tmp1, tmp2, "spectrum {i}");
    }

    let scomp = &exp.spectra[1];
    assert_eq!(scomp.float_data_arrays.len(), 2);
    assert_eq!(scomp.integer_data_arrays.len(), 0);
    assert_eq!(scomp.string_data_arrays.len(), 0);
    let s = cache.spectrum(1).unwrap();
    assert_eq!(s.float_data_arrays.len(), 2);
    assert_eq!(s.integer_data_arrays.len(), 0);
    assert_eq!(s.string_data_arrays.len(), 0);
    assert_eq!(s.float_data_arrays[0].name, scomp.float_data_arrays[0].name);
    assert_eq!(s.float_data_arrays[1].name, scomp.float_data_arrays[1].name);
    assert_eq!(s.float_data_arrays[0].name, "signal to noise array");
    assert_eq!(s.float_data_arrays[1].name, "user-defined name");
    for k in 0..2 {
        let (a, b): (&DataArray<f32>, &DataArray<f32>) =
            (&s.float_data_arrays[k], &scomp.float_data_arrays[k]);
        assert_eq!(a.data.len(), b.data.len());
        for (x, y) in a.data.iter().zip(&b.data) {
            assert!(similar(f64::from(*x), f64::from(*y)));
        }
    }

    let mut cache = CachedmzML::new();
    CachedmzML::load(&tmp_filename, &mut cache).unwrap();
    assert_eq!(cache.nr_chromatograms(), 2);
    for i in 0..2 {
        let c = cache.chromatogram(i).unwrap();
        assert_eq!(c.peaks.len(), exp.chromatograms[i].peaks.len());
        assert_eq!(c.native_id, exp.chromatograms[i].native_id);
        assert_eq!(
            c.instrument_settings,
            exp.chromatograms[i].instrument_settings
        );
        // identical except DataProcessing
        let mut tmp1 = c;
        let mut tmp2 = exp.chromatograms[i].clone();
        tmp1.data_processing.clear();
        tmp2.data_processing.clear();
        assert_eq!(tmp1, tmp2, "chromatogram {i}");
    }
}

// START_SECTION(( size_t getNrSpectra() const ))
// START_SECTION(( size_t getNrChromatograms() const ))
// START_SECTION(( const MSExperiment& getMetaData() const )) (both sections)
#[test]
fn cached_counts_and_metadata() {
    let dir = TempDir::new(false).unwrap();
    let (tmp, exp) = stored_example(&dir);
    let cache_example = CachedmzML::open(&tmp).unwrap();
    assert_eq!(cache_example.nr_spectra(), 4);
    assert_eq!(cache_example.nr_chromatograms(), 2);
    assert_eq!(cache_example.meta_data().spectra.len(), 4);
    assert_eq!(cache_example.meta_data().chromatograms.len(), 2);

    for i in 0..4 {
        // identical except DataProcessing
        let mut tmp1 = cleared_spectrum(&cache_example.meta_data().spectra[i]);
        let mut tmp2 = cleared_spectrum(&exp.spectra[i]);
        tmp1.data_processing.clear();
        tmp2.data_processing.clear();
        assert_eq!(tmp1, tmp2, "spectrum settings {i}");
    }
    for i in 0..2 {
        let mut tmp1 = cleared_chromatogram(&cache_example.meta_data().chromatograms[i]);
        let mut tmp2 = cleared_chromatogram(&exp.chromatograms[i]);
        tmp1.data_processing.clear();
        tmp2.data_processing.clear();
        assert_eq!(tmp1, tmp2, "chromatogram settings {i}");
    }
}

// ---------------------------------------------------------------------------
// Native regressions derived from CachedMzMLHandler.cpp
// ---------------------------------------------------------------------------

#[test]
fn empty_records_keep_their_arrays_and_chromatogram_arrays_keep_values() {
    let mut exp = MSExperiment::new();
    let mut spectrum = MSSpectrum::default();
    spectrum
        .float_data_arrays
        .push(DataArray::new("empty", Vec::new()));
    spectrum
        .integer_data_arrays
        .push(DataArray::new("ints", Vec::new()));
    exp.spectra.push(spectrum);
    let mut second = MSSpectrum::default();
    second.ms_level = 2;
    second.rt = 7.5;
    second.peaks.push(openms::kernel::Peak1D::new(100.0, 5.0));
    second
        .integer_data_arrays
        .push(DataArray::new("charge", vec![3]));
    exp.spectra.push(second);
    let mut chromatogram = MSChromatogram::default();
    chromatogram
        .peaks
        .push(openms::kernel::ChromatogramPeak::new(1.0, 2.0));
    chromatogram
        .float_data_arrays
        .push(DataArray::new("fwhm", vec![0.25]));
    exp.chromatograms.push(chromatogram);

    let dir = TempDir::new(false).unwrap();
    let path = dir.path().join("x.cached");
    let handler = CachedMzMLHandler::new();
    handler.write_memdump(&exp, &path).unwrap();
    let back = handler.read_memdump(&path).unwrap();
    // The empty spectrum's two arrays are written (the source skips them) and
    // come back as float arrays.
    assert_eq!(back.spectra[0].float_data_arrays.len(), 2);
    assert_eq!(back.spectra[1].ms_level, 2);
    assert_eq!(back.spectra[1].rt, 7.5);
    assert_eq!(back.spectra[1].float_data_arrays[0].name, "charge");
    assert_eq!(back.spectra[1].float_data_arrays[0].data, vec![3.0]);
    // The source reads a chromatogram's arrays back empty; the values are kept here.
    assert_eq!(back.chromatograms[0].float_data_arrays[0].data, vec![0.25]);
}

#[test]
fn truncated_and_forged_files_are_refused() {
    let dir = TempDir::new(false).unwrap();
    let handler = CachedMzMLHandler::new();
    let (path, _) = cache_file(&dir, &handler);
    let bytes = std::fs::read(&path).unwrap();
    // Forge a huge spectrum count in the trailer.
    let mut forged = bytes.clone();
    let n = forged.len();
    forged[n - 16..n - 8].copy_from_slice(&u64::MAX.to_le_bytes());
    let forged_path = dir.path().join("forged.cached");
    std::fs::write(&forged_path, &forged).unwrap();
    assert!(matches!(
        handler.read_memdump(&forged_path),
        Err(Error::Parse { .. })
    ));
    let mut indexer = CachedMzMLHandler::new();
    assert!(matches!(
        indexer.create_memdump_index(&forged_path),
        Err(Error::Parse { .. })
    ));
    // Forge a huge peak count in the first spectrum.
    let mut forged = bytes.clone();
    forged[4..12].copy_from_slice(&(1u64 << 40).to_le_bytes());
    std::fs::write(&forged_path, &forged).unwrap();
    assert!(matches!(
        handler.read_memdump(&forged_path),
        Err(Error::Parse { .. })
    ));
    // Truncate.
    std::fs::write(&forged_path, &bytes[..10]).unwrap();
    assert!(matches!(
        handler.read_memdump(&forged_path),
        Err(Error::Parse { .. })
    ));
    // Out-of-range access through CachedmzML.
    let (tmp, _) = stored_example(&dir);
    let mut cache = CachedmzML::open(&tmp).unwrap();
    assert!(matches!(cache.spectrum(4), Err(Error::InvalidValue(_))));
    assert!(matches!(cache.chromatogram(2), Err(Error::InvalidValue(_))));
    let mut copy = cache.clone();
    assert_eq!(copy.spectrum(0).unwrap(), cache.spectrum(0).unwrap());
    assert!(matches!(
        CachedmzML::new().spectrum(0),
        Err(Error::InvalidValue(_))
    ));
}
