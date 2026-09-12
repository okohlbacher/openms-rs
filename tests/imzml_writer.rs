// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Coverage for `FORMAT/HANDLERS/ImzMLWriter.h`.
//!
//! The header has no class test of its own upstream. `ImzMLFile_test.cpp`
//! exercises the writer through `ImzMLFile::store`, and the literals below that
//! come from that suite are transcribed source review (tier 3): the round trip
//! of the two unmodified upstream fixtures in both modes, `MS:1003006` with its
//! single `MS:1002814` unit and `MS:1000786` for a free-text name,
//! `MS:1000821` "pressure array" with `UO:0000110` pascal, `MS:1003007` "raw
//! ion mobility array" written **without** a unit because it allows two, the
//! `binaryDataArrayList count` of 2 and 4, `MS:1000576` on both
//! `referenceableParamGroup`s, the refusal of a spectrum with no `imzml:x`, the
//! toleration of a duplicated pixel, the refusal of an incompatible continuous
//! mode, `FLOAT64` after `setMz32Bit(false)` / `setIntensity32Bit(false)`, the
//! m/z-range filter shrinking a spectrum, and the geometry and instrument
//! metadata round trip with the recomputed 300 µm extents.
//!
//! The strongest evidence here does not depend on that suite: every round trip
//! writes a dataset and reads it back with the stage-1
//! `openms::format::imzml_handler` reader, then compares the decoded m/z,
//! intensities, auxiliary arrays and pixel coordinates against what went in.
//! That closes the loop through the `.ibd` offsets, which is the part a
//! transcribed literal cannot check. The MD5 port is checked against the seven
//! RFC 1321 appendix A.5 test vectors, which are published values and not
//! OpenMS output.
//!
//! The resource ceilings, the atomicity of a rejected store, the deterministic
//! identifier derivation and the exact-versus-source shared-m/z tolerance are
//! independently derived (tier 4): the source has no analogue for any of them.

#![cfg(feature = "mzml")]

use openms::Error;
use openms::format::PeakFileOptions;
use openms::format::imzml_handler::{
    IBD_UUID_BYTES, ImagingMode, ImzMLDataType, ImzMLHandler, UuidStatus, infer_ibd_path,
};
use openms::format::imzml_writer::{
    FloatArraySkipReason, IBD_UUID_NAMESPACE, ImzMLWriteLimits, ImzMLWriteOptions,
    SOURCE_SHARED_MZ_TOLERANCE, StoreReport, apply_store_options, dataset_meta, derive_uuid,
    md5_hex, resolve_float_array_cv, spectra_share_mz, storage_mode, store, store_with_options,
};
use openms::kernel::{DataArray, MSExperiment, MSSpectrum, NumericRange, Peak1D, Precursor};
use openms::metadata::MetaValue;
use openms::system::file::TempDir;
use std::path::{Path, PathBuf};

const CONTINUOUS: &str = "ImzMLFile_1_Example_Continuous.imzML";
const PROCESSED: &str = "ImzMLFile_2_Example_Processed.imzML";
/// The upstream suite's grid for both fixtures.
const GRID: u32 = 3;
/// The upstream suite's spectrum count for both fixtures.
const SPECTRA: usize = 9;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// One peak list plus its pixel, the shape of the upstream suite's
/// `makePixelSpectrum_`.
fn pixel_spectrum(x: i64, y: i64, peaks: Vec<Peak1D>) -> MSSpectrum {
    let mut spectrum = MSSpectrum::from(peaks);
    spectrum
        .metadata
        .insert("imzml:x".into(), MetaValue::from(x));
    spectrum
        .metadata
        .insert("imzml:y".into(), MetaValue::from(y));
    spectrum
        .metadata
        .insert("imzml:z".into(), MetaValue::from(1_i64));
    spectrum
}

fn one_pixel(mz: f64, intensity: f32) -> MSExperiment {
    let mut experiment = MSExperiment::new();
    experiment
        .spectra
        .push(pixel_spectrum(1, 1, vec![Peak1D::new(mz, intensity)]));
    experiment
}

fn set_mode(experiment: &mut MSExperiment, mode: &str) {
    experiment
        .settings
        .metadata
        .insert("imzml:imaging_mode".into(), MetaValue::from(mode));
}

/// What a decoded dataset looks like, for comparing one side of a round trip
/// with the other.
#[derive(Debug, PartialEq)]
struct Decoded {
    coords: Vec<(u32, u32, u32)>,
    mz: Vec<Vec<f64>>,
    intensity: Vec<Vec<f32>>,
    aux: Vec<Vec<(String, Vec<f32>)>>,
}

fn decode(path: &Path) -> Decoded {
    let mut handler = ImzMLHandler::open(path).expect("the written imzML opens");
    let mut decoded = Decoded {
        coords: Vec::new(),
        mz: Vec::new(),
        intensity: Vec::new(),
        aux: Vec::new(),
    };
    for index in 0..handler.len() {
        let entry = handler.index()[index].clone();
        decoded.coords.push((entry.x, entry.y, entry.z));
        let spectrum = handler.spectrum(index).expect("pixel decodes").spectrum;
        decoded
            .mz
            .push(spectrum.peaks.iter().map(|peak| peak.mz).collect());
        decoded
            .intensity
            .push(spectrum.peaks.iter().map(|peak| peak.intensity).collect());
        decoded.aux.push(
            spectrum
                .float_data_arrays
                .iter()
                .map(|array| (array.name.clone(), array.data.clone()))
                .collect(),
        );
    }
    decoded
}

/// Build an in-memory experiment out of an upstream fixture through the stage-1
/// reader, mirroring the dataset metadata onto it the way `ImzMLFile::load`
/// does. This crate's whole-file mzML reader cannot read either fixture yet
/// (see `docs/IMZML_HANDLER_SUPPORT.md`), so this is the reference-data path.
fn load_fixture(name: &str) -> MSExperiment {
    let mut handler = ImzMLHandler::open(data(name)).expect("upstream fixture opens");
    let meta = handler.meta().clone();
    let mut experiment = MSExperiment::new();
    for index in 0..handler.len() {
        experiment
            .spectra
            .push(handler.spectrum(index).expect("pixel decodes").spectrum);
    }
    let settings = &mut experiment.settings.metadata;
    if let Some(mode) = meta.imaging_mode {
        settings.insert("imzml:imaging_mode".into(), MetaValue::from(mode.as_str()));
    }
    settings.insert("imzml:uuid".into(), MetaValue::from(meta.uuid.as_str()));
    settings.insert(
        "imzml:max_count_x".into(),
        MetaValue::from(meta.max_count_x),
    );
    settings.insert(
        "imzml:max_count_y".into(),
        MetaValue::from(meta.max_count_y),
    );
    settings.insert(
        "imzml:max_count_z".into(),
        MetaValue::from(meta.max_count_z),
    );
    for (key, value) in [
        ("imzml:pixel_size_x", meta.pixel_size_x),
        ("imzml:pixel_size_y", meta.pixel_size_y),
        ("imzml:max_dim_x", meta.max_dim_x),
        ("imzml:max_dim_y", meta.max_dim_y),
    ] {
        settings.insert(
            key.into(),
            MetaValue::try_from(value).expect("fixture geometry is finite"),
        );
    }
    for (key, value) in [
        ("imzml:scan_pattern", meta.scan_pattern.as_str()),
        ("imzml:scan_direction", meta.scan_direction.as_str()),
        (
            "imzml:line_scan_direction",
            meta.line_scan_direction.as_str(),
        ),
        ("imzml:polarity", meta.polarity.as_str()),
    ] {
        if !value.is_empty() {
            settings.insert(key.into(), MetaValue::from(value));
        }
    }
    experiment
}

/// Store into a fresh temporary directory and hand back the paths and report.
fn store_into(
    directory: &TempDir,
    name: &str,
    experiment: &MSExperiment,
    options: &PeakFileOptions,
) -> (PathBuf, StoreReport) {
    let path = directory.path().join(name);
    let report = store(&path, experiment, options).expect("store succeeds");
    (path, report)
}

fn xml(path: &Path) -> String {
    std::fs::read_to_string(path).expect("the written imzML is UTF-8")
}

// ---------------------------------------------------------------------------
// Round trips: the writer's primary evidence.
// ---------------------------------------------------------------------------

#[test]
fn continuous_fixture_round_trips_through_the_stage_one_reader() {
    let original = load_fixture(CONTINUOUS);
    assert_eq!(original.spectra.len(), SPECTRA);
    let before = decode(&data(CONTINUOUS));

    let temp = TempDir::new(false).unwrap();
    let (path, report) = store_into(
        &temp,
        "continuous.imzML",
        &original,
        &PeakFileOptions::new(),
    );

    assert_eq!(report.mode, ImagingMode::Continuous);
    assert_eq!(report.spectra_written, SPECTRA);
    assert_eq!(report.spectra_filtered_out, 0);
    assert_eq!(report.aux_arrays_written, 0);
    assert_eq!(report.duplicate_pixel_count, 0);
    assert_eq!(report.dropped_data_array_count, 0);
    assert_eq!(report.skipped_float_array_count, 0);
    assert_eq!(
        (report.meta.max_count_x, report.meta.max_count_y),
        (GRID, GRID)
    );
    assert_eq!(report.meta.max_count_z, 1);
    // Recomputed from the fixture's 100 µm pixels: the upstream suite's 300.
    assert!((report.meta.max_dim_x - 300.0).abs() < 1e-9);
    assert!((report.meta.max_dim_y - 300.0).abs() < 1e-9);
    // The declared identifier is kept verbatim, dashes and all.
    assert_eq!(report.meta.uuid, "12345678-1234-1234-1234-123456789012");

    // The shared axis is stored once, immediately after the UUID header, and
    // every pixel points at it.
    let mut handler = ImzMLHandler::open(&path).unwrap();
    assert_eq!(handler.len(), SPECTRA);
    assert_eq!(handler.meta().imaging_mode, Some(ImagingMode::Continuous));
    assert!(
        handler
            .index()
            .iter()
            .all(|entry| entry.mz_offset == IBD_UUID_BYTES as u64)
    );
    assert_eq!(handler.uuid_status().unwrap(), UuidStatus::Match);
    assert_eq!(decode(&path), before);

    // The .ibd holds the header, one shared float64 m/z array and nine float32
    // intensity arrays; the file is exactly that long.
    let peaks = before.mz[0].len() as u64;
    assert_eq!(
        report.ibd_bytes,
        IBD_UUID_BYTES as u64 + peaks * 8 + SPECTRA as u64 * peaks * 4
    );
    assert_eq!(
        std::fs::metadata(handler.ibd_path()).unwrap().len(),
        report.ibd_bytes
    );
}

#[test]
fn processed_fixture_round_trips_through_the_stage_one_reader() {
    let original = load_fixture(PROCESSED);
    assert_eq!(original.spectra.len(), SPECTRA);
    let before = decode(&data(PROCESSED));

    let temp = TempDir::new(false).unwrap();
    let (path, report) = store_into(&temp, "processed.imzML", &original, &PeakFileOptions::new());

    assert_eq!(report.mode, ImagingMode::Processed);
    assert_eq!(report.spectra_written, SPECTRA);
    let mut handler = ImzMLHandler::open(&path).unwrap();
    assert_eq!(handler.meta().imaging_mode, Some(ImagingMode::Processed));
    assert_eq!(handler.uuid_status().unwrap(), UuidStatus::Match);
    // Each pixel owns its axis, so no two m/z offsets agree.
    let mut offsets: Vec<u64> = handler.index().iter().map(|e| e.mz_offset).collect();
    offsets.sort_unstable();
    offsets.dedup();
    assert_eq!(offsets.len(), SPECTRA);
    assert_eq!(decode(&path), before);
    // The upstream suite's first m/z of processed pixel (1,1).
    assert!((before.mz[0][0] - 100.083336).abs() < 1e-5);
}

#[test]
fn a_written_dataset_declares_a_verifiable_sha1_and_a_well_formed_md5() {
    let temp = TempDir::new(false).unwrap();
    let (path, report) = store_into(
        &temp,
        "checksums.imzML",
        &one_pixel(100.0, 10.0),
        &PeakFileOptions::new(),
    );

    assert_eq!(report.meta.ibd_sha1.len(), 40);
    assert_eq!(report.meta.ibd_md5.len(), 32);
    assert!(report.meta.ibd_md5.bytes().all(|b| b.is_ascii_hexdigit()));

    let mut handler = ImzMLHandler::open(&path).unwrap();
    assert_eq!(handler.meta().ibd_sha1, report.meta.ibd_sha1);
    assert_eq!(handler.meta().ibd_md5, report.meta.ibd_md5);
    assert_eq!(
        handler.verify_ibd_sha1().unwrap(),
        openms::format::imzml_handler::ChecksumStatus::Match
    );
    let bytes = std::fs::read(handler.ibd_path()).unwrap();
    assert_eq!(md5_hex(&bytes), report.meta.ibd_md5);
}

// ---------------------------------------------------------------------------
// Auxiliary float arrays.
// ---------------------------------------------------------------------------

#[test]
fn ion_mobility_and_free_text_float_arrays_round_trip() {
    let mut experiment = MSExperiment::new();
    let mut spectrum = pixel_spectrum(
        1,
        1,
        vec![Peak1D::new(100.0, 10.0), Peak1D::new(200.0, 20.0)],
    );
    spectrum.float_data_arrays.push(DataArray::new(
        "mean inverse reduced ion mobility array",
        vec![0.85, 1.15],
    ));
    spectrum
        .float_data_arrays
        .push(DataArray::new("my custom SNR", vec![3.0, 4.5]));
    experiment.spectra.push(spectrum);
    set_mode(&mut experiment, "processed");

    let temp = TempDir::new(false).unwrap();
    let (path, report) = store_into(&temp, "aux.imzML", &experiment, &PeakFileOptions::new());
    assert_eq!(report.aux_arrays_written, 2);
    assert_eq!(report.skipped_float_array_count, 0);

    // The upstream suite inspects the XML for exactly these.
    let document = xml(&path);
    assert!(document.contains("MS:1003006"));
    assert!(document.contains("mean inverse reduced ion mobility array"));
    assert!(document.contains("MS:1002814"));
    assert!(document.contains("MS:1000786"));
    assert!(document.contains("my custom SNR"));
    assert!(document.contains("binaryDataArrayList count=\"4\""));

    let mut handler = ImzMLHandler::open(&path).unwrap();
    let entry = handler.index()[0].clone();
    assert_eq!(entry.aux.len(), 2);
    assert!(!entry.mz_compressed && !entry.int_compressed);
    let mobility = entry
        .aux
        .iter()
        .find(|array| array.name == "mean inverse reduced ion mobility array")
        .expect("the PSI ion-mobility array is indexed");
    assert_eq!(mobility.unit_accession, "MS:1002814");
    assert_eq!(mobility.accession, "MS:1003006");
    let custom = entry
        .aux
        .iter()
        .find(|array| array.name == "my custom SNR")
        .expect("the free-text array is indexed");
    assert_eq!(custom.accession, "MS:1000786");

    let decoded = handler.spectrum(0).unwrap().spectrum;
    assert_eq!(decoded.float_data_arrays.len(), 2);
    let by_name = |name: &str| -> Vec<f32> {
        decoded
            .float_data_arrays
            .iter()
            .find(|array| array.name == name)
            .expect("array survived the round trip")
            .data
            .clone()
    };
    assert_eq!(
        by_name("mean inverse reduced ion mobility array"),
        [0.85, 1.15]
    );
    assert_eq!(by_name("my custom SNR"), [3.0, 4.5]);
    assert_eq!(
        decoded.float_data_arrays[0]
            .metadata
            .get("unit_accession")
            .map(|value| value.to_string()),
        Some("MS:1002814".to_owned())
    );
}

#[test]
fn a_psi_term_with_one_allowed_unit_gets_it_and_a_term_with_two_gets_none() {
    let mut experiment = MSExperiment::new();
    let mut spectrum = pixel_spectrum(1, 1, vec![Peak1D::new(100.0, 10.0)]);
    spectrum
        .float_data_arrays
        .push(DataArray::new("pressure array", vec![101_325.0]));
    // Two allowed units, milliseconds and seconds: neither may be guessed.
    spectrum
        .float_data_arrays
        .push(DataArray::new("raw ion mobility array", vec![12.5]));
    experiment.spectra.push(spectrum);

    let temp = TempDir::new(false).unwrap();
    let (path, _) = store_into(&temp, "units.imzML", &experiment, &PeakFileOptions::new());
    let document = xml(&path);

    assert!(document.contains("MS:1000821"));
    assert!(document.contains("pressure array"));
    assert!(document.contains("unitAccession=\"UO:0000110\""));
    assert!(document.contains("unitName=\"pascal\""));
    assert!(document.contains("unitCvRef=\"UO\""));

    let raw = document
        .lines()
        .find(|line| line.contains("MS:1003007"))
        .expect("the raw ion mobility term is written");
    assert!(
        !raw.contains("unitAccession"),
        "a term with two allowed units must get no unit: {raw}"
    );
}

#[test]
fn an_unnamed_float_array_is_skipped_and_reported() {
    let mut experiment = MSExperiment::new();
    let mut spectrum = pixel_spectrum(1, 1, vec![Peak1D::new(100.0, 10.0)]);
    spectrum
        .float_data_arrays
        .push(DataArray::new("", vec![1.0]));
    experiment.spectra.push(spectrum);

    let temp = TempDir::new(false).unwrap();
    let (path, report) = store_into(&temp, "unnamed.imzML", &experiment, &PeakFileOptions::new());

    assert_eq!(report.aux_arrays_written, 0);
    assert_eq!(report.skipped_float_array_count, 1);
    assert_eq!(
        report.skipped_float_arrays[0].reason,
        FloatArraySkipReason::Unnamed
    );
    assert!(report.skipped_float_arrays[0].name.is_empty());

    let document = xml(&path);
    assert!(document.contains("binaryDataArrayList count=\"2\""));
    assert!(!document.contains("MS:1000786"));
    let handler = ImzMLHandler::open(&path).unwrap();
    assert!(handler.index()[0].aux.is_empty());
}

#[test]
fn float_arrays_named_after_the_peak_arrays_are_skipped() {
    let mut experiment = MSExperiment::new();
    let mut spectrum = pixel_spectrum(
        1,
        1,
        vec![Peak1D::new(100.0, 10.0), Peak1D::new(200.0, 20.0)],
    );
    spectrum
        .float_data_arrays
        .push(DataArray::new("m/z array", vec![999.0, 998.0]));
    spectrum
        .float_data_arrays
        .push(DataArray::new("intensity array", vec![1.0, 2.0]));
    experiment.spectra.push(spectrum);
    set_mode(&mut experiment, "processed");

    let temp = TempDir::new(false).unwrap();
    let (path, report) = store_into(&temp, "shadow.imzML", &experiment, &PeakFileOptions::new());

    assert_eq!(report.skipped_float_array_count, 2);
    let accessions: Vec<String> = report
        .skipped_float_arrays
        .iter()
        .map(|skipped| match &skipped.reason {
            FloatArraySkipReason::ReservedPeakArrayName { accession } => accession.clone(),
            other => panic!("unexpected reason {other:?}"),
        })
        .collect();
    assert_eq!(accessions, ["MS:1000514", "MS:1000515"]);

    assert!(xml(&path).contains("binaryDataArrayList count=\"2\""));
    let mut handler = ImzMLHandler::open(&path).unwrap();
    assert!(handler.index()[0].aux.is_empty());
    let decoded = handler.spectrum(0).unwrap().spectrum;
    assert_eq!(decoded.peaks.len(), 2);
    assert!((decoded.peaks[0].mz - 100.0).abs() < 1e-9);
    assert!((decoded.peaks[1].mz - 200.0).abs() < 1e-9);
    assert!((decoded.peaks[1].intensity - 20.0).abs() < 1e-6);
    assert!(decoded.float_data_arrays.is_empty());
}

#[test]
fn a_float_array_of_the_wrong_length_or_no_values_is_skipped() {
    let mut experiment = MSExperiment::new();
    let mut spectrum = pixel_spectrum(
        1,
        1,
        vec![Peak1D::new(100.0, 10.0), Peak1D::new(200.0, 20.0)],
    );
    spectrum
        .float_data_arrays
        .push(DataArray::new("temperature array", vec![1.0]));
    spectrum
        .float_data_arrays
        .push(DataArray::new("pressure array", Vec::new()));
    experiment.spectra.push(spectrum);
    set_mode(&mut experiment, "processed");

    let temp = TempDir::new(false).unwrap();
    let (path, report) = store_into(&temp, "lengths.imzML", &experiment, &PeakFileOptions::new());

    assert_eq!(report.aux_arrays_written, 0);
    assert_eq!(report.skipped_float_array_count, 2);
    assert_eq!(
        report.skipped_float_arrays[0].reason,
        FloatArraySkipReason::LengthMismatch {
            length: 1,
            peaks: 2
        }
    );
    assert_eq!(
        report.skipped_float_arrays[1].reason,
        FloatArraySkipReason::Empty
    );
    assert!(xml(&path).contains("binaryDataArrayList count=\"2\""));
}

/// One misaligned annotation array used to be skipped or fatal depending on an
/// unrelated option: with default options nothing in `apply_store_options`
/// touched the spectrum and the array was skipped, while any sort or trimming
/// filter reached `MSSpectrum::sort_by_position` or `MSSpectrum::select`, both
/// of which refuse a data array whose length is neither 0 nor the peak count.
/// The source diverges the same way — `applyStoreOptions_` reaches
/// `sortByPosition` and `select`, and both run `checkDataArraySizes_` — but the
/// writer's own explicit policy, in `appendAndWriteFloatDataArrays_`, is to
/// skip and warn. That is now what every option set does.
#[test]
fn a_misaligned_array_is_skipped_whatever_the_peak_file_options_say() {
    // Deliberately unsorted, so `sort_spectra_by_mz` really sorts. The
    // misaligned arrays are five long against three peaks, so no filter below
    // can leave them accidentally aligned.
    fn experiment() -> MSExperiment {
        let mut experiment = MSExperiment::new();
        let mut spectrum = pixel_spectrum(
            1,
            1,
            vec![
                Peak1D::new(300.0, 30.0),
                Peak1D::new(100.0, 10.0),
                Peak1D::new(200.0, 20.0),
            ],
        );
        spectrum.float_data_arrays.push(DataArray::new(
            "temperature array",
            vec![1.0, 2.0, 3.0, 4.0, 5.0],
        ));
        spectrum
            .float_data_arrays
            .push(DataArray::new("pressure array", vec![7.0, 8.0, 9.0]));
        spectrum
            .integer_data_arrays
            .push(DataArray::new("charge array", vec![1, 2, 3, 4, 5]));
        experiment.spectra.push(spectrum);
        set_mode(&mut experiment, "processed");
        experiment
    }

    // Neither sorts nor filters: nothing in `apply_store_options` touches the
    // spectrum, which is the one option set that always worked.
    let mut untouched = PeakFileOptions::new();
    untouched.sort_spectra_by_mz = false;
    // The library default. `sort_spectra_by_mz` is true, and the spectrum is
    // unsorted, so `MSSpectrum::sort_by_position` runs.
    let sorting = PeakFileOptions::new();
    // Trims a peak, so `MSSpectrum::select` runs.
    let mut trimming = PeakFileOptions::new();
    trimming.sort_spectra_by_mz = false;
    trimming.set_mz_range(NumericRange {
        min: 150.0,
        max: 350.0,
    });

    let temp = TempDir::new(false).unwrap();
    let mut reports = Vec::new();
    for (name, options) in [
        ("untouched.imzML", untouched),
        ("sorted.imzML", sorting),
        ("trimmed.imzML", trimming),
    ] {
        // The point of the test: none of the three is an error.
        let (path, report) = store_into(&temp, name, &experiment(), &options);
        let peaks = if name == "trimmed.imzML" { 2 } else { 3 };

        assert_eq!(report.skipped_float_array_count, 1, "{name}");
        assert_eq!(report.skipped_float_arrays[0].name, "temperature array");
        assert_eq!(
            report.skipped_float_arrays[0].reason,
            FloatArraySkipReason::LengthMismatch { length: 5, peaks },
            "{name}"
        );
        // The misaligned integer array is dropped like any other, not fatal.
        assert_eq!(report.dropped_data_array_count, 1, "{name}");
        assert_eq!(report.dropped_data_array_names, ["charge array"], "{name}");
        // The aligned array is still exported, and survives the reorder.
        assert_eq!(report.aux_arrays_written, 1, "{name}");
        let decoded = decode(&path);
        assert_eq!(decoded.aux[0].len(), 1, "{name}");
        assert_eq!(decoded.aux[0][0].0, "pressure array", "{name}");
        reports.push((report, decoded));
    }

    // The untouched and the sorting set change no peak count, so their skip
    // reports agree exactly — the same array, the same reason, the same numbers.
    assert_eq!(
        reports[0].0.skipped_float_arrays,
        reports[1].0.skipped_float_arrays
    );
    // Untouched: the peaks keep their input order, and so does the annotation.
    assert_eq!(reports[0].1.mz[0], vec![300.0, 100.0, 200.0]);
    assert_eq!(reports[0].1.aux[0][0].1, vec![7.0, 8.0, 9.0]);
    // Sorting reordered the peaks and carried the aligned array with them.
    assert_eq!(reports[1].1.mz[0], vec![100.0, 200.0, 300.0]);
    assert_eq!(reports[1].1.aux[0][0].1, vec![8.0, 9.0, 7.0]);
    // Trimming dropped the 100.0 peak and its annotation with it.
    assert_eq!(reports[2].1.mz[0], vec![300.0, 200.0]);
    assert_eq!(reports[2].1.aux[0][0].1, vec![7.0, 9.0]);
}

/// `apply_store_options` leaves a misaligned array exactly where it was, with
/// the values it had, so a caller that filters and then inspects sees the same
/// spectrum shape the writer reported on.
#[test]
fn apply_store_options_puts_a_misaligned_array_back_in_place() {
    let mut experiment = MSExperiment::new();
    let mut spectrum = MSSpectrum::from(vec![
        Peak1D::new(300.0, 30.0),
        Peak1D::new(100.0, 10.0),
        Peak1D::new(200.0, 20.0),
    ]);
    spectrum
        .float_data_arrays
        .push(DataArray::new("short", vec![1.0]));
    spectrum
        .float_data_arrays
        .push(DataArray::new("aligned", vec![1.0, 2.0, 3.0]));
    spectrum
        .float_data_arrays
        .push(DataArray::new("long", vec![1.0, 2.0, 3.0, 4.0]));
    experiment.spectra.push(spectrum);

    let mut options = PeakFileOptions::new();
    options.sort_spectra_by_mz = true;
    assert_eq!(apply_store_options(&mut experiment, &options).unwrap(), 0);

    let arrays = &experiment.spectra[0].float_data_arrays;
    let names: Vec<&str> = arrays.iter().map(|array| array.name.as_str()).collect();
    assert_eq!(names, ["short", "aligned", "long"]);
    assert_eq!(arrays[0].data, vec![1.0]);
    // Only the aligned array followed the sort.
    assert_eq!(arrays[1].data, vec![2.0, 3.0, 1.0]);
    assert_eq!(arrays[2].data, vec![1.0, 2.0, 3.0, 4.0]);
}

#[test]
fn integer_and_string_data_arrays_are_dropped_and_reported() {
    let mut experiment = MSExperiment::new();
    let mut spectrum = pixel_spectrum(1, 1, vec![Peak1D::new(100.0, 10.0)]);
    spectrum.float_data_arrays.push(DataArray::new(
        "mean inverse reduced ion mobility array",
        vec![0.85],
    ));
    spectrum
        .integer_data_arrays
        .push(DataArray::new("charge array", vec![2]));
    spectrum
        .string_data_arrays
        .push(DataArray::new("annotation array", vec!["peak0".to_owned()]));
    experiment.spectra.push(spectrum);
    set_mode(&mut experiment, "processed");

    let temp = TempDir::new(false).unwrap();
    let (path, report) = store_into(&temp, "dropped.imzML", &experiment, &PeakFileOptions::new());

    assert_eq!(report.dropped_data_array_count, 2);
    assert_eq!(
        report.dropped_data_array_names,
        ["charge array", "annotation array"]
    );
    assert_eq!(report.aux_arrays_written, 1);

    let mut handler = ImzMLHandler::open(&path).unwrap();
    let decoded = handler.spectrum(0).unwrap().spectrum;
    assert!(decoded.integer_data_arrays.is_empty());
    assert!(decoded.string_data_arrays.is_empty());
    assert_eq!(decoded.float_data_arrays.len(), 1);
    assert_eq!(decoded.float_data_arrays[0].data, [0.85]);
}

// ---------------------------------------------------------------------------
// Refusals and tolerations.
// ---------------------------------------------------------------------------

#[test]
fn a_spectrum_without_pixel_coordinates_is_refused() {
    let mut experiment = MSExperiment::new();
    experiment
        .spectra
        .push(MSSpectrum::from(vec![Peak1D::new(100.0, 1000.0)]));
    let temp = TempDir::new(false).unwrap();
    let path = temp.path().join("no-pixel.imzML");
    assert!(matches!(
        store(&path, &experiment, &PeakFileOptions::new()),
        Err(Error::MissingInformation(_))
    ));
    // Nothing was created: the refusal happens before either file is opened.
    assert!(!path.exists());
    assert!(!infer_ibd_path(&path).exists());
}

#[test]
fn a_pixel_coordinate_below_one_or_of_the_wrong_type_is_refused() {
    let temp = TempDir::new(false).unwrap();
    for (x, y) in [(0_i64, 1_i64), (1, 0), (-3, 1)] {
        let mut experiment = MSExperiment::new();
        experiment
            .spectra
            .push(pixel_spectrum(x, y, vec![Peak1D::new(100.0, 1.0)]));
        let path = temp.path().join(format!("pixel-{x}-{y}.imzML"));
        assert!(
            matches!(
                store(&path, &experiment, &PeakFileOptions::new()),
                Err(Error::InvalidValue(_))
            ),
            "({x},{y}) must be refused"
        );
        assert!(!path.exists());
    }

    let mut experiment = one_pixel(100.0, 1.0);
    experiment.spectra[0]
        .metadata
        .insert("imzml:y".into(), MetaValue::from("two"));
    let path = temp.path().join("pixel-type.imzML");
    assert!(matches!(
        store(&path, &experiment, &PeakFileOptions::new()),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn a_z_coordinate_below_one_is_refused_and_an_absent_one_defaults_to_one() {
    let temp = TempDir::new(false).unwrap();
    let mut experiment = one_pixel(100.0, 1.0);
    experiment.spectra[0]
        .metadata
        .insert("imzml:z".into(), MetaValue::from(0_i64));
    assert!(matches!(
        store(
            temp.path().join("z0.imzML"),
            &experiment,
            &PeakFileOptions::new()
        ),
        Err(Error::InvalidValue(_))
    ));

    let mut experiment = one_pixel(100.0, 1.0);
    experiment.spectra[0].metadata.remove("imzml:z");
    let (path, report) = store_into(&temp, "no-z.imzML", &experiment, &PeakFileOptions::new());
    assert_eq!(report.meta.max_count_z, 1);
    assert_eq!(ImzMLHandler::open(&path).unwrap().index()[0].z, 1);
}

#[test]
fn duplicate_pixel_coordinates_are_written_out_with_a_report() {
    let mut experiment = MSExperiment::new();
    experiment
        .spectra
        .push(pixel_spectrum(1, 1, vec![Peak1D::new(100.0, 1000.0)]));
    experiment
        .spectra
        .push(pixel_spectrum(1, 1, vec![Peak1D::new(101.0, 900.0)]));
    set_mode(&mut experiment, "processed");

    let temp = TempDir::new(false).unwrap();
    let (path, report) = store_into(
        &temp,
        "duplicate.imzML",
        &experiment,
        &PeakFileOptions::new(),
    );

    assert_eq!(report.spectra_written, 2);
    assert_eq!(report.duplicate_pixel_count, 1);
    assert_eq!(report.duplicate_pixels.len(), 1);
    assert_eq!(report.duplicate_pixels[0].spectrum, 1);
    assert_eq!(
        (
            report.duplicate_pixels[0].x,
            report.duplicate_pixels[0].y,
            report.duplicate_pixels[0].z
        ),
        (1, 1, 1)
    );

    // Both are readable by index; only the first is reachable by coordinate.
    let handler = ImzMLHandler::open(&path).unwrap();
    assert_eq!(handler.len(), 2);
    assert_eq!(handler.index_at_coord(1, 1, 1), Some(0));
}

#[test]
fn an_incompatible_continuous_mode_is_refused() {
    let mut experiment = MSExperiment::new();
    experiment
        .spectra
        .push(pixel_spectrum(1, 1, vec![Peak1D::new(100.0, 1000.0)]));
    experiment
        .spectra
        .push(pixel_spectrum(2, 1, vec![Peak1D::new(200.0, 800.0)]));
    set_mode(&mut experiment, "continuous");

    let temp = TempDir::new(false).unwrap();
    let path = temp.path().join("bad-continuous.imzML");
    assert!(matches!(
        store(&path, &experiment, &PeakFileOptions::new()),
        Err(Error::InvalidValue(_))
    ));
    assert!(!path.exists());
    assert!(!infer_ibd_path(&path).exists());
}

#[test]
fn an_experiment_with_no_spectra_is_refused() {
    let temp = TempDir::new(false).unwrap();
    assert!(matches!(
        store(
            temp.path().join("empty.imzML"),
            &MSExperiment::new(),
            &PeakFileOptions::new()
        ),
        Err(Error::MissingInformation(_))
    ));
}

#[test]
fn an_experiment_emptied_by_the_filters_is_refused() {
    let mut options = PeakFileOptions::new();
    options.set_ms_levels(&[7]).unwrap();
    let temp = TempDir::new(false).unwrap();
    assert!(matches!(
        store(
            temp.path().join("filtered-empty.imzML"),
            &one_pixel(100.0, 1.0),
            &options
        ),
        Err(Error::MissingInformation(_))
    ));
}

#[test]
fn a_native_identifier_xml_cannot_represent_is_refused() {
    let mut experiment = one_pixel(100.0, 1.0);
    experiment.spectra[0].native_id = "pixel\u{1}one".to_owned();
    let temp = TempDir::new(false).unwrap();
    let path = temp.path().join("bad-id.imzML");
    assert!(matches!(
        store(&path, &experiment, &PeakFileOptions::new()),
        Err(Error::InvalidValue(_))
    ));
    assert!(!path.exists());
}

#[test]
fn a_non_finite_shared_mz_tolerance_is_refused() {
    let temp = TempDir::new(false).unwrap();
    for tolerance in [f64::NAN, f64::INFINITY, -1.0] {
        let options = ImzMLWriteOptions {
            shared_mz_tolerance: tolerance,
            ..ImzMLWriteOptions::default()
        };
        assert!(matches!(
            store_with_options(
                temp.path().join("tolerance.imzML"),
                &one_pixel(100.0, 1.0),
                &PeakFileOptions::new(),
                &options,
                &mut openms::concept::progress_logger::ProgressLogger::new(),
            ),
            Err(Error::InvalidValue(_))
        ));
    }
}

// ---------------------------------------------------------------------------
// Storage-mode choice and the shared-m/z tolerance.
// ---------------------------------------------------------------------------

fn near_identical_axes() -> MSExperiment {
    let mut experiment = MSExperiment::new();
    experiment
        .spectra
        .push(pixel_spectrum(1, 1, vec![Peak1D::new(100.0, 1.0)]));
    experiment
        .spectra
        .push(pixel_spectrum(2, 1, vec![Peak1D::new(100.000_001, 2.0)]));
    experiment
}

#[test]
fn the_native_default_refuses_to_discard_an_almost_shared_axis() {
    let experiment = near_identical_axes();
    assert!(!spectra_share_mz(&experiment, 0.0));
    assert!(spectra_share_mz(&experiment, SOURCE_SHARED_MZ_TOLERANCE));

    let temp = TempDir::new(false).unwrap();
    // Auto-detection: exact comparison keeps both axes, so processed wins.
    let (path, report) = store_into(&temp, "near.imzML", &experiment, &PeakFileOptions::new());
    assert_eq!(report.mode, ImagingMode::Processed);
    let mut handler = ImzMLHandler::open(&path).unwrap();
    assert!((handler.spectrum(1).unwrap().spectrum.peaks[0].mz - 100.000_001).abs() < 1e-12);

    // The source's tolerance collapses the two axes into one.
    let mut logger = openms::concept::progress_logger::ProgressLogger::new();
    let source_path = temp.path().join("near-source.imzML");
    let report = store_with_options(
        &source_path,
        &experiment,
        &PeakFileOptions::new(),
        &ImzMLWriteOptions::source(),
        &mut logger,
    )
    .unwrap();
    assert_eq!(report.mode, ImagingMode::Continuous);
    let mut handler = ImzMLHandler::open(&source_path).unwrap();
    // The second pixel now reads back the first pixel's m/z: that is the loss.
    assert!((handler.spectrum(1).unwrap().spectrum.peaks[0].mz - 100.0).abs() < 1e-12);
}

#[test]
fn an_explicit_continuous_mode_over_an_almost_shared_axis_needs_the_source_tolerance() {
    let mut experiment = near_identical_axes();
    set_mode(&mut experiment, "continuous");
    let temp = TempDir::new(false).unwrap();
    assert!(matches!(
        store(
            temp.path().join("explicit.imzML"),
            &experiment,
            &PeakFileOptions::new()
        ),
        Err(Error::InvalidValue(_))
    ));
    let mut logger = openms::concept::progress_logger::ProgressLogger::new();
    assert!(
        store_with_options(
            temp.path().join("explicit-source.imzML"),
            &experiment,
            &PeakFileOptions::new(),
            &ImzMLWriteOptions::source(),
            &mut logger,
        )
        .is_ok()
    );
}

#[test]
fn storage_mode_honours_a_declaration_and_otherwise_detects() {
    let shared = {
        let mut experiment = MSExperiment::new();
        experiment
            .spectra
            .push(pixel_spectrum(1, 1, vec![Peak1D::new(100.0, 1.0)]));
        experiment
            .spectra
            .push(pixel_spectrum(2, 1, vec![Peak1D::new(100.0, 2.0)]));
        experiment
    };
    let mut meta = dataset_meta(&shared).unwrap();
    assert_eq!(meta.imaging_mode, None);
    assert_eq!(
        storage_mode(&shared, &meta, 0.0).unwrap(),
        ImagingMode::Continuous
    );
    // An explicit "processed" is honoured unconditionally.
    meta.imaging_mode = Some(ImagingMode::Processed);
    assert_eq!(
        storage_mode(&shared, &meta, 0.0).unwrap(),
        ImagingMode::Processed
    );
    meta.imaging_mode = Some(ImagingMode::Continuous);
    assert_eq!(
        storage_mode(&shared, &meta, 0.0).unwrap(),
        ImagingMode::Continuous
    );
}

#[test]
fn spectra_share_mz_needs_a_non_empty_reference_and_equal_lengths() {
    assert!(!spectra_share_mz(&MSExperiment::new(), 0.0));

    // One non-empty spectrum trivially shares its own axis.
    assert!(spectra_share_mz(&one_pixel(100.0, 1.0), 0.0));

    // An all-empty experiment has no reference axis at all.
    let mut empty = MSExperiment::new();
    empty.spectra.push(pixel_spectrum(1, 1, Vec::new()));
    assert!(!spectra_share_mz(&empty, 0.0));

    // A mix of empty and non-empty does not share.
    let mut mixed = MSExperiment::new();
    mixed
        .spectra
        .push(pixel_spectrum(1, 1, vec![Peak1D::new(100.0, 1.0)]));
    mixed.spectra.push(pixel_spectrum(2, 1, Vec::new()));
    assert!(!spectra_share_mz(&mixed, 0.0));

    // Unequal lengths do not share.
    let mut lengths = MSExperiment::new();
    lengths
        .spectra
        .push(pixel_spectrum(1, 1, vec![Peak1D::new(100.0, 1.0)]));
    lengths.spectra.push(pixel_spectrum(
        2,
        1,
        vec![Peak1D::new(100.0, 1.0), Peak1D::new(200.0, 2.0)],
    ));
    assert!(!spectra_share_mz(&lengths, 0.0));
}

#[test]
fn an_all_empty_experiment_stores_as_processed_with_zero_length_arrays() {
    let mut experiment = MSExperiment::new();
    experiment.spectra.push(pixel_spectrum(1, 1, Vec::new()));
    experiment.spectra.push(pixel_spectrum(2, 1, Vec::new()));

    let temp = TempDir::new(false).unwrap();
    let (path, report) = store_into(
        &temp,
        "all-empty.imzML",
        &experiment,
        &PeakFileOptions::new(),
    );
    assert_eq!(report.mode, ImagingMode::Processed);
    assert_eq!(report.ibd_bytes, IBD_UUID_BYTES as u64);

    let mut handler = ImzMLHandler::open(&path).unwrap();
    assert_eq!(handler.len(), 2);
    assert_eq!(handler.index()[0].mz_length, 0);
    assert!(handler.spectrum(1).unwrap().spectrum.peaks.is_empty());
}

// ---------------------------------------------------------------------------
// PeakFileOptions.
// ---------------------------------------------------------------------------

#[test]
fn binary_precision_follows_the_peak_file_options() {
    let temp = TempDir::new(false).unwrap();
    let experiment = one_pixel(100.5, 10.25);

    let mut wide = PeakFileOptions::new();
    wide.mz_32_bit = false;
    wide.intensity_32_bit = false;
    let (path, report) = store_into(&temp, "float64.imzML", &experiment, &wide);
    assert_eq!(report.meta.mz_data_type, ImzMLDataType::Float64);
    assert_eq!(report.meta.int_data_type, ImzMLDataType::Float64);
    let handler = ImzMLHandler::open(&path).unwrap();
    assert_eq!(handler.index()[0].mz_type, ImzMLDataType::Float64);
    assert_eq!(handler.index()[0].int_type, ImzMLDataType::Float64);
    assert_eq!(handler.index()[0].mz_encoded_bytes, 8);
    assert_eq!(handler.index()[0].int_encoded_bytes, 8);

    let mut narrow = PeakFileOptions::new();
    narrow.mz_32_bit = true;
    narrow.intensity_32_bit = true;
    let (path, report) = store_into(&temp, "float32.imzML", &experiment, &narrow);
    assert_eq!(report.meta.mz_data_type, ImzMLDataType::Float32);
    assert_eq!(report.ibd_bytes, IBD_UUID_BYTES as u64 + 4 + 4);
    let handler = ImzMLHandler::open(&path).unwrap();
    assert_eq!(handler.index()[0].mz_type, ImzMLDataType::Float32);
    assert_eq!(handler.index()[0].int_type, ImzMLDataType::Float32);
}

#[test]
fn the_mz_range_filter_shrinks_every_spectrum_before_the_write() {
    let original = load_fixture(CONTINUOUS);
    let full = original.spectra[0].peaks.len();
    let low = original.spectra[0].peaks[0].mz;
    let high = original.spectra[0].peaks[full - 1].mz;
    let middle = (low + high) / 2.0;

    let mut options = PeakFileOptions::new();
    options.set_mz_range(NumericRange {
        min: low,
        max: middle,
    });

    let temp = TempDir::new(false).unwrap();
    let (path, report) = store_into(&temp, "range.imzML", &original, &options);
    assert_eq!(report.spectra_written, SPECTRA);
    assert_eq!(report.spectra_filtered_out, 0);

    let mut handler = ImzMLHandler::open(&path).unwrap();
    let peaks = handler.spectrum(0).unwrap().spectrum.peaks;
    assert!(peaks.len() < full);
    assert!(!peaks.is_empty());
    assert!(peaks.first().unwrap().mz >= low);
    // DRange::encloses is open at the maximum, so every kept m/z is below it.
    assert!(peaks.last().unwrap().mz < middle);
}

#[test]
fn metadata_only_writes_the_geometry_and_no_peaks() {
    let mut options = PeakFileOptions::new();
    options.metadata_only = true;
    let temp = TempDir::new(false).unwrap();
    let (path, report) = store_into(
        &temp,
        "metadata-only.imzML",
        &load_fixture(PROCESSED),
        &options,
    );

    assert_eq!(report.spectra_written, SPECTRA);
    assert_eq!(report.mode, ImagingMode::Processed);
    assert_eq!(report.ibd_bytes, IBD_UUID_BYTES as u64);
    let mut handler = ImzMLHandler::open(&path).unwrap();
    assert_eq!(
        (handler.meta().max_count_x, handler.meta().max_count_y),
        (GRID, GRID)
    );
    assert!(handler.spectrum(0).unwrap().spectrum.peaks.is_empty());
}

/// A finding, not a port defect. `applyStoreOptions_` clears every peak for
/// `metadata_only`, and `isContinuousMode_` then runs `spectraShareMz_` over
/// the emptied spectra, which reports no shared axis because it cannot find a
/// non-empty reference. An experiment that declares `imzml:imaging_mode =
/// "continuous"` — which is what loading any continuous imzML produces —
/// therefore fails its own metadata-only store with
/// `Exception::InvalidParameter`, in the source exactly as here. Recorded in
/// `docs/IMZML_WRITER_SUPPORT.md`; the port reproduces it rather than quietly
/// downgrading the declared mode.
#[test]
fn metadata_only_over_a_declared_continuous_dataset_is_refused_as_in_the_source() {
    let mut options = PeakFileOptions::new();
    options.metadata_only = true;
    let experiment = load_fixture(CONTINUOUS);
    assert_eq!(
        dataset_meta(&experiment).unwrap().imaging_mode,
        Some(ImagingMode::Continuous)
    );
    let temp = TempDir::new(false).unwrap();
    let path = temp.path().join("metadata-only-continuous.imzML");
    assert!(matches!(
        store(&path, &experiment, &options),
        Err(Error::InvalidValue(_))
    ));
    assert!(!path.exists());

    // Dropping the declaration lets the same store succeed as processed.
    let mut relaxed = experiment.clone();
    relaxed.settings.metadata.remove("imzml:imaging_mode");
    let (_, report) = store_into(&temp, "metadata-only-auto.imzML", &relaxed, &options);
    assert_eq!(report.mode, ImagingMode::Processed);
}

#[test]
fn apply_store_options_filters_levels_times_precursors_and_peaks() {
    let mut experiment = MSExperiment::new();
    for (level, rt) in [(1_u32, 10.0), (2, 20.0), (1, 30.0)] {
        let mut spectrum = pixel_spectrum(1, 1, vec![Peak1D::new(100.0, 1.0)]);
        spectrum.ms_level = level;
        spectrum.rt = rt;
        spectrum.precursors.push(Precursor::new(500.0, 2));
        experiment.spectra.push(spectrum);
    }

    let mut levels = experiment.clone();
    let mut options = PeakFileOptions::new();
    options.set_ms_levels(&[1]).unwrap();
    assert_eq!(apply_store_options(&mut levels, &options).unwrap(), 1);
    assert_eq!(levels.spectra.len(), 2);

    // Closed at the minimum, open at the maximum: 10 is kept, 30 is not.
    let mut times = experiment.clone();
    let mut options = PeakFileOptions::new();
    options.set_rt_range(NumericRange {
        min: 10.0,
        max: 30.0,
    });
    assert_eq!(apply_store_options(&mut times, &options).unwrap(), 1);
    assert_eq!(times.spectra.len(), 2);
    assert!((times.spectra[0].rt - 10.0).abs() < 1e-9);

    // Only the MS2 spectrum is tested against the precursor range.
    let mut precursors = experiment.clone();
    let mut options = PeakFileOptions::new();
    options.set_precursor_mz_range(NumericRange {
        min: 600.0,
        max: 700.0,
    });
    assert_eq!(apply_store_options(&mut precursors, &options).unwrap(), 1);
    assert_eq!(precursors.spectra.len(), 2);
    assert!(precursors.spectra.iter().all(|s| s.ms_level == 1));

    // An intensity range drops peaks, not spectra.
    let mut intensities = experiment.clone();
    let mut options = PeakFileOptions::new();
    options.set_intensity_range(NumericRange {
        min: 5.0,
        max: 10.0,
    });
    assert_eq!(apply_store_options(&mut intensities, &options).unwrap(), 0);
    assert!(intensities.spectra.iter().all(|s| s.peaks.is_empty()));
}

#[test]
fn apply_store_options_sorts_by_mz_and_moves_the_annotations_with_it() {
    let mut experiment = MSExperiment::new();
    let mut spectrum = pixel_spectrum(1, 1, vec![Peak1D::new(200.0, 2.0), Peak1D::new(100.0, 1.0)]);
    spectrum
        .float_data_arrays
        .push(DataArray::new("temperature array", vec![20.0, 10.0]));
    experiment.spectra.push(spectrum);

    let mut options = PeakFileOptions::new();
    options.sort_spectra_by_mz = true;
    assert_eq!(apply_store_options(&mut experiment, &options).unwrap(), 0);
    assert!((experiment.spectra[0].peaks[0].mz - 100.0).abs() < 1e-9);
    assert_eq!(
        experiment.spectra[0].float_data_arrays[0].data,
        [10.0, 20.0]
    );

    // With sorting off the stored order is kept.
    let mut unsorted = MSExperiment::new();
    unsorted.spectra.push(pixel_spectrum(
        1,
        1,
        vec![Peak1D::new(200.0, 2.0), Peak1D::new(100.0, 1.0)],
    ));
    let mut options = PeakFileOptions::new();
    options.sort_spectra_by_mz = false;
    apply_store_options(&mut unsorted, &options).unwrap();
    assert!((unsorted.spectra[0].peaks[0].mz - 200.0).abs() < 1e-9);
}

/// This used to assert the opposite: that `apply_store_options` refused an
/// unaligned annotation array, because the sort it runs does. It no longer
/// does, because the writer's policy for such an array is to skip it and say
/// so, and that policy must not turn on whether an unrelated option happens to
/// reach the sort. The array is lifted off for the sort and put back after.
#[test]
fn apply_store_options_sorts_around_an_unaligned_annotation_array() {
    let mut experiment = MSExperiment::new();
    let mut spectrum = pixel_spectrum(1, 1, vec![Peak1D::new(200.0, 2.0), Peak1D::new(100.0, 1.0)]);
    spectrum
        .float_data_arrays
        .push(DataArray::new("temperature array", vec![20.0]));
    experiment.spectra.push(spectrum);

    let mut options = PeakFileOptions::new();
    options.sort_spectra_by_mz = true;
    assert_eq!(apply_store_options(&mut experiment, &options).unwrap(), 0);
    // The peaks are sorted, and the unaligned array is untouched and in place.
    assert!((experiment.spectra[0].peaks[0].mz - 100.0).abs() < 1e-9);
    let arrays = &experiment.spectra[0].float_data_arrays;
    assert_eq!(arrays.len(), 1);
    assert_eq!(arrays[0].name, "temperature array");
    assert_eq!(arrays[0].data, vec![20.0]);

    // A spectrum the sort genuinely refuses is still left exactly as it was.
    let mut nonfinite = MSExperiment::new();
    let mut spectrum = pixel_spectrum(1, 1, vec![Peak1D::new(200.0, 2.0), Peak1D::new(100.0, 1.0)]);
    spectrum.peaks[0].mz = f64::NAN;
    spectrum
        .float_data_arrays
        .push(DataArray::new("temperature array", vec![20.0]));
    nonfinite.spectra.push(spectrum);
    assert!(matches!(
        apply_store_options(&mut nonfinite, &options),
        Err(Error::InvalidValue(_))
    ));
    assert!(nonfinite.spectra[0].peaks[0].mz.is_nan());
    assert_eq!(nonfinite.spectra[0].float_data_arrays.len(), 1);
    assert_eq!(
        nonfinite.spectra[0].float_data_arrays[0].data,
        vec![20.0_f32]
    );
}

// ---------------------------------------------------------------------------
// Dataset metadata.
// ---------------------------------------------------------------------------

#[test]
fn the_acquisition_geometry_and_instrument_model_round_trip() {
    let mut experiment = load_fixture(CONTINUOUS);
    for (key, value) in [
        ("imzml:scan_pattern", "top down"),
        ("imzml:scan_direction", "flyback"),
        ("imzml:line_scan_direction", "left-right"),
        ("imzml:polarity", "positive"),
    ] {
        experiment
            .settings
            .metadata
            .insert(key.into(), MetaValue::from(value));
    }
    experiment.settings.instrument.model = "Test MSI Instrument".to_owned();
    experiment.spectra[0].rt = 12.34;
    experiment.spectra[0].ms_level = 1;

    let temp = TempDir::new(false).unwrap();
    let (path, report) = store_into(&temp, "meta.imzML", &experiment, &PeakFileOptions::new());

    assert_eq!(report.meta.scan_pattern, "top down");
    assert_eq!(report.meta.scan_direction, "flyback");
    assert_eq!(report.meta.line_scan_direction, "left-right");
    assert_eq!(report.meta.polarity, "positive");
    assert!(!report.meta.ibd_sha1.is_empty());
    assert!((report.meta.max_dim_x - 300.0).abs() < 1e-9);
    assert!((report.meta.max_dim_y - 300.0).abs() < 1e-9);

    let handler = ImzMLHandler::open(&path).unwrap();
    let meta = handler.meta();
    assert_eq!(meta.scan_pattern, "top down");
    assert_eq!(meta.scan_direction, "flyback");
    assert_eq!(meta.line_scan_direction, "left-right");
    assert_eq!(meta.polarity, "positive");
    assert!((meta.pixel_size_x - 100.0).abs() < 1e-9);
    assert!((meta.max_dim_x - 300.0).abs() < 1e-9);

    let document = xml(&path);
    assert!(document.contains("value=\"Test MSI Instrument\""));
    assert!(document.contains("IMS:1000401"));
    assert!(document.contains("IMS:1000413"));
    assert!(document.contains("IMS:1000491"));
    assert!(document.contains("MS:1000130"));
    assert!(document.contains("value=\"12.34\""));
}

#[test]
fn an_unrecognised_geometry_term_writes_nothing_and_the_model_falls_back() {
    let mut experiment = one_pixel(100.0, 1.0);
    for (key, value) in [
        ("imzml:scan_pattern", "sideways"),
        ("imzml:scan_direction", "diagonal"),
        ("imzml:line_scan_direction", "up-down"),
        ("imzml:polarity", "sideways"),
    ] {
        experiment
            .settings
            .metadata
            .insert(key.into(), MetaValue::from(value));
    }
    experiment.settings.instrument.name = "Named Only".to_owned();

    let temp = TempDir::new(false).unwrap();
    let (path, _) = store_into(
        &temp,
        "unknown-geom.imzML",
        &experiment,
        &PeakFileOptions::new(),
    );
    let document = xml(&path);
    for accession in [
        "IMS:1000401",
        "IMS:1000402",
        "IMS:1000412",
        "IMS:1000413",
        "IMS:1000480",
        "IMS:1000481",
        "IMS:1000491",
        "IMS:1000492",
        "MS:1000129",
        "MS:1000130",
    ] {
        assert!(!document.contains(accession), "{accession} must not appear");
    }
    // The instrument name stands in for a missing model.
    assert!(document.contains("value=\"Named Only\""));

    // With neither, the source's placeholder is written.
    let temp2 = TempDir::new(false).unwrap();
    let (path, _) = store_into(
        &temp2,
        "no-model.imzML",
        &one_pixel(100.0, 1.0),
        &PeakFileOptions::new(),
    );
    assert!(xml(&path).contains("value=\"OpenMS export\""));
}

#[test]
fn a_declared_grid_larger_than_the_pixels_survives() {
    let mut experiment = one_pixel(100.0, 1.0);
    let settings = &mut experiment.settings.metadata;
    settings.insert("imzml:max_count_x".into(), MetaValue::from(40_u32));
    settings.insert("imzml:max_count_y".into(), MetaValue::from(50_u32));
    settings.insert(
        "imzml:pixel_size_x".into(),
        MetaValue::try_from(2.0_f64).unwrap(),
    );
    // An explicit extent is overwritten by pixel size times pixel count.
    settings.insert(
        "imzml:max_dim_x".into(),
        MetaValue::try_from(9999.0_f64).unwrap(),
    );

    let temp = TempDir::new(false).unwrap();
    let (_, report) = store_into(&temp, "grid.imzML", &experiment, &PeakFileOptions::new());
    assert_eq!((report.meta.max_count_x, report.meta.max_count_y), (40, 50));
    assert!((report.meta.max_dim_x - 80.0).abs() < 1e-9);
    // No pixel size in y, so the declared extent of 0 stays 0 and is not written.
    assert!((report.meta.max_dim_y - 0.0).abs() < 1e-9);
}

#[test]
fn dataset_meta_reads_every_key_and_refuses_the_wrong_type() {
    let mut experiment = MSExperiment::new();
    let settings = &mut experiment.settings.metadata;
    settings.insert("imzml:imaging_mode".into(), MetaValue::from("continuous"));
    settings.insert("imzml:max_count_x".into(), MetaValue::from(4_u32));
    settings.insert("imzml:max_count_y".into(), MetaValue::from(5_u32));
    settings.insert("imzml:max_count_z".into(), MetaValue::from(6_u32));
    settings.insert(
        "imzml:pixel_size_x".into(),
        MetaValue::try_from(1.5_f64).unwrap(),
    );
    settings.insert(
        "imzml:pixel_size_y".into(),
        MetaValue::try_from(2.5_f64).unwrap(),
    );
    settings.insert(
        "imzml:max_dim_x".into(),
        MetaValue::try_from(6.0_f64).unwrap(),
    );
    settings.insert(
        "imzml:max_dim_y".into(),
        MetaValue::try_from(12.5_f64).unwrap(),
    );
    settings.insert("imzml:uuid".into(), MetaValue::from("{ABCD}"));
    settings.insert("imzml:scan_pattern".into(), MetaValue::from("bottom up"));
    settings.insert("imzml:scan_direction".into(), MetaValue::from("meander"));
    settings.insert(
        "imzml:line_scan_direction".into(),
        MetaValue::from("right-left"),
    );
    settings.insert("imzml:polarity".into(), MetaValue::from("negative"));

    let meta = dataset_meta(&experiment).unwrap();
    assert_eq!(meta.imaging_mode, Some(ImagingMode::Continuous));
    assert_eq!(
        (meta.max_count_x, meta.max_count_y, meta.max_count_z),
        (4, 5, 6)
    );
    assert!((meta.pixel_size_x - 1.5).abs() < 1e-12);
    assert!((meta.pixel_size_y - 2.5).abs() < 1e-12);
    assert!((meta.max_dim_x - 6.0).abs() < 1e-12);
    assert!((meta.max_dim_y - 12.5).abs() < 1e-12);
    assert_eq!(meta.uuid, "{ABCD}");
    assert_eq!(meta.scan_pattern, "bottom up");
    assert_eq!(meta.scan_direction, "meander");
    assert_eq!(meta.line_scan_direction, "right-left");
    assert_eq!(meta.polarity, "negative");

    // An unrecognised mode is neither branch, so the mode is auto-detected.
    let mut other = experiment.clone();
    set_mode(&mut other, "sparse");
    assert_eq!(dataset_meta(&other).unwrap().imaging_mode, None);

    // Wrong types are refused.
    let mut wrong = MSExperiment::new();
    wrong
        .settings
        .metadata
        .insert("imzml:max_count_x".into(), MetaValue::from("four"));
    assert!(matches!(dataset_meta(&wrong), Err(Error::InvalidValue(_))));

    let mut wrong = MSExperiment::new();
    wrong
        .settings
        .metadata
        .insert("imzml:pixel_size_x".into(), MetaValue::from("wide"));
    assert!(matches!(dataset_meta(&wrong), Err(Error::InvalidValue(_))));

    // A vocabulary key is NOT refused for its type: `extractMeta_` reads all
    // six of them with the lenient `DataValue::toString()`, which never throws
    // and renders a number as its decimal text. The port used to error here.
    let mut lenient = MSExperiment::new();
    lenient
        .settings
        .metadata
        .insert("imzml:uuid".into(), MetaValue::from(7_i64));
    assert_eq!(dataset_meta(&lenient).unwrap().uuid, "7");

    // A negative count cannot be a pixel count; the source wraps it instead.
    let mut wrong = MSExperiment::new();
    wrong
        .settings
        .metadata
        .insert("imzml:max_count_y".into(), MetaValue::from(-1_i64));
    assert!(matches!(dataset_meta(&wrong), Err(Error::InvalidValue(_))));
}

/// `extractMeta_` reads all six vocabulary keys with `DataValue::toString()`,
/// the stringification `StringUtils.h` documents as the one that never throws.
/// So a non-string value for any of them is rendered, not refused: an integer
/// as its decimal text, a float through the same 15-digit rule the float
/// `cvParam`s use, a list joined as `[a, b]`, and an empty value as `""`.
#[test]
fn the_six_vocabulary_keys_take_any_type_the_way_the_source_does() {
    let mut experiment = MSExperiment::new();
    let settings = &mut experiment.settings.metadata;
    settings.insert("imzml:polarity".into(), MetaValue::from(1_i64));
    settings.insert(
        "imzml:uuid".into(),
        MetaValue::try_from(5.0_f64).expect("finite"),
    );
    settings.insert(
        "imzml:scan_pattern".into(),
        MetaValue::from(vec!["top".to_owned(), "down".to_owned()]),
    );
    settings.insert(
        "imzml:scan_direction".into(),
        MetaValue::from(vec![1_i64, 2]),
    );
    settings.insert(
        "imzml:line_scan_direction".into(),
        MetaValue::try_from(vec![1.5_f64, 100_000.0]).expect("a float list"),
    );
    settings.insert("imzml:imaging_mode".into(), MetaValue::default());

    let meta = dataset_meta(&experiment).expect("no vocabulary key can be refused");
    assert_eq!(meta.polarity, "1");
    // Not Rust's "5": the source's NumericFormatting keeps one fractional digit.
    assert_eq!(meta.uuid, "5.0");
    assert_eq!(meta.scan_pattern, "[top, down]");
    assert_eq!(meta.scan_direction, "[1, 2]");
    assert_eq!(meta.line_scan_direction, "[1.5, 1.0e05]");
    // An empty value stringifies to "", which is neither mode, so the mode is
    // auto-detected rather than refused.
    assert_eq!(meta.imaging_mode, None);

    // And such an experiment stores: "1" is not "positive", so no polarity term
    // is written, but nothing fails.
    let mut storable = experiment.clone();
    storable
        .spectra
        .push(pixel_spectrum(1, 1, vec![Peak1D::new(100.0, 1.0)]));
    let temp = TempDir::new(false).unwrap();
    let (path, report) = store_into(&temp, "lenient.imzML", &storable, &PeakFileOptions::new());
    assert_eq!(report.meta.polarity, "1");
    let document = xml(&path);
    assert!(!document.contains("MS:1000130"));
    assert!(!document.contains("MS:1000129"));
    // "5.0" is not 16 hex bytes, so the usual derivation replaces it rather
    // than the store failing.
    assert_ne!(report.meta.uuid, "5.0");
    assert_eq!(report.meta.uuid.len(), 36);
    assert!(document.contains(&format!(
        "name=\"universally unique identifier\" value=\"{}\"",
        report.meta.uuid
    )));
}

/// The numeric keys keep the source's strictness. `max_count_*` goes through
/// `static_cast<UInt>`, which throws for anything but a non-negative integer;
/// `pixel_size_*` and `max_dim_*` go through `static_cast<double>`, which
/// accepts an integer and throws for an empty value.
#[test]
fn the_numeric_keys_stay_strict_and_a_size_accepts_an_integer() {
    let mut counts = MSExperiment::new();
    counts.settings.metadata.insert(
        "imzml:max_count_x".into(),
        MetaValue::try_from(4.0_f64).expect("finite"),
    );
    assert!(matches!(dataset_meta(&counts), Err(Error::InvalidValue(_))));

    let mut sizes = MSExperiment::new();
    sizes
        .settings
        .metadata
        .insert("imzml:pixel_size_x".into(), MetaValue::from(3_i64));
    assert!((dataset_meta(&sizes).unwrap().pixel_size_x - 3.0).abs() < 1e-12);

    let mut empty = MSExperiment::new();
    empty
        .settings
        .metadata
        .insert("imzml:max_dim_x".into(), MetaValue::default());
    assert!(matches!(dataset_meta(&empty), Err(Error::InvalidValue(_))));
}

// ---------------------------------------------------------------------------
// Float cvParam text.
// ---------------------------------------------------------------------------

/// Every float `cvParam` goes through `StringConversions::toString(double)` →
/// `StringUtils::appendToStr` → `Internal::NumericFormatting::appendNumeric`
/// with 15 digits and `fixed_format = false`: fixed with 15 fractional digits
/// for `|v|` in `[1e-2, 1e4)`, shortest-round-trip scientific with the
/// exponent rewritten to the `e05` form otherwise, trailing zeros trimmed but
/// at least one digit kept after the dot. `std::ostringstream` — six
/// significant digits and no forced `.0` — is what the writer's own self-audit
/// claimed, and is not what the source does.
///
/// `imzml:pixel_size_x` is the probe because it is written verbatim, and
/// `imzml:max_count_x` is pinned to 1 so that `max_dim_x` is the pixel size
/// itself.
#[test]
fn a_float_cv_param_is_written_with_the_sources_numeric_formatting() {
    for (size, expected) in [
        // Just inside the fixed window at the bottom, and just outside it.
        (0.01_f64, "0.01"),
        (0.009_f64, "9.0e-03"),
        (0.009_999_f64, "9.999e-03"),
        // Just inside the fixed window at the top, and just outside it.
        (9999.0_f64, "9999.0"),
        (10000.0_f64, "1.0e04"),
        // A whole number needs the trailing ".0" the source forces.
        (5.0_f64, "5.0"),
        // An exponent that must render as e05, not e5 and not e+05.
        (100_000.0_f64, "1.0e05"),
        (123_456.789_f64, "1.23456789e05"),
        // All 15 fractional digits of the nearest double, which is what
        // separates this rule from a six-significant-digit stream.
        (9999.9_f64, "9999.899999999999636"),
    ] {
        let mut experiment = one_pixel(100.0, 1.0);
        let settings = &mut experiment.settings.metadata;
        settings.insert("imzml:max_count_x".into(), MetaValue::from(1_u32));
        settings.insert(
            "imzml:pixel_size_x".into(),
            MetaValue::try_from(size).expect("finite"),
        );

        let temp = TempDir::new(false).unwrap();
        let (path, _) = store_into(&temp, "floats.imzML", &experiment, &PeakFileOptions::new());
        let document = xml(&path);
        assert!(
            document.contains(&format!("name=\"pixel size x\" value=\"{expected}\"")),
            "{size} must be written as {expected}"
        );
        // max_dim_x is pixel size times a one-pixel grid, so the same text.
        assert!(
            document.contains(&format!("name=\"max dimension x\" value=\"{expected}\"")),
            "the recomputed extent of {size} must be written as {expected}"
        );
    }
}

/// The same rule on `MS:1000016` "scan start time", including the OpenMS unset
/// default of -1, which the source writes as "-1.0" and Rust's own `to_string`
/// would write as "-1".
#[test]
fn the_scan_start_time_is_written_with_the_same_rule() {
    for (rt, expected) in [
        (-1.0_f64, "-1.0"),
        (0.0_f64, "0.0"),
        (12.34_f64, "12.34"),
        (0.000_25_f64, "2.5e-04"),
        (86_400.0_f64, "8.64e04"),
    ] {
        let mut experiment = one_pixel(100.0, 1.0);
        experiment.spectra[0].rt = rt;
        let temp = TempDir::new(false).unwrap();
        let (path, _) = store_into(&temp, "rt.imzML", &experiment, &PeakFileOptions::new());
        assert!(
            xml(&path).contains(&format!("name=\"scan start time\" value=\"{expected}\"")),
            "rt {rt} must be written as {expected}"
        );
    }

    // `MSSpectrum::rt` is a plain public field, so a caller can put a
    // non-finite value there. `appendNumeric` has explicit NaN and infinity
    // branches and so does the port; neither may panic and neither may emit
    // the "inf.0" that once made such a file unreadable.
    for (rt, expected) in [
        (f64::NAN, "NaN"),
        (f64::INFINITY, "inf"),
        (f64::NEG_INFINITY, "-inf"),
    ] {
        let mut experiment = one_pixel(100.0, 1.0);
        experiment.spectra[0].rt = rt;
        let temp = TempDir::new(false).unwrap();
        let path = temp.path().join("nonfinite.imzML");
        match store(&path, &experiment, &PeakFileOptions::new()) {
            Ok(_) => assert!(
                xml(&path).contains(&format!("name=\"scan start time\" value=\"{expected}\"")),
                "rt {rt} must be written as {expected}"
            ),
            // An explicit refusal is acceptable; a panic is not.
            Err(Error::InvalidValue(_)) => {}
            Err(other) => panic!("unexpected error for rt {rt}: {other:?}"),
        }
    }
}

// ---------------------------------------------------------------------------
// Identifiers, checksums and escaping.
// ---------------------------------------------------------------------------

#[test]
fn a_declared_uuid_is_kept_verbatim_and_written_to_the_ibd_header() {
    let mut experiment = one_pixel(100.0, 1.0);
    experiment.settings.metadata.insert(
        "imzml:uuid".into(),
        MetaValue::from("12345678-1234-1234-1234-123456789012"),
    );

    let temp = TempDir::new(false).unwrap();
    let (path, report) = store_into(&temp, "uuid.imzML", &experiment, &PeakFileOptions::new());
    assert_eq!(report.meta.uuid, "12345678-1234-1234-1234-123456789012");

    let mut handler = ImzMLHandler::open(&path).unwrap();
    assert_eq!(handler.uuid_status().unwrap(), UuidStatus::Match);
    let header = std::fs::read(handler.ibd_path()).unwrap();
    assert_eq!(
        &header[..IBD_UUID_BYTES],
        &[
            0x12, 0x34, 0x56, 0x78, 0x12, 0x34, 0x12, 0x34, 0x12, 0x34, 0x12, 0x34, 0x56, 0x78,
            0x90, 0x12
        ]
    );
}

#[test]
fn an_absent_uuid_is_derived_deterministically_from_the_payload() {
    let temp = TempDir::new(false).unwrap();
    let experiment = one_pixel(100.0, 10.0);

    let (first_path, first) = store_into(
        &temp,
        "derive-a.imzML",
        &experiment,
        &PeakFileOptions::new(),
    );
    let (_, second) = store_into(
        &temp,
        "derive-b.imzML",
        &experiment,
        &PeakFileOptions::new(),
    );
    assert_eq!(first.meta.uuid, second.meta.uuid);
    assert!(!first.meta.uuid.is_empty());
    // RFC 4122 shape: 36 characters, version nibble 5, variant 10.
    assert_eq!(first.meta.uuid.len(), 36);
    assert_eq!(first.meta.uuid.as_bytes()[14], b'5');
    assert!(matches!(
        first.meta.uuid.as_bytes()[19],
        b'8' | b'9' | b'a' | b'b'
    ));

    // It is exactly derive_uuid over the bytes after the header.
    let mut handler = ImzMLHandler::open(&first_path).unwrap();
    assert_eq!(handler.uuid_status().unwrap(), UuidStatus::Match);
    let bytes = std::fs::read(handler.ibd_path()).unwrap();
    assert_eq!(
        &bytes[..IBD_UUID_BYTES],
        &derive_uuid(&bytes[IBD_UUID_BYTES..])
    );

    // A different payload gives a different identifier.
    let (_, other) = store_into(
        &temp,
        "derive-c.imzML",
        &one_pixel(100.0, 11.0),
        &PeakFileOptions::new(),
    );
    assert_ne!(first.meta.uuid, other.meta.uuid);
}

#[test]
fn an_unparsable_declared_uuid_falls_back_to_the_derivation() {
    let mut experiment = one_pixel(100.0, 1.0);
    experiment
        .settings
        .metadata
        .insert("imzml:uuid".into(), MetaValue::from("not-a-uuid"));
    let temp = TempDir::new(false).unwrap();
    let (path, report) = store_into(
        &temp,
        "bad-uuid.imzML",
        &experiment,
        &PeakFileOptions::new(),
    );
    assert_ne!(report.meta.uuid, "not-a-uuid");
    assert_eq!(report.meta.uuid.len(), 36);
    let mut handler = ImzMLHandler::open(&path).unwrap();
    assert_eq!(handler.uuid_status().unwrap(), UuidStatus::Match);
}

#[test]
fn derive_uuid_is_the_stamped_sha1_of_the_namespace_and_seed() {
    let first = derive_uuid(b"");
    let second = derive_uuid(b"");
    assert_eq!(first, second);
    assert_ne!(first, derive_uuid(b"x"));
    assert_eq!(first[6] & 0xf0, 0x50);
    assert_eq!(first[8] & 0xc0, 0x80);
    // The namespace is itself a version-5 identifier.
    assert_eq!(IBD_UUID_NAMESPACE[6] & 0xf0, 0x50);
    assert_eq!(IBD_UUID_NAMESPACE[8] & 0xc0, 0x80);
}

/// RFC 1321 appendix A.5, the published MD5 test suite.
#[test]
fn md5_matches_the_rfc_1321_test_suite() {
    for (input, expected) in [
        ("", "d41d8cd98f00b204e9800998ecf8427e"),
        ("a", "0cc175b9c0f1b6a831c399e269772661"),
        ("abc", "900150983cd24fb0d6963f7d28e17f72"),
        ("message digest", "f96b697d7cb7938d525a2f31aaf161d0"),
        (
            "abcdefghijklmnopqrstuvwxyz",
            "c3fcd3d76192e4007dfb496cca67e13b",
        ),
        (
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789",
            "d174ab98d277d9f5a5611c2c9f419d9f",
        ),
        (
            "12345678901234567890123456789012345678901234567890123456789012345678901234567890",
            "57edf4a22be3c955ac49da2e2107b67a",
        ),
    ] {
        assert_eq!(md5_hex(input.as_bytes()), expected, "MD5 of {input:?}");
    }
    // A multi-block input exercises the streaming path, not just the padding.
    let long = vec![b'x'; 4096];
    assert_eq!(md5_hex(&long).len(), 32);
}

#[test]
fn a_native_identifier_needing_escapes_round_trips() {
    let mut experiment = one_pixel(100.0, 1.0);
    experiment.spectra[0].native_id = "pixel=\"1\" & <2>".to_owned();
    let temp = TempDir::new(false).unwrap();
    let (path, _) = store_into(&temp, "escape.imzML", &experiment, &PeakFileOptions::new());
    assert_eq!(
        ImzMLHandler::open(&path).unwrap().index()[0].native_id,
        "pixel=\"1\" & <2>"
    );
}

#[test]
fn an_absent_native_identifier_gets_the_source_placeholder() {
    let mut experiment = MSExperiment::new();
    experiment
        .spectra
        .push(pixel_spectrum(1, 1, vec![Peak1D::new(100.0, 1.0)]));
    experiment
        .spectra
        .push(pixel_spectrum(2, 1, vec![Peak1D::new(100.0, 2.0)]));

    let temp = TempDir::new(false).unwrap();
    let (path, _) = store_into(&temp, "ids.imzML", &experiment, &PeakFileOptions::new());
    let handler = ImzMLHandler::open(&path).unwrap();
    assert_eq!(handler.index()[0].native_id, "spectrum=1");
    assert_eq!(handler.index()[1].native_id, "spectrum=2");
}

#[test]
fn both_reference_groups_declare_no_compression_and_external_data() {
    let temp = TempDir::new(false).unwrap();
    let (path, _) = store_into(
        &temp,
        "groups.imzML",
        &one_pixel(100.0, 10.0),
        &PeakFileOptions::new(),
    );
    let document = xml(&path);
    for id in ["mzArray", "intensityArray"] {
        let open = format!("<referenceableParamGroup id=\"{id}\">");
        let start = document.find(&open).expect("the group is written");
        let end = document[start..]
            .find("</referenceableParamGroup>")
            .expect("the group closes")
            + start;
        let group = &document[start..end];
        assert!(group.contains("MS:1000576"), "{id} needs no-compression");
        assert!(group.contains("IMS:1000101"), "{id} needs external data");
    }
}

#[test]
fn the_ibd_path_follows_the_imzml_name() {
    let temp = TempDir::new(false).unwrap();
    let experiment = one_pixel(100.0, 1.0);

    let (path, report) = store_into(&temp, "case.IMZML", &experiment, &PeakFileOptions::new());
    assert_eq!(report.meta.ibd_file_path, infer_ibd_path(&path));
    assert!(report.meta.ibd_file_path.exists());

    // Any other suffix simply gains .ibd.
    let other = temp.path().join("dataset.txt");
    let report = store(&other, &experiment, &PeakFileOptions::new()).unwrap();
    assert_eq!(
        report.meta.ibd_file_path,
        temp.path().join("dataset.txt.ibd")
    );
    assert!(report.meta.ibd_file_path.exists());
}

// ---------------------------------------------------------------------------
// Vocabulary resolution.
// ---------------------------------------------------------------------------

#[test]
fn resolve_float_array_cv_names_terms_units_and_free_text() {
    let cv = openms::format::controlled_vocabulary::ControlledVocabulary::psi_ms().unwrap();

    let mobility = resolve_float_array_cv("mean inverse reduced ion mobility array", cv).unwrap();
    assert_eq!(mobility.accession, "MS:1003006");
    assert_eq!(mobility.name, "mean inverse reduced ion mobility array");
    assert_eq!(mobility.unit_accession, "MS:1002814");
    assert_eq!(mobility.unit_cv_ref, "MS");
    assert!(!mobility.unit_name.is_empty());
    assert!(!mobility.non_standard);

    let pressure = resolve_float_array_cv("pressure array", cv).unwrap();
    assert_eq!(pressure.accession, "MS:1000821");
    assert_eq!(pressure.unit_accession, "UO:0000110");
    assert_eq!(pressure.unit_cv_ref, "UO");
    assert_eq!(pressure.unit_name, "pascal");

    // Two allowed units: none is written.
    let raw = resolve_float_array_cv("raw ion mobility array", cv).unwrap();
    assert_eq!(raw.accession, "MS:1003007");
    assert!(raw.unit_accession.is_empty());
    assert!(raw.unit_name.is_empty());
    assert!(raw.unit_cv_ref.is_empty());

    // The peak arrays are children of MS:1000513 too, which is why the writer
    // has to check for them by accession.
    assert_eq!(
        resolve_float_array_cv("m/z array", cv).unwrap().accession,
        "MS:1000514"
    );
    assert_eq!(
        resolve_float_array_cv("intensity array", cv)
            .unwrap()
            .accession,
        "MS:1000515"
    );

    let custom = resolve_float_array_cv("my custom SNR", cv).unwrap();
    assert_eq!(custom.accession, "MS:1000786");
    assert_eq!(custom.name, "non-standard data array");
    assert!(custom.non_standard);
    assert!(custom.unit_accession.is_empty());
}

// ---------------------------------------------------------------------------
// Resource ceilings. No upstream fixture reaches any of these.
// ---------------------------------------------------------------------------

fn write_with(limits: ImzMLWriteLimits) -> ImzMLWriteOptions {
    ImzMLWriteOptions {
        limits,
        ..ImzMLWriteOptions::default()
    }
}

fn refuses(name: &str, experiment: &MSExperiment, limits: ImzMLWriteLimits) {
    let temp = TempDir::new(false).unwrap();
    let path = temp.path().join(name);
    let mut logger = openms::concept::progress_logger::ProgressLogger::new();
    let outcome = store_with_options(
        &path,
        experiment,
        &PeakFileOptions::new(),
        &write_with(limits),
        &mut logger,
    );
    assert!(
        matches!(outcome, Err(Error::InvalidValue(_))),
        "{name} must be refused, got {outcome:?}"
    );
    assert!(!path.exists(), "{name} must leave no .imzML behind");
    assert!(
        !infer_ibd_path(&path).exists(),
        "{name} must leave no .ibd behind"
    );
}

#[test]
fn every_ceiling_refuses_before_a_file_is_created() {
    let mut two = MSExperiment::new();
    two.spectra
        .push(pixel_spectrum(1, 1, vec![Peak1D::new(100.0, 1.0)]));
    two.spectra
        .push(pixel_spectrum(2, 1, vec![Peak1D::new(100.0, 2.0)]));

    refuses(
        "spectra.imzML",
        &two,
        ImzMLWriteLimits {
            max_spectra: 1,
            ..ImzMLWriteLimits::default()
        },
    );
    refuses(
        "peaks.imzML",
        &two,
        ImzMLWriteLimits {
            max_peaks_per_spectrum: 0,
            ..ImzMLWriteLimits::default()
        },
    );
    refuses(
        "total-peaks.imzML",
        &two,
        ImzMLWriteLimits {
            max_total_peaks: 1,
            ..ImzMLWriteLimits::default()
        },
    );
    refuses(
        "ibd-bytes.imzML",
        &two,
        ImzMLWriteLimits {
            max_ibd_bytes: IBD_UUID_BYTES as u64,
            ..ImzMLWriteLimits::default()
        },
    );
    refuses(
        "checksum-bytes.imzML",
        &two,
        ImzMLWriteLimits {
            max_checksum_bytes: 1,
            ..ImzMLWriteLimits::default()
        },
    );
    refuses(
        "text-bytes.imzML",
        &two,
        ImzMLWriteLimits {
            max_text_bytes: 4,
            ..ImzMLWriteLimits::default()
        },
    );

    let mut long_id = one_pixel(100.0, 1.0);
    long_id.spectra[0].native_id = "n".repeat(64);
    refuses(
        "name-bytes.imzML",
        &long_id,
        ImzMLWriteLimits {
            max_name_bytes: 32,
            ..ImzMLWriteLimits::default()
        },
    );

    let mut aux = MSExperiment::new();
    let mut spectrum = pixel_spectrum(1, 1, vec![Peak1D::new(100.0, 1.0)]);
    spectrum
        .float_data_arrays
        .push(DataArray::new("pressure array", vec![1.0]));
    spectrum
        .float_data_arrays
        .push(DataArray::new("temperature array", vec![2.0]));
    aux.spectra.push(spectrum);
    refuses(
        "aux-per-spectrum.imzML",
        &aux,
        ImzMLWriteLimits {
            max_aux_arrays_per_spectrum: 1,
            ..ImzMLWriteLimits::default()
        },
    );
    refuses(
        "aux-total.imzML",
        &aux,
        ImzMLWriteLimits {
            max_total_aux_arrays: 1,
            ..ImzMLWriteLimits::default()
        },
    );
    refuses(
        "cv-lookups.imzML",
        &aux,
        ImzMLWriteLimits {
            max_cv_lookups: 1,
            ..ImzMLWriteLimits::default()
        },
    );
}

#[test]
fn the_report_listings_are_capped_while_the_counts_are_not() {
    let mut experiment = MSExperiment::new();
    for _ in 0..5 {
        experiment
            .spectra
            .push(pixel_spectrum(1, 1, vec![Peak1D::new(100.0, 1.0)]));
    }
    set_mode(&mut experiment, "processed");

    let temp = TempDir::new(false).unwrap();
    let mut logger = openms::concept::progress_logger::ProgressLogger::new();
    let report = store_with_options(
        temp.path().join("capped.imzML"),
        &experiment,
        &PeakFileOptions::new(),
        &write_with(ImzMLWriteLimits {
            max_reported_items: 2,
            ..ImzMLWriteLimits::default()
        }),
        &mut logger,
    )
    .unwrap();

    // Four spectra reuse pixel (1,1); only two are listed.
    assert_eq!(report.duplicate_pixel_count, 4);
    assert_eq!(report.duplicate_pixels.len(), 2);
}

#[test]
fn a_large_but_legal_dataset_is_accepted_at_the_default_ceilings() {
    // Nine pixels, 8399 peaks each: the upstream continuous fixture. It has to
    // pass at the defaults, or the ceilings are wrong rather than protective.
    let temp = TempDir::new(false).unwrap();
    let (_, report) = store_into(
        &temp,
        "defaults.imzML",
        &load_fixture(CONTINUOUS),
        &PeakFileOptions::new(),
    );
    assert_eq!(report.spectra_written, SPECTRA);
}
