// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Coverage for `KERNEL/OnDiscImzMLExperiment.h`.
//!
//! The header has no class test of its own upstream; it is exercised through
//! `ImzMLFile_test.cpp` and `ImzMLFile_all_modes_test.cpp`, which construct an
//! `OnDiscImzMLExperiment` twenty-two times between them. The literals
//! transcribed below are that suite's (tier 3): the continuous fixture's nine
//! spectra on a 3 x 3 grid with `getNumberOfPixels() == getNrSpectra()` and
//! `getSpectrumIndex(0, 0) == 0`, the `continuous` and `processed` modes, the
//! duplicate-pixel contract (`size() == 2`, one pixel, spectrum 0 keeps it),
//! and the sorted-external-peaks case whose on-disc ion image at pixel (0,0)
//! is 1310.0 for a spectrum stored as 131.0/121.0 in that order.
//!
//! One check does not depend on that suite: each extracted ion image is
//! compared against a sum computed here directly from
//! `ImzMLHandler::spectrum`, which is the same independent-summation
//! cross-check the C++ section performs against the in-memory
//! `MSImagingExperiment` path — unavailable in Rust, because
//! `IMAGING/MSImagingExperiment.h` is unported.
//!
//! The resource ceilings, the region algebra, the `IonImage` bounds and the
//! error-variant choices are independently derived (tier 4): no upstream
//! fixture reaches them.
//!
//! Section accounting lives in `docs/ON_DISC_IMZML_SUPPORT.md`: of the 22
//! upstream sections that construct an `OnDiscImzMLExperiment`, 17 are ported
//! here, 2 are partial and 3 are writer-driven and cannot be reproduced without
//! `ImzMLWriter`. None is unaccounted for.

#![cfg(feature = "mzml")]

use openms::Error;
use openms::format::imzml_handler::{
    AuxSkipReason, ImagingMode, ImzMLDataType, ImzMLHandler, ImzMLReadLimits, UuidStatus,
};
use openms::kernel::on_disc_imzml_experiment::{
    GeometryReport, ImagingGeometry, ImagingRegion, IonImage, MAX_IMAGE_PIXELS, MAX_REGIONS,
    OnDiscImzMLExperiment, RegionShape, build_imaging_geometry,
};
use openms::system::file::TempDir;
use std::path::{Path, PathBuf};

const CONTINUOUS: &str = "ImzMLFile_1_Example_Continuous.imzML";
const PROCESSED: &str = "ImzMLFile_2_Example_Processed.imzML";

/// UUID of the synthetic datasets below, in the hyphenated spelling the
/// upstream continuous fixture uses.
const UUID_TEXT: &str = "12345678-1234-1234-1234-123456789012";
const UUID_BYTES: [u8; 16] = [
    0x12, 0x34, 0x56, 0x78, 0x12, 0x34, 0x12, 0x34, 0x12, 0x34, 0x12, 0x34, 0x56, 0x78, 0x90, 0x12,
];

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn opened(name: &str) -> OnDiscImzMLExperiment {
    let mut experiment = OnDiscImzMLExperiment::new();
    experiment
        .open(data(name))
        .expect("upstream imzML fixture opens");
    experiment
}

/// A minimal well-formed `.imzML` around `body`, in the shape of the upstream
/// fixtures: the IMS terms the reader looks at and nothing else.
fn document(body: &str) -> String {
    format!(
        concat!(
            "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n",
            "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">",
            "<cvParam accession=\"IMS:1000080\" name=\"universally unique identifier\" value=\"{uuid}\"/>",
            "{body}",
            "</mzML>\n"
        ),
        uuid = UUID_TEXT,
        body = body
    )
}

/// One external float32 `<binaryDataArray>` addressing `.ibd` bytes.
fn array(role: &str, offset: u64, length: u64) -> String {
    typed_array(role, "MS:1000521", "MS:1000576", offset, length, "")
}

/// The same, with the binary data type, the compression term and an optional
/// `unitAccession` under the caller's control.
fn typed_array(
    role: &str,
    data_type: &str,
    compression: &str,
    offset: u64,
    length: u64,
    unit: &str,
) -> String {
    let unit = if unit.is_empty() {
        String::new()
    } else {
        format!(" unitAccession=\"{unit}\"")
    };
    format!(
        concat!(
            "<binaryDataArray encodedLength=\"0\">",
            "<cvParam accession=\"{role}\" name=\"array\"{unit}/>",
            "<cvParam accession=\"{data_type}\" name=\"binary type\"/>",
            "<cvParam accession=\"{compression}\" name=\"compression\"/>",
            "<cvParam accession=\"IMS:1000101\" name=\"external data\" value=\"true\"/>",
            "<cvParam accession=\"IMS:1000102\" name=\"external offset\" value=\"{offset}\"/>",
            "<cvParam accession=\"IMS:1000103\" name=\"external array length\" value=\"{length}\"/>",
            "<binary/>",
            "</binaryDataArray>"
        ),
        role = role,
        unit = unit,
        data_type = data_type,
        compression = compression,
        offset = offset,
        length = length
    )
}

/// One `<spectrum>` at pixel (1,1) with two float32 peaks and one auxiliary
/// external array described by `aux`.
fn spectrum_with_aux(aux: &str) -> String {
    document(&format!(
        concat!(
            "<run><spectrumList count=\"1\">",
            "<spectrum id=\"s=1\" index=\"0\" defaultArrayLength=\"0\">",
            "<scanList count=\"1\"><scan>",
            "<cvParam accession=\"IMS:1000050\" name=\"position x\" value=\"1\"/>",
            "<cvParam accession=\"IMS:1000051\" name=\"position y\" value=\"1\"/>",
            "</scan></scanList>",
            "<binaryDataArrayList count=\"3\">{mz}{int}{aux}</binaryDataArrayList>",
            "</spectrum>",
            "</spectrumList></run>"
        ),
        mz = array("MS:1000514", 16, 2),
        int = array("MS:1000515", 24, 2),
        aux = aux
    ))
}

/// One `<spectrum>` at 1-based pixel `(x, y)` whose peaks live in the `.ibd`.
fn spectrum(id: &str, x: u32, y: u32, mz_offset: u64, int_offset: u64, length: u64) -> String {
    format!(
        concat!(
            "<spectrum id=\"{id}\" index=\"0\" defaultArrayLength=\"0\">",
            "<scanList count=\"1\"><scan>",
            "<cvParam accession=\"IMS:1000050\" name=\"position x\" value=\"{x}\"/>",
            "<cvParam accession=\"IMS:1000051\" name=\"position y\" value=\"{y}\"/>",
            "</scan></scanList>",
            "<binaryDataArrayList count=\"2\">{mz}{int}</binaryDataArrayList>",
            "</spectrum>"
        ),
        id = id,
        x = x,
        y = y,
        mz = array("MS:1000514", mz_offset, length),
        int = array("MS:1000515", int_offset, length)
    )
}

/// Write a synthetic two-file dataset: `xml` plus a `.ibd` holding the UUID
/// header, then `mz` as float32 at offset 16, then `intensity` as float32.
fn write_dataset(dir: &Path, xml: &str, mz: &[f32], intensity: &[f32]) -> PathBuf {
    write_arrays(dir, xml, &[mz, intensity])
}

/// The same, for a dataset with more than the two peak arrays; the arrays are
/// laid out back to back after the 16-byte UUID header.
fn write_arrays(dir: &Path, xml: &str, arrays: &[&[f32]]) -> PathBuf {
    write_arrays_and_raw(dir, xml, arrays, &[])
}

/// The same, with `raw` bytes appended after the float32 arrays, for an array
/// whose declared binary data type is not float32.
fn write_arrays_and_raw(dir: &Path, xml: &str, arrays: &[&[f32]], raw: &[u8]) -> PathBuf {
    let imzml = dir.join("synthetic.imzML");
    std::fs::write(&imzml, xml).unwrap();
    let mut ibd = Vec::from(UUID_BYTES);
    for value in arrays.iter().copied().flatten() {
        ibd.extend_from_slice(&value.to_le_bytes());
    }
    ibd.extend_from_slice(raw);
    std::fs::write(dir.join("synthetic.ibd"), &ibd).unwrap();
    imzml
}

fn temp_dir() -> TempDir {
    TempDir::new_in(std::env::temp_dir(), false).unwrap()
}

/// The sum an ion image pixel must hold, computed straight from the reader
/// rather than from the grid walk under test.
fn reference_sum(name: &str, index: usize, mz: f64, tolerance_ppm: f64) -> f64 {
    let mut handler = ImzMLHandler::open(data(name)).unwrap();
    let decoded = handler.spectrum(index).unwrap();
    let dm = mz * tolerance_ppm * 1e-6;
    decoded
        .spectrum
        .peaks
        .iter()
        .filter(|peak| peak.mz >= mz - dm && peak.mz <= mz + dm)
        .map(|peak| f64::from(peak.intensity))
        .sum()
}

// ---------------------------------------------------------------------------
// Opening, state and metadata
// ---------------------------------------------------------------------------

#[test]
fn a_default_experiment_is_closed_and_empty() {
    let experiment = OnDiscImzMLExperiment::new();
    assert!(!experiment.is_open());
    assert_eq!(experiment.len(), 0);
    assert!(experiment.is_empty());
    assert_eq!(experiment.grid_width(), 0);
    assert_eq!(experiment.grid_height(), 0);
    assert_eq!(experiment.geometry().number_of_pixels(), 0);
    assert!(experiment.uuid_status().is_none());
    assert!(experiment.geometry_report().is_clean());
    assert_eq!(experiment.imzml_path(), Path::new(""));
}

#[test]
fn the_continuous_fixture_opens_with_the_upstream_grid() {
    // ImzMLFile_all_modes_test.cpp mode 5: nine spectra, a 3x3 grid, mode
    // "continuous", and the first pixel reachable both ways.
    let mut experiment = opened(CONTINUOUS);
    assert!(experiment.is_open());
    assert_eq!(experiment.len(), 9);
    assert_eq!(experiment.grid_width(), 3);
    assert_eq!(experiment.grid_height(), 3);
    assert_eq!(
        experiment.meta().imaging_mode,
        Some(ImagingMode::Continuous)
    );
    assert_eq!(experiment.uuid_status(), Some(&UuidStatus::Match));
    assert_eq!(
        experiment.ibd_path(),
        data("ImzMLFile_1_Example_Continuous.ibd")
    );

    let first = experiment.index(0).unwrap().clone();
    let by_index = experiment.spectrum(0).unwrap();
    let by_coord = experiment
        .spectrum_at_coord(first.x, first.y, first.z)
        .unwrap();
    assert!(!by_index.peaks.is_empty());
    assert_eq!(by_coord.peaks.len(), by_index.peaks.len());
}

#[test]
fn the_processed_fixture_opens_and_decodes() {
    // ImzMLFile_all_modes_test.cpp "processed imzML encoding".
    let mut experiment = opened(PROCESSED);
    assert_eq!(experiment.meta().imaging_mode, Some(ImagingMode::Processed));
    assert!(!experiment.is_empty());
    assert!(!experiment.spectrum(0).unwrap().peaks.is_empty());
    // Processed mode declares no IMS:1000052, and an absent param means z == 1,
    // so every pixel still lands on the addressable plane.
    assert_eq!(experiment.index(0).unwrap().z, 1);
    assert!(experiment.geometry_report().is_clean());
}

#[test]
fn an_explicit_ibd_override_is_honoured() {
    // ImzMLFile_test.cpp "[EXTRA] OnDiscImzMLExperiment::open honours an
    // explicit .ibd path override": copy the .imzML somewhere whose inferred
    // sibling does not exist, then open against the real .ibd.
    let dir = temp_dir();
    let copy = dir.path().join("relocated.imzML");
    std::fs::copy(data(PROCESSED), &copy).unwrap();

    let mut inferred = OnDiscImzMLExperiment::new();
    assert!(inferred.open(&copy).is_err());

    let mut experiment = OnDiscImzMLExperiment::new();
    experiment
        .open_with_ibd(&copy, data("ImzMLFile_2_Example_Processed.ibd"))
        .unwrap();
    assert!(!experiment.is_empty());
    assert!(!experiment.spectrum(0).unwrap().peaks.is_empty());
    assert_eq!(experiment.imzml_path(), copy.as_path());
}

#[test]
fn a_failed_open_leaves_the_previous_dataset_in_place() {
    // Native divergence: source open() replaces its Impl on the first line, so
    // a failed open there discards the loaded dataset. This commits only on
    // success.
    let mut experiment = opened(CONTINUOUS);
    let error = experiment.open(data("no_such_dataset.imzML")).unwrap_err();
    assert!(matches!(error, Error::Io(_)));
    assert!(experiment.is_open());
    assert_eq!(experiment.len(), 9);
    assert_eq!(experiment.geometry().number_of_pixels(), 9);
}

#[test]
fn close_releases_the_ibd_but_keeps_the_index() {
    let mut experiment = opened(CONTINUOUS);
    experiment.close();

    assert!(!experiment.is_open());
    // Source close() clears geometry_ and retains index_ and meta_.
    assert_eq!(experiment.len(), 9);
    assert_eq!(experiment.grid_width(), 3);
    assert_eq!(experiment.index(8).unwrap().x, 3);
    assert_eq!(experiment.geometry().number_of_pixels(), 0);
    assert_eq!(experiment.geometry().width(), 0);

    let error = experiment.spectrum(0).unwrap_err();
    match error {
        Error::Io(io) => assert_eq!(io.kind(), std::io::ErrorKind::NotFound),
        other => panic!("expected the source's FileNotFound equivalent, got {other:?}"),
    }
    assert!(matches!(
        experiment.spectrum_at_pixel(1, 1),
        Err(Error::Io(_))
    ));
    assert!(matches!(
        experiment.extract_ion_image(100.0, 10.0),
        Err(Error::Io(_))
    ));

    // close() is idempotent and noexcept-equivalent.
    experiment.close();
    assert_eq!(experiment.len(), 9);
}

#[test]
fn reopening_replaces_the_dataset() {
    let mut experiment = opened(CONTINUOUS);
    experiment.open(data(PROCESSED)).unwrap();
    assert_eq!(experiment.meta().imaging_mode, Some(ImagingMode::Processed));
    assert!(experiment.is_open());
}

// ---------------------------------------------------------------------------
// Index and coordinate access
// ---------------------------------------------------------------------------

#[test]
fn index_entries_carry_offsets_without_reading_the_ibd() {
    let experiment = opened(CONTINUOUS);
    let entry = experiment.index(0).unwrap();
    assert_eq!((entry.x, entry.y, entry.z), (1, 1, 1));
    assert_eq!(entry.mz_offset, 16);
    assert!(entry.mz_length > 0);
    // Continuous mode: every pixel names the one shared m/z array.
    for position in 0..experiment.len() {
        assert_eq!(experiment.index(position).unwrap().mz_offset, 16);
    }
}

#[test]
fn an_out_of_range_spectrum_index_is_an_error() {
    let mut experiment = opened(CONTINUOUS);
    // Source getIndex/getSpectrum throw Exception::IndexOverflow.
    assert!(matches!(experiment.index(9), Err(Error::InvalidValue(_))));
    assert!(matches!(
        experiment.spectrum(9),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        experiment.decoded_spectrum(usize::MAX),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn every_pixel_of_the_upstream_grid_is_addressable() {
    let mut experiment = opened(CONTINUOUS);
    for y in 1..=3 {
        for x in 1..=3 {
            let spectrum = experiment.spectrum_at_pixel(x, y).unwrap();
            assert!(!spectrum.peaks.is_empty());
            let expected = experiment.geometry().spectrum_index(x - 1, y - 1).unwrap();
            assert_eq!(
                spectrum.metadata.get("imzml:x").unwrap().to_string(),
                x.to_string()
            );
            assert_eq!(
                experiment.index(expected).unwrap().x,
                x,
                "grid lookup and index disagree at ({x},{y})"
            );
        }
    }
}

#[test]
fn only_the_first_plane_is_addressable_by_coordinate() {
    let mut experiment = opened(CONTINUOUS);
    // The grid is 2-D; the source answers ElementNotFound for any other z and
    // for a coordinate below 1.
    assert!(matches!(
        experiment.spectrum_at_coord(1, 1, 2),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        experiment.spectrum_at_coord(0, 1, 1),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        experiment.spectrum_at_coord(1, 0, 1),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        experiment.spectrum_at_coord(4, 4, 1),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn a_decoded_spectrum_carries_its_pixel_and_is_sorted() {
    let mut experiment = opened(CONTINUOUS);
    let decoded = experiment.decoded_spectrum(4).unwrap();
    let entry = experiment.index(4).unwrap();
    assert_eq!(
        decoded
            .spectrum
            .metadata
            .get("imzml:x")
            .unwrap()
            .to_string(),
        entry.x.to_string()
    );
    assert_eq!(
        decoded
            .spectrum
            .metadata
            .get("imzml:y")
            .unwrap()
            .to_string(),
        entry.y.to_string()
    );
    assert_eq!(
        decoded
            .spectrum
            .metadata
            .get("imzml:z")
            .unwrap()
            .to_string(),
        "1"
    );
    // Source decodeSpectrum() ends with sortByPosition().
    assert!(decoded.spectrum.is_sorted());
    assert!(!decoded.inline_peaks);
}

// ---------------------------------------------------------------------------
// Geometry construction from the index
// ---------------------------------------------------------------------------

#[test]
fn the_grid_matches_the_dataset_and_maps_every_spectrum() {
    // ImzMLFile_test.cpp "const MSImagingGeometry& getGeometry() const".
    let experiment = opened(CONTINUOUS);
    let geometry = experiment.geometry();
    assert_eq!(geometry.width(), experiment.grid_width());
    assert_eq!(geometry.height(), experiment.grid_height());
    assert_eq!(geometry.number_of_pixels(), experiment.len());

    let first = experiment.index(0).unwrap();
    assert!(geometry.has_pixel(first.x - 1, first.y - 1));
    assert_eq!(geometry.spectrum_index(first.x - 1, first.y - 1), Some(0));
    assert_eq!(geometry.spectrum_index(2, 2), Some(8));
    assert_eq!(geometry.spectrum_index(3, 0), None);
    assert!(!geometry.is_empty());

    // IMS:1000046/47 is 100 um on this fixture, in micrometres.
    assert!((geometry.pixel_size_x() - 100.0).abs() < 1e-9);
    assert!((geometry.pixel_size_y() - 100.0).abs() < 1e-9);
    assert_eq!(geometry.pixel_size_unit(), "micrometer");
}

#[test]
fn a_duplicate_pixel_keeps_the_first_spectrum() {
    // ImzMLFile_test.cpp "OnDiscImzMLExperiment tolerates duplicate pixel
    // coordinates by default": both spectra load, one pixel is mapped, and it
    // belongs to spectrum 0.
    let dir = temp_dir();
    let body = format!(
        "<run><spectrumList count=\"2\">{}{}</spectrumList></run>",
        spectrum("s=1", 1, 1, 16, 24, 2),
        spectrum("s=2", 1, 1, 16, 24, 2)
    );
    let path = write_dataset(dir.path(), &document(&body), &[100.0, 200.0], &[10.0, 20.0]);

    let mut experiment = OnDiscImzMLExperiment::new();
    experiment.open(&path).unwrap();
    assert_eq!(experiment.len(), 2);
    assert_eq!(experiment.geometry().number_of_pixels(), 1);
    assert_eq!(experiment.geometry().spectrum_index(0, 0), Some(0));

    // The dropped spectrum is counted and named rather than logged.
    let report = experiment.geometry_report();
    assert!(!report.is_clean());
    assert_eq!(report.duplicate_count, 1);
    assert_eq!(report.duplicate_pixels, vec![1]);
    assert_eq!(report.out_of_grid_count, 0);
    assert_eq!(report.non_positive_count, 0);

    // It stays reachable by index, which is what the source guarantees.
    assert!(!experiment.spectrum(1).unwrap().peaks.is_empty());
}

#[test]
fn a_pixel_outside_the_declared_grid_is_dropped_not_fatal() {
    let index = vec![
        entry(1, 1, 1),
        entry(9, 9, 1),
        entry(0, 1, 1),
        entry(2, 2, 2),
    ];
    let meta = openms::format::imzml_handler::ImzMLMeta {
        max_count_x: 3,
        max_count_y: 3,
        ..Default::default()
    };
    let (geometry, report) = build_imaging_geometry(&index, &meta).unwrap();
    assert_eq!(geometry.number_of_pixels(), 1);
    assert_eq!(geometry.spectrum_index(0, 0), Some(0));
    assert_eq!(report.out_of_grid, vec![1]);
    assert_eq!(report.non_positive_coordinates, vec![2]);
    assert_eq!(report.other_plane_count, 1);
    assert!(!report.is_clean());
}

#[test]
fn an_undeclared_grid_is_derived_from_the_pixels() {
    // Source: width stays 0 until the loop finishes, then becomes max_x + 1.
    let index = vec![entry(1, 1, 1), entry(4, 2, 1)];
    let (geometry, report) =
        build_imaging_geometry(&index, &openms::format::imzml_handler::ImzMLMeta::default())
            .unwrap();
    assert_eq!((geometry.width(), geometry.height()), (4, 2));
    assert_eq!(geometry.number_of_pixels(), 2);
    assert!(report.is_clean());
    // No IMS:1000046/47, so the default 1.0 um pixel stands.
    assert!((geometry.pixel_size_x() - 1.0).abs() < 1e-9);
}

#[test]
fn an_empty_index_builds_an_empty_grid() {
    let (geometry, report) =
        build_imaging_geometry(&[], &openms::format::imzml_handler::ImzMLMeta::default()).unwrap();
    assert_eq!((geometry.width(), geometry.height()), (0, 0));
    assert!(geometry.is_empty());
    assert!(report.is_clean());
    assert_eq!(GeometryReport::default(), report);
}

#[test]
fn a_declared_grid_above_the_ceiling_is_refused() {
    // Native bound: both dimensions come from the file and size the ion-image
    // allocation, so the product is checked before anything is allocated.
    let meta = openms::format::imzml_handler::ImzMLMeta {
        max_count_x: 100_000,
        max_count_y: 100_000,
        ..Default::default()
    };
    let error = build_imaging_geometry(&[entry(1, 1, 1)], &meta).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)));

    let dir = temp_dir();
    let body = format!(
        concat!(
            "<cvParam accession=\"IMS:1000042\" name=\"max count of pixels x\" value=\"100000\"/>",
            "<cvParam accession=\"IMS:1000043\" name=\"max count of pixels y\" value=\"100000\"/>",
            "<run><spectrumList count=\"1\">{}</spectrumList></run>"
        ),
        spectrum("s=1", 1, 1, 16, 20, 1)
    );
    let path = write_dataset(dir.path(), &document(&body), &[100.0], &[10.0]);
    let mut experiment = OnDiscImzMLExperiment::new();
    let error = experiment.open(&path).unwrap_err();
    match &error {
        Error::InvalidValue(message) => assert!(
            message.contains("above the ceiling"),
            "the grid ceiling should be the reason: {message}"
        ),
        other => panic!("expected Error::InvalidValue, got {other:?}"),
    }
    // The refusal is atomic: nothing was committed.
    assert!(!experiment.is_open());
    assert_eq!(experiment.len(), 0);
}

fn entry(x: u32, y: u32, z: u32) -> openms::format::imzml_handler::ImzMLSpectrumIndex {
    openms::format::imzml_handler::ImzMLSpectrumIndex {
        x,
        y,
        z,
        ..Default::default()
    }
}

// ---------------------------------------------------------------------------
// Ion image extraction
// ---------------------------------------------------------------------------

#[test]
fn a_whole_image_extraction_sums_the_window_at_every_pixel() {
    let mut experiment = opened(CONTINUOUS);
    let first = experiment.spectrum(0).unwrap();
    let mz = first.peaks[0].mz;
    let tolerance_ppm = 1000.0;

    let image = experiment.extract_ion_image(mz, tolerance_ppm).unwrap();
    assert_eq!((image.width(), image.height()), (3, 3));
    assert_eq!(image.data().len(), 9);
    assert_eq!(image.mask().len(), 9);

    let dm = mz * tolerance_ppm * 1e-6;
    assert!((image.mz_range().min().unwrap() - (mz - dm)).abs() < 1e-12);
    assert!((image.mz_range().max().unwrap() - (mz + dm)).abs() < 1e-12);

    // Every pixel of this fixture carries a spectrum, so every cell is valid,
    // and each must equal an independently computed sum.
    for y in 0..3 {
        for x in 0..3 {
            assert!(image.has_pixel(x, y), "pixel ({x},{y}) should be valid");
            let position = experiment.geometry().spectrum_index(x, y).unwrap();
            let expected = reference_sum(CONTINUOUS, position, mz, tolerance_ppm);
            let found = image.intensity(x, y).unwrap();
            assert!(
                (found - expected).abs() <= expected.abs() * 1e-9 + 1e-9,
                "pixel ({x},{y}): {found} != {expected}"
            );
        }
    }
}

#[test]
fn an_extraction_window_with_no_peak_is_valid_and_zero() {
    let mut experiment = opened(CONTINUOUS);
    // A window far above the fixture's m/z range: the pixels exist, so they are
    // marked valid with intensity 0 rather than left invalid.
    let image = experiment.extract_ion_image(50_000.0, 1.0).unwrap();
    for y in 0..3 {
        for x in 0..3 {
            assert!(image.has_pixel(x, y));
            assert_eq!(image.intensity(x, y).unwrap(), 0.0);
        }
    }
}

#[test]
fn a_zero_tolerance_window_is_a_point_query() {
    let mut experiment = opened(CONTINUOUS);
    let first = experiment.spectrum(0).unwrap();
    let mz = first.peaks[0].mz;
    // dm == 0, and the source's MZBegin/MZEnd pair is inclusive at both ends,
    // so the peak exactly at mz still contributes.
    let image = experiment.extract_ion_image(mz, 0.0).unwrap();
    assert!(image.intensity(0, 0).unwrap() > 0.0);
    assert!(image.mz_range().min().unwrap() == image.mz_range().max().unwrap());
}

#[test]
fn the_extraction_is_faithful_to_the_upstream_sorted_peak_case() {
    // ImzMLFile_test.cpp "void load sorts external peaks when
    // getSortSpectraByMZ is true": two peaks stored 131.0/121.0 in that order,
    // and the on-disc ion image at (0,0) for 131.0 +/- 10 ppm is 1310.0.
    let dir = temp_dir();
    let body = format!(
        "<run><spectrumList count=\"1\">{}</spectrumList></run>",
        spectrum("s=1", 1, 1, 16, 24, 2)
    );
    let path = write_dataset(
        dir.path(),
        &document(&body),
        &[131.0, 121.0],
        &[1310.0, 1210.0],
    );

    let mut experiment = OnDiscImzMLExperiment::new();
    experiment.open(&path).unwrap();
    let decoded = experiment.spectrum(0).unwrap();
    assert!(decoded.is_sorted());
    assert!((decoded.peaks[0].mz - 121.0).abs() < 1e-9);
    assert!((decoded.peaks[1].mz - 131.0).abs() < 1e-9);

    let image = experiment.extract_ion_image(131.0, 10.0).unwrap();
    assert!((image.intensity(0, 0).unwrap() - 1310.0).abs() < 1e-6);
    assert_eq!((image.width(), image.height()), (1, 1));
}

#[test]
fn a_negative_or_non_finite_extraction_argument_is_refused() {
    let mut experiment = opened(CONTINUOUS);
    // ImzMLFile_test.cpp asserts Exception::InvalidValue for both.
    for (mz, tolerance) in [
        (-1.0, 1000.0),
        (100.0, -1.0),
        (f64::NAN, 1.0),
        (100.0, f64::INFINITY),
    ] {
        assert!(
            matches!(
                experiment.extract_ion_image(mz, tolerance),
                Err(Error::InvalidValue(_))
            ),
            "mz={mz} tolerance={tolerance} should be refused"
        );
    }
}

#[test]
fn a_region_extraction_covers_only_the_region() {
    // ImzMLFile_test.cpp "IonImage extractIonImage(double, double, Size)":
    // a single-pixel rectangle at the origin.
    let mut experiment = opened(CONTINUOUS);
    let mz = experiment.spectrum(0).unwrap().peaks[0].mz;
    let tolerance_ppm = 1000.0;
    experiment
        .geometry_mut()
        .add_region(ImagingRegion::rectangle(1, "roi", 0, 0, 0, 0).unwrap())
        .unwrap();

    let image = experiment
        .extract_ion_image_in_region(mz, tolerance_ppm, 1)
        .unwrap();
    assert_eq!((image.width(), image.height()), (3, 3));
    assert!(image.has_pixel(0, 0));
    let expected = reference_sum(CONTINUOUS, 0, mz, tolerance_ppm);
    assert!((image.intensity(0, 0).unwrap() - expected).abs() <= expected.abs() * 1e-9 + 1e-9);
    // Every other pixel stays invalid, because the region did not cover it.
    for y in 0..3 {
        for x in 0..3 {
            if (x, y) != (0, 0) {
                assert!(!image.has_pixel(x, y), "pixel ({x},{y}) is outside the roi");
            }
        }
    }

    // An unknown region id must fail, where the source throws ElementNotFound.
    assert!(matches!(
        experiment.extract_ion_image_in_region(mz, tolerance_ppm, 99),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn a_masked_region_extraction_follows_the_bitmask() {
    let mut experiment = opened(CONTINUOUS);
    let mz = experiment.spectrum(0).unwrap().peaks[0].mz;
    // A checkerboard over the top-left 2x2 block.
    let region =
        ImagingRegion::from_mask(7, "mask", 0, 0, 2, 2, vec![true, false, false, true]).unwrap();
    experiment.geometry_mut().add_region(region).unwrap();

    let image = experiment
        .extract_ion_image_in_region(mz, 1000.0, 7)
        .unwrap();
    assert!(image.has_pixel(0, 0));
    assert!(image.has_pixel(1, 1));
    assert!(!image.has_pixel(1, 0));
    assert!(!image.has_pixel(0, 1));
}

#[test]
fn a_region_over_unacquired_pixels_extracts_nothing() {
    let mut experiment = opened(CONTINUOUS);
    experiment
        .geometry_mut()
        .add_region(ImagingRegion::rectangle(3, "outside", 50, 50, 60, 60).unwrap())
        .unwrap();
    let image = experiment
        .extract_ion_image_in_region(100.0, 10.0, 3)
        .unwrap();
    assert_eq!((image.width(), image.height()), (3, 3));
    assert!(image.mask().iter().all(|&valid| !valid));
}

// ---------------------------------------------------------------------------
// ImagingGeometry on its own
// ---------------------------------------------------------------------------

#[test]
fn a_default_geometry_has_unit_micrometer_pixels() {
    let geometry = ImagingGeometry::new();
    assert_eq!((geometry.width(), geometry.height()), (0, 0));
    assert_eq!(geometry.number_of_pixels(), 0);
    assert!((geometry.pixel_size_x() - 1.0).abs() < 1e-9);
    assert!((geometry.pixel_size_y() - 1.0).abs() < 1e-9);
    assert_eq!(geometry.pixel_size_unit(), "micrometer");
    assert_eq!(ImagingGeometry::NO_REGION, usize::MAX);
}

#[test]
fn pixels_keep_insertion_order_and_reject_duplicates_and_strays() {
    let mut geometry = ImagingGeometry::new();
    geometry.add_pixel(2, 0, 10).unwrap();
    geometry.add_pixel(0, 1, 20).unwrap();
    geometry.add_pixel(1, 1, 30).unwrap();
    assert_eq!(geometry.number_of_pixels(), 3);
    assert_eq!(geometry.pixels()[0].x, 2);
    assert_eq!(geometry.pixels()[0].spectrum_index, 10);
    assert_eq!(geometry.pixels()[2].y, 1);
    assert_eq!(geometry.spectrum_index(0, 1), Some(20));
    assert_eq!(geometry.spectrum_index(5, 5), None);

    // Source addPixel throws InvalidValue on a duplicate.
    assert!(matches!(
        geometry.add_pixel(2, 0, 99),
        Err(Error::InvalidValue(_))
    ));

    // The bounds test applies only once both dimensions are set.
    let mut bounded = ImagingGeometry::new();
    bounded.set_dimensions(3, 2).unwrap();
    bounded.add_pixel(2, 1, 0).unwrap();
    assert!(matches!(
        bounded.add_pixel(3, 0, 1),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        bounded.add_pixel(0, 2, 2),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn clear_restores_the_default_geometry() {
    let mut geometry = ImagingGeometry::new();
    geometry.set_dimensions(4, 4).unwrap();
    geometry.set_pixel_size(50.0, 75.0, "nanometer").unwrap();
    geometry.add_pixel(0, 0, 0).unwrap();
    geometry
        .add_region(ImagingRegion::rectangle(1, "r", 0, 0, 1, 1).unwrap())
        .unwrap();
    geometry.clear();

    assert_eq!((geometry.width(), geometry.height()), (0, 0));
    assert_eq!(geometry.number_of_pixels(), 0);
    assert!(!geometry.has_pixel(0, 0));
    assert!((geometry.pixel_size_x() - 1.0).abs() < 1e-9);
    assert_eq!(geometry.pixel_size_unit(), "micrometer");
    assert_eq!(geometry.number_of_regions(), 0);
}

#[test]
fn a_non_finite_pixel_size_is_refused() {
    let mut geometry = ImagingGeometry::new();
    assert!(matches!(
        geometry.set_pixel_size(f64::NAN, 1.0, "micrometer"),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        geometry.set_pixel_size(1.0, f64::INFINITY, "micrometer"),
        Err(Error::InvalidValue(_))
    ));
    // Refused before the mutation, so the default survives.
    assert!((geometry.pixel_size_x() - 1.0).abs() < 1e-9);
}

#[test]
fn regions_must_be_identified_and_disjoint() {
    let mut geometry = ImagingGeometry::new();
    geometry
        .add_region(ImagingRegion::rectangle(1, "left", 0, 0, 1, 3).unwrap())
        .unwrap();
    geometry
        .add_region(ImagingRegion::rectangle(2, "right", 2, 0, 3, 3).unwrap())
        .unwrap();
    assert_eq!(geometry.number_of_regions(), 2);
    assert!(geometry.has_region(2));
    assert_eq!(geometry.region(2).unwrap().name(), "right");
    assert_eq!(geometry.regions()[0].id(), 1);

    // The NO_REGION sentinel, a duplicate id and an overlap all fail.
    assert!(matches!(
        geometry.add_region(
            ImagingRegion::rectangle(ImagingGeometry::NO_REGION, "bad", 9, 9, 9, 9).unwrap()
        ),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        geometry.add_region(ImagingRegion::rectangle(1, "again", 8, 8, 8, 8).unwrap()),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        geometry.add_region(ImagingRegion::rectangle(3, "overlap", 1, 1, 2, 2).unwrap()),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(geometry.number_of_regions(), 2);

    // Removal reindexes the remaining regions.
    geometry.remove_region(1).unwrap();
    assert!(!geometry.has_region(1));
    assert_eq!(geometry.region(2).unwrap().name(), "right");
    assert!(matches!(
        geometry.remove_region(1),
        Err(Error::InvalidValue(_))
    ));
    geometry.clear_regions();
    assert_eq!(geometry.number_of_regions(), 0);
}

#[test]
fn region_membership_needs_an_acquired_pixel() {
    let mut geometry = ImagingGeometry::new();
    geometry.set_dimensions(4, 4).unwrap();
    geometry.add_pixel(0, 0, 5).unwrap();
    geometry.add_pixel(1, 1, 6).unwrap();
    geometry.add_pixel(3, 3, 7).unwrap();
    geometry
        .add_region(ImagingRegion::rectangle(4, "roi", 0, 0, 1, 1).unwrap())
        .unwrap();

    assert_eq!(geometry.region_pixels(4).unwrap(), vec![0, 1]);
    assert_eq!(geometry.region_spectrum_indices(4).unwrap(), vec![5, 6]);
    assert_eq!(geometry.region_of(0, 0), Some(4));
    assert_eq!(geometry.region_of(3, 3), None);
    // (1,0) is inside the footprint but carries no spectrum, so it belongs to
    // no region; the source returns its NO_REGION sentinel here.
    assert!(geometry.region(4).unwrap().contains(1, 0));
    assert_eq!(geometry.region_of(1, 0), None);
    assert!(matches!(
        geometry.region_pixels(99),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        geometry.region_spectrum_indices(99),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn the_region_count_is_bounded() {
    let mut geometry = ImagingGeometry::new();
    for id in 0..MAX_REGIONS {
        let x = u32::try_from(id).unwrap();
        geometry
            .add_region(ImagingRegion::rectangle(id + 1, "r", x, 0, x, 0).unwrap())
            .unwrap();
    }
    let x = u32::try_from(MAX_REGIONS).unwrap();
    assert!(matches!(
        geometry.add_region(ImagingRegion::rectangle(MAX_REGIONS + 1, "r", x, 0, x, 0).unwrap()),
        Err(Error::InvalidValue(_))
    ));
}

// ---------------------------------------------------------------------------
// ImagingRegion on its own
// ---------------------------------------------------------------------------

#[test]
fn a_rectangle_covers_its_inclusive_bounding_box() {
    let region = ImagingRegion::rectangle(1, "roi", 1, 2, 3, 5).unwrap();
    assert_eq!(region.id(), 1);
    assert_eq!(region.name(), "roi");
    assert_eq!(region.shape(), RegionShape::Rectangle);
    assert_eq!(
        (
            region.min_x(),
            region.min_y(),
            region.max_x(),
            region.max_y()
        ),
        (1, 2, 3, 5)
    );
    assert_eq!((region.bbox_width(), region.bbox_height()), (3, 4));
    assert_eq!(region.area(), 12);
    assert!(region.mask().is_empty());
    assert!(region.contains(1, 2));
    assert!(region.contains(3, 5));
    assert!(!region.contains(0, 2));
    assert!(!region.contains(3, 6));

    assert!(matches!(
        ImagingRegion::rectangle(1, "bad", 3, 0, 2, 0),
        Err(Error::InvalidRange(_))
    ));
    assert!(matches!(
        ImagingRegion::rectangle(1, "bad", 0, 3, 0, 2),
        Err(Error::InvalidRange(_))
    ));
}

#[test]
fn a_mask_covers_only_its_set_bits() {
    let region =
        ImagingRegion::from_mask(2, "m", 10, 20, 2, 2, vec![true, false, false, true]).unwrap();
    assert_eq!(region.shape(), RegionShape::Mask);
    assert_eq!(
        (
            region.min_x(),
            region.min_y(),
            region.max_x(),
            region.max_y()
        ),
        (10, 20, 11, 21)
    );
    assert_eq!(region.area(), 2);
    assert_eq!(region.mask().len(), 4);
    assert!(region.contains(10, 20));
    assert!(region.contains(11, 21));
    assert!(!region.contains(11, 20));
    assert!(!region.contains(12, 20));

    assert!(matches!(
        ImagingRegion::from_mask(2, "m", 0, 0, 0, 2, vec![]),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        ImagingRegion::from_mask(2, "m", 0, 0, 2, 2, vec![true, true]),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        ImagingRegion::from_mask(2, "m", 0, 0, 1, 2, vec![false, false]),
        Err(Error::InvalidValue(_))
    ));
    // Native bound: the mask is the allocation, so its cell count is checked.
    assert!(matches!(
        ImagingRegion::from_mask(2, "m", 0, 0, 70_000, 70_000, Vec::new()),
        Err(Error::InvalidValue(_))
    ));
    // Native check: the source computes origin + extent - 1 in wrapping UInt.
    assert!(matches!(
        ImagingRegion::from_mask(2, "m", u32::MAX, 0, 2, 1, vec![true, true]),
        Err(Error::InvalidRange(_))
    ));
}

#[test]
fn intersection_is_symmetric_and_mask_aware() {
    let left = ImagingRegion::rectangle(1, "a", 0, 0, 2, 2).unwrap();
    let right = ImagingRegion::rectangle(2, "b", 2, 2, 4, 4).unwrap();
    let far = ImagingRegion::rectangle(3, "c", 10, 10, 11, 11).unwrap();
    assert!(left.intersects(&right));
    assert!(right.intersects(&left));
    assert!(!left.intersects(&far));
    assert!(!far.intersects(&left));

    // Boxes that touch while the bitmasks do not.
    let checker_a =
        ImagingRegion::from_mask(4, "a", 0, 0, 2, 2, vec![true, false, false, true]).unwrap();
    let checker_b =
        ImagingRegion::from_mask(5, "b", 0, 0, 2, 2, vec![false, true, true, false]).unwrap();
    assert!(!checker_a.intersects(&checker_b));
    assert!(!checker_b.intersects(&checker_a));
    assert!(checker_a.intersects(&checker_a.clone()));
}

// ---------------------------------------------------------------------------
// IonImage on its own
// ---------------------------------------------------------------------------

#[test]
fn a_fresh_image_is_zeroed_and_fully_masked_out() {
    let image = IonImage::new(3, 2).unwrap();
    assert_eq!((image.width(), image.height()), (3, 2));
    assert_eq!(image.data(), &[0.0; 6]);
    assert!(image.mask().iter().all(|&valid| !valid));
    assert!(!image.has_pixel(0, 0));
    assert_eq!(image.intensity(2, 1).unwrap(), 0.0);
    assert!(image.mz_range().is_empty());
    assert_eq!(IonImage::default(), IonImage::new(0, 0).unwrap());
}

#[test]
fn writing_a_cell_marks_it_valid_and_is_row_major() {
    let mut image = IonImage::new(3, 2).unwrap();
    image.set_intensity(2, 1, 7.5).unwrap();
    assert!(image.has_pixel(2, 1));
    assert_eq!(image.intensity(2, 1).unwrap(), 7.5);
    // Row-major: index = y * width + x.
    assert_eq!(image.data()[5], 7.5);
    assert!(image.mask()[5]);
    assert!(!image.has_pixel(0, 0));
}

#[test]
fn image_coordinates_are_bounds_checked() {
    let mut image = IonImage::new(2, 2).unwrap();
    // Source linearIndex_ throws IndexOverflow; hasPixel answers false instead.
    assert!(!image.has_pixel(2, 0));
    assert!(!image.has_pixel(0, 2));
    assert!(matches!(image.intensity(2, 0), Err(Error::InvalidValue(_))));
    assert!(matches!(
        image.set_intensity(0, 2, 1.0),
        Err(Error::InvalidValue(_))
    ));
    // Native check: a non-finite pixel would compare unequal to itself.
    assert!(matches!(
        image.set_intensity(0, 0, f64::NAN),
        Err(Error::InvalidValue(_))
    ));
    let empty = IonImage::new(0, 0).unwrap();
    assert!(!empty.has_pixel(0, 0));
    assert!(matches!(empty.intensity(0, 0), Err(Error::InvalidValue(_))));
}

#[test]
fn resizing_zeroes_and_invalidates_every_cell() {
    let mut image = IonImage::new(2, 2).unwrap();
    image.set_intensity(1, 1, 3.0).unwrap();
    image.resize(1, 3).unwrap();
    assert_eq!((image.width(), image.height()), (1, 3));
    assert_eq!(image.data(), &[0.0; 3]);
    assert!(image.mask().iter().all(|&valid| !valid));
}

#[test]
fn an_image_above_the_pixel_ceiling_is_refused() {
    // Native bound; the source allocates width * height unconditionally.
    let side = u32::try_from(MAX_IMAGE_PIXELS).unwrap();
    assert!(matches!(
        IonImage::new(side, 2),
        Err(Error::InvalidValue(_))
    ));
    let mut image = IonImage::new(2, 2).unwrap();
    image.set_intensity(0, 0, 1.0).unwrap();
    assert!(matches!(image.resize(side, 2), Err(Error::InvalidValue(_))));
    // The refused resize left the image untouched.
    assert_eq!((image.width(), image.height()), (2, 2));
    assert_eq!(image.intensity(0, 0).unwrap(), 1.0);
}

// ---------------------------------------------------------------------------
// Bounds threaded through the façade
// ---------------------------------------------------------------------------

#[test]
fn the_configured_read_limits_reach_the_ibd() {
    let limits = ImzMLReadLimits {
        max_array_elements: 4,
        ..ImzMLReadLimits::default()
    };
    let mut experiment = OnDiscImzMLExperiment::new();
    experiment
        .open_with_limits(
            data(CONTINUOUS),
            data("ImzMLFile_1_Example_Continuous.ibd"),
            limits,
        )
        .unwrap();
    // Opening only parses the XML index; the 8399-element arrays are refused
    // when a decode actually asks for them.
    assert_eq!(experiment.len(), 9);
    assert!(matches!(
        experiment.spectrum(0),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        experiment.extract_ion_image(100.0, 10.0),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn a_hostile_offset_cannot_read_past_the_ibd() {
    let dir = temp_dir();
    let body = format!(
        "<run><spectrumList count=\"1\">{}</spectrumList></run>",
        spectrum("s=1", 1, 1, 16, u64::MAX - 8, 2)
    );
    let path = write_dataset(dir.path(), &document(&body), &[100.0, 200.0], &[1.0, 2.0]);

    let mut experiment = OnDiscImzMLExperiment::new();
    experiment.open(&path).unwrap();
    // The index parses; the range is refused before anything is allocated.
    assert_eq!(experiment.len(), 1);
    let error = experiment.spectrum(0).unwrap_err();
    assert!(matches!(
        error,
        Error::InvalidValue(_) | Error::Parse { .. }
    ));
    assert!(experiment.extract_ion_image(100.0, 1000.0).is_err());
}

#[test]
fn a_length_mismatch_between_the_peak_arrays_is_reported() {
    let dir = temp_dir();
    // Four float32 values follow the header; declare three m/z and one
    // intensity, so the decode finds unequal arrays.
    let body = format!(
        concat!(
            "<spectrum id=\"s=1\" index=\"0\" defaultArrayLength=\"0\">",
            "<scanList count=\"1\"><scan>",
            "<cvParam accession=\"IMS:1000050\" name=\"position x\" value=\"1\"/>",
            "<cvParam accession=\"IMS:1000051\" name=\"position y\" value=\"1\"/>",
            "</scan></scanList>",
            "<binaryDataArrayList count=\"2\">{mz}{int}</binaryDataArrayList>",
            "</spectrum>"
        ),
        mz = array("MS:1000514", 16, 3),
        int = array("MS:1000515", 28, 1)
    );
    let path = write_dataset(
        dir.path(),
        &document(&format!(
            "<run><spectrumList count=\"1\">{body}</spectrumList></run>"
        )),
        &[100.0, 200.0, 300.0],
        &[1.0],
    );

    let mut experiment = OnDiscImzMLExperiment::new();
    experiment.open(&path).unwrap();
    // Source decodePeaksInto_ raises ParseError naming the (x,y,z) pixel, on
    // both the getSpectrum and the extractIonImage path.
    let error = experiment.spectrum(0).unwrap_err();
    match &error {
        Error::Parse { message, .. } => assert!(
            message.contains("(1,1,1)"),
            "message should name the pixel: {message}"
        ),
        other => panic!("expected Error::Parse, got {other:?}"),
    }
    assert!(matches!(
        experiment.extract_ion_image(100.0, 1000.0),
        Err(Error::Parse { .. })
    ));
}

// ---------------------------------------------------------------------------
// Auxiliary external arrays and compression, through the facade
// ---------------------------------------------------------------------------

/// The header's `@note`: `getSpectrum()` returns any indexed auxiliary external
/// array as a float data array, named after its ontology term, so a viewer can
/// call `containsIMData()` without a separate code path.
#[test]
fn an_auxiliary_array_reaches_the_decoded_spectrum() {
    let dir = temp_dir();
    let aux = typed_array(
        "MS:1003006",
        "MS:1000521",
        "MS:1000576",
        32,
        2,
        "MS:1002814",
    );
    let path = write_arrays(
        dir.path(),
        &spectrum_with_aux(&aux),
        &[&[121.0, 131.0], &[1210.0, 1310.0], &[0.85, 0.95]],
    );

    let mut experiment = OnDiscImzMLExperiment::new();
    experiment.open(&path).unwrap();
    assert_eq!(experiment.index(0).unwrap().aux.len(), 1);

    let decoded = experiment.decoded_spectrum(0).unwrap();
    assert!(decoded.skipped_aux.is_empty());
    let arrays = &decoded.spectrum.float_data_arrays;
    assert_eq!(arrays.len(), 1);
    assert_eq!(arrays[0].name, "mean inverse reduced ion mobility array");
    assert_eq!(arrays[0].data.len(), 2);
    assert_eq!(
        arrays[0]
            .metadata
            .get("unit_accession")
            .unwrap()
            .to_string(),
        "MS:1002814"
    );
    assert!(decoded.spectrum.contains_im_data());
    // The sort that getSpectrum performs keeps the array aligned with the peaks:
    // the .ibd stores 121.0/131.0 in order, so 0.85 stays with 121.0.
    assert!((f64::from(arrays[0].data[0]) - 0.85).abs() < 1e-6);
}

/// `void load skips one bad aux length and still returns all spectra`: an
/// auxiliary array whose declared length is not the peak count is skipped with
/// a warning, and the spectrum still decodes.
#[test]
fn an_auxiliary_array_of_the_wrong_length_is_skipped() {
    let dir = temp_dir();
    let aux = typed_array("MS:1003006", "MS:1000521", "MS:1000576", 32, 1, "");
    let path = write_arrays(
        dir.path(),
        &spectrum_with_aux(&aux),
        &[&[121.0, 131.0], &[1210.0, 1310.0], &[0.85]],
    );

    let mut experiment = OnDiscImzMLExperiment::new();
    experiment.open(&path).unwrap();
    let decoded = experiment.decoded_spectrum(0).unwrap();
    assert_eq!(decoded.spectrum.peaks.len(), 2);
    assert!(decoded.spectrum.float_data_arrays.is_empty());
    assert!(!decoded.spectrum.contains_im_data());
    assert_eq!(
        decoded.skipped_aux[0].reason,
        AuxSkipReason::LengthMismatch {
            length: 1,
            peaks: 2
        }
    );
}

/// `void load and OnDisc skip aux array without a supported binary data type`:
/// MS:1000520 (obsolete 16-bit float) is not decodable by either loader.
#[test]
fn an_auxiliary_array_without_a_supported_type_is_skipped() {
    let dir = temp_dir();
    let aux = typed_array("MS:1003006", "MS:1000520", "MS:1000576", 32, 2, "");
    let path = write_arrays(
        dir.path(),
        &spectrum_with_aux(&aux),
        &[&[121.0, 131.0], &[1210.0, 1310.0], &[0.85, 0.95]],
    );

    let mut experiment = OnDiscImzMLExperiment::new();
    experiment.open(&path).unwrap();
    assert_eq!(experiment.index(0).unwrap().aux.len(), 1);
    assert_eq!(
        experiment.index(0).unwrap().aux[0].data_type,
        ImzMLDataType::Unknown
    );
    let decoded = experiment.decoded_spectrum(0).unwrap();
    assert_eq!(decoded.spectrum.peaks.len(), 2);
    assert!(decoded.spectrum.float_data_arrays.is_empty());
    assert_eq!(
        decoded.skipped_aux[0].reason,
        AuxSkipReason::UnknownDataType
    );
}

/// `void load and OnDisc drop zero-length aux without a ghost IM array`.
#[test]
fn a_zero_length_auxiliary_array_leaves_no_array() {
    let dir = temp_dir();
    let aux = typed_array("MS:1003006", "MS:1000521", "MS:1000576", 32, 0, "");
    let path = write_arrays(
        dir.path(),
        &spectrum_with_aux(&aux),
        &[&[121.0, 131.0], &[1210.0, 1310.0]],
    );

    let mut experiment = OnDiscImzMLExperiment::new();
    experiment.open(&path).unwrap();
    let decoded = experiment.decoded_spectrum(0).unwrap();
    assert!(decoded.spectrum.float_data_arrays.is_empty());
    assert!(!decoded.spectrum.contains_im_data());
    assert_eq!(decoded.skipped_aux[0].reason, AuxSkipReason::ZeroLength);
}

/// `OnDiscImzMLExperiment rejects compressed zero-length aux arrays`, and the
/// asymmetry the header documents: the extraction decodes only m/z and
/// intensity, so it never sees the compressed auxiliary array at all.
#[test]
fn a_compressed_auxiliary_array_fails_a_decode_but_not_an_extraction() {
    let dir = temp_dir();
    let aux = typed_array("MS:1003006", "MS:1000521", "MS:1000574", 32, 0, "");
    let path = write_arrays(
        dir.path(),
        &spectrum_with_aux(&aux),
        &[&[121.0, 131.0], &[1210.0, 1310.0]],
    );

    let mut experiment = OnDiscImzMLExperiment::new();
    experiment.open(&path).unwrap();
    assert!(experiment.index(0).unwrap().aux[0].compressed);
    assert_eq!(experiment.index(0).unwrap().aux[0].length, 0);
    // Source decodeSpectrum throws ParseError; a compressed array is the one
    // auxiliary condition that is an error rather than a skip.
    assert!(matches!(experiment.spectrum(0), Err(Error::Unsupported(_))));
    // Extraction sums peaks only, so the same dataset yields an image.
    let image = experiment.extract_ion_image(131.0, 10.0).unwrap();
    assert!((image.intensity(0, 0).unwrap() - 1310.0).abs() < 1e-6);
}

/// `void load and OnDisc reject zlib-compressed external m/z and intensity` and
/// `… reject numpress-compressed external m/z`: MS:1000576 is the only
/// acceptable value, so no compression term has to be enumerated.
#[test]
fn a_compressed_peak_array_is_refused_on_both_paths() {
    for compression in ["MS:1000574", "MS:1002312"] {
        let dir = temp_dir();
        let body = format!(
            concat!(
                "<run><spectrumList count=\"1\">",
                "<spectrum id=\"s=1\" index=\"0\" defaultArrayLength=\"0\">",
                "<scanList count=\"1\"><scan>",
                "<cvParam accession=\"IMS:1000050\" name=\"position x\" value=\"1\"/>",
                "<cvParam accession=\"IMS:1000051\" name=\"position y\" value=\"1\"/>",
                "</scan></scanList>",
                "<binaryDataArrayList count=\"2\">{mz}{int}</binaryDataArrayList>",
                "</spectrum>",
                "</spectrumList></run>"
            ),
            mz = typed_array("MS:1000514", "MS:1000521", compression, 16, 2, ""),
            int = array("MS:1000515", 24, 2)
        );
        let path = write_dataset(
            dir.path(),
            &document(&body),
            &[121.0, 131.0],
            &[1210.0, 1310.0],
        );

        let mut experiment = OnDiscImzMLExperiment::new();
        experiment.open(&path).unwrap();
        assert!(
            experiment.index(0).unwrap().mz_compressed,
            "{compression} should set the compression flag"
        );
        assert!(matches!(experiment.spectrum(0), Err(Error::Unsupported(_))));
        // Unlike an auxiliary array, a compressed peak array also fails the
        // extraction, because the sweep has to decode it.
        assert!(matches!(
            experiment.extract_ion_image(131.0, 10.0),
            Err(Error::Unsupported(_))
        ));
    }
}

/// `void load drops phantom IntegerDataArray for integer-typed aux`: an
/// auxiliary array declared `MS:1000519` (32-bit integer) is widened into a
/// float data array, and no integer data array appears.
#[test]
fn an_integer_typed_auxiliary_array_becomes_a_float_array() {
    let dir = temp_dir();
    let aux = typed_array("MS:1003006", "MS:1000519", "MS:1000576", 32, 2, "");
    let mut raw = Vec::new();
    for value in [7_i32, 9_i32] {
        raw.extend_from_slice(&value.to_le_bytes());
    }
    let path = write_arrays_and_raw(
        dir.path(),
        &spectrum_with_aux(&aux),
        &[&[121.0, 131.0], &[1210.0, 1310.0]],
        &raw,
    );

    let mut experiment = OnDiscImzMLExperiment::new();
    experiment.open(&path).unwrap();
    let decoded = experiment.decoded_spectrum(0).unwrap();
    assert!(decoded.spectrum.integer_data_arrays.is_empty());
    assert_eq!(decoded.spectrum.float_data_arrays.len(), 1);
    assert_eq!(decoded.spectrum.float_data_arrays[0].data, vec![7.0, 9.0]);
}

/// `void load and OnDisc drop inline aux when peaks are external`: an auxiliary
/// array with no `IMS:1000101` keeps its payload inline while the peaks come
/// from the `.ibd`, and neither loader decodes it.
#[test]
fn an_inline_auxiliary_array_is_not_decoded() {
    let dir = temp_dir();
    let aux = concat!(
        "<binaryDataArray encodedLength=\"0\">",
        "<cvParam accession=\"MS:1003006\" name=\"array\"/>",
        "<cvParam accession=\"MS:1000521\" name=\"binary type\"/>",
        "<cvParam accession=\"MS:1000576\" name=\"compression\"/>",
        "<binary/>",
        "</binaryDataArray>"
    );
    let path = write_arrays(
        dir.path(),
        &spectrum_with_aux(aux),
        &[&[121.0, 131.0], &[1210.0, 1310.0]],
    );

    let mut experiment = OnDiscImzMLExperiment::new();
    experiment.open(&path).unwrap();
    assert_eq!(experiment.index(0).unwrap().aux.len(), 0);
    let decoded = experiment.decoded_spectrum(0).unwrap();
    assert_eq!(decoded.spectrum.peaks.len(), 2);
    assert!(decoded.spectrum.float_data_arrays.is_empty());
    assert_eq!(decoded.skipped_aux[0].reason, AuxSkipReason::Inline);
}
