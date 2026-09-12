// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The seven imzML access modes, and their cross-mode agreement.
//!
//! Ports all nine `START_SECTION`s of
//! `src/tests/class_tests/openms/source/ImzMLFile_all_modes_test.cpp`, in
//! upstream order, one Rust test per section. The 41 sections of
//! `ImzMLFile_test.cpp` are in `tests/imzml_file.rs`; the per-section table for
//! both suites is in `docs/IMZML_FILE_SUPPORT.md`.
//!
//! What this suite is for, and why it is worth porting separately from the
//! class test, is the last two sections: one spectrum is fetched through all
//! three access paths — the full in-memory load, the on-disc random access and
//! the RAM pixel lookup — and the three must return the same peaks. That is the
//! property the three-way split of this family exists to guarantee, and it is
//! the only place upstream where it is asserted.
//!
//! The literals are **tier 3 source review**: nine spectra, a 3 x 3 grid, the
//! `continuous` and `processed` modes, pixel (1,1) first in document order, and
//! the pixel (2,3) used for the cross-mode comparison. Whether the three paths
//! agree is not transcribed, it is computed here.

#![cfg(feature = "mzml")]

use openms::Error;
use openms::format::imzml_file::{ImzMLFile, build_imaging_geometry_from_experiment};
use openms::format::imzml_handler::ImagingMode;
use openms::format::{FileHandler, FileType, file_handler, file_types};
use openms::interfaces::MSDataConsumer;
use openms::kernel::on_disc_imzml_experiment::OnDiscImzMLExperiment;
use openms::kernel::{MSChromatogram, MSExperiment, MSSpectrum};
use openms::metadata::ExperimentalSettings;
use std::ops::ControlFlow;
use std::path::PathBuf;

const CONTINUOUS: &str = "ImzMLFile_1_Example_Continuous.imzML";
const PROCESSED: &str = "ImzMLFile_2_Example_Processed.imzML";
/// Upstream `k_continuous_spectra` and `k_grid`.
const CONTINUOUS_SPECTRA: usize = 9;
const GRID: u32 = 3;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// Upstream `loadImzMLExperiment_`.
fn load_experiment(name: &str) -> MSExperiment {
    ImzMLFile::new()
        .load_experiment(data(name))
        .expect("the fixture loads")
        .0
}

fn meta_text(experiment: &MSExperiment, key: &str) -> String {
    experiment.settings.metadata[key]
        .as_str()
        .unwrap_or_default()
        .to_owned()
}

fn meta_count(experiment: &MSExperiment, key: &str) -> u32 {
    u32::try_from(experiment.settings.metadata[key].as_i64().unwrap()).unwrap()
}

/// Upstream `findSpectrumIndexByImzMLCoord_`.
fn spectrum_index_at_imzml_coord(experiment: &MSExperiment, x: i64, y: i64) -> usize {
    experiment
        .spectra
        .iter()
        .position(|spectrum| {
            let metadata = &spectrum.metadata;
            metadata.get("imzml:x").and_then(|v| v.as_i64().ok()) == Some(x)
                && metadata.get("imzml:y").and_then(|v| v.as_i64().ok()) == Some(y)
        })
        .unwrap_or_else(|| panic!("no spectrum at ({x},{y})"))
}

/// Upstream `CollectConsumer`.
#[derive(Default)]
struct CollectConsumer {
    count: usize,
    first_size: usize,
}

impl MSDataConsumer for CollectConsumer {
    fn set_expected_size(&mut self, _: usize, _: usize) -> openms::Result<()> {
        Ok(())
    }
    fn set_experimental_settings(&mut self, _: &ExperimentalSettings) -> openms::Result<()> {
        Ok(())
    }
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> openms::Result<ControlFlow<()>> {
        if self.count == 0 {
            self.first_size = spectrum.peaks.len();
        }
        self.count += 1;
        Ok(ControlFlow::Continue(()))
    }
    fn consume_chromatogram(&mut self, _: &mut MSChromatogram) -> openms::Result<ControlFlow<()>> {
        panic!("imzML carries no chromatograms");
    }
}

// ===========================================================================
// Mode 1 — full load into MSExperiment
// ===========================================================================

#[test]
fn mode_1_a_full_load_gives_every_spectrum_and_the_dataset_metadata() {
    let experiment = load_experiment(CONTINUOUS);
    assert_eq!(experiment.spectra.len(), CONTINUOUS_SPECTRA);
    assert!(
        experiment
            .settings
            .metadata
            .contains_key("imzml:imaging_mode")
    );
    assert_eq!(meta_text(&experiment, "imzml:imaging_mode"), "continuous");
    assert_eq!(meta_count(&experiment, "imzml:max_count_x"), GRID);
    assert_eq!(meta_count(&experiment, "imzml:max_count_y"), GRID);
    assert!(!experiment.spectra[0].peaks.is_empty());
    assert_eq!(
        experiment.spectra[0].metadata["imzml:x"].as_i64().unwrap(),
        1
    );
    assert_eq!(
        experiment.spectra[0].metadata["imzml:y"].as_i64().unwrap(),
        1
    );
}

// ===========================================================================
// Mode 2 — FileHandler and FileTypes
// ===========================================================================

#[test]
fn mode_2_the_file_handler_recognises_imzml_and_refuses_to_load_it() {
    assert_eq!(
        file_types::type_by_file_name(CONTINUOUS),
        FileType::ImzMl,
        "recognised by extension"
    );
    let head = std::fs::read(data(CONTINUOUS)).unwrap();
    assert_eq!(
        file_handler::type_by_content(&head),
        FileType::ImzMl,
        "recognised by the IMS ontology markers in its header"
    );

    // The generic experiment loader must refuse imzML: the peaks are in a
    // second file it knows nothing about, so it would silently produce empty
    // spectra. Source throws Exception::InvalidFileType.
    assert!(!FileHandler::can_read_experiment(FileType::ImzMl));
    assert!(matches!(
        FileHandler::load_experiment(data(CONTINUOUS), &[]),
        Err(Error::Unsupported(_))
    ));
}

// ===========================================================================
// Mode 3 — streaming IMSDataConsumer
// ===========================================================================

#[test]
fn mode_3_a_streaming_load_reaches_the_consumer_once_per_spectrum() {
    let mut consumer = CollectConsumer::default();
    ImzMLFile::new()
        .load_into_consumer(data(CONTINUOUS), &mut consumer)
        .unwrap();
    assert_eq!(consumer.count, CONTINUOUS_SPECTRA);
    assert!(consumer.first_size > 0);
}

// ===========================================================================
// Mode 4 — loadSpectraIndex for on-disc access
// ===========================================================================

#[test]
fn mode_4_the_index_only_load_gives_the_offsets_and_the_dataset_metadata() {
    let index = ImzMLFile::new()
        .load_spectra_index(data(CONTINUOUS))
        .unwrap();
    assert_eq!(index.spectra.len(), CONTINUOUS_SPECTRA);
    assert_eq!(index.meta.imaging_mode, Some(ImagingMode::Continuous));
    assert_eq!(index.meta.max_count_x, GRID);
    assert_eq!(index.meta.max_count_y, GRID);
    assert_eq!(index.spectra[0].x, 1);
    assert_eq!(index.spectra[0].y, 1);
    assert!(index.spectra[0].mz_length > 0);
    // Native: the index records the .ibd it opened, which is what the on-disc
    // reader then reads from.
    assert_eq!(
        index.meta.ibd_file_path,
        openms::format::imzml_handler::infer_ibd_path(data(CONTINUOUS))
    );
}

// ===========================================================================
// Mode 5 — OnDiscImzMLExperiment random access
// ===========================================================================

#[test]
fn mode_5_the_on_disc_reader_serves_one_pixel_at_a_time() {
    let mut on_disc = OnDiscImzMLExperiment::new();
    on_disc.open(data(CONTINUOUS)).unwrap();

    assert!(on_disc.is_open());
    assert_eq!(on_disc.len(), CONTINUOUS_SPECTRA);
    assert_eq!(on_disc.grid_width(), GRID);
    assert_eq!(on_disc.grid_height(), GRID);
    assert_eq!(on_disc.meta().imaging_mode, Some(ImagingMode::Continuous));

    let by_index = on_disc.spectrum(0).unwrap();
    assert!(!by_index.peaks.is_empty());
    let entry = on_disc.index(0).unwrap().clone();
    let by_coord = on_disc
        .spectrum_at_coord(entry.x, entry.y, entry.z)
        .unwrap();
    assert_eq!(by_coord.peaks.len(), by_index.peaks.len());
    assert!(!by_coord.peaks.is_empty());
}

// ===========================================================================
// Mode 6 — MSImagingExperiment in-memory pixel lookup
// ===========================================================================

#[test]
fn mode_6_the_imaging_load_maps_every_spectrum_to_a_pixel() {
    let (imaging, report) = ImzMLFile::new().load(data(CONTINUOUS)).unwrap();
    assert_eq!(imaging.number_of_spectra(), CONTINUOUS_SPECTRA);
    assert_eq!(imaging.number_of_pixels(), CONTINUOUS_SPECTRA);
    assert_eq!(imaging.geometry().width(), GRID);
    assert_eq!(imaging.geometry().height(), GRID);
    assert!(imaging.has_pixel(0, 0));
    assert!(imaging.has_pixel(2, 2));
    assert!(!imaging.spectrum(0, 0).unwrap().peaks.is_empty());
    assert!(!imaging.spectrum(2, 2).unwrap().peaks.is_empty());
    // Native: nothing was left out, which is what makes pixels == spectra.
    assert!(report.geometry.is_clean());
}

// ===========================================================================
// Mode 7 — buildImagingGeometry from MSExperiment
// ===========================================================================

#[test]
fn mode_7_the_meta_value_geometry_builder_agrees_with_document_order() {
    let experiment = load_experiment(CONTINUOUS);
    let (geometry, report) = build_imaging_geometry_from_experiment(&experiment).unwrap();
    assert_eq!(geometry.number_of_pixels(), CONTINUOUS_SPECTRA);
    assert_eq!(geometry.spectrum_index(0, 0), Some(0));
    assert_eq!(
        geometry.spectrum_index(2, 2),
        Some(spectrum_index_at_imzml_coord(&experiment, 3, 3))
    );
    assert!(report.is_clean());
}

// ===========================================================================
// Cross-mode consistency — full load vs on-disc vs RAM lookup
// ===========================================================================

#[test]
fn the_three_access_paths_return_the_same_pixel() {
    let experiment = load_experiment(CONTINUOUS);
    let imaging = ImzMLFile::new().load(data(CONTINUOUS)).unwrap().0;
    let mut on_disc = OnDiscImzMLExperiment::new();
    on_disc.open(data(CONTINUOUS)).unwrap();

    let (px, py) = (2_i64, 3_i64);
    let position = spectrum_index_at_imzml_coord(&experiment, px, py);
    let reference = &experiment.spectra[position];
    let pz = reference
        .metadata
        .get("imzml:z")
        .and_then(|value| value.as_i64().ok())
        .unwrap_or(1);

    let full = reference.peaks.len();
    let ram = imaging
        .spectrum(
            u32::try_from(px - 1).unwrap(),
            u32::try_from(py - 1).unwrap(),
        )
        .unwrap();
    let disc = on_disc
        .spectrum_at_coord(
            u32::try_from(px).unwrap(),
            u32::try_from(py).unwrap(),
            u32::try_from(pz).unwrap(),
        )
        .unwrap();

    assert!(full > 0);
    assert_eq!(ram.peaks.len(), full);
    assert_eq!(disc.peaks.len(), full);

    // Beyond the upstream section, which compares only the sizes: the values
    // must agree too, or three paths of the same size could still disagree.
    for (index, peak) in reference.peaks.iter().enumerate() {
        assert_eq!(ram.peaks[index].mz.to_bits(), peak.mz.to_bits(), "{index}");
        assert_eq!(disc.peaks[index].mz.to_bits(), peak.mz.to_bits(), "{index}");
        assert_eq!(
            ram.peaks[index].intensity.to_bits(),
            peak.intensity.to_bits(),
            "{index}"
        );
        assert_eq!(
            disc.peaks[index].intensity.to_bits(),
            peak.intensity.to_bits(),
            "{index}"
        );
    }
}

// ===========================================================================
// Processed imzML encoding
// ===========================================================================

#[test]
fn a_processed_dataset_works_on_the_in_memory_and_the_on_disc_path() {
    let experiment = load_experiment(PROCESSED);
    assert!(!experiment.spectra.is_empty());
    assert_eq!(meta_text(&experiment, "imzml:imaging_mode"), "processed");
    assert!(!experiment.spectra[0].peaks.is_empty());

    let mut on_disc = OnDiscImzMLExperiment::new();
    on_disc.open(data(PROCESSED)).unwrap();
    assert_eq!(on_disc.meta().imaging_mode, Some(ImagingMode::Processed));
    assert!(!on_disc.spectrum(0).unwrap().peaks.is_empty());

    // Native: in processed mode every pixel owns its m/z array, so no two
    // spectra may share an offset — the property that distinguishes the two
    // modes in the .ibd.
    let mut offsets: Vec<u64> = (0..on_disc.len())
        .map(|index| on_disc.index(index).unwrap().mz_offset)
        .collect();
    offsets.sort_unstable();
    offsets.dedup();
    assert_eq!(offsets.len(), on_disc.len());
}
