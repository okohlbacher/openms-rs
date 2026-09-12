// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Coverage for `FORMAT/MzDataFile.h` and `FORMAT/HANDLERS/MzDataHandler.h`.
//!
//! All fifteen `START_SECTION`s of `MzDataFile_test.cpp` are ported here. Every
//! literal that comes from that file or from the three unmodified upstream
//! fixtures is transcribed source review (tier 3): the three spectra of
//! `MzDataFile_1.mzData` with their retention times 60/120/180 and native IDs
//! `spectrum=10`/`11`/`12`, the eight annotation arrays `area` … `peakShape`,
//! the `lsid` accession number, the two contacts, the two mass analyzers with
//! resolutions 22.33 and 12.3, the sample `MS-Sample` `0-815` in the gas state,
//! the 997530-peak long spectrum, and the 64-bit and big-endian arrays of
//! `MzDataFile_4_64bit.mzData`.
//!
//! One check does not depend on that suite. The m/z payload of
//! `MzDataFile_3_minimal.mzData` and `MzDataFile_4_64bit.mzData` is 36 base64
//! characters, which is 25 bytes — one byte more than the three 64-bit values
//! the `length="3"` attribute declares. `Base64::decodeUncompressed_`
//! (`Base64.h:344`) assigns `s.size() / element_size` elements and so silently
//! drops the stray byte; that the upstream fixtures depend on this tolerance is
//! established here by decoding the payload independently, not by any
//! assertion upstream makes.
//!
//! The synthetic documents — the big-endian peak arrays, the ISO-8859-1 and
//! UTF-8 metadata, the resource ceilings, the hostile `length` attributes and
//! the writer refusals — are independently derived (tier 4), because no
//! upstream fixture reaches them.

#![cfg(feature = "mzml")]

use openms::Error;
use openms::format::PeakFileOptions;
use openms::format::mzdata::{
    Endian, LoadReport, MzDataFile, Precision, ReadLimits, SCHEMA_FILE, SCHEMA_VERSION,
    WriteOptions, load, load_with_options, read, read_with_options, store, store_with_options,
    write_with_options,
};
use openms::kernel::{MSExperiment, MSSpectrum, NumericRange, Peak1D, SpectrumType};
use openms::metadata::{
    ActivationMethod, AnalyzerType, DetectorAcquisitionMode, DetectorType, InletType,
    IonizationMethod, MetaValue, Polarity, ProcessingAction, ReflectronState, ResolutionMethod,
    ResolutionType, SampleState, ScanDirection, ScanLaw, ScanMode,
};
use openms::system::file::TempDir;
use std::path::PathBuf;

const FIXTURE_1: &str = "MzDataFile_1.mzData";
const FIXTURE_3_MINIMAL: &str = "MzDataFile_3_minimal.mzData";
const FIXTURE_4_64BIT: &str = "MzDataFile_4_64bit.mzData";
const FIXTURE_BIG_ENDIAN: &str = "mzdata_big_endian.mzData";
const FIXTURE_LATIN1: &str = "mzdata_latin1.mzData";
const FIXTURE_UTF8: &str = "mzdata_utf8.mzData";

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// `TOLERANCE_ABSOLUTE(0.01)`, the tolerance every `TEST_REAL_SIMILAR` in the
/// ported `load` sections runs under.
fn similar(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 0.01,
        "{actual} is not within 0.01 of {expected}"
    );
}

fn meta(value: Option<&MetaValue>) -> &str {
    value.expect("metadata entry is present").as_str().unwrap()
}

fn range(min: f64, max: f64) -> NumericRange {
    NumericRange { min, max }
}

/// A minimal well-formed mzData document around `body`, shaped like the
/// upstream fixtures' `<description>` block.
fn document(body: &str) -> String {
    format!(
        concat!(
            "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n",
            "<mzData version=\"1.05\" accessionNumber=\"synthetic\">\n",
            "\t<description>\n",
            "\t\t<admin>\n",
            "\t\t\t<sampleName></sampleName>\n",
            "\t\t</admin>\n",
            "\t\t<instrument>\n",
            "\t\t\t<instrumentName></instrumentName>\n",
            "\t\t</instrument>\n",
            "\t\t<dataProcessing>\n",
            "\t\t\t<software>\n",
            "\t\t\t\t<name></name>\n",
            "\t\t\t\t<version></version>\n",
            "\t\t\t</software>\n",
            "\t\t</dataProcessing>\n",
            "\t</description>\n",
            "{}",
            "</mzData>\n"
        ),
        body
    )
}

/// One `<spectrum>` with the two peak arrays spelled out.
fn spectrum_with(arrays: &str) -> String {
    document(&format!(
        concat!(
            "\t<spectrumList count=\"1\">\n",
            "\t\t<spectrum id=\"1\">\n",
            "\t\t\t<spectrumDesc>\n",
            "\t\t\t\t<spectrumSettings>\n",
            "\t\t\t\t\t<spectrumInstrument msLevel=\"1\"/>\n",
            "\t\t\t\t</spectrumSettings>\n",
            "\t\t\t</spectrumDesc>\n",
            "{}",
            "\t\t</spectrum>\n",
            "\t</spectrumList>\n"
        ),
        arrays
    ))
}

fn read_text(text: &str) -> Result<MSExperiment, Error> {
    read(text.as_bytes())
}

// ===========================================================================
// START_SECTION((MzDataFile())) and START_SECTION((~MzDataFile()))
// ===========================================================================

/// `MzDataFile_test.cpp:36-47`. `TEST_NOT_EQUAL(ptr, nullPointer)` asserts the
/// constructor produced an object; `~MzDataFile()` only deletes it. Rust has no
/// null adapter and no explicit destructor, so the equivalent statement is that
/// a default adapter exists, pins the source's schema version and file, and is
/// dropped without a leak at the end of the scope.
#[test]
fn default_construction_and_drop() {
    let file = MzDataFile::new();
    assert_eq!(file, MzDataFile::default());
    assert_eq!(file.version(), "1.05");
    assert_eq!(SCHEMA_VERSION, "1.05");
    assert_eq!(SCHEMA_FILE, "/SCHEMAS/mzData_1_05.xsd");
    assert!(!file.discards_unrepresentable());
    drop(file);
}

// ===========================================================================
// START_SECTION(const PeakFileOptions& getOptions() const)
// START_SECTION(setOptions(const PeakFileOptions & options))
// START_SECTION(PeakFileOptions& getOptions())
// ===========================================================================

/// `MzDataFile_test.cpp:49-79`, all three option sections:
/// `file.getOptions().hasMSLevels()` is false on a fresh adapter, stays false
/// through a const copy, becomes true after `setOptions` of an options object
/// carrying MS level 1, and becomes true when the mutable accessor is used to
/// add it in place.
#[test]
fn option_accessors() {
    let file = MzDataFile::new();
    assert!(!file.options().has_ms_levels());
    let copy = file.options().clone();
    assert!(!copy.has_ms_levels());

    let mut file = MzDataFile::new();
    assert!(!file.options().has_ms_levels());
    let mut options = file.options().clone();
    options.add_ms_level(1).unwrap();
    file.set_options(options);
    assert!(file.options().has_ms_levels());

    let mut file = MzDataFile::new();
    file.options_mut().add_ms_level(1).unwrap();
    assert!(file.options().has_ms_levels());
}

// ===========================================================================
// START_SECTION((template <typename MapType> void load(...)))
// MzDataFile_test.cpp:81-497 — 299 assertion macros, split across the six
// tests below. Each covers one commented block of that section.
// ===========================================================================

/// `MzDataFile_test.cpp:92-107`: the document identifier, the detected file
/// type, and the MS level, retention time and native ID of all three spectra.
#[test]
fn load_document_identity_and_scan_axis() {
    let file = MzDataFile::new();
    let loaded = file.load_report(data(FIXTURE_1)).unwrap();
    let e = &loaded.experiment;
    assert_eq!(
        e.settings.document.loaded_file_path,
        data(FIXTURE_1).to_str().unwrap()
    );
    assert_eq!(e.settings.document.loaded_file_type.name(), "mzData");
    assert_eq!(e.spectra.len(), 3);
    assert_eq!(e.spectra[0].ms_level, 1);
    assert_eq!(e.spectra[1].ms_level, 2);
    assert_eq!(e.spectra[2].ms_level, 1);
    similar(e.spectra[0].rt, 60.0);
    similar(e.spectra[1].rt, 120.0);
    similar(e.spectra[2].rt, 180.0);
    assert_eq!(e.spectra[0].native_id, "spectrum=10");
    assert_eq!(e.spectra[1].native_id, "spectrum=11");
    assert_eq!(e.spectra[2].native_id, "spectrum=12");
    assert_eq!(e.spectra[0].spectrum_type, SpectrumType::Unknown);
    // Every cvParam and userParam of this fixture is one the handler's fixed
    // vocabulary covers, and every declared `length` matches its payload, so
    // nothing is warned about.
    assert_eq!(loaded.report, LoadReport::default());
    assert!(loaded.report.is_clean());
}

/// `MzDataFile_test.cpp:113-136`: the `<supDesc>` descriptions of the eight
/// annotation arrays of the first spectrum, including the `comment` that comes
/// from the `<supDataDesc comment=…>` attribute rather than from a
/// `<userParam>`.
#[test]
fn load_annotation_array_descriptions() {
    let e = load(data(FIXTURE_1)).unwrap();
    let arrays = &e.spectra[0].float_data_arrays;
    assert_eq!(arrays.len(), 8);
    let comments = [
        "Area of the peak",
        "Full width at half max",
        "Left width",
        "Right width",
        "Peak charge",
        "Signal to noise ratio",
        "Correlation value",
        "Peak shape",
    ];
    for (array, expected) in arrays.iter().zip(comments) {
        assert_eq!(meta(array.metadata.get("URL")), "www.open-ms.de");
        assert_eq!(meta(array.metadata.get("Comment")), expected);
    }
    assert_eq!(meta(arrays[0].metadata.get("comment")), "bla|comment|bla");
    assert!(!arrays[1].metadata.contains_key("comment"));
    // The `<arrayName>` children, which the same section reaches through the
    // round trip in the store section.
    let names = [
        "area",
        "fwhm",
        "leftWidth",
        "rightWidth",
        "charge",
        "signalToNoise",
        "rValue",
        "peakShape",
    ];
    for (array, expected) in arrays.iter().zip(names) {
        assert_eq!(array.name, expected);
    }
}

/// `MzDataFile_test.cpp:141-159`: the two precursors of the second spectrum,
/// their m/z, charge, intensity, activation method, activation energy and the
/// `<userParam>`s of both `<ionSelection>` and `<activation>`.
#[test]
fn load_precursors() {
    let e = load(data(FIXTURE_1)).unwrap();
    assert_eq!(e.spectra[0].precursors.len(), 0);
    assert_eq!(e.spectra[1].precursors.len(), 2);
    assert_eq!(e.spectra[2].precursors.len(), 0);

    let first = &e.spectra[1].precursors[0];
    similar(first.mz, 1.2);
    assert_eq!(first.charge, 2);
    similar(f64::from(first.intensity), 2.3);
    assert_eq!(
        meta(first.cv_terms.metadata.get("IonSelectionComment")),
        "selected"
    );
    assert!(first.activation_methods.contains(&ActivationMethod::Cid));
    assert_eq!(first.activation_methods.len(), 1);
    similar(first.activation_energy, 3.4);
    assert_eq!(
        meta(first.cv_terms.metadata.get("ActivationComment")),
        "active"
    );

    let second = &e.spectra[1].precursors[1];
    similar(second.mz, 2.2);
    assert_eq!(second.charge, 3);
    similar(f64::from(second.intensity), 3.3);
    assert_eq!(
        meta(second.cv_terms.metadata.get("IonSelectionComment")),
        "selected2"
    );
    assert!(second.activation_methods.contains(&ActivationMethod::Sid));
    similar(second.activation_energy, 4.4);
    assert_eq!(
        meta(second.cv_terms.metadata.get("ActivationComment")),
        "active2"
    );
}

/// `MzDataFile_test.cpp:164-205`: the instrument settings and the acquisition
/// information of all three spectra, including the scan window whose
/// `mzRangeStop` is absent and therefore zero.
#[test]
fn load_instrument_settings_and_acquisition() {
    let e = load(data(FIXTURE_1)).unwrap();
    let settings: Vec<_> = e
        .spectra
        .iter()
        .map(|spectrum| &spectrum.instrument_settings)
        .collect();
    assert_eq!(meta(settings[0].metadata.get("URL")), "www.open-ms.de");
    assert_eq!(meta(settings[1].metadata.get("URL")), "www.open-ms.de");
    assert!(!settings[2].metadata.contains_key("URL"));
    assert_eq!(meta(settings[0].metadata.get("SpecComment")), "Spectrum 1");
    assert_eq!(meta(settings[1].metadata.get("SpecComment")), "Spectrum 2");
    assert!(!settings[2].metadata.contains_key("SpecComment"));
    assert_eq!(settings[0].scan_mode, ScanMode::MassSpectrum);
    assert_eq!(settings[1].scan_mode, ScanMode::MassSpectrum);
    assert_eq!(settings[2].scan_mode, ScanMode::SelectedIonMonitoring);
    assert_eq!(settings[0].polarity, Polarity::Positive);
    assert_eq!(settings[1].polarity, Polarity::Positive);
    assert_eq!(settings[2].polarity, Polarity::Negative);
    assert_eq!(settings[0].scan_windows.len(), 0);
    assert_eq!(settings[1].scan_windows.len(), 1);
    similar(settings[1].scan_windows[0].begin, 110.0);
    similar(settings[1].scan_windows[0].end, 0.0);
    assert_eq!(settings[2].scan_windows.len(), 1);
    similar(settings[2].scan_windows[0].begin, 100.0);
    similar(settings[2].scan_windows[0].end, 140.0);
    // That window is (110, 0), which `ScanWindow::validate` rejects because
    // its begin exceeds its end. The reader deliberately does not validate:
    // the source has no validation step on load, and refusing here would make
    // this port unable to read its own reference data.
    assert!(e.validate().is_err());

    assert_eq!(e.spectra[0].acquisition_info.acquisitions.len(), 0);
    assert_eq!(e.spectra[1].acquisition_info.acquisitions.len(), 2);
    assert_eq!(e.spectra[1].spectrum_type, SpectrumType::Profile);
    assert_eq!(e.spectra[1].acquisition_info.method_of_combination, "sum");
    let acquisitions = &e.spectra[1].acquisition_info.acquisitions;
    assert_eq!(acquisitions[0].identifier, "501");
    assert_eq!(acquisitions[1].identifier, "502");
    assert_eq!(meta(acquisitions[0].metadata.get("URL")), "www.open-ms.de");
    assert_eq!(meta(acquisitions[1].metadata.get("URL")), "www.open-ms.de");
    assert_eq!(
        meta(acquisitions[0].metadata.get("AcqComment")),
        "Acquisition 1"
    );
    assert_eq!(
        meta(acquisitions[1].metadata.get("AcqComment")),
        "Acquisition 2"
    );

    assert_eq!(e.spectra[2].acquisition_info.acquisitions.len(), 1);
    assert_eq!(e.spectra[2].spectrum_type, SpectrumType::Centroid);
    assert_eq!(
        e.spectra[2].acquisition_info.method_of_combination,
        "average"
    );
    assert_eq!(
        e.spectra[2].acquisition_info.acquisitions[0].identifier,
        "601"
    );
}

/// `MzDataFile_test.cpp:223-324`: the peaks of all three spectra and the eight
/// aligned annotation values behind each of them.
#[test]
fn load_peaks_and_annotation_values() {
    let e = load(data(FIXTURE_1)).unwrap();
    let expected: [&[(f64, f32)]; 3] = [
        &[(120.0, 100.0)],
        &[(110.0, 100.0), (120.0, 200.0), (130.0, 100.0)],
        &[
            (100.0, 100.0),
            (110.0, 200.0),
            (120.0, 300.0),
            (130.0, 200.0),
            (140.0, 100.0),
        ],
    ];
    assert_eq!(e.spectra[0].peaks.len(), 1);
    assert_eq!(e.spectra[1].peaks.len(), 3);
    assert_eq!(e.spectra[2].peaks.len(), 5);
    for (spectrum, points) in e.spectra.iter().zip(expected) {
        for (peak, &(mz, intensity)) in spectrum.peaks.iter().zip(points) {
            similar(peak.mz, mz);
            similar(f64::from(peak.intensity), f64::from(intensity));
        }
        // Every annotation array carries the peak's own intensity, so all
        // eight columns equal the intensity column, as the section asserts
        // value by value.
        for array in &spectrum.float_data_arrays {
            assert_eq!(array.data.len(), points.len());
            for (value, &(_, intensity)) in array.data.iter().zip(points) {
                similar(f64::from(*value), f64::from(intensity));
            }
        }
    }
}

/// `MzDataFile_test.cpp:329-437`: the accession number, the source file, the
/// contacts, the document-wide data processing, the instrument with its two
/// mass analyzers, and the sample.
#[test]
fn load_experiment_metadata() {
    let e = load(data(FIXTURE_1)).unwrap();
    assert_eq!(e.settings.document.identifier, "lsid");

    assert_eq!(e.settings.source_files.len(), 1);
    let file = &e.settings.source_files[0];
    assert_eq!(file.name, "MzDataFile_test_1.raw");
    assert_eq!(file.path, "/share/data/");
    assert_eq!(file.file_type, "MS");
    assert_eq!(file.checksum, "");
    assert_eq!(file.checksum_type, openms::metadata::ChecksumType::Unknown);

    assert_eq!(e.settings.contacts.len(), 2);
    assert_eq!(e.settings.contacts[0].first_name, "John");
    assert_eq!(e.settings.contacts[0].last_name, "Doe");
    assert_eq!(e.settings.contacts[0].institution, "department 1");
    assert_eq!(e.settings.contacts[0].contact_info, "www.john.doe");
    assert_eq!(e.settings.contacts[1].first_name, "Jane");
    assert_eq!(e.settings.contacts[1].last_name, "Doe");
    assert_eq!(e.settings.contacts[1].institution, "department 2");
    assert_eq!(e.settings.contacts[1].contact_info, "www.jane.doe");

    for spectrum in &e.spectra {
        assert_eq!(spectrum.data_processing.len(), 1);
        let processing = &spectrum.data_processing[0];
        assert_eq!(meta(processing.metadata.get("URL")), "www.open-ms.de");
        assert_eq!(
            meta(processing.metadata.get("comment")),
            "ProcessingComment"
        );
        assert_eq!(
            processing.completion_time.unwrap().get(),
            "2001-02-03 04:05:06"
        );
        assert_eq!(processing.software.name, "MS-X");
        assert_eq!(processing.software.version, "1.0");
        assert_eq!(
            meta(processing.software.cv_terms.metadata.get("comment")),
            "SoftwareComment"
        );
        assert!(processing.actions.contains(&ProcessingAction::Deisotoping));
        assert!(
            processing
                .actions
                .contains(&ProcessingAction::ChargeDeconvolution)
        );
    }

    let instrument = &e.settings.instrument;
    assert_eq!(instrument.name, "MS-Instrument");
    assert_eq!(instrument.vendor, "MS-Vendor");
    assert_eq!(instrument.model, "MS 1");
    assert_eq!(instrument.customizations, "tuned");
    assert_eq!(meta(instrument.metadata.get("URL")), "www.open-ms.de");
    assert_eq!(
        meta(instrument.metadata.get("AdditionalComment")),
        "Additional"
    );
    assert_eq!(instrument.ion_sources.len(), 1);
    assert_eq!(
        instrument.ion_sources[0].ionization_method,
        IonizationMethod::Esi
    );
    assert_eq!(instrument.ion_sources[0].inlet_type, InletType::Direct);
    assert_eq!(instrument.ion_sources[0].polarity, Polarity::Negative);
    assert_eq!(
        meta(instrument.ion_sources[0].metadata.get("URL")),
        "www.open-ms.de"
    );
    assert_eq!(
        meta(instrument.ion_sources[0].metadata.get("SourceComment")),
        "Source"
    );
    assert_eq!(instrument.ion_detectors.len(), 1);
    let detector = &instrument.ion_detectors[0];
    assert_eq!(detector.detector_type, DetectorType::FaradayCup);
    assert_eq!(detector.acquisition_mode, DetectorAcquisitionMode::Tdc);
    assert_eq!(detector.resolution, 0.815);
    assert_eq!(detector.adc_sampling_frequency, 11.22);
    assert_eq!(meta(detector.metadata.get("URL")), "www.open-ms.de");
    assert_eq!(meta(detector.metadata.get("DetectorComment")), "Detector");

    assert_eq!(instrument.mass_analyzers.len(), 2);
    let first = &instrument.mass_analyzers[0];
    assert_eq!(first.analyzer_type, AnalyzerType::PaulIonTrap);
    assert_eq!(first.resolution_method, ResolutionMethod::Fwhm);
    assert_eq!(first.resolution_type, ResolutionType::Constant);
    assert_eq!(first.scan_direction, ScanDirection::Up);
    assert_eq!(first.scan_law, ScanLaw::Linear);
    assert_eq!(first.reflectron_state, ReflectronState::Off);
    assert_eq!(first.resolution, 22.33);
    assert_eq!(first.accuracy, 33.44);
    assert_eq!(first.scan_rate, 44.55);
    assert_eq!(first.scan_time, 55.66);
    assert_eq!(first.tof_total_path_length, 66.77);
    assert_eq!(first.isolation_width, 77.88);
    assert_eq!(first.final_ms_exponent, 2);
    assert_eq!(first.magnetic_field_strength, 88.99);
    assert_eq!(meta(first.metadata.get("URL")), "www.open-ms.de");
    assert_eq!(meta(first.metadata.get("AnalyzerComment")), "Analyzer 1");
    let second = &instrument.mass_analyzers[1];
    assert_eq!(second.analyzer_type, AnalyzerType::Quadrupole);
    assert_eq!(second.resolution_method, ResolutionMethod::Baseline);
    assert_eq!(second.resolution_type, ResolutionType::Proportional);
    assert_eq!(second.scan_direction, ScanDirection::Down);
    assert_eq!(second.scan_law, ScanLaw::Exponential);
    assert_eq!(second.reflectron_state, ReflectronState::On);
    assert_eq!(second.resolution, 12.3);
    assert_eq!(second.accuracy, 13.4);
    assert_eq!(second.scan_rate, 14.5);
    assert_eq!(second.scan_time, 15.6);
    assert_eq!(second.tof_total_path_length, 16.7);
    assert_eq!(second.isolation_width, 17.8);
    assert_eq!(second.final_ms_exponent, -2);
    assert_eq!(second.magnetic_field_strength, 18.9);
    assert_eq!(meta(second.metadata.get("URL")), "www.open-ms.de");
    assert_eq!(meta(second.metadata.get("AnalyzerComment")), "Analyzer 2");

    let sample = &e.settings.sample;
    assert_eq!(sample.name, "MS-Sample");
    assert_eq!(sample.number, "0-815");
    assert_eq!(sample.state, SampleState::Gas);
    assert_eq!(sample.mass, 1.01);
    assert_eq!(sample.volume, 2.02);
    assert_eq!(sample.concentration, 3.03);
    assert_eq!(meta(sample.metadata.get("URL")), "www.open-ms.de");
    assert_eq!(meta(sample.metadata.get("SampleComment")), "Sample");
}

/// `MzDataFile_test.cpp:441-495`, the section's four special cases: loading the
/// same document twice yields equal experiments; the minimal fixture with
/// whitespace inside its base64 yields one spectrum of three peaks; the long
/// fixture yields one spectrum of 997530 peaks; and the 64-bit fixture yields
/// three peaks starting at m/z 110 with intensity 100.
///
/// The 10.6 MB upstream `MzDataFile_2_long.mzData` is not copied into this
/// repository — it would nearly double `tests/data` — so the long-spectrum case
/// writes an equivalent 997530-peak document and reads it back. Its hash is
/// recorded in `tests/data/mzdata_provenance.json` so the substitution is
/// visible.
#[test]
fn load_special_cases() {
    let file = MzDataFile::new();
    let first = file.load(data(FIXTURE_1)).unwrap();
    let second = file.load(data(FIXTURE_1)).unwrap();
    assert_eq!(first, second);

    let minimal = file.load_report(data(FIXTURE_3_MINIMAL)).unwrap();
    assert_eq!(minimal.experiment.spectra.len(), 1);
    assert_eq!(minimal.experiment.spectra[0].peaks.len(), 3);
    similar(minimal.experiment.spectra[0].peaks[0].mz, 110.0);
    similar(minimal.experiment.spectra[0].peaks[2].mz, 130.0);
    // Both of this fixture's payloads decode to 25 bytes, one more than the
    // three 64-bit values they declare; the source drops the stray byte
    // silently and this port warns about each.
    assert_eq!(minimal.report.warning_count, 2);
    assert!(
        minimal.report.warnings[0].contains("25 bytes"),
        "{:?}",
        minimal.report.warnings
    );

    let long = long_document(997_530);
    let loaded = read_text(&long).unwrap();
    assert_eq!(loaded.spectra.len(), 1);
    assert_eq!(loaded.spectra[0].peaks.len(), 997_530);

    let sixty_four = file.load_report(data(FIXTURE_4_64BIT)).unwrap();
    let e5 = &sixty_four.experiment;
    assert_eq!(e5.settings.document.identifier, "");
    assert_eq!(e5.spectra.len(), 1);
    assert_eq!(e5.spectra[0].peaks.len(), 3);
    let points = [(110.0, 100.0), (120.0, 200.0), (130.0, 100.0)];
    for (peak, (mz, intensity)) in e5.spectra[0].peaks.iter().zip(points) {
        similar(peak.mz, mz);
        similar(f64::from(peak.intensity), intensity);
    }
    assert_eq!(e5.spectra[0].float_data_arrays.len(), 8);
    for array in &e5.spectra[0].float_data_arrays {
        for (value, (_, intensity)) in array.data.iter().zip(points) {
            similar(f64::from(*value), intensity);
        }
    }
    // Only the 64-bit m/z payload carries the stray trailing byte here.
    assert_eq!(sixty_four.report.warning_count, 1);
}

/// The upstream fixtures' m/z payload really does hold one byte too many.
/// Decoded here without going through the reader, so this does not depend on
/// any assertion the class test makes.
#[test]
fn upstream_payload_holds_a_partial_trailing_element() {
    let text = std::fs::read_to_string(data(FIXTURE_3_MINIMAL)).unwrap();
    let payload: String = text
        .split("<data precision=\"64\" endian=\"little\" length=\"3\">")
        .nth(1)
        .unwrap()
        .split("</data>")
        .next()
        .unwrap()
        .chars()
        .filter(|c| !c.is_whitespace())
        .collect();
    assert_eq!(payload.len(), 36);
    // 36 base64 characters with two '=' pad characters carry 25 bytes, and
    // three 64-bit values need 24.
    assert_eq!(payload.len() / 4 * 3 - 2, 25);
    assert_eq!(25 / 8, 3);
}

/// `<mzArrayBinary>` and `<intenArrayBinary>` for `count` synthetic peaks, in
/// the layout the writer produces.
fn long_document(count: usize) -> String {
    let mut experiment = MSExperiment::new();
    experiment.spectra.push(MSSpectrum {
        rt: 1.0,
        native_id: "spectrum=1".into(),
        peaks: (0..count)
            .map(|index| Peak1D::new(100.0 + index as f64, 1.0))
            .collect(),
        ..Default::default()
    });
    let mut bytes = Vec::new();
    write_with_options(&mut bytes, &experiment, &WriteOptions::default()).unwrap();
    String::from_utf8(bytes).unwrap()
}

// ===========================================================================
// START_SECTION(([EXTRA] load with metadata - only flag))
// ===========================================================================

/// `MzDataFile_test.cpp:499-524`: with `setMetadataOnly(true)` the parse stops
/// at `<spectrumList>`, so no spectrum is retained while the source file, both
/// contacts, the instrument and the sample are.
#[test]
fn load_metadata_only() {
    let mut file = MzDataFile::new();
    file.options_mut().metadata_only = true;
    let e = file.load(data(FIXTURE_1)).unwrap();
    assert_eq!(e.spectra.len(), 0);
    assert_eq!(e.settings.source_files.len(), 1);
    assert_eq!(e.settings.source_files[0].name, "MzDataFile_test_1.raw");
    assert_eq!(e.settings.contacts.len(), 2);
    assert_eq!(e.settings.contacts[0].first_name, "John");
    assert_eq!(e.settings.contacts[0].last_name, "Doe");
    assert_eq!(e.settings.instrument.name, "MS-Instrument");
    assert_eq!(e.settings.instrument.vendor, "MS-Vendor");
    assert_eq!(e.settings.sample.name, "MS-Sample");
    assert_eq!(e.settings.sample.number, "0-815");
}

// ===========================================================================
// START_SECTION(([EXTRA] load with selected MS levels))
// ===========================================================================

/// `MzDataFile_test.cpp:526-555`: selecting MS level 1 keeps only
/// `spectrum=10` and `spectrum=12`, with 1 and 5 peaks; clearing the selection
/// restores all three with 1, 3 and 5 peaks.
#[test]
fn load_selected_ms_levels() {
    let mut file = MzDataFile::new();
    file.options_mut().add_ms_level(1).unwrap();
    let loaded = file.load_report(data(FIXTURE_1)).unwrap();
    let e = &loaded.experiment;
    assert_eq!(e.spectra.len(), 2);
    assert_eq!(e.spectra[0].peaks.len(), 1);
    assert_eq!(e.spectra[0].native_id, "spectrum=10");
    assert_eq!(e.spectra[1].peaks.len(), 5);
    assert_eq!(e.spectra[1].native_id, "spectrum=12");
    assert_eq!(e.spectra[0].ms_level, 1);
    assert_eq!(e.spectra[1].ms_level, 1);
    assert_eq!(loaded.report.spectra_skipped, 1);

    file.options_mut().clear_ms_levels();
    let e = file.load(data(FIXTURE_1)).unwrap();
    assert_eq!(e.spectra.len(), 3);
    assert_eq!(e.spectra[0].peaks.len(), 1);
    assert_eq!(e.spectra[1].peaks.len(), 3);
    assert_eq!(e.spectra[2].peaks.len(), 5);
    assert_eq!(e.spectra[0].ms_level, 1);
    assert_eq!(e.spectra[1].ms_level, 2);
    assert_eq!(e.spectra[2].ms_level, 1);
}

// ===========================================================================
// START_SECTION(([EXTRA] load with RT range))
// ===========================================================================

/// `MzDataFile_test.cpp:557-577`: the retention-time range [100, 200) keeps the
/// spectra at 120 and 180 seconds and drops the one at 60.
#[test]
fn load_rt_range() {
    let mut file = MzDataFile::new();
    file.options_mut().set_rt_range(range(100.0, 200.0));
    let loaded = file.load_report(data(FIXTURE_1)).unwrap();
    let e = &loaded.experiment;
    assert_eq!(e.spectra.len(), 2);
    assert_eq!(e.spectra[0].ms_level, 2);
    assert_eq!(e.spectra[1].ms_level, 1);
    similar(e.spectra[0].rt, 120.0);
    similar(e.spectra[1].rt, 180.0);
    assert_eq!(loaded.report.spectra_skipped, 1);
}

// ===========================================================================
// START_SECTION(([EXTRA] load with MZ range))
// ===========================================================================

/// `MzDataFile_test.cpp:579-614`: the m/z range [115, 135) keeps 1, 2 and 2 of
/// the 1, 3 and 5 peaks, and the retained coordinates are 120; 120 and 130;
/// 120 and 130.
#[test]
fn load_mz_range() {
    let mut file = MzDataFile::new();
    file.options_mut().set_mz_range(range(115.0, 135.0));
    let loaded = file.load_report(data(FIXTURE_1)).unwrap();
    let e = &loaded.experiment;
    assert_eq!(e.spectra.len(), 3);
    assert_eq!(e.spectra[0].peaks.len(), 1);
    assert_eq!(e.spectra[1].peaks.len(), 2);
    assert_eq!(e.spectra[2].peaks.len(), 2);
    similar(e.spectra[0].peaks[0].mz, 120.0);
    similar(f64::from(e.spectra[0].peaks[0].intensity), 100.0);
    similar(e.spectra[1].peaks[0].mz, 120.0);
    similar(f64::from(e.spectra[1].peaks[0].intensity), 200.0);
    similar(e.spectra[1].peaks[1].mz, 130.0);
    similar(f64::from(e.spectra[1].peaks[1].intensity), 100.0);
    similar(e.spectra[2].peaks[0].mz, 120.0);
    similar(f64::from(e.spectra[2].peaks[0].intensity), 300.0);
    similar(e.spectra[2].peaks[1].mz, 130.0);
    similar(f64::from(e.spectra[2].peaks[1].intensity), 200.0);
    assert_eq!(loaded.report.peaks_filtered_out, 4);
    // The annotation arrays are filtered with the same index set.
    for spectrum in &e.spectra {
        for array in &spectrum.float_data_arrays {
            assert_eq!(array.data.len(), spectrum.peaks.len());
        }
    }
}

// ===========================================================================
// START_SECTION(([EXTRA] load with intensity range))
// ===========================================================================

/// `MzDataFile_test.cpp:616-648`: the intensity range [150, 350) leaves the
/// first spectrum empty, keeps one peak of the second and three of the third,
/// the retained coordinates being 120; 110, 120 and 130.
#[test]
fn load_intensity_range() {
    let mut file = MzDataFile::new();
    file.options_mut().set_intensity_range(range(150.0, 350.0));
    let loaded = file.load_report(data(FIXTURE_1)).unwrap();
    let e = &loaded.experiment;
    assert_eq!(e.spectra.len(), 3);
    assert_eq!(e.spectra[0].peaks.len(), 0);
    assert_eq!(e.spectra[1].peaks.len(), 1);
    assert_eq!(e.spectra[2].peaks.len(), 3);
    similar(e.spectra[1].peaks[0].mz, 120.0);
    similar(f64::from(e.spectra[1].peaks[0].intensity), 200.0);
    similar(e.spectra[2].peaks[0].mz, 110.0);
    similar(f64::from(e.spectra[2].peaks[0].intensity), 200.0);
    similar(e.spectra[2].peaks[1].mz, 120.0);
    similar(f64::from(e.spectra[2].peaks[1].intensity), 300.0);
    similar(e.spectra[2].peaks[2].mz, 130.0);
    similar(f64::from(e.spectra[2].peaks[2].intensity), 200.0);
    assert_eq!(loaded.report.peaks_filtered_out, 5);
}

// ===========================================================================
// START_SECTION((template <typename MapType> void store(...)))
// ===========================================================================

/// `MzDataFile_test.cpp:650-667`: the fixture loads to three spectra, stores
/// and reloads to an experiment whose accession number is still `lsid` and
/// which equals the original once the software `comment` is put back.
///
/// Putting it back is what the upstream section does too, on all three spectra
/// (`:662-664`), because `writeTo` emits `<software><name>` and `<version>` and
/// nothing else — the `<comments>` child the reader accepts is never written.
/// This port refuses that loss by default, so the section runs under
/// [`WriteOptions::source`].
#[test]
fn store_round_trip_equals_the_loaded_experiment() {
    let file = MzDataFile::new();
    let original = file.load(data(FIXTURE_1)).unwrap();
    assert_eq!(original.spectra.len(), 3);

    let directory = TempDir::new(false).unwrap();
    let path = directory.path().join("round_trip.mzData");
    store_with_options(&path, &original, &WriteOptions::source()).unwrap();
    let mut reloaded = load(&path).unwrap();
    assert_eq!(reloaded.settings.document.identifier, "lsid");
    for spectrum in &mut reloaded.spectra {
        let processing = std::sync::Arc::make_mut(&mut spectrum.data_processing[0]);
        processing
            .software
            .cv_terms
            .metadata
            .insert("comment".into(), MetaValue::from("SoftwareComment"));
    }
    assert_eq!(original, reloaded);
}

/// The default writer refuses the same store, because it would drop the
/// software comment.
#[test]
fn store_refuses_to_drop_the_software_comment() {
    let original = load(data(FIXTURE_1)).unwrap();
    let directory = TempDir::new(false).unwrap();
    let path = directory.path().join("refused.mzData");
    let error = store(&path, &original).unwrap_err();
    assert!(
        matches!(&error, Error::Unsupported(message) if message.contains("software metadata")),
        "{error}"
    );
    assert!(!path.exists(), "a refused store must not create the file");
}

/// The 64-bit m/z default round-trips coordinates the source's hardcoded
/// `precision="32"` would destroy.
#[test]
fn default_writer_keeps_full_mz_precision() {
    let mut experiment = MSExperiment::new();
    experiment.spectra.push(MSSpectrum {
        rt: 12.5,
        native_id: "spectrum=1".into(),
        peaks: vec![
            Peak1D::new(1234.5678901234567, 1.0),
            Peak1D::new(2345.678901234568, 2.0),
        ],
        ..Default::default()
    });
    let mut bytes = Vec::new();
    write_with_options(&mut bytes, &experiment, &WriteOptions::default()).unwrap();
    let text = String::from_utf8(bytes).unwrap();
    assert!(text.contains("<data precision=\"64\" endian=\"little\" length=\"2\">"));
    let reloaded = read_text(&text).unwrap();
    assert_eq!(reloaded.spectra[0].peaks[0].mz, 1234.5678901234567);
    assert_eq!(reloaded.spectra[0].peaks[1].mz, 2345.678901234568);

    let mut bytes = Vec::new();
    write_with_options(&mut bytes, &experiment, &WriteOptions::source()).unwrap();
    let text = String::from_utf8(bytes).unwrap();
    assert!(text.contains("<data precision=\"32\" endian=\"little\" length=\"2\">"));
    let reloaded = read_text(&text).unwrap();
    assert_ne!(reloaded.spectra[0].peaks[0].mz, 1234.5678901234567);
    similar(reloaded.spectra[0].peaks[0].mz, 1234.5678901234567);
}

// ===========================================================================
// START_SECTION([EXTRA] storing / loading of meta data arrays)
// ===========================================================================

/// The three spectra `MzDataFile_test.cpp:672-711` builds: five peaks at m/z
/// 1 … 5 with intensities 1 … 5, and one, zero and two annotation arrays
/// `MDA1` = 1.1 … 1.5 and `MDA2` = −2.1 … −2.5, at retention times 500, 600
/// and 700.
fn annotation_experiment() -> MSExperiment {
    let mda1 = vec![1.1f32, 1.2, 1.3, 1.4, 1.5];
    let mda2 = vec![-2.1f32, -2.2, -2.3, -2.4, -2.5];
    let base = MSSpectrum {
        peaks: (1..=5)
            .map(|n| Peak1D::new(f64::from(n), n as f32))
            .collect(),
        ..Default::default()
    };
    let mut experiment = MSExperiment::new();

    let mut first = base.clone();
    first.rt = 500.0;
    first
        .float_data_arrays
        .push(openms::kernel::DataArray::new("MDA1", mda1.clone()));
    experiment.spectra.push(first);

    let mut second = base.clone();
    second.rt = 600.0;
    experiment.spectra.push(second);

    let mut third = base;
    third.rt = 700.0;
    third
        .float_data_arrays
        .push(openms::kernel::DataArray::new("MDA1", mda1));
    third
        .float_data_arrays
        .push(openms::kernel::DataArray::new("MDA2", mda2));
    experiment.spectra.push(third);
    experiment
}

/// `MzDataFile_test.cpp:669-828`, the whole section: store the three spectra,
/// reload them and check the array counts, names and values; reload again under
/// an m/z range of [2.5, 7.0) and check that the arrays are filtered with the
/// peaks; then clear the names, store, reload and check that nameless arrays
/// survive.
#[test]
fn store_and_load_annotation_arrays() {
    let experiment = annotation_experiment();
    let directory = TempDir::new(false).unwrap();
    let path = directory.path().join("arrays.mzData");
    store_with_options(&path, &experiment, &WriteOptions::source()).unwrap();

    let reloaded = load(&path).unwrap();
    assert_eq!(reloaded.spectra.len(), 3);
    assert_eq!(reloaded.spectra[0].float_data_arrays.len(), 1);
    assert_eq!(reloaded.spectra[1].float_data_arrays.len(), 0);
    assert_eq!(reloaded.spectra[2].float_data_arrays.len(), 2);
    assert_eq!(reloaded.spectra[0].float_data_arrays[0].name, "MDA1");
    assert_eq!(reloaded.spectra[2].float_data_arrays[0].name, "MDA1");
    assert_eq!(reloaded.spectra[2].float_data_arrays[1].name, "MDA2");
    let mda1 = [1.1, 1.2, 1.3, 1.4, 1.5];
    let mda2 = [-2.1, -2.2, -2.3, -2.4, -2.5];
    assert_eq!(reloaded.spectra[0].float_data_arrays[0].data.len(), 5);
    for (value, expected) in reloaded.spectra[0].float_data_arrays[0]
        .data
        .iter()
        .zip(mda1)
    {
        similar(f64::from(*value), expected);
    }
    assert_eq!(reloaded.spectra[2].float_data_arrays[0].data.len(), 5);
    for (value, expected) in reloaded.spectra[2].float_data_arrays[0]
        .data
        .iter()
        .zip(mda1)
    {
        similar(f64::from(*value), expected);
    }
    assert_eq!(reloaded.spectra[2].float_data_arrays[1].data.len(), 5);
    for (value, expected) in reloaded.spectra[2].float_data_arrays[1]
        .data
        .iter()
        .zip(mda2)
    {
        similar(f64::from(*value), expected);
    }

    let mut options = PeakFileOptions::default();
    options.set_mz_range(range(2.5, 7.0));
    let filtered = load_with_options(&path, &options, &ReadLimits::default())
        .unwrap()
        .experiment;
    assert_eq!(filtered.spectra.len(), 3);
    assert_eq!(filtered.spectra[0].peaks.len(), 3);
    assert_eq!(filtered.spectra[1].peaks.len(), 3);
    assert_eq!(filtered.spectra[2].peaks.len(), 3);
    assert_eq!(filtered.spectra[0].float_data_arrays.len(), 1);
    assert_eq!(filtered.spectra[1].float_data_arrays.len(), 0);
    assert_eq!(filtered.spectra[2].float_data_arrays.len(), 2);
    assert_eq!(filtered.spectra[0].float_data_arrays[0].name, "MDA1");
    assert_eq!(filtered.spectra[2].float_data_arrays[0].name, "MDA1");
    assert_eq!(filtered.spectra[2].float_data_arrays[1].name, "MDA2");
    for (value, expected) in filtered.spectra[0].float_data_arrays[0]
        .data
        .iter()
        .zip([1.3, 1.4, 1.5])
    {
        similar(f64::from(*value), expected);
    }
    for (value, expected) in filtered.spectra[2].float_data_arrays[0]
        .data
        .iter()
        .zip([1.3, 1.4, 1.5])
    {
        similar(f64::from(*value), expected);
    }
    for (value, expected) in filtered.spectra[2].float_data_arrays[1]
        .data
        .iter()
        .zip([-2.3, -2.4, -2.5])
    {
        similar(f64::from(*value), expected);
    }

    let mut nameless = filtered;
    nameless.spectra[0].float_data_arrays[0].name.clear();
    nameless.spectra[2].float_data_arrays[0].name.clear();
    nameless.spectra[2].float_data_arrays[1].name.clear();
    let second_path = directory.path().join("nameless.mzData");
    store_with_options(&second_path, &nameless, &WriteOptions::source()).unwrap();
    let reloaded = load(&second_path).unwrap();
    assert_eq!(reloaded.spectra.len(), 3);
    assert_eq!(reloaded.spectra[0].peaks.len(), 3);
    assert_eq!(reloaded.spectra[1].peaks.len(), 3);
    assert_eq!(reloaded.spectra[2].peaks.len(), 3);
    assert_eq!(reloaded.spectra[0].float_data_arrays.len(), 1);
    assert_eq!(reloaded.spectra[1].float_data_arrays.len(), 0);
    assert_eq!(reloaded.spectra[2].float_data_arrays.len(), 2);
    assert_eq!(reloaded.spectra[0].float_data_arrays[0].name, "");
    assert_eq!(reloaded.spectra[2].float_data_arrays[0].name, "");
    assert_eq!(reloaded.spectra[2].float_data_arrays[1].name, "");
    for (value, expected) in reloaded.spectra[0].float_data_arrays[0]
        .data
        .iter()
        .zip([1.3, 1.4, 1.5])
    {
        similar(f64::from(*value), expected);
    }
    for (value, expected) in reloaded.spectra[2].float_data_arrays[1]
        .data
        .iter()
        .zip([-2.3, -2.4, -2.5])
    {
        similar(f64::from(*value), expected);
    }

    // `write_supplemental_data = false` drops every annotation array, which
    // the default refuses and the source option allows.
    let options = WriteOptions {
        write_supplemental_data: false,
        ..WriteOptions::default()
    };
    let third_path = directory.path().join("no_supplemental.mzData");
    let error = store_with_options(&third_path, &experiment, &options).unwrap_err();
    assert!(
        matches!(&error, Error::Unsupported(message) if message.contains("write_supplemental_data")),
        "{error}"
    );
    let options = WriteOptions {
        write_supplemental_data: false,
        ..WriteOptions::source()
    };
    store_with_options(&third_path, &experiment, &options).unwrap();
    let bare = load(&third_path).unwrap();
    assert_eq!(bare.spectra.len(), 3);
    assert!(
        bare.spectra
            .iter()
            .all(|spectrum| spectrum.float_data_arrays.is_empty())
    );
}

// ===========================================================================
// START_SECTION([EXTRA] static bool isValid(const std::string& filename))
// ===========================================================================

/// `MzDataFile_test.cpp:830-847`, mapped rather than fully ported: both
/// assertions are `f.isValid(tmp, cerr) == true`, a validation of the stored
/// document against `mzData_1_05.xsd`, and that schema does not ship with this
/// crate. What is reproduced is that both documents — the one stored from an
/// empty experiment and the one stored from the loaded fixture — are
/// well-formed and reload, and that the empty-experiment document matches the
/// exact byte layout `writeTo` produces for it, including the schema-mandated
/// placeholder spectrum.
#[test]
fn stored_documents_are_wellformed() {
    let empty = MSExperiment::new();
    let mut bytes = Vec::new();
    write_with_options(&mut bytes, &empty, &WriteOptions::default()).unwrap();
    let text = String::from_utf8(bytes).unwrap();
    assert_eq!(text, EMPTY_DOCUMENT);
    // Reading it back is asymmetric, exactly as in the source: the
    // schema-mandated placeholder contact, analyzer and spectrum become real
    // records.
    let reloaded = read_text(&text).unwrap();
    assert_eq!(reloaded.spectra.len(), 1);
    assert_eq!(reloaded.spectra[0].peaks.len(), 0);
    assert_eq!(reloaded.settings.contacts.len(), 1);
    assert_eq!(reloaded.settings.instrument.mass_analyzers.len(), 1);

    let loaded = load(data(FIXTURE_1)).unwrap();
    let mut bytes = Vec::new();
    write_with_options(&mut bytes, &loaded, &WriteOptions::source()).unwrap();
    let text = String::from_utf8(bytes).unwrap();
    assert!(text.starts_with("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"));
    assert!(text.ends_with("\t</spectrumList>\n</mzData>\n"));
    assert_eq!(read_text(&text).unwrap().spectra.len(), 3);
}

/// `MzDataHandler::writeTo` for an empty experiment, transcribed from
/// `MzDataHandler.cpp:583-1072`. Tier 4: derived from the source text, not
/// from retained C++ output.
const EMPTY_DOCUMENT: &str = concat!(
    "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n",
    "<mzData version=\"1.05\" accessionNumber=\"\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"http://psidev.sourceforge.net/ms/xml/mzdata/mzdata.xsd\">\n",
    "\t<description>\n",
    "\t\t<admin>\n",
    "\t\t\t<sampleName></sampleName>\n",
    "\t\t\t<contact>\n",
    "\t\t\t\t<name></name>\n",
    "\t\t\t\t<institution></institution>\n",
    "\t\t\t</contact>\n",
    "\t\t</admin>\n",
    "\t\t<instrument>\n",
    "\t\t\t<instrumentName></instrumentName>\n",
    "\t\t\t<source>\n",
    "\t\t\t</source>\n",
    "\t\t\t<analyzerList count=\"1\">\n",
    "\t\t\t\t<analyzer>\n",
    "\t\t\t\t</analyzer>\n",
    "\t\t\t</analyzerList>\n",
    "\t\t\t<detector>\n",
    "\t\t\t</detector>\n",
    "\t\t</instrument>\n",
    "\t\t<dataProcessing>\n",
    "\t\t\t<software>\n",
    "\t\t\t\t<name></name>\n",
    "\t\t\t\t<version></version>\n",
    "\t\t\t</software>\n",
    "\t\t</dataProcessing>\n",
    "\t</description>\n",
    "\t<spectrumList count=\"1\">\n",
    "\t\t<spectrum id=\"1\">\n",
    "\t\t\t<spectrumDesc>\n",
    "\t\t\t\t<spectrumSettings>\n",
    "\t\t\t\t\t<spectrumInstrument msLevel=\"1\"/>\n",
    "\t\t\t\t</spectrumSettings>\n",
    "\t\t\t</spectrumDesc>\n",
    "\t\t\t<mzArrayBinary>\n",
    "\t\t\t\t<data length=\"0\" endian=\"little\" precision=\"32\"></data>\n",
    "\t\t\t</mzArrayBinary>\n",
    "\t\t\t<intenArrayBinary>\n",
    "\t\t\t\t<data length=\"0\" endian=\"little\" precision=\"32\"></data>\n",
    "\t\t\t</intenArrayBinary>\n",
    "\t\t</spectrum>\n",
    "\t</spectrumList>\n",
    "</mzData>\n"
);

// ===========================================================================
// START_SECTION(bool isSemanticallyValid(...))
// ===========================================================================

/// `MzDataFile_test.cpp:849-854` is `NOT_TESTABLE`, with the comment that the
/// feature "is not officially supported - the mapping file was hand-crafted".
/// Neither `mzdata-mapping.xml` nor `psi-mzdata.obo` ships with this crate, so
/// the method reports that rather than pretending to validate.
#[test]
fn semantic_validation_is_unsupported() {
    let file = MzDataFile::new();
    let error = file.is_semantically_valid(data(FIXTURE_1)).unwrap_err();
    assert!(
        matches!(&error, Error::Unsupported(message) if message.contains("mzdata-mapping.xml")),
        "{error}"
    );
}

// ===========================================================================
// Native coverage beyond the class test
// ===========================================================================

/// A per-array `endian="big"` is honoured for the m/z array, the intensity
/// array and an annotation array. Dropping this is the single most likely way
/// to break an mzData port, because every upstream peak array is little
/// endian.
#[test]
fn big_endian_arrays_are_honoured() {
    let loaded = read(
        std::fs::File::open(data(FIXTURE_BIG_ENDIAN))
            .map(std::io::BufReader::new)
            .unwrap(),
    )
    .unwrap();
    assert_eq!(loaded.spectra.len(), 1);
    let spectrum = &loaded.spectra[0];
    assert_eq!(spectrum.native_id, "spectrum=7");
    similar(spectrum.rt, 42.0);
    assert_eq!(spectrum.peaks.len(), 3);
    assert_eq!(spectrum.peaks[0].mz, 110.0);
    assert_eq!(spectrum.peaks[1].mz, 120.0);
    assert_eq!(spectrum.peaks[2].mz, 130.0);
    assert_eq!(spectrum.peaks[0].intensity, 100.0);
    assert_eq!(spectrum.peaks[1].intensity, 200.0);
    assert_eq!(spectrum.peaks[2].intensity, 100.0);
    assert_eq!(spectrum.float_data_arrays.len(), 1);
    assert_eq!(spectrum.float_data_arrays[0].name, "widths");
    assert_eq!(spectrum.float_data_arrays[0].data, vec![1.5f32, 2.5, 3.5]);
}

/// The same bytes read as little endian produce different values, which is
/// what proves the byte order is actually used rather than guessed from the
/// host.
#[test]
fn byte_order_changes_the_decoded_values() {
    let big = std::fs::read_to_string(data(FIXTURE_BIG_ENDIAN)).unwrap();
    let little = big.replace("endian=\"big\"", "endian=\"little\"");
    let decoded = read_text(&little).unwrap();
    assert_ne!(decoded.spectra[0].peaks[0].mz, 110.0);
}

/// `precision` and `endian` classification follows the source exactly: only
/// `"32"` selects 32-bit and only `"big"` selects big endian, and any other
/// spelling falls back with a warning.
#[test]
fn transport_attribute_classification() {
    assert_eq!(Precision::parse("32"), (Precision::Bits32, true));
    assert_eq!(Precision::parse("64"), (Precision::Bits64, true));
    assert_eq!(Precision::parse("128"), (Precision::Bits64, false));
    assert_eq!(Precision::parse(""), (Precision::Bits64, false));
    assert_eq!(Precision::Bits32.width(), 4);
    assert_eq!(Precision::Bits64.width(), 8);
    assert_eq!(Precision::Bits32.name(), "32");
    assert_eq!(Endian::parse("big"), (Endian::Big, true));
    assert_eq!(Endian::parse("little"), (Endian::Little, true));
    assert_eq!(Endian::parse("BIG"), (Endian::Little, false));
    assert_eq!(Endian::Big.name(), "big");
    assert_eq!(Endian::default(), Endian::Little);
    assert_eq!(Precision::default(), Precision::Bits64);

    let text = std::fs::read_to_string(data(FIXTURE_BIG_ENDIAN))
        .unwrap()
        .replace("precision=\"64\"", "precision=\"128\"");
    let loaded = read_with_options(
        text.as_bytes(),
        &PeakFileOptions::default(),
        &ReadLimits::default(),
    )
    .unwrap();
    assert!(
        loaded
            .report
            .warnings
            .iter()
            .any(|warning| warning.contains("precision")),
        "{:?}",
        loaded.report.warnings
    );
}

/// ISO-8859-1 metadata is transcoded rather than refused or mangled: the
/// declared encoding is what the source's own writer emits.
#[test]
fn latin1_metadata_is_transcoded() {
    let e = load(data(FIXTURE_LATIN1)).unwrap();
    assert_eq!(e.settings.sample.name, "Prüfmuster");
    assert_eq!(e.settings.instrument.name, "Meßgerät");
    assert_eq!(e.settings.contacts[0].first_name, "François");
    assert_eq!(e.settings.contacts[0].last_name, "Müller");
    assert_eq!(e.settings.contacts[0].institution, "Universität Tübingen");
    assert_eq!(e.spectra[0].peaks.len(), 2);

    // The characters are all inside ISO-8859-1, so a round trip writes them
    // as single bytes and reads them back unchanged.
    let directory = TempDir::new(false).unwrap();
    let path = directory.path().join("latin1.mzData");
    store_with_options(&path, &e, &WriteOptions::source()).unwrap();
    let bytes = std::fs::read(&path).unwrap();
    assert!(bytes.contains(&0xfc), "ü must be one ISO-8859-1 byte");
    assert!(std::str::from_utf8(&bytes).is_err());
    let reloaded = load(&path).unwrap();
    assert_eq!(reloaded.settings.sample.name, "Prüfmuster");
    assert_eq!(reloaded.settings.instrument.name, "Meßgerät");
}

/// Metadata outside ISO-8859-1 survives too: a UTF-8 document reads, and the
/// writer keeps its declared ISO-8859-1 encoding truthful by emitting numeric
/// character references. `writeTo` would emit raw UTF-8 bytes under that
/// declaration.
#[test]
fn non_latin1_metadata_round_trips_as_character_references() {
    let e = load(data(FIXTURE_UTF8)).unwrap();
    assert_eq!(e.settings.sample.name, "日本語");
    assert_eq!(e.settings.instrument.name, "質量分析計");
    assert_eq!(e.settings.contacts[0].institution, "東京");

    let mut bytes = Vec::new();
    write_with_options(&mut bytes, &e, &WriteOptions::source()).unwrap();
    assert!(bytes.is_ascii(), "the output must stay ISO-8859-1");
    let text = String::from_utf8(bytes).unwrap();
    assert!(
        text.contains("<sampleName>&#x65E5;&#x672C;&#x8A9E;</sampleName>"),
        "{text}"
    );
    let reloaded = read_text(&text).unwrap();
    assert_eq!(reloaded.settings.sample.name, "日本語");
    assert_eq!(reloaded.settings.instrument.name, "質量分析計");
}

/// A path whose file name is not ASCII is handled by the path API, not by
/// byte-slicing a string the port did not construct.
#[test]
fn non_ascii_path_loads() {
    let directory = TempDir::new(false).unwrap();
    let path = directory.path().join("日本語.mzData");
    std::fs::copy(data(FIXTURE_3_MINIMAL), &path).unwrap();
    let e = load(&path).unwrap();
    assert_eq!(e.spectra.len(), 1);
    assert!(
        e.settings
            .document
            .loaded_file_path
            .ends_with("日本語.mzData")
    );
}

/// The XML delimiters the source streams unescaped are escaped here, so a
/// metadata value containing them still produces a document that reads back.
#[test]
fn xml_delimiters_in_metadata_are_escaped() {
    let mut experiment = MSExperiment::new();
    experiment.settings.sample.name = "A & B <C> \"D\" 'E'".into();
    experiment.settings.instrument.name = "]]> & <![CDATA[".into();
    let mut bytes = Vec::new();
    write_with_options(&mut bytes, &experiment, &WriteOptions::source()).unwrap();
    let text = String::from_utf8(bytes).unwrap();
    assert!(text.contains("A &amp; B &lt;C&gt;"), "{text}");
    let reloaded = read_text(&text).unwrap();
    assert_eq!(reloaded.settings.sample.name, "A & B <C> \"D\" 'E'");
    assert_eq!(reloaded.settings.instrument.name, "]]> & <![CDATA[");
}

/// Character data split across several events — here by an XML comment inside
/// the base64 — is concatenated, as `data_to_decode_.back() += …` does.
#[test]
fn split_base64_payload_is_concatenated() {
    let text = spectrum_with(concat!(
        "\t\t\t<mzArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADw<!-- split -->Qg==</data>\n",
        "\t\t\t</mzArrayBinary>\n",
        "\t\t\t<intenArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADIQg==</data>\n",
        "\t\t\t</intenArrayBinary>\n"
    ));
    let e = read_text(&text).unwrap();
    assert_eq!(e.spectra[0].peaks.len(), 1);
    assert_eq!(e.spectra[0].peaks[0].mz, 120.0);
    assert_eq!(e.spectra[0].peaks[0].intensity, 100.0);
}

/// Whitespace inside a payload is removed, because "line breaks inside the
/// base64 data are unfortunately no exception".
#[test]
fn whitespace_inside_a_payload_is_removed() {
    let text = spectrum_with(concat!(
        "\t\t\t<mzArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AA\n  Dw\tQg==\n</data>\n",
        "\t\t\t</mzArrayBinary>\n",
        "\t\t\t<intenArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADIQg==</data>\n",
        "\t\t\t</intenArrayBinary>\n"
    ));
    let e = read_text(&text).unwrap();
    assert_eq!(e.spectra[0].peaks[0].mz, 120.0);
}

/// A spectrum whose m/z and intensity arrays disagree in length is refused.
/// `fillData_` only logs it and then indexes the shorter array at every
/// position of the longer one.
#[test]
fn mismatched_peak_array_lengths_are_refused() {
    let text = spectrum_with(concat!(
        "\t\t\t<mzArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"3\">AADcQgAA8EIAAAJD</data>\n",
        "\t\t\t</mzArrayBinary>\n",
        "\t\t\t<intenArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADIQg==</data>\n",
        "\t\t\t</intenArrayBinary>\n"
    ));
    let error = read_text(&text).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("differs from length of intensity data")),
        "{error}"
    );
}

/// An annotation array shorter than the spectrum is refused for the same
/// reason: `fillData_` indexes it at every peak position without checking.
#[test]
fn short_annotation_array_is_refused() {
    let text = spectrum_with(concat!(
        "\t\t\t<mzArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"3\">AADcQgAA8EIAAAJD</data>\n",
        "\t\t\t</mzArrayBinary>\n",
        "\t\t\t<intenArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"3\">AADIQgAASEMAAMhC</data>\n",
        "\t\t\t</intenArrayBinary>\n",
        "\t\t\t<supDataArrayBinary id=\"1\">\n",
        "\t\t\t\t<arrayName>short</arrayName>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADIQg==</data>\n",
        "\t\t\t</supDataArrayBinary>\n"
    ));
    let error = read_text(&text).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("meta data array")),
        "{error}"
    );
}

/// A spectrum with no binary array at all yields no peaks, where `fillData_`
/// reads `precisions_[0]` past the end of an empty vector.
#[test]
fn spectrum_without_binary_arrays_has_no_peaks() {
    let text = document(concat!(
        "\t<spectrumList count=\"1\">\n",
        "\t\t<spectrum id=\"1\">\n",
        "\t\t\t<spectrumDesc>\n",
        "\t\t\t\t<spectrumSettings>\n",
        "\t\t\t\t\t<spectrumInstrument msLevel=\"1\"/>\n",
        "\t\t\t\t</spectrumSettings>\n",
        "\t\t\t</spectrumDesc>\n",
        "\t\t</spectrum>\n",
        "\t</spectrumList>\n"
    ));
    let e = read_text(&text).unwrap();
    assert_eq!(e.spectra.len(), 1);
    assert_eq!(e.spectra[0].peaks.len(), 0);
}

/// A following spectrum does not inherit the previous one's peak count, which
/// `peak_count_` does because it is a handler member assigned only from
/// `<data>` inside `<mzArrayBinary>`.
#[test]
fn peak_count_does_not_leak_between_spectra() {
    let text = document(concat!(
        "\t<spectrumList count=\"2\">\n",
        "\t\t<spectrum id=\"1\">\n",
        "\t\t\t<spectrumDesc>\n",
        "\t\t\t\t<spectrumSettings>\n",
        "\t\t\t\t\t<spectrumInstrument msLevel=\"1\"/>\n",
        "\t\t\t\t</spectrumSettings>\n",
        "\t\t\t</spectrumDesc>\n",
        "\t\t\t<mzArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"3\">AADcQgAA8EIAAAJD</data>\n",
        "\t\t\t</mzArrayBinary>\n",
        "\t\t\t<intenArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"3\">AADIQgAASEMAAMhC</data>\n",
        "\t\t\t</intenArrayBinary>\n",
        "\t\t</spectrum>\n",
        "\t\t<spectrum id=\"2\">\n",
        "\t\t\t<spectrumDesc>\n",
        "\t\t\t\t<spectrumSettings>\n",
        "\t\t\t\t\t<spectrumInstrument msLevel=\"1\"/>\n",
        "\t\t\t\t</spectrumSettings>\n",
        "\t\t\t</spectrumDesc>\n",
        "\t\t</spectrum>\n",
        "\t</spectrumList>\n"
    ));
    let e = read_text(&text).unwrap();
    assert_eq!(e.spectra.len(), 2);
    assert_eq!(e.spectra[0].peaks.len(), 3);
    assert_eq!(e.spectra[1].peaks.len(), 0);
}

/// A declared `length` that disagrees with the payload is a warning, not an
/// error, and the payload wins — as `fillData_` does after its warning.
#[test]
fn declared_length_is_advisory() {
    let text = spectrum_with(concat!(
        "\t\t\t<mzArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"9\">AADcQgAA8EIAAAJD</data>\n",
        "\t\t\t</mzArrayBinary>\n",
        "\t\t\t<intenArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"9\">AADIQgAASEMAAMhC</data>\n",
        "\t\t\t</intenArrayBinary>\n"
    ));
    let loaded = read_with_options(
        text.as_bytes(),
        &PeakFileOptions::default(),
        &ReadLimits::default(),
    )
    .unwrap();
    assert_eq!(loaded.experiment.spectra[0].peaks.len(), 3);
    assert_eq!(loaded.report.warning_count, 1);
    assert!(
        loaded.report.warnings[0].contains("attribute 'length'"),
        "{:?}",
        loaded.report.warnings
    );
}

/// An attacker-controlled `length` is refused before anything is allocated,
/// where `MSExperiment::reserve` and `spec_.reserve` take theirs on trust.
#[test]
fn hostile_declared_lengths_are_refused_before_allocation() {
    let text = spectrum_with(concat!(
        "\t\t\t<mzArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"4000000000\">AADcQg==</data>\n",
        "\t\t\t</mzArrayBinary>\n",
        "\t\t\t<intenArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADIQg==</data>\n",
        "\t\t\t</intenArrayBinary>\n"
    ));
    let error = read_text(&text).unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(message) if message.contains("element ceiling")),
        "{error}"
    );

    let text = document("\t<spectrumList count=\"4000000000\">\n\t</spectrumList>\n");
    let error = read_text(&text).unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(message) if message.contains("spectrum ceiling")),
        "{error}"
    );

    let text = document("\t<spectrumList count=\"-1\">\n\t</spectrumList>\n");
    let error = read_text(&text).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("whole number")),
        "{error}"
    );
}

/// Every ceiling in [`ReadLimits`] refuses rather than allocating.
#[test]
fn read_ceilings_are_enforced() {
    let text = std::fs::read_to_string(data(FIXTURE_1)).unwrap();
    let options = PeakFileOptions::default();
    let tight = |limits: ReadLimits| read_with_options(text.as_bytes(), &options, &limits);

    let error = tight(ReadLimits {
        max_xml_bytes: 100,
        ..ReadLimits::default()
    })
    .unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(m) if m.contains("byte ceiling")),
        "{error}"
    );

    let error = tight(ReadLimits {
        max_xml_depth: 3,
        ..ReadLimits::default()
    })
    .unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(m) if m.contains("depth ceiling")),
        "{error}"
    );

    let error = tight(ReadLimits {
        max_spectra: 2,
        ..ReadLimits::default()
    })
    .unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(m) if m.contains("spectrum ceiling")),
        "{error}"
    );

    let error = tight(ReadLimits {
        max_arrays_per_spectrum: 2,
        ..ReadLimits::default()
    })
    .unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(m) if m.contains("ceiling")),
        "{error}"
    );

    let error = tight(ReadLimits {
        max_array_bytes: 2,
        ..ReadLimits::default()
    })
    .unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(m) if m.contains("byte ceiling")),
        "{error}"
    );

    let error = tight(ReadLimits {
        max_array_elements: 2,
        ..ReadLimits::default()
    })
    .unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(m) if m.contains("ceiling")),
        "{error}"
    );

    let error = tight(ReadLimits {
        max_total_peaks: 2,
        ..ReadLimits::default()
    })
    .unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(m) if m.contains("peak ceiling")),
        "{error}"
    );

    let error = tight(ReadLimits {
        max_text_bytes: 16,
        ..ReadLimits::default()
    })
    .unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(m) if m.contains("character data")),
        "{error}"
    );

    let error = tight(ReadLimits {
        max_metadata_entries: 3,
        ..ReadLimits::default()
    })
    .unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(m) if m.contains("entry ceiling")),
        "{error}"
    );

    // A warning cap truncates the retained list without losing the count.
    let minimal = std::fs::read_to_string(data(FIXTURE_3_MINIMAL)).unwrap();
    let loaded = read_with_options(
        minimal.as_bytes(),
        &options,
        &ReadLimits {
            max_warnings: 1,
            ..ReadLimits::default()
        },
    )
    .unwrap();
    assert_eq!(loaded.report.warnings.len(), 1);
    assert_eq!(loaded.report.warning_count, 2);
}

/// Structurally hostile documents are refused rather than dereferencing the
/// back of an empty vector, which is what the source does in nine places.
#[test]
fn structural_defects_are_refused() {
    let cases: &[(&str, &str)] = &[
        ("<institution>x</institution>", "outside a <contact>"),
        ("<contactInfo>x</contactInfo>", "outside a <contact>"),
        ("<version>x</version>", "outside a <software>"),
        (
            "<sourceFileX><nameOfFile>x</nameOfFile></sourceFileX>",
            "Unhandled character content",
        ),
    ];
    for (body, _) in cases.iter().take(3) {
        let text = format!(
            concat!(
                "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n",
                "<mzData accessionNumber=\"x\"><description><admin>{}</admin></description></mzData>\n"
            ),
            body
        );
        let error = read_text(&text).unwrap_err();
        assert!(matches!(error, Error::Parse { .. }), "{body}: {error}");
    }

    // A missing required attribute is a parse error, as `attributeAsString_`
    // makes it.
    let text = "<?xml version=\"1.0\"?>\n<mzData version=\"1.05\"></mzData>\n";
    let error = read_text(text).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("accessionNumber")),
        "{error}"
    );

    // A malformed payload length is a parse error, as `Base64::decode` makes
    // it.
    let text = spectrum_with(concat!(
        "\t\t\t<mzArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADwQ</data>\n",
        "\t\t\t</mzArrayBinary>\n",
        "\t\t\t<intenArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADIQg==</data>\n",
        "\t\t\t</intenArrayBinary>\n"
    ));
    let error = read_text(&text).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("multiple of 4")),
        "{error}"
    );

    // An unparsable number is refused rather than silently becoming zero, as
    // `asDouble_` makes it.
    let text = document(concat!(
        "\t<spectrumList count=\"1\">\n",
        "\t\t<spectrum id=\"1\">\n",
        "\t\t\t<spectrumDesc>\n\t\t\t\t<spectrumSettings>\n",
        "\t\t\t\t\t<spectrumInstrument msLevel=\"1\">\n",
        "\t\t\t\t\t\t<cvParam cvLabel=\"psi\" accession=\"PSI:1000039\" name=\"TimeInSeconds\" value=\"NaN\"/>\n",
        "\t\t\t\t\t</spectrumInstrument>\n",
        "\t\t\t\t</spectrumSettings>\n\t\t\t</spectrumDesc>\n",
        "\t\t</spectrum>\n\t</spectrumList>\n"
    ));
    let error = read_text(&text).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("Double conversion error")),
        "{error}"
    );

    // A DTD, an unbalanced tree, a mismatched end tag and a foreign root are
    // all refused.
    for text in [
        "<?xml version=\"1.0\"?>\n<!DOCTYPE mzData>\n<mzData accessionNumber=\"x\"></mzData>\n",
        "<?xml version=\"1.0\"?>\n<mzData accessionNumber=\"x\"><description></mzData>\n",
        "<?xml version=\"1.0\"?>\n<mzXML></mzXML>\n",
        "",
    ] {
        assert!(read_text(text).is_err(), "{text}");
    }

    // An undeclared entity is refused rather than expanded.
    let text = document(
        "\t<spectrumList count=\"0\">\n\t\t<spectrum id=\"&external;\">\n\t\t</spectrum>\n\t</spectrumList>\n",
    );
    assert!(read_text(&text).is_err());
}

/// A CDATA section, a character reference and a predefined entity in character
/// data all resolve. Xerces hands the source handler a CDATA section's content
/// through `characters()` like any other text, so refusing it would be
/// stricter than the source.
#[test]
fn cdata_and_entities_in_character_data_resolve() {
    let text = document(concat!(
        "\t<spectrumList count=\"1\">\n",
        "\t\t<spectrum id=\"1\">\n",
        "\t\t\t<spectrumDesc>\n\t\t\t\t<spectrumSettings>\n",
        "\t\t\t\t\t<spectrumInstrument msLevel=\"1\"/>\n",
        "\t\t\t\t</spectrumSettings>\n\t\t\t</spectrumDesc>\n",
        "\t\t\t<mzArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\"><![CDATA[AADwQg==]]></data>\n",
        "\t\t\t</mzArrayBinary>\n",
        "\t\t\t<intenArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADIQg==</data>\n",
        "\t\t\t</intenArrayBinary>\n",
        "\t\t</spectrum>\n\t</spectrumList>\n"
    ));
    let e = read_text(&text).unwrap();
    assert_eq!(e.spectra[0].peaks[0].mz, 120.0);

    let text = concat!(
        "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n",
        "<mzData version=\"1.05\" accessionNumber=\"x\">\n",
        "\t<description>\n\t\t<admin>\n",
        "\t\t\t<sampleName>A &amp; B &#x65E5; &#66;</sampleName>\n",
        "\t\t</admin>\n\t</description>\n",
        "</mzData>\n"
    );
    let e = read_text(text).unwrap();
    assert_eq!(e.settings.sample.name, "A & B 日 B");
}

/// An unrecognised encoding declaration is refused rather than read as UTF-8.
#[test]
fn unknown_encoding_declaration_is_refused() {
    let text = std::fs::read_to_string(data(FIXTURE_3_MINIMAL))
        .unwrap()
        .replace("ISO-8859-1", "Shift_JIS");
    let error = read_text(&text).unwrap_err();
    assert!(
        matches!(&error, Error::Unsupported(message) if message.contains("encoding")),
        "{error}"
    );
}

/// The unknown-scan-mode fallback acts on the spectrum that carries the
/// parameter. `cvParam_` writes it onto `exp_->getSpectra().back()` — the
/// previous spectrum, or the back of an empty vector for the first one.
#[test]
fn unknown_scan_mode_stays_on_its_own_spectrum() {
    let make = |level: u32| {
        document(&format!(
            concat!(
                "\t<spectrumList count=\"1\">\n",
                "\t\t<spectrum id=\"1\">\n",
                "\t\t\t<spectrumDesc>\n\t\t\t\t<spectrumSettings>\n",
                "\t\t\t\t\t<spectrumInstrument msLevel=\"{}\">\n",
                "\t\t\t\t\t\t<cvParam cvLabel=\"psi\" accession=\"PSI:1000036\" name=\"ScanMode\" value=\"PhotodiodeArrayDetector\"/>\n",
                "\t\t\t\t\t</spectrumInstrument>\n",
                "\t\t\t\t</spectrumSettings>\n\t\t\t</spectrumDesc>\n",
                "\t\t</spectrum>\n\t</spectrumList>\n"
            ),
            level
        ))
    };
    let loaded = read_with_options(
        make(1).as_bytes(),
        &PeakFileOptions::default(),
        &ReadLimits::default(),
    )
    .unwrap();
    assert_eq!(
        loaded.experiment.spectra[0].instrument_settings.scan_mode,
        ScanMode::MassSpectrum
    );
    assert!(loaded.report.warnings[0].contains("Assuming full scan"));

    let loaded = read_with_options(
        make(2).as_bytes(),
        &PeakFileOptions::default(),
        &ReadLimits::default(),
    )
    .unwrap();
    assert_eq!(
        loaded.experiment.spectra[0].instrument_settings.scan_mode,
        ScanMode::MsnSpectrum
    );
    assert!(loaded.report.warnings[0].contains("Assuming MSn scan"));
}

/// The whole source vocabulary is exercised through a store/load round trip of
/// one spectrum carrying every enumerated value mzData has a term for.
#[test]
fn every_enumerated_term_round_trips() {
    let mut experiment = MSExperiment::new();
    experiment.settings.sample.state = SampleState::Emulsion;
    experiment.settings.sample.mass = 1.0;
    let source = openms::metadata::IonSource {
        inlet_type: InletType::MovingWire,
        ionization_method: IonizationMethod::Appi,
        polarity: Polarity::Negative,
        ..Default::default()
    };
    experiment.settings.instrument.ion_sources.push(source);
    let detector = openms::metadata::IonDetector {
        detector_type: DetectorType::MultiCollector,
        acquisition_mode: DetectorAcquisitionMode::TransientRecorder,
        ..Default::default()
    };
    experiment.settings.instrument.ion_detectors.push(detector);
    let analyzer = openms::metadata::MassAnalyzer {
        analyzer_type: AnalyzerType::IonStorage,
        resolution_method: ResolutionMethod::TenPercentValley,
        resolution_type: ResolutionType::Proportional,
        scan_direction: ScanDirection::Down,
        scan_law: ScanLaw::Quadratic,
        reflectron_state: ReflectronState::None,
        ..Default::default()
    };
    experiment.settings.instrument.mass_analyzers.push(analyzer);
    let mut spectrum = MSSpectrum {
        native_id: "spectrum=1".into(),
        rt: 3.0,
        ..Default::default()
    };
    spectrum.instrument_settings.scan_mode = ScanMode::ConstantNeutralLoss;
    spectrum.instrument_settings.polarity = Polarity::Negative;
    let mut precursor = openms::kernel::Precursor::new(400.0, 2);
    precursor.activation_energy = 25.0;
    precursor.activation_methods.insert(ActivationMethod::Psd);
    spectrum.precursors.push(precursor);
    experiment.spectra.push(spectrum);

    let mut bytes = Vec::new();
    write_with_options(&mut bytes, &experiment, &WriteOptions::source()).unwrap();
    let text = String::from_utf8(bytes).unwrap();
    for term in [
        "value=\"Emulsion\"",
        "value=\"MovingWire\"",
        "value=\"APPI\"",
        "value=\"NegativeIonMode\"",
        "value=\"Multi-Collector\"",
        "value=\"TransientRecorder\"",
        "value=\"IonStorage\"",
        "value=\"TenPercentValley\"",
        "value=\"Proportional\"",
        "value=\"Down\"",
        "value=\"Quadratic\"",
        "value=\"None\"",
        "value=\"ConstantNeutralLossScan\"",
        "value=\"PSD\"",
    ] {
        assert!(text.contains(term), "{term} missing from\n{text}");
    }
    let reloaded = read_text(&text).unwrap();
    assert_eq!(reloaded.settings.sample.state, SampleState::Emulsion);
    let instrument = &reloaded.settings.instrument;
    assert_eq!(instrument.ion_sources[0].inlet_type, InletType::MovingWire);
    assert_eq!(
        instrument.ion_sources[0].ionization_method,
        IonizationMethod::Appi
    );
    assert_eq!(instrument.ion_sources[0].polarity, Polarity::Negative);
    assert_eq!(
        instrument.ion_detectors[0].detector_type,
        DetectorType::MultiCollector
    );
    assert_eq!(
        instrument.ion_detectors[0].acquisition_mode,
        DetectorAcquisitionMode::TransientRecorder
    );
    assert_eq!(
        instrument.mass_analyzers[0].analyzer_type,
        AnalyzerType::IonStorage
    );
    assert_eq!(
        instrument.mass_analyzers[0].resolution_method,
        ResolutionMethod::TenPercentValley
    );
    assert_eq!(
        instrument.mass_analyzers[0].reflectron_state,
        ReflectronState::None
    );
    assert_eq!(
        reloaded.spectra[0].instrument_settings.scan_mode,
        ScanMode::ConstantNeutralLoss
    );
    assert_eq!(
        reloaded.spectra[0].precursors[0]
            .activation_methods
            .iter()
            .next(),
        Some(&ActivationMethod::Psd)
    );
}

/// `cvStringToEnum_` resolves an unknown term to index 0, which is the unknown
/// value everywhere except `cv_terms_[18]`, where it is `CID`.
#[test]
fn unknown_vocabulary_terms_fall_back_with_a_warning() {
    let text = document(concat!(
        "\t<spectrumList count=\"1\">\n",
        "\t\t<spectrum id=\"1\">\n",
        "\t\t\t<spectrumDesc>\n\t\t\t\t<spectrumSettings>\n",
        "\t\t\t\t\t<spectrumInstrument msLevel=\"2\"/>\n",
        "\t\t\t\t</spectrumSettings>\n",
        "\t\t\t\t<precursorList count=\"1\">\n",
        "\t\t\t\t\t<precursor msLevel=\"1\" spectrumRef=\"0\">\n",
        "\t\t\t\t\t\t<ionSelection>\n",
        "\t\t\t\t\t\t\t<cvParam cvLabel=\"psi\" accession=\"PSI:1000040\" name=\"MassToChargeRatio\" value=\"1\"/>\n",
        "\t\t\t\t\t\t</ionSelection>\n",
        "\t\t\t\t\t\t<activation>\n",
        "\t\t\t\t\t\t\t<cvParam cvLabel=\"psi\" accession=\"PSI:1000044\" name=\"Method\" value=\"ETD\"/>\n",
        "\t\t\t\t\t\t</activation>\n",
        "\t\t\t\t\t</precursor>\n",
        "\t\t\t\t</precursorList>\n",
        "\t\t\t</spectrumDesc>\n",
        "\t\t</spectrum>\n\t</spectrumList>\n"
    ));
    let loaded = read_with_options(
        text.as_bytes(),
        &PeakFileOptions::default(),
        &ReadLimits::default(),
    )
    .unwrap();
    // "ETD" is not in `cv_terms_[18]` = CID;PSD;PD;SID, so the source reads it
    // as CID — the one table whose index 0 is a real value.
    assert_eq!(
        loaded.experiment.spectra[0].precursors[0]
            .activation_methods
            .iter()
            .next(),
        Some(&ActivationMethod::Cid)
    );
    assert!(
        loaded
            .report
            .warnings
            .iter()
            .any(|warning| warning.contains("activation method")),
        "{:?}",
        loaded.report.warnings
    );
}

/// Multiple precursor charges are reported and reset to zero, and a polarity
/// spelling outside the six the source accepts is a warning.
#[test]
fn precursor_and_polarity_tolerances() {
    let text = document(concat!(
        "\t<spectrumList count=\"1\">\n",
        "\t\t<spectrum id=\"1\">\n",
        "\t\t\t<spectrumDesc>\n\t\t\t\t<spectrumSettings>\n",
        "\t\t\t\t\t<spectrumInstrument msLevel=\"2\">\n",
        "\t\t\t\t\t\t<cvParam cvLabel=\"psi\" accession=\"PSI:1000037\" name=\"Polarity\" value=\"sideways\"/>\n",
        "\t\t\t\t\t\t<cvParam cvLabel=\"psi\" accession=\"PSI:1000038\" name=\"TimeInMinutes\" value=\"2\"/>\n",
        "\t\t\t\t\t</spectrumInstrument>\n",
        "\t\t\t\t</spectrumSettings>\n",
        "\t\t\t\t<precursorList count=\"1\">\n",
        "\t\t\t\t\t<precursor msLevel=\"1\" spectrumRef=\"0\">\n",
        "\t\t\t\t\t\t<ionSelection>\n",
        "\t\t\t\t\t\t\t<cvParam cvLabel=\"psi\" accession=\"PSI:1000041\" name=\"ChargeState\" value=\"2\"/>\n",
        "\t\t\t\t\t\t\t<cvParam cvLabel=\"psi\" accession=\"PSI:1000041\" name=\"ChargeState\" value=\"3\"/>\n",
        "\t\t\t\t\t\t</ionSelection>\n",
        "\t\t\t\t\t</precursor>\n",
        "\t\t\t\t</precursorList>\n",
        "\t\t\t</spectrumDesc>\n",
        "\t\t</spectrum>\n\t</spectrumList>\n"
    ));
    let loaded = read_with_options(
        text.as_bytes(),
        &PeakFileOptions::default(),
        &ReadLimits::default(),
    )
    .unwrap();
    assert_eq!(loaded.experiment.spectra[0].precursors[0].charge, 0);
    // PSI:1000038 is minutes; OpenMS stores seconds.
    similar(loaded.experiment.spectra[0].rt, 120.0);
    assert_eq!(
        loaded.experiment.spectra[0].instrument_settings.polarity,
        Polarity::Unknown
    );
    assert!(
        loaded
            .report
            .warnings
            .iter()
            .any(|warning| warning.contains("Multiple precursor charges")),
        "{:?}",
        loaded.report.warnings
    );
    assert!(
        loaded
            .report
            .warnings
            .iter()
            .any(|warning| warning.contains("Invalid scan polarity")),
        "{:?}",
        loaded.report.warnings
    );
}

/// The precursor-m/z filter drops a whole spectrum, as `skip_spectrum_` does.
#[test]
fn precursor_mz_filter_drops_the_spectrum() {
    let mut file = MzDataFile::new();
    file.options_mut().set_precursor_mz_range(range(10.0, 20.0));
    let loaded = file.load_report(data(FIXTURE_1)).unwrap();
    assert_eq!(loaded.experiment.spectra.len(), 2);
    assert_eq!(loaded.report.spectra_skipped, 1);
    assert!(
        loaded
            .experiment
            .spectra
            .iter()
            .all(|spectrum| spectrum.ms_level == 1)
    );
}

/// Every record mzData has no element for is refused by default and named in
/// the error message.
#[test]
fn the_writer_names_what_it_cannot_store() {
    let base = || {
        let mut experiment = MSExperiment::new();
        experiment.spectra.push(MSSpectrum {
            native_id: "spectrum=1".into(),
            rt: 1.0,
            peaks: vec![Peak1D::new(100.0, 1.0)],
            ..Default::default()
        });
        experiment
    };
    let refuses = |experiment: &MSExperiment, fragment: &str| {
        let mut bytes = Vec::new();
        let error = write_with_options(&mut bytes, experiment, &WriteOptions::default())
            .expect_err(fragment);
        assert!(
            matches!(&error, Error::Unsupported(message) if message.contains(fragment)),
            "expected {fragment}, got {error}"
        );
        assert!(bytes.is_empty(), "a refusal must write nothing");
        // The source option accepts it.
        let mut bytes = Vec::new();
        write_with_options(&mut bytes, experiment, &WriteOptions::source()).unwrap();
        assert!(!bytes.is_empty());
    };

    let mut e = base();
    e.chromatograms
        .push(openms::kernel::MSChromatogram::default());
    refuses(&e, "chromatograms");

    let mut e = base();
    e.settings.source_files = vec![Default::default(), Default::default()];
    refuses(&e, "more than one source file");

    let mut e = base();
    e.settings.source_files = vec![openms::metadata::SourceFile {
        checksum: "abc".into(),
        ..Default::default()
    }];
    refuses(&e, "checksum");

    let mut e = base();
    e.settings.contacts = vec![openms::metadata::ContactPerson {
        email: "a@b".into(),
        ..Default::default()
    }];
    refuses(&e, "contact email");

    let mut e = base();
    e.settings.sample.organism = "yeast".into();
    refuses(&e, "sample organism");

    let mut e = base();
    e.settings.instrument.ion_sources = vec![Default::default(), Default::default()];
    refuses(&e, "more than one ion source");

    let mut e = base();
    e.settings.instrument.ion_detectors = vec![Default::default(), Default::default()];
    refuses(&e, "more than one ion detector");

    let mut e = base();
    e.settings.comment = "note".into();
    refuses(&e, "experiment comment");

    let mut e = base();
    e.settings.metadata.insert("k".into(), MetaValue::from("v"));
    refuses(&e, "experiment-level metadata");

    let mut e = base();
    e.spectra[0].name = "scan".into();
    refuses(&e, "a spectrum name");

    let mut e = base();
    e.spectra[0]
        .integer_data_arrays
        .push(openms::kernel::DataArray::new("i", vec![1i32]));
    refuses(&e, "integer or string data arrays");

    let mut e = base();
    e.spectra[0].instrument_settings.scan_windows = vec![Default::default(), Default::default()];
    refuses(&e, "more than one scan window");

    let mut e = base();
    e.spectra[0].instrument_settings.scan_mode = ScanMode::Absorption;
    refuses(&e, "scan mode");

    let mut e = base();
    e.spectra[0].instrument_settings.zoom_scan = true;
    e.spectra[0].instrument_settings.scan_mode = ScanMode::SelectedIonMonitoring;
    refuses(&e, "zoom scan");

    let mut e = base();
    let mut precursor = openms::kernel::Precursor::new(100.0, 1);
    precursor.activation_methods.insert(ActivationMethod::Cid);
    precursor.activation_methods.insert(ActivationMethod::Hcd);
    e.spectra[0].precursors.push(precursor);
    refuses(&e, "more than one activation method");

    let mut e = base();
    let mut precursor = openms::kernel::Precursor::new(100.0, 1);
    precursor.isolation_window_lower_offset = 0.5;
    e.spectra[0].precursors.push(precursor);
    refuses(&e, "isolation window");

    let mut e = base();
    e.spectra[0]
        .metadata
        .insert("k".into(), MetaValue::from("v"));
    refuses(&e, "spectrum-level metadata");

    let mut e = base();
    e.spectra[0].acquisition_info.acquisitions = vec![openms::metadata::Acquisition {
        identifier: "abc".into(),
        ..Default::default()
    }];
    e.spectra[0].spectrum_type = SpectrumType::Centroid;
    refuses(&e, "acquisition identifier");

    // Native IDs that cannot become an integer `id` are refused, and the
    // source option renumbers with a warning.
    let mut e = base();
    e.spectra[0].native_id = "controllerType=0 scan=1".into();
    refuses(&e, "native IDs");
    let mut bytes = Vec::new();
    let report = write_with_options(&mut bytes, &e, &WriteOptions::source()).unwrap();
    assert!(report.renumbered);
    assert!(report.warnings[0].contains("renumbered"));
    assert!(
        String::from_utf8(bytes)
            .unwrap()
            .contains("<spectrum id=\"1\">")
    );
}

/// A spectrum list whose native IDs are all empty is renumbered silently, as
/// `all_empty` makes it.
#[test]
fn empty_native_ids_renumber_without_a_warning() {
    let mut experiment = MSExperiment::new();
    for _ in 0..3 {
        experiment.spectra.push(MSSpectrum {
            rt: 1.0,
            ..Default::default()
        });
    }
    let mut bytes = Vec::new();
    let report = write_with_options(&mut bytes, &experiment, &WriteOptions::default()).unwrap();
    assert!(report.renumbered);
    assert_eq!(report.warning_count, 0);
    let text = String::from_utf8(bytes).unwrap();
    assert!(text.contains("<spectrum id=\"1\">"));
    assert!(text.contains("<spectrum id=\"2\">"));
    assert!(text.contains("<spectrum id=\"3\">"));
}

/// A refused store leaves the destination file untouched.
#[test]
fn a_refused_store_leaves_the_destination_alone() {
    let directory = TempDir::new(false).unwrap();
    let path = directory.path().join("existing.mzData");
    std::fs::write(&path, b"original").unwrap();
    let mut experiment = MSExperiment::new();
    experiment
        .chromatograms
        .push(openms::kernel::MSChromatogram::default());
    assert!(store(&path, &experiment).is_err());
    assert_eq!(std::fs::read(&path).unwrap(), b"original");
}

/// An `<acqInstrument>` element sets the MS level just as `<spectrumInstrument>`
/// does, which is the one alternative spelling the source's dispatch accepts.
#[test]
fn acq_instrument_is_accepted() {
    let text = document(concat!(
        "\t<spectrumList count=\"1\">\n",
        "\t\t<spectrum id=\"1\">\n",
        "\t\t\t<spectrumDesc>\n\t\t\t\t<spectrumSettings>\n",
        "\t\t\t\t\t<acqInstrument msLevel=\"3\"/>\n",
        "\t\t\t\t</spectrumSettings>\n\t\t\t</spectrumDesc>\n",
        "\t\t</spectrum>\n\t</spectrumList>\n"
    ));
    let e = read_text(&text).unwrap();
    assert_eq!(e.spectra[0].ms_level, 3);
}

/// `<supSourceFile>` children are read and dropped, as the three `//ignored`
/// branches do.
#[test]
fn sup_source_file_children_are_dropped() {
    let text = document(concat!(
        "\t<spectrumList count=\"1\">\n",
        "\t\t<spectrum id=\"1\">\n",
        "\t\t\t<spectrumDesc>\n\t\t\t\t<spectrumSettings>\n",
        "\t\t\t\t\t<spectrumInstrument msLevel=\"1\"/>\n",
        "\t\t\t\t</spectrumSettings>\n\t\t\t</spectrumDesc>\n",
        "\t\t\t<supDesc supDataArrayRef=\"1\">\n",
        "\t\t\t\t<supSourceFile>\n",
        "\t\t\t\t\t<nameOfFile>aux.raw</nameOfFile>\n",
        "\t\t\t\t\t<pathToFile>/tmp/</pathToFile>\n",
        "\t\t\t\t\t<fileType>aux</fileType>\n",
        "\t\t\t\t</supSourceFile>\n",
        "\t\t\t</supDesc>\n",
        "\t\t</spectrum>\n\t</spectrumList>\n"
    ));
    let loaded = read_with_options(
        text.as_bytes(),
        &PeakFileOptions::default(),
        &ReadLimits::default(),
    )
    .unwrap();
    assert!(loaded.experiment.settings.source_files.is_empty());
    assert_eq!(loaded.report.warning_count, 0);
}

/// A document with no `<dataProcessing>` leaves the spectra without a
/// processing record, where the source pushes a null `DataProcessingPtr`.
#[test]
fn a_document_without_data_processing_has_no_null_handle() {
    let text = concat!(
        "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n",
        "<mzData version=\"1.05\" accessionNumber=\"x\">\n",
        "\t<spectrumList count=\"1\">\n",
        "\t\t<spectrum id=\"1\">\n",
        "\t\t\t<spectrumDesc>\n\t\t\t\t<spectrumSettings>\n",
        "\t\t\t\t\t<spectrumInstrument msLevel=\"1\"/>\n",
        "\t\t\t\t</spectrumSettings>\n\t\t\t</spectrumDesc>\n",
        "\t\t</spectrum>\n\t</spectrumList>\n",
        "</mzData>\n"
    );
    let e = read_text(text).unwrap();
    assert_eq!(e.spectra.len(), 1);
    assert!(e.spectra[0].data_processing.is_empty());
}

/// The adapter's own read and write settings are plumbed through.
#[test]
fn adapter_settings_are_used() {
    let mut file = MzDataFile::new();
    file.set_limits(ReadLimits {
        max_spectra: 1,
        ..ReadLimits::default()
    });
    assert_eq!(file.limits().max_spectra, 1);
    assert!(file.load(data(FIXTURE_1)).is_err());

    let mut file = MzDataFile::new();
    file.set_discard_unrepresentable(true);
    assert!(file.discards_unrepresentable());
    assert!(file.write_options().discard_unrepresentable);
    file.options_mut().mz_32_bit = true;
    assert!(file.write_options().mz_32_bit);
    file.options_mut().write_supplemental_data = false;
    assert!(!file.write_options().write_supplemental_data);

    let loaded = file.load(data(FIXTURE_1)).unwrap();
    let directory = TempDir::new(false).unwrap();
    let path = directory.path().join("adapter.mzData");
    let report = file.store_report(&path, &loaded).unwrap();
    assert!(!report.renumbered);
    let mut destination = MSExperiment::new();
    file.load_into(&path, &mut destination).unwrap();
    assert_eq!(destination.spectra.len(), 3);
    assert!(
        destination
            .spectra
            .iter()
            .all(|spectrum| spectrum.float_data_arrays.is_empty())
    );
}

/// The m/z and intensity arrays are identified by their parent element, so a
/// document that writes them in the other order is read correctly.
/// `fillData_` identifies them by position among the `<data>` elements, so it
/// would swap them.
#[test]
fn peak_arrays_are_identified_by_their_parent_element() {
    let text = spectrum_with(concat!(
        "\t\t\t<intenArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADIQg==</data>\n",
        "\t\t\t</intenArrayBinary>\n",
        "\t\t\t<mzArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADwQg==</data>\n",
        "\t\t\t</mzArrayBinary>\n"
    ));
    let e = read_text(&text).unwrap();
    assert_eq!(e.spectra[0].peaks[0].mz, 120.0);
    assert_eq!(e.spectra[0].peaks[0].intensity, 100.0);

    // A second array of either kind is refused rather than shifting the
    // annotation arrays by one slot.
    let text = spectrum_with(concat!(
        "\t\t\t<mzArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADwQg==</data>\n",
        "\t\t\t</mzArrayBinary>\n",
        "\t\t\t<mzArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADwQg==</data>\n",
        "\t\t\t</mzArrayBinary>\n",
        "\t\t\t<intenArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADIQg==</data>\n",
        "\t\t\t</intenArrayBinary>\n"
    ));
    let error = read_text(&text).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("second m/z")),
        "{error}"
    );

    // An m/z array with no intensity array is refused, where `fillData_`
    // reads `precisions_[1]` past the end of its vector.
    let text = spectrum_with(concat!(
        "\t\t\t<mzArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADwQg==</data>\n",
        "\t\t\t</mzArrayBinary>\n"
    ));
    let error = read_text(&text).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("without an intensity array")),
        "{error}"
    );
}

/// A `<supDataArrayBinary>` that does not hold exactly one `<data>` child is
/// refused. The source pairs payloads with arrays by position alone, so an
/// array with none shifts every following payload by one slot and the last one
/// reads past the end of `precisions_`.
#[test]
fn annotation_arrays_must_pair_one_to_one_with_their_payloads() {
    let peaks = concat!(
        "\t\t\t<mzArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADwQg==</data>\n",
        "\t\t\t</mzArrayBinary>\n",
        "\t\t\t<intenArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADIQg==</data>\n",
        "\t\t\t</intenArrayBinary>\n"
    );
    // The second annotation array has no payload of its own.
    let text = spectrum_with(&format!(
        concat!(
            "{}",
            "\t\t\t<supDataArrayBinary id=\"1\">\n",
            "\t\t\t\t<arrayName>a</arrayName>\n",
            "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADIQg==</data>\n",
            "\t\t\t</supDataArrayBinary>\n",
            "\t\t\t<supDataArrayBinary id=\"2\">\n",
            "\t\t\t\t<arrayName>b</arrayName>\n",
            "\t\t\t</supDataArrayBinary>\n"
        ),
        peaks
    ));
    let error = read_text(&text).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("<supDataArrayBinary>")),
        "{error}"
    );

    // The first annotation array has two payloads.
    let text = spectrum_with(&format!(
        concat!(
            "{}",
            "\t\t\t<supDataArrayBinary id=\"1\">\n",
            "\t\t\t\t<arrayName>a</arrayName>\n",
            "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADIQg==</data>\n",
            "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADIQg==</data>\n",
            "\t\t\t</supDataArrayBinary>\n"
        ),
        peaks
    ));
    let error = read_text(&text).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("<supDataArrayBinary>")),
        "{error}"
    );

    // An annotation array with no `<arrayName>` keeps an empty name and does
    // not shift anything, where the source pushes the slot from that child.
    let text = spectrum_with(&format!(
        concat!(
            "{}",
            "\t\t\t<supDataArrayBinary id=\"1\">\n",
            "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADIQg==</data>\n",
            "\t\t\t</supDataArrayBinary>\n",
            "\t\t\t<supDataArrayBinary id=\"2\">\n",
            "\t\t\t\t<arrayName>b</arrayName>\n",
            "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AABIQw==</data>\n",
            "\t\t\t</supDataArrayBinary>\n"
        ),
        peaks
    ));
    let e = read_text(&text).unwrap();
    let arrays = &e.spectra[0].float_data_arrays;
    assert_eq!(arrays.len(), 2);
    assert_eq!(arrays[0].name, "");
    assert_eq!(arrays[0].data, vec![100.0f32]);
    assert_eq!(arrays[1].name, "b");
    assert_eq!(arrays[1].data, vec![200.0f32]);
}

/// A nonfinite decoded coordinate or intensity is refused. The source stores
/// whatever the bytes decode to, which makes every later sort and binary
/// search on that spectrum undefined.
#[test]
fn nonfinite_decoded_values_are_refused() {
    // f64 NaN, little endian, then one finite intensity.
    let nan_mz = "AAAAAAAA+H8=";
    let text = spectrum_with(&format!(
        concat!(
            "\t\t\t<mzArrayBinary>\n",
            "\t\t\t\t<data precision=\"64\" endian=\"little\" length=\"1\">{}</data>\n",
            "\t\t\t</mzArrayBinary>\n",
            "\t\t\t<intenArrayBinary>\n",
            "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADIQg==</data>\n",
            "\t\t\t</intenArrayBinary>\n"
        ),
        nan_mz
    ));
    let error = read_text(&text).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("nonfinite")),
        "{error}"
    );

    // f32 +inf intensity.
    let text = spectrum_with(concat!(
        "\t\t\t<mzArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADwQg==</data>\n",
        "\t\t\t</mzArrayBinary>\n",
        "\t\t\t<intenArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AACAfw==</data>\n",
        "\t\t\t</intenArrayBinary>\n"
    ));
    let error = read_text(&text).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("overflows f32")),
        "{error}"
    );
}

/// A `<supDataArrayBinary>` in a spectrum with no `<data>` element at all is
/// refused rather than silently producing an empty annotation array, because
/// the source would read past the end of its decoded lists.
#[test]
fn annotation_array_without_any_payload_is_refused() {
    let text = spectrum_with(concat!(
        "\t\t\t<supDataArrayBinary id=\"1\">\n",
        "\t\t\t\t<arrayName>a</arrayName>\n",
        "\t\t\t</supDataArrayBinary>\n"
    ));
    let error = read_text(&text).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("no <data> at all")),
        "{error}"
    );
}

/// A repeated attribute is not well-formed XML, and Xerces refuses it, so the
/// source parser never sees one; neither does this reader.
#[test]
fn duplicate_attributes_are_refused() {
    let text = spectrum_with(concat!(
        "\t\t\t<mzArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" precision=\"64\" endian=\"little\" length=\"1\">AADwQg==</data>\n",
        "\t\t\t</mzArrayBinary>\n",
        "\t\t\t<intenArrayBinary>\n",
        "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"1\">AADIQg==</data>\n",
        "\t\t\t</intenArrayBinary>\n"
    ));
    assert!(matches!(read_text(&text), Err(Error::Parse { .. })));
}

// ===========================================================================
// [EXTRA] The store/load round trip is closed
// ===========================================================================

/// An annotation array whose length differs from the spectrum's peak count
/// cannot be stored, because the reader refuses a misaligned
/// `<supDataArrayBinary>`.
///
/// `writeTo` writes such an array anyway and only logs
/// (`MzDataHandler.cpp:1032-1037`), which is what made `store` produce a
/// document this crate's own `load` refuses: `fill_data` compares every
/// annotation array with the m/z array, the place where `fillData_` reads past
/// the end of `decoded_list_[2 + i]` (`MzDataHandler.cpp:562-572`). The empty
/// placeholder array that `DataArray`'s documentation permits and
/// `MSExperiment::validate` accepts is the shortest way to reach it.
#[test]
fn annotation_array_length_must_match_the_peak_count() {
    let with_array = |values: Vec<f32>| {
        let mut experiment = MSExperiment::new();
        experiment.spectra.push(MSSpectrum {
            native_id: "spectrum=1".into(),
            rt: 1.0,
            peaks: vec![Peak1D::new(100.0, 1.0), Peak1D::new(200.0, 2.0)],
            float_data_arrays: vec![openms::kernel::DataArray::new("widths", values)],
            ..Default::default()
        });
        experiment
    };
    // The empty placeholder, one value short, and one value too many.
    for values in [vec![], vec![1.5f32], vec![1.5f32, 2.5, 3.5]] {
        let experiment = with_array(values.clone());
        let mut bytes = Vec::new();
        let error =
            write_with_options(&mut bytes, &experiment, &WriteOptions::default()).unwrap_err();
        assert!(
            matches!(&error, Error::Unsupported(message)
                if message.contains("annotation array at index 0")
                    && message.contains("do not match")),
            "{values:?}: {error}"
        );
        assert!(bytes.is_empty(), "a refusal must write nothing");
    }

    // A matching array is written and reloads unchanged.
    let experiment = with_array(vec![1.5f32, 2.5]);
    let mut bytes = Vec::new();
    write_with_options(&mut bytes, &experiment, &WriteOptions::default()).unwrap();
    let text = String::from_utf8(bytes).unwrap();
    let reloaded = read_text(&text).unwrap();
    assert_eq!(
        reloaded.spectra[0].float_data_arrays[0].data,
        vec![1.5f32, 2.5]
    );
}

/// With `discard_unrepresentable` the misaligned array is dropped from the
/// document instead of being written, so the stored file still reloads — and
/// the arrays that survive keep their names, their descriptions and their
/// values, because `<supDesc>` and `<supDataArrayBinary>` are renumbered
/// together.
#[test]
fn misaligned_annotation_arrays_are_dropped_rather_than_written() {
    let mut experiment = MSExperiment::new();
    let mut short = openms::kernel::DataArray::new("short", vec![1.0f32]);
    short
        .metadata
        .insert("Comment".into(), MetaValue::from("dropped"));
    let mut kept = openms::kernel::DataArray::new("kept", vec![7.5f32, 8.5]);
    kept.metadata
        .insert("Comment".into(), MetaValue::from("survives"));
    experiment.spectra.push(MSSpectrum {
        native_id: "spectrum=1".into(),
        rt: 1.0,
        peaks: vec![Peak1D::new(100.0, 1.0), Peak1D::new(200.0, 2.0)],
        float_data_arrays: vec![short, kept],
        ..Default::default()
    });

    let directory = TempDir::new(false).unwrap();
    let path = directory.path().join("dropped.mzData");
    let mut file = MzDataFile::new();
    file.set_discard_unrepresentable(true);
    let report = file.store_report(&path, &experiment).unwrap();
    assert_eq!(report.warning_count, 1);
    assert!(
        report.warnings[0].contains("Length of meta data array (index:'0' name:'short')")
            && report.warnings[0].contains("not stored"),
        "{:?}",
        report.warnings
    );

    let document = std::fs::read_to_string(&path).unwrap();
    assert!(!document.contains("short"), "{document}");
    assert_eq!(
        document.matches("<supDesc supDataArrayRef=\"1\">").count(),
        1
    );
    assert_eq!(document.matches("<supDataArrayBinary id=\"1\">").count(), 1);
    assert!(!document.contains("id=\"2\""), "{document}");

    let reloaded = load(&path).unwrap();
    assert_eq!(reloaded.spectra.len(), 1);
    let arrays = &reloaded.spectra[0].float_data_arrays;
    assert_eq!(arrays.len(), 1);
    assert_eq!(arrays[0].name, "kept");
    assert_eq!(arrays[0].data, vec![7.5f32, 8.5]);
    assert_eq!(meta(arrays[0].metadata.get("Comment")), "survives");
    assert_eq!(reloaded.spectra[0].peaks.len(), 2);
}

/// A nonfinite coordinate, intensity or annotation value is refused in both
/// modes: the reader refuses such an array, so writing one would produce a
/// document this crate cannot load back. `writeBinary_` encodes whatever the
/// `float` holds.
#[test]
fn nonfinite_values_are_refused_by_the_writer() {
    let cases: &[(&str, MSSpectrum)] = &[
        (
            "nonfinite m/z",
            MSSpectrum {
                native_id: "spectrum=1".into(),
                rt: 1.0,
                peaks: vec![Peak1D::new(f64::NAN, 1.0)],
                ..Default::default()
            },
        ),
        (
            "nonfinite intensity",
            MSSpectrum {
                native_id: "spectrum=1".into(),
                rt: 1.0,
                peaks: vec![Peak1D::new(100.0, f32::INFINITY)],
                ..Default::default()
            },
        ),
        (
            "nonfinite annotation",
            MSSpectrum {
                native_id: "spectrum=1".into(),
                rt: 1.0,
                peaks: vec![Peak1D::new(100.0, 1.0)],
                float_data_arrays: vec![openms::kernel::DataArray::new("widths", vec![f32::NAN])],
                ..Default::default()
            },
        ),
    ];
    for (what, spectrum) in cases {
        let mut experiment = MSExperiment::new();
        experiment.spectra.push(spectrum.clone());
        for options in [WriteOptions::default(), WriteOptions::source()] {
            let mut bytes = Vec::new();
            let error = write_with_options(&mut bytes, &experiment, &options).unwrap_err();
            assert!(
                matches!(&error, Error::InvalidValue(message) if message.contains("nonfinite")),
                "{what}: {error}"
            );
            assert!(bytes.is_empty(), "{what}: a refusal must write nothing");
        }
    }
}

/// A big-endian document survives a store followed by a load. The writer emits
/// little endian only, as `writeBinary_` does, so this is the one direction in
/// which the per-array `endian` attribute has to be honoured on reading for the
/// values to come back at all.
#[test]
fn big_endian_input_survives_a_round_trip() {
    let loaded = load(data(FIXTURE_BIG_ENDIAN)).unwrap();
    let directory = TempDir::new(false).unwrap();
    let path = directory.path().join("from_big_endian.mzData");
    store(&path, &loaded).unwrap();
    let document = std::fs::read_to_string(&path).unwrap();
    assert!(!document.contains("endian=\"big\""), "{document}");
    let reloaded = load(&path).unwrap();
    assert_eq!(reloaded, loaded);
    let spectrum = &reloaded.spectra[0];
    assert_eq!(spectrum.peaks[0].mz, 110.0);
    assert_eq!(spectrum.peaks[1].mz, 120.0);
    assert_eq!(spectrum.peaks[2].mz, 130.0);
    assert_eq!(spectrum.peaks[1].intensity, 200.0);
    assert_eq!(spectrum.float_data_arrays[0].name, "widths");
    assert_eq!(spectrum.float_data_arrays[0].data, vec![1.5f32, 2.5, 3.5]);
}

/// An `&…;` reference inside character data is resolved, never dropped.
///
/// quick-xml reports every reference in character data as its own
/// `Event::GeneralRef`, so a reader that takes character data through a
/// catch-all arm silently deletes it — the defect found in the Mascot XML
/// reader, where `1&#46;5` became `15`. The established fix is the explicit
/// arm at `src/format/mzml.rs:2465`, applied to `src/format/imzml_handler.rs`
/// in commit `c218077`. This reader resolves the five predefined entities and
/// every numeric character reference, refuses an undeclared one, and takes
/// CDATA and a DTD through arms of their own. A base64 payload is where
/// dropping a reference would change numbers with no error at all.
#[test]
fn entity_references_in_character_data_never_vanish() {
    // 110/120/130 as three 32-bit little-endian values, with the payload's
    // last character written as the numeric reference for the same character.
    let payload = "AADcQgAA8EIAAAJD";
    let spliced = format!("{}&#x44;", &payload[..payload.len() - 1]);
    let text = spectrum_with(&format!(
        concat!(
            "\t\t\t<mzArrayBinary>\n",
            "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"3\">{}</data>\n",
            "\t\t\t</mzArrayBinary>\n",
            "\t\t\t<intenArrayBinary>\n",
            "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\"3\">AADIQgAASEMAAMhC</data>\n",
            "\t\t\t</intenArrayBinary>\n"
        ),
        spliced
    ));
    let referenced = read_text(&text).unwrap();
    let literal = read_text(&text.replace("&#x44;", "D")).unwrap();
    assert_eq!(referenced.spectra[0].peaks.len(), 3);
    assert_eq!(referenced.spectra[0].peaks[0].mz, 110.0);
    assert_eq!(referenced.spectra[0].peaks[1].mz, 120.0);
    assert_eq!(referenced.spectra[0].peaks[2].mz, 130.0);
    assert_eq!(referenced.spectra[0].peaks[1].intensity, 200.0);
    assert_eq!(referenced, literal);

    // The Mascot shape: a decimal point written as a reference must not turn
    // 1.5 into 15.
    let text = concat!(
        "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n",
        "<mzData version=\"1.05\" accessionNumber=\"x\">\n",
        "\t<description>\n\t\t<admin>\n",
        "\t\t\t<sampleName>1&#46;5</sampleName>\n",
        "\t\t</admin>\n",
        "\t\t<instrument>\n\t\t\t<instrumentName>a&amp;b</instrumentName>\n\t\t</instrument>\n",
        "\t</description>\n",
        "</mzData>\n"
    );
    let loaded = read_text(text).unwrap();
    assert_eq!(loaded.settings.sample.name, "1.5");
    assert_eq!(loaded.settings.instrument.name, "a&b");

    // An undeclared entity in character data is refused rather than expanded
    // or dropped.
    let text = concat!(
        "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n",
        "<mzData version=\"1.05\" accessionNumber=\"x\">\n",
        "\t<description>\n\t\t<admin>\n",
        "\t\t\t<sampleName>&external;</sampleName>\n",
        "\t\t</admin>\n\t</description>\n",
        "</mzData>\n"
    );
    let error = read_text(text).unwrap_err();
    assert!(
        matches!(&error, Error::Unsupported(message) if message.contains("external")),
        "{error}"
    );
}

/// What a store followed by a load still does *not* preserve, for the record.
///
/// Both losses are the source's, reproduced deliberately, and both are cheap
/// to mistake for the defect above. `writeCVS_` writes nothing for a numeric
/// value that is exactly zero (`MzDataHandler.cpp:1442-1448`), so a retention
/// time of exactly 0 s emits no `TimeInSeconds` and reads back as the
/// `MSSpectrum` default of −1; and `writeTo` emits one zero-length
/// placeholder spectrum for an empty experiment, because the schema requires
/// at least one (`:1057-1072`), so an empty experiment reloads with a
/// spectrum in it.
#[test]
fn the_two_round_trip_losses_that_are_the_sources_own() {
    let mut experiment = MSExperiment::new();
    experiment.spectra.push(MSSpectrum {
        native_id: "spectrum=1".into(),
        rt: 0.0,
        peaks: vec![Peak1D::new(100.0, 1.0)],
        ..Default::default()
    });
    let mut bytes = Vec::new();
    write_with_options(&mut bytes, &experiment, &WriteOptions::default()).unwrap();
    let text = String::from_utf8(bytes).unwrap();
    assert!(!text.contains("TimeInSeconds"), "{text}");
    let reloaded = read_text(&text).unwrap();
    assert_eq!(reloaded.spectra[0].rt, -1.0);
    // Any other retention time returns.
    experiment.spectra[0].rt = 60.0;
    let mut bytes = Vec::new();
    write_with_options(&mut bytes, &experiment, &WriteOptions::default()).unwrap();
    let reloaded = read_text(&String::from_utf8(bytes).unwrap()).unwrap();
    assert_eq!(reloaded.spectra[0].rt, 60.0);

    let empty = MSExperiment::new();
    let mut bytes = Vec::new();
    write_with_options(&mut bytes, &empty, &WriteOptions::default()).unwrap();
    let reloaded = read_text(&String::from_utf8(bytes).unwrap()).unwrap();
    assert_eq!(reloaded.spectra.len(), 1);
    assert!(reloaded.spectra[0].peaks.is_empty());
}
