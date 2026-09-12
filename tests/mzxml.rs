// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Coverage for `FORMAT/MzXMLFile.h` and `FORMAT/HANDLERS/MzXMLHandler.h`.
//!
//! Every `START_SECTION` of `MzXMLFile_test.cpp` (652 lines, 15 sections) is
//! represented here. Literals taken from that file and from the four unmodified
//! upstream fixtures are transcribed source review (tier 3): the four scans of
//! `MzXMLFile_1.mzXML` with their 1/3/5/5 peak counts and `scan=10`..`scan=13`
//! native IDs, the peak coordinates 100..140 with intensities 100..300, the two
//! `<parentFile>` records and their SHA-1 digests, the two `<dataProcessing>`
//! steps, the instrument and contact values, the three precursors of scan 13,
//! the 64-bit fixture's retention times 1/121/3661, the 997530-peak long
//! spectrum and the total ion current 2300 the `transform` sections assert.
//!
//! The resource ceilings, the byte-order/precision/content-type rejections, the
//! writer's lossless defaults, the nested-scan metadata attribution, the
//! retention-time sign and the non-ASCII inputs are independently derived
//! (tier 4), because no upstream fixture reaches them. `MzXMLFile_5_nested.mzXML`
//! is hand-written for this suite. No C++ was built or executed, so nothing
//! here is a tier 1 differential.

#![cfg(feature = "mzml")]

use openms::Error;
use openms::format::mzxml::{
    self, MzXMLFile, PeakPrecision, ReadLimits, ReadOptions, TransformOptions, WriteOptions,
};
use openms::interfaces::MSDataConsumer;
use openms::kernel::{MSChromatogram, MSExperiment, MSSpectrum, NumericRange, Peak1D, Precursor};
use openms::metadata::{
    AnalyzerType, ChecksumType, DetectorType, ExperimentalSettings, InletType, IonizationMethod,
    Polarity, ProcessingAction, ResolutionMethod, ResolutionType, ScanMode,
};
use openms::system::file::TempDir;
use std::ops::ControlFlow;
use std::path::{Path, PathBuf};

/// Upstream `TOLERANCE_ABSOLUTE(0.01)` in the `load` sections.
const TOLERANCE: f64 = 0.01;

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}
fn fixture_1() -> MSExperiment {
    mzxml::load(data("MzXMLFile_1.mzXML")).unwrap()
}
fn close(left: f64, right: f64) {
    assert!(
        (left - right).abs() <= TOLERANCE,
        "{left} differs from {right} by more than {TOLERANCE}"
    );
}
fn range(low: f64, high: f64) -> NumericRange {
    NumericRange {
        min: low,
        max: high,
    }
}
fn text(experiment: &MSExperiment, spectrum: usize, key: &str) -> String {
    experiment.spectra[spectrum].metadata[key].to_string()
}

/// Sums intensities and counts spectra, the upstream `TICConsumer` of
/// `MzXMLFile_test.cpp:29-62`.
#[derive(Default)]
struct TicConsumer {
    tic: f64,
    spectra: usize,
    peaks: usize,
    expected: (usize, usize),
    settings_source_files: usize,
    stop_after: Option<usize>,
}
impl MSDataConsumer for TicConsumer {
    fn set_expected_size(&mut self, spectra: usize, chromatograms: usize) -> Result<(), Error> {
        self.expected = (spectra, chromatograms);
        Ok(())
    }
    fn set_experimental_settings(&mut self, settings: &ExperimentalSettings) -> Result<(), Error> {
        self.settings_source_files = settings.source_files.len();
        Ok(())
    }
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<ControlFlow<()>, Error> {
        for peak in &spectrum.peaks {
            self.tic += f64::from(peak.intensity);
        }
        self.peaks += spectrum.peaks.len();
        self.spectra += 1;
        if self.stop_after == Some(self.spectra) {
            return Ok(ControlFlow::Break(()));
        }
        Ok(ControlFlow::Continue(()))
    }
    fn consume_chromatogram(
        &mut self,
        _chromatogram: &mut MSChromatogram,
    ) -> Result<ControlFlow<()>, Error> {
        Ok(ControlFlow::Continue(()))
    }
}

// ===========================================================================
// START_SECTION((MzXMLFile())) and START_SECTION((~MzXMLFile()))
// ===========================================================================

/// Upstream constructs and deletes the adapter and only checks the pointer is
/// non-null (`MzXMLFile_test.cpp:74-85`). Construction cannot fail here and
/// `Drop` is derived, so the observable part is the registered schema version
/// and the default options the constructor installs.
#[test]
fn default_adapter_registers_schema_3_1() {
    let file = MzXMLFile::new();
    assert_eq!(file.version(), "3.1");
    assert_eq!(mzxml::SCHEMA_VERSION, "3.1");
    assert_eq!(mzxml::SCHEMA, "/SCHEMAS/mzXML_idx_3.1.xsd");
    assert!(!file.options.peaks.has_ms_levels());
    // MzXMLFile is dropped here; the upstream destructor section has no assertion.
}

// ===========================================================================
// START_SECTION(const PeakFileOptions& getOptions() const)
// START_SECTION(PeakFileOptions& getOptions())
// ===========================================================================

/// `TEST_EQUAL(file.getOptions().hasMSLevels(), false)`,
/// `MzXMLFile_test.cpp:87-92`.
#[test]
fn default_options_have_no_ms_levels() {
    let file = MzXMLFile::new();
    assert!(!file.options.peaks.has_ms_levels());
}

/// `file.getOptions().addMSLevel(1)` then
/// `TEST_EQUAL(file.getOptions().hasMSLevels(), true)`,
/// `MzXMLFile_test.cpp:94-100`.
#[test]
fn options_are_mutable_through_the_adapter() {
    let mut file = MzXMLFile::new();
    file.options.peaks.add_ms_level(1).unwrap();
    assert!(file.options.peaks.has_ms_levels());
    assert!(file.options.peaks.contains_ms_level(1));
}

// ===========================================================================
// START_SECTION(void load(const std::string& filename, MapType& map))
// MzXMLFile_test.cpp:102-350 — split by the fixture's own section comments.
// ===========================================================================

/// `TEST_EXCEPTION(Exception::FileNotFound, file.load("dummy/dummy.mzXML", e))`,
/// `MzXMLFile_test.cpp:110`. This port reports a missing file as [`Error::Io`].
#[test]
fn loading_a_missing_file_is_an_io_error() {
    let error = mzxml::load("dummy/dummy.mzXML").unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error:?}");
}

/// `TEST_STRING_EQUAL(e.getLoadedFilePath(), ...)` and
/// `TEST_STRING_EQUAL(FileTypes::typeToName(e.getLoadedFileType()), "mzXML")`,
/// `MzXMLFile_test.cpp:116-117`.
#[test]
fn load_records_the_document_identity() {
    let experiment = fixture_1();
    assert!(
        experiment
            .settings
            .document
            .loaded_file_path
            .ends_with("tests/data/MzXMLFile_1.mzXML"),
        "{}",
        experiment.settings.document.loaded_file_path
    );
    assert_eq!(
        experiment.settings.document.loaded_file_type.name(),
        "mzXML"
    );
}

/// The peak block of `MzXMLFile_test.cpp:119-160`: four scans with MS levels
/// 1/1/1/2, sizes 1/3/5/5, native IDs `scan=10`..`scan=13`, the coordinates and
/// intensities listed in the fixture's comment, and scan 10's two `nameValue`
/// entries plus its `<comment>`.
#[test]
fn load_reads_the_four_scans_of_fixture_1() {
    let experiment = fixture_1();
    assert_eq!(experiment.len(), 4);
    let levels: Vec<u32> = experiment.spectra.iter().map(|s| s.ms_level).collect();
    assert_eq!(levels, vec![1, 1, 1, 2]);
    let sizes: Vec<usize> = experiment.spectra.iter().map(|s| s.peaks.len()).collect();
    assert_eq!(sizes, vec![1, 3, 5, 5]);
    let ids: Vec<&str> = experiment
        .spectra
        .iter()
        .map(|s| s.native_id.as_str())
        .collect();
    assert_eq!(ids, vec!["scan=10", "scan=11", "scan=12", "scan=13"]);

    close(experiment.spectra[0].peaks[0].mz, 120.0);
    close(f64::from(experiment.spectra[0].peaks[0].intensity), 100.0);
    let second: Vec<(f64, f64)> = experiment.spectra[1]
        .peaks
        .iter()
        .map(|p| (p.mz, f64::from(p.intensity)))
        .collect();
    assert_eq!(second, vec![(110.0, 100.0), (120.0, 200.0), (130.0, 100.0)]);
    let third: Vec<(f64, f64)> = experiment.spectra[2]
        .peaks
        .iter()
        .map(|p| (p.mz, f64::from(p.intensity)))
        .collect();
    assert_eq!(
        third,
        vec![
            (100.0, 100.0),
            (110.0, 200.0),
            (120.0, 300.0),
            (130.0, 200.0),
            (140.0, 100.0)
        ]
    );

    assert_eq!(text(&experiment, 0, "URL1"), "www.open-ms.de");
    assert_eq!(text(&experiment, 0, "URL2"), "www.uni-tuebingen.de");
    // Upstream reads this through SpectrumSettings::getComment(); this port has
    // no comment field on MSSpectrum and stores it under the internal key.
    assert_eq!(text(&experiment, 0, mzxml::COMMENT_KEY), "Scan Comment");

    // Retention times are not asserted in this upstream section but are in the
    // RT-range one; the durations PT60S/PT120S/PT180S/PT5S parse to these.
    let rts: Vec<f64> = experiment.spectra.iter().map(|s| s.rt).collect();
    assert_eq!(rts, vec![60.0, 120.0, 180.0, 5.0]);
}

/// The source-file block, `MzXMLFile_test.cpp:162-175`.
#[test]
fn load_reads_both_parent_files() {
    let experiment = fixture_1();
    let sources = &experiment.settings.source_files;
    assert_eq!(sources.len(), 2);
    assert_eq!(sources[0].name, "File_test_1.raw");
    assert_eq!(sources[0].path, "");
    assert_eq!(sources[0].file_type, "RAWData");
    assert_eq!(
        sources[0].checksum,
        "2fd4e1c67a2d28fced849ee1bb76e7391b93eb12"
    );
    assert_eq!(sources[0].checksum_type, ChecksumType::Sha1);
    assert_eq!(sources[1].name, "File_test_2.raw");
    assert_eq!(sources[1].path, "");
    assert_eq!(sources[1].file_type, "processedData");
    assert_eq!(
        sources[1].checksum,
        "2fd4e1c67a2d28fced849ee1bb76e7391b93eb13"
    );
    assert_eq!(sources[1].checksum_type, ChecksumType::Sha1);
}

/// The data-processing block, `MzXMLFile_test.cpp:177-202`: both steps are
/// attached to every spectrum, the first carries `#type` `conversion` with a
/// completion time and no actions, the second carries all three actions and the
/// `#intensity_cutoff` 3.4 with an unset completion time.
#[test]
fn load_attaches_both_data_processing_steps_to_every_scan() {
    let experiment = fixture_1();
    for spectrum in &experiment.spectra {
        let processing = &spectrum.data_processing;
        assert_eq!(processing.len(), 2);

        assert_eq!(processing[0].software.name, "MS-X");
        assert_eq!(processing[0].software.version, "1.0");
        assert_eq!(
            processing[0].metadata[mzxml::PROCESSING_TYPE_KEY].to_string(),
            "conversion"
        );
        assert_eq!(processing[0].metadata["processing 1"].to_string(), "done 1");
        assert_eq!(processing[0].metadata["processing 2"].to_string(), "done 2");
        assert_eq!(
            processing[0].completion_time.unwrap().get(),
            "2001-02-03 04:05:06"
        );
        assert!(processing[0].actions.is_empty());

        assert_eq!(processing[1].software.name, "MS-Y");
        assert_eq!(processing[1].software.version, "1.1");
        assert_eq!(
            processing[1].metadata[mzxml::PROCESSING_TYPE_KEY].to_string(),
            "processing"
        );
        close(
            processing[1].metadata[mzxml::INTENSITY_CUTOFF_KEY]
                .as_f64()
                .unwrap(),
            3.4,
        );
        assert_eq!(processing[1].metadata["processing 3"].to_string(), "done 3");
        // Upstream compares getCompletionTime().get() with the null rendering
        // "0000-00-00 00:00:00"; this port uses None for the unset sentinel.
        assert!(processing[1].completion_time.is_none());
        assert_eq!(processing[1].actions.len(), 3);
        assert!(
            processing[1]
                .actions
                .contains(&ProcessingAction::Deisotoping)
        );
        assert!(
            processing[1]
                .actions
                .contains(&ProcessingAction::ChargeDeconvolution)
        );
        assert!(
            processing[1]
                .actions
                .contains(&ProcessingAction::PeakPicking)
        );
    }
}

/// The instrument block, `MzXMLFile_test.cpp:204-240`.
#[test]
fn load_reads_the_instrument() {
    let experiment = fixture_1();
    let instrument = &experiment.settings.instrument;
    assert_eq!(instrument.vendor, "MS-Vendor");
    assert_eq!(instrument.model, "MS 1");
    assert_eq!(instrument.metadata["URL1"].to_string(), "www.open-ms.de");
    assert_eq!(
        instrument.metadata["URL2"].to_string(),
        "www.uni-tuebingen.de"
    );
    assert_eq!(
        instrument.metadata[mzxml::COMMENT_KEY].to_string(),
        "Instrument Comment"
    );
    assert_eq!(instrument.name, "");
    assert_eq!(instrument.customizations, "");

    assert_eq!(instrument.ion_sources.len(), 1);
    assert_eq!(
        instrument.ion_sources[0].ionization_method,
        IonizationMethod::Esi
    );
    assert_eq!(instrument.ion_sources[0].inlet_type, InletType::Unknown);
    assert_eq!(instrument.ion_sources[0].polarity, Polarity::Unknown);

    assert_eq!(instrument.ion_detectors.len(), 1);
    assert_eq!(
        instrument.ion_detectors[0].detector_type,
        DetectorType::FaradayCup
    );
    assert_eq!(instrument.ion_detectors[0].resolution, 0.0);
    assert_eq!(instrument.ion_detectors[0].adc_sampling_frequency, 0.0);

    assert_eq!(instrument.mass_analyzers.len(), 1);
    let analyzer = &instrument.mass_analyzers[0];
    assert_eq!(analyzer.analyzer_type, AnalyzerType::PaulIonTrap);
    assert_eq!(analyzer.resolution_method, ResolutionMethod::Fwhm);
    assert_eq!(analyzer.resolution_type, ResolutionType::Unknown);
    assert_eq!(analyzer.resolution, 0.0);
    assert_eq!(analyzer.accuracy, 0.0);
    assert_eq!(analyzer.scan_rate, 0.0);
    assert_eq!(analyzer.scan_time, 0.0);
    assert_eq!(analyzer.tof_total_path_length, 0.0);
    assert_eq!(analyzer.isolation_width, 0.0);
    assert_eq!(analyzer.final_ms_exponent, 0);
    assert_eq!(analyzer.magnetic_field_strength, 0.0);

    assert_eq!(instrument.software.name, "MS-Z");
    assert_eq!(instrument.software.version, "3.0");
}

/// The contact and sample blocks, `MzXMLFile_test.cpp:242-262`. mzXML carries
/// no sample description, so every sample field stays at its default.
#[test]
fn load_reads_the_operator_and_leaves_the_sample_empty() {
    let experiment = fixture_1();
    let contacts = &experiment.settings.contacts;
    assert_eq!(contacts.len(), 1);
    assert_eq!(contacts[0].first_name, "FirstName");
    assert_eq!(contacts[0].last_name, "LastName");
    assert_eq!(contacts[0].metadata[mzxml::PHONE_KEY].to_string(), "0049");
    assert_eq!(contacts[0].email, "a@b.de");
    assert_eq!(contacts[0].url, "http://bla.de");
    assert_eq!(contacts[0].contact_info, "");

    let sample = &experiment.settings.sample;
    assert_eq!(sample.name, "");
    assert_eq!(sample.number, "");
    assert_eq!(sample.mass, 0.0);
    assert_eq!(sample.volume, 0.0);
    assert_eq!(sample.concentration, 0.0);
}

/// The precursor block, `MzXMLFile_test.cpp:264-288`: only the nested MS2 scan
/// carries precursors, and each `windowWideness` is halved into two offsets.
#[test]
fn load_reads_the_three_precursors_of_the_nested_scan() {
    let experiment = fixture_1();
    assert_eq!(experiment.spectra[0].precursors.len(), 0);
    assert_eq!(experiment.spectra[1].precursors.len(), 0);
    assert_eq!(experiment.spectra[2].precursors.len(), 0);
    let precursors = &experiment.spectra[3].precursors;
    assert_eq!(precursors.len(), 3);
    for (index, (mz, intensity, offset, charge)) in [
        (101.0, 100.0, 5.0, 1),
        (201.0, 200.0, 10.0, 2),
        (301.0, 300.0, 15.0, 3),
    ]
    .into_iter()
    .enumerate()
    {
        close(precursors[index].mz, mz);
        close(f64::from(precursors[index].intensity), intensity);
        close(precursors[index].isolation_window_lower_offset, offset);
        close(precursors[index].isolation_window_upper_offset, offset);
        assert_eq!(precursors[index].charge, charge);
    }
}

/// `TEST_EQUAL(e == e2, true)` after a second load,
/// `MzXMLFile_test.cpp:290-295`.
#[test]
fn loading_twice_gives_the_same_experiment() {
    assert_eq!(fixture_1(), fixture_1());
}

/// 64-bit peak data, `MzXMLFile_test.cpp:297-329`: three scans, retention times
/// 1/121/3661 from `PT1S`, `PT2M1S` and `PT1H61S`, and the same coordinates the
/// 32-bit fixture carries.
#[test]
fn load_reads_64_bit_peaks_and_composite_durations() {
    let experiment = mzxml::load(data("MzXMLFile_3_64bit.mzXML")).unwrap();
    assert_eq!(experiment.len(), 3);
    assert_eq!(
        experiment
            .spectra
            .iter()
            .map(|s| s.ms_level)
            .collect::<Vec<_>>(),
        vec![1, 1, 1]
    );
    close(experiment.spectra[0].rt, 1.0);
    close(experiment.spectra[1].rt, 121.0);
    close(experiment.spectra[2].rt, 3661.0);
    assert_eq!(
        experiment
            .spectra
            .iter()
            .map(|s| s.peaks.len())
            .collect::<Vec<_>>(),
        vec![1, 3, 5]
    );
    close(experiment.spectra[0].peaks[0].mz, 120.0);
    close(f64::from(experiment.spectra[0].peaks[0].intensity), 100.0);
    // The 64-bit scan; upstream asserts these through TEST_REAL_SIMILAR.
    let second: Vec<(f64, f64)> = experiment.spectra[1]
        .peaks
        .iter()
        .map(|p| (p.mz, f64::from(p.intensity)))
        .collect();
    assert_eq!(second, vec![(110.0, 100.0), (120.0, 200.0), (130.0, 100.0)]);
    let third: Vec<(f64, f64)> = experiment.spectra[2]
        .peaks
        .iter()
        .map(|p| (p.mz, f64::from(p.intensity)))
        .collect();
    assert_eq!(
        third,
        vec![
            (100.0, 100.0),
            (110.0, 200.0),
            (120.0, 300.0),
            (130.0, 200.0),
            (140.0, 100.0)
        ]
    );
    // scanType="zoom" on scan 3 and "Full" on scans 1 and 2, with polarity
    // "any"/"+"/"-" — read but not asserted upstream.
    assert_eq!(
        experiment.spectra[0].instrument_settings.scan_mode,
        ScanMode::MassSpectrum
    );
    assert!(experiment.spectra[2].instrument_settings.zoom_scan);
    assert_eq!(
        experiment.spectra[1].instrument_settings.polarity,
        Polarity::Positive
    );
    assert_eq!(
        experiment.spectra[2].instrument_settings.polarity,
        Polarity::Negative
    );
    let window = &experiment.spectra[1].instrument_settings.scan_windows[0];
    close(window.begin, 110.0);
    close(window.end, 130.0);
}

/// The minimal fixture, `MzXMLFile_test.cpp:331-335`: one scan whose base64 is
/// broken by spaces, a newline and a tab.
#[test]
fn load_strips_whitespace_inside_base64() {
    let experiment = mzxml::load(data("MzXMLFile_2_minimal.mzXML")).unwrap();
    assert_eq!(experiment.len(), 1);
    assert_eq!(experiment.spectra[0].peaks.len(), 1);
    close(experiment.spectra[0].peaks[0].mz, 120.0);
    // No retentionTime attribute: the source defaults to zero, not to the
    // MSSpectrum sentinel -1.
    assert_eq!(experiment.spectra[0].rt, 0.0);
}

/// `TEST_EQUAL(e5[0].size(), 997530)` for `MzXMLFile_4_long.mzXML`,
/// `MzXMLFile_test.cpp:337-341`.
///
/// The 10.6 MB upstream fixture is not copied into this repository (its sha256
/// is recorded in `tests/data/mzxml_provenance.json`); an equivalent document
/// with the same 997530 peaks is generated here, with the base64 wrapped at 76
/// characters so the payload arrives as many text events. The peak count is the
/// upstream literal; the coordinates are independently derived.
#[test]
fn load_reads_a_very_long_spectrum_split_across_text_events() {
    const COUNT: usize = 997_530;
    let mut payload = Vec::with_capacity(COUNT * 8);
    for index in 0..COUNT {
        let mz = 100.0 + index as f32 * 0.001;
        payload.extend_from_slice(&mz.to_be_bytes());
        payload.extend_from_slice(&(1.0f32 + (index % 7) as f32).to_be_bytes());
    }
    let encoded = base64_encode(&payload);
    let mut document = String::with_capacity(encoded.len() + 4096);
    document.push_str(
        "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n<mzXML>\n\t<msRun scanCount=\"1\">\n\
         \t\t<parentFile fileName=\"\" fileType=\"processedData\" \
         fileSha1=\"0000000000000000000000000000000000000000\"/>\n",
    );
    document.push_str(&format!(
        "\t\t<scan num=\"1\" msLevel=\"1\" peaksCount=\"{COUNT}\" polarity=\"any\" \
         retentionTime=\"-PT1S\" startMz=\"97.3951\" endMz=\"1999.93\">\n\
         \t\t\t<peaks precision=\"32\" byteOrder=\"network\" contentType=\"m/z-int\">"
    ));
    for (index, chunk) in encoded.as_bytes().chunks(76).enumerate() {
        if index != 0 {
            document.push('\n');
        }
        document.push_str(std::str::from_utf8(chunk).unwrap());
    }
    document.push_str("</peaks>\n\t\t</scan>\n\t</msRun>\n</mzXML>\n");

    let experiment = mzxml::read(std::io::Cursor::new(document)).unwrap();
    assert_eq!(experiment.len(), 1);
    assert_eq!(experiment.spectra[0].peaks.len(), COUNT);
    close(experiment.spectra[0].peaks[0].mz, 100.0);
    close(
        f64::from(experiment.spectra[0].peaks[COUNT - 1].intensity),
        1.0 + ((COUNT - 1) % 7) as f64,
    );
    // The upstream fixture writes retentionTime="-PT1S" for an unset RT. The
    // source drops the leading sign (it takes the suffix after 'T' first) and
    // reads +1; this port honours it, so a store/load cycle is lossless.
    assert_eq!(experiment.spectra[0].rt, -1.0);
    let window = &experiment.spectra[0].instrument_settings.scan_windows[0];
    close(window.begin, 97.3951);
    close(window.end, 1999.93);
}

/// `TEST_EQUAL(zlib == none, true)`, `MzXMLFile_test.cpp:343-348`.
#[test]
fn zlib_compressed_peaks_equal_the_uncompressed_fixture() {
    let plain = fixture_1();
    let compressed = mzxml::load(data("MzXMLFile_1_compressed.mzXML")).unwrap();
    assert_eq!(plain.spectra, compressed.spectra);
    assert_eq!(plain.settings, compressed.settings);
}

// ===========================================================================
// START_SECTION(([EXTRA] load with metadata only flag))
// ===========================================================================

/// `MzXMLFile_test.cpp:352-373`: no spectra, both parent files, the contact and
/// an empty sample.
#[test]
fn metadata_only_stops_at_the_first_scan() {
    let mut file = MzXMLFile::new();
    file.options.peaks.metadata_only = true;
    let experiment = file.load(data("MzXMLFile_1.mzXML")).unwrap();
    assert_eq!(experiment.len(), 0);
    assert_eq!(experiment.settings.source_files.len(), 2);
    assert_eq!(experiment.settings.source_files[0].name, "File_test_1.raw");
    assert_eq!(experiment.settings.source_files[0].path, "");
    assert_eq!(experiment.settings.contacts.len(), 1);
    assert_eq!(experiment.settings.contacts[0].first_name, "FirstName");
    assert_eq!(experiment.settings.contacts[0].last_name, "LastName");
    assert_eq!(experiment.settings.sample.name, "");
    assert_eq!(experiment.settings.sample.number, "");

    // The same path through the stream API.
    let stream = std::io::BufReader::new(std::fs::File::open(data("MzXMLFile_1.mzXML")).unwrap());
    let settings = mzxml::read_metadata(stream, &ReadOptions::default()).unwrap();
    assert_eq!(settings.source_files.len(), 2);
}

// ===========================================================================
// START_SECTION(([EXTRA] load with selected MS levels))
// ===========================================================================

/// `MzXMLFile_test.cpp:375-403`.
#[test]
fn ms_level_selection_drops_the_nested_ms2_scan() {
    let mut file = MzXMLFile::new();
    file.options.peaks.add_ms_level(1).unwrap();
    let experiment = file.load(data("MzXMLFile_1.mzXML")).unwrap();
    assert_eq!(experiment.len(), 3);
    assert_eq!(
        experiment
            .spectra
            .iter()
            .map(|s| s.ms_level)
            .collect::<Vec<_>>(),
        vec![1, 1, 1]
    );
    assert_eq!(
        experiment
            .spectra
            .iter()
            .map(|s| s.peaks.len())
            .collect::<Vec<_>>(),
        vec![1, 3, 5]
    );
    assert_eq!(
        experiment
            .spectra
            .iter()
            .map(|s| s.native_id.as_str())
            .collect::<Vec<_>>(),
        vec!["scan=10", "scan=11", "scan=12"]
    );

    file.options.peaks.clear_ms_levels();
    let all = file.load(data("MzXMLFile_1.mzXML")).unwrap();
    assert_eq!(all.len(), 4);
}

// ===========================================================================
// START_SECTION(([EXTRA] load with selected MZ range))
// ===========================================================================

/// `MzXMLFile_test.cpp:405-435`: m/z 115..135 keeps 120 and 130 in each scan.
#[test]
fn mz_range_selection_filters_individual_peaks() {
    let mut file = MzXMLFile::new();
    file.options.peaks.set_mz_range(range(115.0, 135.0));
    let experiment = file.load(data("MzXMLFile_1.mzXML")).unwrap();
    assert_eq!(experiment.spectra[0].peaks.len(), 1);
    assert_eq!(experiment.spectra[1].peaks.len(), 2);
    assert_eq!(experiment.spectra[2].peaks.len(), 2);
    close(experiment.spectra[0].peaks[0].mz, 120.0);
    close(f64::from(experiment.spectra[0].peaks[0].intensity), 100.0);
    close(experiment.spectra[1].peaks[0].mz, 120.0);
    close(f64::from(experiment.spectra[1].peaks[0].intensity), 200.0);
    close(experiment.spectra[1].peaks[1].mz, 130.0);
    close(f64::from(experiment.spectra[1].peaks[1].intensity), 100.0);
    close(experiment.spectra[2].peaks[0].mz, 120.0);
    close(f64::from(experiment.spectra[2].peaks[0].intensity), 300.0);
    close(experiment.spectra[2].peaks[1].mz, 130.0);
    close(f64::from(experiment.spectra[2].peaks[1].intensity), 200.0);
}

// ===========================================================================
// START_SECTION(([EXTRA] load with RT range))
// ===========================================================================

/// `MzXMLFile_test.cpp:437-470`: RT 100..200 keeps the scans at 120 s and 180 s.
/// The nested MS2 scan at 5 s is dropped even though its parent is kept.
#[test]
fn rt_range_selection_keeps_two_scans() {
    let mut file = MzXMLFile::new();
    file.options.peaks.set_rt_range(range(100.0, 200.0));
    let experiment = file.load(data("MzXMLFile_1.mzXML")).unwrap();
    assert_eq!(experiment.len(), 2);
    assert_eq!(experiment.spectra[0].peaks.len(), 3);
    assert_eq!(experiment.spectra[1].peaks.len(), 5);
    let first: Vec<(f64, f64)> = experiment.spectra[0]
        .peaks
        .iter()
        .map(|p| (p.mz, f64::from(p.intensity)))
        .collect();
    assert_eq!(first, vec![(110.0, 100.0), (120.0, 200.0), (130.0, 100.0)]);
    let second: Vec<(f64, f64)> = experiment.spectra[1]
        .peaks
        .iter()
        .map(|p| (p.mz, f64::from(p.intensity)))
        .collect();
    assert_eq!(
        second,
        vec![
            (100.0, 100.0),
            (110.0, 200.0),
            (120.0, 300.0),
            (130.0, 200.0),
            (140.0, 100.0)
        ]
    );
}

// ===========================================================================
// START_SECTION(([EXTRA] load with intensity range))
// ===========================================================================

/// `MzXMLFile_test.cpp:472-498`: intensity 150..350 empties the first scan.
#[test]
fn intensity_range_selection_can_empty_a_scan() {
    let mut file = MzXMLFile::new();
    file.options.peaks.set_intensity_range(range(150.0, 350.0));
    let experiment = file.load(data("MzXMLFile_1.mzXML")).unwrap();
    assert_eq!(experiment.spectra[0].peaks.len(), 0);
    assert_eq!(experiment.spectra[1].peaks.len(), 1);
    assert_eq!(experiment.spectra[2].peaks.len(), 3);
    close(experiment.spectra[1].peaks[0].mz, 120.0);
    close(f64::from(experiment.spectra[1].peaks[0].intensity), 200.0);
    close(experiment.spectra[2].peaks[0].mz, 110.0);
    close(f64::from(experiment.spectra[2].peaks[0].intensity), 200.0);
    close(experiment.spectra[2].peaks[1].mz, 120.0);
    close(f64::from(experiment.spectra[2].peaks[1].intensity), 300.0);
    close(experiment.spectra[2].peaks[2].mz, 130.0);
    close(f64::from(experiment.spectra[2].peaks[2].intensity), 200.0);
}

// ===========================================================================
// START_SECTION(([EXTRA] load/store for nested scans))
// ===========================================================================

/// `MzXMLFile_test.cpp:500-568`: six MS-level patterns over five spectra, each
/// stored and reloaded with `TEST_EQUAL(e2.size(), 5)`.
#[test]
fn every_ms_level_pattern_survives_a_store_and_load() {
    let directory = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let file = MzXMLFile::new();
    let mut experiment = MSExperiment::new();
    experiment.spectra = vec![MSSpectrum::default(); 5];
    for (case, levels) in [
        [1, 2, 1, 2, 1],
        [1, 2, 1, 2, 2],
        [1, 2, 3, 2, 3],
        [2, 2, 2, 2, 2],
        [2, 2, 3, 2, 3],
        [2, 1, 2, 3, 1],
    ]
    .into_iter()
    .enumerate()
    {
        for (spectrum, level) in experiment.spectra.iter_mut().zip(levels) {
            spectrum.ms_level = level;
        }
        let path = directory.path().join(format!("nested_{case}.mzXML"));
        file.store(&path, &experiment).unwrap();
        file.load_into(&path, &mut experiment).unwrap();
        assert_eq!(experiment.len(), 5, "case {case} levels {levels:?}");
        assert_eq!(
            experiment
                .spectra
                .iter()
                .map(|s| s.ms_level)
                .collect::<Vec<_>>(),
            levels.to_vec(),
            "case {case}"
        );
    }
}

// ===========================================================================
// START_SECTION(void store(const std::string& filename, const MapType& map))
// ===========================================================================

/// `MzXMLFile_test.cpp:570-584`: store `MzXMLFile_1.mzXML` and reload it with
/// `TEST_TRUE(e1 == e2)`.
#[test]
fn store_round_trips_fixture_1_exactly() {
    let directory = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = directory.path().join("round_trip.mzXML");
    let original = fixture_1();
    assert_eq!(original.len(), 4);
    // The loaded experiment is a valid record, not just a parsed one.
    original.validate().unwrap();
    mzxml::store(&path, &original).unwrap();
    let reloaded = mzxml::load(&path).unwrap();
    assert_eq!(original.spectra, reloaded.spectra);
    assert_eq!(original.settings, reloaded.settings);
    assert_eq!(original, reloaded);
}

/// The same round trip with [`WriteOptions::source`], the lossy source writer.
/// 32-bit peaks and six significant digits still reproduce this fixture because
/// every value in it is exactly representable; the m/z narrowing is shown
/// separately in `source_write_options_narrow_mz_to_f32`.
#[test]
fn store_with_source_options_round_trips_fixture_1() {
    let directory = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = directory.path().join("round_trip_source.mzXML");
    let original = fixture_1();
    mzxml::store_with_options(&path, &original, &WriteOptions::source()).unwrap();
    let reloaded = mzxml::load(&path).unwrap();
    assert_eq!(original, reloaded);
    let document = std::fs::read_to_string(&path).unwrap();
    assert!(document.contains("precision=\"32\""), "{document}");
    assert!(document.contains("retentionTime=\"PT60S\""));
}

// ===========================================================================
// START_SECTION([EXTRA] static bool isValid(const std::string& filename))
// ===========================================================================

/// `MzXMLFile_test.cpp:586-600` validates the stored file against
/// `mzXML_idx_3.1.xsd` through Xerces. That XSD is not embedded in this crate
/// and no mzXML schema validator exists here, so the assertion is replaced by
/// the strongest available structural check: the writer's output parses with
/// every element balanced, declares the 3.1 namespace and schema location, and
/// re-reads to the same experiment.
#[test]
fn the_stored_document_is_structurally_valid() {
    let directory = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = directory.path().join("valid.mzXML");
    let original = fixture_1();
    mzxml::store(&path, &original).unwrap();
    let document = std::fs::read_to_string(&path).unwrap();
    assert!(document.starts_with("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>"));
    assert!(document.contains(mzxml::NAMESPACE), "{document}");
    assert!(document.contains("mzXML_idx_3.1.xsd"));
    assert_eq!(
        document.matches("<scan ").count(),
        document.matches("</scan>").count()
    );
    assert!(document.ends_with("</mzXML>\n"));
    assert_eq!(mzxml::load(&path).unwrap(), original);
}

// ===========================================================================
// START_SECTION(void transform(filename_in, consumer, skip_full_count))
// ===========================================================================

/// `MzXMLFile_test.cpp:602-621`: four spectra, 14 peaks, TIC 2300, nothing
/// stored.
#[test]
fn transform_delivers_every_scan_to_the_consumer() {
    let mut file = MzXMLFile::new();
    file.options.peaks.fill_data = true;
    file.options.peaks.skip_xml_checks = true;
    file.options.peaks.max_data_pool_size = 100;
    file.options.peaks.always_append_data = false;
    let mut consumer = TicConsumer::default();
    let report = file
        .transform(data("MzXMLFile_1.mzXML"), &mut consumer, true)
        .unwrap();
    assert_eq!(consumer.spectra, 4);
    assert_eq!(consumer.peaks, 14);
    close(consumer.tic, 2300.0);
    // skip_full_count=true: the first pass stops at the first scan, so the
    // source reports an expected size of zero (MzXMLFile.cpp:100-107).
    assert_eq!(consumer.expected, (0, 0));
    assert_eq!(consumer.settings_source_files, 2);
    assert_eq!(report.delivered, 4);
    assert!(!report.stopped);
    assert_eq!(report.read.scan_count, 4);
}

/// The counting first pass, which the upstream sections do not exercise:
/// `skip_full_count=false` reports the real scan count.
#[test]
fn transform_without_skip_full_count_reports_four_scans() {
    let file = MzXMLFile::new();
    let mut consumer = TicConsumer::default();
    let report = file
        .transform(data("MzXMLFile_1.mzXML"), &mut consumer, false)
        .unwrap();
    assert_eq!(consumer.expected, (4, 0));
    assert_eq!(report.expected_spectra, 4);
    assert_eq!(consumer.spectra, 4);
}

// ===========================================================================
// START_SECTION(void transform(filename_in, consumer, map, skip_full_count))
// ===========================================================================

/// `MzXMLFile_test.cpp:623-645`: the consumer sees four spectra with TIC 2300
/// and the map keeps all four.
#[test]
fn transform_into_a_map_keeps_every_spectrum() {
    let mut file = MzXMLFile::new();
    file.options.peaks.fill_data = true;
    file.options.peaks.skip_xml_checks = true;
    file.options.peaks.max_data_pool_size = 100;
    file.options.peaks.always_append_data = false;
    let mut consumer = TicConsumer::default();
    let mut map = MSExperiment::new();
    file.transform_into(data("MzXMLFile_1.mzXML"), &mut consumer, &mut map, true)
        .unwrap();
    assert_eq!(consumer.spectra, 4);
    assert_eq!(consumer.peaks, 14);
    close(consumer.tic, 2300.0);
    // The source forces always_append_data for this overload, so the false the
    // test sets above does not suppress the map.
    assert_eq!(map.len(), 4);
    assert_eq!(map.settings.source_files.len(), 2);
}

/// A consumer returning `ControlFlow::Break` ends the read successfully. This
/// port's soft stop; the source interface has no such signal.
#[test]
fn a_consumer_can_stop_the_read() {
    let mut file = MzXMLFile::new();
    file.options.peaks.max_data_pool_size = 1;
    let mut consumer = TicConsumer {
        stop_after: Some(2),
        ..TicConsumer::default()
    };
    let report = file
        .transform(data("MzXMLFile_1.mzXML"), &mut consumer, true)
        .unwrap();
    assert!(report.stopped);
    assert_eq!(report.delivered, 2);
    assert_eq!(consumer.spectra, 2);
}

// ===========================================================================
// Native behaviour: nesting, durations, ceilings, rejections, writer options
// ===========================================================================

/// Nested scans become flat spectra in document order, and a child's metadata
/// stays on the child.
///
/// The source attaches `<nameValue>`, `<comment>` and `<peaks>` to
/// `spectrum_data_.back()`, the most recently *started* scan, and nothing pops
/// that on `</scan>`. A parent's `<nameValue>` written after a nested child
/// therefore lands on the child upstream (`MzXMLHandler.cpp:489`, `:605`). This
/// port keeps a stack of open scans, so each level keeps its own metadata.
#[test]
fn nested_scan_metadata_stays_on_its_own_scan() {
    let experiment = mzxml::load(data("MzXMLFile_5_nested.mzXML")).unwrap();
    assert_eq!(experiment.len(), 4);
    assert_eq!(
        experiment
            .spectra
            .iter()
            .map(|s| s.ms_level)
            .collect::<Vec<_>>(),
        vec![1, 2, 3, 1]
    );
    assert_eq!(
        experiment
            .spectra
            .iter()
            .map(|s| s.native_id.as_str())
            .collect::<Vec<_>>(),
        vec!["scan=1", "scan=2", "scan=3", "scan=4"]
    );
    assert_eq!(text(&experiment, 0, "depth"), "one");
    assert_eq!(text(&experiment, 1, "depth"), "two");
    assert_eq!(text(&experiment, 2, "depth"), "three");
    assert_eq!(
        text(&experiment, 0, mzxml::COMMENT_KEY),
        "level one comment"
    );
    assert_eq!(
        text(&experiment, 1, mzxml::COMMENT_KEY),
        "level two comment"
    );
    assert_eq!(
        text(&experiment, 2, mzxml::COMMENT_KEY),
        "level three comment"
    );
    close(experiment.spectra[0].peaks[0].mz, 120.0);
    close(experiment.spectra[1].peaks[0].mz, 121.0);
    close(experiment.spectra[2].peaks[0].mz, 122.0);
    // The xsi:nil peaks element of scan 4 yields no peaks and no error.
    assert!(experiment.spectra[3].peaks.is_empty());
    let precursor = &experiment.spectra[1].precursors[0];
    close(precursor.mz, 120.0);
    close(precursor.isolation_window_lower_offset, 2.0);
    close(precursor.isolation_window_upper_offset, 2.0);
    assert_eq!(precursor.activation_methods.len(), 1);
}

/// The nested fixture also round-trips, which exercises the writer's
/// level-driven `</scan>` closing on a three-level tree.
#[test]
fn the_nested_fixture_round_trips() {
    let directory = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = directory.path().join("nested_round_trip.mzXML");
    let original = mzxml::load(data("MzXMLFile_5_nested.mzXML")).unwrap();
    mzxml::store(&path, &original).unwrap();
    let document = std::fs::read_to_string(&path).unwrap();
    assert_eq!(
        document.matches("<scan ").count(),
        document.matches("</scan>").count()
    );
    assert_eq!(mzxml::load(&path).unwrap(), original);
}

/// `peaksCount` is the loop bound upstream and its `assert` is compiled out in
/// release builds, so a short payload reads past the decoded buffer
/// (`MzXMLHandler.cpp:1175-1187`). This port compares the two and refuses.
#[test]
fn a_peaks_count_that_disagrees_with_the_payload_is_refused() {
    let document = scan_document(
        "peaksCount=\"1000\"",
        "precision=\"32\" byteOrder=\"network\" contentType=\"m/z-int\"",
        "QvAAAELIAAA=",
    );
    let error = mzxml::read(std::io::Cursor::new(document)).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error:?}");
}

/// An empty `<peaks>` payload with a nonzero `peaksCount` follows the source's
/// early return and yields no peaks rather than an error
/// (`MzXMLHandler.cpp:1153-1156`), with a diagnostic.
#[test]
fn an_empty_payload_with_a_nonzero_count_is_only_reported() {
    let document = scan_document(
        "peaksCount=\"5\"",
        "precision=\"32\" byteOrder=\"network\" contentType=\"m/z-int\"",
        "",
    );
    let mut report = mzxml::ReadReport::default();
    let experiment = mzxml::read_with_report(
        std::io::Cursor::new(document),
        &ReadOptions::default(),
        &mut report,
    )
    .unwrap();
    assert!(experiment.spectra[0].peaks.is_empty());
    assert!(
        report
            .diagnostics
            .iter()
            .any(|line| line.contains("peaksCount 5")),
        "{:?}",
        report.diagnostics
    );
}

/// The four `<peaks>` attributes the source checks with a non-fatal `error`
/// (`MzXMLHandler.cpp:167-197`). The source logs and then decodes anyway, which
/// silently misreads little-endian or non-pair data; this port refuses with
/// [`Error::Unsupported`] and still records the source's message.
#[test]
fn undefined_peaks_attributes_are_refused_with_the_source_message() {
    for (attributes, needle) in [
        (
            "precision=\"16\" byteOrder=\"network\" contentType=\"m/z-int\"",
            "Invalid precision '16'",
        ),
        (
            "precision=\"32\" byteOrder=\"little\" contentType=\"m/z-int\"",
            "byte order 'little'",
        ),
        (
            "precision=\"32\" byteOrder=\"network\" contentType=\"int-m/z\"",
            "pair order 'int-m/z'",
        ),
        (
            "precision=\"32\" byteOrder=\"network\" contentType=\"m/z-int\" compressionType=\"bz2\"",
            "Invalid compression type bz2",
        ),
    ] {
        let document = scan_document("peaksCount=\"1\"", attributes, "QvAAAELIAAA=");
        let mut report = mzxml::ReadReport::default();
        let error = mzxml::read_with_report(
            std::io::Cursor::new(document),
            &ReadOptions::default(),
            &mut report,
        )
        .unwrap_err();
        assert!(matches!(error, Error::Unsupported(_)), "{error:?}");
        assert!(
            report.diagnostics.iter().any(|line| line.contains(needle)),
            "{needle} missing from {:?}",
            report.diagnostics
        );
    }
}

/// A `<peaks>` element with no enclosing `<scan>` indexes
/// `spectrum_data_.back()` on an empty vector upstream, which is undefined
/// behaviour (`MzXMLHandler.cpp:170`). This port reports it.
#[test]
fn peaks_outside_a_scan_is_a_parse_error() {
    let document = "<mzXML><msRun><peaks precision=\"32\" byteOrder=\"network\" \
                    contentType=\"m/z-int\">QvAAAELIAAA=</peaks></msRun></mzXML>";
    let error = mzxml::read(std::io::Cursor::new(document)).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error:?}");
    let document = "<mzXML><msRun><precursorMz precursorIntensity=\"1\">100</precursorMz>\
                    </msRun></mzXML>";
    let error = mzxml::read(std::io::Cursor::new(document)).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error:?}");
}

/// `msLevel="0"` warns and is treated as level 1 (`MzXMLHandler.cpp:248-252`).
#[test]
fn ms_level_zero_becomes_one_with_a_warning() {
    let document = "<mzXML><msRun><scan num=\"1\" msLevel=\"0\" peaksCount=\"0\">\
                    <peaks precision=\"32\" byteOrder=\"network\" contentType=\"m/z-int\"/>\
                    </scan></msRun></mzXML>";
    let mut report = mzxml::ReadReport::default();
    let experiment = mzxml::read_with_report(
        std::io::Cursor::new(document),
        &ReadOptions::default(),
        &mut report,
    )
    .unwrap();
    assert_eq!(experiment.spectra[0].ms_level, 1);
    assert!(
        report
            .diagnostics
            .iter()
            .any(|line| line.contains("Assuming ms level 1")),
        "{:?}",
        report.diagnostics
    );
}

/// Every `scanType` the source maps, including the three non-standard ABI
/// Sashimi values (`MzXMLHandler.cpp:328-388`). `EPI` also rewrites the MS
/// level to 2.
#[test]
fn scan_types_map_to_the_source_scan_modes() {
    for (scan_type, mode, zoom, level) in [
        ("", ScanMode::Unknown, false, 1),
        ("zoom", ScanMode::MassSpectrum, true, 1),
        ("Full", ScanMode::MassSpectrum, false, 1),
        ("SIM", ScanMode::SelectedIonMonitoring, false, 1),
        ("SRM", ScanMode::SelectedReactionMonitoring, false, 1),
        ("MRM", ScanMode::SelectedReactionMonitoring, false, 1),
        ("CRM", ScanMode::ConsecutiveReactionMonitoring, false, 1),
        ("Q1", ScanMode::MassSpectrum, false, 1),
        ("Q3", ScanMode::MassSpectrum, false, 1),
        ("EMS", ScanMode::MassSpectrum, false, 1),
        ("EPI", ScanMode::MassSpectrum, false, 2),
        ("ER", ScanMode::MassSpectrum, true, 1),
        ("nonsense", ScanMode::MassSpectrum, false, 1),
    ] {
        let attribute = if scan_type.is_empty() {
            String::new()
        } else {
            format!(" scanType=\"{scan_type}\"")
        };
        let document = format!(
            "<mzXML><msRun><scan num=\"1\" msLevel=\"1\" peaksCount=\"0\"{attribute}>\
             <peaks precision=\"32\" byteOrder=\"network\" contentType=\"m/z-int\"/>\
             </scan></msRun></mzXML>"
        );
        let experiment = mzxml::read(std::io::Cursor::new(document)).unwrap();
        let settings = &experiment.spectra[0].instrument_settings;
        assert_eq!(settings.scan_mode, mode, "scanType {scan_type}");
        assert_eq!(settings.zoom_scan, zoom, "scanType {scan_type}");
        assert_eq!(
            experiment.spectra[0].ms_level, level,
            "scanType {scan_type}"
        );
    }
    // "Full" at MS level 2 selects MSnSpectrum instead.
    let document = "<mzXML><msRun><scan num=\"1\" msLevel=\"2\" peaksCount=\"0\" \
                    scanType=\"Full\"><peaks precision=\"32\" byteOrder=\"network\" \
                    contentType=\"m/z-int\"/></scan></msRun></mzXML>";
    let experiment = mzxml::read(std::io::Cursor::new(document)).unwrap();
    assert_eq!(
        experiment.spectra[0].instrument_settings.scan_mode,
        ScanMode::MsnSpectrum
    );
}

/// The precursor m/z filter drops the scan whose `<precursorMz>` text lies
/// outside the range (`MzXMLHandler.cpp:581-587`).
#[test]
fn the_precursor_mz_filter_drops_the_matching_scan() {
    let mut options = ReadOptions::default();
    options.peaks.set_precursor_mz_range(range(200.0, 250.0));
    let experiment = mzxml::read_with_options(open(data("MzXMLFile_1.mzXML")), &options).unwrap();
    // Only scan 13 has precursors; its first precursor is at 101, outside the
    // range, so the whole scan goes.
    assert_eq!(experiment.len(), 3);
    assert_eq!(
        experiment
            .spectra
            .iter()
            .map(|s| s.native_id.as_str())
            .collect::<Vec<_>>(),
        vec!["scan=10", "scan=11", "scan=12"]
    );
}

/// `fill_data = false` keeps the scan records and discards their peaks.
#[test]
fn fill_data_false_keeps_scans_without_peaks() {
    let mut options = ReadOptions::default();
    options.peaks.fill_data = false;
    let experiment = mzxml::read_with_options(open(data("MzXMLFile_1.mzXML")), &options).unwrap();
    assert_eq!(experiment.len(), 4);
    assert!(experiment.spectra.iter().all(|s| s.peaks.is_empty()));
}

/// `sort_spectra_by_mz` (the source default) sorts a descending payload, and
/// clearing it preserves the stored order.
#[test]
fn sort_spectra_by_mz_orders_a_descending_payload() {
    let payload = {
        let mut bytes = Vec::new();
        for (mz, intensity) in [(130.0f32, 1.0f32), (110.0, 2.0), (120.0, 3.0)] {
            bytes.extend_from_slice(&mz.to_be_bytes());
            bytes.extend_from_slice(&intensity.to_be_bytes());
        }
        base64_encode(&bytes)
    };
    let document = scan_document(
        "peaksCount=\"3\"",
        "precision=\"32\" byteOrder=\"network\" contentType=\"m/z-int\"",
        &payload,
    );
    let sorted = mzxml::read(std::io::Cursor::new(document.clone())).unwrap();
    assert_eq!(
        sorted.spectra[0]
            .peaks
            .iter()
            .map(|p| p.mz)
            .collect::<Vec<_>>(),
        vec![110.0, 120.0, 130.0]
    );
    let mut options = ReadOptions::default();
    options.peaks.sort_spectra_by_mz = false;
    let stored = mzxml::read_with_options(std::io::Cursor::new(document), &options).unwrap();
    assert_eq!(
        stored.spectra[0]
            .peaks
            .iter()
            .map(|p| p.mz)
            .collect::<Vec<_>>(),
        vec![130.0, 110.0, 120.0]
    );
}

/// `read_scan_count` is the source `LD_RAWCOUNTS` pass: it counts every scan,
/// filters included, and stores none.
#[test]
fn raw_counts_counts_scans_without_storing_them() {
    let (count, settings) =
        mzxml::read_scan_count(open(data("MzXMLFile_1.mzXML")), &ReadOptions::default()).unwrap();
    assert_eq!(count, 4);
    assert_eq!(settings.source_files.len(), 2);
    let (nested, _) = mzxml::read_scan_count(
        open(data("MzXMLFile_5_nested.mzXML")),
        &ReadOptions::default(),
    )
    .unwrap();
    assert_eq!(nested, 4);
}

/// Each ceiling refuses before allocating from the declared length.
#[test]
fn resource_ceilings_refuse_hostile_declarations() {
    let base = ReadLimits::default();
    let checks: Vec<(ReadLimits, String)> = vec![
        (
            ReadLimits {
                max_scans: 2,
                ..base
            },
            "scan".into(),
        ),
        (
            ReadLimits {
                max_peaks_per_scan: 2,
                ..base
            },
            "peaksCount".into(),
        ),
        (
            ReadLimits {
                max_total_peaks: 4,
                ..base
            },
            "total peak".into(),
        ),
        (
            ReadLimits {
                max_decoded_bytes: 8,
                ..base
            },
            "decoded array byte".into(),
        ),
        (
            ReadLimits {
                max_encoded_bytes: 4,
                ..base
            },
            "base64 payload".into(),
        ),
        (
            ReadLimits {
                max_xml_bytes: 64,
                ..base
            },
            "XML byte".into(),
        ),
        (
            ReadLimits {
                max_depth: 3,
                ..base
            },
            "XML depth".into(),
        ),
        (
            ReadLimits {
                max_source_files: 1,
                ..base
            },
            "parentFile".into(),
        ),
        (
            ReadLimits {
                max_data_processing: 1,
                ..base
            },
            "dataProcessing".into(),
        ),
        (
            ReadLimits {
                max_metadata_entries: 1,
                ..base
            },
            "metadata entry".into(),
        ),
    ];
    for (limits, needle) in checks {
        let options = ReadOptions {
            peaks: Default::default(),
            limits,
        };
        let error =
            mzxml::read_with_options(open(data("MzXMLFile_1.mzXML")), &options).unwrap_err();
        match &error {
            Error::InvalidValue(message) => assert!(
                message.contains(&needle),
                "expected {needle} in {message:?}"
            ),
            other => panic!("expected a limit error for {needle}, got {other:?}"),
        }
    }
    // A hostile msRun scanCount is rejected before it is used as a hint.
    let document = "<mzXML><msRun scanCount=\"99999999999\"></msRun></mzXML>";
    let error = mzxml::read(std::io::Cursor::new(document)).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    // The precursor ceiling needs a document with several precursors.
    let options = ReadOptions {
        peaks: Default::default(),
        limits: ReadLimits {
            max_precursors_per_scan: 2,
            ..base
        },
    };
    let error = mzxml::read_with_options(open(data("MzXMLFile_1.mzXML")), &options).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
}

/// Non-ASCII input must never panic. A UTF-8 document keeps its multi-byte
/// metadata, the path itself may be non-ASCII, an ISO-8859-1 declaration with a
/// non-ASCII byte is refused, and a multi-byte `completionTime` is clipped by
/// characters rather than bytes.
#[test]
fn non_ascii_input_is_handled_without_slicing_bytes() {
    let directory = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = directory.path().join("日本語.mzXML");
    let document = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
         <mzXML><msRun>\
         <parentFile fileName=\"データ.raw\" fileType=\"RAWData\" fileSha1=\"\"/>\
         <dataProcessing><software type=\"変換\" name=\"日本語\" version=\"1.0\" \
         completionTime=\"日本語日本語日本語日本語\"/></dataProcessing>\
         <scan num=\"1\" msLevel=\"1\" peaksCount=\"1\" filterLine=\"フィルタ\">\
         <peaks precision=\"32\" byteOrder=\"network\" contentType=\"m/z-int\">\
         QvAAAELIAAA=</peaks>\
         <nameValue name=\"日本語\" value=\"値\"/>\
         <comment>日本語のコメント</comment>\
         </scan></msRun></mzXML>";
    std::fs::write(&path, document).unwrap();
    let mut report = mzxml::ReadReport::default();
    let experiment = mzxml::read_with_report(
        std::io::Cursor::new(document),
        &ReadOptions::default(),
        &mut report,
    )
    .unwrap();
    assert_eq!(experiment.len(), 1);
    assert_eq!(text(&experiment, 0, "日本語"), "値");
    assert_eq!(text(&experiment, 0, mzxml::COMMENT_KEY), "日本語のコメント");
    assert_eq!(text(&experiment, 0, mzxml::FILTER_STRING_KEY), "フィルタ");
    assert_eq!(experiment.settings.source_files[0].name, "データ.raw");
    // The unparsable completionTime is reported, not fatal, and the 19-character
    // clip lands on a character boundary.
    assert!(
        report
            .diagnostics
            .iter()
            .any(|line| line.contains("DateTime conversion error")),
        "{:?}",
        report.diagnostics
    );
    assert!(
        experiment.spectra[0].data_processing[0]
            .completion_time
            .is_none()
    );

    // The same document read through a non-ASCII path.
    let loaded = mzxml::load(&path).unwrap();
    assert_eq!(loaded.spectra, experiment.spectra);

    // An ISO-8859-1 declaration with a byte outside ASCII is refused.
    let latin1 = b"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n<mzXML><msRun>\
                   <parentFile fileName=\"\xfc\" fileType=\"\" fileSha1=\"\"/>\
                   </msRun></mzXML>";
    let error = mzxml::read(std::io::Cursor::new(latin1.to_vec())).unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error:?}");

    // And the round trip through the writer keeps the multi-byte text. The
    // source writer hard-codes encoding="ISO-8859-1" whatever the bytes are;
    // this port declares UTF-8 when the document is not ASCII.
    let out = directory.path().join("日本語_out.mzXML");
    mzxml::store(&out, &experiment).unwrap();
    let written = std::fs::read_to_string(&out).unwrap();
    assert!(written.starts_with("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"));
    assert_eq!(mzxml::load(&out).unwrap().spectra, experiment.spectra);
    // An ASCII experiment keeps the source declaration.
    let ascii = directory.path().join("ascii_out.mzXML");
    mzxml::store(&ascii, &fixture_1()).unwrap();
    assert!(
        std::fs::read_to_string(&ascii)
            .unwrap()
            .starts_with("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>")
    );
}

/// Malformed and hostile documents produce errors, never panics.
#[test]
fn malformed_documents_are_rejected() {
    for document in [
        "",
        "<notMzXML/>",
        "<mzXML><msRun><scan num=\"1\" peaksCount=\"1\"/></msRun></mzXML>",
        "<mzXML><msRun><scan num=\"1\" msLevel=\"x\" peaksCount=\"1\"/></msRun></mzXML>",
        "<mzXML><msRun><scan num=\"1\" msLevel=\"1\" peaksCount=\"-1\"/></msRun></mzXML>",
        "<mzXML><msRun><scan num=\"1\" msLevel=\"1\" peaksCount=\"1\"><peaks \
         precision=\"32\" byteOrder=\"network\" contentType=\"m/z-int\">!!!!</peaks>\
         </scan></msRun></mzXML>",
        "<mzXML><msRun><scan msLevel=\"1\" peaksCount=\"0\"><peaks precision=\"32\" \
         byteOrder=\"network\" contentType=\"m/z-int\"/></scan></msRun></mzXML>",
        "<mzXML><msRun><scan num=\"1\" msLevel=\"1\" peaksCount=\"0\" startMz=\"200\" \
         endMz=\"100\"><peaks precision=\"32\" byteOrder=\"network\" \
         contentType=\"m/z-int\"/></scan></msRun></mzXML>",
        "<!DOCTYPE mzXML><mzXML/>",
        "<mzXML><msRun><parentFile/></msRun></mzXML>",
        "<mzXML/><mzXML/>",
    ] {
        let result = mzxml::read(std::io::Cursor::new(document));
        assert!(result.is_err(), "expected an error for {document:?}");
    }
    // A truncated zlib payload and a nonfinite peak value.
    let document = scan_document(
        "peaksCount=\"1\"",
        "precision=\"32\" byteOrder=\"network\" contentType=\"m/z-int\" \
         compressionType=\"zlib\"",
        "eJxz",
    );
    assert!(mzxml::read(std::io::Cursor::new(document)).is_err());
    let nan = base64_encode(&[0x7f, 0xc0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]);
    let document = scan_document(
        "peaksCount=\"1\"",
        "precision=\"32\" byteOrder=\"network\" contentType=\"m/z-int\"",
        &nan,
    );
    assert!(mzxml::read(std::io::Cursor::new(document)).is_err());
    // A 64-bit intensity beyond the f32 the kernel stores. The source narrows
    // it with an implicit conversion and keeps the resulting infinity.
    let mut payload = 100.0f64.to_be_bytes().to_vec();
    payload.extend_from_slice(&1.0e300f64.to_be_bytes());
    let document = scan_document(
        "peaksCount=\"1\"",
        "precision=\"64\" byteOrder=\"network\" contentType=\"m/z-int\"",
        &base64_encode(&payload),
    );
    let error = mzxml::read(std::io::Cursor::new(document)).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error:?}");
}

/// Duplicate `<peaks>` elements inside one scan are refused rather than
/// concatenated into a meaningless payload.
#[test]
fn two_peaks_elements_in_one_scan_are_refused() {
    let document = "<mzXML><msRun><scan num=\"1\" msLevel=\"1\" peaksCount=\"1\">\
                    <peaks precision=\"32\" byteOrder=\"network\" contentType=\"m/z-int\">\
                    QvAAAELIAAA=</peaks>\
                    <peaks precision=\"32\" byteOrder=\"network\" contentType=\"m/z-int\">\
                    QvAAAELIAAA=</peaks></scan></msRun></mzXML>";
    let error = mzxml::read(std::io::Cursor::new(document)).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error:?}");
}

/// The writer's default 64-bit precision keeps an m/z that `f32` cannot hold,
/// and [`WriteOptions::source`] narrows it, as the source writer always does.
#[test]
fn source_write_options_narrow_mz_to_f32() {
    /// An m/z with more digits than `f32` can hold.
    const WIDE_MZ: f64 = 1234.5678901234;
    let mut experiment = MSExperiment::new();
    let mut spectrum = MSSpectrum {
        ms_level: 1,
        rt: 1.5,
        native_id: "scan=1".into(),
        ..MSSpectrum::default()
    };
    spectrum.peaks.push(Peak1D {
        mz: WIDE_MZ,
        intensity: 7.5,
    });
    experiment.spectra.push(spectrum);

    let directory = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let lossless = directory.path().join("lossless.mzXML");
    mzxml::store(&lossless, &experiment).unwrap();
    let reloaded = mzxml::load(&lossless).unwrap();
    assert_eq!(reloaded.spectra[0].peaks[0].mz, WIDE_MZ);
    assert!(
        std::fs::read_to_string(&lossless)
            .unwrap()
            .contains("precision=\"64\"")
    );

    let lossy = directory.path().join("lossy.mzXML");
    mzxml::store_with_options(&lossy, &experiment, &WriteOptions::source()).unwrap();
    let narrowed = mzxml::load(&lossy).unwrap();
    assert_eq!(narrowed.spectra[0].peaks[0].mz, f64::from(WIDE_MZ as f32));
    assert_ne!(narrowed.spectra[0].peaks[0].mz, WIDE_MZ);
    // Six significant digits also round the retention time, and the source
    // truncates the precursor intensity.
    let document = std::fs::read_to_string(&lossy).unwrap();
    assert!(document.contains("retentionTime=\"PT1.5S\""), "{document}");

    let explicit = WriteOptions {
        precision: PeakPrecision::Float32,
        ..WriteOptions::default()
    };
    let path = directory.path().join("explicit32.mzXML");
    mzxml::store_with_options(&path, &experiment, &explicit).unwrap();
    assert_eq!(
        mzxml::load(&path).unwrap().spectra[0].peaks[0].mz,
        f64::from(WIDE_MZ as f32)
    );
}

/// The source writes `(int)precursor.getIntensity()`; the lossless default
/// keeps the fraction.
#[test]
fn precursor_intensity_truncation_is_opt_in() {
    let mut experiment = MSExperiment::new();
    let mut spectrum = MSSpectrum {
        ms_level: 2,
        rt: 1.0,
        native_id: "scan=1".into(),
        ..MSSpectrum::default()
    };
    spectrum.peaks.push(Peak1D {
        mz: 100.0,
        intensity: 1.0,
    });
    spectrum.precursors.push(Precursor {
        mz: 250.0,
        intensity: 12.75,
        charge: 2,
        ..Precursor::default()
    });
    experiment.spectra.push(spectrum);
    let directory = TempDir::new_in(std::env::temp_dir(), false).unwrap();

    let lossless = directory.path().join("precursor.mzXML");
    mzxml::store(&lossless, &experiment).unwrap();
    let reloaded = mzxml::load(&lossless).unwrap();
    assert_eq!(reloaded.spectra[0].precursors[0].intensity, 12.75);

    let lossy = directory.path().join("precursor_source.mzXML");
    mzxml::store_with_options(&lossy, &experiment, &WriteOptions::source()).unwrap();
    let truncated = mzxml::load(&lossy).unwrap();
    assert_eq!(truncated.spectra[0].precursors[0].intensity, 12.0);
}

/// The native zlib writer round-trips through the source-defined reader path.
#[test]
fn zlib_written_peaks_round_trip() {
    let original = fixture_1();
    let directory = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = directory.path().join("compressed.mzXML");
    let options = WriteOptions {
        zlib_compression: true,
        ..WriteOptions::default()
    };
    mzxml::store_with_options(&path, &original, &options).unwrap();
    let document = std::fs::read_to_string(&path).unwrap();
    assert!(document.contains("compressionType=\"zlib\""), "{document}");
    assert!(!document.contains("compressedLen=\"0\""));
    assert_eq!(mzxml::load(&path).unwrap(), original);
}

/// MaxQuant compatibility adds the four statistics attributes, forces
/// `scanType="Full"`, breaks the `<peaks>` tag across lines, falls back to a
/// `CID` activation method, writes an index even when it was not requested and
/// skips empty spectra.
#[test]
fn maxquant_compatibility_changes_the_written_document() {
    let mut experiment = fixture_1();
    experiment.spectra.push(MSSpectrum {
        ms_level: 1,
        rt: 200.0,
        native_id: "scan=14".into(),
        ..MSSpectrum::default()
    });
    let directory = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = directory.path().join("maxquant.mzXML");
    let options = WriteOptions {
        force_mq_compatibility: true,
        write_index: false,
        ..WriteOptions::source()
    };
    mzxml::store_with_options(&path, &experiment, &options).unwrap();
    let document = std::fs::read_to_string(&path).unwrap();
    assert!(document.contains("lowMz="), "{document}");
    assert!(document.contains("highMz="));
    assert!(document.contains("basePeakIntensity="));
    assert!(document.contains("totIonCurrent="));
    assert!(document.contains("scanType=\"Full\""));
    assert!(document.contains("<peaks precision=\"32\"\n byteOrder=\"network\""));
    assert!(document.contains("activationMethod=\"CID\""));
    assert!(document.contains("<indexOffset>"));
    // The empty spectrum is skipped, so only the four original scans remain.
    let reloaded = mzxml::load(&path).unwrap();
    assert_eq!(reloaded.len(), 4);
    // An Xcalibur acquisition software name forces a Thermo manufacturer.
    let mut thermo = MSExperiment::new();
    thermo.settings.instrument.software.name = "Xcalibur".into();
    thermo.spectra.push(MSSpectrum {
        ms_level: 1,
        rt: 0.0,
        native_id: "scan=1".into(),
        ..MSSpectrum::default()
    });
    let path = directory.path().join("thermo.mzXML");
    mzxml::store_with_options(&path, &thermo, &WriteOptions::default()).unwrap();
    assert!(
        std::fs::read_to_string(&path)
            .unwrap()
            .contains("value=\"Thermo Scientific\"")
    );
}

/// MaxQuant compatibility needs sorted m/z for its `lowMz`/`highMz` attributes.
/// The source logs a non-fatal error and writes the first and last stored peak
/// anyway, producing wrong values; this port refuses.
#[test]
fn maxquant_compatibility_refuses_an_unsorted_spectrum() {
    let mut experiment = MSExperiment::new();
    let mut spectrum = MSSpectrum {
        ms_level: 1,
        rt: 0.0,
        native_id: "scan=1".into(),
        ..MSSpectrum::default()
    };
    spectrum.peaks = vec![
        Peak1D {
            mz: 200.0,
            intensity: 1.0,
        },
        Peak1D {
            mz: 100.0,
            intensity: 2.0,
        },
    ];
    experiment.spectra.push(spectrum);
    let mut output = Vec::new();
    let options = WriteOptions {
        force_mq_compatibility: true,
        ..WriteOptions::default()
    };
    let error = mzxml::write_with_options(&mut output, &experiment, &options).unwrap_err();
    assert!(matches!(error, Error::UnsortedData), "{error:?}");
    assert!(output.is_empty());
}

/// The index trailer records the byte offset of each `<scan` tag, and the
/// `<indexOffset>` the position of the `<index>` element itself.
#[test]
fn the_index_trailer_points_at_the_scan_tags() {
    let original = fixture_1();
    let mut output = Vec::new();
    mzxml::write(&mut output, &original).unwrap();
    let document = String::from_utf8(output).unwrap();
    for id in [10, 11, 12, 13] {
        let marker = format!("<offset id = \"{id}\" >");
        let start = document.find(&marker).expect("offset entry") + marker.len();
        let end = start + document[start..].find('<').expect("offset end");
        let offset: usize = document[start..end].parse().unwrap();
        assert!(
            document[offset..].starts_with(&format!("<scan num=\"{id}\"")),
            "offset {offset} for id {id} points at {:?}",
            &document[offset..offset + 24]
        );
    }
    let marker = "<indexOffset>";
    let start = document.find(marker).unwrap() + marker.len();
    let end = start + document[start..].find('<').unwrap();
    let offset: usize = document[start..end].parse().unwrap();
    assert!(document[offset..].starts_with("<index name = \"scan\" >"));
    // Suppressing the index removes both elements.
    let mut output = Vec::new();
    let options = WriteOptions {
        write_index: false,
        ..WriteOptions::default()
    };
    mzxml::write_with_options(&mut output, &original, &options).unwrap();
    let document = String::from_utf8(output).unwrap();
    assert!(!document.contains("<index "));
    assert!(!document.contains("<indexOffset>"));
}

/// An empty experiment still writes a schema-plausible document: the source
/// bumps a zero `scanCount` to one and emits a placeholder `<parentFile>` and
/// `<dataProcessing>`.
#[test]
fn an_empty_experiment_writes_a_scan_count_of_one() {
    let mut output = Vec::new();
    mzxml::write(&mut output, &MSExperiment::new()).unwrap();
    let document = String::from_utf8(output).unwrap();
    assert!(document.contains("scanCount=\"1\""), "{document}");
    assert!(document.contains("<parentFile fileName=\"\" fileType=\"processedData\""));
    assert!(document.contains("<software type=\"processing\" name=\"\" version=\"\"/>"));
    assert!(!document.contains("<msInstrument>"));
    assert_eq!(
        mzxml::read(std::io::Cursor::new(document)).unwrap().len(),
        0
    );
}

/// Native IDs decide the written `num`: `scan=`-prefixed numbers and bare
/// numbers are preserved, anything else is renumbered from one.
#[test]
fn native_ids_are_preserved_or_renumbered() {
    let make = |ids: [&str; 2]| {
        let mut experiment = MSExperiment::new();
        for id in ids {
            let mut spectrum = MSSpectrum {
                ms_level: 1,
                rt: 0.0,
                native_id: id.to_owned(),
                ..MSSpectrum::default()
            };
            spectrum.peaks.push(Peak1D {
                mz: 100.0,
                intensity: 1.0,
            });
            experiment.spectra.push(spectrum);
        }
        experiment
    };
    for (ids, expected) in [
        (["scan=7", "scan=9"], ["scan=7", "scan=9"]),
        (["7", "9"], ["scan=7", "scan=9"]),
        (["spectrum=a", "spectrum=b"], ["scan=1", "scan=2"]),
        (["", ""], ["scan=1", "scan=2"]),
    ] {
        let mut output = Vec::new();
        mzxml::write(&mut output, &make(ids)).unwrap();
        let reloaded = mzxml::read(std::io::Cursor::new(output)).unwrap();
        assert_eq!(
            reloaded
                .spectra
                .iter()
                .map(|s| s.native_id.as_str())
                .collect::<Vec<_>>(),
            expected.to_vec(),
            "{ids:?}"
        );
    }
}

/// A spectrum whose MS level mzXML cannot indent, or a nonfinite retention
/// time, is refused before anything is written.
#[test]
fn the_writer_preflight_refuses_unrepresentable_records() {
    let mut experiment = MSExperiment::new();
    experiment.spectra.push(MSSpectrum {
        ms_level: 1_000,
        rt: 0.0,
        native_id: "scan=1".into(),
        ..MSSpectrum::default()
    });
    let mut output = Vec::new();
    let error = mzxml::write(&mut output, &experiment).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    assert!(output.is_empty());

    experiment.spectra[0].ms_level = 1;
    experiment.spectra[0].rt = f64::NAN;
    let error = mzxml::write(&mut output, &experiment).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error:?}");
    assert!(output.is_empty());

    experiment.spectra[0].rt = 0.0;
    experiment.spectra[0].native_id = "scan=\u{1}".into();
    let error = mzxml::write(&mut output, &experiment).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    assert!(output.is_empty());
}

/// A non-SHA-1 `fileSha1` is kept verbatim but its algorithm is downgraded, so
/// the returned record still passes `SourceFile::validate`, and the writer
/// replaces it with the source's 40-zero placeholder.
#[test]
fn a_malformed_checksum_is_downgraded_and_replaced_on_write() {
    let document = "<mzXML><msRun><parentFile fileName=\"a.raw\" fileType=\"RAWData\" \
                    fileSha1=\"notahash\"/></msRun></mzXML>";
    let mut report = mzxml::ReadReport::default();
    let experiment = mzxml::read_with_report(
        std::io::Cursor::new(document),
        &ReadOptions::default(),
        &mut report,
    )
    .unwrap();
    assert_eq!(experiment.settings.source_files[0].checksum, "notahash");
    assert_eq!(
        experiment.settings.source_files[0].checksum_type,
        ChecksumType::Unknown
    );
    experiment.validate().unwrap();
    assert!(
        report
            .diagnostics
            .iter()
            .any(|line| line.contains("is not a SHA-1 digest")),
        "{:?}",
        report.diagnostics
    );
    let mut output = Vec::new();
    mzxml::write(&mut output, &experiment).unwrap();
    let written = String::from_utf8(output).unwrap();
    assert!(
        written.contains("fileSha1=\"0000000000000000000000000000000000000000\""),
        "{written}"
    );
}

/// `msRun`'s `startTime`/`endTime` are the first and last spectrum's retention
/// times, not the minimum and maximum, and a negative one keeps its sign. The
/// source writes `PT-1S` there, which `xs:duration` does not allow.
#[test]
fn ms_run_times_are_the_first_and_last_scan() {
    let mut experiment = MSExperiment::new();
    for rt in [50.0, 10.0, -1.0] {
        experiment.spectra.push(MSSpectrum {
            ms_level: 1,
            rt,
            native_id: "scan=1".into(),
            ..MSSpectrum::default()
        });
    }
    let mut output = Vec::new();
    mzxml::write(&mut output, &experiment).unwrap();
    let document = String::from_utf8(output).unwrap();
    assert!(document.contains("startTime=\"PT50S\""), "{document}");
    assert!(document.contains("endTime=\"-PT1S\""));
    assert!(!document.contains("PT-1S"));
}

/// A default `PeakFileOptions` carries the two writer switches the source
/// consults; `WriteOptions::from_peak_options` copies them.
#[test]
fn write_options_can_be_derived_from_peak_file_options() {
    let mut peaks = openms::format::PeakFileOptions::new();
    peaks.force_mq_compatibility = true;
    peaks.write_index = false;
    peaks.zlib_compression = true;
    let options = WriteOptions::from_peak_options(&peaks);
    assert!(options.force_mq_compatibility);
    assert!(!options.write_index);
    assert!(options.zlib_compression);
    // The lossless defaults survive.
    assert_eq!(options.precision, PeakPrecision::Float64);
    assert_eq!(options.significant_digits, None);
}

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

fn open(path: PathBuf) -> std::io::BufReader<std::fs::File> {
    std::io::BufReader::new(std::fs::File::open(path).unwrap())
}

fn scan_document(count: &str, peak_attributes: &str, payload: &str) -> String {
    format!(
        "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n<mzXML><msRun>\
         <scan num=\"1\" msLevel=\"1\" {count} polarity=\"any\" retentionTime=\"PT1S\">\
         <peaks {peak_attributes}>{payload}</peaks></scan></msRun></mzXML>"
    )
}

/// Standard base64 without a dependency on the crate's private helpers.
fn base64_encode(bytes: &[u8]) -> String {
    const ALPHABET: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for group in bytes.chunks(3) {
        let b0 = u32::from(group[0]);
        let b1 = group.get(1).copied().map_or(0, u32::from);
        let b2 = group.get(2).copied().map_or(0, u32::from);
        let word = (b0 << 16) | (b1 << 8) | b2;
        out.push(char::from(ALPHABET[(word >> 18) as usize & 63]));
        out.push(char::from(ALPHABET[(word >> 12) as usize & 63]));
        if group.len() > 1 {
            out.push(char::from(ALPHABET[(word >> 6) as usize & 63]));
        } else {
            out.push('=');
        }
        if group.len() > 2 {
            out.push(char::from(ALPHABET[word as usize & 63]));
        } else {
            out.push('=');
        }
    }
    out
}

#[test]
fn the_local_base64_helper_agrees_with_the_fixture_payloads() {
    // The fixture's own scan-10 payload, so the generated long document is
    // encoded the same way the upstream files are.
    assert_eq!(
        base64_encode(&[0x42, 0xf0, 0x00, 0x00, 0x42, 0xc8, 0x00, 0x00]),
        "QvAAAELIAAA="
    );
    assert_eq!(base64_encode(b"a"), "YQ==");
    assert_eq!(base64_encode(b"ab"), "YWI=");
    assert_eq!(base64_encode(b"abc"), "YWJj");
}

/// The `TransformOptions` struct is public API; a default one behaves like the
/// two-pass source call with `skip_full_count = false`.
#[test]
fn default_transform_options_run_both_passes() {
    let mut consumer = TicConsumer::default();
    let report = mzxml::transform_with_options(
        data("MzXMLFile_1.mzXML"),
        &mut consumer,
        &TransformOptions::default(),
    )
    .unwrap();
    assert_eq!(report.expected_spectra, 4);
    assert_eq!(report.delivered, 4);
    assert_eq!(consumer.spectra, 4);
}
