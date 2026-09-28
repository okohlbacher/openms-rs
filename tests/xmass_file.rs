//! `XMassFile`, `AcqusHandler` and `FidHandler` against the upstream class test
//! `XMassFile_test.cpp` (core pin bc9cc12) and its data
//! (`src/tests/class_tests/openms/data/XMassFile_test/{acqus,fid}` and
//! `XMassFile_test_data.txt`), copied unchanged into `tests/data/`.
//!
//! Evidence: the 80478 expected `(m/z, intensity)` pairs are the retained
//! upstream expectation file, compared with the class test's own
//! `TOLERANCE_ABSOLUTE(0.001)` and `TEST_REAL_SIMILAR` rule. The metadata
//! assertions are read from the fixture's `acqus` text and the source
//! assignments in `XMassFile.h` (source review).

use openms::format::acqus_handler::AcqusHandler;
use openms::format::fid_handler::FidHandler;
use openms::format::xmass_file::{LoadOptions, TruncatedFid, XMassFile};
use openms::kernel::SpectrumType;
use openms::metadata::{
    AnalyzerType, InletType, IonizationMethod, Polarity, ProcessingAction, ScanMode,
};
use openms::{Error, MSExperiment, MSSpectrum};
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// `ClassTest::isRealSimilar` with `TOLERANCE_ABSOLUTE(0.001)` and the default
/// relative tolerance of 1e-5: similar when either bound holds.
fn real_similar(actual: f64, expected: f64) -> bool {
    const ABSOLUTE: f64 = 0.001;
    const RATIO: f64 = 1.0 + 1e-5;
    if actual == expected {
        return true;
    }
    if (actual - expected).abs() <= ABSOLUTE {
        return true;
    }
    if actual == 0.0 || expected == 0.0 {
        return false;
    }
    let ratio = actual / expected;
    let ratio = if ratio < 1.0 { 1.0 / ratio } else { ratio };
    ratio > 0.0 && ratio <= RATIO
}

fn scratch(case: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms-xmass-{}-{case}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn load_matches_the_class_test_expectation_file() {
    let spectrum = XMassFile::new().load(data("XMassFile_test/fid")).unwrap();
    // TEST_EQUAL(s.size(), 80478)
    assert_eq!(spectrum.peaks.len(), 80478);

    let expected = std::fs::read_to_string(data("XMassFile_test_data.txt")).unwrap();
    let rows: Vec<&str> = expected.lines().collect();
    // TEST_EQUAL((file.end() - file.begin()), 80478)
    assert_eq!(rows.len(), 80478);
    for (index, (peak, row)) in spectrum.peaks.iter().zip(&rows).enumerate() {
        let values: Vec<f64> = row
            .split(',')
            .map(|item| item.trim().parse().unwrap())
            .collect();
        assert_eq!(values.len(), 2, "row {index}");
        assert!(
            real_similar(peak.mz, values[0]),
            "m/z of point {index}: {} vs {}",
            peak.mz,
            values[0]
        );
        assert!(
            real_similar(f64::from(peak.intensity), values[1]),
            "intensity of point {index}: {} vs {}",
            peak.intensity,
            values[1]
        );
    }
}

#[test]
fn missing_file_is_not_found() {
    // TEST_EXCEPTION(Exception::FileNotFound, f.load("data_Idontexist", s))
    let error = XMassFile::new().load("data_Idontexist").unwrap_err();
    match error {
        Error::Io(error) => assert_eq!(error.kind(), std::io::ErrorKind::NotFound),
        other => panic!("expected a not-found I/O error, got {other:?}"),
    }
}

#[test]
fn store_is_not_implemented() {
    // TEST_EXCEPTION(Exception::NotImplemented, XMassFile().store(...))
    assert!(matches!(
        XMassFile::new().store("", &MSSpectrum::default()),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn load_sets_the_source_metadata() {
    let path = data("XMassFile_test/fid");
    let spectrum = XMassFile::new().load(&path).unwrap();
    assert_eq!(spectrum.rt, 0.0);
    assert_eq!(spectrum.ms_level, 1);
    assert_eq!(
        spectrum.name,
        "Xmass analysis file <779007CD-D8F3-4E54-909F-1E998217AA41>"
    );
    assert_eq!(
        spectrum.native_id,
        "spectrum=xsd:779007CD-D8F3-4E54-909F-1E998217AA41"
    );
    assert_eq!(spectrum.spectrum_type, SpectrumType::Profile);
    assert_eq!(
        spectrum.instrument_settings.scan_mode,
        ScanMode::MassSpectrum
    );
    assert!(!spectrum.instrument_settings.zoom_scan);
    // `##.IONIZATION MODE=  LD+`
    assert_eq!(spectrum.instrument_settings.polarity, Polarity::Positive);
    assert_eq!(
        spectrum.acquisition_info.method_of_combination,
        "Sum of 300 raw spectrum"
    );
    let directory = path.to_str().unwrap();
    let directory = &directory[..directory.len() - 3];
    assert_eq!(spectrum.source_file.name, "fid");
    assert_eq!(spectrum.source_file.path, directory);
    assert_eq!(
        spectrum.source_file.size_mb,
        (4.0 * 80478.0 / 1024.0 / 1024.0) as f32
    );
    assert_eq!(spectrum.source_file.file_type, "Xmass analysis file (fid)");

    assert_eq!(spectrum.data_processing.len(), 1);
    let processing = &spectrum.data_processing[0];
    assert_eq!(processing.software.name, "FlexControl");
    // `##$FCVer= <flexControl 2.4.41.0>`
    assert_eq!(processing.software.version, "2.4.41.0");
    assert_eq!(
        processing
            .software
            .cv_terms
            .metadata
            .get("Acquisition method")
            .map(ToString::to_string)
            .as_deref(),
        Some(r"C:\Methods\flexControlMethods\LP_PepMix.par")
    );
    let actions: Vec<ProcessingAction> = processing.actions.iter().copied().collect();
    assert_eq!(actions.len(), 3);
    assert!(processing.actions.contains(&ProcessingAction::Smoothing));
    assert!(
        processing
            .actions
            .contains(&ProcessingAction::BaselineReduction)
    );
    assert!(
        processing
            .actions
            .contains(&ProcessingAction::MzCalibration)
    );
    assert!(processing.completion_time.is_some());
}

#[test]
fn import_experimental_settings_reads_acqus() {
    let mut exp = MSExperiment::new();
    exp.settings.comment = "kept".into();
    XMassFile::new()
        .import_experimental_settings(data("XMassFile_test/fid"), &mut exp)
        .unwrap();
    let instrument = &exp.settings.instrument;
    assert_eq!(instrument.name, "Bruker Flex Series");
    assert_eq!(
        instrument.vendor,
        "XMASS, Bruker-Daltonics (USA) and Bruker-Daltonik (GER)"
    );
    assert_eq!(instrument.model, "222024.00050");
    assert_eq!(instrument.ion_sources.len(), 1);
    let source = &instrument.ion_sources[0];
    // `##.INLET=  DIRECT`: direct inlet, ionization method left unknown.
    assert_eq!(source.inlet_type, InletType::Direct);
    assert_eq!(source.ionization_method, IonizationMethod::default());
    assert_eq!(source.polarity, Polarity::Positive);
    assert_eq!(source.order, 0);
    assert_eq!(
        source
            .metadata
            .get("MALDI target reference")
            .map(ToString::to_string)
            .as_deref(),
        Some("T_0209518_0000137_11")
    );
    assert_eq!(instrument.mass_analyzers.len(), 1);
    assert_eq!(
        instrument.mass_analyzers[0].analyzer_type,
        AnalyzerType::Tof
    );
    // `##$AQ_DATE= <2009-10-21 12:51:09>`
    assert_eq!(exp.settings.date_time.get(), "2009-10-21 12:51:09");
    assert_eq!(exp.settings.comment, "kept");
}

#[test]
fn handlers_read_the_class_test_files() {
    let acqus = AcqusHandler::open(data("XMassFile_test/acqus")).unwrap();
    assert_eq!(acqus.size(), 80478);
    assert_eq!(acqus.param("$NoSHOTS"), "300");
    assert_eq!(acqus.param("does not exist"), "");
    // The first and last rows of the upstream expectation file.
    assert!(real_similar(acqus.position(0), 500.227));
    assert!(real_similar(
        acqus.param("mzMin").parse::<f64>().unwrap(),
        500.227
    ));
    let expected = std::fs::read_to_string(data("XMassFile_test_data.txt")).unwrap();
    let last_mz: f64 = expected
        .lines()
        .last()
        .unwrap()
        .split(',')
        .next()
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    assert!(real_similar(acqus.position(80477), last_mz));
    assert!(real_similar(
        acqus.param("mzMax").parse::<f64>().unwrap(),
        last_mz
    ));

    let mut fid = FidHandler::open(data("XMassFile_test/fid")).unwrap();
    assert_eq!(fid.index(), 0);
    // First row: `500.227, 317.0`.
    assert_eq!(fid.intensity().unwrap(), 317);
    assert_eq!(fid.index(), 1);
    assert_eq!(fid.intensity().unwrap(), 324);
}

#[test]
fn truncated_fid_is_refused_by_default_and_padded_on_request() {
    let dir = scratch("truncated");
    std::fs::copy(data("XMassFile_test/acqus"), dir.join("acqus")).unwrap();
    let full = std::fs::read(data("XMassFile_test/fid")).unwrap();
    // Ten complete points and two stray bytes.
    std::fs::write(dir.join("fid"), &full[..42]).unwrap();

    let error = XMassFile::new().load(dir.join("fid")).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error:?}");

    let spectrum = XMassFile::new()
        .load_with_options(dir.join("fid"), &LoadOptions::source())
        .unwrap();
    assert_eq!(spectrum.peaks.len(), 80478);
    assert_eq!(spectrum.peaks[0].intensity, 317.0);
    // Point 10 is the two stray bytes in the low-order positions.
    let partial = u32::from(full[40]) | (u32::from(full[41]) << 8);
    assert_eq!(spectrum.peaks[10].intensity, partial as f32);
    assert!(spectrum.peaks[11..].iter().all(|p| p.intensity == 0.0));
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn oversized_td_is_refused_before_allocation() {
    let dir = scratch("ceiling");
    std::fs::copy(data("XMassFile_test/acqus"), dir.join("acqus")).unwrap();
    std::fs::copy(data("XMassFile_test/fid"), dir.join("fid")).unwrap();
    let options = LoadOptions {
        truncated_fid: TruncatedFid::SourcePadding,
        max_points: 1000,
    };
    assert!(matches!(
        XMassFile::new().load_with_options(dir.join("fid"), &options),
        Err(Error::InvalidValue(_))
    ));
    let _ = std::fs::remove_dir_all(&dir);
}
