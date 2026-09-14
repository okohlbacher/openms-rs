// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Spectrum- and scan-level ion mobility, the `MS:1000525` representation reset
//! and unit-bearing ion-mobility arrays in the native mzML reader and writer.
//!
//! Every loaded value is compared with `data/mzml_mobility/a3_format_io_oracle.tsv`,
//! the output of the unmodified product-sdk `FileHandler::loadExperiment` on the
//! same bytes (oracle-generated, tier 1 executed differential; driver, inputs
//! and hashes in `../oracle/a3-format-io/`). Class-test literals come from the
//! pinned `MzMLFile_test.cpp`. `docs/MZML_MOBILITY_SUPPORT.md` records the
//! mapping and the native differences these tests pin.
#![cfg(feature = "mzml")]

use openms::Error;
use openms::format::{FileHandler, FileType, PeakFileOptions, mzml};
use openms::kernel::{DataArray, MSExperiment, MSSpectrum, Peak1D, SpectrumType};
use openms::metadata::{DriftTimeUnit, ImTypes, IonMobilityFormat};
use std::io::Cursor;

const ORACLE: &str = include_str!("data/mzml_mobility/a3_format_io_oracle.tsv");

fn fixture(name: &str) -> String {
    format!(
        "{}/tests/data/mzml_mobility/{name}",
        env!("CARGO_MANIFEST_DIR")
    )
}

/// Oracle records of one kind for one input label, split into fields.
fn rows(record: &str, label: &str) -> Vec<Vec<&'static str>> {
    ORACLE
        .lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|fields| fields[0] == record && fields.get(1) == Some(&label))
        .collect()
}

/// Parse a C `printf("%a")` hexadecimal float exactly.
fn hex(text: &str) -> f64 {
    let (negative, body) = match text.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, text),
    };
    let body = body.strip_prefix("0x").expect("hexadecimal float");
    let (mantissa, exponent) = body.split_once('p').expect("binary exponent");
    let exponent: i32 = exponent.parse().expect("decimal exponent");
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let digits = u64::from_str_radix(&format!("{whole}{fraction}"), 16).expect("hex digits");
    let scale = exponent - 4 * i32::try_from(fraction.len()).expect("short fraction");
    let magnitude = digits as f64 * 2f64.powi(scale);
    if negative { -magnitude } else { magnitude }
}

/// Load as `FileHandler().loadExperiment(path, exp, {MZML})` with source defaults.
fn load(name: &str) -> openms::Result<MSExperiment> {
    FileHandler::load_experiment_with_options(
        fixture(name),
        &[FileType::MzMl],
        &PeakFileOptions::default(),
    )
}

fn read(xml: &str) -> openms::Result<MSExperiment> {
    mzml::read(Cursor::new(xml.as_bytes()))
}

fn unit<T>(array: &DataArray<T>) -> Option<&str> {
    array
        .metadata
        .get("unit_accession")
        .map(|value| value.as_str().expect("string unit accession"))
}

/// The driver's `kind:name|unit_accession|size` rendering of every data array.
fn arrays(spectrum: &MSSpectrum) -> String {
    fn describe<T>(out: &mut Vec<String>, kind: &str, arrays: &[DataArray<T>]) {
        for array in arrays {
            out.push(format!(
                "{kind}:{}|{}|{}",
                array.name,
                unit(array).unwrap_or("-"),
                array.data.len()
            ));
        }
    }
    let mut out = Vec::new();
    describe(&mut out, "f", &spectrum.float_data_arrays);
    describe(&mut out, "i", &spectrum.integer_data_arrays);
    describe(&mut out, "s", &spectrum.string_data_arrays);
    if out.is_empty() {
        "-".into()
    } else {
        out.join(";")
    }
}

/// Compare every loaded spectrum with the oracle. `native_drift` lists spectra
/// whose spectrum-level drift time is a documented native difference.
fn assert_oracle(label: &str, experiment: &MSExperiment, native_drift: &[usize]) {
    let loaded = rows("loaded", label);
    assert_eq!(loaded.len(), 1, "{label}: one oracle load record");
    assert_eq!(
        experiment.spectra.len().to_string(),
        loaded[0][2],
        "{label}"
    );
    assert_eq!(
        experiment.chromatograms.len().to_string(),
        loaded[0][3],
        "{label} chromatograms"
    );
    let expected = rows("spectrum", label);
    assert_eq!(expected.len(), experiment.spectra.len(), "{label}");
    for (index, (spectrum, row)) in experiment.spectra.iter().zip(&expected).enumerate() {
        let context = format!("{label} spectrum {index}");
        assert_eq!(row[2], index.to_string(), "{context}");
        assert_eq!(spectrum.native_id, row[3], "{context}");
        assert_eq!(spectrum.ms_level.to_string(), row[4], "{context}");
        assert_eq!(format!("{:?}", spectrum.spectrum_type), row[5], "{context}");
        if !native_drift.contains(&index) {
            assert_eq!(
                spectrum.drift_time.to_bits(),
                hex(row[6]).to_bits(),
                "{context} drift time"
            );
            assert_eq!(spectrum.drift_time_unit.name(), row[7], "{context} unit");
            assert_eq!(
                ImTypes::determine_im_format(spectrum).name(),
                row[8],
                "{context} format"
            );
        }
        assert_eq!(spectrum.len().to_string(), row[9], "{context} peaks");
        let precursors: Vec<_> = if row[10] == "-" {
            Vec::new()
        } else {
            row[10].split(',').collect()
        };
        assert_eq!(spectrum.precursors.len(), precursors.len(), "{context}");
        for (precursor, text) in spectrum.precursors.iter().zip(precursors) {
            let (drift, unit) = text.split_once(':').expect("drift:unit");
            assert_eq!(
                precursor.drift_time.unwrap_or(-1.0).to_bits(),
                hex(drift).to_bits(),
                "{context} precursor drift time"
            );
            assert_eq!(precursor.drift_time_unit.name(), unit, "{context}");
        }
        assert_eq!(arrays(spectrum), row[11], "{context} arrays");
    }
}

#[test]
fn upstream_faims_interleaved_file_reads_spectrum_level_compensation_voltages() {
    let experiment = load("FAIMS_CV-60C_V-45_Interleaved.mzML").unwrap();
    assert_eq!(experiment.spectra.len(), 12);
    let mut voltages: Vec<f64> = experiment.spectra.iter().map(|s| s.drift_time).collect();
    voltages.sort_by(f64::total_cmp);
    voltages.dedup();
    assert_eq!(voltages, [-60.0, -45.0]);
    assert!(
        experiment
            .spectra
            .iter()
            .all(|s| s.drift_time_unit == DriftTimeUnit::FaimsCompensationVoltage)
    );
    assert_oracle("FAIMS_CV-60C_V-45_Interleaved", &experiment, &[]);
}

#[test]
fn upstream_faims_test_data_reads_the_scan_level_compensation_voltage() {
    let experiment = load("FAIMS_test_data.mzML").unwrap();
    let first = &experiment.spectra[0];
    assert_eq!(first.drift_time, -65.0);
    assert_eq!(
        first.drift_time_unit,
        DriftTimeUnit::FaimsCompensationVoltage
    );
    assert_eq!(
        ImTypes::determine_im_format(first),
        IonMobilityFormat::PerSpectrum
    );
    assert_oracle("FAIMS_test_data", &experiment, &[]);
}

#[test]
fn class_test_fixture_scan_drift_time_and_the_precursor_propagation_gap() {
    let path = format!(
        "{}/tests/data/mzml_validator/MzMLFile_1.mzML",
        env!("CARGO_MANIFEST_DIR")
    );
    let experiment = FileHandler::load_experiment_with_options(
        &path,
        &[FileType::MzMl],
        &PeakFileOptions::default(),
    )
    .unwrap();
    // MzMLFile_test.cpp:368-369, 451-452 and 464-465.
    assert_eq!(experiment.spectra[0].drift_time, 7.1);
    assert_eq!(
        experiment.spectra[0].drift_time_unit,
        DriftTimeUnit::Millisecond
    );
    let precursors = &experiment.spectra[1].precursors;
    assert_eq!(precursors[0].drift_time, Some(8.1));
    assert_eq!(precursors[0].drift_time_unit, DriftTimeUnit::Millisecond);
    assert_eq!(precursors[1].drift_time, None);
    assert_eq!(precursors[1].drift_time_unit, DriftTimeUnit::None);
    // MzMLFile_test.cpp:418-421 also expects the selected-ion drift time 8.1 on
    // the spectrum, and the oracle shows it. This port does not copy it there:
    // tests/precursor_workflow.rs asserts the unpropagated round trip. See
    // "Native differences" in docs/MZML_MOBILITY_SUPPORT.md.
    assert_eq!(hex(rows("spectrum", "MzMLFile_1")[1][6]), 8.1);
    assert_eq!(experiment.spectra[1].drift_time, -1.0);
    assert_eq!(experiment.spectra[1].drift_time_unit, DriftTimeUnit::None);
    assert_oracle("MzMLFile_1", &experiment, &[1]);
    // The scan writer keeps the scan-level value through a round trip.
    let mut bytes = Vec::new();
    mzml::write(&mut bytes, &experiment).unwrap();
    let reread = mzml::read(bytes.as_slice()).unwrap();
    assert_eq!(reread.spectra, experiment.spectra);
}

#[test]
fn scan_and_spectrum_level_mobility_terms_match_the_oracle() {
    let experiment = load("scan_mobility.mzML").unwrap();
    use DriftTimeUnit as D;
    let expected = [
        (7.5, D::Millisecond),
        (0.85, D::InverseReducedMobility),
        (-35.0, D::FaimsCompensationVoltage),
        (350.5, D::CollisionCrossSection),
        (-50.0, D::FaimsCompensationVoltage), // spectrum-level FAIMS voltage
        (-1.0, D::None),                      // spectrum-level MS:1002476 is no source route
        (-1.0, D::None),
        (-1.0, D::FaimsCompensationVoltage), // a voltage equal to the -1 sentinel
        (-40.0, D::FaimsCompensationVoltage), // no unit attributes
        (2.5, D::Millisecond),
        (7.5, D::Millisecond),                // identical repeat
        (-42.0, D::FaimsCompensationVoltage), // second scan
        (-45.0, D::FaimsCompensationVoltage), // legacy UO:000218 spelling
    ];
    let actual: Vec<_> = experiment
        .spectra
        .iter()
        .map(|s| (s.drift_time, s.drift_time_unit))
        .collect();
    assert_eq!(actual, expected);
    assert_eq!(
        ImTypes::determine_im_format(&experiment.spectra[7]),
        IonMobilityFormat::None
    );
    assert_oracle("scan_mobility", &experiment, &[]);
}

#[test]
fn spectrum_representation_reset_matches_the_oracle() {
    let experiment = load("spectrum_representation.mzML").unwrap();
    assert_eq!(
        experiment
            .spectra
            .iter()
            .map(|s| s.spectrum_type)
            .collect::<Vec<_>>(),
        [
            SpectrumType::Unknown,
            SpectrumType::Profile,
            SpectrumType::Unknown,
            SpectrumType::Centroid,
            SpectrumType::Unknown
        ]
    );
    assert_oracle("spectrum_representation", &experiment, &[]);
}

#[test]
fn mobility_arrays_with_and_without_units_are_per_peak() {
    let experiment = load("im_arrays.mzML").unwrap();
    for spectrum in &experiment.spectra {
        assert_eq!(
            ImTypes::determine_im_format(spectrum),
            IonMobilityFormat::PerPeak
        );
    }
    let units: Vec<_> = experiment
        .spectra
        .iter()
        .map(|s| unit(&s.float_data_arrays[0]))
        .collect();
    assert_eq!(
        units,
        [
            None,
            Some("UO:0000028"),
            Some("MS:1002814"),
            Some("UO:0000028"),
            None
        ]
    );
    assert_oracle("im_arrays", &experiment, &[]);
}

#[test]
fn feature_finder_input_matches_the_oracle() {
    let experiment = load("FeatureFinderCentroided_1_input.mzML").unwrap();
    assert_oracle("FeatureFinderCentroided_1_input", &experiment, &[]);
}

#[test]
fn source_lenient_inputs_stay_explicit_native_errors() {
    // The oracle loads all three: a later repeated value wins, units on mobility
    // terms are ignored, and every auxiliary array keeps its unit.
    for (label, message) in [
        (
            "scan_mobility_conflict",
            "conflicting spectrum ion mobility values",
        ),
        (
            "scan_mobility_unit_mismatch",
            "spectrum mobility unit does not match its typed quantity",
        ),
        (
            "im_arrays_non_mobility_unit",
            "units on auxiliary arrays are not represented",
        ),
    ] {
        assert_eq!(rows("loaded", label).len(), 1, "{label} oracle load");
        let error = load(&format!("{label}.mzML")).unwrap_err();
        assert!(
            matches!(&error, Error::Unsupported(text) if text.contains(message)),
            "{label}: {error}"
        );
    }
    let conflict = rows("spectrum", "scan_mobility_conflict");
    assert_eq!(hex(conflict[0][6]), 8.5);
    assert_eq!(hex(conflict[1][6]), -60.0);
    assert_eq!(
        rows("spectrum", "im_arrays_non_mobility_unit")[0][11],
        "f:review array|UO:0000010|3"
    );
}

#[test]
fn mobility_units_values_and_repeats_are_checked() {
    let base = include_str!("data/mzml_mobility/scan_mobility.mzML");
    let faims = r#"accession="MS:1001581" name="FAIMS compensation voltage" value="-35" unitCvRef="UO" unitAccession="UO:0000218""#;
    let millisecond = r#"accession="MS:1002476" name="ion mobility drift time" value="7.5" unitCvRef="UO" unitAccession="UO:0000028""#;
    assert!(base.contains(faims) && base.contains(millisecond));
    // The legacy spelling means volts, for the FAIMS voltage only.
    let legacy = read(&base.replacen(faims, &faims.replace("UO:0000218", "UO:000218"), 1));
    assert_eq!(legacy.unwrap().spectra[2].drift_time, -35.0);
    let wrong = base.replacen(
        millisecond,
        &millisecond.replace("UO:0000028", "UO:000218"),
        1,
    );
    assert!(matches!(read(&wrong), Err(Error::Unsupported(_))));
    for bad in ["NaN", "inf", "", "-3x"] {
        let changed = base.replacen(faims, &faims.replace("\"-35\"", &format!("\"{bad}\"")), 1);
        assert!(read(&changed).is_err(), "value {bad:?} accepted");
    }
    // Spectrum-level then scan-level values: identical is one value, a
    // different value or unit has no single native representation.
    let conflict = include_str!("data/mzml_mobility/scan_mobility_conflict.mzML").replacen(
        "value=\"8.5\"",
        "value=\"7.5\"",
        1,
    );
    assert!(
        matches!(read(&conflict), Err(Error::Unsupported(text)) if text.contains("conflicting"))
    );
    let identical = conflict.replacen("value=\"-60\"", "value=\"-50\"", 1);
    let experiment = read(&identical).unwrap();
    assert_eq!(
        (
            experiment.spectra[0].drift_time,
            experiment.spectra[1].drift_time
        ),
        (7.5, -50.0)
    );
    // The same value under another unit conflicts too.
    let scan_term = r#"accession="MS:1001581" name="FAIMS compensation voltage" value="-50" unitCvRef="UO" unitAccession="UO:0000218""#;
    let last = identical.rfind(scan_term).expect("scan-level voltage");
    let mut changed = identical.clone();
    changed.replace_range(
        last..last + scan_term.len(),
        r#"accession="MS:1002476" name="ion mobility drift time" value="-50" unitCvRef="UO" unitAccession="UO:0000028""#,
    );
    assert!(
        matches!(read(&changed), Err(Error::Unsupported(text)) if text.contains("conflicting"))
    );
}

#[test]
fn mobility_array_units_are_validated_and_kept_as_source_metadata() {
    let base = include_str!("data/mzml_mobility/im_arrays.mzML");
    let unit_term = r#"accession="MS:1002816" name="mean ion mobility array" unitCvRef="UO" unitAccession="UO:0000028" unitName="millisecond""#;
    assert!(base.contains(unit_term));
    for (bad, parse_error) in [
        (
            r#"accession="MS:1002816" name="mean ion mobility array" unitCvRef="UO" unitName="millisecond""#,
            true,
        ),
        (
            r#"accession="MS:1002816" name="mean ion mobility array" unitCvRef="MS" unitAccession="UO:0000028" unitName="millisecond""#,
            true,
        ),
        (
            r#"accession="MS:1002816" name="mean ion mobility array" unitCvRef="XX" unitAccession="XX:0000028" unitName="millisecond""#,
            false,
        ),
    ] {
        let error = read(&base.replacen(unit_term, bad, 1)).unwrap_err();
        assert_eq!(
            matches!(error, Error::Parse { .. }),
            parse_error,
            "{bad}: {error}"
        );
    }
    // A cvParam unit and a unit_accession userParam cannot both own the array.
    let duplicate = base.replacen(
        unit_term,
        &format!("{unit_term}/><userParam name=\"unit_accession\" value=\"UO:0000028\""),
        1,
    );
    assert!(read(&duplicate).is_err());
}

fn spectrum(rt: f64, drift_time: f64, unit: DriftTimeUnit) -> MSSpectrum {
    MSSpectrum {
        native_id: format!("scan={}", rt as i64),
        rt,
        peaks: vec![Peak1D::new(100.0, 10.0), Peak1D::new(101.0, 20.0)],
        drift_time,
        drift_time_unit: unit,
        ..Default::default()
    }
}

#[test]
fn spectrum_mobility_round_trips_through_the_scan_writer() {
    use DriftTimeUnit as D;
    let experiment = MSExperiment {
        spectra: vec![
            spectrum(1.0, 7.5, D::Millisecond),
            spectrum(2.0, 0.85, D::InverseReducedMobility),
            spectrum(3.0, -35.0, D::FaimsCompensationVoltage),
            spectrum(4.0, 350.5, D::CollisionCrossSection),
            spectrum(5.0, -1.0, D::FaimsCompensationVoltage),
            spectrum(6.0, -1.0, D::None),
            // Without a retention time the mobility alone still needs a scan.
            spectrum(-1.0, 3.25, D::Millisecond),
        ],
        ..Default::default()
    };
    for zlib_compression in [false, true] {
        let mut bytes = Vec::new();
        mzml::write_with_options(
            &mut bytes,
            &experiment,
            &mzml::WriteOptions { zlib_compression },
        )
        .unwrap();
        assert_eq!(mzml::read(bytes.as_slice()).unwrap(), experiment);
    }
    let mut bytes = Vec::new();
    mzml::write(&mut bytes, &experiment).unwrap();
    let xml = String::from_utf8(bytes).unwrap();
    // Unit spelling of the pinned writer, MzMLHandler.cpp:5414-5415.
    assert!(xml.contains(r#"<cvParam cvRef="MS" accession="MS:1001581" name="FAIMS compensation voltage" value="-35" unitCvRef="UO" unitAccession="UO:0000218" unitName="volt"/>"#));
    for name in [
        "scan_mobility.mzML",
        "im_arrays.mzML",
        "spectrum_representation.mzML",
        "FAIMS_test_data.mzML",
    ] {
        let loaded = load(name).unwrap();
        let mut bytes = Vec::new();
        mzml::write(&mut bytes, &loaded).unwrap();
        let reread = mzml::read(bytes.as_slice()).unwrap();
        assert_eq!(reread.spectra, loaded.spectra, "{name}");
    }
}

#[test]
fn writer_refuses_lossy_spectrum_mobility_before_any_output() {
    use DriftTimeUnit as D;
    for (drift_time, unit) in [
        (1.0, D::None),
        (-1.0, D::Millisecond),
        (-1.0, D::InverseReducedMobility),
        (-1.0, D::CollisionCrossSection),
        (f64::NAN, D::Millisecond),
        (f64::INFINITY, D::FaimsCompensationVoltage),
    ] {
        let experiment = MSExperiment {
            spectra: vec![spectrum(1.0, drift_time, unit)],
            ..Default::default()
        };
        let mut output = vec![42];
        assert!(
            mzml::write(&mut output, &experiment).is_err(),
            "{drift_time} {unit:?} accepted"
        );
        assert_eq!(output, [42]);
    }
}
