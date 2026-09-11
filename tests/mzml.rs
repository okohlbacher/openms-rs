// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

#![cfg(feature = "mzml")]

use base64::{Engine, engine::general_purpose::STANDARD};
use openms::format::mzml::{self, ReadOptions, WriteOptions};
use openms::kernel::{DataArray, SpectrumType};
use openms::{
    ChromatogramPeak, Error, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, Precursor,
};
use std::io::Cursor;

const INDEPENDENT: &str = include_str!("data/mzml_independent.mzML");
const MINIMAL: &str = include_str!("data/mzml_upstream_minimal.mzML");
const SERUM: &str = include_str!("data/mzml_upstream_serum.mzML");

fn parse(xml: &str) -> openms::Result<MSExperiment> {
    mzml::read(Cursor::new(xml.as_bytes()))
}

fn sample() -> MSExperiment {
    let mut experiment = MSExperiment::new();
    experiment
        .settings
        .metadata
        .insert("sample".into(), "A & B \"quoted\" <values>\nUTF-8 µ".into());
    let mut spectrum =
        MSSpectrum::from_peaks(vec![Peak1D::new(100.25, 5.0), Peak1D::new(200.5, -2.0)]);
    spectrum.native_id = "controllerType=0 controllerNumber=1 scan=7".into();
    spectrum.name = "spectrum & name".into();
    spectrum.rt = 30.25;
    spectrum.ms_level = 2;
    spectrum.spectrum_type = SpectrumType::Centroid;
    spectrum.precursors.push(Precursor {
        mz: 500.25,
        charge: 2,
        intensity: 234.5,
        ..Precursor::default()
    });
    spectrum
        .metadata
        .insert("note".into(), "one\ttwo\rthree".into());
    experiment.spectra.push(spectrum);
    let mut chromatogram = MSChromatogram::from_peaks(vec![
        ChromatogramPeak::new(30.0, 10.0),
        ChromatogramPeak::new(60.0, 20.0),
    ]);
    chromatogram.native_id = "TIC".into();
    chromatogram.name = "signal".into();
    chromatogram.precursor = Precursor {
        mz: 321.5,
        charge: -2,
        intensity: 10.0,
        ..Precursor::default()
    };
    chromatogram
        .metadata
        .insert("source".into(), "integration".into());
    experiment.chromatograms.push(chromatogram);
    experiment
}

fn encoded_xml(experiment: &MSExperiment, zlib: bool) -> String {
    let mut bytes = Vec::new();
    mzml::write_with_options(
        &mut bytes,
        experiment,
        &WriteOptions {
            zlib_compression: zlib,
        },
    )
    .unwrap();
    String::from_utf8(bytes).unwrap()
}

#[test]
fn reads_pinned_upstream_minimal_and_indexed_fixture() {
    assert_eq!(parse(MINIMAL).unwrap(), MSExperiment::new());
    // The original ASCII bytes and Latin-1 declaration are accepted unchanged.
    assert!(SERUM.is_ascii());
    let experiment = parse(SERUM).unwrap();
    assert_eq!(experiment.spectra.len(), 1);
    let spectrum = &experiment.spectra[0];
    assert_eq!(spectrum.native_id, "spectrum=0");
    assert_eq!(spectrum.spectrum_type, SpectrumType::Centroid);
    assert_eq!(spectrum.rt, -1.0);
    let expected = [
        Peak1D::new(109.9828037554159, 7452.0927734375_f64 as f32),
        Peak1D::new(109.99478129476843, 1467.237548828125_f64 as f32),
        Peak1D::new(109.99813009042244, 1007.3275756835938_f64 as f32),
    ];
    assert_eq!(spectrum.peaks, expected);
}

#[test]
fn reads_independently_encoded_mixed_precision_compression_and_minutes() {
    let experiment = parse(INDEPENDENT).unwrap();
    let spectrum = &experiment.spectra[0];
    assert_eq!(
        spectrum.peaks,
        [
            Peak1D::new(100.25, 10.5),
            Peak1D::new(200.5, 0.0),
            Peak1D::new(300.75, 35.25)
        ]
    );
    assert_eq!(spectrum.rt, 90.0);
    assert_eq!(spectrum.ms_level, 2);
    assert_eq!(spectrum.spectrum_type, SpectrumType::Centroid);
    assert_eq!(
        spectrum.precursors,
        [Precursor {
            mz: 500.25,
            charge: 2,
            intensity: 123.5,
            ..Precursor::default()
        }]
    );
    assert_eq!(spectrum.metadata["label"], "A & B");
    let chromatogram = &experiment.chromatograms[0];
    assert_eq!(
        chromatogram.peaks,
        [
            ChromatogramPeak::new(30.0, 1.0),
            ChromatogramPeak::new(60.0, 7.0),
            ChromatogramPeak::new(90.0, 2.0)
        ]
    );
    assert_eq!(chromatogram.native_id, "TIC");
    assert_eq!(chromatogram.metadata["label"], "chromatogram");
}

#[test]
fn round_trips_supported_fields_with_both_compression_modes() {
    for zlib in [false, true] {
        let experiment = sample();
        assert_eq!(parse(&encoded_xml(&experiment, zlib)).unwrap(), experiment);
        assert_eq!(
            parse(&encoded_xml(&MSExperiment::new(), zlib)).unwrap(),
            MSExperiment::new()
        );
        let mut empty = MSExperiment::new();
        empty.spectra.push(MSSpectrum {
            native_id: "scan=1".into(),
            ..Default::default()
        });
        empty.chromatograms.push(MSChromatogram {
            native_id: "TIC".into(),
            ..Default::default()
        });
        assert_eq!(parse(&encoded_xml(&empty, zlib)).unwrap(), empty);
    }
}

#[test]
fn accepts_namespace_prefixes_and_whitespace_in_base64() {
    let xml = encoded_xml(&sample(), false);
    let prefixed = xml
        .replace(
            "xmlns=\"http://psi.hupo.org/ms/mzml\"",
            "xmlns:m=\"http://psi.hupo.org/ms/mzml\"",
        )
        .replace("<mzML ", "<m:mzML ")
        .replace("</mzML>", "</m:mzML>");
    // Other elements retain the namespace through their original default binding.
    let prefixed = prefixed.replace("<m:mzML ", "<m:mzML xmlns=\"http://psi.hupo.org/ms/mzml\" ");
    assert_eq!(parse(&prefixed).unwrap(), sample());
    let spaced = xml
        .replace("<binary>", "<binary>\n ")
        .replace("</binary>", " \n</binary>");
    assert_eq!(parse(&spaced).unwrap(), sample());
}

#[test]
fn rejects_malformed_structure_and_lengths() {
    let source = encoded_xml(&sample(), false);
    let cases = [
        source
            .replace("<mzML ", "<wrong ")
            .replace("</mzML>", "</wrong>"),
        source.replace("http://psi.hupo.org/ms/mzml", "urn:wrong"),
        source.replace("</run></mzML>", "</run>"),
        source.replace("</spectrum>", "</chromatogram>"),
        source.clone() + &source,
        "unwanted text".to_owned() + &source,
        source.replacen("defaultArrayLength=\"2\"", "defaultArrayLength=\"3\"", 1),
        source.replacen("spectrumList count=\"1\"", "spectrumList count=\"2\"", 1),
        source.replacen(
            "binaryDataArrayList count=\"2\"",
            "binaryDataArrayList count=\"3\"",
            1,
        ),
        source.replacen("MS:1000514", "MS:1000595", 1),
        source.replacen("encodedLength=\"24\"", "encodedLength=\"23\"", 1),
        source.replacen("<binary>", "<binary>!", 1),
        source.replacen("<binary>", "<binary/><binary>", 1),
        source.replacen("MS:1000523", "MS:1000521", 1),
        source.replacen("value=\"2\"", "value=\"0\"", 1),
        source.replacen("value=\"30.25\"", "value=\"NaN\"", 1),
        source.replacen("UO:0000010", "UO:0000028", 1),
    ];
    for (i, xml) in cases.into_iter().enumerate() {
        assert!(parse(&xml).is_err(), "invalid case {i} accepted");
    }
    let start = source.find("<binaryDataArrayList").unwrap();
    let end = start
        + source[start..].find("</binaryDataArrayList>").unwrap()
        + "</binaryDataArrayList>".len();
    let missing = format!("{}{}", &source[..start], &source[end..]);
    assert!(parse(&missing).is_err());
}

fn mutate_first_payload(xml: &str, mutation: impl FnOnce(&mut Vec<u8>)) -> String {
    mutate_payload(xml, 0, mutation)
}

fn mutate_payload(xml: &str, index: usize, mutation: impl FnOnce(&mut Vec<u8>)) -> String {
    let start = xml.match_indices("<binary>").nth(index).unwrap().0 + "<binary>".len();
    let end = start + xml[start..].find("</binary>").unwrap();
    let mut bytes = STANDARD.decode(&xml[start..end]).unwrap();
    mutation(&mut bytes);
    let encoded = STANDARD.encode(bytes);
    let array_start = xml[..start]
        .rfind("<binaryDataArray encodedLength=\"")
        .unwrap();
    let length_start = array_start + "<binaryDataArray encodedLength=\"".len();
    let length_end = length_start + xml[length_start..].find('"').unwrap();
    format!(
        "{}{}{}{}{}",
        &xml[..length_start],
        encoded.len(),
        &xml[length_end..start],
        encoded,
        &xml[end..]
    )
}

#[test]
fn rejects_truncated_corrupt_or_excess_compressed_data() {
    let xml = encoded_xml(&sample(), true);
    assert!(
        parse(&mutate_first_payload(&xml, |b| {
            b.truncate(b.len() - 1);
        }))
        .is_err()
    );
    assert!(
        parse(&mutate_first_payload(&xml, |b| {
            *b.last_mut().unwrap() ^= 1;
        }))
        .is_err()
    );
    assert!(parse(&mutate_first_payload(&xml, |b| b.extend([0, 1]))).is_err());
    // Valid compression stream, deliberately incorrect declared decoded size.
    assert!(
        parse(&xml.replacen("defaultArrayLength=\"2\"", "defaultArrayLength=\"1\"", 1)).is_err()
    );
}

#[test]
fn rejects_nonfinite_and_overflowing_peak_values() {
    let xml = encoded_xml(&sample(), false);
    assert!(
        parse(&mutate_first_payload(&xml, |b| b[..8].copy_from_slice(&f64::NAN.to_le_bytes())))
            .is_err()
    );
    let huge = mutate_payload(INDEPENDENT, 1, |bytes| {
        use std::io::Write;
        let mut encoder =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        for _ in 0..3 {
            encoder.write_all(&1e40_f64.to_le_bytes()).unwrap();
        }
        *bytes = encoder.finish().unwrap();
    });
    assert!(
        parse(&huge)
            .unwrap_err()
            .to_string()
            .contains("overflows f32")
    );
}

#[test]
fn rejects_unsupported_encodings_and_array_types() {
    let xml = encoded_xml(&sample(), false);
    for changed in [
        xml.replacen("MS:1000523", "MS:1000522", 1),
        xml.replacen("MS:1000514", "MS:1000516", 1),
        xml.replace("encoding=\"UTF-8\"", "encoding=\"UTF-16\""),
        xml.replace("version=\"1.1.0\"", "version=\"1.0.0\""),
    ] {
        assert!(matches!(parse(&changed), Err(Error::Unsupported(_))));
    }
    // Numpress is now supported, but relabeling ordinary bytes is malformed.
    assert!(parse(&xml.replacen("MS:1000576", "MS:1002312", 1)).is_err());
    let dtd = xml.replace("<mzML ", "<!DOCTYPE mzML [<!ENTITY x 'test'>]><mzML ");
    assert!(matches!(parse(&dtd), Err(Error::Unsupported(_))));
}

#[test]
fn enforces_configured_input_and_output_limits() {
    for options in [
        ReadOptions {
            max_xml_bytes: 128,
            ..Default::default()
        },
        ReadOptions {
            max_array_bytes: 8,
            ..Default::default()
        },
        ReadOptions {
            max_total_peaks: 2,
            ..Default::default()
        },
        ReadOptions {
            max_records: 1,
            ..Default::default()
        },
    ] {
        assert!(mzml::read_with_options(Cursor::new(INDEPENDENT), &options).is_err());
    }
    let options = ReadOptions {
        max_xml_bytes: INDEPENDENT.len() as u64,
        ..Default::default()
    };
    assert!(mzml::read_with_options(Cursor::new(INDEPENDENT), &options).is_ok());
    let too_short = ReadOptions {
        max_xml_bytes: INDEPENDENT.len() as u64 - 1,
        ..Default::default()
    };
    assert!(mzml::read_with_options(Cursor::new(INDEPENDENT), &too_short).is_err());
}

#[test]
fn writing_rejects_data_loss_invalid_ids_and_invalid_xml_before_output() {
    let mut extra = sample();
    extra.spectra[0]
        .float_data_arrays
        .push(DataArray::new("quality", vec![1.0, f32::NAN]));
    let mut invalid_id = sample();
    invalid_id.spectra[0].native_id = "invalid id".into();
    let mut duplicate = sample();
    duplicate.spectra.push(duplicate.spectra[0].clone());
    let mut invalid_text = sample();
    invalid_text
        .settings
        .metadata
        .insert("invalid".into(), "\0".into());
    let mut reserved = sample();
    reserved
        .settings
        .metadata
        .insert("openms-rust:name".into(), "reserved".into());
    let mut bad_value = sample();
    bad_value.spectra[0].peaks[0].mz = f64::INFINITY;
    for experiment in [
        extra,
        invalid_id,
        duplicate,
        invalid_text,
        reserved,
        bad_value,
    ] {
        let mut output = Vec::new();
        assert!(mzml::write(&mut output, &experiment).is_err());
        assert!(output.is_empty());
    }
}

#[test]
fn validates_writer_output_against_pinned_schema_when_xmllint_is_available() {
    if std::process::Command::new("xmllint")
        .arg("--version")
        .output()
        .is_err()
    {
        eprintln!("xmllint unavailable; independent XSD validation is not executed on this host");
        return;
    }
    let schema = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/mzml_1_10.xsd");
    for (index, experiment) in [sample(), MSExperiment::new()].into_iter().enumerate() {
        let path = std::env::temp_dir().join(format!(
            "openms-rust-mzml-{}-{index}.mzML",
            std::process::id()
        ));
        std::fs::write(&path, encoded_xml(&experiment, index == 0)).unwrap();
        let output = std::process::Command::new("xmllint")
            .args(["--nonet", "--noout", "--schema"])
            .arg(&schema)
            .arg(&path)
            .output()
            .unwrap();
        std::fs::remove_file(path).unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
}

#[test]
fn rejects_duplicate_or_inconsistent_precursor_and_scan_lists() {
    let xml = encoded_xml(&sample(), false);
    for changed in [
        xml.replacen("precursorList count=\"1\"", "precursorList count=\"2\"", 1),
        xml.replacen(
            "selectedIonList count=\"1\"",
            "selectedIonList count=\"2\"",
            1,
        ),
        xml.replacen("scanList count=\"1\"", "scanList count=\"2\"", 1),
        xml.replacen("</scanList>", "</scanList><scanList count=\"0\"/>", 1),
        xml.replacen(
            "</precursorList>",
            "</precursorList><precursorList count=\"0\"/>",
            1,
        ),
        xml.replacen(
            "</selectedIonList>",
            "</selectedIonList><selectedIonList count=\"0\"/>",
            1,
        ),
        xml.replacen("<scan>", "<scan><precursorList count=\"0\"/>", 1),
    ] {
        assert!(parse(&changed).is_err());
    }
}

#[test]
fn rejects_malformed_comments_and_non_ascii_bytes_under_ascii_declaration() {
    let xml = encoded_xml(&sample(), false);
    assert!(parse(&xml.replace("</run>", "<!-- invalid -- comment --></run>")).is_err());
    assert!(parse(&xml.replace("encoding=\"UTF-8\"", "encoding=\"US-ASCII\"")).is_err());
    assert!(parse(&MINIMAL.replace("encoding=\"UTF-8\"", "encoding=\"US-ASCII\"")).is_ok());
}

#[test]
fn writer_propagates_buffer_flush_errors() {
    struct FlushFailure;
    impl std::io::Write for FlushFailure {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Err(std::io::Error::other("flush failed"))
        }
    }
    assert!(matches!(
        mzml::write(FlushFailure, &sample()),
        Err(Error::Io(_))
    ));
}

#[test]
fn latin1_declaration_accepts_ascii_and_character_references_without_misdecoding_bytes() {
    let latin = MINIMAL.replace("UTF-8", "iSo-8859-1");
    assert!(parse(&latin).is_ok());
    // ASCII character references are decoded by XML independently of byte encoding.
    let escaped = latin.replace("<run ", "<run harmless=\"caf&#233;\" ");
    assert!(parse(&escaped).is_ok());
    for text in ["<!-- café -->", "<?example café?>"] {
        let xml = latin.replace("<mzML ", &format!("{text}<mzML "));
        assert!(matches!(parse(&xml), Err(Error::Unsupported(_))));
        let mut bytes = xml.as_bytes().to_vec();
        let index = bytes.windows(2).position(|s| s == [0xc3, 0xa9]).unwrap();
        bytes.splice(index..index + 2, [0xe9]);
        assert!(matches!(
            mzml::read(Cursor::new(bytes)),
            Err(Error::Unsupported(_))
        ));
    }
}
