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
    assert_eq!(spectrum.metadata["label"].as_str().unwrap(), "A & B");
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
    assert_eq!(
        chromatogram.metadata["label"].as_str().unwrap(),
        "chromatogram"
    );
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
    // A declared record or binary-array count disagreeing with the actual
    // children is advisory on reading. MzMLHandler.cpp reads `count` only for a
    // progress range and `reserveSpaceSpectra` (:965-979), the chromatogram
    // equivalent (:996-1013), `bin_data_.reserve(...)` (:1015-1017) and the
    // selectedIon warning (:1374); it is compared against nothing. The upstream
    // class-test fixture `MzMLFile_1.mzML` declares
    // `<binaryDataArrayList count="2">` with four arrays and C++ loads it.
    let expected = parse(&source).unwrap();
    for xml in [
        source.replacen("spectrumList count=\"1\"", "spectrumList count=\"2\"", 1),
        source.replacen(
            "chromatogramList count=\"1\"",
            "chromatogramList count=\"9\"",
            1,
        ),
        source.replacen(
            "binaryDataArrayList count=\"2\"",
            "binaryDataArrayList count=\"3\"",
            1,
        ),
    ] {
        assert_eq!(parse(&xml).unwrap(), expected);
    }
    // The attribute itself stays required and numeric at both sites.
    for xml in [
        source.replacen("spectrumList count=\"1\"", "spectrumList", 1),
        source.replacen("spectrumList count=\"1\"", "spectrumList count=\"many\"", 1),
        source.replacen("binaryDataArrayList count=\"2\"", "binaryDataArrayList", 1),
        source.replacen(
            "binaryDataArrayList count=\"2\"",
            "binaryDataArrayList count=\"-1\"",
            1,
        ),
    ] {
        assert!(parse(&xml).is_err());
    }
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

// `Record::finish` skips the per-peak value loop of `MSSpectrum::validate` and
// `MSChromatogram::validate`, because every coordinate it hands them has been
// refused by the binary decoder if nonfinite and every intensity has been
// refused by the f32 narrowing. This pins that premise where it is made, in
// every position of both primary arrays of both record kinds: the refusal must
// carry the reader's own message, never the kernel validator's. If a decode
// path ever stops rejecting nonfinite values, this fails here rather than the
// skipped check silently failing to catch it downstream.
#[test]
fn nonfinite_peak_values_are_refused_by_the_decoder_not_the_validator() {
    let mut experiment = MSExperiment::new();
    let mut spectrum = MSSpectrum::from_peaks(
        (0..5)
            .map(|i| Peak1D::new(100.0 + f64::from(i), 1.0 + i as f32))
            .collect(),
    );
    spectrum.rt = 1.0;
    spectrum.ms_level = 1;
    experiment.spectra.push(spectrum);
    experiment.chromatograms.push(MSChromatogram::from_peaks(
        (0..5)
            .map(|i| ChromatogramPeak::new(f64::from(i), 1.0 + i as f32))
            .collect(),
    ));
    let xml = encoded_xml(&experiment, false);
    // The writer emits, in this order: spectrum m/z as Float64, spectrum
    // intensity as Float32, chromatogram time as Float64, chromatogram
    // intensity as Float32.
    for (array, width) in [(0_usize, 8_usize), (1, 4), (2, 8), (3, 4)] {
        for position in [0_usize, 2, 4] {
            for pattern in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let broken = mutate_payload(&xml, array, |bytes| {
                    let at = position * width;
                    if width == 8 {
                        bytes[at..at + 8].copy_from_slice(&pattern.to_le_bytes());
                    } else {
                        bytes[at..at + 4].copy_from_slice(&(pattern as f32).to_le_bytes());
                    }
                });
                let reported = parse(&broken).unwrap_err().to_string();
                assert!(
                    reported.contains("nonfinite binary value"),
                    "array {array} position {position} pattern {pattern}: {reported}"
                );
                // The kernel validator's wording. Seeing it here would mean the
                // decoder had let the value through.
                assert!(
                    !reported.contains("must be finite"),
                    "array {array} position {position} pattern {pattern}: {reported}"
                );
            }
        }
    }
    // The one nonfinite peak intensity the decoder cannot see is a finite f64
    // with no finite f32. That is the second half of the premise, and
    // `rejects_nonfinite_and_overflowing_peak_values` above pins it: the
    // narrowing refuses it with `intensity overflows f32`, again before the
    // record reaches the validator.
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
fn rejects_duplicate_or_misplaced_precursor_and_scan_lists_despite_advisory_counts() {
    let xml = encoded_xml(&sample(), false);
    let expected = parse(&xml).unwrap();
    // `precursorList` and `scanWindowList` have no open-tag handler in
    // MzMLHandler.cpp, `scanList` is only read as a parent tag (:2259, :3496),
    // and `selectedIonList` only warns when its count exceeds one (:1371-1375).
    // None of them compares the declared count with the children, so a
    // disagreement is advisory on reading.
    for changed in [
        xml.replacen("precursorList count=\"1\"", "precursorList count=\"2\"", 1),
        xml.replacen(
            "selectedIonList count=\"1\"",
            "selectedIonList count=\"2\"",
            1,
        ),
        xml.replacen("scanList count=\"1\"", "scanList count=\"2\"", 1),
    ] {
        assert_eq!(parse(&changed).unwrap(), expected);
    }
    // Each count stays required and numeric, and duplicate or misplaced lists
    // remain structural errors.
    for changed in [
        xml.replacen("precursorList count=\"1\"", "precursorList", 1),
        xml.replacen(
            "selectedIonList count=\"1\"",
            "selectedIonList count=\"\"",
            1,
        ),
        xml.replacen("scanList count=\"1\"", "scanList count=\"two\"", 1),
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

// The reader tests text nodes and attribute values against the XML 1.0 `Char`
// production on their UTF-8 bytes, a fixed-size block at a time. A block that
// holds nothing but ordinary printable ASCII is skipped whole, so every case
// below places the character under test at a series of offsets around that
// block length: at 0, inside the first block, exactly on its boundary and past
// it, so no block length can hide a character from the scan.
#[test]
fn xml_character_production_holds_at_every_offset_in_a_scanned_text() {
    let source = encoded_xml(&sample(), false);
    assert_eq!(source.matches("integration").count(), 1);
    for pad in [0usize, 1, 2, 31, 63, 64, 65, 127, 129] {
        let filler = "z".repeat(pad);
        // Admitted: the last scalar below the surrogate block, the first above
        // it, U+FFFD, and a supplementary-plane scalar. Each has to survive the
        // round trip, not merely be accepted.
        for accepted in ['\u{d7ff}', '\u{e000}', '\u{fffd}', '\u{10000}'] {
            let value = format!("{filler}{accepted}{filler}");
            let xml = source.replacen("integration", &value, 1);
            let mut expected = sample();
            expected.chromatograms[0]
                .metadata
                .insert("source".into(), value.as_str().into());
            assert_eq!(
                parse(&xml).unwrap_or_else(|e| panic!("pad {pad} {accepted:?}: {e}")),
                expected
            );
        }
        // Tab and carriage return are admitted too, but an XML parser is free
        // to normalise them inside an attribute value, so only acceptance is
        // asserted for those.
        for accepted in ['\u{9}', '\u{d}'] {
            let xml = source.replacen("integration", &format!("{filler}{accepted}{filler}"), 1);
            assert!(parse(&xml).is_ok(), "pad {pad} refused {accepted:?}");
        }
        // Refused: the C0 controls other than tab, newline and return, and the
        // two noncharacters at the end of the basic plane.
        for refused in [
            '\u{0}', '\u{1}', '\u{b}', '\u{c}', '\u{1f}', '\u{fffe}', '\u{ffff}',
        ] {
            let xml = source.replacen("integration", &format!("{filler}{refused}{filler}"), 1);
            assert!(
                matches!(parse(&xml), Err(Error::InvalidValue(message))
                    if message.contains("invalid XML 1.0 characters")),
                "pad {pad} accepted {refused:?}"
            );
        }
    }
}

// Base64 text is accumulated in runs of ordinary characters rather than one
// character at a time, and the three ways a `<binary>` node can be malformed
// keep the order they had per character: an invalid XML character anywhere in
// the node outranks a non-ASCII base64 character, which outranks the declared
// `encodedLength` being exceeded by a later character.
#[test]
fn binary_text_keeps_its_whitespace_rule_and_its_error_order() {
    let source = encoded_xml(&sample(), false);
    let payload_start = source.find("<binary>").unwrap() + "<binary>".len();
    let payload_end = payload_start + source[payload_start..].find("</binary>").unwrap();
    let payload = source[payload_start..payload_end].to_owned();
    let expected = parse(&source).unwrap();
    // XML whitespace inside the payload is stripped, at any offset, and in runs.
    for split in 0..=payload.len() {
        for gap in [" ", "\n", "\r\n", "\t \r\n\t"] {
            let spaced = format!("{}{gap}{}", &payload[..split], &payload[split..]);
            let xml = source.replacen(&payload, &spaced, 1);
            assert_eq!(parse(&xml).unwrap(), expected, "split {split} gap {gap:?}");
        }
    }
    // An invalid XML character is rejected as such, before the payload is read.
    let xml = source.replacen(&payload, &format!("\u{1}{payload}"), 1);
    assert!(matches!(parse(&xml), Err(Error::InvalidValue(message))
        if message.contains("invalid XML 1.0 characters")));
    // A non-ASCII character that IS a valid XML character reaches the payload.
    let xml = source.replacen(&payload, &format!("{payload}\u{e9}"), 1);
    assert!(matches!(parse(&xml), Err(Error::Parse { message, .. })
        if message == "non-ASCII base64 text"));
    // One extra base64 character puts the node past its declared length. The
    // non-ASCII character after it is never reached, which pins the order.
    let xml = source.replacen(&payload, &format!("{payload}A\u{e9}"), 1);
    assert!(matches!(parse(&xml), Err(Error::Parse { message, .. })
        if message == "binary text exceeds encodedLength"));
}

// Binary arrays lend their text and decoded bytes from buffers the reader keeps
// across records, so a long array is followed by shorter ones whose payloads
// must not pick up anything the long one left behind, in either direction.
#[test]
fn arrays_of_changing_length_reuse_the_reader_buffers_without_bleeding() {
    let lengths = [512usize, 3, 257, 1, 64, 300, 0, 129];
    let mut experiment = MSExperiment::new();
    for (index, length) in lengths.into_iter().enumerate() {
        let mut spectrum = MSSpectrum::from_peaks(
            (0..length)
                .map(|i| Peak1D::new(100.0 + i as f64, (index + i) as f32))
                .collect(),
        );
        spectrum.native_id = format!("controllerType=0 controllerNumber=1 scan={}", index + 1);
        spectrum.rt = index as f64;
        spectrum.ms_level = 1;
        experiment.spectra.push(spectrum);
    }
    for zlib in [false, true] {
        let xml = encoded_xml(&experiment, zlib);
        assert_eq!(parse(&xml).unwrap(), experiment, "zlib {zlib}");
    }
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
