// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

#![cfg(feature = "mzml")]

//! Independent writer boundary tests. Binary layout follows the pinned
//! MzMLHandler.cpp integer/string branches (source revision 7c029e8).
//! Expected bytes are constructed directly, without the native mzML reader.

use base64::{Engine, engine::general_purpose::STANDARD};
use flate2::read::ZlibDecoder;
use openms::format::mzml::{self, WriteOptions};
use openms::kernel::DataArray;
use openms::{ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D};
use std::io::{Cursor, Read, Write};

fn sample() -> MSExperiment {
    let mut spectrum =
        MSSpectrum::from_peaks(vec![Peak1D::new(100.0, 1.0), Peak1D::new(200.0, 2.0)]);
    spectrum.native_id = "scan=1".into();
    spectrum.float_data_arrays = vec![DataArray::new("float", vec![1.0, -2.0])];
    let mut chromatogram = MSChromatogram::from_peaks(vec![
        ChromatogramPeak::new(1.0, 2.0),
        ChromatogramPeak::new(2.0, 3.0),
    ]);
    chromatogram.native_id = "TIC".into();
    chromatogram.float_data_arrays = vec![DataArray::new("float", vec![1.0, -2.0])];
    chromatogram.integer_data_arrays = vec![DataArray::new("integer", vec![-1, 2])];
    chromatogram.string_data_arrays = vec![DataArray::new(
        "string",
        vec!["left".into(), "right".into()],
    )];
    MSExperiment {
        spectra: vec![spectrum],
        chromatograms: vec![chromatogram],
        ..Default::default()
    }
}

fn encoded(experiment: &MSExperiment, compressed: bool) -> Vec<u8> {
    let mut bytes = Vec::new();
    mzml::write_with_options(
        &mut bytes,
        experiment,
        &WriteOptions {
            zlib_compression: compressed,
        },
    )
    .unwrap();
    bytes
}

#[derive(Default)]
struct CountingWriter {
    writes: usize,
    flushes: usize,
}
impl Write for CountingWriter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.writes += 1;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.flushes += 1;
        Ok(())
    }
}

#[test]
fn every_late_auxiliary_validation_error_precedes_all_output() {
    // Each mutation is in the last record, after a valid spectrum and headers.
    // This proves whole-document preflight instead of validation during output.
    type Mutation = fn(&mut MSExperiment);
    let mutations: &[Mutation] = &[
        |e| e.chromatograms[0].float_data_arrays[0].name.clear(),
        |e| e.chromatograms[0].integer_data_arrays[0].name.clear(),
        |e| e.chromatograms[0].string_data_arrays[0].name.clear(),
        |e| e.chromatograms[0].float_data_arrays[0].name = "bad\u{1}name".into(),
        |e| e.chromatograms[0].integer_data_arrays[0].name = "bad\u{b}name".into(),
        |e| e.chromatograms[0].string_data_arrays[0].name = "bad\u{0}name".into(),
        |e| {
            e.chromatograms[0]
                .metadata
                .insert("bad\u{1}".into(), "value".into());
        },
        |e| {
            e.chromatograms[0]
                .metadata
                .insert("key".into(), "bad\u{b}".into());
        },
        |e| e.chromatograms[0].integer_data_arrays[0].name = "float".into(),
        |e| e.chromatograms[0].string_data_arrays[0].name = "integer".into(),
        |e| {
            e.chromatograms[0]
                .float_data_arrays
                .push(DataArray::new("float", vec![]))
        },
        |e| {
            e.chromatograms[0]
                .integer_data_arrays
                .push(DataArray::new("integer", vec![]))
        },
        |e| {
            e.chromatograms[0]
                .string_data_arrays
                .push(DataArray::new("string", vec![]))
        },
        |e| e.chromatograms[0].float_data_arrays[0].data[1] = f32::NAN,
        |e| e.chromatograms[0].float_data_arrays[0].data[1] = f32::INFINITY,
        |e| e.chromatograms[0].float_data_arrays[0].data[1] = f32::NEG_INFINITY,
        |e| e.chromatograms[0].string_data_arrays[0].data[1] = "embedded\0NUL".into(),
        |e| e.chromatograms[0].string_data_arrays[0].data[1] = "non-ASCII µ".into(),
        |e| e.chromatograms[0].float_data_arrays[0].data.push(1.0),
        |e| {
            e.chromatograms[0].integer_data_arrays[0].data.pop();
        },
        |e| {
            e.chromatograms[0].string_data_arrays[0]
                .data
                .push("extra".into())
        },
    ];
    for (case, change) in mutations.iter().enumerate() {
        for compressed in [false, true] {
            let mut input = sample();
            change(&mut input);
            let mut writer = CountingWriter::default();
            let result = mzml::write_with_options(
                &mut writer,
                &input,
                &WriteOptions {
                    zlib_compression: compressed,
                },
            );
            assert!(
                result.is_err(),
                "accepted invalid case {case}, compressed={compressed}"
            );
            assert_eq!(
                writer.writes, 0,
                "case {case} emitted output before validation"
            );
            assert_eq!(writer.flushes, 0, "case {case} flushed before validation");
        }
    }
}

fn array<'a>(xml: &'a str, name: &str) -> &'a str {
    let name = format!("value=\"{name}\"");
    xml.split("<binaryDataArray ")
        .skip(1)
        .map(|block| block.split_once("</binaryDataArray>").unwrap().0)
        .find(|block| block.contains(&name))
        .unwrap_or_else(|| panic!("missing auxiliary array {name}"))
}
fn payload(block: &str, compressed: bool) -> Vec<u8> {
    let body = block
        .split_once("<binary>")
        .unwrap()
        .1
        .split_once("</binary>")
        .unwrap()
        .0;
    let bytes = STANDARD.decode(body).unwrap();
    if compressed {
        let mut decoded = Vec::new();
        ZlibDecoder::new(bytes.as_slice())
            .read_to_end(&mut decoded)
            .unwrap();
        decoded
    } else {
        bytes
    }
}

#[test]
fn independent_binary_layout_preserves_empty_placeholders_empty_strings_and_integer_signs() {
    let mut input = sample();
    let spectrum = &mut input.spectra[0];
    spectrum.float_data_arrays = vec![
        DataArray::new("float-empty", vec![]),
        DataArray::new("float-bits", vec![-0.0, f32::from_bits(1)]),
    ];
    spectrum.integer_data_arrays = vec![
        DataArray::new("integer-empty", vec![]),
        DataArray::new("integer-extrema", vec![i32::MIN, i32::MAX]),
    ];
    spectrum.string_data_arrays = vec![
        DataArray::new("string-empty", vec![]),
        DataArray::new("two-empty-strings", vec![String::new(), String::new()]),
        DataArray::new(
            "ascii-controls",
            vec!["\u{1}\t\n\r\u{7f}".into(), "<&>\"".into()],
        ),
    ];
    for compressed in [false, true] {
        let bytes = encoded(&input, compressed);
        let xml = std::str::from_utf8(&bytes).unwrap();
        for name in ["float-empty", "integer-empty", "string-empty"] {
            let block = array(xml, name);
            assert!(
                block.contains("arrayLength=\"0\""),
                "placeholder must declare zero length"
            );
            assert!(payload(block, compressed).is_empty());
        }
        let strings = array(xml, "two-empty-strings");
        assert!(strings.contains("MS:1001479"));
        assert!(strings.contains("MS:1000786"));
        assert_eq!(payload(strings, compressed), b"\0\0");
        let integer = array(xml, "integer-extrema");
        assert!(integer.contains("MS:1000522"));
        let expected: Vec<u8> = [i64::from(i32::MIN), i64::from(i32::MAX)]
            .into_iter()
            .flat_map(i64::to_le_bytes)
            .collect();
        assert_eq!(payload(integer, compressed), expected);
        let floats = array(xml, "float-bits");
        assert!(floats.contains("MS:1000521"));
        let expected: Vec<u8> = [-0.0_f32, f32::from_bits(1)]
            .into_iter()
            .flat_map(f32::to_le_bytes)
            .collect();
        assert_eq!(payload(floats, compressed), expected);
        assert_eq!(
            payload(array(xml, "ascii-controls"), compressed),
            b"\x01\t\n\r\x7f\0<&>\"\0"
        );
        let restored = mzml::read(Cursor::new(bytes)).unwrap();
        assert_eq!(
            restored.spectra[0].float_data_arrays[1].data[0].to_bits(),
            (-0.0_f32).to_bits()
        );
        assert_eq!(restored, input);
    }
}

#[test]
fn names_and_within_type_order_survive_xml_escaping_on_both_record_kinds() {
    let mut input = sample();
    input
        .metadata
        .insert("run <&>".into(), "µ\r\n\t<&>\"".into());
    for spectrum in &mut input.spectra {
        spectrum.float_data_arrays = vec![
            DataArray::new("z float", vec![0.0, 1.0]),
            DataArray::new("a float µ\r\n\t<&>\"", vec![]),
        ];
        spectrum.integer_data_arrays = vec![
            DataArray::new("z integer", vec![]),
            DataArray::new("a integer", vec![-1, 2]),
        ];
        spectrum.string_data_arrays = vec![
            DataArray::new("z string", vec!["".into(), "A".into()]),
            DataArray::new("a string", vec![]),
        ];
        spectrum
            .metadata
            .insert("spectrum note".into(), "value & < >".into());
    }
    let spectrum = &input.spectra[0];
    let chromatogram = &mut input.chromatograms[0];
    chromatogram.float_data_arrays = spectrum.float_data_arrays.clone();
    chromatogram.integer_data_arrays = spectrum.integer_data_arrays.clone();
    chromatogram.string_data_arrays = spectrum.string_data_arrays.clone();
    chromatogram
        .metadata
        .insert("chromatogram note".into(), "value & < >".into());
    let before = input.clone();
    for compressed in [false, true] {
        let restored = mzml::read(Cursor::new(encoded(&input, compressed))).unwrap();
        assert_eq!(restored, input);
        assert_eq!(input, before);
    }
}

#[test]
fn auxiliary_writer_output_is_valid_against_pinned_mzml_schema_when_available() {
    use std::process::{Command, Stdio};
    if Command::new("xmllint").arg("--version").output().is_err() {
        eprintln!("xmllint unavailable; auxiliary writer XSD validation not executed on this host");
        return;
    }
    let schema = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/mzml_1_10.xsd");
    let mut populated = sample();
    populated.spectra[0]
        .string_data_arrays
        .push(DataArray::new("empty string placeholder", vec![]));
    let mut empty = MSExperiment::new();
    let mut spectrum = MSSpectrum::default();
    spectrum
        .float_data_arrays
        .push(DataArray::new("empty float", vec![]));
    spectrum
        .integer_data_arrays
        .push(DataArray::new("empty integer", vec![]));
    spectrum
        .string_data_arrays
        .push(DataArray::new("empty string", vec![]));
    empty.spectra.push(spectrum);
    let mut chromatogram = MSChromatogram::default();
    chromatogram
        .string_data_arrays
        .push(DataArray::new("empty chromatogram string", vec![]));
    empty.chromatograms.push(chromatogram);
    for input in [populated, empty] {
        for compressed in [false, true] {
            let mut process = Command::new("xmllint")
                .args(["--nonet", "--noout", "--schema"])
                .arg(&schema)
                .arg("-")
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .unwrap();
            process
                .stdin
                .take()
                .unwrap()
                .write_all(&encoded(&input, compressed))
                .unwrap();
            let result = process.wait_with_output().unwrap();
            assert!(
                result.status.success(),
                "{}",
                String::from_utf8_lossy(&result.stderr)
            );
        }
    }
}
