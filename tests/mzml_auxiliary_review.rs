// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Source conventions inspected at OpenMS4-core 7c029e8cdba6abab503708ecdd56f6ab55e38ce4:
// MzMLHandler.cpp SHA256 8546ea20bbfd478e5464001e0a3e0cb931142eed250c244f4f19113159d9198a
// MzMLHandlerHelper.cpp 339b9f224fe3b142d2af42b65869b838549810d379cd406e85609ad499ae62f9
// Base64.cpp e72ee0f62de04a2a30ece314bfa8db2587950700d0868a34cffe2dc5fd87069e
// psi-ms.obo 1623792d5fd37ab305bc7228ab6b51839cfde06130f9b8727bd8ce99993a27d3

#![cfg(feature = "mzml")]

use base64::{Engine, engine::general_purpose::STANDARD};
use flate2::{Compression, write::ZlibEncoder};
use openms::format::mzml::{self, ReadOptions, WriteOptions};
use openms::kernel::{DataArray, MSExperiment, MSSpectrum, Peak1D};
use std::io::{Cursor, Write};

fn binary_array(precision: &str, kind: &str, bytes: &[u8], count: usize, zlib: bool) -> String {
    let encoded = if zlib {
        let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(bytes).unwrap();
        STANDARD.encode(encoder.finish().unwrap())
    } else {
        STANDARD.encode(bytes)
    };
    let compression = if zlib { "MS:1000574" } else { "MS:1000576" };
    format!(
        "<binaryDataArray arrayLength=\"{count}\" encodedLength=\"{}\"><cvParam accession=\"{precision}\"/><cvParam accession=\"{compression}\"/>{kind}<binary>{encoded}</binary></binaryDataArray>",
        encoded.len()
    )
}
fn auxiliary(precision: &str, bytes: &[u8], count: usize, zlib: bool) -> String {
    binary_array(
        precision,
        "<cvParam accession=\"MS:1000786\" value=\"review array\"/>",
        bytes,
        count,
        zlib,
    )
}
fn document(count: usize, auxiliary_arrays: &[String]) -> String {
    // Independent bytes and XML construction, not the native mzML writer.
    let mz: Vec<_> = (0..count)
        .flat_map(|i| (100.0 + i as f64).to_le_bytes())
        .collect();
    let intensity: Vec<_> = (0..count)
        .flat_map(|i| (i as f32 + 1.0).to_le_bytes())
        .collect();
    let primary_mz = binary_array(
        "MS:1000523",
        "<cvParam accession=\"MS:1000514\"/>",
        &mz,
        count,
        false,
    );
    let primary_intensity = binary_array(
        "MS:1000521",
        "<cvParam accession=\"MS:1000515\"/>",
        &intensity,
        count,
        false,
    );
    format!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?><mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\"><run id=\"run\"><spectrumList count=\"1\"><spectrum id=\"scan=1\" index=\"0\" defaultArrayLength=\"{count}\"><binaryDataArrayList count=\"{}\">{primary_mz}{primary_intensity}{}</binaryDataArrayList></spectrum></spectrumList></run></mzML>",
        auxiliary_arrays.len() + 2,
        auxiliary_arrays.join("")
    )
}
fn read(xml: &str) -> openms::Result<MSExperiment> {
    mzml::read(Cursor::new(xml.as_bytes()))
}
fn read_options(xml: &str, options: &ReadOptions) -> openms::Result<MSExperiment> {
    mzml::read_with_options(Cursor::new(xml.as_bytes()), options)
}

#[test]
fn independent_signed_integer_widths_preserve_native_extremes() {
    let values = [i32::MIN, 0, i32::MAX];
    let int32: Vec<_> = values.into_iter().flat_map(i32::to_le_bytes).collect();
    let int64: Vec<_> = values
        .into_iter()
        .flat_map(|v| i64::from(v).to_le_bytes())
        .collect();
    for (precision, bytes) in [("MS:1000519", int32), ("MS:1000522", int64)] {
        for compressed in [false, true] {
            let result =
                read(&document(3, &[auxiliary(precision, &bytes, 3, compressed)])).unwrap();
            assert_eq!(
                result.spectra[0].integer_data_arrays,
                [DataArray::new("review array", values.to_vec())]
            );
        }
    }
    for invalid in [
        i64::from(i32::MIN) - 1,
        i64::from(i32::MAX) + 1,
        i64::MIN,
        i64::MAX,
    ] {
        assert!(
            read(&document(
                1,
                &[auxiliary("MS:1000522", &invalid.to_le_bytes(), 1, false)]
            ))
            .is_err()
        );
    }
}

#[test]
fn independent_float_widths_preserve_auxiliary_values_and_negative_zero() {
    for compressed in [false, true] {
        for precision in ["MS:1000521", "MS:1000523"] {
            let values = [-0.0_f64, -1.25, 16_777_216.0];
            let bytes: Vec<_> = if precision == "MS:1000521" {
                values
                    .into_iter()
                    .flat_map(|v| (v as f32).to_le_bytes())
                    .collect()
            } else {
                values.into_iter().flat_map(f64::to_le_bytes).collect()
            };
            let result =
                read(&document(3, &[auxiliary(precision, &bytes, 3, compressed)])).unwrap();
            let values = &result.spectra[0].float_data_arrays[0].data;
            assert_eq!(values[0].to_bits(), (-0.0_f32).to_bits());
            assert_eq!(values[1..], [-1.25, 16_777_216.]);
        }
    }
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1e100] {
        assert!(
            read(&document(
                1,
                &[auxiliary("MS:1000523", &invalid.to_le_bytes(), 1, false)]
            ))
            .is_err()
        );
    }
}

#[test]
fn source_string_encoding_preserves_empty_entries_and_ascii_controls() {
    // MS:1001479 explicitly allows zero non-NUL ASCII characters: an empty
    // element is one NUL. Source encodeStrings emits all of these separators.
    for compressed in [false, true] {
        let result = read(&document(
            4,
            &[auxiliary("MS:1001479", b"\0A\t\r\n\0\0z\0", 4, compressed)],
        ))
        .unwrap();
        assert_eq!(
            result.spectra[0].string_data_arrays[0].data,
            ["", "A\t\r\n", "", "z"]
        );
        let empty = read(&document(0, &[auxiliary("MS:1001479", b"", 0, compressed)])).unwrap();
        assert!(empty.spectra[0].string_data_arrays[0].data.is_empty());
    }
}

#[test]
fn string_reader_rejects_missing_terminators_non_ascii_and_wrong_element_counts() {
    for (bytes, count) in [
        (b"ABC".as_slice(), 1),
        (b"A\0tail", 1),
        (b"A\0\0", 1),
        (b"A\0", 2),
        (b"\xc2\xb5\0", 1),
        (b"\xff\0", 1),
        (b"\0", 0),
    ] {
        for compressed in [false, true] {
            assert!(
                read(&document(
                    count,
                    &[auxiliary("MS:1001479", bytes, count, compressed)]
                ))
                .is_err(),
                "bytes={bytes:?}, count={count}, zlib={compressed}"
            );
        }
    }
}

#[test]
fn string_decompression_respects_actual_bytes_with_variable_element_length() {
    let options = ReadOptions {
        max_array_bytes: 64,
        ..Default::default()
    };
    let mut at_limit = vec![b'A'; 63];
    at_limit.push(0);
    let xml = document(1, &[auxiliary("MS:1001479", &at_limit, 1, true)]);
    let result = read_options(&xml, &options).unwrap();
    assert_eq!(result.spectra[0].string_data_arrays[0].data[0].len(), 63);
    for length in [64, 4096] {
        let mut excessive = vec![b'A'; length];
        excessive.push(0);
        let xml = document(1, &[auxiliary("MS:1001479", &excessive, 1, true)]);
        assert!(read_options(&xml, &options).is_err());
    }
}

fn sample() -> MSExperiment {
    let mut spectrum = MSSpectrum::from_peaks(vec![
        Peak1D::new(100., 1.),
        Peak1D::new(200., 2.),
        Peak1D::new(300., 3.),
    ]);
    spectrum.native_id = "scan=1".into();
    spectrum
        .float_data_arrays
        .push(DataArray::new("µ & \"float\"", vec![-0.0, 2.5, -3.]));
    spectrum
        .integer_data_arrays
        .push(DataArray::new("integer", vec![i32::MIN, 0, i32::MAX]));
    spectrum.string_data_arrays.push(DataArray::new(
        "label",
        vec!["".into(), "a\tb\r\nc".into(), "".into()],
    ));
    MSExperiment {
        spectra: vec![spectrum],
        ..Default::default()
    }
}

#[test]
fn writer_roundtrip_retains_empty_strings_and_escapes_unicode_array_names() {
    let input = sample();
    for compressed in [false, true] {
        let mut bytes = Vec::new();
        mzml::write_with_options(
            &mut bytes,
            &input,
            &WriteOptions {
                zlib_compression: compressed,
            },
        )
        .unwrap();
        let text = std::str::from_utf8(&bytes).unwrap();
        assert!(text.contains("accession=\"MS:1000786\""));
        assert!(text.contains("accession=\"MS:1001479\""));
        assert_eq!(mzml::read(Cursor::new(bytes)).unwrap(), input);
    }
}

#[test]
fn writer_rejects_non_ascii_or_embedded_nul_strings_before_any_output() {
    for invalid in ["µ", "a\0b"] {
        let mut input = sample();
        input.spectra[0].string_data_arrays[0].data[1] = invalid.into();
        let mut output = Vec::new();
        assert!(mzml::write(&mut output, &input).is_err());
        assert!(output.is_empty());
    }
}

#[test]
fn cumulative_limits_include_primary_arrays_and_empty_string_elements() {
    // Two f64 coordinates + two f32 intensities + two empty ASCII strings:
    // 26 decoded bytes, six elements, three arrays, regardless of compression.
    for compressed in [false, true] {
        let xml = document(2, &[auxiliary("MS:1001479", b"\0\0", 2, compressed)]);
        let limits = ReadOptions {
            max_total_array_bytes: 26,
            max_total_array_elements: 6,
            max_total_arrays: 3,
            ..Default::default()
        };
        assert!(read_options(&xml, &limits).is_ok());
        for below in [
            ReadOptions {
                max_total_array_bytes: 25,
                ..limits
            },
            ReadOptions {
                max_total_array_elements: 5,
                ..limits
            },
            ReadOptions {
                max_total_arrays: 2,
                ..limits
            },
        ] {
            assert!(read_options(&xml, &below).is_err());
        }
    }
}

#[test]
fn cumulative_array_count_also_bounds_zero_length_placeholders() {
    let extras = [
        auxiliary("MS:1001479", b"", 0, false),
        auxiliary("MS:1001479", b"", 0, false).replace("review array", "second array"),
    ];
    // A populated native record may carry explicit empty annotation placeholders.
    let xml = document(2, &extras);
    let exactly = ReadOptions {
        max_total_arrays: 4,
        max_total_array_bytes: 24,
        max_total_array_elements: 4,
        ..Default::default()
    };
    let result = read_options(&xml, &exactly).unwrap();
    assert_eq!(result.spectra[0].string_data_arrays.len(), 2);
    assert!(
        result.spectra[0]
            .string_data_arrays
            .iter()
            .all(|a| a.data.is_empty())
    );
    assert!(
        read_options(
            &xml,
            &ReadOptions {
                max_total_arrays: 3,
                ..exactly
            }
        )
        .is_err()
    );
}

#[test]
fn auxiliary_units_are_not_silently_dropped_from_any_binary_cv_term() {
    let original = auxiliary("MS:1000521", &1.0_f32.to_le_bytes(), 1, false);
    for accession in ["MS:1000786", "MS:1000521", "MS:1000576"] {
        let changed = original.replace(
            &format!("accession=\"{accession}\""),
            &format!("accession=\"{accession}\" unitAccession=\"UO:0000010\" unitCvRef=\"UO\" unitName=\"second\""),
        );
        assert!(
            read(&document(1, &[changed])).is_err(),
            "accession={accession}"
        );
    }
}
