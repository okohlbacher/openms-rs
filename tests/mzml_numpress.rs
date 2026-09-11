// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "mzml")]
use base64::{Engine, engine::general_purpose::STANDARD};
use openms::format::mzml::{self, LoadOptions, NumpressWriteOptions, ReadOptions};
use openms::format::numpress_coder::{NumpressCompression as Mode, NumpressConfig};
use openms::kernel::{
    ChromatogramPeak, DataArray, MSChromatogram, MSExperiment, MSSpectrum, NumericRange, Peak1D,
};
use std::io::Cursor;

const LITERAL: [f64; 4] = [100., 200., 300.00005, 400.00010];
const LINEAR: &str = "QWR64UAAAADo//8/0P//f1kSgA==";
const PIC: &str = "ZGaMXCFQkQ==";
const SLOF: &str = "QMVagAAAAAAZxX3ivPP8/w==";
fn cv(id: &str) -> String {
    format!("<cvParam accession=\"{id}\" name=\"fixture\"/>")
}
fn array(role: &str, precision: &str, compression: &str, text: &str) -> String {
    format!(
        "<binaryDataArray encodedLength=\"{}\">{}{}{compression}<binary>{text}</binary></binaryDataArray>",
        text.len(),
        cv(role),
        if precision.is_empty() {
            String::new()
        } else {
            cv(precision)
        }
    )
}
fn ordinary(role: &str, values: &[f64]) -> String {
    array(
        role,
        "MS:1000523",
        &cv("MS:1000576"),
        &STANDARD.encode(
            values
                .iter()
                .flat_map(|v| v.to_le_bytes())
                .collect::<Vec<_>>(),
        ),
    )
}
fn spectrum(count: usize, arrays: &[String]) -> String {
    format!(
        "<spectrum id=\"scan=1\" defaultArrayLength=\"{count}\">{}<userParam name=\"transport\" value=\"kept\"/><binaryDataArrayList count=\"{}\">{}</binaryDataArrayList></spectrum>",
        cv("MS:1000127"),
        arrays.len(),
        arrays.join("")
    )
}
fn doc(count: usize, arrays: &[String]) -> String {
    format!(
        "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\"><run id=\"r\"><spectrumList count=\"1\">{}</spectrumList></run></mzML>",
        spectrum(count, arrays)
    )
}
fn read(text: &str) -> openms::Result<MSExperiment> {
    mzml::read(Cursor::new(text))
}
fn config(compression: Mode) -> NumpressConfig {
    NumpressConfig {
        compression,
        ..Default::default()
    }
}
fn all_options() -> NumpressWriteOptions {
    NumpressWriteOptions {
        mass_time: config(Mode::Linear),
        intensity: config(Mode::Pic),
        float_data_array: config(Mode::Slof),
        ..Default::default()
    }
}
fn populated() -> MSExperiment {
    let mut s = MSSpectrum::from_peaks(LITERAL.iter().map(|&v| Peak1D::new(v, v as f32)).collect());
    s.native_id = "scan=1".into();
    s.metadata.insert("tag".into(), "preserved".into());
    s.float_data_arrays.push(DataArray::new(
        "signal to noise array",
        LITERAL.map(|v| v as f32).to_vec(),
    ));
    s.integer_data_arrays
        .push(DataArray::new("charge array", vec![1, 2, 3, 4]));
    s.string_data_arrays.push(DataArray::new(
        "labels",
        vec!["a".into(), "b".into(), "c".into(), "d".into()],
    ));
    let mut c = MSChromatogram::from_peaks(
        LITERAL
            .iter()
            .map(|&v| ChromatogramPeak::new(v, v as f32))
            .collect(),
    );
    c.native_id = "transition".into();
    MSExperiment {
        spectra: vec![s],
        chromatograms: vec![c],
        ..Default::default()
    }
}

#[test]
fn published_bytes_all_six_accessions_and_separate_zlib() {
    for row in include_str!("data/numpress_coder_transport.tsv")
        .lines()
        .skip(1)
    {
        let r: Vec<_> = row.split('\t').collect();
        let (plain, combined, tolerance) = match r[0] {
            "linear" => ("MS:1002312", "MS:1002746", 1e-5),
            "pic" => ("MS:1002313", "MS:1002747", 0.001),
            _ => ("MS:1002314", "MS:1002748", 0.1),
        };
        for (terms, text) in [
            (cv(plain), r[1]),
            (cv(combined), r[2]),
            (format!("{}{}", cv(plain), cv("MS:1000574")), r[2]),
            (format!("{}{}", cv("MS:1000574"), cv(plain)), r[2]),
        ] {
            for precision in ["MS:1000523", "MS:1000521", ""] {
                let xml = doc(
                    4,
                    &[
                        array("MS:1000514", precision, &terms, text),
                        ordinary("MS:1000515", &LITERAL),
                    ],
                );
                let e = read(&xml).unwrap();
                for (p, expected) in e.spectra[0].peaks.iter().zip(LITERAL) {
                    assert!((p.mz - expected).abs() < tolerance);
                }
                assert_eq!(e.spectra[0].metadata["transport"], "kept");
            }
        }
    }
}
#[test]
fn source_real_file_payloads_independent_decode_oracle() {
    let e = read(include_str!("data/mzml_numpress_source_projection.mzML")).unwrap();
    assert_eq!(e.chromatograms.len(), 18);
    let mut count = 0;
    for row in include_str!("data/mzml_numpress_source_values.tsv")
        .lines()
        .skip(1)
    {
        let r: Vec<_> = row.split('\t').collect();
        let i: usize = r[0].parse().unwrap();
        let j: usize = r[1].parse().unwrap();
        let p = e.chromatograms[i].peaks[j];
        assert_eq!(p.rt.to_bits(), u64::from_str_radix(r[3], 16).unwrap());
        let expected = f32::from_bits(u32::from_str_radix(r[5], 16).unwrap());
        assert!((p.intensity - expected).abs() <= 2. * f32::EPSILON * expected.abs().max(1.));
        count += 1;
    }
    assert_eq!(count, 342);
}
#[test]
fn writer_three_configs_and_ordinary_metadata_roundtrip() {
    for compressed in [false, true] {
        let mut options = all_options();
        options.binary.zlib_compression = compressed;
        let mut bytes = Vec::new();
        let report = mzml::write_with_numpress(&mut bytes, &populated(), &options).unwrap();
        assert_eq!(
            (
                report.encoded_arrays,
                report.ordinary_arrays,
                report.fallback_arrays
            ),
            (5, 2, 0)
        );
        let text = String::from_utf8(bytes).unwrap();
        if !compressed {
            for literal in [LINEAR, PIC, SLOF] {
                assert!(text.contains(literal), "{literal}");
            }
        }
        for term in if compressed {
            ["MS:1002746", "MS:1002747", "MS:1002748"]
        } else {
            ["MS:1002312", "MS:1002313", "MS:1002314"]
        } {
            assert!(text.contains(term));
        }
        let e = read(&text).unwrap();
        assert_eq!(e.spectra[0].metadata["tag"], "preserved");
        assert_eq!(
            e.spectra[0].integer_data_arrays,
            populated().spectra[0].integer_data_arrays
        );
        assert_eq!(
            e.spectra[0].string_data_arrays,
            populated().spectra[0].string_data_arrays
        );
        for (p, v) in e.spectra[0].peaks.iter().zip(LITERAL) {
            assert!((p.mz - v).abs() < 1e-5);
            assert!((f64::from(p.intensity) - v).abs() < 0.001);
        }
        for (p, v) in e.spectra[0].float_data_arrays[0].data.iter().zip(LITERAL) {
            assert!((f64::from(*p) - v).abs() < 0.1);
        }
        assert!((e.chromatograms[0].peaks[3].rt - LITERAL[3]).abs() < 1e-5);
    }
}
#[test]
fn default_entrypoint_is_byte_identical_to_ordinary_writer() {
    for compressed in [false, true] {
        let mut options = NumpressWriteOptions::default();
        options.binary.zlib_compression = compressed;
        let mut ordinary = Vec::new();
        mzml::write_with_options(&mut ordinary, &populated(), &options.binary).unwrap();
        let mut new = Vec::new();
        let report = mzml::write_with_numpress(&mut new, &populated(), &options).unwrap();
        assert_eq!(new, ordinary);
        assert_eq!(report.ordinary_arrays, 7);
    }
}
#[test]
fn codec_accuracy_empty_and_canonical_role_fallback_are_ordinary() {
    let mut e = populated();
    e.spectra[0]
        .float_data_arrays
        .push(DataArray::new("mean charge array", vec![1., 2., 3., 4.]));
    e.spectra[0]
        .float_data_arrays
        .push(DataArray::new("empty", vec![]));
    // Finite negative PIC input fails the raw encoder; a coarse linear factor
    // encodes successfully but fails its source error check.
    e.spectra[0].peaks[0].intensity = -2.;
    let mut options = all_options();
    options.mass_time.estimate_fixed_point = false;
    options.mass_time.fixed_point = 0.001;
    let mut bytes = Vec::new();
    let report = mzml::write_with_numpress(&mut bytes, &e, &options).unwrap();
    assert_eq!(report.fallback_arrays, 5); // two coordinates, one intensity, canonical+empty
    let decoded = read(std::str::from_utf8(&bytes).unwrap()).unwrap();
    assert_eq!(decoded.spectra[0].peaks, e.spectra[0].peaks);
    assert_eq!(
        decoded.spectra[0].float_data_arrays[1],
        e.spectra[0].float_data_arrays[1]
    );
    assert!(decoded.spectra[0].float_data_arrays[2].data.is_empty());
}
#[test]
fn precision_repairs_are_numpress_only_and_preserve_canonical_identity() {
    for precision in ["MS:1000519", "MS:1000522"] {
        let e = read(&doc(
            4,
            &[
                array("MS:1000514", precision, &cv("MS:1002313"), PIC),
                ordinary("MS:1000515", &LITERAL),
            ],
        ))
        .unwrap();
        assert_eq!(e.spectra[0].peaks[2].mz, 300.);
        for mode in ["MS:1002312", "MS:1002314"] {
            assert!(
                read(&doc(
                    4,
                    &[
                        array("MS:1000514", precision, &cv(mode), LINEAR),
                        ordinary("MS:1000515", &LITERAL)
                    ]
                ))
                .is_err()
            );
        }
    }
    for (precision, role) in [
        ("MS:1001479", "MS:1000514"),
        ("MS:1000519", "MS:1000516"),
        ("MS:1000521", "MS:1002478"),
    ] {
        assert!(
            read(&doc(
                4,
                &[
                    ordinary("MS:1000514", &LITERAL),
                    ordinary("MS:1000515", &LITERAL),
                    array(role, precision, &cv("MS:1002313"), PIC)
                ]
            ))
            .is_err()
        );
    }
    assert!(
        read(&doc(
            4,
            &[
                array("MS:1000514", "", &cv("MS:1000576"), LINEAR),
                ordinary("MS:1000515", &LITERAL)
            ]
        ))
        .is_err()
    );
}
#[test]
fn conflicting_encoding_and_malformed_payloads_do_not_silently_decode() {
    for terms in [
        format!("{}{}", cv("MS:1002312"), cv("MS:1002313")),
        format!("{}{}", cv("MS:1002312"), cv("MS:1002312")),
        format!("{}{}", cv("MS:1000576"), cv("MS:1002312")),
        format!("{}{}", cv("MS:1002312"), cv("MS:1000576")),
        format!("{}{}", cv("MS:1002746"), cv("MS:1000574")),
        format!("{}{}", cv("MS:1000574"), cv("MS:1002746")),
    ] {
        assert!(
            read(&doc(
                4,
                &[
                    array("MS:1000514", "MS:1000523", &terms, LINEAR),
                    ordinary("MS:1000515", &LITERAL)
                ]
            ))
            .is_err()
        );
    }
    for text in [
        "?",
        "???",
        "!!!!",
        "A===",
        "QWR64UAAAADo//8/0P//f1kSgA=",
        "AAAAAA==",
    ] {
        assert!(
            read(&doc(
                4,
                &[
                    array("MS:1000514", "MS:1000523", &cv("MS:1002312"), text),
                    ordinary("MS:1000515", &LITERAL)
                ]
            ))
            .is_err(),
            "{text}"
        );
    }
    // A syntactically valid zlib member with a damaged checksum cannot be used.
    let mut compressed = STANDARD
        .decode("eJxzTKl66MDAwPDi/3/7C///10cKNQAAWQUJng==")
        .unwrap();
    let last = compressed.len() - 1;
    compressed[last] ^= 1;
    assert!(
        read(&doc(
            4,
            &[
                array(
                    "MS:1000514",
                    "MS:1000523",
                    &cv("MS:1002746"),
                    &STANDARD.encode(compressed)
                ),
                ordinary("MS:1000515", &LITERAL)
            ]
        ))
        .is_err()
    );
}
#[test]
fn exact_counts_and_existing_cumulative_binary_limits_apply() {
    for count in [0, 3, 5] {
        assert!(
            read(&doc(
                count,
                &[
                    array("MS:1000514", "MS:1000523", &cv("MS:1002312"), LINEAR),
                    ordinary("MS:1000515", &LITERAL)
                ]
            ))
            .is_err()
        );
    }
    let xml = doc(
        4,
        &[
            array("MS:1000514", "MS:1000523", &cv("MS:1002312"), LINEAR),
            array("MS:1000515", "MS:1000523", &cv("MS:1002313"), PIC),
        ],
    );
    for limits in [
        ReadOptions {
            max_array_bytes: 31,
            ..Default::default()
        },
        ReadOptions {
            max_total_array_bytes: 63,
            ..Default::default()
        },
        ReadOptions {
            max_total_array_elements: 7,
            ..Default::default()
        },
    ] {
        assert!(mzml::read_with_options(Cursor::new(&xml), &limits).is_err());
    }
    let limits = ReadOptions {
        max_array_bytes: 32,
        max_total_array_bytes: 64,
        max_total_array_elements: 8,
        ..Default::default()
    };
    assert!(mzml::read_with_options(Cursor::new(xml), &limits).is_ok());
    let empty = doc(
        0,
        &[
            array("MS:1000514", "MS:1000523", &cv("MS:1002312"), ""),
            ordinary("MS:1000515", &[]),
        ],
    );
    assert!(read(&empty).is_ok());
    assert!(read(&empty.replace("MS:1002312", "MS:1002746")).is_err());
}
#[test]
fn referenced_compression_matches_inline_and_minutes_are_scaled() {
    let inline = doc(
        4,
        &[
            array("MS:1000514", "MS:1000523", &cv("MS:1002312"), LINEAR),
            ordinary("MS:1000515", &LITERAL),
        ],
    );
    let reference=inline.replace("<run id=\"r\">",&format!("<referenceableParamGroupList count=\"1\"><referenceableParamGroup id=\"np\">{}</referenceableParamGroup></referenceableParamGroupList><run id=\"r\">",cv("MS:1002312"))).replacen(&format!("{}<binary>",cv("MS:1002312")),"<referenceableParamGroupRef ref=\"np\"/><binary>",1);
    assert_eq!(read(&reference).unwrap(), read(&inline).unwrap());
    let minutes = inline
        .replace("spectrumList", "chromatogramList")
        .replace("spectrum ", "chromatogram ")
        .replace("</spectrum>", "</chromatogram>")
        .replace(&cv("MS:1000127"), "")
        .replace(
            &cv("MS:1000514"),
            "<cvParam accession=\"MS:1000595\" unitAccession=\"UO:0000031\"/>",
        );
    let e = read(&minutes).unwrap();
    assert!((e.chromatograms[0].peaks[3].rt - 60. * LITERAL[3]).abs() < 0.001);
}
#[test]
fn selection_filters_decoded_f64_before_primary_and_auxiliary_narrowing() {
    // Independent SLOF words with factor1 decode exp(1)-1 and exp(300)-1.
    let mut raw = 1_f64.to_be_bytes().to_vec();
    raw.extend_from_slice(&1_u16.to_le_bytes());
    raw.extend_from_slice(&300_u16.to_le_bytes());
    let text = STANDARD.encode(raw);
    let xml = doc(
        2,
        &[
            ordinary("MS:1000514", &[100., 200.]),
            array("MS:1000515", "MS:1000523", &cv("MS:1002314"), &text),
            array("MS:1000517", "MS:1000523", &cv("MS:1002314"), &text),
        ],
    );
    assert!(read(&xml).is_err());
    let mut load = LoadOptions::default();
    load.scientific
        .set_intensity_range(NumericRange { min: 0., max: 10. });
    let e =
        mzml::read_with_load_options(Cursor::new(&xml), &load, &ReadOptions::default()).unwrap();
    assert_eq!(e.spectra[0].len(), 1);
    assert_eq!(e.spectra[0].peaks[0].intensity, (1_f64.exp() - 1.) as f32);
    assert_eq!(
        e.spectra[0].float_data_arrays[0].data,
        vec![(1_f64.exp() - 1.) as f32]
    );
    load.skip_spectra = true;
    assert!(
        mzml::read_with_load_options(Cursor::new(&xml), &load, &ReadOptions::default()).is_err()
    );
}
#[test]
fn writer_resource_failures_are_cumulative_and_precede_output() {
    let mut e = populated();
    e.chromatograms.clear();
    e.spectra[0].float_data_arrays.clear();
    e.spectra[0].integer_data_arrays.clear();
    e.spectra[0].string_data_arrays.clear();
    let mut options = all_options();
    options.limits.raw.max_work = 3000;
    let mut bytes = Vec::new();
    assert!(mzml::write_with_numpress(&mut bytes, &e, &options).is_ok());
    let mut second = e.spectra[0].clone();
    second.native_id = "scan=2".into();
    e.spectra.push(second);
    bytes = b"unchanged".to_vec();
    assert!(mzml::write_with_numpress(&mut bytes, &e, &options).is_err());
    assert_eq!(bytes, b"unchanged");
    options.limits.raw.max_work = usize::MAX;
    options.limits.max_total_bytes = 1;
    assert!(mzml::write_with_numpress(&mut bytes, &e, &options).is_err());
    assert_eq!(bytes, b"unchanged");
    options.limits = Default::default();
    e.spectra[0]
        .float_data_arrays
        .push(DataArray::new("late metadata", vec![1.; 4]));
    e.spectra[0].float_data_arrays[0]
        .metadata
        .insert("unrepresented".into(), "value".into());
    assert!(mzml::write_with_numpress(&mut bytes, &e, &options).is_err());
    assert_eq!(bytes, b"unchanged");
}

#[test]
fn compressed_and_plain_numpress_writer_pass_pinned_xsd() {
    if std::process::Command::new("xmllint")
        .arg("--version")
        .output()
        .is_err()
    {
        eprintln!("xmllint unavailable; independent XSD check skipped");
        return;
    }
    for compressed in [false, true] {
        let mut options = all_options();
        options.binary.zlib_compression = compressed;
        let mut bytes = Vec::new();
        mzml::write_with_numpress(&mut bytes, &populated(), &options).unwrap();
        let path = std::env::temp_dir().join(format!(
            "openms-numpress-schema-{}-{compressed}.mzML",
            std::process::id()
        ));
        std::fs::write(&path, bytes).unwrap();
        let output = std::process::Command::new("xmllint")
            .args(["--nonet", "--noout", "--schema"])
            .arg(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/tests/data/mzml_1_10.xsd"
            ))
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
fn numpress_selection_sorts_all_aligned_array_kinds_together() {
    // Independent PIC nibble oracle: 300 => 5,c,2,1; 200 => 6,8,c;
    // 100 => 6,4,6. Concatenated high-nibble-first bytes below.
    let text = STANDARD.encode([0x5c, 0x21, 0x68, 0xc6, 0x46]);
    let ints = array(
        "MS:1000516",
        "MS:1000519",
        &cv("MS:1000576"),
        &STANDARD.encode(
            [3_i32, 2, 1]
                .into_iter()
                .flat_map(i32::to_le_bytes)
                .collect::<Vec<_>>(),
        ),
    );
    let labels = array(
        "MS:1000786",
        "MS:1001479",
        &cv("MS:1000576"),
        &STANDARD.encode(b"c\0b\0a\0"),
    )
    .replace(
        &cv("MS:1000786"),
        "<cvParam accession=\"MS:1000786\" value=\"labels\"/>",
    );
    let xml = doc(
        3,
        &[
            array("MS:1000514", "MS:1000523", &cv("MS:1002313"), &text),
            ordinary("MS:1000515", &[30., 20., 10.]),
            array("MS:1000517", "MS:1000523", &cv("MS:1002313"), &text),
            ints,
            labels,
        ],
    );
    let mut load = LoadOptions::default();
    load.scientific.set_mz_range(NumericRange {
        min: 100.,
        max: 300.,
    });
    let e = mzml::read_with_load_options(Cursor::new(xml), &load, &ReadOptions::default()).unwrap();
    let s = &e.spectra[0];
    assert_eq!(
        s.peaks,
        vec![Peak1D::new(100., 10.), Peak1D::new(200., 20.)]
    );
    assert_eq!(s.float_data_arrays[0].data, vec![100., 200.]);
    assert_eq!(s.integer_data_arrays[0].data, vec![1, 2]);
    assert_eq!(s.string_data_arrays[0].data, vec!["a", "b"]);
}

#[test]
fn binary_count_and_work_limits_precede_native_value_validation() {
    let mut s = MSSpectrum::from_peaks(vec![Peak1D::new(100., 1.); 2048]);
    s.peaks[2047].mz = f64::NAN;
    let e = MSExperiment {
        spectra: vec![s],
        ..Default::default()
    };
    let mut options = all_options();
    options.limits.raw.max_values = 1;
    let mut bytes = b"untouched".to_vec();
    let error = mzml::write_with_numpress(&mut bytes, &e, &options).unwrap_err();
    assert!(error.to_string().contains("value count limit"));
    assert_eq!(bytes, b"untouched");
    options.limits.raw.max_values = 4096;
    options.limits.raw.max_work = 8;
    let error = mzml::write_with_numpress(&mut bytes, &e, &options).unwrap_err();
    assert!(error.to_string().contains("work limit"));
    assert_eq!(bytes, b"untouched");
    let mut e = MSExperiment {
        spectra: vec![MSSpectrum::from_peaks(vec![Peak1D::new(100., 1.)])],
        ..Default::default()
    };
    e.spectra[0]
        .float_data_arrays
        .push(DataArray::new("over limit", vec![f32::NAN; 2048]));
    options.limits = Default::default();
    options.limits.raw.max_values = 1;
    assert!(
        mzml::write_with_numpress(&mut bytes, &e, &options)
            .unwrap_err()
            .to_string()
            .contains("value count limit")
    );
}
