// Copyright (c) 2026 OpenMS Rust contributors
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "numpress")]
use base64::{Engine, engine::general_purpose::STANDARD};
use flate2::{Compression, write::ZlibEncoder};
use openms::format::numpress;
use openms::format::numpress_coder::*;
use std::io::Write;

fn config(mode: NumpressCompression) -> NumpressConfig {
    NumpressConfig {
        compression: mode,
        ..Default::default()
    }
}
fn literal() -> [f64; 4] {
    [100., 200., 300.00005, 400.00010]
}
fn failed(report: &NumpressEncodeReport<impl std::fmt::Debug>) -> &NumpressRejection {
    match &report.status {
        NumpressEncodeStatus::Rejected(reason) => reason,
        _ => panic!("expected rejection: {report:?}"),
    }
}
#[test]
fn literal_source_base64_and_raw_bytes_are_exact_for_all_modes() {
    let coder = MSNumpressCoder::default();
    for row in include_str!("data/numpress_source_bytes.tsv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = row.split('\t').collect();
        let config = config(f[0].parse().unwrap());
        let encoded = coder.encode(&literal(), false, &config).unwrap();
        assert_eq!(encoded.status, NumpressEncodeStatus::Encoded);
        assert_eq!(encoded.output, f[1]);
        let bytes = coder.encode_raw(&literal(), &config).unwrap();
        assert_eq!(STANDARD.encode(&bytes.output), f[1]);
        let decoded = coder.decode(f[1], false, &config).unwrap();
        assert_eq!(decoded, coder.decode_raw(&bytes.output, &config).unwrap());
        for (actual, expected) in decoded.iter().zip(literal()) {
            assert!(
                (actual - expected).abs()
                    <= if config.compression == NumpressCompression::Slof {
                        expected * 1e-4
                    } else {
                        0.001
                    }
            );
        }
    }
}
#[test]
fn existing_configuration_defaults_and_exact_names_are_reused() {
    let cfg = NumpressConfig::default();
    assert_eq!(cfg.fixed_point, 0.);
    assert_eq!(cfg.error_tolerance, 0.0001);
    assert_eq!(cfg.compression, NumpressCompression::None);
    assert!(cfg.estimate_fixed_point);
    assert_eq!(cfg.linear_fp_mass_acc, -1.);
    for mode in NumpressCompression::ALL {
        assert_eq!(mode.name().parse::<NumpressCompression>().unwrap(), mode);
    }
    let mut cfg = cfg;
    cfg.set_compression("linear").unwrap();
    let before = cfg;
    for text in ["LINEAR", " linear", "linear ", "safe", ""] {
        assert!(cfg.set_compression(text).is_err());
        assert_eq!(cfg, before);
    }
}
#[test]
fn empty_none_and_raw_versus_text_destination_semantics() {
    let coder = MSNumpressCoder::default();
    let mut bytes = vec![42];
    let mut text = "old".into();
    let none = NumpressConfig {
        fixed_point: f64::NAN,
        error_tolerance: f64::NAN,
        ..Default::default()
    };
    assert_eq!(
        coder
            .encode_raw_into(&[f64::NAN], &mut bytes, &none)
            .unwrap(),
        NumpressEncodeStatus::Disabled
    );
    assert_eq!(bytes, [42]);
    assert_eq!(
        coder
            .encode_into(&[f64::NAN], &mut text, true, &none)
            .unwrap(),
        NumpressEncodeStatus::Disabled
    );
    assert_eq!(text, "");
    assert_eq!(
        coder.encode_raw(&[], &none).unwrap().status,
        NumpressEncodeStatus::EmptyInput
    );
    let pic = config(NumpressCompression::Pic);
    let status = coder.encode_raw_into(&[-1.], &mut bytes, &pic).unwrap();
    assert!(matches!(status, NumpressEncodeStatus::Rejected(_)));
    assert_eq!(bytes, [42]);
    text = "old".into();
    let status = coder.encode_into(&[-1.], &mut text, false, &pic).unwrap();
    assert!(matches!(status, NumpressEncodeStatus::Rejected(_)));
    assert!(text.is_empty());
    assert!(
        coder
            .decode_raw(b"malformed bytes", &none)
            .unwrap()
            .is_empty()
    );
    assert!(coder.decode("YWJj", false, &none).unwrap().is_empty());
    assert!(coder.decode("!@#$", false, &none).is_err()); // transport still executes before NONE
}
#[test]
fn linear_accuracy_estimation_fallback_and_source_short_vector_quirk() {
    let coder = MSNumpressCoder::default();
    let mut cfg = config(NumpressCompression::Linear);
    cfg.linear_fp_mass_acc = 0.001;
    let report = coder.encode_raw(&literal(), &cfg).unwrap();
    assert!(report.is_encoded());
    assert_eq!(report.fixed_point, Some(500.));
    assert!(!report.used_maximal_fixed_point_fallback);
    cfg.linear_fp_mass_acc = 1e-12;
    let report = coder.encode_raw(&literal(), &cfg).unwrap();
    assert!(report.is_encoded());
    assert!(report.used_maximal_fixed_point_fallback);
    assert_eq!(
        report.fixed_point,
        Some(numpress::optimal_linear_fixed_point(&literal()).unwrap())
    );
    cfg.estimate_fixed_point = false;
    cfg.fixed_point = 1000.;
    let report = coder.encode_raw(&literal(), &cfg).unwrap();
    assert!(report.is_encoded());
    assert_eq!(report.fixed_point, Some(1000.));
    assert!(!report.used_maximal_fixed_point_fallback);
    cfg.estimate_fixed_point = true;
    cfg.linear_fp_mass_acc = 0.001;
    for input in [&[100.][..], &[100., 200.][..]] {
        let report = coder.encode_raw(input, &cfg).unwrap();
        assert_eq!(report.fixed_point, Some(0.));
        assert!(matches!(failed(&report), NumpressRejection::Codec(_)));
        cfg.error_tolerance = 0.;
        let unchecked = coder.encode_raw(input, &cfg).unwrap();
        assert!(unchecked.is_encoded());
        assert!(numpress::decode_linear(&unchecked.output).is_err());
        cfg.error_tolerance = 0.0001;
    }
}
#[test]
fn tolerance_reverse_order_ratio_zero_branch_and_exact_boundary() {
    let coder = MSNumpressCoder::default();
    let mut cfg = NumpressConfig {
        compression: NumpressCompression::Linear,
        estimate_fixed_point: false,
        fixed_point: 1.,
        error_tolerance: 0.15,
        ..Default::default()
    };
    let report = coder.encode_raw(&[1.25, 2.4, 3.4], &cfg).unwrap();
    match failed(&report) {
        NumpressRejection::Accuracy(f) => {
            assert_eq!(f.index, 1);
            assert_eq!(f.original, 2.4);
            assert_eq!(f.decoded, 2.);
        }
        other => panic!("{other:?}"),
    }
    cfg.error_tolerance = 0.25;
    assert!(coder.encode_raw(&[1.5], &cfg).unwrap().is_encoded()); // |1-original/decoded|=.25, inverse ratio differs
    cfg.error_tolerance = f64::from_bits(0.25f64.to_bits() - 1);
    assert!(matches!(
        failed(&coder.encode_raw(&[1.5], &cfg).unwrap()),
        NumpressRejection::Accuracy(_)
    ));
    cfg.error_tolerance = 0.125;
    assert!(coder.encode_raw(&[0., 0.125], &cfg).unwrap().is_encoded());
    cfg.error_tolerance = 0.124;
    match failed(&coder.encode_raw(&[0., 0.125], &cfg).unwrap()) {
        NumpressRejection::Accuracy(f) => assert_eq!(f.index, 1),
        other => panic!("{other:?}"),
    }
}
#[test]
fn source_tolerance_gating_and_pic_magnitude_independence() {
    let coder = MSNumpressCoder::default();
    for tolerance in [0., -1., f64::NEG_INFINITY, f64::NAN] {
        let cfg = NumpressConfig {
            compression: NumpressCompression::Linear,
            estimate_fixed_point: false,
            fixed_point: 1.,
            error_tolerance: tolerance,
            ..Default::default()
        };
        assert!(coder.encode_raw(&[1.4], &cfg).unwrap().is_encoded());
    }
    let cfg = NumpressConfig {
        error_tolerance: f64::INFINITY,
        estimate_fixed_point: false,
        fixed_point: 1.,
        ..config(NumpressCompression::Linear)
    };
    assert!(coder.encode_raw(&[1.4], &cfg).unwrap().is_encoded());
    let cfg = NumpressConfig {
        error_tolerance: f64::MIN_POSITIVE,
        fixed_point: f64::NAN,
        linear_fp_mass_acc: f64::NAN,
        ..config(NumpressCompression::Pic)
    };
    assert!(
        coder
            .encode_raw(&[-0.5, 0.5, 1.5], &cfg)
            .unwrap()
            .is_encoded()
    );
}
#[test]
fn raw_errors_are_reported_without_silent_fallback_and_unused_config_is_ignored() {
    let coder = MSNumpressCoder::default();
    for (mode, input) in [
        (NumpressCompression::Pic, vec![-1.]),
        (NumpressCompression::Slof, vec![-2.]),
        (NumpressCompression::Linear, vec![f64::NAN]),
    ] {
        let report = coder.encode(&input, false, &config(mode)).unwrap();
        assert!(matches!(failed(&report), NumpressRejection::Codec(_)));
        assert!(report.output.is_empty());
    }
    let cfg = NumpressConfig {
        fixed_point: f64::NAN,
        linear_fp_mass_acc: f64::NAN,
        ..config(NumpressCompression::Linear)
    };
    assert!(coder.encode_raw(&literal(), &cfg).unwrap().is_encoded()); // both ignored by maximal estimation path
    let cfg = NumpressConfig {
        estimate_fixed_point: false,
        ..cfg
    };
    assert!(matches!(
        failed(&coder.encode_raw(&literal(), &cfg).unwrap()),
        NumpressRejection::Codec(_)
    ));
}
#[test]
fn float_overload_promotes_before_encoding_without_inventing_f32_decoder() {
    let coder = MSNumpressCoder::default();
    let input = literal().map(|value| value as f32);
    let promoted: Vec<_> = input.iter().map(|&x| f64::from(x)).collect();
    for mode in [
        NumpressCompression::Linear,
        NumpressCompression::Pic,
        NumpressCompression::Slof,
    ] {
        assert_eq!(
            coder.encode_f32(&input, false, &config(mode)).unwrap(),
            coder.encode(&promoted, false, &config(mode)).unwrap()
        );
    }
}
fn compressed(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}
#[test]
fn zlib_transport_is_standard_and_independent_of_source_byte_literals() {
    let coder = MSNumpressCoder::default();
    // Independently compressed by Python zlib, not the Rust writer. These are
    // derived transport fixtures; only the uncompressed strings are source goldens.
    for row in include_str!("data/numpress_coder_transport.tsv")
        .lines()
        .skip(1)
    {
        let fields: Vec<_> = row.split('\t').collect();
        let cfg = config(fields[0].parse().unwrap());
        assert_eq!(
            coder.decode(fields[2], true, &cfg).unwrap(),
            coder.decode(fields[1], false, &cfg).unwrap()
        );
    }
    for mode in [
        NumpressCompression::Linear,
        NumpressCompression::Pic,
        NumpressCompression::Slof,
    ] {
        let cfg = config(mode);
        let raw = coder.encode_raw(&literal(), &cfg).unwrap();
        let zipped = STANDARD.encode(compressed(&raw.output));
        assert_eq!(
            coder.decode(&zipped, true, &cfg).unwrap(),
            coder.decode_raw(&raw.output, &cfg).unwrap()
        );
        let report = coder.encode(&literal(), true, &cfg).unwrap();
        assert!(report.is_encoded());
        assert_eq!(
            coder.decode(&report.output, true, &cfg).unwrap(),
            coder
                .decode(
                    &coder.encode(&literal(), false, &cfg).unwrap().output,
                    false,
                    &cfg
                )
                .unwrap()
        );
    }
    // Multiple chunks, including encoder flush and decoder growth paths.
    let input: Vec<_> = (0..20_000).map(|i| f64::from((i * 97) % 104729)).collect();
    let cfg = config(NumpressCompression::Pic);
    let encoded = coder.encode(&input, true, &cfg).unwrap();
    assert!(encoded.is_encoded());
    assert_eq!(coder.decode(&encoded.output, true, &cfg).unwrap(), input);
}
#[test]
fn short_source_base64_semantics_and_strict_long_form_validation() {
    let coder = MSNumpressCoder::default();
    let cfg = config(NumpressCompression::Pic);
    for short in ["", "?", "??", "???", "λ"] {
        for zlib in [false, true] {
            assert!(coder.decode(short, zlib, &cfg).unwrap().is_empty());
        }
    }
    for malformed in [
        "!!!!", "YQ=", "YQ===", "Y Q=", "YR==", "YQ==\n", "_AAA", "AAAA====",
    ] {
        if malformed.len() < 4 {
            continue;
        }
        assert!(
            coder.decode(malformed, false, &cfg).is_err(),
            "{malformed:?}"
        );
    }
    assert!(
        coder
            .decode_raw(&[], &config(NumpressCompression::Linear))
            .unwrap()
            .is_empty()
    ); // source wrapper empty short-circuit
    assert!(numpress::decode_linear(&[]).is_err()); // raw codec distinction retained
    let mut destination = vec![42.];
    assert!(
        coder
            .decode_into("!!!!", &mut destination, false, &cfg)
            .is_err()
    );
    assert_eq!(destination, [42.]);
}
#[test]
fn zlib_corruption_trailing_members_and_expansion_are_bounded() {
    let coder = MSNumpressCoder::default();
    let cfg = config(NumpressCompression::Pic);
    let valid = compressed(&[0x88; 32]);
    for input in [
        valid[..valid.len() - 1].to_vec(),
        [valid.clone(), vec![0]].concat(),
        [valid.clone(), valid.clone()].concat(),
    ] {
        assert!(coder.decode(&STANDARD.encode(input), true, &cfg).is_err());
    }
    let mut corrupted = valid.clone();
    let n = corrupted.len();
    corrupted[n - 1] ^= 1;
    assert!(
        coder
            .decode(&STANDARD.encode(corrupted), true, &cfg)
            .is_err()
    );
    let bomb = STANDARD.encode(compressed(&vec![0x88; 100_000]));
    let bounded = MSNumpressCoder {
        limits: NumpressCoderLimits {
            raw: numpress::NumpressLimits {
                max_encoded_bytes: 1000,
                ..Default::default()
            },
            ..Default::default()
        },
    };
    assert!(bounded.decode(&bomb, true, &cfg).is_err());
    assert!(
        bounded
            .decode(&bomb, true, &NumpressConfig::default())
            .is_err()
    ); // transport runs even for NONE
}
#[test]
fn cumulative_raw_verification_promotion_and_transport_budgets_do_not_reset() {
    let mut coder = MSNumpressCoder {
        limits: NumpressCoderLimits {
            raw: numpress::NumpressLimits {
                max_work: 128,
                ..Default::default()
            },
            ..Default::default()
        },
    };
    // PIC raw encode64 + raw verify64 fit separately; validation12 exceeds shared128.
    let cfg = config(NumpressCompression::Pic);
    let mut destination = vec![42];
    assert!(
        coder
            .encode_raw_into(&[1.5], &mut destination, &cfg)
            .is_err()
    );
    assert_eq!(destination, [42]);
    coder.limits.raw.max_work = 140;
    assert!(coder.encode_raw(&[1.5], &cfg).unwrap().is_encoded());
    coder.limits.raw.max_work = 1_000_000;
    coder.limits.max_total_bytes = 12;
    assert!(coder.encode_raw(&[1.5], &cfg).is_err()); // raw cap5 + verification8
    coder.limits.max_total_bytes = 13;
    assert!(coder.encode_raw(&[1.5], &cfg).unwrap().is_encoded());
    assert!(coder.encode(&[1.5], false, &cfg).is_err());
    coder.limits.max_total_bytes = 17;
    assert!(coder.encode(&[1.5], false, &cfg).unwrap().is_encoded());
    assert!(coder.encode_f32(&[1.5], false, &cfg).is_err());
    coder.limits.max_total_bytes = 0;
    assert_eq!(
        coder
            .encode(&[f64::NAN], true, &NumpressConfig::default())
            .unwrap()
            .status,
        NumpressEncodeStatus::Disabled
    );
    assert!(
        coder
            .encode_f32(&[1.], false, &NumpressConfig::default())
            .is_err()
    );
}
#[test]
fn decoded_base64_padding_uses_exact_binary_length_and_output_limits() {
    let header = STANDARD.encode(1f64.to_be_bytes());
    let coder = MSNumpressCoder {
        limits: NumpressCoderLimits {
            raw: numpress::NumpressLimits {
                max_encoded_bytes: 8,
                ..Default::default()
            },
            ..Default::default()
        },
    };
    assert!(
        coder
            .decode(&header, false, &config(NumpressCompression::Linear))
            .unwrap()
            .is_empty()
    );
    let coder = MSNumpressCoder {
        limits: NumpressCoderLimits {
            raw: numpress::NumpressLimits {
                max_values: 1,
                ..Default::default()
            },
            ..Default::default()
        },
    };
    assert!(
        coder
            .decode(
                &STANDARD.encode([0x88]),
                false,
                &config(NumpressCompression::Pic)
            )
            .is_err()
    );
}
#[test]
fn source_hundred_value_series_lengths_and_error_bounds_retain_wrapper_policy() {
    let coder = MSNumpressCoder::default();
    let input: Vec<_> = (0..100)
        .map(|i| 400. + f64::from(i) + f64::from(i) * 10f64.powi(-(i % 10) - 1))
        .collect();
    for (mode, length, abs, rel) in [
        (NumpressCompression::Linear, 360, 1e-5, 1e-7),
        (NumpressCompression::Pic, 268, 0.99, 0.01),
        (NumpressCompression::Slof, 280, 0.05, 1e-4),
    ] {
        let cfg = config(mode);
        let report = coder.encode(&input, false, &cfg).unwrap();
        assert!(report.is_encoded());
        assert_eq!(report.output.len(), length);
        for (actual, expected) in coder
            .decode(&report.output, false, &cfg)
            .unwrap()
            .iter()
            .zip(&input)
        {
            assert!((actual - expected).abs() <= abs);
            assert!((actual / expected).max(expected / actual) - 1. <= rel);
        }
    }
}
