// Copyright (c) 2026 OpenMS Rust contributors
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::format::numpress::*;

fn bytes(hex: &str) -> Vec<u8> {
    hex.as_bytes()
        .chunks_exact(2)
        .map(|b| u8::from_str_radix(std::str::from_utf8(b).unwrap(), 16).unwrap())
        .collect()
}
fn values(hex: &str) -> Vec<f64> {
    if hex.is_empty() {
        return vec![];
    }
    hex.split(',')
        .map(|v| f64::from_bits(u64::from_str_radix(v, 16).unwrap()))
        .collect()
}
fn encode(mode: &str, data: &[f64], fp: f64) -> Vec<u8> {
    match mode {
        "linear" => encode_linear(data, fp),
        "pic" => encode_pic(data),
        "slof" => encode_slof(data, fp),
        "safe" => encode_safe(data),
        _ => unreachable!(),
    }
    .unwrap()
}
fn decode(mode: &str, data: &[u8]) -> Vec<f64> {
    match mode {
        "linear" => decode_linear(data),
        "pic" => decode_pic(data),
        "slof" => decode_slof(data),
        "safe" => decode_safe(data),
        _ => unreachable!(),
    }
    .unwrap()
}
fn close(actual: f64, expected: f64, mode: &str) {
    if mode == "slof" {
        assert!(
            (actual - expected).abs() <= 2e-14 * expected.abs().max(1.0),
            "{actual} vs {expected}"
        );
    } else {
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "{mode}: {actual} vs {expected}"
        );
    }
}

#[test]
fn all_three_literal_upstream_byte_fixtures_match_without_transport_wrapper() {
    let input = [100., 200., 300.00005, 400.00010];
    for row in include_str!("data/numpress_source_bytes.tsv")
        .lines()
        .skip(1)
    {
        let fields: Vec<_> = row.split('\t').collect();
        let expected = bytes(fields[2]);
        let fp = match fields[0] {
            "linear" => optimal_linear_fixed_point(&input).unwrap(),
            "slof" => optimal_slof_fixed_point(&input).unwrap(),
            _ => 0.0,
        };
        assert_eq!(encode(fields[0], &input, fp), expected);
        let output = decode(fields[0], &expected);
        assert_eq!(output.len(), 4);
        for (actual, original) in output.iter().zip(input) {
            let tolerance = if fields[0] == "slof" {
                original * 1e-4
            } else {
                0.001
            };
            assert!((actual - original).abs() <= tolerance);
        }
    }
}

#[test]
fn executed_unmodified_cpp_probe_matches_all_raw_bytes_and_decoded_values() {
    let mut cases = 0;
    for row in include_str!("data/numpress_cpp_differential.tsv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = row.split('\t').collect();
        assert_eq!(f.len(), 9);
        let fp = f64::from_bits(u64::from_str_radix(f[2], 16).unwrap());
        let input = values(f[3]);
        let expected = bytes(f[4]);
        assert_eq!(encode(f[1], &input, fp), expected, "case {}", f[0]);
        let output = decode(f[1], &expected);
        let target = values(f[5]);
        assert_eq!(output.len(), target.len());
        for (actual, expected) in output.into_iter().zip(target) {
            close(actual, expected, f[1]);
        }
        let helpers = [
            optimal_linear_fixed_point(&input),
            optimal_linear_fixed_point_mass(&input, 0.001),
            optimal_slof_fixed_point(&input),
        ];
        for (result, literal) in helpers.into_iter().zip(&f[6..]) {
            let expected = f64::from_bits(u64::from_str_radix(literal, 16).unwrap());
            if expected.is_finite() {
                assert_eq!(
                    result.unwrap().to_bits(),
                    expected.to_bits(),
                    "helper case {}",
                    f[0]
                );
            } else {
                assert!(result.is_err());
            }
        }
        cases += 1;
    }
    assert_eq!(cases, 295);
}

// Independent textual hexadecimal oracle: remove high zero/f digits, then
// reverse the retained hex characters. Does not reuse the bitwise production loop.
fn integer_oracle(value: u32) -> Vec<u8> {
    let hex = format!("{value:08x}");
    let leading = if hex.starts_with('0') {
        hex.chars().take_while(|c| *c == '0').count()
    } else if hex.starts_with('f') {
        hex.chars().take_while(|c| *c == 'f').count().min(7)
    } else {
        0
    };
    let head = leading as u8 + if hex.starts_with('f') { 8 } else { 0 };
    let mut digits = vec![head];
    digits.extend(
        hex[leading..]
            .chars()
            .rev()
            .map(|c| c.to_digit(16).unwrap() as u8),
    );
    if digits.len() % 2 != 0 {
        digits.push(0);
    }
    digits.chunks_exact(2).map(|p| p[0] * 16 + p[1]).collect()
}
#[test]
fn leading_nibble_patterns_and_terminal_padding_match_independent_oracles() {
    for (raw, expected) in [
        (vec![0x80], vec![0.]),
        (vec![0x88], vec![0., 0.]),
        (vec![0xff], vec![4294967295.]),
        (vec![0x72], vec![2.]),
        (vec![0x67, 0x10], vec![23.]),
        (vec![0x5f, 0xf7], vec![2047.]),
    ] {
        assert_eq!(decode_pic(&raw).unwrap(), expected);
    }
    let mut state = 0x1257u32;
    for value in [
        0,
        1,
        15,
        16,
        255,
        256,
        4095,
        4096,
        0xfffffff,
        0x10000000,
        0x7ffffffe,
        0x80000000,
        0xf0000000,
        0xfffffff0,
        u32::MAX,
    ]
    .into_iter()
    .chain((0..1000).map(|_| {
        state = state.wrapping_mul(1664525).wrapping_add(1013904223);
        state
    })) {
        let expected = integer_oracle(value);
        assert_eq!(decode_pic(&expected).unwrap(), [f64::from(value)]);
        if value < i32::MAX as u32 {
            assert_eq!(encode_pic(&[f64::from(value)]).unwrap(), expected);
        }
    }
    // Source accepts complete nonminimal encodings; it rejects incomplete ones.
    assert_eq!(decode_pic(&[0, 0, 0, 0, 0]).unwrap(), [0.]);
    for malformed in [&[0][..], &[0x10][..], &[0x80, 0][..], &[0xf0, 0][..]] {
        assert!(decode_pic(malformed).is_err(), "{malformed:x?}");
    }
}

#[test]
fn source_optimal_fixed_points_and_short_accuracy_shortcut() {
    assert_eq!(optimal_linear_fixed_point(&[]).unwrap(), 0.);
    assert_eq!(optimal_linear_fixed_point(&[100.]).unwrap(), 21474836.);
    assert_eq!(
        optimal_linear_fixed_point(&[100., 200.]).unwrap(),
        10737418.
    );
    assert_eq!(
        optimal_linear_fixed_point(&[0., 0., 0.]).unwrap(),
        2147483647.
    );
    assert!(optimal_linear_fixed_point(&[0.]).is_err());
    assert!(optimal_linear_fixed_point(&[0., 0.]).is_err());
    assert_eq!(optimal_linear_fixed_point(&[-1.]).unwrap(), -2147483647.);
    for n in 0..=2 {
        assert_eq!(
            optimal_linear_fixed_point_mass(&[f64::NAN; 2][..n], f64::NAN).unwrap(),
            0.
        );
    }
    let input = [100., 200., 300.];
    assert_eq!(
        optimal_linear_fixed_point_mass(&input, 0.001).unwrap(),
        500.
    );
    assert_eq!(optimal_linear_fixed_point_mass(&input, 1e-12).unwrap(), -1.);
    assert_eq!(optimal_linear_fixed_point_mass(&input, 0.).unwrap(), -1.);
    assert_eq!(
        optimal_linear_fixed_point_mass(&input, -0.001).unwrap(),
        -500.
    );
    assert_eq!(optimal_slof_fixed_point(&[]).unwrap(), 0.);
    assert_eq!(optimal_slof_fixed_point(&[0., -1., -2.]).unwrap(), 65534.);
    assert!(encode_slof(&[-2.], 65534.).is_err());
}

#[test]
fn linear_endianness_signed_residuals_and_initial_low_word_behavior() {
    let first = encode_linear(&[305419896., 305419897.], 1.).unwrap();
    assert_eq!(&first[..8], &1f64.to_be_bytes());
    assert_eq!(&first[8..12], &[0x78, 0x56, 0x34, 0x12]);
    assert_eq!(&first[12..], &[0x79, 0x56, 0x34, 0x12]);
    for residual in [i32::MIN, -100, -1, 0, 1, 100, i32::MAX] {
        let input = [0., 0., f64::from(residual) - 0.5];
        let encoded = encode_linear(&input, 1.).unwrap();
        assert_eq!(&encoded[16..], integer_oracle(residual as u32));
        assert_eq!(
            decode_linear(&encoded).unwrap(),
            [0., 0., f64::from(residual)]
        );
    }
    let encoded = encode_linear(&[4294967296., 4294967297.], 1.).unwrap();
    assert_eq!(decode_linear(&encoded).unwrap(), [0., 1.]);
    let encoded = encode_linear(&[-2.], 1.).unwrap();
    assert_eq!(decode_linear(&encoded).unwrap(), [4294967295.]);
    assert!(encode_linear(&[0., 0., 2147483648.], 1.).is_err());
    assert!(encode_linear(&[9223372036854775808.], 1.).is_err());
    assert!(encode_linear(&[0., -9223372036854775808., 0.], 1.).is_err());
}

#[test]
fn empty_cases_and_unused_fixed_point_payloads_preserve_source_semantics() {
    for fp in [0., -0., -1., f64::NAN, f64::INFINITY] {
        for mode in ["linear", "slof"] {
            let encoded = encode(mode, &[], fp);
            assert_eq!(encoded, fp.to_be_bytes());
            assert!(decode(mode, &encoded).is_empty());
        }
    }
    assert!(encode_pic(&[]).unwrap().is_empty());
    assert!(decode_pic(&[]).unwrap().is_empty());
    assert!(encode_safe(&[]).unwrap().is_empty());
    assert!(decode_safe(&[]).unwrap().is_empty()); // checked extension of undefined C++ empty decoder
}

#[test]
fn pic_uses_source_guard_not_the_broader_header_claim() {
    let input = [-0.5, -0.1, 0., 0.499, 0.5, 1.5, 2147483646.5];
    assert_eq!(
        decode_pic(&encode_pic(&input).unwrap()).unwrap(),
        [0., 0., 0., 0., 1., 2., 2147483647.]
    );
    for value in [
        -0.500001,
        2147483647.,
        4294967294.,
        f64::MAX,
        f64::NAN,
        f64::INFINITY,
    ] {
        assert!(encode_pic(&[value]).is_err());
    }
    assert_eq!(decode_pic(&[0xff]).unwrap(), [4294967295.]);
}

#[test]
fn slof_literal_math_signed_factors_and_cast_boundaries() {
    let encoded = encode_slof(&[0., 1., 2.], 100.).unwrap();
    assert_eq!(&encoded[..8], &100f64.to_be_bytes());
    assert_eq!(&encoded[8..], &[0, 0, 69, 0, 110, 0]);
    assert_eq!(encode_slof(&[-0.1], 1.).unwrap()[8..], [0, 0]);
    assert!(encode_slof(&[-0.9], 1.).is_err());
    assert!(encode_slof(&[-1.], 1.).is_err());
    assert!(encode_slof(&[f64::MAX], 1000.).is_err());
    assert!(decode_slof(&[vec![0; 8], vec![1, 0]].concat()).is_err());
    let negative_zero = [(-0f64).to_be_bytes().to_vec(), vec![1, 0]].concat();
    assert_eq!(decode_slof(&negative_zero).unwrap(), [-1.]);
    assert!(decode_slof(&[(-0f64).to_be_bytes().to_vec(), vec![0, 0]].concat()).is_err());
    let negative = encode_slof(&[1., 1.5, 2.], -1.).unwrap();
    assert_eq!(decode_slof(&negative).unwrap(), [0., 0., 0.]);
    // ln(1+x), rather than ln_1p, rounds this input's logarithm to zero.
    assert_eq!(encode_slof(&[1e-20], 1e20).unwrap()[8..], [0, 0]);
}

#[test]
fn safe_codec_order_and_non_bit_exact_source_boundary() {
    let encoded = encode_safe(&[1., 2., 4.]).unwrap();
    assert_eq!(
        encoded,
        [1f64.to_be_bytes(), 2f64.to_be_bytes(), 1f64.to_be_bytes()].concat()
    );
    assert_eq!(decode_safe(&encoded).unwrap(), [1., 2., 4.]);
    let input = [0.1, 1., 1., 0.1];
    let decoded = decode_safe(&encode_safe(&input).unwrap()).unwrap();
    assert_ne!(decoded[3].to_bits(), input[3].to_bits());
    assert!((decoded[3] - input[3]).abs() < 1e-15);
    assert!(encode_safe(&[f64::MAX, -f64::MAX, 0.]).is_err());
    assert!(decode_safe(&[0; 7]).is_err());
}

#[test]
fn malformed_streams_and_nonfinite_selected_values_return_errors() {
    for length in (0..8).chain(9..12).chain(13..16) {
        assert!(decode_linear(&vec![0; length]).is_err());
    }
    for length in (0..8).chain([9, 11, 13, 15, 17]) {
        assert!(decode_slof(&vec![0; length]).is_err());
    }
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for encode in [encode_linear, encode_slof] {
            assert!(encode(&[1.], bad).is_err());
            assert!(encode(&[bad], 1.).is_err());
        }
        assert!(encode_safe(&[bad]).is_err());
        assert!(decode_safe(&bad.to_be_bytes()).is_err());
        assert!(decode_linear(&[bad.to_be_bytes().to_vec(), vec![0; 4]].concat()).is_err());
        assert!(decode_slof(&[bad.to_be_bytes().to_vec(), vec![0; 2]].concat()).is_err());
    }
    let mut linear = encode_linear(&[1., 2.], 1.).unwrap();
    linear.push(0); // incomplete residual header
    assert!(decode_linear(&linear).is_err());
    let mut extreme = f64::MIN_POSITIVE.to_be_bytes().to_vec();
    extreme.extend_from_slice(&65535u16.to_le_bytes());
    assert!(decode_slof(&extreme).is_err());
}

#[test]
fn exact_output_count_preflight_and_resource_errors_expose_no_partial_values() {
    let one = encode_linear(&[1.], 1.).unwrap();
    let limits = NumpressLimits {
        max_values: 1,
        ..Default::default()
    };
    assert_eq!(decode_linear_with_limits(&one, &limits).unwrap(), [1.]);
    assert!(decode_pic_with_limits(&[0x88], &limits).is_err());
    let before = vec![42.];
    let mut caller = before.clone();
    if let Ok(result) = decode_pic_with_limits(&[0x88], &limits) {
        caller = result;
    }
    assert_eq!(caller, before);
    let no_work = NumpressLimits {
        max_work: 0,
        ..Default::default()
    };
    assert!(encode_pic_with_limits(&[1.], &no_work).is_err());
    assert!(decode_linear_with_limits(&one, &no_work).is_err());
    assert!(optimal_linear_fixed_point_with_limits(&[1.], &no_work).is_err());
    let no_bytes = NumpressLimits {
        max_encoded_bytes: 0,
        ..Default::default()
    };
    assert!(encode_linear_with_limits(&[], 1., &no_bytes).is_err());
    assert!(decode_pic_with_limits(&[0x80], &no_bytes).is_err());
    assert!(encode_safe_with_limits(&[1., 2.], &limits).is_err());
}

#[test]
fn long_linear_prediction_overflow_is_checked_during_validation_pass() {
    // Independent nibble construction: 0x7fffffff is header zero followed by
    // seven low-order f digits and a final 7. Repeated positive residuals make
    // the second-order recurrence exceed i64 within the normal byte/work caps.
    let mut encoded = 1f64.to_be_bytes().to_vec();
    encoded.extend_from_slice(&[0; 8]);
    let mut high = None;
    for _ in 0..100_000 {
        for digit in [0, 15, 15, 15, 15, 15, 15, 15, 7] {
            if let Some(first) = high.take() {
                encoded.push((first << 4) | digit);
            } else {
                high = Some(digit);
            }
        }
    }
    if let Some(first) = high {
        encoded.push(first << 4);
    }
    let error = decode_linear(&encoded).unwrap_err();
    assert!(error.to_string().contains("prediction overflow"));
}

#[test]
fn source_large_series_quantization_bounds_and_encoded_length_envelopes() {
    let input: Vec<_> = (0..100)
        .map(|i| 400.0 + f64::from(i) + f64::from(i) * 10f64.powi(-(i % 10) - 1))
        .collect();
    for (mode, base64_len, absolute, relative) in [
        ("linear", 360, 1e-5, 1e-7),
        ("pic", 268, 0.99, 0.01),
        ("slof", 280, 0.05, 1e-4),
    ] {
        let fp = match mode {
            "linear" => optimal_linear_fixed_point(&input).unwrap(),
            "slof" => optimal_slof_fixed_point(&input).unwrap(),
            _ => 0.,
        };
        let encoded = encode(mode, &input, fp);
        // The source only asserts encoded base64 length; it does not pin all raw bytes.
        assert_eq!(encoded.len().div_ceil(3) * 4, base64_len);
        for (actual, expected) in decode(mode, &encoded).into_iter().zip(&input) {
            assert!((actual - expected).abs() <= absolute);
            assert!((actual / expected).max(expected / actual) - 1. <= relative);
        }
    }
}
