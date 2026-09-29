// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `StringUtils::toInt32` and `StringUtils::toDouble`
//! (`DATASTRUCTURES/StringUtils.h`, `source/DATASTRUCTURES/StringUtils.cpp` at
//! core SDK `bc9cc12`) for the layers below `cli`.
//!
//! The `cli` context module carried the same two conversions privately; this
//! module exists so `metadata` and other low layers can use them without an
//! edge to `cli`. At the group-A integration `concept::version_info`'s own
//! `toInt32` copy was switched to this one, and `to_double` was corrected
//! against the Linux x86_64 Release build (`../oracle/integ-A-stringutils`):
//! it refuses a literal that underflows to zero and keeps the sign of `-nan`,
//! which the copy it was taken from did not.

/// The source's four whitespace characters, skipped from `index` on.
fn skip_whitespace(bytes: &[u8], mut index: usize) -> usize {
    while bytes
        .get(index)
        .is_some_and(|b| matches!(*b, b' ' | b'\t' | b'\n' | b'\r'))
    {
        index += 1;
    }
    index
}

/// Source `StringUtils::toInt32` (`StringUtils.cpp:136-166`), returning the
/// source's `ConversionError` message as the error.
///
/// Leading and trailing space, tab, newline and carriage return are skipped,
/// one `+` may precede the number, and the rest must be a complete decimal
/// `i32`. Like `std::from_chars` after the source strips `+`, a `-` may still
/// follow it.
pub(crate) fn to_int32(text: &str) -> std::result::Result<i32, String> {
    let bytes = text.as_bytes();
    let not_converted = || format!("Could not convert string '{text}' to an integer value");
    let mut cursor = skip_whitespace(bytes, 0);
    if cursor == bytes.len() {
        return Err(not_converted());
    }
    if bytes[cursor] == b'+' {
        cursor += 1;
    }
    let number_start = cursor;
    if bytes.get(cursor) == Some(&b'-') {
        cursor += 1;
    }
    let digits_start = cursor;
    while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
        cursor += 1;
    }
    if cursor == digits_start {
        return Err(not_converted());
    }
    let value = text
        .get(number_start..cursor)
        .and_then(|digits| digits.parse::<i32>().ok())
        .ok_or_else(not_converted)?;
    let after = skip_whitespace(bytes, cursor);
    if after != bytes.len() {
        return Err(format!(
            "Prefix of string '{text}' successfully converted to an int32 value. Additional characters found at position {}",
            after + 1
        ));
    }
    Ok(value)
}

/// Length of a case-insensitive ASCII `word` at `bytes[index..]`, if present.
fn word_at(bytes: &[u8], index: usize, word: &str) -> Option<usize> {
    let candidate = bytes.get(index..index.checked_add(word.len())?)?;
    candidate
        .eq_ignore_ascii_case(word.as_bytes())
        .then_some(word.len())
}

/// Source `StringUtils::toDouble` (`StringUtils.cpp:239-276`), returning the
/// source's `ConversionError` message as the error.
///
/// Whitespace is skipped on both sides and `nan`, optionally followed by a
/// parenthesised payload, is accepted before anything else. Then one `+` may
/// precede a libstdc++ `std::from_chars` general-format number: an optional
/// `-`, digits with an optional fraction and exponent, `inf`/`infinity`, or
/// `nan` with an optional `(n-char-sequence)`, whose sign is kept (`-nan` is
/// `0xfff8000000000000`). A literal whose value is out of range is an error, as
/// `result_out_of_range` is: one that overflows to an infinity, and one with a
/// non-zero digit that rounds to zero (`1e-400`, `2.4e-324`; `0e-400` is zero
/// and accepted). Hexadecimal floats are rejected, as `std::from_chars`'s
/// general format rejects them. Every rule is the Linux x86_64 Release build's,
/// compared case by case with `../oracle/integ-A-stringutils`
/// (`to_double_matches_the_release_build`).
pub(crate) fn to_double(text: &str) -> std::result::Result<f64, String> {
    let bytes = text.as_bytes();
    let not_converted = || format!("Could not convert string '{text}' to a double value");
    let first = skip_whitespace(bytes, 0);
    if first == bytes.len() {
        return Err(not_converted());
    }
    if let Some(length) = word_at(bytes, first, "nan") {
        let mut end = first + length;
        if bytes.get(end) == Some(&b'(') {
            match bytes[end..].iter().position(|b| *b == b')') {
                Some(close) => end += close + 1,
                None => end = first,
            }
        }
        if end != first && skip_whitespace(bytes, end) == bytes.len() {
            return Ok(f64::NAN);
        }
    }
    let mut start = first;
    if bytes[start] == b'+' {
        start += 1;
    }
    let mut cursor = start;
    if bytes.get(cursor) == Some(&b'-') {
        cursor += 1;
    }
    let negative = bytes.get(start) == Some(&b'-');
    let infinity = word_at(bytes, cursor, "infinity").or_else(|| word_at(bytes, cursor, "inf"));
    let nan = word_at(bytes, cursor, "nan");
    let mut nonzero_mantissa = false;
    if let Some(length) = infinity {
        cursor += length;
    } else if let Some(length) = nan {
        cursor += length;
        // `nan(n-char-sequence)`: letters, digits and `_` up to a `)`; without
        // the `)` only `nan` is consumed.
        if bytes.get(cursor) == Some(&b'(') {
            let mut close = cursor + 1;
            while bytes
                .get(close)
                .is_some_and(|b| b.is_ascii_alphanumeric() || *b == b'_')
            {
                close += 1;
            }
            if bytes.get(close) == Some(&b')') {
                cursor = close + 1;
            }
        }
    } else {
        let integer_start = cursor;
        while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
            cursor += 1;
        }
        let mut digits = cursor - integer_start;
        if bytes.get(cursor) == Some(&b'.') {
            let fraction_start = cursor + 1;
            let mut fraction_end = fraction_start;
            while bytes.get(fraction_end).is_some_and(u8::is_ascii_digit) {
                fraction_end += 1;
            }
            if digits > 0 || fraction_end > fraction_start {
                digits += fraction_end - fraction_start;
                cursor = fraction_end;
            }
        }
        if digits == 0 {
            return Err(not_converted());
        }
        nonzero_mantissa = bytes[integer_start..cursor]
            .iter()
            .any(|b| matches!(*b, b'1'..=b'9'));
        if matches!(bytes.get(cursor), Some(b'e' | b'E')) {
            let mut exponent = cursor + 1;
            if matches!(bytes.get(exponent), Some(b'+' | b'-')) {
                exponent += 1;
            }
            let exponent_digits = exponent;
            while bytes.get(exponent).is_some_and(u8::is_ascii_digit) {
                exponent += 1;
            }
            if exponent > exponent_digits {
                cursor = exponent;
            }
        }
    }
    let value: f64 = if nan.is_some() && infinity.is_none() {
        if negative { -f64::NAN } else { f64::NAN }
    } else {
        let literal = text.get(start..cursor).ok_or_else(not_converted)?;
        literal.parse().map_err(|_| not_converted())?
    };
    // `result_out_of_range`: an overflow to infinity, or a non-zero literal
    // that rounds to zero.
    if (value.is_infinite() && infinity.is_none()) || (value == 0.0 && nonzero_mantissa) {
        return Err(not_converted());
    }
    let after = skip_whitespace(bytes, cursor);
    if after != bytes.len() {
        return Err(format!(
            "Prefix of string '{text}' successfully converted to a double value. Additional characters found at position {}",
            after + 1
        ));
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn to_int32_follows_the_source_conversion() {
        assert_eq!(to_int32("42"), Ok(42));
        assert_eq!(to_int32(" \t+42\r\n"), Ok(42));
        assert_eq!(to_int32("+-7"), Ok(-7));
        assert_eq!(to_int32("2147483647"), Ok(i32::MAX));
        assert!(to_int32("2147483648").is_err());
        assert!(to_int32("").is_err());
        assert!(to_int32("  ").is_err());
        assert!(to_int32("12a").is_err());
        assert!(to_int32("++1").is_err());
    }

    #[test]
    fn to_double_follows_the_source_conversion() {
        assert_eq!(to_double(" 2.5 "), Ok(2.5));
        assert_eq!(to_double("+1e3"), Ok(1000.0));
        assert_eq!(to_double("-inf"), Ok(f64::NEG_INFINITY));
        assert!(to_double("nan").unwrap().is_nan());
        assert!(to_double("1e999").is_err());
        assert!(to_double("0x10").is_err());
        assert!(to_double("1.5x").is_err());
        assert!(to_double("").is_err());
    }

    /// The driver's inputs, in its order (`../oracle/integ-A-stringutils`,
    /// `driver.cpp`).
    const DOUBLE_INPUTS: [&str; 38] = [
        "1e-400",
        "-1e-400",
        "1e-320",
        "4.9e-324",
        "2.4e-324",
        "2.5e-324",
        "1e309",
        "-1e309",
        "1.7976931348623157e308",
        "nan",
        "-nan",
        "+nan",
        "NaN(123)",
        "-nan(1)",
        " nan ",
        "inf",
        "-inf",
        "+inf",
        "infinity",
        "-Infinity",
        "  12.5  ",
        "+-3",
        "-+3",
        "0x1p3",
        "1e",
        "1e+",
        ".5",
        "5.",
        ".",
        "-0",
        "-0.0",
        "1.5abc",
        "",
        "   ",
        "1e-400 ",
        "0e-400",
        "123456789012345678901234567890",
        "1_0",
    ];
    const INT_INPUTS: [&str; 14] = [
        "12",
        " +7 ",
        "+-7",
        "-7",
        "2147483647",
        "2147483648",
        "-2147483648",
        "-2147483649",
        "1.5",
        "",
        " ",
        "0x10",
        "12 3",
        "007",
    ];

    /// Tier 1: every `StringUtils::toDouble`/`toInt32` result, bits or the
    /// `ConversionError` text, as the Linux x86_64 Release build returned it
    /// (`../oracle/integ-A-stringutils`, ibminode06, install
    /// openms4-release-bc9cc12-c19e494-174b576, two runs in each of two
    /// invocations, identical; tests/data/string_utils_release_provenance.json).
    #[test]
    fn to_double_and_to_int32_match_the_release_build() {
        let text = std::fs::read_to_string(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/data/string_utils_release/oracle.tsv"),
        )
        .unwrap();
        let mut rows = 0;
        for line in text.lines() {
            let fields: Vec<&str> = line.splitn(4, '\t').collect();
            let (function, index, outcome, value) = (fields[0], fields[1], fields[2], fields[3]);
            let index: usize = index.parse().unwrap();
            let actual = match function {
                "toDouble" => match to_double(DOUBLE_INPUTS[index]) {
                    Ok(number) => format!("ok\t{:#018x}", number.to_bits()),
                    Err(message) => format!("err\t{message}"),
                },
                "toInt32" => match to_int32(INT_INPUTS[index]) {
                    Ok(number) => format!("ok\t{number}"),
                    Err(message) => format!("err\t{message}"),
                },
                other => panic!("unknown function {other}"),
            };
            assert_eq!(
                actual,
                format!("{outcome}\t{value}"),
                "{function} case {index}"
            );
            rows += 1;
        }
        assert_eq!(rows, DOUBLE_INPUTS.len() + INT_INPUTS.len());
    }
}
