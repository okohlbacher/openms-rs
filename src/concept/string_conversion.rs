// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The text-to-number conversions of `StringUtils`
//! (`DATASTRUCTURES/StringUtils.h`, `source/DATASTRUCTURES/StringUtils.cpp` at
//! core SDK `bc9cc12`): `toInt32`, `toInt64`, `toFloat`, `toDouble`,
//! `extractDouble`, `extractInt` and `skipWhitespace`, in the one place every
//! layer can reach.
//!
//! This is the single implementation behind
//! [`crate::data_structures::string_utils`]'s public conversions (which add
//! the crate error type) and behind the low layers (`metadata`, `concept`)
//! that cannot depend on `data_structures` without closing a cycle. The errors
//! here are the source's `ConversionError` message texts.
//!
//! Evidence: `../oracle/integ-A-stringutils` (`to_double_and_to_int32_match_the_release_build`
//! below) and `../oracle/B8a-stringutils` (`tests/string_utils.rs`), both the
//! Linux x86_64 Release install on ibminode06.

/// The four characters every whitespace rule of the header uses: space, tab,
/// line feed and carriage return.
const fn is_blank(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r')
}

/// Index of the first byte at or after `index` that is not one of the four
/// blanks.
pub(crate) fn skip_blanks(bytes: &[u8], mut index: usize) -> usize {
    while bytes.get(index).copied().is_some_and(is_blank) {
        index += 1;
    }
    index
}

/// Length of the case-insensitive ASCII `word` at `bytes[index..]`, if present.
fn word_at(bytes: &[u8], index: usize, word: &[u8]) -> Option<usize> {
    let end = index.checked_add(word.len())?;
    let candidate = bytes.get(index..end)?;
    candidate.eq_ignore_ascii_case(word).then_some(word.len())
}

/// The anonymous `tryParseNaN` (`StringUtils.cpp:35-58`): `nan` in any case at
/// `index`, optionally followed by `(`, anything but `)`, and `)`. Returns the
/// end of the match; an unclosed `(` is no match at all.
fn try_parse_nan(bytes: &[u8], index: usize) -> Option<usize> {
    word_at(bytes, index, b"nan")?;
    let mut end = index + 3;
    if bytes.get(end) == Some(&b'(') {
        let close = bytes.get(end..)?.iter().position(|b| *b == b')')?;
        end += close + 1;
    }
    Some(end)
}

/// What a `std::from_chars` general-format literal denotes.
#[derive(Clone, Copy)]
enum Literal {
    /// Decimal digits; `nonzero` when a mantissa digit is not `0`.
    Decimal {
        nonzero: bool,
    },
    Infinity,
    NotANumber,
}

/// The extent of the `std::from_chars(first, last, value)` general-format
/// literal at `start`: an optional `-`, then `inf` or `infinity`, `nan` with an
/// optional `(n-char-sequence)` of letters, digits and `_`, or digits with an
/// optional fraction and an optional exponent, all case-insensitive. A leading
/// `+`, whitespace and hexadecimal are not part of the grammar. Returns the
/// end, whether the literal is negative, and what it denotes.
fn scan_from_chars(bytes: &[u8], start: usize) -> Option<(usize, bool, Literal)> {
    let mut cursor = start;
    let negative = bytes.get(cursor) == Some(&b'-');
    if negative {
        cursor += 1;
    }
    if let Some(length) =
        word_at(bytes, cursor, b"infinity").or_else(|| word_at(bytes, cursor, b"inf"))
    {
        return Some((cursor + length, negative, Literal::Infinity));
    }
    if word_at(bytes, cursor, b"nan").is_some() {
        let mut end = cursor + 3;
        if bytes.get(end) == Some(&b'(') {
            let mut close = end + 1;
            while bytes
                .get(close)
                .is_some_and(|b| b.is_ascii_alphanumeric() || *b == b'_')
            {
                close += 1;
            }
            if bytes.get(close) == Some(&b')') {
                end = close + 1;
            }
        }
        return Some((end, negative, Literal::NotANumber));
    }
    let mut nonzero = false;
    let mut digits = 0usize;
    while let Some(b) = bytes.get(cursor).filter(|b| b.is_ascii_digit()) {
        nonzero |= *b != b'0';
        digits += 1;
        cursor += 1;
    }
    if bytes.get(cursor) == Some(&b'.') {
        let mut fraction = cursor + 1;
        let mut fraction_digits = 0usize;
        while let Some(b) = bytes.get(fraction).filter(|b| b.is_ascii_digit()) {
            nonzero |= *b != b'0';
            fraction_digits += 1;
            fraction += 1;
        }
        if digits + fraction_digits > 0 {
            digits += fraction_digits;
            cursor = fraction;
        }
    }
    if digits == 0 {
        return None;
    }
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
    Some((cursor, negative, Literal::Decimal { nonzero }))
}

/// A floating-point type `std::from_chars` converts to.
trait FromCharsFloat: Copy + std::str::FromStr + std::ops::Neg<Output = Self> {
    const INFINITY: Self;
    const NAN: Self;
    fn is_infinite(self) -> bool;
    fn is_zero(self) -> bool;
}

impl FromCharsFloat for f64 {
    const INFINITY: Self = f64::INFINITY;
    const NAN: Self = f64::NAN;
    fn is_infinite(self) -> bool {
        f64::is_infinite(self)
    }
    fn is_zero(self) -> bool {
        self == 0.0
    }
}

impl FromCharsFloat for f32 {
    const INFINITY: Self = f32::INFINITY;
    const NAN: Self = f32::NAN;
    fn is_infinite(self) -> bool {
        f32::is_infinite(self)
    }
    fn is_zero(self) -> bool {
        self == 0.0
    }
}

/// `std::from_chars` of a floating-point value at `start` as libstdc++ does it
/// on the reference build: correctly rounded, and `result_out_of_range` (no
/// value) when a finite literal overflows to infinity or when a literal with a
/// nonzero digit rounds to zero. A subnormal result is a value. Returns the
/// value and the end of the literal.
fn from_chars_float<T: FromCharsFloat>(text: &str, start: usize) -> Option<(T, usize)> {
    let bytes = text.as_bytes();
    let (end, negative, literal) = scan_from_chars(bytes, start)?;
    let value = match literal {
        Literal::Infinity => {
            if negative {
                -T::INFINITY
            } else {
                T::INFINITY
            }
        }
        Literal::NotANumber => {
            if negative {
                -T::NAN
            } else {
                T::NAN
            }
        }
        Literal::Decimal { nonzero } => {
            let value: T = text.get(start..end)?.parse().ok()?;
            if value.is_infinite() || (nonzero && value.is_zero()) {
                return None;
            }
            value
        }
    };
    Some((value, end))
}

/// `std::from_chars` of an integer at `start`: an optional `-` and decimal
/// digits, out-of-range values refused.
fn from_chars_int<T: std::str::FromStr>(text: &str, start: usize) -> Option<(T, usize)> {
    let bytes = text.as_bytes();
    let mut cursor = start;
    if bytes.get(cursor) == Some(&b'-') {
        cursor += 1;
    }
    let digits = cursor;
    while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
        cursor += 1;
    }
    if cursor == digits {
        return None;
    }
    let value = text.get(start..cursor)?.parse().ok()?;
    Some((value, cursor))
}

/// The source's `ConversionError` for text with no convertible prefix.
fn not_converted(text: &str, kind: &str) -> String {
    format!("Could not convert string '{text}' to {kind} value")
}

/// The source's `ConversionError` for a converted prefix followed by other
/// characters; `after` is the zero-based index of the first of them.
fn trailing(text: &str, kind: &str, after: usize) -> String {
    format!(
        "Prefix of string '{text}' successfully converted to {kind} value. \
         Additional characters found at position {}",
        after + 1
    )
}

/// `StringUtilsHelper::toInt32` and `toInt64` (`StringUtils.cpp:136-200`).
fn parse_integer<T: std::str::FromStr>(text: &str, prefix_kind: &str) -> Result<T, String> {
    let bytes = text.as_bytes();
    let mut start = skip_blanks(bytes, 0);
    if start == bytes.len() {
        return Err(not_converted(text, "an integer"));
    }
    if bytes.get(start) == Some(&b'+') {
        start += 1;
    }
    let (value, end) =
        from_chars_int(text, start).ok_or_else(|| not_converted(text, "an integer"))?;
    let after = skip_blanks(bytes, end);
    if after != bytes.len() {
        return Err(trailing(text, prefix_kind, after));
    }
    Ok(value)
}

/// `StringUtilsHelper::toFloat` and `toDouble` (`StringUtils.cpp:202-276`).
fn parse_float<T: FromCharsFloat>(text: &str, kind: &str) -> Result<T, String> {
    let bytes = text.as_bytes();
    let first = skip_blanks(bytes, 0);
    if first == bytes.len() {
        return Err(not_converted(text, kind));
    }
    if let Some(end) = try_parse_nan(bytes, first) {
        if skip_blanks(bytes, end) == bytes.len() {
            return Ok(T::NAN);
        }
    }
    let mut start = first;
    if bytes.get(start) == Some(&b'+') {
        start += 1;
    }
    let (value, end) =
        from_chars_float::<T>(text, start).ok_or_else(|| not_converted(text, kind))?;
    let after = skip_blanks(bytes, end);
    if after != bytes.len() {
        return Err(trailing(text, kind, after));
    }
    Ok(value)
}

/// Source `StringUtils::toInt32` (`StringUtils.cpp:136-166`), returning the
/// source's `ConversionError` message as the error.
pub(crate) fn to_int32(text: &str) -> Result<i32, String> {
    parse_integer(text, "an int32")
}

/// Source `StringUtils::toInt64` (`StringUtils.cpp:168-200`).
pub(crate) fn to_int64(text: &str) -> Result<i64, String> {
    parse_integer(text, "an integer")
}

/// Source `StringUtils::toFloat` (`StringUtils.cpp:202-237`).
pub(crate) fn to_float(text: &str) -> Result<f32, String> {
    parse_float(text, "a float")
}

/// Source `StringUtils::toDouble` (`StringUtils.cpp:239-276`).
pub(crate) fn to_double(text: &str) -> Result<f64, String> {
    parse_float(text, "a double")
}

/// Source `StringUtils::extractDouble`: the value at the start of `text` and
/// the bytes it occupies.
pub(crate) fn extract_double(text: &str) -> Option<(f64, usize)> {
    let bytes = text.as_bytes();
    if let Some(end) = try_parse_nan(bytes, 0) {
        return Some((f64::NAN, end));
    }
    let start = usize::from(bytes.first() == Some(&b'+'));
    from_chars_float::<f64>(text, start)
}

/// Source `StringUtils::extractInt`.
pub(crate) fn extract_int(text: &str) -> Option<(i32, usize)> {
    let start = usize::from(text.as_bytes().first() == Some(&b'+'));
    from_chars_int::<i32>(text, start)
}

/// Source `StringUtils::skipWhitespace`: the number of leading blanks.
pub(crate) fn skip_whitespace(text: &str) -> usize {
    skip_blanks(text.as_bytes(), 0)
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
