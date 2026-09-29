// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Source `StringUtils` numeric text conversion shared by the text formats and
//! the command line (`DATASTRUCTURES/StringUtils.h`, `StringUtils.cpp`).
//!
//! These were first written for the command line (`cli/context.rs`) and moved
//! here, below both, so that the text-format adapters convert exactly as the
//! source's `StringUtils::toDouble` / `toInt32` do without a second copy.

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
/// precede a `std::from_chars` general-format number: an optional `-`, digits
/// with an optional fraction and exponent, or `inf`/`infinity`. A finite
/// literal that overflows is an error, as `result_out_of_range` is; underflow
/// rounds, as the oracle's libc++ fallback accepts it. Hexadecimal floats are
/// rejected, following `std::from_chars` rather than the `strtod` fallback the
/// source uses on libc++.
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
    let infinity = word_at(bytes, cursor, "infinity").or_else(|| word_at(bytes, cursor, "inf"));
    let nan = word_at(bytes, cursor, "nan");
    if let Some(length) = infinity.or(nan) {
        cursor += length;
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
    let literal = text.get(start..cursor).ok_or_else(not_converted)?;
    let value: f64 = if nan.is_some() && infinity.is_none() {
        f64::NAN
    } else {
        literal.parse().map_err(|_| not_converted())?
    };
    if value.is_infinite() && infinity.is_none() {
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

/// Source `StringUtils::toFloat` (`StringUtils.cpp:200-237`), returning the
/// source's `ConversionError` message as the error.
///
/// The accepted text is that of [`to_double`]; the number is then rounded once,
/// directly to the nearest `f32`, as `std::from_chars<float>` does, not through
/// a `double`. A finite literal beyond the `f32` range is an error, as
/// `result_out_of_range` is.
pub(crate) fn to_float(text: &str) -> std::result::Result<f32, String> {
    let not_converted = || format!("Could not convert string '{text}' to a float value");
    let wide = to_double(text).map_err(|message| message.replace("a double", "a float"))?;
    if wide.is_nan() {
        return Ok(f32::NAN);
    }
    let literal = text.trim_matches(|c| matches!(c, ' ' | '\t' | '\n' | '\r'));
    let literal = literal.strip_prefix('+').unwrap_or(literal);
    let value: f32 = literal.parse().map_err(|_| not_converted())?;
    if value.is_infinite() && !wide.is_infinite() {
        return Err(not_converted());
    }
    Ok(value)
}

/// The characters C's `isspace` accepts in the "C" locale, which `strtod` and
/// `strtoul` skip before a number.
fn c_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}

/// `std::stod`: `strtod` on the text, which skips leading C whitespace, reads
/// the longest numeric prefix and ignores whatever follows it.
///
/// The prefix is an optional sign and then a decimal number (digits with an
/// optional fraction, at least one digit, and an exponent only when it has
/// digits), `inf`, `infinity` or `nan` (any case), or `nan(...)`.
///
/// # Errors
///
/// Returns `std::stod`'s `std::invalid_argument` ("no conversion") when no
/// prefix converts, and `std::out_of_range` when the value overflows to an
/// infinity or underflows to zero or a subnormal (`ERANGE` in the C library).
/// A hexadecimal number (`0x...`), which `strtod` also reads, is refused with
/// its own message rather than read as its leading `0`.
pub(crate) fn stod(text: &str) -> std::result::Result<f64, String> {
    let bytes = text.as_bytes();
    let invalid = || format!("stod: no conversion of '{text}'");
    let mut start = 0;
    while bytes.get(start).copied().is_some_and(c_space) {
        start += 1;
    }
    let mut cursor = start;
    let negative = bytes.get(cursor) == Some(&b'-');
    if matches!(bytes.get(cursor), Some(b'+' | b'-')) {
        cursor += 1;
    }
    if word_at(bytes, cursor, "nan").is_some() {
        return Ok(f64::NAN);
    }
    if word_at(bytes, cursor, "inf").is_some() {
        return Ok(if negative {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        });
    }
    if bytes.get(cursor) == Some(&b'0')
        && matches!(bytes.get(cursor + 1), Some(b'x' | b'X'))
        && bytes
            .get(cursor + 2)
            .is_some_and(|b| b.is_ascii_hexdigit() || *b == b'.')
    {
        return Err(format!(
            "stod: hexadecimal floating-point text is not supported: '{text}'"
        ));
    }
    let integer_start = cursor;
    while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
        cursor += 1;
    }
    let mut digits = cursor - integer_start;
    if bytes.get(cursor) == Some(&b'.') {
        let mut end = cursor + 1;
        while bytes.get(end).is_some_and(u8::is_ascii_digit) {
            end += 1;
        }
        let fraction = end - cursor - 1;
        if digits > 0 || fraction > 0 {
            digits += fraction;
            cursor = end;
        }
    }
    if digits == 0 {
        return Err(invalid());
    }
    if matches!(bytes.get(cursor), Some(b'e' | b'E')) {
        let mut exponent = cursor + 1;
        if matches!(bytes.get(exponent), Some(b'+' | b'-')) {
            exponent += 1;
        }
        let digits_start = exponent;
        while bytes.get(exponent).is_some_and(u8::is_ascii_digit) {
            exponent += 1;
        }
        if exponent > digits_start {
            cursor = exponent;
        }
    }
    let literal = text.get(start..cursor).ok_or_else(invalid)?;
    let value: f64 = literal.parse().map_err(|_| invalid())?;
    let nonzero_digit = text[integer_start..cursor]
        .bytes()
        .take_while(|b| !matches!(b, b'e' | b'E'))
        .any(|b| matches!(b, b'1'..=b'9'));
    if value.is_infinite() || (nonzero_digit && value.abs() < f64::MIN_POSITIVE) {
        return Err(format!("stod: '{text}' is out of range"));
    }
    Ok(value)
}

/// `std::stoul` with the 64-bit `unsigned long` of the source platforms:
/// `strtoul` in base 10, which skips leading C whitespace, accepts one sign,
/// reads the longest digit prefix and ignores whatever follows. A `-` negates
/// the value in unsigned arithmetic, as C does (`-1` is `u64::MAX`).
///
/// # Errors
///
/// Returns `std::invalid_argument` ("no conversion") when no digit follows,
/// and `std::out_of_range` when the digits exceed `u64`.
pub(crate) fn stoul(text: &str) -> std::result::Result<u64, String> {
    let bytes = text.as_bytes();
    let mut cursor = 0;
    while bytes.get(cursor).copied().is_some_and(c_space) {
        cursor += 1;
    }
    let negative = bytes.get(cursor) == Some(&b'-');
    if matches!(bytes.get(cursor), Some(b'+' | b'-')) {
        cursor += 1;
    }
    let start = cursor;
    while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
        cursor += 1;
    }
    if cursor == start {
        return Err(format!("stoul: no conversion of '{text}'"));
    }
    let value: u64 = text[start..cursor]
        .parse()
        .map_err(|_| format!("stoul: '{text}' is out of range"))?;
    Ok(if negative {
        value.wrapping_neg()
    } else {
        value
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stod_reads_the_longest_prefix() {
        assert_eq!(stod(" 12.5 min"), Ok(12.5));
        assert_eq!(stod("\x0c-3e2x"), Ok(-300.0));
        assert_eq!(stod("7e"), Ok(7.0));
        assert_eq!(stod(".5"), Ok(0.5));
        assert_eq!(stod("5."), Ok(5.0));
        assert_eq!(stod("-inf"), Ok(f64::NEG_INFINITY));
        assert!(stod("NaN").unwrap().is_nan());
        assert!(stod("").is_err());
        assert!(stod(".").is_err());
        assert!(stod("abc").is_err());
        assert!(stod("1e999").is_err());
        assert!(stod("1e-400").is_err());
        assert_eq!(stod("0e-400"), Ok(0.0));
        assert!(stod("0x1p3").is_err());
    }

    #[test]
    fn stoul_follows_strtoul() {
        assert_eq!(stoul("14"), Ok(14));
        assert_eq!(stoul(" +3x"), Ok(3));
        assert_eq!(stoul("-1"), Ok(u64::MAX));
        assert!(stoul("x").is_err());
        assert!(stoul("99999999999999999999").is_err());
    }
}
