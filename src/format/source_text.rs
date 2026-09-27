// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Source `StringUtils` conversions shared by the search-engine file adapters.
//!
//! The search-engine readers (`KroenikFile`, `MsInspectFile`, `SpecArrayFile`,
//! `OMSSACSVFile`, `PepNovoOutfile`, `InspectOutfile`, `SequestOutfile` and the
//! XML readers) convert every field with `StringUtils::toDouble`,
//! `StringUtils::toInt32` and `StringUtils::split`
//! (`DATASTRUCTURES/StringUtils.h`/`.cpp` at the core pin). Their exact
//! acceptance rules are observable: a field that the source refuses raises
//! `Exception::ConversionError`, which each reader turns into its own
//! `ParseError`, so a port that accepted more or less would load different
//! files. The rules are the same as `cli::context::to_double`/`to_int32`, which
//! the format layer cannot import without closing a module cycle.

/// The source's four whitespace characters, skipped from `index` on
/// (`StringUtils.cpp` `skipWS`).
fn skip_whitespace(bytes: &[u8], mut index: usize) -> usize {
    while bytes
        .get(index)
        .is_some_and(|b| matches!(*b, b' ' | b'\t' | b'\n' | b'\r'))
    {
        index += 1;
    }
    index
}

/// Length of a case-insensitive ASCII `word` at `bytes[index..]`, if present.
fn word_at(bytes: &[u8], index: usize, word: &str) -> Option<usize> {
    let candidate = bytes.get(index..index.checked_add(word.len())?)?;
    candidate
        .eq_ignore_ascii_case(word.as_bytes())
        .then_some(word.len())
}

/// Source `StringUtilsHelper::toInt32`, returning the source's
/// `ConversionError` message as the error.
///
/// Leading and trailing space, tab, newline and carriage return are skipped,
/// one `+` may precede the number, and the rest must be a complete decimal
/// `i32`. Like `std::from_chars` after the source strips `+`, a `-` may still
/// follow it. Overflow is a conversion error.
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

/// Source `StringUtilsHelper::toDouble`, returning the source's
/// `ConversionError` message as the error.
///
/// Whitespace is skipped on both sides and `nan`, optionally followed by a
/// parenthesised payload, is accepted before anything else. Then one `+` may
/// precede a `std::from_chars` general-format number: an optional `-`, digits
/// with an optional fraction and exponent, or `inf`/`infinity`. A finite
/// literal that overflows is an error, as `result_out_of_range` is; underflow
/// rounds. Hexadecimal floats are rejected, following `std::from_chars`.
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

/// Source `StringUtils::split(s, splitter, substrings)` without quote
/// protection: an empty string yields no fields at all, anything else every
/// field between separators, empty fields included.
pub(crate) fn split(text: &str, separator: char) -> Vec<&str> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(separator).collect()
    }
}

/// Source `String::trim`: strips space, tab, newline and carriage return from
/// both ends.
pub(crate) fn trim(text: &str) -> &str {
    text.trim_matches([' ', '\t', '\n', '\r'])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conversions_follow_the_source_rules() {
        assert_eq!(to_double(" +1.5 "), Ok(1.5));
        assert!(to_double("1.5x").is_err());
        assert!(to_double("++1").is_err());
        assert!(to_double("1e999").is_err());
        assert!(to_double("").is_err());
        assert!(to_double("nan(x)").unwrap().is_nan());
        assert_eq!(to_double("-inf"), Ok(f64::NEG_INFINITY));
        assert_eq!(to_int32(" +7\t"), Ok(7));
        assert_eq!(to_int32("-7"), Ok(-7));
        assert!(to_int32("7.0").is_err());
        assert!(to_int32("2147483648").is_err());
        assert!(split("", '\t').is_empty());
        assert_eq!(split("a\t\tb", '\t'), ["a", "", "b"]);
        assert_eq!(trim(" \ta b\r\n"), "a b");
    }
}
