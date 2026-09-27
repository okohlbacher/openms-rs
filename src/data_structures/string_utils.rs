// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Source `StringUtils` conversions and splitting, as the file readers and the
//! command line use them.
//!
//! Ports the parts of `DATASTRUCTURES/StringUtils.h`/`.cpp` whose exact
//! acceptance rules are observable: `StringUtilsHelper::toInt32`
//! (`StringUtils.cpp:136-166`), `StringUtilsHelper::toDouble`
//! (`StringUtils.cpp:239-276`), the single-character `StringUtils::split`
//! (`StringUtils.h:596-663`) and `StringUtils::trim` (`StringUtils.h:395-408`),
//! plus the C library conversions `std::stod` and `std::stoul` that some text
//! readers (`MSPGenericFile`) call directly instead of `StringUtils`.
//! A field the source refuses raises `Exception::ConversionError`, which the
//! text readers (`KroenikFile`, `MsInspectFile`, `SpecArrayFile`,
//! `OMSSACSVFile`, the Sequest, Inspect and PepNovo adapters) and the command
//! line (`cli::context`) report in their own way, so a port that accepted more
//! or less than the source would load different files. The rest of the header
//! is not ported here.
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

/// Source `StringUtils::toInt32` (`StringUtils.cpp:136-166`), returning the
/// source's `ConversionError` message as the error.
///
/// Leading and trailing space, tab, newline and carriage return are skipped,
/// one `+` may precede the number, and the rest must be a complete decimal
/// `i32`. Like `std::from_chars` after the source strips `+`, a `-` may still
/// follow it. Overflow is a conversion error.
///
/// # Errors
///
/// The source's `ConversionError` message when the text is not a complete
/// `i32`.
pub fn to_int32(text: &str) -> std::result::Result<i32, String> {
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

/// Source `StringUtils::toDouble` (`StringUtils.cpp:239-276`), returning the
/// source's `ConversionError` message as the error.
///
/// Whitespace is skipped on both sides and `nan`, optionally followed by a
/// parenthesised payload, is accepted before anything else. Then one `+` may
/// precede a `std::from_chars` general-format number: an optional `-`, digits
/// with an optional fraction and exponent, `inf`/`infinity`, or `nan`
/// optionally followed by `(n-char-sequence)` of letters, digits and `_`
/// (`-nan` keeps its sign bit). These are the rules of libstdc++'s
/// `std::from_chars`, which the Linux Release build uses: a finite literal that
/// overflows, and a nonzero literal that rounds to zero, are both
/// `result_out_of_range` and so a conversion error, while a literal that
/// rounds to a subnormal is accepted. Hexadecimal floats are rejected. (The
/// source's `strtod` fallback, compiled only on libc++, accepts underflow; it
/// is not the behaviour of the Linux build. Release records:
/// `tests/data/string_utils_release/cpp_release_records.tsv`.)
///
/// # Errors
///
/// The source's `ConversionError` message when the text is not a complete
/// number.
pub fn to_double(text: &str) -> std::result::Result<f64, String> {
    to_floating::<f64>(text, "double")
}

/// Source `StringUtils::toFloat` (`StringUtils.cpp:200-237`): the rules of
/// [`to_double`], with the literal rounded once, directly to `f32`, as
/// `std::from_chars` into a `float` does.
///
/// # Errors
///
/// The source's `ConversionError` message when the text is not a complete
/// number or a finite literal overflows `f32`.
pub fn to_float(text: &str) -> std::result::Result<f32, String> {
    to_floating::<f32>(text, "float")
}

/// The two floating-point targets of the source's conversions.
trait SourceFloat: std::str::FromStr + Copy {
    const NAN: Self;
    fn infinite(self) -> bool;
    fn zero(self) -> bool;
    fn negated(self) -> Self;
}

impl SourceFloat for f64 {
    const NAN: Self = f64::NAN;
    fn infinite(self) -> bool {
        self.is_infinite()
    }
    fn zero(self) -> bool {
        self == 0.0
    }
    fn negated(self) -> Self {
        -self
    }
}

impl SourceFloat for f32 {
    const NAN: Self = f32::NAN;
    fn infinite(self) -> bool {
        self.is_infinite()
    }
    fn zero(self) -> bool {
        self == 0.0
    }
    fn negated(self) -> Self {
        -self
    }
}

fn to_floating<T: SourceFloat>(text: &str, kind: &str) -> std::result::Result<T, String> {
    let bytes = text.as_bytes();
    let not_converted = || format!("Could not convert string '{text}' to a {kind} value");
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
            return Ok(T::NAN);
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
    let mut nonzero_digit = false;
    if let Some(length) = infinity {
        cursor += length;
    } else if let Some(length) = nan {
        cursor += length;
        // `nan(n-char-sequence)`: consumed only when the parenthesis closes
        // over letters, digits and underscores; otherwise `nan` alone is read.
        if bytes.get(cursor) == Some(&b'(') {
            let payload = cursor + 1;
            let mut end = payload;
            while bytes
                .get(end)
                .is_some_and(|b| b.is_ascii_alphanumeric() || *b == b'_')
            {
                end += 1;
            }
            if bytes.get(end) == Some(&b')') {
                cursor = end + 1;
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
        nonzero_digit = bytes[integer_start..cursor]
            .iter()
            .any(|b| matches!(b, b'1'..=b'9'));
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
    let value: T = if nan.is_some() && infinity.is_none() {
        if negative { T::NAN.negated() } else { T::NAN }
    } else {
        literal.parse().map_err(|_| not_converted())?
    };
    // libstdc++ `from_chars`: out of range on overflow to infinity and on a
    // nonzero literal that rounds to zero.
    if (value.infinite() && infinity.is_none()) || (nonzero_digit && value.zero()) {
        return Err(not_converted());
    }
    let after = skip_whitespace(bytes, cursor);
    if after != bytes.len() {
        return Err(format!(
            "Prefix of string '{text}' successfully converted to a {kind} value. Additional characters found at position {}",
            after + 1
        ));
    }
    Ok(value)
}

/// Source `StringUtils::split(s, splitter, substrings)` without quote
/// protection (`StringUtils.h:596-663`): an empty string yields no fields at
/// all, anything else every field between separators, empty fields included.
pub fn split(text: &str, separator: char) -> Vec<&str> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(separator).collect()
    }
}

/// Source `StringUtils::trim` (`StringUtils.h:395-408`): strips space, tab,
/// newline and carriage return from both ends.
pub fn trim(text: &str) -> &str {
    text.trim_matches([' ', '\t', '\n', '\r'])
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
pub fn stod(text: &str) -> std::result::Result<f64, String> {
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
pub fn stoul(text: &str) -> std::result::Result<u64, String> {
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
    fn conversions_follow_the_source_rules() {
        assert_eq!(to_double(" +1.5 "), Ok(1.5));
        assert!(to_double("1.5x").is_err());
        assert!(to_double("++1").is_err());
        assert!(to_double("1e999").is_err());
        assert!(to_double("").is_err());
        assert!(to_double("nan(x)").unwrap().is_nan());
        assert_eq!(to_double("-inf"), Ok(f64::NEG_INFINITY));
        assert_eq!(to_float("866.606"), Ok(866.606_f32));
        assert!(to_float("1e39").is_err());
        assert!(to_float("1.5x").unwrap_err().contains("to a float value"));
        assert_eq!(to_int32(" +7\t"), Ok(7));
        assert_eq!(to_int32("-7"), Ok(-7));
        assert!(to_int32("7.0").is_err());
        assert!(to_int32("2147483648").is_err());
        assert!(split("", '\t').is_empty());
        assert_eq!(split("a\t\tb", '\t'), ["a", "", "b"]);
        assert_eq!(trim(" \ta b\r\n"), "a b");
    }

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
