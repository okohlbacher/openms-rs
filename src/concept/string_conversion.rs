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
//! `data_structures::string_utils`'s public conversions (which add the crate
//! error type) and behind the low layers (`metadata`, `concept`)
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

/// Source `StringUtils::split(s, splitter, substrings)` without quote
/// protection (`StringUtils.h:596-663`): an empty string yields no fields at
/// all, anything else every field between separators, empty fields included.
pub(crate) fn split(text: &str, separator: char) -> Vec<&str> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(separator).collect()
    }
}

/// Source `StringUtils::trim` (`StringUtils.h:395-408`): strips space, tab,
/// newline and carriage return from both ends.
pub(crate) fn trim(text: &str) -> &str {
    text.trim_matches([' ', '\t', '\n', '\r'])
}

// The C library conversions: std::stod, std::stoll and std::stoul as glibc
// performs them (moved here from format::arrow_io_helpers, lane F1a, whose
// Release oracle measured them, and format::string_utils, lane F4).

/// C locale `isspace`.
fn is_c_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// `std::stoll(text)`: leading whitespace, an optional sign and the longest
/// run of decimal digits, the rest ignored. `None` where `stoll` throws (no
/// digits, or out of range).
#[cfg(feature = "arrow")]
fn c_stoll(text: &str) -> Option<i64> {
    let bytes = text.as_bytes();
    let mut i = bytes.iter().take_while(|b| is_c_space(**b)).count();
    let negative = match bytes.get(i) {
        Some(b'-') => {
            i += 1;
            true
        }
        Some(b'+') => {
            i += 1;
            false
        }
        _ => false,
    };
    let digits = bytes[i..].iter().take_while(|b| b.is_ascii_digit()).count();
    if digits == 0 {
        return None;
    }
    // Accumulate negatively so that i64::MIN is reachable.
    let mut value: i64 = 0;
    for &digit in &bytes[i..i + digits] {
        value = value
            .checked_mul(10)?
            .checked_sub(i64::from(digit - b'0'))?;
    }
    if negative {
        Some(value)
    } else {
        value.checked_neg()
    }
}

/// Whether `text` starts with `prefix`, ASCII case-insensitively.
fn starts_with_ignore_case(text: &[u8], prefix: &str) -> bool {
    text.len() >= prefix.len() && text[..prefix.len()].eq_ignore_ascii_case(prefix.as_bytes())
}

/// Round `mantissa * 2^exponent` (plus a nonzero tail below it when
/// `sticky`) to a double, to nearest with ties to even. `None` where glibc's
/// `strtod` reports `ERANGE`: overflow, or an inexact result below the
/// smallest normal double (tininess after rounding, as on x86-64).
fn assemble_binary(mantissa: u64, exponent: i64, sticky: bool) -> Option<f64> {
    if mantissa == 0 {
        return Some(0.0);
    }
    let bits = 64 - i64::from(mantissa.leading_zeros());
    let top = exponent.saturating_add(bits - 1);
    if top > 1023 {
        return None;
    }
    let lsb = top.saturating_sub(52).max(-1074);
    let shift = lsb.saturating_sub(exponent);
    let wide = u128::from(mantissa);
    let (mut q, inexact) = if shift <= 0 {
        // At most 53 significant bits remain, so the shifted value fits.
        (wide << u32::try_from(-shift).unwrap_or(0).min(64), sticky)
    } else if shift >= 128 {
        (0, true)
    } else {
        let shift = u32::try_from(shift).unwrap_or(127);
        let q = wide >> shift;
        let rest = wide & ((1u128 << shift) - 1);
        let half = 1u128 << (shift - 1);
        let round_up = rest > half || (rest == half && (sticky || q & 1 == 1));
        (q + u128::from(round_up), rest != 0 || sticky)
    };
    let mut e = lsb;
    if q == 1u128 << 53 {
        q >>= 1;
        e += 1;
    }
    let q = u64::try_from(q).ok()?;
    let tiny = q < 1u64 << 52;
    if inexact && tiny {
        return None;
    }
    if tiny {
        return Some(f64::from_bits(q));
    }
    let biased = e + 52 + 1023;
    if biased >= 2047 {
        return None;
    }
    let biased = u64::try_from(biased).ok()?;
    Some(f64::from_bits((biased << 52) | (q & ((1u64 << 52) - 1))))
}

/// The hexadecimal float after `0x`: `None` when no hex digit follows (the
/// caller then reads the leading `0`), `Some(None)` for `ERANGE`.
fn c_hex_float(bytes: &[u8]) -> Option<Option<f64>> {
    let mut mantissa: u64 = 0;
    let mut exponent: i64 = 0;
    let mut sticky = false;
    let mut any = false;
    let mut dot = false;
    let mut i = 0;
    loop {
        match bytes.get(i) {
            Some(&c) if c.is_ascii_hexdigit() => {
                any = true;
                let digit = u64::from(char::from(c).to_digit(16).unwrap_or(0));
                if mantissa >> 60 == 0 {
                    mantissa = (mantissa << 4) | digit;
                    if dot {
                        exponent -= 4;
                    }
                } else {
                    sticky |= digit != 0;
                    if !dot {
                        exponent = exponent.saturating_add(4);
                    }
                }
            }
            Some(b'.') if !dot => dot = true,
            _ => break,
        }
        i += 1;
    }
    if !any {
        return None;
    }
    if matches!(bytes.get(i), Some(b'p' | b'P')) {
        let mut j = i + 1;
        let negative = match bytes.get(j) {
            Some(b'-') => {
                j += 1;
                true
            }
            Some(b'+') => {
                j += 1;
                false
            }
            _ => false,
        };
        let mut value: i64 = 0;
        let mut digits = 0;
        while let Some(digit) = bytes.get(j).filter(|b| b.is_ascii_digit()) {
            value = value
                .saturating_mul(10)
                .saturating_add(i64::from(digit - b'0'))
                .min(1 << 40);
            digits += 1;
            j += 1;
        }
        if digits > 0 {
            exponent = exponent.saturating_add(if negative { -value } else { value });
        }
    }
    Some(assemble_binary(mantissa, exponent, sticky))
}

/// Base-10^9 digits of `k * 5^1074`, most significant first, as a decimal
/// string: the exact digits of the subnormal double `k * 2^-1074` scaled by
/// 10^1074.
fn subnormal_digits(k: u64) -> String {
    const BASE: u64 = 1_000_000_000;
    let mut limbs: Vec<u64> = Vec::new();
    let mut rest = k;
    while rest > 0 {
        limbs.push(rest % BASE);
        rest /= BASE;
    }
    let mut multiply = |factor: u64| {
        let mut carry = 0u64;
        for limb in &mut limbs {
            let product = *limb * factor + carry;
            *limb = product % BASE;
            carry = product / BASE;
        }
        while carry > 0 {
            limbs.push(carry % BASE);
            carry /= BASE;
        }
    };
    // 1074 = 13 * 82 + 8; 5^13 = 1,220,703,125 keeps each product below 2^64.
    for _ in 0..82 {
        multiply(1_220_703_125);
    }
    multiply(390_625);
    let mut text = String::new();
    for (index, limb) in limbs.iter().rev().enumerate() {
        if index == 0 {
            text.push_str(&limb.to_string());
        } else {
            text.push_str(&format!("{limb:09}"));
        }
    }
    text
}

/// Strip leading and trailing zeros of a digit string, adjusting the decimal
/// exponent of its last digit.
fn normalized(digits: &str, exponent: i64) -> (&str, i64) {
    let digits = digits.trim_start_matches('0');
    let trimmed = digits.trim_end_matches('0');
    let dropped = i64::try_from(digits.len() - trimmed.len()).unwrap_or(i64::MAX);
    (trimmed, exponent.saturating_add(dropped))
}

/// `std::stod(text)` as glibc's `strtod` decides it: leading whitespace, an
/// optional sign, then the longest valid prefix — `inf`/`infinity`, `nan`
/// (with or without a parenthesised payload), a hexadecimal float after `0x`
/// or a decimal float — the rest ignored. `None` where `stod` throws: no
/// conversion, overflow, or an inexact result below the smallest normal
/// double (measured: `1e-310` and `0x1.8p-1074` throw, `0x1p-1074` and
/// `0e-400` do not).
fn c_stod(text: &str) -> Option<f64> {
    let bytes = text.as_bytes();
    let mut i = bytes.iter().take_while(|b| is_c_space(**b)).count();
    let negative = match bytes.get(i) {
        Some(b'-') => {
            i += 1;
            true
        }
        Some(b'+') => {
            i += 1;
            false
        }
        _ => false,
    };
    let sign = if negative { -1.0 } else { 1.0 };
    let rest = &bytes[i..];
    if starts_with_ignore_case(rest, "inf") {
        return Some(sign * f64::INFINITY);
    }
    if starts_with_ignore_case(rest, "nan") {
        return Some(if negative { -f64::NAN } else { f64::NAN });
    }
    if rest.len() >= 2 && rest[0] == b'0' && (rest[1] | 0x20) == b'x' {
        if let Some(value) = c_hex_float(&rest[2..]) {
            return value.map(|value| sign * value);
        }
    }
    // Decimal: digits, an optional point and digits, at least one digit, then
    // an exponent only when it has digits.
    let integer = rest.iter().take_while(|b| b.is_ascii_digit()).count();
    let mut end = integer;
    let mut fraction = 0;
    if rest.get(end) == Some(&b'.') {
        fraction = rest[end + 1..]
            .iter()
            .take_while(|b| b.is_ascii_digit())
            .count();
        if integer + fraction > 0 {
            end += 1 + fraction;
        }
    }
    if integer + fraction == 0 {
        return None;
    }
    let mantissa_end = end;
    let mut exponent: i64 = 0;
    if matches!(rest.get(end), Some(b'e' | b'E')) {
        let mut j = end + 1;
        let exponent_negative = match rest.get(j) {
            Some(b'-') => {
                j += 1;
                true
            }
            Some(b'+') => {
                j += 1;
                false
            }
            _ => false,
        };
        let digits = rest[j..].iter().take_while(|b| b.is_ascii_digit()).count();
        if digits > 0 {
            for &digit in &rest[j..j + digits] {
                exponent = exponent
                    .saturating_mul(10)
                    .saturating_add(i64::from(digit - b'0'))
                    .min(1 << 40);
            }
            if exponent_negative {
                exponent = -exponent;
            }
            end = j + digits;
        }
    }
    // Every byte in rest[..end] is ASCII.
    let literal = std::str::from_utf8(&rest[..end]).ok()?;
    let magnitude: f64 = literal.parse().ok()?;
    if magnitude.is_infinite() {
        return None;
    }
    let digits: String = rest[..mantissa_end]
        .iter()
        .filter(|b| b.is_ascii_digit())
        .map(|b| char::from(*b))
        .collect();
    if magnitude == 0.0 {
        return if digits.bytes().any(|b| b != b'0') {
            None
        } else {
            Some(sign * 0.0)
        };
    }
    if magnitude < f64::MIN_POSITIVE {
        let fraction = i64::try_from(fraction).unwrap_or(i64::MAX);
        let input = normalized(&digits, exponent.saturating_sub(fraction));
        let exact_digits = subnormal_digits(magnitude.to_bits());
        let exact = normalized(&exact_digits, -1074);
        if input != exact {
            return None;
        }
    }
    Some(sign * magnitude)
}

/// Whether C's `strtod`/`strtoll` would read any number at the start of
/// `text` (after leading C white space and one sign): a digit, a point
/// followed by a digit, `inf`, `nan` or `0x`. When it would, a `None` from the
/// conversion is `ERANGE`, not "no conversion".
fn has_numeric_prefix(text: &str, float: bool) -> bool {
    let bytes = text.as_bytes();
    let mut i = bytes.iter().take_while(|b| is_c_space(**b)).count();
    if matches!(bytes.get(i), Some(b'+' | b'-')) {
        i += 1;
    }
    let rest = &bytes[i..];
    match rest.first() {
        Some(b) if b.is_ascii_digit() => true,
        Some(b'.') if float => rest.get(1).is_some_and(u8::is_ascii_digit),
        _ => {
            float && (starts_with_ignore_case(rest, "inf") || starts_with_ignore_case(rest, "nan"))
        }
    }
}

/// `std::stod`: glibc's `strtod` on the text, which skips leading C white
/// space, reads the longest numeric prefix and ignores whatever follows it.
///
/// The prefix is an optional sign and then a decimal number (digits with an
/// optional fraction, at least one digit, and an exponent only when it has
/// digits), a hexadecimal float after `0x`, `inf`, `infinity` or `nan` (any
/// case), or `nan(...)`. The rules are the Linux Release build's, measured by
/// lane F1a through `ArrowIOHelpers::readMetaValues`
/// (`tests/data/arrow_parquet/cpp_release_records.tsv`, `d1`-`d32`): `0x1p3` is
/// 8, `1e-310` and `0x1.8p-1074` are out of range (an inexact result below the
/// smallest normal double), `0x1p-1074` and `0e-400` are not.
///
/// # Errors
///
/// `std::stod`'s `std::invalid_argument` ("no conversion") when no prefix
/// converts, and `std::out_of_range` when the value overflows or is an inexact
/// result below the smallest normal double (`ERANGE` in the C library).
pub(crate) fn stod(text: &str) -> std::result::Result<f64, String> {
    c_stod(text).ok_or_else(|| {
        if has_numeric_prefix(text, true) {
            format!("stod: '{text}' is out of range")
        } else {
            format!("stod: no conversion of '{text}'")
        }
    })
}

/// `std::stoll`: `strtoll` in base 10, which skips leading C white space,
/// takes one sign and the longest run of decimal digits, ignoring the rest.
///
/// # Errors
///
/// `std::invalid_argument` ("no conversion") when no digit follows, and
/// `std::out_of_range` when the digits exceed `i64`.
#[cfg(feature = "arrow")]
pub(crate) fn stoll(text: &str) -> std::result::Result<i64, String> {
    c_stoll(text).ok_or_else(|| {
        if has_numeric_prefix(text, false) {
            format!("stoll: '{text}' is out of range")
        } else {
            format!("stoll: no conversion of '{text}'")
        }
    })
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
    while bytes.get(cursor).copied().is_some_and(is_c_space) {
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
    /// invocations, identical; tests/data/string_conversion_release_provenance.json).
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
        // The Release build reads hexadecimal floats (readMetaValues d12,
        // tests/data/arrow_parquet/cpp_release_records.tsv: 0x1p3 -> 8.0);
        // lane F4 had refused them.
        assert_eq!(stod("0x1p3"), Ok(8.0));
    }

    #[test]
    fn stoul_follows_strtoul() {
        assert_eq!(stoul("14"), Ok(14));
        assert_eq!(stoul(" +3x"), Ok(3));
        assert_eq!(stoul("-1"), Ok(u64::MAX));
        assert!(stoul("x").is_err());
        assert!(stoul("99999999999999999999").is_err());
    }

    #[test]
    #[cfg(feature = "arrow")]
    fn stoll_follows_the_c_library() {
        assert_eq!(stoll("12"), Ok(12));
        assert_eq!(stoll("\x0b\x0c 7"), Ok(7));
        assert!(stoll("- 1").is_err());
        assert_eq!(stoll("-9223372036854775808"), Ok(i64::MIN));
        assert!(stoll("9223372036854775808").is_err());
    }

    #[test]
    fn stod_reports_range_errors_as_glibc() {
        assert_eq!(stod("0x1p-1074"), Ok(f64::from_bits(1)));
        assert!(stod("0x1.8p-1074").is_err());
        assert!(stod("1e-310").is_err());
        assert_eq!(stod("0e-400"), Ok(0.0));
        assert_eq!(stod("0x"), Ok(0.0));
        assert_eq!(stod("0x.8"), Ok(0.5));
        assert!(stod("0x1p1024").is_err());
        assert_eq!(stod("0x1.fffffffffffffp1023"), Ok(f64::MAX));
        // (lane F1a) The exact decimal expansion of the smallest subnormal is exact.
        let exact = format!("{}e-1074", subnormal_digits(1));
        assert_eq!(stod(&exact), Ok(f64::from_bits(1)));
    }
}
