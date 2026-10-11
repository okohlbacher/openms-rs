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
//! plus the C library conversions `std::stod`, `std::stoll` and `std::stoul`
//! that some readers (`MSPGenericFile`, `ArrowIOHelpers::readMetaValues`) call
//! directly instead of `StringUtils`, with glibc's rules.
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
pub fn stod(text: &str) -> std::result::Result<f64, String> {
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
pub fn stoll(text: &str) -> std::result::Result<i64, String> {
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
pub fn stoul(text: &str) -> std::result::Result<u64, String> {
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
