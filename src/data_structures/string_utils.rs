// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `DATASTRUCTURES/StringUtils.h`: the string utility namespace of the source.
//!
//! The source header (`src/openms/include/OpenMS/DATASTRUCTURES/StringUtils.h`
//! and `src/openms/source/DATASTRUCTURES/StringUtils.cpp` at core `bc9cc12`)
//! replaced the former `OpenMS::String` class with free functions over
//! `std::string`: number parsing (`StringUtilsHelper`), number formatting
//! (`toStr`, `appendToStr`, `number`, `numberLength`), predicates, accessors,
//! in-place mutators, quoting, splitting and joining, the legacy
//! `StringConversions` forwarding namespace and global `std::string + number`
//! operators. This module ports every public member; the table below maps each
//! one. The support document is `docs/STRING_UTILS_SUPPORT.md`.
//!
//! # Text model
//!
//! A C++ `std::string` is a byte sequence and every function here is defined on
//! bytes; a Rust [`String`] is UTF-8. Every function keeps the source's byte
//! semantics on ASCII text, which is all the source's own callers feed it. The
//! differences are confined to text a C++ `char` cannot address: a `char`
//! argument that is not ASCII is refused where the source's `char` is one byte
//! of the text (quotes, fill characters, splitters), and a byte count that
//! would cut a multi-byte character is refused ([`prefix`], [`suffix`],
//! [`substr`], [`chop`]), because the source's result there is not UTF-8.
//! [`reverse`] refuses non-ASCII text for the same reason.
//!
//! # API mapping
//!
//! | Source member | Rust | Notes |
//! |---|---|---|
//! | `enum class QuotingMethod { NONE, ESCAPE, DOUBLE }` | [`QuotingMethod`] | `ESCAPE`, the source default argument, is [`QuotingMethod::default`] |
//! | `StringUtilsHelper::toInt32`, `StringUtils::toInt32` | [`to_int32`] | |
//! | `StringUtilsHelper::toInt64`, `StringUtils::toInt64` | [`to_int64`] | |
//! | `StringUtilsHelper::toFloat`, `StringUtils::toFloat` | [`to_float`] | |
//! | `StringUtilsHelper::toDouble`, `StringUtils::toDouble` | [`to_double`] | |
//! | `extractDouble` (`const char*` and both `std::string` iterator overloads) | [`extract_double`] | returns the value and the bytes consumed |
//! | `extractInt` (the same three overloads) | [`extract_int`] | |
//! | `StringUtilsHelper()` / `~StringUtilsHelper()` | not needed | a class of static functions; Rust has free functions |
//! | `appendToStr(int / unsigned / short / unsigned short / long / unsigned long / long long / unsigned long long, std::string&)` | [`ToStr::append_to_str`] | |
//! | `toStr` of the same eight integer types | [`ToStr::to_str_text`] | decimal |
//! | `appendToStr(float, std::string&)`, `toStr(float, true)` | [`ToStr`] for `f32`, [`to_str_f32`] | |
//! | `appendToStr(double, std::string&)`, `toStr(double, true)` | [`ToStr`] for `f64`, [`to_str`] | |
//! | `appendToStrLowP(float, std::string&)`, `toStr(float, false)` | [`append_to_str_low_p_f32`], [`to_str_f32_low_p`] | |
//! | `appendToStrLowP(double, std::string&)`, `toStr(double, false)` | [`append_to_str_low_p`], [`to_str_low_p`] | |
//! | `appendToStr(long double, ...)`, `appendToStrLowP(long double, ...)`, `toStr(long double, bool)` | not portable | Rust has no extended-precision type; the `f64` functions cover every value a Rust caller can hold |
//! | `toStr(char)` | [`ToStr`] for `char` | one character |
//! | `toStr(const std::string&)`, `toStr(const char*)`, `toStr(std::string_view)` | [`ToStr`] for `str` and `String` | copied unchanged |
//! | `appendToStr(const DataValue&, bool, std::string&)`, `toStr(const DataValue&, bool)` | [`MetaValue::to_text`](crate::metadata::MetaValue::to_text) | lives beside the value type, which `metadata` owns |
//! | `toStr(const ParamValue&, bool)` | [`ParamValue::to_text`](crate::param::ParamValue::to_text) | |
//! | `skipWhitespace(const char*, const char*)`, `skipWhitespace(const std::string_view&)` | [`skip_whitespace`] | a byte count |
//! | `skipNonWhitespace` (both overloads) | [`skip_non_whitespace`] | |
//! | `number(double, UInt)` | [`number`], [`number_truncated`] | strict and source-compatible |
//! | `numberLength(double, UInt)` | [`number_length`] | |
//! | `random(UInt)` | [`random`], [`random_with`] | |
//! | `hasPrefix` (`std::string`, `char` and `string_view` overloads) | [`str::starts_with`] | the byte test is the same test on UTF-8 |
//! | `hasSuffix` (three overloads) | [`str::ends_with`] | |
//! | `hasSubstring` (three overloads) | [`str::contains`] | |
//! | `has(const std::string&, Byte)`, `has(const std::string&, char)` | `text.as_bytes().contains(&byte)`, [`str::contains`] | |
//! | `prefix(s, size_t)`, `suffix(s, size_t)` | [`prefix`], [`suffix`] | |
//! | `prefix(s, Int)`, `suffix(s, Int)` | [`prefix_signed`], [`suffix_signed`] | |
//! | `prefix(s, char)`, `suffix(s, char)` | [`prefix_until`], [`suffix_after`] | |
//! | `substr` (`std::string` and `string_view` overloads) | [`substr`] | `npos` is `None` |
//! | `chop(s, Size)` | [`chop`] | |
//! | `trim(std::string&)`, `trimmed(std::string)` | [`trim`], [`trimmed`] | |
//! | `toUpper`, `toLower`, `toUppered`, `toLowered` | [`to_upper`], [`to_lower`], [`to_uppered`], [`to_lowered`] | the classic `"C"` locale: ASCII letters only |
//! | `firstToUpper` | [`first_to_upper`] | |
//! | `reverse` | [`reverse`] | |
//! | `simplify` | [`simplify`] | |
//! | `fillLeft`, `fillRight` (lvalue and rvalue overloads) | [`fill_left`], [`fill_right`] | an owned `String` is filled in place |
//! | `substitute(s, char, char)`, `substituted(s, char, char)` | [`substitute_char`], [`substituted_char`] | |
//! | `substitute(s, string, string)`, `substituted(s, string, string)` | [`substitute`], [`substituted`] | |
//! | `remove(s, char)` | [`remove`] | |
//! | `ensureLastChar` (lvalue and rvalue overloads) | [`ensure_last_char`] | |
//! | `removeWhitespaces` | [`remove_whitespaces`] | |
//! | `isQuoted`, `quote`, `unquote` | [`is_quoted`], [`quote`], [`unquote`] | |
//! | `split(s, char, vector&, bool quote_protect)` | [`split`] | returns [`Split`] |
//! | `split(s, const std::string&, vector&)` | [`split_str`] | |
//! | `split_quoted` | [`split_quoted`] | |
//! | `concatenate` (target, iterator and container overloads) | [`concatenate`] | |
//! | `StringConversions::append`, `appendLowP`, `toString`, `toStringLowP` | the functions above | a forwarding namespace |
//! | global `std::string& operator+=(std::string&, number)` | [`ToStr::append_to_str`] | |
//! | global `std::string operator+(std::string, number)` | [`ToStr::append_to_str`] on an owned `String` | |
//!
//! # Checked boundaries
//!
//! Functions whose result can outgrow their input ([`fill_left`],
//! [`fill_right`], [`substitute`], [`quote`], [`concatenate`], [`random`])
//! refuse a result above [`MAX_TEXT_BYTES`] before allocating and reserve
//! fallibly; the source has no ceiling. Splitting refuses more than
//! [`MAX_SPLIT_PARTS`] parts. Every refusal leaves the input unchanged.
//!
//! # Evidence
//!
//! Tier 1: `../oracle/B8a-stringutils` links the Linux x86_64 Release install
//! on ibminode06 and records 5,068 calls over the whole surface (two identical
//! runs); `tests/string_utils.rs` embeds the records and compares every one.
//! The 246 `toInt32`/`toDouble`/`toFloat` rows of the earlier
//! `../oracle/integ-F-stringutils` (same install, same host) use 82 inputs that
//! are all among these records, so they are not compared a second time. The upstream
//! `StringUtils_test.cpp` and `String_test.cpp` literals are ported in the same
//! test (tier 3).

use crate::concept::numeric_formatting::{
    append_numeric_f32, append_numeric_f64, ostream_g, printf_nonfinite, printf_precision,
};
use crate::error::{Error, Result};

/// Largest text, in bytes, that `StringUtils::number` can return.
///
/// The source formats into `char buf[64]` with `snprintf`, which stores at most
/// 63 bytes before the terminating NUL and silently drops the rest.
pub const NUMBER_MAX_TEXT_BYTES: usize = 63;

/// `writtenDigits<float>()`, which is `std::numeric_limits<float>::digits10`:
/// the fraction-digit count of `toStr(float)` and the stream precision
/// FileInfo sets for `float` statistics.
pub const WRITTEN_DIGITS_F32: u32 = 6;

/// `writtenDigits<double>()`, which is `std::numeric_limits<double>::digits10`:
/// the fraction-digit count of `toStr(double)`, the stream precision of
/// `numberLength` and of FileInfo's `double` statistics.
pub const WRITTEN_DIGITS_F64: u32 = 15;

/// Fraction digits of the low-precision `toStr(x, false)` and `appendToStrLowP`.
const LOW_PRECISION_DIGITS: u32 = 3;

/// Native ceiling, in bytes, of any text a function here builds that can grow
/// beyond its input. The source has none.
pub const MAX_TEXT_BYTES: usize = 1 << 30;

/// Native ceiling on the number of parts a split may produce. The source has
/// none.
pub const MAX_SPLIT_PARTS: usize = 1 << 26;

/// Longest [`to_str`] text: sign, 17 significant digits, point and `e-308`.
const F64_TO_STR_MAX_BYTES: usize = 24;
/// Longest [`to_str_f32`] text: sign, 9 significant digits, point and `e-45`.
const F32_TO_STR_MAX_BYTES: usize = 15;
/// Longest decimal 64-bit integer, `-9223372036854775808`.
const INTEGER_TO_STR_MAX_BYTES: usize = 20;

/// The four characters every whitespace rule of the header uses: space, tab,
/// line feed and carriage return. Vertical tab and form feed are not
/// whitespace here.
const fn is_blank(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r')
}

fn blank_char(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\r')
}

/// How embedded quote characters are treated when quoting, unquoting and
/// splitting, the source `enum class QuotingMethod`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum QuotingMethod {
    /// `NONE`: quote characters inside the text are not marked.
    None,
    /// `ESCAPE`: a quote character (and a backslash) is escaped with a
    /// backslash. The source's default argument.
    #[default]
    Escape,
    /// `DOUBLE`: a quote character is doubled.
    Double,
}

// ---------------------------------------------------------------------------
// Parsing (StringUtilsHelper)
// ---------------------------------------------------------------------------

/// Index of the first byte at or after `index` that is not one of the four
/// blanks.
fn skip_blanks(bytes: &[u8], mut index: usize) -> usize {
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
fn not_converted(text: &str, kind: &str) -> Error {
    Error::InvalidValue(format!("Could not convert string '{text}' to {kind} value"))
}

/// The source's `ConversionError` for a converted prefix followed by other
/// characters; `after` is the zero-based index of the first of them.
fn trailing(text: &str, kind: &str, after: usize) -> Error {
    Error::InvalidValue(format!(
        "Prefix of string '{text}' successfully converted to {kind} value. \
         Additional characters found at position {}",
        after + 1
    ))
}

/// `StringUtilsHelper::toInt32` and `toInt64` (`StringUtils.cpp:136-200`).
fn parse_integer<T: std::str::FromStr>(text: &str, prefix_kind: &str) -> Result<T> {
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
fn parse_float<T: FromCharsFloat>(text: &str, kind: &str) -> Result<T> {
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

/// `StringUtils::toInt32`: the whole text as a 32-bit integer.
///
/// Space, tab, line feed and carriage return are skipped on both sides, one `+`
/// may precede the number, and the rest is `std::from_chars`: an optional `-`
/// (so `+-1` is `-1`, while `++1` and `-+1` fail) and decimal digits. Anything
/// after the number other than those blanks is an error.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] with the source's `ConversionError` text:
/// "Could not convert string '…' to an integer value" when there is no number
/// or it is outside `i32`, and "Prefix of string '…' successfully converted to
/// an int32 value. Additional characters found at position N" (one-based) when
/// characters follow it.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::to_int32;
///
/// assert_eq!(to_int32("  2147483647 ").ok(), Some(2_147_483_647));
/// assert_eq!(to_int32("+-1").ok(), Some(-1));
/// assert!(to_int32("2147483648").is_err());
/// assert!(to_int32("1234  moreText").is_err());
/// ```
pub fn to_int32(text: &str) -> Result<i32> {
    parse_integer(text, "an int32")
}

/// `StringUtils::toInt64`: [`to_int32`]'s rules for a 64-bit integer.
///
/// # Errors
///
/// As [`to_int32`]; the text after a converted prefix says "an integer value",
/// as the source's does.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::to_int64;
///
/// assert_eq!(to_int64("9223372036854775807").ok(), Some(i64::MAX));
/// assert!(to_int64("-9223372036854775809").is_err());
/// ```
pub fn to_int64(text: &str) -> Result<i64> {
    parse_integer(text, "an integer")
}

/// `StringUtils::toFloat`: the whole text as a `float`.
///
/// Blanks are skipped on both sides. `nan` in any case, optionally followed by
/// `(` anything `)`, is a positive quiet NaN. Otherwise one `+` may precede a
/// `std::from_chars` literal: an optional `-`, then `inf`, `infinity`, `nan`
/// with an optional `(letters, digits, _)`, or a decimal number with optional
/// fraction and exponent (`1.`, `.5`, `1e3`), converted with correct rounding
/// to the `float` itself. `-nan` keeps its sign bit.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] with the source's `ConversionError` text
/// when there is no literal, when a finite literal overflows (`1e40`) or a
/// literal with a nonzero digit rounds to zero (`7e-46`) — the libstdc++
/// `result_out_of_range` of the reference build; a subnormal (`1e-40`) is a
/// value — and when characters other than blanks follow the literal.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::to_float;
///
/// assert_eq!(to_float("  1234.45 ").ok(), Some(1234.45_f32));
/// assert!(to_float("1e-40").is_ok_and(|v| v > 0.0));
/// assert!(to_float("1e40").is_err());
/// ```
pub fn to_float(text: &str) -> Result<f32> {
    parse_float(text, "a float")
}

/// `StringUtils::toDouble`: [`to_float`]'s rules for a `double`.
///
/// # Errors
///
/// As [`to_float`], at `double` range: `1e400` overflows, `1e-400` and
/// `2e-324` round to zero and are refused, `3e-324` is the smallest subnormal.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::to_double;
///
/// assert_eq!(to_double("2.17388884170148e-321").ok(), Some(2.17388884170148e-321));
/// assert!(to_double("-nan").is_ok_and(|v| v.is_nan() && v.is_sign_negative()));
/// assert!(to_double("1e400").is_err());
/// assert!(to_double(" 1234.45 911.0").is_err());
/// ```
pub fn to_double(text: &str) -> Result<f64> {
    parse_float(text, "a double")
}

/// `StringUtils::extractDouble`: the `double` at the start of `text` and the
/// number of bytes it occupies, or `None` when there is none.
///
/// No whitespace is skipped. `nan` with an optional `(` anything `)` is taken
/// first; otherwise one `+` may precede a `std::from_chars` literal (see
/// [`to_float`]). An out-of-range literal is no value. The source advances its
/// iterator by the returned count and leaves both the iterator and the target
/// unchanged on failure.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::extract_double;
///
/// assert_eq!(extract_double("15.0e6 x"), Some((15.0e6, 6)));
/// assert_eq!(extract_double("+1234.45!"), Some((1234.45, 8)));
/// assert_eq!(extract_double("  -123.45"), None);
/// ```
pub fn extract_double(text: &str) -> Option<(f64, usize)> {
    let bytes = text.as_bytes();
    if let Some(end) = try_parse_nan(bytes, 0) {
        return Some((f64::NAN, end));
    }
    let start = usize::from(bytes.first() == Some(&b'+'));
    from_chars_float::<f64>(text, start)
}

/// `StringUtils::extractInt`: the `int` at the start of `text` and the number
/// of bytes it occupies, or `None` when there is none.
///
/// No whitespace is skipped; one `+` may precede `std::from_chars`'s optional
/// `-` and digits. An out-of-range number is no value.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::extract_int;
///
/// assert_eq!(extract_int("+1234!"), Some((1234, 5)));
/// assert_eq!(extract_int("-5\t9"), Some((-5, 2)));
/// assert_eq!(extract_int("  -123"), None);
/// ```
pub fn extract_int(text: &str) -> Option<(i32, usize)> {
    let start = usize::from(text.as_bytes().first() == Some(&b'+'));
    from_chars_int::<i32>(text, start)
}

// ---------------------------------------------------------------------------
// Whitespace scanning
// ---------------------------------------------------------------------------

/// `StringUtils::skipWhitespace`: the number of leading space, tab, line feed
/// and carriage return bytes.
///
/// The source scans sixteen bytes at a time with SIMD; the result is the same.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::skip_whitespace;
///
/// assert_eq!(skip_whitespace(" \t\n\rx"), 4);
/// assert_eq!(skip_whitespace("x "), 0);
/// ```
pub fn skip_whitespace(text: &str) -> usize {
    skip_blanks(text.as_bytes(), 0)
}

/// `StringUtils::skipNonWhitespace`: the number of leading bytes before the
/// first space, tab, line feed or carriage return (the whole length when there
/// is none).
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::skip_non_whitespace;
///
/// assert_eq!(skip_non_whitespace("0 2"), 1);
/// assert_eq!(skip_non_whitespace("abc"), 3);
/// ```
pub fn skip_non_whitespace(text: &str) -> usize {
    text.as_bytes()
        .iter()
        .position(|b| is_blank(*b))
        .unwrap_or(text.len())
}

// ---------------------------------------------------------------------------
// Number to text
// ---------------------------------------------------------------------------

/// `StringUtils::toStr(double)` (`toStr(d, true)`, `appendToStr(double)`): the
/// full-precision text of `NumericFormatting::appendNumeric(d, s, 15, false)`.
///
/// Zero and magnitudes in `[1e-2, 1e4)` use fixed notation with 15 fraction
/// digits, rounded exactly with ties to even. Trailing zeros are removed but one
/// fraction digit stays (`5.0`, `-0.0`). These are 15 *fraction* digits, not the
/// 15 significant digits the source comment describes, so binary representation
/// digits show: `1234.5678` prints as `1234.567800000000034`. Other magnitudes
/// use the shortest scientific text that reads back as the same double, with at
/// least one fraction digit and an exponent without `+` of at least two digits
/// (`1.0e-05`, `3.6739e04`, `-5.0e-324`). As in `std::to_chars`, the text nearest
/// the exact value wins among equally short ones, and an exact tie takes the even
/// digit (`-1.3516970629882812e05`, not `...813e05`). NaN prints as `NaN`, the
/// infinities as `inf` and `-inf`.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::to_str;
///
/// assert_eq!(to_str(-65.0), "-65.0");
/// assert_eq!(to_str(1e-5), "1.0e-05");
/// assert_eq!(to_str(36739.0), "3.6739e04");
/// assert_eq!(to_str(-135_169.706_298_828_125), "-1.3516970629882812e05");
/// ```
pub fn to_str(value: f64) -> String {
    let mut text = String::new();
    append_numeric_f64(value, &mut text, WRITTEN_DIGITS_F64, false);
    text
}

/// `StringUtils::toStr(float)`: [`to_str`]'s rules applied to the `float`
/// itself with `writtenDigits<float>()`.
///
/// Fixed notation keeps 6 fraction digits of the exact `float` value
/// (`9999.999023`), and scientific notation uses the shortest text that reads
/// back as the same `float` (`9.999999e-03`, `3.4028235e38`), not as the same
/// promoted `double`, with the same tie rule (`2203265.25` prints as
/// `2.2032652e06`).
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::to_str_f32;
///
/// assert_eq!(to_str_f32(50.840_353), "50.840351");
/// assert_eq!(to_str_f32(1e-5), "1.0e-05");
/// assert_eq!(to_str_f32(2_203_265.25), "2.2032652e06");
/// ```
pub fn to_str_f32(value: f32) -> String {
    let mut text = String::new();
    append_numeric_f32(value, &mut text, WRITTEN_DIGITS_F32, false);
    text
}

/// `StringUtils::toStr(double, false)`: three fraction digits.
///
/// `NumericFormatting::appendNumeric(d, s, 3, true)`: zero and magnitudes in
/// `[1e-2, 1e4)` with three fraction digits rounded with ties to even
/// (`0.0625` is `0.062`), trailing zeros removed down to one (`1.0`), other
/// magnitudes as `std::to_chars(scientific, 3)` with the same trimming
/// (`1.235e04`, `5.0e-04`, `4.941e-324`). NaN, `inf` and `-inf` as in
/// [`to_str`].
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::to_str_low_p;
///
/// assert_eq!(to_str_low_p(12345.6789), "1.235e04");
/// assert_eq!(to_str_low_p(0.0625), "0.062");
/// assert_eq!(to_str_low_p(0.0005), "5.0e-04");
/// ```
pub fn to_str_low_p(value: f64) -> String {
    let mut text = String::new();
    append_to_str_low_p(value, &mut text);
    text
}

/// `StringUtils::toStr(float, false)`: [`to_str_low_p`]'s rules applied to the
/// `float` itself.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::to_str_f32_low_p;
///
/// assert_eq!(to_str_f32_low_p(12345.678), "1.235e04");
/// ```
pub fn to_str_f32_low_p(value: f32) -> String {
    let mut text = String::new();
    append_to_str_low_p_f32(value, &mut text);
    text
}

/// `StringUtils::appendToStrLowP(double, std::string&)`: append
/// [`to_str_low_p`]'s text to `target`.
pub fn append_to_str_low_p(value: f64, target: &mut String) {
    append_numeric_f64(value, target, LOW_PRECISION_DIGITS, true);
}

/// `StringUtils::appendToStrLowP(float, std::string&)`: append
/// [`to_str_f32_low_p`]'s text to `target`.
pub fn append_to_str_low_p_f32(value: f32, target: &mut String) {
    append_numeric_f32(value, target, LOW_PRECISION_DIGITS, true);
}

/// The `StringUtils::toStr` / `appendToStr` overload set, and the global
/// `std::string` `operator+=` and `operator+` with a number, which append the
/// same text.
///
/// Implemented for `f64` ([`to_str`]), `f32` ([`to_str_f32`]), the signed and
/// unsigned 16-, 32- and 64-bit and pointer-sized integers (decimal), `char`
/// (the character), `str` and `String` (copied unchanged), and references to
/// any of these. `bool` is not implemented: C++ would promote it to `int` and
/// print `1` or `0`.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::ToStr;
///
/// let mut text = String::from("x=");
/// 2.5_f64.append_to_str(&mut text); // std::string("x=") + 2.5
/// assert_eq!(text, "x=2.5");
/// assert_eq!((-7_i64).to_str_text(), "-7");
/// ```
pub trait ToStr {
    /// An upper bound of the bytes [`ToStr::append_to_str`] appends.
    ///
    /// Callers that build long texts sum these bounds to enforce a ceiling
    /// before they allocate, so an implementation must not return less than it
    /// appends.
    fn text_len_bound(&self) -> usize;

    /// Appends the `StringUtils::toStr` text of `self` to `target`, as
    /// `appendToStr(value, target)` and `target += value` do.
    fn append_to_str(&self, target: &mut String);

    /// The `StringUtils::toStr` text of `self`.
    fn to_str_text(&self) -> String {
        let mut text = String::new();
        self.append_to_str(&mut text);
        text
    }
}

impl ToStr for f64 {
    fn text_len_bound(&self) -> usize {
        F64_TO_STR_MAX_BYTES
    }
    fn append_to_str(&self, target: &mut String) {
        append_numeric_f64(*self, target, WRITTEN_DIGITS_F64, false);
    }
}

impl ToStr for f32 {
    fn text_len_bound(&self) -> usize {
        F32_TO_STR_MAX_BYTES
    }
    fn append_to_str(&self, target: &mut String) {
        append_numeric_f32(*self, target, WRITTEN_DIGITS_F32, false);
    }
}

impl ToStr for char {
    fn text_len_bound(&self) -> usize {
        self.len_utf8()
    }
    fn append_to_str(&self, target: &mut String) {
        target.push(*self);
    }
}

impl ToStr for str {
    fn text_len_bound(&self) -> usize {
        self.len()
    }
    fn append_to_str(&self, target: &mut String) {
        target.push_str(self);
    }
}

impl ToStr for String {
    fn text_len_bound(&self) -> usize {
        self.len()
    }
    fn append_to_str(&self, target: &mut String) {
        target.push_str(self);
    }
}

impl<T: ToStr + ?Sized> ToStr for &T {
    fn text_len_bound(&self) -> usize {
        (**self).text_len_bound()
    }
    fn append_to_str(&self, target: &mut String) {
        (**self).append_to_str(target);
    }
}

macro_rules! integer_to_str {
    ($($integer:ty),*) => {$(
        impl ToStr for $integer {
            fn text_len_bound(&self) -> usize {
                INTEGER_TO_STR_MAX_BYTES
            }
            fn append_to_str(&self, target: &mut String) {
                target.push_str(&self.to_string());
            }
        }
    )*};
}
integer_to_str!(i16, u16, i32, u32, i64, u64, isize, usize);

/// `StringUtils::number(d, n)`: `value` with exactly `digits` fraction digits,
/// as C `snprintf("%.*f")` writes it.
///
/// Rounding uses the exact binary value with ties to even, and the sign of zero
/// is kept (`-0.00`). NaN and the infinities ignore `digits` and print as
/// `printf` writes them: `nan`, `-nan` for a NaN whose sign bit is set, `inf`
/// and `-inf` (glibc, the reference build's C library). The source parameter
/// is a `double`, so a `float` argument is promoted first; pass `f64::from(x)`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] where the source would not print the full
/// number: the text is longer than [`NUMBER_MAX_TEXT_BYTES`], which `number`
/// silently cuts, or `digits` exceeds `i32::MAX` for a finite value, which the
/// source converts to a negative `int` and so prints with six digits. Use
/// [`number_truncated`] where byte parity with those outputs is required.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::number;
///
/// assert_eq!(number(-1.0, 2).ok().as_deref(), Some("-1.00"));
/// assert_eq!(number(0.125, 2).ok().as_deref(), Some("0.12"));
/// assert!(number(1e100, 2).is_err());
/// ```
pub fn number(value: f64, digits: u32) -> Result<String> {
    if !value.is_finite() {
        return Ok(printf_nonfinite(value).to_owned());
    }
    if i32::try_from(digits).is_err() {
        return Err(Error::InvalidValue(format!(
            "StringUtils::number digit count {digits} exceeds INT_MAX; the source prints six digits instead"
        )));
    }
    // "0." and 62 digits already overflow the buffer, so refuse before formatting.
    let digit_count = match usize::try_from(digits) {
        Ok(count) if count <= NUMBER_MAX_TEXT_BYTES - 2 => count,
        _ => return Err(number_truncated_error(value, digits)),
    };
    let text = format!("{value:.digit_count$}");
    if text.len() > NUMBER_MAX_TEXT_BYTES {
        return Err(number_truncated_error(value, digits));
    }
    Ok(text)
}

/// `%f` digit counts at or above this only append zeros: no double has more
/// than 1074 fraction digits.
const MAX_EFFECTIVE_FIXED_DIGITS: usize = 1100;

/// `StringUtils::number(d, n)` including its silent truncation, for byte parity.
///
/// Identical to [`number`] wherever that succeeds. Otherwise it returns what the
/// source returns: the first [`NUMBER_MAX_TEXT_BYTES`] bytes of the `%.*f` text,
/// and six fraction digits when `digits` exceeds `i32::MAX`. The source gives no
/// diagnostic, and a cut result shows a wrong magnitude (`number(1e100, 2)` is a
/// 63-digit integer), so prefer [`number`] outside tool paths that must
/// reproduce the C++ output byte for byte.
///
/// The exception is a digit count at which the `%.*f` text would be `INT_MAX`
/// bytes or longer. There `snprintf` overflows its `int` result and the source
/// text depends on the C library: Apple libc returns an empty string, glibc 2.39
/// the digit prefix or, at `digits == INT_MAX`, 63 spaces. This function still
/// returns the digit prefix.
///
/// Work is bounded: a digit count above 1100 is formatted as 1100. That keeps
/// the returned prefix unchanged, because a double has at most 1074 fraction
/// digits and every further digit is a zero.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::number_truncated;
///
/// assert_eq!(number_truncated(6.0 / 60.0, 1), "0.1");
/// assert_eq!(number_truncated(1e100, 2).len(), 63);
/// assert_eq!(number_truncated(0.125, u32::MAX), "0.125000");
/// ```
pub fn number_truncated(value: f64, digits: u32) -> String {
    if !value.is_finite() {
        return printf_nonfinite(value).to_owned();
    }
    let digit_count = printf_precision(digits).min(MAX_EFFECTIVE_FIXED_DIGITS);
    let text = format!("{value:.digit_count$}");
    match text.get(..NUMBER_MAX_TEXT_BYTES) {
        Some(prefix) => prefix.to_owned(),
        None => text,
    }
}

fn number_truncated_error(value: f64, digits: u32) -> Error {
    Error::InvalidValue(format!(
        "StringUtils::number({value:e}, {digits}) needs more than {NUMBER_MAX_TEXT_BYTES} bytes, \
         which the source silently truncates"
    ))
}

/// `10^exponent` as `std::pow(10.0, Int(exponent))` gives it on the reference
/// build: the correctly rounded power.
fn power_of_ten(exponent: i32) -> f64 {
    format!("1e{exponent}").parse().unwrap_or(f64::NAN)
}

/// `StringUtils::numberLength(d, n)`: `value` as text of at most `length`
/// bytes, switching to a hand-built scientific form when the fixed form would
/// not fit.
///
/// The source (`StringUtils.cpp:533-560`) takes the sign `s` (1 when `value < 0`)
/// and the magnitude `d`. When `d < 10^(n - s - 2)` it streams the sign and `d`
/// at stream precision 15 (`%.15g`). Otherwise it divides `d` by 10 until
/// `d <= 10^(n - s - 4)`, counting the divisions, truncates `d` to an `int`,
/// divides that by 10, counts one more, and streams the sign, the result at the
/// default precision 6, `e`, a `0` when the count is below 10, and the count.
/// The exponents are computed in `unsigned int` and converted to `int`, so a
/// short `length` wraps them negative, and the text is finally cut to `length`
/// bytes. The cut is not rounding: `numberLength(1.0, 2)` is `0e`,
/// `numberLength(3.14159265, 7)` is `3.14159`, and a large value keeps its
/// leading digit and exponent (`1.2e012` style texts, truncated).
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] for a non-finite `value` — the source loops
/// forever on an infinity and converts a NaN to `int`, which is undefined — and
/// when the truncated magnitude does not fit an `int`, which is undefined in
/// the source too (it happens for `length >= 14` with values of `1e12` and
/// above; x86 gives `INT_MIN`).
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::number_length;
///
/// assert_eq!(number_length(3.141_592_65, 7).ok().as_deref(), Some("3.14159"));
/// assert_eq!(number_length(1.0, 2).ok().as_deref(), Some("0e"));
/// assert!(number_length(f64::INFINITY, 5).is_err());
/// ```
pub fn number_length(value: f64, length: u32) -> Result<String> {
    if !value.is_finite() {
        return Err(Error::InvalidValue(format!(
            "StringUtils::numberLength({value}, {length}): the source does not terminate or is undefined for a non-finite value"
        )));
    }
    let sign: u32 = u32::from(value < 0.0);
    let mut magnitude = value.abs();
    // Int(n - sign - 2) with UInt arithmetic: wraps, then reinterprets.
    let fixed_exponent = length.wrapping_sub(sign).wrapping_sub(2) as i32;
    let mut text = String::new();
    if sign == 1 {
        text.push('-');
    }
    if magnitude < power_of_ten(fixed_exponent) {
        text.push_str(&ostream_g(magnitude, WRITTEN_DIGITS_F64));
    } else {
        let limit = power_of_ten(length.wrapping_sub(sign).wrapping_sub(4) as i32);
        let mut exponent: u32 = 0;
        while magnitude > limit {
            magnitude /= 10.0;
            exponent = exponent.wrapping_add(1);
        }
        if magnitude.trunc() > f64::from(i32::MAX) {
            return Err(Error::InvalidValue(format!(
                "StringUtils::numberLength({value}, {length}) converts {magnitude} to int, which is undefined in the source"
            )));
        }
        // static_cast<int>(d) / 10.0
        let truncated = f64::from(magnitude.trunc() as i32) / 10.0;
        exponent = exponent.wrapping_add(1);
        text.push_str(&ostream_g(truncated, 6));
        text.push('e');
        if exponent < 10 {
            text.push('0');
        }
        text.push_str(&exponent.to_string());
    }
    let cut = usize::try_from(length)
        .unwrap_or(usize::MAX)
        .min(text.len());
    text.truncate(cut);
    Ok(text)
}

/// `StringUtils::random(length)` with a caller-supplied source of uniform
/// numbers in `[0, 1)`: `length` characters of `[0-9A-Za-z]`.
///
/// Each character is `floor(u * 62)` of the next number `u`, mapped to `0`-`9`,
/// `A`-`Z`, `a`-`z` in that order, as the source maps
/// `rand() / (RAND_MAX + 1.0)`. A number outside `[0, 1)` is clamped into it.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `length` exceeds [`MAX_TEXT_BYTES`] or
/// cannot be allocated.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::random_with;
///
/// let mut values = [0.0, 0.5, 0.999].into_iter();
/// let text = random_with(3, || values.next().unwrap_or(0.0));
/// assert_eq!(text.ok().as_deref(), Some("0Vz"));
/// ```
pub fn random_with(length: u32, mut next_unit: impl FnMut() -> f64) -> Result<String> {
    let length = checked_text_len(u64::from(length))?;
    let mut text = reserve_text(length)?;
    for _ in 0..length {
        let unit = next_unit();
        let unit = if unit.is_nan() {
            0.0
        } else {
            unit.clamp(0.0, 1.0 - f64::EPSILON / 2.0)
        };
        let index = (unit * 62.0).floor() as u8;
        let byte = match index {
            0..=9 => b'0' + index,
            10..=35 => b'A' + (index - 10),
            _ => b'a' + (index.min(61) - 36),
        };
        text.push(char::from(byte));
    }
    Ok(text)
}

/// `StringUtils::random(length)`: `length` random characters of `[0-9A-Za-z]`.
///
/// The source reseeds the C library generator with `srand(time(nullptr))` on
/// every call, so two calls within the same second return the same text. This
/// port draws from a 32-bit Mersenne Twister seeded from the system clock's
/// nanoseconds and the call count, so consecutive calls differ; neither the
/// source's nor this sequence is reproducible, and no caller may rely on
/// either. Use [`random_with`] for a deterministic sequence.
///
/// # Errors
///
/// As [`random_with`].
pub fn random(length: u32) -> Result<String> {
    use std::sync::atomic::{AtomicU32, Ordering};
    use std::time::{SystemTime, UNIX_EPOCH};
    static CALLS: AtomicU32 = AtomicU32::new(0);
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |elapsed| {
            elapsed.subsec_nanos() ^ (elapsed.as_secs() as u32)
        });
    let mut generator =
        rand_mt::Mt::new(nanos ^ CALLS.fetch_add(1, Ordering::Relaxed).rotate_left(16));
    random_with(length, || f64::from(generator.next_u32()) / 4_294_967_296.0)
}

// ---------------------------------------------------------------------------
// Checked allocation helpers
// ---------------------------------------------------------------------------

/// `bytes` as a text length within [`MAX_TEXT_BYTES`].
fn checked_text_len(bytes: u64) -> Result<usize> {
    usize::try_from(bytes)
        .ok()
        .filter(|bytes| *bytes <= MAX_TEXT_BYTES)
        .ok_or_else(|| {
            Error::InvalidValue(format!(
                "a text of {bytes} bytes exceeds the {MAX_TEXT_BYTES}-byte ceiling"
            ))
        })
}

/// An empty `String` with `capacity` reserved fallibly.
fn reserve_text(capacity: usize) -> Result<String> {
    let mut text = String::new();
    text.try_reserve_exact(capacity)
        .map_err(|_| Error::InvalidValue(format!("cannot allocate {capacity} bytes of text")))?;
    Ok(text)
}

/// `c` as the one byte a C++ `char` would hold.
fn ascii(c: char, role: &str) -> Result<u8> {
    u8::try_from(c).ok().filter(u8::is_ascii).ok_or_else(|| {
        Error::InvalidValue(format!(
            "the {role} {c:?} is not ASCII; the source's char is one byte of the text"
        ))
    })
}

/// A byte range of `text` as a `String`, refused when it would cut a
/// multi-byte character.
fn byte_range(text: &str, start: usize, end: usize) -> Result<String> {
    text.get(start..end).map(str::to_owned).ok_or_else(|| {
        Error::InvalidValue(format!(
            "bytes {start}..{end} of a {}-byte text would split a UTF-8 character",
            text.len()
        ))
    })
}

// ---------------------------------------------------------------------------
// Accessors
// ---------------------------------------------------------------------------

fn index_overflow(index: usize, size: usize) -> Error {
    Error::InvalidRange(format!(
        "the given index was too large: {index} (size = {size})"
    ))
}

fn index_underflow(index: i32) -> Error {
    Error::InvalidRange(format!("the given index was too small: {index} (size = 0)"))
}

/// `StringUtils::prefix(s, size_t length)`: the first `length` bytes.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] with the source's `IndexOverflow` text ("the
/// given index was too large: 10 (size = 3)") when `length` exceeds the text,
/// and [`Error::InvalidValue`] when the cut would split a UTF-8 character.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::prefix;
///
/// assert_eq!(prefix("abc", 2).ok().as_deref(), Some("ab"));
/// assert!(prefix("abc", 4).is_err());
/// ```
pub fn prefix(text: &str, length: usize) -> Result<String> {
    if length > text.len() {
        return Err(index_overflow(length, text.len()));
    }
    byte_range(text, 0, length)
}

/// `StringUtils::suffix(s, size_t length)`: the last `length` bytes.
///
/// # Errors
///
/// As [`prefix`].
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::suffix;
///
/// assert_eq!(suffix("abc", 2).ok().as_deref(), Some("bc"));
/// ```
pub fn suffix(text: &str, length: usize) -> Result<String> {
    if length > text.len() {
        return Err(index_overflow(length, text.len()));
    }
    byte_range(text, text.len() - length, text.len())
}

/// `StringUtils::prefix(s, Int length)`: [`prefix`] for a signed length.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] with the source's `IndexUnderflow` text ("the
/// given index was too small: -1 (size = 0)") for a negative length, otherwise
/// as [`prefix`].
pub fn prefix_signed(text: &str, length: i32) -> Result<String> {
    let length = usize::try_from(length).map_err(|_| index_underflow(length))?;
    prefix(text, length)
}

/// `StringUtils::suffix(s, Int length)`: [`suffix`] for a signed length.
///
/// # Errors
///
/// As [`prefix_signed`].
pub fn suffix_signed(text: &str, length: i32) -> Result<String> {
    let length = usize::try_from(length).map_err(|_| index_underflow(length))?;
    suffix(text, length)
}

/// `StringUtils::prefix(s, char delim)`: the text before the first `delimiter`,
/// or the whole text when there is none.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::prefix_until;
///
/// assert_eq!(prefix_until("x,y,z", ','), "x");
/// assert_eq!(prefix_until("xyz", ','), "xyz");
/// ```
pub fn prefix_until(text: &str, delimiter: char) -> String {
    text.split(delimiter).next().unwrap_or(text).to_owned()
}

/// `StringUtils::suffix(s, char delim)`: the text after the last `delimiter`,
/// or the whole text when there is none.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::suffix_after;
///
/// assert_eq!(suffix_after("x,y,z", ','), "z");
/// ```
pub fn suffix_after(text: &str, delimiter: char) -> String {
    text.rsplit(delimiter).next().unwrap_or(text).to_owned()
}

/// `StringUtils::substr(s, pos, n)`: `std::string::substr` with the start
/// clamped to the text, so a start beyond the end gives an empty text instead
/// of `std::out_of_range`. `None` is `npos`, the rest of the text; a length
/// beyond the end is clamped as `substr` clamps it.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a boundary would split a UTF-8
/// character.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::substr;
///
/// assert_eq!(substr("abc", 1, None).ok().as_deref(), Some("bc"));
/// assert_eq!(substr("abc", 5, Some(2)).ok().as_deref(), Some(""));
/// ```
pub fn substr(text: &str, position: usize, length: Option<usize>) -> Result<String> {
    let start = position.min(text.len());
    let end = match length {
        Some(length) => start.saturating_add(length).min(text.len()),
        None => text.len(),
    };
    byte_range(text, start, end)
}

/// `StringUtils::chop(s, n)`: the text without its last `count` bytes, empty
/// when `count` is at least the length.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the cut would split a UTF-8 character.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::chop;
///
/// assert_eq!(chop("abc", 1).ok().as_deref(), Some("ab"));
/// assert_eq!(chop("abc", 5).ok().as_deref(), Some(""));
/// ```
pub fn chop(text: &str, count: usize) -> Result<String> {
    byte_range(text, 0, text.len().saturating_sub(count))
}

// ---------------------------------------------------------------------------
// In-place mutators
// ---------------------------------------------------------------------------

/// `StringUtils::trim`: remove leading and trailing space, tab, line feed and
/// carriage return, in place.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::trim;
///
/// let mut text = String::from(" \t a b \r\n");
/// trim(&mut text);
/// assert_eq!(text, "a b");
/// ```
pub fn trim(text: &mut String) {
    let trimmed = text.trim_matches(blank_char);
    if trimmed.len() != text.len() {
        *text = trimmed.to_owned();
    }
}

/// `StringUtils::trimmed`: a trimmed copy (see [`trim`]).
pub fn trimmed(text: &str) -> String {
    text.trim_matches(blank_char).to_owned()
}

/// `StringUtils::toUpper`: ASCII letters to upper case in place, as `toupper`
/// in the classic `"C"` locale does; every other byte is unchanged.
pub fn to_upper(text: &mut str) {
    text.make_ascii_uppercase();
}

/// `StringUtils::toLower`: ASCII letters to lower case in place (see
/// [`to_upper`]).
pub fn to_lower(text: &mut str) {
    text.make_ascii_lowercase();
}

/// `StringUtils::toUppered`: an upper-cased copy (see [`to_upper`]).
pub fn to_uppered(text: &str) -> String {
    text.to_ascii_uppercase()
}

/// `StringUtils::toLowered`: a lower-cased copy (see [`to_lower`]).
pub fn to_lowered(text: &str) -> String {
    text.to_ascii_lowercase()
}

/// `StringUtils::firstToUpper`: the first byte to upper case when it is an ASCII
/// letter.
pub fn first_to_upper(text: &mut str) {
    if let Some(first) = text.get_mut(..1) {
        first.make_ascii_uppercase();
    }
}

/// `StringUtils::reverse`: the bytes in reverse order, in place.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`], leaving the text unchanged, when it is not
/// ASCII: reversing the bytes of a multi-byte character, as the source does,
/// does not give text.
pub fn reverse(text: &mut String) -> Result<()> {
    if !text.is_ascii() {
        return Err(Error::InvalidValue(
            "reversing the bytes of non-ASCII text, as the source does, does not give UTF-8".into(),
        ));
    }
    *text = text.chars().rev().collect();
    Ok(())
}

/// `StringUtils::simplify`: every run of space, tab, line feed and carriage
/// return becomes one space, in place. The ends are not trimmed.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::simplify;
///
/// let mut text = String::from(" a \t\n b ");
/// simplify(&mut text);
/// assert_eq!(text, " a b ");
/// ```
pub fn simplify(text: &mut String) {
    let mut result = String::with_capacity(text.len());
    let mut last_blank = false;
    for c in text.chars() {
        if blank_char(c) {
            if !last_blank {
                result.push(' ');
            }
            last_blank = true;
        } else {
            result.push(c);
            last_blank = false;
        }
    }
    *text = result;
}

/// Pad `text` with `fill` to `size` bytes on one side.
fn fill(text: &mut String, fill: char, size: u32, left: bool) -> Result<()> {
    let fill = char::from(ascii(fill, "fill character")?);
    let size = checked_text_len(u64::from(size))?;
    let Some(missing) = size.checked_sub(text.len()).filter(|missing| *missing > 0) else {
        return Ok(());
    };
    let mut result = reserve_text(size)?;
    if left {
        result.extend(std::iter::repeat_n(fill, missing));
        result.push_str(text);
    } else {
        result.push_str(text);
        result.extend(std::iter::repeat_n(fill, missing));
    }
    *text = result;
    Ok(())
}

/// `StringUtils::fillLeft`: prepend `fill` until the text is `size` bytes long;
/// a longer text is unchanged. Covers the rvalue overload too: fill an owned
/// `String` in place.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`], leaving the text unchanged, when `fill` is
/// not ASCII or `size` exceeds [`MAX_TEXT_BYTES`].
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::fill_left;
///
/// let mut text = String::from("7");
/// fill_left(&mut text, '0', 3)?;
/// assert_eq!(text, "007");
/// # Ok::<(), openms::Error>(())
/// ```
pub fn fill_left(text: &mut String, fill_char: char, size: u32) -> Result<()> {
    fill(text, fill_char, size, true)
}

/// `StringUtils::fillRight`: append `fill` until the text is `size` bytes long.
///
/// # Errors
///
/// As [`fill_left`].
pub fn fill_right(text: &mut String, fill_char: char, size: u32) -> Result<()> {
    fill(text, fill_char, size, false)
}

/// `StringUtils::substitute(s, char, char)`: replace every `from` with `to`, in
/// place.
pub fn substitute_char(text: &mut String, from: char, to: char) {
    if text.contains(from) {
        *text = text
            .chars()
            .map(|c| if c == from { to } else { c })
            .collect();
    }
}

/// `StringUtils::substituted(s, char, char)`: [`substitute_char`] on a copy.
pub fn substituted_char(text: &str, from: char, to: char) -> String {
    text.chars()
        .map(|c| if c == from { to } else { c })
        .collect()
}

/// `StringUtils::substitute(s, string, string)`: replace every non-overlapping
/// occurrence of `from`, left to right, with `to`, in place. An empty `from`
/// changes nothing.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`], leaving the text unchanged, when the result
/// would exceed [`MAX_TEXT_BYTES`] or cannot be allocated.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::substitute;
///
/// let mut text = String::from("aaa");
/// substitute(&mut text, "aa", "b")?;
/// assert_eq!(text, "ba");
/// # Ok::<(), openms::Error>(())
/// ```
pub fn substitute(text: &mut String, from: &str, to: &str) -> Result<()> {
    if from.is_empty() || !text.contains(from) {
        return Ok(());
    }
    *text = substituted(text, from, to)?;
    Ok(())
}

/// `StringUtils::substituted(s, string, string)`: [`substitute`] on a copy.
///
/// # Errors
///
/// As [`substitute`].
pub fn substituted(text: &str, from: &str, to: &str) -> Result<String> {
    if from.is_empty() {
        return Ok(text.to_owned());
    }
    let count = text.matches(from).count() as u64;
    let length = (text.len() as u64)
        .checked_sub(count.saturating_mul(from.len() as u64))
        .and_then(|rest| rest.checked_add(count.checked_mul(to.len() as u64)?))
        .unwrap_or(u64::MAX);
    let length = checked_text_len(length)?;
    let mut result = reserve_text(length)?;
    let mut start = 0;
    for (position, _) in text.match_indices(from) {
        result.push_str(text.get(start..position).unwrap_or_default());
        result.push_str(to);
        start = position + from.len();
    }
    result.push_str(text.get(start..).unwrap_or_default());
    Ok(result)
}

/// `StringUtils::remove(s, char)`: delete every `what`, in place.
pub fn remove(text: &mut String, what: char) {
    text.retain(|c| c != what);
}

/// `StringUtils::ensureLastChar`: append `end` unless the text already ends
/// with it (an empty text always gets it). Covers the rvalue overload too.
pub fn ensure_last_char(text: &mut String, end: char) {
    if !text.ends_with(end) {
        text.push(end);
    }
}

/// `StringUtils::removeWhitespaces`: delete every space, tab, line feed and
/// carriage return, in place.
pub fn remove_whitespaces(text: &mut String) {
    text.retain(|c| !blank_char(c));
}

/// `StringUtils::isQuoted`: at least two characters, the first and the last
/// being `quote`.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::is_quoted;
///
/// assert!(is_quoted("\"a\"", '"'));
/// assert!(!is_quoted("\"", '"'));
/// ```
pub fn is_quoted(text: &str, quote: char) -> bool {
    let mut chars = text.chars();
    chars.next() == Some(quote) && chars.next_back() == Some(quote)
}

/// `StringUtils::quote`: mark embedded quotes by `method`, then surround with
/// `quote`, in place.
///
/// `Escape` first doubles every backslash and then puts a backslash before
/// every `quote`; `Double` doubles every `quote`; `None` marks nothing.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`], leaving the text unchanged, when `quote` is
/// not ASCII or the result would exceed [`MAX_TEXT_BYTES`].
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::{quote, QuotingMethod};
///
/// let mut text = String::from(r#"a"b\"#);
/// quote(&mut text, '"', QuotingMethod::Escape)?;
/// assert_eq!(text, r#""a\"b\\""#);
/// # Ok::<(), openms::Error>(())
/// ```
pub fn quote(text: &mut String, quote: char, method: QuotingMethod) -> Result<()> {
    let q = char::from(ascii(quote, "quote character")?);
    checked_text_len((text.len() as u64).saturating_mul(2).saturating_add(2))?;
    let mut result = match method {
        QuotingMethod::Escape => {
            let escaped = substituted(text, "\\", "\\\\")?;
            substituted(&escaped, &q.to_string(), &format!("\\{q}"))?
        }
        QuotingMethod::Double => substituted(text, &q.to_string(), &format!("{q}{q}"))?,
        QuotingMethod::None => text.clone(),
    };
    result.insert(0, q);
    result.push(q);
    *text = result;
    Ok(())
}

/// `StringUtils::unquote`: strip the surrounding `quote` characters and undo
/// `method`'s marking, in place.
///
/// `Escape` replaces backslash-quote by the quote and then a double backslash by
/// one, in that order; `Double` replaces a doubled quote by one.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] with the source's `ConversionError` text
/// ("'…' does not have the expected format of a quoted string") when
/// [`is_quoted`] is false, and when `quote` is not ASCII. The text is unchanged
/// on error.
pub fn unquote(text: &mut String, quote: char, method: QuotingMethod) -> Result<()> {
    let q = char::from(ascii(quote, "quote character")?);
    if !is_quoted(text, q) {
        return Err(Error::InvalidValue(format!(
            "'{text}' does not have the expected format of a quoted string"
        )));
    }
    let inner = text.get(1..text.len() - 1).unwrap_or_default();
    let result = match method {
        QuotingMethod::Escape => {
            let unescaped = substituted(inner, &format!("\\{q}"), &q.to_string())?;
            substituted(&unescaped, "\\\\", "\\")?
        }
        QuotingMethod::Double => substituted(inner, &format!("{q}{q}"), &q.to_string())?,
        QuotingMethod::None => inner.to_owned(),
    };
    *text = result;
    Ok(())
}

// ---------------------------------------------------------------------------
// Split / join
// ---------------------------------------------------------------------------

/// The result of a split: the parts, and the source's return value, whether
/// the text was split at all.
///
/// `split` is not always `parts.len() > 1`: splitting a one-character text on
/// an empty splitter returns one part and `true`, as the source does.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Split {
    /// The substrings, in order.
    pub parts: Vec<String>,
    /// The source's `bool` return value.
    pub split: bool,
}

fn check_parts(count: usize) -> Result<()> {
    if count > MAX_SPLIT_PARTS {
        return Err(Error::InvalidValue(format!(
            "a split into {count} parts exceeds the {MAX_SPLIT_PARTS}-part ceiling"
        )));
    }
    Ok(())
}

/// `trim` and the quote check of `split`'s quote-protected mode.
fn dequote_block(block: &str) -> Result<String> {
    let block = block.trim_matches(blank_char);
    let has_start = block.starts_with('"');
    let has_end = block.ends_with('"');
    if block.len() >= 2 && (has_start ^ has_end) {
        return Err(Error::InvalidValue(format!(
            "Could not dequote string '{block}' due to wrongly placed '\"'."
        )));
    }
    if block.len() >= 2 && has_start && has_end {
        return Ok(block.get(1..block.len() - 1).unwrap_or_default().to_owned());
    }
    Ok(block.to_owned())
}

/// `StringUtils::split(s, char splitter, substrings, quote_protect)`.
///
/// An empty text gives no parts and `false`. Without `quote_protect`, the text
/// is cut at every `splitter` (empty parts kept) and the result is `true` when
/// there was at least one. With `quote_protect`, a `splitter` counts only
/// outside `"`-quoted blocks (an even number of `"` before it); every part is
/// trimmed and loses one pair of surrounding `"`; when no `splitter` counts,
/// the result is the untouched text and `false`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `splitter` is not ASCII, when a
/// quote-protected part of two or more characters has a `"` at only one end
/// (the source's `ConversionError` "Could not dequote string '…' due to wrongly
/// placed '\"'."), and when the parts would exceed [`MAX_SPLIT_PARTS`].
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::split;
///
/// let result = split("x,y,,z", ',', false)?;
/// assert_eq!(result.parts, ["x", "y", "", "z"]);
/// assert!(result.split);
/// let result = split(r#" "a,b" , c"#, ',', true)?;
/// assert_eq!(result.parts, ["a,b", "c"]);
/// # Ok::<(), openms::Error>(())
/// ```
pub fn split(text: &str, splitter: char, quote_protect: bool) -> Result<Split> {
    let splitter = ascii(splitter, "splitter")?;
    if text.is_empty() {
        return Ok(Split::default());
    }
    let bytes = text.as_bytes();
    let count = bytes.iter().filter(|b| **b == splitter).count();
    check_parts(count.saturating_add(1))?;
    if !quote_protect && count == 0 {
        return Ok(Split {
            parts: vec![text.to_owned()],
            split: false,
        });
    }
    let mut parts = Vec::new();
    parts
        .try_reserve(count + 1)
        .map_err(|_| Error::InvalidValue("cannot allocate the split parts".into()))?;
    if quote_protect {
        let mut quotes = 0usize;
        let mut begin = 0usize;
        for (index, byte) in bytes.iter().enumerate() {
            if *byte == b'"' {
                quotes += 1;
            }
            if quotes % 2 == 0 && *byte == splitter {
                parts.push(dequote_block(text.get(begin..index).unwrap_or_default())?);
                begin = index + 1;
            }
        }
        if parts.is_empty() {
            return Ok(Split {
                parts: vec![text.to_owned()],
                split: false,
            });
        }
        parts.push(dequote_block(text.get(begin..).unwrap_or_default())?);
    } else {
        parts.extend(text.split(char::from(splitter)).map(str::to_owned));
    }
    Ok(Split { parts, split: true })
}

/// `StringUtils::split(s, const std::string& splitter, substrings)`.
///
/// An empty text gives no parts and `false`. An empty `splitter` cuts the text
/// into single characters and returns `true` (the source cuts into single
/// bytes, which is the same for ASCII). Otherwise the text is cut at every
/// non-overlapping `splitter`, left to right, and the result is `true` when
/// there is more than one part.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the parts would exceed
/// [`MAX_SPLIT_PARTS`].
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::split_str;
///
/// let result = split_str("::a::b::", "::")?;
/// assert_eq!(result.parts, ["", "a", "b", ""]);
/// # Ok::<(), openms::Error>(())
/// ```
pub fn split_str(text: &str, splitter: &str) -> Result<Split> {
    if text.is_empty() {
        return Ok(Split::default());
    }
    if splitter.is_empty() {
        check_parts(text.chars().count())?;
        return Ok(Split {
            parts: text.chars().map(String::from).collect(),
            split: true,
        });
    }
    check_parts(text.matches(splitter).count().saturating_add(1))?;
    let parts: Vec<String> = text.split(splitter).map(str::to_owned).collect();
    let split = parts.len() > 1;
    Ok(Split { parts, split })
}

/// `StringUtils::split_quoted`: cut at every `splitter` outside quoted blocks.
///
/// A block opens at a `quote` and closes at the next `quote` that `method` does
/// not mark as embedded (`Escape`: a quote after an odd run of backslashes;
/// `Double`: a doubled quote; `None`: none). Quotes stay in the parts. An empty
/// text or an empty `splitter` gives no parts and `false`; otherwise the result
/// is `true` when there is more than one part.
///
/// The source detects an unbalanced quote only when the unmatched `quote` is the
/// **last** character; an unmatched quote earlier silently extends its block to
/// the end of the text (`"a,b` is one part). This port does the same.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] with the source's `ConversionError` text
/// ("unbalanced quotation marks in string '…'") in that one case, when `quote`
/// is not ASCII, and when the parts would exceed [`MAX_SPLIT_PARTS`].
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::{split_quoted, QuotingMethod};
///
/// let result = split_quoted(r#"a,"b,c",d"#, ",", '"', QuotingMethod::Escape)?;
/// assert_eq!(result.parts, ["a", r#""b,c""#, "d"]);
/// # Ok::<(), openms::Error>(())
/// ```
pub fn split_quoted(
    text: &str,
    splitter: &str,
    quote: char,
    method: QuotingMethod,
) -> Result<Split> {
    let q = ascii(quote, "quote character")?;
    let bytes = text.as_bytes();
    let pattern = splitter.as_bytes();
    let Some((&first, rest)) = pattern.split_first() else {
        return Ok(Split::default());
    };
    if bytes.is_empty() {
        return Ok(Split::default());
    }
    let mut parts = Vec::new();
    let mut in_quote = false;
    let mut start = 0usize;
    let mut i = 0usize;
    while i < bytes.len() {
        if in_quote {
            let mut embedded = false;
            match method {
                QuotingMethod::Escape => {
                    while i < bytes.len() {
                        if bytes[i] == b'\\' {
                            embedded = !embedded;
                        } else if bytes[i] == q && !embedded {
                            break;
                        } else {
                            embedded = false;
                        }
                        i += 1;
                    }
                }
                QuotingMethod::None | QuotingMethod::Double => {
                    while i < bytes.len() {
                        if bytes[i] == q {
                            if method == QuotingMethod::None {
                                break;
                            }
                            if bytes.get(i + 1) == Some(&q) {
                                embedded = !embedded;
                            } else if !embedded {
                                break;
                            } else {
                                embedded = false;
                            }
                        }
                        i += 1;
                    }
                }
            }
            in_quote = false;
        } else {
            // s.find_first_of({q, splitter[0]}, i)
            let Some(offset) = bytes
                .get(i..)
                .and_then(|tail| tail.iter().position(|b| *b == q || *b == first))
            else {
                break;
            };
            i += offset;
            if bytes[i] == q {
                in_quote = true;
            } else if bytes.get(i + 1..i + 1 + rest.len()) == Some(rest) {
                check_parts(parts.len() + 2)?;
                parts.push(byte_range(text, start, i)?);
                start = i + pattern.len();
                i = start;
                continue;
            }
        }
        i += 1;
    }
    if in_quote {
        return Err(Error::InvalidValue(format!(
            "unbalanced quotation marks in string '{text}'"
        )));
    }
    parts.push(byte_range(text, start.min(bytes.len()), bytes.len())?);
    let split = parts.len() > 1;
    Ok(Split { parts, split })
}

/// `StringUtils::concatenate`: the items joined with `glue` between them (empty
/// for no items). Covers the target, iterator-pair and container overloads.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the result would exceed
/// [`MAX_TEXT_BYTES`] or cannot be allocated.
///
/// # Examples
///
/// ```
/// use openms::data_structures::string_utils::concatenate;
///
/// assert_eq!(concatenate(["a", "", "b"], ", ").ok().as_deref(), Some("a, , b"));
/// assert_eq!(concatenate(Vec::<String>::new(), "x").ok().as_deref(), Some(""));
/// ```
pub fn concatenate<I>(items: I, glue: &str) -> Result<String>
where
    I: IntoIterator,
    I::Item: AsRef<str>,
{
    let items: Vec<I::Item> = items.into_iter().collect();
    let mut length = 0u64;
    for (index, item) in items.iter().enumerate() {
        if index > 0 {
            length = length.saturating_add(glue.len() as u64);
        }
        length = length.saturating_add(item.as_ref().len() as u64);
    }
    let mut result = reserve_text(checked_text_len(length)?)?;
    for (index, item) in items.iter().enumerate() {
        if index > 0 {
            result.push_str(glue);
        }
        result.push_str(item.as_ref());
    }
    Ok(result)
}
