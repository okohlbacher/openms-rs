// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Locale-independent decimal text for floating-point values:
//! `CONCEPT/Detail/NumericFormatting.h` (in `src/common/include` at core
//! `bc9cc12`).
//!
//! `Internal::NumericFormatting::appendNumeric<T>(value, target, precision,
//! fixed_format)` is the one routine behind every `StringUtils::toStr` of a
//! floating-point value, and so behind `PrecisionWrapper`, `DPosition` output
//! and the `ListUtilsIO` vector operators. It lives in `concept` so that every
//! module can reach it; `format::file_info::text_format` builds its public
//! `toStr` family on it and holds the executed differential evidence
//! (`tests/file_info_text_format.rs`).
//!
//! # API mapping
//!
//! | Source | Rust |
//! |---|---|
//! | `appendNumeric<double>(double, std::string&, int, bool)` | [`append_numeric_f64`](crate::concept::numeric_formatting::append_numeric_f64) |
//! | `appendNumeric<float>(float, std::string&, int, bool)` | [`append_numeric_f32`](crate::concept::numeric_formatting::append_numeric_f32) |
//! | `appendNumeric<long double>(...)` | not ported: Rust has no `long double` |
//! | the local `char buf[64]` | [`NUMERIC_BUFFER_BYTES`](crate::concept::numeric_formatting::NUMERIC_BUFFER_BYTES) |
//!
//! # Rules
//!
//! - NaN writes `NaN` whatever its sign bit; the infinities write `inf` and
//!   `-inf`.
//! - Zero and magnitudes in `[1e-2, 1e4)` (thresholds of the value's own type)
//!   use fixed notation with `precision` fraction digits,
//!   `std::to_chars(fixed, precision)`.
//! - Other magnitudes use scientific notation: with `fixed_format`, three
//!   fraction digits (`std::to_chars(scientific, 3)`, whatever `precision`
//!   says); without it, the shortest text that reads back as the same value
//!   (`std::to_chars(scientific)`), the nearest one of equal length, an exact
//!   tie taking the even digit.
//! - Rounding is of the exact binary value with ties to even.
//! - Trailing fraction zeros are removed but one fraction digit stays (`5.0`,
//!   `1.0e06`); a text without a point gains `.0`. The exponent loses its `+`
//!   and keeps at least two digits (`1.0e-05`, `3.6739e04`, `1.0e300`).
//! - Fixed text longer than the source's 64-byte buffer makes `to_chars` fail,
//!   and the source then appends `std::to_string(static_cast<double>(value))`,
//!   which is `%f`: six fraction digits, untrimmed. Only a `precision` above
//!   about 57 gets there. The fallback is reproduced, and a `precision` that
//!   cannot fit is detected before any text is formatted, so a huge
//!   `precision` costs no work.

/// Size of the source's `char buf[64]` that `std::to_chars` writes into.
pub const NUMERIC_BUFFER_BYTES: usize = 64;

/// `std::to_chars(first, last, value, std::chars_format::scientific)`, in Rust's
/// `{:e}` spelling: the shortest text that reads back as `value`.
///
/// Among the equally short texts that read back, the C++ standard
/// (`[charconv.to.chars]`) takes the one nearest the exact binary value and
/// breaks a remaining tie by round-to-nearest, which gives the even digit.
/// `{:e}` finds the same shortest length but breaks such an exact tie upwards
/// (`-135169.706298828125` gives `...813e5` where `to_chars` gives `...812e+05`).
///
/// Both equally short neighbours of the exact value are candidates, and the
/// nearer one is `{:.N$e}` at that length, which rounds the exact value with ties
/// to even. It is taken whenever it reads back as `value`. When it does not, only
/// the other neighbour reads back, and `{:e}` already chose that one. This
/// happens where the value is a power of two, whose lower half-interval is half as
/// wide as the upper: `2^-24` is `5.9604644775390625e-8` exactly, and only
/// `5.960464477539063e-8` reads back; the `float` `2^-96` is nearer to
/// `1.2621774e-29`, but only `1.2621775e-29` reads back.
fn shortest_round_trip<T>(value: T) -> String
where
    T: Copy + PartialEq + std::fmt::LowerExp + std::str::FromStr,
{
    let shortest = format!("{value:e}");
    let significant = shortest.split_once('e').map_or(0, |(mantissa, _)| {
        mantissa.bytes().filter(u8::is_ascii_digit).count()
    });
    let Some(fraction) = significant.checked_sub(1) else {
        return shortest;
    };
    let nearest = format!("{value:.fraction$e}");
    if nearest != shortest && nearest.parse::<T>().is_ok_and(|parsed| parsed == value) {
        nearest
    } else {
        shortest
    }
}

/// `appendNumeric`'s trim: trailing fraction zeros removed, one digit kept.
fn keep_one_fraction_digit(text: &str) -> &str {
    let Some(point) = text.find('.') else {
        return text;
    };
    let trimmed = text.trim_end_matches('0');
    if trimmed.len() > point + 1 {
        trimmed
    } else {
        text.get(..point + 2).unwrap_or(text)
    }
}

/// Appends a trimmed mantissa, which gains `.0` when it has no point.
fn push_with_fraction(text: &str, target: &mut String) {
    if text.contains('.') {
        target.push_str(keep_one_fraction_digit(text));
    } else {
        target.push_str(text);
        target.push_str(".0");
    }
}

/// Appends `e`, then the exponent without `+` and with at least two digits.
fn push_exponent(exponent: &str, target: &mut String) {
    target.push('e');
    match exponent.parse::<i32>() {
        Ok(exponent) if exponent < 0 => {
            target.push_str(&format!("-{:02}", exponent.unsigned_abs()));
        }
        Ok(exponent) => target.push_str(&format!("{exponent:02}")),
        Err(_) => target.push_str(exponent),
    }
}

/// `static_cast<double>` of a `double`.
fn widen_f64(value: f64) -> f64 {
    value
}

/// `static_cast<double>` of a `float`, which is exact.
fn widen_f32(value: f32) -> f64 {
    f64::from(value)
}

macro_rules! append_numeric {
    ($name:ident, $float:ty, $widen:ident) => {
        fn $name(value: $float, target: &mut String, precision: usize, fixed_format: bool) {
            if value.is_nan() {
                target.push_str("NaN");
                return;
            }
            if value.is_infinite() {
                target.push_str(if value < 0.0 { "-inf" } else { "inf" });
                return;
            }
            let magnitude = value.abs();
            // Source: abs_val >= T(1e4) || abs_val < T(1e-2); NaN returned above.
            if magnitude != 0.0 && !(1e-2..1e4).contains(&magnitude) {
                let text = if fixed_format {
                    // std::to_chars(scientific, 3): exact value, ties to even.
                    format!("{value:.3e}")
                } else {
                    // std::to_chars(scientific) without a precision.
                    shortest_round_trip(value)
                };
                let Some((mantissa, exponent)) = text.split_once('e') else {
                    target.push_str(&text);
                    return;
                };
                push_with_fraction(mantissa, target);
                push_exponent(exponent, target);
                return;
            }
            // std::to_chars(fixed, precision) into char[64]. The text is at
            // least a digit, a point and `precision` digits, so a precision
            // above 62 cannot fit and is not formatted at all.
            let text = if precision <= NUMERIC_BUFFER_BYTES - 2 {
                Some(format!("{value:.precision$}"))
                    .filter(|text| text.len() <= NUMERIC_BUFFER_BYTES)
            } else {
                None
            };
            match text {
                Some(text) => push_with_fraction(&text, target),
                // to_chars fails; the source appends
                // std::to_string(static_cast<double>(value)), i.e. "%f".
                None => target.push_str(&format!("{:.6}", $widen(value))),
            }
        }
    };
}

append_numeric!(append_numeric_double, f64, widen_f64);
append_numeric!(append_numeric_float, f32, widen_f32);

/// Appends `value` as decimal text: source
/// `NumericFormatting::appendNumeric<double>(value, target, precision,
/// fixed_format)` (`NumericFormatting.h:27-140`).
///
/// See the module notes for the rules. `StringUtils::toStr(double, true)` is
/// `precision` 15 (`writtenDigits<double>()`) without `fixed_format`, and
/// `toStr(double, false)` is `precision` 3 with it.
///
/// # Examples
///
/// ```
/// use openms::concept::numeric_formatting::append_numeric_f64;
///
/// let mut text = String::from("x = ");
/// append_numeric_f64(17.012345, &mut text, 15, false);
/// assert_eq!(text, "x = 17.012345");
/// let mut text = String::new();
/// append_numeric_f64(1.0e6, &mut text, 3, true);
/// assert_eq!(text, "1.0e06");
/// ```
pub fn append_numeric_f64(value: f64, target: &mut String, precision: usize, fixed_format: bool) {
    append_numeric_double(value, target, precision, fixed_format);
}

/// Appends `value` as decimal text: source
/// `NumericFormatting::appendNumeric<float>(value, target, precision,
/// fixed_format)`.
///
/// The magnitude thresholds `1e4f` and `1e-2f` are `float` values and rounding
/// is of the exact `float` value, not of the promoted `double`.
/// `StringUtils::toStr(float, true)` is `precision` 6 (`writtenDigits<float>()`)
/// without `fixed_format`.
pub fn append_numeric_f32(value: f32, target: &mut String, precision: usize, fixed_format: bool) {
    append_numeric_float(value, target, precision, fixed_format);
}
