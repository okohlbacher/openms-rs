// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Locale-independent floating-point text: `CONCEPT/Detail/NumericFormatting.h`
//! and the C library conversions the source's text relies on.
//!
//! [`append_numeric_f64`] and [`append_numeric_f32`] port
//! `Internal::NumericFormatting::appendNumeric<T>`, the one routine behind every
//! `StringUtils::toStr` of a floating-point value, both precisions
//! (`toStr(x, true)` and `toStr(x, false)`, `appendToStr` and
//! `appendToStrLowP`), and therefore behind `DataValue` and `ParamValue` text
//! and every XML writer that stores a number through them. [`ostream_g`] is the
//! C++ default-stream insertion of a `double` (`printf("%.*g")`), which
//! `StringUtils::numberLength` and FileInfo's statistics use.
//!
//! This is the numeric core of `crate::data_structures::string_utils` and of
//! `crate::format::file_info::text_format`, which re-export it. It lives in
//! `concept`, as `NumericFormatting.h` lives in `CONCEPT/Detail`, so that every
//! module can reach it without closing a dependency cycle.
//!
//! Sources at core revision `bc9cc12`:
//! `src/common/include/OpenMS/CONCEPT/Detail/NumericFormatting.h` (141 lines,
//! header-only) and the standard-library behaviour it calls: `std::to_chars`
//! (`<charconv>`) and, for [`ostream_g`], `std::ostream::operator<<(double)`.
//! The support document is `docs/STRING_UTILS_SUPPORT.md`.
//!
//! # API mapping
//!
//! | Source | Rust | Notes |
//! |---|---|---|
//! | `appendNumeric<double>(T, std::string&, int precision, bool fixed_format)` | [`append_numeric_f64`] | `precision` is a `u32`; the source only passes 3, 6, 15 and 18 |
//! | `appendNumeric<float>(...)` | [`append_numeric_f32`] | formats the `float` itself, not a promoted `double` |
//! | `appendNumeric<long double>(...)` | not ported | Rust has no extended-precision type; see the support document |
//! | `std::ostream << double` after `precision(p)` | [`ostream_g`] | C `%.*g`, glibc spelling of NaN |
//!
//! # The `appendNumeric` rules
//!
//! - NaN is `NaN` whatever its sign bit; the infinities are `inf` and `-inf`.
//! - Zero and magnitudes in `[1e-2, 1e4)` (the bounds are `T(1e-2)` and
//!   `T(1e4)`, so a `float` compares against the `float` nearest `0.01`) are
//!   written in fixed notation with exactly `precision` fraction digits, rounded
//!   exactly with ties to even. Trailing zeros are removed but one fraction
//!   digit stays (`5.0`, `-0.0`, `0.062` for `0.0625` at precision 3). A text
//!   longer than the source's 64-byte buffer makes `std::to_chars` fail, and the
//!   source then appends `std::to_string(double(value))`, `%f` with six digits
//!   and no trimming; that fallback is reproduced.
//! - Other magnitudes are written in scientific notation: with
//!   `fixed_format = false` as the shortest text that reads back as the same
//!   value (`std::to_chars` without a precision), with `fixed_format = true` as
//!   `std::to_chars(..., scientific, 3)` — three fraction digits whatever
//!   `precision` says. The mantissa keeps at least one fraction digit (`1.0e06`)
//!   and the exponent has no `+` and at least two digits (`1.235e04`,
//!   `5.0e-04`, `4.941e-324`).
//!
//! # Evidence
//!
//! Tier 1: `../oracle/B8a-stringutils` (Linux x86_64 Release install on
//! ibminode06, two identical runs) records `toStr(double|float, true|false)`,
//! `appendToStrLowP` and `std::string + double|float` for 58 doubles and 72
//! floats; `tests/string_utils.rs` compares every row. The shortest-text tie
//! rule is measured by `../oracle/text-format` and `../oracle/a2-textfmt-linux`
//! (`tests/file_info_text_format.rs`).

use std::fmt::LowerExp;
use std::str::FromStr;

/// Size of the `char buf[64]` that `appendNumeric` formats into.
///
/// A fixed-notation text longer than this makes `std::to_chars` fail with
/// `value_too_large`, and the source then falls back to
/// `std::to_string(double(value))`.
pub const NUMERIC_BUFFER_BYTES: usize = 64;

/// `%g` precisions at or above this give the same text: no double has more than
/// 767 significant digits, and no exponent reaches it.
const MAX_EFFECTIVE_G_PRECISION: usize = 800;
/// Precision `printf` applies when the requested one is negative.
const OMITTED_PRECISION: usize = 6;

/// Generates the `appendNumeric` instantiation for one floating-point type,
/// the C++ template `NumericFormatting::appendNumeric<T>`.
macro_rules! append_numeric {
    ($(#[$doc:meta])* $name:ident, $float:ty) => {
        $(#[$doc])*
        pub fn $name(value: $float, target: &mut String, precision: u32, fixed_format: bool) {
            if value.is_nan() {
                target.push_str("NaN");
                return;
            }
            if value.is_infinite() {
                target.push_str(if value < 0.0 { "-inf" } else { "inf" });
                return;
            }
            // Source: abs_val = value < 0 ? -value : value, so -0.0 stays -0.0
            // and compares equal to zero; the bounds are T(1e-2) and T(1e4).
            let magnitude = value.abs();
            let lower: $float = 1e-2;
            let upper: $float = 1e4;
            if magnitude != 0.0 && !(lower..upper).contains(&magnitude) {
                let text = if fixed_format {
                    // std::to_chars(first, last, value, scientific, 3)
                    format!("{value:.3e}")
                } else {
                    // std::to_chars(first, last, value, scientific)
                    shortest_round_trip(value)
                };
                push_scientific(&text, target);
                return;
            }
            // std::to_chars(first, last, value, fixed, precision) into char[64].
            let text = usize::try_from(precision)
                .ok()
                .filter(|digits| *digits <= NUMERIC_BUFFER_BYTES)
                .map(|digits| format!("{value:.digits$}"))
                .filter(|text| text.len() <= NUMERIC_BUFFER_BYTES);
            match text {
                Some(text) if text.contains('.') => {
                    target.push_str(keep_one_fraction_digit(&text));
                }
                // A precision of 0 writes no point; a floating-point value
                // always gets ".0" appended.
                Some(text) => {
                    target.push_str(&text);
                    target.push_str(".0");
                }
                // std::to_chars failed (value_too_large): the source appends
                // std::to_string(static_cast<double>(value)), which is "%f".
                None => target.push_str(&format!("{:.6}", f64::from(value))),
            }
        }
    };
}

append_numeric!(
    /// `NumericFormatting::appendNumeric<double>(value, target, precision, fixed_format)`
    /// (`NumericFormatting.h:27-140`).
    ///
    /// See the module documentation for the rules. `StringUtils::toStr(double)`
    /// passes `precision = 15` and `fixed_format = false`; `toStr(double, false)`
    /// and `appendToStrLowP(double)` pass `3` and `true`.
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::concept::numeric_formatting::append_numeric_f64;
    ///
    /// let mut text = String::from("x=");
    /// append_numeric_f64(12345.6789, &mut text, 3, true);
    /// assert_eq!(text, "x=1.235e04");
    /// let mut text = String::new();
    /// append_numeric_f64(0.0625, &mut text, 3, true);
    /// assert_eq!(text, "0.062");
    /// ```
    append_numeric_f64,
    f64
);

append_numeric!(
    /// `NumericFormatting::appendNumeric<float>(value, target, precision, fixed_format)`
    /// (`NumericFormatting.h:27-140`), applied to the `float` itself.
    ///
    /// Fixed notation writes the exact `float` value's digits (`9999.999023` at
    /// precision 6) and the shortest scientific text is the shortest that reads
    /// back as the same `float` (`2.2032652e06`), not as a promoted `double`.
    /// `StringUtils::toStr(float)` passes `precision = 6` and
    /// `fixed_format = false`; the low-precision overloads pass `3` and `true`.
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::concept::numeric_formatting::append_numeric_f32;
    ///
    /// let mut text = String::new();
    /// append_numeric_f32(50.840_353, &mut text, 6, false);
    /// assert_eq!(text, "50.840351");
    /// ```
    append_numeric_f32,
    f32
);

/// Default `std::ostream` insertion of a `double` after `precision(precision)`,
/// which is C `printf("%.*g")`.
///
/// With `P` the precision (0 counts as 1) and `X` the decimal exponent after
/// rounding to `P` significant digits, the value is written in fixed notation
/// with `P - 1 - X` fraction digits when `-4 <= X < P`, otherwise in scientific
/// notation with `P - 1` fraction digits and a signed exponent of at least two
/// digits. Trailing fraction zeros and a bare point are then removed. FileInfo
/// writes statistics at precision 6 (`3.49692e+06`, `0.00209717`) and
/// consensus coordinate statistics at 15 (`0.277777777777778`);
/// `StringUtils::numberLength` uses 15 and 6. A C++ `float` is promoted to
/// `double` before it is formatted, so pass `f64::from(x)`.
///
/// Negative zero prints as `-0`, a NaN as `nan` or, with its sign bit set, as
/// `-nan` (glibc, the reference build's C library; Apple libc writes `nan`),
/// and the infinities as `inf` and `-inf`. A precision above `i32::MAX`
/// reaches `printf` as a negative `int` in both libc++ and libstdc++ and counts
/// as six.
///
/// Exact decimal ties follow the C standard, rounding half to even and then
/// stripping zeros, as glibc does. Apple libc keeps the zeros of an integer
/// below `1e15` whose tie rounds down to them (`6427405000.0` at precision 6 is
/// `6.42740e+09` there and `6.4274e+09` here); see the native differences of
/// `crate::format::file_info::text_format`.
///
/// Work is bounded: a precision above 800 is formatted as 800, which gives the
/// same text because no double has more than 767 significant digits and no
/// decimal exponent exceeds 308.
///
/// # Examples
///
/// ```
/// use openms::concept::numeric_formatting::ostream_g;
///
/// assert_eq!(ostream_g(3_496_920.99, 6), "3.49692e+06");
/// assert_eq!(ostream_g(1e-5, 6), "1e-05");
/// assert_eq!(ostream_g(2.5 / 9.0, 15), "0.277777777777778");
/// assert_eq!(ostream_g(-0.0, 6), "-0");
/// ```
pub fn ostream_g(value: f64, precision: u32) -> String {
    if !value.is_finite() {
        return printf_nonfinite(value).to_owned();
    }
    let significant = printf_precision(precision).clamp(1, MAX_EFFECTIVE_G_PRECISION);
    let scientific = format!("{value:.fraction$e}", fraction = significant - 1);
    let Some((mantissa, exponent)) = scientific.split_once('e') else {
        // `{:e}` always writes an exponent.
        return scientific;
    };
    let Ok(exponent) = exponent.parse::<i32>() else {
        return scientific;
    };
    let limit = i32::try_from(significant).unwrap_or(i32::MAX);
    if exponent < -4 || exponent >= limit {
        let sign = if exponent < 0 { '-' } else { '+' };
        format!(
            "{}e{sign}{:02}",
            strip_fraction_zeros(mantissa),
            exponent.unsigned_abs()
        )
    } else {
        // -4 <= exponent < limit <= 800, so the difference is small and non-negative.
        let fraction = usize::try_from(limit - 1 - exponent).unwrap_or(0);
        strip_fraction_zeros(&format!("{value:.fraction$}")).to_owned()
    }
}

/// The precision `printf` applies for a `UInt` digit count or `std::streamsize`
/// precision converted to `int`: a value above `i32::MAX` becomes negative and
/// counts as omitted, six.
pub(crate) fn printf_precision(requested: u32) -> usize {
    match i32::try_from(requested) {
        Ok(precision) => usize::try_from(precision).unwrap_or(OMITTED_PRECISION),
        Err(_) => OMITTED_PRECISION,
    }
}

/// `snprintf` and stream spellings of a non-finite value, as the Linux x86_64
/// Release build's glibc writes them.
///
/// Every caller reaches C `printf`: `StringUtils::number(double, UInt)`, whose
/// whole body is `std::snprintf(buf, sizeof(buf), "%.*f", (int)n, d)`
/// (`StringUtils.cpp:526-531` at core `bc9cc12`), and `std::ostringstream`,
/// whose `num_put<char>::_M_insert_float` writes through
/// `std::__convert_from_v` — `__builtin_vsnprintf` under the C locale in
/// libstdc++ 14.4.0's `x86_64-conda-linux-gnu/bits/c++locale.h:74`. glibc's
/// `__printf_fp` writes a NaN whose **sign bit is set** as `-nan`, so this
/// function does.
///
/// The `StringUtils::toStr` path does **not** come here and must not: its first
/// statement is `if (std::isnan(value)) { target += "NaN"; return; }`
/// (`NumericFormatting.h:29`), which never reads the sign bit.
///
/// Measured by `../oracle/a2-textfmt-linux` and `../oracle/a2-textfmt-nan-sweep`
/// on ibminode06 against the Release install: 912 sign-bearing `printf` rows,
/// none of which prints a sign that disagrees with the argument's sign bit.
pub(crate) fn printf_nonfinite(value: f64) -> &'static str {
    if value.is_nan() {
        if value.is_sign_negative() {
            "-nan"
        } else {
            "nan"
        }
    } else if value < 0.0 {
        "-inf"
    } else {
        "inf"
    }
}

/// `%g` without the `#` flag: trailing fraction zeros removed, then a bare point.
fn strip_fraction_zeros(text: &str) -> &str {
    if text.contains('.') {
        text.trim_end_matches('0').trim_end_matches('.')
    } else {
        text
    }
}

/// `appendNumeric`'s scientific post-processing of Rust `{:e}` text: the
/// mantissa keeps one fraction digit at least (a mantissa without a point gets
/// `.0`), then `e`, then the exponent without `+`, zero-padded to two digits.
fn push_scientific(text: &str, target: &mut String) {
    let Some((mantissa, exponent)) = text.split_once('e') else {
        target.push_str(text);
        return;
    };
    if mantissa.contains('.') {
        target.push_str(keep_one_fraction_digit(mantissa));
    } else {
        target.push_str(mantissa);
        target.push_str(".0");
    }
    target.push('e');
    match exponent.parse::<i32>() {
        Ok(exponent) if exponent < 0 => {
            target.push_str(&format!("-{:02}", exponent.unsigned_abs()));
        }
        Ok(exponent) => target.push_str(&format!("{exponent:02}")),
        Err(_) => target.push_str(exponent),
    }
}

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
    T: Copy + PartialEq + LowerExp + FromStr,
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
