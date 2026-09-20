// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! C++ stream and `StringUtils` numeric text formatting for the FileInfo report.
//!
//! `FORMAT/FileInfo.cpp` renders every number through one of four C++ text
//! paths, and the upstream reference reports compare those bytes. Each path is
//! reproduced here as a stateless function:
//!
//! | FileInfo use | Source path | Rust |
//! |---|---|---|
//! | Range bounds `number(x, 2)`; spans `number((max - min) / 60, 1)` | `StringUtils::number(double, UInt)`: `snprintf("%.*f")` into `char[64]` | [`fixed`], [`fixed_truncated`] |
//! | Total ion current, analyzer resolution, `SummaryStatistics` blocks | `std::ostream << double` in the default float field after `os.precision(n)` | [`ostream_g`] |
//! | FAIMS compensation voltages | `StringUtils::toStr(double)` through `NumericFormatting::appendNumeric` | [`to_str`], [`to_str_f32`] |
//! | `IM (FAIMS_CV): [..]` | `operator<<(std::ostream&, const std::vector<T>&)` from `ListUtilsIO.h` | [`list_to_string`], [`ToStr`] |
//!
//! Sources at core revision `bc9cc12`:
//! `src/openms/include/OpenMS/DATASTRUCTURES/StringUtils.h`,
//! `src/openms/source/DATASTRUCTURES/StringUtils.cpp`,
//! `src/common/include/OpenMS/CONCEPT/Detail/NumericFormatting.h`,
//! `src/openms/include/OpenMS/DATASTRUCTURES/ListUtilsIO.h` and
//! `src/openms/include/OpenMS/CONCEPT/Types.h`, plus the C++ standard library's
//! floating-point stream insertion that `FileInfo.cpp` relies on.
//!
//! # Stream state
//!
//! A C++ stream remembers `precision(n)`. `FileInfo.cpp` sets it once per
//! statistics block (lines 2224-2254, 2330-2362 and 2404) and it stays in force
//! for every later number written to that stream. These functions keep no
//! state: callers pass the precision explicitly and track the stream's current
//! precision themselves, starting from [`DEFAULT_STREAM_PRECISION`].
//!
//! # API mapping
//!
//! | Source member | Rust | Notes |
//! |---|---|---|
//! | `StringUtils::number(double, UInt)` | [`fixed`], [`fixed_truncated`] | strict and source-compatible |
//! | `StringUtils::toStr(double, true)`, `appendToStr(double, std::string&)` | [`to_str`] | |
//! | `StringUtils::toStr(float, true)`, `appendToStr(float, std::string&)` | [`to_str_f32`] | |
//! | `toStr`/`appendToStr` for `int`, `unsigned int`, `short`, `unsigned short`, `long`, `unsigned long`, `long long`, `unsigned long long` | [`ToStr`] for `i16`, `u16`, `i32`, `u32`, `i64`, `u64`, `isize`, `usize` | decimal |
//! | `toStr(const std::string&)`, `toStr(const char*)`, `toStr(std::string_view)` | [`ToStr`] for `str` and `String` | copied unchanged |
//! | `toStr(char)` | not ported | FileInfo streams no `char` vector |
//! | `toStr(long double, bool)`, `appendToStr(long double, std::string&)`, `appendToStrLowP(long double, std::string&)` | not ported | Rust has no `long double`; FileInfo uses none |
//! | `toStr(float, false)`, `toStr(double, false)`, `appendToStrLowP(float/double, std::string&)` | not ported | three-fraction-digit form, unused by FileInfo |
//! | `toStr(const DataValue&, bool)`, `toStr(const ParamValue&, bool)`, `appendToStr(const DataValue&, bool, std::string&)` | not ported here | the value types format themselves in `param` and `metadata` |
//! | `StringUtils::numberLength(double, UInt)` | not ported | not on the FileInfo path |
//! | `NumericFormatting::appendNumeric<T>(T, std::string&, int, bool)` | private; reached through [`to_str`] and [`to_str_f32`] with `fixed_format = false` | the `fixed_format = true` and libc++ `long double` branches are unused |
//! | `writtenDigits<float>()`, `writtenDigits<double>()` | [`WRITTEN_DIGITS_F32`], [`WRITTEN_DIGITS_F64`] | the integer and `long double` specialisations are not ported |
//! | `operator<<(std::ostream&, const std::vector<T>&)` | [`list_to_string`] | |
//! | `VecLowPrecision<T>` and its `operator<<` | not ported | built on the low-precision `toStr` |
//! | `operator<<(std::vector<std::string>&, const TString&)` | not ported | `Vec::push` |
//! | `std::ostream::operator<<(double)` and `(float)` in the default float field | [`ostream_g`] | a `float` is promoted to `double` first |
//!
//! # Preserved source conventions
//!
//! - Rounding is exact on the binary value, with exact ties to the even digit,
//!   as `snprintf` and `std::to_chars` do: `number(0.125, 2)` is `0.12` and
//!   `number(2.5, 0)` is `2`, while `number(0.005, 2)` is `0.01` because the
//!   stored double lies just above the tie.
//! - The sign of zero is kept on every path: `-0.00`, `-0` and `-0.0`.
//! - `%g` writes fixed notation when the decimal exponent `X` after rounding to
//!   `P` significant digits satisfies `-4 <= X < P`, otherwise scientific
//!   notation with a sign and at least two exponent digits (`3.49692e+06`).
//!   Trailing fraction zeros and a bare point are removed, and a precision of 0
//!   counts as 1.
//! - `toStr` writes zero and magnitudes in `[1e-2, 1e4)` with 15 (`double`) or
//!   6 (`float`) *fraction* digits, other magnitudes as the shortest scientific
//!   text that reads back as the same value (`std::to_chars` without a
//!   precision). Of the equally short texts that read back, the one nearest the
//!   exact binary value is taken, and an exact tie takes the even digit:
//!   `-135169.706298828125` prints as `-1.3516970629882812e05`. Trailing zeros
//!   are removed but one fraction digit stays, and the exponent has no `+` and at
//!   least two digits (`1.0e-05`, `3.6739e04`).
//! - Spellings: `nan` — `-nan` when the sign bit is set — plus `inf` and `-inf`
//!   from `snprintf` and streams; `NaN` whatever the sign bit, plus `inf` and
//!   `-inf`, from `toStr`.
//! - A digit count or stream precision above `INT_MAX` reaches `printf` as a
//!   negative `int`, which counts as omitted: six digits.
//!
//! # Native differences
//!
//! - [`fixed`] refuses text longer than the 63 bytes the source buffer holds.
//!   `StringUtils::number` silently cuts such text, so `number(1e100, 2)`
//!   returns a 63-digit integer. [`fixed_truncated`] reproduces the cut for
//!   tool paths that need byte parity with the C++ report.
//! - [`fixed`] also refuses a digit count above `INT_MAX` for a finite value,
//!   which the source silently prints with six digits; [`fixed_truncated`]
//!   reproduces that too.
//! - Where the `%.*f` text would be `INT_MAX` bytes or longer (only digit counts
//!   within 311 of `INT_MAX` get there), the source result depends on the C
//!   library and [`fixed_truncated`] does not reproduce it. Apple libc fails with
//!   `EOVERFLOW` and `number` returns an empty string. glibc 2.39 still returns
//!   the 63-byte digit prefix at text lengths `INT_MAX` and `2^31`, and 63 spaces
//!   for a digit count of `INT_MAX`. [`fixed_truncated`] returns the digit prefix
//!   throughout and [`fixed`] refuses. FileInfo uses digit counts 1 and 2.
//! - Default-stream (`%g`) text depends on the C library for one class of exact
//!   decimal ties, and [`ostream_g`] follows the C standard and glibc there. The
//!   class: an integer-valued double below `1e15` lying exactly halfway between
//!   its two `P`-significant-digit neighbours, where the neighbour nearer zero
//!   ends in `0`, so rounding half to even keeps it. Apple libc keeps that
//!   neighbour's trailing zeros in the scientific mantissa, while glibc and the
//!   C standard strip them: `6427405000.0` at precision 6 is `6.42740e+09` on
//!   macOS but `6.4274e+09` with glibc and here, and `205.0` at precision 2 is
//!   `2.0e+02` against `2e+02`. Outside that class both print the same: a tie
//!   that rounds up to a neighbour ending in `0` (`-95.0` at precision 1 is
//!   `-1e+02`), the same tie at or above `1e15` (`5.05e18` at precision 2 is
//!   `5e+18`), a non-integer tie, and a non-tie (`-7723301.0` at precision 6 is
//!   `-7.7233e+06`). A macOS C++ report therefore differs from the port on
//!   such values, for example a total ion current of `1463805` at precision 6,
//!   and macOS comparisons must not count them as port defects.
//! - The `printf` spelling of a NaN follows **glibc**, the reference build's C
//!   library: a NaN whose sign bit is set prints as `-nan` out of [`fixed`],
//!   [`fixed_truncated`] and [`ostream_g`]. Apple libc writes `nan` for the same
//!   bits, so a macOS C++ report differs from the port on such a value and a
//!   macOS comparison must not count that as a port defect — the same caveat the
//!   `%g` tie bullet above carries. `toStr` is not affected on either platform:
//!   `NumericFormatting.h:29` returns `NaN` before the sign bit is ever read.
//!
//!   *What is measured and what is generalised.* `../oracle/a2-textfmt-linux`
//!   re-ran the same driver, `cases.h` and pin probe, byte for byte, against the
//!   Linux x86_64 Release install on ibminode06 (conda-forge GCC 14.4.0,
//!   libstdc++ 6.0.36, glibc 2.39). Of 1018 rows exactly one differs from the
//!   macOS capture: `fff8000000000000`, in its five `printf` columns and not in
//!   its `toStr` column. The corpus holds exactly three NaN bit patterns
//!   (`7ff8000000000000`, `fff8000000000000`, and the `float` `7fc00000`) and
//!   that row is the only sign-bit NaN in it, so what is **measured** is
//!   `number(-NaN, n)` at `n` in `{0, 1, 2}` and `ostream(-NaN, p)` at `p` in
//!   `{6, 15}`. No sign-bit NaN `float` is pinned at all. Every other digit
//!   count, precision and the `float` overload are **generalised** from glibc
//!   writing the sign before `__printf_fp` dispatches on the class, which makes
//!   the spelling independent of both; a wider negative-NaN sweep is being
//!   captured separately.
//! - Only the classic `"C"` locale is modelled; FileInfo never imbues another.
//! - Work is bounded: `%g` precisions above 800 and `%f` digit counts above 1100
//!   are clamped internally. A double's exact decimal expansion has at most 767
//!   significant and 1074 fraction digits, so larger values only add zeros that
//!   `%g` strips and that [`fixed_truncated`] cuts. For `%f` that equivalence
//!   ends at the `INT_MAX`-byte limit above.
//!
//! # Evidence
//!
//! - Oracle-generated, tier 1 executed differential:
//!   `../oracle/a2-textfmt-linux/results/driver.tsv`, produced by linking the
//!   Linux x86_64 **Release** install
//!   `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576` on
//!   ibminode06 (conda-forge GCC 14.4.0, libstdc++ 6.0.36, glibc 2.39), holds
//!   386 doubles, 100 floats, 249 digit counts, 268 stream precisions and 10
//!   vectors. `tests/file_info_text_format.rs` embeds it and compares every row.
//!   It re-runs `../oracle/file-info-text-format`'s `driver.cpp`, `cases.h` and
//!   `pin_probe.cpp` byte for byte against the reference build rather than the
//!   macOS product SDK (core `4fdec46`, Debug, AppleClang 21, macOS arm64). Of
//!   the 1018 rows exactly one differs between the two captures, the sign-bit
//!   NaN of the native-difference note above; the macOS capture stands behind
//!   every other row unchanged.
//! - Executed probe, tier 2: the same oracle run compiles the pinned
//!   `NumericFormatting.h` unchanged. Its `toStr` and `number` rows are
//!   byte-identical to the SDK's, so commit `74526a8` between the oracle build
//!   and the pin leaves this text unchanged.
//! - Tie rule, tier 1 with tier 2 probes: `../oracle/text-format` links the same
//!   SDK `libOpenMS` and prints `std::to_chars` shortest scientific text next to
//!   `StringUtils::toStr` over a tie-dense sweep of 60,457 doubles and 10,973
//!   floats. The sweep holds 6,493 and 1,234 exact decimal ties, every power of
//!   two with both neighbours, and the 41 mismatches the package review found.
//!   At 46 powers of two as `double` (among them the one-sided tie at `2^-24`)
//!   and 3 as `float` (among them `2^-96`), the exact value rounded to the
//!   shortest length does not read back, so `to_chars` takes the other
//!   neighbour; no other value of the sweep does. A standard-library probe of
//!   the same sweep gives identical `to_chars` text on macOS (AppleClang 21,
//!   libc++) and Linux (GCC 13.3, libstdc++, glibc 2.39), so the tie rule does
//!   not depend on the platform. The same run records 15,195 default-stream
//!   rows, 14,950 of them exact decimal ties at their precision. Apple libc and
//!   glibc differ on all 324 rows of the class described under native
//!   differences and on none of the other 14,871; the largest differing value
//!   is `993123135600500`, and the smallest integer round-down tie at or above
//!   `1e15`, `1005000000000000`, prints alike. The run also records
//!   `snprintf("%.*f")` at the `INT_MAX`-byte limit on both platforms. The test
//!   embeds 425 selected rows, among them every sweep row where `to_chars` takes
//!   the other neighbour. The whole sweep was also run against this module on
//!   Linux with Rust 1.96.0 and 1.85.0: every `toStr` row and every
//!   default-stream glibc row matches. Before the fix 3,785 `toStr` rows
//!   differed, every one a two-sided tie, among them all 41 review rows.
//! - Retained upstream outputs: range, statistics and FAIMS lines of
//!   `FileInfo_1`, `_2`, `_3`, `_7`, `_9` and `_19` at test-data `0cb15f2`,
//!   rebuilt in the test from values read out of their inputs.
//!
//! # Class tests
//!
//! The upstream class tests at `bc9cc12`, section by section. Ported literals
//! are the `class_test_*` tests in `tests/file_info_text_format.rs`.
//!
//! | Section | Status |
//! |---|---|
//! | `String_test.cpp` `StringUtils::toStr` for `int`, `unsigned int`, `long int`, `long unsigned int`, `short int`, `short unsigned int`, `long long unsigned int` and `long long signed int` (8 sections) | ported: [`ToStr`] |
//! | `String_test.cpp` `StringUtils::toStr(float, bool full_precision)` | ported for `full_precision = true`: [`to_str_f32`]; the `false` literals are not ported |
//! | `String_test.cpp` `StringUtils::toStr(double, bool full_precision)` | ported for `full_precision = true`: [`to_str`]; the `false` literals are not ported |
//! | `String_test.cpp` `StringUtils::number(double, UInt)` | ported: [`fixed`], [`fixed_truncated`] |
//! | `String_test.cpp` `StringUtils::toStr(long double, bool full_precision)` | not ported: Rust has no `long double` |
//! | `String_test.cpp` `StringUtils::toStr(DataValue)` | not ported here: `DataValue` formats itself in `metadata` |
//! | `String_test.cpp` `StringUtils::numberLength(double, UInt)` | not ported |
//! | `String_test.cpp` `operator+` and `operator+=` of `std::string` with a number (15 sections) | not ported: string concatenation, not on the FileInfo path |
//! | `String_test.cpp`, the other 40 sections (constructors, `random`, `has*`, prefix and suffix, `substr`, `chop`, `reverse`, `trim`, quoting, `simplify`, padding, number parsing, case, `substitute`, `remove`, `ensureLastChar`, `removeWhitespaces`, splitting, `concatenate`) | outside this module |
//! | `StringUtils_test.cpp` `numberLength` and `number` | `NOT_TESTABLE` upstream ("tested in String_test.cpp") |
//! | `StringUtils_test.cpp` `[EXTRA] non-finite values round-trip through toStr/toDouble/toFloat` | formatting literals ported; the `toDouble`/`toFloat` read-back and the `DataValue` list are outside this module |
//! | `StringUtils_test.cpp`, the other 42 sections | outside this module |
//! | `ListUtilsIO_test.cpp` `[EXTRA] StringList& operator<<(StringList&, const StringType&)` | not ported: `Vec::push` |
//! | `NumericFormatting.h`, `Types.h` | no class test at the pin |
//!
//! # Follow-up
//!
//! Private copies of these rules exist in `format/mascot_generic.rs`
//! (`ostream_g`, which drops the sign of negative zero), `format/pepxml.rs`
//! (`general`, likewise), `format/mzxml.rs` (`general_format`),
//! `math/posterior_error_probability.rs` (`format_g`) and `param/value.rs`
//! (`format_float`, a `toStr` port). Consolidating them on this module is left
//! to their owners.
//!
//! `param/value.rs` `format_float` and `format_float32` take the scientific
//! `toStr` text straight from `{:e}`, so they still break exact ties upwards:
//! the defect [`to_str`] and [`to_str_f32`] fix. Executed through
//! `data_structures::list::to_string_list`, they print `1.5000000000000003e15`
//! and `-1.3516970629882813e05` for `1500000000000000.25` and
//! `-135169.706298828125`, and `2.2032653e06` and `2.5000003e06` for the floats
//! `2203265.25` and `2500000.25`. `StringUtils::toStr` and this module print
//! `1.5000000000000002e15`, `-1.3516970629882812e05`, `2.2032652e06` and
//! `2.5000002e06`. Every caller inherits the difference: `ParamValue` text,
//! `data_structures::list` text, ProForma conversion and resolution, and the
//! controlled-vocabulary, imzML, Mascot generic, MSstats, mzTab, mzTab-M,
//! pepXML, percolator and transformation-XML writers. Moving the two functions
//! onto this module fixes it; `param/value.rs` belongs to another package and is
//! unchanged here.
//!
//! [`fixed`]: crate::format::file_info::text_format::fixed
//! [`fixed_truncated`]: crate::format::file_info::text_format::fixed_truncated
//! [`ostream_g`]: crate::format::file_info::text_format::ostream_g
//! [`to_str`]: crate::format::file_info::text_format::to_str
//! [`to_str_f32`]: crate::format::file_info::text_format::to_str_f32
//! [`list_to_string`]: crate::format::file_info::text_format::list_to_string
//! [`ToStr`]: crate::format::file_info::text_format::ToStr
//! [`DEFAULT_STREAM_PRECISION`]: crate::format::file_info::text_format::DEFAULT_STREAM_PRECISION
//! [`WRITTEN_DIGITS_F32`]: crate::format::file_info::text_format::WRITTEN_DIGITS_F32
//! [`WRITTEN_DIGITS_F64`]: crate::format::file_info::text_format::WRITTEN_DIGITS_F64

use crate::error::{Error, Result};
use std::fmt::LowerExp;
use std::str::FromStr;

/// Largest text, in bytes, that `StringUtils::number` can return.
///
/// The source formats into `char buf[64]` with `snprintf`, which stores at most
/// 63 bytes before the terminating NUL and silently drops the rest.
pub const NUMBER_MAX_TEXT_BYTES: usize = 63;

/// Precision of a newly constructed C++ stream, `std::ios_base::precision()`.
///
/// FileInfo writes the total ion current and analyzer resolutions at this
/// precision, before any statistics block changes it.
pub const DEFAULT_STREAM_PRECISION: u32 = 6;

/// `writtenDigits<float>()`, which is `std::numeric_limits<float>::digits10`.
///
/// FileInfo sets this stream precision for feature and peak intensities and for
/// qualities, whose C++ types (`Feature::IntensityType`,
/// `Feature::QualityType`, `Peak1D::IntensityType`) are `float`. `toStr(float)`
/// uses it as its fraction-digit count.
pub const WRITTEN_DIGITS_F32: u32 = 6;

/// `writtenDigits<double>()`, which is `std::numeric_limits<double>::digits10`.
///
/// FileInfo sets this stream precision for consensus coordinate differences
/// (`ConsensusFeature::CoordinateType` is `double`). `toStr(double)` uses it as
/// its fraction-digit count.
pub const WRITTEN_DIGITS_F64: u32 = 15;

/// Output ceiling of [`list_to_string`] in bytes, checked before allocating.
///
/// The source has no ceiling. FileInfo streams a handful of compensation
/// voltages, far below it.
pub const MAX_LIST_TEXT_BYTES: usize = 64 * 1024 * 1024;

/// `%g` precisions at or above this give the same text: no double has more than
/// 767 significant digits, and no exponent reaches it.
const MAX_EFFECTIVE_G_PRECISION: usize = 800;
/// `%f` digit counts at or above this only append zeros: no double has more than
/// 1074 fraction digits.
const MAX_EFFECTIVE_FIXED_DIGITS: usize = 1100;
/// Precision `printf` applies when the requested one is negative.
const OMITTED_PRECISION: usize = 6;
/// Longest [`to_str`] text: sign, 17 significant digits, point and `e-308`.
const F64_TO_STR_MAX_BYTES: usize = 24;
/// Longest [`to_str_f32`] text: sign, 9 significant digits, point and `e-45`.
const F32_TO_STR_MAX_BYTES: usize = 15;
/// Longest decimal 64-bit integer, `-9223372036854775808`.
const INTEGER_TO_STR_MAX_BYTES: usize = 20;

/// `StringUtils::number(d, n)`: `value` with exactly `digits` fraction digits,
/// as C `snprintf("%.*f")` writes it.
///
/// Rounding uses the exact binary value with ties to even, and the sign of zero
/// is kept (`-0.00`). NaN and the infinities ignore `digits` and print as
/// `printf` writes them: `nan`, `-nan` for a NaN whose sign bit is set, `inf`
/// and `-inf` (see [`nonfinite`]). The source parameter is a `double`, so a
/// `float` argument is promoted first; pass `f64::from(x)`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] where the source would not print the full
/// number: the text is longer than [`NUMBER_MAX_TEXT_BYTES`], which `number`
/// silently cuts, or `digits` exceeds `i32::MAX` for a finite value, which the
/// source converts to a negative `int` and so prints with six digits. Use
/// [`fixed_truncated`] where byte parity with those outputs is required.
///
/// # Examples
///
/// ```
/// use openms::format::file_info::text_format::fixed;
///
/// assert_eq!(fixed(-1.0, 2).ok().as_deref(), Some("-1.00"));
/// assert_eq!(fixed(0.125, 2).ok().as_deref(), Some("0.12"));
/// assert!(fixed(1e100, 2).is_err());
/// ```
pub fn fixed(value: f64, digits: u32) -> Result<String> {
    if !value.is_finite() {
        return Ok(nonfinite(value).to_owned());
    }
    if i32::try_from(digits).is_err() {
        return Err(Error::InvalidValue(format!(
            "StringUtils::number digit count {digits} exceeds INT_MAX; the source prints six digits instead"
        )));
    }
    // "0." and 62 digits already overflow the buffer, so refuse before formatting.
    let digit_count = match usize::try_from(digits) {
        Ok(count) if count <= NUMBER_MAX_TEXT_BYTES - 2 => count,
        _ => return Err(truncated(value, digits)),
    };
    let text = format!("{value:.digit_count$}");
    if text.len() > NUMBER_MAX_TEXT_BYTES {
        return Err(truncated(value, digits));
    }
    Ok(text)
}

/// `StringUtils::number(d, n)` including its silent truncation, for byte parity.
///
/// Identical to [`fixed`] wherever that succeeds. Otherwise it returns what the
/// source returns: the first [`NUMBER_MAX_TEXT_BYTES`] bytes of the `%.*f` text,
/// and six fraction digits when `digits` exceeds `i32::MAX`. The source gives no
/// diagnostic, and a cut result shows a wrong magnitude (`number(1e100, 2)` is a
/// 63-digit integer), so prefer [`fixed`] outside tool paths that must
/// reproduce the C++ report byte for byte.
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
/// use openms::format::file_info::text_format::fixed_truncated;
///
/// assert_eq!(fixed_truncated(6.0 / 60.0, 1), "0.1");
/// assert_eq!(fixed_truncated(1e100, 2).len(), 63);
/// assert_eq!(fixed_truncated(0.125, u32::MAX), "0.125000");
/// ```
pub fn fixed_truncated(value: f64, digits: u32) -> String {
    if !value.is_finite() {
        return nonfinite(value).to_owned();
    }
    let digit_count = source_precision(digits).min(MAX_EFFECTIVE_FIXED_DIGITS);
    let text = format!("{value:.digit_count$}");
    match text.get(..NUMBER_MAX_TEXT_BYTES) {
        Some(prefix) => prefix.to_owned(),
        None => text,
    }
}

/// Default `std::ostream` insertion of a `double` after `precision(precision)`,
/// which is C `printf("%.*g")`.
///
/// With `P` the precision (0 counts as 1) and `X` the decimal exponent after
/// rounding to `P` significant digits, the value is written in fixed notation
/// with `P - 1 - X` fraction digits when `-4 <= X < P`, otherwise in scientific
/// notation with `P - 1` fraction digits and a signed exponent of at least two
/// digits. Trailing fraction zeros and a bare point are then removed. FileInfo
/// writes statistics at [`WRITTEN_DIGITS_F32`] (`3.49692e+06`, `0.00209717`)
/// and consensus coordinate statistics at [`WRITTEN_DIGITS_F64`]
/// (`0.277777777777778`). A C++ `float` is promoted to `double` before it is
/// formatted, so pass `f64::from(x)`.
///
/// Negative zero prints as `-0`, a NaN as `nan` or, with its sign bit set, as
/// `-nan` (see [`nonfinite`]), and the infinities as `inf` and `-inf`. A
/// precision above `i32::MAX` reaches `printf` as a negative `int` in both
/// libc++ and libstdc++ and counts as six.
///
/// Exact decimal ties follow the C standard, rounding half to even and then
/// stripping zeros, as glibc does. Apple libc keeps the zeros of an integer below
/// `1e15` whose tie rounds down to them (`6427405000.0` at precision 6 is
/// `6.42740e+09` there and `6.4274e+09` here); see the module's native
/// differences.
///
/// Work is bounded: a precision above 800 is formatted as 800, which gives the
/// same text because no double has more than 767 significant digits and no
/// decimal exponent exceeds 308.
///
/// # Examples
///
/// ```
/// use openms::format::file_info::text_format::ostream_g;
///
/// assert_eq!(ostream_g(3_496_920.99, 6), "3.49692e+06");
/// assert_eq!(ostream_g(1e-5, 6), "1e-05");
/// assert_eq!(ostream_g(2.5 / 9.0, 15), "0.277777777777778");
/// assert_eq!(ostream_g(-0.0, 6), "-0");
/// ```
pub fn ostream_g(value: f64, precision: u32) -> String {
    if !value.is_finite() {
        return nonfinite(value).to_owned();
    }
    let significant = source_precision(precision).clamp(1, MAX_EFFECTIVE_G_PRECISION);
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

/// `StringUtils::toStr(double)`: the full-precision text of
/// `NumericFormatting::appendNumeric(d, s, writtenDigits<double>(), false)`.
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
/// use openms::format::file_info::text_format::to_str;
///
/// assert_eq!(to_str(-65.0), "-65.0");
/// assert_eq!(to_str(1e-5), "1.0e-05");
/// assert_eq!(to_str(36739.0), "3.6739e04");
/// assert_eq!(to_str(-135_169.706_298_828_125), "-1.3516970629882812e05");
/// ```
pub fn to_str(value: f64) -> String {
    let mut text = String::new();
    append_numeric_f64(value, &mut text);
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
/// use openms::format::file_info::text_format::to_str_f32;
///
/// assert_eq!(to_str_f32(50.840_353), "50.840351");
/// assert_eq!(to_str_f32(1e-5), "1.0e-05");
/// assert_eq!(to_str_f32(2_203_265.25), "2.2032652e06");
/// ```
pub fn to_str_f32(value: f32) -> String {
    let mut text = String::new();
    append_numeric_f32(value, &mut text);
    text
}

/// Element text for [`list_to_string`], mirroring the `StringUtils::toStr`
/// overload set that the `ListUtilsIO.h` stream operator calls.
///
/// Implemented for `f64` ([`to_str`]), `f32` ([`to_str_f32`]), the signed and
/// unsigned 16-, 32- and 64-bit and pointer-sized integers (decimal), `str` and
/// `String` (copied unchanged), and references to any of these. `bool` is not
/// implemented: C++ would promote it to `int` and print `1` or `0`.
pub trait ToStr {
    /// An upper bound of the bytes [`ToStr::append_to_str`] appends.
    ///
    /// [`list_to_string`] sums these bounds to enforce
    /// [`MAX_LIST_TEXT_BYTES`] before it allocates, so an implementation must
    /// not return less than it appends.
    fn text_len_bound(&self) -> usize;

    /// Appends the `StringUtils::toStr` text of `self` to `target`.
    fn append_to_str(&self, target: &mut String);
}

impl ToStr for f64 {
    fn text_len_bound(&self) -> usize {
        F64_TO_STR_MAX_BYTES
    }
    fn append_to_str(&self, target: &mut String) {
        append_numeric_f64(*self, target);
    }
}

impl ToStr for f32 {
    fn text_len_bound(&self) -> usize {
        F32_TO_STR_MAX_BYTES
    }
    fn append_to_str(&self, target: &mut String) {
        append_numeric_f32(*self, target);
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

/// `operator<<(std::ostream&, const std::vector<T>&)` from `ListUtilsIO.h`:
/// `[`, each element's `StringUtils::toStr` text separated by `, `, then `]`.
///
/// Elements are neither quoted nor escaped, so a string containing `, ` cannot
/// be told apart from two elements (`["a, b", ""]` prints `[a, b, ]`), exactly
/// as in the source. FileInfo converts its FAIMS compensation voltages with
/// [`to_str`] into a string vector and streams that, giving `[-65.0]`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the text could exceed
/// [`MAX_LIST_TEXT_BYTES`], judged from [`ToStr::text_len_bound`] before any
/// allocation, or when the allocation fails. The source has no ceiling.
///
/// # Examples
///
/// ```
/// use openms::format::file_info::text_format::{list_to_string, to_str};
///
/// let cvs = [to_str(-65.0)];
/// assert_eq!(list_to_string(&cvs).ok().as_deref(), Some("[-65.0]"));
/// assert_eq!(list_to_string::<f64>(&[]).ok().as_deref(), Some("[]"));
/// ```
pub fn list_to_string<T: ToStr>(items: &[T]) -> Result<String> {
    let mut bound = 2usize;
    for (index, item) in items.iter().enumerate() {
        let separator = if index == 0 { 0 } else { 2 };
        bound = bound
            .checked_add(separator)
            .and_then(|total| total.checked_add(item.text_len_bound()))
            .filter(|total| *total <= MAX_LIST_TEXT_BYTES)
            .ok_or_else(|| {
                Error::InvalidValue(format!(
                    "text of a {}-element vector may exceed {MAX_LIST_TEXT_BYTES} bytes",
                    items.len()
                ))
            })?;
    }
    let mut text = String::new();
    text.try_reserve_exact(bound).map_err(|_| {
        Error::InvalidValue(format!("cannot allocate {bound} bytes of vector text"))
    })?;
    text.push('[');
    for (index, item) in items.iter().enumerate() {
        if index > 0 {
            text.push_str(", ");
        }
        item.append_to_str(&mut text);
    }
    text.push(']');
    Ok(text)
}

/// The precision `printf` applies for a `UInt` digit count or `std::streamsize`
/// precision converted to `int`.
fn source_precision(requested: u32) -> usize {
    match i32::try_from(requested) {
        Ok(precision) => usize::try_from(precision).unwrap_or(OMITTED_PRECISION),
        Err(_) => OMITTED_PRECISION,
    }
}

/// `snprintf` and stream spellings of a non-finite value, as the Linux x86_64
/// Release build's glibc writes them.
///
/// Every caller of this function reaches C `printf`: [`fixed`] and
/// [`fixed_truncated`] are `StringUtils::number(double, UInt)`, whose whole body
/// is `std::snprintf(buf, sizeof(buf), "%.*f", (int)n, d)`
/// (`StringUtils.cpp:526-531` at core `bc9cc12`), and [`ostream_g`] is
/// `std::ostringstream`, whose `num_put<char>::_M_insert_float` writes through
/// `std::__convert_from_v` — `__builtin_vsnprintf` under the C locale in
/// libstdc++ 14.4.0's `x86_64-conda-linux-gnu/bits/c++locale.h:74`. glibc's
/// `__printf_fp` writes a NaN whose **sign bit is set** as `-nan`, so this
/// function does.
///
/// The `StringUtils::toStr` path does **not** come here and must not: its first
/// statement is `if (std::isnan(value)) { target += "NaN"; return; }`
/// (`NumericFormatting.h:29`), which never reads the sign bit. That is
/// [`append_numeric_f64`] and [`append_numeric_f32`], which write `NaN` for
/// either sign.
///
/// Measured by `../oracle/a2-textfmt-linux` on ibminode06 against the Release
/// install: of its 1018 rows exactly one differs from the macOS capture, the
/// `fff8000000000000` row, and it differs in the five `printf` columns and not
/// in the `toStr` column. That row pins `number(-NaN, n)` for `n` in `{0, 1, 2}`
/// and `ostream(-NaN, p)` for `p` in `{6, 15}`; no other sign-bit NaN, and no
/// sign-bit NaN `float`, is in the corpus. Every other precision and the `float`
/// overload follow from glibc writing the sign before it dispatches on the
/// class, which the module's native-difference note states as the
/// generalisation it is.
fn nonfinite(value: f64) -> &'static str {
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

fn truncated(value: f64, digits: u32) -> Error {
    Error::InvalidValue(format!(
        "StringUtils::number({value:e}, {digits}) needs more than {NUMBER_MAX_TEXT_BYTES} bytes, \
         which the source silently truncates"
    ))
}

/// `%g` without the `#` flag: trailing fraction zeros removed, then a bare point.
fn strip_fraction_zeros(text: &str) -> &str {
    if text.contains('.') {
        text.trim_end_matches('0').trim_end_matches('.')
    } else {
        text
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

macro_rules! append_numeric {
    ($name:ident, $float:ty, $fraction_digits:expr) => {
        /// `NumericFormatting::appendNumeric(value, target, writtenDigits, false)`
        /// (`NumericFormatting.h:27-140`).
        fn $name(value: $float, target: &mut String) {
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
                // std::to_chars(scientific) without a precision.
                let shortest = shortest_round_trip(value);
                let Some((mantissa, exponent)) = shortest.split_once('e') else {
                    target.push_str(&shortest);
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
            } else {
                let text = format!("{value:.digits$}", digits = $fraction_digits);
                target.push_str(keep_one_fraction_digit(&text));
            }
        }
    };
}
append_numeric!(append_numeric_f64, f64, 15);
append_numeric!(append_numeric_f32, f32, 6);
