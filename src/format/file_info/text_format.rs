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
//! Since B8a the functions themselves live in
//! `crate::data_structures::string_utils` (the `StringUtils.h` port) and
//! `crate::concept::numeric_formatting` (`NumericFormatting.h` and the C stream
//! conversion); this module re-exports them under the names FileInfo uses and
//! keeps the `ListUtilsIO.h` stream operator and the FileInfo evidence.
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
//! | `toStr(char)` | [`ToStr`] for `char` | one character; a `bool` element is promoted to `int` (`1`/`0`) as in C++ |
//! | `toStr(long double, bool)`, `appendToStr(long double, std::string&)`, `appendToStrLowP(long double, std::string&)` | not ported | Rust has no `long double`; FileInfo uses none |
//! | `toStr(double, false)`, `appendToStrLowP(double, std::string&)` | [`to_str_low_precision`] (`string_utils::to_str_low_p`) | three fraction digits; not on the FileInfo path |
//! | `toStr(float, false)`, `appendToStrLowP(float, std::string&)` | [`to_str_low_precision_f32`] (`string_utils::to_str_f32_low_p`) | |
//! | `toStr(const DataValue&, bool)`, `toStr(const ParamValue&, bool)`, `appendToStr(const DataValue&, bool, std::string&)` | `MetaValue::to_text`, `ParamValue::to_text` | the value types format themselves in `metadata` and `param` |
//! | `StringUtils::numberLength(double, UInt)` | `crate::data_structures::string_utils::number_length` | not on the FileInfo path |
//! | `NumericFormatting::appendNumeric<T>(T, std::string&, int, bool)` | `crate::concept::numeric_formatting::append_numeric_f64` and `append_numeric_f32`, reached here through [`to_str`] and [`to_str_f32`] | the libc++ `long double` branch is not portable |
//! | `writtenDigits<float>()`, `writtenDigits<double>()` | [`WRITTEN_DIGITS_F32`], [`WRITTEN_DIGITS_F64`] | the integer and `long double` specialisations are not ported |
//! | `operator<<(std::ostream&, const std::vector<T>&)` | [`list_to_string`] | |
//! | `VecLowPrecision<T>`, its member `value` and constructor | [`VecLowPrecision`], `value`, [`VecLowPrecision::new`] | also `From<&[T]>` and `From<&Vec<T>>` |
//! | `operator<<(std::ostream&, const VecLowPrecision<T>&)` | `Display` for [`VecLowPrecision`]; [`list_to_string_low_precision`] and [`VecLowPrecision::to_text`] into a bounded string | `f32` and `f64` elements, [`ToStrLowPrecision`] |
//! | `operator<<(std::vector<std::string>&, const TString&)` | [`StringListPush::push_item`] | chains like the source operator; `<<` cannot be overloaded on `Vec<String>` |
//! | `PrecisionWrapper<T>`, `precisionWrapper(rhs)` and its `operator<<` (`CONCEPT/PrecisionWrapper.h`) | `concept::precision_wrapper` | writes the [`to_str`] and [`to_str_f32`] text |
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
//!   *What is measured.* `../oracle/a2-textfmt-linux` re-ran the same driver,
//!   `cases.h` and pin probe, byte for byte, against the Linux x86_64 Release
//!   install on ibminode06 (conda-forge GCC 14.4.0, libstdc++ 6.0.36, glibc
//!   2.39). Of 1018 rows exactly one differs from the macOS capture:
//!   `fff8000000000000`, in its five `printf` columns and not in its `toStr`
//!   column. That corpus holds only three NaN bit patterns and that row is its
//!   only sign-bit NaN, so on its own it would pin `number(-NaN, n)` at `n` in
//!   `{0, 1, 2}` and `ostream(-NaN, p)` at `p` in `{6, 15}`, and no sign-bit
//!   NaN `float` at all. The companion sweep `../oracle/a2-textfmt-nan-sweep`
//!   closes that gap against the same install: six NaN shapes (quiet, smallest
//!   signalling, signalling with the payload's top bit set, an arbitrary quiet
//!   payload, the largest payload, and canonical-payload-plus-one) by both
//!   signs by `double` and `float`, over 22 digit counts and 16 precisions.
//!   **912 sign-bearing `printf` rows, of which 0 print a sign that disagrees
//!   with the argument's sign bit**, and 24 `toStr` rows that all print `NaN`.
//!   So both halves of the rule are measured rather than generalised: glibc
//!   keeps the sign everywhere it was swept, and `appendNumeric` discards it
//!   everywhere. One thing the sweep also shows, which the port does not depend
//!   on: `number` takes a `double`, so a `float` argument is promoted first, and
//!   the promotion keeps the sign but quiets a signalling NaN and shifts the
//!   payload.
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
//! [`to_str_low_precision`]: crate::format::file_info::text_format::to_str_low_precision
//! [`to_str_low_precision_f32`]: crate::format::file_info::text_format::to_str_low_precision_f32
//! [`list_to_string_low_precision`]: crate::format::file_info::text_format::list_to_string_low_precision
//! [`ToStrLowPrecision`]: crate::format::file_info::text_format::ToStrLowPrecision
//! [`VecLowPrecision`]: crate::format::file_info::text_format::VecLowPrecision
//! [`VecLowPrecision::new`]: crate::format::file_info::text_format::VecLowPrecision::new
//! [`VecLowPrecision::to_text`]: crate::format::file_info::text_format::VecLowPrecision::to_text
//! [`StringListPush`]: crate::format::file_info::text_format::StringListPush
//! [`StringListPush::push_item`]: crate::format::file_info::text_format::StringListPush::push_item
//! [`DEFAULT_STREAM_PRECISION`]: crate::format::file_info::text_format::DEFAULT_STREAM_PRECISION
//! [`WRITTEN_DIGITS_F32`]: crate::format::file_info::text_format::WRITTEN_DIGITS_F32
//! [`WRITTEN_DIGITS_F64`]: crate::format::file_info::text_format::WRITTEN_DIGITS_F64

use crate::error::{Error, Result};

pub use crate::concept::numeric_formatting::ostream_g;
/// `StringUtils::number(double, UInt)`, strict, under its FileInfo name: see
/// [`crate::data_structures::string_utils::number`].
pub use crate::data_structures::string_utils::number as fixed;
/// `StringUtils::number(double, UInt)` with its silent truncation, under its
/// FileInfo name: see [`crate::data_structures::string_utils::number_truncated`].
pub use crate::data_structures::string_utils::number_truncated as fixed_truncated;
pub use crate::data_structures::string_utils::{
    NUMBER_MAX_TEXT_BYTES, ToStr, WRITTEN_DIGITS_F32, WRITTEN_DIGITS_F64, to_str, to_str_f32,
};
/// `StringUtils::toStr(double, false)` under its FileInfo-era name: see
/// [`crate::data_structures::string_utils::to_str_low_p`].
pub use crate::data_structures::string_utils::to_str_low_p as to_str_low_precision;
/// `StringUtils::toStr(float, false)` under its FileInfo-era name: see
/// [`crate::data_structures::string_utils::to_str_f32_low_p`].
pub use crate::data_structures::string_utils::to_str_f32_low_p as to_str_low_precision_f32;
use crate::data_structures::string_utils::{append_to_str_low_p, append_to_str_low_p_f32};

/// Longest three-fraction-digit `toStr` text: `-9.999e-308` (11 bytes) or
/// `-9999.999` (9 bytes), with one byte to spare.
const LOW_PRECISION_TO_STR_MAX_BYTES: usize = 12;

/// Precision of a newly constructed C++ stream, `std::ios_base::precision()`.
///
/// FileInfo writes the total ion current and analyzer resolutions at this
/// precision, before any statistics block changes it.
pub const DEFAULT_STREAM_PRECISION: u32 = 6;

/// Output ceiling of [`list_to_string`] in bytes, checked before allocating.
///
/// The source has no ceiling. FileInfo streams a handful of compensation
/// voltages, far below it.
pub const MAX_LIST_TEXT_BYTES: usize = 64 * 1024 * 1024;

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

/// A floating-point value with a three-fraction-digit `StringUtils::toStr`
/// form: the element types [`VecLowPrecision`] and
/// [`list_to_string_low_precision`] accept.
///
/// The source's `VecLowPrecision` operator "works for vectors of float, double,
/// long double"; it is implemented here for `f32` ([`to_str_low_precision_f32`])
/// and `f64` ([`to_str_low_precision`]) and references to them. Rust has no
/// `long double`.
pub trait ToStrLowPrecision {
    /// An upper bound of the bytes
    /// [`ToStrLowPrecision::append_to_str_low_precision`] appends, used as
    /// [`ToStr::text_len_bound`] is.
    fn text_len_bound_low_precision(&self) -> usize;

    /// Appends the `StringUtils::toStr(value, false)` text of `self` to
    /// `target`.
    fn append_to_str_low_precision(&self, target: &mut String);
}

impl ToStrLowPrecision for f64 {
    fn text_len_bound_low_precision(&self) -> usize {
        LOW_PRECISION_TO_STR_MAX_BYTES
    }
    fn append_to_str_low_precision(&self, target: &mut String) {
        append_to_str_low_p(*self, target);
    }
}

impl ToStrLowPrecision for f32 {
    fn text_len_bound_low_precision(&self) -> usize {
        LOW_PRECISION_TO_STR_MAX_BYTES
    }
    fn append_to_str_low_precision(&self, target: &mut String) {
        append_to_str_low_p_f32(*self, target);
    }
}

impl<T: ToStrLowPrecision + ?Sized> ToStrLowPrecision for &T {
    fn text_len_bound_low_precision(&self) -> usize {
        (**self).text_len_bound_low_precision()
    }
    fn append_to_str_low_precision(&self, target: &mut String) {
        (**self).append_to_str_low_precision(target);
    }
}

/// `operator<<(std::ostream&, const VecLowPrecision<T>&)` into a string: `[`,
/// each element's `StringUtils::toStr(value, false)` text separated by `, `,
/// then `]`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the text could exceed
/// [`MAX_LIST_TEXT_BYTES`], judged before any allocation, or when the
/// allocation fails, as [`list_to_string`] does. The source has no ceiling.
///
/// # Examples
///
/// ```
/// use openms::format::file_info::text_format::list_to_string_low_precision;
///
/// let text = list_to_string_low_precision(&[17.012345, 2.0]);
/// assert_eq!(text.ok().as_deref(), Some("[17.012, 2.0]"));
/// ```
pub fn list_to_string_low_precision<T: ToStrLowPrecision>(items: &[T]) -> Result<String> {
    let mut bound = 2usize;
    for (index, item) in items.iter().enumerate() {
        let separator = if index == 0 { 0 } else { 2 };
        bound = bound
            .checked_add(separator)
            .and_then(|total| total.checked_add(item.text_len_bound_low_precision()))
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
        item.append_to_str_low_precision(&mut text);
    }
    text.push(']');
    Ok(text)
}

/// A borrowed vector written with three fraction digits: source
/// `VecLowPrecision<T>`.
///
/// The source struct holds `const std::vector<T>& value` and exists only so
/// that `os << VecLowPrecision(my_vec)` selects the low-precision stream
/// operator. Here [`std::fmt::Display`] is that operator: it streams the same
/// text as [`list_to_string_low_precision`] straight into the formatter, so,
/// like the source's `ostream` insertion, it builds no whole-list string and
/// needs no ceiling. [`VecLowPrecision::to_text`] is the bounded, allocating
/// form.
///
/// # Examples
///
/// ```
/// use openms::format::file_info::text_format::VecLowPrecision;
///
/// let values = vec![17.012345f64, 2.0, 1.0e6];
/// assert_eq!(VecLowPrecision::new(&values).to_string(), "[17.012, 2.0, 1.0e06]");
/// assert_eq!(VecLowPrecision::<f32>::new(&[]).to_string(), "[]");
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VecLowPrecision<'a, T> {
    /// The borrowed elements: source member `value`.
    pub value: &'a [T],
}

impl<'a, T> VecLowPrecision<'a, T> {
    /// Wraps `value`: source constructor `VecLowPrecision(const std::vector<T>&)`.
    pub const fn new(value: &'a [T]) -> Self {
        Self { value }
    }
}

impl<T: ToStrLowPrecision> VecLowPrecision<'_, T> {
    /// The text [`std::fmt::Display`] writes, as an owned string.
    ///
    /// # Errors
    ///
    /// As [`list_to_string_low_precision`].
    pub fn to_text(&self) -> Result<String> {
        list_to_string_low_precision(self.value)
    }
}

impl<'a, T> From<&'a [T]> for VecLowPrecision<'a, T> {
    fn from(value: &'a [T]) -> Self {
        Self::new(value)
    }
}

impl<'a, T> From<&'a Vec<T>> for VecLowPrecision<'a, T> {
    fn from(value: &'a Vec<T>) -> Self {
        Self::new(value)
    }
}

impl<T: ToStrLowPrecision> std::fmt::Display for VecLowPrecision<'_, T> {
    /// Writes `[`, the elements' three-fraction-digit text separated by `, `,
    /// and `]`, as the source operator does.
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str("[")?;
        let mut element = String::with_capacity(LOW_PRECISION_TO_STR_MAX_BYTES);
        for (index, item) in self.value.iter().enumerate() {
            if index > 0 {
                formatter.write_str(", ")?;
            }
            element.clear();
            item.append_to_str_low_precision(&mut element);
            formatter.write_str(&element)?;
        }
        formatter.write_str("]")
    }
}

/// Appends entries to a string list "with less code": source
/// `operator<<(std::vector<std::string>&, const TString&)` from
/// `ListUtilsIO.h`.
///
/// The source pushes `string` and returns the list so that insertions chain,
/// `list << "a" << "b"`. Rust cannot overload `<<` on `Vec<String>` (both the
/// trait and the type are foreign), so the operator is the chaining method
/// [`StringListPush::push_item`]; plain [`Vec::push`] is the unchained form.
///
/// # Examples
///
/// ```
/// use openms::format::file_info::text_format::StringListPush;
///
/// let mut list: Vec<String> = Vec::new();
/// list.push_item("a").push_item("b").push_item(String::from("c"));
/// assert_eq!(list, ["a", "b", "c"]);
/// ```
pub trait StringListPush {
    /// Appends `item` and returns the list, for chaining.
    fn push_item<S: Into<String>>(&mut self, item: S) -> &mut Self;
}

impl StringListPush for Vec<String> {
    fn push_item<S: Into<String>>(&mut self, item: S) -> &mut Self {
        self.push(item.into());
        self
    }
}
