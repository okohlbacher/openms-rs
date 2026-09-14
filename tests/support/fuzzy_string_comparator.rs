// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Fuzzy comparison of text that tolerates numeric differences: a port of
//! `OpenMS::FuzzyStringComparator` (`src/testframework`, core bc9cc12) and of the
//! `FuzzyDiff` TOPP tool contract (topp 174b576), as shared **test support**.
//!
//! This is not library code. Integration tests include it with
//!
//! ```text
//! #[path = "support/fuzzy_string_comparator.rs"]
//! mod fuzzy;
//! ```
//!
//! It depends on `std` only. The upstream test suite compares tool output with
//! `FuzzyDiff -test -ini FuzzyDiff.ini [-whitelist ...]`; a Rust test reproduces a
//! registered comparison with [`FuzzyDiffSettings::load_ini`] on the pinned
//! `tests/data/fuzzy_string_comparator/FuzzyDiff.ini`, the registration's
//! whitelist ([`FuzzyDiffSettings::with_whitelist`]) and [`fuzzy_diff`] or
//! [`FuzzyDiffSettings::compare_bytes`]. Tolerances are never widened here; a test
//! that needs other values sets them explicitly and names the reason.
//!
//! The mapping to the C++ members, the preserved quirks and the native
//! differences are listed in `docs/FUZZY_STRING_COMPARATOR_SUPPORT.md`.
#![allow(dead_code)]

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, ErrorKind, Read};
use std::path::Path;

/// Largest input a comparison reads from one stream or file, in bytes.
///
/// Native bound; the source reads without limit. Exceeding it fails the
/// comparison with a log message instead of producing a verdict on partial data.
pub const MAX_INPUT_BYTES: u64 = 1 << 30;

/// Largest log kept in a [`LogDestination::Buffer`], in bytes. Later output is
/// dropped and [`FuzzyStringComparator::log_truncated`] turns true. Native bound.
pub const MAX_LOG_BYTES: usize = 256 << 20;

/// Failure reports written per comparison before further reports are suppressed.
///
/// Only verbose level 3 and above continues after a failure, so only there can
/// the limit be reached; the verdict is unaffected. Native bound.
pub const MAX_FAILURE_REPORTS: usize = 10_000;

/// Largest ParamXML (INI) file [`FuzzyDiffSettings::from_ini`] accepts, in bytes.
pub const MAX_INI_BYTES: usize = 16 << 20;

/// Failure messages of `compareLines_`, verbatim.
pub mod message {
    /// One number is NaN and the other is not.
    pub const ONE_NAN: &str = "one value is NaN and the other is not";
    /// Two infinities of opposite sign.
    pub const INFINITY_SIGNS: &str = "infinities have different signs";
    /// One number is infinite and the other is finite.
    pub const ONE_INFINITE: &str = "one value is infinite and the other is not";
    /// The first number is zero and the second is not, beyond the absolute tolerance.
    pub const FIRST_ZERO: &str = "element_1_.number_ is zero, but element_2_.number_ is not";
    /// The second number is zero and the first is not, beyond the absolute tolerance.
    pub const SECOND_ZERO: &str = "element_1_.number_ is not zero, but element_2_.number_ is";
    /// Opposite signs beyond the absolute tolerance.
    pub const SIGNS: &str = "numbers have different signs";
    /// Ratio beyond both tolerances.
    pub const RATIO: &str = "ratio of numbers is too large";
    /// A number on the first side meets a non-number on the second.
    pub const NUMBER_FIRST: &str = "input_1 is a number, but input_2 is not";
    /// A non-number on the first side meets a number on the second.
    pub const NUMBER_SECOND: &str = "input_1 is not a number, but input_2 is";
    /// Whitespace on the first side meets a letter on the second.
    pub const SPACE_FIRST: &str = "input_1 is whitespace, but input_2 is not";
    /// A letter on the first side meets whitespace on the second.
    pub const SPACE_SECOND: &str = "input_1 is not whitespace, but input_2 is";
    /// Two different letters.
    pub const LETTERS: &str = "different letters";
    /// The second line ended first.
    pub const SECOND_SHORTER: &str = "line from input_2 is shorter than line from input_1";
    /// The first line ended first.
    pub const FIRST_SHORTER: &str = "line from input_1 is shorter than line from input_2";
}

/// Where the comparator writes its log.
///
/// The source holds a `std::ostream*` defaulting to `std::cout`. `Stdout` and
/// `Stderr` go through `print!`/`eprint!` (so the test harness captures them) with
/// non-UTF-8 bytes replaced; `Buffer` keeps the exact bytes, see
/// [`FuzzyStringComparator::log`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum LogDestination {
    /// Standard output, the source default.
    #[default]
    Stdout,
    /// Standard error.
    Stderr,
    /// An in-memory byte buffer, the `std::ostringstream` of the class test.
    Buffer,
}

/// Result of comparing one pair of numbers under the comparator's rule.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NumberComparison {
    /// The source failure message, or `None` when the pair is accepted.
    pub failure: Option<&'static str>,
    /// The absolute difference the source records into `absolute_max`; `None`
    /// for equal or non-finite pairs, where the source returns before computing it.
    pub absdiff: Option<f64>,
    /// The ratio (reciprocal taken when below one) the source considers for the
    /// relative maximum; `None` unless both numbers are non-zero and the quotient
    /// is not negative.
    pub ratio: Option<f64>,
}

/// Compare two numbers exactly as the number branch of `compareLines_` does
/// (`FuzzyStringComparator.cpp` 556-689).
///
/// Equal numbers (including `0 == -0` and equal infinities) are accepted. Two NaNs
/// are accepted; NaN against anything else fails, as do opposite infinities and an
/// infinity against a finite number. Otherwise the pair is accepted when the
/// absolute difference is at most `absdiff_max_allowed`, or when both numbers are
/// non-zero with the same sign and their ratio (taken as at least one) does not
/// exceed `ratio_max_allowed`.
///
/// `ratio_max` is the comparator's running `ratio_max_`. The source reports a
/// ratio failure only when the ratio also exceeds it, and raises it only on such a
/// failure; a fresh comparison passes `1.0`. Two source quirks are kept because
/// they decide verdicts: a quotient that underflows to `-0.0` (opposite signs, e.g.
/// `-1e-300` against `1e300`) passes the sign test and its reciprocal `-inf` never
/// exceeds the maximum, so the pair is accepted; and a NaN tolerance accepts every
/// ratio.
pub fn compare_numbers(
    n1: f64,
    n2: f64,
    ratio_max_allowed: f64,
    absdiff_max_allowed: f64,
    ratio_max: f64,
) -> NumberComparison {
    let mut outcome = NumberComparison {
        failure: None,
        absdiff: None,
        ratio: None,
    };
    if n1 == n2 {
        return outcome;
    }
    if n1.is_nan() || n2.is_nan() {
        if !(n1.is_nan() && n2.is_nan()) {
            outcome.failure = Some(message::ONE_NAN);
        }
        return outcome;
    }
    if n1.is_infinite() || n2.is_infinite() {
        outcome.failure = Some(if n1.is_infinite() && n2.is_infinite() {
            message::INFINITY_SIGNS
        } else {
            message::ONE_INFINITE
        });
        return outcome;
    }
    let absdiff = (n1 - n2).abs();
    outcome.absdiff = Some(absdiff);
    let small = absdiff <= absdiff_max_allowed;
    if n1 == 0.0 {
        if !small {
            outcome.failure = Some(message::FIRST_ZERO);
        }
        return outcome;
    }
    if n2 == 0.0 {
        if !small {
            outcome.failure = Some(message::SECOND_ZERO);
        }
        return outcome;
    }
    let mut ratio = n1 / n2;
    if ratio < 0.0 {
        if !small {
            outcome.failure = Some(message::SIGNS);
        }
        return outcome;
    }
    if ratio < 1.0 {
        ratio = 1.0 / ratio;
    }
    outcome.ratio = Some(ratio);
    if ratio > ratio_max && ratio > ratio_max_allowed && !small {
        outcome.failure = Some(message::RATIO);
    }
    outcome
}

/// Normalise a relative tolerance as `setAcceptableRelative` does: values below
/// one are replaced by their reciprocal.
pub fn normalise_relative(ratio: f64) -> f64 {
    if ratio < 1.0 { 1.0 / ratio } else { ratio }
}

/// Normalise an absolute tolerance as `setAcceptableAbsolute` does: negative
/// values are negated.
pub fn normalise_absolute(absdiff: f64) -> f64 {
    if absdiff < 0.0 { -absdiff } else { absdiff }
}

/// Whitespace in the C locale, the set `isspace` and `std::ws` use in the source.
pub fn is_c_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// Parse the number token at the start of `input` as the source's
/// `extractDouble` does, returning the value and the number of bytes consumed.
///
/// The accepted grammar is the one the source states for `std::from_chars`
/// (general format): an optional single leading `+` (skipped by `extractDouble`),
/// then an optional `-`, then decimal digits with an optional fraction and
/// exponent, `inf`/`infinity`, or `nan` with an optional `(n-char-sequence)`,
/// case-insensitive; `nan(...)` with arbitrary contents up to `)` is accepted
/// first, before the sign handling. An exponent without digits is not consumed.
/// Values that overflow to infinity are rejected (the token then reads as
/// letters); values that underflow are accepted, as the source comment requires.
///
/// Hexadecimal floats (`0x1p3`) are not numbers here, as with `std::from_chars`.
/// The libc++ build of the source (Apple, the product-sdk oracle) parses through
/// a `strtod` fallback that does accept them; see the support document.
pub fn extract_double(input: &[u8]) -> Option<(f64, usize)> {
    if input.is_empty() {
        return None;
    }
    if let Some(consumed) = try_parse_nan(input) {
        return Some((f64::NAN, consumed));
    }
    let offset = usize::from(input.first() == Some(&b'+'));
    let rest = input.get(offset..).unwrap_or(&[]);
    let (value, consumed) = from_chars_general(rest)?;
    Some((value, offset + consumed))
}

fn try_parse_nan(input: &[u8]) -> Option<usize> {
    let head = input.get(..3)?;
    if !head.eq_ignore_ascii_case(b"nan") {
        return None;
    }
    if input.get(3) != Some(&b'(') {
        return Some(3);
    }
    let close = input.iter().skip(4).position(|&b| b == b')')?;
    Some(4 + close + 1)
}

fn from_chars_general(input: &[u8]) -> Option<(f64, usize)> {
    let sign = usize::from(input.first() == Some(&b'-'));
    let body = input.get(sign..).unwrap_or(&[]);
    let negative = sign == 1;
    let lower =
        |len: usize| -> Vec<u8> { body.iter().take(len).map(u8::to_ascii_lowercase).collect() };
    if lower(8) == b"infinity" {
        return Some((signed(f64::INFINITY, negative), sign + 8));
    }
    if lower(3) == b"inf" {
        return Some((signed(f64::INFINITY, negative), sign + 3));
    }
    if lower(3) == b"nan" {
        let mut consumed = 3;
        if body.get(3) == Some(&b'(') {
            let inner = body
                .iter()
                .skip(4)
                .take_while(|b| b.is_ascii_alphanumeric() || **b == b'_')
                .count();
            if body.get(4 + inner) == Some(&b')') {
                consumed = 4 + inner + 1;
            }
        }
        return Some((f64::NAN, sign + consumed));
    }
    let digits = |from: usize| {
        body.iter()
            .skip(from)
            .take_while(|b| b.is_ascii_digit())
            .count()
    };
    let integer = digits(0);
    let mut end = integer;
    let mut fraction = 0;
    if body.get(end) == Some(&b'.') {
        fraction = digits(end + 1);
        end += 1 + fraction;
    }
    if integer + fraction == 0 {
        return None;
    }
    if matches!(body.get(end), Some(b'e' | b'E')) {
        let exponent_sign = usize::from(matches!(body.get(end + 1), Some(b'+' | b'-')));
        let exponent_digits = digits(end + 1 + exponent_sign);
        if exponent_digits > 0 {
            end += 1 + exponent_sign + exponent_digits;
        }
    }
    let text = std::str::from_utf8(input.get(..sign + end)?).ok()?;
    let value: f64 = text.parse().ok()?;
    if value.is_infinite() {
        return None;
    }
    Some((value, sign + end))
}

fn signed(value: f64, negative: bool) -> f64 {
    if negative { -value } else { value }
}

/// Format a double as `std::ostream << double` does by default (`%g`, six
/// significant digits). NaN prints as `nan` whatever its sign, as the Apple libc
/// oracle does.
pub fn format_g(value: f64) -> String {
    if value.is_nan() {
        return "nan".to_owned();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.to_owned();
    }
    if value == 0.0 {
        return if value.is_sign_negative() { "-0" } else { "0" }.to_owned();
    }
    const PRECISION: i32 = 6;
    let scientific = format!("{:.5e}", value);
    let Some((mantissa, exponent)) = scientific.split_once('e') else {
        return scientific;
    };
    let exponent: i32 = exponent.parse().unwrap_or(0);
    if !(-4..PRECISION).contains(&exponent) {
        let sign = if exponent < 0 { '-' } else { '+' };
        format!(
            "{}e{}{:02}",
            trim_fraction(mantissa),
            sign,
            exponent.unsigned_abs()
        )
    } else {
        let decimals = usize::try_from(PRECISION - 1 - exponent).unwrap_or(0);
        trim_fraction(&format!("{:.*}", decimals, value))
    }
}

fn trim_fraction(text: &str) -> String {
    if text.contains('.') {
        text.trim_end_matches('0').trim_end_matches('.').to_owned()
    } else {
        text.to_owned()
    }
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    needle.is_empty()
        || haystack
            .windows(needle.len())
            .any(|window| window == needle)
}

/// Emulates the state of the `std::stringstream` the source reads one line from
/// (`InputLine`): the get position, `eofbit`, `failbit`/`badbit` and the saved
/// position.
#[derive(Clone, Debug, Default)]
struct InputLine {
    line: Vec<u8>,
    pos: usize,
    eof: bool,
    fail: bool,
    line_position: usize,
}

impl InputLine {
    fn set_to_string(&mut self, bytes: &[u8]) {
        self.line.clear();
        self.line.extend_from_slice(bytes);
        self.pos = 0;
        self.eof = false;
        self.fail = false;
        self.line_position = 0;
    }

    fn good(&self) -> bool {
        !self.eof && !self.fail
    }

    fn ok(&self) -> bool {
        !self.fail
    }

    /// `tellg`: a stream that is not good gets `failbit` and reports -1.
    fn tellg(&mut self) -> Option<usize> {
        if self.good() {
            Some(self.pos)
        } else {
            self.fail = true;
            None
        }
    }

    fn update_position(&mut self) {
        self.line_position = self.tellg().unwrap_or(self.line.len());
    }

    /// `line_ >> letter` with `skipws` unset.
    fn read_letter(&mut self, letter: &mut u8) {
        if !self.good() {
            self.fail = true;
            return;
        }
        match self.line.get(self.pos) {
            Some(&byte) => {
                *letter = byte;
                self.pos += 1;
            }
            None => {
                self.eof = true;
                self.fail = true;
            }
        }
    }

    /// `line_ >> std::ws`.
    fn skip_whitespace(&mut self) {
        if !self.good() {
            self.fail = true;
            return;
        }
        while let Some(&byte) = self.line.get(self.pos) {
            if !is_c_space(byte) {
                return;
            }
            self.pos += 1;
        }
        self.eof = true;
    }

    /// `seekg(pos)`: clears `eofbit` first and does nothing on a failed stream.
    fn seekg(&mut self, position: usize) {
        self.eof = false;
        if self.fail {
            return;
        }
        if position <= self.line.len() {
            self.pos = position;
        } else {
            self.fail = true;
        }
    }

    /// `clear()` followed by `seekg(line_position_)`.
    fn seek_to_saved_position(&mut self) {
        self.eof = false;
        self.fail = false;
        self.seekg(self.line_position);
    }
}

/// One character, number or whitespace run read from a line (`StreamElement_`).
#[derive(Clone, Copy, Debug)]
struct StreamElement {
    number: f64,
    letter: u8,
    is_number: bool,
    is_space: bool,
}

impl Default for StreamElement {
    fn default() -> Self {
        Self {
            number: 0.0,
            letter: 0,
            is_number: false,
            is_space: false,
        }
    }
}

impl StreamElement {
    fn reset(&mut self) {
        self.is_number = false;
        self.is_space = false;
        self.letter = 0;
        self.number = f64::NAN;
    }

    /// `fillFromInputLine`. An exhausted line yields a non-space, non-number NUL
    /// letter and leaves the stream failed, exactly as the source's stream does.
    fn fill_from_input_line(&mut self, input: &mut InputLine) {
        self.reset();
        input.update_position();
        input.read_letter(&mut self.letter);
        self.is_space = is_c_space(self.letter);
        if self.is_space {
            input.skip_whitespace();
            return;
        }
        input.seek_to_saved_position();
        let start = input.line_position;
        let parsed = extract_double(input.line.get(start..).unwrap_or(&[]));
        match parsed {
            Some((value, consumed)) => {
                self.number = value;
                self.is_number = true;
                input.seekg(start + consumed);
            }
            None => input.read_letter(&mut self.letter),
        }
    }
}

/// Column bookkeeping for the failure report (`PrefixInfo_`).
struct PrefixInfo {
    prefix: Vec<u8>,
    prefix_whitespaces: Vec<u8>,
    line_column: i32,
}

impl PrefixInfo {
    fn new(input: &InputLine, tab_width: i32, first_column: i32) -> Self {
        let length = input.line_position.min(input.line.len());
        let prefix = input.line.get(..length).unwrap_or(&[]).to_vec();
        let mut prefix_whitespaces = prefix.clone();
        let mut line_column: i32 = 0;
        for byte in &mut prefix_whitespaces {
            if *byte == b'\t' {
                // Source divides by the tab width; a zero width is undefined there
                // (arm64 yields a zero quotient) and gives a zero quotient here.
                let quotient = line_column.checked_div(tab_width).unwrap_or(0);
                line_column = quotient.wrapping_add(1).wrapping_mul(tab_width);
            } else {
                *byte = b' ';
                line_column = line_column.wrapping_add(1);
            }
        }
        line_column = line_column.wrapping_add(first_column);
        Self {
            prefix,
            prefix_whitespaces,
            line_column,
        }
    }
}

/// Signals that a failure report ended the comparison of a line (source
/// `AbortComparison`).
struct AbortComparison;

/// Reads lines as the source's local `getLine` does: `\n`, `\r\n` and a lone
/// `\r` all terminate a line; the terminator is not stored. Tracks `eofbit` and
/// `failbit` so that `while (input_1 || input_2)` behaves as in the source.
struct LineSource<'a, R: BufRead> {
    reader: &'a mut R,
    eof: bool,
    fail: bool,
    total: u64,
    error: Option<String>,
}

impl<'a, R: BufRead> LineSource<'a, R> {
    fn new(reader: &'a mut R) -> Self {
        Self {
            reader,
            eof: false,
            fail: false,
            total: 0,
            error: None,
        }
    }

    fn ok(&self) -> bool {
        !self.fail
    }

    fn fill(&mut self) -> Option<&[u8]> {
        // Retry interrupted reads without holding the buffer across iterations,
        // then borrow the (now buffered) data once.
        loop {
            match self.reader.fill_buf() {
                Ok(_) => break,
                Err(error) if error.kind() == ErrorKind::Interrupted => continue,
                Err(error) => {
                    self.error = Some(format!("Error: reading input failed: {error}"));
                    return None;
                }
            }
        }
        match self.reader.fill_buf() {
            Ok(buffer) => Some(buffer),
            Err(error) => {
                self.error = Some(format!("Error: reading input failed: {error}"));
                None
            }
        }
    }

    fn consume(&mut self, amount: usize) -> bool {
        self.reader.consume(amount);
        self.total = self.total.saturating_add(amount as u64);
        if self.total > MAX_INPUT_BYTES {
            self.error = Some(format!(
                "Error: input exceeds the comparison limit of {MAX_INPUT_BYTES} bytes"
            ));
            return false;
        }
        true
    }

    /// Returns whether the stream is still usable (`!fail()`) after the read.
    fn get_line(&mut self, line: &mut Vec<u8>) -> bool {
        line.clear();
        if self.eof || self.fail || self.error.is_some() {
            self.fail = true;
            return false;
        }
        loop {
            let (consumed, terminator) = {
                let Some(buffer) = self.fill() else {
                    self.fail = true;
                    return false;
                };
                if buffer.is_empty() {
                    self.eof = true;
                    if line.is_empty() {
                        self.fail = true;
                    }
                    return !self.fail;
                }
                match buffer.iter().position(|&b| b == b'\n' || b == b'\r') {
                    Some(index) => {
                        line.extend_from_slice(&buffer[..index]);
                        (index + 1, Some(buffer[index]))
                    }
                    None => {
                        line.extend_from_slice(buffer);
                        (buffer.len(), None)
                    }
                }
            };
            if !self.consume(consumed) {
                self.fail = true;
                return false;
            }
            match terminator {
                Some(b'\r') => {
                    let Some(next) = self.fill() else {
                        self.fail = true;
                        return false;
                    };
                    if next.first() == Some(&b'\n') && !self.consume(1) {
                        self.fail = true;
                        return false;
                    }
                    return true;
                }
                Some(_) => return true,
                None => {}
            }
        }
    }

    /// `readNextLine_`: skip empty and whitespace-only lines, counting every read.
    fn read_next_line(&mut self, line: &mut Vec<u8>, line_number: &mut i64) {
        line.clear();
        loop {
            *line_number = line_number.wrapping_add(1);
            if !self.get_line(line) {
                return;
            }
            if line.iter().any(|&byte| !is_c_space(byte)) {
                return;
            }
        }
    }
}

/// Fuzzy comparison of strings, tolerating numeric differences
/// (`OpenMS::FuzzyStringComparator`).
///
/// Lines are compared element by element: a run of whitespace equals any other
/// run, numbers are compared with [`compare_numbers`], every other byte must be
/// equal. Empty and whitespace-only lines are skipped. A line is skipped when both
/// sides contain the same whitelist entry, or one side contains the first and the
/// other the second string of a matched-whitelist pair.
///
/// The source class is neither copyable nor assignable (both operators are
/// declared but not implemented); this type does not implement `Clone`. Its state
/// is kept between comparisons exactly as in the source: line numbers, the
/// relative and absolute maxima and the whitelist counts are never reset, only the
/// success flag is. Reusing one instance therefore changes later reports and, via
/// the relative maximum, later verdicts; use a fresh instance per comparison.
#[derive(Debug)]
pub struct FuzzyStringComparator {
    log_destination: LogDestination,
    log: Vec<u8>,
    log_truncated: bool,
    reports: usize,
    input_1_name: String,
    input_2_name: String,
    input_line_1: InputLine,
    input_line_2: InputLine,
    line_num_1: i64,
    line_num_2: i64,
    line_num_1_max: i64,
    line_num_2_max: i64,
    line_str_1_max: Vec<u8>,
    line_str_2_max: Vec<u8>,
    ratio_max_allowed: f64,
    ratio_max: f64,
    absdiff_max_allowed: f64,
    absdiff_max: f64,
    element_1: StreamElement,
    element_2: StreamElement,
    verbose_level: i32,
    tab_width: i32,
    first_column: i32,
    is_status_success: bool,
    use_prefix: bool,
    whitelist: Vec<String>,
    whitelist_cases: BTreeMap<String, u32>,
    matched_whitelist: Vec<(String, String)>,
}

impl Default for FuzzyStringComparator {
    fn default() -> Self {
        Self::new()
    }
}

impl FuzzyStringComparator {
    /// A comparator with the source defaults: relative tolerance 1, absolute
    /// tolerance 0, verbose level 2, tab width 8, first column 1, empty whitelists
    /// and the log on standard output.
    pub fn new() -> Self {
        Self {
            log_destination: LogDestination::Stdout,
            log: Vec::new(),
            log_truncated: false,
            reports: 0,
            input_1_name: "input_1".to_owned(),
            input_2_name: "input_2".to_owned(),
            input_line_1: InputLine::default(),
            input_line_2: InputLine::default(),
            line_num_1: 0,
            line_num_2: 0,
            line_num_1_max: -1,
            line_num_2_max: -1,
            line_str_1_max: Vec::new(),
            line_str_2_max: Vec::new(),
            ratio_max_allowed: 1.0,
            ratio_max: 1.0,
            absdiff_max_allowed: 0.0,
            absdiff_max: 0.0,
            element_1: StreamElement::default(),
            element_2: StreamElement::default(),
            verbose_level: 2,
            tab_width: 8,
            first_column: 1,
            is_status_success: true,
            use_prefix: false,
            whitelist: Vec::new(),
            whitelist_cases: BTreeMap::new(),
            matched_whitelist: Vec::new(),
        }
    }

    /// Acceptable relative error, a number of at least one (`getAcceptableRelative`).
    pub fn acceptable_relative(&self) -> f64 {
        self.ratio_max_allowed
    }

    /// Set the acceptable relative error (`setAcceptableRelative`); a value below
    /// one is replaced by its reciprocal. NaN is stored unchanged and then accepts
    /// every ratio, as in the source.
    pub fn set_acceptable_relative(&mut self, ratio: f64) {
        self.ratio_max_allowed = normalise_relative(ratio);
    }

    /// Acceptable absolute difference, a number of at least zero
    /// (`getAcceptableAbsolute`).
    pub fn acceptable_absolute(&self) -> f64 {
        self.absdiff_max_allowed
    }

    /// Set the acceptable absolute difference (`setAcceptableAbsolute`); a negative
    /// value is negated.
    pub fn set_acceptable_absolute(&mut self, absdiff: f64) {
        self.absdiff_max_allowed = normalise_absolute(absdiff);
    }

    /// The whitelist: lines that both contain the same entry are skipped
    /// (`getWhitelist() const`).
    pub fn whitelist(&self) -> &[String] {
        &self.whitelist
    }

    /// Mutable access to the whitelist (`getWhitelist()`).
    pub fn whitelist_mut(&mut self) -> &mut Vec<String> {
        &mut self.whitelist
    }

    /// Replace the whitelist (`setWhitelist`). Entries are tried in order and the
    /// first one found in both lines is counted; an empty entry matches every line.
    pub fn set_whitelist(&mut self, whitelist: Vec<String>) {
        self.whitelist = whitelist;
    }

    /// The matched whitelist (`getMatchedWhitelist`).
    pub fn matched_whitelist(&self) -> &[(String, String)] {
        &self.matched_whitelist
    }

    /// Replace the matched whitelist (`setMatchedWhitelist`): a line pair is
    /// skipped when one line contains the first string of a pair and the other
    /// line the second, in either order. Matches are not counted.
    pub fn set_matched_whitelist(&mut self, matched_whitelist: Vec<(String, String)>) {
        self.matched_whitelist = matched_whitelist;
    }

    /// Verbose level (`getVerboseLevel`): 0 very quiet (no output), 1 quiet (output
    /// only on differences), 2 default (summary at the end), 3 continue after
    /// errors. Any value is stored; below 1 behaves as 0 and above 3 as 3.
    pub fn verbose_level(&self) -> i32 {
        self.verbose_level
    }

    /// Set the verbose level (`setVerboseLevel`).
    pub fn set_verbose_level(&mut self, verbose_level: i32) {
        self.verbose_level = verbose_level;
    }

    /// Tab width used for column numbers in the failure report (`getTabWidth`).
    pub fn tab_width(&self) -> i32 {
        self.tab_width
    }

    /// Set the tab width (`setTabWidth`). Any value is stored.
    pub fn set_tab_width(&mut self, tab_width: i32) {
        self.tab_width = tab_width;
    }

    /// Number of the first column in the failure report (`getFirstColumn`).
    pub fn first_column(&self) -> i32 {
        self.first_column
    }

    /// Set the number of the first column (`setFirstColumn`).
    pub fn set_first_column(&mut self, first_column: i32) {
        self.first_column = first_column;
    }

    /// Where log output is written (`getLogDestination`).
    pub fn log_destination(&self) -> LogDestination {
        self.log_destination
    }

    /// Select where log output is written (`setLogDestination`). Bytes already in
    /// the buffer are kept.
    pub fn set_log_destination(&mut self, destination: LogDestination) {
        self.log_destination = destination;
    }

    /// The bytes written while the destination was [`LogDestination::Buffer`].
    pub fn log(&self) -> &[u8] {
        &self.log
    }

    /// Take and clear the buffered log (the class test's `log.str("")`).
    pub fn take_log(&mut self) -> Vec<u8> {
        self.log_truncated = false;
        std::mem::take(&mut self.log)
    }

    /// Whether buffered output was dropped at [`MAX_LOG_BYTES`].
    pub fn log_truncated(&self) -> bool {
        self.log_truncated
    }

    fn emit(&mut self, bytes: &[u8]) {
        match self.log_destination {
            LogDestination::Buffer => {
                let room = MAX_LOG_BYTES.saturating_sub(self.log.len());
                if bytes.len() > room {
                    self.log.extend_from_slice(&bytes[..room]);
                    self.log_truncated = true;
                } else {
                    self.log.extend_from_slice(bytes);
                }
            }
            LogDestination::Stdout => print!("{}", String::from_utf8_lossy(bytes)),
            LogDestination::Stderr => eprint!("{}", String::from_utf8_lossy(bytes)),
        }
    }

    /// Compare two strings line by line (`compareStrings`); true when no
    /// difference was found.
    pub fn compare_strings(&mut self, lhs: &str, rhs: &str) -> bool {
        self.compare_bytes(lhs.as_bytes(), rhs.as_bytes())
    }

    /// Compare two byte strings line by line, as `compareStrings` does with the
    /// bytes of a `std::string`.
    pub fn compare_bytes(&mut self, lhs: &[u8], rhs: &[u8]) -> bool {
        self.compare_streams(&mut Cursor::new(lhs), &mut Cursor::new(rhs))
    }

    /// Compare two inputs line by line (`compareStreams`); true when no difference
    /// was found.
    ///
    /// Resets only the success flag, then reads one non-blank line from each
    /// input per step until both are exhausted; a missing line compares as an
    /// empty one. Stops after the first failing line pair unless the verbose level
    /// is at least 3, then writes the success report.
    ///
    /// Native difference: a read error or an input beyond [`MAX_INPUT_BYTES`]
    /// fails the comparison with a log message; the source's stream would treat
    /// the error as the end of input.
    pub fn compare_streams<R1: BufRead, R2: BufRead>(
        &mut self,
        input_1: &mut R1,
        input_2: &mut R2,
    ) -> bool {
        self.is_status_success = true;
        self.reports = 0;
        let mut source_1 = LineSource::new(input_1);
        let mut source_2 = LineSource::new(input_2);
        let mut line_1 = Vec::new();
        let mut line_2 = Vec::new();
        while source_1.ok() || source_2.ok() {
            source_1.read_next_line(&mut line_1, &mut self.line_num_1);
            source_2.read_next_line(&mut line_2, &mut self.line_num_2);
            if let Some(error) = source_1.error.take().or_else(|| source_2.error.take()) {
                self.is_status_success = false;
                self.emit(format!("{error}\n").as_bytes());
                return false;
            }
            if !self.compare_lines(&line_1, &line_2) && self.verbose_level < 3 {
                break;
            }
        }
        self.report_success();
        self.is_status_success
    }

    /// Diff-like comparison of two files (`compareFiles`); true when no
    /// difference was found.
    ///
    /// Fails without comparing when both names are the same string ("That's
    /// cheating!"), or when a file cannot be opened; the source reports both open
    /// failures as "Error opening first input file". Files are read in binary mode.
    ///
    /// Native difference: a file larger than [`MAX_INPUT_BYTES`] is refused before
    /// reading.
    pub fn compare_files(&mut self, filename_1: &Path, filename_2: &Path) -> bool {
        self.input_1_name = filename_1.to_string_lossy().into_owned();
        self.input_2_name = filename_2.to_string_lossy().into_owned();
        if filename_1.as_os_str() == filename_2.as_os_str() {
            self.emit(b"Error: first and second input file have the same name. That's cheating!\n");
            return false;
        }
        let Some(file_1) = self.open_input_file(filename_1) else {
            return false;
        };
        let Some(file_2) = self.open_input_file(filename_2) else {
            return false;
        };
        self.compare_streams(&mut BufReader::new(file_1), &mut BufReader::new(file_2));
        self.is_status_success
    }

    fn open_input_file(&mut self, filename: &Path) -> Option<File> {
        let Ok(file) = File::open(filename) else {
            let text = format!(
                "Error opening first input file '{}'.\n",
                filename.to_string_lossy()
            );
            self.emit(text.as_bytes());
            return None;
        };
        let length = file.metadata().map(|m| m.len()).unwrap_or(0);
        if length > MAX_INPUT_BYTES {
            let text = format!(
                "Error: input file '{}' exceeds the comparison limit of {MAX_INPUT_BYTES} bytes.\n",
                filename.to_string_lossy()
            );
            self.emit(text.as_bytes());
            return None;
        }
        Some(file)
    }

    /// Compare two single lines (`compareLines_`, protected in the source and
    /// public here so single-line cases can be tested directly); returns the
    /// success flag after the comparison.
    ///
    /// Identical lines return true immediately without touching any state. A
    /// whitelist or matched-whitelist hit returns the current success flag. A
    /// carriage return that starts a whitespace run on one side is skipped when the
    /// other side holds a letter; `compare_streams` never passes one, because line
    /// reading already treats it as a terminator.
    pub fn compare_lines(&mut self, line_1: &[u8], line_2: &[u8]) -> bool {
        if line_1 == line_2 {
            return true;
        }
        for entry in &self.whitelist {
            if contains(line_1, entry.as_bytes()) && contains(line_2, entry.as_bytes()) {
                let count = self.whitelist_cases.entry(entry.clone()).or_insert(0);
                *count = count.wrapping_add(1);
                return self.is_status_success;
            }
        }
        for (first, second) in &self.matched_whitelist {
            let (first, second) = (first.as_bytes(), second.as_bytes());
            if (contains(line_1, first) && contains(line_2, second))
                || (contains(line_1, second) && contains(line_2, first))
            {
                return self.is_status_success;
            }
        }
        self.input_line_1.set_to_string(line_1);
        self.input_line_2.set_to_string(line_2);
        let mut maximum_on_this_line = false;
        let _ = self.compare_elements(&mut maximum_on_this_line);
        if maximum_on_this_line {
            // The source copies both lines at every new maximum; recording the flag
            // and copying once gives the same report without quadratic copying.
            self.line_str_1_max = line_1.to_vec();
            self.line_str_2_max = line_2.to_vec();
        }
        self.is_status_success
    }

    fn compare_elements(&mut self, maximum_on_this_line: &mut bool) -> Result<(), AbortComparison> {
        while self.input_line_1.ok() && self.input_line_2.ok() {
            self.element_1.fill_from_input_line(&mut self.input_line_1);
            self.element_2.fill_from_input_line(&mut self.input_line_2);
            let (e1, e2) = (self.element_1, self.element_2);
            if e1.is_number {
                if !e2.is_number {
                    self.report_failure(message::NUMBER_FIRST)?;
                    continue;
                }
                let outcome = compare_numbers(
                    e1.number,
                    e2.number,
                    self.ratio_max_allowed,
                    self.absdiff_max_allowed,
                    self.ratio_max,
                );
                if let Some(absdiff) = outcome.absdiff {
                    if absdiff > self.absdiff_max {
                        self.absdiff_max = absdiff;
                    }
                }
                if let Some(ratio) = outcome.ratio {
                    if ratio > self.ratio_max {
                        self.line_num_1_max = self.line_num_1;
                        self.line_num_2_max = self.line_num_2;
                        *maximum_on_this_line = true;
                        if outcome.failure.is_some() {
                            self.ratio_max = ratio;
                        }
                    }
                }
                if let Some(text) = outcome.failure {
                    self.report_failure(text)?;
                }
                continue;
            }
            if e2.is_number {
                self.report_failure(message::NUMBER_SECOND)?;
                continue;
            }
            if e1.is_space {
                if e2.is_space {
                    continue;
                }
                if e1.letter == b'\r' {
                    self.input_line_2.seek_to_saved_position();
                    continue;
                }
                self.report_failure(message::SPACE_FIRST)?;
                continue;
            }
            if e2.is_space {
                if e2.letter == b'\r' {
                    self.input_line_1.seek_to_saved_position();
                    continue;
                }
                self.report_failure(message::SPACE_SECOND)?;
                continue;
            }
            if e1.letter != e2.letter {
                self.report_failure(message::LETTERS)?;
            }
        }
        if self.input_line_1.ok() && !self.input_line_2.ok() {
            self.report_failure(message::SECOND_SHORTER)?;
        }
        if !self.input_line_1.ok() && self.input_line_2.ok() {
            self.report_failure(message::FIRST_SHORTER)?;
        }
        Ok(())
    }

    fn report_failure(&mut self, text: &str) -> Result<(), AbortComparison> {
        self.is_status_success = false;
        if self.verbose_level >= 1 && !self.log_truncated {
            self.reports += 1;
            if self.reports <= MAX_FAILURE_REPORTS {
                let report = self.failure_report(text);
                self.emit(&report);
            } else if self.reports == MAX_FAILURE_REPORTS + 1 {
                self.emit(b"Further failure reports are suppressed.\n");
            }
        }
        if self.verbose_level < 3 {
            return Err(AbortComparison);
        }
        Ok(())
    }

    fn prefix(&self) -> &'static str {
        if self.use_prefix { "   :|:  " } else { "" }
    }

    fn failure_report(&self, text: &str) -> Vec<u8> {
        let p = self.prefix();
        let prefix_1 = PrefixInfo::new(&self.input_line_1, self.tab_width, self.first_column);
        let prefix_2 = PrefixInfo::new(&self.input_line_2, self.tab_width, self.first_column);
        let (e1, e2) = (self.element_1, self.element_2);
        let (pos_1, pos_2) = (
            self.input_line_1.line_position,
            self.input_line_2.line_position,
        );
        let (abs_1, abs_2) = (
            absolute_display(&self.input_1_name),
            absolute_display(&self.input_2_name),
        );
        let mut out = Vec::new();
        let mut put = |s: &str| out.extend_from_slice(s.as_bytes());
        put(&format!(
            "{p}FAILED: '{text}'\n{p}\n{p}  input:\tin1\tin2\n"
        ));
        put(&format!(
            "{p}  line:\t{}\t{}\n",
            self.line_num_1, self.line_num_2
        ));
        put(&format!(
            "{p}  pos/col:\t{pos_1}/{}\t{pos_2}/{}\n{p} --------------------------------\n",
            prefix_1.line_column, prefix_2.line_column
        ));
        put(&format!(
            "{p}  is_number:\t{}\t{}\n",
            e1.is_number, e2.is_number
        ));
        put(&format!(
            "{p}  numbers:\t{}\t{}\n",
            format_g(e1.number),
            format_g(e2.number)
        ));
        put(&format!(
            "{p}  is_space:\t{}\t{}\n",
            e1.is_space, e2.is_space
        ));
        put(&format!(
            "{p}  is_letter:\t{}\t{}\n",
            !e1.is_number && !e1.is_space,
            !e2.is_number && !e2.is_space
        ));
        put(&format!("{p}  letters:\t\""));
        out.push(e1.letter);
        out.extend_from_slice(b"\"\t\"");
        out.push(e2.letter);
        let mut put = |s: &str| out.extend_from_slice(s.as_bytes());
        put("\"\n");
        put(&format!(
            "{p}  char_codes:\t{}\t{}\n{p} --------------------------------\n",
            e1.letter, e2.letter
        ));
        put(&format!(
            "{p}  relative_max:        {}\n{p}  relative_acceptable: {}\n{p} --------------------------------\n",
            format_g(self.ratio_max),
            format_g(self.ratio_max_allowed)
        ));
        put(&format!(
            "{p}  absolute_max:        {}\n{p}  absolute_acceptable: {}\n",
            format_g(self.absdiff_max),
            format_g(self.absdiff_max_allowed)
        ));
        out.extend_from_slice(&self.whitelist_cases_text());
        let mut put = |s: &str| out.extend_from_slice(s.as_bytes());
        put(&format!(
            "{p}\n{p}Offending lines:\t\t\t(tab_width = {}, first_column = {})\n{p}\n",
            self.tab_width, self.first_column
        ));
        for (name, number, position, info, line) in [
            (
                "in1",
                self.line_num_1,
                pos_1,
                &prefix_1,
                &self.input_line_1.line,
            ),
            (
                "in2",
                self.line_num_2,
                pos_2,
                &prefix_2,
                &self.input_line_2.line,
            ),
        ] {
            let absolute = if name == "in1" { &abs_1 } else { &abs_2 };
            out.extend_from_slice(
                format!(
                    "{p}{name}:  {absolute}   (line: {number}, position/column: {position}/{})\n{p}",
                    info.line_column
                )
                .as_bytes(),
            );
            out.extend_from_slice(&info.prefix);
            out.extend_from_slice(format!("!\n{p}").as_bytes());
            out.extend_from_slice(&info.prefix_whitespaces);
            out.extend_from_slice(line.get(info.prefix.len()..).unwrap_or(&[]));
            out.extend_from_slice(format!("\n{p}\n").as_bytes());
        }
        // The loop above wrote "{p}\n" after in2; the source writes "{p}\n\n" there.
        out.extend_from_slice(b"\n");
        out.extend_from_slice(
            format!(
                "Easy Access:\n{abs_1}:{}:{}:\n{abs_2}:{}:{}:\n\ndiff {abs_1} {abs_2}\n",
                self.line_num_1, prefix_1.line_column, self.line_num_2, prefix_2.line_column
            )
            .as_bytes(),
        );
        out
    }

    fn whitelist_cases_text(&self) -> Vec<u8> {
        let mut out = Vec::new();
        if self.whitelist_cases.is_empty() {
            return out;
        }
        let p = self.prefix();
        out.extend_from_slice(format!("{p}\n{p}  whitelist cases:\n").as_bytes());
        let width = self
            .whitelist_cases
            .keys()
            .map(String::len)
            .max()
            .unwrap_or(0)
            + 3;
        for (entry, count) in &self.whitelist_cases {
            let quoted = format!("\"{entry}\"");
            let padding = width.saturating_sub(quoted.len());
            let count = count.to_string();
            let count_padding = 3usize.saturating_sub(count.len());
            out.extend_from_slice(
                format!(
                    "{p}    {quoted}{}{}{count}x\n",
                    " ".repeat(padding),
                    " ".repeat(count_padding)
                )
                .as_bytes(),
            );
        }
        out
    }

    fn report_success(&mut self) {
        if !self.is_status_success || self.verbose_level < 2 {
            return;
        }
        let p = self.prefix();
        let mut out = format!(
            "{p}PASSED.\n{p}\n{p}  relative_max:        {}\n{p}  relative_acceptable: {}\n{p}\n{p}  absolute_max:        {}\n{p}  absolute_acceptable: {}\n",
            format_g(self.ratio_max),
            format_g(self.ratio_max_allowed),
            format_g(self.absdiff_max),
            format_g(self.absdiff_max_allowed)
        )
        .into_bytes();
        out.extend_from_slice(&self.whitelist_cases_text());
        out.extend_from_slice(format!("{p}\n").as_bytes());
        if self.line_num_1_max == -1 && self.line_num_2_max == -1 {
            out.extend_from_slice(
                format!("{p}No numeric differences were found.\n{p}\n").as_bytes(),
            );
        } else {
            out.extend_from_slice(
                format!(
                    "{p}Maximum relative error was attained at these lines, enclosed in \"\":\n{p}\n{}:{}:\n\"",
                    self.input_1_name, self.line_num_1_max
                )
                .as_bytes(),
            );
            out.extend_from_slice(&self.line_str_1_max);
            out.extend_from_slice(
                format!("\"\n\n{}:{}:\n\"", self.input_2_name, self.line_num_2_max).as_bytes(),
            );
            out.extend_from_slice(&self.line_str_2_max);
            out.extend_from_slice(b"\"\n\n");
        }
        self.emit(&out);
    }
}

/// The absolute path the failure report prints for an input name
/// (`std::filesystem::absolute`): the current directory joined with a relative
/// name, without normalisation.
fn absolute_display(name: &str) -> String {
    let path = Path::new(name);
    if !name.is_empty() && path.is_absolute() {
        return name.to_owned();
    }
    match std::env::current_dir() {
        Ok(directory) if name.is_empty() => directory.to_string_lossy().into_owned(),
        Ok(directory) => directory.join(path).to_string_lossy().into_owned(),
        Err(_) => name.to_owned(),
    }
}

/// Exit codes `FuzzyDiff` returns (the used subset of `TOPPBase::ExitCodes`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FuzzyDiffExit {
    /// No difference found (`EXECUTION_OK`).
    ExecutionOk = 0,
    /// An input file does not exist (`INPUT_FILE_NOT_FOUND`).
    InputFileNotFound = 1,
    /// An input file cannot be read (`INPUT_FILE_NOT_READABLE`).
    InputFileNotReadable = 2,
    /// An input file is empty (`INPUT_FILE_EMPTY`).
    InputFileEmpty = 4,
    /// A parameter is out of its registered range or the INI is invalid
    /// (`ILLEGAL_PARAMETERS`).
    IllegalParameters = 6,
    /// `in1` or `in2` is empty (`MISSING_PARAMETERS`).
    MissingParameters = 7,
    /// A malformed `matched_whitelist` entry (`IllegalArgument`, `UNKNOWN_ERROR`).
    UnknownError = 8,
    /// Differences were found (`PARSE_ERROR`; the source notes it should find a
    /// better code).
    ParseError = 10,
}

impl FuzzyDiffExit {
    /// The numeric process exit code.
    pub fn code(self) -> i32 {
        self as i32
    }
}

/// Outcome of a [`fuzzy_diff`] run.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FuzzyDiffOutcome {
    /// The exit code the C++ tool returns for the same invocation.
    pub exit: FuzzyDiffExit,
    /// The comparator log, or a one-line reason for a pre-comparison exit. The
    /// TOPPBase framework messages of the tool are not reproduced.
    pub log: Vec<u8>,
}

impl FuzzyDiffOutcome {
    /// Whether the comparison passed (exit code 0).
    pub fn passed(&self) -> bool {
        self.exit == FuzzyDiffExit::ExecutionOk
    }

    /// The log as text, with invalid UTF-8 replaced.
    pub fn log_text(&self) -> String {
        String::from_utf8_lossy(&self.log).into_owned()
    }
}

/// The parameters of a `FuzzyDiff` invocation (topp 174b576 `FuzzyDiff.cpp`).
#[derive(Clone, Debug, PartialEq)]
pub struct FuzzyDiffSettings {
    /// Acceptable relative error (`-ratio`, at least 1).
    pub ratio: f64,
    /// Acceptable absolute difference (`-absdiff`, at least 0).
    pub absdiff: f64,
    /// Lines containing one of these strings on both sides are skipped
    /// (`-whitelist`).
    pub whitelist: Vec<String>,
    /// Colon-separated pairs `first:second` (`-matched_whitelist`).
    pub matched_whitelist: Vec<String>,
    /// Verbose level 0-3 (`-verbose`).
    pub verbose: i32,
    /// Tab width for column numbers, at least 1 (`-tab_width`).
    pub tab_width: i32,
    /// Number of the first column, at least 0 (`-first_column`).
    pub first_column: i32,
    /// Sort all lines but the first before comparing (`-sort`).
    pub sort: bool,
    /// Problems found while reading an INI file; any entry makes [`fuzzy_diff`]
    /// return [`FuzzyDiffExit::IllegalParameters`], as TOPPBase rejects an invalid
    /// INI during initialisation.
    pub ini_errors: Vec<String>,
}

impl Default for FuzzyDiffSettings {
    fn default() -> Self {
        Self::registered_defaults()
    }
}

impl FuzzyDiffSettings {
    /// The defaults `FuzzyDiff` registers: ratio 1, absdiff 0, whitelist
    /// `<?xml-stylesheet`, no matched whitelist, verbose 2, tab width 8, first
    /// column 1, no sorting.
    pub fn registered_defaults() -> Self {
        Self {
            ratio: 1.0,
            absdiff: 0.0,
            whitelist: vec!["<?xml-stylesheet".to_owned()],
            matched_whitelist: Vec::new(),
            verbose: 2,
            tab_width: 8,
            first_column: 1,
            sort: false,
            ini_errors: Vec::new(),
        }
    }

    /// Path of the pinned upstream `FuzzyDiff.ini` (test-data 0cb15f2
    /// `topp/FuzzyDiff.ini`) copied into this repository's test data.
    pub fn upstream_ini_path() -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/fuzzy_string_comparator/FuzzyDiff.ini")
    }

    /// Settings of the upstream `${DIFF}` command: the registered defaults
    /// overridden by the pinned `FuzzyDiff.ini` (ratio 1.01, absdiff 0.01,
    /// whitelist `<?xml-stylesheet`, verbose 1).
    ///
    /// # Errors
    ///
    /// Returns the reason when the pinned file cannot be read.
    pub fn upstream() -> Result<Self, String> {
        Self::load_ini(&Self::upstream_ini_path())
    }

    /// Read a ParamXML (INI) file; see [`Self::from_ini`].
    ///
    /// # Errors
    ///
    /// Returns the reason when the file cannot be read or is not well-formed
    /// ParamXML.
    pub fn load_ini(path: &Path) -> Result<Self, String> {
        let file =
            File::open(path).map_err(|e| format!("cannot open '{}': {e}", path.display()))?;
        let mut bytes = Vec::new();
        file.take(MAX_INI_BYTES as u64 + 1)
            .read_to_end(&mut bytes)
            .map_err(|e| format!("cannot read '{}': {e}", path.display()))?;
        Self::from_ini(&bytes)
    }

    /// Apply the values of a ParamXML (INI) file to the registered defaults, as
    /// `FuzzyDiff -ini <file>` does.
    ///
    /// Items of the tool instance section `FuzzyDiff:1:` are read: `ratio`,
    /// `absdiff`, `whitelist`, `matched_whitelist`, `verbose`, `tab_width`,
    /// `first_column` and `sort`. The common TOPP items (`in1`, `in2`, `log`,
    /// `debug`, `threads`, `no_progress`, `force`, `test`) and the tool `version`
    /// item are accepted and ignored. Any other item in the instance section, or a
    /// value that does not parse, is recorded in [`Self::ini_errors`]. Items outside
    /// the tool's sections are ignored.
    ///
    /// # Errors
    ///
    /// Returns the reason when the input exceeds [`MAX_INI_BYTES`] or is not
    /// well-formed ParamXML.
    pub fn from_ini(bytes: &[u8]) -> Result<Self, String> {
        let items = param_xml::parse(bytes)?;
        let mut settings = Self::registered_defaults();
        for item in items {
            let in_instance = item.path == ["FuzzyDiff", "1"];
            if !in_instance {
                continue;
            }
            let problem = settings.apply_ini_item(&item);
            if let Err(problem) = problem {
                settings.ini_errors.push(problem);
            }
        }
        Ok(settings)
    }

    fn apply_ini_item(&mut self, item: &param_xml::Item) -> Result<(), String> {
        use param_xml::Value;
        let key = format!("FuzzyDiff:1:{}", item.name);
        let scalar = || match &item.value {
            Value::Scalar(text) => Ok(text.as_str()),
            Value::List(_) => Err(format!("'{key}' must be a single value")),
        };
        let list = || match &item.value {
            Value::List(values) => Ok(values.clone()),
            Value::Scalar(_) => Err(format!("'{key}' must be a list")),
        };
        let double = |text: &str| {
            text.trim()
                .parse::<f64>()
                .map_err(|_| format!("'{key}' is not a number: '{text}'"))
        };
        let int = |text: &str| {
            text.trim()
                .parse::<i32>()
                .map_err(|_| format!("'{key}' is not an integer: '{text}'"))
        };
        match item.name.as_str() {
            "ratio" => self.ratio = double(scalar()?)?,
            "absdiff" => self.absdiff = double(scalar()?)?,
            "whitelist" => self.whitelist = list()?,
            "matched_whitelist" => self.matched_whitelist = list()?,
            "verbose" => self.verbose = int(scalar()?)?,
            "tab_width" => self.tab_width = int(scalar()?)?,
            "first_column" => self.first_column = int(scalar()?)?,
            "sort" => {
                self.sort = match scalar()? {
                    "true" => true,
                    "false" => false,
                    other => return Err(format!("'{key}' is not a flag: '{other}'")),
                }
            }
            "in1" | "in2" | "log" | "debug" | "threads" | "no_progress" | "force" | "test" => {}
            other => {
                return Err(format!(
                    "Unknown (or deprecated) Parameter 'FuzzyDiff:1:{other}'"
                ));
            }
        }
        Ok(())
    }

    /// Replace the whitelist, as a registration's `-whitelist a b` does (the
    /// command-line list replaces the INI list; it is not appended).
    pub fn with_whitelist(mut self, entries: &[&str]) -> Self {
        self.whitelist = entries.iter().map(|&e| e.to_owned()).collect();
        self
    }

    /// Replace the matched whitelist, as `-matched_whitelist a:b` does.
    pub fn with_matched_whitelist(mut self, entries: &[&str]) -> Self {
        self.matched_whitelist = entries.iter().map(|&e| e.to_owned()).collect();
        self
    }

    /// A comparator configured with these settings, logging into a buffer.
    ///
    /// # Errors
    ///
    /// Returns the exit code and reason `FuzzyDiff` would fail with before
    /// comparing: an INI error, a parameter out of range, or a malformed matched
    /// whitelist entry.
    pub fn comparator(&self) -> Result<FuzzyStringComparator, (FuzzyDiffExit, String)> {
        if let Some(problem) = self.ini_errors.first() {
            return Err((FuzzyDiffExit::IllegalParameters, problem.clone()));
        }
        self.check_ranges()?;
        let matched = self.parsed_matched_whitelist()?;
        let mut comparator = FuzzyStringComparator::new();
        comparator.set_log_destination(LogDestination::Buffer);
        comparator.set_acceptable_relative(self.ratio);
        comparator.set_acceptable_absolute(self.absdiff);
        comparator.set_whitelist(self.whitelist.clone());
        comparator.set_matched_whitelist(matched);
        comparator.set_verbose_level(self.verbose);
        comparator.set_tab_width(self.tab_width);
        comparator.set_first_column(self.first_column);
        Ok(comparator)
    }

    /// Compare two in-memory texts with these settings: the `FuzzyDiff` comparison
    /// without its file checks and same-name check. Sorting applies when enabled.
    ///
    /// # Errors
    ///
    /// Returns the comparator log (or the pre-comparison reason) when the texts
    /// differ or the settings are invalid.
    pub fn compare_bytes(&self, actual: &[u8], expected: &[u8]) -> Result<(), String> {
        let mut comparator = self.comparator().map_err(|(_, reason)| reason)?;
        let passed = if self.sort {
            comparator.compare_bytes(&sorted_lines(actual), &sorted_lines(expected))
        } else {
            comparator.compare_bytes(actual, expected)
        };
        if passed {
            Ok(())
        } else {
            Err(String::from_utf8_lossy(comparator.log()).into_owned())
        }
    }

    #[allow(clippy::manual_range_contains)]
    fn check_ranges(&self) -> Result<(), (FuzzyDiffExit, String)> {
        // getDoubleOption_/getIntOption_ check a value only when it differs from
        // the registered default, with plain comparisons: NaN passes, as in C++.
        let defaults = Self::registered_defaults();
        let invalid = |name: &str, value: String| {
            Err((
                FuzzyDiffExit::IllegalParameters,
                format!("Invalid value '{value}' for parameter '{name}' given."),
            ))
        };
        if self.ratio != defaults.ratio && self.ratio < 1.0 {
            return invalid("ratio", self.ratio.to_string());
        }
        if self.absdiff != defaults.absdiff && self.absdiff < 0.0 {
            return invalid("absdiff", self.absdiff.to_string());
        }
        if self.verbose != defaults.verbose && (self.verbose < 0 || self.verbose > 3) {
            return invalid("verbose", self.verbose.to_string());
        }
        if self.tab_width != defaults.tab_width && self.tab_width < 1 {
            return invalid("tab_width", self.tab_width.to_string());
        }
        if self.first_column != defaults.first_column && self.first_column < 0 {
            return invalid("first_column", self.first_column.to_string());
        }
        Ok(())
    }

    fn parsed_matched_whitelist(&self) -> Result<Vec<(String, String)>, (FuzzyDiffExit, String)> {
        let mut pairs = Vec::new();
        for entry in &self.matched_whitelist {
            // StringUtils::split(entry, ":"): an empty entry yields no parts.
            let parts: Vec<&str> = if entry.is_empty() {
                Vec::new()
            } else {
                entry.split(':').collect()
            };
            match parts.as_slice() {
                [first, second] => pairs.push(((*first).to_owned(), (*second).to_owned())),
                _ => {
                    return Err((
                        FuzzyDiffExit::UnknownError,
                        format!("{entry} does not have the format String1:String2"),
                    ));
                }
            }
        }
        Ok(pairs)
    }
}

/// Run the `FuzzyDiff` tool contract on two files and return the exit code the
/// C++ tool gives for the same invocation.
///
/// Order of checks follows TOPPBase and `FuzzyDiff::main_`, as executed on the
/// oracle: during initialisation, an INI error or a parameter outside its
/// registered range, from the INI or the command line (6); then for `in1` and then
/// `in2`, an empty name (7), a missing file (1), an unreadable file (2) or an empty
/// regular file (4); then a matched-whitelist entry that does not split into
/// exactly two parts at `:` (8); then the comparison, 0 when equal and 10
/// otherwise. With `sort` the lines after the first are sorted bytewise (as
/// `std::getline` splits them, at `\n` only) and compared in memory, so the
/// same-name check does not apply, as with the source's temporary files.
pub fn fuzzy_diff(in1: &Path, in2: &Path, settings: &FuzzyDiffSettings) -> FuzzyDiffOutcome {
    let early = |exit: FuzzyDiffExit, reason: String| FuzzyDiffOutcome {
        exit,
        log: format!("{reason}\n").into_bytes(),
    };
    if let Some(problem) = settings.ini_errors.first() {
        return early(FuzzyDiffExit::IllegalParameters, problem.clone());
    }
    if let Err((exit, reason)) = settings.check_ranges() {
        return early(exit, reason);
    }
    for (name, path) in [("in1", in1), ("in2", in2)] {
        if let Err((exit, reason)) = check_input_file(name, path) {
            return early(exit, reason);
        }
    }
    let mut comparator = match settings.comparator() {
        Ok(comparator) => comparator,
        Err((exit, reason)) => return early(exit, reason),
    };
    let passed = if settings.sort {
        match (read_bounded(in1), read_bounded(in2)) {
            (Ok(a), Ok(b)) => {
                comparator.input_1_name = in1.to_string_lossy().into_owned();
                comparator.input_2_name = in2.to_string_lossy().into_owned();
                comparator.compare_bytes(&sorted_lines(&a), &sorted_lines(&b))
            }
            (Err(reason), _) | (_, Err(reason)) => {
                return early(FuzzyDiffExit::InputFileNotReadable, reason);
            }
        }
    } else {
        comparator.compare_files(in1, in2)
    };
    FuzzyDiffOutcome {
        exit: if passed {
            FuzzyDiffExit::ExecutionOk
        } else {
            FuzzyDiffExit::ParseError
        },
        log: comparator.take_log(),
    }
}

fn check_input_file(name: &str, path: &Path) -> Result<(), (FuzzyDiffExit, String)> {
    if path.as_os_str().is_empty() {
        return Err((
            FuzzyDiffExit::MissingParameters,
            format!("Error: The required parameter '{name}' was not given or is empty!"),
        ));
    }
    let Ok(metadata) = std::fs::metadata(path) else {
        return Err((
            FuzzyDiffExit::InputFileNotFound,
            format!("Error: File not found ({})", path.display()),
        ));
    };
    if !metadata.is_dir() && File::open(path).is_err() {
        return Err((
            FuzzyDiffExit::InputFileNotReadable,
            format!("Error: File not readable ({})", path.display()),
        ));
    }
    if !metadata.is_dir() && metadata.len() == 0 {
        return Err((
            FuzzyDiffExit::InputFileEmpty,
            format!("Error: File empty ({})", path.display()),
        ));
    }
    Ok(())
}

fn read_bounded(path: &Path) -> Result<Vec<u8>, String> {
    let file = File::open(path).map_err(|e| format!("cannot open '{}': {e}", path.display()))?;
    let mut bytes = Vec::new();
    file.take(MAX_INPUT_BYTES + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| format!("cannot read '{}': {e}", path.display()))?;
    if bytes.len() as u64 > MAX_INPUT_BYTES {
        return Err(format!(
            "input '{}' exceeds the comparison limit of {MAX_INPUT_BYTES} bytes",
            path.display()
        ));
    }
    Ok(bytes)
}

/// The text `FuzzyDiff -sort` writes to its temporary file: the first line, then
/// the remaining lines sorted bytewise, each followed by `\n`. Lines are split as
/// `std::getline` does, at `\n` only, so a `\r` stays part of its line.
pub fn sorted_lines(text: &[u8]) -> Vec<u8> {
    let mut lines: Vec<&[u8]> = text.split(|&b| b == b'\n').collect();
    if text.is_empty() || text.last() == Some(&b'\n') {
        lines.pop();
    }
    let mut out = Vec::with_capacity(text.len() + 1);
    let mut iter = lines.into_iter();
    out.extend_from_slice(iter.next().unwrap_or(&[]));
    out.push(b'\n');
    let mut rest: Vec<&[u8]> = iter.collect();
    rest.sort();
    for line in rest {
        out.extend_from_slice(line);
        out.push(b'\n');
    }
    out
}

/// A bounded reader for the subset of ParamXML that INI files use: `NODE`,
/// `ITEM`, `ITEMLIST` and `LISTITEM` elements with quoted attributes.
mod param_xml {
    /// Largest element nesting accepted.
    const MAX_DEPTH: usize = 64;
    /// Largest number of items accepted.
    const MAX_ITEMS: usize = 100_000;

    /// A parameter value: a scalar string or a list of strings.
    #[derive(Clone, Debug, PartialEq)]
    pub enum Value {
        /// The `value` attribute of an `ITEM`.
        Scalar(String),
        /// The `value` attributes of the `LISTITEM`s of an `ITEMLIST`.
        List(Vec<String>),
    }

    /// One `ITEM` or `ITEMLIST` with the names of its enclosing `NODE`s.
    #[derive(Clone, Debug, PartialEq)]
    pub struct Item {
        /// Names of the enclosing nodes, outermost first.
        pub path: Vec<String>,
        /// The item name.
        pub name: String,
        /// The item value.
        pub value: Value,
    }

    struct Tag {
        name: String,
        attributes: Vec<(String, String)>,
        closing: bool,
        self_closing: bool,
    }

    /// Parse the items of a ParamXML document.
    pub fn parse(bytes: &[u8]) -> Result<Vec<Item>, String> {
        if bytes.len() > super::MAX_INI_BYTES {
            return Err(format!("INI exceeds {} bytes", super::MAX_INI_BYTES));
        }
        let latin1 = declared_latin1(bytes);
        let mut items = Vec::new();
        let mut path: Vec<String> = Vec::new();
        let mut list: Option<(String, Vec<String>)> = None;
        let mut position = 0;
        while let Some(offset) = bytes
            .get(position..)
            .and_then(|rest| rest.iter().position(|&b| b == b'<'))
        {
            let start = position + offset;
            let rest = &bytes[start..];
            if rest.starts_with(b"<?") {
                position =
                    start + find(rest, b"?>").ok_or("unterminated processing instruction")? + 2;
                continue;
            }
            if rest.starts_with(b"<!--") {
                position = start + find(rest, b"-->").ok_or("unterminated comment")? + 3;
                continue;
            }
            if rest.starts_with(b"<!") {
                position = start + find(rest, b">").ok_or("unterminated declaration")? + 1;
                continue;
            }
            let (tag, length) = parse_tag(rest, latin1)?;
            position = start + length;
            let attribute = |key: &str| {
                tag.attributes
                    .iter()
                    .find(|(k, _)| k == key)
                    .map(|(_, v)| v.clone())
            };
            match (tag.name.as_str(), tag.closing) {
                ("NODE", false) => {
                    if !tag.self_closing {
                        if path.len() >= MAX_DEPTH {
                            return Err("INI nesting is too deep".into());
                        }
                        path.push(attribute("name").ok_or("NODE without name")?);
                    }
                }
                ("NODE", true) => {
                    path.pop().ok_or("unbalanced </NODE>")?;
                }
                ("ITEM", false) => {
                    let name = attribute("name").ok_or("ITEM without name")?;
                    let value = attribute("value").ok_or("ITEM without value")?;
                    items.push(Item {
                        path: path.clone(),
                        name,
                        value: Value::Scalar(value),
                    });
                }
                ("ITEMLIST", false) => {
                    let name = attribute("name").ok_or("ITEMLIST without name")?;
                    if tag.self_closing {
                        items.push(Item {
                            path: path.clone(),
                            name,
                            value: Value::List(Vec::new()),
                        });
                    } else {
                        list = Some((name, Vec::new()));
                    }
                }
                ("LISTITEM", false) => {
                    let value = attribute("value").ok_or("LISTITEM without value")?;
                    let (_, values) = list.as_mut().ok_or("LISTITEM outside ITEMLIST")?;
                    values.push(value);
                }
                ("ITEMLIST", true) => {
                    let (name, values) = list.take().ok_or("unbalanced </ITEMLIST>")?;
                    items.push(Item {
                        path: path.clone(),
                        name,
                        value: Value::List(values),
                    });
                }
                _ => {}
            }
            if items.len() > MAX_ITEMS {
                return Err("INI holds too many items".into());
            }
        }
        if !path.is_empty() || list.is_some() {
            return Err("INI ends inside an element".into());
        }
        Ok(items)
    }

    fn declared_latin1(bytes: &[u8]) -> bool {
        let head = &bytes[..bytes.len().min(200)];
        let lower: Vec<u8> = head.iter().map(u8::to_ascii_lowercase).collect();
        super::contains(&lower, b"encoding=\"iso-8859-1\"")
            || super::contains(&lower, b"encoding='iso-8859-1'")
    }

    fn find(haystack: &[u8], needle: &[u8]) -> Option<usize> {
        haystack.windows(needle.len()).position(|w| w == needle)
    }

    fn parse_tag(bytes: &[u8], latin1: bool) -> Result<(Tag, usize), String> {
        let mut i = 1;
        let closing = bytes.get(i) == Some(&b'/');
        if closing {
            i += 1;
        }
        let name_start = i;
        while bytes
            .get(i)
            .is_some_and(|b| !super::is_c_space(*b) && *b != b'>' && *b != b'/')
        {
            i += 1;
        }
        let name = String::from_utf8_lossy(&bytes[name_start..i]).into_owned();
        let mut attributes = Vec::new();
        loop {
            while bytes.get(i).is_some_and(|b| super::is_c_space(*b)) {
                i += 1;
            }
            match bytes.get(i) {
                None => return Err("unterminated tag".into()),
                Some(b'>') => {
                    return Ok((
                        Tag {
                            name,
                            attributes,
                            closing,
                            self_closing: false,
                        },
                        i + 1,
                    ));
                }
                Some(b'/') if bytes.get(i + 1) == Some(&b'>') => {
                    return Ok((
                        Tag {
                            name,
                            attributes,
                            closing,
                            self_closing: true,
                        },
                        i + 2,
                    ));
                }
                Some(_) => {}
            }
            let key_start = i;
            while bytes
                .get(i)
                .is_some_and(|b| *b != b'=' && !super::is_c_space(*b) && *b != b'>')
            {
                i += 1;
            }
            let key = String::from_utf8_lossy(&bytes[key_start..i]).into_owned();
            while bytes.get(i).is_some_and(|b| super::is_c_space(*b)) {
                i += 1;
            }
            if bytes.get(i) != Some(&b'=') {
                return Err(format!("attribute '{key}' without value"));
            }
            i += 1;
            while bytes.get(i).is_some_and(|b| super::is_c_space(*b)) {
                i += 1;
            }
            let quote = *bytes.get(i).ok_or("unterminated attribute")?;
            if quote != b'"' && quote != b'\'' {
                return Err(format!("attribute '{key}' is not quoted"));
            }
            let value_start = i + 1;
            let length = bytes
                .get(value_start..)
                .and_then(|rest| rest.iter().position(|&b| b == quote))
                .ok_or("unterminated attribute value")?;
            let raw = &bytes[value_start..value_start + length];
            attributes.push((key, decode(raw, latin1)?));
            i = value_start + length + 1;
        }
    }

    fn decode(raw: &[u8], latin1: bool) -> Result<String, String> {
        let text: String = if latin1 {
            raw.iter().map(|&b| char::from(b)).collect()
        } else {
            String::from_utf8(raw.to_vec()).map_err(|_| "attribute is not UTF-8".to_owned())?
        };
        let mut out = String::with_capacity(text.len());
        let mut rest = text.as_str();
        while let Some(amp) = rest.find('&') {
            out.push_str(&rest[..amp]);
            let after = &rest[amp + 1..];
            let semicolon = after.find(';').ok_or("unterminated entity")?;
            let entity = &after[..semicolon];
            let decoded = match entity {
                "lt" => '<',
                "gt" => '>',
                "amp" => '&',
                "quot" => '"',
                "apos" => '\'',
                _ => {
                    let code = if let Some(hex) = entity.strip_prefix("#x") {
                        u32::from_str_radix(hex, 16).ok()
                    } else if let Some(dec) = entity.strip_prefix('#') {
                        dec.parse::<u32>().ok()
                    } else {
                        None
                    };
                    code.and_then(char::from_u32)
                        .ok_or_else(|| format!("unknown entity '&{entity};'"))?
                }
            };
            out.push(decoded);
            rest = &after[semicolon + 1..];
        }
        out.push_str(rest);
        Ok(out)
    }
}
