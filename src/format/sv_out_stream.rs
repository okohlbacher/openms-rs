// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Separated-value output with automatic separators and string quoting.
//!
//! Ports `FORMAT/SVOutStream.h` and `SVOutStream.cpp`. See
//! `docs/SV_OUT_STREAM_SUPPORT.md` for the full API mapping.
//!
//! The source class derives from `std::ostream` and overloads `operator<<` so
//! that a separator appears between items but not at the start of a line. Rust
//! has no `operator<<`, so the overload set becomes named methods on
//! [`SVOutStream`](crate::format::sv_out_stream::SVOutStream): field writers
//! ([`write_field`](crate::format::sv_out_stream::SVOutStream::write_field),
//! [`write_char`](crate::format::sv_out_stream::SVOutStream::write_char),
//! [`write_number`](crate::format::sv_out_stream::SVOutStream::write_number),
//! [`write_display`](crate::format::sv_out_stream::SVOutStream::write_display))
//! and line terminators
//! ([`newline`](crate::format::sv_out_stream::SVOutStream::newline),
//! [`end_line`](crate::format::sv_out_stream::SVOutStream::end_line)).
//!
//! Numeric text is produced by
//! [`SvNumber::sv_text`](crate::format::sv_out_stream::SvNumber::sv_text),
//! which reproduces `StringUtils::toStr` / `NumericFormatting::appendNumeric`
//! at full precision. A general `StringUtils` port is separate work; only the
//! conversion this stream depends on lives here.

use crate::{Error, Result};
use std::fmt::Display;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// How embedded quote characters are handled when a field is quoted.
///
/// Ports `OpenMS::QuotingMethod` from `DATASTRUCTURES/StringUtils.h`, whose
/// three enumerators are the source's own `NONE`, `ESCAPE` and `DOUBLE`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum QuotingMethod {
    /// Do not quote. [`SVOutStream`] then replaces every occurrence of the
    /// separator inside a field with the configured replacement instead.
    None,
    /// Wrap in quote characters, prefixing `\` to backslashes and to the quote
    /// character. The source applies the backslash substitution first, so
    /// `a\"b` becomes `"a\\\"b"`.
    Escape,
    /// Wrap in quote characters and double every embedded quote character, the
    /// source default for this stream.
    #[default]
    Double,
}

/// Whether a numeric value is finite, not-a-number or an infinity.
///
/// [`SVOutStream::write_value_or_nan`] needs this distinction. The source
/// obtains it from `boost::math::isfinite` / `isnan` plus a `thing < 0` test,
/// the comment in `SVOutStream.h` noting the Boost detour exists because
/// `isfinite` was unavailable on macOS.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum NumberClass {
    /// A finite value, including zero and negative zero.
    Finite,
    /// Not-a-number, of either sign.
    NotANumber,
    /// An infinity; `negative` mirrors the source's `thing < 0` branch.
    Infinite {
        /// True for negative infinity.
        negative: bool,
    },
}

/// Administrative ceilings for one stream, checked before anything is written.
///
/// The source has no limits: a field of any length is quoted into a fresh
/// `std::string` and pushed at the stream. These bounds exist because
/// `SVOutStream` renders caller data and quoting can double a field's length.
/// Every ceiling is checked before the rendering allocation, so exceeding one
/// leaves the output byte-for-byte unchanged.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SVLimits {
    /// Maximum bytes of one rendered field or raw chunk, quoting included.
    pub max_field_bytes: usize,
    /// Maximum fields on one line, reset by every line terminator.
    pub max_row_fields: usize,
    /// Maximum line terminators over the life of the stream.
    pub max_rows: usize,
}

impl SVLimits {
    /// Default rendered size of one field: 16 MiB.
    pub const MAX_FIELD_BYTES: usize = 16 * 1024 * 1024;
    /// Default number of fields on one line.
    pub const MAX_ROW_FIELDS: usize = 1 << 20;
    /// Default number of lines in one file.
    pub const MAX_ROWS: usize = 1 << 30;
}

impl Default for SVLimits {
    fn default() -> Self {
        Self {
            max_field_bytes: Self::MAX_FIELD_BYTES,
            max_row_fields: Self::MAX_ROW_FIELDS,
            max_rows: Self::MAX_ROWS,
        }
    }
}

/// Numeric text exactly as the source's `StringUtils::toStr` renders it.
///
/// The source `SVOutStream` has two `operator<<` templates selected by
/// `std::is_arithmetic`. The arithmetic one converts through
/// `StringUtils::toStr` rather than through stream formatting, which is what
/// makes its output independent of the stream's locale and precision; this
/// trait is that conversion. There are no defaulted methods, so a new
/// implementation cannot silently inherit a wrong rendering.
pub trait SvNumber: Copy {
    /// The value's source text. Never contains a newline or the empty string.
    fn sv_text(self) -> String;
    /// Finiteness, as [`SVOutStream::write_value_or_nan`] classifies it.
    fn sv_class(self) -> NumberClass;
}

/// Digits after the decimal point for `f64`, `std::numeric_limits<double>::digits10`.
pub const F64_FIXED_DIGITS: usize = 15;
/// Digits after the decimal point for `f32`, `std::numeric_limits<float>::digits10`.
pub const F32_FIXED_DIGITS: usize = 6;

/// Below this magnitude a nonzero value is written in scientific notation.
pub const SCIENTIFIC_LOWER: f64 = 1e-2;
/// From this magnitude upwards a value is written in scientific notation.
pub const SCIENTIFIC_UPPER: f64 = 1e4;
/// [`SCIENTIFIC_LOWER`] as the source's `T(1e-2)` for `T = float`.
///
/// Not the same comparison: `0.01_f32` is below `1e-2_f64` but equal to
/// `1e-2_f32`, so it is scientific in the `f64` pipeline and fixed in this one.
const SCIENTIFIC_LOWER_F32: f32 = 1e-2;
/// [`SCIENTIFIC_UPPER`] as the source's `T(1e4)` for `T = float`.
const SCIENTIFIC_UPPER_F32: f32 = 1e4;
/// Largest `fixed_digits` the two text renderers honour.
///
/// Native: the source has no clamp and overflows its `char buf[64]` instead,
/// falling through to a `std::to_string(double)` with six fractional digits.
pub const MAX_FIXED_DIGITS: usize = 512;

/// Render a floating-point value as `NumericFormatting::appendNumeric` does at
/// full precision.
///
/// `NaN` becomes `NaN` and an infinity becomes `inf` or `-inf`; the header
/// comment records that the uppercase spelling is kept for compatibility with
/// the removed Karma formatter, and that writing `inf.0` once produced tokens
/// no reader could parse. A nonzero magnitude at or above
/// [`SCIENTIFIC_UPPER`] or below [`SCIENTIFIC_LOWER`] is written in scientific
/// notation using the shortest representation that round-trips; everything
/// else, zero included, is written with `fixed_digits` digits after the
/// decimal point. Trailing fractional zeros are then trimmed, but at least one
/// digit after the point survives, so `5.0` never degrades to `5`. In
/// scientific notation the `+` is dropped from the exponent while its
/// zero padding to two digits is kept, and a mantissa without a decimal point
/// gains `.0`: `1e4` is written `1.0e04`.
///
/// # Arguments
///
/// * `value` - any finite or nonfinite value; no range is rejected.
/// * `fixed_digits` - digits after the decimal point in fixed notation,
///   [`F64_FIXED_DIGITS`] in the source's `double` call. Values above
///   [`MAX_FIXED_DIGITS`] are clamped; the source has no clamp and instead
///   overflows its `char buf[64]`, after which `appendNumeric` falls through to
///   `std::to_string(static_cast<double>(value))`, i.e. exactly six fractional
///   digits (`NumericFormatting.h:136-139`). Neither behaviour is reachable
///   from the source's own call sites, which pass at most
///   `writtenDigits<long double>()` = 18.
///
/// This is the `T = double` instantiation of the source template. An `f32`
/// must go through [`source_f32_text`] instead: promoting it here would change
/// both the shortest round-trip mantissa and the branch boundary.
pub fn source_float_text(value: f64, fixed_digits: usize) -> String {
    if value.is_nan() {
        return "NaN".to_owned();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.to_owned();
    }
    let magnitude = value.abs();
    if magnitude != 0.0 && !(SCIENTIFIC_LOWER..SCIENTIFIC_UPPER).contains(&magnitude) {
        return scientific_text(format!("{value:e}"));
    }
    let digits = fixed_digits.min(MAX_FIXED_DIGITS);
    trim_fraction(format!("{value:.digits$}"))
}

/// [`source_float_text`] for an `f32`: the `T = float` instantiation.
///
/// `appendNumeric` is a template, so for a `float` argument the source compares
/// `abs_val` against `T(1e-2)` and `T(1e4)` in *float* arithmetic and calls
/// `std::to_chars` with the `float` overload, whose shortest round-trip is the
/// shortest decimal that round-trips as a `float`. Both differ from the `f64`
/// pipeline: `1.23e-5_f32` is `1.23e-05` here and `1.2299999980314169e-05`
/// after a promotion to `f64`, and `0.01_f32` — which equals `1e-2_f32` exactly
/// and so takes the fixed branch — is `0.01` here and `9.999999776482582e-03`
/// after a promotion.
///
/// # Arguments
///
/// * `value` - any finite or nonfinite value; no range is rejected.
/// * `fixed_digits` - digits after the decimal point in fixed notation,
///   [`F32_FIXED_DIGITS`] in the source's `float` call, clamped as above.
pub fn source_f32_text(value: f32, fixed_digits: usize) -> String {
    if value.is_nan() {
        return "NaN".to_owned();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.to_owned();
    }
    let magnitude = value.abs();
    if magnitude != 0.0 && !(SCIENTIFIC_LOWER_F32..SCIENTIFIC_UPPER_F32).contains(&magnitude) {
        return scientific_text(format!("{value:e}"));
    }
    let digits = fixed_digits.min(MAX_FIXED_DIGITS);
    trim_fraction(format!("{value:.digits$}"))
}

/// Trim trailing fractional zeros, keeping one digit after the decimal point.
///
/// A rendering without a decimal point gains `.0`, which is the source's
/// `is_floating_point_v` branch. Only text this crate produced with `{:.N}` is
/// inspected, so the search for `.` and the trailing `0` run never leaves a
/// character boundary.
fn trim_fraction(mut text: String) -> String {
    match text.find('.') {
        None => {
            text.push_str(".0");
            text
        }
        Some(dot) => {
            let mut end = text.len();
            while end > dot + 2 && text.as_bytes().get(end - 1) == Some(&b'0') {
                end -= 1;
            }
            text.truncate(end);
            text
        }
    }
}

/// Reshape Rust's `{:e}` into the source's scientific spelling.
///
/// Rust writes `-1.23e45`; `std::to_chars` writes `-1.23e+45` and the source
/// then removes the `+` but keeps the two-digit zero padding, so both agree on
/// the mantissa and differ only in exponent padding, which this restores.
fn scientific_text(text: String) -> String {
    let Some((mantissa, exponent)) = text.split_once('e') else {
        // Rust's LowerExp always emits 'e'; this branch cannot be reached.
        return text;
    };
    let mut out = String::with_capacity(text.len() + 4);
    out.push_str(mantissa);
    if !mantissa.contains('.') {
        out.push_str(".0");
    }
    out.push('e');
    let digits = match exponent.strip_prefix('-') {
        Some(rest) => {
            out.push('-');
            rest
        }
        None => exponent,
    };
    if digits.len() < 2 {
        out.push('0');
    }
    out.push_str(digits);
    out
}

macro_rules! integer_sv_number {
    ($($type:ty),*) => {
        $(impl SvNumber for $type {
            fn sv_text(self) -> String {
                self.to_string()
            }
            fn sv_class(self) -> NumberClass {
                NumberClass::Finite
            }
        })*
    };
}
integer_sv_number!(
    i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize
);

impl SvNumber for f64 {
    fn sv_text(self) -> String {
        source_float_text(self, F64_FIXED_DIGITS)
    }
    fn sv_class(self) -> NumberClass {
        float_class(self.is_nan(), self.is_infinite(), self < 0.0)
    }
}

impl SvNumber for f32 {
    fn sv_text(self) -> String {
        source_f32_text(self, F32_FIXED_DIGITS)
    }
    fn sv_class(self) -> NumberClass {
        float_class(self.is_nan(), self.is_infinite(), self < 0.0)
    }
}

fn float_class(nan: bool, infinite: bool, negative: bool) -> NumberClass {
    if nan {
        NumberClass::NotANumber
    } else if infinite {
        NumberClass::Infinite { negative }
    } else {
        NumberClass::Finite
    }
}

fn limit(what: &str) -> Error {
    Error::InvalidValue(format!("separated-value {what} limit exceeded"))
}

/// A separated-value output stream with automatic separators and quoting.
///
/// Ports `SVOutStream`. Every field writer emits the separator first unless
/// the stream is at the start of a line, and each line terminator puts it back
/// there; [`at_line_start`](Self::at_line_start) exposes that one bit of state,
/// which the source keeps private as `newline_`.
///
/// The source requires `nl` or `std::endl` as the line delimiter and its class
/// documentation states that a literal `"\n"` "won't be accepted": a `"\n"`
/// reaches `operator<<(std::string)`, which throws
/// `Exception::IllegalArgument`. This port keeps that rule as
/// [`Error::InvalidValue`] from every field writer, so the separator
/// bookkeeping cannot desynchronise from the emitted text.
///
/// Unlike the source this type does not derive from an output stream, so there
/// is no way to bypass the separator bookkeeping by casting to the base class.
/// [`write_raw`](Self::write_raw) is the one controlled escape hatch, as
/// `SVOutStream::write` is upstream.
///
/// # Examples
///
/// ```
/// use openms::format::sv_out_stream::{QuotingMethod, SVOutStream};
///
/// let mut out = SVOutStream::with_options(Vec::new(), ",", "_", QuotingMethod::None)?;
/// out.write_field("a")?;
/// out.write_field("d,f")?;
/// out.newline()?;
/// out.write_number(123)?;
/// out.write_number(3.14)?;
/// out.newline()?;
/// assert_eq!(String::from_utf8(out.finish()?).unwrap(), "a,d_f\n123,3.14\n");
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Debug)]
pub struct SVOutStream<W: Write> {
    writer: W,
    separator: String,
    replacement: String,
    nan: String,
    infinity: String,
    quoting: QuotingMethod,
    modify_strings: bool,
    newline: bool,
    limits: SVLimits,
    row_fields: usize,
    rows: usize,
}

impl<W: Write> SVOutStream<W> {
    /// The source default separator, a tab.
    pub const DEFAULT_SEPARATOR: &str = "\t";
    /// The source default replacement for a separator inside an unquoted field.
    pub const DEFAULT_REPLACEMENT: &str = "_";
    /// The source's `nan_`, used by [`Self::write_value_or_nan`].
    pub const DEFAULT_NAN: &str = "nan";
    /// The source's `inf_`, used by [`Self::write_value_or_nan`].
    pub const DEFAULT_INFINITY: &str = "inf";

    /// A stream over `writer` with the source constructor defaults: a tab
    /// separator, `_` as the replacement and [`QuotingMethod::Double`].
    ///
    /// This is the source's `SVOutStream(std::ostream&)` with all three
    /// defaults taken. The source also calls
    /// `precision(std::numeric_limits<double>::digits10)` on itself; that
    /// setting only reached its generic non-arithmetic overload, and this port
    /// has no stream state to configure - see [`Self::write_display`].
    pub fn new(writer: W) -> Self {
        Self {
            writer,
            separator: Self::DEFAULT_SEPARATOR.to_owned(),
            replacement: Self::DEFAULT_REPLACEMENT.to_owned(),
            nan: Self::DEFAULT_NAN.to_owned(),
            infinity: Self::DEFAULT_INFINITY.to_owned(),
            quoting: QuotingMethod::Double,
            modify_strings: true,
            newline: true,
            limits: SVLimits::default(),
            row_fields: 0,
            rows: 0,
        }
    }

    /// A stream over `writer` with an explicit separator, replacement and
    /// quoting method.
    ///
    /// # Arguments
    ///
    /// * `separator` - separator string, typically a comma, semicolon or tab.
    /// * `replacement` - substituted for occurrences of `separator` inside a
    ///   field when `quoting` is [`QuotingMethod::None`]; unused otherwise.
    /// * `quoting` - quoting method for fields, as `String::quote`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `separator` is empty or when
    /// `separator` or `replacement` contains a newline. The source checks
    /// none of the three: an empty separator makes `StringUtils::substitute`
    /// return the field unchanged and emits nothing between fields, so the
    /// file cannot be read back, and a newline in either string moves the
    /// cursor to a new line while the stream still believes it is mid-line.
    pub fn with_options(
        writer: W,
        separator: &str,
        replacement: &str,
        quoting: QuotingMethod,
    ) -> Result<Self> {
        if separator.is_empty() {
            return Err(Error::InvalidValue(
                "separated-value separator must not be empty".into(),
            ));
        }
        for (text, what) in [(separator, "separator"), (replacement, "replacement")] {
            if text.contains('\n') || text.contains('\r') {
                return Err(Error::InvalidValue(format!(
                    "separated-value {what} must not contain a line break"
                )));
            }
        }
        let mut stream = Self::new(writer);
        stream.separator = separator.to_owned();
        stream.replacement = replacement.to_owned();
        stream.quoting = quoting;
        Ok(stream)
    }

    /// Replace the administrative ceilings of [`SVLimits`].
    pub fn with_limits(mut self, limits: SVLimits) -> Self {
        self.limits = limits;
        self
    }

    /// The configured separator string, the source's protected `sep_`.
    pub fn separator(&self) -> &str {
        &self.separator
    }

    /// The configured separator replacement, the source's `replacement_`.
    pub fn replacement(&self) -> &str {
        &self.replacement
    }

    /// The configured quoting method, the source's `quoting_`.
    pub fn quoting(&self) -> QuotingMethod {
        self.quoting
    }

    /// Text written for not-a-number, the source's `nan_`.
    pub fn nan_text(&self) -> &str {
        &self.nan
    }

    /// Text written for an infinity, the source's `inf_`. A negative infinity
    /// is this text with `-` prefixed, as the source builds `"-" + inf_`.
    pub fn infinity_text(&self) -> &str {
        &self.infinity
    }

    /// The administrative ceilings in force.
    pub fn limits(&self) -> SVLimits {
        self.limits
    }

    /// Whether string modification is on, the source's `modify_strings_`.
    pub fn modify_strings(&self) -> bool {
        self.modify_strings
    }

    /// Whether the next field is written without a leading separator.
    ///
    /// The source's `newline_`, true after construction and after every line
    /// terminator. It is protected upstream; this port exposes it because a
    /// caller mixing [`Self::write_raw`] with field writers has to be able to
    /// see it.
    pub fn at_line_start(&self) -> bool {
        self.newline
    }

    /// Fields written on the current line, reset by every line terminator.
    ///
    /// Native accounting with no source counterpart; the ceiling it is checked
    /// against is [`SVLimits::max_row_fields`].
    pub fn row_fields(&self) -> usize {
        self.row_fields
    }

    /// Line terminators written so far, checked against [`SVLimits::max_rows`].
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Switch modification of fields - quoting, or replacing separators - on
    /// or off, returning the previous state.
    ///
    /// Ports `modifyStrings`. With modification off a field is written
    /// verbatim; separators are still inserted between fields, and the
    /// newline rejection still applies. The source uses this internally to
    /// write its `nan_`/`inf_` texts unquoted.
    pub fn set_modify_strings(&mut self, modify: bool) -> bool {
        let previous = self.modify_strings;
        self.modify_strings = modify;
        previous
    }

    /// Write one field, quoting or substituting it as configured.
    ///
    /// Ports `operator<<(std::string)`, and with it the `const char*` overload,
    /// which the source implements as `operator<<(std::string(c_str))`. A
    /// separator is emitted first unless the stream is at the start of a line.
    /// With modification on, [`QuotingMethod::None`] replaces every occurrence
    /// of the separator with the replacement and any other method quotes the
    /// field with `"`; with modification off the field is written verbatim.
    ///
    /// The field and its separator are rendered into one buffer and written
    /// once, so a rejected field emits no bytes and does not advance the
    /// line state.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `text` contains a newline - the source
    /// throws `Exception::IllegalArgument` with the message "argument must not
    /// contain newline characters" - or when a ceiling in [`SVLimits`] would
    /// be exceeded. [`Error::Io`] propagates a failed write.
    ///
    /// A carriage return is rejected as well. The source tests only for `\n`,
    /// so a lone `\r` reaches the file and, on a reader treating it as a line
    /// break, splits a row the writer counted as one.
    pub fn write_field(&mut self, text: &str) -> Result<()> {
        if text.contains('\n') || text.contains('\r') {
            return Err(Error::InvalidValue(
                "separated-value field must not contain newline characters".into(),
            ));
        }
        self.preflight(text.len())?;
        let rendered = if !self.modify_strings {
            text.to_owned()
        } else if self.quoting != QuotingMethod::None {
            quote(text, self.quoting)
        } else {
            self.preflight_substitution(text)?;
            text.replace(self.separator.as_str(), &self.replacement)
        };
        self.emit_field(&rendered)
    }

    /// Write one character as a field.
    ///
    /// Ports `operator<<(const char c)`, which the source routes through
    /// `StringUtils::toStr(c)` into the string overload, so the character is
    /// quoted like any other field. A `char` here is a Unicode scalar value
    /// rather than the source's single byte; a `'\n'` is rejected exactly as
    /// the one-character string would be.
    ///
    /// # Errors
    ///
    /// As [`Self::write_field`].
    pub fn write_char(&mut self, value: char) -> Result<()> {
        self.write_field(&value.to_string())
    }

    /// Write a numeric field using the source's own numeric conversion.
    ///
    /// Ports the arithmetic `operator<<` template, whose body converts through
    /// `StringUtils::toStr` and therefore never quotes and never consults the
    /// stream's locale or precision. A separator is emitted first unless the
    /// stream is at the start of a line.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a ceiling in [`SVLimits`] would be
    /// exceeded; [`Error::Io`] propagates a failed write.
    pub fn write_number<T: SvNumber>(&mut self, value: T) -> Result<()> {
        let text = value.sv_text();
        self.preflight(text.len())?;
        self.emit_field(&text)
    }

    /// Write a numeric field, substituting the not-a-number and infinity texts.
    ///
    /// Ports `writeValueOrNan`. A finite value is written by
    /// [`Self::write_number`]. Otherwise modification is switched off, the
    /// stream writes [`Self::nan_text`], [`Self::infinity_text`] or that text
    /// with `-` prefixed, and the previous modification state is restored -
    /// including when the write fails, which the source's straight-line code
    /// cannot guarantee because it has no failure path.
    ///
    /// The source header notes this "would not be needed for Linux": the
    /// platform difference it works around is in the C++ stream formatting of
    /// nonfinite values, not in this port.
    ///
    /// # Errors
    ///
    /// As [`Self::write_number`].
    pub fn write_value_or_nan<T: SvNumber>(&mut self, value: T) -> Result<()> {
        let text = match value.sv_class() {
            NumberClass::Finite => return self.write_number(value),
            NumberClass::NotANumber => self.nan.clone(),
            NumberClass::Infinite { negative: false } => self.infinity.clone(),
            NumberClass::Infinite { negative: true } => format!("-{}", self.infinity),
        };
        let previous = self.set_modify_strings(false);
        let result = self.write_field(&text);
        self.set_modify_strings(previous);
        result
    }

    /// Write a field from any [`Display`] value, without quoting.
    ///
    /// Ports the generic non-arithmetic `operator<<` template, which inserts
    /// the separator and then streams the value into the underlying
    /// `std::ostream` with no quoting or separator substitution.
    ///
    /// The source's constructor sets the stream precision to
    /// `std::numeric_limits<double>::digits10`, which affected exactly this
    /// overload when it was handed a floating-point value through a
    /// user-defined type. Rust's [`Display`] carries no such setting, so a
    /// value rendering a float here is spelled by Rust rather than by the
    /// source's `%g`-style conversion; [`Self::write_number`] is the
    /// source-equivalent route for numbers.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the rendered text contains a newline or a
    /// carriage return, or when a ceiling in [`SVLimits`] would be exceeded;
    /// [`Error::Io`] propagates a failed write. The source performs no newline
    /// check on this overload at all, so a `Display` implementation emitting
    /// one would silently corrupt the row structure.
    pub fn write_display<T: Display>(&mut self, value: T) -> Result<()> {
        let text = value.to_string();
        if text.contains('\n') || text.contains('\r') {
            return Err(Error::InvalidValue(
                "separated-value field must not contain newline characters".into(),
            ));
        }
        self.preflight(text.len())?;
        self.emit_field(&text)
    }

    /// Write text verbatim: no separator, no quoting, no state change.
    ///
    /// Ports `write(const std::string&)`, whose header comment says
    /// "no quoting: useful for comments, but use only on a line of its own!".
    /// The reason for that warning is that the source's `write` does not touch
    /// `newline_`: after writing a comment ending in `\n` while mid-line, the
    /// next field still receives a leading separator. This port reproduces
    /// that, so [`at_line_start`](Self::at_line_start) is unchanged here and a
    /// caller who needs the flag reset calls [`Self::newline`] first.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the chunk is longer than
    /// [`SVLimits::max_field_bytes`]; [`Error::Io`] propagates a failed write.
    /// The text itself is not inspected, matching the source.
    pub fn write_raw(&mut self, text: &str) -> Result<()> {
        if text.len() > self.limits.max_field_bytes {
            return Err(limit("raw chunk"));
        }
        self.writer.write_all(text.as_bytes())?;
        Ok(())
    }

    /// End the line without flushing.
    ///
    /// Ports `operator<<(enum Newline)`, the source's preferred `nl`
    /// delimiter: it sets the line-start flag and writes `"\n"`. The source
    /// header recommends it over `std::endl` "for improved performance",
    /// which is exactly the missing flush.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when [`SVLimits::max_rows`] would be exceeded;
    /// [`Error::Io`] propagates a failed write.
    pub fn newline(&mut self) -> Result<()> {
        let rows = self.rows.checked_add(1).ok_or_else(|| limit("row"))?;
        if rows > self.limits.max_rows {
            return Err(limit("row"));
        }
        self.writer.write_all(b"\n")?;
        self.rows = rows;
        self.row_fields = 0;
        self.newline = true;
        Ok(())
    }

    /// End the line and flush.
    ///
    /// Ports `operator<<(std::ostream& (*fp)(std::ostream&))`, the manipulator
    /// overload that exists to catch `std::endl`. The source detects it by
    /// applying the function pointer to a scratch `std::stringstream` and
    /// testing whether `"\n"` came out, because comparing the pointer against
    /// `&std::endl` does not work under libc++; any other manipulator is
    /// forwarded without changing the line state. Rust has no manipulators, so
    /// the one case that mattered becomes this method.
    ///
    /// # Errors
    ///
    /// As [`Self::newline`], plus a failed flush as [`Error::Io`].
    pub fn end_line(&mut self) -> Result<()> {
        self.newline()?;
        self.writer.flush()?;
        Ok(())
    }

    /// Flush the underlying writer.
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()?;
        Ok(())
    }

    /// Flush and return the underlying writer.
    ///
    /// The source destructor closes the `std::ofstream` it owns when the
    /// filename constructor was used, and does nothing for a borrowed stream.
    /// Rust drops the writer either way, but a drop cannot report a failed
    /// flush, so this is the checked end of the stream.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the final flush fails.
    pub fn finish(mut self) -> Result<W> {
        self.writer.flush()?;
        Ok(self.writer)
    }

    /// Charge a [`QuotingMethod::None`] substitution before it is allocated.
    ///
    /// [`Self::preflight`]'s worst case bounds both `quote` paths, because
    /// each at most doubles the field. It does not bound the substitution:
    /// [`Self::replacement`] may be arbitrarily longer than the
    /// separator, so a field of separators grows by
    /// `replacement.len() / separator.len()`, which is not 2. The rendered
    /// length is therefore computed exactly here — by counting separators,
    /// which allocates nothing — and charged before `str::replace` runs.
    /// [`Self::emit_field`] checks the same ceiling again on the rendered text.
    fn preflight_substitution(&self, text: &str) -> Result<()> {
        let separator = self.separator.as_str();
        // The constructor refuses an empty separator, so `matches` counts
        // non-overlapping occurrences exactly as `replace` substitutes them.
        let occurrences = text.matches(separator).count();
        let removed = occurrences
            .checked_mul(separator.len())
            .ok_or_else(|| limit("field"))?;
        let substituted = occurrences
            .checked_mul(self.replacement.len())
            .and_then(|grown| grown.checked_add(text.len().saturating_sub(removed)))
            .ok_or_else(|| limit("field"))?;
        let total = substituted
            .checked_add(separator.len())
            .ok_or_else(|| limit("field"))?;
        if total > self.limits.max_field_bytes {
            return Err(limit("field"));
        }
        Ok(())
    }

    /// Charge one field against the ceilings before it is rendered.
    ///
    /// Quoting can at most double a field and add two quote characters, so the
    /// worst case is checked before the rendering allocation happens. The
    /// substitution path is charged separately, by
    /// [`Self::preflight_substitution`].
    fn preflight(&self, raw: usize) -> Result<()> {
        let fields = self
            .row_fields
            .checked_add(1)
            .ok_or_else(|| limit("field"))?;
        if fields > self.limits.max_row_fields {
            return Err(limit("field"));
        }
        let worst = raw
            .checked_mul(2)
            .and_then(|n| n.checked_add(2))
            .and_then(|n| n.checked_add(self.separator.len()))
            .ok_or_else(|| limit("field"))?;
        if worst > self.limits.max_field_bytes {
            return Err(limit("field"));
        }
        Ok(())
    }

    /// Write a rendered field and its separator as one chunk.
    fn emit_field(&mut self, rendered: &str) -> Result<()> {
        if rendered.len() > self.limits.max_field_bytes {
            return Err(limit("field"));
        }
        let fields = self
            .row_fields
            .checked_add(1)
            .ok_or_else(|| limit("field"))?;
        let mut chunk = String::with_capacity(rendered.len().saturating_add(self.separator.len()));
        if !self.newline {
            chunk.push_str(&self.separator);
        }
        chunk.push_str(rendered);
        self.writer.write_all(chunk.as_bytes())?;
        self.row_fields = fields;
        self.newline = false;
        Ok(())
    }
}

impl SVOutStream<BufWriter<File>> {
    /// Create or truncate `path` and write to it, with the source defaults.
    ///
    /// Ports `SVOutStream(const std::string& file_out, ...)`, which opens an
    /// owned `std::ofstream` on the path and overwrites an existing file.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be opened for writing; the source
    /// throws `Exception::FileNotWritable` naming the path.
    pub fn create(path: impl AsRef<Path>) -> Result<Self> {
        Ok(Self::new(BufWriter::new(File::create(path.as_ref())?)))
    }

    /// Create or truncate `path` with an explicit separator, replacement and
    /// quoting method.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be opened for writing, and the
    /// [`Error::InvalidValue`] cases of [`SVOutStream::with_options`].
    pub fn create_with_options(
        path: impl AsRef<Path>,
        separator: &str,
        replacement: &str,
        quoting: QuotingMethod,
    ) -> Result<Self> {
        let file = BufWriter::new(File::create(path.as_ref())?);
        Self::with_options(file, separator, replacement, quoting)
    }
}

/// Quote `text` with `"` as `StringUtils::quote` does.
///
/// [`QuotingMethod::Escape`] prefixes `\` to every backslash and then to every
/// quote character - that order is the source's and is not interchangeable.
/// [`QuotingMethod::Double`] doubles every quote character.
/// [`QuotingMethod::None`] still adds the surrounding quotes, because
/// `StringUtils::quote` skips only the substitution; `SVOutStream` never
/// reaches this case, since it substitutes separators instead.
pub fn quote(text: &str, method: QuotingMethod) -> String {
    let body = match method {
        QuotingMethod::None => text.to_owned(),
        QuotingMethod::Escape => text.replace('\\', "\\\\").replace('"', "\\\""),
        QuotingMethod::Double => text.replace('"', "\"\""),
    };
    let mut out = String::with_capacity(body.len() + 2);
    out.push('"');
    out.push_str(&body);
    out.push('"');
    out
}

// Every behaviour of this module is tested in tests/sv_out_stream.rs, which
// ports all twelve sections of SVOutStream_test.cpp.
