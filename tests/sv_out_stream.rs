// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports every `START_SECTION` of `SVOutStream_test.cpp` (pinned revision
//! bc9cc12, 12 sections) plus native boundary tests. The expected strings are
//! transcribed C++ literals, tier-3 evidence; provenance and section mapping
//! are in `tests/data/sv_out_stream_provenance.json` and
//! `docs/SV_OUT_STREAM_SUPPORT.md`.

// 3.14 is a value transcribed from SVOutStream_test.cpp, not an attempt at pi.
#![allow(clippy::approx_constant)]

use openms::Error;
use openms::format::sv_out_stream::{
    F32_FIXED_DIGITS, F64_FIXED_DIGITS, MAX_FIXED_DIGITS, NumberClass, QuotingMethod, SVLimits,
    SVOutStream, SvNumber, quote, source_f32_text, source_float_text,
};

/// Build a stream over an in-memory buffer, the test's `stringstream`.
fn stream(separator: &str, replacement: &str, quoting: QuotingMethod) -> SVOutStream<Vec<u8>> {
    SVOutStream::with_options(Vec::new(), separator, replacement, quoting)
        .expect("valid separator and replacement")
}

fn text(out: SVOutStream<Vec<u8>>) -> String {
    String::from_utf8(out.finish().expect("flush")).expect("UTF-8 output")
}

// ---------------------------------------------------------------------------
// START_SECTION((SVOutStream(std::ostream& out, const std::string& sep="\t",
//   const std::string& replacement="_",
//   OpenMS::QuotingMethod quoting=OpenMS::QuotingMethod::DOUBLE)))
// The upstream section only asserts the pointer is non-null; the constructor
// defaults it exercises are asserted here instead.
// ---------------------------------------------------------------------------
#[test]
fn stream_constructor_takes_the_source_defaults() {
    let out = SVOutStream::new(Vec::new());
    assert_eq!(out.separator(), "\t");
    assert_eq!(out.replacement(), "_");
    assert_eq!(out.quoting(), QuotingMethod::Double);
    assert_eq!(out.nan_text(), "nan");
    assert_eq!(out.infinity_text(), "inf");
    assert!(out.modify_strings());
    assert!(out.at_line_start());
    assert_eq!(out.row_fields(), 0);
    assert_eq!(out.rows(), 0);
    assert_eq!(out.limits(), SVLimits::default());
    // The explicit form reaches the same state with other values.
    let out = stream(",", "-", QuotingMethod::Escape);
    assert_eq!(out.separator(), ",");
    assert_eq!(out.replacement(), "-");
    assert_eq!(out.quoting(), QuotingMethod::Escape);
}

// ---------------------------------------------------------------------------
// START_SECTION(([EXTRA] ~SVOutStream()))
// The upstream section only deletes the pointer. The destructor's observable
// effect is closing an owned file stream; `finish` is the checked equivalent.
// ---------------------------------------------------------------------------
#[test]
fn destructor_equivalent_flushes_and_returns_the_writer() {
    let mut out = stream(",", "_", QuotingMethod::None);
    out.write_field("a").unwrap();
    let buffer = out.finish().unwrap();
    assert_eq!(buffer, b"a");
    // Dropping without finishing is legal; a Vec sink loses nothing.
    let mut out = stream(",", "_", QuotingMethod::None);
    out.write_field("b").unwrap();
    drop(out);
}

// ---------------------------------------------------------------------------
// START_SECTION((template <typename T> SVOutStream& operator<<(const T& value)))
// Two blocks, two TEST_EQUAL macros. Each accepts three platform spellings of
// -1.23e45; the pinned NumericFormatting produces the third, "-1.23e45".
// ---------------------------------------------------------------------------
#[test]
fn arithmetic_fields_use_the_source_numeric_conversion() {
    let mut out = stream(",", "_", QuotingMethod::Double);
    out.write_number(123).unwrap();
    out.write_number(3.14).unwrap();
    out.write_number(-1.23e45).unwrap();
    out.newline().unwrap();
    out.write_number(456).unwrap();
    out.end_line().unwrap();
    assert_eq!(text(out), "123,3.14,-1.23e45\n456\n");

    let mut out = stream("_/_", "_", QuotingMethod::Double);
    out.write_number(123).unwrap();
    out.write_number(3.14).unwrap();
    out.write_number(-1.23e45).unwrap();
    out.end_line().unwrap();
    out.write_number(456).unwrap();
    out.newline().unwrap();
    assert_eq!(text(out), "123_/_3.14_/_-1.23e45\n456\n");
}

// ---------------------------------------------------------------------------
// START_SECTION((SVOutStream& operator<<(std::string str)))
// Four blocks, four TEST_EQUAL macros: QuotingMethod NONE, ESCAPE, DOUBLE and
// a multi-character separator with a multi-character replacement.
// ---------------------------------------------------------------------------
#[test]
fn string_fields_are_substituted_or_quoted_as_configured() {
    let mut out = stream(",", "_", QuotingMethod::None);
    out.write_field("a").unwrap();
    out.write_field("bc").unwrap();
    out.write_field("d,f").unwrap();
    out.newline().unwrap();
    out.write_field("g\"i\"k").unwrap();
    out.write_char('l').unwrap();
    out.end_line().unwrap();
    assert_eq!(text(out), "a,bc,d_f\ng\"i\"k,l\n");

    let mut out = stream(",", "_", QuotingMethod::Escape);
    out.write_field("a").unwrap();
    out.write_field("bc").unwrap();
    out.write_field("d,f").unwrap();
    out.newline().unwrap();
    out.write_field("g\"i\"k").unwrap();
    out.write_char('l').unwrap();
    out.end_line().unwrap();
    assert_eq!(text(out), "\"a\",\"bc\",\"d,f\"\n\"g\\\"i\\\"k\",\"l\"\n");

    let mut out = stream(",", "_", QuotingMethod::Double);
    out.write_field("a").unwrap();
    out.write_field("bc").unwrap();
    out.write_field("d,f").unwrap();
    out.newline().unwrap();
    out.write_field("g\"i\"k").unwrap();
    out.write_char('l').unwrap();
    out.end_line().unwrap();
    assert_eq!(text(out), "\"a\",\"bc\",\"d,f\"\n\"g\"\"i\"\"k\",\"l\"\n");

    let mut out = stream("; ", ",_", QuotingMethod::None);
    out.write_field("a").unwrap();
    out.write_field("bc").unwrap();
    out.write_field("d; f").unwrap();
    out.newline().unwrap();
    out.write_field("g\"i\"k").unwrap();
    out.write_char('l').unwrap();
    out.end_line().unwrap();
    assert_eq!(text(out), "a; bc; d,_f\ng\"i\"k; l\n");
}

// ---------------------------------------------------------------------------
// START_SECTION((SVOutStream& operator<<(const std::string& str)))
// Upstream NOT_TESTABLE ("tested with operator<<std::string"). Rust has one
// field writer taking &str, so the by-value and by-reference source overloads
// are the same call; this pins that they agree.
// ---------------------------------------------------------------------------
#[test]
fn owned_and_borrowed_text_take_the_same_path() {
    let owned = String::from("d,f");
    let mut out = stream(",", "_", QuotingMethod::None);
    out.write_field(&owned).unwrap();
    out.write_field(owned.as_str()).unwrap();
    assert_eq!(text(out), "d_f,d_f");
}

// ---------------------------------------------------------------------------
// START_SECTION((SVOutStream& operator<<(const char* c_str)))
// Upstream NOT_TESTABLE. The source implements it as
// `operator<<(std::string(c_str))`, so a literal is quoted like any field.
// ---------------------------------------------------------------------------
#[test]
fn string_literals_are_quoted_like_any_field() {
    let mut out = stream(",", "_", QuotingMethod::Escape);
    out.write_field("d,f").unwrap();
    assert_eq!(text(out), "\"d,f\"");
}

// ---------------------------------------------------------------------------
// START_SECTION((SVOutStream& operator<<(const char c)))
// Upstream NOT_TESTABLE. The source routes through StringUtils::toStr(c).
// ---------------------------------------------------------------------------
#[test]
fn characters_are_written_as_one_character_fields() {
    let mut out = stream(",", "_", QuotingMethod::Double);
    out.write_char('l').unwrap();
    out.write_char('"').unwrap();
    out.write_char(',').unwrap();
    assert_eq!(text(out), "\"l\",\"\"\"\",\",\"");
    // A char is a Unicode scalar value here, not the source's single byte.
    let mut out = stream(",", "_", QuotingMethod::None);
    out.write_char('e').unwrap();
    out.write_char('\u{00e4}').unwrap();
    assert_eq!(text(out), "e,\u{00e4}");
}

// ---------------------------------------------------------------------------
// START_SECTION((SVOutStream& operator<<(std::ostream& (*fp)(std::ostream&))))
// One TEST_EQUAL: `out << endl << 123 << endl << "bla"`.
// ---------------------------------------------------------------------------
#[test]
fn end_line_terminates_and_resets_the_line_state() {
    let mut out = stream(",", "_", QuotingMethod::Escape);
    out.end_line().unwrap();
    out.write_number(123).unwrap();
    out.end_line().unwrap();
    out.write_field("bla").unwrap();
    assert_eq!(text(out), "\n123\n\"bla\"");
}

// ---------------------------------------------------------------------------
// START_SECTION((SVOutStream& operator<<(enum Newline)))
// One TEST_EQUAL: the same sequence with `nl`, which differs only by flushing.
// ---------------------------------------------------------------------------
#[test]
fn newline_terminates_without_flushing() {
    let mut out = stream(",", "_", QuotingMethod::Escape);
    out.newline().unwrap();
    out.write_number(123).unwrap();
    out.newline().unwrap();
    out.write_field("bla").unwrap();
    assert_eq!(out.rows(), 2);
    assert_eq!(text(out), "\n123\n\"bla\"");
}

// ---------------------------------------------------------------------------
// START_SECTION((SVOutStream& write(const std::string& str)))
// One TEST_EQUAL. The comment field is written unmodified and, crucially, the
// source's `write` does not touch `newline_`.
// ---------------------------------------------------------------------------
#[test]
fn raw_chunks_bypass_quoting_and_leave_the_line_state_alone() {
    let mut out = stream(",", "_", QuotingMethod::Escape);
    out.write_field("bla").unwrap();
    out.write_number(123).unwrap();
    out.newline().unwrap();
    out.write_raw("#This, is, a, comment\n").unwrap();
    out.write_number(4.56).unwrap();
    out.write_field("test").unwrap();
    out.end_line().unwrap();
    assert_eq!(
        text(out),
        "\"bla\",123\n#This, is, a, comment\n4.56,\"test\"\n"
    );

    // The quirk the source header warns about: a raw chunk written mid-line
    // leaves the stream believing it is mid-line, so the next field still
    // receives a separator even though the cursor is at column zero.
    let mut out = stream(",", "_", QuotingMethod::None);
    out.write_field("a").unwrap();
    out.write_raw("\n").unwrap();
    out.write_field("b").unwrap();
    assert!(!out.at_line_start());
    assert_eq!(text(out), "a\n,b");
}

// ---------------------------------------------------------------------------
// START_SECTION((bool modifyStrings(bool modify)))
// Three TEST_EQUAL macros: the two returned previous states and the output.
// ---------------------------------------------------------------------------
#[test]
fn modify_strings_returns_the_previous_state() {
    let mut out = stream(",", "_", QuotingMethod::Double);
    out.write_field("test").unwrap();
    assert!(out.set_modify_strings(false));
    out.write_field("bla").unwrap();
    assert!(!out.set_modify_strings(true));
    out.write_field("laber").unwrap();
    out.end_line().unwrap();
    assert_eq!(text(out), "\"test\",bla,\"laber\"\n");
}

// ---------------------------------------------------------------------------
// START_SECTION((template <typename NumericT> SVOutStream& writeValueOrNan(NumericT thing)))
// One TEST_EQUAL covering an integer, a finite double, NaN and both infinities.
// ---------------------------------------------------------------------------
#[test]
fn write_value_or_nan_substitutes_the_nonfinite_texts() {
    let mut out = stream(",", "_", QuotingMethod::Double);
    out.write_value_or_nan(123).unwrap();
    out.write_value_or_nan(3.14).unwrap();
    out.newline().unwrap();
    out.write_value_or_nan(456).unwrap();
    out.write_value_or_nan(f64::NAN).unwrap();
    out.newline().unwrap();
    out.write_value_or_nan(f64::INFINITY).unwrap();
    out.write_value_or_nan(f64::NEG_INFINITY).unwrap();
    out.end_line().unwrap();
    assert_eq!(text(out), "123,3.14\n456,nan\ninf,-inf\n");
}

// ---------------------------------------------------------------------------
// Native tests beyond the upstream sections.
// ---------------------------------------------------------------------------

#[test]
fn nonfinite_substitution_restores_the_modification_state() {
    let mut out = stream(",", "_", QuotingMethod::Double);
    out.write_value_or_nan(f64::NAN).unwrap();
    assert!(out.modify_strings());
    out.write_field("x").unwrap();
    assert_eq!(text(out), "nan,\"x\"");
    // With modification already off, the previous state is restored as off.
    let mut out = stream(",", "_", QuotingMethod::Double);
    out.set_modify_strings(false);
    out.write_value_or_nan(f32::NAN).unwrap();
    assert!(!out.modify_strings());
    assert_eq!(text(out), "nan");
}

#[test]
fn a_newline_in_a_field_is_rejected_before_any_byte_is_written() {
    let mut out = stream(",", "_", QuotingMethod::Double);
    out.write_field("a").unwrap();
    let error = out.write_field("b\nc").unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    // No separator was emitted for the rejected field and the state is intact.
    assert_eq!(out.row_fields(), 1);
    out.write_field("d").unwrap();
    assert_eq!(text(out), "\"a\",\"d\"");
}

#[test]
fn a_carriage_return_in_a_field_is_rejected_too() {
    let mut out = stream(",", "_", QuotingMethod::None);
    let error = out.write_field("a\rb").unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    let error = out.write_display("a\rb").unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
}

#[test]
fn a_newline_through_the_display_writer_is_rejected() {
    let mut out = stream(",", "_", QuotingMethod::None);
    let error = out.write_display("a\nb").unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
}

#[test]
fn display_fields_are_never_quoted_or_substituted() {
    let mut out = stream(",", "_", QuotingMethod::Double);
    out.write_display("d,f").unwrap();
    out.write_display(42).unwrap();
    assert_eq!(text(out), "d,f,42");
}

#[test]
fn non_ascii_fields_survive_quoting_and_substitution() {
    // A multi-byte separator and a multi-byte field: nothing may be sliced at
    // a byte offset inside a character.
    let mut out = stream("\u{2014}", "\u{2026}", QuotingMethod::None);
    out.write_field("\u{65e5}\u{672c}\u{8a9e}").unwrap();
    out.write_field("a\u{2014}b").unwrap();
    out.newline().unwrap();
    out.write_field("\u{1f600}").unwrap();
    assert_eq!(
        text(out),
        "\u{65e5}\u{672c}\u{8a9e}\u{2014}a\u{2026}b\n\u{1f600}"
    );

    let mut out = stream("\t", "_", QuotingMethod::Escape);
    out.write_field("\u{e4}\"\u{f6}\\\u{fc}").unwrap();
    assert_eq!(text(out), "\"\u{e4}\\\"\u{f6}\\\\\u{fc}\"");
}

#[test]
fn an_empty_or_broken_separator_is_rejected() {
    let error = SVOutStream::with_options(Vec::new(), "", "_", QuotingMethod::None).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    let error = SVOutStream::with_options(Vec::new(), "\n", "_", QuotingMethod::None).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    let error =
        SVOutStream::with_options(Vec::new(), ",", "a\nb", QuotingMethod::None).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
}

#[test]
fn field_row_and_line_ceilings_are_checked_before_rendering() {
    let limits = SVLimits {
        max_field_bytes: 8,
        max_row_fields: 2,
        max_rows: 1,
    };
    let mut out = stream(",", "_", QuotingMethod::None).with_limits(limits);
    // Worst-case quoting is charged, so a 4-byte field already exceeds 8.
    let error = out.write_field("abcd").unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    out.write_field("ab").unwrap();
    out.write_field("cd").unwrap();
    let error = out.write_field("ef").unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    out.newline().unwrap();
    // The field count resets with the line, the line count does not.
    out.write_field("gh").unwrap();
    let error = out.newline().unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    let error = out.write_raw("123456789").unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    assert_eq!(text(out), "ab,cd\ngh");
}

#[test]
fn a_replacement_longer_than_the_separator_is_charged_before_the_substitution() {
    // The worst case `preflight` charges - twice the field plus two quotes -
    // bounds both quote paths but not the QuotingMethod::None substitution,
    // whose growth factor is replacement.len() / separator.len(). With a
    // 4096-byte replacement, a 400-byte field of commas renders to 1,638,400
    // bytes, so the ceiling has to be charged from the separator count rather
    // than from the raw length.
    let replacement = "R".repeat(4096);
    let limits = SVLimits {
        max_field_bytes: 1024,
        ..SVLimits::default()
    };
    let mut out = SVOutStream::with_options(Vec::new(), ",", &replacement, QuotingMethod::None)
        .expect("valid separator and replacement")
        .with_limits(limits);
    let error = out.write_field(&",".repeat(400)).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    // Nothing was emitted and the line state did not advance.
    assert_eq!(text(out), "");

    // A field that fits after substitution is still written: one separator
    // becomes 4096 bytes, which is under a 8192-byte ceiling.
    let limits = SVLimits {
        max_field_bytes: 8192,
        ..SVLimits::default()
    };
    let mut out = SVOutStream::with_options(Vec::new(), ",", &replacement, QuotingMethod::None)
        .expect("valid separator and replacement")
        .with_limits(limits);
    out.write_field("a,b").expect("3 bytes become 4098");
    assert_eq!(text(out), format!("a{replacement}b"));
}

#[test]
fn the_file_constructor_creates_and_truncates() {
    let dir = std::env::temp_dir().join("openms_sv_out_stream_file");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("table.csv");
    std::fs::write(&path, b"stale contents that must not survive").unwrap();
    let mut out = SVOutStream::create_with_options(&path, ",", "_", QuotingMethod::None).unwrap();
    out.write_field("a").unwrap();
    out.write_number(1).unwrap();
    out.newline().unwrap();
    out.finish().unwrap();
    assert_eq!(std::fs::read_to_string(&path).unwrap(), "a,1\n");

    let mut out = SVOutStream::create(&path).unwrap();
    assert_eq!(out.separator(), "\t");
    out.write_field("x").unwrap();
    out.finish().unwrap();
    assert_eq!(std::fs::read_to_string(&path).unwrap(), "\"x\"");

    let error = SVOutStream::create(dir.join("missing").join("table.csv")).unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error}");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn quote_follows_the_source_substitution_order() {
    // StringUtils::quote escapes backslashes first, then quote characters.
    assert_eq!(quote("a\\\"b", QuotingMethod::Escape), "\"a\\\\\\\"b\"");
    assert_eq!(quote("a\"b", QuotingMethod::Double), "\"a\"\"b\"");
    // QuotingMethod::NONE still adds the surrounding quotes in the source,
    // which is why SVOutStream substitutes separators instead of calling it.
    assert_eq!(quote("a,b", QuotingMethod::None), "\"a,b\"");
    assert_eq!(quote("", QuotingMethod::Double), "\"\"");
}

#[test]
fn numeric_text_matches_the_source_formatter() {
    // Fixed notation with 15 fractional digits, trailing zeros trimmed and
    // one digit always kept after the point.
    assert_eq!(source_float_text(3.14, F64_FIXED_DIGITS), "3.14");
    assert_eq!(source_float_text(4.56, F64_FIXED_DIGITS), "4.56");
    assert_eq!(source_float_text(5.0, F64_FIXED_DIGITS), "5.0");
    assert_eq!(source_float_text(0.0, F64_FIXED_DIGITS), "0.0");
    assert_eq!(source_float_text(-0.0, F64_FIXED_DIGITS), "-0.0");
    assert_eq!(source_float_text(0.01, F64_FIXED_DIGITS), "0.01");
    assert_eq!(source_float_text(9999.5, F64_FIXED_DIGITS), "9999.5");
    // Scientific notation from 1e4 up and below 1e-2, '+' dropped from the
    // exponent, two-digit zero padding kept, mantissa given a ".0".
    assert_eq!(source_float_text(1e4, F64_FIXED_DIGITS), "1.0e04");
    assert_eq!(source_float_text(-1.23e45, F64_FIXED_DIGITS), "-1.23e45");
    assert_eq!(source_float_text(1e-3, F64_FIXED_DIGITS), "1.0e-03");
    assert_eq!(source_float_text(1.5e-10, F64_FIXED_DIGITS), "1.5e-10");
    assert_eq!(source_float_text(1e308, F64_FIXED_DIGITS), "1.0e308");
    // Nonfinite spellings: uppercase NaN, lowercase infinities.
    assert_eq!(source_float_text(f64::NAN, F64_FIXED_DIGITS), "NaN");
    assert_eq!(source_float_text(f64::INFINITY, F64_FIXED_DIGITS), "inf");
    assert_eq!(
        source_float_text(f64::NEG_INFINITY, F64_FIXED_DIGITS),
        "-inf"
    );
    // f32 uses six fractional digits in fixed notation.
    assert_eq!(source_float_text(1.0 / 3.0, F32_FIXED_DIGITS), "0.333333");
    assert_eq!(0.5f32.sv_text(), "0.5");
    assert_eq!(3.14f64.sv_text(), "3.14");
    // A precision above MAX_FIXED_DIGITS is clamped rather than honoured. The
    // source has no clamp: it overflows its char buf[64] and falls through to
    // std::to_string(double), which would print "0.333333"
    // (NumericFormatting.h:136-139). Neither is reachable from the source's own
    // call sites, which pass at most writtenDigits<long double>() = 18.
    assert_eq!(
        source_float_text(1.0 / 3.0, MAX_FIXED_DIGITS + 1),
        source_float_text(1.0 / 3.0, MAX_FIXED_DIGITS)
    );
}

/// `appendNumeric` is a template: for a `float` the branch comparison and
/// `std::to_chars` both run at `float` width, so an `f32` must not be promoted
/// to `f64` before it is formatted.
#[test]
fn f32_text_is_formatted_at_f32_width_not_through_f64() {
    // The third column says whether the promoted `f64` pipeline spells the
    // value differently - the defect this separation fixes.
    for (value, expected, promotion_differs) in [
        (1.23e-5f32, "1.23e-05", true),
        (12345.6f32, "1.23456e04", true),
        (1.23e10f32, "1.23e10", true),
        (0.001f32, "1.0e-03", true),
        // 0.01f32 == float(1e-2) exactly, so `abs_val < T(1e-2)` is false and
        // the fixed branch prints it; the promoted f64 is below 1e-2 and would
        // take the scientific branch.
        (0.01f32, "0.01", true),
        // Both pipelines agree on a short value inside the fixed range.
        (0.5f32, "0.5", false),
        (1.0f32 / 3.0f32, "0.333333", false),
    ] {
        assert_eq!(
            source_f32_text(value, F32_FIXED_DIGITS),
            expected,
            "{value}"
        );
        assert_eq!(value.sv_text(), expected, "{value} through SvNumber");
        let promoted = source_float_text(f64::from(value), F32_FIXED_DIGITS);
        assert_eq!(
            promoted != expected,
            promotion_differs,
            "{value} promoted is {promoted:?}, f32 text is {expected:?}"
        );
    }
    assert_eq!(source_f32_text(f32::NAN, F32_FIXED_DIGITS), "NaN");
    assert_eq!(source_f32_text(f32::INFINITY, F32_FIXED_DIGITS), "inf");
    assert_eq!(source_f32_text(f32::NEG_INFINITY, F32_FIXED_DIGITS), "-inf");
    assert_eq!(source_f32_text(0.0f32, F32_FIXED_DIGITS), "0.0");
    assert_eq!(source_f32_text(-0.0f32, F32_FIXED_DIGITS), "-0.0");
}

#[test]
fn every_integer_width_renders_and_classifies_as_finite() {
    assert_eq!(i8::MIN.sv_text(), "-128");
    assert_eq!(u8::MAX.sv_text(), "255");
    assert_eq!(i16::MIN.sv_text(), "-32768");
    assert_eq!(u16::MAX.sv_text(), "65535");
    assert_eq!(i32::MIN.sv_text(), "-2147483648");
    assert_eq!(u32::MAX.sv_text(), "4294967295");
    assert_eq!(i64::MIN.sv_text(), "-9223372036854775808");
    assert_eq!(u64::MAX.sv_text(), "18446744073709551615");
    assert_eq!(
        i128::MIN.sv_text(),
        "-170141183460469231731687303715884105728"
    );
    assert_eq!(
        u128::MAX.sv_text(),
        "340282366920938463463374607431768211455"
    );
    assert_eq!(7usize.sv_text(), "7");
    assert_eq!((-7isize).sv_text(), "-7");
    assert_eq!(0u8.sv_class(), NumberClass::Finite);
    assert_eq!(0.0f64.sv_class(), NumberClass::Finite);
    assert_eq!(f32::NAN.sv_class(), NumberClass::NotANumber);
    assert_eq!(
        f64::NEG_INFINITY.sv_class(),
        NumberClass::Infinite { negative: true }
    );
    assert_eq!(
        f32::INFINITY.sv_class(),
        NumberClass::Infinite { negative: false }
    );
}

#[test]
fn a_failing_writer_reports_io_and_keeps_the_stream_usable() {
    struct Refusing;
    impl std::io::Write for Refusing {
        fn write(&mut self, _buf: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("refused"))
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut out = SVOutStream::with_options(Refusing, ",", "_", QuotingMethod::None).unwrap();
    let error = out.write_field("a").unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error}");
    // The rejected field did not advance the line state.
    assert!(out.at_line_start());
    assert_eq!(out.row_fields(), 0);
    let error = out.newline().unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error}");
    let error = out.write_raw("x").unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error}");
}

#[test]
fn flush_is_available_without_consuming_the_stream() {
    let mut out = stream(",", "_", QuotingMethod::None);
    out.write_field("a").unwrap();
    out.flush().unwrap();
    out.write_field("b").unwrap();
    assert_eq!(text(out), "a,b");
}
