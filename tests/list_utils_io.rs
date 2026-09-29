// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `DATASTRUCTURES/ListUtilsIO.h`, `CONCEPT/PrecisionWrapper.h`,
//! `CONCEPT/Detail/NumericFormatting.h` and the three-fraction-digit
//! `StringUtils::toStr(value, false)` they build on.
//!
//! Evidence, per test:
//!
//! - tier 3, upstream class-test literals at core `bc9cc12`:
//!   `PrecisionWrapper_test.cpp` (all four sections), `ListUtilsIO_test.cpp`
//!   (its one section) and the `full_precision = false` rows of the
//!   `StringUtils::toStr(float, bool)` and `toStr(double, bool)` sections of
//!   `String_test.cpp`;
//! - tier 4, independently derived: the other low-precision texts were computed
//!   in Python 3 with `'%.3e'` and `'%.3f'` (correctly rounded, ties to even,
//!   on the exact binary value) followed by the header's trimming rule, and the
//!   `float` rows over the `float` value obtained with
//!   `struct.unpack('<f', struct.pack('<f', x))`.

// The source class-test literals are kept verbatim (digit for digit, including
// the `-3.14159` probe), so clippy's literal-shape lints are silenced here.
#![allow(
    clippy::excessive_precision,
    clippy::approx_constant,
    clippy::inconsistent_digit_grouping
)]

use openms::concept::numeric_formatting::{append_numeric_f32, append_numeric_f64};
use openms::concept::precision_wrapper::{PrecisionWrapper, precision_wrapper};
use openms::format::file_info::text_format::{
    MAX_LIST_TEXT_BYTES, StringListPush, VecLowPrecision, list_to_string,
    list_to_string_low_precision, to_str, to_str_f32, to_str_low_precision,
    to_str_low_precision_f32,
};

// ------------------------------------------------------------ PrecisionWrapper

#[test]
fn class_test_precision_wrapper_constructor_stores_the_value() {
    // START_SECTION(explicit PrecisionWrapper(const FloatingPointType rhs))
    let wrapped = PrecisionWrapper::<f64>::new(3.5);
    assert_eq!(wrapped.value, 3.5);
}

#[test]
fn class_test_precision_wrapper_copy_constructor() {
    // START_SECTION(PrecisionWrapper(const PrecisionWrapper& rhs))
    let wrapped = PrecisionWrapper::<f64>::new(2.25);
    let copy = wrapped;
    assert_eq!(copy.value, 2.25);
    assert_eq!(copy, wrapped);
}

#[test]
fn class_test_precision_wrapper_make_function_deduces_the_type() {
    // START_SECTION(precisionWrapper(const FloatingPointType rhs))
    let double: PrecisionWrapper<f64> = precision_wrapper(1.5);
    assert_eq!(double.value, 1.5);
    let single: PrecisionWrapper<f32> = precision_wrapper(0.5f32);
    assert_eq!(single.value, 0.5f32);
}

#[test]
#[allow(clippy::excessive_precision)]
fn class_test_precision_wrapper_output_operator() {
    // START_SECTION(operator<<(std::ostream&, const PrecisionWrapper<T>&))
    assert_eq!(
        precision_wrapper(0.123_456_789_012_345_678_9).to_string(),
        "0.123456789012346"
    );
    assert_eq!(
        precision_wrapper(0.123_456_789_012_345_678_9_f32).to_string(),
        "0.123457"
    );
    assert_eq!(precision_wrapper(1.0).to_string(), "1.0");
    assert_eq!(precision_wrapper(100.5).to_string(), "100.5");
}

#[test]
fn precision_wrapper_writes_to_str_and_ignores_formatter_precision() {
    for value in [
        0.0,
        -0.0,
        17.012345,
        1.0e6,
        1.0e-5,
        -135_169.706_298_828_125,
    ] {
        assert_eq!(precision_wrapper(value).to_string(), to_str(value));
        // The source leaves the stream's precision untouched and does not use
        // it; a Rust formatter's precision and width are likewise ignored.
        assert_eq!(format!("{:.2}", precision_wrapper(value)), to_str(value));
    }
    for value in [
        0.0f32,
        50_254.2,
        1.0e6,
        2_203_265.25,
        f32::MIN_POSITIVE / 10.0,
    ] {
        assert_eq!(precision_wrapper(value).to_string(), to_str_f32(value));
    }
    assert_eq!(precision_wrapper(f64::NAN).to_string(), "NaN");
    assert_eq!(precision_wrapper(f64::NEG_INFINITY).to_string(), "-inf");
    assert_eq!(precision_wrapper(f32::INFINITY).to_string(), "inf");
}

// ------------------------------------------------------------ toStr(x, false)

#[test]
fn class_test_to_str_float_low_precision() {
    // String_test.cpp, StringUtils::toStr(float, bool full_precision)
    assert_eq!(to_str_low_precision_f32(17.0123), "17.012");
    let denorm = f32::MIN_POSITIVE / 10.0;
    assert_eq!(to_str_low_precision_f32(denorm), "1.175e-39");
    assert_eq!(to_str_low_precision_f32(f32::NAN), "NaN");
}

#[test]
fn class_test_to_str_double_low_precision() {
    // String_test.cpp, StringUtils::toStr(double, bool full_precision)
    assert_eq!(to_str_low_precision(17.012345), "17.012");
    let denorm = f64::MIN_POSITIVE / 10.0;
    assert_eq!(to_str_low_precision(denorm), "2.225e-309");
    assert_eq!(to_str_low_precision(f64::NAN), "NaN");
}

#[test]
fn low_precision_double_texts_match_python() {
    let cases = [
        (2.0, "2.0"),
        (0.0, "0.0"),
        (-0.0, "-0.0"),
        (1.0e6, "1.0e06"),
        // Exact tie 1.2345e4 goes to the even digit.
        (12_345.0, "1.234e04"),
        // Exact binary ties in fixed notation go to the even digit.
        (0.0625, "0.062"),
        (0.1875, "0.188"),
        // 0.0125 is stored just above the tie.
        (0.0125, "0.013"),
        (0.0105, "0.011"),
        // Rounding carries into a fifth integer digit but stays fixed.
        (9_999.9996, "10000.0"),
        (0.00999, "9.99e-03"),
        // 1.0005 is stored just below the tie, 2.0005 just above.
        (1.0005, "1.0"),
        (2.0005, "2.001"),
        (123.4565, "123.457"),
        (-3.14159, "-3.142"),
        (1.0e-5, "1.0e-05"),
        (1.0e300, "1.0e300"),
        (5.0e-324, "4.941e-324"),
        (f64::INFINITY, "inf"),
        (f64::NEG_INFINITY, "-inf"),
        (-f64::NAN, "NaN"),
    ];
    for (value, expected) in cases {
        assert_eq!(to_str_low_precision(value), expected, "{value:e}");
    }
}

#[test]
fn low_precision_float_texts_use_float_thresholds_and_values() {
    let cases = [
        (17.0123f32, "17.012"),
        (50_254.199_219f32, "5.025e04"),
        (1.0e6f32, "1.0e06"),
        (0.1f32, "0.1"),
        (2.0005f32, "2.0"),
        (0.0125f32, "0.013"),
        (f32::NEG_INFINITY, "-inf"),
    ];
    for (value, expected) in cases {
        assert_eq!(to_str_low_precision_f32(value), expected, "{value:e}");
    }
    // The float 0.01f is not below the float threshold 1e-2f, so it is fixed;
    // promoted to double it lies below 1e-2 and becomes scientific.
    let hundredth = 0.01f32;
    assert_eq!(to_str_low_precision_f32(hundredth), "0.01");
    assert_eq!(to_str_low_precision(f64::from(hundredth)), "1.0e-02");
}

// ------------------------------------------------------------ ListUtilsIO

#[test]
fn class_test_string_list_push_chains() {
    // ListUtilsIO_test.cpp, [EXTRA] StringList& operator<<(StringList&, const StringType&)
    let mut list: Vec<String> = Vec::new();
    list.push_item("a")
        .push_item("b")
        .push_item("c")
        .push_item("a");
    assert_eq!(list.len(), 4);
    assert_eq!(list[0], "a");
    assert_eq!(list[1], "b");
    assert_eq!(list[2], "c");
    assert_eq!(list[3], "a");
}

#[test]
fn vec_low_precision_writes_the_bracketed_list() {
    // The DataValue(DoubleList) row of String_test.cpp has the same shape:
    // toStr(DataValue(DoubleList({17.012345, 2.0})), false) == "[17.012, 2.0]".
    let values = vec![17.012_345f64, 2.0];
    assert_eq!(VecLowPrecision::new(&values).to_string(), "[17.012, 2.0]");
    assert_eq!(VecLowPrecision::from(&values).to_string(), "[17.012, 2.0]");
    assert_eq!(
        VecLowPrecision::new(&values).to_text().ok().as_deref(),
        Some("[17.012, 2.0]")
    );
    assert_eq!(
        list_to_string_low_precision(&values).ok().as_deref(),
        Some("[17.012, 2.0]")
    );
    // Full precision through the unwrapped operator, for contrast.
    assert_eq!(
        list_to_string(&values).ok().as_deref(),
        Some("[17.012345, 2.0]")
    );

    let empty: [f32; 0] = [];
    assert_eq!(VecLowPrecision::from(&empty[..]).to_string(), "[]");
    assert_eq!(
        list_to_string_low_precision(&empty).ok().as_deref(),
        Some("[]")
    );

    let floats = [17.0123f32, 1.0e6, f32::NAN];
    assert_eq!(
        VecLowPrecision::new(&floats).to_string(),
        "[17.012, 1.0e06, NaN]"
    );
    let references = [&1.0005f64, &2.0005];
    assert_eq!(
        list_to_string_low_precision(&references).ok().as_deref(),
        Some("[1.0, 2.001]")
    );
}

#[test]
fn low_precision_list_ceiling_is_checked_before_allocating() {
    // Each element may need up to 12 bytes plus a 2-byte separator, so a
    // vector of this length cannot be guaranteed to fit.
    let count = MAX_LIST_TEXT_BYTES / 14 + 1;
    let values = vec![0.0f64; count];
    assert!(list_to_string_low_precision(&values).is_err());
    // The streaming Display form has no ceiling; spot-check its start.
    let text = VecLowPrecision::new(&values[..3]).to_string();
    assert_eq!(text, "[0.0, 0.0, 0.0]");
}

// ------------------------------------------------------------ appendNumeric

fn numeric(value: f64, precision: usize, fixed_format: bool) -> String {
    let mut text = String::from("prefix:");
    append_numeric_f64(value, &mut text, precision, fixed_format);
    text
}

#[test]
fn append_numeric_appends_and_matches_the_to_str_paths() {
    for value in [
        0.0,
        -0.0,
        17.012345,
        1.0e6,
        1.0e-5,
        12_345.0,
        0.0625,
        f64::NAN,
    ] {
        assert_eq!(
            numeric(value, 15, false),
            format!("prefix:{}", to_str(value))
        );
        assert_eq!(
            numeric(value, 3, true),
            format!("prefix:{}", to_str_low_precision(value))
        );
    }
    for value in [0.0f32, 17.0123, 50_254.2, 1.0e6, 0.01] {
        let mut full = String::new();
        append_numeric_f32(value, &mut full, 6, false);
        assert_eq!(full, to_str_f32(value));
        let mut low = String::new();
        append_numeric_f32(value, &mut low, 3, true);
        assert_eq!(low, to_str_low_precision_f32(value));
    }
}

#[test]
fn append_numeric_fixed_precision_rules() {
    // Precision 0: to_chars writes no point, so ".0" is appended.
    assert_eq!(numeric(2.5, 0, false), "prefix:2.0");
    assert_eq!(numeric(3.5, 0, false), "prefix:4.0");
    assert_eq!(numeric(-0.4, 0, false), "prefix:-0.0");
    // Other precisions round the exact value, ties to even, then trim.
    assert_eq!(numeric(1.25, 1, false), "prefix:1.2");
    assert_eq!(numeric(1.2, 5, false), "prefix:1.2");
    // Scientific magnitudes ignore the precision: shortest, or three digits.
    assert_eq!(numeric(123_456.0, 0, false), "prefix:1.23456e05");
    assert_eq!(numeric(123_456.0, 9, true), "prefix:1.235e05");
}

#[test]
fn append_numeric_falls_back_to_percent_f_when_the_buffer_overflows() {
    // 62 fraction digits fill the 64-byte buffer exactly and still fit.
    let fits = numeric(1.5, 62, false);
    assert_eq!(fits, "prefix:1.5");
    // One more does not: to_chars fails and std::to_string's "%f" is used,
    // six fraction digits, untrimmed.
    assert_eq!(numeric(1.5, 63, false), "prefix:1.500000");
    assert_eq!(numeric(-1234.5, 60, false), "prefix:-1234.500000");
    assert_eq!(numeric(0.0, usize::MAX, false), "prefix:0.000000");
    let mut text = String::new();
    append_numeric_f32(0.1, &mut text, 100, false);
    // static_cast<double>(0.1f) is 0.100000001490116...
    assert_eq!(text, "0.100000");
}

#[test]
fn vector_operator_covers_char_bool_integer_and_string_elements() {
    // operator<<(std::ostream&, const std::vector<T>&) with toStr(char), a bool
    // promoted to int, integers and strings.
    assert_eq!(list_to_string(&['a', 'b']).ok().as_deref(), Some("[a, b]"));
    assert_eq!(
        list_to_string(&[true, false]).ok().as_deref(),
        Some("[1, 0]")
    );
    assert_eq!(list_to_string(&[-3i32, 4]).ok().as_deref(), Some("[-3, 4]"));
    assert_eq!(
        list_to_string(&["x", "y z"]).ok().as_deref(),
        Some("[x, y z]")
    );
}
