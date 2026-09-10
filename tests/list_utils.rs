// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// Source locations and boundary policies: tests/data/list_utils_provenance.json.
use openms::data_structures::list::{self, CaseSensitivity, ListFormat, ListParse};
use openms::metadata::{MetaValue, MetaValueData};
use openms::param::ParamValue;
use std::borrow::Cow;

#[test]
fn source_literal_create_overloads_and_untrimmed_strings() {
    assert_eq!(
        list::create::<String>("yes,no, maybe", b',').unwrap(),
        ["yes", "no", " maybe"]
    );
    assert_eq!(list::create::<f64>("1.2,3.5", b',').unwrap(), [1.2, 3.5]);
    assert_eq!(list::create::<i32>("1,5", b',').unwrap(), [1, 5]);
    assert_eq!(list::create::<i32>("2", b',').unwrap(), [2]);
    assert!(list::create::<i32>("", b',').unwrap().is_empty());
    assert!(list::create::<String>("", b',').unwrap().is_empty());
    assert_eq!(
        list::create::<String>("yes,no", b',').unwrap(),
        ["yes", "no"]
    );
    assert_eq!(list::create::<String>("no", b',').unwrap(), ["no"]);
    let source = ["test string", "string2", "last string"];
    assert_eq!(
        list::create::<String>("test string,string2,last string", b',').unwrap(),
        source
    );
    assert_eq!(
        list::create::<String>("test string#string2#last string", b'#').unwrap(),
        source
    );
    let values = ["1.2", "1.56", "10.4"];
    assert_eq!(
        list::create_from_strings::<String>(&values).unwrap(),
        values
    );
    assert_eq!(
        list::create_from_strings::<String>(&["1.2", "1.56", "10.4", "a"]).unwrap(),
        ["1.2", "1.56", "10.4", "a"]
    );
    assert_eq!(
        list::create_from_strings::<f64>(&values).unwrap(),
        [1.2, 1.56, 10.4]
    );
    assert!(list::create_from_strings::<f64>(&["1.2", "1.56", "10.4", "a"]).is_err());
    assert_eq!(
        list::create_from_strings::<String>(&[" x ", "\ty\r\n"]).unwrap(),
        [" x ", "\ty\r\n"]
    );
}

#[test]
fn source_contains_index_and_concatenation_literals() {
    for x in 1..=4 {
        assert!(list::contains(&[1, 2, 3, 4], &x));
    }
    for x in [5, 1011] {
        assert!(!list::contains(&[1, 2, 3, 4], &x));
    }
    let strings = [String::from("yes"), String::from("no")];
    for x in ["yes", "no"] {
        assert!(list::contains(&strings, x));
    }
    for x in ["jup", "", "noe"] {
        assert!(!list::contains(&strings, x));
    }
    let integers = [4, 3, 1, 2];
    for (value, expected) in [
        (0, None),
        (1, Some(2)),
        (2, Some(3)),
        (3, Some(1)),
        (4, Some(0)),
        (5, None),
    ] {
        assert_eq!(list::get_index(&integers, &value), expected);
    }
    let words = ["four", "three", "one", "two"].map(String::from);
    for (value, expected) in [
        ("zero", None),
        ("one", Some(2)),
        ("two", Some(3)),
        ("three", Some(1)),
        ("four", Some(0)),
        ("five", None),
    ] {
        assert_eq!(list::get_index(&words, value), expected);
    }
    assert_eq!(list::get_index(&[2, 1, 2], &2), Some(0));
    assert_eq!(
        list::concatenate(["1", "2", "3", "4", "5"], "g").unwrap(),
        "1g2g3g4g5"
    );
    assert_eq!(
        list::concatenate(["1", "2", "3", "4", "5"], "").unwrap(),
        "12345"
    );
    assert_eq!(list::concatenate(Vec::<String>::new(), "g").unwrap(), "");
    assert_eq!(
        list::concatenate(["1\n", "2\n", "3\n"], "").unwrap(),
        "1\n2\n3\n"
    );
    assert_eq!(list::concatenate((1..=3).rev(), "|").unwrap(), "3|2|1");
}

#[test]
fn strict_float_tolerance_and_ascii_case_are_not_normalized() {
    for value in [1.2, 3.4] {
        assert!(list::contains_f64(&[1.2, 3.4], value));
    }
    for value in [1.21, 1.19, 4.2, 2.0, 0.0] {
        assert!(!list::contains_f64(&[1.2, 3.4], value));
    }
    for value in [1.21, 1.19] {
        assert!(list::contains_approx(&[1.2, 3.4], value, 0.02));
    }
    assert!(!list::contains_approx(&[0.0], 0.5, 0.5));
    assert!(list::contains_approx(
        &[0.0],
        0.5,
        f64::from_bits(0.5f64.to_bits() + 1)
    ));
    for tolerance in [0.0, -0.1, f64::NAN] {
        assert!(!list::contains_approx(&[1.0], 1.0, tolerance));
    }
    assert!(list::contains_approx(&[1.0], 2.0, f64::INFINITY));
    assert!(!list::contains_approx(
        &[f64::INFINITY],
        f64::INFINITY,
        f64::INFINITY
    ));
    assert!(!list::contains_approx(&[f64::NAN], 1.0, 1.0));
    assert!(list::contains_string(
        &["YES", "no"],
        "yes",
        CaseSensitivity::Insensitive
    ));
    assert!(!list::contains_string(
        &["YES"],
        "yes",
        CaseSensitivity::Sensitive
    ));
    assert!(!list::contains_string(
        &[" YES "],
        "yes",
        CaseSensitivity::Insensitive
    ));
    assert!(!list::contains_string(
        &["Ä"],
        "ä",
        CaseSensitivity::Insensitive
    ));
    assert!(list::contains_string(
        &["ÄBC"],
        "Äbc",
        CaseSensitivity::Insensitive
    ));
}

#[test]
fn literal_delimiters_empty_fields_quotes_and_utf8() {
    assert_eq!(
        list::create::<String>(",a,,", b',').unwrap(),
        ["", "a", "", ""]
    );
    assert_eq!(
        list::create::<String>("\"a,b\",c", b',').unwrap(),
        ["\"a", "b\"", "c"]
    );
    assert_eq!(list::create::<String>("a\0b", 0).unwrap(), ["a", "b"]);
    assert_eq!(list::create::<String>("α,β", b',').unwrap(), ["α", "β"]);
    assert!(list::create::<String>("α", 0xb1).is_err());
    for bad in ["1,,2", "1,", ",1"] {
        assert!(list::create::<i32>(bad, b',').is_err());
    }
}

#[test]
fn all_numeric_types_trim_only_source_whitespace_and_check_whole_tokens() {
    assert_eq!(
        list::create::<i32>(" \t+1\r\n, -2 ,+-3,2147483647,-2147483648", b',').unwrap(),
        [1, -2, -3, i32::MAX, i32::MIN]
    );
    assert_eq!(
        list::create_from_strings::<f32>(&[" \t1.5\n", "+-2.25\r"]).unwrap(),
        [1.5, -2.25]
    );
    assert_eq!(
        list::create::<f64>(".5,1.,1e-2,+-2", b',').unwrap(),
        [0.5, 1.0, 0.01, -2.0]
    );
    for text in [
        "",
        "++1",
        "1 2",
        "1.5",
        "2147483648",
        "-2147483649",
        "1e2",
        "\u{a0}1",
        "\u{b}1",
    ] {
        assert!(i32::from_list_item(text).is_err(), "{text:?}");
    }
    for text in [
        "", "++1", "1 2", "1e", "1x", "0x1p1", "1e9999", "\u{a0}1", "\u{b}1", "ααα",
    ] {
        assert!(f32::from_list_item(text).is_err(), "{text:?}");
        assert!(f64::from_list_item(text).is_err(), "{text:?}");
    }
    assert!(f32::from_list_item("3.5e38").is_err());
    assert!(f64::from_list_item("3.5e38").is_ok());
    for text in ["1e-9999", "-1e-9999", ".0001e-9999"] {
        assert!(f64::from_list_item(text).is_err());
        assert!(f32::from_list_item(text).is_err());
    }
    assert_eq!(
        f64::from_list_item("-0e-9999").unwrap().to_bits(),
        (-0.0f64).to_bits()
    );
    assert_eq!(
        f32::from_list_item("-0.0e-9999").unwrap().to_bits(),
        (-0.0f32).to_bits()
    );
    assert_eq!(f64::from_list_item("5e-324").unwrap().to_bits(), 1);
    assert_eq!(f32::from_list_item("1e-45").unwrap().to_bits(), 1);
    // Exact decimal just above a binary32 midpoint: direct float parsing must
    // not first round it to the midpoint in binary64.
    assert_eq!(
        f32::from_list_item("1.0000000596046447753906250000000000000001")
            .unwrap()
            .to_bits(),
        1.0f32.to_bits() + 1
    );
}

#[test]
fn special_float_spellings_preserve_source_nan_branches() {
    for text in [
        "nan",
        "NaN",
        "NAN()",
        "nan(payload)",
        "nan(a b)",
        "nan(α)",
        " +NaN ",
        "-nan(foo_1)",
        "+-nan",
    ] {
        assert!(f64::from_list_item(text).unwrap().is_nan(), "{text}");
        assert!(f32::from_list_item(text).unwrap().is_nan(), "{text}");
    }
    for text in ["nan(", "nan(x))", "nan(x)tail", "+nan(a b)", "-nan(α)"] {
        assert!(f64::from_list_item(text).is_err(), "{text}");
    }
    for (text, value) in [
        ("inf", f64::INFINITY),
        ("-Infinity", f64::NEG_INFINITY),
        ("+INF", f64::INFINITY),
        ("+-inf", f64::NEG_INFINITY),
    ] {
        assert_eq!(f64::from_list_item(text).unwrap(), value);
    }
}

#[test]
fn source_primitive_formatting_and_metadata_lists() {
    assert_eq!(
        list::to_string_list([0i32, -1, 42]).unwrap(),
        ["0", "-1", "42"]
    );
    assert_eq!(list::to_string_list([true, false]).unwrap(), ["1", "0"]);
    assert_eq!(list::to_string_list(['x', '\0']).unwrap(), ["x", "\0"]);
    assert!(list::to_string_list(['α']).is_err());
    assert_eq!(
        list::to_string_list([
            0.0f64,
            -0.0,
            1.5,
            0.01,
            10000.0,
            0.001,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN
        ])
        .unwrap(),
        [
            "0.0", "-0.0", "1.5", "0.01", "1.0e04", "1.0e-03", "inf", "-inf", "NaN"
        ]
    );
    assert_eq!(
        list::to_string_list([1.2345678f32, 0.01, 10000.0, 1.0e-6]).unwrap(),
        ["1.234568", "0.01", "1.0e04", "1.0e-06"]
    );
    let values = [
        ParamValue::Empty,
        ParamValue::Integer(42),
        ParamValue::Float(1.0),
        ParamValue::FloatList(vec![1.0, 0.001]),
    ];
    assert_eq!(
        list::to_string_list(values).unwrap(),
        ["", "42", "1.0", "[1.0, 1.0e-03]"]
    );
    let values = [
        MetaValue::default(),
        MetaValue::from("text"),
        MetaValue::from(vec![i64::MIN, i64::MAX]),
        MetaValue::new(MetaValueData::FloatList(vec![1.0, 0.001])).unwrap(),
    ];
    assert_eq!(
        list::to_string_list(values).unwrap(),
        [
            "",
            "text",
            "[-9223372036854775808, 9223372036854775807]",
            "[1.0, 1.0e-03]"
        ]
    );
}

#[test]
fn list_bounds_and_custom_conversion_failures_are_checked() {
    let too_many = ",".repeat(list::MAX_ITEMS);
    assert!(list::create::<String>(&too_many, b',').is_err());
    let long_items = vec![""; list::MAX_ITEMS + 1];
    assert!(list::create_from_strings::<i32>(&long_items).is_err());
    struct Failure;
    impl ListFormat for Failure {
        fn to_list_text(&self) -> openms::Result<Cow<'_, str>> {
            Err(openms::Error::InvalidValue(
                "deliberate formatting failure".into(),
            ))
        }
    }
    assert!(list::to_string_list([Failure]).is_err());
    struct ExcessCapacity;
    impl ListFormat for ExcessCapacity {
        fn to_list_text(&self) -> openms::Result<Cow<'_, str>> {
            let mut text = String::with_capacity(list::MAX_BYTES);
            text.push('x');
            Ok(Cow::Owned(text))
        }
    }
    // Account for retained allocation capacity, including custom formatters.
    assert!(list::to_string_list([ExcessCapacity]).is_err());
    let unchanged = vec!["source".to_string()];
    assert!(list::create_from_strings::<i32>(&unchanged).is_err());
    assert_eq!(unchanged, ["source"]);
}
