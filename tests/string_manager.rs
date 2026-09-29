// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! StringManager and XMLAttributes.
//!
//! The StringManager class test is `XMLHandler_test.cpp` at core `bc9cc12`
//! (its `START_TEST` names StringManager); the UTF-8 section is from
//! `MzMLFile_test.cpp` ("XML metadata uses UTF-8 and a strict ASCII fast
//! path"). Their literals are the expected values. The `compress64` section
//! tests a private function and lives in `src/format/string_manager.rs`.
//! `parseInt` has no upstream test; its cases are Xerces' documented rules
//! (trim, then `strtol` base 10 over the whole rest, then an `(int)` cast).

use openms::Error;
use openms::format::string_manager::StringManager;
use openms::format::xml_attributes::XMLAttributes;

const RUSSIAN_HELLO: [u16; 10] = [
    0x041F, 0x0440, 0x0438, 0x0432, 0x0435, 0x0442, 0x043C, 0x0438, 0x0440, 0x0000,
];
const ASCII: [u16; 13] = [
    0x0048, 0x0065, 0x006C, 0x006C, 0x006F, 0x002C, 0x0057, 0x006F, 0x0072, 0x006C, 0x0064, 0x0021,
    0x0000,
];
const MIXED: [u16; 12] = [
    0x0048, 0x0065, 0x0432, 0x0435, 0x0442, 0x043C, 0x006F, 0x0072, 0x006C, 0x0064, 0x0021, 0x0000,
];
const EMPTY: [u16; 1] = [0];
const UPPER_BOUNDARY: [u16; 3] = [0x00FF, 0x00FF, 0x0000];

fn u(text: &str) -> Vec<u16> {
    text.encode_utf16().collect()
}

#[test]
fn class_test_str_length() {
    // The upstream section compares strLength with itself; the lengths of
    // its arrays up to the terminator are 12, 0, 2 and 9.
    assert_eq!(StringManager::str_length(&ASCII), 12);
    assert_eq!(StringManager::str_length(&EMPTY), 0);
    assert_eq!(StringManager::str_length(&UPPER_BOUNDARY), 2);
    assert_eq!(StringManager::str_length(&RUSSIAN_HELLO), 9);
    assert_eq!(StringManager::str_length(&MIXED), 11);
    // A null pointer has length 0; so has an empty slice. An unterminated
    // slice ends at its end.
    assert_eq!(StringManager::str_length(&[]), 0);
    assert_eq!(StringManager::str_length(&ASCII[..5]), 5);
}

#[test]
fn class_test_is_ascii() {
    let length = |units: &[u16]| StringManager::str_length(units);
    assert!(StringManager::is_ascii(&ASCII[..length(&ASCII)]));
    assert!(!StringManager::is_ascii(
        &RUSSIAN_HELLO[..length(&RUSSIAN_HELLO)]
    ));
    assert!(!StringManager::is_ascii(&MIXED[..length(&MIXED)]));
    assert!(StringManager::is_ascii(&EMPTY[..length(&EMPTY)]));
    // U+00FF is Latin-1, not ASCII.
    assert!(!StringManager::is_ascii(
        &UPPER_BOUNDARY[..length(&UPPER_BOUNDARY)]
    ));
    // MzMLFile_test.cpp
    assert!(StringManager::is_ascii(&u("ASCII")));
    assert!(!StringManager::is_ascii(&u("M\u{fc}ller")));
    assert!(!StringManager::is_ascii(&u("1234567\u{fc}")));
}

#[test]
fn class_test_append_ascii() {
    let ascii = &ASCII[..12];
    let mut o5 = String::new();
    StringManager::append_ascii(ascii, &mut o5).unwrap();
    assert_eq!(o5, "Hello,World!");
    let mut o6 = String::from("Gruess Gott und ");
    StringManager::append_ascii(ascii, &mut o6).unwrap();
    assert_eq!(o6, "Gruess Gott und Hello,World!");
    let mut o7 = String::new();
    StringManager::append_ascii(&EMPTY[..0], &mut o7).unwrap();
    assert_eq!(o7, "");
    // Native: a non-ASCII unit is refused and the string left alone.
    let mut kept = String::from("kept");
    assert!(matches!(
        StringManager::append_ascii(&RUSSIAN_HELLO[..9], &mut kept),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(kept, "kept");
}

#[test]
fn utf8_conversion_both_ways() {
    // MzMLFile_test.cpp: the expected UTF-8 bytes spelled out.
    let text = StringManager::convert_utf16(&u("M\u{fc}ller \u{3b1} \u{1f9ea}")).unwrap();
    assert_eq!(
        text.as_bytes(),
        b"M\xC3\xBCller \xCE\xB1 \xF0\x9F\xA7\xAA".as_slice()
    );
    assert_eq!(
        StringManager::convert_utf16(&RUSSIAN_HELLO).unwrap(),
        "\u{41f}\u{440}\u{438}\u{432}\u{435}\u{442}\u{43c}\u{438}\u{440}"
    );
    assert_eq!(
        StringManager::convert_utf16(&ASCII).unwrap(),
        "Hello,World!"
    );
    // Conversion stops at the terminator.
    assert_eq!(StringManager::convert_utf16(&[0x41, 0, 0x42]).unwrap(), "A");
    assert!(StringManager::convert_utf16(&[0xD800, 0x41]).is_err());
    // Round trip, with a surrogate pair.
    let wide = StringManager::convert("M\u{fc}ller \u{1f9ea}");
    assert_eq!(wide, u("M\u{fc}ller \u{1f9ea}"));
    // "M\u{fc}ller " is seven units, the test tube a surrogate pair.
    assert_eq!(wide.len(), 9);
    assert_eq!(
        StringManager::convert_utf16(&wide).unwrap(),
        "M\u{fc}ller \u{1f9ea}"
    );
    let terminated = StringManager::convert_ptr("ab");
    assert_eq!(terminated, vec![0x61, 0x62, 0]);
    assert_eq!(StringManager::new(), StringManager);
}

#[test]
fn parse_int_follows_xerces() {
    assert_eq!(StringManager::parse_int(&u("42")).unwrap(), 42);
    assert_eq!(StringManager::parse_int(&u("  -17\t\n")).unwrap(), -17);
    assert_eq!(StringManager::parse_int(&u("+0")).unwrap(), 0);
    assert_eq!(StringManager::parse_int(&u("7\u{0}junk")).unwrap(), 7);
    for bad in ["", "   ", "4 2", "1.5", "1e3", "abc", "-", "\u{661}"] {
        assert!(
            matches!(
                StringManager::parse_int(&u(bad)),
                Err(Error::InvalidValue(_))
            ),
            "{bad:?}"
        );
    }
    assert_eq!(
        StringManager::parse_int(&u("-2147483648")).unwrap(),
        i32::MIN
    );
    // Beyond i32: refused natively, truncated by the source's (int) cast.
    assert!(matches!(
        StringManager::parse_int(&u("2147483648")),
        Err(Error::InvalidRange(_))
    ));
    assert_eq!(
        StringManager::parse_int_source(&u("2147483648")).unwrap(),
        i32::MIN
    );
    assert_eq!(
        StringManager::parse_int_source(&u("4294967296")).unwrap(),
        0
    );
    assert!(matches!(
        StringManager::parse_int_source(&u("99999999999999999999")),
        Err(Error::InvalidRange(_))
    ));
}

#[test]
fn xml_attributes_view() {
    let list = vec![
        (
            "xmlns:xsi".to_owned(),
            "http://www.w3.org/2001/XMLSchema-instance".to_owned(),
        ),
        ("id".to_owned(), "run1".to_owned()),
        ("id2".to_owned(), "".to_owned()),
    ];
    let attributes = XMLAttributes::new(&list);
    let copy = attributes;
    assert_eq!(copy, attributes);
    assert_eq!(attributes.len(), 3);
    assert!(!attributes.is_empty());
    assert_eq!(attributes.value("id"), Some("run1"));
    assert_eq!(attributes.value("id2"), Some(""));
    assert_eq!(attributes.value("ID"), None);
    assert_eq!(attributes.value("xsi"), None);
    assert_eq!(
        attributes.value_utf16(&StringManager::convert_ptr("xmlns:xsi")),
        Some("http://www.w3.org/2001/XMLSchema-instance")
    );
    assert_eq!(attributes.index("id2"), Some(2));
    assert_eq!(attributes.index_utf16(&u("id")), Some(1));
    assert_eq!(attributes.index("absent"), None);
    assert_eq!(attributes.index_utf16(&[0xD800]), None);
    assert_eq!(attributes.value_by_index(0).map(str::len), Some(41));
    assert_eq!(attributes.value_by_index(3), None);
    assert_eq!(attributes.handle().len(), 3);
    assert!(XMLAttributes::new(&[]).is_empty());
}
