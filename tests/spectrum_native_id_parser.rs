// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `SpectrumNativeIDParser_test.cpp` (core SDK bc9cc12): one test per
//! `START_SECTION`, cited by source line. Every literal is tier 3, transcribed
//! from the class test; a value derived rather than transcribed says so.

use openms::Error;
use openms::concept::boost_regex::BoostRegex;
use openms::metadata::spectrum_native_id_parser::{
    NATIVE_ID_PREFIXES, extract_scan_number, extract_scan_number_with, is_native_id,
    regex_from_native_id, try_extract_scan_number_with,
};

fn scan(native_id: &str, accession: &str) -> i32 {
    extract_scan_number(native_id, accession)
        .unwrap()
        .source_value()
}

// L24 static bool isNativeID(const std::string& id)
#[test]
fn section_is_native_id() {
    assert!(is_native_id("scan=123"));
    assert!(is_native_id("scanId=456"));
    assert!(is_native_id("scanID=456"));
    assert!(is_native_id("controllerType=0 controllerNumber=1 scan=100"));
    assert!(is_native_id("function=2 process=1 scan=100"));
    assert!(is_native_id("sample=1 period=1 cycle=42 experiment=1"));
    assert!(is_native_id("index=789"));
    assert!(is_native_id("spectrum=101112"));
    assert!(is_native_id("file=42"));

    assert!(!is_native_id("123"));
    assert!(!is_native_id(""));
    assert!(!is_native_id("some_random_string"));
    assert!(!is_native_id("SCAN=123"));

    // Derived: every listed prefix alone is a native ID.
    for prefix in NATIVE_ID_PREFIXES {
        assert!(is_native_id(prefix), "{prefix}");
    }
}

// L46 static std::string getRegExFromNativeID(const std::string& native_id)
#[test]
fn section_regex_from_native_id() {
    let scan = r"scan=(?<GROUP>\d+)";
    assert_eq!(
        regex_from_native_id("controllerType=0 controllerNumber=1 scan=100"),
        scan
    );
    assert_eq!(regex_from_native_id("function=2 process=1 scan=100"), scan);
    assert_eq!(regex_from_native_id("scan=123"), scan);
    assert_eq!(regex_from_native_id("index=456"), r"index=(?<GROUP>\d+)");
    assert_eq!(regex_from_native_id("scanId=789"), r"scanId=(?<GROUP>\d+)");
    assert_eq!(regex_from_native_id("scanID=789"), r"scanID=(?<GROUP>\d+)");
    assert_eq!(
        regex_from_native_id("spectrum=101112"),
        r"spectrum=(?<GROUP>\d+)"
    );
    assert_eq!(regex_from_native_id("file=42"), r"file=(?<GROUP>\d+)");
    assert_eq!(regex_from_native_id("123"), r"(?<GROUP>\d+)");
    assert_eq!(regex_from_native_id("unknown=123"), r"(?<GROUP>\d+)");
    assert_eq!(regex_from_native_id(""), r"(?<GROUP>\d+)");
}

// L82 static Int extractScanNumber(const std::string&, const boost::regex&, bool no_error = false)
#[test]
fn section_extract_scan_number_with_regex() {
    let re_spectrum = BoostRegex::new(r"spectrum=(?<SCAN>\d+)").unwrap();
    assert_eq!(
        extract_scan_number_with("spectrum=42", &re_spectrum).unwrap(),
        42
    );
    assert_eq!(
        extract_scan_number_with("spectrum=0", &re_spectrum).unwrap(),
        0
    );
    assert_eq!(
        extract_scan_number_with("spectrum=99999", &re_spectrum).unwrap(),
        99999
    );

    let re_scan = BoostRegex::new(r"scan=(?<SCAN>\d+)").unwrap();
    assert_eq!(extract_scan_number_with("scan=123", &re_scan).unwrap(), 123);
    assert_eq!(
        extract_scan_number_with("controllerType=0 controllerNumber=1 scan=456", &re_scan).unwrap(),
        456
    );

    // no_error = true: the source's -1 is `None`.
    assert_eq!(
        try_extract_scan_number_with("scan=42", &re_spectrum).unwrap(),
        None
    );
    assert_eq!(
        try_extract_scan_number_with("no_match_here", &re_spectrum).unwrap(),
        None
    );
    assert_eq!(
        try_extract_scan_number_with("", &re_spectrum).unwrap(),
        None
    );

    // no_error = false: the source's Exception::ParseError.
    assert!(matches!(
        extract_scan_number_with("scan=42", &re_spectrum),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        extract_scan_number_with("no_match", &re_spectrum),
        Err(Error::InvalidValue(_))
    ));

    // Several matches: the last one wins.
    let re_multi = BoostRegex::new(r"(?<SCAN>\d+)").unwrap();
    assert_eq!(extract_scan_number_with("1 2 3 42", &re_multi).unwrap(), 42);

    let re_index = BoostRegex::new(r"index=(?<SCAN>\d+)").unwrap();
    assert_eq!(extract_scan_number_with("index=0", &re_index).unwrap(), 0);
}

// L115 static Int extractScanNumber(const std::string&, const std::string& native_id_type_accession)
#[test]
fn section_extract_scan_number_by_accession() {
    // Thermo
    assert_eq!(scan("scan=42", "MS:1000768"), 42);
    assert_eq!(scan("scan=0", "MS:1000768"), 0);
    assert_eq!(
        scan("controllerType=0 controllerNumber=1 scan=123", "MS:1000768"),
        123
    );
    // Waters
    assert_eq!(scan("scan=42", "MS:1000769"), 42);
    assert_eq!(scan("function=1 process=1 scan=99", "MS:1000769"), 99);
    // WIFF: cycle * 1000 + experiment
    assert_eq!(
        scan("sample=1 period=1 cycle=42 experiment=1", "MS:1000770"),
        42001
    );
    assert_eq!(
        scan("sample=1 period=1 cycle=1 experiment=999", "MS:1000770"),
        1999
    );
    assert_eq!(
        scan("sample=1 period=1 cycle=100 experiment=50", "MS:1000770"),
        100050
    );
    // Bruker/Agilent, BAF
    assert_eq!(scan("scan=42", "MS:1000771"), 42);
    assert_eq!(scan("scan=42", "MS:1000772"), 42);
    // Bruker FID
    assert_eq!(scan("file=42", "MS:1000773"), 42);
    assert_eq!(scan("file=0", "MS:1000773"), 0);
    // Index-based: index + 1
    assert_eq!(scan("index=42", "MS:1000774"), 43);
    assert_eq!(scan("index=0", "MS:1000774"), 1);
    // Single peak list
    assert_eq!(scan("file=42", "MS:1000775"), 42);
    // Thermo/Bruker TDF
    assert_eq!(scan("scan=42", "MS:1000776"), 42);
    // Generic spectrum
    assert_eq!(scan("spectrum=42", "MS:1000777"), 42);
    // Agilent MassHunter
    assert_eq!(scan("scanId=42", "MS:1001508"), 42);
    // Numeric only
    assert_eq!(scan("42", "MS:1001530"), 42);
    assert_eq!(scan("0", "MS:1001530"), 0);

    // Unknown accessions.
    assert_eq!(scan("scan=42", "MS:9999999"), -1);
    assert_eq!(scan("scan=42", "invalid"), -1);
    assert_eq!(scan("scan=42", ""), -1);

    // Format mismatch.
    assert_eq!(scan("spectrum=42", "MS:1000768"), -1);
    assert_eq!(scan("", "MS:1000768"), -1);
    assert_eq!(scan("malformed", "MS:1000768"), -1);

    // Merged spectra: the last native ID is the merged scan.
    assert_eq!(
        scan(
            "controllerType=0 controllerNumber=1 scan=100 merged controllerType=0 controllerNumber=1 scan=200",
            "MS:1000768"
        ),
        200
    );
}

/// The warnings the source logs on its three `-1` paths, carried as text
/// (SpectrumNativeIDParser.cpp, `extractScanNumber(native_id, accession)`).
#[test]
fn every_failure_carries_the_sources_warning() {
    let unknown = extract_scan_number("scan=42", "invalid").unwrap();
    assert_eq!(unknown.scan, None);
    assert_eq!(
        unknown.warning.as_deref(),
        Some(
            "native_id: scan=42 accession: invalid Could not extract scan number - no valid native_id_type_accession was provided"
        )
    );
    let no_match = extract_scan_number("malformed", "MS:1000768").unwrap();
    assert_eq!(
        no_match.warning.as_deref(),
        Some("native_id 'malformed' is invalid. Could not extract scan number.")
    );
    // Derived: a digit run beyond i32 fails `toInt32`, the source's
    // ConversionError branch.
    let too_large = extract_scan_number("scan=99999999999", "MS:1000768").unwrap();
    assert_eq!(too_large.scan, None);
    assert_eq!(
        too_large.warning.as_deref(),
        Some(
            "Value: '99999999999' could not be converted to int in string. Native ID='scan=99999999999'"
        )
    );
    let found = extract_scan_number("scan=7", "MS:1000768").unwrap();
    assert_eq!(found.warning, None);
}

/// Where the source throws `Exception::InvalidValue` (experiment of 1000 or
/// more) or overflows `int` (undefined behaviour), derived from the source
/// branches rather than from a class-test literal.
#[test]
fn wiff_and_index_boundaries() {
    assert!(matches!(
        extract_scan_number("sample=1 period=1 cycle=1 experiment=1000", "MS:1000770"),
        Err(Error::InvalidValue(_))
    ));
    // 2147484 * 1000 exceeds i32::MAX; the source's `int` arithmetic
    // overflows, this reports the failure instead.
    let overflow =
        extract_scan_number("sample=1 period=1 cycle=2147484 experiment=1", "MS:1000770").unwrap();
    assert_eq!(overflow.scan, None);
    assert!(overflow.warning.is_some());
    let index_overflow = extract_scan_number("index=2147483647", "MS:1000774").unwrap();
    assert_eq!(index_overflow.scan, None);
    // A WIFF ID with no experiment token has fewer tokens than the two groups.
    assert_eq!(scan("sample=1 period=1 cycle=42", "MS:1000770"), -1);
}

// L178 [EXTRA] Edge cases and error handling
#[test]
fn section_extra_edge_cases() {
    assert_eq!(scan("scan=999999999", "MS:1000768"), 999999999);
    assert_eq!(
        scan("sample=1 period=1 cycle=42   experiment=1", "MS:1000770"),
        42001
    );
    assert!(!is_native_id(" scan=123"));
    assert!(!is_native_id("scan =123"));

    let regex = BoostRegex::new(regex_from_native_id("scan=123")).unwrap();
    assert!(regex.is_search_match(b"scan=456").unwrap());
}

// L199 [EXTRA] Bruker TDF frame= / MS:1002818 native ID
#[test]
fn section_extra_bruker_tdf_frame() {
    assert!(is_native_id("frame=2 scan=529"));
    assert!(is_native_id("frame=10 scan=42 precursor=3"));
    assert_eq!(
        regex_from_native_id("frame=2 scan=529"),
        r"scan=(?<GROUP>\d+)"
    );
    assert_eq!(
        regex_from_native_id("frame=10 scan=42 precursor=3"),
        r"scan=(?<GROUP>\d+)"
    );
    assert_eq!(scan("frame=2 scan=529", "MS:1002818"), 529);
    assert_eq!(scan("frame=10 scan=42 precursor=3", "MS:1002818"), 42);
    assert_eq!(scan("frame=1 scan=5 frame=1 scan=9", "MS:1002818"), 9);
}
