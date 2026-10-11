// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `USI_test.cpp` (core SDK bc9cc12): one test per `START_SECTION`,
//! cited by source line. Every literal is tier 3, transcribed from the class
//! test; a value derived rather than transcribed says so.

use openms::Error;
use openms::metadata::usi::{CV_ACCESSION, CV_NAME, IndexType, Usi};

fn usi(collection: &str, ms_run: &str, index_type: IndexType, index: &str) -> Usi {
    Usi::new(collection, ms_run, index_type, index, "")
}

// L29 USI(), L37 ~USI()
#[test]
fn section_default_constructor() {
    let empty = Usi::default();
    assert!(!empty.is_valid());
    // Derived from the source's member initialisers.
    assert_eq!(empty.index_type, IndexType::Scan);
}

// L43 USI(collection, ms_run, index_type, index, interpretation)
#[test]
fn section_component_constructor() {
    let u = Usi::new(
        "PXD000561",
        "sample.mzML",
        IndexType::Scan,
        "12345",
        "PEPTIDEK/2",
    );
    assert!(u.is_valid());
    assert_eq!(u.collection, "PXD000561");
    assert_eq!(u.ms_run, "sample.mzML");
    assert_eq!(u.index_type, IndexType::Scan);
    assert_eq!(u.index, "12345");
    assert_eq!(u.interpretation, "PEPTIDEK/2");
}

// L55 explicit USI(const std::string& usi_string)
#[test]
fn section_string_constructor() {
    let u = Usi::parse("mzspec:PXD000561:sample.mzML:scan:12345:PEPTIDEK/2").unwrap();
    assert!(u.is_valid());
    assert_eq!(u.collection, "PXD000561");
    assert_eq!(u.ms_run, "sample.mzML");
    assert_eq!(u.index_type, IndexType::Scan);
    assert_eq!(u.index, "12345");
    assert_eq!(u.interpretation, "PEPTIDEK/2");

    let u3 = Usi::parse("mzspec:PXD000561:sample.mzML:scan:12345:EM[UNIMOD:35]K/2").unwrap();
    assert!(u3.is_valid());
    assert_eq!(u3.interpretation, "EM[UNIMOD:35]K/2");
    assert_eq!(
        u3.to_usi_string(),
        "mzspec:PXD000561:sample.mzML:scan:12345:EM[UNIMOD:35]K/2"
    );

    let u4 = Usi::parse("mzspec:PXD000561:sample.mzML:scan:12345").unwrap();
    assert!(u4.is_valid());
    assert_eq!(u4.interpretation, "");

    // Exception::ParseError
    assert!(matches!(
        Usi::parse("invalid_usi_string"),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        Usi::parse("mzspec:PXD000561"),
        Err(Error::InvalidValue(_))
    ));
    // `FromStr` is the same constructor.
    assert_eq!(
        "mzspec:PXD000561:sample.mzML:scan:12345"
            .parse::<Usi>()
            .unwrap(),
        u4
    );
}

// L83 USI(const USI&), L91 USI(USI&&), L99 operator=(const USI&), L108 operator=(USI&&)
#[test]
fn section_copy_and_move() {
    let usi1 = usi("PXD000561", "sample.mzML", IndexType::Scan, "12345");
    let usi2 = usi1.clone();
    assert_eq!(usi1, usi2);
    let mut usi3 = Usi::default();
    assert_ne!(usi3, usi1);
    usi3 = usi1.clone();
    assert_eq!(usi3, usi1);
    let moved = usi1;
    assert_eq!(moved.collection, "PXD000561");
}

// L117 operator==, L128 operator!=
#[test]
fn section_equality() {
    let usi1 = usi("PXD000561", "sample.mzML", IndexType::Scan, "12345");
    let usi2 = usi("PXD000561", "sample.mzML", IndexType::Scan, "12345");
    let usi3 = usi("PXD000562", "sample.mzML", IndexType::Scan, "12345");
    assert!(usi1 == usi2);
    assert!(usi1 != usi3);
}

// L137 operator<
#[test]
fn section_less_than() {
    let usi1 = usi("PXD000561", "sample.mzML", IndexType::Scan, "12345");
    let usi2 = usi("PXD000562", "sample.mzML", IndexType::Scan, "12345");
    assert!(usi1 < usi2);

    // Derived from the source's field-by-field operator<: the index type
    // compares in declaration order and only after collection and MS run,
    // and the interpretation decides last.
    let scan = usi("A", "r", IndexType::Scan, "9");
    let index = usi("A", "r", IndexType::Index, "1");
    let native = usi("A", "r", IndexType::NativeId, "0");
    assert!(scan < index && index < native);
    assert!(usi("A", "s", IndexType::Scan, "0") > native);
    let plain = usi("A", "r", IndexType::Scan, "9");
    let interpreted = Usi::new("A", "r", IndexType::Scan, "9", "P/1");
    assert!(plain < interpreted);
}

// L146 isValid
#[test]
fn section_is_valid() {
    assert!(!Usi::default().is_valid());
    assert!(usi("PXD000561", "sample.mzML", IndexType::Scan, "12345").is_valid());
    assert!(!usi("", "sample.mzML", IndexType::Scan, "12345").is_valid());
}

// L159 static isValidUSI
#[test]
fn section_is_valid_usi() {
    assert!(Usi::is_valid_usi("mzspec:PXD000561:sample.mzML:scan:12345"));
    assert!(Usi::is_valid_usi(
        "mzspec:PXD000561:sample.mzML:scan:12345:PEPTIDEK/2"
    ));
    assert!(Usi::is_valid_usi(
        "mzspec:PXD000561:sample.mzML:scan:12345:EM[UNIMOD:35]K/2"
    ));
    assert!(!Usi::is_valid_usi("invalid_string"));
    assert!(!Usi::is_valid_usi("mzspec:PXD000561"));
    assert!(!Usi::is_valid_usi(""));
}

// L170-L243: the ten accessor sections.
#[test]
fn section_accessors() {
    let u = usi("PXD000561", "sample.mzML", IndexType::Scan, "12345");
    assert_eq!(u.collection, "PXD000561");
    assert_eq!(u.ms_run, "sample.mzML");
    assert_eq!(u.index, "12345");

    let mut u = Usi::default();
    assert!(!u.is_valid());
    u.collection = "PXD000562".into();
    assert_eq!(u.collection, "PXD000562");
    u.ms_run = "test.mzML".into();
    assert_eq!(u.ms_run, "test.mzML");
    u.index_type = IndexType::NativeId;
    assert_eq!(u.index_type, IndexType::NativeId);
    u.index = "99999".into();
    assert_eq!(u.index, "99999");
    u.interpretation = "SEQUENCE/3".into();
    assert_eq!(u.interpretation, "SEQUENCE/3");

    let indexed = usi("PXD000561", "sample.mzML", IndexType::Index, "12345");
    assert_eq!(indexed.index_type, IndexType::Index);
    let interpreted = Usi::new(
        "PXD000561",
        "sample.mzML",
        IndexType::Scan,
        "12345",
        "PEPTIDEK/2",
    );
    assert_eq!(interpreted.interpretation, "PEPTIDEK/2");
}

// L245 hasInterpretation
#[test]
fn section_has_interpretation() {
    let usi1 = Usi::new(
        "PXD000561",
        "sample.mzML",
        IndexType::Scan,
        "12345",
        "PEPTIDEK/2",
    );
    assert!(usi1.has_interpretation());
    let usi2 = usi("PXD000561", "sample.mzML", IndexType::Scan, "12345");
    assert!(!usi2.has_interpretation());
}

// L255 toString
#[test]
fn section_to_string() {
    let usi1 = Usi::new(
        "PXD000561",
        "sample.mzML",
        IndexType::Scan,
        "12345",
        "PEPTIDEK/2",
    );
    assert_eq!(
        usi1.to_usi_string(),
        "mzspec:PXD000561:sample.mzML:scan:12345:PEPTIDEK/2"
    );
    let usi2 = usi("PXD000561", "sample.mzML", IndexType::Scan, "12345");
    assert_eq!(
        usi2.to_usi_string(),
        "mzspec:PXD000561:sample.mzML:scan:12345"
    );
    let usi3 = usi("PXD000561", "sample.mzML", IndexType::Index, "0");
    assert_eq!(usi3.to_usi_string(), "mzspec:PXD000561:sample.mzML:index:0");
    assert_eq!(Usi::default().to_usi_string(), "");
}

// L271 fromString
#[test]
fn section_from_string() {
    let mut u = Usi::default();
    assert!(u.set_from_string("mzspec:PXD000561:sample.mzML:scan:12345:PEPTIDEK/2"));
    assert_eq!(u.collection, "PXD000561");
    assert_eq!(u.ms_run, "sample.mzML");
    assert_eq!(u.index_type, IndexType::Scan);
    assert_eq!(u.index, "12345");
    assert_eq!(u.interpretation, "PEPTIDEK/2");

    assert!(u.set_from_string("mzspec:PXD000561:sample.mzML:index:0"));
    assert_eq!(u.index_type, IndexType::Index);

    assert!(!u.set_from_string("invalid"));
    assert!(!u.set_from_string("mzspec:PXD000561"));
    assert!(!u.set_from_string(""));
}

/// Derived from `USI::fromString`: fields are reset, then assigned as read,
/// so a rejection part-way leaves what was read; a missing field is invalid;
/// a trailing colon is an empty interpretation.
#[test]
fn from_string_boundaries() {
    let mut u = Usi::new("X", "Y", IndexType::NativeId, "Z", "W");
    assert!(!u.set_from_string("mzspec:PXD1:run.mzML:bogus:7"));
    assert_eq!(u.collection, "PXD1");
    assert_eq!(u.ms_run, "run.mzML");
    assert_eq!(u.index_type, IndexType::Scan);
    assert_eq!(u.index, "");
    assert_eq!(u.interpretation, "");

    assert!(!u.set_from_string("mzspec:PXD1::scan:7"));
    assert!(u.set_from_string("mzspec:PXD1:run.mzML:NATIVEID:controllerType=0 scan=3:"));
    assert_eq!(u.index_type, IndexType::NativeId);
    assert_eq!(u.index, "controllerType=0 scan=3");
    assert!(!u.has_interpretation());
    // The index type is written back in its canonical spelling.
    assert_eq!(
        u.to_usi_string(),
        "mzspec:PXD1:run.mzML:nativeId:controllerType=0 scan=3"
    );
}

// L293 createFromScanNumber
#[test]
fn section_create_from_scan_number() {
    let usi1 = Usi::from_scan_number("PXD000561", "sample.mzML", 12345, "");
    assert_eq!(
        usi1.to_usi_string(),
        "mzspec:PXD000561:sample.mzML:scan:12345"
    );
    let usi2 = Usi::from_scan_number("PXD000561", "sample.mzML", 12345, "PEPTIDEK/2");
    assert_eq!(
        usi2.to_usi_string(),
        "mzspec:PXD000561:sample.mzML:scan:12345:PEPTIDEK/2"
    );
    let usi3 = Usi::from_scan_number("PXD000561", "sample.mzML", 12345, "EM[UNIMOD:35]K/2");
    assert_eq!(
        usi3.to_usi_string(),
        "mzspec:PXD000561:sample.mzML:scan:12345:EM[UNIMOD:35]K/2"
    );
}

// L306 createFromNativeID
#[test]
fn section_create_from_native_id() {
    let usi1 = Usi::from_native_id("PXD000561", "sample.mzML", "scan=12345").unwrap();
    assert_eq!(
        usi1.to_usi_string(),
        "mzspec:PXD000561:sample.mzML:scan:12345"
    );
    let usi2 = Usi::from_native_id(
        "PXD000561",
        "sample.mzML",
        "controllerType=0 controllerNumber=1 scan=54321",
    )
    .unwrap();
    assert_eq!(usi2.index, "54321");
    let usi3 = Usi::from_native_id("PXD000561", "sample.mzML", "custom_id_format").unwrap();
    assert_eq!(usi3.index_type, IndexType::NativeId);
    assert_eq!(usi3.index, "custom_id_format");
}

// L322 extractScanNumberFromNativeID
#[test]
fn section_extract_scan_number_from_native_id() {
    let extract = |id: &str| Usi::extract_scan_number_from_native_id(id).unwrap();
    assert_eq!(extract("scan=12345"), Some(12345));
    assert_eq!(
        extract("controllerType=0 controllerNumber=1 scan=54321"),
        Some(54321)
    );
    assert_eq!(extract("scanId=99999"), Some(99999));
    assert_eq!(extract("index=100"), Some(100));
    assert_eq!(extract("12345"), Some(12345));
    assert_eq!(extract("custom_format_no_number"), None);

    // Derived from the source: a WIFF ID has no prefix-specific expression,
    // so the generic `(?<GROUP>\d+)` picks the last number, the experiment.
    assert_eq!(extract("sample=1 period=1 cycle=42 experiment=7"), Some(7));
    // A native-ID prefix without a number falls through to toInt32, which
    // fails on the whole string.
    assert_eq!(extract("scan=abc"), None);
}

// L350 indexTypeToString, L358 indexTypeFromString
#[test]
fn section_index_type_strings() {
    assert_eq!(IndexType::Scan.as_str(), "scan");
    assert_eq!(IndexType::Index.as_str(), "index");
    assert_eq!(IndexType::NativeId.as_str(), "nativeId");
    assert_eq!(IndexType::NativeId.to_string(), "nativeId");

    assert_eq!(IndexType::parse("scan").unwrap(), IndexType::Scan);
    assert_eq!(IndexType::parse("SCAN").unwrap(), IndexType::Scan);
    assert_eq!(IndexType::parse("index").unwrap(), IndexType::Index);
    assert_eq!(IndexType::parse("INDEX").unwrap(), IndexType::Index);
    assert_eq!(IndexType::parse("nativeId").unwrap(), IndexType::NativeId);
    assert_eq!(IndexType::parse("NATIVEID").unwrap(), IndexType::NativeId);
    assert!(matches!(
        IndexType::parse("invalid_type"),
        Err(Error::InvalidValue(_))
    ));
}

// L371 extractBasename
#[test]
fn section_extract_basename() {
    assert_eq!(Usi::extract_basename("/path/to/sample.mzML"), "sample.mzML");
    assert_eq!(
        Usi::extract_basename("/data/experiments/run1.mzML"),
        "run1.mzML"
    );
    assert_eq!(Usi::extract_basename("/sample.mzML"), "sample.mzML");
    assert_eq!(
        Usi::extract_basename("C:\\data\\sample.mzML"),
        "sample.mzML"
    );
    assert_eq!(
        Usi::extract_basename("C:\\Users\\test\\data\\run1.mzML"),
        "run1.mzML"
    );
    assert_eq!(
        Usi::extract_basename("file:///path/to/sample.mzML"),
        "sample.mzML"
    );
    assert_eq!(
        Usi::extract_basename("file:///C:/Users/bielow/data/sample.mzML"),
        "sample.mzML"
    );
    assert_eq!(
        Usi::extract_basename("file:///C:/data/experiments/ES-0014b_2.mzML"),
        "ES-0014b_2.mzML"
    );
    assert_eq!(Usi::extract_basename("sample.mzML"), "sample.mzML");
    assert_eq!(Usi::extract_basename("run1.raw"), "run1.raw");
    assert_eq!(Usi::extract_basename(""), "");
    assert_eq!(
        Usi::extract_basename("/data/sub\\folder/sample.mzML"),
        "sample.mzML"
    );
}

// L399 getCVAccession, L405 getCVName
#[test]
fn section_cv_term() {
    assert_eq!(CV_ACCESSION, "MS:1003063");
    assert_eq!(CV_NAME, "universal spectrum identifier");
}

// L411 operator<<
#[test]
fn section_stream_operator() {
    let u = usi("PXD000561", "sample.mzML", IndexType::Scan, "12345");
    assert_eq!(format!("{u}"), "mzspec:PXD000561:sample.mzML:scan:12345");
}

// L421 [EXTRA] USI construction with file paths
#[test]
fn section_extra_file_paths() {
    let usi1 = Usi::from_scan_number("PXD000561", "/path/to/sample.mzML", 12345, "");
    assert_eq!(usi1.ms_run, "/path/to/sample.mzML");

    let full_path = "file:///C:/Users/bielow/AppData/Local/Temp/20190911_110348_8204_1/Untitled_workflow/002_FileFilter/out/ES-0014b_2.mzML";
    let basename = Usi::extract_basename(full_path);
    let usi2 = Usi::from_scan_number("PXD000561", &basename, 219, "");
    assert_eq!(usi2.ms_run, "ES-0014b_2.mzML");
    assert_eq!(
        usi2.to_usi_string(),
        "mzspec:PXD000561:ES-0014b_2.mzML:scan:219"
    );

    let win_path = "C:\\data\\experiments\\sample.mzML";
    let usi3 = Usi::from_scan_number("PXD000561", &Usi::extract_basename(win_path), 100, "");
    assert_eq!(usi3.ms_run, "sample.mzML");
}
