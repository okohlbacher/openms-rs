// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `SpectrumLookup_test.cpp` (core SDK bc9cc12): one test per
//! `START_SECTION`, cited by source line. The class test shares one `lookup`
//! across sections; [`read`] rebuilds that state. `Exception::ElementNotFound`
//! is `Error::InvalidRange` here, `IllegalArgument` and `ParseError` are
//! `Error::InvalidValue`. Literals are tier 3; derived values say so.

use openms::concept::boost_regex::BoostRegex;
use openms::metadata::spectrum_lookup::{DEFAULT_SCAN_REGEXP, SpectrumLookup};
use openms::{Error, MSSpectrum};

/// The class test's three spectra (SpectrumLookup_test.cpp:41-50).
fn spectra() -> Vec<MSSpectrum> {
    (0..3)
        .map(|i| MSSpectrum {
            native_id: format!("spectrum={i}"),
            rt: f64::from(i + 1),
            ..MSSpectrum::default()
        })
        .collect()
}

fn read() -> SpectrumLookup {
    let mut lookup = SpectrumLookup::new();
    let warnings = lookup
        .read_spectra(&spectra(), DEFAULT_SCAN_REGEXP)
        .unwrap();
    assert!(warnings.is_empty());
    lookup
}

// L27 SpectrumLookup(), L35 ~SpectrumLookup()
#[test]
fn section_constructor() {
    let lookup = SpectrumLookup::new();
    assert_eq!(lookup.rt_tolerance, 0.01);
}

// L55 bool empty() const
#[test]
fn section_empty() {
    assert!(SpectrumLookup::new().is_empty());
}

// L62 readSpectra
#[test]
fn section_read_spectra() {
    assert!(!read().is_empty());
}

// L70 Size findByRT(double) const
#[test]
fn section_find_by_rt() {
    let lookup = read();
    assert_eq!(lookup.find_by_rt(2.0).unwrap(), 1);
    assert!(matches!(
        lookup.find_by_rt(5.0),
        Err(Error::InvalidRange(_))
    ));
}

// L79 Size findByNativeID(const std::string&) const
#[test]
fn section_find_by_native_id() {
    let lookup = read();
    assert_eq!(lookup.find_by_native_id("spectrum=1").unwrap(), 1);
    assert!(matches!(
        lookup.find_by_native_id("spectrum=3"),
        Err(Error::InvalidRange(_))
    ));
}

// L89 Size findByIndex(Size, bool) const
#[test]
fn section_find_by_index() {
    let lookup = read();
    assert_eq!(lookup.find_by_index(1, false).unwrap(), 1);
    assert_eq!(lookup.find_by_index(1, true).unwrap(), 0);
    assert!(matches!(
        lookup.find_by_index(0, true),
        Err(Error::InvalidRange(_))
    ));
}

// L99 Size findByScanNumber(Size) const
#[test]
fn section_find_by_scan_number() {
    let lookup = read();
    assert_eq!(lookup.find_by_scan_number(1).unwrap(), 1);
    assert!(matches!(
        lookup.find_by_scan_number(5),
        Err(Error::InvalidRange(_))
    ));
}

// L108 void addReferenceFormat(const std::string&), L119 Size findByReference(const std::string&) const
#[test]
fn section_add_reference_format_and_find_by_reference() {
    let mut lookup = read();
    assert!(matches!(
        lookup.add_reference_format("XXX"),
        Err(Error::InvalidValue(_))
    ));
    assert!(lookup.reference_formats.is_empty());
    lookup
        .add_reference_format(r"scan_number=(?<SCAN>\d+)")
        .unwrap();
    lookup.add_reference_format(r"(?<ID>spectrum=\d+)").unwrap();

    assert_eq!(lookup.find_by_reference("scan_number=1").unwrap(), 1);
    assert_eq!(lookup.find_by_reference("name=bla,spectrum=0").unwrap(), 0);
    assert!(matches!(
        lookup.find_by_reference("test123"),
        Err(Error::InvalidValue(_))
    ));
}

// L129 static Int extractScanNumber(const std::string&, const boost::regex&) (deprecated forwarder)
#[test]
#[allow(deprecated)]
fn section_deprecated_extract_scan_number_with_regex() {
    let re = BoostRegex::new(r"spectrum=(?<SCAN>\d+)").unwrap();
    assert_eq!(
        SpectrumLookup::extract_scan_number_with("spectrum=42", &re).unwrap(),
        42
    );
    assert_eq!(
        SpectrumLookup::try_extract_scan_number_with("scan=42", &re).unwrap(),
        None
    );
    assert!(matches!(
        SpectrumLookup::extract_scan_number_with("scan=42", &re),
        Err(Error::InvalidValue(_))
    ));
    assert!(SpectrumLookup::is_native_id("scan=42"));
    assert_eq!(
        SpectrumLookup::regex_from_native_id("scan=42"),
        r"scan=(?<GROUP>\d+)"
    );
}

// L140 static Int extractScanNumber(const std::string&, const std::string&) (deprecated forwarder)
#[test]
#[allow(deprecated)]
fn section_deprecated_extract_scan_number_by_accession() {
    let scan = |id: &str, accession: &str| {
        SpectrumLookup::extract_scan_number(id, accession)
            .unwrap()
            .source_value()
    };
    assert_eq!(scan("scan=42", "MS:1000768"), 42);
    assert_eq!(scan("scan=42", "MS:1000769"), 42);
    assert_eq!(scan("scan=42", "MS:1000771"), 42);
    assert_eq!(scan("scan=42", "MS:1000772"), 42);
    assert_eq!(scan("scan=42", "MS:1000776"), 42);
    assert_eq!(
        scan("sample=1 period=1 cycle=42 experiment=1", "MS:1000770"),
        42001
    );
    assert_eq!(scan("file=42", "MS:1000773"), 42);
    assert_eq!(scan("file=42", "MS:1000775"), 42);
    assert_eq!(scan("index=42", "MS:1000774"), 43);
    assert_eq!(scan("scanId=42", "MS:1001508"), 42);
    assert_eq!(scan("spectrum=42", "MS:1000777"), 42);
    assert_eq!(scan("42", "MS:1001530"), 42);
}

/// Derived from `SpectrumLookup.cpp`: the RT tie rule, the tolerance edge,
/// last-wins duplicates, the scan-regexp check, reference groups in their
/// fixed order, and a matched format without usable information.
#[test]
fn source_rules_beyond_the_class_test() {
    let mut duplicated = spectra();
    duplicated.push(MSSpectrum {
        native_id: "spectrum=1".into(),
        rt: 2.0,
        ..MSSpectrum::default()
    });
    let mut lookup = SpectrumLookup::new();
    lookup
        .read_spectra(&duplicated, DEFAULT_SCAN_REGEXP)
        .unwrap();
    // Duplicate RT and native ID: the later spectrum wins.
    assert_eq!(lookup.find_by_rt(2.0).unwrap(), 3);
    assert_eq!(lookup.find_by_native_id("spectrum=1").unwrap(), 3);
    assert_eq!(lookup.find_by_scan_number(1).unwrap(), 3);
    // Halfway between 1.0 and 2.0 is outside 0.01; a wider tolerance then
    // prefers the upper neighbour on a tie.
    assert!(lookup.find_by_rt(1.5).is_err());
    lookup.rt_tolerance = 0.5;
    assert_eq!(lookup.find_by_rt(1.5).unwrap(), 3);
    assert_eq!(lookup.find_by_rt(1.4).unwrap(), 0);
    assert!(lookup.find_by_rt(f64::NAN).is_err());

    // A scan expression must name SCAN; the old index survives the refusal.
    assert!(matches!(
        lookup.read_spectra(&duplicated, r"=(\d+)$"),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(lookup.find_by_rt(2.0).unwrap(), 3);
    // An empty expression disables scan numbers without a warning.
    let warnings = lookup.read_spectra(&duplicated, "").unwrap();
    assert!(warnings.is_empty());
    assert!(lookup.find_by_scan_number(1).is_err());
    // A native ID the expression cannot parse warns.
    let odd = [MSSpectrum {
        native_id: "controllerType=0 scan".into(),
        rt: 1.0,
        ..MSSpectrum::default()
    }];
    let warnings = lookup.read_spectra(&odd, DEFAULT_SCAN_REGEXP).unwrap();
    assert_eq!(warnings.len(), 1);

    let mut lookup = SpectrumLookup::new();
    lookup
        .read_spectra(&spectra(), DEFAULT_SCAN_REGEXP)
        .unwrap();
    lookup
        .add_reference_format(r"idx=(?<INDEX1>\d*),rt=(?<RT>[\d.]+)")
        .unwrap();
    // INDEX1 is consulted before RT; an empty INDEX1 falls through to RT.
    assert_eq!(lookup.find_by_reference("idx=3,rt=1.0").unwrap(), 2);
    assert_eq!(lookup.find_by_reference("idx=,rt=1.0").unwrap(), 0);
    // A matched format that captures nothing usable.
    let mut empty = SpectrumLookup::new();
    empty.read_spectra(&spectra(), DEFAULT_SCAN_REGEXP).unwrap();
    empty.add_reference_format(r"x(?<ID>\d*)").unwrap();
    assert!(matches!(
        empty.find_by_reference("x"),
        Err(Error::MissingInformation(_))
    ));
}
