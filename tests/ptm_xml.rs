// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FORMAT/PTMXMLFile.h` and `FORMAT/HANDLERS/PTMXMLHandler.h`.
//!
//! Expectations come from `PTMXMLFile_test.cpp` and its unchanged fixture
//! `PTMs.xml` (here `tests/data/PTMXMLFile_PTMs.xml`); the byte layout of the
//! writer comes from `PTMXMLHandler::writeTo` (`PTMXMLHandler.cpp:27-39`).
//! See `tests/data/ptm_xml_provenance.json`.

#![cfg(feature = "kb-xml")]

use openms::Error;
use openms::format::ptm_xml::{self, PtmInformation, ReadLimits};
use std::path::PathBuf;

fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/PTMXMLFile_PTMs.xml")
}

fn temporary(name: &str) -> PathBuf {
    let directory = std::env::temp_dir().join(format!(
        "openms-ptm-xml-{}-{name}",
        std::process::id()
    ));
    std::fs::create_dir_all(&directory).unwrap();
    directory.join(name)
}

/// `PTMXMLFile_test.cpp`, section `load`.
#[test]
fn load_reads_the_class_test_fixture() {
    let ptm_informations = ptm_xml::load(fixture()).unwrap();
    assert_eq!(ptm_informations["TEST"].0, "N2O2-CH3");
    assert_eq!(ptm_informations["TEST"].1, "KLR");
    assert_eq!(ptm_informations.len(), 1);
}

/// `PTMXMLFile_test.cpp`, section `store`: store, clear, load again.
#[test]
fn store_round_trips_the_class_test_fixture() {
    let ptm_informations = ptm_xml::load(fixture()).unwrap();
    let path = temporary("store.xml");
    ptm_xml::store(&path, &ptm_informations).unwrap();
    let reloaded = ptm_xml::load(&path).unwrap();
    assert_eq!(reloaded["TEST"].0, "N2O2-CH3");
    assert_eq!(reloaded["TEST"].1, "KLR");
    // The upstream fixture is laid out exactly as writeTo writes it: no
    // declaration, tab indentation, one line per element. So storing what it
    // holds reproduces it byte for byte.
    assert_eq!(
        std::fs::read(&path).unwrap(),
        std::fs::read(fixture()).unwrap()
    );
    std::fs::remove_file(&path).unwrap();
}

/// The writeTo layout for several entries, from `PTMXMLHandler.cpp:29-38`:
/// entries in `std::map` key order.
#[test]
fn writer_emits_entries_in_key_order() {
    let map = PtmInformation::from([
        ("b".to_owned(), ("H2O".to_owned(), "S".to_owned())),
        ("a".to_owned(), ("C2".to_owned(), "KR".to_owned())),
    ]);
    let mut bytes = Vec::new();
    ptm_xml::write(&mut bytes, &map).unwrap();
    assert_eq!(
        String::from_utf8(bytes).unwrap(),
        "<PTMs>\n\
         \t<PTM>\n\t\t<name>a</name>\n\t\t<composition>C2</composition>\n\
         \t\t<possible_amino_acids>KR</possible_amino_acids>\n\t</PTM>\n\
         \t<PTM>\n\t\t<name>b</name>\n\t\t<composition>H2O</composition>\n\
         \t\t<possible_amino_acids>S</possible_amino_acids>\n\t</PTM>\n\
         </PTMs>\n"
    );
}

/// Source `writeTo` does not escape; this does, so markup characters survive.
#[test]
fn escaped_names_round_trip() {
    let map = PtmInformation::from([(
        "a<b>&\"c'".to_owned(),
        ("H2O".to_owned(), "S".to_owned()),
    )]);
    let mut bytes = Vec::new();
    ptm_xml::write(&mut bytes, &map).unwrap();
    assert!(
        String::from_utf8(bytes.clone())
            .unwrap()
            .contains("<name>a&lt;b&gt;&amp;&quot;c&apos;</name>")
    );
    assert_eq!(ptm_xml::read(bytes.as_slice()).unwrap(), map);
}

/// `PTMXMLFile::store` goes through `XMLFile::save_`, which compresses by the
/// `.gz` / `.bz2` suffix; `parse_` detects compression by magic bytes.
#[test]
fn compressed_store_and_load() {
    let map = ptm_xml::load(fixture()).unwrap();
    for name in ["c.xml.gz", "c.xml.bz2"] {
        let path = temporary(name);
        ptm_xml::store(&path, &map).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        assert!(bytes.starts_with(&[0x1f, 0x8b]) || bytes.starts_with(b"BZh"));
        assert_eq!(ptm_xml::load(&path).unwrap(), map);
        std::fs::remove_file(&path).unwrap();
    }
}

#[test]
fn failures_leave_no_partial_result() {
    assert!(matches!(
        ptm_xml::load(temporary("missing.xml")),
        Err(Error::Io(_))
    ));
    assert!(matches!(
        ptm_xml::read(b"<PTMs><PTM><name>A</name>".as_slice()),
        Err(Error::Parse { .. })
    ));
    let limits = ReadLimits {
        max_entries: 0,
        ..ReadLimits::default()
    };
    assert!(matches!(
        ptm_xml::load_with_limits(fixture(), &limits),
        Err(Error::InvalidValue(_))
    ));
}
