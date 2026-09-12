// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Declared mzML list `count` attributes are advisory on reading.
//!
//! `MzMLHandler.cpp` reads the attribute in exactly four places — a progress
//! range plus `reserveSpaceSpectra` for `spectrumList` (:965-979), the
//! chromatogram equivalent (:996-1013), `bin_data_.reserve(...)` for
//! `binaryDataArrayList` (:1015-1017) and a warning when `selectedIonList`
//! exceeds one ion (:1371-1375) — and never compares it with the actual number
//! of children. `precursorList`, `productList` and `scanWindowList` have no
//! open-tag handler at all. So a wrong count loads upstream, and this port must
//! load it too; only the attribute's presence and numeric form stay checked,
//! and declared values are still spent as resource ceilings.
#![cfg(feature = "mzml")]

use openms::format::mzml::{self, ReadOptions};
use std::io::Cursor;

/// The upstream `MzMLFile_test.cpp` input, unmodified.
const SOURCE: &str = include_str!("data/mzml_validator/MzMLFile_1.mzML");

fn read(xml: &str) -> openms::Result<openms::MSExperiment> {
    mzml::read(Cursor::new(xml.as_bytes()))
}
/// One spectrum carrying `body`, inside a one-record `spectrumList`.
fn doc(body: &str) -> String {
    format!(
        "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\"><run><spectrumList count=\"1\"><spectrum id=\"scan=1\" defaultArrayLength=\"0\">{body}</spectrum></spectrumList></run></mzML>"
    )
}
const SCANS: &str = "<scanList count=\"1\"><scan><scanWindowList count=\"1\"><scanWindow><cvParam accession=\"MS:1000501\" value=\"1\"/><cvParam accession=\"MS:1000500\" value=\"2\"/></scanWindow></scanWindowList></scan></scanList>";
const PRECURSORS: &str = "<precursorList count=\"1\"><precursor><selectedIonList count=\"1\"><selectedIon><cvParam accession=\"MS:1000744\" value=\"500\"/></selectedIon></selectedIonList><activation/></precursor></precursorList>";
const PRODUCTS: &str = "<productList count=\"1\"><product><isolationWindow><cvParam accession=\"MS:1000827\" value=\"18.88\"/></isolationWindow></product></productList>";

/// Declared count and actual child count for every `list` instance in `xml`.
fn declared_and_actual(xml: &str, list: &str, child: &str) -> Vec<(usize, usize)> {
    let open = format!("<{list} ");
    let close = format!("</{list}>");
    let mut found = Vec::new();
    let mut rest = xml;
    while let Some(start) = rest.find(&open) {
        let head = &rest[start + open.len()..];
        let attributes = &head[..head.find('>').expect("unterminated list element")];
        let value = attributes
            .split_once("count=\"")
            .expect("list without a count attribute")
            .1;
        let declared = value[..value.find('"').expect("unterminated count")]
            .parse()
            .expect("non-numeric emitted count");
        let body_start = start + open.len() + attributes.len() + 1;
        let body_end = body_start
            + rest[body_start..]
                .find(&close)
                .expect("unterminated list body");
        let body = &rest[body_start..body_end];
        let actual = [
            format!("<{child}>"),
            format!("<{child} "),
            format!("<{child}/>"),
        ]
        .iter()
        .map(|pattern| body.matches(pattern.as_str()).count())
        .sum();
        found.push((declared, actual));
        rest = &rest[body_end + close.len()..];
    }
    found
}

#[test]
fn the_upstream_class_fixture_with_wrong_declared_counts_loads() {
    // The mismatch this test exists for: `MzMLFile_1.mzML` declares two binary
    // data arrays on its second spectrum and carries four, and declares one
    // product on its third and carries two. `MzMLFile_test.cpp` loads the file
    // regardless, so rejecting it left this port unable to read its own
    // reference data.
    assert_eq!(
        declared_and_actual(SOURCE, "binaryDataArrayList", "binaryDataArray")[1],
        (2, 4)
    );
    assert_eq!(
        declared_and_actual(SOURCE, "productList", "product"),
        vec![(1, 2)]
    );
    let e = read(SOURCE).unwrap();
    assert_eq!((e.spectra.len(), e.chromatograms.len()), (4, 2));
    // All four of the under-declared arrays are present: two primary plus the
    // signal-to-noise and non-standard auxiliary arrays.
    assert_eq!(e.spectra[1].peaks.len(), 10);
    assert_eq!(
        e.spectra[1]
            .float_data_arrays
            .iter()
            .map(|a| a.name.as_str())
            .collect::<Vec<_>>(),
        ["signal to noise array", "user-defined name"]
    );
    // Both of the under-declared products are present.
    assert_eq!(
        e.spectra[2]
            .products
            .iter()
            .map(|p| p.mz)
            .collect::<Vec<_>>(),
        [18.88, 19.99]
    );
}

#[test]
fn writing_emits_the_true_count_for_every_list() {
    let mut bytes = Vec::new();
    mzml::write(&mut bytes, &read(SOURCE).unwrap()).unwrap();
    let xml = String::from_utf8(bytes).unwrap();
    for (list, child) in [
        ("spectrumList", "spectrum"),
        ("chromatogramList", "chromatogram"),
        ("binaryDataArrayList", "binaryDataArray"),
        ("precursorList", "precursor"),
        ("productList", "product"),
        ("selectedIonList", "selectedIon"),
        ("scanList", "scan"),
        ("scanWindowList", "scanWindow"),
    ] {
        let pairs = declared_and_actual(&xml, list, child);
        assert!(!pairs.is_empty(), "{list} never emitted");
        for (declared, actual) in pairs {
            assert_eq!(declared, actual, "{list} emitted a wrong count");
        }
    }
    // The output therefore corrects the input: four arrays and two products.
    assert_eq!(
        declared_and_actual(&xml, "binaryDataArrayList", "binaryDataArray")[1],
        (4, 4)
    );
    assert_eq!(
        declared_and_actual(&xml, "productList", "product"),
        vec![(2, 2)]
    );
    // And the corrected document reads back to the same experiment.
    assert_eq!(read(&xml).unwrap(), read(SOURCE).unwrap());
}

#[test]
fn a_wrong_count_is_accepted_at_every_reading_site() {
    let body = format!("{SCANS}{PRECURSORS}{PRODUCTS}");
    let source = doc(&body);
    let expected = read(&source).unwrap();
    for (from, to) in [
        ("spectrumList count=\"1\"", "spectrumList count=\"7\""),
        ("scanList count=\"1\"", "scanList count=\"0\""),
        ("scanWindowList count=\"1\"", "scanWindowList count=\"3\""),
        ("precursorList count=\"1\"", "precursorList count=\"2\""),
        ("selectedIonList count=\"1\"", "selectedIonList count=\"5\""),
        ("productList count=\"1\"", "productList count=\"4\""),
    ] {
        let changed = source.replacen(from, to, 1);
        assert_ne!(changed, source);
        assert_eq!(read(&changed).unwrap(), expected, "{to}");
    }
    // Records with binary arrays, and the chromatogram list, behave the same.
    let arrays = SOURCE.replacen(
        "<binaryDataArrayList count=\"2\">",
        "<binaryDataArrayList count=\"11\">",
        1,
    );
    let both = arrays.replacen(
        "<chromatogramList count=\"2\"",
        "<chromatogramList count=\"6\"",
        1,
    );
    assert_eq!(read(&both).unwrap(), read(SOURCE).unwrap());
}

#[test]
fn the_count_attribute_stays_required_and_numeric() {
    let source = doc(&format!("{SCANS}{PRECURSORS}{PRODUCTS}"));
    for (from, to) in [
        ("spectrumList count=\"1\"", "spectrumList"),
        ("spectrumList count=\"1\"", "spectrumList count=\"\""),
        ("spectrumList count=\"1\"", "spectrumList count=\"1.0\""),
        ("spectrumList count=\"1\"", "spectrumList count=\"-1\""),
        ("scanList count=\"1\"", "scanList"),
        ("scanList count=\"1\"", "scanList count=\"many\""),
        ("scanWindowList count=\"1\"", "scanWindowList"),
        ("scanWindowList count=\"1\"", "scanWindowList count=\" 1\""),
        ("precursorList count=\"1\"", "precursorList"),
        ("precursorList count=\"1\"", "precursorList count=\"1 \""),
        ("selectedIonList count=\"1\"", "selectedIonList"),
        (
            "selectedIonList count=\"1\"",
            "selectedIonList count=\"0x1\"",
        ),
        ("productList count=\"1\"", "productList"),
        ("productList count=\"1\"", "productList count=\"one\""),
    ] {
        let changed = source.replacen(from, to, 1);
        assert_ne!(changed, source);
        assert!(read(&changed).is_err(), "{to} accepted");
    }
    for to in [
        "<binaryDataArrayList>",
        "<binaryDataArrayList count=\"two\">",
    ] {
        let changed = SOURCE.replacen("<binaryDataArrayList count=\"2\">", to, 1);
        assert_ne!(changed, SOURCE);
        assert!(read(&changed).is_err(), "{to} accepted");
    }
    for to in ["<chromatogramList", "<chromatogramList count=\"nope\""] {
        let changed = SOURCE.replacen("<chromatogramList count=\"2\"", to, 1);
        assert_ne!(changed, SOURCE);
        assert!(read(&changed).is_err(), "{to} accepted");
    }
}

#[test]
fn declared_counts_remain_resource_ceilings_before_allocation() {
    let options = ReadOptions {
        max_total_params: 1000,
        ..Default::default()
    };
    let read_bounded =
        |xml: &str| mzml::read_with_options(Cursor::new(xml.as_bytes()), &options).map(|_| ());
    for body in [
        format!("{SCANS}{PRECURSORS}{PRODUCTS}"),
        format!("{SCANS}{PRODUCTS}"),
    ] {
        assert!(read_bounded(&doc(&body)).is_ok(), "{body}");
    }
    // A hostile declared count is refused before the parameter arena grows,
    // exactly as before this fix: only the comparison with the children went.
    for ((from, to), message) in [
        (
            ("productList count=\"1\"", "productList count=\"100000000\""),
            "product count exceeds parameter limit",
        ),
        (
            (
                "scanWindowList count=\"1\"",
                "scanWindowList count=\"100000000\"",
            ),
            "scan window count exceeds parameter limit",
        ),
        (
            ("scanList count=\"1\"", "scanList count=\"100000000\""),
            "scan count exceeds parameter limit",
        ),
    ] {
        let body = format!("{SCANS}{PRECURSORS}{PRODUCTS}").replacen(from, to, 1);
        let error = read_bounded(&doc(&body)).unwrap_err().to_string();
        assert!(error.contains(message), "{to} gave {error}");
    }
    // The parameter-group list keeps its configured ceiling too.
    let groups = "<referenceableParamGroupList count=\"2\"><referenceableParamGroup id=\"a\"><userParam name=\"x\" value=\"1\"/></referenceableParamGroup><referenceableParamGroup id=\"b\"><userParam name=\"y\" value=\"2\"/></referenceableParamGroup></referenceableParamGroupList>";
    let with_groups = doc("").replacen("<run>", &format!("{groups}<run>"), 1);
    assert!(read(&with_groups).is_ok());
    let error = mzml::read_with_options(
        Cursor::new(with_groups.as_bytes()),
        &ReadOptions {
            max_param_groups: 1,
            ..Default::default()
        },
    )
    .unwrap_err()
    .to_string();
    assert!(
        error.contains("parameter group count is zero or exceeds configured limit"),
        "{error}"
    );
    // Record and binary-array counts drive no reserve, so an absurd declared
    // value is accepted without allocating for it.
    let absurd = doc(&format!("{SCANS}{PRODUCTS}")).replacen(
        "spectrumList count=\"1\"",
        "spectrumList count=\"1000000000\"",
        1,
    );
    assert_eq!(read(&absurd).unwrap().spectra.len(), 1);
    let absurd = SOURCE.replacen(
        "<binaryDataArrayList count=\"2\">",
        "<binaryDataArrayList count=\"1000000000\">",
        1,
    );
    assert_eq!(read(&absurd).unwrap(), read(SOURCE).unwrap());
}

#[test]
fn the_parameter_group_list_count_stays_strict_as_a_documented_divergence() {
    // Upstream has no `referenceableParamGroupList` handler either, so this
    // check is not source behaviour. It is kept so that `read` and `read_size`
    // agree on the same document: `src/format/mzml_counts.rs` rejects the same
    // mismatch and that reader is outside this fix. The fixture declares two
    // groups and carries two, so both readers accept it unchanged.
    assert!(read(SOURCE).is_ok());
    assert!(mzml::read_size(Cursor::new(SOURCE.as_bytes())).is_ok());
    let broken = SOURCE.replacen(
        "<referenceableParamGroupList count=\"2\">",
        "<referenceableParamGroupList count=\"3\">",
        1,
    );
    assert_ne!(broken, SOURCE);
    assert!(read(&broken).is_err());
    assert!(mzml::read_size(Cursor::new(broken.as_bytes())).is_err());
}

#[test]
fn advisory_counts_do_not_relax_list_structure() {
    let source = doc(&format!("{SCANS}{PRECURSORS}{PRODUCTS}"));
    for (from, to) in [
        (
            "</scanList>",
            "</scanList><scanList count=\"1\"><scan/></scanList>",
        ),
        (
            "</precursorList>",
            "</precursorList><precursorList count=\"0\"/>",
        ),
        ("</productList>", "</productList><productList count=\"0\"/>"),
        ("<scan>", "<scan><precursorList count=\"0\"/>"),
        ("<product>", "<product><scanWindowList count=\"0\"/>"),
    ] {
        let changed = source.replacen(from, to, 1);
        assert_ne!(changed, source);
        assert!(read(&changed).is_err(), "{to} accepted");
    }
    let duplicate = SOURCE.replacen(
        "</binaryDataArrayList>",
        "</binaryDataArrayList><binaryDataArrayList count=\"0\"/>",
        1,
    );
    assert!(read(&duplicate).is_err());
}
