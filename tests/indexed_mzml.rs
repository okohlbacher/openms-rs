// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

#![cfg(feature = "mzml")]

use openms::format::indexed_mzml::{IndexOffsets, IndexedMzMLDecoder, has_index};
use openms::system::file::TempDir;
use std::path::PathBuf;

fn file(bytes: &[u8]) -> (TempDir, PathBuf) {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("input.mzML");
    std::fs::write(&path, bytes).unwrap();
    (dir, path)
}

fn index(sections: &str) -> String {
    format!(
        "<indexList count=\"2\">{sections}</indexList><indexListOffset>0</indexListOffset><fileChecksum>0</fileChecksum></indexedmzML>"
    )
}

#[test]
fn unchanged_source_file_has_literal_offsets() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/indexed_mzml/IndexedmzMLFile_1.mzML");
    let decoder = IndexedMzMLDecoder::default();
    assert!(has_index(&path).unwrap());
    assert_eq!(decoder.find_index_list_offset(&path).unwrap(), Some(667742));
    let offsets = decoder.parse_offsets(&path, 667742).unwrap();
    assert_eq!(
        offsets.spectra,
        vec![
            ("controllerType=0 controllerNumber=1 scan=1".into(), 24146),
            ("controllerType=0 controllerNumber=1 scan=2".into(), 345745),
        ]
    );
    assert_eq!(offsets.chromatograms, vec![("TIC".into(), 665563)]);
}

#[test]
fn probe_matches_source_prefix_and_first_match_rules() {
    let decoder = IndexedMzMLDecoder::default();
    for (bytes, expected) in [
        (b"".as_slice(), None),
        (b"<indexedmzML/>".as_slice(), None),
        (b"<indexListOffset>0</indexListOffset>".as_slice(), Some(0)),
        (b"<x:indexListOffset \t>\n 123tail".as_slice(), Some(123)),
        (b"<anyindexListOffset>7".as_slice(), Some(7)),
        (
            b"</indexListOffset><indexListOffset>19".as_slice(),
            Some(19),
        ),
        (
            b"<indexListOffset/><indexListOffset>21".as_slice(),
            Some(21),
        ),
        (b"<indexListOffset> <indexListOffset>13".as_slice(), None),
        (b"\0<indexListOffset>25".as_slice(), None),
        (b"<indexListOffset>-1".as_slice(), None),
    ] {
        let (_dir, path) = file(bytes);
        assert_eq!(
            decoder.find_index_list_offset(&path).unwrap(),
            expected,
            "{bytes:?}"
        );
    }
}

#[test]
fn probe_window_and_overflow_are_bounded() {
    let mut bytes = b"<indexListOffset>123".to_vec();
    bytes.extend(vec![b' '; 1023]);
    let (_dir, path) = file(&bytes);
    let mut decoder = IndexedMzMLDecoder::default();
    assert_eq!(decoder.find_index_list_offset(&path).unwrap(), None);
    assert_eq!(
        decoder
            .find_index_list_offset_with_buffer_size(&path, bytes.len())
            .unwrap(),
        Some(123)
    );
    assert_eq!(
        decoder
            .find_index_list_offset_with_buffer_size(&path, 0)
            .unwrap(),
        None
    );
    decoder.limits.max_footer_bytes = 10;
    assert!(decoder.find_index_list_offset(&path).is_err());
    for value in ["9223372036854775808", "18446744073709551616"] {
        let (_dir, path) = file(format!("<indexListOffset>{value}").as_bytes());
        assert!(
            IndexedMzMLDecoder::default()
                .find_index_list_offset(path)
                .is_err()
        );
    }
}

#[test]
fn compact_xml_keeps_first_offset_and_decodes_ids_and_numbers() {
    let xml = index(
        "<index name=\"spectrum\"><offset idRef=\"α&amp;β&#32;x\">&#49;<![CDATA[23]]></offset><offset idRef=\"same\">+4</offset><offset idRef=\"same\">5</offset></index><index name=\"chromatogram\"><offset>0</offset></index>",
    );
    let (_dir, path) = file(xml.as_bytes());
    let offsets = IndexedMzMLDecoder::default()
        .parse_offsets(path, 0)
        .unwrap();
    assert_eq!(
        offsets.spectra,
        vec![
            ("α&β x".into(), 123),
            ("same".into(), 4),
            ("same".into(), 5)
        ]
    );
    assert_eq!(offsets.chromatograms, vec![(String::new(), 0)]);
}

#[test]
fn repeated_sections_replace_and_absent_sections_are_empty() {
    let xml = index(
        "<index name=\"spectrum\"><offset idRef=\"old\">8</offset></index><index name=\"spectrum\"><offset idRef=\"new\">9</offset></index>",
    );
    let (_dir, path) = file(xml.as_bytes());
    let mut decoder = IndexedMzMLDecoder::default();
    assert_eq!(
        decoder.parse_offsets(&path, 0).unwrap(),
        IndexOffsets {
            spectra: vec![("new".into(), 9)],
            chromatograms: vec![],
        }
    );
    decoder.limits.max_offsets = 1;
    assert!(decoder.parse_offsets(path, 0).is_err());
    let (_dir, path) = file(index("").as_bytes());
    decoder.limits.max_offsets = 0;
    assert_eq!(
        decoder.parse_offsets(path, 0).unwrap(),
        IndexOffsets::default()
    );
}

#[test]
fn invalid_offsets_and_index_xml_return_no_partial_output() {
    for section in [
        "<index name=\"other\"/>",
        "<index><offset>1</offset></index>",
        "<offset>1</offset>",
        "<index name=\"spectrum\"><offset> 1</offset></index>",
        "<index name=\"spectrum\"><offset>-1</offset></index>",
        "<index name=\"spectrum\"><offset>1.0</offset></index>",
        "<index name=\"spectrum\"><offset>9223372036854775808</offset></index>",
        "<index name=\"spectrum\"><offset>&unknown;</offset></index>",
        "<index name=\"spectrum\"><offset>&#0;</offset></index>",
        "<index name=\"spectrum\"><offset idRef=\"x\" idRef=\"y\">1</offset></index>",
        "<index name=\"spectrum\"><offset><nested>1</nested></offset></index>",
        "<index name=\"spectrum\"><offset>1</offset></wrong>",
    ] {
        let (_dir, path) = file(index(section).as_bytes());
        assert!(
            IndexedMzMLDecoder::default()
                .parse_offsets(path, 0)
                .is_err(),
            "{section}"
        );
    }
    for xml in [
        "</indexedmzML>",
        "<indexList/><indexList/></indexedmzML>",
        "<indexList/>",
        "<!DOCTYPE x><indexList/></indexedmzML>",
        "<indexList/></indexedmzML><extra/>",
        "<indexList/></indexedmzML>&#32;",
        "<indexList/></indexedmzML><![CDATA[ ]]>",
        "<indexList/><!-----></indexedmzML>",
    ] {
        let (_dir, path) = file(xml.as_bytes());
        assert!(
            IndexedMzMLDecoder::default()
                .parse_offsets(path, 0)
                .is_err(),
            "{xml}"
        );
    }
}

#[test]
fn byte_depth_identifier_and_record_limits_precede_results() {
    let xml = index("<index name=\"spectrum\"><offset idRef=\"alpha\">1</offset></index>");
    let (_dir, path) = file(xml.as_bytes());
    let original = IndexedMzMLDecoder::default();
    let mut cases = [original; 4];
    cases[0].limits.max_index_bytes = xml.len() - 1;
    cases[1].limits.max_depth = 3;
    cases[2].limits.max_id_bytes = 4;
    cases[3].limits.max_offsets = 0;
    for decoder in cases {
        assert!(decoder.parse_offsets(&path, 0).is_err());
    }
    let mut exact = original;
    exact.limits.max_index_bytes = xml.len();
    exact.limits.max_depth = 4;
    exact.limits.max_id_bytes = 5;
    exact.limits.max_offsets = 1;
    assert_eq!(exact.parse_offsets(path, 0).unwrap().spectra.len(), 1);
}

#[test]
fn io_and_address_errors_are_distinct_from_no_index() {
    let (_dir, path) = file(index("").as_bytes());
    let decoder = IndexedMzMLDecoder::default();
    let length = std::fs::metadata(&path).unwrap().len();
    assert!(decoder.parse_offsets(&path, length).is_err());
    assert!(decoder.parse_offsets(&path, length + 1).is_err());
    std::fs::remove_file(&path).unwrap();
    assert!(matches!(
        decoder.find_index_list_offset(&path),
        Err(openms::Error::Io(_))
    ));
    assert!(matches!(
        decoder.parse_offsets(path, 0),
        Err(openms::Error::Io(_))
    ));
}

#[test]
fn xml_lexical_boundaries_and_attribute_work_are_checked() {
    for xml in [
        "<indexList/><1bad/></indexedmzML>",
        "<indexList/><other 1bad=\"x\"/></indexedmzML>",
        "<indexList count=\"0\"other=\"x\"/></indexedmzML>",
        "<indexList/></indexedmzML>\u{a0}",
        "<indexList/><?XML bad?></indexedmzML>",
    ] {
        let (_dir, path) = file(xml.as_bytes());
        assert!(
            IndexedMzMLDecoder::default()
                .parse_offsets(path, 0)
                .is_err(),
            "{xml}"
        );
    }
    let xml = index("<index name=\"spectrum\"><offset idRef=\"a\r\nb&#10;c\">1</offset></index>");
    let (_dir, path) = file(xml.as_bytes());
    assert_eq!(
        IndexedMzMLDecoder::default()
            .parse_offsets(path, 0)
            .unwrap()
            .spectra[0]
            .0,
        "a b\nc"
    );
    let attributes = (0..1000)
        .map(|n| format!("a{n}=\"x\" "))
        .collect::<String>();
    let (_dir, path) = file(format!("<indexList {attributes}/></indexedmzML>").as_bytes());
    assert!(
        IndexedMzMLDecoder::default()
            .parse_offsets(path, 0)
            .is_err()
    );
}
