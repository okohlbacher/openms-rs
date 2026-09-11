// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::{
    format::{FileType, file_handler::type_by_content},
    metadata::DocumentIdentifier,
    system::file::TempDir,
};
use std::fs;

#[test]
fn source_identity_copy_assignment_swap_and_file_independent_equality() {
    let mut a = DocumentIdentifier::new();
    assert_eq!(a.identifier, "");
    assert_eq!(a.loaded_file_path, "");
    assert_eq!(a.loaded_file_type, FileType::Unknown);
    a.identifier = "this is a test".into(); // literal source class test
    a.loaded_file_path = "/source/file".into();
    a.loaded_file_type = FileType::MzMl;
    let mut b = a.clone();
    b.loaded_file_path = "different".into();
    b.loaded_file_type = FileType::Fasta;
    assert_eq!(a, b);
    b.identifier.push('!');
    assert_ne!(a, b);
    let mut empty = DocumentIdentifier::default();
    std::mem::swap(&mut a, &mut empty);
    assert_eq!(a.identifier, "");
    assert_eq!(a.loaded_file_path, "");
    assert_eq!(a.loaded_file_type, FileType::Unknown);
    assert_eq!(empty.identifier, "this is a test");
    assert_eq!(empty.loaded_file_path, "/source/file");
    assert_eq!(empty.loaded_file_type, FileType::MzMl);
    a = empty.clone();
    assert_eq!(a.loaded_file_path, empty.loaded_file_path);
    assert_eq!(a.loaded_file_type, empty.loaded_file_type);
}

#[test]
fn absolute_input_preserves_case_dot_segments_and_separators_without_io() {
    let mut value = DocumentIdentifier::default();
    for text in ["/", "/Case/./NotPresent/../Δ//file.mzML", "/missing/"] {
        value.set_loaded_file_path(text).unwrap();
        assert_eq!(value.loaded_file_path, text);
    }
    #[cfg(windows)]
    for text in [r"C:\Case\..\File", r"\\server\share\Case"] {
        value.set_loaded_file_path(text).unwrap();
        assert_eq!(value.loaded_file_path, text);
    }
}

#[test]
fn relative_and_empty_paths_use_current_directory_with_atomic_bounds() {
    let cwd = std::env::current_dir().unwrap();
    let mut value = DocumentIdentifier::default();
    for text in ["missing/../Δ/file", "."] {
        value.set_loaded_file_path(text).unwrap();
        let expected = cwd.join(text).to_str().unwrap().to_owned();
        #[cfg(windows)]
        let expected = expected.replace('\\', "/");
        assert_eq!(value.loaded_file_path, expected);
    }
    value.set_loaded_file_path("").unwrap();
    let expected = cwd.to_str().unwrap().to_owned();
    #[cfg(windows)]
    let expected = expected.replace('\\', "/");
    assert_eq!(value.loaded_file_path, expected);
    let old = value.loaded_file_path.clone();
    for text in ["bad\0path".into(), "x".repeat(1024 * 1024 + 1)] {
        assert!(value.set_loaded_file_path(&text).is_err());
        assert_eq!(value.loaded_file_path, old);
    }
    assert_eq!(value.loaded_file_type, FileType::Unknown);
}

#[test]
fn source_empty_file_and_content_selection_ignore_filename_and_stored_path() {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("deliberately-wrong.fasta");
    let mut value = DocumentIdentifier {
        identifier: "unchanged".into(),
        loaded_file_path: "unchanged".into(),
        loaded_file_type: FileType::Fasta,
    };
    fs::write(&path, b"").unwrap();
    value.set_loaded_file_type(&path).unwrap();
    assert_eq!(value.loaded_file_type, FileType::Unknown);
    fs::write(&path, b"<mzML>").unwrap();
    value.set_loaded_file_type(&path).unwrap();
    assert_eq!(value.loaded_file_type, FileType::MzMl);
    assert_eq!(value.identifier, "unchanged");
    assert_eq!(value.loaded_file_path, "unchanged");
    assert!(
        value
            .set_loaded_file_type(dir.path().join("missing"))
            .is_err()
    );
    assert_eq!(value.loaded_file_type, FileType::MzMl);
    fs::write(&path, b"PK\x03\x04").unwrap();
    assert!(value.set_loaded_file_type(&path).is_err());
    assert_eq!(value.loaded_file_type, FileType::MzMl);
}

#[test]
fn source_tabular_markers_and_detection_precedence() {
    for (text, kind) in [
        ("MTD\tmzTab-version\t1.0.0", FileType::MzTab),
        (
            "scan\ttime\tmz\taccurateMZ\tmass\tintensity\tcharge\tchargeStates\tkl\tbackground\tmedian\tpeaks\tscanFirst\tscanLast\tscanCount\ttotalIntensity\tsumSquaresDist\tdescription",
            FileType::Tsv,
        ),
        (
            "       m/z\t     rt(min)\t       snr\t      charge\t   intensity",
            FileType::Peplist,
        ),
        (
            "File\tFirst Scan\tLast Scan\tNum of Scans\tCharge\tMonoisotopic Mass\tBase Isotope Peak\tBest Intensity\tSummed Intensity\tFirst RTime\tLast RTime\tBest RTime\tBest Correlation\tModifications",
            FileType::Kroenik,
        ),
        (
            "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds",
            FileType::Psms,
        ),
    ] {
        assert_eq!(type_by_content(text.as_bytes()), kind, "{text}");
        assert_eq!(
            type_by_content(format!("<mzML>\n{text}").as_bytes()),
            FileType::MzMl
        );
    }
    assert_eq!(type_by_content(b"\nMTD\tmzTab-version"), FileType::MzTab);
    assert_eq!(
        type_by_content(b"\n\n\n\n\nMTD\tmzTab-version"),
        FileType::Unknown
    );
    assert_eq!(
        type_by_content(b"notPSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds"),
        FileType::Unknown
    );
}

#[test]
fn content_preview_bound_and_compression_feature_boundary() {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("content");
    let mut text = vec![b' '; 65_536];
    text.extend_from_slice(b"<mzML>");
    fs::write(&path, text).unwrap();
    let mut value = DocumentIdentifier::default();
    value.set_loaded_file_type(&path).unwrap();
    assert_eq!(value.loaded_file_type, FileType::Unknown);
    #[cfg(feature = "file-compression")]
    {
        use std::io::Write;
        let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gzip.write_all(b"<mzML>").unwrap();
        fs::write(&path, gzip.finish().unwrap()).unwrap();
        value.set_loaded_file_type(&path).unwrap();
        assert_eq!(value.loaded_file_type, FileType::MzMl);
        let mut bzip = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
        bzip.write_all(b">source\nPEPTIDE\n").unwrap();
        fs::write(&path, bzip.finish().unwrap()).unwrap();
        value.set_loaded_file_type(&path).unwrap();
        assert_eq!(value.loaded_file_type, FileType::Fasta);
    }
    #[cfg(not(feature = "file-compression"))]
    {
        fs::write(&path, b"\x1f\x8b000").unwrap();
        assert!(value.set_loaded_file_type(&path).is_err());
        assert_eq!(value.loaded_file_type, FileType::Unknown);
    }
}
