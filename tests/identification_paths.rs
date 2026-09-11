// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "idxml")]

use openms::chemistry::ModificationsDB;
use openms::format::idxml::{self, IdXmlDocument, ReadOptions, WriteOptions};
use openms::format::{FileHandler, FileType};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

const SOURCE: &[u8] = include_bytes!("data/idxml_upstream_whole.idXML");
struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "openms-id-paths-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
    fn write(&self, name: &str, bytes: &[u8]) -> PathBuf {
        let path = self.path(name);
        std::fs::write(&path, bytes).unwrap();
        path
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
fn document() -> IdXmlDocument {
    idxml::read(SOURCE).unwrap()
}
fn assert_plain(path: &Path, document: &IdXmlDocument) {
    let bytes = std::fs::read(path).unwrap();
    assert!(bytes.starts_with(b"<?xml"));
    assert_eq!(idxml::read(bytes.as_slice()).unwrap(), *document);
}

#[test]
fn path_loading_preserves_source_document_and_uses_content_magic() {
    let directory = Directory::new();
    let expected = document();
    for name in ["source.idXML", "source.idXML.gz", "unknown.bin"] {
        let path = directory.write(name, SOURCE);
        assert_eq!(idxml::load(&path).unwrap(), expected);
        assert_eq!(
            FileHandler::load_identifications(&path, &[]).unwrap(),
            expected
        );
    }
    assert_eq!(
        idxml::load_with_registry(
            directory.path("source.idXML"),
            &ReadOptions::default(),
            ModificationsDB::global()
        )
        .unwrap(),
        expected
    );
    assert_eq!(
        FileHandler::read_identifications(SOURCE, FileType::IdXml).unwrap(),
        expected
    );
}

#[test]
fn load_failure_keeps_destination_and_checks_decoded_limits() {
    let directory = Directory::new();
    let path = directory.write("bad.idXML", b"<IdXML><broken></IdXML>");
    let mut destination = document();
    destination.document_id = "existing".into();
    let before = destination.clone();
    assert!(idxml::load_into(&path, &mut destination).is_err());
    assert_eq!(destination, before);
    assert!(idxml::load_into(directory.path("missing.idXML"), &mut destination).is_err());
    assert_eq!(destination, before);
    std::fs::write(&path, SOURCE).unwrap();
    assert!(
        idxml::load_with_options(
            &path,
            &ReadOptions {
                max_xml_bytes: 64,
                ..Default::default()
            }
        )
        .is_err()
    );
    idxml::load_into(&path, &mut destination).unwrap();
    assert_eq!(destination, document());
}

#[test]
fn source_plain_output_rule_and_atomic_publication() {
    let directory = Directory::new();
    let expected = document();
    for name in [
        "out.idXML",
        "out.idXML.gz",
        "out.idXML.bz2",
        "out.idXML.zip",
    ] {
        let path = directory.path(name);
        idxml::store(&path, &expected).unwrap();
        assert_plain(&path, &expected);
        let before = std::fs::read(&path).unwrap();
        assert!(idxml::store(&path, &IdXmlDocument::default()).is_err());
        assert_eq!(std::fs::read(&path).unwrap(), before);
        assert!(
            idxml::store_with_options(
                &path,
                &expected,
                &WriteOptions {
                    max_xml_bytes: 1,
                    ..Default::default()
                }
            )
            .is_err()
        );
        assert_eq!(std::fs::read(&path).unwrap(), before);
    }
    for name in ["wrong.mzML", "wrong.mzML.gz"] {
        let path = directory.write(name, b"existing destination");
        assert!(idxml::store(&path, &expected).is_err());
        assert_eq!(std::fs::read(&path).unwrap(), b"existing destination");
    }
    let unknown = directory.path("unknown.extension");
    idxml::store(&unknown, &expected).unwrap();
    assert_plain(&unknown, &expected);
    let path = directory.path("registry.idXML");
    idxml::store_with_registry(
        &path,
        &expected,
        &WriteOptions::default(),
        ModificationsDB::global(),
    )
    .unwrap();
    assert_plain(&path, &expected);
    assert_eq!(std::fs::read_dir(&directory.0).unwrap().count(), 8);
}

#[test]
fn source_output_allowlist_fallback_and_input_recognition() {
    let directory = Directory::new();
    let expected = document();
    assert!(FileHandler::can_read_identifications(FileType::IdXml));
    assert!(FileHandler::can_write_identifications(FileType::IdXml));
    for kind in [FileType::MzMl, FileType::MzIdentMl, FileType::Unknown] {
        assert!(!FileHandler::can_read_identifications(kind));
        assert!(FileHandler::read_identifications(SOURCE, kind).is_err());
        let mut output = Vec::new();
        assert!(FileHandler::write_identifications(&mut output, &expected, kind).is_err());
        assert!(output.is_empty());
    }
    let path = directory.path("output.unknown");
    assert!(FileHandler::store_identifications(&path, &expected, &[]).is_err());
    assert!(!path.exists());
    assert!(
        FileHandler::store_identifications(
            &path,
            &expected,
            &[FileType::IdXml, FileType::MzIdentMl]
        )
        .is_err()
    );
    FileHandler::store_identifications(&path, &expected, &[FileType::IdXml]).unwrap();
    assert_eq!(
        FileHandler::load_identifications(&path, &[FileType::IdXml]).unwrap(),
        expected
    );
    assert!(FileHandler::load_identifications(&path, &[FileType::MzIdentMl]).is_err());
    let misleading = directory.write("source.mzML", SOURCE);
    assert!(FileHandler::load_identifications(&misleading, &[]).is_err());
    let before = std::fs::read(&misleading).unwrap();
    assert!(
        FileHandler::store_identifications(&misleading, &expected, &[FileType::IdXml]).is_err()
    );
    assert_eq!(std::fs::read(&misleading).unwrap(), before);
    let compressed_name = directory.path("stored.idXML.gz");
    FileHandler::store_identifications(&compressed_name, &expected, &[]).unwrap();
    assert_plain(&compressed_name, &expected);
}

#[test]
#[cfg(feature = "file-compression")]
fn gzip_and_bzip_input_are_detected_independently_of_extension() {
    use std::io::Write;
    let directory = Directory::new();
    let expected = document();
    let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gzip.write_all(SOURCE).unwrap();
    let mut bzip = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
    bzip.write_all(SOURCE).unwrap();
    for (name, encoded) in [
        ("gzip.data", gzip.finish().unwrap()),
        ("bzip.idXML", bzip.finish().unwrap()),
    ] {
        let path = directory.write(name, &encoded);
        assert_eq!(idxml::load(&path).unwrap(), expected);
        assert_eq!(
            FileHandler::load_identifications(&path, &[]).unwrap(),
            expected
        );
        assert!(
            idxml::load_with_options(
                &path,
                &ReadOptions {
                    max_xml_bytes: 64,
                    ..Default::default()
                }
            )
            .is_err()
        );
        let mut truncated = encoded.clone();
        truncated.truncate(encoded.len() / 2);
        std::fs::write(&path, truncated).unwrap();
        let mut destination = expected.clone();
        assert!(idxml::load_into(&path, &mut destination).is_err());
        assert_eq!(destination, expected);
    }
}

#[test]
#[cfg(not(feature = "file-compression"))]
fn compression_feature_errors_preserve_plain_loading() {
    let directory = Directory::new();
    for (name, bytes) in [
        ("gzip.idXML", b"\x1f\x8b".as_slice()),
        ("bzip.idXML", b"BZh".as_slice()),
    ] {
        let path = directory.write(name, bytes);
        assert!(matches!(
            idxml::load(path),
            Err(openms::Error::Unsupported(_))
        ));
    }
    assert_eq!(
        idxml::load(directory.write("plain.idXML.gz", SOURCE)).unwrap(),
        document()
    );
}
