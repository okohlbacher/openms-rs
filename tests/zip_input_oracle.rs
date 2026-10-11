// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! ZIP input through the loaders (SHARED-02), against the executed Release
//! build.
//!
//! `tests/data/compressed_streams/b8a_oracle/records.tsv` is the output of
//! `../oracle/B8a-zip-input` (lane B8a): `FeatureXMLFile::load`/`loadSize` and
//! `MzMLFile::load`/`loadSize` of the Linux x86_64 Release install
//! `openms4-release-bc9cc12-c19e494-174b576`, run twice on ibminode06 with
//! identical records, over ZIP archives stored, deflated, bzip2-compressed,
//! LZMA-compressed (which that libzip cannot read), with a directory entry,
//! with two entries, empty, directory-only, with a bad CRC and with junk
//! after the signature, plus the two plain inputs.

#![cfg(all(feature = "featurexml", feature = "mzml"))]

use openms::Error;
use openms::format::{featurexml, mzml};
use std::collections::BTreeMap;
use std::io::ErrorKind;
use std::path::{Path, PathBuf};

fn dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/compressed_streams/b8a_oracle")
}

fn records() -> BTreeMap<(String, String), String> {
    let text = std::fs::read_to_string(dir().join("records.tsv")).unwrap();
    text.lines()
        .map(|line| {
            let mut fields = line.splitn(3, '\t');
            let function = fields.next().unwrap().to_owned();
            let file = fields.next().unwrap().to_owned();
            ((function, file), fields.next().unwrap_or("").to_owned())
        })
        .collect()
}

/// The exception class of an oracle error as thrown, and the class it wraps:
/// `MzMLFile::loadSize` (and `MzMLFile::load` of an archive that holds only
/// its end record) wraps the exception of `ZipIfstream::open` in a
/// `ParseError` ("due to that error of type X").
fn oracle_class(result: &str) -> Option<(&str, Option<&str>)> {
    let rest = result.strip_prefix("ERR:")?;
    let outer = rest.split(':').next()?;
    let inner = rest
        .split("due to that error of type ")
        .nth(1)
        .and_then(|inner| inner.split(" in:").next());
    Some((outer, inner))
}

fn port_class(error: &Error) -> &'static str {
    match error {
        Error::Parse { .. } => "Parse Error",
        Error::Io(e) if e.kind() == ErrorKind::NotFound => "FileNotFound",
        Error::Io(e) if e.kind() == ErrorKind::UnexpectedEof => "FileEmpty",
        Error::Io(e)
            if matches!(
                e.kind(),
                ErrorKind::InvalidData | ErrorKind::PermissionDenied
            ) =>
        {
            "FileNotReadable"
        }
        Error::InvalidValue(_) => "InvalidValue",
        Error::Unsupported(_) => "Unsupported",
        Error::MissingInformation(_) => "MissingInformation",
        _ => "other",
    }
}

/// Errors that come from the archive rather than from the XML inside it.
const ZIP_ERRORS: &[&str] = &["FileEmpty", "FileNotReadable", "FileNotFound"];

fn check(
    oracle: &BTreeMap<(String, String), String>,
    function: &str,
    file: &str,
    port: Result<String, Error>,
) -> bool {
    let expected = &oracle[&(function.to_owned(), file.to_owned())];
    match (&port, oracle_class(expected)) {
        (Ok(text), None) => {
            assert_eq!(text, expected, "{function} {file}");
            true
        }
        (Err(error), Some((outer, inner))) => {
            let class = inner.unwrap_or(outer);
            let archive_refusal =
                ZIP_ERRORS.contains(&class) || expected.contains("expected exactly 1:");
            if archive_refusal {
                if inner.is_some() && function == "MzMLFile::loadSize" {
                    // `MzMLFile::loadSize` reports the wrapper, and so does
                    // `mzml::load_size`.
                    assert_eq!(
                        port_class(error),
                        "Parse Error",
                        "{function} {file}: oracle {expected}, port {error:?}"
                    );
                } else if inner.is_some() {
                    // The one row where `MzMLFile::load` wraps (an archive that
                    // holds only its end record): the port reports the class
                    // of the exception itself, which `FeatureXMLFile::load`
                    // and the other `MzMLFile::load` rows do too.
                    assert!(
                        port_class(error) == class || port_class(error) == "Parse Error",
                        "{function} {file}: oracle {expected}, port {error:?}"
                    );
                } else {
                    assert_eq!(
                        port_class(error),
                        outer,
                        "{function} {file}: oracle {expected}, port {error:?}"
                    );
                }
                if expected.contains("expected exactly 1:") {
                    // The whole `what()`: the message, " in: ", and no expression.
                    let Error::Parse { message, .. } = error else {
                        unreachable!()
                    };
                    let whole = expected.strip_prefix("ERR:Parse Error:").unwrap();
                    assert_eq!(message, whole, "{function} {file}");
                }
            }
            // An XML-level refusal (the wrong format inside the archive) is
            // an error on both sides; its class is the loader's, not ZIP's.
            true
        }
        (Ok(text), Some(_)) => panic!("{function} {file}: oracle {expected}, port Ok({text})"),
        (Err(error), None) => panic!("{function} {file}: oracle {expected}, port {error:?}"),
    }
}

fn bits(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

#[test]
fn loaders_read_zip_archives_as_the_release_build_does() {
    let oracle = records();
    let mut files: Vec<String> = oracle.keys().map(|(_, file)| file.clone()).collect();
    files.sort();
    files.dedup();
    assert_eq!(files.len(), 13);
    let mut compared = 0;
    for file in &files {
        let path = dir().join(file);
        let features = featurexml::load(&path).map(|map| {
            let mut text = map.features.len().to_string();
            if let Some(first) = map.features.first() {
                text += &format!(
                    " rt={} mz={} uid={}",
                    bits(first.base.rt),
                    bits(first.base.mz),
                    first.base.unique_id
                );
            }
            text
        });
        compared += usize::from(check(&oracle, "FeatureXMLFile::load", file, features));
        let size = featurexml::load_size(&path, &featurexml::ReadOptions::default())
            .map(|n| n.to_string());
        compared += usize::from(check(&oracle, "FeatureXMLFile::loadSize", file, size));
        let experiment = mzml::load(&path).map(|exp| {
            let mut text = format!(
                "{} spectra {} chromatograms",
                exp.spectra.len(),
                exp.chromatograms.len()
            );
            if let Some(first) = exp.spectra.first() {
                text += &format!(
                    " rt0={} peaks0={} id0={}",
                    bits(first.rt),
                    first.peaks.len(),
                    first.native_id
                );
            }
            text
        });
        compared += usize::from(check(&oracle, "MzMLFile::load", file, experiment));
        let counts = mzml::load_size(&path)
            .map(|counts| format!("{} {}", counts.spectra, counts.chromatograms));
        compared += usize::from(check(&oracle, "MzMLFile::loadSize", file, counts));
    }
    assert_eq!(compared, 52);
}
