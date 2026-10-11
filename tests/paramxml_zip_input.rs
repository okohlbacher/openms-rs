// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `ParamXMLFile::load` of one-entry ZIP input, against the Linux x86_64
//! Release build.
//!
//! `ParamXMLFile::load` (`ParamXMLFile.cpp:357-361` at core `bc9cc12`) calls
//! `XMLFile::parse_`, whose `PK` branch (`XMLFile.cpp:141-165`) streams the
//! archive's single file entry; the port's `paramxml::load` opens its path
//! through `path_io::open`, which takes the same branch
//! (`src/format/zip_input.rs`).
//!
//! The archives in `tests/data/paramxml_zip_input/` were written by
//! `../oracle/integ-B-paramxml-zip/make_inputs.py` from three retained
//! documents of `tests/data/paramxml/` (the upstream
//! `ParamXMLFile_test_writeXMLToStream.xml`, `Param_advanced_no_required.ini`
//! and `Param_pre16_update.ini`): one entry deflated, stored and
//! bzip2-compressed, one entry below a directory entry, two entries, and an
//! empty archive. `oracle.tsv` is that oracle's `results/records.tsv` and
//! `release/<name>.xml` is `ParamXMLFile::writeXMLToStream` of what the
//! Release build loaded from `<name>`: ibminode06, the Release install
//! `openms4-release-bc9cc12-c19e494-174b576`, two runs byte-identical.
#![cfg(feature = "paramxml")]

use openms::format::paramxml::{self, Limits, WriteOptions};
use std::path::PathBuf;

const ORACLE: &str = include_str!("data/paramxml_zip_input/oracle.tsv");

fn data(path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(path)
}

fn fixture(name: &str) -> PathBuf {
    if name.ends_with(".zip") {
        data(&format!("paramxml_zip_input/{name}"))
    } else {
        data(&format!("paramxml/{name}"))
    }
}

/// The plain document an archive holds (`make_inputs.py`).
fn plain_of(name: &str) -> &'static str {
    match name {
        "advanced_deflate.zip" => "Param_advanced_no_required.ini",
        "pre16_deflate.zip" => "Param_pre16_update.ini",
        _ => "ParamXMLFile_test_writeXMLToStream.xml",
    }
}

#[test]
fn zip_input_matches_the_release_build() {
    let mut loaded = 0;
    let mut refused = 0;
    for line in ORACLE.lines() {
        let (name, expected) = line.split_once('\t').unwrap();
        let result = paramxml::load(fixture(name));
        if let Some(error) = expected.strip_prefix("ERR:") {
            let refusal = result.expect_err(name);
            if let Some(what) = error.strip_prefix("Parse Error:") {
                // The source's whole `what()`, " in: " and its empty
                // expression included, is the message of a line-0 parse error.
                match &refusal {
                    openms::Error::Parse { line: 0, message } => {
                        assert_eq!(message, what, "{name}");
                    }
                    other => panic!("{name}: C++ ParseError {what:?}, port {other:?}"),
                }
            } else {
                // `FileEmpty` names its message as the file: "the file
                // '<message>' is empty". The port has no file-empty error and
                // refuses with the message itself (src/format/zip_input.rs).
                let what = error
                    .strip_prefix("FileEmpty:the file '")
                    .and_then(|text| text.strip_suffix("' is empty"))
                    .unwrap();
                let message = refusal.to_string();
                assert!(
                    message.contains(what),
                    "{name}: C++ {what:?}, port {message:?}"
                );
            }
            refused += 1;
            continue;
        }
        let param = result.unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(expected, format!("OK {} entries", param.size()), "{name}");
        let mut written = Vec::new();
        paramxml::write_with_options(
            &mut written,
            &param,
            Limits::default(),
            WriteOptions::source(),
        )
        .unwrap();
        let release =
            std::fs::read(data(&format!("paramxml_zip_input/release/{name}.xml"))).unwrap();
        assert!(
            written == release,
            "{name}: port\n{}\nRelease\n{}",
            String::from_utf8_lossy(&written),
            String::from_utf8_lossy(&release)
        );
        if name.ends_with(".zip") {
            assert_eq!(
                param,
                paramxml::load(fixture(plain_of(name))).unwrap(),
                "{name}"
            );
        }
        loaded += 1;
    }
    assert_eq!((loaded, refused), (9, 2));
}
