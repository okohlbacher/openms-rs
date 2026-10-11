// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The retained inputs of the FeatureFinderIdentification and
//! FeatureFinderMetaboIdent tests, as shared **test support**.
//!
//! The upstream input `FeatureFinderIdentification_1_input.idXML` cannot be read
//! by the native idXML reader as it is, for two reasons that are the reader's,
//! not the algorithm's (both are recorded in the lane report):
//!
//! 1. Two hits carry the legacy attribute `prot="PH_14"` (lines 81 and 135).
//!    The C++ `IdXMLFile` handler never reads it and ignores it; the native
//!    reader refuses an unknown attribute. It is removed here, which is all the
//!    C++ reader does with it.
//! 2. The reader spends a flat work allowance of 50 million units, and every
//!    annotated sequence costs a scan of the whole registry. With the global
//!    registry of about 7,000 modifications the 21 `(Carbamidomethyl)` hits of
//!    this 23 KB document exhaust it ("peptide parent parsing resource limit
//!    exceeded"). The registry passed here holds the one named modification the
//!    document uses; the two mass annotations stay mass tags, as in the source
//!    (which logs "unknown modification ... adding it to the database").
#![allow(dead_code)]

use openms::chemistry::ModificationsDB;
use openms::format::featurexml;
use openms::format::idxml::{self, IdXmlDocument};
use openms::kernel::FeatureMap;
use openms::system::file::TempDir;
use std::fs;
use std::io::BufReader;
use std::path::{Path, PathBuf};

/// The directory of the retained fixtures.
pub fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/feature_finder_identification")
        .join(name)
}

/// A registry holding the one named modification the fixtures use (see the
/// module documentation for why the global registry is not passed).
pub fn registry() -> ModificationsDB {
    let records: Vec<_> = ModificationsDB::global()
        .entries()
        .iter()
        .filter(|m| m.name() == "Carbamidomethyl")
        .map(|m| (**m).clone())
        .collect();
    assert!(!records.is_empty());
    ModificationsDB::from_records(records).unwrap()
}

/// A featureXML file, read against [`registry`] (the featureXML reader spends
/// the same flat allowance per annotated sequence as the idXML one).
pub fn load_feature_map(path: &Path) -> FeatureMap {
    let file = fs::File::open(path).unwrap();
    featurexml::read_with_registry(
        BufReader::new(file),
        &featurexml::ReadOptions::default(),
        &registry(),
    )
    .unwrap()
}

/// A feature map as the tools leave it in their output file: written, then read.
pub fn through_file(map: &FeatureMap, work: &TempDir, name: &str) -> FeatureMap {
    let path = work.path().join(name);
    let mut bytes = Vec::new();
    featurexml::write_with_registry(
        &mut bytes,
        map,
        &featurexml::WriteOptions::default(),
        &registry(),
    )
    .unwrap();
    fs::write(&path, bytes).unwrap();
    load_feature_map(&path)
}

/// The identifications of `FeatureFinderIdentification_1_input.idXML`.
pub fn load_identifications(work: &TempDir) -> IdXmlDocument {
    const ATTRIBUTE: &str = " prot=\"";
    let original = fs::read_to_string(data("FeatureFinderIdentification_1_input.idXML")).unwrap();
    let mut cleaned = String::with_capacity(original.len());
    let mut rest = original.as_str();
    let mut removed = 0;
    while let Some(at) = rest.find(ATTRIBUTE) {
        cleaned.push_str(&rest[..at]);
        let after = &rest[at + ATTRIBUTE.len()..];
        let end = after.find('"').unwrap();
        rest = &after[end + 1..];
        removed += 1;
    }
    cleaned.push_str(rest);
    assert_eq!(removed, 2, "the two legacy `prot` attributes");
    let path = work
        .path()
        .join("FeatureFinderIdentification_1_input.idXML");
    fs::write(&path, cleaned).unwrap();

    let file = fs::File::open(&path).unwrap();
    idxml::read_with_registry(
        BufReader::new(file),
        &idxml::ReadOptions::default(),
        &registry(),
    )
    .unwrap()
}
