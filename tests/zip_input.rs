// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! One-entry ZIP input through the path readers, against the Linux x86_64
//! Release build.
//!
//! The fixtures in `tests/data/zip_input/` were written with Python's `zipfile`
//! from `tests/data/featurexml_source_1.featureXML` (the upstream
//! `FeatureXMLFile_1.featureXML`) and `tests/data/mzml_upstream_serum.mzML`:
//! one entry stored, deflated, bzip2- and LZMA-compressed, one entry below a
//! directory entry, two entries, a directory only, an empty archive (its end
//! record alone), a deflated entry whose CRC-32 is wrong in both headers, and
//! eight bytes of a local header followed by junk. `oracle.tsv` is
//! `results/records.tsv` of `../oracle/B8a-zip-input`: the Release build's
//! `FeatureXMLFile::load`/`loadSize` and `MzMLFile::load`/`loadSize` on every
//! fixture and on the two plain inputs, two runs byte-identical on ibminode06.
#![cfg(all(feature = "featurexml", feature = "mzml"))]

use openms::format::featurexml::{self, ReadOptions};
use openms::format::mzml;
use openms::kernel::{Feature, FeatureMap};
use std::path::PathBuf;

const ORACLE: &str = include_str!("data/zip_input/oracle.tsv");

fn fixture(name: &str) -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    match name {
        "featurexml_source_1.featureXML" | "mzml_upstream_serum.mzML" => dir.join(name),
        _ => dir.join("zip_input").join(name),
    }
}

fn featurexml_load(name: &str) -> Result<String, String> {
    let map = featurexml::load(fixture(name)).map_err(|e| e.to_string())?;
    let mut text = map.features.len().to_string();
    if let Some(first) = map.features.first() {
        text.push_str(&format!(
            " rt={:016x} mz={:016x} uid={}",
            first.rt.to_bits(),
            first.mz.to_bits(),
            first.unique_id
        ));
    }
    Ok(text)
}

fn featurexml_load_size(name: &str) -> Result<String, String> {
    featurexml::load_size(fixture(name), &ReadOptions::default())
        .map(|n| n.to_string())
        .map_err(|e| e.to_string())
}

fn mzml_load(name: &str) -> Result<String, String> {
    let experiment = mzml::load(fixture(name)).map_err(|e| e.to_string())?;
    let mut text = format!(
        "{} spectra {} chromatograms",
        experiment.spectra.len(),
        experiment.chromatograms.len()
    );
    if let Some(first) = experiment.spectra.first() {
        text.push_str(&format!(
            " rt0={:016x} peaks0={} id0={}",
            first.rt.to_bits(),
            first.peaks.len(),
            first.native_id
        ));
    }
    Ok(text)
}

fn mzml_load_size(name: &str) -> Result<String, String> {
    mzml::load_size(fixture(name))
        .map(|counts| format!("{} {}", counts.spectra, counts.chromatograms))
        .map_err(|e| e.to_string())
}

#[test]
fn zip_input_matches_the_release_build() {
    let mut compared = 0;
    for line in ORACLE.lines() {
        let mut columns = line.splitn(3, '\t');
        let function = columns.next().unwrap();
        let name = columns.next().unwrap();
        let expected = columns.next().unwrap();
        let actual = match function {
            "FeatureXMLFile::load" => featurexml_load(name),
            "FeatureXMLFile::loadSize" => featurexml_load_size(name),
            "MzMLFile::load" => mzml_load(name),
            "MzMLFile::loadSize" => mzml_load_size(name),
            other => panic!("oracle function {other}"),
        };
        if expected.starts_with("ERR:") {
            assert!(
                actual.is_err(),
                "{function} {name}: C++ {expected}, port {actual:?}"
            );
        } else {
            assert_eq!(actual.as_deref(), Ok(expected), "{function} {name}");
        }
        compared += 1;
    }
    assert_eq!(compared, 52);
}

#[test]
fn zip_refusals_carry_the_source_messages() {
    let two = featurexml::load(fixture("two_entries.zip"))
        .unwrap_err()
        .to_string();
    assert!(
        two.contains(
            "ZIP archive contains 2 file entries; expected exactly 1: a.featureXML b.featureXML in: "
        ),
        "{two}"
    );
    for empty in ["directory_only.zip", "empty_archive.zip"] {
        let error = featurexml::load(fixture(empty)).unwrap_err().to_string();
        assert!(
            error.contains("ZIP archive contains no file entries"),
            "{error}"
        );
    }
    // LZMA: the reference build's libzip has no decoder either, and the
    // oracle's class is FileNotReadable (`ZipIfstream.cpp:94`).
    assert!(matches!(
        featurexml::load(fixture("featurexml_lzma.zip")),
        Err(openms::Error::Io(error)) if error.kind() == std::io::ErrorKind::InvalidData
    ));
    // A wrong CRC-32 is found at the end of the entry, as libzip finds it:
    // `Parse Error:zip_fread returned error in: error reading from zip entry`.
    let crc = featurexml::load(fixture("featurexml_bad_crc.zip"))
        .unwrap_err()
        .to_string();
    assert!(crc.contains("error reading from zip entry"), "{crc}");
}

/// A loaded map with every identification-run identifier replaced by the
/// run's position.
///
/// A run without an `identifier` meta value gets a generated one that ends in
/// a per-process counter (`map_xml::read_run`), so two loads of the same
/// document never agree on it. The position keeps the link from every peptide
/// identification to its run, which is all the identifier encodes.
fn with_run_positions(mut map: FeatureMap) -> FeatureMap {
    fn relabel(features: &mut [Feature], position: &dyn Fn(&mut String)) {
        for feature in features {
            for peptide in &mut feature.peptide_identifications {
                position(&mut peptide.identifier);
            }
            relabel(&mut feature.subordinates, position);
        }
    }
    let names: Vec<String> = map
        .protein_identifications
        .iter()
        .map(|run| run.identifier.clone())
        .collect();
    let position = |identifier: &mut String| {
        let index = names
            .iter()
            .position(|name| name == identifier)
            .unwrap_or_else(|| panic!("peptide identification names no run: {identifier}"));
        *identifier = format!("run {index}");
    };
    for (index, run) in map.protein_identifications.iter_mut().enumerate() {
        run.identifier = format!("run {index}");
    }
    relabel(&mut map.features, &position);
    for peptide in &mut map.unassigned_peptide_identifications {
        position(&mut peptide.identifier);
    }
    map
}

#[test]
fn zip_entry_equals_the_plain_document() {
    let plain =
        with_run_positions(featurexml::load(fixture("featurexml_source_1.featureXML")).unwrap());
    assert!(!plain.protein_identifications.is_empty());
    for name in [
        "featurexml_deflate.zip",
        "featurexml_stored.zip",
        "featurexml_bzip2.zip",
        "featurexml_with_directory.zip",
    ] {
        let zipped = with_run_positions(featurexml::load(fixture(name)).unwrap());
        assert_eq!(zipped.features, plain.features, "{name}");
        assert_eq!(
            zipped.protein_identifications, plain.protein_identifications,
            "{name}"
        );
        assert_eq!(
            zipped.unassigned_peptide_identifications, plain.unassigned_peptide_identifications,
            "{name}"
        );
    }
    let plain = mzml::load(fixture("mzml_upstream_serum.mzML")).unwrap();
    let zipped = mzml::load(fixture("mzml_deflate.zip")).unwrap();
    assert_eq!(zipped.spectra, plain.spectra);
}
