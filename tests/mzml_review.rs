// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

#![cfg(feature = "mzml")]

use openms::format::mzml;
use openms::kernel::SpectrumType;
use openms::{MSExperiment, MSSpectrum, Peak1D, Precursor};

fn fixture() -> String {
    let experiment = MSExperiment {
        spectra: vec![MSSpectrum {
            peaks: vec![Peak1D::new(100.0, 1.0)],
            spectrum_type: SpectrumType::Centroid,
            precursors: vec![Precursor {
                mz: 500.0,
                intensity: 10.0,
                charge: 2,
                ..Precursor::default()
            }],
            ..MSSpectrum::default()
        }],
        ..MSExperiment::default()
    };
    let mut bytes = Vec::new();
    mzml::write(&mut bytes, &experiment).unwrap();
    String::from_utf8(bytes).unwrap()
}

fn after_cv(xml: &str, original_accession: &str, content: &str) -> String {
    let start = xml
        .find(&format!("accession=\"{original_accession}\""))
        .unwrap();
    let end = start + xml[start..].find("/>").unwrap() + 2;
    format!("{}{}{}", &xml[..end], content, &xml[end..])
}

fn inside_run(xml: &str, content: &str) -> String {
    let start = xml.find("<run ").unwrap();
    let end = start + xml[start..].find('>').unwrap() + 1;
    format!("{}{}{}", &xml[..end], content, &xml[end..])
}

#[test]
fn mzml_rejects_conflicting_duplicate_scientific_cv_values() {
    let xml = fixture();
    for (original, duplicate, value) in [
        ("MS:1000511", "MS:1000511", "2"),
        ("MS:1000127", "MS:1000128", ""),
        ("MS:1000744", "MS:1000744", "600"),
        ("MS:1000744", "MS:1000040", "600"),
        ("MS:1000041", "MS:1000041", "3"),
        ("MS:1000042", "MS:1000042", "100"),
    ] {
        let added = format!(
            "<cvParam cvRef=\"MS\" accession=\"{duplicate}\" name=\"review value\" value=\"{value}\"/>"
        );
        let malformed = after_cv(&xml, original, &added);
        assert!(
            mzml::read(malformed.as_bytes()).is_err(),
            "overwrote {original} with {duplicate}"
        );
    }
}

#[test]
fn mzml_reader_and_writer_reserve_internal_name_key_consistently() {
    let xml = inside_run(
        &fixture(),
        "<userParam name=\"openms-rust:name\" value=\"reserved\"/>",
    );
    assert!(mzml::read(xml.as_bytes()).is_err());
}

#[test]
fn mzml_rejects_forbidden_xml_characters_in_ignored_text() {
    for character in ['\u{1}', '\u{b}', '\u{fffe}', '\u{ffff}'] {
        let xml = inside_run(
            &fixture(),
            &format!("<userParam name=\"ignored\">{character}</userParam>"),
        );
        assert!(
            mzml::read(xml.as_bytes()).is_err(),
            "accepted forbidden XML character {character:?}"
        );
    }
}
