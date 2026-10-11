// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `MascotXMLFile::initializeLookup` with a caller-supplied `scan_regex`
//! (`src/format/mascot_xml.rs`, `TitleReferenceFormat::Regex`).
//!
//! Tier 1, executed differential: an oracle driver linked against the Linux
//! x86_64 Release libOpenMS of openms4-release-bc9cc12-c19e494-174b576 ran on
//! ibminode06, twice, identically (`../oracle/F3a-mascot-regex`). Each case
//! builds an experiment from `inputs/spectra_*.tsv`, calls `initializeLookup`
//! with the case's expression (or the default), loads an upstream Mascot XML
//! fixture and prints the retention time and m/z of every identification;
//! the handler's warnings went to stderr. `tests/data/mascot_xml_regex/`
//! holds the records; see `tests/data/mascot_xml_regex_provenance.json`.
#![cfg(feature = "idxml")]

use openms::Error;
use openms::format::file_info::text_format::ostream_g;
use openms::format::mascot_xml::{MascotXmlFile, TitleReferenceFormat};
use openms::{MSExperiment, MSSpectrum, Precursor};
use std::fs;
use std::path::PathBuf;

fn data() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/mascot_xml_regex")
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// The driver's experiment: one spectrum per line of native ID, RT,
/// precursor m/z, precursor charge and MS level, parsed with `std::stod` and
/// `std::stoi` (plain decimal text here).
fn experiment(name: &str) -> MSExperiment {
    let mut experiment = MSExperiment::default();
    if name == "-" {
        return experiment;
    }
    let text = fs::read_to_string(data().join("inputs").join(name)).unwrap();
    for line in text.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let mut spectrum = MSSpectrum {
            native_id: fields[0].to_owned(),
            rt: fields[1].parse().unwrap(),
            ms_level: fields[4].parse().unwrap(),
            ..MSSpectrum::default()
        };
        spectrum.precursors.push(Precursor {
            mz: fields[2].parse().unwrap(),
            charge: fields[3].parse().unwrap(),
            ..Precursor::default()
        });
        experiment.spectra.push(spectrum);
    }
    experiment
}

fn g17(value: Option<f64>) -> String {
    value.map_or_else(|| "none".to_owned(), |value| ostream_g(value, 17))
}

#[test]
fn release_scan_regex_differential() {
    let cases = fs::read_to_string(data().join("inputs/cases.tsv")).unwrap();
    let mut count = 0;
    for row in cases.lines() {
        let fields: Vec<&str> = row.split('\t').collect();
        let (name, mascot, spectra, regex) = (fields[0], fields[1], fields[2], fields[3]);
        let regex = (regex != "<default>").then_some(regex);
        let recorded =
            fs::read_to_string(data().join("expected").join(format!("{name}.out"))).unwrap();
        let stderr =
            fs::read_to_string(data().join("expected").join(format!("{name}.err"))).unwrap();
        let lines: Vec<&str> = recorded.lines().collect();
        count += 1;

        let lookup = MascotXmlFile::initialize_lookup(&experiment(spectra), regex);
        let (lookup, _) = match lookup {
            Ok(lookup) => {
                assert!(lines[0].starts_with("lookup ok"), "{name}: {}", lines[0]);
                assert_eq!(
                    lines[0],
                    format!("lookup ok empty={}", u8::from(lookup.0.is_empty())),
                    "{name}"
                );
                lookup
            }
            Err(error) => {
                // IllegalArgument (no named group) and boost::regex_error.
                assert!(
                    lines[0].starts_with("lookup exception IllegalArgument")
                        || lines[0].starts_with("lookup std_exception"),
                    "{name}: {}",
                    lines[0]
                );
                assert!(matches!(error, Error::InvalidValue(_)), "{name}: {error:?}");
                continue;
            }
        };
        if let Some(expression) = regex {
            assert_eq!(
                lookup.reference_formats(),
                [TitleReferenceFormat::Regex(0)],
                "{name}: {expression}"
            );
        }
        match MascotXmlFile::new().load(fixture(mascot), &lookup) {
            Ok(result) => {
                assert_eq!(
                    lines[1],
                    format!("load ok {}", result.peptide_identifications.len()),
                    "{name}"
                );
                let rendered: Vec<String> = result
                    .peptide_identifications
                    .iter()
                    .map(|id| format!("rt={} mz={}", g17(id.rt), g17(id.mz)))
                    .collect();
                assert_eq!(rendered, lines[2..], "{name}");
                // One non-fatal error per title the lookup could not resolve.
                let unresolved = |text: &str| text.contains("> element has unexpected format '");
                // The source's log stream collapses a run of identical lines
                // into the first and `<line> occurred N times`, which stands
                // for N - 1 more.
                let recorded: usize = stderr
                    .lines()
                    .filter(|line| unresolved(line))
                    .map(|line| {
                        line.strip_suffix(" times")
                            .and_then(|rest| rest.rsplit_once("> occurred "))
                            .map_or(1, |(_, n)| n.parse::<usize>().unwrap() - 1)
                    })
                    .sum();
                assert_eq!(
                    result.warnings.iter().filter(|w| unresolved(w)).count(),
                    recorded,
                    "{name}"
                );
            }
            Err(error) => {
                assert!(
                    lines[1].starts_with("load exception MissingInformation"),
                    "{name}: {} vs {error:?}",
                    lines[1]
                );
                assert!(
                    matches!(error, Error::MissingInformation(_)),
                    "{name}: {error:?}"
                );
            }
        }
    }
    assert_eq!(count, 14);
}

#[test]
fn a_regex_format_is_registered_after_the_builtin_ones_in_call_order() {
    let mut lookup = MascotXmlFile::initialize_lookup(&MSExperiment::default(), None)
        .unwrap()
        .0;
    let format = lookup
        .add_regex_reference_format(r"_(?<RT>\d+\.\d+)$")
        .unwrap();
    assert_eq!(format, TitleReferenceFormat::Regex(0));
    assert_eq!(
        lookup.reference_formats(),
        [
            TitleReferenceFormat::MzThenRt,
            TitleReferenceFormat::Regex(0)
        ]
    );
    // Nothing is registered when the expression is refused.
    assert!(lookup.add_regex_reference_format("(?<SCAN>").is_err());
    assert_eq!(lookup.reference_formats().len(), 2);
}
