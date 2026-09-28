// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `src/format/pepxml_mascot.rs` (`FORMAT/PepXMLFileMascot.h`).
//!
//! Evidence, strongest first (see `docs/PEPXML_MASCOT_SUPPORT.md` and
//! `tests/data/pepxml_mascot_provenance.json`):
//!
//! - tier 1, executed differential: `../oracle/F3a-pepxml-mascot` loaded the
//!   upstream class-test fixture and 27 generated edge cases with the Release
//!   `libOpenMS` of `openms4-release-bc9cc12-c19e494-174b576` on ibminode06,
//!   twice (identical); `tests/data/pepxml_mascot/expected/` holds its
//!   stdout, stderr and exit codes, compared here title by title, sequence by
//!   sequence (`toString` and `toUnmodifiedString`), log line by log line,
//!   and exception by exception;
//! - tier 3: the one `START_SECTION` of `PepXMLFileMascot_test.cpp`, with its
//!   literals;
//! - tier 4: the resource ceilings.
#![cfg(feature = "idxml")]

use openms::Error;
use openms::format::pepxml_mascot::{MascotPepXmlDocument, PepXmlFileMascot, ReadLimits};
use std::fs;
use std::path::{Path, PathBuf};

fn data() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/pepxml_mascot")
}

fn input(name: &str) -> PathBuf {
    data().join("inputs").join(format!("{name}.pepXML"))
}

fn sequences(document: &MascotPepXmlDocument, title: &str) -> Vec<(String, String)> {
    document.peptides[title]
        .iter()
        .map(|sequence| (sequence.to_string(), sequence.as_str().to_owned()))
        .collect()
}

/// `START_SECTION(void load(const std::string& filename, std::map<...>& peptides))`
/// of `PepXMLFileMascot_test.cpp`, literal by literal.
#[test]
fn class_test_load() {
    let document = PepXmlFileMascot::new().load(input("m01_upstream")).unwrap();
    assert_eq!(document.peptides.len(), 3);
    let titles: Vec<&str> = document.peptides.keys().map(String::as_str).collect();
    assert_eq!(
        titles,
        [
            "402.42_14_3688.53",
            "404.875_14_333.442",
            "411.766_14_3724.98"
        ]
    );
    let expect = |title: &str, expected: &[(&str, &str)]| {
        let got = sequences(&document, title);
        let expected: Vec<(String, String)> = expected
            .iter()
            .map(|(full, bare)| ((*full).to_owned(), (*bare).to_owned()))
            .collect();
        assert_eq!(got, expected, "{title}");
    };
    expect(
        "402.42_14_3688.53",
        &[
            ("LAPSAAEDGAFR", "LAPSAAEDGAFR"),
            ("GQPLGQLEAHR", "GQPLGQLEAHR"),
            ("SPAPLVPVRLR", "SPAPLVPVRLR"),
            ("Q(Deamidated)PYLHPSFSK", "QPYLHPSFSK"),
        ],
    );
    expect(
        "404.875_14_333.442",
        &[
            ("AQAVAAEIRER", "AQAVAAEIRER"),
            ("AALNAADIVTVR", "AALNAADIVTVR"),
            ("AEPAAELALEAK", "AEPAAELALEAK"),
            ("LANAASPAITQ(Deamidated)R", "LANAASPAITQR"),
            ("SGHSAVLLQ(Deamidated)DGK", "SGHSAVLLQDGK"),
            ("AGAGIAN(Deamidated)VQ(Deamidated)AAIR", "AGAGIANVQAAIR"),
            ("VTAPAARSAALGK", "VTAPAARSAALGK"),
        ],
    );
    expect(
        "411.766_14_3724.98",
        &[
            ("LLAWMGRTER", "LLAWMGRTER"),
            ("VLALYRAAQ(Deamidated)AR", "VLALYRAAQAR"),
            ("RTLLMSLTGLK", "RTLLMSLTGLK"),
            ("LLGLSRFGLQ(Deamidated)K", "LLGLSRFGLQK"),
            ("M(Oxidation)GGIALLDEIGK", "MGGIALLDEIGK"),
            ("DQMDN(Deamidated)ALRIR", "DQMDNALRIR"),
            ("Q(Deamidated)TLAGRMVVQ(Deamidated)K", "QTLAGRMVVQK"),
        ],
    );
    assert!(document.log.is_empty());
}

/// The oracle's record of one load: the lines after `load <file>` up to the
/// next load.
struct Recorded {
    status: String,
    body: Vec<String>,
}

fn recorded(case: &str) -> Vec<Recorded> {
    let text = fs::read_to_string(data().join("expected").join(format!("{case}.out"))).unwrap();
    let mut loads: Vec<Recorded> = Vec::new();
    for line in text.lines() {
        if line.starts_with("load ") {
            loads.push(Recorded {
                status: String::new(),
                body: Vec::new(),
            });
        } else if let Some(status) = line.strip_prefix("status ") {
            loads.last_mut().unwrap().status = status.to_owned();
        } else {
            loads.last_mut().unwrap().body.push(line.to_owned());
        }
    }
    let exit = fs::read_to_string(data().join("expected").join(format!("{case}.exit"))).unwrap();
    assert_eq!(exit.trim(), "0", "{case}: the driver itself failed");
    loads
}

/// The handler's own log lines in the oracle's stderr (other lines are the
/// fatal-error echo and the modification registry's own diagnostics).
fn recorded_log(case: &str, file: &str) -> Vec<String> {
    let text = fs::read_to_string(data().join("expected").join(format!("{case}.err"))).unwrap();
    text.lines()
        .filter(|line| line.starts_with("Non-fatal error while loading"))
        .map(|line| line.replace(&format!("'{file}'"), "'<file>'"))
        .collect()
}

/// The oracle driver's rendering of a successful load.
fn render(document: &MascotPepXmlDocument) -> Vec<String> {
    let mut out = vec![format!("peptides {}", document.peptides.len())];
    for (title, list) in &document.peptides {
        out.push(format!("title [{title}] {}", list.len()));
        for sequence in list {
            out.push(format!("seq [{sequence}] [{}]", sequence.as_str()));
        }
    }
    out
}

/// The crate's error kind for each C++ exception the oracle recorded.
fn same_kind(cpp: &str, error: &Error) -> bool {
    let name = cpp.split(':').next().unwrap_or_default();
    match name {
        // `fatalError`, `fromString` and `StringUtils::toDouble`.
        "exception Parse Error" | "exception ConversionError" => {
            matches!(error, Error::Parse { .. })
        }
        "exception IndexOverflow" => matches!(error, Error::InvalidRange(_)),
        "exception InvalidValue" => matches!(error, Error::InvalidValue(_)),
        _ => false,
    }
}

/// Where the port is stricter than the source, with the source's result
/// recorded in the expected files: an empty peptide cannot hold a terminal
/// modification in [`openms::chemistry::AASequence`].
const STRICTER: &[&str] = &["m19_empty_peptide_fixed_nterm"];

fn check_case(case: &str) {
    let file = format!("{case}.pepXML");
    let loads = recorded(case);
    assert_eq!(loads.len(), 1, "{case}");
    let expected = &loads[0];
    let result = PepXmlFileMascot::new().load(input(case));
    if STRICTER.contains(&case) {
        assert_eq!(expected.status, "ok", "{case}");
        assert!(
            matches!(result, Err(Error::InvalidValue(_))),
            "{case}: {result:?}"
        );
        return;
    }
    match result {
        Ok(document) => {
            assert_eq!(expected.status, "ok", "{case}: the source threw");
            assert_eq!(render(&document), expected.body, "{case}");
            let path = input(case).display().to_string();
            let log: Vec<String> = document
                .log
                .iter()
                .map(|line| line.replace(&format!("'{path}'"), "'<file>'"))
                .collect();
            assert_eq!(log, recorded_log(case, &file), "{case}");
            assert_eq!(document.suppressed_log_lines, 0);
        }
        Err(error) => {
            assert!(
                same_kind(&expected.status, &error),
                "{case}: source {:?}, port {error:?}",
                expected.status
            );
            // The source's out-parameter is cleared before parsing and only
            // filled by completed queries; a throw leaves it empty in every
            // recorded case, which is what an error return means here.
            assert_eq!(expected.body, ["peptides 0"], "{case}");
        }
    }
}

const CASES: &[&str] = &[
    "m01_upstream",
    "m02_var_nterm",
    "m03_var_cterm",
    "m04_fixed",
    "m05_fixed_bad",
    "m06_unmatched",
    "m07_duplicates_and_orphans",
    "m08_position_zero",
    "m09_position_past_end",
    "m10_position_negative",
    "m11_unknown_name",
    "m12_missing_peptide",
    "m13_bad_mass",
    "m14_both_missing",
    "m15_missing_variable",
    "m16_three_words",
    "m17_first_match_and_window",
    "m18_empty_peptide",
    "m19_empty_peptide_fixed_nterm",
    "m20_permissive",
    "m21_invalid_sequence",
    "m22_lowercase_variable",
    "m23_variable_then_fixed",
    "m24_empty",
    "m25_position_text",
    "m26_var_protein_nterm",
    "m27_title_order",
];

#[test]
fn release_differential() {
    for case in CASES {
        check_case(case);
    }
}

/// The one missing-attribute case where two attributes are missing: the
/// Release build names `mass`, the argument its compiler evaluates first.
#[test]
fn both_missing_names_mass_first() {
    let error = PepXmlFileMascot::new()
        .load(input("m14_both_missing"))
        .unwrap_err();
    let recorded = &recorded("m14_both_missing")[0].status;
    assert!(recorded.contains("Required attribute 'mass' not present!"));
    match error {
        Error::Parse { message, .. } => {
            assert_eq!(message, "Required attribute 'mass' not present!");
        }
        other => panic!("{other:?}"),
    }
}

/// The source keeps its modification definitions after a load that threw,
/// so its second load matched against the first file's `Oxidation (M)`
/// (`m28_reuse`); every read here starts empty and gives what a fresh object
/// gives (`m28_fresh`).
#[test]
fn a_read_after_a_failed_read_starts_empty() {
    let reader = PepXmlFileMascot::new();
    assert!(reader.load(input("m28_reuse_first_throws")).is_err());
    let document = reader.load(input("m28_reuse_second")).unwrap();
    let fresh = &recorded("m28_fresh")[0];
    assert_eq!(fresh.status, "ok");
    assert_eq!(render(&document), fresh.body);
    let reused = recorded("m28_reuse");
    assert_eq!(reused.len(), 2);
    assert_ne!(reused[1].body, fresh.body, "the recorded source defect");
    assert_eq!(document.log.len(), 1);
    assert!(document.log[0].ends_with("Cannot parse modification '@1'"));
}

#[test]
fn read_from_memory_logs_an_empty_file_name() {
    let text = fs::read(input("m06_unmatched")).unwrap();
    let document = PepXmlFileMascot::new().read(text.as_slice()).unwrap();
    assert_eq!(
        document.log,
        ["Non-fatal error while loading '': Cannot parse modification '@2'"]
    );
}

#[test]
fn a_missing_file_is_an_io_error() {
    let error = PepXmlFileMascot::new()
        .load(data().join("no_such_file.pepXML"))
        .unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error:?}");
}

#[test]
fn ceilings_refuse_instead_of_truncating() {
    let reader = PepXmlFileMascot::new();
    let upstream = fs::read(input("m01_upstream")).unwrap();
    let refused = |limits: ReadLimits| {
        let error = reader
            .read_with_limits(upstream.as_slice(), &limits)
            .unwrap_err();
        assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    };
    let defaults = ReadLimits::default();
    refused(ReadLimits {
        max_bytes: 1_000,
        ..defaults
    });
    refused(ReadLimits {
        max_events: 50,
        ..defaults
    });
    refused(ReadLimits {
        max_depth: 3,
        ..defaults
    });
    refused(ReadLimits {
        max_modification_definitions: 3,
        ..defaults
    });
    refused(ReadLimits {
        max_modifications_per_hit: 1,
        ..defaults
    });
    refused(ReadLimits {
        max_residues: 11,
        ..defaults
    });
    refused(ReadLimits {
        max_sequences: 17,
        ..defaults
    });
    // Exactly at every ceiling the fixture reads (18 hits, 4 definitions,
    // at most 2 modifications per hit, 13 residues).
    let exact = ReadLimits {
        max_modification_definitions: 4,
        max_modifications_per_hit: 2,
        max_residues: 13,
        max_sequences: 18,
        ..defaults
    };
    let document = reader
        .read_with_limits(upstream.as_slice(), &exact)
        .unwrap();
    assert_eq!(document.peptides.values().map(Vec::len).sum::<usize>(), 18);
}

#[test]
fn log_lines_past_the_ceiling_are_counted() {
    let text = fs::read(input("m05_fixed_bad")).unwrap();
    let limits = ReadLimits {
        max_log_lines: 1,
        ..ReadLimits::default()
    };
    let document = PepXmlFileMascot::new()
        .read_with_limits(text.as_slice(), &limits)
        .unwrap();
    assert_eq!(document.log.len(), 1);
    assert_eq!(document.suppressed_log_lines, 2);
}
