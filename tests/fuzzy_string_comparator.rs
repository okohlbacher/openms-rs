// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! C3-FUZZY: the `FuzzyStringComparator` test-support port and the decoded
//! comparator.
//!
//! Evidence, in order of strength (see `tests/data/fuzzy_string_comparator_provenance.json`):
//! - an executed C++ differential over 137 comparator cases, run on the
//!   unmodified `FuzzyStringComparator` object from product-sdk (verdicts and
//!   complete logs);
//! - executed `FuzzyDiff` exit codes for 39 invocations: TOPP_FuzzyDiff_1..4 as
//!   registered and 35 more, including the retained-vs-current-C++ pairs of
//!   FileInfo_3, FileInfo_17 and FeatureFinderCentroided_1;
//! - the 25 `START_SECTION`s of `FuzzyStringComparator_test.cpp`, transcribed or
//!   accounted for.

#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;

#[path = "support/decoded_compare.rs"]
mod decoded;

use fuzzy::{
    FuzzyDiffExit, FuzzyDiffSettings, FuzzyStringComparator, LogDestination, compare_numbers,
    extract_double, format_g, fuzzy_diff, message,
};
use std::collections::BTreeMap;
use std::io::Cursor;
use std::path::{Path, PathBuf};

fn data(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/fuzzy_string_comparator")
        .join(relative)
}

fn read(relative: &str) -> Vec<u8> {
    std::fs::read(data(relative)).unwrap_or_else(|e| panic!("{relative}: {e}"))
}

fn temp_dir(case: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms-fuzzy-{}-{case}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// `StringUtils::split(log, '\n', substrings)`: nothing for an empty string,
/// otherwise every piece including the one after the last newline.
fn split_lines(log: &[u8]) -> Vec<&[u8]> {
    if log.is_empty() {
        Vec::new()
    } else {
        log.split(|&b| b == b'\n').collect()
    }
}

fn buffered(verbose: i32, ratio: f64, absdiff: f64) -> FuzzyStringComparator {
    let mut comparator = FuzzyStringComparator::new();
    comparator.set_log_destination(LogDestination::Buffer);
    comparator.set_verbose_level(verbose);
    comparator.set_acceptable_relative(ratio);
    comparator.set_acceptable_absolute(absdiff);
    comparator
}

// ------------------------------------------------------------------------
// FuzzyStringComparator_test.cpp (core bc9cc12), section by section.
// ------------------------------------------------------------------------

// START_SECTION((FuzzyStringComparator())) and ((virtual ~FuzzyStringComparator()))
#[test]
fn class_test_constructor_and_destructor() {
    let comparator = FuzzyStringComparator::new();
    assert_eq!(comparator.verbose_level(), 2);
    drop(comparator);
}

// START_SECTION((FuzzyStringComparator& operator=(...))) and the copy constructor:
// NOT_TESTABLE in the source (declared, not implemented). The Rust type does not
// implement Clone, so neither operation exists.

// START_SECTIONs getAcceptableAbsolute, getAcceptableRelative, getVerboseLevel,
// getTabWidth, getFirstColumn, getLogDestination: NOT_TESTABLE, "tested along with
// set-method"; they are exercised by the setter tests below.

// START_SECTION((void setAcceptableAbsolute(const double rhs)))
#[test]
fn class_test_set_acceptable_absolute() {
    let mut comparator = FuzzyStringComparator::new();
    comparator.set_acceptable_absolute(2345.6789);
    assert_eq!(comparator.acceptable_absolute(), 2345.6789);
}

// START_SECTION((void setAcceptableRelative(const double rhs)))
#[test]
fn class_test_set_acceptable_relative() {
    let mut comparator = FuzzyStringComparator::new();
    comparator.set_acceptable_relative(6789.2345);
    assert_eq!(comparator.acceptable_relative(), 6789.2345);
}

// START_SECTION((void setTabWidth(const int rhs)))
#[test]
fn class_test_set_tab_width() {
    let mut comparator = FuzzyStringComparator::new();
    comparator.set_tab_width(1452);
    assert_eq!(comparator.tab_width(), 1452);
}

// START_SECTION((void setFirstColumn(const int rhs)))
#[test]
fn class_test_set_first_column() {
    let mut comparator = FuzzyStringComparator::new();
    comparator.set_first_column(4321235);
    assert_eq!(comparator.first_column(), 4321235);
}

// START_SECTION((void setLogDestination(std::ostream & rhs)))
#[test]
fn class_test_set_log_destination() {
    let mut comparator = FuzzyStringComparator::new();
    assert_eq!(comparator.log_destination(), LogDestination::Stdout);
    comparator.set_log_destination(LogDestination::Stderr);
    assert_eq!(comparator.log_destination(), LogDestination::Stderr);
    assert_ne!(comparator.log_destination(), LogDestination::Stdout);
    comparator.set_log_destination(LogDestination::Stdout);
    assert_ne!(comparator.log_destination(), LogDestination::Stderr);
    assert_eq!(comparator.log_destination(), LogDestination::Stdout);
}

// START_SECTION((void setVerboseLevel(const int rhs)))
#[test]
fn class_test_set_verbose_level() {
    let mut comparator = FuzzyStringComparator::new();
    assert_eq!(comparator.verbose_level(), 2);
    comparator.set_verbose_level(88);
    assert_eq!(comparator.verbose_level(), 88);
    comparator.set_verbose_level(-21);
    assert_eq!(comparator.verbose_level(), -21);
}

// START_SECTIONs ((const StringList& getWhitelist() const)), ((StringList&
// getWhitelist())) and ((void setWhitelist(const StringList &rhs))), which share one
// comparator in the source.
#[test]
fn class_test_whitelist_accessors() {
    let mut comparator = FuzzyStringComparator::new();
    assert!(comparator.whitelist().is_empty());
    assert!(comparator.whitelist_mut().is_empty());
    let list = |text: &str| text.split(',').map(str::to_owned).collect::<Vec<_>>();
    comparator.set_whitelist(list("null,eins,zwei,drei"));
    assert_eq!(comparator.whitelist_mut()[0], "null");
    assert_eq!(comparator.whitelist()[1], "eins");
    assert_eq!(comparator.whitelist().len(), 4);
    comparator.set_whitelist(list("zero,one,two,three,four"));
    assert_eq!(comparator.whitelist_mut()[0], "zero");
    assert_eq!(comparator.whitelist()[1], "one");
    assert_eq!(comparator.whitelist().len(), 5);
}

const V1_LHS: &str = "1 \n \t\t   2\t\n 3";
const V1_RHS: &str = "1.01 \n \n\t\t\n\n  \t\t\t\t\t  \t0002.01000 \n 3";
const V3_RHS: &str = "1.11 \n \n\t\t\n\n  \t\t\t\t\t  \t0004.01000 \n 3";
const V4_LHS: &str = "1 \n xx\n 2.008\t\n 3";
const V4_RHS: &str = "1.11 \nU\n\t\t\n\n  q\t\t\t\t\t  \t0002.04000 \n 3";

/// The seven failure reports of the verbose-3 comparison, 35 lines each plus the
/// piece after the final newline.
fn assert_verbose_3_report(log: &[u8]) {
    let lines = split_lines(log);
    assert_eq!(lines.len(), 246, "{}", String::from_utf8_lossy(log));
    let expected = [
        (0, "FAILED: 'ratio of numbers is too large'"),
        (35, "FAILED: 'input_1 is whitespace, but input_2 is not'"),
        (70, "FAILED: 'different letters'"),
        (
            105,
            "FAILED: 'line from input_2 is shorter than line from input_1'",
        ),
        (140, "FAILED: 'input_1 is a number, but input_2 is not'"),
        (175, "FAILED: 'input_1 is not a number, but input_2 is'"),
        (
            210,
            "FAILED: 'line from input_1 is shorter than line from input_2'",
        ),
    ];
    for (index, text) in expected {
        assert_eq!(lines[index], text.as_bytes(), "line {index}");
    }
}

// START_SECTION((bool compareStrings(std::string const &lhs, std::string const &rhs)))
#[test]
fn class_test_compare_strings() {
    // What regular expressions could not do.
    assert!(buffered(2, 1.00021, 0.0).compare_strings("0.9999E4", "1.0001E4"));
    assert!(buffered(2, 1.0, 2.0).compare_strings("0.9999E4", "1.0001E4"));
    // Mixing letters, whitespace and numbers.
    assert!(
        buffered(1, 1.01, 0.001)
            .compare_strings("bl   a b 00.0022 asdfdf", "bl a  b 0.00225 asdfdf")
    );
    assert!(!buffered(1, 1.01, 0.01).compare_strings("bl   a 1.2   b", "bl a 1.25 b"));
    assert!(buffered(1, 2.0, 0.01).compare_strings("bl   a 1.2   b", "bl a 1.25 b"));
    assert!(buffered(1, 1.01, 0.0).compare_strings("bl   a 1.002   b", "bl a 1.0025 b"));

    // The impact of the verbose level.
    let mut quiet = buffered(1, 1.03, 0.01);
    assert!(quiet.compare_strings(V1_LHS, V1_RHS));
    assert!(split_lines(quiet.log()).is_empty());

    let mut summary = buffered(2, 1.03, 0.01);
    assert!(summary.compare_strings(V1_LHS, V1_RHS));
    let lines = split_lines(summary.log());
    assert_eq!(
        lines.len(),
        17,
        "{}",
        String::from_utf8_lossy(summary.log())
    );
    assert_eq!(lines[0], b"PASSED.");

    let mut failure = buffered(1, 1.01, 0.01);
    failure.compare_strings(V1_LHS, V3_RHS);
    let lines = split_lines(failure.log());
    assert_eq!(
        lines.len(),
        36,
        "{}",
        String::from_utf8_lossy(failure.log())
    );
    assert_eq!(lines[0], b"FAILED: 'ratio of numbers is too large'");

    let mut all = buffered(3, 1.01, 0.01);
    all.compare_strings(V4_LHS, V4_RHS);
    assert_verbose_3_report(all.log());

    assert!(buffered(2, 1.0, 2.0).compare_strings("0.9999X", "1.0001X"));
}

// START_SECTION((bool compareStreams(std::istream &input_1, std::istream &input_2)))
#[test]
fn class_test_compare_streams() {
    let mut comparator = buffered(3, 1.01, 0.01);
    comparator.compare_streams(&mut Cursor::new(V4_LHS), &mut Cursor::new(V4_RHS));
    assert_verbose_3_report(comparator.log());
}

// START_SECTION((bool compareFiles(const std::string &filename_1, const std::string &filename_2)))
#[test]
fn class_test_compare_files() {
    let dir = temp_dir("class-files");
    let (file_1, file_2) = (dir.join("1.tmp"), dir.join("2.tmp"));
    std::fs::write(&file_1, V4_LHS).unwrap();
    std::fs::write(&file_2, V4_RHS).unwrap();
    let mut comparator = buffered(3, 1.01, 0.01);
    comparator.compare_files(&file_1, &file_2);
    assert_verbose_3_report(comparator.log());
    std::fs::remove_dir_all(dir).unwrap();
}

// START_SECTION(([EXTRA] non-finite numbers are compared by category, not by arithmetic))
#[test]
fn class_test_non_finite_numbers() {
    let equal = |lhs: &str, rhs: &str| buffered(2, 1.01, 0.01).compare_strings(lhs, rhs);
    assert!(equal("x nan", "x nan"));
    assert!(equal("x inf", "x inf"));
    assert!(equal("x -inf", "x -inf"));
    assert!(!equal("x nan", "x 5.0"));
    assert!(!equal("x 5.0", "x nan"));
    assert!(!equal("x nan", "x inf"));
    assert!(!equal("x nan", "x -inf"));
    assert!(!equal("x nan", "x 0.0"));
    assert!(!equal("x inf", "x -inf"));
    assert!(!equal("x -inf", "x inf"));
    assert!(!equal("x inf", "x 5.0"));
    assert!(!equal("x -inf", "x 5.0"));
    assert!(equal("x 1.000 nan", "x 1.001 nan"));
    assert!(!equal("x 1.0 nan", "x 2.0 nan"));
    assert!(!equal("x nan", "x abc"));
}

// The two remaining START_SECTIONs, reportFailure_ and reportSuccess_, are commented
// out in the source ("Tested in compare...() methods"); both reports are checked
// line for line by the executed differential below.

// ------------------------------------------------------------------------
// Executed C++ differential: 137 comparator cases (../oracle/fuzzy-string-comparator).
// ------------------------------------------------------------------------

struct OracleCase {
    name: String,
    kind: String,
    mode: String,
    ratio: f64,
    absdiff: f64,
    verbose: i32,
    tab_width: i32,
    first_column: i32,
    whitelist: Vec<String>,
    matched: Vec<(String, String)>,
    lhs: Vec<u8>,
    rhs: Vec<u8>,
    lhs2: Vec<u8>,
    rhs2: Vec<u8>,
}

struct OracleResult {
    verdict_1: Option<bool>,
    verdict_2: Option<bool>,
    cwd: String,
    log: Vec<u8>,
}

fn unhex(text: &str) -> Vec<u8> {
    assert_eq!(text.len() % 2, 0, "odd hex length");
    (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).unwrap())
        .collect()
}

fn hex_list(field: &str) -> Vec<String> {
    let mut parts = field.split('|');
    let count: usize = parts.next().unwrap().parse().unwrap();
    let items: Vec<String> = parts
        .map(|h| String::from_utf8(unhex(h)).unwrap())
        .collect();
    assert_eq!(items.len(), count);
    items
}

fn tsv_rows(relative: &str) -> Vec<Vec<String>> {
    String::from_utf8(read(relative))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect()
}

fn oracle_cases() -> Vec<OracleCase> {
    tsv_rows("oracle/comparator_cases.tsv")
        .into_iter()
        .map(|f| {
            assert_eq!(f.len(), 14);
            let matched = {
                let mut parts = f[9].split('|');
                let count: usize = parts.next().unwrap().parse().unwrap();
                let pairs: Vec<(String, String)> = parts
                    .map(|item| {
                        let (a, b) = item.split_once(':').unwrap();
                        (
                            String::from_utf8(unhex(a)).unwrap(),
                            String::from_utf8(unhex(b)).unwrap(),
                        )
                    })
                    .collect();
                assert_eq!(pairs.len(), count);
                pairs
            };
            OracleCase {
                name: f[0].clone(),
                kind: f[1].clone(),
                mode: f[2].clone(),
                ratio: f[3].parse().unwrap(),
                absdiff: f[4].parse().unwrap(),
                verbose: f[5].parse().unwrap(),
                tab_width: f[6].parse().unwrap(),
                first_column: f[7].parse().unwrap(),
                whitelist: hex_list(&f[8]),
                matched,
                lhs: unhex(&f[10]),
                rhs: unhex(&f[11]),
                lhs2: unhex(&f[12]),
                rhs2: unhex(&f[13]),
            }
        })
        .collect()
}

fn oracle_results() -> BTreeMap<String, OracleResult> {
    let verdict = |text: &str| match text {
        "1" => Some(true),
        "0" => Some(false),
        "-" => None,
        other => panic!("bad verdict {other}"),
    };
    tsv_rows("oracle/comparator_results.tsv")
        .into_iter()
        .map(|f| {
            assert_eq!(f.len(), 5);
            (
                f[0].clone(),
                OracleResult {
                    verdict_1: verdict(&f[1]),
                    verdict_2: verdict(&f[2]),
                    cwd: f[3].clone(),
                    log: unhex(&f[4]),
                },
            )
        })
        .collect()
}

fn replace_all(haystack: &[u8], needle: &[u8], replacement: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(haystack.len());
    let mut i = 0;
    while i < haystack.len() {
        if !needle.is_empty() && haystack[i..].starts_with(needle) {
            out.extend_from_slice(replacement);
            i += needle.len();
        } else {
            out.push(haystack[i]);
            i += 1;
        }
    }
    out
}

/// Replace the `numbers:` value of every element whose `is_number:` is false. The
/// value of a non-number is meaningless: the source resets it to NaN, and the libc++
/// build's `strtod` fallback then overwrites it with `strtod`'s result for the
/// letter (usually 0); `std::from_chars`, whose contract this port follows, leaves NaN.
fn mask_non_number_values(log: &[u8]) -> Vec<u8> {
    let lines: Vec<&[u8]> = log.split(|&b| b == b'\n').collect();
    let mut out: Vec<Vec<u8>> = Vec::with_capacity(lines.len());
    let mut is_number = [true, true];
    for line in lines {
        if let Some(rest) = line.strip_prefix(b"  is_number:\t") {
            let mut fields = rest.split(|&b| b == b'\t');
            is_number = [
                fields.next() == Some(b"true"),
                fields.next() == Some(b"true"),
            ];
            out.push(line.to_vec());
        } else if let Some(rest) = line.strip_prefix(b"  numbers:\t") {
            let fields: Vec<&[u8]> = rest.split(|&b| b == b'\t').collect();
            let mut masked = b"  numbers:".to_vec();
            for (index, field) in fields.iter().enumerate() {
                masked.push(b'\t');
                if is_number.get(index).copied().unwrap_or(true) {
                    masked.extend_from_slice(field);
                } else {
                    masked.extend_from_slice(b"<not a number>");
                }
            }
            out.push(masked);
        } else {
            out.push(line.to_vec());
        }
    }
    out.join(&b'\n')
}

/// Cases where this port deliberately follows the documented `std::from_chars`
/// contract and the libc++ oracle's `strtod` fallback does not.
const HEX_FLOAT_DIVERGENCE: &str = "tok_hex_vs_decimal";

#[test]
fn executed_cpp_comparator_corpus_verdicts_and_logs_match() {
    let cases = oracle_cases();
    let results = oracle_results();
    assert_eq!(cases.len(), 137);
    assert_eq!(results.len(), cases.len());
    let root = temp_dir("oracle-corpus");
    let rust_cwd = std::env::current_dir()
        .unwrap()
        .to_string_lossy()
        .into_owned();
    let mut failures = Vec::new();
    for case in &cases {
        let expected = results
            .get(&case.name)
            .unwrap_or_else(|| panic!("no result for {}", case.name));
        let mut comparator = FuzzyStringComparator::new();
        comparator.set_log_destination(LogDestination::Buffer);
        comparator.set_verbose_level(case.verbose);
        comparator.set_acceptable_relative(case.ratio);
        comparator.set_acceptable_absolute(case.absdiff);
        comparator.set_tab_width(case.tab_width);
        comparator.set_first_column(case.first_column);
        comparator.set_whitelist(case.whitelist.clone());
        comparator.set_matched_whitelist(case.matched.clone());
        let mut verdict_2 = None;
        let (verdict_1, rust_prefix, cpp_prefix, replacement): (bool, String, String, &[u8]) =
            match case.kind.as_str() {
                "S" => (
                    comparator.compare_bytes(&case.lhs, &case.rhs),
                    format!("{rust_cwd}/"),
                    format!("{}/", expected.cwd),
                    b"<CWD>/",
                ),
                "L" => (
                    comparator.compare_lines(&case.lhs, &case.rhs),
                    format!("{rust_cwd}/"),
                    format!("{}/", expected.cwd),
                    b"<CWD>/",
                ),
                "T" => {
                    let first = comparator.compare_bytes(&case.lhs, &case.rhs);
                    verdict_2 = Some(comparator.compare_bytes(&case.lhs2, &case.rhs2));
                    (
                        first,
                        format!("{rust_cwd}/"),
                        format!("{}/", expected.cwd),
                        b"<CWD>/",
                    )
                }
                "F" => {
                    let dir = root.join(&case.name);
                    std::fs::create_dir_all(&dir).unwrap();
                    std::fs::write(dir.join("a.txt"), &case.lhs).unwrap();
                    std::fs::write(dir.join("b.txt"), &case.rhs).unwrap();
                    let (a, b) = match case.mode.as_str() {
                        "files" => ("a.txt", "b.txt"),
                        "same" => ("a.txt", "a.txt"),
                        "missing1" => ("missing.txt", "b.txt"),
                        "missing2" => ("a.txt", "missing.txt"),
                        other => panic!("mode {other}"),
                    };
                    let verdict = comparator.compare_files(&dir.join(a), &dir.join(b));
                    (
                        verdict,
                        format!("{}/", dir.to_string_lossy()),
                        format!("{}/", expected.cwd),
                        b"",
                    )
                }
                other => panic!("kind {other}"),
            };
        if case.name == HEX_FLOAT_DIVERGENCE {
            assert_eq!(
                expected.verdict_1,
                Some(true),
                "libc++ strtod reads 0x10 as 16"
            );
            assert!(
                !verdict_1,
                "std::from_chars reads 0x10 as 0 followed by letters"
            );
            continue;
        }
        if Some(verdict_1) != expected.verdict_1 || verdict_2 != expected.verdict_2 {
            failures.push(format!(
                "{}: verdict {verdict_1}/{verdict_2:?}, C++ {:?}/{:?}",
                case.name, expected.verdict_1, expected.verdict_2
            ));
            continue;
        }
        let rust_log = mask_non_number_values(&replace_all(
            comparator.log(),
            rust_prefix.as_bytes(),
            replacement,
        ));
        let cpp_log = mask_non_number_values(&replace_all(
            &expected.log,
            cpp_prefix.as_bytes(),
            replacement,
        ));
        if rust_log != cpp_log {
            failures.push(format!(
                "{}: log differs\n--- Rust\n{}\n--- C++\n{}",
                case.name,
                String::from_utf8_lossy(&rust_log),
                String::from_utf8_lossy(&cpp_log)
            ));
        }
    }
    std::fs::remove_dir_all(root).unwrap();
    assert!(
        failures.is_empty(),
        "{} of {} cases differ:\n{}",
        failures.len(),
        cases.len(),
        failures.join("\n")
    );
}

#[test]
fn decoded_number_rule_agrees_with_executed_single_number_cases() {
    let results = oracle_results();
    let mut checked = 0;
    for case in oracle_cases() {
        if case.kind != "S" || !case.whitelist.is_empty() || !case.matched.is_empty() {
            continue;
        }
        let whole = |bytes: &[u8]| {
            extract_double(bytes)
                .filter(|(_, n)| *n == bytes.len())
                .map(|(v, _)| v)
        };
        let (Some(a), Some(b)) = (whole(&case.lhs), whole(&case.rhs)) else {
            continue;
        };
        let tolerance = decoded::Tolerance::new(case.ratio, case.absdiff);
        assert_eq!(
            Some(tolerance.check(a, b).is_ok()),
            results[&case.name].verdict_1,
            "{}",
            case.name
        );
        checked += 1;
    }
    assert!(checked >= 25, "only {checked} single-number cases");
}

// ------------------------------------------------------------------------
// FuzzyDiff tool contract: executed exit codes (TOPP_FuzzyDiff_1..4 and more).
// ------------------------------------------------------------------------

fn upstream() -> FuzzyDiffSettings {
    let settings = FuzzyDiffSettings::upstream().unwrap();
    assert!(settings.ini_errors.is_empty());
    settings
}

fn invocation(case: &str) -> (FuzzyDiffSettings, PathBuf, PathBuf) {
    let fd = |name: &str| data(&format!("fuzzydiff/{name}"));
    let retained = |name: &str| data(&format!("retained/{name}"));
    let current = |name: &str| data(&format!("oracle/{name}"));
    let ini = |name: &str| FuzzyDiffSettings::load_ini(&fd(name)).unwrap();
    let file_name = || upstream().with_whitelist(&["File name"]);
    let id = || upstream().with_whitelist(&["id="]);
    let with = |mut settings: FuzzyDiffSettings, edit: fn(&mut FuzzyDiffSettings)| {
        edit(&mut settings);
        settings
    };
    match case {
        "topp_fuzzydiff_1" => (
            upstream(),
            retained("FuzzyDiff_1_in1.featureXML"),
            retained("FuzzyDiff_1_in1.featureXML"),
        ),
        "topp_fuzzydiff_2" => (
            upstream(),
            retained("FuzzyDiff_1_in1.featureXML"),
            retained("FuzzyDiff_1_in2.featureXML"),
        ),
        "topp_fuzzydiff_3" => (
            upstream(),
            retained("FuzzyDiff_3_in1.featureXML"),
            retained("FuzzyDiff_3_in2.featureXML"),
        ),
        "topp_fuzzydiff_4" => (
            upstream(),
            retained("lorem_ipsum.featureXML"),
            retained("FuzzyDiff_3_in2.featureXML"),
        ),
        "fileinfo_3_current_vs_retained" => (
            file_name(),
            current("FileInfo_3.current.txt"),
            retained("FileInfo_3_output.txt"),
        ),
        "fileinfo_17_current_vs_retained" => (
            file_name(),
            current("FileInfo_17.current.txt"),
            retained("FileInfo_17_output.txt"),
        ),
        "ffc_1_current_vs_retained" => (
            id(),
            current("FeatureFinderCentroided_1.current.featureXML"),
            retained("FeatureFinderCentroided_1_1_output.featureXML"),
        ),
        "fileinfo_3_without_whitelist" => (
            upstream(),
            current("FileInfo_3.current.txt"),
            retained("FileInfo_3_output.txt"),
        ),
        "ffc_1_without_whitelist" => (
            upstream(),
            current("FeatureFinderCentroided_1.current.featureXML"),
            retained("FeatureFinderCentroided_1_1_output.featureXML"),
        ),
        "fileinfo_3_beyond" => (
            file_name(),
            current("FileInfo_3.current.txt"),
            fd("FileInfo_3_output.beyond.txt"),
        ),
        "fileinfo_3_within" => (
            file_name(),
            current("FileInfo_3.current.txt"),
            fd("FileInfo_3_output.within.txt"),
        ),
        "fileinfo_17_beyond" => (
            file_name(),
            current("FileInfo_17.current.txt"),
            fd("FileInfo_17_output.beyond.txt"),
        ),
        "ffc_1_beyond" => (
            id(),
            current("FeatureFinderCentroided_1.current.featureXML"),
            fd("FeatureFinderCentroided_1_1_output.beyond.featureXML"),
        ),
        "ffc_1_within" => (
            id(),
            current("FeatureFinderCentroided_1.current.featureXML"),
            fd("FeatureFinderCentroided_1_1_output.within.featureXML"),
        ),
        "ffc_1_id_change" => (
            id(),
            current("FeatureFinderCentroided_1.current.featureXML"),
            fd("FeatureFinderCentroided_1_1_output.id_change.featureXML"),
        ),
        "defaults_without_ini" => (
            FuzzyDiffSettings::registered_defaults(),
            fd("x_1.txt"),
            fd("x_1001.txt"),
        ),
        "ini_tolerance" => (upstream(), fd("x_1.txt"), fd("x_1001.txt")),
        "empty_in1" => (upstream(), fd("empty.txt"), fd("x_1.txt")),
        "empty_in2" => (upstream(), fd("x_1.txt"), fd("empty.txt")),
        "missing_in1" => (upstream(), fd("missing.txt"), fd("x_1.txt")),
        "missing_in1_empty_in2" => (upstream(), fd("missing.txt"), fd("empty.txt")),
        "matched_single_token" => (
            upstream().with_matched_whitelist(&["abc"]),
            fd("alpha.txt"),
            fd("beta.txt"),
        ),
        "matched_three_tokens" => (
            upstream().with_matched_whitelist(&["a:b:c"]),
            fd("alpha.txt"),
            fd("beta.txt"),
        ),
        "matched_pair" => (
            upstream().with_matched_whitelist(&["alpha:beta"]),
            fd("alpha.txt"),
            fd("beta.txt"),
        ),
        "matched_trailing_colon" => (
            upstream().with_matched_whitelist(&["a:"]),
            fd("a_1.txt"),
            fd("z_2.txt"),
        ),
        "ratio_below_min" => (
            with(upstream(), |s| s.ratio = 0.5),
            fd("x_1.txt"),
            fd("x_1001.txt"),
        ),
        "verbose_above_max" => (
            with(upstream(), |s| s.verbose = 4),
            fd("x_1.txt"),
            fd("x_1001.txt"),
        ),
        "tab_width_zero" => (
            with(upstream(), |s| s.tab_width = 0),
            fd("x_1.txt"),
            fd("x_1001.txt"),
        ),
        "ini_ratio_below_min" => (
            ini("FuzzyDiff_ratio_below_min.ini"),
            fd("x_1.txt"),
            fd("x_1001.txt"),
        ),
        "ini_absdiff_negative" => (
            ini("FuzzyDiff_absdiff_negative.ini"),
            fd("x_1.txt"),
            fd("x_1001.txt"),
        ),
        "ini_first_column_negative" => (
            ini("FuzzyDiff_first_column_negative.ini"),
            fd("x_1.txt"),
            fd("x_1001.txt"),
        ),
        "ini_unknown_item" => (
            ini("FuzzyDiff_unknown_item.ini"),
            fd("x_1.txt"),
            fd("x_1001.txt"),
        ),
        "ini_verbose_3" => (ini("FuzzyDiff_verbose_3.ini"), fd("x_1.txt"), fd("a_1.txt")),
        "sort_reordered_rows" => (
            with(upstream(), |s| s.sort = true),
            fd("table_a.tsv"),
            fd("table_b.tsv"),
        ),
        "unsorted_reordered_rows" => (upstream(), fd("table_a.tsv"), fd("table_b.tsv")),
        "sort_same_file" => (
            with(upstream(), |s| s.sort = true),
            fd("table_a.tsv"),
            fd("table_a.tsv"),
        ),
        "ini_whitelist_stylesheet" => (upstream(), fd("stylesheet_a.xml"), fd("stylesheet_b.xml")),
        "cli_whitelist_replaces_ini" => (id(), fd("stylesheet_a.xml"), fd("stylesheet_b.xml")),
        "trailing_space" => (
            upstream(),
            fd("trailing_space.txt"),
            fd("no_trailing_space.txt"),
        ),
        other => panic!("oracle FuzzyDiff case '{other}' has no Rust invocation"),
    }
}

#[test]
fn executed_fuzzydiff_exit_codes_match() {
    let rows = tsv_rows("oracle/tool_runs.tsv");
    let mut checked = 0;
    for row in &rows {
        let Some(case) = row[0].strip_prefix("fd_") else {
            // FileInfo and FeatureFinderCentroided runs that produced the current outputs.
            assert_eq!(row[1], "0", "{}", row[0]);
            continue;
        };
        let expected: i32 = row[1].parse().unwrap();
        let (settings, in1, in2) = invocation(case);
        let outcome = fuzzy_diff(&in1, &in2, &settings);
        assert_eq!(
            outcome.exit.code(),
            expected,
            "{case}: {}",
            outcome.log_text()
        );
        checked += 1;
    }
    assert_eq!(checked, 39);
}

#[test]
fn topp_fuzzydiff_1_to_4_verdict_parity() {
    // CMakeLists.txt 138-144 (test-data 0cb15f2): WILL_FAIL on 1, 2 and 4.
    let verdict = |case| {
        let (settings, in1, in2) = invocation(case);
        fuzzy_diff(&in1, &in2, &settings)
    };
    let first = verdict("topp_fuzzydiff_1");
    assert_eq!(first.exit, FuzzyDiffExit::ParseError);
    assert!(first.log_text().contains("That's cheating!"));
    assert_eq!(verdict("topp_fuzzydiff_2").exit, FuzzyDiffExit::ParseError);
    assert!(verdict("topp_fuzzydiff_3").passed());
    assert_eq!(
        verdict("topp_fuzzydiff_4").exit,
        FuzzyDiffExit::InputFileNotFound
    );
}

#[test]
fn pinned_fuzzydiff_ini_is_loaded() {
    let settings = upstream();
    assert_eq!(settings.ratio, 1.01);
    assert_eq!(settings.absdiff, 0.01);
    assert_eq!(settings.whitelist, ["<?xml-stylesheet"]);
    assert!(settings.matched_whitelist.is_empty());
    assert_eq!(settings.verbose, 1);
    assert_eq!(settings.tab_width, 8);
    assert_eq!(settings.first_column, 1);
    assert!(!settings.sort);
    // A registration's -whitelist replaces the INI list.
    assert_eq!(settings.clone().with_whitelist(&["id="]).whitelist, ["id="]);
    let defaults = FuzzyDiffSettings::registered_defaults();
    assert_eq!(
        (defaults.ratio, defaults.absdiff, defaults.verbose),
        (1.0, 0.0, 2)
    );
}

#[cfg(feature = "paramxml")]
#[test]
fn ini_reader_agrees_with_the_crate_paramxml_reader() {
    let param = openms::format::paramxml::load(FuzzyDiffSettings::upstream_ini_path()).unwrap();
    let value = |key: &str| param.value(&format!("FuzzyDiff:1:{key}")).unwrap();
    let settings = upstream();
    assert_eq!(value("ratio").to_f64().unwrap(), settings.ratio);
    assert_eq!(value("absdiff").to_f64().unwrap(), settings.absdiff);
    assert_eq!(
        value("whitelist").as_string_list().unwrap(),
        settings.whitelist.as_slice()
    );
    assert_eq!(
        value("matched_whitelist").as_string_list().unwrap(),
        settings.matched_whitelist.as_slice()
    );
    assert_eq!(
        value("verbose").to_i64().unwrap(),
        i64::from(settings.verbose)
    );
    assert_eq!(
        value("tab_width").to_i64().unwrap(),
        i64::from(settings.tab_width)
    );
    assert_eq!(
        value("first_column").to_i64().unwrap(),
        i64::from(settings.first_column)
    );
}

#[test]
fn ini_reader_rejects_malformed_input_and_records_bad_items() {
    assert!(FuzzyDiffSettings::from_ini(b"<PARAMETERS><NODE name=\"FuzzyDiff\">").is_err());
    assert!(
        FuzzyDiffSettings::from_ini(
            b"<NODE name=\"FuzzyDiff\"><ITEM name=\"ratio\" value=\"1</NODE>"
        )
        .is_err()
    );
    assert!(FuzzyDiffSettings::from_ini(b"<ITEM name=\"a\" value=\"&bogus;\"/>").is_err());
    assert!(FuzzyDiffSettings::from_ini(&vec![b' '; fuzzy::MAX_INI_BYTES + 1]).is_err());
    let settings = FuzzyDiffSettings::from_ini(
        b"<?xml version=\"1.0\"?><!-- c --><PARAMETERS><NODE name=\"FuzzyDiff\"><NODE name=\"1\">\
          <ITEM name=\"ratio\" value=\"abc\" type=\"double\"/>\
          <ITEMLIST name=\"whitelist\" type=\"string\"><LISTITEM value=\"a&amp;b &#x41;&#66; &apos;&quot;&gt;\"/></ITEMLIST>\
          </NODE><NODE name=\"Other\"><ITEM name=\"bogus\" value=\"1\"/></NODE></NODE></PARAMETERS>",
    )
    .unwrap();
    assert_eq!(settings.whitelist, ["a&b AB '\">"]);
    assert_eq!(settings.ini_errors.len(), 1, "{:?}", settings.ini_errors);
    let dir = temp_dir("ini-error");
    std::fs::write(dir.join("a.txt"), "x\n").unwrap();
    let outcome = fuzzy_diff(&dir.join("missing"), &dir.join("a.txt"), &settings);
    assert_eq!(
        outcome.exit,
        FuzzyDiffExit::IllegalParameters,
        "INI errors precede file checks"
    );
    std::fs::remove_dir_all(dir).unwrap();
}

// ------------------------------------------------------------------------
// Retained upstream expectations against current C++ output (acceptance 5).
// ------------------------------------------------------------------------

#[test]
fn retained_expectations_match_current_cpp_output_under_registered_whitelists() {
    let cases = [
        (
            "oracle/FileInfo_3.current.txt",
            "retained/FileInfo_3_output.txt",
            "File name",
        ),
        (
            "oracle/FileInfo_17.current.txt",
            "retained/FileInfo_17_output.txt",
            "File name",
        ),
        (
            "oracle/FeatureFinderCentroided_1.current.featureXML",
            "retained/FeatureFinderCentroided_1_1_output.featureXML",
            "id=",
        ),
    ];
    for (current, retained, whitelist) in cases {
        let settings = upstream().with_whitelist(&[whitelist]);
        if let Err(log) = settings.compare_bytes(&read(current), &read(retained)) {
            panic!("{current} vs {retained}:\n{log}");
        }
        assert!(
            upstream()
                .compare_bytes(&read(current), &read(retained))
                .is_err(),
            "{current} needs its whitelist"
        );
    }
    // FileInfo_3 differs from current C++ only by the width of the intensity line;
    // FileInfo_17 is retained with CRLF line ends.
    assert!(
        !read("retained/FileInfo_3_output.txt")
            .windows(12)
            .any(|w| w == b"intensity: 1")
    );
    assert!(
        read("oracle/FileInfo_3.current.txt")
            .windows(12)
            .any(|w| w == b"intensity: 1")
    );
    assert!(
        read("retained/FileInfo_17_output.txt")
            .windows(2)
            .any(|w| w == b"\r\n")
    );
    assert!(!read("oracle/FileInfo_17.current.txt").contains(&b'\r'));
}

#[test]
fn one_digit_changes_beyond_tolerance_fail_and_within_tolerance_pass() {
    let differing = |a: &[u8], b: &[u8]| {
        a.iter().zip(b).filter(|(x, y)| x != y).count() + a.len().abs_diff(b.len())
    };
    let cases = [
        (
            "FileInfo_3_output.beyond.txt",
            "FileInfo_3_output.txt",
            "oracle/FileInfo_3.current.txt",
            "File name",
            false,
        ),
        (
            "FileInfo_3_output.within.txt",
            "FileInfo_3_output.txt",
            "oracle/FileInfo_3.current.txt",
            "File name",
            true,
        ),
        (
            "FileInfo_17_output.beyond.txt",
            "FileInfo_17_output.txt",
            "oracle/FileInfo_17.current.txt",
            "File name",
            false,
        ),
        (
            "FeatureFinderCentroided_1_1_output.beyond.featureXML",
            "FeatureFinderCentroided_1_1_output.featureXML",
            "oracle/FeatureFinderCentroided_1.current.featureXML",
            "id=",
            false,
        ),
        (
            "FeatureFinderCentroided_1_1_output.within.featureXML",
            "FeatureFinderCentroided_1_1_output.featureXML",
            "oracle/FeatureFinderCentroided_1.current.featureXML",
            "id=",
            true,
        ),
        (
            "FeatureFinderCentroided_1_1_output.id_change.featureXML",
            "FeatureFinderCentroided_1_1_output.featureXML",
            "oracle/FeatureFinderCentroided_1.current.featureXML",
            "id=",
            true,
        ),
    ];
    for (mutated, original, current, whitelist, passes) in cases {
        let mutated_bytes = read(&format!("fuzzydiff/{mutated}"));
        assert_eq!(
            differing(&mutated_bytes, &read(&format!("retained/{original}"))),
            1,
            "{mutated}"
        );
        let outcome = upstream()
            .with_whitelist(&[whitelist])
            .compare_bytes(&read(current), &mutated_bytes);
        assert_eq!(outcome.is_ok(), passes, "{mutated}");
    }
}

// ------------------------------------------------------------------------
// Rule details, reports and bounds (source review and Rust-only checks).
// ------------------------------------------------------------------------

#[test]
fn number_rule_keeps_the_source_quirks() {
    let fresh = |a, b, ratio, absdiff| compare_numbers(a, b, ratio, absdiff, 1.0).failure;
    assert_eq!(fresh(1.0, 1.0, 1.0, 0.0), None);
    assert_eq!(fresh(0.0, -0.0, 1.0, 0.0), None);
    assert_eq!(fresh(f64::NAN, f64::NAN, 1.0, 0.0), None);
    assert_eq!(fresh(f64::NAN, 1.0, 1.0, 0.0), Some(message::ONE_NAN));
    assert_eq!(
        fresh(f64::INFINITY, f64::NEG_INFINITY, 1.0, 0.0),
        Some(message::INFINITY_SIGNS)
    );
    assert_eq!(
        fresh(f64::INFINITY, 1.0, 1e300, 1e300),
        Some(message::ONE_INFINITE)
    );
    assert_eq!(fresh(0.0, 0.02, 1.01, 0.01), Some(message::FIRST_ZERO));
    assert_eq!(fresh(0.02, 0.0, 1.01, 0.01), Some(message::SECOND_ZERO));
    assert_eq!(fresh(-0.001, 0.001, 1.01, 0.0), Some(message::SIGNS));
    assert_eq!(fresh(-0.001, 0.001, 1.01, 0.01), None);
    assert_eq!(fresh(1.0, 1.02, 1.01, 0.01), Some(message::RATIO));
    assert_eq!(
        fresh(100.0, 101.0, 1.01, 0.01),
        None,
        "one percent passes although the integers differ"
    );
    // A negative quotient that underflows to -0.0 passes the sign test.
    assert_eq!(fresh(-1e-300, 1e300, 1.01, 0.0), None);
    // A NaN tolerance accepts every ratio; a running maximum hides smaller failures.
    assert_eq!(fresh(1.0, 5.0, f64::NAN, 0.0), None);
    assert_eq!(compare_numbers(1.0, 1.2, 1.01, 0.0, 1.5).failure, None);
}

#[test]
fn reused_comparator_carries_its_relative_maximum() {
    let mut comparator = buffered(1, 1.01, 0.0);
    assert!(!comparator.compare_strings("1", "1.5"));
    assert!(
        comparator.compare_strings("1", "1.2"),
        "source quirk: 1.2 < previous maximum 1.5"
    );
    assert!(!buffered(1, 1.01, 0.0).compare_strings("1", "1.2"));
}

#[test]
fn number_tokens_follow_the_from_chars_contract() {
    let parse = |text: &str| extract_double(text.as_bytes()).map(|(v, n)| (v.to_bits(), n));
    let bits = |value: f64| value.to_bits();
    assert_eq!(parse("1.5e3x"), Some((bits(1500.0), 5)));
    assert_eq!(parse("+-5"), Some((bits(-5.0), 3)));
    assert_eq!(parse("++5"), None);
    assert_eq!(parse("0x10"), Some((bits(0.0), 1)));
    assert_eq!(parse("1e"), Some((bits(1.0), 1)));
    assert_eq!(parse("1e+"), Some((bits(1.0), 1)));
    assert_eq!(parse(".5"), Some((bits(0.5), 2)));
    assert_eq!(parse("5."), Some((bits(5.0), 2)));
    assert_eq!(parse("."), None);
    assert_eq!(parse("-"), None);
    assert_eq!(parse(" 1"), None);
    assert_eq!(parse("1e999"), None, "overflow is rejected");
    assert_eq!(
        parse("1e-400"),
        Some((bits(0.0), 6)),
        "underflow is accepted"
    );
    assert_eq!(parse("infinity"), Some((bits(f64::INFINITY), 8)));
    assert_eq!(parse("-INF"), Some((bits(f64::NEG_INFINITY), 4)));
    assert_eq!(extract_double(b"information").map(|(_, n)| n), Some(3));
    assert_eq!(extract_double(b"nan(x y)z").map(|(_, n)| n), Some(8));
    assert_eq!(extract_double(b"-nan(ab)").map(|(_, n)| n), Some(8));
    assert_eq!(extract_double(b"-nan( x)").map(|(_, n)| n), Some(4));
    assert_eq!(extract_double(b"Nancy").map(|(_, n)| n), Some(3));
    assert_eq!(extract_double(b""), None);
}

#[test]
fn ostream_double_formatting() {
    for (value, text) in [
        (1.0, "1"),
        (1.01, "1.01"),
        (0.11, "0.11"),
        (123_456_789.0, "1.23457e+08"),
        (1e-7, "1e-07"),
        (0.0001, "0.0001"),
        (0.00001, "1e-05"),
        (100_000.0, "100000"),
        (1_000_000.0, "1e+06"),
        (3.49692e6, "3.49692e+06"),
        (-2.5, "-2.5"),
        (-0.0, "-0"),
        (f64::NAN, "nan"),
        (f64::NEG_INFINITY, "-inf"),
    ] {
        assert_eq!(format_g(value), text, "{value:?}");
    }
}

#[test]
fn line_reading_treats_every_carriage_return_as_a_terminator() {
    assert!(buffered(0, 1.0, 0.0).compare_bytes(b"a 1\r\nb 2\r\n", b"a 1\nb 2\n"));
    assert!(buffered(0, 1.0, 0.0).compare_bytes(b"a\rb", b"a\nb"));
    assert!(!buffered(0, 1.0, 0.0).compare_bytes(b"a\rb", b"ab"));
    assert!(buffered(0, 1.0, 0.0).compare_bytes(b"a\n\n \t\x0b\x0c\n", b"a"));
    assert!(
        !buffered(0, 1.0, 0.0).compare_bytes(b"a ", b"a"),
        "trailing whitespace on one side fails"
    );
}

#[test]
fn sort_matches_the_temporary_file_fuzzydiff_writes() {
    assert_eq!(fuzzy::sorted_lines(b"h\nc\na\nb"), b"h\na\nb\nc\n");
    assert_eq!(fuzzy::sorted_lines(b"h\r\nc\r\na\r\n"), b"h\r\na\r\nc\r\n");
    assert_eq!(fuzzy::sorted_lines(b""), b"\n");
    assert_eq!(fuzzy::sorted_lines(b"h\n\nb\n"), b"h\n\nb\n");
}

#[test]
fn read_errors_and_oversized_inputs_fail_instead_of_passing() {
    struct Failing;
    impl std::io::Read for Failing {
        fn read(&mut self, _: &mut [u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("boom"))
        }
    }
    let mut comparator = buffered(1, 1.0, 0.0);
    let mut broken = std::io::BufReader::new(Failing);
    assert!(!comparator.compare_streams(&mut broken, &mut Cursor::new(b"")));
    assert!(String::from_utf8_lossy(comparator.log()).contains("reading input failed"));

    let dir = temp_dir("directory-input");
    let outcome = fuzzy_diff(&dir, &dir.join("missing"), &upstream());
    assert_eq!(outcome.exit, FuzzyDiffExit::InputFileNotFound);
    assert_eq!(
        fuzzy_diff(Path::new(""), &dir, &upstream()).exit,
        FuzzyDiffExit::MissingParameters
    );
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn verbose_3_failure_reports_are_bounded() {
    // One failing letter per short line keeps every report small, so the report
    // count limit is reached before the log byte limit.
    let lhs = "a\n".repeat(fuzzy::MAX_FAILURE_REPORTS + 10);
    let rhs = "b\n".repeat(fuzzy::MAX_FAILURE_REPORTS + 10);
    let mut comparator = buffered(3, 1.0, 0.0);
    assert!(!comparator.compare_strings(&lhs, &rhs));
    let log = String::from_utf8_lossy(comparator.log()).into_owned();
    assert_eq!(
        log.matches("FAILED: 'different letters'").count(),
        fuzzy::MAX_FAILURE_REPORTS
    );
    assert!(log.ends_with("Further failure reports are suppressed.\n"));
    assert!(!comparator.log_truncated());
}

#[test]
#[ignore = "fills the 256 MiB log buffer"]
fn verbose_3_log_bytes_are_bounded() {
    // Every report repeats the whole 20 kB line, so the byte limit comes first.
    let lhs = "a ".repeat(fuzzy::MAX_FAILURE_REPORTS);
    let rhs = "b ".repeat(fuzzy::MAX_FAILURE_REPORTS);
    let mut comparator = buffered(3, 1.0, 0.0);
    assert!(!comparator.compare_strings(&lhs, &rhs));
    assert!(comparator.log_truncated());
    assert_eq!(comparator.log().len(), fuzzy::MAX_LOG_BYTES);
}

// ------------------------------------------------------------------------
// Decoded comparison (decision D6).
// ------------------------------------------------------------------------

use decoded::{DecodedOptions, Tolerance, compare_experiments, compare_feature_maps};
use openms::kernel::{
    ChromatogramPeak, ConvexHull2D, DataArray, Feature, FeatureMap, MSChromatogram, MSExperiment,
    MSSpectrum, Peak1D, Point2D, Precursor,
};
use openms::metadata::{MetaValue, MetaValueData};

fn upstream_tolerance() -> Tolerance {
    Tolerance::from_settings(&FuzzyDiffSettings::upstream().unwrap())
}

fn feature_fixture() -> FeatureMap {
    let mut feature = Feature::new(4389.1, 646.24, 44767.48);
    feature.unique_id = 17;
    feature.charge = 2;
    feature.quality = 0.88;
    feature.metadata.insert(
        "label".into(),
        MetaValue::new(MetaValueData::Integer(0)).unwrap(),
    );
    feature.metadata.insert(
        "score_fit".into(),
        MetaValue::new(MetaValueData::Float(0.78)).unwrap(),
    );
    let mut hull = ConvexHull2D::new();
    hull.set_hull_points(&[
        Point2D::new(4374.19, 646.229),
        Point2D::new(4443.42, 646.229),
        Point2D::new(4443.42, 646.258),
        Point2D::new(4374.19, 646.258),
    ])
    .unwrap();
    feature.convex_hulls = vec![hull.clone(), hull];
    let mut map = FeatureMap::from_features(vec![feature.clone(), feature]);
    map.unique_id = 99;
    map
}

#[test]
fn decoded_feature_maps_report_the_first_mismatch_with_a_path() {
    let options = DecodedOptions::new(upstream_tolerance());
    let expected = feature_fixture();
    compare_feature_maps(&expected, &expected.clone(), &options).unwrap();
    let path_of = |edit: &dyn Fn(&mut FeatureMap)| {
        let mut actual = expected.clone();
        edit(&mut actual);
        compare_feature_maps(&actual, &expected, &options).map_err(|m| m.path)
    };
    assert_eq!(path_of(&|m| m.unique_id = 1), Err("unique_id".into()));
    assert_eq!(
        path_of(&|m| m.features[1].base.rt = 4489.1),
        Err("features[1].rt".into())
    );
    assert_eq!(
        path_of(&|m| m.features[1].base.rt = 4389.2),
        Ok(()),
        "within the ratio"
    );
    assert_eq!(
        path_of(&|m| m.features[0].base.charge = 3),
        Err("features[0].charge".into())
    );
    assert_eq!(
        path_of(&|m| m.features[0].base.intensity = f32::NAN),
        Err("features[0].intensity".into())
    );
    assert_eq!(
        path_of(&|m| {
            m.features.pop();
        }),
        Err("features".into())
    );
    assert_eq!(
        path_of(&|m| {
            m.features[0].base.metadata.remove("label");
        }),
        Err("features[0].metadata".into())
    );
    assert_eq!(
        path_of(&|m| {
            m.features[0].base.metadata.insert(
                "label".into(),
                MetaValue::new(MetaValueData::Float(0.0)).unwrap(),
            );
        }),
        Err("features[0].metadata[\"label\"].value".into())
    );
    assert_eq!(
        path_of(&|m| {
            m.features[0].base.metadata.insert(
                "score_fit".into(),
                MetaValue::new(MetaValueData::Float(0.9)).unwrap(),
            );
        }),
        Err("features[0].metadata[\"score_fit\"].value".into())
    );
    assert_eq!(
        path_of(&|m| {
            let mut points = m.features[1].convex_hulls[1].hull_points();
            points[3].mz = 700.0;
            m.features[1].convex_hulls[1]
                .set_hull_points(&points)
                .unwrap();
        }),
        Err("features[1].convex_hulls[1].points[3].mz".into())
    );
    assert_eq!(
        path_of(&|m| {
            let mut points = m.features[0].convex_hulls[0].hull_points();
            points.swap(0, 1);
            m.features[0].convex_hulls[0]
                .set_hull_points(&points)
                .unwrap();
        }),
        Err("features[0].convex_hulls[0].points[0].rt".into()),
        "hull point order is part of the contract"
    );
    assert_eq!(
        path_of(&|m| {
            m.features[0].convex_hulls.pop();
        }),
        Err("features[0].convex_hulls".into())
    );
    assert_eq!(
        path_of(&|m| m.features[0].subordinates.push(Feature::new(1.0, 2.0, 3.0))),
        Err("features[0].subordinates".into())
    );
    let mut with_subordinate = expected.clone();
    with_subordinate.features[0]
        .subordinates
        .push(Feature::new(1.0, 2.0, 3.0));
    let mut changed = with_subordinate.clone();
    changed.features[0].subordinates[0].base.mz = 3.0;
    assert_eq!(
        compare_feature_maps(&changed, &with_subordinate, &options)
            .unwrap_err()
            .path,
        "features[0].subordinates[0].mz"
    );
    assert_eq!(
        path_of(&|m| m.features[0].base.unique_id = 5),
        Err("features[0].unique_id".into())
    );
    let ignoring = options.ignoring_unique_ids();
    let mut renumbered = expected.clone();
    renumbered.unique_id = 5;
    renumbered.features[0].base.unique_id = 6;
    renumbered.loaded_file_path = "elsewhere.featureXML".into();
    compare_feature_maps(&renumbered, &expected, &ignoring).unwrap();
}

#[test]
fn decoded_experiments_report_the_first_mismatch_with_a_path() {
    let options = DecodedOptions::new(upstream_tolerance());
    let mut spectrum = MSSpectrum {
        rt: 12.5,
        ms_level: 2,
        native_id: "scan=100".into(),
        peaks: vec![Peak1D::new(100.0, 10.0), Peak1D::new(200.0, 20.0)],
        ..MSSpectrum::default()
    };
    spectrum.precursors.push(Precursor::new(500.25, 2));
    spectrum
        .float_data_arrays
        .push(DataArray::new("Ion Mobility", vec![1.0, 2.0]));
    let chromatogram = MSChromatogram {
        native_id: "TIC".into(),
        peaks: vec![
            ChromatogramPeak::new(1.0, 5.0),
            ChromatogramPeak::new(2.0, 6.0),
            ChromatogramPeak::new(3.0, 7.0),
        ],
        ..MSChromatogram::default()
    };
    let expected = MSExperiment {
        spectra: vec![MSSpectrum::default(), spectrum],
        chromatograms: vec![chromatogram],
        ..MSExperiment::default()
    };
    compare_experiments(&expected, &expected.clone(), &options).unwrap();
    let path_of = |edit: &dyn Fn(&mut MSExperiment)| {
        let mut actual = expected.clone();
        edit(&mut actual);
        compare_experiments(&actual, &expected, &options).map_err(|m| m.path)
    };
    assert_eq!(
        path_of(&|e| e.spectra[1].peaks[1].intensity = 30.0),
        Err("spectra[1].peaks[1].intensity".into())
    );
    assert_eq!(path_of(&|e| e.spectra[1].peaks[1].intensity = 20.1), Ok(()));
    assert_eq!(
        path_of(&|e| {
            e.spectra[1].peaks.pop();
        }),
        Err("spectra[1].peaks".into())
    );
    assert_eq!(
        path_of(&|e| e.spectra[1].ms_level = 1),
        Err("spectra[1].ms_level".into())
    );
    assert_eq!(
        path_of(&|e| e.spectra[1].native_id = "scan=101".into()),
        Err("spectra[1].native_id".into())
    );
    assert_eq!(
        path_of(&|e| e.spectra[1].precursors[0].charge = 3),
        Err("spectra[1].precursors[0].charge".into())
    );
    assert_eq!(
        path_of(&|e| e.spectra[1].precursors[0].isolation_target_mz = Some(500.0)),
        Err("spectra[1].precursors[0].isolation_target_mz".into())
    );
    assert_eq!(
        path_of(&|e| e.spectra[1].float_data_arrays[0].data[1] = 3.0),
        Err("spectra[1].float_data_arrays[0].data[1]".into())
    );
    assert_eq!(
        path_of(&|e| e.spectra[1].float_data_arrays[0].name = "other".into()),
        Err("spectra[1].float_data_arrays[0].name".into())
    );
    assert_eq!(
        path_of(&|e| e.chromatograms[0].peaks[2].rt = 4.0),
        Err("chromatograms[0].peaks[2].rt".into())
    );
    assert_eq!(
        path_of(&|e| e.settings.comment = "changed".into()),
        Err("settings".into())
    );
    assert_eq!(
        path_of(&|e| e.settings.document.loaded_file_path = "x.mzML".into()),
        Ok(())
    );
    let ignoring = options.ignoring_unique_ids();
    let mut renamed = expected.clone();
    renamed.spectra[1].native_id = "index=1".into();
    renamed.chromatograms[0].native_id = "chromatogram=0".into();
    compare_experiments(&renamed, &expected, &ignoring).unwrap();
    let mismatch = compare_experiments(&renamed, &expected, &options).unwrap_err();
    assert_eq!(
        mismatch.to_string(),
        "spectra[1].native_id: \"index=1\" vs \"scan=100\""
    );
}

#[cfg(feature = "featurexml")]
#[test]
fn decoded_ffc_1_current_cpp_output_matches_the_retained_expectation() {
    let load = |relative: &str| openms::format::featurexml::load(data(relative)).unwrap();
    let current = load("oracle/FeatureFinderCentroided_1.current.featureXML");
    let retained = load("retained/FeatureFinderCentroided_1_1_output.featureXML");
    assert_eq!(current.features.len(), 8);
    let options = DecodedOptions::new(upstream_tolerance()).ignoring_unique_ids();
    compare_feature_maps(&current, &retained, &options).unwrap();
    // Without the id exclusion the generated map ids differ first.
    let strict = DecodedOptions::new(upstream_tolerance());
    assert_eq!(
        compare_feature_maps(&current, &retained, &strict)
            .unwrap_err()
            .path,
        "unique_id"
    );
    // One-digit mutations of the retained file, as run through C++ FuzzyDiff.
    let beyond = load("fuzzydiff/FeatureFinderCentroided_1_1_output.beyond.featureXML");
    let mismatch = compare_feature_maps(&current, &beyond, &options).unwrap_err();
    assert_eq!(mismatch.path, "features[1].rt", "{mismatch}");
    let within = load("fuzzydiff/FeatureFinderCentroided_1_1_output.within.featureXML");
    compare_feature_maps(&current, &within, &options).unwrap();
    let id_change = load("fuzzydiff/FeatureFinderCentroided_1_1_output.id_change.featureXML");
    compare_feature_maps(&retained, &id_change, &options).unwrap();
    assert_eq!(
        compare_feature_maps(&retained, &id_change, &strict)
            .unwrap_err()
            .path,
        "features[0].unique_id"
    );
    // The retained file against itself is exact, and a tighter tolerance than the
    // upstream one exposes the ~1e-10 drift between retained and current C++.
    compare_feature_maps(
        &retained,
        &retained.clone(),
        &DecodedOptions::new(Tolerance::exact()),
    )
    .unwrap();
    assert!(
        compare_feature_maps(
            &current,
            &retained,
            &DecodedOptions::new(Tolerance::exact()).ignoring_unique_ids()
        )
        .is_err()
    );
}

#[cfg(feature = "mzml")]
#[test]
fn decoded_mzml_experiment_matches_itself_and_reports_changes() {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/baseline_filter_tool_input.mzML");
    let expected = openms::format::mzml::load(&path).unwrap();
    let options = DecodedOptions::new(upstream_tolerance());
    let mut reloaded = openms::format::mzml::load(&path).unwrap();
    reloaded.settings.document.loaded_file_path = "other/place.mzML".into();
    compare_experiments(&reloaded, &expected, &options).unwrap();
    let index = expected
        .spectra
        .iter()
        .position(|s| !s.peaks.is_empty())
        .expect("a spectrum with peaks");
    let mut changed = expected.clone();
    changed.spectra[index].peaks[0].intensity =
        changed.spectra[index].peaks[0].intensity * 1.5 + 1.0;
    assert_eq!(
        compare_experiments(&changed, &expected, &options)
            .unwrap_err()
            .path,
        format!("spectra[{index}].peaks[0].intensity")
    );
}
