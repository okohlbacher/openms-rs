// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! A5-FILEINFO-TOOL: the FileInfo TOPP tool (`OpenMS4-topp/src/FileInfo.cpp`,
//! topp `174b576`), its executable and the preview workflows.
//!
//! Evidence, in order of strength (see `tests/data/topp_file_info_provenance.json`
//! and `docs/TOPP_FILE_INFO_SUPPORT.md`):
//! - tier 1, retained upstream outputs: TOPP_FileInfo_1, _2, _3 and _9, and
//!   since A8 _4, _5 and _6, which also match the Release build's `-out` byte
//!   for byte (`../oracle/a8-fileinfo`)
//!   (test-data `0cb15f2`, `topp/CMakeLists.txt:881-904`) run through the tool
//!   and compared with FuzzyDiff (`FuzzyDiff.ini`, ratio 1.01, absdiff 0.01) and
//!   the registered whitelist `File name`;
//! - tier 1, executed differential: the product-SDK FileInfo (Debug, core
//!   `4fdec46`) from the C1 oracle (`../oracle/topp-early-bundle`, run1) and
//!   from this package's oracle (`../oracle/topp-file-info-tool`, run1, both
//!   runs reproduced): exit codes, diagnostics, usage text, `-write_ini`, and
//!   reports compared byte for byte with only the `File name` lines normalised;
//! - tier 4, the explicit refusal of every branch and flag the preview does not
//!   run, and native checks of the output routing.
//!
//! Every case runs in its own temporary directory, because the tests run in
//! parallel.

#![cfg(all(feature = "mzml", feature = "paramxml", feature = "featurexml"))]

#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;
#[path = "support/release_runs.rs"]
mod release_runs;
#[path = "support/took_line.rs"]
mod took_line;

use openms::cli::tools::FileInfo;
use openms::cli::{ExitCode, Tool, run_with, verbose_version};
use openms::format::file_info::model::Options;
use openms::format::file_info::report::FileInfo as FileInfoLibrary;
use openms::system::file::TempDir;
use std::fs;
use std::path::{Path, PathBuf};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// The C++ usage version line: the product SDK names its build revision.
const CPP_VERSION: &str = "1.0.0 (OpenMS core 4.0.0, revision 4fdec46)";

/// The console-width probe of the C++ usage text, absent without a terminal
/// in this port (`docs/TOPP_CLI_SUPPORT.md`, *Usage text*).
const STTY_LINE: &str = "stty: stdin isn't a terminal\n";

fn data(relative: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(relative)
        .to_string_lossy()
        .into_owned()
}

/// A4's library fixtures, read in place.
fn library(relative: &str) -> String {
    data(&format!("file_info/{relative}"))
}

/// This package's fixtures.
fn tool(relative: &str) -> String {
    data(&format!("topp_file_info/{relative}"))
}

fn read(path: impl AsRef<Path>) -> String {
    let path = path.as_ref();
    fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

struct Outcome {
    code: ExitCode,
    /// The output stream without the closing `FileInfo took …` line.
    out: String,
    err: String,
    /// The closing line: present when the source's `main_` returns, absent
    /// when it throws.
    took: Option<String>,
}

/// Run the tool as `FileInfo args...` through `run_with`. The closing
/// `FileInfo took …` line of a completed run is checked and kept apart from
/// the output stream (`support/took_line.rs`).
fn run(args: &[&str]) -> Outcome {
    let arguments: Vec<String> = std::iter::once(FileInfo::NAME.to_owned())
        .chain(args.iter().map(|a| (*a).to_owned()))
        .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<FileInfo>(&arguments, &mut out, &mut err);
    let (out, took) = took_line::split_took_line(
        "FileInfo",
        &String::from_utf8(out).expect("UTF-8 output stream"),
    );
    Outcome {
        code,
        out,
        err: String::from_utf8(err).expect("UTF-8 error stream"),
        took,
    }
}

/// A fresh, uniquely named directory per case, removed when the case ends.
struct Workdir(TempDir);
impl Workdir {
    fn new() -> Self {
        Self(TempDir::new_in(std::env::temp_dir(), false).unwrap())
    }
    fn file(&self, name: &str) -> String {
        self.0.path().join(name).to_string_lossy().into_owned()
    }
}

/// Replace the value of every `File name: ` text line and `general: file name`
/// TSV line, which embed the path the report ran on.
fn normalise_file_name(report: &str) -> String {
    report
        .split_inclusive('\n')
        .map(|line| {
            if line.starts_with("File name: ") {
                "File name: <input>\n".to_owned()
            } else if line.starts_with("general: file name\t") {
                "general: file name\t<input>\n".to_owned()
            } else {
                line.to_owned()
            }
        })
        .collect()
}

/// The first differing line, for a readable failure.
fn first_difference(actual: &str, expected: &str) -> String {
    let mut expected_lines = expected.split_inclusive('\n');
    for (index, line) in actual.split_inclusive('\n').enumerate() {
        match expected_lines.next() {
            Some(other) if other == line => {}
            other => return format!("line {}: actual {line:?}, expected {other:?}", index + 1),
        }
    }
    match expected_lines.next() {
        Some(line) => format!("actual ends early; expected next {line:?}"),
        None => "no line differs".to_owned(),
    }
}

/// A report equal to the C++ file byte for byte, apart from the file name.
fn assert_report(actual: &str, expected_file: &str) {
    assert_report_text(actual, &read(expected_file), expected_file);
}

/// A report equal to the C++ text byte for byte, apart from the file name.
fn assert_report_text(actual: &str, expected: &str, label: &str) {
    let actual = normalise_file_name(actual);
    let expected = normalise_file_name(expected);
    assert!(
        actual == expected,
        "{label}: {}",
        first_difference(&actual, &expected)
    );
}

/// The C++ stdout of a report run without `-out`: the report, then the
/// `FileInfo took ...` timing line, which is taken off here as the port's is
/// taken off its own output.
fn cpp_stdout_report(file: &str) -> String {
    let stdout = read(file);
    stdout
        .strip_suffix('\n')
        .and_then(|text| text.rsplit_once('\n'))
        .map(|(body, last)| {
            assert!(last.starts_with("FileInfo took "), "{last}");
            format!("{body}\n")
        })
        .unwrap_or_else(|| panic!("{file}: no timing line"))
}

/// A C++ stderr with the usage text: the `stty` probe line removed and the
/// build revision replaced by this port's.
fn cpp_usage_stream(file: &str) -> String {
    let text = read(file);
    assert_eq!(text.matches(STTY_LINE).count(), 1, "{file}");
    text.replacen(STTY_LINE, "", 1)
        .replace(CPP_VERSION, &verbose_version::<FileInfo>())
}

fn fuzzy_diff_against_retained(produced: &str, retained: &str) {
    let settings = fuzzy::FuzzyDiffSettings::load_ini(Path::new(&data(
        "fuzzy_string_comparator/FuzzyDiff.ini",
    )))
    .unwrap()
    .with_whitelist(&["File name"]);
    assert_eq!(settings.ratio, 1.01);
    assert_eq!(settings.absdiff, 0.01);
    let expected = fs::read(retained).unwrap();
    if let Err(log) = settings.compare_bytes(read(produced).as_bytes(), &expected) {
        panic!("{retained}: FuzzyDiff failed:\n{log}");
    }
}

fn assert_code(outcome: &Outcome, code: ExitCode) {
    assert_eq!(
        outcome.code, code,
        "out:\n{}\nerr:\n{}",
        outcome.out, outcome.err
    );
}

// ---------------------------------------------------------------------------
// Registration and usage (tier 1: ../oracle/topp-file-info-tool help,
// helphelp and write_ini; C1 FileInfo_no_args)
// ---------------------------------------------------------------------------

/// `--help` and `--helphelp` against the C++ usage text, byte for byte on the
/// error stream apart from the `stty` probe and the revision, with nothing on
/// the output stream. The text carries every registered option, its verbatim
/// description, valid formats and strings, and the advanced `-out_tsv`.
#[test]
fn usage_text_matches_the_cpp_text() {
    for (flag, file) in [
        ("--help", "expected/help.stderr.txt"),
        ("--helphelp", "expected/helphelp.stderr.txt"),
    ] {
        let outcome = run(&[flag]);
        assert_code(&outcome, ExitCode::ExecutionOk);
        assert!(outcome.out.is_empty(), "{}", outcome.out);
        assert_eq!(outcome.err, cpp_usage_stream(&tool(file)), "{flag}");
    }
}

/// C1 `FileInfo_no_args`: usage, `No options given. Aborting!`, exit 6.
#[test]
fn a_bare_invocation_prints_usage_and_exits_6() {
    let arguments = vec![FileInfo::NAME.to_owned()];
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<FileInfo>(&arguments, &mut out, &mut err);
    assert_eq!(code, ExitCode::IllegalParameters);
    assert!(out.is_empty());
    assert_eq!(
        String::from_utf8(err).unwrap(),
        cpp_usage_stream(&tool("expected/FileInfo_no_args.stderr.txt"))
    );
}

/// Oracle `write_ini`: the file equals the C++ `-test -write_ini` file line by
/// line with exact numbers, `version` lines skipped, and as a decoded
/// parameter tree with descriptions, tags, restrictions and formats.
#[test]
fn write_ini_matches_the_cpp_file() {
    let dir = Workdir::new();
    let written = dir.file("FileInfo.tmp.ini");
    let outcome = run(&["-test", "-write_ini", &written]);
    assert_code(&outcome, ExitCode::ExecutionOk);
    let expected = tool("expected/write_ini.ini");

    let mut comparator = fuzzy::FuzzyStringComparator::new();
    comparator.set_acceptable_relative(1.0);
    comparator.set_acceptable_absolute(0.0);
    comparator.set_whitelist(vec!["version".to_owned()]);
    comparator.set_log_destination(fuzzy::LogDestination::Buffer);
    assert!(
        comparator.compare_files(Path::new(&written), Path::new(&expected)),
        "{}",
        String::from_utf8_lossy(comparator.log())
    );
    let produced = openms::format::paramxml::load(&written).unwrap();
    let oracle = openms::format::paramxml::load(&expected).unwrap();
    assert_eq!(produced, oracle);
}

/// Oracle `ini_with_flags`: the INI the C++ tool writes, with `test`,
/// `no_progress`, `m`, `p` and `s` set, is accepted as it is, and the run
/// equals the C++ run on the same INI. A benchmark that hands both tools one
/// C++-written INI depends on this.
#[test]
fn an_ini_written_by_the_cpp_tool_is_accepted() {
    let dir = Workdir::new();
    let (text, tsv) = (dir.file("ini.tmp.txt"), dir.file("ini.tmp.tsv"));
    let outcome = run(&[
        "-ini",
        &tool("inputs/FileInfo_flags.ini"),
        "-in",
        &library("inputs/FileInfo_3_input.featureXML"),
        "-out",
        &text,
        "-out_tsv",
        &tsv,
    ]);
    assert_code(&outcome, ExitCode::ExecutionOk);
    assert!(outcome.err.is_empty(), "{}", outcome.err);
    assert!(outcome.out.is_empty(), "{}", outcome.out);
    assert_report(&read(&text), &tool("expected/ini_with_flags.txt"));
    assert_report(&read(&tsv), &tool("expected/ini_with_flags.tsv"));
}

// ---------------------------------------------------------------------------
// Upstream registrations (tier 1: retained outputs through FuzzyDiff, and the
// C1 oracle outputs byte for byte)
// ---------------------------------------------------------------------------

/// Run a registration into `<dir>/<name>.tmp.txt` and return that path.
fn registration(dir: &Workdir, name: &str, args: &[&str]) -> String {
    let out = dir.file(&format!("{name}.tmp.txt"));
    let mut full: Vec<&str> = args.to_vec();
    full.extend(["-out", &out]);
    let outcome = run(&full);
    assert_code(&outcome, ExitCode::ExecutionOk);
    assert!(outcome.out.is_empty(), "{}", outcome.out);
    out
}

/// TOPP_FileInfo_1 (`CMakeLists.txt:881-883`).
#[test]
fn topp_file_info_1() {
    let dir = Workdir::new();
    let input = library("inputs/FileInfo_1_input.dta");
    let out = registration(
        &dir,
        "FileInfo_1",
        &["-test", "-in", &input, "-in_type", "dta", "-no_progress"],
    );
    fuzzy_diff_against_retained(&out, &library("retained/FileInfo_1_output.txt"));
    // C1 TOPP_FileInfo_1 wrote the same bytes as C1 FileInfo_1_tsv's -out.
    assert_report(&read(&out), &library("expected/FileInfo_1_tsv.txt"));
}

/// TOPP_FileInfo_7 (`CMakeLists.txt:899-901`), reproduced by A7 with the
/// registration's own `-s -m -p`. The retained output is the upstream
/// `FileInfo_7_output.txt`, compared through FuzzyDiff as the registration
/// does; `tests/file_info_a7.rs` compares the same report with the Release C++
/// output byte for byte.
#[cfg(feature = "consensusxml")]
#[test]
fn topp_file_info_7() {
    let dir = Workdir::new();
    let input = tool("inputs/FileInfo_7_input.consensusXML");
    let out = registration(
        &dir,
        "FileInfo_7",
        &["-test", "-in", &input, "-s", "-m", "-p", "-no_progress"],
    );
    fuzzy_diff_against_retained(&out, &library("retained/FileInfo_7_output.txt"));
}

/// TOPP_FileInfo_10 (`CMakeLists.txt:905-907`), reproduced by A7.
#[cfg(feature = "idxml")]
#[test]
fn topp_file_info_10() {
    let dir = Workdir::new();
    let input = tool("inputs/FileInfo_10_input.idXML");
    let out = registration(
        &dir,
        "FileInfo_10",
        &["-test", "-in", &input, "-no_progress"],
    );
    fuzzy_diff_against_retained(&out, &library("retained/FileInfo_10_output.txt"));
}

/// TOPP_FileInfo_13 (`CMakeLists.txt:912`), reproduced by A7. The registration
/// has no retained output and no comparison: it exists only to show that an
/// empty consensusXML does not crash the tool, so the exit code is the whole
/// assertion.
#[cfg(feature = "consensusxml")]
#[test]
fn topp_file_info_13() {
    let input = tool("inputs/FileInfo_13_input.consensusXML");
    let outcome = run(&["-test", "-in", &input, "-no_progress"]);
    assert_code(&outcome, ExitCode::ExecutionOk);
    assert!(
        outcome
            .out
            .contains("No consensus features found, map is empty!\n"),
        "{}",
        outcome.out
    );
}

/// TOPP_FileInfo_17, _18 and _20 (`CMakeLists.txt:922-924`, `:925-927`,
/// `:931-933`), reproduced by A7.
#[test]
fn topp_file_info_17_18_and_20() {
    for name in ["FileInfo_17", "FileInfo_18", "FileInfo_20"] {
        let dir = Workdir::new();
        let input = tool(&format!("inputs/{name}_input.fasta"));
        let out = registration(&dir, name, &["-test", "-in", &input, "-no_progress"]);
        fuzzy_diff_against_retained(&out, &library(&format!("retained/{name}_output.txt")));
    }
}

/// TOPP_FileInfo_4, _5 and _6 (`CMakeLists.txt:890-892`, `:893-895`,
/// `:896-898`), reproduced by A8 with the registrations' own flags: mzXML with
/// `-m`, a `.mzDat` file forced to mzData with `-m -s`, and mzData with
/// `-d -s`. The retained outputs are compared through FuzzyDiff as the
/// registrations do, and each `-out` byte for byte with the Release build's
/// (`../oracle/a8-fileinfo`, cases `x4`, `d5` and `d6`, which
/// `tests/file_info_a8.rs` runs through the library as well).
#[test]
fn topp_file_info_4_5_and_6() {
    let fi4 = tool("inputs/FileInfo_4_input.mzXML");
    let fi5 = data("mzml_mobility/FileInfo_5_input.mzDat");
    let fi6 = tool("inputs/FileInfo_6_input.mzData");
    let cases: [(&str, &str, Vec<&str>); 3] = [
        ("FileInfo_4", "x4", vec!["-in", &fi4, "-m"]),
        (
            "FileInfo_5",
            "d5",
            vec!["-in", &fi5, "-in_type", "mzData", "-m", "-s"],
        ),
        ("FileInfo_6", "d6", vec!["-in", &fi6, "-d", "-s"]),
    ];
    for (name, case, flags) in cases {
        let dir = Workdir::new();
        let mut args = vec!["-test"];
        args.extend(flags);
        args.push("-no_progress");
        let out = registration(&dir, name, &args);
        fuzzy_diff_against_retained(&out, &library(&format!("retained/{name}_output.txt")));
        assert_report(
            &read(&out),
            &data(&format!("file_info_a8/expected/{case}.txt")),
        );
    }
}

/// Oracle `x4_bare`: TOPP_FileInfo_4 without `-out`, so the report goes to the
/// output stream; byte for byte with the Release build's standard output
/// apart from its `FileInfo took` line.
#[test]
fn topp_file_info_4_on_the_output_stream() {
    let fi4 = tool("inputs/FileInfo_4_input.mzXML");
    let outcome = run(&["-test", "-in", &fi4, "-no_progress", "-m"]);
    assert_code(&outcome, ExitCode::ExecutionOk);
    assert!(outcome.err.is_empty(), "{}", outcome.err);
    let expected = data("file_info_a8/expected/x4_bare.stdout.txt");
    assert_report_text(&outcome.out, &cpp_stdout_report(&expected), &expected);
}

/// Oracle `x4_i`: `-i` on an mzXML file is refused by the tool before the
/// class runs (`OpenMS4-topp/src/FileInfo.cpp:118-121`), exit 6, with the same
/// error line and usage text as on DTA; the Release build's stderr starts with
/// that line too. The class itself runs the check (`tests/file_info_a8.rs`,
/// driver case `lib_x4_i`).
#[test]
fn the_index_check_on_mzxml_exits_6_with_usage() {
    let input = tool("inputs/FileInfo_4_input.mzXML");
    let outcome = run(&["-test", "-in", &input, "-i", "-no_progress"]);
    assert_code(&outcome, ExitCode::IllegalParameters);
    assert!(outcome.out.is_empty(), "{}", outcome.out);
    assert!(outcome.took.is_some(), "main_ returned 6");
    assert!(
        outcome
            .err
            .starts_with("Error: Can only validate indices for mzML files\n"),
        "{}",
        outcome.err
    );
    assert_eq!(
        outcome.err,
        cpp_usage_stream(&tool("expected/FileInfo_index_on_dta.stderr.txt"))
    );
}

/// Oracle `m_test`: the tool's `-in` lists the formats it accepts, and `ms2`
/// is not among them, so an MS2 file is refused before the library runs, in
/// the Release build as here: exit 6 and the same message. The library's MS2
/// branch is reached only through the class (`tests/file_info_a8.rs`).
#[test]
fn an_ms2_input_is_refused_by_the_input_format_check() {
    let input = data("text_peak_lists/MS2File_test_spectra.ms2");
    let outcome = run(&["-test", "-in", &input, "-no_progress"]);
    assert_code(&outcome, ExitCode::IllegalParameters);
    assert!(outcome.out.is_empty(), "{}", outcome.out);
    assert_eq!(
        outcome.err,
        format!(
            "Invalid parameter: Input file '{input}' has invalid format 'ms2'. Valid formats are: \
             'mzData','mzXML','mzML','sqMass','dta','dta2d','mgf','featureXML','consensusXML',\
             'idXML','pepXML','mzTab','fid','mzid','trafoXML','fasta','pqp'.\n"
        )
    );
}

/// A7 implemented the consensusXML, idXML, mzIdentML and FASTA branches, so the
/// six rows the not-ported table used to hold for them are gone. Pinned here is
/// only that the tool no longer answers those inputs with a not-ported refusal,
/// so the table cannot quietly regain them.
#[test]
fn the_branches_a7_implemented_are_no_longer_refused() {
    // Each branch is behind the feature that carries its format, and this file
    // is built in slices that lack some of them (CI's minimum-rust line builds
    // `mzml paramxml featurexml` only). A branch whose feature is absent is
    // refused with "this build lacks the <feature> feature", which is the
    // build's own answer and not the not-ported refusal this test guards.
    let mut cases: Vec<(String, &[&str])> = Vec::new();
    #[cfg(feature = "consensusxml")]
    {
        cases.push((tool("inputs/FileInfo_7_input.consensusXML"), &["-s"][..]));
        cases.push((tool("inputs/FileInfo_13_input.consensusXML"), &[][..]));
    }
    #[cfg(feature = "idxml")]
    cases.push((tool("inputs/FileInfo_10_input.idXML"), &[][..]));
    cases.push((tool("inputs/FileInfo_17_input.fasta"), &[][..]));
    cases.push((tool("inputs/FileInfo_18_input.fasta"), &[][..]));
    cases.push((tool("inputs/FileInfo_20_input.fasta"), &[][..]));
    for (input, flags) in cases {
        let dir = Workdir::new();
        let out = dir.file("ported.tmp.txt");
        let mut full: Vec<&str> = vec!["-test", "-no_progress", "-in", &input, "-out", &out];
        full.extend(flags.iter().copied());
        let outcome = run(&full);
        assert_code(&outcome, ExitCode::ExecutionOk);
        assert!(
            !outcome.err.contains("is not ported"),
            "{input}: {}",
            outcome.err
        );
        assert!(read(&out).contains("-- General information --"), "{input}");
    }
}

/// TOPP_FileInfo_2 (`CMakeLists.txt:884-886`).
#[test]
fn topp_file_info_2() {
    let dir = Workdir::new();
    let input = library("inputs/FileInfo_2_input.dta2d");
    let out = registration(
        &dir,
        "FileInfo_2",
        &["-test", "-in", &input, "-no_progress"],
    );
    fuzzy_diff_against_retained(&out, &library("retained/FileInfo_2_output.txt"));
    assert_report(&read(&out), &library("expected/FileInfo_2_tsv.txt"));
}

/// TOPP_FileInfo_3 (`CMakeLists.txt:887-889`). The retained file writes six
/// spaces after `intensity:`; the source at the pin and the port write one,
/// which FuzzyDiff's whitespace rule accepts.
#[test]
fn topp_file_info_3() {
    let dir = Workdir::new();
    let input = library("inputs/FileInfo_3_input.featureXML");
    let out = registration(
        &dir,
        "FileInfo_3",
        &["-test", "-in", &input, "-m", "-s", "-p", "-no_progress"],
    );
    assert!(read(&out).contains("  intensity: 1376.00 .. 6712.00\n"));
    fuzzy_diff_against_retained(&out, &library("retained/FileInfo_3_output.txt"));
    assert_report(&read(&out), &library("expected/FileInfo_3_tsv.txt"));
}

/// TOPP_FileInfo_9 (`CMakeLists.txt:902-904`) on its registered input: the
/// strict native mzML reader still refuses it, for three gaps outside this
/// package (`docs/FILE_INFO_SUPPORT.md`, *Known reader gaps*), where C++ exits
/// 0. The dangling-reference option the tool enables does not cover them. When
/// a reader change makes the input load, this test fails; then compare it as
/// `topp_file_info_9_on_the_derived_input` does.
#[test]
fn topp_file_info_9_registered_input_is_refused_by_the_mzml_reader() {
    let dir = Workdir::new();
    let out = dir.file("FileInfo_9.tmp.txt");
    let input = library("inputs/FileInfo_9_input.mzML");
    let outcome = run(&[
        "-test",
        "-in",
        &input,
        "-m",
        "-p",
        "-s",
        "-no_progress",
        "-out",
        &out,
    ]);
    assert_code(&outcome, ExitCode::IncompatibleInputData);
    assert_eq!(
        outcome.err,
        "Error: unsupported: duplicate userParam name name\n"
    );
    assert_eq!(read(&out), "");
}

/// TOPP_FileInfo_9 on the derived input whose C++ report equals the
/// registered input's apart from the file name (A4 oracle
/// `fileinfo_9_strict_reader_mps`).
#[test]
fn topp_file_info_9_on_the_derived_input() {
    let dir = Workdir::new();
    let input = library("inputs/FileInfo_9_strict_reader.mzML");
    let out = registration(
        &dir,
        "FileInfo_9",
        &["-test", "-in", &input, "-m", "-p", "-s", "-no_progress"],
    );
    fuzzy_diff_against_retained(&out, &library("retained/FileInfo_9_output.txt"));
    assert_report(&read(&out), &library("expected/FileInfo_9_tsv.txt"));
}

// ---------------------------------------------------------------------------
// C1 oracle regressions (tier 1)
// ---------------------------------------------------------------------------

/// Run with `-out` and `-out_tsv` in a fresh directory and compare both files
/// with the C1 pair.
fn check_text_and_tsv(args: &[&str], expected_text: &str, expected_tsv: &str) {
    let dir = Workdir::new();
    let (text, tsv) = (dir.file("report.tmp.txt"), dir.file("report.tmp.tsv"));
    let mut full: Vec<&str> = args.to_vec();
    full.extend(["-out", &text, "-out_tsv", &tsv]);
    let outcome = run(&full);
    assert_code(&outcome, ExitCode::ExecutionOk);
    assert!(outcome.out.is_empty(), "{}", outcome.out);
    assert_report(&read(&text), expected_text);
    assert_report(&read(&tsv), expected_tsv);
}

/// C1 `FileInfo_1_tsv`, `_2_tsv`, `_3_tsv` and `_9_tsv` (`_9` on the derived
/// input): text and TSV exact apart from the file-name lines.
#[test]
fn c1_text_and_tsv_of_the_registrations() {
    let fi1 = library("inputs/FileInfo_1_input.dta");
    check_text_and_tsv(
        &["-test", "-in", &fi1, "-in_type", "dta", "-no_progress"],
        &library("expected/FileInfo_1_tsv.txt"),
        &library("expected/FileInfo_1_tsv.tsv"),
    );
    let fi2 = library("inputs/FileInfo_2_input.dta2d");
    check_text_and_tsv(
        &["-test", "-in", &fi2, "-no_progress"],
        &library("expected/FileInfo_2_tsv.txt"),
        &library("expected/FileInfo_2_tsv.tsv"),
    );
    let fi3 = library("inputs/FileInfo_3_input.featureXML");
    check_text_and_tsv(
        &["-test", "-in", &fi3, "-m", "-s", "-p", "-no_progress"],
        &library("expected/FileInfo_3_tsv.txt"),
        &library("expected/FileInfo_3_tsv.tsv"),
    );
    let fi9 = library("inputs/FileInfo_9_strict_reader.mzML");
    check_text_and_tsv(
        &["-test", "-in", &fi9, "-m", "-p", "-s", "-no_progress"],
        &library("expected/FileInfo_9_tsv.txt"),
        &library("expected/FileInfo_9_tsv.tsv"),
    );
}

/// C1 `FileInfo_empty_featureXML`.
#[test]
fn c1_empty_featurexml() {
    let input = library("inputs/empty.featureXML");
    check_text_and_tsv(
        &["-test", "-in", &input, "-m", "-p", "-s", "-no_progress"],
        &library("expected/FileInfo_empty_featureXML.txt"),
        &library("expected/FileInfo_empty_featureXML.tsv"),
    );
}

/// C1 `FileInfo_empty_mzML`, on the original input with its dangling
/// `defaultDataProcessingRef`. The tool enables the source-compatible reading
/// (decision D10); the library default stays strict and refuses the same file.
#[test]
fn c1_empty_mzml_with_a_dangling_reference() {
    let input = library("inputs/empty.mzML");
    check_text_and_tsv(
        &["-test", "-in", &input, "-m", "-p", "-s", "-no_progress"],
        &library("expected/FileInfo_empty_mzML.txt"),
        &library("expected/FileInfo_empty_mzML.tsv"),
    );
    let strict = Options {
        meta: true,
        processing: true,
        statistics: true,
        ..Options::default()
    };
    assert!(!strict.source_dangling_references);
    let error = FileInfoLibrary::new().run(&input, &strict).unwrap_err();
    assert!(
        error.to_string().contains("unresolved dataProcessingRef"),
        "{error}"
    );
    let lenient = Options {
        source_dangling_references: true,
        ..strict
    };
    let result = FileInfoLibrary::new().run(&input, &lenient).unwrap();
    assert_report(&result.text, &library("expected/FileInfo_empty_mzML.txt"));
}

/// C1 `FileInfo_on_FFC_1_oracle_output` and `FileInfo_on_FFC_1_retained_output`:
/// `-m -p -s` on the C++ FeatureFinderCentroided_1 outputs.
#[test]
fn c1_feature_finder_centroided_1_outputs() {
    let oracle = tool("inputs/FeatureFinderCentroided_1_oracle_output.featureXML");
    check_text_and_tsv(
        &["-test", "-in", &oracle, "-m", "-p", "-s", "-no_progress"],
        &tool("expected/FileInfo_on_FFC_1_oracle_output.txt"),
        &tool("expected/FileInfo_on_FFC_1_oracle_output.tsv"),
    );
    let retained = data("mzml_mobility/FeatureFinderCentroided_1_1_output.featureXML");
    check_text_and_tsv(
        &["-test", "-in", &retained, "-m", "-p", "-s", "-no_progress"],
        &library("expected/FileInfo_on_FFC_1_retained_output.txt"),
        &library("expected/FileInfo_on_FFC_1_retained_output.tsv"),
    );
}

/// C1 `FileInfo_out_dat`: an extension no type claims is accepted for `-out`.
#[test]
fn c1_out_with_an_unclaimed_extension() {
    let dir = Workdir::new();
    let out = dir.file("rep.dat");
    let input = library("inputs/FileInfo_1_input.dta");
    let outcome = run(&[
        "-test",
        "-in",
        &input,
        "-in_type",
        "dta",
        "-no_progress",
        "-out",
        &out,
    ]);
    assert_code(&outcome, ExitCode::ExecutionOk);
    assert_report(&read(&out), &library("expected/FileInfo_1_tsv.txt"));
}

/// C1 `FileInfo_notype`: an undetermined type warns in the input check, then
/// fails with exit 10; `-out` was opened first and stays empty.
#[test]
fn c1_unknown_type_exits_10() {
    let dir = Workdir::new();
    let out = dir.file("notype.tmp.txt");
    let input = data("mzml_mobility/FileInfo_8_input.notype");
    let outcome = run(&["-test", "-in", &input, "-no_progress", "-out", &out]);
    assert_code(&outcome, ExitCode::ParseError);
    assert_eq!(
        outcome.err,
        format!(
            "Warning: Could not determine format of input file '{input}'!\nError: Could not determine input file type!\n"
        )
    );
    assert!(outcome.out.is_empty());
    assert_eq!(read(&out), "");
}

/// Oracle `notype_prefilled_out`: the source opens `-out` before resolving the
/// type, so an existing file is truncated even though the run fails.
#[test]
fn an_existing_out_is_truncated_by_a_failing_run() {
    let dir = Workdir::new();
    let out = dir.file("prefilled.txt");
    fs::write(&out, "previous content\n").unwrap();
    let input = data("mzml_mobility/FileInfo_8_input.notype");
    let outcome = run(&["-test", "-in", &input, "-no_progress", "-out", &out]);
    assert_code(&outcome, ExitCode::ParseError);
    assert_eq!(read(&out), "");
}

/// C1 `FileInfo_index_on_dta`: `-i` on a DTA file prints the error line, then
/// the usage text, and exits 6 with nothing on the output stream but the
/// closing line: `outputTo_` returns `ILLEGAL_PARAMETERS`, so `TOPPBase`
/// prints `FileInfo took …`, in the C1 product-SDK run as in the Release
/// build (`release_fileinfo_refusals_end_as_in_the_release_build`).
#[test]
fn c1_index_check_on_a_non_mzml_file_exits_6_with_usage() {
    let input = library("inputs/FileInfo_1_input.dta");
    let outcome = run(&[
        "-test",
        "-in",
        &input,
        "-in_type",
        "dta",
        "-i",
        "-no_progress",
    ]);
    assert_code(&outcome, ExitCode::IllegalParameters);
    assert!(outcome.out.is_empty(), "{}", outcome.out);
    assert!(outcome.took.is_some(), "main_ returned 6");
    assert_eq!(
        outcome.err,
        cpp_usage_stream(&tool("expected/FileInfo_index_on_dta.stderr.txt"))
    );
}

/// C1 `FileInfo_in_type_foo`: an invalid `-in_type` fails the strict update.
#[test]
fn c1_invalid_in_type_exits_6() {
    let dir = Workdir::new();
    let out = dir.file("in_type_foo.tmp.txt");
    let input = library("inputs/FileInfo_1_input.dta");
    let outcome = run(&[
        "-test",
        "-in",
        &input,
        "-in_type",
        "foo",
        "-no_progress",
        "-out",
        &out,
    ]);
    assert_code(&outcome, ExitCode::IllegalParameters);
    for line in [
        "Parameters passed to 'FileInfo' are invalid. To prevent usage of wrong defaults, please update/fix the parameters!\n",
        "Invalid string parameter value 'foo' for parameter 'in_type' given! Valid values are: 'mzData,mzXML,mzML,sqMass,dta,dta2d,mgf,featureXML,consensusXML,idXML,pepXML,mzTab,fid,mzid,trafoXML,fasta,pqp'. Updating failed!\n",
    ] {
        assert!(outcome.err.contains(line), "{}", outcome.err);
    }
    assert!(!Path::new(&out).exists());
}

/// C1 `FileInfo_out_csv` and `FileInfo_out_tsv_txt`: a claimed extension that
/// is not the registered one exits 6.
#[test]
fn c1_wrong_output_extensions_exit_6() {
    let dir = Workdir::new();
    let input = library("inputs/FileInfo_1_input.dta");
    let csv = dir.file("report.csv");
    let outcome = run(&[
        "-test",
        "-in",
        &input,
        "-in_type",
        "dta",
        "-no_progress",
        "-out",
        &csv,
    ]);
    assert_code(&outcome, ExitCode::IllegalParameters);
    assert_eq!(
        outcome.err,
        format!(
            "Invalid parameter: Invalid output file extension for file '{csv}'. Valid file extensions are: 'txt'.\n"
        )
    );
    let (text, tsv_txt) = (dir.file("report.txt"), dir.file("report_tsv.txt"));
    let outcome = run(&[
        "-test",
        "-in",
        &input,
        "-in_type",
        "dta",
        "-no_progress",
        "-out",
        &text,
        "-out_tsv",
        &tsv_txt,
    ]);
    assert_code(&outcome, ExitCode::IllegalParameters);
    assert_eq!(
        outcome.err,
        format!(
            "Invalid parameter: Invalid output file extension for file '{tsv_txt}'. Valid file extensions are: 'tsv'.\n"
        )
    );
}

/// C1 `FileInfo_out_unwritable`: exit 5.
#[test]
fn c1_unwritable_out_exits_5() {
    let input = library("inputs/FileInfo_1_input.dta");
    let outcome = run(&[
        "-test",
        "-in",
        &input,
        "-in_type",
        "dta",
        "-no_progress",
        "-out",
        "/nonexistent_dir/report.txt",
    ]);
    assert_code(&outcome, ExitCode::CannotWriteOutputFile);
    assert_eq!(
        outcome.err,
        "Cannot write output file given from parameter '-out'!\nError: Unable to write file (the file '/nonexistent_dir/report.txt' could not be created. )\n"
    );
}

/// C1 `FileInfo_in_not_found` (exit 1) and `FileInfo_missing_in` (exit 7).
///
/// The source reads `-in` only inside `outputTo_`, after opening `-out`, so its
/// missing-`-in` run leaves an empty `-out`; the framework checks every option
/// before the tool body, so no `-out` is created here
/// (`docs/TOPP_CLI_SUPPORT.md`, *Validation order*).
#[test]
fn c1_missing_input_file_and_parameter() {
    let dir = Workdir::new();
    let missing = dir.file("does_not_exist.mzML");
    let outcome = run(&["-test", "-in", &missing, "-no_progress"]);
    assert_code(&outcome, ExitCode::InputFileNotFound);
    assert!(
        outcome.err.starts_with(&format!(
            "Cannot read input file given from parameter '-in'!\nError: File not found (the file '{missing}' "
        )),
        "{}",
        outcome.err
    );
    let out = dir.file("missing_in.tmp.txt");
    let outcome = run(&["-test", "-no_progress", "-out", &out]);
    assert_code(&outcome, ExitCode::MissingParameters);
    assert_eq!(
        outcome.err,
        "Error: The required parameter 'in' [valid: mzData, mzXML, mzML, sqMass, dta, dta2d, mgf, featureXML, consensusXML, idXML, pepXML, mzTab, fid, mzid, trafoXML, fasta, pqp] was not given or is empty!\n"
    );
    assert!(!Path::new(&out).exists());
}

/// C1 `FileInfo_truncated_featureXML` and `FileInfo_truncated_mzML`: the first
/// 3000 bytes of the FileInfo_3 and FileInfo_9 inputs, as C1 derived them,
/// exit 3; `-out` was opened first and stays empty.
#[test]
fn c1_truncated_inputs_exit_3() {
    // The provenance manifest records the sha256 of both prefixes, equal to
    // C1's cases/derived/corrupt.featureXML and corrupt.mzML.
    for (source, name) in [
        ("inputs/FileInfo_3_input.featureXML", "corrupt.featureXML"),
        ("inputs/FileInfo_9_input.mzML", "corrupt.mzML"),
    ] {
        let dir = Workdir::new();
        let input = dir.file(name);
        let bytes = fs::read(library(source)).unwrap();
        fs::write(&input, &bytes[..3000]).unwrap();
        let out = dir.file("corrupt.tmp.txt");
        let outcome = run(&["-test", "-in", &input, "-no_progress", "-out", &out]);
        assert_code(&outcome, ExitCode::InputFileCorrupt);
        assert!(
            outcome.err.starts_with("Error: Unable to read file ("),
            "{}",
            outcome.err
        );
        assert_eq!(read(&out), "");
    }
}

// ---------------------------------------------------------------------------
// Output routing, forced types and other tool-level cases (tier 1:
// ../oracle/topp-file-info-tool)
// ---------------------------------------------------------------------------

/// Oracle `stdout_dta`: without `-out` the report goes to the output stream,
/// exactly the C++ stdout before its timing line.
#[test]
fn without_out_the_report_goes_to_the_output_stream() {
    let input = library("inputs/FileInfo_1_input.dta");
    let outcome = run(&["-test", "-in", &input, "-in_type", "dta", "-no_progress"]);
    assert_code(&outcome, ExitCode::ExecutionOk);
    assert!(outcome.err.is_empty(), "{}", outcome.err);
    let expected = tool("expected/stdout_dta.stdout.txt");
    assert_report_text(&outcome.out, &cpp_stdout_report(&expected), &expected);
}

/// Oracle `stdout_featurexml_mps_tsv_only`: `-out_tsv` without `-out` writes the
/// TSV file and the text report to the output stream.
#[test]
fn out_tsv_without_out() {
    let dir = Workdir::new();
    let tsv = dir.file("only.tmp.tsv");
    let input = library("inputs/FileInfo_3_input.featureXML");
    let outcome = run(&[
        "-test",
        "-in",
        &input,
        "-m",
        "-p",
        "-s",
        "-no_progress",
        "-out_tsv",
        &tsv,
    ]);
    assert_code(&outcome, ExitCode::ExecutionOk);
    let expected = tool("expected/stdout_featurexml_mps_tsv_only.stdout.txt");
    assert_report_text(&outcome.out, &cpp_stdout_report(&expected), &expected);
    assert_report(
        &read(&tsv),
        &tool("expected/stdout_featurexml_mps_tsv_only.tsv"),
    );
}

/// Oracle `forced_featurexml_tmp` and `detected_featurexml_tmp`: a featureXML
/// map under a `.tmp` name runs the featureXML branch when forced, and when
/// its type comes from its content.
#[test]
fn a_featurexml_map_under_a_tmp_name() {
    for (case, name, forced) in [
        ("forced_featurexml_tmp", "forced.tmp", true),
        ("detected_featurexml_tmp", "detected.tmp", false),
    ] {
        let dir = Workdir::new();
        let input = dir.file(name);
        fs::copy(library("inputs/FileInfo_3_input.featureXML"), &input).unwrap();
        let mut args = vec!["-test", "-in", &input];
        if forced {
            args.extend(["-in_type", "featureXML"]);
        }
        args.extend(["-m", "-p", "-s", "-no_progress"]);
        let (text, tsv) = (dir.file("report.tmp.txt"), dir.file("report.tmp.tsv"));
        args.extend(["-out", &text, "-out_tsv", &tsv]);
        let outcome = run(&args);
        assert_code(&outcome, ExitCode::ExecutionOk);
        assert!(
            !outcome.err.contains("Could not determine format"),
            "{}",
            outcome.err
        );
        assert_report(&read(&text), &tool(&format!("expected/{case}.txt")));
        assert_report(&read(&tsv), &tool(&format!("expected/{case}.tsv")));
    }
}

/// Oracle `threads_4`: the tool is serial, as the source; `-threads` is
/// accepted and every thread count writes the same bytes.
#[test]
fn threads_do_not_change_the_reports() {
    let input = library("inputs/FileInfo_3_input.featureXML");
    let mut reports = Vec::new();
    for threads in ["1", "4", "0"] {
        let dir = Workdir::new();
        let (text, tsv) = (dir.file("threads.tmp.txt"), dir.file("threads.tmp.tsv"));
        let outcome = run(&[
            "-test",
            "-in",
            &input,
            "-m",
            "-p",
            "-s",
            "-no_progress",
            "-threads",
            threads,
            "-out",
            &text,
            "-out_tsv",
            &tsv,
        ]);
        assert_code(&outcome, ExitCode::ExecutionOk);
        reports.push((read(&text), read(&tsv)));
    }
    assert_report(&reports[1].0, &tool("expected/threads_4.txt"));
    assert_report(&reports[1].1, &tool("expected/threads_4.tsv"));
    assert!(reports.iter().all(|report| *report == reports[0]));
}

/// Oracle `d_c_on_featurexml`: `-d` and `-c` have no effect on a featureXML map.
#[test]
fn detailed_and_corrupt_flags_on_featurexml() {
    let dir = Workdir::new();
    let out = dir.file("dc.tmp.txt");
    let input = library("inputs/FileInfo_3_input.featureXML");
    let outcome = run(&[
        "-test",
        "-in",
        &input,
        "-d",
        "-c",
        "-no_progress",
        "-out",
        &out,
    ]);
    assert_code(&outcome, ExitCode::ExecutionOk);
    assert_report(&read(&out), &tool("expected/d_c_on_featurexml.txt"));
}

/// Oracle `out_is_directory`: `-out` naming an existing directory passes the
/// writability check, and the open fails as the source's `FileNotWritable`,
/// exit 8. The exception unwinds past `TOPPBase`'s closing line, so the
/// output stream is empty, in the product-SDK oracle as in the Release build
/// (`release_fileinfo_refusals_end_as_in_the_release_build`).
#[test]
fn out_naming_a_directory_exits_8() {
    let dir = Workdir::new();
    let out = dir.file("dir.txt");
    fs::create_dir(&out).unwrap();
    let input = library("inputs/FileInfo_1_input.dta");
    let outcome = run(&[
        "-test",
        "-in",
        &input,
        "-in_type",
        "dta",
        "-no_progress",
        "-out",
        &out,
    ]);
    assert_code(&outcome, ExitCode::UnknownError);
    assert_eq!(
        outcome.err,
        format!(
            "Error: Unexpected internal error (the file '{out}' is not writable for the current user)\n"
        )
    );
    assert_eq!(outcome.out, "");
    assert_eq!(outcome.took, None);
}

/// Oracle `in_is_directory`: a directory as `-in` has no type; exit 10, an
/// exit code `outputTo_` returns, so the closing line follows, and the
/// output stream holds nothing else. That oracle ran the macOS product SDK,
/// whose libc++ reads a directory as an empty file; the Linux Release build's
/// libstdc++ throws instead (see
/// `a_directory_input_is_a_recorded_difference_from_the_linux_release_build`).
#[test]
fn in_naming_a_directory_exits_10() {
    let dir = Workdir::new();
    let input = dir.file("somedir");
    fs::create_dir(&input).unwrap();
    let outcome = run(&["-test", "-in", &input, "-no_progress"]);
    assert_code(&outcome, ExitCode::ParseError);
    assert_eq!(
        outcome.err,
        format!(
            "Warning: Could not determine format of input file '{input}'!\nError: Could not determine input file type!\n"
        )
    );
    assert_eq!(outcome.out, "");
    assert!(outcome.took.is_some(), "main_ returned 10");
}

/// The Release build (`../oracle/topp-exception-exits`, retained in
/// `tests/data/topp_exception_exits`) on the tool's refusals, with `-log`:
///
/// * `fi_index_on_dta_log`: `-i` on a DTA file. `outputTo_` writes the error
///   with `writeLogError_` (`FileInfo.cpp:118-121`) and returns 6, so the log
///   file holds the line, the usage text follows on standard error, and the
///   closing line on standard output.
/// * `fi_notype_log`: a file of no type. The framework's warning and the
///   tool's `writeLogError_` line reach the log, and `outputTo_` returns 10:
///   the closing line follows.
/// * `fi_out_is_directory_log`: `-out` an existing directory. The source
///   throws `FileNotWritable`: the `BaseException` arm writes its line to the
///   error stream and the log, exit 8, and no closing line.
#[test]
fn release_fileinfo_refusals_end_as_in_the_release_build() {
    use release_runs::ReleaseRun;
    ReleaseRun::new("fi_index_on_dta_log").assert_replayed::<FileInfo>(&[], |_| {});
    ReleaseRun::new("fi_notype_log").assert_replayed::<FileInfo>(&[], |_| {});
    ReleaseRun::new("fi_out_is_directory_log").assert_replayed::<FileInfo>(&[], |cwd| {
        fs::create_dir(cwd.join("dir.txt")).unwrap();
    });
}

/// Oracle `fi_debug2_log`: at debug level 2 the tool's own `writeDebug_` line
/// `Input file type: dta` (`FileInfo.cpp:111`) reaches the log file, where
/// the Release build writes it, right after the `in_type` read that decides
/// the type is detected. The framework's lines up to the tool body are
/// compared in order; the body's `Value of … option` lines follow the order
/// in which the port reads its options, and the input check runs before the
/// body (`docs/TOPP_CLI_SUPPORT.md`, *Log file and debug levels*), so the
/// rest is compared as a collection. The report on standard output is the
/// Release build's, and so is the closing line.
#[test]
fn release_fileinfo_debug_lines_reach_the_log() {
    use release_runs::ReleaseRun;
    let case = ReleaseRun::new("fi_debug2_log");
    let replay = case.replay::<FileInfo>(&[], |_| {});
    let release = case.release("FileInfo", &replay, &[]);
    assert_eq!(replay.code.as_i32(), release.exit, "{}", replay.err);
    assert_eq!(replay.out, release.out);
    assert_eq!(replay.took.is_some(), release.took.is_some());
    assert_eq!(replay.err, release.err);
    let (actual, expected) = (replay.log.unwrap(), release.log.unwrap());
    let framework = 1 + expected
        .iter()
        .position(|line| line.contains("Value of string option 'no_progress'"))
        .unwrap();
    assert_eq!(actual[..framework], expected[..framework]);
    let line = "<time> FileInfo:1:: Input file type: dta\n";
    assert!(expected.iter().any(|l| l == line));
    let (mut rest_actual, mut rest_expected) =
        (actual[framework..].to_vec(), expected[framework..].to_vec());
    rest_actual.sort();
    rest_expected.sort();
    assert_eq!(rest_actual, rest_expected);
}

/// Oracle `fi_undetermined_log`: a directory as `-in` on the Linux Release
/// build. `TOPPBase`'s input check asks `FileHandler::getType`, whose content
/// sniffing reads the directory; libstdc++ throws `std::ios_base::failure`,
/// which only the initialisation catch handles: `Unable to initialize or run
/// FileInfo: basic_filebuf::underflow error reading the file: Is a
/// directory`, exit 12, nothing logged. The macOS product SDK (libc++) reads
/// the directory as an empty file instead and ends with the tool's own
/// refusal, exit 10 (`in_naming_a_directory_exits_10`), which is what this
/// port does. Recorded here, and in `docs/TOPP_CLI_SUPPORT.md` (*Input
/// checks on a directory*), as a difference the port has not closed.
#[test]
fn a_directory_input_is_a_recorded_difference_from_the_linux_release_build() {
    use release_runs::ReleaseRun;
    let case = ReleaseRun::new("fi_undetermined_log");
    let replay = case.replay::<FileInfo>(&[], |cwd| {
        fs::create_dir(cwd.join("somedir")).unwrap();
    });
    let release = case.release("FileInfo", &replay, &[]);
    assert_eq!(release.exit, ExitCode::InternalError.as_i32());
    assert_eq!(
        release.err,
        "Unable to initialize or run FileInfo: basic_filebuf::underflow error reading the file: Is a directory\n"
    );
    assert_eq!((release.out.as_str(), release.took.is_some()), ("", false));
    assert_eq!(release.log, None);
    // The port follows the product SDK: the framework warns, the tool refuses.
    assert_eq!(replay.code, ExitCode::ParseError);
    assert!(replay.took.is_some());
}

/// Oracle `zero_byte_input`: exit 4, from the framework's input check. The
/// source reads `-in` after opening `-out` and leaves an empty `-out`; the
/// framework checks first, so none is created (see
/// `c1_missing_input_file_and_parameter`).
#[test]
fn a_zero_byte_input_exits_4() {
    let dir = Workdir::new();
    let input = dir.file("zero.dta");
    fs::write(&input, "").unwrap();
    let out = dir.file("zero.tmp.txt");
    let outcome = run(&["-test", "-in", &input, "-no_progress", "-out", &out]);
    assert_code(&outcome, ExitCode::InputFileEmpty);
    assert_eq!(
        outcome.err,
        format!(
            "Cannot read input file given from parameter '-in'!\nError: File empty (the file '{input}' is empty)\n"
        )
    );
    assert!(!Path::new(&out).exists());
}

/// Oracle `forced_dta_on_featurexml`, `forced_mzml_on_featurexml` and
/// `forced_featurexml_on_dta`: a forced type the loader's own detection
/// contradicts. The source loader throws: `ParseError` for the experiment
/// loader (exit 3) and `InvalidFileType` for the feature loader (exit 8). The
/// library maps the first to `Error::InvalidValue` (exit 6) and lets the
/// feature loader refuse the second itself; both are documented native
/// differences (`docs/TOPP_FILE_INFO_SUPPORT.md`). No report is written.
#[test]
fn a_forced_type_the_file_contradicts_is_refused() {
    let featurexml = library("inputs/FileInfo_3_input.featureXML");
    let dta = library("inputs/FileInfo_1_input.dta");
    for (input, forced, code) in [
        (&featurexml, "dta", ExitCode::IllegalParameters),
        (&featurexml, "mzML", ExitCode::IllegalParameters),
        (&dta, "featureXML", ExitCode::IncompatibleInputData),
    ] {
        let dir = Workdir::new();
        let out = dir.file("mismatch.tmp.txt");
        let outcome = run(&[
            "-test",
            "-in",
            input,
            "-in_type",
            forced,
            "-no_progress",
            "-out",
            &out,
        ]);
        assert_code(&outcome, code);
        assert!(outcome.out.is_empty(), "{}", outcome.out);
        assert_eq!(read(&out), "", "{forced}");
    }
}

// ---------------------------------------------------------------------------
// Branches the preview does not run (tier 4): explicit refusal, exit 11
// ---------------------------------------------------------------------------

/// Each unported registration exits `INCOMPATIBLE_INPUT_DATA` with a message
/// naming the branch or flag, writes nothing to the output stream and leaves
/// `-out` empty. The C++ tool exits 0 on all of them (C1 `TOPP_FileInfo_14`
/// and `_16`). A8 wired mzXML and mzData, so the three rows TOPP_FileInfo_4, _5
/// and _6 held here are gone; `topp_file_info_4_5_and_6` reproduces them.
#[test]
fn unported_branches_are_refused_explicitly() {
    let cases: Vec<(&str, Vec<String>, &str)> = vec![
        (
            "TOPP_FileInfo_14",
            args(&tool("inputs/FileInfo_14_input.mzid"), &["-v"]),
            "FileInfo schema and semantic validation (-v) is not ported",
        ),
        (
            "TOPP_FileInfo_16",
            args(&tool("inputs/FileInfo_16_input.trafoXML"), &[]),
            "FileInfo trafoXML branch is not ported",
        ),
        (
            "validation of a featureXML map",
            args(&library("inputs/FileInfo_3_input.featureXML"), &["-v"]),
            "FileInfo schema and semantic validation (-v) is not ported",
        ),
    ];
    for (case, case_args, message) in cases {
        let dir = Workdir::new();
        let out = dir.file("unported.tmp.txt");
        let mut full: Vec<&str> = vec!["-test", "-no_progress", "-out", &out];
        full.extend(case_args.iter().map(String::as_str));
        let outcome = run(&full);
        assert_code(&outcome, ExitCode::IncompatibleInputData);
        assert!(
            outcome
                .err
                .ends_with(&format!("Error: unsupported: {message}\n")),
            "{case}: {}",
            outcome.err
        );
        assert!(outcome.out.is_empty(), "{case}: {}", outcome.out);
        assert_eq!(read(&out), "", "{case}");
    }
}

/// A6 implemented `-i`, `-d` and `-c`, so the four rows this table used to hold
/// for them are gone. What they do instead is compared with the Release C++
/// output in `tests/file_info_checks.rs`; pinned here is only that the tool no
/// longer answers them with a not-ported refusal, so the table cannot quietly
/// regain them.
#[test]
fn the_flags_a6_implemented_are_no_longer_refused() {
    let faims = data("mzml_mobility/FAIMS_test_data.mzML");
    // TOPP_FileInfo_11 runs FileInfo_11_input.mzML, byte-identical to the
    // FileInfo_9 input (sha256 e14e087c...).
    let fi11 = library("inputs/FileInfo_9_input.mzML");
    let fi12 = library("inputs/FileInfo_12_input.mzML");
    let fi1 = library("inputs/FileInfo_1_input.dta");
    for (case, case_args) in [
        ("TOPP_FileInfo_11", args(&fi11, &["-i"])),
        ("TOPP_FileInfo_12", args(&fi12, &["-i"])),
        ("TOPP_FileInfo_19", args(&faims, &["-d"])),
        (
            "corrupt-data check on a peak file",
            args(&fi1, &["-in_type", "dta", "-c"]),
        ),
    ] {
        let dir = Workdir::new();
        let out = dir.file("a6.tmp.txt");
        let mut full: Vec<&str> = vec!["-test", "-no_progress", "-out", &out];
        full.extend(case_args.iter().map(String::as_str));
        let outcome = run(&full);
        for flag in ["(-i)", "(-d)", "(-c)"] {
            assert!(
                !outcome.err.contains("is not ported") || !outcome.err.contains(flag),
                "{case}: still refused: {}",
                outcome.err
            );
        }
    }
}

/// `-in <input> flags...` for the refusal table.
fn args(input: &str, flags: &[&str]) -> Vec<String> {
    ["-in", input]
        .iter()
        .chain(flags)
        .map(|part| (*part).to_owned())
        .collect()
}

// ---------------------------------------------------------------------------
// The executable (tier 4)
// ---------------------------------------------------------------------------

/// The shipped `FileInfo` binary routes the report to the process's standard
/// output and returns the tool's exit status.
#[test]
fn the_executable_writes_the_report_to_standard_output() {
    let binary = PathBuf::from(env!("CARGO_BIN_EXE_FileInfo"));
    let input = library("inputs/FileInfo_1_input.dta");
    let output = std::process::Command::new(&binary)
        .args(["-test", "-in", &input, "-in_type", "dta", "-no_progress"])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{output:?}");
    let (stdout, took) =
        took_line::split_took_line("FileInfo", &String::from_utf8(output.stdout).unwrap());
    assert!(took.is_some(), "the closing line ends the output: {stdout}");
    let expected = tool("expected/stdout_dta.stdout.txt");
    assert_report_text(&stdout, &cpp_stdout_report(&expected), &expected);

    // A7 ported the consensusXML branch, so the executable now writes its
    // report to standard output like any other branch instead of refusing.
    // Only when this build carries the format: without the feature the tool
    // answers "this build lacks the consensusxml feature", which is the
    // build's own refusal and is asserted by the feature-sliced tests.
    #[cfg(feature = "consensusxml")]
    {
        let consensus = tool("inputs/FileInfo_7_input.consensusXML");
        let output = std::process::Command::new(&binary)
            .args(["-test", "-in", &consensus, "-no_progress"])
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(0), "{output:?}");
        let stdout = String::from_utf8(output.stdout).unwrap();
        assert!(stdout.contains("File type: consensusXML\n"), "{stdout}");
        assert!(
            stdout.contains("Number of consensus features:\n"),
            "{stdout}"
        );
        assert!(
            stdout.contains("Assigned peptide identifications: 0\n"),
            "{stdout}"
        );
    }
}
