// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The FuzzyDiff TOPP tool (`OpenMS4-topp/src/FuzzyDiff.cpp`, topp `174b576`),
//! its executable, and the test-support emulation of its contract.
//!
//! Evidence, in order of strength (see `tests/data/topp_fuzzy_diff_provenance.json`
//! and `docs/TOPP_FUZZY_DIFF_SUPPORT.md`):
//!
//! - tier 1, executed differential: the pinned C++ Release build
//!   `openms4-release-bc9cc12-c19e494-174b576` ran `FuzzyDiff` on ibminode06
//!   for 80 invocations (`../oracle/fuzzy-diff-tool`, two runs, identical after
//!   normalisation): the four upstream registrations TOPP_FuzzyDiff_1..4
//!   (test-data `0cb15f2`, `topp/CMakeLists.txt:138-144`) on their pinned
//!   inputs, the 35 further invocations the product-SDK oracle judged by exit
//!   code only, and 41 that cover the usage text, `-write_ini`, parameter
//!   errors, every verbose level, tab width and first column, both whitelists,
//!   `-sort`, number tokens, raw bytes, relative paths and directories. Each
//!   case is run through the tool here and compared on the exit code, the
//!   output stream byte for byte and the error stream byte for byte, after the
//!   normalisations [`cpp_streams_for_port`] names;
//! - tier 4, the native refusals (an unreadable `-sort` input, the input bound)
//!   and the agreement of the test-support emulation with the tool.
//!
//! The four upstream registrations have no retained output: `WILL_FAIL` on 1, 2
//! and 4 is their whole expectation (`CMakeLists.txt:142-144`). The executed
//! run pins more than that - the exit code and both streams - and the tests
//! below check the registered expectation first and the executed streams second.
//!
//! Every case runs in its own temporary directory where it writes anything.

#![cfg(feature = "paramxml")]

#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;

use openms::cli::tools::FuzzyDiff;
use openms::cli::{ExitCode, Tool, run_with};
use openms::system::file::TempDir;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

// ---------------------------------------------------------------------------
// The oracle
// ---------------------------------------------------------------------------

/// The console-width probe the C++ usage text runs, which fails without a
/// terminal on Linux; the port has none (`docs/TOPP_CLI_SUPPORT.md`, *Usage
/// text*).
const LINUX_STTY_LINE: &str = "stty: 'standard input': Inappropriate ioctl for device\n";

/// The `TOPPBase` timing line after `main_` returns, as the oracle masks it.
/// The port's framework does not write it (`docs/TOPP_CLI_SUPPORT.md`).
const TIMING_LINE: &[u8] = b"FuzzyDiff took <T>.\n";

/// One invocation of `../oracle/fuzzy-diff-tool/cases.tsv`.
#[derive(Clone, Debug)]
struct Case {
    name: String,
    /// `-` for a fresh directory, otherwise a directory relative to `<D>`.
    cwd: String,
    args: Vec<String>,
}

/// What the Release build did, normalised by the oracle: the data directory
/// is `<D>`, the case directory `<CASE>`, the timing line's values `<T>` and
/// the unique part of a `-sort` temporary file `<UNIQUE>`.
#[derive(Clone, Debug)]
struct Executed {
    exit: i32,
    stdout: Vec<u8>,
    stderr: Vec<u8>,
}

fn data_dir() -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .to_string_lossy()
        .into_owned()
}

fn oracle_file(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/topp_fuzzy_diff/oracle")
        .join(name)
}

fn cases() -> Vec<Case> {
    let text = fs::read_to_string(oracle_file("cases.tsv")).unwrap();
    text.lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let mut fields = line.split('\t').map(str::to_owned);
            Case {
                name: fields.next().unwrap(),
                cwd: fields.next().unwrap(),
                args: fields.collect(),
            }
        })
        .collect()
}

fn case(name: &str) -> Case {
    cases()
        .into_iter()
        .find(|case| case.name == name)
        .unwrap_or_else(|| panic!("no oracle case {name}"))
}

fn hex(text: &str) -> Vec<u8> {
    assert!(text.len() % 2 == 0, "odd hex length");
    (0..text.len())
        .step_by(2)
        .map(|at| u8::from_str_radix(&text[at..at + 2], 16).unwrap())
        .collect()
}

fn executed() -> BTreeMap<String, Executed> {
    let text = fs::read_to_string(oracle_file("results.tsv")).unwrap();
    text.lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(fields.len(), 4, "{line}");
            (
                fields[0].to_owned(),
                Executed {
                    exit: fields[1].parse().unwrap(),
                    stdout: hex(fields[2]),
                    stderr: hex(fields[3]),
                },
            )
        })
        .collect()
}

fn executed_case(name: &str) -> Executed {
    executed()
        .remove(name)
        .unwrap_or_else(|| panic!("no executed result for {name}"))
}

// ---------------------------------------------------------------------------
// Running the port
// ---------------------------------------------------------------------------

struct Run {
    code: i32,
    out: Vec<u8>,
    err: Vec<u8>,
}

fn replace(haystack: &[u8], needle: &[u8], replacement: &[u8]) -> Vec<u8> {
    if needle.is_empty() {
        return haystack.to_vec();
    }
    let mut out = Vec::with_capacity(haystack.len());
    let mut rest = haystack;
    while let Some(at) = rest.windows(needle.len()).position(|w| w == needle) {
        out.extend_from_slice(&rest[..at]);
        out.extend_from_slice(replacement);
        rest = &rest[at + needle.len()..];
    }
    out.extend_from_slice(rest);
    out
}

/// The arguments of a case with its placeholders resolved: `<D>` is this
/// repository's `tests/data`, `<OUT>` the case's temporary directory.
fn arguments(case: &Case, out_dir: &Path) -> Vec<String> {
    let data = data_dir();
    let out = out_dir.to_string_lossy();
    std::iter::once(FuzzyDiff::NAME.to_owned())
        .chain(case.args.iter().map(|arg| {
            if arg == "<EMPTY>" {
                String::new()
            } else {
                arg.replace("<D>", &data).replace("<OUT>", &out)
            }
        }))
        .collect()
}

/// Run a case through `run_with`, in process, and put the oracle's
/// placeholders back into both streams.
fn run_case(case: &Case, out_dir: &Path) -> Run {
    assert_eq!(case.cwd, "-", "{}: needs the executable", case.name);
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<FuzzyDiff>(&arguments(case, out_dir), &mut out, &mut err).as_i32();
    Run {
        code,
        out: normalise_port(&out, out_dir),
        err: normalise_port(&err, out_dir),
    }
}

fn normalise_port(bytes: &[u8], out_dir: &Path) -> Vec<u8> {
    let bytes = replace(bytes, data_dir().as_bytes(), b"<D>");
    replace(&bytes, out_dir.to_string_lossy().as_bytes(), b"<CASE>/out")
}

/// The Release build's streams as this port writes them, applying exactly
/// the framework's documented differences and nothing of the tool's:
///
/// * the timing line `TOPPBase::main` writes after `main_` returns is dropped
///   from the output stream; whether it was there is returned, so a caller can
///   hold it to the cases where `main_` returned;
/// * the `stty` probe line of the usage text is dropped;
/// * the framework's missing-file parenthetical reads `does not exist` where
///   the source says `could not be found` (`docs/TOPP_CLI_SUPPORT.md`,
///   *Diagnostics*); only the framework's two-line input-check diagnostic is
///   rewritten, not the tool's own `FileNotFound` text;
/// * when the strict parameter update fails, the C++ error stream carries
///   `Parameters passed to 'FuzzyDiff' are invalid...` *before* the diagnostic
///   that caused it, although `TOPPBase` writes the diagnostic first: the two
///   go through different log streams (`getGlobalLogWarn()` and
///   `OPENMS_LOG_ERROR`, `TOPPBase.cpp:338-343` at cli c19e494) and the error
///   stream is flushed first. The port's framework writes them in code order;
///   the C++ order is restored to that here, the lines themselves unchanged;
/// * a `-sort` temporary file `<CASE>/tmp/<name>.sorted.<UNIQUE>.tmp` is named
///   by the input it copies, as the port compares in memory and names the
///   inputs (`docs/TOPP_FUZZY_DIFF_SUPPORT.md`, native difference 1).
fn cpp_streams_for_port(case: &Case, cpp: &Executed) -> (Vec<u8>, Vec<u8>, bool) {
    let (mut stdout, timed) = match cpp.stdout.strip_suffix(TIMING_LINE) {
        Some(body) => (body.to_vec(), true),
        None => (cpp.stdout.clone(), false),
    };
    assert!(
        !stdout
            .windows(b"FuzzyDiff took".len())
            .any(|w| w == b"FuzzyDiff took"),
        "{}: a timing line that is not the last one",
        case.name
    );
    let mut stderr = replace(&cpp.stderr, LINUX_STTY_LINE.as_bytes(), b"");
    if stderr.starts_with(b"Cannot read input file given from parameter '-") {
        stderr = replace(&stderr, b"' could not be found)\n", b"' does not exist)\n");
    }
    const INVALID: &[u8] = b"Parameters passed to 'FuzzyDiff' are invalid. To prevent usage of wrong defaults, please update/fix the parameters!\n";
    if let Some(diagnostics) = stderr.strip_prefix(INVALID) {
        stderr = [diagnostics, INVALID].concat();
    }
    if case.args.iter().any(|arg| arg == "-sort") {
        for input in input_arguments(case) {
            let name = Path::new(&input).file_name().unwrap().to_string_lossy();
            let temporary = format!("<CASE>/tmp/{name}.sorted.<UNIQUE>.tmp");
            stdout = replace(&stdout, temporary.as_bytes(), input.as_bytes());
        }
    }
    (stdout, stderr, timed)
}

/// The `-in1` and `-in2` values of a case, placeholders kept.
fn input_arguments(case: &Case) -> Vec<String> {
    case.args
        .windows(2)
        .filter(|pair| pair[0] == "-in1" || pair[0] == "-in2")
        .map(|pair| pair[1].clone())
        .collect()
}

/// Cases where the port deliberately differs from the Release build; each is
/// asserted on its own below.
const KNOWN_DIVERGENCES: [&str; 3] = [
    "token_underflow",
    "token_underflow_absdiff",
    "sort_directory",
];

/// Cases that run with a working directory of their own, through the executable.
const PROCESS_CASES: [&str; 2] = ["relative_paths", "relative_same_name"];

/// The first differing line of two streams, for a readable failure.
fn first_difference(actual: &[u8], expected: &[u8]) -> String {
    let actual_lines: Vec<&[u8]> = actual.split_inclusive(|&b| b == b'\n').collect();
    let expected_lines: Vec<&[u8]> = expected.split_inclusive(|&b| b == b'\n').collect();
    for (index, line) in actual_lines.iter().enumerate() {
        match expected_lines.get(index) {
            Some(other) if other == line => {}
            other => {
                return format!(
                    "line {}: port {:?}, C++ {:?}",
                    index + 1,
                    String::from_utf8_lossy(line),
                    other.map(|o| String::from_utf8_lossy(o))
                );
            }
        }
    }
    match expected_lines.get(actual_lines.len()) {
        Some(line) => format!(
            "port ends early; C++ continues with {:?}",
            String::from_utf8_lossy(line)
        ),
        None => "no line differs".to_owned(),
    }
}

/// Compare one case with the Release build; returns what differs.
fn compare_case(case: &Case, cpp: &Executed) -> Vec<String> {
    let dir = TempDir::new(false).unwrap();
    let run = run_case(case, dir.path());
    let (stdout, stderr, _) = cpp_streams_for_port(case, cpp);
    let mut problems = Vec::new();
    if run.code != cpp.exit {
        problems.push(format!(
            "{}: exit {} where C++ exits {}\nport stdout:\n{}\nport stderr:\n{}",
            case.name,
            run.code,
            cpp.exit,
            String::from_utf8_lossy(&run.out),
            String::from_utf8_lossy(&run.err)
        ));
    }
    if run.out != stdout {
        problems.push(format!(
            "{}: stdout {}",
            case.name,
            first_difference(&run.out, &stdout)
        ));
    }
    if run.err != stderr {
        problems.push(format!(
            "{}: stderr {}",
            case.name,
            first_difference(&run.err, &stderr)
        ));
    }
    problems
}

fn assert_case_matches(name: &str) {
    let problems = compare_case(&case(name), &executed_case(name));
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}

// ---------------------------------------------------------------------------
// The four upstream registrations (test-data 0cb15f2, topp/CMakeLists.txt:138-144)
// ---------------------------------------------------------------------------

/// Run an upstream registration and check its registered expectation: a
/// `WILL_FAIL` test must exit non-zero, any other must exit 0.
fn upstream(name: &str, will_fail: bool) -> Run {
    let dir = TempDir::new(false).unwrap();
    let run = run_case(&case(name), dir.path());
    assert_eq!(
        run.code != 0,
        will_fail,
        "{name}: exit {}\n{}{}",
        run.code,
        String::from_utf8_lossy(&run.out),
        String::from_utf8_lossy(&run.err)
    );
    run
}

/// TOPP_FuzzyDiff_1 compares `FuzzyDiff_1_in1.featureXML` with itself and is
/// registered `WILL_FAIL`: the same name twice is refused as cheating. (The
/// comments at `CMakeLists.txt:142-143` have the reasons of 1 and 2 swapped.)
#[test]
fn topp_fuzzydiff_1_refuses_the_same_file_twice() {
    let run = upstream("TOPP_FuzzyDiff_1", true);
    assert_eq!(run.code, ExitCode::ParseError.as_i32());
    assert!(String::from_utf8_lossy(&run.out).contains("That's cheating!"));
    assert_case_matches("TOPP_FuzzyDiff_1");
}

/// TOPP_FuzzyDiff_2 is registered `WILL_FAIL`: a number differs beyond the
/// pinned `FuzzyDiff.ini` tolerances (ratio 1.01, absdiff 0.01).
#[test]
fn topp_fuzzydiff_2_fails_on_a_number() {
    let run = upstream("TOPP_FuzzyDiff_2", true);
    assert_eq!(run.code, ExitCode::ParseError.as_i32());
    assert!(String::from_utf8_lossy(&run.out).contains("FAILED: 'ratio of numbers is too large'"));
    assert_case_matches("TOPP_FuzzyDiff_2");
}

/// TOPP_FuzzyDiff_3 passes: two byte-identical files under different names.
#[test]
fn topp_fuzzydiff_3_passes() {
    let run = upstream("TOPP_FuzzyDiff_3", false);
    assert_eq!(run.code, ExitCode::ExecutionOk.as_i32());
    assert_case_matches("TOPP_FuzzyDiff_3");
}

/// TOPP_FuzzyDiff_4 is registered `WILL_FAIL`: `lorem_ipsum.featureXML` does
/// not exist at the pin, so the framework's input check exits 1.
#[test]
fn topp_fuzzydiff_4_exits_1_for_a_missing_input() {
    let run = upstream("TOPP_FuzzyDiff_4", true);
    assert_eq!(run.code, ExitCode::InputFileNotFound.as_i32());
    assert_case_matches("TOPP_FuzzyDiff_4");
}

// ---------------------------------------------------------------------------
// Every executed case
// ---------------------------------------------------------------------------

/// Every in-process oracle case against the Release build: exit code and both
/// streams, byte for byte after [`cpp_streams_for_port`]. The deliberate
/// divergences and the cases that need a working directory are checked in
/// their own tests; the count pins the case list.
#[test]
fn every_executed_case_matches_the_release_build() {
    let executed = executed();
    let all = cases();
    assert_eq!(all.len(), 80);
    assert_eq!(executed.len(), all.len());
    let mut problems = Vec::new();
    let mut compared = 0;
    for case in &all {
        if KNOWN_DIVERGENCES.contains(&case.name.as_str())
            || PROCESS_CASES.contains(&case.name.as_str())
        {
            continue;
        }
        problems.extend(compare_case(case, &executed[&case.name]));
        compared += 1;
    }
    assert_eq!(compared, 75);
    assert!(
        problems.is_empty(),
        "{} problem(s):\n{}",
        problems.len(),
        problems.join("\n\n")
    );
}

/// The oracle's own consistency: the timing line is written exactly when
/// `main_` returned, which is every exit 0 or 10 except usage and `-write_ini`,
/// which end before `main_`. This is what licenses dropping it.
#[test]
fn the_timing_line_marks_exactly_the_runs_whose_main_returned() {
    for case in cases() {
        let cpp = &executed()[&case.name];
        let (_, _, timed) = cpp_streams_for_port(&case, cpp);
        let before_main = ["help", "helphelp", "write_ini"].contains(&case.name.as_str());
        let returned = (cpp.exit == 0 || cpp.exit == 10) && !before_main;
        assert_eq!(timed, returned, "{}", case.name);
    }
}

/// `relative_paths` and `relative_same_name` run with a working directory:
/// the report prints relative names as absolute paths
/// (`std::filesystem::absolute`), and two equal relative names are cheating.
/// Run through the executable, whose streams are the process's.
#[test]
fn relative_names_resolve_against_the_working_directory() {
    let binary = env!("CARGO_BIN_EXE_FuzzyDiff");
    for name in PROCESS_CASES {
        let case = case(name);
        let cpp = executed_case(name);
        let dir = TempDir::new(false).unwrap();
        let arguments = arguments(&case, dir.path());
        let done = std::process::Command::new(binary)
            .args(&arguments[1..])
            .current_dir(Path::new(&data_dir()).join(&case.cwd))
            .output()
            .unwrap();
        let (stdout, stderr, _) = cpp_streams_for_port(&case, &cpp);
        assert_eq!(done.status.code(), Some(cpp.exit), "{name}");
        let out = normalise_port(&done.stdout, dir.path());
        let err = normalise_port(&done.stderr, dir.path());
        assert!(out == stdout, "{name}: {}", first_difference(&out, &stdout));
        assert!(err == stderr, "{name}: {}", first_difference(&err, &stderr));
    }
}

/// `-write_ini`: the file equals the Release build's, line by line with exact
/// numbers apart from nothing (the `version` item is 1.0.0 in both), and as a
/// decoded parameter tree with descriptions, tags and restrictions.
#[test]
fn write_ini_matches_the_release_build() {
    let dir = TempDir::new(false).unwrap();
    let written = dir.path().join("written.ini");
    let arguments = vec![
        FuzzyDiff::NAME.to_owned(),
        "-write_ini".to_owned(),
        written.to_string_lossy().into_owned(),
    ];
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<FuzzyDiff>(&arguments, &mut out, &mut err);
    assert_eq!(code, ExitCode::ExecutionOk);
    assert!(out.is_empty() && err.is_empty());
    let expected = oracle_file("write_ini.ini");

    let mut comparator = fuzzy::FuzzyStringComparator::new();
    comparator.set_acceptable_relative(1.0);
    comparator.set_acceptable_absolute(0.0);
    comparator.set_whitelist(Vec::new());
    comparator.set_log_destination(fuzzy::LogDestination::Buffer);
    assert!(
        comparator.compare_files(&written, &expected),
        "{}",
        String::from_utf8_lossy(comparator.log())
    );
    let produced = openms::format::paramxml::load(&written).unwrap();
    let oracle = openms::format::paramxml::load(&expected).unwrap();
    assert_eq!(produced, oracle);
}

// ---------------------------------------------------------------------------
// Deliberate divergences
// ---------------------------------------------------------------------------

/// `token_underflow` and `token_underflow_absdiff`: `value 1e-400` against
/// `value 0`, with absdiff 0 and 1. The Linux Release build's `std::from_chars`
/// (libstdc++) reports a value that rounds to zero as out of range, so `1e-400`
/// reads as the letter `1` and the comparison fails (exit 10) whatever the
/// tolerance. The source's own comment says underflow is accepted, and its
/// libc++ build accepts it (`FuzzyStringComparator.cpp:205-208`); the
/// comparator follows that stated contract, as it did before it was promoted,
/// so the port reads 0 against 0 and passes. The verdict therefore depends on
/// the C++ standard library. Only a value that rounds to zero differs: a
/// subnormal (`token_subnormal`), the smallest normal (`token_min_normal`) and
/// an overflow (`token_overflow`) are compared above, with the Release build's
/// verdicts. Recorded in `docs/TOPP_FUZZY_DIFF_SUPPORT.md`.
#[test]
fn an_underflowing_number_is_a_number_as_the_source_comment_says() {
    for (name, absdiff) in [("token_underflow", "0"), ("token_underflow_absdiff", "1")] {
        let cpp = executed_case(name);
        assert_eq!(cpp.exit, 10, "{name}");
        assert!(
            String::from_utf8_lossy(&cpp.stdout)
                .contains("input_1 is not a number, but input_2 is"),
            "{name}"
        );
        let dir = TempDir::new(false).unwrap();
        let run = run_case(&case(name), dir.path());
        assert_eq!(run.code, ExitCode::ExecutionOk.as_i32(), "{name}");
        assert_eq!(
            String::from_utf8_lossy(&run.out),
            format!(
                "PASSED.\n\n  relative_max:        1\n  relative_acceptable: 1\n\n  absolute_max:        0\n  absolute_acceptable: {absdiff}\n\nNo numeric differences were found.\n\n"
            ),
            "{name}"
        );
        assert!(run.err.is_empty(), "{name}");
    }
    for name in ["token_subnormal", "token_min_normal", "token_overflow"] {
        assert!(
            !KNOWN_DIVERGENCES.contains(&name),
            "{name} is compared in full"
        );
    }
}

/// `sort_directory`: a directory with `-sort`. The source's `std::getline`
/// swallows the read error, so the directory sorts as an empty text and is
/// compared (exit 10). The port does not compare a text it could not read:
/// `INTERNAL_ERROR` (12), the code the source gives the same directory without
/// `-sort` (`directory_vs_file`).
#[test]
fn a_sort_input_that_cannot_be_read_is_refused() {
    let cpp = executed_case("sort_directory");
    assert_eq!(cpp.exit, 10);
    let dir = TempDir::new(false).unwrap();
    let case = case("sort_directory");
    let run = run_case(&case, dir.path());
    assert_eq!(run.code, ExitCode::InternalError.as_i32());
    assert!(run.out.is_empty());
    assert_eq!(
        String::from_utf8_lossy(&run.err),
        "Unable to initialize or run FuzzyDiff: error reading the file '<D>/topp_fuzzy_diff/inputs/directory_a': Is a directory\n"
    );
    assert_eq!(executed_case("directory_vs_file").exit, 12);
}

// ---------------------------------------------------------------------------
// The test-support emulation
// ---------------------------------------------------------------------------

/// Settings and inputs of a case in the emulation's terms, for the options the
/// emulation models; `None` for a case that uses anything else.
fn emulated(case: &Case) -> Option<(fuzzy::FuzzyDiffSettings, PathBuf, PathBuf)> {
    let data = data_dir();
    let args: Vec<String> = case
        .args
        .iter()
        .map(|arg| {
            if arg == "<EMPTY>" {
                String::new()
            } else {
                arg.replace("<D>", &data)
            }
        })
        .collect();
    let is_option = |arg: &str| {
        arg.len() > 1 && arg.starts_with('-') && arg[1..2].chars().all(char::is_alphabetic)
    };
    let mut settings = fuzzy::FuzzyDiffSettings::registered_defaults();
    if let Some(at) = args.iter().position(|arg| arg == "-ini") {
        settings = fuzzy::FuzzyDiffSettings::load_ini(Path::new(&args[at + 1])).unwrap();
    }
    let (mut in1, mut in2) = (None, None);
    let mut index = 0;
    while index < args.len() {
        let option = args[index].as_str();
        let values: Vec<&str> = args[index + 1..]
            .iter()
            .take_while(|arg| !is_option(arg))
            .map(String::as_str)
            .collect();
        let one = || values.first().copied();
        match option {
            "-test" | "-ini" => {}
            "-in1" => in1 = one().map(PathBuf::from),
            "-in2" => in2 = one().map(PathBuf::from),
            "-ratio" => settings.ratio = one()?.parse().ok()?,
            "-absdiff" => settings.absdiff = one()?.parse().ok()?,
            "-verbose" => settings.verbose = one()?.parse().ok()?,
            "-tab_width" => settings.tab_width = one()?.parse().ok()?,
            "-first_column" => settings.first_column = one()?.parse().ok()?,
            "-sort" => settings.sort = true,
            "-whitelist" => settings = settings.with_whitelist(&values),
            "-matched_whitelist" => settings = settings.with_matched_whitelist(&values),
            _ => return None,
        }
        index += 1 + values.len();
    }
    Some((settings, in1?, in2?))
}

/// The emulation that tests without the `paramxml` feature use,
/// `tests/support/fuzzy_string_comparator.rs::fuzzy_diff`, gives the tool's
/// exit code, and the Release build's, on the 39 invocations it was built
/// against (the product-SDK oracle's) and on every newer case it models,
/// except where it keeps its own documented answer for an unreadable input.
#[test]
fn the_test_support_emulation_agrees_with_the_tool() {
    let executed = executed();
    let mut checked = 0;
    for case in cases() {
        if case.cwd != "-" || KNOWN_DIVERGENCES.contains(&case.name.as_str()) {
            continue;
        }
        let Some((settings, in1, in2)) = emulated(&case) else {
            continue;
        };
        let outcome = fuzzy::fuzzy_diff(&in1, &in2, &settings);
        let dir = TempDir::new(false).unwrap();
        let tool = run_case(&case, dir.path()).code;
        let cpp = executed[&case.name].exit;
        assert_eq!(tool, cpp, "{}", case.name);
        if cpp == ExitCode::InternalError.as_i32() {
            // A directory: the emulation's comparator stops with a log line
            // and a failed verdict, as before the tool existed.
            assert_eq!(
                outcome.exit,
                fuzzy::FuzzyDiffExit::ParseError,
                "{}",
                case.name
            );
        } else {
            assert_eq!(
                outcome.exit.code(),
                cpp,
                "{}: {}",
                case.name,
                outcome.log_text()
            );
        }
        checked += 1;
    }
    assert!(checked >= 60, "{checked}");
}

// ---------------------------------------------------------------------------
// Native bounds and refusals (tier 4)
// ---------------------------------------------------------------------------

/// `-sort` reads each input into memory; a file beyond the comparator's
/// `MAX_INPUT_BYTES` would be refused with 11 before reading. A sparse file of
/// that size costs no disk; where the filesystem cannot make one, the check is
/// skipped rather than faked.
#[test]
fn a_sort_input_beyond_the_bound_is_refused_before_reading() {
    let dir = TempDir::new(false).unwrap();
    let big = dir.path().join("big.tsv");
    let file = fs::File::create(&big).unwrap();
    if file.set_len(fuzzy::MAX_INPUT_BYTES + 1).is_err() {
        return;
    }
    let small = Path::new(&data_dir()).join("topp_fuzzy_diff/inputs/sort_a.tsv");
    let arguments = vec![
        FuzzyDiff::NAME.to_owned(),
        "-sort".to_owned(),
        "-in1".to_owned(),
        big.to_string_lossy().into_owned(),
        "-in2".to_owned(),
        small.to_string_lossy().into_owned(),
    ];
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<FuzzyDiff>(&arguments, &mut out, &mut err);
    assert_eq!(code, ExitCode::IncompatibleInputData);
    assert!(out.is_empty());
    assert_eq!(
        String::from_utf8_lossy(&err),
        format!(
            "Error: input file '{}' exceeds the comparison limit of {} bytes.\n",
            big.display(),
            fuzzy::MAX_INPUT_BYTES
        )
    );
}

/// Without `-sort` the comparator refuses the same file itself, and the tool
/// reports the refusal as 11 rather than as a difference (10).
#[test]
fn an_input_beyond_the_bound_is_not_reported_as_a_difference() {
    let dir = TempDir::new(false).unwrap();
    let big = dir.path().join("big.txt");
    let file = fs::File::create(&big).unwrap();
    if file.set_len(fuzzy::MAX_INPUT_BYTES + 1).is_err() {
        return;
    }
    let small = Path::new(&data_dir()).join("fuzzy_string_comparator/fuzzydiff/x_1.txt");
    let arguments = vec![
        FuzzyDiff::NAME.to_owned(),
        "-in1".to_owned(),
        small.to_string_lossy().into_owned(),
        "-in2".to_owned(),
        big.to_string_lossy().into_owned(),
    ];
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<FuzzyDiff>(&arguments, &mut out, &mut err);
    assert_eq!(code, ExitCode::IncompatibleInputData);
    assert!(out.is_empty());
    assert_eq!(
        String::from_utf8_lossy(&err),
        format!(
            "Error: input file '{}' exceeds the comparison limit of {} bytes.\n",
            big.display(),
            fuzzy::MAX_INPUT_BYTES
        )
    );
}

/// An INI `tab_width` beyond the `Int` the source's `getIntOption_` returns
/// never reaches the tool: the framework's ParamXML reader refuses the file
/// (exit 3) before the tool body runs, as the command line's `to_int32` does
/// for a command-line value. The tool's own `i32` conversion is therefore a
/// defensive check that no input reaches. Native (tier 4); the C++ response
/// to such an INI was not executed.
#[test]
fn an_ini_integer_beyond_i32_never_reaches_the_tool() {
    let dir = TempDir::new(false).unwrap();
    let ini = dir.path().join("wide.ini");
    let pinned =
        fs::read_to_string(Path::new(&data_dir()).join("fuzzy_string_comparator/FuzzyDiff.ini"))
            .unwrap();
    let wide = pinned.replace(
        r#"<ITEM name="tab_width" value="8""#,
        r#"<ITEM name="tab_width" value="4294967296""#,
    );
    assert_ne!(wide, pinned, "the pinned INI has a tab_width item");
    fs::write(&ini, wide).unwrap();
    let x = Path::new(&data_dir()).join("fuzzy_string_comparator/fuzzydiff/x_1.txt");
    let y = Path::new(&data_dir()).join("fuzzy_string_comparator/fuzzydiff/x_1001.txt");
    let arguments: Vec<String> = [
        FuzzyDiff::NAME,
        "-ini",
        &ini.to_string_lossy(),
        "-in1",
        &x.to_string_lossy(),
        "-in2",
        &y.to_string_lossy(),
    ]
    .iter()
    .map(|s| (*s).to_owned())
    .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<FuzzyDiff>(&arguments, &mut out, &mut err);
    assert_eq!(
        code,
        ExitCode::InputFileCorrupt,
        "{}",
        String::from_utf8_lossy(&err)
    );
    assert!(out.is_empty());
    assert!(String::from_utf8_lossy(&err).contains("4294967296"));
}
