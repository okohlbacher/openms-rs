// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! JSON `-ini` files against the C++ Release build: the retained runs of
//! `../oracle/json-ini` (`tests/data/topp_json_ini/`, provenance in
//! `tests/data/topp_json_ini_provenance.json`), executed on ibminode06 against
//! the Release installation at core `bc9cc12`, cli `c19e494`, topp `174b576`.
//!
//! The source reads an `-ini` file whose name is of type JSON through
//! `ParamJSONFile::load` onto the tool's defaults (`TOPPBase.cpp:283-296`
//! for a run, `2619-2629` for `-write_ini`). Each case is replayed through
//! [`run_with_registry`] in a fresh directory with the Release run's paths
//! mapped to the replay's, and compared on the exit code, standard output
//! (the closing `<tool> took …` line by shape), standard error, the INI file
//! `-write_ini` wrote and the `-log` file with its time stamps masked. What is
//! not compared as recorded is named where it is applied, each a difference
//! `docs/TOPP_CLI_SUPPORT.md` records: the source's debug line that carries a
//! C++ file and line; at debug level 2, the order of the lifecycle's file
//! checks and the tool body's `Value of` lines; the text `serde_json` gives a
//! malformed document, where `nlohmann::json` words the same failure
//! differently; and the text of a floating-point value in the written INI,
//! which is the parameter XML writer's (the parameters themselves are
//! compared whole). Two cases are not replayed at all ([`OPEN`]).
#![cfg(all(feature = "mzml", feature = "paramxml"))]

#[path = "support/release_runs.rs"]
mod release_runs;
#[path = "support/took_line.rs"]
mod took_line;

use openms::cli::tools::{FuzzyDiff, PeakPickerHiRes};
use openms::cli::{ExitCode, Tool, ToolHandler, ToolRegistrySources, run_with_registry};
use openms::system::file::TempDir;
use release_runs::{RELEASE_PREFIX, masked_log};
use std::fs;
use std::path::{Path, PathBuf};

/// Cases whose Release behaviour the port does not reach yet, outside the
/// command line and the JSON reader: an integer parameter that the JSON
/// document sets beyond 32 bits (`1e30` becomes `-9223372036854775808`). The
/// source's `ParamEntry::isValid` and `getIntOption_` cast it to `int`, and
/// its INI writer writes it; this port's `Param` refuses to cast it
/// (`src/param.rs`, `ParamEntry::validation_error`) and its parameter XML
/// writer refuses to write it, so the run ends in the initialisation catch.
/// The loaded value itself is the Release build's (`tests/param_json.rs`).
const OPEN: [&str; 2] = ["r_pphr_huge", "w_pphr_int_huge"];

/// Where the Release run kept its case directories.
const ORACLE_RESULTS: &str = "/ceph/ibmi/abi/oliver/oracle-runs/json-ini/results2";
/// Where the Release run read its inputs.
const ORACLE_INPUTS: &str = "/ceph/ibmi/abi/oliver/oracle-runs/json-ini/inputs";

fn repository(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

fn data(case: &str) -> PathBuf {
    repository("tests/data/topp_json_ini").join(case)
}

fn text(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

fn read(case: &str, file: &str) -> Option<String> {
    fs::read_to_string(data(case).join(file)).ok()
}

/// One replay of a retained case.
struct Replay {
    code: ExitCode,
    out: String,
    took: Option<String>,
    err: String,
    ini: Option<String>,
    log: Option<Vec<String>>,
    cwd: PathBuf,
    inputs: PathBuf,
    _work: TempDir,
}

impl Replay {
    /// The Release run's `file`, with its paths mapped to this replay's.
    fn release(&self, case: &str, file: &str) -> Option<String> {
        Some(
            read(case, file)?
                .replace(&format!("{ORACLE_RESULTS}/{case}/cwd"), &text(&self.cwd))
                .replace(ORACLE_INPUTS, &text(&self.inputs))
                .replace(&format!("{RELEASE_PREFIX}/bin/"), ""),
        )
    }
}

/// Runs the retained case `case` with tool `T` in a fresh directory: the
/// inputs are copied to `inputs/`, and a `*.json.in` document the case
/// wrote before it ran is written to `cwd/run.json` as the oracle's
/// `templated` did.
fn replay<T: Tool>(case: &str) -> Replay {
    let work = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let (home, cwd, inputs) = (
        work.path().join("home"),
        work.path().join("cwd"),
        work.path().join("inputs"),
    );
    for dir in [&home, &cwd, &inputs] {
        fs::create_dir_all(dir).unwrap();
    }
    for entry in fs::read_dir(data("inputs")).unwrap() {
        let entry = entry.unwrap();
        fs::copy(entry.path(), inputs.join(entry.file_name())).unwrap();
    }
    fs::copy(
        repository("tests/data/file_info/inputs/FileInfo_1_input.dta"),
        inputs.join("FileInfo_1_input.dta"),
    )
    .unwrap();
    let template = match case {
        "r_fi_run_override" => Some("fi_run_mzml_type.json.in"),
        _ if case.starts_with("r_fi_") => Some("fi_run.json.in"),
        _ => None,
    };
    if let Some(template) = template {
        let document = fs::read_to_string(inputs.join(template)).unwrap();
        fs::write(cwd.join("run.json"), document.replace("@IN@", &text(&inputs))).unwrap();
    }
    if case == "w_directory" {
        fs::create_dir_all(cwd.join("dir.json")).unwrap();
    }
    let mut replay = Replay {
        code: ExitCode::ExecutionOk,
        out: String::new(),
        took: None,
        err: String::new(),
        ini: None,
        log: None,
        cwd,
        inputs,
        _work: work,
    };
    let arguments: Vec<String> = replay
        .release(case, "argv.txt")
        .unwrap()
        .lines()
        .map(str::to_owned)
        .collect();
    assert_eq!(arguments[0], T::NAME, "{case}");
    let mut sources = ToolRegistrySources::from_environment().unwrap();
    sources.prefixes.clear();
    sources.internal_tools_path = None;
    sources.ttd_internal_path = None;
    sources.file_context.user_override = Some(home);
    let registry = ToolHandler::new(sources);
    let (mut out, mut err) = (Vec::new(), Vec::new());
    replay.code = run_with_registry::<T>(&arguments, &mut out, &mut err, &registry);
    let (out, took) = took_line::split_took_line(T::NAME, &String::from_utf8_lossy(&out));
    replay.out = out;
    replay.took = took;
    replay.err = String::from_utf8_lossy(&err).into_owned();
    replay.ini = fs::read_to_string(replay.cwd.join("out.ini")).ok();
    replay.log = fs::read_to_string(replay.cwd.join("log.txt"))
        .ok()
        .map(|content| masked_log(&content));
    replay
}

/// The Release log without its source-location debug lines, which carry the
/// C++ file, line and function of a caught exception and have no counterpart.
fn release_log(replay: &Replay, case: &str) -> Option<Vec<String>> {
    replay.release(case, "log.txt").map(|content| {
        masked_log(&content)
            .into_iter()
            .filter(|line| !line.contains(":: Error occurred in line "))
            .collect()
    })
}

/// The written INI: the same parameters, descriptions, tags and restrictions
/// (`Param`'s full equality after reading both), and the same text on every
/// line but those of a floating-point item, whose number the Release build
/// writes as `4.0`, `1.0e20` and `0.0:` where this port's writer writes `4`,
/// `100000000000000000000` and `0:` (the parameter XML writer's text, a gap
/// of `ParamXMLFile`, not of the command line).
fn assert_ini(case: &str, port: Option<&str>, release: Option<&str>) {
    let (Some(port), Some(release)) = (port, release) else {
        assert_eq!(port, release, "{case}: the written INI");
        return;
    };
    let read = |text: &str| openms::format::paramxml::read(text.as_bytes()).unwrap();
    assert_eq!(read(port), read(release), "{case}: the written parameters");
    let (port_lines, release_lines): (Vec<&str>, Vec<&str>) =
        (port.lines().collect(), release.lines().collect());
    assert_eq!(port_lines.len(), release_lines.len(), "{case}: INI lines");
    for (port_line, release_line) in port_lines.iter().zip(&release_lines) {
        if port_line != release_line {
            assert!(
                port_line.contains("type=\"double\"") && release_line.contains("type=\"double\""),
                "{case}:\n{port_line}\n{release_line}"
            );
        }
    }
}

/// The log file. At debug level 2 the lifecycle checks the input and output
/// files before the tool body reads them, where the source checks each when
/// the body first reads it (`docs/TOPP_CLI_SUPPORT.md`, "-log and -debug"):
/// the lines up to the `-no_progress` value, which hold every parameter
/// dump, are compared in order and the rest as a set of lines. When the input
/// check fails, the body never runs, so the Release build's `Value of string
/// option 'in'` line, which its body writes before the check, has no
/// counterpart.
fn assert_log(case: &str, port: Option<Vec<String>>, release: Option<Vec<String>>) {
    let (Some(port), Some(mut release)) = (port.clone(), release.clone()) else {
        assert_eq!(port, release, "{case}: the log file");
        return;
    };
    if !case.contains("debug2") {
        assert_eq!(port, release, "{case}: the log file");
        return;
    }
    let framework = 1 + release
        .iter()
        .position(|line| line.contains(":: Value of string option 'no_progress'"))
        .unwrap_or_else(|| panic!("{case}: no -no_progress line"));
    assert_eq!(port[..framework], release[..framework], "{case}: the log file");
    if case == "r_pphr_debug2_log" {
        let before = release.len();
        release.retain(|line| !line.ends_with(":: Value of string option 'in': /x/in.mzML\n"));
        assert_eq!(release.len(), before - 1, "{case}");
    }
    let (mut port_rest, mut release_rest) = (port[framework..].to_vec(), release[framework..].to_vec());
    port_rest.sort();
    release_rest.sort();
    assert_eq!(port_rest, release_rest, "{case}: the log file after the parameters");
}

/// Compares a replay with the Release run; `err` maps the Release standard
/// error before the comparison.
fn assert_case<T: Tool>(case: &str, err: impl Fn(String) -> String) {
    let replay = replay::<T>(case);
    let exit: i32 = read(case, "exit_code.txt").unwrap().trim().parse().unwrap();
    assert_eq!(
        replay.code.as_i32(),
        exit,
        "{case}: exit\nstdout:\n{}\nstderr:\n{}",
        replay.out,
        replay.err
    );
    let stdout = replay.release(case, "stdout.txt").unwrap();
    let (out, took) = took_line::split_took_line(T::NAME, &stdout);
    assert_eq!(replay.out, out, "{case}: stdout");
    assert_eq!(replay.took.is_some(), took.is_some(), "{case}: the closing line");
    assert_eq!(
        replay.err,
        err(replay.release(case, "stderr.txt").unwrap()),
        "{case}: stderr"
    );
    assert_ini(case, replay.ini.as_deref(), replay.release(case, "out.ini").as_deref());
    assert_log(case, replay.log.clone(), release_log(&replay, case));
}

/// The retained cases whose names start with `prefix` and that ran `tool`.
fn cases(prefix: &str, tool: &str) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(repository("tests/data/topp_json_ini"))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|name| name.starts_with(prefix))
        .filter(|name| {
            read(name, "argv.txt").is_some_and(|argv| argv.lines().next() == Some(tool))
        })
        .collect();
    names.sort();
    names
}

/// The cases whose document `serde_json` refuses with its own words: the
/// source's `nlohmann::json` names the same failure as
/// `[json.exception.parse_error.101] parse error at line <l>, column <c>:
/// syntax error while parsing <context> - <reason>`.
const SYNTAX_ERRORS: [&str; 5] = [
    "w_pphr_syntax",
    "w_pphr_trailing_comma",
    "w_pphr_ws",
    "w_pphr_trailing_garbage",
    "r_pphr_syntax",
];

/// Standard error with the text of a syntax error replaced by the port's:
/// the line keeps the source's frame, `Error: Unable to read file (… in: )`,
/// and the Release text must be `nlohmann::json`'s parse error 101.
fn with_serde_wording(port_err: &str) -> impl Fn(String) -> String + '_ {
    move |release: String| {
        let frame = |text: &str| -> (usize, usize) {
            let start = text.find("Error: Unable to read file (").expect("the frame") + 28;
            let end = text.rfind(" in: )\n").expect("the frame's end");
            (start, end)
        };
        let (start, end) = frame(&release);
        assert!(
            release[start..end].starts_with("[json.exception.parse_error.101] parse error at line "),
            "{release}"
        );
        let (port_start, port_end) = frame(port_err);
        format!(
            "{}{}{}",
            &release[..start],
            &port_err[port_start..port_end],
            &release[end..]
        )
    }
}

/// Every `-write_ini` case of PeakPickerHiRes: the written INI shows each
/// value the load set, and each failure is the Release build's, exit code,
/// wording and all. That includes the conversions of the `nlohmann::json`
/// accessors (`1e30` to `-9223372036854775808`, `18446744073709551615` to
/// `-1`, `true` in an integer list to 1, `1e10` to `-2147483648`), `-0` as the
/// integer 0, the byte order mark, the names in byte order, the instance that
/// is always `1`, the file type that ignores case and a compression suffix, a
/// missing file, and a directory that fails in `basic_filebuf::underflow`.
#[test]
fn write_ini_from_a_json_ini_is_the_release_builds() {
    let names = cases("w_", "PeakPickerHiRes");
    assert_eq!(names.len(), 45);
    let mut syntax = 0;
    for case in &names {
        if OPEN.contains(&case.as_str()) {
            continue;
        }
        if SYNTAX_ERRORS.contains(&case.as_str()) {
            let replay = replay::<PeakPickerHiRes>(case);
            assert_case::<PeakPickerHiRes>(case, with_serde_wording(&replay.err));
            syntax += 1;
        } else {
            assert_case::<PeakPickerHiRes>(case, |text| text);
        }
    }
    assert_eq!(syntax, 4);
}

/// FuzzyDiff's string lists, integers, floats and two input files, and its
/// list type errors.
#[test]
fn write_ini_with_lists_is_the_release_builds() {
    let names = cases("w_", "FuzzyDiff");
    assert_eq!(names, ["w_fd_list_number", "w_fd_list_string", "w_fd_lists"]);
    for case in &names {
        assert_case::<FuzzyDiff>(case, |text| text);
    }
}

/// Runs that read a JSON `-ini`: the strict update after the load, an unknown
/// name in the run phase, with `-log` and `-debug` (`Assuming INI is a cwl
/// file` at level 1, the parameter dumps at level 2), a missing file and the
/// instance.
#[test]
fn a_run_with_a_json_ini_is_the_release_builds() {
    let names = cases("r_", "PeakPickerHiRes");
    assert_eq!(names.len(), 9);
    for case in &names {
        if OPEN.contains(&case.as_str()) {
            continue;
        }
        if SYNTAX_ERRORS.contains(&case.as_str()) {
            let replay = replay::<PeakPickerHiRes>(case);
            assert_case::<PeakPickerHiRes>(case, with_serde_wording(&replay.err));
        } else {
            assert_case::<PeakPickerHiRes>(case, |text| text);
        }
    }
}

/// FileInfo run from a JSON `-ini` that names its input as a file object:
/// the report, the command line winning over the INI, and the parameter dump
/// at debug level 2.
#[cfg(feature = "featurexml")]
#[test]
fn file_info_runs_from_a_json_ini_as_the_release_build_does() {
    let names = cases("r_fi_", "FileInfo");
    assert_eq!(names.len(), 5);
    for case in &names {
        assert_case::<openms::cli::tools::FileInfo>(case, |text| text);
    }
}

/// The retained Release runs are the ones the oracle made: every case of its
/// summary is here with the exit code it recorded.
#[test]
fn the_retained_runs_are_the_oracles() {
    let summary = read("", "summary.tsv").unwrap();
    let mut count = 0;
    for line in summary.lines().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            read(fields[0], "exit_code.txt").unwrap().trim(),
            fields[2],
            "{line}"
        );
        count += 1;
    }
    assert_eq!(count, 62);
    let host = read("", "host.txt").unwrap();
    assert!(host.starts_with("ibminode06\n"), "{host}");
}
