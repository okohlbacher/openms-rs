// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The tool executables on a console, against the C++ Release build: the
//! usage text shaped to the width `COLUMNS` or `stty size` reports
//! (`ConsoleUtils::readConsoleSize_`, `breakString_`, `IndentedStream`, core
//! `bc9cc12`; `TOPPBase::printUsage_`, cli `c19e494`).
//!
//! Evidence (tier 1, executed differential):
//! `../oracle/toppbase-completion/console.sh` ran the Release build of
//! `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576` on
//! ibminode06, each case in a clean environment (`env -i`, its own `HOME`,
//! `OPENMS_HOME_PATH` and working directory, standard input `/dev/null`). The
//! 115 case directories are retained under `tests/data/topp_cli_console` with
//! their hashes in `fixtures.sha256.json`. The executables here run the same
//! way (with none of the variables a tool reads but the case's own, their
//! own `HOME` and working directory, standard input `/dev/null`) and are
//! compared byte for byte on both streams and in their exit code;
//! only the Release run's paths are mapped to this run's, the closing
//! `<tool> took …` line is compared by shape, and the line `stty size` prints
//! when standard input is not a terminal is the local `stty`'s own (GNU's on
//! the oracle host, BSD's on macOS), taken from running it.
//!
//! The `tty_*` cases ran inside a pseudo-terminal, which a test cannot open
//! without unsafe code; their usage text and log lines are compared in the
//! crate's unit tests (`src/cli/usage.rs`, `src/cli/console.rs`) and the
//! executables in a pseudo-terminal by
//! `../oracle/toppbase-completion/compare_tty.py` (see
//! `docs/TOPP_CLI_SUPPORT.md`).
#![cfg(all(feature = "mzml", feature = "paramxml"))]

#[path = "support/took_line.rs"]
mod took_line;

use openms::system::file::TempDir;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// Where the Release run kept its case directories and inputs.
const ORACLE_RESULTS: &str = "/scratch/kohlbach/toppbase-oracle/run1/console_results";
const ORACLE_INPUTS: &str = "/scratch/kohlbach/toppbase-oracle/run1/inputs";
/// The GNU `stty` complaint the Release run recorded.
const GNU_STTY_LINE: &str = "stty: 'standard input': Inappropriate ioctl for device\n";

fn data(case: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/topp_cli_console")
        .join(case)
}

fn lines(case: &str, file: &str) -> Vec<String> {
    fs::read_to_string(data(case).join(file))
        .unwrap_or_else(|e| panic!("{case}/{file}: {e}"))
        .lines()
        .map(str::to_owned)
        .collect()
}

/// The executable of a ported tool, when this build has it.
fn executable(tool: &str) -> Option<&'static str> {
    match tool {
        "BaselineFilter" => Some(env!("CARGO_BIN_EXE_BaselineFilter")),
        "DTAExtractor" => Some(env!("CARGO_BIN_EXE_DTAExtractor")),
        "MapNormalizer" => Some(env!("CARGO_BIN_EXE_MapNormalizer")),
        "MzMLSplitter" => Some(env!("CARGO_BIN_EXE_MzMLSplitter")),
        "PeakPickerHiRes" => Some(env!("CARGO_BIN_EXE_PeakPickerHiRes")),
        "SpectraFilterWindowMower" => Some(env!("CARGO_BIN_EXE_SpectraFilterWindowMower")),
        #[cfg(feature = "featurexml")]
        "FileInfo" => Some(env!("CARGO_BIN_EXE_FileInfo")),
        #[cfg(feature = "featurexml")]
        "FeatureFinderCentroided" => Some(env!("CARGO_BIN_EXE_FeatureFinderCentroided")),
        _ => None,
    }
}

/// What `stty size` prints on standard error when standard input is
/// `/dev/null`, here, as the executable lets it through.
#[cfg(unix)]
fn local_stty_line() -> String {
    let output = Command::new("/bin/sh")
        .args(["-c", "stty size"])
        .stdin(Stdio::null())
        .output()
        .expect("/bin/sh runs");
    String::from_utf8(output.stderr).expect("UTF-8 stty message")
}

/// Off Unix the port has no `stty` to ask, so it prints nothing.
#[cfg(not(unix))]
fn local_stty_line() -> String {
    String::new()
}

struct Run {
    code: i32,
    out: String,
    err: String,
}

/// Run a retained case's command line on this port's executable, as the
/// oracle ran it; `None` when this build lacks the tool.
fn run_case(case: &str) -> Option<(Run, TempDir)> {
    let argv = lines(case, "argv.txt");
    let binary = executable(&argv[0])?;
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let home = dir.path().join("home");
    let cwd = dir.path().join("cwd");
    fs::create_dir_all(&home).unwrap();
    fs::create_dir_all(&cwd).unwrap();
    let arguments: Vec<String> = argv[1..]
        .iter()
        .map(|argument| map_paths(case, argument, &dir))
        .collect();
    // The oracle's `env -i` stands for "none of the variables a tool reads";
    // the rest of the environment stays, because the executable may need it
    // to start (a library path, for one).
    let mut command = Command::new(binary);
    command
        .args(&arguments)
        .env_remove("COLUMNS")
        .env_remove("OPENMS_TOOL_PREFIX_PATH")
        .env_remove("OPENMS_TTD_INTERNAL_PATH")
        .env_remove("OPENMS_DATA_PATH")
        .env("HOME", &home)
        .env("OPENMS_HOME_PATH", &home)
        .current_dir(&cwd)
        .stdin(Stdio::null());
    for variable in lines(case, "env.txt") {
        if variable.is_empty() {
            continue;
        }
        let (name, value) = variable.split_once('=').expect("NAME=value");
        command.env(name, value);
    }
    let output = command.output().expect("the executable runs");
    Some((
        Run {
            code: output.status.code().expect("an exit code"),
            out: String::from_utf8(output.stdout).expect("UTF-8 output"),
            err: String::from_utf8(output.stderr).expect("UTF-8 errors"),
        },
        dir,
    ))
}

/// The Release run's paths mapped to this run's: its working directory, and
/// the one input the cases name, `tests/data/baseline_filter_tool_input.mzML`.
fn map_paths(case: &str, text: &str, dir: &TempDir) -> String {
    let input =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/baseline_filter_tool_input.mzML");
    text.replace(
        &format!("{ORACLE_RESULTS}/{case}/cwd"),
        &dir.path().join("cwd").to_string_lossy(),
    )
    .replace(
        &format!("{ORACLE_INPUTS}/baseline_filter_tool_input.mzML"),
        &input.to_string_lossy(),
    )
}

/// The Release run's text with its paths mapped and its `stty` line
/// replaced by the local one.
fn expected(case: &str, file: &str, dir: &TempDir) -> String {
    let text =
        fs::read_to_string(data(case).join(file)).unwrap_or_else(|e| panic!("{case}/{file}: {e}"));
    map_paths(case, &text, dir).replace(GNU_STTY_LINE, &local_stty_line())
}

/// Compare one case on both streams and the exit code; the closing line is
/// required on standard output exactly where the Release run printed one.
fn assert_case(case: &str, transform: impl Fn(String) -> String) -> bool {
    let Some((run, dir)) = run_case(case) else {
        return false;
    };
    let want_code: i32 = lines(case, "exit_code.txt")[0].parse().unwrap();
    assert_eq!(run.code, want_code, "{case}: {}", run.err);
    assert_eq!(
        run.err,
        transform(expected(case, "stderr.txt", &dir)),
        "{case}: standard error"
    );
    let want_out = expected(case, "stdout.txt", &dir);
    let (want_body, want_took) = took_line::split_took_line(&lines(case, "argv.txt")[0], &want_out);
    let (got_body, got_took) = took_line::split_took_line(&lines(case, "argv.txt")[0], &run.out);
    assert_eq!(got_body, want_body, "{case}: standard output");
    assert_eq!(
        got_took.is_some(),
        want_took.is_some(),
        "{case}: closing line"
    );
    true
}

/// The retained cases whose names start with `prefix`.
fn cases(prefix: &str) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(data(""))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|name| name.starts_with(prefix))
        .collect();
    names.sort();
    names
}

/// `--help` and `--helphelp` of every tool at `COLUMNS` 20, 28, 40, 60, 80
/// and 120: words carried to the next line, option descriptions continued
/// under their column, items longer than ten lines cut with `...`, the
/// subsection list broken twice, and a width narrower than the option column,
/// where the source's unsigned arithmetic stops shortening later lines.
#[test]
fn usage_is_shaped_to_columns_as_the_release_build_shapes_it() {
    let mut compared = 0;
    let mut expected_count = 0;
    for case in cases("help") {
        let tool = lines(&case, "argv.txt")[0].clone();
        if tool == "FeatureFinderCentroided" {
            continue;
        }
        if executable(&tool).is_some() {
            expected_count += 1;
        }
        if assert_case(&case, |text| text) {
            compared += 1;
        }
    }
    assert_eq!(compared, expected_count);
    let tools = if cfg!(feature = "featurexml") { 7 } else { 6 };
    assert_eq!(compared, tools * 10);
}

/// `COLUMNS` as `StringUtils::toInt32` reads it (surrounding space and one
/// `+` allowed, nothing else), less one, and no shaping below 10: `abc`,
/// empty, `0`, `9`, `10`, `-5`, `45x`, `4 5` and out-of-range values leave the
/// text unshaped; `11`, `12`, ` 45 `, `+45` and `2147483647` shape it. Unset,
/// the executable asks `stty size`, whose complaint about a standard input
/// that is not a terminal reaches standard error first.
#[test]
fn columns_values_are_read_as_the_release_build_reads_them() {
    let names = cases("columns_");
    assert_eq!(names.len(), 16);
    for case in &names {
        assert!(assert_case(case, |text| text), "{case}");
    }
    assert!(
        fs::read_to_string(data("columns_unset").join("stderr.txt"))
            .unwrap()
            .starts_with(GNU_STTY_LINE)
    );
}

/// The usage text after a failure is shaped as `--help`'s is: a missing
/// required option, an unknown option, and the usage `FileInfo` prints
/// itself when `-i` meets a file that is not mzML.
#[test]
fn usage_after_a_failure_is_shaped_as_the_release_build_shapes_it() {
    let mut compared = 0;
    for case in [
        "missing_in_c50",
        "unknown_opt_c50",
        "fileinfo_i_dta_c50",
        "fileinfo_i_dta_c50_helphelp",
    ] {
        if assert_case(case, |text| text) {
            compared += 1;
        }
    }
    assert_eq!(compared, if cfg!(feature = "featurexml") { 4 } else { 2 });
}

/// `FeatureFinderCentroided` does not yet declare the two citations the
/// source's tool registers (`FeatureFinderCentroided.cpp:124-136`; its tool
/// lives in `src/cli/tools`), so its usage lacks the `To cite
/// FeatureFinderCentroided:` block. Everything else is shaped as the Release
/// build shapes it: the Release text from that heading to `Usage:` is
/// replaced by what the source writes there without citations, and the
/// heading is required, so that this test fails once the citations are
/// declared.
#[cfg(feature = "featurexml")]
#[test]
fn feature_finder_centroided_differs_only_by_its_undeclared_citations() {
    let names: Vec<String> = cases("help")
        .into_iter()
        .filter(|case| case.contains("FeatureFinderCentroided"))
        .collect();
    assert_eq!(names.len(), 10);
    for case in &names {
        // Without citations the source writes the one empty line
        // (`is << is.indent(0) << "\n"`) between the OpenMS citation and
        // `Usage:`. With them, the block and what its last line leaves behind:
        // a last line that fills the width exactly makes the next line break
        // an extra empty line (`helphelp_…_c28`).
        let without_citations = |text: String| {
            // The narrowest widths break the heading itself.
            let start = text
                .find("To cite Feature")
                .unwrap_or_else(|| panic!("{case}: the Release text cites the tool"));
            let usage = text.find("Usage:").expect("the usage heading");
            assert!(start < usage, "{case}");
            format!("{}\n{}", &text[..start], &text[usage..])
        };
        assert!(assert_case(case, without_citations), "{case}");
    }
}
