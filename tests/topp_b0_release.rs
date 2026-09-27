// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The TOPP cases of the B0 integration lane's Release oracle
//! (`../oracle/B0-attacks`, 2026-09-27), retained under
//! `tests/data/b0_release/topp` with provenance in
//! `tests/data/b0_release_provenance.json`: the Linux x86_64 C++ Release
//! build at the pins (core `bc9cc12`, cli `c19e494`, topp `174b576`) on
//! ibminode06, each case in a clean environment, run twice with identical
//! results once paths and run times are masked (`bf_dp_negative`'s standard
//! error is the shell's report of the crash, which names a process id).
//!
//! Each case is replayed in process with the Release run's paths mapped to the
//! replay's, and compared with what the Release build did.
#![cfg(all(feature = "mzml", feature = "featurexml", feature = "paramxml"))]

use openms::cli::tools::{
    BaselineFilter, DTAExtractor, FeatureFinderCentroided, MapNormalizer, MzMLSplitter,
    PeakPickerHiRes, SpectraFilterWindowMower,
};
use openms::cli::{ExitCode, Tool, run_with};
use openms::system::file::TempDir;
use std::fs;
use std::path::{Path, PathBuf};

/// Where the Release run kept its inputs and its per-case working directories.
const RELEASE_INPUTS: &str = "/scratch/kohlbach/b0-oracle/inputs";
const RELEASE_RUN: &str = "/scratch/kohlbach/b0-oracle/run1";

fn release(case: &str, file: &str) -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/b0_release/topp")
        .join(case)
        .join(file);
    fs::read_to_string(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

struct Replay {
    code: ExitCode,
    out: String,
    err: String,
    /// The Release run's text with its paths replaced by the replay's.
    map: Box<dyn Fn(&str) -> String>,
    _work: TempDir,
}

/// Run `T` with the Release case's arguments, the inputs staged in a fresh
/// directory under the name the Release run gave them.
fn replay<T: Tool>(case: &str, inputs: &[(&str, PathBuf)]) -> Replay {
    let work = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let (input_dir, cwd) = (work.path().join("inputs"), work.path().join("cwd"));
    fs::create_dir_all(&input_dir).unwrap();
    fs::create_dir_all(&cwd).unwrap();
    for (name, source) in inputs {
        fs::copy(source, input_dir.join(name)).unwrap();
    }
    let (input_text, cwd_text) = (
        input_dir.to_str().unwrap().to_owned(),
        cwd.to_str().unwrap().to_owned(),
    );
    let release_cwd = format!("{RELEASE_RUN}/{case}/cwd");
    let map = move |text: &str| {
        text.replace(RELEASE_INPUTS, &input_text)
            .replace(&release_cwd, &cwd_text)
    };
    let arguments: Vec<String> = release(case, "argv.txt").lines().map(&map).collect();
    assert_eq!(arguments[0], T::NAME, "{case}");
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<T>(&arguments, &mut out, &mut err);
    Replay {
        code,
        out: String::from_utf8(out).unwrap(),
        err: String::from_utf8(err).unwrap(),
        map: Box::new(map),
        _work: work,
    }
}

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn unknown_input() -> Vec<(&'static str, PathBuf)> {
    vec![(
        "unknown_name_and_content",
        data("topp_cli_completion/unknown_name_and_content"),
    )]
}

/// An input whose name and content no type claims, for every ported tool that
/// loads with an allowed list: the framework's format check only warns
/// (`TOPPBase.cpp:1580-1583`), and the load refuses the type with
/// `ParseError` (`FileHandler.cpp:859-863`), exit 3, before anything is
/// written. The exit code, standard output and the warning line are the
/// Release build's; the error line is the Release build's text inside the
/// `parse error on line 0: ` with which the lifecycle renders every
/// `Error::Parse` (open item 3 of the TOPPBase ledger entry).
#[test]
fn an_undetermined_input_format_is_the_loaders_parse_error_in_every_tool() {
    /// A retained case and the tool that replays it.
    type Case = (&'static str, fn(&str) -> Replay);
    let cases: [Case; 7] = [
        ("unknown_bf", |case| {
            replay::<BaselineFilter>(case, &unknown_input())
        }),
        ("unknown_pphr", |case| {
            replay::<PeakPickerHiRes>(case, &unknown_input())
        }),
        ("unknown_ffc", |case| {
            replay::<FeatureFinderCentroided>(case, &unknown_input())
        }),
        ("unknown_dta", |case| {
            replay::<DTAExtractor>(case, &unknown_input())
        }),
        ("unknown_sfwm", |case| {
            replay::<SpectraFilterWindowMower>(case, &unknown_input())
        }),
        ("unknown_ms", |case| {
            replay::<MzMLSplitter>(case, &unknown_input())
        }),
        ("unknown_mn", |case| {
            replay::<MapNormalizer>(case, &unknown_input())
        }),
    ];
    for (case, run) in cases {
        let replay = run(case);
        let exit: i32 = release(case, "exit_code.txt").trim().parse().unwrap();
        assert_eq!(exit, 3, "{case}: the Release build's exit code");
        assert_eq!(replay.code.as_i32(), exit, "{case}: {}", replay.err);
        assert_eq!(
            replay.out,
            (replay.map)(&release(case, "stdout.txt")),
            "{case}"
        );
        let expected = (replay.map)(&release(case, "stderr.txt")).replacen(
            "Error: Unable to read file (",
            "Error: Unable to read file (parse error on line 0: ",
            1,
        );
        assert_eq!(replay.err, expected, "{case}");
        assert!(
            !replay
                ._work
                .path()
                .join("cwd")
                .read_dir()
                .unwrap()
                .any(|_| true),
            "{case}: nothing is written"
        );
    }
}

/// BaselineFilter with a `DataPoints` structuring element below one sample.
/// The source's `MorphologicalFilter` has no restriction on
/// `struc_elem_length` and converts it with `(UInt)(double)`
/// (`MorphologicalFilter.h:139`, `:177`, `:281`): a length in (-1, 1) becomes
/// a zero-sample element, which the source widens to one sample, and the
/// Release run completes (exit 0); a negative length is undefined behaviour
/// and the Release build ends with SIGSEGV (exit 139). This port refuses both
/// as illegal parameters, exit 6, and writes nothing
/// (`docs/MORPHOLOGICAL_FILTER_SUPPORT.md`, `docs/TOPP_CLI_SUPPORT.md`,
/// *Errors*).
#[test]
fn a_structuring_element_below_one_sample_is_refused_as_an_illegal_parameter() {
    let input = || {
        vec![(
            "baseline_filter_tool_input.mzML",
            data("baseline_filter_tool_input.mzML"),
        )]
    };
    for (case, release_exit, message) in [
        (
            "bf_dp_negative",
            139,
            "Invalid parameter: invalid value: struc_elem_length must not be negative",
        ),
        (
            "bf_dp_negative_fraction",
            0,
            "Invalid parameter: invalid value: structuring element must contain at least one point",
        ),
    ] {
        let exit: i32 = release(case, "exit_code.txt").trim().parse().unwrap();
        assert_eq!(exit, release_exit, "{case}: the Release build's exit code");
        let replay = replay::<BaselineFilter>(case, &input());
        assert_eq!(
            replay.code,
            ExitCode::IllegalParameters,
            "{case}: {}",
            replay.err
        );
        assert!(
            replay.err.lines().any(|line| line == message),
            "{case}: {}",
            replay.err
        );
        assert!(
            !replay._work.path().join("cwd/out.mzML").exists(),
            "{case}: nothing is written"
        );
    }
}
