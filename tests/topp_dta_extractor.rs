// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The three upstream `TOPP_DTAExtractor_*` tests, reproduced against the
//! retained C++ reference outputs. These fixtures were produced by the C++
//! tool, so agreement here is executed differential evidence for the whole
//! chain: command line, parameter validation, mzML reading and DTA writing.

// The TOPP framework lives behind `paramxml` (every tool supports -ini) and these
// tools read mzML, so the whole file is inert without both features. Without this
// gate `cargo test --no-default-features` fails to compile.
#![cfg(all(feature = "mzml", feature = "paramxml"))]

use openms::cli::tools::DTAExtractor;
use openms::cli::{ExitCode, run_with};
use std::fs;
use std::path::{Path, PathBuf};

fn fixture(name: &str) -> PathBuf {
    Path::new("tests/data").join(name)
}

/// Run the tool into a fresh directory and compare one produced file with its
/// retained C++ output, byte for byte.
fn check(case: &str, args: &[&str], produced: &str, expected: &str) {
    // Cases share produced file names, so the directory is keyed by the case.
    let dir = std::env::temp_dir().join(format!("openms-dta-{}-{case}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let base = dir.join("DTAExtractor");

    let mut arguments = vec![
        "DTAExtractor".to_string(),
        "-test".to_string(),
        "-in".to_string(),
        fixture("dta_extractor_input.mzML")
            .to_string_lossy()
            .into_owned(),
        "-out".to_string(),
        base.to_string_lossy().into_owned(),
    ];
    arguments.extend(args.iter().map(|a| (*a).to_string()));

    let mut out = Vec::new();
    let mut err = Vec::new();
    let code = run_with::<DTAExtractor>(&arguments, &mut out, &mut err);
    assert_eq!(
        code,
        ExitCode::ExecutionOk,
        "tool failed: {}",
        String::from_utf8_lossy(&err)
    );

    let actual = fs::read(dir.join(produced))
        .unwrap_or_else(|_| panic!("{produced} was not written; produced: {:?}", listing(&dir)));
    let reference = fs::read(fixture(expected)).unwrap();
    assert_eq!(
        String::from_utf8_lossy(&actual),
        String::from_utf8_lossy(&reference),
        "{produced} differs from the C++ output {expected}"
    );
    let _ = fs::remove_dir_all(&dir);
}

fn listing(dir: &Path) -> Vec<String> {
    fs::read_dir(dir)
        .map(|entries| {
            entries
                .filter_map(|e| e.ok())
                .map(|e| e.file_name().to_string_lossy().into_owned())
                .collect()
        })
        .unwrap_or_default()
}

#[test]
fn topp_dta_extractor_1_retention_time_range() {
    check(
        "1",
        &["-rt", ":61"],
        "DTAExtractor_RT60.0.dta",
        "dta_extractor_1_output.dta",
    );
}

#[test]
fn topp_dta_extractor_2_ms_level_one() {
    check(
        "2",
        &["-level", "1"],
        "DTAExtractor_RT60.0.dta",
        "dta_extractor_2_output.dta",
    );
}

#[test]
fn topp_dta_extractor_3_precursor_mz_range() {
    check(
        "3",
        &["-level", "2", "-mz", ":1000"],
        "DTAExtractor_RT140.0_MZ5.0.dta",
        "dta_extractor_3_output.dta",
    );
}

#[test]
fn usage_and_exit_codes_follow_the_source_contract() {
    let run = |args: &[&str]| -> (ExitCode, String, String) {
        let arguments: Vec<String> = std::iter::once("DTAExtractor".to_string())
            .chain(args.iter().map(|a| (*a).to_string()))
            .collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = run_with::<DTAExtractor>(&arguments, &mut out, &mut err);
        (
            code,
            String::from_utf8_lossy(&out).into_owned(),
            String::from_utf8_lossy(&err).into_owned(),
        )
    };

    // Usage short-circuits before any validation, so required values are not demanded.
    let (code, out, _) = run(&["--help"]);
    assert_eq!(code, ExitCode::ExecutionOk);
    assert!(out.contains("DTAExtractor --"), "{out}");
    assert!(out.contains("-in <file>*"), "{out}");
    // Advanced parameters are hidden until asked for.
    assert!(!out.contains("-no_progress"), "{out}");
    assert!(run(&["--helphelp"]).1.contains("-no_progress"));

    // A bare invocation is ILLEGAL_PARAMETERS (TOPPBase.cpp:227-232; oracle
    // no_arguments in ../oracle/topp-cli-lifecycle).
    assert_eq!(run(&[]).0, ExitCode::IllegalParameters);

    // A missing input file is INPUT_FILE_NOT_FOUND.
    let (code, _, err) = run(&["-in", "absent.mzML", "-out", "x"]);
    assert_eq!(code, ExitCode::InputFileNotFound);
    assert!(err.contains("does not exist"), "{err}");

    // A registered format is enforced on the input path.
    let (code, _, err) = run(&[
        "-in",
        &fixture("dta_extractor_1_output.dta").to_string_lossy(),
        "-out",
        "x",
    ]);
    assert_eq!(code, ExitCode::IllegalParameters);
    assert!(err.contains("unsupported format"), "{err}");
}

#[test]
fn write_ini_round_trips_through_the_parameter_file() {
    let dir = std::env::temp_dir().join(format!("openms-dta-ini-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let ini = dir.join("DTAExtractor.ini");

    let arguments: Vec<String> = ["DTAExtractor", "-write_ini", &ini.to_string_lossy()]
        .iter()
        .map(|a| (*a).to_string())
        .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    assert_eq!(
        run_with::<DTAExtractor>(&arguments, &mut out, &mut err),
        ExitCode::ExecutionOk,
        "{}",
        String::from_utf8_lossy(&err)
    );

    // The written file is a loadable parameter tree holding the tool defaults,
    // and does not carry the write_ini request itself.
    let loaded = openms::format::paramxml::load(&ini).unwrap();
    assert_eq!(
        loaded
            .value("DTAExtractor:1:level")
            .unwrap()
            .as_str()
            .unwrap(),
        "1,2,3"
    );
    assert!(!loaded.exists("DTAExtractor:1:write_ini").unwrap());

    // Feeding it back supplies the same defaults, and the command line still wins.
    let arguments: Vec<String> = [
        "DTAExtractor",
        "-ini",
        &ini.to_string_lossy(),
        "-in",
        &fixture("dta_extractor_input.mzML").to_string_lossy(),
        "-out",
        &dir.join("out").to_string_lossy(),
        "-level",
        "1",
    ]
    .iter()
    .map(|a| (*a).to_string())
    .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    assert_eq!(
        run_with::<DTAExtractor>(&arguments, &mut out, &mut err),
        ExitCode::ExecutionOk,
        "{}",
        String::from_utf8_lossy(&err)
    );
    assert!(dir.join("out_RT60.0.dta").exists());
    let _ = fs::remove_dir_all(&dir);
}
