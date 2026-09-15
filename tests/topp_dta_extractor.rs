// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The three upstream `TOPP_DTAExtractor_*` tests plus a real-instrument slice,
//! reproduced against the retained C++ reference outputs. Every `.dta` fixture
//! here was written by the pinned C++ `DTAExtractor` itself, so agreement is
//! executed differential evidence for the whole chain: command line, parameter
//! validation, mzML reading, DTA file naming and DTA writing.
//!
//! # Provenance
//!
//! Produced on 2026-09-15 on `ibminode06` with the Release build at
//! `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`
//! (core `bc9cc12`, cli `c19e494`, topp `174b576`):
//!
//! - `dta_extractor_{1,2,3}_output.dta`: `DTAExtractor -test -in
//!   dta_extractor_input.mzML -out <dir>/DTAExtractor` with the same flags the
//!   cases below pass. They replace the upstream `TOPP_DTAExtractor_*_output`
//!   files, which predate the C++ move from Boost.Karma to `std::to_chars`
//!   (`NumericFormatting.h`) and spell the peak m/z `120` where the pinned tool
//!   writes `120.0`. The upstream TOPP suite compares with `FuzzyDiff`, so the
//!   stale spelling still passes there; a byte comparison needs the real bytes.
//! - `dta_extractor_velos_slice.mzML`: `FileFilter -in
//!   /ceph/ibmi/abi/oliver/bench/openms4/inputs/derived/sub_centroid_velos_50amol_r1_first4000.mzML
//!   -rt 2.0:4.0 -mz 350:450`, three spectra of the LTQ Orbitrap Velos
//!   centroided run the TOPP benchmark uses.
//! - `dta_extractor_velos_*_output.dta`: the pinned `DTAExtractor` on that
//!   slice. These carry the numbers the round ones cannot: 15-fraction-digit
//!   m/z (`350.127610656737602`), 15-significant-digit intensities
//!   (`583.898498535156`), zero intensities that print as `0` rather than `0.0`
//!   because they take the other formatter, and a 15-significant-digit
//!   precursor header (`1015.72283935547 2`).

// The TOPP framework lives behind `paramxml` (every tool supports -ini) and these
// tools read mzML, so the whole file is inert without both features. Without this
// gate `cargo test --no-default-features` fails to compile.
#![cfg(all(feature = "mzml", feature = "paramxml"))]

use openms::cli::tools::DTAExtractor;
use openms::cli::{ExitCode, run_with};
use openms::system::file::TempDir;
use std::fs;
use std::path::{Path, PathBuf};

fn fixture(name: &str) -> PathBuf {
    Path::new("tests/data").join(name)
}

/// Run the tool into a fresh directory and compare one produced file with its
/// retained C++ output, byte for byte.
fn check(case: &str, args: &[&str], produced: &str, expected: &str) {
    // Cases share produced file names, so each case writes into its own
    // uniquely named directory.
    let temp = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let dir = temp.path();
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
        .unwrap_or_else(|_| panic!("{produced} was not written; produced: {:?}", listing(dir)));
    let reference = fs::read(fixture(expected)).unwrap();
    assert_eq!(
        String::from_utf8_lossy(&actual),
        String::from_utf8_lossy(&reference),
        "{produced} differs from the C++ output {expected} (case {case})"
    );
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

/// The real-instrument slice: every produced file name and every produced byte
/// against the pinned C++ tool's own output.
///
/// This is the case that pins the numeric text. The round values of the three
/// upstream cases agree under several formatting rules; these do not. Both the
/// file names and the peak m/z go through `StringUtils::toStr(double)`, and the
/// precursor header and the intensities through the stream's default float
/// field at precision 15, so a single wrong rule shows up here immediately.
#[test]
fn topp_dta_extractor_real_velos_slice_matches_cpp_bytes() {
    let temp = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let dir = temp.path();

    let arguments: Vec<String> = [
        "DTAExtractor",
        "-test",
        "-in",
        &fixture("dta_extractor_velos_slice.mzML").to_string_lossy(),
        "-out",
        &dir.join("DTAExtractor").to_string_lossy(),
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

    // The names the C++ tool produced, in the order the spectra appear. An
    // MS1 name carries the retention time alone; an MS2 name also carries the
    // precursor m/z, and both are `StringUtils::toStr(double)` text.
    let expected: [(&str, &str); 3] = [
        (
            "DTAExtractor_RT2.5674_MZ508.361419677733977.dta",
            "dta_extractor_velos_ms2_output.dta",
        ),
        (
            "DTAExtractor_RT2.82019999998.dta",
            "dta_extractor_velos_ms1_first_output.dta",
        ),
        (
            "DTAExtractor_RT3.53560000002.dta",
            "dta_extractor_velos_ms1_second_output.dta",
        ),
    ];

    let mut produced = listing(dir);
    produced.sort();
    let mut wanted: Vec<String> = expected
        .iter()
        .map(|(name, _)| (*name).to_owned())
        .collect();
    wanted.sort();
    assert_eq!(
        produced, wanted,
        "produced file names differ from the C++ run"
    );

    for (name, reference) in expected {
        let actual = fs::read(dir.join(name)).unwrap();
        let expected_bytes = fs::read(fixture(reference)).unwrap();
        // Compared as text so that a mismatch reports the differing line
        // rather than a byte offset.
        assert_eq!(
            String::from_utf8_lossy(&actual),
            String::from_utf8_lossy(&expected_bytes),
            "{name} differs from the C++ output {reference}"
        );
    }
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

    // Usage short-circuits before any validation, so required values are not
    // demanded. The source prints usage to standard error, for --help too
    // (TOPPBase.cpp:637; oracle help in ../oracle/topp-cli-lifecycle, whose
    // stdout is empty), so it is asserted on the error stream.
    let (code, out, err) = run(&["--help"]);
    assert_eq!(code, ExitCode::ExecutionOk);
    assert!(out.is_empty(), "{out}");
    assert!(err.contains("DTAExtractor --"), "{err}");
    assert!(err.contains("-in <file>*"), "{err}");
    // Advanced parameters are hidden until asked for.
    assert!(!err.contains("-no_progress"), "{err}");
    assert!(run(&["--helphelp"]).2.contains("-no_progress"));

    // A bare invocation is ILLEGAL_PARAMETERS (TOPPBase.cpp:227-232; oracle
    // no_arguments in ../oracle/topp-cli-lifecycle).
    assert_eq!(run(&[]).0, ExitCode::IllegalParameters);

    // A missing input file is INPUT_FILE_NOT_FOUND.
    let (code, _, err) = run(&["-in", "absent.mzML", "-out", "x"]);
    assert_eq!(code, ExitCode::InputFileNotFound);
    assert!(err.contains("does not exist"), "{err}");

    // A registered format is enforced on the input path, with the source's
    // InvalidParameter wording (TOPPBase.cpp:1584-1591; the same message for a
    // .dta input of SpectraFilterWindowMower in oracle in_dta_extension,
    // ../oracle/topp-cli-lifecycle/cli2/manifest.json).
    let dta = fixture("dta_extractor_1_output.dta")
        .to_string_lossy()
        .into_owned();
    let (code, _, err) = run(&["-in", &dta, "-out", "x"]);
    assert_eq!(code, ExitCode::IllegalParameters);
    assert!(
        err.contains(&format!(
            "Invalid parameter: Input file '{dta}' has invalid format 'dta'. Valid formats are: 'mzML'."
        )),
        "{err}"
    );
}

#[test]
fn write_ini_round_trips_through_the_parameter_file() {
    let temp = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let dir = temp.path();
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
}
