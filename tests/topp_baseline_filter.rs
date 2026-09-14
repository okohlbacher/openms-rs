// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The upstream `TOPP_BaselineFilter_1` test against the retained C++ output.

// The TOPP framework lives behind `paramxml` (every tool supports -ini) and these
// tools read mzML, so the whole file is inert without both features. Without this
// gate `cargo test --no-default-features` fails to compile.
#![cfg(all(feature = "mzml", feature = "paramxml"))]

use openms::cli::tools::BaselineFilter;
use openms::cli::{ExitCode, run_with};
use openms::format::file_handler::FileHandler;
use openms::format::file_types::FileType;
use openms::kernel::MSExperiment;
use std::fs;
use std::path::{Path, PathBuf};

fn fixture(name: &str) -> PathBuf {
    Path::new("tests/data").join(name)
}
fn load(path: impl AsRef<Path>) -> MSExperiment {
    FileHandler::load_experiment(path, &[FileType::MzMl]).unwrap()
}
fn run(args: &[&str]) -> (ExitCode, String) {
    let arguments: Vec<String> = std::iter::once("BaselineFilter".to_string())
        .chain(args.iter().map(|a| (*a).to_string()))
        .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<BaselineFilter>(&arguments, &mut out, &mut err);
    (code, String::from_utf8_lossy(&err).into_owned())
}

#[test]
fn topp_baseline_filter_1_matches_the_retained_output() {
    let dir = std::env::temp_dir().join(format!("openms-bl-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let out = dir.join("out.mzML");
    let (code, err) = run(&[
        "-test",
        "-in",
        &fixture("baseline_filter_tool_input.mzML").to_string_lossy(),
        "-out",
        &out.to_string_lossy(),
        "-struc_elem_length",
        "1.5",
    ]);
    assert_eq!(code, ExitCode::ExecutionOk, "{err}");

    let produced = load(&out);
    let reference = load(fixture("baseline_filter_tool_output.mzML"));
    assert_eq!(produced.spectra.len(), reference.spectra.len());
    let mut compared = 0usize;
    for (index, (a, e)) in produced
        .spectra
        .iter()
        .zip(reference.spectra.iter())
        .enumerate()
    {
        assert_eq!(a.peaks.len(), e.peaks.len(), "spectrum {index} peak count");
        for (i, (pa, pe)) in a.peaks.iter().zip(e.peaks.iter()).enumerate() {
            assert!(
                (pa.mz - pe.mz).abs() <= 1e-9,
                "spectrum {index} peak {i} m/z"
            );
            assert!(
                (f64::from(pa.intensity) - f64::from(pe.intensity)).abs() <= 1e-4,
                "spectrum {index} peak {i} intensity"
            );
            compared += 1;
        }
    }
    assert!(compared > 100, "expected the fixture's full peak set");
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn unsorted_input_and_bad_parameters_are_refused() {
    // A registered restriction rejects an unknown method outright.
    assert_eq!(
        run(&[
            "-in",
            &fixture("baseline_filter_tool_input.mzML").to_string_lossy(),
            "-out",
            "x.mzML",
            "-method",
            "sharpen"
        ])
        .0,
        ExitCode::IllegalParameters
    );
    // A bare invocation is ILLEGAL_PARAMETERS (TOPPBase.cpp:227-232).
    assert_eq!(run(&[]).0, ExitCode::IllegalParameters);
}
