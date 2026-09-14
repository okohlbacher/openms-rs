// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The upstream `TOPP_BaselineFilter_1` test against the retained C++ output.

// The TOPP framework lives behind `paramxml` (every tool supports -ini) and these
// tools read mzML, so the whole file is inert without both features. Without this
// gate `cargo test --no-default-features` fails to compile.
#![cfg(all(feature = "mzml", feature = "paramxml"))]

use openms::cli::tools::BaselineFilter;
use openms::cli::{ExitCode, TEST_MODE_COMPLETION_TIME, run_with};
use openms::data_structures::DateTime;
use openms::format::file_handler::FileHandler;
use openms::format::file_types::FileType;
use openms::kernel::MSExperiment;
use openms::system::file::TempDir;
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

/// The processing records of every spectrum and chromatogram agree with the
/// retained C++ output: the same records in the same order, each with the same
/// software name and version, actions, completion time and metadata.
///
/// Decision D4: the source tool attaches `getProcessingInfo_` to its output
/// (`BaselineFilter.cpp`, `addDataProcessing_(ms_exp,
/// getProcessingInfo_(DataProcessing::BASELINE_REDUCTION))`), and the retained
/// `BaselineFilter_output.mzML` carries that record after the two the input
/// already had.
///
/// Completion times are compared to the minute, the precision the C++ mzML
/// writer keeps (`MzMLHandler.cpp:3947` writes `yyyy-MM-dd+hh:mm`): the record
/// the tool builds under `-test` says 23:59:59, and the retained file says
/// 23:59. This port's writer keeps the seconds.
fn assert_same_processing(produced: &MSExperiment, expected: &MSExperiment) {
    let minutes = |time: Option<DateTime>| time.map(|t| t.format("yyyy-MM-dd+hh:mm").unwrap());
    let records = |experiment: &MSExperiment| {
        experiment
            .spectra
            .iter()
            .map(|s| s.data_processing.clone())
            .chain(
                experiment
                    .chromatograms
                    .iter()
                    .map(|c| c.data_processing.clone()),
            )
            .collect::<Vec<_>>()
    };
    let (produced, expected) = (records(produced), records(expected));
    assert_eq!(produced.len(), expected.len(), "record holders");
    for (index, (a, e)) in produced.iter().zip(&expected).enumerate() {
        assert_eq!(a.len(), e.len(), "holder {index}: processing records");
        for (i, (pa, pe)) in a.iter().zip(e).enumerate() {
            let at = format!("holder {index}, record {i}");
            assert_eq!(pa.software.name, pe.software.name, "{at}: software");
            assert_eq!(pa.software.version, pe.software.version, "{at}: version");
            assert_eq!(pa.actions, pe.actions, "{at}: actions");
            assert_eq!(
                minutes(pa.completion_time),
                minutes(pe.completion_time),
                "{at}: completion time"
            );
            assert_eq!(pa.metadata, pe.metadata, "{at}: metadata");
        }
    }
}

#[test]
fn topp_baseline_filter_1_matches_the_retained_output() {
    let temp = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let out = temp.path().join("out.mzML");
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
    assert_same_processing(&produced, &reference);
    // The seconds the C++ writer drops are kept here: the -test record's time.
    for spectrum in &produced.spectra {
        assert_eq!(
            spectrum.data_processing.last().unwrap().completion_time,
            Some(DateTime::parse(TEST_MODE_COMPLETION_TIME).unwrap())
        );
    }
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
