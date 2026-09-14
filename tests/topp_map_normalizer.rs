// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The upstream `TOPP_MapNormalizer_1` test, compared canonically against the
//! retained C++ output. The upstream `${DIFF}` runs with an index whitelist, so
//! byte equality was never the contract; this compares decoded content.

// The TOPP framework lives behind `paramxml` (every tool supports -ini) and these
// tools read mzML, so the whole file is inert without both features. Without this
// gate `cargo test --no-default-features` fails to compile.
#![cfg(all(feature = "mzml", feature = "paramxml"))]

use openms::cli::tools::MapNormalizer;
use openms::cli::{ExitCode, TEST_MODE_COMPLETION_TIME, run_with};
use openms::data_structures::DateTime;
use openms::format::file_handler::FileHandler;
use openms::format::file_types::FileType;
use openms::kernel::MSExperiment;
use openms::system::file::TempDir;
use std::path::{Path, PathBuf};

/// The processing records of every spectrum and chromatogram agree with the
/// retained C++ output: the same records in the same order, each with the same
/// software name and version, actions, completion time and metadata.
///
/// Decision D4: the source tool attaches `getProcessingInfo_` to its output
/// (`MapNormalizer.cpp`, `addDataProcessing_(exp,
/// getProcessingInfo_(DataProcessing::NORMALIZATION))`), and the retained
/// `MapNormalizer_output.mzML` carries that `intensity normalization` record
/// after the two the input already had.
///
/// Completion times are compared to the minute, the precision the C++ mzML
/// writer keeps (`MzMLHandler.cpp:3947` writes `yyyy-MM-dd+hh:mm`).
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

fn fixture(name: &str) -> PathBuf {
    Path::new("tests/data").join(name)
}

fn load(path: impl AsRef<Path>) -> MSExperiment {
    FileHandler::load_experiment(path, &[FileType::MzMl]).unwrap()
}

fn run(args: &[&str]) -> (ExitCode, String) {
    let arguments: Vec<String> = std::iter::once("MapNormalizer".to_string())
        .chain(args.iter().map(|a| (*a).to_string()))
        .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<MapNormalizer>(&arguments, &mut out, &mut err);
    (code, String::from_utf8_lossy(&err).into_owned())
}

#[test]
fn topp_map_normalizer_1_scales_ms1_to_percent_of_the_run_maximum() {
    let temp = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let out = temp.path().join("out.mzML");

    let (code, err) = run(&[
        "-test",
        "-in",
        &fixture("map_normalizer_input.mzML").to_string_lossy(),
        "-out",
        &out.to_string_lossy(),
    ]);
    assert_eq!(code, ExitCode::ExecutionOk, "{err}");

    let produced = load(&out);
    let reference = load(fixture("map_normalizer_output.mzML"));
    assert_eq!(produced.spectra.len(), reference.spectra.len());

    let mut scaled = 0usize;
    for (index, (a, e)) in produced
        .spectra
        .iter()
        .zip(reference.spectra.iter())
        .enumerate()
    {
        assert_eq!(a.ms_level, e.ms_level, "spectrum {index} MS level");
        assert_eq!(a.peaks.len(), e.peaks.len(), "spectrum {index} peak count");
        for (i, (pa, pe)) in a.peaks.iter().zip(e.peaks.iter()).enumerate() {
            assert!(
                (pa.mz - pe.mz).abs() <= 1e-9,
                "spectrum {index} peak {i} m/z"
            );
            // Intensities are stored as f32 after an f64 division, so compare
            // relatively rather than exactly.
            let (x, y) = (f64::from(pa.intensity), f64::from(pe.intensity));
            assert!(
                (x - y).abs() <= 1e-3 * y.abs().max(1.0),
                "spectrum {index} peak {i} intensity {x} vs {y}"
            );
        }
        if a.ms_level < 2 {
            scaled += 1;
        }
    }
    assert!(scaled > 0, "the fixture must contain MS1 spectra");

    // The most intense MS1 peak becomes 100, which is what "normalize" means here.
    let peak = produced
        .spectra
        .iter()
        .filter(|s| s.ms_level < 2)
        .flat_map(|s| s.peaks.iter())
        .fold(0.0f32, |acc, p| acc.max(p.intensity));
    assert!(
        (f64::from(peak) - 100.0).abs() < 1e-3,
        "MS1 maximum should normalise to 100, got {peak}"
    );

    assert_same_processing(&produced, &reference);
    for spectrum in &produced.spectra {
        assert_eq!(
            spectrum.data_processing.last().unwrap().completion_time,
            Some(DateTime::parse(TEST_MODE_COMPLETION_TIME).unwrap())
        );
    }
}

#[test]
fn required_parameters_and_formats_are_enforced() {
    // A bare invocation is ILLEGAL_PARAMETERS (TOPPBase.cpp:227-232).
    assert_eq!(run(&[]).0, ExitCode::IllegalParameters);
    assert_eq!(
        run(&["-in", "absent.mzML", "-out", "x.mzML"]).0,
        ExitCode::InputFileNotFound
    );
    // 'out' is registered as mzML only.
    assert_eq!(
        run(&[
            "-in",
            &fixture("map_normalizer_input.mzML").to_string_lossy(),
            "-out",
            "x.dta"
        ])
        .0,
        ExitCode::IllegalParameters
    );
}
