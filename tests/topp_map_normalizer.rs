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

    // The most intense MS1 peak becomes 100 *for this fixture*, which holds only
    // because it carries no chromatogram and its most intense peak is an MS1 one.
    // The general rule is `combined maximum / 100`; see
    // `a_chromatogram_above_the_spectrum_maximum_sets_the_scale`.
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

/// Run the tool on `input` and compare every decoded intensity — spectra and
/// chromatograms alike — against the retained C++ output.
///
/// The C++ side was executed at the pinned Release build
/// (`openms4-release-bc9cc12-c19e494-174b576`); the command and the digests are
/// in `tests/data/topp_map_normalizer_provenance.json`. Both sides divide an
/// `f64` by the same `f64` scale and narrow the quotient to `f32`, so the
/// intensities must agree bit for bit — a tolerance here would hide exactly the
/// kind of scale error this case exists to catch.
///
/// Every other array is compared exactly for the same reason. The source only
/// ever calls `pk.setIntensity(...)`, so m/z and the chromatogram time array
/// must survive the round trip untouched; the executed evidence shows they do,
/// bitwise, on both fixtures and on all 88,434,492 m/z points and 43,745 time
/// points of the full 1.2 GB run. A tolerance on them would only hide a future
/// regression that rewrites them.
fn assert_matches_cpp(input: &str, expected: &str) {
    let temp = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let out = temp.path().join("out.mzML");
    let (code, err) = run(&[
        "-test",
        "-in",
        &fixture(input).to_string_lossy(),
        "-out",
        &out.to_string_lossy(),
    ]);
    assert_eq!(code, ExitCode::ExecutionOk, "{err}");

    let produced = load(&out);
    let reference = load(fixture(expected));
    assert_eq!(produced.spectra.len(), reference.spectra.len(), "spectra");
    assert_eq!(
        produced.chromatograms.len(),
        reference.chromatograms.len(),
        "chromatograms"
    );
    for (index, (a, e)) in produced
        .spectra
        .iter()
        .zip(reference.spectra.iter())
        .enumerate()
    {
        assert_eq!(a.ms_level, e.ms_level, "spectrum {index} MS level");
        assert_eq!(a.peaks.len(), e.peaks.len(), "spectrum {index} peak count");
        for (i, (pa, pe)) in a.peaks.iter().zip(e.peaks.iter()).enumerate() {
            assert_eq!(pa.mz, pe.mz, "spectrum {index} peak {i} m/z");
            assert_eq!(
                pa.intensity, pe.intensity,
                "spectrum {index} peak {i} intensity"
            );
        }
    }
    // The source never touches chromatogram intensities: its chromatogram
    // branch is commented out. They still set the scale, so they are the one
    // place a wrong scale could also be written back.
    for (index, (a, e)) in produced
        .chromatograms
        .iter()
        .zip(reference.chromatograms.iter())
        .enumerate()
    {
        assert_eq!(a.peaks.len(), e.peaks.len(), "chromatogram {index} points");
        for (i, (pa, pe)) in a.peaks.iter().zip(e.peaks.iter()).enumerate() {
            assert_eq!(pa.rt, pe.rt, "chromatogram {index} point {i} RT");
            assert_eq!(
                pa.intensity, pe.intensity,
                "chromatogram {index} point {i} intensity"
            );
        }
    }
    assert_same_processing(&produced, &reference);
}

/// The scale is `getMaxIntensity() / 100` over the *combined* range, so a
/// chromatogram more intense than every peak sets it on its own.
///
/// This is the shape the 1.2 GB Velos benchmark run has: its TIC chromatogram
/// peaks at 1,788,496,256 against a spectrum maximum of 135,038,560, and taking
/// the spectrum maximum alone scaled every MS1 peak 13.2443x too high. The
/// fixture reproduces it in miniature — chromatogram maximum 654,321 against a
/// spectrum maximum of 50,000 — so the MS1 peak of 3,000 must come out as
/// 3000 / 6543.21 = 0.4584906, not 3000 / 500 = 6.
#[test]
fn a_chromatogram_above_the_spectrum_maximum_sets_the_scale() {
    assert_matches_cpp(
        "map_normalizer_chromatogram_above_input.mzML",
        "map_normalizer_chromatogram_above_output.mzML",
    );
}

/// The other side of the same maximum: a chromatogram below every peak changes
/// nothing, and the scale stays the spectrum maximum — which here lives in the
/// MS2 spectrum the tool never rescales, so MS1 normalizes to 6, not to 100.
#[test]
fn a_chromatogram_below_the_spectrum_maximum_leaves_the_scale_alone() {
    assert_matches_cpp(
        "map_normalizer_chromatogram_below_input.mzML",
        "map_normalizer_chromatogram_below_output.mzML",
    );
}

/// A run whose *combined intensity* range is empty is refused, and no output
/// file is written — as the source refuses it.
///
/// `main_` asks `exp.getMaxIntensity()` unconditionally, one line before the
/// peak loop (`MapNormalizer.cpp:93-94`), and `RangeBase::getMax()`
/// (`RangeManager.h:139-146` at core bc9cc12) throws `Exception::InvalidRange`
/// on an empty range with no assertion guard, so a Release build throws too.
/// Executed at the pinned Release build on both fixtures below: exit 8, *Empty
/// or uninitialized range object. Did you forget to call updateRanges()?*, no
/// output file. The port refuses both and writes nothing, but exits 6, because
/// `Error::InvalidValue` maps to `ILLEGAL_PARAMETERS` crate-wide where the
/// source's same-named exception reaches its `UNKNOWN_ERROR` arm; see
/// `run_failure` in `src/cli.rs`. The C++ FileInfo prints
/// `intensity: <none> .. <none>` under *Combined Ranges* for both.
///
/// The two shapes are the two ways to get there: scans that carry a retention
/// time but no point at all (the combined *RT* range is non-empty, the
/// intensity range is not), and a run with no spectra and no chromatograms.
#[test]
fn an_empty_combined_intensity_range_is_refused() {
    for name in [
        "map_normalizer_empty_range_input.mzML",
        "map_normalizer_empty_run_input.mzML",
    ] {
        let temp = TempDir::new_in(std::env::temp_dir(), false).unwrap();
        let out = temp.path().join("out.mzML");
        let (code, err) = run(&[
            "-test",
            "-in",
            &fixture(name).to_string_lossy(),
            "-out",
            &out.to_string_lossy(),
        ]);
        assert_eq!(code, ExitCode::IllegalParameters, "{name}: {err}");
        assert!(
            err.contains("run has no intensities to normalize"),
            "{name}: unexpected diagnostic {err}"
        );
        assert!(
            !out.exists(),
            "{name}: the C++ writes no output file here, and neither may the port"
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
