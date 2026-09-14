// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The upstream `TOPP_SpectraFilterWindowMower_1` test, compared canonically
//! against the retained C++ output, plus the algorithm-subsection plumbing that
//! every algorithm-wrapping TOPP tool depends on.

// The TOPP framework lives behind `paramxml` (every tool supports -ini) and these
// tools read mzML, so the whole file is inert without both features. Without this
// gate `cargo test --no-default-features` fails to compile.
#![cfg(all(feature = "mzml", feature = "paramxml"))]

use openms::cli::tools::SpectraFilterWindowMower;
use openms::cli::{ExitCode, TEST_MODE_COMPLETION_TIME, TOPP_PRODUCT_VERSION, run_with};
use openms::data_structures::DateTime;
use openms::format::controlled_vocabulary::ControlledVocabulary;
use openms::format::file_handler::FileHandler;
use openms::format::file_types::FileType;
use openms::kernel::MSExperiment;
use openms::metadata::{MetaValue, ProcessingAction};
use openms::system::file::TempDir;
use std::path::{Path, PathBuf};

/// Completion times as the C++ mzML writer keeps them: to the minute
/// (`MzMLHandler.cpp:3947` writes `yyyy-MM-dd+hh:mm`).
fn minutes(time: Option<DateTime>) -> Option<String> {
    time.map(|t| t.format("yyyy-MM-dd+hh:mm").unwrap())
}

/// The software name a record named `name` has after the C++ mzML writer and
/// reader: the writer stores the first PSI-MS software term below `MS:1000531`
/// named `name`, `name software` or `TOPP name`, and a custom software with
/// the name as its value otherwise (`MzMLHandler.cpp:3763-3787`); the reader
/// takes the term's name. So the record `SpectraFilterWindowMower` reloads as
/// `TOPP SpectraFilterWindowMower` (`MS:1002146`). This port's writer
/// transports the exact name instead (`src/format/mzml_header/write.rs`).
fn source_written_software_name(name: &str) -> String {
    let cv = ControlledVocabulary::psi_ms().unwrap();
    for candidate in [
        name.to_owned(),
        format!("{name} software"),
        format!("TOPP {name}"),
    ] {
        if let Some(term) = cv.first_child_with_name("MS:1000531", &candidate).unwrap() {
            return if term.id == "MS:1000799" {
                String::new()
            } else {
                term.name.clone()
            };
        }
    }
    name.to_owned()
}

/// The processing records of every spectrum and chromatogram agree with the
/// retained C++ output: the same records in the same order, each with the same
/// software name and version, actions, completion time and metadata.
///
/// Decision D4: the source tool attaches `getProcessingInfo_` to its output
/// (`SpectraFilterWindowMower.cpp`, `addDataProcessing_(exp,
/// getProcessingInfo_(DataProcessing::FILTERING))`), and the retained
/// `SpectraFilterWindowMower_1_output.mzML` carries that `data filtering` record
/// after the two the input already had. Completion times are compared to the
/// minute (see [`minutes`]) and software names as the C++ writer stores them
/// (see [`source_written_software_name`]).
fn assert_same_processing(produced: &MSExperiment, expected: &MSExperiment) {
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
            assert_eq!(
                source_written_software_name(&pa.software.name),
                pe.software.name,
                "{at}: software"
            );
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
fn run(args: &[&str]) -> (ExitCode, String, String) {
    let arguments: Vec<String> = std::iter::once("SpectraFilterWindowMower".to_string())
        .chain(args.iter().map(|a| (*a).to_string()))
        .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<SpectraFilterWindowMower>(&arguments, &mut out, &mut err);
    (
        code,
        String::from_utf8_lossy(&out).into_owned(),
        String::from_utf8_lossy(&err).into_owned(),
    )
}
/// A uniquely named directory for one case, removed when the guard drops.
fn workdir() -> TempDir {
    TempDir::new_in(std::env::temp_dir(), false).unwrap()
}

#[test]
fn topp_spectra_filter_window_mower_1_matches_the_retained_output() {
    let temp = workdir();
    let out = temp.path().join("out.mzML");
    let (code, _, err) = run(&[
        "-test",
        "-in",
        &fixture("window_mower_tool_input.mzML").to_string_lossy(),
        "-out",
        &out.to_string_lossy(),
    ]);
    assert_eq!(code, ExitCode::ExecutionOk, "{err}");

    let produced = load(&out);
    let reference = load(fixture("window_mower_tool_output.mzML"));
    assert_eq!(produced.spectra.len(), reference.spectra.len());
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
                (f64::from(pa.intensity) - f64::from(pe.intensity)).abs() <= 1e-6,
                "spectrum {index} peak {i} intensity"
            );
        }
    }
    assert_same_processing(&produced, &reference);
    for spectrum in &produced.spectra {
        assert_eq!(
            spectrum.data_processing.last().unwrap().completion_time,
            Some(DateTime::parse(TEST_MODE_COMPLETION_TIME).unwrap())
        );
    }
}

/// The processing record outside `-test`, against the C++ product SDK's output
/// for the same command (oracle `window_mower_notest` in
/// `../oracle/topp-cli-lifecycle/cli2/manifest.json`, retained as
/// `tests/data/topp_cli_lifecycle/swm_notest_output.mzML`).
///
/// The source records the product version, the time of the run and every
/// resolved parameter as `parameter: <name>` with its typed value
/// (`TOPPBase.cpp:556-568`). Compared: the records the input already carried
/// in full, and in the new record the software, version, action, the metadata
/// keys in order and every value except the two paths, which name the oracle's
/// and this run's own files; the completion time is only required to be set.
#[test]
fn the_processing_record_outside_test_mode_matches_the_cpp_output() {
    let temp = workdir();
    let out = temp.path().join("notest.mzML");
    let input = fixture("window_mower_tool_input.mzML")
        .to_string_lossy()
        .into_owned();
    let (code, _, err) = run(&["-in", &input, "-out", &out.to_string_lossy()]);
    assert_eq!(code, ExitCode::ExecutionOk, "{err}");

    let produced = load(&out);
    let oracle = load("tests/data/topp_cli_lifecycle/swm_notest_output.mzML");
    assert_eq!(produced.spectra.len(), oracle.spectra.len());
    for (index, (a, e)) in produced.spectra.iter().zip(&oracle.spectra).enumerate() {
        assert_eq!(a.data_processing.len(), 3, "spectrum {index}");
        assert_eq!(e.data_processing.len(), 3, "spectrum {index}");
        let (earlier_a, earlier_e) = (&a.data_processing[..2], &e.data_processing[..2]);
        for (pa, pe) in earlier_a.iter().zip(earlier_e) {
            assert_eq!(
                source_written_software_name(&pa.software.name),
                pe.software.name
            );
            assert_eq!(pa.software.version, pe.software.version);
            assert_eq!(pa.actions, pe.actions);
            assert_eq!(minutes(pa.completion_time), minutes(pe.completion_time));
            assert_eq!(pa.metadata, pe.metadata);
        }
        let (new_a, new_e) = (&a.data_processing[2], &e.data_processing[2]);
        assert_eq!(new_a.software.name, "SpectraFilterWindowMower");
        assert_eq!(
            source_written_software_name(&new_a.software.name),
            new_e.software.name
        );
        assert_eq!(new_a.software.version, TOPP_PRODUCT_VERSION);
        assert_eq!(new_a.software.version, new_e.software.version);
        assert_eq!(
            new_a.actions.iter().copied().collect::<Vec<_>>(),
            vec![ProcessingAction::DataFiltering]
        );
        assert_eq!(new_a.actions, new_e.actions);
        assert!(new_a.completion_time.is_some());
        let keys = |record: &openms::metadata::DataProcessing| {
            record.metadata.keys().cloned().collect::<Vec<_>>()
        };
        assert_eq!(keys(new_a), keys(new_e));
        for (key, value) in &new_a.metadata {
            match key.as_str() {
                "parameter: in" => assert_eq!(value, &MetaValue::from(input.as_str())),
                "parameter: out" => {
                    assert_eq!(value, &MetaValue::from(out.to_string_lossy().into_owned()))
                }
                _ => assert_eq!(Some(value), new_e.metadata.get(key), "{key}"),
            }
        }
    }
}

/// A non-finite window size never reaches the processing record: this port's
/// window mower refuses it first, exit 6, and nothing is written. The C++ tool
/// exits 0 and records `parameter: algorithm:windowsize` as `xsd:double` `inf`
/// (oracle `window_mower_notest_windowsize_inf` in
/// `../oracle/topp-cli-lifecycle/cli2/manifest.json`); `MetaValue` cannot hold
/// it. Both are documented native differences in `docs/TOPP_CLI_SUPPORT.md`
/// (*Non-finite values in processing records*).
#[test]
fn a_non_finite_window_size_is_refused_before_the_processing_record() {
    let temp = workdir();
    let out = temp.path().join("inf.mzML");
    let (code, _, err) = run(&[
        "-in",
        &fixture("window_mower_tool_input.mzML").to_string_lossy(),
        "-out",
        &out.to_string_lossy(),
        "-algorithm:windowsize",
        "inf",
    ]);
    assert_eq!(code, ExitCode::IllegalParameters, "{err}");
    assert!(err.contains("window width must be finite"), "{err}");
    assert!(!out.exists());
}

#[test]
fn the_algorithm_subsection_carries_the_source_defaults() {
    let temp = workdir();
    let dir = temp.path();
    let ini = dir.join("tool.ini");
    let (code, _, err) = run(&["-write_ini", &ini.to_string_lossy()]);
    assert_eq!(code, ExitCode::ExecutionOk, "{err}");

    // The subsection must reach the INI with the C++ WindowMower defaults.
    let loaded = openms::format::paramxml::load(&ini).unwrap();
    let key = |k: &str| format!("SpectraFilterWindowMower:1:algorithm:{k}");
    assert_eq!(
        loaded.value(&key("windowsize")).unwrap().to_f64().unwrap(),
        50.0
    );
    assert_eq!(
        loaded.value(&key("peakcount")).unwrap().to_i64().unwrap(),
        2
    );
    assert_eq!(
        loaded.value(&key("movetype")).unwrap().as_str().unwrap(),
        "slide"
    );
}

#[test]
fn a_subsection_value_from_an_ini_changes_the_result() {
    let temp = workdir();
    let dir = temp.path();
    let ini = dir.join("tool.ini");
    assert_eq!(
        run(&["-write_ini", &ini.to_string_lossy()]).0,
        ExitCode::ExecutionOk
    );

    // Keeping one peak per window must retain no more peaks than the default of two.
    let mut param = openms::format::paramxml::load(&ini).unwrap();
    param
        .set_value(
            "SpectraFilterWindowMower:1:algorithm:peakcount",
            openms::param::ParamValue::Integer(1),
            "",
            &[],
        )
        .unwrap();
    openms::format::paramxml::store(&ini, &param).unwrap();

    let out = dir.join("one.mzML");
    let (code, _, err) = run(&[
        "-test",
        "-ini",
        &ini.to_string_lossy(),
        "-in",
        &fixture("window_mower_tool_input.mzML").to_string_lossy(),
        "-out",
        &out.to_string_lossy(),
    ]);
    assert_eq!(code, ExitCode::ExecutionOk, "{err}");

    let one: usize = load(&out).spectra.iter().map(|s| s.peaks.len()).sum();
    let two: usize = load(fixture("window_mower_tool_output.mzML"))
        .spectra
        .iter()
        .map(|s| s.peaks.len())
        .sum();
    assert!(
        one <= two,
        "peakcount=1 kept {one} peaks, default kept {two}"
    );
}

/// An INI value that violates a registered restriction is **rejected**. The
/// source applies INI and command-line values with
/// `param_.update(finalParam, false, false, true, true, ...)`, that is with
/// `fail_on_invalid_values` and `fail_on_unknown_parameters` set
/// (`TOPPBase.cpp:339-342`), and exits with ILLEGAL_PARAMETERS instead of
/// falling back to the default. The C++ product SDK does exactly that for this
/// INI value (`../oracle/topp-cli-lifecycle`, case
/// `ini_invalid_subsection_value`).
#[test]
fn an_invalid_subsection_value_is_rejected() {
    let temp = workdir();
    let dir = temp.path();
    let ini = dir.join("tool.ini");
    assert_eq!(
        run(&["-write_ini", &ini.to_string_lossy()]).0,
        ExitCode::ExecutionOk
    );
    let mut param = openms::format::paramxml::load(&ini).unwrap();
    param
        .set_value(
            "SpectraFilterWindowMower:1:algorithm:movetype",
            openms::param::ParamValue::String("sideways".into()),
            "",
            &[],
        )
        .unwrap();
    openms::format::paramxml::store(&ini, &param).unwrap();
    let out = dir.join("x.mzML");
    let (code, _, err) = run(&[
        "-test",
        "-ini",
        &ini.to_string_lossy(),
        "-in",
        &fixture("window_mower_tool_input.mzML").to_string_lossy(),
        "-out",
        &out.to_string_lossy(),
    ]);
    assert_eq!(code, ExitCode::IllegalParameters, "{err}");
    assert!(
        err.contains("Parameters passed to 'SpectraFilterWindowMower' are invalid"),
        "{err}"
    );
    assert!(!out.exists(), "no output is written after a refused update");
}
