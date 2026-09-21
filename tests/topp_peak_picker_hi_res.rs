// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The `PeakPickerHiRes` TOPP tool against the C++ outputs.
//!
//! * The upstream registrations `TOPP_PeakPickerHiRes_1`, `_2`, `_5` and `_6`
//!   (test-data `0cb15f2` `topp/CMakeLists.txt:2523-2549`), compared by decoded
//!   content against the retained outputs (decision D6).
//! * The six parameter-failure registrations and `TOPPWRITEINI_OVERWRITE`
//!   (`CMakeLists.txt:104-131`), with the diagnostics `ExpectToolFailure.cmake`
//!   requires and the line-level `FuzzyDiff` comparison of the INI.
//! * The C1 oracle regressions (`../oracle/topp-early-bundle`): the centroided
//!   refusal without `-force`, `SignalToNoise:auto_mode 1` and the command-line
//!   override `-algorithm:signal_to_noise 2`.
//! * The P3 oracle (`../oracle/topp-peak-picker-tool`): the INI the C++ tool
//!   writes, fed back unchanged at several thread counts, an empty input,
//!   unsorted records, per-peak ion mobility and a run outside `-test`.
//! * The tool's load options against a synthetic document beyond the fixed
//!   ceilings the tool shipped with and beyond the library's strictness. The
//!   executed run on the 2.3 GB benchmark input, which no test suite can carry,
//!   is in `../oracle/topp-peak-picker-scale`.
//!
//! Every C++ output was produced by the product SDK (Debug, core `4fdec46`,
//! decision D7). Hashes, command lines and derivation rules are in
//! `tests/data/topp_peak_picker_hi_res_provenance.json`.

#![cfg(all(feature = "mzml", feature = "paramxml"))]

#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;

#[path = "support/decoded_compare.rs"]
mod decoded;

use base64::{Engine, engine::general_purpose::STANDARD};
use decoded::{DecodedOptions, Tolerance, compare_experiments};
use openms::cli::tools::PeakPickerHiRes;
use openms::cli::{ExitCode, TEST_MODE_COMPLETION_TIME, TOPP_PRODUCT_VERSION, run_with};
use openms::data_structures::DateTime;
use openms::format::PeakFileOptions;
use openms::format::controlled_vocabulary::ControlledVocabulary;
use openms::format::file_handler::FileHandler;
use openms::format::file_types::FileType;
use openms::format::mzml::{InputScaling, ReadOptions};
use openms::kernel::{MSExperiment, SpectrumType};
use openms::metadata::{DataProcessing, MetaValue, ProcessingAction};
use openms::system::file::TempDir;
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// A file under `tests/data`.
fn data(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(relative)
}

/// A file of this group's fixture directory.
fn fixture(name: &str) -> PathBuf {
    data("topp_peak_picker_hi_res").join(name)
}

/// The retained workflow inputs, shared with P1 (`tests/data/peak_picking/`)
/// and P2 (`tests/data/mzml_header_leniency/`) instead of copied.
fn workflow_input(workflow: u8) -> PathBuf {
    match workflow {
        1 => data("peak_picking/PeakPickerHiRes_input.mzML"),
        2 => data("peak_picking/PeakPickerHiRes_2_input.mzML"),
        5 => data("mzml_header_leniency/PeakPickerHiRes_5_input.mzML"),
        6 => data("peak_picking/PeakPickerHiRes_6_input.mzML"),
        _ => unreachable!("no such workflow"),
    }
}

fn text(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

/// Load an mzML file the way the tool does, through the tool's own load
/// options, so the two cannot drift apart.
fn load(path: impl AsRef<Path>) -> MSExperiment {
    FileHandler::load_experiment_with_read_options(
        path,
        &[FileType::MzMl],
        &PeakFileOptions::default(),
        &PeakPickerHiRes::read_options(),
    )
    .unwrap()
}

/// The outcome of one in-process run.
struct Run {
    code: ExitCode,
    out: String,
    err: String,
}

fn run(args: &[&str]) -> Run {
    let arguments: Vec<String> = std::iter::once("PeakPickerHiRes".to_owned())
        .chain(args.iter().map(|a| (*a).to_owned()))
        .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<PeakPickerHiRes>(&arguments, &mut out, &mut err);
    Run {
        code,
        out: String::from_utf8_lossy(&out).into_owned(),
        err: String::from_utf8_lossy(&err).into_owned(),
    }
}

/// A uniquely named directory for one case, removed when the guard drops.
fn workdir() -> TempDir {
    TempDir::new_in(std::env::temp_dir(), false).unwrap()
}

/// The software name a record named `name` has after the C++ mzML writer and
/// reader: the writer stores the first PSI-MS software term below `MS:1000531`
/// named `name`, `name software` or `TOPP name`, and a custom software with
/// the name as its value otherwise (`MzMLHandler.cpp:3763-3787`); the reader
/// takes the term's name. So `PeakPickerHiRes` reloads as
/// `TOPP PeakPickerHiRes` (`MS:1002135`). This port's writer transports the
/// exact name instead (`src/format/mzml_header/write.rs`).
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

/// Completion times as the C++ mzML writer keeps them: to the minute
/// (`MzMLHandler.cpp:3947` writes `yyyy-MM-dd+hh:mm`).
fn minutes(time: Option<DateTime>) -> Option<String> {
    time.map(|t| t.format("yyyy-MM-dd+hh:mm").unwrap())
}

/// A processing record as both mzML containers can carry it: the software name
/// as the C++ writer stores it, the completion time to the minute and no
/// software CV terms (the C++ writer's `MS:1002135` alias term is container
/// detail, documented and not compared).
fn container_neutral(record: &DataProcessing, rename: bool) -> Arc<DataProcessing> {
    let mut record = record.clone();
    if rename {
        record.software.name = source_written_software_name(&record.software.name);
    }
    record.completion_time = minutes(record.completion_time)
        .map(|text| DateTime::parse(&format!("{}:00", text.replace('+', " "))).unwrap());
    record.software.cv_terms = Default::default();
    Arc::new(record)
}

/// Apply [`container_neutral`] to every processing history of an experiment.
fn neutral(experiment: &MSExperiment, rename: bool) -> MSExperiment {
    let mut copy = experiment.clone();
    let convert = |history: &mut Vec<Arc<DataProcessing>>| {
        *history = history
            .iter()
            .map(|record| container_neutral(record, rename))
            .collect();
    };
    for spectrum in &mut copy.spectra {
        convert(&mut spectrum.data_processing);
        for array in &mut spectrum.float_data_arrays {
            convert(&mut array.data_processing);
        }
    }
    for chromatogram in &mut copy.chromatograms {
        convert(&mut chromatogram.data_processing);
        for array in &mut chromatogram.float_data_arrays {
            convert(&mut array.data_processing);
        }
    }
    copy
}

/// The decoded comparison of a Rust output with a C++ output (decision D6).
///
/// First the explicit contract of the package: record counts, native ids, MS
/// levels, spectrum types, float array names, every m/z, retention time,
/// intensity and float array value **exactly**, and the length of every
/// processing history. Then the complete decoded walk of `compare_experiments`
/// with exact numbers over both experiments with container-neutral processing
/// records (see [`container_neutral`]).
fn assert_decoded_equal(produced: &MSExperiment, expected: &MSExperiment) {
    assert_eq!(produced.spectra.len(), expected.spectra.len(), "spectra");
    assert_eq!(
        produced.chromatograms.len(),
        expected.chromatograms.len(),
        "chromatograms"
    );
    for (index, (a, e)) in produced.spectra.iter().zip(&expected.spectra).enumerate() {
        let at = format!("spectrum {index}");
        assert_eq!(a.native_id, e.native_id, "{at}: native id");
        assert_eq!(a.ms_level, e.ms_level, "{at}: MS level");
        assert_eq!(a.spectrum_type, e.spectrum_type, "{at}: type");
        assert_eq!(a.peaks.len(), e.peaks.len(), "{at}: peaks");
        for (i, (pa, pe)) in a.peaks.iter().zip(&e.peaks).enumerate() {
            assert_eq!(pa.mz.to_bits(), pe.mz.to_bits(), "{at}, peak {i}: m/z");
            assert_eq!(
                pa.intensity.to_bits(),
                pe.intensity.to_bits(),
                "{at}, peak {i}: intensity"
            );
        }
        let names = |s: &openms::kernel::MSSpectrum| {
            s.float_data_arrays
                .iter()
                .map(|array| array.name.clone())
                .collect::<Vec<_>>()
        };
        assert_eq!(names(a), names(e), "{at}: float arrays");
        for (aa, ae) in a.float_data_arrays.iter().zip(&e.float_data_arrays) {
            let bits = |values: &[f32]| values.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
            assert_eq!(bits(&aa.data), bits(&ae.data), "{at}: {}", aa.name);
        }
        assert_eq!(
            a.data_processing.len(),
            e.data_processing.len(),
            "{at}: processing history"
        );
    }
    for (index, (a, e)) in produced
        .chromatograms
        .iter()
        .zip(&expected.chromatograms)
        .enumerate()
    {
        let at = format!("chromatogram {index}");
        assert_eq!(a.native_id, e.native_id, "{at}: native id");
        assert_eq!(a.peaks.len(), e.peaks.len(), "{at}: peaks");
        for (i, (pa, pe)) in a.peaks.iter().zip(&e.peaks).enumerate() {
            assert_eq!(pa.rt.to_bits(), pe.rt.to_bits(), "{at}, peak {i}: RT");
            assert_eq!(
                pa.intensity.to_bits(),
                pe.intensity.to_bits(),
                "{at}, peak {i}: intensity"
            );
        }
        let names = |c: &openms::kernel::MSChromatogram| {
            c.float_data_arrays
                .iter()
                .map(|array| array.name.clone())
                .collect::<Vec<_>>()
        };
        assert_eq!(names(a), names(e), "{at}: float arrays");
        for (aa, ae) in a.float_data_arrays.iter().zip(&e.float_data_arrays) {
            let bits = |values: &[f32]| values.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
            assert_eq!(bits(&aa.data), bits(&ae.data), "{at}: {}", aa.name);
        }
        assert_eq!(
            a.data_processing.len(),
            e.data_processing.len(),
            "{at}: processing history"
        );
    }
    compare_experiments(
        &neutral(produced, true),
        &neutral(expected, false),
        &DecodedOptions::new(Tolerance::exact()),
    )
    .unwrap();
}

/// Every spectrum and chromatogram ends with the same shared test-mode
/// `peak picking` record (`addDataProcessing_` under `-test`).
fn assert_final_test_mode_record(produced: &MSExperiment) {
    let histories = produced
        .spectra
        .iter()
        .map(|s| &s.data_processing)
        .chain(produced.chromatograms.iter().map(|c| &c.data_processing));
    let mut first: Option<&Arc<DataProcessing>> = None;
    for history in histories {
        let last = history.last().expect("a processing record");
        assert_eq!(last.software.name, "PeakPickerHiRes");
        assert_eq!(last.software.version, "version_string");
        assert_eq!(
            last.actions.iter().copied().collect::<Vec<_>>(),
            vec![ProcessingAction::PeakPicking]
        );
        assert_eq!(
            minutes(last.completion_time),
            minutes(Some(DateTime::parse(TEST_MODE_COMPLETION_TIME).unwrap()))
        );
        assert_eq!(
            last.metadata.iter().collect::<Vec<_>>(),
            vec![(&"parameter: mode".to_owned(), &MetaValue::from("test_mode"))]
        );
        if let Some(first) = first {
            assert_eq!(**first, **last);
        }
        first = Some(last);
    }
}

/// The lengths of every record's processing history, spectra then
/// chromatograms.
fn history_lengths(experiment: &MSExperiment) -> Vec<usize> {
    experiment
        .spectra
        .iter()
        .map(|s| s.data_processing.len())
        .chain(
            experiment
                .chromatograms
                .iter()
                .map(|c| c.data_processing.len()),
        )
        .collect()
}

const VERSION_WARNING_3_6_0: &str = "Warning: Parameters file version (3.6.0) does not match the version of this tool (1.0.0).\nYour current parameters are still valid, but there might be new valid values or even new parameters. Upgrading the INI might be useful.\n";

/// Run a registered workflow into a temporary directory and return the run,
/// the produced experiment and the guard.
fn workflow(ini: Option<&str>, input: &Path, extra: &[&str]) -> (Run, MSExperiment, TempDir) {
    let temp = workdir();
    let out = temp.path().join("PeakPickerHiRes.tmp.mzML");
    let ini_path = ini.map(|name| text(&fixture(name)));
    let (input, out_text) = (text(input), text(&out));
    let mut args = vec!["-test"];
    if let Some(ini) = &ini_path {
        args.extend(["-ini", ini]);
    }
    args.extend(["-in", &input, "-out", &out_text]);
    args.extend_from_slice(extra);
    let outcome = run(&args);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let produced = load(&out);
    (outcome, produced, temp)
}

/// `TOPP_PeakPickerHiRes_1` (`CMakeLists.txt:2523-2525`): two MS1 profile
/// spectra picked with the FWHM in ppm, three MS2 spectra copied.
#[test]
fn topp_peak_picker_hi_res_1_matches_the_retained_output() {
    let (outcome, produced, _temp) = workflow(
        Some("PeakPickerHiRes_parameters.ini"),
        &workflow_input(1),
        &[],
    );
    assert_eq!(
        outcome.out,
        format!(
            "{VERSION_WARNING_3_6_0}#Spectra that needed to and could be picked by MS-level:\n  MS-level 1: 2 / 2\n  MS-level 2: 0 / 3\n"
        )
    );
    assert_eq!(outcome.err, "");
    let expected = load(fixture("PeakPickerHiRes_output.mzML"));
    assert_decoded_equal(&produced, &expected);
    assert_eq!(history_lengths(&produced), vec![6; 5]);
    assert_eq!(
        produced
            .spectra
            .iter()
            .map(|s| s.peaks.len())
            .collect::<Vec<_>>(),
        expected
            .spectra
            .iter()
            .map(|s| s.peaks.len())
            .collect::<Vec<_>>()
    );
    assert!(
        produced
            .spectra
            .iter()
            .filter(|s| s.ms_level == 1)
            .all(|s| {
                s.float_data_arrays.len() == 1 && s.float_data_arrays[0].name == "FWHM_ppm"
            })
    );
    assert_final_test_mode_record(&produced);
}

/// `TOPP_PeakPickerHiRes_2` (`CMakeLists.txt:2527-2529`): five chromatograms,
/// no spectra.
#[test]
fn topp_peak_picker_hi_res_2_matches_the_retained_output() {
    let (outcome, produced, _temp) = workflow(
        Some("PeakPickerHiRes_parameters.ini"),
        &workflow_input(2),
        &[],
    );
    assert_eq!(
        outcome.out,
        format!(
            "{VERSION_WARNING_3_6_0}#Spectra that needed to and could be picked by MS-level:\n"
        )
    );
    assert_eq!(outcome.err, "");
    assert_decoded_equal(&produced, &load(fixture("PeakPickerHiRes_2_output.mzML")));
    assert_eq!(history_lengths(&produced), vec![2; 5]);
    assert_eq!(
        produced
            .chromatograms
            .iter()
            .map(|c| c.peaks.len())
            .collect::<Vec<_>>(),
        vec![2, 4, 1, 2, 3]
    );
    assert_final_test_mode_record(&produced);
}

/// `TOPP_PeakPickerHiRes_5` (`CMakeLists.txt:2543-2545`): no INI, an input with
/// a dangling instrument `softwareRef` and chromatogram
/// `defaultDataProcessingRef` that the source reader drops (decision D10).
#[test]
fn topp_peak_picker_hi_res_5_matches_the_retained_output() {
    let (outcome, produced, _temp) = workflow(None, &workflow_input(5), &[]);
    assert_eq!(
        outcome.out,
        "#Spectra that needed to and could be picked by MS-level:\n"
    );
    let expected = load(fixture("PeakPickerHiRes_5_output.mzML"));
    assert_decoded_equal(&produced, &expected);
    assert_eq!(
        history_lengths(&produced),
        vec![1; expected.chromatograms.len()]
    );
    assert!(
        produced
            .chromatograms
            .iter()
            .all(|c| c.peaks.len() == 1 && c.float_data_arrays.is_empty())
    );
    assert_final_test_mode_record(&produced);
}

/// `TOPP_PeakPickerHiRes_6` (`CMakeLists.txt:2547-2549`): a 1.8.0 INI with
/// `force` set as `type="bool"`, manual mode on MS level 1 and a spectrum the
/// type estimator calls centroided.
#[test]
fn topp_peak_picker_hi_res_6_matches_the_retained_output() {
    let (outcome, produced, _temp) =
        workflow(Some("PeakPickerHiRes_6.ini"), &workflow_input(6), &[]);
    assert!(
        outcome.out.ends_with(
            "#Spectra that needed to and could be picked by MS-level:\n  MS-level 1: 1 / 1\n"
        ),
        "{}",
        outcome.out
    );
    assert_decoded_equal(&produced, &load(fixture("PeakPickerHiRes_6_output.mzML")));
    assert_eq!(history_lengths(&produced), vec![3]);
    assert_eq!(produced.spectra[0].peaks.len(), 5);
    assert_final_test_mode_record(&produced);
}

/// The six registered parameter failures (`CMakeLists.txt:107-131`) exit 6
/// with every diagnostic `ExpectToolFailure.cmake` requires as a substring of
/// the combined output.
#[test]
fn registered_parameter_failures_exit_6_with_the_registered_diagnostics() {
    let invalid_ini = |name: &str| text(&fixture(name));
    let cases: [(&str, Vec<String>, &[&str]); 6] = [
        (
            "TOPP_INI_INVALIDVALUE",
            vec!["-ini".into(), invalid_ini("PPHiRes_invalidValue.ini")],
            &[
                "Invalid string parameter value 'THISVALUEISINVALID' for parameter 'processOption'",
                "Parameters passed to 'PeakPickerHiRes' are invalid",
            ],
        ),
        (
            "TOPP_CLI_INVALIDVALUE",
            vec!["-processOption".into(), "THISVALUEISINVALID".into()],
            &[
                "Invalid string parameter value 'THISVALUEISINVALID' for parameter 'processOption'",
                "Parameters passed to 'PeakPickerHiRes' are invalid",
            ],
        ),
        (
            "TOPP_INI_INVALIDVALUE_SECTION",
            vec![
                "-ini".into(),
                invalid_ini("PPHiRes_invalidValueSection.ini"),
            ],
            &[
                "Invalid string parameter value 'THISVALUEISINVALID' for parameter 'report_FWHM_unit'",
                "Parameters passed to 'PeakPickerHiRes' are invalid",
            ],
        ),
        (
            "TOPP_CLI_INVALIDVALUE_SECTION",
            vec![
                "-algorithm:report_FWHM_unit".into(),
                "THISVALUEISINVALID".into(),
            ],
            &[
                "Invalid string parameter value 'THISVALUEISINVALID' for parameter 'report_FWHM_unit'",
                "Parameters passed to 'PeakPickerHiRes' are invalid",
            ],
        ),
        (
            "TOPP_INI_INVALIDNAME",
            vec!["-ini".into(), invalid_ini("PPHiRes_invalidParamName.ini")],
            &[
                "Unknown (or deprecated) Parameter 'INVALIDPARAMNAME_BASE'",
                "Unknown (or deprecated) Parameter 'algorithm:invalidParamName'",
                "Parameters passed to 'PeakPickerHiRes' are invalid",
            ],
        ),
        (
            "TOPP_CLI_INVALIDNAME",
            vec!["-algorithm:invalidParamName".into(), "somevalue".into()],
            &[
                "Unknown option(s)",
                "-algorithm:invalidParamName",
                "given. Aborting!",
            ],
        ),
    ];
    for (name, args, diagnostics) in cases {
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        let outcome = run(&args);
        assert_eq!(outcome.code, ExitCode::IllegalParameters, "{name}");
        let combined = format!("{}\n{}", outcome.out, outcome.err);
        for diagnostic in diagnostics {
            assert!(
                combined.contains(diagnostic),
                "{name}: missing {diagnostic:?} in {combined}"
            );
        }
    }
}

/// `TOPPWRITEINI_OVERWRITE` (`CMakeLists.txt:104-106`): `-write_ini` updated
/// leniently from an INI of the retired `PeakPickerWavelet`, compared line by
/// line with the registered `FuzzyDiff` whitelist `version` against the
/// retained `WRITE_INI_OUT.ini`, which P1 keeps in `tests/data/peak_picking/`.
#[test]
fn toppwriteini_overwrite_matches_the_retained_ini() {
    let temp = workdir();
    let written = temp.path().join("WRITE_INI.tmp.ini");
    let outcome = run(&[
        "-test",
        "-write_ini",
        &text(&written),
        "-ini",
        &text(&fixture("WRITE_INI_IN.ini")),
    ]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    for line in [
        "Warning: The provided INI file does not contain any parameters specific for this tool (expected in 'PeakPickerHiRes:1:').",
        "Default-Parameter 'PeakPickerHiRes:1:debug' overridden: '0' --> '4'!",
        "Default-Parameter 'PeakPickerHiRes:1:algorithm:signal_to_noise' overridden: '0' --> '1'!",
        "Unknown (or deprecated) Parameter 'PeakPickerWavelet:1:in_OLDNAME' given in outdated parameter file! Ignoring parameter. ",
        "Found 'PeakPickerWavelet:1:algorithm:SignalToNoiseEstimationParameter:auto_mode' as 'PeakPickerHiRes:1:algorithm:SignalToNoise:auto_mode' in new param.",
    ] {
        assert!(
            outcome.err.contains(line),
            "missing {line:?} in {}",
            outcome.err
        );
    }
    let settings = fuzzy::FuzzyDiffSettings::upstream()
        .unwrap()
        .with_whitelist(&["version"]);
    let verdict = fuzzy::fuzzy_diff(&written, &data("peak_picking/WRITE_INI_OUT.ini"), &settings);
    assert!(verdict.passed(), "{}", verdict.log_text());
}

/// `-write_ini` without `-ini` writes what the C++ tool writes (P3 oracle
/// `write_ini`), line by line under the upstream `FuzzyDiff` settings without
/// a `version` whitelist: both report product version 1.0.0.
#[test]
fn write_ini_equals_the_ini_the_cpp_tool_writes() {
    let temp = workdir();
    let written = temp.path().join("PeakPickerHiRes.ini");
    let outcome = run(&["-write_ini", &text(&written)]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_eq!(outcome.err, "");
    let verdict = fuzzy::fuzzy_diff(
        &written,
        &fixture("oracle_write_ini.ini"),
        &fuzzy::FuzzyDiffSettings::upstream().unwrap(),
    );
    assert!(verdict.passed(), "{}", verdict.log_text());
    assert!(
        std::fs::read_to_string(&written)
            .unwrap()
            .contains(&format!("value=\"{TOPP_PRODUCT_VERSION}\""))
    );
}

/// The benchmark contract: the INI the C++ tool writes with `-write_ini` loads
/// unchanged (exit 0, no diagnostics), and the output at `-threads 1`, `16` and
/// `0` is byte-identical and decoded-equal to the C++ output for the same INI
/// (P3 oracle `cpp_ini_threads_1`, `_16` and `_0`, which are themselves
/// byte-identical).
#[test]
fn the_cpp_written_ini_is_accepted_unchanged_at_every_thread_count() {
    let temp = workdir();
    let ini = text(&fixture("oracle_write_ini.ini"));
    let input = text(&workflow_input(1));
    let expected = load(fixture("oracle_cpp_ini.mzML"));
    let mut bytes: Vec<Vec<u8>> = Vec::new();
    for threads in ["1", "2", "8", "16", "32", "0"] {
        let out = temp.path().join(format!("threads_{threads}.mzML"));
        let outcome = run(&[
            "-test",
            "-ini",
            &ini,
            "-in",
            &input,
            "-out",
            &text(&out),
            "-threads",
            threads,
        ]);
        assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
        assert_eq!(outcome.err, "", "-threads {threads}");
        assert_eq!(
            outcome.out,
            "#Spectra that needed to and could be picked by MS-level:\n  MS-level 1: 2 / 2\n  MS-level 2: 0 / 3\n"
        );
        assert_decoded_equal(&load(&out), &expected);
        bytes.push(std::fs::read(&out).unwrap());
    }
    assert!(
        bytes.iter().all(|b| *b == bytes[0]),
        "outputs differ by thread count"
    );

    // The benchmark harness runs without `-test`: the same INI, the same input
    // and the same centroids, with a processing record carrying every parameter.
    let out = temp.path().join("benchmark.mzML");
    let outcome = run(&[
        "-ini",
        &ini,
        "-in",
        &input,
        "-out",
        &text(&out),
        "-threads",
        "16",
    ]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let produced = load(&out);
    for (a, e) in produced.spectra.iter().zip(&expected.spectra) {
        assert_eq!(a.peaks, e.peaks);
    }
}

/// `-threads` on a registered workflow changes no output byte.
///
/// The tool body now runs on the worker pool `-threads` sizes
/// (`ToolContext::in_thread_pool`) and the picker's spectrum loop is parallel
/// inside it, so this is the tool-level half of the determinism contract: the
/// written file is fixed by the input and the parameters. The counts span one
/// worker, two, eight, thirty-two — more than an ordinary development machine
/// has, so the schedule genuinely differs — and `0`, every available processor.
/// The library half, over every fixture and both entry points, is
/// `tests/peak_picking_experiment.rs`; the whole-file comparison on the 2.3 GB
/// benchmark input is in `../oracle/topp-peak-picker-scale`.
#[test]
fn threads_do_not_change_the_workflow_1_output() {
    let temp = workdir();
    let ini = text(&fixture("PeakPickerHiRes_parameters.ini"));
    let input = text(&workflow_input(1));
    let mut bytes: Vec<Vec<u8>> = Vec::new();
    for threads in ["1", "2", "8", "16", "32", "0"] {
        let out = temp.path().join(format!("threads_{threads}.mzML"));
        let outcome = run(&[
            "-test",
            "-ini",
            &ini,
            "-in",
            &input,
            "-out",
            &text(&out),
            "-threads",
            threads,
        ]);
        assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
        bytes.push(std::fs::read(&out).unwrap());
    }
    assert!(
        bytes.iter().all(|b| *b == bytes[0]),
        "outputs differ by thread count"
    );
}

/// C1 oracle `PPHR_6_noforce`: workflow 6 with `force` false refuses the
/// centroided spectrum in manual mode with the source's message and
/// `UNKNOWN_ERROR`, and writes nothing.
#[test]
fn a_centroided_spectrum_without_force_exits_8_with_the_source_message() {
    let temp = workdir();
    let out = temp.path().join("PeakPickerHiRes_6.tmp.mzML");
    let outcome = run(&[
        "-test",
        "-ini",
        &text(&fixture("PeakPickerHiRes_6_noforce.ini")),
        "-in",
        &text(&workflow_input(6)),
        "-out",
        &text(&out),
    ]);
    assert_eq!(outcome.code, ExitCode::UnknownError);
    assert_eq!(
        outcome.err,
        "Error: Unexpected internal error (Error: Centroided data provided but profile spectra expected.)\n"
    );
    assert!(!out.exists());
    // The same refusal through a command-line subsection value.
    let outcome = run(&[
        "-test",
        "-in",
        &text(&workflow_input(6)),
        "-out",
        &text(&out),
        "-algorithm:ms_levels",
        "1",
    ]);
    assert_eq!(outcome.code, ExitCode::UnknownError, "{}", outcome.err);
    assert!(!out.exists());
}

/// C1 oracle `PPHR_auto_mode_1`: with noise estimation enabled the C++ tool
/// crashes in `AUTOMAXBYPERCENT` (SIGBUS or SIGSEGV); the port refuses
/// explicitly with `INCOMPATIBLE_INPUT_DATA` and writes nothing. With
/// `signal_to_noise` 0 the estimator never runs, and both implementations
/// succeed with the default output (P3 oracle `auto_mode_1_without_estimation`).
#[test]
fn signal_to_noise_auto_mode_1_is_refused_where_the_source_crashes() {
    let temp = workdir();
    let out = temp.path().join("PeakPickerHiRes_1.tmp.mzML");
    let outcome = run(&[
        "-test",
        "-ini",
        &text(&fixture("PeakPickerHiRes_parameters.ini")),
        "-in",
        &text(&workflow_input(1)),
        "-out",
        &text(&out),
        "-algorithm:SignalToNoise:auto_mode",
        "1",
    ]);
    assert_eq!(
        outcome.code,
        ExitCode::IncompatibleInputData,
        "{}",
        outcome.err
    );
    assert!(
        outcome.err.starts_with("Error: unsupported: ")
            && outcome.err.contains("auto_mode 1")
            && outcome.err.ends_with('\n'),
        "{}",
        outcome.err
    );
    assert!(!out.exists());

    let (_, produced, _guard) = workflow(
        None,
        &workflow_input(1),
        &["-algorithm:SignalToNoise:auto_mode", "1"],
    );
    assert_decoded_equal(&produced, &load(fixture("oracle_cpp_ini.mzML")));
}

/// C1 oracle `PPHR_cli_signal_to_noise_2`: command-line subsection values reach
/// the picker; `scan=12663` keeps 102 centroids, and the whole output is
/// decoded-equal to the C++ output.
#[test]
fn command_line_subsection_values_reach_the_picker() {
    let (_, produced, _temp) = workflow(
        None,
        &workflow_input(1),
        &[
            "-algorithm:signal_to_noise",
            "2",
            "-algorithm:ms_levels",
            "1",
        ],
    );
    let scan = produced
        .spectra
        .iter()
        .find(|s| s.native_id == "scan=12663")
        .expect("scan=12663");
    assert_eq!(scan.peaks.len(), 102);
    assert_decoded_equal(
        &produced,
        &load(fixture("oracle_cli_signal_to_noise_2.mzML")),
    );
}

/// One `-processOption lowmemory` run into a temporary directory: the run, the
/// bytes written and the guard. The bytes, not only the decoded experiment,
/// because the mode's container is part of what it produces.
fn low_memory(ini: Option<&str>, input: &Path, extra: &[&str]) -> (Run, Vec<u8>, LowMemoryRun) {
    let temp = workdir();
    let out = temp.path().join("PeakPickerHiRes_lowmem.tmp.mzML");
    let ini_path = ini.map(|name| text(&fixture(name)));
    let (input, out_text) = (text(input), text(&out));
    let mut args = vec!["-test"];
    if let Some(ini) = &ini_path {
        args.extend(["-ini", ini]);
    }
    args.extend([
        "-in",
        &input,
        "-out",
        &out_text,
        "-processOption",
        "lowmemory",
    ]);
    args.extend_from_slice(extra);
    let outcome = run(&args);
    let bytes = std::fs::read(&out).unwrap_or_default();
    (outcome, bytes, LowMemoryRun { out, _temp: temp })
}

/// The output path of a [`low_memory`] run and the guard that removes it.
struct LowMemoryRun {
    out: PathBuf,
    _temp: TempDir,
}

/// The same run in the in-memory mode, for the mode-to-mode comparisons.
fn in_memory_bytes(ini: Option<&str>, input: &Path, extra: &[&str]) -> Vec<u8> {
    let temp = workdir();
    let out = temp.path().join("PeakPickerHiRes_inmem.tmp.mzML");
    let ini_path = ini.map(|name| text(&fixture(name)));
    let (input, out_text) = (text(input), text(&out));
    let mut args = vec!["-test"];
    if let Some(ini) = &ini_path {
        args.extend(["-ini", ini]);
    }
    args.extend(["-in", &input, "-out", &out_text]);
    args.extend_from_slice(extra);
    assert_eq!(run(&args).code, ExitCode::ExecutionOk);
    std::fs::read(&out).unwrap()
}

/// `TOPP_PeakPickerHiRes_3` (test-data `0cb15f2` `topp/CMakeLists.txt:2533-2536`):
/// workflow 1 through `-processOption lowmemory`, compared against the retained
/// `PeakPickerHiRes_output_lowMem.mzML`.
///
/// The registration exists because of one attribute. Upstream's own comment on
/// it reads: the output "SHOULD be identical to 'PeakPickerHiRes_output.mzML',
/// but due to a missing 'Dataprocessing' entry (which is not known when writing
/// the mzML header), we need an extra output file". The two retained C++ files
/// differ in exactly one byte-length-preserving place, `dataProcessingList
/// count="3"` against `count="2"`, which
/// [`the_two_retained_cpp_outputs_differ_only_in_the_data_processing_count`]
/// pins; the C++ count is `max(1, dps.size() + float data arrays of the whole
/// experiment)` (`MzMLHandler.cpp:5160-5170`) and the consumer's header sees a
/// one-record dummy map, so it counts one record's arrays. This port's writer
/// counts the processing histories it writes, so that artifact has no analogue
/// here and both modes write the same count - a native difference of the mzML
/// writer, not of this mode.
///
/// Nothing reaches standard output: the per-MS-level summary belongs to
/// `pickExperiment`, which this path never calls.
#[test]
fn topp_peak_picker_hi_res_3_matches_the_retained_low_memory_output() {
    let (outcome, bytes, produced_at) = low_memory(
        Some("PeakPickerHiRes_parameters.ini"),
        &workflow_input(1),
        &[],
    );
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_eq!(outcome.out, VERSION_WARNING_3_6_0);
    assert_eq!(outcome.err, "");
    let produced = load(&produced_at.out);
    // The C++ Release build at the pins reproduces this retained file byte for
    // byte on this command line (oracle `p4-lowmemory`, case `w1_lowmem`,
    // sha256 `91b0fdb0...`), so comparing against it is comparing against that
    // build's own output.
    let expected = load(fixture("PeakPickerHiRes_output_lowMem.mzML"));
    assert_decoded_equal(&produced, &expected);
    assert_eq!(history_lengths(&produced), vec![6; 5]);
    assert_final_test_mode_record(&produced);
    // And, as upstream's comment wanted, identical to the in-memory run.
    assert_eq!(
        bytes,
        in_memory_bytes(
            Some("PeakPickerHiRes_parameters.ini"),
            &workflow_input(1),
            &[]
        )
    );
}

/// `TOPP_PeakPickerHiRes_4` (`CMakeLists.txt:2538-2540`): workflow 2, five
/// chromatograms and no spectra, through the low-memory mode. Upstream compares
/// this one against the **in-memory** output file, so there the two modes agree
/// in the source as well; `processChromatogram_` picks every chromatogram
/// unconditionally, exactly as `pickExperiment` does.
#[test]
fn topp_peak_picker_hi_res_4_matches_the_in_memory_retained_output() {
    let (outcome, bytes, produced_at) = low_memory(
        Some("PeakPickerHiRes_parameters.ini"),
        &workflow_input(2),
        &[],
    );
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_eq!(outcome.out, VERSION_WARNING_3_6_0);
    assert_eq!(outcome.err, "");
    let produced = load(&produced_at.out);
    assert_decoded_equal(&produced, &load(fixture("PeakPickerHiRes_2_output.mzML")));
    assert_eq!(history_lengths(&produced), vec![2; 5]);
    assert_eq!(
        produced
            .chromatograms
            .iter()
            .map(|c| c.peaks.len())
            .collect::<Vec<_>>(),
        vec![2, 4, 1, 2, 3]
    );
    assert_eq!(
        bytes,
        in_memory_bytes(
            Some("PeakPickerHiRes_parameters.ini"),
            &workflow_input(2),
            &[]
        )
    );
}

/// The retained C++ pair, byte for byte: the low-memory output differs from the
/// in-memory one in the `dataProcessingList count` attribute and nowhere else.
///
/// This pins the source divergence the registration exists for, from the two
/// files upstream retained, without depending on either implementation.
#[test]
fn the_two_retained_cpp_outputs_differ_only_in_the_data_processing_count() {
    let in_memory = std::fs::read(fixture("PeakPickerHiRes_output.mzML")).unwrap();
    let low_memory = std::fs::read(fixture("PeakPickerHiRes_output_lowMem.mzML")).unwrap();
    assert_eq!(in_memory.len(), low_memory.len());
    let differing: Vec<usize> = in_memory
        .iter()
        .zip(&low_memory)
        .enumerate()
        .filter(|(_, (a, b))| a != b)
        .map(|(i, _)| i)
        .collect();
    assert_eq!(differing.len(), 1, "{differing:?}");
    let at = differing[0];
    assert_eq!((in_memory[at], low_memory[at]), (b'3', b'2'));
    let context = &in_memory[at.saturating_sub(30)..at + 2];
    assert!(
        String::from_utf8_lossy(context).contains("dataProcessingList count=\""),
        "{}",
        String::from_utf8_lossy(context)
    );
}

/// The mode writes indexed mzML, as `doCleanup_` does through
/// `MzMLHandlerHelper::writeFooter_` with the consumer's inherited
/// `write_index_`: the retained `PeakPickerHiRes_output_lowMem.mzML` is an
/// `indexedmzML`, and so is this.
///
/// Every offset in the index is checked to land on the record element it names,
/// which is the property a downstream random-access reader depends on and the
/// reason the mode is worth using on a file too large to hold.
#[test]
fn the_low_memory_output_is_indexed_and_its_offsets_land_on_the_records() {
    let retained = std::fs::read(fixture("PeakPickerHiRes_output_lowMem.mzML")).unwrap();
    assert!(String::from_utf8_lossy(&retained).contains("<indexedmzML "));

    let (outcome, bytes, _produced_at) = low_memory(
        Some("PeakPickerHiRes_parameters.ini"),
        &workflow_input(1),
        &[],
    );
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let document = String::from_utf8(bytes).unwrap();
    assert!(document.contains("<indexedmzML "), "{document}");
    assert!(
        document.ends_with("</fileChecksum>\n</indexedmzML>\n"),
        "{}",
        &document[document.len() - 120..]
    );
    let entries: Vec<(&str, usize)> = document
        .lines()
        .filter_map(|line| {
            let rest = line.strip_prefix("<offset idRef=\"")?;
            let (id, rest) = rest.split_once("\">")?;
            let offset = rest.strip_suffix("</offset>")?.parse().ok()?;
            Some((id, offset))
        })
        .collect();
    assert_eq!(entries.len(), 5, "{entries:?}");
    for (id, offset) in entries {
        let at = &document[offset..];
        assert!(
            at.starts_with("<spectrum "),
            "{id} at {offset}: {:?}",
            &at[..40]
        );
        assert!(at.contains(&format!("id=\"{id}\"")), "{id} at {offset}");
    }
}

/// **The divergence that matters.** Automatic mode in the low-memory path tests
/// the *stored* spectrum type only - `s.getType()`, the `SpectrumSettings`
/// accessor `MSSpectrum` re-exposes (`MSSpectrum.h:655`) - where `pickExperiment`
/// tests `getType(true)`, which falls through to the data-processing history and
/// then to `PeakTypeEstimator` (`OpenMS4-topp/src/PeakPickerHiRes.cpp:119-124` against `510`).
///
/// Workflow 6's input carries `MS:1000525` and neither `MS:1000127` nor
/// `MS:1000128`, so its stored type is unknown while the estimator calls its 33
/// samples centroided. The in-memory mode therefore copies it and reports
/// `0 / 1`; the low-memory mode picks it, to four centroids. Two different
/// files from one input, and the source specifies both.
#[test]
fn low_memory_automatic_mode_tests_only_the_stored_spectrum_type() {
    let (outcome, in_memory, _temp) = {
        let (outcome, produced, temp) = workflow(None, &workflow_input(6), &[]);
        (outcome, produced, temp)
    };
    assert_eq!(
        outcome.out,
        "#Spectra that needed to and could be picked by MS-level:\n  MS-level 1: 0 / 1\n"
    );
    assert_eq!(in_memory.spectra.len(), 1);
    assert_eq!(in_memory.spectra[0].peaks.len(), 33);
    assert_eq!(in_memory.spectra[0].spectrum_type, SpectrumType::Unknown);

    let (outcome, _bytes, produced_at) = low_memory(None, &workflow_input(6), &[]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_eq!(outcome.out, "");
    let low = load(&produced_at.out);
    assert_eq!(low.spectra.len(), 1);
    assert_eq!(low.spectra[0].peaks.len(), 4);
    assert_eq!(low.spectra[0].spectrum_type, SpectrumType::Centroid);
    // The C++ Release build's own output on this command line, which is what
    // makes this a divergence of the source rather than a claim about it.
    assert_decoded_equal(&low, &load(fixture("oracle_lowmem_w6_auto.mzML")));
}

/// The low-memory path runs `pp.pick` with no spectrum-type check at all, so a
/// centroided spectrum on a selected MS level is picked instead of refused and
/// **`-force` is inert**: `check_spectrum_type` is read only by the experiment
/// entry points. The in-memory mode on the same input and parameters exits 8
/// with the source's message.
#[test]
fn low_memory_never_refuses_centroided_data_and_force_is_inert() {
    let refusal = run(&[
        "-test",
        "-in",
        &text(&workflow_input(6)),
        "-out",
        &text(&workdir().path().join("unused.mzML")),
        "-algorithm:ms_levels",
        "1",
    ]);
    assert_eq!(refusal.code, ExitCode::UnknownError);
    assert_eq!(
        refusal.err,
        "Error: Unexpected internal error (Error: Centroided data provided but profile spectra expected.)\n"
    );

    let (outcome, without_force, without_force_at) =
        low_memory(None, &workflow_input(6), &["-algorithm:ms_levels", "1"]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_eq!(outcome.err, "");
    let picked = load(&without_force_at.out);
    assert_eq!(picked.spectra[0].peaks.len(), 4);
    assert_decoded_equal(&picked, &load(fixture("oracle_lowmem_w6_auto.mzML")));

    let (outcome, with_force, _with_force_at) = low_memory(
        None,
        &workflow_input(6),
        &["-algorithm:ms_levels", "1", "-force"],
    );
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_eq!(with_force, without_force);
}

/// None of the in-memory mode's input phases exists on this path, in the source
/// or here: the per-peak ion mobility warning, the
/// [`ExitCode::IncompatibleInputData`] for an input without spectra and
/// chromatograms, and the two sortedness refusals are all `main_` code that
/// `doLowMemAlgorithm` returns before. Each input below is one the in-memory
/// mode reacts to; the low-memory mode writes a file and exits 0 on every one.
#[test]
fn low_memory_runs_none_of_the_in_memory_input_checks() {
    // Per-peak ion mobility: the in-memory mode warns once on standard error.
    let warned = run(&[
        "-test",
        "-in",
        &text(&fixture("p3_im_peak.mzML")),
        "-out",
        &text(&workdir().path().join("im.mzML")),
    ]);
    assert!(warned.err.contains("IM_PEAK"), "{}", warned.err);
    let (outcome, bytes, im_at) = low_memory(None, &fixture("p3_im_peak.mzML"), &[]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_eq!(outcome.err, "");
    assert!(!bytes.is_empty());
    // Silent, and the four centroids the C++ Release build writes here - the
    // automatic-mode divergence again, on a second input: the in-memory mode
    // copies this spectrum's 33 samples.
    assert_decoded_equal(
        &load(&im_at.out),
        &load(fixture("oracle_lowmem_im_peak.mzML")),
    );

    // An input without spectra and chromatograms: the in-memory mode exits 11.
    let empty = run(&[
        "-test",
        "-in",
        &text(&fixture("empty.mzML")),
        "-out",
        &text(&workdir().path().join("empty.mzML")),
    ]);
    assert_eq!(empty.code, ExitCode::IncompatibleInputData);
    let (outcome, bytes, _produced_at) = low_memory(None, &fixture("empty.mzML"), &[]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_eq!(outcome.err, "");
    // No record ever reached the consumer, so no header was written either: the
    // source's `doCleanup_` writes a footer only when `started_writing_` is set.
    assert!(bytes.is_empty(), "{}", String::from_utf8_lossy(&bytes));

    // Unsorted records: the in-memory mode's two sortedness refusals are
    // unreachable through its loader, which sorts every record by position
    // (`unsorted_records_are_sorted_on_load_and_picked`). The low-memory path
    // does not even contain the checks, and its reader sorts the same way, so
    // the two modes agree here rather than differing.
    for (name, ini, oracle) in [
        (
            "p3_unsorted_chromatogram.mzML",
            "PeakPickerHiRes_parameters.ini",
            "PeakPickerHiRes_2_output.mzML",
        ),
        (
            "p3_unsorted_spectrum.mzML",
            "PeakPickerHiRes_6.ini",
            "oracle_lowmem_unsorted_spectrum.mzML",
        ),
    ] {
        let (outcome, bytes, unsorted_at) = low_memory(Some(ini), &fixture(name), &[]);
        assert_eq!(
            outcome.code,
            ExitCode::ExecutionOk,
            "{name}: {}",
            outcome.err
        );
        assert_eq!(outcome.err, "", "{name}");
        assert_eq!(
            bytes,
            in_memory_bytes(Some(ini), &fixture(name), &[]),
            "{name}"
        );
        assert_decoded_equal(&load(&unsorted_at.out), &load(fixture(oracle)));
    }
}

/// The mode is serial - the source's consumer dispatch loop hands over one
/// record at a time and this port's reader does the same - so `-threads` reaches
/// nothing on this path and the written bytes cannot depend on it. Pinned at 1,
/// 8 and 32, on an input with spectra and on one with chromatograms.
#[test]
fn the_low_memory_output_is_bit_identical_at_every_thread_count() {
    for input in [workflow_input(1), workflow_input(2)] {
        let mut reference: Option<Vec<u8>> = None;
        for threads in ["1", "8", "32"] {
            let (outcome, bytes, _produced_at) = low_memory(
                Some("PeakPickerHiRes_parameters.ini"),
                &input,
                &["-threads", threads],
            );
            assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
            match &reference {
                None => reference = Some(bytes),
                Some(first) => assert_eq!(first, &bytes, "{input:?} at -threads {threads}"),
            }
        }
    }
}

/// A file derived from a retained input by one stated rule, in a directory of
/// its own that the returned guard removes.
///
/// The rules below are byte for byte the ones the C++ differential ran on
/// (`../oracle/p4-lowmemory/fixdiff_06.sh`), so the two implementations are
/// compared on the same bytes without a new fixture. The retained inputs are
/// pure ASCII despite their `ISO-8859-1` declaration, which is why a text rule
/// is a byte rule here.
fn derived(name: &str, from: &Path, rule: impl FnOnce(&str) -> String) -> (PathBuf, TempDir) {
    let temp = workdir();
    let path = temp.path().join(name);
    let source = std::fs::read(from).unwrap();
    let source = String::from_utf8(source).expect("the retained inputs are ASCII");
    std::fs::write(&path, rule(&source)).unwrap();
    (path, temp)
}

/// A corrupt input is a corrupt input in **both** process options.
///
/// `doLowMemAlgorithm` catches nothing (`OpenMS4-topp/src/PeakPickerHiRes.cpp:170-186`), so a
/// reader failure leaves `MzMLFile::transform` and reaches `TOPPBase::main`,
/// whose `ParseError` arm writes `Error: Unable to read file (…)` and returns
/// `INPUT_FILE_CORRUPT` (`TOPPBase.cpp:460-465`) — the same arm the in-memory
/// mode's `loadExperiment` failures take. The exit code and the diagnostic are
/// therefore the mode's, not the reader's caller's.
///
/// Executed against the C++ Release build on `ibminode06`
/// (`../oracle/p4-lowmemory`, `logs/fixdiff_06.log`, cases `trunc`, `garbage`
/// and `badb64`): on each input below the C++ tool exits 3 in both modes. The
/// low-memory runs leave a zero-byte output behind — the consumer's constructor
/// created the file before anything was read — and the in-memory runs leave
/// none, on both sides.
#[test]
fn a_corrupt_input_exits_3_in_both_modes() {
    for (name, rule) in [
        // Truncated inside the third record: neither pass can finish the XML.
        (
            "trunc.mzML",
            (|s: &str| s[..215_000].to_owned()) as fn(&str) -> String,
        ),
        // Not XML at all.
        ("garbage.mzML", |_: &str| {
            "this is not xml at all\n".repeat(100)
        }),
        // Well-formed XML whose third record carries base64 of a length no
        // decoder accepts. The C++ first pass steps over record contents
        // (`skip_spectrum_`, `MzMLHandler.cpp:966-974`), so the C++ low-memory
        // run reaches this only in its second pass and still writes nothing:
        // the refusal comes before the first record is handed over.
        ("bad_base64.mzML", |s: &str| {
            let record = s.find("<spectrum id=\"scan=12665\" index=\"2\"").unwrap();
            let binary = s[record..].find("<binary>").unwrap() + record + "<binary>".len();
            format!("{}!!!{}", &s[..binary], &s[binary..])
        }),
    ] {
        let (input, _temp) = derived(name, &workflow_input(1), rule);
        let (outcome, bytes, produced_at) = low_memory(None, &input, &[]);
        assert_eq!(
            outcome.code,
            ExitCode::InputFileCorrupt,
            "{name}: {}",
            outcome.err
        );
        assert!(
            outcome.err.starts_with("Error: Unable to read file ("),
            "{name}: {}",
            outcome.err
        );
        // The file the consumer's constructor created is still there and still
        // empty: no record reached the writer, so `doCleanup_` wrote no footer.
        assert!(produced_at.out.exists(), "{name}");
        assert!(bytes.is_empty(), "{name}");

        let temp = workdir();
        let out = temp.path().join("in_memory.tmp.mzML");
        let in_memory = run(&["-test", "-in", &text(&input), "-out", &text(&out)]);
        assert_eq!(in_memory.code, ExitCode::InputFileCorrupt, "{name}");
        assert!(
            in_memory.err.starts_with("Error: Unable to read file ("),
            "{name}: {}",
            in_memory.err
        );
        assert!(!out.exists(), "{name}");
    }
}

/// The one place this port is stricter than the source, pinned.
///
/// A `spectrumList count` that overstates the records that follow is written
/// silently by the source: the header carries the first pass's count, nothing
/// re-checks it, and its own class note says a wrong count "will lead to an
/// inconsistent mzML". The consumer's `CountPolicy::Checked` default is kept
/// here, so the run ends with an error — **after** the document has been closed
/// and indexed, so the file left behind is complete and readable, carrying the
/// same lying count the source writes.
///
/// Executed on `ibminode06` against the C++ Release build on this exact input
/// (`../oracle/p4-lowmemory`, `logs/fixdiff_06.log`, case `badcount`): the C++
/// low-memory run exits 0 and writes `<spectrumList count="9">` over five
/// records. Both in-memory runs are unaffected.
#[test]
fn a_low_memory_run_reports_a_lying_list_count_over_a_closed_document() {
    let (input, _temp) = derived("bad_count.mzML", &workflow_input(1), |s| {
        s.replacen("<spectrumList count=\"5\"", "<spectrumList count=\"9\"", 1)
    });
    let (outcome, bytes, produced_at) = low_memory(None, &input, &[]);
    assert_eq!(outcome.code, ExitCode::UnknownError, "{}", outcome.err);
    assert_eq!(
        outcome.err,
        "Error: Unexpected internal error (invalid value: mzML list counts announce (9, 0) \
         records but (5, 0) were written)\n"
    );
    // Closed, indexed, complete - and announcing the count the source announces.
    let written = String::from_utf8(bytes).unwrap();
    assert!(
        written.contains("<spectrumList count=\"9\""),
        "{written:.400}"
    );
    assert!(written.contains("</spectrumList>\n</run></mzML>\n"));
    assert!(written.ends_with("</indexedmzML>\n"));
    let reloaded = load(&produced_at.out);
    assert_eq!(reloaded.spectra.len(), 5);

    // The in-memory mode never sees the declared count as a promise.
    let temp = workdir();
    let out = temp.path().join("in_memory.tmp.mzML");
    let in_memory = run(&["-test", "-in", &text(&input), "-out", &text(&out)]);
    assert_eq!(in_memory.code, ExitCode::ExecutionOk, "{}", in_memory.err);
    assert_eq!(load(&out).spectra.len(), 5);
}

/// A failure after the first record still leaves a closed document.
///
/// The source's `~MSDataWritingConsumer` calls `doCleanup_` on every path
/// (`MSDataWritingConsumer.cpp:37-40`), which closes the open list and writes
/// the footer whenever writing started (`:151-173`), so a source run that
/// throws after the first record still leaves a closed, indexed document.
/// `MSDataWritingConsumer::finish` is this port's destructor, and
/// `run_low_memory` calls it on the failing path too.
///
/// `SignalToNoise:auto_mode 1` with `ms_levels 2` is the reachable case: the
/// MS1 spectrum of workflow 1 is copied and written, and the first MS2 spectrum
/// is the first record the picker touches, so the refusal of native difference
/// 2 arrives with one record already on disc. The source cannot be compared
/// here: the same command line on `ibminode06` writes all five records and then
/// dies of SIGSEGV (exit 139), leaving 404,915 bytes with no index and no
/// footer, because a signal runs no destructor (oracle `PPHR_auto_mode_1`, and
/// `../oracle/p4-lowmemory/logs/fixdiff_06.log` case `am1`). The assertions
/// below are therefore on the source's `doCleanup_` contract rather than on
/// executed C++ bytes.
#[test]
fn a_failure_after_the_first_record_still_closes_the_document() {
    let (outcome, bytes, produced_at) = low_memory(
        None,
        &workflow_input(1),
        &[
            "-algorithm:ms_levels",
            "2",
            "-algorithm:signal_to_noise",
            "1",
            "-algorithm:SignalToNoise:auto_mode",
            "1",
        ],
    );
    assert_eq!(
        outcome.code,
        ExitCode::IncompatibleInputData,
        "{}",
        outcome.err
    );
    assert!(outcome.err.contains("auto_mode 1"), "{}", outcome.err);
    let written = String::from_utf8(bytes).unwrap();
    assert!(
        written.contains("</spectrumList>\n</run></mzML>\n"),
        "truncated output"
    );
    assert!(written.ends_with("</indexedmzML>\n"), "no footer");
    // One record written, under the count the first pass declared: a streaming
    // writer cannot take back what it has already sent.
    assert!(written.contains("<spectrumList count=\"5\""));
    let reloaded = load(&produced_at.out);
    assert_eq!(reloaded.spectra.len(), 1);
    assert_eq!(reloaded.spectra[0].native_id, "scan=12663");
}

/// P3 oracle `empty`: an mzML without spectra and chromatograms (the C1 derived
/// `empty.mzML`) exits 11 with the source warning on standard error and
/// nothing on standard output, and writes nothing.
#[test]
fn an_input_without_spectra_and_chromatograms_exits_11() {
    let temp = workdir();
    let out = temp.path().join("empty.tmp.mzML");
    let outcome = run(&[
        "-test",
        "-in",
        &text(&fixture("empty.mzML")),
        "-out",
        &text(&out),
    ]);
    assert_eq!(outcome.code, ExitCode::IncompatibleInputData);
    assert_eq!(
        outcome.err,
        std::fs::read_to_string(fixture("oracle_empty.stderr.txt")).unwrap()
    );
    assert_eq!(outcome.out, "");
    assert!(!out.exists());
}

/// P3 oracle `unsorted_spectrum` and `unsorted_chromatogram`: the loader sorts
/// each record, in both implementations, before the tool's own sortedness
/// checks, so an unsorted record is picked. The derived inputs swap one pair
/// of samples of workflow 6 and workflow 2, so the outputs are those workflows'
/// retained outputs.
#[test]
fn unsorted_records_are_sorted_on_load_and_picked() {
    let (_, produced, _temp) = workflow(
        Some("PeakPickerHiRes_6.ini"),
        &fixture("p3_unsorted_spectrum.mzML"),
        &[],
    );
    assert_decoded_equal(&produced, &load(fixture("PeakPickerHiRes_6_output.mzML")));
    let (_, produced, _temp) = workflow(
        Some("PeakPickerHiRes_parameters.ini"),
        &fixture("p3_unsorted_chromatogram.mzML"),
        &[],
    );
    assert_decoded_equal(&produced, &load(fixture("PeakPickerHiRes_2_output.mzML")));
}

/// P3 oracle `im_peak`: per-peak ion mobility warns once on standard error
/// with the source text, and the intensity-weighted mean ion mobility array of
/// the picked spectrum equals the C++ output bit for bit.
#[test]
fn per_peak_ion_mobility_warns_and_reports_the_weighted_mobility() {
    let (outcome, produced, _temp) = workflow(
        Some("PeakPickerHiRes_6.ini"),
        &fixture("p3_im_peak.mzML"),
        &[],
    );
    assert_eq!(
        outcome.err,
        std::fs::read_to_string(fixture("oracle_im_peak.stderr.txt")).unwrap()
    );
    let expected = load(fixture("oracle_im_peak.mzML"));
    assert_eq!(
        produced.spectra[0].float_data_arrays[0].name,
        "mean ion mobility array"
    );
    assert_eq!(produced.spectra[0].float_data_arrays[0].data.len(), 5);
    assert_decoded_equal(&produced, &expected);
}

/// A diagnostic the input checks write survives a refusal raised later in the
/// run, because it reaches the real stream where the source writes it
/// (`OpenMS4-topp/src/PeakPickerHiRes.cpp:222-226`, before picking) rather than being collected
/// and written at the end.
///
/// The per-peak ion mobility input of `im_peak` with the `auto_mode 1` refusal
/// of `PPHR_auto_mode_1` on top: the source writes its warning and then dies in
/// `AUTOMAXBYPERCENT`, so the warning is what a user is left with. The port
/// must not swallow it.
#[test]
fn a_refusal_after_the_input_checks_keeps_the_ion_mobility_warning() {
    let temp = workdir();
    let out = temp.path().join("im_auto_mode.tmp.mzML");
    let outcome = run(&[
        "-test",
        "-ini",
        &text(&fixture("PeakPickerHiRes_6.ini")),
        "-in",
        &text(&fixture("p3_im_peak.mzML")),
        "-out",
        &text(&out),
        "-algorithm:signal_to_noise",
        "2",
        "-algorithm:SignalToNoise:auto_mode",
        "1",
    ]);
    assert_eq!(
        outcome.code,
        ExitCode::IncompatibleInputData,
        "{}",
        outcome.err
    );
    let warning = std::fs::read_to_string(fixture("oracle_im_peak.stderr.txt")).unwrap();
    let refusal = outcome
        .err
        .strip_prefix(&warning)
        .unwrap_or_else(|| panic!("the warning is written first: {}", outcome.err));
    assert!(
        refusal.starts_with("Error: unsupported: ")
            && refusal.contains("auto_mode 1")
            && refusal.ends_with('\n'),
        "{refusal}"
    );
    assert!(!out.exists());
}

/// The per-MS-level summary is written where the source writes it — after
/// picking and **before** the output is stored (`CENTROIDING/PeakPickerHiRes.cpp:559-563`,
/// then `:571`) — so a run whose store fails still reports what it picked.
///
/// The store is made to fail with an `-out` that is an existing directory: the
/// framework's writability pre-check (`src/cli/context.rs`, `file::writable`)
/// answers yes for a directory, so the run reaches the store, which cannot
/// write a file there.
#[test]
fn the_summary_is_written_before_the_output_is_stored() {
    let temp = workdir();
    let out = temp.path().join("occupied.mzML");
    std::fs::create_dir(&out).unwrap();
    let outcome = run(&[
        "-test",
        "-ini",
        &text(&fixture("PeakPickerHiRes_parameters.ini")),
        "-in",
        &text(&workflow_input(1)),
        "-out",
        &text(&out),
    ]);
    assert_ne!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_eq!(
        outcome.out,
        format!(
            "{VERSION_WARNING_3_6_0}#Spectra that needed to and could be picked by MS-level:\n  MS-level 1: 2 / 2\n  MS-level 2: 0 / 3\n"
        )
    );
    assert!(out.is_dir());
}

/// P3 oracle `notest`: outside `-test` the record carries the product version,
/// a completion time and every resolved parameter as `parameter: <name>`, in
/// the C++ order and with the C++ values, except the two paths, which name each
/// run's own files. The empty `algorithm:ms_levels` list reaches the output as
/// the source's `[]` text, which the tool renders because the native mzML
/// writer refuses list metadata.
#[test]
fn the_processing_record_outside_test_mode_matches_the_cpp_output() {
    let temp = workdir();
    let out = temp.path().join("notest.tmp.mzML");
    let input = text(&workflow_input(2));
    let outcome = run(&["-in", &input, "-out", &text(&out)]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let produced = load(&out);
    let oracle = load(fixture("oracle_notest.mzML"));
    assert_eq!(produced.chromatograms.len(), oracle.chromatograms.len());
    for (index, (a, e)) in produced
        .chromatograms
        .iter()
        .zip(&oracle.chromatograms)
        .enumerate()
    {
        assert_eq!(a.peaks, e.peaks, "chromatogram {index}");
        assert_eq!(a.data_processing.len(), 2, "chromatogram {index}");
        assert_eq!(e.data_processing.len(), 2, "chromatogram {index}");
        let (new_a, new_e) = (&a.data_processing[1], &e.data_processing[1]);
        assert_eq!(new_a.software.name, "PeakPickerHiRes");
        assert_eq!(
            source_written_software_name(&new_a.software.name),
            new_e.software.name
        );
        assert_eq!(new_a.software.version, TOPP_PRODUCT_VERSION);
        assert_eq!(new_a.software.version, new_e.software.version);
        assert_eq!(new_a.actions, new_e.actions);
        assert!(new_a.completion_time.is_some());
        let keys = |record: &DataProcessing| record.metadata.keys().cloned().collect::<Vec<_>>();
        assert_eq!(keys(new_a), keys(new_e));
        for (key, value) in &new_a.metadata {
            match key.as_str() {
                "parameter: in" => assert_eq!(value, &MetaValue::from(input.as_str())),
                "parameter: out" => {
                    assert_eq!(value, &MetaValue::from(text(&out)))
                }
                _ => assert_eq!(Some(value), new_e.metadata.get(key), "{key}"),
            }
        }
        assert_eq!(
            new_a.metadata.get("parameter: algorithm:ms_levels"),
            Some(&MetaValue::from("[]"))
        );
    }
}

/// Base64 of `bytes`, uncompressed, as the synthetic input below encodes its
/// binary arrays.
fn base64(bytes: &[u8]) -> String {
    STANDARD.encode(bytes)
}

/// One profile peak sampled five times at 0.01 Th spacing around its centre,
/// repeated `peaks` times at 10 Th intervals from 400 Th.
///
/// The samples of one peak are a Gaussian of 0.008 Th standard deviation, so
/// the picker finds one centroid per peak and no sample is zero.
fn profile_samples(peaks: usize) -> (Vec<f64>, Vec<f32>) {
    let (mut mz, mut intensity) = (Vec::new(), Vec::new());
    for peak in 0..peaks {
        let centre = 400.0 + peak as f64 * 10.0;
        for step in 0..5_i32 {
            let offset = f64::from(step - 2) * 0.01;
            mz.push(centre + offset);
            intensity.push((1_000.0 * (-(offset * offset) / (2.0 * 0.008 * 0.008)).exp()) as f32);
        }
    }
    (mz, intensity)
}

/// A synthetic profile mzML of `spectra` MS1 spectra, each carrying three
/// five-sample profile peaks and the parameter payload a converter writes, with
/// a dangling instrument `softwareRef` in the header.
///
/// The shape is the one `tests/mzml_reader_scale.rs` measures the former
/// fixed reader ceilings against: at 20,000 spectra the cumulative parameter
/// storage charge passes the fixed 512 MiB floor, while the size-derived
/// allowance of the same option set admits it. The dangling reference is what
/// the library default refuses and a tool path accepts (decision D10), so one
/// document exercises both differences between the tool's load options and the
/// library's. Returns the path and the document size in bytes.
fn instrument_scale_input(dir: &Path, spectra: usize) -> (PathBuf, u64) {
    let (mz, intensity) = profile_samples(3);
    let mz_text = base64(&mz.iter().flat_map(|v| v.to_le_bytes()).collect::<Vec<u8>>());
    let intensity_text = base64(
        &intensity
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect::<Vec<u8>>(),
    );
    let array = |accession: &str, bits: &str, text: &str| {
        format!(
            "<binaryDataArray encodedLength=\"{}\">\
             <cvParam cvRef=\"MS\" accession=\"{bits}\" name=\"bits\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>\
             <cvParam cvRef=\"MS\" accession=\"{accession}\" name=\"array\"/>\
             <binary>{text}</binary></binaryDataArray>",
            text.len()
        )
    };
    let arrays = format!(
        "<binaryDataArrayList count=\"2\">{}{}</binaryDataArrayList>",
        array("MS:1000514", "MS:1000523", &mz_text),
        array("MS:1000515", "MS:1000521", &intensity_text),
    );
    let mut xml = String::from(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\
         <mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">\
         <cvList count=\"1\"><cv id=\"MS\" fullName=\"PSI-MS\" URI=\"https://purl.obolibrary.org/obo/ms.obo\"/></cvList>\
         <fileDescription><fileContent/></fileDescription>\
         <referenceableParamGroupList count=\"1\"><referenceableParamGroup id=\"common\">\
         <cvParam cvRef=\"MS\" accession=\"MS:1000579\" name=\"MS1 spectrum\"/>\
         <cvParam cvRef=\"MS\" accession=\"MS:1000130\" name=\"positive scan\"/>\
         </referenceableParamGroup></referenceableParamGroupList>\
         <softwareList count=\"1\"><software id=\"sw\" version=\"1.0\"/></softwareList>\
         <instrumentConfigurationList count=\"1\"><instrumentConfiguration id=\"ic\">\
         <softwareRef ref=\"missing\"/></instrumentConfiguration></instrumentConfigurationList>\
         <dataProcessingList count=\"1\"><dataProcessing id=\"dp\">\
         <processingMethod order=\"0\" softwareRef=\"sw\">\
         <cvParam cvRef=\"MS\" accession=\"MS:1000544\" name=\"Conversion to mzML\"/>\
         </processingMethod></dataProcessing></dataProcessingList>\
         <run id=\"run\" defaultInstrumentConfigurationRef=\"ic\" startTimeStamp=\"2016-11-18T23:31:16\">",
    );
    xml += &format!("<spectrumList count=\"{spectra}\" defaultDataProcessingRef=\"dp\">");
    for index in 0..spectra {
        xml += &format!(
            "<spectrum id=\"controllerType=0 controllerNumber=1 scan={}\" index=\"{index}\" defaultArrayLength=\"{}\">\
             <referenceableParamGroupRef ref=\"common\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"1\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000128\" name=\"profile spectrum\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000504\" name=\"base peak m/z\" value=\"400.0\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000505\" name=\"base peak intensity\" value=\"1000.0\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000285\" name=\"total ion current\" value=\"12345.0\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000528\" name=\"lowest observed m/z\" value=\"399.98\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000527\" name=\"highest observed m/z\" value=\"420.02\"/>\
             <userParam name=\"filter string\" value=\"FTMS + p NSI Full ms [300.00-1500.00]\"/>\
             <userParam name=\"preset scan configuration\" value=\"1\"/>\
             <scanList count=\"1\"><cvParam cvRef=\"MS\" accession=\"MS:1000795\" name=\"no combination\"/>\
             <scan><cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"{}\" unitCvRef=\"UO\" unitAccession=\"UO:0000010\" unitName=\"second\"/>\
             <scanWindowList count=\"1\"><scanWindow>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000501\" name=\"scan window lower limit\" value=\"300.0\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000500\" name=\"scan window upper limit\" value=\"1500.0\"/>\
             </scanWindow></scanWindowList></scan></scanList>{arrays}</spectrum>",
            index + 1,
            mz.len(),
            index as f64 * 0.5,
        );
    }
    xml += "</spectrumList></run></mzML>";
    let path = dir.join("instrument_scale.tmp.mzML");
    std::fs::write(&path, &xml).unwrap();
    (path, xml.len() as u64)
}

/// The tool's load options admit an input that the library's own default
/// options refuse, in both ways they differ.
///
/// The scaffold this tool was written against passed the library defaults with
/// one switch flipped, and those defaults were fixed ceilings: 10,000,000 raw
/// points and 512 MiB of XML. No instrument-sized profile run fits under
/// either — the benchmark input `UK222.mzML` is 2.3 GB with 197,765,338 raw
/// points — so the Rust tool could not read the data the C++ tool reads.
/// [`PeakPickerHiRes::read_options`] is now `ReadOptions::source()`, whose
/// ceilings grow with the consumed input (`InputScaling`).
///
/// A 2.3 GB document cannot live in a test suite; the executed proof on that
/// input is in `docs/TOPP_PEAK_PICKER_HI_RES_SUPPORT.md`. What is checked here
/// is that the tool's options are the ones that admit a document beyond the
/// fixed floors and beyond the library's strictness, and that the run through
/// the tool is complete: every spectrum read, picked and written back.
#[test]
fn the_tool_load_options_admit_an_input_the_library_defaults_refuse() {
    let temp = workdir();
    let (input, bytes) = instrument_scale_input(temp.path(), 20_000);
    // 42,776,703 bytes: the size a converter writes for 20,000 sparse profile
    // spectra, and four times the largest fixture in this suite.
    assert!((40..44).contains(&(bytes / (1 << 20))), "{bytes} bytes");

    let read = |options: &ReadOptions| {
        FileHandler::load_experiment_with_read_options(
            &input,
            &[FileType::MzMl],
            &PeakFileOptions::default(),
            options,
        )
    };
    // The library default refuses the dangling instrument `softwareRef`, which
    // source `MzMLHandler` drops; a tool path accepts it (decision D10).
    let refused = read(&ReadOptions::default()).unwrap_err();
    assert!(
        refused.to_string().contains("unresolved softwareRef"),
        "{refused}"
    );
    // The fixed ceilings this tool shipped with refuse the size, with the
    // source-compatibility switches unchanged: only the scaling differs. The
    // floor reached first here is the cumulative parameter storage (512 MiB),
    // not the 10,000,000-point one, because a document that declares ten
    // million points must carry them: that is 120 MB of binary before base64,
    // which is why the point ceiling is measured on the real input instead.
    let fixed = ReadOptions {
        scaling: InputScaling::default().fixed(),
        ..PeakPickerHiRes::read_options()
    };
    let refused = read(&fixed).unwrap_err();
    assert!(
        refused
            .to_string()
            .contains("parameter bytes exceed configured limit"),
        "{refused}"
    );

    // The tool's own options read every spectrum.
    let raw = read(&PeakPickerHiRes::read_options()).unwrap();
    assert_eq!(raw.spectra.len(), 20_000);
    assert_eq!(raw.spectra[19_999].peaks.len(), 15);

    // And the tool picks the whole document, one centroid per profile peak.
    let out = temp.path().join("instrument_scale_picked.tmp.mzML");
    let outcome = run(&["-test", "-in", &text(&input), "-out", &text(&out)]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_eq!(
        outcome.out,
        "#Spectra that needed to and could be picked by MS-level:\n  MS-level 1: 20000 / 20000\n"
    );
    let picked = load(&out);
    assert_eq!(picked.spectra.len(), 20_000);
    for spectrum in [&picked.spectra[0], &picked.spectra[19_999]] {
        let centroids: Vec<f64> = spectrum.peaks.iter().map(|peak| peak.mz).collect();
        assert_eq!(centroids.len(), 3);
        for (centroid, expected) in centroids.iter().zip([400.0, 410.0, 420.0]) {
            assert!((centroid - expected).abs() < 1e-3, "{centroid} {expected}");
        }
    }
}

/// The new fixtures this round adds, in `tests/data/peak_picking/`.
fn close_fixture(name: &str) -> PathBuf {
    data("peak_picking").join(name)
}

/// The record start tags of an mzML document, in order.
fn start_tags(text: &str) -> Vec<String> {
    text.match_indices("<spectrum id=")
        .map(|(at, _)| {
            let rest = &text[at..];
            rest[..rest.find('>').expect("a start tag")].to_owned()
        })
        .collect()
}

/// Native difference 12, now reproduced rather than refused: the low-memory
/// mode writes the source's dangling header references.
///
/// Only the first record reaches the header, so a later record's
/// `dataProcessing` cannot be numbered against it. The source numbers the
/// reference by the record's position in the stream instead
/// (`MzMLHandler.cpp:5258-5272`, `dps_` holding one entry), and a
/// `sourceFileRef` by the same number for **every** record after the first
/// that carries one (`:5252-5255`, which never consults the header).
/// [`ReferencePolicy::SourceDangling`] reproduces both, because refusing would
/// stop the mode on any `FileMerger` output.
///
/// Executed against the C++ Release build on `ibminode06` on this exact
/// fixture (`../oracle/p4-lowmemory`, `logs/closediff2_06.log` section A and
/// `logs/closediff3_06.log` section A). The C++ low-memory output's five start
/// tags carry, in order:
///
/// | record | `sourceFileRef` | `dataProcessingRef` |
/// | --- | --- | --- |
/// | 0 | `sf_sp_0`, declared | `dp_sp_0`, declared |
/// | 1 | `sf_sp_1` | `dp_sp_1` |
/// | 2 | `sf_sp_2` | `dp_sp_2` |
/// | 3 | `sf_sp_3` | none |
/// | 4 | `sf_sp_4` | none |
///
/// against a header declaring one record `sourceFile` and one
/// `dataProcessing` — six dangling references. This port writes the same
/// references for records 1 to 4, in the source's own spelling; the first
/// record's point at the header entries this writer's own scheme names.
#[test]
fn the_low_memory_mode_writes_the_sources_dangling_references() {
    let input = close_fixture("PeakPickerHiRes_refs_input.mzML");
    let (outcome, bytes, produced_at) = low_memory(None, &input, &[]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let written = String::from_utf8(bytes).unwrap();
    let tags = start_tags(&written);
    assert_eq!(tags.len(), 5, "{written:.800}");
    // Record 0 establishes the header, so its references are the ones the
    // header declares and it needs no `dataProcessingRef` of its own.
    assert!(
        tags[0].contains(" sourceFileRef=\"sf_00000000000000000003\""),
        "{}",
        tags[0]
    );
    assert!(!tags[0].contains("dataProcessingRef="), "{}", tags[0]);
    for (index, processing) in [(1usize, true), (2, true), (3, false), (4, false)] {
        let tag = &tags[index];
        assert!(
            tag.contains(&format!(" sourceFileRef=\"sf_sp_{index}\"")),
            "{tag}"
        );
        assert_eq!(
            tag.contains(&format!(" dataProcessingRef=\"dp_sp_{index}\"")),
            processing,
            "{tag}"
        );
        assert_eq!(tag.contains("dataProcessingRef="), processing, "{tag}");
    }
    // Every one of those identifiers dangles: nothing declares it.
    for dangling in [
        "sf_sp_1", "sf_sp_2", "sf_sp_3", "sf_sp_4", "dp_sp_1", "dp_sp_2",
    ] {
        assert!(
            !written.contains(&format!(" id=\"{dangling}\"")),
            "{dangling} is declared"
        );
    }
    // The document is complete, and the records carry the peaks the in-memory
    // mode produces — compared through the text, which is the stronger
    // comparison here because it also covers the encoding, not only the
    // decoded values the round trip below checks.
    assert!(written.ends_with("</indexedmzML>\n"));
    let temp = workdir();
    let out = temp.path().join("refs_in_memory.tmp.mzML");
    let in_memory = run(&["-test", "-in", &text(&input), "-out", &text(&out)]);
    assert_eq!(in_memory.code, ExitCode::ExecutionOk, "{}", in_memory.err);
    let mem = std::fs::read_to_string(&out).unwrap();
    assert_eq!(
        binaries(&written),
        binaries(&mem),
        "the encoded arrays differ"
    );

    // And the tool reads its own output back. Until decision D14 the reader
    // refused an unregistered `sourceFileRef` under either policy, so the port
    // wrote, on this fixture, a file neither it nor the C++ reader's strict
    // equivalent would take — while the C++ reader took both its own output
    // and this one. `PeakPickerHiRes::read_options()` is
    // `mzml::ReadOptions::source()`, which now covers `sourceFileRef` as it
    // covers `dataProcessingRef`, and the records survive the round trip.
    // Measured on `ibminode06` in both directions
    // (`../oracle/reader-roundtrip/logs/roundtrip_06.log`, section `refs`):
    // each implementation's `FileInfo` exits 0 on the other's low-memory
    // output, and reading either file back and writing it out again gives the
    // same decoded content, `902f49fd94e5a4a8`, in all four combinations.
    let back = FileHandler::load_experiment_with_read_options(
        &produced_at.out,
        &[FileType::MzMl],
        &PeakFileOptions::default(),
        &PeakPickerHiRes::read_options(),
    )
    .unwrap();
    assert_eq!(back.spectra.len(), 5);
    // The first record's source file is the declared one; the four renumbered
    // references named nothing and are dropped, as `MzMLHandler.cpp:896-906`
    // drops them.
    assert_eq!(back.spectra[0].source_file.name, "part_one.mzML");
    assert!(
        back.spectra[1..]
            .iter()
            .all(|s| s.source_file == Default::default())
    );
    // The strict default still refuses the file, which is what makes the
    // reference genuinely dangling rather than merely unusual.
    let error = FileHandler::load_experiment_with_read_options(
        &produced_at.out,
        &[FileType::MzMl],
        &PeakFileOptions::default(),
        &ReadOptions::default(),
    )
    .unwrap_err();
    assert!(error.to_string().contains("unresolved"), "{error}");
}

/// The pointer rule, at the tool level: two textually identical
/// `dataProcessing` entries under different identifiers are two histories.
///
/// `PeakPickerHiRes_dupdp_input.mzML` is `PeakPickerHiRes_refs_input.mzML`
/// with `dp_sp_1`'s `softwareRef` repointed from `so_dp_1` to `so_dp_0`, so
/// that `dp_sp_0` and `dp_sp_1` render identically and only their identifiers
/// differ. The source compares `spec.getDataProcessing() != dps[0]`
/// (`MzMLHandler.cpp:5258`) by pointer, so it still writes `dp_sp_1` on record
/// 1 and `dp_sp_2` on record 2.
///
/// Measured on `ibminode06` against the Release build at the pins
/// (`../oracle/reader-roundtrip/logs/roundtrip_06.log`, section `dupdp`): the
/// C++ low-memory output of this fixture is **byte-identical** to its
/// low-memory output of the unmodified `refs` fixture, sha256
/// `7c75908440e1c154…`. This port's two low-memory outputs are likewise
/// byte-identical to each other, sha256 `6ef28364d62ea7a7…`, and carry the
/// same six dangling identifiers as the C++ output. Before this round the port
/// compared the rendered declaration text instead and wrote no
/// `dataProcessingRef` at all here; that was native difference 12, and it is
/// closed.
#[test]
fn the_low_memory_mode_decides_a_duplicate_history_by_pointer() {
    let refs = low_memory(None, &close_fixture("PeakPickerHiRes_refs_input.mzML"), &[]);
    let dupdp = low_memory(
        None,
        &close_fixture("PeakPickerHiRes_dupdp_input.mzML"),
        &[],
    );
    assert_eq!(refs.0.code, ExitCode::ExecutionOk, "{}", refs.0.err);
    assert_eq!(dupdp.0.code, ExitCode::ExecutionOk, "{}", dupdp.0.err);
    // The two inputs differ, in exactly one attribute value.
    assert_ne!(
        std::fs::read(close_fixture("PeakPickerHiRes_refs_input.mzML")).unwrap(),
        std::fs::read(close_fixture("PeakPickerHiRes_dupdp_input.mzML")).unwrap()
    );
    // The outputs do not, on either side.
    assert_eq!(refs.1, dupdp.1);
    let written = String::from_utf8(dupdp.1).unwrap();
    let tags = start_tags(&written);
    assert_eq!(tags.len(), 5, "{written:.800}");
    assert!(
        tags[1].contains(" dataProcessingRef=\"dp_sp_1\""),
        "{}",
        tags[1]
    );
    assert!(
        tags[2].contains(" dataProcessingRef=\"dp_sp_2\""),
        "{}",
        tags[2]
    );
    assert!(!tags[3].contains("dataProcessingRef="), "{}", tags[3]);
    assert!(!tags[4].contains("dataProcessingRef="), "{}", tags[4]);
    for dangling in ["dp_sp_1", "dp_sp_2"] {
        assert!(
            !written.contains(&format!(" id=\"{dangling}\"")),
            "{dangling} is declared"
        );
    }
}

/// Every `<binary>` payload of a document, in order: the encoded arrays, which
/// are the same in both modes whatever the header around them says.
fn binaries(text: &str) -> Vec<&str> {
    text.match_indices("<binary>")
        .map(|(at, _)| {
            let rest = &text[at + "<binary>".len()..];
            &rest[..rest.find("</binary>").expect("a closed array")]
        })
        .collect()
}

/// Where a low-memory run fails, it leaves the batches it had already sent.
///
/// Both implementations read in batches of `max_data_pool_size`, the 100 of
/// `PeakFileOptions.h:248`, and neither counting pass reads record contents —
/// the source's runs with `LD_RAWCOUNTS` and `skip_spectrum_`
/// (`MzMLHandler.cpp:966-974`), this port's sets `state.raw` and a
/// `skip_depth` at the list tag (`src/format/mzml_counts.rs:859`, `:465-469`).
/// A record that is well-formed XML but wrong inside is therefore discovered
/// only in the second pass, with `floor(index / 100) * 100` records already
/// written, and what stays on disc is a closed, indexed, **reloadable**
/// document announcing the count the first pass declared.
///
/// Executed on `ibminode06` on this fixture with four corruptions at index 4
/// and at index 104 (`../oracle/p4-lowmemory`, `logs/closediff2_06.log`
/// section D and `logs/closediff3_06.log` section C). Malformed base64 fails
/// on both sides: the C++ low-memory run leaves 0 bytes at index 4 and 100
/// records at index 104. The other three — a `defaultArrayLength` one too
/// large, a non-numeric `scan start time`, two records sharing a native id —
/// the C++ accepts, exiting 0 with all 110 records in **both** modes, while
/// this port's reader refuses them; so at index 104 this port leaves a
/// truncated 100-record document where the C++ writes the file in full.
#[test]
fn a_low_memory_failure_leaves_the_batches_already_written() {
    let batches = close_fixture("PeakPickerHiRes_batches_input.mzML");
    let corrupt = |index: usize, rule: &dyn Fn(&str) -> String| -> (PathBuf, TempDir) {
        derived(&format!("batch_{index}.mzML"), &batches, |s| {
            let at = s
                .match_indices("<spectrum ")
                .nth(index)
                .expect("the record")
                .0;
            let end = at + s[at..].find("</spectrum>").expect("the record end");
            format!("{}{}{}", &s[..at], rule(&s[at..end]), &s[end..])
        })
    };
    let bad_base64 = |record: &str| {
        let at = record.find("<binary>").expect("an array") + "<binary>".len();
        format!("{}!!!{}", &record[..at], &record[at..])
    };
    let duplicate_id = |record: &str| {
        let at = record.find("id=\"").expect("the id") + "id=\"".len();
        let end = at + record[at..].find('"').expect("the id end");
        format!("{}spectrum=0{}", &record[..at], &record[end..])
    };

    for (index, records) in [(4usize, 0usize), (104, 100)] {
        for (kind, rule) in [
            ("base64", &bad_base64 as &dyn Fn(&str) -> String),
            ("duplicate id", &duplicate_id),
        ] {
            let (input, _temp) = corrupt(index, rule);
            let (outcome, bytes, produced_at) = low_memory(None, &input, &[]);
            let at = format!("{kind} at {index}");
            assert_eq!(
                outcome.code,
                ExitCode::InputFileCorrupt,
                "{at}: {}",
                outcome.err
            );
            assert!(
                outcome.err.starts_with("Error: Unable to read file ("),
                "{at}: {}",
                outcome.err
            );
            let written = String::from_utf8(bytes).unwrap();
            assert_eq!(start_tags(&written).len(), records, "{at}");
            if records == 0 {
                // Nothing reached the writer, so there is no document at all.
                assert!(written.is_empty(), "{at}");
                continue;
            }
            // Closed, indexed, and announcing the whole input's count.
            assert!(written.contains("<spectrumList count=\"110\""), "{at}");
            assert!(written.contains("</spectrumList>\n</run></mzML>\n"), "{at}");
            assert!(written.ends_with("</indexedmzML>\n"), "{at}");
            assert_eq!(load(&produced_at.out).spectra.len(), records, "{at}");

            // The in-memory mode leaves no file behind on the same input.
            let temp = workdir();
            let out = temp.path().join("batch_in_memory.tmp.mzML");
            let in_memory = run(&["-test", "-in", &text(&input), "-out", &text(&out)]);
            assert_eq!(in_memory.code, ExitCode::InputFileCorrupt, "{at}");
            assert!(!out.exists(), "{at}");
        }
    }

    // The healthy fixture crosses the same boundary with both modes agreeing.
    let (outcome, bytes, _at) = low_memory(None, &batches, &[]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_eq!(
        start_tags(&String::from_utf8(bytes.clone()).unwrap()).len(),
        110
    );
    assert_eq!(bytes, in_memory_bytes(None, &batches, &[]));
}

/// An `-out` that names an existing directory: the one measured case where the
/// C++ low-memory run reports success on a run that produced nothing.
///
/// `MSDataWritingConsumer`'s constructor never checks its `std::ofstream`
/// (`MSDataWritingConsumer.cpp:33`), and `doLowMemAlgorithm` returns
/// `EXECUTION_OK` regardless, so the C++ low-memory run exits 0 with an empty
/// standard error having written nothing, where its in-memory run exits 5
/// `Error: Unable to write file`. Measured on `ibminode06`
/// (`../oracle/p4-lowmemory`, `logs/closediff1_06.log` section F,
/// `logs/closediff3_06.log` section D), with the two controls that do **not**
/// reach the consumer — a read-only `-out` and an `-out` under a missing
/// directory — exiting 5 with `Cannot write output file given from parameter
/// '-out'!` in both implementations and both modes.
///
/// This port reports the operating system's refusal instead, in both modes.
#[test]
fn an_out_that_names_a_directory_is_reported_in_both_modes() {
    let temp = workdir();
    let out = temp.path().join("isdir.mzML");
    std::fs::create_dir(&out).unwrap();
    let input = text(&workflow_input(6));
    for extra in [&["-processOption", "lowmemory"][..], &[]] {
        let mut args = vec!["-test", "-in", &input, "-out"];
        let out_text = text(&out);
        args.push(&out_text);
        args.extend_from_slice(extra);
        let outcome = run(&args);
        assert_eq!(
            outcome.code,
            ExitCode::UnknownError,
            "{:?}: {}",
            extra,
            outcome.err
        );
        assert!(
            outcome
                .err
                .starts_with("Error: Unexpected internal error (")
                && outcome.err.contains("directory"),
            "{extra:?}: {}",
            outcome.err
        );
        // Nothing was written into the directory either way.
        assert_eq!(std::fs::read_dir(&out).unwrap().count(), 0, "{extra:?}");
    }
}
