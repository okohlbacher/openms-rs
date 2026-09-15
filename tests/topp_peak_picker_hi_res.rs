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
use openms::kernel::MSExperiment;
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
    for threads in ["1", "16", "0"] {
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
#[test]
fn threads_do_not_change_the_workflow_1_output() {
    let temp = workdir();
    let ini = text(&fixture("PeakPickerHiRes_parameters.ini"));
    let input = text(&workflow_input(1));
    let mut bytes: Vec<Vec<u8>> = Vec::new();
    for threads in ["1", "2", "16", "0"] {
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

/// `-processOption lowmemory` is refused explicitly until package P4 ports the
/// source's transforming consumer, and nothing is written.
#[test]
fn low_memory_processing_is_refused_explicitly() {
    let temp = workdir();
    let out = temp.path().join("PeakPickerHiRes_3.tmp.mzML");
    let outcome = run(&[
        "-test",
        "-ini",
        &text(&fixture("PeakPickerHiRes_parameters.ini")),
        "-in",
        &text(&workflow_input(1)),
        "-out",
        &text(&out),
        "-processOption",
        "lowmemory",
    ]);
    assert_eq!(outcome.code, ExitCode::IncompatibleInputData);
    assert!(
        outcome.err.contains(
            "Error: unsupported: PeakPickerHiRes -processOption lowmemory is not ported yet"
        ),
        "{}",
        outcome.err
    );
    assert!(!out.exists());
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
