// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Executable specification of the TOPPBase lifecycle: package C4-CLI-SPEC,
//! turned green by the CLI closure W1.4.
//!
//! Two kinds of evidence back these cases.
//!
//! * **Oracle cases** run DTAExtractor and SpectraFilterWindowMower exactly as
//!   `../oracle/topp-cli-lifecycle/run.sh` ran the C++ product SDK (core
//!   4fdec46, Debug build). Each asserts that run's exit code and a diagnostic
//!   line from its output; the case name quoted in each comment is the one in
//!   the oracle's `manifest.json`. Evidence: oracle-generated (tier 1 executed
//!   differential). None of these cases reaches a Debug-only precondition.
//! * **Upstream class-test cases** transcribe the `TOPPBase_test.cpp` sections
//!   (cli c19e494) that apply to the native API, with their literals unchanged
//!   (tier 3 source review). Sections that cannot apply are named in
//!   `docs/TOPP_CLI_SUPPORT.md`.
//!
//! A few native cases cover services the source has no test for: run-phase
//! error mapping, `run_io` stream routing and the context services a tool uses.
//!
//! Every case runs in its own temporary directory, because the tests run in
//! parallel. `tests/data/topp_cli_lifecycle/topp_cli_lifecycle_provenance.json`
//! records the fixtures, their hashes and the oracle manifest.

// The TOPP framework lives behind `paramxml` (every tool supports -ini) and the
// stand-in tools read mzML, so the whole file is inert without both features.
#![cfg(all(feature = "mzml", feature = "paramxml"))]

use openms::cli::tools::{DTAExtractor, SpectraFilterWindowMower};
use openms::cli::{
    ExitCode, MAX_ARGUMENTS, TEST_MODE_COMPLETION_TIME, TEST_MODE_PARAMETER_KEY,
    TEST_MODE_PARAMETER_VALUE, TEST_MODE_UNIQUE_ID_SEED, TEST_MODE_VERSION, Tool, ToolContext,
    ToolSpec, input_file_readable, output_file_writable, parse_range, run_with,
};
use openms::concept::parallel::Threads;
use openms::concept::progress_logger::ProgressLogType;
use openms::data_structures::DateTime;
use openms::format::file_handler::FileHandler;
use openms::format::file_types::FileType;
use openms::kernel::{
    ColumnHeader, ConsensusMap, FeatureMap, MSChromatogram, MSExperiment, MSSpectrum,
};
use openms::metadata::{MetaValue, ProcessingAction};
use openms::param::{Param, ParamValue};
use openms::{Error, Result};
use std::cell::RefCell;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn fixture(name: &str) -> PathBuf {
    Path::new("tests/data").join(name)
}
fn lifecycle(name: &str) -> String {
    text(Path::new("tests/data/topp_cli_lifecycle").join(name))
}
fn text(path: impl AsRef<Path>) -> String {
    path.as_ref().to_string_lossy().into_owned()
}
fn swm_input() -> String {
    text(fixture("window_mower_tool_input.mzML"))
}

struct Outcome {
    code: ExitCode,
    out: String,
    err: String,
}

fn run_arguments<T: Tool>(arguments: &[String]) -> Outcome {
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<T>(arguments, &mut out, &mut err);
    Outcome {
        code,
        out: String::from_utf8_lossy(&out).into_owned(),
        err: String::from_utf8_lossy(&err).into_owned(),
    }
}

/// Run a tool as `<tool> args...`.
fn run<T: Tool>(args: &[&str]) -> Outcome {
    let arguments: Vec<String> = std::iter::once(T::NAME.to_owned())
        .chain(args.iter().map(|a| (*a).to_owned()))
        .collect();
    run_arguments::<T>(&arguments)
}

/// A fresh directory per case, removed when the case ends.
struct Workdir(PathBuf);
impl Workdir {
    fn new(case: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "openms-cli-lifecycle-{}-{case}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }
    fn file(&self, name: &str) -> String {
        text(self.0.join(name))
    }
}
impl Drop for Workdir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn load(path: impl AsRef<Path>) -> MSExperiment {
    FileHandler::load_experiment(path, &[FileType::MzMl]).unwrap()
}

/// Canonical peak comparison with the tolerances of the retained
/// `TOPP_SpectraFilterWindowMower_1` comparison.
fn assert_same_peaks(produced: impl AsRef<Path>, expected: impl AsRef<Path>) {
    let produced = load(produced);
    let expected = load(expected);
    assert_eq!(produced.spectra.len(), expected.spectra.len());
    for (index, (a, e)) in produced.spectra.iter().zip(&expected.spectra).enumerate() {
        assert_eq!(a.peaks.len(), e.peaks.len(), "spectrum {index} peak count");
        for (i, (pa, pe)) in a.peaks.iter().zip(&e.peaks).enumerate() {
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
}

fn total_peaks(path: impl AsRef<Path>) -> usize {
    load(path).spectra.iter().map(|s| s.peaks.len()).sum()
}

thread_local! {
    /// The context a test tool's body last received on this thread.
    static CAPTURED: RefCell<Option<ToolContext>> = const { RefCell::new(None) };
}
fn capture(ctx: &ToolContext) {
    CAPTURED.with(|slot| *slot.borrow_mut() = Some(ctx.clone()));
}
fn captured() -> ToolContext {
    CAPTURED
        .with(|slot| slot.borrow_mut().take())
        .expect("the tool body ran")
}

// ---------------------------------------------------------------------------
// Upstream test tools (TOPPBase_test.cpp:33-377)
// ---------------------------------------------------------------------------

/// `TOPPBaseTest`: optional parameters of every type.
struct ToppBaseTest;
impl Tool for ToppBaseTest {
    const NAME: &'static str = "TOPPBaseTest";
    const DESCRIPTION: &'static str = "A test class";
    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_string_option(
            "stringoption",
            "<string>",
            "string default",
            "string description",
            false,
            false,
        )?;
        spec.register_int_option("intoption", "<int>", 4711, "int description", false, false)?;
        spec.register_double_option(
            "doubleoption",
            "<double>",
            0.4711,
            "double description",
            false,
            false,
        )?;
        spec.register_int_list(
            "intlist",
            "<intlist>",
            &[1, 2, 3, 4],
            "intlist description",
            false,
            false,
        )?;
        spec.register_double_list(
            "doublelist",
            "<doublelist>",
            &[0.4711, 1.022, 4.0],
            "doubelist description",
            false,
            false,
        )?;
        spec.register_string_list(
            "stringlist",
            "<stringlist>",
            &["abc", "def", "ghi", "jkl"],
            "stringlist description",
            false,
            false,
        )?;
        spec.register_flag("flag", "flag description", false)?;
        spec.register_string_list(
            "stringlist2",
            "<stringlist>",
            &["hopla", "dude"],
            "stringlist with restrictions",
            false,
            false,
        )?;
        spec.set_valid_strings("stringlist2", &["hopla", "dude"])?;
        spec.register_int_list(
            "intlist2",
            "<int>",
            &[3, 4, 5],
            "intlist with restrictions",
            false,
            false,
        )?;
        spec.set_min_int("intlist2", 2)?;
        spec.set_max_int("intlist2", 6)?;
        spec.register_double_list(
            "doublelist2",
            "<double>",
            &[1.2, 2.33],
            "doublelist with restrictions",
            false,
            false,
        )?;
        spec.set_min_float("doublelist2", 0.2)?;
        spec.set_max_float("doublelist2", 5.4)?;
        Ok(())
    }
    fn run(ctx: &ToolContext) -> Result<ExitCode> {
        capture(ctx);
        Ok(ExitCode::ExecutionOk)
    }
}

/// `TOPPBaseTestNOP`: required parameters. `registerStringOption_` and the
/// list registrations default to `required = true` in `TOPPBase.h`.
struct ToppBaseTestNop;
impl Tool for ToppBaseTestNop {
    const NAME: &'static str = "TOPPBaseTestNOP";
    const DESCRIPTION: &'static str = "A test class with non-optional parameters";
    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_string_option(
            "stringoption",
            "<string>",
            "",
            "string description",
            true,
            false,
        )?;
        spec.register_int_option("intoption", "<int>", 0, "int description", false, false)?;
        spec.register_double_option(
            "doubleoption",
            "<double>",
            -1.0,
            "double description",
            false,
            false,
        )?;
        spec.register_flag("flag", "flag description", false)?;
        spec.register_string_list(
            "stringlist",
            "<stringlist>",
            &[],
            "stringlist description",
            true,
            false,
        )?;
        spec.register_int_list(
            "intlist",
            "<intlist>",
            &[],
            "intlist description",
            true,
            false,
        )?;
        spec.register_double_list(
            "doublelist",
            "<doublelist>",
            &[],
            "doubelist description",
            true,
            false,
        )?;
        Ok(())
    }
    fn run(ctx: &ToolContext) -> Result<ExitCode> {
        capture(ctx);
        Ok(ExitCode::ExecutionOk)
    }
}

/// `TOPPBaseTestParam`: parameters registered from a whole `Param`.
struct ToppBaseTestParam;
fn full_param() -> Param {
    let mut param = Param::new();
    param
        .set_value(
            "param_int",
            ParamValue::Integer(123),
            "param int description",
            &[],
        )
        .unwrap();
    param
        .set_value(
            "param_double",
            ParamValue::Float(-4.56),
            "param double description",
            &[],
        )
        .unwrap();
    param
        .set_value(
            "param_string",
            ParamValue::String("test".into()),
            "param string description",
            &[],
        )
        .unwrap();
    param
        .set_value(
            "param_stringlist",
            ParamValue::StringList(vec!["this".into(), "is".into(), "a".into(), "test".into()]),
            "param stringlist description",
            &[],
        )
        .unwrap();
    param
        .set_value(
            "param_intlist",
            ParamValue::IntegerList(vec![7, -8, 9]),
            "param intlist description",
            &[],
        )
        .unwrap();
    param
        .set_value(
            "param_doublelist",
            ParamValue::FloatList(vec![123.0, -4.56, 0.789]),
            "param doublelist description",
            &[],
        )
        .unwrap();
    param
        .set_value(
            "param_flag",
            ParamValue::String("true".into()),
            "param flag description",
            &[],
        )
        .unwrap();
    param
        .set_valid_strings("param_flag", &["true".into(), "false".into()])
        .unwrap();
    param
}
impl Tool for ToppBaseTestParam {
    const NAME: &'static str = "TOPPBaseTestParam";
    const DESCRIPTION: &'static str = "A test class with parameters derived from Param";
    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_full_param(&full_param())
    }
    fn run(ctx: &ToolContext) -> Result<ExitCode> {
        capture(ctx);
        Ok(ExitCode::ExecutionOk)
    }
}

/// `TOPPBaseCmdParseTest`: no parameters of its own.
struct ToppBaseCmdParseTest;
impl Tool for ToppBaseCmdParseTest {
    const NAME: &'static str = "TOPPBaseCmdParseTest";
    const DESCRIPTION: &'static str = "A test class to test parts of the cmd parser functionality";
    fn register(_spec: &mut ToolSpec) -> Result<()> {
        Ok(())
    }
    fn run(_ctx: &ToolContext) -> Result<ExitCode> {
        Ok(ExitCode::ExecutionOk)
    }
}

/// `TOPPBaseCmdParseSubsectionsTest`: two algorithm subsections.
struct ToppBaseCmdParseSubsectionsTest;
impl Tool for ToppBaseCmdParseSubsectionsTest {
    const NAME: &'static str = "TOPPBaseCmdParseSubsectionsTest";
    const DESCRIPTION: &'static str = "A test class to test parts of the cmd parser functionality";
    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_string_option(
            "stringoption",
            "<string>",
            "",
            "string description",
            true,
            false,
        )?;
        spec.register_subsection("algorithm", "Algorithm parameters section")?;
        spec.register_subsection("other", "Other parameters section")?;
        Ok(())
    }
    fn subsection_defaults(section: &str) -> Result<Option<Param>> {
        let value = |text: &str| ParamValue::String(text.to_owned());
        let flag_strings = ["true".to_owned(), "false".to_owned()];
        let mut param = Param::new();
        if section == "algorithm" {
            param.set_value("param1", value("param1_value"), "param1_description", &[])?;
            param.set_value("param2", value("param2_value"), "param2_description", &[])?;
        } else {
            param.set_value("param3", value("param3_value"), "param3_description", &[])?;
            param.set_value("param4", value("param4_value"), "param4_description", &[])?;
            param.set_value("flagparam", value("false"), "this will be a flag", &[])?;
            param.set_valid_strings("flagparam", &flag_strings)?;
            param.set_value(
                "nonflagparam",
                value("true"),
                "this will be a string param with true/false",
                &[],
            )?;
            param.set_valid_strings("nonflagparam", &flag_strings)?;
        }
        Ok(Some(param))
    }
    fn run(ctx: &ToolContext) -> Result<ExitCode> {
        capture(ctx);
        Ok(ExitCode::ExecutionOk)
    }
}

fn string_value(ctx: &ToolContext, key: &str) -> String {
    ctx.param().value(key).unwrap().as_str().unwrap().to_owned()
}

// ---------------------------------------------------------------------------
// Oracle cases: the command-line parse closure (TOPPBase.cpp:182-255, 2288-2470)
// ---------------------------------------------------------------------------

/// Oracle `no_arguments`: exit 6, "No options given. Aborting!", usage on stderr.
#[test]
#[ignore = "INT-W1.4 d: a bare invocation exits 6 at TOPPBase.cpp:227-232 (oracle no_arguments), but tests/topp_dta_extractor.rs:131, tests/topp_baseline_filter.rs:93 and tests/topp_map_normalizer.rs:102 still assert MISSING_PARAMETERS and package CLI-1 may not change them"]
fn a_bare_invocation_is_refused() {
    let outcome = run::<DTAExtractor>(&[]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome.err.contains("No options given. Aborting!"),
        "{}",
        outcome.err
    );
    assert!(outcome.err.contains("DTAExtractor --"), "{}", outcome.err);
}

/// Oracle `unknown_option`: exit 6, usage on stderr.
#[test]
fn an_unknown_option_is_refused() {
    let dir = Workdir::new("unknown-option");
    let outcome = run::<SpectraFilterWindowMower>(&[
        "-in",
        &swm_input(),
        "-out",
        &dir.file("out.mzML"),
        "-bogus",
        "1",
    ]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome
            .err
            .contains("Unknown option(s) '[-bogus]' given. Aborting!"),
        "{}",
        outcome.err
    );
    assert!(
        outcome.err.contains("SpectraFilterWindowMower --"),
        "{}",
        outcome.err
    );
    assert!(outcome.out.is_empty(), "{}", outcome.out);
}

/// Oracle `trailing_text`: exit 6.
#[test]
fn trailing_text_is_refused() {
    let dir = Workdir::new("trailing-text");
    let outcome = run::<SpectraFilterWindowMower>(&[
        "-in",
        &swm_input(),
        "-out",
        &dir.file("out.mzML"),
        "extra",
    ]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome
            .err
            .contains("Trailing text argument(s) '[extra]' given. Aborting!"),
        "{}",
        outcome.err
    );
}

/// Oracle `flag_with_trailing_text`: a flag followed by text is a parse error.
#[test]
fn a_flag_followed_by_text_is_refused() {
    let dir = Workdir::new("flag-text");
    let outcome = run::<SpectraFilterWindowMower>(&[
        "-in",
        &swm_input(),
        "-out",
        &dir.file("out.mzML"),
        "-test",
        "extra",
    ]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome.err.contains(
            "Invalid parameter values (InvalidParameter): Command line error: Trailing arguments after flag '-test': extra. Aborting!"
        ),
        "{}",
        outcome.err
    );
}

/// Oracle `int_option_not_numeric`: conversion happens while parsing.
#[test]
fn a_non_numeric_integer_is_refused() {
    let dir = Workdir::new("int-not-numeric");
    let outcome = run::<SpectraFilterWindowMower>(&[
        "-in",
        &swm_input(),
        "-out",
        &dir.file("out.mzML"),
        "-threads",
        "abc",
    ]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome.err.contains(
            "Invalid parameter values (ConversionError): Could not convert string 'abc' to an integer value. Aborting!"
        ),
        "{}",
        outcome.err
    );
}

/// Oracle `instance_on_command_line`: `-instance` is registered but left out of
/// the defaults (TOPPBase.cpp:2104), so the strict update rejects it.
#[test]
fn instance_on_the_command_line_is_refused() {
    let dir = Workdir::new("instance");
    let outcome = run::<SpectraFilterWindowMower>(&[
        "-test",
        "-instance",
        "1",
        "-in",
        &swm_input(),
        "-out",
        &dir.file("out.mzML"),
    ]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome.err.contains(
            "Unknown (or deprecated) Parameter 'instance' given in outdated parameter file!"
        ),
        "{}",
        outcome.err
    );
    assert!(
        outcome
            .err
            .contains("Parameters passed to 'SpectraFilterWindowMower' are invalid."),
        "{}",
        outcome.err
    );
}

// ---------------------------------------------------------------------------
// Oracle cases: file checks (TOPPBase.cpp:1968-2013) and run-phase errors
// ---------------------------------------------------------------------------

/// Oracle `missing_input`: exit 1.
#[test]
fn a_missing_input_is_input_file_not_found() {
    let dir = Workdir::new("missing-input");
    let outcome = run::<SpectraFilterWindowMower>(&[
        "-in",
        &dir.file("absent.mzML"),
        "-out",
        &dir.file("out.mzML"),
    ]);
    assert_eq!(outcome.code, ExitCode::InputFileNotFound, "{}", outcome.err);
    assert!(
        outcome
            .err
            .contains("Cannot read input file given from parameter '-in'!"),
        "{}",
        outcome.err
    );
    assert!(
        outcome.err.contains("Error: File not found"),
        "{}",
        outcome.err
    );
}

/// Oracle `zero_byte_input`: exit 4, not a parse error.
#[test]
fn a_zero_byte_input_is_input_file_empty() {
    let dir = Workdir::new("zero-byte");
    let empty = dir.file("empty.mzML");
    fs::write(&empty, b"").unwrap();
    let outcome = run::<SpectraFilterWindowMower>(&["-in", &empty, "-out", &dir.file("out.mzML")]);
    assert_eq!(outcome.code, ExitCode::InputFileEmpty, "{}", outcome.err);
    assert!(
        outcome
            .err
            .contains("Cannot read input file given from parameter '-in'!"),
        "{}",
        outcome.err
    );
    assert!(outcome.err.contains("Error: File empty"), "{}", outcome.err);
}

/// Oracle `unreadable_input`: exit 2. Skipped when permissions do not restrict
/// the process, as for root.
#[cfg(unix)]
#[test]
fn an_unreadable_input_is_input_file_not_readable() {
    use std::os::unix::fs::PermissionsExt;
    let dir = Workdir::new("unreadable");
    let path = dir.file("unreadable.mzML");
    fs::copy(fixture("window_mower_tool_input.mzML"), &path).unwrap();
    fs::set_permissions(&path, fs::Permissions::from_mode(0o000)).unwrap();
    if fs::File::open(&path).is_ok() {
        fs::set_permissions(&path, fs::Permissions::from_mode(0o644)).unwrap();
        return;
    }
    let outcome = run::<SpectraFilterWindowMower>(&["-in", &path, "-out", &dir.file("out.mzML")]);
    fs::set_permissions(&path, fs::Permissions::from_mode(0o644)).unwrap();
    assert_eq!(
        outcome.code,
        ExitCode::InputFileNotReadable,
        "{}",
        outcome.err
    );
    assert!(
        outcome.err.contains("Error: File not readable"),
        "{}",
        outcome.err
    );
}

/// Oracle `unwritable_output`: exit 5.
#[test]
fn an_unwritable_output_is_cannot_write_output_file() {
    let outcome =
        run::<SpectraFilterWindowMower>(&["-in", &swm_input(), "-out", "/nonexistent_dir/x.mzML"]);
    assert_eq!(
        outcome.code,
        ExitCode::CannotWriteOutputFile,
        "{}",
        outcome.err
    );
    assert!(
        outcome
            .err
            .contains("Cannot write output file given from parameter '-out'!"),
        "{}",
        outcome.err
    );
    assert!(
        outcome.err.contains("Error: Unable to write file"),
        "{}",
        outcome.err
    );
}

/// Oracle `output_wrong_extension`: exit 6.
#[test]
fn a_wrong_output_extension_is_refused() {
    let dir = Workdir::new("wrong-extension");
    let outcome =
        run::<SpectraFilterWindowMower>(&["-in", &swm_input(), "-out", &dir.file("out.dta")]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome
            .err
            .contains("Invalid output file extension for file"),
        "{}",
        outcome.err
    );
    assert!(
        outcome.err.contains("Valid file extensions are: 'mzML'."),
        "{}",
        outcome.err
    );
}

/// Oracle `missing_required_output`: exit 7.
#[test]
fn a_missing_required_output_is_missing_parameters() {
    let outcome = run::<SpectraFilterWindowMower>(&["-in", &swm_input()]);
    assert_eq!(outcome.code, ExitCode::MissingParameters, "{}", outcome.err);
    assert!(
        outcome.err.contains(
            "Error: The required parameter 'out' [valid: mzML] was not given or is empty!"
        ),
        "{}",
        outcome.err
    );
}

/// Oracle `corrupt_input`: a parse failure inside the tool is exit 3
/// (TOPPBase.cpp:460-465), not PARSE_ERROR.
#[test]
fn a_corrupt_input_is_input_file_corrupt() {
    let dir = Workdir::new("corrupt");
    let outcome = run::<SpectraFilterWindowMower>(&[
        "-in",
        &lifecycle("corrupt.mzML"),
        "-out",
        &dir.file("out.mzML"),
    ]);
    assert_eq!(outcome.code, ExitCode::InputFileCorrupt, "{}", outcome.err);
    assert!(
        outcome.err.contains("Error: Unable to read file ("),
        "{}",
        outcome.err
    );
}

// ---------------------------------------------------------------------------
// Oracle cases: subsection parameters on the command line
// ---------------------------------------------------------------------------

/// Oracle `algorithm_peakcount_1`: the override reaches the algorithm and the
/// output equals the C++ output.
#[test]
fn a_subsection_value_on_the_command_line_is_applied() {
    let dir = Workdir::new("peakcount-1");
    let out = dir.file("out.mzML");
    let outcome = run::<SpectraFilterWindowMower>(&[
        "-test",
        "-in",
        &swm_input(),
        "-out",
        &out,
        "-algorithm:peakcount",
        "1",
    ]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_same_peaks(&out, lifecycle("swm_algorithm_peakcount_1.mzML"));
    assert!(total_peaks(&out) < total_peaks(fixture("window_mower_tool_output.mzML")));
}

/// Oracle `algorithm_peakcount_twice`: the last occurrence wins, with a warning.
#[test]
fn a_repeated_option_keeps_the_last_value() {
    let dir = Workdir::new("peakcount-twice");
    let out = dir.file("out.mzML");
    let outcome = run::<SpectraFilterWindowMower>(&[
        "-test",
        "-in",
        &swm_input(),
        "-out",
        &out,
        "-algorithm:peakcount",
        "5",
        "-algorithm:peakcount",
        "1",
    ]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert!(
        outcome.err.contains(
            "Warning: Duplicate parameter '-algorithm:peakcount' given. Using last occurrence with value '1' (ignoring '5')."
        ),
        "{}",
        outcome.err
    );
    assert_same_peaks(&out, lifecycle("swm_algorithm_peakcount_1.mzML"));
}

/// Oracle `algorithm_bogus`: an unknown subsection name is an unknown option.
#[test]
fn an_unknown_subsection_parameter_is_refused() {
    let dir = Workdir::new("algorithm-bogus");
    let outcome = run::<SpectraFilterWindowMower>(&[
        "-test",
        "-in",
        &swm_input(),
        "-out",
        &dir.file("out.mzML"),
        "-algorithm:bogus",
        "1",
    ]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome
            .err
            .contains("Unknown option(s) '[-algorithm:bogus]' given. Aborting!"),
        "{}",
        outcome.err
    );
}

/// Oracle `algorithm_movetype_sideways`: the strict update rejects a value
/// outside the valid strings.
#[test]
fn an_invalid_subsection_value_on_the_command_line_is_refused() {
    let dir = Workdir::new("movetype-sideways");
    let outcome = run::<SpectraFilterWindowMower>(&[
        "-test",
        "-in",
        &swm_input(),
        "-out",
        &dir.file("out.mzML"),
        "-algorithm:movetype",
        "sideways",
    ]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome.err.contains(
            "Invalid string parameter value 'sideways' for parameter 'movetype' given! Valid values are: 'slide,jump'. Updating failed!"
        ),
        "{}",
        outcome.err
    );
    assert!(
        outcome.err.contains(
            "Parameters passed to 'SpectraFilterWindowMower' are invalid. To prevent usage of wrong defaults, please update/fix the parameters!"
        ),
        "{}",
        outcome.err
    );
}

/// Oracle `algorithm_peakcount_not_int`: subsection values convert by type.
#[test]
fn a_non_numeric_subsection_integer_is_refused() {
    let dir = Workdir::new("peakcount-abc");
    let outcome = run::<SpectraFilterWindowMower>(&[
        "-test",
        "-in",
        &swm_input(),
        "-out",
        &dir.file("out.mzML"),
        "-algorithm:peakcount",
        "abc",
    ]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome.err.contains(
            "Invalid parameter values (ConversionError): Could not convert string 'abc' to an integer value. Aborting!"
        ),
        "{}",
        outcome.err
    );
}

// ---------------------------------------------------------------------------
// Oracle cases: INI handling (TOPPBase.cpp:274-367)
// ---------------------------------------------------------------------------

fn swm_with_ini(dir: &Workdir, ini: &str) -> (Outcome, String) {
    let out = dir.file("out.mzML");
    let outcome = run::<SpectraFilterWindowMower>(&[
        "-test",
        "-ini",
        &lifecycle(ini),
        "-in",
        &swm_input(),
        "-out",
        &out,
    ]);
    (outcome, out)
}

/// Oracle `ini_unknown_item`: exit 6.
#[test]
fn an_unknown_ini_item_is_refused() {
    let dir = Workdir::new("ini-unknown");
    let (outcome, _) = swm_with_ini(&dir, "unknown_item.ini");
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome.err.contains(
            "Unknown (or deprecated) Parameter 'bogus_item' given in outdated parameter file!"
        ),
        "{}",
        outcome.err
    );
    assert!(
        outcome
            .err
            .contains("Parameters passed to 'SpectraFilterWindowMower' are invalid."),
        "{}",
        outcome.err
    );
}

/// Oracle `ini_foreign_section`: exit 0 with a warning; defaults apply.
#[test]
fn an_ini_for_another_tool_warns_and_applies_defaults() {
    let dir = Workdir::new("ini-foreign");
    let (outcome, out) = swm_with_ini(&dir, "foreign_section.ini");
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert!(
        outcome.err.contains(
            "Warning: The provided INI file does not contain any parameters specific for this tool (expected in 'SpectraFilterWindowMower:1:')."
        ),
        "{}",
        outcome.err
    );
    assert_same_peaks(&out, fixture("window_mower_tool_output.mzML"));
}

/// Oracle `ini_invalid_subsection_value`: exit 6 (TOPPBase.cpp:339-342).
#[test]
fn an_invalid_subsection_value_in_an_ini_is_refused() {
    let dir = Workdir::new("ini-sideways");
    let (outcome, out) = swm_with_ini(&dir, "movetype_sideways.ini");
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome.err.contains(
            "Invalid string parameter value 'sideways' for parameter 'movetype' given! Valid values are: 'slide,jump'. Updating failed!"
        ),
        "{}",
        outcome.err
    );
    assert!(
        outcome
            .err
            .contains("Parameters passed to 'SpectraFilterWindowMower' are invalid."),
        "{}",
        outcome.err
    );
    assert!(!Path::new(&out).exists());
}

/// Oracle `ini_value_type_mismatch`: exit 6.
#[test]
fn an_ini_value_of_the_wrong_type_is_refused() {
    let dir = Workdir::new("ini-type");
    let (outcome, _) = swm_with_ini(&dir, "peakcount_string.ini");
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome
            .err
            .contains("Parameter 'algorithm:peakcount' has changed value type!\n Updating failed!"),
        "{}",
        outcome.err
    );
}

/// Oracle `ini_version_mismatch`: exit 0 with a notice on stdout.
#[test]
fn an_ini_from_another_version_is_noted() {
    let dir = Workdir::new("ini-version");
    let (outcome, out) = swm_with_ini(&dir, "version_mismatch.ini");
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert!(
        outcome.out.contains(
            "Warning: Parameters file version (0.0.1) does not match the version of this tool (1.0.0)."
        ),
        "{}",
        outcome.out
    );
    assert_same_peaks(&out, fixture("window_mower_tool_output.mzML"));
}

/// Oracle `ini_malformed`: exit 3. The INI is loaded inside the source's
/// run-phase `try` (TOPPBase.cpp:258, 296), so its `ParseError` is
/// INPUT_FILE_CORRUPT; it is not an initialisation error.
#[test]
fn a_malformed_ini_is_input_file_corrupt() {
    let dir = Workdir::new("ini-malformed");
    let (outcome, _) = swm_with_ini(&dir, "malformed.ini");
    assert_eq!(outcome.code, ExitCode::InputFileCorrupt, "{}", outcome.err);
    assert!(
        outcome.err.contains("Error: Unable to read file ("),
        "{}",
        outcome.err
    );
}

/// Oracle `ini_missing`: exit 1.
#[test]
fn a_missing_ini_is_input_file_not_found() {
    let dir = Workdir::new("ini-missing");
    let outcome = run::<SpectraFilterWindowMower>(&[
        "-test",
        "-ini",
        &dir.file("absent.ini"),
        "-in",
        &swm_input(),
        "-out",
        &dir.file("out.mzML"),
    ]);
    assert_eq!(outcome.code, ExitCode::InputFileNotFound, "{}", outcome.err);
    assert!(
        outcome.err.contains("Error: File not found"),
        "{}",
        outcome.err
    );
}

/// Oracle `ini_common_tool_section`: a `common:<tool>:` value is found by leaf
/// name and applied.
#[test]
fn a_common_tool_section_is_applied() {
    let dir = Workdir::new("ini-common-tool");
    let (outcome, out) = swm_with_ini(&dir, "common_tool_section.ini");
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert!(
        outcome
            .err
            .contains("does not contain any parameters specific for this tool"),
        "{}",
        outcome.err
    );
    assert!(
        outcome.err.contains(
            "Found 'SpectraFilterWindowMower:algorithm:peakcount' as 'algorithm:peakcount' in new param."
        ),
        "{}",
        outcome.err
    );
    assert_same_peaks(&out, lifecycle("swm_algorithm_peakcount_1.mzML"));
}

/// Oracle `ini_instance_and_common`: the source applies the `common:` copy of a
/// subsection value after the instance value, so `peakcount` 1 from the common
/// section wins over 2 from the instance section. Reproduced, not endorsed; see
/// the C++ issue candidate in `docs/TOPP_CLI_SUPPORT.md`.
#[test]
fn a_common_subsection_value_overrides_the_instance_value_as_in_the_source() {
    let dir = Workdir::new("ini-instance-common");
    let (outcome, out) = swm_with_ini(&dir, "instance_and_common.ini");
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert!(
        outcome.err.contains(
            "Found 'SpectraFilterWindowMower:algorithm:peakcount' as 'algorithm:peakcount' in new param."
        ),
        "{}",
        outcome.err
    );
    assert_same_peaks(&out, lifecycle("swm_algorithm_peakcount_1.mzML"));
}

/// Oracle `ini_common_top_level`: a `common:<tool>:` value of a top-level
/// parameter has no nested match, so the strict update rejects it.
#[test]
fn a_common_top_level_value_is_refused_as_in_the_source() {
    let dir = Workdir::new("ini-common-top");
    let (outcome, _) = swm_with_ini(&dir, "common_top_level.ini");
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome.err.contains(
            "Unknown (or deprecated) Parameter 'SpectraFilterWindowMower:threads' given in outdated parameter file!"
        ),
        "{}",
        outcome.err
    );
}

// ---------------------------------------------------------------------------
// Oracle cases: write commands (TOPPBase.cpp:2546-2686)
// ---------------------------------------------------------------------------

/// Oracle `write_ini_with_cli_value`: the INI holds defaults, never
/// command-line values.
#[test]
fn write_ini_ignores_command_line_values() {
    let dir = Workdir::new("write-ini-level");
    let ini = dir.file("written.ini");
    let outcome = run::<DTAExtractor>(&["-write_ini", &ini, "-level", "2"]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let written = openms::format::paramxml::load(&ini).unwrap();
    assert_eq!(
        written
            .value("DTAExtractor:1:level")
            .unwrap()
            .as_str()
            .unwrap(),
        "1,2,3"
    );
    assert_eq!(
        written
            .value("DTAExtractor:version")
            .unwrap()
            .as_str()
            .unwrap(),
        "1.0.0"
    );
    assert!(!written.exists("DTAExtractor:1:instance").unwrap());
}

/// Oracle `write_ini_with_ini`: `-ini` updates the written defaults leniently.
#[test]
fn write_ini_with_an_invalid_ini_value_keeps_the_default() {
    let dir = Workdir::new("write-ini-ini");
    let ini = dir.file("written.ini");
    let outcome = run::<SpectraFilterWindowMower>(&[
        "-write_ini",
        &ini,
        "-ini",
        &lifecycle("movetype_sideways.ini"),
    ]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert!(
        outcome.err.contains(
            "Invalid string parameter value 'sideways' for parameter 'movetype' given! Valid values are: 'slide,jump'. Ignoring invalid value (using new default 'slide')!"
        ),
        "{}",
        outcome.err
    );
    let written = openms::format::paramxml::load(&ini).unwrap();
    assert_eq!(
        written
            .value("SpectraFilterWindowMower:1:algorithm:movetype")
            .unwrap()
            .as_str()
            .unwrap(),
        "slide"
    );
}

/// Oracle `write_ini_unwritable`: exit 5.
#[test]
fn write_ini_to_an_unwritable_path_is_refused() {
    let outcome = run::<DTAExtractor>(&["-write_ini", "/nonexistent_dir/x.ini"]);
    assert_eq!(
        outcome.code,
        ExitCode::CannotWriteOutputFile,
        "{}",
        outcome.err
    );
    assert!(
        outcome
            .err
            .contains("Cannot write output file given from parameter '-write_ini'!"),
        "{}",
        outcome.err
    );
}

/// Oracle `write_cwl`, `write_nested_cwl`, `write_json` and `write_nested_json`
/// exit 12 in a build without TDL support; `write_ctd` exits 0 there. The
/// writers are not ported, so all five are refused explicitly with the same
/// INTERNAL_ERROR rather than silently doing nothing.
#[test]
fn tool_description_writers_are_refused_explicitly() {
    let dir = Workdir::new("description-writers");
    for writer in [
        "-write_ctd",
        "-write_cwl",
        "-write_nested_cwl",
        "-write_json",
        "-write_nested_json",
    ] {
        let outcome = run::<DTAExtractor>(&[writer, &text(&dir.0)]);
        assert_eq!(
            outcome.code,
            ExitCode::InternalError,
            "{writer}: {}",
            outcome.err
        );
        assert!(
            outcome
                .err
                .contains(&format!("'{writer}' is not supported")),
            "{writer}: {}",
            outcome.err
        );
    }
}

// ---------------------------------------------------------------------------
// Native cases: phases, streams and context services (W1.4 a, b, f)
// ---------------------------------------------------------------------------

/// A tool whose body writes to both streams and fails in a chosen way.
struct StreamTool;
thread_local! {
    static STREAM_RESULT: RefCell<Option<Result<ExitCode>>> = const { RefCell::new(None) };
}
impl Tool for StreamTool {
    const NAME: &'static str = "StreamTool";
    const DESCRIPTION: &'static str = "Writes to its streams";
    fn register(_spec: &mut ToolSpec) -> Result<()> {
        Ok(())
    }
    fn run(_ctx: &ToolContext) -> Result<ExitCode> {
        unreachable!("run_io is overridden")
    }
    fn run_io(_ctx: &ToolContext, out: &mut dyn Write, err: &mut dyn Write) -> Result<ExitCode> {
        writeln!(out, "report line")?;
        writeln!(err, "diagnostic line")?;
        STREAM_RESULT
            .with(|slot| slot.borrow_mut().take())
            .unwrap_or(Ok(ExitCode::ExecutionOk))
    }
}
fn stream_tool(result: Result<ExitCode>) -> Outcome {
    STREAM_RESULT.with(|slot| *slot.borrow_mut() = Some(result));
    run::<StreamTool>(&["-test"])
}

#[test]
fn run_io_output_reaches_the_callers_streams() {
    let outcome = stream_tool(Ok(ExitCode::InputFileEmpty));
    assert_eq!(outcome.code, ExitCode::InputFileEmpty);
    assert_eq!(outcome.out, "report line\n");
    assert_eq!(outcome.err, "diagnostic line\n");
}

#[test]
fn run_phase_errors_map_like_the_source_inner_catch() {
    let parse = stream_tool(Err(Error::Parse {
        line: 3,
        message: "unexpected end".into(),
    }));
    assert_eq!(parse.code, ExitCode::InputFileCorrupt, "{}", parse.err);
    assert!(
        parse.err.contains("Error: Unable to read file ("),
        "{}",
        parse.err
    );

    let missing = stream_tool(Err(Error::Io(std::io::Error::from(
        std::io::ErrorKind::NotFound,
    ))));
    assert_eq!(missing.code, ExitCode::InputFileNotFound, "{}", missing.err);

    let invalid = stream_tool(Err(Error::InvalidValue("bad".into())));
    assert_eq!(invalid.code, ExitCode::IllegalParameters, "{}", invalid.err);

    let unsupported = stream_tool(Err(Error::Unsupported("not here".into())));
    assert_eq!(
        unsupported.code,
        ExitCode::IncompatibleInputData,
        "{}",
        unsupported.err
    );
}

struct FailingRegistration;
impl Tool for FailingRegistration {
    const NAME: &'static str = "FailingRegistration";
    const DESCRIPTION: &'static str = "Registers a parameter twice";
    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_flag("twice", "first", false)?;
        spec.register_flag("twice", "second", false)
    }
    fn run(_ctx: &ToolContext) -> Result<ExitCode> {
        Ok(ExitCode::ExecutionOk)
    }
}

/// An initialisation failure is ILLEGAL_PARAMETERS (TOPPBase.cpp:505-508).
#[test]
fn an_initialisation_failure_is_illegal_parameters() {
    let outcome = run::<FailingRegistration>(&["-test"]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome
            .err
            .contains("Unable to initialize or run FailingRegistration:"),
        "{}",
        outcome.err
    );
}

#[test]
fn an_oversized_command_line_is_refused_before_parsing() {
    let arguments = vec!["-test".to_owned(); MAX_ARGUMENTS + 1];
    let outcome = run_arguments::<ToppBaseTest>(&arguments);
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
}

/// `-test` seeds the unique-id generator (TOPPBase.cpp:369-376). The two draws
/// are those the cli tracer recorded for MT19937-64 seeded 19991231235959.
#[test]
fn test_mode_seeds_the_unique_id_generator() {
    assert_eq!(run::<ToppBaseTest>(&["-test"]).code, ExitCode::ExecutionOk);
    let ctx = captured();
    assert!(ctx.test_mode());
    let mut generator = ctx.unique_id_generator();
    assert_eq!(generator.seed(), TEST_MODE_UNIQUE_ID_SEED);
    assert_eq!(generator.get_unique_id(), 5_233_264_595_117_471_314);
    assert_eq!(generator.get_unique_id(), 4_835_329_514_588_776_807);

    assert_eq!(run::<ToppBaseTest>(&["-flag"]).code, ExitCode::ExecutionOk);
    let ctx = captured();
    assert!(!ctx.test_mode());
    assert_ne!(ctx.unique_id_generator().seed(), TEST_MODE_UNIQUE_ID_SEED);
}

/// `-no_progress`, `-threads`, `-debug` and `-force` reach the context
/// (TOPPBase.cpp:393-408).
#[test]
fn common_options_become_context_services() {
    assert_eq!(run::<ToppBaseTest>(&["-flag"]).code, ExitCode::ExecutionOk);
    let ctx = captured();
    assert_eq!(ctx.progress_log_type(), ProgressLogType::Cmd);
    assert_eq!(ctx.thread_policy(), Threads::serial());
    assert_eq!(ctx.debug_level(), 0);
    assert!(!ctx.force());
    assert_eq!(ctx.ini_location(), "TOPPBaseTest:1:");
    assert_eq!(ctx.tool_name(), "TOPPBaseTest");
    assert_eq!(ctx.version(), "1.0.0");

    let outcome = run::<ToppBaseTest>(&["-no_progress", "-threads", "0", "-debug", "3", "-force"]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let ctx = captured();
    assert_eq!(ctx.progress_log_type(), ProgressLogType::None);
    assert_eq!(ctx.thread_policy(), Threads::all());
    assert_eq!(ctx.threads(), 0);
    assert_eq!(ctx.debug_level(), 3);
    assert!(ctx.force());
    assert!(ctx.no_progress());

    assert_eq!(
        run::<ToppBaseTest>(&["-threads", "4"]).code,
        ExitCode::ExecutionOk
    );
    assert_eq!(captured().thread_policy(), Threads::from_cli(4));
}

/// Test-mode processing record (TOPPBase.cpp:543-555).
#[test]
fn test_mode_processing_info_is_machine_independent() {
    assert_eq!(run::<ToppBaseTest>(&["-test"]).code, ExitCode::ExecutionOk);
    let ctx = captured();
    let processing = ctx
        .processing_info(&[ProcessingAction::Quantitation])
        .unwrap();
    assert_eq!(processing.software.name, "TOPPBaseTest");
    assert_eq!(processing.software.version, TEST_MODE_VERSION);
    assert_eq!(
        processing.completion_time,
        Some(DateTime::parse(TEST_MODE_COMPLETION_TIME).unwrap())
    );
    assert_eq!(processing.metadata.len(), 1);
    assert_eq!(
        processing.metadata.get(TEST_MODE_PARAMETER_KEY),
        Some(&MetaValue::from(TEST_MODE_PARAMETER_VALUE))
    );
    assert_eq!(
        processing.actions.iter().copied().collect::<Vec<_>>(),
        vec![ProcessingAction::Quantitation]
    );
}

/// Processing record outside test mode (TOPPBase.cpp:556-568).
#[test]
fn processing_info_records_the_version_time_and_parameters() {
    assert_eq!(
        run::<ToppBaseTest>(&["-stringoption", "commandline"]).code,
        ExitCode::ExecutionOk
    );
    let ctx = captured();
    let processing = ctx
        .processing_info(&[ProcessingAction::PeakPicking])
        .unwrap();
    assert_eq!(processing.software.version, "1.0.0");
    let completion = processing.completion_time.expect("completion time is set");
    assert_ne!(completion, DateTime::default());
    assert_eq!(
        processing.metadata.get("parameter: stringoption"),
        Some(&MetaValue::from("commandline"))
    );
    assert_eq!(
        processing.metadata.get("parameter: intoption"),
        Some(&MetaValue::from(4711_i64))
    );
    assert!(!processing.metadata.contains_key(TEST_MODE_PARAMETER_KEY));
}

/// A consensus map keeps base names only under `-test` (TOPPBase.cpp:573-585).
#[test]
fn test_mode_strips_consensus_column_paths() {
    assert_eq!(run::<ToppBaseTest>(&["-test"]).code, ExitCode::ExecutionOk);
    let ctx = captured();
    let processing = ctx
        .processing_info(&[ProcessingAction::FeatureGrouping])
        .unwrap();
    let mut map = ConsensusMap::new();
    map.column_headers.insert(
        0,
        ColumnHeader {
            filename: "/data/run/input_1.featureXML".into(),
            ..ColumnHeader::default()
        },
    );
    ctx.add_data_processing(&mut map, &processing);
    assert_eq!(map.data_processing.len(), 1);
    assert_eq!(map.column_headers[&0].filename, "input_1.featureXML");
}

// ---------------------------------------------------------------------------
// Upstream TOPPBase_test.cpp sections
// ---------------------------------------------------------------------------

/// `[EXTRA] getIniLocation_` (TOPPBase_test.cpp:432-440), default part. The
/// `-instance 5` part cannot apply: the source rejects `-instance` in the
/// strict update (oracle `instance_on_command_line`).
#[test]
fn upstream_ini_location() {
    assert_eq!(
        run_arguments::<ToppBaseTest>(&[]).code,
        ExitCode::ExecutionOk
    );
    assert_eq!(captured().ini_location(), "TOPPBaseTest:1:");
}

/// `[EXTRA] getStringOption_` (TOPPBase_test.cpp:443-481). The INI cases read
/// values the source leaves in `param_` after its update has failed (a
/// `common:` value of a top-level parameter is rejected, oracle
/// `ini_common_top_level`), and the write_ini comparison is W2.2; neither is
/// transcribed here.
#[test]
fn upstream_get_string_option() {
    assert_eq!(
        run_arguments::<ToppBaseTest>(&[]).code,
        ExitCode::ExecutionOk
    );
    assert_eq!(captured().string("stringoption").unwrap(), "string default");

    assert_eq!(
        run::<ToppBaseTest>(&["-stringoption", "commandline"]).code,
        ExitCode::ExecutionOk
    );
    let ctx = captured();
    assert_eq!(ctx.string("stringoption").unwrap(), "commandline");
    assert!(ctx.string("doubleoption").is_err());
    assert!(ctx.string("imleeewenit").is_err());

    // Missing required parameter: RequiredParameterNotGiven is MISSING_PARAMETERS.
    let outcome = run::<ToppBaseTestNop>(&["-flag"]);
    assert_eq!(outcome.code, ExitCode::MissingParameters, "{}", outcome.err);
    assert!(outcome.err.contains("'stringoption'"), "{}", outcome.err);
}

/// `[EXTRA] getIntOption_` (TOPPBase_test.cpp:539-553).
#[test]
fn upstream_get_int_option() {
    assert_eq!(
        run_arguments::<ToppBaseTest>(&[]).code,
        ExitCode::ExecutionOk
    );
    assert_eq!(captured().int("intoption").unwrap(), 4711);
    assert_eq!(
        run::<ToppBaseTest>(&["-intoption", "6"]).code,
        ExitCode::ExecutionOk
    );
    let ctx = captured();
    assert_eq!(ctx.int("intoption").unwrap(), 6);
    assert!(ctx.int("doubleoption").is_err());
    assert!(ctx.int("imleeewenit").is_err());
}

/// `[EXTRA] getDoubleOption_` (TOPPBase_test.cpp:555-566).
#[test]
fn upstream_get_double_option() {
    assert_eq!(
        run_arguments::<ToppBaseTest>(&[]).code,
        ExitCode::ExecutionOk
    );
    assert!((captured().double("doubleoption").unwrap() - 0.4711).abs() < 1e-5);
    assert_eq!(
        run::<ToppBaseTest>(&["-doubleoption", "4.5"]).code,
        ExitCode::ExecutionOk
    );
    let ctx = captured();
    assert!((ctx.double("doubleoption").unwrap() - 4.5).abs() < 1e-5);
    assert!(ctx.double("intoption").is_err());
    assert!(ctx.double("imleeewenit").is_err());
}

/// `[EXTRA] getIntList_` (TOPPBase_test.cpp:568-588).
#[test]
fn upstream_get_int_list() {
    assert_eq!(
        run_arguments::<ToppBaseTest>(&[]).code,
        ExitCode::ExecutionOk
    );
    assert_eq!(captured().int_list("intlist").unwrap(), &[1, 2, 3, 4]);
    assert_eq!(
        run::<ToppBaseTest>(&["-intlist", "6", "5", "4711"]).code,
        ExitCode::ExecutionOk
    );
    let ctx = captured();
    assert_eq!(ctx.int_list("intlist").unwrap(), &[6, 5, 4711]);
    assert!(ctx.int_list("intoption").is_err());
    assert!(ctx.int_list("imleeewenit").is_err());
    assert_eq!(
        run::<ToppBaseTest>(&["-intlist", "6"]).code,
        ExitCode::ExecutionOk
    );
    assert_eq!(captured().int_list("intlist").unwrap(), &[6]);

    let outcome = run::<ToppBaseTestNop>(&["-flag", "-stringoption", "x", "-stringlist", "a"]);
    assert_eq!(outcome.code, ExitCode::MissingParameters, "{}", outcome.err);
    assert!(outcome.err.contains("'intlist'"), "{}", outcome.err);
}

/// `[EXTRA] getDoubleList_` (TOPPBase_test.cpp:590-613).
#[test]
fn upstream_get_double_list() {
    assert_eq!(
        run_arguments::<ToppBaseTest>(&[]).code,
        ExitCode::ExecutionOk
    );
    assert_eq!(
        captured().double_list("doublelist").unwrap(),
        &[0.4711, 1.022, 4.0]
    );
    assert_eq!(
        run::<ToppBaseTest>(&["-doublelist", "0.411"]).code,
        ExitCode::ExecutionOk
    );
    let ctx = captured();
    assert_eq!(ctx.double_list("doublelist").unwrap(), &[0.411]);
    assert!(ctx.double_list("intoption").is_err());
    assert!(ctx.double_list("imleeewenit").is_err());
    assert_eq!(
        run::<ToppBaseTest>(&["-doublelist", "0.411", "4.5", "4.0"]).code,
        ExitCode::ExecutionOk
    );
    assert_eq!(
        captured().double_list("doublelist").unwrap(),
        &[0.411, 4.5, 4.0]
    );
    assert_eq!(
        run::<ToppBaseTest>(&["-doublelist", "0.411", "4.5"]).code,
        ExitCode::ExecutionOk
    );
    assert_eq!(captured().double_list("doublelist").unwrap(), &[0.411, 4.5]);

    let outcome = run::<ToppBaseTestNop>(&[
        "-flag",
        "-stringoption",
        "x",
        "-stringlist",
        "a",
        "-intlist",
        "1",
    ]);
    assert_eq!(outcome.code, ExitCode::MissingParameters, "{}", outcome.err);
    assert!(outcome.err.contains("'doublelist'"), "{}", outcome.err);
}

/// `[EXTRA] getStringList_` (TOPPBase_test.cpp:615-640).
#[test]
fn upstream_get_string_list() {
    assert_eq!(
        run_arguments::<ToppBaseTest>(&[]).code,
        ExitCode::ExecutionOk
    );
    assert_eq!(
        captured().string_list("stringlist").unwrap(),
        &["abc", "def", "ghi", "jkl"]
    );
    assert_eq!(
        run::<ToppBaseTest>(&["-stringlist", "commandline"]).code,
        ExitCode::ExecutionOk
    );
    let ctx = captured();
    assert_eq!(ctx.string_list("stringlist").unwrap(), &["commandline"]);
    assert!(ctx.string_list("intoption").is_err());
    assert!(ctx.string_list("imleeewenit").is_err());

    let toolcommon = lifecycle("TOPPBase_toolcommon.ini");
    let common = lifecycle("TOPPBase_common.ini");
    assert_eq!(
        run::<ToppBaseTest>(&["-stringlist", "commandline", &toolcommon, &common]).code,
        ExitCode::ExecutionOk
    );
    assert_eq!(
        captured().string_list("stringlist").unwrap(),
        &["commandline".to_owned(), toolcommon, common]
    );

    let outcome = run::<ToppBaseTestNop>(&["-flag", "-stringoption", "x"]);
    assert_eq!(outcome.code, ExitCode::MissingParameters, "{}", outcome.err);
    assert!(outcome.err.contains("'stringlist'"), "{}", outcome.err);
}

/// `[EXTRA] getFlag_` (TOPPBase_test.cpp:642-653).
#[test]
fn upstream_get_flag() {
    assert_eq!(
        run_arguments::<ToppBaseTest>(&[]).code,
        ExitCode::ExecutionOk
    );
    assert!(!captured().flag("flag").unwrap());
    assert_eq!(run::<ToppBaseTest>(&["-flag"]).code, ExitCode::ExecutionOk);
    let ctx = captured();
    assert!(ctx.flag("flag").unwrap());
    assert!(ctx.flag("doubleoption").is_err());
    assert!(ctx.flag("imleeewenit").is_err());
}

/// `[EXTRA] inputFileReadable_` (TOPPBase_test.cpp:655-660).
#[test]
fn upstream_input_file_readable() {
    let mut sink = Vec::new();
    assert_eq!(
        input_file_readable("/this/file/does/not/exist.txt", "someparam", &mut sink),
        Some(ExitCode::InputFileNotFound)
    );
    assert_eq!(
        input_file_readable(&lifecycle("TOPPBase_empty.txt"), "someparam", &mut sink),
        Some(ExitCode::InputFileEmpty)
    );
    assert_eq!(
        input_file_readable(&lifecycle("TOPPBase_common.ini"), "ini", &mut sink),
        None
    );
}

/// `[EXTRA] outputFileWritable_` (TOPPBase_test.cpp:662-673).
#[test]
fn upstream_output_file_writable() {
    let mut sink = Vec::new();
    assert_eq!(
        output_file_writable(
            "/this/file/cannot/be/written/does_not_exists.txt",
            "someparam",
            &mut sink
        ),
        Some(ExitCode::CannotWriteOutputFile)
    );
    let dir = Workdir::new("output-writable");
    let file = dir.file("source.tmp");
    assert_eq!(output_file_writable(&file, "", &mut sink), None);
    // Asking must not create the file under the caller's name.
    assert!(!Path::new(&file).exists());
}

/// `[EXTRA] parseRange_` (TOPPBase_test.cpp:675-710).
#[test]
fn upstream_parse_range() {
    let (mut a, mut b) = (-1.0, -1.0);
    assert!(!parse_range(":", &mut a, &mut b).unwrap());
    assert_eq!((a, b), (-1.0, -1.0));

    assert!(parse_range("4.5:", &mut a, &mut b).unwrap());
    assert!((a - 4.5).abs() < 1e-5 && (b + 1.0).abs() < 1e-5);

    assert!(parse_range(":5.5", &mut a, &mut b).unwrap());
    assert!((a - 4.5).abs() < 1e-5 && (b - 5.5).abs() < 1e-5);

    assert!(parse_range("6.5:7.5", &mut a, &mut b).unwrap());
    assert!((a - 6.5).abs() < 1e-5 && (b - 7.5).abs() < 1e-5);

    // A colon-less range is malformed and must fail loudly.
    assert!(parse_range("400", &mut a, &mut b).is_err());
}

/// `[EXTRA] data processing methods` (TOPPBase_test.cpp:715-731), with the
/// feature map and consensus map the source helper also fills.
#[test]
fn upstream_data_processing_methods() {
    assert_eq!(
        run_arguments::<ToppBaseTest>(&[]).code,
        ExitCode::ExecutionOk
    );
    let ctx = captured();
    let processing = ctx
        .processing_info(&[ProcessingAction::RetentionTimeAlignment])
        .unwrap();

    let mut experiment = MSExperiment::new();
    experiment.spectra = vec![MSSpectrum::default(), MSSpectrum::default()];
    experiment.chromatograms = vec![MSChromatogram::default()];
    ctx.add_data_processing(&mut experiment, &processing);
    for spectrum in &experiment.spectra {
        assert_eq!(spectrum.data_processing.len(), 1);
        let entry = &spectrum.data_processing[0];
        assert_eq!(entry.software.name, "TOPPBaseTest");
        assert_ne!(entry.software.version, "1.1");
        assert!(
            entry
                .completion_time
                .is_some_and(|t| t != DateTime::default())
        );
        assert_eq!(entry.actions.len(), 1);
        assert_eq!(
            entry.actions.iter().next(),
            Some(&ProcessingAction::RetentionTimeAlignment)
        );
    }
    // One shared record, as the source pushes one shared pointer everywhere.
    assert!(Arc::ptr_eq(
        &experiment.spectra[0].data_processing[0],
        &experiment.spectra[1].data_processing[0]
    ));
    assert!(Arc::ptr_eq(
        &experiment.spectra[0].data_processing[0],
        &experiment.chromatograms[0].data_processing[0]
    ));

    let mut features = FeatureMap::default();
    ctx.add_data_processing(&mut features, &processing);
    assert_eq!(features.data_processing, vec![processing.clone()]);
    let mut consensus = ConsensusMap::new();
    ctx.add_data_processing(&mut consensus, &processing);
    assert_eq!(consensus.data_processing, vec![processing]);
}

/// `[EXTRA] const Param& getParam_()` (TOPPBase_test.cpp:733-752).
#[test]
fn upstream_get_param() {
    assert_eq!(
        run_arguments::<ToppBaseTestParam>(&[]).code,
        ExitCode::ExecutionOk
    );
    let result = captured();
    let expected = full_param();
    for item in expected.iter().unwrap() {
        let entry = result.param().entry(&item.key).unwrap();
        // Source ParamEntry::operator== compares name and value.
        assert_eq!(entry.name, item.entry.name, "{}", item.key);
        assert_eq!(entry.value, item.entry.value, "{}", item.key);
    }
}

/// `setMaxNumberOfThreads` is NOT_TESTABLE upstream (TOPPBase_test.cpp:754-761);
/// the native policy it becomes is covered by
/// `common_options_become_context_services`.
#[test]
fn upstream_misc_options_on_command_line() {
    // TOPPBase_test.cpp:763-777: trailing text, then an unknown option.
    let outcome = run_arguments::<ToppBaseCmdParseTest>(&[
        "TOPPBaseTest".into(),
        "commandline".into(),
        "-test".into(),
    ]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters);
    let outcome = run_arguments::<ToppBaseCmdParseTest>(&[
        "TOPPBaseTest".into(),
        "-stringoption".into(),
        "commandline".into(),
        "-test".into(),
    ]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters);
}

/// `[EXTRA] test subsection parameters` (TOPPBase_test.cpp:790-824).
#[test]
fn upstream_subsection_parameters() {
    let outcome =
        run::<ToppBaseCmdParseSubsectionsTest>(&["-stringoption", "commandline", "-test"]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let ctx = captured();
    assert_eq!(ctx.string("stringoption").unwrap(), "commandline");
    assert_eq!(string_value(&ctx, "algorithm:param1"), "param1_value");
    assert_eq!(string_value(&ctx, "algorithm:param2"), "param2_value");
    assert_eq!(string_value(&ctx, "other:param3"), "param3_value");
    assert_eq!(string_value(&ctx, "other:param4"), "param4_value");

    let outcome = run::<ToppBaseCmdParseSubsectionsTest>(&[
        "-stringoption",
        "commandline",
        "-algorithm:param1",
        "val1",
        "-algorithm:param2",
        "val2",
        "-other:param3",
        "val3",
        "-other:param4",
        "val4",
        "-test",
    ]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let ctx = captured();
    assert_eq!(ctx.string("stringoption").unwrap(), "commandline");
    assert_eq!(string_value(&ctx, "algorithm:param1"), "val1");
    assert_eq!(string_value(&ctx, "algorithm:param2"), "val2");
    assert_eq!(string_value(&ctx, "other:param3"), "val3");
    assert_eq!(string_value(&ctx, "other:param4"), "val4");

    let outcome = run::<ToppBaseCmdParseSubsectionsTest>(&[
        "-ini",
        &lifecycle("TOPPBaseCmdParseSubsectionsTest.ini"),
        "-algorithm:param1",
        "val1",
        "-other:param4",
        "val4",
        "-stringoption",
        "commandline",
        "-test",
    ]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let ctx = captured();
    assert_eq!(ctx.string("stringoption").unwrap(), "commandline");
    assert_eq!(string_value(&ctx, "algorithm:param1"), "val1");
    assert_eq!(string_value(&ctx, "algorithm:param2"), "param2_ini_value");
    assert_eq!(string_value(&ctx, "other:param3"), "param3_ini_value");
    assert_eq!(string_value(&ctx, "other:param4"), "val4");

    // Native: a subsection entry restricted to true/false with default false is
    // a flag token (TOPPBase.cpp:900-908); one defaulting to true is a string.
    let outcome = run::<ToppBaseCmdParseSubsectionsTest>(&[
        "-stringoption",
        "x",
        "-other:flagparam",
        "-other:nonflagparam",
        "false",
    ]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let ctx = captured();
    assert_eq!(string_value(&ctx, "other:flagparam"), "true");
    assert_eq!(string_value(&ctx, "other:nonflagparam"), "false");
}

/// `[EXTRA] test duplicate parameters` (TOPPBase_test.cpp:826-852).
#[test]
fn upstream_duplicate_parameters() {
    let outcome = run::<ToppBaseTest>(&["-stringoption", "commandline", "-stringoption", "4711"]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_eq!(captured().string("stringoption").unwrap(), "4711");
    assert!(
        outcome.err.contains(
            "Warning: Duplicate parameter '-stringoption' given. Using last occurrence with value '4711' (ignoring 'commandline')."
        ),
        "{}",
        outcome.err
    );

    let outcome = run::<ToppBaseTest>(&["-intoption", "5", "-intoption", "4711"]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_eq!(captured().int("intoption").unwrap(), 4711);

    let outcome = run::<ToppBaseTest>(&["-doubleoption", "0.411", "-doubleoption", "4.5"]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert!((captured().double("doubleoption").unwrap() - 4.5).abs() < 1e-5);
}

/// `[EXTRA] test flag with trailing arguments` (TOPPBase_test.cpp:854-868).
#[test]
fn upstream_flag_with_trailing_arguments() {
    assert_eq!(
        run::<ToppBaseTest>(&["-flag", "commandline", "-test"]).code,
        ExitCode::IllegalParameters
    );
    assert_eq!(
        run::<ToppBaseTest>(&["-flag", "commandline", "4711", "4.5", "-test"]).code,
        ExitCode::IllegalParameters
    );
}

/// Source `StringUtils::toInt32` and `toDouble` messages (StringUtils.cpp:136-276)
/// reach the parse failure unchanged.
#[test]
fn number_conversion_follows_string_utils() {
    let outcome = run::<ToppBaseTest>(&["-intoption", "4.5"]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters);
    assert!(
        outcome.err.contains(
            "Prefix of string '4.5' successfully converted to an int32 value. Additional characters found at position 2"
        ),
        "{}",
        outcome.err
    );
    let outcome = run::<ToppBaseTest>(&["-doubleoption", "1e999"]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters);
    assert!(
        outcome
            .err
            .contains("Could not convert string '1e999' to a double value"),
        "{}",
        outcome.err
    );
    assert_eq!(
        run::<ToppBaseTest>(&["-intoption", " +12 ", "-doubleoption", "+.5e1"]).code,
        ExitCode::ExecutionOk
    );
    let ctx = captured();
    assert_eq!(ctx.int("intoption").unwrap(), 12);
    assert_eq!(ctx.double("doubleoption").unwrap(), 5.0);
}
