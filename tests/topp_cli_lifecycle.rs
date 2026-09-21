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
//! parallel. `tests/data/topp_cli_lifecycle_provenance.json` records the
//! fixtures, their hashes and the oracle manifest.

// The TOPP framework lives behind `paramxml` (every tool supports -ini) and the
// stand-in tools read mzML, so the whole file is inert without both features.
#![cfg(all(feature = "mzml", feature = "paramxml"))]

#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;

use openms::cli::tools::{
    BaselineFilter, DTAExtractor, MapNormalizer, MzMLSplitter, SpectraFilterWindowMower,
};
use openms::cli::{
    ExitCode, MAX_ARGUMENTS, TEST_MODE_COMPLETION_TIME, TEST_MODE_PARAMETER_KEY,
    TEST_MODE_PARAMETER_VALUE, TEST_MODE_UNIQUE_ID_SEED, TEST_MODE_VERSION, Tool, ToolContext,
    ToolSpec, input_file_readable, output_file_writable, parse_range, run_with, verbose_version,
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
use openms::system::file::TempDir;
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

/// A fresh, uniquely named directory per case, removed when the case ends.
///
/// The case name is kept only to make a failing assertion easier to place; the
/// directory itself comes from the crate's `TempDir`, so concurrent runs never
/// share it.
struct Workdir(TempDir);
impl Workdir {
    fn new(_case: &str) -> Self {
        Self(TempDir::new_in(std::env::temp_dir(), false).unwrap())
    }
    fn path(&self) -> &Path {
        self.0.path()
    }
    fn file(&self, name: &str) -> String {
        text(self.path().join(name))
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

/// Oracle `ini_unreadable` and `write_ini_ini_unreadable`
/// (`../oracle/topp-cli-lifecycle/ini_read_failures/manifest.json`): an INI
/// that exists but cannot be read is exit 2 with the source's `FileNotReadable`
/// wording, before a run and with `-write_ini`, and nothing is written. As in
/// the oracle's `ini_readable_control`, the same INI runs cleanly while it is
/// readable. Skipped when permissions do not restrict the process, as for root.
#[cfg(unix)]
#[test]
fn an_unreadable_ini_is_input_file_not_readable() {
    use std::os::unix::fs::PermissionsExt;
    let dir = Workdir::new("ini-unreadable");
    let ini = dir.file("unreadable.ini");
    let written = run::<SpectraFilterWindowMower>(&["-write_ini", &ini]);
    assert_eq!(written.code, ExitCode::ExecutionOk, "{}", written.err);
    let control = run::<SpectraFilterWindowMower>(&[
        "-test",
        "-ini",
        &ini,
        "-in",
        &swm_input(),
        "-out",
        &dir.file("control.mzML"),
    ]);
    assert_eq!(control.code, ExitCode::ExecutionOk, "{}", control.err);

    fs::set_permissions(&ini, fs::Permissions::from_mode(0o000)).unwrap();
    if fs::File::open(&ini).is_ok() {
        fs::set_permissions(&ini, fs::Permissions::from_mode(0o644)).unwrap();
        return;
    }
    let out = dir.file("out.mzML");
    let before_run = run::<SpectraFilterWindowMower>(&[
        "-test",
        "-ini",
        &ini,
        "-in",
        &swm_input(),
        "-out",
        &out,
    ]);
    let written_ini = dir.file("written.ini");
    let with_write_ini =
        run::<SpectraFilterWindowMower>(&["-write_ini", &written_ini, "-ini", &ini]);
    fs::set_permissions(&ini, fs::Permissions::from_mode(0o644)).unwrap();

    let expected =
        format!("Error: File not readable (the file '{ini}' is not readable for the current user)");
    for outcome in [&before_run, &with_write_ini] {
        assert_eq!(
            outcome.code,
            ExitCode::InputFileNotReadable,
            "{}",
            outcome.err
        );
        assert!(outcome.err.contains(&expected), "{}", outcome.err);
    }
    assert!(!Path::new(&out).exists());
    assert!(!Path::new(&written_ini).exists());
}

/// Oracle `ini_directory` and `write_ini_ini_directory`
/// (`../oracle/topp-cli-lifecycle/ini_read_failures/manifest.json`): a
/// directory given as `-ini` is the source's `ParseError`, exit 3, before a run
/// and with `-write_ini`, and nothing is written.
#[test]
fn an_ini_directory_is_input_file_corrupt() {
    let dir = Workdir::new("ini-directory");
    let ini = dir.file("directory.ini");
    fs::create_dir(&ini).unwrap();
    let out = dir.file("out.mzML");
    let before_run = run::<SpectraFilterWindowMower>(&[
        "-test",
        "-ini",
        &ini,
        "-in",
        &swm_input(),
        "-out",
        &out,
    ]);
    let written_ini = dir.file("written.ini");
    let with_write_ini =
        run::<SpectraFilterWindowMower>(&["-write_ini", &written_ini, "-ini", &ini]);

    let expected =
        format!("Error: Unable to read file (While loading '{ini}': unable to read data from file");
    for outcome in [&before_run, &with_write_ini] {
        assert_eq!(outcome.code, ExitCode::InputFileCorrupt, "{}", outcome.err);
        assert!(outcome.err.contains(&expected), "{}", outcome.err);
    }
    assert!(!Path::new(&out).exists());
    assert!(!Path::new(&written_ini).exists());
}

/// Oracle `ini_dev_null` and `write_ini_ini_dev_null`
/// (`../oracle/topp-cli-lifecycle/ini_read_failures/manifest.json`): the
/// character device `/dev/null` given as `-ini` is neither a regular file nor a
/// directory, but it is readable, and the source reads it as an empty document.
/// That is the source's `ParseError`, exit 3 with an `Error: Unable to read
/// file (` line, before a run and with `-write_ini`; it is not reported as
/// unreadable, and nothing is written.
#[cfg(unix)]
#[test]
fn a_character_device_ini_is_input_file_corrupt() {
    let dir = Workdir::new("ini-dev-null");
    let out = dir.file("out.mzML");
    let before_run = run::<SpectraFilterWindowMower>(&[
        "-test",
        "-ini",
        "/dev/null",
        "-in",
        &swm_input(),
        "-out",
        &out,
    ]);
    let written_ini = dir.file("written.ini");
    let with_write_ini =
        run::<SpectraFilterWindowMower>(&["-write_ini", &written_ini, "-ini", "/dev/null"]);

    for outcome in [&before_run, &with_write_ini] {
        assert_eq!(outcome.code, ExitCode::InputFileCorrupt, "{}", outcome.err);
        assert!(
            outcome
                .err
                .lines()
                .any(|line| line.starts_with("Error: Unable to read file (")),
            "{}",
            outcome.err
        );
        assert!(!outcome.err.contains("not readable"), "{}", outcome.err);
    }
    assert!(!Path::new(&out).exists());
    assert!(!Path::new(&written_ini).exists());
}

/// Oracle `ini_fifo_denied` and `write_ini_ini_fifo_denied`
/// (`../oracle/topp-cli-lifecycle/ini_read_failures/manifest.json`): a FIFO
/// this user cannot open (mode 000) is neither a regular file nor a directory,
/// so it is not asked for readability before the load; the load's open is
/// refused, exit 2 with the source's `FileNotReadable` wording, before a run
/// and with `-write_ini`, and nothing is written. The open is refused before it
/// would wait for a writer, so neither tool blocks. Skipped when `mkfifo` is
/// unavailable or permissions do not restrict the process, as for root, where
/// the open would wait instead.
#[cfg(unix)]
#[test]
fn a_fifo_ini_this_user_cannot_open_is_input_file_not_readable() {
    use std::os::unix::fs::PermissionsExt;
    let dir = Workdir::new("ini-fifo-denied");
    let probe = dir.file("probe.ini");
    fs::write(&probe, b"").unwrap();
    fs::set_permissions(&probe, fs::Permissions::from_mode(0o000)).unwrap();
    if fs::File::open(&probe).is_ok() {
        return;
    }
    let fifo = dir.file("denied.fifo");
    let made = std::process::Command::new("mkfifo").arg(&fifo).status();
    if !made.is_ok_and(|status| status.success()) {
        return;
    }
    fs::set_permissions(&fifo, fs::Permissions::from_mode(0o000)).unwrap();

    let out = dir.file("out.mzML");
    let before_run = run::<SpectraFilterWindowMower>(&[
        "-test",
        "-ini",
        &fifo,
        "-in",
        &swm_input(),
        "-out",
        &out,
    ]);
    let written_ini = dir.file("written.ini");
    let with_write_ini =
        run::<SpectraFilterWindowMower>(&["-write_ini", &written_ini, "-ini", &fifo]);

    let expected = format!(
        "Error: File not readable (the file '{fifo}' is not readable for the current user)"
    );
    for outcome in [&before_run, &with_write_ini] {
        assert_eq!(
            outcome.code,
            ExitCode::InputFileNotReadable,
            "{}",
            outcome.err
        );
        assert!(outcome.err.contains(&expected), "{}", outcome.err);
    }
    assert!(!Path::new(&out).exists());
    assert!(!Path::new(&written_ini).exists());
}

/// A FIFO given as `-ini` is read like a file once a writer opens it: an INI
/// written by the tool and sent through the FIFO by a writer that opens it
/// once is loaded by `-write_ini`, exit 0, and is not refused as unreadable.
///
/// A deliberate difference (tier 4, native): the C++ tool opens the INI twice,
/// once to look for compression (`XMLFile.cpp:141-147`) and once for xerces
/// (`:166`), so after the first open closes no writer is left and the second
/// open waits. The oracle observation `write_ini_ini_fifo_single_writer`
/// (`../oracle/topp-cli-lifecycle/ini_read_failures/manifest.json`) records
/// that: the same single writer sends the whole INI, the C++ tool blocks until
/// its 10-second alarm ends it (exit 142) and writes nothing. This port opens
/// the FIFO once. The tool runs on its own thread with a timeout, so a
/// regression that blocks fails instead of hanging the suite. Skipped when
/// `mkfifo` is unavailable.
#[cfg(unix)]
#[test]
fn an_ini_fifo_is_read_once_a_writer_opens_it() {
    let dir = Workdir::new("ini-fifo");
    let source = dir.file("source.ini");
    let written = run::<SpectraFilterWindowMower>(&["-write_ini", &source]);
    assert_eq!(written.code, ExitCode::ExecutionOk, "{}", written.err);
    let fifo = dir.file("fifo.ini");
    let made = std::process::Command::new("mkfifo").arg(&fifo).status();
    if !made.is_ok_and(|status| status.success()) {
        return;
    }
    let bytes = fs::read(&source).unwrap();

    // Opening the FIFO for writing waits until the tool opens it for reading.
    let writer_path = fifo.clone();
    let writer = std::thread::spawn(move || -> std::io::Result<()> {
        fs::OpenOptions::new()
            .write(true)
            .open(&writer_path)?
            .write_all(&bytes)
    });
    let target = dir.file("written.ini");
    let (tool_target, tool_fifo) = (target.clone(), fifo.clone());
    let (sender, receiver) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let outcome =
            run::<SpectraFilterWindowMower>(&["-write_ini", &tool_target, "-ini", &tool_fifo]);
        let _ = sender.send(outcome);
    });
    let outcome = receiver
        .recv_timeout(std::time::Duration::from_secs(120))
        .expect("loading the FIFO as -ini did not finish");

    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert!(!outcome.err.contains("not readable"), "{}", outcome.err);
    writer.join().unwrap().unwrap();
    assert!(Path::new(&target).exists());
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

/// Release build, `../oracle/toppbase-completion` cases `write_cwl`,
/// `write_nested_cwl`, `write_json` and `write_nested_json`: the four CWL and
/// JSON writers end with INTERNAL_ERROR and, on standard error, exactly the
/// line the Release build prints, because it is built without TDL
/// (`ParamCWLFile.cpp:331`, `ParamJSONFile.cpp:326`). `-write_ctd` writes the
/// CTD; `tests/topp_cli_completion.rs` compares it with the Release build's.
///
/// The Release build opens the target before it fails and so leaves an empty
/// `<dir>/<tool>.cwl` or `.json` behind, emptying one that existed; this port
/// leaves the target alone (a documented native difference).
#[test]
fn cwl_and_json_writers_are_refused_as_the_release_build_refuses_them() {
    let dir = Workdir::new("description-writers");
    for (writer, extension) in [
        ("-write_cwl", "cwl"),
        ("-write_nested_cwl", "cwl"),
        ("-write_json", "json"),
        ("-write_nested_json", "json"),
    ] {
        let outcome = run::<DTAExtractor>(&["-test", writer, &text(dir.path())]);
        assert_eq!(
            outcome.code,
            ExitCode::InternalError,
            "{writer}: {}",
            outcome.err
        );
        assert_eq!(
            outcome.err,
            "Unable to initialize or run DTAExtractor: TDL support is not available. Rebuild with -DENABLE_TDL=ON to enable this feature.\n",
            "{writer}"
        );
        assert!(outcome.out.is_empty(), "{writer}: {}", outcome.out);
        assert!(
            !dir.path()
                .join(format!("DTAExtractor.{extension}"))
                .exists()
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

/// The body's streams are the caller's. After the body returns, whatever its
/// exit code, the source prints its run-time line on the info log, standard
/// output (`TOPPBase.cpp:413-424`; oracle `plain_run` in
/// `../oracle/toppbase-completion` shows its shape).
#[test]
fn run_io_output_reaches_the_callers_streams() {
    let outcome = stream_tool(Ok(ExitCode::InputFileEmpty));
    assert_eq!(outcome.code, ExitCode::InputFileEmpty);
    let took = outcome
        .out
        .strip_prefix("report line\n")
        .unwrap_or_else(|| panic!("{}", outcome.out));
    assert_took_line("StreamTool", took);
    assert_eq!(outcome.err, "diagnostic line\n");
}

/// The run-time line `<tool> took <t> (wall), <t> (CPU), <t> (system), <t>
/// (user); Peak Memory Usage: <n> MB.` and its newline, with each `<t>` in
/// `StopWatch::toString`'s seconds form and the memory part only where the
/// platform reports it.
fn assert_took_line(tool: &str, line: &str) {
    let rest = line
        .strip_prefix(&format!("{tool} took "))
        .unwrap_or_else(|| panic!("{line:?}"));
    let rest = rest
        .strip_suffix(".\n")
        .unwrap_or_else(|| panic!("{line:?}"));
    let (times, memory) = match rest.split_once("; Peak Memory Usage: ") {
        Some((times, memory)) => (times, Some(memory)),
        None => (rest, None),
    };
    let parts: Vec<&str> = times.split(", ").collect();
    assert_eq!(parts.len(), 4, "{line:?}");
    for (part, label) in parts.iter().zip(["(wall)", "(CPU)", "(system)", "(user)"]) {
        // `StopWatch::summary` renders a component the platform cannot report
        // as `n/a` (src/system/stop_watch.rs); the source always has a value.
        if *part == format!("n/a {label}") {
            continue;
        }
        let seconds = part
            .strip_suffix(&format!(" s {label}"))
            .unwrap_or_else(|| panic!("{line:?}"));
        let (whole, fraction) = seconds
            .split_once('.')
            .unwrap_or_else(|| panic!("{line:?}"));
        assert!(
            !whole.is_empty() && whole.bytes().all(|b| b.is_ascii_digit()),
            "{line:?}"
        );
        assert!(
            fraction.len() == 2 && fraction.bytes().all(|b| b.is_ascii_digit()),
            "{line:?}"
        );
    }
    if let Some(memory) = memory {
        let megabytes = memory
            .strip_suffix(" MB")
            .unwrap_or_else(|| panic!("{line:?}"));
        assert!(megabytes.bytes().all(|b| b.is_ascii_digit()), "{line:?}");
    }
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

    // Native mappings (tier 4). A std::io::Error does not say whether it read
    // or wrote, and inputs are checked before the body runs, so a denied
    // permission is taken as a failed write. The source's own InvalidRange and
    // MissingInformation exceptions would reach its BaseException arm (exit 8);
    // these follow InvalidParameter and RequiredParameterNotGiven instead.
    let denied = stream_tool(Err(Error::Io(std::io::Error::from(
        std::io::ErrorKind::PermissionDenied,
    ))));
    assert_eq!(
        denied.code,
        ExitCode::CannotWriteOutputFile,
        "{}",
        denied.err
    );
    assert!(
        denied.err.contains("Error: Unable to write file ("),
        "{}",
        denied.err
    );

    let other = stream_tool(Err(Error::Io(std::io::Error::other("device failure"))));
    assert_eq!(other.code, ExitCode::UnknownError, "{}", other.err);
    assert!(
        other.err.contains("Error: Unexpected internal error ("),
        "{}",
        other.err
    );

    let range = stream_tool(Err(Error::InvalidRange("inverted".into())));
    assert_eq!(range.code, ExitCode::IllegalParameters, "{}", range.err);

    let information = stream_tool(Err(Error::MissingInformation("no charge".into())));
    assert_eq!(
        information.code,
        ExitCode::MissingParameters,
        "{}",
        information.err
    );

    let unsorted = stream_tool(Err(Error::UnsortedData));
    assert_eq!(
        unsorted.code,
        ExitCode::IncompatibleInputData,
        "{}",
        unsorted.err
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

/// A command line at the `MAX_ARGUMENTS` bound is parsed in full, and the text
/// left after several options keeps command-line order. The source inserts each
/// such chunk at the front of its list (`TOPPBase.cpp:2436-2444`), which is
/// quadratic in the token count; the port gathers the chunks in reverse and
/// orders them once. No timing is asserted.
#[test]
fn a_command_line_at_the_argument_bound_is_parsed() {
    let outcome = run::<DTAExtractor>(&["t1", "-in", "a", "t2", "t3", "-out", "b", "t4"]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome
            .err
            .contains("Trailing text argument(s) '[t1, t2, t3, t4]' given. Aborting!"),
        "{}",
        outcome.err
    );

    let mut arguments = vec![DTAExtractor::NAME.to_owned()];
    while arguments.len() + 2 <= MAX_ARGUMENTS {
        arguments.push("-bogus".to_owned());
        arguments.push("x".to_owned());
    }
    let outcome = run_arguments::<DTAExtractor>(&arguments);
    // The diagnostic lists half a million options; keep it out of the message.
    assert_eq!(outcome.code, ExitCode::IllegalParameters);
    assert!(
        outcome
            .err
            .starts_with("Unknown option(s) '[-bogus, -bogus, ")
    );
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
    // A tool in no tool manifest reports the core version (TOPPBase.cpp:119,
    // 144-150), as ToolManifest_test.cpp's "product version is written to INI
    // with a core-version fallback" section asserts.
    assert_eq!(ctx.version(), openms::CORE_SDK_VERSION);

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
    // The version of a tool in no manifest: the core version (TOPPBase.cpp:119,
    // 559).
    assert_eq!(processing.software.version, openms::CORE_SDK_VERSION);
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
/// `ini_common_top_level`), so they are not transcribed. The section's
/// `-write_ini` comparison (483-536) is `upstream_write_ini_matches_the_retained_files`.
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

// ---------------------------------------------------------------------------
// CLI part 2 (W2.2): -write_ini parity, usage text, format checks and INI open
// failures. Oracle cases named here are in
// ../oracle/topp-cli-lifecycle/cli2/manifest.json unless marked otherwise.
// ---------------------------------------------------------------------------

fn lifecycle_path(name: &str) -> PathBuf {
    Path::new("tests/data/topp_cli_lifecycle").join(name)
}

/// Line-level comparison with the upstream comparator: `whitelist` lines are
/// skipped, numbers are compared with the given tolerances, and the log of the
/// first difference is the failure message.
fn assert_similar(produced: &Path, expected: &Path, ratio: f64, absdiff: f64, whitelist: &[&str]) {
    let mut comparator = fuzzy::FuzzyStringComparator::new();
    comparator.set_acceptable_relative(ratio);
    comparator.set_acceptable_absolute(absdiff);
    comparator.set_whitelist(whitelist.iter().map(|line| (*line).to_owned()).collect());
    comparator.set_log_destination(fuzzy::LogDestination::Buffer);
    let similar = comparator.compare_files(produced, expected);
    assert!(
        similar,
        "{} differs from {}:\n{}",
        produced.display(),
        expected.display(),
        String::from_utf8_lossy(comparator.log())
    );
}

/// `[EXTRA] getStringOption_`, option `write_ini` (TOPPBase_test.cpp:483-536).
///
/// `TEST_EQUAL(p1, p2)` is transcribed with the class test's keys and values
/// and `Param::operator==` semantics (names and values), including the
/// source's `VersionInfo::getVersion()` for `TOPPBaseTest:version`: the
/// class-test tool is in no tool manifest, so it reports the core version
/// (TOPPBase.cpp:119, 144-150), `openms::CORE_SDK_VERSION`. The file comparisons run the
/// retained C++ files `TOPPBase_test_write_ini_out.ini` and
/// `TOPPBase_test_write_ini_subsec_out.ini` (cli c19e494) through the upstream
/// comparator with `TEST_FILE_SIMILAR`'s default tolerances (absolute 1e-5,
/// relative 1 + 1e-5, `ClassTest.cpp:35-38`) and the `WHITELIST("version")`
/// the section sets, which stays in force for the second comparison.
/// Evidence: tier 1 (retained C++ output).
#[test]
fn upstream_write_ini_matches_the_retained_files() {
    let dir = Workdir::new("upstream-write-ini");
    let written = dir.file("TOPPBaseTest.ini");
    let outcome = run::<ToppBaseTest>(&["-write_ini", &written]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);

    let p1 = openms::format::paramxml::load(&written).unwrap();
    let mut p2 = Param::new();
    let strings = |values: &[&str]| values.iter().map(|v| (*v).to_owned()).collect::<Vec<_>>();
    let expected = [
        (
            "TOPPBaseTest:version",
            ParamValue::String(openms::CORE_SDK_VERSION.to_owned()),
        ),
        (
            "TOPPBaseTest:1:stringoption",
            ParamValue::String("string default".into()),
        ),
        ("TOPPBaseTest:1:intoption", ParamValue::Integer(4711)),
        ("TOPPBaseTest:1:doubleoption", ParamValue::Float(0.4711)),
        (
            "TOPPBaseTest:1:intlist",
            ParamValue::IntegerList(vec![1, 2, 3, 4]),
        ),
        (
            "TOPPBaseTest:1:doublelist",
            ParamValue::FloatList(vec![0.4711, 1.022, 4.0]),
        ),
        (
            "TOPPBaseTest:1:stringlist",
            ParamValue::StringList(strings(&["abc", "def", "ghi", "jkl"])),
        ),
        ("TOPPBaseTest:1:flag", ParamValue::String("false".into())),
        ("TOPPBaseTest:1:log", ParamValue::String(String::new())),
        ("TOPPBaseTest:1:debug", ParamValue::Integer(0)),
        ("TOPPBaseTest:1:threads", ParamValue::Integer(1)),
        (
            "TOPPBaseTest:1:no_progress",
            ParamValue::String("false".into()),
        ),
        ("TOPPBaseTest:1:force", ParamValue::String("false".into())),
        ("TOPPBaseTest:1:test", ParamValue::String("false".into())),
        (
            "TOPPBaseTest:1:stringlist2",
            ParamValue::StringList(strings(&["hopla", "dude"])),
        ),
        (
            "TOPPBaseTest:1:intlist2",
            ParamValue::IntegerList(vec![3, 4, 5]),
        ),
        (
            "TOPPBaseTest:1:doublelist2",
            ParamValue::FloatList(vec![1.2, 2.33]),
        ),
    ];
    for (key, value) in expected {
        p2.set_value(key, value, "", &[]).unwrap();
    }
    assert!(
        p1.source_equal(&p2).unwrap(),
        "written: {p1:?}\nexpected: {p2:?}"
    );
    assert_similar(
        Path::new(&written),
        &lifecycle_path("TOPPBase_test_write_ini_out.ini"),
        1.0 + 1e-5,
        1e-5,
        &["version"],
    );

    let written = dir.file("TOPPBaseCmdParseSubsectionsTest.ini");
    let outcome = run::<ToppBaseCmdParseSubsectionsTest>(&["-write_ini", &written]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_similar(
        Path::new(&written),
        &lifecycle_path("TOPPBase_test_write_ini_subsec_out.ini"),
        1.0 + 1e-5,
        1e-5,
        &["version"],
    );
}

/// `-test -write_ini` of one ported tool against the C++ file.
///
/// The registered `TOPPWRITEINI_<tool>` test runs the tool, and
/// `TOPPWRITEINI_<tool>_SectionName` (test-data `topp/CMakeLists.txt:83-85`,
/// `check_ini.cmake`) requires the first line matching `^  <NODE name="…"` to
/// name the tool. Beyond that, the file must match the product SDK's file for
/// the same command (oracle `write_ini_<tool>`): line by line with exact numbers
/// and only `version` lines skipped, as `TOPPWRITEINI_OVERWRITE` compares
/// (CMakeLists.txt:104-106), and as decoded parameter trees, entry for entry,
/// with descriptions, tags, restrictions and supported formats.
fn assert_write_ini_matches_the_cpp_file<T: Tool>() {
    let dir = Workdir::new("write-ini-parity");
    let written = dir.file(&format!("{}.tmp.ini", T::NAME));
    let outcome = run::<T>(&["-test", "-write_ini", &written]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);

    let bytes = fs::read(&written).unwrap();
    let content = String::from_utf8(bytes).unwrap();
    assert!(
        content.starts_with("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"),
        "{content}"
    );
    let section = content
        .lines()
        .find_map(|line| line.strip_prefix("  <NODE name=\""))
        .and_then(|rest| rest.split('"').next());
    assert_eq!(section, Some(T::NAME), "{content}");

    let expected = lifecycle_path(&format!("write_ini_{}.ini", T::NAME));
    assert_similar(Path::new(&written), &expected, 1.0, 0.0, &["version"]);
    let produced = openms::format::paramxml::load(&written).unwrap();
    let oracle = openms::format::paramxml::load(&expected).unwrap();
    assert_eq!(produced, oracle, "{}", T::NAME);
}

#[test]
fn write_ini_matches_the_cpp_file_for_every_ported_tool() {
    assert_write_ini_matches_the_cpp_file::<BaselineFilter>();
    assert_write_ini_matches_the_cpp_file::<DTAExtractor>();
    assert_write_ini_matches_the_cpp_file::<MapNormalizer>();
    assert_write_ini_matches_the_cpp_file::<MzMLSplitter>();
    assert_write_ini_matches_the_cpp_file::<SpectraFilterWindowMower>();
}

/// `--help` and `--helphelp` against the C++ usage text (oracle `help_<tool>`
/// and `helphelp_<tool>`), byte for byte on the error stream, with nothing on
/// the output stream. The oracle ran without a terminal and without `COLUMNS`,
/// so the source wrote plain text without line shaping; its first line, the
/// `stty: stdin isn't a terminal` message of the console-width probe, is not
/// part of the retained text. The version line names the revision of the
/// build: the product SDK prints `4fdec46`, this port its pinned core revision.
fn assert_usage_matches_the_cpp_text<T: Tool>() {
    let cpp_version = "1.0.0 (OpenMS core 4.0.0, revision 4fdec46)";
    for (flag, prefix) in [("--help", "help"), ("--helphelp", "helphelp")] {
        let outcome = run::<T>(&[flag]);
        assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
        assert!(outcome.out.is_empty(), "{}", outcome.out);
        let expected = fs::read_to_string(lifecycle_path(&format!("{prefix}_{}.txt", T::NAME)))
            .unwrap()
            .replace(cpp_version, &verbose_version::<T>());
        assert_eq!(outcome.err, expected, "{} {flag}", T::NAME);
    }
}

#[test]
fn usage_text_matches_the_cpp_text_for_every_ported_tool() {
    assert_eq!(
        verbose_version::<DTAExtractor>(),
        format!(
            "1.0.0 (OpenMS core 4.0.0, revision {})",
            &openms::CORE_SDK_REVISION[..7]
        )
    );
    assert_usage_matches_the_cpp_text::<BaselineFilter>();
    assert_usage_matches_the_cpp_text::<DTAExtractor>();
    assert_usage_matches_the_cpp_text::<MapNormalizer>();
    assert_usage_matches_the_cpp_text::<MzMLSplitter>();
    assert_usage_matches_the_cpp_text::<SpectraFilterWindowMower>();
}

/// SpectraFilterWindowMower on `input` into `output`, both inside `dir`.
fn swm_between(dir: &Workdir, input: &str, output: &str) -> (Outcome, String) {
    let out = dir.file(output);
    let outcome = run::<SpectraFilterWindowMower>(&["-test", "-in", input, "-out", &out]);
    (outcome, out)
}

/// A copy of the window mower input under another name.
fn swm_input_as(dir: &Workdir, name: &str) -> String {
    let path = dir.file(name);
    fs::copy(fixture("window_mower_tool_input.mzML"), &path).unwrap();
    path
}

const UNDETERMINED_FORMAT: &str = "Warning: Could not determine format of input file";

/// The input format comes from `FileHandler::get_type`, by name and then by
/// content (TOPPBase.cpp:1575-1593). Oracle `in_no_extension_mzml_content`,
/// `in_unknown_extension_mzml_content` and `in_uppercase_extension`: an mzML
/// file without an extension, with an unknown one, or with `.MZML` runs, exit
/// 0, with no warning. Oracle `in_txt_extension_mzml_content` and
/// `in_dta_extension`: the name decides a known type first, so mzML content
/// named `.txt` or `.dta` is refused, exit 6.
#[test]
fn the_input_format_is_detected_by_name_then_content() {
    let dir = Workdir::new("input-format");
    for name in [
        "window_mower_input",
        "window_mower_input.foo",
        "window_mower_input.MZML",
    ] {
        let input = swm_input_as(&dir, name);
        let (outcome, out) = swm_between(&dir, &input, &format!("{name}.out.mzML"));
        assert_eq!(
            outcome.code,
            ExitCode::ExecutionOk,
            "{name}: {}",
            outcome.err
        );
        assert!(
            !outcome.err.contains(UNDETERMINED_FORMAT),
            "{}",
            outcome.err
        );
        assert_same_peaks(&out, fixture("window_mower_tool_output.mzML"));
    }
    for (name, format) in [
        ("window_mower_input.txt", "txt"),
        ("window_mower_input.dta", "dta"),
    ] {
        let input = swm_input_as(&dir, name);
        let (outcome, out) = swm_between(&dir, &input, &format!("{name}.out.mzML"));
        assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
        assert!(
            outcome.err.contains(&format!(
                "Invalid parameter: Input file '{input}' has invalid format '{format}'. Valid formats are: 'mzML'."
            )),
            "{}",
            outcome.err
        );
        assert!(!Path::new(&out).exists());
    }
}

/// Oracle `in_unknown_name_and_content`: a file whose name and content are both
/// unknown only warns in the format check (TOPPBase.cpp:1580-1583), and the run
/// continues to the load, which refuses it; nothing is written. The C++ load
/// failure is a `ParseError`, exit 3 ("type: unknown is not allowed for loading
/// an experiment"); this port's `FileHandler::load_experiment` reports
/// `Error::InvalidValue`, exit 6. That difference belongs to the loader
/// (`src/format/file_handler.rs`) and is recorded in `docs/TOPP_CLI_SUPPORT.md`;
/// the exit code is therefore only asserted to be a failure.
#[test]
fn an_undetermined_input_format_only_warns() {
    let dir = Workdir::new("input-unknown");
    let input = dir.file("unknown.foo");
    fs::write(&input, b"neither a known name nor known content\n").unwrap();
    let (outcome, out) = swm_between(&dir, &input, "out.mzML");
    assert!(
        outcome
            .err
            .contains(&format!("{UNDETERMINED_FORMAT} '{input}'!")),
        "{}",
        outcome.err
    );
    assert_ne!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_ne!(outcome.code, ExitCode::InputFileNotFound, "{}", outcome.err);
    assert!(!Path::new(&out).exists());
}

/// The output format comes from the file name alone (TOPPBase.cpp:1595-1608).
/// Oracle `out_unknown_extension`, `out_no_extension` and
/// `out_uppercase_extension`: a name no file type claims, no extension, or
/// `.MZML` is accepted and written, exit 0. Oracle `out_txt_extension`: a known
/// type the parameter does not accept is refused, exit 6.
#[test]
fn the_output_format_is_checked_by_name_only() {
    let dir = Workdir::new("output-format");
    for name in ["out.foo", "out", "out.MZML"] {
        let (outcome, out) = swm_between(&dir, &swm_input(), name);
        assert_eq!(
            outcome.code,
            ExitCode::ExecutionOk,
            "{name}: {}",
            outcome.err
        );
        assert_same_peaks(&out, fixture("window_mower_tool_output.mzML"));
    }
    let (outcome, out) = swm_between(&dir, &swm_input(), "out.txt");
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    assert!(
        outcome.err.contains(&format!(
            "Invalid parameter: Invalid output file extension for file '{out}'. Valid file extensions are: 'mzML'."
        )),
        "{}",
        outcome.err
    );
    assert!(!Path::new(&out).exists());
}

/// Run SpectraFilterWindowMower with `ini` before a run and with `-write_ini`,
/// and assert the source's `IOException` outcome on both paths.
fn assert_ini_open_failure_is_an_unexpected_internal_error(dir: &Workdir, ini: &str) {
    let out = dir.file("out.mzML");
    let before_run =
        run::<SpectraFilterWindowMower>(&["-test", "-ini", ini, "-in", &swm_input(), "-out", &out]);
    let written_ini = dir.file("written.ini");
    let with_write_ini =
        run::<SpectraFilterWindowMower>(&["-write_ini", &written_ini, "-ini", ini]);
    let expected = format!("Error: Unexpected internal error (IO error for file '{ini}')");
    for outcome in [&before_run, &with_write_ini] {
        assert_eq!(outcome.code, ExitCode::UnknownError, "{}", outcome.err);
        assert!(outcome.err.contains(&expected), "{}", outcome.err);
    }
    assert!(!Path::new(&out).exists());
    assert!(!Path::new(&written_ini).exists());
}

/// Oracle `ini_socket` and `write_ini_ini_socket`
/// (`../oracle/topp-cli-lifecycle/ini_read_failures/manifest.json`): a Unix
/// socket given as `-ini` exists and is readable by mode, but opening it fails.
/// `TextFile` throws `IOException` for such a file (TextFile.cpp:44-47), which
/// is exit 8 with `Error: Unexpected internal error (IO error for file '…')`
/// (TOPPBase.cpp:495-499), before a run and with `-write_ini`. Skipped where
/// the socket cannot be bound or opens for reading.
#[cfg(unix)]
#[test]
fn a_socket_ini_is_an_unexpected_internal_error() {
    use std::os::unix::net::UnixListener;
    let dir = Workdir::new("ini-socket");
    let socket = dir.file("ini.sock");
    match UnixListener::bind(&socket) {
        Ok(listener) => drop(listener),
        Err(error) => {
            eprintln!("skipped: cannot bind {socket}: {error}");
            return;
        }
    }
    if fs::File::open(&socket).is_ok() {
        eprintln!("skipped: {socket} opens for reading here");
        return;
    }
    assert_ini_open_failure_is_an_unexpected_internal_error(&dir, &socket);
    eprintln!("ran: socket -ini {socket}");
}

/// Oracle `ini_tty` and `write_ini_ini_tty`
/// (`../oracle/topp-cli-lifecycle/ini_read_failures/manifest.json`): `/dev/tty`
/// opened by a process without a controlling terminal fails, so it is exit 8 as
/// for the socket, on both paths. The oracle started the C++ tool in a new
/// session; this case runs only when the test process itself has no
/// controlling terminal, which holds under the gate and in CI, and is skipped
/// otherwise, where the load would wait for terminal input.
#[cfg(unix)]
#[test]
fn a_terminal_ini_without_a_controlling_terminal_is_an_unexpected_internal_error() {
    let tty = "/dev/tty";
    match fs::File::open(tty) {
        Ok(_) => {
            eprintln!("skipped: {tty} opens, this process has a controlling terminal");
            return;
        }
        Err(error)
            if matches!(
                error.kind(),
                std::io::ErrorKind::NotFound | std::io::ErrorKind::PermissionDenied
            ) =>
        {
            eprintln!("skipped: {tty}: {error}");
            return;
        }
        Err(_) => {}
    }
    let dir = Workdir::new("ini-tty");
    assert_ini_open_failure_is_an_unexpected_internal_error(&dir, tty);
    eprintln!("ran: terminal -ini {tty}");
}
