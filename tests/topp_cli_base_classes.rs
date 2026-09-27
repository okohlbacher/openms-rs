// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The cli base classes against the C++ Release build: `TOPPExternalToolBase`,
//! `SearchEngineBase`, `TOPPMapAlignerBase`, `INIUpdater` and the
//! `TOPPBase_defs.h` exceptions of `OpenMS4-cli` `c19e494`.
//!
//! Evidence tier 1. The oracle `../oracle/B7-cli` (lane B7, 2026-09-27) links a
//! driver with one small tool per base class against the Linux x86_64 Release
//! install at the pins (core `bc9cc12`, cli `c19e494`, topp `174b576`) on
//! ibminode06, and runs 65 cases twice in clean environments; both runs agree
//! once run times are masked. Its cases are retained under
//! `tests/data/topp_cli_base_classes/cases` with the inputs under `inputs`,
//! and each is replayed here in process with the same tool defined natively,
//! the Release run's paths mapped to the replay's. Provenance:
//! `tests/data/topp_cli_base_classes_provenance.json`.
//!
//! The class tests `TOPPExternalToolBase_test.cpp`, `SearchEngineBase_test.cpp`
//! and `INIUpdater_test.cpp` are transcribed as well (tier 3); their cases
//! are among the oracle's.
#![cfg(feature = "paramxml")]

use openms::cli::{
    ExitCode, INIUpdater, ReferenceParameterKind, TOPPMapAlignerBase, Tool, ToolContext, ToolError,
    ToolHandler, ToolResult, ToolSpec, exit_code_of, run_with,
};
use openms::data_structures::ToolDescriptionInternal;
use openms::param::{Param, ParamValue};
use openms::system::external_process::ReturnState;
use openms::system::file::TempDir;
use std::collections::BTreeMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Where the Release run kept its inputs, its per-case directories and the
/// shared data it was given as `OPENMS_DATA_PATH`.
const RELEASE_INPUTS: &str = "/scratch/kohlbach/oracle/B7-cli/data";
const RELEASE_RUN: &str = "/scratch/kohlbach/oracle/B7-cli/runs/run1";
#[cfg(feature = "mzml")]
const RELEASE_DATA: &str =
    "/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576/share/OpenMS/4.0.0";

fn retained() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/topp_cli_base_classes")
}

fn release(case: &str, file: &str) -> String {
    let path = retained().join("cases").join(case).join(file);
    fs::read_to_string(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// The closing `<tool> took … .` line with its figures masked, as the retained
/// standard output has it.
fn mask_took(tool: &str, out: &str) -> String {
    let prefix = format!("{tool} took ");
    out.lines()
        .map(|line| {
            if line.starts_with(&prefix) && line.ends_with('.') {
                format!("{tool} took <T>.")
            } else {
                line.to_owned()
            }
        })
        .map(|line| line + "\n")
        .collect()
}

/// One replayed case.
struct Replay {
    code: ExitCode,
    out: String,
    err: String,
    cwd: PathBuf,
    /// The Release run's text with its paths replaced by the replay's.
    map: Box<dyn Fn(&str) -> String>,
    _work: TempDir,
    /// A further directory the replay needs as long as it lives.
    _extra: Option<TempDir>,
}

impl Replay {
    /// The Release run's exit status.
    fn expected_code(&self, case: &str) -> i32 {
        release(case, "exit_code.txt").trim().parse().unwrap()
    }
    fn expected(&self, case: &str, file: &str) -> String {
        (self.map)(&release(case, file))
    }
    /// Exit code, standard output (closing line masked) and standard error
    /// equal the Release run's.
    fn assert_streams(&self, case: &str, tool: &str) {
        assert_eq!(
            self.code.as_i32(),
            self.expected_code(case),
            "{case}: exit code\n{}",
            self.err
        );
        assert_eq!(
            mask_took(tool, &self.out),
            self.expected(case, "stdout.txt"),
            "{case}: standard output"
        );
        assert_eq!(
            self.err,
            self.expected(case, "stderr.txt"),
            "{case}: standard error"
        );
    }
}

/// Run `T` with the Release case's arguments, the inputs staged in a fresh
/// directory. `prepare` sees the replay's case directory and path map before
/// the run.
fn replay<T: Tool>(case: &str, prepare: impl FnOnce(&Path, &dyn Fn(&str) -> String)) -> Replay {
    let work = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let inputs = work.path().join("inputs");
    let case_dir = work.path().join("case");
    let cwd = case_dir.join("cwd");
    fs::create_dir_all(&cwd).unwrap();
    fs::create_dir_all(case_dir.join("home")).unwrap();
    copy_dir(&retained().join("inputs"), &inputs);
    let (input_text, case_text) = (
        inputs.to_str().unwrap().to_owned(),
        case_dir.to_str().unwrap().to_owned(),
    );
    let release_case = format!("{RELEASE_RUN}/{case}");
    let map = move |text: &str| {
        text.replace(RELEASE_INPUTS, &input_text)
            .replace(&release_case, &case_text)
    };
    prepare(&case_dir, &map);
    let mut arguments: Vec<String> = release(case, "argv.txt").lines().map(&map).collect();
    arguments[0] = T::NAME.to_owned();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<T>(&arguments, &mut out, &mut err);
    Replay {
        code,
        out: String::from_utf8(out).unwrap(),
        err: String::from_utf8(err).unwrap(),
        cwd,
        map: Box::new(map),
        _work: work,
        _extra: None,
    }
}

fn copy_dir(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_dir(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), &target).unwrap();
        }
    }
}

fn no_preparation(_: &Path, _: &dyn Fn(&str) -> String) {}

// ---------------------------------------------------------------------------
// TOPPExternalToolBase
// ---------------------------------------------------------------------------

/// The oracle's `B7ExternalTool`: runs `-exe` with `-args` (an `@` in front of
/// an argument is dropped, so that arguments starting with `-` reach the
/// child) in `-workdir` with the `-env` overlay, writes the captured streams
/// to `<-captured>.stdout` and `.stderr`, and returns the run's exit code.
struct B7ExternalTool;

impl Tool for B7ExternalTool {
    const NAME: &'static str = "B7ExternalTool";
    const DESCRIPTION: &'static str = "B7 oracle for TOPPExternalToolBase";

    fn register(spec: &mut ToolSpec) -> openms::Result<()> {
        spec.register_string_option("exe", "<path>", "", "executable", true, false)?;
        spec.register_string_list("args", "<list>", &[], "arguments", false, false)?;
        spec.register_string_option("workdir", "<dir>", "", "working directory", false, false)?;
        spec.register_string_list("env", "<K=V>", &[], "environment overlay", false, false)?;
        spec.register_string_option(
            "captured",
            "<prefix>",
            "",
            "where to write the captured streams",
            false,
            false,
        )?;
        spec.register_flag(
            "first_overload",
            "use the overload without captured streams",
            false,
        )
    }

    fn run(_ctx: &ToolContext) -> ToolResult {
        unreachable!("run_io is overridden")
    }

    fn run_io(ctx: &ToolContext, out: &mut dyn Write, err: &mut dyn Write) -> ToolResult {
        let env: BTreeMap<String, String> = ctx
            .string_list("env")?
            .iter()
            .map(|kv| {
                let (key, value) = kv.split_once('=').unwrap_or((kv, ""));
                (key.to_owned(), value.to_owned())
            })
            .collect();
        let args: Vec<String> = ctx
            .string_list("args")?
            .iter()
            .map(|a| a.strip_prefix('@').unwrap_or(a).to_owned())
            .collect();
        let exe = ctx.string("exe")?;
        let workdir = ctx.string("workdir")?;
        if ctx.flag("first_overload")? {
            return Ok(ctx.run_external_process(out, err, exe, &args, workdir, &env)?);
        }
        let run = ctx.run_external_process_output(out, err, exe, &args, workdir, &env)?;
        let captured = ctx.string("captured")?;
        if !captured.is_empty() {
            fs::write(format!("{captured}.stdout"), &run.stdout)?;
            fs::write(format!("{captured}.stderr"), &run.stderr)?;
        }
        Ok(run.code)
    }
}

fn ext(case: &str) -> Replay {
    replay::<B7ExternalTool>(case, no_preparation)
}

/// The captured streams equal the Release run's.
fn assert_captured(case: &str, replay: &Replay) {
    for stream in ["stdout", "stderr"] {
        let produced = fs::read_to_string(replay.cwd.join(format!("cap.{stream}"))).unwrap();
        assert_eq!(
            produced,
            replay.expected(case, &format!("cwd/cap.{stream}")),
            "{case}: captured {stream}"
        );
    }
}

/// `TOPPExternalToolBase_test.cpp`, section `runExternalProcess_` (tier 3):
/// a program that does not exist is `EXTERNAL_PROGRAM_NOTFOUND`, one that
/// fails is `EXTERNAL_PROGRAM_ERROR`, one that succeeds `EXECUTION_OK`.
#[cfg(unix)]
#[test]
fn class_test_run_external_process_maps_the_three_outcomes() {
    for (case, code) in [
        ("ext_missing_first", ExitCode::ExternalProgramNotFound),
        ("ext_broken_first", ExitCode::ExternalProgramError),
        ("ext_ok_first", ExitCode::ExecutionOk),
    ] {
        let replay = ext(case);
        assert_eq!(replay.code, code, "{case}: {}", replay.err);
        assert_eq!(replay.code.as_i32(), replay.expected_code(case), "{case}");
    }
}

/// Every outcome of a child whose output is the same on every platform, at
/// debug levels 0 and 4, byte for byte as the Release build: a program that
/// cannot start, a non-zero exit with output on both streams, a crash, output
/// without a final line break, an environment overlay, a working directory,
/// and a repeated line, which the info stream counts at debug level 4.
#[cfg(unix)]
#[test]
fn runs_equal_the_release_build() {
    for case in [
        "ext_missing",
        "ext_missing_debug4",
        "ext_echo_debug4",
        "ext_repeat_debug4",
        "ext_env",
        "ext_workdir",
        "ext_both_streams",
        "ext_crash",
        "ext_no_newline",
        "ext_repeat_fail",
    ] {
        let replay = ext(case);
        replay.assert_streams(case, B7ExternalTool::NAME);
        assert_captured(case, &replay);
    }
    let replay = ext("ext_missing_first");
    replay.assert_streams("ext_missing_first", B7ExternalTool::NAME);
}

/// `ls -0` and `ls -ld .`: the child's own words differ between platforms
/// (GNU and BSD `ls`) and, on Linux, in the program name the child prints,
/// because the source starts the program by the path `PATH` resolves it to
/// (`bp::search_path`) and this port by its name. The structure around them
/// is the Release build's: the banner, the two `Standard …:` records, the
/// process message and the exit code; at debug level 4 everything goes to
/// standard output instead, and with `-log` the two records reach the log.
#[cfg(unix)]
#[test]
fn a_failing_platform_program_is_reported_as_the_release_build_reports_it() {
    for case in ["ext_broken", "ext_broken_first", "ext_broken_log"] {
        let replay = ext(case);
        assert_eq!(replay.code.as_i32(), replay.expected_code(case), "{case}");
        assert!(
            replay
                .err
                .starts_with("Standard output: Running: ls -0\n\nStandard error: "),
            "{case}: {}",
            replay.err
        );
        assert!(
            replay
                .err
                .contains("\nProcess 'ls' did not finish successfully (exit code: ")
                && replay.err.ends_with("). Please check the log.\n\n"),
            "{case}: {}",
            replay.err
        );
        assert_eq!(
            mask_took(B7ExternalTool::NAME, &replay.out),
            replay.expected(case, "stdout.txt"),
            "{case}"
        );
    }
    let replay = ext("ext_broken_log");
    let log = fs::read_to_string(replay.cwd.join("run.log")).unwrap();
    let records: Vec<&str> = log
        .lines()
        .filter(|line| line.contains(" B7ExternalTool:1:: "))
        .map(|line| line.split_once(" B7ExternalTool:1:: ").unwrap().1)
        .collect();
    assert_eq!(
        records.first().copied(),
        Some("Standard output: Running: ls -0"),
        "{log}"
    );
    assert!(
        records[1..]
            .iter()
            .any(|r| r.starts_with("Standard error: ")),
        "{log}"
    );
    let expected = replay.expected("ext_broken_log", "run.log.masked");
    assert!(expected.starts_with("<DATE> B7ExternalTool:1:: Standard output: Running: ls -0\n\n<DATE> B7ExternalTool:1:: Standard error: "));

    let replay = ext("ext_broken_debug4");
    assert_eq!(replay.code, ExitCode::ExternalProgramError);
    assert_eq!(
        replay.err, "",
        "at debug level 4 nothing reaches standard error"
    );
    assert!(replay.out.starts_with("Running: ls -0\n"), "{}", replay.out);
    assert!(
        mask_took(B7ExternalTool::NAME, &replay.out)
            .ends_with("). Please check the log.\nB7ExternalTool took <T>.\n"),
        "{}",
        replay.out
    );

    let replay = ext("ext_ok");
    assert_eq!(replay.code, ExitCode::ExecutionOk, "{}", replay.err);
    assert_eq!(replay.err, "");
    let captured = fs::read_to_string(replay.cwd.join("cap.stdout")).unwrap();
    assert!(captured.starts_with("Running: ls -ld .\n"), "{captured}");
    assert!(
        captured.ends_with("Executed 'ls' successfully!\n"),
        "{captured}"
    );
}

#[test]
fn every_return_state_is_the_source_exit_code() {
    assert_eq!(exit_code_of(ReturnState::Success), ExitCode::ExecutionOk);
    assert_eq!(
        exit_code_of(ReturnState::NonzeroExit),
        ExitCode::ExternalProgramError
    );
    assert_eq!(
        exit_code_of(ReturnState::Crash),
        ExitCode::ExternalProgramError
    );
    assert_eq!(
        exit_code_of(ReturnState::FailedToStart),
        ExitCode::ExternalProgramNotFound
    );
}

// ---------------------------------------------------------------------------
// SearchEngineBase
// ---------------------------------------------------------------------------

#[cfg(feature = "mzml")]
mod search_engine {
    use super::*;
    use openms::cli::SearchEngineBase;
    use openms::format::paramxml;
    use openms::system::file::FileContext;
    use std::cell::RefCell;

    thread_local! {
        /// The file context the replayed tool resolves its database with.
        static FILES: RefCell<Option<FileContext>> = const { RefCell::new(None) };
    }

    /// The indexer defaults the Release build's `PeptideIndexing` has,
    /// retained from the oracle (`peptide_indexing_defaults.ini`); the native
    /// `PeptideIndexing` has no parameter surface yet.
    pub(super) fn peptide_indexing_defaults() -> Param {
        paramxml::load(retained().join("peptide_indexing_defaults.ini")).unwrap()
    }

    /// The oracle's `B7SearchEngine`: `SearchEngineBase_test.cpp`'s test tool
    /// with `-in` also accepting mzXML and MGF, an `-ms_level` option and the
    /// `PeptideIndexing` parameters registered.
    pub(super) struct B7SearchEngine;

    impl Tool for B7SearchEngine {
        const NAME: &'static str = "B7SearchEngine";
        const DESCRIPTION: &'static str = "B7 oracle for SearchEngineBase";

        fn register(spec: &mut ToolSpec) -> openms::Result<()> {
            spec.register_input_file("in", "<file>", "", "Input file", true, false, &[])?;
            spec.set_valid_formats("in", &["mzML", "mzXML", "mgf"])?;
            spec.register_input_file(
                "database",
                "<file>",
                "",
                "FASTA file",
                true,
                false,
                &["skipexists"],
            )?;
            spec.set_valid_formats("database", &["FASTA"])?;
            spec.register_int_option("ms_level", "<n>", 2, "MS level to check", false, false)?;
            SearchEngineBase::register_peptide_indexing_parameter(spec, peptide_indexing_defaults())
        }

        fn run(_ctx: &ToolContext) -> ToolResult {
            unreachable!("run_io is overridden")
        }

        fn run_io(ctx: &ToolContext, out: &mut dyn Write, err: &mut dyn Write) -> ToolResult {
            let ms_level = i32::try_from(ctx.int("ms_level")?).unwrap();
            let input = SearchEngineBase::get_rawfile_name(ctx, err, ms_level)?;
            let db = FILES.with(|files| {
                let files = files.borrow();
                SearchEngineBase::get_db_filename(ctx, out, err, files.as_ref().unwrap(), "")
            })?;
            writeln!(out, "in={input}\ndb={db}")?;
            Ok(ExitCode::ExecutionOk)
        }
    }

    /// Replay a search-engine case with a file context whose shared data is
    /// a fresh directory named as `OPENMS_DATA_PATH`, and whose home holds
    /// the case's `OpenMS.ini` when the Release run had one.
    pub(super) fn se(case: &str) -> Replay {
        let data = TempDir::new_in(std::env::temp_dir(), false).unwrap();
        let data_text = data.path().to_str().unwrap().to_owned();
        let mut replay = replay::<B7SearchEngine>(case, |case_dir, map| {
            let mut files = FileContext::new(
                case_dir.join("bin"),
                case_dir.join("home"),
                std::env::temp_dir(),
            )
            .unwrap();
            files.data_override = Some(data.path().to_path_buf());
            let ini = retained()
                .join("cases")
                .join(case)
                .join("home/.config/OpenMS/OpenMS.ini");
            if ini.exists() {
                let config = files.get_openms_config_dir().unwrap();
                fs::create_dir_all(&config).unwrap();
                fs::write(
                    config.join("OpenMS.ini"),
                    map(&fs::read_to_string(ini).unwrap()),
                )
                .unwrap();
            }
            FILES.with(|slot| *slot.borrow_mut() = Some(files));
        });
        FILES.with(|slot| *slot.borrow_mut() = None);
        let map = std::mem::replace(&mut replay.map, Box::new(|text: &str| text.to_owned()));
        replay.map = Box::new(move |text: &str| map(text).replace(RELEASE_DATA, &data_text));
        // The shared-data directory lives as long as the replay.
        replay._extra = Some(data);
        replay
    }

    /// `SearchEngineBase_test.cpp`, section `getRawfileName` (tier 3): an mzML
    /// file with only an unknown-type spectrum and one with profile spectra
    /// end with `UNKNOWN_ERROR`, and both run with `-force`.
    #[test]
    fn class_test_get_rawfile_name() {
        for (case, code) in [
            ("se_empty", ExitCode::UnknownError),
            ("se_empty_force", ExitCode::ExecutionOk),
            ("se_profile", ExitCode::UnknownError),
            ("se_profile_force", ExitCode::ExecutionOk),
        ] {
            let replay = se(case);
            assert_eq!(replay.code, code, "{case}: {}", replay.err);
        }
    }

    /// Every spectra check and database resolution, byte for byte as the
    /// Release build: no spectrum of the level (`FileEmpty`, levels 1, 0 and
    /// -1), profile or no centroided spectra with and without `-force`,
    /// centroided MS1 data, MGF without a check, mzXML with the warning the
    /// source leaves unflushed (at the end of standard error, after a failed
    /// run's error lines), a database found through `OpenMS.ini`'s
    /// `id_db_dir` (also in a subdirectory) and one found nowhere.
    #[test]
    fn runs_equal_the_release_build() {
        for case in [
            "se_empty",
            "se_empty_force",
            "se_profile",
            "se_profile_force",
            "se_profile_ms1",
            "se_profile_ms0",
            "se_profile_msneg",
            "se_mower_ms1",
            "se_mower_ms2",
            "se_dta_ms1",
            "se_dta_ms2",
            "se_mgf",
            "se_mzxml",
            "se_mzxml_twice",
            "se_db_missing",
            "se_db_iddir",
            "se_db_iddir_sub",
            "se_db_iddir_missing",
        ] {
            se(case).assert_streams(case, B7SearchEngine::NAME);
        }
    }

    /// The `-log` file receives the catch block's line, and not the lines the
    /// source's `findDatabase` writes to its error log stream.
    #[test]
    fn the_log_file_holds_the_catch_block_line_only() {
        let case = "se_db_missing_log";
        let replay = se(case);
        replay.assert_streams(case, B7SearchEngine::NAME);
        let log = fs::read_to_string(replay.cwd.join("run.log")).unwrap();
        let masked: String = log
            .lines()
            .map(|line| match line.split_once(" B7SearchEngine:1:: ") {
                Some((_, text)) => format!("<DATE> B7SearchEngine:1:: {text}\n"),
                None => format!("{line}\n"),
            })
            .collect();
        assert_eq!(masked, replay.expected(case, "run.log.masked"));
    }

    /// A blank `-database` and one of another format fail in the framework's
    /// file checks, before the tool body: this port validates every
    /// parameter before the body where the source checks each when the body
    /// first reads it (`docs/TOPP_CLI_SUPPORT.md`, *Native differences*), so
    /// the body's `-force` warning, which the Release build prints first, is
    /// not printed. The error line and the exit code are the Release build's.
    #[test]
    fn database_values_the_framework_refuses_end_as_in_the_release_build() {
        for case in ["se_db_blank", "se_db_data", "se_db_data_dotdot"] {
            let replay = se(case);
            assert_eq!(
                replay.code.as_i32(),
                replay.expected_code(case),
                "{case}: {}",
                replay.err
            );
            let expected = replay.expected(case, "stderr.txt");
            let (warning, error) = expected.split_once('\n').unwrap();
            assert!(warning.starts_with("Warning: No centroided MS2"), "{case}");
            assert_eq!(replay.err, error, "{case}");
            assert_eq!(replay.out, replay.expected(case, "stdout.txt"), "{case}");
        }
    }

    /// `registerPeptideIndexingParameter_` with the Release build's indexer
    /// defaults writes the Release build's INI file: `-reindex`, and the
    /// `PeptideIndexing:` section with `missing_decoy_action` reset to a bare
    /// `warn` entry and the twelve advanced tags.
    #[test]
    fn the_peptide_indexing_parameters_are_the_release_builds() {
        let case = "se_write_ini";
        let replay = se(case);
        assert_eq!(replay.code, ExitCode::ExecutionOk, "{}", replay.err);
        let produced = fs::read_to_string(replay.cwd.join("out.ini")).unwrap();
        assert_eq!(produced, replay.expected(case, "cwd/out.ini"));
    }

    /// The indexer defaults are patched as the source patches them, whatever
    /// the caller passed, and a tree without one of the twelve names is
    /// refused (the source's `Param::addTag` throws `ElementNotFound`).
    #[test]
    fn register_peptide_indexing_parameter_patches_and_checks_the_tree() {
        let mut spec = ToolSpec::new();
        SearchEngineBase::register_peptide_indexing_parameter(
            &mut spec,
            peptide_indexing_defaults(),
        )
        .unwrap();
        let reindex = spec.find("reindex").unwrap();
        assert_eq!(reindex.valid_strings, vec!["true", "false"]);
        let action = spec.find("PeptideIndexing:missing_decoy_action").unwrap();
        assert_eq!(action.default_value, ParamValue::String("warn".into()));
        assert!(action.advanced);
        assert!(
            action.valid_strings.is_empty(),
            "setValue drops the restrictions"
        );
        assert!(spec.find("PeptideIndexing:enzyme:name").unwrap().advanced);
        assert!(
            !spec
                .find("PeptideIndexing:allow_nterm_protein_cleavage")
                .unwrap()
                .advanced
        );

        let mut short = peptide_indexing_defaults();
        short.remove("IL_equivalent").unwrap();
        let mut spec = ToolSpec::new();
        assert!(SearchEngineBase::register_peptide_indexing_parameter(&mut spec, short).is_err());
    }
}

// ---------------------------------------------------------------------------
// TOPPMapAlignerBase
// ---------------------------------------------------------------------------

/// The oracle's map-aligner tools: the common registration for featureXML and
/// consensusXML with the kind of reference parameters `KIND`,
/// `checkParameters_`, then `transformSpectraFiles_` with `y = x + 10` for
/// every input.
trait AlignerKind {
    const NAME: &'static str;
    const BASE: TOPPMapAlignerBase;
}

struct B7MapAligner<K>(std::marker::PhantomData<K>);
struct NoneKind;
struct RestrictedKind;
struct FlexibleKind;
impl AlignerKind for NoneKind {
    const NAME: &'static str = "B7MapAlignerNone";
    const BASE: TOPPMapAlignerBase = TOPPMapAlignerBase::new(ReferenceParameterKind::None);
}
impl AlignerKind for RestrictedKind {
    const NAME: &'static str = "B7MapAlignerRestricted";
    const BASE: TOPPMapAlignerBase = TOPPMapAlignerBase::new(ReferenceParameterKind::Restricted);
}
impl AlignerKind for FlexibleKind {
    const NAME: &'static str = "B7MapAlignerFlexible";
    const BASE: TOPPMapAlignerBase = TOPPMapAlignerBase::new(ReferenceParameterKind::Flexible);
}

impl<K: AlignerKind> Tool for B7MapAligner<K> {
    const NAME: &'static str = K::NAME;
    const DESCRIPTION: &'static str = "B7 oracle for TOPPMapAlignerBase";

    fn register(spec: &mut ToolSpec) -> openms::Result<()> {
        K::BASE.register_options_and_flags_map_aligners(spec, &["featureXML", "consensusXML"])?;
        spec.register_flag(
            "store_original_rt",
            "keep the original retention times",
            false,
        )
    }

    fn run(_ctx: &ToolContext) -> ToolResult {
        unreachable!("run_io is overridden")
    }

    fn run_io(ctx: &ToolContext, _out: &mut dyn Write, err: &mut dyn Write) -> ToolResult {
        let code = K::BASE.check_parameters(ctx, err)?;
        if code != ExitCode::ExecutionOk {
            return Ok(code);
        }
        #[cfg(feature = "mzml")]
        {
            use openms::analysis::transformations::{
                DataPoint, LinearOptions, ModelConfig, TransformationDescription,
            };
            let mut trafos = Vec::new();
            for _ in ctx.string_list("in")? {
                let mut trafo = TransformationDescription::new(vec![
                    DataPoint::new(0.0, 10.0),
                    DataPoint::new(100.0, 110.0),
                ])?;
                trafo.fit_model(ModelConfig::Linear(LinearOptions::default()))?;
                trafos.push(trafo);
            }
            K::BASE.transform_spectra_files(
                ctx,
                ctx.string_list("in_spectra_files")?,
                ctx.string_list("out_spectra_files")?,
                &trafos,
                ctx.flag("store_original_rt")?,
            )?;
        }
        Ok(ExitCode::ExecutionOk)
    }
}

type AlignerNone = B7MapAligner<NoneKind>;
type AlignerRestricted = B7MapAligner<RestrictedKind>;
type AlignerFlexible = B7MapAligner<FlexibleKind>;

/// Every check of `checkParameters_` in the source's order, each with the
/// Release build's line and exit code, and the reference checks of the
/// restricted and flexible kinds; a tool without reference parameters does
/// not know `-reference:index`, and a negative index fails the registered
/// minimum before the body.
#[test]
fn check_parameters_equals_the_release_build() {
    for case in [
        "ma_none_no_out",
        "ma_none_out_count",
        "ma_none_trafo_count",
        "ma_none_mixed",
        "ma_none_spectra_in",
        "ma_none_spectra_out",
        "ma_none_spectra_only_out",
        "ma_none_ok",
    ] {
        replay::<AlignerNone>(case, no_preparation).assert_streams(case, AlignerNone::NAME);
    }
    for case in [
        "ma_restricted_index_high",
        "ma_restricted_index_eq",
        "ma_restricted_index_neg",
        "ma_restricted_both",
        "ma_restricted_type",
        "ma_restricted_file_ok",
    ] {
        replay::<AlignerRestricted>(case, no_preparation)
            .assert_streams(case, AlignerRestricted::NAME);
    }
    let case = "ma_flexible_type";
    replay::<AlignerFlexible>(case, no_preparation).assert_streams(case, AlignerFlexible::NAME);
}

/// The registration of each kind writes the Release build's INI file,
/// including the `reference` section and its description.
#[test]
fn the_registration_writes_the_release_ini_files() {
    fn check<T: Tool>(case: &str) {
        let replay = replay::<T>(case, no_preparation);
        assert_eq!(replay.code, ExitCode::ExecutionOk, "{case}: {}", replay.err);
        let produced = fs::read_to_string(replay.cwd.join("out.ini")).unwrap();
        assert_eq!(produced, replay.expected(case, "cwd/out.ini"), "{case}");
    }
    check::<AlignerNone>("ma_none_write_ini");
    check::<AlignerRestricted>("ma_restricted_write_ini");
    check::<AlignerFlexible>("ma_flexible_write_ini");
}

/// `transformSpectraFiles_` writes the Release build's spectra: every
/// retention time moved by the transformation, the original kept as meta
/// value when asked, and one `Retention time alignment` record of the tool.
#[cfg(feature = "mzml")]
#[test]
fn transformed_spectra_files_equal_the_release_build() {
    use openms::format::file_handler::FileHandler;
    use openms::format::file_types::FileType;
    for case in ["ma_none_transform", "ma_none_transform_store"] {
        let replay = replay::<AlignerNone>(case, no_preparation);
        assert_eq!(replay.code, ExitCode::ExecutionOk, "{case}: {}", replay.err);
        assert_eq!(replay.err, "");
        for name in ["a.mzML", "b.mzML"] {
            let produced =
                FileHandler::load_experiment(replay.cwd.join(name), &[FileType::MzMl]).unwrap();
            let expected = FileHandler::load_experiment(
                retained().join("cases").join(case).join("cwd").join(name),
                &[FileType::MzMl],
            )
            .unwrap();
            assert_eq!(
                produced.spectra.len(),
                expected.spectra.len(),
                "{case} {name}"
            );
            for (i, (a, e)) in produced.spectra.iter().zip(&expected.spectra).enumerate() {
                let at = format!("{case} {name} spectrum {i}");
                assert_eq!(a.native_id, e.native_id, "{at}");
                assert_eq!(a.rt, e.rt, "{at}: retention time");
                assert_eq!(a.metadata, e.metadata, "{at}: metadata");
                assert_eq!(a.peaks, e.peaks, "{at}: peaks");
                assert_eq!(a.data_processing.len(), e.data_processing.len(), "{at}");
                for (pa, pe) in a.data_processing.iter().zip(&e.data_processing) {
                    assert_eq!(pa.software.name, pe.software.name, "{at}");
                    assert_eq!(pa.software.version, pe.software.version, "{at}");
                    assert_eq!(pa.actions, pe.actions, "{at}");
                    assert_eq!(pa.completion_time, pe.completion_time, "{at}");
                    assert_eq!(pa.metadata, pe.metadata, "{at}");
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// INIUpdater
// ---------------------------------------------------------------------------

/// The registry of the built-in product manifest alone: the pinned `topp`
/// package's `topp.tools.tsv`, byte-identical to the Release install's.
fn builtin_registry() -> ToolHandler {
    let mut sources = openms::cli::ToolRegistrySources::only_prefixes(Vec::new()).unwrap();
    sources.builtin_manifest = true;
    ToolHandler::new(sources)
}

/// `INIUpdater_test.cpp` (tier 3): the constructor, `getToolNamesFromINI`
/// on an empty tree, two `:version` entries and a `:version2` one,
/// `getNameMapping`'s `FeatureFinder -type centroided` entry, and
/// `getNewToolName` for mapped, registered, typed and unknown names.
#[test]
fn class_test_ini_updater() {
    let updater = INIUpdater::new();
    let mut p = Param::new();
    assert!(updater.get_tool_names_from_ini(&p).unwrap().is_empty());
    p.set_value(
        "FeatureFinder:version",
        ParamValue::String("1.9".into()),
        "",
        &[],
    )
    .unwrap();
    p.set_value(
        "SomeTool:version",
        ParamValue::String("whatever".into()),
        "",
        &[],
    )
    .unwrap();
    assert_eq!(updater.get_tool_names_from_ini(&p).unwrap().len(), 2);
    p.set_value(
        "BrokenTool:version2",
        ParamValue::String("1.9".into()),
        "",
        &[],
    )
    .unwrap();
    assert_eq!(updater.get_tool_names_from_ini(&p).unwrap().len(), 2);

    let map = updater.get_name_mapping();
    assert_ne!(map.len(), 0);
    assert_eq!(
        map.get(&ToolDescriptionInternal::with_types(
            "FeatureFinder",
            &["centroided".to_owned()]
        )),
        Some(&ToolDescriptionInternal::with_types(
            "FeatureFinderCentroided",
            &[]
        ))
    );

    let registry = builtin_registry();
    let name = |old: &str, kind: &str| updater.get_new_tool_name(old, kind, &registry).unwrap();
    assert_eq!(
        name("FeatureFinder", "centroided").as_deref(),
        Some("FeatureFinderCentroided")
    );
    assert_eq!(
        name("PeakPicker", "wavelet").as_deref(),
        Some("PeakPickerWavelet")
    );
    assert_eq!(name("FileInfo", "").as_deref(), Some("FileInfo"));
    assert_eq!(name("FileInfo", "bogus type").as_deref(), Some("FileInfo"));
    assert_eq!(name("UNKNOWNTOOL", "bogus type"), None);
}

/// Every answer of the Release build's `INIUpdater`
/// (`tests/data/topp_cli_base_classes/iniupdater.tsv`): the tool names of
/// four trees, the whole table in map order with its type-list lengths, and
/// nineteen lookups, with the built-in product manifest as the registry.
#[test]
fn ini_updater_equals_the_release_build() {
    let updater = INIUpdater::new();
    let registry = builtin_registry();
    let mut produced = String::new();
    let mut p = Param::new();
    let mut names = |label: &str, p: &Param| {
        produced.push_str(&format!("names\t{label}"));
        for name in updater.get_tool_names_from_ini(p).unwrap() {
            produced.push_str(&format!("\t{name}"));
        }
        produced.push('\n');
    };
    names("empty", &p);
    let text = |v: &str| ParamValue::String(v.into());
    p.set_value("FeatureFinder:version", text("1.9"), "", &[])
        .unwrap();
    p.set_value("SomeTool:version", text("whatever"), "", &[])
        .unwrap();
    names("two", &p);
    p.set_value("BrokenTool:version2", text("1.9"), "", &[])
        .unwrap();
    names("broken", &p);
    p.set_value("Nested:1:version", text("1.9"), "", &[])
        .unwrap();
    p.set_value("Aaa:version", text("1"), "", &[]).unwrap();
    p.set_value("version", text("top"), "", &[]).unwrap();
    p.set_value("x:version:y", text("no"), "", &[]).unwrap();
    names("more", &p);
    for entry in updater.get_name_mapping().iter() {
        produced.push_str(&format!(
            "map\t{}\t{}\t{}\t{}\t{}\t{}\n",
            entry.old.name,
            entry.old.types.join(","),
            entry.old.types.len(),
            entry.new.name,
            entry.new.types.join(","),
            entry.new.types.len()
        ));
    }
    let expected = fs::read_to_string(retained().join("iniupdater.tsv")).unwrap();
    for line in expected.lines().filter(|line| line.starts_with("new\t")) {
        let fields: Vec<&str> = line.split('\t').collect();
        let (old, kind) = (fields[1], fields[2]);
        let answer = updater.get_new_tool_name(old, kind, &registry).unwrap();
        produced.push_str(&format!(
            "new\t{old}\t{kind}\t{}\t{}\n",
            u8::from(answer.is_some()),
            answer.unwrap_or_default()
        ));
    }
    assert_eq!(produced, expected);
}

// ---------------------------------------------------------------------------
// TOPPBase_defs.h
// ---------------------------------------------------------------------------

/// The three exceptions a tool throws itself end the run as the source's
/// catch blocks end it (`TOPPBase.cpp:466-493`): `RequiredParameterNotGiven`
/// quotes its text unless it starts with a quote and returns
/// `MISSING_PARAMETERS`; `UnregisteredParameter` and `WrongParameterType`
/// return `INTERNAL_ERROR`. None prints the closing line.
#[test]
fn the_defs_exceptions_end_the_run_as_their_catch_blocks() {
    struct Thrower;
    thread_local! {
        static THROW: std::cell::Cell<u8> = const { std::cell::Cell::new(0) };
    }
    impl Tool for Thrower {
        const NAME: &'static str = "B7Thrower";
        const DESCRIPTION: &'static str = "throws a TOPPBase_defs exception";
        fn register(_spec: &mut ToolSpec) -> openms::Result<()> {
            Ok(())
        }
        fn run(_ctx: &ToolContext) -> ToolResult {
            Err(match THROW.with(std::cell::Cell::get) {
                0 => ToolError::required_parameter_not_given("in"),
                1 => ToolError::required_parameter_not_given("'in' [valid: mzML]"),
                2 => ToolError::unregistered_parameter("nope"),
                _ => ToolError::wrong_parameter_type("count"),
            })
        }
    }
    for (which, code, line) in [
        (
            0,
            ExitCode::MissingParameters,
            "Error: The required parameter 'in' was not given or is empty!\n",
        ),
        (
            1,
            ExitCode::MissingParameters,
            "Error: The required parameter 'in' [valid: mzML] was not given or is empty!\n",
        ),
        (
            2,
            ExitCode::InternalError,
            "Internal error: Request for unregistered parameter 'nope'\n",
        ),
        (
            3,
            ExitCode::InternalError,
            "Internal error: Request for parameter with wrong type 'count'\n",
        ),
    ] {
        THROW.with(|slot| slot.set(which));
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let result = run_with::<Thrower>(&["B7Thrower".to_owned()], &mut out, &mut err);
        assert_eq!(result, code, "{which}");
        assert_eq!(String::from_utf8(err).unwrap(), line, "{which}");
        assert!(
            !String::from_utf8(out).unwrap().contains(" took "),
            "{which}"
        );
    }
}
