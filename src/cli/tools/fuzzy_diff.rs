// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Compares two files, tolerating numeric differences.
//!
//! Ports `OpenMS4-topp/src/FuzzyDiff.cpp` (topp `174b576`) on the TOPP
//! lifecycle, over the library comparator
//! [`crate::concept::fuzzy_string_comparator`], exactly as the source tool
//! drives `OpenMS::FuzzyStringComparator`. See [`super`] on why the tool lives
//! in the library rather than in the binary, and `docs/TOPP_FUZZY_DIFF_SUPPORT.md`
//! for the parameter table, the exit codes and the evidence.
//!
//! The source's page documentation, carried over: in the diff output,
//! "position" refers to the characters in the string, whereas "column" is meant
//! for the text editor. Only one of `ratio` or `absdiff` has to be satisfied;
//! use `absdiff` to deal with cases like "zero vs. epsilon".
//!
//! The tool body follows the source `main_` (`OpenMS4-topp/src/FuzzyDiff.cpp:88-208`):
//!
//! 1. The parameters are read. The framework has already checked both inputs
//!    as `getStringOption_` does (a missing file exits 1, an unreadable one 2,
//!    an empty regular file 4) and the registered ranges (6).
//! 2. Each `-matched_whitelist` entry must split at `:` into exactly two
//!    parts; the first that does not is the source's `IllegalArgument`,
//!    reported as `Error: Unexpected internal error (<entry> does not have the
//!    format String1:String2)` with `UNKNOWN_ERROR` (8).
//! 3. The comparator is configured with the tolerances, both whitelists, the
//!    verbose level, the tab width and the first column.
//! 4. The two files are compared, or with `-sort` their texts with every line
//!    but the first sorted ([`sorted_lines`]). The comparator's report goes to
//!    the output stream, as the source's goes to `std::cout`.
//! 5. The tool exits `EXECUTION_OK` (0) when no difference was found and
//!    `PARSE_ERROR` (10) otherwise; the source's `TODO` about better exit
//!    codes stands.
//!
//! Three native differences, each documented at the item that makes it:
//! `-sort` compares in memory instead of through temporary files, so the report
//! names the inputs rather than the deleted copies; an input that cannot be
//! read exits `INTERNAL_ERROR` (12) in both modes, where the source does so only
//! without `-sort`; and an input beyond the comparator's
//! [`MAX_INPUT_BYTES`] is refused with `INCOMPATIBLE_INPUT_DATA` (11).

use crate::cli::{ExitCode, Tool, ToolContext, ToolSpec};
use crate::concept::fuzzy_string_comparator::{
    FuzzyStringComparator, InputFailure, LogDestination, MAX_INPUT_BYTES, parse_matched_whitelist,
    sorted_lines,
};
use crate::{Error, Result};
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

/// The `FuzzyDiff` TOPP tool.
///
/// Registration and the comparison run; every comparison rule is the library
/// comparator's.
pub struct FuzzyDiff;

impl Tool for FuzzyDiff {
    const NAME: &'static str = "FuzzyDiff";
    const DESCRIPTION: &'static str = "Compares two files, tolerating numeric differences.";

    /// Source `registerOptionsAndFlags_` (`OpenMS4-topp/src/FuzzyDiff.cpp:54-86`),
    /// in source order with the source texts.
    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.add_empty_line()?;
        spec.register_input_file("in1", "<file>", "", "first input file", true, false, &[])?;
        spec.register_input_file("in2", "<file>", "", "second input file", true, false, &[])?;
        spec.add_empty_line()?;
        spec.register_double_option(
            "ratio",
            "<double>",
            1.0,
            r#"acceptable relative error. Only one of 'ratio' or 'absdiff' has to be satisfied.  Use "absdiff" to deal with cases like "zero vs. epsilon"."#,
            false,
            false,
        )?;
        spec.set_min_float("ratio", 1.0)?;
        spec.register_double_option(
            "absdiff",
            "<double>",
            0.0,
            "acceptable absolute difference. Only one of 'ratio' or 'absdiff' has to be satisfied. ",
            false,
            false,
        )?;
        spec.set_min_float("absdiff", 0.0)?;
        spec.add_empty_line()?;
        spec.register_string_list(
            "whitelist",
            "<string list>",
            &["<?xml-stylesheet"],
            "Lines containing one of these strings are skipped",
            false,
            true,
        )?;
        // Source default `ListUtils::create<std::string>("")`, which splits the
        // empty string into no elements: an empty list.
        spec.register_string_list(
            "matched_whitelist",
            "<string list>",
            &[],
            "Lines where one file contains one string and the other file another string are skipped. Input is given as list of colon separated tuples, e.g. String1:String2 String3:String4",
            false,
            true,
        )?;
        spec.register_int_option(
            "verbose",
            "<int>",
            2,
            "set verbose level:\n0 = very quiet mode (absolutely no output)\n1 = quiet mode (no output unless differences detected)\n2 = default (include summary at end)\n3 = continue after errors\n",
            false,
            false,
        )?;
        spec.set_min_int("verbose", 0)?;
        spec.set_max_int("verbose", 3)?;
        spec.register_int_option(
            "tab_width",
            "<int>",
            8,
            "tabulator width, used for calculation of column numbers",
            false,
            false,
        )?;
        spec.set_min_int("tab_width", 1)?;
        spec.register_int_option(
            "first_column",
            "<int>",
            1,
            "number of first column, used for calculation of column numbers",
            false,
            false,
        )?;
        spec.set_min_int("first_column", 0)?;
        spec.add_empty_line()?;
        spec.register_flag(
            "sort",
            "sort the input files before comparison (useful for tabular files where row order may vary). The first line of each file is assumed to be a header and is not sorted.",
            false,
        )?;
        Ok(())
    }

    /// Run against the process streams; see `run_io`.
    fn run(ctx: &ToolContext) -> Result<ExitCode> {
        Self::run_io(ctx, &mut std::io::stdout(), &mut std::io::stderr())
    }

    /// Source `main_` (`OpenMS4-topp/src/FuzzyDiff.cpp:88-208`).
    ///
    /// `out` receives the comparator's report, which the source writes to
    /// `std::cout`; `err` receives the diagnostics. The body runs on the
    /// calling thread: the comparison is serial in the source too, so
    /// `-threads` has nothing to size, and no worker pool is built for it
    /// (see [`ToolContext::in_thread_pool`]).
    ///
    /// The source's two debug lines (`writeDebug_` of both whitelists at
    /// `-debug 1`) are not written, because the framework does not port
    /// `writeDebug_`, and neither is its trailing timing line.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when a stream cannot be written, and
    /// [`Error::InvalidValue`] for an integer parameter outside the `i32`
    /// range the source's `getIntOption_` returns. Every other outcome is an
    /// exit code.
    fn run_io(ctx: &ToolContext, out: &mut dyn Write, err: &mut dyn Write) -> Result<ExitCode> {
        let in1 = ctx.string("in1")?.to_owned();
        let in2 = ctx.string("in2")?.to_owned();
        let ratio = ctx.double("ratio")?;
        let absdiff = ctx.double("absdiff")?;
        let whitelist = ctx.string_list("whitelist")?.to_vec();
        let raw_matched_whitelist = ctx.string_list("matched_whitelist")?.to_vec();
        let verbose_level = int_option(ctx, "verbose")?;
        let tab_width = int_option(ctx, "tab_width")?;
        let first_column = int_option(ctx, "first_column")?;
        let do_sort = ctx.flag("sort")?;

        // The source throws IllegalArgument, which TOPPBase's BaseException
        // arm reports (TOPPBase.cpp:495-499 at cli c19e494).
        let matched_whitelist = match parse_matched_whitelist(&raw_matched_whitelist) {
            Ok(pairs) => pairs,
            Err(Error::InvalidValue(message)) => {
                writeln!(err, "Error: Unexpected internal error ({message})")?;
                return Ok(ExitCode::UnknownError);
            }
            Err(other) => return Err(other),
        };

        let mut fsc = FuzzyStringComparator::new();
        fsc.set_log_destination(LogDestination::Buffer);
        fsc.set_acceptable_relative(ratio);
        fsc.set_acceptable_absolute(absdiff);
        fsc.set_whitelist(whitelist);
        fsc.set_matched_whitelist(matched_whitelist);
        fsc.set_verbose_level(verbose_level);
        fsc.set_tab_width(tab_width);
        fsc.set_first_column(first_column);

        let result = if do_sort {
            let text_1 = match read_for_sort(&in1, err)? {
                Ok(text) => text,
                Err(code) => return Ok(code),
            };
            let text_2 = match read_for_sort(&in2, err)? {
                Ok(text) => text,
                Err(code) => return Ok(code),
            };
            fsc.set_input_names(&in1, &in2);
            fsc.compare_bytes(&sorted_lines(&text_1), &sorted_lines(&text_2))
        } else {
            fsc.compare_files(Path::new(&in1), Path::new(&in2))
        };

        if let Some(failure) = fsc.input_failure().cloned() {
            return report_input_failure(&fsc, &failure, out, err);
        }
        out.write_all(fsc.log())?;
        warn_if_truncated(&fsc, err)?;
        Ok(if result {
            ExitCode::ExecutionOk
        } else {
            // Source: "TODO think about better exit codes."
            ExitCode::ParseError
        })
    }
}

/// An integer parameter as the source's `getIntOption_` returns it, an `Int`.
///
/// The registered ranges bound `verbose`; `tab_width` and `first_column` have
/// no upper bound, and an INI value can exceed what a command-line value
/// can. Such a value is refused rather than wrapped.
fn int_option(ctx: &ToolContext, name: &str) -> Result<i32> {
    let value = ctx.int(name)?;
    i32::try_from(value).map_err(|_| {
        Error::InvalidValue(format!(
            "Invalid value '{value}' for integer parameter '{name}' given. Out of the range of a 32-bit integer."
        ))
    })
}

/// Read one input for `-sort`, within the comparator's [`MAX_INPUT_BYTES`].
///
/// Source `sortFile` opens the file with `std::ifstream` and throws
/// `FileNotFound` when that fails, which `TOPPBase` reports as `Error: File
/// not found (the file '<name>' could not be found)` with
/// `INPUT_FILE_NOT_FOUND` (1), whatever the reason; the framework's input
/// check makes that a race here, as there.
///
/// Native differences: a file larger than [`MAX_INPUT_BYTES`] is refused with
/// `INCOMPATIBLE_INPUT_DATA` (11), because the text is held in memory to be
/// sorted, where the source has no limit. A read that fails is refused with
/// `INTERNAL_ERROR` (12), the code the source gives a failed read without
/// `-sort`: the source's `std::getline` swallows the error and sorts what it
/// read, so a directory compares as an empty text (executed on the Release
/// build, `sort_directory`, exit 10), and so would a file cut short by an I/O
/// error. This port does not compare a text it could not read in full.
fn read_for_sort(
    path: &str,
    err: &mut dyn Write,
) -> Result<std::result::Result<Vec<u8>, ExitCode>> {
    let Ok(file) = File::open(path) else {
        writeln!(
            err,
            "Error: File not found (the file '{path}' could not be found)"
        )?;
        return Ok(Err(ExitCode::InputFileNotFound));
    };
    let too_large = || {
        format!(
            "Error: input file '{path}' exceeds the comparison limit of {MAX_INPUT_BYTES} bytes."
        )
    };
    if file.metadata().map(|m| m.len()).unwrap_or(0) > MAX_INPUT_BYTES {
        writeln!(err, "{}", too_large())?;
        return Ok(Err(ExitCode::IncompatibleInputData));
    }
    let mut text = Vec::new();
    if let Err(error) = file.take(MAX_INPUT_BYTES + 1).read_to_end(&mut text) {
        let description = match InputFailure::from_io_error(&error) {
            InputFailure::Read { description, .. } => description,
            InputFailure::TooLarge => error.to_string(),
        };
        writeln!(
            err,
            "Unable to initialize or run FuzzyDiff: error reading the file '{path}': {description}"
        )?;
        return Ok(Err(ExitCode::InternalError));
    }
    if text.len() as u64 > MAX_INPUT_BYTES {
        writeln!(err, "{}", too_large())?;
        return Ok(Err(ExitCode::IncompatibleInputData));
    }
    Ok(Ok(text))
}

/// Report a comparison that stopped before it had read all of its input.
///
/// The reports written before the failure go to `out`, as the source has
/// already written them to `std::cout` when its read fails.
///
/// * A failed read: the source's line reader lets `std::ios_base::failure`
///   escape `compareFiles`, and `TOPPBase`'s last catch reports `Unable to
///   initialize or run FuzzyDiff: <what>` with `INTERNAL_ERROR` (12)
///   (`TOPPBase.cpp:510-513` at cli c19e494). `<what>` is libstdc++'s text,
///   `basic_filebuf::underflow error reading the file: <strerror>`, which is
///   what the Linux Release build prints for a directory given as an input
///   (executed: `directory_vs_file`, `file_vs_directory`,
///   `directory_vs_directory`); it is reproduced so the stream matches.
/// * The comparator's [`MAX_INPUT_BYTES`]: native, the source reads without a
///   limit. The comparator's own message goes to `err` and the tool exits
///   `INCOMPATIBLE_INPUT_DATA` (11), not the 10 that would claim a difference.
fn report_input_failure(
    fsc: &FuzzyStringComparator,
    failure: &InputFailure,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<ExitCode> {
    let kept = fsc.log_without_input_failure();
    out.write_all(kept)?;
    warn_if_truncated(fsc, err)?;
    match failure {
        InputFailure::Read { description, .. } => {
            writeln!(
                err,
                "Unable to initialize or run FuzzyDiff: basic_filebuf::underflow error reading the file: {description}"
            )?;
            Ok(ExitCode::InternalError)
        }
        InputFailure::TooLarge => {
            let log = fsc.log();
            err.write_all(log.get(kept.len()..).unwrap_or_default())?;
            Ok(ExitCode::IncompatibleInputData)
        }
    }
}

/// Say so on `err` when the comparator's report reached its buffer limit.
///
/// Native: the source writes its report straight to `std::cout` without a
/// limit. Only verbose level 3 reports every difference, so only there can the
/// report reach [`MAX_LOG_BYTES`](crate::concept::fuzzy_string_comparator::MAX_LOG_BYTES);
/// the verdict and the exit code are unaffected.
fn warn_if_truncated(fsc: &FuzzyStringComparator, err: &mut dyn Write) -> Result<()> {
    if fsc.log_truncated() {
        writeln!(
            err,
            "Warning: the comparison report was cut at {} bytes; the verdict is unaffected.",
            crate::concept::fuzzy_string_comparator::MAX_LOG_BYTES
        )?;
    }
    Ok(())
}
