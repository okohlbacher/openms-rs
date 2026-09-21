// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Shows basic information about the file, such as data ranges and file type.
//!
//! Ports `OpenMS4-topp/src/FileInfo.cpp` (topp `174b576`) as a thin wrapper
//! over the FileInfo library, `crate::format::file_info::report::FileInfo`,
//! exactly as the source tool delegates to `OpenMS::FileInfo`. See [`super`]
//! on why the tool lives in the library rather than in the binary, and
//! `docs/TOPP_FILE_INFO_SUPPORT.md` for the capability table, the exit codes
//! and the evidence.
//!
//! The tool body follows the source `main_` and `outputTo_`:
//!
//! 1. `-out` and `-out_tsv` are opened, creating or truncating them, before
//!    anything is read, so a run that fails later leaves them empty; without
//!    `-out` the report goes to the output stream, and without `-out_tsv` the
//!    TSV report is discarded.
//! 2. The type is `-in_type`, or else detected from the file name and then the
//!    content; an undetermined type prints `Error: Could not determine input
//!    file type!` and exits `PARSE_ERROR` (10).
//! 3. `-i` on anything but mzML prints `Error: Can only validate indices for
//!    mzML files` and the usage text, and exits `ILLEGAL_PARAMETERS` (6).
//! 4. The flags become library options with the resolved type forced, and the
//!    library's text and TSV reports are written.
//! 5. A failed `-i` index check exits `ILLEGAL_PARAMETERS`; everything else
//!    exits `EXECUTION_OK`.
//!
//! A branch or flag the library does not run yet fails with
//! [`Error::Unsupported`] naming it, which the framework reports as `Error:
//! unsupported: ...` with `INCOMPATIBLE_INPUT_DATA` (11), before anything is
//! written to the reports (decisions D3 and D5 of the early TOPP bundle). The
//! source reports those branches and exits 0.

use crate::cli::{ExitCode, Tool, ToolContext, ToolError, ToolResult, ToolSpec};
use crate::format::file_info::model::Options;
use crate::format::file_info::report::FileInfo as FileInfoLibrary;
use crate::format::{FileHandler, FileType};
use crate::system::file;
use crate::{Error, Result};
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// The `FileInfo` TOPP tool.
///
/// Registration, option mapping and output routing only; every report branch
/// is the library's, `crate::format::file_info::report::FileInfo`.
pub struct FileInfo;

/// The input formats `-in` accepts and the types `-in_type` names, in source
/// order (`OpenMS4-topp/src/FileInfo.cpp:85`).
const INPUT_TYPES: [&str; 17] = [
    "mzData",
    "mzXML",
    "mzML",
    "sqMass",
    "dta",
    "dta2d",
    "mgf",
    "featureXML",
    "consensusXML",
    "idXML",
    "pepXML",
    "mzTab",
    "fid",
    "mzid",
    "trafoXML",
    "fasta",
    "pqp",
];

impl Tool for FileInfo {
    const NAME: &'static str = "FileInfo";
    const DESCRIPTION: &'static str =
        "Shows basic information about the file, such as data ranges and file type.";

    /// Source `registerOptionsAndFlags_` (`OpenMS4-topp/src/FileInfo.cpp:83-100`), with the
    /// descriptions verbatim; the usage text capitalises their first letter.
    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_input_file("in", "<file>", "", "input file", true, false, &[])?;
        spec.set_valid_formats("in", &INPUT_TYPES)?;
        spec.register_string_option(
            "in_type",
            "<type>",
            "",
            "input file type -- default: determined from file extension or content",
            false,
            false,
        )?;
        spec.set_valid_strings("in_type", &INPUT_TYPES)?;
        spec.register_output_file(
            "out",
            "<file>",
            "",
            "Optional output file. If left out, the output is written to the command line.",
            false,
            false,
        )?;
        spec.set_valid_formats("out", &["txt"])?;
        spec.register_output_file(
            "out_tsv",
            "<file>",
            "",
            "Second optional output file. Tab separated flat text file.",
            false,
            true,
        )?;
        spec.set_valid_formats("out_tsv", &["tsv"])?;
        spec.register_flag(
            "m",
            "Show meta information about the whole experiment",
            false,
        )?;
        spec.register_flag("p", "Shows data processing information", false)?;
        spec.register_flag(
            "s",
            "Computes a five-number statistics of intensities, qualities, and widths",
            false,
        )?;
        spec.register_flag(
            "d",
            "Show detailed listing of all spectra and chromatograms (peak files only)",
            false,
        )?;
        spec.register_flag(
            "c",
            "Check for corrupt data in the file (peak files only)",
            false,
        )?;
        spec.register_flag(
            "v",
            "Validate the file only (for mzML, mzData, mzXML, featureXML, idXML, consensusXML, pepXML)",
            false,
        )?;
        spec.register_flag(
            "i",
            "Check whether a given mzML file contains valid indices (conforming to the indexedmzML standard)",
            false,
        )?;
        Ok(())
    }

    /// Run against the process streams; see `run_io`.
    fn run(ctx: &ToolContext) -> ToolResult {
        Self::run_io(ctx, &mut std::io::stdout(), &mut std::io::stderr())
    }

    /// Source `main_` and `outputTo_` (`OpenMS4-topp/src/FileInfo.cpp:102-199`).
    ///
    /// `out` receives the report when `-out` is not given, as the source's
    /// info log does; `err` receives the diagnostics, the usage text of the
    /// `-i` refusal and the library's warnings, which the source writes with
    /// `OPENMS_LOG_WARN`. The two refusals the source writes with
    /// `writeLogError_` — an undetermined type and `-i` on anything but mzML —
    /// also reach the `-log` file. Both are exit codes `outputTo_` returns, so
    /// the lifecycle's closing `FileInfo took … .` line follows them, as in
    /// the Release build (oracle `in_is_directory`).
    ///
    /// # Errors
    ///
    /// A report file that cannot be opened is the source's `FileNotWritable`,
    /// [`ToolError::unexpected`] (8), which ends the run without the closing
    /// line (oracle `out_is_directory`). The library's errors, which the
    /// framework maps to exit codes: [`Error::Parse`] for a corrupt input (3),
    /// [`Error::Unsupported`] for an unported branch or flag (11),
    /// [`Error::InvalidValue`] for a detected type other than the forced one
    /// (6, where the source exits 3 or 8; see the support document), and
    /// [`Error::Io`] when a report cannot be written.
    fn run_io(ctx: &ToolContext, out: &mut dyn Write, err: &mut dyn Write) -> ToolResult {
        let input = ctx.string("in")?.to_owned();
        let mut text_file = open_output(ctx.string("out")?)?;
        let mut tsv_file = open_output(ctx.string("out_tsv")?)?;

        // Source `FileTypes::nameToType`: the empty default names no type.
        let mut in_type = FileType::from_name(ctx.string("in_type")?);
        if in_type == FileType::Unknown {
            in_type = detect_type(&input)?;
            ctx.write_debug(&format!("Input file type: {}", in_type.name()), 2);
        }
        if in_type == FileType::Unknown {
            ctx.write_log_error(err, "Error: Could not determine input file type!")?;
            return Ok(ExitCode::ParseError);
        }
        if ctx.flag("i")? && in_type != FileType::MzMl {
            ctx.write_log_error(err, "Error: Can only validate indices for mzML files")?;
            let spec = crate::cli::tool_spec::<Self>()?;
            let subsections = crate::cli::subsection_defaults::<Self>(&spec)?;
            crate::cli::print_usage::<Self>(err, &spec, &subsections, false)?;
            return Ok(ExitCode::IllegalParameters);
        }

        let options = Options {
            forced_type: in_type,
            meta: ctx.flag("m")?,
            processing: ctx.flag("p")?,
            statistics: ctx.flag("s")?,
            detailed: ctx.flag("d")?,
            check_corrupt: ctx.flag("c")?,
            validate: ctx.flag("v")?,
            check_index: ctx.flag("i")?,
            log_type: ctx.progress_log_type(),
            // The source loader reads dangling mzML header references
            // (decision D10); the library default stays strict.
            source_dangling_references: true,
        };
        let result = FileInfoLibrary::new().run(&input, &options)?;
        // The library's `OPENMS_LOG_WARN` lines: the warning log stream,
        // yellow on a terminal, and not the -log file.
        for warning in &result.warnings {
            crate::cli::log_warning(err, warning)?;
        }

        match text_file.as_mut() {
            Some(handle) => handle.write_all(FileInfoLibrary::to_text(&result).as_bytes())?,
            None => out.write_all(FileInfoLibrary::to_text(&result).as_bytes())?,
        }
        if let Some(handle) = tsv_file.as_mut() {
            handle.write_all(FileInfoLibrary::to_tsv(&result).as_bytes())?;
        }

        if options.check_index && result.validation.index_checked && !result.validation.index_valid
        {
            return Ok(ExitCode::IllegalParameters);
        }
        Ok(ExitCode::ExecutionOk)
    }
}

/// Open one report file as the source's `std::ofstream::open` does, creating
/// or truncating it; an empty `path` opens nothing.
///
/// The framework has already confirmed that the path is writable, as the
/// source's `outputFileWritable_` has; an open that still fails, such as on an
/// existing directory, is the source's thrown `FileNotWritable`, which
/// reaches the `BaseException` arm of `TOPPBase::main`: `UNKNOWN_ERROR` (8)
/// with `Error: Unexpected internal error (the file '<path>' is not writable
/// for the current user)`, and no closing line (oracle `out_is_directory`).
fn open_output(path: &str) -> std::result::Result<Option<File>, ToolError> {
    if path.is_empty() {
        return Ok(None);
    }
    File::create(path).map(Some).map_err(|_| {
        ToolError::unexpected(format!(
            "the file '{path}' is not writable for the current user"
        ))
    })
}

/// Source `FileHandler::getType(in)`, by file name and then by content.
///
/// The source content check reads a directory as an empty file and returns
/// `UNKNOWN` (oracle `in_is_directory`, exit 10), where
/// [`FileHandler::get_type`] returns the I/O error of reading it; only that
/// case is mapped, as the library's own type detection does. Every other error
/// is returned.
fn detect_type(input: &str) -> Result<FileType> {
    match FileHandler::get_type(Path::new(input)) {
        Err(Error::Io(_)) if file::is_directory(input) => Ok(FileType::Unknown),
        other => other,
    }
}
