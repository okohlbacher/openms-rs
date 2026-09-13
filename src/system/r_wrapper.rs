// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Calling R scripts, mainly to produce plots, from `SYSTEM/RWrapper.h`.
//!
//! Three entry points, in the order a caller uses them: is the interpreter
//! present ([`find_r`](crate::system::r_wrapper::find_r)), where does the
//! script live ([`find_script`](crate::system::r_wrapper::find_script)), and
//! run it ([`run_script`](crate::system::r_wrapper::run_script)).
//!
//! The graceful-degradation contract the class test pins is preserved: a
//! missing interpreter or a missing script is reported as a `false` outcome or
//! a not-found error, never as a panic.
//!
//! Similar modules exist for other external tools, for example
//! [`java_info`](crate::system::java_info) and
//! [`python_info`](crate::system::python_info).
//!
//! See `docs/R_WRAPPER_SUPPORT.md` for the full API mapping.

use crate::system::external_process::{Invocation, IoMode, ReturnState, capture};
use crate::system::file::FileContext;
use crate::{Error, Result};
use std::path::PathBuf;

/// The interpreter the source defaults to, `Rscript`.
///
/// Rust has no default arguments, so the source's `executable = "Rscript"` is
/// this constant and every caller passes it explicitly.
pub const DEFAULT_EXECUTABLE: &str = "Rscript";

/// Subdirectory of the OpenMS shared data that holds bundled R scripts.
pub const SCRIPTS_DIRECTORY: &str = "SCRIPTS";

/// Arguments the source hands `Rscript` to probe for a working interpreter:
/// `--vanilla -e sessionInfo()`.
pub const PROBE_ARGUMENTS: [&str; 3] = ["--vanilla", "-e", "sessionInfo()"];

/// Arguments the source puts in front of a script: `--vanilla --quiet`.
pub const SCRIPT_ARGUMENTS: [&str; 2] = ["--vanilla", "--quiet"];

/// Largest number of command line arguments [`run_script`] passes on to a
/// script.
///
/// The source passes the caller's vector through unbounded; this preflights.
pub const MAX_SCRIPT_ARGUMENTS: usize = 1024;

/// Outcome of probing for the R interpreter.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RCheck {
    /// The source's return value: `true` when the interpreter ran
    /// `sessionInfo()` to a zero exit code.
    pub found: bool,
    /// Standard output followed by standard error, as the source merges them to
    /// imitate Qt's `MergedChannels`.
    pub output: String,
    /// The diagnosis the source writes to `OPENMS_LOG_ERROR`, empty on success.
    pub message: String,
    /// The resolved interpreter, when one was found. Native: the source keeps
    /// the resolution inside `bp::search_path` and never reports it.
    pub executable: Option<PathBuf>,
}

/// Outcome of running an R script.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ScriptOutcome {
    /// The source's return value.
    pub success: bool,
    /// The script that was run, absent when it could not be located.
    pub script: Option<PathBuf>,
    /// The script's standard output, the source's `stdout_content`.
    pub stdout: String,
    /// The script's standard error, the source's `stderr_content`.
    pub stderr: String,
    /// The diagnosis the source writes to `OPENMS_LOG_INFO` / `OPENMS_LOG_ERROR`,
    /// empty on success.
    pub message: String,
}

/// Rebuild a captured stream the way the source's drain loop does: read it line
/// by line with `std::getline` and append each line plus `'\n'`.
///
/// The only observable effect is that a final line without a newline gains one,
/// and that a completely empty stream stays empty. It is reproduced so that the
/// strings a caller inspects are the ones the C++ assembled.
///
/// # Examples
///
/// ```
/// use openms::system::r_wrapper::reassemble_lines;
///
/// assert_eq!(reassemble_lines("a\nb\n"), "a\nb\n");
/// assert_eq!(reassemble_lines("a\nb"), "a\nb\n");
/// assert_eq!(reassemble_lines(""), "");
/// ```
pub fn reassemble_lines(text: &str) -> String {
    if text.is_empty() {
        return String::new();
    }
    let trimmed = text.strip_suffix('\n').unwrap_or(text);
    let mut value = String::with_capacity(text.len() + 1);
    for line in trimmed.split('\n') {
        value.push_str(line);
        value.push('\n');
    }
    value
}

/// Look for an R script in the `share/OpenMS/SCRIPTS` folder.
///
/// The script filename can be absolute, in which case an existing file is
/// returned unchanged; otherwise it is searched in the OpenMS `SCRIPTS` path
/// and the full filename is returned.
///
/// # Arguments
///
/// * `context` — runtime locations; `get_openms_data_path` supplies the base
///   the source reads from `File::getOpenMSDataPath()`.
/// * `script_file` — name of the R script.
///
/// # Errors
///
/// Returns [`Error::Io`] with [`std::io::ErrorKind::NotFound`] when the file
/// cannot be found, which is this crate's counterpart of the source's
/// `Exception::FileNotFound`, and whatever
/// [`FileContext::get_openms_data_path`] reports when the shared data itself
/// cannot be located — the source turns that into the same `FileNotFound`,
/// losing the distinction between "no script" and "no OpenMS data".
///
/// The source's `verbose` parameter only chose whether to also print the
/// message; here the message is in the error.
pub fn find_script(context: &FileContext, script_file: &str) -> Result<PathBuf> {
    let scripts = context.get_openms_data_path()?.join(SCRIPTS_DIRECTORY);
    context.find(script_file, &[scripts])
}

fn probe_message(executable: &str, merged: &str) -> String {
    let arguments = PROBE_ARGUMENTS.join(" ");
    format!(
        "Error: '{executable}' executable returned with error (command: '{executable} {arguments}')\nOutput was:\n------>\n{merged}\n<------\nMake sure '{executable}' is installed properly.\n"
    )
}

fn missing_interpreter_message(executable: &str) -> String {
    format!(
        "Error: Could not find or run '{executable}' executable (FailedToStart).\nPlease install '{executable}', make sure it's in PATH and is flagged as executable.\n"
    )
}

/// Check for the presence of `Rscript`.
///
/// Runs `<executable> --vanilla -e sessionInfo()` and reports whether it
/// finished with a zero exit code. This is the "is R available?" signal callers
/// use to skip optional R-based steps, and it never fails loudly: an absent or
/// broken interpreter comes back as `found == false` with a `message`.
///
/// # Arguments
///
/// * `context` — runtime locations, whose `search_path` stands in for the
///   source's `bp::search_path`.
/// * `executable` — name of the R interpreter; pass [`DEFAULT_EXECUTABLE`] for
///   the source's default.
///
/// # Errors
///
/// Returns [`Error::Io`] when the filesystem search or the wait on the child
/// fails, and [`Error::InvalidValue`] when `executable` exceeds a ceiling.
///
/// # Notes
///
/// The source gives this probe no time budget, unlike the Java and Python
/// probes, and that is preserved: a wedged interpreter blocks the caller here
/// exactly as it does in C++. Set [`Invocation::timeout`] yourself through
/// [`capture`] if you need one.
///
/// The source's failure message names the literal string `Rscript` even when a
/// different interpreter was requested; this interpolates the requested name in
/// both places, so the diagnosis matches the command that actually ran.
pub fn find_r(context: &FileContext, executable: &str) -> Result<RCheck> {
    let Some(resolved) = context.find_executable(executable)? else {
        return Ok(RCheck {
            found: false,
            output: String::new(),
            message: missing_interpreter_message(executable),
            executable: None,
        });
    };
    let run = capture(
        &Invocation::new(&resolved)
            .with_arguments(PROBE_ARGUMENTS)
            .with_io_mode(IoMode::ReadWrite),
    )?;
    // Qt's MergedChannels put standard output first; the source imitates that
    // by concatenating the two drained buffers in this order.
    let merged = format!(
        "{}{}",
        reassemble_lines(&run.stdout),
        reassemble_lines(&run.stderr)
    );
    if run.report.state == ReturnState::FailedToStart {
        return Ok(RCheck {
            found: false,
            output: merged,
            message: missing_interpreter_message(executable),
            executable: Some(resolved),
        });
    }
    if !run.report.is_success() {
        let message = probe_message(executable, &merged);
        return Ok(RCheck {
            found: false,
            output: merged,
            message,
            executable: Some(resolved),
        });
    }
    Ok(RCheck {
        found: true,
        output: merged,
        message: String::new(),
        executable: Some(resolved),
    })
}

/// Run an R script with the given arguments on the command line.
///
/// The following checks are done before running the script:
///
/// 1. optionally, the interpreter is searched with [`find_r`] — set
///    `find_r_first`;
/// 2. `script_file` is searched with [`find_script`];
/// 3. the script is run as `<executable> --vanilla --quiet <script> <args…>`.
///
/// If any of those steps fails, the diagnosis is in
/// [`ScriptOutcome::message`] and `success` is `false`. `cmd_args` are passed
/// on the command line and should be read by the R script through R's
/// `commandArgs()`; usually they are input and output filenames.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `cmd_args` exceeds
/// [`MAX_SCRIPT_ARGUMENTS`] or an argument violates
/// [`Invocation::preflight`], and [`Error::Io`] when the wait on the child
/// fails. Neither a missing interpreter nor a missing script is an error: the
/// source swallows its own `FileNotFound` here and returns `false`, and so does
/// this.
///
/// # Notes
///
/// Arguments are passed as a vector and never interpolated into a shell string,
/// so a filename containing spaces or shell metacharacters arrives at the
/// script intact. That matches the source, which also uses `bp::args`.
///
/// As in [`find_r`], there is no time budget: a plotting script may legitimately
/// run for a long time, and adding one would change observable behaviour.
pub fn run_script(
    context: &FileContext,
    script_file: &str,
    cmd_args: &[String],
    executable: &str,
    find_r_first: bool,
) -> Result<ScriptOutcome> {
    if cmd_args.len() > MAX_SCRIPT_ARGUMENTS {
        return Err(Error::InvalidValue(format!(
            "R script argument count {} exceeds the ceiling of {MAX_SCRIPT_ARGUMENTS}",
            cmd_args.len()
        )));
    }
    if find_r_first {
        let check = find_r(context, executable)?;
        if !check.found {
            return Ok(ScriptOutcome {
                success: false,
                script: None,
                stdout: String::new(),
                stderr: String::new(),
                message: check.message,
            });
        }
    }
    let script = match find_script(context, script_file) {
        Ok(script) => script,
        Err(_) => {
            return Ok(ScriptOutcome {
                success: false,
                script: None,
                stdout: String::new(),
                stderr: String::new(),
                message: format!("\n\nCould not find R script '{script_file}'!\n\n"),
            });
        }
    };

    let mut arguments: Vec<String> = SCRIPT_ARGUMENTS.iter().map(|a| (*a).to_owned()).collect();
    arguments.push(script.display().to_string());
    arguments.extend(cmd_args.iter().cloned());

    let Some(interpreter) = context.find_executable(executable)? else {
        return Ok(ScriptOutcome {
            success: false,
            script: Some(script),
            stdout: String::new(),
            stderr: String::new(),
            message: format!("Error: Could not run '{executable}'. Is it installed and in PATH?\n"),
        });
    };
    let run = capture(
        &Invocation::new(interpreter)
            .with_arguments(arguments)
            .with_io_mode(IoMode::ReadWrite),
    )?;
    let stdout = reassemble_lines(&run.stdout);
    let stderr = reassemble_lines(&run.stderr);
    if run.report.state == ReturnState::FailedToStart {
        return Ok(ScriptOutcome {
            success: false,
            script: Some(script),
            stdout,
            stderr,
            message: format!("Error: Could not run '{executable}'. Is it installed and in PATH?\n"),
        });
    }
    if !run.report.is_success() {
        let message = format!(
            "\n--- ERROR MESSAGES ---\n{stderr}\n--- OTHER MESSAGES ---\n{stdout}\n\nScript failed. See above for an error description. \n"
        );
        return Ok(ScriptOutcome {
            success: false,
            script: Some(script),
            stdout,
            stderr,
            message,
        });
    }
    Ok(ScriptOutcome {
        success: true,
        script: Some(script),
        stdout,
        stderr,
        message: String::new(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn context() -> FileContext {
        let mut context = FileContext::new(
            std::env::temp_dir(),
            std::env::temp_dir(),
            std::env::temp_dir(),
        )
        .unwrap();
        context.search_path.clear();
        context.data_candidates.clear();
        context
    }

    #[test]
    fn a_bogus_interpreter_is_reported_cleanly() {
        let check = find_r(&context(), "this_is_not_a_real_R_interpreter_xyz").unwrap();
        assert!(!check.found);
        assert!(check.message.contains("FailedToStart"));
        assert_eq!(check.executable, None);
    }

    #[test]
    fn a_missing_script_is_a_not_found_error() {
        let error = find_script(&context(), "definitely_nonexistent_script_qwerty.R").unwrap_err();
        assert!(matches!(error, Error::Io(_)));
    }

    #[test]
    fn run_script_swallows_a_missing_interpreter_and_a_missing_script() {
        let context = context();
        let outcome = run_script(
            &context,
            "any.R",
            &[],
            "this_is_not_a_real_R_interpreter_xyz",
            true,
        )
        .unwrap();
        assert!(!outcome.success);
        assert_eq!(outcome.script, None);

        let outcome = run_script(
            &context,
            "definitely_nonexistent_script_qwerty.R",
            &[],
            DEFAULT_EXECUTABLE,
            false,
        )
        .unwrap();
        assert!(!outcome.success);
        assert!(outcome.message.contains("Could not find R script"));
    }

    #[test]
    fn too_many_script_arguments_are_refused_before_anything_runs() {
        let many = vec![String::from("x"); MAX_SCRIPT_ARGUMENTS + 1];
        assert!(run_script(&context(), "any.R", &many, DEFAULT_EXECUTABLE, false).is_err());
    }
}
