// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Detecting a Python installation, from `SYSTEM/PythonInfo.h`.
//!
//! Three probes: can the interpreter be run
//! ([`can_run`](crate::system::python_info::can_run)), is a package importable
//! ([`is_package_installed`](crate::system::python_info::is_package_installed)),
//! and what does `--version` say
//! ([`version`](crate::system::python_info::version)). Each runs the real
//! executable under the source's 30 second budget.
//!
//! The banner parsing is separated out as
//! [`version_from_output`](crate::system::python_info::version_from_output) so
//! it can be exercised on a recorded string rather than a live interpreter.
//!
//! Similar modules exist for other external tools, for example
//! [`java_info`](crate::system::java_info) and
//! [`r_wrapper`](crate::system::r_wrapper).
//!
//! See `docs/PYTHON_INFO_SUPPORT.md` for the full API mapping.

use crate::system::external_process::{Invocation, IoMode, PROBE_TIMEOUT, ReturnState, capture};
use crate::system::file::FileContext;
use crate::system::path_utils::to_path;
use crate::{Error, Result};
use std::path::PathBuf;

/// The argument the source passes to read a version banner: `python --version`.
pub const VERSION_ARGUMENT: &str = "--version";

/// The flag the source uses to test an import: `python -c "import <package>"`.
pub const COMMAND_ARGUMENT: &str = "-c";

/// Longest package name [`is_package_installed`] will accept.
pub const MAX_PACKAGE_NAME_BYTES: usize = 256;

/// Outcome of a Python probe.
///
/// The source's `canRun` takes `python_executable` by non-`const` reference and
/// rewrites it to the resolved absolute path, and fills an `error_msg`
/// out-parameter that also carries a purely informational "resolved to" note on
/// success. Both out-parameters are fields here, so nothing is hidden in an
/// argument the caller might not have expected to be modified.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PythonCheck {
    /// The source's return value: `true` when the interpreter ran to a zero
    /// exit code.
    pub can_run: bool,
    /// The resolved interpreter — the value the source writes back into its
    /// `python_executable` argument. Unchanged from the input when resolution
    /// failed.
    pub executable: PathBuf,
    /// The source's `error_msg`. Non-empty on success when the name resolved to
    /// a different path, because the source appends that note unconditionally.
    pub message: String,
    /// Exit code of `python --version`, when it ran and exited normally.
    pub exit_code: Option<i32>,
    /// Whether the probe exceeded [`PROBE_TIMEOUT`] and Python was terminated.
    pub timed_out: bool,
}

fn not_found_message(python_executable: &str, relative: bool) -> String {
    let mut text = format!("  Python not found at '{python_executable}'!\n");
    text.push_str("  Make sure Python is installed and this location is correct.\n");
    if relative {
        let path = std::env::var_os("PATH").unwrap_or_default();
        text.push_str("  You might need to add the Python binary to your PATH variable\n");
        text.push_str("  or use an absolute path+filename pointing to Python.\n");
        text.push_str(&format!(
            "  The current SYSTEM PATH is: '{}'.\n\n",
            path.to_string_lossy()
        ));
        if cfg!(target_os = "macos") {
            text.push_str(
                "  On MacOSX, application bundles change the system PATH; Open your executable (e.g. KNIME/TOPPAS/TOPPView) from within the bundle (e.g. ./TOPPAS.app/Contents/MacOS/TOPPAS) to preserve the system PATH or use an absolute path to Python!\n",
            );
        }
    }
    text
}

/// Determine whether Python is installed and executable.
///
/// The call fails if Python is not installed, or if a relative location is
/// given and Python is not on the search `PATH`. When Python is found, the
/// resolved absolute path is reported in [`PythonCheck::executable`]; when it
/// is not, the diagnosis is in [`PythonCheck::message`].
///
/// # Arguments
///
/// * `context` — runtime locations; its `search_path` and
///   [`FileContext::find_executable`] stand in for the source's
///   `File::findExecutable`.
/// * `python_executable` — path to the Python executable. May be absolute,
///   relative or just a filename.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `python_executable` exceeds the path
/// ceiling or contains a NUL byte, and [`Error::Io`] when the filesystem search
/// or the wait on the child fails. An absent or unusable interpreter is not an
/// error; it is `can_run == false`.
///
/// # Notes
///
/// After resolving the name to an absolute path, the source hands that absolute
/// path back to `bp::search_path`, which only ever walks the `PATH` entries. A
/// perfectly good absolute interpreter therefore fails the probe whenever
/// `PATH` is empty or unset. This port runs the resolved path directly, so the
/// answer no longer depends on `PATH` once resolution has succeeded.
///
/// The source's `PATH` hint is built with `getenv("PATH")` assigned into a
/// `std::string`, which is undefined behaviour when `PATH` is unset, and cached
/// in a function-local `static`. This reads the variable on every call.
pub fn can_run(context: &FileContext, python_executable: &str) -> Result<PythonCheck> {
    let relative = to_path(python_executable)?.is_relative();
    let Some(resolved) = context.find_executable(python_executable)? else {
        return Ok(PythonCheck {
            can_run: false,
            executable: PathBuf::from(python_executable),
            message: not_found_message(python_executable, relative),
            exit_code: None,
            timed_out: false,
        });
    };

    let mut message = String::new();
    if resolved.as_os_str() != python_executable {
        message.push_str(&format!(
            "Python executable ('{python_executable}') resolved to '{}'\n",
            resolved.display()
        ));
    }
    let shown = resolved.display().to_string();
    let run = capture(
        &Invocation::new(&resolved)
            .with_arguments([VERSION_ARGUMENT])
            .with_io_mode(IoMode::ReadWrite)
            .with_timeout(PROBE_TIMEOUT),
    )?;

    if run.report.timed_out {
        message.push_str(&format!(
            "  Python was found at '{shown}' but the process timed out (can happen on very busy systems).\n  Please free some resources or if you want to run the TOPP tool nevertheless set the TOPP tools 'force' flag in order to avoid this check.\n"
        ));
        return Ok(PythonCheck {
            can_run: false,
            executable: resolved,
            message,
            exit_code: run.report.exit_code,
            timed_out: true,
        });
    }
    if run.report.state == ReturnState::FailedToStart {
        message.push_str(&format!(
            "  Python found at '{shown}' but failed to run!\n  Make sure you have the rights to execute this binary file.\n"
        ));
        return Ok(PythonCheck {
            can_run: false,
            executable: resolved,
            message,
            exit_code: None,
            timed_out: false,
        });
    }
    Ok(PythonCheck {
        can_run: run.report.is_success(),
        executable: resolved,
        message,
        exit_code: run.report.exit_code,
        timed_out: false,
    })
}

/// Whether `name` is a Python module path this port is willing to interpolate
/// into `python -c "import <name>"`.
///
/// Accepts a dotted sequence of identifiers — each component starting with an
/// ASCII letter or underscore and continuing with ASCII letters, digits or
/// underscores — of at most [`MAX_PACKAGE_NAME_BYTES`] bytes.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] for anything else. This check has no
/// counterpart in the source, which concatenates the caller's string into a
/// Python program and executes it: `isPackageInstalled(py, "os, os.system('…')")`
/// runs that command. Non-ASCII identifiers, which Python 3 permits, are
/// refused here; the narrower rule is deliberate, since the probe exists to
/// answer a yes/no question about ordinary package names.
pub fn validate_package_name(name: &str) -> Result<()> {
    if name.is_empty() || name.len() > MAX_PACKAGE_NAME_BYTES {
        return Err(Error::InvalidValue(format!(
            "Python package name must be 1..={MAX_PACKAGE_NAME_BYTES} bytes"
        )));
    }
    for component in name.split('.') {
        let mut characters = component.chars();
        let Some(first) = characters.next() else {
            return Err(Error::InvalidValue(
                "Python package name has an empty component".into(),
            ));
        };
        if !(first.is_ascii_alphabetic() || first == '_') {
            return Err(Error::InvalidValue(format!(
                "Python package component '{component}' must start with a letter or underscore"
            )));
        }
        if !characters.all(|c| c.is_ascii_alphanumeric() || c == '_') {
            return Err(Error::InvalidValue(format!(
                "Python package component '{component}' must be alphanumeric or underscore"
            )));
        }
    }
    Ok(())
}

/// Determine whether the Python at `python_executable` already has the package
/// `package_name` installed, by running `python -c "import <package_name>"`.
///
/// If Python cannot be found the function returns `false`, so make sure
/// [`can_run`] succeeds before calling this. Package names are case sensitive.
///
/// A `package_name` that is not a syntactically valid module path returns
/// `Ok(false)` **without starting an interpreter**: such a name cannot name an
/// installed package, so `false` is the right answer and is what the source
/// reports too, having executed the string to find out. Call
/// [`validate_package_name`] first when the distinction matters.
///
/// # Errors
///
/// Returns [`Error::Io`] when the wait on the child fails, and whatever
/// [`Invocation::preflight`] rejects.
///
/// # Notes
///
/// The source pipes both of the child's streams to the null device here, so
/// this uses [`IoMode::NoIo`].
pub fn is_package_installed(
    context: &FileContext,
    python_executable: &str,
    package_name: &str,
) -> Result<bool> {
    if validate_package_name(package_name).is_err() {
        return Ok(false);
    }
    let Some(resolved) = context.find_executable(python_executable)? else {
        return Ok(false);
    };
    let run = capture(
        &Invocation::new(resolved)
            .with_arguments([
                COMMAND_ARGUMENT.to_owned(),
                format!("import {package_name}"),
            ])
            .with_io_mode(IoMode::NoIo)
            .with_timeout(PROBE_TIMEOUT),
    )?;
    Ok(run.report.is_success())
}

/// Reassemble a version banner the way the source does.
///
/// The source reads `python --version` line by line with `std::getline`,
/// appending each line to one string **without a separator**, first from
/// standard output and then from standard error — some interpreters report the
/// version on standard error — and finally trims spaces, tabs, newlines and
/// carriage returns from both ends.
///
/// The missing separator is the source's behaviour, not an oversight of this
/// port: a two-line banner comes back with its lines run together. It is
/// reproduced because a caller matching on the string would otherwise see a
/// different value than the C++ produced.
///
/// # Examples
///
/// ```
/// use openms::system::python_info::version_from_output;
///
/// assert_eq!(version_from_output("Python 3.11.4\n", ""), "Python 3.11.4");
/// // Python 2 printed its banner on standard error.
/// assert_eq!(version_from_output("", "Python 2.7.18\n"), "Python 2.7.18");
/// // Two lines are concatenated with nothing between them.
/// assert_eq!(version_from_output("first\nsecond\n", ""), "firstsecond");
/// ```
pub fn version_from_output(stdout: &str, stderr: &str) -> String {
    let mut value = String::new();
    for text in [stdout, stderr] {
        for line in getline_sequence(text) {
            value.push_str(line);
        }
    }
    value
        .trim_matches(|c| c == ' ' || c == '\t' || c == '\n' || c == '\r')
        .to_owned()
}

/// The lines `std::getline` would yield from `text`: split on `'\n'`, with the
/// delimiter dropped and no empty final element for a trailing newline.
fn getline_sequence(text: &str) -> impl Iterator<Item = &str> {
    let empty = text.is_empty();
    let trimmed = text.strip_suffix('\n').unwrap_or(text);
    trimmed.split('\n').filter(move |_| !empty)
}

/// Determine the version of the Python at `python_executable` by calling
/// `--version`.
///
/// Returns the empty string when Python cannot be found, when it could not be
/// started, when it timed out, or when it exited non-zero — so make sure
/// [`can_run`] succeeds before calling this. The banner itself is assembled by
/// [`version_from_output`].
///
/// # Errors
///
/// Returns [`Error::Io`] when the filesystem search or the wait on the child
/// fails, and [`Error::InvalidValue`] when `python_executable` exceeds a
/// ceiling.
pub fn version(context: &FileContext, python_executable: &str) -> Result<String> {
    let Some(resolved) = context.find_executable(python_executable)? else {
        return Ok(String::new());
    };
    let run = capture(
        &Invocation::new(resolved)
            .with_arguments([VERSION_ARGUMENT])
            .with_io_mode(IoMode::ReadWrite)
            .with_timeout(PROBE_TIMEOUT),
    )?;
    if !run.report.is_success() {
        return Ok(String::new());
    }
    Ok(version_from_output(&run.stdout, &run.stderr))
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
    fn an_absent_interpreter_reports_python_not_found() {
        let check = can_run(&context(), "does_not_exist_@@").unwrap();
        assert!(!check.can_run);
        assert!(check.message.contains("Python not found at"));
        assert_eq!(check.executable, PathBuf::from("does_not_exist_@@"));
    }

    #[test]
    fn a_banner_on_either_stream_is_trimmed() {
        assert_eq!(
            version_from_output("  Python 3.9.7  \n", ""),
            "Python 3.9.7"
        );
        assert_eq!(version_from_output("Python 3.9.7\r\n", ""), "Python 3.9.7");
        assert_eq!(version_from_output("", ""), "");
        assert_eq!(version_from_output("\n\n", ""), "");
    }

    #[test]
    fn a_package_name_that_could_inject_code_is_refused() {
        assert!(validate_package_name("math").is_ok());
        assert!(validate_package_name("os.path").is_ok());
        assert!(validate_package_name("_private1").is_ok());
        assert!(validate_package_name("os; import sys").is_err());
        assert!(validate_package_name("veryWeirdPackage___@@__@").is_err());
        assert!(validate_package_name("").is_err());
        assert!(validate_package_name("a..b").is_err());
        assert!(validate_package_name(&"a".repeat(MAX_PACKAGE_NAME_BYTES + 1)).is_err());
    }

    #[test]
    fn an_unusable_package_name_answers_false_without_an_interpreter() {
        assert!(
            !is_package_installed(&context(), "does_not_exist_@@", "veryWeirdPackage___@@__@")
                .unwrap()
        );
    }

    #[test]
    fn the_version_of_an_absent_interpreter_is_empty() {
        assert_eq!(version(&context(), "does_not_exist_@@").unwrap(), "");
    }
}
