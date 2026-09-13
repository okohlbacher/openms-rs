// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Detecting a Java installation, from `SYSTEM/JavaInfo.h`.
//!
//! The source runs `java -version` under a 30 second budget and reports one
//! `bool`, writing its diagnosis to `OPENMS_LOG_ERROR`. This port returns the
//! diagnosis instead of printing it, so the caller decides where it goes; see
//! [`JavaCheck`](crate::system::java_info::JavaCheck) and
//! [`can_run`](crate::system::java_info::can_run).
//!
//! Similar modules exist for other external tools, for example
//! [`python_info`](crate::system::python_info) and
//! [`r_wrapper`](crate::system::r_wrapper).
//!
//! See `docs/JAVA_INFO_SUPPORT.md` for the full API mapping.

use crate::Result;
use crate::system::external_process::{self, Invocation, IoMode, PROBE_TIMEOUT, ReturnState};
use crate::system::file::FileContext;
use crate::system::path_utils::to_path;

/// The single argument the source passes: `java -version`.
pub const VERSION_ARGUMENT: &str = "-version";

/// Outcome of a Java probe.
///
/// The source returns a bare `bool` and sends everything else to
/// `OPENMS_LOG_ERROR`. This port has no global log stream, so the text the
/// source would have logged is returned in `message` and the caller prints it,
/// which is why the source's `verbose_on_error` parameter has no counterpart.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct JavaCheck {
    /// The source's return value: `true` when the Java executable was found and
    /// ran to a zero exit code.
    pub can_run: bool,
    /// The diagnosis, empty when `can_run` is true. Reproduces the source's
    /// `OPENMS_LOG_ERROR` text, including its leading `"Java-Check:\n"` line.
    pub message: String,
    /// Exit code of `java -version`, when it ran and exited normally.
    pub exit_code: Option<i32>,
    /// Whether the probe exceeded [`PROBE_TIMEOUT`] and Java was terminated.
    pub timed_out: bool,
    /// Everything `java -version` wrote, standard output then standard error.
    /// Java reports its banner on standard error, which is why the source pipes
    /// both. Native: the source discards this.
    pub output: String,
}

fn not_found_message(java_executable: &str, relative: bool) -> String {
    let mut text = String::from("Java-Check:\n");
    text.push_str(&format!("  Java not found at '{java_executable}'!\n"));
    text.push_str("  Make sure Java is installed and this location is correct.\n");
    if relative {
        let path = std::env::var_os("PATH").unwrap_or_default();
        text.push_str("  You might need to add the Java binary to your PATH variable\n");
        text.push_str("  or use an absolute path+filename pointing to Java.\n");
        text.push_str(&format!(
            "  The current SYSTEM PATH is: '{}'.\n\n",
            path.to_string_lossy()
        ));
        if cfg!(target_os = "macos") {
            text.push_str(
                "  On MacOSX, application bundles change the system PATH; Open your executable (e.g. KNIME/TOPPAS/TOPPView) from within the bundle (e.g. ./TOPPAS.app/Contents/MacOS/TOPPAS) to preserve the system PATH or use an absolute path to Java!\n",
            );
        }
        text.push('\n');
    } else {
        text.push_str("  You gave an absolute path to Java. Please check if it's correct.\n");
        text.push_str("  You can also try 'java' if your system path is correctly configured.\n\n");
    }
    text
}

/// Determine whether Java is installed and reachable.
///
/// The call fails if Java is not installed, or if a relative location is given
/// and Java is not on the search `PATH`.
///
/// # Arguments
///
/// * `context` — runtime locations, whose `search_path` stands in for the
///   source's `bp::search_path`. Native: the source reads the process `PATH`
///   directly, which is what makes it untestable without a Java installation.
/// * `java_executable` — path to the Java executable. May be absolute, relative
///   or just a filename.
///
/// Returns a [`JavaCheck`] whose `can_run` is the source's `bool`: `false` when
/// the executable cannot be called, `true` when it runs to a zero exit code.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`](crate::Error::InvalidValue) when
/// `java_executable` exceeds the path ceiling or contains a NUL byte, and
/// [`Error::Io`](crate::Error::Io) when the filesystem search or the wait on
/// the child fails. An absent or unusable Java is *not* an error: it is
/// `can_run == false` with a `message`, as the source returns `false`.
///
/// # Notes
///
/// The source builds its `PATH` hint with `path = getenv("PATH")`, which
/// constructs a `std::string` from a null pointer — undefined behaviour — when
/// `PATH` is unset, and caches the first value it reads in a function-local
/// `static` so a later change is never seen. This reads the variable on every
/// call and treats "unset" as empty.
///
/// # Examples
///
/// ```
/// use openms::system::file::FileContext;
/// use openms::system::java_info;
///
/// // An empty name can never resolve, which is what the class test asserts.
/// let context = FileContext::from_environment()?;
/// let check = java_info::can_run(&context, "")?;
/// assert!(!check.can_run);
/// assert!(check.message.starts_with("Java-Check:\n"));
/// # Ok::<(), openms::Error>(())
/// ```
pub fn can_run(context: &FileContext, java_executable: &str) -> Result<JavaCheck> {
    let relative = to_path(java_executable)?.is_relative();
    let Some(resolved) = context.find_executable(java_executable)? else {
        return Ok(JavaCheck {
            can_run: false,
            message: not_found_message(java_executable, relative),
            exit_code: None,
            timed_out: false,
            output: String::new(),
        });
    };

    let run = external_process::capture(
        &Invocation::new(resolved)
            .with_arguments([VERSION_ARGUMENT])
            .with_io_mode(IoMode::ReadWrite)
            .with_timeout(PROBE_TIMEOUT),
    )?;
    let output = format!("{}{}", run.stdout, run.stderr);

    if run.report.timed_out {
        return Ok(JavaCheck {
            can_run: false,
            message: format!(
                "Java-Check:\n  Java was found at '{java_executable}' but the process timed out (can happen on very busy systems).\n  Please free some resources or if you want to run the TOPP tool nevertheless set the TOPP tools 'force' flag in order to avoid this check.\n"
            ),
            exit_code: run.report.exit_code,
            timed_out: true,
            output,
        });
    }
    if run.report.state == ReturnState::FailedToStart {
        return Ok(JavaCheck {
            can_run: false,
            message: not_found_message(java_executable, relative),
            exit_code: None,
            timed_out: false,
            output,
        });
    }
    if !run.report.is_success() {
        return Ok(JavaCheck {
            can_run: false,
            message: format!(
                "Java-Check:\n  Error executing '{java_executable}'!\n  Java returned a non-zero exit code ({}).\n",
                run.report.exit_code.unwrap_or_default()
            ),
            exit_code: run.report.exit_code,
            timed_out: false,
            output,
        });
    }
    Ok(JavaCheck {
        can_run: true,
        message: String::new(),
        exit_code: run.report.exit_code,
        timed_out: false,
        output,
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
    fn an_empty_executable_can_never_run() {
        let check = can_run(&context(), "").unwrap();
        assert!(!check.can_run);
        assert!(check.message.contains("Java not found at ''!"));
    }

    #[test]
    fn a_relative_name_that_is_absent_mentions_the_path_variable() {
        let check = can_run(&context(), "no_such_java_@@").unwrap();
        assert!(!check.can_run);
        assert!(check.message.contains("The current SYSTEM PATH is:"));
    }

    #[test]
    fn an_absolute_name_that_is_absent_does_not_mention_the_path_variable() {
        let absent = std::env::temp_dir().join("no_such_java_@@");
        let check = can_run(&context(), absent.to_str().unwrap()).unwrap();
        assert!(!check.can_run);
        assert!(check.message.contains("You gave an absolute path to Java."));
        assert!(!check.message.contains("The current SYSTEM PATH is:"));
    }

    #[test]
    fn a_path_with_a_nul_byte_is_refused_rather_than_probed() {
        assert!(can_run(&context(), "java\0").is_err());
    }
}
