// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Tools that invoke an external executable: the native form of the
//! `OpenMS4-cli` package's `APPLICATIONS/TOPPExternalToolBase.h` and
//! `TOPPExternalToolBase.cpp` (cli `c19e494`).
//!
//! The source class derives from `TOPPBase` and adds one protected method,
//! `runExternalProcess_` (two overloads), "a thin adapter over
//! `ExternalProcess` that maps its `RETURNSTATE` to `TOPPBase::ExitCodes` and
//! logs stdout/stderr on failure (or when the debug level is >= 4)". It adds
//! behaviour and no state, so here it is two methods of the run context,
//! [`ToolContext::run_external_process`] and
//! [`ToolContext::run_external_process_output`], which every tool can call;
//! a tool that shells out implements [`Tool::run_io`](super::Tool::run_io) and
//! passes its output and error streams through. The source's deleted default
//! and copy constructors and its inherited constructor have no counterpart,
//! because a native tool is a type implementing [`Tool`](super::Tool), not an
//! object built from a name and a description; its defaulted destructor is
//! ordinary ownership.
//!
//! See `docs/TOPP_CLI_SUPPORT.md` (*Tool base classes*) for the member mapping
//! and the evidence.

use super::context::ToolContext;
use super::parameter::ExitCode;
use crate::system::external_process::{ExternalProcess, Invocation, IoMode, ReturnState};
use crate::{Error, Result};
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::io::Write;

/// The debug level from which every chunk of the child's output is echoed on
/// the output stream as it arrives (`TOPPExternalToolBase.cpp:36-37`).
pub const EXTERNAL_OUTPUT_DEBUG_LEVEL: i64 = 4;

/// What an external program run left behind: the exit code the tool returns
/// and everything the child wrote.
///
/// The source's second `runExternalProcess_` overload returns the exit code
/// and fills two out-parameters, `proc_stdout` and `proc_stderr`; the three
/// are one value here.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExternalProcessOutput {
    /// The run's outcome as a tool exit code; see
    /// [`ToolContext::run_external_process`].
    pub code: ExitCode,
    /// Everything that reached the standard-output callback, starting with the
    /// `Running: <executable> <arguments>` banner and ending, on success, with
    /// `Executed '<executable>' successfully!`.
    pub stdout: String,
    /// Everything that reached the standard-error callback, ending, on failure,
    /// with the `Process '<executable>' …` message.
    pub stderr: String,
}

/// The tool exit code for an external program's outcome
/// (`TOPPExternalToolBase.cpp:45-58`): success is
/// [`ExitCode::ExecutionOk`], a non-zero exit and a crash are both
/// [`ExitCode::ExternalProgramError`], and a program that could not be
/// started is [`ExitCode::ExternalProgramNotFound`].
///
/// The source's `default` arm throws `InternalToolError` ("Unknown return
/// state of external process."); [`ReturnState`] has exactly the four
/// states, so this match is exhaustive and that arm cannot be reached.
pub fn exit_code_of(state: ReturnState) -> ExitCode {
    match state {
        ReturnState::Success => ExitCode::ExecutionOk,
        ReturnState::NonzeroExit | ReturnState::Crash => ExitCode::ExternalProgramError,
        ReturnState::FailedToStart => ExitCode::ExternalProgramNotFound,
    }
}

impl ToolContext {
    /// Run an external program and map its outcome to the tool's exit code,
    /// printing what it wrote when it fails or when the debug level is at
    /// least [`EXTERNAL_OUTPUT_DEBUG_LEVEL`] (source `runExternalProcess_`,
    /// first overload).
    ///
    /// As [`run_external_process_output`](Self::run_external_process_output),
    /// which the source overload calls too, with the captured output dropped.
    ///
    /// # Errors
    ///
    /// As [`run_external_process_output`](Self::run_external_process_output).
    pub fn run_external_process(
        &self,
        out: &mut dyn Write,
        err: &mut dyn Write,
        executable: &str,
        arguments: &[String],
        workdir: &str,
        env: &BTreeMap<String, String>,
    ) -> Result<ExitCode> {
        self.run_external_process_output(out, err, executable, arguments, workdir, env)
            .map(|run| run.code)
    }

    /// Run an external program and return its exit code together with
    /// everything it wrote (source `runExternalProcess_`, second overload,
    /// `TOPPExternalToolBase.cpp:28-59`).
    ///
    /// The program runs as `ExternalProcess::run(executable, arguments,
    /// workdir, verbose = true, READ_WRITE, env)` runs it: `executable` is
    /// looked up on `PATH` unless it is a path, the child inherits this
    /// process's environment with `env` laid over it, an empty `workdir`
    /// means the current directory, and its standard output and standard
    /// error are read as they arrive. Because the run is verbose, the
    /// standard-output text starts with the `Running: …` banner and ends with
    /// `Executed '<executable>' successfully!` on success, and on failure the
    /// standard-error text ends with the process's `Process '<executable>' …`
    /// message (see [`ExternalProcess::run`]).
    ///
    /// * At debug level [`EXTERNAL_OUTPUT_DEBUG_LEVEL`] or above, every chunk
    ///   of either stream goes to the run's info log stream as it arrives
    ///   ([`write_info_stream`](Self::write_info_stream)), as the source's
    ///   callbacks send it to `OPENMS_LOG_INFO`: on `out`, with no line break
    ///   added, a line repeated within the stream's two-line cache counted
    ///   rather than printed, and nothing in the `-log` file.
    /// * Below that level, a failed run writes `Standard output: <stdout>`
    ///   and then `Standard error: <stderr>` through
    ///   [`write_log_error`](Self::write_log_error): on `err`, red on a
    ///   terminal, and in the `-log` file.
    ///
    /// The exit code is [`exit_code_of`] the run's [`ReturnState`].
    ///
    /// # Errors
    ///
    /// Returns what [`Invocation::preflight`] refuses — more arguments or
    /// command-line bytes than its ceilings, an interior NUL byte, an empty
    /// environment key or one containing `=` — before anything is started,
    /// [`Error::Io`] when waiting on the child fails, and [`Error::Io`] when
    /// `out` or `err` fails. The source checks none of the ceilings; an
    /// interior NUL would silently shorten the argument it hands to `exec`.
    /// A program that cannot be started is not an error: it is
    /// [`ExitCode::ExternalProgramNotFound`], as in the source.
    pub fn run_external_process_output(
        &self,
        out: &mut dyn Write,
        err: &mut dyn Write,
        executable: &str,
        arguments: &[String],
        workdir: &str,
        env: &BTreeMap<String, String>,
    ) -> Result<ExternalProcessOutput> {
        let mut invocation = Invocation::new(executable)
            .with_arguments(arguments.iter())
            .with_environment(env.clone())
            .with_io_mode(IoMode::ReadWrite);
        if !workdir.is_empty() {
            invocation = invocation.with_working_directory(workdir);
        }
        invocation.preflight()?;

        let echo = self.debug_level() >= EXTERNAL_OUTPUT_DEBUG_LEVEL;
        let stdout = RefCell::new(String::new());
        let stderr = RefCell::new(String::new());
        let console = RefCell::new(&mut *out);
        let failure: RefCell<Option<Error>> = RefCell::new(None);
        let forward = |target: &RefCell<String>, text: &str| {
            target.borrow_mut().push_str(text);
            if echo {
                let mut console = console.borrow_mut();
                if let Err(error) = self.write_info_stream(&mut **console, text) {
                    failure.borrow_mut().get_or_insert(error);
                }
            }
        };
        let report = {
            let mut process = ExternalProcess::with_callbacks(
                |text| forward(&stdout, text),
                |text| forward(&stderr, text),
            );
            process.run(&invocation, true)?
        };
        if let Some(error) = failure.into_inner() {
            return Err(error);
        }
        let stdout = stdout.into_inner();
        let stderr = stderr.into_inner();
        if !echo && report.state != ReturnState::Success {
            self.write_log_error(err, &format!("Standard output: {stdout}"))?;
            self.write_log_error(err, &format!("Standard error: {stderr}"))?;
        }
        Ok(ExternalProcessOutput {
            code: exit_code_of(report.state),
            stdout,
            stderr,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_return_state_maps_to_the_source_exit_code() {
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
}
