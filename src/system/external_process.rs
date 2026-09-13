// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Starting an external program and forwarding its output, from
//! `SYSTEM/ExternalProcess.h`.
//!
//! The source wraps `boost::process` and `boost::asio`: it resolves the program
//! on `PATH`, inherits the caller's environment plus a supplied overlay, streams
//! the child's two output pipes into callbacks while polling a 50 ms loop, and
//! finally reports one of four outcomes. `std::process::Command` covers all of
//! that, so this module takes no dependency.
//!
//! [`ExternalProcess`](crate::system::external_process::ExternalProcess) is the
//! callback-driven form and
//! [`capture`](crate::system::external_process::capture) the accumulating one
//! used by the interpreter probes in
//! [`java_info`](crate::system::java_info),
//! [`python_info`](crate::system::python_info) and
//! [`r_wrapper`](crate::system::r_wrapper).
//!
//! See `docs/EXTERNAL_PROCESS_SUPPORT.md` for the full API mapping.

use crate::{Error, Result};
use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, ExitStatus, Stdio};
use std::sync::mpsc::{self, Receiver, RecvTimeoutError, Sender};
use std::thread::{self, JoinHandle};
use std::time::{Duration, Instant};

/// Largest number of arguments a single invocation may carry.
///
/// The source passes `args` straight to `boost::process` with no ceiling; this
/// port preflights, so an oversized request fails before anything is spawned.
pub const MAX_ARGUMENTS: usize = 4096;

/// Largest combined size of the program name and every argument, in bytes.
///
/// Comfortably below the usual `ARG_MAX`, so the ceiling reports a clear error
/// instead of leaving the kernel to refuse the `exec`.
pub const MAX_COMMAND_BYTES: usize = 1024 * 1024;

/// Largest number of environment overrides a single invocation may carry.
pub const MAX_ENVIRONMENT_ENTRIES: usize = 4096;

/// Largest amount of output [`capture`] accumulates per stream, in bytes.
///
/// [`ExternalProcess::run`] has no such ceiling because it never accumulates:
/// every chunk goes straight to a callback.
pub const MAX_CAPTURED_BYTES: usize = 64 * 1024 * 1024;

/// Poll cadence of the run loop, matching the source's
/// `io_ctx.run_for(milliseconds(50))` and its `sleep_for(milliseconds(50))`.
pub const POLL_INTERVAL: Duration = Duration::from_millis(50);

/// Read buffer size, matching the source's `std::array<char, 4096>`.
pub const READ_BUFFER_BYTES: usize = 4096;

/// Budget the source's interpreter probes give a child before terminating it,
/// from `Internal::waitForProcess(child, std::chrono::seconds(30))`.
///
/// `ExternalProcess` itself imposes no budget, and neither does this port's
/// [`Invocation::timeout`] unless a caller sets one.
pub const PROBE_TIMEOUT: Duration = Duration::from_secs(30);

/// Which of a child's two output streams a chunk of bytes came from.
///
/// Native: the source keeps the streams apart by having two separate callbacks
/// rather than by tagging the data.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Stream {
    /// The child's standard output.
    Standard,
    /// The child's standard error.
    Error,
}

/// Result of calling an external executable, the source's `RETURNSTATE`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ReturnState {
    /// Everything went smoothly (exit code 0).
    Success,
    /// Finished, but returned with an exit code other than 0.
    NonzeroExit,
    /// Ran, but crashed (segfault and similar).
    Crash,
    /// Executable not found, or not enough access rights for the user.
    FailedToStart,
}

impl ReturnState {
    /// True only for [`ReturnState::Success`], the source's `RETURNSTATE::SUCCESS`.
    pub fn is_success(self) -> bool {
        matches!(self, Self::Success)
    }
}

/// Open mode for the process, the source's `IO_MODE`.
///
/// # Notes
///
/// The source distinguishes four modes but acts on one question only — may the
/// parent read? `NO_IO` and `WRITE_ONLY` both send the child's output to the
/// null device, `READ_ONLY` and `READ_WRITE` both capture it, and neither
/// branch ever wires the child's standard *input*, so the "write" half of the
/// name has no effect at all. That is reproduced here rather than tidied away,
/// and [`IoMode::reads`] is the predicate the implementation actually asks.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum IoMode {
    /// No read nor write access; output is discarded.
    NoIo,
    /// Read access; output is forwarded.
    ReadOnly,
    /// Write access. Identical to [`IoMode::NoIo`] in the source and here.
    WriteOnly,
    /// Read and write access. The source's default argument, and this port's
    /// [`Default`].
    #[default]
    ReadWrite,
}

impl IoMode {
    /// Whether the parent captures the child's output, the source's `can_read`.
    pub fn reads(self) -> bool {
        matches!(self, Self::ReadOnly | Self::ReadWrite)
    }
}

/// One external program call: what to run, where, and with which environment.
///
/// The source spreads these across `run`'s eight parameters. Collecting them
/// into a value keeps the two `run` overloads from becoming four Rust methods
/// and gives [`Invocation::preflight`] one place to check its ceilings.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Invocation {
    /// The program to call. May contain spaces; a bare name is looked up on
    /// `PATH`, replacing the source's explicit `bp::search_path` call.
    pub program: PathBuf,
    /// Extra arguments, passed as a vector and never as one shell string.
    pub arguments: Vec<String>,
    /// Directory to execute in. `None` means the source's `"."`, which is the
    /// current working directory.
    pub working_directory: Option<PathBuf>,
    /// Environment variables added on top of the inherited environment, as the
    /// source's `env` map. A [`BTreeMap`] keeps the application order
    /// deterministic where the source's `std::map` does the same.
    pub environment: BTreeMap<String, String>,
    /// Open mode for the process.
    pub io_mode: IoMode,
    /// Optional budget after which the child is killed. `None` — the default,
    /// and what [`ExternalProcess`] uses — reproduces the source, which waits
    /// indefinitely; the interpreter probes set [`PROBE_TIMEOUT`].
    pub timeout: Option<Duration>,
}

impl Invocation {
    /// An invocation of `program` with no arguments, in the current directory,
    /// inheriting the environment unchanged, reading output, without a budget.
    pub fn new(program: impl Into<PathBuf>) -> Self {
        Self {
            program: program.into(),
            ..Self::default()
        }
    }

    /// Replace the argument list.
    pub fn with_arguments<I, S>(mut self, arguments: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.arguments = arguments.into_iter().map(Into::into).collect();
        self
    }

    /// Run in `directory` rather than the current working directory.
    pub fn with_working_directory(mut self, directory: impl Into<PathBuf>) -> Self {
        self.working_directory = Some(directory.into());
        self
    }

    /// Add environment entries on top of the inherited environment.
    pub fn with_environment(mut self, environment: BTreeMap<String, String>) -> Self {
        self.environment = environment;
        self
    }

    /// Select the open mode; see [`IoMode`].
    pub fn with_io_mode(mut self, io_mode: IoMode) -> Self {
        self.io_mode = io_mode;
        self
    }

    /// Terminate the child if it has not finished within `timeout`.
    pub fn with_timeout(mut self, timeout: Duration) -> Self {
        self.timeout = Some(timeout);
        self
    }

    /// The verbose banner the source writes to the stdout callback before
    /// starting: `"Running: "`, the program, then each argument separated by a
    /// single space, then a newline.
    ///
    /// This is a display string only. It is never handed to a shell, and the
    /// arguments it shows unquoted are passed to the child as a vector.
    pub fn command_line(&self) -> String {
        let mut text = format!("Running: {}", self.program.display());
        for argument in &self.arguments {
            text.push(' ');
            text.push_str(argument);
        }
        text.push('\n');
        text
    }

    /// Check every ceiling before anything is spawned or allocated.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the argument count exceeds
    /// [`MAX_ARGUMENTS`], the environment overlay exceeds
    /// [`MAX_ENVIRONMENT_ENTRIES`], the program name and arguments together
    /// exceed [`MAX_COMMAND_BYTES`], an environment key is empty or contains
    /// `=`, or any of those strings contains an interior NUL byte. The source
    /// checks none of this: an interior NUL silently truncates the C string it
    /// hands to `exec`, so the child would receive a different argument than
    /// the caller wrote.
    pub fn preflight(&self) -> Result<()> {
        if self.arguments.len() > MAX_ARGUMENTS {
            return Err(Error::InvalidValue(format!(
                "process argument count {} exceeds the ceiling of {MAX_ARGUMENTS}",
                self.arguments.len()
            )));
        }
        if self.environment.len() > MAX_ENVIRONMENT_ENTRIES {
            return Err(Error::InvalidValue(format!(
                "process environment entry count {} exceeds the ceiling of {MAX_ENVIRONMENT_ENTRIES}",
                self.environment.len()
            )));
        }
        let mut total = self.program.as_os_str().len();
        reject_nul("program name", self.program.as_os_str().as_encoded_bytes())?;
        for argument in &self.arguments {
            reject_nul("process argument", argument.as_bytes())?;
            total = total
                .checked_add(argument.len())
                .and_then(|value| value.checked_add(1))
                .ok_or_else(|| Error::InvalidValue("process command line size overflow".into()))?;
        }
        if total > MAX_COMMAND_BYTES {
            return Err(Error::InvalidValue(format!(
                "process command line of {total} bytes exceeds the ceiling of {MAX_COMMAND_BYTES}"
            )));
        }
        for (key, value) in &self.environment {
            if key.is_empty() || key.contains('=') {
                return Err(Error::InvalidValue(format!(
                    "environment key '{key}' must be non-empty and free of '='"
                )));
            }
            reject_nul("environment key", key.as_bytes())?;
            reject_nul("environment value", value.as_bytes())?;
        }
        Ok(())
    }
}

fn reject_nul(what: &str, bytes: &[u8]) -> Result<()> {
    if bytes.contains(&0) {
        return Err(Error::InvalidValue(format!(
            "{what} must not contain an interior NUL byte"
        )));
    }
    Ok(())
}

/// What happened to one child process.
///
/// The source returns `RETURNSTATE` and fills an `error_msg` out-parameter; the
/// two are one value here. `exit_code`, `terminating_signal` and `timed_out`
/// are native additions, so that the detail behind [`ReturnState::Crash`] is
/// available rather than only summarised in prose.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RunReport {
    /// The source's `RETURNSTATE`.
    pub state: ReturnState,
    /// The source's `error_msg`: empty on success, and a message to display to
    /// the user otherwise.
    pub error_message: String,
    /// Exit code, when the child exited normally. `None` when it was killed by
    /// a signal or never started.
    pub exit_code: Option<i32>,
    /// Unix signal that terminated the child, when one did.
    pub terminating_signal: Option<i32>,
    /// Whether [`Invocation::timeout`] elapsed and the child was killed.
    pub timed_out: bool,
}

impl RunReport {
    /// True only when the child ran and exited with code 0.
    pub fn is_success(&self) -> bool {
        self.state.is_success()
    }
}

/// A finished run together with everything it wrote.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CapturedRun {
    /// Outcome and diagnostics, exactly as [`ExternalProcess::run`] reports them.
    pub report: RunReport,
    /// Everything the child wrote to standard output, empty unless
    /// [`IoMode::reads`] held.
    pub stdout: String,
    /// Everything the child wrote to standard error, empty unless
    /// [`IoMode::reads`] held.
    pub stderr: String,
}

#[derive(Clone, Copy, Debug)]
struct Outcome {
    state: ReturnState,
    exit_code: Option<i32>,
    terminating_signal: Option<i32>,
    timed_out: bool,
}

impl Outcome {
    fn failed_to_start() -> Self {
        Self {
            state: ReturnState::FailedToStart,
            exit_code: None,
            terminating_signal: None,
            timed_out: false,
        }
    }
}

#[cfg(unix)]
fn status_parts(status: ExitStatus) -> (Option<i32>, Option<i32>) {
    use std::os::unix::process::ExitStatusExt;
    (status.code(), status.signal())
}

#[cfg(not(unix))]
fn status_parts(status: ExitStatus) -> (Option<i32>, Option<i32>) {
    (status.code(), None)
}

/// Whether a finished child died abnormally.
///
/// Unix: the source asks `WIFSIGNALED(native_exit_code)`, which is exactly
/// `ExitStatus::signal().is_some()`. Windows: the source asks
/// `exit_code < 0 || static_cast<unsigned int>(exit_code) > 0x80000000u`. The
/// second half of that disjunction can never fire on its own — for a
/// non-negative `i32` the unsigned value is at most `0x7fffffff`, and for a
/// negative one the first half already holds — so the whole condition reduces
/// to `exit_code < 0`, which is what this writes.
fn crashed(exit_code: Option<i32>, terminating_signal: Option<i32>) -> bool {
    if terminating_signal.is_some() {
        return true;
    }
    cfg!(windows) && exit_code.is_some_and(|code| code < 0)
}

fn pump<R: Read + Send + 'static>(
    stream: Stream,
    mut reader: R,
    sender: Sender<(Stream, Vec<u8>)>,
) -> std::io::Result<JoinHandle<()>> {
    // Builder rather than thread::spawn: the latter panics when the operating
    // system refuses a thread, and this port does not panic on resource limits.
    thread::Builder::new().spawn(move || {
        let mut buffer = [0_u8; READ_BUFFER_BYTES];
        loop {
            match reader.read(&mut buffer) {
                Ok(0) => break,
                Ok(count) => match buffer.get(..count) {
                    Some(chunk) => {
                        if sender.send((stream, chunk.to_vec())).is_err() {
                            break;
                        }
                    }
                    None => break,
                },
                Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
                Err(_) => break,
            }
        }
    })
}

fn spawn(invocation: &Invocation) -> std::io::Result<Child> {
    let mut command = Command::new(&invocation.program);
    command.args(&invocation.arguments);
    command.current_dir(
        invocation
            .working_directory
            .as_deref()
            .unwrap_or_else(|| Path::new(".")),
    );
    for (key, value) in &invocation.environment {
        command.env(key, value);
    }
    if invocation.io_mode.reads() {
        command.stdout(Stdio::piped()).stderr(Stdio::piped());
    } else {
        command.stdout(Stdio::null()).stderr(Stdio::null());
    }
    command.spawn()
}

fn poll_loop(
    child: &mut Child,
    receiver: &Receiver<(Stream, Vec<u8>)>,
    deadline: Option<Instant>,
    idle: &mut Option<&mut dyn FnMut()>,
    sink: &mut dyn FnMut(Stream, &[u8]) -> Result<()>,
) -> Result<bool> {
    let mut connected = true;
    loop {
        if connected {
            match receiver.recv_timeout(POLL_INTERVAL) {
                Ok((stream, bytes)) => sink(stream, &bytes)?,
                Err(RecvTimeoutError::Timeout) => {}
                Err(RecvTimeoutError::Disconnected) => connected = false,
            }
        } else {
            thread::sleep(POLL_INTERVAL);
        }
        if let Some(callback) = idle.as_deref_mut() {
            callback();
        }
        let finished = child.try_wait()?.is_some();
        if finished && !connected {
            return Ok(false);
        }
        if let Some(limit) = deadline {
            if Instant::now() >= limit && !finished {
                return Ok(true);
            }
        }
    }
}

fn execute(
    invocation: &Invocation,
    mut idle: Option<&mut dyn FnMut()>,
    sink: &mut dyn FnMut(Stream, &[u8]) -> Result<()>,
) -> Result<Outcome> {
    invocation.preflight()?;
    if invocation.program.as_os_str().is_empty() {
        // The source reaches the same answer by a longer route: search_path("")
        // yields an empty path, bp::child then throws process_error, and the
        // catch block reports FAILED_TO_START.
        return Ok(Outcome::failed_to_start());
    }
    let mut child = match spawn(invocation) {
        Ok(child) => child,
        Err(_) => return Ok(Outcome::failed_to_start()),
    };

    let (sender, receiver) = mpsc::channel::<(Stream, Vec<u8>)>();
    let mut pumps = Vec::new();
    let mut refused = None;
    if let Some(pipe) = child.stdout.take() {
        match pump(Stream::Standard, pipe, sender.clone()) {
            Ok(handle) => pumps.push(handle),
            Err(error) => refused = Some(error),
        }
    }
    if refused.is_none() {
        if let Some(pipe) = child.stderr.take() {
            match pump(Stream::Error, pipe, sender.clone()) {
                Ok(handle) => pumps.push(handle),
                Err(error) => refused = Some(error),
            }
        }
    }
    drop(sender);
    if let Some(error) = refused {
        let _ = child.kill();
        let _ = child.wait();
        join_all(pumps);
        return Err(error.into());
    }

    // `Instant + Duration` panics on overflow; a budget that cannot be
    // represented is treated as no budget at all, which is the source's
    // behaviour for `ExternalProcess`.
    let deadline = invocation
        .timeout
        .and_then(|budget| Instant::now().checked_add(budget));
    let timed_out = match poll_loop(&mut child, &receiver, deadline, &mut idle, sink) {
        Ok(value) => value,
        Err(error) => {
            let _ = child.kill();
            let _ = child.wait();
            join_all(pumps);
            return Err(error);
        }
    };
    if timed_out {
        let _ = child.kill();
    }
    let status = child.wait()?;
    join_all(pumps);

    let (exit_code, terminating_signal) = status_parts(status);
    let state = if crashed(exit_code, terminating_signal) {
        ReturnState::Crash
    } else if exit_code.is_some_and(|code| code != 0) {
        ReturnState::NonzeroExit
    } else {
        ReturnState::Success
    };
    Ok(Outcome {
        state,
        exit_code,
        terminating_signal,
        timed_out,
    })
}

fn join_all(pumps: Vec<JoinHandle<()>>) {
    for handle in pumps {
        let _ = handle.join();
    }
}

fn report(invocation: &Invocation, outcome: Outcome) -> RunReport {
    let program = invocation.program.display();
    let error_message = if outcome.timed_out {
        format!("Process '{program}' did not finish in time and was terminated.")
    } else {
        match outcome.state {
            ReturnState::Success => String::new(),
            ReturnState::FailedToStart => {
                format!("Process '{program}' failed to start. Does it exist? Is it executable?")
            }
            ReturnState::Crash => {
                format!("Process '{program}' crashed hard (segfault-like). Please check the log.")
            }
            ReturnState::NonzeroExit => format!(
                "Process '{program}' did not finish successfully (exit code: {}). Please check the log.",
                outcome.exit_code.unwrap_or_default()
            ),
        }
    };
    RunReport {
        state: outcome.state,
        error_message,
        exit_code: outcome.exit_code,
        terminating_signal: outcome.terminating_signal,
        timed_out: outcome.timed_out,
    }
}

/// Run a program and collect everything it writes.
///
/// This is the shape the source's interpreter probes use, where `boost::process`
/// pipes into an `ipstream` that is drained after the wait. The probes here call
/// it with [`PROBE_TIMEOUT`].
///
/// # Errors
///
/// Returns whatever [`Invocation::preflight`] rejects, [`Error::InvalidValue`]
/// when either stream exceeds [`MAX_CAPTURED_BYTES`] — at which point the child
/// is killed and reaped — and [`Error::Io`] when waiting on the child fails. A
/// program that cannot be started is *not* an error: it is
/// [`ReturnState::FailedToStart`], as in the source.
///
/// # Notes
///
/// Output that is not valid UTF-8 is replaced lossily rather than rejected,
/// because the source stores raw bytes in a `std::string` and a version banner
/// is not worth failing a probe over.
pub fn capture(invocation: &Invocation) -> Result<CapturedRun> {
    let mut stdout: Vec<u8> = Vec::new();
    let mut stderr: Vec<u8> = Vec::new();
    let outcome = {
        let mut sink = |stream: Stream, bytes: &[u8]| -> Result<()> {
            let target = match stream {
                Stream::Standard => &mut stdout,
                Stream::Error => &mut stderr,
            };
            let total = target
                .len()
                .checked_add(bytes.len())
                .ok_or_else(|| Error::InvalidValue("captured output size overflow".into()))?;
            if total > MAX_CAPTURED_BYTES {
                return Err(Error::InvalidValue(format!(
                    "captured process output of {total} bytes exceeds the ceiling of {MAX_CAPTURED_BYTES}"
                )));
            }
            target.extend_from_slice(bytes);
            Ok(())
        };
        execute(invocation, None, &mut sink)?
    };
    Ok(CapturedRun {
        report: report(invocation, outcome),
        stdout: String::from_utf8_lossy(&stdout).into_owned(),
        stderr: String::from_utf8_lossy(&stderr).into_owned(),
    })
}

/// A wrapper around [`std::process::Command`] to start an external program and
/// forward its output, the source's `ExternalProcess`.
///
/// Provide callbacks for standard output and standard error through
/// [`ExternalProcess::with_callbacks`], or set them later with
/// [`ExternalProcess::set_callbacks`].
///
/// Running an external program blocks the caller, so do not use this on a main
/// GUI thread unless you have some other means of telling the user that no
/// interaction is possible at the moment.
///
/// The source additionally points at `ExternalProcessMBox`, a Qt wrapper that
/// raises message boxes on failure; it lives outside the scientific SDK and has
/// no counterpart here.
///
/// # Examples
///
/// ```
/// use openms::system::external_process::{ExternalProcess, Invocation, ReturnState};
///
/// let mut collected = String::new();
/// let mut process = ExternalProcess::with_callbacks(|text| collected.push_str(text), |_| {});
/// let report = process.run(&Invocation::new("no_such_program_@@"), false)?;
/// assert_eq!(report.state, ReturnState::FailedToStart);
/// assert!(report.error_message.contains("failed to start"));
/// # Ok::<(), openms::Error>(())
/// ```
pub struct ExternalProcess<'a> {
    standard_callback: Box<dyn FnMut(&str) + 'a>,
    error_callback: Box<dyn FnMut(&str) + 'a>,
}

impl std::fmt::Debug for ExternalProcess<'_> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.debug_struct("ExternalProcess").finish()
    }
}

impl Default for ExternalProcess<'_> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> ExternalProcess<'a> {
    /// Default construction; the callbacks for standard output and standard
    /// error do nothing, as the source's default constructor installs two empty
    /// lambdas.
    pub fn new() -> Self {
        Self {
            standard_callback: Box::new(|_| {}),
            error_callback: Box::new(|_| {}),
        }
    }

    /// Set the callbacks that process standard output and standard error as the
    /// external process generates them.
    pub fn with_callbacks(standard: impl FnMut(&str) + 'a, error: impl FnMut(&str) + 'a) -> Self {
        Self {
            standard_callback: Box::new(standard),
            error_callback: Box::new(error),
        }
    }

    /// Re-wire the callbacks used during [`ExternalProcess::run`].
    pub fn set_callbacks(&mut self, standard: impl FnMut(&str) + 'a, error: impl FnMut(&str) + 'a) {
        self.standard_callback = Box::new(standard);
        self.error_callback = Box::new(error);
    }

    /// Run a program, calling the callbacks whenever output is available.
    ///
    /// `verbose` makes the call itself visible: before starting, the standard
    /// output callback receives [`Invocation::command_line`]; afterwards it
    /// receives `Executed '<program>' successfully!` on success, and the error
    /// callback receives the report's `error_message` otherwise. Both carry a
    /// trailing newline, as in the source.
    ///
    /// The returned [`RunReport`] carries both of the source's outputs: its
    /// `state` is the `RETURNSTATE` that both `run` overloads return, and its
    /// `error_message` is the `error_msg` that only the first fills. The second
    /// overload therefore has no separate counterpart — ignore the field.
    ///
    /// # Errors
    ///
    /// Returns whatever [`Invocation::preflight`] rejects, and [`Error::Io`]
    /// when waiting on the child fails. A program that could not be started is
    /// reported as [`ReturnState::FailedToStart`], not as an error, exactly as
    /// the source turns `bp::process_error` into that state.
    ///
    /// # Notes
    ///
    /// Chunks arrive in the order each stream produced them, but the two
    /// streams interleave nondeterministically — as they do in the source,
    /// which reads both pipes asynchronously.
    pub fn run(&mut self, invocation: &Invocation, verbose: bool) -> Result<RunReport> {
        self.run_with_idle(invocation, verbose, &mut || {})
    }

    /// [`ExternalProcess::run`] with the source's `idle_callback`, invoked once
    /// per [`POLL_INTERVAL`] while the child runs so a caller can pump a GUI
    /// event loop.
    ///
    /// # Errors
    ///
    /// As [`ExternalProcess::run`].
    pub fn run_with_idle(
        &mut self,
        invocation: &Invocation,
        verbose: bool,
        idle: &mut dyn FnMut(),
    ) -> Result<RunReport> {
        if verbose {
            (self.standard_callback)(&invocation.command_line());
        }
        let outcome = {
            let standard = &mut self.standard_callback;
            let error = &mut self.error_callback;
            let mut sink = |stream: Stream, bytes: &[u8]| -> Result<()> {
                let text = String::from_utf8_lossy(bytes);
                match stream {
                    Stream::Standard => standard(&text),
                    Stream::Error => error(&text),
                }
                Ok(())
            };
            execute(invocation, Some(idle), &mut sink)?
        };
        let report = report(invocation, outcome);
        if verbose {
            if report.state.is_success() {
                (self.standard_callback)(&format!(
                    "Executed '{}' successfully!\n",
                    invocation.program.display()
                ));
            } else {
                (self.error_callback)(&format!("{}\n", report.error_message));
            }
        }
        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_empty_program_fails_to_start_without_spawning() {
        let mut process = ExternalProcess::new();
        let report = process.run(&Invocation::new(""), false).unwrap();
        assert_eq!(report.state, ReturnState::FailedToStart);
        assert_eq!(report.exit_code, None);
        assert!(!report.timed_out);
    }

    #[test]
    fn the_verbose_banner_separates_arguments_with_single_spaces() {
        let invocation = Invocation::new("ls").with_arguments(["-l", "a b"]);
        assert_eq!(invocation.command_line(), "Running: ls -l a b\n");
    }

    #[test]
    fn preflight_rejects_an_interior_nul_in_an_argument() {
        let invocation = Invocation::new("ls").with_arguments(["a\0b"]);
        assert!(invocation.preflight().is_err());
    }

    #[test]
    fn a_crash_is_recognised_from_a_signal_on_every_platform() {
        assert!(crashed(None, Some(11)));
        assert!(!crashed(Some(0), None));
        assert!(!crashed(Some(2), None));
    }

    #[test]
    fn read_modes_are_exactly_the_two_the_source_captures_for() {
        assert!(!IoMode::NoIo.reads());
        assert!(!IoMode::WriteOnly.reads());
        assert!(IoMode::ReadOnly.reads());
        assert!(IoMode::ReadWrite.reads());
        assert_eq!(IoMode::default(), IoMode::ReadWrite);
    }
}
