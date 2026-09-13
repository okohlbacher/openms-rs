// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `SYSTEM/ExternalProcess.h`, `SYSTEM/JavaInfo.h`, `SYSTEM/PythonInfo.h` and
//! `SYSTEM/RWrapper.h`: starting external programs and probing for interpreters.
//!
//! All fourteen class-test sections across the four headers are covered here.
//!
//! Two deliberate substitutions keep the suite hermetic:
//!
//! * The `ExternalProcess` class test drives `ls -l` and `ls -0` (or `cmd /C`
//!   on Windows) because it needs *some* tool present on every box. Upstream
//!   issue #9948 then had to give it a private working directory, because under
//!   `ctest --parallel` other tests' temporary files appeared and vanished while
//!   `ls -l` walked the shared directory, making it exit non-zero. Here the
//!   three outcomes — zero exit with output on standard output, non-zero exit
//!   with output on standard error, and an executable that does not exist — come
//!   from scripts written into a private temporary directory, so the asserted
//!   shapes are the class test's and nothing depends on the host's `ls`.
//! * The interpreter probes must not require Java, Python or R. Absence is
//!   therefore tested directly on every platform, and the "interpreter present"
//!   half is exercised against a recorded stand-in script on Unix. Banner
//!   parsing is tested on recorded strings alone.
//!
//! The tests in this binary run one at a time. On Linux a child forked by one
//! test inherits every descriptor open at that instant, including the write
//! handle another test holds while it writes a stand-in script, and keeps it
//! until its own `exec`. Executing that script in the meantime fails with
//! `ETXTBSY` ("Text file busy"), which the port reports as
//! `ReturnState::FailedToStart`, as the source would. Measured on the Linux
//! build node with the spawn error printed: 18 of 25 parallel runs failed that
//! way under both Rust 1.85 and 1.96, and 0 of 20 single-threaded runs did. The
//! race is in the harness, not in `ExternalProcess`, so every test here takes
//! `serial` first; a test that only spawns can still hand a write handle to its
//! child, so the guard covers spawning tests as well as script writers.

use openms::system::external_process::{
    self, ExternalProcess, Invocation, IoMode, ReturnState, capture,
};
use openms::system::file::{FileContext, TempDir};
use openms::system::{java_info, python_info, r_wrapper};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

fn temp() -> TempDir {
    TempDir::new_in(std::env::temp_dir(), false).unwrap()
}

/// Holds the lock that keeps this binary's tests from running concurrently.
///
/// A panicking test poisons the mutex; the next test takes the lock anyway,
/// because the guarded state is `()` and the failure is already reported.
fn serial() -> std::sync::MutexGuard<'static, ()> {
    static LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
    LOCK.lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

/// A context that reaches nothing outside `base`: no `PATH`, no shared data.
fn context(base: &Path) -> FileContext {
    let mut context =
        FileContext::new(base.join("bin"), base.join("home"), base.join("temp")).unwrap();
    fs::create_dir_all(&context.executable_directory).unwrap();
    fs::create_dir_all(&context.temporary_directory).unwrap();
    context.search_path.clear();
    context.data_candidates.clear();
    context
}

#[cfg(unix)]
fn script(directory: &Path, name: &str, body: &str) -> PathBuf {
    use std::os::unix::fs::PermissionsExt;
    fs::create_dir_all(directory).unwrap();
    let path = directory.join(name);
    fs::write(&path, body).unwrap();
    fs::set_permissions(&path, fs::Permissions::from_mode(0o755)).unwrap();
    path
}

#[cfg(unix)]
const OK_SCRIPT: &str = "#!/bin/sh\necho MARKER_STDOUT_OK\nexit 0\n";
#[cfg(unix)]
const FAIL_SCRIPT: &str = "#!/bin/sh\necho MARKER_STDERR_FAIL 1>&2\nexit 2\n";

// ---------------------------------------------------------------------------
// ExternalProcess
// ---------------------------------------------------------------------------

/// Class-test section `ExternalProcess()`.
///
/// The source marks it `NOT_TESTABLE` ("tested below"); the observable property
/// of the default constructor is that its two callbacks are empty, so a child
/// that writes to both streams produces no observable effect and the state is
/// still reported.
#[cfg(unix)]
#[test]
fn a_default_constructed_process_swallows_output_and_still_reports_state() {
    let _serial = serial();
    let directory = temp();
    let ok = script(directory.path(), "ok.sh", OK_SCRIPT);
    let mut process = ExternalProcess::default();
    let report = process.run(&Invocation::new(&ok), false).unwrap();
    assert_eq!(report.state, ReturnState::Success);
    assert_eq!(report.error_message, "");
    assert_eq!(report.exit_code, Some(0));
}

/// Class-test section
/// `ExternalProcess(std::function<void(const std::string&)>, std::function<void(const std::string&)>)`.
///
/// Reproduces the class test's "with callbacks" block: on success the standard
/// output callback receives something and the error callback nothing; on a
/// non-zero exit it is the other way round.
#[cfg(unix)]
#[test]
fn constructed_callbacks_receive_the_two_streams_separately() {
    let _serial = serial();
    let directory = temp();
    let ok = script(directory.path(), "ok.sh", OK_SCRIPT);
    let fail = script(directory.path(), "fail.sh", FAIL_SCRIPT);
    let mut out = String::new();
    let mut err = String::new();
    {
        let mut process =
            ExternalProcess::with_callbacks(|text| out.push_str(text), |text| err.push_str(text));
        let report = process.run(&Invocation::new(&ok), false).unwrap();
        assert_eq!(report.state, ReturnState::Success);
    }
    assert!(out.contains("MARKER_STDOUT_OK"));
    assert!(err.is_empty());

    out.clear();
    err.clear();
    {
        let mut process =
            ExternalProcess::with_callbacks(|text| out.push_str(text), |text| err.push_str(text));
        let report = process.run(&Invocation::new(&fail), false).unwrap();
        assert_eq!(report.state, ReturnState::NonzeroExit);
        assert!(!report.error_message.is_empty());
    }
    assert!(out.is_empty());
    assert!(err.contains("MARKER_STDERR_FAIL"));
}

/// Class-test section `~ExternalProcess()`.
///
/// The source's destructor is `= default`; what a caller can observe is that
/// everything the callbacks accumulated survives the process value, because the
/// process borrows the callbacks rather than owning their state.
#[cfg(unix)]
#[test]
fn what_the_callbacks_collected_outlives_the_process_value() {
    let _serial = serial();
    let directory = temp();
    let ok = script(directory.path(), "ok.sh", OK_SCRIPT);
    let mut collected = String::new();
    {
        let mut process = ExternalProcess::with_callbacks(|t| collected.push_str(t), |_| {});
        process.run(&Invocation::new(&ok), false).unwrap();
        drop(process);
    }
    assert!(collected.contains("MARKER_STDOUT_OK"));
}

/// Class-test section
/// `void setCallbacks(std::function<void(const std::string&)>, std::function<void(const std::string&)>)`.
///
/// The class test swaps the two lambdas and re-runs the failing command: the
/// buffer that was empty now holds the standard error text and vice versa.
#[cfg(unix)]
#[test]
fn set_callbacks_swaps_which_buffer_receives_standard_error() {
    let _serial = serial();
    let directory = temp();
    let fail = script(directory.path(), "fail.sh", FAIL_SCRIPT);
    let mut first = String::new();
    let mut second = String::new();
    {
        // Swapped relative to the natural order: standard error now lands in
        // `first`, as the class test's `ep.setCallbacks(l_err, l_out)` arranges.
        let mut process = ExternalProcess::new();
        process.set_callbacks(|t| second.push_str(t), |t| first.push_str(t));
        let report = process.run(&Invocation::new(&fail), false).unwrap();
        assert_eq!(report.state, ReturnState::NonzeroExit);
    }
    assert!(first.contains("MARKER_STDERR_FAIL"));
    assert!(second.is_empty());
}

/// Class-test section
/// `RETURNSTATE run(exe, args, working_dir, verbose, error_msg, io_mode, env, idle_callback)`.
///
/// The three asserted outcomes of the class test's "without callbacks" block,
/// with the same pairing of state and `error_msg` emptiness:
/// `SUCCESS` with an empty message, `FAILED_TO_START` with a non-empty one, and
/// `NONZERO_EXIT` with a non-empty one.
#[cfg(unix)]
#[test]
fn run_reports_success_failure_to_start_and_a_nonzero_exit() {
    let _serial = serial();
    let directory = temp();
    let ok = script(directory.path(), "ok.sh", OK_SCRIPT);
    let fail = script(directory.path(), "fail.sh", FAIL_SCRIPT);
    let mut process = ExternalProcess::new();

    let report = process
        .run(
            &Invocation::new(&ok).with_working_directory(directory.path()),
            true,
        )
        .unwrap();
    assert_eq!(report.state, ReturnState::Success);
    assert_eq!(report.error_message.len(), 0);

    let report = process
        .run(
            &Invocation::new("this_exe_does_not_exist").with_working_directory(directory.path()),
            true,
        )
        .unwrap();
    assert_eq!(report.state, ReturnState::FailedToStart);
    assert_ne!(report.error_message.len(), 0);
    assert_eq!(
        report.error_message,
        "Process 'this_exe_does_not_exist' failed to start. Does it exist? Is it executable?"
    );

    let report = process
        .run(
            &Invocation::new(&fail).with_working_directory(directory.path()),
            true,
        )
        .unwrap();
    assert_eq!(report.state, ReturnState::NonzeroExit);
    assert_ne!(report.error_message.len(), 0);
    assert!(report.error_message.contains("(exit code: 2)"));
}

/// Class-test section
/// `RETURNSTATE run(exe, args, working_dir, verbose, io_mode, env, idle_callback)`.
///
/// The source marks the message-less overload `NOT_TESTABLE` ("tested above").
/// It has no separate counterpart here because one [`openms::system::external_process::RunReport`]
/// carries both of the first overload's outputs, which this pins: a caller that
/// ignores `error_message` sees exactly the state the other overload returned.
#[cfg(unix)]
#[test]
fn one_report_carries_both_overloads_outputs() {
    let _serial = serial();
    let directory = temp();
    let fail = script(directory.path(), "fail.sh", FAIL_SCRIPT);
    let mut process = ExternalProcess::new();
    let report = process.run(&Invocation::new(&fail), false).unwrap();
    assert_eq!(report.state, ReturnState::NonzeroExit);
    assert!(!report.is_success());
    assert!(!report.error_message.is_empty());
}

/// Class-test section `[EXTRA] run with spaces in the executable path and arguments`.
///
/// An executable whose path contains spaces must launch, its standard output
/// must be captured intact, and a single argument that itself contains spaces
/// must arrive as one argument rather than being split by a shell. Reproduces
/// the class test's two literals, `MARKER_STDOUT_OK` and `arg=[one two three]`.
///
/// The section's own call is
/// `ep.run(script.string(), spaced_args, "", true, error_msg)` — an **empty**
/// working directory, which the source folds into `"."` — and expects
/// `SUCCESS`, so the empty string is passed here rather than left unset.
#[cfg(unix)]
#[test]
fn spaces_in_the_executable_path_and_in_one_argument_survive() {
    let _serial = serial();
    let directory = temp();
    let spaced = directory.path().join("open ms space dir");
    let path = script(
        &spaced,
        "my script.sh",
        "#!/bin/sh\necho MARKER_STDOUT_OK\necho \"arg=[$1]\"\n",
    );
    let mut out = String::new();
    {
        let mut process = ExternalProcess::with_callbacks(|t| out.push_str(t), |_| {});
        let report = process
            .run(
                &Invocation::new(&path)
                    .with_arguments(["one two three"])
                    .with_working_directory(""),
                true,
            )
            .unwrap();
        assert_eq!(report.state, ReturnState::Success);
    }
    assert!(out.contains("MARKER_STDOUT_OK"));
    assert!(out.contains("arg=[one two three]"));
}

/// Native: the source's crash branch, which the class test never reaches.
///
/// `WIFSIGNALED` is `ExitStatus::signal().is_some()`, so a child killed by a
/// signal is `CRASH` rather than a non-zero exit. It is reported for every
/// [`openms::system::external_process::IoMode`] here; the source only asks the
/// question in its reading branch, so under `NO_IO` or `WRITE_ONLY` it reports
/// the signal number as if it were an exit code.
#[cfg(unix)]
#[test]
fn a_child_killed_by_a_signal_is_a_crash_in_every_io_mode() {
    let _serial = serial();
    let directory = temp();
    let path = script(directory.path(), "crash.sh", "#!/bin/sh\nkill -9 $$\n");
    for mode in [IoMode::ReadWrite, IoMode::NoIo] {
        let run = capture(&Invocation::new(&path).with_io_mode(mode)).unwrap();
        assert_eq!(run.report.state, ReturnState::Crash);
        assert_eq!(run.report.terminating_signal, Some(9));
        assert_eq!(run.report.exit_code, None);
        assert!(run.report.error_message.contains("crashed hard"));
    }
}

/// Native: [`openms::system::external_process::Invocation::timeout`], which the
/// source's `ExternalProcess` lacks and its interpreter probes implement
/// separately as `Internal::waitForProcess`.
#[cfg(unix)]
#[test]
fn a_child_that_overruns_its_budget_is_terminated() {
    let _serial = serial();
    let directory = temp();
    let path = script(directory.path(), "sleep.sh", "#!/bin/sh\nsleep 30\n");
    let run = capture(
        &Invocation::new(&path)
            .with_io_mode(IoMode::NoIo)
            .with_timeout(Duration::from_millis(250)),
    )
    .unwrap();
    assert!(run.report.timed_out);
    assert!(!run.report.is_success());
    assert!(run.report.error_message.contains("did not finish in time"));
}

/// The source's `env` parameter: entries are *added* to the inherited
/// environment rather than replacing it.
#[cfg(unix)]
#[test]
fn the_environment_overlay_is_added_to_the_inherited_environment() {
    let _serial = serial();
    let directory = temp();
    let path = script(
        directory.path(),
        "env.sh",
        "#!/bin/sh\necho \"added=[$OPENMS_PORT_TEST]\"\necho \"inherited=[${PATH:+yes}]\"\n",
    );
    let mut environment = BTreeMap::new();
    environment.insert("OPENMS_PORT_TEST".to_owned(), "value with space".to_owned());
    let run = capture(&Invocation::new(&path).with_environment(environment)).unwrap();
    assert!(run.stdout.contains("added=[value with space]"));
    assert!(run.stdout.contains("inherited=[yes]"));
}

/// The source's `working_dir` parameter, and its "leave empty for the current
/// working directory" rule, which this port spells `None` or an empty path.
#[cfg(unix)]
#[test]
fn the_working_directory_is_where_relative_paths_resolve() {
    let _serial = serial();
    let directory = temp();
    let path = script(directory.path(), "cat.sh", "#!/bin/sh\ncat ./marker.txt\n");
    fs::write(
        directory.path().join("marker.txt"),
        "in the working directory",
    )
    .unwrap();

    let run = capture(&Invocation::new(&path).with_working_directory(directory.path())).unwrap();
    assert_eq!(run.report.state, ReturnState::Success);
    assert_eq!(run.stdout, "in the working directory");

    // Without a working directory the relative `marker.txt` is not found, so
    // `cat` fails; that is the source's `"."` default.
    let run = capture(&Invocation::new(&path)).unwrap();
    assert_eq!(run.report.state, ReturnState::NonzeroExit);
}

/// The other half of that rule: `working_dir.empty() ? "." : working_dir`
/// (`ExternalProcess.cpp:131`, `:233`). An empty directory is not a directory
/// that cannot be entered — it is the current one, so the child must run.
#[cfg(unix)]
#[test]
fn an_empty_working_directory_is_the_current_one() {
    let _serial = serial();
    let directory = temp();
    let path = script(directory.path(), "pwd.sh", "#!/bin/sh\npwd\n");

    let empty = capture(&Invocation::new(&path).with_working_directory("")).unwrap();
    assert_eq!(empty.report.state, ReturnState::Success);

    let unset = capture(&Invocation::new(&path)).unwrap();
    assert_eq!(empty.stdout, unset.stdout);
}

/// The source's `idle_callback`, invoked from the poll loop so a caller can
/// pump a GUI event loop while the child runs.
#[cfg(unix)]
#[test]
fn the_idle_callback_runs_while_the_child_runs() {
    let _serial = serial();
    let directory = temp();
    let path = script(directory.path(), "slow.sh", "#!/bin/sh\nsleep 0.4\n");
    let mut ticks = 0_usize;
    let mut process = ExternalProcess::new();
    let report = process
        .run_with_idle(
            &Invocation::new(&path).with_io_mode(IoMode::NoIo),
            false,
            &mut || ticks += 1,
        )
        .unwrap();
    assert_eq!(report.state, ReturnState::Success);
    // At a 50 ms cadence a 400 ms child yields several ticks; one is enough to
    // prove the callback is reached from inside the loop.
    assert!(ticks >= 1, "idle callback never ran");
}

/// The cadence itself. The source reaches `idle_callback` once per
/// `io_ctx.run_for(milliseconds(50))`, however much the child writes in
/// between, so the callback rate is bounded by elapsed time and not by the
/// number of output chunks. A child that floods both is the case that tells
/// the two apart.
#[cfg(unix)]
#[test]
fn the_idle_cadence_is_set_by_the_clock_and_not_by_the_output() {
    let _serial = serial();
    let directory = temp();
    let path = script(
        directory.path(),
        "chatty.sh",
        "#!/bin/sh\nyes 0123456789012345678901234567890123456789 | head -n 60000\nsleep 0.3\n",
    );
    let mut ticks = 0_usize;
    let mut process = ExternalProcess::new();
    let started = std::time::Instant::now();
    let report = process
        .run_with_idle(&Invocation::new(&path), false, &mut || ticks += 1)
        .unwrap();
    let elapsed = started.elapsed();
    assert_eq!(report.state, ReturnState::Success);

    // Every tick needs a full POLL_INTERVAL since the previous one, so the
    // count is at most one more than the intervals that fit in the run. The
    // spare tick absorbs the partial interval at each end.
    let allowed = elapsed.as_millis() / external_process::POLL_INTERVAL.as_millis() + 2;
    assert!(
        u128::try_from(ticks).unwrap() <= allowed,
        "idle callback ran {ticks} times in {elapsed:?}, which is more often than the \
         source's 50 ms poll allows ({allowed}); it is following the output, not the clock"
    );
}

/// A budget has to bound the *call*, not merely the child. A child that exits
/// at once but leaves a descendant behind hands that descendant the write ends
/// of its pipes, so they never reach end of file; without this rule the drain
/// waits for the descendant and the budget can never fire, because the child it
/// would have killed is already gone.
///
/// The call is made on a worker thread so that the failure is a failed
/// assertion rather than a suite that never finishes.
#[cfg(unix)]
#[test]
fn a_budget_ends_the_call_when_a_descendant_still_holds_the_pipes() {
    let _serial = serial();
    let directory = temp();
    let path = script(
        directory.path(),
        "orphan.sh",
        "#!/bin/sh\nsleep 30 &\nexit 0\n",
    );
    let (sender, receiver) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let run = capture(&Invocation::new(&path).with_timeout(Duration::from_millis(250)));
        let _ = sender.send(run.map(|run| run.report));
    });

    let report = receiver
        .recv_timeout(Duration::from_secs(20))
        .expect(
            "the call never returned: the budget cannot end a drain a descendant is holding open",
        )
        .unwrap();
    // The child itself exited cleanly, and it is the child's outcome that is
    // reported: nothing was killed, so `timed_out` stays false.
    assert_eq!(report.state, ReturnState::Success);
    assert_eq!(report.exit_code, Some(0));
    assert!(!report.timed_out);
}

/// Verbose mode writes the source's two banners: the command line before the
/// child starts, and the success line after it finishes.
#[cfg(unix)]
#[test]
fn verbose_mode_brackets_the_call_with_the_sources_two_banners() {
    let _serial = serial();
    let directory = temp();
    let ok = script(directory.path(), "ok.sh", OK_SCRIPT);
    let mut out = String::new();
    {
        let mut process = ExternalProcess::with_callbacks(|t| out.push_str(t), |_| {});
        process
            .run(&Invocation::new(&ok).with_arguments(["a b"]), true)
            .unwrap();
    }
    assert!(out.starts_with(&format!("Running: {} a b\n", ok.display())));
    assert!(out.ends_with(&format!("Executed '{}' successfully!\n", ok.display())));
}

/// A program that does not exist is a state, not an error, on every platform —
/// the branch the class test reaches with `"this_exe_does_not_exist"`.
#[test]
fn a_missing_program_is_a_state_and_never_an_error() {
    let _serial = serial();
    let mut process = ExternalProcess::new();
    let report = process
        .run(&Invocation::new("this_exe_does_not_exist_@@"), false)
        .unwrap();
    assert_eq!(report.state, ReturnState::FailedToStart);
    assert!(report.error_message.contains("failed to start"));
    assert_eq!(report.exit_code, None);
    assert!(!report.timed_out);
}

/// Bounded work: every ceiling is checked before anything is spawned, so an
/// oversized request leaves the system untouched.
#[test]
fn every_ceiling_is_checked_before_a_child_exists() {
    let _serial = serial();
    let too_many = Invocation::new("true")
        .with_arguments(vec!["x".to_owned(); external_process::MAX_ARGUMENTS + 1]);
    assert!(too_many.preflight().is_err());
    assert!(external_process::capture(&too_many).is_err());

    let too_long = Invocation::new("true")
        .with_arguments(["x".repeat(external_process::MAX_COMMAND_BYTES + 1)]);
    assert!(too_long.preflight().is_err());

    let mut environment = BTreeMap::new();
    environment.insert("BAD=KEY".to_owned(), "v".to_owned());
    assert!(
        Invocation::new("true")
            .with_environment(environment)
            .preflight()
            .is_err()
    );

    let mut environment = BTreeMap::new();
    for index in 0..=external_process::MAX_ENVIRONMENT_ENTRIES {
        environment.insert(format!("K{index}"), String::new());
    }
    assert!(
        Invocation::new("true")
            .with_environment(environment)
            .preflight()
            .is_err()
    );
}

// ---------------------------------------------------------------------------
// JavaInfo
// ---------------------------------------------------------------------------

/// Class-test section `static bool canRun(const std::string& file)`.
///
/// The single assertion is `JavaInfo::canRun("") == false`. Reproduced here on
/// every platform, together with the diagnosis the source would have logged.
#[test]
fn java_cannot_run_from_an_empty_executable_name() {
    let _serial = serial();
    let directory = temp();
    let check = java_info::can_run(&context(directory.path()), "").unwrap();
    assert!(!check.can_run);
    assert!(check.message.starts_with("Java-Check:\n"));
    assert!(check.message.contains("Java not found at ''!"));
    assert_eq!(check.exit_code, None);
    assert!(!check.timed_out);
}

/// A Java that *is* present, without requiring Java: a recorded stand-in on the
/// context's search path answers `-version` on standard error, as a real JVM
/// does, and a second stand-in exits non-zero to reach the source's
/// "returned a non-zero exit code" branch.
#[cfg(unix)]
#[test]
fn java_is_probed_through_a_recorded_stand_in_on_the_search_path() {
    let _serial = serial();
    let directory = temp();
    let good = directory.path().join("good");
    let bad = directory.path().join("bad");
    script(
        &good,
        "java",
        "#!/bin/sh\necho 'openjdk version \"21.0.2\" 2024-01-16' 1>&2\nexit 0\n",
    );
    script(&bad, "java", "#!/bin/sh\nexit 3\n");

    let mut ok = context(directory.path());
    ok.search_path.push(good);
    let check = java_info::can_run(&ok, "java").unwrap();
    assert!(check.can_run);
    assert_eq!(check.message, "");
    assert_eq!(check.exit_code, Some(0));
    assert!(check.output.contains("openjdk version \"21.0.2\""));

    let mut broken = context(directory.path());
    broken.search_path.push(bad);
    let check = java_info::can_run(&broken, "java").unwrap();
    assert!(!check.can_run);
    assert_eq!(check.exit_code, Some(3));
    assert!(
        check
            .message
            .contains("Java returned a non-zero exit code (3).")
    );
}

// ---------------------------------------------------------------------------
// PythonInfo
// ---------------------------------------------------------------------------

#[cfg(unix)]
const FAKE_PYTHON: &str = "#!/bin/sh\ncase \"$1\" in\n  --version) echo 'Python 3.11.4'; exit 0 ;;\n  -c) case \"$2\" in 'import math') exit 0 ;; *) echo 'ModuleNotFoundError' 1>&2; exit 1 ;; esac ;;\nesac\nexit 2\n";

/// Class-test section `static bool canRun(std::string& python_executable, std::string& error_msg)`.
///
/// Both of the class test's message literals are reproduced: `"Python not found
/// at"` for a name that resolves nowhere, and `"failed to run"` for a file that
/// exists but cannot be executed. The third block — run only when a real Python
/// is present — asserted that the executable was rewritten to an existing,
/// absolute path; the stand-in interpreter lets that be asserted unconditionally
/// on Unix.
#[test]
fn python_can_run_reproduces_the_two_message_literals() {
    let _serial = serial();
    let directory = temp();
    let context = context(directory.path());

    let check = python_info::can_run(&context, "does_not_exist_@@").unwrap();
    assert!(!check.can_run);
    assert!(check.message.contains("Python not found at"));

    let not_executable = directory.path().join("empty_file");
    fs::write(&not_executable, "").unwrap();
    let check = python_info::can_run(&context, not_executable.to_str().unwrap()).unwrap();
    assert!(!check.can_run);
    assert!(check.message.contains("failed to run"));
}

/// The "Python is present" half of the same section, without Python.
#[cfg(unix)]
#[test]
fn python_can_run_resolves_a_bare_name_to_an_absolute_path() {
    let _serial = serial();
    let directory = temp();
    let bin = directory.path().join("bin_python");
    script(&bin, "python", FAKE_PYTHON);
    let mut context = context(directory.path());
    context.search_path.push(bin.clone());

    let check = python_info::can_run(&context, "python").unwrap();
    assert!(check.can_run);
    assert_eq!(check.executable, bin.join("python"));
    assert!(check.executable.is_absolute());
    assert!(check.executable.exists());
    // The source appends this note to `error_msg` even on success.
    assert!(
        check
            .message
            .starts_with("Python executable ('python') resolved to '")
    );
}

/// Class-test section
/// `bool PythonInfo::isPackageInstalled(const std::string&, const std::string&)`.
///
/// The class test asserts `false` for `"veryWeirdPackage___@@__@"` and `true`
/// for `"math"`. Both are reproduced: the weird name is refused as a module
/// path before any interpreter starts, which is the same answer by a safer
/// route, and `math` is answered by the stand-in interpreter.
#[test]
fn python_package_detection_answers_the_class_tests_two_values() {
    let _serial = serial();
    let directory = temp();
    let context = context(directory.path());
    assert!(
        !python_info::is_package_installed(
            &context,
            "does_not_exist_@@",
            "veryWeirdPackage___@@__@"
        )
        .unwrap()
    );
    assert!(python_info::validate_package_name("veryWeirdPackage___@@__@").is_err());
}

/// The `true` half of the same section, plus the hardening the source lacks.
#[cfg(unix)]
#[test]
fn python_package_detection_uses_the_stand_in_interpreter() {
    let _serial = serial();
    let directory = temp();
    let bin = directory.path().join("bin_python");
    script(&bin, "python", FAKE_PYTHON);
    let mut context = context(directory.path());
    context.search_path.push(bin);

    assert!(python_info::is_package_installed(&context, "python", "math").unwrap());
    assert!(!python_info::is_package_installed(&context, "python", "absent_module").unwrap());
    // The source would have executed this as Python source; here it never runs.
    assert!(!python_info::is_package_installed(&context, "python", "os; import sys").unwrap());
}

/// Class-test section `static std::string getVersion(const std::string& python_executable)`.
///
/// The class test only asserts that the version is non-empty when Python is
/// present. The stand-in interpreter lets the exact banner be asserted, and the
/// recorded strings pin the source's concatenate-then-trim rule independently.
#[test]
fn python_version_parsing_follows_the_sources_concatenate_then_trim() {
    let _serial = serial();
    assert_eq!(
        python_info::version_from_output("Python 3.11.4\n", ""),
        "Python 3.11.4"
    );
    assert_eq!(
        python_info::version_from_output("", "Python 2.7.18\r\n"),
        "Python 2.7.18"
    );
    // Two lines are joined with nothing between them: the source appends each
    // `std::getline` result without a separator.
    assert_eq!(
        python_info::version_from_output("Python 3.11.4\ntrailing\n", ""),
        "Python 3.11.4trailing"
    );
    assert_eq!(python_info::version_from_output("", ""), "");
}

/// The live half of `getVersion`, against the stand-in, plus the empty-string
/// contract for an interpreter that is not there.
#[cfg(unix)]
#[test]
fn python_version_is_read_from_the_stand_in_and_empty_when_absent() {
    let _serial = serial();
    let directory = temp();
    let bin = directory.path().join("bin_python");
    script(&bin, "python", FAKE_PYTHON);
    let mut context = context(directory.path());
    context.search_path.push(bin);
    assert_eq!(
        python_info::version(&context, "python").unwrap(),
        "Python 3.11.4"
    );
    assert_eq!(
        python_info::version(&context, "does_not_exist_@@").unwrap(),
        ""
    );
}

// ---------------------------------------------------------------------------
// RWrapper
// ---------------------------------------------------------------------------

#[cfg(unix)]
const FAKE_RSCRIPT: &str = "#!/bin/sh\nfor a in \"$@\"; do echo \"[$a]\"; done\nexit 0\n";

/// Class-test section `static bool findR(const std::string& executable, bool verbose)`.
///
/// A non-existent interpreter is reported as "not found" cleanly — the source
/// catches `boost::process`'s `process_error` and returns `false`, never letting
/// an exception escape. Probing for the real `Rscript` must likewise return a
/// clean answer either way, which here is `Ok` rather than a panic.
#[test]
fn find_r_reports_a_bogus_interpreter_cleanly() {
    let _serial = serial();
    let directory = temp();
    let context = context(directory.path());
    let check = r_wrapper::find_r(&context, "this_is_not_a_real_R_interpreter_xyz").unwrap();
    assert!(!check.found);
    assert_eq!(check.executable, None);
    assert!(check.message.contains("FailedToStart"));

    // The "is R available?" signal itself: an answer, not an error.
    let check = r_wrapper::find_r(&context, r_wrapper::DEFAULT_EXECUTABLE).unwrap();
    assert!(!check.found);
}

/// The interpreter-present half of `findR`, without R.
#[cfg(unix)]
#[test]
fn find_r_accepts_a_recorded_stand_in_interpreter() {
    let _serial = serial();
    let directory = temp();
    let bin = directory.path().join("bin_r");
    script(&bin, "Rscript", FAKE_RSCRIPT);
    let mut context = context(directory.path());
    context.search_path.push(bin.clone());

    let check = r_wrapper::find_r(&context, r_wrapper::DEFAULT_EXECUTABLE).unwrap();
    assert!(check.found);
    assert_eq!(check.executable, Some(bin.join("Rscript")));
    // The source probes with exactly these three arguments.
    assert_eq!(check.output, "[--vanilla]\n[-e]\n[sessionInfo()]\n");
    assert_eq!(check.message, "");
}

/// Class-test section `static std::string findScript(const std::string& script_file, bool verbose)`.
///
/// A non-existent script raises a clear not-found error — the source's
/// `Exception::FileNotFound`, which is this crate's `Error::Io` with
/// `ErrorKind::NotFound`. This is the only `RWrapper` entry point that fails
/// loudly; `runScript` swallows it.
#[test]
fn find_script_reports_a_missing_script_as_not_found() {
    let _serial = serial();
    let directory = temp();
    let error = r_wrapper::find_script(
        &context(directory.path()),
        "definitely_nonexistent_script_qwerty.R",
    )
    .unwrap_err();
    let openms::Error::Io(io) = error else {
        panic!("expected an I/O error");
    };
    assert_eq!(io.kind(), std::io::ErrorKind::NotFound);
}

/// The success half of `findScript`: a bundled script resolves under
/// `share/OpenMS/SCRIPTS`.
#[test]
fn find_script_resolves_a_bundled_script_under_the_scripts_directory() {
    let _serial = serial();
    let directory = temp();
    let share = directory.path().join("share");
    fs::create_dir_all(share.join("CHEMISTRY")).unwrap();
    fs::write(share.join("CHEMISTRY/unimod.xml"), "marker").unwrap();
    fs::create_dir_all(share.join(r_wrapper::SCRIPTS_DIRECTORY)).unwrap();
    fs::write(share.join("SCRIPTS/plot.R"), "cat('hi')\n").unwrap();

    let mut context = context(directory.path());
    context.data_override = Some(share.clone());
    let found = r_wrapper::find_script(&context, "plot.R").unwrap();
    assert_eq!(found, share.join("SCRIPTS/plot.R"));
}

/// Class-test section
/// `static bool runScript(const std::string&, const std::vector<std::string>&, const std::string&, bool, bool)`.
///
/// Both degradation paths the class test pins: a bogus interpreter with
/// `find_R = true` returns `false` cleanly, and a missing script returns `false`
/// cleanly because the internal not-found error is swallowed. Neither escapes as
/// an error.
#[test]
fn run_script_degrades_to_false_without_raising() {
    let _serial = serial();
    let directory = temp();
    let context = context(directory.path());

    let outcome = r_wrapper::run_script(
        &context,
        "any.R",
        &[],
        "this_is_not_a_real_R_interpreter_xyz",
        true,
    )
    .unwrap();
    assert!(!outcome.success);
    assert_eq!(outcome.script, None);

    let outcome = r_wrapper::run_script(
        &context,
        "definitely_nonexistent_script_qwerty.R",
        &[],
        r_wrapper::DEFAULT_EXECUTABLE,
        false,
    )
    .unwrap();
    assert!(!outcome.success);
    assert!(outcome.message.contains("Could not find R script"));
}

/// The success path of `runScript`, against the stand-in interpreter: the
/// source's fixed `--vanilla --quiet` prefix, the resolved script, and the
/// caller's arguments, each passed as one argument even when it contains
/// spaces.
#[cfg(unix)]
#[test]
fn run_script_passes_the_sources_argument_vector_intact() {
    let _serial = serial();
    let directory = temp();
    let share = directory.path().join("share");
    fs::create_dir_all(share.join("CHEMISTRY")).unwrap();
    fs::write(share.join("CHEMISTRY/unimod.xml"), "marker").unwrap();
    fs::create_dir_all(share.join(r_wrapper::SCRIPTS_DIRECTORY)).unwrap();
    fs::write(share.join("SCRIPTS/plot.R"), "cat('hi')\n").unwrap();
    let bin = directory.path().join("bin_r");
    script(&bin, "Rscript", FAKE_RSCRIPT);

    let mut context = context(directory.path());
    context.data_override = Some(share.clone());
    context.search_path.push(bin);

    let outcome = r_wrapper::run_script(
        &context,
        "plot.R",
        &["out with space.png".to_owned()],
        r_wrapper::DEFAULT_EXECUTABLE,
        true,
    )
    .unwrap();
    assert!(outcome.success);
    assert_eq!(outcome.script, Some(share.join("SCRIPTS/plot.R")));
    assert_eq!(
        outcome.stdout,
        format!(
            "[--vanilla]\n[--quiet]\n[{}]\n[out with space.png]\n",
            share.join("SCRIPTS/plot.R").display()
        )
    );
    assert_eq!(outcome.message, "");
}

/// The source drains its pipes line by line and appends `'\n'` to each line;
/// `reassemble_lines` reproduces that, so the strings a caller inspects are the
/// ones the C++ assembled.
#[test]
fn captured_output_is_reassembled_the_way_the_source_drains_it() {
    let _serial = serial();
    assert_eq!(r_wrapper::reassemble_lines("a\nb\n"), "a\nb\n");
    assert_eq!(r_wrapper::reassemble_lines("a\nb"), "a\nb\n");
    assert_eq!(r_wrapper::reassemble_lines("\n"), "\n");
    assert_eq!(r_wrapper::reassemble_lines(""), "");
}

/// Bounded work in `runScript`: the argument ceiling is checked before the
/// interpreter is looked up.
#[test]
fn run_script_refuses_an_oversized_argument_list() {
    let _serial = serial();
    let directory = temp();
    let many = vec![String::from("x"); r_wrapper::MAX_SCRIPT_ARGUMENTS + 1];
    assert!(
        r_wrapper::run_script(
            &context(directory.path()),
            "any.R",
            &many,
            r_wrapper::DEFAULT_EXECUTABLE,
            false,
        )
        .is_err()
    );
}
