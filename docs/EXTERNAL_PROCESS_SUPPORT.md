# Starting an external program

`system::external_process` covers `SYSTEM/ExternalProcess.h` and
`SYSTEM/ExternalProcess.cpp` at SDK `bc9cc12`, and the private
`source/SYSTEM/ProcessWait.h` helper that its sibling probes use.

The C++ wraps `boost::process` and `boost::asio`. Everything it does —
resolving a program on `PATH`, spawning with an argument vector, setting a
working directory, overlaying environment variables, piping both output
streams, polling for exit, and distinguishing a crash from a non-zero exit —
is available from `std::process::Command` plus `std::thread` and
`std::sync::mpsc`, so **no dependency is added**.

## API mapping

Every public member of the header, with its Rust counterpart.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `enum class RETURNSTATE` | `external_process::ReturnState` | |
| `RETURNSTATE::SUCCESS` | `ReturnState::Success` | Also `ReturnState::is_success`, `RunReport::is_success`. |
| `RETURNSTATE::NONZERO_EXIT` | `ReturnState::NonzeroExit` | |
| `RETURNSTATE::CRASH` | `ReturnState::Crash` | |
| `RETURNSTATE::FAILED_TO_START` | `ReturnState::FailedToStart` | |
| `enum class IO_MODE` | `external_process::IoMode` | |
| `IO_MODE::NO_IO` | `IoMode::NoIo` | |
| `IO_MODE::READ_ONLY` | `IoMode::ReadOnly` | |
| `IO_MODE::WRITE_ONLY` | `IoMode::WriteOnly` | |
| `IO_MODE::READ_WRITE` | `IoMode::ReadWrite` | The source's default argument; also this port's `Default`. |
| `ExternalProcess()` | `ExternalProcess::new` / `ExternalProcess::default` | Two callbacks that do nothing, as the source's delegating constructor installs. |
| `ExternalProcess(std::function<…> out, std::function<…> err)` | `ExternalProcess::with_callbacks` | `std::function<void(const std::string&)>` → `impl FnMut(&str) + 'a`. |
| `~ExternalProcess()` | — | Implicit on both sides; the source's destructor is `= default`. |
| `void setCallbacks(out, err)` | `ExternalProcess::set_callbacks` | |
| `RETURNSTATE run(exe, args, working_dir, verbose, error_msg, io_mode, env, idle_callback)` | `ExternalProcess::run_with_idle` | The eight parameters become `&Invocation` + `verbose` + `idle`; `error_msg` is `RunReport::error_message`. |
| `RETURNSTATE run(exe, args, working_dir, verbose, io_mode, env, idle_callback)` | `ExternalProcess::run` | The message-less overload. One `RunReport` carries both outputs, so the only difference is whether the caller reads `error_message`. |
| `callbackStdOut_`, `callbackStdErr_` (private) | `ExternalProcess` private fields | |
| `Internal::waitForProcess(child, timeout)` (`source/SYSTEM/ProcessWait.h`) | `Invocation::timeout` + the run loop | The source's polling deadline helper, folded into the one loop instead of a second one. Its 30 s budget is `PROBE_TIMEOUT`. |

Native additions, each documented at the item:

| Rust item | Why |
|---|---|
| `Invocation` (+ `new`, `with_arguments`, `with_working_directory`, `with_environment`, `with_io_mode`, `with_timeout`, `command_line`, `preflight`) | The eight `run` parameters as one value, so the two overloads do not become four methods and the ceilings have one place to live. |
| `RunReport` (+ `is_success`), fields `state`, `error_message`, `exit_code`, `terminating_signal`, `timed_out` | `RETURNSTATE` plus the `error_msg` out-parameter, plus the detail behind `CRASH` that the source computes and discards. |
| `CapturedRun`, `capture` | The accumulating shape the source's three interpreter probes each re-implement with `bp::ipstream`; here it is written once and they call it. |
| `Stream` | Tags a chunk with its origin so one sink can serve both the callback and the accumulating form. |
| `ReturnState::is_success`, `IoMode::reads` | `IoMode::reads` is the source's local `can_read` predicate, promoted so the collapse of four modes into two is visible rather than buried. |
| `MAX_ARGUMENTS`, `MAX_COMMAND_BYTES`, `MAX_ENVIRONMENT_ENTRIES`, `MAX_CAPTURED_BYTES` | Bounded work; the source has no ceilings. |
| `POLL_INTERVAL`, `READ_BUFFER_BYTES`, `PROBE_TIMEOUT` | The source's three magic numbers — 50 ms, 4096 bytes, 30 s — named. |

`ExternalProcessMBox`, which the header's documentation points at as "a
convenient wrapper" that raises Qt message boxes, lives outside the scientific
SDK and has no counterpart here.

## Preserved source conventions

**The four outcomes and their messages are exact.** `error_msg` is cleared at
entry and set only on failure, and the three failure strings are transcribed
character for character:

- `Process '<exe>' failed to start. Does it exist? Is it executable?`
- `Process '<exe>' crashed hard (segfault-like). Please check the log.`
- `Process '<exe>' did not finish successfully (exit code: <n>). Please check the log.`

**The verbose banners are exact**, including their trailing newlines: the
command line is `Running: <exe><space><arg>…` on the stdout callback before the
child starts, and `Executed '<exe>' successfully!` on the same callback after
it succeeds; a failure sends `error_msg + "\n"` to the stderr callback.

**Arguments are a vector, never a shell string.** The source passes
`bp::args(std_args)` and never builds a command string for execution; this
passes `Command::args`. A path or argument containing spaces therefore arrives
intact, which the class test's `[EXTRA]` section asserts with
`arg=[one two three]` and which `spaces_in_the_executable_path_and_in_one_argument_survive`
reproduces. `Invocation::command_line` builds a display string with the same
naive space-joining the source uses for its verbose banner, and that string is
never executed.

**The environment is inherited and overlaid, not replaced**
(`ExternalProcess.cpp:113`–`116`). A `BTreeMap` keeps application order
deterministic, matching `std::map`.

**An empty working directory means `"."`** (`ExternalProcess.cpp:131`, `:233`),
which `Invocation::working_directory == None` spells.

**The poll cadence is 50 ms** and the read buffer is 4096 bytes, so a caller's
`idle_callback` is reached at the source's rate.

**`IO_MODE` collapses to one question.** The source computes
`can_read = (io_mode == READ_ONLY || io_mode == READ_WRITE)`
(`ExternalProcess.cpp:122`) and branches on nothing else; neither branch ever
wires the child's standard input, so `WRITE_ONLY` is indistinguishable from
`NO_IO`. `IoMode::reads` reproduces that exactly and says so.

## Native differences

**`PATH` resolution comes from the standard library.** The source calls
`bp::search_path(exe)` and falls back to the unresolved string when that returns
empty (`ExternalProcess.cpp:106`); `Command::new` performs the platform's own
lookup, which also handles an absolute or separator-bearing path when `PATH` is
empty. The observable difference is confined to that last case, where the source
would fail to start; see `OpenMS_CPP_ISSUES.md` entry *empty `PATH` defeats an
absolute executable*.

**A crash is detected in every `IO_MODE`.** The source asks `WIFSIGNALED` only
inside its reading branch (`ExternalProcess.cpp:209`); the non-reading branch
tests `exit_code != 0` alone, and since Boost's `exit_code()` yields `WTERMSIG`
for a signalled child, a segfault under `NO_IO` is reported as
`NONZERO_EXIT (exit code: 11)`. This port classifies identically in both modes,
because reporting a crash as an ordinary non-zero exit is a wrong answer rather
than a different one. `a_child_killed_by_a_signal_is_a_crash_in_every_io_mode`
pins both modes.

**The Windows crash predicate is simplified, provably.** The source writes
`exit_code < 0 || static_cast<unsigned int>(exit_code) > 0x80000000u`
(`ExternalProcess.cpp:205`). For a non-negative `i32` the unsigned value is at
most `0x7fffffff`, and for a negative one the first disjunct already holds, so
the second can never decide the outcome; the port writes `exit_code < 0`.

**Output reaches callbacks as UTF-8, lossily.** The source hands raw bytes to a
`std::string`; `String::from_utf8_lossy` replaces invalid sequences rather than
failing a probe over a byte in a version banner.

**Every ceiling is checked before a child exists.** `Invocation::preflight`
rejects too many arguments, too many environment entries, an oversized command
line, an empty or `=`-bearing environment key, and any interior NUL byte. The
source checks none of these, and an interior NUL in particular would silently
truncate the string handed to `exec`. A rejected invocation spawns nothing.

**`run` returns `Result`, and a program that will not start is not an error.**
`Err` is reserved for a rejected preflight, a capture ceiling, and a failed
wait; `FAILED_TO_START` stays a `ReturnState`, as in the source.

**Threading.** Neither file carries `#pragma omp`; there is no data parallelism
to mirror. The source uses `boost::asio` to read both pipes without blocking,
and `RWrapper` uses one `std::thread` for the same reason. This port spawns one
reader thread per pipe and delivers their chunks to the caller's thread through
a channel, so callbacks are still invoked from the calling thread only. Within a
stream, chunk order is deterministic; between the two streams it is not, exactly
as in the source.

**Timeouts are opt-in.** `ExternalProcess::run` has no budget, as the source
has none. `Invocation::timeout` exists because the source's own probes need one
(`ProcessWait.h`), and killing on expiry mirrors their `terminate()` + `wait()`.

## Checked boundaries and evidence

| Boundary | Behaviour |
|---|---|
| Empty program name | `FailedToStart` without spawning; the source reaches the same answer through `search_path("") → bp::child("") → process_error`. |
| Program not on `PATH` | `FailedToStart` with the source's message. |
| Interior NUL in program, argument, environment key or value | `Error::InvalidValue` from `preflight`. |
| `arguments.len() > MAX_ARGUMENTS` | `Error::InvalidValue`; nothing spawned. |
| Command line over `MAX_COMMAND_BYTES` | `Error::InvalidValue`; nothing spawned. |
| `environment.len() > MAX_ENVIRONMENT_ENTRIES`, or a key that is empty or contains `=` | `Error::InvalidValue`; nothing spawned. |
| Either stream over `MAX_CAPTURED_BYTES` in `capture` | `Error::InvalidValue`; the child is killed and reaped first. |
| Child killed by a signal | `ReturnState::Crash`, `terminating_signal` set, `exit_code` `None`. |
| `Invocation::timeout` elapsed | child killed, `timed_out` true, message names the overrun. |
| Non-UTF-8 output | replaced lossily; never a panic and never an error. |

Evidence is **tier 3** (transcribed class-test literals) for the outcome/message
pairing, the `[EXTRA]` section's two markers and the `error_msg` emptiness
assertions, and **tier 4** for everything derived here: the Windows-predicate
simplification (a closed-form argument over `i32`, not a transcription), the
crash-in-every-mode classification, the ceilings, and the timeout. No C++ was
executed; the header has no retained output fixture, so tier 1 is not available
for this group.

Tests: `tests/system_process.rs` and `src/system/external_process.rs::tests`.
