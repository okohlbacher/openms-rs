# Log streams

The native `concept::log_stream` module implements the retained OpenMS `LogStream`,
`LogStreamBuf`, sink guard and notifier behavior through owned writers and callbacks.
The [reference tests](../tests/log_stream.rs) use four unchanged source output files
and independent boundary cases. [Provenance](../tests/data/log_stream_provenance.json)
records the exact source revision, hashes and adaptations. One C++ program was built
and executed, for the process-exit behaviour only (see *Process exit* below); every
other expectation comes from the source files or from independent boundary cases.

## Public use

`LogStream::new(level)` creates a bound, initially unrouted logger; `Default` is the
source's unbound stream. It implements `std::io::Write`, including formatted writes.
`LogSink::new(writer)` owns a `Write + Send` destination. Cloning a sink shares its
identity and writer. `stdout()` and `stderr()` return stable process identities.
Repeated `insert` calls for the same identity are harmless.

The logger provides `level`/`set_level`, `insert`, `remove`, `remove_all_streams`,
`has_stream`, `insert_notification`, `set_prefix`, `set_all_prefixes`, `set_color`,
`set_clock`, `clear_cache`, `flush_incomplete`, `finish`, and `clone_configuration`.
A notification is an owned `Arc` callback receiving each complete rendered record,
including its newline, after output locks are released. A captured shared writer
can retain the accumulated text provided by the source's stringstream notifier.
Removing its sink unregisters the callback; ownership prevents dangling references.
`LogSinkGuard` temporarily removes a destination and provides mutable access to the
logger. Explicit `restore` exposes flush errors; scope exit also restores it.

`with_global_log` configures one of five `LogLevel` routes. The thread-local accessor
copies that route configuration on first use; each thread then owns its buffers and
cache. Later global changes do not change existing thread-local loggers. Concurrent
configuration or recursive access to a borrowed logger returns `WouldBlock`.
The six `openms_log_*` macros, `log_message` and `log_message_at` write and flush one
line. Raw `Write` supports source-style incremental expressions. Fatal logs include
the supplied source path and line but do not terminate the process. Debug locations
use a basename. Debug has no destination by default; info uses stdout, and fatal,
error and warning use stderr.

## API mapping

Every public member of `LogStream.h` at core `bc9cc12`, and the process-level
behaviour of `LogStream.cpp` and `Colorizer.cpp` it relies on.

| Source | Rust |
|---|---|
| `LogStreamBuf::MAX_TIME` | Not ported: nothing in the source reads it. |
| `LogStreamBuf::UNKNOWN_LOG_LEVEL` | `UNKNOWN_LOG_LEVEL` |
| `LogStreamBuf(log_level, Colorizer*)` | `LogStream::new(level)` and `LogStream::set_color` |
| `LogStreamBuf(source_buf, Colorizer*)` (the thread-local copy) | `LogStream::clone_configuration` |
| `~LogStreamBuf` | `Drop for LogStream`; `LogStream::finish` returns its errors |
| `sync` | `Write::flush` |
| `overflow` | The put buffer of `LOG_BUFFER_BYTES` synchronises when full, inside `Write::write` |
| `LogStreamBuf::setLevel`, `getLevel` | Declared without a definition in the source; `LogStream::set_level`, `level` |
| `LogStreamBuf::StreamStruct` (`stream`, `prefix`, `target`) | Private route: a `LogSink`, its prefix, an optional `LogNotification` |
| `clearCache` | `LogStream::clear_cache` |
| Protected `distribute_`, `expandPrefix_`, `getStreamList_`, `syncLF_`, `isInCache_`, `addToCache_`, `getNextLogCounter_`, `LogCacheStruct`, the buffer, level, stream list, incomplete line, colorizer and cache fields | Private implementation; no subclass hook |
| `LogStreamNotifier` (`logNotify`, `registerAt`, `unregister`, `stream_`, `registered_at_`) | A `LogNotification` callback given to `insert_notification`; `remove` unregisters it; a captured writer keeps the text `stream_` accumulates |
| `LogStream(buf, delete_buf, stream)` | `LogStream::new` then `insert`; `LogStream::default()` is the unbound stream |
| `~LogStream` | `Drop for LogStream` |
| `rdbuf`, `operator->` | Not ported: the buffer is not a separate object; its operations are methods of `LogStream` |
| `setLevel`, `getLevel` | `LogStream::set_level`, `LogStream::level` |
| `insert`, `remove`, `removeAllStreams`, `hasStream` | `insert`, `remove`, `remove_all_streams`, `has_stream` |
| `insertNotification` | `insert_notification` |
| `setPrefix(stream, prefix)`, `setPrefix(prefix)` | `set_prefix`, `set_all_prefixes` |
| `flush`, `flushIncomplete` | `Write::flush`, `flush_incomplete` |
| Private `findStream_`, `hasStream_`, `bound_`, `delete_buffer_` | Private implementation |
| `LogSinkGuard` (constructor, destructor; copy and move deleted) | `LogSinkGuard::new`, `Drop`, `restore`; it borrows the stream, so it cannot be copied or moved away from it |
| `getThreadLocalLogFatal` ... `getThreadLocalLogDebug` | `with_thread_local_log(LogLevel::Fatal ... Debug, ...)` |
| `getGlobalLogFatal` ... `getGlobalLogDebug` | `with_global_log(LogLevel::Fatal ... Debug, ...)` |
| `OPENMS_LOG_FATAL_ERROR`, `OPENMS_LOG_ERROR`, `OPENMS_LOG_WARN`, `OPENMS_LOG_INFO`, `OPENMS_LOG_DEBUG`, `OPENMS_LOG_DEBUG_NOFILE` | `openms_log_fatal_error!`, `openms_log_error!`, `openms_log_warn!`, `openms_log_info!`, `openms_log_debug!`, `openms_log_debug_nofile!` (`log_message`, `log_message_at`) |
| Destructors of the thread-local streams at thread end, and of the global streams at process exit (`LogStream.cpp:667-672`, `:694-722`) | `Drop` of the thread's streams at thread end; `finish_process_logs` and `ProcessLogs` for the process (*Process exit* below) |
| `InitConsole` and `Colorizer::isTTY` (`Colorizer.cpp:26-81`, `:155-170`), which the coloured streams rely on | `enable_console_colors`; the console test inside `LogSink` (*Console colours on Windows* below) |

## Preserved behavior

- The put buffer is 32,768 bytes. Writing an LF does not itself flush a smaller
  buffer. A full buffer or explicit `flush` processes complete lines and retains
  the final partial line. Only LF terminates lines; CR and arbitrary payload bytes
  remain unchanged, including UTF-8 split between writes.
- The two-entry cache suppresses repeated nonempty lines and refreshes their last
  use. Eviction emits `<message> occurred N times`, counting the original line.
  Clearing the cache emits summaries in bytewise lexical key order. Empty lines
  are never cached. Explicit partial-line flush bypasses the cache.
- `remove` flushes complete lines but retains partial text and cache. With no
  destinations, flushed put-buffer text is discarded while an earlier partial
  line remains. `remove_all_streams` flushes partial text and sinks but preserves
  cache contents. Reattaching a sink can therefore deliver a summary whose
  original line was sent elsewhere. `finish` emits complete lines, cache
  summaries, then the final partial line, matching buffer destruction.
- Prefixes belong to existing routes only. `%y` is the arbitrary current level;
  `%T`, `%t`, `%D`, `%d`, `%S`, `%s` have the source's local calendar formats;
  `%%` escapes a percent. Unknown escapes and a trailing percent are dropped.
  Summaries use the current prefix/level, not those present at first occurrence.
  Default local time uses pinned `chrono 0.4.45` with `clock`, without its default
  features. An injected `LogClock` provides deterministic calendar fields.
- All nine source single colors/styles and their corresponding reset sequences
  are supported. Colors are suppressed for non-terminal stdout/stderr, checked
  for every emission. Arbitrary owned writers receive the exact ANSI sequences.
- Sink guards flush before removal and before restoration, preventing partial
  text from leaking across the guarded interval. Restoration inserts a fresh
  route last, losing its old prefix and notifier exactly as in the source.
- One output mutex serializes complete records across logger instances and
  threads. Notifications run afterward and can log through another logger.

The source header's prose mentions numeric level filtering, temporary priorities
and stored message history, but its actual public implementation supplies none of
those operations. The native module does not invent them. C++ `streambuf` pointer
ownership, `rdbuf`/arrow access, raw notifier pointers and the unused `MAX_TIME`
constant are replaced by safe `Write`, owned configuration and callbacks. The
header declares buffer-level level setters/getters without corresponding source
implementations; the functioning `LogStream` operations are supported.

## Checked boundaries and platform differences

Levels and unexpanded prefixes are limited to 65,536 bytes, logical lines to 1 MiB,
and attached sinks to 1,024. Aggregate stored prefixes are limited to 16 MiB.
Each write, formatted write, flush or finalization has a shared conservative
16 MiB allocation/input budget and 50 million work-unit budget. Expansion and
multi-sink amplification are charged before allocation or output. A formatted
write shares one budget across all of its internal small writes. Arbitrary user
`Display`, writer, clock and callback implementations control their own work and
storage. Writers must not recursively log while their `write`/`flush` executes;
use the notification callback for logging after output instead.

I/O and callback errors are returned as `io::Error`; earlier external writes
cannot be rolled back. A failed output/formatted operation disables further
output to prevent destructor replay. `finish` and guard `restore` expose errors
that destruction can only ignore. Counter overflow and invalid injected calendar
fields are checked. Calendar years are restricted to 0000–9999; leap seconds up
to 60 can be supplied by a clock.

Whether an injected day exists is decided by `chrono::NaiveDate::from_ymd_opt`,
which uses the proleptic Gregorian calendar, so year 0000 is a leap year. The
time fields are checked in this crate, because
`chrono::NaiveTime::from_hms_opt` refuses the leap second 60. The chrono call
replaced a hand-written month-length table on 2026-09-13 without changing which
fields are accepted. The crate survey found no difference over years 0–9999,
months 0–13 and days 0–32. The unit test
`injected_calendar_days_follow_the_proleptic_gregorian_calendar` checks the same
against an independent leap-year rule. It tries every month and day in 14
years, including 10000 and 65535, which the year ceiling refuses, plus the
February and 30/31-day boundaries of every year through 10000. It passed before
and after the change. The check is integer arithmetic, so results are the same
on every machine.

General mutable Colorizer composition and its separate insertion API are outside
this logging module. Local-time zone interpretation follows the host's
chrono-supported zone facilities.

## Process exit

At process exit the source's destructors deliver what the streams still hold
(`~LogStreamBuf`, `LogStream.cpp:111-128`): complete buffered lines, then the
repetition summaries of the cache, then an unterminated line. `exit()` first
destroys the calling thread's `thread_local` streams (`LogStream.cpp:694-722`),
most recently created first, and then the static global streams in reverse
order of definition (`LogStream.cpp:667-672`): debug, info, warning, error,
fatal. Another thread's streams are destroyed when that thread ends.

Rust drops thread-local values when their thread ends, in the same order here,
but never drops a static. `finish_process_logs` therefore does what the
destructors do, in the source's order: the calling thread's streams, most
recently created first, then the five global streams. Each stays usable
afterwards. A stream whose output failed earlier is skipped; the others are
finished even when one fails, and the first error is returned. `ProcessLogs` is
a guard that calls it when dropped, at the end of a scope, on an early return,
and while a panic unwinds; `ProcessLogs::finish` returns the error instead.

Every tool executable holds a `ProcessLogs` for the whole run
(`run_from_environment` in `src/cli.rs`): whatever exit code the tool returns,
and also when it panics, the streams are finished after the tool's own output
and before the terminal reset. The crate has no `std::process::exit` call.

A program of its own needs neither. When a thread ends, its streams are
finished, most recently created first, and then, if the thread used
`with_global_log`, the global streams: the thread that configured or wrote to
them owns their exit-time flush. For the main thread that happens at process
exit, because Rust runs the main thread's thread-local destructors there, for a
return from `main` and for `std::process::exit` alike; measured on macOS and on
Linux with glibc (dax). A platform that does not run them (the Rust
documentation names Unix targets whose thread-local storage uses pthread keys)
loses what the main thread's and the global streams still hold unless the
program holds a `ProcessLogs` or calls `finish_process_logs`. A worker thread
that uses `with_global_log`, which the source documents as configuration to do
before threads start, finishes the global streams when it ends, earlier than
the source would.

The expected bytes of `process_exit_paths_deliver_what_the_release_build_delivers`
come from the C++ Release build at the pins (`../oracle/logstream-exit`,
`driver.cpp`, run twice on ibminode06, Linux x86_64): pending thread-local
state, pending global state, both at once, the order of two thread-local
streams, and a thread that ends before `main` logs again, each ended by a
return from `main` and by `std::exit`, which print the same bytes. The test runs
each case in a child process that ends by a return, an error return, a panic and
`std::process::exit`, and once without the guard, where the ending thread
finishes the streams.

One divergence: the source constructs its global debug stream with `delete_buf =
false` (`LogStream.cpp:672`, with the comment that `true` "segfaults"), so the
stream's buffer is never destroyed and what it holds is lost at exit. The port
delivers it, first among the global streams (CPP-355).

## Console colours on Windows

When libOpenMS is loaded, the source's `InitConsole` (`Colorizer.cpp:26-81`) adds
`ENABLE_VIRTUAL_TERMINAL_PROCESSING` to the console mode of standard output and,
independently, of standard error, ignoring a stream that is not a console, so
that the ANSI colour codes render; it does not restore the mode at exit.
`enable_console_colors` does the same through `winapi-util`
(`console::mode`/`console::set_mode`), once per process: the first time a record
is coloured for standard output or standard error, and when a tool executable
starts. On other platforms it does nothing. The console test behind colouring,
the source's `Colorizer::isTTY` (`Colorizer.cpp:155-170`), is `isatty` on Unix and, on Windows, whether
`GetConsoleMode` succeeds, which is false for the pipe a MSYS or Cygwin terminal
provides; the port asks the same questions.

No test executes the Windows branch. The GitHub Windows runners give a test
process pipes, not a console, so `GetConsoleMode` fails there and neither the
mode change nor colouring can be observed; this machine has no Windows console
either. The branch is compile-checked for `x86_64-pc-windows-msvc` (Rust 1.96);
the `schema-minimum-platforms` workflow job builds it with Rust 1.85 on
`windows-latest` when it runs.
