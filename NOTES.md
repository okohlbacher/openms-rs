# Lane notes: close CONCEPT/LogStream.h (do not ship; delete in final commit)

Pins: core = ../OpenMS4-tests/packages/core @ bc9cc12 (abs: /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-tests/packages/core)

Ledger (docs/core-sdk-reviewed-apis.json, lead-owned; report via .planning/ledger-fragments/<lane>.json):
- LogStream.h row ~l.1235: partial. Remaining: automatic Windows console VT mode; process-exit
  flushing of directly buffered global streams. Rust: src/concept/log_stream.rs; tests/log_stream.rs;
  docs/LOG_STREAM_SUPPORT.md
- Colorizer.h row ~l.1174: SEPARATE row, partial (class Colorizer object model, 9 globals, undoAll \033[0m,
  IndentedStream hook, Windows VT mode + exit reset not done).

Source facts (LogStream.cpp @ bc9cc12):
- ~LogStreamBuf (l.111-127): syncLF_(); clearCache(); distribute incomplete_line_; runs for the five
  file-static globals g_log_* (l.666-672) at static teardown AND for thread_local TLS streams
  (getThreadLocalLog*, l.690-718) at thread exit (main thread TLS at exit).
- logSinkMutex_ leaked on purpose so teardown can lock (l.49-53).

Findings (explored):
- Windows VT: Colorizer.cpp @bc9cc12 l.~25-75: static InitConsole windows_console_prep_and_restore;
  ctor (WIN only) GetStdHandle(STD_OUTPUT/ERROR) -> GetConsoleMode -> |= ENABLE_VIRTUAL_TERMINAL_PROCESSING
  -> SetConsoleMode, errors ignored; dtor (all OS) writes undoAll \033[0m to cout/cerr if isTTY.
  isattyWin = GetConsoleMode succeeds. -> crate candidates: anstyle-query::windows::enable_ansi_colors
  (check Cargo.lock), enable-ansi-support.
- Exit reset \033[0m already in src/cli/console.rs with_process_streams/reset_terminals (other lane).
- Rust globals: log_stream.rs `fn global()` static OnceLock<RwLock<LogStream>> x5 -> NEVER dropped.
  TLS: thread_local LOCAL [RefCell<Option<LogStream>>;5]; Drop for LogStream calls finish()
  (sync, clear_cache, partial). Main-thread TLS drop at exit is platform dependent (glibc/macOS run it
  via exit(); musl/pthread-key TLS do not).
- NO lifecycle call finishes the log streams: src/cli.rs:2322 `session.log.finish()` is ToolLog (-log file).
  Process path: src/cli.rs run() -> run_from_environment() (~l.2455) -> console::with_process_streams(||
  run_with(..)). src/bin/*.rs main returns ExitCode (no process::exit). No catch_unwind in cli; no
  panic=abort in Cargo.toml.
- Library users of log macros: format/mzml_header/read.rs, kernel/im_data_converter.rs,
  processing/peak_picking/noise.rs, cli/tools/feature_finder_centroided.rs.

Plan:
- log_stream.rs: pub fn finish_process_logs() -> io::Result<()>: calling thread's TLS streams in reverse
  first-access order (C++ thread_local dtor order), then the 5 globals in reverse construction order
  (debug, info, warn, error, fatal). Skip already-failed streams; attempt all; return first error.
- pub struct ProcessLogs guard (must_use) whose Drop calls it (library mains; unwinding).
- src/cli.rs run_from_environment: call inside with_process_streams closure after run_with (TLS summary
  precedes the \033[0m reset as in C++: thread_local dtors before statics).
- Tests: self re-exec child process per exit path (return, process::exit, panic, error exit code).

Windows VT crate decision (explored):
- anstyle-query 1.1.5 enable_virtual_terminal_processing: stdout then stderr with `?` -> skips stderr when
  stdout is redirected (C++ does each handle independently) -> REJECT.
- winapi-util 0.1.11 (BurntSushi, Unlicense/MIT, edition 2021, windows-sys >=0.48,<=0.61, cached in
  ~/.cargo/registry): Console::stdout()/stderr() + set_virtual_terminal_processing(true); also pub
  console::mode(HandleRef::stderr()) == exact isattyWin. -> USE (cfg(windows) target dep).
- Local target x86_64-pc-windows-msvc installed: `cargo check --target x86_64-pc-windows-msvc
  --no-default-features` can compile-check. CI windows job (cross-platform) runs only on tags/dispatch,
  and GH runners have pipes, not a console -> VT mode not observable in CI.

More source facts:
- g_log_debug is built with delete_buf=false (LogStream.cpp l.672 "should be 'true', but segfaults") ->
  its buffer is leaked, never flushed at exit. Port flushes it anyway (deliberate divergence, document).
- C++ exit order: main-thread thread_locals first (reverse construction = reverse first access), then
  statics in reverse definition order: debug(leaked), info, warn, error, fatal.
- clearCache prints "<line> occurred N times" in std::map key order (Rust already sorts; ok).
- Rust LOCAL is an array -> drops Fatal..Debug in index order, NOT reverse first-access -> fix.
- Macros (C++) do not flush; Rust log_message flushes per line (pending empty); summaries remain.

Decisions:
- No automatic main-thread hook (thread-name heuristic rejected). Explicit finish_process_logs() +
  #[must_use] guard (Drop = best effort, covers unwinding) + lifecycle call in src/cli.rs
  run_from_environment (not src/cli/**). Library main without guard: TLS finished by TLS dtor where the
  platform runs it; direct global writes need the guard -> document.
- VT: pub fn in log_stream.rs, OnceLock-once, called lazily from LogSink colour path for Stdout/Stderr
  and at start of run_from_environment. isTTY mintty divergence (IsTerminal vs GetConsoleMode) is a
  Colorizer.h member -> note only.

Work order (commit after each):
A. log_stream.rs: TLS container w/ order + Drop; finish_process_logs; guard.
B. tests/log_stream.rs: order + child re-exec per exit path (return, process::exit, panic, error code).
C. src/cli.rs run_from_environment: guard.
D. VT: Cargo.toml [target.'cfg(windows)'.dependencies] winapi-util =0.1.11; Cargo.lock; decisions doc;
   cargo check --target x86_64-pc-windows-msvc --no-default-features.
E. docs/LOG_STREAM_SUPPORT.md. F. gates, fragment .planning/ledger-fragments/<lane>.json, sweep.

Oracle: ../oracle/logstream-exit/{driver.cpp,remote_run.sh,run.sh} (written; run.sh launched in bg,
  results -> ../oracle/logstream-exit/results/<case>/{stdout,stderr,combined}.bin + exit, FETCHED marker).
  Cases: tls_return tls_exit tls_order global_return global_exit both thread. If not fetched, rerun run.sh.
- FYI (other lane): src/cli/tools/feature_finder_centroided.rs LogStreamLines duplicates the LogStream
  cache incl. an admitted exit-order divergence -> report as follow-up, do not touch.

ORACLE DONE (ibminode06, Linux Release build at pins; run1==run2, differ.txt empty). Combined (2>&1):
 tls_return/tls_exit(exit 3): 'same\nerr\ninfo partial\n<err> occurred 2 times\n<same> occurred 3 times\n'
 tls_order: 'e\nw\n<w> occurred 2 times\n<e> occurred 2 times\n'
 global_return/global_exit: 'global line\nglobal partial\ngw\n<gw> occurred 2 times\nge\ngf\n' (debug 'gd' LOST)
 both: 'w\ninfo partial\n<w> occurred 2 times\ngw\n<gw> occurred 2 times\n'
 thread: 'te\ntw\n<tw> occurred 2 times\ntw partial\n<te> occurred 2 times\nafter join\n'
 -> confirms: TLS reverse-creation order, then globals info,warn,error,fatal; per stream sync,cache,partial.
 Decision: port flushes the global debug stream too = deliberate divergence, provisional CPP-355.

DONE: A (c80f5ef), B test passes 18/18 (--no-default-features --test log_stream), mutation (forward
 order) fails it; C src/cli.rs run_from_environment guard; D Cargo.toml winapi-util =0.1.11 (lock unified
 with windows-sys 0.52.0), enable_console_colors + is_console (cfg windows), cargo check --target
 x86_64-pc-windows-msvc --no-default-features ok (only pre-existing src/system/file.rs:1493 unused mut warning).
 MSRV 1.85 on windows NOT measured (1.85 toolchain lacks windows target; not downloading).
E done: docs/LOG_STREAM_SUPPORT.md Process exit + Windows sections, THIRD_PARTY row; citations fixed
 (LogStream.cpp:111-128, 667-672, 694-722; Colorizer.cpp:26-81, 87-110, 139-153, 155-170).
 Python gates pass (core_sdk, citations, cycles; doc coverage improved, floor not yet written).
Next step (for complete): 1 rustdoc 100% for log_stream.rs (was 47.8%, 22/46) + check_doc_coverage --write;
 2 API mapping table (every LogStream.h member) in LOG_STREAM_SUPPORT.md; 3 tests/data/log_stream_provenance.json
 add oracle hashes (../oracle/logstream-exit/results/hashes.txt); 4 CHANGELOG/VALIDATION if convention;
 5 clippy all-features + sweep on dax; 6 fragment .planning/ledger-fragments/<lane>.json.
