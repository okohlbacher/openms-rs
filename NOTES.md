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

Next step:
1. (done) Read Colorizer.cpp @ pin for Windows VT code.
2. Read src/concept/log_stream.rs (970 lines): find global streams, "finish", buffering.
3. Find the TOPP lifecycle's explicit global finish in src/cli (read only) and all exit paths.
