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

Next step:
1. Read Colorizer.cpp @ pin for Windows VT code.
2. Read src/concept/log_stream.rs (970 lines): find global streams, "finish", buffering.
3. Find the TOPP lifecycle's explicit global finish in src/cli (read only) and all exit paths.
