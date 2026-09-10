# Log streams

The native `concept::log_stream` module implements the retained OpenMS `LogStream`,
`LogStreamBuf`, sink guard and notifier behavior through owned writers and callbacks.
The [reference tests](../tests/log_stream.rs) use four unchanged source output files
and independent boundary cases. [Provenance](../tests/data/log_stream_provenance.json)
records the exact source revision, hashes and adaptations. No C++ program was built
or executed to generate expectations.

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

The Windows C++ Colorizer also enables virtual-terminal console mode and resets
console state at process exit. This module emits the same per-record ANSI bytes
but does not alter host console mode; a Windows terminal must already support
those sequences. General mutable Colorizer composition and its separate insertion
API are outside this logging module. Local-time zone interpretation follows the
host's chrono-supported zone facilities. The process's global route objects live
for its lifetime; call `finish` explicitly for any directly buffered global text.
Thread-local and ordinary owned logger destruction flush their buffers normally.
