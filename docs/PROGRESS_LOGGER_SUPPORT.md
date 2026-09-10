# ProgressLogger

The native module `concept::progress_logger` implements the pinned Core SDK
`CONCEPT/ProgressLogger.h` and `.cpp` dispatch and command-output behavior. It is
independent of a GUI toolkit and logging framework. Its only direct dependency
is the pinned `cpu-time =1.0.0` platform CPU-clock adapter on Unix and Windows.

## API mapping

| Source API | Native API |
| --- | --- |
| `ProgressLogger()`, destruction | `ProgressLogger::new()` / `Default`; ordinary drop, with no implicit end |
| copy constructor, assignment | `Clone`, `clone_from`; fresh backend, copied type and throttle timestamp |
| `LogType::{CMD,GUI,NONE}` | `ProgressLogType::{Cmd,Gui,None}`, with the same discriminants 0/1/2 |
| `setLogType`, `getLogType` | `set_log_type`, `log_type` |
| `setLogger` | `set_logger(Box<dyn ProgressBackend>)`, transferring ownership |
| `startProgress`, `setProgress`, `nextProgress`, `endProgress` | `start_progress`, `set_progress`, `next_progress`, `end_progress(bytes_processed)` |
| four virtual `ProgressLoggerImpl` operations | four corresponding methods of `ProgressBackend`, returning `Result` |
| `make_gui_progress_logger` | per-logger `set_gui_factory`, retained by copies; GUI defaults to a no-op |
| source process-static recursion depth | shared `ProgressNesting::global()` by default; `Default` creates an isolated nesting context |

`CommandProgressLogger<W: std::io::Write>` is also directly usable as a backend.
`new(writer)` uses the system clock. `with_clock(writer, clock)` permits supplied
wall and CPU samples; `into_inner()` returns the writer. Use the
`ProgressBackend` trait to call its four methods. `ProgressLogger` uses stdout
when selecting `Cmd`; replacing its backend permits other writers or callbacks.

```rust
use openms::concept::progress_logger::{ProgressLogType, ProgressLogger};

let mut logger = ProgressLogger::new();
logger.set_log_type(ProgressLogType::Cmd);
logger.start_progress(0, 3, "processing")?;
for _ in 0..3 {
    logger.next_progress()?;
}
logger.end_progress(0)?;
# Ok::<(), openms::Error>(())
```

## Preserved behavior

The default is `None`. Selecting a log type always creates a new backend, even
when that type was already selected. Backend replacement does not change the
reported type. Copying a logger discards its active/custom backend, constructing
a fresh backend for the copied type. Backend counter/timer state is not copied.
Changing the GUI factory does not replace an already selected GUI backend.

Every successful start dispatches at the current nesting depth and increments
that depth afterward, including starts with logging disabled. Every end first
decrements the depth if nonzero, then dispatches, even without a matching start.
Dropping a logger does not end progress or unwind nesting. A second start on an active command backend prints its new header, replaces
the range/counter and resets the timer, then errors: the current source
`StopWatch::reset()` restarts a running timer, so the following `start()` throws.
The failed wrapper call does not increment nesting. The caller balances
successful starts and ends.

Start records the current whole wall-clock second. Set suppresses output if that
second equals the last recorded second. This compares calendar-second buckets,
not elapsed durations: moving the wall clock backward to a different second
also permits output. An unsuppressed set reads the clock twice, matching the
source comparison and assignment. Next increments the backend counter before
this throttle. Explicit set never repositions the command backend's next counter.
The no-op backend's next operation returns zero.

The command backend retains the literal labels, two-space nesting indentation,
carriage returns, dots for an equal begin/end range, and out-of-range diagnostic.
For valid ranges it calculates `(value - begin) as f32 / (end - begin) as f32 *
100_f32` and formats two fractional digits. An invalid value is printed as a
diagnostic, not returned as a range error. Equal endpoints print a dot regardless
of the supplied value. Source output is rendered in the classic numeric locale.

Duration formatting follows the source second/minute/hour/day branches; branch
selection truncates elapsed seconds first. Thus 59.999 seconds renders as
`60.00 s`, whereas 60 seconds renders as `01:00 m`. Throughput first truncates the
floating bytes-per-second quotient to `u64`, then selects binary units and four
significant digits. The source fallback text above the PiB range is retained.

## Native timing and checked boundaries

`ProgressClock` is a shared callable returning `ProgressTime`: an independent
`wall_second` for throttling, absolute `wall_seconds` for the timer, and optional
absolute `cpu_seconds`. Only differences of the latter two are displayed.
`with_clock_and_nesting` provides deterministic timing and independent nesting.
The standard clock uses `SystemTime` seconds for the throttle and `Instant`
elapsed time for duration. On Unix and Windows it samples total process CPU time
through the safe [`ProcessTime::try_now` API](https://docs.rs/cpu-time/1.0.0/cpu_time/struct.ProcessTime.html)
and propagates operating-system errors. CPU sums the process's threads rather
than measuring only the calling thread. Other platforms explicitly report CPU
as unavailable; the command summary then says `unavailable (CPU)`. Injected
clocks can also omit CPU samples. CPU is never inferred from wall time or
represented as a false zero.

The GUI factory is owned per logger instead of being a mutable global function
pointer. Standard clock, factory, and nesting handles are shared by clones.
Backends are `Send`, permitting `Arc<Mutex<ProgressLogger>>` for shared
parallel progress. Methods take mutable Rust references instead of mutating
through C++ const methods. Shared nesting uses checked atomic operations, avoiding source data
races; related operations should still be serialized when exact cross-logger
ordering matters. Isolated contexts avoid coupling independent jobs.

Labels are bounded by `MAX_PROGRESS_LABEL_BYTES` (1 MiB) and nesting/indentation
by `MAX_PROGRESS_DEPTH` (1024). These checks apply before output allocation or
backend dispatch. There is no cumulative lifetime operation/output cap; each
built-in operation uses bounded scratch and fixed-size state. The caller owns
writer buffering, and custom clocks/backends/factories are caller code outside
these internal work bounds.

Inverted ranges, counter/difference overflow, invalid timer values, elapsed
intervals above `i32::MAX` seconds, and nonfinite/out-of-`u64` throughput return
errors. Zero elapsed time is valid without throughput but errors when a nonzero
byte count requests a rate. Negative or nonfinite elapsed samples are not
silently clamped. The command timer must have started before end. GUI/None end
remains a no-op even without start. Checks avoid undefined source signed
arithmetic and floating-to-integer conversions.

Clock/backend/writer errors propagate. Output streams may already contain a
partial write when an I/O failure occurs. A command start updates its range/counter before header output and samples
timing after the header flush, excluding slow header I/O. A late start failure
can therefore follow observable state/output changes. A rejected command end
remains retryable. Wrapper throttle
assignment happens before backend set/start, and wrapper end decrements nesting
before dispatch, retaining the source operation order. Failed start dispatch
does not increment nesting. These operations do not pretend to roll back user
callback side effects or partially written bytes.

## Evidence

`tests/progress_logger.rs` includes the source mode/copy/NONE smoke cases plus
independent deterministic clock/backend traces and literal command-output
oracles. It checks same-second and backward-second updates, two clock reads,
copy timestamps, next-counter behavior, shared and bounded nesting, GUI
replacement, duration/rate boundaries, CPU-unavailable output, and checked I/O,
clock, overflow, and zero-duration failures. Output literals are derived from
source expressions; they are not claimed as captured C++ fixture output.

`tests/data/progress_logger_provenance.json` pins the original source/header/test
and formatting dependencies. No C++ code was built or executed to generate the
native tests. A live-clock sanity test checks finite nondecreasing total process
CPU samples without requiring a minimum elapsed time or workload.
The main-crate checks are recorded in [validation results](VALIDATION.md).
