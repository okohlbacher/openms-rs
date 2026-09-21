# ProgressLogger

The native module `concept::progress_logger` implements the pinned Core SDK
`CONCEPT/ProgressLogger.h` and `.cpp` dispatch and command-output behavior. It is
independent of a GUI toolkit and logging framework. Its only direct dependency
is the pinned `cpu-time =1.0.0` platform CPU-clock adapter on Unix and Windows.

The reference is the Linux x86_64 **Release** build (a project decision). That
build compiles `OPENMS_PRECONDITION`/`OPENMS_POSTCONDITION` out, so a check that
exists only in a Debug build is not source behavior for this port.

## API mapping

| Source API | Native API |
| --- | --- |
| `ProgressLogger()`, destruction | `ProgressLogger::new()` / `Default`; ordinary drop, with no implicit end |
| copy constructor, assignment | `Clone`, `clone_from`; fresh backend, copied type and throttle timestamp |
| `LogType::{CMD,GUI,NONE}` | `ProgressLogType::{Cmd,Gui,None}`, with the same discriminants 0/1/2 |
| `setLogType`, `getLogType` | `set_log_type`, `log_type` |
| `setLogger` | `set_logger(Box<dyn ProgressBackend>)`, transferring ownership |
| `startProgress`, `setProgress`, `nextProgress`, `endProgress` | `start_progress`, `set_progress`, `next_progress`, `end_progress(bytes_processed)` |
| `startProgress`'s Debug-only `OPENMS_PRECONDITION(begin <= end)` | not ported: absent from the Release reference build, so any range is accepted |
| four virtual `ProgressLoggerImpl` operations | four corresponding methods of `ProgressBackend`, returning `Result` |
| `make_gui_progress_logger` | per-logger `set_gui_factory`, retained by copies; GUI defaults to a no-op |
| source process-static recursion depth | shared `ProgressNesting::global()` by default; `Default` creates an isolated nesting context |
| the `ProgressLogger` *base* of an algorithm class (`class GaussFilter : public ProgressLogger`) | the algorithm's `*_with_progress` entry point, which borrows a caller's `ProgressLogger` for the call; see [Consumers](#consumers) |
| — | `ProgressReporter`: one run's progress calls sent to an optional logger (`start`, `start_count`, `set`, `set_count`, `next_progress`, `end`, `end_with_bytes`; the two `_count` forms convert a record count only when reporting, so a silent run cannot fail on one), with `section`, which ends the section on failure too |
| — | `progress_value`: a `usize` count as the source's `SignedSize` value, refusing what would wrap |

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

Every range is accepted, including `begin > end`, and reaches the backend
unchanged. The source's only range check,
`OPENMS_PRECONDITION(begin <= end, "ProgressLogger::init : invalid range!")`
(`ProgressLogger.cpp:235`), is Debug-only; the Release build has none. Its
command backend stores the range as given (`:36-40`). A set then prints a dot
when begin equals end (`:48`), the diagnostic
`ProgressLogger: Invalid progress value '<v>'. Should be between '<begin>' and '<end>'!`
when the value lies below begin or above end (`:52-55`), and the percentage
otherwise (`:60`). For an inverted range every value, including both endpoints,
takes the diagnostic branch. The percentage is therefore only computed with
`begin < end` and never divides by zero. The no-op backend and the default GUI
factory ignore the range. This is the executed Release behavior (see
[Evidence](#evidence)).

Every successful start dispatches at the current nesting depth and increments
that depth afterward, including starts with logging disabled. Every end first
decrements the depth if nonzero, then dispatches, even without a matching start.
Dropping a logger does not end progress or unwind nesting. A second start on an active command backend prints its new header, replaces
the range/counter and resets the timer, then errors: the current source
`StopWatch::reset()` restarts a running timer, so the following `start()` throws
`Exception::Precondition` ("StopWatch is already started!", `StopWatch.cpp:43`).
That throw is unconditional, not a Debug-only macro, and the Release build
executes it. An end on a command backend that was never started likewise fails
in Release (`StopWatch::stop`, `StopWatch.cpp:55`) before printing anything.
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

Counter/difference overflow, invalid timer values, elapsed intervals above
`i32::MAX` seconds, and nonfinite/out-of-`u64` throughput return errors. Zero
elapsed time is valid without throughput but errors when a nonzero byte count
requests a rate. Negative or nonfinite elapsed samples are not silently clamped.
GUI/None end remains a no-op even without start. Checks avoid undefined source
signed arithmetic and floating-to-integer conversions. An inverted range is not
an error (see [Preserved behavior](#preserved-behavior)).

Clock/backend/writer errors propagate. Output streams may already contain a
partial write when an I/O failure occurs. A command start updates its range/counter before header output and samples
timing after the header flush, excluding slow header I/O. A late start failure
can therefore follow observable state/output changes. A rejected command end
remains retryable. Wrapper throttle
assignment happens before backend set/start, and wrapper end decrements nesting
before dispatch, retaining the source operation order. Failed start dispatch
does not increment nesting. These operations do not pretend to roll back user
callback side effects or partially written bytes.

### Native differences: every refusal and its source status

Each refusal in `src/concept/progress_logger.rs` was checked against the pinned
source for a Debug-only `OPENMS_PRECONDITION`/`OPENMS_POSTCONDITION` ported as a
Release refusal. The only one found, `begin > end`, was removed from the wrapper
and the command backend. `ProgressLogger.cpp:235` is the only use of either
macro in `ProgressLogger.{h,cpp}`, `StopWatch.{h,cpp}`, `SysInfo.{h,cpp}` and
the `StringUtils` formatting they call. Every remaining refusal is either a
native bound or a check the Release build also executes:

| Refusal | Where | Status |
| --- | --- | --- |
| label longer than `MAX_PROGRESS_LABEL_BYTES` | wrapper and command start | native bound: the source accepts any label |
| depth at `MAX_PROGRESS_DEPTH` on start, above it on a backend call | wrapper start; command start/set/end | native bound: source `static int recursion_depth_` is unbounded and overflows at `INT_MAX` |
| second start on a running command backend | command start | source, Release: `StopWatch::start` throws (`StopWatch.cpp:43`); captured |
| end without a running command timer | command end | source, Release: `StopWatch::stop` throws (`StopWatch.cpp:55`); captured |
| next-counter overflow | command next | native: source `++current_` is undefined signed overflow |
| `value - begin` or `end - begin` overflows `i64` | command set, percentage branch | native: undefined signed overflow in the source (`:60`) |
| nonfinite or negative clock sample | command start/end | native: injected clocks have no source counterpart |
| elapsed time above `i32::MAX` seconds | command end | native: source day arithmetic multiplies `int`s |
| nonfinite or out-of-`u64` throughput, including a byte count over zero elapsed time | command end | native: undefined floating-to-integer conversion in the source (`:78`) |
| clock, writer and custom-backend errors | all | native: `Result` propagation where the source reports no failure |
| system clock before the epoch or beyond `i64` seconds | `system_progress_clock` | native |

## Consumers

The header's own public surface is complete (the table above, with the
executed Release evidence below). What the ledger entries of its consumers
recorded as "the `ProgressLogger` base is not ported" was, for each of them,
the missing *use* of it: the source's algorithm calls `startProgress`,
`setProgress` and `endProgress` on itself, and the Rust algorithm made no call.

The source class inherits the logger, and its `const` members mutate it through
`mutable` fields. A Rust algorithm keeps its methods `&self` and borrows the
logger from its caller for the one call instead: a `*_with_progress` entry
point takes `&mut ProgressLogger`, and the caller chooses the type with
`set_log_type` (a TOPP tool passes `ToolContext::progress_log_type`, which is
`Cmd` unless `-no_progress`, `TOPPBase.cpp:400-403`) or a backend with
`set_logger`, exactly as a source caller does on the algorithm object. The
entry points without a logger report nothing, which is the source's default
type `NONE` without its cost. Every consumer below builds its calls with
`ProgressReporter`, so the silent and the reporting entry point share one body,
and a test (`progress_changes_no_result`) checks that they return the same
result.

`ProgressReporter::section` ends a section whether its body succeeds or fails.
The source does not: when an algorithm throws inside its section, `endProgress`
is never reached, the command output stops after the last percentage, the
process-static depth stays one level deeper for every later section, and the
object's next command `startProgress` throws `StopWatch is already started!`.
The executed Release run shows the first two for `PeakPickerHiRes` in manual
mode on a centroided spectrum and for `GaussFilter` with a ppm width on a
chromatogram. The port prints the `-- done` line instead and leaves the
logger reusable; that is the one difference in output on a path the source also
takes, and the replay asserts it against the captured C++.

| Consumer (ledger entry) | Rust entry point | Section |
|---|---|---|
| `PeakPickerHiRes::pickExperiment` | `PeakPickerHiRes::pick_experiment_with_progress`, `pick_experiment_in_place_with_progress` | `picking peaks`, spectra + chromatograms, `1..=n` |
| `PeakPickerIterative::pickExperiment` | `PeakPickerIterative::pick_experiment_with_progress` | `picking peaks`, spectra, `0..n` |
| `LinearResamplerAlign::rasterExperiment` | `LinearResamplerAlign::raster_experiment_with_progress` (and the newly ported `raster_experiment`) | `resampling of data`, spectra, `0..n` |
| `GaussFilter::filterExperiment` | `GaussFilter::filter_experiment_with_progress` | `smoothing data`, spectra + chromatograms, `1..=n` |
| `SavitzkyGolayFilter::filterExperiment` | `SavitzkyGolayFilter::filter_experiment_with_progress` | `smoothing data`, spectra + chromatograms, `1..=n` |
| `MorphologicalFilter::filterExperiment` | `MorphologicalFilter::filter_experiment_with_progress` | `filtering baseline`, spectra, `0..n` |
| `PeptideIndexing::run` | `PeptideIndexing::run_with_progress` | `Load first DB chunk` `0..1`, then `Aho-Corasick` over the proteins, `1..=n` |
| `SignalToNoiseEstimatorMedian::init` (earlier) | `SignalToNoiseEstimatorMedian::estimate_with_progress` | `noise estimation of data` |
| `FeatureFinderAlgorithmPicked::run` (earlier) | its `set_log_type` / `set_progress_logger` | as in `docs/FEATURE_FINDER_PICKED_SUPPORT.md` |

`CONCEPT/Macros.h` also names this row, but only as the precedent for never
turning a Debug-only `OPENMS_PRECONDITION` into a Rust refusal; it owes no
progress call.

## FORMAT readers

The file adapters derive from `ProgressLogger` too, and their handlers make
the calls. Each Rust reader gains `*_with_progress` entry points that take the
caller's logger; they run the same code as the silent entry points, whose calls
go to `ProgressReporter::silent()`,
so a result, a written byte or an error cannot differ between the two
(`progress_changes_no_result` and `progress_changes_no_error` in
`tests/progress_format_readers.rs` check it on every reader). `ImzMLFile` was
wired earlier (`docs/IMZML_FILE_SUPPORT.md`), and `FileInfo` reports once it
passes its log type to these loaders.

Three rules follow the Release build rather than `ProgressReporter::section`:

1. **A failed section stays open.** When a reader fails after a start, the
   source's exception bypasses `endProgress`: no `-- done` line, the static
   depth stays one level deeper, and the object's command backend refuses its
   next start with `StopWatch is already started!` (`StopWatch.cpp:43`). The
   readers make their calls one by one and end a section only on success, so
   all three hold here too; the replay's `mzml_reuse_after_failure` case loads
   a truncated mzML and then a good one through one logger, and both builds
   refuse the second load's list start in command mode.
2. **The calls reach the logger the source uses.** `DTA2DFile`,
   `MascotGenericFile`, `MzDataFile`, `MzXMLFile` and `MzMLFile` hand their
   handler the file object itself, so the calls go to the caller's logger.
   `ConsensusXMLFile` and `FeatureXMLFile` hand theirs only
   `setLogType(getLogType())`, and `MzMLHandler` reports the whole document on
   `pg_outer`, a thread-local copy of the file's logger; those calls go to a
   copy of the caller's logger, made as `ProgressLogger::clone` makes one: a
   fresh backend of its type (its GUI factory for `Gui`), its nesting shared. A
   backend installed with `set_logger` therefore sees only the calls made on
   the logger itself, as a source file's `setLogger` backend does.
3. **A store opens its destination first.** `XMLFile::save_` opens the file
   before its handler's `writeTo` makes a call, so a destination that cannot be
   created makes none. The streaming writers already build inside the atomic
   publication; the ones that build a whole document in memory first
   (featureXML, consensusXML, mzData) build it inside the publication when
   reporting (`path_io::store_reporting`), and if the destination cannot be
   prepared they build it again without reporting, so a refusal still wins over
   the destination's error as it does for the silent writer.

| Source member | Rust entry point | Calls, as the Release build makes them | Logger |
|---|---|---|---|
| `DTA2DFile::load` | `dta2d::load_with_progress` | `startProgress(0, 0, "loading DTA2D file")` before the file opens, `setProgress(0)` as each spectrum begins, `endProgress()` (`DTA2DFile.h:72-248`) | the caller's |
| `DTA2DFile::store`, `storeTIC` | `dta2d::store_with_progress`, `store_tic_with_progress` | `startProgress(0, spectra, "storing DTA2D file")` before the file is created, `setProgress(i)` per spectrum (not for the TIC), `endProgress()` (`:258-287`, `:297-320`) | the caller's |
| `MS2File::load` | none needed | no call: the only one is commented out (`MS2File.h:52`) | — |
| `MascotGenericFile::load` | `mascot_generic::load_with_progress`, `MascotGenericFile::load_with_progress` | `startProgress(0, file size, "loading MGF")`, `setProgress(is.tellg())` per block, which is -1 after a last `END IONS` without a newline, `endProgress()` (`MascotGenericFile.h:74-104`); none for a missing file | the caller's |
| `MascotGenericFile::store` (both overloads) | `MascotGenericFile::store_with_progress`, `store_to_with_progress` | `startProgress(0, spectra, "storing mascot generic file")` after the header, `setProgress(i)`, `endProgress()` (`MascotGenericFile.cpp:458-476`) | the caller's |
| `MzIdentMLFile::load`, `store` | none needed | no call: the handler makes none | — |
| `ConsensusXMLFile::load` | `consensusxml::load_with_progress` | `startProgress(0, 0, "loading consensusXML file")`, then `setProgress(1)`, `setProgress(2)`, … for the root and every `map`, `consensusElement`, `IdentificationRun`, `ProteinHit`, `PeptideHit` and `dataProcessing` element, `endProgress()` (`ConsensusXMLHandler.cpp:130-133`, `:254-256`) | a copy |
| `ConsensusXMLFile::store` | `consensusxml::store_with_progress` | `startProgress(0, 0, "storing consensusXML file")`, `setProgress(1)` … `setProgress(5 + runs + column headers + features)`, `endProgress()` (`:606-837`) | a copy |
| `FeatureXMLFile::load` | `featurexml::load_with_progress` | `startProgress(0, count, "Loading featureXML file")`, `setProgress(features kept)` as each top-level feature begins, `endProgress()` (`FeatureXMLHandler.cpp:319`, `:1047`, `:838`) | a copy |
| `FeatureXMLFile::loadSize` | `featurexml::load_size` | no call: the load stops before the section (`:306-316`) | — |
| `FeatureXMLFile::store` | `featurexml::store_with_progress` | `startProgress(0, features, "Storing featureXML file")`, `setProgress(i)` after feature `i`, `endProgress()` (`:215-222`) | a copy |
| `MzDataFile::load` | `mzdata::load_with_progress`, `MzDataFile::load_with_progress` | `startProgress(0, count, "loading mzData file")` at `<spectrumList>` (`MzDataHandler.cpp:347-356`), a set per spectrum with a process-wide counter, incremented first (`:436`, `:452`), `endProgress()` at `</mzData>`, which resets the counter (`:460-464`) | the caller's |
| `MzDataFile::store` | `mzdata::store_with_progress`, `MzDataFile::store_with_progress` | `startProgress(0, spectra, "storing mzData file")`, `setProgress(s)`, `endProgress()` (`:579-1073`) | the caller's |
| `MzXMLFile::load` | `mzxml::load_with_progress`, `MzXMLFile::load_with_progress` | `startProgress(0, scanCount, "loading mzXML file")` at `<msRun>`, `setProgress(n)` as each scan begins, nested or filtered, `endProgress()` at `</mzXML>` (`MzXMLHandler.cpp:130-136`, `:281-282`, `:524-531`) | the caller's |
| `MzXMLFile::store` | `mzxml::store_with_progress`, `MzXMLFile::store_with_progress` | `startProgress(0, spectra, "storing mzXML file")`, `setProgress(s)`, `endProgress()` (`:636`, `:864`, `:1119`) | the caller's |
| `MzMLFile::load` | `mzml::load_with_progress` | on a copy: `startProgress(0, 1, "loading mzML")` at `<mzML>` and `endProgress(file size)` at `</mzML>` (`MzMLHandler.cpp:1203`, `:1524`); on the caller's: `startProgress(0, count, "loading spectra list")` or `"loading chromatogram list"`, `nextProgress()` per record, `endProgress()` per list (`:966`, `:997`, `:1443`, `:1483`, `:1491`, `:1497`) | both |
| `MzMLFile::store` | `mzml::store_with_progress` | `startProgress(0, spectra + chromatograms, "storing mzML file")`, `setProgress(n)` per record, `endProgress(bytes written)` (`:4763-4838`) | the caller's |

Where the port's calls differ, the replay asserts the difference against the
captured calls:

- **consensusXML parses before it reports.** The reader parses the whole
  document before it converts any of it, so it makes its calls after a
  successful parse: a document that is not well-formed makes none, where the
  Release build made those for the elements before the defect
  (`consensus_load_truncated`), and one the conversion refuses has made every
  set and no end. featureXML hands each feature over as it is parsed, so only a
  failure inside a feature, or inside `<featureList>` before its first feature
  is complete, shows: the port makes a feature's set when the feature ends.
- **mzML decodes arrays as they close.** The source decodes a pool of spectra
  when it is flushed, by default at `</mzML>` (`MzMLHandler.cpp:1425-1428`,
  `:1522-1523`), so a document with an undecodable array fails after fewer
  calls here.
- **The mzML store counts its own bytes.** `endProgress(os.tellp())` reports
  the size of the document written; the port's document is not the source's,
  so the count is the port's file size. With a `.gz` or `.bz2` suffix it is the
  count before compression, where the source's compressing stream has no
  position and passes -1.
- **Native refusals come first where the port checks before it writes.** The
  source has no ceilings and refuses almost nothing; a refusal the port makes
  before a section starts makes no call.
- **Skipping chromatograms follows the corrected reader (CPP-017).** With
  `PeakFileOptions::setSkipChromatograms` the source handler ignores every
  element until `</chromatogramList>` (`MzMLHandler.cpp:149`, `:870-873`): the
  Release build starts no section, advances on every record end, ends two
  sections it never began and, in command mode, fails the load when the
  backend refuses that end (`StopWatch.cpp:55`). The port skips only the
  chromatograms, and its calls are those of the ordinary load
  (`mzml_load_skip_chromatograms`).

The mzData handler's scan counter is a function-local `static UInt`
(`MzDataHandler.cpp:436`) that only `</mzData>` resets, so after a load that
failed inside the spectrum list the next load's values continue past the
range. The port keeps the same process-wide counter, atomically, and every
mzData load advances it, reporting or not (`mzdata_static_counter`). The
source's other file entry points that report progress are not wired here and
report nothing: `MzMLFile::loadSize`, `loadBuffer`, `storeBuffer` and both
`transform` overloads, and `MzXMLFile`'s `transform` overloads.

## Evidence

**Executed Release differential (tier 1).** An oracle driver,
`../oracle/progress-logger-release-range/driver.cpp`, links the Linux x86_64
Release install `openms4-release-bc9cc12-c19e494-174b576` (core `bc9cc12`, cli
`c19e494`, topp `174b576`). It was built with that install's g++ 14.4.0 and
`-O2`, and run three times on ibminode06. The installed `config.h` leaves
`OPENMS_ASSERTIONS` undefined, and `libOpenMS.so` contains no copy of the
`invalid range` message. The driver makes 60 calls:

- `startProgress(5, 0)` with values below, on, between and above both ends, from
  `INT64_MIN` to `INT64_MAX`, plus two `nextProgress`;
- `startProgress(2, 1)` with values 0 to 3;
- `startProgress(0, 0)` with values 0 and 1;
- nested ranges on two loggers, including an inverted inner range;
- a second start on an active command logger;
- an end without start;
- NONE and the default GUI backend with an inverted range.

`setProgress` is throttled to one dispatch per wall-clock second
(`ProgressLogger.cpp:244`). The driver therefore waits before every set/next
until `time(nullptr)` differs from the logger's `last_invoke_`, and the op log
confirms that all 31 dispatched. Each call's stdout bytes are cut out by file
offset. Only the two timing texts of the summary line are masked. All three runs
exit 0 with empty stderr and are identical after masking. The masked table is
`tests/data/progress_logger_release_range.tsv`.
`release_range_fixture_replays_call_for_call` replays every row in order with a
manual clock that advances its second at the same calls. It requires, for each
call, the same outcome (success, or the refusal matching the two StopWatch
throws), the same nesting depth afterwards and the same output bytes.
`checked_bounds_and_failures_leave_dispatch_state_usable` asserts the
`startProgress(2, 1, "bad")` rows with literal strings copied from that capture.
`gui_and_none_accept_an_inverted_range` checks that a custom GUI backend
receives the inverted range unchanged.

**Independent oracles.** `tests/progress_logger.rs` also includes the source
mode/copy/NONE smoke cases plus independent deterministic clock/backend traces
and literal command-output oracles. It checks same-second and backward-second
updates, two clock reads, copy timestamps, next-counter behavior, shared and
bounded nesting, GUI replacement, duration/rate boundaries, CPU-unavailable
output, and checked I/O, clock, overflow, and zero-duration failures. Apart from
the Release rows named above, output literals are derived from source
expressions and are not claimed as captured C++ output.

`tests/data/progress_logger_provenance.json` pins the original source/header/test
and formatting dependencies. Its `release_oracle` section records the install,
its manifest and library hashes, the pins, the method, the masking and the raw
capture hashes. `external_reference_artifacts` holds the driver and scripts,
which stay outside this repository. A live-clock sanity test checks finite
nondecreasing total process CPU samples without requiring a minimum elapsed
time or workload.
**Consumer calls (tier 1).** `../oracle/progress-consumers/driver.cpp` links
the same Release install and runs the seven consumers above twice per case on a
fresh object: once with a recording backend installed through `setLogger`,
which captures every call that reaches it with its arguments and the depth the
wrapper passes, and once with `setLogType(CMD)`, whose stdout bytes are cut out
by file offset. The driver defines `time()` itself (linked with `-rdynamic`, so
`libOpenMS.so` binds to it) and returns a new second on every call, so the
whole-second throttle never suppresses a set. Three runs on kim are
byte-identical; `tests/data/progress_consumers_release.tsv` is the masked table
and `tests/progress_consumers.rs` replays all 18 cases in both modes. See
`tests/data/progress_consumers_provenance.json`.

**FORMAT reader calls (tier 1).** `../oracle/progress-format-readers/driver.cpp`
links the same Release install and runs 31 cases over DTA2D, MS2, MGF,
mzIdentML, consensusXML, featureXML, mzData, mzXML and mzML on ibminode06, each
twice on a fresh file object: once with `setLogType(GUI)` and
`make_gui_progress_logger` replaced by a factory for a recording backend, so
the fresh backends the files make from their log type and `MzMLHandler`'s
`pg_outer` copy are recorded as well as the file's own, and once with
`setLogType(CMD)`. The inputs are the class-test files at `bc9cc12` and
derived ones (`make_inputs.py`): truncated XML documents, a DTA2D line with a
bad number, and an MGF whose last `END IONS` has no newline. Three runs are
identical after masking the timing and throughput texts;
`tests/data/progress_format_readers_release.tsv` is the masked table and
`tests/progress_format_readers.rs` replays it with the port's GUI factory
installing the same kinds of backend, matching every call, depth and stdout
byte except the differences listed under [FORMAT readers](#format-readers).
The replay also found a reader defect that is not a progress one: the mzXML
reader accepts a document truncated after a complete `</scan>`, where the
Release build throws `ParseError`; the case records it. See
`tests/data/progress_format_readers_provenance.json`.

The main-crate checks are recorded in [validation results](VALIDATION.md).
