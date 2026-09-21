// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Replaceable progress reporting with source-compatible whole-second throttling.
//!
//! Port of `CONCEPT/ProgressLogger.h` and `ProgressLogger.cpp`. Behaviour
//! follows the reference Linux x86_64 **Release** build, which compiles the
//! source's `OPENMS_PRECONDITION` checks out; see
//! `docs/PROGRESS_LOGGER_SUPPORT.md` for the API mapping, the native
//! differences and the executed Release evidence.

use crate::{Error, Result};
use std::io::{self, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, MutexGuard, OnceLock, PoisonError};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

/// Possible log types (source `ProgressLogger::LogType`), with the source
/// discriminants. The default is `None`, as in the source.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum ProgressLogType {
    /// Command-line progress (source `CMD`), written to stdout.
    Cmd = 0,
    /// Progress dialog (source `GUI`). The backend comes from the logger's GUI
    /// factory, which defaults to the no-op backend as in the core library.
    Gui = 1,
    /// No progress logging (source `NONE`); every operation is a no-op.
    #[default]
    None = 2,
}

/// Absolute clock samples. Only differences of wall/CPU seconds are displayed.
/// `wall_second` independently drives the source's wall-clock bucket throttle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProgressTime {
    /// Whole civil seconds, the counterpart of the source's `time(nullptr)`.
    pub wall_second: i64,
    /// Absolute elapsed wall time in seconds for the command timer.
    pub wall_seconds: f64,
    /// `None` means that process CPU timing is unavailable, not zero CPU usage.
    pub cpu_seconds: Option<f64>,
}

/// Injectable timing for deterministic callbacks and command-line output.
pub type ProgressClock = Arc<dyn Fn() -> Result<ProgressTime> + Send + Sync>;
/// Owned GUI backend factory. Changing it affects subsequent GUI selection/copies.
pub type ProgressBackendFactory = Arc<dyn Fn() -> Box<dyn ProgressBackend> + Send + Sync>;

/// Civil seconds for throttling and monotonic elapsed time for the command
/// timer. Unix/Windows process CPU time uses the checked `cpu-time` API;
/// other platforms explicitly report CPU timing as unavailable.
pub fn system_progress_clock() -> ProgressClock {
    let origin = Instant::now();
    Arc::new(move || {
        let wall_second = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_err(|_| invalid("system clock precedes the Unix epoch"))?
            .as_secs();
        #[cfg(any(unix, windows))]
        let cpu_seconds = Some(
            cpu_time::ProcessTime::try_now()?
                .as_duration()
                .as_secs_f64(),
        );
        #[cfg(not(any(unix, windows)))]
        let cpu_seconds = None;
        Ok(ProgressTime {
            wall_second: i64::try_from(wall_second)
                .map_err(|_| invalid("system clock exceeds signed seconds"))?,
            wall_seconds: origin.elapsed().as_secs_f64(),
            cpu_seconds,
        })
    })
}

/// Four operations of the source backend (`ProgressLogger::ProgressLoggerImpl`),
/// with owned replacement and checked errors. An explicit `set_progress` need
/// not update `next_progress`'s counter.
pub trait ProgressBackend: Send {
    /// Starts a section. `begin` and `end` arrive exactly as the caller passed
    /// them, including `begin > end`: the source's only range check is a
    /// Debug-only precondition, absent from the Release build.
    fn start_progress(&mut self, begin: i64, end: i64, label: &str, depth: usize) -> Result<()>;
    /// Shows `value`; `depth` is the nesting depth after the section's start.
    fn set_progress(&mut self, value: i64, depth: usize) -> Result<()>;
    /// Advances the backend's own counter and returns it; shows nothing.
    fn next_progress(&mut self) -> Result<i64>;
    /// Finalizes a section; `bytes_processed` (0 = none) requests a rate.
    fn end_progress(&mut self, depth: usize, bytes_processed: u64) -> Result<()>;
}

/// Native bound on the sections open at once, and so on indentation (two
/// spaces per level). The source's `static int recursion_depth_` has no
/// limit; this bound keeps caller-controlled indentation allocations finite.
/// It is not a source check. Sections abandoned by a finished call (see
/// [`ProgressReporter`]) do not count against it, and no backend is ever
/// handed a deeper depth.
pub const MAX_PROGRESS_DEPTH: usize = 1024;
/// Native bound on a label's length in bytes (1 MiB), limiting the
/// caller-controlled header allocation. The source accepts any label; this
/// bound is not a source check.
pub const MAX_PROGRESS_LABEL_BYTES: usize = 1024 * 1024;

/// Shared source-style nesting. `Default` creates an isolated context; ordinary
/// loggers use `global`, retaining nesting across separate logger instances.
///
/// The source's depth is one `static int` (`ProgressLogger.h:105`) that only
/// `endProgress` decrements (`ProgressLogger.cpp:266-269`), so a section an
/// exception leaves open stays in it for the rest of the process and indents
/// every later section of every object. Here such a section, once the call
/// that started it has finished (see [`ProgressReporter`]), is *abandoned*:
/// it stays in [`depth`](Self::depth) and indents the later calls of the
/// logger that started it, and of that logger's copies, as the source's
/// static does, until the last of them is dropped. It does not indent the
/// calls of any other logger and does not count against
/// [`MAX_PROGRESS_DEPTH`], so a failed call cannot change what a later call
/// through another logger does, and no number of them makes a later call fail.
#[derive(Clone, Debug, Default)]
pub struct ProgressNesting(Arc<Mutex<Levels>>);

/// The sections of one nesting context.
#[derive(Clone, Copy, Debug, Default)]
struct Levels {
    /// Started and neither ended nor abandoned. Every logger on the context
    /// sees them, as every source object sees the static depth.
    open: usize,
    /// Abandoned, summed over the loggers that abandoned them; each logger
    /// sees only its own (see [`Owner::abandoned`]).
    abandoned: usize,
}

impl ProgressNesting {
    /// The process-wide context, the counterpart of the source's static
    /// `ProgressLogger::recursion_depth_`.
    pub fn global() -> Self {
        static GLOBAL: OnceLock<ProgressNesting> = OnceLock::new();
        GLOBAL.get_or_init(Self::default).clone()
    }
    /// The current nesting depth: successful starts minus ends, never below 0,
    /// counting the sections abandoned by loggers that are still alive.
    pub fn depth(&self) -> usize {
        let levels = self.levels();
        levels.open.saturating_add(levels.abandoned)
    }
    /// Every count changes under this lock. No backend runs while it is held.
    fn levels(&self) -> MutexGuard<'_, Levels> {
        // Nothing panics while the lock is held, and the counts are valid at
        // every step, so a poisoned lock still holds consistent counts.
        self.0.lock().unwrap_or_else(PoisonError::into_inner)
    }
}

/// What one logger and its copies share: the counterpart of one source file
/// object, whose sections its handlers' copies of its logger also report
/// (`MzMLHandler.cpp:135`, `FeatureXMLFile.cpp:54`). Its counts change only
/// under the lock of its nesting context, which orders them, so relaxed
/// atomics suffice.
struct Owner {
    nesting: ProgressNesting,
    /// Of the context's open sections, those this owner's loggers started.
    /// Counts, not identities: an end closes the innermost open section,
    /// whoever started it, so this can exceed what is left, and is capped by
    /// the context's count whenever it is used.
    open: AtomicUsize,
    /// The sections this owner abandoned, included in the context's
    /// `abandoned` until the owner is dropped.
    abandoned: AtomicUsize,
    /// The live [`ProgressReporter`]s on this owner's loggers.
    calls: AtomicUsize,
    /// `open` when the first of them began: sections the caller started on
    /// the logger before its call, which the call does not abandon.
    before_call: AtomicUsize,
}

impl Owner {
    fn new(nesting: ProgressNesting) -> Arc<Self> {
        Arc::new(Self {
            nesting,
            open: AtomicUsize::new(0),
            abandoned: AtomicUsize::new(0),
            calls: AtomicUsize::new(0),
            before_call: AtomicUsize::new(0),
        })
    }
    /// The depth this owner's calls are dispatched at: every open section plus
    /// its own abandoned ones, at most [`MAX_PROGRESS_DEPTH`].
    fn view(&self, levels: &Levels) -> usize {
        levels
            .open
            .saturating_add(self.abandoned.load(Ordering::Relaxed))
            .min(MAX_PROGRESS_DEPTH)
    }
    /// The dispatch depth, and whether another section may open.
    fn depth_and_room(&self) -> (usize, bool) {
        let levels = self.nesting.levels();
        (self.view(&levels), levels.open < MAX_PROGRESS_DEPTH)
    }
    fn depth(&self) -> usize {
        self.view(&self.nesting.levels())
    }
    /// Counts a started section.
    fn open_section(&self) -> Result<()> {
        let mut levels = self.nesting.levels();
        if levels.open >= MAX_PROGRESS_DEPTH {
            return Err(invalid("progress nesting limit exceeded"));
        }
        levels.open += 1;
        self.open.fetch_add(1, Ordering::Relaxed);
        Ok(())
    }
    /// Counts an end, returning the depth after it. Like the source's, an end
    /// with nothing open changes nothing; the innermost section is an open
    /// one if any is, else one this owner abandoned.
    fn close_section(&self) -> usize {
        let mut levels = self.nesting.levels();
        if levels.open > 0 {
            levels.open -= 1;
            decrement(&self.open);
        } else if self.abandoned.load(Ordering::Relaxed) > 0 {
            levels.abandoned = levels.abandoned.saturating_sub(1);
            decrement(&self.abandoned);
        }
        self.view(&levels)
    }
    /// A reporter on one of this owner's loggers begins.
    fn enter(&self) {
        let _levels = self.nesting.levels();
        if self.calls.fetch_add(1, Ordering::Relaxed) == 0 {
            self.before_call
                .store(self.open.load(Ordering::Relaxed), Ordering::Relaxed);
        }
    }
    /// A reporter ends. When it was the last, the call is over, and the
    /// sections it started and left open are abandoned.
    fn leave(&self) {
        let mut levels = self.nesting.levels();
        if self.calls.fetch_sub(1, Ordering::Relaxed) != 1 {
            return;
        }
        let open = self.open.load(Ordering::Relaxed);
        let left = open
            .saturating_sub(self.before_call.load(Ordering::Relaxed))
            .min(levels.open);
        self.open.store(open - left, Ordering::Relaxed);
        levels.open -= left;
        levels.abandoned = levels.abandoned.saturating_add(left);
        self.abandoned.fetch_add(left, Ordering::Relaxed);
    }
}

impl Drop for Owner {
    /// The last of a logger and its copies is gone: no call can end its
    /// sections any more, so they leave the context. Their backends are not
    /// called, so nothing is printed and no `-- done` line appears.
    fn drop(&mut self) {
        let mut levels = self.nesting.levels();
        let open = self.open.load(Ordering::Relaxed).min(levels.open);
        levels.open -= open;
        let abandoned = self.abandoned.load(Ordering::Relaxed).min(levels.abandoned);
        levels.abandoned -= abandoned;
    }
}

/// Decrements a count that is not already zero.
fn decrement(count: &AtomicUsize) {
    // The closure never declines, so the result carries no information.
    let _ = count.fetch_update(Ordering::Relaxed, Ordering::Relaxed, |value| {
        Some(value.saturating_sub(1))
    });
}

#[derive(Default)]
struct NoProgress;
impl ProgressBackend for NoProgress {
    fn start_progress(&mut self, _: i64, _: i64, _: &str, _: usize) -> Result<()> {
        Ok(())
    }
    fn set_progress(&mut self, _: i64, _: usize) -> Result<()> {
        Ok(())
    }
    fn next_progress(&mut self) -> Result<i64> {
        Ok(0)
    }
    fn end_progress(&mut self, _: usize, _: u64) -> Result<()> {
        Ok(())
    }
}

/// Progress dispatch. Replacing a backend does not change its reported log type;
/// selecting any type always creates a new backend, even when the type is equal.
///
/// Clone and `clone_from` copy the type and throttle timestamp, but create a
/// fresh backend using that type. They never clone an active/custom backend.
/// A copy shares the original's sections, as a source handler's copy of its
/// file object's logger reports that object's sections: sections abandoned by
/// either indent the calls of both (see [`ProgressNesting`]).
///
/// Dropping a logger never calls its backend: no section is ended and no
/// `-- done` line is printed, as the source's destructor ends nothing
/// (`ProgressLogger.cpp:192-195`). Unlike the source's, dropping the last of
/// a logger and its copies takes the sections it left open out of the
/// nesting depth, since no call can end them any more.
pub struct ProgressLogger {
    log_type: ProgressLogType,
    last_invoke: i64,
    backend: Box<dyn ProgressBackend>,
    clock: ProgressClock,
    owner: Arc<Owner>,
    gui_factory: ProgressBackendFactory,
}
impl Default for ProgressLogger {
    fn default() -> Self {
        Self::with_clock_and_nesting(system_progress_clock(), ProgressNesting::global())
    }
}
impl Clone for ProgressLogger {
    fn clone(&self) -> Self {
        Self {
            log_type: self.log_type,
            last_invoke: self.last_invoke,
            backend: self.make_backend(self.log_type),
            clock: self.clock.clone(),
            owner: self.owner.clone(),
            gui_factory: self.gui_factory.clone(),
        }
    }
}
impl ProgressLogger {
    /// A logger of type `None` using the system clock and the global nesting
    /// context (source default constructor).
    pub fn new() -> Self {
        Self::default()
    }
    /// A logger of type `None` with an injected `clock` and `nesting` context,
    /// for deterministic throttling and timing or isolated nesting.
    pub fn with_clock_and_nesting(clock: ProgressClock, nesting: ProgressNesting) -> Self {
        Self {
            log_type: ProgressLogType::None,
            last_invoke: 0,
            backend: Box::new(NoProgress),
            clock,
            owner: Owner::new(nesting),
            gui_factory: Arc::new(|| Box::new(NoProgress)),
        }
    }
    /// The type of progress log being used (source `getLogType`).
    pub fn log_type(&self) -> ProgressLogType {
        self.log_type
    }
    /// Selects the progress log type (source `setLogType`); the default is
    /// `None`. Always creates a fresh backend, even for the current type.
    pub fn set_log_type(&mut self, log_type: ProgressLogType) {
        self.backend = self.make_backend(log_type);
        self.log_type = log_type;
    }
    /// Replaces the backend used for progress logging (source `setLogger`),
    /// taking ownership. The reported log type is unchanged.
    pub fn set_logger(&mut self, backend: Box<dyn ProgressBackend>) {
        self.backend = backend;
    }
    /// Per-logger counterpart of the source global GUI factory. Copies retain
    /// this factory; existing backends are not replaced until type selection.
    pub fn set_gui_factory(&mut self, factory: ProgressBackendFactory) {
        self.gui_factory = factory;
    }
    fn make_backend(&self, log_type: ProgressLogType) -> Box<dyn ProgressBackend> {
        match log_type {
            ProgressLogType::None => Box::new(NoProgress),
            ProgressLogType::Gui => (self.gui_factory)(),
            ProgressLogType::Cmd => Box::new(CommandProgressLogger::with_clock(
                io::stdout(),
                self.clock.clone(),
            )),
        }
    }
    /// Initializes the progress display (source `startProgress`).
    ///
    /// Sets the progress range from `begin` to `end` and the label to `label`.
    /// If `begin` equals `end`, [`set_progress`](Self::set_progress) only
    /// indicates that the program is still running, without an absolute
    /// progress value. As the source notes, select a type with
    /// [`set_log_type`](Self::set_log_type) first; the default `None` shows
    /// nothing.
    ///
    /// Every range is accepted and reaches the backend unchanged, including
    /// `begin > end`. The source's `OPENMS_PRECONDITION(begin <= end, ...)`
    /// (`ProgressLogger.cpp:235`) exists only in Debug builds; the reference
    /// Release build has no range check. The command backend stores an
    /// inverted range as given and then reports every value as invalid,
    /// because no value lies inside it.
    ///
    /// Records the current wall-clock second for the throttle, dispatches at
    /// the current nesting depth, then increments the depth. The depth
    /// dispatched is this logger's (see [`ProgressNesting`]): every open
    /// section, and the sections this logger or a copy of it abandoned.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `label` is longer than
    /// [`MAX_PROGRESS_LABEL_BYTES`] or [`MAX_PROGRESS_DEPTH`] sections are
    /// open. These are native bounds, not source checks, and are tested
    /// before any state changes. Clock and backend errors propagate. After a
    /// backend error the depth is not incremented, as in the source when its
    /// backend throws (the command backend's `StopWatch is already started!`
    /// on a second start).
    pub fn start_progress(&mut self, begin: i64, end: i64, label: &str) -> Result<()> {
        let (depth, room) = self.owner.depth_and_room();
        check_start(label, depth)?;
        if !room {
            return Err(invalid("progress nesting limit exceeded"));
        }
        self.last_invoke = (self.clock)()?.wall_second;
        self.backend.start_progress(begin, end, label, depth)?;
        self.owner.open_section()
    }
    /// Sets the current progress (source `setProgress`).
    ///
    /// Does nothing when the current wall-clock second equals the recorded
    /// one; otherwise records the second and dispatches `value` at the current
    /// depth. Any `value` is accepted; the command backend reports one outside
    /// the range as a printed diagnostic, not as an error.
    ///
    /// # Errors
    ///
    /// Clock and backend errors propagate.
    pub fn set_progress(&mut self, value: i64) -> Result<()> {
        // The source calls time() twice when dispatching, once when suppressed.
        // A backward clock adjustment to another bucket also dispatches.
        if self.last_invoke == (self.clock)()?.wall_second {
            return Ok(());
        }
        self.last_invoke = (self.clock)()?.wall_second;
        self.backend.set_progress(value, self.owner.depth())
    }
    /// Increments progress by one within the range (source `nextProgress`).
    ///
    /// The backend counter advances before the same-second throttle, so a
    /// suppressed call still counts; a dispatched call shows the new count.
    ///
    /// # Errors
    ///
    /// Backend (including counter overflow), clock and display errors
    /// propagate.
    pub fn next_progress(&mut self) -> Result<()> {
        let value = self.backend.next_progress()?;
        self.set_progress(value)
    }
    /// Ends the progress display (source `endProgress`).
    ///
    /// `bytes_processed` optionally requests a bytes-per-second estimate; 0
    /// requests none. The depth is decremented first when nonzero, then the
    /// backend is called, even without a matching start.
    ///
    /// # Errors
    ///
    /// Backend errors propagate. The command backend refuses an end without a
    /// running timer, as the Release build's `StopWatch::stop` throws
    /// `StopWatch cannot be stopped if not running!`.
    pub fn end_progress(&mut self, bytes_processed: u64) -> Result<()> {
        let depth = self.owner.close_section();
        self.backend.end_progress(depth, bytes_processed)
    }
}

/// The progress calls of one algorithm run, sent to an optional logger.
///
/// Source algorithms such as `PeakPickerHiRes` or `GaussFilter` *derive* from
/// `ProgressLogger`, so every one of their objects carries a logger, and their
/// `const` methods mutate it through the source's `mutable` members. A Rust
/// algorithm keeps its methods `&self` and takes the logger from its caller
/// instead: its `*_with_progress` entry point borrows a [`ProgressLogger`] for
/// the call, and the caller selects its type with
/// [`ProgressLogger::set_log_type`] or replaces its backend with
/// [`ProgressLogger::set_logger`], as a source caller does on the algorithm
/// object. Its other entry points report nothing, which is the source's default
/// type `NONE` (`ProgressLogger.cpp:127-128`) without its cost.
///
/// `None` makes every call a no-op: no clock is read and the nesting depth is
/// not touched. A source object of type `NONE` still increments and decrements
/// the static depth around each section, which no output can observe, because
/// the no-op backend prints nothing and the section is balanced.
///
/// A reporter on a logger also marks one call of an entry point, which it
/// lasts for. A call can finish with sections it never ended: a reader that
/// fails inside its section returns without ending it, as the source's
/// exception bypasses `endProgress`, and a metadata-only mzML load stops inside
/// its document section. The call lasts while any reporter on the logger or
/// its copies lives. When the last is dropped, the sections the call started
/// and did not end are *abandoned* (see [`ProgressNesting`]): they stay in the
/// depth for that logger alone and no longer count against
/// [`MAX_PROGRESS_DEPTH`]. Sections the caller had started on the logger before
/// the call are not the call's and stay open. No backend is called, so an
/// abandoned section prints no `-- done` line, as in the source.
pub struct ProgressReporter<'a> {
    logger: Option<&'a mut ProgressLogger>,
    /// The call this reporter marks; `None` when silent.
    _call: Option<Call>,
}

/// One live reporter on a logger. It holds the owner rather than the logger,
/// so it needs no borrow: the reporter's borrow of the logger can end at its
/// last use while the call lasts until the reporter is dropped.
struct Call(Arc<Owner>);
impl Call {
    fn new(owner: &Arc<Owner>) -> Self {
        owner.enter();
        Self(owner.clone())
    }
}
impl Drop for Call {
    fn drop(&mut self) {
        self.0.leave();
    }
}

impl ProgressReporter<'static> {
    /// Calls that go nowhere, for the entry points that report no progress.
    pub fn silent() -> Self {
        Self {
            logger: None,
            _call: None,
        }
    }
}

impl<'a> ProgressReporter<'a> {
    /// Calls that go to `logger`, if any.
    pub fn new(logger: Option<&'a mut ProgressLogger>) -> Self {
        let call = logger.as_deref().map(|logger| Call::new(&logger.owner));
        Self {
            logger,
            _call: call,
        }
    }

    /// Whether the calls reach a logger.
    pub fn is_reporting(&self) -> bool {
        self.logger.is_some()
    }

    /// Source `startProgress(begin, end, label)`; see
    /// [`ProgressLogger::start_progress`].
    ///
    /// # Errors
    ///
    /// The errors of [`ProgressLogger::start_progress`]; never when silent.
    pub fn start(&mut self, begin: i64, end: i64, label: &str) -> Result<()> {
        match self.logger.as_deref_mut() {
            Some(logger) => logger.start_progress(begin, end, label),
            None => Ok(()),
        }
    }

    /// Source `startProgress(0, count, label)` for a record count, which the
    /// source passes as its signed `SignedSize`. The count is converted only
    /// when the calls reach a logger, so a silent reporter can never fail here.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `count` exceeds `i64::MAX`, and the errors
    /// of [`ProgressLogger::start_progress`]; never when silent.
    pub fn start_count(&mut self, count: usize, label: &str) -> Result<()> {
        if self.logger.is_none() {
            return Ok(());
        }
        self.start(0, progress_value(count)?, label)
    }

    /// Source `setProgress(value)`; see [`ProgressLogger::set_progress`].
    ///
    /// # Errors
    ///
    /// The errors of [`ProgressLogger::set_progress`]; never when silent.
    pub fn set(&mut self, value: i64) -> Result<()> {
        match self.logger.as_deref_mut() {
            Some(logger) => logger.set_progress(value),
            None => Ok(()),
        }
    }

    /// Source `setProgress(value)` for a record count or index, which the
    /// source passes as its signed `SignedSize`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `value` exceeds `i64::MAX`, which no
    /// in-memory record count reaches, and the errors of
    /// [`ProgressReporter::set`]; never when silent.
    pub fn set_count(&mut self, value: usize) -> Result<()> {
        if self.logger.is_none() {
            return Ok(());
        }
        self.set(progress_value(value)?)
    }

    /// Source `nextProgress()`; see [`ProgressLogger::next_progress`].
    ///
    /// # Errors
    ///
    /// The errors of [`ProgressLogger::next_progress`]; never when silent.
    pub fn next_progress(&mut self) -> Result<()> {
        match self.logger.as_deref_mut() {
            Some(logger) => logger.next_progress(),
            None => Ok(()),
        }
    }

    /// Source `endProgress()`, without a byte count; see
    /// [`ProgressLogger::end_progress`].
    ///
    /// # Errors
    ///
    /// The errors of [`ProgressLogger::end_progress`]; never when silent.
    pub fn end(&mut self) -> Result<()> {
        self.end_with_bytes(0)
    }

    /// Source `endProgress(bytes_processed)`, whose nonzero byte count asks
    /// the command backend for a throughput; see
    /// [`ProgressLogger::end_progress`].
    ///
    /// # Errors
    ///
    /// The errors of [`ProgressLogger::end_progress`]; never when silent.
    pub fn end_with_bytes(&mut self, bytes_processed: u64) -> Result<()> {
        match self.logger.as_deref_mut() {
            Some(logger) => logger.end_progress(bytes_processed),
            None => Ok(()),
        }
    }

    /// One source `startProgress(begin, end, label)` ... `endProgress()`
    /// section around `body`, which makes the section's `setProgress` calls
    /// through the reporter it is handed.
    ///
    /// The section is ended whether `body` succeeds or fails, so a failed run
    /// leaves the logger's nesting depth where it found it and its command
    /// timer stopped, and the logger can report the next run. **This differs
    /// from the source on its failure paths:** when an algorithm throws inside
    /// its section, the source never reaches `endProgress`, so the command
    /// backend prints no `-- done` line, the static depth stays one level
    /// deeper for the rest of the process, and the object's next
    /// `startProgress` throws `StopWatch is already started!`
    /// (`StopWatch.cpp:43`). Here the `-- done` line is printed and the error
    /// is returned.
    ///
    /// # Errors
    ///
    /// The error of `start`, in which case `body` does not run; otherwise the
    /// error of `body`, or, when `body` succeeded, the error of ending the
    /// section. When both `body` and the end fail, the body's error wins.
    pub fn section<T>(
        &mut self,
        begin: i64,
        end: i64,
        label: &str,
        body: impl FnOnce(&mut Self) -> Result<T>,
    ) -> Result<T> {
        self.start(begin, end, label)?;
        let outcome = body(self);
        let ended = self.end();
        let value = outcome?;
        ended?;
        Ok(value)
    }
}

/// A record count or index as the source's signed progress value
/// (`SignedSize`), which the source converts implicitly from `Size`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `count` exceeds `i64::MAX`, where the source's
/// implicit conversion would wrap.
pub fn progress_value(count: usize) -> Result<i64> {
    i64::try_from(count).map_err(|_| invalid("progress value exceeds signed range"))
}

/// Command backend writing through `std::io::Write` (source
/// `CMDProgressLoggerImpl`). It preserves source f32 percentage arithmetic,
/// labels, indentation, diagnostics, and timer formatting. Missing native
/// process CPU timing is displayed as `unavailable (CPU)`.
///
/// Like the Release build, it stores any range as given
/// (`ProgressLogger.cpp:36-40`). Setting a value then prints a dot when begin
/// equals end (`:48-51`), the invalid-value diagnostic when the value lies
/// below begin or above end (`:52-56`), and otherwise the percentage (`:57-63`).
/// For an inverted range every value takes the diagnostic branch, and the
/// percentage branch is only reached with `begin < end`, so it never divides
/// by zero.
pub struct CommandProgressLogger<W: Write> {
    writer: W,
    clock: ProgressClock,
    begin: i64,
    end: i64,
    current: i64,
    started: Option<ProgressTime>,
}
impl<W: Write> CommandProgressLogger<W> {
    /// A backend writing to `writer`, timed by the system clock.
    pub fn new(writer: W) -> Self {
        Self::with_clock(writer, system_progress_clock())
    }
    /// A backend writing to `writer`, timed by the injected `clock`.
    pub fn with_clock(writer: W, clock: ProgressClock) -> Self {
        Self {
            writer,
            clock,
            begin: 0,
            end: 0,
            current: 0,
            started: None,
        }
    }
    /// Returns the writer, with everything written so far.
    pub fn into_inner(self) -> W {
        self.writer
    }
}
impl<W: Write + Send> ProgressBackend for CommandProgressLogger<W> {
    fn start_progress(&mut self, begin: i64, end: i64, label: &str, depth: usize) -> Result<()> {
        // Native label/depth bounds only; any begin/end pair is stored, as in
        // the Release build.
        check_start(label, depth)?;
        let already_running = self.started.is_some();
        self.begin = begin;
        self.end = end;
        self.current = begin;
        writeln!(
            self.writer,
            "{}Progress of '{}':",
            " ".repeat(depth * 2),
            label
        )?;
        self.writer.flush()?;
        // Source excludes header I/O from the timer. Its reset() restarts an
        // active StopWatch, so the subsequent start() throws after that reset
        // (`StopWatch.cpp:43`, an unconditional throw that the Release build
        // executes; not a Debug-only precondition).
        let started = (self.clock)()?;
        check_time(started)?;
        self.started = Some(started);
        if already_running {
            return Err(invalid("progress timer is already running"));
        }
        Ok(())
    }
    fn set_progress(&mut self, value: i64, depth: usize) -> Result<()> {
        check_depth(depth)?;
        if self.begin == self.end {
            write!(self.writer, ".")?;
        } else if value < self.begin || value > self.end {
            writeln!(
                self.writer,
                "ProgressLogger: Invalid progress value '{}'. Should be between '{}' and '{}'!",
                value, self.begin, self.end
            )?;
        } else {
            // Reached only with begin < end. The source subtracts without a
            // check; an i64 difference that overflows is undefined there and an
            // error here.
            let distance = value
                .checked_sub(self.begin)
                .ok_or_else(|| invalid("progress difference exceeds signed range"))?;
            let range = self
                .end
                .checked_sub(self.begin)
                .ok_or_else(|| invalid("progress range exceeds signed arithmetic"))?;
            let percent = distance as f32 / range as f32 * 100.0_f32;
            write!(
                self.writer,
                "\r{}{percent:.2} %               ",
                " ".repeat(depth * 2)
            )?;
        }
        self.writer.flush()?;
        Ok(())
    }
    fn next_progress(&mut self) -> Result<i64> {
        self.current = self
            .current
            .checked_add(1)
            .ok_or_else(|| invalid("progress counter overflow"))?;
        Ok(self.current)
    }
    fn end_progress(&mut self, depth: usize, bytes_processed: u64) -> Result<()> {
        check_depth(depth)?;
        // Release `StopWatch::stop` throws unconditionally here
        // (`StopWatch.cpp:55`) before anything is printed.
        let start = self
            .started
            .ok_or_else(|| invalid("progress timer is not running"))?;
        let end = (self.clock)()?;
        check_time(end)?;
        let wall = end.wall_seconds - start.wall_seconds;
        let wall_text = duration_text(wall)?;
        let cpu = match (start.cpu_seconds, end.cpu_seconds) {
            (Some(start), Some(end)) => duration_text(end - start)?,
            _ => "unavailable".into(),
        };
        let throughput = if bytes_processed == 0 {
            String::new()
        } else {
            let rate = bytes_processed as f64 / wall;
            // Floating-to-unsigned conversion is defined in the source only in
            // this interval. Zero wall time must not invent an infinite rate.
            if !rate.is_finite() || !(0.0..18_446_744_073_709_551_616.0).contains(&rate) {
                return Err(invalid("progress throughput is outside finite u64 range"));
            }
            format!(" @ {}/s", bytes_text(rate as u64))
        };
        writeln!(
            self.writer,
            "\r{}-- done [took {} (CPU), {} (Wall){}] -- ",
            " ".repeat(depth * 2),
            cpu,
            wall_text,
            throughput
        )?;
        self.writer.flush()?;
        self.started = None;
        Ok(())
    }
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
/// Native indentation bound for backend calls; not a source check.
fn check_depth(depth: usize) -> Result<()> {
    if depth > MAX_PROGRESS_DEPTH {
        Err(invalid("progress nesting limit exceeded"))
    } else {
        Ok(())
    }
}
/// Native start bounds. There is deliberately no range check: the source's
/// `begin <= end` precondition is Debug-only and absent from the Release build.
fn check_start(label: &str, depth: usize) -> Result<()> {
    if label.len() > MAX_PROGRESS_LABEL_BYTES {
        return Err(invalid("progress label limit exceeded"));
    }
    check_depth(depth)
}
fn check_time(time: ProgressTime) -> Result<()> {
    if !time.wall_seconds.is_finite()
        || time.wall_seconds < 0.0
        || time
            .cpu_seconds
            .is_some_and(|value| !value.is_finite() || value < 0.0)
    {
        return Err(invalid("progress timing must be finite and nonnegative"));
    }
    Ok(())
}
fn duration_text(seconds: f64) -> Result<String> {
    // Source day decomposition multiplies signed ints. Bound the arithmetic
    // rather than reproducing overflow for multi-decade timer intervals.
    if !seconds.is_finite() || !(0.0..=f64::from(i32::MAX)).contains(&seconds) {
        return Err(invalid("progress elapsed duration exceeds checked range"));
    }
    let whole = seconds as u64;
    let (days, hours, minutes, secs) = (
        whole / 86400,
        whole / 3600 % 24,
        whole / 60 % 60,
        whole % 60,
    );
    Ok(if days > 0 {
        format!("{days}d {hours:02}:{minutes:02}:{secs:02} h")
    } else if hours > 0 {
        format!("{hours:02}:{minutes:02}:{secs:02} h")
    } else if minutes > 0 {
        format!("{minutes:02}:{secs:02} m")
    } else {
        format!("{seconds:.2} s")
    })
}
fn bytes_text(bytes: u64) -> String {
    let mut value = bytes as f64;
    for unit in ["byte", "KiB", "MiB", "GiB", "TiB", "PiB"] {
        if value < 1024.0 {
            // Every nonzero value is >=1 after binary-unit selection, so four
            // significant digits never requires exponent notation here.
            let decimals = if value == 0.0 {
                0
            } else {
                (3 - value.log10().floor() as i32).max(0) as usize
            };
            let mut number = format!("{value:.decimals$}");
            if number.contains('.') {
                while number.ends_with('0') {
                    number.pop();
                }
                if number.ends_with('.') {
                    number.pop();
                }
            }
            return format!("{number} {unit}");
        }
        value /= 1024.0;
    }
    format!("Congrats. That's a lot of bytes: {bytes}")
}
