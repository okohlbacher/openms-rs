// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Replaceable progress reporting with source-compatible whole-second throttling.

use crate::{Error, Result};
use std::io::{self, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, OnceLock};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

/// Source enum values, including GUI's default no-op backend.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum ProgressLogType {
    Cmd = 0,
    Gui = 1,
    #[default]
    None = 2,
}

/// Absolute clock samples. Only differences of wall/CPU seconds are displayed.
/// `wall_second` independently drives the source's wall-clock bucket throttle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProgressTime {
    pub wall_second: i64,
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

/// Four operations of the source backend, with owned replacement and checked
/// errors. An explicit `set_progress` need not update `next_progress`'s counter.
pub trait ProgressBackend: Send {
    fn start_progress(&mut self, begin: i64, end: i64, label: &str, depth: usize) -> Result<()>;
    fn set_progress(&mut self, value: i64, depth: usize) -> Result<()>;
    fn next_progress(&mut self) -> Result<i64>;
    fn end_progress(&mut self, depth: usize, bytes_processed: u64) -> Result<()>;
}

/// Bounded indentation and labels prevent caller-controlled output allocations.
pub const MAX_PROGRESS_DEPTH: usize = 1024;
pub const MAX_PROGRESS_LABEL_BYTES: usize = 1024 * 1024;

/// Shared source-style nesting. `Default` creates an isolated context; ordinary
/// loggers use `global`, retaining nesting across separate logger instances.
#[derive(Clone, Debug, Default)]
pub struct ProgressNesting(Arc<AtomicUsize>);
impl ProgressNesting {
    pub fn global() -> Self {
        static GLOBAL: OnceLock<ProgressNesting> = OnceLock::new();
        GLOBAL.get_or_init(Self::default).clone()
    }
    pub fn depth(&self) -> usize {
        self.0.load(Ordering::Relaxed)
    }
    fn increment(&self) -> Result<()> {
        self.0
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |value| {
                (value < MAX_PROGRESS_DEPTH).then(|| value + 1)
            })
            .map(|_| ())
            .map_err(|_| invalid("progress nesting limit exceeded"))
    }
    fn decrement(&self) -> usize {
        self.0
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |value| {
                Some(value.saturating_sub(1))
            })
            .expect("the saturating update always succeeds")
            .saturating_sub(1)
    }
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
pub struct ProgressLogger {
    log_type: ProgressLogType,
    last_invoke: i64,
    backend: Box<dyn ProgressBackend>,
    clock: ProgressClock,
    nesting: ProgressNesting,
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
            nesting: self.nesting.clone(),
            gui_factory: self.gui_factory.clone(),
        }
    }
}
impl ProgressLogger {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn with_clock_and_nesting(clock: ProgressClock, nesting: ProgressNesting) -> Self {
        Self {
            log_type: ProgressLogType::None,
            last_invoke: 0,
            backend: Box::new(NoProgress),
            clock,
            nesting,
            gui_factory: Arc::new(|| Box::new(NoProgress)),
        }
    }
    pub fn log_type(&self) -> ProgressLogType {
        self.log_type
    }
    pub fn set_log_type(&mut self, log_type: ProgressLogType) {
        self.backend = self.make_backend(log_type);
        self.log_type = log_type;
    }
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
    pub fn start_progress(&mut self, begin: i64, end: i64, label: &str) -> Result<()> {
        check_start(begin, end, label, self.nesting.depth())?;
        if self.nesting.depth() >= MAX_PROGRESS_DEPTH {
            return Err(invalid("progress nesting limit exceeded"));
        }
        self.last_invoke = (self.clock)()?.wall_second;
        self.backend
            .start_progress(begin, end, label, self.nesting.depth())?;
        self.nesting.increment()
    }
    pub fn set_progress(&mut self, value: i64) -> Result<()> {
        // The source calls time() twice when dispatching, once when suppressed.
        // A backward clock adjustment to another bucket also dispatches.
        if self.last_invoke == (self.clock)()?.wall_second {
            return Ok(());
        }
        self.last_invoke = (self.clock)()?.wall_second;
        self.backend.set_progress(value, self.nesting.depth())
    }
    pub fn next_progress(&mut self) -> Result<()> {
        let value = self.backend.next_progress()?;
        self.set_progress(value)
    }
    pub fn end_progress(&mut self, bytes_processed: u64) -> Result<()> {
        let depth = self.nesting.decrement();
        self.backend.end_progress(depth, bytes_processed)
    }
}

/// Command backend writing through `std::io::Write`. It preserves source f32
/// percentage arithmetic, labels, indentation, diagnostics, and timer formatting.
/// Missing native process CPU timing is displayed as `unavailable (CPU)`.
pub struct CommandProgressLogger<W: Write> {
    writer: W,
    clock: ProgressClock,
    begin: i64,
    end: i64,
    current: i64,
    started: Option<ProgressTime>,
}
impl<W: Write> CommandProgressLogger<W> {
    pub fn new(writer: W) -> Self {
        Self::with_clock(writer, system_progress_clock())
    }
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
    pub fn into_inner(self) -> W {
        self.writer
    }
}
impl<W: Write + Send> ProgressBackend for CommandProgressLogger<W> {
    fn start_progress(&mut self, begin: i64, end: i64, label: &str, depth: usize) -> Result<()> {
        check_start(begin, end, label, depth)?;
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
        // active StopWatch, so the subsequent start() throws after that reset.
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
fn check_depth(depth: usize) -> Result<()> {
    if depth > MAX_PROGRESS_DEPTH {
        Err(invalid("progress nesting limit exceeded"))
    } else {
        Ok(())
    }
}
fn check_start(begin: i64, end: i64, label: &str, depth: usize) -> Result<()> {
    if begin > end {
        return Err(invalid("invalid progress range"));
    }
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
