// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Wall-clock and process CPU timing of the Core SDK `SYSTEM/StopWatch.h`.
//!
//! [`StopWatch`](crate::system::stop_watch::StopWatch) accumulates the time
//! spent between `start`/`resume` and `stop`, so intermediate steps can be
//! excluded from a measurement. CPU time is the sum over all threads of the
//! current process. A watch must be started before a time can be read, but it
//! need not be stopped: reading a running watch adds the interval since the
//! last start to the accumulated total.
//!
//! Two things differ from the source and both are visible in the API.
//! [`StopWatch::clock_time`](crate::system::stop_watch::StopWatch::clock_time)
//! is monotonic here, because the port samples [`std::time::Instant`] where the
//! source samples `gettimeofday`, which an administrator or NTP can step
//! backwards. And CPU time is an [`Option`]: the source reaches the per-process
//! user/kernel split through the POSIX `times` system call on every Unix, which
//! this port cannot do without a libc dependency, so the split is read from
//! `/proc/self/stat` on Linux and is reported as unavailable — never as zero —
//! elsewhere. See `docs/STOP_WATCH_SUPPORT.md` for the platform matrix.

use crate::{Error, Result};
use std::sync::OnceLock;
use std::time::Instant;

/// Microseconds in one second; the source's `1000000L` normalisation constant.
const MICROS_PER_SECOND: i64 = 1_000_000;

/// Microseconds in one `USER_HZ` tick, the unit of `/proc/[pid]/stat`.
///
/// The Linux kernel reports `utime`/`stime` in `USER_HZ`, which is fixed at 100
/// for that interface regardless of the configured timer frequency, so this is
/// a constant rather than a `sysconf(_SC_CLK_TCK)` query. The source divides by
/// `sysconf(_SC_CLK_TCK)`, which glibc also answers with 100 on Linux.
#[cfg(target_os = "linux")]
const MICROS_PER_TICK: i64 = MICROS_PER_SECOND / 100;

/// Largest `/proc` record this module will read, in bytes.
#[cfg(target_os = "linux")]
const MAX_PROC_BYTES: u64 = 64 * 1024;

/// Largest magnitude accepted by
/// [`StopWatch::format_seconds`](crate::system::stop_watch::StopWatch::format_seconds).
///
/// The source casts its `double` argument to an integer type without a range
/// check, which is undefined behaviour for a value that does not fit; this port
/// refuses instead. The bound is well inside [`i64`] and far beyond any
/// measurable duration — 9e18 seconds is about 285 billion years.
pub const MAX_FORMATTABLE_SECONDS: f64 = 9.0e18;

/// Process CPU accounting as one platform reports it.
///
/// The source's `TimeDiff_` always carries a user/kernel pair because `times`
/// and `GetProcessTimes` both provide one. Here the pair exists only where the
/// port can read it without a libc dependency, so the variant records which of
/// the two shapes a sample has. Mixing shapes is not meaningful, and the
/// arithmetic below reports such a mix as unavailable rather than guessing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CpuSample {
    /// User and kernel time reported separately, in microseconds.
    ///
    /// Read from `/proc/self/stat` on Linux, which aggregates over every thread
    /// of the process, as the source's `times` does.
    Split {
        /// Time executing in user mode, in microseconds.
        user_micros: i64,
        /// Time executing in kernel mode on the process's behalf, in microseconds.
        kernel_micros: i64,
    },
    /// Only the aggregate process CPU total is available, in microseconds.
    ///
    /// Used on non-Linux Unix and on Windows, where the port reads the
    /// process CPU clock but has no portable access to the split.
    Total {
        /// Total process CPU time, user plus kernel, in microseconds.
        micros: i64,
    },
}

impl CpuSample {
    /// Time in user mode, or `None` when the platform reports no split.
    pub fn user_micros(self) -> Option<i64> {
        match self {
            Self::Split { user_micros, .. } => Some(user_micros),
            Self::Total { .. } => None,
        }
    }

    /// Time in kernel mode, or `None` when the platform reports no split.
    pub fn kernel_micros(self) -> Option<i64> {
        match self {
            Self::Split { kernel_micros, .. } => Some(kernel_micros),
            Self::Total { .. } => None,
        }
    }

    /// Total process CPU time, user plus kernel.
    pub fn total_micros(self) -> i64 {
        match self {
            Self::Split {
                user_micros,
                kernel_micros,
            } => user_micros.saturating_add(kernel_micros),
            Self::Total { micros } => micros,
        }
    }

    fn difference(self, earlier: Self) -> Option<Self> {
        match (self, earlier) {
            (
                Self::Split {
                    user_micros: user,
                    kernel_micros: kernel,
                },
                Self::Split {
                    user_micros: earlier_user,
                    kernel_micros: earlier_kernel,
                },
            ) => Some(Self::Split {
                user_micros: user.saturating_sub(earlier_user),
                kernel_micros: kernel.saturating_sub(earlier_kernel),
            }),
            (Self::Total { micros }, Self::Total { micros: earlier }) => Some(Self::Total {
                micros: micros.saturating_sub(earlier),
            }),
            _ => None,
        }
    }

    fn sum(self, other: Self) -> Option<Self> {
        match (self, other) {
            (
                Self::Split {
                    user_micros: user,
                    kernel_micros: kernel,
                },
                Self::Split {
                    user_micros: other_user,
                    kernel_micros: other_kernel,
                },
            ) => Some(Self::Split {
                user_micros: user.saturating_add(other_user),
                kernel_micros: kernel.saturating_add(other_kernel),
            }),
            (Self::Total { micros }, Self::Total { micros: other }) => Some(Self::Total {
                micros: micros.saturating_add(other),
            }),
            _ => None,
        }
    }
}

/// One clock reading, or one accumulated interval: the source's `TimeDiff_`.
///
/// Wall time is kept as whole seconds plus an *additional* microsecond count,
/// exactly as the source keeps the two members `gettimeofday` fills in, so that
/// the carry and borrow rules below — and their observable effect on equality —
/// are the source's own. Absolute wall readings are measured from a
/// process-wide monotonic origin rather than the Unix epoch; only differences
/// are ever reported, so the origin is not observable.
///
/// All arithmetic saturates at the [`i64`] bounds instead of wrapping. Reaching
/// them requires roughly 292 billion years of accumulated time, so it cannot
/// occur, but the port does not depend on that to stay panic-free.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TimeSample {
    /// Whole seconds of wall-clock time.
    pub wall_seconds: i64,
    /// Microseconds of wall-clock time *in addition to* [`Self::wall_seconds`].
    pub wall_micros: i64,
    /// Process CPU accounting, or `None` where the platform reports none.
    pub cpu: Option<CpuSample>,
}

impl Default for TimeSample {
    /// A zero reading, with the CPU shape this platform would report.
    ///
    /// The CPU field is a zeroed [`CpuSample`] rather than `None` wherever the
    /// platform has an implementation, so a freshly cleared watch reports
    /// `Some(0.0)` — the source's `0.0` — instead of claiming the platform
    /// cannot measure CPU time.
    fn default() -> Self {
        Self {
            wall_seconds: 0,
            wall_micros: 0,
            cpu: zero_cpu(),
        }
    }
}

impl TimeSample {
    /// Read the current wall clock and process CPU time.
    ///
    /// This is the source's private `snapShot_`. It never fails: a platform
    /// that cannot report CPU time, or a `/proc` read that does not succeed,
    /// yields `cpu: None`, which propagates through the arithmetic so that the
    /// watch reports CPU time as unavailable rather than as zero.
    pub fn now() -> Self {
        let elapsed = origin().elapsed();
        Self {
            wall_seconds: i64::try_from(elapsed.as_secs()).unwrap_or(i64::MAX),
            wall_micros: i64::from(elapsed.subsec_micros()),
            cpu: sample_cpu(),
        }
    }

    /// Accumulated wall-clock time, in seconds.
    ///
    /// Reproduces the source's `(double)start_time + (double)start_time_usec / 1e6`
    /// term by term, including the order of the two conversions.
    pub fn clock_time(self) -> f64 {
        self.wall_seconds as f64 + self.wall_micros as f64 / 1e6
    }

    /// Accumulated user time in seconds, or `None` when the split is unavailable.
    pub fn user_time(self) -> Option<f64> {
        self.cpu?.user_micros().map(micros_to_seconds)
    }

    /// Accumulated kernel time in seconds, or `None` when the split is unavailable.
    pub fn kernel_time(self) -> Option<f64> {
        self.cpu?.kernel_micros().map(micros_to_seconds)
    }

    /// Accumulated process CPU time in seconds, or `None` when unavailable.
    ///
    /// Where the platform reports a split this is the sum of the two integer
    /// microsecond counts; [`StopWatch::cpu_time`] prefers the source's own
    /// `user + kernel` sum of the two already-converted `f64` values and only
    /// falls back to this on a platform that reports no split.
    pub fn cpu_time(self) -> Option<f64> {
        self.cpu.map(|cpu| micros_to_seconds(cpu.total_micros()))
    }

    /// This reading minus an earlier one: the source's `TimeDiff_::operator-`.
    ///
    /// Borrows a whole second for every negative microsecond remainder, as the
    /// source's `while (diff.start_time_usec < 0L)` loop does, in constant time
    /// rather than one iteration per second.
    pub fn difference(self, earlier: Self) -> Self {
        let seconds = self.wall_seconds.saturating_sub(earlier.wall_seconds);
        let micros = self.wall_micros.saturating_sub(earlier.wall_micros);
        let (wall_seconds, wall_micros) = borrow_negative_micros(seconds, micros);
        Self {
            wall_seconds,
            wall_micros,
            cpu: match (self.cpu, earlier.cpu) {
                (Some(later), Some(earlier)) => later.difference(earlier),
                _ => None,
            },
        }
    }

    /// This reading plus another: the source's `TimeDiff_::operator+=`.
    ///
    /// The source carries a whole second only while the microsecond remainder
    /// is *strictly greater* than one million, so a remainder of exactly
    /// 1 000 000 stays unnormalised. That is preserved: it cannot change any
    /// reported duration, because one million microseconds contribute the same
    /// second either way, but it is observable through [`PartialEq`], and a
    /// port that silently normalised would make two watches compare equal that
    /// the source considers different.
    pub fn sum(self, other: Self) -> Self {
        let seconds = self.wall_seconds.saturating_add(other.wall_seconds);
        let micros = self.wall_micros.saturating_add(other.wall_micros);
        let (wall_seconds, wall_micros) = carry_excess_micros(seconds, micros);
        Self {
            wall_seconds,
            wall_micros,
            cpu: match (self.cpu, other.cpu) {
                (Some(left), Some(right)) => left.sum(right),
                _ => None,
            },
        }
    }
}

/// Wall, user and kernel timing of the current process, with pause and resume.
///
/// CPU time is the sum over all threads of the process. A watch must be started
/// before a reading means anything, but it need not be stopped: reading a
/// running watch adds the interval since the last start to the accumulated
/// total. [`stop`](Self::stop) and [`resume`](Self::resume) bracket the steps
/// that should not count towards the measurement.
///
/// Equality compares the accumulated interval, the last start reading and the
/// running flag, exactly as the source's `operator==` does. The source's
/// `operator<` compares something else entirely — only CPU time — so it is
/// deliberately *not* exposed as [`PartialOrd`], which would then disagree with
/// [`PartialEq`]; see [`cpu_time_cmp`](Self::cpu_time_cmp).
///
/// ```
/// use openms::system::stop_watch::StopWatch;
///
/// let mut watch = StopWatch::default();
/// assert_eq!(watch.clock_time(), 0.0);
/// watch.start()?;
/// assert!(watch.is_running());
/// watch.stop()?;
/// assert!(!watch.is_running());
/// assert!(watch.clock_time() >= 0.0);
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct StopWatch {
    accumulated: TimeSample,
    last_start: TimeSample,
    is_running: bool,
}

impl StopWatch {
    /// A cleared, stopped watch.
    pub fn new() -> Self {
        Self::default()
    }

    /// Start the stop watch.
    ///
    /// Data from previous measurements is discarded first, so `start`, `stop`,
    /// `start` cannot be used to resume — use [`resume`](Self::resume),
    /// [`stop`](Self::stop), [`resume`](Self::resume) for that.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the watch is already running, where
    /// the source throws `Exception::Precondition`.
    pub fn start(&mut self) -> Result<()> {
        if self.is_running {
            return Err(precondition("StopWatch is already started!"));
        }
        self.clear();
        self.last_start = TimeSample::now();
        self.is_running = true;
        Ok(())
    }

    /// Stop the stop watch; it can be resumed later.
    ///
    /// The interval since the last start is added to the accumulated total.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the watch is not running, where the
    /// source throws `Exception::Precondition`.
    pub fn stop(&mut self) -> Result<()> {
        if !self.is_running {
            return Err(precondition("StopWatch cannot be stopped if not running!"));
        }
        let difference = TimeSample::now().difference(self.last_start);
        self.accumulated = self.accumulated.sum(difference);
        self.is_running = false;
        Ok(())
    }

    /// Resume a stopped stop watch, keeping what it has already accumulated.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the watch is already running, where
    /// the source throws `Exception::Precondition`.
    pub fn resume(&mut self) -> Result<()> {
        if self.is_running {
            return Err(precondition(
                "StopWatch cannot be resumed if already running!",
            ));
        }
        self.last_start = TimeSample::now();
        self.is_running = true;
        Ok(())
    }

    /// Set the accumulated time to zero but keep running if the watch is running.
    ///
    /// A stopped watch is left stopped and cleared, which is
    /// [`clear`](Self::clear); a running watch is cleared and restarted, so its
    /// next reading measures from this call. Unlike the source's `reset`, which
    /// routes through a `start` that can throw, this cannot fail: the clear
    /// always precedes the restart, so the precondition is never violated.
    pub fn reset(&mut self) {
        let was_running = self.is_running;
        self.clear();
        if was_running {
            self.last_start = TimeSample::now();
            self.is_running = true;
        }
    }

    /// Set the accumulated time to zero and stop the watch.
    pub fn clear(&mut self) {
        *self = Self::default();
    }

    /// Accumulated wall-clock (real) time, in seconds.
    ///
    /// A watch that has never been started reports `0.0`, as the source does; a
    /// stopped watch reports only its accumulated total; a running watch adds
    /// the interval since its last start, each converted to `f64` separately
    /// and then summed, in the source's order.
    ///
    /// Unlike the source, which samples `gettimeofday`, this is monotonic: a
    /// backwards step of the system clock cannot shorten a measured interval.
    pub fn clock_time(&self) -> f64 {
        match self.running_interval() {
            None => self.accumulated.clock_time(),
            Some(interval) => self.accumulated.clock_time() + interval.clock_time(),
        }
    }

    /// Accumulated user time in seconds, summed over all threads.
    ///
    /// `None` means the platform does not report the user/kernel split, not
    /// that no user time was spent. Linux reports it; other platforms do not.
    pub fn user_time(&self) -> Option<f64> {
        match self.running_interval() {
            None => self.accumulated.user_time(),
            Some(interval) => Some(self.accumulated.user_time()? + interval.user_time()?),
        }
    }

    /// Accumulated system (kernel) time in seconds, summed over all threads.
    ///
    /// `None` means the platform does not report the user/kernel split, not
    /// that no system time was spent. Linux reports it; other platforms do not.
    pub fn system_time(&self) -> Option<f64> {
        match self.running_interval() {
            None => self.accumulated.kernel_time(),
            Some(interval) => Some(self.accumulated.kernel_time()? + interval.kernel_time()?),
        }
    }

    /// Accumulated CPU time in seconds: user time plus system time.
    ///
    /// Where the split is available this is literally
    /// [`user_time`](Self::user_time) `+` [`system_time`](Self::system_time),
    /// the source's `getCPUTime`, with the same two clock samples and the same
    /// addition order. Where it is not, the platform's aggregate process CPU
    /// clock answers instead, so a caller that only needs a total still gets
    /// one; `None` then means no process CPU clock at all.
    pub fn cpu_time(&self) -> Option<f64> {
        match (self.user_time(), self.system_time()) {
            (Some(user), Some(system)) => Some(user + system),
            _ => match self.running_interval() {
                None => self.accumulated.cpu_time(),
                Some(interval) => Some(self.accumulated.cpu_time()? + interval.cpu_time()?),
            },
        }
    }

    /// Whether the watch is currently running.
    pub fn is_running(&self) -> bool {
        self.is_running
    }

    /// Order two watches by CPU time, the source's `operator<` comparison.
    ///
    /// The source documents `<`, `<=`, `>` and `>=` as comparing clock, user
    /// *and* system time, but implements all four in terms of `getCPUTime()`
    /// alone; the implementation is what this reproduces, and the four
    /// operators map onto this one result — `a < b` is `Ordering::Less`,
    /// `a <= b` is "not `Greater`", and so on.
    ///
    /// Returns `None` when either watch cannot report CPU time on this
    /// platform. Because both operands are finite, a `Some` result is always a
    /// total order.
    pub fn cpu_time_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.cpu_time()?.partial_cmp(&other.cpu_time()?)
    }

    /// A compact report of wall, CPU, system and user time.
    ///
    /// Reproduces the source's `toString()`, for example
    /// `2.10 s (wall), 1.67 s (CPU), 0.12 s (system), 1.54 s (user)`. A
    /// component the platform cannot report renders as `n/a`, where the source
    /// — which always has a value — cannot produce that form.
    pub fn summary(&self) -> String {
        format!(
            "{} (wall), {} (CPU), {} (system), {} (user)",
            render(Some(self.clock_time())),
            render(self.cpu_time()),
            render(self.system_time()),
            render(self.user_time())
        )
    }

    /// Format a duration using only the units it needs.
    ///
    /// Below a minute the value is rendered in seconds with two decimals
    /// (`0.00 s`, `1.50 s`); from a minute, as zero-padded `MM:SS m`; from an
    /// hour, as `HH:MM:SS h`; from a day, as `Dd HH:MM:SS h` with an unpadded
    /// day count. Seconds are truncated, not rounded, before the split, so
    /// `100.5` renders as `01:40 m`. A negative duration takes the seconds form
    /// unchanged, because none of the three unit branches can then be positive.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `seconds` is not finite or exceeds
    /// [`MAX_FORMATTABLE_SECONDS`] in magnitude. The source casts its argument
    /// to an integer type with no check, which C++ leaves undefined for exactly
    /// those inputs.
    ///
    /// ```
    /// use openms::system::stop_watch::StopWatch;
    /// assert_eq!(StopWatch::format_seconds(0.0)?, "0.00 s");
    /// assert_eq!(StopWatch::format_seconds(100.5)?, "01:40 m");
    /// assert_eq!(StopWatch::format_seconds(3600.0 * 23.0 + 160.5)?, "23:02:40 h");
    /// # Ok::<(), openms::Error>(())
    /// ```
    pub fn format_seconds(seconds: f64) -> Result<String> {
        if !seconds.is_finite() {
            return Err(precondition("duration to format is not finite"));
        }
        if seconds.abs() > MAX_FORMATTABLE_SECONDS {
            return Err(precondition("duration to format exceeds the checked range"));
        }
        // The source truncates towards zero with a C cast; `as` matches that for
        // every value the range check above admits.
        let mut remaining = seconds as i64;
        let days = remaining / (3600 * 24);
        remaining -= days * (3600 * 24);
        let hours = remaining / 3600;
        remaining -= hours * 3600;
        let minutes = remaining / 60;
        remaining -= minutes * 60;
        let whole_seconds = remaining;
        Ok(if days > 0 {
            format!("{days}d {hours:02}:{minutes:02}:{whole_seconds:02} h")
        } else if hours > 0 {
            format!("{hours:02}:{minutes:02}:{whole_seconds:02} h")
        } else if minutes > 0 {
            format!("{minutes:02}:{whole_seconds:02} m")
        } else {
            format!("{seconds:.2} s")
        })
    }

    /// The interval since the last start, or `None` when the watch is stopped.
    fn running_interval(&self) -> Option<TimeSample> {
        self.is_running
            .then(|| TimeSample::now().difference(self.last_start))
    }
}

fn render(seconds: Option<f64>) -> String {
    seconds
        .and_then(|seconds| StopWatch::format_seconds(seconds).ok())
        .unwrap_or_else(|| "n/a".to_string())
}

fn precondition(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

fn micros_to_seconds(micros: i64) -> f64 {
    micros as f64 / 1e6
}

/// Monotonic origin for absolute wall readings, fixed at first use.
fn origin() -> Instant {
    static ORIGIN: OnceLock<Instant> = OnceLock::new();
    *ORIGIN.get_or_init(Instant::now)
}

/// The source's `while (start_time_usec > 1000000L)` carry, in constant time.
fn carry_excess_micros(seconds: i64, micros: i64) -> (i64, i64) {
    if micros > MICROS_PER_SECOND {
        let carry = (micros - 1) / MICROS_PER_SECOND;
        (
            seconds.saturating_add(carry),
            micros - carry * MICROS_PER_SECOND,
        )
    } else {
        (seconds, micros)
    }
}

/// The source's `while (diff.start_time_usec < 0L)` borrow, in constant time.
fn borrow_negative_micros(seconds: i64, micros: i64) -> (i64, i64) {
    if micros < 0 {
        let borrow = (MICROS_PER_SECOND - 1).saturating_sub(micros) / MICROS_PER_SECOND;
        (
            seconds.saturating_sub(borrow),
            micros.saturating_add(borrow.saturating_mul(MICROS_PER_SECOND)),
        )
    } else {
        (seconds, micros)
    }
}

#[cfg(target_os = "linux")]
fn zero_cpu() -> Option<CpuSample> {
    Some(CpuSample::Split {
        user_micros: 0,
        kernel_micros: 0,
    })
}

#[cfg(all(not(target_os = "linux"), any(unix, windows)))]
fn zero_cpu() -> Option<CpuSample> {
    Some(CpuSample::Total { micros: 0 })
}

#[cfg(not(any(unix, windows)))]
fn zero_cpu() -> Option<CpuSample> {
    None
}

/// Read `utime` and `stime` from `/proc/self/stat`.
///
/// `/proc/self` names the thread group, so both values already aggregate over
/// every thread of the process, which is what the source's `times` reports. The
/// second field is the executable name in parentheses and may itself contain
/// spaces and parentheses, so the scan starts after the last `)`.
#[cfg(target_os = "linux")]
fn sample_cpu() -> Option<CpuSample> {
    use std::io::Read;
    let mut text = String::new();
    std::fs::File::open("/proc/self/stat")
        .ok()?
        .take(MAX_PROC_BYTES)
        .read_to_string(&mut text)
        .ok()?;
    // After the command field the remainder starts at field 3 (`state`), so
    // `utime` (field 14) is the twelfth entry and `stime` (field 15) the next.
    let mut fields = text.rsplit_once(')')?.1.split_ascii_whitespace();
    let user: i64 = fields.nth(11)?.parse().ok()?;
    let kernel: i64 = fields.next()?.parse().ok()?;
    Some(CpuSample::Split {
        user_micros: user.checked_mul(MICROS_PER_TICK)?,
        kernel_micros: kernel.checked_mul(MICROS_PER_TICK)?,
    })
}

#[cfg(all(not(target_os = "linux"), any(unix, windows)))]
fn sample_cpu() -> Option<CpuSample> {
    let elapsed = cpu_time::ProcessTime::try_now().ok()?.as_duration();
    i64::try_from(elapsed.as_micros())
        .ok()
        .map(|micros| CpuSample::Total { micros })
}

#[cfg(not(any(unix, windows)))]
fn sample_cpu() -> Option<CpuSample> {
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(wall_seconds: i64, wall_micros: i64) -> TimeSample {
        TimeSample {
            wall_seconds,
            wall_micros,
            cpu: Some(CpuSample::Split {
                user_micros: 0,
                kernel_micros: 0,
            }),
        }
    }

    #[test]
    fn difference_borrows_one_second_per_negative_microsecond_remainder() {
        let later = sample(10, 250_000);
        let earlier = sample(9, 750_000);
        let difference = later.difference(earlier);
        assert_eq!(difference.wall_seconds, 0);
        assert_eq!(difference.wall_micros, 500_000);
        assert!((difference.clock_time() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn sum_leaves_exactly_one_million_microseconds_unnormalised() {
        // The source's carry loop tests `> 1000000`, not `>=`.
        let total = sample(1, 400_000).sum(sample(1, 600_000));
        assert_eq!((total.wall_seconds, total.wall_micros), (2, 1_000_000));
        assert!((total.clock_time() - 3.0).abs() < 1e-12);
        let carried = sample(1, 400_001).sum(sample(1, 600_000));
        assert_eq!((carried.wall_seconds, carried.wall_micros), (3, 1));
    }

    #[test]
    fn mixed_cpu_shapes_report_unavailable_rather_than_a_guess() {
        let split = sample(0, 0);
        let total = TimeSample {
            wall_seconds: 0,
            wall_micros: 0,
            cpu: Some(CpuSample::Total { micros: 5 }),
        };
        assert!(total.difference(split).cpu.is_none());
        assert!(split.sum(total).cpu.is_none());
    }

    #[test]
    fn cpu_sample_accessors_follow_the_reported_shape() {
        let split = CpuSample::Split {
            user_micros: 7,
            kernel_micros: 3,
        };
        assert_eq!(split.user_micros(), Some(7));
        assert_eq!(split.kernel_micros(), Some(3));
        assert_eq!(split.total_micros(), 10);
        let total = CpuSample::Total { micros: 10 };
        assert_eq!(total.user_micros(), None);
        assert_eq!(total.kernel_micros(), None);
        assert_eq!(total.total_micros(), 10);
    }
}
