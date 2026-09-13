// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `SYSTEM/StopWatch.h`: start/stop/resume state and duration formatting.
//!
//! The class test's 21 sections are covered here. Its timing sections assert
//! ranges rather than values, because neither it nor this port controls how
//! much CPU the host grants; the state machine, the frozen-while-stopped
//! invariant and the whole of `toString` are exact.

use openms::system::stop_watch::{CpuSample, MAX_FORMATTABLE_SECONDS, StopWatch, TimeSample};
use std::cmp::Ordering;
use std::thread::sleep;
use std::time::Duration;

/// The source class test busy-waits; sleeping costs no CPU and is enough for
/// every wall-clock assertion made here.
fn pause(millis: u64) {
    sleep(Duration::from_millis(millis));
}

/// Class-test section `StopWatch()`.
#[test]
fn a_fresh_watch_is_stopped_and_reports_zero() {
    let watch = StopWatch::new();
    assert!(!watch.is_running());
    assert_eq!(watch.clock_time(), 0.0);
    assert_eq!(watch, StopWatch::default());
    // Before the first start every reading is zero, not an error and not a
    // missing value: the platform can measure, there is simply nothing yet.
    assert_eq!(watch.cpu_time(), Some(0.0));
    #[cfg(target_os = "linux")]
    {
        assert_eq!(watch.user_time(), Some(0.0));
        assert_eq!(watch.system_time(), Some(0.0));
    }
}

/// Class-test section `StopWatch& operator=(const StopWatch& stop_watch)`.
///
/// The source test only says "tested below"; assignment there is a whole-object
/// copy, and the same holds here, where the type is [`Copy`].
#[test]
fn assignment_copies_the_complete_state() {
    let mut source = StopWatch::new();
    source.start().unwrap();
    pause(10);
    source.stop().unwrap();

    let assigned = source;
    assert_eq!(assigned, source);
    assert_eq!(assigned.clock_time(), source.clock_time());
    assert!(assigned.clock_time() > 0.0);

    // A copy is independent: restarting one does not disturb the other.
    let mut restarted = assigned;
    restarted.start().unwrap();
    assert!(restarted.is_running());
    assert!(!assigned.is_running());
}

/// Class-test section `StopWatch(const StopWatch& stop_watch)`.
///
/// Reproduces the source sequence: a running watch differs from a default one,
/// so does a stopped watch that accumulated time, a copy is equal, and `reset`
/// on a stopped watch returns it to the default value.
#[test]
fn copying_preserves_equality_and_reset_restores_the_default() {
    let mut first = StopWatch::new();
    let second = StopWatch::new();
    first.start().unwrap();
    pause(10);
    assert!(first != second); // running differs from stopped
    first.stop().unwrap();
    assert!(first != second); // accumulated time differs from none

    let assigned = first;
    assert_eq!(first, assigned);
    let copy_constructed = assigned;
    assert_eq!(first, copy_constructed);

    first.reset();
    assert_eq!(first, StopWatch::default());
}

/// Class-test section `bool isRunning() const`.
#[test]
fn is_running_follows_start_and_stop() {
    let mut watch = StopWatch::new();
    assert!(!watch.is_running());
    watch.start().unwrap();
    assert!(watch.is_running());
    watch.stop().unwrap();
    assert!(!watch.is_running());
    watch.resume().unwrap();
    assert!(watch.is_running());
}

/// Class-test section `bool operator==(const StopWatch& stop_watch) const`.
///
/// Equality compares the accumulated interval, the last start reading and the
/// running flag — not the elapsed times two running watches would report.
#[test]
fn equality_compares_state_not_elapsed_time() {
    let mut running = StopWatch::new();
    running.start().unwrap();
    let copy = running;
    assert_eq!(running, copy);

    let mut other = StopWatch::new();
    pause(5);
    other.start().unwrap();
    // Two watches started at different moments hold different start readings.
    assert!(running != other);
}

/// Class-test section `bool operator!=(const StopWatch& stop_watch) const`.
#[test]
fn inequality_is_the_negation_of_equality() {
    let mut watch = StopWatch::new();
    let fresh = StopWatch::new();
    assert_eq!(watch, fresh);
    watch.start().unwrap();
    assert!(watch != fresh);
    watch.clear();
    assert_eq!(watch, fresh);
}

/// Class-test section `bool operator<(const StopWatch& stop_watch) const`.
///
/// The source calls this untestable because it does not control system time.
/// The relation itself is testable: it compares CPU time alone, so two watches
/// that have been cleared compare equal however they got there.
#[test]
fn ordering_compares_cpu_time_only() {
    let fresh = StopWatch::new();
    let cleared = {
        let mut watch = StopWatch::new();
        watch.start().unwrap();
        pause(5);
        watch.stop().unwrap();
        watch.clear();
        watch
    };
    // Different histories, both back to zero CPU time.
    assert_eq!(fresh.cpu_time_cmp(&cleared), Some(Ordering::Equal));
    assert_eq!(cleared.cpu_time_cmp(&fresh), Some(Ordering::Equal));
}

/// Class-test section `bool operator<=(const StopWatch& stop_watch) const`.
#[test]
fn less_or_equal_is_not_greater() {
    let first = StopWatch::new();
    let second = StopWatch::new();
    assert_ne!(first.cpu_time_cmp(&second), Some(Ordering::Greater));
    assert_ne!(second.cpu_time_cmp(&first), Some(Ordering::Greater));
}

/// Class-test section `bool operator>=(const StopWatch& stop_watch) const`.
#[test]
fn greater_or_equal_is_not_less() {
    let first = StopWatch::new();
    let second = StopWatch::new();
    assert_ne!(first.cpu_time_cmp(&second), Some(Ordering::Less));
    assert_ne!(second.cpu_time_cmp(&first), Some(Ordering::Less));
}

/// Class-test section `bool operator>(const StopWatch& stop_watch) const`.
///
/// A watch that only ran while the process slept cannot be strictly ahead on
/// CPU time, which is the one thing the relation compares.
#[test]
fn greater_needs_more_cpu_time_not_more_wall_time() {
    let mut waited = StopWatch::new();
    waited.start().unwrap();
    pause(30);
    waited.stop().unwrap();
    let fresh = StopWatch::new();

    assert!(waited.clock_time() > fresh.clock_time());
    assert_ne!(fresh.cpu_time_cmp(&waited), Some(Ordering::Greater));
}

/// Class-test section `bool start()`.
///
/// Starting twice is refused, and a start discards what was accumulated before
/// it — which is exactly how it differs from `resume`.
#[test]
fn starting_twice_is_refused_and_a_start_discards_earlier_data() {
    let mut watch = StopWatch::new();
    watch.start().unwrap();
    assert!(watch.start().is_err());
    pause(200);
    watch.stop().unwrap();
    let accumulated = watch.clock_time();
    assert!(accumulated >= 0.2);

    watch.start().unwrap();
    watch.stop().unwrap();
    assert!(watch.clock_time() < accumulated); // the earlier interval is gone

    let mut resumed = StopWatch::new();
    resumed.resume().unwrap();
    pause(20);
    resumed.stop().unwrap();
    let first_interval = resumed.clock_time();
    resumed.resume().unwrap();
    pause(20);
    resumed.stop().unwrap();
    assert!(resumed.clock_time() > first_interval); // resume keeps it
}

/// Class-test section `bool stop()`.
///
/// The source section carries most of the class's behaviour: stopping twice is
/// refused, a stopped watch is frozen, `reset` keeps a running watch running,
/// `resume` keeps accumulating, and a watch that never stopped stays ahead.
#[test]
fn stopping_freezes_the_reading_and_resume_keeps_accumulating() {
    let wait = 200;
    let wait_more = 100;
    let mut stopped = StopWatch::new();
    let mut never_stopped = StopWatch::new();
    let mut restarted = StopWatch::new();
    let mut resumed = StopWatch::new();
    stopped.start().unwrap();
    never_stopped.start().unwrap();
    restarted.start().unwrap();
    resumed.resume().unwrap();
    pause(wait);
    stopped.stop().unwrap();
    resumed.stop().unwrap();
    assert!(stopped.stop().is_err()); // cannot stop twice

    assert!(stopped.clock_time() > 0.1);

    let cpu = stopped.cpu_time();
    let clock = stopped.clock_time();
    let system = stopped.system_time();
    let user = stopped.user_time();
    restarted.reset();
    assert!(restarted.is_running()); // reset keeps it running
    resumed.resume().unwrap();
    pause(wait_more);

    // A stopped watch does not move.
    assert_eq!(stopped.cpu_time(), cpu);
    assert_eq!(stopped.clock_time(), clock);
    assert_eq!(stopped.system_time(), system);
    assert_eq!(stopped.user_time(), user);

    assert!(stopped.clock_time() > (wait as f64 / 1000.0) * 0.95);

    // The watch that never stopped is ahead on wall time and not behind on CPU.
    assert!(stopped.clock_time() < never_stopped.clock_time());
    assert_ne!(
        stopped.cpu_time_cmp(&never_stopped),
        Some(Ordering::Greater)
    );

    stopped.reset(); // was stopped, so stays stopped
    assert!(!stopped.is_running());
    assert_eq!(stopped, StopWatch::default());

    // Kept running across the reset, so it accumulated again.
    assert!(restarted.clock_time() > 0.0);

    // Never stopped after the second resume: queried on the fly.
    assert!(resumed.clock_time() > ((wait + wait_more) as f64 / 1000.0) * 0.95);
}

/// Class-test section `void clear()`.
#[test]
fn clear_stops_and_zeroes_the_watch() {
    let mut watch = StopWatch::new();
    watch.start().unwrap();
    pause(10);
    watch.clear();
    assert!(!watch.is_running());
    assert_eq!(watch, StopWatch::default());
    assert_eq!(watch.clock_time(), 0.0);
}

/// Class-test section `void reset()`.
///
/// The source calls this "done above"; both branches are asserted here, because
/// whether the watch keeps running is the one thing `reset` decides.
#[test]
fn reset_keeps_a_running_watch_running_and_a_stopped_one_stopped() {
    let mut kept = StopWatch::new();
    let mut restarted = StopWatch::new();
    kept.start().unwrap();
    restarted.start().unwrap();
    pause(50);
    restarted.reset();
    assert!(restarted.is_running());
    pause(10);
    // Both ran for the same wall time, but the reset one measures from the reset.
    assert!(restarted.clock_time() < kept.clock_time());

    let mut stopped = StopWatch::new();
    stopped.start().unwrap();
    pause(10);
    stopped.stop().unwrap();
    stopped.reset();
    assert!(!stopped.is_running());
    assert_eq!(stopped.clock_time(), 0.0);
}

/// Class-test section `void resume()`.
#[test]
fn resuming_a_running_watch_is_refused() {
    let mut watch = StopWatch::new();
    watch.start().unwrap();
    assert!(watch.resume().is_err());
    watch.stop().unwrap();
    assert!(watch.resume().is_ok());
}

/// Class-test section `double getClockTime() const`.
#[test]
fn clock_time_is_zero_before_the_first_start_and_monotonic_afterwards() {
    let mut watch = StopWatch::new();
    assert_eq!(watch.clock_time(), 0.0);
    watch.start().unwrap();
    let first = watch.clock_time();
    pause(20);
    let second = watch.clock_time();
    assert!(second >= first);
    assert!(second >= 0.015);
}

/// Class-test section `double getUserTime() const`.
#[test]
fn user_time_is_reported_or_explicitly_unavailable() {
    let mut watch = StopWatch::new();
    watch.start().unwrap();
    pause(10);
    watch.stop().unwrap();
    match watch.user_time() {
        Some(seconds) => assert!(seconds >= 0.0),
        // Only a platform without the user/kernel split may answer this way.
        None => assert_eq!(watch.system_time(), None),
    }
    #[cfg(target_os = "linux")]
    assert!(watch.user_time().is_some());
}

/// Class-test section `double getSystemTime() const`.
#[test]
fn system_time_is_reported_or_explicitly_unavailable() {
    let mut watch = StopWatch::new();
    watch.start().unwrap();
    pause(10);
    watch.stop().unwrap();
    match watch.system_time() {
        Some(seconds) => assert!(seconds >= 0.0),
        None => assert_eq!(watch.user_time(), None),
    }
    #[cfg(target_os = "linux")]
    assert!(watch.system_time().is_some());
}

/// Class-test section `double getCPUTime() const`.
///
/// Where the split is available, CPU time is exactly user plus system, which is
/// the source's own definition rather than an independent measurement.
#[test]
fn cpu_time_is_user_plus_system_where_the_split_exists() {
    let mut watch = StopWatch::new();
    watch.start().unwrap();
    pause(10);
    watch.stop().unwrap();
    match (watch.user_time(), watch.system_time()) {
        (Some(user), Some(system)) => assert_eq!(watch.cpu_time(), Some(user + system)),
        // No split: an aggregate is still reported where a CPU clock exists.
        _ => assert!(watch.cpu_time().is_some() || watch.user_time().is_none()),
    }
}

/// Class-test section `virtual ~StopWatch()`.
///
/// The source's destructor is implicit and releases nothing; the port's type
/// owns no resources either, so it is [`Copy`] and a copy outlives its origin's
/// scope unchanged.
#[test]
fn a_watch_owns_no_resources_and_survives_its_origin_scope() {
    let escaped = {
        let mut watch = StopWatch::new();
        watch.start().unwrap();
        watch
    };
    assert!(escaped.is_running());
}

/// Class-test section `static std::string toString(double time)`.
///
/// Every literal is the source test's own expectation, transcribed.
#[test]
fn format_seconds_reproduces_every_class_test_literal() {
    assert_eq!(StopWatch::format_seconds(0.0).unwrap(), "0.00 s");
    assert_eq!(StopWatch::format_seconds(1.0).unwrap(), "1.00 s");
    assert_eq!(StopWatch::format_seconds(1.5).unwrap(), "1.50 s");
    assert_eq!(StopWatch::format_seconds(100.5).unwrap(), "01:40 m");
    assert_eq!(
        StopWatch::format_seconds(3600.0 * 24.0 * 5.0 + 3600.0 * 9.0 + 5.0).unwrap(),
        "5d 09:00:05 h"
    );
    assert_eq!(StopWatch::format_seconds(160.5).unwrap(), "02:40 m");
    assert_eq!(
        StopWatch::format_seconds(3600.0 * 23.0 + 160.5).unwrap(),
        "23:02:40 h"
    );
}

/// Boundaries around the unit thresholds, derived from the source expression.
///
/// The seconds branch keeps two decimals of the *unrounded* argument, while the
/// minute, hour and day branches truncate to whole seconds first, so 59.999 s
/// prints as `60.00 s` and 60.0 s as the first minute form.
#[test]
fn format_seconds_switches_units_exactly_where_the_source_does() {
    assert_eq!(StopWatch::format_seconds(59.999).unwrap(), "60.00 s");
    assert_eq!(StopWatch::format_seconds(60.0).unwrap(), "01:00 m");
    assert_eq!(StopWatch::format_seconds(3599.0).unwrap(), "59:59 m");
    assert_eq!(StopWatch::format_seconds(3600.0).unwrap(), "01:00:00 h");
    assert_eq!(StopWatch::format_seconds(86399.0).unwrap(), "23:59:59 h");
    assert_eq!(StopWatch::format_seconds(86400.0).unwrap(), "1d 00:00:00 h");
    assert_eq!(
        StopWatch::format_seconds(100.0 * 86400.0).unwrap(),
        "100d 00:00:00 h"
    );
    // Negative durations cannot satisfy any unit branch and stay in seconds.
    assert_eq!(StopWatch::format_seconds(-100.5).unwrap(), "-100.50 s");
    assert_eq!(StopWatch::format_seconds(-0.004).unwrap(), "-0.00 s");
}

/// The range check the source lacks.
#[test]
fn format_seconds_refuses_what_the_source_casts_unchecked() {
    assert!(StopWatch::format_seconds(f64::NAN).is_err());
    assert!(StopWatch::format_seconds(f64::INFINITY).is_err());
    assert!(StopWatch::format_seconds(f64::NEG_INFINITY).is_err());
    assert!(StopWatch::format_seconds(MAX_FORMATTABLE_SECONDS * 2.0).is_err());
    assert!(StopWatch::format_seconds(MAX_FORMATTABLE_SECONDS).is_ok());
}

/// The instance `toString()`: the four components in the source's order.
#[test]
fn summary_lists_wall_cpu_system_and_user_in_the_source_order() {
    let watch = StopWatch::new();
    let summary = watch.summary();
    assert!(summary.starts_with("0.00 s (wall), "));
    assert!(summary.contains(" (CPU), "));
    assert!(summary.contains(" (system), "));
    assert!(summary.ends_with(" (user)"));
    #[cfg(target_os = "linux")]
    assert_eq!(
        summary,
        "0.00 s (wall), 0.00 s (CPU), 0.00 s (system), 0.00 s (user)"
    );
}

/// `TimeDiff_::operator-` and `operator+=`, exposed here as [`TimeSample`].
///
/// The carry loop tests `> 1000000`, so exactly one million microseconds stays
/// unnormalised; the borrow loop tests `< 0`. Both are reproduced, and the
/// difference is invisible in the reported duration but visible in the fields.
#[test]
fn time_sample_arithmetic_follows_the_source_normalisation() {
    let cpu = Some(CpuSample::Split {
        user_micros: 0,
        kernel_micros: 0,
    });
    let sample = |wall_seconds, wall_micros| TimeSample {
        wall_seconds,
        wall_micros,
        cpu,
    };

    let difference = sample(10, 250_000).difference(sample(9, 750_000));
    assert_eq!(difference.wall_seconds, 0);
    assert_eq!(difference.wall_micros, 500_000);
    assert_eq!(difference.clock_time(), 0.5);

    let exact = sample(1, 400_000).sum(sample(1, 600_000));
    assert_eq!((exact.wall_seconds, exact.wall_micros), (2, 1_000_000));
    assert_eq!(exact.clock_time(), 3.0);

    let carried = sample(1, 400_001).sum(sample(1, 600_000));
    assert_eq!((carried.wall_seconds, carried.wall_micros), (3, 1));

    // A reading taken now is never before the process origin.
    let now = TimeSample::now();
    assert!(now.wall_seconds >= 0);
    assert!(now.wall_micros >= 0 && now.wall_micros < 1_000_000);
}

/// A mixed pair of CPU shapes reports unavailable instead of a wrong number.
#[test]
fn mixing_split_and_aggregate_cpu_samples_yields_no_cpu_time() {
    let split = TimeSample {
        wall_seconds: 0,
        wall_micros: 0,
        cpu: Some(CpuSample::Split {
            user_micros: 10,
            kernel_micros: 4,
        }),
    };
    let total = TimeSample {
        wall_seconds: 0,
        wall_micros: 0,
        cpu: Some(CpuSample::Total { micros: 14 }),
    };
    assert_eq!(split.cpu_time(), Some(1.4e-5));
    assert_eq!(total.cpu_time(), Some(1.4e-5));
    assert_eq!(split.user_time(), Some(1.0e-5));
    assert_eq!(total.user_time(), None);
    assert_eq!(split.sum(total).cpu_time(), None);
    assert_eq!(total.difference(split).cpu_time(), None);
}
