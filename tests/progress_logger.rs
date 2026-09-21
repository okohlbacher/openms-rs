// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::concept::progress_logger::{
    CommandProgressLogger, MAX_PROGRESS_DEPTH, MAX_PROGRESS_LABEL_BYTES, ProgressBackend,
    ProgressClock, ProgressLogType, ProgressLogger, ProgressNesting, ProgressReporter,
    ProgressTime,
};
use openms::{Error, Result};
use std::collections::{BTreeMap, VecDeque};
use std::io::{self, Write};
use std::sync::{Arc, Mutex};

/// Captured OpenMS4 Release output, one row per call; see
/// `tests/data/progress_logger_provenance.json` (`release_oracle`).
const RELEASE_FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/data/progress_logger_release_range.tsv"
);

/// A cloneable writer, so a command backend owned by a `ProgressLogger` can be
/// read back per call.
#[derive(Clone, Default)]
struct SharedOutput(Arc<Mutex<Vec<u8>>>);
impl SharedOutput {
    fn take(&self) -> String {
        String::from_utf8(std::mem::take(&mut *self.0.lock().unwrap())).unwrap()
    }
}
impl Write for SharedOutput {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// The fixture's mask: only the two timing texts of each command summary
/// line become `<TIME>`, exactly as the oracle's `extract.py` masks them.
/// Whether `text` has one of the four shapes `StopWatch::toString(double)`
/// prints (StopWatch.cpp:231-234): `Nd HH:MM:SS h`, `HH:MM:SS h`, `MM:SS m`, or
/// `StringUtils::number(seconds, 2)` followed by ` s`. The mask accepts only
/// these, so text the Release build never prints (the port's
/// `unavailable (CPU)`) cannot pass as Release output.
fn is_release_time(text: &str) -> bool {
    let all_digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    let two_digit_fields = |s: &str, n: usize| {
        s.split(':').count() == n && s.split(':').all(|f| f.len() == 2 && all_digits(f))
    };
    if let Some(seconds) = text.strip_suffix(" s") {
        return match seconds.split_once('.') {
            Some((whole, fraction)) => all_digits(whole) && all_digits(fraction),
            None => all_digits(seconds),
        };
    }
    if let Some(minutes) = text.strip_suffix(" m") {
        return two_digit_fields(minutes, 2);
    }
    if let Some(hours) = text.strip_suffix(" h") {
        return match hours.split_once("d ") {
            Some((days, rest)) => all_digits(days) && two_digit_fields(rest, 3),
            None => two_digit_fields(hours, 3),
        };
    }
    false
}

fn mask_timing(text: &str) -> String {
    const HEAD: &str = "-- done [took ";
    let mut out = String::new();
    let mut rest = text;
    while let Some(start) = rest.find(HEAD) {
        let after = &rest[start + HEAD.len()..];
        let cpu = after.find(" (CPU), ").expect("CPU timing text");
        let wall_start = cpu + " (CPU), ".len();
        let wall = after[wall_start..]
            .find(" (Wall)")
            .expect("wall timing text");
        for timing in [&after[..cpu], &after[wall_start..wall_start + wall]] {
            assert!(
                is_release_time(timing),
                "{timing:?} is not a time the Release build prints"
            );
        }
        out.push_str(&rest[..start]);
        out.push_str("-- done [took <TIME> (CPU), <TIME> (Wall)");
        rest = &after[wall_start + wall + " (Wall)".len()..];
    }
    out.push_str(rest);
    out
}

fn unescape(field: &str) -> String {
    let mut out = String::new();
    let mut chars = field.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('r') => out.push('\r'),
            Some('t') => out.push('\t'),
            Some('\\') => out.push('\\'),
            other => panic!("bad escape {other:?} in fixture field {field:?}"),
        }
    }
    out
}

fn sample(second: i64, wall: f64, cpu: Option<f64>) -> ProgressTime {
    ProgressTime {
        wall_second: second,
        wall_seconds: wall,
        cpu_seconds: cpu,
    }
}
fn manual_clock(time: ProgressTime) -> (ProgressClock, Arc<Mutex<ProgressTime>>) {
    let state = Arc::new(Mutex::new(time));
    let clock_state = state.clone();
    (Arc::new(move || Ok(*clock_state.lock().unwrap())), state)
}
#[derive(Clone, Debug, PartialEq, Eq)]
enum Event {
    Start(i64, i64, String, usize),
    Set(i64, usize),
    Next(i64),
    End(usize, u64),
}
struct Recorder {
    events: Arc<Mutex<Vec<Event>>>,
    current: i64,
}
impl ProgressBackend for Recorder {
    fn start_progress(&mut self, begin: i64, end: i64, label: &str, depth: usize) -> Result<()> {
        self.events
            .lock()
            .unwrap()
            .push(Event::Start(begin, end, label.into(), depth));
        self.current = begin;
        Ok(())
    }
    fn set_progress(&mut self, value: i64, depth: usize) -> Result<()> {
        self.events.lock().unwrap().push(Event::Set(value, depth));
        Ok(())
    }
    fn next_progress(&mut self) -> Result<i64> {
        self.current += 1;
        self.events.lock().unwrap().push(Event::Next(self.current));
        Ok(self.current)
    }
    fn end_progress(&mut self, depth: usize, bytes: u64) -> Result<()> {
        self.events.lock().unwrap().push(Event::End(depth, bytes));
        Ok(())
    }
}
fn recorder(events: &Arc<Mutex<Vec<Event>>>) -> Box<dyn ProgressBackend> {
    Box::new(Recorder {
        events: events.clone(),
        current: 0,
    })
}

#[test]
fn source_mode_copy_and_disabled_calls() {
    let (clock, _) = manual_clock(sample(100, 0.0, None));
    let nesting = ProgressNesting::default();
    let mut logger = ProgressLogger::with_clock_and_nesting(clock, nesting.clone());
    assert_eq!(logger.log_type(), ProgressLogType::None);
    for kind in [
        ProgressLogType::Cmd,
        ProgressLogType::Gui,
        ProgressLogType::None,
    ] {
        logger.set_log_type(kind);
        assert_eq!(logger.log_type(), kind);
        let copy = logger.clone();
        assert_eq!(copy.log_type(), kind);
        let mut assigned = ProgressLogger::new();
        assigned.clone_from(&copy);
        assert_eq!(assigned.log_type(), kind);
    }
    logger.start_progress(0, 10, "disabled").unwrap();
    for value in [0, 5, 10] {
        logger.set_progress(value).unwrap();
    }
    logger.next_progress().unwrap();
    logger.end_progress(1024).unwrap();
    logger.end_progress(0).unwrap();
    assert_eq!(nesting.depth(), 0);
}

#[test]
fn increment_precedes_throttling_and_explicit_set_does_not_move_counter() {
    let events = Arc::new(Mutex::new(Vec::new()));
    let (clock, time) = manual_clock(sample(42, 0.0, None));
    let mut logger = ProgressLogger::with_clock_and_nesting(clock, ProgressNesting::default());
    logger.set_logger(recorder(&events));
    assert_eq!(logger.log_type(), ProgressLogType::None);
    logger.start_progress(10, 20, "range").unwrap();
    logger.set_progress(19).unwrap(); // same second is suppressed
    logger.next_progress().unwrap();
    logger.next_progress().unwrap();
    time.lock().unwrap().wall_second = 43;
    logger.set_progress(19).unwrap();
    logger.next_progress().unwrap(); // increment13; output remains suppressed
    time.lock().unwrap().wall_second = 41; // backward time is another bucket
    logger.next_progress().unwrap();
    logger.end_progress(2048).unwrap();
    assert_eq!(
        *events.lock().unwrap(),
        vec![
            Event::Start(10, 20, "range".into(), 0),
            Event::Next(11),
            Event::Next(12),
            Event::Set(19, 1),
            Event::Next(13),
            Event::Next(14),
            Event::Set(14, 1),
            Event::End(0, 2048),
        ]
    );
}

#[test]
fn two_clock_reads_on_dispatch_and_cloned_timestamp_are_observable() {
    let samples = Arc::new(Mutex::new(VecDeque::from([7, 8, 9, 9, 10, 11])));
    let used = samples.clone();
    let clock: ProgressClock = Arc::new(move || {
        let second = used
            .lock()
            .unwrap()
            .pop_front()
            .expect("unexpected clock read");
        Ok(sample(second, 0.0, None))
    });
    let events = Arc::new(Mutex::new(Vec::new()));
    let mut logger = ProgressLogger::with_clock_and_nesting(clock, ProgressNesting::default());
    logger.set_logger(recorder(&events));
    logger.start_progress(0, 5, "timing").unwrap(); // records7
    logger.set_progress(1).unwrap(); // compares8, records9
    let mut copy = logger.clone();
    copy.set_logger(recorder(&events));
    copy.set_progress(2).unwrap(); // equal9 suppresses
    copy.set_progress(3).unwrap(); // compares10, records11
    logger.end_progress(0).unwrap();
    assert!(samples.lock().unwrap().is_empty());
    assert_eq!(
        *events.lock().unwrap(),
        vec![
            Event::Start(0, 5, "timing".into(), 0),
            Event::Set(1, 1),
            Event::Set(3, 1),
            Event::End(0, 0)
        ]
    );
}

#[test]
fn nesting_is_shared_including_none_and_end_without_start_dispatches() {
    let events = Arc::new(Mutex::new(Vec::new()));
    let (clock, time) = manual_clock(sample(1, 0.0, None));
    let nesting = ProgressNesting::default();
    let mut outer = ProgressLogger::with_clock_and_nesting(clock.clone(), nesting.clone());
    let mut inner = ProgressLogger::with_clock_and_nesting(clock, nesting.clone());
    inner.set_logger(recorder(&events));
    outer.start_progress(0, 0, "none").unwrap();
    inner.start_progress(-2, 2, "inner").unwrap();
    time.lock().unwrap().wall_second = 2;
    inner.set_progress(0).unwrap();
    inner.end_progress(0).unwrap();
    outer.end_progress(0).unwrap();
    inner.end_progress(5).unwrap();
    assert_eq!(
        *events.lock().unwrap(),
        vec![
            Event::Start(-2, 2, "inner".into(), 1),
            Event::Set(0, 2),
            Event::End(1, 0),
            Event::End(0, 5)
        ]
    );
    assert_eq!(nesting.depth(), 0);
}

#[test]
fn gui_factory_reselection_copy_and_drop_follow_source_lifetimes() {
    let events = Arc::new(Mutex::new(Vec::new()));
    let created = Arc::new(Mutex::new(0));
    let (clock, time) = manual_clock(sample(100, 0.0, None));
    let nesting = ProgressNesting::default();
    let mut logger = ProgressLogger::with_clock_and_nesting(clock, nesting.clone());
    let factory_events = events.clone();
    let factory_created = created.clone();
    logger.set_gui_factory(Arc::new(move || {
        *factory_created.lock().unwrap() += 1;
        recorder(&factory_events)
    }));
    logger.set_log_type(ProgressLogType::Gui);
    logger.start_progress(8, 12, "active").unwrap();
    let mut copy = logger.clone();
    assert_eq!(*created.lock().unwrap(), 2);
    time.lock().unwrap().wall_second = 101;
    copy.next_progress().unwrap(); // fresh backend, current0 ->1
    copy.set_log_type(ProgressLogType::Gui); // equal type resets again
    copy.next_progress().unwrap(); // fresh current0 ->1, now suppressed
    assert_eq!(*created.lock().unwrap(), 3);
    logger.end_progress(0).unwrap();
    let length = events.lock().unwrap().len();
    drop(copy); // no implicit end
    assert_eq!(events.lock().unwrap().len(), length);
    assert_eq!(
        *events.lock().unwrap(),
        vec![
            Event::Start(8, 12, "active".into(), 0),
            Event::Next(1),
            Event::Set(1, 1),
            Event::Next(1),
            Event::End(0, 0)
        ]
    );
    assert_eq!(nesting.depth(), 0);
    logger.set_logger(recorder(&events));
    let mut none = ProgressLogger::with_clock_and_nesting(
        Arc::new(|| Ok(sample(9, 0.0, None))),
        ProgressNesting::default(),
    );
    none.set_logger(recorder(&events));
    let mut copy = none.clone(); // NONE mode discards even a custom recording backend
    let length = events.lock().unwrap().len();
    copy.next_progress().unwrap();
    assert_eq!(events.lock().unwrap().len(), length);
}

#[test]
fn command_output_preserves_f32_percentage_dots_and_range_diagnostics() {
    let (clock, time) = manual_clock(sample(1, 0.0, Some(3.0)));
    let mut cmd = CommandProgressLogger::with_clock(Vec::new(), clock);
    cmd.start_progress(-1, 2, "test", 0).unwrap();
    cmd.set_progress(0, 1).unwrap();
    cmd.set_progress(9, 1).unwrap();
    assert_eq!(cmd.next_progress().unwrap(), 0);
    assert_eq!(cmd.next_progress().unwrap(), 1);
    *time.lock().unwrap() = sample(3, 2.0, Some(3.25));
    cmd.end_progress(0, 2049).unwrap(); // double rate truncation precedes unit choice
    assert_eq!(
        String::from_utf8(cmd.into_inner()).unwrap(),
        "Progress of 'test':\n\r  33.33 %               ProgressLogger: Invalid progress value '9'. Should be between '-1' and '2'!\n\r-- done [took 0.25 s (CPU), 2.00 s (Wall) @ 1 KiB/s] -- \n"
    );
    let (clock, time) = manual_clock(sample(4, 1.0, None));
    let mut cmd = CommandProgressLogger::with_clock(Vec::new(), clock);
    cmd.start_progress(4, 4, "unknown", 1).unwrap();
    cmd.set_progress(i64::MIN, 2).unwrap(); // zero-span ignores value
    time.lock().unwrap().wall_seconds = 1.0;
    cmd.end_progress(1, 0).unwrap(); // zero wall is valid without throughput
    assert_eq!(
        String::from_utf8(cmd.into_inner()).unwrap(),
        "  Progress of 'unknown':\n.\r  -- done [took unavailable (CPU), 0.00 s (Wall)] -- \n"
    );
}

#[test]
fn command_duration_and_rate_boundaries_are_source_expression_oracles() {
    for (seconds, expected) in [
        (59.999, "60.00 s"),
        (60.0, "01:00 m"),
        (61.9, "01:01 m"),
        (3661.5, "01:01:01 h"),
        (90061.9, "1d 01:01:01 h"),
    ] {
        let (clock, time) = manual_clock(sample(0, 0.0, Some(0.0)));
        let mut cmd = CommandProgressLogger::with_clock(Vec::new(), clock);
        cmd.start_progress(0, 1, "duration", 0).unwrap();
        *time.lock().unwrap() = sample(1, seconds, Some(seconds));
        cmd.end_progress(0, 0).unwrap();
        let text = String::from_utf8(cmd.into_inner()).unwrap();
        assert!(
            text.contains(&format!("{expected} (CPU), {expected} (Wall)")),
            "{text}"
        );
    }
    for (bytes, expected) in [
        (1, "1 byte"),
        (1023, "1023 byte"),
        (1024, "1 KiB"),
        (1536, "1.5 KiB"),
        (1_048_575, "1024 KiB"),
        (1_048_576, "1 MiB"),
        (
            1_u64 << 60,
            "Congrats. That's a lot of bytes: 1152921504606846976",
        ),
    ] {
        let (clock, time) = manual_clock(sample(0, 0.0, None));
        let mut cmd = CommandProgressLogger::with_clock(Vec::new(), clock);
        cmd.start_progress(0, 1, "rate", 0).unwrap();
        time.lock().unwrap().wall_seconds = 1.0;
        cmd.end_progress(0, bytes).unwrap();
        let text = String::from_utf8(cmd.into_inner()).unwrap();
        assert!(text.contains(&format!(" @ {expected}/s")), "{text}");
    }
}

#[test]
fn checked_bounds_and_failures_leave_dispatch_state_usable() {
    let (clock, time) = manual_clock(sample(1, 0.0, None));
    let nesting = ProgressNesting::default();
    let events = Arc::new(Mutex::new(Vec::new()));
    let mut logger = ProgressLogger::with_clock_and_nesting(clock.clone(), nesting.clone());
    logger.set_logger(recorder(&events));
    assert!(
        logger
            .start_progress(0, 1, &"x".repeat(MAX_PROGRESS_LABEL_BYTES + 1))
            .is_err()
    );
    assert_eq!(nesting.depth(), 0);
    assert!(events.lock().unwrap().is_empty());

    // An inverted range is not an error. The source's only range check,
    // `OPENMS_PRECONDITION(begin <= end)` (ProgressLogger.cpp:235), is compiled
    // out of the reference Release build, which starts the section and hands
    // the range to its backend unchanged (:237).
    logger.start_progress(2, 1, "bad").unwrap();
    assert_eq!(nesting.depth(), 1);
    logger.end_progress(0).unwrap();
    assert_eq!(nesting.depth(), 0);
    assert_eq!(
        *events.lock().unwrap(),
        vec![Event::Start(2, 1, "bad".into(), 0), Event::End(0, 0)]
    );
    // Through the command backend, every value is out of the inverted range.
    // Each expected string is Release output, rows 15-20 (case `bad`) of
    // tests/data/progress_logger_release_range.tsv. The Release build always
    // has a CPU time, so this block's clock reports one too.
    time.lock().unwrap().cpu_seconds = Some(0.0);
    let output = SharedOutput::default();
    let mut logger = ProgressLogger::with_clock_and_nesting(clock.clone(), nesting.clone());
    logger.set_log_type(ProgressLogType::Cmd);
    logger.set_logger(Box::new(CommandProgressLogger::with_clock(
        output.clone(),
        clock.clone(),
    )));
    logger.start_progress(2, 1, "bad").unwrap();
    assert_eq!(nesting.depth(), 1);
    assert_eq!(output.take(), "Progress of 'bad':\n");
    for (value, expected) in [
        (
            0,
            "ProgressLogger: Invalid progress value '0'. Should be between '2' and '1'!\n",
        ),
        (
            1,
            "ProgressLogger: Invalid progress value '1'. Should be between '2' and '1'!\n",
        ),
        (
            2,
            "ProgressLogger: Invalid progress value '2'. Should be between '2' and '1'!\n",
        ),
        (
            3,
            "ProgressLogger: Invalid progress value '3'. Should be between '2' and '1'!\n",
        ),
    ] {
        time.lock().unwrap().wall_second += 1; // the driver waited for a new second
        logger.set_progress(value).unwrap();
        assert_eq!(output.take(), expected, "set_progress({value})");
    }
    logger.end_progress(0).unwrap();
    assert_eq!(nesting.depth(), 0);
    assert_eq!(
        mask_timing(&output.take()),
        "\r-- done [took <TIME> (CPU), <TIME> (Wall)] -- \n"
    );
    *time.lock().unwrap() = sample(1, 0.0, None);

    // End without start: the Release build's StopWatch::stop throws
    // unconditionally (StopWatch.cpp:55; fixture row 50), not a Debug check.
    let mut cmd = CommandProgressLogger::with_clock(Vec::new(), clock);
    assert!(cmd.end_progress(0, 0).is_err());
    cmd.start_progress(i64::MAX, i64::MAX, "counter", 0)
        .unwrap();
    // Native bounds, not source checks: the source's `++current_` and the
    // differences below are undefined signed overflow, and the depth bound
    // limits indentation.
    assert!(cmd.next_progress().is_err());
    assert!(cmd.next_progress().is_err()); // error does not wrap counter
    assert!(cmd.set_progress(0, MAX_PROGRESS_DEPTH + 1).is_err());
    cmd.end_progress(0, 0).unwrap();
    cmd.start_progress(i64::MIN, i64::MAX, "arithmetic", 0)
        .unwrap();
    assert!(cmd.set_progress(0, 0).is_err());
    assert!(cmd.end_progress(0, 1).is_err()); // division by zero
    time.lock().unwrap().wall_seconds = 1.0;
    assert!(cmd.end_progress(0, u64::MAX).is_err()); // rounded f64 becomes2^64
    cmd.end_progress(0, 1).unwrap(); // rejected end remains retryable
    assert!(cmd.end_progress(0, 0).is_err());
}

#[test]
fn nesting_limit_is_checked_without_wrapping_and_none_can_unwind() {
    let (clock, _) = manual_clock(sample(1, 0.0, None));
    let nesting = ProgressNesting::default();
    let mut logger = ProgressLogger::with_clock_and_nesting(clock, nesting.clone());
    for _ in 0..MAX_PROGRESS_DEPTH {
        logger.start_progress(0, 0, "").unwrap();
    }
    // The zero-width range is valid (it prints dots); this start fails only on
    // the native MAX_PROGRESS_DEPTH bound, which has no source counterpart.
    assert!(logger.start_progress(0, 0, "").is_err());
    assert_eq!(nesting.depth(), MAX_PROGRESS_DEPTH);
    for _ in 0..MAX_PROGRESS_DEPTH {
        logger.end_progress(0).unwrap();
    }
    assert_eq!(nesting.depth(), 0);
}

#[test]
fn writer_clock_and_backend_errors_propagate() {
    struct Broken;
    impl Write for Broken {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            Err(io::Error::other("writer failed"))
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }
    let (clock, time) = manual_clock(sample(1, 0.0, None));
    let mut cmd = CommandProgressLogger::with_clock(Broken, clock.clone());
    assert!(matches!(
        cmd.start_progress(0, 1, "io", 0),
        Err(Error::Io(_))
    ));
    // Injected clocks are native: a NaN sample, a failing clock and a failing
    // backend have no source counterpart, and all use valid ranges.
    let mut cmd = CommandProgressLogger::with_clock(Vec::new(), clock);
    time.lock().unwrap().wall_seconds = f64::NAN;
    assert!(cmd.start_progress(0, 1, "bad clock", 0).is_err());
    assert_eq!(cmd.into_inner(), b"Progress of 'bad clock':\n");
    let bad_clock: ProgressClock = Arc::new(|| Err(Error::InvalidValue("clock failed".into())));
    let nesting = ProgressNesting::default();
    let mut logger = ProgressLogger::with_clock_and_nesting(bad_clock, nesting.clone());
    assert!(logger.start_progress(0, 1, "bad").is_err());
    assert_eq!(nesting.depth(), 0);
    struct BadBackend;
    impl ProgressBackend for BadBackend {
        fn start_progress(&mut self, _: i64, _: i64, _: &str, _: usize) -> Result<()> {
            Err(Error::Unsupported("start".into()))
        }
        fn set_progress(&mut self, _: i64, _: usize) -> Result<()> {
            Err(Error::Unsupported("set".into()))
        }
        fn next_progress(&mut self) -> Result<i64> {
            Err(Error::Unsupported("next".into()))
        }
        fn end_progress(&mut self, _: usize, _: u64) -> Result<()> {
            Err(Error::Unsupported("end".into()))
        }
    }
    let mut logger = ProgressLogger::with_clock_and_nesting(
        Arc::new(|| Ok(sample(1, 0.0, None))),
        nesting.clone(),
    );
    logger.set_logger(Box::new(BadBackend));
    assert!(logger.start_progress(0, 1, "bad").is_err());
    assert_eq!(nesting.depth(), 0);
    assert!(logger.next_progress().is_err());
    assert!(logger.end_progress(0).is_err());
}

#[test]
fn parallel_next_calls_keep_every_increment_under_same_second_throttling() {
    let (clock, time) = manual_clock(sample(10, 0.0, None));
    let events = Arc::new(Mutex::new(Vec::new()));
    let mut logger = ProgressLogger::with_clock_and_nesting(clock, ProgressNesting::default());
    logger.set_logger(recorder(&events));
    logger.start_progress(0, 100, "threads").unwrap();
    let logger = Arc::new(Mutex::new(logger));
    let workers: Vec<_> = (0..4)
        .map(|_| {
            let logger = logger.clone();
            std::thread::spawn(move || {
                for _ in 0..25 {
                    logger.lock().unwrap().next_progress().unwrap();
                }
            })
        })
        .collect();
    for worker in workers {
        worker.join().unwrap();
    }
    time.lock().unwrap().wall_second = 11;
    logger.lock().unwrap().next_progress().unwrap();
    logger.lock().unwrap().end_progress(0).unwrap();
    let events = events.lock().unwrap();
    assert_eq!(events.len(), 104);
    for (index, event) in events[1..102].iter().enumerate() {
        assert_eq!(event, &Event::Next(index as i64 + 1));
    }
    assert_eq!(events[102], Event::Set(101, 1));
    assert_eq!(events[103], Event::End(0, 0));
}

#[test]
fn command_timer_excludes_header_io_and_active_restart_matches_current_stopwatch() {
    struct SlowWriter {
        output: Vec<u8>,
        time: Arc<Mutex<ProgressTime>>,
    }
    impl Write for SlowWriter {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            self.output.extend_from_slice(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            let mut time = self.time.lock().unwrap();
            time.wall_seconds += 5.0;
            *time.cpu_seconds.as_mut().unwrap() += 2.0;
            Ok(())
        }
    }
    let (clock, time) = manual_clock(sample(1, 0.0, Some(0.0)));
    let mut cmd = CommandProgressLogger::with_clock(
        SlowWriter {
            output: Vec::new(),
            time: time.clone(),
        },
        clock,
    );
    cmd.start_progress(0, 2, "first", 0).unwrap(); // sample wall5/cpu2 after flush
    *time.lock().unwrap() = sample(2, 8.0, Some(3.0));
    cmd.end_progress(0, 0).unwrap();
    let output = String::from_utf8(cmd.into_inner().output).unwrap();
    assert!(output.contains("1.00 s (CPU), 3.00 s (Wall)"), "{output}");

    let (clock, time) = manual_clock(sample(1, 0.0, None));
    let mut cmd = CommandProgressLogger::with_clock(Vec::new(), clock);
    cmd.start_progress(0, 10, "first", 0).unwrap();
    time.lock().unwrap().wall_seconds = 10.0;
    // A second start on a running command backend fails in the Release build
    // too: StopWatch::start throws unconditionally (StopWatch.cpp:43) after the
    // header is printed and the range replaced (fixture row 45).
    assert!(cmd.start_progress(5, 20, "restart", 1).is_err());
    assert_eq!(cmd.next_progress().unwrap(), 6);
    time.lock().unwrap().wall_seconds = 12.0;
    cmd.end_progress(0, 0).unwrap();
    let output = String::from_utf8(cmd.into_inner()).unwrap();
    assert_eq!(
        output,
        "Progress of 'first':\n  Progress of 'restart':\n\r-- done [took unavailable (CPU), 2.00 s (Wall)] -- \n"
    );
}

#[test]
fn system_clock_reports_finite_nondecreasing_process_cpu_without_timing_assumptions() {
    let clock = openms::concept::progress_logger::system_progress_clock();
    let first = clock().unwrap();
    let second = clock().unwrap();
    assert!(first.wall_seconds.is_finite());
    assert!(first.wall_seconds >= 0.0);
    assert!(second.wall_seconds >= first.wall_seconds);
    #[cfg(any(unix, windows))]
    {
        let first_cpu = first.cpu_seconds.expect("supported process CPU clock");
        let second_cpu = second.cpu_seconds.expect("supported process CPU clock");
        assert!(first_cpu.is_finite() && first_cpu >= 0.0);
        assert!(second_cpu.is_finite() && second_cpu >= first_cpu);
    }
    #[cfg(not(any(unix, windows)))]
    {
        assert_eq!(first.cpu_seconds, None);
        assert_eq!(second.cpu_seconds, None);
    }
}

/// Tier-1 differential: every call of the Release oracle driver, replayed in
/// order against the port with one shared nesting context (the source's static
/// depth). The driver waited for a new wall-clock second before each
/// set/next so that every one dispatched; the manual clock advances its second
/// at the same points. Outcome, depth after the call and the bytes the call
/// wrote (timing texts masked) must equal the captured Release row.
#[test]
fn release_range_fixture_replays_call_for_call() {
    let fixture = std::fs::read_to_string(RELEASE_FIXTURE).unwrap();
    let mut lines = fixture.lines().filter(|line| !line.starts_with('#'));
    assert_eq!(
        lines.next(),
        Some("seq\tcase\tlogger\top\ta\tb\tlabel\toutcome\tdepth_after\tdispatch\toutput")
    );
    let (clock, time) = manual_clock(sample(1_000, 0.0, Some(0.0)));
    let nesting = ProgressNesting::default();
    let output = SharedOutput::default();
    let mut loggers: BTreeMap<String, ProgressLogger> = BTreeMap::new();
    let mut rows = 0;
    let mut exceptions = 0;
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 11, "{line}");
        let (seq, name, op) = (fields[0], fields[2], fields[3]);
        let a: i64 = fields[4].parse().unwrap();
        let b: i64 = fields[5].parse().unwrap();
        let label = unescape(fields[6]);
        let (outcome, depth_after, dispatch) = (fields[7], fields[8], fields[9]);
        let expected = unescape(fields[10]);
        let logger = loggers.entry(name.to_owned()).or_insert_with(|| {
            ProgressLogger::with_clock_and_nesting(clock.clone(), nesting.clone())
        });
        {
            // Every sample moves on, so the command timer sees real intervals.
            let mut now = time.lock().unwrap();
            now.wall_seconds += 0.25;
            now.cpu_seconds = now.cpu_seconds.map(|cpu| cpu + 0.01);
            if matches!(op, "set" | "next") {
                assert_eq!(dispatch, "forced", "row {seq}");
                now.wall_second += 1;
            } else {
                assert_eq!(dispatch, "-", "row {seq}");
            }
        }
        let result = match op {
            "type" => {
                let kind = match a {
                    0 => ProgressLogType::Cmd,
                    1 => ProgressLogType::Gui,
                    2 => ProgressLogType::None,
                    other => panic!("row {seq}: log type {other}"),
                };
                logger.set_log_type(kind);
                assert_eq!(logger.log_type(), kind, "row {seq}");
                if kind == ProgressLogType::Cmd {
                    // setLogType(CMD), with stdout redirected to the capture.
                    logger.set_logger(Box::new(CommandProgressLogger::with_clock(
                        output.clone(),
                        clock.clone(),
                    )));
                }
                Ok(())
            }
            "start" => logger.start_progress(a, b, &label),
            "set" => logger.set_progress(a),
            "next" => logger.next_progress(),
            "end" => logger.end_progress(0),
            other => panic!("row {seq}: operation {other}"),
        };
        match outcome {
            "ok" => assert!(result.is_ok(), "row {seq}: {result:?}"),
            thrown => {
                // Release throws Exception::Precondition from StopWatch, which is
                // not a Debug-only check; the port reports the same condition.
                let fields: Vec<&str> = thrown.split('|').collect();
                assert_eq!(
                    fields[..2],
                    ["exception", "Precondition failed"],
                    "row {seq}"
                );
                let message = match fields[3] {
                    "StopWatch.cpp:43" => "progress timer is already running",
                    "StopWatch.cpp:55" => "progress timer is not running",
                    other => panic!("row {seq}: unexpected source location {other}"),
                };
                assert!(
                    matches!(&result, Err(Error::InvalidValue(text)) if text == message),
                    "row {seq}: {result:?}"
                );
                exceptions += 1;
            }
        }
        assert_eq!(nesting.depth().to_string(), depth_after, "row {seq}: depth");
        assert_eq!(mask_timing(&output.take()), expected, "row {seq}: output");
        rows += 1;
    }
    assert_eq!((rows, exceptions), (60, 2));
    assert_eq!(nesting.depth(), 0);
}

/// GUI and NONE: the source's `NoProgressLoggerImpl` and default GUI factory
/// ignore every argument, and the wrapper has no Release range check, so an
/// inverted range is accepted, reaches a custom GUI backend unchanged and
/// still nests (fixture rows 51-60 cover the default backends).
#[test]
fn gui_and_none_accept_an_inverted_range() {
    let events = Arc::new(Mutex::new(Vec::new()));
    let (clock, time) = manual_clock(sample(5, 0.0, None));
    let nesting = ProgressNesting::default();
    let mut none = ProgressLogger::with_clock_and_nesting(clock.clone(), nesting.clone());
    none.start_progress(5, 0, "none inverted").unwrap();
    assert_eq!(nesting.depth(), 1);
    let mut gui = ProgressLogger::with_clock_and_nesting(clock, nesting.clone());
    let factory_events = events.clone();
    gui.set_gui_factory(Arc::new(move || recorder(&factory_events)));
    gui.set_log_type(ProgressLogType::Gui);
    gui.start_progress(i64::MAX, i64::MIN, "gui inverted")
        .unwrap();
    assert_eq!(nesting.depth(), 2);
    time.lock().unwrap().wall_second = 6;
    gui.set_progress(3).unwrap();
    gui.end_progress(0).unwrap();
    none.set_progress(3).unwrap();
    none.end_progress(0).unwrap();
    assert_eq!(nesting.depth(), 0);
    assert_eq!(
        *events.lock().unwrap(),
        vec![
            Event::Start(i64::MAX, i64::MIN, "gui inverted".into(), 1),
            Event::Set(3, 2),
            Event::End(1, 0),
        ]
    );
}

// Sections a finished call left open. The source keeps them in its static
// depth for the rest of the process (`ProgressLogger.h:105`, `.cpp:266-269`);
// the expectations below follow from the documented contract of
// `ProgressNesting` and `ProgressReporter`, which bounds that to the logger
// that left them, not from Rust output.

/// A logger whose calls go to a recording backend.
fn recording(
    clock: &ProgressClock,
    nesting: &ProgressNesting,
    events: &Arc<Mutex<Vec<Event>>>,
) -> ProgressLogger {
    let mut logger = ProgressLogger::with_clock_and_nesting(clock.clone(), nesting.clone());
    logger.set_logger(recorder(events));
    logger
}

#[test]
fn an_abandoned_section_indents_only_its_own_logger_until_it_is_dropped() {
    let (clock, _) = manual_clock(sample(1, 0.0, None));
    let nesting = ProgressNesting::default();
    let events = Arc::new(Mutex::new(Vec::new()));
    let mut failed = recording(&clock, &nesting, &events);
    let mut other = recording(&clock, &nesting, &events);
    // A call that returns inside its section, as a reader that fails does.
    ProgressReporter::new(Some(&mut failed))
        .start(0, 2, "left open")
        .unwrap();
    assert_eq!(nesting.depth(), 1);
    other.start_progress(0, 2, "other").unwrap();
    other.end_progress(0).unwrap();
    failed.start_progress(0, 2, "again").unwrap();
    failed.end_progress(0).unwrap();
    assert_eq!(nesting.depth(), 1);
    let length = events.lock().unwrap().len();
    drop(failed);
    assert_eq!(nesting.depth(), 0);
    // Dropping the logger ended nothing.
    assert_eq!(events.lock().unwrap().len(), length);
    assert_eq!(
        *events.lock().unwrap(),
        vec![
            Event::Start(0, 2, "left open".into(), 0),
            Event::Start(0, 2, "other".into(), 0),
            Event::End(0, 0),
            Event::Start(0, 2, "again".into(), 1),
            Event::End(1, 0),
        ]
    );
}

#[test]
fn abandoned_sections_never_reach_the_nesting_bound() {
    let (clock, _) = manual_clock(sample(1, 0.0, None));
    let nesting = ProgressNesting::default();
    let events = Arc::new(Mutex::new(Vec::new()));
    let mut logger = recording(&clock, &nesting, &events);
    for _ in 0..=MAX_PROGRESS_DEPTH {
        ProgressReporter::new(Some(&mut logger))
            .start(0, 0, "left open")
            .unwrap();
    }
    assert_eq!(nesting.depth(), MAX_PROGRESS_DEPTH + 1);
    events.lock().unwrap().clear();
    // No section is open, so one may start; its backend is never handed a
    // depth beyond the bound, which only limits indentation.
    logger.start_progress(0, 0, "after").unwrap();
    logger.end_progress(0).unwrap();
    assert_eq!(
        *events.lock().unwrap(),
        vec![
            Event::Start(0, 0, "after".into(), MAX_PROGRESS_DEPTH),
            Event::End(MAX_PROGRESS_DEPTH, 0),
        ]
    );
    // The bound still limits the sections open at once.
    for _ in 0..MAX_PROGRESS_DEPTH {
        logger.start_progress(0, 0, "open").unwrap();
    }
    assert!(logger.start_progress(0, 0, "one too many").is_err());
    drop(logger);
    assert_eq!(nesting.depth(), 0);
}

/// The mzML reader reports its document section through short-lived
/// reporters on a copy of the caller's logger while the reporter of its lists
/// lives for the whole load. While any reporter on a logger or its copies
/// lives, the call is not over and nothing is abandoned: the document section
/// stays open for every logger.
#[test]
fn a_section_stays_open_while_any_reporter_on_its_logger_lives() {
    let (clock, _) = manual_clock(sample(1, 0.0, None));
    let nesting = ProgressNesting::default();
    let events = Arc::new(Mutex::new(Vec::new()));
    let mut logger = recording(&clock, &nesting, &events);
    let mut other = recording(&clock, &nesting, &events);
    let mut copy = logger.clone();
    let mut lists = ProgressReporter::new(Some(&mut logger));
    ProgressReporter::new(Some(&mut copy))
        .start(0, 1, "document")
        .unwrap();
    other.start_progress(0, 0, "another logger").unwrap();
    other.end_progress(0).unwrap();
    lists.start(0, 2, "list").unwrap();
    lists.end().unwrap();
    ProgressReporter::new(Some(&mut copy)).end().unwrap();
    drop(lists);
    assert_eq!(nesting.depth(), 0);
    assert_eq!(
        *events.lock().unwrap(),
        vec![
            Event::Start(0, 0, "another logger".into(), 1),
            Event::End(1, 0),
            Event::Start(0, 2, "list".into(), 1),
            Event::End(1, 0),
        ]
    );
}

/// A section started directly on a logger, with no reporter, is not
/// abandoned when a call ends, but it leaves the depth when the last of the
/// logger and its copies is dropped, without an end.
#[test]
fn dropping_the_last_copy_of_a_logger_takes_its_open_sections_out_of_the_depth() {
    let (clock, _) = manual_clock(sample(1, 0.0, None));
    let nesting = ProgressNesting::default();
    let events = Arc::new(Mutex::new(Vec::new()));
    let mut logger = recording(&clock, &nesting, &events);
    logger.start_progress(0, 1, "direct").unwrap();
    let copy = logger.clone();
    drop(logger);
    assert_eq!(nesting.depth(), 1);
    drop(copy);
    assert_eq!(nesting.depth(), 0);
    assert_eq!(
        *events.lock().unwrap(),
        vec![Event::Start(0, 1, "direct".into(), 0)]
    );
}
