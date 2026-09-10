// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::concept::progress_logger::{
    CommandProgressLogger, MAX_PROGRESS_DEPTH, MAX_PROGRESS_LABEL_BYTES, ProgressBackend,
    ProgressClock, ProgressLogType, ProgressLogger, ProgressNesting, ProgressTime,
};
use openms::{Error, Result};
use std::collections::VecDeque;
use std::io::{self, Write};
use std::sync::{Arc, Mutex};

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
    assert!(logger.start_progress(2, 1, "bad").is_err());
    assert!(
        logger
            .start_progress(0, 1, &"x".repeat(MAX_PROGRESS_LABEL_BYTES + 1))
            .is_err()
    );
    assert_eq!(nesting.depth(), 0);
    assert!(events.lock().unwrap().is_empty());
    let mut cmd = CommandProgressLogger::with_clock(Vec::new(), clock);
    assert!(cmd.end_progress(0, 0).is_err());
    cmd.start_progress(i64::MAX, i64::MAX, "counter", 0)
        .unwrap();
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
