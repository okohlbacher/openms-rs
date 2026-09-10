// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::concept::log_stream::*;
use std::io::{self, Write};
use std::sync::{
    Arc, Mutex,
    atomic::{AtomicUsize, Ordering},
};

#[derive(Clone, Default)]
struct Capture(Arc<Mutex<Vec<u8>>>);
impl Capture {
    fn bytes(&self) -> Vec<u8> {
        self.0.lock().unwrap().clone()
    }
    fn text(&self) -> String {
        String::from_utf8(self.bytes()).unwrap()
    }
}
impl Write for Capture {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn captured() -> (LogStream, LogSink, Capture) {
    let capture = Capture::default();
    let sink = LogSink::new(capture.clone());
    let mut logger = LogStream::new("UNKNOWN_LOG_LEVEL").unwrap();
    logger.insert(&sink).unwrap();
    (logger, sink, capture)
}
fn line(log: &mut LogStream, text: &str) {
    writeln!(log, "{text}").unwrap();
    log.flush().unwrap();
}

#[test]
fn literal_source_general_and_color_files() {
    for (color, expected) in [
        (
            None,
            include_bytes!("data/log_stream_general.txt").as_slice(),
        ),
        (
            Some(LogColor::Red),
            include_bytes!("data/log_stream_general_red.txt").as_slice(),
        ),
        (
            Some(LogColor::Yellow),
            include_bytes!("data/log_stream_general_yellow.txt").as_slice(),
        ),
    ] {
        let (mut log, _, capture) = captured();
        log.set_color(color);
        line(&mut log, "1");
        line(&mut log, "2");
        log.finish().unwrap();
        assert_eq!(capture.bytes(), expected);
    }
}

#[test]
fn source_cache_fixture_and_empty_lines() {
    let (mut log, _, capture) = captured();
    for _ in 0..3 {
        line(&mut log, "This is a repeptitive message");
        line(&mut log, "This is another repeptitive message");
    }
    line(&mut log, "This is a non-repetitive message");
    log.finish().unwrap();
    assert_eq!(
        capture.bytes(),
        include_bytes!("data/log_stream_caching.txt")
    );
    let (mut log, _, capture) = captured();
    line(&mut log, "No caching for the following empty lines");
    log.write_all(b"\n\n\n\n").unwrap();
    log.finish().unwrap();
    assert_eq!(
        capture.text(),
        "No caching for the following empty lines\n\n\n\n\n"
    );
}

#[test]
fn source_flush_waits_for_sync_and_preserves_cr_and_raw_bytes() {
    let (mut log, _, capture) = captured();
    line(&mut log, "flushtest");
    log.write_all(b"unfinishedline...\n").unwrap();
    assert_eq!(capture.text(), "flushtest\n");
    log.flush().unwrap();
    assert_eq!(capture.text(), "flushtest\nunfinishedline...\n");
    log.write_all(&[0xff, b'\r', b'\n', b'p']).unwrap();
    log.flush().unwrap();
    assert!(capture.bytes().ends_with(&[0xff, b'\r', b'\n']));
    log.flush_incomplete().unwrap();
    assert!(capture.bytes().ends_with(b"p\n"));
}

#[test]
fn exact_prefix_tokens_and_per_route_clock_queries() {
    let (mut log, first, capture) = captured();
    let second_capture = Capture::default();
    let second = LogSink::new(second_capture.clone());
    log.insert(&second).unwrap();
    log.set_level("DEVELOPMENT").unwrap();
    let calls = Arc::new(AtomicUsize::new(0));
    let calls_copy = calls.clone();
    log.set_clock(Some(Arc::new(move || {
        calls_copy.fetch_add(1, Ordering::Relaxed);
        Ok(LogTime {
            year: 2007,
            month: 9,
            day: 8,
            hour: 1,
            minute: 2,
            second: 3,
        })
    })));
    log.set_prefix(&first, "%y|%T|%t|%D|%d|%S|%s|%%|%q|%")
        .unwrap();
    log.set_prefix(&second, "%T ").unwrap();
    line(&mut log, "payload");
    assert_eq!(
        capture.text(),
        "DEVELOPMENT|01:02:03|01:02|2007/09/08|09/08|2007/09/08, 01:02:03|09/08, 01:02|%||payload\n"
    );
    assert_eq!(second_capture.text(), "01:02:03 payload\n");
    assert_eq!(calls.load(Ordering::Relaxed), 2);
}

#[test]
fn prefixes_affect_existing_routes_and_cache_summaries_use_current_level() {
    let (mut log, first, capture) = captured();
    log.set_all_prefixes("%y: ").unwrap();
    log.insert(&first).unwrap(); // identity dedup does not discard original prefix
    line(&mut log, "a");
    line(&mut log, "a");
    log.set_level("CHANGED").unwrap();
    let other = Capture::default();
    log.insert(&LogSink::new(other.clone())).unwrap();
    log.clear_cache().unwrap();
    assert_eq!(
        capture.text(),
        "UNKNOWN_LOG_LEVEL: a\nCHANGED: <a> occurred 2 times\n"
    );
    assert_eq!(other.text(), "<a> occurred 2 times\n");
}

#[test]
fn lru_refresh_lexical_cache_clear_and_partial_bypass() {
    let (mut log, _, capture) = captured();
    for text in ["b", "a", "b", "c", "c", "b"] {
        line(&mut log, text);
    }
    log.write_all(b"b").unwrap();
    log.flush_incomplete().unwrap();
    log.finish().unwrap();
    assert_eq!(
        capture.text(),
        "b\na\nc\nb\n<b> occurred 3 times\n<c> occurred 2 times\n"
    );
}

#[test]
fn removal_preserves_partial_and_cache_but_remove_all_flushes_partial() {
    let (mut log, sink, capture) = captured();
    log.write_all(b"old").unwrap();
    log.remove(&sink).unwrap();
    log.write_all(b"discarded\n").unwrap();
    log.flush().unwrap();
    log.insert(&sink).unwrap();
    line(&mut log, "new");
    line(&mut log, "oldnew");
    log.write_all(b"tail").unwrap();
    log.remove_all_streams().unwrap();
    assert_eq!(capture.text(), "oldnew\ntail\n");
    log.insert(&sink).unwrap();
    log.finish().unwrap();
    assert_eq!(capture.text(), "oldnew\ntail\n<oldnew> occurred 2 times\n");
}

#[test]
fn guards_flush_boundaries_nest_restore_and_drop_route_metadata() {
    let (mut log, sink, capture) = captured();
    let other = Capture::default();
    log.insert(&LogSink::new(other.clone())).unwrap();
    log.set_prefix(&sink, "P:").unwrap();
    log.write_all(b"pending_before_guard").unwrap();
    {
        let mut outer = LogSinkGuard::new(&mut log, &sink).unwrap();
        {
            let mut inner = LogSinkGuard::new(&mut outer, &sink).unwrap();
            inner.write_all(b"unterminated_while_guarded").unwrap();
        }
        outer.restore().unwrap();
    }
    line(&mut log, "after_guard");
    assert_eq!(capture.text(), "P:pending_before_guard\nafter_guard\n");
    assert_eq!(
        other.text(),
        "pending_before_guard\nunterminated_while_guarded\nafter_guard\n"
    );
    let caught = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _guard = LogSinkGuard::new(&mut log, &sink).unwrap();
        panic!("source exception scenario");
    }));
    assert!(caught.is_err());
    assert!(log.has_stream(&sink));
}

#[test]
fn notifications_receive_rendered_lines_after_sink_lock_release() {
    let (mut log, sink, capture) = captured();
    let (nested, _, nested_capture) = captured();
    let nested = Arc::new(Mutex::new(nested));
    let notifications = Capture::default();
    let notified = notifications.clone();
    log.insert_notification(
        &sink,
        Arc::new(move |bytes| {
            let mut notified = notified.clone();
            notified.write_all(bytes)?;
            line(&mut nested.lock().unwrap(), "callback");
            Ok(())
        }),
    )
    .unwrap();
    log.set_prefix(&sink, "[x]").unwrap();
    line(&mut log, "line");
    line(&mut log, "line");
    assert_eq!(notifications.text(), "[x]line\n");
    assert_eq!(capture.text(), "[x]line\n");
    assert_eq!(nested_capture.text(), "callback\n");
    log.remove(&sink).unwrap();
    line(&mut log, "ignored");
}

#[test]
fn configuration_copy_has_fresh_buffers_cache_and_shared_owned_sinks() {
    let (mut original, _, capture) = captured();
    line(&mut original, "a");
    original.write_all(b"pending").unwrap();
    let mut copied = original.clone_configuration();
    line(&mut copied, "a");
    copied.finish().unwrap();
    assert_eq!(capture.text(), "a\na\n");
    original.finish().unwrap();
    assert_eq!(capture.text(), "a\na\npending\n");
    let mut unbound = LogStream::default();
    unbound.insert(&LogSink::new(capture.clone())).unwrap();
    unbound.set_level("changed").unwrap();
    line(&mut unbound, "not printed");
    assert_eq!(unbound.level(), UNKNOWN_LOG_LEVEL);
}

#[test]
fn exact_32k_buffer_boundary_and_chunked_lines() {
    let (mut log, _, capture) = captured();
    let mut bytes = vec![b'x'; LOG_BUFFER_BYTES - 1];
    bytes.push(b'\n');
    log.write_all(&bytes[..LOG_BUFFER_BYTES - 1]).unwrap();
    assert!(capture.bytes().is_empty());
    log.write_all(b"\n").unwrap();
    assert_eq!(capture.bytes(), bytes);
    let large = vec![b'y'; LOG_BUFFER_BYTES * 3];
    log.write_all(&large).unwrap();
    assert_eq!(capture.bytes(), bytes);
    log.flush_incomplete().unwrap();
    assert_eq!(capture.bytes().len(), bytes.len() + large.len() + 1);
}

#[test]
fn limits_fail_before_writing_and_errors_do_not_replay_on_drop() {
    let (mut log, sink, capture) = captured();
    assert!(
        log.set_prefix(&sink, "x".repeat(MAX_LOG_TEXT_BYTES + 1))
            .is_err()
    );
    assert!(
        log.write_all(&vec![b'x'; MAX_LOG_OPERATION_BYTES + 1])
            .is_err()
    );
    assert!(capture.bytes().is_empty());
    // A small stored prefix can expand dramatically through a long level.
    log.set_level("x".repeat(MAX_LOG_TEXT_BYTES)).unwrap();
    log.set_prefix(&sink, "%y".repeat(1024)).unwrap();
    log.write_all(b"line\n").unwrap();
    assert!(log.flush().is_err());
    assert!(capture.bytes().is_empty());
    assert_eq!(log.flush().unwrap_err().kind(), io::ErrorKind::BrokenPipe);
    struct Failing(Capture);
    impl Write for Failing {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            if self.0.bytes().is_empty() {
                self.0.write(&bytes[..1])
            } else {
                Err(io::Error::other("injected write error"))
            }
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }
    let capture = Capture::default();
    {
        let mut log = LogStream::new("x").unwrap();
        log.insert(&LogSink::new(Failing(capture.clone()))).unwrap();
        log.write_all(b"payload\n").unwrap();
        assert!(log.flush().is_err());
    }
    assert_eq!(capture.text(), "p");
}

#[test]
fn local_clock_default_and_checked_injected_calendar() {
    let (mut log, sink, capture) = captured();
    log.set_prefix(&sink, "%D %T ").unwrap();
    line(&mut log, "now");
    let text = capture.text();
    assert_eq!(text.len(), 24);
    assert_eq!(&text[4..5], "/");
    assert_eq!(&text[10..11], " ");
    assert_eq!(&text[13..14], ":");
    let (mut log, sink, capture) = captured();
    log.set_prefix(&sink, "%D").unwrap();
    log.set_clock(Some(Arc::new(|| {
        Ok(LogTime {
            year: 2025,
            month: 2,
            day: 29,
            hour: 0,
            minute: 0,
            second: 0,
        })
    })));
    log.write_all(b"bad\n").unwrap();
    assert!(log.flush().is_err());
    assert!(capture.bytes().is_empty());
}

#[test]
fn global_tls_snapshots_thread_safety_and_location_helpers() {
    // This is the only test touching global configuration, so other tests stay isolated.
    let first = Capture::default();
    let first_sink = LogSink::new(first.clone());
    with_global_log(LogLevel::Debug, |log| {
        log.remove_all_streams()?;
        log.clear_cache()?;
        log.set_color(None);
        log.insert(&first_sink)
    })
    .unwrap();
    std::thread::scope(|scope| {
        let barrier = Arc::new(std::sync::Barrier::new(2));
        let child_barrier = barrier.clone();
        let child = scope.spawn(move || {
            log_message(LogLevel::Debug, format_args!("first snapshot")).unwrap();
            child_barrier.wait();
            child_barrier.wait();
            log_message_at(
                LogLevel::Debug,
                "/source/tree/tool.rs",
                42,
                format_args!("still first"),
            )
            .unwrap();
        });
        barrier.wait();
        let second = Capture::default();
        with_global_log(LogLevel::Debug, |log| {
            log.remove_all_streams()?;
            log.insert(&LogSink::new(second.clone()))
        })
        .unwrap();
        barrier.wait();
        child.join().unwrap();
        assert_eq!(first.text(), "first snapshot\ntool.rs(42): still first\n");
        let mut threads = Vec::new();
        for thread in 0..8 {
            threads.push(scope.spawn(move || {
                for index in 0..625 {
                    openms::openms_log_debug_nofile!("racing_line_{}", thread * 625 + index)
                        .unwrap();
                }
            }));
        }
        for thread in threads {
            thread.join().unwrap();
        }
        let mut values: Vec<usize> = second
            .text()
            .lines()
            .map(|line| line.strip_prefix("racing_line_").unwrap().parse().unwrap())
            .collect();
        values.sort_unstable();
        assert_eq!(values, (0..5000).collect::<Vec<_>>());
    });
    with_global_log(LogLevel::Debug, |log| log.remove_all_streams()).unwrap();
}

#[test]
fn all_source_single_styles_restore_their_matching_ansi_state() {
    for (style, on, off) in [
        (LogColor::Red, "91", "39"),
        (LogColor::Green, "92", "39"),
        (LogColor::Yellow, "93", "39"),
        (LogColor::Blue, "94", "39"),
        (LogColor::Magenta, "95", "39"),
        (LogColor::Cyan, "96", "39"),
        (LogColor::Underline, "4", "24"),
        (LogColor::Bright, "1", "22"),
        (LogColor::Invert, "7", "27"),
    ] {
        let (mut log, _, capture) = captured();
        log.set_color(Some(style));
        line(&mut log, "value");
        assert_eq!(capture.text(), format!("\x1b[{on}mvalue\x1b[{off}m\n"));
    }
}

#[test]
fn line_limits_and_shared_formatting_budget() {
    let (mut log, _, capture) = captured();
    let chunk = vec![b'x'; LOG_BUFFER_BYTES];
    for _ in 0..MAX_LOG_LINE_BYTES / LOG_BUFFER_BYTES {
        log.write_all(&chunk).unwrap();
    }
    log.write_all(b"overflow").unwrap();
    assert!(log.flush().is_err());
    assert!(capture.bytes().is_empty());
    // Individual fmt calls are small but together must consume one budget.
    struct Repeated;
    impl std::fmt::Display for Repeated {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            for _ in 0..MAX_LOG_OPERATION_BYTES / 4096 + 1 {
                f.write_str(&"x".repeat(4096))?;
            }
            Ok(())
        }
    }
    let mut no_sinks = LogStream::new("ignored").unwrap();
    assert!(write!(&mut no_sinks, "{Repeated}").is_err());
}

#[test]
fn removal_flush_failure_disables_output_and_destructor_retry() {
    struct FailFirstFlush {
        capture: Capture,
        flushes: Arc<AtomicUsize>,
    }
    impl Write for FailFirstFlush {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            self.capture.write(bytes)
        }
        fn flush(&mut self) -> io::Result<()> {
            if self.flushes.fetch_add(1, Ordering::Relaxed) == 0 {
                Err(io::Error::other("one-shot flush failure"))
            } else {
                Ok(())
            }
        }
    }
    let capture = Capture::default();
    let flushes = Arc::new(AtomicUsize::new(0));
    let sink = LogSink::new(FailFirstFlush {
        capture: capture.clone(),
        flushes: flushes.clone(),
    });
    {
        let mut log = LogStream::new("test").unwrap();
        log.insert(&sink).unwrap();
        assert_eq!(
            log.remove_all_streams().unwrap_err().to_string(),
            "one-shot flush failure"
        );
        assert!(log.has_stream(&sink));
        assert_eq!(
            log.write_all(b"must not be delivered\n")
                .unwrap_err()
                .kind(),
            io::ErrorKind::BrokenPipe
        );
        assert_eq!(log.flush().unwrap_err().kind(), io::ErrorKind::BrokenPipe);
        assert_eq!(log.finish().unwrap_err().kind(), io::ErrorKind::BrokenPipe);
        assert_eq!(
            log.remove_all_streams().unwrap_err().kind(),
            io::ErrorKind::BrokenPipe
        );
    }
    assert!(capture.bytes().is_empty());
    assert_eq!(flushes.load(Ordering::Relaxed), 1);
}
