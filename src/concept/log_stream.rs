// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Owned, bounded logging streams with OpenMS line, route and duplicate-cache semantics.
//!
//! Writers are shared explicitly through [`LogSink`]. Callbacks run after the shared
//! output lock is released. Errors propagate; output already written to a sink cannot
//! be rolled back. After an output failure the stream rejects further output so its
//! destructor cannot replay a partially delivered record. Use [`LogStream::finish`]
//! to observe final flush errors; destruction can only make a best-effort flush.

use chrono::{Datelike, Timelike};
use std::cell::RefCell;
use std::fmt;
use std::io::{self, IsTerminal, Write};
use std::ops::{Deref, DerefMut};
use std::sync::{Arc, Mutex, OnceLock, RwLock};

pub const UNKNOWN_LOG_LEVEL: &str = "UNKNOWN_LOG_LEVEL";
pub const LOG_BUFFER_BYTES: usize = 32_768;
pub const MAX_LOG_LINE_BYTES: usize = 1_048_576;
pub const MAX_LOG_TEXT_BYTES: usize = 65_536;
pub const MAX_LOG_SINKS: usize = 1024;
pub const MAX_LOG_OPERATION_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_LOG_OPERATION_WORK: usize = 50_000_000;

fn invalid(message: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message)
}
fn poisoned() -> io::Error {
    io::Error::other("logging lock poisoned")
}

#[derive(Default)]
struct Budget {
    work: usize,
    bytes: usize,
}
impl Budget {
    fn spend(&mut self, work: usize, bytes: usize) -> io::Result<()> {
        self.work = self
            .work
            .checked_add(work)
            .ok_or_else(|| invalid("logging work overflow"))?;
        self.bytes = self
            .bytes
            .checked_add(bytes)
            .ok_or_else(|| invalid("logging allocation overflow"))?;
        if self.work > MAX_LOG_OPERATION_WORK || self.bytes > MAX_LOG_OPERATION_BYTES {
            return Err(invalid("logging operation resource limit exceeded"));
        }
        Ok(())
    }
}

/// A stable, cloneable output identity. Cloning shares the writer, not its contents.
#[derive(Clone)]
pub struct LogSink {
    writer: Arc<Mutex<Box<dyn Write + Send>>>,
    kind: SinkKind,
}
#[derive(Clone, Copy, Debug)]
enum SinkKind {
    Other,
    Stdout,
    Stderr,
}
impl fmt::Debug for LogSink {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LogSink")
            .field("kind", &self.kind)
            .finish_non_exhaustive()
    }
}
impl PartialEq for LogSink {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.writer, &other.writer)
    }
}
impl Eq for LogSink {}
impl LogSink {
    fn color_allowed(&self) -> bool {
        match self.kind {
            SinkKind::Other => true,
            SinkKind::Stdout => io::stdout().is_terminal(),
            SinkKind::Stderr => io::stderr().is_terminal(),
        }
    }
    /// Arbitrary writers receive ANSI colors when a colored logger uses them, as in C++.
    pub fn new(writer: impl Write + Send + 'static) -> Self {
        Self {
            writer: Arc::new(Mutex::new(Box::new(writer))),
            kind: SinkKind::Other,
        }
    }
    /// The process stdout identity. ANSI coloring is enabled only for a terminal.
    pub fn stdout() -> Self {
        static SINK: OnceLock<LogSink> = OnceLock::new();
        SINK.get_or_init(|| Self {
            writer: Arc::new(Mutex::new(Box::new(io::stdout()))),
            kind: SinkKind::Stdout,
        })
        .clone()
    }
    /// The process stderr identity. ANSI coloring is enabled only for a terminal.
    pub fn stderr() -> Self {
        static SINK: OnceLock<LogSink> = OnceLock::new();
        SINK.get_or_init(|| Self {
            writer: Arc::new(Mutex::new(Box::new(io::stderr()))),
            kind: SinkKind::Stderr,
        })
        .clone()
    }
}

/// A notification receives the complete rendered record, including its trailing LF.
/// Use a captured shared writer when an accumulating source-style notifier buffer is needed.
pub type LogNotification = Arc<dyn Fn(&[u8]) -> io::Result<()> + Send + Sync>;

/// Calendar fields used by the source's local-time prefix tokens.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct LogTime {
    pub year: u16,
    pub month: u8,
    pub day: u8,
    pub hour: u8,
    pub minute: u8,
    pub second: u8,
}
impl LogTime {
    fn validate(self) -> io::Result<Self> {
        let leap = self.year % 4 == 0 && (self.year % 100 != 0 || self.year % 400 == 0);
        let days = match self.month {
            2 => {
                if leap {
                    29
                } else {
                    28
                }
            }
            4 | 6 | 9 | 11 => 30,
            1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
            _ => 0,
        };
        if self.year > 9999
            || self.day == 0
            || self.day > days
            || self.hour > 23
            || self.minute > 59
            || self.second > 60
        {
            return Err(invalid("invalid logging calendar fields"));
        }
        Ok(self)
    }
}
/// Supplies local calendar fields. No implicit UTC substitution is made for local time.
pub type LogClock = Arc<dyn Fn() -> io::Result<LogTime> + Send + Sync>;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LogColor {
    Red,
    Green,
    Yellow,
    Blue,
    Magenta,
    Cyan,
    Underline,
    Bright,
    Invert,
}
impl LogColor {
    fn enable(self) -> &'static [u8] {
        match self {
            Self::Red => b"\x1b[91m",
            Self::Green => b"\x1b[92m",
            Self::Yellow => b"\x1b[93m",
            Self::Blue => b"\x1b[94m",
            Self::Magenta => b"\x1b[95m",
            Self::Cyan => b"\x1b[96m",
            Self::Underline => b"\x1b[4m",
            Self::Bright => b"\x1b[1m",
            Self::Invert => b"\x1b[7m",
        }
    }
    fn disable(self) -> &'static [u8] {
        match self {
            Self::Underline => b"\x1b[24m",
            Self::Bright => b"\x1b[22m",
            Self::Invert => b"\x1b[27m",
            _ => b"\x1b[39m",
        }
    }
}

#[derive(Clone)]
struct Route {
    sink: LogSink,
    prefix: String,
    notification: Option<LogNotification>,
}
struct CachedLine {
    text: Vec<u8>,
    occurrences: u64,
}

/// Source-compatible buffered line logger. `Default` is the source's unbound stream.
/// Use `new` to create a bound stream, initially without destinations.
pub struct LogStream {
    bound: bool,
    level: String,
    routes: Vec<Route>,
    pending: Vec<u8>,
    incomplete: Vec<u8>,
    // At most two entries, oldest last-use first. No general cache machinery is needed.
    cache: Vec<CachedLine>,
    color: Option<LogColor>,
    clock: Option<LogClock>,
    failed: bool,
}
impl fmt::Debug for LogStream {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LogStream")
            .field("bound", &self.bound)
            .field("level", &self.level)
            .field("sinks", &self.routes.len())
            .field("failed", &self.failed)
            .finish_non_exhaustive()
    }
}
impl Default for LogStream {
    fn default() -> Self {
        Self {
            bound: false,
            level: String::new(),
            routes: Vec::new(),
            pending: Vec::new(),
            incomplete: Vec::new(),
            cache: Vec::new(),
            color: None,
            clock: Some(Arc::new(local_time)),
            failed: false,
        }
    }
}
impl LogStream {
    pub fn new(level: impl AsRef<str>) -> io::Result<Self> {
        let level = level.as_ref();
        check_text(level)?;
        let mut logger = Self::default();
        logger.bound = true;
        logger.level = level.to_owned();
        Ok(logger)
    }
    pub fn level(&self) -> &str {
        if self.bound {
            &self.level
        } else {
            UNKNOWN_LOG_LEVEL
        }
    }
    pub fn set_level(&mut self, level: impl AsRef<str>) -> io::Result<()> {
        if self.bound {
            check_text(level.as_ref())?;
            self.level = level.as_ref().to_owned();
        }
        Ok(())
    }
    pub fn set_color(&mut self, color: Option<LogColor>) {
        self.color = color;
    }
    pub fn set_clock(&mut self, clock: Option<LogClock>) {
        self.clock = clock;
    }
    /// Copy only configuration, with fresh line buffers and cache (the source TLS rule).
    pub fn clone_configuration(&self) -> Self {
        Self {
            bound: self.bound,
            level: self.level.clone(),
            routes: self.routes.clone(),
            color: self.color,
            clock: self.clock.clone(),
            pending: Vec::new(),
            incomplete: Vec::new(),
            cache: Vec::new(),
            failed: false,
        }
    }
    pub fn has_stream(&self, sink: &LogSink) -> bool {
        self.routes.iter().any(|route| route.sink == *sink)
    }
    pub fn insert(&mut self, sink: &LogSink) -> io::Result<()> {
        if !self.bound || self.has_stream(sink) {
            return Ok(());
        }
        if self.routes.len() == MAX_LOG_SINKS {
            return Err(invalid("too many logging sinks"));
        }
        self.routes.push(Route {
            sink: sink.clone(),
            prefix: String::new(),
            notification: None,
        });
        Ok(())
    }
    pub fn insert_notification(
        &mut self,
        sink: &LogSink,
        callback: LogNotification,
    ) -> io::Result<()> {
        self.insert(sink)?;
        if let Some(route) = self.routes.iter_mut().find(|r| r.sink == *sink) {
            route.notification = Some(callback);
        }
        Ok(())
    }
    pub fn remove(&mut self, sink: &LogSink) -> io::Result<()> {
        if let Some(index) = self.routes.iter().position(|r| r.sink == *sink) {
            self.flush()?; // Source keeps incomplete text and duplicate cache across removal.
            self.routes.remove(index);
        }
        Ok(())
    }
    pub fn remove_all_streams(&mut self) -> io::Result<()> {
        let result = (|| {
            self.flush_incomplete()?;
            for route in &self.routes {
                route.sink.writer.lock().map_err(|_| poisoned())?.flush()?;
            }
            self.routes.clear();
            Ok(())
        })();
        self.finish_operation(result)
    }
    pub fn set_prefix(&mut self, sink: &LogSink, prefix: impl AsRef<str>) -> io::Result<()> {
        let Some(index) = self.routes.iter().position(|r| r.sink == *sink) else {
            return Ok(());
        };
        let prefix = prefix.as_ref();
        check_text(prefix)?;
        let bytes = self.routes.iter().map(|r| r.prefix.len()).sum::<usize>()
            - self.routes[index].prefix.len()
            + prefix.len();
        if bytes > MAX_LOG_OPERATION_BYTES {
            return Err(invalid("logging prefix storage limit exceeded"));
        }
        self.routes[index].prefix = prefix.to_owned();
        Ok(())
    }
    /// Affects existing routes only; subsequently inserted destinations start without a prefix.
    pub fn set_all_prefixes(&mut self, prefix: impl AsRef<str>) -> io::Result<()> {
        if !self.bound {
            return Ok(());
        }
        let prefix = prefix.as_ref();
        check_text(prefix)?;
        if prefix.len().saturating_mul(self.routes.len()) > MAX_LOG_OPERATION_BYTES {
            return Err(invalid("logging prefix storage limit exceeded"));
        }
        for route in &mut self.routes {
            route.prefix = prefix.to_owned();
        }
        Ok(())
    }
    pub fn clear_cache(&mut self) -> io::Result<()> {
        self.ready()?;
        let result = self.clear_cache_inner(&mut Budget::default());
        self.finish_operation(result)
    }
    pub fn flush_incomplete(&mut self) -> io::Result<()> {
        self.ready()?;
        let mut budget = Budget::default();
        let result = self
            .sync(&mut budget)
            .and_then(|()| self.partial(&mut budget));
        self.finish_operation(result)
    }
    /// Flush complete lines, then cache summaries, then partial text, as in destruction.
    pub fn finish(&mut self) -> io::Result<()> {
        self.ready()?;
        let mut budget = Budget::default();
        let result = self
            .sync(&mut budget)
            .and_then(|()| self.clear_cache_inner(&mut budget))
            .and_then(|()| self.partial(&mut budget));
        self.finish_operation(result)
    }
    fn ready(&self) -> io::Result<()> {
        if self.failed {
            Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "logging stream previously failed",
            ))
        } else {
            Ok(())
        }
    }
    fn finish_operation<T>(&mut self, result: io::Result<T>) -> io::Result<T> {
        if result.is_err() {
            self.failed = true;
        }
        result
    }
    fn partial(&mut self, budget: &mut Budget) -> io::Result<()> {
        if !self.incomplete.is_empty() {
            let line = std::mem::take(&mut self.incomplete);
            self.distribute(&line, budget)?;
        }
        Ok(())
    }
    fn clear_cache_inner(&mut self, budget: &mut Budget) -> io::Result<()> {
        // The C++ map traverses keys lexically, independently of LRU eviction order.
        budget.spend(self.cache.iter().map(|line| line.text.len()).sum(), 0)?;
        self.cache.sort_by(|a, b| a.text.cmp(&b.text));
        let cache = std::mem::take(&mut self.cache);
        for line in cache {
            self.summary(&line, budget)?;
        }
        Ok(())
    }
    fn summary(&self, line: &CachedLine, budget: &mut Budget) -> io::Result<()> {
        if line.occurrences <= 1 {
            return Ok(());
        }
        budget.spend(line.text.len(), line.text.len() + 64)?;
        let mut message = Vec::with_capacity(line.text.len() + 64);
        message.push(b'<');
        message.extend_from_slice(&line.text);
        write!(&mut message, "> occurred {} times", line.occurrences)?;
        self.distribute(&message, budget)
    }
    fn complete(&mut self, text: Vec<u8>, budget: &mut Budget) -> io::Result<()> {
        if text.is_empty() {
            return self.distribute(&text, budget);
        }
        budget.spend(text.len().saturating_mul(self.cache.len()), 0)?;
        if let Some(index) = self.cache.iter().position(|line| line.text == text) {
            let mut line = self.cache.remove(index);
            line.occurrences = line
                .occurrences
                .checked_add(1)
                .ok_or_else(|| invalid("log repetition count overflow"))?;
            self.cache.push(line);
            return Ok(());
        }
        if self.cache.len() == 2 {
            let oldest = self.cache.remove(0);
            self.summary(&oldest, budget)?;
        }
        self.distribute(&text, budget)?;
        self.cache.push(CachedLine {
            text,
            occurrences: 1,
        });
        Ok(())
    }
    fn sync(&mut self, budget: &mut Budget) -> io::Result<()> {
        if !self.bound || self.routes.is_empty() {
            self.pending.clear();
            return Ok(());
        }
        let pending = std::mem::take(&mut self.pending);
        budget.spend(pending.len(), 0)?;
        for segment in pending.split_inclusive(|&byte| byte == b'\n') {
            let complete = segment.last() == Some(&b'\n');
            let bytes = if complete {
                &segment[..segment.len() - 1]
            } else {
                segment
            };
            let needed = self.incomplete.len() + bytes.len();
            if needed > MAX_LOG_LINE_BYTES {
                return Err(invalid("logging line length limit exceeded"));
            }
            if needed > self.incomplete.capacity() {
                budget.spend(needed, needed.saturating_mul(2))?;
            } else {
                budget.spend(bytes.len(), 0)?;
            }
            self.incomplete.extend_from_slice(bytes);
            if complete {
                let line = std::mem::take(&mut self.incomplete);
                self.complete(line, budget)?;
            }
        }
        Ok(())
    }
    fn distribute(&self, text: &[u8], budget: &mut Budget) -> io::Result<()> {
        static OUTPUT: Mutex<()> = Mutex::new(());
        budget.spend(
            self.routes.len(),
            self.routes.len() * std::mem::size_of::<(Vec<u8>, Option<LogNotification>)>(),
        )?;
        let mut records = Vec::with_capacity(self.routes.len());
        for route in &self.routes {
            let prefix = expand_prefix(&route.prefix, &self.level, self.clock.as_ref(), budget)?;
            let colored = self.color.filter(|_| route.sink.color_allowed());
            let bytes = prefix
                .len()
                .checked_add(text.len())
                .and_then(|v| v.checked_add(16))
                .ok_or_else(|| invalid("logging record overflow"))?;
            budget.spend(bytes, bytes)?;
            let mut record = Vec::with_capacity(bytes);
            if let Some(color) = colored {
                record.extend_from_slice(color.enable());
            }
            record.extend_from_slice(&prefix);
            record.extend_from_slice(text);
            if let Some(color) = colored {
                record.extend_from_slice(color.disable());
            }
            record.push(b'\n');
            records.push((record, route.notification.clone()));
        }
        {
            // ponytail: one process-wide output lock matches source record serialization.
            let _output = OUTPUT.lock().map_err(|_| poisoned())?;
            for (route, (record, _)) in self.routes.iter().zip(&records) {
                let mut writer = route.sink.writer.lock().map_err(|_| poisoned())?;
                writer.write_all(record)?;
                writer.flush()?;
            }
        }
        for (record, notification) in records {
            if let Some(notify) = notification {
                notify(&record)?;
            }
        }
        Ok(())
    }
}
impl LogStream {
    fn write_budgeted(&mut self, bytes: &[u8], budget: &mut Budget) -> io::Result<usize> {
        self.ready()?;
        if !self.bound {
            return Ok(bytes.len());
        }
        if bytes.len() > MAX_LOG_OPERATION_BYTES {
            return Err(invalid("logging input limit exceeded"));
        }
        budget.spend(bytes.len(), bytes.len())?;
        let result = (|| {
            let mut remaining = bytes;
            while !remaining.is_empty() {
                let size = remaining.len().min(LOG_BUFFER_BYTES - self.pending.len());
                self.pending.extend_from_slice(&remaining[..size]);
                remaining = &remaining[size..];
                if self.pending.len() == LOG_BUFFER_BYTES {
                    self.sync(budget)?;
                }
            }
            Ok(bytes.len())
        })();
        self.finish_operation(result)
    }
}
impl Write for LogStream {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.write_budgeted(bytes, &mut Budget::default())
    }
    fn write_fmt(&mut self, arguments: fmt::Arguments<'_>) -> io::Result<()> {
        struct SharedBudget<'a> {
            logger: &'a mut LogStream,
            budget: Budget,
        }
        impl Write for SharedBudget<'_> {
            fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
                self.logger.write_budgeted(bytes, &mut self.budget)
            }
            fn flush(&mut self) -> io::Result<()> {
                self.logger.sync(&mut self.budget)
            }
        }
        let result = SharedBudget {
            logger: self,
            budget: Budget::default(),
        }
        .write_fmt(arguments);
        self.finish_operation(result)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.ready()?;
        let result = self.sync(&mut Budget::default());
        self.finish_operation(result)
    }
}
impl Drop for LogStream {
    fn drop(&mut self) {
        if !self.failed {
            let _ = self.finish();
        }
    }
}

fn local_time() -> io::Result<LogTime> {
    let now = chrono::Local::now();
    Ok(LogTime {
        year: u16::try_from(now.year()).map_err(|_| invalid("local logging year out of range"))?,
        month: now.month() as u8,
        day: now.day() as u8,
        hour: now.hour() as u8,
        minute: now.minute() as u8,
        second: now.second() as u8,
    })
}

fn check_text(text: &str) -> io::Result<()> {
    if text.len() > MAX_LOG_TEXT_BYTES {
        Err(invalid("logging text length limit exceeded"))
    } else {
        Ok(())
    }
}
fn expand_prefix(
    prefix: &str,
    level: &str,
    clock: Option<&LogClock>,
    budget: &mut Budget,
) -> io::Result<Vec<u8>> {
    let bytes = prefix.as_bytes();
    budget.spend(bytes.len(), 0)?;
    let mut at = 0;
    let mut length = 0usize;
    while at < bytes.len() {
        let added = if bytes[at] != b'%' {
            at += 1;
            1
        } else {
            let Some(&token) = bytes.get(at + 1) else {
                break;
            };
            at += 2;
            match token {
                b'%' => 1,
                b'y' => level.len(),
                b'T' => 8,
                b't' | b'd' => 5,
                b'D' => 10,
                b'S' => 20,
                b's' => 12,
                _ => 0,
            }
        };
        length = length
            .checked_add(added)
            .ok_or_else(|| invalid("log prefix length overflow"))?;
    }
    budget.spend(length + bytes.len(), length)?;
    let mut output = Vec::with_capacity(length);
    let mut at = 0;
    let mut time = None;
    while at < bytes.len() {
        if bytes[at] != b'%' {
            output.push(bytes[at]);
            at += 1;
            continue;
        }
        let Some(&token) = bytes.get(at + 1) else {
            break;
        };
        at += 2;
        match token {
            b'%' => output.push(b'%'),
            b'y' => output.extend_from_slice(level.as_bytes()),
            b'T' | b't' | b'D' | b'd' | b'S' | b's' => {
                let t = match time {
                    Some(t) => t,
                    None => {
                        let provider = clock.ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::Unsupported,
                                "local-time log prefixes require a LogClock provider",
                            )
                        })?;
                        let t = provider()?.validate()?;
                        time = Some(t);
                        t
                    }
                };
                match token {
                    b'T' => write!(&mut output, "{:02}:{:02}:{:02}", t.hour, t.minute, t.second)?,
                    b't' => write!(&mut output, "{:02}:{:02}", t.hour, t.minute)?,
                    b'D' => write!(&mut output, "{:04}/{:02}/{:02}", t.year, t.month, t.day)?,
                    b'd' => write!(&mut output, "{:02}/{:02}", t.month, t.day)?,
                    b'S' => write!(
                        &mut output,
                        "{:04}/{:02}/{:02}, {:02}:{:02}:{:02}",
                        t.year, t.month, t.day, t.hour, t.minute, t.second
                    )?,
                    _ => write!(
                        &mut output,
                        "{:02}/{:02}, {:02}:{:02}",
                        t.month, t.day, t.hour, t.minute
                    )?,
                }
            }
            _ => (),
        }
    }
    Ok(output)
}

/// Temporarily detaches one sink. Like C++, restoration inserts a fresh route at the
/// end, with no prefix or notifier. Dereference the guard to write while suppressed.
pub struct LogSinkGuard<'a> {
    logger: &'a mut LogStream,
    sink: LogSink,
    removed: bool,
}
impl<'a> LogSinkGuard<'a> {
    pub fn new(logger: &'a mut LogStream, sink: &LogSink) -> io::Result<Self> {
        let removed = logger.has_stream(sink);
        if removed {
            logger.flush_incomplete()?;
            logger.remove(sink)?;
        }
        Ok(Self {
            logger,
            sink: sink.clone(),
            removed,
        })
    }
    pub fn restore(mut self) -> io::Result<()> {
        self.restore_inner()
    }
    fn restore_inner(&mut self) -> io::Result<()> {
        if !self.removed {
            return Ok(());
        }
        self.removed = false;
        let flushed = self.logger.flush_incomplete();
        let inserted = self.logger.insert(&self.sink);
        flushed.and(inserted)
    }
}
impl Deref for LogSinkGuard<'_> {
    type Target = LogStream;
    fn deref(&self) -> &Self::Target {
        self.logger
    }
}
impl DerefMut for LogSinkGuard<'_> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.logger
    }
}
impl Drop for LogSinkGuard<'_> {
    fn drop(&mut self) {
        let _ = self.restore_inner();
    }
}

/// Five independently routed levels; their labels are not numeric filtering priorities.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LogLevel {
    Fatal,
    Error,
    Warn,
    Info,
    Debug,
}
impl LogLevel {
    fn index(self) -> usize {
        self as usize
    }
    pub fn label(self) -> &'static str {
        match self {
            Self::Fatal => "FATAL_ERROR",
            Self::Error => "ERROR",
            Self::Warn => "WARNING",
            Self::Info => "INFO",
            Self::Debug => "DEBUG",
        }
    }
    fn configured(self) -> LogStream {
        let mut logger = LogStream::new(self.label()).expect("constant log level");
        logger.color = match self {
            Self::Fatal | Self::Error => Some(LogColor::Red),
            Self::Warn => Some(LogColor::Yellow),
            Self::Debug => Some(LogColor::Magenta),
            Self::Info => None,
        };
        let sink = match self {
            Self::Fatal | Self::Error | Self::Warn => Some(LogSink::stderr()),
            Self::Info => Some(LogSink::stdout()),
            Self::Debug => None,
        };
        if let Some(sink) = sink {
            logger.insert(&sink).expect("one logging sink");
        }
        logger
    }
}
fn global(level: LogLevel) -> &'static RwLock<LogStream> {
    static LOGS: [OnceLock<RwLock<LogStream>>; 5] = [const { OnceLock::new() }; 5];
    LOGS[level.index()].get_or_init(|| RwLock::new(level.configured()))
}
/// Configure a global route before threads first access it. Existing TLS streams retain
/// their earlier snapshot. Recursive access to the same stream returns `WouldBlock`.
pub fn with_global_log<T>(
    level: LogLevel,
    action: impl FnOnce(&mut LogStream) -> io::Result<T>,
) -> io::Result<T> {
    let mut logger = global(level).try_write().map_err(|error| match error {
        std::sync::TryLockError::WouldBlock => {
            io::Error::new(io::ErrorKind::WouldBlock, "global logger is in use")
        }
        std::sync::TryLockError::Poisoned(_) => poisoned(),
    })?;
    action(&mut logger)
}
thread_local! {
    static LOCAL: [RefCell<Option<LogStream>>; 5] = std::array::from_fn(|_| RefCell::new(None));
}
/// Access the current thread's independently buffered and cached global-route snapshot.
pub fn with_thread_local_log<T>(
    level: LogLevel,
    action: impl FnOnce(&mut LogStream) -> io::Result<T>,
) -> io::Result<T> {
    LOCAL
        .try_with(|logs| {
            let mut slot = logs[level.index()].try_borrow_mut().map_err(|_| {
                io::Error::new(io::ErrorKind::WouldBlock, "thread-local logger is in use")
            })?;
            if slot.is_none() {
                let configuration = global(level).try_read().map_err(|error| match error {
                    std::sync::TryLockError::WouldBlock => io::Error::new(
                        io::ErrorKind::WouldBlock,
                        "global logger is being configured",
                    ),
                    std::sync::TryLockError::Poisoned(_) => poisoned(),
                })?;
                *slot = Some(configuration.clone_configuration());
            }
            action(slot.as_mut().expect("initialized thread-local logger"))
        })
        .map_err(|_| {
            io::Error::new(
                io::ErrorKind::BrokenPipe,
                "thread-local logging is being destroyed",
            )
        })?
}
/// Write and flush one formatted line. Fatal logging does not exit the process.
pub fn log_message(level: LogLevel, message: fmt::Arguments<'_>) -> io::Result<()> {
    with_thread_local_log(level, |logger| {
        logger.write_fmt(message)?;
        logger.write_all(b"\n")?;
        logger.flush()
    })
}
/// Prefix a source location. Debug uses the basename; fatal retains the supplied path.
pub fn log_message_at(
    level: LogLevel,
    file: &str,
    line: u32,
    message: fmt::Arguments<'_>,
) -> io::Result<()> {
    let file = if level == LogLevel::Debug {
        file.rsplit(['/', '\\']).next().unwrap_or(file)
    } else {
        file
    };
    with_thread_local_log(level, |logger| {
        writeln!(logger, "{file}({line}): {message}")?;
        logger.flush()
    })
}

/// Log and flush one informational line, returning any writer/callback error.
#[macro_export]
macro_rules! openms_log_info {
    ($($arg:tt)*) => { $crate::concept::log_stream::log_message($crate::concept::log_stream::LogLevel::Info, format_args!($($arg)*)) };
}
/// Log and flush one warning line.
#[macro_export]
macro_rules! openms_log_warn {
    ($($arg:tt)*) => { $crate::concept::log_stream::log_message($crate::concept::log_stream::LogLevel::Warn, format_args!($($arg)*)) };
}
/// Log and flush one error line.
#[macro_export]
macro_rules! openms_log_error {
    ($($arg:tt)*) => { $crate::concept::log_stream::log_message($crate::concept::log_stream::LogLevel::Error, format_args!($($arg)*)) };
}
/// Log and flush one fatal line with its source location; does not exit the process.
#[macro_export]
macro_rules! openms_log_fatal_error {
    ($($arg:tt)*) => { $crate::concept::log_stream::log_message_at($crate::concept::log_stream::LogLevel::Fatal, file!(), line!(), format_args!($($arg)*)) };
}
/// Log and flush one debug line with its source basename and line.
#[macro_export]
macro_rules! openms_log_debug {
    ($($arg:tt)*) => { $crate::concept::log_stream::log_message_at($crate::concept::log_stream::LogLevel::Debug, file!(), line!(), format_args!($($arg)*)) };
}
/// Log and flush one debug line without a source location.
#[macro_export]
macro_rules! openms_log_debug_nofile {
    ($($arg:tt)*) => { $crate::concept::log_stream::log_message($crate::concept::log_stream::LogLevel::Debug, format_args!($($arg)*)) };
}
