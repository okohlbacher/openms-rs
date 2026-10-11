// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The info log stream of one tool run, as far as text reaches it without a
//! line break: the source's `OPENMS_LOG_INFO << text` (core
//! `CONCEPT/LogStream.h`, `LogStream.cpp:180-300` at `bc9cc12`).
//!
//! The source's stream collects text until it is flushed, splits it into
//! lines, and passes each through a two-line cache: a line equal to one of the
//! two most recently printed is counted instead of printed, and when a new
//! line pushes the older of the two out, `<line> occurred N times` is printed
//! first if it was repeated. Empty lines bypass the cache. At the end of the
//! process, every cached line that was repeated is reported, in lexicographic
//! order, and text without a final line break is printed as a line.
//!
//! `TOPPBase` writes the closing `<tool> took … .` line to the same stream, so
//! it too passes the cache, and joins a line the tool left unfinished. This
//! type models that for the text a tool writes with
//! [`ToolContext::write_info_stream`](super::ToolContext::write_info_stream);
//! the lifecycle closes it with the closing line. A run that never writes to
//! it prints the closing line directly, as before.
//!
//! The stream is the crate's [`LogStream`] (the port of `LogStreamBuf`, with
//! its cache and its ceilings), bound to a buffer that each call moves on to
//! the run's `out`.

use crate::concept::log_stream::{LogSink, LogStream};
use std::io::{self, Write};
use std::sync::{Arc, Mutex};

fn poisoned() -> io::Error {
    io::Error::other("the info log stream is poisoned")
}

/// The one destination of the run's [`LogStream`]: a buffer shared with the
/// [`InfoStream`], which moves every record the stream delivers on to the
/// run's `out`.
#[derive(Clone, Debug, Default)]
struct Records(Arc<Mutex<Vec<u8>>>);

impl Records {
    fn take(&self) -> io::Result<Vec<u8>> {
        let mut records = self.0.lock().map_err(|_| poisoned())?;
        Ok(std::mem::take(&mut *records))
    }
}

impl Write for Records {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0
            .lock()
            .map_err(|_| poisoned())?
            .extend_from_slice(bytes);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// The run's info log stream: a [`LogStream`] of level `INFO` without a
/// prefix, whose records go to the `out` of the call that completes them.
#[derive(Debug, Default)]
pub(crate) struct InfoStream {
    used: bool,
    logger: Option<LogStream>,
    records: Records,
}

impl InfoStream {
    /// The stream, bound on first use.
    fn logger(&mut self) -> io::Result<&mut LogStream> {
        if self.logger.is_none() {
            let mut logger = LogStream::new("INFO")?;
            logger.insert(&LogSink::new(self.records.clone()))?;
            self.logger = Some(logger);
        }
        self.logger
            .as_mut()
            .ok_or_else(|| io::Error::other("the info log stream is not bound"))
    }

    /// Append `text`; every line it completes goes through the cache to
    /// `out`.
    pub(crate) fn write(&mut self, out: &mut dyn Write, text: &str) -> io::Result<()> {
        self.used = true;
        let logger = self.logger()?;
        logger.write_all(text.as_bytes())?;
        logger.flush()?;
        out.write_all(&self.records.take()?)
    }

    /// Whether anything was written to the stream in this run.
    pub(crate) fn used(&self) -> bool {
        self.used
    }

    /// The end of the run: `closing` (the lifecycle's `<tool> took … .` line,
    /// when the body returned) as the last line written to the stream, then
    /// the stream's end-of-process output ([`LogStream::finish`]): the
    /// repeated cached lines in lexicographic order, then unfinished text as a
    /// line.
    pub(crate) fn finish(&mut self, out: &mut dyn Write, closing: Option<&str>) -> io::Result<()> {
        if let Some(closing) = closing {
            self.write(out, closing)?;
            self.write(out, "\n")?;
        }
        if let Some(logger) = self.logger.as_mut() {
            logger.finish()?;
        }
        out.write_all(&self.records.take()?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(chunks: &[&str], closing: Option<&str>) -> String {
        let mut stream = InfoStream::default();
        let mut out = Vec::new();
        for chunk in chunks {
            stream.write(&mut out, chunk).unwrap();
        }
        stream.finish(&mut out, closing).unwrap();
        String::from_utf8(out).unwrap()
    }

    /// Release build, oracle `ext_repeat_debug4` of `../oracle/B7-cli`: the
    /// child's `a`, `a`, `b`, `a` and the banner lines through the info
    /// stream, then the closing line.
    #[test]
    fn repeated_lines_are_counted_and_reported_on_eviction() {
        let text = run(
            &[
                "Running: sh -c echo a; echo a; echo b; echo a\n",
                "a\na\nb\na\n",
                "Executed 'sh' successfully!\n",
            ],
            Some("T took <T>."),
        );
        assert_eq!(
            text,
            "Running: sh -c echo a; echo a; echo b; echo a\na\nb\nExecuted 'sh' successfully!\n<a> occurred 3 times\nT took <T>.\n"
        );
    }

    #[test]
    fn unfinished_text_joins_the_next_line_and_ends_as_a_line() {
        assert_eq!(run(&["abc", "def\n", "x"], None), "abcdef\nx\n");
        assert_eq!(run(&["abc"], Some("T took.")), "abcT took.\n");
        assert_eq!(run(&["\n", "\n"], None), "\n\n");
    }

    #[test]
    fn remaining_repeats_are_reported_in_lexical_order_at_the_end() {
        assert_eq!(
            run(&["b\nb\na\na\n"], None),
            "b\na\n<a> occurred 2 times\n<b> occurred 2 times\n"
        );
    }
}
