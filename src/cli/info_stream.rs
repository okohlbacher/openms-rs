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

use std::io::{self, Write};

/// One cached line: its text, how often it occurred, and when it was last
/// seen.
#[derive(Clone, Debug)]
struct Cached {
    text: String,
    occurrences: u64,
    last_use: u64,
}

/// The pending text and line cache of the run's info stream.
#[derive(Clone, Debug, Default)]
pub(crate) struct InfoStream {
    used: bool,
    pending: String,
    cache: Vec<Cached>,
    clock: u64,
}

impl InfoStream {
    /// Append `text`; every line it completes goes through the cache to
    /// `out`.
    pub(crate) fn write(&mut self, out: &mut dyn Write, text: &str) -> io::Result<()> {
        self.used = true;
        self.pending.push_str(text);
        while let Some(end) = self.pending.find('\n') {
            let line: String = self.pending.drain(..=end).collect();
            self.line(out, &line[..line.len() - 1])?;
        }
        Ok(())
    }

    /// Whether anything was written to the stream in this run.
    pub(crate) fn used(&self) -> bool {
        self.used
    }

    /// One complete line through the cache (`LogStreamBuf::sync` and
    /// `addToCache_`).
    fn line(&mut self, out: &mut dyn Write, text: &str) -> io::Result<()> {
        if text.is_empty() {
            return writeln!(out);
        }
        self.clock += 1;
        if let Some(entry) = self.cache.iter_mut().find(|entry| entry.text == text) {
            entry.occurrences += 1;
            entry.last_use = self.clock;
            return Ok(());
        }
        if self.cache.len() >= 2 {
            let oldest = (0..self.cache.len())
                .min_by_key(|&index| self.cache[index].last_use)
                .unwrap_or(0);
            let evicted = self.cache.remove(oldest);
            if evicted.occurrences > 1 {
                writeln!(
                    out,
                    "<{}> occurred {} times",
                    evicted.text, evicted.occurrences
                )?;
            }
        }
        writeln!(out, "{text}")?;
        self.cache.push(Cached {
            text: text.to_owned(),
            occurrences: 1,
            last_use: self.clock,
        });
        Ok(())
    }

    /// The end of the run: `closing` (the lifecycle's `<tool> took … .` line,
    /// when the body returned) as the last line written to the stream, then
    /// the stream's end-of-process output — the repeated cached lines in
    /// lexicographic order, then unfinished text as a line.
    pub(crate) fn finish(&mut self, out: &mut dyn Write, closing: Option<&str>) -> io::Result<()> {
        if let Some(closing) = closing {
            self.write(out, closing)?;
            self.write(out, "\n")?;
        }
        self.cache.sort_by(|a, b| a.text.cmp(&b.text));
        for entry in std::mem::take(&mut self.cache) {
            if entry.occurrences > 1 {
                writeln!(out, "<{}> occurred {} times", entry.text, entry.occurrences)?;
            }
        }
        if !self.pending.is_empty() {
            let rest = std::mem::take(&mut self.pending);
            writeln!(out, "{rest}")?;
        }
        Ok(())
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
