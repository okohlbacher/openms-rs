// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! What a tool's text looks like on a console: the width the usage text is
//! shaped to (`ConsoleUtils::readConsoleSize_`, `breakString_`,
//! `ConsoleUtils.cpp:39-196`, core `bc9cc12`) and the colours
//! (`Colorizer.cpp`, and the red and yellow of the error and warning log
//! streams, `LogStream.cpp:256-287`, `:667-669`).
//!
//! The source asks the process's own streams: `COLUMNS`, else `stty size`
//! with the process's standard input, and `isatty` on standard output and
//! error. This port asks them only while a tool runs as its own executable
//! ([`with_process_streams`], which [`run`](crate::cli::run) uses). A caller
//! that drives a tool in process through
//! [`run_with`](crate::cli::run_with) passes streams that are not a console,
//! so the text is the one the source writes when standard error is not a
//! terminal and `COLUMNS` is unset: no colours and no shaping.

use crate::concept::log_stream::LogColor;
use std::cell::Cell;
use std::io::{BufRead, IsTerminal, Read, Result as IoResult, Write};
use std::sync::OnceLock;

/// The width that turns shaping off, the source's
/// `std::numeric_limits<int>::max()`.
pub(crate) const UNSHAPED: i64 = i32::MAX as i64;

/// Most bytes of `stty size` output the source reads, `fgets(buff, 100)`.
const STTY_LINE_BYTES: u64 = 99;

thread_local! {
    /// Whether the tool on this thread writes to the process's own standard
    /// output and error.
    static PROCESS_STREAMS: Cell<bool> = const { Cell::new(false) };
}

/// The width the process's console reports, measured once, when the first
/// usage text is written, as the source's `ConsoleUtils` singleton does.
static PROCESS_WIDTH: OnceLock<i64> = OnceLock::new();

/// Run `body` as a tool executable runs: its streams are the process's
/// standard output and error, so the usage text is shaped to the console and
/// coloured on a terminal. Afterwards, as the source's `InitConsole`
/// destructor at process exit, every terminal among the two streams gets the
/// reset `\x1b[0m`.
pub(crate) fn with_process_streams<R>(body: impl FnOnce() -> R) -> R {
    let previous = PROCESS_STREAMS.with(|flag| flag.replace(true));
    let result = body();
    PROCESS_STREAMS.with(|flag| flag.set(previous));
    if !previous {
        reset_terminals();
    }
    result
}

/// Source `InitConsole::~InitConsole` (`Colorizer.cpp:30-40`): `undoAll` on
/// `cout`, then on `cerr`, each written only to a terminal.
fn reset_terminals() {
    if std::io::stdout().is_terminal() {
        let mut out = std::io::stdout().lock();
        let _ = out.write_all(b"\x1b[0m");
        let _ = out.flush();
    }
    if std::io::stderr().is_terminal() {
        let _ = std::io::stderr().lock().write_all(b"\x1b[0m");
    }
}

fn process_streams() -> bool {
    PROCESS_STREAMS.with(Cell::get)
}

/// Whether text written to the error stream is coloured: the tool runs on
/// the process's streams and standard error is a terminal, as the source's
/// `Colorizer::isTTY(std::cerr)`.
fn error_stream_coloured() -> bool {
    process_streams() && std::io::stderr().is_terminal()
}

/// How the usage text is laid out: the console width (at least 10, or
/// [`UNSHAPED`]) and whether it is coloured.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Console {
    pub(crate) width: i64,
    pub(crate) colour: bool,
}

impl Console {
    /// The layout of text written through explicit streams: no shaping, no
    /// colours.
    pub(crate) const PLAIN: Self = Self {
        width: UNSHAPED,
        colour: false,
    };

    /// The layout of the usage text of the running tool. On the process's
    /// streams this measures the console the first time it is asked, which
    /// may run `stty size` (and let it print its complaint on standard error,
    /// as the source's does).
    pub(crate) fn for_usage() -> Self {
        if !process_streams() {
            return Self::PLAIN;
        }
        Self {
            width: process_width(),
            colour: error_stream_coloured(),
        }
    }
}

/// Source `ConsoleUtils::readConsoleSize_` for a `COLUMNS` value and a way to
/// run `stty size`.
///
/// `COLUMNS`, when set, must be an integer as `StringUtils::toInt32` reads
/// one; otherwise the second of exactly two space-separated fields of the
/// first line `stty size` prints. The width is that number less one, and any
/// width below 10, including an unreadable one, turns shaping off.
pub(crate) fn width_from(
    columns: Option<&std::ffi::OsStr>,
    stty: impl FnOnce() -> Option<Vec<u8>>,
) -> i64 {
    let measured = match columns {
        Some(value) => value
            .to_str()
            .and_then(|text| super::context::to_int32(text).ok()),
        None => stty().and_then(|line| {
            let line = String::from_utf8_lossy(&line);
            let fields: Vec<&str> = line.split(' ').collect();
            match fields.as_slice() {
                [_, columns] => super::context::to_int32(columns).ok(),
                _ => None,
            }
        }),
    };
    // The source decrements what it read; `INT_MIN - 1` is undefined there
    // and wraps here.
    let width = measured.map_or(-1, |columns| i64::from(columns.wrapping_sub(1)));
    if width < 10 { UNSHAPED } else { width }
}

/// Source `popen("stty size", "r")` and one `fgets`: run `stty size` through
/// `/bin/sh` with this process's standard input and standard error, and
/// return the first line it prints, at most 99 bytes.
#[cfg(unix)]
fn stty_size() -> Option<Vec<u8>> {
    use std::process::{Command, Stdio};
    let mut child = Command::new("/bin/sh")
        .args(["-c", "stty size"])
        .stdin(Stdio::inherit())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .spawn()
        .ok()?;
    let mut line = Vec::new();
    if let Some(output) = child.stdout.take() {
        let mut reader = std::io::BufReader::new(output.take(STTY_LINE_BYTES));
        let _ = reader.read_until(b'\n', &mut line);
    }
    let _ = child.wait();
    (!line.is_empty()).then_some(line)
}

/// The source asks the Windows console for its buffer width, which the
/// standard library cannot; without `COLUMNS` the text is not shaped there.
#[cfg(not(unix))]
fn stty_size() -> Option<Vec<u8>> {
    None
}

/// Source `ConsoleUtils::breakString_` (`ConsoleUtils.cpp:122-194`) on bytes,
/// for a console `width` ([`UNSHAPED`] or at least 10).
///
/// The first line takes what is left of the console after
/// `first_line_prefill` bytes, every later line `width - indentation` bytes
/// after `indentation` spaces, and a line break in `input` ends a line early.
/// A line that fills its width and ends in a word shorter than four bytes
/// gives that word to the next line. An input ending in a line break ends in
/// an indented empty line, and more than `max_lines` lines keep the first
/// `max_lines - 2`, an indented `...` and the last one. When the indentation
/// equals the width, the input is returned unbroken; when it exceeds it, the
/// source's unsigned difference wraps and later lines are not shortened.
pub(crate) fn break_string(
    input: &[u8],
    indentation: usize,
    max_lines: usize,
    first_line_prefill: usize,
    width: i64,
) -> Vec<Vec<u8>> {
    let mut result: Vec<Vec<u8>> = Vec::new();
    if input.is_empty() {
        return result;
    }
    // `Size short_line_len = console_width_ - indentation;`
    let width_bytes = u64::try_from(width).unwrap_or(0);
    let short_line_len = width_bytes.wrapping_sub(indentation as u64);
    if short_line_len < 1 {
        result.push(input.to_vec());
        return result;
    }
    let mut prefill = first_line_prefill as u64;
    // `(int)first_line_prefill > console_width_`
    if i64::from(prefill as i32) > width && width_bytes != 0 {
        prefill %= width_bytes;
    }
    let mut position = 0usize;
    while position < input.len() {
        let remaining = if result.is_empty() {
            width_bytes.wrapping_sub(prefill)
        } else {
            short_line_len
        };
        let indent = if result.is_empty() { 0 } else { indentation };
        let rest = &input[position..];
        let take = usize::try_from(remaining).map_or(rest.len(), |n| n.min(rest.len()));
        let mut line = &rest[..take];
        let mut advance = line.len();
        if let Some(at) = line.iter().position(|&byte| byte == b'\n') {
            line = &line[..at];
            advance = at + 1;
        }
        if line.len() as u64 == remaining && short_line_len > 8 {
            if let Some(space) = line.iter().rposition(|&byte| byte == b' ') {
                let last_word = line.len() - space - 1;
                if last_word < 4 {
                    line = &line[..line.len() - last_word];
                    advance -= last_word;
                }
            }
        }
        position += advance;
        let mut shaped = vec![b' '; indent];
        shaped.extend_from_slice(line);
        result.push(shaped);
    }
    if input.last() == Some(&b'\n') {
        result.push(vec![b' '; indentation]);
    }
    if result.len() > max_lines && max_lines >= 2 {
        let last = result.pop().unwrap_or_default();
        result.truncate(max_lines - 2);
        let mut ellipsis = vec![b' '; indentation];
        ellipsis.extend_from_slice(b"...");
        result.push(ellipsis);
        result.push(last);
    }
    result
}

/// Terminal-width-aware string wrapping for console output: source
/// `APPLICATIONS/ConsoleUtils.h`.
///
/// The source is a singleton that detects the console width once per process,
/// on construction, and wraps long strings to it with a hanging indentation,
/// the pattern TOPP tools use for `--help` output, parameter descriptions and
/// error messages. The width is `COLUMNS` when that is set (read as
/// `StringUtils::toInt32` reads it), otherwise the second field of `stty size`
/// run with the process's standard input; one less than that; and
/// `i32::MAX`, which disables wrapping, when it is below 10 or cannot be
/// determined. Later resizes of the terminal are not seen.
///
/// | Source | Rust |
/// |---|---|
/// | `static const ConsoleUtils& getInstance()` | [`ConsoleUtils::get_instance`] |
/// | private `ConsoleUtils()`, deleted copy and assignment | private construction, no `Clone` |
/// | `int getConsoleWidth() const` | [`ConsoleUtils::get_console_width`] |
/// | `static StringList breakStringList(input, indentation, max_lines, first_line_prefill = 0)` | [`ConsoleUtils::break_string_list`] (bytes: [`ConsoleUtils::break_bytes_list`]) |
/// | `static std::string breakString(input, indentation, max_lines, first_line_prefill = 0)` | [`ConsoleUtils::break_string`] |
/// | `friend struct ConsoleWidthTest` (pins the width in the class test) | [`ConsoleUtils::with_console_width`] and the instance methods |
/// | private `readConsoleSize_`, `breakString_`, `console_width_` | crate-private `width_from` and `break_string`; a private field |
///
/// # Windows
///
/// On Windows the source reads the console buffer width with
/// `GetConsoleScreenBufferInfo`. The standard library has no safe
/// counterpart and the crate forbids `unsafe`, so there only `COLUMNS` is
/// read and wrapping is otherwise off.
///
/// # Examples
///
/// The upstream class test pins the width to 9 columns:
///
/// ```
/// use openms::cli::ConsoleUtils;
///
/// let console = ConsoleUtils::with_console_width(9);
/// let x20 = "x".repeat(19);
/// assert_eq!(console.break_text(&x20, 3, 10, 0), "xxxxxxxxx\n   xxxxxx\n   xxxx");
/// ```
#[derive(Debug, PartialEq, Eq)]
pub struct ConsoleUtils {
    /// Source `console_width_`: at least 10, or `i32::MAX` for no wrapping.
    console_width: i32,
}

impl ConsoleUtils {
    /// The process-wide instance, whose width is measured on the first call:
    /// source `getInstance()`.
    ///
    /// Measuring may run `stty size` through `/bin/sh` (on Unix), which prints
    /// its own complaint on standard error when standard input is not a
    /// terminal, as the source's does. Tools running as their own executable
    /// share this measurement with their usage text.
    pub fn get_instance() -> &'static Self {
        static INSTANCE: OnceLock<ConsoleUtils> = OnceLock::new();
        INSTANCE.get_or_init(|| Self {
            console_width: i32::try_from(process_width()).unwrap_or(i32::MAX),
        })
    }

    /// An instance with a fixed width, as the class test's `ConsoleWidthTest`
    /// hook sets one: the width is taken as given, while a measured width is
    /// never below 10. Line breaking follows the source's `breakString_` for
    /// any width, including its unsigned wrap-around when `indentation`
    /// exceeds the width.
    pub fn with_console_width(width: i32) -> Self {
        Self {
            console_width: width,
        }
    }

    /// The cached console width: source `getConsoleWidth()`. `i32::MAX`, the
    /// source's `std::numeric_limits<int>::max()`, means that the width could
    /// not be detected and wrapping is disabled.
    pub fn get_console_width(&self) -> i32 {
        self.console_width
    }

    /// Wraps `input` to the console width, one element per output line:
    /// source `breakStringList`.
    ///
    /// Lines 2 to N are prefixed by `indentation` spaces, giving a
    /// left-aligned hanging block, and lines wider than the console are broken
    /// at its width. An `indentation` of 0 falls back to the native console
    /// behaviour: a raw break at the right edge with no prefix. More than
    /// `max_lines` lines keep the first `max_lines - 2`, a `...` line and the
    /// final line, so the tail of the message stays visible.
    /// `first_line_prefill` is the number of characters already on the
    /// current console line, e.g. a label this string follows, so the first
    /// break falls correctly. A line that fills the width and ends in a word
    /// shorter than four bytes gives that word to the next line, and an input
    /// ending in a line break ends in an indented empty line.
    ///
    /// Breaking is on bytes, as the source's `std::string` is; a multi-byte
    /// UTF-8 character split at a line end is replaced by U+FFFD here (use
    /// [`Self::break_bytes_list`] for the exact bytes).
    pub fn break_string_list(
        input: &str,
        indentation: usize,
        max_lines: usize,
        first_line_prefill: usize,
    ) -> Vec<String> {
        Self::get_instance().break_lines(input, indentation, max_lines, first_line_prefill)
    }

    /// [`Self::break_string_list`] over bytes, exactly as the source breaks a
    /// `std::string`.
    pub fn break_bytes_list(
        input: &[u8],
        indentation: usize,
        max_lines: usize,
        first_line_prefill: usize,
    ) -> Vec<Vec<u8>> {
        Self::get_instance().break_byte_lines(input, indentation, max_lines, first_line_prefill)
    }

    /// [`Self::break_string_list`] joined with `\n`: source `breakString`.
    pub fn break_string(
        input: &str,
        indentation: usize,
        max_lines: usize,
        first_line_prefill: usize,
    ) -> String {
        Self::get_instance().break_text(input, indentation, max_lines, first_line_prefill)
    }

    /// [`Self::break_string_list`] at this instance's width.
    pub fn break_lines(
        &self,
        input: &str,
        indentation: usize,
        max_lines: usize,
        first_line_prefill: usize,
    ) -> Vec<String> {
        self.break_byte_lines(input.as_bytes(), indentation, max_lines, first_line_prefill)
            .into_iter()
            .map(|line| String::from_utf8_lossy(&line).into_owned())
            .collect()
    }

    /// [`Self::break_bytes_list`] at this instance's width.
    pub fn break_byte_lines(
        &self,
        input: &[u8],
        indentation: usize,
        max_lines: usize,
        first_line_prefill: usize,
    ) -> Vec<Vec<u8>> {
        break_string(
            input,
            indentation,
            max_lines,
            first_line_prefill,
            i64::from(self.console_width),
        )
    }

    /// [`Self::break_string`] at this instance's width.
    pub fn break_text(
        &self,
        input: &str,
        indentation: usize,
        max_lines: usize,
        first_line_prefill: usize,
    ) -> String {
        let lines =
            self.break_byte_lines(input.as_bytes(), indentation, max_lines, first_line_prefill);
        String::from_utf8_lossy(&lines.join(&b'\n')).into_owned()
    }
}

/// The process's console width, measured once: `COLUMNS`, else `stty size`.
fn process_width() -> i64 {
    *PROCESS_WIDTH.get_or_init(|| width_from(std::env::var_os("COLUMNS").as_deref(), stty_size))
}

/// Write `text` as the source's error log stream writes to `std::cerr`: each
/// line in red on a terminal (`\x1b[91m` … `\x1b[39m`), as it is otherwise.
///
/// This is what `writeLogError_` and `OPENMS_LOG_ERROR` print. The colour
/// applies only while the tool runs as its own executable and standard error
/// is a terminal.
///
/// # Errors
///
/// When `err` fails.
pub fn log_error(err: &mut dyn Write, text: &str) -> IoResult<()> {
    log_line(err, LogColor::Red, text)
}

/// Write `text` as the source's warning log stream writes to `std::cerr`:
/// each line in yellow on a terminal (`\x1b[93m` … `\x1b[39m`); see
/// [`log_error`].
///
/// # Errors
///
/// When `err` fails.
pub fn log_warning(err: &mut dyn Write, text: &str) -> IoResult<()> {
    log_line(err, LogColor::Yellow, text)
}

fn log_line(err: &mut dyn Write, colour: LogColor, text: &str) -> IoResult<()> {
    write_log_line(err, error_stream_coloured().then_some(colour), text)
}

/// One record as `LogStreamBuf::distribute_` writes it to a stream: each line
/// of `text` wrapped in `colour`'s codes when there is one, then a line break.
fn write_log_line(err: &mut dyn Write, colour: Option<LogColor>, text: &str) -> IoResult<()> {
    let Some(colour) = colour else {
        return writeln!(err, "{text}");
    };
    for line in text.split('\n') {
        err.write_all(colour.enable())?;
        err.write_all(line.as_bytes())?;
        err.write_all(colour.disable())?;
        err.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    //! The pieces of the console layout a test process can reach: the width
    //! rule, the line breaking and the log-line colours. The executables are
    //! compared with the Release build in `tests/topp_cli_console.rs`; the
    //! usage text on a terminal in `super::super::usage`'s tests.
    use super::*;
    use std::ffi::OsStr;

    fn tty_case(case: &str, file: &str) -> Vec<u8> {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/topp_cli_console")
            .join(case)
            .join(file);
        std::fs::read(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
    }

    /// `ConsoleUtils::readConsoleSize_`: COLUMNS wins without running
    /// `stty`; otherwise `stty size` of a 40x70 terminal gives 69 (oracle
    /// `tty_help_BaselineFilter` breaks at 69) and of a 0x0 one no shaping
    /// (`tty_help_nosize`); anything unreadable or below 10 turns shaping off.
    #[test]
    fn the_width_follows_the_source_rule() {
        let never = || -> Option<Vec<u8>> { panic!("stty must not run when COLUMNS is set") };
        assert_eq!(width_from(Some(OsStr::new("70")), never), 69);
        assert_eq!(width_from(Some(OsStr::new(" +45 ")), never), 44);
        assert_eq!(width_from(Some(OsStr::new("11")), never), 10);
        assert_eq!(width_from(Some(OsStr::new("10")), never), UNSHAPED);
        assert_eq!(width_from(Some(OsStr::new("")), never), UNSHAPED);
        assert_eq!(width_from(Some(OsStr::new("45x")), never), UNSHAPED);
        assert_eq!(width_from(Some(OsStr::new("-2147483648")), never), UNSHAPED);
        assert_eq!(width_from(None, || Some(b"40 70\n".to_vec())), 69);
        assert_eq!(width_from(None, || Some(b"0 0\n".to_vec())), UNSHAPED);
        assert_eq!(width_from(None, || Some(b"40  70\n".to_vec())), UNSHAPED);
        assert_eq!(width_from(None, || None), UNSHAPED);
    }

    /// `breakString_` on the source's own edge cases, read from its code:
    /// a short last word moves to the next line, a trailing line break leaves
    /// an indented empty line, and an eleventh line becomes `...`.
    #[test]
    fn breaking_follows_the_source_rules() {
        let lines = |input: &str, indentation, prefill, width| -> Vec<String> {
            break_string(input.as_bytes(), indentation, 10, prefill, width)
                .into_iter()
                .map(|line| String::from_utf8(line).unwrap())
                .collect()
        };
        assert_eq!(lines("", 0, 0, 20), Vec::<String>::new());
        assert_eq!(lines("test this break", 0, 0, 12), ["test this ", "break"]);
        assert_eq!(lines("test thisbreak", 2, 0, 12), ["test thisbre", "  ak"]);
        assert_eq!(lines("a\n", 3, 0, UNSHAPED), ["a", "   "]);
        assert_eq!(lines("abcdef", 4, 8, 10), ["ab", "    cdef"]);
        // An indentation equal to the width leaves the input whole.
        assert_eq!(lines("abc\ndef", 12, 0, 12), ["abc\ndef"]);
        let many = "1\n2\n3\n4\n5\n6\n7\n8\n9\n10\n11";
        assert_eq!(
            lines(many, 1, 0, UNSHAPED),
            ["1", " 2", " 3", " 4", " 5", " 6", " 7", " 8", " ...", " 11"]
        );
    }

    /// Error and warning records on a terminal, byte for byte as the Release
    /// build wrote them (oracles `tty_missing_in` and `tty_ini_unknown`).
    #[test]
    fn log_records_are_coloured_as_the_release_build_colours_them() {
        let mut written = Vec::new();
        write_log_line(
            &mut written,
            Some(LogColor::Red),
            "Error: The required parameter 'in' [valid: mzML] was not given or is empty!",
        )
        .unwrap();
        written.extend_from_slice(b"\x1b[0m");
        assert_eq!(written, tty_case("tty_missing_in", "tty.bin"));

        let mut written = Vec::new();
        write_log_line(
            &mut written,
            Some(LogColor::Red),
            "Parameters passed to 'BaselineFilter' are invalid. To prevent usage of wrong defaults, please update/fix the parameters!",
        )
        .unwrap();
        write_log_line(
            &mut written,
            Some(LogColor::Yellow),
            "Unknown (or deprecated) Parameter 'bogus_item' given in outdated parameter file!",
        )
        .unwrap();
        written.extend_from_slice(b"\x1b[0m");
        assert_eq!(written, tty_case("tty_ini_unknown", "tty.bin"));

        let mut plain = Vec::new();
        write_log_line(&mut plain, None, "a\nb").unwrap();
        assert_eq!(plain, b"a\nb\n");
    }
}
