// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Colored and styled console text: `CONCEPT/Colorizer.h`.
//!
//! A [`Colorizer`](crate::concept::colorizer::Colorizer) colors the console font's foreground (a new color, or the
//! current one brighter), underlines it or inverts it. Styles combine, with
//! the limitation that undoing does not nest: [`Colorizer::undo`](crate::concept::colorizer::Colorizer::undo) does not
//! remember a previous color but resets to the default. When standard output
//! or standard error is redirected away from a terminal, no ANSI codes are
//! written to it, so they never end up in a file; any other writer always
//! receives them.
//!
//! # API mapping
//!
//! | Source (`Colorizer.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `enum class ConsoleColor` (nine values) | [`ConsoleColor`](crate::concept::colorizer::ConsoleColor), the `concept::log_stream` color enum |
//! | `Colorizer(const ConsoleColor)` | [`Colorizer::new`](crate::concept::colorizer::Colorizer::new) |
//! | `Colorizer()`, `Colorizer(const Colorizer&)` (both deleted) | no `Default`, no `Clone` |
//! | `operator()()` | [`Colorizer::apply`](crate::concept::colorizer::Colorizer::apply) |
//! | `template <T> operator()(T s)` | [`Colorizer::apply_to`](crate::concept::colorizer::Colorizer::apply_to) |
//! | `undo()`, `undoAll()` | [`Colorizer::undo`](crate::concept::colorizer::Colorizer::undo), [`Colorizer::undo_all`](crate::concept::colorizer::Colorizer::undo_all) |
//! | `static isTTY(const std::ostream&)` | [`Colorizer::is_tty`](crate::concept::colorizer::Colorizer::is_tty) over [`ConsoleStream`](crate::concept::colorizer::ConsoleStream) |
//! | `friend operator<<(std::ostream&, Colorizer&)` | [`Colorizer::write_to`](crate::concept::colorizer::Colorizer::write_to); `Display` for any other stream |
//! | protected `outputToStream_`, `colorStream_` | inside [`Colorizer::write_to`](crate::concept::colorizer::Colorizer::write_to) |
//! | protected `colors_` table and `color_undo_all_` | [`Colorizer::enable_code`](crate::concept::colorizer::Colorizer::enable_code), [`Colorizer::disable_code`](crate::concept::colorizer::Colorizer::disable_code), [`COLOR_UNDO_ALL`](crate::concept::colorizer::COLOR_UNDO_ALL) |
//! | protected `getInternalChars_`, `setInternalChars_` (for `IndentedStream`) | [`Colorizer::internal_chars`](crate::concept::colorizer::Colorizer::internal_chars), [`Colorizer::set_internal_chars`](crate::concept::colorizer::Colorizer::set_internal_chars) |
//! | protected members `color_`, `undo_`, `undo_all_`, `undos_only`, `input_` | private fields; [`Colorizer::color`](crate::concept::colorizer::Colorizer::color) |
//! | `extern Colorizer red, green, yellow, blue, magenta, cyan, bright, underline, invert` | [`RED`](crate::concept::colorizer::RED), [`GREEN`](crate::concept::colorizer::GREEN), [`YELLOW`](crate::concept::colorizer::YELLOW), [`BLUE`](crate::concept::colorizer::BLUE), [`MAGENTA`](crate::concept::colorizer::MAGENTA), [`CYAN`](crate::concept::colorizer::CYAN), [`BRIGHT`](crate::concept::colorizer::BRIGHT), [`UNDERLINE`](crate::concept::colorizer::UNDERLINE), [`INVERT`](crate::concept::colorizer::INVERT) |
//!
//! # Preserved source conventions
//!
//! - A new colorizer has both undo flags set, as the source's member
//!   initialisers do. Written before any `apply`, it therefore emits its color,
//!   no text, and the reset-all code `\x1b[0m`.
//! - [`Colorizer::undo`](crate::concept::colorizer::Colorizer::undo) does not clear the undo-all flag. After
//!   [`Colorizer::undo_all`](crate::concept::colorizer::Colorizer::undo_all), or on a colorizer never applied, `undo` emits
//!   `\x1b[0m`, not the color's own disable code. This is why the global
//!   colorizers keep their state between uses.
//! - Writing a colorizer does not change it: writing it twice emits the same
//!   codes and text twice.
//!
//! # Native differences
//!
//! - The source's globals are mutable objects shared by the whole program;
//!   here they are `Mutex<Colorizer>` statics, locked for each use.
//! - The source's file-local `InitConsole` object switches Windows consoles to
//!   virtual-terminal processing at start-up and writes the reset-all code to
//!   `std::cout` and `std::cerr` at exit. Rust runs no static constructors or
//!   destructors. The reset is performed when a TOPP tool runs as its own
//!   executable (`cli::run` writes `\x1b[0m` to every terminal among its two
//!   streams afterwards); a library caller that leaves a color on should write
//!   [`COLOR_UNDO_ALL`](crate::concept::colorizer::COLOR_UNDO_ALL) itself. The Windows mode switch is not performed: a
//!   Windows terminal must already accept ANSI sequences (Windows Terminal and
//!   current consoles do).
//! - `isTTY` recognises the standard streams by object identity; here the
//!   caller names the stream with [`ConsoleStream`](crate::concept::colorizer::ConsoleStream).

use std::fmt;
use std::io::{self, IsTerminal, Write};
use std::sync::Mutex;

/// Text colors and styles a [`Colorizer`] supports: source
/// `enum class ConsoleColor`.
///
/// The values are, in source order, `RED`, `GREEN`, `YELLOW`, `BLUE`,
/// `MAGENTA`, `CYAN`, `UNDERLINE`, `BRIGHT` (keeps the foreground color but
/// makes it brighter) and `INVERT` (swaps foreground and background;
/// inverting twice stays inverted). The logging module colors its levels with
/// the same enum.
pub use super::log_stream::LogColor as ConsoleColor;

/// ANSI code that resets every color and style to the default: source
/// `color_undo_all_`.
pub const COLOR_UNDO_ALL: &str = "\x1b[0m";

/// The stream a [`Colorizer`] is written to, as far as coloring is concerned.
///
/// The source compares the stream's address with `std::cout` and `std::cerr`
/// and checks `isatty` for those two; every other stream is colored.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ConsoleStream {
    /// The process's standard output, colored only when it is a terminal.
    Stdout,
    /// The process's standard error, colored only when it is a terminal.
    Stderr,
    /// Any other writer, always colored ("mostly for testing purposes").
    Other,
}

/// Colors and styles console text: source `class Colorizer`.
///
/// # Examples
///
/// ```
/// use openms::concept::colorizer::{Colorizer, ConsoleColor, ConsoleStream};
///
/// let mut red = Colorizer::new(ConsoleColor::Red);
/// let mut out = Vec::new();
/// red.apply_to("make this red!").write_to(&mut out, ConsoleStream::Other).ok();
/// assert_eq!(out, b"\x1b[91mmake this red!\x1b[39m");
///
/// // Into a string, as to a stream that is not the console.
/// assert_eq!(red.apply().to_string(), "\x1b[91m");
/// assert_eq!(red.undo().to_string(), "\x1b[39m");
/// ```
#[derive(Debug)]
pub struct Colorizer {
    /// Source `color_`, fixed at construction.
    color: ConsoleColor,
    /// Source `undo_`: clear this colorizer's style after writing.
    undo: bool,
    /// Source `undo_all_`: clear all styles after writing.
    undo_all: bool,
    /// Source `undos_only`: write only the undo code, neither color nor text.
    undos_only: bool,
    /// Source `input_`: the text to write in color.
    input: String,
}

impl Colorizer {
    /// A colorizer for `color`: source `Colorizer(const ConsoleColor)`.
    ///
    /// Both undo flags start set and the text empty, as in the source.
    pub const fn new(color: ConsoleColor) -> Self {
        Self {
            color,
            undo: true,
            undo_all: true,
            undos_only: false,
            input: String::new(),
        }
    }

    /// The color or style: source member `color_`.
    pub const fn color(&self) -> ConsoleColor {
        self.color
    }

    /// Prepares to color the next stream written to until reset downstream:
    /// source `operator()()`, as in `red()`.
    ///
    /// Clears the text and all undo flags.
    pub fn apply(&mut self) -> &mut Self {
        self.input.clear();
        self.undo = false;
        self.undo_all = false;
        self.undos_only = false;
        self
    }

    /// Prepares to write `item` in color and reset the color right after:
    /// source `operator()(T s)`, as in `red("make this red!")`.
    ///
    /// `item` is converted to text as the source's `std::stringstream << s`
    /// does, here through `Display`.
    pub fn apply_to<T: fmt::Display>(&mut self, item: T) -> &mut Self {
        self.input.clear();
        // Writing to a String cannot fail; a failing Display impl leaves
        // whatever it wrote so far, as a failed stream insertion does.
        let _ = fmt::Write::write_fmt(&mut self.input, format_args!("{item}"));
        self.undo = true;
        self.undo_all = false;
        self.undos_only = false;
        self
    }

    /// Prepares to undo this colorizer's effect on the next stream written
    /// to: source `undo()`, as in `red.undo()`.
    ///
    /// Clears the text and sets the undo flag, but leaves the undo-all flag
    /// as it is, as the source does.
    pub fn undo(&mut self) -> &mut Self {
        self.input.clear();
        self.undo = true;
        self.undos_only = true;
        self
    }

    /// Prepares to reset the console to its default colors and styles: source
    /// `undoAll()`.
    pub fn undo_all(&mut self) -> &mut Self {
        self.input.clear();
        self.undo_all = true;
        self.undos_only = true;
        self
    }

    /// Whether `stream` is shown on a terminal rather than redirected to a
    /// file or the null device: source `isTTY`. Only standard output and
    /// standard error can be terminals; [`ConsoleStream::Other`] is `false`,
    /// as any other stream is in the source.
    pub fn is_tty(stream: ConsoleStream) -> bool {
        match stream {
            ConsoleStream::Stdout => io::stdout().is_terminal(),
            ConsoleStream::Stderr => io::stderr().is_terminal(),
            ConsoleStream::Other => false,
        }
    }

    /// ANSI code that switches `color` on: an entry of the source's `colors_`.
    pub fn enable_code(color: ConsoleColor) -> &'static str {
        std::str::from_utf8(color.enable()).unwrap_or_default()
    }

    /// ANSI code that switches `color` off again: an entry of the source's
    /// `colors_`. The six colors share the default-foreground code
    /// `\x1b[39m`, which is why undoing does not nest.
    pub fn disable_code(color: ConsoleColor) -> &'static str {
        std::str::from_utf8(color.disable()).unwrap_or_default()
    }

    /// Writes the colorizer to `out`: source `operator<<(std::ostream&,
    /// Colorizer&)`.
    ///
    /// Unless only an undo was prepared, this writes the color's enable code
    /// and the text; then the reset-all code if the undo-all flag is set,
    /// otherwise the color's disable code if the undo flag is set. For
    /// [`ConsoleStream::Stdout`] and [`ConsoleStream::Stderr`] the codes are
    /// written only when that stream is a terminal; the text always is.
    ///
    /// # Errors
    ///
    /// Any error of `out`.
    pub fn write_to<W: Write + ?Sized>(
        &self,
        out: &mut W,
        stream: ConsoleStream,
    ) -> io::Result<()> {
        let colored = match stream {
            ConsoleStream::Other => true,
            _ => Self::is_tty(stream),
        };
        self.output(colored, |text| out.write_all(text.as_bytes()))
    }

    /// The text to be written in color: the source's friend accessor
    /// `getInternalChars_`, which `IndentedStream` uses to re-indent it.
    pub fn internal_chars(&self) -> &str {
        &self.input
    }

    /// Replaces the text to be written in color: the source's friend accessor
    /// `setInternalChars_`.
    pub fn set_internal_chars(&mut self, data: &str) {
        self.input.clear();
        self.input.push_str(data);
    }

    /// Source `outputToStream_`, with `colorStream_`'s redirection check
    /// reduced to `colored`.
    fn output<E>(
        &self,
        colored: bool,
        mut write: impl FnMut(&str) -> Result<(), E>,
    ) -> Result<(), E> {
        if !self.undos_only {
            if colored {
                write(Self::enable_code(self.color))?;
            }
            write(&self.input)?;
        }
        if colored {
            if self.undo_all {
                write(COLOR_UNDO_ALL)?;
            } else if self.undo {
                write(Self::disable_code(self.color))?;
            }
        }
        Ok(())
    }
}

impl fmt::Display for Colorizer {
    /// Writes the colorizer as to a stream that is neither standard output
    /// nor standard error, so the codes are always included.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.output(true, |text| formatter.write_str(text))
    }
}

/// The global red colorizer: source `extern Colorizer red`.
pub static RED: Mutex<Colorizer> = Mutex::new(Colorizer::new(ConsoleColor::Red));
/// The global green colorizer: source `extern Colorizer green`.
pub static GREEN: Mutex<Colorizer> = Mutex::new(Colorizer::new(ConsoleColor::Green));
/// The global yellow colorizer: source `extern Colorizer yellow`.
pub static YELLOW: Mutex<Colorizer> = Mutex::new(Colorizer::new(ConsoleColor::Yellow));
/// The global blue colorizer: source `extern Colorizer blue`.
pub static BLUE: Mutex<Colorizer> = Mutex::new(Colorizer::new(ConsoleColor::Blue));
/// The global magenta colorizer: source `extern Colorizer magenta`.
pub static MAGENTA: Mutex<Colorizer> = Mutex::new(Colorizer::new(ConsoleColor::Magenta));
/// The global cyan colorizer: source `extern Colorizer cyan`.
pub static CYAN: Mutex<Colorizer> = Mutex::new(Colorizer::new(ConsoleColor::Cyan));
/// The global brightening colorizer: source `extern Colorizer bright`.
pub static BRIGHT: Mutex<Colorizer> = Mutex::new(Colorizer::new(ConsoleColor::Bright));
/// The global underlining colorizer: source `extern Colorizer underline`.
pub static UNDERLINE: Mutex<Colorizer> = Mutex::new(Colorizer::new(ConsoleColor::Underline));
/// The global inverting colorizer: source `extern Colorizer invert`.
pub static INVERT: Mutex<Colorizer> = Mutex::new(Colorizer::new(ConsoleColor::Invert));
