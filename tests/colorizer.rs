// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `Colorizer_test.cpp` (core `bc9cc12`), one test per testable
//! `START_SECTION`, literals transcribed (tier 3), plus the state rules read
//! from `Colorizer.cpp` (tier 3, source review): the initial undo flags, `undo`
//! leaving the undo-all flag, and the redirection check on the standard
//! streams. The `[EXTRA] visual inspection` section writes to the console and
//! asserts nothing; it has no counterpart.

use openms::concept::colorizer::{
    BLUE, BRIGHT, COLOR_UNDO_ALL, CYAN, Colorizer, ConsoleColor, ConsoleStream, GREEN, INVERT,
    MAGENTA, RED, UNDERLINE, YELLOW,
};
use std::fmt::Write as _;

const TEST_INT: i32 = 15;
const RED_ANSI: &str = "\x1b[91m";
const GREEN_ANSI: &str = "\x1b[92m";
const YELLOW_ANSI: &str = "\x1b[93m";
const BLUE_ANSI: &str = "\x1b[94m";
const MAGENTA_ANSI: &str = "\x1b[95m";
const CYAN_ANSI: &str = "\x1b[96m";
const RESET_COLOR_ANSI: &str = "\x1b[39m";

#[test]
fn class_test_constructor() {
    let mut test = Colorizer::new(ConsoleColor::Blue);
    let mut s = String::from("-");
    write!(s, "{}", test.apply_to("test")).ok();
    s.push('-');
    assert_eq!(s, format!("-{BLUE_ANSI}test{RESET_COLOR_ANSI}-"));
}

#[test]
fn class_test_bracket_operator() {
    let mut test = Colorizer::new(ConsoleColor::Blue);
    let mut s = String::from("-");
    write!(s, "{}test-{TEST_INT}", test.apply()).ok();
    assert_eq!(s, format!("-{BLUE_ANSI}test-{TEST_INT}"));
}

#[test]
fn class_test_bracket_operator_with_an_item() {
    let mut test = Colorizer::new(ConsoleColor::Magenta);
    let mut s = String::from("-");
    write!(s, "{}test-{TEST_INT}", test.apply_to(TEST_INT)).ok();
    assert_eq!(
        s,
        format!("-{MAGENTA_ANSI}{TEST_INT}{RESET_COLOR_ANSI}test-{TEST_INT}")
    );
}

#[test]
fn class_test_undo() {
    let mut test = Colorizer::new(ConsoleColor::Cyan);
    let mut s = String::from("-");
    write!(s, "{}", test.apply()).ok();
    write!(s, "test{TEST_INT}").ok();
    write!(s, "{}", test.undo()).ok();
    s.push_str("nocol");
    assert_eq!(
        s,
        format!("-{CYAN_ANSI}test{TEST_INT}{RESET_COLOR_ANSI}nocol")
    );

    // double coloring + reset (using any Colorizer)
    let mut yellow_test = Colorizer::new(ConsoleColor::Yellow);
    let mut s = String::from("-");
    write!(s, "{}test", test.apply()).ok();
    write!(s, "{}{TEST_INT}", yellow_test.apply()).ok();
    write!(s, "{}nocol", yellow_test.undo()).ok();
    assert_eq!(
        s,
        format!("-{CYAN_ANSI}test{YELLOW_ANSI}{TEST_INT}{RESET_COLOR_ANSI}nocol")
    );
}

#[test]
fn the_color_table_matches_the_header() {
    let expected = [
        (ConsoleColor::Red, RED_ANSI, RESET_COLOR_ANSI),
        (ConsoleColor::Green, GREEN_ANSI, RESET_COLOR_ANSI),
        (ConsoleColor::Yellow, YELLOW_ANSI, RESET_COLOR_ANSI),
        (ConsoleColor::Blue, BLUE_ANSI, RESET_COLOR_ANSI),
        (ConsoleColor::Magenta, MAGENTA_ANSI, RESET_COLOR_ANSI),
        (ConsoleColor::Cyan, CYAN_ANSI, RESET_COLOR_ANSI),
        (ConsoleColor::Underline, "\x1b[4m", "\x1b[24m"),
        (ConsoleColor::Bright, "\x1b[1m", "\x1b[22m"),
        (ConsoleColor::Invert, "\x1b[7m", "\x1b[27m"),
    ];
    for (color, enable, disable) in expected {
        assert_eq!(Colorizer::enable_code(color), enable, "{color:?}");
        assert_eq!(Colorizer::disable_code(color), disable, "{color:?}");
    }
    assert_eq!(COLOR_UNDO_ALL, "\x1b[0m");
}

#[test]
fn a_new_colorizer_carries_both_undo_flags() {
    // Member initialisers undo_ = true, undo_all_ = true: written before any
    // operator(), a colorizer emits its color, no text and the reset-all code.
    let fresh = Colorizer::new(ConsoleColor::Green);
    assert_eq!(fresh.to_string(), format!("{GREEN_ANSI}{COLOR_UNDO_ALL}"));
    // undo() leaves undo_all_ set, so a never-applied colorizer resets all.
    let mut fresh = Colorizer::new(ConsoleColor::Green);
    assert_eq!(fresh.undo().to_string(), COLOR_UNDO_ALL);
}

#[test]
fn undo_after_undo_all_still_resets_everything() {
    let mut red = Colorizer::new(ConsoleColor::Red);
    assert_eq!(red.undo_all().to_string(), COLOR_UNDO_ALL);
    assert_eq!(red.undo().to_string(), COLOR_UNDO_ALL);
    // operator() clears the flag again.
    assert_eq!(red.apply().to_string(), RED_ANSI);
    assert_eq!(red.undo().to_string(), RESET_COLOR_ANSI);
}

#[test]
fn writing_does_not_change_the_colorizer() {
    let mut blue = Colorizer::new(ConsoleColor::Blue);
    blue.apply_to("x");
    let once = blue.to_string();
    assert_eq!(blue.to_string(), once);
    assert_eq!(blue.internal_chars(), "x");
}

#[test]
fn internal_chars_can_be_replaced() {
    let mut underline = Colorizer::new(ConsoleColor::Underline);
    underline.apply_to("first");
    underline.set_internal_chars("second\n  line");
    assert_eq!(underline.internal_chars(), "second\n  line");
    assert_eq!(underline.to_string(), "\x1b[4msecond\n  line\x1b[24m");
}

#[test]
fn write_to_colors_other_streams_and_checks_the_console() {
    let mut red = Colorizer::new(ConsoleColor::Red);
    red.apply_to("text");
    let mut out = Vec::new();
    assert!(red.write_to(&mut out, ConsoleStream::Other).is_ok());
    assert_eq!(
        out,
        format!("{RED_ANSI}text{RESET_COLOR_ANSI}").into_bytes()
    );
    assert!(!Colorizer::is_tty(ConsoleStream::Other));

    // Written "as" standard output, the codes depend on whether it is a
    // terminal; the text is always written.
    let mut out = Vec::new();
    assert!(red.write_to(&mut out, ConsoleStream::Stdout).is_ok());
    let expected = if Colorizer::is_tty(ConsoleStream::Stdout) {
        format!("{RED_ANSI}text{RESET_COLOR_ANSI}")
    } else {
        "text".to_owned()
    };
    assert_eq!(out, expected.into_bytes());

    // An undo alone writes no text.
    let mut out = Vec::new();
    assert!(red.undo().write_to(&mut out, ConsoleStream::Other).is_ok());
    assert_eq!(out, RESET_COLOR_ANSI.as_bytes());
}

#[test]
fn global_colorizers_have_their_colors_and_keep_state() {
    let globals = [
        (&RED, ConsoleColor::Red),
        (&GREEN, ConsoleColor::Green),
        (&YELLOW, ConsoleColor::Yellow),
        (&BLUE, ConsoleColor::Blue),
        (&MAGENTA, ConsoleColor::Magenta),
        (&CYAN, ConsoleColor::Cyan),
        (&BRIGHT, ConsoleColor::Bright),
        (&UNDERLINE, ConsoleColor::Underline),
        (&INVERT, ConsoleColor::Invert),
    ];
    for (global, color) in globals {
        let colorizer = global.lock().unwrap_or_else(|poison| poison.into_inner());
        assert_eq!(colorizer.color(), color);
    }
    let mut cyan = CYAN.lock().unwrap_or_else(|poison| poison.into_inner());
    cyan.apply();
    assert_eq!(cyan.undo().to_string(), RESET_COLOR_ANSI);
    cyan.undo_all();
    assert_eq!(cyan.undo().to_string(), COLOR_UNDO_ALL);
}
