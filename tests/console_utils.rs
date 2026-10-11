// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// `cli` lives behind `paramxml`; without this gate `cargo test
// --no-default-features` fails to compile.
#![cfg(feature = "paramxml")]

//! Port of `ConsoleUtils_test.cpp` (core `bc9cc12`), one test per testable
//! `START_SECTION`, literals transcribed (tier 3).
//!
//! The source test pins the singleton's width to 9 through its
//! `ConsoleWidthTest` friend; [`ConsoleUtils::with_console_width`] is that hook
//! and the instance methods break at the pinned width. The measured width and
//! the byte-exact line breaking are also compared with the executed Release
//! build in `tests/topp_cli_console.rs`.

use openms::cli::ConsoleUtils;

const TEST_CONSOLE_WIDTH: usize = 9;

fn console() -> ConsoleUtils {
    ConsoleUtils::with_console_width(TEST_CONSOLE_WIDTH as i32)
}

fn x(count: usize) -> String {
    "x".repeat(count)
}

fn spaces(count: usize) -> String {
    " ".repeat(count)
}

#[test]
fn get_console_width() {
    assert_eq!(console().get_console_width(), TEST_CONSOLE_WIDTH as i32);
}

#[test]
fn get_instance_measures_once() {
    let first = ConsoleUtils::get_instance();
    let second = ConsoleUtils::get_instance();
    assert!(std::ptr::eq(first, second));
    let width = first.get_console_width();
    assert!(width == i32::MAX || width >= 10, "{width}");
    // The static wrappers use the singleton's width.
    let text = "a short text";
    assert_eq!(
        ConsoleUtils::break_string(text, 2, 10, 0),
        first.break_text(text, 2, 10, 0)
    );
    assert_eq!(
        ConsoleUtils::break_string_list(text, 2, 10, 0),
        first.break_lines(text, 2, 10, 0)
    );
    assert_eq!(
        ConsoleUtils::break_bytes_list(text.as_bytes(), 2, 10, 0),
        first.break_byte_lines(text.as_bytes(), 2, 10, 0)
    );
}

#[test]
fn break_string_list() {
    let console = console();
    let x20 = x(TEST_CONSOLE_WIDTH * 2 + 1); // two full lines plus one 'x'
    let x_c = x(TEST_CONSOLE_WIDTH); // full console width of 'x'

    // indent = 0, twice
    for _ in 0..2 {
        assert_eq!(
            console.break_text(&x20, 0, 10, 0),
            format!("{x_c}\n{x_c}\nx")
        );
    }
    // indent = 3
    let indent = 3;
    let short_x = x(TEST_CONSOLE_WIDTH - indent);
    let s_indent = spaces(indent);
    assert_eq!(
        console.break_text(&x20, indent, 10, 0),
        format!("{x_c}\n{s_indent}{short_x}\n{s_indent}xxxx")
    );
    // prefilled first line
    let prefill = 5;
    let first_x = x(TEST_CONSOLE_WIDTH - prefill);
    assert_eq!(
        console.break_text(&x20, indent, 10, prefill),
        format!("{first_x}\n{s_indent}{short_x}\n{s_indent}{short_x}\n{s_indent}xxx")
    );
    // just a line break
    assert_eq!(console.break_text("\n", 0, 10, 0), "\n");
    // just a line break with indent
    assert_eq!(
        console.break_text("\n", 3, 10, 0),
        format!("\n{}", spaces(3))
    );
    // prefilled line break with indent (makes no difference)
    assert_eq!(
        console.break_text("\n", 3, 10, 5),
        format!("\n{}", spaces(3))
    );
    // text with a line break, indent and prefill: one byte left on line one
    assert_eq!(
        console.break_text("xxx\n", 3, 10, TEST_CONSOLE_WIDTH - 1),
        format!("x\n{0}xx\n{0}", spaces(3))
    );
    // corner case: one byte per line
    let indent = TEST_CONSOLE_WIDTH - 1;
    let pad = spaces(indent);
    assert_eq!(
        console.break_text("xxx\n", indent, 10, indent),
        format!("x\n{pad}x\n{pad}x\n{pad}\n{pad}")
    );
    // max_lines: two bytes per line, three lines kept
    let indent = TEST_CONSOLE_WIDTH - 2;
    let pad = spaces(indent);
    assert_eq!(
        console.break_text(&x(99), indent, 3, indent),
        format!("xx\n{pad}...\n{pad}x")
    );
}

#[test]
fn break_string_list_returns_one_element_per_line() {
    let console = console();
    let lines = console.break_lines(&x(19), 3, 10, 0);
    assert_eq!(lines, [x(9), format!("   {}", x(6)), "   xxxx".to_owned()]);
    let bytes = console.break_byte_lines(x(19).as_bytes(), 3, 10, 0);
    assert_eq!(bytes.len(), 3);
    assert!(console.break_lines("", 3, 10, 0).is_empty());
}

#[test]
fn a_split_multibyte_character_is_replaced_only_in_the_string_form() {
    // Nine bytes per line: "aaaaaaaa" plus the first byte of "é".
    let console = console();
    let text = "aaaaaaaa\u{e9}b";
    let bytes = console.break_byte_lines(text.as_bytes(), 0, 10, 0);
    assert_eq!(bytes[0], b"aaaaaaaa\xc3");
    assert_eq!(bytes[1], b"\xa9b");
    let lines = console.break_lines(text, 0, 10, 0);
    assert_eq!(lines, ["aaaaaaaa\u{fffd}", "\u{fffd}b"]);
}
