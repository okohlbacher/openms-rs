// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The closing line `TOPPBase::main` prints on standard output once the tool
//! body has returned (`TOPPBase.cpp:413-424`):
//! `<tool> took <t> (wall), <t> (CPU), <t> (system), <t> (user); Peak Memory
//! Usage: <n> MB.`
//!
//! Its figures differ from run to run, so a test that compares a tool's
//! standard output takes the line off with [`split_took_line`], which checks
//! its shape, and compares the rest.

/// Split `out` into the output before the closing line and the closing line,
/// which is checked: the tool's name, four times in `StopWatch::toString`'s
/// seconds form (or `n/a` where this port's `StopWatch` cannot read a
/// component on the platform), and the peak-memory part where the platform
/// reports it. Returns `out` unchanged and `None` when its last line is not a
/// closing line.
#[allow(dead_code)]
pub fn split_took_line(tool: &str, out: &str) -> (String, Option<String>) {
    let body = out.strip_suffix('\n').unwrap_or(out);
    let (before, last) = match body.rsplit_once('\n') {
        Some((before, last)) => (format!("{before}\n"), last),
        None => (String::new(), body),
    };
    let prefix = format!("{tool} took ");
    let Some(rest) = last.strip_prefix(&prefix) else {
        return (out.to_owned(), None);
    };
    let rest = rest.strip_suffix('.').unwrap_or_else(|| panic!("{last:?}"));
    let (times, memory) = match rest.split_once("; Peak Memory Usage: ") {
        Some((times, memory)) => (times, Some(memory)),
        None => (rest, None),
    };
    let parts: Vec<&str> = times.split(", ").collect();
    assert_eq!(parts.len(), 4, "{last:?}");
    for (part, label) in parts.iter().zip(["(wall)", "(CPU)", "(system)", "(user)"]) {
        if *part == format!("n/a {label}") {
            continue;
        }
        let seconds = part
            .strip_suffix(&format!(" s {label}"))
            .unwrap_or_else(|| panic!("{last:?}"));
        let (whole, fraction) = seconds
            .split_once('.')
            .unwrap_or_else(|| panic!("{last:?}"));
        assert!(
            !whole.is_empty() && whole.bytes().all(|b| b.is_ascii_digit()),
            "{last:?}"
        );
        assert!(
            fraction.len() == 2 && fraction.bytes().all(|b| b.is_ascii_digit()),
            "{last:?}"
        );
    }
    if let Some(memory) = memory {
        let megabytes = memory
            .strip_suffix(" MB")
            .unwrap_or_else(|| panic!("{last:?}"));
        assert!(megabytes.bytes().all(|b| b.is_ascii_digit()), "{last:?}");
    }
    (before, Some(last.to_owned()))
}

/// The output before the closing line, with the line checked as
/// [`split_took_line`] checks it.
#[allow(dead_code)]
pub fn strip_took_line(tool: &str, out: &str) -> String {
    split_took_line(tool, out).0
}
