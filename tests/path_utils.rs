// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `SYSTEM/PathUtils.h`: `PathUtils::basename` and `to_path`.
//!
//! The class test has a single section, `to_path`, whose UTF-8 round-trip
//! literals are reproduced here byte for byte. `PathUtils::basename` has no
//! section at all; its cases are derived from the source expression
//! `file.substr(file.find_last_of("\\/") + 1)`.

use openms::system::path_utils::{MAX_PATH_BYTES, basename, to_path};

/// Class-test section `std::filesystem::path to_path(const std::string& s)`.
///
/// The source test asserts a byte-exact round trip through `u8string()` for
/// every one of these inputs; here the round trip is through `to_str`, which is
/// the same bytes on Unix and the same scalar values on Windows.
#[test]
fn to_path_round_trips_every_class_test_literal() {
    // "ASCII_only.txt" — identical on all platforms.
    assert_eq!(
        to_path("ASCII_only.txt").unwrap().to_str(),
        Some("ASCII_only.txt")
    );

    // U+00E4 as UTF-8, 0xC3 0xA4.
    let ae = std::str::from_utf8(&[0xC3, 0xA4]).unwrap().to_string() + ".mzML";
    assert_eq!(to_path(&ae).unwrap().to_str(), Some(ae.as_str()));

    // CJK: Japanese 日本語.
    let cjk = std::str::from_utf8(&[0xE6, 0x97, 0xA5, 0xE6, 0x9C, 0xAC, 0xE8, 0xAA, 0x9E])
        .unwrap()
        .to_string()
        + ".mzML";
    assert_eq!(cjk, "日本語.mzML");
    assert_eq!(to_path(&cjk).unwrap().to_str(), Some(cjk.as_str()));

    // Embedded spaces and parentheses.
    let spaced = "my data file (run 1).mzML";
    assert_eq!(to_path(spaced).unwrap().to_str(), Some(spaced));

    // Well beyond the historical Windows MAX_PATH of 260: no truncation.
    let long_name = "a".repeat(300) + ".mzML";
    let long_path = to_path(&long_name).unwrap();
    assert_eq!(long_path.to_str(), Some(long_name.as_str()));
    assert_eq!(long_path.to_str().unwrap().len(), 305);

    // Mixed CJK and Latin-1 diacritics with a space: 测试 äö.
    let mixed = std::str::from_utf8(&[
        0xE6, 0xB5, 0x8B, 0xE8, 0xAF, 0x95, b' ', 0xC3, 0xA4, 0xC3, 0xB6,
    ])
    .unwrap()
    .to_string()
        + ".featureXML";
    assert_eq!(mixed, "测试 äö.featureXML");
    assert_eq!(to_path(&mixed).unwrap().to_str(), Some(mixed.as_str()));
}

/// The source's Windows-only ANSI fallback has no counterpart.
///
/// The class test wraps the `_WIN32` cases in a `try`/`catch` because the lone
/// byte 0xE4 is not valid UTF-8 and the unguarded source threw
/// `std::system_error` on it. A Rust `&str` cannot hold that byte sequence at
/// all, so the failure mode the fallback exists for is unrepresentable and the
/// only remaining refusals are the port's own bounds.
#[test]
fn to_path_refuses_only_what_no_filesystem_could_use() {
    assert!(to_path("\\\\server\\share\\data.mzML").is_ok());
    assert!(to_path("").is_ok());
    assert!(to_path("a\0b").is_err());
    assert!(to_path(&"a".repeat(MAX_PATH_BYTES + 1)).is_err());
}

/// `PathUtils::basename`, which the class test does not exercise.
///
/// A path without a separator yields itself (the source reaches this through
/// `npos + 1` wrapping to zero), and a path ending in a separator — including a
/// path that is only a separator — yields the empty string.
#[test]
fn basename_takes_everything_after_the_last_separator() {
    assert_eq!(basename("/data/run1.mzML"), "run1.mzML");
    assert_eq!(basename("C:\\data\\run1.mzML"), "run1.mzML");
    assert_eq!(basename("mixed/separators\\run1.mzML"), "run1.mzML");
    assert_eq!(basename("run1.mzML"), "run1.mzML");
    assert_eq!(basename("/data/"), "");
    assert_eq!(basename("C:\\data\\"), "");
    assert_eq!(basename("/"), "");
    assert_eq!(basename("\\"), "");
    assert_eq!(basename(""), "");
    // Purely lexical: no normalisation of dot segments or repeated separators.
    assert_eq!(basename("/data/./run1.mzML"), "run1.mzML");
    assert_eq!(basename("/data//"), "");
    assert_eq!(basename("/data/.."), "..");
}

/// The lexical helpers never touch the filesystem.
#[test]
fn basename_and_to_path_ignore_whether_the_path_exists() {
    let missing = "/nonexistent-directory-for-openms-tests/run1.mzML";
    assert_eq!(basename(missing), "run1.mzML");
    assert_eq!(to_path(missing).unwrap().to_str(), Some(missing));
}
