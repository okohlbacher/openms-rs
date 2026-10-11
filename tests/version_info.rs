// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `VersionInfo_test.cpp` (core `bc9cc12`), one test per testable
//! `START_SECTION`, literals transcribed (tier 3). `getTime` is
//! `NOT_TESTABLE` upstream and not ported. The update-check use of the parser
//! is covered by `tests/update_check.rs` (feature `network`).

use openms::concept::version_info::{VersionDetails, get_version, get_version_struct};

fn version(major: i32, minor: i32, patch: i32) -> VersionDetails {
    VersionDetails {
        major,
        minor,
        patch,
        ..VersionDetails::default()
    }
}

#[test]
fn get_version_is_the_trimmed_package_version() {
    assert_eq!(get_version(), env!("CARGO_PKG_VERSION").trim());
}

#[test]
fn get_version_struct_parses_the_package_version() {
    // Parsed independently of VersionDetails::create, as the source test does.
    let mut parts = env!("CARGO_PKG_VERSION").split(['.', '-']);
    let mut next = || parts.next().and_then(|part| part.parse::<i32>().ok());
    let (major, minor, patch) = (next(), next(), next());
    let parsed = get_version_struct();
    assert_eq!(Some(parsed.major), major);
    assert_eq!(Some(parsed.minor), minor);
    assert_eq!(Some(parsed.patch), patch);
}

#[test]
fn default_constructor_is_empty() {
    let detail = VersionDetails::default();
    assert!(detail == VersionDetails::EMPTY);
    assert!(detail.is_empty());
}

#[test]
fn copy_constructor() {
    let detail = VersionDetails::create("1.9.2");
    let c = detail.clone();
    assert_eq!(c.major, detail.major);
    assert_eq!(c.minor, detail.minor);
    assert_eq!(c.patch, detail.patch);
}

#[test]
fn less_than() {
    let detail = VersionDetails::create("1.9.2");
    let mut c = version(1, 9, 2);
    assert!(!detail.is_less_than(&c));
    c.patch = 3;
    assert!(detail.is_less_than(&c));
    c.patch = 1;
    assert!(!detail.is_less_than(&c));
    c.major = 2;
    assert!(detail.is_less_than(&c));
}

#[test]
fn equality() {
    let detail = VersionDetails::create("1.9.2");
    let mut c = version(1, 9, 2);
    assert!(detail == c);
    c.patch = 3;
    assert!(detail != c);
    c.patch = 1;
    assert!(detail != c);
    c.major = 2;
    assert!(detail != c);
}

#[test]
fn greater_than() {
    let detail = VersionDetails::create("1.9.2");
    let mut c = version(1, 9, 2);
    assert!(!detail.is_greater_than(&c));
    c.patch = 3;
    assert!(!detail.is_greater_than(&c));
    c.patch = 1;
    assert!(detail.is_greater_than(&c));
    c.patch = 11;
    assert!(detail.is_less_than(&c));
    c.patch = 2;
    assert!(!detail.is_greater_than(&c));
    // Any version with a pre-release identifier is "less than" the release.
    c.pre_release = "alpha".into();
    assert!(detail.is_greater_than(&c));
}

#[test]
fn create() {
    assert!(VersionDetails::create("1.9.2") == version(1, 9, 2));
    assert!(VersionDetails::create("1.9") == version(1, 9, 0));
    assert!(VersionDetails::create("1.0") == version(1, 0, 0));
    assert!(VersionDetails::create("somestring") == version(0, 0, 0));
    assert!(VersionDetails::create("1.2a.bla") == version(0, 0, 0));
    let detail = VersionDetails::create("1.2.1-bla");
    let mut c = version(1, 2, 1);
    c.pre_release = "bla".into();
    assert_eq!(detail.major, c.major);
    assert_eq!(detail.minor, c.minor);
    assert_eq!(detail.patch, c.patch);
    assert_eq!(detail.pre_release, c.pre_release);
    assert!(detail == c);
}

#[test]
fn create_cannot_tell_a_failure_from_zero_but_parse_can() {
    assert!(VersionDetails::create("0.0.0").is_empty());
    assert!(VersionDetails::create("x").is_empty());
    assert_eq!(VersionDetails::parse("0.0.0"), Some(VersionDetails::EMPTY));
    assert_eq!(VersionDetails::parse("x"), None);
}
