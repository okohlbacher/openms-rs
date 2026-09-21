// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! `APPLICATIONS/ParameterInformation.h`: `ParameterInformation_test.cpp`
//! (cli `c19e494`), transcribed with its literals (tier 3).
//!
//! The source's restriction members hold `±numeric_limits::max()` when unset;
//! this port holds `None`. The class test's sentinel comparisons are
//! transcribed as the value the port's consumers read, `unwrap_or` the
//! sentinel, and as `None`.
#![cfg(feature = "paramxml")]

use openms::cli::{ParameterInformation, ParameterType};
use openms::param::ParamValue;

fn assert_unrestricted(p: &ParameterInformation) {
    assert!(p.valid_strings.is_empty());
    assert_eq!(p.min_int, None);
    assert_eq!(p.max_int, None);
    assert_eq!(p.min_float, None);
    assert_eq!(p.max_float, None);
    assert_eq!(p.min_int.unwrap_or(-i32::MAX), -i32::MAX);
    assert_eq!(p.max_int.unwrap_or(i32::MAX), i32::MAX);
    assert_eq!(p.min_float.unwrap_or(-f64::MAX), -f64::MAX);
    assert_eq!(p.max_float.unwrap_or(f64::MAX), f64::MAX);
}

/// `START_SECTION(ParameterInformation())` (`ParameterInformation_test.cpp:28-46`).
#[test]
fn upstream_default_constructor() {
    let p = ParameterInformation::default();
    assert_eq!(p.name, "");
    assert_eq!(p.kind, ParameterType::None);
    assert_eq!(p.default_value, ParamValue::Empty);
    assert_eq!(p.description, "");
    assert_eq!(p.argument, "");
    assert!(p.required);
    assert!(!p.advanced);
    assert!(p.tags.is_empty());
    assert_unrestricted(&p);
}

fn pi1() -> ParameterInformation {
    ParameterInformation::new(
        "pi1_name",
        ParameterType::String,
        "<STRING>",
        ParamValue::String("def_value".into()),
        "this is a description",
        false,
        true,
    )
    .with_tags(&["tag1", "tag2"])
}

fn assert_pi1(p: &ParameterInformation) {
    assert_eq!(p.name, "pi1_name");
    assert_eq!(p.kind, ParameterType::String);
    assert_eq!(p.default_value, ParamValue::String("def_value".into()));
    assert_eq!(p.description, "this is a description");
    assert_eq!(p.argument, "<STRING>");
    assert!(!p.required);
    assert!(p.advanced);
    assert_eq!(p.tags.len(), 2);
    assert_eq!(p.tags[0], "tag1");
    assert_eq!(p.tags[1], "tag2");
}

/// `START_SECTION((ParameterInformation(n, t, arg, def, desc, req, adv,
/// tag_values)))` (`ParameterInformation_test.cpp:54-77`).
#[test]
fn upstream_detailed_constructor() {
    let p = pi1();
    assert_pi1(&p);
    assert_unrestricted(&p);
}

/// `START_SECTION((ParameterInformation& operator=(...)))`
/// (`ParameterInformation_test.cpp:79-117`): assignment is `Clone`.
#[test]
fn upstream_assignment() {
    let p = pi1();
    assert_pi1(&p);
    assert_unrestricted(&p);
    let mut p2 = ParameterInformation::default();
    assert!(p2.required);
    p2 = p.clone();
    assert_pi1(&p2);
}
