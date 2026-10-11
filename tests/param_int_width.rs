// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! A 64-bit integer parameter against the C++ Release build: the B8c lane's
//! oracle (`../oracle/B8c-numbers`, the Linux x86_64 Release install at core
//! `bc9cc12`, cli `c19e494`, topp `174b576` on ibminode06, gcc 14.4.0, two
//! identical runs; retained in `tests/data/b8c_release`, provenance
//! `tests/data/b8c_release_provenance.json`).
//!
//! `ParamValue` holds 64 bits, and the source narrows to `int` only where it
//! reads the value as one (`ParamValue::operator int`, `ParamValue.cpp:445-451`:
//! the low 32 bits): `ParamEntry::isValid` (`Param.cpp:122`) and therefore
//! `Param::update` check that narrowed value against the restriction and keep
//! the 64-bit value, and `GaussTraceFitter::updateMembers_` keeps it whole.
//! `width.tsv` records `isValid` and `update`, `trace.tsv` the fitter.

use openms::analysis::feature_finder_picked::trace_fitter::TraceFitterParams;
use openms::param::{Param, ParamEntry, ParamUpdateOptions, ParamValue};

const WIDTH: &str = include_str!("data/b8c_release/width.tsv");
const TRACE: &str = include_str!("data/b8c_release/trace.tsv");

fn records<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|fields| fields[0] == kind)
        .collect()
}

/// `isValid\t<case>\t<value>\t<min_int>\t<max_int>\t<valid>\t<message>`: the
/// entry's verdict for every case, and the value the source's message names,
/// which is the narrowed one (`-1` for `i64::MAX`).
#[test]
fn is_valid_checks_the_low_32_bits_as_the_release_build() {
    let cases = records(WIDTH, "isValid");
    assert_eq!(cases.len(), 8);
    for fields in cases {
        let (case, value, min_int, max_int, valid, message) = (
            fields[1],
            fields[2].parse::<i64>().unwrap(),
            fields[3].parse::<i32>().unwrap(),
            fields[4].parse::<i32>().unwrap(),
            fields[5] == "1",
            fields.get(6).copied().unwrap_or(""),
        );
        let mut entry = ParamEntry::new("i", ParamValue::Integer(value), "", &[]).unwrap();
        entry.min_int = min_int;
        entry.max_int = max_int;
        assert_eq!(entry.is_valid().unwrap(), valid, "{case}");
        if !valid {
            let narrowed = message
                .strip_prefix("Invalid integer parameter value '")
                .and_then(|rest| rest.split('\'').next())
                .unwrap();
            assert_eq!(narrowed, (value as i32).to_string(), "{case}");
        }
    }
}

/// `update\t<case>\t<success>\t<stored>\t<log>`: the current value 9 with a
/// maximum of 0, updated from `old` leniently or strictly
/// (`update(old, verbose, add_unknown = false, fail_on_invalid, false)`).
#[test]
fn update_keeps_the_64_bit_value_the_narrowed_check_admits() {
    let cases = records(WIDTH, "update");
    assert_eq!(cases.len(), 4);
    for fields in cases {
        let (case, success, stored) = (
            fields[1],
            fields[2] == "1",
            fields[3].parse::<i64>().unwrap(),
        );
        let (old_value, fail_on_invalid_values) = match case {
            "equal_invalid" => (9, true),
            "i64_max_lenient" => (i64::MAX, false),
            "2_32_plus_1_lenient" => ((1 << 32) + 1, false),
            "2_32_plus_1_strict" => ((1 << 32) + 1, true),
            other => panic!("{other}"),
        };
        let mut current = Param::new();
        current
            .set_value("same", ParamValue::Integer(9), "", &[])
            .unwrap();
        current.set_max_int("same", 0).unwrap();
        let mut old = Param::new();
        old.set_value("same", ParamValue::Integer(old_value), "", &[])
            .unwrap();
        let report = current
            .update_with_options(
                &old,
                ParamUpdateOptions {
                    fail_on_invalid_values,
                    ..ParamUpdateOptions::default()
                },
            )
            .unwrap();
        assert_eq!(report.success, success, "{case}");
        assert_eq!(
            current.value("same").unwrap(),
            &ParamValue::Integer(stored),
            "{case}"
        );
    }
}

/// `case\t<value>\taccepted\tparam\t<stored>\tmember\t<member>`:
/// `GaussTraceFitter::setParameters` takes every `max_iteration`, beyond
/// `int` too, and keeps it whole in `param_` and in `max_iterations_`.
#[test]
fn trace_fitter_takes_max_iteration_beyond_int_as_the_release_build() {
    let cases = records(TRACE, "case");
    assert_eq!(cases.len(), 6);
    for fields in cases {
        assert_eq!(fields[2], "accepted", "{fields:?}");
        let value = fields[1].parse::<i64>().unwrap();
        let (stored, member) = (
            fields[4].parse::<i64>().unwrap(),
            fields[6].parse::<i64>().unwrap(),
        );
        let mut param = TraceFitterParams::defaults().unwrap();
        param
            .set_value("max_iteration", ParamValue::Integer(value), "", &[])
            .unwrap();
        let (fitter, warnings) = TraceFitterParams::from_param(&param).unwrap();
        assert!(warnings.is_empty(), "{value}: {warnings:?}");
        assert_eq!(fitter.max_iteration, member, "{value}");
        assert_eq!(
            fitter.to_param().unwrap().value("max_iteration").unwrap(),
            &ParamValue::Integer(stored),
            "{value}"
        );
    }
}
