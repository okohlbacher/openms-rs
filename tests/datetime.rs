// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::data_structures::{
    DateTime,
    datetime::{DATETIME_FORMATS, MAX_DATETIME_INPUT_BYTES},
};
use std::hash::{Hash, Hasher};

fn unhex(value: &str) -> String {
    String::from_utf8(
        value
            .as_bytes()
            .chunks_exact(2)
            .map(|pair| u8::from_str_radix(std::str::from_utf8(pair).unwrap(), 16).unwrap())
            .collect(),
    )
    .unwrap()
}
fn seed(value: &str) -> DateTime {
    if value.is_empty() {
        return DateTime::default();
    }
    if let Some(time) = value.strip_prefix("time:") {
        let mut d = DateTime::default();
        d.set_time(time).unwrap();
        d
    } else {
        DateTime::parse(value).unwrap()
    }
}
#[test]
fn all_19_published_class_test_string_literals() {
    let mut count = 0;
    for row in include_str!("data/datetime_class_literals.tsv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = row.split('\t').collect();
        let value = DateTime::parse(f[1]);
        assert_eq!(value.is_ok(), f[2] == "ok", "source line {}", f[0]);
        if let Ok(value) = value {
            assert_eq!(value.get(), f[3], "source line {}", f[0]);
        }
        count += 1;
    }
    assert_eq!(count, 19);
}
#[test]
fn all_301_executed_cpp_probe_rows_and_seven_renderings() {
    let corrections: std::collections::BTreeMap<_, _> =
        include_str!("data/datetime_calendar_corrections.tsv")
            .lines()
            .skip(1)
            .map(|row| {
                let fields: Vec<_> = row.split('\t').collect();
                (fields[0], fields)
            })
            .collect();
    let mut corrected = 0;
    let mut count = 0;
    for row in include_str!("data/datetime_cpp_probe.tsv").lines().skip(1) {
        let f: Vec<_> = row.split('\t').collect();
        assert_eq!(f.len(), 19);
        let mut d = seed(&unhex(f[1]));
        let input = unhex(f[3]);
        let format = unhex(f[4]);
        let n: Vec<u32> = f[5]
            .split(',')
            .filter(|s| !s.is_empty())
            .filter_map(|s| s.parse().ok())
            .collect();
        let result = match f[2] {
            "set" => d.set(&input),
            "date" => d.set_date(&input),
            "time" => d.set_time(&input),
            "from" => DateTime::from_format(&input, &format).map(|value| d = value),
            "add" => d.add_seconds(f[5].parse().unwrap()).map(|_| ()),
            "clear" => {
                d.clear();
                Ok(())
            }
            "none" => Ok(()),
            "parts" => d.set_components(n[0], n[1], n[2], n[3], n[4], n[5]),
            "dateparts" => d.set_date_components(n[0], n[1], n[2]),
            "timeparts" => d.set_time_components(n[0], n[1], n[2]),
            _ => panic!("unknown operation"),
        };
        assert_eq!(result.is_ok(), f[6] == "ok", "{} {:?}", f[0], input);
        assert_eq!(d.is_valid(), f[7] == "1", "{}", f[0]);
        // Original source output is retained. For the 35 macOS timegm failures,
        // independent Python calendar month-stepping supplies corrected fields.
        let expected = if let Some(row) = corrections.get(f[0]) {
            corrected += 1;
            assert_eq!(f[2], "add");
            &row[1..]
        } else {
            &f[8..]
        };
        let c: Vec<i32> = expected[0].split(',').map(|n| n.parse().unwrap()).collect();
        assert_eq!(
            d.components(),
            (c[0], c[1], c[2], c[3], c[4], c[5]),
            "{} {:?}",
            f[0],
            input
        );
        for (index, format) in DATETIME_FORMATS.iter().enumerate() {
            assert_eq!(
                d.format(format).unwrap(),
                unhex(expected[index + 1]),
                "{} {:?} {}",
                f[0],
                input,
                format
            );
        }
        assert_eq!(d.get(), unhex(expected[8]), "{}", f[0]);
        assert_eq!(d.date_string(), unhex(expected[9]), "{}", f[0]);
        assert_eq!(d.time_string(), unhex(expected[10]), "{}", f[0]);
        count += 1;
    }
    assert_eq!(count, 301);
    assert_eq!(corrected, 35);
    assert_eq!(count - corrected, 266);
}
#[test]
fn numeric_getters_copy_and_source_default_lifecycle() {
    let mut d = DateTime::default();
    assert!(d.is_null());
    assert!(!d.is_valid());
    assert_eq!(d.components(), (0, 0, 0, 0, 0, 0));
    assert_eq!(d.millisecond(), 0);
    assert_eq!(d.iso_string(), "");
    assert_eq!(d.get(), "0000-00-00 00:00:00");
    d.set_components(12, 14, 2006, 11, 59, 58).unwrap();
    assert_eq!(d.date_components(), (12, 14, 2006));
    assert_eq!(d.time_components(), (11, 59, 58));
    assert_eq!(d.to_string(), "2006-12-14 11:59:58");
    assert_eq!(d.iso_string(), "2006-12-14T11:59:58");
    let copied = d;
    assert_eq!(copied, d);
    d.clear();
    assert_eq!(d, DateTime::default());
    assert_ne!(copied, d);
    d.set_components(5, 4, 666, 3, 2, 1).unwrap();
    assert_eq!(d.get(), "0666-05-04 03:02:01");
}
#[test]
fn independent_validity_partial_setters_and_millisecond_retention() {
    let mut d = DateTime::default();
    d.set_time("11:59:58").unwrap();
    assert!(d.is_valid());
    assert_eq!(d.get(), "0000-00-00 11:59:58");
    assert_eq!(d, DateTime::from_format("11:59:58", "hh:mm:ss").unwrap());
    d.set("2011-08-05T15:32:07.468").unwrap();
    d.set_date("02/03/2000").unwrap();
    d.set_time_components(1, 2, 3).unwrap();
    assert_eq!(
        d.format(DATETIME_FORMATS[1]).unwrap(),
        "2000-02-03T01:02:03.468"
    );
    d.set_components(2, 3, 2000, 1, 2, 3).unwrap();
    assert_eq!(d.millisecond(), 0);
    let null = DateTime::default();
    let mut partial = null;
    partial.set_time_components(0, 0, 0).unwrap();
    assert_ne!(partial, null);
    assert!(!partial.source_less(&null));
    assert!(!null.source_less(&partial));
}
#[test]
fn failed_operations_preserve_source_distinct_publication_policy() {
    let original = DateTime::parse("2011-08-05T15:32:07.468").unwrap();
    let mut d = original;
    assert!(d.set_date("2000-02-30").is_err());
    assert_eq!(d, original);
    assert!(d.set_time("23:59:60").is_err());
    assert_eq!(d, original);
    assert!(d.set_date_components(2, 30, 2000).is_err());
    assert_eq!(d, original);
    assert!(d.set_time_components(24, 0, 0).is_err());
    assert_eq!(d, original);
    assert!(d.set_components(1, 1, 0, 0, 0, 0).is_err());
    assert_eq!(d, original);
    let large = "x".repeat(MAX_DATETIME_INPUT_BYTES + 1);
    assert!(d.set(&large).is_err());
    assert_eq!(d, original);
    assert!(d.set_date(&large).is_err());
    assert_eq!(d, original);
    assert!(d.set_time(&large).is_err());
    assert_eq!(d, original);
    assert!(DateTime::from_format(&large, "unknown").is_err());
    assert!(d.set("invalid").is_err());
    assert_eq!(d, DateTime::default());
    assert_eq!(
        DateTime::from_format("invalid", DATETIME_FORMATS[0]).unwrap(),
        d
    );
}
#[test]
fn checked_integer_and_fraction_bounds_preserve_defined_domain() {
    for value in [
        "2147483648-01-01T00:00:00",
        "-2147483649-01-01T00:00:00",
        "2000-01-01T00:00:2147483648",
        "2000-01-01T00:00:00.2147483648",
        "2000-01-01T00:00:00. 2147483647",
    ] {
        assert!(DateTime::parse(value).is_err(), "{value}");
        assert!(
            DateTime::from_format(value, DATETIME_FORMATS[1])
                .unwrap()
                .is_null()
        );
    }
    let mut d = DateTime::default();
    for value in [i32::MAX as u32 + 1, u32::MAX] {
        assert!(d.set_components(1, 1, value, 0, 0, 0).is_err());
        assert!(d.set_date_components(value, 1, 2000).is_err());
        assert!(d.set_time_components(value, 0, 0).is_err());
        assert!(d.is_null());
    }
    let mut d = DateTime::parse("2147483647-12-31T23:59:59.999").unwrap();
    let before = d;
    assert!(d.add_seconds(1).is_err());
    assert_eq!(d, before);
    assert!(d.format(DATETIME_FORMATS[1]).unwrap().len() < 64);
}
#[test]
fn invalid_format_dispatch_and_bounded_text_scanning() {
    assert_eq!(DateTime::default().format("invalid").unwrap(), "");
    let d = DateTime::parse("2000-01-02T03:04:05").unwrap();
    assert!(d.format("invalid").is_err());
    assert!(
        DateTime::from_format("2000-01-02", "invalid")
            .unwrap()
            .is_null()
    );
    // A valid scanf prefix can coexist with bounded ignored trailing text.
    let mut text = String::from("2000-01-02T03:04:05");
    text.extend(std::iter::repeat_n(
        'x',
        MAX_DATETIME_INPUT_BYTES - text.len(),
    ));
    assert_eq!(DateTime::parse(&text).unwrap(), d);
    assert!(DateTime::parse("２０００-01-02T03:04:05").is_err());
    // Source C whitespace includes vertical tab; ordinary Unicode whitespace does not.
    assert_eq!(DateTime::parse("\u{b}2000-01-02T03:04:05").unwrap(), d);
    assert!(DateTime::parse("\u{a0}2000-01-02T03:04:05").is_err());
}
#[test]
fn naive_gregorian_arithmetic_preserves_flags_and_fraction() {
    let mut invalid = DateTime::default();
    invalid.add_seconds(0).unwrap();
    assert_eq!(invalid.components(), (11, 30, -1, 0, 0, 0));
    assert!(invalid.is_null());
    assert_ne!(invalid, DateTime::default());
    assert_eq!(invalid.get(), DateTime::default().get());
    let mut d = DateTime::parse("0001-01-01T00:00:00.123").unwrap();
    d.add_seconds(-1).unwrap();
    assert_eq!(d.get(), "0000-12-31 23:59:59");
    assert!(d.is_valid());
    assert_eq!(d.millisecond(), 123);
    for year in [1, 4, 100, 400, 1900, 2000, 2100, 123456, 2147483646] {
        for month in 1..=12 {
            let mut d = DateTime::default();
            d.set_components(month, 1, year, 12, 34, 56).unwrap();
            let before = d;
            d.add_seconds(1_000_000)
                .unwrap()
                .add_seconds(-1_000_000)
                .unwrap();
            assert_eq!(d, before);
        }
    }
}
#[test]
fn hashing_follows_rendered_millisecond_identity_not_validity_order() {
    fn hash(d: DateTime) -> u64 {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        d.hash(&mut h);
        h.finish()
    }
    let first = DateTime::parse("2000-01-02T03:04:05.1").unwrap();
    let second = DateTime::parse("2000-01-02T03:04:05.2").unwrap();
    assert_eq!(hash(first), hash(first));
    assert_ne!(hash(first), hash(second));
    assert!(first.source_less(&second));
    let mut invalid = DateTime::default();
    invalid.add_seconds(1).unwrap();
    assert_eq!(hash(invalid), hash(DateTime::default()));
}
#[test]
fn current_local_and_utc_clock_values_have_source_seconds_precision() {
    for d in [DateTime::now(), DateTime::now_utc()] {
        assert!(d.is_valid());
        assert!(!d.is_null());
        assert_eq!(d.millisecond(), 0);
        let (m, day, y, h, min, s) = d.components();
        assert!((1..=12).contains(&m));
        assert!((1..=31).contains(&day));
        assert!(y >= 2026);
        assert!((0..24).contains(&h));
        assert!((0..60).contains(&min));
        assert!((0..60).contains(&s));
    }
}
