// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
use openms::data_structures::string_list::*;
use openms::format::text::{ReadOptions, TextFile};
use std::io::Cursor;

fn source_lines() -> Vec<String> {
    TextFile::from_reader(
        Cursor::new(include_bytes!("data/text_file_source.txt")),
        &ReadOptions::default(),
    )
    .unwrap()
    .lines()
    .to_vec()
}
#[test]
fn source_literal_prefix_queries_and_iterator_ranges() {
    let lines = source_lines();
    for (query, expected) in [
        ("first_line", Some(0)),
        ("middle_line", Some(3)),
        ("space_line", None),
        ("tab_line", None),
        ("last_line", Some(10)),
        ("invented_line", None),
        (" ", Some(5)),
        ("\t", Some(6)),
    ] {
        assert_eq!(search_prefix(&lines, query, false), expected);
        assert_eq!(
            search_prefix_in(&lines, 0..lines.len(), query, false).unwrap(),
            expected
        );
    }
    assert_eq!(
        search_prefix_in(&lines, 1..lines.len(), "first_line", false).unwrap(),
        None
    );
    assert_eq!(
        search_prefix_in(&lines, 9..lines.len(), "\t", false).unwrap(),
        None
    );
    for (query, expected) in [
        ("first_line", Some(0)),
        ("space_line", Some(5)),
        ("tab_line", Some(6)),
        ("invented_line", None),
    ] {
        assert_eq!(search_prefix(&lines, query, true), expected);
    }
    assert_eq!(
        search_prefix_in(&lines, 1..lines.len(), "first_line", true).unwrap(),
        None
    );
    let trimmed: Vec<_> = lines
        .iter()
        .map(|x| x.trim_matches([' ', '\t', '\r', '\n']).to_owned())
        .collect();
    assert_eq!(search_prefix(&trimmed, "space_line", false), Some(5));
    assert_eq!(search_prefix(&trimmed, "tab_line", false), Some(6));
    assert_eq!(source_lines(), lines);
}
#[test]
fn source_literal_suffix_queries_and_nonmutating_trim() {
    let lines = source_lines();
    for (query, expected) in [
        ("invented_line", None),
        ("back_space_line", Some(7)),
        ("back_tab_line", Some(8)),
    ] {
        assert_eq!(search_suffix(&lines, query, true), expected);
        assert_eq!(
            search_suffix_in(&lines, 0..lines.len(), query, true).unwrap(),
            expected
        );
        assert_eq!(search_suffix(&lines, query, false), None);
    }
    assert_eq!(
        search_suffix_in(&lines, 8..lines.len(), "back_space_line", true).unwrap(),
        None
    );
    assert_eq!(
        search_suffix(&lines, " \tback_space_line\r\n", true),
        Some(7)
    );
    assert_eq!(search_prefix(&lines, " \ttab_line\r\n", true), Some(6));
    assert_eq!(source_lines(), lines);
}
#[test]
fn source_case_literals_and_native_ascii_locale() {
    let mut values = vec!["yes".into(), "no".into()];
    to_upper(&mut values);
    assert_eq!(values, ["YES", "NO"]);
    let mut values = vec!["yES".into(), "nO".into(), "Äbc αß".into()];
    to_lower(&mut values);
    assert_eq!(values, ["yes", "no", "Äbc αß"]);
    to_upper(&mut values);
    assert_eq!(values, ["YES", "NO", "ÄBC αß"]);
}
#[test]
fn empty_query_trim_range_and_utf8_boundaries() {
    let lines = ["x", "", "\u{a0}y", "\u{b}z", "αβ"];
    assert_eq!(search_prefix(&lines, "", false), Some(0));
    assert_eq!(search_suffix(&lines, " \t\r\n", true), Some(0));
    assert_eq!(search_prefix_in(&lines, 1..5, "", false).unwrap(), Some(1));
    assert_eq!(search_prefix_in(&lines, 5..5, "", false).unwrap(), None);
    assert_eq!(search_suffix_in(&lines, 2..2, "", false).unwrap(), None);
    assert!(search_prefix_in(&lines, 0..6, "", false).is_err());
    let reversed = std::ops::Range { start: 3, end: 2 };
    assert!(search_suffix_in(&lines, reversed, "", false).is_err());
    assert_eq!(search_prefix(&lines, "y", true), None);
    assert_eq!(search_prefix(&lines, "z", true), None);
    assert_eq!(search_suffix(&lines, "β", false), Some(4));
    assert_eq!(search_prefix(&[] as &[String], "", true), None);
}
