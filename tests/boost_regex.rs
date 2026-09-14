// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Tests for `src/concept/boost_regex.rs`, the Boost.Regex facade over
//! `fancy-regex`.
//!
//! The main test is tier 1, an executed differential: `boost_regex_corpus.txt`
//! lists patterns, construction flags, input sets and operations, and
//! `boost_regex_boost.txt` is what Boost.Regex 1.92 printed for them, produced by
//! the oracle driver `../oracle/boost-regex/driver.cpp` (hashes in
//! `tests/data/boost_regex_provenance.json`). Every pattern of the pinned OpenMS
//! sources that uses Boost.Regex is in the corpus, with source-derived and
//! adversarial inputs, plus syntax probes, the review probes and grammar-generated
//! expressions. The facade must reproduce every compile outcome, every match, every
//! group span, every named-group lookup and every token, byte for byte. Which
//! patterns it refuses, and how many cases it compares, is not taken from the
//! facade: `../oracle/boost-regex/refusals.py` derives both from Boost's output and
//! the refusal rules documented in `docs/BOOST_REGEX_SUPPORT.md`.
//!
//! The case-insensitive range, negated-class, nullable-repeat, leading-repeat,
//! open-group backreference, engine-rewrite, atomic-alternation, start-map and
//! escape tests are tier 1 as well, with Boost's answers transcribed from the same
//! oracle output. The class-test section is tier 3 (literals transcribed from the pinned
//! class tests); the limit, work-bound, error and robustness sections are tier 4.

use openms::Error;
use openms::concept::boost_regex::{
    BACKTRACK_LIMIT_BYTES, BoostRegex, MAX_AUTOMATON_ATOMS, MAX_BACKTRACK_SCALE, MAX_GROUP_DEPTH,
    MAX_LOOKBEHIND_WIDTH, MAX_PATTERN_BYTES, MAX_TRANSLATED_BYTES, RegexOptions, SubMatch,
};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::ops::Range;
use std::path::PathBuf;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

fn data(name: &str) -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name);
    std::fs::read_to_string(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// Text field of the corpus: a leading quote, `\\` for a backslash and `\xHH`
/// for every byte outside `0x21..=0x7e`.
fn unescape(field: &str) -> Vec<u8> {
    let bytes = field.as_bytes();
    assert_eq!(bytes.first(), Some(&b'\''), "bad field {field}");
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 1;
    while i < bytes.len() {
        if bytes[i] == b'\\' && bytes.get(i + 1) == Some(&b'\\') {
            out.push(b'\\');
            i += 2;
        } else if bytes[i] == b'\\' && bytes.get(i + 1) == Some(&b'x') {
            out.push(u8::from_str_radix(&field[i + 2..i + 4], 16).expect("hex escape"));
            i += 4;
        } else {
            out.push(bytes[i]);
            i += 1;
        }
    }
    out
}

struct Pattern {
    family: String,
    flags: String,
    names: Vec<String>,
    text: String,
}

struct Run {
    pattern: usize,
    op: String,
    set: String,
}

struct Corpus {
    sets: BTreeMap<String, Vec<Vec<u8>>>,
    patterns: Vec<Pattern>,
    runs: Vec<Run>,
}

impl Corpus {
    fn parse(text: &str) -> Self {
        let mut corpus = Self {
            sets: BTreeMap::new(),
            patterns: Vec::new(),
            runs: Vec::new(),
        };
        let mut family = String::new();
        for line in text.lines() {
            if let Some(name) = line.strip_prefix("# family ") {
                family = name.to_string();
                continue;
            }
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split(' ').collect();
            match fields[0] {
                "S" => {
                    corpus.sets.insert(
                        fields[1].to_string(),
                        fields[2..].iter().map(|field| unescape(field)).collect(),
                    );
                }
                "P" => {
                    assert_eq!(fields[1].parse::<usize>().unwrap(), corpus.patterns.len());
                    corpus.patterns.push(Pattern {
                        family: family.clone(),
                        flags: fields[2].to_string(),
                        names: if fields[3] == "-" {
                            Vec::new()
                        } else {
                            fields[3].split(',').map(str::to_string).collect()
                        },
                        text: String::from_utf8(unescape(fields[4])).expect("patterns are UTF-8"),
                    });
                }
                "R" => corpus.runs.push(Run {
                    pattern: fields[1].parse().unwrap(),
                    op: fields[2].to_string(),
                    set: fields[3].to_string(),
                }),
                other => panic!("unknown corpus line kind {other}"),
            }
        }
        corpus
    }
}

fn render_sub(token: &SubMatch) -> String {
    if token.matched {
        format!("{},{}", token.range.start, token.range.end)
    } else if token.range.is_empty() {
        format!("~{}", token.range.start)
    } else {
        format!("~{}!{}", token.range.start, token.range.end)
    }
}

/// The facade's answer in the driver's output format.
fn render(regex: &BoostRegex, op: &str, input: &[u8], names: &[String]) -> String {
    if op == "s" || op == "m" {
        let (found, flag) = if op == "s" {
            (regex.search(input), regex.is_search_match(input))
        } else {
            (regex.full_match(input), regex.is_full_match(input))
        };
        let captures = match found {
            Err(error) => return format!("RTERR {error}"),
            Ok(captures) => captures,
        };
        assert_eq!(
            flag.ok(),
            Some(captures.is_some()),
            "boolean and capturing calls disagree for {:?} on {input:?}",
            regex.as_str()
        );
        let Some(captures) = captures else {
            return "0".to_string();
        };
        let mut out = String::from("1 ");
        for group in 0..captures.len() {
            if group > 0 {
                out.push(';');
            }
            match captures.get(group) {
                Some(range) => {
                    let _ = write!(out, "{},{}", range.start, range.end);
                }
                None => {
                    let _ = write!(out, "~{}", input.len());
                }
            }
        }
        if !names.is_empty() {
            out.push(' ');
            for (index, name) in names.iter().enumerate() {
                if index > 0 {
                    out.push(';');
                }
                match captures.name(name) {
                    Some(range) => {
                        let _ = write!(out, "{},{}", range.start, range.end);
                    }
                    None => {
                        let _ = write!(out, "~{}", captures.range().end);
                    }
                }
            }
        }
        return out;
    }
    let submatches: Vec<i32> = op[1..]
        .split(',')
        .map(|index| index.parse().unwrap())
        .collect();
    let tokens = match regex.tokens(input, &submatches) {
        Err(error) => return format!("RTERR {error}"),
        Ok(tokens) => tokens,
    };
    let mut rendered = Vec::new();
    for token in tokens {
        match token {
            Ok(token) => rendered.push(render_sub(&token)),
            Err(error) => return format!("RTERR {error}"),
        }
    }
    if rendered.is_empty() {
        "-".to_string()
    } else {
        rendered.join(";")
    }
}

/// Patterns outside the fuzz family that Boost compiles and the facade refuses
/// with `Error::Unsupported`, as `(flags, pattern)`: the output of
/// `../oracle/boost-regex/refusals.py` for the committed fixture, which applies the
/// documented rules to Boost's output. `docs/BOOST_REGEX_SUPPORT.md` documents each
/// construct.
const EXPECTED_UNSUPPORTED: &[(&str, &str)] = &[
    ("-", "((?=a))*"),
    ("-", "((a)|\\1b)x"),
    ("-", "()*"),
    ("-", "()+"),
    ("-", "()a{1,3}?\\b"),
    ("-", "(*ACCEPT)"),
    ("-", "(*FAIL)"),
    ("-", "(?!(?=a+))b"),
    ("-", "(?!(a))b"),
    ("-", "(?!.*b)a"),
    ("-", "(?!a*$)x"),
    ("-", "(?!a+b)a"),
    ("-", "(?!a{2})a"),
    ("-", "(?!x)a{1,3}?\\b"),
    ("-", "(?#c)a{1,3}?\\b"),
    ("-", "(?(1)a|b)"),
    ("-", "(?-i)a{1,3}?\\b"),
    ("-", "(?-s).{1,3}?\\b"),
    ("-", "(?:$)+"),
    ("-", "(?:$|a){2}"),
    (
        "-",
        "(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:a+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+)+",
    ),
    (
        "-",
        "(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:a+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+$)+",
    ),
    (
        "-",
        "(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:(?:a{2,}){2,}){2,}){2,}){2,}){2,}){2,}){2,}){2,}){2,}){2,}){2,}){2,}){2,}){2,}){2,}){2,}){2,}){2,}){2,}){2,}",
    ),
    ("-", "(?:(?:(?=a)(?=a)){999999999}){999999999}"),
    ("-", "(?:(?:(?=a)a{999999999}){999999999}){999999999}"),
    ("-", "(?:(?:\\bb{999999999}){999999999}){999}"),
    ("-", "(?:(?:ab)?|a)*"),
    ("-", "(?:(?:a{0}(?=a)){99999}){99999}"),
    ("-", "(?:(?:a{0}\\b){999999999}){999999999}"),
    ("-", "(?:(?:a{999999999}){999999999}){999999999}"),
    ("-", "(?:(?:a|b?)*?)+b"),
    ("-", "(?:(?=(a))|b)*"),
    ("-", "(?:(?=a)(?=a))+b"),
    ("-", "(?:(?=a)(?=a)){99999999}"),
    ("-", "(?:(?=a)|a)+?b"),
    ("-", "(?:(?=a)|a){1,}"),
    ("-", "(?:(?=a)|a){2,3}"),
    ("-", "(?:(?=a)|a){2,}"),
    ("-", "(?:(?=a)|a){2}"),
    ("-", "(?:(?=a)|b){2}"),
    ("-", "(?:(?=b)|a)+"),
    ("-", "(?:(?i)x|y+)*C"),
    ("-", "(?:(?i:a|b+)c)+"),
    ("-", "(?:(a|\\1b)c)+"),
    ("-", "(?:(a|b\\1)c)+"),
    ("-", "(?:)a{1,3}?\\b"),
    ("-", "(?:.{0,9999})*?b"),
    ("-", "(?:.{0,999})*?b"),
    ("-", "(?:\\.|\\<b+)*c"),
    ("-", "(?:\\B|a){2}"),
    ("-", "(?:\\W\\W??|a)+\\>"),
    ("-", "(?:\\W\\W??|a){2}\\>"),
    ("-", "(?:\\b)*"),
    ("-", "(?:\\b|a)+?b"),
    ("-", "(?:\\b|a){1,2}b"),
    ("-", "(?:\\b|a){2,}b"),
    ("-", "(?:\\b|a){2}b"),
    ("-", "(?:\\b|a){3}"),
    ("-", "(?:\\w\\w+){2}\\>"),
    ("-", "(?:\\w{2,}){2,}\\>"),
    ("-", "(?:^)*"),
    ("-", "(?:^|a){2,}b"),
    ("-", "(?:a*){2}b"),
    ("-", "(?:a*?b?)*"),
    ("-", "(?:a?){0,2}"),
    ("-", "(?:a?){2}?b"),
    ("-", "(?:a?){2}b"),
    ("-", "(?:a?\\b){2}"),
    ("-", "(?:a?b?)+?c"),
    ("-", "(?:ab){999999999}"),
    ("-", "(?:a{0}$){99999999}"),
    ("-", "(?:a{0}\\b)*a"),
    ("-", "(?:a{0}\\b)+a"),
    ("-", "(?:a{0}\\b){0,999999999}"),
    ("-", "(?:a{0}\\b){2}a"),
    ("-", "(?:a{0}\\b){999999999,}"),
    ("-", "(?:a{0}\\b){999999999}"),
    ("-", "(?:a|(?=b)){2}b"),
    ("-", "(?:a|){2}b"),
    ("-", "(?:a|\\b)*"),
    ("-", "(?:a|\\b){1,}"),
    ("-", "(?:a|\\b){2}b"),
    ("-", "(?:a|ab)+\\<"),
    ("-", "(?:b?|a)*"),
    ("-", "(?:b?|a)*?"),
    ("-", "(?:b?|a)+"),
    ("-", "(?:b?|a){0,}"),
    ("-", "(?:b?|a){1,}"),
    ("-", "(?:x{0}$)+"),
    ("-", "(?:x{0}(?=a)){3}a"),
    ("-", "(?:y|(?i:b.+))+C"),
    ("-", "(?:|a){2}b"),
    ("-", "(?<!(?=.*b)){1000,}"),
    ("-", "(?<!(?=.*b)\\w)[ab](?<=a)b"),
    ("-", "(?<!(?=.*b)a)b"),
    ("-", "(?<!(?=a+)a)b"),
    ("-", "(?<!a{255})b"),
    ("-", "(?<!a{256})b"),
    ("-", "(?<!a{2})b"),
    ("-", "(?<=(a))b"),
    ("-", "(?<=(a)|b)c"),
    ("-", "(?<=\\Z)"),
    ("-", "(?<=a)\\Z"),
    ("-", "(?<=a{100}b{156})c"),
    ("-", "(?<=a{256})b"),
    ("-", "(?<n-x>a)"),
    ("-", "(?<n>a)(?P>n)"),
    ("-", "(?<n>a)\\k<n>"),
    ("-", "(?<n>a|\\1b)x"),
    ("-", "(?=(\\d+))\\d"),
    ("-", "(?=(a))*"),
    ("-", "(?=(a))*?b"),
    ("-", "(?=(a)){0}"),
    ("-", "(?=(a)){2}"),
    ("-", "(?=a)a{1,3}?\\b"),
    ("-", "(?>(?:ab)+)c"),
    ("-", "(?>a*$)b"),
    ("-", "(?>a*)b"),
    ("-", "(?>a+)b"),
    ("-", "(?>a?){2}"),
    ("-", "(?>a{2,})b"),
    ("-", "(?>b?|a)*"),
    ("-", "(?Pim-sx:a|b)"),
    ("-", "(?Px)a b"),
    ("-", "(?i)(?-i:b.+)?[a-z]"),
    ("-", "(?i)(?:a(?-i:b.+))+c"),
    ("-", "(?i)*"),
    ("-", "(?i)+a"),
    ("-", "(?i)\\<A"),
    ("-", "(?i)\\<A|b"),
    ("-", "(?i)\\<[A-Z]"),
    ("-", "(?i:b.+)*C"),
    ("-", "(?i:b.+){1}C"),
    ("-", "(?i:bc+)+C"),
    ("-", "(?s)a{1,3}?\\b"),
    ("-", "(?x)a b"),
    ("-", "(?|(a)|(b))"),
    ("-", "(\\w+?)+\\>"),
    ("-", "(a(?=\\1))b"),
    ("-", "(a)(?(1)b|c)"),
    ("-", "(a)(?1)"),
    ("-", "(a)(?:\\1){2}"),
    ("-", "(a)(b)(c)(d)(e)(f)(g)(h)(i)(j)\\10"),
    ("-", "(a)(b|\\2a)x"),
    ("-", "(a)\\10"),
    ("-", "(a)\\g1"),
    ("-", "(a)\\g{-1}"),
    ("-", "(a)\\g{1}"),
    ("-", "(a*)*"),
    ("-", "(a*)+"),
    ("-", "(a?)\\1{0,999999999}"),
    ("-", "(a?)\\1{2}"),
    ("-", "(a?\\1*)\\b"),
    ("-", "(a\\1)"),
    ("-", "(a\\1?)+"),
    ("-", "(a\\1{2})"),
    ("-", "(ab|\\1b)x"),
    ("-", "(a{1,3}?)\\b"),
    ("-", "(a{999999999})\\1{999999999}"),
    ("-", "(a|)\\1{1,3}b"),
    ("-", "(a|\\1*)$"),
    ("-", "(a|\\1+)b"),
    ("-", "(a|\\1a)b"),
    ("-", "(a|\\1b)x"),
    ("-", "(a|b\\1)+"),
    ("-", "(x|a|\\1\\1)b"),
    ("-", "(|a)*"),
    ("-", ".{1,3}?(?=b)"),
    ("-", ".{3,5}?\\b"),
    ("-", "[!-[.-.][.-.]]"),
    ("-", "[!-[.].]]"),
    ("-", "[-[.a.]]"),
    ("-", "[A-[.a.]]"),
    ("-", "[A-[.a.]x]"),
    ("-", "[[.a.]-z]"),
    ("-", "[[.a.]]"),
    ("-", "[[:<:]]"),
    ("-", "[[:>:]]"),
    ("-", "[[:UNICODE:]]"),
    ("-", "[[:Unicode:]]"),
    ("-", "[[:^unicode:]]"),
    ("-", "[[:unicode:][:alpha:]]"),
    ("-", "[[:unicode:]]"),
    ("-", "[[=a=]]"),
    ("-", "[\\0]"),
    ("-", "[\\n-[.-.]]"),
    ("-", "[\\x 4]"),
    ("-", "[\\x41-[.a.]]"),
    ("-", "[\\y]"),
    ("-", "[^A-[.a.]]"),
    ("-", "[^[:unicode:]]"),
    ("-", "[^b]{0,2}?b"),
    ("-", "[a-[.a.]]"),
    ("-", "[a-\\d]"),
    ("-", "[a[.b.]c]"),
    ("-", "[a[:unicode:]]"),
    ("-", "[q[.a.]-[.z.]]"),
    ("-", "[é]"),
    ("-", "\\0"),
    ("-", "\\012"),
    ("-", "\\1{2}(a)"),
    ("-", "\\<(?:\\w\\w*)+"),
    ("-", "\\<a(?i)b"),
    ("-", "\\C"),
    ("-", "\\E"),
    ("-", "\\G"),
    ("-", "\\K"),
    ("-", "\\L"),
    ("-", "\\P{L}"),
    ("-", "\\Q"),
    ("-", "\\Qa.b\\E"),
    ("-", "\\Qab"),
    ("-", "\\R"),
    ("-", "\\U"),
    ("-", "\\X"),
    ("-", "\\X41"),
    ("-", "\\Z"),
    ("-", "\\Z\\n?"),
    ("-", "\\ba{1,3}?\\b"),
    ("-", "\\cA"),
    ("-", "\\cZ"),
    ("-", "\\ca"),
    ("-", "\\j"),
    ("-", "\\l"),
    ("-", "\\pL"),
    ("-", "\\p{L}"),
    ("-", "\\q"),
    ("-", "\\u"),
    ("-", "\\w{3,5}?\\b"),
    ("-", "\\x 4"),
    ("-", "\\x+4"),
    ("-", "\\x-0"),
    ("-", "\\x80"),
    ("-", "\\xc3"),
    ("-", "\\xff"),
    ("-", "\\x{	41}"),
    ("-", "\\x{ 41}"),
    ("-", "\\x{+41}"),
    ("-", "\\x{-0}"),
    ("-", "\\x{0X41}"),
    ("-", "\\x{0x41}"),
    ("-", "\\y"),
    ("-", "^\\Z"),
    ("-", "^a{1,3}?\\b"),
    (
        "-",
        "a$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$",
    ),
    ("-", "a(?i)*"),
    ("-", "a*+"),
    ("-", "a++"),
    ("-", "a?+"),
    ("-", "a\\Kb"),
    ("-", "a\\Z"),
    ("-", "a\\Z\\n"),
    ("-", "a{0,2}?\\b"),
    ("-", "a{1,3}?(?:\\b|x)"),
    ("-", "a{1,3}?\\b"),
    ("-", "a{1000000000}"),
    ("-", "a{2}+"),
    ("-", "a٣"),
    (
        "-",
        "x*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$",
    ),
    ("-", "é"),
    ("-", "é+"),
    ("i", "(?-i)(?:\\<)A"),
    ("i", "(a|\\1b)x"),
    ("i", "[A-[.a.]]"),
    ("i", "[ab]{3,5}?(?!a)"),
    ("i", "a{1,3}?\\b"),
];

/// Refusals over the whole corpus, fuzz family included, by the construct the
/// facade names in its `Error::Unsupported` message; also from `refusals.py`.
const EXPECTED_REFUSALS: &[(&str, usize)] = &[
    (
        "\\< in an expression that switches case sensitivity (Boost builds its start map with the wrong case)",
        11,
    ),
    (
        "\\< or \\> in an expression with a repeated group holding a repeat or alternation (Boost's start map drops the bytes that start another iteration)",
        17,
    ),
    ("\\G", 1),
    ("\\K", 2),
    ("\\Q...\\E quoting", 4),
    (
        "\\Z (Boost never tries it at a form feed when it starts the expression)",
        7,
    ),
    (
        "a backreference inside the group it refers to (Boost compares the span an abandoned attempt at the group left)",
        29,
    ),
    ("a backreference with more than one digit", 2),
    ("a backtracking control verb", 2),
    ("a branch-reset group", 1),
    (
        "a capturing group inside a lookaround or atomic group (Boost keeps its capture when the surrounding match backtracks)",
        84,
    ),
    ("a character class name that is not a POSIX class", 2),
    ("a character property escape", 3),
    ("a collating element", 14),
    ("a conditional expression", 2),
    ("a control-character escape", 3),
    (
        "a counted repeat of a backreference that can match the empty string (Boost ends a repeat after an empty iteration)",
        7,
    ),
    ("a group name outside [A-Za-z0-9_]", 1),
    ("a hexadecimal escape above 0x7F", 3),
    (
        "a hexadecimal escape whose digits start with a sign, white space or 0x (Boost reads them with std::istream)",
        10,
    ),
    (
        "a lazy repeat with a finite maximum of a one-byte atom that starts the expression (Boost's leading-repeat optimization skips start positions)",
        27,
    ),
    ("a lookbehind wider than MAX_LOOKBEHIND_WIDTH bytes", 2),
    ("a named or relative backreference", 4),
    ("a named-character, line-ending or grapheme escape", 4),
    ("a non-ASCII byte", 4),
    ("a possessive quantifier", 4),
    ("a recursive sub-expression", 2),
    ("a repeat bound above MAX_REPEAT", 1),
    (
        "a repeat inside an atomic group or a negative lookaround (the engine discards its work there without counting it)",
        93,
    ),
    (
        "a repeat of a group holding a repeat or alternation under other case sensitivity (Boost builds its start map with the wrong case)",
        20,
    ),
    (
        "a repeat of a group that can match the empty string and captures (Boost records a final empty iteration)",
        26,
    ),
    ("a repeat of a group that only asserts a position", 11),
    (
        "a repeat of a modifier group that switches case sensitivity (Boost undoes the switch when the empty repetition is abandoned)",
        3,
    ),
    (
        "a repeat of more than one iteration of a group that can match the empty string (Boost ends a repeat after an empty iteration)",
        73,
    ),
    (
        "a repeat whose shortest match is longer than MAX_REPEAT bytes",
        5,
    ),
    ("an equivalence class", 1),
    (
        "an escape inside a character class without a translated meaning",
        3,
    ),
    ("an escape letter without a translated meaning", 7),
    ("an octal escape", 2),
    (
        "groups and repeats nested deeper than the engine's parser allows once translated",
        2,
    ),
    (
        "more line anchors, word-boundary assertions, alternations and repeated groups than Boost's start-map recursion limit allows (Boost throws error_complexity)",
        3,
    ),
    ("the [:unicode:] class", 7),
    ("the x (extended) modifier", 3),
];

/// Cases compared against Boost: the fixture's size, so a truncated fixture fails
/// (`refusals.py`).
const EXPECTED_CASES: usize = 162_714;

const NULLABLE_GROUP: &str = "a repeat of more than one iteration of a group that can match the empty string (Boost ends a repeat after an empty iteration)";
const LEADING_LAZY: &str = "a lazy repeat with a finite maximum of a one-byte atom that starts the expression (Boost's leading-repeat optimization skips start positions)";
const DISCARDED: &str = "a repeat inside an atomic group or a negative lookaround (the engine discards its work there without counting it)";
const LONG_SHORTEST_MATCH: &str = "a repeat whose shortest match is longer than MAX_REPEAT bytes";
const OPEN_GROUP_BACKREFERENCE: &str = "a backreference inside the group it refers to (Boost compares the span an abandoned attempt at the group left)";

/// The construct an `Error::Unsupported` message names.
fn refusal_category(message: &str) -> String {
    match message.rfind(" uses ") {
        Some(index) => {
            let rest = &message[index + " uses ".len()..];
            rest.split(" at byte ").next().unwrap_or(rest).to_string()
        }
        None => "an engine refusal".to_string(),
    }
}

/// The category of a refusal, or a panic naming what happened instead.
fn refused(pattern: &str, result: Result<BoostRegex, Error>) -> String {
    match result {
        Err(Error::Unsupported(message)) => refusal_category(&message),
        Err(other) => panic!("{pattern:?}: expected Unsupported, got {other}"),
        Ok(_) => panic!("{pattern:?}: expected Unsupported, but it compiled"),
    }
}

/// Run `work` on its own thread and fail if it has not finished after `seconds`,
/// so a hang fails this test instead of stalling the whole run.
fn within<T: Send + 'static>(
    seconds: u64,
    what: &str,
    work: impl FnOnce() -> T + Send + 'static,
) -> T {
    let (sender, receiver) = mpsc::channel();
    thread::spawn(move || {
        let _ = sender.send(work());
    });
    match receiver.recv_timeout(Duration::from_secs(seconds)) {
        Ok(value) => value,
        Err(mpsc::RecvTimeoutError::Timeout) => {
            panic!("{what} did not finish within {seconds} s")
        }
        Err(mpsc::RecvTimeoutError::Disconnected) => panic!("{what} panicked"),
    }
}

#[test]
fn boost_differential_corpus_has_no_mismatch() {
    let corpus = Corpus::parse(&data("boost_regex_corpus.txt"));
    let boost_text = data("boost_regex_boost.txt");
    let mut boost = boost_text.lines().filter(|line| !line.starts_with('#'));

    let mut mismatches: Vec<String> = Vec::new();
    let mut unsupported: BTreeSet<(String, String)> = BTreeSet::new();
    let mut refusals: BTreeMap<String, usize> = BTreeMap::new();
    let mut compiled: Vec<Option<BoostRegex>> = Vec::with_capacity(corpus.patterns.len());
    let mut boost_compiled: Vec<bool> = Vec::with_capacity(corpus.patterns.len());
    let mut families: BTreeMap<&str, [usize; 4]> = BTreeMap::new();

    for pattern in &corpus.patterns {
        let expected = boost
            .next()
            .expect("Boost output ends before the patterns do");
        let options = RegexOptions {
            icase: pattern.flags.contains('i'),
            no_mod_s: pattern.flags.contains('S'),
            ..RegexOptions::default()
        };
        let result = BoostRegex::with_options(&pattern.text, options);
        let stats = families.entry(pattern.family.as_str()).or_default();
        stats[0] += 1;
        let boost_ok = expected.starts_with("OK ");
        boost_compiled.push(boost_ok);
        match (boost_ok, result) {
            (true, Ok(regex)) => {
                let marks: usize = expected[3..].parse().unwrap();
                if regex.mark_count() != marks {
                    mismatches.push(format!(
                        "{:?} flags={}: mark_count {} vs Boost {marks}",
                        pattern.text,
                        pattern.flags,
                        regex.mark_count()
                    ));
                }
                compiled.push(Some(regex));
            }
            (true, Err(Error::Unsupported(message))) => {
                *refusals.entry(refusal_category(&message)).or_default() += 1;
                if pattern.family != "FUZZ" {
                    unsupported.insert((pattern.flags.clone(), pattern.text.clone()));
                }
                stats[1] += 1;
                compiled.push(None);
            }
            (false, Err(_)) => {
                stats[2] += 1;
                compiled.push(None);
            }
            (true, Err(error)) => {
                mismatches.push(format!(
                    "{:?} flags={}: Boost compiles, facade: {error}",
                    pattern.text, pattern.flags
                ));
                compiled.push(None);
            }
            (false, Ok(_)) => {
                mismatches.push(format!(
                    "{:?} flags={}: Boost refuses ({expected}), facade compiles",
                    pattern.text, pattern.flags
                ));
                compiled.push(None);
            }
        }
    }

    let mut cases = 0usize;
    let mut grouped: BTreeMap<String, (usize, String)> = BTreeMap::new();
    for run in &corpus.runs {
        if !boost_compiled[run.pattern] {
            continue;
        }
        let pattern = &corpus.patterns[run.pattern];
        for input in &corpus.sets[&run.set] {
            let expected = boost.next().expect("Boost output ends before the cases do");
            let Some(regex) = &compiled[run.pattern] else {
                continue;
            };
            cases += 1;
            families.get_mut(pattern.family.as_str()).unwrap()[3] += 1;
            let got = render(regex, &run.op, input, &pattern.names);
            if got != expected {
                let key = format!(
                    "{} {:?} flags={} op={}",
                    pattern.family, pattern.text, pattern.flags, run.op
                );
                let entry = grouped.entry(key).or_insert_with(|| {
                    (
                        0,
                        format!(
                            "input={:?}: facade {got} vs Boost {expected}",
                            String::from_utf8_lossy(input)
                        ),
                    )
                });
                entry.0 += 1;
                mismatches.push(String::new());
            }
        }
    }
    for (key, (count, example)) in grouped.iter().take(400) {
        eprintln!("MISMATCH {count:6} {key}  e.g. {example}");
    }
    assert!(
        boost.next().is_none(),
        "Boost output has lines the corpus does not account for"
    );

    for (family, [patterns, refused, both_error, compared]) in &families {
        eprintln!(
            "{family:8} patterns {patterns:5}  unsupported {refused:3}  both-error {both_error:4}  compared cases {compared}"
        );
    }
    eprintln!("total compared cases {cases}");

    let shown: Vec<&String> = mismatches.iter().filter(|line| !line.is_empty()).collect();
    assert!(
        mismatches.is_empty(),
        "{} mismatches against Boost, first ones:\n{}",
        mismatches.len(),
        shown
            .iter()
            .map(|line| line.as_str())
            .collect::<Vec<_>>()
            .join("\n")
    );

    let expected: BTreeSet<(String, String)> = EXPECTED_UNSUPPORTED
        .iter()
        .map(|(flags, pattern)| ((*flags).to_string(), (*pattern).to_string()))
        .collect();
    let listing: String = unsupported
        .iter()
        .map(|(flags, pattern)| format!("    ({flags:?}, {pattern:?}),\n"))
        .collect();
    let categories: String = refusals
        .iter()
        .map(|(category, count)| format!("    ({category:?}, {count}),\n"))
        .collect();
    assert_eq!(
        unsupported, expected,
        "refused non-fuzz patterns differ from the oracle's list:\n{listing}"
    );
    let expected_refusals: BTreeMap<String, usize> = EXPECTED_REFUSALS
        .iter()
        .map(|(category, count)| ((*category).to_string(), *count))
        .collect();
    assert_eq!(
        refusals, expected_refusals,
        "refusal categories differ from the oracle's counts:\n{categories}"
    );
    assert_eq!(cases, EXPECTED_CASES, "the fixture lost or gained cases");
}

/// Under `icase` Boost lower-cases both endpoints of a range before ordering them
/// and lower-cases each input byte before its set lookup
/// (`basic_regex_creator::append_set`, `perl_matcher::match_set`), so a range that
/// spans non-letters matches differently than a folded byte range would. Compile
/// outcomes and `regex_search` answers transcribed from Boost's output for the
/// `ADV` family of the fixture.
#[test]
fn case_insensitive_ranges_follow_boost() {
    let icase = RegexOptions {
        icase: true,
        ..RegexOptions::default()
    };
    let plain = RegexOptions::default();
    // (options, pattern, one-byte or short haystack, whether Boost matches all of it)
    let cases: &[(RegexOptions, &str, &[u8], bool)] = &[
        (icase, "[@-Z]", b"[", true),
        (icase, "[@-Z]", b"\\", true),
        (icase, "[@-Z]", b"]", true),
        (icase, "[@-Z]", b"@", true),
        (icase, "[@-Z]", b"a", true),
        (icase, "[@-Z]", b"_", true),
        (icase, "[@-Z]", b"`", true),
        (icase, "[A-z]", b"[", false),
        (icase, "[A-z]", b"_", false),
        (icase, "[A-z]", b"a", true),
        (icase, "[A-z]", b"Z", true),
        (icase, "[A-z]", b"`", false),
        (icase, "[a-Z]", b"a", true),
        (icase, "[a-Z]", b"Z", true),
        (icase, "[a-Z]", b"[", false),
        (icase, "[a-Z]", b"k", true),
        (icase, "[k-K]", b"k", true),
        (icase, "[k-K]", b"K", true),
        (icase, "[k-K]", b"z", false),
        (plain, "[Z-a]", b"Z", true),
        (plain, "[Z-a]", b"[", true),
        (plain, "[Z-a]", b"a", true),
        (plain, "(?i)[@-Z]", b"[", true),
        (plain, "(?i)[@-Z]", b"z", true),
        (icase, "x[@-Z]y", b"x_y", true),
        (icase, "x[@-Z]y", b"x[y", true),
        (icase, "[^A-z]", b"[", true),
        (icase, "[^A-z]", b"a", false),
        (icase, "[^A-z]", b"1", true),
    ];
    for &(options, pattern, haystack, matches) in cases {
        let regex = BoostRegex::with_options(pattern, options)
            .unwrap_or_else(|error| panic!("{pattern:?}: {error}"));
        let found = regex
            .search(haystack)
            .unwrap()
            .map(|captures| captures.range());
        assert_eq!(
            found,
            matches.then_some(0..haystack.len()),
            "{pattern:?} icase={} on {haystack:?}",
            options.icase
        );
    }
    // Boost's error_range (code 4): the lower-cased endpoints are out of order.
    for (options, pattern) in [
        (icase, "[Z-a]"),
        (icase, "[A-_]"),
        (icase, r"[\x41-\x5b]"),
        (plain, "[a-Z]"),
    ] {
        assert!(
            matches!(
                BoostRegex::with_options(pattern, options),
                Err(Error::InvalidValue(_))
            ),
            "{pattern:?} icase={}",
            options.icase
        );
    }
}

/// Boost ORs every negated class of one bracket expression into one mask and
/// complements it once, so `[\S\D]` matches only bytes that are neither space nor
/// digit (`basic_regex_creator::append_set`). Answers transcribed from Boost's
/// output for the `ADV` family of the fixture.
#[test]
fn negated_classes_are_complemented_together() {
    // (pattern, icase, haystack, whether Boost matches it)
    let cases: &[(&str, bool, &[u8], bool)] = &[
        (r"[\S\D]", false, b"1", false),
        (r"[\S\D]", false, b" ", false),
        (r"[\S\D]", false, b"\t", false),
        (r"[\S\D]", false, b"a", true),
        (r"[\S\D]", true, b"1", false),
        (r"[\S\D]", true, b"a", true),
        ("[[:^lower:][:^upper:]]", false, b"a", false),
        ("[[:^lower:][:^upper:]]", false, b"A", false),
        ("[[:^lower:][:^upper:]]", false, b"1", true),
        (r"[\W\S]", false, b" ", false),
        (r"[\W\S]", false, b"a", false),
        (r"[\W\S]", false, b"_", false),
        (r"[^\S\D]", false, b" ", true),
        (r"[^\S\D]", false, b"1", true),
        (r"[^\S\D]", false, b"a", false),
        (r"[\D[:^alpha:]x]", false, b"1", false),
        (r"[\D[:^alpha:]x]", false, b"a", false),
        (r"[\D[:^alpha:]x]", false, b"_", true),
        // A single negated class is complemented on its own, as before.
        (r"[\S\d]", false, b"1", true),
        (r"[\S\d]", false, b" ", false),
        (r"[\S\d]", false, b"a", true),
    ];
    for &(pattern, icase, haystack, matches) in cases {
        let options = RegexOptions {
            icase,
            ..RegexOptions::default()
        };
        let regex = BoostRegex::with_options(pattern, options).unwrap();
        let found = regex
            .search(haystack)
            .unwrap()
            .map(|captures| captures.range());
        assert_eq!(
            found,
            matches.then_some(0..haystack.len()),
            "{pattern:?} icase={icase} on {haystack:?}"
        );
    }
}

fn tokens_text(regex: &BoostRegex, haystack: &str, submatches: &[i32]) -> Vec<String> {
    regex
        .tokens(haystack.as_bytes(), submatches)
        .unwrap()
        .map(|token| {
            String::from_utf8(token.unwrap().as_bytes(haystack.as_bytes()).to_vec()).unwrap()
        })
        .collect()
}

/// `SpectrumNativeIDParser_test`, `SpectrumLookup_test` and
/// `SpectrumMetaDataLookup_test`: the regular-expression half of the expected
/// values. `extractScanNumber` takes the last group-1 token and converts it.
#[test]
fn class_test_expressions() {
    let spectrum = BoostRegex::new(r"spectrum=(?<SCAN>\d+)").unwrap();
    assert_eq!(tokens_text(&spectrum, "spectrum=42", &[1]), ["42"]);
    assert_eq!(tokens_text(&spectrum, "spectrum=0", &[1]), ["0"]);
    assert_eq!(tokens_text(&spectrum, "spectrum=99999", &[1]), ["99999"]);
    assert!(tokens_text(&spectrum, "scan=42", &[1]).is_empty());
    assert!(tokens_text(&spectrum, "no_match_here", &[1]).is_empty());
    assert!(tokens_text(&spectrum, "", &[1]).is_empty());

    let scan = BoostRegex::new(r"scan=(?<SCAN>\d+)").unwrap();
    assert_eq!(tokens_text(&scan, "scan=123", &[1]), ["123"]);
    assert_eq!(
        tokens_text(&scan, "controllerType=0 controllerNumber=1 scan=456", &[1]),
        ["456"]
    );
    let multi = BoostRegex::new(r"(?<SCAN>\d+)").unwrap();
    assert_eq!(tokens_text(&multi, "1 2 3 42", &[1]).last().unwrap(), "42");
    let index = BoostRegex::new(r"index=(?<SCAN>\d+)").unwrap();
    assert_eq!(tokens_text(&index, "index=0", &[1]), ["0"]);

    let from_native_id = BoostRegex::new(r"scan=(?<GROUP>\d+)").unwrap();
    assert!(from_native_id.search(b"scan=456").unwrap().is_some());

    let default_scan = BoostRegex::new(r"=(?<SCAN>\d+)$").unwrap();
    assert_eq!(tokens_text(&default_scan, "spectrum=2", &[1]), ["2"]);

    let rt_mz = BoostRegex::new(r"rt=(?<RT>\d+(\.\d+)?),mz=(?<MZ>\d+(\.\d+)?)").unwrap();
    let text = b"rt=5.0,mz=1000.0";
    let captures = rt_mz.search(text).unwrap().unwrap();
    assert_eq!(&text[captures.name("RT").unwrap()], b"5.0");
    assert_eq!(&text[captures.name("MZ").unwrap()], b"1000.0");
    assert!(captures.name("SCAN").is_none());
}

#[test]
fn line_separators_follow_boost() {
    // `$` before \r, \n and \f, but never between \r and \n.
    let dot_end = BoostRegex::new(".$").unwrap();
    assert_eq!(dot_end.search(b"\r\n").unwrap().unwrap().range(), 1..2);
    let scan = BoostRegex::new(r"=(?<SCAN>\d+)$").unwrap();
    for text in [
        "scanId=21\r",
        "scanId=21\n",
        "scanId=21\x0c",
        "scanId=21\r\n",
    ] {
        assert_eq!(tokens_text(&scan, text, &[1]), ["21"], "{text:?}");
    }
    // ASCII digits only: U+0663 is not \d.
    assert!(!scan.is_search_match("=\u{663}".as_bytes()).unwrap());
    // `.` matches a newline unless no_mod_s or (?-s).
    assert!(
        BoostRegex::new("a.b")
            .unwrap()
            .is_full_match(b"a\nb")
            .unwrap()
    );
    let no_mod_s = RegexOptions {
        no_mod_s: true,
        ..RegexOptions::default()
    };
    assert!(
        !BoostRegex::with_options("a.b", no_mod_s)
            .unwrap()
            .is_full_match(b"a\nb")
            .unwrap()
    );
}

#[test]
fn token_iteration_after_empty_matches() {
    // Enzyme expressions are zero-width: Boost yields an empty first token when
    // the cut is at position 0, and never re-matches the same empty position.
    let trypsin = BoostRegex::new(r"(?<=[KRX])(?!P)").unwrap();
    let tokens: Vec<_> = trypsin
        .tokens(b"AKRPDKA", &[-1])
        .unwrap()
        .map(Result::unwrap)
        .collect();
    let lengths: Vec<usize> = tokens.iter().map(|token| token.range.len()).collect();
    assert_eq!(lengths, [2, 4, 1]);
    let asp_n = BoostRegex::new(r"(?=[DBX])").unwrap();
    let first = asp_n
        .tokens(b"DAD", &[-1])
        .unwrap()
        .next()
        .unwrap()
        .unwrap();
    assert_eq!(
        first,
        SubMatch {
            range: 0..0,
            matched: false
        }
    );
    // After an empty match a non-empty alternative at the same position still wins.
    let mixed = BoostRegex::new("(?=K)|K").unwrap();
    let spans: Vec<_> = mixed
        .tokens(b"AK", &[0])
        .unwrap()
        .map(|token| token.unwrap().range)
        .collect();
    assert_eq!(spans, [1..1, 1..2]);
}

#[test]
fn refuses_what_boost_refuses() {
    for pattern in [
        "(?<=K|RR)",
        "(?<=[KR]+)",
        "(?<=KR?)",
        "(?<=a(?:b|cd))e",
        "(?<=(?:a){2})b",
        "(a)(?<=\\1)b",
        "(",
        ")",
        "[a",
        "a**",
        "^*",
        "a{3,2}",
        "(a)\\2",
        "(?P<n>a)",
        "[z-a]",
        "[[:foo:]]",
        "(?=)",
        "(?)",
        "\\",
    ] {
        match BoostRegex::new(pattern) {
            Err(Error::InvalidValue(_)) => {}
            other => panic!("{pattern:?}: expected InvalidValue, got {other:?}"),
        }
    }
}

#[test]
fn refuses_untranslated_syntax() {
    for pattern in [
        "a++",
        "a{2}+",
        "\\Qa.b\\E",
        "(a)(?1)",
        "(?R)",
        "(a)(?(1)b|c)",
        "(?|(a)|(b))",
        "a\\Kb",
        "\\Ga",
        "a\\Z",
        "(?x)a b",
        "\\p{L}",
        "\\012",
        "\u{e9}",
        "[\u{e9}]",
        "\\x80",
        "(*FAIL)",
        "[[.a.]]",
        "(a)\\10",
        "(?>a+)b",
        "(?!a*b)a",
        "(?<=a{256})b",
    ] {
        match BoostRegex::new(pattern) {
            Err(Error::Unsupported(_)) => {}
            other => panic!("{pattern:?}: expected Unsupported, got {other:?}"),
        }
    }
}

#[test]
fn limits_are_errors() {
    let long = "a".repeat(MAX_PATTERN_BYTES + 1);
    assert!(matches!(
        BoostRegex::new(&long),
        Err(Error::InvalidValue(_))
    ));
    let deep = format!(
        "{}a{}",
        "(".repeat(MAX_GROUP_DEPTH + 1),
        ")".repeat(MAX_GROUP_DEPTH + 1)
    );
    assert!(matches!(BoostRegex::new(&deep), Err(Error::Unsupported(_))));
    let fine = format!(
        "{}a{}",
        "(".repeat(MAX_GROUP_DEPTH),
        ")".repeat(MAX_GROUP_DEPTH)
    );
    assert_eq!(
        BoostRegex::new(&fine).unwrap().mark_count(),
        MAX_GROUP_DEPTH
    );
    let zero = RegexOptions {
        backtrack_limit: 0,
        ..RegexOptions::default()
    };
    assert!(matches!(
        BoostRegex::with_options("a", zero),
        Err(Error::InvalidValue(_))
    ));

    // A lookbehind of MAX_LOOKBEHIND_WIDTH bytes is translated, one byte wider is not.
    let widest = BoostRegex::new(&format!("(?<=a{{{MAX_LOOKBEHIND_WIDTH}}})b")).unwrap();
    let mut haystack = vec![b'a'; MAX_LOOKBEHIND_WIDTH];
    haystack.push(b'b');
    assert_eq!(
        widest.search(&haystack).unwrap().unwrap().range(),
        MAX_LOOKBEHIND_WIDTH..MAX_LOOKBEHIND_WIDTH + 1
    );
    assert_eq!(widest.search(&haystack[1..]).unwrap(), None);
    let wider = format!("(?<=a{{{}}})b", MAX_LOOKBEHIND_WIDTH + 1);
    assert_eq!(
        refused(&wider, BoostRegex::new(&wider)),
        "a lookbehind wider than MAX_LOOKBEHIND_WIDTH bytes"
    );

    // A case-insensitive letter translates to the 9 bytes `(?i:\x61)`, and the
    // counted spelling adds 12 around the whole, `(?:` and `)(?=)(?=)`: the largest
    // such pattern whose engine pattern fits MAX_TRANSLATED_BYTES compiles, one
    // letter more does not.
    let icase = RegexOptions {
        icase: true,
        ..RegexOptions::default()
    };
    let fitting = (MAX_TRANSLATED_BYTES - 12) / 9;
    assert!(BoostRegex::with_options(&"a".repeat(fitting), icase).is_ok());
    let too_long = "a".repeat(fitting + 1);
    assert_eq!(
        refused("icase letters", BoostRegex::with_options(&too_long, icase)),
        "more than MAX_TRANSLATED_BYTES bytes of engine pattern"
    );

    let regex = BoostRegex::new("(a)(b)").unwrap();
    assert!(matches!(
        regex.tokens(b"ab", &[]),
        Err(Error::InvalidValue(_))
    ));
    // Boost's match_results::operator[]: -2 is the suffix, and an index past
    // the groups is the null sub-match at the end of the match.
    let tokens: Vec<SubMatch> = regex
        .tokens(b"xabyab", &[-2, 3, -1])
        .unwrap()
        .map(Result::unwrap)
        .collect();
    assert_eq!(
        tokens,
        [
            SubMatch {
                range: 3..6,
                matched: true
            },
            SubMatch {
                range: 3..3,
                matched: false
            },
            SubMatch {
                range: 0..1,
                matched: true
            },
            SubMatch {
                range: 6..6,
                matched: false
            },
            SubMatch {
                range: 6..6,
                matched: false
            },
            SubMatch {
                range: 3..4,
                matched: true
            },
        ]
    );
}

/// `fancy-regex`'s analyzer multiplies a repeated sub-expression's shortest match
/// by the repeat's minimum without an overflow check. In this test profile, with
/// overflow checks on, constructing these patterns panicked inside the engine; the
/// facade now refuses them before the engine sees them.
#[test]
fn nested_repeat_bounds_do_not_overflow_the_engine() {
    // The review calls this the 43-byte pattern; it is 42 bytes.
    let nested = "(?:(?:a{999999999}){999999999}){999999999}";
    assert_eq!(nested.len(), 42);
    for pattern in [
        nested,
        "(?:(?:(?=a)a{999999999}){999999999}){999999999}",
        r"(?:(?:\bb{999999999}){999999999}){999}",
        "(?:ab){999999999}",
        r"(a{999999999})\1{999999999}",
    ] {
        let category = within(60, pattern, move || {
            refused(pattern, BoostRegex::new(pattern))
        });
        assert_eq!(category, LONG_SHORTEST_MATCH, "{pattern:?}");
    }
    // Below the bound the engine gets the pattern, and needs 999^3 bytes to match.
    let below = BoostRegex::new(r"(?:(?:\ba{999}){999}){999}").unwrap();
    assert_eq!(below.search(b"aaa").unwrap(), None);
}

/// The first review's hang: a bounded repeat of a group that can match the empty
/// string ran every one of its empty iterations without spending the backtracking
/// budget (`(?:a{0}\b){999999999}` took 16.9 s on one byte), and gave answers Boost
/// does not, because Boost ends a repeat after an empty iteration. The second review
/// found that the unbounded forms differ from Boost as well: the engine fails an
/// empty iteration and backtracks into the group's other alternatives, where Boost
/// accepts it and leaves the repeat. Each is refused at once, nested forms included.
#[test]
fn nullable_group_repeats_are_refused_within_a_time_bound() {
    for pattern in [
        r"(?:a{0}\b){999999999}",
        r"(?:a{0}\b){999999999,}",
        "(?:(?=a)(?=a)){99999999}",
        r"(?:(?:a{0}\b){999999999}){999999999}",
        "(?:(?:(?=a)(?=a)){999999999}){999999999}",
        "(?:(?:a{0}(?=a)){99999}){99999}",
        r"(?:a{0}\b){0,999999999}",
        "(?:a{0}$){99999999}",
        r"(?:\b|a){2}b",
        "(?:(?=a)|a){2}",
        r"(?:\b|a){3}",
        "(?:^|a){2,}b",
        "(?:(?=a)|a){2,}",
        // The second review: Boost gives 0..1 on "ba", 0..1 on "ba", 0..2 on "aba"
        // and 0..2 on "abab", where the engine gave 0..2, 0..2, 0..3 and 0..4.
        "(?:b?|a)*",
        "(?:a*?b?)*",
        "(?:(?:ab)?|a)*",
        "(?:(?:a|b?)*?)+b",
        "(?:b?|a)+",
        "(?:b?|a){1,}",
        "(?:b?|a)*?",
        r"(?:a{0}\b)+a",
        r"(?:\b|a)+?b",
        "(?:(?=a)|a)+?b",
        "(?:(?:(?:b?|a)*)*)*",
    ] {
        let category = within(60, pattern, move || {
            refused(pattern, BoostRegex::new(pattern))
        });
        assert_eq!(category, NULLABLE_GROUP, "{pattern:?}");
    }
    // At most one iteration, a group that cannot match empty, or a backreference to
    // a group that is not open, which matches the same text in every iteration:
    // compiled, with Boost's answers (transcribed from the oracle output for the
    // `ADV` family).
    // (pattern, haystack, Boost's groups of the match)
    type Case<'a> = (&'a str, &'a [u8], Option<Vec<Option<Range<usize>>>>);
    let cases: &[Case] = &[
        ("(?:b?|a)?", b"ba", Some(vec![Some(0..1)])),
        ("(?:b|a)*", b"bab", Some(vec![Some(0..3)])),
        (r"(a?)\1+", b"aba", Some(vec![Some(0..0), Some(0..0)])),
        (r"(b?)\1*a", b"bab", Some(vec![Some(0..2), Some(0..1)])),
    ];
    for (pattern, haystack, expected) in cases {
        let found = BoostRegex::new(pattern)
            .unwrap()
            .search(haystack)
            .unwrap()
            .map(|captures| {
                (0..captures.len())
                    .map(|group| captures.get(group))
                    .collect()
            });
        assert_eq!(&found, expected, "{pattern:?} on {haystack:?}");
    }
    // Over a long haystack each of them answers or reports a limit, within the bound.
    let haystack: &'static [u8] = Box::leak(vec![b'a'; 16 * 1024].into_boxed_slice());
    for pattern in ["(?:a?)?b", r"(a?)\1+b", "(?:b|a)+?c", r"(b?)\1*c"] {
        let outcome = within(120, pattern, move || {
            let regex = BoostRegex::new(pattern).unwrap();
            regex
                .search(haystack)
                .map(|found| found.map(|captures| captures.range()))
        });
        assert!(
            matches!(outcome, Ok(None) | Err(Error::InvalidValue(_))),
            "{pattern:?}: {outcome:?}"
        );
    }
}

/// Boost marks a repeat of a one-byte atom as leading when only group starts and
/// ends, flag groups that keep case sensitivity, `^ $ \b \B \< \> \A \z` and whole
/// lookarounds precede it, in an expression without backreferences
/// (`basic_regex_creator::probe_leading_repeat`). When a lazy leading repeat with a
/// finite maximum extends, Boost records where it got to, and a failed start
/// position resumes the search behind that point (`perl_matcher::unwind_*_repeat`,
/// `match_prefix`), skipping start positions that match. The facade refuses such
/// a repeat; the others are compiled and give Boost's answers (transcribed from the
/// oracle output for the `ADV` family).
#[test]
fn leading_lazy_repeats_are_refused() {
    let icase = RegexOptions {
        icase: true,
        ..RegexOptions::default()
    };
    let plain = RegexOptions::default();
    // The second review's examples: Boost answers 3..4, 3..4, no match, no match and
    // no match, where a search that tries every start position finds 1..4, 1..4,
    // 1..6, 2..7 and 1..6.
    for (options, pattern) in [
        (plain, r"a{1,3}?\b"),
        (plain, ".{1,3}?(?=b)"),
        (plain, r".{3,5}?\b"),
        (plain, r"\w{3,5}?\b"),
        (icase, "[ab]{3,5}?(?!a)"),
        (plain, r"a{0,2}?\b"),
        (plain, r"^a{1,3}?\b"),
        (plain, r"(?=a)(?!x)()(?:(a{1,3}?))(?#c)(?-i)\b"),
        (icase, r"a{1,3}?\b"),
        (plain, "[^b]{0,2}?b"),
    ] {
        let category = within(60, pattern, move || {
            refused(pattern, BoostRegex::with_options(pattern, options))
        });
        assert_eq!(
            category, LEADING_LAZY,
            "{pattern:?} icase={}",
            options.icase
        );
    }
    // (options, pattern, haystack, Boost's match)
    type Case<'a> = (RegexOptions, &'a str, &'a [u8], Option<Range<usize>>);
    let cases: &[Case] = &[
        // An alternation state in front of the repeat.
        (plain, r"x|a{1,3}?\b", b"aaaa", Some(1..4)),
        (plain, r"a{1,3}?\b|x", b"aaaa", Some(1..4)),
        // A toggle_case state in front of it.
        (plain, r"(?i)a{1,3}?\b", b"aaaa", Some(1..4)),
        (plain, r"(?i:)a{1,3}?\b", b"aaaa", Some(1..4)),
        (icase, r"(?-i)a{1,3}?\b", b"aaaa", Some(1..4)),
        // A repeat state in front of it.
        (plain, r"(?=x)?a{1,3}?\b", b"aaaa", Some(1..4)),
        (plain, r"(?:a{1,3}?)+\b", b"aaaa", Some(0..4)),
        // No room for two iterations above the minimum, or no maximum.
        (plain, r"a{1,2}?\b", b"aaaa", Some(2..4)),
        (plain, r"a{1,}?\b", b"aaaa", Some(0..4)),
        // A backreference anywhere in the expression.
        (plain, r"a{1,3}?\b(a)\1", b"aaaa", None),
    ];
    for &(options, pattern, haystack, ref expected) in cases {
        let found = BoostRegex::with_options(pattern, options)
            .unwrap_or_else(|error| panic!("{pattern:?}: {error}"))
            .search(haystack)
            .unwrap()
            .map(|captures| captures.range());
        assert_eq!(&found, expected, "{pattern:?} on {haystack:?}");
    }
}

/// The third review: Boost saves a group's span when the group starts and restores
/// it only when the match backtracks past that start, so a backreference inside the
/// group compares against what an abandoned alternative or iteration left
/// (`perl_matcher::match_startmark`). Boost matches all of `abx` with `(a|\1b)x`,
/// `aab` with `(a|\1a)b` and `abax` with `(a)(b|\2a)x`, where the engine found no
/// match, `1..3` and no match; on `(?:(a|\1b)c)+` and `acbc` the engine panicked.
/// Each is refused within a watchdog. Backreferences to closed, later, unset,
/// repeated, duplicate-named, zero-repeated and case-insensitive groups, and into
/// lookarounds and atomic groups, are compiled and give Boost's answers
/// (transcribed from the oracle output for the `ADV` family).
#[test]
fn backreferences_inside_their_group_are_refused() {
    let icase = RegexOptions {
        icase: true,
        ..RegexOptions::default()
    };
    let plain = RegexOptions::default();
    for (options, pattern) in [
        (plain, r"(a|\1b)x"),
        (plain, r"(a|\1a)b"),
        (plain, r"(a)(b|\2a)x"),
        (plain, r"(ab|\1b)x"),
        (plain, r"(x|a|\1\1)b"),
        (plain, r"(a|\1+)b"),
        (plain, r"(a|\1*)$"),
        (plain, r"(a?\1*)\b"),
        (plain, r"(?:(a|\1b)c)+"),
        (plain, r"(?:(a|b\1)c)+"),
        (plain, r"(a|b\1)+"),
        (plain, r"((a)|\1b)x"),
        (plain, r"(?<n>a|\1b)x"),
        (plain, r"(a(?=\1))b"),
        (plain, r"(a\1?)+"),
        (plain, r"(a\1)"),
        (icase, r"(a|\1b)x"),
    ] {
        let category = within(60, pattern, move || {
            refused(pattern, BoostRegex::with_options(pattern, options))
        });
        assert_eq!(
            category, OPEN_GROUP_BACKREFERENCE,
            "{pattern:?} icase={}",
            options.icase
        );
    }
    // (options, pattern, haystack, Boost's groups of the match)
    type Case<'a> = (
        RegexOptions,
        &'a str,
        &'a [u8],
        Option<Vec<Option<Range<usize>>>>,
    );
    let cases: &[Case] = &[
        // Group 2 is closed inside the open group 1, and an abandoned branch
        // restores it.
        (
            plain,
            r"((a)|\2b)x",
            b"abax",
            Some(vec![Some(2..4), Some(2..3), Some(2..3)]),
        ),
        (plain, r"((a)|\2b)x", b"abx", None),
        (
            plain,
            r"(a)(?:x|\1b)",
            b"aab",
            Some(vec![Some(0..3), Some(0..1)]),
        ),
        (
            plain,
            r"(a)(b|\1a)x",
            b"abx",
            Some(vec![Some(0..3), Some(0..1), Some(1..2)]),
        ),
        (plain, r"(a)(b|\1a)x", b"abax", None),
        // Unset in the current alternative, later in the pattern, or repeated.
        (
            plain,
            r"(?:(a)|b)\1",
            b"aab",
            Some(vec![Some(0..2), Some(0..1)]),
        ),
        (
            plain,
            r"(?:\1b|(a))+",
            b"aaab",
            Some(vec![Some(0..4), Some(1..2)]),
        ),
        (
            plain,
            r"(?:(a)|b\1)+",
            b"abab",
            Some(vec![Some(0..3), Some(0..1)]),
        ),
        (
            plain,
            r"(a|b)+\1",
            b"abb",
            Some(vec![Some(0..3), Some(1..2)]),
        ),
        (
            plain,
            r"(a|ab)*c\1",
            b"acacbc",
            Some(vec![Some(0..3), Some(0..1)]),
        ),
        (plain, r"\1(a)", b"aab", None),
        (plain, r"(a){0}\1", b"aab", None),
        (
            plain,
            r"(a?)\1*b",
            b"acbc",
            Some(vec![Some(2..3), Some(2..2)]),
        ),
        // Into lookarounds and atomic groups.
        (
            plain,
            r"(a)(?!\1)\w",
            b"aab",
            Some(vec![Some(1..3), Some(1..2)]),
        ),
        (
            plain,
            r"(a)(?>\1)b",
            b"aaab",
            Some(vec![Some(1..4), Some(1..2)]),
        ),
        (
            plain,
            r"(a)(?<=(?=\1).)",
            b"bab",
            Some(vec![Some(1..2), Some(1..2)]),
        ),
        // Duplicate names.
        (
            plain,
            r"(?:(?<n>a)|(?<n>b))\1",
            b"aab",
            Some(vec![Some(0..2), Some(0..1), None]),
        ),
        // Case-insensitive, with bytes above 0x7F, which have no case.
        (icase, r"(.)\1", b"aA", Some(vec![Some(0..2), Some(0..1)])),
        (icase, r"(.)\1", b"\xc3\xa9\xc3\x89", None),
        (icase, r"(.)\1", b"\xe9\xc9", None),
        (
            plain,
            r"(a)(?i)\1",
            b"aA",
            Some(vec![Some(0..2), Some(0..1)]),
        ),
        (plain, r"(a)(?i)\1", b"Aa", None),
        (plain, r"(?i:(a))\1", b"aA", None),
        (
            icase,
            r"(.+)\1",
            b"abab",
            Some(vec![Some(0..4), Some(0..2)]),
        ),
    ];
    for (options, pattern, haystack, expected) in cases {
        let found = BoostRegex::with_options(pattern, *options)
            .unwrap_or_else(|error| panic!("{pattern:?}: {error}"))
            .search(haystack)
            .unwrap()
            .map(|captures| {
                (0..captures.len())
                    .map(|group| captures.get(group))
                    .collect()
            });
        assert_eq!(&found, expected, "{pattern:?} on {haystack:?}");
    }
    let duplicate = BoostRegex::new(r"(?:(?<n>a)|(?<n>b))\1").unwrap();
    let captures = duplicate.search(b"aab").unwrap().unwrap();
    assert_eq!(captures.name("n"), Some(0..1));
}

/// The third round's fuzzing found rewrites inside the engine crates that change
/// answers. `fancy-regex`'s optimizer turned `X+Y?X+` into `X+(?:YX+)?` (so
/// `\w{1,}b?\w{1,}` matched `a`), `(X+)+` into `(X+)` under a backreference, and
/// `(X)*` into `(X)?` (which captures differently when `X` is lazy); `regex-syntax`
/// factored the common prefix out of `XA|XB`, which loses leftmost-first priority
/// when `X` can match in several ways (`[ab]+b|[ab]+c` matched `0..4` of `abbc`).
/// The facade now spells these expressions so that neither rewrite applies, and
/// gives Boost's answers (transcribed from the oracle output for the `ADV` family).
#[test]
fn engine_rewrites_keep_boost_answers() {
    // (pattern, whole-haystack match instead of search, haystack, Boost's groups)
    type Case<'a> = (&'a str, bool, &'a [u8], Option<Vec<Option<Range<usize>>>>);
    let cases: &[Case] = &[
        (r"\w{1,}b?\w{1,}", false, b"a", None),
        ("a+b?a+", false, b"aab", Some(vec![Some(0..2)])),
        ("a+b?a+", true, b"aab", None),
        ("a+b*a+", true, b"abba", Some(vec![Some(0..4)])),
        (
            "(a)+b?(a)+",
            false,
            b"aab",
            Some(vec![Some(0..2), Some(0..1), Some(1..2)]),
        ),
        (
            r"(.{1,})+\1+",
            false,
            b"abb",
            Some(vec![Some(0..3), Some(1..2)]),
        ),
        (
            r"(.{1,})+\1+",
            true,
            b"abab",
            Some(vec![Some(0..4), Some(0..2)]),
        ),
        (
            r"(\w+?)*?(?<=b)",
            true,
            b"ab",
            Some(vec![Some(0..2), Some(1..2)]),
        ),
        (
            r"(\w+?)*?(?<=b)",
            true,
            b"abab",
            Some(vec![Some(0..4), Some(3..4)]),
        ),
        (
            r"|(a+?)*?\w|\w()",
            true,
            b"aab",
            Some(vec![Some(0..3), Some(1..2), None]),
        ),
        ("(a+?)*b", false, b"aab", Some(vec![Some(0..3), Some(1..2)])),
        (
            "(a{2,}?)*b",
            false,
            b"aab",
            Some(vec![Some(0..3), Some(0..2)]),
        ),
        ("(?:a+(?:ba*)?)+$", false, b"abba", Some(vec![Some(3..4)])),
        ("(?:a+(?:ba*)?)+$", true, b"abba", None),
        ("(?:a+b?a*)+$", false, b"abba", Some(vec![Some(3..4)])),
        (
            r"\w+?a{1,3}?|\w{1,}?()",
            false,
            b"aba",
            Some(vec![Some(0..3), None]),
        ),
        (
            r"\w+?a{1,3}?|\w{1,}?()",
            false,
            b"abba",
            Some(vec![Some(0..4), None]),
        ),
        ("[ab]+b|[ab]+c", false, b"abbc", Some(vec![Some(0..3)])),
        ("[ab]+b|[ab]+c", true, b"abbc", Some(vec![Some(0..4)])),
        (
            "(?:[ab]+?a|[ab]+?())c",
            false,
            b"abbc",
            Some(vec![Some(0..4), Some(3..3)]),
        ),
        ("(?:a|ab)b|(?:a|ab)c", true, b"abb", Some(vec![Some(0..3)])),
        ("a*b|a*()", false, b"a", Some(vec![Some(0..1), Some(1..1)])),
        (
            "(a|b)*c|(a|b)*b",
            false,
            b"abbc",
            Some(vec![Some(0..4), Some(2..3), None]),
        ),
    ];
    for (pattern, whole, haystack, expected) in cases {
        let regex = BoostRegex::new(pattern).unwrap_or_else(|error| panic!("{pattern:?}: {error}"));
        let found = if *whole {
            regex.full_match(haystack)
        } else {
            regex.search(haystack)
        }
        .unwrap()
        .map(|captures| {
            (0..captures.len())
                .map(|group| captures.get(group))
                .collect()
        });
        assert_eq!(
            &found, expected,
            "{pattern:?} whole={whole} on {haystack:?}"
        );
    }
}

/// The fourth review: `fancy-regex` hands the body of an atomic group that needs no
/// backtracking (or a branch or trailing run of it) to `regex-automata`, where
/// `regex-syntax` factored a common prefix out of an alternation, so the group kept a
/// different match than Boost's first one. At the previous commit
/// `(?>[ab]?b|[ab]?c)` found `0..2` in `bc`, `(?>(?:a|ab)c|(?:a|ab)b)` did not match
/// all of `abc`, `(?=(?>(?:a?|ab)c|(?:a?|ab)b)c)` matched at 0 in `abc`, and
/// `(?!(?>[ab]?b|[ab]?c)c)\w` found `0..1` in `bc`. Every alternation inside an atomic
/// group is now guarded, except where it is part of a lookbehind's width, and the
/// answers are Boost's (transcribed from the oracle output for the `ADV` family).
#[test]
fn atomic_alternations_keep_boost_answers() {
    let icase = RegexOptions {
        icase: true,
        ..RegexOptions::default()
    };
    let plain = RegexOptions::default();
    // (options, pattern, whole-haystack match instead of search, haystack, Boost's groups)
    type Case<'a> = (
        RegexOptions,
        &'a str,
        bool,
        &'a [u8],
        Option<Vec<Option<Range<usize>>>>,
    );
    let cases: &[Case] = &[
        (
            plain,
            "(?>[ab]?b|[ab]?c)",
            false,
            b"bc",
            Some(vec![Some(0..1)]),
        ),
        (
            plain,
            "(?>[ab]?b|[ab]?c)",
            false,
            b"\xc3\xa9bc",
            Some(vec![Some(2..3)]),
        ),
        (
            plain,
            "(?>(?:a|ab)c|(?:a|ab)b)",
            false,
            b"abc",
            Some(vec![Some(0..3)]),
        ),
        (
            plain,
            "(?>(?:a|ab)c|(?:a|ab)b)",
            true,
            b"abc",
            Some(vec![Some(0..3)]),
        ),
        (
            icase,
            "(?>(?:a|ab)c|(?:a|ab)b)",
            false,
            b"Abc",
            Some(vec![Some(0..3)]),
        ),
        (icase, "(?>(?:a|ab)c|(?:a|ab)b)", true, b"aAbc", None),
        (
            plain,
            "(?>(?:a|ab)c|(?:a|ab)b)c?",
            false,
            b"abcc",
            Some(vec![Some(0..4)]),
        ),
        (
            plain,
            "(?=(?>(?:a?|ab)c|(?:a?|ab)b)c)",
            false,
            b"abc",
            Some(vec![Some(1..1)]),
        ),
        (
            plain,
            r"\b(?>(?:a|ab)c|(?:a|ab)b)",
            false,
            b"abc",
            Some(vec![Some(0..3)]),
        ),
        (
            plain,
            "(a)(?>(?:a|ab)c|(?:a|ab)b)",
            false,
            b"aabc",
            Some(vec![Some(0..4), Some(0..1)]),
        ),
        (
            plain,
            r"(?>\b(?:[ab]?b|[ab]?c))",
            false,
            b"bc",
            Some(vec![Some(0..1)]),
        ),
        (
            plain,
            r"(?>(?:[ab]?b|[ab]?c)\b)",
            false,
            b"abc",
            Some(vec![Some(1..3)]),
        ),
        (
            plain,
            r"(?!(?>[ab]?b|[ab]?c)c)\w",
            false,
            b"bc",
            Some(vec![Some(1..2)]),
        ),
        (
            plain,
            "(?>a??b|a??c|a??)",
            false,
            b"ac",
            Some(vec![Some(0..2)]),
        ),
        // A lookahead inside a lookbehind is guarded; an alternation of one width in
        // a lookbehind is not, and needs no guard.
        (
            plain,
            "(?<=(?=(?>[ab]?b|[ab]?c)c)..)",
            false,
            b"bcc",
            Some(vec![Some(2..2)]),
        ),
        (
            plain,
            "(?<=(?=(?>(?:a|ab)c|(?:a|ab)b)).)c",
            false,
            b"abc",
            None,
        ),
        (
            plain,
            "(?<=(?>ab|a[bc]))c",
            false,
            b"abc",
            Some(vec![Some(2..3)]),
        ),
        (plain, "(?<!(?>ab|a.))c", false, b"abc", None),
        // Nested and repeated atomic groups, and backreferences beside them.
        (
            plain,
            "(?>(?>[ab]?b|[ab]?c)|b)c",
            false,
            b"bcc",
            Some(vec![Some(0..2)]),
        ),
        (
            plain,
            "(?:x|(?>[ab]?b|[ab]?c))+c",
            false,
            b"abbc",
            Some(vec![Some(0..4)]),
        ),
        (
            plain,
            "(?>(?:a|ab)c|(?:a|ab)b){1,2}c",
            false,
            b"abcc",
            Some(vec![Some(0..4)]),
        ),
        (
            plain,
            r"(a|ab)(?>\1?c|\1?b)",
            false,
            b"abc",
            Some(vec![Some(0..2), Some(0..1)]),
        ),
        (
            plain,
            r"(.)(?i)(?>\1?[ab]?b|\1?[ab]?c)",
            false,
            b"aAbc",
            Some(vec![Some(0..3), Some(0..1)]),
        ),
        (plain, r"(?>[ab]?b|[ab]?c)(a)?\1", false, b"bc", None),
    ];
    for (options, pattern, whole, haystack, expected) in cases {
        let regex = BoostRegex::with_options(pattern, *options)
            .unwrap_or_else(|error| panic!("{pattern:?}: {error}"));
        let found = if *whole {
            regex.full_match(haystack)
        } else {
            regex.search(haystack)
        }
        .unwrap()
        .map(|captures| {
            (0..captures.len())
                .map(|group| captures.get(group))
                .collect()
        });
        assert_eq!(
            &found, expected,
            "{pattern:?} icase={} whole={whole} on {haystack:?}",
            options.icase
        );
    }
    // Boost's tokens -1 and 0 over `abcc`: `ab`, then `c` twice, each after an
    // empty unmatched gap.
    let tokens: Vec<SubMatch> = BoostRegex::new("(?>[ab]?b|[ab]?c)")
        .unwrap()
        .tokens(b"abcc", &[-1, 0])
        .unwrap()
        .map(Result::unwrap)
        .collect();
    let token = |range: Range<usize>, matched: bool| SubMatch { range, matched };
    assert_eq!(
        tokens,
        [
            token(0..0, false),
            token(0..2, true),
            token(2..2, false),
            token(2..3, true),
            token(3..3, false),
            token(3..4, true),
        ]
    );
}

/// The fifth review: before matching, Boost builds a map of the bytes that can start
/// each way on of every repeat and alternation (`basic_regex_creator::create_startmap`),
/// and the map is wrong in two ways the engine does not share. When the walk that
/// builds it recurses, at `\<`, `\>` or a repeat it loops back to, it restarts from
/// the case sensitivity of the state whose map it builds, so after a case switch a
/// letter is looked up with the wrong case: `(?i)\<A` does not match `A`, `(?i:b.+)*C`
/// finds `2..3` in `bcC` and `(?i:bc+)+C` nothing, where the engine found `0..1`,
/// `0..3` and `0..3`. And a recursion at `\<` or `\>` removes bytes from the whole
/// map it fills, which after a loop back already holds the bytes that start another
/// iteration: `(?:\w\w+){2}\>` does not match `aaaa`, and `(\w+?)+\>` captures group
/// 1 at `0..4`, where the engine found `0..4` and `3..4`. The walk also stops with
/// `error_complexity` beyond 100 nested recursions, which `$` chains, alternations in
/// repeated groups and nested repeated groups reach. Each of these is refused within
/// a watchdog; the neighbours are compiled and give Boost's answers (transcribed from
/// the oracle output for the `ADV` family).
#[test]
fn start_map_shapes_are_refused() {
    const CASE_SWITCHED_REPEAT: &str = "a repeat of a group holding a repeat or alternation under other case sensitivity (Boost builds its start map with the wrong case)";
    const WORD_START_CASE_SWITCH: &str = "\\< in an expression that switches case sensitivity (Boost builds its start map with the wrong case)";
    const WORD_BOUNDARY_AFTER_LOOP: &str = "\\< or \\> in an expression with a repeated group holding a repeat or alternation (Boost's start map drops the bytes that start another iteration)";
    const START_MAP_RECURSION: &str = "more line anchors, word-boundary assertions, alternations and repeated groups than Boost's start-map recursion limit allows (Boost throws error_complexity)";
    let icase = RegexOptions {
        icase: true,
        ..RegexOptions::default()
    };
    let plain = RegexOptions::default();
    let alternatives = |count: usize| {
        "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"[..count]
            .chars()
            .map(String::from)
            .collect::<Vec<_>>()
            .join("|")
    };
    let refusals: Vec<(RegexOptions, String, &str)> = vec![
        (plain, r"(?i:b.+)*C".into(), CASE_SWITCHED_REPEAT),
        (plain, r"(?i:b.+){1}C".into(), CASE_SWITCHED_REPEAT),
        (plain, r"(?i:bc+)+C".into(), CASE_SWITCHED_REPEAT),
        (plain, r"(?i)(?-i:b.+)?[a-z]".into(), CASE_SWITCHED_REPEAT),
        (plain, r"(?:y|(?i:b.+))+C".into(), CASE_SWITCHED_REPEAT),
        (plain, r"(?i)(?:a(?-i:b.+))+c".into(), CASE_SWITCHED_REPEAT),
        (plain, r"(?:(?i)x|y+)*C".into(), CASE_SWITCHED_REPEAT),
        (plain, r"(?:(?i:a|b+)c)+".into(), CASE_SWITCHED_REPEAT),
        (plain, r"(?i)\<A".into(), WORD_START_CASE_SWITCH),
        (plain, r"(?i)\<[A-Z]".into(), WORD_START_CASE_SWITCH),
        (plain, r"(?i)\<A|b".into(), WORD_START_CASE_SWITCH),
        (icase, r"(?-i)(?:\<)A".into(), WORD_START_CASE_SWITCH),
        (plain, r"\<a(?i)b".into(), WORD_START_CASE_SWITCH),
        (plain, r"(?:\w\w+){2}\>".into(), WORD_BOUNDARY_AFTER_LOOP),
        (plain, r"(?:\w{2,}){2,}\>".into(), WORD_BOUNDARY_AFTER_LOOP),
        (plain, r"(?:\W\W??|a){2}\>".into(), WORD_BOUNDARY_AFTER_LOOP),
        (plain, r"(?:\W\W??|a)+\>".into(), WORD_BOUNDARY_AFTER_LOOP),
        (plain, r"(\w+?)+\>".into(), WORD_BOUNDARY_AFTER_LOOP),
        (plain, r"(?:a|ab)+\<".into(), WORD_BOUNDARY_AFTER_LOOP),
        (plain, r"\<(?:\w\w*)+".into(), WORD_BOUNDARY_AFTER_LOOP),
        (plain, r"(?:\.|\<b+)*c".into(), WORD_BOUNDARY_AFTER_LOOP),
        // Boost: error_complexity.
        (plain, format!("{}a", "$".repeat(101)), START_MAP_RECURSION),
        (
            plain,
            format!("(?:a+)*{}", "$".repeat(99)),
            START_MAP_RECURSION,
        ),
        (
            plain,
            format!("(?:(?:{})x?)+", alternatives(51)),
            START_MAP_RECURSION,
        ),
        (
            plain,
            format!("{}{}\\z", "$".repeat(80), r"(?:\b|\B)".repeat(12)),
            START_MAP_RECURSION,
        ),
        // Boost compiles these; the bound counts more than Boost's walk reaches.
        (plain, format!("a{}", "$".repeat(101)), START_MAP_RECURSION),
        (plain, format!("x*{}", "$".repeat(99)), START_MAP_RECURSION),
        (
            plain,
            format!("{}a+{}", "(?:".repeat(30), "$)+".repeat(30)),
            START_MAP_RECURSION,
        ),
    ];
    for (options, pattern, expected) in refusals {
        let shown = pattern.clone();
        let category = within(60, &shown, move || {
            refused(&pattern, BoostRegex::with_options(&pattern, options))
        });
        assert_eq!(category, expected, "{shown:?} icase={}", options.icase);
    }
    // (options, pattern, haystack, Boost's groups of the first match)
    type Case = (
        RegexOptions,
        String,
        &'static [u8],
        Option<Vec<Option<Range<usize>>>>,
    );
    let cases: &[Case] = &[
        // The same groups without a repeat, or with the switch where Boost's walk
        // keeps it.
        (plain, r"(?i:b.+)C".into(), b"bcC", Some(vec![Some(0..3)])),
        (plain, r"(?i:bc)+C".into(), b"aBcC", Some(vec![Some(1..4)])),
        (
            icase,
            r"(?i:b(?:c|d))+C".into(),
            b"bCc",
            Some(vec![Some(0..3)]),
        ),
        (
            plain,
            r"(?:b.+(?i))*C".into(),
            b"bCc",
            Some(vec![Some(1..2)]),
        ),
        (
            plain,
            r"(?:b.+(?i))*C".into(),
            b"aBcC",
            Some(vec![Some(3..4)]),
        ),
        (
            plain,
            r"(?i)(?:b.+)*C".into(),
            b"bc",
            Some(vec![Some(1..2)]),
        ),
        (
            plain,
            r"(?:(?i:a|b)c+)+C".into(),
            b"BcC",
            Some(vec![Some(0..3)]),
        ),
        (
            plain,
            r"(?:(?i)a|b)*C".into(),
            b"bCc",
            Some(vec![Some(0..2)]),
        ),
        (icase, r"(?i)\<A".into(), b" A", Some(vec![Some(1..2)])),
        (icase, r"\<A".into(), b"\xe9A", Some(vec![Some(1..2)])),
        (plain, r"(?i)a\>".into(), b"aa aa", Some(vec![Some(1..2)])),
        (plain, r"(?i)\bA".into(), b"..a", Some(vec![Some(2..3)])),
        // `\b`, `\B` and `$` after the same groups, and the groups without a repeat
        // or without an inner repeat.
        (
            plain,
            r"(?:\w\w+){2}\b".into(),
            b"abab",
            Some(vec![Some(0..4)]),
        ),
        (
            plain,
            r"(?:\W\W??|a)+\B".into(),
            b"aaaa",
            Some(vec![Some(0..3)]),
        ),
        (
            plain,
            r"(\w+?)+$".into(),
            b"aaaa",
            Some(vec![Some(0..4), Some(3..4)]),
        ),
        (
            plain,
            r"(\w+?)\>".into(),
            b"ab ab",
            Some(vec![Some(0..2), Some(0..2)]),
        ),
        (
            plain,
            r"(?:\w\w){2}\>".into(),
            b"abab",
            Some(vec![Some(0..4)]),
        ),
        (plain, r"(?:ab)+\>".into(), b"ab ab", Some(vec![Some(0..2)])),
        (plain, r"\<\w+\>".into(), b"..a", Some(vec![Some(2..3)])),
        // Below the recursion limit.
        (plain, format!("{}a", "$".repeat(98)), b"a", None),
        (
            plain,
            format!("(?:a+)*{}", "$".repeat(96)),
            b"aaaa",
            Some(vec![Some(0..4)]),
        ),
        (
            plain,
            format!("(?:(?:{})x?)+", alternatives(48)),
            b"abab",
            Some(vec![Some(0..4)]),
        ),
        (
            plain,
            format!("{}a+{}", "(?:".repeat(20), "$)+".repeat(20)),
            b"aaaa",
            Some(vec![Some(0..4)]),
        ),
        (
            plain,
            format!("{}{}\\z", "$".repeat(60), r"(?:\b|\B)".repeat(12)),
            b"A\n",
            Some(vec![Some(2..2)]),
        ),
    ];
    for (options, pattern, haystack, expected) in cases {
        let found = BoostRegex::with_options(pattern, *options)
            .unwrap_or_else(|error| panic!("{pattern:?}: {error}"))
            .search(haystack)
            .unwrap()
            .map(|captures| {
                (0..captures.len())
                    .map(|group| captures.get(group))
                    .collect()
            });
        assert_eq!(&found, expected, "{pattern:?} on {haystack:?}");
    }
}

/// The fifth review's minor findings. Boost reads the digits of `\x` with
/// `std::istream` (`cpp_regex_traits::toi`), which skips white space and accepts a
/// sign and a `0x` prefix: `\x{+41}`, `\x{ 41}`, `\x{0x41}` and `\x+4` are valid
/// Boost patterns, where the facade reported a syntax error, and `\x0x` is not, where
/// the facade compiled it. Such escapes are refused; the others read as before
/// (compile outcomes and answers from the oracle output). Boost's lookbehind width
/// calculation rejects a lookbehind with more than 1,024 alternations on one path,
/// and a translation nested as deep as `fancy-regex`'s parser refuses is refused
/// before the engine sees it.
#[test]
fn escapes_and_limits_follow_boost() {
    const HEX_STREAM: &str = "a hexadecimal escape whose digits start with a sign, white space or 0x (Boost reads them with std::istream)";
    const LOOKBEHIND_ALTERNATIONS: &str =
        "a lookbehind with more than 1,024 alternations (Boost's backstep calculation gives up)";
    const TOO_DEEP: &str =
        "groups and repeats nested deeper than the engine's parser allows once translated";
    // Boost compiles the first ten and rejects the other four.
    for pattern in [
        r"\x{+41}",
        r"\x{ 41}",
        "\\x{\t41}",
        r"\x{0x41}",
        r"\x{0X41}",
        r"\x{-0}",
        r"\x+4",
        r"\x 4",
        r"[\x 4]",
        r"\x-0",
        r"\x0x",
        r"\x0X",
        r"[\x0x]",
        r"\x{0x}",
    ] {
        assert_eq!(refused(pattern, BoostRegex::new(pattern)), HEX_STREAM);
    }
    // (pattern, haystack, Boost's match)
    for (pattern, haystack, expected) in [
        (r"\x{004}", &b"\x04"[..], Some(0..1)),
        (r"\x04", b"\x04", Some(0..1)),
        (r"[\x4-\x{41}]", b"ab ab", Some(2..3)),
        (r"[\x4-\x{41}]", b"a", None),
        (r"\x4x", b"4x", None),
    ] {
        let regex = BoostRegex::new(pattern).unwrap();
        assert_eq!(
            regex.search(haystack).unwrap().map(|found| found.range()),
            expected,
            "{pattern:?} on {haystack:?}"
        );
    }
    // Boost rejects these (error_badbrace, error_escape), and so does the facade.
    for pattern in [r"\x{4 }", r"\x{,4}", r"\xg", r"\x{}"] {
        assert!(
            matches!(BoostRegex::new(pattern), Err(Error::InvalidValue(_))),
            "{pattern:?}"
        );
    }
    // 1,026 alternations on one path: Boost reports an invalid lookbehind. 1,030 in
    // one alternation: Boost compiles it; refused as well.
    for pattern in [
        format!("(?<={})a", "(?:|)".repeat(1026)),
        format!("(?<!{})a", "|".repeat(1030)),
    ] {
        assert_eq!(
            refused(&pattern, BoostRegex::new(&pattern)),
            LOOKBEHIND_ALTERNATIONS
        );
    }
    assert!(BoostRegex::new(&format!("(?<={})a", "(?:|)".repeat(16))).is_ok());
    // A repeat nested 48 deep in repeats, which the engine's optimizer check spells
    // with a shield at every level, would reach the engine's 64 levels.
    let deep = format!("{}a+{}", "(?:".repeat(48), ")+".repeat(48));
    assert_eq!(refused(&deep, BoostRegex::new(&deep)), TOO_DEEP);
    let too_deep = format!("{}a{{2,}}{}", "(?:".repeat(20), "){2,}".repeat(20));
    assert_eq!(refused(&too_deep, BoostRegex::new(&too_deep)), TOO_DEEP);
    // One level less compiles, and does not match `aaaa`, as in Boost.
    let shallow = format!("{}a{{2,}}{}", "(?:".repeat(19), "){2,}".repeat(19));
    let shallow = BoostRegex::new(&shallow).unwrap();
    assert_eq!(
        shallow.search(b"aaaa").unwrap().map(|found| found.range()),
        None
    );
}

/// The bytes the class-name test reports on, in the order the transcribed lists
/// use.
const CLASS_PROBE: [u8; 16] = [
    b'0', b'9', b'a', b'z', b'A', b'Z', b'_', b'-', b' ', b'\t', b'\n', 0x0B, b'\r', 0x7F, 0x00,
    0xE9,
];

/// The sixth review's findings, both in the bracket-expression parser.
/// `get_next_set_literal` reads a collating element at *either* endpoint of a
/// range, so `[A-[.a.]]` is the range `A` to `a` and only a `[` that is not
/// followed by `.` is a literal; and `cpp_regex_traits::lookup_classname` retries
/// its lookup with the name lower-cased and knows the one-letter aliases
/// `d h l s u v w` beside a `unicode` class that no `char` belongs to. Collating
/// elements and `[[:unicode:]]` are refused; the other spellings are translated.
/// Compile outcomes and answers transcribed from the oracle output for the `ADV`
/// family of the fixture.
#[test]
fn bracket_elements_and_class_names_follow_boost() {
    const COLLATING: &str = "a collating element";
    const UNICODE_CLASS: &str = "the [:unicode:] class";
    const NOT_A_POSIX_CLASS: &str = "a character class name that is not a POSIX class";
    // Boost compiles every one of these, with a range whose endpoints come from
    // its collating sequence.
    for pattern in [
        "[A-[.a.]]",
        "[a-[.a.]]",
        r"[\n-[.-.]]",
        "[!-[.].]]",
        "[^A-[.a.]]",
        "[A-[.a.]x]",
        "[q[.a.]-[.z.]]",
        r"[\x41-[.a.]]",
        "[!-[.-.][.-.]]",
        "[[.a.]-z]",
        "[[.a.]]",
        "[a[.b.]c]",
    ] {
        assert_eq!(
            refused(pattern, BoostRegex::new(pattern)),
            COLLATING,
            "{pattern:?}"
        );
    }
    // Boost rejects these (error_collate, error_range, error_ctype); the facade
    // refuses the collating element before it reaches the same fault.
    for pattern in ["[A-[.ab.]]", "[A-[.a.]-z]", "[A-[.tab.]]", "[A-[.]"] {
        assert_eq!(
            refused(pattern, BoostRegex::new(pattern)),
            COLLATING,
            "{pattern:?}"
        );
    }
    // A `[` that no `.` follows is a literal at a range end, so `[A-[x]` is
    // `A` to `[` plus `x`. (pattern, haystack, Boost's match)
    for (pattern, haystack, expected) in [
        ("[A-[x]", &b"A"[..], Some(0..1)),
        ("[A-[x]", b"[", Some(0..1)),
        ("[A-[x]", b"Z", Some(0..1)),
        ("[A-[x]", b"x", Some(0..1)),
        ("[A-[x]", b"\\", None),
        ("[A-[x]", b"a", None),
        ("[A-[]", b"[", Some(0..1)),
        ("[A-[]", b"\\", None),
        ("[A-[=a=]]x", b"=]x", Some(0..3)),
        ("[A-[=a=]]x", b"[]x", Some(0..3)),
        ("[A-[:alpha:]]x", b":]x", Some(0..3)),
        ("[A-[:alpha:]]x", b"p]x", Some(0..3)),
        (r"[A-\x5b]", b"[", Some(0..1)),
        (r"[A-\x5b]", b"\\", None),
    ] {
        let regex = BoostRegex::new(pattern).unwrap_or_else(|error| panic!("{pattern:?}: {error}"));
        assert_eq!(
            regex.search(haystack).unwrap().map(|found| found.range()),
            expected,
            "{pattern:?} on {haystack:?}"
        );
    }
    // (pattern, the bytes of CLASS_PROBE Boost matches, in that order)
    for (pattern, matched) in [
        ("[[:ALPHA:]]", &b"azAZ"[..]),
        ("[[:Alpha:]]", b"azAZ"),
        ("[[:aLPHA:]]", b"azAZ"),
        ("[[:ALNUM:]]", b"09azAZ"),
        ("[[:Word:]]", b"09azAZ_"),
        ("[[:XDIGIT:]]", b"09aA"),
        ("[[:Blank:]]", b" \t\x0b"),
        ("[[:CNTRL:]]", b"\t\n\x0b\r\x7f\x00"),
        ("[[:GRAPH:]]", b"09azAZ_-"),
        ("[[:PRINT:]]", b"09azAZ_- "),
        ("[[:PUNCT:]]", b"_-"),
        ("[[:SPACE:]]", b" \t\n\x0b\r"),
        ("[[:UPPER:]]", b"AZ"),
        ("[[:LOWER:]]", b"az"),
        // Boost's one-letter aliases, which the facade reported as unknown names.
        ("[[:d:]]", b"09"),
        ("[[:D:]]", b"09"),
        ("[[:h:]]", b" \t"),
        ("[[:H:]]", b" \t"),
        ("[[:l:]]", b"az"),
        ("[[:L:]]", b"az"),
        ("[[:s:]]", b" \t\n\x0b\r"),
        ("[[:S:]]", b" \t\n\x0b\r"),
        ("[[:u:]]", b"AZ"),
        ("[[:U:]]", b"AZ"),
        ("[[:v:]]", b"\n\x0b\r"),
        ("[[:V:]]", b"\n\x0b\r"),
        ("[[:w:]]", b"09azAZ_"),
        ("[[:W:]]", b"09azAZ_"),
        ("[[:^D:]]", b"azAZ_- \t\n\x0b\r\x7f\x00\xe9"),
        ("[[:^ALPHA:]]", b"09_- \t\n\x0b\r\x7f\x00\xe9"),
    ] {
        let regex = BoostRegex::new(pattern).unwrap_or_else(|error| panic!("{pattern:?}: {error}"));
        let found: Vec<u8> = CLASS_PROBE
            .iter()
            .copied()
            .filter(|byte| regex.is_search_match(&[*byte]).unwrap())
            .collect();
        assert_eq!(found, matched, "{pattern:?}");
    }
    // No `char` belongs to Boost's own `unicode` class, in any spelling.
    for pattern in [
        "[[:unicode:]]",
        "[[:UNICODE:]]",
        "[[:Unicode:]]",
        "[[:^unicode:]]",
        "[^[:unicode:]]",
        "[a[:unicode:]]",
        "[[:unicode:][:alpha:]]",
    ] {
        assert_eq!(
            refused(pattern, BoostRegex::new(pattern)),
            UNICODE_CLASS,
            "{pattern:?}"
        );
    }
    // `[[:<:]]` and `[[:>:]]` are Boost's word-boundary spellings, still refused.
    for pattern in ["[[:<:]]", "[[:>:]]"] {
        assert_eq!(
            refused(pattern, BoostRegex::new(pattern)),
            NOT_A_POSIX_CLASS,
            "{pattern:?}"
        );
    }
    // A name neither the exact nor the folded lookup finds stays an error, as
    // Boost's error_ctype.
    for pattern in ["[[:FOO:]]", "[[:Ascii:]]", "[[:ANY:]]", "[[:B:]]"] {
        assert!(
            matches!(BoostRegex::new(pattern), Err(Error::InvalidValue(_))),
            "{pattern:?}"
        );
    }
}

/// The sixth review's third finding: `parse_perl_extension` consumes the `P` of
/// `(?P`, reads `(?P>name)` as a recursion and hands everything else to the
/// option-group parser, so `(?Pi)` is `(?i)`, `(?P:a|b)` is `(?:a|b)` and `(?P)`
/// is an empty option group. Only the Python spellings `(?P<name>...)` and
/// `(?P=name)` are errors, because the option parser rejects their `<` and `=`.
/// Answers transcribed from the oracle output for the `ADV` family of the fixture.
#[test]
fn python_style_group_openings_follow_boost() {
    let plain = &RegexOptions::default();
    let icase = &RegexOptions {
        icase: true,
        ..RegexOptions::default()
    };
    // (options, pattern, haystack, Boost's match)
    for (options, pattern, haystack, expected) in [
        (plain, "(?Pi)A", &b"a"[..], Some(0..1)),
        (plain, "(?Pi)A", b"A", Some(0..1)),
        (plain, "(?Pi)A", b"b", None),
        (plain, "(?P:a|b)c", b"ac", Some(0..2)),
        (plain, "(?P:a|b)c", b"bc", Some(0..2)),
        (plain, "(?P:a|b)c", b"abc", Some(1..3)),
        (plain, "(?P:a|b)c", b"c", None),
        (plain, "(?P)", b"a", Some(0..0)),
        (plain, "(?P)", b"", Some(0..0)),
        (plain, "(?P)a", b"a", Some(0..1)),
        (plain, "(?P)a", b"b", None),
        (plain, "(?P:)", b"a", Some(0..0)),
        (icase, "(?P-i)A", b"a", None),
        (icase, "(?P-i)A", b"A", Some(0..1)),
        (icase, "(?Pi)a", b"A", Some(0..1)),
        (plain, "(?Pm)^a", b"b\na", Some(2..3)),
        (plain, "(?Pm)^a", b"ba", None),
        (plain, "(?Ps).", b"\n", Some(0..1)),
        (plain, "(?Pi-s).", b"\n", None),
        (plain, "(?Pi-s).", b"A", Some(0..1)),
        (plain, "a(?Pi)b", b"ab", Some(0..2)),
        (plain, "a(?Pi)b", b"aB", Some(0..2)),
        (plain, "a(?Pi)b", b"Ab", None),
        (plain, "(?Pi:a)", b"A", Some(0..1)),
        (plain, "(?Pi:a)", b"b", None),
        (icase, "(?P-i:A)", b"A", Some(0..1)),
        (icase, "(?P-i:A)", b"a", None),
        (plain, "(?Pi)(?P-i)a", b"a", Some(0..1)),
        (plain, "(?Pi)(?P-i)a", b"A", None),
    ] {
        let regex = BoostRegex::with_options(pattern, *options)
            .unwrap_or_else(|error| panic!("{pattern:?}: {error}"));
        assert_eq!(
            regex.search(haystack).unwrap().map(|found| found.range()),
            expected,
            "{pattern:?} icase={} on {haystack:?}",
            options.icase
        );
    }
    // `(?P>name)` is a recursion by name, and `x` in the option list is the
    // extended modifier: both refused, as their `(?...` spellings are.
    assert_eq!(
        refused("(?<n>a)(?P>n)", BoostRegex::new("(?<n>a)(?P>n)")),
        "a recursive sub-expression"
    );
    for pattern in ["(?Pim-sx:a|b)", "(?Px)a b"] {
        assert_eq!(
            refused(pattern, BoostRegex::new(pattern)),
            "the x (extended) modifier",
            "{pattern:?}"
        );
    }
    // Everything Boost rejects here the facade rejects too, as an invalid pattern
    // rather than an unsupported construct.
    for pattern in [
        "(?P<n>a)",
        "(?P=n)",
        "(?<n>a)(?P=n)",
        "(?Pz)",
        "(?P",
        "(?Pi",
        "(?PP)",
        "(?PPi)",
        "(?P&n)",
        "(?P1)",
        "(?P?i)",
    ] {
        assert!(
            matches!(BoostRegex::new(pattern), Err(Error::InvalidValue(_))),
            "{pattern:?}"
        );
    }
}

/// An expression that needs no backtracking was handed whole to an automaton,
/// whose slowest mode does work proportional to the haystack length times the
/// expression's atoms, which repeat bounds multiply: `\w{0,9999}b` took 8 s over
/// 100 KB, where Boost reports `error_complexity` in 80 ms. Above
/// `MAX_AUTOMATON_ATOMS`, divided by the budget's scale factor for a long haystack,
/// the search runs on the backtracking machine and spends the budget.
#[test]
fn large_automata_spend_the_budget() {
    // Boost's answer on a short haystack, from the oracle output.
    let wide = BoostRegex::new(r"\w{0,9999}b").unwrap();
    assert_eq!(
        wide.search(b"aab").unwrap().map(|found| found.range()),
        Some(0..3)
    );
    let mut haystack = vec![b'a'; 100_000];
    haystack.push(b'c');
    let haystack: &'static [u8] = Box::leak(haystack.into_boxed_slice());
    match within(120, r"\w{0,9999}b", move || wide.search(haystack)) {
        Err(Error::InvalidValue(message)) => assert!(
            message.contains("exceeded the backtracking limit of 4000000 steps"),
            "{message}"
        ),
        other => panic!("expected the backtracking limit, got {other:?}"),
    }

    // The limit behaviour at the first two tiers, with a budget of one step, which
    // only a search on the backtracking machine can exceed. `\w{0,4095}b` has exactly
    // MAX_AUTOMATON_ATOMS atoms; its match in `ab` followed by twenty `a` takes more
    // than four steps to backtrack to, and the spaces that follow end every run.
    let one_step = RegexOptions {
        backtrack_limit: 1,
        ..RegexOptions::default()
    };
    let at_limit = format!(r"\w{{0,{}}}b", MAX_AUTOMATON_ATOMS - 1);
    let above_limit = format!(r"\w{{0,{MAX_AUTOMATON_ATOMS}}}b");
    let text = |length: usize| {
        let mut text = b"ab".to_vec();
        text.extend(std::iter::repeat_n(b'a', 20));
        text.resize(length, b' ');
        text
    };
    let at_limit = BoostRegex::with_options(&at_limit, one_step).unwrap();
    assert_eq!(
        at_limit
            .search(&text(BACKTRACK_LIMIT_BYTES))
            .unwrap()
            .map(|found| found.range()),
        Some(0..2)
    );
    for (regex, length, limit) in [
        (&at_limit, BACKTRACK_LIMIT_BYTES + 1, 4),
        (
            &BoostRegex::with_options(&above_limit, one_step).unwrap(),
            64,
            1,
        ),
    ] {
        match regex.search(&text(length)) {
            Err(Error::InvalidValue(message)) => assert!(
                message.contains(&format!("exceeded the backtracking limit of {limit} steps")),
                "{length} bytes: {message}"
            ),
            other => panic!("{length} bytes: expected the backtracking limit, got {other:?}"),
        }
    }
}

/// The counted spelling ends in `(?=)(?=)` so that `fancy-regex` runs the expression
/// on its backtracking machine: one trailing `(?=)` is removed by the engine's
/// `optimize_trailing_lookahead`, after which a trailing run that needs no
/// backtracking, or a whole expression that needs none, is handed to an automaton
/// that spends no budget. Were a crate update to undo this, these searches would
/// answer (slowly) instead of stopping with the budget error.
#[test]
fn forced_backtracking_spends_the_budget() {
    let upper: &'static [u8] = Box::leak(vec![b'A'; 64 * 1024].into_boxed_slice());
    for pattern in ["(?=A).*B", ".{0,5000}B", "(?:A|B){0,4096}C"] {
        let outcome = within(120, pattern, move || {
            BoostRegex::new(pattern).unwrap().search(upper)
        });
        match outcome {
            Err(Error::InvalidValue(message)) => assert!(
                message.contains("exceeded the backtracking limit of 1000000 steps"),
                "{pattern:?}: {message}"
            ),
            other => panic!("{pattern:?}: expected the backtracking limit, got {other:?}"),
        }
    }
}

/// The engine's positive lookaheads are not atomic: when what follows fails, it
/// backtracks into the lookahead's body, which Boost never does. A body with
/// ambiguous repeats then exhausts the budget on a short haystack where Boost
/// answers at once (no match for each of these); the search still ends within the
/// bound, with Boost's answer or the budget error.
#[test]
fn lookahead_bodies_backtrack_within_the_budget() {
    for (pattern, haystack) in [
        ("(?=(?:a|aa)*)b", format!("{}c", "a".repeat(60))),
        ("(?=a*a*a*)b", format!("{}c", "a".repeat(300))),
        ("(?=.*a).*b", "a".repeat(4000)),
        ("(?=[A-Z]*K)[A-Z]+R", "A".repeat(12_000)),
    ] {
        let outcome = within(120, pattern, move || {
            BoostRegex::new(pattern)
                .unwrap()
                .search(haystack.as_bytes())
        });
        match outcome {
            Ok(None) => {}
            Err(Error::InvalidValue(message)) if message.contains("backtracking limit") => {}
            other => panic!("{pattern:?}: {other:?}"),
        }
    }
}

/// `fancy-regex` counts backtracking steps, not the work between them: an
/// automaton it runs anchored can scan to the end of the haystack, a counted repeat
/// pushes no backtracking branch below its minimum, and the branches of an atomic
/// group or a negative lookahead are discarded without being counted. Repeated from
/// every start position, each of these ran for seconds to minutes on 64 KiB
/// (`(?:(?=A)A){999999}` 26.5 s, `(?!A*$)x` 35 s on 100 KB). Now the budget counts
/// that work, so each search here stops with the budget error, or the pattern is
/// refused.
#[test]
fn uncounted_engine_work_is_bounded() {
    let limit = "exceeded the backtracking limit of 1000000 steps";
    let upper: &'static [u8] = Box::leak(vec![b'A'; 16 * 1024].into_boxed_slice());
    let lower: &'static [u8] = Box::leak(vec![b'a'; 16 * 1024].into_boxed_slice());
    let lines: &'static [u8] = Box::leak(b"a\n".repeat(8 * 1024).into_boxed_slice());
    for (pattern, haystack) in [
        ("(?:(?=A)A){999999}", upper),
        (r"(A)\1{999999}", upper),
        ("(?=A).*B", upper),
        ("(?=.*B)A", upper),
        ("^(.+): (.+)", lines),
        ("(?=A*$)x", upper),
        ("(?<=a{255})b", lower),
    ] {
        let outcome = within(120, pattern, move || {
            BoostRegex::new(pattern).unwrap().search(haystack)
        });
        match outcome {
            Err(Error::InvalidValue(message)) => {
                assert!(message.contains(limit), "{pattern:?}: {message}")
            }
            other => panic!("{pattern:?}: expected the backtracking limit, got {other:?}"),
        }
    }
    // A negative lookbehind discards its branches too, including those of a
    // lookahead nested in it that scans forward (found by fuzzing: 12 s on 64 KiB).
    for pattern in [
        "(?>A*)B",
        "(?>A*$)B",
        "(?!A*$)x",
        "(?<!(?=.*b)a)b",
        r"(?<!(?=.*b)\w)[ab](?<=a)b",
    ] {
        let category = within(60, pattern, move || {
            refused(pattern, BoostRegex::new(pattern))
        });
        assert_eq!(category, DISCARDED, "{pattern:?}");
    }
}

/// The limit behaviour the facade chose. One search shares its backtracking
/// budget across all its start positions, so the budget grows with the haystack:
/// by 1 up to `BACKTRACK_LIMIT_BYTES`, then by 4, 16 and at most
/// `MAX_BACKTRACK_SCALE`. The engine's fixed stack of 1,000,000 backtracking
/// branches is a second limit. Both are errors, never a different answer.
#[test]
fn backtracking_budget_scales_with_the_haystack() {
    // A line-anchored scan takes a few steps per line: 500 KB exceed one unscaled
    // budget but not the scaled one, and Boost finds no match.
    let caret = BoostRegex::new("^x").unwrap();
    let lines = b"a\n".repeat(250_000);
    assert_eq!(
        within(120, "^x", move || caret.search(&lines)).unwrap(),
        None
    );

    let options = RegexOptions {
        backtrack_limit: 10_000,
        ..RegexOptions::default()
    };
    // Twenty two-way choices before a failing tail: the stack stays a few dozen
    // branches deep however long the haystack, so the budget, not the stack, ends
    // every search, and the message names the scaled budget.
    let catastrophic = BoostRegex::with_options(r"(?:a|a){20}(?!x)b", options).unwrap();
    for (length, scale) in [
        (BACKTRACK_LIMIT_BYTES, 1),
        (BACKTRACK_LIMIT_BYTES + 1, 4),
        (4 * BACKTRACK_LIMIT_BYTES, 4),
        (4 * BACKTRACK_LIMIT_BYTES + 1, 16),
        (16 * BACKTRACK_LIMIT_BYTES, 16),
        (16 * BACKTRACK_LIMIT_BYTES + 1, MAX_BACKTRACK_SCALE),
    ] {
        let mut haystack = vec![b'a'; length - 1];
        haystack.push(b'c');
        match catastrophic.search(&haystack) {
            Err(Error::InvalidValue(message)) => assert!(
                message.contains(&format!(
                    "exceeded the backtracking limit of {} steps",
                    10_000 * scale
                )),
                "{length} bytes: {message}"
            ),
            other => panic!("{length} bytes: expected the backtracking limit, got {other:?}"),
        }
    }

    // A greedy repeat before a line anchor pushes one branch per byte it takes.
    let scan = BoostRegex::new(r"=(?<SCAN>\d+)$").unwrap();
    let digits = |count: usize| {
        let mut haystack = vec![b'='];
        haystack.extend(std::iter::repeat_n(b'1', count));
        haystack.push(b'x');
        haystack
    };
    assert_eq!(scan.search(&digits(100_000)).unwrap(), None);
    match scan.search(&digits(1_000_000)) {
        Err(Error::InvalidValue(message)) => assert!(
            message.contains("exceeded the engine's backtracking stack of 1000000 entries"),
            "{message}"
        ),
        other => panic!("expected the stack limit, got {other:?}"),
    }
}

#[test]
fn backtracking_limit_is_an_error_not_a_hang() {
    let options = RegexOptions {
        backtrack_limit: 10_000,
        ..RegexOptions::default()
    };
    let regex = BoostRegex::with_options(r"(?:a|a)+(?!x)b", options).unwrap();
    let haystack = format!("{}c", "a".repeat(40));
    match regex.search(haystack.as_bytes()) {
        Err(Error::InvalidValue(message)) => {
            assert!(message.contains("backtracking limit"), "{message}")
        }
        other => panic!("expected the backtracking limit, got {other:?}"),
    }
    let mut tokens = regex.tokens(haystack.as_bytes(), &[-1]).unwrap();
    assert!(matches!(tokens.next(), Some(Err(Error::InvalidValue(_)))));
    assert!(tokens.next().is_none());
    assert!(regex.search(b"aab").unwrap().is_some());
}

/// Arbitrary patterns, including malformed and non-ASCII ones, return a value
/// or an error, never a panic, and every search of them answers or reports a limit.
#[test]
fn arbitrary_patterns_never_panic() {
    const PIECES: &[&str] = &[
        "a", "b", "0", "_", "-", ",", ".", "^", "$", "|", "(", ")", "[", "]", "{", "}", "*", "+",
        "?", "\\", ":", "=", "!", "<", ">", "'", "#", "i", "s", "m", "x", "d", "w", "1", "2",
        "\u{e9}", " ", "\\d", "(?", "(?<", "(?<=", "[:", ":]", "{2}", "{1,3}", "\\x", "\\b", "P",
        "R", "&", "Q", "E", "\\n",
    ];
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let inputs: [&[u8]; 5] = [b"", b"ab\r\n0_", b"aaab", b"\xff\x00a", b"x-y,z"];
    let mut compiled = 0;
    for _ in 0..3000 {
        let length = (next() % 16) as usize;
        let pattern: String = (0..length)
            .map(|_| PIECES[(next() % PIECES.len() as u64) as usize])
            .collect();
        let Ok(regex) = BoostRegex::new(&pattern) else {
            continue;
        };
        compiled += 1;
        for input in inputs {
            let _ = regex.search(input);
            let _ = regex.full_match(input);
            if let Ok(tokens) = regex.tokens(input, &[-1, 0]) {
                assert!(tokens.take(4 * input.len() + 8).count() <= 4 * input.len() + 8);
            }
        }
    }
    assert!(compiled > 100, "only {compiled} random patterns compiled");
}

#[test]
fn regex_is_send_and_sync() {
    fn check<T: Send + Sync>() {}
    check::<BoostRegex>();
}
