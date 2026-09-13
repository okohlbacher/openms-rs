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
//! adversarial inputs, plus syntax probes and grammar-generated expressions. The
//! facade must reproduce every compile outcome, every match, every group span,
//! every named-group lookup and every token, byte for byte.
//!
//! The class-test section is tier 3 (literals transcribed from the pinned class
//! tests); the limit, error and robustness sections are tier 4.

use openms::Error;
use openms::concept::boost_regex::{
    BoostRegex, MAX_GROUP_DEPTH, MAX_PATTERN_BYTES, RegexOptions, SubMatch,
};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::PathBuf;

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
/// with `Error::Unsupported`, as `(flags, pattern)`. `docs/BOOST_REGEX_SUPPORT.md`
/// documents each construct.
const EXPECTED_UNSUPPORTED: &[(&str, &str)] = &[
    ("-", "((?=a))*"),
    ("-", "()*"),
    ("-", "()+"),
    ("-", "(*ACCEPT)"),
    ("-", "(*FAIL)"),
    ("-", "(?!(a))b"),
    ("-", "(?(1)a|b)"),
    ("-", "(?:$)+"),
    ("-", "(?:(?=(a))|b)*"),
    ("-", "(?:\\b)*"),
    ("-", "(?:^)*"),
    ("-", "(?<=(a))b"),
    ("-", "(?<=(a)|b)c"),
    ("-", "(?<=\\Z)"),
    ("-", "(?<=a)\\Z"),
    ("-", "(?<n-x>a)"),
    ("-", "(?<n>a)\\k<n>"),
    ("-", "(?=(\\d+))\\d"),
    ("-", "(?=(a))*"),
    ("-", "(?=(a))*?b"),
    ("-", "(?=(a)){0}"),
    ("-", "(?=(a)){2}"),
    ("-", "(?i)*"),
    ("-", "(?i)+a"),
    ("-", "(?x)a b"),
    ("-", "(?|(a)|(b))"),
    ("-", "(a)(?(1)b|c)"),
    ("-", "(a)(?1)"),
    ("-", "(a)(b)(c)(d)(e)(f)(g)(h)(i)(j)\\10"),
    ("-", "(a)\\10"),
    ("-", "(a)\\g1"),
    ("-", "(a)\\g{-1}"),
    ("-", "(a)\\g{1}"),
    ("-", "(a*)*"),
    ("-", "(a*)+"),
    ("-", "(|a)*"),
    ("-", "[[.a.]]"),
    ("-", "[[=a=]]"),
    ("-", "[\\0]"),
    ("-", "[\\y]"),
    ("-", "[a-\\d]"),
    ("-", "[é]"),
    ("-", "\\0"),
    ("-", "\\012"),
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
    ("-", "\\cA"),
    ("-", "\\cZ"),
    ("-", "\\ca"),
    ("-", "\\j"),
    ("-", "\\l"),
    ("-", "\\pL"),
    ("-", "\\p{L}"),
    ("-", "\\q"),
    ("-", "\\u"),
    ("-", "\\x80"),
    ("-", "\\xc3"),
    ("-", "\\xff"),
    ("-", "\\y"),
    ("-", "^\\Z"),
    ("-", "a(?i)*"),
    ("-", "a*+"),
    ("-", "a++"),
    ("-", "a?+"),
    ("-", "a\\Kb"),
    ("-", "a\\Z"),
    ("-", "a\\Z\\n"),
    ("-", "a{2}+"),
    ("-", "a٣"),
    ("-", "é"),
    ("-", "é+"),
];

/// Refusals over the whole corpus, fuzz family included, by the construct the
/// facade names in its `Error::Unsupported` message.
const EXPECTED_REFUSALS: &[(&str, usize)] = &[
    (r"\G", 1),
    (r"\K", 2),
    (r"\Q...\E quoting", 4),
    (
        r"\Z (Boost never tries it at a form feed when it starts the expression)",
        7,
    ),
    ("a backreference with more than one digit", 2),
    ("a backtracking control verb", 2),
    ("a branch-reset group", 1),
    (
        "a capturing group inside a lookaround or atomic group (Boost keeps its capture when the surrounding match backtracks)",
        111,
    ),
    ("a character property escape", 3),
    ("a collating element", 1),
    ("a conditional expression", 2),
    ("a control-character escape", 3),
    ("a group name outside [A-Za-z0-9_]", 1),
    ("a hexadecimal escape above 0x7F", 3),
    ("a named or relative backreference", 4),
    ("a named-character, line-ending or grapheme escape", 4),
    ("a non-ASCII byte", 4),
    ("a possessive quantifier", 4),
    ("a recursive sub-expression", 1),
    (
        "a repeat of a group that can match the empty string and captures (Boost records a final empty iteration)",
        34,
    ),
    ("a repeat of a group that only asserts a position", 14),
    (
        "a repeat of a modifier group that switches case sensitivity (Boost undoes the switch when the empty repetition is abandoned)",
        3,
    ),
    ("an equivalence class", 1),
    (
        "an escape inside a character class without a translated meaning",
        3,
    ),
    ("an escape letter without a translated meaning", 7),
    ("an octal escape", 2),
    ("the x (extended) modifier", 1),
];

/// Cases compared against Boost: the fixture's size, so a truncated fixture fails.
const EXPECTED_CASES: usize = 146_437;

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
        "refused non-fuzz patterns differ from the documented list:\n{listing}"
    );
    let expected_refusals: BTreeMap<String, usize> = EXPECTED_REFUSALS
        .iter()
        .map(|(category, count)| ((*category).to_string(), *count))
        .collect();
    assert_eq!(
        refusals, expected_refusals,
        "refusal categories differ from the documented counts:\n{categories}"
    );
    assert_eq!(cases, EXPECTED_CASES, "the fixture lost or gained cases");
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
/// or an error, never a panic, and matching them stays bounded.
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
