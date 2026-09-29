// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `DATASTRUCTURES/StringUtils.h` against the Linux x86_64 Release build.
//!
//! `tests/data/string_utils_oracle.tsv` is `results/records.tsv` of
//! `../oracle/B8a-stringutils` (driver `driver.cpp`, linked against
//! `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576` on
//! ibminode06, two runs byte-identical): 5,068 calls over the whole header, one
//! record per call, `function <TAB> arguments <TAB> result`. Every string is
//! escaped (`\\`, `\|`, `\xNN` for every byte outside printable ASCII), lists
//! are joined with `|`, doubles are `d:` and floats `f:` plus their bits, and a
//! thrown exception is `ERR:<name>:<what()>`. [`oracle_records_match`]
//! recomputes every record with the port and compares the text.
//!
//! The class-test sections below port the literals of
//! `src/tests/class_tests/openms/source/StringUtils_test.cpp` (`class_test_*`)
//! and `String_test.cpp` (`string_test_*`) at core `bc9cc12` (tier 3); the
//! oracle is tier 1.

use openms::Error;
use openms::data_structures::string_utils::{
    self as su, QuotingMethod, Split, ToStr, extract_double, extract_int, number, number_length,
    number_truncated, skip_non_whitespace, skip_whitespace, to_double, to_float, to_int32,
    to_int64, to_str, to_str_f32, to_str_f32_low_p, to_str_low_p,
};
use openms::metadata::MetaValue;
use openms::param::ParamValue;

const ORACLE: &str = include_str!("data/string_utils_oracle.tsv");

/// The driver's `esc()`.
fn esc(bytes: &[u8]) -> String {
    let mut out = String::new();
    for &b in bytes {
        match b {
            b'\\' => out.push_str("\\\\"),
            b'|' => out.push_str("\\|"),
            0x20..=0x7e => out.push(char::from(b)),
            _ => out.push_str(&format!("\\x{b:02x}")),
        }
    }
    out
}

/// Inverse of [`esc`].
fn unesc(text: &str) -> Vec<u8> {
    let bytes = text.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'\\' {
            match bytes.get(i + 1) {
                Some(b'x') => {
                    let hex = std::str::from_utf8(&bytes[i + 2..i + 4]).unwrap();
                    out.push(u8::from_str_radix(hex, 16).unwrap());
                    i += 4;
                }
                Some(&c) => {
                    out.push(c);
                    i += 2;
                }
                None => panic!("dangling escape in {text:?}"),
            }
        } else {
            out.push(bytes[i]);
            i += 1;
        }
    }
    out
}

/// Split an escaped field on its unescaped `|` separators.
fn fields(text: &str) -> Vec<String> {
    let mut out = vec![String::new()];
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        match c {
            '\\' => {
                let last = out.last_mut().unwrap();
                last.push('\\');
                last.push(chars.next().unwrap());
            }
            '|' => out.push(String::new()),
            _ => out.last_mut().unwrap().push(c),
        }
    }
    out
}

fn text(field: &str) -> Option<String> {
    String::from_utf8(unesc(field)).ok()
}

fn f64_bits(field: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(field.strip_prefix("d:").unwrap(), 16).unwrap())
}

fn f32_bits(field: &str) -> f32 {
    f32::from_bits(u32::from_str_radix(field.strip_prefix("f:").unwrap(), 16).unwrap())
}

fn dbits(value: f64) -> String {
    format!("d:{:016x}", value.to_bits())
}

fn fbits(value: f32) -> String {
    format!("f:{:08x}", value.to_bits())
}

/// A port error in the driver's `ERR:<name>:<what()>` form, the exception name
/// taken from the error the port documents for it.
fn err(error: &Error) -> String {
    let (name, message) = match error {
        Error::InvalidRange(m) if m.contains("too large") => ("IndexOverflow", m.as_str()),
        Error::InvalidRange(m) if m.contains("too small") => ("IndexUnderflow", m.as_str()),
        Error::InvalidValue(m) => ("ConversionError", m.as_str()),
        other => panic!("unexpected error kind {other:?}"),
    };
    format!("ERR:{name}:{}", esc(message.as_bytes()))
}

/// Whether the port refused because the source's result is not UTF-8.
fn utf8_refusal(error: &Error) -> bool {
    matches!(error, Error::InvalidValue(m) if m.contains("UTF-8"))
}

fn text_or_err(result: openms::Result<String>) -> std::result::Result<String, Error> {
    result.map(|t| esc(t.as_bytes()))
}

fn split_text(split: &Split, with_count: bool) -> String {
    let parts: Vec<String> = split.parts.iter().map(|p| esc(p.as_bytes())).collect();
    if with_count {
        format!(
            "{}#{}#{}",
            u8::from(split.split),
            split.parts.len(),
            parts.join("|")
        )
    } else {
        format!("{}#{}", u8::from(split.split), parts.join("|"))
    }
}

fn method(name: &str) -> QuotingMethod {
    match name {
        "NONE" => QuotingMethod::None,
        "ESCAPE" => QuotingMethod::Escape,
        "DOUBLE" => QuotingMethod::Double,
        other => panic!("method {other}"),
    }
}

fn meta_value(name: &str) -> Option<MetaValue> {
    Some(match name {
        "empty" => MetaValue::default(),
        "string" => MetaValue::from("a, b"),
        "int" => MetaValue::from(-7_i64),
        "int64" => MetaValue::from(9_007_199_254_740_993_i64),
        "double1" => MetaValue::try_from(3.25).unwrap(),
        "double2" => MetaValue::try_from(1.23456789e-7).unwrap(),
        "double3" => MetaValue::try_from(12345.6789).unwrap(),
        "double4" => MetaValue::try_from(0.0625).unwrap(),
        "stringlist" => MetaValue::from(vec!["a".to_owned(), "b, c".to_owned(), String::new()]),
        "intlist" => MetaValue::from(vec![1_i64, -2, 3]),
        "emptystringlist" => MetaValue::from(Vec::<String>::new()),
        "emptyintlist" => MetaValue::from(Vec::<i64>::new()),
        "emptydoublelist" => MetaValue::try_from(Vec::<f64>::new()).unwrap(),
        // NaN and the infinities are refused by the public constructors
        // (MetaValue::validate); their text is the `toStr` rows'.
        "double5" | "double6" | "doublelist" => return None,
        other => panic!("DataValue case {other}"),
    })
}

fn param_value(name: &str) -> ParamValue {
    match name {
        "empty" => ParamValue::Empty,
        "string" => ParamValue::String("a, b".into()),
        "int" => ParamValue::Integer(-7),
        "double1" => ParamValue::Float(3.25),
        "double2" => ParamValue::Float(1.23456789e-7),
        "double3" => ParamValue::Float(12345.6789),
        "double4" => ParamValue::Float(0.0625),
        "stringlist" => ParamValue::StringList(vec!["a".into(), "b, c".into(), String::new()]),
        "intlist" => ParamValue::IntegerList(vec![1, -2, 3]),
        "doublelist" => ParamValue::FloatList(vec![1.5, 1e-5, 12345.678, 0.0625, -0.0]),
        other => panic!("ParamValue case {other}"),
    }
}

/// The port's text for one oracle record, or `None` when the record's input is
/// not UTF-8 (a Rust `str` cannot hold it).
fn port_result(function: &str, args: &str) -> Option<std::result::Result<String, Error>> {
    let f = fields(args);
    let ok = |s: String| Some(Ok(s));
    match function {
        "toInt32" => Some(to_int32(&text(args)?).map(|v| v.to_string())),
        "toInt64" => Some(to_int64(&text(args)?).map(|v| v.to_string())),
        "toFloat" => Some(to_float(&text(args)?).map(fbits)),
        "toDouble" => Some(to_double(&text(args)?).map(dbits)),
        "extractDouble" => ok(match extract_double(&text(args)?) {
            Some((value, used)) => format!("1|{used}|{}", dbits(value)),
            None => "0|0|-".into(),
        }),
        "extractInt" => ok(match extract_int(&text(args)?) {
            Some((value, used)) => format!("1|{used}|{value}"),
            None => "0|0|-".into(),
        }),
        "toStr(double,true)" => ok(esc(to_str(f64_bits(args)).as_bytes())),
        "toStr(double,false)" => ok(esc(to_str_low_p(f64_bits(args)).as_bytes())),
        "toStr(float,true)" => ok(esc(to_str_f32(f32_bits(args)).as_bytes())),
        "toStr(float,false)" => ok(esc(to_str_f32_low_p(f32_bits(args)).as_bytes())),
        "operator+(double)" => {
            let mut s = String::from("x=");
            f64_bits(args).append_to_str(&mut s);
            ok(esc(s.as_bytes()))
        }
        "operator+(float)" => {
            let mut s = String::from("x=");
            f32_bits(args).append_to_str(&mut s);
            ok(esc(s.as_bytes()))
        }
        "appendToStrLowP(double)" => {
            let mut s = String::from("p");
            su::append_to_str_low_p(f64_bits(args), &mut s);
            ok(esc(s.as_bytes()))
        }
        "appendToStrLowP(float)" => {
            let mut s = String::from("p");
            su::append_to_str_low_p_f32(f32_bits(args), &mut s);
            ok(esc(s.as_bytes()))
        }
        "toStr(long long)" => ok(args.parse::<i64>().unwrap().to_str_text()),
        "operator+(long long)" => {
            let mut s = String::from("n");
            args.parse::<i64>().unwrap().append_to_str(&mut s);
            ok(s)
        }
        "toStr(unsigned long long)" => ok(args.parse::<u64>().unwrap().to_str_text()),
        "toStr(short)" => ok(args.parse::<i16>().unwrap().to_str_text()),
        "toStr(unsigned short)" => ok(args.parse::<u16>().unwrap().to_str_text()),
        "toStr(unsigned int)" => ok(args.parse::<u32>().unwrap().to_str_text()),
        "toStr(char)" | "toStr(const char*)" => {
            let input = text(args)?;
            let rendered = match input.chars().next() {
                Some(c) if function == "toStr(char)" => c.to_str_text(),
                _ => input.as_str().to_str_text(),
            };
            ok(esc(rendered.as_bytes()))
        }
        "numberLength" => {
            let value = f64_bits(&f[0]);
            let length: u32 = f[1].parse().unwrap();
            Some(text_or_err(number_length(value, length)))
        }
        "number" => {
            let value = f64_bits(&f[0]);
            let digits: u32 = f[1].parse().unwrap();
            let truncated = number_truncated(value, digits);
            // The strict form agrees wherever it succeeds.
            if let Ok(strict) = number(value, digits) {
                assert_eq!(strict, truncated, "number({value}, {digits})");
            } else {
                assert!(truncated.len() == 63, "number({value}, {digits}) refused");
            }
            ok(esc(truncated.as_bytes()))
        }
        "random" => {
            let length: u32 = args.parse().unwrap();
            let r = su::random(length).unwrap();
            assert_eq!(r.len(), length as usize);
            assert!(r.bytes().all(|b| b.is_ascii_alphanumeric()));
            ok(if length == 64 {
                "len=64 alnum".into()
            } else {
                esc(r.as_bytes())
            })
        }
        "concatenate" | "concatenate(target)" => {
            let items = ["a", "", "b c", "|"];
            ok(esc(match args {
                "a||b c|\\|;glue=, " => su::concatenate(items, ", ").unwrap(),
                "a||b c|\\|;glue=" => su::concatenate(items, "").unwrap(),
                ";glue=x" => su::concatenate(Vec::<String>::new(), "x").unwrap(),
                other => panic!("concatenate {other}"),
            }
            .as_bytes()))
        }
        "toStr(DataValue,true)" | "toStr(DataValue,false)" | "appendToStr(DataValue,false)" => {
            let value = meta_value(args)?;
            let full = function == "toStr(DataValue,true)";
            let mut rendered = value.to_text(full).unwrap();
            if function.starts_with("appendToStr") {
                rendered.insert(0, 'p');
            }
            ok(esc(rendered.as_bytes()))
        }
        "toStr(ParamValue,true)" | "toStr(ParamValue,false)" => {
            let full = function == "toStr(ParamValue,true)";
            ok(esc(param_value(args).to_text(full).unwrap().as_bytes()))
        }
        _ => string_function(function, &f),
    }
}

/// The records whose first argument is a text.
fn string_function(function: &str, f: &[String]) -> Option<std::result::Result<String, Error>> {
    let s = text(&f[0])?;
    let ok = |s: String| Some(Ok(esc(s.as_bytes())));
    let bit = |b: bool| Some(Ok(String::from(if b { "1" } else { "0" })));
    let mut owned = s.clone();
    match function {
        "trim" => {
            su::trim(&mut owned);
            assert_eq!(su::trimmed(&s), owned);
            ok(owned)
        }
        "simplify" => {
            su::simplify(&mut owned);
            ok(owned)
        }
        "removeWhitespaces" => {
            su::remove_whitespaces(&mut owned);
            ok(owned)
        }
        "toUpper" => {
            su::to_upper(&mut owned);
            assert_eq!(su::to_uppered(&s), owned);
            ok(owned)
        }
        "toLower" => {
            su::to_lower(&mut owned);
            assert_eq!(su::to_lowered(&s), owned);
            ok(owned)
        }
        "firstToUpper" => {
            su::first_to_upper(&mut owned);
            ok(owned)
        }
        "reverse" => Some(su::reverse(&mut owned).map(|()| esc(owned.as_bytes()))),
        "skipWhitespace" => Some(Ok(skip_whitespace(&s).to_string())),
        "skipNonWhitespace" => Some(Ok(skip_non_whitespace(&s).to_string())),
        "isQuoted" => bit(su::is_quoted(&s, '"')),
        "hasPrefix" => bit(s.starts_with("a")),
        "hasSuffix" => bit(s.ends_with("c")),
        "hasSubstring" => bit(s.contains(",")),
        "has(Byte)" => bit(s.as_bytes().contains(&b' ')),
        "remove" => {
            su::remove(&mut owned, ' ');
            ok(owned)
        }
        "ensureLastChar" => {
            su::ensure_last_char(&mut owned, '/');
            ok(owned)
        }
        "substitute(char)" => {
            su::substitute_char(&mut owned, 'a', 'Z');
            assert_eq!(su::substituted_char(&s, 'a', 'Z'), owned);
            ok(owned)
        }
        "substitute(string)" => {
            let from = text(&f[1])?;
            let to = text(&f[2])?;
            su::substitute(&mut owned, &from, &to).unwrap();
            assert_eq!(su::substituted(&s, &from, &to).unwrap(), owned);
            ok(owned)
        }
        "chop" | "prefix(size_t)" | "suffix(size_t)" | "prefix(Int)" | "suffix(Int)"
        | "fillLeft" | "fillRight" | "substr(pos)" | "substr(pos,2)" => {
            let n: i64 = f[1].parse().unwrap();
            let size = usize::try_from(n).unwrap_or(0);
            let result = match function {
                "chop" => su::chop(&s, size),
                "prefix(size_t)" => su::prefix(&s, size),
                "suffix(size_t)" => su::suffix(&s, size),
                "prefix(Int)" => su::prefix_signed(&s, i32::try_from(n).unwrap()),
                "suffix(Int)" => su::suffix_signed(&s, i32::try_from(n).unwrap()),
                "fillLeft" => su::fill_left(&mut owned, 'x', u32::try_from(n).unwrap())
                    .map(|()| owned.clone()),
                "fillRight" => su::fill_right(&mut owned, 'x', u32::try_from(n).unwrap())
                    .map(|()| owned.clone()),
                "substr(pos)" => su::substr(&s, size, None),
                _ => su::substr(&s, size, Some(2)),
            };
            Some(text_or_err(result))
        }
        "prefix(char)" => ok(su::prefix_until(&s, ',')),
        "suffix(char)" => ok(su::suffix_after(&s, ',')),
        "quote" | "unquote" => {
            let q = text(&f[1])?.chars().next().unwrap();
            let m = method(&f[2]);
            let result = if function == "quote" {
                su::quote(&mut owned, q, m)
            } else {
                su::unquote(&mut owned, q, m)
            };
            Some(result.map(|()| esc(owned.as_bytes())))
        }
        "split_quoted" => {
            let q = text(&f[1])?.chars().next().unwrap();
            let m = method(&f[2]);
            let splitter = text(&f[3])?;
            let with_count = splitter.is_empty();
            Some(su::split_quoted(&s, &splitter, q, m).map(|r| split_text(&r, with_count)))
        }
        "split(char)" => {
            let splitter = text(&f[1])?.chars().next().unwrap();
            Some(su::split(&s, splitter, f[2] == "1").map(|r| split_text(&r, true)))
        }
        "split(string)" => {
            let splitter = text(&f[1])?;
            if splitter.is_empty() && !s.is_ascii() {
                // The source cuts into single bytes, which are not UTF-8; the
                // port cuts into characters (native difference).
                return Some(Err(Error::InvalidValue(
                    "one part per byte does not give UTF-8".into(),
                )));
            }
            Some(su::split_str(&s, &splitter).map(|r| split_text(&r, true)))
        }
        other => panic!("no port mapping for oracle function {other}"),
    }
}

#[test]
fn oracle_records_match() {
    let mut compared = 0usize;
    let mut not_utf8_inputs = 0usize;
    let mut not_utf8_results = 0usize;
    let mut undefined_in_source = 0usize;
    let mut skipped_nonfinite_values = 0usize;
    let mut mismatches = Vec::new();
    for line in ORACLE.lines() {
        let mut columns = line.splitn(3, '\t');
        let function = columns.next().unwrap();
        let args = columns.next().unwrap();
        let expected = columns.next().unwrap();
        let Some(actual) = port_result(function, args) else {
            if function.contains("DataValue") {
                skipped_nonfinite_values += 1;
            } else {
                not_utf8_inputs += 1;
                assert!(args.contains("\\xff"), "{function} {args} skipped");
            }
            continue;
        };
        match actual {
            Ok(actual) => {
                compared += 1;
                if actual != expected {
                    mismatches.push(format!(
                        "{function}\t{args}\n  C++:  {expected}\n  Rust: {actual}"
                    ));
                }
            }
            Err(error)
                if utf8_refusal(&error)
                    || matches!(&error, Error::InvalidValue(m) if m.contains("not give UTF-8")) =>
            {
                // The source's result is bytes that are not UTF-8.
                let bytes = unesc(expected);
                assert!(
                    std::str::from_utf8(&bytes).is_err(),
                    "{function} {args}: refused, but C++ gave text {expected}"
                );
                not_utf8_results += 1;
            }
            Err(error) if function == "numberLength" => {
                // Undefined in the source: static_cast<int> of a value above
                // INT_MAX, which x86 turns into INT_MIN (-214748364.8 after /10).
                assert!(matches!(error, Error::InvalidValue(_)));
                assert!(
                    expected.starts_with("-2.14748e+08") || expected.starts_with("--2.14748e+08"),
                    "{function} {args}: refused, C++ gave {expected}"
                );
                undefined_in_source += 1;
            }
            Err(error) => {
                compared += 1;
                let actual = err(&error);
                if actual != expected {
                    mismatches.push(format!(
                        "{function}\t{args}\n  C++:  {expected}\n  Rust: {actual}"
                    ));
                }
            }
        }
    }
    assert!(
        mismatches.is_empty(),
        "{} mismatches:\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
    // Every record is accounted for; the counts pin the classification.
    assert_eq!(
        compared
            + not_utf8_inputs
            + not_utf8_results
            + undefined_in_source
            + skipped_nonfinite_values,
        5068
    );
    assert_eq!(skipped_nonfinite_values, 9);
    eprintln!(
        "compared {compared}, non-UTF-8 inputs {not_utf8_inputs}, non-UTF-8 results {not_utf8_results}, undefined in source {undefined_in_source}"
    );
    assert!(compared > 4700, "only {compared} records compared");
}

// --------------------------------------------------------------------------
// StringUtils_test.cpp at core bc9cc12, section by section (tier 3)
// --------------------------------------------------------------------------

#[test]
fn class_test_skip_whitespace() {
    // StringUtils_test.cpp:29-54, for each of the four whitespace characters.
    for ws in ['\t', '\r', '\n', ' '] {
        let x16 = "xxxxxxxxxxxxxxxx";
        let s16 = "                ";
        let at1 = format!("0 2  3456789101112{x16}").replace(' ', &ws.to_string());
        assert_eq!(skip_whitespace(&at1), 0);
        assert_eq!(skip_whitespace(&at1[1..]), 1);
        assert_eq!(skip_whitespace(&at1[2..]), 0);
        assert_eq!(skip_whitespace(&at1[3..]), 2);
        let at2 = format!("{s16}{s16}1{x16}").replace(' ', &ws.to_string());
        assert_eq!(skip_whitespace(&at2), 32);
        assert_eq!(skip_whitespace(&at2[2..]), 30);
        let short = "0 2  34".replace(' ', &ws.to_string());
        assert_eq!(skip_whitespace(&short), 0);
        assert_eq!(skip_whitespace(&short[1..]), 1);
        assert_eq!(skip_whitespace(&short[2..]), 0);
        assert_eq!(skip_whitespace(&short[3..]), 2);
    }
}

#[test]
fn class_test_skip_non_whitespace() {
    // StringUtils_test.cpp:56-90.
    for ws in ['\t', '\r', '\n', ' '] {
        let x16 = "xxxxxxxxxxxxxxxx";
        let at1 = format!("0 2  3456789101112{x16}").replace(' ', &ws.to_string());
        assert_eq!(skip_non_whitespace(&at1), 1);
        assert_eq!(skip_non_whitespace(&at1[1..]), 0);
        assert_eq!(skip_non_whitespace(&at1[2..]), 1);
        assert_eq!(skip_non_whitespace(&at1[3..]), 0);
        assert_eq!(skip_non_whitespace(&at1[5..]), 13 + 16);
        let at2 = format!("{x16}{x16} {x16}").replace(' ', &ws.to_string());
        assert_eq!(skip_non_whitespace(&at2), 32);
        assert_eq!(skip_non_whitespace(&at2[31..]), 1);
        assert_eq!(skip_non_whitespace(&at2[33..]), 16);
        let short = "0 2  34".replace(' ', &ws.to_string());
        assert_eq!(skip_non_whitespace(&short), 1);
        assert_eq!(skip_non_whitespace(&short[1..]), 0);
        assert_eq!(skip_non_whitespace(&short[2..]), 1);
        assert_eq!(skip_non_whitespace(&short[3..]), 0);
        assert_eq!(skip_non_whitespace(&short[5..]), 2);
    }
}

#[test]
fn class_test_to_int32_and_to_int64() {
    // StringUtils_test.cpp:263-305.
    assert_eq!(to_int32("2147483647").unwrap(), 2_147_483_647);
    assert_eq!(to_int32("  2147483647").unwrap(), 2_147_483_647);
    assert_eq!(to_int32("2147483647 ").unwrap(), 2_147_483_647);
    assert_eq!(to_int32("   2147483647  ").unwrap(), 2_147_483_647);
    for bad in [
        "2147483648",
        "-2147483649",
        "1234  moreText",
        " 1234 911.0",
        " abc ",
        " 123.45 ",
    ] {
        assert!(
            matches!(to_int32(bad), Err(Error::InvalidValue(_))),
            "{bad}"
        );
    }
    assert_eq!(to_int64("9223372036854775807").unwrap(), i64::MAX);
    assert_eq!(to_int64("  9223372036854775807").unwrap(), i64::MAX);
    assert_eq!(to_int64("9223372036854775807 ").unwrap(), i64::MAX);
    assert_eq!(to_int64("   9223372036854775807  ").unwrap(), i64::MAX);
    for bad in [
        "9223372036854775808",
        "-9223372036854775809",
        "1234  moreText",
        " 1234 911.0",
        " abc ",
        " 123.45 ",
    ] {
        assert!(
            matches!(to_int64(bad), Err(Error::InvalidValue(_))),
            "{bad}"
        );
    }
}

#[test]
fn class_test_to_float_and_to_double() {
    // StringUtils_test.cpp:307-347; TEST_REAL_SIMILAR's default relative
    // tolerance is 1e-5.
    let similar = |a: f64, b: f64| ((a - b) / b).abs() <= 1e-5;
    for s in ["1234.45", "  1234.45", "1234.45 ", "   1234.45  "] {
        assert!(similar(f64::from(to_float(s).unwrap()), 1234.45));
        assert!(similar(to_double(s).unwrap(), 1234.45));
    }
    for bad in ["1234.45  moreText", " 1234.45 911.0", " abc "] {
        assert!(to_float(bad).is_err(), "{bad}");
        assert!(to_double(bad).is_err(), "{bad}");
    }
    assert!(to_float("1e-40").unwrap() > 0.0);
    assert!(to_float("1e40").is_err());
    assert_eq!(
        to_double("2.17388884170148e-321").unwrap(),
        2.17388884170148e-321
    );
    assert!(to_double("4.9e-324").unwrap() > 0.0);
    assert!(to_double("1e400").is_err());
}

#[test]
fn class_test_extract_double() {
    // StringUtils_test.cpp:349-424.
    let similar = |a: f64, b: f64| ((a - b) / b).abs() <= 1e-5;
    let (d, n) = extract_double("12345.45  ").unwrap();
    assert!(similar(d, 12345.45));
    assert_eq!(n, 8);
    let (d, n) = extract_double("+1234.45!").unwrap();
    assert!(similar(d, 1234.45));
    assert_eq!(n, 8);
    assert_eq!(extract_double("  -123.45"), None);
    assert_eq!(extract_double("15.0e6"), Some((15.0e6, 6)));
    let s = "-5.0\t9.1";
    assert_eq!(extract_double(s), Some((-5.0, 4)));
    let (d, n) = extract_double(&s[5..]).unwrap();
    assert!(similar(d, 9.1));
    assert_eq!(5 + n, 8);
    assert_eq!(extract_double("15.0e6 x"), Some((15.0e6, 6)));
    assert_eq!(extract_double("16e6!"), Some((16e6, 4)));
    assert_eq!(extract_double("!noNumber"), None);
    assert_eq!(
        extract_double("2.17388884170148e-321 x"),
        Some((2.17388884170148e-321, 21))
    );
}

#[test]
fn class_test_extract_int() {
    // StringUtils_test.cpp:426-479.
    assert_eq!(extract_int("12345  "), Some((12345, 5)));
    assert_eq!(extract_int("+1234!"), Some((1234, 5)));
    assert_eq!(extract_int("  -123"), None);
    assert_eq!(extract_int(&"y15++"[1..]), Some((15, 2)));
    let s = "-5\t9";
    assert_eq!(extract_int(s), Some((-5, 2)));
    assert_eq!(extract_int(&s[3..]), Some((9, 1)));
    assert_eq!(extract_int("!noNumber"), None);
}

#[test]
fn class_test_non_finite_values_round_trip() {
    // StringUtils_test.cpp:530-574 ([EXTRA] non-finite values round-trip).
    assert_eq!(to_str(f64::INFINITY), "inf");
    assert_eq!(to_str(f64::NEG_INFINITY), "-inf");
    assert_eq!(to_str_f32(f32::INFINITY), "inf");
    assert_eq!(to_str_f32(f32::NEG_INFINITY), "-inf");
    assert_eq!(to_str(f64::NAN), "NaN");
    assert!(to_double(&to_str(f64::INFINITY)).unwrap() == f64::INFINITY);
    assert!(to_double(&to_str(f64::NEG_INFINITY)).unwrap() == f64::NEG_INFINITY);
    assert!(to_float(&to_str_f32(f32::INFINITY)).unwrap() == f32::INFINITY);
    assert!(to_float(&to_str_f32(f32::NEG_INFINITY)).unwrap() < 0.0);
    assert!(to_double(&to_str(f64::NAN)).unwrap().is_nan());
    assert_eq!(to_str(5.0), "5.0");
    assert_eq!(to_str(-0.5), "-0.5");
    assert_eq!(to_str(0.0), "0.0");
    // The DataValue list half: the list text split on ',' reads back.
    let text = ["1.5", "inf", "-inf", "NaN", "2.5"]
        .iter()
        .map(|t| to_double(t).unwrap())
        .map(to_str)
        .collect::<Vec<_>>()
        .join(", ");
    let parts = su::split(&text, ',', false).unwrap().parts;
    assert_eq!(parts.len(), 5);
    let back: Vec<f64> = parts.iter().map(|p| to_double(p).unwrap()).collect();
    assert!((back[0] - 1.5).abs() < 1e-9);
    assert!(back[1].is_infinite() && back[1] > 0.0);
    assert!(back[2].is_infinite() && back[2] < 0.0);
    assert!(back[3].is_nan());
    assert!((back[4] - 2.5).abs() < 1e-9);
}

#[test]
fn native_boundaries() {
    // A byte count that cuts a multi-byte character is refused, not cut.
    assert!(su::prefix("\u{e4}b", 1).is_err());
    assert!(su::chop("a\u{e4}", 1).is_err());
    // C++ chars are one byte: a non-ASCII quote, fill or splitter is refused.
    let mut s = String::from("a");
    assert!(su::quote(&mut s, '\u{e4}', QuotingMethod::Escape).is_err());
    assert!(su::fill_left(&mut s, '\u{e4}', 4).is_err());
    assert!(su::split("a,b", '\u{e4}', false).is_err());
    assert_eq!(s, "a");
    // Ceilings are checked before allocating and leave the input unchanged.
    assert!(su::fill_left(&mut s, 'x', u32::MAX).is_err());
    assert_eq!(s, "a");
    assert!(su::random(u32::MAX).is_err());
    // numberLength refuses what the source cannot finish.
    assert!(number_length(f64::NAN, 5).is_err());
    assert!(number_length(f64::NEG_INFINITY, 5).is_err());
}

// --------------------------------------------------------------------------
// String_test.cpp at core bc9cc12: the class test of the former String class,
// which now tests the StringUtils free functions directly (tier 3). The
// std::string constructor sections (String_test.cpp:40-89) test the standard
// library and have no port counterpart; the long double sections
// (:163-186, :657-659, :685-696) have no Rust type to test.
// --------------------------------------------------------------------------

const AMINO: &str = "ACDEFGHIKLMNPQRSTVWY";

fn relative_similar(a: f64, b: f64) -> bool {
    // TEST_REAL_SIMILAR's default relative tolerance.
    ((a - b) / b).abs() <= 1e-5
}

#[test]
fn string_test_integer_to_str() {
    // String_test.cpp:94-124.
    assert_eq!((-17_i32).to_str_text(), "-17");
    assert_eq!(17_u32.to_str_text(), "17");
    assert_eq!((-17_i64).to_str_text(), "-17");
    assert_eq!(17_u64.to_str_text(), "17");
    assert_eq!((-17_i16).to_str_text(), "-17");
    assert_eq!(17_u16.to_str_text(), "17");
    assert_eq!(12_345_678_u64.to_str_text(), "12345678");
    assert_eq!((-12_345_678_i64).to_str_text(), "-12345678");
}

#[test]
// The float literals are the class test's own spellings (17.0123456f,
// 50254.199219f), which carry more digits than an f32 holds.
#[allow(clippy::excessive_precision)]
fn string_test_float_to_str() {
    // String_test.cpp:126-140.
    assert_eq!(to_str_f32(17.012_345_6), "17.012346");
    assert_eq!(to_str_f32_low_p(17.0123), "17.012");
    assert_eq!(to_str_f32(50_254.199_219), "5.02542e04");
    assert_eq!(to_str_f32(1_000_000.0), "1.0e06");
    assert_eq!(to_str_f32(0.000_001), "1.0e-06");
    let denorm = f32::MIN_POSITIVE / 10.0;
    assert_eq!(to_str_f32(denorm), "1.175495e-39");
    assert_eq!(to_str_f32_low_p(denorm), "1.175e-39");
    assert_eq!(to_str_f32(f32::NAN), "NaN");
    assert_eq!(to_str_f32_low_p(f32::NAN), "NaN");
}

#[test]
fn string_test_double_to_str() {
    // String_test.cpp:142-161.
    assert_eq!(to_str(17.012345), "17.012345");
    assert_eq!(to_str_low_p(17.012345), "17.012");
    assert_eq!(to_str(1_000_000.0), "1.0e06");
    assert_eq!(to_str(-1_000_000.0), "-1.0e06");
    assert_eq!(to_str(10_000.0), "1.0e04");
    assert_eq!(to_str(0.000_001), "1.0e-06");
    assert_eq!(to_str(1.5e10), "1.5e10");
    let denorm = f64::MIN_POSITIVE / 10.0;
    assert_eq!(to_str(denorm), "2.225073858507203e-309");
    assert_eq!(to_str_low_p(denorm), "2.225e-309");
    assert_eq!(to_str(f64::NAN), "NaN");
    assert_eq!(to_str_low_p(f64::NAN), "NaN");
}

#[test]
fn string_test_data_value_to_str() {
    // String_test.cpp:188-195.
    let value = MetaValue::try_from(17.012345).unwrap();
    assert_eq!(value.to_text(true).unwrap(), "17.012345");
    assert_eq!(value.to_text(false).unwrap(), "17.012");
    let list = MetaValue::try_from(vec![17.012345, 2.0]).unwrap();
    assert_eq!(list.to_text(true).unwrap(), "[17.012345, 2.0]");
    assert_eq!(list.to_text(false).unwrap(), "[17.012, 2.0]");
    assert_eq!(MetaValue::from("bla").to_text(true).unwrap(), "bla");
    assert_eq!(MetaValue::from(4711_i64).to_text(true).unwrap(), "4711");
}

#[test]
fn string_test_number_length_number_random() {
    // String_test.cpp:201-223.
    assert_eq!(number_length(12345678.9123, 11).unwrap(), "12345678.91");
    assert_eq!(number_length(-12345678.9123, 11).unwrap(), "-12345678.9");
    assert_eq!(number_length(12345678.9123, 10).unwrap(), "12345678.9");
    assert_eq!(number_length(-12345678.9123, 10).unwrap(), "-1234.5e04");
    assert_eq!(number_length(12345678.9123, 9).unwrap(), "1234.5e04");
    assert_eq!(number_length(-12345678.9123, 9).unwrap(), "-123.4e05");
    for (digits, expected) in [
        (0, "123"),
        (1, "123.1"),
        (2, "123.12"),
        (3, "123.123"),
        (4, "123.1234"),
        (5, "123.12340"),
    ] {
        assert_eq!(number(123.1234, digits).unwrap(), expected);
    }
    assert_eq!(number(0.0, 5).unwrap(), "0.00000");
    assert_eq!(su::random(10).unwrap().len(), 10);
}

#[test]
fn string_test_predicates() {
    // String_test.cpp:231-262; hasPrefix, hasSuffix, hasSubstring and has map
    // onto str::starts_with, ends_with and contains (see the module table).
    let s = AMINO;
    assert!(s.starts_with(""));
    assert!(s.starts_with("ACDEF"));
    assert!(s.starts_with(AMINO));
    assert!(!s.starts_with("ABCDEF"));
    assert!(!s.starts_with("ACDEFGHIKLMNPQRSTVWYACDEF"));
    assert!(s.ends_with(""));
    assert!(s.ends_with("TVWY"));
    assert!(s.ends_with(AMINO));
    assert!(!s.ends_with("WXYZ"));
    assert!(!s.ends_with("ACDEFACDEFGHIKLMNPQRSTVWY"));
    assert!(s.contains(""));
    assert!(s.contains("GHIKLM"));
    assert!(s.contains(AMINO));
    assert!(!s.contains("MLKIGH"));
    assert!(!s.contains("ACDEFGHIKLMNPQRSTVWYACDEF"));
    assert!(s.as_bytes().contains(&b'A'));
    assert!(!s.as_bytes().contains(&b'O'));
}

fn is_overflow(result: openms::Result<String>) -> bool {
    matches!(result, Err(Error::InvalidRange(ref m)) if m.contains("too large"))
}

fn is_underflow(result: openms::Result<String>) -> bool {
    matches!(result, Err(Error::InvalidRange(ref m)) if m.contains("too small"))
}

#[test]
fn string_test_accessors() {
    // String_test.cpp:264-333.
    let s = AMINO;
    assert_eq!(su::prefix_signed(s, 4).unwrap(), "ACDE");
    assert_eq!(su::prefix_signed(s, 0).unwrap(), "");
    assert!(is_overflow(su::prefix_signed(s, s.len() as i32 + 1)));
    assert!(is_underflow(su::prefix_signed(s, -1)));
    assert_eq!(su::suffix_signed(s, 4).unwrap(), "TVWY");
    assert_eq!(su::suffix_signed(s, 0).unwrap(), "");
    assert!(is_overflow(su::suffix_signed(s, s.len() as i32 + 1)));
    assert!(is_underflow(su::suffix_signed(s, -1)));
    assert_eq!(su::prefix(s, 4).unwrap(), "ACDE");
    assert_eq!(su::prefix(s, 0).unwrap(), "");
    assert!(is_overflow(su::prefix(s, s.len() + 1)));
    assert_eq!(su::suffix(s, 4).unwrap(), "TVWY");
    assert_eq!(su::suffix(s, 0).unwrap(), "");
    assert!(is_overflow(su::suffix(s, s.len() + 1)));
    assert_eq!(su::prefix_until(s, 'F'), "ACDE");
    assert_eq!(su::prefix_until(s, 'A'), "");
    assert_eq!(su::prefix_until(s, 'Z'), s);
    assert_eq!(su::suffix_after(s, 'S'), "TVWY");
    assert_eq!(su::suffix_after(s, 'Y'), "");
    assert_eq!(su::suffix_after(s, 'Z'), s);
    let s2 = "abcdef";
    for (position, length, expected) in [
        (0, Some(4), "abcd"),
        (1, Some(1), "b"),
        (1, Some(3), "bcd"),
        (0, Some(6), "abcdef"),
        (5, Some(1), "f"),
        (6, Some(1), ""),
        (0, Some(7), "abcdef"),
        (0, None, "abcdef"),
        (1, None, "bcdef"),
        (5, None, "f"),
        (6, None, ""),
    ] {
        assert_eq!(su::substr(s2, position, length).unwrap(), expected);
    }
    for (count, expected) in [
        (0, "abcdef"),
        (1, "abcde"),
        (2, "abcd"),
        (3, "abc"),
        (4, "ab"),
        (5, "a"),
        (6, ""),
        (9, ""),
        (usize::MAX, ""),
    ] {
        assert_eq!(su::chop(s2, count).unwrap(), expected);
    }
}

#[test]
fn string_test_reverse_trim_simplify_fill() {
    // String_test.cpp:335-438.
    let mut r = String::from(AMINO);
    su::reverse(&mut r).unwrap();
    assert_eq!(r, "YWVTSRQPNMLKIHGFEDCA");
    r.clear();
    su::reverse(&mut r).unwrap();
    assert_eq!(r, "");

    let mut t = String::from("\n\r\t test \n\r\t");
    su::trim(&mut t);
    assert_eq!(t, "test");
    su::trim(&mut t);
    assert_eq!(t, "test");
    for (input, expected) in [("", ""), (" t", "t"), ("t ", "t"), ("\t\r\n ", "")] {
        let mut t = String::from(input);
        su::trim(&mut t);
        assert_eq!(t, expected);
    }

    let mut si = String::from("\n\r\t te\tst \n\r\t");
    su::simplify(&mut si);
    assert_eq!(si, " te st ");
    su::simplify(&mut si);
    assert_eq!(si, " te st ");
    for (input, expected) in [("", ""), (" t", " t"), ("t ", "t "), ("\t\r\n ", " ")] {
        let mut si = String::from(input);
        su::simplify(&mut si);
        assert_eq!(si, expected);
    }

    let mut fl = String::from("TEST");
    su::fill_left(&mut fl, 'x', 4).unwrap();
    assert_eq!(fl, "TEST");
    su::fill_left(&mut fl, 'y', 5).unwrap();
    assert_eq!(fl, "yTEST");
    su::fill_left(&mut fl, 'z', 7).unwrap();
    assert_eq!(fl, "zzyTEST");
    let mut fr = String::from("TEST");
    su::fill_right(&mut fr, 'x', 4).unwrap();
    assert_eq!(fr, "TEST");
    su::fill_right(&mut fr, 'y', 5).unwrap();
    assert_eq!(fr, "TESTy");
    su::fill_right(&mut fr, 'z', 7).unwrap();
    assert_eq!(fr, "TESTyzz");
}

#[test]
fn string_test_quote_unquote() {
    // String_test.cpp:364-398.
    let mut q = String::new();
    su::quote(&mut q, '\'', QuotingMethod::None).unwrap();
    assert_eq!(q, "''");
    su::quote(&mut q, '\'', QuotingMethod::Escape).unwrap();
    assert_eq!(q, "'\\'\\''");
    let mut q = String::from("ab\"cd\\ef");
    su::quote(&mut q, '"', QuotingMethod::None).unwrap();
    assert_eq!(q, "\"ab\"cd\\ef\"");
    su::quote(&mut q, '"', QuotingMethod::Escape).unwrap();
    assert_eq!(q, "\"\\\"ab\\\"cd\\\\ef\\\"\"");
    let mut q = String::from("ab\"cd\\ef");
    su::quote(&mut q, '"', QuotingMethod::Double).unwrap();
    assert_eq!(q, "\"ab\"\"cd\\ef\"");

    let mut u = String::new();
    assert!(matches!(
        su::unquote(&mut u, '"', QuotingMethod::default()),
        Err(Error::InvalidValue(_))
    ));
    let mut u = String::from("''");
    su::unquote(&mut u, '\'', QuotingMethod::None).unwrap();
    assert_eq!(u, "");
    let mut u = String::from("'\\'\\''");
    su::unquote(&mut u, '\'', QuotingMethod::Escape).unwrap();
    assert_eq!(u, "''");
    let mut u = String::from(r#""ab"cd\ef""#);
    su::unquote(&mut u, '"', QuotingMethod::None).unwrap();
    assert_eq!(u, "ab\"cd\\ef");
    let mut u = String::from(r#""\"ab\"cd\\ef\"""#);
    su::unquote(&mut u, '"', QuotingMethod::Escape).unwrap();
    assert_eq!(u, "\"ab\"cd\\ef\"");
    let mut u = String::from(r#""ab""cd\ef""#);
    su::unquote(&mut u, '"', QuotingMethod::Double).unwrap();
    assert_eq!(u, "ab\"cd\\ef");
}

fn conversion_message<T: std::fmt::Debug>(result: openms::Result<T>) -> String {
    match result {
        Err(Error::InvalidValue(message)) => message,
        other => panic!("expected a ConversionError, got {other:?}"),
    }
}

#[test]
fn string_test_converters() {
    // String_test.cpp:444-495.
    assert_eq!(to_int32("123").unwrap(), 123);
    assert_eq!(to_int32("-123").unwrap(), -123);
    assert_eq!(to_int32("  -123  ").unwrap(), -123);
    conversion_message(to_int32("524 starts with an int"));
    assert_eq!(
        conversion_message(to_int32("not an int")),
        "Could not convert string 'not an int' to an integer value"
    );

    assert!(relative_similar(
        f64::from(to_float("123.456").unwrap()),
        f64::from(123.456_f32)
    ));
    assert!(relative_similar(
        f64::from(to_float("-123.456").unwrap()),
        f64::from(-123.456_f32)
    ));
    assert_eq!(to_str_f32(to_float("73629.9").unwrap()), "7.36299e04");
    assert!(to_float("nan").unwrap().is_nan());
    assert!(to_float("NaN").unwrap().is_nan());
    assert_eq!(
        conversion_message(to_float("not a number")),
        "Could not convert string 'not a number' to a float value"
    );

    assert!(relative_similar(to_double("123.456").unwrap(), 123.456));
    assert!(relative_similar(
        to_double("-123.4567890123").unwrap(),
        -123.4567890123
    ));
    assert_eq!(
        to_str(to_double("73629.980123").unwrap()),
        "7.3629980123e04"
    );
    assert!(to_double("nan").unwrap().is_nan());
    assert!(to_double("NaN").unwrap().is_nan());
    assert_eq!(
        conversion_message(to_double("not a number")),
        "Could not convert string 'not a number' to a double value"
    );
}

#[test]
fn string_test_case_substitute_remove() {
    // String_test.cpp:497-619.
    let mut tu = String::from("test45%#.,");
    su::to_upper(&mut tu);
    assert_eq!(tu, "TEST45%#.,");
    let mut empty = String::new();
    su::to_upper(&mut empty);
    assert_eq!(empty, "");
    let mut tl = String::from("TEST45%#.,");
    su::to_lower(&mut tl);
    assert_eq!(tl, "test45%#.,");
    su::to_lower(&mut empty);
    assert_eq!(empty, "");
    for (input, expected) in [("test45%#.,", "Test45%#.,"), (" ", " "), ("", "")] {
        let mut fu = String::from(input);
        su::first_to_upper(&mut fu);
        assert_eq!(fu, expected);
    }

    let mut sub = String::from("abcdefg");
    su::substitute_char(&mut sub, 'a', 'x');
    assert_eq!(sub, "xbcdefg");
    su::substitute_char(&mut sub, 'g', 'y');
    assert_eq!(sub, "xbcdefy");
    su::substitute_char(&mut sub, 'c', '-');
    assert_eq!(sub, "xb-defy");
    let mut dots = String::from(".....");
    su::substitute_char(&mut dots, '.', ',');
    assert_eq!(dots, ",,,,,");

    let mut sub = String::from("abcdefg");
    for (from, to, expected) in [
        ("a", "x", "xbcdefg"),
        ("bcd", "y", "xyefg"),
        ("fg", "", "xye"),
        ("e", "z!", "xyz!"),
        ("u", "blblblblbl", "xyz!"),
        ("", "blblblblbl", "xyz!"),
    ] {
        su::substitute(&mut sub, from, to).unwrap();
        assert_eq!(sub, expected);
    }
    let mut sub = String::from("abcdefgabcdefgabcdefgab");
    su::substitute(&mut sub, "ab", "x").unwrap();
    assert_eq!(sub, "xcdefgxcdefgxcdefgx");
    su::substitute(&mut sub, "x", "").unwrap();
    assert_eq!(sub, "cdefgcdefgcdefg");

    let mut rm = String::from("abcabc");
    su::remove(&mut rm, 'a');
    assert_eq!(rm, "bcbc");
    su::remove(&mut rm, 'c');
    assert_eq!(rm, "bb");
    su::remove(&mut rm, 'b');
    assert_eq!(rm, "");

    let mut elc = String::from("/");
    for (end, expected) in [('/', "/"), ('\\', "/\\"), ('\\', "/\\"), ('/', "/\\/")] {
        su::ensure_last_char(&mut elc, end);
        assert_eq!(elc, expected);
    }

    for (input, expected) in [
        ("", ""),
        (" \t \n ", ""),
        ("test", "test"),
        ("\n\r\t test \n\r\t", "test"),
        ("\n\r\t t\ne \ts\rt \n\r\t", "test"),
    ] {
        let mut rw = String::from(input);
        su::remove_whitespaces(&mut rw);
        assert_eq!(rw, expected);
    }
    // The source shuffles 100 spaces into 16 KiB + 1 'x' with a randomly
    // seeded engine; a fixed spread of the same 100 spaces stands in for it.
    let test = "x".repeat(16 * 1024 + 1);
    let mut rw = String::with_capacity(test.len() + 100);
    let mut spaces = 0;
    for (index, c) in test.chars().enumerate() {
        if index % 164 == 0 && spaces < 100 {
            rw.push(' ');
            spaces += 1;
        }
        rw.push(c);
    }
    assert_eq!(spaces, 100);
    su::remove_whitespaces(&mut rw);
    assert_eq!(rw, test);
}

#[test]
fn string_test_operators() {
    // String_test.cpp:621-683: std::string + number and std::string += number.
    fn plus(value: &dyn ToStr) -> String {
        let mut text = String::from("test");
        value.append_to_str(&mut text);
        text
    }
    assert_eq!(plus(&4_i32), "test4");
    assert_eq!(plus(&4_u32), "test4");
    assert_eq!(plus(&4_i16), "test4");
    assert_eq!(plus(&4_u16), "test4");
    assert_eq!(plus(&4_i64), "test4");
    assert_eq!(plus(&4_u64), "test4");
    assert_eq!(plus(&4.0_f32), "test4.0");
    assert_eq!(plus(&4.0_f64), "test4.0");
    assert_eq!(plus(&7_i32), "test7");
    assert_eq!(plus(&7_u32), "test7");
    assert_eq!(plus(&7.4_f32), "test7.4");
    assert_eq!(plus(&7.4_f64), "test7.4");
}

#[test]
fn string_test_split_and_concatenate() {
    // String_test.cpp:698-811.
    let r = su::split(";1;2;3;4;5;", ';', false).unwrap();
    assert!(r.split);
    assert_eq!(r.parts.len(), 7);
    assert_eq!(r.parts[0], "");
    assert_eq!(r.parts[1], "1");
    assert_eq!(r.parts[6], "");
    let r = su::split("1;2;3;4;5", ';', false).unwrap();
    assert!(r.split);
    assert_eq!(r.parts.len(), 5);
    assert_eq!(su::split("", ',', false).unwrap(), Split::default());
    let r = su::split(";", ';', false).unwrap();
    assert!(r.split);
    assert_eq!(r.parts.len(), 2);
    let r = su::split("nodelim", ';', false).unwrap();
    assert!(!r.split);
    assert_eq!(r.parts.len(), 1);
    let r = su::split(" \"hello\", world, 23.3", ',', true).unwrap();
    assert!(r.split);
    assert_eq!(r.parts, ["hello", "world", "23.3"]);
    let r = su::split(r#" "hello", " donot,splitthis ", "23.4 " "#, ',', true).unwrap();
    assert!(r.split);
    assert_eq!(r.parts.len(), 3);
    assert_eq!(r.parts[1], " donot,splitthis ");
    assert!(matches!(
        su::split(
            r#" "first", "seconds"<thisshouldnotbehere>, third"#,
            ',',
            true
        ),
        Err(Error::InvalidValue(_))
    ));

    let sp2 = "abcdebcfghbc";
    let r = su::split_str(sp2, "bc").unwrap();
    assert!(r.split);
    assert_eq!(r.parts, ["a", "de", "fgh", ""]);
    let r = su::split_str(sp2, "xy").unwrap();
    assert!(!r.split);
    assert_eq!(r.parts.len(), 1);
    let r = su::split_str(sp2, "").unwrap();
    assert!(r.split);
    assert_eq!(r.parts.len(), sp2.len());
    let r = su::split_str("", ",").unwrap();
    assert!(!r.split);
    assert!(r.parts.is_empty());

    let r = su::split_quoted(sp2, "bc", '"', QuotingMethod::default()).unwrap();
    assert!(r.split);
    assert_eq!(r.parts.len(), 4);
    let r = su::split_quoted(r#""a,b,c","d,\",f","""#, ",", '"', QuotingMethod::Escape).unwrap();
    assert!(r.split);
    assert_eq!(r.parts, ["\"a,b,c\"", "\"d,\\\",f\"", "\"\""]);
    assert!(matches!(
        su::split_quoted(r#""a,"b""#, ",", '"', QuotingMethod::Escape),
        Err(Error::InvalidValue(_))
    ));
    let r = su::split_quoted(r#""ab"___"cd""ef""#, "___", '"', QuotingMethod::Double).unwrap();
    assert!(r.split);
    assert_eq!(r.parts, ["\"ab\"", "\"cd\"\"ef\""]);

    let parts = su::split("1;2;3;4;5", ';', false).unwrap().parts;
    assert_eq!(su::concatenate(&parts, "g").unwrap(), "1g2g3g4g5");
    assert_eq!(su::concatenate(&parts, "").unwrap(), "12345");
    let parts = su::split("", ';', false).unwrap().parts;
    assert_eq!(su::concatenate(&parts, "").unwrap(), "");
}
