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
//! `src/tests/class_tests/openms/source/StringUtils_test.cpp` at core
//! `bc9cc12` (tier 3); the oracle is tier 1.

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
