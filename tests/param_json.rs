// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! `ParamJSONFile::load` (`src/format/param_json.rs`): the upstream class test
//! `ParamJSONFile_test.cpp` (core pin `bc9cc12`), transcribed, and the rules
//! of `ParamJSONFile.cpp:38-167` and of the `nlohmann::json` 3.11.2 accessors
//! it calls, each on a parameter tree built for the case. The TOPP tools'
//! JSON `-ini` files are compared with the Release build in
//! `tests/topp_json_ini.rs`.
#![cfg(feature = "paramxml")]

use openms::Error;
use openms::format::param_json::{self, MAX_JSON_BYTES};
use openms::param::{Param, ParamValue};
use openms::system::file::TempDir;
use std::fs;

fn tags(tags: &[&str]) -> Vec<String> {
    tags.iter().map(|tag| (*tag).to_owned()).collect()
}

/// The parameter tree of `ParamJSONFile_test.cpp:26-49`.
fn upstream_param() -> Param {
    let mut param = Param::new();
    param
        .set_value("test:1:value", ParamValue::Integer(1), "description", &[])
        .unwrap();
    for (key, value) in [
        ("bool1", "false"),
        ("bool2", "false"),
        ("bool3", "true"),
        ("bool4", "true"),
    ] {
        let key = format!("test:1:{key}");
        param.set_value(&key, value.into(), "", &[]).unwrap();
    }
    let true_false = tags(&["true", "false"]);
    let false_true = tags(&["false", "true"]);
    param
        .set_valid_strings("test:1:bool1", &true_false)
        .unwrap();
    param
        .set_valid_strings("test:1:bool2", &true_false)
        .unwrap();
    param
        .set_valid_strings("test:1:bool3", &false_true)
        .unwrap();
    param
        .set_valid_strings("test:1:bool4", &false_true)
        .unwrap();
    param
        .set_value("test:1:int", ParamValue::Integer(0), "", &[])
        .unwrap();
    param
        .set_value("test:1:double", ParamValue::Float(0.0), "", &[])
        .unwrap();
    param
        .set_value("test:1:string", "".into(), "", &[])
        .unwrap();
    param
        .set_value("test:1:int_list", ParamValue::IntegerList(vec![]), "", &[])
        .unwrap();
    param
        .set_value("test:1:double_list", ParamValue::FloatList(vec![]), "", &[])
        .unwrap();
    param
        .set_value(
            "test:1:string_list",
            ParamValue::StringList(vec![]),
            "",
            &[],
        )
        .unwrap();
    param
        .set_value(
            "test:1:file_output",
            "".into(),
            "some description",
            &tags(&["output file"]),
        )
        .unwrap();
    param
        .set_value(
            "test:1:is_executable_v1",
            "".into(),
            "test is executable tag, giving a string",
            &tags(&["is_executable", "input file"]),
        )
        .unwrap();
    param
        .set_value(
            "test:1:is_executable_v2",
            "".into(),
            "test is executable tag, giving a type: File",
            &tags(&["is_executable", "input file"]),
        )
        .unwrap();
    param
}

fn load(document: &str, param: &mut Param) -> openms::Result<Vec<String>> {
    param_json::load_from_slice(document.as_bytes(), "doc.json", param)
}

fn parse_message(result: openms::Result<Vec<String>>) -> String {
    match result {
        Err(Error::Parse { message, .. }) => message,
        other => panic!("expected a parse error, got {other:?}"),
    }
}

/// `START_SECTION((bool ParamJSONFile::load(...)))`
/// (`ParamJSONFile_test.cpp:21-107`).
#[test]
fn upstream_load() {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let mut param = upstream_param();

    // FileNotFound (`29`): the port's I/O error of kind NotFound.
    match param_json::load("/does/not/exist/FileDoesNotExist.json", &mut param) {
        Err(Error::Io(error)) => assert_eq!(error.kind(), std::io::ErrorKind::NotFound),
        other => panic!("{other:?}"),
    }

    // ParseError (`32-37`).
    let filename = dir.path().join("test.json");
    fs::write(&filename, "not a json").unwrap();
    assert!(matches!(
        param_json::load(&filename, &mut param),
        Err(Error::Parse { .. })
    ));

    // All types (`56-106`).
    fs::write(
        &filename,
        "{\n  \"bool1\": true,\n  \"bool2\": false,\n  \"bool3\": true,\n  \"bool4\": false,\n  \"int\": 5,\n  \"double\": 6.1,\n  \"string\": \"Hello OpenMS\",\n  \"int_list\": [10, 11, 12],\n  \"double_list\": [13.25, 15.125],\n  \"string_list\": [\"SeqAn\", \"rocks\"],\n  \"file_output\": \"/some/made/up/path\",\n  \"is_executable_v1\": \"/some/made/up/path\",\n  \"is_executable_v2\": {\n        \"class\": \"File\",\n        \"path\": \"/some/made/up/path\"\n  }\n}\n",
    )
    .unwrap();
    assert_eq!(
        param_json::load(&filename, &mut param).unwrap(),
        Vec::<String>::new()
    );
    let value = |key: &str| param.value(&format!("test:1:{key}")).unwrap().clone();
    assert!(value("bool1").to_bool().unwrap());
    assert!(!value("bool2").to_bool().unwrap());
    assert!(value("bool3").to_bool().unwrap());
    assert!(!value("bool4").to_bool().unwrap());
    assert_eq!(value("int"), ParamValue::Integer(5));
    assert_eq!(value("double"), ParamValue::Float(6.1));
    assert_eq!(value("string"), ParamValue::String("Hello OpenMS".into()));
    assert_eq!(value("int_list"), ParamValue::IntegerList(vec![10, 11, 12]));
    assert_eq!(
        value("double_list"),
        ParamValue::FloatList(vec![13.25, 15.125])
    );
    assert_eq!(
        value("string_list"),
        ParamValue::StringList(vec!["SeqAn".into(), "rocks".into()])
    );
    // Not asserted upstream; the values the three branches of
    // `ParamJSONFile.cpp:100-123` take.
    assert_eq!(
        value("file_output"),
        ParamValue::String("/some/made/up/path".into())
    );
    assert_eq!(
        value("is_executable_v1"),
        ParamValue::String("/some/made/up/path".into())
    );
    assert_eq!(
        value("is_executable_v2"),
        ParamValue::String("/some/made/up/path".into())
    );
}

/// A string list tagged `input file` is read from the `path` member of a file
/// object (`ParamJSONFile.cpp:135-138`), and nothing else.
#[test]
fn an_input_file_list_comes_from_a_file_object() {
    let files = || {
        let mut param = Param::new();
        param
            .set_value(
                "t:1:in",
                ParamValue::StringList(vec![]),
                "",
                &tags(&["input file"]),
            )
            .unwrap();
        param
    };
    let mut param = files();
    load(
        r#"{"in": {"class": "File", "path": ["a", "b"]}}"#,
        &mut param,
    )
    .unwrap();
    assert_eq!(
        param.value("t:1:in").unwrap(),
        &ParamValue::StringList(vec!["a".into(), "b".into()])
    );
    // `Param::setValue` replaces the tags with its empty ones
    // (`Param.cpp:466-469`, `408-411`), so the loaded entry is no longer an
    // input file.
    assert!(param.tags("t:1:in").unwrap().is_empty());
    assert_eq!(
        parse_message(load(r#"{"in": ["a"]}"#, &mut files())),
        "[json.exception.type_error.305] cannot use operator[] with a string argument with array"
    );
}

/// `get<std::vector<double>>()`: numbers only; a boolean is refused, as for
/// `get<double>()` (json.hpp `get_arithmetic_value`, 4605-4640).
#[test]
fn a_float_list_takes_numbers_only() {
    let mut param = Param::new();
    param
        .set_value("t:1:d", ParamValue::FloatList(vec![]), "", &[])
        .unwrap();
    load(r#"{"d": [1, -2, 0.5, 18446744073709551615]}"#, &mut param).unwrap();
    assert_eq!(
        param.value("t:1:d").unwrap(),
        &ParamValue::FloatList(vec![1.0, -2.0, 0.5, 18446744073709551615.0])
    );
    assert_eq!(
        parse_message(load(r#"{"d": [true]}"#, &mut param)),
        "[json.exception.type_error.302] type must be number, but is boolean"
    );
    assert_eq!(
        parse_message(load(r#"{"d": 1.5}"#, &mut param)),
        "[json.exception.type_error.302] type must be array, but is number"
    );
}

/// `get<int>()` for an integer list element takes a boolean as 0 or 1 and
/// casts the rest to 32 bits (json.hpp 4891-4934), where `get<int64_t>()`
/// refuses a boolean.
#[test]
fn integer_conversions_follow_the_accessor_the_source_calls() {
    let mut param = Param::new();
    param
        .set_value("t:1:i", ParamValue::Integer(0), "", &[])
        .unwrap();
    param
        .set_value("t:1:l", ParamValue::IntegerList(vec![]), "", &[])
        .unwrap();
    load(
        r#"{"i": -7.9, "l": [true, false, 2.9, -1.5, 4294967298, -4294967297, 1e10, -1e10, -0]}"#,
        &mut param,
    )
    .unwrap();
    assert_eq!(param.value("t:1:i").unwrap(), &ParamValue::Integer(-7));
    assert_eq!(
        param.value("t:1:l").unwrap(),
        &ParamValue::IntegerList(vec![1, 0, 2, -1, 2, -1, i32::MIN, i32::MIN, 0])
    );
    for (document, value) in [
        (r#"{"i": 1e30}"#, i64::MIN),
        (r#"{"i": -1e30}"#, i64::MIN),
        (r#"{"i": 18446744073709551615}"#, -1),
        (r#"{"i": 9223372036854775808}"#, i64::MIN),
        (r#"{"i": -9223372036854775808}"#, i64::MIN),
    ] {
        load(document, &mut param).unwrap();
        assert_eq!(
            param.value("t:1:i").unwrap(),
            &ParamValue::Integer(value),
            "{document}"
        );
    }
    assert_eq!(
        parse_message(load(r#"{"i": true}"#, &mut param)),
        "[json.exception.type_error.302] type must be number, but is boolean"
    );
    assert_eq!(
        parse_message(load(r#"{"l": [[1]]}"#, &mut param)),
        "[json.exception.type_error.302] type must be number, but is array"
    );
}

/// `-0` is the integer 0 in `nlohmann::json` and `-0.0` the float -0.0.
#[test]
fn negative_zero_is_zero_unless_written_as_a_float() {
    let mut param = Param::new();
    param
        .set_value("t:1:d", ParamValue::Float(1.0), "", &[])
        .unwrap();
    param.set_value("t:1:s", "".into(), "", &[]).unwrap();
    for (document, negative) in [
        (r#"{"d": -0}"#, false),
        (r#"{"d":-0 }"#, false),
        (r#"{"d": -0.0}"#, true),
        (r#"{"d": -0e0}"#, true),
        (r#"{"d": 1e-0, "s": "-0"}"#, false),
    ] {
        load(document, &mut param).unwrap();
        let ParamValue::Float(value) = param.value("t:1:d").unwrap() else {
            panic!()
        };
        assert_eq!(value.is_sign_negative(), negative, "{document}");
    }
    assert_eq!(
        param.value("t:1:s").unwrap(),
        &ParamValue::String("-0".into())
    );
}

/// The walk (`ParamJSONFile.cpp:64-160`): `__` for `:`, sections, a `null`
/// skipped, the first failure in the byte order of the names, the instance
/// always `1`, and the root an object.
#[test]
fn names_sections_and_the_order_of_failures() {
    let mut param = Param::new();
    param
        .set_value("tool:version", "1.0.0".into(), "", &[])
        .unwrap();
    param
        .set_value("tool:1:a:b", ParamValue::Integer(0), "", &[])
        .unwrap();
    param
        .set_value("tool:1:a:c", ParamValue::Integer(0), "", &[])
        .unwrap();
    load(r#"{"a__b": 1, "a": {"c": 2}, "z": null}"#, &mut param).unwrap();
    assert_eq!(param.value("tool:1:a:b").unwrap(), &ParamValue::Integer(1));
    assert_eq!(param.value("tool:1:a:c").unwrap(), &ParamValue::Integer(2));

    // `a:b` sorts before `a__b` and is applied first; the last value wins.
    load(r#"{"a__b": 4, "a:b": 3}"#, &mut param).unwrap();
    assert_eq!(param.value("tool:1:a:b").unwrap(), &ParamValue::Integer(4));

    let unknown = |key: &str| {
        format!(
            "Parameter {key} passed to 'tool' is invalid. To prevent usage of wrong defaults, please update/fix the parameters!"
        )
    };
    assert_eq!(
        parse_message(load(r#"{"zz": 1, "b": 1}"#, &mut param)),
        unknown("tool:1:b")
    );
    assert_eq!(
        parse_message(load(r#"{"version": "2"}"#, &mut param)),
        unknown("tool:1:version")
    );
    assert_eq!(
        parse_message(load(r#"{"a___b": 1}"#, &mut param)),
        unknown("tool:1:a:_b")
    );
    // An object with a non-string `class` is a section.
    assert_eq!(
        parse_message(load(r#"{"a": {"class": 1}}"#, &mut param)),
        unknown("tool:1:a:class")
    );
    // A failed load changes nothing.
    assert!(load(r#"{"a": {"b": 9, "x": 1}}"#, &mut param).is_err());
    assert_eq!(param.value("tool:1:a:b").unwrap(), &ParamValue::Integer(4));

    // The first section names the tool, and the instance is always 1.
    let mut second = Param::new();
    second
        .set_value("tool:2:a", ParamValue::Integer(0), "", &[])
        .unwrap();
    assert_eq!(
        parse_message(load(r#"{"a": 1}"#, &mut second)),
        unknown("tool:1:a")
    );

    for document in ["[1]", "\"x\"", "3", "true", "null"] {
        assert_eq!(
            parse_message(load(document, &mut param)),
            "Ignoring JSON file 'doc.json' because of unexpected data type. Expecting a dictionary as type.",
            "{document}"
        );
    }
}

/// The value checks of `ParamJSONFile.cpp:93-123` and their `nlohmann::json`
/// messages (json.hpp 4644-4674 and 21353).
#[test]
fn type_errors_carry_the_accessors_messages() {
    let mut param = Param::new();
    param
        .set_value("t:1:flag", "false".into(), "", &[])
        .unwrap();
    param
        .set_valid_strings("t:1:flag", &tags(&["true", "false"]))
        .unwrap();
    param.set_value("t:1:choice", "a".into(), "", &[]).unwrap();
    param
        .set_valid_strings("t:1:choice", &tags(&["a", "b", "c"]))
        .unwrap();
    param
        .set_value("t:1:in", "".into(), "", &tags(&["input file"]))
        .unwrap();
    param
        .set_value(
            "t:1:exe",
            "".into(),
            "",
            &tags(&["input file", "is_executable"]),
        )
        .unwrap();
    param.set_value("t:1:s", "".into(), "", &[]).unwrap();
    param
        .set_value("t:1:l", ParamValue::StringList(vec![]), "", &[])
        .unwrap();
    for (document, message) in [
        (
            r#"{"flag": "true"}"#,
            "[json.exception.type_error.302] type must be boolean, but is string",
        ),
        (
            r#"{"flag": 1}"#,
            "[json.exception.type_error.302] type must be boolean, but is number",
        ),
        (
            r#"{"choice": true}"#,
            "[json.exception.type_error.302] type must be string, but is boolean",
        ),
        (
            r#"{"in": "x"}"#,
            "[json.exception.type_error.305] cannot use operator[] with a string argument with string",
        ),
        (
            r#"{"in": 1}"#,
            "[json.exception.type_error.305] cannot use operator[] with a string argument with number",
        ),
        (
            r#"{"in": {"class": "File"}}"#,
            "[json.exception.type_error.302] type must be string, but is null",
        ),
        (
            r#"{"in": {"class": "File", "path": 1}}"#,
            "[json.exception.type_error.302] type must be string, but is number",
        ),
        (
            r#"{"exe": 1}"#,
            "[json.exception.type_error.302] type must be string, but is number",
        ),
        (
            r#"{"exe": {"class": "File"}}"#,
            "[json.exception.type_error.302] type must be string, but is null",
        ),
        (
            r#"{"s": {"class": "File", "path": "x"}}"#,
            "[json.exception.type_error.302] type must be string, but is object",
        ),
        (
            r#"{"l": "x"}"#,
            "[json.exception.type_error.302] type must be array, but is string",
        ),
        (
            r#"{"l": ["x", null]}"#,
            "[json.exception.type_error.302] type must be string, but is null",
        ),
    ] {
        assert_eq!(
            parse_message(load(document, &mut param)),
            message,
            "{document}"
        );
    }
    // A three-valued choice is a string; an executable is a string or a file.
    load(
        r#"{"choice": "c", "exe": "run.sh", "in": {"class": "File", "path": "p"}}"#,
        &mut param,
    )
    .unwrap();
    assert_eq!(
        param.value("t:1:choice").unwrap(),
        &ParamValue::String("c".into())
    );
    assert_eq!(
        param.value("t:1:exe").unwrap(),
        &ParamValue::String("run.sh".into())
    );
    assert_eq!(
        param.value("t:1:in").unwrap(),
        &ParamValue::String("p".into())
    );
}

/// An entry without a value is left empty with the source's warning
/// (`ParamJSONFile.cpp:152-156`).
#[test]
fn an_entry_without_a_value_is_reported_and_kept() {
    let mut param = Param::new();
    param
        .set_value("t:1:e", ParamValue::Empty, "", &[])
        .unwrap();
    let warnings = load(r#"{"e": [1, "anything"]}"#, &mut param).unwrap();
    assert_eq!(
        warnings,
        vec!["Ignoring entry 't:1:e' because of unknown type 'EMPTY_VALUE'.".to_owned()]
    );
    assert_eq!(param.value("t:1:e").unwrap(), &ParamValue::Empty);
}

/// What `json::parse` accepts and refuses: a byte order mark, and the
/// messages of an empty document and a float beyond `f64` (Release oracle
/// `../oracle/json-ini`, cases `w_pphr_bom`, `w_pphr_empty` and `w_pphr_inf`).
#[test]
fn parser_edges() {
    let mut param = Param::new();
    param
        .set_value("t:1:i", ParamValue::Integer(0), "", &[])
        .unwrap();
    param
        .set_value("t:1:d", ParamValue::Float(0.0), "", &[])
        .unwrap();
    load("\u{feff}{\"i\": 3}", &mut param).unwrap();
    assert_eq!(param.value("t:1:i").unwrap(), &ParamValue::Integer(3));
    assert_eq!(
        parse_message(load("", &mut param)),
        "[json.exception.parse_error.101] parse error at line 1, column 1: syntax error while parsing value - unexpected end of input; expected '[', '{', or a literal"
    );
    assert_eq!(
        parse_message(load("{\"d\": 1e400}", &mut param)),
        "[json.exception.out_of_range.406] number overflow parsing '1e400'"
    );
    assert_eq!(
        parse_message(load("{\n \"i\": 1,\n \"d\": -2.5E+999 }", &mut param)),
        "[json.exception.out_of_range.406] number overflow parsing '-2.5E+999'"
    );
    // Other syntax errors are refused with serde_json's wording.
    for document in ["not a json", "{\"i\": 2,}", "  \n", "{\"i\": 2} x"] {
        assert!(
            matches!(load(document, &mut param), Err(Error::Parse { .. })),
            "{document:?}"
        );
    }
    // A document beyond the bound is refused before it is parsed.
    let oversized = vec![b' '; MAX_JSON_BYTES + 1];
    assert!(matches!(
        param_json::load_from_slice(&oversized, "big.json", &mut param),
        Err(Error::InvalidValue(_))
    ));
    // A tree without a section names no tool.
    assert!(matches!(
        load("{}", &mut Param::new()),
        Err(Error::InvalidValue(_))
    ));
}

/// The store side as the Release build has it, without TDL: `store` and
/// `writeToStream` throw the `std::runtime_error` of `ParamJSONFile.cpp:326`
/// (Release oracle `write_json` of `../oracle/toppbase-completion`, through
/// the command line). This port leaves the target untouched (CPP-350).
#[test]
fn the_store_side_is_the_release_builds_without_tdl() {
    use openms::data_structures::ToolInfo;
    use openms::format::param_json::{ParamJsonFile, TDL_UNAVAILABLE};
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let target = dir.path().join("out.json");
    let param = upstream_param();
    for file in [
        ParamJsonFile::default(),
        ParamJsonFile {
            flat_hierarchy: true,
        },
    ] {
        match file.store(&target.to_string_lossy(), &param, &ToolInfo::default()) {
            Err(Error::Unsupported(message)) => assert_eq!(message, TDL_UNAVAILABLE),
            other => panic!("{other:?}"),
        }
        let mut sink = Vec::new();
        assert!(matches!(
            file.write_to_stream(&mut sink, &param),
            Err(Error::Unsupported(_))
        ));
        assert!(sink.is_empty());
    }
    assert!(!target.exists());
    assert_eq!(
        TDL_UNAVAILABLE,
        "TDL support is not available. Rebuild with -DENABLE_TDL=ON to enable this feature."
    );
}

/// Every number token as the Release build stores it, bit for bit: the B0
/// integration lane's oracle (`../oracle/B0-attacks`, the Linux x86_64 C++
/// Release build at core `bc9cc12` on ibminode06, gcc 14.4.0, two identical
/// runs; `tests/data/b0_release/json_release.tsv`, provenance
/// `tests/data/b0_release_provenance.json`). A driver loads each document of
/// `tests/data/b0_release/json` with `ParamJSONFile::load` onto a tree with a
/// floating-point, a 64-bit integer, a floating-point list and an integer
/// list parameter and prints the stored values, doubles as their IEEE-754
/// bits. The tokens probe what `serde_json` with `float_roundtrip` and the
/// `nlohmann::json` 3.11.2 lexer (`strtod`, `strtoull`, `strtoll`) must agree
/// on: correctly rounded halfway cases and long mantissas, subnormals and
/// underflow, the largest double, integers beyond 2^53 and beyond 64 bits,
/// `-0` against `-0.0`, and the integer accessors' casts of each.
#[test]
fn number_tokens_store_the_release_builds_bits() {
    let root = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/b0_release");
    let records = fs::read_to_string(format!("{root}/json_release.tsv")).unwrap();
    let float = |param: &Param, key: &str| match param.value(key).unwrap() {
        ParamValue::Float(value) => format!("{:016x}", value.to_bits()),
        other => panic!("{key}: {other:?}"),
    };
    let mut compared = 0;
    for line in records.lines().filter(|line| *line != "DONE") {
        let fields: Vec<&str> = line.split('\t').collect();
        let case = fields[0];
        let mut param = Param::new();
        for (key, value) in [
            ("T:1:d", ParamValue::Float(0.0)),
            ("T:1:i", ParamValue::Integer(0)),
            ("T:1:dl", ParamValue::FloatList(Vec::new())),
            ("T:1:il", ParamValue::IntegerList(Vec::new())),
        ] {
            param.set_value(key, value, "", &[]).unwrap();
        }
        let result = param_json::load(format!("{root}/json/{case}.json"), &mut param);
        match fields[1] {
            "ok" => {
                result.unwrap_or_else(|error| panic!("{case}: {error}"));
                assert_eq!(float(&param, "T:1:d"), fields[2], "{case}: d");
                let ParamValue::Integer(integer) = param.value("T:1:i").unwrap() else {
                    panic!("{case}: i");
                };
                assert_eq!(integer.to_string(), fields[3], "{case}: i");
                let ParamValue::FloatList(floats) = param.value("T:1:dl").unwrap() else {
                    panic!("{case}: dl");
                };
                let floats: Vec<String> = floats
                    .iter()
                    .map(|value| format!("{:016x}", value.to_bits()))
                    .collect();
                assert_eq!(floats.join(","), fields[4], "{case}: dl");
                let ParamValue::IntegerList(integers) = param.value("T:1:il").unwrap() else {
                    panic!("{case}: il");
                };
                let integers: Vec<String> = integers.iter().map(i32::to_string).collect();
                assert_eq!(integers.join(","), fields[5], "{case}: il");
            }
            "X" => {
                // `ParseError` with an empty expression: `<message> in: `.
                assert_eq!(fields[2], "Parse Error", "{case}");
                assert_eq!(
                    format!("{} in: ", parse_message(result)),
                    fields[3],
                    "{case}"
                );
            }
            other => panic!("{case}: unexpected Release record {other}"),
        }
        compared += 1;
    }
    assert_eq!(compared, 51);
}
