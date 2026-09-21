// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! `FORMAT/ParamCTDFile.h`: the CTD writer behind `-write_ctd`.
//!
//! * `ParamCTDFile_test.cpp` (core `bc9cc12`) is transcribed section by
//!   section. Its `writeCTDToStream` section compares with the retained C++
//!   output `ParamCTDFile_test_writeCTDToStream.ctd` (tier 1); the other
//!   sections read the written file back through the parameter-XML reader
//!   (tier 3, the class test's own literals).
//! * `../oracle/toppbase-completion/ctd_driver.cpp` writes five documents
//!   through the Release build's `ParamCTDFile::writeCTDToStream`, built to
//!   reach every branch and the source's escaping defects; each is compared
//!   byte for byte (tier 1, `tests/data/param_ctd/oracle/`).
#![cfg(feature = "paramxml")]

use openms::cli::ParamCtdFile;
use openms::data_structures::ToolInfo;
use openms::format::paramxml;
use openms::param::{Param, ParamValue};
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/param_ctd")
        .join(name)
}

fn strings(values: &[&str]) -> Vec<String> {
    values.iter().map(|v| (*v).to_owned()).collect()
}

fn info(values: [&str; 5], citations: &[&str]) -> ToolInfo {
    ToolInfo {
        version: values[0].into(),
        name: values[1].into(),
        docurl: values[2].into(),
        category: values[3].into(),
        description: values[4].into(),
        citations: strings(citations),
    }
}

fn set(p: &mut Param, key: &str, value: ParamValue, description: &str, tags: &[&str]) {
    p.set_value(key, value, description, &strings(tags))
        .unwrap();
}

fn text(value: &str) -> ParamValue {
    ParamValue::String(value.into())
}

/// A fresh file name in the system temporary directory.
fn temporary(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms-param-ctd-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir.join("out.ctd")
}

/// Read a written CTD back as `ParamXMLFile::load` does in the class test.
///
/// This port's parameter-XML reader accepts `<PARAMETERS>` as the root only,
/// where the source's handler also reads it inside the CTD's `<tool>` wrapper;
/// the test hands the reader the `<PARAMETERS>` element of the written file.
fn load_ctd(path: &Path) -> Param {
    let written = std::fs::read_to_string(path).unwrap();
    let start = written.find("<PARAMETERS").unwrap();
    let end = written.find("</PARAMETERS>").unwrap() + "</PARAMETERS>".len();
    paramxml::read(&written.as_bytes()[start..end]).unwrap()
}

/// The shared fixture of the class test (`ParamCTDFile_test.cpp:46-54`).
fn fixture() -> Param {
    let mut p = Param::new();
    set(
        &mut p,
        "test:float",
        ParamValue::Float(f64::from(17.4f32)),
        "floatdesc",
        &[],
    );
    set(
        &mut p,
        "test:string",
        text("test,test,test"),
        "stringdesc",
        &[],
    );
    set(&mut p, "test:int", ParamValue::Integer(17), "intdesc", &[]);
    set(
        &mut p,
        "test2:float",
        ParamValue::Float(f64::from(17.5f32)),
        "",
        &[],
    );
    set(&mut p, "test2:string", text("test2"), "", &[]);
    set(&mut p, "test2:int", ParamValue::Integer(18), "", &[]);
    p.set_section_description("test", "sectiondesc").unwrap();
    p.add_tags("test:float", &strings(&["a", "b", "c"]))
        .unwrap();
    p
}

/// `START_SECTION((void store(...)))`, the round trip and the missing
/// directory (`ParamCTDFile_test.cpp:56-93`).
#[test]
fn upstream_store_round_trips_values_descriptions_and_tags() {
    let file = ParamCtdFile;
    let info_a = info(["a"; 5], &[]);
    // TEST_EXCEPTION(std::ios::failure, store("/does/not/exist/...")): the
    // source's message is `Unable to create file: <name>`.
    let error = file
        .store(
            "/does/not/exist/FileDoesNotExist.ctd",
            &Param::new(),
            &info_a,
        )
        .unwrap_err();
    assert_eq!(
        error.to_string(),
        "Unable to create file: /does/not/exist/FileDoesNotExist.ctd"
    );

    let mut p2 = fixture();
    set(
        &mut p2,
        "test:a:a1",
        ParamValue::Float(47.1),
        "a1desc\"<>\nnewline",
        &[],
    );
    set(&mut p2, "test:b:b1", ParamValue::Float(47.1), "", &[]);
    p2.set_section_description("test:b", "bdesc\"<>\nnewline")
        .unwrap();
    set(&mut p2, "test2:a:a1", ParamValue::Float(47.1), "", &[]);
    set(
        &mut p2,
        "test2:b:b1",
        ParamValue::Float(47.1),
        "",
        &["advanced"],
    );
    p2.set_section_description("test2:a", "adesc").unwrap();
    let path = temporary("store");
    file.store(&path, &p2, &info_a).unwrap();
    let p3 = load_ctd(&path);
    let float = |p: &Param, key: &str| p.value(key).unwrap().to_f64().unwrap() as f32;
    assert!((float(&p2, "test:float") - float(&p3, "test:float")).abs() < 1e-5);
    assert_eq!(
        p2.value("test:string").unwrap(),
        p3.value("test:string").unwrap()
    );
    assert_eq!(p2.value("test:int").unwrap(), p3.value("test:int").unwrap());
    assert!((float(&p2, "test2:float") - float(&p3, "test2:float")).abs() < 1e-5);
    assert_eq!(
        p2.value("test2:string").unwrap(),
        p3.value("test2:string").unwrap()
    );
    assert_eq!(
        p2.value("test2:int").unwrap(),
        p3.value("test2:int").unwrap()
    );
    for key in ["test:float", "test:string", "test:int"] {
        assert_eq!(p2.description(key).unwrap(), p3.description(key).unwrap());
    }
    assert_eq!(p3.section_description("test").unwrap(), "sectiondesc");
    assert_eq!(p3.description("test:a:a1").unwrap(), "a1desc\"<>\nnewline");
    assert_eq!(
        p3.section_description("test:b").unwrap(),
        "bdesc\"<>\nnewline"
    );
    assert_eq!(p3.section_description("test2:a").unwrap(), "adesc");
    assert!(p3.has_tag("test2:b:b1", "advanced").unwrap());
    assert!(!p3.has_tag("test2:a:a1", "advanced").unwrap());
}

/// The `advanced` and restriction parts of the store section
/// (`ParamCTDFile_test.cpp:95-213`).
#[test]
fn upstream_store_keeps_advanced_flags_and_restrictions() {
    let file = ParamCtdFile;
    let info_a = info(["a"; 5], &[]);
    let mut p7 = Param::new();
    set(&mut p7, "true", ParamValue::Integer(5), "", &["advanced"]);
    set(&mut p7, "false", ParamValue::Integer(5), "", &[]);
    let path = temporary("advanced");
    file.store(&path, &p7, &info_a).unwrap();
    let p8 = load_ctd(&path);
    assert!(p8.entry("true").unwrap().tags.contains("advanced"));
    assert!(!p8.entry("false").unwrap().tags.contains("advanced"));

    let mut p5 = Param::new();
    for key in ["int", "int_min", "int_max", "int_min_max"] {
        set(&mut p5, key, ParamValue::Integer(5), "", &[]);
    }
    p5.set_min_int("int_min", 4).unwrap();
    p5.set_max_int("int_max", 6).unwrap();
    p5.set_min_int("int_min_max", 0).unwrap();
    p5.set_max_int("int_min_max", 10).unwrap();
    for key in ["float", "float_min", "float_max", "float_min_max"] {
        set(&mut p5, key, ParamValue::Float(5.1), "", &[]);
    }
    p5.set_min_float("float_min", 4.1).unwrap();
    p5.set_max_float("float_max", 6.1).unwrap();
    p5.set_min_float("float_min_max", 0.1).unwrap();
    p5.set_max_float("float_min_max", 10.1).unwrap();
    set(&mut p5, "string", text("bli"), "", &[]);
    set(&mut p5, "string_2", text("bla"), "", &[]);
    p5.set_valid_strings("string_2", &strings(&["bla", "bluff"]))
        .unwrap();
    set(
        &mut p5,
        "stringlist2",
        ParamValue::StringList(strings(&["a.txt", "b.xml", "c.pdf"])),
        "",
        &[],
    );
    set(
        &mut p5,
        "stringlist",
        ParamValue::StringList(strings(&["aa.C", "bb.h", "c.doxygen"])),
        "",
        &[],
    );
    p5.set_valid_strings("stringlist2", &strings(&["xml", "txt"]))
        .unwrap();
    for key in ["intlist", "intlist2", "intlist3", "intlist4"] {
        set(
            &mut p5,
            key,
            ParamValue::IntegerList(vec![2, 5, 10]),
            "",
            &[],
        );
    }
    p5.set_min_int("intlist2", 1).unwrap();
    p5.set_max_int("intlist3", 11).unwrap();
    p5.set_min_int("intlist4", 0).unwrap();
    p5.set_max_int("intlist4", 15).unwrap();
    for key in ["doublelist", "doublelist2", "doublelist3", "doublelist4"] {
        set(
            &mut p5,
            key,
            ParamValue::FloatList(vec![1.2, 3.33, 4.44]),
            "",
            &[],
        );
    }
    p5.set_min_float("doublelist2", 1.1).unwrap();
    p5.set_max_float("doublelist3", 4.45).unwrap();
    p5.set_min_float("doublelist4", 0.1).unwrap();
    p5.set_max_float("doublelist4", 5.8).unwrap();
    let path = temporary("restrictions");
    file.store(&path, &p5, &info_a).unwrap();
    let p6 = load_ctd(&path);
    let ints = |key: &str| {
        let e = p6.entry(key).unwrap();
        (e.min_int, e.max_int)
    };
    assert_eq!(ints("int"), (-i32::MAX, i32::MAX));
    assert_eq!(ints("int_min"), (4, i32::MAX));
    assert_eq!(ints("int_max"), (-i32::MAX, 6));
    assert_eq!(ints("int_min_max"), (0, 10));
    let floats = |key: &str| {
        let e = p6.entry(key).unwrap();
        (e.min_float, e.max_float)
    };
    let similar = |(a, b): (f64, f64), (c, d): (f64, f64)| {
        let close = |x: f64, y: f64| x == y || ((x - y) / y).abs() < 1e-5;
        assert!(close(a, c) && close(b, d), "({a}, {b}) != ({c}, {d})");
    };
    similar(floats("float"), (-f64::MAX, f64::MAX));
    similar(floats("float_min"), (4.1, f64::MAX));
    similar(floats("float_max"), (-f64::MAX, 6.1));
    similar(floats("float_min_max"), (0.1, 10.1));
    assert!(p6.entry("string").unwrap().valid_strings.is_empty());
    assert_eq!(
        p6.entry("string_2").unwrap().valid_strings,
        strings(&["bla", "bluff"])
    );
    assert!(p6.entry("stringlist").unwrap().valid_strings.is_empty());
    assert_eq!(
        p6.entry("stringlist2").unwrap().valid_strings,
        strings(&["xml", "txt"])
    );
    assert_eq!(ints("intlist"), (-i32::MAX, i32::MAX));
    assert_eq!(ints("intlist2"), (1, i32::MAX));
    assert_eq!(ints("intlist3"), (-i32::MAX, 11));
    assert_eq!(ints("intlist4"), (0, 15));
    similar(floats("doublelist"), (-f64::MAX, f64::MAX));
    similar(floats("doublelist2"), (1.1, f64::MAX));
    similar(floats("doublelist3"), (-f64::MAX, 4.45));
    similar(floats("doublelist4"), (0.1, 5.8));
}

/// The NaN part of the store section (`ParamCTDFile_test.cpp:215-229`): a NaN
/// reads back as NaN, and the file spells it `NaN` on its eighth and ninth
/// lines.
#[test]
fn upstream_store_writes_nan_as_capitalised_nan() {
    let mut p = Param::new();
    set(
        &mut p,
        "float_nan",
        ParamValue::Float(f64::from(f32::NAN)),
        "",
        &[],
    );
    set(&mut p, "double_nan", ParamValue::Float(f64::NAN), "", &[]);
    let path = temporary("nan");
    ParamCtdFile.store(&path, &p, &info(["a"; 5], &[])).unwrap();
    let lines: Vec<String> = std::fs::read_to_string(&path)
        .unwrap()
        .lines()
        .map(str::to_owned)
        .collect();
    assert!(lines[7].contains("value=\"NaN\""), "{}", lines[7]);
    assert!(lines[8].contains("value=\"NaN\""), "{}", lines[8]);
    let back = load_ctd(&path);
    assert!(back.value("float_nan").unwrap().to_f64().unwrap().is_nan());
    assert!(back.value("double_nan").unwrap().to_f64().unwrap().is_nan());
}

/// `START_SECTION((void writeCTDToStream(...)))` (`ParamCTDFile_test.cpp:233-272`):
/// `TEST_FILE_EQUAL` against the retained C++ output. `TEST_FILE_EQUAL`
/// compares line by line with surrounding white space trimmed; the retained
/// file was committed with CRLF line ends, so the comparison here is exact
/// once those are removed, which is stricter.
#[test]
fn upstream_write_ctd_to_stream_matches_the_retained_file() {
    let mut p = Param::new();
    set(
        &mut p,
        "stringlist",
        ParamValue::StringList(strings(&["a", "bb", "ccc"])),
        "StringList Description",
        &[],
    );
    set(
        &mut p,
        "intlist",
        ParamValue::IntegerList(vec![1, 22, 333]),
        "",
        &[],
    );
    set(&mut p, "item", text("bla"), "", &[]);
    set(
        &mut p,
        "stringlist2",
        ParamValue::StringList(Vec::new()),
        "",
        &[],
    );
    set(
        &mut p,
        "intlist2",
        ParamValue::IntegerList(Vec::new()),
        "",
        &[],
    );
    set(&mut p, "item1", ParamValue::Integer(7), "", &[]);
    set(
        &mut p,
        "intlist3",
        ParamValue::IntegerList(vec![1]),
        "",
        &[],
    );
    set(
        &mut p,
        "stringlist3",
        ParamValue::StringList(strings(&["1"])),
        "",
        &[],
    );
    set(&mut p, "item3", ParamValue::Float(7.6), "", &[]);
    set(
        &mut p,
        "doublelist",
        ParamValue::FloatList(vec![1.22, 2.33, 4.55]),
        "",
        &[],
    );
    set(
        &mut p,
        "doublelist3",
        ParamValue::FloatList(vec![1.4]),
        "",
        &[],
    );
    set(
        &mut p,
        "file_parameter",
        text(""),
        "This is a file parameter.",
        &[],
    );
    p.add_tag("file_parameter", "input file").unwrap();
    p.set_valid_strings("file_parameter", &strings(&["*.mzML", "*.mzXML"]))
        .unwrap();
    set(
        &mut p,
        "outdir_parameter",
        text(""),
        "This is a outdir parameter.",
        &[],
    );
    p.add_tag("outdir_parameter", "output dir").unwrap();
    set(
        &mut p,
        "advanced_parameter",
        text(""),
        "This is an advanced parameter.",
        &["advanced"],
    );
    set(
        &mut p,
        "flag",
        text("false"),
        "This is a flag i.e. in a command line input it does not need a value.",
        &[],
    );
    p.set_valid_strings("flag", &strings(&["true", "false"]))
        .unwrap();
    set(
        &mut p,
        "noflagJustTrueFalse",
        text("true"),
        "This is not a flag but has a boolean meaning.",
        &[],
    );
    p.set_valid_strings("noflagJustTrueFalse", &strings(&["true", "false"]))
        .unwrap();
    let tool = info(
        [
            "2.6.0-pre-STL-ParamCTD-2021-06-02",
            "AccurateMassSearch",
            "http://www.openms.de/doxygen/nightly/html/TOPP_AccurateMassSearch.html",
            "Utilities",
            "Match MS signals to molecules from a database by mass.",
        ],
        &["10.1038/s41592-024-02197-7"],
    );
    let mut written = Vec::new();
    ParamCtdFile
        .write_ctd_to_stream(&mut written, &p, &tool)
        .unwrap();
    let expected = std::fs::read_to_string(data("ParamCTDFile_test_writeCTDToStream.ctd"))
        .unwrap()
        .replace("\r\n", "\n");
    assert_eq!(String::from_utf8(written).unwrap(), expected);
}

/// `START_SECTION([EXTRA] storing of lists)` (`ParamCTDFile_test.cpp:274-350`).
#[test]
fn upstream_storing_of_lists_round_trips() {
    let mut p = Param::new();
    set(
        &mut p,
        "stringlist",
        ParamValue::StringList(strings(&["a", "bb", "ccc"])),
        "",
        &[],
    );
    set(
        &mut p,
        "intlist",
        ParamValue::IntegerList(vec![1, 22, 333]),
        "",
        &[],
    );
    set(&mut p, "item", text("bla"), "", &[]);
    set(
        &mut p,
        "stringlist2",
        ParamValue::StringList(Vec::new()),
        "",
        &[],
    );
    set(
        &mut p,
        "intlist2",
        ParamValue::IntegerList(Vec::new()),
        "",
        &[],
    );
    set(&mut p, "item1", ParamValue::Integer(7), "", &[]);
    set(
        &mut p,
        "intlist3",
        ParamValue::IntegerList(vec![1]),
        "",
        &[],
    );
    set(
        &mut p,
        "stringlist3",
        ParamValue::StringList(strings(&["1"])),
        "",
        &[],
    );
    set(&mut p, "item3", ParamValue::Float(7.6), "", &[]);
    set(
        &mut p,
        "doublelist",
        ParamValue::FloatList(vec![1.22, 2.33, 4.55]),
        "",
        &[],
    );
    set(
        &mut p,
        "doublelist2",
        ParamValue::FloatList(Vec::new()),
        "",
        &[],
    );
    set(
        &mut p,
        "doublelist3",
        ParamValue::FloatList(vec![1.4]),
        "",
        &[],
    );
    let path = temporary("lists");
    ParamCtdFile
        .store(&path, &p, &info(["a", "b", "c", "d", "e"], &["f"]))
        .unwrap();
    let p2 = load_ctd(&path);
    assert_eq!(p2.size(), 12);
    assert_eq!(
        p2.value("stringlist").unwrap(),
        &ParamValue::StringList(strings(&["a", "bb", "ccc"]))
    );
    assert_eq!(
        p2.value("stringlist2").unwrap(),
        &ParamValue::StringList(Vec::new())
    );
    assert_eq!(
        p2.value("stringlist3").unwrap(),
        &ParamValue::StringList(strings(&["1"]))
    );
    assert_eq!(
        p2.value("intlist").unwrap(),
        &ParamValue::IntegerList(vec![1, 22, 333])
    );
    assert_eq!(
        p2.value("intlist2").unwrap(),
        &ParamValue::IntegerList(Vec::new())
    );
    assert_eq!(
        p2.value("intlist3").unwrap(),
        &ParamValue::IntegerList(vec![1])
    );
    assert_eq!(
        p2.value("doublelist").unwrap(),
        &ParamValue::FloatList(vec![1.22, 2.33, 4.55])
    );
    assert_eq!(
        p2.value("doublelist2").unwrap(),
        &ParamValue::FloatList(Vec::new())
    );
    assert_eq!(
        p2.value("doublelist3").unwrap(),
        &ParamValue::FloatList(vec![1.4])
    );
}

/// `START_SECTION([EXTRA] Escaping of characters)` (`ParamCTDFile_test.cpp:353-391`).
///
/// The class test checks the reloaded tree for one description only and its
/// other assertions against the tree it wrote; both are transcribed, and the
/// reloaded values are checked too.
#[test]
fn upstream_escaping_of_characters_round_trips() {
    let entries = [
        ("string", "bla", "string"),
        (
            "string_with_ampersand",
            "bla2&blubb",
            "string with ampersand",
        ),
        (
            "string_with_ampersand_in_descr",
            "blaxx",
            "std::string with & in description",
        ),
        (
            "string_with_single_quote",
            "bla'xxx",
            "std::string with single quotes",
        ),
        (
            "string_with_single_quote_in_descr",
            "blaxxx",
            "std::string with ' quote in description",
        ),
        (
            "string_with_double_quote",
            "bla\"xxx",
            "std::string with double quote",
        ),
        (
            "string_with_double_quote_in_descr",
            "bla\"xxx",
            "std::string with \" description",
        ),
        (
            "string_with_greater_sign",
            "bla>xxx",
            "std::string with greater sign",
        ),
        (
            "string_with_greater_sign_in_descr",
            "bla greater xxx",
            "std::string with >",
        ),
        (
            "string_with_less_sign",
            "bla<xxx",
            "std::string with less sign",
        ),
        (
            "string_with_less_sign_in_descr",
            "bla less sign_xxx",
            "std::string with less sign <",
        ),
    ];
    let mut p = Param::new();
    for (key, value, description) in entries {
        set(&mut p, key, text(value), description, &[]);
    }
    let path = temporary("escaping");
    ParamCtdFile
        .store(&path, &p, &info(["a"; 5], &["a"]))
        .unwrap();
    let p2 = load_ctd(&path);
    assert_eq!(p2.description("string").unwrap(), "string");
    for (key, value, description) in entries {
        assert_eq!(p.value(key).unwrap(), &text(value));
        assert_eq!(p.description(key).unwrap(), description);
        assert_eq!(p2.value(key).unwrap(), &text(value));
        assert_eq!(p2.description(key).unwrap(), description);
    }
}

/// The metadata the oracle driver writes for its `escaping` and `empty`
/// documents: none of it is escaped by the source.
fn driver_info() -> ToolInfo {
    info(
        [
            "1.2.3",
            "Probe",
            "http://example.org/doc?a=1&b=2",
            "Cat & \"Dog\" <x>",
            "Desc with <b>markup</b> & ]]> end",
        ],
        &["10.1/abc", "https://doi.org/x?y=1&z=2"],
    )
}

fn oracle(name: &str) -> String {
    std::fs::read_to_string(data(&format!("oracle/{name}.ctd"))).unwrap()
}

fn written(p: &Param, tool: &ToolInfo) -> String {
    ParamCtdFile.to_ctd_string(p, tool).unwrap()
}

/// Oracle `escaping`: the source's escaping, including its defects — of two
/// adjacent special characters only the first is escaped, a second
/// consecutive line break stays raw, a tab becomes `&amp;#x9;`, and the
/// `<tool>` attributes, the description and the citations are not escaped.
#[test]
fn escaping_matches_the_release_build_byte_for_byte() {
    let mut p = Param::new();
    set(&mut p, "amp", text("a&b"), "one & two", &[]);
    set(&mut p, "amp2", text("a&&b"), "two && adjacent", &[]);
    set(&mut p, "amp3", text("&&&"), "three &&&", &[]);
    set(&mut p, "mixed", text("<>\"'&"), "all <>\"'& five", &[]);
    set(&mut p, "gtgt", text(">>"), ">>", &[]);
    set(&mut p, "quotes", text("\"\""), "''", &[]);
    set(&mut p, "tab", text("a\tb"), "tab in value", &[]);
    set(&mut p, "tab2", text("a\t\tb"), "two tabs", &[]);
    set(&mut p, "newline", text("x"), "line1\nline2", &[]);
    set(&mut p, "newline2", text("x"), "line1\n\nline3", &[]);
    set(&mut p, "newline3", text("x"), "\n\n\n", &[]);
    set(
        &mut p,
        "utf8",
        text("\u{e4}\u{f6}\u{fc}"),
        "\u{20ac} & \u{df}",
        &[],
    );
    set(
        &mut p,
        "list:tabs",
        ParamValue::StringList(strings(&["a\tb", "c\t\td", "e&f", "g&&h"])),
        "list desc",
        &[],
    );
    p.set_section_description("list", "section <desc>\n\nwith & newlines")
        .unwrap();
    set(
        &mut p,
        "names:a&b",
        ParamValue::Integer(1),
        "name with amp",
        &[],
    );
    assert_eq!(written(&p, &driver_info()), oracle("escaping"));
}

/// Oracle `numbers`: numbers use `ParamValue::toString`, list items the
/// stream's `%.15g`, restrictions `std::to_string`'s `%f`.
#[test]
fn numbers_and_restrictions_match_the_release_build_byte_for_byte() {
    let mut p = Param::new();
    let int = |p: &mut Param, key: &str, v: i64| set(p, key, ParamValue::Integer(v), "", &[]);
    let float = |p: &mut Param, key: &str, v: f64| set(p, key, ParamValue::Float(v), "", &[]);
    int(&mut p, "i", 5);
    int(&mut p, "i_min", 5);
    p.set_min_int("i_min", -3).unwrap();
    int(&mut p, "i_max", 5);
    p.set_max_int("i_max", 7).unwrap();
    int(&mut p, "i_both", 5);
    p.set_min_int("i_both", 0).unwrap();
    p.set_max_int("i_both", 10).unwrap();
    int(&mut p, "i_neg", -2_147_483_647);
    float(&mut p, "d", 3.0);
    float(&mut p, "d_small", 1e-5);
    float(&mut p, "d_big", 123_456.789);
    float(&mut p, "d_huge", 1e300);
    float(&mut p, "d_third", 1.0 / 3.0);
    float(&mut p, "d_nan", f64::NAN);
    float(&mut p, "d_inf", f64::INFINITY);
    float(&mut p, "d_ninf", f64::NEG_INFINITY);
    float(&mut p, "d_negzero", -0.0);
    float(&mut p, "d_min", 5.1);
    p.set_min_float("d_min", 4.1).unwrap();
    float(&mut p, "d_max", 5.1);
    p.set_max_float("d_max", 1e20).unwrap();
    float(&mut p, "d_both", 5.1);
    p.set_min_float("d_both", -0.5).unwrap();
    p.set_max_float("d_both", 0.123_456_789).unwrap();
    float(&mut p, "d_tiny_min", 5.1);
    p.set_min_float("d_tiny_min", 1e-10).unwrap();
    set(
        &mut p,
        "il",
        ParamValue::IntegerList(vec![1, -22, 333]),
        "",
        &[],
    );
    set(&mut p, "il_r", ParamValue::IntegerList(vec![1, 2]), "", &[]);
    p.set_min_int("il_r", 1).unwrap();
    p.set_max_int("il_r", 2).unwrap();
    set(
        &mut p,
        "dl",
        ParamValue::FloatList(vec![1.22, 2.33, 4.55]),
        "",
        &[],
    );
    set(
        &mut p,
        "dl_g",
        ParamValue::FloatList(vec![
            0.1 + 0.2,
            1e-5,
            1e20,
            123_456_789_012_345_678.0,
            1.0 / 3.0,
            100.0,
            1e15,
            1e16,
            -0.0,
            f64::NAN,
            f64::INFINITY,
        ]),
        "",
        &[],
    );
    set(&mut p, "dl_r", ParamValue::FloatList(vec![1.5]), "", &[]);
    p.set_min_float("dl_r", 0.25).unwrap();
    set(
        &mut p,
        "empty_dl",
        ParamValue::FloatList(Vec::new()),
        "",
        &[],
    );
    set(
        &mut p,
        "empty_il",
        ParamValue::IntegerList(Vec::new()),
        "",
        &[],
    );
    set(
        &mut p,
        "empty_sl",
        ParamValue::StringList(Vec::new()),
        "",
        &[],
    );
    assert_eq!(written(&p, &ToolInfo::default()), oracle("numbers"));
}

/// Oracle `types`: file types from tags, flags only for `false` restricted to
/// `true,false` in that order, the remaining tags sorted, and nodes after a
/// section's own entries.
#[test]
fn types_tags_flags_and_nesting_match_the_release_build_byte_for_byte() {
    let mut p = Param::new();
    set(&mut p, "in", text(""), "in", &["input file", "required"]);
    p.set_valid_strings("in", &strings(&["*.mzML", "*.mzXML"]))
        .unwrap();
    set(
        &mut p,
        "exe",
        text("java"),
        "exe",
        &["input file", "is_executable"],
    );
    set(
        &mut p,
        "out",
        text("x.mzML"),
        "out",
        &["output file", "advanced"],
    );
    p.set_valid_strings("out", &strings(&["*.mzML"])).unwrap();
    set(&mut p, "prefix", text("p"), "prefix", &["output prefix"]);
    p.set_valid_strings("prefix", &strings(&["*.mgf"])).unwrap();
    set(&mut p, "dir", text("d"), "dir", &["output dir"]);
    set(
        &mut p,
        "ins",
        ParamValue::StringList(strings(&["a.mzML", "b.mzML"])),
        "ins",
        &["input file"],
    );
    p.set_valid_strings("ins", &strings(&["*.mzML"])).unwrap();
    set(
        &mut p,
        "outs",
        ParamValue::StringList(Vec::new()),
        "outs",
        &["output file", "required"],
    );
    set(
        &mut p,
        "plain_list",
        ParamValue::StringList(strings(&["x"])),
        "plain",
        &["zeta", "alpha", "advanced"],
    );
    p.set_valid_strings("plain_list", &strings(&["x", "y"]))
        .unwrap();
    set(&mut p, "flag", text("false"), "flag", &[]);
    p.set_valid_strings("flag", &strings(&["true", "false"]))
        .unwrap();
    set(&mut p, "flag_true", text("true"), "not a flag", &[]);
    p.set_valid_strings("flag_true", &strings(&["true", "false"]))
        .unwrap();
    set(&mut p, "flag_order", text("false"), "reversed order", &[]);
    p.set_valid_strings("flag_order", &strings(&["false", "true"]))
        .unwrap();
    set(
        &mut p,
        "tagged",
        text("v"),
        "tagged",
        &["b", "a", "required", "advanced"],
    );
    set(&mut p, "a:b:c:deep", ParamValue::Integer(1), "deep", &[]);
    set(
        &mut p,
        "a:b:sibling",
        ParamValue::Integer(2),
        "sibling",
        &[],
    );
    set(&mut p, "a:after", ParamValue::Integer(3), "after", &[]);
    p.set_section_description("a", "A").unwrap();
    p.set_section_description("a:b", "B").unwrap();
    set(&mut p, "z_last", ParamValue::Integer(4), "last", &[]);
    assert_eq!(written(&p, &ToolInfo::default()), oracle("types"));
}

/// Oracles `nested_only` and `empty`: sections still open at the end are
/// closed from the end trace, and an empty tree writes an empty
/// `<PARAMETERS>` element.
#[test]
fn open_sections_and_an_empty_tree_match_the_release_build_byte_for_byte() {
    let mut p = Param::new();
    set(&mut p, "x:y:z", ParamValue::Integer(1), "deep only", &[]);
    assert_eq!(written(&p, &ToolInfo::default()), oracle("nested_only"));
    assert_eq!(written(&Param::new(), &driver_info()), oracle("empty"));
}
