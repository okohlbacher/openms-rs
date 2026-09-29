// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! `-write_ini` byte for byte against every retained C++ `-write_ini` file of
//! a ported tool, where the tools' own tests compare numbers by value (`4`
//! against `4.0`): the text of the parameter XML writer's numbers,
//! [`paramxml::NumberText::Source`], which `-write_ini` selects.
//!
//! The retained files, each written by the tool of the same name with the
//! same arguments (the provenance is the named test's):
//!
//! * `tests/data/topp_cli_lifecycle/write_ini_<tool>.ini`, oracle
//!   `write_ini_<tool>` (`tests/topp_cli_lifecycle.rs`,
//!   `write_ini_matches_the_cpp_file_for_every_ported_tool`), `-test -write_ini`;
//! * `tests/data/topp_peak_picker_hi_res/oracle_write_ini.ini`, P3 oracle
//!   `write_ini` (`tests/topp_peak_picker_hi_res.rs`), `-write_ini`;
//! * `tests/data/topp_file_info/expected/write_ini.ini`
//!   (`tests/topp_file_info.rs`), `-test -write_ini`;
//! * `tests/data/topp_fuzzy_diff/oracle/write_ini.ini`
//!   (`tests/topp_fuzzy_diff.rs`), `-write_ini`;
//! * `tests/data/feature_finder_picked/FeatureFinderCentroided_defaults.ini`,
//!   oracles `FFC_write_ini` and `FFC_write_ini_test`
//!   (`tests/topp_feature_finder_centroided.rs`), with and without `-test`.
//!
//! Between them they hold every kind of number the writer spells: fixed
//! floating-point items with and without a fraction (`3.0`, `0.333`), one
//! whose 15 fractional digits show its binary value (`99.632000000000005`),
//! one in scientific notation (`1.0e20`), open and closed floating-point
//! restrictions (`0.0:`, `0.0:1.0`), integer items and restrictions, and
//! integer and floating-point list items. The `version` item's value is the
//! only text not compared: it names the build that wrote the file.
#![cfg(feature = "paramxml")]

use openms::cli::{ExitCode, Tool, run_with};
use openms::format::paramxml;
use openms::system::file::TempDir;
use std::fs;
use std::path::Path;

fn repository(relative: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(relative)
        .to_string_lossy()
        .into_owned()
}

/// The INI `<tool> [-test] -write_ini <file>` writes.
fn written<T: Tool>(test: bool) -> String {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join(format!("{}.tmp.ini", T::NAME));
    let mut arguments = vec![T::NAME.to_owned()];
    if test {
        arguments.push("-test".to_owned());
    }
    arguments.push("-write_ini".to_owned());
    arguments.push(path.to_string_lossy().into_owned());
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<T>(&arguments, &mut out, &mut err);
    assert_eq!(
        code,
        ExitCode::ExecutionOk,
        "{}: {}",
        T::NAME,
        String::from_utf8_lossy(&err)
    );
    fs::read_to_string(&path).unwrap()
}

/// `text` with the value of the tool's `version` item replaced by `*`.
fn without_version(text: &str) -> String {
    const ITEM: &str = "<ITEM name=\"version\" value=\"";
    text.lines()
        .map(|line| match line.find(ITEM) {
            Some(start) => {
                let value = start + ITEM.len();
                let end = line[value..].find('"').map_or(line.len(), |at| value + at);
                format!("{}*{}\n", &line[..value], &line[end..])
            }
            None => format!("{line}\n"),
        })
        .collect()
}

/// The INI `T` writes equals the retained file `expected`, byte for byte
/// apart from the `version` item's value, and the two trees read back equal.
fn assert_written<T: Tool>(expected: &str, test: bool) {
    let port = written::<T>(test);
    let release = fs::read_to_string(repository(expected)).unwrap();
    assert!(
        release.contains("<ITEM name=\"version\" value=\""),
        "{expected}"
    );
    assert_eq!(
        without_version(&port),
        without_version(&release),
        "{} against {expected}",
        T::NAME
    );
    let tree = |text: &str| {
        let mut param = paramxml::read(text.as_bytes()).unwrap();
        let version = format!("{}:version", T::NAME);
        param.remove(&version).unwrap();
        param
    };
    assert_eq!(tree(&port), tree(&release), "{}", T::NAME);
}

#[test]
fn write_ini_is_the_retained_text_for_the_framework_tools() {
    use openms::cli::tools::{
        BaselineFilter, DTAExtractor, MapNormalizer, MzMLSplitter, SpectraFilterWindowMower,
    };
    let lifecycle = |tool: &str| format!("tests/data/topp_cli_lifecycle/write_ini_{tool}.ini");
    assert_written::<BaselineFilter>(&lifecycle(BaselineFilter::NAME), true);
    assert_written::<DTAExtractor>(&lifecycle(DTAExtractor::NAME), true);
    assert_written::<MapNormalizer>(&lifecycle(MapNormalizer::NAME), true);
    assert_written::<MzMLSplitter>(&lifecycle(MzMLSplitter::NAME), true);
    assert_written::<SpectraFilterWindowMower>(&lifecycle(SpectraFilterWindowMower::NAME), true);
}

#[test]
fn write_ini_is_the_retained_text_for_fuzzy_diff() {
    assert_written::<openms::cli::tools::FuzzyDiff>(
        "tests/data/topp_fuzzy_diff/oracle/write_ini.ini",
        false,
    );
}

#[cfg(feature = "mzml")]
#[test]
fn write_ini_is_the_retained_text_for_peak_picker_hi_res() {
    assert_written::<openms::cli::tools::PeakPickerHiRes>(
        "tests/data/topp_peak_picker_hi_res/oracle_write_ini.ini",
        false,
    );
}

#[cfg(all(feature = "mzml", feature = "featurexml"))]
#[test]
fn write_ini_is_the_retained_text_for_file_info_and_feature_finder_centroided() {
    use openms::cli::tools::{FeatureFinderCentroided, FileInfo};
    assert_written::<FileInfo>("tests/data/topp_file_info/expected/write_ini.ini", true);
    let defaults = "tests/data/feature_finder_picked/FeatureFinderCentroided_defaults.ini";
    assert_written::<FeatureFinderCentroided>(defaults, false);
    assert_written::<FeatureFinderCentroided>(defaults, true);
}

/// The trees of the B8c lane's oracle (`../oracle/B8c-numbers`, `driver.cpp`,
/// mode `numbers`): `ParamXMLFile::store` on the C++ Release build at the pins
/// (ibminode06, gcc 14.4.0, two identical runs; `tests/data/b8c_release/ini`,
/// provenance `tests/data/b8c_release_provenance.json`) wrote them.
fn oracle_trees() -> [(&'static str, openms::param::Param); 2] {
    use openms::param::{Param, ParamValue};
    let set = |param: &mut Param, key: &str, value: ParamValue| {
        param.set_value(key, value, "", &[]).unwrap();
    };
    let mut numbers = Param::new();
    set(&mut numbers, "four", ParamValue::Float(4.0));
    set(&mut numbers, "big", ParamValue::Float(1e20));
    set(&mut numbers, "small", ParamValue::Float(1.5e-5));
    set(&mut numbers, "binary", ParamValue::Float(99.632));
    set(&mut numbers, "negative_zero", ParamValue::Float(-0.0));
    set(
        &mut numbers,
        "infinite",
        ParamValue::Float(f64::NEG_INFINITY),
    );
    set(&mut numbers, "huge", ParamValue::Integer(i64::MIN));
    set(
        &mut numbers,
        "list",
        ParamValue::FloatList(vec![4.0, 1e20, 1.5e-5, 0.1, 99.632, 123_456.789, -0.0]),
    );
    numbers.set_min_float("four", 0.0).unwrap();
    numbers.set_max_float("four", 1e20).unwrap();

    let mut specials = Param::new();
    set(&mut specials, "nan", ParamValue::Float(f64::NAN));
    set(&mut specials, "minus_nan", ParamValue::Float(-f64::NAN));
    set(&mut specials, "plus_inf", ParamValue::Float(f64::INFINITY));
    set(&mut specials, "tiny", ParamValue::Float(5e-324));
    set(&mut specials, "largest", ParamValue::Float(f64::MAX));
    set(&mut specials, "third", ParamValue::Float(1.0 / 3.0));
    set(&mut specials, "ten_thousand", ParamValue::Float(1e4));
    set(
        &mut specials,
        "below_ten_thousand",
        ParamValue::Float(9999.999999999999),
    );
    set(&mut specials, "hundredth", ParamValue::Float(0.01));
    set(
        &mut specials,
        "below_hundredth",
        ParamValue::Float(0.009999999999999998),
    );
    set(&mut specials, "half_ulp", ParamValue::Float(0.1 + 0.2));
    set(&mut specials, "int_max", ParamValue::Integer(i64::MAX));
    set(&mut specials, "int_2_32", ParamValue::Integer(1 << 32));
    set(
        &mut specials,
        "specials",
        ParamValue::FloatList(vec![
            f64::NAN,
            -f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            5e-324,
            f64::MAX,
            1.0 / 3.0,
            1e15,
            1e14,
            123_456_789_012_345.6,
            1e-4,
            1e-5,
            0.30000000000000004,
            1e100,
            -2.5e-300,
        ]),
    );
    set(
        &mut specials,
        "ints",
        ParamValue::IntegerList(vec![i32::MIN, -1, 0, i32::MAX]),
    );
    set(
        &mut specials,
        "restricted_list",
        ParamValue::FloatList(vec![0.5]),
    );
    specials.set_min_float("restricted_list", -1e-5).unwrap();
    specials
        .set_max_float("restricted_list", 12345.678)
        .unwrap();
    set(&mut specials, "restricted_int", ParamValue::Integer(3));
    specials.set_min_int("restricted_int", -7).unwrap();
    specials.set_max_int("restricted_int", 9).unwrap();
    set(&mut specials, "open_max", ParamValue::Float(2.0));
    specials.set_max_float("open_max", 1e-3).unwrap();
    [("numbers.ini", numbers), ("specials.ini", specials)]
}

/// Every kind of number the writer spells, byte for byte against the Release
/// build's `ParamXMLFile::store`: fixed items with 15 fractional digits
/// (`99.632000000000005`, `9999.999999999998181`, `0.333333333333333`) and
/// the scientific ones outside `[1e-2, 1e4)` (`1.0e04`, `9.999999999999998e-03`,
/// `5.0e-324`, `1.7976931348623157e308`), `NaN` for either sign, `inf` and
/// `-inf`; list items as `%.15g` (`1e+20`, `4.94065645841247e-324`,
/// `100000000000000`, `123456789012346`, `nan`, `-nan`); restrictions in the
/// item's text (`0.0:1.0e20`, `-1.0e-05:1.2345678e04`, `:1.0e-03`); and
/// 64-bit integers written whole.
#[test]
fn source_number_text_is_the_release_builds() {
    use openms::format::paramxml::{Limits, NumberText, WriteOptions};
    for (file, param) in oracle_trees() {
        let release =
            fs::read_to_string(repository(&format!("tests/data/b8c_release/ini/{file}"))).unwrap();
        let mut bytes = Vec::new();
        paramxml::write_with_options(
            &mut bytes,
            &param,
            Limits::default(),
            WriteOptions::source(),
        )
        .unwrap();
        assert_eq!(String::from_utf8(bytes).unwrap(), release, "{file}");
        // The native text keeps every value, and refuses the integer beyond
        // 32 bits that the source's reader cannot read back.
        let native = WriteOptions {
            numbers: NumberText::RoundTrip,
            ..WriteOptions::source()
        };
        let mut native_bytes = Vec::new();
        assert!(
            paramxml::write_with_options(&mut native_bytes, &param, Limits::default(), native)
                .is_err(),
            "{file}"
        );
    }
    let [(_, mut numbers), _] = oracle_trees();
    numbers.remove("huge").unwrap();
    let mut bytes = Vec::new();
    paramxml::write_with_options(
        &mut bytes,
        &numbers,
        Limits::default(),
        WriteOptions {
            numbers: NumberText::RoundTrip,
            ..WriteOptions::source()
        },
    )
    .unwrap();
    assert_eq!(paramxml::read(bytes.as_slice()).unwrap(), numbers);
}
