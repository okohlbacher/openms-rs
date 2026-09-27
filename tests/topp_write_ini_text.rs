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

/// The text of each kind of number, on one parameter tree: the source's
/// (`WriteOptions::source()`) against the lossless native text.
#[test]
fn source_number_text_against_the_native_text() {
    use openms::format::paramxml::{Limits, NumberText, WriteOptions};
    use openms::param::{Param, ParamValue};
    let mut param = Param::new();
    let mut set = |key: &str, value: ParamValue| param.set_value(key, value, "", &[]).unwrap();
    set("four", ParamValue::Float(4.0));
    set("big", ParamValue::Float(1e20));
    set("small", ParamValue::Float(1.5e-5));
    set("binary", ParamValue::Float(99.632));
    set("negative_zero", ParamValue::Float(-0.0));
    set("infinite", ParamValue::Float(f64::NEG_INFINITY));
    set("huge", ParamValue::Integer(i64::MIN));
    set(
        "list",
        ParamValue::FloatList(vec![4.0, 1e20, 1.5e-5, 0.1, 99.632, 123_456.789, -0.0]),
    );
    param.set_min_float("four", 0.0).unwrap();
    param.set_max_float("four", 1e20).unwrap();
    let text = |param: &Param, options: WriteOptions| {
        let mut bytes = Vec::new();
        paramxml::write_with_options(&mut bytes, param, Limits::default(), options)
            .map(|()| String::from_utf8(bytes).unwrap())
    };
    let source = text(&param, WriteOptions::source()).unwrap();
    for expected in [
        "<ITEM name=\"four\" value=\"4.0\" type=\"double\" description=\"\" required=\"false\" advanced=\"false\" restrictions=\"0.0:1.0e20\" />",
        "<ITEM name=\"big\" value=\"1.0e20\" type=\"double\"",
        "<ITEM name=\"small\" value=\"1.5e-05\" type=\"double\"",
        "<ITEM name=\"binary\" value=\"99.632000000000005\" type=\"double\"",
        "<ITEM name=\"negative_zero\" value=\"-0.0\" type=\"double\"",
        "<ITEM name=\"infinite\" value=\"-inf\" type=\"double\"",
        "<ITEM name=\"huge\" value=\"-9223372036854775808\" type=\"int\"",
        "<LISTITEM value=\"4\"/>\n",
        "<LISTITEM value=\"1e+20\"/>\n",
        "<LISTITEM value=\"1.5e-05\"/>\n",
        "<LISTITEM value=\"0.1\"/>\n",
        "<LISTITEM value=\"99.632\"/>\n",
        "<LISTITEM value=\"123456.789\"/>\n",
        "<LISTITEM value=\"-0\"/>\n",
    ] {
        assert!(source.contains(expected), "{expected}\n{source}");
    }
    // The native text keeps every value, and refuses the integer the source's
    // reader cannot read back.
    let native = WriteOptions {
        numbers: NumberText::RoundTrip,
        ..WriteOptions::source()
    };
    assert!(text(&param, native).is_err());
    param.remove("huge").unwrap();
    let native = text(&param, native).unwrap();
    assert_eq!(paramxml::read(native.as_bytes()).unwrap(), param);
}
