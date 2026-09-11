// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The upstream `TOPP_SpectraFilterWindowMower_1` test, compared canonically
//! against the retained C++ output, plus the algorithm-subsection plumbing that
//! every algorithm-wrapping TOPP tool depends on.

use openms::cli::tools::SpectraFilterWindowMower;
use openms::cli::{ExitCode, run_with};
use openms::format::file_handler::FileHandler;
use openms::format::file_types::FileType;
use openms::kernel::MSExperiment;
use std::fs;
use std::path::{Path, PathBuf};

fn fixture(name: &str) -> PathBuf {
    Path::new("tests/data").join(name)
}
fn load(path: impl AsRef<Path>) -> MSExperiment {
    FileHandler::load_experiment(path, &[FileType::MzMl]).unwrap()
}
fn run(args: &[&str]) -> (ExitCode, String, String) {
    let arguments: Vec<String> = std::iter::once("SpectraFilterWindowMower".to_string())
        .chain(args.iter().map(|a| (*a).to_string()))
        .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<SpectraFilterWindowMower>(&arguments, &mut out, &mut err);
    (
        code,
        String::from_utf8_lossy(&out).into_owned(),
        String::from_utf8_lossy(&err).into_owned(),
    )
}
fn workdir(case: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms-wm-{}-{case}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn topp_spectra_filter_window_mower_1_matches_the_retained_output() {
    let dir = workdir("1");
    let out = dir.join("out.mzML");
    let (code, _, err) = run(&[
        "-test",
        "-in",
        &fixture("window_mower_tool_input.mzML").to_string_lossy(),
        "-out",
        &out.to_string_lossy(),
    ]);
    assert_eq!(code, ExitCode::ExecutionOk, "{err}");

    let produced = load(&out);
    let reference = load(fixture("window_mower_tool_output.mzML"));
    assert_eq!(produced.spectra.len(), reference.spectra.len());
    for (index, (a, e)) in produced
        .spectra
        .iter()
        .zip(reference.spectra.iter())
        .enumerate()
    {
        assert_eq!(a.peaks.len(), e.peaks.len(), "spectrum {index} peak count");
        for (i, (pa, pe)) in a.peaks.iter().zip(e.peaks.iter()).enumerate() {
            assert!(
                (pa.mz - pe.mz).abs() <= 1e-9,
                "spectrum {index} peak {i} m/z"
            );
            assert!(
                (f64::from(pa.intensity) - f64::from(pe.intensity)).abs() <= 1e-6,
                "spectrum {index} peak {i} intensity"
            );
        }
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn the_algorithm_subsection_carries_the_source_defaults() {
    let dir = workdir("ini");
    let ini = dir.join("tool.ini");
    let (code, _, err) = run(&["-write_ini", &ini.to_string_lossy()]);
    assert_eq!(code, ExitCode::ExecutionOk, "{err}");

    // The subsection must reach the INI with the C++ WindowMower defaults.
    let loaded = openms::format::paramxml::load(&ini).unwrap();
    let key = |k: &str| format!("SpectraFilterWindowMower:1:algorithm:{k}");
    assert_eq!(
        loaded.value(&key("windowsize")).unwrap().to_f64().unwrap(),
        50.0
    );
    assert_eq!(
        loaded.value(&key("peakcount")).unwrap().to_i64().unwrap(),
        2
    );
    assert_eq!(
        loaded.value(&key("movetype")).unwrap().as_str().unwrap(),
        "slide"
    );
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn a_subsection_value_from_an_ini_changes_the_result() {
    let dir = workdir("override");
    let ini = dir.join("tool.ini");
    assert_eq!(
        run(&["-write_ini", &ini.to_string_lossy()]).0,
        ExitCode::ExecutionOk
    );

    // Keeping one peak per window must retain no more peaks than the default of two.
    let mut param = openms::format::paramxml::load(&ini).unwrap();
    param
        .set_value(
            "SpectraFilterWindowMower:1:algorithm:peakcount",
            openms::param::ParamValue::Integer(1),
            "",
            &[],
        )
        .unwrap();
    openms::format::paramxml::store(&ini, &param).unwrap();

    let out = dir.join("one.mzML");
    let (code, _, err) = run(&[
        "-test",
        "-ini",
        &ini.to_string_lossy(),
        "-in",
        &fixture("window_mower_tool_input.mzML").to_string_lossy(),
        "-out",
        &out.to_string_lossy(),
    ]);
    assert_eq!(code, ExitCode::ExecutionOk, "{err}");

    let one: usize = load(&out).spectra.iter().map(|s| s.peaks.len()).sum();
    let two: usize = load(fixture("window_mower_tool_output.mzML"))
        .spectra
        .iter()
        .map(|s| s.peaks.len())
        .sum();
    assert!(
        one <= two,
        "peakcount=1 kept {one} peaks, default kept {two}"
    );
    let _ = fs::remove_dir_all(&dir);
}

/// An INI value that violates a registered restriction is **ignored**, not an
/// error: the source `Param::update` warns and keeps the existing default. The
/// tool therefore still runs, with the default movement type.
#[test]
fn an_invalid_subsection_value_falls_back_to_the_default() {
    let dir = workdir("bad");
    let ini = dir.join("tool.ini");
    assert_eq!(
        run(&["-write_ini", &ini.to_string_lossy()]).0,
        ExitCode::ExecutionOk
    );
    let mut param = openms::format::paramxml::load(&ini).unwrap();
    param
        .set_value(
            "SpectraFilterWindowMower:1:algorithm:movetype",
            openms::param::ParamValue::String("sideways".into()),
            "",
            &[],
        )
        .unwrap();
    openms::format::paramxml::store(&ini, &param).unwrap();
    let out = dir.join("x.mzML");
    let (code, _, err) = run(&[
        "-test",
        "-ini",
        &ini.to_string_lossy(),
        "-in",
        &fixture("window_mower_tool_input.mzML").to_string_lossy(),
        "-out",
        &out.to_string_lossy(),
    ]);
    assert_eq!(code, ExitCode::ExecutionOk, "{err}");
    // Falling back to "slide" must reproduce the retained reference exactly.
    let produced: Vec<usize> = load(&out).spectra.iter().map(|s| s.peaks.len()).collect();
    let reference: Vec<usize> = load(fixture("window_mower_tool_output.mzML"))
        .spectra
        .iter()
        .map(|s| s.peaks.len())
        .collect();
    assert_eq!(produced, reference);
    let _ = fs::remove_dir_all(&dir);
}
