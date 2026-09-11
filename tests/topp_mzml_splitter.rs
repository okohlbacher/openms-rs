// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The upstream `TOPP_MzMLSplitter_*` tests, compared canonically against the
//! retained C++ outputs.
//!
//! The upstream test uses FuzzyDiff, not a byte comparison, so byte equality is
//! not the contract for mzML here either. Per `docs/DIFFERENTIAL_VALIDATION.md`
//! this compares decoded content: how many records each part received, which
//! records they are, and their peak values.

use openms::cli::tools::MzMLSplitter;
use openms::cli::{ExitCode, run_with};
use openms::format::file_handler::FileHandler;
use openms::format::file_types::FileType;
use openms::kernel::{MSExperiment, MSSpectrum};
use std::fs;
use std::path::{Path, PathBuf};

fn fixture(name: &str) -> PathBuf {
    Path::new("tests/data").join(name)
}

fn split(case: &str, args: &[&str]) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms-split-{}-{case}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let mut arguments = vec![
        "MzMLSplitter".to_string(),
        "-test".to_string(),
        "-in".to_string(),
        fixture("mzml_splitter_input.mzML")
            .to_string_lossy()
            .into_owned(),
        "-out".to_string(),
        dir.join("out").to_string_lossy().into_owned(),
    ];
    arguments.extend(args.iter().map(|a| (*a).to_string()));
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<MzMLSplitter>(&arguments, &mut out, &mut err);
    assert_eq!(
        code,
        ExitCode::ExecutionOk,
        "tool failed: {}",
        String::from_utf8_lossy(&err)
    );
    dir
}

/// Records of two runs agree when they hold the same spectra in the same order,
/// with equal native identifiers, MS levels, retention times and peak values.
fn assert_same_records(actual: &MSExperiment, expected: &MSExperiment, what: &str) {
    assert_eq!(
        actual.spectra.len(),
        expected.spectra.len(),
        "{what}: spectrum count"
    );
    assert_eq!(
        actual.chromatograms.len(),
        expected.chromatograms.len(),
        "{what}: chromatogram count"
    );
    for (i, (a, e)) in actual
        .spectra
        .iter()
        .zip(expected.spectra.iter())
        .enumerate()
    {
        assert_eq!(a.native_id, e.native_id, "{what}: spectrum {i} identifier");
        assert_eq!(a.ms_level, e.ms_level, "{what}: spectrum {i} MS level");
        assert!(
            (a.rt - e.rt).abs() <= 1e-9,
            "{what}: spectrum {i} retention time {} vs {}",
            a.rt,
            e.rt
        );
        assert_peaks(a, e, &format!("{what}: spectrum {i}"));
    }
}

fn assert_peaks(a: &MSSpectrum, e: &MSSpectrum, what: &str) {
    assert_eq!(a.peaks.len(), e.peaks.len(), "{what}: peak count");
    for (i, (pa, pe)) in a.peaks.iter().zip(e.peaks.iter()).enumerate() {
        assert!(
            (pa.mz - pe.mz).abs() <= 1e-9,
            "{what}: peak {i} m/z {} vs {}",
            pa.mz,
            pe.mz
        );
        assert!(
            (f64::from(pa.intensity) - f64::from(pe.intensity)).abs() <= 1e-6,
            "{what}: peak {i} intensity {} vs {}",
            pa.intensity,
            pe.intensity
        );
    }
}

fn load(path: impl AsRef<Path>) -> MSExperiment {
    FileHandler::load_experiment(path, &[FileType::MzMl]).unwrap()
}

#[test]
fn topp_mzml_splitter_1_splits_into_a_requested_number_of_parts() {
    let dir = split("1", &["-parts", "2"]);
    for part in 1..=2 {
        let produced = load(dir.join(format!("out_part{part}of2.mzML")));
        let reference = load(fixture(&format!("mzml_splitter_output_part{part}.mzML")));
        assert_same_records(&produced, &reference, &format!("part {part}"));
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn topp_mzml_splitter_2_derives_the_part_count_from_a_size_limit() {
    // 40 KB over a ~59 KB input rounds up to the same two parts as case 1.
    let dir = split("2", &["-size", "40", "-unit", "KB"]);
    for part in 1..=2 {
        let produced = load(dir.join(format!("out_part{part}of2.mzML")));
        let reference = load(fixture(&format!("mzml_splitter_output_part{part}.mzML")));
        assert_same_records(&produced, &reference, &format!("part {part}"));
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn every_record_is_placed_exactly_once() {
    // The remainder is spread over the parts still to come, so no record is
    // dropped or duplicated for a part count that does not divide evenly.
    let whole = load(fixture("mzml_splitter_input.mzML"));
    // parts=1 with no size limit is refused by the source, so it starts at 2.
    for parts in [2usize, 3, 4, 7, 11] {
        let dir = split(&format!("n{parts}"), &["-parts", &parts.to_string()]);
        let width = parts.to_string().len();
        let mut seen = Vec::new();
        for part in 1..=parts {
            let path = dir.join(format!("out_part{part:0width$}of{parts}.mzML"));
            seen.extend(load(&path).spectra.into_iter().map(|s| s.native_id));
        }
        let expected: Vec<String> = whole.spectra.iter().map(|s| s.native_id.clone()).collect();
        assert_eq!(seen, expected, "parts={parts} lost or reordered records");
        let _ = fs::remove_dir_all(&dir);
    }
}

#[test]
fn conflicting_and_missing_options_are_rejected() {
    let run = |args: &[&str]| -> ExitCode {
        let arguments: Vec<String> = std::iter::once("MzMLSplitter".to_string())
            .chain(args.iter().map(|a| (*a).to_string()))
            .collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        run_with::<MzMLSplitter>(&arguments, &mut out, &mut err)
    };
    let input = fixture("mzml_splitter_input.mzML")
        .to_string_lossy()
        .into_owned();
    // Source refuses both record filters at once.
    assert_eq!(
        run(&["-in", &input, "-no_chrom", "-no_spec"]),
        ExitCode::IllegalParameters
    );
    // One part and no size limit leaves nothing to do.
    assert_eq!(run(&["-in", &input]), ExitCode::IllegalParameters);
    // The unit is restricted to the registered set.
    assert_eq!(
        run(&["-in", &input, "-size", "1", "-unit", "TB"]),
        ExitCode::IllegalParameters
    );
    // parts has a registered minimum.
    assert_eq!(
        run(&["-in", &input, "-parts", "0"]),
        ExitCode::IllegalParameters
    );
}
