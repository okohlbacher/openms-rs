// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The `FeatureFinderLFQ` TOPP tool (`OpenMS4-topp/src/FeatureFinderLFQ.cpp`,
//! topp `174b576`), the wrapper of `Biosaur2Algorithm`.
//!
//! Evidence (`tests/data/topp_feature_finder_lfq_provenance.json`):
//!
//! * **tier 1, retained output**: `TOPP_FeatureFinderLFQ_1`
//!   (test-data `0cb15f2`, `topp/CMakeLists.txt:429-431`) against the retained
//!   `FeatureFinderLFQ_1_output.featureXML`, decoded, under the upstream
//!   `FuzzyDiff` settings with generated ids skipped (decision D6);
//! * **tier 1, executed C++**: the Linux x86_64 Release
//!   `FeatureFinderLFQ` of `openms4-release-bc9cc12-c19e494-174b576`, run twice
//!   on ibminode06 by `../oracle/FF4-biosaur2/remote_run.sh` (identical
//!   outputs): the same invocation, one with both TSV exports, one with
//!   `-write_hills` and default parameters, a FAIMS input and a `mass_traces`
//!   run on a second mzML, and `-write_ini`. Each output is compared exactly:
//!   the featureXML decoded with no tolerance (unique ids included, which the
//!   test-mode generator reproduces), the TSV files byte for byte;
//! * **tier 3, source review**: the error branches, cited to the source.
//!
//! Every case runs the built executable in its own temporary directory.

#![cfg(all(feature = "mzml", feature = "paramxml", feature = "featurexml"))]

#[path = "support/decoded_compare.rs"]
mod decoded;
#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;

use openms::cli::ExitCode;
use openms::format::file_handler::FileHandler;
use openms::format::file_types::FileType;
use openms::kernel::FeatureMap;
use openms::metadata::ProcessingAction;
use openms::system::file::TempDir;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn repository(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

fn fixture(name: &str) -> PathBuf {
    repository("tests/data/topp_feature_finder_lfq").join(name)
}

/// The retained upstream input of `TOPP_FeatureFinderLFQ_1`
/// (`FeatureFinderCentroided_1_input.mzML`), reused from `tests/data/mzml_mobility`.
fn ffc1_input() -> PathBuf {
    repository("tests/data/mzml_mobility/FeatureFinderCentroided_1_input.mzML")
}

struct Workdir(TempDir);

impl Workdir {
    fn new() -> Self {
        Self(TempDir::new_in(std::env::temp_dir(), false).unwrap())
    }
    fn path(&self) -> &Path {
        self.0.path()
    }
}

struct Outcome {
    code: i32,
    out: String,
    err: String,
}

impl Outcome {
    fn assert_exit(&self, expected: ExitCode) {
        assert_eq!(
            self.code,
            expected.as_i32(),
            "expected {}\nstdout:\n{}\nstderr:\n{}",
            expected.name(),
            self.out,
            self.err
        );
    }
}

fn run_in(dir: &Workdir, args: &[&str]) -> Outcome {
    let output = Command::new(env!("CARGO_BIN_EXE_FeatureFinderLFQ"))
        .args(args)
        .current_dir(dir.path())
        .output()
        .expect("the FeatureFinderLFQ executable runs");
    Outcome {
        code: output.status.code().expect("the tool exits"),
        out: String::from_utf8_lossy(&output.stdout).into_owned(),
        err: String::from_utf8_lossy(&output.stderr).into_owned(),
    }
}

fn load(path: &Path) -> FeatureMap {
    FileHandler::load_feature_map(path, &[FileType::FeatureXml]).unwrap()
}

fn assert_exactly_release(actual: &Path, release: &str) {
    let options = decoded::DecodedOptions::new(decoded::Tolerance::exact());
    if let Err(mismatch) =
        decoded::compare_feature_maps(&load(actual), &load(&fixture(release)), &options)
    {
        panic!("{}: {mismatch}", actual.display());
    }
}

fn assert_same_bytes(actual: &Path, release: &str) {
    assert_eq!(
        fs::read_to_string(actual).unwrap(),
        fs::read_to_string(fixture(release)).unwrap(),
        "{}",
        actual.display()
    );
}

const TOPP_1_PARAMETERS: [&str; 12] = [
    "-mini", "20", "-htol", "100", "-itol", "100", "-minlh", "2", "-cmin", "2", "-cmax", "2",
];

#[test]
fn topp_feature_finder_lfq_1_matches_the_retained_expectation() {
    // topp/CMakeLists.txt:429-431
    let dir = Workdir::new();
    let input = ffc1_input();
    let mut args = vec![
        "-test",
        "-in",
        input.to_str().unwrap(),
        "-out",
        "FeatureFinderLFQ_1.tmp.featureXML",
    ];
    args.extend(TOPP_1_PARAMETERS);
    let outcome = run_in(&dir, &args);
    outcome.assert_exit(ExitCode::ExecutionOk);
    let written = dir.path().join("FeatureFinderLFQ_1.tmp.featureXML");

    // The upstream comparison: FuzzyDiff with `-whitelist "id="`, decoded.
    let options = decoded::DecodedOptions::new(decoded::Tolerance::from_settings(
        &fuzzy::FuzzyDiffSettings::upstream().unwrap(),
    ))
    .ignoring_unique_ids();
    if let Err(mismatch) = decoded::compare_feature_maps(
        &load(&written),
        &load(&fixture("FeatureFinderLFQ_1_output.featureXML")),
        &options,
    ) {
        panic!("TOPP_FeatureFinderLFQ_1: {mismatch}");
    }
    // And exactly what the Release build wrote for the same call.
    assert_exactly_release(&written, "release_lfq1.featureXML");

    let map = load(&written);
    assert_eq!(map.features.len(), 20);
    assert_eq!(
        map.primary_ms_run_path().unwrap(),
        ["file://FeatureFinderCentroided_1_input.mzML"]
    );
    assert_eq!(map.data_processing.len(), 1);
    assert_eq!(map.data_processing[0].software.name, "FeatureFinderLFQ");
    assert!(
        map.data_processing[0]
            .actions
            .contains(&ProcessingAction::Quantitation)
    );
    assert!(
        outcome
            .out
            .contains("Wrote 20 features to: FeatureFinderLFQ_1.tmp.featureXML")
    );
    assert!(
        outcome
            .out
            .contains("Hill splitting for group (FAIMS CV=nan) produced 133 hills")
    );
    assert!(
        outcome
            .out
            .contains("After RT apex filter (hrttol=10 s), 46 potential isotope patterns remain.")
    );
}

#[test]
fn both_tsv_exports_match_the_release_build() {
    let dir = Workdir::new();
    let input = ffc1_input();
    let mut args = vec![
        "-test",
        "-in",
        input.to_str().unwrap(),
        "-out",
        "lfq1b.featureXML",
        "-out_tsv",
        "lfq1b.tsv",
        "-out_hills",
        "lfq1b.hills.tsv",
    ];
    args.extend(TOPP_1_PARAMETERS);
    let outcome = run_in(&dir, &args);
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert_exactly_release(
        &dir.path().join("lfq1b.featureXML"),
        "release_lfq1.featureXML",
    );
    assert_same_bytes(&dir.path().join("lfq1b.tsv"), "release_lfq1b.tsv");
    assert_same_bytes(
        &dir.path().join("lfq1b.hills.tsv"),
        "release_lfq1b.hills.tsv",
    );
    assert!(outcome.out.contains("Wrote 133 hills to: lfq1b.hills.tsv"));
    assert!(
        outcome
            .out
            .contains("Wrote 20 features to TSV file: lfq1b.tsv")
    );
}

#[test]
fn write_hills_places_the_hills_next_to_the_output() {
    // FeatureFinderLFQ.cpp:174-187: the hills file is `out` up to its last
    // `.` followed by `.hills.tsv`.
    let dir = Workdir::new();
    let input = ffc1_input();
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-in",
            input.to_str().unwrap(),
            "-out",
            "lfq_def.featureXML",
            "-write_hills",
        ],
    );
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert_exactly_release(
        &dir.path().join("lfq_def.featureXML"),
        "release_lfq_def.featureXML",
    );
    assert_same_bytes(
        &dir.path().join("lfq_def.hills.tsv"),
        "release_lfq_def.hills.tsv",
    );
}

#[test]
fn faims_input_matches_the_release_build() {
    let dir = Workdir::new();
    let input = repository("tests/data/mzml_mobility/FAIMS_test_data.mzML");
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-in",
            input.to_str().unwrap(),
            "-out",
            "lfq_faims.featureXML",
            "-out_tsv",
            "lfq_faims.tsv",
            "-minmz",
            "100",
            "-maxmz",
            "3000",
        ],
    );
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert_exactly_release(
        &dir.path().join("lfq_faims.featureXML"),
        "release_lfq_faims.featureXML",
    );
    assert_same_bytes(&dir.path().join("lfq_faims.tsv"), "release_lfq_faims.tsv");
    assert!(
        outcome
            .out
            .contains("Group (FAIMS CV=-65) has 1 / 1 spectra with ion-mobility data.")
    );
}

#[test]
fn mass_trace_hulls_match_the_release_build() {
    let dir = Workdir::new();
    let input = repository("tests/data/feature_finder_picked/FeatureFinderAlgorithmPicked.mzML");
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-in",
            input.to_str().unwrap(),
            "-out",
            "lfq_ffap.featureXML",
            "-out_tsv",
            "lfq_ffap.tsv",
            "-convex_hulls",
            "mass_traces",
            "-minmz",
            "100",
            "-maxmz",
            "3000",
        ],
    );
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert_exactly_release(
        &dir.path().join("lfq_ffap.featureXML"),
        "release_lfq_ffap.featureXML",
    );
    assert_same_bytes(&dir.path().join("lfq_ffap.tsv"), "release_lfq_ffap.tsv");
}

#[test]
fn write_ini_matches_the_release_build() {
    let dir = Workdir::new();
    let outcome = run_in(&dir, &["-write_ini", "lfq.ini"]);
    outcome.assert_exit(ExitCode::ExecutionOk);
    let written = dir.path().join("lfq.ini");
    let expected = fixture("release_lfq_defaults.ini");
    let mut comparator = fuzzy::FuzzyStringComparator::new();
    comparator.set_acceptable_relative(1.0);
    comparator.set_acceptable_absolute(0.0);
    comparator.set_whitelist(vec!["version".to_owned()]);
    comparator.set_log_destination(fuzzy::LogDestination::Buffer);
    assert!(
        comparator.compare_files(&written, &expected),
        "{}",
        String::from_utf8_lossy(comparator.log())
    );
    assert_eq!(
        openms::format::paramxml::load(&written).unwrap(),
        openms::format::paramxml::load(&expected).unwrap()
    );
}

#[test]
fn an_input_without_ms1_spectra_is_the_algorithms_exception() {
    // Biosaur2Algorithm.cpp:221-224 throws InvalidValue, which TOPPBase reports
    // through its BaseException arm (TOPPBase.cpp:495-499 at cli c19e494).
    let dir = Workdir::new();
    let mzml = r#"<?xml version="1.0" encoding="UTF-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">
  <cvList count="1"><cv id="MS" fullName="PSI-MS" URI="https://example.org/psi-ms.obo"/></cvList>
  <run id="r">
    <spectrumList count="0"/>
  </run>
</mzML>
"#;
    fs::write(dir.path().join("empty.mzML"), mzml).unwrap();
    let outcome = run_in(
        &dir,
        &["-test", "-in", "empty.mzML", "-out", "x.featureXML"],
    );
    outcome.assert_exit(ExitCode::UnknownError);
    assert!(
        outcome.err.contains(
            "Error: Unexpected internal error (the value '' was used but is not valid; No MS1 spectra found in input experiment!)"
        ),
        "{}",
        outcome.err
    );
    assert!(!dir.path().join("x.featureXML").exists());
}

#[test]
fn an_unwritable_tsv_is_a_write_failure() {
    // FeatureFinderLFQ.cpp:197-207 via Biosaur2Algorithm::writeTSV's
    // UnableToCreateFile (Biosaur2Algorithm.cpp:3299-3303).
    let dir = Workdir::new();
    let input = ffc1_input();
    let mut args = vec![
        "-test",
        "-in",
        input.to_str().unwrap(),
        "-out",
        "ok.featureXML",
        "-out_tsv",
        "missing_dir/x.tsv",
    ];
    args.extend(TOPP_1_PARAMETERS);
    let outcome = run_in(&dir, &args);
    outcome.assert_exit(ExitCode::CannotWriteOutputFile);
}
