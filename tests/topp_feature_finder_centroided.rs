// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The FeatureFinderCentroided TOPP wrapper: registration, `-write_ini`,
//! loading, every error branch, the FAIMS refusal and the output annotation
//! (package C5-FFC-WRAPPER of the early TOPP bundle).
//!
//! Evidence:
//!
//! * **Oracle cases.** Every case whose name is quoted in a comment as
//!   `FFC_...`, `TOPP_FeatureFinderCentroided_1` or `c5_...` was executed
//!   against the product SDK (Debug, OpenMS core 4fdec46): the planned wrapper
//!   regressions in `../oracle/topp-early-bundle` (C1, manifest
//!   `64e98543...`) and the branch-order cases in `../oracle/ffc-wrapper-c5`
//!   (this package's supplement). Each asserts that run's exit code and its
//!   `Error:` line. Evidence: oracle-generated (tier 1 executed differential).
//!   Two C1 cases exit through a Debug-only precondition
//!   (`FFC_FileFilter_44_force`, `FFC_im_arrays_ms2_only`); their exit codes are
//!   never asserted, only that this wrapper does not refuse the input.
//! * **Source review.** Registration, the `-write_ini` defaults and the
//!   annotation steps are transcribed from `FeatureFinderCentroided.cpp` at
//!   topp `174b576` (tier 3), and checked against the executed `-write_ini`
//!   output where one exists.
//!
//! The algorithm stops after seed selection until package B7 lands, so every
//! input that passes the wrapper ends in the documented `Error: unsupported:`
//! line with exit 11 and no output; `TOPP_FeatureFinderCentroided_1` itself is
//! package B10's. What the wrapper does after the algorithm is tested directly
//! through `FeatureFinderCentroided::finish_features`.
//!
//! Synthetic inputs are not committed: each is derived here from the retained
//! `FeatureFinderCentroided_1_input.mzML` by the rule the oracle manifests
//! record, and the derived bytes are checked against the SHA-1 of the file the
//! C++ tool actually read. `tests/data/topp_feature_finder_centroided_provenance.json`
//! records every hash.
//!
//! Every case runs the built executable in its own temporary working
//! directory, because the `FeatureFinderCentroided_1` INI sets `log=TOPP.log`.

// The tool exists under the features its executable requires.
#![cfg(all(feature = "mzml", feature = "paramxml", feature = "featurexml"))]

#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;

use base64::Engine;
use openms::cli::tools::FeatureFinderCentroided;
use openms::cli::{
    ExitCode, TEST_MODE_COMPLETION_TIME, TEST_MODE_PARAMETER_KEY, TEST_MODE_PARAMETER_VALUE,
    TEST_MODE_UNIQUE_ID_SEED, TEST_MODE_VERSION, Tool, ToolContext, ToolSpec, run_with, tool_spec,
};
use openms::concept::UniqueIdGenerator;
use openms::format::file_handler::FileHandler;
use openms::format::file_types::FileType;
use openms::kernel::features::{BaseFeature, Feature, FeatureMap};
use openms::kernel::geometry::{ConvexHull2D, Point2D};
use openms::kernel::{MSExperiment, MSSpectrum};
use openms::metadata::{MetaValue, ProcessingAction};
use openms::param::ParamValue;
use openms::system::file::TempDir;
use openms::{Error, Result};
use std::cell::RefCell;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

// ---------------------------------------------------------------------------
// Paths, running and small helpers
// ---------------------------------------------------------------------------

fn repository(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

/// An A3-owned mobility fixture, reused read-only (`tests/data/mzml_mobility`).
fn mobility(name: &str) -> PathBuf {
    repository("tests/data/mzml_mobility").join(name)
}

/// A B6-owned feature-finder fixture, reused read-only.
fn picked(name: &str) -> PathBuf {
    repository("tests/data/feature_finder_picked").join(name)
}

/// A fixture of this package.
fn fixture(name: &str) -> PathBuf {
    repository("tests/data/topp_feature_finder_centroided").join(name)
}

/// The retained upstream input of `TOPP_FeatureFinderCentroided_1`
/// (test-data `0cb15f2`, sha256 `a3dfae63...`), reused from A3's fixtures.
fn ffc1_input() -> PathBuf {
    mobility("FeatureFinderCentroided_1_input.mzML")
}

/// The retained `FeatureFinderCentroided_1` INI (`log=TOPP.log`, `threads=1`,
/// version 3.6.0), reused from B6's fixtures.
fn ffc1_ini() -> PathBuf {
    picked("FeatureFinderCentroided_1_parameters.ini")
}

fn text(path: impl AsRef<Path>) -> String {
    path.as_ref().to_string_lossy().into_owned()
}

/// A fresh, uniquely named working directory per case, removed when it ends.
struct Workdir(TempDir);

impl Workdir {
    fn new() -> Self {
        Self(TempDir::new_in(std::env::temp_dir(), false).unwrap())
    }
    fn path(&self) -> &Path {
        self.0.path()
    }
    fn file(&self, name: &str) -> String {
        text(self.path().join(name))
    }
    /// Write `data` into the directory and return its path.
    fn put(&self, name: &str, data: &[u8]) -> String {
        let path = self.path().join(name);
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).unwrap();
        }
        fs::write(&path, data).unwrap();
        text(path)
    }
}

/// What one executable run produced.
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
    fn assert_err_contains(&self, needle: &str) {
        assert!(
            self.err.contains(needle),
            "stderr does not contain {needle:?}:\n{}",
            self.err
        );
    }
    fn assert_out_contains(&self, needle: &str) {
        assert!(
            self.out.contains(needle),
            "stdout does not contain {needle:?}:\n{}",
            self.out
        );
    }
}

/// Run the built `FeatureFinderCentroided` executable inside `dir`.
///
/// The working directory matters: the `FeatureFinderCentroided_1` INI sets
/// `log=TOPP.log`, which the C++ tool writes into the current directory.
fn run_in(dir: &Workdir, args: &[&str]) -> Outcome {
    let output = Command::new(env!("CARGO_BIN_EXE_FeatureFinderCentroided"))
        .args(args)
        .current_dir(dir.path())
        .output()
        .expect("the FeatureFinderCentroided executable runs");
    Outcome {
        code: output
            .status
            .code()
            .expect("the tool exits, it is not signalled"),
        out: String::from_utf8_lossy(&output.stdout).into_owned(),
        err: String::from_utf8_lossy(&output.stderr).into_owned(),
    }
}

/// The message the picked algorithm returns until package B7 lands, as the
/// framework reports an [`Error::Unsupported`].
const ALGORITHM_NOT_PORTED: &str = "Error: unsupported: FeatureFinderAlgorithmPicked seed extension, trace fitting and feature resolution (source step 3.3 onward) are not ported yet";

/// Assert that the wrapper accepted the input and handed it to the algorithm,
/// which is as far as this package's port runs.
fn assert_reached_the_algorithm(outcome: &Outcome) {
    outcome.assert_out_contains(FeatureFinderCentroided::NO_FAIMS_MESSAGE);
    assert!(
        !outcome.err.contains("Profile data provided"),
        "the profile check refused the input:\n{}",
        outcome.err
    );
    assert!(
        !outcome.err.contains("per-peak ion mobility"),
        "the ion-mobility check refused the input:\n{}",
        outcome.err
    );
}

// ---------------------------------------------------------------------------
// Derived inputs
//
// The C++ tool read files derived from the retained FFC_1 input by the rules
// recorded in ../oracle/topp-early-bundle/cases/derived.json and
// ../oracle/ffc-wrapper-c5/manifest.json. The same rules are applied here and
// the result is checked against the SHA-1 of the file that was executed, so no
// derived megabyte enters the repository and the identity is still pinned.
// ---------------------------------------------------------------------------

const REPRESENTATION_CV: &[u8] =
    br#"<cvParam cvRef="MS" accession="MS:1000525" name="spectrum representation" />"#;
const PROFILE_CV: &[u8] =
    br#"<cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum" />"#;
const SPECTRA: usize = 112;

fn find(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    haystack
        .windows(needle.len())
        .position(|window| window == needle)
}

fn count(haystack: &[u8], needle: &[u8]) -> usize {
    haystack
        .windows(needle.len())
        .filter(|window| *window == needle)
        .count()
}

fn lines(source: &[u8]) -> Vec<&[u8]> {
    source.split(|byte| *byte == b'\n').collect()
}

fn join(lines: Vec<Vec<u8>>) -> Vec<u8> {
    let mut out = Vec::new();
    for (index, line) in lines.into_iter().enumerate() {
        if index != 0 {
            out.push(b'\n');
        }
        out.extend_from_slice(&line);
    }
    out
}

fn trimmed(line: &[u8]) -> &[u8] {
    let start = line
        .iter()
        .position(|b| !b.is_ascii_whitespace())
        .unwrap_or(line.len());
    let end = line
        .iter()
        .rposition(|b| !b.is_ascii_whitespace())
        .map_or(start, |index| index + 1);
    &line[start..end]
}

/// SHA-1 of `data`, used only to pin a derived input to the file the C++ tool
/// read; the SHA-256 of both is in the provenance manifest.
fn digest(data: &[u8]) -> String {
    use sha1::{Digest, Sha1};
    let mut hasher = Sha1::new();
    hasher.update(data);
    hasher
        .finalize()
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

fn assert_digest(data: &[u8], expected: &str, case: &str) {
    assert_eq!(digest(data), expected, "derived input {case}");
}

/// Replace the `MS:1000525` term of the spectra `select` accepts (0-based) by
/// the `MS:1000128` profile term.
fn with_profile_terms(source: &[u8], select: impl Fn(usize) -> bool) -> Vec<u8> {
    assert_eq!(count(source, REPRESENTATION_CV), SPECTRA);
    let mut out = Vec::with_capacity(source.len());
    let mut rest = source;
    let mut index = 0;
    while let Some(at) = find(rest, REPRESENTATION_CV) {
        out.extend_from_slice(&rest[..at]);
        out.extend_from_slice(if select(index) {
            PROFILE_CV
        } else {
            REPRESENTATION_CV
        });
        rest = &rest[at + REPRESENTATION_CV.len()..];
        index += 1;
    }
    out.extend_from_slice(rest);
    out
}

/// Oracle `profile/FeatureFinderCentroided_1_input.mzML` (C1): every spectrum
/// stored as profile.
fn derive_profile(source: &[u8]) -> Vec<u8> {
    let derived = with_profile_terms(source, |_| true);
    assert_digest(
        &derived,
        "db09cf03b27a8381f3c88ca019bdf7e1db81ca27",
        "profile",
    );
    derived
}

/// Oracle `profile_before_1000525/FeatureFinderCentroided_1_input.mzML` (C1):
/// an unindented profile term before every `MS:1000525` term, which the source
/// reader resets to `UNKNOWN`.
fn derive_profile_then_representation(source: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    for line in lines(source) {
        if trimmed(line) == REPRESENTATION_CV {
            out.push(PROFILE_CV.to_vec());
        }
        out.push(line.to_vec());
    }
    let derived = join(out);
    assert_digest(
        &derived,
        "54288a97441c37647209ef56395c231c8c04aea9",
        "profile_before_1000525",
    );
    derived
}

/// The third binary data array of oracle `im_peak_ms1.mzML` (C1):
/// `count` uncompressed 32-bit values of 0.5 in a mean ion mobility array.
fn im_block(count: usize, with_units: bool) -> Vec<Vec<u8>> {
    let raw: Vec<u8> = std::iter::repeat_n([0x00, 0x00, 0x00, 0x3f], count)
        .flatten()
        .collect();
    let payload = base64::engine::general_purpose::STANDARD.encode(&raw);
    let im_cv: &[u8] = if with_units {
        br#"<cvParam cvRef="MS" accession="MS:1002816" name="mean ion mobility array" unitAccession="UO:0000028" unitName="millisecond" unitCvRef="UO"/>"#
    } else {
        br#"<cvParam cvRef="MS" accession="MS:1002816" name="mean ion mobility array"/>"#
    };
    vec![
        format!(
            "\t\t\t\t<binaryDataArray encodedLength=\"{}\">",
            payload.len()
        )
        .into_bytes(),
        br#"<cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" />"#.to_vec(),
        br#"<cvParam cvRef="MS" accession="MS:1000576" name="no compression" />"#.to_vec(),
        im_cv.to_vec(),
        format!("<binary>{payload}</binary>").into_bytes(),
        b"</binaryDataArray>".to_vec(),
        b"</binaryDataArrayList>".to_vec(),
    ]
}

/// Oracle `im_peak_ms1.mzML` and `im_peak_ms1_nounit.mzML` (C1): a third
/// binary data array with per-peak ion mobility on every spectrum.
fn derive_im_peak(source: &[u8], with_units: bool) -> Vec<u8> {
    let mut out: Vec<Vec<u8>> = Vec::new();
    let mut length = 0usize;
    let mut lists = 0usize;
    for line in lines(source) {
        if let Some(value) = default_array_length(line) {
            length = value;
        }
        if line == b"\t\t\t\t<binaryDataArrayList count=\"2\">" {
            out.push(b"\t\t\t\t<binaryDataArrayList count=\"3\">".to_vec());
            continue;
        }
        if line == b"\t\t\t\t</binaryDataArrayList>" {
            out.extend(im_block(length, with_units));
            lists += 1;
            continue;
        }
        out.push(line.to_vec());
    }
    assert_eq!(lists, SPECTRA);
    let derived = join(out);
    let expected = if with_units {
        "c7dbb4e5a1948ddf1295c49025655cbb6651b5f0"
    } else {
        "7c3d5bea5e87986b3efacc046c536730929721b9"
    };
    assert_digest(&derived, expected, "im_peak_ms1");
    derived
}

/// `defaultArrayLength` of a `<spectrum>` line, if this is one.
fn default_array_length(line: &[u8]) -> Option<usize> {
    if !trimmed(line).starts_with(b"<spectrum ") {
        return None;
    }
    let marker = br#"defaultArrayLength=""#;
    let at = find(line, marker)? + marker.len();
    let digits: Vec<u8> = line[at..]
        .iter()
        .copied()
        .take_while(u8::is_ascii_digit)
        .collect();
    String::from_utf8(digits).ok()?.parse().ok()
}

/// Oracle `faims_one_cv.mzML` and `faims_two_cv.mzML` (C1) and
/// `faims_partial_cv.mzML` (C5): a scan-level FAIMS compensation voltage after
/// every `<scan>` line, taken from `values` in turn; an empty entry adds none.
fn derive_faims(source: &[u8], values: &[&str]) -> Vec<u8> {
    let mut out: Vec<Vec<u8>> = Vec::new();
    let mut scans = 0usize;
    for line in lines(source) {
        out.push(line.to_vec());
        if trimmed(line) == b"<scan>" {
            let value = values[scans % values.len()];
            if !value.is_empty() {
                out.push(format!(
                    r#"<cvParam cvRef="MS" accession="MS:1001581" name="FAIMS compensation voltage" value="{value}" unitAccession="UO:0000218" unitName="volt" unitCvRef="UO"/>"#
                ).into_bytes());
            }
            scans += 1;
        }
    }
    assert_eq!(scans, SPECTRA);
    join(out)
}

fn derive_faims_one_cv(source: &[u8]) -> Vec<u8> {
    let derived = derive_faims(source, &["-45"]);
    assert_digest(
        &derived,
        "44aba90786460e1e167db1b549e04a960db9edd9",
        "faims_one_cv",
    );
    derived
}

fn derive_faims_two_cv(source: &[u8]) -> Vec<u8> {
    let derived = derive_faims(source, &["-45", "-60"]);
    assert_digest(
        &derived,
        "9c0e2c828f30123468799c10169d61e6d728e6e4",
        "faims_two_cv",
    );
    derived
}

/// Oracle `negative_intensities.mzML` (C5): the sign bit of every 32-bit
/// intensity set, so the source load filter drops every peak.
fn derive_negative_intensities(source: &[u8]) -> Vec<u8> {
    let engine = base64::engine::general_purpose::STANDARD;
    let mut out: Vec<Vec<u8>> = Vec::new();
    let mut in_intensity = false;
    let mut touched = 0usize;
    for line in lines(source) {
        if find(line, br#"accession="MS:1000515""#).is_some() {
            in_intensity = true;
        }
        let content = trimmed(line);
        if in_intensity && content.starts_with(b"<binary>") && content.ends_with(b"</binary>") {
            let payload = &content[b"<binary>".len()..content.len() - b"</binary>".len()];
            let mut raw = engine.decode(payload).unwrap();
            assert_eq!(raw.len() % 4, 0);
            for value in raw.chunks_mut(4) {
                assert_eq!(
                    value[3] & 0x80,
                    0,
                    "the FFC_1 input has no negative intensity"
                );
                value[3] |= 0x80;
            }
            let encoded = engine.encode(&raw);
            assert_eq!(encoded.len(), payload.len());
            let indent_end = line
                .iter()
                .position(|byte| !byte.is_ascii_whitespace())
                .unwrap_or(line.len());
            let mut replaced = line[..indent_end].to_vec();
            replaced.extend_from_slice(format!("<binary>{encoded}</binary>").as_bytes());
            out.push(replaced);
            in_intensity = false;
            touched += 1;
            continue;
        }
        out.push(line.to_vec());
    }
    assert_eq!(touched, SPECTRA);
    let derived = join(out);
    assert_digest(
        &derived,
        "8adf8b8a0ec464d69e57ca84ff3316e941cba08e",
        "negative_intensities",
    );
    derived
}

// ---------------------------------------------------------------------------
// Registration and -write_ini
// ---------------------------------------------------------------------------

/// Source `registerOptionsAndFlags_` (`FeatureFinderCentroided.cpp:140-163`)
/// and `getSubsectionDefaults_` (166-169), parameter by parameter.
///
/// Evidence: tier 3 (source review); the values are pinned again by the
/// `-write_ini` comparison below, which is tier 1.
#[test]
fn registration_matches_the_source() {
    let spec = tool_spec::<FeatureFinderCentroided>().unwrap();
    let parameter = |name: &str| {
        spec.find(name)
            .unwrap_or_else(|| panic!("parameter {name} is registered"))
    };

    let input = parameter("in");
    assert_eq!(input.description, "input file");
    assert_eq!(input.argument, "<file>");
    assert!(input.required && !input.advanced);
    assert_eq!(input.valid_formats, ["mzML"]);

    let output = parameter("out");
    assert_eq!(output.description, "output file");
    assert!(output.required && !output.advanced);
    assert_eq!(output.valid_formats, ["featureXML"]);

    let seeds = parameter("seeds");
    assert_eq!(seeds.description, "User specified seed list");
    assert!(!seeds.required && !seeds.advanced);
    assert_eq!(seeds.valid_formats, ["featureXML"]);

    let merge = parameter("faims_merge_features");
    assert_eq!(merge.argument, "<true/false>");
    assert_eq!(merge.default_value, ParamValue::String("true".into()));
    assert_eq!(merge.valid_strings, ["true", "false"]);
    assert!(!merge.required && !merge.advanced);
    assert!(merge.description.starts_with(
        "For FAIMS data with multiple compensation voltages: Merge features representing the same analyte"
    ));
    assert!(
        merge
            .description
            .ends_with("Has no effect on non-FAIMS data.")
    );

    // The source registers exactly one subsection, whose defaults are the
    // picked algorithm's.
    assert_eq!(
        spec.subsections(),
        [("algorithm".to_owned(), "Algorithm section".to_owned())]
    );
    let defaults = FeatureFinderCentroided::subsection_defaults("algorithm")
        .unwrap()
        .expect("the algorithm subsection has defaults");
    assert_eq!(
        defaults,
        openms::analysis::feature_finder_picked::algorithm::default_parameters().unwrap()
    );
    // Any other name returns the same tree, as the source ignores the name.
    assert_eq!(
        FeatureFinderCentroided::subsection_defaults("anything").unwrap(),
        Some(defaults)
    );
}

/// Oracle `FFC_write_ini` and `FFC_write_ini_test`: both wrote the same file
/// (sha256 `2869134a...`), which is B6's retained
/// `FeatureFinderCentroided_defaults.ini`.
///
/// It is compared as the upstream `TOPPWRITEINI_OVERWRITE_out` comparison does,
/// line by line with exact numbers and `version` lines skipped, and then as a
/// decoded parameter tree, entry for entry with descriptions, tags and
/// restrictions. `TOPPWRITEINI_<tool>_SectionName` (test-data
/// `topp/CMakeLists.txt:83-85`) requires the first `  <NODE name="` line to name
/// the tool. Evidence: tier 1 (executed differential).
#[test]
fn write_ini_matches_the_cpp_file() {
    let expected = picked("FeatureFinderCentroided_defaults.ini");
    for arguments in [
        vec!["-write_ini", "FeatureFinderCentroided.tmp.ini"],
        vec!["-test", "-write_ini", "FeatureFinderCentroided.tmp.ini"],
    ] {
        let dir = Workdir::new();
        let outcome = run_in(&dir, &arguments);
        outcome.assert_exit(ExitCode::ExecutionOk);
        let written = dir.path().join("FeatureFinderCentroided.tmp.ini");
        let content = fs::read_to_string(&written).unwrap();
        assert!(
            content.starts_with("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"),
            "{content}"
        );
        let section = content
            .lines()
            .find_map(|line| line.strip_prefix("  <NODE name=\""))
            .and_then(|rest| rest.split('"').next());
        assert_eq!(section, Some("FeatureFinderCentroided"));

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
}

// ---------------------------------------------------------------------------
// Loading and the error branches
// ---------------------------------------------------------------------------

/// The loading options of `main_` (`FeatureFinderCentroided.cpp:179-192`).
///
/// MS level 1 only and the executed intensity range `[0, f64::MAX)`: the source
/// comment intends `DBL_MIN` but `std::numeric_limits<DPosition<1>>::min()` is
/// zero, so a zero intensity is kept and a negative one dropped. The retained
/// input has neither, and gives 112 spectra with 3084 peaks either way (A3's
/// executed loader oracle). Evidence: tier 3 for the option values, tier 1 for
/// the loaded counts.
#[test]
fn the_load_options_are_the_executed_ones() {
    let options = FeatureFinderCentroided::peak_file_options().unwrap();
    assert_eq!(options.ms_levels(), [1]);
    assert_eq!(options.intensity_range().min, 0.0);
    assert_eq!(options.intensity_range().max, f64::MAX);

    let loaded = FileHandler::load_experiment_with_options(
        ffc1_input(),
        &[FileType::MzMl, FileType::Raw],
        &options,
    )
    .unwrap();
    assert_eq!(loaded.spectra.len(), 112);
    assert_eq!(
        loaded.spectra.iter().map(MSSpectrum::len).sum::<usize>(),
        3084
    );

    // A zero intensity survives the executed range and a negative one does not.
    let dir = Workdir::new();
    let path = dir.file("intensities.mzML");
    let mut experiment = MSExperiment::new();
    let mut spectrum = MSSpectrum {
        ms_level: 1,
        rt: 1.0,
        ..MSSpectrum::default()
    };
    for (mz, intensity) in [(100.0, -1.0f32), (200.0, 0.0), (300.0, 7.0)] {
        spectrum
            .peaks
            .push(openms::kernel::Peak1D { mz, intensity });
    }
    experiment.spectra.push(spectrum);
    FileHandler::store_experiment(&path, &experiment, Some(FileType::MzMl)).unwrap();
    let filtered =
        FileHandler::load_experiment_with_options(&path, &[FileType::MzMl], &options).unwrap();
    assert_eq!(
        filtered.spectra[0]
            .peaks
            .iter()
            .map(|peak| peak.mz)
            .collect::<Vec<_>>(),
        [200.0, 300.0]
    );
}

/// Oracle `FFC_ms2_only`: the MS-level filter leaves no spectrum, so the
/// source throws `FileEmpty` and `TOPPBase` exits 4 with the message inside its
/// `FileEmpty` wording; no output is written.
#[test]
fn an_input_without_ms1_spectra_is_input_file_empty() {
    let dir = Workdir::new();
    let out = dir.file("FFC_ms2_only.tmp.featureXML");
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-in",
            &text(fixture("SimpleSearchEngine_1.mzML")),
            "-out",
            &out,
        ],
    );
    outcome.assert_exit(ExitCode::InputFileEmpty);
    outcome.assert_err_contains(FeatureFinderCentroided::NO_MS1_SPECTRA_MESSAGE);
    assert!(!Path::new(&out).exists());
}

/// Oracle `FFC_profile_noforce` and `FFC_profile_force`: the first spectrum's
/// stored type decides. Without `-force` the source's `IllegalArgument` reaches
/// `TOPPBase`'s catch-all, exit 8, and nothing is written; with `-force` the run
/// continues (the C++ run then finds the FFC_1 features, this port stops in the
/// algorithm).
#[test]
fn profile_data_is_refused_without_force() {
    let dir = Workdir::new();
    let input = dir.put(
        "profile/FeatureFinderCentroided_1_input.mzML",
        &derive_profile(&fs::read(ffc1_input()).unwrap()),
    );
    let out = dir.file("FeatureFinderCentroided_1.tmp.featureXML");
    let arguments = [
        "-test",
        "-ini",
        &text(ffc1_ini()),
        "-in",
        &input,
        "-out",
        &out,
    ];
    let outcome = run_in(&dir, &arguments);
    outcome.assert_exit(ExitCode::UnknownError);
    outcome.assert_err_contains(FeatureFinderCentroided::PROFILE_DATA_MESSAGE);
    assert!(!Path::new(&out).exists());

    let mut forced: Vec<&str> = arguments.to_vec();
    forced.push("-force");
    let outcome = run_in(&dir, &forced);
    assert_reached_the_algorithm(&outcome);
    outcome.assert_exit(ExitCode::IncompatibleInputData);
    outcome.assert_err_contains(ALGORITHM_NOT_PORTED);
    assert!(!Path::new(&out).exists());
}

/// Oracle `FFC_profile_then_spectrum_representation`: a profile term followed
/// by `MS:1000525` resets the stored type to unknown
/// (`MzMLHandler.cpp:1642-1645`), so the check does not fire and the C++ run
/// produces the FFC_1 output.
#[test]
fn a_spectrum_representation_term_after_a_profile_term_is_not_profile() {
    let dir = Workdir::new();
    let input = dir.put(
        "profile_before_1000525/FeatureFinderCentroided_1_input.mzML",
        &derive_profile_then_representation(&fs::read(ffc1_input()).unwrap()),
    );
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-ini",
            &text(ffc1_ini()),
            "-in",
            &input,
            "-out",
            &dir.file("FeatureFinderCentroided_1.tmp.featureXML"),
        ],
    );
    assert_reached_the_algorithm(&outcome);
    outcome.assert_exit(ExitCode::IncompatibleInputData);
    outcome.assert_err_contains(ALGORITHM_NOT_PORTED);
}

/// Oracle `c5_first_profile_only` and `c5_later_profile_only`: only `exp[0]` is
/// examined (`FeatureFinderCentroided.cpp:214`). One profile spectrum at the
/// front is refused with exit 8; 111 profile spectra behind a centroided first
/// one are not (the C++ run writes the FFC_1 output).
#[test]
fn only_the_first_spectrum_decides_the_profile_check() {
    let source = fs::read(ffc1_input()).unwrap();
    let dir = Workdir::new();

    let first = with_profile_terms(&source, |index| index == 0);
    assert_digest(
        &first,
        "446344da1e94460d68e884186a626fdc2f14db4e",
        "first_profile_only",
    );
    let path = dir.put("first/FeatureFinderCentroided_1_input.mzML", &first);
    let out = dir.file("first.featureXML");
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-ini",
            &text(ffc1_ini()),
            "-in",
            &path,
            "-out",
            &out,
        ],
    );
    outcome.assert_exit(ExitCode::UnknownError);
    outcome.assert_err_contains(FeatureFinderCentroided::PROFILE_DATA_MESSAGE);
    assert!(!Path::new(&out).exists());

    let later = with_profile_terms(&source, |index| index > 0);
    assert_digest(
        &later,
        "1fa22f4184f1af9b378c30ca074cf89950a8aba3",
        "later_profile_only",
    );
    let path = dir.put("later/FeatureFinderCentroided_1_input.mzML", &later);
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-ini",
            &text(ffc1_ini()),
            "-in",
            &path,
            "-out",
            &dir.file("later.featureXML"),
        ],
    );
    assert_reached_the_algorithm(&outcome);
    outcome.assert_exit(ExitCode::IncompatibleInputData);
}

/// Oracle `FFC_im_peak_with_units` and `FFC_im_peak_without_units`: a per-peak
/// ion-mobility array on the MS1 spectra is refused with exit 11 and the
/// source's message, whether or not the array carries a unit.
#[test]
fn per_peak_ion_mobility_is_incompatible_input_data() {
    let source = fs::read(ffc1_input()).unwrap();
    for with_units in [true, false] {
        let dir = Workdir::new();
        let input = dir.put("im_peak.mzML", &derive_im_peak(&source, with_units));
        let out = dir.file("FFC_im_peak.tmp.featureXML");
        let outcome = run_in(
            &dir,
            &[
                "-test",
                "-ini",
                &text(ffc1_ini()),
                "-in",
                &input,
                "-out",
                &out,
            ],
        );
        outcome.assert_exit(ExitCode::IncompatibleInputData);
        outcome.assert_err_contains(&FeatureFinderCentroided::im_peak_message());
        assert!(!Path::new(&out).exists());
    }
}

/// Oracle `c5_im_peak_profile_noforce`: the ion-mobility check runs before the
/// profile check (`FeatureFinderCentroided.cpp:200-222`), so an input that is
/// both exits 11 with the ion-mobility message, not 8.
#[test]
fn the_ion_mobility_check_precedes_the_profile_check() {
    let dir = Workdir::new();
    let source = fs::read(ffc1_input()).unwrap();
    let derived = with_profile_terms(&derive_im_peak(&source, true), |_| true);
    assert_digest(
        &derived,
        "117e6a477f1a18a103e466c4567379cabb6baa78",
        "im_peak_profile",
    );
    let input = dir.put("im_peak_profile.mzML", &derived);
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-ini",
            &text(ffc1_ini()),
            "-in",
            &input,
            "-out",
            &dir.file("c5.tmp.featureXML"),
        ],
    );
    outcome.assert_exit(ExitCode::IncompatibleInputData);
    outcome.assert_err_contains(&FeatureFinderCentroided::im_peak_message());
    assert!(!outcome.err.contains("Profile data provided"));
}

/// Oracle `FFC_im_arrays_ms2_only`: ion-mobility arrays that exist only on MS2
/// spectra are removed by the MS-level filter, so the wrapper does not refuse
/// the input. The C++ run then fails inside the algorithm through a Debug-only
/// precondition (`ProgressLogger::init : invalid range!`), which is why no exit
/// code is asserted here.
#[test]
fn ion_mobility_arrays_on_ms2_spectra_only_are_not_refused() {
    let dir = Workdir::new();
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-in",
            &text(fixture("FileConverter_31_output.mzML")),
            "-out",
            &dir.file("FFC_im_ms2.tmp.featureXML"),
        ],
    );
    assert_reached_the_algorithm(&outcome);
}

/// Oracle `FFC_FileFilter_44_noforce`: a registered upstream profile fixture,
/// numpress-compressed, is refused by the profile check with exit 8. The
/// `-force` companion (`FFC_FileFilter_44_force`) exits 8 through a Debug-only
/// precondition in the C++ build, so only the refusal is an expectation.
#[test]
fn the_numpress_profile_fixture_is_refused_without_force() {
    let dir = Workdir::new();
    let out = dir.file("FFC_FileFilter_44.tmp.featureXML");
    let input = text(fixture("FileFilter_44_input.mzML"));
    let outcome = run_in(&dir, &["-test", "-in", &input, "-out", &out]);
    outcome.assert_exit(ExitCode::UnknownError);
    outcome.assert_err_contains(FeatureFinderCentroided::PROFILE_DATA_MESSAGE);
    assert!(!Path::new(&out).exists());

    let outcome = run_in(&dir, &["-test", "-in", &input, "-out", &out, "-force"]);
    assert_reached_the_algorithm(&outcome);
}

/// Oracle `c5_negative_intensities`: with every MS1 peak negative the load
/// filter empties the spectra, and the algorithm's `getSize() == 0` check
/// throws `IllegalArgument`, which `TOPPBase` reports as an unexpected internal
/// error with exit 8. The port reports the same message and code.
#[test]
fn an_input_whose_peaks_are_all_filtered_is_an_unexpected_internal_error() {
    let dir = Workdir::new();
    let input = dir.put(
        "negative_intensities.mzML",
        &derive_negative_intensities(&fs::read(ffc1_input()).unwrap()),
    );
    let out = dir.file("c5.tmp.featureXML");
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-ini",
            &text(ffc1_ini()),
            "-in",
            &input,
            "-out",
            &out,
        ],
    );
    outcome.assert_exit(ExitCode::UnknownError);
    outcome.assert_err_contains(
        "Error: Unexpected internal error (FeatureFinder needs updated ranges on input map. Aborting.)",
    );
    assert!(!Path::new(&out).exists());
}

/// Oracle `FFC_out_no_extension`: an output name without an extension is
/// accepted, because the format check only refuses an extension another type
/// claims. The C++ run writes the FFC_1 featureXML into it.
#[test]
fn an_output_name_without_an_extension_is_accepted() {
    let dir = Workdir::new();
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-ini",
            &text(ffc1_ini()),
            "-in",
            &text(ffc1_input()),
            "-out",
            &dir.file("FFC_noext"),
        ],
    );
    assert!(
        !outcome.err.contains("Invalid output file extension"),
        "{}",
        outcome.err
    );
    assert_reached_the_algorithm(&outcome);
    outcome.assert_exit(ExitCode::IncompatibleInputData);
}

/// Oracle `FFC_invalid_rt_shape`: an invalid value of a subsection parameter is
/// refused by the strict INI/command-line update with exit 6, naming the value
/// and the valid strings.
#[test]
fn an_invalid_algorithm_value_is_illegal_parameters() {
    let dir = Workdir::new();
    let out = dir.file("FeatureFinderCentroided_1.tmp.featureXML");
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-ini",
            &text(ffc1_ini()),
            "-in",
            &text(ffc1_input()),
            "-out",
            &out,
            "-algorithm:feature:rt_shape",
            "bogus",
        ],
    );
    outcome.assert_exit(ExitCode::IllegalParameters);
    outcome.assert_err_contains(
        "Invalid string parameter value 'bogus' for parameter 'rt_shape' given! Valid values are: 'symmetric,asymmetric'.",
    );
    outcome.assert_err_contains(
        "Parameters passed to 'FeatureFinderCentroided' are invalid. To prevent usage of wrong defaults, please update/fix the parameters!",
    );
    assert!(!Path::new(&out).exists());
}

/// Oracle `c5_seeds_not_featurexml`: `-seeds` accepts featureXML only, and the
/// framework's input-format check refuses an mzML file with exit 6 before the
/// tool body runs.
#[test]
fn seeds_must_be_featurexml() {
    let dir = Workdir::new();
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-ini",
            &text(ffc1_ini()),
            "-in",
            &text(ffc1_input()),
            "-seeds",
            &text(ffc1_input()),
            "-out",
            &dir.file("c5.tmp.featureXML"),
        ],
    );
    outcome.assert_exit(ExitCode::IllegalParameters);
    outcome.assert_err_contains("has invalid format 'mzML'. Valid formats are: 'featureXML'.");
}

/// A `-seeds` map is loaded, and the run continues into the algorithm. The
/// retained FFC_1 output is the seed list of the oracle case `FFC_seeds`, whose
/// C++ run reported 24 seeds and 8 features; the port stops in the algorithm,
/// so only the load is asserted here.
#[test]
fn a_featurexml_seed_list_is_loaded() {
    let dir = Workdir::new();
    let seeds = mobility("FeatureFinderCentroided_1_1_output.featureXML");
    assert_eq!(
        FileHandler::load_feature_map(&seeds, &[FileType::FeatureXml])
            .unwrap()
            .features
            .len(),
        8
    );
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-no_progress",
            "-ini",
            &text(ffc1_ini()),
            "-seeds",
            &text(seeds),
            "-in",
            &text(ffc1_input()),
            "-out",
            &dir.file("g.featureXML"),
        ],
    );
    assert_reached_the_algorithm(&outcome);
    outcome.assert_exit(ExitCode::IncompatibleInputData);
}

/// Oracle `c5_faims_corrupt_seeds`: the seed list is loaded before the FAIMS
/// split (`FeatureFinderCentroided.cpp:224-241`), so an unreadable seed file on
/// FAIMS input fails as a parse error with exit 3, not as the FAIMS refusal.
#[test]
fn seeds_are_loaded_before_the_faims_check() {
    let dir = Workdir::new();
    let source = fs::read(ffc1_input()).unwrap();
    let input = dir.put("faims_two_cv.mzML", &derive_faims_two_cv(&source));
    let retained = fs::read(mobility("FeatureFinderCentroided_1_1_output.featureXML")).unwrap();
    let seeds = dir.put("corrupt_seeds.featureXML", &retained[..3000]);
    let out = dir.file("c5.tmp.featureXML");
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-ini",
            &text(ffc1_ini()),
            "-in",
            &input,
            "-seeds",
            &seeds,
            "-out",
            &out,
        ],
    );
    outcome.assert_exit(ExitCode::InputFileCorrupt);
    outcome.assert_err_contains("Error: Unable to read file");
    assert!(!outcome.err.contains("FAIMS"), "{}", outcome.err);
    assert!(!Path::new(&out).exists());
}

// ---------------------------------------------------------------------------
// The FAIMS refusal (decision D5)
// ---------------------------------------------------------------------------

/// The refusal message and the informational line, for one FAIMS input.
fn assert_faims_refusal(outcome: &Outcome, voltages: &[f64], out_path: &str) {
    outcome.assert_exit(ExitCode::IncompatibleInputData);
    outcome.assert_out_contains(&FeatureFinderCentroided::faims_detected_message(
        voltages.len(),
    ));
    outcome.assert_err_contains(&FeatureFinderCentroided::faims_refusal_message(voltages));
    assert!(
        !outcome
            .out
            .contains(FeatureFinderCentroided::NO_FAIMS_MESSAGE),
        "{}",
        outcome.out
    );
    assert!(!Path::new(out_path).exists(), "no output is written");
}

/// Oracle `FFC_faims_test_data`, `FFC_faims_interleaved_force`,
/// `FFC_faims_interleaved_force_nomerge`, `FFC_faims_one_cv`,
/// `FFC_faims_two_cv` and `FFC_faims_two_cv_nomerge`: the C++ tool reaches the
/// algorithm and fails with exit 8 (`the value '1' was used but is not valid;
/// No ranges for this MS level`) on every FAIMS input, whatever
/// `-faims_merge_features` says. Decision D5 defers the closure, so this port
/// refuses such input explicitly with exit 11 and writes nothing.
#[test]
fn faims_input_is_refused() {
    let source = fs::read(ffc1_input()).unwrap();
    let dir = Workdir::new();
    let one_cv = dir.put("faims_one_cv.mzML", &derive_faims_one_cv(&source));
    let two_cv = dir.put("faims_two_cv.mzML", &derive_faims_two_cv(&source));
    let interleaved = text(mobility("FAIMS_CV-60C_V-45_Interleaved.mzML"));
    let test_data = text(mobility("FAIMS_test_data.mzML"));

    let ini = text(ffc1_ini());
    for (case, input, voltages, extra) in [
        ("FFC_faims_one_cv", one_cv.as_str(), &[-45.0][..], &[][..]),
        (
            "FFC_faims_two_cv",
            two_cv.as_str(),
            &[-60.0, -45.0][..],
            &[][..],
        ),
        (
            "FFC_faims_two_cv_nomerge",
            two_cv.as_str(),
            &[-60.0, -45.0][..],
            &["-faims_merge_features", "false"][..],
        ),
    ] {
        let out = dir.file(&format!("{case}.featureXML"));
        let mut arguments = vec!["-test", "-ini", &ini, "-in", input, "-out", &out];
        arguments.extend_from_slice(extra);
        let outcome = run_in(&dir, &arguments);
        assert_faims_refusal(&outcome, voltages, &out);
    }

    // The interleaved upstream file is profile data, so without -force the
    // profile check fires first, exactly as in oracle
    // FFC_faims_interleaved_noforce (exit 8).
    let out = dir.file("interleaved.featureXML");
    let outcome = run_in(&dir, &["-test", "-in", &interleaved, "-out", &out]);
    outcome.assert_exit(ExitCode::UnknownError);
    outcome.assert_err_contains(FeatureFinderCentroided::PROFILE_DATA_MESSAGE);
    let outcome = run_in(
        &dir,
        &["-test", "-in", &interleaved, "-out", &out, "-force"],
    );
    assert_faims_refusal(&outcome, &[-60.0, -45.0], &out);

    let out = dir.file("test_data.featureXML");
    let outcome = run_in(&dir, &["-test", "-in", &test_data, "-out", &out]);
    assert_faims_refusal(&outcome, &[-65.0], &out);
}

/// Oracle `c5_faims_partial_cv`: a compensation voltage on half of the MS1
/// spectra is still FAIMS input. The C++ tool assigns the annotated spectra to
/// the voltage group, skips the others with a warning and fails with exit 8;
/// this port refuses the input.
#[test]
fn a_faims_voltage_on_some_spectra_is_refused() {
    let dir = Workdir::new();
    let source = fs::read(ffc1_input()).unwrap();
    let derived = derive_faims(&source, &["-45", ""]);
    assert_digest(
        &derived,
        "bf274d61a77be9b2af64e51ab4f5bf266452aab9",
        "faims_partial_cv",
    );
    let input = dir.put("faims_partial_cv.mzML", &derived);
    let out = dir.file("c5.tmp.featureXML");
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-ini",
            &text(ffc1_ini()),
            "-in",
            &input,
            "-out",
            &out,
        ],
    );
    assert_faims_refusal(&outcome, &[-45.0], &out);
}

// ---------------------------------------------------------------------------
// TOPP_FeatureFinderCentroided_1 itself
// ---------------------------------------------------------------------------

/// The registered upstream workflow (test-data `topp/CMakeLists.txt:425-428`,
/// oracle `TOPP_FeatureFinderCentroided_1`): the C++ tool exits 0 and writes
/// the eight features of `FeatureFinderCentroided_1_1_output.featureXML`.
///
/// Until package B7 ports seed extension and fitting, the run stops in the
/// algorithm with the documented message and exit 11 and writes nothing; the
/// output comparison is package B10's.
#[test]
fn the_upstream_workflow_stops_in_the_unported_algorithm() {
    let dir = Workdir::new();
    let out = dir.file("FeatureFinderCentroided_1.tmp.featureXML");
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-ini",
            &text(ffc1_ini()),
            "-in",
            &text(ffc1_input()),
            "-out",
            &out,
        ],
    );
    assert_reached_the_algorithm(&outcome);
    outcome.assert_exit(ExitCode::IncompatibleInputData);
    outcome.assert_err_contains(ALGORITHM_NOT_PORTED);
    assert!(!Path::new(&out).exists());
}

// ---------------------------------------------------------------------------
// The annotation and clean-up steps
// ---------------------------------------------------------------------------

thread_local! {
    /// The map `FinishProbe` produced, for the assertions below.
    static FINISHED: RefCell<Option<FeatureMap>> = const { RefCell::new(None) };
}

/// A stand-in tool that runs only the wrapper's post-algorithm steps.
///
/// It registers the same parameters under the same name as the real tool, so
/// its [`ToolContext`] is the one `FeatureFinderCentroided::run_io` would build,
/// and its processing record names the same software.
struct FinishProbe;

impl Tool for FinishProbe {
    const NAME: &'static str = "FeatureFinderCentroided";
    const DESCRIPTION: &'static str = "Detects two-dimensional features in LC-MS data.";

    fn register(spec: &mut ToolSpec) -> Result<()> {
        <FeatureFinderCentroided as Tool>::register(spec)
    }

    fn subsection_defaults(section: &str) -> Result<Option<openms::param::Param>> {
        <FeatureFinderCentroided as Tool>::subsection_defaults(section)
    }

    fn run(_ctx: &ToolContext) -> Result<ExitCode> {
        Err(Error::Unsupported("the probe needs its streams".into()))
    }

    fn run_io(ctx: &ToolContext, out: &mut dyn Write, _err: &mut dyn Write) -> Result<ExitCode> {
        let features =
            FeatureFinderCentroided::finish_features(ctx, ctx.string("in")?, sample_map(), out)?;
        FINISHED.with(|cell| *cell.borrow_mut() = Some(features));
        Ok(ExitCode::ExecutionOk)
    }
}

/// A hull with four corners at the given bounding box, in the order
/// `ConvexHull2D::expandToBoundingBox` adds them.
fn box_hull(rt: (f64, f64), mz: (f64, f64)) -> Vec<Point2D> {
    ConvexHull2D::from_points(&[
        Point2D::new(rt.0, mz.0),
        Point2D::new(rt.0, mz.1),
        Point2D::new(rt.1, mz.0),
        Point2D::new(rt.1, mz.1),
    ])
    .unwrap()
    .hull_points()
}

/// Two features shaped like the algorithm's output: mass-trace hulls with more
/// than four points, metadata, and one subordinate feature.
fn sample_map() -> FeatureMap {
    let hull = |rt: f64| {
        ConvexHull2D::from_points(&[
            Point2D::new(rt, 300.0),
            Point2D::new(rt + 1.0, 300.5),
            Point2D::new(rt + 2.0, 301.0),
            Point2D::new(rt + 3.0, 300.25),
        ])
        .unwrap()
    };
    let feature = |rt: f64, mz: f64, hulls: Vec<ConvexHull2D>| Feature {
        base: BaseFeature {
            rt,
            mz,
            intensity: 1000.0,
            charge: 2,
            ..BaseFeature::default()
        },
        convex_hulls: hulls,
        ..Feature::default()
    };
    let mut first = feature(100.0, 300.0, vec![hull(100.0), hull(110.0)]);
    first
        .metadata
        .insert("label".into(), MetaValue::from("feature one"));
    first
        .metadata
        .insert("score_fit".into(), MetaValue::try_from(0.75).unwrap());
    first.subordinates = vec![feature(101.0, 300.5, vec![hull(101.0)])];

    let second = feature(200.0, 400.0, vec![hull(200.0)]);

    let mut map = FeatureMap::from_features(vec![first, second]);
    map.unique_id = 0;
    map
}

/// Run the probe with the given extra arguments and return its map and output.
fn finish(extra: &[&str]) -> (FeatureMap, String) {
    let dir = Workdir::new();
    let mut arguments = vec![
        "FeatureFinderCentroided".to_owned(),
        "-in".to_owned(),
        text(ffc1_input()),
        "-out".to_owned(),
        dir.file("probe.featureXML"),
    ];
    arguments.extend(extra.iter().map(|argument| (*argument).to_owned()));
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<FinishProbe>(&arguments, &mut out, &mut err);
    assert_eq!(
        code,
        ExitCode::ExecutionOk,
        "{}",
        String::from_utf8_lossy(&err)
    );
    let map = FINISHED.with(|cell| cell.borrow_mut().take()).unwrap();
    (map, String::from_utf8(out).unwrap())
}

/// Source steps 10 to 13 under `-test` (`FeatureFinderCentroided.cpp:318-373`):
/// the primary run path is the base name behind `file://`, the ids come from
/// the seeded generator in the source's visiting order, the processing record
/// is the test-mode `Quantitation` one, and below debug level 5 the hulls
/// become bounding boxes and the subordinates go.
///
/// Evidence: tier 3 (source review) for the steps; the unique-id sequence is
/// the seeded generator's, which `tests/topp_cli_lifecycle.rs` pins against the
/// C++ seed.
#[test]
fn the_output_is_annotated_and_cleaned_in_test_mode() {
    let (map, out) = finish(&["-test"]);
    assert!(out.is_empty(), "{out}");

    assert_eq!(
        map.primary_ms_run_path().unwrap(),
        ["file://FeatureFinderCentroided_1_input.mzML"]
    );

    // ensureUniqueId draws one id for the map, which applyMemberFunction then
    // overwrites; after that the map, both features and the subordinate are
    // drawn depth first.
    let mut generator = UniqueIdGenerator::from_seed(TEST_MODE_UNIQUE_ID_SEED);
    let draws: Vec<u64> = (0..5).map(|_| generator.get_unique_id()).collect();
    assert_eq!(map.unique_id, draws[1]);
    assert_eq!(map.features[0].unique_id, draws[2]);
    assert_eq!(map.features[1].unique_id, draws[4]);

    assert_eq!(map.data_processing.len(), 1);
    let processing = &map.data_processing[0];
    assert!(processing.actions.contains(&ProcessingAction::Quantitation));
    assert_eq!(processing.software.name, "FeatureFinderCentroided");
    assert_eq!(processing.software.version, TEST_MODE_VERSION);
    assert_eq!(
        processing.completion_time.unwrap().to_string(),
        TEST_MODE_COMPLETION_TIME
    );
    assert_eq!(
        processing.metadata.get(TEST_MODE_PARAMETER_KEY),
        Some(&MetaValue::from(TEST_MODE_PARAMETER_VALUE))
    );

    // Debug level 0: every hull is its bounding box and no subordinate remains.
    assert_eq!(map.features[0].convex_hulls.len(), 2);
    assert_eq!(
        map.features[0].convex_hulls[0].hull_points(),
        box_hull((100.0, 103.0), (300.0, 301.0))
    );
    assert_eq!(
        map.features[0].convex_hulls[1].hull_points(),
        box_hull((110.0, 113.0), (300.0, 301.0))
    );
    assert!(map.features[0].subordinates.is_empty());
    assert!(map.features[1].subordinates.is_empty());
}

/// Without `-test` the primary run path is the input path as given, the
/// processing record carries the product version and every resolved parameter,
/// and the ids come from an unseeded generator.
#[test]
fn the_output_records_the_input_path_and_parameters_without_test_mode() {
    let (map, _) = finish(&[]);
    assert_eq!(map.primary_ms_run_path().unwrap(), [text(ffc1_input())]);

    let processing = &map.data_processing[0];
    assert_eq!(processing.software.version, "1.0.0");
    assert_eq!(
        processing.metadata.get("parameter: faims_merge_features"),
        Some(&MetaValue::from("true"))
    );
    assert_eq!(
        processing
            .metadata
            .get("parameter: algorithm:feature:rt_shape"),
        Some(&MetaValue::from("symmetric"))
    );
    assert!(!processing.metadata.contains_key(TEST_MODE_PARAMETER_KEY));
    assert_ne!(map.unique_id, 0);
}

/// From `-debug 5` on, the source keeps the mass-trace hulls and the
/// subordinate features; the oracle case `FFC_debug5` writes 1054 hull points
/// where the default run writes 120.
#[test]
fn debug_level_five_keeps_the_hulls_and_subordinates() {
    let (map, _) = finish(&["-test", "-debug", "5"]);
    assert_eq!(map.features[0].convex_hulls[0].hull_points().len(), 4 + 2);
    assert_eq!(map.features[0].subordinates.len(), 1);
    assert_eq!(
        map.features[0].subordinates[0].convex_hulls[0]
            .hull_points()
            .len(),
        4 + 2
    );
}

/// Above `-debug 10` the source lists the metadata of every feature that has
/// any; its own ordering and number formatting differ (see the item
/// documentation).
#[test]
fn debug_level_eleven_lists_the_feature_metadata() {
    let (map, out) = finish(&["-test", "-debug", "11"]);
    let expected = format!(
        "Feature {}\n  label = feature one\n  score_fit = 0.75\n",
        map.features[0].unique_id
    );
    assert_eq!(out, expected, "{out}");
}

/// The hulls of a feature without any, and a map without features, pass
/// through the clean-up unchanged; the map still receives its annotations.
#[test]
fn an_empty_feature_map_is_annotated() {
    let dir = Workdir::new();
    let arguments = [
        "FeatureFinderCentroided".to_owned(),
        "-in".to_owned(),
        text(ffc1_input()),
        "-out".to_owned(),
        dir.file("probe.featureXML"),
        "-test".to_owned(),
    ];
    struct EmptyProbe;
    impl Tool for EmptyProbe {
        const NAME: &'static str = "FeatureFinderCentroided";
        const DESCRIPTION: &'static str = "Detects two-dimensional features in LC-MS data.";
        fn register(spec: &mut ToolSpec) -> Result<()> {
            <FeatureFinderCentroided as Tool>::register(spec)
        }
        fn subsection_defaults(section: &str) -> Result<Option<openms::param::Param>> {
            <FeatureFinderCentroided as Tool>::subsection_defaults(section)
        }
        fn run(_ctx: &ToolContext) -> Result<ExitCode> {
            Err(Error::Unsupported("the probe needs its streams".into()))
        }
        fn run_io(
            ctx: &ToolContext,
            out: &mut dyn Write,
            _err: &mut dyn Write,
        ) -> Result<ExitCode> {
            let map = FeatureFinderCentroided::finish_features(
                ctx,
                ctx.string("in")?,
                FeatureMap::new(),
                out,
            )?;
            FINISHED.with(|cell| *cell.borrow_mut() = Some(map));
            Ok(ExitCode::ExecutionOk)
        }
    }
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<EmptyProbe>(&arguments, &mut out, &mut err);
    assert_eq!(
        code,
        ExitCode::ExecutionOk,
        "{}",
        String::from_utf8_lossy(&err)
    );
    let map = FINISHED.with(|cell| cell.borrow_mut().take()).unwrap();
    assert!(map.features.is_empty());
    assert_eq!(map.data_processing.len(), 1);
    assert_eq!(
        map.primary_ms_run_path().unwrap(),
        ["file://FeatureFinderCentroided_1_input.mzML"]
    );
}
