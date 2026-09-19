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
//! * **The C++ Release build.** The cases that run the whole chain are also
//!   measured against
//!   `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576/bin/FeatureFinderCentroided`,
//!   which is the build the Debug oracle's featureXML reproduces to within
//!   `9e-11` relative. That run is what the numeric expectations below are
//!   quoted from; the two oracle cases whose Debug exit comes from an
//!   `OPENMS_PRECONDITION` (`FFC_FileFilter_44_force`, `FFC_im_arrays_ms2_only`)
//!   have no Debug exit-code expectation (decision D7), and the Release
//!   behaviour is recorded with them instead.
//!
//! The algorithm now runs end to end (package B7), so an input that passes the
//! wrapper produces a feature map. The exact acceptance of
//! `TOPP_FeatureFinderCentroided_1` — the `1e-9` tight comparison against the
//! C1 oracle output and the `-algorithm:fit:max_iterations` sweep — stays
//! package B10's; what this file asserts is the decoded comparison against the
//! retained expectation plus the structure the C++ Release run showed. What the
//! wrapper does after the algorithm is additionally tested on its own through
//! `FeatureFinderCentroided::finish_features`.
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

#[path = "support/decoded_compare.rs"]
mod decoded;
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

/// Assert that the wrapper accepted the input and handed it to the algorithm.
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

/// The algorithm's progress lines of the `TOPP_FeatureFinderCentroided_1` run,
/// in the order the C++ Release build printed them (its `stdout` also carries
/// the framework's loading progress and timing lines, which are not compared).
const FFC1_ALGORITHM_LINES: &[&str] = &[
    "Not FAIMS compensation voltages found in the data. Returning PeakMap as CV NaN.",
    "Found 25 seeds for charge 2.",
    "Found 8 feature candidates for charge 2.",
    "Removed 0 overlapping features.",
    "",
    "Info: reasons for not finalizing a feature during its construction:",
    " - Invalid fit: Fitted model is bigger than 'max_rt_span': 1 times",
    "",
    "8 features found.",
];

/// The same lines for the `-seeds` run, whose 24 given seeds give one candidate
/// more and one overlap removal (oracle `FFC_seeds`, C++ Release confirmed).
const FFC_SEEDS_ALGORITHM_LINES: &[&str] = &[
    "Not FAIMS compensation voltages found in the data. Returning PeakMap as CV NaN.",
    "Found 24 seeds for charge 2.",
    "Found 9 feature candidates for charge 2.",
    "Removed 1 overlapping features.",
    "",
    "Info: reasons for not finalizing a feature during its construction:",
    " - Could not extend seed: 1 times",
    "",
    "8 features found.",
];

/// Assert that `lines` appear in `outcome.out` consecutively and in order.
/// One whole line of stdout, so that a message is not matched inside a longer
/// one.
fn assert_out_contains_line(outcome: &Outcome, line: &str) {
    assert!(
        outcome.out.lines().any(|written| written == line),
        "stdout has no line {line:?}:\n{}",
        outcome.out
    );
}

fn assert_out_block(outcome: &Outcome, lines: &[&str]) {
    let actual: Vec<&str> = outcome.out.lines().collect();
    let found = actual.windows(lines.len()).any(|window| window == lines);
    assert!(
        found,
        "stdout does not contain\n{lines:#?}\ngot\n{}",
        outcome.out
    );
}

/// The retained upstream expectation of `TOPP_FeatureFinderCentroided_1`
/// (`FeatureFinderCentroided_1_1_output.featureXML`), loaded as a feature map.
fn ffc1_expected_map() -> FeatureMap {
    FileHandler::load_feature_map(
        mobility("FeatureFinderCentroided_1_1_output.featureXML"),
        &[FileType::FeatureXml],
    )
    .unwrap()
}

/// Compare a written featureXML with the retained `TOPP_FeatureFinderCentroided_1`
/// expectation under decision D6: decoded content, the upstream `FuzzyDiff` rule
/// (ratio `1.01` or absolute difference `0.01`) and generated identifiers
/// skipped, which is the decoded form of the upstream `-whitelist "id="`.
///
/// The loose rule is the upstream one and is never relied on alone: every
/// caller also pins the structure, and
/// [`the_upstream_workflow_matches_the_retained_expectation`] pins the numbers
/// far more tightly than `FuzzyDiff` would.
fn assert_matches_ffc1_expectation(path: &str) {
    let actual = FileHandler::load_feature_map(path, &[FileType::FeatureXml]).unwrap();
    let options =
        decoded::DecodedOptions::new(decoded::Tolerance::new(1.01, 0.01)).ignoring_unique_ids();
    if let Err(mismatch) = decoded::compare_feature_maps(&actual, &ffc1_expected_map(), &options) {
        panic!("{path}: {mismatch}");
    }
}

/// The structure every `TOPP_FeatureFinderCentroided_1`-shaped run wrote, read
/// off the C++ Release output: eight features of charge 2, each with four-point
/// mass-trace hulls and no subordinates, 30 hulls and 120 hull points in all,
/// and one `UserParam` key set.
fn assert_ffc1_structure(path: &str, input_basename: &str) {
    let map = FileHandler::load_feature_map(path, &[FileType::FeatureXml]).unwrap();
    assert_eq!(map.features.len(), 8, "{path}: feature count");
    assert_eq!(
        map.primary_ms_run_path().unwrap(),
        [format!("file://{input_basename}")],
        "{path}: spectra_data"
    );
    let mut hulls = 0;
    let mut points = 0;
    for (index, feature) in map.features.iter().enumerate() {
        assert_eq!(feature.base.charge, 2, "{path}: features[{index}].charge");
        assert!(
            feature.subordinates.is_empty(),
            "{path}: features[{index}] has subordinates"
        );
        hulls += feature.convex_hulls.len();
        for hull in &feature.convex_hulls {
            assert_eq!(
                hull.hull_points().len(),
                4,
                "{path}: features[{index}] hull is not a bounding box"
            );
            points += hull.hull_points().len();
        }
        // `MetaInfo` is a `BTreeMap`, so its keys already come out in order.
        let keys: Vec<&str> = feature.metadata.keys().map(String::as_str).collect();
        assert_eq!(
            keys,
            [
                "FWHM",
                "label",
                "num_of_datapoints",
                "score_correlation",
                "score_fit",
                "spectrum_index",
                "spectrum_native_id",
            ],
            "{path}: features[{index}] metadata keys"
        );
    }
    assert_eq!(hulls, 30, "{path}: hull count");
    assert_eq!(points, 120, "{path}: hull point count");

    assert_eq!(map.data_processing.len(), 1, "{path}: processing entries");
    let processing = &map.data_processing[0];
    assert_eq!(processing.software.name, "FeatureFinderCentroided");
    assert_eq!(processing.software.version, TEST_MODE_VERSION);
    assert_eq!(
        processing.completion_time.unwrap().to_string(),
        TEST_MODE_COMPLETION_TIME
    );
    assert_eq!(processing.actions.len(), 1);
    assert!(processing.actions.contains(&ProcessingAction::Quantitation));
    assert_eq!(
        processing.metadata.get(TEST_MODE_PARAMETER_KEY),
        Some(&MetaValue::from(TEST_MODE_PARAMETER_VALUE))
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
/// continues and finds the `TOPP_FeatureFinderCentroided_1` features — the
/// oracle's `-force` output is byte-identical to its `FFC_1` output, and the
/// C++ Release run reproduces that.
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
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert_out_block(&outcome, FFC1_ALGORITHM_LINES);
    assert_ffc1_structure(&out, "FeatureFinderCentroided_1_input.mzML");
    assert_matches_ffc1_expectation(&out);
}

/// Oracle `FFC_profile_then_spectrum_representation`: a profile term followed
/// by `MS:1000525` resets the stored type to unknown
/// (`MzMLHandler.cpp:1642-1645`), so the check does not fire and the run
/// produces the FFC_1 output — the oracle wrote a file identical to its FFC_1
/// output, and the C++ Release run reproduces that.
#[test]
fn a_spectrum_representation_term_after_a_profile_term_is_not_profile() {
    let dir = Workdir::new();
    let input = dir.put(
        "profile_before_1000525/FeatureFinderCentroided_1_input.mzML",
        &derive_profile_then_representation(&fs::read(ffc1_input()).unwrap()),
    );
    let out = dir.file("FeatureFinderCentroided_1.tmp.featureXML");
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
    assert_reached_the_algorithm(&outcome);
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert_out_block(&outcome, FFC1_ALGORITHM_LINES);
    assert_ffc1_structure(&out, "FeatureFinderCentroided_1_input.mzML");
    assert_matches_ffc1_expectation(&out);
}

/// Oracle `c5_first_profile_only` and `c5_later_profile_only`: only `exp[0]` is
/// examined (`FeatureFinderCentroided.cpp:214`). One profile spectrum at the
/// front is refused with exit 8; 111 profile spectra behind a centroided first
/// one are not, and that run wrote the FFC_1 output (oracle
/// `c5_later_profile_only`, exit 0, 15234 bytes).
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
    let out = dir.file("later.featureXML");
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
    assert_reached_the_algorithm(&outcome);
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert_out_block(&outcome, FFC1_ALGORITHM_LINES);
    assert_ffc1_structure(&out, "FeatureFinderCentroided_1_input.mzML");
    assert_matches_ffc1_expectation(&out);
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
/// the input. The Debug oracle then fails inside the algorithm through an
/// `OPENMS_PRECONDITION` (`ProgressLogger::init : invalid range!`), so its exit
/// 8 is `debug_only` (decision D7) and is not an expectation.
///
/// The C++ **Release** build, which has no preconditions, runs the input to the
/// end: it reports no seed and no candidate for charges 1 to 4, `0 features
/// found.`, exits 0 and writes an empty feature map. That is what this port
/// does, and what is asserted here.
#[test]
fn ion_mobility_arrays_on_ms2_spectra_only_are_not_refused() {
    let dir = Workdir::new();
    let out = dir.file("FFC_im_ms2.tmp.featureXML");
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-in",
            &text(fixture("FileConverter_31_output.mzML")),
            "-out",
            &out,
        ],
    );
    assert_reached_the_algorithm(&outcome);
    outcome.assert_exit(ExitCode::ExecutionOk);
    for charge in 1..=4 {
        outcome.assert_out_contains(&format!("Found 0 seeds for charge {charge}."));
        outcome.assert_out_contains(&format!("Found 0 feature candidates for charge {charge}."));
    }
    outcome.assert_out_contains("0 features found.");
    let map = FileHandler::load_feature_map(&out, &[FileType::FeatureXml]).unwrap();
    assert!(map.features.is_empty());
    assert_eq!(
        map.primary_ms_run_path().unwrap(),
        ["file://FileConverter_31_output.mzML"]
    );
}

/// Oracle `FFC_FileFilter_44_noforce`: a registered upstream profile fixture,
/// numpress-compressed, is refused by the profile check with exit 8. The C++
/// Release build agrees (exit 8, the same message; `../oracle/ffap-sem-completion`
/// case `filefilter_44_noforce`, two repetitions).
///
/// The `-force` companion (`FFC_FileFilter_44_force`) exits 8 in the Debug
/// oracle through `OPENMS_PRECONDITION(begin <= end)` in
/// `ProgressLogger::startProgress` (`ProgressLogger.cpp:235`), called with
/// `(5, 0)` for the mass-trace scores (`FeatureFinderAlgorithmPicked.cpp:298`)
/// because the input is shorter than the seed loop; that exit is `debug_only`
/// (decision D7). What `-force` does in the Release build, and here, is
/// [`a_short_input_never_reaches_the_seed_loop_as_in_the_cpp_release_build`].
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
    outcome.assert_exit(ExitCode::ExecutionOk);
}

/// The algorithm lines of a run in which no seed of any of the default
/// charges 1 to 4 is found, as the C++ Release build printed them.
const NOTHING_FOUND_FOR_CHARGES_1_TO_4: &[&str] = &[
    "Not FAIMS compensation voltages found in the data. Returning PeakMap as CV NaN.",
    "Found 0 seeds for charge 1.",
    "Found 0 feature candidates for charge 1.",
    "Found 0 seeds for charge 2.",
    "Found 0 feature candidates for charge 2.",
    "Found 0 seeds for charge 3.",
    "Found 0 feature candidates for charge 3.",
    "Found 0 seeds for charge 4.",
    "Found 0 feature candidates for charge 4.",
    "Removed 0 overlapping features.",
    "",
    "Info: reasons for not finalizing a feature during its construction:",
    "",
    "0 features found.",
];

/// The same for the FeatureFinderCentroided_1 INI, which searches charge 2
/// only.
const NOTHING_FOUND_FOR_CHARGE_2: &[&str] = &[
    "Not FAIMS compensation voltages found in the data. Returning PeakMap as CV NaN.",
    "Found 0 seeds for charge 2.",
    "Found 0 feature candidates for charge 2.",
    "Removed 0 overlapping features.",
    "",
    "Info: reasons for not finalizing a feature during its construction:",
    "",
    "0 features found.",
];

/// Assert an exit 0 with an empty feature map whose `spectra_data` names
/// `input_basename`, as each of the C++ Release runs quoted below wrote it.
fn assert_empty_release_map(outcome: &Outcome, out: &str, input_basename: &str) {
    outcome.assert_exit(ExitCode::ExecutionOk);
    let map = FileHandler::load_feature_map(out, &[FileType::FeatureXml]).unwrap();
    assert!(map.features.is_empty());
    assert_eq!(
        map.primary_ms_run_path().unwrap(),
        [format!("file://{input_basename}")]
    );
}

/// **The short input: the seed loop is empty.** Measured against the C++
/// Release build (`../oracle/ffap-sem-completion` case `filefilter_44_force`,
/// three repetitions, identical).
///
/// `FileFilter_44_input.mzML` has two MS1 spectra, both at retention time
/// `0.273` s. Without an INI, `mass_trace:min_spectra` is 10, so the source's
/// `min_spectra_` is 5 and its seed loop runs over the scans `5 .. n - min(5,
/// n)` (`FeatureFinderAlgorithmPicked.cpp:493-498`), which is empty for every
/// input of at most 10 scans. No seed can be found whatever the scores are:
/// the C++ Release build prints no seed and no candidate for charges 1 to 4,
/// exits 0 and writes an empty feature map, and so does this port. The zero
/// retention-time extent is incidental here: the executed intensity scores of
/// these peaks are NaN, but nothing reads them
/// (`degenerate_bin_steps_match_the_linux_release_build`,
/// `tests/feature_finder_picked_seeds.rs`), and
/// `FileConverter_31_output.mzML`, four scans at retention times 5 to 8 s
/// ([`ion_mobility_arrays_on_ms2_spectra_only_are_not_refused`]), gives the
/// same result with a non-zero extent. Until this port reproduced the Release
/// build's scoring of a zero extent it refused this input, and this test was
/// ignored under a name that blamed the extent.
#[test]
fn a_short_input_never_reaches_the_seed_loop_as_in_the_cpp_release_build() {
    let dir = Workdir::new();
    let out = dir.file("FFC_FileFilter_44_force.tmp.featureXML");
    let input = text(fixture("FileFilter_44_input.mzML"));
    let outcome = run_in(&dir, &["-test", "-in", &input, "-out", &out, "-force"]);
    outcome.assert_exit(ExitCode::ExecutionOk);
    outcome.assert_out_contains("0 features found.");
    let map = FileHandler::load_feature_map(&out, &[FileType::FeatureXml]).unwrap();
    assert!(map.features.is_empty());
    assert_out_block(&outcome, NOTHING_FOUND_FOR_CHARGES_1_TO_4);
    assert_empty_release_map(&outcome, &out, "FileFilter_44_input.mzML");
}

/// FeatureFinderCentroided_1's input with every `scan start time` replaced by
/// the first one, `4114.53`: a zero retention-time extent (oracle input
/// `zero_rt_ffc1.mzML` of `../oracle/ffap-sem-completion/make_inputs.py`; its
/// sha256 and those of the two m/z inputs are in
/// `tests/data/feature_finder_picked_provenance.json`, `oracle.release_build`).
fn derive_zero_rt(source: &[u8]) -> Vec<u8> {
    const MARKER: &[u8] = br#"name="scan start time" value=""#;
    let mut out = Vec::with_capacity(source.len());
    let mut rest = source;
    let mut replaced = 0;
    while let Some(at) = find(rest, MARKER) {
        let value_start = at + MARKER.len();
        out.extend_from_slice(&rest[..value_start]);
        rest = &rest[value_start..];
        let end = rest.iter().position(|byte| *byte == b'"').unwrap();
        if replaced == 0 {
            assert_eq!(&rest[..end], b"4114.53");
        }
        out.extend_from_slice(b"4114.53");
        rest = &rest[end..];
        replaced += 1;
    }
    out.extend_from_slice(rest);
    assert_eq!(replaced, SPECTRA);
    assert_digest(&out, "c58ee032b11af4ae864fb3016047bf248cbe42ea", "zero_rt");
    out
}

/// FeatureFinderCentroided_1's input with every m/z array payload replaced by
/// `defaultArrayLength` little-endian doubles of 500.0 (oracle input
/// `zero_mz_ffc1.mzML`), and with `first` as the first m/z of the first
/// spectrum (`zero_mz_control_ffc1.mzML`, 499.0).
fn derive_mz(source: &[u8], first: Option<f64>) -> Vec<u8> {
    let engine = base64::engine::general_purpose::STANDARD;
    let mut out: Vec<Vec<u8>> = Vec::new();
    let mut length = 0usize;
    let mut spectrum = None::<usize>;
    let mut pending = false;
    let mut arrays = 0usize;
    for line in lines(source) {
        if let Some(value) = default_array_length(line) {
            length = value;
            spectrum = Some(spectrum.map_or(0, |index| index + 1));
        }
        if find(line, br#"name="m/z array""#).is_some() {
            pending = true;
        }
        if pending && trimmed(line).starts_with(b"<binary>") {
            let mut values = vec![500.0f64; length];
            if spectrum == Some(0) {
                if let Some(value) = first {
                    values[0] = value;
                }
            }
            let raw: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
            let payload = engine.encode(&raw);
            let prefix = &line[..find(line, b"<binary>").unwrap()];
            let old = &trimmed(line)[b"<binary>".len()..trimmed(line).len() - b"</binary>".len()];
            assert_eq!(old.len(), payload.len());
            let mut replaced = prefix.to_vec();
            replaced.extend_from_slice(b"<binary>");
            replaced.extend_from_slice(payload.as_bytes());
            replaced.extend_from_slice(b"</binary>");
            out.push(replaced);
            pending = false;
            arrays += 1;
            continue;
        }
        out.push(line.to_vec());
    }
    assert_eq!(arrays, SPECTRA);
    let derived = join(out);
    let (expected, case) = match first {
        None => ("8311af59499e10154cafab2ff56acbbcde3d9f60", "zero_mz"),
        Some(_) => (
            "5df2a4e1a7fbbf47a2bb7d6ea142bfab61c20d5e",
            "zero_mz_control",
        ),
    };
    assert_digest(&derived, expected, case);
    derived
}

/// FeatureFinderCentroided_1's input with every `scan start time` value `v`
/// written as `v` followed by `suffix` (`e36` or `e39`): the retention times
/// scaled by that power of ten in the text (oracle inputs `rt_e36.mzML` and
/// `rt_e39.mzML` of `../oracle/ffap-complete-fix4/node/run_tool.sh`).
fn derive_rt_suffix(source: &[u8], suffix: &str) -> Vec<u8> {
    const MARKER: &[u8] = br#"name="scan start time" value=""#;
    let mut out = Vec::with_capacity(source.len() + SPECTRA * suffix.len());
    let mut rest = source;
    let mut replaced = 0;
    while let Some(at) = find(rest, MARKER) {
        let value_start = at + MARKER.len();
        let end = value_start
            + rest[value_start..]
                .iter()
                .position(|byte| *byte == b'"')
                .unwrap();
        out.extend_from_slice(&rest[..end]);
        out.extend_from_slice(suffix.as_bytes());
        rest = &rest[end..];
        replaced += 1;
    }
    out.extend_from_slice(rest);
    assert_eq!(replaced, SPECTRA);
    let expected = match suffix {
        "e36" => "baadf62aa5ef060d07c8cf010b768da00adf376c",
        "e39" => "f2acc3a566d722abadf685851c953dd8229c18e7",
        other => panic!("{other}"),
    };
    assert_digest(&out, expected, suffix);
    out
}

/// **Features of infinite width and intensity are written**, measured against
/// the C++ Release build (`../oracle/featurexml-inf`, cases `rt_e36`, `rt_e39`
/// and `rt_e39_egh`, two runs each whose output featureXML is byte-identical;
/// `rt_scaled_release.tsv` and `rt_scaled_release_features.tsv`).
///
/// With every retention time of FeatureFinderCentroided_1 scaled by `1e36` the
/// fitted features' `float` intensities overflow, and with `1e39` their widths
/// too; the Release tool prints its usual lines, exits 0 and writes
/// `<intensity>inf</intensity>` and `FWHM` `inf` into the featureXML, because
/// `precisionWrapper` and `writeUserParam_` spell a non-finite value rather
/// than refusing it. This port now does the same, and the written document
/// carries every value the Release build's does: the fixture holds all 1263
/// of them, taken from the executed output, and each is compared after
/// rendering this port's value the way the source renders it — six fractional
/// digits for a `float` field, fifteen for a `double`
/// (`writtenDigits`, `NumericFormatting.h:27-135`), which is the only
/// difference left, and a difference of text rather than of value.
///
/// This closes TOPP native difference 16, which recorded the earlier refusal.
#[test]
fn infinite_feature_values_are_written_as_the_release_build_writes_them() {
    let source = fs::read(ffc1_input()).unwrap();
    let release = fs::read_to_string(fixture("rt_scaled_release.tsv")).unwrap();
    let values = fs::read_to_string(fixture("rt_scaled_release_features.tsv")).unwrap();
    let ini = text(ffc1_ini());
    for (case, suffix, extra) in [
        ("rt_e36", "e36", &[][..]),
        ("rt_e39", "e39", &[][..]),
        (
            "rt_e39_egh",
            "e39",
            &["-algorithm:feature:rt_shape", "asymmetric"][..],
        ),
    ] {
        let rows: Vec<Vec<&str>> = release
            .lines()
            .skip(1)
            .map(|line| line.split('\t').collect())
            .filter(|row: &Vec<&str>| row[0] == case)
            .collect();
        assert_eq!(rows[0][1], "0", "{case}: the Release tool exits 0");
        let count: usize = rows[0][2].parse().unwrap();
        let intensities: Vec<&str> = rows
            .iter()
            .filter(|row| row[3] == "intensity")
            .map(|row| row[4])
            .collect();
        assert_eq!(intensities.len(), count, "{case}");
        assert!(intensities.iter().all(|value| *value == "inf"), "{case}");
        let widths: Vec<&str> = rows
            .iter()
            .filter(|row| row[3] == "FWHM")
            .map(|row| row[4])
            .collect();
        assert_eq!(widths.len(), count, "{case}");
        assert_eq!(
            widths.iter().all(|value| *value == "inf"),
            suffix == "e39",
            "{case}: {widths:?}"
        );
        let block: Vec<&str> = rows
            .iter()
            .filter(|row| row[3] == "stdout")
            .map(|row| row[4])
            .collect();

        let dir = Workdir::new();
        let input = dir.put(
            &format!("rt_{suffix}.mzML"),
            &derive_rt_suffix(&source, suffix),
        );
        let out = dir.file(&format!("{case}.tmp.featureXML"));
        let mut args = vec!["-test", "-ini", &ini, "-in", &input, "-out", &out];
        args.extend_from_slice(extra);
        let outcome = run_in(&dir, &args);
        assert_out_block(&outcome, &block);
        outcome.assert_exit(ExitCode::ExecutionOk);
        assert_eq!(outcome.err, "", "{case}");

        // The document on disk still spells the infinities the way the
        // Release build spells them.
        let written = fs::read_to_string(&out).unwrap();
        assert!(written.contains("<intensity>inf</intensity>"), "{case}");
        assert_eq!(
            written.matches("<intensity>inf</intensity>").count(),
            count,
            "{case}"
        );
        assert_eq!(
            written
                .matches("name=\"FWHM\" type=\"float\" value=\"inf\"")
                .count(),
            if suffix == "e39" { count } else { 0 },
            "{case}"
        );

        let map = openms::format::featurexml::load(&out).unwrap();
        assert_eq!(map.len(), count, "{case}");
        assert_release_values(&values, case, &map);
    }
}

/// Every recorded value of the Release build's own output document, compared
/// with `map`.
///
/// The fixture holds the Release text; this port's value is rendered the way
/// `NumericFormatting::appendNumeric` renders it before the comparison, so a
/// row passes only when the two agree at the source's own precision. A
/// `float` field is rendered by the `float` instantiation and everything else
/// by the `double` one, exactly as `precisionWrapper` picks them.
fn assert_release_values(fixture: &str, case: &str, map: &FeatureMap) {
    use openms::format::sv_out_stream::{
        F32_FIXED_DIGITS, F64_FIXED_DIGITS, source_f32_text, source_float_text,
    };
    let wide = |value: f64| source_float_text(value, F64_FIXED_DIGITS);
    let narrow = |value: f32| source_f32_text(value, F32_FIXED_DIGITS);
    let mut checked = 0usize;
    for line in fixture.lines().skip(1) {
        let row: Vec<&str> = line.split('\t').collect();
        if row[0] != case {
            continue;
        }
        let (feature, detail, kind, want) = (row[1], row[2], row[3], row[4]);
        if feature == "map" {
            assert_eq!(kind, "count");
            assert_eq!(map.len().to_string(), want, "{case} count");
            checked += 1;
            continue;
        }
        let f = &map.features[feature.parse::<usize>().unwrap()];
        let got = match kind {
            "id" => format!("f_{}", f.unique_id),
            "position0" => wide(f.rt),
            "position1" => wide(f.mz),
            "intensity" => narrow(f.intensity),
            "quality0" => narrow(f.quality_rt),
            "quality1" => narrow(f.quality_mz),
            "overallquality" => narrow(f.quality),
            "charge" => f.charge.to_string(),
            "hull_points" => f.convex_hulls[detail.parse::<usize>().unwrap()]
                .hull_points()
                .len()
                .to_string(),
            other if other.starts_with("pt") => {
                let hull = f.convex_hulls[detail.parse::<usize>().unwrap()].hull_points();
                let (index, axis) = other[2..].split_once('_').unwrap();
                let point = hull[index.parse::<usize>().unwrap()];
                wide(if axis == "x" { point.rt } else { point.mz })
            }
            other => {
                let value = &f.metadata[other.strip_prefix("meta:").unwrap()];
                match detail {
                    "float" => wide(value.as_f64().unwrap()),
                    "int" => value.as_i64().unwrap().to_string(),
                    _ => value.as_str().unwrap().to_owned(),
                }
            }
        };
        assert_eq!(got, want, "{case} feature {feature} {kind}");
        checked += 1;
    }
    assert!(checked > 400, "{case}: only {checked} values compared");
}

/// **A store that fails is a write failure, not a read failure.**
///
/// The source's run-phase catch has one write-side arm: `UnableToCreateFile`
/// becomes `Error: Unable to write file (<what>)` with
/// `CANNOT_WRITE_OUTPUT_FILE` (`TOPPBase.cpp:430-435`, cli pin `c19e494`), and
/// a featureXML store raises exactly that exception when it cannot produce the
/// file — for a stream it cannot open (`XMLFile.cpp:366-372`) and for a name
/// whose extension it does not accept (`FeatureXMLFile.cpp:74-77`), both
/// recorded in `../oracle/featurexml-inf/results/store_failures.tsv`. This
/// port used to report such a failure through the `ParseError` arm instead.
///
/// The case below is the one that reaches the store rather than
/// `outputFileWritable_`: `-out` is an existing **directory** whose name
/// carries the featureXML extension, so the writability check passes.
/// Executed against the Release build on `ibminode06`
/// (`../oracle/featurexml-inf/results/err_dir.txt`, `out_directory.status`):
/// exit 5 and the single line below. With `-out` inside a directory that does
/// not exist the check fires first and both builds print two lines and exit 5
/// (`err_missing.txt`), which is what the port already did.
#[test]
fn a_store_that_fails_is_the_sources_write_failure() {
    let dir = Workdir::new();
    let out = dir.file("adir.featureXML");
    fs::create_dir(&out).unwrap();
    let ini = text(ffc1_ini());
    let input = text(ffc1_input());
    let outcome = run_in(&dir, &["-test", "-ini", &ini, "-in", &input, "-out", &out]);
    outcome.assert_exit(ExitCode::CannotWriteOutputFile);
    assert_eq!(
        outcome.err.trim_end(),
        format!("Error: Unable to write file (the file '{out}' could not be created. )"),
        "stdout:\n{}",
        outcome.out
    );

    let missing = dir.file("nosuch/out.featureXML");
    let outcome = run_in(
        &dir,
        &["-test", "-ini", &ini, "-in", &input, "-out", &missing],
    );
    outcome.assert_exit(ExitCode::CannotWriteOutputFile);
    assert_eq!(
        outcome.err.lines().collect::<Vec<_>>(),
        vec![
            "Cannot write output file given from parameter '-out'!",
            &format!("Error: Unable to write file (the file '{missing}' could not be created. )"),
        ]
    );
}

/// **A zero retention-time extent that reaches the seed loop**, measured
/// against the C++ Release build (`../oracle/ffap-sem-completion` cases
/// `zero_rt`, `zero_rt_threads4` and `zero_rt_min_score_0`, three repetitions
/// each, identical).
///
/// Step 1 divides the zero extent by `intensity:bins` and `intensityScore_`
/// converts `floor(0 / 0)` to `UInt` for every peak
/// (`FeatureFinderAlgorithmPicked.cpp:1837`), which is undefined behaviour.
/// The Release build computes it as `cvttsd2si` does, every intensity score is
/// NaN, and no peak becomes a seed, even with `seed:min_score` 0: the run prints
/// no seed, exits 0 and writes an empty map, at one and at four threads. This
/// port reproduces that (`DegenerateBinStep::Source`, the tool's setting).
#[test]
fn a_zero_width_retention_time_range_follows_the_cpp_release_build() {
    let source = fs::read(ffc1_input()).unwrap();
    let derived = derive_zero_rt(&source);
    for extra in [
        &[][..],
        &["-threads", "4"][..],
        &["-algorithm:seed:min_score", "0"][..],
    ] {
        let dir = Workdir::new();
        let input = dir.put("zero_rt_ffc1.mzML", &derived);
        let out = dir.file("zero_rt.tmp.featureXML");
        let ini = text(ffc1_ini());
        let mut args = vec!["-test", "-ini", &ini, "-in", &input, "-out", &out];
        args.extend_from_slice(extra);
        let outcome = run_in(&dir, &args);
        assert_out_block(&outcome, NOTHING_FOUND_FOR_CHARGE_2);
        assert_empty_release_map(&outcome, &out, "zero_rt_ffc1.mzML");
    }
}

/// **A zero m/z extent that reaches the seed loop**, measured against the C++
/// Release build (cases `zero_mz`, `zero_mz_threads4` and
/// `zero_mz_min_score_0`, three repetitions each, identical): the same
/// undefined conversion at `FeatureFinderAlgorithmPicked.cpp:1838`, the same
/// empty result.
///
/// The control moves one peak to m/z 499 (case
/// `zero_mz_control_min_score_0`, two repetitions): the extent is no longer
/// zero, the intensity scores are defined, and with `seed:min_score` 0 the
/// Release build finds 735 seeds, none of which has an isotope pattern.
#[test]
fn a_zero_width_mz_range_follows_the_cpp_release_build() {
    let source = fs::read(ffc1_input()).unwrap();
    let derived = derive_mz(&source, None);
    let ini = text(ffc1_ini());
    for extra in [
        &[][..],
        &["-threads", "4"][..],
        &["-algorithm:seed:min_score", "0"][..],
    ] {
        let dir = Workdir::new();
        let input = dir.put("zero_mz_ffc1.mzML", &derived);
        let out = dir.file("zero_mz.tmp.featureXML");
        let mut args = vec!["-test", "-ini", &ini, "-in", &input, "-out", &out];
        args.extend_from_slice(extra);
        let outcome = run_in(&dir, &args);
        assert_out_block(&outcome, NOTHING_FOUND_FOR_CHARGE_2);
        assert_empty_release_map(&outcome, &out, "zero_mz_ffc1.mzML");
    }

    let dir = Workdir::new();
    let input = dir.put(
        "zero_mz_control_ffc1.mzML",
        &derive_mz(&source, Some(499.0)),
    );
    let out = dir.file("zero_mz_control.tmp.featureXML");
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-ini",
            &ini,
            "-in",
            &input,
            "-out",
            &out,
            "-algorithm:seed:min_score",
            "0",
        ],
    );
    assert_out_block(
        &outcome,
        &[
            "Not FAIMS compensation voltages found in the data. Returning PeakMap as CV NaN.",
            "Found 735 seeds for charge 2.",
            "Found 0 feature candidates for charge 2.",
            "Removed 0 overlapping features.",
            "",
            "Info: reasons for not finalizing a feature during its construction:",
            " - Could not find good enough isotope pattern containing the seed: 735 times",
            "",
            "0 features found.",
        ],
    );
    assert_empty_release_map(&outcome, &out, "zero_mz_control_ffc1.mzML");
}

/// **A changed isotope abundance: the one designed difference** (`CPP-247`).
///
/// The executed C++ Release tool (`../oracle/ffap-sem-completion` case
/// `ffc1_abundance_12C_90`, three repetitions, identical) builds the override
/// with a stray `(0, 1)` peak and finds nothing on FeatureFinderCentroided_1
/// with `-algorithm:isotopic_pattern:abundance_12C 90`: no seed, no candidate,
/// no feature, exit 0. The tool uses the library default
/// `AbundanceOverride::Intended` (`Options::default()`), which computes the
/// override the source intends and does find features. The expected values of
/// that result are not this port's: they are the adapted Release replay of the
/// intended override (driver `intended_abundance.cpp`, configuration
/// `ffc1_12C_90`, recorded in
/// `tests/data/feature_finder_picked/intended_abundance.tsv`): 18 seeds, one
/// candidate, one feature of charge 2 at RT `0x40b12596f2a04222` and m/z
/// `0x4084420ded67bc0c`, with intensity `81362.57` and quality `0.7551466` as
/// `f32`, 65 data points, and three abort reasons.
#[test]
fn a_changed_abundance_finds_the_intended_features_where_the_cpp_release_build_finds_none() {
    // The executed C++ tool's lines, for the record: this port differs here by
    // design.
    const CPP_RELEASE: &[&str] = &[
        "Not FAIMS compensation voltages found in the data. Returning PeakMap as CV NaN.",
        "Found 0 seeds for charge 2.",
        "Found 0 feature candidates for charge 2.",
        "Removed 0 overlapping features.",
        "",
        "Info: reasons for not finalizing a feature during its construction:",
        "",
        "0 features found.",
    ];
    let dir = Workdir::new();
    let out = dir.file("abundance_12C_90.tmp.featureXML");
    let ini = text(ffc1_ini());
    let input = text(ffc1_input());
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-ini",
            &ini,
            "-in",
            &input,
            "-out",
            &out,
            "-algorithm:isotopic_pattern:abundance_12C",
            "90",
        ],
    );
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert!(
        !outcome.out.contains(CPP_RELEASE[1]),
        "the tool follows the executed stray-peak override:\n{}",
        outcome.out
    );
    assert_out_block(
        &outcome,
        &[
            "Not FAIMS compensation voltages found in the data. Returning PeakMap as CV NaN.",
            "Found 18 seeds for charge 2.",
            "Found 1 feature candidates for charge 2.",
            "Removed 0 overlapping features.",
            "",
            "Info: reasons for not finalizing a feature during its construction:",
            " - Could not extend seed: 2 times",
            " - Could not find good enough isotope pattern containing the seed: 8 times",
            " - Feature quality too low after fit: 6 times",
            "",
            "1 features found.",
        ],
    );
    let map = FileHandler::load_feature_map(&out, &[FileType::FeatureXml]).unwrap();
    assert_eq!(map.features.len(), 1);
    let feature = &map.features[0];
    assert_eq!(feature.charge, 2);
    // The retention time is the fitted centre. It is bit-identical on every
    // platform: the Gaussian fit calls the reference build's glibc `exp` and
    // `log`, ported (lead decision D10; `tests/feature_finder_picked.rs`,
    // `tolerance`).
    let rt = f64::from_bits(0x40b1_2596_f2a0_4222);
    assert_eq!(feature.rt.to_bits(), rt.to_bits());
    assert_eq!(feature.mz.to_bits(), 0x4084_420d_ed67_bc0c);
    assert_eq!(feature.intensity.to_bits(), 0x479e_e949);
    assert_eq!(feature.quality.to_bits(), 0x3f41_5149);
    assert_eq!(
        feature
            .metadata
            .get("num_of_datapoints")
            .map(MetaValue::to_string),
        Some("65".to_owned())
    );
    assert_eq!(
        map.primary_ms_run_path().unwrap(),
        ["file://FeatureFinderCentroided_1_input.mzML"]
    );
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
/// claims, and the FFC_1 featureXML is written into it (the oracle's `FFC_noext`
/// is identical to its FFC_1 output).
#[test]
fn an_output_name_without_an_extension_is_accepted() {
    let dir = Workdir::new();
    let out = dir.file("FFC_noext");
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
    assert!(
        !outcome.err.contains("Invalid output file extension"),
        "{}",
        outcome.err
    );
    assert_reached_the_algorithm(&outcome);
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert_out_block(&outcome, FFC1_ALGORITHM_LINES);
    assert_ffc1_structure(&out, "FeatureFinderCentroided_1_input.mzML");
    assert_matches_ffc1_expectation(&out);
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

/// Oracle `FFC_seeds`: the retained FFC_1 output is used as the seed list. Its
/// eight features give the 24 seeds the source keeps after the charge filter,
/// which produce one candidate more than the computed seeds and one overlap
/// removal, and one seed that cannot be extended: `24 seeds`, `9 feature
/// candidates`, `Removed 1 overlapping features.`, `Could not extend seed: 1
/// times`, `8 features found.` The C++ Release run prints the same lines and
/// writes a feature map the decoded comparison cannot tell from the FFC_1
/// expectation.
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
    let out = dir.file("g.featureXML");
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
            &out,
        ],
    );
    assert_reached_the_algorithm(&outcome);
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert_out_block(&outcome, FFC_SEEDS_ALGORITHM_LINES);
    assert_ffc1_structure(&out, "FeatureFinderCentroided_1_input.mzML");
    assert_matches_ffc1_expectation(&out);
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
// The FAIMS closure (package B11)
//
// The corrected path has no whole-tool C++ oracle: the executed C++ tool exits
// 8 on every FAIMS input, because its voltage groups carry no per-MS-level
// ranges (`CPP-278`; re-executed in `../oracle/b11-faims`, eight runs on six
// FAIMS inputs, all rc 8 with `the value '1' was used but is not valid; No
// ranges for this MS level`). The oracle is built from the parts instead: each compensation
// voltage group is written as its own single-voltage mzML, with the FAIMS
// cvParam removed so that the C++ tool takes its non-FAIMS path, and the C++
// Release build is run on that file. The port's features for that group must
// equal what the Release build found (decision D6). Only the cross-voltage
// merge has no executed counterpart; it is pinned against the specification
// derived in `FaimsMergeFidelity` and against hand-derived cases whose numbers
// are written out below.
// ---------------------------------------------------------------------------

/// The features the C++ Release build found on one compensation-voltage group,
/// run on that group written as its own single-voltage mzML
/// (`../oracle/b11-faims`, `run.sh`, `cases.sh`; `OMP_NUM_THREADS=1`, `-test`).
fn faims_group_map(name: &str) -> FeatureMap {
    FileHandler::load_feature_map(fixture(name), &[FileType::FeatureXml]).unwrap()
}

/// The feature map the tool must write for a FAIMS input whose voltage groups
/// the C++ Release build produced one by one: the groups concatenated in
/// ascending voltage order, each feature annotated with its voltage, under the
/// FAIMS input's own `spectra_data`.
fn expected_faims_map(groups: &[(&str, f64)], input_basename: &str) -> FeatureMap {
    let mut expected = FeatureMap::new();
    for (name, volts) in groups {
        let group = faims_group_map(name);
        if expected.data_processing.is_empty() {
            expected.data_processing = group.data_processing.clone();
        }
        for mut feature in group.features {
            feature
                .metadata
                .insert("FAIMS_CV".into(), MetaValue::try_from(*volts).unwrap());
            expected.features.push(feature);
        }
    }
    expected
        .set_primary_ms_run_path(&[format!("file://{input_basename}")])
        .unwrap();
    expected
}

/// Compare a written featureXML with an expected map under decision D6:
/// decoded content, the upstream `FuzzyDiff` rule (ratio `1.01` or absolute
/// difference `0.01`) and generated identifiers skipped. The loose rule is the
/// upstream one; every caller also pins counts and the exact numbers that
/// matter.
fn assert_decoded_matches(path: &str, expected: &FeatureMap) {
    let actual = FileHandler::load_feature_map(path, &[FileType::FeatureXml]).unwrap();
    let options =
        decoded::DecodedOptions::new(decoded::Tolerance::new(1.01, 0.01)).ignoring_unique_ids();
    if let Err(mismatch) = decoded::compare_feature_maps(&actual, expected, &options) {
        panic!("{path}: {mismatch}");
    }
}

/// The `FAIMS_CV` of a feature, or `None` once a merge replaced it with the
/// merged-voltage list.
fn faims_cv(feature: &Feature) -> Option<f64> {
    feature
        .metadata
        .get("FAIMS_CV")
        .map(|v| v.as_f64().unwrap())
}

fn float_list(feature: &Feature, key: &str) -> Vec<f64> {
    feature.metadata[key].as_float_list().unwrap().to_vec()
}

/// A single compensation voltage: the split hands the whole input to one
/// group, so the algorithm sees exactly what it sees without FAIMS and the
/// features are `TOPP_FeatureFinderCentroided_1`'s, each annotated with the
/// voltage. Executed: the C++ Release build on the same 112 spectra is
/// `TOPP_FeatureFinderCentroided_1` itself, whose retained output this
/// compares against; the C++ tool on the FAIMS file itself exits 8 after
/// printing the two lines asserted here (oracle `faims_one_cv`).
///
/// The merge runs — there are eight FAIMS features — and merges nothing,
/// because every feature carries the same voltage.
#[test]
fn a_single_faims_voltage_finds_the_features_of_the_plain_input() {
    let dir = Workdir::new();
    let source = fs::read(ffc1_input()).unwrap();
    let input = dir.put("faims_one_cv.mzML", &derive_faims_one_cv(&source));
    let out = dir.file("one_cv.featureXML");
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
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert!(
        !outcome
            .out
            .contains(FeatureFinderCentroided::NO_FAIMS_MESSAGE),
        "the split reported no voltages:\n{}",
        outcome.out
    );
    let detected = FeatureFinderCentroided::faims_detected_message(1);
    let group = FeatureFinderCentroided::processing_group_message(-45.0, 112);
    let combined = FeatureFinderCentroided::combined_features_message(8);
    let merge = FeatureFinderCentroided::faims_merge_message(8, 8);
    assert_out_block(
        &outcome,
        &[
            detected.as_str(),
            group.as_str(),
            "Found 25 seeds for charge 2.",
            "Found 8 feature candidates for charge 2.",
            "Removed 0 overlapping features.",
            "",
            "Info: reasons for not finalizing a feature during its construction:",
            " - Invalid fit: Fitted model is bigger than 'max_rt_span': 1 times",
            "",
            "8 features found.",
            combined.as_str(),
            merge.as_str(),
        ],
    );

    let mut expected = ffc1_expected_map();
    for feature in &mut expected.features {
        feature
            .metadata
            .insert("FAIMS_CV".into(), MetaValue::try_from(-45.0).unwrap());
    }
    expected
        .set_primary_ms_run_path(&["file://faims_one_cv.mzML".to_owned()])
        .unwrap();
    assert_decoded_matches(&out, &expected);

    let written = FileHandler::load_feature_map(&out, &[FileType::FeatureXml]).unwrap();
    assert_eq!(written.features.len(), 8);
    for feature in &written.features {
        assert_eq!(faims_cv(feature), Some(-45.0));
        assert!(!feature.metadata.contains_key("merged_centroid_IMs"));
        assert!(!feature.metadata.contains_key("FAIMS_merge_count"));
    }
}

/// Two voltages, one algorithm run each. With `-faims_merge_features false`
/// the output is exactly the two groups, in ascending voltage order, so it can
/// be compared feature by feature with what the C++ Release build found on
/// each group on its own (oracle cases `group_m60` and `group_m45`: three and
/// two features, and the console lines asserted here).
#[test]
fn two_faims_voltages_reproduce_the_release_build_group_by_group() {
    let dir = Workdir::new();
    let source = fs::read(ffc1_input()).unwrap();
    let input = dir.put("faims_two_cv.mzML", &derive_faims_two_cv(&source));
    let out = dir.file("two_cv.featureXML");
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
            "-faims_merge_features",
            "false",
        ],
    );
    outcome.assert_exit(ExitCode::ExecutionOk);
    let detected = FeatureFinderCentroided::faims_detected_message(2);
    let first = FeatureFinderCentroided::processing_group_message(-60.0, 56);
    let second = FeatureFinderCentroided::processing_group_message(-45.0, 56);
    assert_out_block(
        &outcome,
        &[
            detected.as_str(),
            first.as_str(),
            "Found 14 seeds for charge 2.",
            "Found 3 feature candidates for charge 2.",
            "Removed 0 overlapping features.",
            "",
            "Info: reasons for not finalizing a feature during its construction:",
            " - Could not extend seed: 6 times",
            "",
            "3 features found.",
            second.as_str(),
            "Found 14 seeds for charge 2.",
            "Found 2 feature candidates for charge 2.",
        ],
    );
    assert_out_contains_line(&outcome, "2 features found.");
    assert_out_contains_line(
        &outcome,
        &FeatureFinderCentroided::combined_features_message(5),
    );
    assert!(
        !outcome.out.contains("FAIMS feature merge:"),
        "-faims_merge_features false still merged:\n{}",
        outcome.out
    );

    assert_decoded_matches(
        &out,
        &expected_faims_map(
            &[
                ("faims_group_m60.featureXML", -60.0),
                ("faims_group_m45.featureXML", -45.0),
            ],
            "faims_two_cv.mzML",
        ),
    );
    let written = FileHandler::load_feature_map(&out, &[FileType::FeatureXml]).unwrap();
    let voltages: Vec<Option<f64>> = written.features.iter().map(faims_cv).collect();
    assert_eq!(
        voltages,
        [
            Some(-60.0),
            Some(-60.0),
            Some(-60.0),
            Some(-45.0),
            Some(-45.0)
        ]
    );
}

/// The merge of the same two voltages. Hand-derived from the executed group
/// features (`../oracle/b11-faims`, `results/out/group_m45`, `group_m60`) and
/// the specification of `FaimsMergeFidelity::Corrected`:
///
/// | analyte | −45 V | −60 V | survivor | merged intensity |
/// |---|---|---|---|---|
/// | 4389.11 s, 648.257 | 45181.723 | 44601.215 | −45 | 89782.9375 |
/// | 4300.95 s, 651.760 | 35109.383 | 34660.066 | −45 | 69769.453125 |
/// | 4278.07 s, 653.770 | — | 19216.809 | −60 | 19216.80859375 |
///
/// Every pair is within 5 s and 0.05 Da and has charge 2, the third feature is
/// 23 s and 2 Da away from the second and merges with nothing. The survivor is
/// the member of higher intensity, and the sum is `f32(f64(a) + f64(b))`, the
/// source's `setIntensity(double + double)` into a `float`. The survivors keep
/// the descending-intensity order of the merge.
// The literals below are written with the digits the C++ featureXML writer
// printed, and the `f32` sums with every digit the value has, so that each one
// can be checked against the oracle file by eye. Clippy would shorten them to
// the fewest digits that name the same float.
#[allow(clippy::excessive_precision)]
#[test]
fn the_faims_merge_joins_the_two_voltages_of_each_analyte() {
    let dir = Workdir::new();
    let source = fs::read(ffc1_input()).unwrap();
    let input = dir.put("faims_two_cv.mzML", &derive_faims_two_cv(&source));
    let out = dir.file("two_cv_merged.featureXML");
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
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert_out_contains_line(
        &outcome,
        &FeatureFinderCentroided::combined_features_message(5),
    );
    assert_out_contains_line(
        &outcome,
        &FeatureFinderCentroided::faims_merge_message(5, 3),
    );

    let written = FileHandler::load_feature_map(&out, &[FileType::FeatureXml]).unwrap();
    assert_eq!(written.features.len(), 3);
    let intensities: Vec<u32> = written
        .features
        .iter()
        .map(|f| f.base.intensity.to_bits())
        .collect();
    assert_eq!(
        intensities,
        [
            89782.9375f32.to_bits(),
            69769.453125f32.to_bits(),
            19216.80859375f32.to_bits()
        ]
    );
    for (index, feature) in written.features.iter().take(2).enumerate() {
        assert_eq!(faims_cv(feature), None, "features[{index}] kept FAIMS_CV");
        assert_eq!(
            float_list(feature, "merged_centroid_IMs"),
            [-45.0, -60.0],
            "features[{index}] merged voltages"
        );
        assert_eq!(
            feature.metadata["FAIMS_merge_count"].as_i64().unwrap(),
            2,
            "features[{index}] merge count"
        );
    }
    assert_eq!(
        float_list(&written.features[0], "merged_centroid_rts"),
        [4389.113192580538453, 4389.204616887192060]
    );
    assert_eq!(
        float_list(&written.features[0], "merged_centroid_mzs"),
        [648.257466504270724, 648.255764126700569]
    );
    // The third analyte was found at one voltage only, so it is untouched.
    assert_eq!(faims_cv(&written.features[2]), Some(-60.0));
    assert!(
        !written.features[2]
            .metadata
            .contains_key("merged_centroid_IMs")
    );
}

/// Three voltages, with `mass_trace:min_spectra 5` so that every group finds
/// the same analytes: the case that separates the corrected merge from the
/// source's, because six clusters hold **three** features each. The source
/// merges a survivor once and then refuses (`CPP-283`), which would leave those
/// six clusters as twelve features; the corrected merge collapses each to one.
///
/// The three groups are executed (`../oracle/b11-faims`, `cases2.sh`,
/// `group3s5_m45`, `group3s5_m60`, `group3s5_m70`: 8, 8 and 7 features), and
/// `-faims_merge_features false` reproduces them exactly. The merged
/// expectation is hand-derived from those 23 features by the specification and
/// written out in the table below; the intensity of each survivor is the `f32`
/// running sum of its cluster.
// The literals below are written with the digits the C++ featureXML writer
// printed, and the `f32` sums with every digit the value has, so that each one
// can be checked against the oracle file by eye. Clippy would shorten them to
// the fewest digits that name the same float.
#[allow(clippy::excessive_precision)]
#[test]
fn three_faims_voltages_collapse_a_cluster_of_three_into_one_feature() {
    let dir = Workdir::new();
    let source = fs::read(ffc1_input()).unwrap();
    let derived = derive_faims(&source, &["-45", "-60", "-70"]);
    assert_digest(
        &derived,
        "7edd8c78aebe6f6553498bbce3c0586bb6817fa7",
        "faims_three_cv",
    );
    let input = dir.put("faims_three_cv.mzML", &derived);
    let ini = text(ffc1_ini());
    let short_trace = ["-algorithm:mass_trace:min_spectra", "5"];

    // Without the merge: the three executed groups, in ascending voltage order.
    let unmerged = dir.file("three_cv_unmerged.featureXML");
    let mut arguments = vec![
        "-test",
        "-ini",
        &ini,
        "-in",
        &input,
        "-out",
        &unmerged,
        "-faims_merge_features",
        "false",
    ];
    arguments.extend_from_slice(&short_trace);
    let outcome = run_in(&dir, &arguments);
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert_out_contains_line(
        &outcome,
        &FeatureFinderCentroided::faims_detected_message(3),
    );
    for (volts, spectra) in [(-70.0, 37), (-60.0, 37), (-45.0, 38)] {
        assert_out_contains_line(
            &outcome,
            &FeatureFinderCentroided::processing_group_message(volts, spectra),
        );
    }
    assert_out_contains_line(
        &outcome,
        &FeatureFinderCentroided::combined_features_message(23),
    );
    assert_decoded_matches(
        &unmerged,
        &expected_faims_map(
            &[
                ("faims_group3s5_m70.featureXML", -70.0),
                ("faims_group3s5_m60.featureXML", -60.0),
                ("faims_group3s5_m45.featureXML", -45.0),
            ],
            "faims_three_cv.mzML",
        ),
    );

    // With the merge: ten survivors.
    let merged = dir.file("three_cv_merged.featureXML");
    let mut arguments = vec!["-test", "-ini", &ini, "-in", &input, "-out", &merged];
    arguments.extend_from_slice(&short_trace);
    let outcome = run_in(&dir, &arguments);
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert_out_contains_line(
        &outcome,
        &FeatureFinderCentroided::faims_merge_message(23, 10),
    );

    let written = FileHandler::load_feature_map(&merged, &[FileType::FeatureXml]).unwrap();
    assert_eq!(written.features.len(), 10);
    // survivor voltage, retention time, m/z, merged voltages, intensity
    let expected: &[(f64, f64, f64, &[f64], f32)] = &[
        (
            -70.0,
            4407.335644965823121,
            646.237161048441294,
            &[-70.0, -60.0, -45.0],
            152660.390625,
        ),
        (
            -45.0,
            4388.999406220842502,
            648.257121711850118,
            &[-45.0, -70.0, -60.0],
            134546.578125,
        ),
        (
            -70.0,
            4301.190310864986714,
            651.757210503288547,
            &[-70.0, -60.0, -45.0],
            103111.703125,
        ),
        // The only cluster of this package whose `f32` running sum could
        // depend on the order in which the survivor absorbs its two partners.
        // In the order below it is exact twice over: 20089.396484375 +
        // 19694.748046875 = 39784.14453125 and + 18861.0546875 =
        // 58645.19921875, both representable. The other order rounds twice —
        // 20089.396484375 + 18861.0546875 = 38950.451171875 is a tie, rounded
        // to even as 38950.453125, and the second sum ties again — and would
        // give 58645.203125. Which order the merge takes is the quadtree's
        // query order, **derived** rather than measured: the derivation is in
        // the note below this table.
        (
            -70.0,
            4278.163376914249966,
            653.775964531109253,
            &[-70.0, -45.0, -60.0],
            58645.19921875,
        ),
        (
            -60.0,
            4201.935576947686059,
            652.764241594621922,
            &[-60.0, -70.0, -45.0],
            38816.55859375,
        ),
        (
            -45.0,
            4221.610088053402251,
            646.766440532254705,
            &[-45.0, -60.0, -70.0],
            25299.296875,
        ),
        (
            -60.0,
            4183.539151442289949,
            654.782096916037176,
            &[-60.0, -70.0],
            14590.9267578125,
        ),
    ];
    let by_position = |rt: f64, mz: f64| {
        written
            .features
            .iter()
            .find(|f| (f.base.rt - rt).abs() < 1e-6 && (f.base.mz - mz).abs() < 1e-9)
            .unwrap_or_else(|| panic!("no survivor at {rt} s, {mz} Da"))
    };
    for &(_, rt, mz, voltages, intensity) in expected {
        let feature = by_position(rt, mz);
        assert_eq!(faims_cv(feature), None, "{rt}: kept FAIMS_CV");
        assert_eq!(
            float_list(feature, "merged_centroid_IMs"),
            voltages,
            "{rt}: merged voltages"
        );
        assert_eq!(
            usize::try_from(feature.metadata["FAIMS_merge_count"].as_i64().unwrap()).unwrap(),
            voltages.len(),
            "{rt}: merge count"
        );
        assert_eq!(
            feature.base.intensity.to_bits(),
            intensity.to_bits(),
            "{rt}: merged intensity"
        );
    }

    // The absorption order of the 653.776 Da cluster, on which its running sum
    // hangs, is derived and not read off this run. `../oracle/b11-faims/quadorder.py`
    // re-implements `Plan::new`, `Plan::feature_box`, the stable sort by
    // intensity and the quadtree's `add`/`split`/`quadrant`/`query_node` in
    // `f32`, reads the three C++ Release group fixtures and executes no Rust.
    // Its answer for that cluster is `[-70, -45, -60]`, the row above. The same
    // run reproduces the other six rows of the table, the survivor count and the
    // three single-voltage intensities below, which is what shows the
    // re-implementation to be of this algorithm; `../oracle/b11-faims/results/quadorder.txt`
    // holds its output.

    // Three analytes were found at one voltage only and keep their `FAIMS_CV`.
    let untouched: Vec<f64> = written.features.iter().filter_map(faims_cv).collect();
    assert_eq!(untouched.len(), 3);
    let single: Vec<u32> = written
        .features
        .iter()
        .filter(|f| faims_cv(f).is_some())
        .map(|f| f.base.intensity.to_bits())
        .collect();
    assert_eq!(
        single,
        [
            4783.8017578125f32.to_bits(),
            4171.59326171875f32.to_bits(),
            3829.223388671875f32.to_bits()
        ]
    );
}

/// Oracle `c5_faims_partial_cv`: a compensation voltage on half of the MS1
/// spectra is still FAIMS input. The split assigns the annotated spectra to the
/// voltage group and skips the others, each with the warning the executed C++
/// wrote to **stderr**, folded by the log stream's line cache into one line and
/// one `occurred 56 times` report. The group is exactly the 56 spectra the
/// oracle case `group_m45` was run on, so its features are that case's.
#[test]
fn a_faims_voltage_on_some_spectra_leaves_one_group_and_a_warning_per_skip() {
    let dir = Workdir::new();
    let source = fs::read(ffc1_input()).unwrap();
    let derived = derive_faims(&source, &["-45", ""]);
    assert_digest(
        &derived,
        "bf274d61a77be9b2af64e51ab4f5bf266452aab9",
        "faims_partial_cv",
    );
    let input = dir.put("faims_partial_cv.mzML", &derived);
    let out = dir.file("partial_cv.featureXML");
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
    outcome.assert_exit(ExitCode::ExecutionOk);
    let skipped =
        "Skipping spectrum without FAIMS CV (no prior FAIMS CV context or unexpected layout).";
    assert_eq!(
        outcome.err.lines().collect::<Vec<&str>>(),
        [skipped, &format!("<{skipped}> occurred 56 times")]
    );
    assert_out_contains_line(
        &outcome,
        &FeatureFinderCentroided::faims_detected_message(1),
    );
    assert_out_contains_line(
        &outcome,
        &FeatureFinderCentroided::processing_group_message(-45.0, 56),
    );
    assert_decoded_matches(
        &out,
        &expected_faims_map(
            &[("faims_group_m45.featureXML", -45.0)],
            "faims_partial_cv.mzML",
        ),
    );
}

/// The two upstream FAIMS fixtures (test-data `0cb15f2`), which spell the volt
/// unit `UO:000218` (`CPP-240`) and are too short for the default
/// `mass_trace:min_spectra` of 10, so every group gives an empty map — as the
/// C++ Release build does on each group written out on its own (oracle
/// `testdata_cvm65`, `interleaved_cvm45`, `interleaved_cvm60`: exit 0, `0
/// features found.`).
///
/// `FAIMS_test_data.mzML` holds one MS1 and one MS2 spectrum at −65 V, and the
/// loader keeps MS level 1 only, so the group holds one spectrum — the count
/// the executed C++ printed. The interleaved file is profile data, so without
/// `-force` the profile check still fires before the split
/// (oracle `FFC_faims_interleaved_noforce`).
#[test]
fn the_upstream_faims_fixtures_run_to_an_empty_map() {
    let dir = Workdir::new();
    let test_data = text(mobility("FAIMS_test_data.mzML"));
    let interleaved = text(mobility("FAIMS_CV-60C_V-45_Interleaved.mzML"));

    let out = dir.file("test_data.featureXML");
    let outcome = run_in(&dir, &["-test", "-in", &test_data, "-out", &out]);
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert_out_contains_line(
        &outcome,
        &FeatureFinderCentroided::faims_detected_message(1),
    );
    assert_out_contains_line(
        &outcome,
        &FeatureFinderCentroided::processing_group_message(-65.0, 1),
    );
    assert_out_contains_line(&outcome, "0 features found.");
    assert_out_contains_line(
        &outcome,
        &FeatureFinderCentroided::combined_features_message(0),
    );
    assert!(
        FileHandler::load_feature_map(&out, &[FileType::FeatureXml])
            .unwrap()
            .features
            .is_empty()
    );

    let out = dir.file("interleaved.featureXML");
    let outcome = run_in(&dir, &["-test", "-in", &interleaved, "-out", &out]);
    outcome.assert_exit(ExitCode::UnknownError);
    outcome.assert_err_contains(FeatureFinderCentroided::PROFILE_DATA_MESSAGE);
    assert!(!Path::new(&out).exists());

    let outcome = run_in(
        &dir,
        &["-test", "-in", &interleaved, "-out", &out, "-force"],
    );
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert_out_contains_line(
        &outcome,
        &FeatureFinderCentroided::faims_detected_message(2),
    );
    for volts in [-60.0, -45.0] {
        assert_out_contains_line(
            &outcome,
            &FeatureFinderCentroided::processing_group_message(volts, 6),
        );
    }
    assert_out_contains_line(
        &outcome,
        &FeatureFinderCentroided::combined_features_message(0),
    );
    assert!(
        FileHandler::load_feature_map(&out, &[FileType::FeatureXml])
            .unwrap()
            .features
            .is_empty()
    );
}

/// The per-group seed filter (`FeatureFinderCentroided.cpp:258-281`): a seed
/// with a `FAIMS_CV` within 0.01 V of the group joins it, a seed without one
/// joins every group, and without FAIMS input the list is used as it is.
#[test]
fn seeds_are_filtered_by_compensation_voltage() {
    let mut seeds = FeatureMap::new();
    let mut annotated = |cv: Option<f64>, rt: f64| {
        let mut feature = Feature::new(rt, 0.0, 0.0);
        if let Some(cv) = cv {
            feature
                .metadata
                .insert("FAIMS_CV".into(), MetaValue::try_from(cv).unwrap());
        }
        seeds.features.push(feature);
    };
    annotated(Some(-45.0), 1.0);
    annotated(Some(-60.0), 2.0);
    annotated(None, 3.0);
    annotated(Some(-45.005), 4.0);
    annotated(Some(-45.02), 5.0);

    let times = |map: &FeatureMap| -> Vec<f64> { map.features.iter().map(|f| f.base.rt).collect() };
    assert_eq!(
        times(&FeatureFinderCentroided::seeds_of_group(&seeds, true, -45.0).unwrap()),
        [1.0, 3.0, 4.0]
    );
    assert_eq!(
        times(&FeatureFinderCentroided::seeds_of_group(&seeds, true, -60.0).unwrap()),
        [2.0, 3.0]
    );
    // Not FAIMS input: the whole list, whatever the voltage argument says.
    assert_eq!(
        times(&FeatureFinderCentroided::seeds_of_group(&seeds, false, f64::NAN).unwrap()),
        [1.0, 2.0, 3.0, 4.0, 5.0]
    );
    // An empty list stays empty rather than being filtered.
    assert!(
        FeatureFinderCentroided::seeds_of_group(&FeatureMap::new(), true, -45.0)
            .unwrap()
            .features
            .is_empty()
    );
    // A non-numeric FAIMS_CV is a ConversionError in the source.
    let mut bad = FeatureMap::new();
    let mut feature = Feature::new(0.0, 0.0, 0.0);
    feature
        .metadata
        .insert("FAIMS_CV".into(), MetaValue::from("text"));
    bad.features.push(feature);
    assert!(FeatureFinderCentroided::seeds_of_group(&bad, true, -45.0).is_err());
}

// ---------------------------------------------------------------------------
// TOPP_FeatureFinderCentroided_1 itself
// ---------------------------------------------------------------------------

/// The registered upstream workflow (test-data `topp/CMakeLists.txt:425-428`,
/// oracle `TOPP_FeatureFinderCentroided_1`): exit 0 and the eight features of
/// `FeatureFinderCentroided_1_1_output.featureXML`.
///
/// The decoded comparison (decision D6) uses the upstream `FuzzyDiff` rule, but
/// the loose rule is not what this case rests on: the numbers below are the
/// measured gap to the C++ **Release** build
/// `openms4-release-bc9cc12-c19e494-174b576` on the same input and INI, which
/// is far tighter than `FuzzyDiff` and tighter than the `1e-9` relative bound
/// package B10 asks for on the fitted fields:
///
/// | field | worst gap to C++ Release | C++ Debug vs C++ Release |
/// |---|---|---|
/// | convex-hull `rt` and `mz` | 0 (bit-identical) | 0 |
/// | feature `mz` | 0 (bit-identical) | 0 |
/// | feature `rt` | `5.5e-13` relative | `2.2e-13` |
/// | `intensity`, `FWHM` | 0 (identical as `f32`) | 0 |
/// | `score_fit` | `2.2e-10` relative | `9.1e-11` |
/// | `score_correlation` | `7.7e-12` relative | `3.1e-12` |
/// | `overallquality` | agrees to the six decimals C++ prints | — |
///
/// So the remaining difference is the last bits of the Levenberg-Marquardt fit,
/// of the same order as the C++ build's own Debug-to-Release spread, and every
/// integral and structural field agrees exactly. The retained expectation is
/// itself not bit-reproducible by current C++ (its `intensity` is printed with
/// one digit fewer), which is why the tight numbers are quoted against the
/// executed build and the in-repo assertion is the decoded `FuzzyDiff` one.
///
/// The `1e-9` comparison against the C1 oracle output itself, the thread sweep
/// and the `-algorithm:fit:max_iterations` boundary stay package B10's.
#[test]
fn the_upstream_workflow_matches_the_retained_expectation() {
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
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert_out_block(&outcome, FFC1_ALGORITHM_LINES);
    assert!(Path::new(&out).exists(), "the output file is written");
    assert_ffc1_structure(&out, "FeatureFinderCentroided_1_input.mzML");
    assert_matches_ffc1_expectation(&out);

    // The fitted fields, against the retained expectation, at the tolerance the
    // measurement above supports rather than the FuzzyDiff one. The expectation
    // prints `intensity` with seven significant digits, so intensity is checked
    // as the `f32` it is stored as.
    let actual = FileHandler::load_feature_map(&out, &[FileType::FeatureXml]).unwrap();
    let expected = ffc1_expected_map();
    for (index, (a, e)) in actual.features.iter().zip(&expected.features).enumerate() {
        assert!(
            (a.base.rt - e.base.rt).abs() <= 1e-9 * e.base.rt.abs(),
            "features[{index}].rt: {} vs {}",
            a.base.rt,
            e.base.rt
        );
        assert_eq!(a.base.mz, e.base.mz, "features[{index}].mz");
        // The retained file prints seven significant digits of the `f32`
        // intensity where current C++ prints its full decimal expansion (B10's
        // note: the fixture is not bit-reproducible even by C++), so the two
        // reconstructed `f32` values differ by up to 1.6e-7 relative. The
        // assertion is therefore that both print the same seven digits — the
        // whole information the fixture carries. Against the executed C++
        // Release run the intensity is bit-identical as an `f32`.
        assert_eq!(
            format!("{:.6e}", a.base.intensity),
            format!("{:.6e}", e.base.intensity),
            "features[{index}].intensity"
        );
        let (a_fit, e_fit) = (meta_f64(a, "score_fit"), meta_f64(e, "score_fit"));
        assert!(
            (a_fit - e_fit).abs() <= 1e-9 * e_fit.abs(),
            "features[{index}].score_fit: {a_fit} vs {e_fit}"
        );
        let (a_cor, e_cor) = (
            meta_f64(a, "score_correlation"),
            meta_f64(e, "score_correlation"),
        );
        assert!(
            (a_cor - e_cor).abs() <= 1e-9 * e_cor.abs(),
            "features[{index}].score_correlation: {a_cor} vs {e_cor}"
        );
        assert_eq!(
            a.convex_hulls
                .iter()
                .map(ConvexHull2D::hull_points)
                .collect::<Vec<_>>(),
            e.convex_hulls
                .iter()
                .map(ConvexHull2D::hull_points)
                .collect::<Vec<_>>(),
            "features[{index}] hull points are not bit-identical"
        );
    }
}

/// A `float`-valued metadata entry of a feature.
fn meta_f64(feature: &Feature, key: &str) -> f64 {
    feature
        .metadata
        .get(key)
        .unwrap_or_else(|| panic!("{key} is present"))
        .as_f64()
        .unwrap_or_else(|error| panic!("{key} is a number: {error}"))
}

/// `-threads` reaches the seed loop, and the determinism contract holds through
/// the tool: 1, 2, 4, 8 and 0 (every core) write byte-identical output, unique
/// ids included, because the loop returns its results in seed order and every
/// later step is serial. The C++ oracle asserts the same across
/// `FFC_1_threads_0/1/2/4/8`.
#[test]
fn the_output_is_byte_identical_at_every_thread_count() {
    let dir = Workdir::new();
    let mut reference: Option<Vec<u8>> = None;
    for threads in ["1", "2", "4", "8", "0"] {
        let out = dir.file(&format!("threads_{threads}.featureXML"));
        let outcome = run_in(
            &dir,
            &[
                "-test",
                "-ini",
                &text(ffc1_ini()),
                "-in",
                &text(ffc1_input()),
                "-threads",
                threads,
                "-out",
                &out,
            ],
        );
        outcome.assert_exit(ExitCode::ExecutionOk);
        assert_out_block(&outcome, FFC1_ALGORITHM_LINES);
        let written = fs::read(&out).unwrap();
        match &reference {
            None => reference = Some(written),
            Some(first) => assert!(
                *first == written,
                "-threads {threads} wrote a different file"
            ),
        }
    }
}

/// The same contract on the path this package adds, which the case above does
/// not reach: three compensation voltages are three algorithm runs, three seed
/// filters and a cross-voltage merge whose intensity is an `f32` running sum
/// over a quadtree query order. Nothing there may depend on the thread count,
/// so 1, 2, 4, 8 and 0 must again write one byte-identical file — unique ids
/// included, since the generator is drawn from in the merge as well.
#[test]
fn the_faims_output_is_byte_identical_at_every_thread_count() {
    let dir = Workdir::new();
    let source = fs::read(ffc1_input()).unwrap();
    let derived = derive_faims(&source, &["-45", "-60", "-70"]);
    assert_digest(
        &derived,
        "7edd8c78aebe6f6553498bbce3c0586bb6817fa7",
        "faims_three_cv",
    );
    let input = dir.put("faims_three_cv.mzML", &derived);
    let ini = text(ffc1_ini());
    let mut reference: Option<Vec<u8>> = None;
    for threads in ["1", "2", "4", "8", "0"] {
        let out = dir.file(&format!("faims_threads_{threads}.featureXML"));
        let outcome = run_in(
            &dir,
            &[
                "-test",
                "-ini",
                &ini,
                "-in",
                &input,
                "-threads",
                threads,
                "-out",
                &out,
                "-algorithm:mass_trace:min_spectra",
                "5",
            ],
        );
        outcome.assert_exit(ExitCode::ExecutionOk);
        assert_out_contains_line(
            &outcome,
            &FeatureFinderCentroided::faims_merge_message(23, 10),
        );
        let written = fs::read(&out).unwrap();
        match &reference {
            None => reference = Some(written),
            Some(first) => assert!(
                *first == written,
                "-threads {threads} wrote a different file"
            ),
        }
    }
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

// ---------------------------------------------------------------------------
// -algorithm:write_debug (package ffap-instrumentation)
// ---------------------------------------------------------------------------
//
// Every case below was executed with the C++ Release FeatureFinderCentroided
// (`../oracle/ffap-instr-completion`, `tool_cases.py`, three repetitions each
// with OMP_NUM_THREADS=1): the exit status, the console output and every file
// the run left behind are compared, log.txt byte for byte (by SHA-1) and the
// featureXML and mzML files by decoded content (decision D6), with the exact
// unique ids the seeded -test generator draws.

/// A fixture of the instrumentation package.
fn instrumentation(name: &str) -> PathBuf {
    repository("tests/data/feature_finder_picked_instrumentation").join(name)
}

/// The executed size and SHA-1 of one debug file (`debug_digests.tsv`).
fn executed_digest(case: &str, file: &str) -> (usize, String) {
    let table = fs::read_to_string(instrumentation("debug_digests.tsv")).unwrap();
    for line in table.lines().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields[0] == case && fields[1] == file {
            return (fields[2].parse().unwrap(), fields[3].to_owned());
        }
    }
    panic!("no executed digest for {case} {file}");
}

fn assert_executed_file(path: &Path, case: &str, file: &str) {
    let data = fs::read(path).unwrap();
    let (bytes, sha1) = executed_digest(case, file);
    assert_eq!(data.len(), bytes, "{case} {file}: size");
    assert_eq!(digest(&data), sha1, "{case} {file}: content");
}

fn debug_d6() -> decoded::DecodedOptions {
    decoded::DecodedOptions::new(decoded::Tolerance::from_settings(
        &fuzzy::FuzzyDiffSettings::upstream().unwrap(),
    ))
    .ignoring_unique_ids()
}

/// A featureXML file, decoded; `.gz` fixtures included.
fn decoded_features(path: &Path) -> FeatureMap {
    openms::format::featurexml::load(path).unwrap()
}

fn assert_same_features(actual: &Path, expected: &Path) {
    if let Err(mismatch) = decoded::compare_feature_maps(
        &decoded_features(actual),
        &decoded_features(expected),
        &debug_d6(),
    ) {
        panic!("{}: {mismatch}", actual.display());
    }
}

/// An mzML file with NaN scores allowed, as the source reads it; `.gz`
/// fixtures are decompressed first.
fn decoded_debug_input(path: &Path) -> MSExperiment {
    let mut bytes = fs::read(path).unwrap();
    if path.extension().is_some_and(|e| e == "gz") {
        use std::io::Read;
        let mut xml = Vec::new();
        flate2::read::GzDecoder::new(bytes.as_slice())
            .read_to_end(&mut xml)
            .unwrap();
        bytes = xml;
    }
    let options = openms::format::mzml::ReadOptions {
        source_nonfinite_float_arrays: true,
        ..Default::default()
    };
    openms::format::mzml::read_with_options(bytes.as_slice(), &options).unwrap()
}

/// The debug input written by the port against the executed one: the same
/// spectra, peaks and score arrays, bit for bit (NaN bits and the overall
/// scores of the reference build's `powf` included).
fn assert_same_debug_input(actual: &Path, expected: &Path) {
    let a = decoded_debug_input(actual);
    let e = decoded_debug_input(expected);
    assert_eq!(a.spectra.len(), e.spectra.len());
    for (x, y) in a.spectra.iter().zip(&e.spectra) {
        assert_eq!(x.native_id, y.native_id);
        assert_eq!(x.rt.to_bits(), y.rt.to_bits());
        assert_eq!(x.peaks, y.peaks);
        let names = |s: &MSSpectrum| -> Vec<String> {
            s.float_data_arrays.iter().map(|a| a.name.clone()).collect()
        };
        assert_eq!(names(x), names(y));
        for (u, v) in x.float_data_arrays.iter().zip(&y.float_data_arrays) {
            assert_eq!(u.data.len(), v.data.len(), "{} {}", x.native_id, u.name);
            for (p, q) in u.data.iter().zip(&v.data) {
                assert_eq!(
                    p.to_bits(),
                    q.to_bits(),
                    "{} {}: {p} against {q}",
                    x.native_id,
                    u.name
                );
            }
        }
    }
}

/// The executed console block of a debug run: from the FAIMS line to the end,
/// without TOPPBase's closing timing line and the fatal block of a terminated
/// run.
fn executed_block(name: &str) -> Vec<String> {
    let text = fs::read_to_string(instrumentation(name)).unwrap();
    let start = text
        .find(FeatureFinderCentroided::NO_FAIMS_MESSAGE)
        .unwrap();
    let mut lines = Vec::new();
    for line in text[start..].lines() {
        if line.starts_with("FeatureFinderCentroided took ")
            || line.starts_with("-----------------")
        {
            break;
        }
        lines.push(line.to_owned());
    }
    while lines.last().is_some_and(String::is_empty) {
        lines.pop();
    }
    lines
}

fn port_block(outcome: &Outcome) -> Vec<String> {
    let start = outcome
        .out
        .find(FeatureFinderCentroided::NO_FAIMS_MESSAGE)
        .unwrap();
    let mut lines: Vec<String> = outcome.out[start..].lines().map(str::to_owned).collect();
    while lines.last().is_some_and(String::is_empty) {
        lines.pop();
    }
    lines
}

fn debug_args<'a>(dir: &'a Workdir, extra: &[&'a str]) -> (String, String, Vec<String>) {
    let ini = text(ffc1_ini());
    let input = text(ffc1_input());
    let mut args = vec![
        "-test".to_owned(),
        "-ini".to_owned(),
        ini.clone(),
        "-in".to_owned(),
        input.clone(),
        "-out".to_owned(),
        dir.file("out.featureXML"),
        "-algorithm:write_debug".to_owned(),
    ];
    args.extend(extra.iter().map(|s| (*s).to_owned()));
    (ini, input, args)
}

fn run_args(dir: &Workdir, args: &[String]) -> Outcome {
    let refs: Vec<&str> = args.iter().map(String::as_str).collect();
    run_in(dir, &refs)
}

/// The id `-test` makes the source draw first: the debug abort map's.
const FIRST_TEST_MODE_ID: u64 = 5_233_264_595_117_471_314;
/// The output map's id after the abort map's draw (the executed a1 and a2
/// outputs); without debug it is the second draw, 4835329514588776807.
const OUTPUT_ID_AFTER_DEBUG: u64 = 17_749_660_155_506_638_460;

/// Executed case x0: `write_debug` is a true/false string, which TOPPBase
/// registers as a flag, so a value after it is a command-line error.
#[test]
fn write_debug_is_a_flag_on_the_command_line() {
    let dir = Workdir::new();
    let (_, _, mut args) = debug_args(&dir, &[]);
    args.push("true".into());
    let outcome = run_args(&dir, &args);
    outcome.assert_exit(ExitCode::IllegalParameters);
    let executed = fs::read_to_string(instrumentation("tool_x0_stderr.txt")).unwrap();
    let first = executed.lines().next().unwrap();
    assert_eq!(
        first,
        "Invalid parameter values (InvalidParameter): Command line error: Trailing arguments after flag '-algorithm:write_debug': true. Aborting!"
    );
    outcome.assert_err_contains(first);
    assert!(!dir.path().join("debug").exists());
}

/// Executed cases a1 and a2: no seed reaches the fit, so the debug run
/// completes. The console, log.txt, the seed map, the abort map, the input
/// with its score arrays and the output match the Release build.
#[test]
fn a_completed_debug_run_writes_the_executed_files() {
    for (case, extra, stdout, input, seeds, aborts, out) in [
        (
            "a1",
            ["-algorithm:mass_trace:min_spectra", "1"],
            "tool_a1_stdout.txt",
            "a1_input.mzML.gz",
            "empty_seed_map.featureXML",
            "empty_abort_map.featureXML",
            "tool_a1_out.featureXML",
        ),
        (
            "a2",
            ["-algorithm:feature:min_isotope_fit", "1.0"],
            "tool_a2_stdout.txt",
            "a2_input.mzML.gz",
            "ffc1_seed_map.featureXML.gz",
            "a2_abort_map.featureXML.gz",
            "tool_a2_out.featureXML",
        ),
    ] {
        let dir = Workdir::new();
        let (_, _, args) = debug_args(&dir, &extra);
        let outcome = run_args(&dir, &args);
        outcome.assert_exit(ExitCode::ExecutionOk);
        assert_eq!(port_block(&outcome), executed_block(stdout), "{case}");
        let debug = dir.path().join("debug");
        assert!(debug.join("features").is_dir());
        assert_eq!(fs::read_dir(debug.join("features")).unwrap().count(), 0);
        assert_executed_file(&debug.join("log.txt"), case, "log.txt");
        assert_same_features(&debug.join("seeds_2.featureXML"), &instrumentation(seeds));
        let abort_path = debug.join("abort_reasons.featureXML");
        assert_same_features(&abort_path, &instrumentation(aborts));
        assert_eq!(decoded_features(&abort_path).unique_id, FIRST_TEST_MODE_ID);
        assert_same_debug_input(&debug.join("input.mzML"), &instrumentation(input));
        let output = PathBuf::from(dir.file("out.featureXML"));
        assert_same_features(&output, &instrumentation(out));
        assert_eq!(decoded_features(&output).unique_id, OUTPUT_ID_AFTER_DEBUG);
        let mut entries: Vec<String> = fs::read_dir(&debug)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        entries.sort();
        assert_eq!(
            entries,
            [
                "abort_reasons.featureXML",
                "features",
                "input.mzML",
                "log.txt",
                "seeds_2.featureXML"
            ]
        );
    }
}

/// Executed case a3: four scans and all four charges, no seed: one seed map
/// per charge, and the repeated store line the OpenMS log stream reports
/// once, with its count, when a later line pushes it out.
#[test]
fn a_debug_run_on_a_short_input_writes_the_executed_files() {
    let dir = Workdir::new();
    let input = text(fixture("FileConverter_31_output.mzML"));
    let out = dir.file("out.featureXML");
    let outcome = run_in(
        &dir,
        &[
            "-test",
            "-in",
            &input,
            "-out",
            &out,
            "-force",
            "-algorithm:write_debug",
        ],
    );
    outcome.assert_exit(ExitCode::ExecutionOk);
    let block = port_block(&outcome);
    assert_eq!(block, executed_block("tool_a3_stdout.txt"));
    assert!(block.contains(
        &"<FeatureXMLHandler::store():  found 1 invalid unique ids> occurred 4 times".to_owned()
    ));
    let debug = dir.path().join("debug");
    assert_executed_file(&debug.join("log.txt"), "a3", "log.txt");
    for charge in 1..=4 {
        assert_same_features(
            &debug.join(format!("seeds_{charge}.featureXML")),
            &instrumentation("empty_seed_map.featureXML"),
        );
    }
    assert_same_features(
        &debug.join("abort_reasons.featureXML"),
        &instrumentation("empty_abort_map.featureXML"),
    );
    assert_same_debug_input(
        &debug.join("input.mzML"),
        &instrumentation("a3_input.mzML.gz"),
    );
    assert_same_features(Path::new(&out), &instrumentation("tool_a3_out.featureXML"));
}

/// Executed case b1: the first seed reaches the fit, `writeFeatureDebugInfo_`
/// throws `ElementNotFound` inside the OpenMP region, and the executed tool
/// prints OpenMS's fatal-exception block and is killed by SIGABRT (shell
/// status 134) in all three repetitions. A safe port cannot end that way: it
/// exits 8 with the message TOPPBase gives the exception where it can catch
/// it. Everything the executed run wrote before it died is written the same:
/// the console lines, the seed map, and log.txt up to the last byte the file
/// buffer had flushed; no output and no later debug file.
#[test]
fn a_debug_run_that_reaches_the_fit_ends_where_the_release_build_terminates() {
    let dir = Workdir::new();
    let (_, _, args) = debug_args(&dir, &[]);
    let outcome = run_args(&dir, &args);
    outcome.assert_exit(ExitCode::UnknownError);
    outcome.assert_err_contains(
        "Error: Unexpected internal error (the element 'debug:pseudo_rt_shift' could not be found)",
    );
    let executed = fs::read_to_string(instrumentation("tool_b1_stdout.txt")).unwrap();
    assert!(executed.contains("FATAL: uncaught exception!"));
    assert!(
        executed.contains("error message: the element 'debug:pseudo_rt_shift' could not be found")
    );
    assert_eq!(port_block(&outcome), executed_block("tool_b1_stdout.txt"));
    let debug = dir.path().join("debug");
    assert_executed_file(&debug.join("log.txt"), "b1", "log.txt");
    assert_same_features(
        &debug.join("seeds_2.featureXML"),
        &instrumentation("ffc1_seed_map.featureXML.gz"),
    );
    assert!(debug.join("features").is_dir());
    assert!(!debug.join("abort_reasons.featureXML").exists());
    assert!(!debug.join("input.mzML").exists());
    assert!(!dir.path().join("out.featureXML").exists());
}

/// Executed tool case avg0 (`../oracle/ffap-complete-fix3`): FFC_1 with
/// `reported_mz` average and the trace, seed, feature and isotope-fit
/// thresholds at 0, without `write_debug`. A seed whose best isotope pattern
/// stayed empty reaches `extendMassTraces_`, which reads the pattern's first
/// entry: the executed tool dies with SIGSEGV (shell status 139) in both
/// repetitions, with nothing on stderr and no output. The port refuses at
/// that seed (lead decision D1) and reports it as TOPPBase reports the
/// algorithm's other errors. Its console block starts with everything the
/// executed process wrote; after it, the port also prints the `std::cout`
/// lines of `run_` that the executed process still held in its buffer when it
/// died, the per-charge seed counts (TOPP native difference: a crash loses
/// unflushed output).
#[test]
fn a_run_that_reaches_an_empty_best_pattern_is_refused_where_the_release_build_crashes() {
    let dir = Workdir::new();
    let args: Vec<String> = [
        "-test",
        "-ini",
        &text(ffc1_ini()),
        "-in",
        &text(ffc1_input()),
        "-out",
        &dir.file("out.featureXML"),
        "-algorithm:feature:reported_mz",
        "average",
        "-algorithm:feature:min_trace_score",
        "0",
        "-algorithm:seed:min_score",
        "0",
        "-algorithm:feature:min_score",
        "0",
        "-algorithm:feature:min_isotope_fit",
        "0",
    ]
    .iter()
    .map(|arg| (*arg).to_owned())
    .collect();
    let outcome = run_args(&dir, &args);
    outcome.assert_exit(ExitCode::UnknownError);
    outcome.assert_err_contains(
        "Error: Unexpected internal error (FeatureFinderAlgorithmPicked seed extension: the \
         isotope pattern matched no peak; the source reads its first entry here)",
    );
    let executed = executed_block("tool_avg0_stdout.txt");
    assert_eq!(
        executed,
        vec![FeatureFinderCentroided::NO_FAIMS_MESSAGE.to_owned()]
    );
    let port = port_block(&outcome);
    assert_eq!(port[..executed.len()], executed[..]);
    assert!(!port[executed.len()..].is_empty());
    assert!(
        port[executed.len()..]
            .iter()
            .all(|line| line.starts_with("Found ") && line.contains(" seeds for charge ")),
        "{port:?}"
    );
    assert!(!dir.path().join("out.featureXML").exists());
    assert!(!dir.path().join("debug").exists());
}

/// Executed cases c1 and c2 ran the C++ tool at `-threads 4`: its debug log
/// then differs from run to run (c2's three logs are pairwise different,
/// c1's too), because `abort_` and the log writes race; that output is
/// undefined and not compared. The port writes the single-thread files at
/// every thread count.
#[test]
fn debug_files_do_not_depend_on_the_thread_count() {
    let mut reference: Option<Vec<(String, Vec<u8>)>> = None;
    for threads in ["1", "4"] {
        let dir = Workdir::new();
        let (_, _, args) = debug_args(
            &dir,
            &[
                "-algorithm:feature:min_isotope_fit",
                "1.0",
                "-threads",
                threads,
            ],
        );
        run_args(&dir, &args).assert_exit(ExitCode::ExecutionOk);
        let mut files = Vec::new();
        for name in [
            "debug/log.txt",
            "debug/seeds_2.featureXML",
            "debug/abort_reasons.featureXML",
            "debug/input.mzML",
            "out.featureXML",
        ] {
            files.push((name.to_owned(), fs::read(dir.path().join(name)).unwrap()));
        }
        match &reference {
            None => reference = Some(files),
            Some(expected) => assert_eq!(&files, expected),
        }
    }
}

/// Executed case c3: case a1 at `-threads 4`. No seed is found, so nothing is
/// written from inside the parallel region and the run is defined: both
/// executed repetitions wrote exactly a1's files (`make_fixtures.py` checks
/// them byte for byte against a1's). The port at four threads writes them
/// too; `log.txt` is compared with c3's own digest.
#[test]
fn a_debug_run_without_seeds_at_four_threads_writes_the_executed_files() {
    let dir = Workdir::new();
    let (_, _, args) = debug_args(
        &dir,
        &["-algorithm:mass_trace:min_spectra", "1", "-threads", "4"],
    );
    let outcome = run_args(&dir, &args);
    outcome.assert_exit(ExitCode::ExecutionOk);
    assert_eq!(port_block(&outcome), executed_block("tool_a1_stdout.txt"));
    let debug = dir.path().join("debug");
    assert_eq!(fs::read_dir(debug.join("features")).unwrap().count(), 0);
    assert_executed_file(&debug.join("log.txt"), "c3", "log.txt");
    assert_same_features(
        &debug.join("seeds_2.featureXML"),
        &instrumentation("empty_seed_map.featureXML"),
    );
    let abort_path = debug.join("abort_reasons.featureXML");
    assert_same_features(&abort_path, &instrumentation("empty_abort_map.featureXML"));
    assert_eq!(decoded_features(&abort_path).unique_id, FIRST_TEST_MODE_ID);
    assert_same_debug_input(
        &debug.join("input.mzML"),
        &instrumentation("a1_input.mzML.gz"),
    );
    let output = PathBuf::from(dir.file("out.featureXML"));
    assert_same_features(&output, &instrumentation("tool_a1_out.featureXML"));
    assert_eq!(decoded_features(&output).unique_id, OUTPUT_ID_AFTER_DEBUG);
}

/// FeatureFinderCentroided_1's input with the last m/z of its last spectrum
/// replaced by `value` (oracle inputs `huge_mz_1e19.mzML` and
/// `huge_mz_2e18.mzML` of `../oracle/ffap-complete-fix2/make_inputs.py`; the
/// SHA-1 is that of the file the C++ tool read).
fn derive_last_mz(source: &[u8], value: f64, expected_sha1: &str) -> Vec<u8> {
    let engine = base64::engine::general_purpose::STANDARD;
    let all = lines(source);
    let mut pending = false;
    let mut last = None;
    for (index, line) in all.iter().enumerate() {
        if find(line, br#"name="m/z array""#).is_some() {
            pending = true;
        }
        if pending && trimmed(line).starts_with(b"<binary>") {
            last = Some(index);
            pending = false;
        }
    }
    let last = last.unwrap();
    let mut out: Vec<Vec<u8>> = all.iter().map(|line| line.to_vec()).collect();
    let line = all[last];
    let prefix = &line[..find(line, b"<binary>").unwrap()];
    let old = &trimmed(line)[b"<binary>".len()..trimmed(line).len() - b"</binary>".len()];
    let mut raw = engine.decode(old).unwrap();
    assert_eq!(raw.len(), 24 * 8);
    let at = raw.len() - 8;
    let previous = f64::from_le_bytes(raw[at..].try_into().unwrap());
    assert!(previous < value);
    raw[at..].copy_from_slice(&value.to_le_bytes());
    let payload = engine.encode(&raw);
    assert_eq!(payload.len(), old.len());
    let mut replaced = prefix.to_vec();
    replaced.extend_from_slice(b"<binary>");
    replaced.extend_from_slice(payload.as_bytes());
    replaced.extend_from_slice(b"</binary>");
    out[last] = replaced;
    let derived = join(out);
    assert_digest(&derived, expected_sha1, "huge m/z");
    derived
}

/// The executed rows of `length_error.tsv` for one input and kind of the tool.
fn length_error_tool(input: &str, kind: &str) -> Vec<String> {
    fs::read_to_string(instrumentation("length_error.tsv"))
        .unwrap()
        .lines()
        .skip(2)
        .filter_map(|line| {
            let fields: Vec<&str> = line.splitn(4, '\t').collect();
            (fields[0] == input && fields[1] == "tool" && fields[2] == kind)
                .then(|| fields[3].to_owned())
        })
        .collect()
}

/// Executed cases `tool_1e19` and `tool_2e18` (`../oracle/ffap-complete-fix2`,
/// two identical runs each): a debug run on an input whose maximum m/z asks
/// step 2.5 for more isotope windows than it can allocate.
///
/// At m/z `1e19` (`2e17 + 1` windows at the INI's charge 2 and width 100,
/// above `vector::max_size()`) the source's `resize` throws
/// `std::length_error`, which is no OpenMS exception:
/// TOPPBase's outer handler prints `Unable to initialize or run
/// FeatureFinderCentroided: vector::_M_default_append` and returns 12
/// (`INTERNAL_ERROR`), as this port does. The run has created
/// `debug/features` and written the first log line, which the unwinding
/// flushes: 40 bytes. No output is written.
///
/// At m/z `2e18` (`4e16 + 1` windows, below the bound) the executed
/// allocation throws `std::bad_alloc`, again exit 12. This port refuses the
/// count with its native window ceiling and exits 8 with that message (a
/// recorded native difference); the debug side effects are the executed ones.
#[test]
fn a_debug_run_beyond_the_isotope_window_limit_exits_as_the_release_build() {
    let source = fs::read(ffc1_input()).unwrap();
    let log = fs::read(instrumentation("length_error_log.txt")).unwrap();
    for (tag, value, sha1) in [
        ("1e19", 1e19, "76954c288ddadebd6d4d9fd579ba3be764c85e26"),
        ("2e18", 2e18, "fee5f77d5582d01388065004ab621116840a8be8"),
    ] {
        let dir = Workdir::new();
        let input = dir.put(
            &format!("huge_mz_{tag}.mzML"),
            &derive_last_mz(&source, value, sha1),
        );
        let ini = text(ffc1_ini());
        let out = dir.file("out.featureXML");
        let outcome = run_in(
            &dir,
            &[
                "-test",
                "-ini",
                &ini,
                "-in",
                &input,
                "-out",
                &out,
                "-algorithm:write_debug",
                "-algorithm:feature:min_isotope_fit",
                "1.0",
            ],
        );
        let executed_status = length_error_tool(tag, "status");
        let executed_stderr = length_error_tool(tag, "stderr");
        assert_eq!(executed_status, ["12"]);
        if tag == "1e19" {
            outcome.assert_exit(ExitCode::InternalError);
            assert_eq!(
                executed_stderr,
                ["Unable to initialize or run FeatureFinderCentroided: vector::_M_default_append"]
            );
            assert_eq!(outcome.err.lines().collect::<Vec<_>>(), executed_stderr);
        } else {
            assert_eq!(
                executed_stderr,
                ["Unable to initialize or run FeatureFinderCentroided: std::bad_alloc"]
            );
            outcome.assert_exit(ExitCode::UnknownError);
            outcome.assert_err_contains("Error: Unexpected internal error (");
            outcome.assert_err_contains("exceed the limit");
        }
        assert_eq!(
            port_block(&outcome),
            length_error_tool(tag, "stdout_block"),
            "{tag}"
        );
        let debug = dir.path().join("debug");
        assert!(debug.join("features").is_dir());
        assert_eq!(fs::read_dir(debug.join("features")).unwrap().count(), 0);
        assert_eq!(fs::read(debug.join("log.txt")).unwrap(), log);
        assert_eq!(length_error_tool(tag, "log_bytes"), [log.len().to_string()]);
        let mut entries: Vec<String> = fs::read_dir(&debug)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        entries.sort();
        assert_eq!(entries, ["features", "log.txt"]);
        let tree = length_error_tool(tag, "tree");
        assert!(tree.contains(&"./debug/log.txt".to_owned()));
        assert!(!tree.iter().any(|path| path.contains("out.featureXML")));
        assert!(!Path::new(&out).exists());
    }
}
