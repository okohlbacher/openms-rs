// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The instance, its side channels and its debug mode
//! (`FEATUREFINDER/FeatureFinderAlgorithmPicked.h`): reuse of one object, a
//! caller's non-empty feature map, `aborts_`, `abort_reasons_`, the
//! `DefaultParamHandler` and `ProgressLogger` bases, `write_debug` with its
//! log, seed maps, abort map and debug input, and `writeFeatureDebugInfo_`.
//!
//! Evidence (`tests/data/feature_finder_picked_instrumentation_provenance.json`):
//!
//! - tier 1, executed C++: the Linux x86_64 **Release** build
//!   `openms4-release-bc9cc12-c19e494-174b576`, through the oracle driver
//!   `../oracle/ffap-instr-completion/drivers/ffap_instr_driver.cpp`, which
//!   subclasses the algorithm to read its protected members and prints every
//!   number as its IEEE-754 bit pattern, and through the Release
//!   `FeatureFinderCentroided`. Every case ran twice with `OMP_NUM_THREADS=1`
//!   and reproduced (apart from timing text and the random map id the driver
//!   draws outside test mode). Text outputs are compared byte for byte,
//!   through their SHA-1 where the file is large; featureXML and mzML outputs
//!   by decoded content (decision D6);
//! - tier 4: the undefined cases, which have no executed answer.
//!
//! Fitted values (retention time, `score_fit`, `score_correlation`, the
//! `EGH_*` parameters) are compared bit for bit on every platform: both trace
//! fitters call the reference build's glibc `exp` and `log`, ported (lead
//! decision D10). The one exception is the intensity of the asymmetric case on
//! a host without the GNU C Library, whose `atan` is not the reference's
//! (`fit_bound`). Everything else is compared exactly.

#![cfg(all(feature = "mzml", feature = "paramxml", feature = "featurexml"))]

#[path = "support/decoded_compare.rs"]
mod decoded;
#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use openms::analysis::feature_finder_picked::algorithm::{
    Limits, Options, PseudoRtShiftKey, RejectedParameters, run_with_options,
};
use openms::analysis::feature_finder_picked::debug::{
    DebugOutput, FeatureDebugInput, HEAP_ADDRESS_END, PseudoRtShift, ReportLine, TerminationKind,
    TerminationPoint, write_feature_debug_info,
};
use openms::analysis::feature_finder_picked::instance::FeatureFinderAlgorithmPicked;
use openms::analysis::feature_finder_picked::source_sort::source_sort_permutation;
use openms::concept::parallel::Threads;
use openms::concept::progress_logger::{
    CommandProgressLogger, ProgressBackend, ProgressLogType, ProgressLogger, ProgressNesting,
    ProgressTime,
};
use openms::format::{FileHandler, FileType, PeakFileOptions, featurexml, paramxml};
use openms::kernel::{ConvexHull2D, Feature, FeatureMap, MSExperiment, NumericRange, Point2D};
use openms::metadata::{MetaValue, MetaValueData};
use openms::param::{Param, ParamValue};

// ---------------------------------------------------------------------------
// Inputs, fixtures and the dump format of the oracle driver
// ---------------------------------------------------------------------------

fn repository(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

fn fixture_path(name: &str) -> PathBuf {
    repository("tests/data/feature_finder_picked_instrumentation").join(name)
}

/// A fixture as text; a `.gz` fixture is decompressed.
fn fixture(name: &str) -> String {
    let bytes = std::fs::read(fixture_path(name)).unwrap();
    if name.ends_with(".gz") {
        let mut text = String::new();
        flate2::read::GzDecoder::new(bytes.as_slice())
            .read_to_string(&mut text)
            .unwrap();
        text
    } else {
        String::from_utf8(bytes).unwrap()
    }
}

/// An input as the driver loads it: MS1, intensities in `[0, f64::MAX)`.
fn load_input(path: PathBuf) -> MSExperiment {
    let mut options = PeakFileOptions::default();
    options.add_ms_level(1).unwrap();
    options.set_intensity_range(NumericRange {
        min: 0.0,
        max: f64::MAX,
    });
    FileHandler::load_experiment_with_options(path, &[FileType::MzMl], &options).unwrap()
}

/// FeatureFinderCentroided_1_input.mzML (test-data 0cb15f2, sha256 a3dfae63...).
fn ffc1_input() -> MSExperiment {
    load_input(repository(
        "tests/data/mzml_mobility/FeatureFinderCentroided_1_input.mzML",
    ))
}

/// FileConverter_31_output.mzML: four MS1 scans, fewer than the default seed
/// window needs.
fn short_input() -> MSExperiment {
    load_input(repository(
        "tests/data/topp_feature_finder_centroided/FileConverter_31_output.mzML",
    ))
}

/// The `algorithm` section of the FeatureFinderCentroided_1 INI.
fn ffc1_parameters() -> Param {
    paramxml::load(repository(
        "tests/data/feature_finder_picked/FeatureFinderCentroided_1_parameters.ini",
    ))
    .unwrap()
    .copy("FeatureFinderCentroided:1:algorithm:", true)
    .unwrap()
}

fn set(parameters: &mut Param, key: &str, value: ParamValue) {
    parameters
        .set_value(key, value, "", &[] as &[String])
        .unwrap();
}

fn hex64(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

fn hex32(value: f32) -> String {
    format!("{:08x}", value.to_bits())
}

/// The driver's `escape`.
fn escape(text: &str) -> String {
    text.replace('\\', "\\\\")
        .replace('\n', "\\n")
        .replace('\t', "\\t")
}

fn meta_text(value: &MetaValue) -> String {
    match value.data() {
        MetaValueData::String(text) => format!("s {}", escape(text)),
        MetaValueData::Integer(value) => format!("i {value}"),
        MetaValueData::Float(value) => format!("d {}", hex64(*value)),
        MetaValueData::StringList(values) => {
            let mut out = "sl".to_owned();
            for value in values {
                out.push(' ');
                out.push_str(&escape(value));
            }
            out
        }
        MetaValueData::IntegerList(values) => {
            let mut out = "il".to_owned();
            for value in values {
                out.push_str(&format!(" {value}"));
            }
            out
        }
        MetaValueData::FloatList(values) => {
            let mut out = "dl".to_owned();
            for value in values {
                out.push_str(&format!(" {}", hex64(*value)));
            }
            out
        }
        MetaValueData::Empty => "e".to_owned(),
    }
}

/// The driver's `dumpFeature`.
fn dump_feature(out: &mut String, feature: &Feature, depth: usize, index: usize) {
    let indent = " ".repeat(2 * depth);
    out.push_str(&format!(
        "{indent}F {depth} {index} rt={} mz={} int={} charge={} q={} qrt={} qmz={} w={} uid={} hulls={} subs={}\n",
        hex64(feature.rt),
        hex64(feature.mz),
        hex32(feature.intensity),
        feature.charge,
        hex32(feature.quality),
        hex32(feature.quality_rt),
        hex32(feature.quality_mz),
        hex32(feature.width),
        feature.unique_id,
        feature.convex_hulls.len(),
        feature.subordinates.len()
    ));
    for (key, value) in &feature.metadata {
        out.push_str(&format!("{indent}  M {key} {}\n", meta_text(value)));
    }
    for (h, hull) in feature.convex_hulls.iter().enumerate() {
        let points = hull.hull_points();
        out.push_str(&format!("{indent}  H {h} {}", points.len()));
        for point in points {
            out.push_str(&format!(" {},{}", hex64(point.rt), hex64(point.mz)));
        }
        out.push('\n');
    }
    for (s, sub) in feature.subordinates.iter().enumerate() {
        dump_feature(out, sub, depth + 1, s);
    }
}

/// The driver's `dumpMap`: the map and the instance's `aborts_`.
fn dump_map(map: &FeatureMap, aborts: &BTreeMap<String, u32>) -> String {
    let mut out = format!("MAP size={} uid={}\n", map.len(), map.unique_id);
    for (index, feature) in map.features.iter().enumerate() {
        dump_feature(&mut out, feature, 0, index);
    }
    for (reason, count) in aborts {
        out.push_str(&format!("ABORT {count} {}\n", escape(reason)));
    }
    out
}

/// The driver's `dumpParam`.
fn dump_param(parameters: &Param) -> String {
    let mut out = String::new();
    let mut sections = std::collections::BTreeSet::new();
    for item in parameters.iter().unwrap() {
        let key = item.key.clone();
        for (position, _) in key.match_indices(':') {
            sections.insert(key[..position].to_owned());
        }
        let entry = item.entry;
        let value = match &entry.value {
            ParamValue::String(value) => format!("s {}", escape(value)),
            ParamValue::Integer(value) => format!("i {value}"),
            ParamValue::Float(value) => format!("d {}", hex64(*value)),
            ParamValue::StringList(values) => {
                let mut out = "sl".to_owned();
                for value in values {
                    out.push(' ');
                    out.push_str(&escape(value));
                }
                out
            }
            ParamValue::IntegerList(values) => {
                let mut out = "il".to_owned();
                for value in values {
                    out.push_str(&format!(" {value}"));
                }
                out
            }
            ParamValue::FloatList(values) => {
                let mut out = "dl".to_owned();
                for value in values {
                    out.push_str(&format!(" {}", hex64(*value)));
                }
                out
            }
            ParamValue::Empty => "e".to_owned(),
        };
        out.push_str(&format!("P {key} {value}\n"));
        out.push_str("  tags");
        for tag in &entry.tags {
            out.push(' ');
            out.push_str(tag);
        }
        out.push('\n');
        out.push_str(&format!(
            "  int {} {} float {} {}\n",
            entry.min_int,
            entry.max_int,
            hex64(entry.min_float),
            hex64(entry.max_float)
        ));
        out.push_str("  valid");
        for valid in &entry.valid_strings {
            out.push(' ');
            out.push_str(valid);
        }
        out.push('\n');
        out.push_str(&format!("  desc {}\n", escape(&entry.description)));
    }
    for section in sections {
        out.push_str(&format!(
            "S {section} {}\n",
            escape(parameters.section_description(&section).unwrap_or(""))
        ));
    }
    out
}

/// Fitted values admit no difference: both trace fitters call the reference
/// build's glibc `exp` and `log`, ported (lead decision D10), so every fitted
/// value is the Linux x86_64 Release build's on every platform, NaN bits
/// included.
const BIT_FOR_BIT: f64 = 0.0;

/// The relative bound of a differing fitted token of `case`: [`BIT_FOR_BIT`],
/// except the intensity (`int=`) of the asymmetric case on a host other than
/// x86_64 Linux with the GNU C Library, whose EGH area calls the `libm`
/// crate's `atan` instead of the reference's (D10's fallback;
/// [`EGH_ATAN_GAP`]).
fn fit_bound(case: &str, token: &str) -> f64 {
    if case == "egh" && token.starts_with("int=") {
        EGH_ATAN_GAP
    } else {
        BIT_FOR_BIT
    }
}

/// The EGH intensity bound: exact with the host's `atan` on x86_64 Linux with
/// the GNU C Library, which is the reference's `__atan_fma` where the host
/// library selects it (GNU C Library 2.39 on a CPU with FMA, as on the
/// reference node and the gate hosts; other versions and CPUs are not
/// measured).
#[cfg(all(target_os = "linux", target_env = "gnu", target_arch = "x86_64"))]
const EGH_ATAN_GAP: f64 = BIT_FOR_BIT;

/// The EGH intensity bound elsewhere: the largest relative departure measured
/// over these cases on macOS arm64 with the `libm` crate's `atan`, a measured
/// maximum, not a guarantee (the crate's `atan` itself differs from the
/// reference for 1.6 to 6.2 % of arguments in the fits' range; the `float`
/// intensity hides nearly all of it).
#[cfg(not(all(target_os = "linux", target_env = "gnu", target_arch = "x86_64")))]
const EGH_ATAN_GAP: f64 = 0.0;

/// Whether a dump token holds a fitted value, compared within [`fit_bound`];
/// everything else is compared exactly.
fn fitted_field(token: &str) -> bool {
    ["rt=", "q=", "w=", "int="]
        .iter()
        .any(|prefix| token.starts_with(prefix))
}

fn fitted_meta(key: &str) -> bool {
    matches!(
        key,
        "score_fit" | "score_correlation" | "FWHM" | "EGH_tau" | "EGH_height" | "EGH_sigma"
    )
}

fn value_of(token: &str) -> f64 {
    let hex = token.split('=').nth(1).unwrap_or(token);
    if hex.len() == 8 {
        f64::from(f32::from_bits(u32::from_str_radix(hex, 16).unwrap()))
    } else {
        f64::from_bits(u64::from_str_radix(hex, 16).unwrap())
    }
}

/// Whether two differing fitted tokens are within `bound`, and their relative
/// departure. A zero bound admits nothing, NaN bits included; elsewhere two
/// NaNs of any bits are alike, since only the reference platform reproduces
/// x86_64's NaN bits inside the fit.
fn close(expected: f64, actual: f64, bound: f64) -> (bool, f64) {
    if expected.is_nan() && actual.is_nan() {
        return (bound != 0.0, 0.0);
    }
    let departure = (expected - actual).abs() / expected.abs().max(actual.abs());
    (bound != 0.0 && departure <= bound, departure)
}

/// Compare two dumps: identical line structure, exact tokens, and fitted
/// values within [`fit_bound`] of `case`. Returns the number of fitted values
/// that were not bit-identical.
fn assert_dumps_match(expected: &str, actual: &str, case: &str) -> usize {
    let mut largest = 0.0f64;
    let expected_lines: Vec<&str> = expected.lines().collect();
    let actual_lines: Vec<&str> = actual.lines().collect();
    assert_eq!(
        expected_lines.len(),
        actual_lines.len(),
        "{case}: line count\nexpected:\n{expected}\nactual:\n{actual}"
    );
    let mut inexact = 0usize;
    for (number, (e, a)) in expected_lines.iter().zip(&actual_lines).enumerate() {
        if e == a {
            continue;
        }
        let et: Vec<&str> = e.split_whitespace().collect();
        let at: Vec<&str> = a.split_whitespace().collect();
        assert_eq!(et.len(), at.len(), "{case}: line {number}\n{e}\n{a}");
        let meta_key = (et.first() == Some(&"M")).then(|| et[1]);
        for (index, (x, y)) in et.iter().zip(&at).enumerate() {
            if x == y {
                continue;
            }
            let fitted = match meta_key {
                Some(key) => fitted_meta(key) && index == 3,
                None => et.first() == Some(&"F") && fitted_field(x),
            };
            let bound = fit_bound(case, x);
            let (within, departure) = close(value_of(x), value_of(y), bound);
            assert!(
                fitted && within,
                "{case}: line {number} token {index}: expected {x}, got {y} (departure \
                 {departure:e}, bound {bound:e})\n{e}\n{a}"
            );
            largest = largest.max(departure);
            inexact += 1;
        }
    }
    if inexact > 0 {
        eprintln!(
            "{case}: {inexact} fitted values within their bound, largest departure {largest:e}"
        );
    }
    inexact
}

// ---------------------------------------------------------------------------
// The C++ Release build's std::sort
// ---------------------------------------------------------------------------

/// The driver's key generator (Knuth's MMIX LCG, upper 31 bits).
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }
}

fn permutation_text(name: &str, order: &[usize]) -> String {
    let mut out = name.to_owned();
    for index in order {
        out.push_str(&format!(" {index}"));
    }
    out
}

/// `FeatureMap::sortByMZ` and `sortByIntensity(true)` of the Release
/// `libOpenMS.so` on 130 generated inputs with few distinct keys and on the
/// eight key sequences of `sort_keys.txt` (three built by McIlroy's adversary
/// against the emulation so that the depth budget runs out and the heapsort
/// takes over, three with NaN keys, which are no strict weak ordering but kept
/// the executed sort inside the range, and two with signed zeros): every
/// permutation, ties included, is the one the executed `std::sort` produced.
#[test]
fn sorts_leave_ties_where_the_release_build_leaves_them() {
    let expected = fixture("sort.txt.gz");
    let mut actual = Vec::new();
    let mut lcg = Lcg(12345);
    for n in [1usize, 2, 3, 5, 16, 17, 18, 31, 64, 100, 257, 1000, 2048] {
        for k in [1u64, 2, 3, 7, 1_000_000] {
            let mut mz = Vec::with_capacity(n);
            let mut intensity = Vec::with_capacity(n);
            for _ in 0..n {
                mz.push(100.0 + (lcg.next() % k) as f64);
                intensity.push((lcg.next() % k) as f32);
            }
            let by_mz = source_sort_permutation(n, |a, b| mz[a] < mz[b]).unwrap();
            actual.push(permutation_text(&format!("mz n={n} k={k}"), &by_mz));
            let by_int = source_sort_permutation(n, |a, b| intensity[b] < intensity[a]).unwrap();
            actual.push(permutation_text(&format!("int n={n} k={k}"), &by_int));
        }
    }
    for (case, line) in fixture("sort_keys.txt")
        .lines()
        .filter(|line| !line.is_empty())
        .enumerate()
    {
        let keys: Vec<f64> = line
            .split(' ')
            .map(|key| key.parse::<f64>().unwrap())
            .collect();
        let intensities: Vec<f32> = keys.iter().map(|&key| key as f32).collect();
        let by_mz = source_sort_permutation(keys.len(), |a, b| keys[a] < keys[b]).unwrap();
        actual.push(permutation_text(&format!("file-mz {case}"), &by_mz));
        let by_int =
            source_sort_permutation(keys.len(), |a, b| intensities[b] < intensities[a]).unwrap();
        actual.push(permutation_text(&format!("file-int {case}"), &by_int));
    }
    let expected: Vec<&str> = expected.lines().collect();
    assert_eq!(expected.len(), actual.len());
    for (e, a) in expected.iter().zip(&actual) {
        assert_eq!(*e, a.as_str());
    }
}

// ---------------------------------------------------------------------------
// The DefaultParamHandler surface
// ---------------------------------------------------------------------------

/// `getName`, `getDefaults`, `getDefaultParameters` and `getParameters` after
/// construction, with an unknown key, after a run and after an empty-input
/// run, against the executed driver (`params_*.txt`).
#[test]
fn the_parameter_surface_matches_the_release_build() {
    let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
    let info = fixture("params_info.txt");
    assert!(info.starts_with(&format!("name {}\n", algorithm.name())));
    let defaults = fixture("params_defaults.txt");
    assert_eq!(dump_param(algorithm.defaults()), defaults);
    assert_eq!(dump_param(&algorithm.default_parameters()), defaults);
    assert_eq!(dump_param(algorithm.parameters()), defaults);

    // An unknown key is kept, with the source's warning.
    let mut extra = ffc1_parameters();
    set(
        &mut extra,
        "debug:pseudo_rt_shift",
        ParamValue::Float(250.0),
    );
    let warnings = algorithm.set_parameters(&extra).unwrap();
    assert_eq!(
        warnings,
        [
            "Warning: FeatureFinderAlgorithmPicked received the unknown parameter 'debug:pseudo_rt_shift'!"
        ]
    );
    let with_unknown = fixture("params_with_unknown_key.txt");
    assert_eq!(dump_param(algorithm.parameters()), with_unknown);

    // A restriction violation and a type violation are refused. The executed
    // messages are `InvalidParameter` texts; the port reports the same
    // violation. The source has already assigned the refused set, so the
    // executed `getParameters()` then shows it (`bins 0`), and so does the
    // port by default (`RejectedParameters::Shown`).
    let mut bad = Param::new();
    set(&mut bad, "intensity:bins", ParamValue::Integer(0));
    let error = algorithm.set_parameters(&bad).unwrap_err().to_string();
    assert!(error.contains("bins"), "{error}");
    assert!(info.contains(
        "bad_bins InvalidParameter FeatureFinderAlgorithmPicked: Invalid integer parameter value '0' for parameter 'bins' given!"
    ));
    let failed = fixture("params_after_failed_set.txt");
    assert!(failed.contains("P intensity:bins i 0\n"));
    assert_eq!(dump_param(algorithm.parameters()), failed);
    // The settings a run would use are still the accepted set's.
    assert_eq!(algorithm.settings().intensity_bins, 1);
    let mut bad_type = Param::new();
    set(
        &mut bad_type,
        "intensity:bins",
        ParamValue::String("ten".into()),
    );
    assert!(algorithm.set_parameters(&bad_type).is_err());
    assert!(info.contains("bad_type InvalidParameter"));

    // After a run the run's parameters, merged with the defaults.
    let mut features = FeatureMap::new();
    algorithm
        .run(
            ffc1_input(),
            &mut features,
            &ffc1_parameters(),
            &FeatureMap::new(),
        )
        .unwrap();
    assert!(info.contains("features_after_run 8\n"));
    assert_eq!(features.len(), 8);
    let after_run = fixture("params_after_run.txt");
    assert_eq!(dump_param(algorithm.parameters()), after_run);

    // An empty input clears the caller's map completely and leaves the
    // parameters alone.
    let mut other = Param::new();
    set(&mut other, "intensity:bins", ParamValue::Integer(7));
    let mut cleared = prefilled();
    algorithm
        .run(
            MSExperiment::new(),
            &mut cleared,
            &other,
            &FeatureMap::new(),
        )
        .unwrap();
    assert_eq!(dump_param(algorithm.parameters()), after_run);
    assert_eq!(
        dump_map(&cleared, algorithm.aborts()),
        fixture("params_empty_run_map.txt")
    );
    assert!(info.contains("empty_run_map_uid 0 meta_empty 1\n"));
    assert!(cleared.metadata.is_empty());
}

/// The sections of `rejected_stdout.txt` (`== <name>` headers) and of
/// `rejected_stderr.txt`, by name.
fn rejected_sections(text: &str) -> BTreeMap<String, String> {
    let mut sections = BTreeMap::new();
    let mut current: Option<(String, String)> = None;
    for line in text.split_inclusive('\n') {
        if let Some(name) = line.strip_prefix("== ") {
            if let Some((name, body)) = current.take() {
                sections.insert(name, body);
            }
            current = Some((name.trim_end().to_owned(), String::new()));
        } else if let Some((_, body)) = current.as_mut() {
            body.push_str(line);
        }
    }
    if let Some((name, body)) = current {
        sections.insert(name, body);
    }
    sections
}

/// The driver case `ffap_progress_driver rejected`: an accepted set with
/// unknown keys (the warnings in `checkDefaults` order), a refused set with
/// unknown keys before and after the refused entry (only the ones before are
/// logged, then `InvalidParameter`, and `getParameters()` shows the refused
/// set), and a run with the refused set (the same warnings, which the log
/// stream folds into its `occurred 2 times` lines, the same exception, and
/// the caller's map untouched).
#[test]
fn a_refused_parameter_set_matches_the_release_build() {
    let out = rejected_sections(&fixture("rejected_stdout.txt.gz"));
    let err = rejected_sections(&fixture("rejected_stderr.txt"));
    let lines = |text: &str| -> Vec<String> { text.lines().map(str::to_owned).collect() };
    let mut good = ffc1_parameters();
    set(&mut good, "zz:late_unknown", ParamValue::Integer(1));
    set(&mut good, "a_unknown", ParamValue::Float(1.5));
    set(
        &mut good,
        "intensity:c_unknown",
        ParamValue::String("x".into()),
    );
    let refused = || {
        let mut p = Param::new();
        set(&mut p, "a_unknown", ParamValue::Integer(1));
        set(&mut p, "intensity:bins", ParamValue::Integer(0));
        set(&mut p, "intensity:b_unknown", ParamValue::Integer(2));
        set(&mut p, "z_unknown", ParamValue::Integer(3));
        p
    };

    let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
    let warnings = algorithm.set_parameters(&good).unwrap();
    assert_eq!(warnings, lines(&err["set good"]));
    assert_eq!(dump_param(algorithm.parameters()), out["params after good"]);

    let mut warnings = Vec::new();
    let error = algorithm
        .set_parameters_logged(&refused(), &mut warnings)
        .unwrap_err();
    assert!(matches!(&error, openms::Error::InvalidValue(m) if m.contains("bins")));
    assert!(out["set rejected"].starts_with("exception InvalidParameter "));
    assert_eq!(warnings, lines(&err["set rejected"]));
    assert_eq!(
        dump_param(algorithm.parameters()),
        out["params after rejected"]
    );

    let mut features = FeatureMap::from_features(vec![Feature::new(0.0, 500.0, 0.0)]);
    let before = features.clone();
    let result = algorithm.run(ffc1_input(), &mut features, &refused(), &FeatureMap::new());
    assert!(matches!(result, Err(openms::Error::InvalidValue(_))));
    assert!(out["run rejected"].starts_with("exception InvalidParameter "));
    assert!(out["run rejected"].ends_with("features 1\n"));
    assert_eq!(features, before);
    let warned: Vec<String> = algorithm
        .report()
        .iter()
        .filter_map(|line| match line {
            ReportLine::Warn(text) => Some(format!("<{text}> occurred 2 times")),
            _ => None,
        })
        .collect();
    assert_eq!(warned, lines(&err["end"]));
    assert_eq!(dump_param(algorithm.parameters()), out["params after run"]);

    // The native option keeps the accepted set.
    let mut atomic = FeatureFinderAlgorithmPicked::with_options(Options {
        rejected_parameters: RejectedParameters::Discarded,
        ..Options::default()
    })
    .unwrap();
    atomic.set_parameters(&good).unwrap();
    assert!(atomic.set_parameters(&refused()).is_err());
    assert_eq!(dump_param(atomic.parameters()), out["params after good"]);
}

// ---------------------------------------------------------------------------
// Reusing an instance and extending a caller's map
// ---------------------------------------------------------------------------

/// A hull from `(rt, mz)` points, as `ConvexHull2D::addPoint` builds it.
fn hull(points: &[(f64, f64)]) -> ConvexHull2D {
    let points: Vec<Point2D> = points
        .iter()
        .map(|&(rt, mz)| Point2D::new(rt, mz))
        .collect();
    ConvexHull2D::from_points(&points).unwrap()
}

fn synthetic(rt: f64, mz: f64, intensity: f32, charge: i32, quality: f32, label: &str) -> Feature {
    let mut feature = Feature::new(rt, mz, intensity);
    feature.charge = charge;
    feature.quality = quality;
    feature
        .metadata
        .insert("label".into(), MetaValue::from(label.to_owned()));
    feature
}

/// The driver's `prefilled()` map, feature for feature.
fn prefilled() -> FeatureMap {
    let mut features = Vec::new();
    let mut f = synthetic(1720.0, 445.2, 5.0e6, 2, 0.9, "caller-0");
    f.convex_hulls.push(hull(&[
        (1650.0, 445.20),
        (1700.0, 445.18),
        (1790.0, 445.22),
    ]));
    f.convex_hulls
        .push(hull(&[(1650.0, 445.70), (1790.0, 445.72)]));
    features.push(f);
    features.push(synthetic(1500.0, 500.0, 1000.0, 1, 0.5, "caller-1-no-hull"));
    let mut f = synthetic(1600.0, 600.0, 2000.0, 3, 0.4, "caller-2-empty-hull");
    f.convex_hulls.push(ConvexHull2D::new());
    f.convex_hulls
        .push(hull(&[(1590.0, 600.0), (1610.0, 600.05)]));
    features.push(f);
    let mut f = synthetic(1600.0, 700.0, 0.0, 2, 0.8, "caller-3-zero");
    f.convex_hulls
        .push(hull(&[(1590.0, 700.0), (1610.0, 700.05)]));
    features.push(f);
    let mut f = synthetic(1.0e6, 800.0, 3000.0, 2, 0.8, "caller-4-late");
    f.metadata
        .insert("spectrum_index".into(), MetaValue::from(12345i64));
    f.convex_hulls
        .push(hull(&[(999_990.0, 800.0), (1_000_010.0, 800.05)]));
    features.push(f);
    let mut f = synthetic(1400.0, 900.0, 4000.0, 2, 0.7, "caller-5-sub");
    f.unique_id = 4242;
    f.convex_hulls
        .push(hull(&[(1390.0, 900.0), (1410.0, 900.05)]));
    f.subordinates
        .push(synthetic(1405.0, 900.01, 10.0, 2, 0.1, "caller-5-sub-sub"));
    features.push(f);
    let mut f = synthetic(1402.0, 900.0, 4000.0, 2, 0.7, "caller-6-tie");
    f.unique_id = 4242;
    f.convex_hulls
        .push(hull(&[(1392.0, 900.0), (1412.0, 900.05)]));
    features.push(f);
    let mut f = synthetic(1300.0, 950.0, 100.0, 4, 0.2, "caller-7-z4");
    f.convex_hulls
        .push(hull(&[(1290.0, 950.0), (1310.0, 950.02)]));
    features.push(f);
    let mut f = synthetic(1301.0, 950.01, 200_000.0, 2, 0.9, "caller-8-z2");
    f.convex_hulls
        .push(hull(&[(1291.0, 950.0), (1311.0, 950.02)]));
    features.push(f);
    let mut f = synthetic(1200.0, 975.0, 100.0, 3, 0.6, "caller-9-z3");
    f.convex_hulls
        .push(hull(&[(1190.0, 975.0), (1210.0, 975.02)]));
    features.push(f);
    let mut f = synthetic(1200.5, 975.0, 50.0, 2, 0.6, "caller-10-z2");
    f.convex_hulls
        .push(hull(&[(1190.0, 975.0), (1210.0, 975.02)]));
    features.push(f);
    let mut map = FeatureMap::from_features(features);
    map.metadata
        .insert("caller_map_meta".into(), MetaValue::from("kept".to_owned()));
    map.unique_id = 777;
    map
}

/// One instance, three runs: FeatureFinderCentroided_1 twice into the same
/// map (old and new features are resolved together, labels restart, the
/// abort counts accumulate), then charges 1 to 3 into a fresh map (the counts
/// keep accumulating) with the new parameters in `getParameters()`.
#[test]
fn a_reused_instance_matches_the_release_build() {
    let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
    let mut features = FeatureMap::new();
    let mut inexact = 0;
    algorithm
        .run(
            ffc1_input(),
            &mut features,
            &ffc1_parameters(),
            &FeatureMap::new(),
        )
        .unwrap();
    inexact += assert_dumps_match(
        &fixture("reuse_run1.txt.gz"),
        &dump_map(&features, algorithm.aborts()),
        "run 1",
    );
    algorithm
        .run(
            ffc1_input(),
            &mut features,
            &ffc1_parameters(),
            &FeatureMap::new(),
        )
        .unwrap();
    inexact += assert_dumps_match(
        &fixture("reuse_run2.txt.gz"),
        &dump_map(&features, algorithm.aborts()),
        "run 2",
    );
    // The second run works with the isotope windows the first one left
    // (`isotope_distributions_` is never cleared), finds other features and
    // aborts no seed: the executed count stays at the first run's one.
    assert!(
        fixture("reuse_run2.txt.gz")
            .ends_with("ABORT 1 Invalid fit: Fitted model is bigger than 'max_rt_span'\n")
    );
    assert_eq!(
        algorithm
            .aborts()
            .get("Invalid fit: Fitted model is bigger than 'max_rt_span'"),
        Some(&1)
    );
    assert!(algorithm.isotope_windows().is_some());
    let mut changed = ffc1_parameters();
    set(
        &mut changed,
        "isotopic_pattern:charge_low",
        ParamValue::Integer(1),
    );
    set(
        &mut changed,
        "isotopic_pattern:charge_high",
        ParamValue::Integer(3),
    );
    let mut fresh = FeatureMap::new();
    algorithm
        .run(ffc1_input(), &mut fresh, &changed, &FeatureMap::new())
        .unwrap();
    inexact += assert_dumps_match(
        &fixture("reuse_run3.txt.gz"),
        &dump_map(&fresh, algorithm.aborts()),
        "run 3",
    );
    assert_eq!(
        dump_param(algorithm.parameters()),
        fixture("reuse_run3_params.txt")
    );
    eprintln!("reuse: {inexact} fitted values not bit-identical");
}

/// A caller's map with features that overlap new ones, no hull, an empty
/// hull, zero intensity, an RT past the last scan, a subordinate, equal m/z,
/// and charge pairs for the multiple-of and quality rules: the resulting map
/// equals the executed one, for the FeatureFinderCentroided_1 parameters and
/// for the defaults (all four charges).
#[test]
fn a_caller_map_is_extended_as_the_release_build_extends_it() {
    let mut inexact = 0;
    for (parameters, expected, case) in [
        (
            ffc1_parameters(),
            "prefilled_run.txt.gz",
            "FFC_1 parameters",
        ),
        (Param::new(), "prefilled_defaults_run.txt.gz", "defaults"),
    ] {
        let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
        let mut features = prefilled();
        algorithm
            .run(ffc1_input(), &mut features, &parameters, &FeatureMap::new())
            .unwrap();
        inexact += assert_dumps_match(
            &fixture(expected),
            &dump_map(&features, algorithm.aborts()),
            case,
        );
        assert_eq!(features.unique_id, 777);
        assert_eq!(
            features.metadata.get("caller_map_meta"),
            Some(&MetaValue::from("kept".to_owned()))
        );
        assert!(algorithm.report().contains(&ReportLine::Warn(
            "Could not assign 'spectrum_native_id' for 1 feature(s), because the computed apex \
             spectrum index was out of range."
                .into()
        )));
    }
    assert_eq!(fixture("prefilled_info.txt"), "uid 777 meta kept\n");
    eprintln!("prefilled: {inexact} fitted values not bit-identical");
}

// ---------------------------------------------------------------------------
// Debug mode
// ---------------------------------------------------------------------------

fn sha1_hex(data: &[u8]) -> String {
    use sha1::{Digest, Sha1};
    let mut hasher = Sha1::new();
    hasher.update(data);
    hasher
        .finalize()
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

/// The executed file sizes and SHA-1 digests of `debug_digests.tsv`, by case
/// and file name.
fn debug_digests() -> BTreeMap<(String, String), (usize, String)> {
    // `shift_nonfinite_digests.tsv`: the combined fix round's executed runs
    // with an infinite and a NaN `debug:pseudo_rt_shift`
    // (`../oracle/ffap-complete-fix1/node/run_shift.sh`, two runs each,
    // identical).
    // `crash_digests.tsv`: fix round 3's negative-intensity library runs that
    // die with SIGSEGV (`../oracle/ffap-complete-fix3`, `vfi2_driver neg`,
    // two runs each, identical).
    let text = fixture("debug_digests.tsv")
        + &fixture("shift_nonfinite_digests.tsv")
        + &fixture("crash_digests.tsv");
    text.lines()
        .filter(|line| !line.starts_with("case\t"))
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            (
                (fields[0].to_owned(), fields[1].to_owned()),
                (fields[2].parse().unwrap(), fields[3].to_owned()),
            )
        })
        .collect()
}

/// Assert that `data` is byte for byte the executed file.
fn assert_executed_bytes(case: &str, file: &str, data: &[u8]) {
    let digests = debug_digests();
    let (bytes, sha1) = digests
        .get(&(case.to_owned(), file.to_owned()))
        .unwrap_or_else(|| panic!("no executed digest for {case} {file}"));
    assert_eq!(data.len(), *bytes, "{case} {file}: size");
    assert_eq!(&sha1_hex(data), sha1, "{case} {file}: content");
}

/// The upstream FuzzyDiff tolerance with generated ids excluded (decision D6).
fn d6() -> decoded::DecodedOptions {
    decoded::DecodedOptions::new(decoded::Tolerance::from_settings(
        &fuzzy::FuzzyDiffSettings::upstream().unwrap(),
    ))
    .ignoring_unique_ids()
}

/// A featureXML fixture, decoded.
fn feature_fixture(name: &str) -> FeatureMap {
    featurexml::load(fixture_path(name)).unwrap()
}

/// An mzML fixture, decoded.
fn experiment_fixture(name: &str) -> MSExperiment {
    let bytes = std::fs::read(fixture_path(name)).unwrap();
    let mut xml = Vec::new();
    flate2::read::GzDecoder::new(bytes.as_slice())
        .read_to_end(&mut xml)
        .unwrap();
    let options = openms::format::mzml::ReadOptions {
        source_nonfinite_float_arrays: true,
        ..Default::default()
    };
    openms::format::mzml::read_with_options(xml.as_slice(), &options).unwrap()
}

fn assert_maps_decoded_equal(actual: &FeatureMap, expected: &FeatureMap, case: &str) {
    if let Err(mismatch) = decoded::compare_feature_maps(actual, expected, &d6()) {
        panic!("{case}: {mismatch}");
    }
}

fn assert_experiments_decoded_equal(actual: &MSExperiment, expected: &MSExperiment, case: &str) {
    assert_eq!(actual.spectra.len(), expected.spectra.len(), "{case}");
    for (index, (a, e)) in actual.spectra.iter().zip(&expected.spectra).enumerate() {
        assert_eq!(a.native_id, e.native_id, "{case} spectrum {index}");
        assert_eq!(a.rt.to_bits(), e.rt.to_bits(), "{case} spectrum {index} rt");
        assert_eq!(a.peaks, e.peaks, "{case} spectrum {index} peaks");
        assert_eq!(
            a.float_data_arrays.len(),
            e.float_data_arrays.len(),
            "{case} spectrum {index} arrays"
        );
        for (x, y) in a.float_data_arrays.iter().zip(&e.float_data_arrays) {
            assert_eq!(x.name, y.name, "{case} spectrum {index}");
            assert_eq!(x.data.len(), y.data.len(), "{case} spectrum {index}");
            for (u, v) in x.data.iter().zip(&y.data) {
                // Every array, the overall scores and NaN bits included, is
                // the executed one: the port computes the Release build's
                // `powf` (`seeds::overall_score`) and its SSE NaN rules.
                assert_eq!(
                    u.to_bits(),
                    v.to_bits(),
                    "{case} spectrum {index} array {}: {u} against {v}",
                    x.name
                );
            }
        }
    }
}

/// Run one debug case through a fresh instance.
fn debug_run(
    experiment: MSExperiment,
    parameters: &Param,
    options: Options,
) -> (openms::Result<()>, FeatureFinderAlgorithmPicked, FeatureMap) {
    let mut algorithm = FeatureFinderAlgorithmPicked::with_options(options).unwrap();
    let mut features = FeatureMap::new();
    let result = algorithm.run(experiment, &mut features, parameters, &FeatureMap::new());
    (result, algorithm, features)
}

fn with_debug(mut parameters: Param, overrides: &[(&str, ParamValue)]) -> Param {
    set(
        &mut parameters,
        "write_debug",
        ParamValue::String("true".into()),
    );
    for (key, value) in overrides {
        set(&mut parameters, key, value.clone());
    }
    parameters
}

fn store_points(report: &[ReportLine]) -> Vec<String> {
    report
        .iter()
        .map(|line| match line {
            ReportLine::Out(text) => format!("out {text}"),
            ReportLine::Info(text) => format!("info {text}"),
            ReportLine::Warn(text) => format!("warn {text}"),
            ReportLine::StoreSeedMap(index) => format!("store seeds {index}"),
            ReportLine::StoreAbortReasons => "store aborts".to_owned(),
            ReportLine::StoreInput => "store input".to_owned(),
        })
        .collect()
}

/// Executed case a1: `mass_trace:min_spectra = 1` finds no seed, so no seed
/// reaches the fit and the debug run completes: the log, the empty seed map,
/// the empty abort map and the input with its (NaN) trace scores.
#[test]
fn a_debug_run_without_seeds_matches_the_release_build() {
    let parameters = with_debug(
        ffc1_parameters(),
        &[("mass_trace:min_spectra", ParamValue::Integer(1))],
    );
    let (result, algorithm, features) = debug_run(ffc1_input(), &parameters, Options::default());
    result.unwrap();
    assert!(features.is_empty());
    let out = algorithm.debug_output().unwrap();
    assert!(out.log_opened);
    assert!(out.termination.is_none());
    assert_executed_bytes("a1", "log.txt", out.log.text().as_bytes());
    assert_eq!(out.seed_maps.len(), 1);
    assert_eq!(out.seed_maps[0].charge, 2);
    assert_maps_decoded_equal(
        &out.seed_maps[0].map,
        &feature_fixture("empty_seed_map.featureXML"),
        "a1 seeds",
    );
    assert_maps_decoded_equal(
        out.abort_reasons.as_ref().unwrap(),
        &feature_fixture("empty_abort_map.featureXML"),
        "a1 aborts",
    );
    assert_experiments_decoded_equal(
        out.input.as_ref().unwrap(),
        &experiment_fixture("a1_input.mzML.gz"),
        "a1 input",
    );
    assert!(out.feature_files.is_empty());
    assert_eq!(
        store_points(algorithm.report()),
        [
            "store seeds 0",
            "out Found 0 seeds for charge 2.",
            "out Found 0 feature candidates for charge 2.",
            "info Removed 0 overlapping features.",
            "info ",
            "info Info: reasons for not finalizing a feature during its construction:",
            "info ",
            "info 0 features found.",
            "store aborts",
            "store input",
        ]
    );
}

/// Executed case a2: `feature:min_isotope_fit = 1` aborts all 25 seeds before
/// the fit; the abort map holds them by ascending intensity.
#[test]
fn a_debug_run_whose_seeds_all_abort_matches_the_release_build() {
    let parameters = with_debug(
        ffc1_parameters(),
        &[("feature:min_isotope_fit", ParamValue::Float(1.0))],
    );
    let (result, algorithm, _) = debug_run(ffc1_input(), &parameters, Options::default());
    result.unwrap();
    let out = algorithm.debug_output().unwrap();
    assert_executed_bytes("a2", "log.txt", out.log.text().as_bytes());
    assert_maps_decoded_equal(
        &out.seed_maps[0].map,
        &feature_fixture("ffc1_seed_map.featureXML.gz"),
        "a2 seeds",
    );
    let aborts = out.abort_reasons.as_ref().unwrap();
    assert_maps_decoded_equal(
        aborts,
        &feature_fixture("a2_abort_map.featureXML.gz"),
        "a2 aborts",
    );
    // The source numbers the abort features 0, 1, ...; id 0 is the invalid id.
    let ids: Vec<u64> = aborts.features.iter().map(|f| f.unique_id).collect();
    assert_eq!(ids, (0..25).collect::<Vec<u64>>());
    assert_eq!(algorithm.abort_reasons().len(), 25);
    assert_experiments_decoded_equal(
        out.input.as_ref().unwrap(),
        &experiment_fixture("a2_input.mzML.gz"),
        "a2 input",
    );
}

/// Executed case a3: four scans, all four charges, no seed.
#[test]
fn a_debug_run_on_a_short_input_matches_the_release_build() {
    let parameters = with_debug(Param::new(), &[]);
    let (result, algorithm, _) = debug_run(short_input(), &parameters, Options::default());
    result.unwrap();
    let out = algorithm.debug_output().unwrap();
    assert_executed_bytes("a3", "log.txt", out.log.text().as_bytes());
    assert_eq!(
        out.seed_maps.iter().map(|m| m.charge).collect::<Vec<_>>(),
        [1, 2, 3, 4]
    );
    for seeds in &out.seed_maps {
        assert_maps_decoded_equal(
            &seeds.map,
            &feature_fixture("empty_seed_map.featureXML"),
            "a3 seeds",
        );
    }
    assert_experiments_decoded_equal(
        out.input.as_ref().unwrap(),
        &experiment_fixture("a3_input.mzML.gz"),
        "a3 input",
    );
}

/// Executed cases b1 and a4: a seed reaches the fit, `writeFeatureDebugInfo_`
/// throws `ElementNotFound` inside the OpenMP region and the executed process
/// is killed by SIGABRT. The port refuses at that seed. What the executed
/// build wrote before it died is reproduced: the seed map of the charge being
/// extended and the log up to the last byte the file buffer had flushed,
/// which the port's buffer model predicts exactly.
#[test]
fn a_debug_run_that_reaches_the_fit_stops_where_the_release_build_terminates() {
    for (case, experiment, parameters, charge, seed_map) in [
        (
            "b1",
            ffc1_input(),
            with_debug(ffc1_parameters(), &[]),
            2,
            "ffc1_seed_map.featureXML.gz",
        ),
        (
            "a4",
            ffc1_input(),
            with_debug(
                Param::new(),
                &[("feature:min_isotope_fit", ParamValue::Float(1.0))],
            ),
            1,
            "a4_seed_map_1.featureXML.gz",
        ),
    ] {
        let (result, algorithm, _) = debug_run(experiment, &parameters, Options::default());
        let error = result.unwrap_err();
        assert!(
            matches!(error, openms::Error::Unsupported(_)),
            "{case}: {error}"
        );
        assert!(
            error
                .to_string()
                .contains("the element 'debug:pseudo_rt_shift' could not be found")
        );
        let out = algorithm.debug_output().unwrap();
        let termination = out.termination.as_ref().unwrap();
        // The first seed that reaches the fit, which receives plot 0.
        assert!(
            matches!(
                termination.point,
                TerminationPoint::Seed { charge: c, plot_nr: 0, .. } if c == charge
            ),
            "{case}: {:?}",
            termination.point
        );
        assert_eq!(termination.kind, TerminationKind::Exception, "{case}");
        assert_eq!(
            termination.log_file_bytes,
            Some(out.log.flushed_bytes()),
            "{case}"
        );
        assert_eq!(algorithm.termination(), Some(termination), "{case}");
        assert_eq!(termination.exception, "ElementNotFound");
        assert_eq!(
            termination.message,
            "the element 'debug:pseudo_rt_shift' could not be found"
        );
        // `what.txt`: the executed exception text.
        assert!(fixture("what.txt").contains(&format!(
            "missing ElementNotFound {}\n",
            termination.message
        )));
        let flushed = &out.log.text().as_bytes()[..out.log.flushed_bytes()];
        assert_executed_bytes(case, "log.txt", flushed);
        assert!(out.log.len() > out.log.flushed_bytes());
        assert_eq!(out.seed_maps.len(), 1, "{case}");
        assert_eq!(out.seed_maps[0].charge, charge);
        assert_maps_decoded_equal(&out.seed_maps[0].map, &feature_fixture(seed_map), case);
        assert!(out.abort_reasons.is_none());
        assert!(out.input.is_none());
        assert!(out.feature_files.is_empty());
        let report = store_points(algorithm.report());
        assert_eq!(
            report,
            [
                "store seeds 0".to_owned(),
                format!(
                    "out Found {} seeds for charge {charge}.",
                    out.seed_maps[0].map.len()
                ),
            ]
        );
    }
    // `ParamValue::EMPTY` under the key: the executed ConversionError.
    let mut parameters = with_debug(ffc1_parameters(), &[]);
    set(&mut parameters, "debug:pseudo_rt_shift", ParamValue::Empty);
    let (result, algorithm, _) = debug_run(ffc1_input(), &parameters, Options::default());
    assert!(result.is_err());
    let termination = algorithm
        .debug_output()
        .unwrap()
        .termination
        .clone()
        .unwrap();
    assert_eq!(termination.exception, "ConversionError");
    assert!(
        fixture("what.txt").contains(&format!("empty ConversionError {}\n", termination.message))
    );
    // A string under the key does not terminate the source; see
    // `a_heap_address_shift_matches_the_release_build_where_it_is_reproducible`.
}

/// A driver case `declared-*`: its name, its parameter overrides and the
/// fixture of its executed feature map.
type DeclaredCase = (&'static str, Vec<(&'static str, ParamValue)>, &'static str);

/// The executed driver cases `declared-*`: with `debug:pseudo_rt_shift` in the
/// parameters the source writes its feature files, and so does the port under
/// either key policy. Every `.dta`, `_cropped.dta` and `.plot` file and the
/// log are byte-identical to the executed ones; the feature map is the plain
/// run's.
///
/// The `shift_*` cases hold a non-finite shift (`shift_nonfinite_digests.tsv`,
/// under [`PseudoRtShiftKey::Source`] only, since the declared key is
/// restricted to at least 1): `inf * 0.0` is x86_64's default NaN, whose sign
/// bit is set, so trace 0 of every plot is centred at a negative NaN, and the
/// `.plot` formulas spell it `-nan`, as glibc's `operator<<` does; a negative
/// NaN shift writes `-nan` in every trace, a positive one `nan`.
#[test]
fn feature_debug_files_match_the_release_build() {
    let digests = debug_digests();
    let cases: [DeclaredCase; 8] = [
        (
            "shift500",
            vec![("debug:pseudo_rt_shift", ParamValue::Float(500.0))],
            "reuse_run1.txt.gz",
        ),
        (
            "shift123",
            vec![
                ("debug:pseudo_rt_shift", ParamValue::Float(123.25)),
                ("advanced:pseudo_rt_shift", ParamValue::Float(123.25)),
            ],
            "reuse_run1.txt.gz",
        ),
        (
            "int250",
            vec![
                ("debug:pseudo_rt_shift", ParamValue::Integer(250)),
                ("advanced:pseudo_rt_shift", ParamValue::Float(250.0)),
            ],
            "reuse_run1.txt.gz",
        ),
        (
            "egh",
            vec![
                ("debug:pseudo_rt_shift", ParamValue::Float(500.0)),
                ("feature:rt_shape", ParamValue::String("asymmetric".into())),
            ],
            "declared_egh_map.txt.gz",
        ),
        (
            "shift_inf",
            vec![("debug:pseudo_rt_shift", ParamValue::Float(f64::INFINITY))],
            "reuse_run1.txt.gz",
        ),
        (
            "shift_neginf",
            vec![(
                "debug:pseudo_rt_shift",
                ParamValue::Float(f64::NEG_INFINITY),
            )],
            "reuse_run1.txt.gz",
        ),
        (
            "shift_negnan",
            vec![(
                "debug:pseudo_rt_shift",
                ParamValue::Float(f64::from_bits(0xfff8_0000_0000_0000)),
            )],
            "reuse_run1.txt.gz",
        ),
        (
            "shift_nan",
            vec![(
                "debug:pseudo_rt_shift",
                ParamValue::Float(f64::from_bits(0x7ff8_0000_0000_0000)),
            )],
            "reuse_run1.txt.gz",
        ),
    ];
    let mut inexact = 0;
    for (case, overrides, expected_map) in cases {
        let policies: &[PseudoRtShiftKey] = if case.starts_with("shift_") {
            &[PseudoRtShiftKey::Source]
        } else {
            &[PseudoRtShiftKey::Source, PseudoRtShiftKey::Declared]
        };
        for &policy in policies {
            let mut overrides = overrides.clone();
            if policy == PseudoRtShiftKey::Declared {
                // The declared key alone; the undeclared one must not matter.
                // shift500 and egh leave advanced:pseudo_rt_shift at its
                // default 500.
                overrides.retain(|(key, _)| *key != "debug:pseudo_rt_shift");
            }
            let parameters = with_debug(ffc1_parameters(), &overrides);
            let options = Options {
                pseudo_rt_shift: policy,
                ..Options::default()
            };
            let (result, algorithm, features) = debug_run(ffc1_input(), &parameters, options);
            result.unwrap();
            let out = algorithm.debug_output().unwrap();
            assert_executed_bytes(case, "log.txt", out.log.text().as_bytes());
            let expected_files: Vec<&String> = digests
                .keys()
                .filter(|(c, f)| c == case && f.starts_with("features/"))
                .map(|(_, f)| f)
                .collect();
            let mut produced = 0;
            for files in &out.feature_files {
                assert_eq!(files.path, "debug/features/");
                let name = |full: String| full.trim_start_matches("debug/").to_owned();
                assert_executed_bytes(case, &name(files.dta_name()), files.dta.as_bytes());
                assert_executed_bytes(case, &name(files.plot_name()), &files.plot);
                produced += 2;
                if let Some(cropped) = &files.cropped_dta {
                    assert_executed_bytes(
                        case,
                        &name(files.cropped_dta_name()),
                        cropped.as_bytes(),
                    );
                    produced += 1;
                }
            }
            assert_eq!(produced, expected_files.len(), "{case}: file count");
            // The executed maps of shift500, shift123 and int250 are the plain
            // first run's, byte for byte (make_fixtures.py checks it).
            inexact += assert_dumps_match(
                &fixture(expected_map),
                &dump_map(&features, algorithm.aborts()),
                case,
            );
            let abort_fixture = match case {
                "egh" => "declared_egh_abort_map.featureXML",
                _ => "declared_shift500_abort_map.featureXML",
            };
            assert_maps_decoded_equal(
                out.abort_reasons.as_ref().unwrap(),
                &feature_fixture(abort_fixture),
                case,
            );
        }
    }
    eprintln!("declared: {inexact} fitted values not bit-identical");
}

/// The executed driver case `declared-prefilled`: the feature files number the
/// accepted features after the caller's eleven.
#[test]
fn feature_debug_files_count_the_callers_features() {
    let parameters = with_debug(
        ffc1_parameters(),
        &[("debug:pseudo_rt_shift", ParamValue::Float(500.0))],
    );
    let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
    let mut features = prefilled();
    algorithm
        .run(ffc1_input(), &mut features, &parameters, &FeatureMap::new())
        .unwrap();
    let out = algorithm.debug_output().unwrap();
    assert_executed_bytes("prefilled", "log.txt", out.log.text().as_bytes());
    for files in &out.feature_files {
        let name = |full: String| full.trim_start_matches("debug/").to_owned();
        assert_executed_bytes("prefilled", &name(files.dta_name()), files.dta.as_bytes());
        assert_executed_bytes("prefilled", &name(files.plot_name()), &files.plot);
        if let Some(cropped) = &files.cropped_dta {
            assert_executed_bytes(
                "prefilled",
                &name(files.cropped_dta_name()),
                cropped.as_bytes(),
            );
        }
    }
    assert!(
        String::from_utf8_lossy(&out.feature_files[0].plot).contains("title 'feature 12 (score: ")
    );
    assert_dumps_match(
        &fixture("prefilled_run.txt.gz"),
        &dump_map(&features, algorithm.aborts()),
        "declared prefilled",
    );
    assert_maps_decoded_equal(
        out.abort_reasons.as_ref().unwrap(),
        &feature_fixture("declared_prefilled_abort_map.featureXML"),
        "declared prefilled",
    );
}

/// The executed driver case `debug_twice`: the second debug run of one
/// object cannot reopen the log, so the file keeps the first run's text; the
/// abort map accumulates over both runs; a third run without debug changes no
/// debug file. While the object lived, the file held the flushed prefix.
#[test]
fn a_second_debug_run_of_one_object_matches_the_release_build() {
    let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
    let first = with_debug(
        ffc1_parameters(),
        &[("feature:min_isotope_fit", ParamValue::Float(1.0))],
    );
    let mut features = FeatureMap::new();
    algorithm
        .run(ffc1_input(), &mut features, &first, &FeatureMap::new())
        .unwrap();
    let run1 = algorithm.take_debug_output().unwrap();
    assert!(run1.log_opened);
    assert_executed_bytes("twice", "log.txt", run1.log.text().as_bytes());
    assert_executed_bytes(
        "twice",
        "log_after_run1.txt",
        &run1.log.text().as_bytes()[..run1.log.flushed_bytes()],
    );
    assert_eq!(algorithm.abort_reasons().len(), 25);

    let mut second = first.clone();
    set(
        &mut second,
        "isotopic_pattern:charge_low",
        ParamValue::Integer(1),
    );
    set(
        &mut second,
        "isotopic_pattern:charge_high",
        ParamValue::Integer(3),
    );
    let mut features = FeatureMap::new();
    algorithm
        .run(ffc1_input(), &mut features, &second, &FeatureMap::new())
        .unwrap();
    let run2 = algorithm.take_debug_output().unwrap();
    assert!(!run2.log_opened);
    assert!(run2.log.is_empty());
    assert_eq!(algorithm.abort_reasons().len(), 25);
    for (index, name) in [
        "twice_seed_map_1.featureXML.gz",
        "twice_seed_map_2.featureXML.gz",
        "empty_seed_map.featureXML",
    ]
    .iter()
    .enumerate()
    {
        assert_maps_decoded_equal(&run2.seed_maps[index].map, &feature_fixture(name), name);
    }
    assert_maps_decoded_equal(
        run2.abort_reasons.as_ref().unwrap(),
        &feature_fixture("twice_abort_map.featureXML.gz"),
        "twice aborts",
    );
    assert_experiments_decoded_equal(
        run2.input.as_ref().unwrap(),
        &experiment_fixture("twice_input.mzML.gz"),
        "twice input",
    );
    assert_eq!(
        algorithm
            .aborts()
            .get("Could not find good enough isotope pattern containing the seed"),
        Some(&55)
    );

    let mut third = ffc1_parameters();
    set(
        &mut third,
        "feature:min_isotope_fit",
        ParamValue::Float(1.0),
    );
    let mut features = FeatureMap::new();
    algorithm
        .run(ffc1_input(), &mut features, &third, &FeatureMap::new())
        .unwrap();
    assert!(algorithm.debug_output().is_none());
    assert_eq!(
        algorithm
            .aborts()
            .get("Could not find good enough isotope pattern containing the seed"),
        Some(&60)
    );
    assert!(
        algorithm
            .report()
            .contains(&ReportLine::Out("Found 5 seeds for charge 2.".into()))
    );
}

/// The debug output does not depend on the worker count: the port collects
/// every seed's lines and aborts in seed order, which is the source's
/// single-thread output. The executed C++ at four threads is undefined there
/// (its logs differed from run to run: cases c1 and c2).
#[test]
fn debug_output_is_identical_at_every_thread_count() {
    let parameters = with_debug(
        ffc1_parameters(),
        &[("debug:pseudo_rt_shift", ParamValue::Float(500.0))],
    );
    let reference = debug_run(
        ffc1_input(),
        &parameters,
        Options {
            threads: Threads::from_cli(1),
            ..Options::default()
        },
    );
    let reference_out: DebugOutput = reference.1.debug_output().unwrap().clone();
    for threads in [2, 8] {
        let run = debug_run(
            ffc1_input(),
            &parameters,
            Options {
                threads: Threads::from_cli(threads),
                ..Options::default()
            },
        );
        run.0.unwrap();
        assert_eq!(run.1.debug_output().unwrap(), &reference_out);
        assert_eq!(run.2, reference.2);
    }
}

// ---------------------------------------------------------------------------
// The ProgressLogger base
// ---------------------------------------------------------------------------

/// A writer the progress logger and the console sink share.
#[derive(Clone, Default)]
struct Shared(Arc<Mutex<Vec<u8>>>);

impl std::io::Write for Shared {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// The deterministic part of a `CMD` progress transcript: the percentage
/// updates, which the logger prints only when a wall-clock second has passed,
/// are removed, and the timings of the `done` lines are masked.
fn mask_progress(text: &str) -> Vec<String> {
    let mut lines = Vec::new();
    for raw in text.split('\n') {
        let mut line = raw.to_owned();
        // Drop every `\r<indent>NN.NN %<15 spaces>` update.
        while let Some(start) = line.find('\r') {
            let rest = &line[start + 1..];
            let trimmed = rest.trim_start_matches(' ');
            let digits = trimmed.find(" %               ").filter(|&end| {
                trimmed[..end]
                    .chars()
                    .all(|c| c.is_ascii_digit() || c == '.' || c == '-')
            });
            match digits {
                Some(end) => {
                    let cut =
                        start + 1 + (rest.len() - trimmed.len()) + end + " %               ".len();
                    line.replace_range(start..cut, "");
                }
                None => {
                    line.replace_range(start..start + 1, "");
                }
            }
        }
        if let (Some(open), Some(close)) = (line.find("[took "), line.find("] -- ")) {
            line.replace_range(open..close + 1, "[took <timing>]");
        }
        lines.push(line);
    }
    while lines.last().is_some_and(String::is_empty) {
        lines.pop();
    }
    lines
}

/// `setLogType(CMD)`: the executed progress transcript of three runs (the
/// FeatureFinderCentroided_1 run, a four-scan input whose seed ranges are
/// inverted, and a run into the caller's map), interleaved with the
/// `std::cout` seed and candidate lines, has the same structure here: the
/// same labels in the same order, each closed by its `done` line.
#[test]
fn the_progress_transcript_matches_the_release_build() {
    let executed = fixture("progress_stdout.txt");
    let mut sections: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut current: Option<(String, String)> = None;
    for line in executed.split_inclusive('\n') {
        if let Some(name) = line.strip_prefix("== CASE ") {
            current = Some((name.trim().to_owned(), String::new()));
        } else if line.starts_with("== END ") {
            let (name, text) = current.take().unwrap();
            sections.insert(name, mask_progress(&text));
        } else if let Some((_, text)) = current.as_mut() {
            text.push_str(line);
        }
    }
    assert_eq!(sections.len(), 3);
    for (name, experiment, parameters, mut features) in [
        ("ffc1", ffc1_input(), ffc1_parameters(), FeatureMap::new()),
        ("short", short_input(), Param::new(), FeatureMap::new()),
        ("prefilled", ffc1_input(), ffc1_parameters(), prefilled()),
    ] {
        let shared = Shared::default();
        let mut logger = ProgressLogger::new();
        logger.set_log_type(ProgressLogType::Cmd);
        logger.set_logger(Box::new(CommandProgressLogger::new(shared.clone())));
        let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
        assert_eq!(algorithm.log_type(), ProgressLogType::None);
        algorithm.set_progress_logger(Some(logger));
        assert_eq!(algorithm.log_type(), ProgressLogType::Cmd);
        let sink = shared.clone();
        algorithm.set_console(Some(Box::new(move |line: &ReportLine| {
            // A library program prints only `std::cout` lines on standard
            // output; the OpenMS info log is not attached to it.
            if let ReportLine::Out(text) = line {
                let mut sink = sink.clone();
                use std::io::Write;
                writeln!(sink, "{text}").unwrap();
            }
        })));
        algorithm
            .run(experiment, &mut features, &parameters, &FeatureMap::new())
            .unwrap();
        let produced = String::from_utf8(shared.0.lock().unwrap().clone()).unwrap();
        assert_eq!(mask_progress(&produced), sections[name], "{name}");
    }
}

/// The default progress type is `NONE`, and a `NONE` logger is removed
/// rather than called, so it costs nothing.
#[test]
fn progress_is_silent_by_default() {
    let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
    assert_eq!(algorithm.log_type(), ProgressLogType::None);
    algorithm.set_log_type(ProgressLogType::Cmd);
    assert_eq!(algorithm.log_type(), ProgressLogType::Cmd);
    algorithm.set_log_type(ProgressLogType::None);
    assert!(algorithm.progress_logger_mut().is_none());
    let mut features = FeatureMap::new();
    algorithm
        .run(
            short_input(),
            &mut features,
            &Param::new(),
            &FeatureMap::new(),
        )
        .unwrap();
    assert!(features.is_empty());
}

/// The fresh-object convenience entry point returns the debug output with the
/// features, and equals a fresh instance run.
#[test]
fn the_convenience_run_returns_the_debug_output() {
    let parameters = with_debug(
        ffc1_parameters(),
        &[("feature:min_isotope_fit", ParamValue::Float(1.0))],
    );
    let output = run_with_options(
        ffc1_input(),
        &FeatureMap::new(),
        &parameters,
        &Options::default(),
    )
    .unwrap();
    let (result, algorithm, features) = debug_run(ffc1_input(), &parameters, Options::default());
    result.unwrap();
    assert_eq!(output.features, features);
    assert_eq!(output.debug.as_ref(), algorithm.debug_output());
    assert_eq!(
        output
            .aborts
            .get("Could not find good enough isotope pattern containing the seed"),
        Some(&25)
    );
    let texts: Vec<String> = algorithm
        .report()
        .iter()
        .filter_map(|line| line.text().map(str::to_owned))
        .collect();
    assert_eq!(output.log, texts);
    // The convenience run refuses where the source terminates, too.
    let terminated = run_with_options(
        ffc1_input(),
        &FeatureMap::new(),
        &with_debug(ffc1_parameters(), &[]),
        &Options::default(),
    );
    assert!(matches!(terminated, Err(openms::Error::Unsupported(_))));
}

// ---------------------------------------------------------------------------
// A caller's map that overlaps the new features, and its undefined variants
// ---------------------------------------------------------------------------

/// The driver's `prefilledOverlapping(variant)` map.
fn overlapping(variant: &str) -> FeatureMap {
    let mut features = Vec::new();
    let mut f = synthetic(4407.0, 646.24, 1.0e7, 2, 0.9, "over-0-same-charge");
    f.convex_hulls
        .push(hull(&[(4374.19, 646.229), (4443.42, 646.2585)]));
    f.convex_hulls
        .push(hull(&[(4370.78, 646.7377), (4443.42, 646.7672)]));
    features.push(f);
    let mut f = synthetic(4389.0, 648.25, 10.0, 4, 0.1, "over-1-z4");
    f.convex_hulls
        .push(hull(&[(4360.0, 648.25), (4420.0, 648.26)]));
    f.convex_hulls
        .push(hull(&[(4360.0, 648.75), (4420.0, 648.76)]));
    features.push(f);
    let mut f = synthetic(4301.0, 651.75, 5.0e4, 3, 0.99, "over-2-z3");
    f.convex_hulls
        .push(hull(&[(4280.0, 651.75), (4320.0, 651.77)]));
    f.convex_hulls
        .push(hull(&[(4280.0, 652.08), (4320.0, 652.10)]));
    features.push(f);
    let mut f = synthetic(4278.0, 653.70, 100.0, 2, 0.5, "over-3-empty-hull");
    f.convex_hulls.push(ConvexHull2D::new());
    f.convex_hulls
        .push(hull(&[(4250.0, 653.77), (4300.0, 653.78)]));
    features.push(f);
    features.push(synthetic(4201.0, 652.70, 500.0, 2, 0.5, "over-4-no-hull"));
    let mut extra =
        |rt: f64, mz: f64, intensity: f32, charge: i32, label: &str, points: &[(f64, f64)]| {
            let mut f = synthetic(rt, mz, intensity, charge, 0.5, label);
            f.convex_hulls.push(hull(points));
            features.push(f);
        };
    match variant {
        "nan_mz" => extra(
            4200.0,
            f64::NAN,
            700.0,
            2,
            "ub-nan-mz",
            &[(4190.0, 660.0), (4210.0, 660.01)],
        ),
        "nan_intensity" => extra(
            4200.0,
            660.0,
            f32::NAN,
            2,
            "ub-nan-intensity",
            &[(4190.0, 660.0), (4210.0, 660.01)],
        ),
        "zero_charge" => extra(
            4201.8,
            652.70,
            700.0,
            0,
            "ub-zero-charge",
            &[(4150.0, 652.70), (4260.0, 652.80)],
        ),
        "odd_rt" => {
            extra(
                f64::INFINITY,
                700.0,
                700.0,
                2,
                "rt-inf",
                &[(4190.0, 700.0), (4210.0, 700.01)],
            );
            extra(
                f64::NAN,
                701.0,
                700.0,
                2,
                "rt-nan",
                &[(4190.0, 701.0), (4210.0, 701.01)],
            );
            extra(
                f64::NEG_INFINITY,
                702.0,
                700.0,
                2,
                "rt-minus-inf",
                &[(4190.0, 702.0), (4210.0, 702.01)],
            );
        }
        _ => {}
    }
    FeatureMap::from_features(features)
}

/// The driver's `overlap` mode: the caller's features overlap new ones
/// under all three resolution rules (same charge, multiple of the charge,
/// quality), an empty hull makes a full-plane box that meets every feature,
/// and `feature:max_intersection = 0` resolves the `-0.0` overlap that box
/// yields. NaN m/z, NaN intensity and infinite or NaN retention times are
/// defined inputs here (the sorts stay inside the map, `lower_bound` is
/// defined for them), and the port matches the executed maps.
#[test]
fn an_overlapping_caller_map_matches_the_release_build() {
    let mut inexact = 0;
    for variant in ["ffc1", "zero", "nan_mz", "nan_intensity", "odd_rt"] {
        let mut parameters = ffc1_parameters();
        if variant == "zero" {
            set(
                &mut parameters,
                "feature:max_intersection",
                ParamValue::Float(0.0),
            );
        }
        let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
        let mut features = overlapping(variant);
        algorithm
            .run(ffc1_input(), &mut features, &parameters, &FeatureMap::new())
            .unwrap_or_else(|error| panic!("{variant}: {error}"));
        inexact += assert_dumps_match(
            &fixture(&format!("overlap_{variant}.txt.gz")),
            &dump_map(&features, algorithm.aborts()),
            variant,
        );
    }
    eprintln!("overlap: {inexact} fitted values not bit-identical");
}

/// The driver's `overlap zero_charge`: a caller's feature of charge 0 meets a
/// charge-2 feature, and the source's `2 % 0` traps; the executed driver is
/// killed by SIGFPE (signal 8) in both repetitions. The port refuses at that
/// pair.
#[test]
fn a_charge_zero_overlap_is_refused_where_the_release_build_traps() {
    let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
    let mut features = overlapping("zero_charge");
    let error = algorithm
        .run(
            ffc1_input(),
            &mut features,
            &ffc1_parameters(),
            &FeatureMap::new(),
        )
        .unwrap_err();
    assert!(
        matches!(&error, openms::Error::InvalidValue(message) if message.contains("% 0") && message.contains("SIGFPE")),
        "{error}"
    );
    // The process ends there: a termination without debug output, and no
    // debug run of this instance has opened `debug/log.txt`.
    let termination = algorithm.termination().unwrap();
    assert_eq!(termination.kind, TerminationKind::ArithmeticTrap);
    assert_eq!(termination.exception, "SIGFPE");
    assert_eq!(termination.message, error.to_string());
    assert_eq!(termination.log_file_bytes, None);
    let TerminationPoint::OverlapResolution { first, second } = termination.point else {
        panic!("{:?}", termination.point);
    };
    assert!(first < second && second < features.len());
    assert!(features.features[first].charge == 0 || features.features[second].charge == 0);
    assert!(algorithm.debug_output().is_none());
    // The source's state at the trap: the new features are in the map, which
    // step 4 has sorted by m/z.
    assert!(features.len() > overlapping("zero_charge").len());
    assert!(
        features
            .features
            .windows(2)
            .all(|pair| pair[0].mz <= pair[1].mz)
    );
}

// ---------------------------------------------------------------------------
// Stale abort-reason seeds
// ---------------------------------------------------------------------------

fn scaled_ffc1_input() -> MSExperiment {
    let mut experiment = ffc1_input();
    for spectrum in &mut experiment.spectra {
        for peak in &mut spectrum.peaks {
            peak.intensity *= 2.0;
        }
    }
    experiment
}

/// The driver's `stale scaled`: a second debug run of one object on the same
/// scans with doubled intensities. The 25 seeds of the first run stay in
/// `abort_reasons_` under their old intensities, the second run adds its own
/// under the doubled ones, and the abort map reads all of them from the
/// second input.
#[test]
fn stale_abort_seeds_are_read_from_the_current_input() {
    let parameters = with_debug(
        ffc1_parameters(),
        &[("feature:min_isotope_fit", ParamValue::Float(1.0))],
    );
    let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
    let mut first = FeatureMap::new();
    algorithm
        .run(ffc1_input(), &mut first, &parameters, &FeatureMap::new())
        .unwrap();
    assert_eq!(
        format!("{}\n", algorithm.abort_reasons().len()),
        fixture("stale_scaled_run1_count.txt")
    );
    let mut second = FeatureMap::new();
    algorithm
        .run(
            scaled_ffc1_input(),
            &mut second,
            &parameters,
            &FeatureMap::new(),
        )
        .unwrap();
    assert_eq!(
        format!("{}\n", algorithm.abort_reasons().len()),
        fixture("stale_scaled_run2_count.txt")
    );
    assert_dumps_match(
        &fixture("stale_scaled_run2_map.txt.gz"),
        &dump_map(&second, algorithm.aborts()),
        "stale scaled",
    );
    let out = algorithm.debug_output().unwrap();
    assert!(!out.log_opened);
    assert_maps_decoded_equal(
        &out.seed_maps[0].map,
        &feature_fixture("stale_scaled_seed_map_2.featureXML.gz"),
        "stale seeds",
    );
    assert_maps_decoded_equal(
        out.abort_reasons.as_ref().unwrap(),
        &feature_fixture("stale_scaled_abort_map.featureXML.gz"),
        "stale aborts",
    );
}

/// The driver's `stale oob`: the second debug run reads four scans, while
/// the first run's seeds address scans up to about 100. The source reads
/// them without a bounds check; the executed driver is killed by SIGSEGV in
/// all five repetitions. The port refuses at the abort map, after everything
/// the source does before it.
#[test]
fn stale_abort_seeds_outside_the_input_are_refused() {
    let first_parameters = with_debug(
        ffc1_parameters(),
        &[("feature:min_isotope_fit", ParamValue::Float(1.0))],
    );
    let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
    let mut first = FeatureMap::new();
    algorithm
        .run(
            ffc1_input(),
            &mut first,
            &first_parameters,
            &FeatureMap::new(),
        )
        .unwrap();
    let mut second = FeatureMap::new();
    let error = algorithm
        .run(
            short_input(),
            &mut second,
            &with_debug(Param::new(), &[]),
            &FeatureMap::new(),
        )
        .unwrap_err();
    assert!(
        matches!(&error, openms::Error::InvalidValue(message) if message.contains("abort_reasons_")),
        "{error}"
    );
    let out = algorithm.debug_output().unwrap();
    // The out-of-bounds read of the first stored seed ends the process, and
    // the file keeps what the first run's buffer had flushed.
    let termination = out.termination.as_ref().unwrap();
    assert_eq!(algorithm.termination(), Some(termination));
    assert_eq!(termination.point, TerminationPoint::AbortMap { entry: 0 });
    assert_eq!(termination.kind, TerminationKind::OutOfBounds);
    assert_eq!(termination.exception, "SIGSEGV");
    assert_eq!(termination.message, error.to_string());
    let stream = algorithm.debug_log_file().unwrap();
    assert!(stream.failed);
    assert_eq!(termination.log_file_bytes, Some(stream.flushed));
    // The executed run wrote the four seed maps before it crashed, and no
    // abort map or input of its own.
    assert_eq!(out.seed_maps.len(), 4);
    assert!(out.abort_reasons.is_none());
    assert!(out.input.is_none());
    assert!(
        algorithm
            .report()
            .contains(&ReportLine::Info("0 features found.".into()))
    );
}

// ---------------------------------------------------------------------------
// Where the process ends outside the seed loop, and the never-closed stream
// ---------------------------------------------------------------------------

/// One case of `termination_digests.tsv.gz`: the executed exit status, the
/// driver's report lines, and every file below `debug/` with its size and
/// SHA-1 (`None` for the abort map, whose random unique id changes its bytes).
#[derive(Default)]
struct Executed {
    status: String,
    report: Vec<String>,
    files: BTreeMap<String, Option<(usize, String)>>,
}

impl Executed {
    /// The number after `<prefix>` in a report line, if the driver printed it.
    fn reported(&self, prefix: &str) -> Option<usize> {
        self.report
            .iter()
            .find_map(|line| line.strip_prefix(prefix)?.parse().ok())
    }

    fn printed(&self, line: &str) -> bool {
        self.report.iter().any(|printed| printed == line)
    }
}

/// `termination_digests.tsv.gz` (`../oracle/ffap-complete-fix5`, `fix5_driver`
/// and the unchanged `ffap_instr_driver`, every case twice, identical but for
/// the abort map's unique id; `extract/make_fixtures.py`).
fn termination_digests() -> BTreeMap<String, Executed> {
    let mut cases: BTreeMap<String, Executed> = BTreeMap::new();
    for line in fixture("termination_digests.tsv.gz").lines().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();
        let case = cases.entry(fields[0].to_owned()).or_default();
        match fields[1] {
            "status" => case.status = fields[4].to_owned(),
            "report" => case.report.push(fields[2].to_owned()),
            "file" => {
                let digest =
                    (fields[3] != "-").then(|| (fields[3].parse().unwrap(), fields[4].to_owned()));
                case.files.insert(fields[2].to_owned(), digest);
            }
            other => panic!("unknown row kind {other}"),
        }
    }
    cases
}

/// The files below `debug/` that a caller writes from the port's output,
/// following the `DebugOutput` protocol: a run that opened the stream
/// replaces `log.txt` with its log, every stored map and feature file
/// replaces its namesake, and after a termination `log.txt` is cut to
/// `DebugTermination::log_file_bytes`. The featureXML and mzML stores are
/// recorded by name only (their bytes are the shared writers', compared
/// decoded elsewhere).
#[derive(Default)]
struct Disk {
    files: BTreeMap<String, Option<Vec<u8>>>,
}

impl Disk {
    fn write_run(
        &mut self,
        out: Option<&DebugOutput>,
        termination: Option<&openms::analysis::feature_finder_picked::debug::DebugTermination>,
    ) {
        if let Some(out) = out {
            if out.log_opened {
                self.files
                    .insert("log.txt".into(), Some(out.log.text().as_bytes().to_vec()));
            }
            for seeds in &out.seed_maps {
                self.files
                    .insert(format!("seeds_{}.featureXML", seeds.charge), None);
            }
            for files in &out.feature_files {
                let name = |full: String| full.trim_start_matches("debug/").to_owned();
                self.files
                    .insert(name(files.dta_name()), Some(files.dta.as_bytes().to_vec()));
                if let Some(cropped) = &files.cropped_dta {
                    self.files.insert(
                        name(files.cropped_dta_name()),
                        Some(cropped.as_bytes().to_vec()),
                    );
                }
                self.files
                    .insert(name(files.plot_name()), Some(files.plot.clone()));
            }
            if out.abort_reasons.is_some() {
                self.files.insert("abort_reasons.featureXML".into(), None);
            }
            if out.input.is_some() {
                self.files.insert("input.mzML".into(), None);
            }
        }
        if let Some(bytes) = termination.and_then(|termination| termination.log_file_bytes) {
            let log = self
                .files
                .get_mut("log.txt")
                .and_then(Option::as_mut)
                .expect("a termination with a stream length follows a written log");
            assert!(bytes <= log.len());
            log.truncate(bytes);
        }
    }

    /// Every executed file exists here and nothing else; the log and the
    /// feature files are byte for byte the executed ones.
    fn assert_executed(&self, executed: &Executed, case: &str) {
        assert_eq!(
            self.files.keys().collect::<Vec<_>>(),
            executed.files.keys().collect::<Vec<_>>(),
            "{case}: files"
        );
        for (name, data) in &self.files {
            if name.ends_with(".featureXML") || name.ends_with(".mzML") {
                continue;
            }
            let data = data.as_ref().unwrap();
            let (bytes, sha1) = executed.files[name].as_ref().unwrap();
            assert_eq!(data.len(), *bytes, "{case} {name}: size");
            assert_eq!(&sha1_hex(data), sha1, "{case} {name}: content");
        }
    }
}

/// The FFC_1 section with `write_debug` and `debug:pseudo_rt_shift` 500, as
/// `fix5_driver` sets them.
fn ffc1_debug_parameters() -> Param {
    with_debug(
        ffc1_parameters(),
        &[("debug:pseudo_rt_shift", ParamValue::Float(500.0))],
    )
}

/// `fix5_driver`'s caller map: the common part of [`overlapping`] and the
/// variant feature labelled `ub-<variant>`.
fn overlapping_v4(variant: &str) -> FeatureMap {
    let mut map = overlapping("none");
    let (mz, charge) = match variant {
        "zero_charge" => (652.70, 0),
        "zero_charge_hi" => (652.79, 0),
        _ => return map,
    };
    let mut feature = synthetic(4201.8, mz, 700.0, charge, 0.5, &format!("ub-{variant}"));
    feature
        .convex_hulls
        .push(hull(&[(4150.0, 652.70), (4260.0, 652.80)]));
    map.features.push(feature);
    map
}

/// Executed `fix5_driver single` and `ffap_instr_driver stale` runs
/// (`../oracle/ffap-complete-fix5`): the two process-ending points after the
/// seed loop and the one before it, each with a defined neighbour.
///
/// - A caller's charge-0 feature below the found charge-2 feature at
///   m/z 652.766 makes step 4 compute `2 % 0`: SIGFPE (status 136), with and
///   without `write_debug`. The debug run left `debug/log.txt` at 1,163,782
///   bytes, `seeds_2.featureXML` and the files of 25 plots, and no abort map
///   or input. The port refuses at that pair and records an
///   `ArithmeticTrap`; the file protocol reproduces every file. The same
///   feature above it (`0 % 2 == 0`) and no extra feature return with the
///   complete log.
/// - A reused object whose abort seeds lie outside a four-scan input reads
///   them when it builds the abort map: SIGSEGV (status 139), leaving the
///   first run's flushed prefix, 1,114,578 bytes, and the second run's four
///   seed maps. The port records an `OutOfBounds` termination at entry 0
///   with the first run's stream length. On the same scans with doubled
///   intensities the run returns and the object's destructor completes the
///   log (1,118,230 bytes).
/// - `charge_low` 4 with `charge_high` 2 wraps the score-array count to one
///   array, which the source writes past before it creates `debug/`:
///   SIGSEGV and no file. The port records an `OutOfBounds` termination at
///   the score arrays without a stream length.
#[test]
fn process_ending_refusals_outside_the_seed_loop_record_their_termination() {
    let digests = termination_digests();

    // Step 4.
    for (case, debug) in [("zero_charge_dbg", true), ("zero_charge_nd", false)] {
        let executed = &digests[case];
        assert_eq!(executed.status, "136", "{case}: SIGFPE");
        let parameters = if debug {
            ffc1_debug_parameters()
        } else {
            ffc1_parameters()
        };
        let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
        let mut features = overlapping_v4("zero_charge");
        let error = algorithm
            .run(ffc1_input(), &mut features, &parameters, &FeatureMap::new())
            .unwrap_err();
        assert!(error.to_string().contains("2 % 0"), "{case}: {error}");
        let termination = algorithm.termination().unwrap();
        assert!(
            matches!(termination.point, TerminationPoint::OverlapResolution { first, second } if first < second),
            "{case}: {:?}",
            termination.point
        );
        assert_eq!(termination.kind, TerminationKind::ArithmeticTrap, "{case}");
        assert_eq!(termination.exception, "SIGFPE", "{case}");
        assert_eq!(termination.message, error.to_string(), "{case}");
        let mut disk = Disk::default();
        let out = algorithm.debug_output();
        assert_eq!(out.is_some(), debug, "{case}");
        if let Some(out) = out {
            assert_eq!(out.termination.as_ref(), Some(termination), "{case}");
            assert!(out.log_opened);
            assert_eq!(termination.log_file_bytes, Some(out.log.flushed_bytes()));
            assert!(out.log.flushed_bytes() < out.log.len());
            assert_eq!(out.seed_maps.len(), 1);
            assert_maps_decoded_equal(
                &out.seed_maps[0].map,
                &feature_fixture("zero_charge_seed_map_2.featureXML.gz"),
                case,
            );
            assert!(out.abort_reasons.is_none() && out.input.is_none());
            let plots: std::collections::BTreeSet<i64> = out
                .feature_files
                .iter()
                .map(|files| files.plot_nr)
                .collect();
            assert_eq!(plots.len(), 25, "{case}");
        } else {
            assert_eq!(termination.log_file_bytes, None, "{case}");
        }
        disk.write_run(out, Some(termination));
        disk.assert_executed(executed, case);
    }
    for case in ["zero_charge_hi_dbg", "none_dbg"] {
        let executed = &digests[case];
        assert_eq!(executed.status, "0", "{case}");
        let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
        let mut features = overlapping_v4(case.trim_end_matches("_dbg"));
        algorithm
            .run(
                ffc1_input(),
                &mut features,
                &ffc1_debug_parameters(),
                &FeatureMap::new(),
            )
            .unwrap();
        assert!(algorithm.termination().is_none(), "{case}");
        assert_eq!(executed.reported("run features "), Some(features.len()));
        let mut disk = Disk::default();
        disk.write_run(algorithm.debug_output(), None);
        disk.assert_executed(executed, case);
    }

    // The abort map of a reused object.
    for (case, scaled) in [("stale_oob", false), ("stale_scaled", true)] {
        let executed = &digests[case];
        let first_parameters = with_debug(
            ffc1_parameters(),
            &[("feature:min_isotope_fit", ParamValue::Float(1.0))],
        );
        let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
        let mut disk = Disk::default();
        algorithm
            .run(
                ffc1_input(),
                &mut FeatureMap::new(),
                &first_parameters,
                &FeatureMap::new(),
            )
            .unwrap();
        let first = algorithm.debug_output().unwrap().clone();
        disk.write_run(Some(&first), None);
        let (input, parameters) = if scaled {
            (scaled_ffc1_input(), first_parameters)
        } else {
            (short_input(), with_debug(Param::new(), &[]))
        };
        let result = algorithm.run(
            input,
            &mut FeatureMap::new(),
            &parameters,
            &FeatureMap::new(),
        );
        let stream = algorithm.debug_log_file().unwrap();
        assert_eq!(stream.written, first.log.len(), "{case}");
        assert_eq!(stream.flushed, first.log.flushed_bytes(), "{case}");
        assert!(stream.failed, "{case}");
        let termination = algorithm.termination();
        if scaled {
            assert_eq!(executed.status, "0", "{case}");
            result.unwrap();
            assert!(termination.is_none(), "{case}");
        } else {
            assert_eq!(executed.status, "139", "{case}: SIGSEGV");
            let error = result.unwrap_err();
            let termination = termination.unwrap();
            assert_eq!(termination.point, TerminationPoint::AbortMap { entry: 0 });
            assert_eq!(termination.kind, TerminationKind::OutOfBounds);
            assert_eq!(termination.message, error.to_string());
            assert_eq!(termination.log_file_bytes, Some(first.log.flushed_bytes()));
            assert_eq!(first.log.flushed_bytes(), 1_114_578);
        }
        disk.write_run(algorithm.debug_output(), termination);
        disk.assert_executed(executed, case);
    }

    // The score arrays, before `debug/` exists.
    let executed = &digests["wrap42_fresh_dbg"];
    assert_eq!(executed.status, "139", "SIGSEGV");
    assert!(executed.files.is_empty());
    let mut parameters = ffc1_debug_parameters();
    set(
        &mut parameters,
        "isotopic_pattern:charge_low",
        ParamValue::Integer(4),
    );
    set(
        &mut parameters,
        "isotopic_pattern:charge_high",
        ParamValue::Integer(2),
    );
    let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
    let error = algorithm
        .run(
            ffc1_input(),
            &mut FeatureMap::new(),
            &parameters,
            &FeatureMap::new(),
        )
        .unwrap_err();
    let termination = algorithm.termination().unwrap();
    assert_eq!(termination.point, TerminationPoint::ScoreArrays);
    assert_eq!(termination.kind, TerminationKind::OutOfBounds);
    assert_eq!(termination.exception, "SIGSEGV");
    assert_eq!(termination.message, error.to_string());
    assert_eq!(termination.log_file_bytes, None);
    assert!(algorithm.debug_output().is_none());
    assert!(algorithm.debug_log_file().is_none());
}

/// Executed `fix5_driver reuse` runs (`../oracle/ffap-complete-fix5`): one
/// object runs FeatureFinderCentroided_1 with `write_debug`, which returns and
/// leaves `debug/log.txt` open with an unflushed tail (the driver measured the
/// file at the flushed length), then a second run, with or without
/// `write_debug`, that ends the process or returns:
///
/// - a charge-0 caller feature (SIGFPE in step 4). After a plain first run
///   the reused windows change the second run's features and nothing traps;
///   after a first run with a NaN `intensity_percentage_optional`, whose
///   empty windows leave the second run a fresh object's windows, it traps;
/// - an empty best isotope pattern (`vfi2_driver neg`'s `avg0` section:
///   SIGSEGV in the seed loop);
/// - a score-array count that wraps to 1 (4/2) or to 1003 arrays
///   (`INT_MAX`/498): SIGSEGV before the stream is touched; one that wraps to
///   `2^32 - 5` arrays (7/2) throws `std::bad_alloc` under the 16 GB address
///   space, which the driver catches;
/// - the same run again, which returns.
///
/// Wherever the process ended, `debug/log.txt` is the first run's flushed
/// prefix, whatever the second run is; where it returned, the destroyed
/// object has written the complete first-run log. The port records every
/// termination with the first run's stream length, and the file protocol
/// reproduces every executed file.
#[test]
fn a_reused_instance_leaves_the_executed_log_at_every_later_termination() {
    let digests = termination_digests();
    let cases: Vec<(&str, &str, bool, bool)> = vec![
        ("reuse_zero_charge_dbg", "zero_charge", true, false),
        ("reuse_zero_charge_nd", "zero_charge", false, false),
        ("reuse_nanipo_zero_charge_dbg", "zero_charge", true, true),
        ("reuse_nanipo_zero_charge_nd", "zero_charge", false, true),
        ("reuse_nanipo_ok_dbg", "ok", true, true),
        ("reuse_nanipo_ok_nd", "ok", false, true),
        ("reuse_avg0_dbg", "avg0", true, false),
        ("reuse_avg0_nd", "avg0", false, false),
        ("reuse_wrap42_dbg", "wrap42", true, false),
        ("reuse_wrap42_nd", "wrap42", false, false),
        ("reuse_wrapmax498_dbg", "wrapmax498", true, false),
        ("reuse_wrapmax498_nd", "wrapmax498", false, false),
        ("reuse_wrap72_dbg", "wrap72", true, false),
        ("reuse_wrap72_nd", "wrap72", false, false),
        ("reuse_ok_dbg", "ok", true, false),
        ("reuse_ok_nd", "ok", false, false),
    ];
    assert_eq!(
        digests
            .keys()
            .filter(|case| case.starts_with("reuse_"))
            .count(),
        cases.len()
    );
    let replay = |&(case, second, debug, nan_first): &(&str, &str, bool, bool)| {
        let executed = &digests[case];
        let mut first_parameters = ffc1_debug_parameters();
        if nan_first {
            set(
                &mut first_parameters,
                "isotopic_pattern:intensity_percentage_optional",
                ParamValue::Float(f64::NAN),
            );
        }
        let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
        let mut features = FeatureMap::new();
        algorithm
            .run(
                ffc1_input(),
                &mut features,
                &first_parameters,
                &FeatureMap::new(),
            )
            .unwrap();
        assert_eq!(
            executed.reported("run1 features "),
            Some(features.len()),
            "{case}"
        );
        let first = algorithm.debug_output().unwrap().clone();
        assert!(first.log_opened && first.termination.is_none(), "{case}");
        assert_eq!(
            executed.reported("run1 file_bytes "),
            Some(first.log.flushed_bytes()),
            "{case}: the file while the object lives"
        );
        let mut disk = Disk::default();
        disk.write_run(Some(&first), None);

        let (input, mut parameters, mut features) = match second {
            "zero_charge" => (ffc1_input(), ffc1_parameters(), overlapping_v4(second)),
            "avg0" => {
                let (input, parameters) = negated_band(1.0, 0.0, 1.0, 0.0, false);
                (input, parameters, FeatureMap::new())
            }
            "wrap42" | "wrapmax498" | "wrap72" => {
                let (low, high) = match second {
                    "wrap42" => (4, 2),
                    "wrap72" => (7, 2),
                    _ => (i64::from(i32::MAX), 498),
                };
                let mut parameters = ffc1_parameters();
                set(
                    &mut parameters,
                    "isotopic_pattern:charge_low",
                    ParamValue::Integer(low),
                );
                set(
                    &mut parameters,
                    "isotopic_pattern:charge_high",
                    ParamValue::Integer(high),
                );
                (ffc1_input(), parameters, FeatureMap::new())
            }
            _ => (ffc1_input(), ffc1_parameters(), FeatureMap::new()),
        };
        if debug {
            parameters = with_debug(
                parameters,
                &[("debug:pseudo_rt_shift", ParamValue::Float(500.0))],
            );
        }
        let result = algorithm.run(input, &mut features, &parameters, &FeatureMap::new());
        // A debug run fails to reopen the stream, unless it ends before the
        // source's `log_.open` (the score arrays come first).
        let reopened = debug && !second.starts_with("wrap");
        let stream = algorithm.debug_log_file().unwrap();
        assert_eq!(
            (stream.written, stream.flushed, stream.failed),
            (first.log.len(), first.log.flushed_bytes(), reopened),
            "{case}: the stream"
        );
        let out = algorithm.debug_output();
        let termination = algorithm.termination();
        if let Some(out) = out {
            assert!(!out.log_opened && out.log.is_empty(), "{case}");
            assert_eq!(out.termination.as_ref(), termination, "{case}");
        }
        if executed.printed("run2 returned") {
            assert_eq!(executed.status, "0", "{case}");
            result.unwrap_or_else(|error| panic!("{case}: {error}"));
            assert!(termination.is_none(), "{case}");
            assert_eq!(
                executed.reported("run2 features "),
                Some(features.len()),
                "{case}"
            );
        } else if executed.printed("run2 threw St9bad_alloc what=std::bad_alloc") {
            // Caught by the driver: no termination, and the destroyed object
            // writes the complete log.
            assert_eq!(executed.status, "0", "{case}");
            assert!(result.is_err(), "{case}");
            assert!(termination.is_none(), "{case}");
        } else {
            assert!(!executed.printed("run2 returned"), "{case}");
            let error = result.unwrap_err();
            let termination = termination.unwrap_or_else(|| panic!("{case}: {error}"));
            let (status, kind, exception) = match second {
                "zero_charge" => ("136", TerminationKind::ArithmeticTrap, "SIGFPE"),
                _ => ("139", TerminationKind::OutOfBounds, "SIGSEGV"),
            };
            assert_eq!(executed.status, status, "{case}");
            assert_eq!(termination.kind, kind, "{case}");
            assert_eq!(termination.exception, exception, "{case}");
            assert_eq!(termination.message, error.to_string(), "{case}");
            assert_eq!(
                termination.log_file_bytes,
                Some(first.log.flushed_bytes()),
                "{case}"
            );
            match second {
                "zero_charge" => assert!(
                    matches!(
                        termination.point,
                        TerminationPoint::OverlapResolution { .. }
                    ),
                    "{case}"
                ),
                "avg0" => assert!(
                    matches!(
                        termination.point,
                        TerminationPoint::Seed { plot_nr: -1, .. }
                    ),
                    "{case}: {:?}",
                    termination.point
                ),
                _ => assert_eq!(termination.point, TerminationPoint::ScoreArrays, "{case}"),
            }
        }
        if executed.status == "0" {
            assert_eq!(
                executed.reported("run2 file_bytes "),
                Some(first.log.flushed_bytes()),
                "{case}"
            );
            assert_eq!(
                executed.reported("destroyed file_bytes "),
                Some(first.log.len()),
                "{case}"
            );
        }
        disk.write_run(out, termination);
        disk.assert_executed(executed, case);
    };
    std::thread::scope(|scope| {
        let replay = &replay;
        let handles: Vec<_> = cases
            .chunks(4)
            .map(|chunk| scope.spawn(move || chunk.iter().for_each(replay)))
            .collect();
        for handle in handles {
            handle
                .join()
                .unwrap_or_else(|panic| std::panic::resume_unwind(panic));
        }
    });
}

// ---------------------------------------------------------------------------
// The complete progress sequence (second driver, counting clock)
// ---------------------------------------------------------------------------

/// A progress backend that records every call with its raw arguments, in the
/// driver's `Recorder` format.
struct Recorder(Shared);

impl ProgressBackend for Recorder {
    fn start_progress(
        &mut self,
        begin: i64,
        end: i64,
        label: &str,
        depth: usize,
    ) -> openms::Result<()> {
        use std::io::Write;
        writeln!(self.0, "S {begin} {end} {depth} {label}")?;
        Ok(())
    }
    fn set_progress(&mut self, value: i64, depth: usize) -> openms::Result<()> {
        use std::io::Write;
        writeln!(self.0, "V {value} {depth}")?;
        Ok(())
    }
    fn next_progress(&mut self) -> openms::Result<i64> {
        use std::io::Write;
        writeln!(self.0, "N")?;
        Ok(0)
    }
    fn end_progress(&mut self, depth: usize, bytes_processed: u64) -> openms::Result<()> {
        use std::io::Write;
        writeln!(self.0, "E {depth} {bytes_processed}")?;
        Ok(())
    }
}

/// The driver's caller map of the `caller` progress case.
fn progress_caller_map() -> FeatureMap {
    let feature = |rt: f64, mz: f64, intensity: f32, charge: i32, quality: f32| {
        let mut f = Feature::new(rt, mz, intensity);
        f.charge = charge;
        f.quality = quality;
        f
    };
    let mut a = feature(4407.0, 646.24, 1.0e7, 2, 0.9);
    a.convex_hulls
        .push(hull(&[(4374.19, 646.229), (4443.42, 646.2585)]));
    let mut b = feature(4301.0, 651.75, 5.0e4, 3, 0.99);
    b.convex_hulls
        .push(hull(&[(4280.0, 651.75), (4320.0, 651.77)]));
    let c = feature(1500.0, 500.0, 1000.0, 1, 0.5);
    FeatureMap::from_features(vec![a, b, c])
}

/// The `FeatureXMLHandler::store()` line of a debug seed or abort map: the
/// map's own id and every feature's that is 0 (`FeatureXMLFile.cpp:80-88`).
/// The source draws the abort map's id before it stores it.
fn store_line(map: &FeatureMap, map_id_drawn: bool) -> Option<String> {
    let map_invalid = !map_id_drawn && map.unique_id == 0;
    let invalid =
        usize::from(map_invalid) + map.features.iter().filter(|f| f.unique_id == 0).count();
    (invalid > 0)
        .then(|| format!("FeatureXMLHandler::store():  found {invalid} invalid unique ids"))
}

/// The executed driver `ffap_progress_driver progress <case>`: its time() is a
/// counter, so `ProgressLogger::setProgress` forwards every call, and a
/// recording backend prints each call with its raw arguments. With a counting
/// clock the port's `ProgressLogger` forwards every call too, and every event
/// and every `std::cout` line matches in order, and the `OPENMS_LOG_INFO`
/// lines in their own order. The inverted ranges the source starts in steps 2
/// and 3.2 on fewer than `2 * min_spectra` scans (case `short`: `S 5 0`, once
/// for step 2 and once per charge for step 3.2) are passed unchanged.
#[test]
fn the_progress_event_sequence_matches_the_release_build() {
    for case in ["ffc1", "defaults", "short", "caller", "debug"] {
        let executed = fixture(&format!("progress_events_{case}.txt.gz"));
        let mut events = Vec::new();
        let mut info = Vec::new();
        let mut in_case = false;
        for line in executed.lines() {
            if line.starts_with("== CASE ") {
                in_case = true;
            } else if line.starts_with("== END ") {
                in_case = false;
            } else if let Some(calls) = line.strip_prefix("TIMECALLS ") {
                let calls: u64 = calls.parse().unwrap();
                assert!(calls > 100, "{case}: time() was not interposed");
            } else if in_case
                && (line.starts_with("S ")
                    || line.starts_with("V ")
                    || line.starts_with("N")
                    || line.starts_with("E ")
                    || line.starts_with("Found "))
            {
                events.push(line.to_owned());
            } else {
                info.push(line.to_owned());
            }
        }

        let (experiment, parameters, mut features) = match case {
            "ffc1" => (ffc1_input(), ffc1_parameters(), FeatureMap::new()),
            "defaults" => (ffc1_input(), Param::new(), FeatureMap::new()),
            "short" => (short_input(), Param::new(), FeatureMap::new()),
            "caller" => (ffc1_input(), ffc1_parameters(), progress_caller_map()),
            _ => (
                ffc1_input(),
                with_debug(
                    ffc1_parameters(),
                    &[("feature:min_isotope_fit", ParamValue::Float(1.0))],
                ),
                FeatureMap::new(),
            ),
        };
        let shared = Shared::default();
        let ticks = Arc::new(std::sync::atomic::AtomicI64::new(0));
        let clock_ticks = ticks.clone();
        let mut logger = ProgressLogger::with_clock_and_nesting(
            Arc::new(move || {
                Ok(ProgressTime {
                    wall_second: clock_ticks.fetch_add(1, std::sync::atomic::Ordering::Relaxed),
                    wall_seconds: 0.0,
                    cpu_seconds: None,
                })
            }),
            ProgressNesting::default(),
        );
        logger.set_logger(Box::new(Recorder(shared.clone())));
        let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
        algorithm.set_progress_logger(Some(logger));
        let sink = shared.clone();
        algorithm.set_console(Some(Box::new(move |line: &ReportLine| {
            if let ReportLine::Out(text) = line {
                use std::io::Write;
                writeln!(sink.clone(), "{text}").unwrap();
            }
        })));
        algorithm
            .run(experiment, &mut features, &parameters, &FeatureMap::new())
            .unwrap();
        let produced = String::from_utf8(shared.0.lock().unwrap().clone()).unwrap();
        let produced: Vec<&str> = produced.lines().collect();
        assert_eq!(produced.len(), events.len(), "{case}: event count");
        for (index, (ours, theirs)) in produced.iter().zip(&events).enumerate() {
            assert_eq!(ours, theirs, "{case} event {index}");
        }
        assert_eq!(produced, events, "{case}: events");
        // `S <begin> <end> <depth> <label>` with begin > end: the executed
        // `short` case has five (`S 5 0`: step 2, and step 3.2 for each of the
        // four default charges), every other case none.
        let inverted: Vec<&str> = produced
            .iter()
            .copied()
            .filter(|line| {
                line.strip_prefix("S ").is_some_and(|rest| {
                    let mut fields = rest.splitn(3, ' ');
                    let begin: i64 = fields.next().unwrap().parse().unwrap();
                    let end: i64 = fields.next().unwrap().parse().unwrap();
                    begin > end
                })
            })
            .collect();
        if case == "short" {
            assert_eq!(inverted.len(), 5, "{case}: inverted ranges");
            assert!(
                inverted.iter().all(|line| line.starts_with("S 5 0 0 ")),
                "{case}: {inverted:?}"
            );
        } else {
            assert!(inverted.is_empty(), "{case}: {inverted:?}");
        }

        // The OPENMS_LOG_INFO lines, with the store lines of a debug run.
        let debug = algorithm.debug_output();
        let mut ours_info = Vec::new();
        for line in algorithm.report() {
            match line {
                ReportLine::Info(text) => ours_info.push(text.clone()),
                ReportLine::StoreSeedMap(index) => {
                    let map = &debug.unwrap().seed_maps[*index].map;
                    ours_info.extend(store_line(map, false));
                }
                ReportLine::StoreAbortReasons => {
                    let map = debug.unwrap().abort_reasons.as_ref().unwrap();
                    ours_info.extend(store_line(map, true));
                }
                ReportLine::Out(_) | ReportLine::Warn(_) | ReportLine::StoreInput => {}
            }
        }
        while info.last().is_some_and(String::is_empty) {
            info.pop();
        }
        assert_eq!(ours_info, info, "{case}: info lines");
        assert!(
            algorithm
                .report()
                .iter()
                .all(|l| !matches!(l, ReportLine::Warn(_)))
        );
    }
}

// ---------------------------------------------------------------------------
// A string or list under debug:pseudo_rt_shift (second driver)
// ---------------------------------------------------------------------------

/// `ffap_progress_driver pun_value`: the `double` of a string and a list
/// `ParamValue` is the bit pattern of a heap pointer, a subnormal number below
/// `2^-1027`, different in each of three processes. The pointers lie below
/// `DEFAULT_MAP_WINDOW = (1 << 47) - 4096` of the reference host's kernel
/// headers (`shift_band/address_space.txt`), the port's bound.
#[test]
fn a_string_shift_is_a_heap_address_in_the_release_build() {
    assert_eq!(HEAP_ADDRESS_END, (1 << 47) - 4096);
    let values = fixture("pun_values.txt");
    let mut seen = std::collections::BTreeSet::new();
    for line in values.lines() {
        let bits = u64::from_str_radix(line.rsplit(' ').next().unwrap(), 16).unwrap();
        let value = f64::from_bits(bits);
        // `2^-1027`, the subnormal of integer `2^47`.
        assert!(value > 0.0 && value < f64::from_bits(1 << 47), "{line}");
        assert!(bits < HEAP_ADDRESS_END, "{line}");
        seen.insert(bits);
    }
    assert_eq!(seen.len(), 9, "every value differs");
}

/// `ffap_progress_driver pun string` and `pun list`: FeatureFinderCentroided_1
/// with `debug:pseudo_rt_shift` a string or a string list. Every retention
/// time is far from zero, so the address vanishes in `rt + k * shift`: the 75
/// feature files and the log were byte-identical in all six processes, and the
/// port writes the same bytes. The feature map is the plain run's.
#[test]
fn a_heap_address_shift_matches_the_release_build_where_it_is_reproducible() {
    let digests = debug_digests();
    for value in [
        ParamValue::String("500".into()),
        ParamValue::StringList(vec!["500".into()]),
    ] {
        let parameters = with_debug(ffc1_parameters(), &[("debug:pseudo_rt_shift", value)]);
        let (result, algorithm, features) =
            debug_run(ffc1_input(), &parameters, Options::default());
        result.unwrap();
        let out = algorithm.debug_output().unwrap();
        assert!(out.termination.is_none());
        assert_executed_bytes("pun", "log.txt", out.log.text().as_bytes());
        let mut produced = 0;
        for files in &out.feature_files {
            let name = |full: String| full.trim_start_matches("debug/").to_owned();
            assert_executed_bytes("pun", &name(files.dta_name()), files.dta.as_bytes());
            assert_executed_bytes("pun", &name(files.plot_name()), &files.plot);
            produced += 2;
            if let Some(cropped) = &files.cropped_dta {
                assert_executed_bytes("pun", &name(files.cropped_dta_name()), cropped.as_bytes());
                produced += 1;
            }
        }
        let expected = digests
            .keys()
            .filter(|(case, file)| case == "pun" && file.starts_with("features/"))
            .count();
        assert_eq!(produced, expected);
        assert_dumps_match(
            &fixture("reuse_run1.txt.gz"),
            &dump_map(&features, algorithm.aborts()),
            "pun",
        );
    }
}

/// FeatureFinderCentroided_1 as `vfi2_driver neg` loads it: MS1, every
/// intensity kept.
fn raw_ffc1_input() -> MSExperiment {
    let mut options = PeakFileOptions::default();
    options.add_ms_level(1).unwrap();
    FileHandler::load_experiment_with_options(
        repository("tests/data/mzml_mobility/FeatureFinderCentroided_1_input.mzML"),
        &[FileType::MzMl],
        &options,
    )
    .unwrap()
}

/// `vfi2_driver neg`'s input and parameters: every intensity whose m/z lies in
/// `[lo, hi]` multiplied by `-factor` (a `float` times a `double`, narrowed),
/// the FFC_1 algorithm section with `reported_mz` average, `rt_shape`
/// symmetric, the four feature thresholds 0 and `seed:min_score`, and for a
/// debug run `write_debug` with `debug:pseudo_rt_shift` 500.
fn negated_band(
    lo: f64,
    hi: f64,
    factor: f64,
    seed_min: f64,
    debug: bool,
) -> (MSExperiment, Param) {
    let mut experiment = raw_ffc1_input();
    for spectrum in &mut experiment.spectra {
        for peak in &mut spectrum.peaks {
            if peak.mz >= lo && peak.mz <= hi {
                peak.intensity = (f64::from(peak.intensity) * -factor) as f32;
            }
        }
    }
    let mut parameters = ffc1_parameters();
    set(
        &mut parameters,
        "feature:reported_mz",
        ParamValue::String("average".into()),
    );
    set(
        &mut parameters,
        "feature:rt_shape",
        ParamValue::String("symmetric".into()),
    );
    for key in [
        "feature:min_score",
        "feature:min_isotope_fit",
        "feature:min_trace_score",
        "feature:min_rt_span",
    ] {
        set(&mut parameters, key, ParamValue::Float(0.0));
    }
    set(
        &mut parameters,
        "seed:min_score",
        ParamValue::Float(seed_min),
    );
    if debug {
        parameters = with_debug(
            parameters,
            &[("debug:pseudo_rt_shift", ParamValue::Float(500.0))],
        );
    }
    (experiment, parameters)
}

/// Executed `vfi2_driver neg` runs (`../oracle/ffap-complete-fix3`, after
/// `../oracle/ffc-instrumentation-v2`): with `feature:min_isotope_fit` 0 a
/// seed whose best isotope pattern stayed empty reaches `extendMassTraces_`,
/// which reads the pattern's first entry, and the Release build dies with
/// SIGSEGV, with and without `write_debug`, two runs of each, identical. The
/// three inputs negate an m/z band around 646.74 with `seed:min_score` 0.3 or
/// 0.35, or nothing with `seed:min_score` 0.
///
/// The port refuses at that seed (lead decision D1) and records the
/// termination. The executed `debug/log.txt` is exactly the prefix of the
/// port's log that its model of the file buffer had flushed when the seed was
/// refused, the refused seed's own lines counted, and the executed feature
/// files are the port's, byte for byte, no more and no fewer.
#[test]
fn a_seed_loop_crash_keeps_what_the_executed_process_had_written() {
    let digests = debug_digests();
    for (case, lo, hi, factor, seed_min, plots) in [
        ("neg_oob1", 646.72686, 646.75686, 0.05, 0.3, 53),
        ("neg_oob_seed035", 646.72686, 646.75686, 0.05, 0.35, 51),
        ("neg_none_avg0", 1.0, 0.0, 1.0, 0.0, 26),
    ] {
        let (experiment, parameters) = negated_band(lo, hi, factor, seed_min, false);
        let (result, algorithm, _) = debug_run(experiment, &parameters, Options::default());
        let refusal = result.unwrap_err().to_string();
        assert!(
            refusal.contains("the isotope pattern matched no peak"),
            "{case}: {refusal}"
        );
        assert!(algorithm.debug_output().is_none(), "{case}");

        let (experiment, parameters) = negated_band(lo, hi, factor, seed_min, true);
        let (result, algorithm, _) = debug_run(experiment, &parameters, Options::default());
        assert_eq!(result.unwrap_err().to_string(), refusal, "{case}");
        let out = algorithm.debug_output().unwrap();
        let termination = out.termination.as_ref().unwrap();
        assert_eq!(termination.kind, TerminationKind::OutOfBounds, "{case}");
        assert_eq!(termination.exception, "SIGSEGV", "{case}");
        assert!(
            matches!(
                termination.point,
                TerminationPoint::Seed {
                    charge: 2,
                    plot_nr: -1,
                    ..
                }
            ),
            "{case}: the seed ends before the fit: {:?}",
            termination.point
        );
        assert_eq!(termination.message, refusal, "{case}");
        assert_eq!(
            termination.log_file_bytes,
            Some(out.log.flushed_bytes()),
            "{case}"
        );
        let log = out.log.text().as_bytes();
        let flushed = out.log.flushed_bytes();
        assert!(
            flushed < log.len(),
            "{case}: the refused seed's lines are buffered"
        );
        assert_executed_bytes(case, "log.txt", &log[..flushed]);
        let mut produced = 0;
        let mut plot_nrs = std::collections::BTreeSet::new();
        for files in &out.feature_files {
            let name = |full: String| full.trim_start_matches("debug/").to_owned();
            assert_executed_bytes(case, &name(files.dta_name()), files.dta.as_bytes());
            assert_executed_bytes(case, &name(files.plot_name()), &files.plot);
            produced += 2;
            if let Some(cropped) = &files.cropped_dta {
                assert_executed_bytes(case, &name(files.cropped_dta_name()), cropped.as_bytes());
                produced += 1;
            }
            plot_nrs.insert(files.plot_nr);
        }
        let expected = digests
            .keys()
            .filter(|(c, file)| c == case && file.starts_with("features/"))
            .count();
        assert_eq!(produced, expected, "{case}: feature files");
        assert_eq!(plot_nrs.len(), plots, "{case}: plots");
        assert_eq!(plot_nrs.last().copied(), Some(plots as i64 - 1), "{case}");
    }
}

/// Executed file sizes and SHA-1 digests, by case and file name.
type FileDigests = BTreeMap<(String, String), (usize, String)>;

/// The executed digests of `never_returns_digests.tsv`, by case and file,
/// and each case's exit status.
fn never_returns_digests() -> (FileDigests, BTreeMap<String, String>) {
    let mut files = BTreeMap::new();
    let mut statuses = BTreeMap::new();
    for line in fixture("never_returns_digests.tsv").lines().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields[1] == "status" {
            statuses.insert(fields[0].to_owned(), fields[3].to_owned());
        } else {
            files.insert(
                (fields[0].to_owned(), fields[1].to_owned()),
                (fields[2].parse().unwrap(), fields[3].to_owned()),
            );
        }
    }
    (files, statuses)
}

/// Executed `fix4_vfi` runs (`../oracle/ffap-complete-fix4`, the round-3
/// instrumentation verifier's `vfi3_driver`, two runs each, identical;
/// `never_returns_digests.tsv`): FeatureFinderCentroided_1 with the retention
/// time of scan 50 (or 20) set to NaN and `write_debug` with
/// `debug:pseudo_rt_shift` 500. A seed whose mass traces include the NaN
/// scan reaches `MassTraces::computeIntensityProfile` in the fitter's start
/// values, whose merge never ends (`CPP-242`): the Release process was killed
/// after 30 s (status 137), Gaussian and EGH alike, with and without
/// `write_debug`.
///
/// The port refuses at that merge (lead decision D1) and records a
/// `NeverReturns` termination with the plot number the hanging seed had
/// received. The executed `debug/log.txt` is the prefix of the port's log that
/// its model of the file buffer had flushed, the executed seed map is the
/// port's, and the executed feature files (those of the seeds fitted before)
/// are the port's, byte for byte, no more and no fewer.
#[test]
fn a_seed_loop_that_never_returns_keeps_what_the_executed_process_had_written() {
    let (digests, statuses) = never_returns_digests();
    let nan_rt = |scan: usize, egh: bool, debug: bool| {
        // `vfi3_driver filt`: MS1 with the intensity range
        // `[numeric_limits<DPosition<1>>::min(), maxPositive()]`, which is
        // `[0, DBL_MAX]` (`DPosition<1>()` is 0), as `ffc1_input` loads it.
        let mut experiment = ffc1_input();
        experiment.spectra[scan].rt = f64::NAN;
        let mut parameters = ffc1_parameters();
        if egh {
            set(
                &mut parameters,
                "feature:rt_shape",
                ParamValue::String("asymmetric".into()),
            );
        }
        if debug {
            parameters = with_debug(
                parameters,
                &[("debug:pseudo_rt_shift", ParamValue::Float(500.0))],
            );
        }
        (experiment, parameters)
    };

    assert_eq!(statuses["nanrt50_nd"], "137");
    let (experiment, parameters) = nan_rt(50, false, false);
    let (result, algorithm, _) = debug_run(experiment, &parameters, Options::default());
    let refusal = result.unwrap_err().to_string();
    assert!(
        refusal.contains("a NaN retention time cannot be merged into an intensity profile"),
        "{refusal}"
    );
    assert!(algorithm.debug_output().is_none());

    for (case, scan, egh, plots) in [
        ("nanrt50_dbg", 50, false, 3),
        ("nanrt20_dbg", 20, false, 10),
        ("nanrt50_dbg_egh", 50, true, 3),
    ] {
        assert_eq!(statuses[case], "137", "{case}: killed after 30 s");
        let (experiment, parameters) = nan_rt(scan, egh, true);
        let (result, algorithm, _) = debug_run(experiment, &parameters, Options::default());
        let refusal = result.unwrap_err().to_string();
        assert!(
            refusal.contains("a NaN retention time cannot be merged into an intensity profile"),
            "{case}: {refusal}"
        );
        let out = algorithm.debug_output().unwrap();
        let termination = out.termination.as_ref().unwrap();
        assert_eq!(termination.kind, TerminationKind::NeverReturns, "{case}");
        assert_eq!(termination.exception, "", "{case}");
        assert!(
            matches!(
                termination.point,
                TerminationPoint::Seed { charge: 2, plot_nr, .. } if plot_nr == plots as i64
            ),
            "{case}: the hanging seed's plot number: {:?}",
            termination.point
        );
        assert_eq!(termination.message, refusal, "{case}");
        assert_eq!(algorithm.termination(), Some(termination), "{case}");
        let log = out.log.text().as_bytes();
        let flushed = out.log.flushed_bytes();
        let (bytes, sha1) = &digests[&(case.to_owned(), "log.txt".to_owned())];
        assert_eq!(flushed, *bytes, "{case}: flushed log bytes");
        assert_eq!(termination.log_file_bytes, Some(*bytes), "{case}");
        // The killed process had stored the charge's seed map
        // (`debug/seeds_2.featureXML`, `../oracle/ffap-complete-fix4`, both
        // runs identical): the algorithm's seeds, decoded under D6.
        assert_eq!(out.seed_maps.len(), 1, "{case}");
        assert_eq!(out.seed_maps[0].charge, 2, "{case}");
        assert_maps_decoded_equal(
            &out.seed_maps[0].map,
            &feature_fixture(&format!("never_returns_seed_map_{case}.featureXML.gz")),
            case,
        );
        assert!(
            flushed < log.len(),
            "{case}: the hanging seed's lines are buffered"
        );
        assert_eq!(&sha1_hex(&log[..flushed]), sha1, "{case}: log.txt");
        let mut produced = 0;
        let mut plot_nrs = std::collections::BTreeSet::new();
        let check = |file: String, data: &[u8]| {
            let name = file.trim_start_matches("debug/").to_owned();
            let (bytes, sha1) = digests
                .get(&(case.to_owned(), name.clone()))
                .unwrap_or_else(|| panic!("{case}: {name} was not written by the executed run"));
            assert_eq!(data.len(), *bytes, "{case} {name}: size");
            assert_eq!(&sha1_hex(data), sha1, "{case} {name}: content");
        };
        for files in &out.feature_files {
            check(files.dta_name(), files.dta.as_bytes());
            check(files.plot_name(), &files.plot);
            produced += 2;
            if let Some(cropped) = &files.cropped_dta {
                check(files.cropped_dta_name(), cropped.as_bytes());
                produced += 1;
            }
            plot_nrs.insert(files.plot_nr);
        }
        let expected = digests
            .keys()
            .filter(|(c, file)| c == case && file.starts_with("features/"))
            .count();
        assert_eq!(produced, expected, "{case}: feature files");
        assert_eq!(plot_nrs.len(), plots, "{case}: plots");
        assert_eq!(plot_nrs.last().copied(), Some(plots as i64 - 1), "{case}");
    }
}

/// `fix3_driver progress` (`../oracle/ffap-complete-fix3`, two runs each,
/// identical): step 1 calls `startProgress(0, intensity_bins_ *
/// intensity_bins_)` with a `UInt` product, which wraps modulo 2^32 (65,536
/// bins give `S 0 0`), and the driver ends the process right after that event.
/// The port logs the same start event (lead decision D12: an in-bounds wrap is
/// reproduced). Its intensity-bin and work ceilings are raised or lowered so
/// that the run fails right after the event instead of computing 2^32 cells;
/// `intensity:bins` 2^32 narrows to 0, which the parameter check refuses with
/// the executed text.
#[test]
fn the_step_one_progress_range_wraps_as_the_release_build_computes_it() {
    for line in fixture("progress_bins.tsv").lines().skip(1) {
        let (bins, outcome) = line.split_once('\t').unwrap();
        let bins: i64 = bins.parse().unwrap();
        let mut parameters = ffc1_parameters();
        set(&mut parameters, "intensity:bins", ParamValue::Integer(bins));
        let shared = Shared::default();
        // A fresh nesting depth per case, as each executed case is a process
        // of its own; a failed run leaves its section open.
        let mut logger = ProgressLogger::with_clock_and_nesting(
            Arc::new(|| {
                Ok(ProgressTime {
                    wall_second: 0,
                    wall_seconds: 0.0,
                    cpu_seconds: None,
                })
            }),
            ProgressNesting::default(),
        );
        logger.set_logger(Box::new(Recorder(shared.clone())));
        let mut algorithm = FeatureFinderAlgorithmPicked::with_options(Options {
            limits: Limits {
                max_intensity_bins: 1 << 17,
                max_work: 1,
                ..Limits::default()
            },
            ..Options::default()
        })
        .unwrap();
        algorithm.set_progress_logger(Some(logger));
        let result = algorithm.run(
            ffc1_input(),
            &mut FeatureMap::new(),
            &parameters,
            &FeatureMap::new(),
        );
        if let Some(what) = outcome.strip_prefix("threw InvalidParameter: ") {
            assert_eq!(
                result.unwrap_err().to_string(),
                format!("invalid value: {what}"),
                "{bins}"
            );
            continue;
        }
        let events = String::from_utf8(shared.0.lock().unwrap().clone()).unwrap();
        let start = events
            .lines()
            .find(|line| line.ends_with(" Precalculating intensity scores"))
            .unwrap_or_else(|| panic!("{bins}: no step-1 event in {events}"));
        assert_eq!(start, outcome, "{bins}");
        if bins > 2000 {
            // The raised bin ceiling lets the run reach step 1; the work
            // ceiling ends it there.
            assert!(result.is_err(), "{bins}");
        }
    }
}

/// `ffap_progress_driver pun zero_rt`: the same input with every retention
/// time reduced by 4404.89 s, so the scan of the first seed that reaches the
/// fit sits at RT 0 and `0 + k * shift` prints the address: `0.dta` differed in
/// each of three processes (while the log and the seed map did not). The port
/// refuses at that seed's files, after the seed map and the log up to that
/// point, which match the executed ones.
#[test]
fn a_heap_address_shift_is_refused_where_the_written_text_depends_on_the_address() {
    let digests = debug_digests();
    let executed: std::collections::BTreeSet<&String> = (1..=3)
        .map(|rep| &digests[&(format!("zero_rt.{rep}"), "features/0.dta".to_owned())].1)
        .collect();
    assert_eq!(
        executed.len(),
        3,
        "the executed 0.dta differs from process to process"
    );

    let mut experiment = ffc1_input();
    for spectrum in &mut experiment.spectra {
        spectrum.rt -= 4404.89;
    }
    let parameters = with_debug(
        ffc1_parameters(),
        &[("debug:pseudo_rt_shift", ParamValue::String("500".into()))],
    );
    let (result, algorithm, features) = debug_run(experiment, &parameters, Options::default());
    let error = result.unwrap_err();
    assert!(
        matches!(&error, openms::Error::Unsupported(message)
            if message.contains("in 0.dta trace 1") && message.contains("retention time 0e0")),
        "{error}"
    );
    assert!(features.is_empty());
    let out = algorithm.debug_output().unwrap();
    assert!(out.termination.is_none());
    assert!(out.feature_files.is_empty());
    assert_eq!(out.seed_maps.len(), 1);
    assert_maps_decoded_equal(
        &out.seed_maps[0].map,
        &feature_fixture("pun_zero_rt_seed_map_2.featureXML.gz"),
        "zero_rt seeds",
    );
    assert_executed_bytes("zero_rt", "log_to_plot0.txt", out.log.text().as_bytes());
}

/// What `ffap_shift_band_driver pun_rt <value>` printed in each process: the
/// moved retention time's bits, the shift's bits and the feature count, by
/// value.
struct BandRun {
    value: f64,
    shifts: Vec<u64>,
    features: usize,
}

fn band_run(value: &str) -> BandRun {
    let text = fixture("band_stdout.txt");
    let mut values = std::collections::BTreeSet::new();
    let mut shifts = Vec::new();
    let mut features = std::collections::BTreeSet::new();
    for line in text.lines() {
        let (process, rest) = line.split_once(' ').unwrap();
        if process.rsplit_once('.').unwrap().0 != value {
            continue;
        }
        let fields: Vec<&str> = rest.split(' ').collect();
        match fields[0] {
            "moved" => {
                assert_eq!(fields[1], "1", "{line}");
                values.insert(u64::from_str_radix(fields[3], 16).unwrap());
            }
            "shift_bits" => shifts.push(u64::from_str_radix(fields[1], 16).unwrap()),
            "features" => {
                features.insert(fields[1].parse::<usize>().unwrap());
            }
            _ => {}
        }
    }
    assert_eq!(values.len(), 1, "{value}");
    assert_eq!(features.len(), 1, "{value}");
    assert_eq!(shifts.len(), 3, "{value}");
    BandRun {
        value: f64::from_bits(values.into_iter().next().unwrap()),
        shifts,
        features: features.into_iter().next().unwrap(),
    }
}

/// The zero_rt input (every retention time reduced by 4404.89 s) with the one
/// scan at RT 0 moved to `value`, as the driver builds it.
fn band_input(value: f64) -> MSExperiment {
    let mut experiment = ffc1_input();
    let mut moved = 0;
    for spectrum in &mut experiment.spectra {
        spectrum.rt -= 4404.89;
        if spectrum.rt == 0.0 {
            spectrum.rt = value;
            moved += 1;
        }
    }
    assert_eq!(moved, 1);
    experiment
}

/// `ffap_shift_band_driver pun_rt 5e-275` and `pun_rt 1e-289`: the seed's scan
/// at a retention time of that magnitude, `debug:pseudo_rt_shift` a string.
/// Three processes with three different shifts wrote the same log, seed map,
/// input and 75 feature files, because no shift below `HEAP_ADDRESS_END`
/// changes `k * shift + rt` there. The port writes those bytes. `1e-289` lies
/// below the bound an address below `2^56` would give (`1e-289 + 2 * 2^-1007`
/// is not `1e-289`), which would refuse it.
#[test]
fn a_heap_address_shift_matches_the_release_build_down_to_the_address_bound() {
    for value in ["5e-275", "1e-289"] {
        let executed = band_run(value);
        let distinct: std::collections::BTreeSet<&u64> = executed.shifts.iter().collect();
        assert_eq!(distinct.len(), 3, "{value}: three processes, three shifts");
        assert!(executed.shifts.iter().all(|&bits| bits < HEAP_ADDRESS_END));
        let case = format!("band_{value}");
        let parameters = with_debug(
            ffc1_parameters(),
            &[("debug:pseudo_rt_shift", ParamValue::String("500".into()))],
        );
        let (result, algorithm, features) =
            debug_run(band_input(executed.value), &parameters, Options::default());
        result.unwrap();
        assert_eq!(features.len(), executed.features, "{value}");
        let out = algorithm.debug_output().unwrap();
        assert!(out.termination.is_none());
        assert_executed_bytes(&case, "log.txt", out.log.text().as_bytes());
        let mut produced = 0;
        for files in &out.feature_files {
            let name = |full: String| full.trim_start_matches("debug/").to_owned();
            assert_executed_bytes(&case, &name(files.dta_name()), files.dta.as_bytes());
            assert_executed_bytes(&case, &name(files.plot_name()), &files.plot);
            produced += 2;
            if let Some(cropped) = &files.cropped_dta {
                assert_executed_bytes(&case, &name(files.cropped_dta_name()), cropped.as_bytes());
                produced += 1;
            }
        }
        assert_eq!(produced, 75, "{value}");
        assert_eq!(out.seed_maps.len(), 1);
        assert_maps_decoded_equal(
            &out.seed_maps[0].map,
            &feature_fixture(&format!("band_{value}_seed_map_2.featureXML.gz")),
            &format!("{case} seeds"),
        );
        assert_maps_decoded_equal(
            out.abort_reasons.as_ref().unwrap(),
            &feature_fixture(&format!("band_{value}_abort_map.featureXML.gz")),
            &format!("{case} aborts"),
        );
    }
}

/// `ffap_shift_band_driver pun_rt 1e-295` and `pun_rt 1e-300`: there the
/// shift changes the scan's retention time in trace 1 of plot 0, and each of
/// the three processes wrote a different `0.dta` (the log and the seed map
/// were the same). The port refuses at that seed's files, after the seed map
/// and the log up to that point, which match the executed ones.
#[test]
fn a_heap_address_shift_is_refused_below_the_address_bound() {
    let digests = debug_digests();
    for value in ["1e-295", "1e-300"] {
        let executed = band_run(value);
        let files: std::collections::BTreeSet<&String> = (1..=3)
            .map(|rep| &digests[&(format!("band_{value}.{rep}"), "features/0.dta".to_owned())].1)
            .collect();
        assert_eq!(
            files.len(),
            3,
            "{value}: 0.dta differs from process to process"
        );
        let parameters = with_debug(
            ffc1_parameters(),
            &[("debug:pseudo_rt_shift", ParamValue::String("500".into()))],
        );
        let (result, algorithm, features) =
            debug_run(band_input(executed.value), &parameters, Options::default());
        let error = result.unwrap_err();
        let expected = format!("in 0.dta trace 1 has the retention time {value}, ");
        assert!(
            matches!(&error, openms::Error::Unsupported(message) if message.contains(&expected)),
            "{error}"
        );
        assert!(features.is_empty());
        let out = algorithm.debug_output().unwrap();
        assert!(out.termination.is_none());
        assert!(out.feature_files.is_empty());
        assert_eq!(out.seed_maps.len(), 1);
        assert_maps_decoded_equal(
            &out.seed_maps[0].map,
            &feature_fixture(&format!("band_{value}_seed_map_2.featureXML.gz")),
            &format!("band_{value} seeds"),
        );
        assert_executed_bytes(
            &format!("band_{value}"),
            "log_to_plot0.txt",
            out.log.text().as_bytes(),
        );
    }
}

/// The check behind the refusal, on its own, at its boundaries. With
/// `L = (2^47 - 4096) * 2^-1074`, just below `2^-1027`, trace `k` keeps a
/// retention time `v` exactly when `v + k * L` rounds to `v`, which needs half
/// the spacing of the doubles next to `v` (on the side of the shift) above
/// `k * L`:
///
/// - `2^-974` has the spacing `2^-1026` above it, so `2^-974 + L` rounds back
///   (`L < 2^-1027`), and `2 * L` does not; `2^-973` keeps `2 * L`.
/// - The double below `2^-974` has the spacing `2^-1027`, so `L` moves it.
/// - `-2^-974` moves towards zero, where the spacing is `2^-1027`: moved.
///   The double below it, `-(2^-974 + 2^-1026)`, has the spacing `2^-1026`
///   on both sides: kept.
/// - `0`, `-0`, `1e-300` move; NaN and the infinities stay.
///
/// The `.plot` prints the centre with six significant digits: `1e-303 + k * L`
/// prints as `1e-303` up to `k = 7` (`7 * L / 1e-303 < 4.9e-6`) and as
/// `1.00001e-303` from `k = 8` (`5.6e-6`). A centre of `0` changes at `k = 1`.
/// Trace 0 is never shifted. Where every value is kept, the files are those of
/// shift `0`.
#[test]
fn a_heap_address_shift_is_refused_only_where_an_address_changes_the_text() {
    use openms::analysis::feature_finder_picked::gauss_trace_fitter::GaussTraceFitter;
    use openms::analysis::feature_finder_picked::helper_structs::{
        MassTrace, MassTraces, TracePeak,
    };
    let trace = |rt: f64| {
        let mut trace = MassTrace::default();
        trace.peaks.push(TracePeak::new(0, 0, rt, 500.0, 10.0));
        trace.theoretical_int = 1.0;
        trace
    };
    let traces = |rts: &[f64]| {
        let mut traces = MassTraces::new();
        for &rt in rts {
            traces.push(trace(rt));
        }
        traces
    };
    let fitter = |center: f64| {
        let mut fitter = GaussTraceFitter::default();
        fitter.set_optimized_parameters([100.0, center, 2.0]);
        fitter
    };
    let write = |center: f64, extended: &MassTraces, cropped: &MassTraces, shift: PseudoRtShift| {
        let fitter = fitter(center);
        write_feature_debug_info(FeatureDebugInput {
            fitter: &fitter,
            traces: extended,
            new_traces: cropped,
            feature_ok: true,
            error_msg: "",
            final_score: 0.5,
            plot_nr: 3,
            seed_mz: 500.0,
            features_len: 0,
            pseudo_rt_shift: shift,
            path: "debug/features/",
        })
    };
    let kept = |center: f64, extended: &MassTraces, cropped: &MassTraces| {
        let files = write(center, extended, cropped, PseudoRtShift::HeapAddress).unwrap();
        assert_eq!(
            files,
            write(center, extended, cropped, PseudoRtShift::Value(0.0)).unwrap()
        );
    };
    let refused = |center: f64, extended: &MassTraces, cropped: &MassTraces, place: &str| {
        let result = write(center, extended, cropped, PseudoRtShift::HeapAddress);
        assert!(
            matches!(&result, Err(openms::Error::Unsupported(message))
                if message.contains(place) && !message.contains("  ")),
            "{place}: {result:?}"
        );
    };
    let empty = MassTraces::new();
    // `2^-974` and `2^-973`: biased exponents 49 and 50.
    let low = f64::from_bits(49 << 52);
    let high = f64::from_bits(50 << 52);
    let below_low = f64::from_bits(low.to_bits() - 1);
    let beyond_minus_low = f64::from_bits((-low).to_bits() + 1);
    assert!(beyond_minus_low < -low);

    // Trace 0 is never shifted, whatever its values.
    kept(0.0, &traces(&[0.0]), &traces(&[-0.0]));
    // The `.dta` boundaries, with a centre that keeps its text.
    kept(100.0, &traces(&[10.0, low]), &empty);
    kept(100.0, &traces(&[10.0, beyond_minus_low]), &empty);
    kept(100.0, &traces(&[10.0, 10.0, high]), &empty);
    kept(
        100.0,
        &traces(&[10.0, f64::NAN, f64::INFINITY]),
        &traces(&[1.0, f64::NEG_INFINITY]),
    );
    for rt in [below_low, -low, 0.0, -0.0, 1e-300] {
        refused(
            100.0,
            &traces(&[10.0, rt]),
            &empty,
            "3.dta trace 1 has the retention time",
        );
        refused(
            100.0,
            &traces(&[10.0, 20.0]),
            &traces(&[10.0, rt]),
            "3_cropped.dta trace 1",
        );
    }
    refused(100.0, &traces(&[10.0, 10.0, low]), &empty, "3.dta trace 2");
    // The `.plot` centre.
    let eight = traces(&[10.0; 8]);
    kept(1e-303, &eight, &eight);
    refused(
        1e-303,
        &traces(&[10.0; 9]),
        &empty,
        "3.plot trace 8 has the fitted centre 1e-303",
    );
    refused(
        0.0,
        &traces(&[10.0, 10.0]),
        &empty,
        "3.plot trace 1 has the fitted centre 0e0",
    );
    kept(f64::NAN, &traces(&[10.0, 10.0]), &empty);
    // A `.dta` value is checked before the centre.
    refused(0.0, &traces(&[10.0, 0.0]), &empty, "3.dta trace 1");
}

/// The rest of the `DefaultParamHandler` base the source object inherits
/// (source-reviewed: `DefaultParamHandler.cpp:32-40`, `:65`, `:110-123`;
/// `Param.cpp:1085`): `setName` renames the handler in the unknown-parameter
/// warnings, `getSubsections` is empty, and `operator==` compares the handler
/// part only, with the refused set the source keeps in `param_`.
#[test]
fn the_handler_base_renames_compares_and_has_no_subsections() {
    let mut a = FeatureFinderAlgorithmPicked::new().unwrap();
    let mut b = FeatureFinderAlgorithmPicked::new().unwrap();
    assert!(a.subsections().is_empty());
    assert!(a.handler_equal(&b).unwrap());
    // State outside the handler is not compared.
    a.set_seeds(&prefilled());
    assert!(a.handler_equal(&b).unwrap());

    let mut unknown = Param::new();
    set(&mut unknown, "zz_unknown", ParamValue::Integer(1));
    assert_eq!(
        a.set_parameters(&unknown).unwrap(),
        ["Warning: FeatureFinderAlgorithmPicked received the unknown parameter 'zz_unknown'!"]
    );
    assert!(!a.handler_equal(&b).unwrap());
    b.set_parameters(&unknown).unwrap();
    assert!(a.handler_equal(&b).unwrap());

    a.set_name("Renamed").unwrap();
    assert_eq!(a.name(), "Renamed");
    assert!(!a.handler_equal(&b).unwrap());
    assert_eq!(
        a.set_parameters(&unknown).unwrap(),
        ["Warning: Renamed received the unknown parameter 'zz_unknown'!"]
    );
    b.set_name("Renamed").unwrap();
    assert!(a.handler_equal(&b).unwrap());

    let mut refused = Param::new();
    set(&mut refused, "intensity:bins", ParamValue::Integer(0));
    assert!(a.set_parameters(&refused).is_err());
    assert!(!a.handler_equal(&b).unwrap());
    assert!(b.set_parameters(&refused).is_err());
    assert!(a.handler_equal(&b).unwrap());
}

// ---------------------------------------------------------------------------
// Integer parameters beyond the int range, a failed step 2.5 in a debug run,
// and the gnuplot formulas on non-finite operands (combined fix round 2,
// ../oracle/ffap-complete-fix2, every case twice and identical)
// ---------------------------------------------------------------------------

/// The executed rows of one case of `param_narrowing.tsv`: kind, key, value.
fn narrowing_rows() -> BTreeMap<String, Vec<(String, String, String)>> {
    let mut cases: BTreeMap<String, Vec<(String, String, String)>> = BTreeMap::new();
    for line in fixture("param_narrowing.tsv").lines().skip(2) {
        let fields: Vec<&str> = line.splitn(4, '\t').collect();
        cases.entry(fields[0].to_owned()).or_default().push((
            fields[1].to_owned(),
            fields[2].to_owned(),
            fields[3].to_owned(),
        ));
    }
    cases
}

fn rows_of<'a>(rows: &'a [(String, String, String)], kind: &str) -> Vec<(&'a str, &'a str)> {
    rows.iter()
        .filter(|(k, _, _)| k == kind)
        .map(|(_, key, value)| (key.as_str(), value.as_str()))
        .collect()
}

fn row_value<'a>(rows: &'a [(String, String, String)], kind: &str, key: &str) -> Option<&'a str> {
    rows_of(rows, kind)
        .into_iter()
        .find(|(k, _)| *k == key)
        .map(|(_, v)| v)
}

/// The typed members against the driver's dump of the protected members.
fn assert_members(
    settings: &openms::analysis::feature_finder_picked::algorithm::Settings,
    rows: &[(String, String, String)],
    tag: &str,
    case: &str,
) {
    use openms::analysis::feature_finder_picked::algorithm::ReportedMz;
    let members = rows_of(rows, tag);
    assert_eq!(members.len(), 16, "{case} {tag}");
    for (name, value) in members {
        let double = |actual: f64| {
            assert_eq!(
                format!("{:016x}", actual.to_bits()),
                value,
                "{case} {tag} {name}"
            );
        };
        match name {
            "pattern_tolerance" => double(settings.pattern_tolerance),
            "trace_tolerance" => double(settings.trace_tolerance),
            "slope_bound" => double(settings.slope_bound),
            "intensity_percentage" => double(settings.intensity_percentage),
            "intensity_percentage_optional" => double(settings.intensity_percentage_optional),
            "optional_fit_improvement" => double(settings.optional_fit_improvement),
            "mass_window_width" => double(settings.mass_window_width),
            "min_isotope_fit" => double(settings.min_isotope_fit),
            "min_trace_score" => double(settings.min_trace_score),
            "min_rt_span" => double(settings.min_rt_span),
            "max_rt_span" => double(settings.max_rt_span),
            "max_feature_intersection" => double(settings.max_feature_intersection),
            "min_spectra" => assert_eq!(settings.min_spectra.to_string(), value, "{case} {name}"),
            "max_missing" => assert_eq!(
                settings.max_missing_trace_peaks.to_string(),
                value,
                "{case} {name}"
            ),
            "intensity_bins" => {
                assert_eq!(settings.intensity_bins.to_string(), value, "{case} {name}")
            }
            "reported_mz" => {
                let expected = match value {
                    "maximum" => ReportedMz::Maximum,
                    "average" => ReportedMz::Average,
                    "monoisotopic" => ReportedMz::Monoisotopic,
                    other => panic!("{other}"),
                };
                assert_eq!(settings.reported_mz, expected, "{case} {name}");
            }
            other => panic!("unknown member {other}"),
        }
    }
}

/// A stored parameter value as the driver printed it (`ParamValue::toString`),
/// read back into its type.
fn assert_stored(parameters: &Param, key: &str, text: &str, case: &str) {
    match parameters.value(key).unwrap() {
        ParamValue::Integer(value) => assert_eq!(value.to_string(), text, "{case} {key}"),
        ParamValue::Float(value) => {
            assert_eq!(
                value.to_bits(),
                text.parse::<f64>().unwrap().to_bits(),
                "{case} {key}"
            )
        }
        other => panic!("{case} {key}: {other:?}"),
    }
}

/// The source narrows a 64-bit integer parameter to `int` for its restriction
/// check (`Param::ParamEntry::isValid`) and converts the members as the
/// Release build does: `(Int)` and `operator unsigned int` keep the low 32
/// bits, a negative value throws `ConversionError` from `operator unsigned
/// int` (half way through `updateMembers_`, or at the start of `run_` for
/// `fit:max_iterations`), and `min_spectra_` is the `cvttsd2si` of
/// `floor(value * 0.5)`. Every executed case: the outcome and its text, the
/// features bit for bit, the abort counts, the console lines, the typed
/// members afterwards, and the stored value.
#[test]
fn integer_parameters_beyond_the_int_range_follow_the_release_build() {
    let cases = narrowing_rows();
    assert_eq!(cases.len(), 21);
    let key_of = |case: &str| -> Option<(&'static str, i64)> {
        Some(match case {
            "bins_2p32_plus_10" => ("intensity:bins", 4_294_967_306),
            "bins_2p32" => ("intensity:bins", 4_294_967_296),
            "bins_2p31" => ("intensity:bins", 2_147_483_648),
            "bins_m2p32_plus_3" => ("intensity:bins", -4_294_967_293),
            "charge_high_2p32_plus_2" => ("isotopic_pattern:charge_high", 4_294_967_298),
            "charge_high_m2p32_plus_3" => ("isotopic_pattern:charge_high", -4_294_967_293),
            "charge_low_2p32_plus_2" => ("isotopic_pattern:charge_low", 4_294_967_298),
            "max_missing_2p32_plus_1" => ("mass_trace:max_missing", 4_294_967_297),
            "max_missing_2p32_minus_1" => ("mass_trace:max_missing", 4_294_967_295),
            "max_missing_m2p32_plus_1" => ("mass_trace:max_missing", -4_294_967_295),
            "max_iterations_2p32_plus_500" => ("fit:max_iterations", 4_294_967_796),
            "max_iterations_m2p32_plus_500" => ("fit:max_iterations", -4_294_966_796),
            "min_spectra_2p32_plus_10" => ("mass_trace:min_spectra", 4_294_967_306),
            "min_spectra_2p33_plus_22" => ("mass_trace:min_spectra", 8_589_934_614),
            "min_spectra_m2p32_plus_4" => ("mass_trace:min_spectra", -4_294_967_292),
            "min_spectra_2p62_plus_30" => ("mass_trace:min_spectra", 4_611_686_018_427_387_934),
            "min_spectra_2p63_minus_2p32_plus_8" => {
                ("mass_trace:min_spectra", 9_223_372_032_559_808_520)
            }
            _ => return None,
        })
    };
    let mut run_cases = 0;
    for (case, rows) in &cases {
        if case.starts_with("partial_") || case == "run_max_iterations_negative" {
            continue;
        }
        run_cases += 1;
        let mut parameters = ffc1_parameters();
        let stored_key = match key_of(case) {
            Some((key, value)) => {
                set(&mut parameters, key, ParamValue::Integer(value));
                key
            }
            None => {
                assert_eq!(case, "control");
                "intensity:bins"
            }
        };
        let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
        let mut features = FeatureMap::new();
        let result = algorithm.run(ffc1_input(), &mut features, &parameters, &FeatureMap::new());
        match row_value(rows, "outcome", "status").unwrap() {
            "returned" => result.unwrap_or_else(|error| panic!("{case}: {error}")),
            "threw" => {
                let what = row_value(rows, "outcome", "what").unwrap();
                match result {
                    Err(openms::Error::InvalidValue(message)) => {
                        assert_eq!(message, what, "{case}")
                    }
                    other => panic!("{case}: {other:?}"),
                }
            }
            other => panic!("{other}"),
        }
        let count: usize = row_value(rows, "count", "features")
            .unwrap()
            .parse()
            .unwrap();
        assert_eq!(features.len(), count, "{case}");
        let expected_features = rows_of(rows, "feature");
        assert_eq!(expected_features.len(), count, "{case}");
        for (feature, (rt, rest)) in features.features.iter().zip(expected_features) {
            let rest: Vec<&str> = rest.split(' ').collect();
            // The retention time and the intensity are fitted, and Gaussian:
            // bit for bit on every platform.
            let fitted = |actual: String, expected: &str, what: &str| {
                if actual != expected {
                    let (within, departure) =
                        close(value_of(expected), value_of(&actual), BIT_FOR_BIT);
                    assert!(
                        within,
                        "{case} {what}: {actual} against {expected} ({departure:e})"
                    );
                }
            };
            fitted(format!("{:016x}", feature.rt.to_bits()), rt, "rt");
            assert_eq!(
                format!("{:016x}", feature.mz.to_bits()),
                rest[0],
                "{case} mz"
            );
            assert_eq!(feature.charge.to_string(), rest[1], "{case} charge");
            fitted(
                format!("{:08x}", feature.intensity.to_bits()),
                rest[2],
                "intensity",
            );
        }
        let aborts: Vec<(String, String)> = algorithm
            .aborts()
            .iter()
            .map(|(reason, count)| (reason.clone(), count.to_string()))
            .collect();
        let expected_aborts: Vec<(String, String)> = rows_of(rows, "abort")
            .into_iter()
            .map(|(reason, count)| (reason.to_owned(), count.to_owned()))
            .collect();
        if result_returned(rows) {
            assert_eq!(aborts, expected_aborts, "{case}");
        }
        let console: Vec<String> = algorithm
            .report()
            .iter()
            .filter_map(|line| line.text().map(str::to_owned))
            .collect();
        let expected_console: Vec<String> = rows_of(rows, "stdout")
            .into_iter()
            .map(|(_, line)| line.to_owned())
            .collect();
        assert_eq!(console, expected_console, "{case}");
        assert_members(algorithm.settings(), rows, "members_after", case);
        assert_stored(
            algorithm.parameters(),
            stored_key,
            row_value(rows, "stored", stored_key).unwrap(),
            case,
        );
    }
    assert_eq!(run_cases, 18);

    // A refused set after an accepted one: updateMembers_ threw half way.
    for case in ["partial_max_missing", "partial_bins"] {
        let rows = &cases[case];
        let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
        let mut first = ffc1_parameters();
        set(
            &mut first,
            "mass_trace:mz_tolerance",
            ParamValue::Float(0.05),
        );
        set(&mut first, "mass_trace:min_spectra", ParamValue::Integer(4));
        set(&mut first, "mass_trace:max_missing", ParamValue::Integer(3));
        set(
            &mut first,
            "isotopic_pattern:mass_window_width",
            ParamValue::Float(20.0),
        );
        set(&mut first, "intensity:bins", ParamValue::Integer(7));
        set(
            &mut first,
            "feature:min_isotope_fit",
            ParamValue::Float(0.7),
        );
        algorithm.set_parameters(&first).unwrap();
        assert_members(algorithm.settings(), rows, "members_first", case);
        let mut second = ffc1_parameters();
        set(
            &mut second,
            "mass_trace:mz_tolerance",
            ParamValue::Float(0.02),
        );
        set(
            &mut second,
            "mass_trace:min_spectra",
            ParamValue::Integer(8),
        );
        set(
            &mut second,
            "isotopic_pattern:mass_window_width",
            ParamValue::Float(30.0),
        );
        set(
            &mut second,
            "feature:min_isotope_fit",
            ParamValue::Float(0.9),
        );
        if case == "partial_max_missing" {
            set(
                &mut second,
                "mass_trace:max_missing",
                ParamValue::Integer(-4_294_967_295),
            );
            set(&mut second, "intensity:bins", ParamValue::Integer(9));
        } else {
            set(
                &mut second,
                "mass_trace:max_missing",
                ParamValue::Integer(2),
            );
            set(
                &mut second,
                "intensity:bins",
                ParamValue::Integer(-4_294_967_293),
            );
        }
        assert_eq!(
            row_value(rows, "outcome", "exception"),
            Some("ConversionError")
        );
        match algorithm.set_parameters(&second) {
            Err(openms::Error::InvalidValue(message)) => {
                assert_eq!(
                    message,
                    row_value(rows, "outcome", "what").unwrap(),
                    "{case}"
                )
            }
            other => panic!("{case}: {other:?}"),
        }
        assert_members(algorithm.settings(), rows, "members_second", case);
        for (key, text) in rows_of(rows, "stored") {
            assert_stored(algorithm.parameters(), key, text, case);
        }
    }

    // `fit:max_iterations` is read at the start of run_, before the user
    // seeds are sorted and before anything touches the caller's map.
    let rows = &cases["run_max_iterations_negative"];
    let mut parameters = ffc1_parameters();
    set(
        &mut parameters,
        "fit:max_iterations",
        ParamValue::Integer(-4_294_966_796),
    );
    let mut features = FeatureMap::new();
    let mut kept = Feature::new(1.0, 2.0, 0.0);
    kept.intensity = 0.0;
    features.features.push(kept);
    let mut seeds = FeatureMap::new();
    for mz in [700.0, 500.0, 600.0] {
        seeds.features.push(Feature::new(4200.0, mz, 0.0));
    }
    let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
    match algorithm.run(ffc1_input(), &mut features, &parameters, &seeds) {
        Err(openms::Error::InvalidValue(message)) => {
            assert_eq!(message, row_value(rows, "outcome", "what").unwrap())
        }
        other => panic!("{other:?}"),
    }
    assert_eq!(
        features.len().to_string(),
        row_value(rows, "count", "features").unwrap()
    );
    assert_eq!(row_value(rows, "count", "windows"), Some("0"));
    assert!(algorithm.isotope_windows().is_none());
    let order: Vec<String> = algorithm
        .seeds()
        .features
        .iter()
        .map(|seed| seed.mz.to_string())
        .collect();
    assert_eq!(order.join(" "), row_value(rows, "seeds", "mz").unwrap());
    assert_members(
        algorithm.settings(),
        rows,
        "members_after",
        "run_max_iterations_negative",
    );
}

fn result_returned(rows: &[(String, String, String)]) -> bool {
    row_value(rows, "outcome", "status") == Some("returned")
}

/// The executed rows of `length_error.tsv` for one input and mode.
fn length_error_rows(input: &str, mode: &str, kind: &str) -> Vec<String> {
    fixture("length_error.tsv")
        .lines()
        .skip(2)
        .filter_map(|line| {
            let fields: Vec<&str> = line.splitn(4, '\t').collect();
            (fields[0] == input && fields[1] == mode && fields[2] == kind)
                .then(|| fields[3].to_owned())
        })
        .collect()
}

/// FFC_1 with the last m/z of the last spectrum set to `mz`: the oracle's
/// `huge_mz_1e19.mzML` and `huge_mz_2e18.mzML`
/// (`../oracle/ffap-complete-fix2/make_inputs.py`), which the tool test
/// derives byte for byte.
fn huge_mz_input(mz: f64) -> MSExperiment {
    let mut experiment = ffc1_input();
    let last = experiment.spectra.last_mut().unwrap();
    let peak = last.peaks.last_mut().unwrap();
    assert!(peak.mz < mz);
    peak.mz = mz;
    experiment
}

/// A debug run whose step 2.5 fails: the Release build has opened
/// `debug/log.txt` and created `debug/features` before step 1 and written the
/// step-1 line, which its destructor flushes (the executed log is empty while
/// the object lives and 40 bytes after it). Above `vector::max_size()` the
/// source throws `std::length_error` (`what()` `vector::_M_default_append`);
/// below it the allocation of `4e16` windows fails with `std::bad_alloc`,
/// where the port's native window ceiling refuses instead, with the same debug
/// output. The object keeps its (empty) isotope windows, and its stream stays
/// open: a second debug run on FFC_1 writes the seed map, the abort map and the
/// input exactly as a fresh object does (the executed files are identical),
/// but no log.
#[test]
fn a_debug_run_that_fails_in_step_two_point_five_keeps_the_executed_debug_output() {
    use openms::analysis::feature_finder_picked::seeds::{LENGTH_ERROR_WHAT, is_length_error};
    let log = std::fs::read(fixture_path("length_error_log.txt")).unwrap();
    assert_eq!(log, b"Precalculating intensity thresholds ...\n");
    let debug_parameters = with_debug(
        ffc1_parameters(),
        &[("feature:min_isotope_fit", ParamValue::Float(1.0))],
    );
    for (tag, mz) in [("1e19", 1e19), ("2e18", 2e18)] {
        let report = length_error_rows(tag, "lenerr_single", "report");
        let mut algorithm = FeatureFinderAlgorithmPicked::new().unwrap();
        let mut features = FeatureMap::new();
        let error = algorithm
            .run(
                huge_mz_input(mz),
                &mut features,
                &debug_parameters,
                &FeatureMap::new(),
            )
            .unwrap_err();
        if tag == "1e19" {
            assert_eq!(
                report[0],
                format!("run1 threw St12length_error what={LENGTH_ERROR_WHAT}")
            );
            assert!(is_length_error(&error), "{error}");
        } else {
            assert_eq!(report[0], "run1 threw St9bad_alloc what=std::bad_alloc");
            assert!(!is_length_error(&error), "{error}");
            assert!(error.to_string().contains("exceed the limit"), "{error}");
        }
        assert_eq!(
            report[1],
            "run1 features=0 aborts=0 abort_reasons=0 windows=0"
        );
        assert!(features.is_empty());
        assert!(algorithm.aborts().is_empty());
        assert!(algorithm.abort_reasons().is_empty());
        assert!(algorithm.isotope_windows().is_none());
        let out = algorithm.debug_output().unwrap();
        assert!(out.log_opened);
        assert!(out.termination.is_none());
        assert_eq!(out.log.text().as_bytes(), log.as_slice());
        assert_eq!(
            length_error_rows(tag, "lenerr_single", "log_bytes_destroyed"),
            [log.len().to_string()]
        );
        assert_eq!(
            length_error_rows(tag, "lenerr_single", "tree"),
            ["debug", "debug/features", "debug/log.txt"]
        );
        assert!(out.seed_maps.is_empty() && out.feature_files.is_empty());
        assert!(out.abort_reasons.is_none() && out.input.is_none());

        // The reused object: its stream is open, so the second run's log is
        // dropped, and its files are a fresh object's.
        let mut second = FeatureMap::new();
        algorithm
            .run(
                ffc1_input(),
                &mut second,
                &debug_parameters,
                &FeatureMap::new(),
            )
            .unwrap();
        let (fresh_result, fresh, fresh_features) =
            debug_run(ffc1_input(), &debug_parameters, Options::default());
        fresh_result.unwrap();
        let reuse_report = length_error_rows(tag, "lenerr_reuse", "report");
        assert_eq!(
            reuse_report[3],
            "reuse_run2 features=0 aborts=1 abort_reasons=25 windows=15"
        );
        assert_eq!(
            reuse_report[6],
            "fresh_run features=0 aborts=1 abort_reasons=25 windows=15"
        );
        for (object, map) in [(&algorithm, &second), (&fresh, &fresh_features)] {
            assert!(map.is_empty());
            assert_eq!(object.abort_reasons().len(), 25);
            assert_eq!(object.isotope_windows().unwrap().patterns().len(), 15);
            assert_eq!(
                object.aborts().iter().collect::<Vec<_>>(),
                [(
                    &"Could not find good enough isotope pattern containing the seed".to_owned(),
                    &25
                )]
            );
        }
        assert_eq!(
            length_error_rows(tag, "lenerr_reuse", "compare"),
            [
                "seeds_2.featureXML same",
                "abort_reasons.featureXML same",
                "input.mzML same"
            ]
        );
        let reused = algorithm.debug_output().unwrap();
        let fresh_out = fresh.debug_output().unwrap();
        assert!(!reused.log_opened);
        assert!(reused.log.text().is_empty());
        assert_eq!(
            length_error_rows(tag, "lenerr_reuse", "log_bytes_reuse"),
            [log.len().to_string()]
        );
        assert!(fresh_out.log_opened);
        assert_eq!(
            length_error_rows(tag, "lenerr_reuse", "log_bytes_fresh"),
            [fresh_out.log.text().len().to_string()]
        );
        assert_eq!(reused.seed_maps, fresh_out.seed_maps);
        assert_eq!(reused.abort_reasons, fresh_out.abort_reasons);
        assert_eq!(reused.input, fresh_out.input);
        assert!(reused.feature_files.is_empty() && fresh_out.feature_files.is_empty());
        assert_eq!(
            length_error_rows(tag, "lenerr_reuse", "tree_reuse"),
            [
                "debug",
                "debug/abort_reasons.featureXML",
                "debug/features",
                "debug/input.mzML",
                "debug/log.txt",
                "debug/seeds_2.featureXML"
            ]
        );
        // The fresh object's files are the executed a2 files.
        assert_executed_bytes("a2", "log.txt", fresh_out.log.text().as_bytes());
    }
}

/// `getGnuplotFormula` of both fitters on non-finite operands, against the
/// executed Release build (`gnuplot_formula_nonfinite.tsv`, 648 rows): the
/// sum `rt_shift + centre` keeps `rt_shift`'s NaN when both are NaN and gives
/// SSE's default NaN (`-nan`) for `inf + -inf`, the product
/// `theoretical_int * height` likewise, and EGH's `2 * sigma * sigma` is
/// `(sigma + sigma) * sigma`; the Gaussian stores `|sigma|`, a NaN's sign
/// included.
#[test]
fn gnuplot_formulas_print_the_executed_nan_signs() {
    use openms::analysis::feature_finder_picked::egh_trace_fitter::EGHTraceFitter;
    use openms::analysis::feature_finder_picked::gauss_trace_fitter::GaussTraceFitter;
    use openms::analysis::feature_finder_picked::helper_structs::MassTrace;
    use openms::analysis::feature_finder_picked::trace_fitter::TraceFitter;
    let table = fixture("gnuplot_formula_nonfinite.tsv");
    let bits = |text: &str| f64::from_bits(u64::from_str_radix(text, 16).unwrap());
    let mut rows = 0;
    let mut negative_nan = 0;
    for line in table.lines().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 10, "{line}");
        let (theoretical_int, height, rt_shift, centre, sigma, tau, baseline) = (
            bits(fields[2]),
            bits(fields[3]),
            bits(fields[4]),
            bits(fields[5]),
            bits(fields[6]),
            bits(fields[7]),
            bits(fields[8]),
        );
        let trace = MassTrace {
            theoretical_int,
            ..MassTrace::default()
        };
        let formula = match fields[0] {
            "gauss" => {
                let mut fitter = GaussTraceFitter::new();
                fitter.set_optimized_parameters([height, centre, sigma]);
                fitter.gnuplot_formula(&trace, 'f', baseline, rt_shift)
            }
            "egh" => {
                let mut fitter = EGHTraceFitter::new();
                fitter.set_optimized_parameters([height, centre, sigma, tau]);
                fitter.gnuplot_formula(&trace, 'f', baseline, rt_shift)
            }
            other => panic!("{other}"),
        };
        assert_eq!(formula, fields[9], "{line}");
        rows += 1;
        negative_nan += usize::from(fields[9].contains("-nan"));
    }
    assert_eq!(rows, 648);
    assert!(negative_nan > 100);
}
