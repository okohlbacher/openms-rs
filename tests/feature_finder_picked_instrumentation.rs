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
//! `EGH_*` parameters) differ from the Release build's Eigen in their last
//! bits, as `docs/TRACE_FITTER_SUPPORT.md` records; they are compared within
//! `FIT_RELATIVE` and everything else exactly.

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
    Options, PseudoRtShiftKey, run_with_options,
};
use openms::analysis::feature_finder_picked::debug::{DebugOutput, ReportLine};
use openms::analysis::feature_finder_picked::instance::FeatureFinderAlgorithmPicked;
use openms::analysis::feature_finder_picked::source_sort::source_sort_permutation;
use openms::concept::parallel::Threads;
use openms::concept::progress_logger::{CommandProgressLogger, ProgressLogType, ProgressLogger};
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

/// The relative tolerance of fitted values: the Levenberg-Marquardt
/// transcription departs from the Release build's Eigen in the last bits.
const FIT_RELATIVE: f64 = 1e-9;

/// Whether a dump token holds a fitted value, compared within
/// [`FIT_RELATIVE`]; everything else is compared exactly.
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

fn close(expected: f64, actual: f64) -> bool {
    expected == actual
        || (expected.is_nan() && actual.is_nan())
        || (expected - actual).abs() <= FIT_RELATIVE * expected.abs().max(actual.abs())
}

/// Compare two dumps: identical line structure, exact tokens, and fitted
/// values within [`FIT_RELATIVE`]. Returns the number of fitted values that
/// were not bit-identical.
fn assert_dumps_match(expected: &str, actual: &str, case: &str) -> usize {
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
            assert!(
                fitted && close(value_of(x), value_of(y)),
                "{case}: line {number} token {index}: expected {x}, got {y}\n{e}\n{a}"
            );
            inexact += 1;
        }
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
    // violation, and leaves the parameters as they were, where the source has
    // already assigned the rejected set (native difference, recorded in the
    // support document): the executed `getParameters()` then shows `bins 0`.
    let mut bad = Param::new();
    set(&mut bad, "intensity:bins", ParamValue::Integer(0));
    let error = algorithm.set_parameters(&bad).unwrap_err().to_string();
    assert!(error.contains("bins"), "{error}");
    assert!(info.contains(
        "bad_bins InvalidParameter FeatureFinderAlgorithmPicked: Invalid integer parameter value '0' for parameter 'bins' given!"
    ));
    assert_eq!(dump_param(algorithm.parameters()), with_unknown);
    let failed = fixture("params_after_failed_set.txt");
    assert!(failed.contains("P intensity:bins i 0\n"));
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
    eprintln!("reuse: {inexact} fitted values within {FIT_RELATIVE} but not bit-identical");
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
    eprintln!("prefilled: {inexact} fitted values within {FIT_RELATIVE} but not bit-identical");
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
    fixture("debug_digests.tsv")
        .lines()
        .skip(1)
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
    let powf_steps = POWF_STEPS.get_or_init(|| Mutex::new(0));
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
                // A NaN equals a NaN (decision D6). Its sign is the
                // hardware's default NaN of `0.0 / 0.0`, negative on x86_64
                // as in the executed build and positive on arm64.
                if u.to_bits() == v.to_bits() || (u.is_nan() && v.is_nan()) {
                    continue;
                }
                // The overall score is the C++ `std::pow(float, float)`,
                // glibc's `powf` in the Release build, where the port computes
                // the correctly rounded power (`seeds::overall_score`,
                // CPP-272); the two differ by one binary32 step on a few
                // scores. Every other array is bit-identical.
                let step = (i64::from(u.to_bits()) - i64::from(v.to_bits())).abs();
                assert!(
                    x.name.starts_with("overall_score_") && step == 1,
                    "{case} spectrum {index} array {}: {u} against {v}",
                    x.name
                );
                *powf_steps.lock().unwrap() += 1;
            }
        }
    }
}

/// How many overall scores differed from the Release build by one binary32
/// step in the input comparisons of this process.
static POWF_STEPS: std::sync::OnceLock<Mutex<usize>> = std::sync::OnceLock::new();

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
        assert_eq!(termination.charge, charge, "{case}");
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
    // A string under the key has no reproducible value.
    let mut parameters = with_debug(ffc1_parameters(), &[]);
    set(
        &mut parameters,
        "debug:pseudo_rt_shift",
        ParamValue::String("500".into()),
    );
    let (result, algorithm, _) = debug_run(ffc1_input(), &parameters, Options::default());
    assert!(matches!(result, Err(openms::Error::Unsupported(_))));
    assert!(algorithm.debug_output().unwrap().termination.is_none());
}

/// A driver case `declared-*`: its name, its parameter overrides and the
/// fixture of its executed feature map.
type DeclaredCase = (&'static str, Vec<(&'static str, ParamValue)>, &'static str);

/// The executed driver cases `declared-*`: with `debug:pseudo_rt_shift` in the
/// parameters the source writes its feature files, and so does the port under
/// either key policy. Every `.dta`, `_cropped.dta` and `.plot` file and the
/// log are byte-identical to the executed ones; the feature map is the plain
/// run's.
#[test]
fn feature_debug_files_match_the_release_build() {
    let digests = debug_digests();
    let cases: [DeclaredCase; 4] = [
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
    ];
    let mut inexact = 0;
    for (case, overrides, expected_map) in cases {
        for policy in [PseudoRtShiftKey::Source, PseudoRtShiftKey::Declared] {
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
    eprintln!("declared: {inexact} fitted values within {FIT_RELATIVE} but not bit-identical");
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
/// rather than called, so an inverted range never reaches the progress
/// logger's own range check.
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
    eprintln!("overlap: {inexact} fitted values within {FIT_RELATIVE} but not bit-identical");
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
