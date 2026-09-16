// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The feature stage of `FeatureFinderAlgorithmPicked`: isotope fit, mass-trace
//! extension, trace fitting, cropping, quality checks, the parallel seed loop
//! and overlap resolution.
//!
//! Evidence (see `docs/FEATURE_FINDER_PICKED_SUPPORT.md` and
//! `tests/data/feature_finder_picked_provenance.json`):
//!
//! - tier 1, executed C++ (the Linux x86_64 Release build
//!   `openms4-release-bc9cc12-c19e494-174b576`, the reference platform): the
//!   final `FeatureMap` of `FeatureFinderAlgorithmPicked::run` for six
//!   configurations (FFC_1 symmetric, FFC_1 asymmetric, FFC_1 with user seeds,
//!   the class-test input and the two `#9247` tolerance swaps), its printed
//!   seed and candidate counts and its `aborts_` map, in
//!   `b7_feature_records.tsv`. The fixtures were re-extracted from that build
//!   by `../oracle/ffap-sem-completion/extract/extract_linux.py`, which runs the
//!   B7 extraction unchanged; they replaced the macOS arm64 product-SDK (Debug)
//!   capture of package B7;
//! - adapted: the per-seed intermediate state of `b7_seed_records.tsv`, which
//!   the C2 driver produced by replaying the protected library steps
//!   (`findBestIsotopeFit_`, `extendMassTraces_`, the chosen fitter,
//!   `cropFeature_`, `checkFeatureQuality_`) outside the source's OpenMP
//!   region; the driver checks its own step-4 replica against the library
//!   output;
//! - tier 3: the literals of `FeatureFinderAlgorithmPicked_test.cpp`;
//! - tier 4: hand-derived cases for `intersection_`, the preserved
//!   `extendMassTraces_` defect and the resource ceilings.
//!
//! Peak identities, counts, charges, labels and abort reasons are compared
//! exactly. Coordinates, intensities, qualities and fitted parameters are
//! compared bit for bit on Linux x86_64 with glibc, except for the asymmetric
//! (EGH) configuration, and within measured platform bounds elsewhere; see
//! [`tolerance`].

#![cfg(all(feature = "mzml", feature = "paramxml"))]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use openms::analysis::feature_finder_picked::algorithm::{
    AbundanceOverride, Limits, Options, RtShape, default_parameters, run, run_with_options,
};
use openms::analysis::feature_finder_picked::extension::{
    OverallScores, extend_mass_traces, find_best_isotope_fit,
};
use openms::analysis::feature_finder_picked::fitting::{
    FittedModel, QualityOutcome, check_feature_quality, crop_feature,
};
use openms::analysis::feature_finder_picked::helper_structs::{
    MassTrace, MassTraces, PatternPeak, TracePeak,
};
use openms::analysis::feature_finder_picked::resolution::intersection;
use openms::analysis::feature_finder_picked::seeds::SeedStage;
use openms::analysis::feature_finder_picked::trace_fitter::TraceFitterParams;
use openms::concept::parallel::Threads;
use openms::format::{FileHandler, FileType, PeakFileOptions, featurexml, paramxml};
use openms::kernel::{ConvexHull2D, Feature, FeatureMap, MSExperiment, NumericRange, Point2D};
use openms::metadata::{MetaValue, MetaValueData};
use openms::param::{Param, ParamValue};

// ---------------------------------------------------------------------------
// Inputs, parameters and fixtures
// ---------------------------------------------------------------------------

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/feature_finder_picked")
        .join(name)
}

/// The FeatureFinderCentroided_1 input, reused read-only from the A3 fixtures
/// (sha256 a3dfae63..., test-data 0cb15f2).
fn ffc1_input() -> MSExperiment {
    let mut options = PeakFileOptions::default();
    options.add_ms_level(1).unwrap();
    options.set_intensity_range(NumericRange {
        min: 0.0,
        max: f64::MAX,
    });
    FileHandler::load_experiment_with_options(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/mzml_mobility/FeatureFinderCentroided_1_input.mzML"),
        &[FileType::MzMl],
        &options,
    )
    .unwrap()
}

/// `FeatureFinderAlgorithmPicked_test.cpp` loading: MS1 only.
fn class_test_input() -> MSExperiment {
    let mut options = PeakFileOptions::default();
    options.add_ms_level(1).unwrap();
    FileHandler::load_experiment_with_options(
        data("FeatureFinderAlgorithmPicked.mzML"),
        &[FileType::MzMl],
        &options,
    )
    .unwrap()
}

fn ffc1_parameters() -> Param {
    paramxml::load(data("FeatureFinderCentroided_1_parameters.ini"))
        .unwrap()
        .copy("FeatureFinderCentroided:1:algorithm:", true)
        .unwrap()
}

fn class_test_parameters() -> Param {
    paramxml::load(data("FeatureFinderAlgorithmPicked.ini"))
        .unwrap()
        .copy("FeatureFinder:1:algorithm:", true)
        .unwrap()
}

/// The eight features of the retained FeatureFinderCentroided_1 output, used as
/// user seeds by the `ffc1_user_seeds` configuration.
fn ffc1_user_seeds() -> FeatureMap {
    featurexml::load(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/mzml_mobility/FeatureFinderCentroided_1_1_output.featureXML"),
    )
    .unwrap()
}

fn set(param: &mut Param, key: &str, value: ParamValue) {
    param.set_value(key, value, "", &[]).unwrap();
}

fn f64_hex(text: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(text, 16).unwrap())
}

fn f32_hex(text: &str) -> f32 {
    f32::from_bits(u32::from_str_radix(text, 16).unwrap())
}

fn rows(file: &str) -> Vec<Vec<String>> {
    std::fs::read_to_string(data(file))
        .unwrap()
        .lines()
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect()
}

/// Rows of one kind and configuration, without those two columns.
fn records(file: &str, kind: &str, config: &str) -> Vec<Vec<String>> {
    rows(file)
        .into_iter()
        .filter(|row| row[0] == kind && row[1] == config)
        .map(|row| row[2..].to_vec())
        .collect()
}

/// Bit equality, for the quantities measured to agree exactly with the executed
/// C++ on every platform: every isotope-fit score, every isotope-pattern
/// intensity and m/z score, and every mass trace (peak identity, theoretical
/// intensity and baseline).
const BITWISE: f64 = 0.0;

/// The comparison bound of the fitted parameters, the qualities and the
/// feature coordinates of one seed or feature (`index`, `None` for a final
/// feature) of one configuration.
///
/// The fixtures are the **Linux x86_64 Release** build
/// (`openms4-release-bc9cc12-c19e494-174b576`, the C2 driver `ffap_stages` run on
/// ibminode06, AMD EPYC 7763, glibc 2.39; `../oracle/ffap-sem-completion`), the
/// reference platform the user chose on 2026-09-15. Since lane B3b the port's
/// Levenberg-Marquardt solver follows that build's Eigen kernels, and its
/// Gaussian fit calls the platform `exp` and `log`, as the source does, so the
/// bound depends on the platform the test runs on:
///
/// - Linux x86_64 with glibc: every Gaussian configuration is **bit for bit**,
///   measured on dax (AMD EPYC 9654). glibc selects FMA variants of `exp` and
///   `log` on CPUs that have FMA, as both measured hosts do; the exact
///   comparison assumes such a CPU. The asymmetric (EGH) configuration departs
///   by at most `2.3038e-12` relative (seed 24's lower retention-time bound):
///   `EGHTraceFitter` in this port calls the `libm` crate's `exp`, `log` and
///   `atan` where the source calls glibc's (`docs/EGH_TRACE_FITTER_SUPPORT.md`),
///   so [`EGH_LIBM_GAP`] bounds it.
/// - macOS arm64 (Apple libm), measured against the same Linux capture: the
///   Gaussian fits depart by at most `5.355e-13` relative, except seeds 11 and 12
///   of `classtest_9247_tight_pattern` ([`KNOWN_FIT_GAP_MACOS`]); the EGH
///   configuration departs by the same `2.3038e-12` as on Linux.
/// - Any other platform is unmeasured: the bound is the work package's `1e-9`
///   contract, and the two ill-conditioned seeds keep the largest departure
///   recorded before the Linux capture existed ([`KNOWN_FIT_GAP_UNMEASURED`]).
///
/// The two ill-conditioned seeds are rejected by `checkFeatureQuality_` in the
/// executed C++ and here, with the same reason, on every measured platform, so
/// no feature changes (`every_seed_matches_the_executed_intermediate_state`
/// asserts that a seed with a platform gap never becomes a feature).
fn tolerance(config: &str, index: Option<usize>) -> f64 {
    if config == "ffc1_asymmetric" {
        return EGH_LIBM_GAP;
    }
    platform_tolerance(config, index)
}

/// The EGH configuration's measured departure, `2.3038102266706174e-12`
/// relative on both Linux x86_64 and macOS arm64, rounded up at the second
/// significant digit.
const EGH_LIBM_GAP: f64 = 2.4e-12;

#[cfg(all(target_os = "linux", target_arch = "x86_64", target_env = "gnu"))]
fn platform_tolerance(_config: &str, _index: Option<usize>) -> f64 {
    0.0
}

/// macOS arm64: the two seeds whose Gaussian fit is ill-conditioned enough for
/// Apple's `exp` to move the fitted area by `1.0698e-3` relative (seed 11; seed
/// 12 fits the same traces), and `5.355e-13` for every other fitted quantity.
#[cfg(all(target_os = "macos", target_arch = "aarch64"))]
const KNOWN_FIT_GAP_MACOS: [(&str, usize, f64); 2] = [
    ("classtest_9247_tight_pattern", 11, 1.1e-3),
    ("classtest_9247_tight_pattern", 12, 1.1e-3),
];

#[cfg(all(target_os = "macos", target_arch = "aarch64"))]
fn platform_tolerance(config: &str, index: Option<usize>) -> f64 {
    KNOWN_FIT_GAP_MACOS
        .iter()
        .find(|(c, i, _)| *c == config && Some(*i) == index)
        .map_or(5.4e-13, |(_, _, bound)| *bound)
}

/// Unmeasured platforms: the work package's `1e-9` contract, and for the two
/// ill-conditioned seeds the largest departure measured before the Linux
/// capture, `2.25e-3` (Linux x86_64 against the macOS arm64 product SDK).
#[cfg(not(any(
    all(target_os = "linux", target_arch = "x86_64", target_env = "gnu"),
    all(target_os = "macos", target_arch = "aarch64")
)))]
const KNOWN_FIT_GAP_UNMEASURED: [(&str, usize, f64); 2] = [
    ("classtest_9247_tight_pattern", 11, 2.3e-3),
    ("classtest_9247_tight_pattern", 12, 2.3e-3),
];

#[cfg(not(any(
    all(target_os = "linux", target_arch = "x86_64", target_env = "gnu"),
    all(target_os = "macos", target_arch = "aarch64")
)))]
fn platform_tolerance(config: &str, index: Option<usize>) -> f64 {
    KNOWN_FIT_GAP_UNMEASURED
        .iter()
        .find(|(c, i, _)| *c == config && Some(*i) == index)
        .map_or(1e-9, |(_, _, bound)| *bound)
}

/// Whether a seed's fit has a platform gap larger than the configuration's
/// general bound on the platform running the test.
fn has_platform_gap(config: &str, index: usize) -> bool {
    tolerance(config, Some(index)) > tolerance(config, None)
}

#[track_caller]
fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    if actual.to_bits() == expected.to_bits() {
        return;
    }
    let deviation = if expected == 0.0 {
        actual.abs()
    } else {
        ((actual - expected) / expected).abs()
    };
    assert!(
        deviation <= tolerance,
        "{what}: {actual:e} differs from the executed {expected:e} by {deviation:e} relative, \
         above {tolerance:e}"
    );
}

/// One executed configuration and how to reproduce it.
struct Case {
    config: &'static str,
    experiment: fn() -> MSExperiment,
    parameters: fn() -> Param,
    seeds: fn() -> FeatureMap,
}

fn no_seeds() -> FeatureMap {
    FeatureMap::new()
}

fn tight_pattern_parameters() -> Param {
    let mut p = class_test_parameters();
    set(
        &mut p,
        "isotopic_pattern:mz_tolerance",
        ParamValue::Float(0.005),
    );
    set(&mut p, "mass_trace:mz_tolerance", ParamValue::Float(0.5));
    p
}

fn tight_trace_parameters() -> Param {
    let mut p = class_test_parameters();
    set(
        &mut p,
        "isotopic_pattern:mz_tolerance",
        ParamValue::Float(0.5),
    );
    set(&mut p, "mass_trace:mz_tolerance", ParamValue::Float(0.005));
    p
}

fn asymmetric_parameters() -> Param {
    let mut p = ffc1_parameters();
    set(
        &mut p,
        "feature:rt_shape",
        ParamValue::String("asymmetric".into()),
    );
    p
}

const CASES: [Case; 6] = [
    Case {
        config: "ffc1_symmetric",
        experiment: ffc1_input,
        parameters: ffc1_parameters,
        seeds: no_seeds,
    },
    Case {
        config: "ffc1_asymmetric",
        experiment: ffc1_input,
        parameters: asymmetric_parameters,
        seeds: no_seeds,
    },
    Case {
        config: "ffc1_user_seeds",
        experiment: ffc1_input,
        parameters: ffc1_parameters,
        seeds: ffc1_user_seeds,
    },
    Case {
        config: "classtest",
        experiment: class_test_input,
        parameters: class_test_parameters,
        seeds: no_seeds,
    },
    Case {
        config: "classtest_9247_tight_pattern",
        experiment: class_test_input,
        parameters: tight_pattern_parameters,
        seeds: no_seeds,
    },
    Case {
        config: "classtest_9247_tight_trace",
        experiment: class_test_input,
        parameters: tight_trace_parameters,
        seeds: no_seeds,
    },
];

// ---------------------------------------------------------------------------
// The class test (tier 3)
// ---------------------------------------------------------------------------

/// `FeatureFinderAlgorithmPicked_test.cpp`, `run()` section.
#[test]
fn class_test_run_finds_the_expected_eight_features() {
    let output = run(
        class_test_input(),
        &FeatureMap::new(),
        &class_test_parameters(),
    )
    .unwrap();
    let features = &output.features.features;
    assert_eq!(features.len(), 8);
    for (index, expected) in [(0usize, 88i64), (3, 71), (7, 47)] {
        assert_eq!(
            features[index].metadata["num_of_datapoints"].data(),
            &MetaValueData::Integer(expected)
        );
    }
    let qualities = [
        0.8826, 0.8680, 0.9077, 0.9270, 0.9398, 0.9098, 0.9403, 0.9245,
    ];
    for (feature, expected) in features.iter().zip(qualities) {
        assert!(
            (f64::from(feature.quality) - expected).abs() <= 0.001,
            "quality {} vs {expected}",
            feature.quality
        );
    }
    let intensities = [
        51366.2, 44767.6, 34731.1, 19494.2, 12570.2, 8532.26, 7318.62, 5038.81,
    ];
    for (feature, expected) in features.iter().zip(intensities) {
        assert!(
            (f64::from(feature.intensity) - expected).abs() <= 20.0,
            "intensity {} vs {expected}",
            feature.intensity
        );
    }
}

/// `FeatureFinderAlgorithmPicked_test.cpp`, the `[EXTRA]` #9247 section: the
/// isotope-pattern and mass-trace tolerances are not interchangeable.
#[test]
fn class_test_tolerance_swap_is_directional() {
    let first = run(
        class_test_input(),
        &FeatureMap::new(),
        &tight_pattern_parameters(),
    )
    .unwrap();
    let second = run(
        class_test_input(),
        &FeatureMap::new(),
        &tight_trace_parameters(),
    )
    .unwrap();
    assert_eq!(first.features.len(), 1);
    assert_eq!(second.features.len(), 0);
    let feature = &first.features.features[0];
    assert_eq!(
        feature.metadata["num_of_datapoints"].data(),
        &MetaValueData::Integer(33)
    );
    assert!((feature.rt - 4278.1601).abs() <= 0.001, "{}", feature.rt);
    assert!((feature.mz - 653.7722).abs() <= 0.001, "{}", feature.mz);
    assert!(
        (f64::from(feature.quality) - 0.9609).abs() <= 0.001,
        "{}",
        feature.quality
    );
    assert!(
        (f64::from(feature.intensity) - 18467.8).abs() <= 20.0,
        "{}",
        feature.intensity
    );
}

// ---------------------------------------------------------------------------
// The executed library output (tier 1)
// ---------------------------------------------------------------------------

/// Every executed configuration: counts, printed lines, abort reasons and the
/// features themselves.
#[test]
fn every_configuration_matches_the_executed_library() {
    for case in &CASES {
        let output = run_with_options(
            (case.experiment)(),
            &(case.seeds)(),
            &(case.parameters)(),
            &Options {
                threads: Threads::serial(),
                ..Options::default()
            },
        )
        .unwrap_or_else(|error| panic!("{}: {error}", case.config));
        check_run(case.config, &output.features, &output.log, &output.aborts);
    }
}

fn check_run(config: &str, map: &FeatureMap, log: &[String], aborts: &BTreeMap<String, usize>) {
    let run_row = &records("b7_feature_records.tsv", "run", config)[0];
    assert_eq!(
        map.len(),
        run_row[3].parse::<usize>().unwrap(),
        "{config}: feature count"
    );
    // The two `std::cout` lines of the source, adjacent and in its order; the
    // class-test INI's unknown `debug` entry puts a parameter warning first.
    let seed_line = log
        .iter()
        .position(|line| line == &run_row[4])
        .unwrap_or_else(|| panic!("{config}: no seed line in {log:?}"));
    assert_eq!(&log[seed_line + 1], &run_row[5], "{config}: candidate line");
    assert!(
        log.contains(&format!(
            "Removed {} overlapping features.",
            run_row[2].parse::<usize>().unwrap()
        )),
        "{config}: overlap count, log {log:?}"
    );
    assert!(
        log.contains(&format!("{} features found.", map.len())),
        "{config}: feature count line"
    );

    let expected_aborts: BTreeMap<String, usize> =
        records("b7_feature_records.tsv", "abort", config)
            .into_iter()
            .map(|row| (row[1].clone(), row[0].parse().unwrap()))
            .collect();
    assert_eq!(aborts, &expected_aborts, "{config}: abort reasons");

    check_features("b7_feature_records.tsv", config, map);
}

/// The recorded `feature`, `meta` and `hull` rows of `config` in `file`
/// against `map`.
fn check_features(file: &str, config: &str, map: &FeatureMap) {
    let expected: Vec<Vec<String>> = records(file, "feature", config);
    assert_eq!(
        expected.len(),
        map.len(),
        "{config}: recorded feature count"
    );
    for (index, row) in expected.iter().enumerate() {
        let feature = &map.features[index];
        let what = |field: &str| format!("{config}[{index}].{field}");
        assert_eq!(row[0].parse::<usize>().unwrap(), index);
        let relative = tolerance(config, None);
        close(feature.rt, f64_hex(&row[1]), relative, &what("rt"));
        close(feature.mz, f64_hex(&row[2]), relative, &what("mz"));
        close(
            f64::from(feature.intensity),
            f64::from(f32_hex(&row[3])),
            relative,
            &what("intensity"),
        );
        assert_eq!(
            feature.charge,
            row[4].parse::<i32>().unwrap(),
            "{}",
            what("charge")
        );
        close(
            f64::from(feature.quality),
            f64::from(f32_hex(&row[5])),
            relative,
            &what("quality"),
        );
        assert_eq!(
            feature.quality_rt,
            f32_hex(&row[6]),
            "{}",
            what("quality_rt")
        );
        assert_eq!(
            feature.quality_mz,
            f32_hex(&row[7]),
            "{}",
            what("quality_mz")
        );
        close(
            f64::from(feature.width),
            f64::from(f32_hex(&row[8])),
            relative,
            &what("width"),
        );
        assert_eq!(
            feature.subordinates.len(),
            row[9].parse::<usize>().unwrap(),
            "{}",
            what("subordinates")
        );
        assert_eq!(
            feature.convex_hulls.len(),
            row[10].parse::<usize>().unwrap(),
            "{}",
            what("hull count")
        );
    }
    check_meta(file, config, map);
    check_hulls(file, config, map);
}

fn check_meta(file: &str, config: &str, map: &FeatureMap) {
    let mut expected: BTreeMap<usize, BTreeMap<String, (String, String)>> = BTreeMap::new();
    for row in records(file, "meta", config) {
        expected
            .entry(row[0].parse().unwrap())
            .or_default()
            .insert(row[1].clone(), (row[2].clone(), row[3].clone()));
    }
    for (index, keys) in expected {
        let feature = &map.features[index];
        assert_eq!(
            feature.metadata.len(),
            keys.len(),
            "{config}[{index}]: meta key count, {:?} vs {:?}",
            feature.metadata.keys().collect::<Vec<_>>(),
            keys.keys().collect::<Vec<_>>()
        );
        for (key, (kind, value)) in keys {
            let actual = feature
                .metadata
                .get(&key)
                .unwrap_or_else(|| panic!("{config}[{index}]: missing meta {key}"));
            match (kind.as_str(), actual.data()) {
                ("int", MetaValueData::Integer(got)) => {
                    assert_eq!(got.to_string(), value, "{config}[{index}].{key}");
                }
                ("string", MetaValueData::String(got)) => {
                    assert_eq!(got.as_str(), value, "{config}[{index}].{key}");
                }
                ("double", MetaValueData::Float(got)) => close(
                    *got,
                    f64_hex(&value),
                    tolerance(config, None),
                    &format!("{config}[{index}].{key}"),
                ),
                other => panic!("{config}[{index}].{key}: unexpected {other:?} for {kind}"),
            }
        }
    }
}

fn check_hulls(file: &str, config: &str, map: &FeatureMap) {
    for row in records(file, "hull", config) {
        let index: usize = row[0].parse().unwrap();
        let hull: usize = row[1].parse().unwrap();
        let count: usize = row[2].parse().unwrap();
        let points = map.features[index].convex_hulls[hull].hull_points();
        assert_eq!(points.len(), count, "{config}[{index}] hull {hull}: points");
        for (position, point) in row[3].split(' ').enumerate() {
            let (rt, mz) = point.split_once(',').unwrap();
            // Hull points are input coordinates, so they are bit-identical.
            assert_eq!(
                points[position].rt.to_bits(),
                f64_hex(rt).to_bits(),
                "{config}[{index}] hull {hull} point {position} rt"
            );
            assert_eq!(
                points[position].mz.to_bits(),
                f64_hex(mz).to_bits(),
                "{config}[{index}] hull {hull} point {position} mz"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// The per-seed intermediate state (adapted)
// ---------------------------------------------------------------------------

/// Replay step 3.3 seed by seed, as the C2 driver replays it, and compare the
/// isotope fit, the extended and cropped traces, the fitted model and the
/// quality with the executed library's.
#[test]
fn every_seed_matches_the_executed_intermediate_state() {
    for case in &CASES {
        let stage = SeedStage::run((case.experiment)(), &(case.seeds)(), &(case.parameters)())
            .unwrap()
            .unwrap();
        let settings = stage.settings().clone();
        let spectra = &stage.experiment().spectra;
        let seed_rows = records("b7_seed_records.tsv", "seed", case.config);
        let mut row_index = 0usize;
        for (charge_index, charge_seeds) in stage.charges().iter().enumerate() {
            let overall = OverallScores::new(stage.scores(), charge_index);
            for (index, seed) in charge_seeds.seeds.iter().enumerate() {
                let row = &seed_rows[row_index];
                row_index += 1;
                let what = |field: &str| format!("{} seed {index}: {field}", case.config);
                assert_eq!(row[0].parse::<i32>().unwrap(), charge_seeds.charge);
                assert_eq!(row[1].parse::<usize>().unwrap(), index);
                assert_eq!(
                    row[2].parse::<usize>().unwrap(),
                    seed.spectrum,
                    "{}",
                    what("spectrum")
                );
                assert_eq!(
                    row[3].parse::<usize>().unwrap(),
                    seed.peak,
                    "{}",
                    what("peak")
                );

                let (quality, pattern) = find_best_isotope_fit(
                    spectra,
                    stage.windows(),
                    &settings,
                    *seed,
                    charge_seeds.charge,
                )
                .unwrap();
                close(quality, f64_hex(&row[4]), BITWISE, &what("isotope fit"));
                check_pattern(case.config, index, &pattern);
                if quality < settings.min_isotope_fit {
                    assert_eq!(
                        row[6],
                        "Could not find good enough isotope pattern containing the seed",
                        "{}",
                        what("abort")
                    );
                    continue;
                }

                let mut traces = extend_mass_traces(spectra, overall, &settings, &pattern).unwrap();
                check_traces(case.config, index, "extended", &traces);
                let seed_mz = spectra[seed.spectrum].peaks[seed.peak].mz;
                if !traces.is_valid(seed_mz, settings.trace_tolerance) {
                    assert_eq!(row[6], "Could not extend seed", "{}", what("abort"));
                    continue;
                }
                traces.update_baseline();
                traces.baseline *= 0.75;
                traces.get_mut(traces.max_trace).unwrap().update_maximum();
                check_traces(case.config, index, "fit_input", &traces);

                let mut model = FittedModel::new(
                    settings.rt_shape,
                    TraceFitterParams {
                        max_iteration: i64::from(settings.max_iterations),
                        weighted: false,
                    },
                );
                model.fit(&traces).unwrap();
                check_fitter(case.config, index, &model);
                let cropped =
                    crop_feature(model.as_fitter(), &traces, settings.min_trace_score).unwrap();
                check_traces(case.config, index, "cropped", &cropped);
                match check_feature_quality(model.as_fitter(), &cropped, seed_mz, &settings)
                    .unwrap()
                {
                    QualityOutcome::Rejected(reason) => {
                        assert_eq!(row[5], "false", "{}", what("feature_ok"));
                        assert_eq!(row[6], reason, "{}", what("abort"));
                    }
                    QualityOutcome::Accepted(q) => {
                        assert_eq!(row[5], "true", "{}", what("feature_ok"));
                        assert!(
                            !has_platform_gap(case.config, index),
                            "{}: a seed whose fit departs from the executed Eigen became a \
                             feature; the known gap must never change an output",
                            what("known gap")
                        );
                        let relative = tolerance(case.config, Some(index));
                        close(q.fit_score, f64_hex(&row[8]), relative, &what("fit_score"));
                        close(
                            q.correlation,
                            f64_hex(&row[9]),
                            relative,
                            &what("correlation"),
                        );
                        close(
                            q.final_score,
                            f64_hex(&row[10]),
                            relative,
                            &what("final_score"),
                        );
                    }
                }
            }
        }
        assert_eq!(row_index, seed_rows.len(), "{}: seed count", case.config);
    }
}

fn seed_record(file: &str, config: &str, kind: &str, index: usize) -> Vec<Vec<String>> {
    records(file, kind, config)
        .into_iter()
        .filter(|row| row[1].parse::<usize>().unwrap() == index)
        .map(|row| row[2..].to_vec())
        .collect()
}

fn check_pattern(
    config: &str,
    index: usize,
    pattern: &openms::analysis::feature_finder_picked::helper_structs::IsotopePattern,
) {
    let row = &seed_record("b7_seed_records.tsv", config, "pattern", index)[0];
    let codes: Vec<String> = pattern
        .peak
        .iter()
        .map(|peak| match peak {
            PatternPeak::NotFound => "-1".to_string(),
            PatternPeak::Removed => "-2".to_string(),
            PatternPeak::Found(i) => i.to_string(),
        })
        .collect();
    assert_eq!(
        codes.join(" "),
        row[0],
        "{config} seed {index}: pattern peaks"
    );
    let spectra: Vec<String> = pattern.spectrum.iter().map(usize::to_string).collect();
    assert_eq!(
        spectra.join(" "),
        row[1],
        "{config} seed {index}: pattern spectra"
    );
    for (position, expected) in row[2].split(' ').enumerate() {
        if expected.is_empty() {
            continue;
        }
        close(
            pattern.intensity[position],
            f64_hex(expected),
            BITWISE,
            &format!("{config} seed {index}: pattern intensity {position}"),
        );
    }
    for (position, expected) in row[3].split(' ').enumerate() {
        if expected.is_empty() {
            continue;
        }
        close(
            pattern.mz_score[position],
            f64_hex(expected),
            BITWISE,
            &format!("{config} seed {index}: pattern m/z score {position}"),
        );
    }
}

fn check_traces(config: &str, index: usize, stage: &str, traces: &MassTraces) {
    let summary = seed_record("b7_seed_records.tsv", config, "traces", index)
        .into_iter()
        .find(|row| row[0] == stage)
        .unwrap_or_else(|| panic!("{config} seed {index}: no {stage} record"));
    let what = |field: &str| format!("{config} seed {index} {stage}: {field}");
    assert_eq!(
        traces.len(),
        summary[1].parse::<usize>().unwrap(),
        "{}",
        what("size")
    );
    assert_eq!(
        traces.max_trace,
        summary[2].parse::<usize>().unwrap(),
        "{}",
        what("max_trace")
    );
    assert_eq!(
        traces.peak_count(),
        summary[3].parse::<usize>().unwrap(),
        "{}",
        what("peak count")
    );
    if summary[4] != "none" {
        close(
            traces.baseline,
            f64_hex(&summary[4]),
            BITWISE,
            &what("baseline"),
        );
    }
    for row in seed_record("b7_seed_records.tsv", config, "trace", index) {
        if row[0] != stage {
            continue;
        }
        let position: usize = row[1].parse().unwrap();
        let trace = &traces[position];
        close(
            trace.theoretical_int,
            f64_hex(&row[2]),
            BITWISE,
            &what(&format!("trace {position} theoretical_int")),
        );
        let pairs: Vec<String> = trace
            .peaks
            .iter()
            .map(|peak| format!("{}.{}", peak.spectrum, peak.peak))
            .collect();
        assert_eq!(
            pairs.join(" "),
            row[3],
            "{}",
            what(&format!("trace {position} peaks"))
        );
    }
}

fn check_fitter(config: &str, index: usize, model: &FittedModel) {
    let row = &seed_record("b7_seed_records.tsv", config, "fitter", index)[0];
    let fitter = model.as_fitter();
    let relative = tolerance(config, Some(index));
    let what = |field: &str| format!("{config} seed {index} fit: {field}");
    let expected_shape = match model {
        FittedModel::Gauss(_) => "GaussTraceFitter",
        FittedModel::Egh(_) => "EGHTraceFitter",
    };
    assert_eq!(row[0], expected_shape, "{}", what("shape"));
    close(fitter.center(), f64_hex(&row[1]), relative, &what("center"));
    close(fitter.height(), f64_hex(&row[2]), relative, &what("height"));
    close(fitter.fwhm(), f64_hex(&row[3]), relative, &what("fwhm"));
    close(fitter.area(), f64_hex(&row[4]), relative, &what("area"));
    close(
        fitter.lower_rt_bound(),
        f64_hex(&row[5]),
        relative,
        &what("lower bound"),
    );
    close(
        fitter.upper_rt_bound(),
        f64_hex(&row[6]),
        relative,
        &what("upper bound"),
    );
    match model {
        FittedModel::Egh(egh) => {
            close(egh.sigma(), f64_hex(&row[7]), relative, &what("sigma"));
            close(egh.tau(), f64_hex(&row[8]), relative, &what("tau"));
        }
        FittedModel::Gauss(gauss) => {
            close(gauss.sigma(), f64_hex(&row[7]), relative, &what("sigma"));
        }
    }
}

/// The intended isotope-abundance override against the Linux x86_64 Release
/// library, replayed with the override the source intends (adapted).
///
/// The library cannot compute the intended override, so the driver
/// `../oracle/ffap-sem-completion/drivers/intended_abundance.cpp` runs
/// `FeatureFinderAlgorithmPicked::run` with the changed abundance (the executed
/// result: no seed, no feature), then recomputes step 2.5 with an override
/// distribution that is cleared before its two isotopes are inserted, assigns
/// those windows to the protected `isotope_distributions_`, and replays steps
/// 3.1 to 4 with the library's own protected functions on the library's own
/// arrays (two repetitions at one and four threads, identical). Three
/// configurations: FFC_1 with `abundance_12C` 90 and 99 and with
/// `abundance_14N` 95. This port's default, `AbundanceOverride::Intended`,
/// reproduces the replay: every window (bit for bit, which fixes the
/// override's `f32` weights), the seeds with their pattern and overall scores,
/// the candidate counts, the abort reasons and every feature.
#[test]
fn the_intended_abundance_override_matches_the_adapted_release_replay() {
    const FILE: &str = "intended_abundance.tsv";
    for (config, key, value) in [
        ("ffc1_12C_90", "isotopic_pattern:abundance_12C", 90.0),
        ("ffc1_14N_95", "isotopic_pattern:abundance_14N", 95.0),
        ("ffc1_12C_99", "isotopic_pattern:abundance_12C", 99.0),
    ] {
        let mut parameters = ffc1_parameters();
        set(&mut parameters, key, ParamValue::Float(value));
        // The executed library run with the stray-peak override finds nothing.
        let executed = &records(FILE, "run", config)[0];
        assert_eq!(
            executed.as_slice(),
            [
                "Found 0 seeds for charge 2.",
                "Found 0 feature candidates for charge 2."
            ]
        );
        assert_eq!(records(FILE, "executed", config)[0][0], "0");

        let stage = SeedStage::run(ffc1_input(), &FeatureMap::new(), &parameters)
            .unwrap()
            .unwrap();
        let windows = records(FILE, "window", config);
        let patterns = stage.windows().patterns();
        assert_eq!(patterns.len(), windows.len(), "{config}: windows");
        for (pattern, row) in patterns.iter().zip(&windows) {
            let what = format!("{config}: window {}", row[0]);
            assert_eq!(pattern.len().to_string(), row[1], "{what}");
            assert_eq!(pattern.optional_begin.to_string(), row[2], "{what}");
            assert_eq!(pattern.optional_end.to_string(), row[3], "{what}");
            assert_eq!(pattern.max.to_bits(), f64_hex(&row[4]).to_bits(), "{what}");
            assert_eq!(pattern.trimmed_left.to_string(), row[5], "{what}");
            let bits: Vec<u64> = pattern.intensity.iter().map(|v| v.to_bits()).collect();
            let expected: Vec<u64> = row[6..].iter().map(|v| f64_hex(v).to_bits()).collect();
            assert_eq!(bits, expected, "{what}");
        }
        let seeds = records(FILE, "seed", config);
        let charge = &stage.charges()[0];
        assert_eq!(charge.seeds.len(), seeds.len(), "{config}: seed count");
        for (seed, row) in charge.seeds.iter().zip(&seeds) {
            let scores = stage.scores();
            assert_eq!(
                [
                    seed.spectrum.to_string(),
                    seed.peak.to_string(),
                    format!("{:08x}", seed.intensity.to_bits()),
                    format!(
                        "{:08x}",
                        scores.pattern(0, seed.spectrum).unwrap()[seed.peak].to_bits()
                    ),
                    format!(
                        "{:08x}",
                        scores.overall(0, seed.spectrum).unwrap()[seed.peak].to_bits()
                    ),
                ]
                .as_slice(),
                &row[2..],
                "{config}: seed {}",
                row[1]
            );
        }

        let output = run_with_options(
            ffc1_input(),
            &FeatureMap::new(),
            &parameters,
            &Options {
                threads: Threads::serial(),
                ..Options::default()
            },
        )
        .unwrap();
        let candidates = &records(FILE, "candidates", config)[0];
        assert!(output.log.contains(&format!(
            "Found {} feature candidates for charge {}.",
            candidates[1], candidates[0]
        )));
        let intended = &records(FILE, "intended", config)[0];
        assert!(
            output
                .log
                .contains(&format!("Removed {} overlapping features.", intended[0]))
        );
        assert_eq!(output.features.len().to_string(), intended[2], "{config}");
        let aborts: BTreeMap<String, usize> = records(FILE, "abort", config)
            .into_iter()
            .map(|row| (row[1].clone(), row[0].parse().unwrap()))
            .collect();
        assert_eq!(output.aborts, aborts, "{config}: abort reasons");
        check_features(FILE, config, &output.features);
    }
}

/// The one deliberate divergence from the executed C++: a changed isotope
/// abundance.
///
/// The source builds the override from a default `IsotopeDistribution` that
/// already holds `(0, 1)`, so its patterns grow and FeatureFinderCentroided_1
/// with `abundance_12C = 90` finds no seed, no candidate and no feature (C2
/// `ffap_ffc1_abundance_12C_90`, the `run` row of the fixture). The port
/// computes the *intended* two-isotope override instead (`CPP-247`, lead
/// decision of 2026-09-15), so it does find features. This test states both
/// sides, so the divergence cannot become invisible.
#[test]
fn the_abundance_override_deliberately_differs_from_the_executed_library() {
    let executed = &records("b7_feature_records.tsv", "run", "ffc1_abundance_12C_90")[0];
    assert_eq!(executed[0], "0", "executed seeds");
    assert_eq!(executed[1], "0", "executed candidates");
    assert_eq!(executed[3], "0", "executed features");
    assert_eq!(executed[4], "Found 0 seeds for charge 2.");

    let mut parameters = ffc1_parameters();
    set(
        &mut parameters,
        "isotopic_pattern:abundance_12C",
        ParamValue::Float(90.0),
    );
    let output = run(ffc1_input(), &FeatureMap::new(), &parameters).unwrap();
    assert!(
        !output.features.is_empty(),
        "the intended override should find features where the source's stray (0, 1) peak finds none"
    );
    // Opting out refuses rather than differing.
    assert!(matches!(
        run_with_options(
            ffc1_input(),
            &FeatureMap::new(),
            &parameters,
            &Options {
                abundance_override: AbundanceOverride::Refuse,
                ..Options::default()
            },
        ),
        Err(openms::Error::Unsupported(_))
    ));
}

// ---------------------------------------------------------------------------
// Determinism
// ---------------------------------------------------------------------------

/// The seed loop's result does not depend on the thread count: the source's
/// results are schedule-independent and this port's ordered `map_collect` makes
/// that structural.
#[test]
fn the_seed_loop_is_bit_identical_across_thread_counts() {
    let reference = run_with_options(
        ffc1_input(),
        &FeatureMap::new(),
        &ffc1_parameters(),
        &Options {
            threads: Threads::serial(),
            ..Options::default()
        },
    )
    .unwrap();
    for threads in [2i64, 8] {
        let output = run_with_options(
            ffc1_input(),
            &FeatureMap::new(),
            &ffc1_parameters(),
            &Options {
                threads: Threads::from_cli(threads),
                ..Options::default()
            },
        )
        .unwrap();
        assert_eq!(output.log, reference.log, "{threads} threads: log");
        assert_eq!(output.aborts, reference.aborts, "{threads} threads: aborts");
        assert_eq!(
            output.features.len(),
            reference.features.len(),
            "{threads} threads: count"
        );
        for (actual, expected) in output
            .features
            .features
            .iter()
            .zip(&reference.features.features)
        {
            assert_eq!(actual.rt.to_bits(), expected.rt.to_bits());
            assert_eq!(actual.mz.to_bits(), expected.mz.to_bits());
            assert_eq!(actual.intensity.to_bits(), expected.intensity.to_bits());
            assert_eq!(actual.quality.to_bits(), expected.quality.to_bits());
            assert_eq!(actual.width.to_bits(), expected.width.to_bits());
            assert_eq!(actual.metadata, expected.metadata);
            assert_eq!(actual.convex_hulls, expected.convex_hulls);
            assert_eq!(actual.subordinates, expected.subordinates);
        }
    }
}

// ---------------------------------------------------------------------------
// Overlap resolution (tier 1 through the user-seed configuration, tier 4 here)
// ---------------------------------------------------------------------------

fn hull(points: &[(f64, f64)]) -> ConvexHull2D {
    let points: Vec<Point2D> = points
        .iter()
        .map(|&(rt, mz)| Point2D::new(rt, mz))
        .collect();
    ConvexHull2D::from_points(&points).unwrap()
}

fn feature_with(hulls: Vec<ConvexHull2D>, charge: i32, intensity: f32, quality: f32) -> Feature {
    let mut feature = Feature::new(0.0, 0.0, intensity);
    feature.charge = charge;
    feature.quality = quality;
    feature.convex_hulls = hulls;
    feature
}

/// The four containment and partial-overlap cases of `intersection_`, and its
/// division by the smaller total width.
#[test]
fn intersection_follows_the_source_cases() {
    // bb1 contains bb2: the overlap is bb2's width, the denominator bb2's total.
    let f1 = feature_with(vec![hull(&[(0.0, 100.0), (10.0, 100.0)])], 1, 1.0, 1.0);
    let f2 = feature_with(vec![hull(&[(2.0, 100.0), (6.0, 100.0)])], 1, 1.0, 1.0);
    assert_eq!(intersection(&f1, &f2), 1.0);
    assert_eq!(intersection(&f2, &f1), 1.0);
    // Partial overlap: 8..10 of a width-10 and a width-8 feature.
    let f3 = feature_with(vec![hull(&[(8.0, 100.0), (16.0, 100.0)])], 1, 1.0, 1.0);
    assert_eq!(intersection(&f1, &f3), 2.0 / 8.0);
    assert_eq!(intersection(&f3, &f1), 2.0 / 8.0);
    // Disjoint boxes contribute nothing; touching ones intersect inclusively.
    let f4 = feature_with(vec![hull(&[(20.0, 100.0), (30.0, 100.0)])], 1, 1.0, 1.0);
    assert_eq!(intersection(&f1, &f4), 0.0);
    let f5 = feature_with(vec![hull(&[(10.0, 100.0), (20.0, 100.0)])], 1, 1.0, 1.0);
    assert_eq!(intersection(&f1, &f5), 0.0);
    // Several hulls sum their widths, but only the hulls whose *boxes* intersect
    // contribute: the second hull sits at a different m/z, so the overlap stays
    // one hull's while the denominator is still the smaller feature's total.
    let two = feature_with(
        vec![
            hull(&[(0.0, 100.0), (10.0, 100.0)]),
            hull(&[(0.0, 101.0), (10.0, 101.0)]),
        ],
        1,
        1.0,
        1.0,
    );
    assert_eq!(intersection(&two, &f2), 4.0 / 4.0);
    // Two hulls at the same m/z both overlap, and both are counted.
    let same = feature_with(
        vec![
            hull(&[(0.0, 100.0), (10.0, 100.0)]),
            hull(&[(1.0, 100.0), (11.0, 100.0)]),
        ],
        1,
        1.0,
        1.0,
    );
    assert_eq!(intersection(&same, &f2), (4.0 + 4.0) / 4.0);
}

// ---------------------------------------------------------------------------
// Preserved defects and native refusals (tier 4)
// ---------------------------------------------------------------------------

fn tiny_stage() -> SeedStage {
    // Twelve scans of two peaks each, an isotope pair that survives the default
    // trace search with min_spectra 6.
    let mut spectra = Vec::new();
    for scan in 0..12 {
        let intensity = 100.0 + 10.0 * f32::from(6 - (scan as i8 - 6).abs());
        spectra.push(openms::MSSpectrum {
            rt: 10.0 * f64::from(scan),
            ms_level: 1,
            native_id: format!("scan={scan}"),
            peaks: vec![
                openms::Peak1D::new(500.0, intensity),
                openms::Peak1D::new(500.5, intensity * 0.5),
            ],
            ..openms::MSSpectrum::default()
        });
    }
    let experiment = MSExperiment {
        spectra,
        ..MSExperiment::default()
    };
    let mut parameters = default_parameters().unwrap();
    set(&mut parameters, "intensity:bins", ParamValue::Integer(1));
    set(
        &mut parameters,
        "isotopic_pattern:charge_high",
        ParamValue::Integer(2),
    );
    SeedStage::run(experiment, &FeatureMap::new(), &parameters)
        .unwrap()
        .unwrap()
}

/// `extendMassTraces_` refuses a pattern that matched no peak, where the source
/// dereferences its first entry.
#[test]
fn an_empty_pattern_is_refused_instead_of_dereferenced() {
    let stage = tiny_stage();
    let settings = stage.settings().clone();
    let overall = OverallScores::new(stage.scores(), 0);
    let empty =
        openms::analysis::feature_finder_picked::helper_structs::IsotopePattern::new(3).unwrap();
    assert!(matches!(
        extend_mass_traces(&stage.experiment().spectra, overall, &settings, &empty),
        Err(openms::Error::InvalidValue(_))
    ));
}

/// The resource ceilings of the seed loop are checked before it runs.
#[test]
fn the_seed_loop_ceilings_are_checked_first() {
    let limited = |limits: Limits| {
        run_with_options(
            ffc1_input(),
            &FeatureMap::new(),
            &ffc1_parameters(),
            &Options {
                limits,
                ..Options::default()
            },
        )
    };
    assert!(matches!(
        limited(Limits {
            max_seeds: 4,
            ..Limits::default()
        }),
        Err(openms::Error::InvalidValue(_))
    ));
    assert!(matches!(
        limited(Limits {
            max_seed_work: 10,
            ..Limits::default()
        }),
        Err(openms::Error::InvalidValue(_))
    ));
    // The real workload passes with room to spare.
    assert!(limited(Limits::default()).is_ok());
}

/// A trace whose peaks all lie outside the fitted bounds is dropped, and the
/// source's position rules decide what happens to the traces around it.
#[test]
fn cropping_follows_the_source_position_rules() {
    // Three traces; the model keeps only the middle retention times.
    let make = |offset: f64| {
        let mut trace = MassTrace::default();
        for k in 0..5 {
            trace.peaks.push(TracePeak::new(
                k,
                0,
                offset + f64::from(k as i32),
                500.0,
                100.0,
            ));
        }
        trace.theoretical_int = 1.0;
        trace
    };
    let mut traces = MassTraces::new();
    traces.push(make(0.0));
    traces.push(make(100.0));
    traces.max_trace = 0;
    traces.baseline = 0.0;
    let mut model = FittedModel::new(RtShape::Symmetric, TraceFitterParams::default());
    model.fit(&traces).unwrap();
    let cropped = crop_feature(model.as_fitter(), &traces, 0.5).unwrap();
    // The far trace lies beyond the model's bounds, so it is dropped; it comes
    // after `max_trace`, so the cropping simply stops there.
    assert!(cropped.len() <= 1, "{}", cropped.len());
    assert_eq!(cropped.baseline, traces.baseline);
}

/// The port's label is the source's `MetaInfoRegistry` index 3, the key
/// `label`, and it carries the feature number after the containment pass.
#[test]
fn labels_are_the_feature_numbers_in_order() {
    let output = run(ffc1_input(), &FeatureMap::new(), &ffc1_parameters()).unwrap();
    let mut labels: Vec<i64> = output
        .features
        .features
        .iter()
        .map(|feature| match feature.metadata["label"].data() {
            MetaValueData::Integer(value) => *value,
            other => panic!("label is {other:?}"),
        })
        .collect();
    labels.sort_unstable();
    assert_eq!(labels, (0..8).collect::<Vec<_>>());
    assert!(
        output
            .features
            .features
            .iter()
            .all(|f| f.metadata.contains_key("spectrum_index")
                && f.metadata.contains_key("spectrum_native_id"))
    );
    let _ = MetaValue::from(0i64);
}
