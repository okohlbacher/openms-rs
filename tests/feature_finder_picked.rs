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
//! compared bit for bit on every platform, since both trace fitters call the
//! reference build's glibc `exp` and `log` (lead decision D10); only an EGH
//! area on a host without glibc has a measured bound; see [`tolerance`] and
//! [`area_tolerance`].

#![cfg(all(feature = "mzml", feature = "paramxml", feature = "featurexml"))]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use openms::analysis::feature_finder_picked::algorithm::{
    AbundanceOverride, DegenerateBinStep, Limits, Options, RtShape, default_parameters, run,
    run_with_options,
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
use openms::kernel::{
    ChromatogramPeak, ConvexHull2D, DataArray, Feature, FeatureMap, MSChromatogram, MSExperiment,
    NumericRange, Point2D,
};
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
/// feature) of one configuration: bit for bit, NaN bits included.
///
/// The fixtures are the **Linux x86_64 Release** build
/// (`openms4-release-bc9cc12-c19e494-174b576`, the C2 driver `ffap_stages` run on
/// ibminode06, AMD EPYC 7763, glibc 2.39; `../oracle/ffap-sem-completion`), the
/// reference platform the user chose on 2026-09-15. Since lane B3b the port's
/// Levenberg-Marquardt solver follows that build's Eigen kernels, and since
/// lead decision D10 both trace fitters call that build's glibc `exp` and
/// `log`, ported (the FMA variants `__ieee754_exp_fma` and
/// `__ieee754_log_fma`, equal to the executed library on every probed input,
/// `../oracle/ffap-complete-fix3`). No fitted quantity depends on the host
/// any more, so every platform compares bit for bit, except the area of an
/// asymmetric fit ([`area_tolerance`]).
fn tolerance(_config: &str, _index: Option<usize>) -> f64 {
    BITWISE
}

/// The bound of an EGH fit's area and of the intensity derived from it.
///
/// `EGHTraceFitter::getArea` calls `atan`, whose reference implementation (the
/// IBM Accurate Mathematical Library in glibc) has no licence-clean upstream,
/// so the port calls the host's `atan` on x86_64 Linux with glibc, exact where
/// that library selects the reference's `__atan_fma` (glibc 2.39 on a CPU
/// with FMA, as on the reference node and the gate hosts), and the `libm`
/// crate's elsewhere (lead decision D10's fallback). On those other hosts the
/// bound is [`EGH_ATAN_GAP`], the largest departure measured over these
/// fixtures on macOS arm64.
///
/// That measured maximum is `0.0`, which is [`BITWISE`], so this function
/// relaxes nothing anywhere: the area and the intensity derived from it must
/// be bit for bit the executed value on every host, and a departure fails with
/// [`close`]'s bitwise message. What the measurement bounds is how often the
/// two `atan` implementations may differ at all without being seen here: the
/// crate's `atan` differs from the reference for 1.6 to 6.2 % of the arguments
/// in the fits' range, and only the `float` narrowing of the intensity hides
/// it over these fixtures, so a failure on a non-reference host means that
/// host's `atan` is not the reference's on an argument this file reaches, not
/// that the port regressed. The `libm` crate is pure Rust, so the value is the
/// same on every such host.
fn area_tolerance(config: &str) -> f64 {
    #[cfg(not(all(target_os = "linux", target_env = "gnu", target_arch = "x86_64")))]
    if config == "ffc1_asymmetric" {
        return EGH_ATAN_GAP;
    }
    let _ = config;
    BITWISE
}

/// The largest relative departure of an EGH area or intensity from the Linux
/// capture measured on macOS arm64 over the fixtures of this file (the
/// `libm` crate's `atan`), rounded up at the second significant digit.
///
/// The measurement is `0.0`, that is [`BITWISE`], so the comparison is exact
/// on every host; [`area_tolerance`] says what that does and does not promise.
/// Only a new measurement on a host whose `atan` departs may change it, and
/// only upwards from a re-measured maximum, never to pass a failing run.
#[cfg(not(all(target_os = "linux", target_env = "gnu", target_arch = "x86_64")))]
const EGH_ATAN_GAP: f64 = 0.0;

/// Whether `actual` is `expected` within `tolerance` relative (absolute
/// against an executed zero). A [`BITWISE`] tolerance compares the bits, so
/// the sign of a zero and the payload of a NaN must be the executed ones too.
#[track_caller]
fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    if actual.to_bits() == expected.to_bits() {
        return;
    }
    assert!(
        tolerance != BITWISE,
        "{what}: {actual:e} ({:016x}) is not bit for bit the executed {expected:e} ({:016x})",
        actual.to_bits(),
        expected.to_bits()
    );
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
            area_tolerance(config),
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
    close(
        fitter.area(),
        f64_hex(&row[4]),
        area_tolerance(config),
        &what("area"),
    );
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

/// The `UInt` score-array count `3 + 2 * charge_count` (lead decision D12):
/// every count that wraps is refused whatever the [`Limits`], where the
/// executed build writes out of bounds (`charge_high - charge_low + 1` of
/// `2^31 - 1`, `-1` and `-4` or less: SIGSEGV where the wrapped count is small,
/// `boundary_stage.tsv.gz`) or allocates about 2^32 arrays per spectrum (`-2`,
/// `-3`: `std::bad_alloc`). A count of 0 runs and finds nothing, and the
/// counts in between stay behind the native charge ceiling, which the caller
/// can raise.
#[test]
fn charge_count_wraps_are_refused_whatever_the_limits() {
    let unlimited = Options {
        limits: Limits {
            max_charges: usize::MAX,
            ..Limits::default()
        },
        ..Options::default()
    };
    let run_charges = |low: i64, high: i64, options: &Options| {
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
        run_with_options(ffc1_input(), &FeatureMap::new(), &parameters, options)
    };
    let int_max = i64::from(i32::MAX);
    for (low, high, text) in [
        (
            1,
            int_max,
            "wraps to 1 and the source writes past the arrays",
        ),
        (4, 2, "wraps and the source writes past the arrays"),
        (7, 2, "wraps and the source writes past the arrays"),
        (int_max, 1, "wraps and the source writes past the arrays"),
        (int_max, 498, "wraps and the source writes past the arrays"),
        (5, 2, "wraps to 4294967295 arrays per spectrum"),
        (6, 2, "wraps to 4294967293 arrays per spectrum"),
    ] {
        match run_charges(low, high, &unlimited) {
            Err(openms::Error::InvalidValue(message)) => {
                assert!(message.contains(text), "{low}..{high}: {message}");
            }
            other => panic!("{low}..{high}: {other:?}"),
        }
    }
    // Count 0: the source's three arrays and no charge (executed: no feature).
    let output = run_charges(3, 2, &Options::default()).unwrap();
    assert!(output.features.is_empty());
    // Count 2^31 - 2 does not wrap: the native ceiling refuses it (executed:
    // `std::bad_alloc`).
    match run_charges(2, int_max, &Options::default()) {
        Err(openms::Error::InvalidValue(message)) => {
            assert!(message.contains("charges exceed the limit"), "{message}");
        }
        other => panic!("{other:?}"),
    }
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

// ---------------------------------------------------------------------------
// Non-finite input: the Linux x86_64 Release build (tier 1)
// ---------------------------------------------------------------------------

/// The rows of a stage fixture in the format of
/// `../oracle/ffap-sem-completion/extract/extract_nonfinite.py`.
///
/// `nonfinite_stage.tsv.gz`: the driver `nonfinite_stage` run against
/// `openms4-release-bc9cc12-c19e494-174b576` on 189 modified
/// FeatureFinderCentroided_1 inputs, twice each, identical, five of them also
/// twice at four threads, identical apart from the one-thread abort rows.
/// `sort_mobility_stage.tsv.gz`: its variant with drift times
/// (`nonfinite_stage_dt`) on the sort, drift-time and step-2.5 cases of
/// `../oracle/ffap-complete-fix1/node/run_stage.sh`, checked and written by
/// `../oracle/ffap-complete-fix1/extract/extract_stage.py` in the same way.
fn stage_rows(file: &str) -> Vec<Vec<String>> {
    use std::io::Read;
    let bytes = std::fs::read(data(file)).unwrap();
    let mut text = String::new();
    flate2::read::GzDecoder::new(bytes.as_slice())
        .read_to_string(&mut text)
        .unwrap();
    text.lines()
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect()
}

/// FNV-1a over little-endian bytes, as the driver and the extraction digest.
struct Fnv(u64);

impl Fnv {
    fn new() -> Self {
        Self(0xcbf2_9ce4_8422_2325)
    }

    fn bytes(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.0 ^= u64::from(byte);
            self.0 = self.0.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
}

/// A spectrum index of an option: a number or `last`.
fn spectrum_index(experiment: &MSExperiment, text: &str) -> usize {
    if text == "last" {
        experiment.spectra.len() - 1
    } else {
        text.parse().unwrap()
    }
}

/// The input, parameters and user seeds of one case, built as the driver
/// built them: the mzML (`input=class`: the class test's, loaded as the
/// driver loads every input), `keep=` first, then every modification in
/// option order.
fn nonfinite_input(options: &[String]) -> (MSExperiment, Param, FeatureMap) {
    let mut experiment = if options.iter().any(|option| option == "input=class") {
        let mut load = PeakFileOptions::default();
        load.add_ms_level(1).unwrap();
        load.set_intensity_range(NumericRange {
            min: 0.0,
            max: f64::MAX,
        });
        FileHandler::load_experiment_with_options(
            data("FeatureFinderAlgorithmPicked.mzML"),
            &[FileType::MzMl],
            &load,
        )
        .unwrap()
    } else {
        ffc1_input()
    };
    let mut parameters = ffc1_parameters();
    let mut seeds = FeatureMap::new();
    for option in options {
        if let Some(keep) = option.strip_prefix("keep=") {
            experiment.spectra.truncate(keep.parse().unwrap());
        }
    }
    for option in options {
        let (lhs, value) = option.split_once('=').unwrap();
        let parts: Vec<&str> = lhs.split(':').collect();
        match parts[0] {
            "keep" | "scores" | "input" | "isowin" => {}
            // The whole-input modifications of `fix4_stage`
            // (`../oracle/ffap-complete-fix4`, the round-3 numerics
            // verifier's `v3_stage`), with the driver's arithmetic.
            "rtscale" => {
                let factor = f64_hex(value);
                for spectrum in &mut experiment.spectra {
                    spectrum.rt *= factor;
                }
            }
            "rtshift" => {
                let shift = f64_hex(value);
                for spectrum in &mut experiment.spectra {
                    spectrum.rt += shift;
                }
            }
            "inscale" => {
                let factor = f32_hex(value);
                for spectrum in &mut experiment.spectra {
                    for peak in &mut spectrum.peaks {
                        peak.intensity *= factor;
                    }
                }
            }
            "skew" => {
                let a = f64_hex(value);
                let n = experiment.spectra.len() as f64;
                for (k, spectrum) in experiment.spectra.iter_mut().enumerate() {
                    let factor = (a * k as f64) / n + 1.0;
                    for peak in &mut spectrum.peaks {
                        peak.intensity = (f64::from(peak.intensity) * factor) as f32;
                    }
                }
            }
            "jit" => {
                // splitmix64's finalizer over `(seed << 40) + (k << 20) + q`.
                let seed: u64 = parts[1].parse().unwrap();
                let amplitude = f64_hex(value);
                for (k, spectrum) in experiment.spectra.iter_mut().enumerate() {
                    for (q, peak) in spectrum.peaks.iter_mut().enumerate() {
                        let mut z = (seed << 40)
                            .wrapping_add((k as u64) << 20)
                            .wrapping_add(q as u64);
                        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
                        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
                        z ^= z >> 31;
                        let u = (z >> 11) as f64 * f64::from_bits(0x3ca0_0000_0000_0000);
                        let factor = (u - 0.5) * amplitude + 1.0;
                        peak.intensity = (f64::from(peak.intensity) * factor) as f32;
                    }
                }
            }
            "seeds" => {
                assert_eq!(
                    value.rsplit('/').next().unwrap(),
                    "FeatureFinderCentroided_1_1_output.featureXML"
                );
                seeds = ffc1_user_seeds();
            }
            "seedkeep" => seeds.features.truncate(value.parse().unwrap()),
            "seedrt" => seeds.features[parts[1].parse::<usize>().unwrap()].rt = f64_hex(value),
            "seedmz" => seeds.features[parts[1].parse::<usize>().unwrap()].mz = f64_hex(value),
            "dt" => {
                let s = spectrum_index(&experiment, parts[1]);
                experiment.spectra[s].drift_time = f64_hex(value);
            }
            "rtall" => {
                for spectrum in &mut experiment.spectra {
                    spectrum.rt = f64_hex(value);
                }
            }
            "empty" => {
                let s = spectrum_index(&experiment, value);
                experiment.spectra[s].peaks.clear();
            }
            "rt" => {
                let s = spectrum_index(&experiment, parts[1]);
                experiment.spectra[s].rt = f64_hex(value);
            }
            "trim" => {
                let s = spectrum_index(&experiment, parts[1]);
                experiment.spectra[s].peaks.truncate(value.parse().unwrap());
            }
            "mzall" => {
                let s = spectrum_index(&experiment, parts[1]);
                for peak in &mut experiment.spectra[s].peaks {
                    peak.mz = f64_hex(value);
                }
            }
            "mz" | "in" => {
                let s = spectrum_index(&experiment, parts[1]);
                let peaks = &mut experiment.spectra[s].peaks;
                let p = if parts[2] == "last" {
                    peaks.len() - 1
                } else {
                    parts[2].parse().unwrap()
                };
                if parts[0] == "mz" {
                    peaks[p].mz = f64_hex(value);
                } else {
                    peaks[p].intensity = f32_hex(value);
                }
            }
            "innear" => {
                // `MSSpectrum::findNearest` on the unmodified, sorted target.
                let s = spectrum_index(&experiment, parts[1]);
                let p: usize = parts[2].parse().unwrap();
                let d: isize = parts[3].parse().unwrap();
                let mz = experiment.spectra[s].peaks[p].mz;
                let target = &mut experiment.spectra[s.checked_add_signed(d).unwrap()].peaks;
                let above = target.partition_point(|peak| peak.mz < mz);
                let nearest = if above == 0 {
                    0
                } else if above == target.len() {
                    above - 1
                } else if (target[above].mz - mz).abs() < (target[above - 1].mz - mz).abs() {
                    above
                } else {
                    above - 1
                };
                target[nearest].intensity = f32_hex(value);
            }
            "swap" => {
                let a = spectrum_index(&experiment, parts[1]);
                let b = spectrum_index(&experiment, value);
                experiment.spectra.swap(a, b);
            }
            "rev" => {
                let s = spectrum_index(&experiment, parts[1]);
                experiment.spectra[s].peaks.reverse();
            }
            "fa" | "sa" | "ia" => {
                // `fix3_stage.cpp`: array `i` of spectrum `s` with `n` entries,
                // `n` a count or `size`, `size+k`, `size-k`; lower arrays are
                // created empty.
                let s = spectrum_index(&experiment, parts[1]);
                let i: usize = parts[2].parse().unwrap();
                let spectrum = &mut experiment.spectra[s];
                let n = match value.strip_prefix("size") {
                    Some("") => spectrum.peaks.len(),
                    Some(delta) => spectrum
                        .peaks
                        .len()
                        .checked_add_signed(delta.parse::<isize>().unwrap())
                        .unwrap(),
                    None => value.parse().unwrap(),
                };
                match parts[0] {
                    "fa" => {
                        let arrays = &mut spectrum.float_data_arrays;
                        while arrays.len() <= i {
                            arrays.push(DataArray::new("", Vec::new()));
                        }
                        arrays[i] = DataArray::new(
                            format!("fa{i}"),
                            (0..n).map(|k| 0.5 * k as f32).collect(),
                        );
                    }
                    "sa" => {
                        let arrays = &mut spectrum.string_data_arrays;
                        while arrays.len() <= i {
                            arrays.push(DataArray::new("", Vec::new()));
                        }
                        arrays[i] = DataArray::new(
                            format!("sa{i}"),
                            (0..n).map(|k| format!("s{k}")).collect(),
                        );
                    }
                    _ => {
                        let arrays = &mut spectrum.integer_data_arrays;
                        while arrays.len() <= i {
                            arrays.push(DataArray::new("", Vec::new()));
                        }
                        arrays[i] =
                            DataArray::new(format!("ia{i}"), (0..n).map(|k| k as i32).collect());
                    }
                }
            }
            "chrom" => {
                // A chromatogram of `n` peaks with descending retention times
                // `n - 1 .. 0`, intensities `1 ..= n`, product m/z 500.
                let n: usize = value.parse().unwrap();
                let mut chromatogram = MSChromatogram::default();
                chromatogram.product.mz = 500.0;
                chromatogram.peaks = (0..n)
                    .map(|k| ChromatogramPeak {
                        rt: (n - 1 - k) as f64,
                        intensity: (k + 1) as f32,
                    })
                    .collect();
                experiment.chromatograms.push(chromatogram);
            }
            "cfa" => {
                let n: usize = value.parse().unwrap();
                let chromatogram = experiment.chromatograms.last_mut().unwrap();
                let arrays = &mut chromatogram.float_data_arrays;
                if arrays.is_empty() {
                    arrays.push(DataArray::new("", Vec::new()));
                }
                arrays[0] = DataArray::new("cfa0", (0..n).map(|k| 0.5 * k as f32).collect());
            }
            "i" => set(
                &mut parameters,
                &lhs[2..],
                ParamValue::Integer(value.parse().unwrap()),
            ),
            "d" => set(
                &mut parameters,
                &lhs[2..],
                ParamValue::Float(value.parse().unwrap()),
            ),
            "s" => set(&mut parameters, &lhs[2..], ParamValue::String(value.into())),
            other => panic!("unknown option {other}"),
        }
    }
    (experiment, parameters, seeds)
}

/// Replay every case of a stage fixture ([`stage_rows`]) and return how many
/// cases ended in each outcome.
///
/// For every case the port gives the executed outcome:
///
/// - the executed exception, as [`openms::Error`] with the same `what()` text:
///   `-inf` m/z fails the positive-m/z check after the sort; `+inf` or `1e300`
///   m/z leaves no isotope window (the `Size` conversion of `ceil(inf) + 1` is
///   0), and a NaN m/z asks for window `2^63`, both at the first pattern
///   lookup of step 3.1; every retention time or every m/z NaN leaves an empty
///   range. A window count above `vector::max_size()` =
///   164,703,072,086,692,425 makes the source's `resize` throw
///   `std::length_error`, whose text the port returns; the executed
///   `std::bad_alloc` of a count just below it is where the port's native
///   window ceiling refuses instead;
/// - status 137, the executed run killed after 30 s: a NaN retention time in a
///   mass trace makes `computeIntensityProfile` loop forever
///   (`FeatureFinderAlgorithmPickedHelperStructs.cpp:210-236`), and the port
///   refuses at exactly that merge;
/// - otherwise the printed lines, the feature count, the abort reasons, the
///   bin steps and the window count, the seeds, a digest of every quantile
///   and every per-peak score (NaN bits included; the full rows where few
///   differ from the unmodified input), and every feature with its meta values
///   and a digest of its convex hulls. Every overall score is the executed
///   one, including those the Release build's `powf` rounds one binary32 step
///   away from the correctly rounded value (the `rounding` rows, `CPP-272`).
///   Fitted quantities use the platform bound of [`tolerance`].
fn replay_stage_fixture(rows: &[Vec<String>]) -> BTreeMap<&'static str, usize> {
    use openms::Error;
    use openms::analysis::feature_finder_picked::algorithm::feature_stage;
    let cases: Vec<&Vec<String>> = rows.iter().filter(|row| row[0] == "case").collect();
    let mut outcomes: BTreeMap<&'static str, usize> = BTreeMap::new();
    for case in cases {
        let name = case[1].as_str();
        let of = |kind: &str| -> Vec<&[String]> {
            rows.iter()
                .filter(|row| row[0] == kind && row[1] == name)
                .map(|row| &row[2..])
                .collect()
        };
        let (experiment, parameters, seeds) = nonfinite_input(&case[2..]);
        let input = of("input")[0];
        assert_eq!(experiment.spectra.len().to_string(), input[0], "{name}");
        // `MSExperiment::getSize` counts chromatogram peaks too.
        let peaks: usize = experiment
            .spectra
            .iter()
            .map(|s| s.peaks.len())
            .sum::<usize>()
            + experiment
                .chromatograms
                .iter()
                .map(|c| c.peaks.len())
                .sum::<usize>();
        assert_eq!(peaks.to_string(), input[1], "{name}");
        let options = Options {
            threads: Threads::serial(),
            ..Options::default()
        };
        let stage = SeedStage::run_with_options(experiment.clone(), &seeds, &parameters, &options);
        let status = of("status")[0][0].as_str();
        let rt_config = if case[2..]
            .iter()
            .any(|o| o == "s:feature:rt_shape=asymmetric")
        {
            "ffc1_asymmetric"
        } else {
            "ffc1_symmetric"
        };
        if status == "137" {
            let error = stage
                .and_then(|stage| feature_stage(&stage.unwrap(), &options))
                .unwrap_err();
            assert!(
                matches!(&error, Error::InvalidValue(m) if m.contains("NaN retention time cannot be merged")),
                "{name}: {error}"
            );
            *outcomes.entry("hang").or_default() += 1;
            continue;
        }
        if status == "139" {
            // SIGSEGV, twice: the port refuses where the source reads or writes
            // out of bounds, at an empty best isotope pattern
            // (`extendMassTraces_`) or at the wrapped score-array count.
            let error = stage
                .and_then(|stage| feature_stage(&stage.unwrap(), &options))
                .unwrap_err();
            assert!(
                matches!(&error, Error::InvalidValue(m)
                    if m.contains("the isotope pattern matched no peak")
                        || (m.contains("score-array count") && m.contains("undefined"))),
                "{name}: {error}"
            );
            *outcomes.entry("crash").or_default() += 1;
            continue;
        }
        assert_eq!(status, "0", "{name}");
        if let Some(threw) = of("threw").first() {
            let (kind, text) = threw[0].split_once(": ").unwrap();
            let error = stage
                .and_then(|stage| feature_stage(&stage.unwrap(), &options))
                .unwrap_err();
            match (kind, text, &error) {
                ("IllegalArgument" | "InvalidValue", _, Error::InvalidValue(message))
                | ("InvalidRange", _, Error::InvalidRange(message))
                | ("std::exception", "vector::_M_default_append", Error::InvalidValue(message)) => {
                    assert_eq!(message, text, "{name}");
                }
                ("Precondition failed", _, Error::InvalidValue(message)) => {
                    // `MSSpectrum::sort` or `MSChromatogram::sort` of an
                    // unsorted input with a mis-sized data array.
                    assert_eq!(message, text, "{name}");
                }
                ("std::exception", "std::bad_alloc", Error::InvalidValue(message)) => {
                    // Below `vector::max_size()` the source allocates, which
                    // fails or not depending on memory; the port's window
                    // ceiling refuses first. So does its charge ceiling in
                    // front of `3 + 2 * charge_count` score arrays per
                    // spectrum, and a wrapped count near 2^32 is refused
                    // whatever the limits.
                    assert!(
                        (message.contains("isotope windows")
                            && message.contains("exceed the limit"))
                            || (message.contains("charges exceed the limit"))
                            || (message.contains("score-array count wraps")
                                && message.contains("depends on memory"))
                            || (message.contains("score-array count wraps")
                                && message.contains("undefined")),
                        "{name}: {message}"
                    );
                }
                _ => panic!("{name}: executed {kind}: {text}, port {error}"),
            }
            *outcomes.entry("threw").or_default() += 1;
            continue;
        }
        let stage = stage
            .unwrap_or_else(|error| panic!("{name}: {error}"))
            .unwrap();
        let output =
            feature_stage(&stage, &options).unwrap_or_else(|error| panic!("{name}: {error}"));

        let printed: Vec<&String> = output
            .log
            .iter()
            .filter(|line| line.starts_with("Found "))
            .collect();
        let stdout: Vec<&String> = of("stdout").into_iter().map(|row| &row[0]).collect();
        assert_eq!(printed, stdout, "{name}: printed lines");
        assert_eq!(
            output.features.len().to_string(),
            of("features")[0][0],
            "{name}: features"
        );
        let aborts: BTreeMap<String, usize> = of("abort")
            .into_iter()
            .map(|row| (row[1].clone(), row[0].parse().unwrap()))
            .collect();
        assert_eq!(output.aborts, aborts, "{name}: abort reasons");

        // Step 1 and step 2.5.
        let thresholds = stage.thresholds();
        let bins = of("bins")[0];
        assert_eq!(thresholds.bins().to_string(), bins[0], "{name}");
        for (value, expected) in [
            thresholds.rt_start(),
            thresholds.mz_start(),
            thresholds.rt_step(),
            thresholds.mz_step(),
        ]
        .iter()
        .zip(&bins[1..])
        {
            assert_eq!(value.to_bits(), f64_hex(expected).to_bits(), "{name}: bins");
        }
        assert_eq!(
            stage.windows().patterns().len().to_string(),
            of("windows")[0][0],
            "{name}: windows"
        );
        // `isowin` rows: the executed `isotope_distributions_` entries.
        for row in of("isowin") {
            let index: usize = row[0].parse().unwrap();
            let pattern = &stage.windows().patterns()[index];
            let mut actual = vec![
                pattern.trimmed_left.to_string(),
                pattern.optional_begin.to_string(),
                pattern.optional_end.to_string(),
                format!("{:016x}", pattern.max.to_bits()),
                pattern.intensity.len().to_string(),
            ];
            actual.extend(
                pattern
                    .intensity
                    .iter()
                    .map(|value| format!("{:016x}", value.to_bits())),
            );
            assert_eq!(
                actual.as_slice(),
                &row[1..],
                "{name}: isotope window {index}"
            );
        }

        // Every score, through the digest and the listed rows.
        let scores = stage.scores();
        let charges = scores.charge_count();
        let spectra = &stage.experiment().spectra;
        let arrays_of = |s: usize, p: usize| -> Vec<u32> {
            let mut values = vec![
                scores.trace(s).unwrap()[p].to_bits(),
                scores.intensity(s).unwrap()[p].to_bits(),
                scores.local_max(s).unwrap()[p].to_bits(),
            ];
            for c in 0..charges {
                values.push(scores.pattern(c, s).unwrap()[p].to_bits());
            }
            for c in 0..charges {
                values.push(scores.overall(c, s).unwrap()[p].to_bits());
            }
            values
        };
        // The executed `powf` misrounds these scores; the port computes the
        // same misrounded values.
        for row in of("rounding") {
            let (s, p, c): (usize, usize, usize) = (
                row[0].parse().unwrap(),
                row[1].parse().unwrap(),
                row[2].parse().unwrap(),
            );
            let executed = u32::from_str_radix(&row[3], 16).unwrap();
            let correct = u32::from_str_radix(&row[4], 16).unwrap();
            assert_eq!(executed.abs_diff(correct), 1, "{name}: rounding row");
            assert_eq!(
                scores.overall(c, s).unwrap()[p].to_bits(),
                executed,
                "{name}: overall score s{s} p{p} c{c}"
            );
        }
        for row in of("score") {
            let s: usize = row[0].parse().unwrap();
            let p: usize = row[1].parse().unwrap();
            assert_eq!(
                spectra[s].peaks[p].mz.to_bits(),
                f64_hex(&row[2]).to_bits(),
                "{name}"
            );
            let expected: Vec<u32> = row[4..row.len() - 1]
                .iter()
                .map(|v| u32::from_str_radix(v, 16).unwrap())
                .collect();
            assert_eq!(arrays_of(s, p), expected, "{name}: spectrum {s} peak {p}");
            let score = thresholds
                .score(
                    spectra[s].rt,
                    spectra[s].peaks[p].mz,
                    f64::from(spectra[s].peaks[p].intensity),
                )
                .unwrap();
            assert_eq!(
                score.to_bits(),
                f64_hex(&row[row.len() - 1]).to_bits(),
                "{name}: intensityScore_({s}, {p})"
            );
        }
        for row in of("quantiles") {
            let actual = thresholds
                .quantiles(row[0].parse().unwrap(), row[1].parse().unwrap())
                .unwrap();
            let expected: Vec<u64> = row[2..].iter().map(|q| f64_hex(q).to_bits()).collect();
            let actual: Vec<u64> = actual.iter().map(|q| q.to_bits()).collect();
            assert_eq!(actual, expected, "{name}: quantiles");
        }

        let mut digest = Fnv::new();
        for rt in 0..thresholds.bins() {
            for mz in 0..thresholds.bins() {
                for q in thresholds.quantiles(rt, mz).unwrap() {
                    digest.bytes(&q.to_bits().to_le_bytes());
                }
            }
        }
        for (s, spectrum) in spectra.iter().enumerate() {
            for p in 0..spectrum.peaks.len() {
                for value in arrays_of(s, p) {
                    digest.bytes(&value.to_le_bytes());
                }
            }
        }
        assert_eq!(
            format!("{:016x}", digest.0),
            of("digest")[0][0],
            "{name}: score digest"
        );

        // Seeds (automatic seeds only), with the executed overall score.
        if seeds_are_automatic(&case[2..]) {
            let mut actual = Vec::new();
            for charge in stage.charges() {
                let index = (charge.charge - stage.settings().charge_low) as usize;
                for (rank, seed) in charge.seeds.iter().enumerate() {
                    let overall = arrays_of(seed.spectrum, seed.peak)[3 + charges + index];
                    actual.push(vec![
                        charge.charge.to_string(),
                        rank.to_string(),
                        seed.spectrum.to_string(),
                        seed.peak.to_string(),
                        format!("{:08x}", seed.intensity.to_bits()),
                        format!("{overall:08x}"),
                    ]);
                }
            }
            let expected = of("seed");
            assert_eq!(actual.len(), expected.len(), "{name}: seed count");
            for (a, e) in actual.iter().zip(&expected) {
                assert_eq!(a.as_slice(), *e, "{name}: seed");
            }
        }

        check_nonfinite_features(name, rt_config, &of, &output.features);
        *outcomes.entry("features").or_default() += 1;

        // `DegenerateBinStep::Refuse` refuses exactly the runs whose executed
        // bin steps are zero or infinite (an infinite coordinate makes them
        // so) and whose seed loop visits a scan; the others run unchanged.
        let degenerate = [&bins[3], &bins[4]].iter().any(|step| {
            let step = f64_hex(step);
            step == 0.0 || step.is_infinite()
        });
        let scans = spectra.len();
        let min_spectra = stage.settings().min_spectra;
        let read = min_spectra < scans - min_spectra.min(scans);
        let refusing = SeedStage::run_with_options(
            experiment,
            &seeds,
            &parameters,
            &Options {
                degenerate_bin_step: DegenerateBinStep::Refuse,
                ..options
            },
        );
        if degenerate && read {
            assert!(
                matches!(&refusing, Err(Error::InvalidValue(m)) if m.contains("DegenerateBinStep::Refuse")),
                "{name}"
            );
            *outcomes
                .entry("refused under DegenerateBinStep::Refuse")
                .or_default() += 1;
        } else {
            let refusing = refusing
                .unwrap_or_else(|error| panic!("{name}: {error}"))
                .unwrap();
            assert_eq!(refusing.log(), stage.log(), "{name}");
        }
    }
    outcomes
}

/// Infinite and NaN retention times, m/z values, intensities and user-seed
/// positions, against the executed Linux x86_64 Release build
/// (`nonfinite_stage.tsv.gz`, [`replay_stage_fixture`]).
///
/// Infinite retention times make every step infinite and every intensity score
/// NaN (no feature); infinite intensities shift the quantiles, score NaN at
/// their own peak, join mass traces and are cut off again by the slope check,
/// as in the source. NaN keys of the source's sorts are sorted as the Release
/// build's libstdc++ sorts them (the introsort of `std::sort`, the merge sort
/// of `std::stable_sort`).
#[test]
fn non_finite_inputs_match_the_linux_release_build() {
    let rows = stage_rows("nonfinite_stage.tsv.gz");
    assert_eq!(rows.iter().filter(|row| row[0] == "case").count(), 189);
    let outcomes = replay_stage_fixture(&rows);
    // Of the 170 executed runs that returned features, all are reproduced,
    // among them the nine whose sorts see a NaN key (a NaN intensity in a
    // step-1 cell, and `seeds_mz_nan`); the 16 that threw are reproduced; the
    // 3 that never returned are refused at the endless profile merge, where
    // `rt_nan_mid_unsorted` gets after its NaN retention-time sort. Of the
    // returned runs, the opt-out refuses the seven with an infinite step and a
    // non-empty seed loop: `rt_posinf_last`, `rt_posinf_last_min0`,
    // `rt_neginf_first`, `rt_posinf_all`, `rt_neginf_all`,
    // `rt_posinf_last_bins3` and `mz_posinf_last_nocharge`; `rt_neginf_short`
    // (10 scans, `min_spectra_` 7) has an empty seed loop.
    assert_eq!(
        outcomes,
        BTreeMap::from([
            ("features", 170),
            ("hang", 3),
            ("refused under DegenerateBinStep::Refuse", 7),
            ("threw", 16)
        ])
    );
}

/// The sorts, the drift-time filter and the step-2.5 bound, against the
/// executed Linux x86_64 Release build (`sort_mobility_stage.tsv.gz`,
/// [`replay_stage_fixture`]):
///
/// - `std::sort` of a step-1 cell holding NaN intensities, NaN bits of both
///   signs and signed zeros, whether or not the keys are strictly weakly
///   ordered (`v2_cell_*`, `v3_cell_*`, `v3_ms1_nan_mixed_40`);
/// - `MSExperiment::sortSpectra` of an unsorted input whose retention times
///   tie or hold a NaN (`v2_rt_*`, `v3_rt_*`: `v2_rt_tie3_unsorted` finds 26
///   seeds only in the introsort order; `v3_rt_nan_first_unsorted` never
///   returns);
/// - `MSSpectrum::sortByPosition`'s `std::stable_sort` of peaks holding NaN
///   m/z values (`v3_mz_nan_*`; the `nocharge` cases show the order in the
///   trace scores, the others end in the NaN window lookup);
/// - `FeatureMap::sortByMZ` of user seeds with equal and NaN m/z values
///   (`v2_seeds_*`, `v3_seeds_*`);
/// - the area iterator's drift-time filter: NaN, `+inf` and `-inf` drift times
///   leave their scan out of every step-1 cell, `f64::MAX`, `f64::MIN` and
///   finite ones keep it (`v2_dt_*`, `v3_dt_*`);
/// - window counts around `vector::max_size()` (`v2_mz_2p64`,
///   `v2_mz_below_2p64`, `v3_mz_count_above_max_size`,
///   `v3_mz_count_below_max_size`) and a `-0.0` m/z (`v2_mz_negzero`).
#[test]
fn sort_and_mobility_cases_match_the_linux_release_build() {
    let rows = stage_rows("sort_mobility_stage.tsv.gz");
    assert_eq!(rows.iter().filter(|row| row[0] == "case").count(), 52);
    let outcomes = replay_stage_fixture(&rows);
    // 43 executed runs returned, 8 threw (four NaN window lookups, three
    // length errors, one allocation failure) and `v3_rt_nan_first_unsorted`
    // never returned. The two sixteen-scan inputs with one retention time
    // have a zero step and a non-empty seed loop.
    assert_eq!(
        outcomes,
        BTreeMap::from([
            ("features", 43),
            ("hang", 1),
            ("refused under DegenerateBinStep::Refuse", 2),
            ("threw", 8)
        ])
    );
}

/// Boundaries of the combined fix round 3, against the executed Linux x86_64
/// Release build (`boundary_stage.tsv.gz`, `../oracle/ffap-complete-fix3`,
/// [`replay_stage_fixture`]):
///
/// - averagine windows whose 20 binary32 bins all underflow, at width 100
///   from window 2738 (centre 273,850 Da; the executed mass boundary lies
///   between 273,769.5 and 273,770.5 Da, see
///   [`extended_cases_match_the_linux_release_build`]) on: the source's NaN
///   weights empty those windows and the run
///   continues (`u_mz136850` has no such window, `u_mz136850_5` the first;
///   `u_mz100k` to `u_mz1e6`, `u_ch1000_keep20`, with `seed:min_score` 0, the
///   EGH fit and a zero cutoff);
/// - a NaN `isotopic_pattern:intensity_percentage_optional`, which empties
///   every window (`p_ipo_nan*`), and the controls at 100 and `-0.0`;
/// - EGH and Gaussian fits over 34 configurations, compared bit for bit on
///   every platform (lead decision D10);
/// - the empty best isotope pattern of `extendMassTraces_` at
///   `feature:min_isotope_fit` 0, where the source crashes (`g_avg_trace0`,
///   `g_iso0_seed0`, `p_ipo_*_seed0_iso0`) and a bound of `1e-300` returns;
/// - the `UInt` count `3 + 2 * charge_count` (lead decision D12): SIGSEGV at
///   the wraps that write out of bounds, `std::bad_alloc` where the wrapped
///   count is near 2^32 arrays per spectrum, and the count 0 of
///   `charge_low = charge_high + 1`;
/// - unsorted input with mis-sized float, string and integer data arrays and a
///   chromatogram: the source's `Exception::Precondition` text for the first
///   such spectrum in its introsort order, then chromatograms; sorted
///   spectra and exact arrays run through.
#[test]
fn boundary_cases_match_the_linux_release_build() {
    let rows = stage_rows("boundary_stage.tsv.gz");
    assert_eq!(rows.iter().filter(|row| row[0] == "case").count(), 85);
    let outcomes = replay_stage_fixture(&rows);
    // 63 runs returned (`a_ties_exact`, every retention time equal, has a
    // zero step and a non-empty seed loop), 8 crashed and 14 threw: four
    // `std::bad_alloc` of the wrapped charge count and ten
    // `Exception::Precondition`.
    assert_eq!(
        outcomes,
        BTreeMap::from([
            ("crash", 8),
            ("features", 63),
            ("refused under DegenerateBinStep::Refuse", 1),
            ("threw", 14)
        ])
    );
}

/// Cases beyond the earlier fixtures, against the executed Linux x86_64
/// Release build (`extended_stage.tsv.gz`, `../oracle/ffap-complete-fix4`,
/// the round-3 numerics verifier's case list re-executed, every case twice
/// and identical, and equal to the verifier's own capture;
/// [`replay_stage_fixture`]):
///
/// - EGH and Gaussian fits on jittered, skewed, rescaled and shifted inputs,
///   tiny, huge and infinite intensities, the class test's input, user seeds,
///   `fit:max_iterations` 2 and 100, a NaN intensity and a NaN drift time
///   (`ve_*`, `vg_*`), bit for bit on every platform;
/// - isotope windows at 1-Da resolution across the binary32 underflow of the
///   averagine bins (`vw_*`: the window centred at 273,769.5 Da keeps one
///   bin, 273,770.5 Da is the first empty one; at width 100, windows 2734 to
///   2737 hold a single bin and 2738 is empty), and under NaN, tiny and 100 %
///   cutoffs (`isowin` rows);
/// - both sides of every charge-count wrap (`vc_*`, lead decision D12) and
///   of the empty best pattern (`vi_*`: `feature:min_isotope_fit` 0, `-0.0`
///   and NaN crash, `5e-324` returns);
/// - the `Exception::Precondition` of mis-sized data arrays under all-equal
///   retention times and 20 chromatograms of equal product m/z (`vp_*`,
///   the introsort's order of equal keys);
/// - extreme retention-time and intensity scales (`vx_*`, `vy_*`). Once the
///   fitted `sigma` passes about `1.44e38` the `float` width overflows: the
///   Release build returns those features with an infinite width, `FWHM` meta
///   value and intensity (7 of 9 at `vy_rt_1e37`, all from `vy_rt_1e38` to
///   `vy_rt_1e150`, `vx_rt_1e150` to `vx_rt_1e300`, the FFC_1 parameters
///   unchanged or all thresholds 0), and so does the port, whose
///   [`openms::kernel::BaseFeature::validate`] then refuses them;
///   `vy_rt_1e36` still has finite widths and infinite intensities, and
///   `vy_rt_1e33` finite intensities (the onset in detail:
///   [`width_onset_cases_match_the_linux_release_build`]).
#[test]
fn extended_cases_match_the_linux_release_build() {
    let rows = stage_rows("extended_stage.tsv.gz");
    assert_eq!(rows.iter().filter(|row| row[0] == "case").count(), 134);
    let outcomes = replay_stage_fixture_in_parallel(&rows);
    // 8 crashed (SIGSEGV: three wrapped charge counts, five empty best
    // patterns); 9 threw (five `std::bad_alloc` of the charge count, one of
    // them the native ceiling's, and four `Exception::Precondition`); the
    // other 117 returned, none of them with a zero or infinite bin step.
    assert_eq!(
        outcomes,
        BTreeMap::from([("crash", 8), ("features", 117), ("threw", 9)])
    );

    // The infinite widths, as stored.
    let mut infinite = 0;
    for name in ["vy_rt_1e37", "vy_rt_1e39", "vy_rt_1e39_egh", "vx_rt_1e300"] {
        let case = rows
            .iter()
            .find(|row| row[0] == "case" && row[1] == name)
            .unwrap();
        let (experiment, parameters, seeds) = nonfinite_input(&case[2..]);
        let options = Options {
            threads: Threads::serial(),
            ..Options::default()
        };
        let stage = SeedStage::run_with_options(experiment, &seeds, &parameters, &options)
            .unwrap()
            .unwrap();
        let output =
            openms::analysis::feature_finder_picked::algorithm::feature_stage(&stage, &options)
                .unwrap();
        for feature in &output.features.features {
            if feature.width.is_infinite() {
                infinite += 1;
                let fwhm = feature.metadata.get("FWHM").unwrap();
                assert_eq!(fwhm.as_f64().unwrap(), f64::INFINITY, "{name}");
                assert!(fwhm.validate().is_err(), "{name}");
                assert!(feature.validate().is_err(), "{name}");
            }
        }
    }
    assert!(infinite >= 30, "{infinite}");
}

/// The onset of the `float` width overflow on FeatureFinderCentroided_1
/// (`width_onset_stage.tsv.gz`, `../oracle/ffap-complete-fix5`, `run_onset.sh`:
/// the round-4 numerics verifier's cases re-executed with `fix4_stage`, every
/// case twice and identical; [`replay_stage_fixture`]). With every retention
/// time scaled and the FFC_1 parameters unchanged, the executed features have
/// infinite widths (and `FWHM` values) as follows, of 9 Gaussian and 8 EGH
/// features: none at `2e36` and `4e36`; 1 and 0 at `6e36`; 3 and 2 at `8e36`;
/// all from `1.5e37` (`extended_stage.tsv.gz` adds 7 and 6 at `1e37`, and
/// finite widths at `1e36`). Every intensity is already infinite from `1e36`
/// on and finite at `1e33`. Three jittered inputs at `1e37` (`seed:min_score`
/// 0) have 12 of 14, 14 of 15 and 17 of 17 infinite Gaussian widths and 12 of
/// 14, 12 of 14 and 16 of 16 EGH ones. The port reproduces every case bit for
/// bit.
#[test]
fn width_onset_cases_match_the_linux_release_build() {
    let rows = stage_rows("width_onset_stage.tsv.gz");
    assert_eq!(rows.iter().filter(|row| row[0] == "case").count(), 21);
    let outcomes = replay_stage_fixture_in_parallel(&rows);
    assert_eq!(outcomes, BTreeMap::from([("features", 21)]));

    // The executed counts the documentation states: features, infinite
    // intensities (column 5) and infinite widths (column 10), as bits.
    let infinite = f32::INFINITY.to_bits();
    let count = |case: &str| {
        let features: Vec<&Vec<String>> = rows
            .iter()
            .filter(|row| row[0] == "feature" && row[1] == case)
            .collect();
        let inf = |column: usize| {
            features
                .iter()
                .filter(|row| u32::from_str_radix(&row[column], 16).unwrap() == infinite)
                .count()
        };
        (features.len(), inf(5), inf(10))
    };
    for (case, expected) in [
        ("nb_rt_2e36", (9, 9, 0)),
        ("nb_rt_2e36_egh", (8, 8, 0)),
        ("nb_rt_4e36", (9, 9, 0)),
        ("nb_rt_4e36_egh", (8, 8, 0)),
        ("nb_rt_6e36", (9, 9, 1)),
        ("nb_rt_6e36_egh", (8, 8, 0)),
        ("nb_rt_8e36", (9, 9, 3)),
        ("nb_rt_8e36_egh", (8, 8, 2)),
        ("nb_rt_1p5e37", (9, 9, 9)),
        ("nb_rt_1p5e37_egh", (8, 8, 8)),
        ("nb_rt_3e37", (9, 9, 9)),
        ("nb_rt_3e37_egh", (8, 8, 8)),
        ("nb_rt_5e37", (9, 9, 9)),
        ("nb_rt_5e37_egh", (8, 8, 8)),
        ("nb_rt_1e37_jit81", (14, 14, 12)),
        ("nb_rt_1e37_jit81_egh", (14, 14, 12)),
        ("nb_rt_1e37_jit82", (15, 15, 14)),
        ("nb_rt_1e37_jit82_egh", (14, 14, 12)),
        ("nb_rt_1e37_jit83", (17, 17, 17)),
        ("nb_rt_1e37_jit83_egh", (16, 16, 16)),
    ] {
        assert_eq!(count(case), expected, "{case}");
    }
    let extended = stage_rows("extended_stage.tsv.gz");
    let count_extended = |case: &str| {
        let features: Vec<&Vec<String>> = extended
            .iter()
            .filter(|row| row[0] == "feature" && row[1] == case)
            .collect();
        let inf = |column: usize| {
            features
                .iter()
                .filter(|row| u32::from_str_radix(&row[column], 16).unwrap() == infinite)
                .count()
        };
        (features.len(), inf(5), inf(10))
    };
    for (case, expected) in [
        ("vy_rt_1e33", (9, 0, 0)),
        ("vy_rt_1e33_egh", (8, 0, 0)),
        ("vy_rt_1e36", (9, 9, 0)),
        ("vy_rt_1e36_egh", (8, 8, 0)),
        ("vy_rt_1e37", (9, 9, 7)),
        ("vy_rt_1e37_egh", (8, 8, 6)),
    ] {
        assert_eq!(count_extended(case), expected, "{case}");
    }
}

/// [`replay_stage_fixture`] over the cases of `rows`, spread over worker
/// threads: every case is replayed and asserted exactly as in one call, and
/// the outcome counts are summed. The cases are independent (each builds its
/// own input and runs serially), and a few of them compute a quarter of a
/// million isotope windows, which a single thread of an unoptimised test
/// build takes minutes for.
fn replay_stage_fixture_in_parallel(rows: &[Vec<String>]) -> BTreeMap<&'static str, usize> {
    let names: Vec<&str> = rows
        .iter()
        .filter(|row| row[0] == "case")
        .map(|row| row[1].as_str())
        .collect();
    let workers = std::thread::available_parallelism()
        .map_or(1, |n| n.get())
        .clamp(1, 8);
    let mut groups: Vec<Vec<Vec<String>>> = vec![Vec::new(); workers];
    for (index, name) in names.iter().enumerate() {
        groups[index % workers].extend(
            rows.iter()
                .filter(|row| row.len() > 1 && row[1] == *name)
                .cloned(),
        );
    }
    let mut outcomes: BTreeMap<&'static str, usize> = BTreeMap::new();
    std::thread::scope(|scope| {
        let handles: Vec<_> = groups
            .iter()
            .map(|group| scope.spawn(move || replay_stage_fixture(group)))
            .collect();
        for handle in handles {
            let part = handle
                .join()
                .unwrap_or_else(|panic| std::panic::resume_unwind(panic));
            for (kind, count) in part {
                *outcomes.entry(kind).or_default() += count;
            }
        }
    });
    let replayed: usize = groups
        .iter()
        .map(|group| group.iter().filter(|row| row[0] == "case").count())
        .sum();
    assert_eq!(replayed, names.len(), "every case is replayed once");
    outcomes
}

/// The bound of one feature of a stage fixture: [`tolerance`] of the
/// configuration, and [`area_tolerance`] for its intensity.
fn nonfinite_tolerance(_case: &str, config: &str, _index: usize) -> f64 {
    tolerance(config, None)
}

fn seeds_are_automatic(options: &[String]) -> bool {
    !options.iter().any(|option| option.starts_with("seeds="))
}

/// The `feature`, `meta` and `hulls` rows of one case against `map`.
fn check_nonfinite_features<'a>(
    name: &str,
    config: &str,
    of: &dyn Fn(&str) -> Vec<&'a [String]>,
    map: &FeatureMap,
) {
    let expected = of("feature");
    assert_eq!(expected.len(), map.len(), "{name}: feature rows");
    for (index, row) in expected.iter().enumerate() {
        let relative = nonfinite_tolerance(name, config, index);
        let feature = &map.features[index];
        let what = |field: &str| format!("{name}[{index}].{field}");
        assert_eq!(row[0].parse::<usize>().unwrap(), index);
        close(feature.rt, f64_hex(&row[1]), relative, &what("rt"));
        close(feature.mz, f64_hex(&row[2]), relative, &what("mz"));
        close(
            f64::from(feature.intensity),
            f64::from(f32_hex(&row[3])),
            area_tolerance(config),
            &what("intensity"),
        );
        assert_eq!(feature.charge.to_string(), row[4], "{}", what("charge"));
        close(
            f64::from(feature.quality),
            f64::from(f32_hex(&row[5])),
            relative,
            &what("quality"),
        );
        assert_eq!(feature.quality_rt.to_bits(), f32_hex(&row[6]).to_bits());
        assert_eq!(feature.quality_mz.to_bits(), f32_hex(&row[7]).to_bits());
        close(
            f64::from(feature.width),
            f64::from(f32_hex(&row[8])),
            relative,
            &what("width"),
        );
        assert_eq!(
            feature.subordinates.len().to_string(),
            row[9],
            "{}",
            what("subordinates")
        );
        assert_eq!(
            feature.convex_hulls.len().to_string(),
            row[10],
            "{}",
            what("hulls")
        );
    }
    let mut metas: BTreeMap<usize, BTreeMap<String, (String, String)>> = BTreeMap::new();
    for row in of("meta") {
        metas
            .entry(row[0].parse().unwrap())
            .or_default()
            .insert(row[1].clone(), (row[2].clone(), row[3].clone()));
    }
    assert_eq!(metas.len(), map.len(), "{name}: features with meta values");
    for (index, keys) in metas {
        let feature = &map.features[index];
        assert_eq!(
            feature.metadata.len(),
            keys.len(),
            "{name}[{index}]: meta keys"
        );
        for (key, (kind, value)) in keys {
            let actual = feature
                .metadata
                .get(&key)
                .unwrap_or_else(|| panic!("{name}[{index}]: missing meta {key}"));
            match (kind.as_str(), actual.data()) {
                ("int", MetaValueData::Integer(got)) => {
                    assert_eq!(got.to_string(), value, "{name}[{index}].{key}");
                }
                ("string", MetaValueData::String(got)) => {
                    assert_eq!(got.as_str(), value, "{name}[{index}].{key}");
                }
                ("double", MetaValueData::Float(got)) => close(
                    *got,
                    f64_hex(&value),
                    nonfinite_tolerance(name, config, index),
                    &format!("{name}[{index}].{key}"),
                ),
                other => panic!("{name}[{index}].{key}: unexpected {other:?} for {kind}"),
            }
        }
    }
    let hulls = of("hulls");
    assert_eq!(hulls.len(), map.len(), "{name}: hull rows");
    for row in hulls {
        let index: usize = row[0].parse().unwrap();
        let feature = &map.features[index];
        let counts: Vec<String> = feature
            .convex_hulls
            .iter()
            .map(|hull| hull.hull_points().len().to_string())
            .collect();
        assert_eq!(counts.join(","), row[1], "{name}[{index}]: hull sizes");
        let mut digest = Fnv::new();
        for hull in &feature.convex_hulls {
            let points = hull.hull_points();
            digest.bytes(&(points.len() as u64).to_le_bytes());
            for point in points {
                digest.bytes(&point.rt.to_bits().to_le_bytes());
                digest.bytes(&point.mz.to_bits().to_le_bytes());
            }
        }
        // Hull points are input coordinates, so they are bit-identical.
        assert_eq!(
            format!("{:016x}", digest.0),
            row[2],
            "{name}[{index}]: hulls"
        );
    }
}

// ---------------------------------------------------------------------------
// The source's sorts against the executed library
// ---------------------------------------------------------------------------

/// The sorts `run` reaches, against the executed Linux x86_64 Release build
/// (`sort_probe.tsv.gz`: `../oracle/ffap-complete-fix1/drivers/sort_probe.cpp`
/// on ibminode06, two runs, identical; 2,272 inputs with ties, signed zeros,
/// infinities and NaN keys of four bit patterns, drawn from a 32-value
/// palette).
///
/// - `spec`, `specda`, `chrom`, `chromda`: `MSSpectrum::sortByPosition` and
///   `MSChromatogram::sortByPosition` of one spectrum or chromatogram, without
///   and with a float data array. For every buffer limit the probe applied
///   (`full`; `0`, where every `operator new(nothrow)` fails; and two
///   partial buffers) the libstdc++ `std::stable_sort` port asks for the
///   executed sequence of buffer sizes and leaves the executed order. With the
///   full buffer, what the port runs, `validate_input` leaves the executed
///   order too, with the data array aligned.
/// - `spectra`, `chroms`: `MSExperiment::sortSpectra(true)` and
///   `sortChromatograms(true)` through `validate_input`.
/// - `features`: `FeatureMap::sortByMZ`, the user-seed sort, through
///   `source_sort_by` with `Feature::MZLess`.
#[test]
fn every_source_sort_matches_the_executed_library() {
    use openms::analysis::feature_finder_picked::algorithm::validate_input;
    use openms::kernel::{ChromatogramPeak, DataArray, MSChromatogram, MSSpectrum, Peak1D};
    use openms::math::source_sort::{
        TemporaryBuffer, source_sort_by, source_stable_sort_permutation,
    };
    let rows = stage_rows("sort_probe.tsv.gz");
    let palette: Vec<f64> = rows[0][1..].iter().map(|bits| f64_hex(bits)).collect();
    assert_eq!(rows[0][0], "palette");
    assert_eq!(palette.len(), 32);
    assert_eq!(rows[1], ["sizes", "16", "16", "8"]);
    // An MS1 spectrum that `isSorted(true)` rejects, so that `validate_input`
    // sorts; its own two peaks do not take part in the sorts under test.
    let unsorted_spectrum = || MSSpectrum {
        rt: 0.0,
        peaks: vec![Peak1D::new(200.0, 1.0), Peak1D::new(100.0, 1.0)],
        ..MSSpectrum::default()
    };
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for row in &rows[2..] {
        assert_eq!(row[0], "case");
        let (id, site, limit, codes, requests) = (&row[1], &row[2], &row[3], &row[4], &row[5]);
        let keys: Vec<f64> = (0..codes.len() / 2)
            .map(|i| palette[usize::from_str_radix(&codes[2 * i..2 * i + 2], 16).unwrap()])
            .collect();
        let n = keys.len();
        let executed: Vec<usize> = if row[6].is_empty() {
            Vec::new()
        } else {
            row[6].split(',').map(|v| v.parse().unwrap()).collect()
        };
        let label = format!("{id} {site} {limit}");
        let tag = |i: usize| i as f32;
        let order_of = |intensities: &mut dyn Iterator<Item = f32>| -> Vec<usize> {
            intensities.map(|v| v as usize).collect()
        };
        match site.as_str() {
            "spec" | "specda" | "chrom" | "chromda" => {
                assert!(row.len() == 7 || row[7] == "ok", "{label}");
                let size = if site.ends_with("da") { 8 } else { 16 };
                // `isSorted`: no key less than its predecessor.
                let sorted = !keys.windows(2).any(|pair| pair[1] < pair[0]);
                let mut asked = Vec::new();
                let permutation = if sorted {
                    (0..n).collect()
                } else {
                    let cap = limit.parse::<usize>().ok();
                    let mut grant = |count: usize| {
                        let ok = cap.is_none_or(|cap| count * size <= cap);
                        asked.push(format!("{}:{}", count * size, u8::from(ok)));
                        ok
                    };
                    source_stable_sort_permutation(
                        n,
                        |a, b| keys[a] < keys[b],
                        TemporaryBuffer::Model(&mut grant),
                    )
                    .unwrap()
                };
                let asked = if asked.is_empty() {
                    "-".to_owned()
                } else {
                    asked.join(",")
                };
                assert_eq!(&asked, requests, "{label}: buffer requests");
                assert_eq!(permutation, executed, "{label}");
                if limit == "full" {
                    let array = if site.ends_with("da") {
                        vec![DataArray::new("tag", (0..n).map(tag).collect())]
                    } else {
                        Vec::new()
                    };
                    let mut experiment = MSExperiment::default();
                    if site.starts_with("spec") {
                        experiment.spectra.push(MSSpectrum {
                            rt: 0.0,
                            peaks: keys
                                .iter()
                                .enumerate()
                                .map(|(i, &k)| Peak1D::new(k, tag(i)))
                                .collect(),
                            float_data_arrays: array,
                            ..MSSpectrum::default()
                        });
                        let _ = validate_input(&mut experiment, &mut Vec::new());
                        let spectrum = &experiment.spectra[0];
                        let got = order_of(&mut spectrum.peaks.iter().map(|p| p.intensity));
                        assert_eq!(got, executed, "{label}: validate_input");
                        if let Some(array) = spectrum.float_data_arrays.first() {
                            assert_eq!(
                                order_of(&mut array.data.iter().copied()),
                                executed,
                                "{label}"
                            );
                        }
                    } else {
                        experiment.spectra.push(unsorted_spectrum());
                        experiment.chromatograms.push(MSChromatogram {
                            peaks: keys
                                .iter()
                                .enumerate()
                                .map(|(i, &k)| ChromatogramPeak::new(k, tag(i)))
                                .collect(),
                            float_data_arrays: array,
                            ..MSChromatogram::default()
                        });
                        validate_input(&mut experiment, &mut Vec::new()).unwrap();
                        let chromatogram = &experiment.chromatograms[0];
                        let got = order_of(&mut chromatogram.peaks.iter().map(|p| p.intensity));
                        assert_eq!(got, executed, "{label}: validate_input");
                        if let Some(array) = chromatogram.float_data_arrays.first() {
                            assert_eq!(
                                order_of(&mut array.data.iter().copied()),
                                executed,
                                "{label}"
                            );
                        }
                    }
                }
            }
            "spectra" | "chroms" => {
                assert_eq!(requests, "-", "{label}");
                let mut experiment = MSExperiment::default();
                experiment.spectra.push(unsorted_spectrum());
                if site == "spectra" {
                    // The probe's spectra hold no peaks; the first carries the
                    // two unsorted ones, which no RT comparison sees.
                    experiment.spectra.clear();
                    for (i, &k) in keys.iter().enumerate() {
                        let mut spectrum = if i == 0 {
                            unsorted_spectrum()
                        } else {
                            MSSpectrum::default()
                        };
                        spectrum.rt = k;
                        spectrum.native_id = i.to_string();
                        experiment.spectra.push(spectrum);
                    }
                } else {
                    for (i, &k) in keys.iter().enumerate() {
                        let mut chromatogram = MSChromatogram {
                            native_id: i.to_string(),
                            ..MSChromatogram::default()
                        };
                        chromatogram.product.mz = k;
                        experiment.chromatograms.push(chromatogram);
                    }
                }
                let _ = validate_input(&mut experiment, &mut Vec::new());
                let got: Vec<usize> = if site == "spectra" {
                    experiment
                        .spectra
                        .iter()
                        .map(|s| s.native_id.parse().unwrap())
                        .collect()
                } else {
                    experiment
                        .chromatograms
                        .iter()
                        .map(|c| c.native_id.parse().unwrap())
                        .collect()
                };
                assert_eq!(got, executed, "{label}");
            }
            "features" => {
                assert_eq!(requests, "-", "{label}");
                let mut features: Vec<Feature> = keys
                    .iter()
                    .enumerate()
                    .map(|(i, &k)| Feature::new(0.0, k, tag(i)))
                    .collect();
                source_sort_by(&mut features, |a, b| a.mz < b.mz).unwrap();
                let got = order_of(&mut features.iter().map(|f| f.intensity));
                assert_eq!(got, executed, "{label}");
            }
            other => panic!("unknown site {other}"),
        }
        *counts
            .entry(format!(
                "{site} {}",
                if limit == "full" { "full" } else { "limited" }
            ))
            .or_default() += 1;
    }
    assert_eq!(counts.values().sum::<usize>(), 2272);
}
