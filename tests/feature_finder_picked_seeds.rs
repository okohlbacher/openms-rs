// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The seed stage of `FeatureFinderAlgorithmPicked`: parameters, input
//! validation, intensity/trace/isotope scores, isotope windows and seeds.
//!
//! Evidence (see `docs/FEATURE_FINDER_PICKED_SUPPORT.md` and
//! `tests/data/feature_finder_picked_provenance.json`):
//!
//! - tier 1, executed C++: the `-write_ini` algorithm section of
//!   FeatureFinderCentroided (C1, product SDK), and, from the Linux x86_64
//!   Release build `openms4-release-bc9cc12-c19e494-174b576` (the reference
//!   platform), the library state of the C2 driver `ffap_stages` (effective
//!   members, intensity quantiles, isotope windows, per-peak score arrays,
//!   printed seed counts) and of the B6 driver `seed_stage` (default
//!   parameters, 7 bins and charges 1 to 3, min_spectra 1), re-extracted by
//!   `../oracle/ffap-sem-completion/extract` with the B6 extraction unchanged;
//!   the degenerate intensity bins of the same build (driver
//!   `degenerate_stage`, the generalised `seed_stage`);
//! - adapted: the ordered seed lists, which both drivers re-derive from the
//!   library's score arrays with the source's selection code;
//! - tier 3 and 4: source-review and hand-derived cases for validation, the
//!   conversions and the scoring functions.
//!
//! Floats are compared bit for bit.

#![cfg(all(feature = "mzml", feature = "paramxml"))]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use openms::analysis::feature_finder_picked::algorithm::{
    AbundanceOverride, DegenerateBinStep, HANDLER_NAME, Limits, Options, ReportedMz, RtShape,
    Settings, UNSORTED_WARNING, default_parameters, run, run_with_options, validate_input,
};
use openms::analysis::feature_finder_picked::helper_structs::{
    IsotopePattern, PatternPeak, TheoreticalIsotopePattern,
};
use openms::analysis::feature_finder_picked::scoring::{
    IntensityThresholds, QUANTILE_COUNT, find_isotope, isotope_score, nearest_from, position_score,
};
use openms::analysis::feature_finder_picked::seeds::{SeedStage, overall_score};
use openms::format::{FileHandler, FileType, PeakFileOptions, paramxml};
use openms::kernel::{FeatureMap, NumericRange};
use openms::param::{Param, ParamValue};
use openms::{Error, MSExperiment, MSSpectrum, Peak1D};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/feature_finder_picked")
        .join(name)
}

/// The FeatureFinderCentroided_1 input, reused read-only from the A3 fixtures
/// (sha256 a3dfae63..., test-data 0cb15f2).
fn ffc1_input_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/mzml_mobility/FeatureFinderCentroided_1_input.mzML")
}

/// FeatureFinderCentroided `main_` loading: MS1 only, intensity range
/// `[0, DBL_MAX)` as executed.
fn ffc1_input() -> MSExperiment {
    let mut options = PeakFileOptions::default();
    options.add_ms_level(1).unwrap();
    options.set_intensity_range(NumericRange {
        min: 0.0,
        max: f64::MAX,
    });
    FileHandler::load_experiment_with_options(ffc1_input_path(), &[FileType::MzMl], &options)
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

/// `overall_rounding.tsv` rows of one score table, keyed by spectrum, peak and
/// charge index: (executed score, correctly rounded score). The extraction
/// checks each pair against 60-digit decimal arithmetic and that they are one
/// binary32 step apart; this re-checks the step.
fn rounding_rows(table: &str) -> BTreeMap<(usize, usize, usize), (f32, f32)> {
    rows("overall_rounding.tsv")
        .into_iter()
        .filter(|row| row[0] == table)
        .map(|row| {
            let misrounded = f32_hex(&row[4]);
            let correct = f32_hex(&row[5]);
            assert_eq!(misrounded.to_bits().abs_diff(correct.to_bits()), 1);
            (
                (
                    row[1].parse().unwrap(),
                    row[2].parse().unwrap(),
                    row[3].parse().unwrap(),
                ),
                (misrounded, correct),
            )
        })
        .collect()
}

struct Case {
    records: &'static str,
    config: &'static str,
    scores: &'static str,
    /// First score column: C2 tables carry m/z and intensity first.
    first_array: usize,
}

fn check_settings(settings: &Settings, case: &Case) {
    let members: BTreeMap<String, String> = records(case.records, "member", case.config)
        .into_iter()
        .map(|row| (row[0].clone(), row[1].clone()))
        .collect();
    let float = |name: &str| f64_hex(&members[name]).to_bits();
    assert_eq!(
        settings.pattern_tolerance.to_bits(),
        float("pattern_tolerance")
    );
    assert_eq!(settings.trace_tolerance.to_bits(), float("trace_tolerance"));
    assert_eq!(settings.min_spectra.to_string(), members["min_spectra"]);
    assert_eq!(
        settings.max_missing_trace_peaks.to_string(),
        members["max_missing_trace_peaks"]
    );
    assert_eq!(settings.slope_bound.to_bits(), float("slope_bound"));
    assert_eq!(
        settings.intensity_percentage.to_bits(),
        float("intensity_percentage")
    );
    assert_eq!(
        settings.intensity_percentage_optional.to_bits(),
        float("intensity_percentage_optional")
    );
    assert_eq!(
        settings.optional_fit_improvement.to_bits(),
        float("optional_fit_improvement")
    );
    assert_eq!(
        settings.mass_window_width.to_bits(),
        float("mass_window_width")
    );
    assert_eq!(
        settings.intensity_bins.to_string(),
        members["intensity_bins"]
    );
    assert_eq!(settings.min_isotope_fit.to_bits(), float("min_isotope_fit"));
    assert_eq!(settings.min_trace_score.to_bits(), float("min_trace_score"));
    assert_eq!(settings.min_rt_span.to_bits(), float("min_rt_span"));
    assert_eq!(settings.max_rt_span.to_bits(), float("max_rt_span"));
    assert_eq!(
        settings.max_feature_intersection.to_bits(),
        float("max_feature_intersection")
    );
    let reported = match settings.reported_mz {
        ReportedMz::Maximum => "maximum",
        ReportedMz::Average => "average",
        ReportedMz::Monoisotopic => "monoisotopic",
    };
    assert_eq!(reported, members["reported_mz"]);
}

/// Compare a seed stage with the executed C++ state of one configuration.
fn check_stage(stage: &SeedStage, case: &Case) {
    check_settings(stage.settings(), case);
    let experiment = stage.experiment();

    let input = &records(case.records, "input", case.config)[0];
    assert_eq!(experiment.spectra.len().to_string(), input[0]);
    let peaks: usize = experiment.spectra.iter().map(|s| s.peaks.len()).sum();
    assert_eq!(peaks.to_string(), input[1]);

    let thresholds = stage.thresholds();
    let bins = &records(case.records, "bins", case.config)[0];
    assert_eq!(thresholds.bins().to_string(), bins[0]);
    assert_eq!(thresholds.rt_start().to_bits(), f64_hex(&bins[1]).to_bits());
    assert_eq!(thresholds.mz_start().to_bits(), f64_hex(&bins[2]).to_bits());
    assert_eq!(thresholds.rt_step().to_bits(), f64_hex(&bins[3]).to_bits());
    assert_eq!(thresholds.mz_step().to_bits(), f64_hex(&bins[4]).to_bits());
    let quantiles = records(case.records, "quantiles", case.config);
    assert_eq!(quantiles.len(), thresholds.bins() * thresholds.bins());
    for row in &quantiles {
        let rt: usize = row[0].parse().unwrap();
        let mz: usize = row[1].parse().unwrap();
        let actual = thresholds.quantiles(rt, mz).unwrap();
        for (i, expected) in row[2..].iter().enumerate() {
            assert_eq!(
                actual[i].to_bits(),
                f64_hex(expected).to_bits(),
                "{}: quantile {i} of bin {rt}/{mz}",
                case.config
            );
        }
    }

    let windows = records(case.records, "window", case.config);
    let patterns = stage.windows().patterns();
    assert_eq!(
        patterns.len(),
        windows.len(),
        "{}: window count",
        case.config
    );
    for (pattern, row) in patterns.iter().zip(&windows) {
        let label = format!("{}: window {}", case.config, row[0]);
        assert_eq!(pattern.len().to_string(), row[2], "{label}");
        assert_eq!(pattern.optional_begin.to_string(), row[3], "{label}");
        assert_eq!(pattern.optional_end.to_string(), row[4], "{label}");
        assert_eq!(pattern.max.to_bits(), f64_hex(&row[5]).to_bits(), "{label}");
        assert_eq!(pattern.trimmed_left.to_string(), row[6], "{label}");
        let bits: Vec<u64> = pattern.intensity.iter().map(|v| v.to_bits()).collect();
        let expected: Vec<u64> = row[7..].iter().map(|v| f64_hex(v).to_bits()).collect();
        assert_eq!(bits, expected, "{label}");
    }

    let scores = stage.scores();
    let charges = scores.charge_count();
    let table = rows(case.scores);
    assert_eq!(table.len(), peaks, "{}: score rows", case.config);
    // Overall scores that the executed glibc powf rounds one binary32 step away
    // from the correctly rounded power, each with the correctly rounded value
    // that the port produces (8 of 30,840; the macOS arm64 product SDK's Apple
    // powf misrounded 99).
    let rounding = rounding_rows(case.scores);
    let mut rounded = 0;
    let mut mismatches = Vec::new();
    for row in &table {
        let s: usize = row[0].parse().unwrap();
        let p: usize = row[1].parse().unwrap();
        let arrays = &row[case.first_array..];
        assert_eq!(
            arrays.len(),
            3 + 2 * charges,
            "{}: array count",
            case.config
        );
        let mut actual = vec![
            scores.trace(s).unwrap()[p],
            scores.intensity(s).unwrap()[p],
            scores.local_max(s).unwrap()[p],
        ];
        for c in 0..charges {
            actual.push(scores.pattern(c, s).unwrap()[p]);
        }
        for c in 0..charges {
            actual.push(scores.overall(c, s).unwrap()[p]);
        }
        for (column, (value, oracle)) in actual.iter().zip(arrays).enumerate() {
            let mut expected = f32_hex(oracle);
            if column >= 3 + charges {
                if let Some(&(misrounded, correct)) = rounding.get(&(s, p, column - 3 - charges)) {
                    assert_eq!(misrounded.to_bits(), expected.to_bits(), "{}", case.config);
                    expected = correct;
                    rounded += 1;
                }
            }
            if value.to_bits() != expected.to_bits() {
                mismatches.push((s, p, column, *value, expected));
            }
        }
    }
    assert_eq!(rounded, rounding.len(), "{}: rounding rows", case.config);
    assert!(
        mismatches.is_empty(),
        "{}: {} score mismatches, first {:?}",
        case.config,
        mismatches.len(),
        &mismatches[..mismatches.len().min(5)]
    );

    let stdout: Vec<String> = records(case.records, "stdout", case.config)
        .into_iter()
        .map(|row| row[0].clone())
        .filter(|line| line.contains(" seeds for charge "))
        .collect();
    let logged: Vec<&String> = stage
        .log()
        .iter()
        .filter(|line| line.contains(" seeds for charge "))
        .collect();
    assert_eq!(logged, stdout.iter().collect::<Vec<_>>(), "{}", case.config);
    let seeds = records(case.records, "seed", case.config);
    for charge in stage.charges() {
        let expected: Vec<(String, String, String, String)> = seeds
            .iter()
            .filter(|row| row[0] == charge.charge.to_string())
            .map(|row| {
                (
                    row[2].clone(),
                    row[3].clone(),
                    row[4].clone(),
                    row[5].clone(),
                )
            })
            .collect();
        let charge_index = (charge.charge - stage.settings().charge_low) as usize;
        let actual: Vec<(String, String, String, String)> = charge
            .seeds
            .iter()
            .map(|seed| {
                (
                    seed.spectrum.to_string(),
                    seed.peak.to_string(),
                    format!("{:08x}", seed.intensity.to_bits()),
                    // The drivers record the executed (possibly misrounded) score.
                    format!(
                        "{:08x}",
                        rounding
                            .get(&(seed.spectrum, seed.peak, charge_index))
                            .map_or(
                                scores.overall(charge_index, seed.spectrum).unwrap()[seed.peak],
                                |&(misrounded, _)| misrounded
                            )
                            .to_bits()
                    ),
                )
            })
            .collect();
        assert_eq!(
            actual, expected,
            "{}: seeds of charge {}",
            case.config, charge.charge
        );
    }
    for row in records(case.records, "rt", case.config) {
        let s: usize = row[0].parse().unwrap();
        assert_eq!(
            experiment.spectra[s].rt.to_bits(),
            f64_hex(&row[1]).to_bits()
        );
        assert_eq!(experiment.spectra[s].native_id, row[2]);
    }
}

fn stage(experiment: MSExperiment, seeds: &FeatureMap, parameters: &Param) -> SeedStage {
    SeedStage::run(experiment, seeds, parameters)
        .unwrap()
        .unwrap()
}

// ---------------------------------------------------------------------------
// FeatureFinderAlgorithmPicked_test.cpp sections and the parameter contract
// ---------------------------------------------------------------------------

/// START_SECTION((FeatureFinderAlgorithmPicked())) and the destructor section:
/// the handler constructs with its defaults; the settings and the stage own no
/// external resource, so the destructor section has nothing to release.
#[test]
fn constructor_and_destructor() {
    let defaults = default_parameters().unwrap();
    assert_eq!(defaults.size(), 29);
    let (settings, warnings) = Settings::from_parameters(&Param::new()).unwrap();
    assert!(warnings.is_empty());
    assert_eq!(settings.intensity_bins, 10);
}

/// The defaults equal the algorithm section that the executed
/// FeatureFinderCentroided writes with `-write_ini` (C1), entry by entry and as
/// a whole, including descriptions, restrictions, tags and section descriptions.
#[test]
fn default_parameters_equal_the_executed_write_ini_section() {
    let written = paramxml::load(data("FeatureFinderCentroided_defaults.ini"))
        .unwrap()
        .copy("FeatureFinderCentroided:1:algorithm:", true)
        .unwrap();
    let defaults = default_parameters().unwrap();
    let written_items: Vec<_> = written.iter().unwrap().collect();
    let default_items: Vec<_> = defaults.iter().unwrap().collect();
    assert_eq!(written_items.len(), 29);
    assert_eq!(default_items.len(), 29);
    for (w, d) in written_items.iter().zip(&default_items) {
        assert_eq!(w.key, d.key);
        assert_eq!(w.entry, d.entry, "{}", d.key);
    }
    for section in [
        "intensity",
        "mass_trace",
        "isotopic_pattern",
        "seed",
        "fit",
        "feature",
        "user-seed",
        "advanced",
    ] {
        assert_eq!(
            written.section_description(section).unwrap(),
            defaults.section_description(section).unwrap(),
            "{section}"
        );
    }
    assert_eq!(written, defaults);
}

/// The 27 algorithm items of `FeatureFinderCentroided_1_parameters.ini` are
/// known defaults of the same value type; only the two abundances are missing,
/// so the strict tool path accepts the file without an unknown-parameter warning.
#[test]
fn ffc_1_parameters_are_a_subset_of_the_defaults() {
    let parameters = ffc1_parameters();
    let defaults = default_parameters().unwrap();
    let keys: Vec<String> = parameters.iter().unwrap().map(|item| item.key).collect();
    assert_eq!(keys.len(), 27);
    for item in parameters.iter().unwrap() {
        let default = defaults.entry(&item.key).unwrap();
        assert_eq!(
            default.value.value_type(),
            item.entry.value.value_type(),
            "{}",
            item.key
        );
    }
    let missing: Vec<String> = defaults
        .iter()
        .unwrap()
        .map(|item| item.key)
        .filter(|key| !keys.contains(key))
        .collect();
    assert_eq!(
        missing,
        [
            "isotopic_pattern:abundance_12C",
            "isotopic_pattern:abundance_14N"
        ]
    );
    let (_, warnings) = Settings::from_parameters(&parameters).unwrap();
    assert!(warnings.is_empty(), "{warnings:?}");
}

/// The class-test INI carries a legacy `debug` item; the source warns about it
/// and runs.
#[test]
fn unknown_parameters_are_warnings() {
    let (_, warnings) = Settings::from_parameters(&class_test_parameters()).unwrap();
    assert_eq!(
        warnings,
        [format!("{HANDLER_NAME}: unknown parameter 'debug'")]
    );
}

/// `updateMembers_` and the `run_` reads, beyond the executed configurations.
#[test]
fn settings_follow_update_members() {
    let (defaults, _) = Settings::from_parameters(&Param::new()).unwrap();
    assert_eq!(defaults.min_spectra, 5);
    assert_eq!(defaults.charge_low, 1);
    assert_eq!(defaults.charge_high, 4);
    assert_eq!(defaults.charge_count().unwrap(), 4);
    assert_eq!(defaults.max_iterations, 500);
    assert_eq!(defaults.max_isotopes(), 20);
    assert_eq!(defaults.rt_shape, RtShape::Symmetric);
    assert_eq!(defaults.reported_mz, ReportedMz::Monoisotopic);
    assert!(!defaults.write_debug);
    assert!(!defaults.abundance_12c_changed && !defaults.abundance_14n_changed);
    assert_eq!(defaults.seed_min_score, 0.8);
    assert_eq!(defaults.min_feature_score, 0.7);
    assert_eq!(defaults.user_seed_rt_tolerance, 5.0);
    assert_eq!(defaults.user_seed_mz_tolerance, 1.1);
    assert_eq!(defaults.user_seed_min_score, 0.5);
    assert_eq!(defaults.intensity_percentage, 10.0 / 100.0);
    assert_eq!(defaults.intensity_percentage_optional, 0.1 / 100.0);
    assert_eq!(defaults.optional_fit_improvement, 2.0 / 100.0);

    // min_spectra_ = floor(mass_trace:min_spectra * 0.5)
    for (given, half) in [(1, 0), (2, 1), (3, 1), (14, 7), (15, 7)] {
        let mut p = Param::new();
        set(&mut p, "mass_trace:min_spectra", ParamValue::Integer(given));
        assert_eq!(Settings::from_parameters(&p).unwrap().0.min_spectra, half);
    }

    // Each changed abundance adds 1000 isotopes; equality is exact.
    let mut p = Param::new();
    set(
        &mut p,
        "isotopic_pattern:abundance_12C",
        ParamValue::Float(98.93),
    );
    assert!(
        !Settings::from_parameters(&p)
            .unwrap()
            .0
            .abundance_12c_changed
    );
    set(
        &mut p,
        "isotopic_pattern:abundance_12C",
        ParamValue::Float(90.0),
    );
    let (one, _) = Settings::from_parameters(&p).unwrap();
    assert!(one.abundance_12c_changed);
    assert_eq!(one.max_isotopes(), 1020);
    set(
        &mut p,
        "isotopic_pattern:abundance_14N",
        ParamValue::Float(99.0),
    );
    assert_eq!(
        Settings::from_parameters(&p).unwrap().0.max_isotopes(),
        2020
    );

    let mut p = Param::new();
    set(
        &mut p,
        "feature:rt_shape",
        ParamValue::String("asymmetric".into()),
    );
    set(
        &mut p,
        "feature:reported_mz",
        ParamValue::String("average".into()),
    );
    set(&mut p, "write_debug", ParamValue::String("true".into()));
    set(&mut p, "fit:max_iterations", ParamValue::Integer(7));
    let (s, _) = Settings::from_parameters(&p).unwrap();
    assert_eq!(s.rt_shape, RtShape::Asymmetric);
    assert_eq!(s.reported_mz, ReportedMz::Average);
    assert!(s.write_debug);
    assert_eq!(s.max_iterations, 7);

    // charge_low == charge_high + 1 is zero charges; more is refused.
    let mut p = Param::new();
    set(
        &mut p,
        "isotopic_pattern:charge_low",
        ParamValue::Integer(5),
    );
    let (s, _) = Settings::from_parameters(&p).unwrap();
    assert_eq!(s.charge_count().unwrap(), 0);
    set(
        &mut p,
        "isotopic_pattern:charge_low",
        ParamValue::Integer(6),
    );
    let (s, _) = Settings::from_parameters(&p).unwrap();
    assert!(matches!(s.charge_count(), Err(Error::InvalidValue(_))));
}

/// Restriction and type violations are errors (source
/// `Exception::InvalidParameter`).
#[test]
fn invalid_parameters_are_errors() {
    for (key, value) in [
        ("intensity:bins", ParamValue::Integer(0)),
        ("mass_trace:min_spectra", ParamValue::Integer(0)),
        (
            "isotopic_pattern:intensity_percentage",
            ParamValue::Float(100.5),
        ),
        ("isotopic_pattern:mass_window_width", ParamValue::Float(0.5)),
        ("feature:rt_shape", ParamValue::String("skewed".into())),
        ("write_debug", ParamValue::String("yes".into())),
        ("intensity:bins", ParamValue::Float(10.0)),
    ] {
        let mut p = Param::new();
        set(&mut p, key, value);
        assert!(
            matches!(Settings::from_parameters(&p), Err(Error::InvalidValue(_))),
            "{key}"
        );
    }
}

// ---------------------------------------------------------------------------
// Executed C++ state, steps 0 to 3.2
// ---------------------------------------------------------------------------

#[test]
fn ffc_1_stage_matches_the_executed_library() {
    let s = stage(ffc1_input(), &FeatureMap::new(), &ffc1_parameters());
    check_stage(
        &s,
        &Case {
            records: "c2_stage_records.tsv",
            config: "ffc1_symmetric",
            scores: "scores_ffc1.tsv",
            first_array: 4,
        },
    );
    // The score table also pins the loaded input peak for peak.
    for row in rows("scores_ffc1.tsv") {
        let s_index: usize = row[0].parse().unwrap();
        let p: usize = row[1].parse().unwrap();
        let peak = s.experiment().spectra[s_index].peaks[p];
        assert_eq!(peak.mz.to_bits(), f64_hex(&row[2]).to_bits());
        assert_eq!(peak.intensity.to_bits(), f32_hex(&row[3]).to_bits());
    }
    assert_eq!(s.charges().len(), 1);
    assert_eq!(s.charges()[0].seeds.len(), 25);
    assert_eq!(s.log(), ["Found 25 seeds for charge 2."]);
}

#[test]
fn class_test_stage_matches_the_executed_library() {
    let s = stage(
        class_test_input(),
        &FeatureMap::new(),
        &class_test_parameters(),
    );
    check_stage(
        &s,
        &Case {
            records: "c2_stage_records.tsv",
            config: "classtest",
            scores: "scores_ffc1.tsv",
            first_array: 4,
        },
    );
    assert_eq!(
        s.log(),
        [
            format!("{HANDLER_NAME}: unknown parameter 'debug'"),
            "Found 25 seeds for charge 2.".to_string()
        ]
    );
}

/// [EXTRA] #9247: the isotope-pattern and mass-trace tolerances feed different
/// scores; both asymmetric configurations reproduce the executed arrays.
#[test]
fn tolerance_swap_stages_match_the_executed_library() {
    for (config, scores, pattern, trace, count) in [
        (
            "classtest_9247_tight_pattern",
            "scores_tight_pattern.tsv",
            0.005,
            0.5,
            15,
        ),
        (
            "classtest_9247_tight_trace",
            "scores_tight_trace.tsv",
            0.5,
            0.005,
            18,
        ),
    ] {
        let mut parameters = class_test_parameters();
        set(
            &mut parameters,
            "isotopic_pattern:mz_tolerance",
            ParamValue::Float(pattern),
        );
        set(
            &mut parameters,
            "mass_trace:mz_tolerance",
            ParamValue::Float(trace),
        );
        let s = stage(class_test_input(), &FeatureMap::new(), &parameters);
        check_stage(
            &s,
            &Case {
                records: "c2_stage_records.tsv",
                config,
                scores,
                first_array: 4,
            },
        );
        assert_eq!(s.charges()[0].seeds.len(), count);
    }
}

/// User seeds from the retained FFC_1 output (reused read-only from the A3
/// fixtures, sha256 4eaa53c8...): 24 seeds near the 8 user seeds.
#[cfg(feature = "featurexml")]
#[test]
fn user_seed_stage_matches_the_executed_library() {
    let reader = std::io::BufReader::new(
        std::fs::File::open(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/data/mzml_mobility/FeatureFinderCentroided_1_1_output.featureXML"),
        )
        .unwrap(),
    );
    let user_seeds = FileHandler::read_feature_map(reader, FileType::FeatureXml).unwrap();
    assert_eq!(user_seeds.len(), 8);
    let s = stage(ffc1_input(), &user_seeds, &ffc1_parameters());
    assert_eq!(s.user_seeds().len(), 8);
    assert!(s.user_seeds().windows(2).all(|w| w[0].mz <= w[1].mz));
    check_stage(
        &s,
        &Case {
            records: "c2_stage_records.tsv",
            config: "ffc1_user_seeds",
            scores: "scores_ffc1.tsv",
            first_array: 4,
        },
    );
    assert_eq!(s.charges()[0].seeds.len(), 24);
}

/// Default parameters: 10 intensity bins (the four-bin interpolation), charges 1
/// to 4, min_spectra 5, window width 25 (B6 driver).
#[test]
fn default_parameter_stage_matches_the_executed_library() {
    let s = stage(ffc1_input(), &FeatureMap::new(), &Param::new());
    check_stage(
        &s,
        &Case {
            records: "b6_stage_records.tsv",
            config: "ffc1_defaults",
            scores: "scores_defaults.tsv",
            first_array: 2,
        },
    );
    let counts: Vec<usize> = s.charges().iter().map(|c| c.seeds.len()).collect();
    assert_eq!(counts, [42, 44, 0, 0]);
}

/// An odd bin count and three charges on the FFC_1 parameters (B6 driver).
#[test]
fn seven_bin_three_charge_stage_matches_the_executed_library() {
    let mut parameters = ffc1_parameters();
    set(&mut parameters, "intensity:bins", ParamValue::Integer(7));
    set(
        &mut parameters,
        "isotopic_pattern:charge_low",
        ParamValue::Integer(1),
    );
    set(
        &mut parameters,
        "isotopic_pattern:charge_high",
        ParamValue::Integer(3),
    );
    let s = stage(ffc1_input(), &FeatureMap::new(), &parameters);
    check_stage(
        &s,
        &Case {
            records: "b6_stage_records.tsv",
            config: "ffc1_bins7_charge13",
            scores: "scores_bins7_charge13.tsv",
            first_array: 2,
        },
    );
    let counts: Vec<usize> = s.charges().iter().map(|c| c.seeds.len()).collect();
    assert_eq!(counts, [16, 23, 0]);
}

/// The inclusive cell traversal equals the kernel's area iterator
/// (`areaBeginConst`) for every cell of a seven-bin grid.
#[test]
fn intensity_bins_equal_the_area_iterator() {
    let experiment = ffc1_input();
    for bins in [1, 7, 10] {
        let thresholds = IntensityThresholds::compute(&experiment, bins).unwrap();
        for rt in 0..bins {
            let min_rt = thresholds.rt_start() + rt as f64 * thresholds.rt_step();
            let max_rt = thresholds.rt_start() + (rt + 1) as f64 * thresholds.rt_step();
            for mz in 0..bins {
                let min_mz = thresholds.mz_start() + mz as f64 * thresholds.mz_step();
                let max_mz = thresholds.mz_start() + (mz + 1) as f64 * thresholds.mz_step();
                let mut values: Vec<f64> = experiment
                    .area_begin(min_rt, max_rt, min_mz, max_mz, 1)
                    .unwrap()
                    .map(|point| f64::from(point.peak.intensity))
                    .collect();
                let mut expected = [0.0; QUANTILE_COUNT];
                if !values.is_empty() {
                    values.sort_by(f64::total_cmp);
                    for (i, q) in expected.iter_mut().enumerate() {
                        *q = values[(0.05 * i as f64 * (values.len() - 1) as f64).floor() as usize];
                    }
                }
                assert_eq!(thresholds.quantiles(rt, mz).unwrap(), &expected);
            }
        }
    }
}

/// mass_trace:min_spectra 1 gives min_spectra_ 0. The executed C++ divides every
/// trace score by zero (NaN), finds no seed and exits normally (B6 driver), and
/// the port follows it (lead decision of 2026-09-15, `CPP-271`).
#[test]
fn min_spectra_one_follows_the_source_and_finds_no_seed() {
    assert_eq!(
        records("b6_stage_records.tsv", "stdout", "ffc1_min_spectra_1")[0][0],
        "Found 0 seeds for charge 2."
    );
    assert_eq!(
        records("b6_stage_records.tsv", "member", "ffc1_min_spectra_1")
            .iter()
            .find(|row| row[0] == "min_spectra")
            .unwrap()[1],
        "0"
    );
    let mut parameters = ffc1_parameters();
    set(
        &mut parameters,
        "mass_trace:min_spectra",
        ParamValue::Integer(1),
    );
    let stage = SeedStage::run(ffc1_input(), &FeatureMap::new(), &parameters)
        .unwrap()
        .unwrap();
    assert_eq!(stage.settings().min_spectra, 0);
    assert_eq!(stage.charges().len(), 1);
    assert!(stage.charges()[0].seeds.is_empty());
    assert_eq!(stage.log(), ["Found 0 seeds for charge 2."]);
    // Every trace score is the NaN of 0.0 / 0, so no peak can reach the seed
    // threshold; the whole run therefore ends with an empty map.
    assert!(stage.scores().trace(0).unwrap().iter().all(|s| s.is_nan()));
    let output = run(ffc1_input(), &FeatureMap::new(), &parameters).unwrap();
    assert!(output.features.is_empty());
    // The two empty entries are the blank lines the source prints before the
    // abort block and before the feature count
    // (`FeatureFinderAlgorithmPicked.cpp:1019` and `1026`), which it prints even
    // when there is no abort reason: the C++ Release build
    // `openms4-release-bc9cc12-c19e494-174b576` writes exactly this shape on a
    // run that finds nothing (`FileConverter_31_output.mzML`).
    assert_eq!(
        output.log,
        [
            "Found 0 seeds for charge 2.",
            "Found 0 feature candidates for charge 2.",
            "Removed 0 overlapping features.",
            "",
            "Info: reasons for not finalizing a feature during its construction:",
            "",
            "0 features found.",
        ]
    );
}

// ---------------------------------------------------------------------------
// run() validation (source order) and the unported remainder
// ---------------------------------------------------------------------------

fn spectrum(rt: f64, level: u32, peaks: &[(f64, f32)]) -> MSSpectrum {
    MSSpectrum {
        rt,
        ms_level: level,
        peaks: peaks.iter().map(|&(mz, i)| Peak1D::new(mz, i)).collect(),
        ..MSSpectrum::default()
    }
}

fn experiment(spectra: Vec<MSSpectrum>) -> MSExperiment {
    MSExperiment {
        spectra,
        ..MSExperiment::default()
    }
}

#[test]
fn empty_input_returns_an_empty_map_before_the_parameters_are_read() {
    let mut invalid = Param::new();
    set(&mut invalid, "intensity:bins", ParamValue::Integer(0));
    let output = run(MSExperiment::new(), &FeatureMap::new(), &invalid).unwrap();
    assert!(output.features.is_empty());
    assert!(output.log.is_empty());
    assert!(
        SeedStage::run(MSExperiment::new(), &FeatureMap::new(), &invalid)
            .unwrap()
            .is_none()
    );
}

#[test]
fn input_checks_follow_the_source() {
    let message = |result: openms::Result<Option<SeedStage>>| match result {
        Err(Error::InvalidValue(message)) => message,
        other => panic!("expected InvalidValue, got {other:?}"),
    };
    let none = FeatureMap::new();
    let p = Param::new();
    // No peak at all.
    let e = experiment(vec![spectrum(1.0, 1, &[]), spectrum(2.0, 1, &[])]);
    assert!(message(SeedStage::run(e, &none, &p)).contains("needs updated ranges"));
    // MS2 data, alone or mixed.
    let e = experiment(vec![spectrum(1.0, 2, &[(500.0, 1.0)])]);
    assert!(message(SeedStage::run(e, &none, &p)).contains("MS level 1"));
    let e = experiment(vec![
        spectrum(1.0, 1, &[(500.0, 1.0)]),
        spectrum(2.0, 2, &[(500.0, 1.0)]),
    ]);
    assert!(message(SeedStage::run(e, &none, &p)).contains("MS level 1"));
    // A negative first m/z, checked after sorting.
    let e = experiment(vec![spectrum(1.0, 1, &[(500.0, 1.0), (-1.0, 1.0)])]);
    assert!(message(SeedStage::run(e, &none, &p)).contains("positive m/z"));
    // Non-finite values (native).
    let e = experiment(vec![spectrum(1.0, 1, &[(500.0, f32::NAN)])]);
    assert!(message(SeedStage::run(e, &none, &p)).contains("finite"));
    // The parameters are checked after the input.
    let mut invalid = Param::new();
    set(&mut invalid, "intensity:bins", ParamValue::Integer(0));
    let e = experiment(vec![spectrum(1.0, 2, &[(500.0, 1.0)])]);
    assert!(message(SeedStage::run(e, &none, &invalid)).contains("MS level 1"));
    let e = experiment(vec![
        spectrum(1.0, 1, &[(500.0, 1.0)]),
        spectrum(2.0, 1, &[(501.0, 2.0)]),
    ]);
    assert!(matches!(
        SeedStage::run(e, &none, &invalid),
        Err(Error::InvalidValue(_))
    ));
}

/// Unsorted input is sorted with the source warning and scored like the sorted
/// input.
#[test]
fn unsorted_input_is_sorted_with_a_warning() {
    let sorted = stage(ffc1_input(), &FeatureMap::new(), &ffc1_parameters());
    let mut reversed = ffc1_input();
    reversed.spectra.reverse();
    reversed.spectra[3].peaks.reverse();
    let mut log = Vec::new();
    let mut copy = reversed.clone();
    assert!(validate_input(&mut copy, &mut log).unwrap());
    assert_eq!(log, [UNSORTED_WARNING]);
    let s = stage(reversed, &FeatureMap::new(), &ffc1_parameters());
    assert_eq!(s.log()[0], UNSORTED_WARNING);
    assert_eq!(s.experiment(), sorted.experiment());
    assert_eq!(s.scores(), sorted.scores());
    assert_eq!(s.charges(), sorted.charges());
    assert_eq!(s.windows(), sorted.windows());
    assert_eq!(s.thresholds(), sorted.thresholds());
}

/// Past seed selection the run continues into the feature stage, which
/// `tests/feature_finder_picked.rs` compares with the executed C++ in detail;
/// here only that it runs and reports the FeatureFinderCentroided_1 counts.
#[test]
fn run_continues_past_seed_selection() {
    let output = run(ffc1_input(), &FeatureMap::new(), &ffc1_parameters()).unwrap();
    assert_eq!(output.features.len(), 8);
    assert_eq!(output.log[0], "Found 25 seeds for charge 2.");
    assert_eq!(output.log[1], "Found 8 feature candidates for charge 2.");
}

#[test]
fn undefined_source_configurations_are_refused() {
    let none = FeatureMap::new();
    // write_debug reads an undeclared parameter in the source.
    let mut p = ffc1_parameters();
    set(&mut p, "write_debug", ParamValue::String("true".into()));
    assert!(matches!(
        SeedStage::run(ffc1_input(), &none, &p),
        Err(Error::Unsupported(_))
    ));
    // charge_low more than one above charge_high.
    let mut p = ffc1_parameters();
    set(
        &mut p,
        "isotopic_pattern:charge_low",
        ParamValue::Integer(4),
    );
    assert!(matches!(
        SeedStage::run(ffc1_input(), &none, &p),
        Err(Error::InvalidValue(_))
    ));
    // charge_low one above charge_high: no charge, no seed.
    let mut p = ffc1_parameters();
    set(
        &mut p,
        "isotopic_pattern:charge_low",
        ParamValue::Integer(3),
    );
    let s = stage(ffc1_input(), &none, &p);
    assert!(s.charges().is_empty());
    assert_eq!(s.scores().charge_count(), 0);
    // A single retention time or a single m/z: a zero bin step, for which the
    // source converts floor(NaN) to UInt. The default follows the Linux x86_64
    // Release build (every intensity score NaN, no seed; executed in
    // `degenerate_bin_steps_match_the_linux_release_build`), and
    // DegenerateBinStep::Refuse refuses once the seed loop reads the scores:
    // with the default min_spectra 10 (min_spectra_ 5), from 11 scans on. Two
    // scans never reach the seed loop and are not refused.
    let refuse = Options {
        degenerate_bin_step: DegenerateBinStep::Refuse,
        ..Options::default()
    };
    let one_rt = |n: usize| {
        experiment(
            (0..n)
                .map(|_| spectrum(1.0, 1, &[(500.0, 1.0), (501.0, 2.0)]))
                .collect(),
        )
    };
    let one_mz = |n: usize| {
        experiment(
            (0..n)
                .map(|s| spectrum(s as f64, 1, &[(500.0, 1.0 + s as f32)]))
                .collect(),
        )
    };
    for e in [one_rt(11), one_mz(11)] {
        assert!(matches!(
            SeedStage::run_with_options(e.clone(), &none, &Param::new(), &refuse),
            Err(Error::InvalidValue(_))
        ));
        let s = stage(e, &none, &Param::new());
        assert!(
            s.scores()
                .intensity(5)
                .unwrap()
                .iter()
                .all(|score| score.to_bits() == 0xffc0_0000)
        );
        assert!(s.charges().iter().all(|charge| charge.seeds.is_empty()));
    }
    for e in [one_rt(2), one_mz(2), one_rt(10), one_mz(10)] {
        let refusing = SeedStage::run_with_options(e.clone(), &none, &Param::new(), &refuse)
            .unwrap()
            .unwrap();
        // NaN scores make the stages unequal under `PartialEq`; compare bits.
        let default = stage(e, &none, &Param::new());
        let bits = |s: &SeedStage| -> Vec<u32> {
            (0..s.experiment().spectra.len())
                .flat_map(|index| s.scores().intensity(index).unwrap().to_vec())
                .map(f32::to_bits)
                .collect()
        };
        assert_eq!(bits(&refusing), bits(&default));
        assert!(bits(&refusing).iter().all(|&b| b == 0xffc0_0000));
        assert_eq!(refusing.log(), default.log());
        assert_eq!(refusing.charges(), default.charges());
        assert!(
            refusing
                .charges()
                .iter()
                .all(|charge| charge.seeds.is_empty())
        );
    }
    // A non-finite user seed.
    let mut seeds = FeatureMap::new();
    seeds
        .features
        .push(openms::kernel::Feature::new(100.0, f64::NAN, 1.0));
    assert!(matches!(
        SeedStage::run(ffc1_input(), &seeds, &ffc1_parameters()),
        Err(Error::InvalidValue(_))
    ));
}

/// A changed abundance computes the intended two-isotope override by default
/// (lead decision of 2026-09-15, `CPP-247`); it differs from the executed C++,
/// whose override keeps a stray (0, 1) peak, so that its first window has 27
/// bins and the run finds no seed (C2 `ffap_ffc1_abundance_12C_90`).
/// [`AbundanceOverride::Refuse`] is the opt-in that refuses instead of
/// differing.
#[test]
fn abundance_overrides_use_the_intended_override_unless_refusal_is_selected() {
    let mut p = ffc1_parameters();
    set(
        &mut p,
        "isotopic_pattern:abundance_12C",
        ParamValue::Float(90.0),
    );
    let s = SeedStage::run(ffc1_input(), &FeatureMap::new(), &p)
        .unwrap()
        .unwrap();
    assert_eq!(s.settings().max_isotopes(), 1020);
    let first = &s.windows().patterns()[0];
    assert!(first.len() <= 6, "{}", first.len());
    assert_eq!(first.intensity.iter().copied().fold(0.0, f64::max), 1.0);
    // The executed C++ found no seed here; the intended override does.
    assert!(!s.charges()[0].seeds.is_empty());

    let refusing = Options {
        abundance_override: AbundanceOverride::Refuse,
        ..Options::default()
    };
    assert!(matches!(
        SeedStage::run_with_options(ffc1_input(), &FeatureMap::new(), &p, &refusing),
        Err(Error::Unsupported(_))
    ));
    assert!(matches!(
        run_with_options(ffc1_input(), &FeatureMap::new(), &p, &refusing),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn limits_are_checked_before_the_work() {
    let p = ffc1_parameters();
    let none = FeatureMap::new();
    let limited = |limits: Limits| {
        SeedStage::run_with_options(
            ffc1_input(),
            &none,
            &p,
            &Options {
                limits,
                ..Options::default()
            },
        )
    };
    for limits in [
        Limits {
            max_spectra: 111,
            ..Limits::default()
        },
        Limits {
            max_peaks: 3083,
            ..Limits::default()
        },
        Limits {
            max_charges: 0,
            ..Limits::default()
        },
        Limits {
            max_intensity_bins: 0,
            ..Limits::default()
        },
        Limits {
            max_isotope_windows: 14,
            ..Limits::default()
        },
        Limits {
            max_pattern_values: 15 * 20 - 1,
            ..Limits::default()
        },
        Limits {
            max_score_bytes: 3084 * 5 * 4 - 1,
            ..Limits::default()
        },
        Limits {
            max_work: 1000,
            ..Limits::default()
        },
    ] {
        assert!(
            matches!(limited(limits), Err(Error::InvalidValue(_))),
            "{limits:?}"
        );
    }
    assert!(
        limited(Limits {
            max_spectra: 112,
            max_peaks: 3084,
            max_charges: 1,
            max_intensity_bins: 1,
            max_isotope_windows: 15,
            max_pattern_values: 15 * 20,
            max_score_bytes: 3084 * 5 * 4,
            ..Limits::default()
        })
        .is_ok()
    );
}

// ---------------------------------------------------------------------------
// Scoring functions (source review and hand-derived values)
// ---------------------------------------------------------------------------

#[test]
fn position_score_follows_the_source_formula() {
    assert_eq!(
        position_score(100.0, 100.0, 0.02),
        0.1 * (0.01 - 0.0) / 0.01 + 0.9
    );
    assert_eq!(position_score(100.0, 100.0, 0.02), 1.0);
    let d = (100.0f64 - 100.005).abs();
    assert_eq!(
        position_score(100.0, 100.005, 0.02),
        0.1 * (0.5 * 0.02 - d) / (0.5 * 0.02) + 0.9
    );
    let d = (100.0f64 - 100.015).abs();
    assert_eq!(
        position_score(100.0, 100.015, 0.02),
        0.9 * (0.02 - d) / (0.5 * 0.02)
    );
    assert_eq!(position_score(100.0, 100.0201, 0.02), 0.0);
    // Zero tolerance: 0/0 at equal positions, as in the source.
    assert!(position_score(100.0, 100.0, 0.0).is_nan());
    assert_eq!(position_score(100.0, 100.1, 0.0), 0.0);
}

#[test]
fn nearest_from_walks_upwards_to_the_first_local_minimum() {
    let peaks: Vec<Peak1D> = [1.0, 2.0, 3.0, 4.0, 5.0]
        .iter()
        .map(|&mz| Peak1D::new(mz, 1.0))
        .collect();
    assert_eq!(nearest_from(&peaks, 3.1, 0), Some((2, 2)));
    assert_eq!(nearest_from(&peaks, 3.1, 2), Some((2, 0)));
    // Never walks down, and ties keep the lower index.
    assert_eq!(nearest_from(&peaks, 1.0, 3), Some((3, 0)));
    assert_eq!(nearest_from(&peaks, 2.5, 1), Some((1, 0)));
    assert_eq!(nearest_from(&peaks, 3.0, 5), None);
    // Stops at the first local minimum even when a later peak is closer.
    let unsorted: Vec<Peak1D> = [1.0, 3.0, 4.0, 2.95]
        .iter()
        .map(|&mz| Peak1D::new(mz, 1.0))
        .collect();
    assert_eq!(nearest_from(&unsorted, 2.95, 0), Some((1, 1)));
    assert_eq!(nearest_from(&[], 3.0, 0), None);
}

fn theoretical(
    intensity: &[f64],
    optional_begin: usize,
    optional_end: usize,
) -> TheoreticalIsotopePattern {
    TheoreticalIsotopePattern {
        intensity: intensity.to_vec(),
        optional_begin,
        optional_end,
        max: 1.0,
        trimmed_left: 0,
    }
}

fn found_pattern(intensity: &[f64], mz_score: &[f64]) -> IsotopePattern {
    let mut pattern = IsotopePattern::new(intensity.len()).unwrap();
    for i in 0..intensity.len() {
        pattern.peak[i] = PatternPeak::Found(i);
        pattern.intensity[i] = intensity[i];
        pattern.mz_score[i] = mz_score[i];
    }
    pattern
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    openms::math::statistic_functions::pearson_correlation_coefficient(a, b).unwrap()
}

#[test]
fn isotope_score_requires_the_core_peaks() {
    let isotopes = theoretical(&[0.05, 1.0, 0.6, 0.25, 0.05], 1, 1);
    let mut pattern = found_pattern(&[0.1, 1.0, 0.5, 0.3, 0.1], &[1.0; 5]);
    pattern.peak[2] = PatternPeak::NotFound;
    assert_eq!(
        isotope_score(&isotopes, &mut pattern, true, 0.8, 0.02).unwrap(),
        0.0
    );
    // A missing optional peak is fine.
    let mut pattern = found_pattern(&[0.0, 1.0, 0.5, 0.3, 0.1], &[1.0; 5]);
    pattern.peak[0] = PatternPeak::NotFound;
    let score = isotope_score(&isotopes, &mut pattern, false, 0.8, 0.02).unwrap();
    let expected = pearson(&isotopes.intensity[1..5], &pattern.intensity[1..5]);
    assert_eq!(score, expected);
    // Mismatched lengths are an error, not undefined behaviour.
    let mut short = found_pattern(&[1.0, 1.0], &[1.0, 1.0]);
    assert!(isotope_score(&isotopes, &mut short, false, 0.8, 0.02).is_err());
}

/// The source re-reads the best trailing count at the start of every inner
/// loop: after (b 0, e 1) became the best fit, (b 1, e 0) is never evaluated,
/// although its correlation (0.850) would win. The last peak is removed, the
/// first kept.
#[test]
fn isotope_score_narrows_later_candidates_after_a_new_best_fit() {
    let isotopes = theoretical(&[0.05, 1.0, 0.6, 0.25, 0.05], 1, 1);
    let found = [0.05, 0.93, 0.96, 0.82, 0.79];
    let s00 = pearson(&isotopes.intensity, &found);
    let s01 = pearson(&isotopes.intensity[..4], &found[..4]);
    let s10 = pearson(&isotopes.intensity[1..], &found[1..]);
    let s11 = pearson(&isotopes.intensity[1..4], &found[1..4]);
    assert!(s00 / 0.01 >= 1.02 && s01 / s00 >= 1.02 && s10 / s01 >= 1.02);
    assert!(s11 / s01 < 1.02);
    let mut pattern = found_pattern(&found, &[1.0; 5]);
    let score = isotope_score(&isotopes, &mut pattern, true, 0.8, 0.02).unwrap();
    assert_eq!(score, s01 * (4.0 / 4.0));
    assert_eq!(pattern.peak[0], PatternPeak::Found(0));
    assert_eq!(pattern.peak[4], PatternPeak::Removed);
    assert_eq!(pattern.intensity[4], 0.0);
    assert_eq!(pattern.mz_score[4], 0.0);
}

/// Two remaining isotopes are only tried first and are capped at
/// min_isotope_fit; the m/z factor is the mean m/z score of the kept isotopes.
#[test]
fn isotope_score_caps_two_isotope_fits_and_weights_by_mz() {
    let isotopes = theoretical(&[1.0, 0.5, 0.05], 0, 1);
    let mut pattern = found_pattern(&[1.0, 0.4, 0.2], &[0.9, 0.8, 0.7]);
    pattern.peak[2] = PatternPeak::NotFound;
    // best_end starts at 1 (missing optional peak): only (0, 1) with two
    // isotopes is evaluated; its correlation 1 is capped at 0.8.
    let score = isotope_score(&isotopes, &mut pattern, true, 0.8, 0.02).unwrap();
    assert_eq!(score, 0.8 * ((0.9 + 0.8) / 2.0));
    assert_eq!(pattern.peak[2], PatternPeak::Removed);
    // No candidate: the starting best score 0.01 is returned.
    let isotopes = theoretical(&[1.0, 0.5], 0, 0);
    let mut pattern = found_pattern(&[1.0, 1.0], &[1.0, 1.0]);
    assert_eq!(
        isotope_score(&isotopes, &mut pattern, false, 0.8, 0.02).unwrap(),
        0.01
    );
}

#[test]
fn find_isotope_matches_the_centre_and_neighbour_scans() {
    let spectra = vec![
        spectrum(1.0, 1, &[(100.0, 10.0), (101.0, 11.0)]),
        spectrum(2.0, 1, &[(99.0, 1.0), (100.005, 20.0), (102.0, 3.0)]),
        spectrum(3.0, 1, &[]),
    ];
    let mut pattern = IsotopePattern::new(2).unwrap();
    let mut peak_index = 0;
    let work = find_isotope(&spectra, 100.0, 1, &mut pattern, 0, &mut peak_index, 0.02).unwrap();
    assert_eq!(peak_index, 1);
    assert!(work >= 2);
    assert_eq!(pattern.peak[0], PatternPeak::Found(1));
    assert_eq!(pattern.spectrum[0], 1);
    assert_eq!(pattern.theoretical_mz[0], 100.0);
    let centre = position_score(100.0, 100.005, 0.02);
    let before = position_score(100.0, 100.0, 0.02);
    assert_eq!(pattern.mz_score[0], (centre + before) / 2.0);
    assert_eq!(pattern.intensity[0], (20.0 + 10.0) / 2.0);
    // Nothing near: not found, zeros.
    let mut peak_index = 1;
    find_isotope(&spectra, 100.5, 1, &mut pattern, 1, &mut peak_index, 0.02).unwrap();
    assert_eq!(pattern.peak[1], PatternPeak::NotFound);
    assert_eq!(pattern.intensity[1], 0.0);
    // Out-of-range indices are errors.
    assert!(find_isotope(&spectra, 100.0, 3, &mut pattern, 0, &mut 0, 0.02).is_err());
    assert!(find_isotope(&spectra, 100.0, 1, &mut pattern, 2, &mut 0, 0.02).is_err());
    assert!(find_isotope(&spectra, 100.0, 2, &mut pattern, 0, &mut 0, 0.02).is_err());
}

#[test]
fn intensity_bin_score_interpolates_the_quantiles() {
    // Two spectra, one bin: 21 quantiles of the four intensities.
    let e = experiment(vec![
        spectrum(1.0, 1, &[(100.0, 1.0), (200.0, 2.0)]),
        spectrum(2.0, 1, &[(100.0, 4.0), (200.0, 8.0)]),
    ]);
    let t = IntensityThresholds::compute(&e, 1).unwrap();
    let q = t.quantiles(0, 0).unwrap();
    let values = [1.0, 2.0, 4.0, 8.0];
    for (i, quantile) in q.iter().enumerate() {
        assert_eq!(*quantile, values[(0.05 * i as f64 * 3.0).floor() as usize]);
    }
    assert_eq!(t.bin_score(0, 0, 9.0), Some(1.0));
    // Below the first quantile: 0.05 * 0.5 / 1 - 0.05 < 0, clamped.
    assert_eq!(t.bin_score(0, 0, 0.5), Some(0.0));
    let position = q.partition_point(|&v| v < 3.0);
    let expected = 0.05 * (3.0 - q[position - 1]) / (q[position] - q[position - 1])
        + 0.05 * (position as f64 - 1.0);
    assert_eq!(t.bin_score(0, 0, 3.0), Some(expected));
    assert_eq!(t.bin_score(1, 0, 3.0), None);
    // A single bin interpolates the one cell with itself; at the grid centre all
    // four distances are 1 and all four weights sqrt(2) / (4 sqrt(2)).
    let d = 2.0f64.sqrt();
    let weight = d / (d + d + d + d);
    assert_eq!(
        t.score(1.5, 150.0, 3.0).unwrap(),
        expected * weight + expected * weight + expected * weight + expected * weight
    );
    // Below the range, the position floor(-1.0) converts to UInt as the Linux
    // x86_64 Release build's cvttsd2si does: 0xffffffff, capped at half-bin 1,
    // the last. Both RT neighbours are bin 0 at distance 1, both m/z neighbours
    // bin 0 at distance 0, so each weight is 1 / 4 (the probe test covers such
    // positions on executed grids).
    let quarter = 0.25 * expected;
    assert_eq!(
        t.score(0.5, 150.0, 3.0).unwrap(),
        quarter + quarter + quarter + quarter
    );
    // A NaN retention time selects half-bin 0, and its NaN, sign cleared by the
    // distance's absolute value, reaches the result first.
    assert_eq!(
        t.score(-f64::NAN, 150.0, 3.0).unwrap().to_bits(),
        0x7ff8_0000_0000_0000
    );
    assert!(IntensityThresholds::compute(&e, 0).is_err());
}

#[test]
fn overall_score_is_the_float_cube_root_of_the_product() {
    assert_eq!(overall_score(1.0, 1.0, 1.0), 1.0);
    assert_eq!(overall_score(0.0, 0.5, 1.0), 0.0);
    let product = 0.9f32 * 0.8f32 * 0.7f32;
    assert_eq!(
        overall_score(0.9, 0.8, 0.7),
        libm::pow(f64::from(product), f64::from(1.0f32 / 3.0f32)) as f32
    );
    assert!(overall_score(f32::NAN, 1.0, 1.0).is_nan());
}

// ---------------------------------------------------------------------------
// Degenerate intensity bins: the Linux x86_64 Release build (tier 1)
// ---------------------------------------------------------------------------

/// The rows of `degenerate_stage.tsv.gz`: the driver `degenerate_stage` (the
/// generalised B6 `seed_stage`) run against
/// `openms4-release-bc9cc12-c19e494-174b576` on 26 configurations, three
/// repetitions at one and at four threads, identical
/// (`../oracle/ffap-sem-completion/extract/extract_degenerate.py`).
fn degenerate_rows() -> Vec<Vec<String>> {
    use std::io::Read;
    let bytes = std::fs::read(data("degenerate_stage.tsv.gz")).unwrap();
    let mut text = String::new();
    flate2::read::GzDecoder::new(bytes.as_slice())
        .read_to_string(&mut text)
        .unwrap();
    text.lines()
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect()
}

/// The FeatureFinderCentroided loading of a fixture of the tool's tests.
fn tool_fixture(name: &str) -> MSExperiment {
    let mut options = PeakFileOptions::default();
    options.add_ms_level(1).unwrap();
    options.set_intensity_range(NumericRange {
        min: 0.0,
        max: f64::MAX,
    });
    FileHandler::load_experiment_with_options(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/topp_feature_finder_centroided")
            .join(name),
        &[FileType::MzMl],
        &options,
    )
    .unwrap()
}

/// The experiment and parameters of one `case` row, rebuilt as the driver
/// built them: the mzML (the derived FFC_1 inputs in memory, by the rules of
/// `make_inputs.py`, which the score rows then pin peak for peak), then
/// `keep=`, `rt=`, `rtlast=` and `empty=` in that order, and the parameter
/// overrides on the INI's algorithm section or on no parameters.
fn degenerate_input(case: &[String]) -> (MSExperiment, Param) {
    let mut experiment = match case[0].as_str() {
        "FeatureFinderCentroided_1_input.mzML" => ffc1_input(),
        "zero_rt_ffc1.mzML" => {
            let mut e = ffc1_input();
            for spectrum in &mut e.spectra {
                spectrum.rt = 4114.53;
            }
            e
        }
        name @ ("zero_mz_ffc1.mzML" | "zero_mz_control_ffc1.mzML") => {
            let mut e = ffc1_input();
            for spectrum in &mut e.spectra {
                for peak in &mut spectrum.peaks {
                    peak.mz = 500.0;
                }
            }
            if name == "zero_mz_control_ffc1.mzML" {
                e.spectra[0].peaks[0].mz = 499.0;
            }
            e
        }
        name @ ("FileFilter_44_input.mzML" | "FileConverter_31_output.mzML") => tool_fixture(name),
        other => panic!("unknown input {other}"),
    };
    let mut parameters = match case[1].as_str() {
        "-" => Param::new(),
        "FeatureFinderCentroided_1_parameters.ini" => ffc1_parameters(),
        other => panic!("unknown INI {other}"),
    };
    let option = |prefix: &str| {
        case[2..]
            .iter()
            .filter_map(|o| o.strip_prefix(prefix))
            .collect::<Vec<_>>()
    };
    if let Some(keep) = option("keep=").first() {
        experiment.spectra.truncate(keep.parse().unwrap());
    }
    if let Some(bits) = option("rt=").first() {
        for spectrum in &mut experiment.spectra {
            spectrum.rt = f64_hex(bits);
        }
    }
    if let Some(bits) = option("rtlast=").first() {
        experiment.spectra.last_mut().unwrap().rt = f64_hex(bits);
    }
    for index in option("empty=") {
        experiment.spectra[index.parse::<usize>().unwrap()]
            .peaks
            .clear();
    }
    for assignment in option("i:") {
        let (key, value) = assignment.split_once('=').unwrap();
        set(
            &mut parameters,
            key,
            ParamValue::Integer(value.parse().unwrap()),
        );
    }
    for assignment in option("d:") {
        let (key, value) = assignment.split_once('=').unwrap();
        set(
            &mut parameters,
            key,
            ParamValue::Float(value.parse().unwrap()),
        );
    }
    (experiment, parameters)
}

/// Every configuration of the degenerate captures, compared with the executed
/// library: the input as loaded and modified, the effective members, the bin
/// steps and quantiles, the isotope windows, every per-peak float array with
/// its NaN bits, the double result of `intensityScore_(spectrum, peak)` for
/// every peak, the seeds, the printed lines, the feature count and the abort
/// reasons.
///
/// A zero RT or m/z bin step (every retention time equal, every m/z equal, or
/// a subnormal extent that underflows in the division) and an infinite RT step
/// (an overflowing extent) make every intensity score the default NaN
/// `0xfff8000000000000` (stored as `0xffc00000`), every overall score of the
/// scans the seed loop visits NaN, and the run find nothing, because the source
/// converts `floor(NaN)` and `floor(inf)` to `UInt` as `cvttsd2si` does; the
/// port reproduces that ([`DegenerateBinStep::Source`], the default). The
/// subnormal step with one bin (`tiny_rt`) is not degenerate and finds the 25
/// seeds of FFC_1, which the extreme retention times then make fail the
/// quality check. The short inputs (at most `2 * min_spectra_` scans:
/// `filefilter_44`, `fileconverter_31`, `*_keep10`, `*_keep14`) never reach the
/// seed loop.
///
/// [`DegenerateBinStep::Refuse`] refuses exactly the configurations whose steps
/// the executed build computed as zero or infinite and whose seed loop is not
/// empty.
#[test]
fn degenerate_bin_steps_match_the_linux_release_build() {
    use openms::concept::parallel::Threads;
    let rows = degenerate_rows();
    let cases: Vec<&Vec<String>> = rows.iter().filter(|row| row[0] == "case").collect();
    assert_eq!(cases.len(), 26);
    let mut refused = 0;
    let mut nan_configurations = 0;
    for case in cases {
        let config = case[1].as_str();
        let of = |kind: &str| -> Vec<&[String]> {
            rows.iter()
                .filter(|row| row[0] == kind && row[1] == config)
                .map(|row| &row[2..])
                .collect()
        };
        let (experiment, parameters) = degenerate_input(&case[2..]);

        // The input as the executed build saw it.
        let input = of("input")[0];
        let peaks: usize = experiment.spectra.iter().map(|s| s.peaks.len()).sum();
        assert_eq!(experiment.spectra.len().to_string(), input[0], "{config}");
        assert_eq!(peaks.to_string(), input[1], "{config}");
        assert_eq!(
            experiment.spectra[0].rt.to_bits(),
            f64_hex(&input[2]).to_bits(),
            "{config}"
        );
        assert_eq!(
            experiment.spectra.last().unwrap().rt.to_bits(),
            f64_hex(&input[3]).to_bits(),
            "{config}"
        );
        for row in of("rt") {
            let s: usize = row[0].parse().unwrap();
            assert_eq!(
                experiment.spectra[s].rt.to_bits(),
                f64_hex(&row[1]).to_bits(),
                "{config}: rt of {s}"
            );
            assert_eq!(experiment.spectra[s].native_id, row[2], "{config}");
        }

        let stage = SeedStage::run(experiment.clone(), &FeatureMap::new(), &parameters)
            .unwrap_or_else(|error| panic!("{config}: {error}"))
            .unwrap();
        let settings = stage.settings();
        let members: BTreeMap<&str, &str> = of("member")
            .into_iter()
            .map(|row| (row[0].as_str(), row[1].as_str()))
            .collect();
        assert_eq!(
            settings.min_spectra.to_string(),
            members["min_spectra"],
            "{config}"
        );
        assert_eq!(
            settings.intensity_bins.to_string(),
            members["intensity_bins"],
            "{config}"
        );
        assert_eq!(
            settings.pattern_tolerance.to_bits(),
            f64_hex(members["pattern_tolerance"]).to_bits(),
            "{config}"
        );
        assert_eq!(
            settings.trace_tolerance.to_bits(),
            f64_hex(members["trace_tolerance"]).to_bits(),
            "{config}"
        );

        // Step 1: the bins, as the executed build computed them.
        let thresholds = stage.thresholds();
        let bins = of("bins")[0];
        assert_eq!(thresholds.bins().to_string(), bins[0], "{config}");
        assert_eq!(thresholds.rt_start().to_bits(), f64_hex(&bins[1]).to_bits());
        assert_eq!(thresholds.mz_start().to_bits(), f64_hex(&bins[2]).to_bits());
        assert_eq!(
            thresholds.rt_step().to_bits(),
            f64_hex(&bins[3]).to_bits(),
            "{config}: rt step"
        );
        assert_eq!(
            thresholds.mz_step().to_bits(),
            f64_hex(&bins[4]).to_bits(),
            "{config}: m/z step"
        );
        let quantiles = of("quantiles");
        assert_eq!(quantiles.len(), thresholds.bins() * thresholds.bins());
        for row in quantiles {
            let actual = thresholds
                .quantiles(row[0].parse().unwrap(), row[1].parse().unwrap())
                .unwrap();
            let expected: Vec<u64> = row[2..].iter().map(|q| f64_hex(q).to_bits()).collect();
            let actual: Vec<u64> = actual.iter().map(|q| q.to_bits()).collect();
            assert_eq!(
                actual, expected,
                "{config}: quantiles {}/{}",
                row[0], row[1]
            );
        }

        // Step 2.5.
        let windows = of("window");
        let patterns = stage.windows().patterns();
        assert_eq!(patterns.len(), windows.len(), "{config}: windows");
        for (pattern, row) in patterns.iter().zip(&windows) {
            let bits: Vec<u64> = pattern.intensity.iter().map(|v| v.to_bits()).collect();
            let expected: Vec<u64> = row[7..].iter().map(|v| f64_hex(v).to_bits()).collect();
            assert_eq!(bits, expected, "{config}: window {}", row[0]);
            assert_eq!(pattern.max.to_bits(), f64_hex(&row[5]).to_bits());
        }

        // Every per-peak array, NaN bits included, and intensityScore_.
        let scores = stage.scores();
        let charges = scores.charge_count();
        let rounding: BTreeMap<(usize, usize, usize), (u32, u32)> = of("rounding")
            .into_iter()
            .map(|row| {
                let oracle = u32::from_str_radix(&row[3], 16).unwrap();
                let correct = u32::from_str_radix(&row[4], 16).unwrap();
                assert_eq!(oracle.abs_diff(correct), 1, "{config}: rounding row");
                (
                    (
                        row[0].parse().unwrap(),
                        row[1].parse().unwrap(),
                        row[2].parse().unwrap(),
                    ),
                    (oracle, correct),
                )
            })
            .collect();
        let score_rows = of("score");
        let with_scores = !score_rows.is_empty();
        if with_scores {
            assert_eq!(score_rows.len(), peaks, "{config}: score rows");
        }
        let mut rounded = 0;
        let mut nan_scores = 0;
        for row in score_rows {
            let s: usize = row[0].parse().unwrap();
            let p: usize = row[1].parse().unwrap();
            let spectrum = &stage.experiment().spectra[s];
            let peak = spectrum.peaks[p];
            assert_eq!(peak.mz.to_bits(), f64_hex(&row[2]).to_bits(), "{config}");
            assert_eq!(
                peak.intensity.to_bits(),
                f32_hex(&row[3]).to_bits(),
                "{config}"
            );
            let arrays = &row[4..row.len() - 1];
            assert_eq!(arrays.len(), 3 + 2 * charges, "{config}: arrays");
            let mut actual = vec![
                scores.trace(s).unwrap()[p],
                scores.intensity(s).unwrap()[p],
                scores.local_max(s).unwrap()[p],
            ];
            for c in 0..charges {
                actual.push(scores.pattern(c, s).unwrap()[p]);
            }
            for c in 0..charges {
                actual.push(scores.overall(c, s).unwrap()[p]);
            }
            for (column, (value, oracle)) in actual.iter().zip(arrays).enumerate() {
                let mut expected = u32::from_str_radix(oracle, 16).unwrap();
                if column >= 3 + charges {
                    if let Some(&(misrounded, correct)) =
                        rounding.get(&(s, p, column - 3 - charges))
                    {
                        assert_eq!(misrounded, expected);
                        expected = correct;
                        rounded += 1;
                    }
                }
                assert_eq!(
                    value.to_bits(),
                    expected,
                    "{config}: spectrum {s} peak {p} array {column}"
                );
            }
            let score = thresholds
                .score(spectrum.rt, peak.mz, f64::from(peak.intensity))
                .unwrap();
            assert_eq!(
                score.to_bits(),
                f64_hex(&row[row.len() - 1]).to_bits(),
                "{config}: intensityScore_({s}, {p})"
            );
            nan_scores += usize::from(score.is_nan());
        }
        let degenerate = [&bins[3], &bins[4]].iter().any(|step| {
            let step = f64_hex(step);
            step == 0.0 || step.is_infinite()
        });
        if with_scores {
            assert_eq!(rounded, rounding.len(), "{config}: rounding rows");
            // Every intensity score is NaN exactly when a step is degenerate.
            assert_eq!(nan_scores, if degenerate { peaks } else { 0 }, "{config}");
            nan_configurations += usize::from(degenerate);
        }

        // Seeds, with the executed (possibly misrounded) overall score.
        let seeds = of("seed");
        let mut actual_seeds = Vec::new();
        for charge in stage.charges() {
            let index = (charge.charge - settings.charge_low) as usize;
            for (rank, seed) in charge.seeds.iter().enumerate() {
                let overall = scores.overall(index, seed.spectrum).unwrap()[seed.peak].to_bits();
                let executed = rounding
                    .get(&(seed.spectrum, seed.peak, index))
                    .map_or(overall, |&(misrounded, _)| misrounded);
                actual_seeds.push(vec![
                    charge.charge.to_string(),
                    rank.to_string(),
                    seed.spectrum.to_string(),
                    seed.peak.to_string(),
                    format!("{:08x}", seed.intensity.to_bits()),
                    format!("{executed:08x}"),
                ]);
            }
        }
        assert_eq!(actual_seeds.len(), seeds.len(), "{config}: seed count");
        for (actual, expected) in actual_seeds.iter().zip(&seeds) {
            assert_eq!(actual.as_slice(), *expected, "{config}: seed");
        }

        // The whole run: printed lines, feature count, abort reasons.
        let output = run_with_options(
            experiment.clone(),
            &FeatureMap::new(),
            &parameters,
            &Options {
                threads: Threads::serial(),
                ..Options::default()
            },
        )
        .unwrap_or_else(|error| panic!("{config}: {error}"));
        let printed: Vec<&String> = output
            .log
            .iter()
            .filter(|line| line.starts_with("Found "))
            .collect();
        let stdout: Vec<&String> = of("stdout").into_iter().map(|row| &row[0]).collect();
        assert_eq!(printed, stdout, "{config}: printed lines");
        assert_eq!(
            output.features.len().to_string(),
            of("features")[0][0],
            "{config}: features"
        );
        let aborts: BTreeMap<String, usize> = of("abort")
            .into_iter()
            .map(|row| (row[1].clone(), row[0].parse().unwrap()))
            .collect();
        assert_eq!(output.aborts, aborts, "{config}: abort reasons");

        // The opt-out refuses exactly the undefined cases the seed loop reads.
        let spectra = experiment.spectra.len();
        let loop_end = spectra - settings.min_spectra.min(spectra);
        let read = settings.min_spectra < loop_end;
        let outcome = run_with_options(
            experiment,
            &FeatureMap::new(),
            &parameters,
            &Options {
                threads: Threads::serial(),
                degenerate_bin_step: DegenerateBinStep::Refuse,
                ..Options::default()
            },
        );
        if degenerate && read {
            refused += 1;
            assert!(
                matches!(&outcome, Err(Error::InvalidValue(message)) if message.contains("DegenerateBinStep::Refuse")),
                "{config}: {:?}",
                outcome.map(|o| o.features.len())
            );
        } else {
            let refusing = outcome.unwrap_or_else(|error| panic!("{config}: {error}"));
            assert_eq!(refusing.log, output.log, "{config}");
            assert_eq!(refusing.features.len(), output.features.len(), "{config}");
        }
    }
    // zero_rt, its min_score 0 variant, its defaults, keep15 and keep15_empty7
    // (the FFC_1 INI's min_spectra_ is 7, so keep10 to keep14 are short);
    // zero_mz, its min_score 0 variant and keep15; tiny_rt with the default
    // 10 bins (three configurations) and with 2 bins; huge_rt, its min_score 0
    // variant and its defaults.
    assert_eq!(refused, 15);
    // The score rows of 17 configurations: all but tiny_rt, fileconverter_31
    // and zero_mz_control_min_score_0 are degenerate.
    assert_eq!(nan_configurations, 17);
}

/// Driver `iscore_probe`: `intensityScore_(spectrum, peak)` of the Linux x86_64
/// Release build at 57 positions on four grids of 11 spectra by 11 peaks
/// (default parameters, 10 bins), including positions outside the binned range
/// (negative, beyond `2^31`, `2^32` and `2^63`, infinite and NaN) whose
/// `UInt` conversion is undefined, NaN and infinite intensities, and the grids
/// with a zero RT step, a zero m/z step and a subnormal RT extent that
/// underflows to a zero step. The probe moved spectrum 0, peak 0 to each
/// position without updating the ranges. Every result is compared bit for bit,
/// NaN sign and payload included.
#[test]
fn intensity_scores_outside_the_bins_match_the_linux_release_build() {
    use std::io::Read;
    let bytes = std::fs::read(data("iscore_probe.tsv.gz")).unwrap();
    let mut text = String::new();
    flate2::read::GzDecoder::new(bytes.as_slice())
        .read_to_string(&mut text)
        .unwrap();
    let rows: Vec<Vec<&str>> = text
        .lines()
        .map(|line| line.split('\t').collect())
        .collect();
    let mut queries = 0;
    for grid in ["grid", "zero_rt", "zero_mz", "tiny_rt"] {
        // iscore_probe.cpp make_grid.
        let spectra = (0..=10)
            .map(|s| {
                let rt = match grid {
                    "zero_rt" => 5.0,
                    "tiny_rt" if s == 10 => f64::from_bits(1),
                    "tiny_rt" => 0.0,
                    _ => f64::from(s),
                };
                let peaks: Vec<(f64, f32)> = (0..=10)
                    .map(|p| {
                        let mz = if grid == "zero_mz" {
                            150.0
                        } else {
                            100.0 + 10.0 * f64::from(p)
                        };
                        (mz, ((s + 1) * (p + 1)) as f32)
                    })
                    .collect();
                spectrum(rt, 1, &peaks)
            })
            .collect();
        let stage = SeedStage::run(experiment(spectra), &FeatureMap::new(), &Param::new())
            .unwrap()
            .unwrap();
        let thresholds = stage.thresholds();
        let of = |kind: &'static str| {
            rows.iter()
                .filter(move |row| row[0] == kind && row[1] == grid)
        };
        let header = of("grid").next().unwrap();
        assert_eq!(thresholds.bins().to_string(), header[2]);
        assert_eq!(
            thresholds.rt_step().to_bits(),
            f64_hex(header[3]).to_bits(),
            "{grid}"
        );
        assert_eq!(
            thresholds.mz_step().to_bits(),
            f64_hex(header[4]).to_bits(),
            "{grid}"
        );
        for row in of("quantiles") {
            let actual = thresholds
                .quantiles(row[2].parse().unwrap(), row[3].parse().unwrap())
                .unwrap();
            for (value, expected) in actual.iter().zip(&row[4..]) {
                assert_eq!(value.to_bits(), f64_hex(expected).to_bits(), "{grid}");
            }
        }
        for row in of("query") {
            let rt = f64_hex(row[2]);
            let mz = f64_hex(row[3]);
            let intensity = f32_hex(row[4]);
            let score = thresholds.score(rt, mz, f64::from(intensity)).unwrap();
            assert_eq!(
                format!("{:016x}", score.to_bits()),
                row[5],
                "{grid}: rt {rt:e}, m/z {mz:e}, intensity {intensity:e}"
            );
            queries += 1;
        }
    }
    assert_eq!(queries, 4 * 57);
}
