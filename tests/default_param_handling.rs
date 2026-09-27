// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The `DefaultParamHandler` contract of every type implementing
//! [`DefaultParamHandling`], against the Linux x86-64 Release build.
//!
//! Tier 1: `tests/data/param_defaults/*.defaults.ini` is `getDefaults()` of a
//! default-constructed object written through `ParamXMLFile`, and
//! `<case>.ini` is `getParameters()` after the case's `setParameters` (or
//! constructor), both produced by `../oracle/B6-param-defaults` on ibminode06
//! against the Release install `openms4-release-bc9cc12-c19e494-174b576`, every
//! case run twice in fresh processes with identical output.
//! `oracle_outcomes.json` holds each case's stdout (exception name and
//! message), stderr (the warnings) and exit code. See
//! `tests/data/param_defaults_provenance.json`.

#![cfg(feature = "paramxml")]

use openms::Error;
use openms::format::paramxml;
use openms::param::{DefaultParamHandling, Param, ParamValue};
use openms::processing::iterative::PeakPickerIterative;
use openms::processing::mean_noise::{
    MeanNoiseHistogramRange, SignalToNoiseEstimatorMeanIterative,
};
use openms::processing::peak_picking::{PeakPickerHiRes, SignalToNoiseEstimatorMedian};
use openms::processing::smoothing::{GaussFilter, GaussianWidth, SavitzkyGolayFilter};
use openms::processing::window_mower::{WindowMower, WindowMowerMethod};
use openms::processing::{
    LinearResamplerAlign, NLargest, NormalizationMethod, Normalizer, RankScaler, SqrtScaler,
    ThresholdMower,
};
use std::collections::BTreeMap;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/param_defaults")
        .join(name)
}

fn oracle(name: &str) -> Param {
    paramxml::load(data(name)).unwrap_or_else(|e| panic!("{name}: {e}"))
}

/// One oracle case: exit code, stdout and stderr.
struct Outcome {
    stdout: String,
    stderr: String,
}

fn outcomes() -> BTreeMap<String, Outcome> {
    let text = std::fs::read_to_string(data("oracle_outcomes.json")).unwrap();
    let json: serde_json_lite::Value = serde_json_lite::parse(&text);
    json.entries()
        .into_iter()
        .map(|(case, value)| {
            assert_eq!(value.get("exit").number(), 0.0, "{case} exit");
            (
                case,
                Outcome {
                    stdout: value.get("stdout").string(),
                    stderr: value.get("stderr").string(),
                },
            )
        })
        .collect()
}

fn outcome(case: &str) -> Outcome {
    outcomes()
        .remove(case)
        .unwrap_or_else(|| panic!("no oracle case {case}"))
}

fn param(entries: &[(&str, ParamValue)]) -> Param {
    let mut param = Param::new();
    for (key, value) in entries {
        param.set_value(key, value.clone(), "", &[]).unwrap();
    }
    param
}

/// The whole contract of one type against its oracle dumps.
fn contract<T: DefaultParamHandling + Default>(name: &str) {
    let defaults = T::defaults().unwrap();
    // Tier 1: names, values, types, descriptions, restrictions, tags, order.
    assert_eq!(defaults, oracle(&format!("{name}.defaults.ini")), "{name}");
    assert_eq!(T::NAME, name);
    // `defaultsToParam_()`: the Release build's `getParameters()` of a
    // default-constructed object equals its `getDefaults()`.
    assert!(
        outcome(&format!("defaults_{name}"))
            .stdout
            .contains("defaults_equal_parameters\t1")
    );
    assert_eq!(T::default().to_param().unwrap(), defaults, "{name}");
    let rebuilt = T::from_param(&defaults).unwrap();
    assert_eq!(rebuilt.to_param().unwrap(), defaults, "{name}");
    let (_, handler, warnings) = T::from_param_with_handler(&Param::new()).unwrap();
    assert_eq!(handler.defaults(), &defaults);
    assert_eq!(handler.name(), name);
    let empty = defaults.is_empty();
    assert_eq!(warnings.len(), usize::from(empty), "{name}: {warnings:?}");

    // setParameters with one unknown key: the source warns and keeps it.
    let unknown = param(&[("unknown_key", ParamValue::Integer(1))]);
    let (_, handler, warnings) = T::from_param_with_handler(&unknown).unwrap();
    assert_eq!(
        handler.parameters(),
        &oracle(&format!("unknown_{name}.ini")),
        "{name}: getParameters() after setParameters"
    );
    let expected = outcome(&format!("unknown_{name}"));
    let lines: Vec<&str> = expected.stderr.lines().collect();
    assert_eq!(warnings.len(), lines.len(), "{name}: {warnings:?}");
    for (warning, line) in warnings.iter().zip(&lines) {
        if line.contains("No default parameters") {
            assert_eq!(warning, line);
        } else {
            // The source's `Warning: <name> received the unknown parameter
            // '<key>'!`; the handler words it `<name>: unknown parameter
            // '<key>'` (see DefaultParamHandler::set_parameters).
            assert!(line.contains("'unknown_key'"), "{line}");
            assert_eq!(warning, &format!("{name}: unknown parameter 'unknown_key'"));
        }
    }
    assert!(expected.stdout.contains("END setParameters ok"));
}

#[test]
fn nlargest_contract_and_constructor() {
    contract::<NLargest>("NLargest");
    // `NLargest(7)` sets `n` after `defaultsToParam_()`.
    assert_eq!(
        NLargest { n: 7 }.to_param().unwrap(),
        oracle("ctor_NLargest7.ini")
    );
    let filter = NLargest::from_param(&param(&[("n", ParamValue::Integer(3))])).unwrap();
    assert_eq!(filter.n, 3);
    // No minimum, so setParameters accepts it; `(UInt)` then throws
    // ConversionError in updateMembers_.
    assert!(matches!(
        NLargest::from_param(&param(&[("n", ParamValue::Integer(-1))])),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        NLargest {
            n: i32::MAX as usize + 1
        }
        .to_param(),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn threshold_mower_contract() {
    contract::<ThresholdMower>("ThresholdMower");
    let mower =
        ThresholdMower::from_param(&param(&[("threshold", ParamValue::Float(2.5))])).unwrap();
    assert_eq!(mower.threshold, 2.5);
    assert_eq!(
        mower.to_param().unwrap().value("threshold").unwrap(),
        &ParamValue::Float(2.5)
    );
    // An integer where the default is a double is the wrong type.
    assert!(ThresholdMower::from_param(&param(&[("threshold", ParamValue::Integer(1))])).is_err());
}

#[test]
fn window_mower_contract() {
    contract::<WindowMower>("WindowMower");
    let changed = param(&[
        ("windowsize", ParamValue::Float(10.0)),
        ("peakcount", ParamValue::Integer(4)),
        ("movetype", ParamValue::String("jump".into())),
    ]);
    let mower = WindowMower::from_param(&changed).unwrap();
    assert_eq!(
        (mower.window_size, mower.peak_count, mower.method),
        (10.0, 4, WindowMowerMethod::Jumping)
    );
    let mut expected = WindowMower::defaults().unwrap();
    expected
        .replace_value("windowsize", ParamValue::Float(10.0))
        .unwrap();
    expected
        .replace_value("peakcount", ParamValue::Integer(4))
        .unwrap();
    expected
        .replace_value("movetype", ParamValue::String("jump".into()))
        .unwrap();
    assert_eq!(mower.to_param().unwrap(), expected);
    assert!(
        WindowMower::from_param(&param(&[("movetype", ParamValue::String("walk".into()))]))
            .is_err()
    );
    // set_parameters keeps the native ceilings and resets unnamed parameters
    // to their defaults, as the source's setParameters does.
    let mut native = WindowMower {
        max_points: 17,
        window_size: 3.0,
        ..WindowMower::default()
    };
    native
        .set_parameters(&param(&[("peakcount", ParamValue::Integer(9))]))
        .unwrap();
    assert_eq!(
        (native.max_points, native.window_size, native.peak_count),
        (17, 50.0, 9)
    );
}

#[test]
fn normalizer_contract() {
    contract::<Normalizer>("Normalizer");
    let tic =
        Normalizer::from_param(&param(&[("method", ParamValue::String("to_TIC".into()))])).unwrap();
    assert_eq!(tic.method, NormalizationMethod::ToTic);
    assert_eq!(
        tic.to_param().unwrap().value("method").unwrap(),
        &ParamValue::String("to_TIC".into())
    );
    // The Release build throws InvalidParameter from the restriction check.
    let bad = param(&[("method", ParamValue::String("to_max".into()))]);
    assert!(matches!(
        Normalizer::from_param(&bad),
        Err(Error::InvalidValue(_))
    ));
    let expected = outcome("set_Normalizer_bad_method");
    assert!(
        expected
            .stdout
            .contains("InvalidParameter\tNormalizer: Invalid string parameter value 'to_max'")
    );
}

#[test]
fn scalers_declare_no_parameters() {
    contract::<RankScaler>("RankScaler");
    contract::<SqrtScaler>("SqrtScaler");
    for (name, warnings) in [
        (
            "RankScaler",
            RankScaler::from_param_with_warnings(&Param::new())
                .unwrap()
                .1,
        ),
        (
            "SqrtScaler",
            SqrtScaler::from_param_with_warnings(&Param::new())
                .unwrap()
                .1,
        ),
    ] {
        let expected = outcome(&format!("empty_{name}"));
        assert_eq!(warnings, expected.stderr.lines().collect::<Vec<_>>());
        assert_eq!(
            oracle(&format!("empty_{name}.ini")),
            Param::new(),
            "{name}: getParameters() stays empty"
        );
    }
}

#[test]
fn gauss_filter_contract() {
    contract::<GaussFilter>("GaussFilter");
    let ppm = param(&[
        ("use_ppm_tolerance", ParamValue::String("true".into())),
        ("ppm_tolerance", ParamValue::Float(5.0)),
        ("gaussian_width", ParamValue::Float(0.4)),
        ("write_log_messages", ParamValue::String("true".into())),
    ]);
    let filter = GaussFilter::from_param(&ppm).unwrap();
    assert_eq!(filter.algorithm.width, GaussianWidth::Ppm(5.0));
    assert_eq!(filter.algorithm.kernel_spacing, 0.01);
    assert_eq!(filter.inactive_width, 0.4);
    assert!(filter.write_log_messages);
    let round = filter.to_param().unwrap();
    for (key, value) in [
        ("use_ppm_tolerance", ParamValue::String("true".into())),
        ("ppm_tolerance", ParamValue::Float(5.0)),
        ("gaussian_width", ParamValue::Float(0.4)),
        ("write_log_messages", ParamValue::String("true".into())),
    ] {
        assert_eq!(round.value(key).unwrap(), &value, "{key}");
    }
    assert_eq!(
        GaussFilter::from_param(&round).unwrap().to_param().unwrap(),
        round
    );
    // Restriction `gaussian_width >= 0`: the Release build throws
    // InvalidParameter and leaves the refused value in param_.
    let negative = param(&[("gaussian_width", ParamValue::Float(-1.0))]);
    assert!(matches!(
        GaussFilter::from_param(&negative),
        Err(Error::InvalidValue(_))
    ));
    assert!(
        outcome("set_GaussFilter_negative")
            .stdout
            .contains("InvalidParameter\tGaussFilter: Invalid double parameter value '-1.000000'")
    );
    assert_eq!(
        oracle("set_GaussFilter_negative.ini")
            .value("gaussian_width")
            .unwrap(),
        &ParamValue::Float(-1.0)
    );
    // Zero passes the restriction; the source divides by it when filtering.
    let zero = GaussFilter::from_param(&param(&[("gaussian_width", ParamValue::Float(0.0))]));
    assert_eq!(zero.unwrap().algorithm.width, GaussianWidth::Absolute(0.0));
    // set_parameters keeps the native lookup spacing and coefficient ceiling.
    let mut native = GaussFilter::default();
    native.algorithm.max_coefficients = 99;
    native.set_parameters(&ppm).unwrap();
    assert_eq!(native.algorithm.max_coefficients, 99);
    assert_eq!(native.algorithm.width, GaussianWidth::Ppm(5.0));
}

#[test]
fn savitzky_golay_contract() {
    contract::<SavitzkyGolayFilter>("SavitzkyGolayFilter");
    // An even frame length is incremented for the filter and kept in param_.
    let even = param(&[("frame_length", ParamValue::Integer(10))]);
    let (filter, handler, warnings) = SavitzkyGolayFilter::from_param_with_handler(&even).unwrap();
    assert!(warnings.is_empty());
    assert_eq!(filter.frame_length(), 11);
    assert_eq!(filter.polynomial_order(), 4);
    assert_eq!(
        handler.parameters(),
        &oracle("set_SavitzkyGolayFilter_even.ini")
    );
    assert_eq!(
        filter.to_param().unwrap().value("frame_length").unwrap(),
        &ParamValue::Integer(10)
    );
    assert!(
        outcome("set_SavitzkyGolayFilter_even")
            .stdout
            .contains("END setParameters ok")
    );
    // updateMembers_ throws InvalidValue when the degree is not below the
    // frame length, with this what() in the Release build.
    let order = param(&[
        ("frame_length", ParamValue::Integer(4)),
        ("polynomial_order", ParamValue::Integer(5)),
    ]);
    let expected = outcome("set_SavitzkyGolayFilter_order");
    let message = expected
        .stdout
        .lines()
        .find_map(|line| line.strip_prefix("END setParameters exception\tInvalidValue\t"))
        .unwrap()
        .to_owned();
    match SavitzkyGolayFilter::from_param(&order) {
        Err(Error::InvalidValue(text)) => assert_eq!(text, message),
        other => panic!("{other:?}"),
    }
    assert!(
        SavitzkyGolayFilter::from_param(&param(&[("polynomial_order", ParamValue::Integer(-2))]))
            .is_err()
    );
}

#[test]
fn linear_resampler_align_contract() {
    contract::<LinearResamplerAlign>("LinearResamplerAlign");
    let resampler =
        LinearResamplerAlign::from_param(&param(&[("spacing", ParamValue::Float(0.5))])).unwrap();
    assert_eq!(resampler.spacing, 0.5);
    // `ppm` has no valid-strings restriction; `toBool()` throws for anything
    // but true/false, and the ppm grid is not ported.
    assert!(matches!(
        LinearResamplerAlign::from_param(&param(&[("ppm", ParamValue::String("yes".into()))])),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        LinearResamplerAlign::from_param(&param(&[("ppm", ParamValue::String("true".into()))])),
        Err(Error::Unsupported(_))
    ));
    let mut native = LinearResamplerAlign {
        max_points: 5,
        ..LinearResamplerAlign::default()
    };
    native
        .set_parameters(&param(&[("spacing", ParamValue::Float(2.0))]))
        .unwrap();
    assert_eq!((native.spacing, native.max_points), (2.0, 5));
}

#[test]
fn mean_iterative_contract() {
    contract::<SignalToNoiseEstimatorMeanIterative>("SignalToNoiseEstimatorMeanIterative");
    for (mode, expected) in [
        (
            -1,
            MeanNoiseHistogramRange::Manual {
                max_intensity: 250.0,
            },
        ),
        (
            0,
            MeanNoiseHistogramRange::StandardDeviation { factor: 2.5 },
        ),
        (
            1,
            MeanNoiseHistogramRange::LegacyPercentile { percentile: 80 },
        ),
    ] {
        let given = param(&[
            ("auto_mode", ParamValue::Integer(mode)),
            ("max_intensity", ParamValue::Integer(250)),
            ("auto_max_stdev_factor", ParamValue::Float(2.5)),
            ("auto_max_percentile", ParamValue::Integer(80)),
            ("win_len", ParamValue::Float(40.0)),
            ("bin_count", ParamValue::Integer(12)),
            ("min_required_elements", ParamValue::Integer(4)),
        ]);
        let estimator = SignalToNoiseEstimatorMeanIterative::from_param(&given).unwrap();
        assert_eq!(estimator.histogram_range, expected);
        assert_eq!(
            (
                estimator.window_length,
                estimator.bin_count,
                estimator.min_required_elements
            ),
            (40.0, 12, 4)
        );
        // Every value survives, including the two the mode does not use.
        let round = estimator.to_param().unwrap();
        for item in [
            "auto_mode",
            "max_intensity",
            "auto_max_stdev_factor",
            "auto_max_percentile",
            "win_len",
            "bin_count",
            "min_required_elements",
        ] {
            assert_eq!(
                round.value(item).unwrap(),
                given.value(item).unwrap(),
                "{item}"
            );
        }
    }
    // auto_mode 2 breaks the restriction -1..1.
    assert!(
        SignalToNoiseEstimatorMeanIterative::from_param(&param(&[(
            "auto_mode",
            ParamValue::Integer(2)
        )]))
        .is_err()
    );
    let mut native = SignalToNoiseEstimatorMeanIterative {
        max_bins: 64,
        ..Default::default()
    };
    native
        .set_parameters(&param(&[("bin_count", ParamValue::Integer(40))]))
        .unwrap();
    assert_eq!((native.max_bins, native.bin_count), (64, 40));
}

#[test]
fn peak_picker_iterative_contract() {
    contract::<PeakPickerIterative>("PeakPickerIterative");
    let given = param(&[
        ("signal_to_noise_", ParamValue::Float(0.0)),
        ("peak_width", ParamValue::Float(0.02)),
        ("spacing_difference", ParamValue::Float(2.5)),
        ("sn_bin_count_", ParamValue::Integer(40)),
        ("nr_iterations_", ParamValue::Integer(2)),
        ("sn_win_len_", ParamValue::Float(35.0)),
        ("check_width_internally", ParamValue::String("true".into())),
        ("ms1_only", ParamValue::String("true".into())),
        ("clear_meta_data", ParamValue::String("true".into())),
    ]);
    let picker = PeakPickerIterative::from_param(&given).unwrap();
    assert_eq!(picker.signal_to_noise, 0.0);
    assert_eq!(picker.peak_width, 0.02);
    assert_eq!(picker.spacing_difference, 2.5);
    assert_eq!(picker.noise_estimator.bin_count, 40);
    assert_eq!(picker.noise_estimator.window_length, 35.0);
    assert_eq!(picker.iterations, 2);
    assert!(picker.check_width_internally && picker.ms1_only && picker.clear_meta_data);
    // The rest of the refinement estimator is a default median estimator.
    let median = SignalToNoiseEstimatorMedian::default();
    assert_eq!(
        picker.noise_estimator.min_required_elements,
        median.min_required_elements
    );
    let round = picker.to_param().unwrap();
    let mut expected = PeakPickerIterative::defaults().unwrap();
    for key in [
        "signal_to_noise_",
        "peak_width",
        "spacing_difference",
        "sn_bin_count_",
        "nr_iterations_",
        "sn_win_len_",
        "check_width_internally",
        "ms1_only",
        "clear_meta_data",
    ] {
        expected
            .replace_value(key, given.value(key).unwrap().clone())
            .unwrap();
    }
    assert_eq!(round, expected);
    // `nr_iterations_ >= 1`.
    assert!(matches!(
        PeakPickerIterative::from_param(&param(&[("nr_iterations_", ParamValue::Integer(0))])),
        Err(Error::InvalidValue(_))
    ));
    assert!(
        outcome("set_PeakPickerIterative_zero_iterations")
            .stdout
            .contains("InvalidParameter\tPeakPickerIterative: Invalid integer parameter value '0'")
    );
    // set_parameters keeps the estimator's non-parameter members.
    let mut native = PeakPickerIterative::default();
    native.noise_estimator.min_required_elements = 3;
    native.max_points = 11;
    native.set_parameters(&given).unwrap();
    assert_eq!(native.noise_estimator.min_required_elements, 3);
    assert_eq!(native.noise_estimator.bin_count, 40);
    assert_eq!(native.max_points, 11);
}

#[test]
fn adopted_reference_types_meet_the_contract() {
    contract::<PeakPickerHiRes>("PeakPickerHiRes");
    contract::<SignalToNoiseEstimatorMedian>("SignalToNoiseEstimatorMedian");
    // The inherent entry points and the trait agree.
    let given = param(&[("signal_to_noise", ParamValue::Float(2.0))]);
    let inherent = PeakPickerHiRes::from_param_with_warnings(&given).unwrap();
    let trait_path =
        <PeakPickerHiRes as DefaultParamHandling>::from_param_with_warnings(&given).unwrap();
    assert_eq!(inherent, trait_path);
    // set_parameters keeps the native members of the picker and its estimator.
    let mut picker = PeakPickerHiRes {
        check_spectrum_type: false,
        max_points: 12,
        ..PeakPickerHiRes::default()
    };
    picker.noise_estimator.max_bins = 77;
    DefaultParamHandling::set_parameters(&mut picker, &given).unwrap();
    assert_eq!(picker.signal_to_noise, 2.0);
    assert!(!picker.check_spectrum_type);
    assert_eq!(picker.max_points, 12);
    assert_eq!(picker.noise_estimator.max_bins, 77);
}

/// The handler follows the source's atomicity only where it is safe: a
/// refused tree changes nothing, where the source keeps the refused values in
/// `param_` (the oracle's dump after the failed call shows them).
#[test]
fn refused_parameters_leave_the_value_unchanged() {
    let mut mower = WindowMower {
        peak_count: 7,
        ..WindowMower::default()
    };
    let before = mower;
    let refused = param(&[("movetype", ParamValue::String("walk".into()))]);
    assert!(mower.set_parameters(&refused).is_err());
    assert_eq!(mower, before);
    assert_eq!(
        oracle("set_PeakPickerIterative_zero_iterations.ini")
            .value("nr_iterations_")
            .unwrap(),
        &ParamValue::Integer(0)
    );
}

/// A minimal JSON reader for the oracle outcome file, so that this binary
/// needs no feature beyond `paramxml`.
mod serde_json_lite {
    pub enum Value {
        Object(Vec<(String, Value)>),
        String(String),
        Number(f64),
    }

    impl Value {
        pub fn entries(self) -> Vec<(String, Value)> {
            match self {
                Value::Object(entries) => entries,
                _ => panic!("not an object"),
            }
        }
        pub fn get(&self, key: &str) -> &Value {
            match self {
                Value::Object(entries) => {
                    &entries
                        .iter()
                        .find(|(name, _)| name == key)
                        .unwrap_or_else(|| panic!("no key {key}"))
                        .1
                }
                _ => panic!("not an object"),
            }
        }
        pub fn string(&self) -> String {
            match self {
                Value::String(text) => text.clone(),
                _ => panic!("not a string"),
            }
        }
        pub fn number(&self) -> f64 {
            match self {
                Value::Number(value) => *value,
                _ => panic!("not a number"),
            }
        }
    }

    pub fn parse(text: &str) -> Value {
        let chars: Vec<char> = text.chars().collect();
        let mut at = 0;
        let value = value(&chars, &mut at);
        skip(&chars, &mut at);
        assert_eq!(at, chars.len(), "trailing JSON");
        value
    }

    fn skip(chars: &[char], at: &mut usize) {
        while *at < chars.len() && chars[*at].is_whitespace() {
            *at += 1;
        }
    }

    fn value(chars: &[char], at: &mut usize) -> Value {
        skip(chars, at);
        match chars[*at] {
            '{' => {
                *at += 1;
                let mut entries = Vec::new();
                loop {
                    skip(chars, at);
                    if chars[*at] == '}' {
                        *at += 1;
                        return Value::Object(entries);
                    }
                    let key = match value(chars, at) {
                        Value::String(key) => key,
                        _ => panic!("object key"),
                    };
                    skip(chars, at);
                    assert_eq!(chars[*at], ':');
                    *at += 1;
                    entries.push((key, value(chars, at)));
                    skip(chars, at);
                    if chars[*at] == ',' {
                        *at += 1;
                    }
                }
            }
            '"' => {
                *at += 1;
                let mut text = String::new();
                loop {
                    let c = chars[*at];
                    *at += 1;
                    match c {
                        '"' => return Value::String(text),
                        '\\' => {
                            let escaped = chars[*at];
                            *at += 1;
                            match escaped {
                                'n' => text.push('\n'),
                                't' => text.push('\t'),
                                'r' => text.push('\r'),
                                '"' => text.push('"'),
                                '\\' => text.push('\\'),
                                '/' => text.push('/'),
                                'u' => {
                                    let hex: String = chars[*at..*at + 4].iter().collect();
                                    *at += 4;
                                    let code = u32::from_str_radix(&hex, 16).unwrap();
                                    text.push(char::from_u32(code).unwrap());
                                }
                                other => panic!("escape {other}"),
                            }
                        }
                        c => text.push(c),
                    }
                }
            }
            _ => {
                let start = *at;
                while *at < chars.len() && (chars[*at] == '-' || chars[*at].is_ascii_digit()) {
                    *at += 1;
                }
                let number: String = chars[start..*at].iter().collect();
                Value::Number(number.parse().unwrap())
            }
        }
    }
}
