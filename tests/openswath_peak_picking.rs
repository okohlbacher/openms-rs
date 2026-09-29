// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Chromatogram peak picking on the OpenSWATH path: the `DefaultParamHandler`
//! surfaces of `PeakPickerChromatogram` and `PeakIntegrator`, the integrator's
//! index-bounded (`ConstIterator&`) overloads, and `MRMTransitionGroupPicker`.
//!
//! Expected values are transcribed from the class tests at OpenMS4-core
//! `bc9cc12` (`PeakPickerChromatogram_test.cpp`, `PeakIntegrator_test.cpp`,
//! `MRMTransitionGroupPicker_test.cpp`) or read from the source constructors;
//! provenance is in `tests/data/openswath_peak_picking_provenance.json`. The
//! shared chromatogram literals are the exact-bit TSVs of
//! `tests/data/chromatogram_processing_*`.

use openms::analysis::openswath::masstrace_correlator::{
    MASSTRACE_CORRELATOR_NAME, MasstraceCorrelator,
};
use openms::analysis::openswath::mrm_transition_group_picker::{
    BackgroundSubtraction, BoundarySelection, MRMTransitionGroupPicker, PeakIntegrationSignal,
};
use openms::analysis::openswath::peak_picker_mobilogram::{PeakPickerMobilogram, PeakPositions};
use openms::analysis::openswath::scoring::{calculate_cross_correlation, standardize_data};
use openms::analysis::peak_integrator::{
    BaselineType, IntegrationMethod, PEAK_INTEGRATOR_NAME, PeakIntegrator,
};
use openms::concept::unique_id::UniqueIdGenerator;
use openms::kernel::mrm::{MRMTransitionGroup, SimpleTransition};
use openms::kernel::{
    ChromatogramPeak, DataArray, Feature, MSChromatogram, MSSpectrum, MobilityPeak1D, Mobilogram,
    Peak1D,
};
use openms::metadata::{MetaValue, MetaValueData};
use openms::param::{Param, ParamValue};
use openms::processing::chromatogram::{
    ChromatogramPickingMethod, ChromatogramSmoothing, PEAK_PICKER_CHROMATOGRAM_NAME,
    PeakPickerChromatogram,
};

fn rows(text: &str) -> impl Iterator<Item = Vec<&str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .skip(1)
        .map(|line| line.split('\t').collect())
}

/// One exact-bit trace of `tests/data/chromatogram_processing_traces.tsv`.
fn trace(name: &str) -> MSChromatogram {
    let peaks = rows(include_str!("data/chromatogram_processing_traces.tsv"))
        .filter(|row| row[0] == name)
        .map(|row| {
            let position: f64 = row[2].parse().unwrap();
            let intensity: f32 = row[3].parse().unwrap();
            assert_eq!(position.to_bits(), u64::from_str_radix(row[4], 16).unwrap());
            assert_eq!(
                intensity.to_bits(),
                u32::from_str_radix(row[5], 16).unwrap()
            );
            ChromatogramPeak::new(position, intensity)
        })
        .collect();
    MSChromatogram {
        peaks,
        ..Default::default()
    }
}

fn spectrum(chrom: &MSChromatogram) -> MSSpectrum {
    MSSpectrum::from_peaks(
        chrom
            .peaks
            .iter()
            .map(|p| Peak1D::new(p.rt, p.intensity))
            .collect(),
    )
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16}, expected={expected:.16}, tolerance={tolerance}"
    );
}

fn string(param: &Param, key: &str) -> String {
    param.value(key).unwrap().as_str().unwrap().to_owned()
}

fn set(param: &mut Param, key: &str, value: ParamValue) {
    param.set_value(key, value, "", &[]).unwrap();
}

// ---------------------------------------------------------------------------
// PeakPickerChromatogram parameters
// ---------------------------------------------------------------------------

/// `PeakPickerChromatogram.cpp:23-45`: names, order, types, values and
/// restrictions of the eleven source defaults.
#[test]
fn chromatogram_picker_defaults_follow_the_source_constructor() {
    let defaults = PeakPickerChromatogram::defaults().unwrap();
    let keys: Vec<_> = defaults
        .iter()
        .unwrap()
        .map(|item| item.key.clone())
        .collect();
    assert_eq!(
        keys,
        [
            "sgolay_frame_length",
            "sgolay_polynomial_order",
            "gauss_width",
            "use_gauss",
            "peak_width",
            "signal_to_noise",
            "sn_win_len",
            "sn_bin_count",
            "write_sn_log_messages",
            "report_sn",
            "remove_overlapping_peaks",
            "method",
        ]
    );
    assert_eq!(
        defaults.value("sgolay_frame_length").unwrap(),
        &ParamValue::Integer(15)
    );
    assert_eq!(
        defaults.value("sgolay_polynomial_order").unwrap(),
        &ParamValue::Integer(3)
    );
    assert_eq!(
        defaults.value("gauss_width").unwrap(),
        &ParamValue::Float(50.0)
    );
    assert_eq!(string(&defaults, "use_gauss"), "true");
    assert_eq!(
        defaults.value("peak_width").unwrap(),
        &ParamValue::Float(-1.0)
    );
    assert_eq!(
        defaults.value("signal_to_noise").unwrap(),
        &ParamValue::Float(1.0)
    );
    assert_eq!(
        defaults.value("sn_win_len").unwrap(),
        &ParamValue::Float(1000.0)
    );
    assert_eq!(
        defaults.value("sn_bin_count").unwrap(),
        &ParamValue::Integer(30)
    );
    assert_eq!(string(&defaults, "write_sn_log_messages"), "false");
    assert_eq!(string(&defaults, "report_sn"), "false");
    assert_eq!(string(&defaults, "remove_overlapping_peaks"), "false");
    assert_eq!(string(&defaults, "method"), "corrected");
    let entry = |key: &str| defaults.entry(key).unwrap().clone();
    assert_eq!(entry("use_gauss").valid_strings, ["false", "true"]);
    assert_eq!(
        entry("write_sn_log_messages").valid_strings,
        ["true", "false"]
    );
    assert_eq!(entry("report_sn").valid_strings, ["true", "false"]);
    assert_eq!(
        entry("remove_overlapping_peaks").valid_strings,
        ["false", "true"]
    );
    assert_eq!(
        entry("method").valid_strings,
        ["legacy", "corrected", "crawdad"]
    );
    assert_eq!(entry("signal_to_noise").min_float, 0.0);
    assert_eq!(
        defaults.description("sgolay_frame_length").unwrap(),
        "The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added."
    );
    for item in defaults.iter().unwrap() {
        assert!(!item.entry.tags.contains("advanced"), "{}", item.key);
    }

    // The defaults build the default picker.
    let picker = PeakPickerChromatogram::from_param(&defaults).unwrap();
    let native = PeakPickerChromatogram::default();
    assert_eq!(picker.method, native.method);
    assert_eq!(picker.smoothing, native.smoothing);
    assert_eq!(picker.inactive_smoothing, native.inactive_smoothing);
    assert_eq!(picker.peak_width, None);
    assert_eq!(picker.signal_to_noise, 1.0);
    assert_eq!(picker.seed_signal_to_noise, 1.0);
    assert_eq!(picker.noise_estimator.window_length, 1000.0);
    assert_eq!(picker.noise_estimator.bin_count, 30);
    assert!(!picker.noise_estimator.write_log_messages);
    assert!(!picker.report_sn && !picker.remove_overlapping_peaks);
    assert!(picker.to_param().unwrap().source_equal(&defaults).unwrap());
}

/// `setParameters` then `getParameters` returns every value, including the
/// parameters of the smoother `use_gauss` does not select.
#[test]
fn chromatogram_picker_parameters_round_trip() {
    let mut param = PeakPickerChromatogram::defaults().unwrap();
    set(&mut param, "sgolay_frame_length", ParamValue::Integer(9));
    set(
        &mut param,
        "sgolay_polynomial_order",
        ParamValue::Integer(2),
    );
    set(&mut param, "gauss_width", ParamValue::Float(30.0));
    set(&mut param, "use_gauss", ParamValue::String("false".into()));
    set(&mut param, "peak_width", ParamValue::Float(12.5));
    set(&mut param, "signal_to_noise", ParamValue::Float(0.5));
    set(&mut param, "sn_win_len", ParamValue::Float(500.0));
    set(&mut param, "sn_bin_count", ParamValue::Integer(12));
    set(
        &mut param,
        "write_sn_log_messages",
        ParamValue::String("true".into()),
    );
    set(&mut param, "report_sn", ParamValue::String("true".into()));
    set(
        &mut param,
        "remove_overlapping_peaks",
        ParamValue::String("true".into()),
    );
    set(&mut param, "method", ParamValue::String("legacy".into()));
    let (picker, warnings) = PeakPickerChromatogram::from_param_with_warnings(&param).unwrap();
    assert!(warnings.is_empty(), "{warnings:?}");
    assert_eq!(picker.method, ChromatogramPickingMethod::Legacy);
    assert_eq!(
        picker.smoothing,
        ChromatogramSmoothing::SavitzkyGolay {
            frame_length: 9,
            polynomial_order: 2
        }
    );
    assert_eq!(
        picker.inactive_smoothing,
        ChromatogramSmoothing::Gaussian { width: 30.0 }
    );
    assert_eq!(picker.peak_width, Some(12.5));
    assert_eq!(picker.signal_to_noise, 0.5);
    assert_eq!(picker.noise_estimator.window_length, 500.0);
    assert_eq!(picker.noise_estimator.bin_count, 12);
    assert!(picker.noise_estimator.write_log_messages);
    assert!(picker.report_sn && picker.remove_overlapping_peaks);
    // The seed picker is configured once, in the source constructor.
    assert_eq!(picker.seed_signal_to_noise, 1.0);
    let written = picker.to_param().unwrap();
    for item in param.iter().unwrap() {
        assert_eq!(
            written.value(&item.key).unwrap(),
            &item.entry.value,
            "{}",
            item.key
        );
    }
}

#[test]
fn chromatogram_picker_parameter_refusals_follow_the_source() {
    let with = |key: &str, value: ParamValue| {
        let mut param = PeakPickerChromatogram::defaults().unwrap();
        set(&mut param, key, value);
        PeakPickerChromatogram::from_param(&param)
    };
    // The default source build has no Crawdad (`WITH_CRAWDAD` off), and
    // `updateMembers_` throws `IllegalArgument` with this text.
    let error = with("method", ParamValue::String("crawdad".into()))
        .unwrap_err()
        .to_string();
    assert!(
        error.contains(
            "PeakPickerChromatogram was not compiled with crawdad, please choose a different algorithm!"
        ),
        "{error}"
    );
    assert!(with("method", ParamValue::String("other".into())).is_err());
    assert!(with("signal_to_noise", ParamValue::Float(-0.5)).is_err());
    // `snt_.setParameters`: bin_count at least three, win_len at least one.
    assert!(with("sn_bin_count", ParamValue::Integer(2)).is_err());
    assert!(with("sn_win_len", ParamValue::Float(0.5)).is_err());
    // `sgolay_.setParameters`: the order must be below the frame length, and
    // it is configured even when the Gaussian filter is selected.
    assert!(with("sgolay_polynomial_order", ParamValue::Integer(15)).is_err());
    assert!(with("use_gauss", ParamValue::String("maybe".into())).is_err());
    assert!(with("gauss_width", ParamValue::Integer(50)).is_err());
    // Every peak_width not above zero disables the forced width.
    for disabled in [0.0, -1.0, -7.0, f64::NAN] {
        assert_eq!(
            with("peak_width", ParamValue::Float(disabled))
                .unwrap()
                .peak_width,
            None
        );
    }
    // An unknown parameter is a warning, as in the source.
    let mut param = PeakPickerChromatogram::defaults().unwrap();
    set(&mut param, "not_a_parameter", ParamValue::Integer(1));
    let (_, warnings) = PeakPickerChromatogram::from_param_with_warnings(&param).unwrap();
    assert_eq!(warnings.len(), 1);
    assert!(
        warnings[0].contains(PEAK_PICKER_CHROMATOGRAM_NAME),
        "{warnings:?}"
    );
}

/// `PeakPickerChromatogram_test.cpp:75-137`: the class test configures the
/// picker through `getDefaults`, `setValue("method" / "peak_width")` and
/// `setParameters`, then picks both SRM traces.
#[test]
fn chromatogram_picker_class_test_through_parameters() {
    for row in rows(include_str!("data/chromatogram_processing_picked.tsv")) {
        let mut param = PeakPickerChromatogram::defaults().unwrap();
        set(&mut param, "method", ParamValue::String(row[1].into()));
        set(
            &mut param,
            "peak_width",
            ParamValue::Float(row[2].parse().unwrap()),
        );
        let picker = PeakPickerChromatogram::from_param(&param).unwrap();
        let input = trace(row[0]);
        let result = picker.pick_chromatogram(&input).unwrap();
        let picked = &result.picked.chromatogram;
        assert_eq!(picked.len(), 1);
        assert_eq!(
            picked.float_data_arrays.len(),
            PeakPickerChromatogram::SIZE_OF_FLOATINDICES
        );
        let names: Vec<_> = picked
            .float_data_arrays
            .iter()
            .map(|a| a.name.as_str())
            .collect();
        assert_eq!(
            names,
            [
                "FWHM",
                "IntegratedIntensity",
                "leftWidth",
                "rightWidth",
                "SN"
            ]
        );
        let array = |index: usize| picked.float_data_arrays[index].data[0];
        close(picked.peaks[0].rt, row[3].parse().unwrap(), 0.001);
        close(
            f64::from(picked.peaks[0].intensity),
            row[4].parse().unwrap(),
            0.02,
        );
        close(
            f64::from(array(PeakPickerChromatogram::IDX_ABUNDANCE)),
            row[5].parse().unwrap(),
            0.6,
        );
        assert_eq!(
            array(PeakPickerChromatogram::IDX_LEFTBORDER),
            row[6].parse::<f32>().unwrap()
        );
        assert_eq!(
            array(PeakPickerChromatogram::IDX_RIGHTBORDER),
            row[7].parse::<f32>().unwrap()
        );
        // `signal_to_noise` is 1 > 0, so the SN array holds estimates, not -1.
        assert_ne!(array(PeakPickerChromatogram::IDX_SN), -1.0);
        assert!(result.noise_log.is_empty());
    }
}

/// `write_sn_log_messages` reaches the boundary noise estimator, whose
/// warnings the result carries instead of the source's log stream.
#[test]
fn chromatogram_picker_returns_noise_warnings_when_asked() {
    // Samples about 2.8 s apart in a 3 s window leave every window below the
    // estimator's ten-point minimum, which the estimator reports as sparse.
    let input = trace("srm_2");
    let mut param = PeakPickerChromatogram::defaults().unwrap();
    set(&mut param, "sn_win_len", ParamValue::Float(3.0));
    let quiet = PeakPickerChromatogram::from_param(&param)
        .unwrap()
        .pick_chromatogram(&input)
        .unwrap();
    assert!(quiet.noise_log.is_empty());
    set(
        &mut param,
        "write_sn_log_messages",
        ParamValue::String("true".into()),
    );
    let loud = PeakPickerChromatogram::from_param(&param)
        .unwrap()
        .pick_chromatogram(&input)
        .unwrap();
    assert!(!loud.noise_log.is_empty());
    assert_eq!(loud.picked, quiet.picked);
}

// ---------------------------------------------------------------------------
// PeakIntegrator parameters and index-bounded overloads
// ---------------------------------------------------------------------------

/// `PeakIntegrator.cpp:82-94` and `PeakIntegrator_test.cpp:165-171`.
#[test]
fn integrator_defaults_follow_get_default_parameters() {
    let defaults = PeakIntegrator::defaults().unwrap();
    let keys: Vec<_> = defaults
        .iter()
        .unwrap()
        .map(|item| item.key.clone())
        .collect();
    assert_eq!(keys, ["integration_type", "baseline_type", "fit_EMG"]);
    assert_eq!(
        string(&defaults, "integration_type"),
        PeakIntegrator::INTEGRATION_TYPE_INTENSITYSUM
    );
    assert_eq!(
        string(&defaults, "baseline_type"),
        PeakIntegrator::BASELINE_TYPE_BASETOBASE
    );
    assert_eq!(string(&defaults, "fit_EMG"), "false");
    assert_eq!(
        defaults.entry("integration_type").unwrap().valid_strings,
        ["intensity_sum", "simpson", "trapezoid"]
    );
    assert_eq!(
        defaults.entry("baseline_type").unwrap().valid_strings,
        [
            "base_to_base",
            "vertical_division",
            "vertical_division_min",
            "vertical_division_max"
        ]
    );
    assert_eq!(
        defaults.entry("fit_EMG").unwrap().valid_strings,
        ["false", "true"]
    );
    assert_eq!(
        PeakIntegrator::from_param(&defaults).unwrap(),
        PeakIntegrator::default()
    );
    assert_eq!(
        [
            PeakIntegrator::INTEGRATION_TYPE_TRAPEZOID,
            PeakIntegrator::INTEGRATION_TYPE_SIMPSON,
            PeakIntegrator::BASELINE_TYPE_VERTICALDIVISION,
            PeakIntegrator::BASELINE_TYPE_VERTICALDIVISION_MIN,
            PeakIntegrator::BASELINE_TYPE_VERTICALDIVISION_MAX,
        ],
        [
            "trapezoid",
            "simpson",
            "vertical_division",
            "vertical_division_min",
            "vertical_division_max"
        ]
    );
}

#[test]
fn integrator_parameters_round_trip_and_refuse_unknown_types() {
    let mut param = PeakIntegrator::defaults().unwrap();
    set(
        &mut param,
        "integration_type",
        ParamValue::String("simpson".into()),
    );
    set(
        &mut param,
        "baseline_type",
        ParamValue::String("vertical_division".into()),
    );
    set(&mut param, "fit_EMG", ParamValue::String("true".into()));
    let integrator = PeakIntegrator::from_param(&param).unwrap();
    assert_eq!(integrator.integration_method, IntegrationMethod::Simpson);
    assert_eq!(integrator.baseline_type, BaselineType::VerticalDivisionMin);
    assert!(integrator.emg.is_some());
    let written = integrator.to_param().unwrap();
    assert_eq!(string(&written, "integration_type"), "simpson");
    // The source alias is written back under the name it stands for.
    assert_eq!(string(&written, "baseline_type"), "vertical_division_min");
    assert_eq!(string(&written, "fit_EMG"), "true");

    let mut bad = PeakIntegrator::defaults().unwrap();
    set(
        &mut bad,
        "integration_type",
        ParamValue::String("midpoint".into()),
    );
    assert!(PeakIntegrator::from_param(&bad).is_err());
    let mut unknown = PeakIntegrator::defaults().unwrap();
    set(&mut unknown, "surplus", ParamValue::Integer(1));
    let (_, warnings) = PeakIntegrator::from_param_with_warnings(&unknown).unwrap();
    assert_eq!(warnings.len(), 1);
    assert!(warnings[0].contains(PEAK_INTEGRATOR_NAME), "{warnings:?}");
}

/// The source iterators of `PeakIntegrator_test.cpp:111-118`: `RTBegin(left)`
/// and `RTEnd(right) - 1`, as sample indices.
fn begin(input: &MSChromatogram, position: f64) -> usize {
    input.peaks.partition_point(|p| p.rt < position)
}
fn last(input: &MSChromatogram, position: f64) -> usize {
    input.peaks.partition_point(|p| p.rt <= position) - 1
}

/// `PeakIntegrator_test.cpp:506-586` and the spectrum twin at `669-749`.
#[test]
fn integrator_index_overloads_reproduce_the_class_test() {
    let input = trace("glutamate");
    let spec = spectrum(&input);
    let left = begin(&input, 2.472833334);
    let right = last(&input, 3.022891666);
    let cases = [
        // (method, right index, area, height, apex)
        (
            IntegrationMethod::IntensitySum,
            right,
            6768778.0,
            966489.0,
            2.7045,
        ),
        (
            IntegrationMethod::Trapezoid,
            right,
            71540.2,
            966489.0,
            2.7045,
        ),
        (
            IntegrationMethod::Trapezoid,
            last(&input, 2.477966667),
            0.0,
            881.0,
            2.477966667,
        ),
        (
            IntegrationMethod::Simpson,
            right,
            71720.443144994,
            966489.0,
            2.7045,
        ),
        (
            IntegrationMethod::Simpson,
            right - 1,
            71515.0792609335,
            966489.0,
            2.7045,
        ),
        (
            IntegrationMethod::Simpson,
            last(&input, 2.477966667),
            0.0,
            881.0,
            2.477966667,
        ),
        (
            IntegrationMethod::Simpson,
            last(&input, 2.488216667),
            11.6081250000001,
            1384.0,
            2.488216667,
        ),
    ];
    for (method, right, area, height, apex) in cases {
        let integrator = PeakIntegrator {
            integration_method: method,
            ..Default::default()
        };
        let chrom = integrator
            .integrate_chromatogram_between(&input, left, right)
            .unwrap();
        // The six-significant-digit trapezoid literal, as in the neighbouring
        // position-bounded test.
        let tolerance = if area == 71540.2 { 0.01 } else { 1e-6 };
        close(chrom.area, area, tolerance);
        close(chrom.height, height, 0.0);
        close(chrom.apex_pos, apex, 0.0);
        assert_eq!(chrom.hull_points.len(), right - left + 1);
        for (point, sample) in chrom.hull_points.iter().zip(&input.peaks[left..=right]) {
            assert_eq!(*point, [sample.rt, f64::from(sample.intensity)]);
        }
        assert_eq!(
            integrator
                .integrate_spectrum_between(&spec, left, right)
                .unwrap(),
            chrom
        );
        assert_eq!(
            integrator
                .integrate_chromatogram(&input, input.peaks[left].rt, input.peaks[right].rt)
                .unwrap(),
            chrom
        );
    }
}

/// `PeakIntegrator_test.cpp:232-289` and the spectrum twin at `350-407`.
#[test]
fn integrator_index_overloads_reproduce_background_class_test() {
    let input = trace("glutamate");
    let spec = spectrum(&input);
    let left = begin(&input, 2.472833334);
    let right = last(&input, 3.022891666);
    for row in rows(include_str!("data/chromatogram_processing_backgrounds.tsv")) {
        let integrator = PeakIntegrator {
            integration_method: row[0].parse().unwrap(),
            baseline_type: row[1].parse().unwrap(),
            ..Default::default()
        };
        let apex = integrator
            .integrate_chromatogram_between(&input, left, right)
            .unwrap()
            .apex_pos;
        let chrom = integrator
            .estimate_background_chromatogram_between(&input, left, right, apex)
            .unwrap();
        close(chrom.area, row[2].parse().unwrap(), 1e-6);
        close(chrom.height, row[3].parse().unwrap(), 1e-8);
        assert_eq!(
            integrator
                .estimate_background_spectrum_between(&spec, left, right, apex)
                .unwrap(),
            chrom
        );
    }
}

/// `PeakIntegrator_test.cpp:801-849` and the spectrum twin at `901-949`.
#[test]
fn integrator_index_overloads_reproduce_shape_class_test() {
    let input = trace("glutamate");
    let spec = spectrum(&input);
    let integrator = PeakIntegrator::default();
    let left = begin(&input, 2.472833334);
    let right = last(&input, 3.022891666);
    let shape = integrator
        .calculate_shape_metrics_chromatogram_between(&input, left, right, 966489.0, 2.7045)
        .unwrap();
    for row in rows(include_str!("data/chromatogram_processing_shape.tsv")) {
        let value = match row[0] {
            "width_at_5" => shape.width_at_5,
            "width_at_10" => shape.width_at_10,
            "width_at_50" => shape.width_at_50,
            "start_position_at_5" => shape.start_position_at_5,
            "start_position_at_10" => shape.start_position_at_10,
            "start_position_at_50" => shape.start_position_at_50,
            "end_position_at_5" => shape.end_position_at_5,
            "end_position_at_10" => shape.end_position_at_10,
            "end_position_at_50" => shape.end_position_at_50,
            "total_width" => shape.total_width,
            "tailing_factor" => shape.tailing_factor,
            "asymmetry_factor" => shape.asymmetry_factor,
            "slope_of_baseline" => shape.slope_of_baseline,
            "baseline_delta_2_height" => shape.baseline_delta_2_height,
            "points_across_baseline" => shape.points_across_baseline as f64,
            "points_across_half_height" => shape.points_across_half_height as f64,
            other => panic!("unknown shape field {other}"),
        };
        close(value, row[1].parse().unwrap(), 1e-10);
    }
    assert_eq!(
        integrator
            .calculate_shape_metrics_spectrum_between(&spec, left, right, 966489.0, 2.7045)
            .unwrap(),
        shape
    );
    // The source dereferences the iterators unchecked; an index that names no
    // sample is refused here.
    let past = input.len();
    assert!(
        integrator
            .integrate_chromatogram_between(&input, left, past)
            .is_err()
    );
    assert!(
        integrator
            .integrate_spectrum_between(&spec, past, right)
            .is_err()
    );
    assert!(
        integrator
            .estimate_background_chromatogram_between(&input, left, past, 2.7045)
            .is_err()
    );
    assert!(
        integrator
            .calculate_shape_metrics_chromatogram_between(
                &MSChromatogram::default(),
                0,
                0,
                1.0,
                0.0
            )
            .is_err()
    );
}

// ---------------------------------------------------------------------------
// MRMTransitionGroupPicker
// ---------------------------------------------------------------------------

/// `TEST_REAL_SIMILAR` at its defaults (`ClassTest.cpp:35-38`): absolute
/// difference within 1e-5, or ratio within 1 + 1e-5.
fn real_similar(actual: f64, expected: f64) {
    let absdiff = (actual - expected).abs();
    let ratio = if actual.abs() > expected.abs() {
        actual / expected
    } else {
        expected / actual
    };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{actual} !~ {expected} (absdiff {absdiff}, ratio {ratio})"
    );
}

type Group = MRMTransitionGroup<MSChromatogram, SimpleTransition>;

fn srm_chromatogram(id: &str, rts: &[f64], intensities: &[f64]) -> MSChromatogram {
    MSChromatogram {
        native_id: id.into(),
        peaks: rts
            .iter()
            .zip(intensities)
            .map(|(&rt, &intensity)| ChromatogramPeak::new(rt, intensity as f32))
            .collect(),
        ..Default::default()
    }
}

fn with_product_mz(mut chromatogram: MSChromatogram, mz: f64) -> MSChromatogram {
    chromatogram
        .metadata
        .insert("product_mz".into(), MetaValue::try_from(mz).unwrap());
    chromatogram
}

/// `setup_transition_group` of `MRMTransitionGroupPicker_test.cpp:34-117`.
fn transition_group_1() -> Group {
    let rt_1 = [
        1474.34, 1477.11, 1479.88, 1482.64, 1485.41, 1488.19, 1490.95, 1493.72, 1496.48, 1499.25,
        1502.03, 1504.8, 1507.56, 1510.33, 1513.09, 1515.87, 1518.64, 1521.42,
    ];
    let rt_2 = [
        1473.55, 1476.31, 1479.08, 1481.84, 1484.61, 1487.39, 1490.15, 1492.92, 1495.69, 1498.45,
        1501.23, 1504.0, 1506.76, 1509.53, 1512.29, 1515.07, 1517.84, 1520.62,
    ];
    let int_1 = [
        3.26958, 3.74189, 3.31075, 86.1901, 3.47528, 387.864, 13281.0, 6375.84, 39852.6, 2.66726,
        612.747, 3.34313, 793.12, 3.29156, 4.00586, 4.1591, 3.23035, 3.90591,
    ];
    let int_2 = [
        3.44054, 2142.31, 3.58763, 3076.97, 6663.55, 45681.0, 157694.0, 122844.0, 86034.7, 85391.1,
        15992.8, 2293.94, 6934.85, 2735.18, 459.413, 3.93863, 3.36564, 3.44005,
    ];
    let ms1 = [
        900.0, 2600.0, 4100.0, 5400.0, 6500.0, 7400.0, 8100.0, 8600.0, 8900.0, 9000.0, 8900.0,
        8600.0, 8100.0, 7400.0, 6500.0, 5400.0, 4100.0, 2600.0,
    ];
    let mut group = Group::new();
    for (id, rts, intensities, mz) in [("1", &rt_1, &int_1, 618.31), ("2", &rt_2, &int_2, 619.31)] {
        group
            .add_chromatogram(
                with_product_mz(srm_chromatogram(id, rts, intensities), mz),
                id,
            )
            .unwrap();
        group
            .add_transition(SimpleTransition::new(id, -1.0), id)
            .unwrap();
    }
    let shifted: Vec<f64> = rt_2.iter().map(|rt| rt + 0.5).collect();
    group
        .add_precursor_chromatogram(
            srm_chromatogram("Precursor_i0", &shifted, &ms1),
            "Precursor_i0",
        )
        .unwrap();
    group
}

/// `setup_transition_group2` of `MRMTransitionGroupPicker_test.cpp:119-186`.
fn transition_group_2() -> Group {
    let rts: Vec<f64> = (0..30).map(f64::from).collect();
    let traces: [(&str, [f64; 30], f64); 3] = [
        (
            "1",
            [
                3., 3., 3., 3., 3., 3., 3., 3., 4., 6., 8., 12., 14., 15., 16., 15., 14., 12., 8.,
                6., 4., 3., 3., 3., 3., 3., 3., 3., 3., 3.,
            ],
            618.31,
        ),
        (
            "2",
            [
                3., 3., 4., 6., 8., 10., 11., 11., 10., 8., 6., 4., 4., 6., 8., 10., 11., 11., 10.,
                8., 6., 4., 3., 3., 3., 3., 3., 3., 3., 3.,
            ],
            619.31,
        ),
        ("3", [3.0; 30], 812.2),
    ];
    let mut group = Group::new();
    for (id, intensities, mz) in traces {
        group
            .add_chromatogram(
                with_product_mz(srm_chromatogram(id, &rts, &intensities), mz),
                id,
            )
            .unwrap();
        group
            .add_transition(SimpleTransition::new(id, -1.0), id)
            .unwrap();
    }
    group
}

fn picker_with(settings: &[(&str, ParamValue)]) -> MRMTransitionGroupPicker {
    let mut param = MRMTransitionGroupPicker::defaults().unwrap();
    for (key, value) in settings {
        set(&mut param, key, value.clone());
    }
    let (picker, warnings) = MRMTransitionGroupPicker::from_param_with_warnings(&param).unwrap();
    assert!(warnings.is_empty(), "{warnings:?}");
    picker
}

fn legacy_picker() -> MRMTransitionGroupPicker {
    picker_with(&[
        (
            "PeakPickerChromatogram:method",
            ParamValue::String("legacy".into()),
        ),
        ("PeakPickerChromatogram:peak_width", ParamValue::Float(40.0)),
    ])
}

fn corrected_picker(use_consensus: bool) -> MRMTransitionGroupPicker {
    picker_with(&[
        ("resample_boundary", ParamValue::Float(1.0)),
        (
            "PeakPickerChromatogram:gauss_width",
            ParamValue::Float(10.0),
        ),
        ("PeakPickerChromatogram:peak_width", ParamValue::Float(-1.0)),
        (
            "PeakPickerChromatogram:signal_to_noise",
            ParamValue::Float(1.0),
        ),
        (
            "PeakPickerChromatogram:method",
            ParamValue::String("corrected".into()),
        ),
        (
            "use_consensus",
            ParamValue::String(if use_consensus { "true" } else { "false" }.into()),
        ),
    ])
}

fn meta(feature: &Feature, key: &str) -> f64 {
    feature.base.metadata[key].as_f64().unwrap()
}

fn hull_sum(feature: &Feature) -> f64 {
    feature.convex_hulls[0]
        .hull_points()
        .iter()
        .fold(0.0, |sum, point| sum + point.mz)
}

fn ids() -> UniqueIdGenerator {
    UniqueIdGenerator::from_seed(1)
}

/// `MRMTransitionGroupPicker.cpp:27-65`.
#[test]
fn group_picker_defaults_follow_the_source_constructor() {
    let defaults = MRMTransitionGroupPicker::defaults().unwrap();
    let top: Vec<_> = defaults
        .iter()
        .unwrap()
        .map(|item| item.key.clone())
        .filter(|key| !key.contains(':'))
        .collect();
    assert_eq!(
        top,
        [
            "stop_after_feature",
            "stop_after_intensity_ratio",
            "min_peak_width",
            "peak_integration",
            "background_subtraction",
            "recalculate_peaks",
            "use_precursors",
            "use_consensus",
            "recalculate_peaks_max_z",
            "minimal_quality",
            "resample_boundary",
            "compute_peak_quality",
            "compute_peak_shape_metrics",
            "compute_total_mi",
            "boundary_selection_method",
        ]
    );
    for key in &top[2..] {
        assert!(
            defaults.entry(key).unwrap().tags.contains("advanced"),
            "{key}"
        );
    }
    for key in &top[..2] {
        assert!(
            !defaults.entry(key).unwrap().tags.contains("advanced"),
            "{key}"
        );
    }
    assert_eq!(
        defaults.value("stop_after_feature").unwrap(),
        &ParamValue::Integer(-1)
    );
    assert_eq!(
        defaults.value("stop_after_intensity_ratio").unwrap(),
        &ParamValue::Float(0.0001)
    );
    assert_eq!(
        defaults.value("min_peak_width").unwrap(),
        &ParamValue::Float(0.001)
    );
    assert_eq!(string(&defaults, "peak_integration"), "original");
    assert_eq!(string(&defaults, "background_subtraction"), "none");
    assert_eq!(string(&defaults, "use_consensus"), "true");
    assert_eq!(
        defaults.value("recalculate_peaks_max_z").unwrap(),
        &ParamValue::Float(1.0)
    );
    assert_eq!(
        defaults.value("minimal_quality").unwrap(),
        &ParamValue::Float(-10000.0)
    );
    assert_eq!(
        defaults.value("resample_boundary").unwrap(),
        &ParamValue::Float(15.0)
    );
    assert_eq!(string(&defaults, "boundary_selection_method"), "largest");
    assert_eq!(
        defaults
            .entry("background_subtraction")
            .unwrap()
            .valid_strings,
        ["none", "original", "exact"]
    );
    // The two subsections are the nested classes' own defaults.
    assert!(
        defaults
            .copy("PeakPickerChromatogram:", true)
            .unwrap()
            .source_equal(&PeakPickerChromatogram::defaults().unwrap())
            .unwrap()
    );
    assert!(
        defaults
            .copy("PeakIntegrator:", true)
            .unwrap()
            .source_equal(&PeakIntegrator::defaults().unwrap())
            .unwrap()
    );
    let picker = MRMTransitionGroupPicker::from_param(&defaults).unwrap();
    assert!(picker.use_consensus && !picker.use_precursors);
    assert_eq!(picker.stop_after_feature, -1);
    assert_eq!(picker.boundary_selection_method, BoundarySelection::Largest);
    assert_eq!(picker.peak_integration, PeakIntegrationSignal::Original);
    assert_eq!(picker.background_subtraction, BackgroundSubtraction::None);
}

#[test]
fn group_picker_parameters_reach_the_nested_classes() {
    let picker = picker_with(&[
        (
            "PeakPickerChromatogram:method",
            ParamValue::String("legacy".into()),
        ),
        (
            "PeakIntegrator:integration_type",
            ParamValue::String("simpson".into()),
        ),
        ("background_subtraction", ParamValue::String("exact".into())),
        (
            "boundary_selection_method",
            ParamValue::String("widest".into()),
        ),
        ("stop_after_feature", ParamValue::Integer(3)),
    ]);
    assert_eq!(picker.picker.method, ChromatogramPickingMethod::Legacy);
    assert_eq!(
        picker.integrator.integration_method,
        IntegrationMethod::Simpson
    );
    assert_eq!(picker.background_subtraction, BackgroundSubtraction::Exact);
    assert_eq!(picker.boundary_selection_method, BoundarySelection::Widest);
    assert_eq!(picker.stop_after_feature, 3);
    let written = picker.to_param().unwrap();
    assert_eq!(string(&written, "PeakPickerChromatogram:method"), "legacy");
    assert_eq!(
        string(&written, "PeakIntegrator:integration_type"),
        "simpson"
    );
    assert_eq!(string(&written, "boundary_selection_method"), "widest");

    // An unknown subsection key warns twice, as the source's two
    // `setParameters` calls do; a bad value is refused.
    let mut param = MRMTransitionGroupPicker::defaults().unwrap();
    set(&mut param, "PeakIntegrator:surplus", ParamValue::Integer(1));
    let (_, warnings) = MRMTransitionGroupPicker::from_param_with_warnings(&param).unwrap();
    assert_eq!(warnings.len(), 2, "{warnings:?}");
    let mut bad = MRMTransitionGroupPicker::defaults().unwrap();
    set(
        &mut bad,
        "peak_integration",
        ParamValue::String("both".into()),
    );
    assert!(MRMTransitionGroupPicker::from_param(&bad).is_err());
}

/// `MRMTransitionGroupPicker_test.cpp:249-311`.
#[test]
fn group_picker_class_test_transition_group_1() {
    let mut group = transition_group_1();
    legacy_picker()
        .pick_transition_group(&mut group, &mut ids())
        .unwrap();
    assert_eq!(group.features().len(), 1);
    let feature = &group.features()[0];
    real_similar(feature.base.rt, 1492.83060);
    real_similar(f64::from(feature.base.intensity), 564236.0);
    real_similar(meta(feature, "leftWidth"), 1484.60998535156);
    real_similar(meta(feature, "rightWidth"), 1501.23);
    assert_eq!(
        feature.feature("1").unwrap().convex_hulls[0]
            .hull_points()
            .len(),
        6
    );
    assert_eq!(
        feature.feature("2").unwrap().convex_hulls[0]
            .hull_points()
            .len(),
        6
    );
    real_similar(hull_sum(feature.feature("2").unwrap()), 504308.354492188);
    real_similar(
        f64::from(feature.feature("2").unwrap().base.intensity),
        504308.0,
    );
    real_similar(hull_sum(feature.feature("1").unwrap()), 59927.5729999542);
    real_similar(
        f64::from(feature.feature("1").unwrap().base.intensity),
        59927.6,
    );
    assert_eq!(
        feature.precursor_feature_ids().collect::<Vec<_>>(),
        ["Precursor_i0"]
    );
    let resampling_loss = 6043.9760742188;
    let precursor = feature.precursor_feature("Precursor_i0").unwrap();
    real_similar(hull_sum(precursor), 53900.0 - resampling_loss);
    real_similar(
        f64::from(precursor.base.intensity),
        53900.0 - resampling_loss,
    );
    // Native: the unique ID is drawn from the caller's generator.
    assert_ne!(feature.base.unique_id, 0);
}

/// `MRMTransitionGroupPicker_test.cpp:313-346`.
#[test]
fn group_picker_class_test_transition_group_2_consensus() {
    let mut group = transition_group_2();
    corrected_picker(true)
        .pick_transition_group(&mut group, &mut ids())
        .unwrap();
    assert_eq!(group.features().len(), 1);
    let feature = &group.features()[0];
    real_similar(feature.base.rt, 14.0);
    real_similar(meta(feature, "leftWidth"), 8.0);
    real_similar(meta(feature, "rightWidth"), 20.0);
    for (id, intensity, apex_int, apex_pos) in [
        ("1", 134.0, 16.0, 14.0),
        ("2", 102.0, 11.0, 16.0),
        ("3", 39.0, 3.0, 8.0),
    ] {
        let sub = feature.feature(id).unwrap();
        real_similar(f64::from(sub.base.intensity), intensity);
        real_similar(meta(sub, "peak_apex_int"), apex_int);
        real_similar(meta(sub, "peak_apex_position"), apex_pos);
    }
}

/// `MRMTransitionGroupPicker_test.cpp:348-377` and `432-498`.
#[test]
fn group_picker_class_test_transition_group_2_without_consensus() {
    let mut group = transition_group_2();
    corrected_picker(false)
        .pick_transition_group(&mut group, &mut ids())
        .unwrap();
    assert_eq!(group.features().len(), 2);
    let feature = &group.features()[0];
    real_similar(f64::from(feature.base.intensity), 237.0);
    real_similar(feature.base.rt, 14.0);
    real_similar(meta(feature, "leftWidth"), 8.0);
    real_similar(meta(feature, "rightWidth"), 20.0);
    for (id, intensity, apex_int, apex_pos) in [
        ("1", 134.0, 16.0, 14.0),
        ("2", 64.0, 11.0, 16.0),
        ("3", 39.0, 3.0, 8.0),
    ] {
        let sub = feature.feature(id).unwrap();
        real_similar(f64::from(sub.base.intensity), intensity);
        real_similar(meta(sub, "peak_apex_int"), apex_int);
        real_similar(meta(sub, "peak_apex_position"), apex_pos);
    }
    // Source issue 9138: equally long hulls outside consensus mode.
    let hulls: Vec<_> = ["1", "2", "3"]
        .iter()
        .map(|id| {
            feature.feature(id).unwrap().convex_hulls[0]
                .hull_points()
                .len()
        })
        .collect();
    assert!(hulls[0] > 0);
    assert_eq!(hulls, [hulls[0]; 3]);

    // A non-detecting transition is an error outside consensus mode, and the
    // group is left unchanged.
    let mut group = transition_group_2();
    group.transitions_mut()[0].detecting = false;
    assert!(
        corrected_picker(false)
            .pick_transition_group(&mut group, &mut ids())
            .is_err()
    );
    assert!(group.features().is_empty());
}

/// `MRMTransitionGroupPicker_test.cpp:379-430`.
#[test]
fn group_picker_class_test_only_quantifying_transitions_count() {
    let picker = legacy_picker();
    for (quantifying, expected) in [
        ([true, true], Some(564236.0)),
        ([false, true], Some(504308.0)),
        ([true, false], Some(59927.6)),
        ([false, false], None),
    ] {
        let mut group = transition_group_1();
        for (transition, flag) in group.transitions_mut().iter_mut().zip(quantifying) {
            transition.quantifying = flag;
        }
        picker
            .pick_transition_group(&mut group, &mut ids())
            .unwrap();
        match expected {
            Some(intensity) => {
                assert_eq!(group.features().len(), 1);
                real_similar(f64::from(group.features()[0].base.intensity), intensity);
                real_similar(group.features()[0].base.rt, 1492.83060);
            }
            // Zero-intensity peak groups are dropped.
            None => assert!(group.features().is_empty()),
        }
    }
}

fn empty_arrays() -> Vec<DataArray<f32>> {
    (0..PeakPickerChromatogram::SIZE_OF_FLOATINDICES)
        .map(|_| DataArray::new("", Vec::new()))
        .collect()
}

/// The hand-made picked chromatograms of
/// `MRMTransitionGroupPicker_test.cpp:506-534`: one peak at 1490 s with
/// boundaries 1481.84 and 1512.29 in each fragment chromatogram.
fn hand_picked(group: &Group) -> Vec<MSChromatogram> {
    group
        .chromatograms()
        .iter()
        .map(|chromatogram| {
            let mut arrays = empty_arrays();
            arrays[PeakPickerChromatogram::IDX_ABUNDANCE] =
                DataArray::new("IntegratedIntensity", vec![1000.0]);
            arrays[PeakPickerChromatogram::IDX_LEFTBORDER] =
                DataArray::new("leftWidth", vec![1481.840_f64 as f32]);
            arrays[PeakPickerChromatogram::IDX_RIGHTBORDER] =
                DataArray::new("rightWidth", vec![1512.290_f64 as f32]);
            MSChromatogram {
                native_id: chromatogram.native_id.clone(),
                peaks: vec![ChromatogramPeak::new(1490.0, 170.0)],
                float_data_arrays: arrays,
                ..Default::default()
            }
        })
        .collect()
}

/// Settings and expected `[area 2, apex 2, area 1, apex 1, ratio]` of one
/// `createMRMFeature` class-test block.
type CreateCase<'a> = (Vec<(&'a str, ParamValue)>, [f64; 5]);

/// `MRMTransitionGroupPicker_test.cpp:504-634`. The class test spells its
/// picker keys `PeakPickerChromatogram::method`, which match no parameter; the
/// chromatogram picker is not used by `createMRMFeature` either way.
#[test]
fn group_picker_class_test_create_mrm_feature() {
    let group = transition_group_1();
    // (settings, [area 2, apex 2, area 1, apex 1 (NaN: not asserted), ratio])
    let cases: [CreateCase<'_>; 4] = [
        (
            Vec::new(),
            [535801.503, 157694.0, 61405.95106, 30286.9130513267, 8.72556],
        ),
        (
            vec![(
                "background_subtraction",
                ParamValue::String("original".into()),
            )],
            [514583.0, 155925.8085, 60913.8, f64::NAN, 8.44773],
        ),
        (
            vec![
                ("background_subtraction", ParamValue::String("exact".into())),
                (
                    "PeakIntegrator:integration_type",
                    ParamValue::String("simpson".into()),
                ),
            ],
            [
                1.42194e+06,
                155331.37806798,
                168685.0,
                30251.2414481247,
                8.42957,
            ],
        ),
        (
            vec![("background_subtraction", ParamValue::String("exact".into()))],
            [
                514590.0,
                155331.37806798,
                61009.7,
                30251.2414481247,
                8.43456040596823,
            ],
        ),
    ];
    for (index, (settings, [area_2, apex_2, area_1, apex_1, ratio])) in cases.iter().enumerate() {
        let picker = picker_with(settings);
        let mut picked = hand_picked(&group);
        let feature = picker
            .create_mrm_feature(&group, &mut picked, &[], 1, 0, &mut ids())
            .unwrap();
        let one = feature.feature("1").unwrap();
        let two = feature.feature("2").unwrap();
        real_similar(f64::from(two.base.intensity), *area_2);
        real_similar(meta(two, "peak_apex_int"), *apex_2);
        real_similar(f64::from(one.base.intensity), *area_1);
        if !apex_1.is_nan() {
            real_similar(meta(one, "peak_apex_int"), *apex_1);
        }
        real_similar(
            f64::from(two.base.intensity) / f64::from(one.base.intensity),
            *ratio,
        );
        if index == 0 {
            assert_eq!(feature.base.rt, 1490.0);
            assert_eq!(one.convex_hulls[0].hull_points().len(), 12);
            assert_eq!(two.convex_hulls[0].hull_points().len(), 12);
            real_similar(hull_sum(two), 535801.503);
            real_similar(hull_sum(one), 61405.95106);
            real_similar(meta(&feature, "leftWidth"), 1481.840);
            real_similar(meta(&feature, "rightWidth"), 1512.290);
            // The seed is consumed.
            assert_eq!(picked[1].peaks[0].intensity, 0.0);
        }
    }
}

/// `MRMTransitionGroupPicker_test.cpp:641-677`.
#[test]
fn group_picker_class_test_find_largest_peak() {
    let picked: Vec<MSChromatogram> = (0..3_u8)
        .map(|k| MSChromatogram {
            peaks: vec![
                ChromatogramPeak::new(3120.0, 100.0 + f32::from(k)),
                ChromatogramPeak::new(4120.0, 200.0 + f32::from(k)),
            ],
            ..Default::default()
        })
        .collect();
    assert_eq!(
        MRMTransitionGroupPicker::find_largest_peak(&picked),
        Some((2, 1))
    );
}

fn with_borders(peaks: &[(f64, f32)], borders: &[(f32, f32)]) -> MSChromatogram {
    let mut arrays = empty_arrays();
    arrays[PeakPickerChromatogram::IDX_LEFTBORDER] =
        DataArray::new("leftWidth", borders.iter().map(|b| b.0).collect());
    arrays[PeakPickerChromatogram::IDX_RIGHTBORDER] =
        DataArray::new("rightWidth", borders.iter().map(|b| b.1).collect());
    MSChromatogram {
        peaks: peaks
            .iter()
            .map(|&(rt, intensity)| ChromatogramPeak::new(rt, intensity))
            .collect(),
        float_data_arrays: arrays,
        ..Default::default()
    }
}

/// `MRMTransitionGroupPicker_test.cpp:679-727`.
#[test]
fn group_picker_class_test_find_widest_peak_indices() {
    let picked = [
        with_borders(
            &[(110.0, 2000.0), (150.0, 6000.0), (190.0, 2000.0)],
            &[(100.0, 120.0), (120.0, 180.0), (180.0, 200.0)],
        ),
        with_borders(
            &[(150.0, 5500.0), (190.0, 2000.0)],
            &[(100.0, 180.0), (180.0, 200.0)],
        ),
        with_borders(
            &[(110.0, 2000.0), (150.0, 7000.0), (190.0, 2000.0)],
            &[(105.0, 115.0), (125.0, 175.0), (185.0, 195.0)],
        ),
    ];
    assert_eq!(
        MRMTransitionGroupPicker::find_widest_peak_indices(&picked).unwrap(),
        Some((1, 0))
    );
}

/// `MRMTransitionGroupPicker_test.cpp:730-841`.
#[test]
fn group_picker_class_test_remove_overlapping_features() {
    let original = [
        with_borders(
            &[(3120.0, 170.0), (3090.0, 170.0), (3060.0, 170.0)],
            &[(3100.0, 3140.0), (3070.0, 3120.0), (3050.0, 3090.0)],
        ),
        with_borders(
            &[(3120.0, 170.0), (3060.0, 170.0)],
            &[(3100.0, 3140.0), (3050.0, 3090.0)],
        ),
    ];
    let intensities =
        |c: &MSChromatogram| -> Vec<f32> { c.peaks.iter().map(|p| p.intensity).collect() };
    for (left, right, first, second) in [
        (3100.0, 3140.0, vec![0.0, 0.0, 170.0], vec![0.0, 170.0]),
        (3070.0, 3120.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0]),
        (3050.0, 3090.0, vec![170.0, 0.0, 0.0], vec![170.0, 0.0]),
    ] {
        let mut picked = original.clone();
        MRMTransitionGroupPicker::remove_overlapping_features(&mut picked, left, right);
        assert_eq!(intensities(&picked[0]), first);
        assert_eq!(intensities(&picked[1]), second);
    }
}

// ---------------------------------------------------------------------------
// Tier 1: MRMTransitionGroupPicker against the executed Release build
// ---------------------------------------------------------------------------

/// A C99 `%a` hex float as glibc prints it (`0x1.8p+1`, `-0x0p+0`, `inf`,
/// `nan`), parsed exactly.
fn parse_hex_float(text: &str) -> f64 {
    let (negative, body) = match text.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, text),
    };
    let value = match body {
        "inf" => f64::INFINITY,
        "nan" => f64::NAN,
        _ => {
            let body = body.strip_prefix("0x").expect("hex float");
            let (mantissa, exponent) = body.split_once('p').expect("hex exponent");
            let (int_part, frac_part) = mantissa.split_once('.').unwrap_or((mantissa, ""));
            let digits = format!("{int_part}{frac_part}");
            let bits = u64::from_str_radix(&digits, 16).expect("hex mantissa");
            let exponent: i32 = exponent.parse().expect("decimal exponent");
            let scale = exponent - 4 * i32::try_from(frac_part.len()).unwrap();
            // At most 53 significant bits, so the conversion is exact; the
            // scaling by a power of two is exact for the normal values here.
            (bits as f64) * 2f64.powi(scale)
        }
    };
    if negative { -value } else { value }
}

#[derive(Debug)]
enum Record {
    Text(String),
    Float(f64),
    Hull(Vec<(f64, f64)>),
}

fn push_feature(out: &mut Vec<(String, Record)>, prefix: &str, feature: &Feature) {
    let key = |name: &str| format!("{prefix}{name}");
    out.push((key("rt"), Record::Float(feature.base.rt)));
    out.push((key("mz"), Record::Float(feature.base.mz)));
    out.push((
        key("intensity"),
        Record::Float(f64::from(feature.base.intensity)),
    ));
    out.push((
        key("quality"),
        Record::Float(f64::from(feature.base.quality)),
    ));
    for (name, value) in &feature.base.metadata {
        let record = match value.data() {
            MetaValueData::Float(x) => Record::Float(*x),
            MetaValueData::Integer(i) => Record::Text(format!("i:{i}")),
            MetaValueData::String(s) => Record::Text(format!("s:{s}")),
            other => Record::Text(format!("?:{other:?}")),
        };
        out.push((key(&format!("meta.{name}")), record));
    }
    out.push((
        key("hulls"),
        Record::Text(feature.convex_hulls.len().to_string()),
    ));
    for (h, hull) in feature.convex_hulls.iter().enumerate() {
        out.push((
            key(&format!("hull{h}")),
            Record::Hull(hull.hull_points().iter().map(|p| (p.rt, p.mz)).collect()),
        ));
    }
}

fn oracle_group(name: &str) -> Group {
    let mut group = Group::new();
    group.set_transition_group_id(name);
    for row in rows_all(include_str!("data/openswath_peak_picking_groups.tsv")) {
        if row[0] != name {
            continue;
        }
        let numbers =
            |csv: &str| -> Vec<f64> { csv.split(',').map(|t| t.parse::<f64>().unwrap()).collect() };
        let rts = numbers(row[7]);
        let intensities = numbers(row[8]);
        let mut chromatogram = srm_chromatogram(row[2], &rts, &intensities);
        if row[5] != "-" {
            chromatogram = with_product_mz(chromatogram, row[5].parse().unwrap());
        }
        if row[6] != "-" {
            chromatogram.metadata.insert(
                "precursor_mz".into(),
                MetaValue::try_from(row[6].parse::<f64>().unwrap()).unwrap(),
            );
        }
        if row[1] == "fragment" {
            group.add_chromatogram(chromatogram, row[2]).unwrap();
            let transition = SimpleTransition {
                detecting: row[3] == "1",
                quantifying: row[4] == "1",
                ..SimpleTransition::new(row[2], -1.0)
            };
            group.add_transition(transition, row[2]).unwrap();
        } else {
            group
                .add_precursor_chromatogram(chromatogram, row[2])
                .unwrap();
        }
    }
    group
}

/// Every data line of a `#`-commented TSV without a header row.
fn rows_all(text: &str) -> impl Iterator<Item = Vec<&str>> {
    text.lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| line.split('\t').collect())
}

fn oracle_settings(settings: &str) -> Vec<(String, ParamValue)> {
    settings
        .split(';')
        .filter(|kv| !kv.is_empty())
        .map(|kv| {
            let (key, typed) = kv.split_once('=').unwrap();
            let (kind, value) = typed.split_once(':').unwrap();
            let value = match kind {
                "i" => ParamValue::Integer(value.parse().unwrap()),
                "f" => ParamValue::Float(value.parse().unwrap()),
                _ => ParamValue::String(value.into()),
            };
            (key.to_owned(), value)
        })
        .collect()
}

fn run_oracle_case(group_name: &str, settings: &str) -> Vec<(String, Record)> {
    let mut group = oracle_group(group_name);
    let mut param = MRMTransitionGroupPicker::defaults().unwrap();
    for (key, value) in oracle_settings(settings) {
        set(&mut param, &key, value);
    }
    let picker = MRMTransitionGroupPicker::from_param(&param).unwrap();
    let mut out = Vec::new();
    if let Err(error) = picker.pick_transition_group(&mut group, &mut ids()) {
        out.push(("error".into(), Record::Text(error.to_string())));
        assert!(group.features().is_empty(), "failed pick changed the group");
        return out;
    }
    out.push((
        "features".into(),
        Record::Text(group.features().len().to_string()),
    ));
    for (i, feature) in group.features().iter().enumerate() {
        let prefix = format!("f{i}.");
        push_feature(&mut out, &prefix, feature);
        out.push((
            format!("{prefix}subs"),
            Record::Text(feature.features().len().to_string()),
        ));
        for (k, sub) in feature.features().iter().enumerate() {
            push_feature(&mut out, &format!("{prefix}s{k}."), sub);
        }
        let ids: Vec<&str> = feature.precursor_feature_ids().collect();
        out.push((
            format!("{prefix}precs"),
            Record::Text(ids.len().to_string()),
        ));
        for (k, id) in ids.iter().enumerate() {
            push_feature(
                &mut out,
                &format!("{prefix}p{k}."),
                feature.precursor_feature(id).unwrap(),
            );
        }
    }
    out
}

fn same_float(actual: f64, expected: f64) -> bool {
    actual.to_bits() == expected.to_bits() || (actual.is_nan() && expected.is_nan())
}

/// Every case of `tests/data/openswath_peak_picking_configs.tsv`, replayed
/// against the records the Release build wrote on ibminode06
/// (`tests/data/openswath_peak_picking_oracle.tsv`): the same peak groups,
/// sub-features, metadata keys and hull points, every floating-point value
/// bit for bit, and an error wherever the source threw.
#[test]
fn group_picker_matches_the_release_build_bit_for_bit() {
    let oracle = include_str!("data/openswath_peak_picking_oracle.tsv");
    let mut cases = 0;
    let mut mismatches = Vec::new();
    for config in rows_all(include_str!("data/openswath_peak_picking_configs.tsv")) {
        let case = config[0];
        let expected: Vec<(&str, &str)> = oracle
            .lines()
            .filter_map(|line| {
                let mut fields = line.splitn(3, '\t');
                (fields.next() == Some(case))
                    .then(|| (fields.next().unwrap(), fields.next().unwrap()))
            })
            .collect();
        assert!(!expected.is_empty(), "{case}: no oracle records");
        let actual = run_oracle_case(config[1], config.get(2).copied().unwrap_or(""));
        if expected[0].0 == "error" {
            assert_eq!(
                actual.len(),
                1,
                "{case}: expected {}, got {actual:?}",
                expected[0].1
            );
            assert_eq!(
                actual[0].0, "error",
                "{case}: the source threw {}",
                expected[0].1
            );
            cases += 1;
            continue;
        }
        let actual_keys: Vec<_> = actual.iter().map(|(k, _)| k.as_str()).collect();
        let expected_keys: Vec<_> = expected.iter().map(|(k, _)| *k).collect();
        if actual_keys != expected_keys {
            mismatches.push(format!(
                "{case}: record keys differ: Rust {actual:?}, Release build {expected:?}"
            ));
            cases += 1;
            continue;
        }
        for ((key, record), (_, text)) in actual.iter().zip(&expected) {
            let ok = match record {
                Record::Text(value) => value == text,
                Record::Float(value) => {
                    let text = text.strip_prefix("f:").unwrap_or(text);
                    same_float(*value, parse_hex_float(text))
                }
                Record::Hull(points) => {
                    let expected: Vec<(f64, f64)> = if text.is_empty() {
                        Vec::new()
                    } else {
                        text.split(',')
                            .map(|pair| {
                                let (x, y) = pair.split_once(':').unwrap();
                                (parse_hex_float(x), parse_hex_float(y))
                            })
                            .collect()
                    };
                    points.len() == expected.len()
                        && points
                            .iter()
                            .zip(&expected)
                            .all(|(a, b)| same_float(a.0, b.0) && same_float(a.1, b.1))
                }
            };
            if !ok {
                mismatches.push(format!(
                    "{case}: {key}: Rust {record:?}, Release build {text}"
                ));
            }
        }
        cases += 1;
    }
    assert_eq!(cases, 32);
    assert!(mismatches.is_empty(), "{}", mismatches.join("\n"));
}

#[test]
fn hex_float_parser_reads_glibc_output() {
    assert_eq!(parse_hex_float("0x1.8p+1"), 3.0);
    assert_eq!(parse_hex_float("-0x1p-2"), -0.25);
    assert_eq!(parse_hex_float("0x0p+0"), 0.0);
    assert_eq!(parse_hex_float("0x1.753527c57eb84p+10"), 1492.8305524575708);
    assert!(parse_hex_float("nan").is_nan());
}

// ---------------------------------------------------------------------------
// PeakPickerMobilogram
// ---------------------------------------------------------------------------

/// `get_mobilogram(0)` of `PeakPickerMobilogram_test.cpp:21-44`: the first 60
/// of the 61 listed mobilities (the loop stops at 60).
fn class_test_mobilogram() -> Mobilogram {
    let mobilities = [
        0.93734956, 0.93835497, 0.9393905, 0.9404265, 0.9414631, 0.9425003, 0.9435083, 0.9445466,
        0.9455854, 0.94662476, 0.947635, 0.94867545, 0.9496867, 0.95072824, 0.95177037, 0.9527832,
        0.9538264, 0.95484036, 0.9558847, 0.9568997, 0.9579451, 0.9589612, 0.96000767, 0.9610248,
        0.96207243, 0.96309066, 0.9641094, 0.96515864, 0.9661785, 0.96719885, 0.96824974,
        0.9692712, 0.9702931, 0.97131556, 0.9723687, 0.97339225, 0.9744164, 0.975441, 0.9764661,
        0.977522, 0.9785482, 0.979575, 0.98060226, 0.98163015, 0.9826585, 0.9836874, 0.98471683,
        0.9857468, 0.9867773, 0.98780835, 0.98883986, 0.989872, 0.9909046, 0.9919378, 0.99297154,
        0.99397534, 0.99501014, 0.9960454, 0.9970813, 0.9981177,
    ];
    let intensities = [
        341.995757,
        200.000741,
        325.0076900000001,
        79.002655,
        227.00538,
        487.00108,
        576.0032550000001,
        526.996605,
        715.9864020000001,
        778.010235,
        686.003534,
        548.0023570000001,
        481.00395999999995,
        547.00563,
        950.006949,
        1453.9987189999997,
        831.026847,
        1839.009966,
        2306.021277,
        2588.9953999999993,
        2126.0273610000004,
        3254.0281579999996,
        3961.0108749999995,
        4708.020719000002,
        4991.072095,
        5230.084562000001,
        6559.103413000001,
        6486.9137040000005,
        9512.035332,
        14461.755780000001,
        16570.804681,
        19887.687351000004,
        22368.220470000004,
        31683.985819,
        37679.199531000006,
        46487.67955,
        48958.673168000016,
        52820.41709300001,
        58070.474128999995,
        57655.34146700001,
        58849.449743000005,
        56964.34470800001,
        56100.919934,
        48964.276839,
        43019.251419,
        39263.106511000005,
        28199.275315000003,
        25018.205981999996,
        20248.883416,
        14727.949587000001,
        10372.823519,
        6309.068338999999,
        5370.119231000001,
        2253.037892,
        2025.998686,
        868.9924269999999,
        209.0077,
        44.99693,
        317.999914,
        57.00195,
    ];
    Mobilogram::from_peaks(
        mobilities
            .iter()
            .zip(intensities)
            .map(|(&mobility, intensity)| MobilityPeak1D::new(mobility, intensity as f32))
            .collect(),
    )
}

/// `PeakPickerMobilogram.cpp:26-45`.
#[test]
fn mobilogram_picker_defaults_follow_the_source_constructor() {
    let defaults = PeakPickerMobilogram::defaults().unwrap();
    let keys: Vec<_> = defaults
        .iter()
        .unwrap()
        .map(|item| item.key.clone())
        .collect();
    assert_eq!(
        keys,
        [
            "sgolay_frame_length",
            "sgolay_polynomial_order",
            "gauss_width",
            "use_gauss",
            "peak_width",
            "signal_to_noise",
            "sn_win_len",
            "sn_bin_count",
            "write_sn_log_messages",
            "remove_overlapping_peaks",
            "method",
        ]
    );
    assert_eq!(
        defaults.value("sgolay_frame_length").unwrap(),
        &ParamValue::Integer(9)
    );
    assert_eq!(
        defaults.value("gauss_width").unwrap(),
        &ParamValue::Float(0.002)
    );
    assert_eq!(string(&defaults, "use_gauss"), "false");
    // `defaults_.setValue("sn_win_len", 1, ...)`: an integer parameter.
    assert_eq!(
        defaults.value("sn_win_len").unwrap(),
        &ParamValue::Integer(1)
    );
    assert_eq!(
        defaults.value("sn_bin_count").unwrap(),
        &ParamValue::Integer(4)
    );
    assert_eq!(
        defaults.entry("method").unwrap().valid_strings,
        ["legacy", "corrected", "crawdad"]
    );
    assert_eq!(
        defaults.description("method").unwrap(),
        "Which method to choose for mobilogram peak-picking (OpenSWATH legacy on raw data, corrected picking on smoothed mobilogram)."
    );
    let picker = PeakPickerMobilogram::from_param(&defaults).unwrap();
    assert_eq!(
        picker.smoothing,
        ChromatogramSmoothing::SavitzkyGolay {
            frame_length: 9,
            polynomial_order: 3
        }
    );
    assert!(picker.to_param().unwrap().source_equal(&defaults).unwrap());
    let mut param = defaults.clone();
    set(&mut param, "method", ParamValue::String("crawdad".into()));
    assert!(PeakPickerMobilogram::from_param(&param).is_err());
    let mut param = defaults;
    set(
        &mut param,
        "remove_overlapping_peaks",
        ParamValue::String("true".into()),
    );
    set(&mut param, "sn_bin_count", ParamValue::Integer(-3));
    let picker = PeakPickerMobilogram::from_param(&param).unwrap();
    assert!(picker.remove_overlapping_peaks);
    // Stored without effect, as in the source.
    assert_eq!(picker.sn_bin_count, -3);
}

/// `PeakPickerMobilogram_test.cpp:67-84`: `getParameters`, `method`
/// `corrected`, `use_gauss` `false`, `setParameters`, then
/// `pickMobilogram(mobilogram, picked, smoothed)`.
#[test]
fn mobilogram_picker_class_test() {
    let mut param = PeakPickerMobilogram::default().to_param().unwrap();
    set(&mut param, "method", ParamValue::String("corrected".into()));
    set(&mut param, "use_gauss", ParamValue::String("false".into()));
    let picker = PeakPickerMobilogram::from_param(&param).unwrap();
    let input = class_test_mobilogram();
    let result = picker.pick_mobilogram(&input).unwrap();
    let picked = &result.picked;
    real_similar(f64::from(picked.peaks[0].intensity), 58956.1);
    real_similar(picked.peaks[0].mobility, 0.978364);
    let array = |index: usize| f64::from(picked.float_data_arrays[index].data[0]);
    real_similar(array(PeakPickerMobilogram::IDX_ABUNDANCE), 884145.0);
    real_similar(array(PeakPickerMobilogram::IDX_LEFTBORDER), 0.948675);
    real_similar(array(PeakPickerMobilogram::IDX_RIGHTBORDER), 0.997081);
    let names: Vec<_> = picked
        .float_data_arrays
        .iter()
        .map(|a| a.name.as_str())
        .collect();
    assert_eq!(
        names,
        [
            "FWHM",
            "IntegratedIntensity",
            "leftWidth",
            "rightWidth",
            "leftWidthIndex",
            "rightWidthIndex"
        ]
    );
    // The index arrays name the border samples; the public per-peak vectors
    // hold the same values.
    let left = result.left_widths[0];
    let right = result.right_widths[0];
    assert_eq!(
        array(PeakPickerMobilogram::IDX_OF_LEFTBORDER_IDX),
        left as f64
    );
    assert_eq!(
        array(PeakPickerMobilogram::IDX_OF_RIGHTBORDER_IDX),
        right as f64
    );
    assert_eq!(
        array(PeakPickerMobilogram::IDX_LEFTBORDER),
        f64::from(input.peaks[left].mobility as f32)
    );
    let sum = input.peaks[left..=right]
        .iter()
        .fold(0.0, |s, p| s + f64::from(p.intensity));
    assert_eq!(result.integrated_intensities[0], sum);
    assert_eq!(result.smoothed.len(), input.len());
    assert!(
        PeakPickerMobilogram::default()
            .pick_mobilogram(&Mobilogram::default())
            .unwrap()
            .picked
            .peaks
            .is_empty()
    );
}

/// `filterTopPeak` and `findHighestPeak_` (`PeakPickerMobilogram.cpp:242-436`).
#[test]
fn mobilogram_top_peak_filters() {
    let input = class_test_mobilogram();
    let result = PeakPickerMobilogram::default()
        .pick_mobilogram(&input)
        .unwrap();
    let (left, right) = (result.left_widths[0], result.right_widths[0]);
    let mut many = vec![input.clone(), input.clone()];
    let positions = PeakPickerMobilogram::filter_top_peak(&result.picked, &mut many).unwrap();
    // `apex` is the picked peak's index, as the source returns it.
    assert_eq!(
        positions,
        PeakPositions {
            left,
            apex: 0,
            right
        }
    );
    for cut in &many {
        assert_eq!(cut.len(), right - left + 1);
        assert_eq!(cut.peaks[0], input.peaks[left]);
    }
    let mut one = input.clone();
    assert_eq!(
        PeakPickerMobilogram::filter_top_peak_single(&result.picked, &mut one).unwrap(),
        positions
    );
    assert_eq!(one, many[0]);
    // No picked peaks: the whole first mobilogram, apex in its middle.
    let mut empty_picked = result.picked.clone();
    for array in &mut empty_picked.float_data_arrays {
        array.data.clear();
    }
    empty_picked.peaks.clear();
    let mut whole = vec![input.clone()];
    assert_eq!(
        PeakPickerMobilogram::filter_top_peak(&empty_picked, &mut whole).unwrap(),
        PeakPositions {
            left: 0,
            apex: 30,
            right: 59
        }
    );
    assert_eq!(whole[0], input);
    // Missing arrays or no mobilograms: all zeros, nothing cut.
    assert_eq!(
        PeakPickerMobilogram::filter_top_peak(&Mobilogram::default(), &mut many).unwrap(),
        PeakPositions::default()
    );
    assert_eq!(
        PeakPickerMobilogram::filter_top_peak(&result.picked, &mut []).unwrap(),
        PeakPositions::default()
    );
}

// ---------------------------------------------------------------------------
// MasstraceCorrelator
// ---------------------------------------------------------------------------

/// `MasstraceCorrelator.cpp:34-45`: three stored-only parameters under the
/// source's (misnamed) handler.
#[test]
fn masstrace_correlator_parameters() {
    let defaults = MasstraceCorrelator::defaults().unwrap();
    assert_eq!(
        defaults.value("sgolay_frame_length").unwrap(),
        &ParamValue::Integer(15)
    );
    assert_eq!(
        defaults.value("sgolay_polynomial_order").unwrap(),
        &ParamValue::Integer(3)
    );
    // `defaults_.setValue("gauss_width", 50, ...)`: an integer parameter.
    assert_eq!(
        defaults.value("gauss_width").unwrap(),
        &ParamValue::Integer(50)
    );
    assert_eq!(
        defaults.description("sgolay_polynomial_order").unwrap(),
        "Order or the polynomial that is fitted."
    );
    let mut param = defaults;
    set(&mut param, "gauss_width", ParamValue::Integer(7));
    set(&mut param, "unknown", ParamValue::Integer(1));
    let (correlator, warnings) = MasstraceCorrelator::from_param_with_warnings(&param).unwrap();
    assert_eq!(correlator.gauss_width, 7);
    assert_eq!(warnings.len(), 1);
    assert!(
        warnings[0].contains(MASSTRACE_CORRELATOR_NAME),
        "{warnings:?}"
    );
    assert_eq!(MASSTRACE_CORRELATOR_NAME, "MRMFeatureFinderScoring");
}

fn traces(rts: &[f64], values: &[f64]) -> Vec<(f64, f64)> {
    rts.iter().copied().zip(values.iter().copied()).collect()
}

/// `MasstraceCorrelator_test.cpp:211-283`.
#[test]
fn masstrace_correlator_class_test_score_hullpoints() {
    let mut data1 = vec![0.0, 1.0, 3.0, 5.0, 2.0, 0.0];
    let mut data2 = vec![1.0, 3.0, 5.0, 2.0, 0.0, 0.0];
    let rt1 = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let rt2 = [0.0, 1.0, 2.0, 4.0, 5.0, 6.0];
    standardize_data(&mut data1);
    standardize_data(&mut data2);
    let mut result = calculate_cross_correlation(&data1, &data2, 2, 1).unwrap();
    for entry in result.iter_mut() {
        entry.1 /= 6.0;
    }
    let lags: Vec<i32> = result.iter().map(|e| e.0).collect();
    assert_eq!(lags, [-2, -1, 0, 1, 2]);
    for (index, expected) in [
        (4, -0.7374631),
        (3, -0.567846),
        (2, 0.4159292),
        (1, 0.8215339),
        (0, 0.15634218),
    ] {
        real_similar(result.data[index].1, expected);
    }

    let correlator = MasstraceCorrelator::default();
    let min_pearson_score = -1.1;
    let max_lag = 6;
    let same_rt = correlator
        .score_hullpoints(
            &traces(&rt1, &data1),
            &traces(&rt1, &data2),
            min_pearson_score,
            max_lag,
            0.1,
        )
        .unwrap();
    let (lag, lag_intensity) = same_rt.lag.unwrap();
    assert_eq!(lag, -1);
    real_similar(lag_intensity, 0.821534);
    real_similar(same_rt.pearson_score, 0.41593);

    let wide = correlator
        .score_hullpoints(
            &traces(&rt1, &data1),
            &traces(&rt2, &data2),
            min_pearson_score,
            max_lag,
            1.5,
        )
        .unwrap();
    let (lag, lag_intensity) = wide.lag.unwrap();
    assert_eq!(lag, -1);
    real_similar(lag_intensity, 0.821534);
    real_similar(wide.pearson_score, 0.41593);

    let narrow = correlator
        .score_hullpoints(
            &traces(&rt1, &data1),
            &traces(&rt2, &data2),
            min_pearson_score,
            max_lag,
            0.1,
        )
        .unwrap();
    let (lag, lag_intensity) = narrow.lag.unwrap();
    assert_eq!(lag, -1);
    real_similar(lag_intensity, 0.625368);
    real_similar(narrow.pearson_score, 0.405604);

    // Below the minimum Pearson score the lag is not computed.
    let skipped = correlator
        .score_hullpoints(
            &traces(&rt1, &data1),
            &traces(&rt1, &data2),
            0.9,
            max_lag,
            0.1,
        )
        .unwrap();
    assert_eq!(skipped.lag, None);
    assert!(correlator.score_hullpoints(&[], &[], 0.0, 1, 0.1).is_err());
}

/// `MasstraceCorrelator_test.cpp:285-315`, on the retained
/// `Masstraces_Testdata.consensusXML`, sorted by descending intensity.
#[cfg(feature = "consensusxml")]
#[test]
fn masstrace_correlator_class_test_create_pseudo_spectra() {
    let mut map = openms::format::consensusxml::read(
        include_bytes!("data/openswath_peak_picking_masstraces.consensusXML").as_slice(),
    )
    .unwrap();
    map.sort_by_intensity(true).unwrap();
    let correlator = MasstraceCorrelator::default();

    let spectra = correlator
        .create_pseudo_spectra(&map, 0, 0.7, 1, 3.0)
        .unwrap();
    assert_eq!(spectra.len(), 3);
    assert_eq!(spectra[1].spectrum.peaks.len(), 1);
    real_similar(spectra[1].spectrum.rt, 4203.0);
    real_similar(spectra[1].spectrum.peaks[0].mz, 668.5);
    for spectrum in &spectra {
        assert_eq!(spectrum.spectrum.ms_level, 2);
        assert_eq!(
            spectrum.correlations.len(),
            spectrum.spectrum.peaks.len() - 1
        );
    }

    let spectra = correlator
        .create_pseudo_spectra(&map, 1, 0.7, 1, 3.0)
        .unwrap();
    assert_eq!(spectra.len(), 2);
    assert_eq!(spectra[0].spectrum.peaks.len(), 2);
    assert_eq!(spectra[1].spectrum.peaks.len(), 2);
    real_similar(spectra[0].spectrum.rt, 5201.0);
    real_similar(spectra[0].spectrum.peaks[0].mz, 568.5);
    real_similar(spectra[1].spectrum.rt, 5203.0);
    real_similar(spectra[1].spectrum.peaks[0].mz, 768.5);

    // The cache: one RT-sorted trace, maximum and RT per feature.
    let cache = correlator.create_consensus_map_cache(&map).unwrap();
    assert_eq!(cache.feature_points.len(), map.features.len());
    for (points, feature) in cache.feature_points.iter().zip(&map.features) {
        assert_eq!(points.len(), feature.handles().len());
        assert!(points.windows(2).all(|w| w[0].0 <= w[1].0));
    }
    // A work budget that cannot cover one pair is refused.
    let starved = MasstraceCorrelator {
        max_work: 1,
        ..MasstraceCorrelator::default()
    };
    assert!(starved.create_pseudo_spectra(&map, 0, 0.7, 1, 3.0).is_err());
}
