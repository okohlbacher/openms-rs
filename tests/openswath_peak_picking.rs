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

use openms::analysis::peak_integrator::{
    BaselineType, IntegrationMethod, PEAK_INTEGRATOR_NAME, PeakIntegrator,
};
use openms::kernel::{ChromatogramPeak, MSChromatogram, MSSpectrum, Peak1D};
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

/// `PeakIntegrator_test.cpp:506-585` and the spectrum twin at `669-749`.
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

/// `PeakIntegrator_test.cpp:801-850` and the spectrum twin at `901-950`.
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
