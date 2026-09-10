// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Integration of fitted native traces. The cutoff trace and 28-point fitted
//! length below are literal golden data from EmgGradientDescent_test.cpp at
//! OpenMS4-core 7c029e8cdba6abab503708ecdd56f6ab55e38ce4. No C++ was executed.

use openms::analysis::emg::EmgGradientDescent;
use openms::analysis::peak_integrator::{
    BaselineType, IntegrationMethod, PeakIntegrator, PeakShapeMetrics,
};
use openms::kernel::DataArray;
use openms::{ChromatogramPeak, MSChromatogram, MSSpectrum, Peak1D};

fn source_trace(name: &str) -> MSChromatogram {
    let peaks = include_str!("data/emg_traces.tsv")
        .lines()
        .skip(1)
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|row| row[0] == name)
        .map(|row| {
            ChromatogramPeak::new(
                f64::from_bits(u64::from_str_radix(row[4], 16).unwrap()),
                f32::from_bits(u32::from_str_radix(row[6], 16).unwrap()),
            )
        })
        .collect();
    MSChromatogram::from_peaks(peaks)
}

fn cutoff() -> MSChromatogram {
    let mut input = source_trace("cutoff_min");
    input.name = "source cutoff peak".into();
    input
        .metadata
        .insert("origin".into(), "pinned class test".into());
    input
        .float_data_arrays
        .push(DataArray::new("observed quality", vec![1.0; input.len()]));
    input.integer_data_arrays.push(DataArray::new(
        "original index",
        (0..input.len() as i32).collect(),
    ));
    input.string_data_arrays.push(DataArray::new(
        "observed label",
        vec!["sample".into(); input.len()],
    ));
    input
}
fn bounds(input: &MSChromatogram) -> (f64, f64) {
    (input.peaks[0].rt, input.peaks.last().unwrap().rt)
}
fn spectrum(input: &MSChromatogram) -> MSSpectrum {
    MSSpectrum {
        peaks: input
            .peaks
            .iter()
            .map(|p| Peak1D::new(p.rt, p.intensity))
            .collect(),
        metadata: input.metadata.clone(),
        float_data_arrays: input.float_data_arrays.clone(),
        integer_data_arrays: input.integer_data_arrays.clone(),
        string_data_arrays: input.string_data_arrays.clone(),
        ..Default::default()
    }
}
fn enabled(method: IntegrationMethod, fit: &EmgGradientDescent) -> PeakIntegrator {
    PeakIntegrator {
        integration_method: method,
        emg: Some(fit.clone()),
        ..Default::default()
    }
}

#[test]
fn all_integration_methods_include_the_source_extrapolated_span_on_both_containers() {
    let input = cutoff();
    let original = input.clone();
    let fitter = EmgGradientDescent::default();
    let (left, right) = bounds(&input);
    let fitted = fitter
        .fit_chromatogram(&input, Some(left), Some(right))
        .unwrap();
    assert_eq!(fitted.chromatogram.len(), 28); // Literal upstream fit output size.
    let (fitted_left, fitted_right) = bounds(&fitted.chromatogram);
    assert!(fitted_right > right);
    assert!(fitted.estimate.parameters.h.is_finite());
    // Parameters are typed diagnostics, not an invalid four-element peak array.
    fitted.chromatogram.validate().unwrap();
    let input_spectrum = spectrum(&input);
    let before_spectrum = input_spectrum.clone();
    for method in [
        IntegrationMethod::IntensitySum,
        IntegrationMethod::Trapezoid,
        IntegrationMethod::Simpson,
    ] {
        let native = PeakIntegrator {
            integration_method: method,
            ..Default::default()
        };
        let expected = native
            .integrate_chromatogram(&fitted.chromatogram, fitted_left, fitted_right)
            .unwrap();
        let actual = enabled(method, &fitter)
            .integrate_chromatogram(&input, left, right)
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(actual.hull_points.len(), 28);
        assert!(actual.area > 0.0);
        assert_eq!(
            actual,
            enabled(method, &fitter)
                .integrate_spectrum(&input_spectrum, left, right)
                .unwrap()
        );
        // A crop of the fitted result to original bounds omits the restored tail.
        assert!(
            actual.area
                > native
                    .integrate_chromatogram(&fitted.chromatogram, left, right)
                    .unwrap()
                    .area
        );
    }
    assert_eq!(input, original);
    assert_eq!(input_spectrum, before_spectrum);
}

#[test]
fn fitted_background_and_shape_use_expanded_endpoints_but_keep_supplied_height_and_apex() {
    let input = cutoff();
    let input_spectrum = spectrum(&input);
    let fitter = EmgGradientDescent::default();
    let (left, right) = bounds(&input);
    let fitted = fitter
        .fit_chromatogram(&input, Some(left), Some(right))
        .unwrap()
        .chromatogram;
    let (fit_left, fit_right) = bounds(&fitted);
    let supplied_apex = input.peaks[6].rt;
    let supplied_height = f64::from(input.peaks[7].intensity) / 5.0;
    for method in [
        IntegrationMethod::IntensitySum,
        IntegrationMethod::Trapezoid,
        IntegrationMethod::Simpson,
    ] {
        for baseline in [
            BaselineType::BaseToBase,
            BaselineType::VerticalDivisionMin,
            BaselineType::VerticalDivisionMax,
        ] {
            let native = PeakIntegrator {
                integration_method: method,
                baseline_type: baseline,
                ..Default::default()
            };
            let emg = PeakIntegrator {
                emg: Some(fitter.clone()),
                ..native.clone()
            };
            let expected = native
                .estimate_background_chromatogram(&fitted, fit_left, fit_right, supplied_apex)
                .unwrap();
            assert_eq!(
                emg.estimate_background_chromatogram(&input, left, right, supplied_apex)
                    .unwrap(),
                expected
            );
            assert_eq!(
                emg.estimate_background_spectrum(&input_spectrum, left, right, supplied_apex)
                    .unwrap(),
                expected
            );
        }
    }
    let native = PeakIntegrator::default();
    let expected = native
        .calculate_shape_metrics_chromatogram(
            &fitted,
            fit_left,
            fit_right,
            supplied_height,
            supplied_apex,
        )
        .unwrap();
    let emg = enabled(IntegrationMethod::IntensitySum, &fitter);
    assert_eq!(
        emg.calculate_shape_metrics_chromatogram(
            &input,
            left,
            right,
            supplied_height,
            supplied_apex
        )
        .unwrap(),
        expected
    );
    assert_eq!(
        emg.calculate_shape_metrics_spectrum(
            &input_spectrum,
            left,
            right,
            supplied_height,
            supplied_apex
        )
        .unwrap(),
        expected
    );
    assert_eq!(expected.points_across_baseline, 28);
    assert_eq!(expected.total_width, fit_right - fit_left);
    let recomputed = native
        .integrate_chromatogram(&fitted, fit_left, fit_right)
        .unwrap();
    assert_ne!(
        expected,
        native
            .calculate_shape_metrics_chromatogram(
                &fitted,
                fit_left,
                fit_right,
                recomputed.height,
                recomputed.apex_pos
            )
            .unwrap()
    );
}

#[test]
fn cropped_fit_excludes_outside_measurements_and_can_disable_extrapolation() {
    let mut input = cutoff();
    let (left, right) = bounds(&input);
    input
        .peaks
        .insert(0, ChromatogramPeak::new(left - 1.0, 1e8));
    input.peaks.push(ChromatogramPeak::new(right + 1.0, 1e8));
    input.float_data_arrays.clear();
    input.integer_data_arrays.clear();
    input.string_data_arrays.clear();
    let before = input.clone();
    let fitter = EmgGradientDescent {
        compute_additional_points: false,
        ..Default::default()
    };
    let emg = enabled(IntegrationMethod::IntensitySum, &fitter);
    let output = emg.integrate_chromatogram(&input, left, right).unwrap();
    assert_eq!(output.hull_points.len(), 12);
    assert_eq!(output.hull_points[0][0], left);
    assert_eq!(output.hull_points.last().unwrap()[0], right);
    assert!(output.height < 1e8);
    assert_eq!(
        output,
        emg.integrate_chromatogram(&cutoff(), left, right).unwrap()
    );
    assert_eq!(input, before);
    assert!(PeakIntegrator::default().emg.is_none());
}

#[test]
fn zero_is_a_literal_fit_boundary_and_never_the_cpp_unbounded_sentinel() {
    let mut input = cutoff();
    // Unrelated negative-coordinate observations must not enter a zero-bounded fit.
    // Retain the pinned peak's coordinate scale and numerically supported fit.
    input.peaks.insert(0, ChromatogramPeak::new(-1.0, 1e8));
    input.peaks.insert(0, ChromatogramPeak::new(-2.0, 1e8));
    input.float_data_arrays.clear();
    input.integer_data_arrays.clear();
    input.string_data_arrays.clear();
    let fitter = EmgGradientDescent {
        compute_additional_points: false,
        ..Default::default()
    };
    let right = input.peaks.last().unwrap().rt;
    let fit = fitter
        .fit_chromatogram(&input, Some(0.0), Some(right))
        .unwrap();
    let output = enabled(IntegrationMethod::IntensitySum, &fitter)
        .integrate_chromatogram(&input, 0.0, right)
        .unwrap();
    assert!(output.hull_points.iter().all(|point| point[0] >= 0.0));
    assert!(output.hull_points.len() < input.len());
    assert_eq!(output.hull_points.len(), fit.chromatogram.len());
    let (a, b) = bounds(&fit.chromatogram);
    assert_eq!(
        output,
        PeakIntegrator::default()
            .integrate_chromatogram(&fit.chromatogram, a, b)
            .unwrap()
    );
}

#[test]
fn emg_errors_and_budgets_leave_input_and_annotations_unchanged() {
    let input = cutoff();
    let original = input.clone();
    let (left, right) = bounds(&input);
    for fitter in [
        EmgGradientDescent {
            max_iterations: 0,
            ..Default::default()
        },
        EmgGradientDescent {
            max_evaluations: 1,
            ..Default::default()
        },
        EmgGradientDescent {
            max_points: 1,
            ..Default::default()
        },
    ] {
        let emg = enabled(IntegrationMethod::Simpson, &fitter);
        assert!(emg.integrate_chromatogram(&input, left, right).is_err());
        assert!(
            emg.estimate_background_chromatogram(&input, left, right, input.peaks[7].rt)
                .is_err()
        );
        assert!(
            emg.calculate_shape_metrics_chromatogram(
                &input,
                left,
                right,
                2000.0,
                input.peaks[7].rt
            )
            .is_err()
        );
        assert_eq!(input, original);
    }
    // Original12 points fit the integration limit; restored28 points do not.
    let limited = PeakIntegrator {
        max_points: input.len(),
        emg: Some(EmgGradientDescent::default()),
        ..Default::default()
    };
    assert!(limited.integrate_chromatogram(&input, left, right).is_err());
    assert_eq!(input, original);
    let empty = MSChromatogram::default();
    let emg = enabled(
        IntegrationMethod::IntensitySum,
        &EmgGradientDescent::default(),
    );
    assert!(emg.integrate_chromatogram(&empty, 0.0, 1.0).is_err());
    assert!(
        emg.estimate_background_chromatogram(&empty, 0.0, 1.0, 0.5)
            .is_err()
    );
    // Source shape returns its empty result before invoking the fitter.
    assert_eq!(
        emg.calculate_shape_metrics_chromatogram(&empty, 0.0, 1.0, 0.0, 0.5)
            .unwrap(),
        PeakShapeMetrics::default()
    );
}

#[test]
fn left_extrapolation_changes_baseline_endpoints_and_shape_span() {
    let input = source_trace("saturated_cutoff_min");
    let (left, right) = bounds(&input);
    let fitter = EmgGradientDescent::default();
    let fitted = fitter
        .fit_chromatogram(&input, Some(left), Some(right))
        .unwrap()
        .chromatogram;
    assert_eq!(fitted.len(), 71); // Source's 66-point left-cutoff fit golden.
    let (fit_left, fit_right) = bounds(&fitted);
    assert!(fit_left < left);
    assert_eq!(fit_right, right);
    let emg = PeakIntegrator {
        emg: Some(fitter),
        baseline_type: BaselineType::VerticalDivisionMin,
        ..Default::default()
    };
    let area = emg.integrate_chromatogram(&input, left, right).unwrap();
    assert_eq!(area.hull_points[0][0], fit_left);
    let apex = input.peaks[8].rt;
    let minimum = f64::from(
        fitted.peaks[0]
            .intensity
            .min(fitted.peaks.last().unwrap().intensity),
    );
    let background = emg
        .estimate_background_chromatogram(&input, left, right, apex)
        .unwrap();
    assert_eq!(background.height, minimum);
    assert_eq!(background.area, minimum * 71.0);
    let shape = emg
        .calculate_shape_metrics_chromatogram(
            &input,
            left,
            right,
            f64::from(input.peaks[8].intensity),
            apex,
        )
        .unwrap();
    assert_eq!(shape.total_width, fit_right - fit_left);
    assert_eq!(shape.points_across_baseline, 71);
}
