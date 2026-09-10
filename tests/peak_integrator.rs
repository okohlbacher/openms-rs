// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::analysis::peak_integrator::{
    BaselineType, IntegrationMethod, PeakIntegrator, PeakShapeMetrics,
};
use openms::error::Error;
use openms::kernel::{ChromatogramPeak, MSChromatogram, MSSpectrum, Peak1D};

fn chromatogram(points: &[(f64, f32)]) -> MSChromatogram {
    MSChromatogram {
        peaks: points
            .iter()
            .map(|&(x, y)| ChromatogramPeak::new(x, y))
            .collect(),
        ..Default::default()
    }
}
fn integrator(method: IntegrationMethod) -> PeakIntegrator {
    PeakIntegrator {
        integration_method: method,
        ..Default::default()
    }
}
fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 1e-10, "{actual} != {expected}");
}

#[test]
fn source_defaults_names_and_legacy_baseline_alias() {
    let default = PeakIntegrator::default();
    assert_eq!(default.integration_method, IntegrationMethod::IntensitySum);
    assert_eq!(default.baseline_type, BaselineType::BaseToBase);
    for method in [
        IntegrationMethod::IntensitySum,
        IntegrationMethod::Trapezoid,
        IntegrationMethod::Simpson,
    ] {
        assert_eq!(method.name().parse::<IntegrationMethod>().unwrap(), method);
        assert_eq!(method.to_string(), method.name());
    }
    for baseline in [
        BaselineType::BaseToBase,
        BaselineType::VerticalDivisionMin,
        BaselineType::VerticalDivisionMax,
    ] {
        assert_eq!(baseline.name().parse::<BaselineType>().unwrap(), baseline);
        assert_eq!(baseline.to_string(), baseline.name());
    }
    assert_eq!(
        "vertical_division".parse::<BaselineType>().unwrap(),
        BaselineType::VerticalDivisionMin
    );
    assert!("trapezoidal".parse::<IntegrationMethod>().is_err());
    assert!("unknown".parse::<BaselineType>().is_err());
}

#[test]
fn bounds_are_inclusive_samples_with_first_positive_maximum() {
    let input = chromatogram(&[(0.0, 99.0), (1.0, 2.0), (2.0, 7.0), (3.0, 7.0), (4.0, 99.0)]);
    let before = input.clone();
    let algorithm = PeakIntegrator::default();
    let output = algorithm.integrate_chromatogram(&input, 0.1, 3.9).unwrap();
    assert_eq!(output.area, 16.0);
    assert_eq!(output.height, 7.0);
    assert_eq!(output.apex_pos, 2.0);
    assert_eq!(output.hull_points, vec![[1.0, 2.0], [2.0, 7.0], [3.0, 7.0]]);
    assert_eq!(
        output,
        algorithm.integrate_chromatogram(&input, 1.0, 3.0).unwrap()
    );
    let spec = MSSpectrum {
        peaks: input
            .peaks
            .iter()
            .map(|p| Peak1D::new(p.rt, p.intensity))
            .collect(),
        ..Default::default()
    };
    assert_eq!(
        output,
        algorithm.integrate_spectrum(&spec, 1.0, 3.0).unwrap()
    );
    assert_eq!(input, before);
}

#[test]
fn empty_singleton_and_nonpositive_peak_conventions() {
    for method in [
        IntegrationMethod::IntensitySum,
        IntegrationMethod::Trapezoid,
        IntegrationMethod::Simpson,
    ] {
        let algorithm = integrator(method);
        let empty = algorithm
            .integrate_chromatogram(&MSChromatogram::default(), 2.0, 4.0)
            .unwrap();
        assert_eq!((empty.area, empty.height, empty.apex_pos), (0.0, 0.0, 3.0));
        assert!(empty.hull_points.is_empty());
        let input = chromatogram(&[(3.0, -7.0)]);
        let output = algorithm.integrate_chromatogram(&input, 0.0, 8.0).unwrap();
        assert_eq!(
            output.area,
            if method == IntegrationMethod::IntensitySum {
                -7.0
            } else {
                0.0
            }
        );
        assert_eq!((output.height, output.apex_pos), (0.0, 4.0));
        let positive = algorithm
            .integrate_chromatogram(&chromatogram(&[(3.0, 7.0)]), 0.0, 8.0)
            .unwrap();
        assert_eq!((positive.height, positive.apex_pos), (7.0, 3.0));
    }
    let input = chromatogram(&[(0.0, -2.0), (2.0, -4.0)]);
    for method in [IntegrationMethod::Trapezoid, IntegrationMethod::Simpson] {
        assert_eq!(
            integrator(method)
                .integrate_chromatogram(&input, 0.0, 2.0)
                .unwrap()
                .area,
            -6.0
        );
    }
}

#[test]
fn trapezoid_preserves_float_intensity_addition_and_two_point_simpson_fallback() {
    let input = chromatogram(&[(0.0, 16_777_216.0), (1.0, 1.0)]);
    for method in [IntegrationMethod::Trapezoid, IntegrationMethod::Simpson] {
        assert_eq!(
            integrator(method)
                .integrate_chromatogram(&input, 0.0, 1.0)
                .unwrap()
                .area,
            8_388_608.0
        );
    }
    let input = chromatogram(&[(0.0, f32::MAX), (1.0, f32::MAX)]);
    assert!(
        integrator(IntegrationMethod::Trapezoid)
            .integrate_chromatogram(&input, 0.0, 1.0)
            .is_err()
    );
    assert!(
        PeakIntegrator::default()
            .integrate_chromatogram(&input, 0.0, 1.0)
            .unwrap()
            .area
            .is_finite()
    );
}

#[test]
fn nonuniform_simpson_integrates_a_quadratic_without_clamping_negative_results() {
    let input = chromatogram(&[(0.0, 0.0), (1.0, 1.0), (3.0, 9.0)]);
    close(
        integrator(IntegrationMethod::Simpson)
            .integrate_chromatogram(&input, 0.0, 3.0)
            .unwrap()
            .area,
        9.0,
    );
    // Algebraic negative three-point integral: h=1,k=3,left intensity=1.5.
    let input = chromatogram(&[(0.0, 1.5), (1.0, 0.0), (4.0, 0.0)]);
    assert_eq!(
        integrator(IntegrationMethod::Simpson)
            .integrate_chromatogram(&input, 0.0, 4.0)
            .unwrap()
            .area,
        -1.0
    );
}

#[test]
fn all_baseline_methods_follow_sampling_and_endpoint_direction() {
    // The straight baseline is y=2+2x, observed at nonuniform x=0,1,3.
    let input = chromatogram(&[(0.0, 2.0), (1.0, 20.0), (3.0, 8.0)]);
    let falling = chromatogram(&[(0.0, 8.0), (1.0, 20.0), (3.0, 2.0)]);
    for method in [
        IntegrationMethod::IntensitySum,
        IntegrationMethod::Trapezoid,
        IntegrationMethod::Simpson,
    ] {
        for (baseline, height, sum_area, weighted_area) in [
            (BaselineType::BaseToBase, 4.0, 14.0, 15.0),
            (BaselineType::VerticalDivisionMin, 2.0, 6.0, 6.0),
            (BaselineType::VerticalDivisionMax, 8.0, 24.0, 24.0),
        ] {
            let algorithm = PeakIntegrator {
                integration_method: method,
                baseline_type: baseline,
                ..Default::default()
            };
            let result = algorithm
                .estimate_background_chromatogram(&input, -1.0, 4.0, 1.0)
                .unwrap();
            assert_eq!(result.height, height);
            assert_eq!(
                result.area,
                if method == IntegrationMethod::IntensitySum {
                    sum_area
                } else {
                    weighted_area
                }
            );
        }
        let algorithm = integrator(method);
        let result = algorithm
            .estimate_background_chromatogram(&falling, 0.0, 3.0, 1.0)
            .unwrap();
        assert_eq!(result.height, 6.0);
        assert_eq!(
            result.area,
            if method == IntegrationMethod::IntensitySum {
                16.0
            } else {
                15.0
            }
        );
    }
}

#[test]
fn baseline_empty_and_zero_span_are_explicit() {
    let empty = MSChromatogram::default();
    let singleton = chromatogram(&[(1.0, 4.0)]);
    let duplicates = chromatogram(&[(1.0, 4.0), (1.0, 6.0)]);
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
            let algorithm = PeakIntegrator {
                integration_method: method,
                baseline_type: baseline,
                ..Default::default()
            };
            assert!(
                algorithm
                    .estimate_background_chromatogram(&empty, 0.0, 2.0, 1.0)
                    .is_err()
            );
            for input in [&singleton, &duplicates] {
                let result = algorithm.estimate_background_chromatogram(input, 0.0, 2.0, 1.0);
                if baseline == BaselineType::BaseToBase {
                    assert!(result.is_err());
                } else {
                    let result = result.unwrap();
                    assert_eq!(
                        result.area,
                        if method == IntegrationMethod::IntensitySum {
                            result.height * input.peaks.len() as f64
                        } else {
                            0.0
                        }
                    );
                }
            }
        }
    }
}

#[test]
fn sampled_shape_thresholds_are_inclusive_and_clipped_without_interpolation() {
    let input = chromatogram(&[
        (0.0, 0.0),
        (1.0, 5.0),
        (2.0, 10.0),
        (3.0, 50.0),
        (4.0, 100.0),
        (5.0, 50.0),
        (6.0, 10.0),
        (7.0, 5.0),
        (8.0, 0.0),
    ]);
    let algorithm = PeakIntegrator::default();
    let shape = algorithm
        .calculate_shape_metrics_chromatogram(&input, 0.0, 8.0, 100.0, 4.0)
        .unwrap();
    assert_eq!(
        (shape.width_at_5, shape.width_at_10, shape.width_at_50),
        (6.0, 4.0, 2.0)
    );
    assert_eq!(
        (
            shape.start_position_at_5,
            shape.start_position_at_10,
            shape.start_position_at_50
        ),
        (1.0, 2.0, 3.0)
    );
    assert_eq!(
        (
            shape.end_position_at_5,
            shape.end_position_at_10,
            shape.end_position_at_50
        ),
        (7.0, 6.0, 5.0)
    );
    assert_eq!((shape.tailing_factor, shape.asymmetry_factor), (1.0, 1.0));
    assert_eq!(
        (
            shape.points_across_baseline,
            shape.points_across_half_height
        ),
        (9, 3)
    );
    let clipped = algorithm
        .calculate_shape_metrics_chromatogram(&input, 3.0, 5.0, 100.0, 4.0)
        .unwrap();
    assert_eq!(
        (clipped.width_at_5, clipped.width_at_10, clipped.width_at_50),
        (2.0, 2.0, 2.0)
    );
    let spectrum = MSSpectrum {
        peaks: input
            .peaks
            .iter()
            .map(|p| Peak1D::new(p.rt, p.intensity))
            .collect(),
        ..Default::default()
    };
    assert_eq!(
        shape,
        algorithm
            .calculate_shape_metrics_spectrum(&spectrum, 0.0, 8.0, 100.0, 4.0)
            .unwrap()
    );
    assert_eq!(
        algorithm
            .estimate_background_spectrum(&spectrum, 0.0, 8.0, 4.0)
            .unwrap(),
        algorithm
            .estimate_background_chromatogram(&input, 0.0, 8.0, 4.0)
            .unwrap()
    );
}

#[test]
fn shape_zero_height_singleton_empty_and_endpoint_apex_guards() {
    let algorithm = PeakIntegrator::default();
    assert_eq!(
        algorithm
            .calculate_shape_metrics_chromatogram(&MSChromatogram::default(), 0.0, 1.0, 0.0, 0.5)
            .unwrap(),
        PeakShapeMetrics::default()
    );
    let input = chromatogram(&[(2.0, 0.0)]);
    let shape = algorithm
        .calculate_shape_metrics_chromatogram(&input, 0.0, 3.0, 0.0, 2.0)
        .unwrap();
    assert_eq!(
        (
            shape.total_width,
            shape.tailing_factor,
            shape.asymmetry_factor,
            shape.baseline_delta_2_height
        ),
        (0.0, 0.0, 0.0, 0.0)
    );
    assert_eq!(
        (shape.start_position_at_5, shape.end_position_at_5),
        (2.0, 2.0)
    );
    assert_eq!(
        (
            shape.points_across_baseline,
            shape.points_across_half_height
        ),
        (1, 1)
    );
    assert!(
        algorithm
            .calculate_shape_metrics_chromatogram(&input, 0.0, 3.0, 0.0, 2.5)
            .is_err()
    );
    assert!(
        algorithm
            .calculate_shape_metrics_chromatogram(&input, 0.0, 1.0, 0.0, 0.5)
            .is_err()
    );
    for height in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(
            algorithm
                .calculate_shape_metrics_chromatogram(&input, 0.0, 3.0, height, 2.0)
                .is_err()
        );
    }
}

#[test]
fn duplicate_spacing_and_nonfinite_calculation_errors_preserve_inputs() {
    let input = chromatogram(&[(0.0, 1.0), (0.0, 2.0), (1.0, 3.0)]);
    let before = input.clone();
    assert_eq!(
        PeakIntegrator::default()
            .integrate_chromatogram(&input, 0.0, 1.0)
            .unwrap()
            .area,
        6.0
    );
    assert_eq!(
        integrator(IntegrationMethod::Trapezoid)
            .integrate_chromatogram(&input, 0.0, 1.0)
            .unwrap()
            .area,
        2.5
    );
    assert!(
        integrator(IntegrationMethod::Simpson)
            .integrate_chromatogram(&input, 0.0, 1.0)
            .is_err()
    );
    assert_eq!(input, before);
    let tiny = chromatogram(&[(0.0, 1.0), (1e-200, 1.0), (2e-200, 1.0)]);
    assert!(
        integrator(IntegrationMethod::Simpson)
            .integrate_chromatogram(&tiny, 0.0, 2e-200)
            .is_err()
    );
    let huge = chromatogram(&[(-f64::MAX, 1.0), (f64::MAX, 1.0)]);
    assert!(
        integrator(IntegrationMethod::Trapezoid)
            .integrate_chromatogram(&huge, -f64::MAX, f64::MAX)
            .is_err()
    );
}

#[test]
fn preflight_rejects_unsorted_nonfinite_invalid_bounds_and_whole_input_limit() {
    let algorithm = PeakIntegrator::default();
    let input = chromatogram(&[(2.0, 1.0), (1.0, 1.0)]);
    assert!(matches!(
        algorithm.integrate_chromatogram(&input, 0.0, 3.0),
        Err(Error::UnsortedData)
    ));
    for input in [
        chromatogram(&[(f64::NAN, 1.0)]),
        chromatogram(&[(0.0, f32::INFINITY)]),
    ] {
        assert!(algorithm.integrate_chromatogram(&input, 0.0, 3.0).is_err());
    }
    let input = chromatogram(&[(0.0, 1.0), (1.0, 2.0)]);
    for (left, right) in [(2.0, 1.0), (f64::NAN, 2.0), (0.0, f64::INFINITY)] {
        assert!(
            algorithm
                .integrate_chromatogram(&input, left, right)
                .is_err()
        );
    }
    for max_points in [0, 1] {
        let limited = PeakIntegrator {
            max_points,
            ..Default::default()
        };
        assert!(limited.integrate_chromatogram(&input, 0.0, 0.0).is_err());
        assert!(
            limited
                .estimate_background_chromatogram(&input, 0.0, 1.0, 0.0)
                .is_err()
        );
        assert!(
            limited
                .calculate_shape_metrics_chromatogram(&input, 0.0, 1.0, 2.0, 1.0)
                .is_err()
        );
    }
    assert!(
        algorithm
            .estimate_background_chromatogram(&input, 0.0, 1.0, 2.0)
            .is_err()
    );
    assert!(
        algorithm
            .calculate_shape_metrics_chromatogram(&input, 0.0, 1.0, 2.0, f64::NAN)
            .is_err()
    );
}

#[test]
fn shape_rejects_overflowing_denominator_even_when_division_would_return_zero() {
    let input = chromatogram(&[(-1e308, 0.0), (0.0, 100.0), (6e307, 0.0)]);
    assert!(
        PeakIntegrator::default()
            .calculate_shape_metrics_chromatogram(&input, -1e308, 6e307, 100.0, 0.0)
            .is_err()
    );
}
