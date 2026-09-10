// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::analysis::transformations::*;
fn data(points: &[(f64, f64)]) -> Vec<DataPoint> {
    points.iter().copied().map(DataPoint::from).collect()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}, tolerance {tolerance}"
    );
}
fn table(text: &str) -> Vec<Vec<f64>> {
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split_whitespace().map(|v| v.parse().unwrap()).collect())
        .collect()
}
#[test]
fn upstream_linear_fits_weights_and_explicit_inverse() {
    let points = data(&[(0., 1.), (1., 2.), (1., 4.)]);
    let model = LinearModel::fit(&points, LinearOptions::default()).unwrap();
    for x in [-0.5, 0., 0.5, 1., 1.5] {
        close(model.apply(x).unwrap(), 2. * x + 1., 1e-14);
    }
    let one = LinearModel::fit(&data(&[(4., 8.)]), Default::default()).unwrap();
    assert_eq!(
        one.coefficients(),
        LinearCoefficients {
            slope: 1.,
            intercept: 4.
        }
    );
    let two = LinearModel::fit(&data(&[(0., 1.), (1e-20, 2.)]), Default::default()).unwrap();
    close(two.apply(1e-20).unwrap(), 2., 1e-15);
    let log = CoordinateWeight {
        function: WeightFunction::Log,
        ..Default::default()
    };
    let weighted = LinearModel::fit(
        &data(&[(1., 2.), (2., 4.), (4., 8.)]),
        LinearOptions {
            x_weight: log,
            y_weight: log,
            ..Default::default()
        },
    )
    .unwrap();
    close(weighted.apply(2.).unwrap(), 4., 1e-14);
    close(weighted.inverse().unwrap().apply(4.).unwrap(), 2., 1e-14);
    let reciprocal = CoordinateWeight {
        function: WeightFunction::Reciprocal,
        ..Default::default()
    };
    let weighted = LinearModel::fit(
        &data(&[(1., 2.), (2., 4.), (4., 8.)]),
        LinearOptions {
            x_weight: reciprocal,
            ..Default::default()
        },
    )
    .unwrap();
    close(weighted.apply(2.).unwrap(), 5.285714286, 1e-9);
    close(
        weighted.inverse().unwrap().apply(5.285714286).unwrap(),
        2.,
        1e-9,
    );
    let explicit = LinearModel::from_coefficients(12.3, -45.6).unwrap();
    let inverse = explicit.inverse().unwrap();
    for x in [-1., 0., 10.] {
        close(inverse.apply(explicit.apply(x).unwrap()).unwrap(), x, 1e-14);
    }
    let points = data(&[(0., 1.), (1., 2.), (1., 4.), (2., 2.)]);
    let fit = LinearModel::fit(
        &points,
        LinearOptions {
            x_weight: CoordinateWeight {
                function: WeightFunction::Log,
                min: 1e-4,
                max: 1e15,
            },
            y_weight: CoordinateWeight {
                function: WeightFunction::Log,
                min: 1e-7,
                max: 1e15,
            },
            ..Default::default()
        },
    )
    .unwrap();
    close(fit.coefficients().slope, 0.095_036_911_971_605_03, 1e-15);
    close(fit.coefficients().intercept, 0.895_509_115_454_389_9, 1e-15);
}
#[test]
fn upstream_interpolation_table_linear_and_natural_cubic() {
    let input = table(include_str!("data/transformations_interpolation_input.tsv"));
    let points: Vec<_> = input.iter().map(|r| DataPoint::new(r[0], r[1])).collect();
    let linear = InterpolatedModel::fit(
        &points,
        InterpolationOptions {
            interpolation: Interpolation::Linear,
            ..Default::default()
        },
    )
    .unwrap();
    let cubic = InterpolatedModel::fit(&points, Default::default()).unwrap();
    for row in table(include_str!(
        "data/transformations_interpolation_golden.tsv"
    )) {
        close(linear.apply(row[0]).unwrap(), row[1], 1e-5);
        close(cubic.apply(row[0]).unwrap(), row[3], 1e-5);
    }
}
#[test]
fn upstream_duplicate_x_averaging_and_all_extrapolation_modes() {
    let points = data(&[(0., 1.), (0.5, 4.), (1., 2.), (1., 4.)]);
    for (mode, left, right) in [
        (Extrapolation::TwoPointLinear, 0., 4.),
        (Extrapolation::FourPointLinear, -2., 2.),
        (
            Extrapolation::GlobalLinear,
            0.909090909090909,
            4.181818181818182,
        ),
    ] {
        let model = InterpolatedModel::fit(
            &points,
            InterpolationOptions {
                extrapolation: mode,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            model.knots().collect::<Vec<_>>(),
            [(0., 1.), (0.5, 4.), (1., 3.)]
        );
        close(model.apply(-0.5).unwrap(), left, 1e-14);
        close(model.apply(1.5).unwrap(), right, 1e-14);
        assert_eq!(model.apply(1.).unwrap(), 3.);
    }
}
#[test]
fn upstream_lowess_sine_and_extrapolation_goldens() {
    let points: Vec<_> = table(include_str!("data/transformations_lowess_input.tsv"))
        .iter()
        .map(|r| DataPoint::new(r[0], r[1]))
        .collect();
    let model = LowessModel::fit(
        &points,
        LowessOptions {
            span: 0.3,
            ..Default::default()
        },
    )
    .unwrap();
    for row in table(include_str!("data/transformations_lowess_golden.tsv")) {
        close(model.apply(row[0]).unwrap(), row[1], 1e-6);
    }
    for (mode, left, right) in [
        (
            Extrapolation::FourPointLinear,
            0.815490292172986,
            -0.571905836956494,
        ),
        (Extrapolation::TwoPointLinear, -0.04240732863, 0.046870277),
        (Extrapolation::GlobalLinear, -0.9501004, 1.08486397),
    ] {
        let model = LowessModel::fit(
            &points,
            LowessOptions {
                span: 0.3,
                extrapolation: mode,
                ..Default::default()
            },
        )
        .unwrap();
        close(model.apply(-4.).unwrap(), left, 1e-7);
        close(model.apply(4.).unwrap(), right, 1e-7);
    }
}
#[test]
fn upstream_lowess_cars_robust_fit_and_tied_coordinates() {
    let rows = table(include_str!("data/transformations_lowess_cars.tsv"));
    let x: Vec<_> = rows.iter().map(|r| r[0]).collect();
    let y: Vec<_> = rows.iter().map(|r| r[1]).collect();
    for (span, column) in [(2. / 3., 2), (0.2, 3)] {
        let actual = lowess(
            &x,
            &y,
            LowessOptions {
                span,
                ..Default::default()
            },
        )
        .unwrap();
        for (value, row) in actual.iter().zip(&rows) {
            close(*value, row[column], 1e-10);
        }
    }
}
#[test]
fn description_identity_lock_notes_reset_and_source_deviations() {
    let mut points = data(&[(0., 1.), (0.25, 1.0625), (0.5, 1.25), (1., 2.)]);
    points[1].note = "peptide A".into();
    let mut description = TransformationDescription::new(points.clone()).unwrap();
    description
        .fit_model(ModelConfig::Linear(Default::default()))
        .unwrap();
    assert_eq!(
        description.deviations(false, true).unwrap(),
        [0.75, 0.8125, 1., 1.]
    );
    for (a, b) in description.deviations(true, false).unwrap().iter().zip([
        0.125,
        0.0714285714285714,
        0.1428571428571428,
        0.0892857142857142,
    ]) {
        close(*a, b, 1e-14);
    }
    let stats = description.statistics().unwrap();
    assert_eq!(
        stats
            .percentiles
            .iter()
            .map(|p| p.percent)
            .collect::<Vec<_>>(),
        [100, 99, 95, 90, 75, 50, 25]
    );
    close(stats.percentiles[0].after, 0.1428571428571428, 1e-14);
    close(stats.percentiles[1].after, 0.125, 1e-14);
    description.invert().unwrap();
    assert_eq!(description.data_points()[1].note, "peptide A");
    description.invert().unwrap();
    assert_eq!(description.data_points(), points);
    description.fit_model(ModelConfig::Identity).unwrap();
    description
        .fit_model(ModelConfig::Linear(Default::default()))
        .unwrap();
    assert_eq!(description.model().name(), "identity");
    assert_eq!(description.apply(8.).unwrap(), 8.);
    description.set_data_points(points).unwrap();
    assert_eq!(description.model().name(), "none");
}

#[test]
fn upstream_original_lowess_delta_and_robust_pass_goldens() {
    let rows = table(include_str!("data/transformations_lowess_original.tsv"));
    let x: Vec<_> = rows.iter().map(|r| r[0]).collect();
    let y: Vec<_> = rows.iter().map(|r| r[1]).collect();
    for (iterations, delta, column) in [(0, 0., 2), (0, 3., 3), (2, 0., 4)] {
        let actual = lowess(
            &x,
            &y,
            LowessOptions {
                span: 0.25,
                iterations,
                delta: Some(delta),
                ..Default::default()
            },
        )
        .unwrap();
        for (value, row) in actual.iter().zip(&rows) {
            close(*value, row[column], 1e-3);
        }
    }
}

#[test]
fn inverse_description_refits_but_explicit_linear_inverse_is_algebraic() {
    let mut description =
        TransformationDescription::new(data(&[(0., 0.), (1., 1.), (2., 3.)])).unwrap();
    description
        .fit_model(ModelConfig::Linear(Default::default()))
        .unwrap();
    let inverse = description.inverse().unwrap();
    // Inverse least squares is not the reciprocal of the forward least squares line.
    close(inverse.apply(1.).unwrap(), 11. / 14., 1e-14);
    let TransformationModel::Linear(model) = description.model() else {
        panic!()
    };
    close(model.inverse().unwrap().apply(1.).unwrap(), 7. / 9., 1e-14);
    let mut explicit = TransformationDescription::default();
    explicit
        .fit_model(ModelConfig::Linear(LinearOptions {
            coefficients: Some(LinearCoefficients {
                slope: 3.,
                intercept: -4.,
            }),
            ..Default::default()
        }))
        .unwrap();
    explicit.invert().unwrap();
    close(explicit.apply(2.).unwrap(), 2., 1e-15);
    let mut constant =
        TransformationDescription::new(data(&[(0., 1.), (1., 1.), (2., 1.)])).unwrap();
    constant
        .fit_model(ModelConfig::Linear(Default::default()))
        .unwrap();
    let anchors = constant.data_points().to_vec();
    assert!(constant.invert().is_err());
    assert_eq!(constant.data_points(), anchors);
    assert_eq!(constant.apply(9.).unwrap(), 1.);
    assert!(
        LinearModel::from_coefficients(0., 2.)
            .unwrap()
            .inverse()
            .is_err()
    );
}

#[test]
fn adaptive_residual_windows_sparse_dense_small_and_inverse_units() {
    let mut points: Vec<_> = (0..99).map(|y| DataPoint::new(0., f64::from(y))).collect();
    points.push(DataPoint::new(0., 100_000.));
    let description = TransformationDescription::new(points).unwrap();
    let options = WindowOptions {
        quantile: 1.,
        inverse: false,
        full_window: false,
        padding_factor: 1.,
    };
    close(description.estimate_window(options).unwrap(), 148.5, 1e-14);
    close(
        description
            .estimate_window(WindowOptions {
                full_window: true,
                padding_factor: 1.5,
                ..options
            })
            .unwrap(),
        445.5,
        1e-14,
    );
    let mut points: Vec<_> = (0..90).map(|y| DataPoint::new(0., f64::from(y))).collect();
    points.extend((0..10).map(|_| DataPoint::new(0., 10_000.)));
    assert_eq!(
        TransformationDescription::new(points)
            .unwrap()
            .estimate_window(options)
            .unwrap(),
        10_000.
    );
    assert_eq!(
        TransformationDescription::new(data(&[(0., 1.), (0., 100.)]))
            .unwrap()
            .estimate_window(options)
            .unwrap(),
        100.
    );
    let mut description =
        TransformationDescription::new(data(&[(0., 0.), (1., 1.), (2., 3.)])).unwrap();
    description
        .fit_model(ModelConfig::Linear(Default::default()))
        .unwrap();
    close(
        description.estimate_window(options).unwrap(),
        1. / 3.,
        1e-14,
    );
    close(
        description
            .estimate_window(WindowOptions {
                inverse: true,
                ..options
            })
            .unwrap(),
        3. / 14.,
        1e-14,
    );
}

#[test]
fn empty_and_single_anchor_statistics_have_checked_small_sample_behavior() {
    let empty = TransformationDescription::default();
    assert_eq!(
        empty.statistics().unwrap(),
        TransformationStatistics::default()
    );
    assert_eq!(empty.estimate_window(Default::default()).unwrap(), 0.);
    assert_eq!(empty.apply(-4.).unwrap(), -4.);
    let description =
        TransformationDescription::new(vec![DataPoint::with_note(2., 5., "one")]).unwrap();
    let stats = description.statistics().unwrap();
    assert_eq!(stats.x_range.unwrap().min, 2.);
    assert_eq!(stats.y_range.unwrap().max, 5.);
    assert!(
        stats
            .percentiles
            .iter()
            .all(|p| p.before == 3. && p.after == 3.)
    );
}

#[test]
fn transformation_errors_and_atomic_batch_fit_and_data_updates() {
    let weight = CoordinateWeight {
        function: WeightFunction::Log,
        ..Default::default()
    };
    close(
        weight.transform_training(0.).unwrap(),
        1e-15_f64.ln(),
        1e-14,
    );
    assert!(weight.transform(0.).is_err());
    assert!(
        CoordinateWeight {
            function: WeightFunction::ReciprocalSquared,
            ..Default::default()
        }
        .transform(f64::MAX)
        .is_err()
    );
    assert!(
        CoordinateWeight {
            min: 2.,
            max: 1.,
            ..Default::default()
        }
        .transform(1.)
        .is_err()
    );
    let mut description = TransformationDescription::default();
    description
        .fit_model(ModelConfig::Linear(LinearOptions {
            x_weight: weight,
            coefficients: Some(LinearCoefficients {
                slope: 1.,
                intercept: 0.,
            }),
            ..Default::default()
        }))
        .unwrap();
    let mut batch = [1., 0., 2.];
    assert!(description.apply_values(&mut batch).is_err());
    assert_eq!(batch, [1., 0., 2.]);
    assert!(
        description
            .fit_model(ModelConfig::Linear(Default::default()))
            .is_err()
    );
    assert_eq!(description.apply(1.).unwrap(), 0.);
    assert!(
        description
            .set_data_points(data(&[(f64::NAN, 2.)]))
            .is_err()
    );
    assert!(description.data_points().is_empty());
    assert_eq!(description.apply(1.).unwrap(), 0.);
    assert!(TransformationDescription::new(data(&[(1., f64::INFINITY)])).is_err());
    assert!(LinearModel::fit(&data(&[(1., 2.), (1., 4.)]), Default::default()).is_err());
    assert!(
        LinearModel::fit(
            &data(&[(-f64::MAX, 0.), (f64::MAX, 1.)]),
            Default::default()
        )
        .is_err()
    );
    assert!(
        LinearModel::from_coefficients(f64::MAX, 1.)
            .unwrap()
            .apply(2.)
            .is_err()
    );
    assert!(
        InterpolatedModel::fit(&data(&[(0., 1.), (1., 2.), (1., 3.)]), Default::default()).is_err()
    );
    assert!(
        InterpolatedModel::fit(&data(&[(0., 1.), (1., 2.), (2., 3.)]), Default::default())
            .unwrap()
            .apply(f64::NAN)
            .is_err()
    );
    for options in [
        WindowOptions {
            quantile: 0.,
            ..Default::default()
        },
        WindowOptions {
            quantile: f64::NAN,
            ..Default::default()
        },
        WindowOptions {
            padding_factor: -1.,
            ..Default::default()
        },
    ] {
        assert!(description.estimate_window(options).is_err());
    }
}

#[test]
fn lowess_limits_sorting_invalid_parameters_and_degenerate_inputs() {
    let x = [0., 1., 2., 3.];
    let y = [2., 3., 5., 7.];
    for options in [
        LowessOptions {
            span: 0.,
            ..Default::default()
        },
        LowessOptions {
            span: 1.1,
            ..Default::default()
        },
        LowessOptions {
            iterations: 65,
            ..Default::default()
        },
        LowessOptions {
            delta: Some(-1.),
            ..Default::default()
        },
        LowessOptions {
            max_work: 1,
            ..Default::default()
        },
        LowessOptions {
            max_points: 3,
            ..Default::default()
        },
    ] {
        assert!(lowess(&x, &y, options).is_err());
    }
    assert!(lowess(&[1., 0.], &[1., 2.], Default::default()).is_err());
    assert!(lowess(&[1., 2.], &[1.], Default::default()).is_err());
    assert!(lowess(&[1., 2.], &[1., f64::NAN], Default::default()).is_err());
    assert!(
        lowess(
            &[f64::MAX / 2., f64::MAX],
            &[1., 2.],
            LowessOptions {
                delta: Some(f64::MAX),
                ..Default::default()
            }
        )
        .is_err()
    );
    let two = lowess(&[1., 2.], &[3., 5.], Default::default()).unwrap();
    assert_eq!(two, [3., 5.]);
    assert!(LowessModel::fit(&data(&[(1., 3.), (2., 5.)]), Default::default()).is_err());
    let repeated = lowess(&[1., 1., 1.], &[1., 2., 3.], Default::default()).unwrap();
    assert_eq!(repeated, [2., 2., 2.]);
    let mut limited = TransformationDescription::with_max_points(data(&[(0., 1.)]), 1).unwrap();
    assert!(
        limited
            .set_data_points(data(&[(0., 1.), (1., 2.)]))
            .is_err()
    );
    assert!(limited.apply_values(&mut [0., 1.]).is_err());
    assert!(
        LinearModel::fit(
            &data(&[(0., 1.), (1., 2.)]),
            LinearOptions {
                max_points: 1,
                ..Default::default()
            }
        )
        .is_err()
    );
}
