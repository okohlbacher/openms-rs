// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::kernel::{ChromatogramPeak, DataArray, MSChromatogram, MSSpectrum, Peak1D};
use openms::processing::mean_noise::{
    MeanNoiseHistogramRange as Range, SignalToNoiseEstimatorMeanIterative as Estimator,
};

fn points(data: &str) -> MSSpectrum {
    MSSpectrum::from_peaks(
        data.lines()
            .skip(1)
            .map(|line| {
                let mut fields = line.split_whitespace();
                Peak1D::new(
                    fields.next().unwrap().parse().unwrap(),
                    fields.next().unwrap().parse().unwrap(),
                )
            })
            .collect(),
    )
}
fn small() -> Estimator {
    Estimator {
        histogram_range: Range::Manual {
            max_intensity: 30.0,
        },
        bin_count: 3,
        min_required_elements: 1,
        ..Default::default()
    }
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-12, "{a} != {b}");
}

#[test]
fn complete_historical_source_fixture_matches_every_signal_to_noise_value() {
    let input = points(include_str!("data/peak_picking_noise_input.dta"));
    let expected = points(include_str!("data/mean_noise_historical.dta"));
    assert_eq!(input.len(), 2526);
    assert_eq!(input.len(), expected.len());
    let estimator = Estimator {
        window_length: 40.1,
        noise_for_empty_window: 2.0,
        ..Default::default()
    };
    let actual = estimator.estimate_spectrum(&input).unwrap();
    for ((peak, expected), actual) in input
        .peaks
        .iter()
        .zip(expected.peaks)
        .zip(actual.signal_to_noise)
    {
        assert_eq!(peak.mz, expected.mz);
        // Exact source class-test tolerance; independent analytical cases below
        // exercise operation details more tightly than this historical fixture.
        assert!(
            (actual - f64::from(expected.intensity)).abs() <= 0.5,
            "m/z {}: {actual} != {}",
            peak.mz,
            expected.intensity
        );
    }
}

#[test]
fn clipping_keeps_original_count_and_uses_third_pass_mean() {
    let estimator = Estimator {
        stdev_multiplier: 0.5,
        ..small()
    };
    let result = estimator.estimate(&[0., 1., 2.], &[4., 12., 28.]).unwrap();
    // Histogram centers5,15,25; passes retain3,2,1 bins. Every mean still
    // divides by3, so the third mean is5/3, not a renormalized5.
    for noise in result.noise {
        close(noise, 5.0 / 3.0);
    }
    for (actual, expected) in result.signal_to_noise.into_iter().zip([2.4, 7.2, 16.8]) {
        close(actual, expected);
    }
    assert_eq!(result.sparse_window_percent, 0.0);
}

#[test]
fn windows_are_left_inclusive_right_exclusive_and_exclude_above_histogram_range() {
    let estimator = Estimator {
        window_length: 2.0,
        min_required_elements: 2,
        noise_for_empty_window: 2.0,
        ..small()
    };
    let result = estimator
        .estimate(&[0., 1., 2., 3.], &[2., 8., 14., 20.])
        .unwrap();
    assert_eq!(result.noise, [2., 5., 10., 20.]);
    assert_eq!(result.sparse_window_percent, 25.0);
    // Effective width is at least1: with3bins and max2, intensities in[2,3)
    // are retained even though they exceed the configured maximum.
    let estimator = Estimator {
        histogram_range: Range::Manual { max_intensity: 2.0 },
        noise_for_empty_window: 17.0,
        min_required_elements: 2,
        ..small()
    };
    let result = estimator.estimate(&[1., 1., 1.], &[2., 2.5, 3.]).unwrap();
    // Mean2.5 and zero deviation give a truncated rightmost-bin value2,
    // excluding the sole occupied bin on later passes; the final floor is1.
    assert_eq!(result.noise, [1.; 3]);
    let estimator = Estimator {
        min_required_elements: 3,
        ..estimator
    };
    let result = estimator.estimate(&[1., 1., 1.], &[2., 2.5, 3.]).unwrap();
    assert_eq!(result.noise, [17.; 3]);
    assert_eq!(result.sparse_window_percent, 100.0);
}

#[test]
fn signed_intensity_empty_and_zero_signals_have_explicit_results() {
    let estimator = Estimator {
        stdev_multiplier: 0.01,
        histogram_range: Range::Manual { max_intensity: 3.0 },
        ..small()
    };
    let result = estimator.estimate(&[-1., 0., 1.], &[-5., 0., 2.]).unwrap();
    assert_eq!(result.noise, [1.; 3]);
    assert_eq!(result.signal_to_noise, [-5., 0., 2.]);
    let empty = Estimator::default().estimate(&[], &[]).unwrap();
    assert!(empty.noise.is_empty());
    assert_eq!(empty.max_intensity, 0.0);
    assert_eq!(empty.sparse_window_percent, 0.0);
    let result = Estimator {
        min_required_elements: 1,
        ..Default::default()
    }
    .estimate(&[0., 1., 2.], &[0.; 3])
    .unwrap();
    assert_eq!(result.noise, [1.; 3]);
    assert_eq!(result.signal_to_noise, [0.; 3]);
    assert!(
        Estimator::default()
            .estimate(&[0., 1.], &[-1., -1.])
            .is_err()
    );
}

#[test]
fn legacy_percentile_retains_defined_source_behavior_and_checks_unsafe_cases() {
    let estimator = Estimator {
        histogram_range: Range::LegacyPercentile { percentile: 95 },
        noise_for_empty_window: 4.0,
        ..Default::default()
    };
    // Source minimum100 gives f32bin_size1. Its bin scan stops after10
    // input points even though all histogram entries are in bin99.
    let result = estimator.estimate(&[0.; 10], &[100.; 10]).unwrap();
    assert_eq!(result.max_intensity, 9.5);
    assert_eq!(result.signal_to_noise, [25.; 10]);
    let result = estimator.estimate(&[0.; 100], &[100.; 100]).unwrap();
    assert_eq!(result.max_intensity, 99.5);
    assert_eq!(result.sparse_window_percent, 100.0);
    // The reversed maximum comparator picks1, causing the value100 to index
    // outside the source100-entry histogram. The port rejects it before access.
    assert!(estimator.estimate(&[0., 1.], &[1., 100.]).is_err());
    assert!(estimator.estimate(&[0., 1.], &[0., 0.]).is_err());
    assert!(estimator.estimate(&[], &[]).is_err());
    assert!(
        Estimator {
            histogram_range: Range::LegacyPercentile { percentile: 0 },
            ..estimator
        }
        .estimate(&[0.; 10], &[100.; 10])
        .is_err()
    );
}

#[test]
fn spectrum_and_chromatogram_conveniences_preserve_input_and_annotations() {
    let mut input = MSSpectrum::from_peaks(vec![
        Peak1D::new(10., 4.),
        Peak1D::new(11., 12.),
        Peak1D::new(12., 28.),
    ]);
    input.metadata.insert("sample".into(), "mean noise".into());
    input
        .integer_data_arrays
        .push(DataArray::new("indices", vec![0, 1, 2]));
    let original = input.clone();
    let chromatogram = MSChromatogram {
        peaks: input
            .peaks
            .iter()
            .map(|p| ChromatogramPeak::new(p.mz, p.intensity))
            .collect(),
        ..Default::default()
    };
    let estimator = small();
    assert_eq!(
        estimator.estimate_spectrum(&input).unwrap(),
        estimator.estimate_chromatogram(&chromatogram).unwrap()
    );
    assert_eq!(input, original);
    input.integer_data_arrays[0].data.pop();
    assert!(estimator.estimate_spectrum(&input).is_err());
}

#[test]
fn resource_limits_are_inclusive_and_numerical_errors_are_checked() {
    let estimator = Estimator {
        stdev_multiplier: 0.5,
        max_work: 45,
        ..small()
    };
    assert!(estimator.estimate(&[0., 1., 2.], &[4., 12., 28.]).is_ok());
    assert!(
        Estimator {
            max_work: 44,
            ..estimator.clone()
        }
        .estimate(&[0., 1., 2.], &[4., 12., 28.])
        .is_err()
    );
    assert!(
        Estimator {
            max_points: 2,
            ..estimator.clone()
        }
        .estimate(&[0., 1., 2.], &[4., 12., 28.])
        .is_err()
    );
    assert!(
        Estimator {
            max_bins: 2,
            ..estimator.clone()
        }
        .estimate(&[0., 1., 2.], &[4., 12., 28.])
        .is_err()
    );
    assert!(estimator.estimate(&[0., 1.], &[4.]).is_err());
    assert!(estimator.estimate(&[1., 0.], &[4., 12.]).is_err());
    assert!(estimator.estimate(&[0., f64::NAN], &[4., 12.]).is_err());
    assert!(estimator.estimate(&[0., 1.], &[4., f64::INFINITY]).is_err());
    assert!(
        Estimator::default()
            .estimate(&[0., 1.], &[1e308, 1e308])
            .is_err()
    );
    assert!(estimator.estimate(&[1e308], &[4.]).is_err());
    assert!(
        Estimator {
            noise_for_empty_window: 0.,
            ..estimator
        }
        .estimate(&[0.], &[4.])
        .is_err()
    );
}
