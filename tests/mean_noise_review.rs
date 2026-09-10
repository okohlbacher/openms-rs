// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Independent checks of pinned 7c029e8 SignalToNoiseEstimatorMeanIterative.h.
//! Rebuild each histogram from the whole trace to test the incremental window
//! implementation, then check source f32 percentile and signed global statistics.

use openms::processing::mean_noise::{
    MeanNoiseHistogramRange as Range, SignalToNoiseEstimatorMeanIterative as Estimator,
};

// Deliberately rescan all samples at every center: no sliding indices or updates.
fn rescan_noise(estimator: &Estimator, positions: &[f64], y: &[f64], max: f64) -> Vec<f64> {
    let width = (max / estimator.bin_count as f64).max(1.0);
    positions
        .iter()
        .map(|&center| {
            let members: Vec<_> = positions
                .iter()
                .zip(y)
                .filter(|(x, value)| {
                    **x >= center - estimator.window_length / 2.0
                        && **x < center + estimator.window_length / 2.0
                        && value.max(0.0) / width < estimator.bin_count as f64
                })
                .map(|(_, value)| (value.max(0.0) / width) as usize)
                .collect();
            if members.len() < estimator.min_required_elements {
                return estimator.noise_for_empty_window;
            }
            let histogram: Vec<_> = (0..estimator.bin_count)
                .map(|bin| members.iter().filter(|&&value| value == bin).count())
                .collect();
            let mut considered = estimator.bin_count;
            let mut mean = 0.0;
            for _ in 0..3 {
                mean = histogram[..considered]
                    .iter()
                    .enumerate()
                    .map(|(bin, &count)| {
                        count as f64 / members.len() as f64 * ((bin as f64 + 0.5) * width)
                    })
                    .sum();
                let variance: f64 = histogram[..considered]
                    .iter()
                    .enumerate()
                    .map(|(bin, &count)| {
                        let delta = (bin as f64 + 0.5) * width - mean;
                        count as f64 / members.len() as f64 * delta * delta
                    })
                    .sum();
                considered = (((mean + variance.sqrt() * estimator.stdev_multiplier - 1.0) / width
                    + 1.0) as usize)
                    .min(estimator.bin_count);
            }
            mean.max(1.0)
        })
        .collect()
}

#[test]
fn incremental_windows_match_independent_rescans_with_duplicates_gaps_and_signed_data() {
    let positions = [-4., -4., -3., -2., 0., 0., 0.5, 1., 2., 8., 9., 9., 10.];
    let y = [-3., 0., 0.5, 3., 2.5, 6., 15., 30., 29.5, 70., 1., 14., 7.];
    for max in [2.0, 30.0] {
        for bin_count in [3, 7] {
            for window_length in [1.0, 2.0, 6.0, 20.0] {
                for stdev_multiplier in [0.01, 0.5, 3.0] {
                    for min_required_elements in [1, 3] {
                        let estimator = Estimator {
                            histogram_range: Range::Manual { max_intensity: max },
                            bin_count,
                            window_length,
                            stdev_multiplier,
                            min_required_elements,
                            noise_for_empty_window: 19.0,
                            ..Default::default()
                        };
                        let expected = rescan_noise(&estimator, &positions, &y, max);
                        let actual = estimator.estimate(&positions, &y).unwrap();
                        assert_eq!(actual.noise, expected, "{estimator:?}");
                        for ((noise, observed), value) in
                            actual.noise.iter().zip(actual.signal_to_noise).zip(y)
                        {
                            assert_eq!(observed, value / noise);
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn global_standard_deviation_uses_signed_raw_intensities_before_histogram_clamping() {
    let estimator = Estimator {
        histogram_range: Range::StandardDeviation { factor: 1.0 },
        min_required_elements: 1,
        ..Default::default()
    };
    let actual = estimator.estimate(&[0., 1., 2.], &[-10., 0., 10.]).unwrap();
    // Signed population mean is zero, variance is (100+0+100)/3. Clamping
    // negatives before calculating this range would produce a different result.
    assert_eq!(actual.max_intensity, (200.0_f64 / 3.0).sqrt());
    assert!(actual.signal_to_noise[0] < 0.0);
    assert_eq!(actual.signal_to_noise[1], 0.0);
}

#[test]
fn historical_percentile_uses_f32_division_and_truncates_negative_fraction_to_zero() {
    let estimator = Estimator {
        histogram_range: Range::LegacyPercentile { percentile: 95 },
        ..Default::default()
    };
    let value = f32::from_bits(1.0_f32.to_bits() - 1);
    // Source (int)((value-1)/bin_size) is zero for this negative fraction,
    // not -1. All samples occupy bin0; its center is half the f32 spacing.
    let actual = estimator
        .estimate(&[0.; 100], &[f64::from(value); 100])
        .unwrap();
    let source_spacing = f64::from(value / 100.0_f32);
    assert_eq!(actual.max_intensity, 0.5 * source_spacing);
    assert_ne!(source_spacing, f64::from(value) / 100.0);
    assert_eq!(actual.noise, [1.; 100]);
    assert_eq!(actual.signal_to_noise, [f64::from(value); 100]);
}
