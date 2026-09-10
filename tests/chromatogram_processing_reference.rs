// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Independent pinned-source chromatogram and integration oracles.
//! Source literals, expected assertions and SHA-256 provenance are in
//! tests/data/chromatogram_processing_*. No C++ code was executed.

use openms::analysis::peak_integrator::{BaselineType, IntegrationMethod, PeakIntegrator};
use openms::kernel::{ChromatogramPeak, MSChromatogram, MSSpectrum, Peak1D};
use openms::processing::chromatogram::{
    ChromatogramPickingMethod, ChromatogramSmoothing, PeakPickerChromatogram,
};

fn rows(text: &str) -> impl Iterator<Item = Vec<&str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .skip(1)
        .map(|line| line.split('\t').collect())
}
fn trace(name: &str) -> MSChromatogram {
    let peaks = rows(include_str!("data/chromatogram_processing_traces.tsv"))
        .filter(|row| row[0] == name)
        .enumerate()
        .map(|(index, row)| {
            assert_eq!(row[1].parse::<usize>().unwrap(), index);
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
fn chromatogram(values: &[(f64, f32)]) -> MSChromatogram {
    MSChromatogram {
        peaks: values
            .iter()
            .map(|&(x, y)| ChromatogramPeak::new(x, y))
            .collect(),
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
fn method(name: &str) -> IntegrationMethod {
    match name {
        "intensity_sum" => IntegrationMethod::IntensitySum,
        "trapezoid" => IntegrationMethod::Trapezoid,
        "simpson" => IntegrationMethod::Simpson,
        _ => panic!("unknown source method {name}"),
    }
}

#[test]
fn pinned_integration_areas_apices_and_sample_outlines() {
    for row in rows(include_str!("data/chromatogram_processing_integrals.tsv")) {
        let input = trace(row[0]);
        let config = PeakIntegrator {
            integration_method: method(row[1]),
            ..Default::default()
        };
        let left: f64 = row[2].parse().unwrap();
        let right: f64 = row[3].parse().unwrap();
        let chrom = config.integrate_chromatogram(&input, left, right).unwrap();
        let spec = config
            .integrate_spectrum(&spectrum(&input), left, right)
            .unwrap();
        // Source's full-interval trapezoid assertion has only six significant digits.
        let tolerance = if row[1] == "trapezoid" { 0.01 } else { 1e-6 };
        close(chrom.area, row[4].parse().unwrap(), tolerance);
        close(chrom.height, row[5].parse().unwrap(), 0.0);
        close(chrom.apex_pos, row[6].parse().unwrap(), 0.0);
        assert_eq!(spec, chrom);
        let sampled: Vec<_> = input
            .peaks
            .iter()
            .filter(|p| left <= p.rt && p.rt <= right)
            .map(|p| [p.rt, f64::from(p.intensity)])
            .collect();
        assert_eq!(chrom.hull_points, sampled);
    }
}

#[test]
fn pinned_background_values_use_selected_sample_endpoints() {
    let input = trace("glutamate");
    for row in rows(include_str!("data/chromatogram_processing_backgrounds.tsv")) {
        let baseline_type = match row[1] {
            "base_to_base" => BaselineType::BaseToBase,
            "vertical_division_min" => BaselineType::VerticalDivisionMin,
            "vertical_division_max" => BaselineType::VerticalDivisionMax,
            _ => unreachable!(),
        };
        let config = PeakIntegrator {
            integration_method: method(row[0]),
            baseline_type,
            ..Default::default()
        };
        let actual = config
            .estimate_background_chromatogram(&input, 2.472833334, 3.022891666, 2.7045)
            .unwrap();
        close(actual.area, row[2].parse().unwrap(), 1e-6);
        close(actual.height, row[3].parse().unwrap(), 1e-8);
        assert_eq!(
            actual,
            config
                .estimate_background_spectrum(&spectrum(&input), 2.472833334, 3.022891666, 2.7045)
                .unwrap()
        );
    }
}

#[test]
fn pinned_shape_metrics_and_cropped_sample_crossings() {
    let input = trace("glutamate");
    let config = PeakIntegrator::default();
    let actual = config
        .calculate_shape_metrics_chromatogram(&input, 2.472833334, 3.022891666, 966489.0, 2.7045)
        .unwrap();
    assert_eq!(
        actual,
        config
            .calculate_shape_metrics_spectrum(
                &spectrum(&input),
                2.472833334,
                3.022891666,
                966489.0,
                2.7045
            )
            .unwrap()
    );
    for row in rows(include_str!("data/chromatogram_processing_shape.tsv")) {
        let value = match row[0] {
            "width_at_5" => actual.width_at_5,
            "width_at_10" => actual.width_at_10,
            "width_at_50" => actual.width_at_50,
            "start_position_at_5" => actual.start_position_at_5,
            "start_position_at_10" => actual.start_position_at_10,
            "start_position_at_50" => actual.start_position_at_50,
            "end_position_at_5" => actual.end_position_at_5,
            "end_position_at_10" => actual.end_position_at_10,
            "end_position_at_50" => actual.end_position_at_50,
            "total_width" => actual.total_width,
            "tailing_factor" => actual.tailing_factor,
            "asymmetry_factor" => actual.asymmetry_factor,
            "slope_of_baseline" => actual.slope_of_baseline,
            "baseline_delta_2_height" => actual.baseline_delta_2_height,
            "points_across_baseline" => actual.points_across_baseline as f64,
            "points_across_half_height" => actual.points_across_half_height as f64,
            _ => panic!("unknown shape field {}", row[0]),
        };
        close(value, row[1].parse().unwrap(), 1e-10);
    }
    // Exact crop indices from the source tests. Missing crossings remain at
    // the selected sampled edges rather than extrapolating the peak flanks.
    for (left_index, right_index, cutoff) in [(41, 54, 5), (42, 53, 10), (44, 49, 50), (46, 48, 5)]
    {
        let left = input.peaks[left_index].rt;
        let right = input.peaks[right_index].rt;
        let shape = config
            .calculate_shape_metrics_chromatogram(&input, left, right, 966489.0, 2.7045)
            .unwrap();
        assert_eq!(
            (shape.start_position_at_5, shape.end_position_at_5),
            (left, right)
        );
        if cutoff >= 10 {
            assert_eq!(
                (shape.start_position_at_10, shape.end_position_at_10),
                (left, right)
            );
        }
        if cutoff >= 50 {
            assert_eq!(
                (shape.start_position_at_50, shape.end_position_at_50),
                (left, right)
            );
        }
    }
}

#[test]
fn nonuniform_even_simpson_averages_source_neighbor_intervals() {
    // A quadratic is integrated exactly by every three-point unequal-grid
    // Simpson segment. Analytic integrals, not a second numerical routine,
    // therefore expose which neighboring intervals are included in the mean.
    let input = chromatogram(&[
        (0., 0.),
        (1., 1.),
        (2., 4.),
        (4., 16.),
        (7., 49.),
        (11., 121.),
    ]);
    let config = PeakIntegrator {
        integration_method: IntegrationMethod::Simpson,
        ..Default::default()
    };
    for (left, right, expected) in [(1., 7., 2071. / 12.), (0., 4., 46.), (2., 11., 2932. / 9.)] {
        close(
            config
                .integrate_chromatogram(&input, left, right)
                .unwrap()
                .area,
            expected,
            1e-12,
        );
    }
}

#[test]
fn trapezoid_preserves_f32_addition_before_f64_scaling() {
    let input = chromatogram(&[(0., 16_777_216.), (2., 1.)]);
    let config = PeakIntegrator {
        integration_method: IntegrationMethod::Trapezoid,
        ..Default::default()
    };
    // f32(2^24 + 1) rounds to 2^24; a premature f64 promotion changes this by 1.
    assert_eq!(
        config.integrate_chromatogram(&input, 0., 2.).unwrap().area,
        16_777_216.
    );
}

#[test]
fn source_even_simpson_minus_one_sentinel_differs_from_odd_integration() {
    let input = chromatogram(&[(0., 1.5), (1., 0.), (4., 0.), (5., 0.)]);
    let config = PeakIntegrator {
        integration_method: IntegrationMethod::Simpson,
        ..Default::default()
    };
    // h=1,k=3 gives coefficient -2/3 for the first intensity: -1 exactly.
    assert_eq!(
        config.integrate_chromatogram(&input, 0., 4.).unwrap().area,
        -1.
    );
    // In even mode source excludes this valid -1 term using its missing-slot
    // sentinel. The remaining odd segment contributes zero; retain finite parity.
    assert_eq!(
        config.integrate_chromatogram(&input, 0., 5.).unwrap().area,
        0.
    );
    let no_finite_average = chromatogram(&[(0., 1.5), (1., 0.), (4., 0.), (5., 1.5)]);
    // Both available odd subareas equal the sentinel. Source divides zero by
    // zero here; the native checked API reports an error instead.
    assert!(
        config
            .integrate_chromatogram(&no_finite_average, 0., 5.)
            .is_err()
    );
}

#[test]
fn sampled_linear_background_and_vertical_division_alias() {
    let input = chromatogram(&[(0., 2.), (1., 3.), (4., 6.)]);
    let alias: BaselineType = "vertical_division".parse().unwrap();
    assert_eq!(alias, BaselineType::VerticalDivisionMin);
    let base = PeakIntegrator::default();
    let background = base
        .estimate_background_chromatogram(&input, -0.5, 4.5, 1.)
        .unwrap();
    assert_eq!((background.area, background.height), (11., 3.));
    // Uneven spacing means sum of sampled baseline differs from n*mean(endpoints).
    let trapezoid = PeakIntegrator {
        integration_method: IntegrationMethod::Trapezoid,
        ..base
    };
    let continuous = trapezoid
        .estimate_background_chromatogram(&input, -0.5, 4.5, 1.)
        .unwrap();
    assert_eq!((continuous.area, continuous.height), (16., 3.));
}

#[test]
fn shape_baseline_delta_and_background_have_distinct_source_precision() {
    let input = chromatogram(&[(0., 1.), (1., 33_554_432.), (2., 33_554_432.)]);
    let config = PeakIntegrator::default();
    let shape = config
        .calculate_shape_metrics_chromatogram(&input, 0., 2., 33_554_432., 1.)
        .unwrap();
    // Shape's f32 subtraction rounds 2^25 - 1 upward before promotion.
    assert_eq!(shape.slope_of_baseline, 33_554_432.);
    assert_eq!(shape.baseline_delta_2_height, 1.);
    // Background promotes each endpoint first; its line retains that unit.
    let background = config
        .estimate_background_chromatogram(&input, 0., 2., 1.)
        .unwrap();
    assert_eq!(background.height, 16_777_216.5);
    assert_eq!(background.area, 50_331_649.5);
}

#[test]
fn both_source_srm_traces_match_legacy_and_corrected_peak_goldens() {
    for row in rows(include_str!("data/chromatogram_processing_picked.tsv")) {
        let input = trace(row[0]);
        let picker = PeakPickerChromatogram {
            method: if row[1] == "legacy" {
                ChromatogramPickingMethod::Legacy
            } else {
                ChromatogramPickingMethod::Corrected
            },
            peak_width: if row[2] == "-1" {
                None
            } else {
                Some(row[2].parse().unwrap())
            },
            ..Default::default()
        };
        let result = picker.pick_chromatogram(&input).unwrap();
        let picked = &result.picked.chromatogram;
        assert_eq!(picked.len(), 1);
        assert_eq!(picked.float_data_arrays.len(), 5);
        close(picked.peaks[0].rt, row[3].parse().unwrap(), 0.001);
        close(
            f64::from(picked.peaks[0].intensity),
            row[4].parse().unwrap(),
            0.02,
        );
        let array = |name: &str| {
            picked
                .float_data_arrays
                .iter()
                .find(|a| a.name == name)
                .unwrap()
                .data[0]
        };
        close(
            f64::from(array("IntegratedIntensity")),
            row[5].parse().unwrap(),
            0.6, // Upstream areas are approximate six-significant-digit goldens.
        );
        assert_eq!(array("leftWidth"), row[6].parse::<f32>().unwrap());
        assert_eq!(array("rightWidth"), row[7].parse::<f32>().unwrap());
        let region = &result.regions[0];
        assert_eq!(
            input.peaks[region.left_index].rt,
            row[6].parse::<f64>().unwrap()
        );
        assert_eq!(
            input.peaks[region.right_index].rt,
            row[7].parse::<f64>().unwrap()
        );
        // Independently sum source raw samples with the recorded boundaries.
        let raw_sum: f64 = input.peaks[region.left_index..=region.right_index]
            .iter()
            .map(|p| f64::from(p.intensity))
            .sum();
        assert_eq!(array("IntegratedIntensity"), raw_sum as f32);
        assert_eq!(result.smoothed.len(), input.len());
    }
}

#[test]
fn source_seed_noise_is_separate_from_boundary_threshold_and_reporting() {
    let input = chromatogram(&[
        (0., 0.),
        (1., 1.),
        (2., 4.),
        (3., 10.),
        (4., 4.),
        (5., 1.),
        (6., 0.),
    ]);
    let config = PeakPickerChromatogram {
        smoothing: ChromatogramSmoothing::SavitzkyGolay {
            frame_length: 3,
            polynomial_order: 2,
        },
        signal_to_noise: 0.,
        ..Default::default()
    };
    // Fewer than the noise estimator's required ten samples makes S/N tiny.
    // Source still applies its internal seed threshold1 when outer threshold0.
    assert!(
        config
            .pick_chromatogram(&input)
            .unwrap()
            .picked
            .chromatogram
            .is_empty()
    );
    let config = PeakPickerChromatogram {
        seed_signal_to_noise: 0.,
        ..config
    };
    let result = config.pick_chromatogram(&input).unwrap();
    assert_eq!(result.picked.chromatogram.len(), 1);
    assert_eq!(
        (result.regions[0].left_index, result.regions[0].right_index),
        (0, 6)
    );
    let array = |result: &openms::processing::chromatogram::ChromatogramPickingResult,
                 name: &str| {
        result
            .picked
            .chromatogram
            .float_data_arrays
            .iter()
            .find(|a| a.name == name)
            .unwrap()
            .data[0]
    };
    assert_eq!(array(&result, "IntegratedIntensity"), 20.);
    assert_eq!(array(&result, "SN"), -1.);
    let reported = PeakPickerChromatogram {
        report_sn: true,
        ..config.clone()
    }
    .pick_chromatogram(&input)
    .unwrap();
    assert_eq!(array(&reported, "SN"), (10.0_f64 / 1e20) as f32);
    let stopped = PeakPickerChromatogram {
        signal_to_noise: 1.,
        ..config
    }
    .pick_chromatogram(&input)
    .unwrap();
    assert_eq!(
        (
            stopped.regions[0].left_index,
            stopped.regions[0].right_index
        ),
        (2, 4)
    );
    assert_eq!(array(&stopped, "IntegratedIntensity"), 18.);
}

#[test]
fn overlap_resolution_retains_both_peaks_and_sums_the_shared_valley_twice() {
    let input = chromatogram(&[
        (0., 0.),
        (1., 1.),
        (2., 4.),
        (3., 10.),
        (4., 4.),
        (5., 1.),
        (6., 5.),
        (7., 12.),
        (8., 5.),
        (9., 1.),
        (10., 0.),
    ]);
    let picker = PeakPickerChromatogram {
        smoothing: ChromatogramSmoothing::SavitzkyGolay {
            frame_length: 3,
            polynomial_order: 2,
        },
        peak_width: Some(100.),
        signal_to_noise: 0.,
        seed_signal_to_noise: 0.,
        ..Default::default()
    };
    let original = picker.pick_chromatogram(&input).unwrap();
    assert_eq!(original.regions.len(), 2);
    for region in &original.regions {
        assert_eq!((region.left_index, region.right_index), (0, 10));
    }
    let corrected = PeakPickerChromatogram {
        remove_overlapping_peaks: true,
        ..picker
    }
    .pick_chromatogram(&input)
    .unwrap();
    assert_eq!(corrected.regions.len(), 2);
    assert_eq!(
        (
            corrected.regions[0].left_index,
            corrected.regions[0].right_index
        ),
        (0, 5)
    );
    assert_eq!(
        (
            corrected.regions[1].left_index,
            corrected.regions[1].right_index
        ),
        (5, 10)
    );
    let areas = &corrected
        .picked
        .chromatogram
        .float_data_arrays
        .iter()
        .find(|a| a.name == "IntegratedIntensity")
        .unwrap()
        .data;
    assert_eq!(areas, &[20., 24.]);
    // Source keeps both inclusive valley endpoints. Their total 44 includes the
    // shared intensity 1 twice; the original sampled signal sums to 43.
}

#[test]
fn gaussian_table_setup_is_charged_even_for_an_empty_chromatogram() {
    // Width 2000 and the fixed 0.01 table spacing require 100001 coefficients.
    // Empty sampling does not remove that setup work in the shared smoother.
    for limit in [1, 100_000] {
        let picker = PeakPickerChromatogram {
            smoothing: ChromatogramSmoothing::Gaussian { width: 2000. },
            max_work: limit,
            ..Default::default()
        };
        let mut input = MSChromatogram::default();
        let before = input.clone();
        assert!(picker.filter_chromatogram(&mut input).is_err());
        assert_eq!(input, before);
    }
    let exact = PeakPickerChromatogram {
        smoothing: ChromatogramSmoothing::Gaussian { width: 2000. },
        max_work: 100_001,
        ..Default::default()
    };
    assert!(
        exact
            .pick_chromatogram(&MSChromatogram::default())
            .unwrap()
            .picked
            .chromatogram
            .is_empty()
    );
}
