// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Source boundary semantics plus independent metadata and failure regressions.
//! Literal SRM goldens are independently checked in chromatogram_processing_reference.

use openms::Error;
use openms::kernel::{ChromatogramPeak, DataArray, MSChromatogram};
use openms::processing::chromatogram::{
    ChromatogramPickingMethod, ChromatogramSmoothing, PeakPickerChromatogram,
};
use openms::processing::smoothing::{GaussFilter, GaussianWidth, SavitzkyGolayFilter};

fn signal(intensities: &[f32]) -> MSChromatogram {
    MSChromatogram {
        peaks: intensities
            .iter()
            .enumerate()
            .map(|(i, &y)| ChromatogramPeak::new(i as f64, y))
            .collect(),
        ..Default::default()
    }
}
fn options() -> PeakPickerChromatogram {
    PeakPickerChromatogram {
        smoothing: ChromatogramSmoothing::SavitzkyGolay {
            frame_length: 1,
            polynomial_order: 0,
        },
        seed_signal_to_noise: 0.0,
        signal_to_noise: 0.0,
        ..Default::default()
    }
}
fn array<'a>(chromatogram: &'a MSChromatogram, name: &str) -> &'a [f32] {
    &chromatogram
        .float_data_arrays
        .iter()
        .find(|array| array.name == name)
        .unwrap()
        .data
}

#[test]
fn defaults_keep_seed_detection_and_boundary_noise_independent() {
    let config = PeakPickerChromatogram::default();
    assert_eq!(config.method, ChromatogramPickingMethod::Corrected);
    assert_eq!(
        config.smoothing,
        ChromatogramSmoothing::Gaussian { width: 50.0 }
    );
    assert_eq!(config.peak_width, None);
    assert_eq!(
        (config.signal_to_noise, config.seed_signal_to_noise),
        (1.0, 1.0)
    );
    assert_eq!(config.noise_estimator.window_length, 1000.0);
    assert_eq!(config.seed_noise_estimator.window_length, 200.0);
    assert_eq!(
        (
            config.noise_estimator.bin_count,
            config.seed_noise_estimator.bin_count
        ),
        (30, 30)
    );
    assert!(!config.report_sn);
    assert!(!config.remove_overlapping_peaks);
}

#[test]
fn regions_preserve_exact_sampling_and_raw_sums_are_not_time_integrals() {
    let mut input = signal(&[0., 1., 4., 10., 4., 1., 0.]);
    for (peak, rt) in input.peaks.iter_mut().zip([0., 1., 3., 6., 10., 15., 21.]) {
        peak.rt = rt;
    }
    for method in [
        ChromatogramPickingMethod::Legacy,
        ChromatogramPickingMethod::Corrected,
    ] {
        let result = PeakPickerChromatogram {
            method,
            ..options()
        }
        .pick_chromatogram(&input)
        .unwrap();
        assert_eq!(result.regions.len(), 1);
        let region = result.regions[0];
        assert_eq!(
            (region.apex_index, region.left_index, region.right_index),
            (3, 0, 6)
        );
        let picked = &result.picked.chromatogram;
        assert_eq!(array(picked, "IntegratedIntensity"), [20.0]);
        assert_eq!(array(picked, "leftWidth"), [0.0]);
        assert_eq!(array(picked, "rightWidth"), [21.0]);
        assert_eq!(array(picked, "SN"), [-1.0]);
        assert!(array(picked, "FWHM")[0] > 0.0);
        assert_eq!(result.smoothed, input);
    }
}

#[test]
fn forced_width_extends_nonmonotonic_flanks_but_never_overrides_noise_gate() {
    let input = signal(&[0., 1., 4., 10., 4., 5., 3., 2., 0.]);
    let ordinary = options().pick_chromatogram(&input).unwrap();
    let first = ordinary.regions.iter().find(|r| r.apex_index == 3).unwrap();
    assert_eq!(first.right_index, 4); // rise to 5 stops the descent
    let forced = PeakPickerChromatogram {
        peak_width: Some(20.0),
        ..options()
    }
    .pick_chromatogram(&input)
    .unwrap();
    let first = forced.regions.iter().find(|r| r.apex_index == 3).unwrap();
    assert_eq!((first.left_index, first.right_index), (0, 8));
    let gated = PeakPickerChromatogram {
        peak_width: Some(20.0),
        signal_to_noise: 1.0,
        ..options()
    }
    .pick_chromatogram(&input)
    .unwrap();
    let first = gated.regions.iter().find(|r| r.apex_index == 3).unwrap();
    // The nine-point sparse window has noise 1e20. Immediate neighbors remain
    // unconditional; second-away candidates fail S/N despite the forced width.
    assert_eq!((first.left_index, first.right_index), (2, 4));
    assert!(array(&gated.picked.chromatogram, "SN")[0] > 0.0);
}

#[test]
fn overlap_adjustment_uses_valleys_and_inclusive_shared_endpoints() {
    let input = signal(&[0., 1., 8., 5., 2., 4., 10., 6., 1., 0., 0.]);
    let forced = PeakPickerChromatogram {
        peak_width: Some(30.0),
        ..options()
    };
    let original = forced.pick_chromatogram(&input).unwrap();
    assert_eq!(original.regions.len(), 2);
    assert!(
        original
            .regions
            .iter()
            .all(|r| r.left_index == 0 && r.right_index == 10)
    );
    let adjusted = PeakPickerChromatogram {
        remove_overlapping_peaks: true,
        ..forced
    }
    .pick_chromatogram(&input)
    .unwrap();
    assert_eq!(
        adjusted.picked.chromatogram.peaks,
        original.picked.chromatogram.peaks
    );
    assert_eq!(
        (
            adjusted.regions[0].left_index,
            adjusted.regions[0].right_index
        ),
        (0, 4)
    );
    assert_eq!(
        (
            adjusted.regions[1].left_index,
            adjusted.regions[1].right_index
        ),
        (4, 10)
    );
    assert_eq!(
        array(&adjusted.picked.chromatogram, "IntegratedIntensity"),
        [16.0, 23.0]
    );
    // Both inclusive intervals contain the valley intensity 2.
    assert_eq!(input.peaks.iter().map(|p| p.intensity).sum::<f32>(), 37.0);
}

#[test]
fn both_smoothing_choices_reuse_existing_filters_and_preserve_the_sample_grid() {
    let input = signal(&[0., 1., 4., 10., 4., 1., 0., 1., 5., 1., 0.]);
    for smoothing in [
        ChromatogramSmoothing::Gaussian { width: 8.0 },
        ChromatogramSmoothing::SavitzkyGolay {
            frame_length: 4,
            polynomial_order: 2,
        },
    ] {
        let mut expected = input.clone();
        match smoothing {
            ChromatogramSmoothing::Gaussian { width } => {
                GaussFilter::new(GaussianWidth::Absolute(width))
                    .unwrap()
                    .filter_chromatogram(&mut expected)
                    .unwrap()
            }
            ChromatogramSmoothing::SavitzkyGolay {
                frame_length,
                polynomial_order,
            } => SavitzkyGolayFilter::new(frame_length, polynomial_order)
                .unwrap()
                .filter_chromatogram(&mut expected)
                .unwrap(),
        }
        let result = PeakPickerChromatogram {
            smoothing,
            ..options()
        }
        .pick_chromatogram(&input)
        .unwrap();
        assert_eq!(result.smoothed, expected);
        assert_eq!(
            result
                .smoothed
                .peaks
                .iter()
                .map(|p| p.rt)
                .collect::<Vec<_>>(),
            input.peaks.iter().map(|p| p.rt).collect::<Vec<_>>()
        );
    }
}

#[test]
fn metadata_survives_and_omitted_profile_arrays_are_explicit() {
    let mut input = signal(&[0., 1., 4., 10., 4., 1., 0.]);
    input.native_id = "transition-1".into();
    input.name = "peak trace".into();
    input.metadata.insert("sample".into(), "A".into());
    input.precursor.mz = 456.7;
    input.precursor.charge = 2;
    input
        .float_data_arrays
        .push(DataArray::new("profile float", vec![1.; 7]));
    input
        .integer_data_arrays
        .push(DataArray::new("profile integer", (0..7).collect()));
    input
        .string_data_arrays
        .push(DataArray::new("profile label", vec!["x".into(); 7]));
    let original = input.clone();
    let result = options().pick_chromatogram(&input).unwrap();
    assert_eq!(input, original);
    assert_eq!(result.smoothed, input);
    assert_eq!(
        result.picked.omitted_arrays,
        ["profile float", "profile integer", "profile label"]
    );
    let picked = &result.picked.chromatogram;
    assert_eq!(picked.native_id, input.native_id);
    assert_eq!(picked.name, input.name);
    assert_eq!(picked.metadata, input.metadata);
    assert_eq!(picked.precursor, input.precursor);
    assert_eq!(
        picked
            .float_data_arrays
            .iter()
            .map(|a| a.name.as_str())
            .collect::<Vec<_>>(),
        [
            "FWHM",
            "IntegratedIntensity",
            "leftWidth",
            "rightWidth",
            "SN"
        ]
    );
    assert!(picked.integer_data_arrays.is_empty());
    assert!(picked.string_data_arrays.is_empty());
    assert_eq!(result.picked.boundaries.len(), picked.len());
    let mut filtered = input.clone();
    options().filter_chromatogram(&mut filtered).unwrap();
    assert_eq!(&filtered, picked);
    picked.validate().unwrap();
}

#[test]
fn empty_short_constant_and_zero_signals_have_well_formed_empty_results() {
    for input in [
        MSChromatogram::default(),
        signal(&[1., 2., 1.]),
        signal(&[1.; 9]),
        signal(&[0.; 9]),
    ] {
        let result = options().pick_chromatogram(&input).unwrap();
        assert!(result.picked.chromatogram.is_empty());
        assert!(result.regions.is_empty());
        assert_eq!(result.picked.chromatogram.float_data_arrays.len(), 5);
        assert!(
            result
                .picked
                .chromatogram
                .float_data_arrays
                .iter()
                .all(|a| a.data.is_empty())
        );
        assert_eq!(result.smoothed, input);
        result.picked.chromatogram.validate().unwrap();
    }
}

#[test]
fn invalid_inputs_options_and_exhausted_work_leave_mutable_input_unchanged() {
    let input = signal(&[0., 1., 4., 10., 4., 1., 0.]);
    for config in [
        PeakPickerChromatogram {
            signal_to_noise: -1.,
            ..options()
        },
        PeakPickerChromatogram {
            seed_signal_to_noise: f64::NAN,
            ..options()
        },
        PeakPickerChromatogram {
            peak_width: Some(0.0),
            ..options()
        },
        PeakPickerChromatogram {
            smoothing: ChromatogramSmoothing::Gaussian {
                width: f64::INFINITY,
            },
            ..options()
        },
        PeakPickerChromatogram {
            smoothing: ChromatogramSmoothing::SavitzkyGolay {
                frame_length: 3,
                polynomial_order: 3,
            },
            ..options()
        },
        PeakPickerChromatogram {
            max_points: 6,
            ..options()
        },
        PeakPickerChromatogram {
            max_work: 1,
            ..options()
        },
    ] {
        let mut unchanged = input.clone();
        assert!(config.filter_chromatogram(&mut unchanged).is_err());
        assert_eq!(unchanged, input);
    }
    for kind in 0..4 {
        let mut invalid = input.clone();
        match kind {
            0 => invalid.peaks.reverse(),
            1 => invalid.peaks[3].rt = invalid.peaks[2].rt,
            2 => invalid.peaks[3].intensity = -1.0,
            _ => invalid
                .float_data_arrays
                .push(DataArray::new("bad", vec![1.0])),
        }
        let before = invalid.clone();
        let error = options().filter_chromatogram(&mut invalid).unwrap_err();
        if kind == 0 {
            assert!(matches!(error, Error::UnsortedData));
        }
        assert_eq!(invalid, before);
    }
}

#[test]
fn source_f32_output_metadata_overflow_is_checked_atomically() {
    let mut input = signal(&[0., 1., 4., 10., 4., 1., 0.]);
    for (index, peak) in input.peaks.iter_mut().enumerate() {
        peak.rt = 1e39 + index as f64 * 1e33;
    }
    let original = input.clone();
    assert!(options().filter_chromatogram(&mut input).is_err());
    assert_eq!(input, original);
    let mut input = signal(&[0., 1., f32::MAX / 2., f32::MAX, f32::MAX / 2., 1., 0.]);
    let original = input.clone();
    assert!(options().filter_chromatogram(&mut input).is_err());
    assert_eq!(input, original);
}
