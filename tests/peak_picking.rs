// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::kernel::{
    ChromatogramPeak, DataArray, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, SpectrumType,
};
use openms::processing::SpectrumFilter;
use openms::processing::peak_picking::*;
fn signal(x: &[f64], y: &[f32]) -> MSSpectrum {
    MSSpectrum::from_peaks(x.iter().zip(y).map(|(&x, &y)| Peak1D::new(x, y)).collect())
}
fn fixture(text: &str) -> MSSpectrum {
    MSSpectrum::from_peaks(
        text.lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let a: Vec<_> = l.split_whitespace().collect();
                Peak1D::new(a[0].parse().unwrap(), a[1].parse().unwrap())
            })
            .collect(),
    )
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a} != {b}; tolerance {tol}");
}
#[test]
fn upstream_natural_spline_values_and_derivatives() {
    let x = [
        486.784, 486.787, 486.790, 486.793, 486.795, 486.797, 486.800, 486.802, 486.805, 486.808,
        486.811,
    ];
    let y = [
        0., 154683.17, 620386.5, 1701390.12, 2848879.25, 3564045.5, 2744585.7, 1605583., 1518984.,
        1591352.21, 1691345.1,
    ];
    let s = CubicSpline2d::new(&x, &y).unwrap();
    close(s.eval(486.785).unwrap(), 35173.1841778984, 1e-5);
    close(s.eval(486.794).unwrap(), 2271426.93316241, 1e-5);
    close(s.derivative(486.785, 1).unwrap(), 39270152.2996247, 0.001);
    close(s.derivative(486.794, 2).unwrap(), 7415503644.8958, 0.01);
    close(s.derivative(x[0], 2).unwrap(), 0.0, 1e-9);
    close(s.derivative(x[x.len() - 1], 2).unwrap(), 0.0, 0.001);
    for (&x, &y) in x.iter().zip(&y) {
        close(s.eval(x).unwrap(), y, 1e-8);
    }
}
#[test]
fn upstream_dummy_peak_and_weighted_mobility() {
    let mut input = signal(
        &[100., 100.01, 100.02, 100.03, 100.04],
        &[200., 250., 450., 250., 200.],
    );
    input.name = "acquisition".into();
    input.metadata.insert("sample".into(), "A".into());
    input.float_data_arrays.push(DataArray::new(
        "Ion Mobility",
        vec![100., 150., 150., 150., 100.],
    ));
    input
        .integer_data_arrays
        .push(DataArray::new("labels", vec![1, 2, 3, 4, 5]));
    let result = PeakPickerHiRes::default().pick_spectrum(&input).unwrap();
    assert_eq!(result.spectrum.len(), 1);
    close(result.spectrum.peaks[0].mz, 100.02, 1e-6);
    close(f64::from(result.spectrum.peaks[0].intensity), 450., 1e-3);
    close(
        f64::from(result.spectrum.float_data_arrays[0].data[0]),
        135.1852,
        1e-4,
    );
    assert_eq!(result.omitted_arrays, ["labels"]);
    assert_eq!(result.spectrum.metadata, input.metadata);
    assert_eq!(result.spectrum.name, input.name);
    assert_eq!(result.spectrum.spectrum_type, SpectrumType::Centroid);
    assert_eq!(
        result.boundaries,
        [PeakBoundary {
            min: 100.,
            max: 100.04
        }]
    );
}
#[test]
fn upstream_real_data_centroids_match_source_fixtures() {
    let cases = [
        (
            include_str!("data/peak_picking_orbitrap.tsv"),
            include_str!("data/peak_picking_orbitrap_sn1_out.tsv"),
            1.0,
        ),
        (
            include_str!("data/peak_picking_orbitrap.tsv"),
            include_str!("data/peak_picking_orbitrap_sn4_out.tsv"),
            4.0,
        ),
        (
            include_str!("data/peak_picking_ftms.tsv"),
            include_str!("data/peak_picking_ftms_sn1_out.tsv"),
            1.0,
        ),
        (
            include_str!("data/peak_picking_ftms.tsv"),
            include_str!("data/peak_picking_ftms_sn4_out.tsv"),
            4.0,
        ),
    ];
    for (source, expected, sn) in cases {
        let picker = PeakPickerHiRes {
            signal_to_noise: sn,
            ..Default::default()
        };
        let actual = picker.pick_spectrum(&fixture(source)).unwrap();
        let expected = fixture(expected);
        assert_eq!(actual.spectrum.len(), expected.len());
        for (a, b) in actual.spectrum.peaks.iter().zip(&expected.peaks) {
            close(a.mz, b.mz, 1e-6);
            close(
                f64::from(a.intensity),
                f64::from(b.intensity),
                f64::from(b.intensity) * 2e-5 + 1e-4,
            );
        }
    }
}
#[test]
fn upstream_noise_values_match_source_fixture() {
    let input =
        openms::format::dta::read(include_bytes!("data/peak_picking_noise_input.dta").as_slice())
            .unwrap();
    let historical = openms::format::dta::read(
        include_bytes!("data/peak_picking_noise_historical.dta").as_slice(),
    )
    .unwrap();
    let config = SignalToNoiseEstimatorMedian {
        window_length: 40.,
        noise_for_empty_window: 2.,
        ..Default::default()
    };
    let x: Vec<_> = input.peaks.iter().map(|p| p.mz).collect();
    let y: Vec<_> = input.peaks.iter().map(|p| f64::from(p.intensity)).collect();
    let actual = config.estimate(&x, &y).unwrap();
    assert_eq!(actual.signal_to_noise.len(), historical.len());
    for (a, b) in actual.signal_to_noise.iter().zip(&historical.peaks) {
        close(
            *a,
            f64::from(b.intensity),
            f64::from(b.intensity).abs() * 1e-5 + 1e-7,
        );
    }
}
#[test]
fn source_boundaries_and_missing_flank_mobility() {
    let input = fixture(include_str!("data/peak_picking_orbitrap.tsv"));
    let out = PeakPickerHiRes {
        signal_to_noise: 1.0,
        ..Default::default()
    }
    .pick_spectrum(&input)
    .unwrap();
    close(out.boundaries[25].min, 367.206604003906, 1e-10);
    close(out.boundaries[25].max, 367.214569091797, 1e-10);
    close(out.boundaries[26].min, 369.042205810547, 1e-10);
    close(out.boundaries[26].max, 369.051574707031, 1e-10);
    let mut input = signal(
        &[100., 100.03, 100.06, 100.07, 100.08],
        &[200., 250., 450., 250., 200.],
    );
    input
        .float_data_arrays
        .push(DataArray::new("Ion Mobility", vec![1., 1.1, 1.2, 1.3, 1.4]));
    assert!(
        PeakPickerHiRes::default()
            .pick_spectrum(&input)
            .unwrap()
            .spectrum
            .is_empty()
    );
    let picker = PeakPickerHiRes {
        allow_missing_flank: true,
        ..Default::default()
    };
    let out = picker.pick_spectrum(&input).unwrap();
    assert_eq!(out.spectrum.len(), 1);
    close(
        f64::from(out.spectrum.float_data_arrays[0].data[0]),
        1.27222,
        1e-5,
    );
    let input = signal(
        &[100., 100.01, 100.02, 100.05, 100.08],
        &[200., 250., 450., 250., 200.],
    );
    assert!(
        PeakPickerHiRes::default()
            .pick_spectrum(&input)
            .unwrap()
            .spectrum
            .is_empty()
    );
    assert_eq!(picker.pick_spectrum(&input).unwrap().spectrum.len(), 1);
}
#[test]
fn histogram_interpolates_median_and_reports_sparse_and_clipped_windows() {
    let config = SignalToNoiseEstimatorMedian {
        histogram_range: NoiseHistogramRange::Manual { max_intensity: 30. },
        bin_count: 3,
        min_required_elements: 1,
        ..Default::default()
    };
    let out = config
        .estimate(&[0., 1., 2., 3., 4.], &[1., 2., 3., 4., 5.])
        .unwrap();
    assert_eq!(out.noise, vec![6.; 5]);
    close(out.signal_to_noise[4], 5. / 6., 1e-15);
    assert_eq!(out.sparse_window_percent, 0.);
    let out = config.estimate(&[0., 1., 2.], &[30., 40., 50.]).unwrap();
    close(out.noise[0], 20. + 20. / 3., 1e-14);
    assert_eq!(out.histogram_rightmost_percent, 100.);
    let sparse = SignalToNoiseEstimatorMedian {
        window_length: 1.,
        noise_for_empty_window: 2.,
        ..Default::default()
    };
    let out = sparse.estimate(&[0., 1000.], &[2., 4.]).unwrap();
    assert_eq!(out.signal_to_noise, [1., 2.]);
    assert_eq!(out.sparse_window_percent, 100.);
    let empty = config.estimate(&[], &[]).unwrap();
    assert!(empty.signal_to_noise.is_empty());
    assert_eq!(empty.sparse_window_percent, 0.);
    let zero = SignalToNoiseEstimatorMedian {
        min_required_elements: 1,
        ..Default::default()
    }
    .estimate(&[0., 1.], &[0., 0.])
    .unwrap();
    assert_eq!(zero.signal_to_noise, [0., 0.]);
    assert_eq!(zero.noise, [1., 1.]);
}
#[test]
fn fwhm_uses_spline_half_height_and_clips_to_available_support() {
    let input = signal(
        &[98., 99., 100., 101., 102.],
        &[200., 250., 450., 250., 200.],
    );
    let mut picker = PeakPickerHiRes {
        report_fwhm: Some(FwhmUnit::Absolute),
        ..Default::default()
    };
    let absolute = picker.pick_spectrum(&input).unwrap();
    let width = absolute.spectrum.float_data_arrays[0].data[0];
    assert_eq!(absolute.spectrum.float_data_arrays[0].name, "FWHM");
    assert!(width > 2. && width < 4.);
    picker.report_fwhm = Some(FwhmUnit::Ppm);
    let ppm = picker.pick_spectrum(&input).unwrap();
    assert_eq!(ppm.spectrum.float_data_arrays[0].name, "FWHM_ppm");
    close(
        f64::from(ppm.spectrum.float_data_arrays[0].data[0]),
        f64::from(width) * 1e4,
        0.01,
    );
    let clipped = signal(
        &[98., 99., 100., 101., 102.],
        &[300., 350., 450., 350., 300.],
    );
    picker.report_fwhm = Some(FwhmUnit::Absolute);
    assert_eq!(
        picker
            .pick_spectrum(&clipped)
            .unwrap()
            .spectrum
            .float_data_arrays[0]
            .data,
        [4.]
    );
}
#[test]
fn chromatograms_disable_spacing_and_keep_acquisition_metadata() {
    let input = MSChromatogram {
        peaks: [(0., 200.), (3., 250.), (6., 450.), (7., 250.), (8., 200.)]
            .into_iter()
            .map(|(x, y)| ChromatogramPeak::new(x, y))
            .collect(),
        name: "TIC".into(),
        native_id: "TIC=1".into(),
        ..Default::default()
    };
    let picker = PeakPickerHiRes {
        report_fwhm: Some(FwhmUnit::Absolute),
        ..Default::default()
    };
    let out = picker.pick_chromatogram(&input).unwrap();
    assert_eq!(out.chromatogram.len(), 1);
    assert_eq!(out.chromatogram.name, "TIC");
    assert_eq!(out.chromatogram.native_id, "TIC=1");
    assert_eq!(out.boundaries, [PeakBoundary { min: 0., max: 8. }]);
    assert!(
        picker
            .pick_chromatogram_with_spacing(&input, true)
            .unwrap()
            .chromatogram
            .is_empty()
    );
}
#[test]
fn source_peak_core_rejections_and_gap_missing_boundaries() {
    let picker = PeakPickerHiRes::default();
    for input in [
        signal(&[0., 1., 2., 3.], &[1., 2., 1., 0.]),
        signal(&[0., 1., 2., 3., 4.], &[0., 0., 4., 2., 0.]),
        signal(&[0., 1., 2., 3., 4.], &[1., 2., 2., 2., 1.]),
        signal(&[0., 1., 2., 3., 4.], &[5., 2., 4., 2., 5.]),
    ] {
        assert!(picker.pick_spectrum(&input).unwrap().spectrum.is_empty());
    }
    // Gap at the first left extension rejects that sample. Right support remains intact.
    let input = signal(&[0., 5., 6., 7., 8.], &[1., 2., 4., 2., 1.]);
    let out = picker.pick_spectrum(&input).unwrap();
    assert_eq!(out.spectrum.len(), 1);
    assert_eq!(out.boundaries[0].min, 5.);
    // A rejected missing sample is nevertheless reported as the source boundary.
    let input = signal(&[0., 2., 3., 4., 5.], &[1., 2., 4., 2., 1.]);
    let no_missing = PeakPickerHiRes {
        missing: 0,
        ..Default::default()
    };
    let out = no_missing.pick_spectrum(&input).unwrap();
    assert_eq!(out.boundaries[0].min, 0.);
    let disabled = PeakPickerHiRes {
        spacing_difference: 0.,
        spacing_difference_gap: 0.,
        ..Default::default()
    };
    assert_eq!(
        disabled
            .pick_spectrum(&signal(
                &[0., 3., 6., 7., 8.],
                &[200., 250., 450., 250., 200.]
            ))
            .unwrap()
            .spectrum
            .len(),
        1
    );
}
#[test]
fn experiment_auto_manual_selection_and_atomic_errors() {
    let profile = signal(
        &[100., 100.01, 100.02, 100.03, 100.04, 100.05, 100.06],
        &[1., 10., 30., 50., 30., 10., 1.],
    );
    assert_eq!(
        estimate_spectrum_type(&profile).unwrap(),
        SpectrumType::Profile
    );
    assert_eq!(
        estimate_spectrum_type(&signal(&[1., 2., 3., 4., 5.], &[1., 1., 1., 1., 1.])).unwrap(),
        SpectrumType::Centroid
    );
    assert_eq!(
        estimate_spectrum_type(&signal(&[1.], &[1.])).unwrap(),
        SpectrumType::Unknown
    );
    let mut centroid = profile.clone();
    centroid.spectrum_type = SpectrumType::Centroid;
    centroid.ms_level = 2;
    let mut experiment = MSExperiment {
        spectra: vec![profile, centroid.clone()],
        ..Default::default()
    };
    experiment
        .settings
        .metadata
        .insert("study".into(), "test".into());
    let out = PeakPickerHiRes::default()
        .pick_experiment(&experiment)
        .unwrap();
    assert_eq!(out.experiment.spectra[0].len(), 1);
    assert_eq!(out.experiment.spectra[1], centroid);
    assert!(out.spectrum_boundaries[0].is_some());
    assert!(out.spectrum_boundaries[1].is_none());
    assert_eq!(
        out.experiment.settings.metadata,
        experiment.settings.metadata
    );
    let manual = PeakPickerHiRes {
        ms_levels: vec![2],
        ..Default::default()
    };
    let before = experiment.clone();
    assert!(manual.filter_experiment(&mut experiment).is_err());
    assert_eq!(experiment, before);
    let unchecked = PeakPickerHiRes {
        check_spectrum_type: false,
        ..manual
    };
    let out = unchecked.pick_experiment(&experiment).unwrap();
    assert_eq!(out.experiment.spectra[0], experiment.spectra[0]);
    assert_eq!(out.experiment.spectra[1].len(), 1);
}
#[test]
fn checked_resources_invalid_data_and_spline_boundaries() {
    let valid = signal(&[0., 1., 2., 3., 4.], &[1., 2., 4., 2., 1.]);
    for picker in [
        PeakPickerHiRes {
            max_points: 4,
            ..Default::default()
        },
        PeakPickerHiRes {
            max_work: 1,
            ..Default::default()
        },
        PeakPickerHiRes {
            signal_to_noise: f64::NAN,
            ..Default::default()
        },
        PeakPickerHiRes {
            spacing_difference: -1.,
            ..Default::default()
        },
        PeakPickerHiRes {
            ms_levels: vec![0],
            ..Default::default()
        },
    ] {
        let mut data = valid.clone();
        assert!(picker.filter_spectrum(&mut data).is_err());
        assert_eq!(data, valid);
    }
    for input in [
        signal(&[0., 1., 1., 2., 3.], &[1., 2., 4., 2., 1.]),
        signal(&[0., 1., 3., 2., 4.], &[1., 2., 4., 2., 1.]),
        signal(&[0., 1., 2., 3., 4.], &[-1., 2., 4., 2., 1.]),
        signal(&[0., 1., 2., 3., f64::INFINITY], &[1., 2., 4., 2., 1.]),
    ] {
        assert!(PeakPickerHiRes::default().pick_spectrum(&input).is_err());
    }
    let no_work = SignalToNoiseEstimatorMedian {
        max_work: 1,
        ..Default::default()
    };
    assert!(no_work.estimate(&[0., 1.], &[1., 2.]).is_err());
    let no_bins = SignalToNoiseEstimatorMedian {
        max_bins: 2,
        ..Default::default()
    };
    assert!(no_bins.estimate(&[], &[]).is_err());
    assert!(CubicSpline2d::new(&[0., 0.], &[1., 2.]).is_err());
    assert!(CubicSpline2d::new(&[-f64::MAX, f64::MAX], &[1., 2.]).is_err());
    assert!(CubicSpline2d::new(&[0., 1e-300, 2e-300], &[1., 2., 1.]).is_err());
    let spline = CubicSpline2d::new(&[0., 1., 2.], &[1., 2., 1.]).unwrap();
    assert!(spline.eval(-1.).is_err());
    assert!(spline.eval(f64::NAN).is_err());
    assert!(spline.derivative(1., 4).is_err());
    assert!(spline.peak_maximum(0., 2., 0.).is_err());
    close(spline.peak_maximum(0., 2., 1e-6).unwrap().0, 1., 1e-12);
}

#[test]
fn auxiliary_selection_and_nonrepresentable_outputs_fail_atomically() {
    let valid = signal(&[0., 1., 2., 3., 4.], &[1., 2., 4., 2., 1.]);
    let missing = PeakPickerHiRes {
        ion_mobility_array: Some("custom mobility".into()),
        ..Default::default()
    };
    assert!(missing.pick_spectrum(&valid).is_err());
    let mut annotated = valid.clone();
    annotated
        .float_data_arrays
        .push(DataArray::new("custom mobility", vec![2.; 5]));
    assert_eq!(
        missing
            .pick_spectrum(&annotated)
            .unwrap()
            .spectrum
            .float_data_arrays[0]
            .data,
        [2.]
    );
    annotated.float_data_arrays[0].data[0] = f32::NAN;
    assert!(missing.pick_spectrum(&annotated).is_err());
    let mut duplicate = valid.clone();
    duplicate.float_data_arrays = vec![
        DataArray::new("Ion Mobility", vec![2.; 5]),
        DataArray::new("Ion Mobility", vec![2.; 5]),
    ];
    assert!(
        PeakPickerHiRes::default()
            .pick_spectrum(&duplicate)
            .is_err()
    );
    let mut overflow = signal(
        &[0., 1., 2., 3., 4.],
        &[0., f32::MAX / 2., f32::MAX, f32::MAX * 0.75, 0.],
    );
    let before = overflow.clone();
    assert!(
        PeakPickerHiRes::default()
            .filter_spectrum(&mut overflow)
            .is_err()
    );
    assert_eq!(overflow, before);
    let mut zero_center = signal(&[-2., -1., 0., 1., 2.], &[1., 2., 4., 2., 1.]);
    let before = zero_center.clone();
    assert!(
        PeakPickerHiRes {
            report_fwhm: Some(FwhmUnit::Ppm),
            ..Default::default()
        }
        .filter_spectrum(&mut zero_center)
        .is_err()
    );
    assert_eq!(zero_center, before);
}

#[test]
fn each_compatibility_flag_lifts_only_its_own_refusal() {
    let base = PickingCompatibility::default();
    let cases = [
        (
            signal(&[0., 1., 1., 2., 3., 4.], &[1., 2., 3., 5., 2., 1.]),
            PickingCompatibility {
                allow_duplicate_positions: true,
                ..base
            },
        ),
        (
            signal(&[0., 1., 3., 2., 4.], &[1., 2., 4., 2., 1.]),
            PickingCompatibility {
                allow_unsorted_positions: true,
                ..base
            },
        ),
        (
            signal(&[0., 1., 2., 3., 4.], &[-1., 2., 4., 2., 1.]),
            PickingCompatibility {
                allow_negative_intensities: true,
                ..base
            },
        ),
    ];
    for (input, flag) in cases {
        assert!(PeakPickerHiRes::default().pick_spectrum(&input).is_err());
        let lifted = PeakPickerHiRes {
            compatibility: flag,
            ..Default::default()
        };
        assert!(lifted.pick_spectrum(&input).is_ok(), "{flag:?}");
        // Every other flag together does not lift it.
        let others = PeakPickerHiRes {
            compatibility: PickingCompatibility {
                allow_duplicate_positions: !flag.allow_duplicate_positions,
                allow_unsorted_positions: !flag.allow_unsorted_positions,
                allow_negative_intensities: !flag.allow_negative_intensities,
                allow_nonpositive_maximum: !flag.allow_nonpositive_maximum,
                allow_nonpositive_fwhm_position: !flag.allow_nonpositive_fwhm_position,
                source_mobility_arrays: !flag.source_mobility_arrays,
            },
            ..Default::default()
        };
        assert!(others.pick_spectrum(&input).is_err(), "{flag:?}");
    }
    // A non-positive spline maximum needs its own flag on top of negatives.
    let negative = signal(
        &[0., 1., 2., 3., 4., 5., 6.],
        &[-100., -90., -50., -20., -50., -90., -100.],
    );
    let negatives_only = PeakPickerHiRes {
        compatibility: PickingCompatibility {
            allow_negative_intensities: true,
            ..base
        },
        ..Default::default()
    };
    assert!(negatives_only.pick_spectrum(&negative).is_err());
    let maximum = PeakPickerHiRes {
        compatibility: PickingCompatibility {
            allow_negative_intensities: true,
            allow_nonpositive_maximum: true,
            ..base
        },
        ..Default::default()
    };
    let out = maximum.pick_spectrum(&negative).unwrap();
    assert_eq!(out.spectrum.len(), 1);
    assert!(out.spectrum.peaks[0].intensity < 0.0);
    // With FWHM reporting the source's half-height loop never terminates for a
    // negative maximum; the port reports that instead of looping.
    let fwhm = PeakPickerHiRes {
        report_fwhm: Some(FwhmUnit::Absolute),
        ..maximum
    };
    assert!(fwhm.pick_spectrum(&negative).is_err());
    // A ppm width at a non-positive position needs its own flag.
    let mirrored = signal(
        &[-102., -101., -100., -99., -98.],
        &[200., 250., 450., 250., 200.],
    );
    let ppm = PeakPickerHiRes {
        report_fwhm: Some(FwhmUnit::Ppm),
        ..Default::default()
    };
    assert!(ppm.pick_spectrum(&mirrored).is_err());
    let lifted = PeakPickerHiRes {
        compatibility: PickingCompatibility {
            allow_nonpositive_fwhm_position: true,
            ..base
        },
        ..ppm
    };
    let width = lifted
        .pick_spectrum(&mirrored)
        .unwrap()
        .spectrum
        .float_data_arrays[0]
        .data[0];
    assert!(width < 0.0);
    // The first of two ion mobility arrays is used only with the source flag.
    let mut two = signal(
        &[100., 100.01, 100.02, 100.03, 100.04],
        &[200., 250., 450., 250., 200.],
    );
    two.float_data_arrays = vec![
        DataArray::new("raw ion mobility array", vec![1.; 5]),
        DataArray::new("Ion Mobility", vec![2.; 5]),
    ];
    assert!(PeakPickerHiRes::default().pick_spectrum(&two).is_err());
    let first = PeakPickerHiRes {
        compatibility: PickingCompatibility {
            source_mobility_arrays: true,
            ..base
        },
        ..Default::default()
    }
    .pick_spectrum(&two)
    .unwrap();
    assert_eq!(first.spectrum.float_data_arrays.len(), 1);
    assert_eq!(
        first.spectrum.float_data_arrays[0].name,
        "raw ion mobility array"
    );
    assert_eq!(first.spectrum.float_data_arrays[0].data, [1.]);
    assert_eq!(first.omitted_arrays, ["Ion Mobility"]);
    assert_eq!(
        PickingCompatibility::source(),
        PickingCompatibility {
            allow_duplicate_positions: true,
            allow_unsorted_positions: true,
            allow_negative_intensities: true,
            allow_nonpositive_maximum: true,
            allow_nonpositive_fwhm_position: true,
            source_mobility_arrays: true,
        }
    );
}

#[test]
fn noise_estimator_source_edge_cases() {
    // A negative automatic range: the source returns before estimating and
    // leaves every ratio at zero.
    let estimator = SignalToNoiseEstimatorMedian::default();
    let x = [0., 1., 2., 3., 4.];
    let y = [-100., -99., -95., -98., -100.];
    assert!(estimator.estimate(&x, &y).is_err());
    let out = estimator
        .estimate_with_compatibility(&x, &y, &PickingCompatibility::source())
        .unwrap();
    assert!(out.max_intensity < 0.0);
    assert_eq!(out.signal_to_noise, [0.; 5]);
    assert!(out.signal_to_noise.iter().all(|v| v.is_sign_positive()));
    assert_eq!(out.noise, [f64::INFINITY; 5]);
    assert_eq!(out.sparse_window_percent, 0.0);
    // Percentages use the source's `count * 100 / n`.
    let sparse = SignalToNoiseEstimatorMedian {
        window_length: 1.0,
        min_required_elements: 2,
        ..Default::default()
    };
    let out = sparse
        .estimate(&[0., 10., 20., 20.5, 30., 40., 50.], &[1.; 7])
        .unwrap();
    assert_eq!(out.sparse_window_percent, 5.0 * 100.0 / 7.0);
    // Manual mode refuses a non-positive upper end only when estimating; the
    // picker never estimates with signal_to_noise = 0.
    let manual = SignalToNoiseEstimatorMedian {
        histogram_range: NoiseHistogramRange::Manual {
            max_intensity: -1.0,
        },
        ..Default::default()
    };
    assert!(manual.estimate(&[0., 1.], &[1., 2.]).is_err());
    let profile = signal(
        &[100., 100.01, 100.02, 100.03, 100.04],
        &[200., 250., 450., 250., 200.],
    );
    let quiet = PeakPickerHiRes {
        noise_estimator: manual.clone(),
        ..Default::default()
    };
    assert_eq!(quiet.pick_spectrum(&profile).unwrap().spectrum.len(), 1);
    assert!(
        PeakPickerHiRes {
            signal_to_noise: 1.0,
            ..quiet
        }
        .pick_spectrum(&profile)
        .is_err()
    );
}

#[test]
fn spacing_constraints_accept_infinity_as_disabled() {
    // The source maps zero to infinity; an explicit infinity behaves the same.
    let input = signal(&[0., 3., 6., 7., 8.], &[200., 250., 450., 250., 200.]);
    let zero = PeakPickerHiRes {
        spacing_difference: 0.0,
        spacing_difference_gap: 0.0,
        ..Default::default()
    };
    let infinite = PeakPickerHiRes {
        spacing_difference: f64::INFINITY,
        spacing_difference_gap: f64::INFINITY,
        ..Default::default()
    };
    assert_eq!(
        zero.pick_spectrum(&input).unwrap(),
        infinite.pick_spectrum(&input).unwrap()
    );
    for invalid in [f64::NAN, -0.5, f64::NEG_INFINITY] {
        let picker = PeakPickerHiRes {
            spacing_difference_gap: invalid,
            ..Default::default()
        };
        assert!(picker.pick_spectrum(&input).is_err());
    }
}
