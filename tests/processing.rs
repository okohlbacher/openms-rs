// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::format::dta;
use openms::kernel::DataArray;
use openms::processing::*;
use openms::{ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D};

fn fixture() -> MSSpectrum {
    dta::read(&include_bytes!("data/Transformers_tests.dta")[..]).unwrap()
}
fn spectrum(values: &[f32]) -> MSSpectrum {
    MSSpectrum::from_peaks(
        values
            .iter()
            .enumerate()
            .map(|(i, &x)| Peak1D::new(i as f64, x))
            .collect(),
    )
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-5, "{a} != {b}");
}

#[test]
fn upstream_normalizer_golden() {
    // Normalizer_test.cpp: the original 121-peak DTA has base intensity 46.
    let mut spec = fixture();
    assert_eq!(spec.len(), 121);
    assert_eq!(spec.base_peak().unwrap().intensity, 46.0);
    Normalizer::default().filter_spectrum(&mut spec).unwrap();
    assert_eq!(spec.base_peak().unwrap().intensity, 1.0);
    Normalizer {
        method: NormalizationMethod::ToTic,
    }
    .filter_spectrum(&mut spec)
    .unwrap();
    close(spec.peaks.iter().map(|p| f64::from(p.intensity)).sum(), 1.0);
}

#[test]
fn upstream_threshold_and_largest_golden() {
    // ThresholdMower.cpp sets a double threshold of 0.05. Check the adjacent
    // f32 observations on either side and preserve equality at explicit zero.
    let below = f32::from_bits(0.05_f32.to_bits() - 1);
    let mut boundary = spectrum(&[0.0, below, 0.05, 1.0]);
    ThresholdMower::default()
        .filter_spectrum(&mut boundary)
        .unwrap();
    assert_eq!(
        boundary.peaks,
        [Peak1D::new(2.0, 0.05), Peak1D::new(3.0, 1.0)]
    );
    let mut zero = spectrum(&[0.0]);
    ThresholdMower { threshold: 0.0 }
        .filter_spectrum(&mut zero)
        .unwrap();
    assert_eq!(zero.len(), 1);
    let mut spec = fixture();
    ThresholdMower { threshold: 1.0 }
        .filter_spectrum(&mut spec)
        .unwrap();
    assert_eq!(spec.len(), 121);
    ThresholdMower { threshold: 10.0 }
        .filter_spectrum(&mut spec)
        .unwrap();
    assert_eq!(spec.len(), 14);
    NLargest { n: 10 }.filter_spectrum(&mut spec).unwrap();
    assert_eq!(spec.len(), 10);
    assert!(
        spec.peaks
            .windows(2)
            .all(|p| p[0].intensity >= p[1].intensity)
    );
}

#[test]
fn upstream_rank_golden() {
    // RankScaler_test.cpp: dense descending ranks offset by peak count.
    let mut spec = fixture();
    RankScaler.filter_spectrum(&mut spec).unwrap();
    assert_eq!(spec.peaks[0].intensity, 96.0);
    assert_eq!(spec.peaks[120].intensity, 121.0);
    close(spec.peaks[120].mz, 136.0765);
}

#[test]
fn upstream_nlargest_preserves_triangle_annotations() {
    // NLargest_test.cpp's triangle-shaped 100-peak regression.
    let mut spec = MSSpectrum::from_peaks(
        (0..100)
            .map(|i| {
                Peak1D::new(
                    i as f64,
                    if i < 50 {
                        i as f32 + 0.1
                    } else {
                        (100 - i) as f32 + 0.2
                    },
                )
            })
            .collect(),
    );
    spec.integer_data_arrays
        .push(DataArray::new("index", (0..100).collect()));
    spec.string_data_arrays.push(DataArray::new(
        "slope",
        (0..100)
            .map(|i| if i < 50 { "up".into() } else { "down".into() })
            .collect(),
    ));
    NLargest { n: 10 }.filter_spectrum(&mut spec).unwrap();
    assert_eq!(
        spec.integer_data_arrays[0].data,
        [50, 51, 49, 52, 48, 53, 47, 54, 46, 55]
    );
    assert_eq!(
        &spec.string_data_arrays[0].data[..3],
        ["down", "down", "up"]
    );
    close(f64::from(spec.peaks[9].intensity), 45.2);
}

#[test]
fn normalization_zero_and_signed_policy() {
    let mut zero = spectrum(&[0.0, 0.0]);
    Normalizer::default().filter_spectrum(&mut zero).unwrap();
    assert_eq!(zero, spectrum(&[0.0, 0.0]));
    let mut signed = spectrum(&[-1.0, 1.0]);
    let before = signed.clone();
    assert!(
        Normalizer {
            method: NormalizationMethod::ToTic
        }
        .filter_spectrum(&mut signed)
        .is_err()
    );
    assert_eq!(before, signed);
    let mut negative = spectrum(&[-4.0, -2.0]);
    Normalizer::default()
        .filter_spectrum(&mut negative)
        .unwrap();
    assert_eq!(negative, spectrum(&[2.0, 1.0]));
}

#[test]
fn nonfinite_and_overflow_errors_do_not_mutate() {
    let mut spec = spectrum(&[1.0, 2.0]);
    let before = spec.clone();
    assert!(
        ThresholdMower {
            threshold: f64::NAN
        }
        .filter_spectrum(&mut spec)
        .is_err()
    );
    assert_eq!(before, spec);
    let mut overflow = spectrum(&[-f32::MAX, f32::MIN_POSITIVE]);
    let before = overflow.clone();
    assert!(
        Normalizer::default()
            .filter_spectrum(&mut overflow)
            .is_err()
    );
    assert_eq!(before, overflow);
}

#[test]
fn sqrt_negative_clamping_and_dense_rank_ties() {
    let mut spec = spectrum(&[-9.0, 0.0, 4.0, 9.0]);
    SqrtScaler.filter_spectrum(&mut spec).unwrap();
    assert_eq!(spec, spectrum(&[0.0, 0.0, 2.0, 3.0]));
    RankScaler.filter_spectrum(&mut spec).unwrap();
    assert_eq!(spec, spectrum(&[2.0, 2.0, 3.0, 4.0]));
    let mut zero = spectrum(&[0.0, 0.0]);
    RankScaler.filter_spectrum(&mut zero).unwrap();
    assert_eq!(zero, spectrum(&[3.0, 3.0]));
}

#[test]
fn largest_noop_and_empty_selection() {
    let mut spec = spectrum(&[1.0, 3.0, 2.0]);
    NLargest { n: 3 }.filter_spectrum(&mut spec).unwrap();
    assert_eq!(spec, spectrum(&[1.0, 3.0, 2.0]));
    spec.float_data_arrays
        .push(DataArray::new("signal", vec![10.0, 30.0, 20.0]));
    NLargest { n: 0 }.filter_spectrum(&mut spec).unwrap();
    assert!(spec.is_empty());
    assert!(spec.float_data_arrays[0].data.is_empty());
}

#[test]
fn experiment_transform_is_atomic() {
    let mut exp = MSExperiment {
        spectra: vec![spectrum(&[1.0, 2.0]), spectrum(&[-1.0, 1.0])],
        ..Default::default()
    };
    let before = exp.clone();
    assert!(
        Normalizer {
            method: NormalizationMethod::ToTic
        }
        .filter_experiment(&mut exp)
        .is_err()
    );
    assert_eq!(before, exp);
}

fn resample_fixture() -> MSSpectrum {
    MSSpectrum::from_peaks(vec![
        Peak1D::new(0.0, 3.0),
        Peak1D::new(0.5, 6.0),
        Peak1D::new(1.0, 8.0),
        Peak1D::new(1.6, 2.0),
        Peak1D::new(1.8, 1.0),
    ])
}

#[test]
fn upstream_linear_resampler_golden() {
    // LinearResamplerAlign_test.cpp: five raw points -> grid spacing 0.75.
    let mut spec = resample_fixture();
    spec.name = "kept".into();
    LinearResamplerAlign::new(0.75)
        .unwrap()
        .raster(&mut spec)
        .unwrap();
    assert_eq!(
        spec.peaks.iter().map(|p| p.mz).collect::<Vec<_>>(),
        [0.0, 0.75, 1.5, 2.25]
    );
    for (p, expected) in spec.peaks.iter().zip([
        5.0,
        4.0 + 16.0 / 3.0,
        8.0 / 3.0 + 2.0 + 1.0 / 3.0,
        2.0 / 3.0,
    ]) {
        close(f64::from(p.intensity), expected);
    }
    close(
        spec.peaks.iter().map(|p| f64::from(p.intensity)).sum(),
        20.0,
    );
    assert_eq!(spec.name, "kept");
}

#[test]
fn resample_chromatogram_and_boundaries() {
    let mut chrom = MSChromatogram {
        peaks: resample_fixture()
            .peaks
            .iter()
            .map(|p| ChromatogramPeak::new(p.mz, p.intensity))
            .collect(),
        ..Default::default()
    };
    LinearResamplerAlign::new(0.75)
        .unwrap()
        .raster_chromatogram(&mut chrom)
        .unwrap();
    close(
        chrom.peaks.iter().map(|p| f64::from(p.intensity)).sum(),
        20.0,
    );
    let output = resample_to_grid(
        &[
            Peak1D::new(-1.0, 2.0),
            Peak1D::new(0.5, 4.0),
            Peak1D::new(9.0, 3.0),
        ],
        &[0.0, 1.0],
    )
    .unwrap();
    assert_eq!(output, [Peak1D::new(0.0, 4.0), Peak1D::new(1.0, 5.0)]);
    let output = resample_to_grid(&resample_fixture().peaks, &[1.0]).unwrap();
    assert_eq!(output, [Peak1D::new(1.0, 20.0)]);
}

#[test]
fn aligned_grid_excludes_outside_peaks() {
    let mut spec = resample_fixture();
    LinearResamplerAlign::new(0.5)
        .unwrap()
        .raster_align(&mut spec, 0.5, 1.0)
        .unwrap();
    assert_eq!(spec.peaks, [Peak1D::new(0.5, 6.0), Peak1D::new(1.0, 8.0)]);
}

#[test]
fn resample_rejects_invalid_grids_and_preserves_arrays() {
    for spacing in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(LinearResamplerAlign::new(spacing).is_err());
    }
    for grid in [
        vec![],
        vec![1.0, 1.0],
        vec![2.0, 1.0],
        vec![f64::NAN],
        vec![-f64::MAX, f64::MAX],
    ] {
        assert!(resample_to_grid(&[], &grid).is_err());
    }
    let mut spec = resample_fixture();
    spec.integer_data_arrays
        .push(DataArray::new("index", vec![0, 1, 2, 3, 4]));
    let before = spec.clone();
    assert!(LinearResamplerAlign::default().raster(&mut spec).is_err());
    assert_eq!(before, spec);
    let mut spec = resample_fixture();
    assert!(
        LinearResamplerAlign {
            spacing: 1e-9,
            max_points: 10
        }
        .raster(&mut spec)
        .is_err()
    );
    assert!(
        resample_to_grid(&[Peak1D::new(2.0, 1.0), Peak1D::new(1.0, 1.0)], &[1.0, 2.0]).is_err()
    );
}

#[test]
fn empty_spectra_are_valid_filter_inputs() {
    let filters: Vec<Box<dyn SpectrumFilter>> = vec![
        Box::new(Normalizer::default()),
        Box::new(ThresholdMower::default()),
        Box::new(NLargest::default()),
        Box::new(SqrtScaler),
        Box::new(RankScaler),
    ];
    for filter in filters {
        let mut spec = MSSpectrum::default();
        filter.filter_spectrum(&mut spec).unwrap();
        assert!(spec.is_empty());
    }
}
