// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::kernel::{DataArray, SpectrumType};
use openms::processing::{
    SpectrumFilter,
    smoothing::{GaussFilter, GaussFilterAlgorithm, GaussianWidth, SavitzkyGolayFilter},
};
use openms::{ChromatogramPeak, Error, MSChromatogram, MSExperiment, MSSpectrum, Peak1D};

fn fixture(text: &str) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let rows: Vec<Vec<f64>> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            line.split_whitespace()
                .map(|v| v.parse().unwrap())
                .collect()
        })
        .collect();
    (
        rows.iter().map(|r| r[0]).collect(),
        rows.iter().map(|r| r[1]).collect(),
        rows.iter().map(|r| r[2]).collect(),
    )
}
fn near(actual: &[f64], expected: &[f64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (index, (&a, &e)) in actual.iter().zip(expected).enumerate() {
        assert!((a - e).abs() <= tolerance, "index {index}: {a} != {e}");
    }
}
fn spectrum(x: &[f64], y: &[f64]) -> MSSpectrum {
    let mut spectrum = MSSpectrum::from_peaks(
        x.iter()
            .zip(y)
            .map(|(&x, &y)| Peak1D::new(x, y as f32))
            .collect(),
    );
    spectrum.name = "preserved".into();
    spectrum.rt = 42.0;
    spectrum.spectrum_type = SpectrumType::Centroid;
    spectrum.metadata.insert("note".into(), "untouched".into());
    spectrum
        .integer_data_arrays
        .push(DataArray::new("index", (0..x.len() as i32).collect()));
    spectrum
}
fn chromatogram(x: &[f64], y: &[f64]) -> MSChromatogram {
    let mut chromatogram = MSChromatogram::from_peaks(
        x.iter()
            .zip(y)
            .map(|(&x, &y)| ChromatogramPeak::new(x, y as f32))
            .collect(),
    );
    chromatogram.name = "chromatogram".into();
    chromatogram
        .string_data_arrays
        .push(DataArray::new("label", vec!["kept".into(); x.len()]));
    chromatogram
}
fn signal(spectrum: &MSSpectrum) -> Vec<f64> {
    spectrum
        .peaks
        .iter()
        .map(|p| f64::from(p.intensity))
        .collect()
}

#[test]
fn gaussian_matches_upstream_numeric_golden_for_spectra_and_chromatograms() {
    let (x, y, expected) = fixture(include_str!("data/smoothing_gaussian_golden.tsv"));
    let algorithm = GaussFilterAlgorithm::new(GaussianWidth::Absolute(0.2), 0.01).unwrap();
    let output = algorithm.filter(&x, &y).unwrap();
    assert!(output.found_signal);
    near(&output.intensities, &expected, 1e-6);
    let mut spectrum = spectrum(&x, &y);
    let before = spectrum.clone();
    GaussFilter::default()
        .filter_spectrum(&mut spectrum)
        .unwrap();
    near(&signal(&spectrum), &expected, 1e-6);
    assert_eq!(spectrum.spectrum_type, SpectrumType::Profile);
    assert_eq!(spectrum.name, before.name);
    assert_eq!(spectrum.metadata, before.metadata);
    assert_eq!(spectrum.integer_data_arrays, before.integer_data_arrays);
    assert_eq!(spectrum.peaks.iter().map(|p| p.mz).collect::<Vec<_>>(), x);
    let mut chromatogram = chromatogram(&x, &y);
    let labels = chromatogram.string_data_arrays.clone();
    GaussFilter::default()
        .filter_chromatogram(&mut chromatogram)
        .unwrap();
    near(
        &chromatogram
            .peaks
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect::<Vec<_>>(),
        &expected,
        1e-6,
    );
    assert_eq!(chromatogram.string_data_arrays, labels);
}

#[test]
fn gaussian_preserves_upstream_constant_signal_and_strict_boundary_rule() {
    let x: Vec<_> = (0..5).map(|i| 500.0 + i as f64 * 0.2).collect();
    let algorithm = GaussFilterAlgorithm::new(GaussianWidth::Absolute(8.0), 0.01).unwrap();
    near(
        &algorithm.filter(&x, &[1.0; 5]).unwrap().intensities,
        &[1.0; 5],
        1e-12,
    );
    // The only intervals touch a global endpoint and are excluded at the
    // center. At an endpoint, only its inward interval is integrated.
    let algorithm = GaussFilterAlgorithm::new(GaussianWidth::Absolute(8.0), 0.25).unwrap();
    let result = algorithm
        .filter(&[0.0, 1.0, 2.0], &[2.0, 2.0, 2.0])
        .unwrap();
    near(&result.intensities, &[2.0, 0.0, 2.0], 1e-12);
    assert_eq!(
        algorithm
            .filter(&[0.0, 1.0], &[2.0, 2.0])
            .unwrap()
            .intensities,
        [0.0, 0.0]
    );
}

#[test]
fn gaussian_tabulated_trapezoids_handle_nonuniform_intervals() {
    let algorithm = GaussFilterAlgorithm::new(GaussianWidth::Absolute(8.0), 0.25).unwrap();
    let positions = [0.0, 0.125, 0.75, 1.25, 2.0];
    let values = [0.0, 2.0, 3.0, 4.0, 0.0];
    let result = algorithm.filter(&positions, &values).unwrap();
    // Center integrates intervals [0.125,0.75] and [0.75,1.25];
    // kernel at distance 0.625 interpolates coefficients at 0.5 and 0.75.
    let a = ((-0.5_f64.powi(2) / 2.0).exp() + (-0.75_f64.powi(2) / 2.0).exp()) / 2.0;
    let b = (-0.5_f64.powi(2) / 2.0).exp();
    let expected = (0.625 / 2.0 * (2.0 * a + 3.0) + 0.5 / 2.0 * (3.0 + 4.0 * b))
        / (0.625 / 2.0 * (a + 1.0) + 0.5 / 2.0 * (1.0 + b));
    near(&[result.intensities[2]], &[expected], 1e-12);
}

#[test]
fn gaussian_ppm_width_is_recalculated_per_position() {
    let x: Vec<_> = (0..9).map(|i| 500.0 + i as f64 * 0.03).collect();
    let y = [0.0, 0.0, 0.0, 1.0, 0.8, 1.2, 0.0, 0.0, 0.0];
    let ppm = GaussFilterAlgorithm::new(GaussianWidth::Ppm(400.0), 0.01).unwrap();
    let result = ppm.filter(&x, &y).unwrap();
    for i in 0..x.len() {
        let absolute =
            GaussFilterAlgorithm::new(GaussianWidth::Absolute(400.0 / 1e6 * x[i]), 0.01).unwrap();
        near(
            &[result.intensities[i]],
            &[absolute.filter(&x, &y).unwrap().intensities[i]],
            1e-12,
        );
    }
    let mut c = chromatogram(&x, &y);
    let before = c.clone();
    assert!(
        GaussFilter::new(GaussianWidth::Ppm(10.0))
            .unwrap()
            .filter_chromatogram(&mut c)
            .is_err()
    );
    assert_eq!(c, before);
}

#[test]
fn gaussian_wrapper_preserves_all_zero_results_for_three_or_more_points() {
    let filter = GaussFilter::new(GaussianWidth::Absolute(0.01)).unwrap();
    let mut s = spectrum(&[1.0, 2.0, 3.0], &[1.0, 2.0, 3.0]);
    filter.filter_spectrum(&mut s).unwrap();
    assert_eq!(signal(&s), [1.0, 2.0, 3.0]);
    assert_eq!(s.spectrum_type, SpectrumType::Profile);
    let mut single = spectrum(&[1.0], &[7.0]);
    filter.filter_spectrum(&mut single).unwrap();
    assert_eq!(signal(&single), [0.0]);
    let algorithm = GaussFilterAlgorithm::default();
    assert!(!algorithm.filter(&[], &[]).unwrap().found_signal);
    assert_eq!(
        algorithm
            .filter(&[1.0, 1.0, 1.0], &[2.0, 2.0, 2.0])
            .unwrap()
            .intensities,
        [0.0, 0.0, 0.0]
    );
}

#[test]
fn gaussian_errors_leave_spectra_and_whole_experiments_unchanged() {
    let mut s = spectrum(&[2.0, 1.0, 3.0], &[1.0, 2.0, 3.0]);
    let before = s.clone();
    assert!(matches!(
        GaussFilter::default().filter_spectrum(&mut s),
        Err(Error::UnsortedData)
    ));
    assert_eq!(s, before);
    for width in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e-200] {
        let filter = GaussFilter {
            algorithm: GaussFilterAlgorithm {
                width: GaussianWidth::Absolute(width),
                ..Default::default()
            },
        };
        let mut input = spectrum(&[1.0, 2.0, 3.0], &[1.0, 2.0, 3.0]);
        let before = input.clone();
        assert!(filter.filter_spectrum(&mut input).is_err());
        assert_eq!(input, before);
    }
    let mut experiment = MSExperiment {
        spectra: vec![spectrum(&[1.0, 1.03, 1.06], &[1.0, 2.0, 3.0]), before],
        ..Default::default()
    };
    let original = experiment.clone();
    assert!(
        GaussFilter::default()
            .filter_experiment(&mut experiment)
            .is_err()
    );
    assert_eq!(experiment, original);
    let limited = GaussFilterAlgorithm {
        max_coefficients: 2,
        ..Default::default()
    };
    assert!(limited.filter(&[1.0], &[1.0]).is_err());
    assert!(GaussFilterAlgorithm::default().filter(&[1.0], &[]).is_err());
    assert!(
        GaussFilterAlgorithm::default()
            .filter(&[1.0], &[f64::NAN])
            .is_err()
    );
}

#[test]
fn savitzky_golay_matches_upstream_golden_with_asymmetric_edges_and_even_rounding() {
    let (x, y, expected) = fixture(include_str!("data/smoothing_savgol_golden.tsv"));
    let filter = SavitzkyGolayFilter::new(4, 2).unwrap();
    assert_eq!(filter.frame_length(), 5);
    assert_eq!(filter.polynomial_order(), 2);
    near(&filter.filter(&x, &y).unwrap(), &expected, 4e-6);
    let mut s = spectrum(&x, &y);
    let before = s.clone();
    filter.filter_spectrum(&mut s).unwrap();
    near(&signal(&s), &expected, 4e-6);
    assert_eq!(s.integer_data_arrays, before.integer_data_arrays);
    assert_eq!(s.metadata, before.metadata);
    assert_eq!(s.spectrum_type, SpectrumType::Centroid);
    let mut c = chromatogram(&x, &y);
    let labels = c.string_data_arrays.clone();
    filter.filter_chromatogram(&mut c).unwrap();
    near(
        &c.peaks
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect::<Vec<_>>(),
        &expected,
        4e-6,
    );
    assert_eq!(c.string_data_arrays, labels);
}

#[test]
fn savitzky_golay_reproduces_polynomials_including_edges() {
    for (frame, degree) in [(3, 2), (5, 2), (11, 4), (31, 8)] {
        let filter = SavitzkyGolayFilter::new(frame, degree).unwrap();
        let x: Vec<_> = (0..51).map(|i| i as f64).collect();
        let y: Vec<_> = x.iter().map(|v| 2.0 + 0.05 * v + 0.001 * v * v).collect();
        near(&filter.filter(&x, &y).unwrap(), &y, 1e-10);
    }
    let filter = SavitzkyGolayFilter::new(3, 2).unwrap();
    near(
        &filter
            .filter(&[0.0; 5], &[0.0, 0.0, 1.0, 0.0, 0.0])
            .unwrap(),
        &[0.0, 0.0, 1.0, 0.0, 0.0],
        1e-12,
    );
}

#[test]
fn savitzky_golay_short_inputs_clamps_and_validation() {
    let filter = SavitzkyGolayFilter::new(5, 2).unwrap();
    assert_eq!(
        filter.filter(&[0.0, 1.0], &[-2.0, 3.0]).unwrap(),
        [-2.0, 3.0]
    );
    near(
        &filter
            .filter(&[0.0, 1.0, 2.0, 3.0, 4.0], &[-1.0; 5])
            .unwrap(),
        &[0.0; 5],
        1e-12,
    );
    for (frame, degree) in [(1, 1), (5, 5), (1025, 2), (usize::MAX, 2), (33, 33)] {
        assert!(SavitzkyGolayFilter::new(frame, degree).is_err());
    }
    assert_eq!(SavitzkyGolayFilter::new(0, 0).unwrap().frame_length(), 1);
    let mut invalid = spectrum(&[2.0, 1.0], &[1.0, 2.0]);
    let original = invalid.clone();
    assert!(filter.filter_spectrum(&mut invalid).is_err());
    assert_eq!(invalid, original);
}

#[test]
fn smoothing_experiments_include_chromatograms_and_are_atomic() {
    let (x, y, expected) = fixture(include_str!("data/smoothing_savgol_golden.tsv"));
    let filter = SavitzkyGolayFilter::new(5, 2).unwrap();
    let mut e = MSExperiment {
        spectra: vec![spectrum(&x, &y)],
        chromatograms: vec![chromatogram(&x, &y)],
        ..Default::default()
    };
    filter.filter_experiment(&mut e).unwrap();
    near(&signal(&e.spectra[0]), &expected, 4e-6);
    near(
        &e.chromatograms[0]
            .peaks
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect::<Vec<_>>(),
        &expected,
        4e-6,
    );
    let mut e = MSExperiment {
        spectra: vec![spectrum(&x, &y)],
        chromatograms: vec![chromatogram(&[1.0, 0.0], &[1.0, 2.0])],
        ..Default::default()
    };
    let before = e.clone();
    assert!(filter.filter_experiment(&mut e).is_err());
    assert_eq!(e, before);
    let ppm = GaussFilter::new(GaussianWidth::Ppm(400.0)).unwrap();
    assert!(ppm.filter_experiment(&mut e).is_err());
    assert_eq!(e, before);
}
