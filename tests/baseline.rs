// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::kernel::{DataArray, SpectrumType};
use openms::processing::{
    SpectrumFilter,
    baseline::{MorphologicalFilter, MorphologicalMethod as Method, StructuringElement as Element},
};
use openms::{ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D};

fn golden() -> Vec<Vec<f32>> {
    include_str!("data/baseline_morphological_golden.txt")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            line.split_whitespace()
                .map(|v| v.parse().unwrap())
                .collect()
        })
        .collect()
}
fn spectrum(values: &[f32], spacing: f64) -> MSSpectrum {
    let mut s = MSSpectrum::from_peaks(
        values
            .iter()
            .enumerate()
            .map(|(i, &y)| Peak1D::new(i as f64 * spacing, y))
            .collect(),
    );
    s.name = "preserved".into();
    s.metadata.insert("comment".into(), "untouched".into());
    s.spectrum_type = SpectrumType::Centroid;
    s.integer_data_arrays
        .push(DataArray::new("index", (0..values.len() as i32).collect()));
    s
}
fn values(s: &MSSpectrum) -> Vec<f32> {
    s.peaks.iter().map(|p| p.intensity).collect()
}

#[test]
fn all_operations_match_the_pinned_upstream_golden_table() {
    let rows = golden();
    let input: Vec<_> = rows.iter().map(|r| r[0]).collect();
    for (method, column) in [
        (Method::Identity, 0),
        (Method::Erosion, 1),
        (Method::Opening, 2),
        (Method::Dilation, 3),
        (Method::Closing, 4),
        (Method::Gradient, 5),
        (Method::TopHat, 6),
        (Method::BottomHat, 7),
        (Method::ErosionSimple, 1),
        (Method::DilationSimple, 3),
    ] {
        let filter = MorphologicalFilter::new(method, Element::DataPoints(3)).unwrap();
        let expected: Vec<_> = rows.iter().map(|r| r[column]).collect();
        assert_eq!(filter.filter_range(&input).unwrap(), expected, "{method:?}");
        let mut s = spectrum(&input, 0.25);
        let original = s.clone();
        filter.filter_spectrum(&mut s).unwrap();
        assert_eq!(values(&s), expected);
        assert_eq!(s.spectrum_type, SpectrumType::Profile);
        assert_eq!(s.metadata, original.metadata);
        assert_eq!(s.integer_data_arrays, original.integer_data_arrays);
        assert_eq!(
            s.peaks.iter().map(|p| p.mz).collect::<Vec<_>>(),
            original.peaks.iter().map(|p| p.mz).collect::<Vec<_>>()
        );
    }
}

fn reference(values: &[f32], points: usize, maximum: bool) -> Vec<f32> {
    (0..values.len())
        .map(|i| {
            let start = i.saturating_sub(points / 2);
            let end = (i + points / 2 + 1).min(values.len());
            values[start..end]
                .iter()
                .copied()
                .reduce(|a, b| if maximum { a.max(b) } else { a.min(b) })
                .unwrap()
        })
        .collect()
}

#[test]
fn monotonic_and_irregular_signals_match_upstream_style_boundary_reference() {
    for size in 0..50 {
        for input in [
            (0..size)
                .map(|i| i as f32 - size as f32 / 2.0)
                .collect::<Vec<_>>(),
            (0..size).map(|i| size as f32 / 2.0 - i as f32).collect(),
            (0..size).map(|i| ((i * 17) % 11) as f32 - 5.0).collect(),
        ] {
            for points in (1..=2 * size + 3).step_by(2) {
                for (method, maximum) in [(Method::Erosion, false), (Method::Dilation, true)] {
                    let filter =
                        MorphologicalFilter::new(method, Element::DataPoints(points)).unwrap();
                    assert_eq!(
                        filter.filter_range(&input).unwrap(),
                        reference(&input, points, maximum),
                        "size={size},points={points}"
                    );
                }
            }
        }
    }
}

#[test]
fn thomson_conversion_uses_average_spacing_then_ceil_and_odd_rounding() {
    let input: Vec<_> = golden().iter().map(|r| r[0]).collect();
    for tenth in 5..=20 {
        let width = tenth as f64 / 10.0;
        let filter = MorphologicalFilter::new(Method::Dilation, Element::Thomson(width)).unwrap();
        let mut s = spectrum(&input, 0.25);
        let mut points = (width / 0.25).ceil() as usize;
        if points % 2 == 0 {
            points += 1;
        }
        assert_eq!(
            filter
                .effective_window(&s.peaks.iter().map(|p| p.mz).collect::<Vec<_>>())
                .unwrap(),
            points
        );
        filter.filter_spectrum(&mut s).unwrap();
        assert_eq!(values(&s), reference(&input, points, true));
    }
    let filter = MorphologicalFilter::new(Method::Erosion, Element::Thomson(1.0)).unwrap();
    assert_eq!(filter.effective_window(&[0.0, 0.1, 0.9, 1.0]).unwrap(), 3);
    assert_eq!(
        MorphologicalFilter::new(Method::Dilation, Element::DataPoints(4))
            .unwrap()
            .filter_range(&input)
            .unwrap(),
        reference(&input, 5, true)
    );
}

#[test]
fn range_and_spectrum_singletons_preserve_the_upstream_wrapper_difference() {
    let filter = MorphologicalFilter::new(Method::TopHat, Element::DataPoints(3)).unwrap();
    assert_eq!(filter.filter_range(&[9.0]).unwrap(), [0.0]);
    let mut s = spectrum(&[9.0], 1.0);
    filter.filter_spectrum(&mut s).unwrap();
    assert_eq!(values(&s), [9.0]);
    assert_eq!(s.spectrum_type, SpectrumType::Profile);
    assert!(filter.filter_range(&[]).unwrap().is_empty());
    let mut empty = MSSpectrum::new();
    filter.filter_spectrum(&mut empty).unwrap();
    assert_eq!(empty.spectrum_type, SpectrumType::Profile);
    let huge = MorphologicalFilter::new(Method::Dilation, Element::DataPoints(usize::MAX)).unwrap();
    assert_eq!(huge.filter_range(&[1.0, 5.0, 3.0]).unwrap(), [5.0; 3]);
}

#[test]
fn signed_bottom_hat_and_opening_order_are_preserved() {
    let input = [0.0, 3.0, 0.0];
    let opening = MorphologicalFilter::new(Method::Opening, Element::DataPoints(3)).unwrap();
    let closing = MorphologicalFilter::new(Method::Closing, Element::DataPoints(3)).unwrap();
    let bottom = MorphologicalFilter::new(Method::BottomHat, Element::DataPoints(3)).unwrap();
    assert_eq!(opening.filter_range(&input).unwrap(), [0.0; 3]);
    assert_eq!(closing.filter_range(&input).unwrap(), [3.0; 3]);
    assert_eq!(bottom.filter_range(&input).unwrap(), [-3.0, 0.0, -3.0]);
}

#[test]
fn chromatogram_extension_uses_seconds_and_preserves_annotations() {
    let mut c = MSChromatogram::from_peaks(
        [1.0, 4.0, 1.0, 1.0]
            .into_iter()
            .enumerate()
            .map(|(i, y)| ChromatogramPeak::new(i as f64 * 0.5, y))
            .collect(),
    );
    c.name = "TIC".into();
    c.string_data_arrays
        .push(DataArray::new("labels", vec!["a".into(); 4]));
    let before = c.clone();
    let filter = MorphologicalFilter::new(Method::TopHat, Element::Seconds(1.0)).unwrap();
    filter.filter_chromatogram(&mut c).unwrap();
    assert_eq!(
        c.peaks.iter().map(|p| p.intensity).collect::<Vec<_>>(),
        [0.0, 3.0, 0.0, 0.0]
    );
    assert_eq!(c.name, before.name);
    assert_eq!(c.string_data_arrays, before.string_data_arrays);
    assert!(
        MorphologicalFilter::default()
            .filter_chromatogram(&mut c)
            .is_err()
    );
}

#[test]
fn invalid_parameters_and_overflow_leave_data_unchanged() {
    for element in [
        Element::DataPoints(0),
        Element::Thomson(0.0),
        Element::Thomson(-1.0),
        Element::Seconds(f64::NAN),
    ] {
        assert!(MorphologicalFilter::new(Method::TopHat, element).is_err());
    }
    let mut unsorted = spectrum(&[1.0, 2.0, 3.0], 1.0);
    unsorted.peaks.swap(0, 1);
    let original = unsorted.clone();
    assert!(
        MorphologicalFilter::default()
            .filter_spectrum(&mut unsorted)
            .is_err()
    );
    assert_eq!(unsorted, original);
    let mut degenerate = spectrum(&[1.0, 2.0], 0.0);
    let original = degenerate.clone();
    assert!(
        MorphologicalFilter::default()
            .filter_spectrum(&mut degenerate)
            .is_err()
    );
    assert_eq!(degenerate, original);
    let filter = MorphologicalFilter::new(Method::Gradient, Element::DataPoints(3)).unwrap();
    let mut s = spectrum(&[-f32::MAX, f32::MAX], 1.0);
    let original = s.clone();
    assert!(filter.filter_spectrum(&mut s).is_err());
    assert_eq!(s, original);
    assert!(filter.filter_range(&[f32::NAN]).is_err());
    assert!(MorphologicalFilter::default().filter_range(&[1.0]).is_err());
}

#[test]
fn experiment_filter_is_atomic_and_leaves_chromatograms_unchanged() {
    let filter = MorphologicalFilter::new(Method::TopHat, Element::DataPoints(3)).unwrap();
    let c = MSChromatogram::from_peaks(vec![ChromatogramPeak::new(0.0, 7.0)]);
    let mut e = MSExperiment {
        spectra: vec![spectrum(&[1.0, 4.0, 1.0], 1.0)],
        chromatograms: vec![c.clone()],
        ..Default::default()
    };
    filter.filter_experiment(&mut e).unwrap();
    assert_eq!(values(&e.spectra[0]), [0.0, 3.0, 0.0]);
    assert_eq!(e.chromatograms[0], c);
    let mut bad = spectrum(&[1.0, 2.0], 1.0);
    bad.peaks.swap(0, 1);
    e.spectra.push(bad);
    let original = e.clone();
    assert!(filter.filter_experiment(&mut e).is_err());
    assert_eq!(e, original);
}
