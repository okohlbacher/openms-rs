// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::kernel::{DataArray, NumericRange, SpectrumType};
use openms::{
    ChromatogramPeak, Error, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, Precursor,
};

// Ported from OpenMS MSSpectrum_test.cpp, spec_test (revision 7c029e8).
fn reference_spectrum() -> MSSpectrum {
    MSSpectrum::from_peaks(vec![
        Peak1D::new(412.321, 29.0),
        Peak1D::new(412.824, 60.0),
        Peak1D::new(413.8, 34.0),
        Peak1D::new(414.301, 29.0),
        Peak1D::new(415.287, 37.0),
        Peak1D::new(416.293, 31.0),
        Peak1D::new(418.232, 31.0),
        Peak1D::new(419.113, 31.0),
        Peak1D::new(420.13, 201.0),
        Peak1D::new(423.269, 56.0),
        Peak1D::new(426.292, 34.0),
        Peak1D::new(427.28, 82.0),
        Peak1D::new(428.322, 87.0),
        Peak1D::new(430.269, 30.0),
        Peak1D::new(431.246, 29.0),
        Peak1D::new(432.289, 42.0),
        Peak1D::new(436.161, 32.0),
        Peak1D::new(437.219, 54.0),
        Peak1D::new(439.186, 40.0),
        Peak1D::new(440.27, 40.0),
        Peak1D::new(441.224, 23.0),
    ])
}

fn annotated_spectrum() -> MSSpectrum {
    MSSpectrum {
        peaks: vec![
            Peak1D::new(30.0, 4.0),
            Peak1D::new(10.0, 2.0),
            Peak1D::new(20.0, 4.0),
        ],
        rt: 42.0,
        name: "example".into(),
        float_data_arrays: vec![
            DataArray::new("signal to noise", vec![3.0, 1.0, 2.0]),
            DataArray::new("empty", vec![]),
        ],
        integer_data_arrays: vec![DataArray::new("charge", vec![3, 1, 2])],
        string_data_arrays: vec![DataArray::new(
            "labels",
            vec!["c".into(), "a".into(), "b".into()],
        )],
        ..MSSpectrum::default()
    }
}

#[test]
fn defaults_match_openms() {
    assert_eq!(Peak1D::default(), Peak1D::new(0.0, 0.0));
    assert_eq!(ChromatogramPeak::default(), ChromatogramPeak::new(0.0, 0.0));
    assert_eq!(
        Precursor::new(500.0, 2),
        Precursor {
            mz: 500.0,
            intensity: 0.0,
            charge: 2,
            ..Precursor::default()
        }
    );
    let spectrum = MSSpectrum::new();
    assert_eq!(spectrum.rt, -1.0);
    assert_eq!(spectrum.ms_level, 1);
    assert_eq!(spectrum.spectrum_type, SpectrumType::Unknown);
    assert!(spectrum.validate().is_ok());
    assert_eq!(spectrum.clone(), spectrum);
}

#[test]
fn upstream_nearest_peak_reference_cases() {
    let spectrum = reference_spectrum();
    for (query, index) in [
        (400.0, 0),
        (500.0, 20),
        (412.4, 0),
        (441.224, 20),
        (426.29, 10),
        (426.3, 10),
        (427.2, 11),
        (427.3, 11),
    ] {
        assert_eq!(spectrum.find_nearest(query).unwrap(), Some(index));
    }
    for (query, tolerance, index) in [
        (400.0, 1.0, None),
        (500.0, 1.0, None),
        (412.4, 0.01, None),
        (412.4, 0.1, Some(0)),
        (441.3, 0.01, None),
        (441.3, 0.1, Some(20)),
        (427.3, 0.001, None),
    ] {
        assert_eq!(
            spectrum
                .find_nearest_with_tolerance(query, tolerance)
                .unwrap(),
            index
        );
    }
    assert_eq!(
        spectrum.find_nearest_in_window(427.3, 0.1, 0.001).unwrap(),
        Some(11)
    );
    assert_eq!(
        spectrum.find_nearest_in_window(427.3, 0.001, 1.01).unwrap(),
        None
    );
    assert_eq!(
        spectrum.find_nearest_in_window(427.3, 0.001, 1.1).unwrap(),
        Some(12)
    );
}

#[test]
fn binary_bounds_include_duplicates_and_exact_endpoints() {
    let spectrum = MSSpectrum::from_peaks(vec![
        Peak1D::new(1.0, 1.0),
        Peak1D::new(2.0, 2.0),
        Peak1D::new(2.0, 3.0),
        Peak1D::new(4.0, 4.0),
    ]);
    assert_eq!(spectrum.mz_begin(-1.0).unwrap(), 0);
    assert_eq!(spectrum.mz_begin(2.0).unwrap(), 1);
    assert_eq!(spectrum.mz_end(2.0).unwrap(), 3);
    assert_eq!(spectrum.mz_begin(3.0).unwrap(), 3);
    assert_eq!(spectrum.mz_end(4.0).unwrap(), 4);
    assert_eq!(spectrum.mz_begin(8.0).unwrap(), 4);
    assert_eq!(spectrum.find_nearest(2.0).unwrap(), Some(1));
    assert_eq!(spectrum.find_nearest(3.0).unwrap(), Some(2));
    assert_eq!(
        spectrum.find_nearest_with_tolerance(3.0, 1.0).unwrap(),
        Some(2)
    );
    assert_eq!(
        spectrum.find_nearest_in_window(3.0, 0.1, 1.0).unwrap(),
        Some(3)
    );
    assert_eq!(
        spectrum.find_nearest_in_window(3.5, 1.5, 0.1).unwrap(),
        Some(2)
    );
}

#[test]
fn empty_peak_searches_are_safe() {
    let spectrum = MSSpectrum::new();
    assert_eq!(spectrum.mz_begin(50.0).unwrap(), 0);
    assert_eq!(spectrum.mz_end(50.0).unwrap(), 0);
    assert_eq!(spectrum.find_nearest(50.0).unwrap(), None);
    assert_eq!(
        spectrum.find_nearest_in_window(50.0, 1.0, 1.0).unwrap(),
        None
    );
    assert_eq!(
        spectrum.find_highest_in_window(50.0, 1.0, 1.0).unwrap(),
        None
    );
    assert_eq!(spectrum.base_peak(), None);
    assert_eq!(spectrum.calculate_tic(), 0.0);
    assert_eq!(spectrum.ranges().unwrap().mz, None);
}

#[test]
fn searches_reject_unsorted_nonfinite_and_negative_tolerances() {
    let unsorted = annotated_spectrum();
    assert!(matches!(unsorted.mz_begin(20.0), Err(Error::UnsortedData)));
    assert!(matches!(
        unsorted.find_nearest(20.0),
        Err(Error::UnsortedData)
    ));
    let sorted = reference_spectrum();
    for query in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            sorted.mz_begin(query),
            Err(Error::InvalidValue(_))
        ));
        assert!(matches!(
            sorted.find_nearest(query),
            Err(Error::InvalidValue(_))
        ));
    }
    assert!(sorted.find_nearest_with_tolerance(420.0, -1.0).is_err());
    assert!(sorted.find_nearest_with_tolerance(420.0, f64::NAN).is_err());
    assert!(
        sorted
            .find_nearest_in_window(420.0, 1.0, f64::INFINITY)
            .is_err()
    );
    for coordinate in [f64::NAN, f64::INFINITY] {
        let invalid = MSSpectrum::from_peaks(vec![Peak1D::new(coordinate, 1.0)]);
        assert!(!invalid.is_sorted());
        assert!(invalid.find_nearest(1.0).is_err());
    }
}

#[test]
fn upstream_base_peak_and_tic_values() {
    let spectrum = reference_spectrum();
    assert_eq!(spectrum.calculate_tic(), 1032.0);
    assert_eq!(spectrum.base_peak(), Some(&spectrum.peaks[8]));
    let ties = MSSpectrum::from_peaks(vec![
        Peak1D::new(1.0, -2.0),
        Peak1D::new(2.0, -1.0),
        Peak1D::new(3.0, -1.0),
    ]);
    assert_eq!(ties.base_peak(), Some(&ties.peaks[1]));
    assert_eq!(ties.calculate_tic(), -4.0);
}

#[test]
fn highest_in_window_selects_first_maximum_in_inclusive_window() {
    let spectrum = MSSpectrum::from_peaks(vec![
        Peak1D::new(1.0, 10.0),
        Peak1D::new(2.0, 20.0),
        Peak1D::new(2.0, 20.0),
        Peak1D::new(3.0, 50.0),
    ]);
    assert_eq!(
        spectrum.find_highest_in_window(2.0, 1.0, 0.0).unwrap(),
        Some(1)
    );
    assert_eq!(
        spectrum.find_highest_in_window(2.0, 0.0, 1.0).unwrap(),
        Some(3)
    );
    assert_eq!(
        spectrum.find_highest_in_window(10.0, 0.5, 0.5).unwrap(),
        None
    );
    assert!(spectrum.find_highest_in_window(2.0, -1.0, 1.0).is_err());
}

#[test]
fn sorting_permutes_all_arrays_stably() {
    let mut spectrum = annotated_spectrum();
    spectrum.sort_by_position().unwrap();
    assert_eq!(
        spectrum
            .peaks
            .iter()
            .map(|peak| peak.mz)
            .collect::<Vec<_>>(),
        vec![10.0, 20.0, 30.0]
    );
    assert_eq!(spectrum.float_data_arrays[0].data, vec![1.0, 2.0, 3.0]);
    assert_eq!(spectrum.integer_data_arrays[0].data, vec![1, 2, 3]);
    assert_eq!(spectrum.string_data_arrays[0].data, vec!["a", "b", "c"]);
    assert!(spectrum.float_data_arrays[1].data.is_empty());
    assert_eq!(spectrum.name, "example");
    assert_eq!(spectrum.rt, 42.0);
    spectrum.sort_by_intensity(true).unwrap();
    assert_eq!(spectrum.integer_data_arrays[0].data, vec![2, 3, 1]);
    spectrum.sort_by_intensity(false).unwrap();
    assert_eq!(spectrum.integer_data_arrays[0].data, vec![1, 2, 3]);
    let mut equal_mz = MSSpectrum::from_peaks(vec![
        Peak1D::new(2.0, 1.0),
        Peak1D::new(1.0, 2.0),
        Peak1D::new(2.0, 3.0),
    ]);
    equal_mz.sort_by_position().unwrap();
    assert_eq!(
        equal_mz
            .peaks
            .iter()
            .map(|peak| peak.intensity)
            .collect::<Vec<_>>(),
        vec![2.0, 1.0, 3.0]
    );
}

#[test]
fn invalid_selection_or_alignment_is_transactional() {
    let mut spectrum = annotated_spectrum();
    let before = spectrum.clone();
    assert!(spectrum.select(&[0, 7]).is_err());
    assert_eq!(spectrum, before);
    assert!(spectrum.select(&[1, 1]).is_err());
    assert_eq!(spectrum, before);
    spectrum.integer_data_arrays[0].data.push(99);
    let before = spectrum.clone();
    assert!(spectrum.select(&[2, 0]).is_err());
    assert_eq!(spectrum, before);
    assert!(spectrum.sort_by_position().is_err());
    assert_eq!(spectrum, before);
    assert!(spectrum.sort_by_intensity(true).is_err());
    assert_eq!(spectrum, before);
}

#[test]
fn selection_and_retention_preserve_annotations_and_metadata() {
    let mut spectrum = annotated_spectrum();
    spectrum.select(&[2, 0]).unwrap();
    assert_eq!(
        spectrum
            .peaks
            .iter()
            .map(|peak| peak.mz)
            .collect::<Vec<_>>(),
        vec![20.0, 30.0]
    );
    assert_eq!(spectrum.string_data_arrays[0].data, vec!["b", "c"]);
    assert_eq!(spectrum.float_data_arrays[0].data, vec![2.0, 3.0]);
    spectrum.retain_peaks(|peak| peak.mz > 25.0).unwrap();
    assert_eq!(spectrum.integer_data_arrays[0].data, vec![3]);
    spectrum.select(&[]).unwrap();
    assert!(spectrum.is_empty());
    assert!(spectrum.string_data_arrays[0].data.is_empty());
    assert_eq!(spectrum.string_data_arrays[0].name, "labels");
    assert_eq!(spectrum.name, "example");
    assert!(spectrum.select(&[]).is_ok());
}

#[test]
fn ranges_follow_mutation_without_update_call() {
    let mut spectrum = annotated_spectrum();
    let ranges = spectrum.ranges().unwrap();
    assert_eq!(
        ranges.mz,
        Some(NumericRange {
            min: 10.0,
            max: 30.0
        })
    );
    assert_eq!(ranges.intensity, Some(NumericRange { min: 2.0, max: 4.0 }));
    spectrum.peaks[0].mz = 100.0;
    spectrum.peaks[0].intensity = -10.0;
    assert_eq!(spectrum.ranges().unwrap().mz.unwrap().max, 100.0);
    assert_eq!(spectrum.ranges().unwrap().intensity.unwrap().min, -10.0);
}

#[test]
fn invalid_numeric_data_is_rejected_before_sorting() {
    let mut spectrum =
        MSSpectrum::from_peaks(vec![Peak1D::new(2.0, f32::INFINITY), Peak1D::new(1.0, 1.0)]);
    let before = spectrum.clone();
    assert!(spectrum.validate().is_err());
    assert!(spectrum.sort_by_position().is_err());
    assert_eq!(spectrum, before);
    spectrum.peaks[0].intensity = 1.0;
    spectrum.ms_level = 0;
    assert!(spectrum.validate().is_err());
    spectrum.ms_level = 1;
    spectrum.rt = f64::NAN;
    assert!(spectrum.validate().is_err());
    spectrum.rt = -1.0;
    spectrum.precursors.push(Precursor::new(f64::INFINITY, 2));
    assert!(spectrum.validate().is_err());
}

#[test]
fn clear_option_preserves_or_resets_metadata() {
    let mut spectrum = annotated_spectrum();
    spectrum.precursors.push(Precursor::new(500.0, 2));
    spectrum.clear(false);
    assert!(spectrum.is_empty());
    assert!(spectrum.string_data_arrays.is_empty());
    assert_eq!(spectrum.name, "example");
    assert_eq!(spectrum.precursors.len(), 1);
    spectrum.clear(true);
    assert_eq!(spectrum, MSSpectrum::new());
}

#[test]
fn chromatogram_searches_and_permutations_match_peak_semantics() {
    let mut chromatogram = MSChromatogram {
        peaks: vec![
            ChromatogramPeak::new(3.0, 5.0),
            ChromatogramPeak::new(1.0, 2.0),
            ChromatogramPeak::new(2.0, 5.0),
        ],
        string_data_arrays: vec![DataArray::new(
            "labels",
            vec!["c".into(), "a".into(), "b".into()],
        )],
        ..MSChromatogram::default()
    };
    assert!(matches!(
        chromatogram.find_nearest(2.0),
        Err(Error::UnsortedData)
    ));
    chromatogram.sort_by_position().unwrap();
    assert_eq!(chromatogram.string_data_arrays[0].data, vec!["a", "b", "c"]);
    assert_eq!(chromatogram.rt_begin(2.0).unwrap(), 1);
    assert_eq!(chromatogram.rt_end(2.0).unwrap(), 2);
    assert_eq!(chromatogram.find_nearest(2.5).unwrap(), Some(1));
    assert_eq!(chromatogram.find_nearest(0.0).unwrap(), Some(0));
    assert_eq!(chromatogram.find_nearest(10.0).unwrap(), Some(2));
    assert_eq!(
        chromatogram.ranges().unwrap().rt,
        Some(NumericRange { min: 1.0, max: 3.0 })
    );
    chromatogram.sort_by_intensity(true).unwrap();
    assert_eq!(chromatogram.string_data_arrays[0].data, vec!["b", "c", "a"]);
    let before = chromatogram.clone();
    assert!(chromatogram.select(&[99]).is_err());
    assert_eq!(chromatogram, before);
    chromatogram.select(&[2, 0]).unwrap();
    assert_eq!(chromatogram.string_data_arrays[0].data, vec!["a", "b"]);
    assert_eq!(MSChromatogram::new().find_nearest(1.0).unwrap(), None);
}

// Ported from the no-binning fixture in MSExperiment_test.cpp::calculateTIC.
fn reference_experiment() -> MSExperiment {
    let scans = [
        (0.0, 1, vec![3.0, 5.0]),
        (2.0, 1, vec![2.0]),
        (2.0, 2, vec![0.5]),
        (5.0, 1, vec![2.0, 3.0, 4.0]),
    ];
    MSExperiment {
        spectra: scans
            .into_iter()
            .map(|(rt, ms_level, intensities)| MSSpectrum {
                peaks: intensities
                    .into_iter()
                    .enumerate()
                    .map(|(i, intensity)| Peak1D::new((i + 1) as f64 * 5.0, intensity))
                    .collect(),
                rt,
                ms_level,
                ..MSSpectrum::default()
            })
            .collect(),
        ..MSExperiment::default()
    }
}

#[test]
fn upstream_experiment_tic_filters_levels_and_preserves_scan_order() {
    let experiment = reference_experiment();
    assert_eq!(
        experiment.calculate_tic(1).peaks,
        vec![
            ChromatogramPeak::new(0.0, 8.0),
            ChromatogramPeak::new(2.0, 2.0),
            ChromatogramPeak::new(5.0, 9.0)
        ]
    );
    assert_eq!(
        experiment.calculate_tic(2).peaks,
        vec![ChromatogramPeak::new(2.0, 0.5)]
    );
    assert_eq!(experiment.calculate_tic(0).len(), 4);
    assert!(experiment.calculate_tic(3).is_empty());
    assert!(MSExperiment::new().calculate_tic(1).is_empty());
    assert_eq!(experiment.ms_levels(), vec![1, 2]);
}

#[test]
fn experiment_rt_searches_use_higher_rt_for_midpoint_ties() {
    let experiment = reference_experiment();
    assert_eq!(experiment.rt_begin(2.0).unwrap(), 1);
    assert_eq!(experiment.rt_end(2.0).unwrap(), 3);
    assert_eq!(experiment.closest_spectrum_in_rt(1.0, 0).unwrap(), Some(1));
    assert_eq!(experiment.closest_spectrum_in_rt(2.0, 1).unwrap(), Some(1));
    assert_eq!(experiment.closest_spectrum_in_rt(2.0, 2).unwrap(), Some(2));
    assert_eq!(experiment.closest_spectrum_in_rt(3.5, 1).unwrap(), Some(3));
    assert_eq!(experiment.closest_spectrum_in_rt(9.0, 2).unwrap(), Some(2));
    assert_eq!(experiment.closest_spectrum_in_rt(-5.0, 1).unwrap(), Some(0));
    assert_eq!(experiment.closest_spectrum_in_rt(9.0, 3).unwrap(), None);
    assert_eq!(
        MSExperiment::new().closest_spectrum_in_rt(2.0, 0).unwrap(),
        None
    );
    assert_eq!(
        MSExperiment::new().closest_spectrum_in_rt(2.0, 2).unwrap(),
        None
    );
}

#[test]
fn experiment_rt_selection_is_inclusive_and_level_aware() {
    let experiment = reference_experiment();
    assert_eq!(
        experiment.spectra_in_rt_range(2.0, 2.0, 0).unwrap().len(),
        2
    );
    assert_eq!(
        experiment.spectra_in_rt_range(2.0, 5.0, 1).unwrap().len(),
        2
    );
    assert!(
        experiment
            .spectra_in_rt_range(3.0, 4.0, 0)
            .unwrap()
            .is_empty()
    );
    assert!(experiment.spectra_in_rt_range(3.0, 2.0, 0).is_err());
    assert!(experiment.spectra_in_rt_range(f64::NAN, 2.0, 0).is_err());
}

#[test]
fn experiment_sorting_and_ranges_are_current_and_transactional() {
    let mut experiment = reference_experiment();
    experiment.spectra.reverse();
    experiment.spectra[0].peaks.reverse();
    assert!(!experiment.is_sorted(false));
    assert!(matches!(experiment.rt_begin(1.0), Err(Error::UnsortedData)));
    experiment.sort_spectra(false).unwrap();
    assert!(experiment.is_sorted(false));
    assert!(!experiment.is_sorted(true));
    experiment.sort_spectra(true).unwrap();
    assert!(experiment.is_sorted(true));
    let ranges = experiment.ranges(1).unwrap();
    assert_eq!(ranges.rt, Some(NumericRange { min: 0.0, max: 5.0 }));
    assert_eq!(
        ranges.mz,
        Some(NumericRange {
            min: 5.0,
            max: 15.0
        })
    );
    assert_eq!(ranges.intensity, Some(NumericRange { min: 2.0, max: 5.0 }));
    assert_eq!(experiment.ranges(3).unwrap().rt, None);
    experiment.spectra[0].peaks[0].mz = 999.0;
    assert_eq!(experiment.ranges(1).unwrap().mz.unwrap().max, 999.0);
    experiment.spectra[1]
        .integer_data_arrays
        .push(DataArray::new("malformed", vec![1, 2, 3]));
    let before = experiment.clone();
    assert!(experiment.sort_spectra(true).is_err());
    assert_eq!(experiment, before);
}

#[test]
fn closest_peak_matches_linear_reference_across_many_queries() {
    let positions = [0.0, 1.5, 2.0, 5.5, 10.0, 11.0, 45.0];
    let spectrum = MSSpectrum::from_peaks(
        positions
            .into_iter()
            .map(|mz| Peak1D::new(mz, 1.0))
            .collect(),
    );
    for i in -20..500 {
        let query = f64::from(i) / 10.0;
        let expected = positions
            .iter()
            .enumerate()
            .min_by(|a, b| {
                (a.1 - query)
                    .abs()
                    .partial_cmp(&(b.1 - query).abs())
                    .unwrap()
            })
            .map(|(i, _)| i);
        assert_eq!(
            spectrum.find_nearest(query).unwrap(),
            expected,
            "query {query}"
        );
    }
}

// The peak checks of `validate` accumulate over the whole peak list without an
// early exit and rescan only when that accumulation failed. The next three
// tests pin what the rescan must reproduce: the source's message for each
// field, its m/z-before-intensity order within a peak, and its first-peak-first
// order across peaks. They also cover a list long enough to run the vectorised
// body and its scalar tail rather than the tail alone.
#[test]
fn nonfinite_peak_values_name_the_first_offending_field() {
    let message = |spectrum: &MSSpectrum| spectrum.validate().unwrap_err().to_string();
    let bad_mz =
        MSSpectrum::from_peaks(vec![Peak1D::new(1.0, 1.0), Peak1D::new(f64::NAN, f32::NAN)]);
    let reported = message(&bad_mz);
    assert!(reported.contains("peak m/z must be finite"), "{reported}");
    // An earlier peak's intensity outranks a later peak's m/z.
    let bad_intensity = MSSpectrum::from_peaks(vec![
        Peak1D::new(1.0, f32::INFINITY),
        Peak1D::new(f64::NEG_INFINITY, 1.0),
    ]);
    let reported = message(&bad_intensity);
    assert!(
        reported.contains("peak intensity must be finite"),
        "{reported}"
    );
    for mz in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            MSSpectrum::from_peaks(vec![Peak1D::new(mz, 1.0)])
                .validate()
                .is_err()
        );
    }
    for intensity in [f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
        assert!(
            MSSpectrum::from_peaks(vec![Peak1D::new(1.0, intensity)])
                .validate()
                .is_err()
        );
    }
}

#[test]
fn a_long_peak_list_is_checked_to_its_last_peak() {
    let peaks = |n: usize| {
        (0..n)
            .map(|i| Peak1D::new(i as f64, 1.0))
            .collect::<Vec<_>>()
    };
    // Lengths either side of any unrolled or vectorised block width, so a
    // nonfinite value in the tail cannot be skipped.
    for length in [1_usize, 2, 3, 7, 8, 9, 15, 16, 17, 31, 33, 1_000] {
        let good = MSSpectrum::from_peaks(peaks(length));
        assert!(good.validate().is_ok(), "length {length}");
        for position in [0, length / 2, length - 1] {
            let mut bad = MSSpectrum::from_peaks(peaks(length));
            bad.peaks[position].intensity = f32::NAN;
            let reported = bad.validate().unwrap_err().to_string();
            assert!(
                reported.contains("peak intensity must be finite"),
                "length {length} position {position}: {reported}"
            );
            let mut bad = MSSpectrum::from_peaks(peaks(length));
            bad.peaks[position].mz = f64::INFINITY;
            assert!(
                bad.validate().is_err(),
                "length {length} position {position}"
            );
        }
    }
}

#[test]
fn nonfinite_chromatogram_values_keep_their_own_messages() {
    let bad_rt = MSChromatogram::from_peaks(vec![ChromatogramPeak::new(f64::NAN, 1.0)]);
    let reported = bad_rt.validate().unwrap_err().to_string();
    assert!(
        reported.contains("chromatogram retention time must be finite"),
        "{reported}"
    );
    let bad_intensity =
        MSChromatogram::from_peaks(vec![ChromatogramPeak::new(1.0, f32::NEG_INFINITY)]);
    let reported = bad_intensity.validate().unwrap_err().to_string();
    assert!(
        reported.contains("chromatogram intensity must be finite"),
        "{reported}"
    );
}

// Finiteness and order are accumulated in one pass over the coordinates, so
// this pins that a nonfinite coordinate is still reported ahead of an unsorted
// pair, and that a container which is merely unsorted still reports that.
#[test]
fn nonfinite_coordinates_outrank_unsorted_coordinates() {
    let nonfinite_and_unsorted = MSSpectrum::from_peaks(vec![
        Peak1D::new(3.0, 1.0),
        Peak1D::new(f64::NAN, 1.0),
        Peak1D::new(1.0, 1.0),
    ]);
    assert!(matches!(
        nonfinite_and_unsorted.mz_begin(2.0),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        nonfinite_and_unsorted.mz_end(2.0),
        Err(Error::InvalidValue(_))
    ));
    assert!(!nonfinite_and_unsorted.is_sorted());
    let unsorted = MSSpectrum::from_peaks(vec![Peak1D::new(3.0, 1.0), Peak1D::new(1.0, 1.0)]);
    assert!(matches!(unsorted.mz_begin(2.0), Err(Error::UnsortedData)));
    assert!(matches!(
        unsorted.find_highest_in_window(2.0, 1.0, 1.0),
        Err(Error::UnsortedData)
    ));
    // An infinite coordinate is nonfinite, not merely out of order, even when
    // it is the first one and so is never the right-hand side of a comparison.
    let infinite_first = MSSpectrum::from_peaks(vec![
        Peak1D::new(f64::NEG_INFINITY, 1.0),
        Peak1D::new(1.0, 1.0),
    ]);
    assert!(matches!(
        infinite_first.mz_begin(2.0),
        Err(Error::InvalidValue(_))
    ));
    // Sorted, finite and empty containers are unaffected.
    assert!(MSSpectrum::new().mz_begin(2.0).is_ok());
    assert!(reference_spectrum().mz_begin(420.0).is_ok());
    let chromatogram = MSChromatogram::from_peaks(vec![
        ChromatogramPeak::new(2.0, 1.0),
        ChromatogramPeak::new(f64::INFINITY, 1.0),
    ]);
    assert!(matches!(
        chromatogram.rt_begin(1.0),
        Err(Error::InvalidValue(_))
    ));
}

// `sort_spectra` validates every spectrum once and then sorts with the checked
// permutation path; the annotation arrays must still travel with the peaks.
#[test]
fn experiment_sorting_keeps_annotation_arrays_aligned() {
    let mut experiment = MSExperiment::new();
    experiment.spectra.push(annotated_spectrum());
    experiment.spectra[0].rt = 2.0;
    let mut second = annotated_spectrum();
    second.rt = 1.0;
    experiment.spectra.push(second);
    experiment.sort_spectra(true).unwrap();
    assert_eq!(experiment.spectra[0].rt, 1.0);
    for spectrum in &experiment.spectra {
        assert_eq!(
            spectrum
                .peaks
                .iter()
                .map(|peak| peak.mz)
                .collect::<Vec<_>>(),
            vec![10.0, 20.0, 30.0]
        );
        assert_eq!(spectrum.float_data_arrays[0].data, vec![1.0, 2.0, 3.0]);
        assert_eq!(spectrum.integer_data_arrays[0].data, vec![1, 2, 3]);
        assert_eq!(spectrum.string_data_arrays[0].data, vec!["a", "b", "c"]);
        assert!(spectrum.float_data_arrays[1].data.is_empty());
    }
}
