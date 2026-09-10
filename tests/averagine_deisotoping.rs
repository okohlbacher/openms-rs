// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::chemistry::{C13C12_MASSDIFF_U, PROTON_MASS_U};
use openms::comparison::Tolerance;
use openms::kernel::DataArray;
use openms::processing::{SpectrumFilter, deisotoping::AveragineDeisotoper};
use openms::{MSExperiment, MSSpectrum, Peak1D, Precursor};

fn options() -> AveragineDeisotoper {
    AveragineDeisotoper {
        min_charge: 1,
        max_charge: 1,
        min_isotope_peaks: 2,
        max_isotope_peaks: 2,
        top_n: None,
        make_single_charged: false,
        annotate_charge: true,
        annotate_isotope_peak_count: true,
        annotate_features: true,
        allow_shared_isotopes: false,
        ..Default::default()
    }
}
fn spectrum(peaks: &[(f64, f32)]) -> MSSpectrum {
    MSSpectrum::from_peaks(
        peaks
            .iter()
            .map(|&(mz, intensity)| Peak1D::new(mz, intensity))
            .collect(),
    )
}
fn ladder() -> MSSpectrum {
    let mz = 1800. + PROTON_MASS_U;
    spectrum(&[(mz, 120.), (mz + C13C12_MASSDIFF_U, 120.)])
}
fn array<'a>(spectrum: &'a MSSpectrum, name: &str) -> &'a [i32] {
    &spectrum
        .integer_data_arrays
        .iter()
        .find(|a| a.name == name)
        .unwrap()
        .data
}

#[test]
fn defaults_follow_source_and_include_fixed_threshold() {
    let config = AveragineDeisotoper::default();
    assert_eq!(config.top_n, Some(5000));
    assert!(config.allow_shared_isotopes);
    assert_eq!((config.min_charge, config.max_charge), (1, 3));
    assert_eq!(
        (config.min_isotope_peaks, config.max_isotope_peaks),
        (2, 10)
    );
    assert!(!config.keep_only_deisotoped);
    assert!(config.make_single_charged);
    let threshold = 0.05_f32;
    let below = f32::from_bits(threshold.to_bits() - 1);
    let mut input = spectrum(&[(100., 0.), (200., below), (300., threshold), (400., 1.)]);
    input
        .integer_data_arrays
        .push(DataArray::new("origin", vec![0, 1, 2, 3]));
    let output = config.deisotope(&input).unwrap().spectrum;
    assert_eq!(
        output.peaks,
        spectrum(&[(300., threshold), (400., 1.)]).peaks
    );
    assert_eq!(array(&output, "origin"), [2, 3]);
    assert_eq!(input.len(), 4);
}

#[test]
fn top_n_is_preprocessing_and_preserves_array_alignment() {
    let mz = 1800. + PROTON_MASS_U;
    let mut input = spectrum(&[(mz, 120.), (mz + C13C12_MASSDIFF_U, 120.), (3000., 200.)]);
    input
        .integer_data_arrays
        .push(DataArray::new("origin", vec![0, 1, 2]));
    input.string_data_arrays.push(DataArray::new(
        "label",
        vec!["mono".into(), "isotope".into(), "noise".into()],
    ));
    input
        .float_data_arrays
        .push(DataArray::new("signal", vec![1., 2., 3.]));
    input.metadata.insert("preserve".into(), "metadata".into());
    let all = options().deisotope(&input).unwrap();
    assert_eq!(all.clusters.len(), 1);
    let filtered = AveragineDeisotoper {
        top_n: Some(2),
        ..options()
    }
    .deisotope(&input)
    .unwrap();
    assert!(filtered.clusters.is_empty());
    assert_eq!(array(&filtered.spectrum, "origin"), [0, 2]); // equal intensities retain input order
    assert_eq!(
        filtered.spectrum.string_data_arrays[0].data,
        ["mono", "noise"]
    );
    assert_eq!(filtered.spectrum.float_data_arrays[0].data, [1., 3.]);
    assert_eq!(filtered.spectrum.metadata, input.metadata);
    let only = AveragineDeisotoper {
        top_n: Some(2),
        keep_only_deisotoped: true,
        ..options()
    }
    .deisotope(&input)
    .unwrap();
    assert!(only.spectrum.is_empty());
}

#[test]
fn cluster_indices_reference_original_input_after_both_selections() {
    let mz = 1800. + PROTON_MASS_U;
    let mut input = spectrum(&[
        (10., 0.),
        (20., 0.1),
        (mz, 120.),
        (mz + C13C12_MASSDIFF_U, 120.),
        (3000., 50.),
    ]);
    input
        .integer_data_arrays
        .push(DataArray::new("original", vec![0, 1, 2, 3, 4]));
    let result = AveragineDeisotoper {
        top_n: Some(2),
        ..options()
    }
    .deisotope(&input)
    .unwrap();
    assert_eq!(result.clusters[0].peak_indices, [2, 3]);
    assert_eq!(array(&result.spectrum, "original"), [2]);
    assert_eq!(array(&result.spectrum, "charge"), [1]);
    assert_eq!(array(&result.spectrum, "iso_peak_count"), [2]);
    assert_eq!(array(&result.spectrum, "feature_number"), [0]);
}

#[test]
fn accepted_clusters_are_disjoint_and_members_advance() {
    let mz = 1800. + PROTON_MASS_U;
    let input = spectrum(&[
        (mz, 120.),
        (mz + 0.1, 120.),
        (mz + C13C12_MASSDIFF_U + 0.05, 120.),
    ]);
    let result = AveragineDeisotoper {
        tolerance: Tolerance::Absolute(0.1),
        ..options()
    }
    .deisotope(&input)
    .unwrap();
    assert_eq!(result.clusters.len(), 1);
    assert_eq!(result.clusters[0].peak_indices, [0, 2]);
    assert_eq!(array(&result.spectrum, "charge"), [1, 0]);
    assert_eq!(array(&result.spectrum, "iso_peak_count"), [2, 1]);
    assert_eq!(array(&result.spectrum, "feature_number"), [0, -1]);
    // The source's wide-tolerance/high-charge search can find its own seed.
    let single = spectrum(&[(1800. / 255. + PROTON_MASS_U, 120.)]);
    let result = AveragineDeisotoper {
        min_charge: 255,
        max_charge: 255,
        tolerance: Tolerance::Absolute(0.1),
        ..options()
    }
    .deisotope(&single)
    .unwrap();
    assert!(result.clusters.is_empty());
    assert_eq!(result.spectrum.peaks, single.peaks);
}

#[test]
fn precursor_unknown_charge_and_multiple_precursors_do_not_limit_mapping() {
    for precursors in [
        vec![],
        vec![Precursor::new(100., 0)],
        vec![Precursor::new(10_000., 2)],
        vec![Precursor::new(100., 1), Precursor::new(100., 1)],
    ] {
        let mut input = ladder();
        input.precursors = precursors.clone();
        let result = options().deisotope(&input).unwrap();
        assert_eq!(result.clusters.len(), 1);
        assert_eq!(result.spectrum.precursors, precursors);
    }
    let mut input = ladder();
    input.precursors = vec![Precursor::new(100., 1)];
    assert!(options().deisotope(&input).unwrap().clusters.is_empty());
}

#[test]
fn zero_model_probability_rejects_instead_of_accepting_nan() {
    let input = spectrum(&[
        (PROTON_MASS_U, 120.),
        (PROTON_MASS_U + C13C12_MASSDIFF_U, 120.),
    ]);
    let result = options().deisotope(&input).unwrap();
    assert!(result.clusters.is_empty());
    assert_eq!(result.spectrum.len(), 2);
    let negative = spectrum(&[(0.5, 120.), (0.5 + C13C12_MASSDIFF_U, 120.)]);
    assert!(options().deisotope(&negative).unwrap().clusters.is_empty());
    let mut mixed = negative.clone();
    mixed.peaks.extend(ladder().peaks);
    let result = options().deisotope(&mixed).unwrap();
    assert_eq!(result.clusters.len(), 1);
    assert_eq!(result.clusters[0].peak_indices, [2, 3]);
    assert_eq!(result.spectrum.len(), 3);
    // Model-generation log fallback remains checked when only its far tail has
    // representable probability. Underflowed prefix probabilities cannot pass KL.
    let huge = spectrum(&[(1e12, 120.), (1e12 + C13C12_MASSDIFF_U, 120.)]);
    assert!(options().deisotope(&huge).unwrap().clusters.is_empty());
}

#[test]
fn summed_intensity_and_converted_coordinates_keep_annotations_aligned() {
    let mz = 900. + PROTON_MASS_U;
    let mut input = spectrum(&[
        (1000., 50.),
        (mz, 120.),
        (mz + C13C12_MASSDIFF_U / 2., 120.),
    ]);
    input.sort_by_position().unwrap();
    input
        .integer_data_arrays
        .push(DataArray::new("origin", vec![0, 1, 2]));
    let result = AveragineDeisotoper {
        min_charge: 2,
        max_charge: 2,
        make_single_charged: true,
        add_up_intensity: true,
        ..options()
    }
    .deisotope(&input)
    .unwrap();
    assert_eq!(result.clusters[0].peak_indices, [0, 1]);
    assert_eq!(result.spectrum.peaks[0].mz, 1000.);
    assert!((result.spectrum.peaks[1].mz - (1800. + PROTON_MASS_U)).abs() < 1e-12);
    assert_eq!(result.spectrum.peaks[1].intensity, 240.);
    assert_eq!(array(&result.spectrum, "origin"), [2, 0]);
    assert_eq!(array(&result.spectrum, "charge"), [0, 2]);
    assert_eq!(array(&result.spectrum, "iso_peak_count"), [1, 2]);
    assert_eq!(array(&result.spectrum, "feature_number"), [-1, 0]);
}

#[test]
fn resource_budget_counts_kl_terms_and_filter_errors_are_atomic() {
    let config = AveragineDeisotoper {
        max_work: 8,
        ..options()
    };
    let mut input = ladder();
    let original = input.clone();
    // This two-peak run uses eight work units before counting its two KL terms.
    assert!(config.filter_spectrum(&mut input).is_err());
    assert_eq!(input, original);
    let mut first = ladder();
    first.metadata.insert("keep".into(), "first".into());
    let mut second = ladder();
    second.peaks[0].intensity = -1.;
    let mut experiment = MSExperiment {
        spectra: vec![first, second],
        ..Default::default()
    };
    let original = experiment.clone();
    assert!(options().filter_experiment(&mut experiment).is_err());
    assert_eq!(experiment, original);
}

#[test]
fn invalid_options_arrays_and_numeric_overflow_are_rejected() {
    for config in [
        AveragineDeisotoper {
            top_n: Some(0),
            ..options()
        },
        AveragineDeisotoper {
            min_charge: 0,
            ..options()
        },
        AveragineDeisotoper {
            min_isotope_peaks: 1,
            ..options()
        },
        AveragineDeisotoper {
            max_isotope_peaks: 1,
            ..options()
        },
        AveragineDeisotoper {
            max_isotope_peaks: usize::MAX,
            ..options()
        },
        AveragineDeisotoper {
            max_work: 0,
            ..options()
        },
        AveragineDeisotoper {
            tolerance: Tolerance::Absolute(0.1001),
            ..options()
        },
        AveragineDeisotoper {
            tolerance: Tolerance::Ppm(f64::NAN),
            ..options()
        },
    ] {
        assert!(config.deisotope(&ladder()).is_err());
    }
    for name in ["charge", "iso_peak_count", "feature_number"] {
        let mut input = ladder();
        input
            .string_data_arrays
            .push(DataArray::new(name, vec!["a".into(), "b".into()]));
        let original = input.clone();
        assert!(options().filter_spectrum(&mut input).is_err());
        assert_eq!(input, original);
    }
    let mut malformed = ladder();
    malformed
        .float_data_arrays
        .push(DataArray::new("bad", vec![1.]));
    assert!(options().deisotope(&malformed).is_err());
    let mut unsorted = ladder();
    unsorted.peaks.reverse();
    assert!(options().deisotope(&unsorted).is_err());
    let mut overflow = ladder();
    for peak in &mut overflow.peaks {
        peak.intensity = f32::MAX;
    }
    let original = overflow.clone();
    assert!(
        AveragineDeisotoper {
            add_up_intensity: true,
            ..options()
        }
        .filter_spectrum(&mut overflow)
        .is_err()
    );
    assert_eq!(overflow, original);
    let overflow = spectrum(&[(f64::MAX, 1.)]);
    assert!(
        AveragineDeisotoper {
            min_charge: 2,
            max_charge: 2,
            ..options()
        }
        .deisotope(&overflow)
        .is_err()
    );
}

#[test]
fn empty_and_entirely_thresholded_spectra_remain_well_formed() {
    let empty = MSSpectrum::new();
    assert_eq!(options().deisotope(&empty).unwrap().spectrum, empty);
    let mut input = spectrum(&[(1., 0.), (2., 0.04)]);
    input
        .string_data_arrays
        .push(DataArray::new("labels", vec!["a".into(), "b".into()]));
    let result = options().deisotope(&input).unwrap();
    assert!(result.spectrum.is_empty());
    assert!(array(&result.spectrum, "charge").is_empty());
    assert!(array(&result.spectrum, "iso_peak_count").is_empty());
    result.spectrum.validate().unwrap();
}

#[test]
fn shared_isotope_policy_is_available_on_both_algorithms() {
    use openms::processing::deisotoping::Deisotoper;
    let mz = 1800. + PROTON_MASS_U;
    let input = spectrum(&[
        (mz, 120.),
        (mz + 0.1, 120.),
        (mz + C13C12_MASSDIFF_U + 0.05, 120.),
    ]);
    assert!(!Deisotoper::default().allow_shared_isotopes);
    for shared in [false, true] {
        let simple = Deisotoper {
            tolerance: Tolerance::Absolute(0.1),
            min_charge: 1,
            max_charge: 1,
            min_isotope_peaks: 2,
            max_isotope_peaks: 2,
            make_single_charged: false,
            allow_shared_isotopes: shared,
            add_up_intensity: true,
            keep_only_deisotoped: true,
            ..Default::default()
        }
        .deisotope(&input)
        .unwrap();
        let averaged = AveragineDeisotoper {
            tolerance: Tolerance::Absolute(0.1),
            allow_shared_isotopes: shared,
            add_up_intensity: true,
            keep_only_deisotoped: true,
            ..options()
        }
        .deisotope(&input)
        .unwrap();
        assert_eq!(simple.clusters, averaged.clusters);
        assert_eq!(simple.spectrum.peaks, averaged.spectrum.peaks);
        assert_eq!(simple.clusters.len(), if shared { 2 } else { 1 });
        if shared {
            assert_eq!(simple.clusters[0].peak_indices, [0, 2]);
            assert_eq!(simple.clusters[1].peak_indices, [1, 2]);
            // C++ counts a shared isotope in the summed signal of each cluster.
            assert_eq!(
                simple
                    .spectrum
                    .peaks
                    .iter()
                    .map(|p| f64::from(p.intensity))
                    .sum::<f64>(),
                480.
            );
            assert_eq!(
                input
                    .peaks
                    .iter()
                    .map(|p| f64::from(p.intensity))
                    .sum::<f64>(),
                360.
            );
        }
    }
}
