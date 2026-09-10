// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Includes source regressions from pinned Deisotoper_test.cpp (#10067 and
// unequal-charge precursor neutral-mass constraints).

use openms::chemistry::{C13C12_MASSDIFF_U, PROTON_MASS_U};
use openms::comparison::Tolerance;
use openms::kernel::DataArray;
use openms::processing::{SpectrumFilter, deisotoping::Deisotoper};
use openms::{MSExperiment, MSSpectrum, Peak1D, Precursor};

fn ladder(mz: f64, charge: u8, intensities: &[f32]) -> MSSpectrum {
    MSSpectrum::from_peaks(
        intensities
            .iter()
            .enumerate()
            .map(|(i, &intensity)| {
                Peak1D::new(
                    mz + i as f64 * C13C12_MASSDIFF_U / f64::from(charge),
                    intensity,
                )
            })
            .collect(),
    )
}
fn options() -> Deisotoper {
    Deisotoper {
        make_single_charged: false,
        annotate_charge: true,
        annotate_isotope_peak_count: true,
        annotate_features: true,
        ..Default::default()
    }
}
fn integers<'a>(spectrum: &'a MSSpectrum, name: &str) -> &'a [i32] {
    &spectrum
        .integer_data_arrays
        .iter()
        .find(|a| a.name == name)
        .unwrap()
        .data
}

#[test]
fn source_unknown_precursor_charge_does_not_remove_clusters() {
    let base = ladder(200.0, 2, &[1.0; 3]);
    let algorithm = Deisotoper {
        max_charge: 2,
        min_isotope_peaks: 2,
        keep_only_deisotoped: true,
        ..options()
    };
    let expected = algorithm.deisotope(&base).unwrap().spectrum;
    assert_eq!(expected.len(), 1);
    assert_eq!(integers(&expected, "charge"), [2]);
    for precursor in [Precursor::new(200.0, 0), Precursor::new(2000.0, 2)] {
        let mut input = base.clone();
        input.precursors = vec![precursor];
        let output = algorithm.deisotope(&input).unwrap().spectrum;
        assert_eq!(output.peaks, expected.peaks);
        assert_eq!(output.integer_data_arrays, expected.integer_data_arrays);
    }
    let mut input = base;
    input.precursors = vec![Precursor::new(100.0, 1)];
    assert!(algorithm.deisotope(&input).unwrap().spectrum.is_empty());
}

#[test]
fn source_precursor_mass_uses_atomic_units_for_unequal_charges() {
    for mz in [998.5, 999.5] {
        let mut input = ladder(mz, 1, &[100.0, 50.0, 100.0 / 3.0]);
        input.precursors = vec![Precursor::new(500.0, 2)];
        for tolerance in [Tolerance::Absolute(0.01), Tolerance::Ppm(10.0)] {
            for keep_only_deisotoped in [false, true] {
                let algorithm = Deisotoper {
                    tolerance,
                    max_charge: 1,
                    keep_only_deisotoped,
                    ..options()
                };
                let output = algorithm.deisotope(&input).unwrap().spectrum;
                let expected = if mz == 998.5 {
                    1
                } else if keep_only_deisotoped {
                    0
                } else {
                    3
                };
                assert_eq!(output.len(), expected);
                assert!(
                    integers(&output, "charge")
                        .iter()
                        .all(|&q| q == if mz == 998.5 { 1 } else { 0 })
                );
            }
        }
    }
}

#[test]
fn charge_conversion_sum_and_sort_keep_all_annotations_aligned() {
    let mut input = ladder(200.0, 2, &[10.0, 5.0, 2.0]);
    input.peaks.push(Peak1D::new(300.0, 9.0));
    input.float_data_arrays.push(DataArray {
        name: "error".into(),
        data: vec![1.0, 2.0, 3.0, 4.0],
    });
    input.string_data_arrays.push(DataArray {
        name: "identity".into(),
        data: vec!["mono".into(), "iso1".into(), "iso2".into(), "other".into()],
    });
    let result = Deisotoper {
        make_single_charged: true,
        add_up_intensity: true,
        ..options()
    }
    .deisotope(&input)
    .unwrap();
    assert_eq!(result.clusters[0].peak_indices, [0, 1, 2]);
    let output = result.spectrum;
    assert_eq!(output.len(), 2);
    assert_eq!(output.peaks[0], input.peaks[3]);
    assert!((output.peaks[1].mz - (400.0 - PROTON_MASS_U)).abs() < 1e-12);
    assert_eq!(output.peaks[1].intensity, 17.0);
    assert_eq!(integers(&output, "charge"), [0, 2]);
    assert_eq!(integers(&output, "iso_peak_count"), [1, 3]);
    assert_eq!(integers(&output, "feature_number"), [-1, 0]);
    assert_eq!(output.float_data_arrays[0].data, [4.0, 1.0]);
    assert_eq!(output.string_data_arrays[0].data, ["other", "mono"]);
    assert_eq!(input.len(), 4);
}

#[test]
fn decreasing_model_starts_at_configured_isotope_and_partial_counts_do_not_leak() {
    let input = ladder(500.0, 1, &[1.0, 10.0, 5.0, 6.0]);
    let result = Deisotoper {
        min_charge: 1,
        max_charge: 1,
        ..options()
    }
    .deisotope(&input)
    .unwrap();
    assert_eq!(result.clusters.len(), 1);
    assert_eq!(result.clusters[0].peak_indices, [0, 1, 2]);
    assert_eq!(integers(&result.spectrum, "iso_peak_count"), [3, 1]);
    let all = Deisotoper {
        max_charge: 1,
        use_decreasing_model: false,
        ..options()
    }
    .deisotope(&input)
    .unwrap();
    assert_eq!(all.clusters[0].peak_indices, [0, 1, 2, 3]);
    let two = ladder(500.0, 1, &[10.0, 5.0]);
    let output = options().deisotope(&two).unwrap().spectrum;
    assert_eq!(integers(&output, "iso_peak_count"), [1, 1]);
    assert_eq!(integers(&output, "feature_number"), [-1, -1]);
}

#[test]
fn missing_isotopes_tolerance_endpoints_and_high_charge_priority() {
    let full = ladder(100.0, 2, &[10.0, 9.0, 8.0, 7.0, 6.0]);
    let output = options().deisotope(&full).unwrap();
    assert_eq!(output.clusters[0].charge, 2);
    let mut missing = full.clone();
    missing.peaks.remove(1);
    let output = Deisotoper {
        min_charge: 2,
        max_charge: 2,
        ..options()
    }
    .deisotope(&missing)
    .unwrap();
    assert_eq!(output.clusters[0].peak_indices, [1, 2, 3]);
    let mut boundary = ladder(100.0, 1, &[10.0, 9.0]);
    boundary.peaks[1].mz += 0.1;
    assert_eq!(
        Deisotoper {
            tolerance: Tolerance::Absolute(0.1),
            max_charge: 1,
            min_isotope_peaks: 2,
            ..options()
        }
        .deisotope(&boundary)
        .unwrap()
        .clusters
        .len(),
        1
    );
}

#[test]
fn repeated_or_already_consumed_peaks_cannot_form_new_clusters() {
    // At very high m/z, even supported ppm windows can contain the seed itself.
    // A single observed peak cannot stand in for several isotopologues.
    let input = ladder(1e6, 3, &[10.0]);
    let output = Deisotoper {
        tolerance: Tolerance::Ppm(100.0),
        keep_only_deisotoped: true,
        ..options()
    }
    .deisotope(&input)
    .unwrap();
    assert!(output.clusters.is_empty());
    assert!(output.spectrum.is_empty());
    // Every accepted input peak occurs in exactly one returned cluster.
    let mut mixed = ladder(200.0, 3, &[10.0; 6]);
    mixed.peaks.extend(ladder(201.0, 2, &[8.0; 4]).peaks);
    mixed.sort_by_position().unwrap();
    let output = options().deisotope(&mixed).unwrap();
    let members: Vec<_> = output
        .clusters
        .iter()
        .flat_map(|c| &c.peak_indices)
        .collect();
    let unique: std::collections::BTreeSet<_> = members.iter().copied().collect();
    assert_eq!(members.len(), unique.len());
}

#[test]
fn invalid_options_data_resource_and_annotation_collisions_are_atomic() {
    for tolerance in [
        Tolerance::Ppm(-1.0),
        Tolerance::Ppm(100.01),
        Tolerance::Absolute(0.10001),
        Tolerance::Absolute(f64::NAN),
    ] {
        assert!(!Deisotoper::is_tolerance_supported(tolerance));
        assert!(
            Deisotoper {
                tolerance,
                ..options()
            }
            .deisotope(&MSSpectrum::default())
            .is_err()
        );
    }
    assert!(Deisotoper::is_tolerance_supported(Tolerance::Ppm(100.0)));
    assert!(Deisotoper::is_tolerance_supported(Tolerance::Absolute(0.1)));
    let input = ladder(500.0, 1, &[f32::MAX; 3]);
    for algorithm in [
        Deisotoper {
            min_charge: 0,
            ..options()
        },
        Deisotoper {
            max_work: 1,
            ..options()
        },
        Deisotoper {
            min_isotope_peaks: 1,
            ..options()
        },
        Deisotoper {
            add_up_intensity: true,
            ..options()
        },
    ] {
        let mut copy = input.clone();
        assert!(algorithm.filter_spectrum(&mut copy).is_err());
        assert_eq!(copy, input);
    }
    let mut collision = input.clone();
    collision.integer_data_arrays.push(DataArray {
        name: "charge".into(),
        data: vec![],
    });
    assert!(options().deisotope(&collision).is_err());
    let mut unsorted = input.clone();
    unsorted.peaks.reverse();
    assert!(options().deisotope(&unsorted).is_err());
    let mut negative = input.clone();
    negative.peaks[0].intensity = -1.0;
    assert!(options().deisotope(&negative).is_err());
    let mut experiment = MSExperiment {
        spectra: vec![input, unsorted],
        ..Default::default()
    };
    let before = experiment.clone();
    assert!(options().filter_experiment(&mut experiment).is_err());
    assert_eq!(experiment, before);
}

#[test]
fn empty_input_and_unassigned_zero_peaks_follow_source_simple_implementation() {
    let empty = MSSpectrum::default();
    assert_eq!(options().deisotope(&empty).unwrap().spectrum, empty);
    let input = ladder(500.0, 1, &[0.0]);
    let output = options().deisotope(&input).unwrap().spectrum;
    assert_eq!(output.peaks, input.peaks);
    assert_eq!(integers(&output, "charge"), [0]);
}
