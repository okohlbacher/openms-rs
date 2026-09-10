// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Independent source-derived numeric oracles; no C++ execution.
// OpenMS4-core 7c029e8cdba6abab503708ecdd56f6ab55e38ce4:
// src/openms/source/PROCESSING/DEISOTOPING/Deisotoper.cpp SHA-256
// 026cb13ddf15362c5ddf7a7f736e6067b057cee6701f598de93154dfdb96e2e2
// src/openms/source/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.cpp
// 0974fabc36a93f741c210d89e940215fa773306b99b5b45309a35cb4b43a2f8c

use openms::chemistry::{C13C12_MASSDIFF_U, PROTON_MASS_U};
use openms::comparison::Tolerance;
use openms::kernel::{MSSpectrum, Peak1D, Precursor};
use openms::processing::deisotoping::{AveragineDeisotoper, Deisotoper};

fn settings(size: usize) -> AveragineDeisotoper {
    AveragineDeisotoper {
        tolerance: Tolerance::Ppm(10.0),
        min_charge: 1,
        max_charge: 1,
        min_isotope_peaks: size,
        max_isotope_peaks: size,
        keep_only_deisotoped: true,
        make_single_charged: false,
        top_n: None,
        ..Default::default()
    }
}
fn ladder(seed: f64, charge: u8, intensities: &[f32]) -> MSSpectrum {
    MSSpectrum::from_peaks(
        intensities
            .iter()
            .enumerate()
            .map(|(i, &intensity)| {
                Peak1D::new(
                    seed + i as f64 * C13C12_MASSDIFF_U / f64::from(charge),
                    intensity,
                )
            })
            .collect(),
    )
}

#[test]
fn source_mixed_precision_kl_boundaries_for_two_through_seven_peaks() {
    // Neutral mass 1800 gives lambda = 1. Earlier observed prefixes are 120/k!,
    // hence all pass. Vary only the final f32 intensity across adjacent values.
    // Expected bits were derived independently using Python math.log/f64 and
    // struct.pack('f') after EACH accumulator update, matching C++ float += double.
    // For size 4 the accepted value gives KL = 0.20000000298023224. Casting each
    // term to f32 first instead gives 0.20000001788139343 and wrongly rejects it.
    let preceding = [120.0, 120.0, 60.0, 20.0, 5.0, 1.0];
    let cases = [
        (2, 0x4365_a1b1_u32, 0x4365_a1b2_u32),
        (3, 0x431c_d519, 0x431c_d51a),
        (4, 0x42d5_b9ff, 0x42d5_ba00),
        (5, 0x42b3_38eb, 0x42b3_38ec),
        (6, 0x4294_de7a, 0x4294_de7b),
        (7, 0x4243_0c42, 0x4243_0c43),
    ];
    for (size, accepted, rejected) in cases {
        for (bits, expected_clusters) in [(accepted, 1), (rejected, 0)] {
            let mut intensities = preceding[..size - 1].to_vec();
            intensities.push(f32::from_bits(bits));
            let input = ladder(1800.0 + PROTON_MASS_U, 1, &intensities);
            let result = settings(size).deisotope(&input).unwrap();
            assert_eq!(
                result.clusters.len(),
                expected_clusters,
                "size={size}, final intensity bits={bits:08x}"
            );
            assert_eq!(result.spectrum.len(), expected_clusters);
            if let Some(cluster) = result.clusters.first() {
                assert_eq!(cluster.peak_indices, (0..size).collect::<Vec<_>>());
            }
        }
    }
}

fn competing_charge_patterns(include_fourth_q1_peak: bool) -> MSSpectrum {
    // For charge 1, lambda 1 fits offsets 0,1,2,3 with intensities 120,120,60,20.
    // For charge 2, lambda 2 accepts offsets 0,0.5,1 with intensities 120,240,120;
    // its 3-point KL is about 0.04986 (<0.1), then offset 1.5 is missing.
    let mut values = vec![(0.0, 120.0), (0.5, 240.0), (1.0, 120.0), (2.0, 60.0)];
    if include_fourth_q1_peak {
        values.push((3.0, 20.0));
    }
    let seed = 1800.0 + PROTON_MASS_U;
    MSSpectrum::from_peaks(
        values
            .into_iter()
            .map(|(offset, intensity)| Peak1D::new(seed + offset * C13C12_MASSDIFF_U, intensity))
            .collect(),
    )
}

#[test]
fn longest_cluster_beats_higher_charge_and_owns_count_and_intensity() {
    let result = AveragineDeisotoper {
        max_charge: 2,
        min_isotope_peaks: 2,
        annotate_isotope_peak_count: true,
        add_up_intensity: true,
        ..settings(4)
    }
    .deisotope(&competing_charge_patterns(true))
    .unwrap();
    assert_eq!(result.clusters.len(), 1);
    assert_eq!(result.clusters[0].charge, 1);
    assert_eq!(result.clusters[0].peak_indices, [0, 2, 3, 4]);
    assert_eq!(result.spectrum.peaks[0].intensity, 320.0);
    // Source writes this count during each charge trial; the native correction
    // must describe the selected four-member charge-1 ladder, not later charge 2 length 3.
    assert_eq!(result.spectrum.integer_data_arrays[0].data, [4]);
}

#[test]
fn equal_cluster_lengths_choose_the_higher_charge() {
    let result = AveragineDeisotoper {
        max_charge: 2,
        min_isotope_peaks: 2,
        annotate_charge: true,
        add_up_intensity: true,
        ..settings(4)
    }
    .deisotope(&competing_charge_patterns(false))
    .unwrap();
    assert_eq!(result.clusters.len(), 1);
    assert_eq!(result.clusters[0].charge, 2);
    assert_eq!(result.clusters[0].peak_indices, [0, 1, 2]);
    assert_eq!(result.spectrum.peaks[0].intensity, 480.0);
    assert_eq!(result.spectrum.integer_data_arrays[0].data, [2]);
}

#[test]
fn seed_ppm_window_is_inclusive_at_both_f64_endpoints() {
    let seed = 1800.0 + PROTON_MASS_U;
    let expected = seed + C13C12_MASSDIFF_U;
    let tolerance = (10.0 / 1e6) * seed;
    for (endpoint, outward) in [(expected - tolerance, -1_i64), (expected + tolerance, 1)] {
        for (observed, accepted) in [
            (endpoint, true),
            (
                f64::from_bits(endpoint.to_bits().checked_add_signed(outward).unwrap()),
                false,
            ),
        ] {
            let input = MSSpectrum::from_peaks(vec![
                Peak1D::new(seed, 120.0),
                Peak1D::new(observed, 120.0),
            ]);
            assert_eq!(
                settings(2).deisotope(&input).unwrap().clusters.len(),
                usize::from(accepted),
                "observed={observed:.16}, endpoint={endpoint:.16}"
            );
        }
    }
}

#[test]
fn nearest_distance_tie_prefers_the_lower_mz_model_candidate() {
    let seed = 1800.0 + PROTON_MASS_U;
    let expected = seed + C13C12_MASSDIFF_U;
    let delta = 1.0 / 4096.0; // exact binary fraction, symmetric at this magnitude
    let input = MSSpectrum::from_peaks(vec![
        Peak1D::new(seed, 120.0),
        Peak1D::new(expected - delta, 120.0),
        Peak1D::new(expected + delta, 1200.0),
    ]);
    let result = AveragineDeisotoper {
        tolerance: Tolerance::Absolute(0.001),
        ..settings(2)
    }
    .deisotope(&input)
    .unwrap();
    assert_eq!(result.clusters.len(), 1);
    assert_eq!(result.clusters[0].peak_indices, [0, 1]);
}

#[test]
fn precursor_constraint_preserves_source_operation_order_at_zero_tolerance() {
    // For seed 1234.56789 and charge 3 the source constraint is 3700.681840599687,
    // whereas the Poisson mass expression q*(seed-PROTON) is 3700.6818405996873.
    // Using that latter expression for rejection would lose the exact-boundary cluster.
    let generator = AveragineDeisotoper {
        tolerance: Tolerance::Absolute(0.0),
        min_charge: 3,
        max_charge: 3,
        ..settings(2)
    };
    let simple = Deisotoper {
        tolerance: generator.tolerance,
        min_charge: 3,
        max_charge: 3,
        min_isotope_peaks: 2,
        max_isotope_peaks: 2,
        keep_only_deisotoped: true,
        make_single_charged: false,
        use_decreasing_model: false,
        ..Default::default()
    };
    for (precursor_bits, accepted) in [
        (0x40ac_eb60_d3f3_bf1a, false),
        (0x40ac_eb60_d3f3_bf1b, true),
        (0x40ac_eb60_d3f3_bf1c, true),
    ] {
        let mut input = ladder(1234.56789, 3, &[120.0, 246.71213]);
        input.precursors = vec![Precursor::new(f64::from_bits(precursor_bits), 1)];
        assert_eq!(
            generator.deisotope(&input).unwrap().clusters.len(),
            usize::from(accepted),
            "precursor bits={precursor_bits:016x}"
        );
        assert_eq!(
            simple.deisotope(&input).unwrap().clusters.len(),
            usize::from(accepted),
            "simple method: precursor bits={precursor_bits:016x}"
        );
    }
}

#[test]
fn shared_isotope_belongs_to_both_clusters_and_is_summed_in_each() {
    // Two source-valid seeds: charge 1 uses positions 0,2, while charge 3 uses 1,2.
    // Lambda 1 fits 120:120 and lambda ~3 fits 40:120. The second seed is
    // unassigned when visited; only its extension was claimed previously.
    let seed = 1800.0 + PROTON_MASS_U;
    let input = MSSpectrum::from_peaks(vec![
        Peak1D::new(seed, 120.0),
        Peak1D::new(seed + 2.0 * C13C12_MASSDIFF_U / 3.0, 40.0),
        Peak1D::new(seed + C13C12_MASSDIFF_U, 120.0),
    ]);
    let config = AveragineDeisotoper {
        max_charge: 3,
        annotate_charge: true,
        annotate_isotope_peak_count: true,
        annotate_features: true,
        add_up_intensity: true,
        ..settings(2)
    };
    assert!(config.allow_shared_isotopes);
    let result = config.deisotope(&input).unwrap();
    assert_eq!(result.clusters.len(), 2);
    assert_eq!(result.clusters[0].charge, 1);
    assert_eq!(result.clusters[0].peak_indices, [0, 2]);
    assert_eq!(result.clusters[1].charge, 3);
    assert_eq!(result.clusters[1].peak_indices, [1, 2]);
    assert_eq!(result.spectrum.peaks.len(), 2);
    assert_eq!(result.spectrum.peaks[0].mz, input.peaks[0].mz);
    assert_eq!(result.spectrum.peaks[1].mz, input.peaks[1].mz);
    assert_eq!(result.spectrum.peaks[0].intensity, 240.0);
    assert_eq!(result.spectrum.peaks[1].intensity, 160.0);
    // Source sums the shared intensity 120 independently for each accepted cluster.
    // Thus total output 400 exceeds input 280; membership is not a partition.
    for (name, expected) in [
        ("charge", [1, 3]),
        ("iso_peak_count", [2, 2]),
        ("feature_number", [0, 1]),
    ] {
        let array = result
            .spectrum
            .integer_data_arrays
            .iter()
            .find(|a| a.name == name)
            .unwrap();
        assert_eq!(array.data, expected);
    }
    let disjoint = AveragineDeisotoper {
        allow_shared_isotopes: false,
        ..config
    }
    .deisotope(&input)
    .unwrap();
    assert_eq!(disjoint.clusters, result.clusters[..1]);
    assert_eq!(disjoint.spectrum.peaks, result.spectrum.peaks[..1]);
}

#[test]
fn nonpositive_model_mass_does_not_suppress_a_later_valid_ladder() {
    // A present first extension makes both zero and negative Poisson model
    // masses observable. Native rejection applies to this candidate alone.
    for seed in [0.0, 0.5, PROTON_MASS_U] {
        let mut input = ladder(seed, 1, &[120.0, 120.0]);
        input
            .peaks
            .extend(ladder(1800.0 + PROTON_MASS_U, 1, &[120.0, 120.0]).peaks);
        let result = AveragineDeisotoper {
            keep_only_deisotoped: false,
            ..settings(2)
        }
        .deisotope(&input)
        .unwrap();
        assert_eq!(result.clusters.len(), 1, "low-m/z seed={seed}");
        assert_eq!(result.clusters[0].peak_indices, [2, 3]);
        assert_eq!(result.spectrum.peaks, input.peaks[..3]);
    }
}
