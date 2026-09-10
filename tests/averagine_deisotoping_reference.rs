// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Complete real-spectrum golden from Deisotoper_test.cpp at pinned 7c029e8.
//! The input/output binary peaks were independently decoded without C++ or
//! Rust processing. Provenance and source/fixture hashes accompany the TSVs.

use openms::chemistry::{C13C12_MASSDIFF_U, PROTON_MASS_U};
use openms::comparison::Tolerance;
use openms::kernel::{DataArray, MSSpectrum, Peak1D, SpectrumType};
use openms::processing::deisotoping::AveragineDeisotoper;
use std::collections::{BTreeMap, BTreeSet};

fn rows(text: &str) -> impl Iterator<Item = Vec<&str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .skip(1)
        .map(|line| line.split('\t').collect())
}

fn peaks(text: &str) -> Vec<Peak1D> {
    rows(text)
        .enumerate()
        .map(|(index, row)| {
            assert_eq!(row.len(), 3);
            assert_eq!(row[0].parse::<usize>().unwrap(), index);
            Peak1D::new(row[1].parse().unwrap(), row[2].parse().unwrap())
        })
        .collect()
}

fn input() -> MSSpectrum {
    MSSpectrum {
        peaks: peaks(include_str!("data/averagine_deisotoping_in.tsv")),
        rt: 0.505_660_614,
        ms_level: 1,
        native_id: "controllerType=0 controllerNumber=1 scan=2".into(),
        spectrum_type: SpectrumType::Centroid,
        ..Default::default()
    }
}

fn input_with_markers() -> MSSpectrum {
    let mut input = input();
    input.integer_data_arrays.push(DataArray::new(
        "original_index",
        (0..input.len()).map(|i| i as i32).collect(),
    ));
    input.float_data_arrays.push(DataArray::new(
        "original_index_float",
        (0..input.len()).map(|i| i as f32 + 0.25).collect(),
    ));
    input.string_data_arrays.push(DataArray::new(
        "original_label",
        (0..input.len()).map(|i| format!("peak-{i}")).collect(),
    ));
    input
        .metadata
        .insert("source".into(), "pinned fixture".into());
    input
}

fn expected() -> Vec<Peak1D> {
    peaks(include_str!("data/averagine_deisotoping_out.tsv"))
}

fn assert_peaks(actual: &[Peak1D], expected: &[Peak1D]) {
    assert!(
        actual == expected,
        "ordered source peaks differ: actual count {}, expected count {}; unexpected {:?}; missing {:?}",
        actual.len(),
        expected.len(),
        actual
            .iter()
            .filter(|p| !expected.contains(p))
            .collect::<Vec<_>>(),
        expected
            .iter()
            .filter(|p| !actual.contains(p))
            .collect::<Vec<_>>(),
    );
}

/// Unique input seed and charge for each ordered source output, inferred from
/// unchanged intensity and the source charge-conversion equation, not Rust.
fn seeds() -> Vec<(usize, u8)> {
    rows(include_str!("data/averagine_deisotoping_seeds.tsv"))
        .enumerate()
        .map(|(index, row)| {
            assert_eq!(row.len(), 3);
            assert_eq!(row[0].parse::<usize>().unwrap(), index);
            (row[1].parse().unwrap(), row[2].parse().unwrap())
        })
        .collect()
}

fn options(top_n: Option<usize>) -> AveragineDeisotoper {
    // Exact source call: (..., 10.0, true, 5000, 1, 3, true); subsequent
    // defaults are min=2/max=10, single-charge=true and summed intensity=false.
    AveragineDeisotoper {
        tolerance: Tolerance::Ppm(10.0),
        top_n,
        min_charge: 1,
        max_charge: 3,
        keep_only_deisotoped: true,
        min_isotope_peaks: 2,
        max_isotope_peaks: 10,
        make_single_charged: true,
        ..Default::default()
    }
}

#[test]
fn every_real_spectrum_peak_matches_the_pinned_source_output() {
    let input = input();
    let before = input.clone();
    let expected = expected();
    let seeds = seeds();
    assert_eq!(input.len(), 5407);
    assert_eq!(
        input.peaks.iter().filter(|p| p.intensity == 0.0).count(),
        4514
    );
    assert_eq!(
        input.peaks.iter().filter(|p| p.intensity >= 0.05).count(),
        893
    );
    assert_eq!(expected.len(), 104);
    assert_eq!(seeds.len(), expected.len());
    assert!(input.precursors.is_empty());

    // Independently check that each recorded mapping is the sole source input
    // and charge reproducing BOTH source output fields exactly.
    for (output, &(seed, charge)) in expected.iter().zip(&seeds) {
        let candidates: Vec<_> = input
            .peaks
            .iter()
            .enumerate()
            .flat_map(|(index, peak)| {
                (1..=3_u8).filter_map(move |charge| {
                    let mz = peak.mz * f64::from(charge) - f64::from(charge - 1) * PROTON_MASS_U;
                    (peak.intensity == output.intensity && mz == output.mz)
                        .then_some((index, charge))
                })
            })
            .collect();
        assert_eq!(candidates, [(seed, charge)]);
    }

    // Both settings remove the source ThresholdMower's low peaks. Since only
    // 893 peaks survive, 5000 and no peak limit must produce identical spectra.
    for top_n in [Some(5000), None] {
        let result = options(top_n).deisotope(&input).unwrap();
        assert_peaks(&result.spectrum.peaks, &expected);
        assert_eq!(result.spectrum.rt, input.rt);
        assert_eq!(result.spectrum.ms_level, input.ms_level);
        assert_eq!(result.spectrum.native_id, input.native_id);
        assert_eq!(result.spectrum.spectrum_type, input.spectrum_type);
        let actual_seeds: Vec<_> = result
            .clusters
            .iter()
            .map(|cluster| (cluster.peak_indices[0], cluster.charge))
            .collect();
        let mut expected_seeds = seeds.clone();
        expected_seeds.sort_unstable();
        assert_eq!(actual_seeds, expected_seeds);
        result.spectrum.validate().unwrap();
    }
    assert_eq!(input, before);
}

#[test]
fn source_seeds_keep_original_arrays_aligned_after_threshold_and_charge_sort() {
    // The source files contain no auxiliary arrays. These independent identity
    // markers exercise the documented source select/sort behavior across all
    // 4514 removed zeros and the subsequent charge-induced permutation.
    let input = input_with_markers();
    let before = input.clone();
    let result = AveragineDeisotoper {
        annotate_charge: true,
        annotate_isotope_peak_count: true,
        annotate_features: true,
        ..options(Some(5000))
    }
    .deisotope(&input)
    .unwrap();
    assert_eq!(input, before);
    assert_peaks(&result.spectrum.peaks, &expected());
    assert_eq!(result.spectrum.metadata, input.metadata);
    let arrays: BTreeMap<_, _> = result
        .spectrum
        .integer_data_arrays
        .iter()
        .map(|a| (a.name.as_str(), a.data.as_slice()))
        .collect();
    let clusters: BTreeMap<_, _> = result
        .clusters
        .iter()
        .enumerate()
        .map(|(index, cluster)| (cluster.peak_indices[0], (index, cluster)))
        .collect();
    assert_eq!(clusters.len(), 104);
    // The real source fixture contains overlapping ladders: these seeds at
    // charges 1 and 3 both use peak 4578. Their two-peak KL divergences are
    // 0.0182303842 and 0.0250655282 (<0.05), and neither has another extension
    // within 10 ppm. Preventing reuse loses the source's m/z 1931.072690310334.
    assert_eq!(clusters[&4550].1.peak_indices, [4550, 4578]);
    assert_eq!(clusters[&4564].1.peak_indices, [4564, 4578]);
    for (output, (seed, charge)) in seeds().into_iter().enumerate() {
        assert_eq!(arrays["original_index"][output], seed as i32);
        assert_eq!(arrays["charge"][output], i32::from(charge));
        assert_eq!(
            result.spectrum.float_data_arrays[0].data[output],
            seed as f32 + 0.25
        );
        assert_eq!(
            result.spectrum.string_data_arrays[0].data[output],
            format!("peak-{seed}")
        );
        let (feature, cluster) = clusters[&seed];
        assert_eq!(cluster.charge, charge);
        assert_eq!(arrays["feature_number"][output], feature as i32);
        // Full memberships are not stored by the source fixture: check their
        // original-coordinate, spacing and annotation consistency separately.
        assert!((2..=10).contains(&cluster.peak_indices.len()));
        assert_eq!(
            arrays["iso_peak_count"][output],
            cluster.peak_indices.len() as i32
        );
        assert!(cluster.peak_indices.windows(2).all(|p| p[0] < p[1]));
        for (isotope, &index) in cluster.peak_indices.iter().enumerate() {
            let observed = input.peaks[index];
            let expected =
                input.peaks[seed].mz + isotope as f64 * C13C12_MASSDIFF_U / f64::from(charge);
            assert!(observed.intensity >= 0.05);
            assert!((observed.mz - expected).abs() <= input.peaks[seed].mz * 10.0 / 1e6);
        }
    }
    result.spectrum.validate().unwrap();
}

#[test]
fn disjoint_policy_removes_only_the_independently_identified_shared_ladder() {
    let input = input_with_markers();
    let before = input.clone();
    let result = AveragineDeisotoper {
        allow_shared_isotopes: false,
        annotate_charge: true,
        ..options(Some(5000))
    }
    .deisotope(&input)
    .unwrap();
    // This is an explicit native policy comparison, not a changed source
    // golden: seed 4564 cannot take isotope 4578 from the earlier seed 4550.
    let seeds = seeds();
    let expected: Vec<_> = expected()
        .into_iter()
        .zip(&seeds)
        .filter_map(|(peak, &(seed, _))| (seed != 4564).then_some(peak))
        .collect();
    assert_eq!(expected.len(), 103);
    assert_peaks(&result.spectrum.peaks, &expected);
    let arrays: BTreeMap<_, _> = result
        .spectrum
        .integer_data_arrays
        .iter()
        .map(|a| (a.name.as_str(), a.data.as_slice()))
        .collect();
    for (output, &(seed, charge)) in seeds.iter().filter(|(seed, _)| *seed != 4564).enumerate() {
        assert_eq!(arrays["original_index"][output], seed as i32);
        assert_eq!(arrays["charge"][output], i32::from(charge));
        assert_eq!(
            result.spectrum.float_data_arrays[0].data[output],
            seed as f32 + 0.25
        );
        assert_eq!(
            result.spectrum.string_data_arrays[0].data[output],
            format!("peak-{seed}")
        );
    }
    let mut used = BTreeSet::new();
    for cluster in &result.clusters {
        assert_ne!(cluster.peak_indices[0], 4564);
        for &index in &cluster.peak_indices {
            assert!(used.insert(index), "isotope {index} was assigned twice");
        }
    }
    assert_eq!(result.spectrum.metadata, input.metadata);
    assert_eq!(input, before);
    result.spectrum.validate().unwrap();
}
