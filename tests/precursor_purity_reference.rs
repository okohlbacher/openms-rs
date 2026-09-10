// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Independent source fixture and arithmetic references. No C++ runtime is used.
//! See data/precursor_purity_provenance.json for extraction and source hashes.

use openms::analysis::precursor_purity::{PrecursorPurity, PurityScores};
use openms::chemistry::{AASequence, TheoreticalSpectrumGenerator};
use openms::comparison::Tolerance;
use openms::kernel::{MSExperiment, MSSpectrum, Peak1D, Precursor};

const PEAKS: &str = include_str!("data/precursor_purity_peaks.tsv");
const SPECTRA: &str = include_str!("data/precursor_purity_spectra.tsv");
const SCORES: &str = include_str!("data/precursor_purity_scores.tsv");
const SPS: &str = include_str!("data/precursor_purity_sps.tsv");

fn rows(table: &str) -> impl Iterator<Item = Vec<&str>> {
    table
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .skip(1)
        .map(|line| line.split('\t').collect())
}
fn number(text: &str) -> f64 {
    text.parse().unwrap()
}
fn index(text: &str) -> usize {
    text.parse().unwrap()
}
fn f64_bits(text: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(text, 16).unwrap())
}
fn f32_bits(text: &str) -> f32 {
    f32::from_bits(u32::from_str_radix(text, 16).unwrap())
}
fn tolerance(value: &str, unit: &str) -> Tolerance {
    match unit {
        "Da" => Tolerance::Absolute(number(value)),
        "ppm" => Tolerance::Ppm(number(value)),
        _ => panic!("unknown fixture tolerance unit"),
    }
}
fn precursor(mz: f64, charge: i32, lower: f64, upper: f64) -> Precursor {
    Precursor {
        isolation_window_lower_offset: lower,
        isolation_window_upper_offset: upper,
        ..Precursor::new(mz, charge)
    }
}
fn fixture() -> MSExperiment {
    let mut experiment = MSExperiment::default();
    for row in rows(SPECTRA) {
        assert_eq!(index(row[0]), experiment.spectra.len());
        let mut spectrum = MSSpectrum {
            native_id: row[1].into(),
            ms_level: row[2].parse().unwrap(),
            rt: f64_bits(row[4]),
            ..Default::default()
        };
        assert_eq!(spectrum.rt, number(row[3]));
        if row[6] != "-" {
            let mut selected = precursor(
                f64_bits(row[7]),
                row[8].parse().unwrap(),
                number(row[11]),
                number(row[12]),
            );
            assert_eq!(selected.mz, number(row[6]));
            selected.intensity = f32_bits(row[10]);
            assert_eq!(f64::from(selected.intensity), number(row[9]));
            spectrum.precursors.push(selected);
        }
        experiment.spectra.push(spectrum);
    }
    for row in rows(PEAKS) {
        let peak = Peak1D::new(f64_bits(row[4]), f32_bits(row[5]));
        assert_eq!(peak.mz, number(row[2]));
        assert_eq!(f64::from(peak.intensity), number(row[3]));
        experiment.spectra[index(row[0])].peaks.push(peak);
    }
    assert_eq!(experiment.spectra[0].len(), 37);
    assert_eq!(experiment.spectra[6].len(), 35);
    experiment.validate().unwrap();
    experiment
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16} != {expected:.16} (absolute tolerance {tolerance})"
    );
}
fn check_score(actual: &PurityScores, row: &[&str]) {
    // Preserve truncated upstream decimal literals separately from exact decoded
    // f32 intensity sums and the independently calculated binary64 ratio.
    near(actual.total_intensity, number(row[5]), 0.000_01);
    near(actual.target_intensity, number(row[6]), 0.000_01);
    near(actual.signal_proportion, number(row[7]), 0.000_01);
    assert_eq!(actual.total_intensity, number(row[8]));
    assert_eq!(actual.target_intensity, number(row[9]));
    near(actual.signal_proportion, number(row[10]), 1e-14);
    assert_eq!(actual.target_peak_count, index(row[11]));
    assert_eq!(actual.interfering_peak_count, index(row[12]));
    let source_indices: Vec<usize> = if row[14] == "-" {
        vec![]
    } else {
        row[14].split(',').map(index).collect()
    };
    let expected: Vec<Peak1D> = rows(PEAKS)
        .filter(|peak| index(peak[0]) == index(row[1]) && source_indices.contains(&index(peak[1])))
        .map(|peak| Peak1D::new(f64_bits(peak[4]), f32_bits(peak[5])))
        .collect();
    assert_eq!(actual.interfering_peaks.peaks, expected);
    assert!(actual.interfering_peaks.native_id.is_empty());
}

#[test]
fn source_scalar_totals_counts_and_interfering_peaks_match_decoded_fixture() {
    let experiment = fixture();
    for row in rows(SCORES) {
        let score = PrecursorPurity::compute(
            &experiment.spectra[index(row[1])],
            &experiment.spectra[index(row[2])].precursors[0],
            tolerance(row[3], row[4]),
        )
        .unwrap();
        check_score(&score, &row);
    }
}

#[test]
fn source_batch_returns_five_ms2_entries_using_the_earlier_ms1() {
    let experiment = fixture();
    let scores =
        PrecursorPurity::compute_all(&experiment, Tolerance::Absolute(0.1), false).unwrap();
    assert_eq!(scores.len(), 5);
    for row in rows(SCORES).filter(|row| row[0].starts_with("experiment_")) {
        check_score(&scores[&experiment.spectra[index(row[2])].native_id], &row);
    }
    // C++ unordered_map::operator[] in the source test inserts MS1/random keys;
    // those zero-valued insertions are not actual algorithm output.
    assert!(!scores.contains_key(&experiment.spectra[0].native_id));
    assert!(!scores.contains_key(&experiment.spectra[6].native_id));
}

#[test]
fn source_sps_literal_fragment_references_and_charge_cases() {
    for row in rows(SPS) {
        let peptide: AASequence = if row[1] == "-" { "" } else { row[1] }.parse().unwrap();
        let selected: Vec<Precursor> = if row[2] == "-" {
            vec![]
        } else {
            row[2]
                .split(',')
                .map(|mz| Precursor::new(number(mz), 0))
                .collect()
        };
        let count = PrecursorPurity::count_sps_matches(
            &selected,
            &peptide,
            tolerance(row[3], row[4]),
            row[5].parse().unwrap(),
        )
        .unwrap();
        assert_eq!(count, index(row[6]), "{}", row[0]);
    }
}

#[test]
fn sps_zero_tolerance_uses_binary32_bounds_and_counts_duplicate_precursors() {
    // Independently calculated A b1: proton + (C3H7NO2 - H2O), with pinned
    // elemental masses, gives 72.044390626271 in f64 and bits 429016ba in f32.
    // A quarter-ULP perturbation still rounds to that same source f32 mass.
    let rounded = f32::from_bits(0x4290_16ba);
    let next = f32::from_bits(rounded.to_bits() + 1);
    let perturbed = f64::from(rounded) + f64::from(next - rounded) / 4.0;
    assert_ne!(perturbed, f64::from(rounded));
    let peptide: AASequence = "AG".parse().unwrap();
    let selected = [Precursor::new(perturbed, 7), Precursor::new(perturbed, 0)];
    assert_eq!(
        PrecursorPurity::count_sps_matches(&selected, &peptide, Tolerance::Absolute(0.0), 1)
            .unwrap(),
        2
    );
}

#[test]
fn scalar_inclusive_window_doubled_tolerance_and_lower_nearest_tie() {
    let spectrum = MSSpectrum::from_peaks(vec![
        Peak1D::new(499.75, 1.0),
        Peak1D::new(499.8, 7.0),
        Peak1D::new(500.2, 13.0),
        Peak1D::new(500.25, 2.0),
    ]);
    let score = PrecursorPurity::compute(
        &spectrum,
        &precursor(500.0, 0, 0.25, 0.25),
        Tolerance::Absolute(0.1),
    )
    .unwrap();
    assert_eq!(score.total_intensity, 23.0);
    assert_eq!(score.target_intensity, 7.0);
    assert_eq!(score.target_peak_count, 1);
    assert_eq!(score.interfering_peak_count, 3);
    assert_eq!(score.interfering_peaks.peaks[0].mz, 499.75);
    assert_eq!(score.interfering_peaks.peaks[2].mz, 500.25);
}

#[test]
fn scalar_can_match_an_isotope_without_a_monoisotopic_peak() {
    // The source header's absent-precursor statement contradicts its actual loop.
    let spectrum = MSSpectrum::from_peaks(vec![
        Peak1D::new(500.4, 20.0),
        Peak1D::new(500.0 + 1.003_354_837_8, 10.0),
    ]);
    for charge in [1, -1] {
        let score = PrecursorPurity::compute(
            &spectrum,
            &precursor(500.0, charge, 0.0, 1.1),
            Tolerance::Absolute(0.01),
        )
        .unwrap();
        assert_eq!(score.total_intensity, 30.0);
        assert_eq!(score.target_intensity, 10.0);
        assert_eq!(score.signal_proportion, 1.0 / 3.0);
        assert_eq!(score.target_peak_count, 1);
    }
}

fn single_input(peaks: Vec<Peak1D>, selected: Precursor) -> MSExperiment {
    MSExperiment {
        spectra: vec![
            MSSpectrum::from_peaks(peaks),
            MSSpectrum {
                ms_level: 2,
                precursors: vec![selected],
                ..Default::default()
            },
        ],
        ..Default::default()
    }
}

#[test]
fn fuzzy_source_accumulates_and_divides_in_binary32() {
    let input = single_input(
        vec![
            Peak1D::new(499.9, 1.0),
            Peak1D::new(500.0, 16_777_216.0),
            Peak1D::new(500.1, 1.0),
        ],
        precursor(500.0, 1, 0.2, 0.2),
    );
    assert_eq!(
        PrecursorPurity::compute_single_scan(&input, 1, 0, 20.0).unwrap(),
        [1.0]
    );
    let input = single_input(
        vec![Peak1D::new(499.9, 2.0), Peak1D::new(500.0, 1.0)],
        precursor(500.0, 1, 0.2, 0.2),
    );
    assert_eq!(
        PrecursorPurity::compute_single_scan(&input, 1, 0, 20.0).unwrap(),
        [f64::from(f32::from_bits(0x3eaa_aaab))]
    );
}

#[test]
fn fuzzy_source_strict_borders_are_half_weighted_and_outer_borders_excluded() {
    let low = 500.0 - 0.2;
    let high = 500.0 + 0.2;
    let ppm = 100.0;
    let input = single_input(
        vec![
            Peak1D::new(low * (1.0 - ppm / 1e6), 1_000.0),
            Peak1D::new(low, 10.0),
            Peak1D::new(500.0, 20.0),
            Peak1D::new(high, 30.0),
            Peak1D::new(high * (1.0 + ppm / 1e6), 1_000.0),
        ],
        precursor(500.0, 1, 0.2, 0.2),
    );
    assert_eq!(
        PrecursorPurity::compute_single_scan(&input, 1, 0, ppm).unwrap(),
        [0.5]
    );
}

#[test]
fn fuzzy_source_ignores_a_closer_isotope_below_its_lower_bound_candidate() {
    let expected_left = 500.0 - 1.008_664_915_66;
    let input = single_input(
        vec![
            Peak1D::new(expected_left - 0.0001, 100.0),
            Peak1D::new(499.5, 10.0),
            Peak1D::new(500.0, 20.0),
            Peak1D::new(501.1, 20.0),
            Peak1D::new(501.15, 10.0),
        ],
        precursor(500.0, 1, 1.2, 1.2),
    );
    // The close isotope is the predecessor. Source compares lower_bound and
    // its successor instead, so neither side accepts an isotope: 20 / 160.
    assert_eq!(
        PrecursorPurity::compute_single_scan(&input, 1, 0, 20.0).unwrap(),
        [0.125]
    );
}

#[test]
fn source_interpolation_fallback_assertions_hold_on_compact_fixture() {
    let experiment = fixture();
    let early = PrecursorPurity::compute_single_scan(&experiment, 1, 0, 20.0).unwrap();
    for next in [None, Some(999), Some(2)] {
        assert_eq!(
            PrecursorPurity::compute_interpolated(&experiment, 1, 0, next, 20.0).unwrap(),
            early
        );
    }
    let interpolated =
        PrecursorPurity::compute_interpolated(&experiment, 1, 0, Some(6), 20.0).unwrap();
    assert_eq!(interpolated.len(), early.len());
    assert!(interpolated.iter().all(|v| (0.0..=1.0).contains(v)));
}

#[test]
fn batch_reference_prefers_an_older_named_parent_then_falls_back_backward() {
    let mut experiment = MSExperiment {
        spectra: vec![
            MSSpectrum {
                native_id: "older".into(),
                peaks: vec![Peak1D::new(500.0, 10.0)],
                ..Default::default()
            },
            MSSpectrum {
                native_id: "nearest".into(),
                peaks: vec![Peak1D::new(500.0, 20.0)],
                ..Default::default()
            },
            MSSpectrum {
                native_id: "product".into(),
                ms_level: 2,
                precursors: vec![precursor(500.0, 1, 0.25, 0.25)],
                ..Default::default()
            },
            MSSpectrum {
                native_id: "later".into(),
                peaks: vec![Peak1D::new(500.0, 30.0)],
                ..Default::default()
            },
        ],
        ..Default::default()
    };
    for (reference, expected) in [
        (Some("older"), 10.0),
        (Some("missing"), 20.0),
        (Some("later"), 20.0),
        (None, 20.0),
    ] {
        experiment.spectra[2].precursors[0].spectrum_reference = reference.map(str::to_owned);
        let scores =
            PrecursorPurity::compute_all(&experiment, Tolerance::Absolute(0.0), false).unwrap();
        assert_eq!(scores.len(), 1);
        assert_eq!(scores["product"].target_intensity, expected);
    }
}

#[test]
fn source_sps_helper_self_consistency_has_32_full_length_b_y_values() {
    // Upstream uses this same self-consistency check; the independent count is
    // eight lengths times two series times two charges, including full length.
    let peptide: AASequence = "PEPTIDER".parse().unwrap();
    let mut mz = Vec::new();
    TheoreticalSpectrumGenerator::default()
        .append_mass_spectrum(&mut mz, &peptide, 2)
        .unwrap();
    assert_eq!(mz.len(), 32);
    let selected: Vec<_> = mz
        .into_iter()
        .map(|mz| Precursor::new(f64::from(mz), 0))
        .collect();
    assert_eq!(
        PrecursorPurity::count_sps_matches(&selected, &peptide, Tolerance::Absolute(0.001), 2)
            .unwrap(),
        32
    );
}

#[test]
fn source_interpolation_preserves_finite_extrapolation_outside_zero_to_one() {
    let mut input = single_input(
        vec![Peak1D::new(500.0, 20.0)],
        precursor(500.0, 1, 0.2, 0.2),
    );
    input.spectra[0].rt = 0.0;
    input.spectra[1].rt = 3.0;
    input.spectra.push(MSSpectrum {
        rt: 1.0,
        peaks: vec![Peak1D::new(499.9, 20.0), Peak1D::new(500.0, 20.0)],
        ..Default::default()
    });
    // Early purity 1, late purity 1/2, RT factor 3 => 1 + 3*(1/2 - 1).
    assert_eq!(
        PrecursorPurity::compute_interpolated(&input, 1, 0, Some(2), 20.0).unwrap(),
        [-0.5]
    );
    input.spectra[2].rt = input.spectra[0].rt;
    assert_eq!(
        PrecursorPurity::compute_interpolated(&input, 1, 0, Some(2), 20.0).unwrap(),
        [1.0]
    );
}
