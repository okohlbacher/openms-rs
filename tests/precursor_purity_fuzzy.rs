// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Independent analytical expectations for the defined PrecursorPurity.cpp loops.

use openms::analysis::precursor_purity::PrecursorPurity;
use openms::kernel::{MSExperiment, MSSpectrum, Peak1D, Precursor};

const NEUTRON: f64 = 1.008_664_915_66;

fn precursor(mz: f64, charge: i32, lower: f64, upper: f64) -> Precursor {
    Precursor {
        isolation_window_lower_offset: lower,
        isolation_window_upper_offset: upper,
        ..Precursor::new(mz, charge)
    }
}
fn spectrum(peaks: &[(f64, f32)]) -> MSSpectrum {
    MSSpectrum::from_peaks(peaks.iter().map(|&(mz, i)| Peak1D::new(mz, i)).collect())
}
fn experiment(peaks: &[(f64, f32)], precursors: Vec<Precursor>) -> MSExperiment {
    MSExperiment {
        spectra: vec![
            spectrum(peaks),
            MSSpectrum {
                ms_level: 2,
                precursors,
                ..MSSpectrum::default()
            },
        ],
        ..MSExperiment::default()
    }
}
fn compute(experiment: &MSExperiment, ppm: f64) -> Vec<f64> {
    PrecursorPurity::compute_single_scan(experiment, 1, 0, ppm).unwrap()
}

#[test]
fn strict_borders_have_half_weight_and_fuzzy_borders_are_excluded() {
    let lower = 100.0 - NEUTRON;
    let upper = 100.0 + NEUTRON;
    let e = experiment(
        &[
            (lower * (1.0 - 20.0 / 1e6), 100.0),
            (lower, 4.0),
            (99.5, 5.0),
            (100.0, 8.0),
            (upper, 6.0),
            (upper * (1.0 + 20.0 / 1e6), 100.0),
        ],
        vec![precursor(100.0, 1, NEUTRON, NEUTRON)],
    );
    // Target =8+4/2+6/2; total additionally includes interference5.
    assert_eq!(compute(&e, 20.0), [f64::from(13.0_f32 / 18.0)]);
    assert_eq!(e.spectra[0].peaks[0].intensity, 100.0);
}

#[test]
fn source_lower_bound_ignores_a_closer_predecessor() {
    let expected = 100.0 + NEUTRON;
    let e = experiment(
        &[
            (100.0, 1.0),
            (expected - 0.00001, 9.0),
            (expected + 0.01, 7.0),
        ],
        vec![precursor(100.0, 1, 0.0, 1.5)],
    );
    // The point just below the isotope is within20ppm, but the source tests
    // the first point ABOVE it instead. Neither point enters target intensity.
    assert_eq!(compute(&e, 20.0), [f64::from(1.0_f32 / 17.0)]);
}

#[test]
fn physical_end_fallback_uses_a_lone_candidate_and_handles_no_candidate() {
    let last_isotope = experiment(
        &[(100.0, 1.0), (100.0 + NEUTRON, 3.0)],
        vec![precursor(100.0, 1, 0.0, 1.1)],
    );
    assert_eq!(compute(&last_isotope, 20.0), [1.0]);
    let no_candidate = experiment(&[(100.0, 1.0)], vec![precursor(100.0, 1, 0.0, 2.0)]);
    assert_eq!(compute(&no_candidate, 20.0), [1.0]);
}

#[test]
fn duplicate_successor_ties_choose_the_second_peak_and_zero_deviation_is_strict() {
    let mut e = experiment(
        &[
            (100.0, 1.0),
            (100.0 + NEUTRON, 10.0),
            (100.0 + NEUTRON, 20.0),
        ],
        vec![precursor(100.0, 1, 0.0, 1.1)],
    );
    assert_eq!(compute(&e, 20.0), [f64::from(21.0_f32 / 31.0)]);
    e.spectra[1].precursors[0].charge = 0;
    assert_eq!(compute(&e, 20.0), [f64::from(21.0_f32 / 31.0)]);
    // Even an exact mass match fails source's strict ppm<deviation test at0.
    assert_eq!(compute(&e, 0.0), [f64::from(1.0_f32 / 31.0)]);
}

#[test]
fn a_physically_valid_candidate_past_the_logical_window_keeps_source_overcount() {
    let e = experiment(
        &[(100.0, 1.0), (101.02, 8.0), (102.0, 1.0)],
        vec![precursor(100.0, 1, 0.0, 0.999)],
    );
    // Fuzzy upper101.0191998 excludes101.02 from total. Source still reads
    // this valid physical iterator at logical end, matches it within200ppm of
    // 101.00866491566, and adds half its intensity to the target.
    assert_eq!(compute(&e, 200.0), [5.0]);
}

#[test]
fn zero_windows_stop_remaining_precursors_and_empty_parents_return_ones() {
    let p = precursor(100.0, 1, 0.0, 0.5);
    let mut e = experiment(
        &[(100.0, 1.0), (100.2, 3.0)],
        vec![p.clone(), precursor(100.0, 1, 0.0, 0.0), p],
    );
    assert_eq!(compute(&e, 0.0), [0.25, 1.0, 1.0]);
    e.spectra[0].peaks.clear();
    assert_eq!(compute(&e, 0.0), [1.0; 3]);
    e.spectra[1].precursors.clear();
    assert!(compute(&e, 0.0).is_empty());
}

#[test]
fn global_seed_and_f32_arithmetic_are_preserved() {
    let outside = experiment(
        &[(100.4, 1.0), (100.5, 2.0)],
        vec![precursor(100.0, 0, 0.1, 0.1)],
    );
    assert_eq!(compute(&outside, 0.0), [1.0]);
    let quotient = experiment(
        &[(100.0, 1.0), (100.1, 2.0)],
        vec![precursor(100.0, 1, 0.0, 0.2)],
    );
    assert_eq!(compute(&quotient, 0.0), [f64::from(1.0_f32 / 3.0)]);
    assert_ne!(compute(&quotient, 0.0)[0], 1.0_f64 / 3.0);
    let rounded = experiment(
        &[(99.9, 1.0), (100.0, 16_777_216.0), (100.1, 1.0)],
        vec![precursor(100.0, 1, 0.2, 0.2)],
    );
    assert_eq!(compute(&rounded, 0.0), [1.0]);
}

fn interpolation_experiment() -> MSExperiment {
    let mut e = experiment(
        &[(100.0, 1.0), (100.1, 3.0)],
        vec![precursor(100.0, 1, 0.0, 0.2)],
    );
    e.spectra[0].rt = 10.0;
    e.spectra[1].rt = 20.0;
    let mut later = spectrum(&[(100.0, 3.0), (100.1, 1.0)]);
    later.rt = 30.0;
    e.spectra.push(later);
    e
}

#[test]
fn interpolation_uses_absolute_rt_distance_and_keeps_source_extrapolation() {
    let mut e = interpolation_experiment();
    assert_eq!(
        PrecursorPurity::compute_interpolated(&e, 1, 0, Some(2), 0.0).unwrap(),
        [0.5]
    );
    e.spectra[1].rt = 50.0;
    assert_eq!(
        PrecursorPurity::compute_interpolated(&e, 1, 0, Some(2), 0.0).unwrap(),
        [1.25]
    );
    e.spectra[1].rt = -10.0;
    assert_eq!(
        PrecursorPurity::compute_interpolated(&e, 1, 0, Some(2), 0.0).unwrap(),
        [0.75]
    );
}

#[test]
fn interpolation_fallbacks_skip_an_unusable_late_scan() {
    let e = interpolation_experiment();
    for next in [None, Some(99), Some(1), Some(0)] {
        assert_eq!(
            PrecursorPurity::compute_interpolated(&e, 1, 0, next, 0.0).unwrap(),
            [0.25]
        );
    }
    let mut invalid_rt = e.clone();
    invalid_rt.spectra[2].rt = f64::NAN;
    assert_eq!(
        PrecursorPurity::compute_interpolated(&invalid_rt, 1, 0, Some(2), 0.0).unwrap(),
        [0.25]
    );
    invalid_rt.spectra[0].rt = -f64::MAX;
    invalid_rt.spectra[2].rt = f64::MAX;
    assert_eq!(
        PrecursorPurity::compute_interpolated(&invalid_rt, 1, 0, Some(2), 0.0).unwrap(),
        [0.25]
    );
}

#[test]
fn invalid_precursor_scalars_are_rejected_before_empty_parent_shortcut() {
    for empty_parent in [false, true] {
        for field in 0..5 {
            for invalid in [f64::NAN, f64::INFINITY, -1.0] {
                let mut e = interpolation_experiment();
                if empty_parent {
                    e.spectra[0].peaks.clear();
                }
                let precursor = &mut e.spectra[1].precursors[0];
                match field {
                    0 => precursor.mz = invalid,
                    1 => precursor.intensity = invalid as f32,
                    2 => precursor.isolation_window_lower_offset = invalid,
                    3 => precursor.isolation_window_upper_offset = invalid,
                    _ => precursor.isolation_target_mz = Some(invalid),
                }
                assert!(
                    PrecursorPurity::compute_single_scan(&e, 1, 0, 20.0).is_err(),
                    "empty_parent={empty_parent}, field={field}, value={invalid}"
                );
            }
        }
    }
}

#[test]
fn malformed_inputs_and_nonterminating_source_cases_fail_without_mutation() {
    let e = interpolation_experiment();
    assert!(PrecursorPurity::compute_single_scan(&e, 9, 0, 20.0).is_err());
    assert!(PrecursorPurity::compute_single_scan(&e, 1, 9, 20.0).is_err());
    for deviation in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(PrecursorPurity::compute_single_scan(&e, 1, 0, deviation).is_err());
    }
    let mut cases = Vec::new();
    let mut negative_charge = e.clone();
    negative_charge.spectra[1].precursors[0].charge = -1;
    cases.push(negative_charge);
    let mut negative_intensity = e.clone();
    negative_intensity.spectra[0].peaks[0].intensity = -1.0;
    cases.push(negative_intensity);
    let mut unsorted = e.clone();
    unsorted.spectra[0].peaks.reverse();
    cases.push(unsorted);
    let mut zero_total = e.clone();
    zero_total.spectra[0]
        .peaks
        .iter_mut()
        .for_each(|p| p.intensity = 0.0);
    cases.push(zero_total);
    let mut negative_bound = e.clone();
    negative_bound.spectra[1].precursors[0].isolation_window_lower_offset = 101.0;
    cases.push(negative_bound);
    let mut overflow = e.clone();
    overflow.spectra[0]
        .peaks
        .iter_mut()
        .for_each(|p| p.intensity = f32::MAX);
    cases.push(overflow);
    for case in cases {
        let before = case.clone();
        assert!(PrecursorPurity::compute_single_scan(&case, 1, 0, 20.0).is_err());
        assert_eq!(case, before);
    }
    let stalled = experiment(
        &[(1000.0, 1.0), (1001.0, 1.0)],
        vec![precursor(1000.0, 1, 2.0, 2.0)],
    );
    assert!(PrecursorPurity::compute_single_scan(&stalled, 1, 0, 5000.0).is_err());
    let huge_ladder = experiment(
        &[(1000.0, 1.0)],
        vec![precursor(1000.0, i32::MAX, 10.0, 10.0)],
    );
    assert!(PrecursorPurity::compute_single_scan(&huge_ladder, 1, 0, 20.0).is_err());
}
