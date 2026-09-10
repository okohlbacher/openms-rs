// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Source: PrecursorPurity.cpp/class tests at OpenMS4-core7c029e8.

use openms::analysis::precursor_purity::{PrecursorPurity as Purity, PurityScores};
use openms::chemistry::{AASequence, C13C12_MASSDIFF_U, TheoreticalSpectrumGenerator};
use openms::comparison::Tolerance;
use openms::kernel::{DataArray, MSExperiment, MSSpectrum, Peak1D, Precursor};

fn spectrum(peaks: &[(f64, f32)]) -> MSSpectrum {
    MSSpectrum {
        peaks: peaks
            .iter()
            .map(|&(mz, intensity)| Peak1D::new(mz, intensity))
            .collect(),
        ..Default::default()
    }
}
fn precursor(mz: f64, charge: i32, lower: f64, upper: f64) -> Precursor {
    Precursor {
        isolation_window_lower_offset: lower,
        isolation_window_upper_offset: upper,
        ..Precursor::new(mz, charge)
    }
}

#[test]
fn exact_pinned_scalar_goldens_from_four_decoded_ms1_peaks() {
    // Exact binary64 coordinates and binary32 intensities independently decoded
    // from PrecursorPurity_input.mzML, spectrum6 originalindices273..276.
    let rows = [
        (0x407a58f47dec136a, 0x4b082b0b),
        (0x407a59334f8c8667, 0x489d2eec),
        (0x407a59a094ac7928, 0x495f5f29),
        (0x407a5c377ee6caba, 0x49aa88e2),
    ];
    let s = MSSpectrum {
        peaks: rows
            .into_iter()
            .map(|(mz, intensity)| Peak1D::new(f64::from_bits(mz), f32::from_bits(intensity)))
            .collect(),
        ..Default::default()
    };
    let p = precursor(421.558584883603, 3, 0.25, 0.25);
    let result = Purity::compute(&s, &p, Tolerance::Ppm(10.)).unwrap();
    assert_eq!(result.total_intensity, 11557777.1875);
    assert_eq!(result.target_intensity, 8923915.);
    assert!((result.signal_proportion - 0.77211).abs() < 0.000005);
    assert_eq!(result.target_peak_count, 1);
    assert_eq!(result.interfering_peak_count, 3);
    assert_eq!(result.interfering_peaks.peaks, s.peaks[1..]);
}

#[test]
fn missing_monoisotope_still_matches_isotopes_using_charge_magnitude() {
    let shift = C13C12_MASSDIFF_U / 2.;
    let s = spectrum(&[(100. - shift, 2.), (100.25, 5.), (100. + shift, 3.)]);
    for charge in [2, -2] {
        let result = Purity::compute(
            &s,
            &precursor(100., charge, 1.1, 1.1),
            Tolerance::Absolute(0.),
        )
        .unwrap();
        assert_eq!(result.total_intensity, 10.);
        assert_eq!(result.target_intensity, 5.);
        assert_eq!(result.signal_proportion, 0.5);
        assert_eq!(result.target_peak_count, 2);
        assert_eq!(result.interfering_peaks.peaks, [Peak1D::new(100.25, 5.)]);
    }
    for charge in [0, 1, i32::MIN] {
        let result = Purity::compute(
            &spectrum(&[(100., 7.)]),
            &precursor(100., charge, 0., 0.),
            Tolerance::Absolute(0.),
        )
        .unwrap();
        assert_eq!(result.signal_proportion, 1.);
        assert_eq!(result.target_peak_count, 1);
    }
}

#[test]
fn doubled_tolerance_and_isolation_boundaries_are_inclusive() {
    let p = precursor(2., 1, 0., 0.6);
    let at = Purity::compute(&spectrum(&[(2.5, 3.)]), &p, Tolerance::Absolute(0.25)).unwrap();
    assert_eq!(at.target_intensity, 3.);
    let outside = f64::from_bits(2.5f64.to_bits() + 1);
    let out = Purity::compute(&spectrum(&[(outside, 3.)]), &p, Tolerance::Absolute(0.25)).unwrap();
    assert_eq!(out.total_intensity, 3.);
    assert_eq!(out.target_peak_count, 0);
    assert_eq!(out.signal_proportion, 0.);
    let s = spectrum(&[(1.5, 1.), (2., 2.), (2.5, 3.)]);
    let result = Purity::compute(&s, &precursor(2., 1, 0.5, 0.5), Tolerance::Absolute(0.)).unwrap();
    assert_eq!(result.total_intensity, 6.);
    assert_eq!(result.target_intensity, 2.);
    // A large f64 ppm value is accepted when the actual arithmetic is finite.
    let huge_ppm = Purity::compute(
        &spectrum(&[(0., 2.)]),
        &precursor(0., 1, 0., 0.),
        Tolerance::Ppm(1e100),
    )
    .unwrap();
    assert_eq!(huge_ppm.signal_proportion, 1.);
}

#[test]
fn greedy_nearest_ties_and_removal_preserve_source_peak_selection() {
    let s = spectrum(&[(1.5, 2.), (2.5, 9.)]);
    let tied = Purity::compute(&s, &precursor(2., 1, 0.5, 0.5), Tolerance::Absolute(0.25)).unwrap();
    assert_eq!(tied.target_intensity, 2.);
    assert_eq!(tied.interfering_peaks.peaks, [Peak1D::new(2.5, 9.)]);
    // The first expected isotope consumes the nearer heavy observation;
    // the next isotope cannot reuse it even though its tolerance overlaps.
    let s = spectrum(&[(100.5, 4.)]);
    let consumed =
        Purity::compute(&s, &precursor(100., 1, 0., 2.), Tolerance::Absolute(0.3)).unwrap();
    assert_eq!(consumed.target_peak_count, 1);
    assert_eq!(consumed.target_intensity, 4.);
    // Exact duplicates use the first observation at equality, as lower_bound.
    let duplicates = spectrum(&[(100., 1.), (100., 9.)]);
    let matched = Purity::compute(
        &duplicates,
        &precursor(100., 1, 0., 0.),
        Tolerance::Absolute(0.),
    )
    .unwrap();
    assert_eq!(matched.target_intensity, 1.);
    assert_eq!(matched.interfering_peaks.peaks, [Peak1D::new(100., 9.)]);
}

#[test]
fn zero_windows_and_signed_zero_have_defined_scores_and_no_input_mutation() {
    let mut s = spectrum(&[(100., -0.), (100.2, 0.)]);
    s.native_id = "original".into();
    s.metadata.insert("instrument".into(), "test".into());
    s.float_data_arrays
        .push(DataArray::new("quality", vec![0.2, 0.3]));
    let p = precursor(100., 1, 0., 0.25);
    let before = s.clone();
    let before_p = p.clone();
    let result = Purity::compute(&s, &p, Tolerance::Absolute(0.)).unwrap();
    assert_eq!(result.total_intensity.to_bits(), 0.0f64.to_bits());
    assert_eq!(result.target_intensity.to_bits(), 0.0f64.to_bits());
    assert_eq!(result.signal_proportion, 0.);
    assert_eq!(result.target_peak_count, 1);
    assert_eq!(result.interfering_peak_count, 1);
    assert!(result.interfering_peaks.metadata.is_empty());
    assert!(result.interfering_peaks.float_data_arrays.is_empty());
    assert_eq!(s, before);
    assert_eq!(p, before_p);
    assert_eq!(
        Purity::compute(&s, &precursor(200., 1, 0., 0.5), Tolerance::Absolute(0.)).unwrap(),
        PurityScores::default()
    );
}

#[test]
fn invalid_inputs_and_unbounded_or_nonadvancing_isotope_work_fail_checked() {
    let p = precursor(100., 1, 0.5, 0.5);
    for s in [
        spectrum(&[(100., -1.)]),
        spectrum(&[(100., f32::NAN)]),
        spectrum(&[(f64::INFINITY, 1.)]),
        spectrum(&[(101., 1.), (100., 1.)]),
    ] {
        assert!(Purity::compute(&s, &p, Tolerance::Absolute(0.)).is_err());
    }
    let s = spectrum(&[(100., 1.)]);
    let before = s.clone();
    for invalid in [
        MSSpectrum {
            rt: f64::NAN,
            ..s.clone()
        },
        MSSpectrum {
            ms_level: 0,
            ..s.clone()
        },
    ] {
        assert!(Purity::compute(&invalid, &p, Tolerance::Absolute(0.)).is_err());
    }
    for invalid in [
        Tolerance::Absolute(-1.),
        Tolerance::Ppm(f64::INFINITY),
        Tolerance::Absolute(f64::MAX),
    ] {
        assert!(Purity::compute(&s, &p, invalid).is_err());
    }
    for bad in [
        precursor(100., 1, -1., 0.),
        precursor(100., 1, 0., 1e20),
        precursor(100., 1, f64::from(i32::MAX) + 1., 0.),
        precursor(100., 1, 0., f64::INFINITY),
        Precursor {
            isolation_target_mz: Some(f64::NAN),
            ..p.clone()
        },
        Precursor {
            intensity: f32::NAN,
            ..p.clone()
        },
    ] {
        assert!(Purity::compute(&s, &bad, Tolerance::Absolute(0.)).is_err());
    }
    let huge = spectrum(&[(1e100, 1.)]);
    assert!(Purity::compute(&huge, &precursor(1e100, 1, 0., 0.), Tolerance::Absolute(0.)).is_err());
    assert_eq!(s, before);
}

#[test]
fn repeated_parent_validation_does_not_traverse_unused_annotations() {
    let clean = spectrum(&[(100., 1.)]);
    let mut parent = clean.clone();
    parent.native_id = "parent".into();
    // These placeholders are unconsumed; repeatedly walking all 10,000 arrays
    // for 10,000 children would exceed the shared 50-million-unit work budget.
    parent.float_data_arrays = (0..10_000)
        .map(|i| DataArray::new(format!("unused{i}"), Vec::new()))
        .collect();
    // General kernel validation still rejects malformed auxiliary payloads.
    // Purity intentionally does not inspect or copy them.
    parent.float_data_arrays[0].data = vec![f32::NAN; 2];
    assert!(parent.validate().is_err());
    let mut p = precursor(100., 1, 0., 0.);
    p.activation_energy = f64::NAN;
    assert!(p.validate().is_err());
    let expected =
        Purity::compute(&clean, &precursor(100., 1, 0., 0.), Tolerance::Absolute(0.)).unwrap();
    assert_eq!(
        Purity::compute(&parent, &p, Tolerance::Absolute(0.)).unwrap(),
        expected
    );
    let mut experiment = MSExperiment::default();
    experiment.spectra.push(parent);
    for i in 0..10_000 {
        experiment.spectra.push(MSSpectrum {
            native_id: format!("child{i}"),
            ms_level: 2,
            precursors: vec![p.clone()],
            ..Default::default()
        });
    }
    let result = Purity::compute_all(&experiment, Tolerance::Absolute(0.), false).unwrap();
    assert_eq!(result.len(), 10_000);
    assert!(result.values().all(|score| score == &expected));
    assert_eq!(experiment.spectra[0].float_data_arrays.len(), 10_000);
    assert!(experiment.spectra[0].float_data_arrays[0].data[0].is_nan());
    assert!(
        experiment.spectra[1].precursors[0]
            .activation_energy
            .is_nan()
    );
}

#[test]
fn source_sps_goldens_charge_zero_and_duplicate_windows() {
    let peptide = AASequence::parse("PEPTIDER").unwrap();
    let precursors: Vec<_> = [425.20308, 419.18849, 730.37299, 500.]
        .into_iter()
        .map(|mz| Precursor::new(mz, 0))
        .collect();
    assert_eq!(
        Purity::count_sps_matches(&precursors, &peptide, Tolerance::Absolute(0.02), 1).unwrap(),
        3
    );
    assert_eq!(
        Purity::count_sps_matches(&precursors[..1], &peptide, Tolerance::Absolute(0.02), 0)
            .unwrap(),
        1
    );
    let doubly = [Precursor::new(430.21143, 0)];
    assert_eq!(
        Purity::count_sps_matches(&doubly, &peptide, Tolerance::Absolute(0.02), 1).unwrap(),
        0
    );
    assert_eq!(
        Purity::count_sps_matches(&doubly, &peptide, Tolerance::Absolute(0.02), 2).unwrap(),
        1
    );
    let ppm = [
        Precursor::new(425.20308 + 0.003, 1),
        Precursor::new(425.20308 + 0.05, 1),
    ];
    assert_eq!(
        Purity::count_sps_matches(&ppm, &peptide, Tolerance::Ppm(20.), 1).unwrap(),
        1
    );
    let repeated = vec![precursors[0].clone(); 3];
    assert_eq!(
        Purity::count_sps_matches(&repeated, &peptide, Tolerance::Absolute(0.02), 1).unwrap(),
        3
    );
    assert_eq!(
        Purity::count_sps_matches(
            &[],
            &AASequence::parse("BZX").unwrap(),
            Tolerance::Absolute(0.),
            1
        )
        .unwrap(),
        0
    );
    assert_eq!(
        Purity::count_sps_matches(
            &precursors,
            &AASequence::default(),
            Tolerance::Absolute(0.),
            1
        )
        .unwrap(),
        0
    );
}

#[test]
fn sps_bounds_narrow_to_f32_before_matching_and_theoretical_limits_propagate() {
    let p = AASequence::parse("A").unwrap();
    let mut masses = Vec::new();
    TheoreticalSpectrumGenerator::default()
        .append_mass_spectrum(&mut masses, &p, 1)
        .unwrap();
    let mass = f64::from(masses[0]);
    let ulp = f64::from(f32::from_bits(masses[0].to_bits() + 1)) - mass;
    let rounds_equal = [Precursor::new(mass + ulp * 0.25, 0)];
    let rounds_next = [Precursor::new(mass + ulp * 0.75, 0)];
    assert_eq!(
        Purity::count_sps_matches(&rounds_equal, &p, Tolerance::Absolute(0.), 1).unwrap(),
        1
    );
    assert_eq!(
        Purity::count_sps_matches(&rounds_next, &p, Tolerance::Absolute(0.), 1).unwrap(),
        0
    );
    assert!(
        Purity::count_sps_matches(
            &rounds_equal,
            &AASequence::parse("B").unwrap(),
            Tolerance::Absolute(0.),
            1
        )
        .is_err()
    );
    assert!(
        Purity::count_sps_matches(&[Precursor::new(1e100, 1)], &p, Tolerance::Absolute(0.), 1)
            .is_err()
    );
}
