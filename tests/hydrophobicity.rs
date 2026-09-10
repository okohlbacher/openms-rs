// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{
    AAIndex, AASequence, HydrophobicityProfile as Profile, HydrophobicityScale as Scale,
};

fn sequence(input: &str) -> AASequence {
    AASequence::parse(input).unwrap()
}
fn close(actual: &[f64], expected: &[f64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (&actual, &expected) in actual.iter().zip(expected) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual} != {expected}"
        );
    }
}

#[test]
fn source_scalar_profile_and_window_goldens() {
    // Literal HydrophobicityProfile_test.cpp assertions at pinned rev7c029e8.
    assert!((Profile::compute_gravy(&sequence("ACDE")).unwrap() + 0.675).abs() < 1e-15);
    assert_eq!(
        Profile::compute_profile(&sequence("ACDE"), Scale::default()).unwrap(),
        [1.8, 2.5, -3.5, -3.5]
    );
    let peptide = sequence("ACDEF");
    close(
        &Profile::compute_windowed_profile(&peptide, 3, Scale::KyteDoolittle).unwrap(),
        &[0.266666666666667, -1.5, -1.4],
        1e-14,
    );
    close(
        &Profile::compute_windowed_profile(&peptide, 6, Scale::KyteDoolittle).unwrap(),
        &[0.02],
        1e-14,
    );
    close(
        &Profile::compute_hydrophobic_moment(&peptide, 3, 100.0).unwrap(),
        &[0.511576803, 0.435170599, 0.734926405],
        1e-9,
    );
}

#[test]
fn every_scale_is_available_and_unknown_numeric_codes_error() {
    assert_eq!(Scale::ALL.len(), 7);
    assert_eq!(Scale::default(), Scale::KyteDoolittle);
    let alanine = [1.8, 0.62, -0.5, 0.61, 0.616, 0.1, 0.25];
    for (scale, value) in Scale::ALL.into_iter().zip(alanine) {
        assert_eq!(scale.value('A').unwrap(), value);
        for residue in ['B', 'J', 'O', 'U', 'X', 'Z', 'a', 'é', '\0'] {
            assert!(scale.value(residue).is_err());
        }
    }
    assert_eq!(Scale::Eisenberg.value('R').unwrap(), -2.53);
    assert_eq!(Scale::EisenbergConsensus.value('R').unwrap(), -1.76);
}

#[test]
fn rolling_update_order_is_preserved_and_oversize_windows_clamp_first() {
    let peptide = sequence("ACDEFGHIKLMNPQRSTVWYACDEFGHIK");
    for scale in Scale::ALL {
        let values = Profile::compute_profile(&peptide, scale).unwrap();
        assert_eq!(
            Profile::compute_windowed_profile(&peptide, 1, scale).unwrap(),
            values
        );
        let width = 7;
        let mut sum = values[..width]
            .iter()
            .copied()
            .fold(0.0, |sum, value| sum + value);
        let mut expected = vec![sum / width as f64];
        for i in width..values.len() {
            sum -= values[i - width];
            sum += values[i];
            expected.push(sum / width as f64);
        }
        assert_eq!(
            Profile::compute_windowed_profile(&peptide, width, scale).unwrap(),
            expected
        );
        assert_eq!(
            Profile::compute_windowed_profile(&peptide, usize::MAX, scale).unwrap(),
            Profile::compute_windowed_profile(&peptide, peptide.len(), scale).unwrap()
        );
    }
    assert_eq!(
        Profile::compute_hydrophobic_moment(&peptide, usize::MAX, 160.0).unwrap(),
        Profile::compute_hydrophobic_moment(&peptide, peptide.len(), 160.0).unwrap()
    );
}

#[test]
fn moment_phase_resets_and_eisenberg_scale_are_observable() {
    let peptide = sequence("ARAR");
    let moment = Profile::compute_hydrophobic_moment(&peptide, 2, 0.0).unwrap();
    // Mean normalized Eisenberg: (0.62 - 2.53)/2, absolute value.
    close(&moment, &[0.955; 3], 1e-15);
    let first = Profile::compute_hydrophobic_moment(&peptide, 2, 100.0).unwrap();
    assert_eq!(first[0], first[2]); // Identical windows have identical phase.
    close(
        &first,
        &Profile::compute_hydrophobic_moment(&peptide, 2, -100.0).unwrap(),
        1e-15,
    );
    close(
        &Profile::compute_hydrophobic_moment(&peptide, 1, 160.0).unwrap(),
        &[0.62, 2.53, 0.62, 2.53],
        1e-15,
    );
    // Source uses angle*PI/180, without remapping finite >360-degree inputs.
    assert!(
        Profile::compute_hydrophobic_moment(&peptide, 2, 460.0)
            .unwrap()
            .iter()
            .all(|x| x.is_finite())
    );
}

#[test]
fn annotations_are_ignored_even_when_composition_is_unknown() {
    let plain = sequence("ACMK");
    let annotated = sequence("(Acetyl)AC(Carbamidomethyl)M[+12.3456789]K.[+0.123456789]");
    assert!(annotated.formula().is_err());
    let before = annotated.clone();
    assert_eq!(
        Profile::compute_gravy(&annotated).unwrap(),
        Profile::compute_gravy(&plain).unwrap()
    );
    for scale in Scale::ALL {
        assert_eq!(
            Profile::compute_profile(&annotated, scale).unwrap(),
            Profile::compute_profile(&plain, scale).unwrap()
        );
        assert_eq!(
            Profile::compute_windowed_profile(&annotated, 2, scale).unwrap(),
            Profile::compute_windowed_profile(&plain, 2, scale).unwrap()
        );
    }
    assert_eq!(
        Profile::compute_hydrophobic_moment(&annotated, 2, 100.0).unwrap(),
        Profile::compute_hydrophobic_moment(&plain, 2, 100.0).unwrap()
    );
    assert_eq!(annotated, before);
}

#[test]
fn empty_zero_window_unsupported_codes_and_nonfinite_angles_are_errors() {
    let empty = sequence("");
    assert!(Profile::compute_gravy(&empty).is_err());
    assert!(Profile::compute_profile(&empty, Scale::Guy).is_err());
    assert!(Profile::compute_windowed_profile(&empty, 1, Scale::Guy).is_err());
    assert!(Profile::compute_hydrophobic_moment(&empty, 1, 100.0).is_err());
    let peptide = sequence("ACDE");
    assert!(Profile::compute_windowed_profile(&peptide, 0, Scale::Guy).is_err());
    assert!(Profile::compute_hydrophobic_moment(&peptide, 0, 100.0).is_err());
    for code in ['B', 'J', 'O', 'U', 'X', 'Z'] {
        let peptide = sequence(&format!("ACDE{code}"));
        assert!(Profile::compute_gravy(&peptide).is_err());
        assert!(Profile::compute_profile(&peptide, Scale::Guy).is_err());
        assert!(Profile::compute_windowed_profile(&peptide, 2, Scale::Guy).is_err());
        assert!(Profile::compute_hydrophobic_moment(&peptide, 2, 100.0).is_err());
    }
    for angle in [
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::MAX,
        -f64::MAX,
    ] {
        assert!(Profile::compute_hydrophobic_moment(&peptide, 1, angle).is_err());
    }
    let longer = sequence(&"A".repeat(1000));
    assert!(Profile::compute_hydrophobic_moment(&longer, 1000, f64::MAX / 4.0).is_err());
}

#[test]
fn quadratic_moment_work_and_linear_input_caps_are_checked() {
    let moderate = sequence(&"A".repeat(20_000));
    assert!(Profile::compute_hydrophobic_moment(&moderate, 10_000, 100.0).is_err());
    // The same sequence is valid for linear work and a clamped single window.
    assert_eq!(
        Profile::compute_profile(&moderate, Scale::Guy)
            .unwrap()
            .len(),
        20_000
    );
    assert_eq!(
        Profile::compute_hydrophobic_moment(&moderate, usize::MAX, 100.0)
            .unwrap()
            .len(),
        1
    );
    // One constructed oversized input exercises every entry point; limits are
    // checked before property traversal/output allocation and leave it intact.
    let oversized = sequence(&"A".repeat(1_000_001));
    assert!(AAIndex::calculate_gb(&oversized, 500.0).is_err());
    assert!(Profile::compute_gravy(&oversized).is_err());
    assert!(Profile::compute_profile(&oversized, Scale::Guy).is_err());
    assert!(Profile::compute_windowed_profile(&oversized, 7, Scale::Guy).is_err());
    assert!(Profile::compute_hydrophobic_moment(&oversized, 1, 100.0).is_err());
    assert_eq!(oversized.len(), 1_000_001);
}
