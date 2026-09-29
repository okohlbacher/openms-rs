// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `CoarseIsotopePatternGenerator.h`: `setMaxIsotope`, `setRoundMasses`,
//! `getIsotopeOverrides` and `estimateForFragmentFromPeptideWeightAndS`.
//!
//! Expected values are the literals of `CoarseIsotopeDistribution_test.cpp` at
//! core `bc9cc12`, sections at lines 84-104 and 355-397 (tier 3: class-test
//! literals, compared with the class test's `TEST_REAL_SIMILAR` defaults of
//! absolute 1e-5 or ratio 1 + 1e-5), checked in both probability precisions.

use openms::chemistry::isotopes::{
    AveragineComposition, CoarseIsotopePatternGenerator as Generator, CoarseMassMode,
    IsotopeDistribution, IsotopePeak, ProbabilityPrecision,
};

/// `TEST_REAL_SIMILAR` with the class-test defaults: absolute difference
/// within 1e-5, or ratio within 1 + 1e-5.
fn similar(actual: f64, expected: f64) {
    let absdiff = (actual - expected).abs();
    let ratio = if actual.abs() > expected.abs() {
        actual / expected
    } else {
        expected / actual
    };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{actual} is not similar to {expected}"
    );
}

fn generator(precision: ProbabilityPrecision) -> Generator {
    Generator::new(None, CoarseMassMode::Approximate)
        .unwrap()
        .with_precision(precision)
}

const PRECISIONS: [ProbabilityPrecision; 2] = [
    ProbabilityPrecision::Double,
    ProbabilityPrecision::SourceSingle,
];

#[test]
fn set_mass_mode_and_set_max_peaks() {
    let mut solver = Generator::default();
    assert_eq!(solver.mass_mode(), CoarseMassMode::Approximate);
    solver.set_mass_mode(CoarseMassMode::Nominal);
    assert_eq!(solver.mass_mode(), CoarseMassMode::Nominal);

    let mut solver = generator(ProbabilityPrecision::Double);
    assert_eq!(solver.max_peaks(), None);
    assert_eq!(
        solver.estimate_from_peptide_weight(1234.2).unwrap().len(),
        317
    );
    solver.set_max_peaks(Some(117)).unwrap();
    assert_eq!(solver.max_peaks(), Some(117));
    assert!(solver.set_max_peaks(Some(0)).is_err());
    assert_eq!(solver.max_peaks(), Some(117));
    solver.set_max_peaks(None).unwrap();
    assert_eq!(solver.max_peaks(), None);
}

#[test]
fn isotope_overrides_lists_each_override_by_symbol() {
    let mut solver = Generator::default();
    assert_eq!(solver.isotope_overrides().len(), 0);
    let carbon = IsotopeDistribution::from_peaks(vec![
        IsotopePeak {
            mass: 12.0,
            probability: 0.5,
        },
        IsotopePeak {
            mass: 13.003355,
            probability: 0.5,
        },
    ])
    .unwrap();
    solver.set_isotope_override("C", carbon.clone()).unwrap();
    let listed: Vec<_> = solver.isotope_overrides().collect();
    assert_eq!(listed.len(), 1);
    assert_eq!(listed[0].0, "C");
    assert_eq!(listed[0].1, &carbon);
    solver.clear_isotope_overrides();
    assert_eq!(solver.isotope_overrides().len(), 0);
}

#[test]
fn sulfur_fragment_estimate_matches_the_class_test() {
    for precision in PRECISIONS {
        let solver = generator(precision);
        // Isolating the M+2 precursor isotope.
        for (sulfur, tail) in [
            (0, 0.355445559123552),
            (1, 0.900804974056174),
            (2, 0.947862830751023),
            (3, 0.969454586761089),
        ] {
            let mut iso = solver
                .estimate_fragment_from_peptide_weight_and_sulfur(
                    200.0,
                    sulfur,
                    100.0,
                    sulfur,
                    &[2],
                )
                .unwrap();
            iso.renormalize_with(precision).unwrap();
            similar(iso.peaks().last().unwrap().probability, tail);
        }

        // With no sulfur the formula matches the plain peptide averagine
        // estimate, so the masses agree exactly.
        let mut with_sulfur = solver
            .estimate_fragment_from_peptide_weight_and_sulfur(200.0, 0, 100.0, 0, &[2])
            .unwrap();
        with_sulfur.renormalize_with(precision).unwrap();
        let mut plain = solver
            .estimate_fragment_from_weights(200.0, 100.0, &[2], AveragineComposition::PEPTIDE)
            .unwrap();
        plain.renormalize_with(precision).unwrap();
        assert_eq!(with_sulfur.len(), plain.len());
        for (left, right) in with_sulfur.peaks().iter().zip(plain.peaks()) {
            assert_eq!(left.mass, right.mass);
        }
    }
}

#[test]
fn sulfur_fragment_estimate_refusals() {
    let solver = generator(ProbabilityPrecision::Double);
    // The source subtracts the unsigned sulfur counts and wraps.
    assert!(
        solver
            .estimate_fragment_from_peptide_weight_and_sulfur(200.0, 0, 100.0, 1, &[2])
            .is_err()
    );
    assert!(
        solver
            .estimate_fragment_from_peptide_weight_and_sulfur(200.0, 0, 100.0, 0, &[])
            .is_err()
    );
    assert!(
        solver
            .estimate_fragment_from_peptide_weight_and_sulfur(100.0, 0, 200.0, 0, &[0])
            .is_err()
    );
}
