// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `IsotopeDistribution.h`: `merge`, `operator<`, `std::hash`, `Sorted` and the
//! validated in-place access.
//!
//! The ordering cases are the `operator <` section of
//! `IsotopeDistribution_test.cpp` at core `bc9cc12` (tier 2). The pinned class
//! test has no `merge` or `std::hash` section; the merge values are derived by
//! hand from `IsotopeDistribution.cpp:240-273` at the same pin on inputs chosen
//! to be exact in binary (tier 4), and the hash cases check the documented
//! consistency with `==`.

use openms::chemistry::EmpiricalFormula;
use openms::chemistry::isotopes::{
    CoarseIsotopePatternGenerator as Generator, CoarseMassMode, IsotopeDistribution, IsotopePeak,
    Sorted,
};
use std::collections::HashSet;

fn run(formula: &str, max: usize, mode: CoarseMassMode) -> IsotopeDistribution {
    Generator::new(Some(max), mode)
        .unwrap()
        .run(&EmpiricalFormula::parse(formula).unwrap())
        .unwrap()
}

fn distribution(peaks: &[(f64, f64)]) -> IsotopeDistribution {
    IsotopeDistribution::from_peaks(
        peaks
            .iter()
            .map(|&(mass, probability)| IsotopePeak { mass, probability })
            .collect(),
    )
    .unwrap()
}

fn pairs(distribution: &IsotopeDistribution) -> Vec<(f64, f64)> {
    distribution
        .peaks()
        .iter()
        .map(|p| (p.mass, p.probability))
        .collect()
}

#[test]
fn ordering_matches_the_source_cases() {
    let approximate = CoarseMassMode::Approximate;
    assert!(IsotopeDistribution::default() >= IsotopeDistribution::default());
    assert!(run("C4", 11, approximate) < run("C5", 11, approximate));
    assert!(run("C5", 1, approximate) < run("C5", 1000, approximate));
    // The second isotope's mass is 61 in nominal mode and 61.003 otherwise.
    assert!(run("C5", 11, CoarseMassMode::Nominal) < run("C5", 11, approximate));
}

#[test]
fn ordering_compares_size_then_mass_then_probability() {
    let small = distribution(&[(5.0, 1.0)]);
    let large = distribution(&[(1.0, 0.1), (2.0, 0.1)]);
    assert!(small < large);
    let lighter = distribution(&[(1.0, 0.9), (2.0, 0.1)]);
    let heavier = distribution(&[(1.0, 0.1), (3.0, 0.1)]);
    assert!(lighter < heavier);
    let weaker = distribution(&[(1.0, 0.1), (2.0, 0.1)]);
    assert!(weaker < lighter);
    assert_eq!(
        distribution(&[(0.0, 1.0)]).cmp(&distribution(&[(-0.0, 1.0)])),
        std::cmp::Ordering::Equal
    );
}

#[test]
fn hashes_follow_equality() {
    let approximate = CoarseMassMode::Approximate;
    let c4 = run("C4", 11, approximate);
    assert_eq!(
        c4.source_hash64(),
        run("C4", 11, approximate).source_hash64()
    );
    assert_ne!(
        c4.source_hash64(),
        run("C5", 11, approximate).source_hash64()
    );
    assert_ne!(
        c4.source_hash64(),
        run("C4", 11, CoarseMassMode::Nominal).source_hash64()
    );
    let mut set = HashSet::new();
    set.insert(c4.clone());
    set.insert(run("C4", 11, approximate));
    set.insert(distribution(&[(-0.0, 1.0)]));
    set.insert(distribution(&[(0.0, 1.0)]));
    assert_eq!(set.len(), 2);
}

#[test]
fn merge_bins_by_the_source_rule() {
    // range 1.5 at resolution 1: ceil -> 2 bins of width 0.75 from 100.
    // 100.4 rounds to bin 0; 101.0 to bin 1; 101.5 rounds to bin 2 and is dropped.
    let mut merged = distribution(&[(101.0, 0.125), (100.0, 0.5), (101.5, 0.125), (100.4, 0.25)]);
    merged.merge(1.0, 0.0).unwrap();
    assert_eq!(pairs(&merged), [(100.0, 0.75), (100.75, 0.125)]);

    // An exact multiple drops the heaviest peak: range 2, 2 bins, the peak at
    // 102 lands in bin 2.
    let mut exact = distribution(&[(100.0, 0.5), (101.0, 0.25), (101.25, 0.125), (102.0, 0.125)]);
    exact.merge(1.0, 0.0).unwrap();
    assert_eq!(pairs(&exact), [(100.0, 0.5), (101.0, 0.375)]);

    // Tails below min_prob are trimmed first, empty bins are removed by the
    // final trimIntensities.
    let mut trimmed = distribution(&[
        (99.0, 0.0625),
        (100.0, 0.5),
        (100.5, 0.25),
        (103.0, 0.25),
        (104.0, 0.0625),
    ]);
    // After trimming: 100, 100.5, 103; range 3 at resolution 1.5 -> 2 bins of
    // 1.5; 100 -> 0, 100.5 -> round(1/3) = 0, 103 -> 2 (dropped).
    trimmed.merge(1.5, 0.1).unwrap();
    assert_eq!(pairs(&trimmed), [(100.0, 0.75)]);

    // One remaining mass gives zero bins and an empty distribution.
    let mut single = distribution(&[(100.0, 1.0)]);
    single.merge(1.0, 0.0).unwrap();
    assert!(single.is_empty());
}

#[test]
fn merge_refusals_leave_the_distribution_unchanged() {
    let original = distribution(&[(100.0, 0.5), (100.5, 0.5)]);
    // range 0.5 at resolution 0.1: 5 bins for 2 peaks.
    let mut refined = original.clone();
    assert!(refined.merge(0.1, 0.0).is_err());
    assert_eq!(refined, original);
    for resolution in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(refined.merge(resolution, 0.0).is_err());
    }
    assert!(refined.merge(1.0, -0.5).is_err());
    assert!(IsotopeDistribution::empty().merge(1.0, 0.0).is_err());
    assert_eq!(refined, original);
}

#[test]
fn in_place_access_is_validated() {
    let mut edited = distribution(&[(1.0, 0.5), (2.0, 0.5)]);
    edited
        .set_peak(
            1,
            IsotopePeak {
                mass: 3.0,
                probability: 0.25,
            },
        )
        .unwrap();
    assert_eq!(pairs(&edited), [(1.0, 0.5), (3.0, 0.25)]);
    assert!(edited.set_peak(2, IsotopePeak::default()).is_err());
    assert!(
        edited
            .set_peak(
                0,
                IsotopePeak {
                    mass: -1.0,
                    probability: 0.1
                }
            )
            .is_err()
    );

    edited
        .try_for_each_peak_mut(|p| p.probability *= 2.0)
        .unwrap();
    assert_eq!(pairs(&edited), [(1.0, 1.0), (3.0, 0.5)]);
    let before = edited.clone();
    assert!(edited.try_for_each_peak_mut(|p| p.mass -= 2.0).is_err());
    assert_eq!(edited, before);

    assert_eq!(Sorted::default(), Sorted::Undefined);
    assert_ne!(Sorted::Intensity, Sorted::Mass);
}
