// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Golden cases: OpenMS4-core 7c029e8 CoarseIsotopeDistribution_test.cpp and
// IsotopeDistribution_test.cpp. Additional checks exercise native error handling.

use openms::chemistry::isotopes::{
    AveragineComposition, CoarseIsotopePatternGenerator as Generator, CoarseMassMode,
    IsotopeDistribution, IsotopePeak, MAX_ISOTOPE_PEAKS,
};
use openms::chemistry::{C13C12_MASSDIFF_U, EmpiricalFormula, PROTON_MASS_U, element};

fn formula(text: &str) -> EmpiricalFormula {
    text.parse().unwrap()
}
fn solver(max: Option<usize>) -> Generator {
    Generator::new(max, CoarseMassMode::Approximate).unwrap()
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
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.15} != {expected:.15} (tol {tolerance})"
    );
}
fn probabilities(distribution: &IsotopeDistribution) -> Vec<f64> {
    distribution.peaks().iter().map(|p| p.probability).collect()
}

#[test]
fn distribution_identity_queries_and_validated_container() {
    let mut identity = IsotopeDistribution::default();
    assert_eq!(
        identity.peaks(),
        &[IsotopePeak {
            mass: 0.0,
            probability: 1.0
        }]
    );
    identity.resize(3).unwrap();
    assert_eq!(identity.len(), 3);
    assert_eq!(identity.peaks()[2], IsotopePeak::default());
    identity.clear();
    assert!(identity.is_empty());
    assert_eq!(identity.min_mass(), None);
    assert_eq!(identity.max_mass(), None);
    assert_eq!(identity.most_abundant(), None);
    assert!(identity.average_mass().is_err());
    identity.renormalize().unwrap();
    identity
        .insert(IsotopePeak {
            mass: 10.0,
            probability: 3.0,
        })
        .unwrap();
    identity
        .insert(IsotopePeak {
            mass: 12.0,
            probability: 1.0,
        })
        .unwrap();
    near(identity.average_mass().unwrap(), 10.5, 1e-12);
    identity.renormalize().unwrap();
    assert_eq!(probabilities(&identity), [0.75, 0.25]);
    assert_eq!(identity.min_mass(), Some(10.0));
    assert_eq!(identity.max_mass(), Some(12.0));
    assert_eq!(identity.most_abundant().unwrap().mass, 10.0);
    let old = identity.clone();
    assert!(
        identity
            .insert(IsotopePeak {
                mass: f64::NAN,
                probability: 1.0
            })
            .is_err()
    );
    assert_eq!(identity, old);
    assert!(identity.resize(MAX_ISOTOPE_PEAKS + 1).is_err());
    assert_eq!(identity, old);
    for peak in [
        IsotopePeak {
            mass: -1.0,
            probability: 1.0,
        },
        IsotopePeak {
            mass: 1.0,
            probability: -1.0,
        },
        IsotopePeak {
            mass: 1.0,
            probability: f64::INFINITY,
        },
    ] {
        assert!(IsotopeDistribution::from_peaks(vec![peak]).is_err());
    }
}

#[test]
fn trimming_is_inclusive_and_preserves_interior_gaps() {
    let raw = distribution(&[(0.0, 0.01), (1.0, 0.2), (2.0, 0.0), (3.0, 0.1), (4.0, 0.01)]);
    let mut both = raw.clone();
    both.trim_left(0.1).unwrap();
    both.trim_right(0.1).unwrap();
    assert_eq!(probabilities(&both), [0.2, 0.0, 0.1]);
    near(both.probability_sum(), 0.3, 1e-12);
    both.trim_intensities(0.1).unwrap();
    assert_eq!(probabilities(&both), [0.2, 0.1]);
    both.sort_by_probability();
    assert_eq!(both.peaks()[0].mass, 1.0);
    both.sort_by_mass();
    assert_eq!(both.peaks()[1].mass, 3.0);
    let old = both.clone();
    assert!(both.trim_right(f64::NAN).is_err());
    assert_eq!(both, old);
    for mut all in [raw.clone(), raw.clone()] {
        all.trim_left(1.0).unwrap();
        assert!(all.is_empty());
    }
    let mut all = raw;
    all.trim_right(1.0).unwrap();
    assert!(all.is_empty());
    let mut zero = distribution(&[(0.0, 0.0)]);
    let old = zero.clone();
    assert!(zero.renormalize().is_err());
    assert_eq!(zero, old);
    assert!(zero.average_mass().is_err());
}

#[test]
fn upstream_glucose_golden_and_positive_charge_convention() {
    let generator = solver(Some(3));
    let neutral = generator.run(&formula("C6H12O6")).unwrap();
    assert_eq!(neutral.len(), 3);
    near(neutral.peaks()[0].mass, 180.0633903828, 1e-10);
    near(neutral.peaks()[0].probability, 0.923456, 1e-6);
    near(neutral.peaks()[2].mass, 182.0701, 1e-4);
    near(neutral.peaks()[2].probability, 0.013232, 1e-6);
    near(neutral.probability_sum(), 1.0, 1e-14);
    let charged = generator.run(&formula("C6H12O6+2")).unwrap();
    near(charged.peaks()[0].mass, 182.077943, 1e-6);
    near(charged.peaks()[0].probability, 0.923246, 1e-6);
    near(charged.peaks()[2].mass, 184.0846529, 1e-6);
    near(charged.peaks()[2].probability, 0.0132435, 1e-6);
    let explicit = generator.run(&formula("C6H14O6")).unwrap();
    for (a, b) in charged.peaks().iter().zip(explicit.peaks()) {
        near(a.probability, b.probability, 1e-14);
        near(
            b.mass - a.mass,
            2.0 * (element("H").unwrap().mono_mass() - PROTON_MASS_U),
            1e-12,
        );
    }
    assert!(generator.run(&formula("C6H12O6-2")).is_err());
    assert!(generator.run(&formula("C-1H2")).is_err());
    assert_eq!(
        generator.run(&formula("")).unwrap(),
        IsotopeDistribution::default()
    );
}

#[test]
fn upstream_heavy_formula_and_bromine_gap_goldens() {
    let large = solver(Some(11)).run(&formula("C222N190O110")).unwrap();
    let expected = [
        0.0349429, 0.109888, 0.180185, 0.204395, 0.179765, 0.130358, 0.0809864, 0.0442441,
        0.0216593, 0.00963707, 0.0039406,
    ];
    for (index, (&expected, peak)) in expected.iter().zip(large.peaks()).enumerate() {
        assert_eq!(peak.mass.round(), 7084.0 + index as f64);
        near(peak.probability, expected, 1e-6);
    }
    for (text, expected) in [
        ("Br2", vec![0.25694761, 0.0, 0.49990478, 0.0, 0.24314761]),
        (
            "CBr2",
            vec![
                0.254198270573,
                0.002749339427,
                0.494555798854,
                0.005348981146,
                0.240545930573,
                0.002601679427,
            ],
        ),
    ] {
        let actual = solver(None).run(&formula(text)).unwrap();
        assert_eq!(actual.len(), expected.len());
        for (peak, value) in actual.peaks().iter().zip(expected) {
            near(peak.probability, value, 1e-10);
        }
    }
}

#[test]
fn lightest_mass_labels_and_nominal_mode_follow_source_conventions() {
    let approximate = solver(None);
    let nominal = Generator::new(None, CoarseMassMode::Nominal).unwrap();
    let selenium = approximate.run(&formula("Se")).unwrap();
    let lightest = element("Se").unwrap().isotopes()[0].mass;
    near(selenium.peaks()[0].mass, lightest, 1e-12);
    assert!(lightest < element("Se").unwrap().mono_mass());
    let labeled = approximate.run(&formula("(13)C2D2")).unwrap();
    assert_eq!(labeled.len(), 1);
    near(
        labeled.peaks()[0].mass,
        formula("(13)C2D2").mono_mass(),
        1e-12,
    );
    assert_eq!(labeled.peaks()[0].probability, 1.0);
    let a = approximate.run(&formula("C160")).unwrap();
    let n = nominal.run(&formula("C160")).unwrap();
    for (a, n) in a.peaks().iter().zip(n.peaks()) {
        assert_eq!(n.mass, a.mass.round());
        assert_eq!(n.probability, a.probability);
    }
    // Nominal labels round the corrected carbon-spacing mass rather than simply
    // adding integer indices to the rounded lightest mass.
    assert_eq!(n.peaks()[160].mass, 2081.0);
    assert_eq!(approximate.run(&formula("H2")).unwrap().len(), 5); // includes zero tritium tail
}

#[test]
fn raw_convolution_power_and_truncation_probability_conservation() {
    let carbon = distribution(&[(12.0, 0.9893), (13.003355, 0.0107)]);
    let full = solver(None);
    let squared = full.convolve(&carbon, &carbon).unwrap();
    assert_eq!(squared, full.convolve_power(&carbon, 2).unwrap());
    assert_eq!(
        full.convolve_power(&carbon, 0).unwrap(),
        IsotopeDistribution::default()
    );
    near(squared.probability_sum(), 1.0, 1e-12);
    near(squared.peaks()[0].probability, 0.9893_f64.powi(2), 1e-14);
    near(squared.peaks()[1].probability, 2.0 * 0.9893 * 0.0107, 1e-14);
    let truncated = solver(Some(2)).convolve(&carbon, &carbon).unwrap();
    assert!(truncated.probability_sum() < 1.0);
    assert_eq!(truncated.peaks(), &squared.peaks()[..2]);
    let generated = solver(Some(2)).run(&formula("C2")).unwrap();
    for (raw, normalized) in truncated.peaks().iter().zip(generated.peaks()) {
        near(
            normalized.probability,
            raw.probability / truncated.probability_sum(),
            1e-12,
        );
    }
    let doubled = distribution(&[(0.0, 2.0), (1.0, 1.0)]);
    near(
        full.convolve(&doubled, &doubled).unwrap().probability_sum(),
        9.0,
        1e-12,
    );
    let empty = IsotopeDistribution::empty();
    assert!(full.convolve(&empty, &carbon).unwrap().is_empty());
    assert!(full.convolve_power(&empty, 0).unwrap().is_empty());
}

#[test]
fn local_isotope_overrides_preserve_bins_and_do_not_change_global_elements() {
    let mut generator = solver(None);
    let natural = generator.run(&formula("C2")).unwrap();
    generator
        .set_isotope_override("C", distribution(&[(12.0, 0.0), (13.0, 1.0)]))
        .unwrap();
    let heavy = generator.run(&formula("C2")).unwrap();
    assert_eq!(probabilities(&heavy), [0.0, 0.0, 1.0]);
    near(heavy.peaks()[2].mass, 24.0 + 2.0 * C13C12_MASSDIFF_U, 1e-12);
    assert_eq!(solver(None).run(&formula("C2")).unwrap(), natural);
    // Explicit labeled isotopes are separate override keys.
    assert_eq!(
        generator.run(&formula("(13)C2")).unwrap().peaks()[0].probability,
        1.0
    );
    assert!(
        generator
            .set_isotope_override("C", distribution(&[(13.0, 1.0)]))
            .is_err()
    );
    assert!(
        generator
            .set_isotope_override("C", distribution(&[(12.0, 0.5), (14.0, 0.5)]))
            .is_err()
    );
    assert!(
        generator
            .set_isotope_override("C", IsotopeDistribution::empty())
            .is_err()
    );
    assert!(
        generator
            .set_isotope_override("C2", natural.clone())
            .is_err()
    );
    generator.clear_isotope_overrides();
    assert_eq!(generator.run(&formula("C2")).unwrap(), natural);
    let mut truncated = solver(Some(1));
    truncated
        .set_isotope_override("C", distribution(&[(12.0, 0.0), (13.0, 1.0)]))
        .unwrap();
    assert!(truncated.run(&formula("C2")).is_err()); // retained total is zero
}

#[test]
fn upstream_averagine_peptide_rna_dna_goldens() {
    let generator = solver(Some(3));
    for (mass, probability, first_mass) in [
        (100.0, 0.949735, 100.170),
        (1000.0, 0.586906, 999.714),
        (10000.0, 0.046495, 9994.041),
    ] {
        let pattern = generator.estimate_from_peptide_weight(mass).unwrap();
        near(pattern.peaks()[0].probability, probability, 2e-6);
        near(pattern.peaks()[0].mass, first_mass, 0.001);
    }
    for (mass, rna, dna) in [
        (100.0, 0.958166, 0.958166),
        (1000.0, 0.668538, 0.657083),
        (10000.0, 0.080505, 0.075138),
    ] {
        near(
            generator.estimate_from_rna_weight(mass).unwrap().peaks()[0].probability,
            rna,
            2e-6,
        );
        near(
            generator.estimate_from_dna_weight(mass).unwrap().peaks()[0].probability,
            dna,
            2e-6,
        );
    }
    for mass in [100.0, 1000.0, 10000.0] {
        let estimate = AveragineComposition::PEPTIDE
            .estimate_mono_mass(mass)
            .unwrap();
        assert!(estimate.hydrogen_adjustment_succeeded);
        assert!(
            (estimate.formula.mono_mass() - mass).abs() <= element("H").unwrap().mono_mass() / 2.0
        );
        assert_eq!(
            generator.estimate_from_peptide_mono_weight(mass).unwrap(),
            generator.run(&estimate.formula).unwrap()
        );
    }
    let pure_carbon = AveragineComposition {
        carbon: 1.0,
        hydrogen: 0.0,
        nitrogen: 0.0,
        oxygen: 0.0,
        sulfur: 0.0,
        phosphorus: 0.0,
    };
    let too_small = pure_carbon.estimate_mono_mass(6.0).unwrap();
    assert!(!too_small.hydrogen_adjustment_succeeded);
    assert_eq!(too_small.formula, formula("C1"));
    let zero = AveragineComposition::PEPTIDE
        .estimate_average_mass(0.0)
        .unwrap();
    assert!(zero.formula.is_empty());
}

#[test]
fn upstream_fixed_sulfur_goldens() {
    let generator = solver(Some(3));
    for (sulfur, tail) in [
        (0, 0.00290370998965918),
        (1, 0.0439547771832361),
        (2, 0.0804989104418586),
        (3, 0.117023432503842),
    ] {
        let pattern = generator
            .estimate_from_peptide_weight_and_sulfur(100.0, sulfur)
            .unwrap();
        near(pattern.peaks()[2].probability, tail, 1e-6);
        let estimate = AveragineComposition::PEPTIDE
            .estimate_average_mass_with_sulfur(100.0, sulfur)
            .unwrap();
        assert_eq!(estimate.formula.count("S").unwrap(), sulfur as i32);
    }
    assert!(
        generator
            .estimate_from_peptide_weight_and_sulfur(100.0, 4)
            .is_err()
    );
}

#[test]
fn fragment_isolation_matches_independent_conditional_probability_and_bounds() {
    let fragment = distribution(&[(12.0, 0.9), (13.0, 0.1)]);
    let complementary = distribution(&[(24.0, 0.8), (25.0, 0.2)]);
    let mut result = solver(None)
        .calc_fragment_isotope_dist(&fragment, &complementary, &[1], 12.0)
        .unwrap();
    near(result.peaks()[0].probability, 0.18, 1e-14);
    near(result.peaks()[1].probability, 0.08, 1e-14);
    result.renormalize().unwrap();
    near(result.peaks()[0].probability, 9.0 / 13.0, 1e-14);
    near(result.peaks()[1].probability, 4.0 / 13.0, 1e-14);
    let bounded = solver(Some(1))
        .calc_fragment_isotope_dist(&fragment, &complementary, &[1, 1], 12.0)
        .unwrap();
    assert_eq!(bounded.len(), 1);
    near(bounded.peaks()[0].probability, 0.18, 1e-14);
    assert!(
        solver(None)
            .calc_fragment_isotope_dist(&fragment, &complementary, &[], 12.0)
            .is_err()
    );
    let mut impossible = solver(None)
        .calc_fragment_isotope_dist(&fragment, &complementary, &[9], 12.0)
        .unwrap();
    assert!(impossible.renormalize().is_err());
    // Source C/C2 conditional golden.
    let generator = solver(None);
    let c = generator.run(&formula("C")).unwrap();
    let c2 = generator.run(&formula("C2")).unwrap();
    let mut conditioned = generator
        .calc_fragment_isotope_dist(&c, &c2, &[0, 1], 12.0)
        .unwrap();
    conditioned.renormalize().unwrap();
    near(conditioned.peaks()[0].probability, 0.989524, 4e-6);
    near(conditioned.peaks()[1].probability, 0.010479, 4e-6);
}

#[test]
fn fragment_averagine_estimates_match_upstream_goldens() {
    let generator = solver(None);
    for (precursor, fragment, expected) in [
        (200.0, 100.0, 0.954654801320083),
        (2000.0, 100.0, 0.975984866212216),
        (20000.0, 100.0, 0.995783521351781),
        (2000.0, 1000.0, 0.741290977639283),
        (20000.0, 1000.0, 0.95467154987681),
        (20000.0, 10000.0, 0.542260764523188),
    ] {
        let mut distribution = generator
            .estimate_fragment_from_weights(
                precursor,
                fragment,
                &[0, 1],
                AveragineComposition::PEPTIDE,
            )
            .unwrap();
        distribution.renormalize().unwrap();
        near(distribution.peaks()[0].probability, expected, 1e-6);
    }
    let mut whole = generator
        .estimate_fragment_from_weights(200.0, 200.0, &[0, 1], AveragineComposition::PEPTIDE)
        .unwrap();
    whole.renormalize().unwrap();
    let precursor = solver(Some(2)).estimate_from_peptide_weight(200.0).unwrap();
    for (whole, precursor) in whole.peaks().iter().zip(precursor.peaks()) {
        near(whole.mass, precursor.mass, 1e-12);
        near(whole.probability, precursor.probability, 1e-12);
    }
}

#[test]
fn poisson_approximation_is_normalized_and_stable() {
    // Independently evaluated binary64 recurrence from pinned
    // CoarseIsotopePatternGenerator.cpp: factor/k, then multiplication, then
    // running normalization sum. Log-space evaluation differs in its last bits.
    assert_eq!(
        Generator::approximate_intensities(1234.56789, 6)
            .unwrap()
            .into_iter()
            .map(f64::to_bits)
            .collect::<Vec<_>>(),
        [
            0x3fe01e3e8eb96f3b,
            0x3fd61c24a2e850b1,
            0x3fbe5441b0587550,
            0x3f9bbc5f2c12fd31,
            0x3f7305f07f0e5727,
            0x3f44e03da294ab0a,
        ]
    );
    let weights = Generator::approximate_intensities(1800.0, 4).unwrap();
    let normalizer = 1.0 + 1.0 + 0.5 + 1.0 / 6.0;
    for (&actual, expected) in weights.iter().zip([1.0, 1.0, 0.5, 1.0 / 6.0]) {
        near(actual, expected / normalizer, 1e-14);
    }
    let peaks = Generator::approximate_from_peptide_weight(1800.0, 4, 2).unwrap();
    near(peaks.peaks()[0].mass, 1800.0, 1e-12);
    near(peaks.peaks()[1].mass, 1800.0 + 1.00866491566 / 2.0, 1e-12);
    assert_eq!(
        Generator::approximate_intensities(0.0, 3).unwrap(),
        [1.0, 0.0, 0.0]
    );
    let enormous = Generator::approximate_intensities(f64::MAX, 20).unwrap();
    near(enormous.iter().sum(), 1.0, 1e-14);
    assert!(enormous.iter().all(|p| p.is_finite()));
    assert!(Generator::approximate_intensities(1.0, 0).is_err());
    assert!(Generator::approximate_from_peptide_weight(1.0, 3, 0).is_err());
    for mass in [20.0, 300.0, 1000.0, 2500.0] {
        let approximation = Generator::approximate_intensities(mass, 20).unwrap();
        let coarse = solver(Some(20)).estimate_from_peptide_weight(mass).unwrap();
        let kl: f64 = coarse
            .peaks()
            .iter()
            .zip(&approximation)
            .filter(|(p, _)| p.probability > 0.0)
            .map(|(p, q)| p.probability * (p.probability / q).ln())
            .sum();
        assert!(kl > 0.0 && kl < 0.05, "mass {mass}: KL {kl}");
    }
}

#[test]
fn invalid_inputs_and_computational_limits_fail_cleanly() {
    assert!(Generator::new(Some(0), CoarseMassMode::Approximate).is_err());
    assert!(Generator::new(Some(MAX_ISOTOPE_PEAKS + 1), CoarseMassMode::Approximate).is_err());
    let identity = IsotopeDistribution::default();
    for input in [
        distribution(&[(13.0, 0.5), (12.0, 0.5)]),
        distribution(&[(12.0, 0.5), (12.1, 0.5)]),
        distribution(&[(1e20, 1.0)]),
    ] {
        assert!(solver(None).convolve(&input, &identity).is_err());
    }
    let wide = distribution(&[(0.0, 0.5), (9999.0, 0.5)]);
    assert!(solver(None).convolve(&wide, &wide).is_err()); // 100 million products
    assert_eq!(
        probabilities(&solver(Some(3)).convolve(&wide, &wide).unwrap()),
        [0.25, 0.0, 0.0]
    );
    let huge_gap = distribution(&[(0.0, 0.5), (1e12, 0.5)]);
    assert!(solver(None).convolve(&huge_gap, &identity).is_err());
    assert!(solver(Some(3)).convolve(&huge_gap, &identity).is_ok());
    let explosive = distribution(&[(0.0, f64::MAX)]);
    assert!(solver(None).convolve(&explosive, &explosive).is_err());
    for mass in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            AveragineComposition::PEPTIDE
                .estimate_average_mass(mass)
                .is_err()
        );
        assert!(Generator::approximate_intensities(mass, 3).is_err());
    }
    let invalid_comp = AveragineComposition {
        carbon: 0.0,
        hydrogen: 0.0,
        nitrogen: 0.0,
        oxygen: 0.0,
        sulfur: 0.0,
        phosphorus: 0.0,
    };
    assert!(invalid_comp.estimate_mono_mass(100.0).is_err());
    assert!(
        AveragineComposition::PEPTIDE
            .estimate_mono_mass(f64::MAX)
            .is_err()
    );
    assert!(
        solver(None)
            .estimate_fragment_from_weights(10.0, 20.0, &[0], AveragineComposition::PEPTIDE)
            .is_err()
    );
}
