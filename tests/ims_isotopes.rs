use openms::chemistry::{
    IMSIsotopeDistribution as Distribution, IMSIsotopeOptions as Options, IMSIsotopePeak as Peak,
};

fn distribution(nominal: u32, values: &[(f64, f64)]) -> Distribution {
    Distribution::from_peaks(
        values
            .iter()
            .map(|&(mass, abundance)| Peak { mass, abundance })
            .collect(),
        nominal,
    )
    .unwrap()
}
fn options(size: usize) -> Options {
    Options {
        size,
        abundances_sum_error: 0.0,
    }
}

#[test]
fn constructors_accessors_and_hidden_tail_follow_source() {
    let empty = Distribution::new(12);
    assert_eq!(empty.nominal_mass(), 12);
    assert!(empty.is_empty());
    assert_eq!(empty.average_mass().unwrap(), 0.0);
    assert!(empty.mass(0).is_err());
    let single = Distribution::from_mass(15.9994).unwrap();
    assert_eq!(single.nominal_mass(), 0);
    assert_eq!(single.mass(0).unwrap(), 15.9994);
    assert_eq!(single.abundance(0).unwrap(), 1.0);
    assert_eq!(single.size(Options::default()).unwrap(), 0);
    assert!(!single.is_empty());
    let mut d = distribution(10, &[(0.125, 0.5), (0.25, 0.5), (0.5, 1.0)]);
    assert_eq!(d.size(options(2)).unwrap(), 2);
    assert_eq!(d.mass(2).unwrap(), 12.5);
    assert_eq!(d.average_mass().unwrap(), 23.1875);
    assert_eq!(d.masses(options(2)).unwrap(), vec![10.125, 11.25]);
    assert_eq!(d.abundances(options(2)).unwrap(), vec![0.5, 0.5]);
    assert_ne!(d, distribution(10, &[(0.125, 0.5), (0.25, 0.5)]));
    d.set_nominal_mass(20);
    assert_eq!(d.mass(2).unwrap(), 22.5);
    let large = distribution(
        1,
        &[
            (9_007_199_254_740_992.0, 1.0),
            (9_007_199_254_740_992.0, 1.0),
        ],
    );
    // (defect + nominal) + index: combining the integer offsets first would differ.
    assert_eq!(large.mass(1).unwrap(), 9_007_199_254_740_992.0);
}

#[test]
fn normalization_uses_all_bins_and_literal_signed_sum_condition() {
    let mut d = distribution(0, &[(0.0, 1.0), (0.0, 3.0)]);
    d.normalize(options(0)).unwrap();
    assert_eq!(d.abundances(options(2)).unwrap(), vec![0.25, 0.75]);
    let mut signed = distribution(0, &[(0.0, -1.0), (0.0, 3.0)]);
    signed.normalize(options(1)).unwrap();
    assert_eq!(signed.abundances(options(2)).unwrap(), vec![-0.5, 1.5]);
    for values in [[-2.0, 1.0], [-1.0, 1.0]] {
        let mut d = distribution(0, &[(0.0, values[0]), (0.0, values[1])]);
        let before = d.clone();
        d.normalize(options(2)).unwrap();
        assert_eq!(d, before);
    }
    let mut d = distribution(0, &[(0.0, 0.25), (0.0, 0.5)]);
    let before = d.clone();
    d.normalize(Options {
        size: usize::MAX,
        abundances_sum_error: 0.25,
    })
    .unwrap();
    assert_eq!(d, before); // strict > threshold, size unused by normalize
    d.normalize(Options {
        size: 0,
        abundances_sum_error: -1.0,
    })
    .unwrap();
    assert_eq!(d.abundance(0).unwrap(), 1.0 / 3.0);
}

#[test]
fn normalization_preserves_finite_overflow_results_and_rolls_back_invalid_results() {
    let mut d = distribution(0, &[(0.0, f64::MAX), (0.0, f64::MAX), (0.0, -1.0)]);
    d.normalize(options(3)).unwrap();
    assert_eq!(d.abundances(options(3)).unwrap(), vec![0.0, 0.0, -0.0]);
    assert!(d.abundance(2).unwrap().is_sign_negative());
    let tiny = f64::from_bits(1);
    let mut d = distribution(0, &[(0.0, tiny)]);
    let before = d.clone();
    assert!(d.normalize(options(1)).is_err());
    assert_eq!(d, before);
    assert!(
        d.normalize(Options {
            size: 1,
            abundances_sum_error: f64::NAN
        })
        .is_err()
    );
    assert_eq!(d, before);
    // Once the finite-input running sum overflows, later negative finite terms leave +Inf.
    let mut d = distribution(
        0,
        &[
            (0.0, f64::MAX),
            (0.0, f64::MAX),
            (0.0, -f64::MAX),
            (0.0, -f64::MAX),
        ],
    );
    d.normalize(options(4)).unwrap();
    assert_eq!(d.abundance(0).unwrap(), 0.0); // sum remains +Inf in this ordering
}

#[test]
fn independently_derived_dyadic_convolution_and_truncated_renormalization() {
    let lhs = distribution(10, &[(0.125, 0.5), (0.25, 0.5)]);
    let rhs = distribution(20, &[(0.5, 0.25), (0.75, 0.75)]);
    let folded = lhs.convolve(&rhs, options(3)).unwrap();
    assert_eq!(
        folded,
        distribution(30, &[(0.625, 0.125), (0.84375, 0.5), (1.0, 0.375)])
    );
    assert_eq!(
        folded.masses(options(3)).unwrap(),
        vec![30.625, 31.84375, 33.0]
    );
    let truncated = lhs.convolve(&rhs, options(2)).unwrap();
    assert_eq!(truncated, distribution(30, &[(0.625, 0.2), (0.84375, 0.8)]));
    assert_eq!(rhs.stored_len(), 2);
    let mut assigned = lhs.clone();
    assigned.convolve_assign(&rhs, options(3)).unwrap();
    assert_eq!(assigned, folded);
    assert_eq!(rhs.stored_len(), 2);
}

#[test]
fn explicit_rhs_padding_and_hidden_tail_are_atomic() {
    let mut lhs = distribution(10, &[(0.125, 0.5), (0.25, 0.5)]);
    let mut rhs = distribution(20, &[(0.5, 0.25), (0.75, 0.75)]);
    let expected = lhs.convolve(&rhs, options(3)).unwrap();
    lhs.convolve_assign_with_padded_rhs(&mut rhs, options(3))
        .unwrap();
    assert_eq!(lhs, expected);
    assert_eq!(rhs.stored_len(), 3);
    assert_eq!(rhs.peaks()[2], Peak::default());
    assert_eq!(rhs.mass(2).unwrap(), 22.0);
    let before = rhs.clone();
    lhs.convolve_assign_with_padded_rhs(&mut rhs, options(1))
        .unwrap();
    assert_eq!(rhs, before);
    lhs.set_nominal_mass(u32::MAX);
    let before = (lhs.clone(), rhs.clone());
    assert!(
        lhs.convolve_assign_with_padded_rhs(&mut rhs, options(4))
            .is_err()
    );
    assert_eq!((lhs, rhs), before);
}

#[test]
fn empty_shortcuts_precede_options_and_nominal_arithmetic() {
    let invalid = Options {
        size: usize::MAX,
        abundances_sum_error: f64::NAN,
    };
    let empty = Distribution::new(u32::MAX);
    let full = distribution(5, &[(0.0, 2.0), (0.0, 3.0)]);
    assert_eq!(full.convolve(&empty, invalid).unwrap(), full);
    assert_eq!(empty.convolve(&full, invalid).unwrap(), full);
    assert_eq!(
        empty.convolve(&Distribution::new(7), invalid).unwrap(),
        empty
    );
    let mut lhs = empty.clone();
    let mut rhs = full.clone();
    lhs.convolve_assign_with_padded_rhs(&mut rhs, options(10))
        .unwrap();
    assert_eq!(lhs, full);
    assert_eq!(rhs, full);
    let zero_size = full.convolve(&full, options(0)).unwrap();
    assert!(zero_size.is_empty());
    assert_eq!(zero_size.nominal_mass(), 10);
}

#[test]
fn source_power_zero_and_binary_empty_size_zero_rules() {
    let full = distribution(3, &[(0.125, 0.5), (0.375, 0.5)]);
    let invalid = Options {
        size: usize::MAX,
        abundances_sum_error: f64::NAN,
    };
    for power in [0, 1] {
        assert_eq!(full.pow(power, invalid).unwrap(), full);
        let mut d = full.clone();
        d.pow_assign(power, invalid).unwrap();
        assert_eq!(d, full);
    }
    assert_eq!(full.pow(2, options(0)).unwrap(), Distribution::default());
    assert_eq!(full.pow(3, options(0)).unwrap(), full);
    let empty = Distribution::new(9);
    assert_eq!(empty.pow(2, invalid).unwrap(), Distribution::default());
    assert_eq!(empty.pow(3, invalid).unwrap(), empty);
    // Once a square is empty, source collection is a no-op. Redundant clones
    // of this large retained result must not spend the remaining byte budget.
    let large = Distribution::from_peaks(vec![Peak::default(); 300_000], 0).unwrap();
    assert_eq!(large.pow(u32::MAX, options(0)).unwrap(), large);
}

#[test]
fn full_bernoulli_power_matches_independent_binomial_probabilities_and_mass_defects() {
    let base = distribution(3, &[(0.125, 0.5), (0.375, 0.5)]);
    for power in 2..=12u32 {
        let actual = base.pow(power, options(power as usize + 1)).unwrap();
        let mut coefficient = 1u64;
        for k in 0..=power {
            let probability = coefficient as f64 / 2f64.powi(power as i32);
            let defect = f64::from(power) * 0.125 + f64::from(k) * 0.25;
            assert_eq!(actual.abundance(k as usize).unwrap(), probability);
            assert_eq!(actual.peaks()[k as usize].mass, defect);
            if k < power {
                coefficient = coefficient * u64::from(power - k) / u64::from(k + 1);
            }
        }
    }
}

#[test]
fn exhaustive_small_dyadic_fold_matches_direct_species_enumeration() {
    for a in 1..=3 {
        for b in 1..=3 {
            let left: Vec<_> = (0..a)
                .map(|i| (i as f64 / 8.0, (i + 1) as f64 / 8.0))
                .collect();
            let right: Vec<_> = (0..b)
                .map(|i| (-(i as f64) / 4.0, (i + 1) as f64 / 4.0))
                .collect();
            let actual = distribution(2, &left)
                .convolve(&distribution(5, &right), options(a + b - 1))
                .unwrap();
            let mut probability = vec![0.0; a + b - 1];
            let mut moment = vec![0.0; a + b - 1];
            // Cartesian species enumeration, distinct from the destination-bin fold.
            for (i, &(mass_a, p_a)) in left.iter().enumerate() {
                for (j, &(mass_b, p_b)) in right.iter().enumerate() {
                    probability[i + j] += p_a * p_b;
                    moment[i + j] += (mass_a + mass_b) * p_a * p_b;
                }
            }
            let total: f64 = probability.iter().sum();
            for k in 0..probability.len() {
                assert_eq!(actual.peaks()[k].mass, moment[k] / probability[k]);
                assert_eq!(actual.peaks()[k].abundance, probability[k] * (1.0 / total));
            }
        }
    }
}

#[test]
fn zero_abundance_discards_nonfinite_weighted_numerator_but_other_failures_are_atomic() {
    let mut left = distribution(0, &[(f64::MAX, 0.0)]);
    let right = left.clone();
    left.convolve_assign(&right, options(1)).unwrap();
    assert_eq!(left, distribution(0, &[(0.0, 0.0)]));
    let mut left = distribution(1, &[(f64::MAX, 1.0)]);
    let before = left.clone();
    assert!(left.convolve_assign(&before, options(1)).is_err());
    assert_eq!(left, before);
    assert!(left.pow_assign(4, options(1)).is_err());
    assert_eq!(left, before);
}

#[test]
fn invalid_inputs_and_resource_limits_are_checked_before_expensive_work() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(Distribution::from_mass(value).is_err());
        assert!(
            Distribution::from_peaks(
                vec![Peak {
                    mass: 0.0,
                    abundance: value
                }],
                0
            )
            .is_err()
        );
    }
    assert!(Distribution::from_peaks(vec![Peak::default(); 1_000_001], 0).is_err());
    let mut d = Distribution::from_mass(1.0).unwrap();
    let before = d.clone();
    assert!(d.convolve_assign(&before, options(1_000_000)).is_err());
    assert_eq!(d, before);
    assert!(d.size(options(1_000_001)).is_err());
    assert!(d.mass(1).is_err());
    assert!(d.abundance(usize::MAX).is_err());
    let large = Distribution::from_peaks(vec![Peak::default(); 300_000], 0).unwrap();
    assert!(large.to_text(options(300_000)).is_err());
    let mut d = distribution(0, &[(0.0, 1.0)]);
    let before = d.clone();
    assert!(d.pow_assign(u32::MAX, options(1_000)).is_err());
    assert_eq!(d, before);
}

#[test]
fn classic_stream_text_applies_size_and_six_significant_digits() {
    let d = distribution(
        1,
        &[
            (0.0078250319, 0.999885),
            (0.01410178, 0.000115),
            (0.01604927, 0.0),
        ],
    );
    assert_eq!(
        d.to_text(options(3)).unwrap(),
        "1.00783 0.999885\n2.0141 0.000115\n3.01605 0\n"
    );
    assert_eq!(d.to_text(Options::default()).unwrap(), "");
    assert_eq!(
        distribution(0, &[(1e-20, -0.0)])
            .to_text(options(1))
            .unwrap(),
        "1e-20 -0\n"
    );
}

#[test]
fn scalar_formatting_cost_is_shared_across_the_whole_stream() {
    let d = Distribution::from_peaks(vec![Peak::default(); 50_000], 0).unwrap();
    // Destination alone is only 1.6MB; cumulative scalar formatting exceeds
    // 50M work units and must fail before entering the 50,000-row loop.
    assert!(d.to_text(options(50_000)).is_err());
    assert_eq!(d.to_text(options(2)).unwrap(), "0 0\n1 0\n");
}
