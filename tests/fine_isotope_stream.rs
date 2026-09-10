// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::fine_isotopes::MAX_FINE_CUSTOM_CELLS;
use openms::chemistry::{
    EmpiricalFormula, FineIsotopeConfiguration as Configuration, FineIsotopeIterator as Stream,
    FineIsotopePatternGenerator as Fine, FineIsotopeStop as Stop, IsotopeDistribution, IsotopePeak,
};
use std::iter::FusedIterator;

fn formula() -> EmpiricalFormula {
    EmpiricalFormula::parse("C3H4O2").unwrap()
}
fn stored(mut values: Vec<Configuration>) -> IsotopeDistribution {
    values.sort_by(|a, b| a.mass.total_cmp(&b.mass));
    IsotopeDistribution::from_peaks(
        values
            .into_iter()
            .map(|configuration| {
                let peak = configuration.to_peak().unwrap();
                IsotopePeak {
                    mass: peak.mz,
                    probability: f64::from(peak.intensity),
                }
            })
            .collect(),
    )
    .unwrap()
}

#[test]
fn formula_stream_collectors_reproduce_all_materialized_selection_modes() {
    let formula = formula();
    for threshold in [0.0, 1e-6, 0.01, 2.0] {
        let absolute = Stream::from_formula(&formula)
            .unwrap()
            .with_absolute_threshold(threshold)
            .unwrap()
            .collect::<openms::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(
            stored(absolute),
            Fine {
                stop: Stop::AbsoluteThreshold(threshold)
            }
            .run(&formula)
            .unwrap()
        );
        let relative = Stream::from_formula(&formula)
            .unwrap()
            .with_relative_threshold(threshold)
            .unwrap()
            .collect::<openms::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(
            stored(relative),
            Fine {
                stop: Stop::RelativeThreshold(threshold)
            }
            .run(&formula)
            .unwrap()
        );
    }
    for unexplained in [0.0, 0.01, 0.1] {
        let mut values = Vec::new();
        let mut sum = 0.0;
        for item in Stream::from_formula(&formula).unwrap() {
            let configuration = item.unwrap();
            sum += f64::from(configuration.to_peak().unwrap().intensity);
            values.push(configuration);
            if sum >= 1.0 - unexplained {
                break;
            }
        }
        assert_eq!(
            stored(values),
            Fine {
                stop: Stop::UnexplainedProbability(unexplained)
            }
            .run(&formula)
            .unwrap()
        );
    }
}

#[test]
fn thresholds_replace_previous_settings_and_use_original_mode_after_a_prefix() {
    let input = || Stream::from_isotopes(&[2], &[vec![10., 11.]], &[vec![0.25, 0.75]]).unwrap();
    let reopened = input()
        .with_absolute_threshold(2.0)
        .unwrap()
        .with_absolute_threshold(0.0)
        .unwrap();
    assert_eq!(
        reopened.collect::<openms::Result<Vec<_>>>().unwrap().len(),
        3
    );
    let mut partial = input();
    let mode = partial.next().unwrap().unwrap();
    assert_eq!(mode.mass, 22.);
    // Remaining maximum=.375; .8 of the ORIGINAL .5625 mode excludes it.
    let mut excluded = partial.with_relative_threshold(0.8).unwrap();
    assert!(excluded.next().is_none());
    assert!(
        excluded
            .with_absolute_threshold(0.0)
            .unwrap()
            .next()
            .is_none()
    );
    let mut partial = input();
    partial.next().unwrap().unwrap();
    let included = partial
        .with_relative_threshold(0.5)
        .unwrap()
        .collect::<openms::Result<Vec<_>>>()
        .unwrap();
    assert_eq!(included.len(), 1);
    assert_eq!(included[0].mass, 21.);
}

#[test]
fn raw_formula_charge_is_ignored_and_owned_inputs_can_be_changed_or_dropped() {
    let neutral = formula();
    let expected = Stream::from_formula(&neutral)
        .unwrap()
        .collect::<openms::Result<Vec<_>>>()
        .unwrap();
    for charge in [2, -2, i32::MIN, i32::MAX] {
        let charged = neutral.clone().with_charge(charge);
        assert_eq!(
            Stream::from_formula(&charged)
                .unwrap()
                .collect::<openms::Result<Vec<_>>>()
                .unwrap(),
            expected
        );
    }
    let mut masses = vec![vec![10., 20.]];
    let mut probabilities = vec![vec![0.2, 0.8]];
    let mut stream = Stream::from_isotopes(&[1], &masses, &probabilities).unwrap();
    masses[0][1] = 999.;
    probabilities[0][1] = 0.01;
    drop(masses);
    drop(probabilities);
    let first = stream.next().unwrap().unwrap();
    assert_eq!(first.mass, 20.);
    assert!((first.probability - 0.8).abs() < 1e-15);
    assert_eq!(stream.next().unwrap().unwrap().mass, 10.);
    drop(stream);
    assert_eq!(first.to_peak().unwrap().mz, 20.);
}

#[test]
fn stream_is_fused_after_normal_exhaustion_or_an_immediate_numeric_error() {
    fn requires_fused<T: FusedIterator>(_: &T) {}
    let mut identity = Stream::from_isotopes(&[], &[], &[]).unwrap();
    requires_fused(&identity);
    assert_eq!(
        identity.next().unwrap().unwrap(),
        Configuration {
            mass: 0.,
            probability: 1.,
            log_probability: 0.
        }
    );
    for _ in 0..3 {
        assert!(identity.next().is_none());
    }
    let mut overflow = Stream::from_isotopes(&[2], &[vec![1.]], &[vec![1e308]]).unwrap();
    assert!(overflow.next().unwrap().is_err());
    for _ in 0..3 {
        assert!(overflow.next().is_none());
    }
    assert!(
        overflow
            .with_relative_threshold(0.0)
            .unwrap()
            .next()
            .is_none()
    );
}

#[test]
fn public_configuration_conversion_checks_each_field_and_allows_underflow() {
    let valid = Configuration {
        mass: 12.,
        probability: 0.5,
        log_probability: 0.5f64.ln(),
    };
    assert_eq!(valid.to_peak().unwrap().intensity, 0.5);
    for bad in [
        Configuration { mass: -1., ..valid },
        Configuration {
            mass: f64::NAN,
            ..valid
        },
        Configuration {
            probability: -1.,
            ..valid
        },
        Configuration {
            probability: f64::INFINITY,
            ..valid
        },
        Configuration {
            probability: f64::MAX,
            ..valid
        },
        Configuration {
            log_probability: f64::NEG_INFINITY,
            ..valid
        },
    ] {
        assert!(bad.to_peak().is_err());
    }
    let underflow = Configuration {
        mass: 12.,
        probability: f64::MIN_POSITIVE,
        log_probability: f64::MIN_POSITIVE.ln(),
    };
    assert_eq!(underflow.to_peak().unwrap().intensity, 0.);
}

#[test]
fn zero_count_rows_are_validated_even_for_empty_materialized_coverage() {
    let none = Fine {
        stop: Stop::UnexplainedProbability(1.),
    };
    for (counts, masses, weights) in [
        (vec![0], vec![], vec![]),
        (vec![0], vec![vec![]], vec![vec![]]),
        (vec![0], vec![vec![1., 2.]], vec![vec![1.]]),
        (vec![0], vec![vec![f64::NAN]], vec![vec![1.]]),
        (vec![0], vec![vec![1.]], vec![vec![0.]]),
        (vec![0], vec![vec![1.]], vec![vec![f64::INFINITY]]),
    ] {
        assert!(Stream::from_isotopes(&counts, &masses, &weights).is_err());
        assert!(none.run_with_isotopes(&counts, &masses, &weights).is_err());
    }
    assert!(
        none.run_with_isotopes(&[0], &[vec![12.]], &[vec![0.5]])
            .unwrap()
            .is_empty()
    );
    let identity = Stream::from_isotopes(&[0], &[vec![12.]], &[vec![0.5]])
        .unwrap()
        .next()
        .unwrap()
        .unwrap();
    assert_eq!(identity.probability, 1.);
}

#[test]
fn custom_input_cap_is_checked_even_when_every_population_is_empty() {
    let masses = vec![vec![1.; MAX_FINE_CUSTOM_CELLS + 1]];
    let probabilities = vec![vec![1.; MAX_FINE_CUSTOM_CELLS + 1]];
    assert!(Stream::from_isotopes(&[0], &masses, &probabilities).is_err());
}

#[test]
fn large_support_can_be_streamed_as_a_small_prefix_without_materialization() {
    let formula = EmpiricalFormula::parse("Sn100").unwrap();
    assert!(
        Fine {
            stop: Stop::AbsoluteThreshold(0.)
        }
        .run(&formula)
        .is_err()
    );
    let prefix = Stream::from_formula(&formula)
        .unwrap()
        .take(3)
        .collect::<openms::Result<Vec<_>>>()
        .unwrap();
    assert_eq!(prefix.len(), 3);
    assert!(
        prefix
            .windows(2)
            .all(|pair| pair[0].log_probability >= pair[1].log_probability)
    );
}
