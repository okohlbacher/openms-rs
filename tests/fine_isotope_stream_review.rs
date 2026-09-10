// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Independent adversarial and analytical checks for the native ordered stream.
//! These are hand-derived cases, not claimed IsoSpec backend runtime goldens.

use openms::chemistry::fine_isotopes::{
    FineIsotopeConfiguration as Configuration, FineIsotopeIterator as Stream,
    FineIsotopePatternGenerator as Fine, FineIsotopeStop as Stop, MAX_FINE_CUSTOM_CELLS,
};

fn one_group(count: u32, masses: &[f64], weights: &[f64]) -> Stream {
    Stream::from_isotopes(&[count], &[masses.to_vec()], &[weights.to_vec()]).unwrap()
}
fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        ((actual - expected) / expected).abs() <= tolerance,
        "{actual:e} != {expected:e} (relative tolerance {tolerance})"
    );
}
fn assert_fused(stream: &mut Stream) {
    fn fused_trait<T: std::iter::FusedIterator>(_iterator: &T) {}
    fused_trait(stream);
    for _ in 0..4 {
        assert!(stream.next().is_none());
    }
}

#[test]
fn borrowed_custom_tables_and_returned_values_have_independent_owned_lifetimes() {
    let mut masses = vec![vec![10.0, 11.0]];
    let mut probabilities = vec![vec![0.25, 0.75]];
    let mut stream = Stream::from_isotopes(&[2], &masses, &probabilities).unwrap();
    masses[0].fill(f64::NAN);
    probabilities[0].clear();
    drop(masses);
    drop(probabilities);
    let first = stream.next().unwrap().unwrap();
    let preserved = first;
    for (mass, probability) in [(21.0, 0.375), (20.0, 0.0625)] {
        let item = stream.next().unwrap().unwrap();
        assert_eq!(item.mass, mass);
        relative(item.probability, probability, 2e-15);
    }
    assert_fused(&mut stream);
    assert_eq!(first, preserved);
    assert_eq!(first.mass, 22.0);
    relative(first.probability, 0.5625, 2e-15);
    drop(stream);
    let peak = first.to_peak().unwrap();
    assert_eq!(peak.mz, 22.0);
    assert_eq!(peak.intensity, 0.5625);
}

#[test]
fn finite_extreme_weights_are_scaled_only_for_mode_selection() {
    // The unscaled sum overflows, but every one-atom configuration is finite.
    let mut stream = one_group(1, &[1.0, 2.0], &[f64::MAX, f64::MAX]);
    let first = stream.next().unwrap().unwrap();
    let second = stream.next().unwrap().unwrap();
    for item in [first, second] {
        relative(item.probability, f64::MAX, 1e-12);
        relative(item.log_probability, f64::MAX.ln(), 1e-14);
        assert!(item.to_peak().is_err()); // raw f64 value cannot become f32.
    }
    assert_ne!(first.mass, second.mass);
    assert_fused(&mut stream);

    // Scaling this ratio underflows. The small category must still be retained
    // in the raw model and returned later with its original log probability.
    let smallest = f64::from_bits(1);
    let mut stream = one_group(1, &[1.0, 2.0], &[f64::MAX, smallest]);
    let large = stream.next().unwrap().unwrap();
    let small = stream.next().unwrap().unwrap();
    assert_eq!(large.mass, 1.0);
    assert_eq!(small.mass, 2.0);
    assert_eq!(small.probability, smallest);
    assert_eq!(small.log_probability, smallest.ln());
    assert_fused(&mut stream);
}

#[test]
fn underflow_is_a_valid_configuration_with_a_finite_log_not_exhaustion() {
    let mut stream = one_group(2000, &[10.0], &[0.5]);
    let item = stream.next().unwrap().unwrap();
    assert_eq!(item.mass, 20_000.0);
    assert_eq!(item.probability, 0.0);
    assert_eq!(item.log_probability, 2000.0 * 0.5_f64.ln());
    assert_eq!(item.to_peak().unwrap().intensity, 0.0);
    assert_fused(&mut stream);
    // Conversely, non-unit weights above one are not normalized away.
    let item = one_group(2, &[10.0], &[2.0]).next().unwrap().unwrap();
    assert_eq!(item.mass, 20.0);
    assert_eq!(item.probability, 4.0);

    // Both raw probabilities underflow, but their ratio is still ten. Relative
    // thresholding must use informative logs instead of comparing two zeros.
    let selected: Vec<_> = Stream::from_isotopes(
        &[1, 2],
        &[vec![1.0, 2.0], vec![10.0]],
        &[vec![1e-300, 1e-301], vec![1e-100]],
    )
    .unwrap()
    .with_relative_threshold(0.5)
    .unwrap()
    .collect::<openms::Result<_>>()
    .unwrap();
    assert_eq!(selected.len(), 1);
    assert_eq!(selected[0].mass, 21.0);
    assert_eq!(selected[0].probability, 0.0);
}

#[test]
fn a_late_mass_overflow_returns_one_error_then_permanent_exhaustion() {
    let mut stream = one_group(2, &[1.0, 1e308], &[0.99, 0.01]);
    let first = stream.next().unwrap().unwrap();
    let second = stream.next().unwrap().unwrap();
    assert_eq!(first.mass, 2.0);
    assert_eq!(second.mass, 1e308);
    assert!(stream.next().unwrap().is_err()); // two heavy atoms overflow mass.
    assert_fused(&mut stream);
    let mut stream = stream.with_absolute_threshold(0.0).unwrap();
    assert_fused(&mut stream); // changing a threshold never revives failure.
    assert_eq!(first.mass, 2.0);
    relative(first.probability, 0.9801, 2e-15);
    relative(second.probability, 0.0198, 2e-15);
}

#[test]
fn changed_thresholds_filter_remaining_states_using_the_original_mode() {
    let mut stream = one_group(2, &[10.0, 11.0], &[0.25, 0.75]);
    assert_eq!(stream.next().unwrap().unwrap().mass, 22.0);
    assert_eq!(stream.next().unwrap().unwrap().mass, 21.0);
    // Original mode .5625 * .15 = .084375 rejects the last .0625. Using the
    // last returned .375 (or the remaining .0625) as mode would wrongly retain it.
    let mut stream = stream.with_relative_threshold(0.15).unwrap();
    assert_fused(&mut stream);

    let mut stream = one_group(2, &[10.0, 11.0], &[0.25, 0.75])
        .with_absolute_threshold(0.5)
        .unwrap();
    assert_eq!(stream.next().unwrap().unwrap().mass, 22.0);
    let remaining: Vec<_> = stream
        .with_absolute_threshold(0.0)
        .unwrap()
        .map(|item| item.unwrap().mass)
        .collect();
    assert_eq!(remaining, [21.0, 20.0]);

    let mut stream = one_group(2, &[10.0, 11.0], &[0.25, 0.75])
        .with_absolute_threshold(0.5)
        .unwrap();
    stream.next().unwrap().unwrap();
    assert_fused(&mut stream);
    let mut stream = stream.with_absolute_threshold(0.0).unwrap();
    assert_fused(&mut stream); // discarded suffix cannot be restarted.
}

#[test]
fn threshold_equality_is_inclusive_and_invalid_cutoffs_are_checked() {
    let items: Vec<_> = one_group(1, &[10.0, 11.0], &[0.75, 0.25])
        .with_absolute_threshold(0.25)
        .unwrap()
        .collect::<openms::Result<_>>()
        .unwrap();
    assert_eq!(items.len(), 2);
    assert_eq!(items[1].probability, 0.25);
    let items: Vec<_> = one_group(1, &[10.0, 11.0], &[0.75, 0.25])
        .with_relative_threshold(1.0)
        .unwrap()
        .collect::<openms::Result<_>>()
        .unwrap();
    assert_eq!(items.len(), 1);
    assert_eq!(items[0].mass, 10.0);
    // A cutoff obtained from the raw accessor must also retain equality; the
    // round trip ln(exp(log p)) is not always bit-identical to the stored log p.
    let configuration = one_group(2, &[10.0], &[0.63]).next().unwrap().unwrap();
    for (threshold, expected) in [
        (f64::from_bits(configuration.probability.to_bits() - 1), 1),
        (configuration.probability, 1),
        (f64::from_bits(configuration.probability.to_bits() + 1), 0),
    ] {
        let items: Vec<_> = one_group(2, &[10.0], &[0.63])
            .with_absolute_threshold(threshold)
            .unwrap()
            .collect::<openms::Result<_>>()
            .unwrap();
        assert_eq!(items.len(), expected, "absolute cutoff {threshold}");
    }
    // Preserve the existing materialized collector's log-space cutoff. Its
    // historical formula counts must not be changed to fix the raw stream API.
    // Whether this log/exp round trip crosses the cutoff can vary by platform.
    let materializer = Fine {
        stop: Stop::AbsoluteThreshold(configuration.probability),
    };
    let expected = usize::from(2.0 * 0.63_f64.ln() >= configuration.probability.ln());
    assert_eq!(
        materializer
            .run_with_isotopes(&[2], &[vec![10.0]], &[vec![0.63]])
            .unwrap()
            .len(),
        expected
    );

    let prefix: Vec<_> = one_group(2, &[10.0, 11.0], &[0.25, 0.75])
        .take(2)
        .collect::<openms::Result<_>>()
        .unwrap();
    let ratio = prefix[1].probability / prefix[0].probability;
    for (threshold, expected) in [(ratio, 2), (f64::from_bits(ratio.to_bits() + 1), 1)] {
        let items: Vec<_> = one_group(2, &[10.0, 11.0], &[0.25, 0.75])
            .with_relative_threshold(threshold)
            .unwrap()
            .collect::<openms::Result<_>>()
            .unwrap();
        assert_eq!(items.len(), expected, "relative cutoff {threshold}");
    }
    // Adding ln(next_up(1)) to a very negative mode log rounds back to that
    // same log. A relative cutoff strictly above one still excludes the mode.
    for (threshold, expected) in [(1.0, 1), (f64::from_bits(1.0_f64.to_bits() + 1), 0)] {
        let items: Vec<_> = one_group(1_000_000, &[1.0], &[1e-300])
            .with_relative_threshold(threshold)
            .unwrap()
            .collect::<openms::Result<_>>()
            .unwrap();
        assert_eq!(items.len(), expected);
    }
    for value in [f64::NAN, f64::INFINITY, -0.1] {
        assert!(
            one_group(1, &[10.0], &[1.0])
                .with_absolute_threshold(value)
                .is_err()
        );
        assert!(
            one_group(1, &[10.0], &[1.0])
                .with_relative_threshold(value)
                .is_err()
        );
    }
}

#[test]
fn zero_count_and_zero_coverage_still_validate_supplied_custom_data() {
    let none = Fine {
        stop: Stop::UnexplainedProbability(1.0),
    };
    for bad in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        let masses = [vec![10.0]];
        let weights = [vec![bad]];
        assert!(Stream::from_isotopes(&[0], &masses, &weights).is_err());
        assert!(none.run_with_isotopes(&[0], &masses, &weights).is_err());
    }
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        let masses = [vec![bad]];
        let weights = [vec![1.0]];
        assert!(Stream::from_isotopes(&[0], &masses, &weights).is_err());
        assert!(none.run_with_isotopes(&[0], &masses, &weights).is_err());
    }
    for (counts, masses, weights) in [
        (vec![0], vec![], vec![]),
        (vec![0], vec![vec![1.0]], vec![]),
        (vec![0], vec![vec![]], vec![vec![]]),
        (vec![0], vec![vec![1.0, 2.0]], vec![vec![1.0]]),
    ] {
        assert!(Stream::from_isotopes(&counts, &masses, &weights).is_err());
        assert!(none.run_with_isotopes(&counts, &masses, &weights).is_err());
    }
    let mut stream = one_group(0, &[f64::MAX], &[f64::MAX]);
    assert_eq!(
        stream.next().unwrap().unwrap(),
        Configuration {
            mass: 0.0,
            probability: 1.0,
            log_probability: 0.0,
        }
    );
    assert_fused(&mut stream);
    assert!(
        none.run_with_isotopes(&[0], &[vec![10.0]], &[vec![1.0]])
            .unwrap()
            .is_empty()
    );
}

#[test]
fn custom_cell_limit_is_cumulative_even_for_zero_count_rows() {
    // Each row is individually below the cap. Their combined cells must be
    // checked before skipping zero populations or copying custom model tables.
    let row_length = MAX_FINE_CUSTOM_CELLS / 2 + 1;
    let masses = vec![vec![10.0; row_length], vec![11.0; row_length]];
    let weights = vec![vec![0.5; row_length], vec![0.5; row_length]];
    let result = Stream::from_isotopes(&[0, 0], &masses, &weights);
    let error = match result {
        Err(error) => error,
        Ok(_) => panic!("over-budget zero-count input was accepted"),
    };
    assert!(error.to_string().contains("input cell limit"));
    let none = Fine {
        stop: Stop::UnexplainedProbability(1.0),
    };
    let error = none
        .run_with_isotopes(&[0, 0], &masses, &weights)
        .unwrap_err();
    assert!(error.to_string().contains("input cell limit"));
}

#[test]
fn high_dimension_coverage_stops_before_unused_expansion_and_error_fuses() {
    let masses = vec![(1..=10_000).map(f64::from).collect::<Vec<_>>()];
    let mut weights = vec![vec![1e-12; 10_000]];
    weights[0][0] = 1.0;
    let fine = Fine {
        stop: Stop::UnexplainedProbability(0.01),
    };
    let result = fine.run_with_isotopes(&[1], &masses, &weights).unwrap();
    assert_eq!(result.len(), 1);
    assert_eq!(result.peaks()[0].mass, 1.0);
    assert_eq!(result.peaks()[0].probability, 1.0);

    let mut stream = Stream::from_isotopes(&[1], &masses, &weights).unwrap();
    assert_eq!(stream.next().unwrap().unwrap().mass, 1.0);
    let error = stream.next().unwrap().unwrap_err();
    assert!(error.to_string().contains("limit"));
    assert_fused(&mut stream);

    let mut stream = Stream::from_isotopes(&[1], &masses, &weights).unwrap();
    stream.next().unwrap().unwrap();
    let mut stream = stream.with_absolute_threshold(2.0).unwrap();
    assert_fused(&mut stream); // no excluded suffix expansion or limit error.
}

#[test]
fn ordered_stream_can_drain_more_than_the_materialized_peak_cap() {
    // Two 320-atom binary components have exactly 321² configurations. This
    // exceeds the materialized100,000 cap but fits streaming state/work bounds.
    let masses = [vec![10.0, 11.0], vec![1000.0, 1001.0]];
    let weights = [vec![0.5, 0.5], vec![0.5, 0.5]];
    let mut stream = Stream::from_isotopes(&[320, 320], &masses, &weights).unwrap();
    let mut count = 0;
    let mut last_log = f64::INFINITY;
    for item in stream.by_ref() {
        let item = item.unwrap();
        assert!(item.log_probability <= last_log + 1e-10);
        last_log = item.log_probability;
        count += 1;
    }
    assert_eq!(count, 103_041);
    assert_fused(&mut stream);
    let full = Fine {
        stop: Stop::AbsoluteThreshold(0.0),
    };
    let error = full
        .run_with_isotopes(&[320, 320], &masses, &weights)
        .unwrap_err();
    assert!(error.to_string().contains("peak limit"));
}
