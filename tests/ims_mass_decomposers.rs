use openms::chemistry::{
    IMSIntegerMassDecomposer as Integer, IMSMassDecomposer, IMSRealMassDecomposer as Real,
    IMSWeights,
};
use std::collections::BTreeMap;

// Independent exhaustive Cartesian traversal, with no residue table/sentinel.
fn brute(weights: &[u64], mass: u64) -> Vec<Vec<u32>> {
    fn visit(weights: &[u64], mass: u64, prefix: &mut Vec<u32>, output: &mut Vec<Vec<u32>>) {
        if weights.is_empty() {
            if mass == 0 {
                output.push(prefix.clone());
            }
            return;
        }
        for n in 0..=mass / weights[0] {
            prefix.push(n as u32);
            visit(&weights[1..], mass - n * weights[0], prefix, output);
            prefix.pop();
        }
    }
    let mut output = Vec::new();
    visit(weights, mass, &mut Vec::new(), &mut output);
    output
}
fn sum(weights: &[u64], counts: &[u32]) -> u64 {
    weights
        .iter()
        .zip(counts)
        .map(|(&w, &n)| w * u64::from(n))
        .sum()
}
fn real(masses: &[f64], precision: f64) -> Real {
    Real::new(&IMSWeights::from_masses(masses, precision).unwrap()).unwrap()
}

#[test]
fn integer_queries_match_independent_exhaustive_small_alphabets() {
    for weights in [
        vec![3],
        vec![2, 3],
        vec![6, 9, 20],
        vec![4, 5, 6],
        vec![10, 6, 15],
        vec![15, 21, 10, 13],
        vec![7, 4, 9],
        vec![6, 6, 6],
        vec![8, 12],
    ] {
        let solver = Integer::from_integer_weights(&weights).unwrap();
        let generic: &dyn IMSMassDecomposer = &solver;
        for mass in 0..=100 {
            let mut expected = brute(&weights, mass);
            let mut actual = generic.decompositions(mass).unwrap();
            expected.sort();
            actual.sort();
            assert_eq!(actual, expected, "all {weights:?} {mass}");
            assert_eq!(
                generic.number_of_decompositions(mass).unwrap() as usize,
                expected.len()
            );
            assert_eq!(
                generic.exists(mass).unwrap(),
                !expected.is_empty(),
                "exists {weights:?} {mass}"
            );
            // Source gcd-5 cache loop writes witness[3]=(2,0), see the
            // dedicated recurrence regression. All/count/exists above remain
            // checked against the independent mathematical oracle.
            if weights == [10, 6, 15] && mass >= 33 && mass % 10 == 3 {
                assert!(generic.decomposition(mass).is_err());
                continue;
            }
            match generic
                .decomposition(mass)
                .unwrap_or_else(|e| panic!("one {weights:?} {mass}: {e}"))
            {
                Some(value) => {
                    assert_eq!(value.len(), weights.len());
                    assert_eq!(
                        sum(&weights, &value),
                        mass,
                        "one {weights:?} {mass}: {value:?}"
                    );
                }
                None => assert!(expected.is_empty()),
            }
        }
    }
}

#[test]
fn source_loop_order_and_witness_are_distinct_from_first_enumeration() {
    let solver = Integer::from_integer_weights(&[6, 9, 20]).unwrap();
    assert_eq!(
        solver.decompositions(60).unwrap(),
        [
            vec![10, 0, 0],
            vec![7, 2, 0],
            vec![4, 4, 0],
            vec![1, 6, 0],
            vec![0, 0, 3]
        ]
    );
    let solver = Integer::from_integer_weights(&[4, 5, 6]).unwrap();
    assert_eq!(
        solver.decompositions(10).unwrap(),
        [vec![0, 2, 0], vec![1, 0, 1]]
    );
    assert_eq!(solver.decomposition(10).unwrap(), Some(vec![1, 0, 1]));
    // No implicit sorting or GCD division: counts remain in supplied index space.
    let permuted = Integer::from_integer_weights(&[6, 4, 5]).unwrap();
    assert_eq!(permuted.decomposition(10).unwrap(), Some(vec![1, 1, 0]));
    let mut values = permuted.decompositions(10).unwrap();
    values.sort();
    assert_eq!(values, [vec![0, 0, 2], vec![1, 1, 0]]);
}

#[test]
fn finite_source_sentinel_quirk_is_not_hidden_by_brute_force_oracle() {
    // These are deductions from source table recurrences, NOT executed C++
    // class-test goldens. Those source numerical test sections are TODOs.
    let weights = [6, 9, 2];
    let solver = Integer::from_integer_weights(&weights).unwrap();
    assert_eq!(brute(&weights, 13), [vec![0, 1, 2]]);
    assert_eq!(solver.decompositions(13).unwrap(), [vec![0, 1, 2]]);
    assert!(!solver.exists(13).unwrap());
    assert_eq!(solver.decomposition(13).unwrap(), None);
    assert_eq!(solver.number_of_decompositions(13).unwrap(), 1);
    let weights = [6, 9, 2, 2];
    let solver = Integer::from_integer_weights(&weights).unwrap();
    assert_eq!(
        brute(&weights, 13),
        [vec![0, 1, 0, 2], vec![0, 1, 1, 1], vec![0, 1, 2, 0]]
    );
    assert_eq!(
        solver.decompositions(13).unwrap(),
        [vec![0, 1, 1, 1], vec![0, 1, 0, 2]]
    );
    assert_eq!(solver.number_of_decompositions(13).unwrap(), 2);
}

#[test]
fn empty_singleton_zero_duplicate_and_integer_width_boundaries() {
    assert!(Integer::new(&IMSWeights::new()).is_err());
    assert!(Integer::from_integer_weights(&[]).is_err());
    for weights in [vec![0], vec![2, 0], vec![0, 2], vec![1; 129]] {
        assert!(Integer::from_integer_weights(&weights).is_err());
    }
    let singleton = Integer::from_integer_weights(&[u64::MAX]).unwrap();
    assert!(singleton.exists(u64::MAX).unwrap());
    assert_eq!(singleton.decomposition(u64::MAX).unwrap(), Some(vec![1]));
    assert_eq!(singleton.decompositions(0).unwrap(), [vec![0]]);
    assert_eq!(singleton.number_of_decompositions(1).unwrap(), 0);
    let ones = Integer::from_integer_weights(&[1]).unwrap();
    assert_eq!(
        ones.decomposition(u64::from(u32::MAX)).unwrap(),
        Some(vec![u32::MAX])
    );
    assert!(ones.decomposition(u64::from(u32::MAX) + 1).is_err());
    assert!(ones.decompositions(u64::from(u32::MAX) + 1).is_err());
    let deepest = Integer::from_integer_weights(&[1; 128]).unwrap();
    assert_eq!(deepest.number_of_decompositions(1).unwrap(), 128);
    let duplicate = Integer::from_integer_weights(&[6, 6, 6]).unwrap();
    assert_eq!(duplicate.decomposition(0).unwrap(), Some(vec![0, 0, 0]));
    assert_eq!(
        duplicate.decompositions(6).unwrap(),
        [vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]]
    );
}

#[test]
fn constructor_and_query_caps_are_checked_and_count_does_not_materialize() {
    assert!(Integer::from_integer_weights(&[2_000_001, 2_000_001]).is_err());
    assert!(Integer::from_integer_weights(&[2, u64::MAX]).is_err());
    let solver = Integer::from_integer_weights(&[1, 1]).unwrap();
    assert!(solver.decompositions(100_000).is_err()); // 100,001 outputs
    assert_eq!(solver.number_of_decompositions(100_000).unwrap(), 100_001);
    let real = real(&[1.0, 1.0], 1.0);
    assert!(real.decompositions(100_000.0, 1.0).is_err());
    // Source integer window 99,999..100,001 contains n+1 count pairs at each n.
    assert_eq!(
        real.number_of_decompositions(100_000.0, 1.0).unwrap(),
        200_001
    );
    // A failed materialization does not alter the immutable table.
    assert_eq!(
        solver.decompositions(2).unwrap(),
        [vec![2, 0], vec![1, 1], vec![0, 2]]
    );
    let values = solver.clone();
    drop(solver);
    assert_eq!(values.number_of_decompositions(2).unwrap(), 3);
}

#[test]
fn real_rounding_gcd_bounds_and_negative_error_follow_source_order() {
    let masses = [2.25, 3.75];
    let rounded = real(&masses, 0.5);
    // Rounded weights are [5,8]. The source exclusive endpoint can exclude
    // even an exact physical target with a small tolerance (start=end=5).
    assert!(rounded.decompositions(2.25, 0.1).unwrap().is_empty());
    assert_eq!(rounded.decompositions(2.25, 0.5).unwrap(), [vec![1, 0]]);
    assert_eq!(rounded.decompositions(3.75, 0.5).unwrap(), [vec![0, 1]]);
    let solver = real(&[3.0, 5.0, 8.0], 1.0);
    assert!(solver.decompositions(8.0, 0.0).unwrap().is_empty());
    assert_eq!(
        solver.decompositions(8.0, 1.0).unwrap(),
        [vec![1, 1, 0], vec![0, 0, 1]]
    );
    assert_eq!(solver.number_of_decompositions(8.0, 1.0).unwrap(), 2);
    assert!(solver.decompositions(8.0, -1.0).unwrap().is_empty());
    assert_eq!(solver.number_of_decompositions(8.0, -1.0).unwrap(), 0);
    let mut weights = IMSWeights::from_masses(&[3.0, 5.0, 8.0], 0.1).unwrap();
    let before = Real::new(&weights).unwrap();
    assert!(weights.divide_by_gcd().unwrap());
    let after = Real::new(&weights).unwrap();
    assert_eq!(weights.precision(), Some(1.0));
    assert_eq!(
        before.decompositions(8.0, 1.0).unwrap(),
        after.decompositions(8.0, 1.0).unwrap()
    );
}

#[test]
fn real_candidates_and_constraints_match_independent_dyadic_brute_force() {
    let masses = [1.5, 2.5, 4.0];
    let solver = real(&masses, 0.5);
    let integer_weights = [3, 5, 8];
    for target_quarters in 4..=100 {
        let target = f64::from(target_quarters) / 4.0;
        let tolerance = 0.375;
        // Independent enumerate all count tuples through the maximum physical
        // mass, then apply exact dyadic physical mass and source integer range.
        let start = ((target - tolerance) / 0.5).ceil() as u64;
        let end = ((target + tolerance) / 0.5).floor() as u64;
        let mut expected = Vec::new();
        for integer in start..end {
            for counts in brute(&integer_weights, integer) {
                let mass: f64 = masses
                    .iter()
                    .zip(&counts)
                    .map(|(m, n)| m * f64::from(*n))
                    .sum();
                if (mass - target).abs() <= tolerance {
                    expected.push(counts);
                }
            }
        }
        expected.sort();
        let mut actual = solver.decompositions(target, tolerance).unwrap();
        actual.sort();
        assert_eq!(actual, expected);
        assert_eq!(
            solver.number_of_decompositions(target, tolerance).unwrap() as usize,
            expected.len()
        );
        let constraints = BTreeMap::from([(0, (1, 2)), (2, (0, 1))]);
        expected.retain(|c| c[0] >= 1 && c[0] <= 2 && c[2] <= 1);
        let mut actual = solver
            .decompositions_with_constraints(target, tolerance, &constraints)
            .unwrap();
        actual.sort();
        assert_eq!(actual, expected);
    }
}

#[test]
fn real_zero_count_endpoint_distinction_and_checked_invalid_inputs() {
    let solver = real(&[2.0, 3.0], 1.0);
    assert!(solver.decompositions(0.0, 0.0).unwrap().is_empty());
    assert_eq!(solver.number_of_decompositions(0.0, 0.0).unwrap(), 0);
    // ceil of a small negative fraction is signed zero, a valid unsigned endpoint.
    assert!(solver.decompositions(0.0, 0.25).unwrap().is_empty());
    assert_eq!(solver.decompositions(0.5, 0.5).unwrap(), [vec![0, 0]]);
    assert_eq!(solver.number_of_decompositions(0.5, 0.5).unwrap(), 0);
    // Source enumeration would cast a negative endpoint outside unsigned range;
    // its count overload skips that calculation and uses one instead.
    assert!(solver.decompositions(0.0, 1.0).is_err());
    assert_eq!(solver.number_of_decompositions(0.0, 1.0).unwrap(), 0);
    for (mass, error) in [
        (f64::NAN, 1.0),
        (1.0, f64::INFINITY),
        (f64::MAX, 1.0),
        (-3.0, 0.0),
    ] {
        assert!(solver.decompositions(mass, error).is_err());
    }
    assert!(solver.decompositions(10.0, 100_000_000.0).is_err());
    assert!(
        solver
            .decompositions_with_constraints(6.0, 1.0, &BTreeMap::from([(2, (0, 1))]))
            .is_err()
    );
    assert!(
        solver
            .decompositions_with_constraints(6.0, 1.0, &BTreeMap::from([(0, (2, 1))]))
            .unwrap()
            .is_empty()
    );
    assert!(Real::new(&IMSWeights::new()).is_err());
    assert!(Real::new(&IMSWeights::from_masses(&[0.0, 1.0], 1.0).unwrap()).is_err());
}

#[test]
fn owned_weights_and_finite_signed_real_fields_remain_independent() {
    let mut weights = IMSWeights::from_masses(&[2.0, 3.0], 1.0).unwrap();
    let solver = Real::new(&weights).unwrap();
    let expected = solver.decompositions(6.0, 1.0).unwrap();
    weights.swap(0, 1).unwrap();
    weights.set_precision(0.5).unwrap();
    drop(weights);
    assert_eq!(solver.decompositions(6.0, 1.0).unwrap(), expected);
    // IMSWeights accepts signed mass/precision pairs yielding positive integers.
    // The real source range can be empty for positive tolerance in this case.
    let negative = real(&[-2.0, -3.0], -1.0);
    assert!(negative.decompositions(-6.0, 1.0).unwrap().is_empty());
    assert!(negative.decompositions(-6.0, -1.0).unwrap().is_empty());
}

#[test]
fn source_class_constructor_uses_natural19_precision_point_zero_one_and_gcd() {
    use openms::chemistry::{AASequence, PeptideFragmentType};
    // The upstream Integer/Real class tests assert successful construction from
    // these 19 sorted residue symbols at 0.01 precision after divideByGCD.
    let masses: Vec<_> = "ACDEFGHKLMNPQRSTVWY"
        .chars()
        .map(|symbol| {
            AASequence::parse(&symbol.to_string())
                .unwrap()
                .mono_mass_for(PeptideFragmentType::Internal, 0)
                .unwrap()
        })
        .collect();
    let mut weights = IMSWeights::from_masses(&masses, 0.01).unwrap();
    weights.divide_by_gcd().unwrap();
    let integer = Integer::new(&weights).unwrap();
    let real = Real::new(&weights).unwrap();
    assert_eq!(weights.len(), 19);
    assert_eq!(integer.decompositions(0).unwrap(), [vec![0; 19]]);
    assert_eq!(real.number_of_decompositions(0.0, 0.0).unwrap(), 0);
}
