// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//
// Class-test sections of EmgGradientDescent_test.cpp that the earlier
// emg.rs / emg_reference.rs / emg_integration.rs suites do not reach:
// construction, the DefaultParamHandler defaults, and the mean constraint that
// computeMuMaxDistance feeds. Source literals come from the pinned test at
// revision bc9cc12514c768385ce121d6ca4bb710fe1983c4 (tier 3, source review).

use openms::analysis::emg::EmgGradientDescent;

// Literal cutoff trace from the class test (cutoff_pos_min / cutoff_int).
const X: [f64; 12] = [
    15.34253311,
    15.35624981,
    15.36995029,
    15.38366699,
    15.39736652,
    15.41156673,
    15.42574978,
    15.44018364,
    15.45436668,
    15.46856689,
    15.48274994,
    15.49695015,
];
const Y: [f64; 12] = [
    3.48297429,
    15.54384613,
    50.31319046,
    151.8971405,
    411.25631714,
    946.44311523,
    1642.56152344,
    2118.89526367,
    2055.13647461,
    1665.13232422,
    1275.53015137,
    1009.70056152,
];

// START_SECTION(EmgGradientDescent())
// START_SECTION(~EmgGradientDescent())
// START_SECTION(getParameters())
#[test]
fn default_construction_reproduces_the_source_parameter_defaults() {
    // The source constructor calls getDefaultParameters and writes them into
    // its Param; getParameters() then reports print_debug 0, max_gd_iter 100000
    // and compute_additional_points "true".
    let fitter = EmgGradientDescent::default();
    assert_eq!(fitter.max_iterations, 100_000);
    assert!(fitter.compute_additional_points);
    // `print_debug` has no counterpart: the port returns diagnostics in its
    // estimate instead of writing to stdout, so there is no level to default.

    // The source's destructor is `= default`; in Rust the value owns nothing
    // that needs one, so there is no observable destruction step and clippy
    // rejects an explicit `drop` of it. What is worth pinning instead is that
    // a clone is an independent, equal value.
    let copy = fitter.clone();
    assert_eq!(copy, fitter);
    assert_eq!(copy.max_iterations, fitter.max_iterations);

    // Resource ceilings have no source counterpart; they are native additions
    // and they are positive by default so that the defaults fit anything.
    assert!(copy.max_points > 0 && copy.max_evaluations > 0);
    copy.validate().unwrap();
}

// START_SECTION(double computeMuMaxDistance(const std::vector<double>& xs) const)
#[test]
fn the_fitted_mean_stays_within_thirty_five_percent_of_the_training_span() {
    // computeMuMaxDistance is private in C++ and reachable only through the
    // EmgGradientDescent_friend shim; its numeric section is mapped by the
    // private unit test in src/analysis/emg.rs. What a caller can observe is
    // the constraint it feeds: mu never leaves initial_mu +/- 0.35 * span.
    let fitter = EmgGradientDescent {
        compute_additional_points: false,
        ..Default::default()
    };
    // Independent reference: the initial mean of this trace, from the
    // derived-oracle row `cutoff_min` in tests/data/emg_training.tsv.
    let initial_mean = 15.448341685833332;
    let radius = (X[11] - X[0]) * 0.35;
    assert!((radius - 0.054045964).abs() < 1e-9);

    for iterations in [1_usize, 2, 10, 100, 100_000] {
        let estimate = EmgGradientDescent {
            max_iterations: iterations,
            ..fitter.clone()
        }
        .estimate_parameters(&X, &Y)
        .unwrap();
        assert!(
            (estimate.parameters.mu - initial_mean).abs() <= radius,
            "mu {} left the constraint after {iterations} iterations",
            estimate.parameters.mu
        );
    }
    // The source's own fitted mean for this trace, 15.4227, lies inside it.
    assert!((15.4227_f64 - initial_mean).abs() < radius);
}
