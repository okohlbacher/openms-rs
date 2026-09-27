// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The transformation models' `Param` surface, the B-spline model, Akima
//! interpolation, LOWESS automatic span selection and
//! `TransformationModelDefaults`, against the literals of the upstream class
//! tests at core pin `bc9cc12` (evidence tier 3, source review):
//! `TransformationModel_test.cpp`, `TransformationModelBSpline_test.cpp`,
//! `TransformationModelInterpolated_test.cpp`,
//! `TransformationModelLowess_test.cpp`, `TransformationModelDefaults_test.cpp`
//! and `TransformationDescription_test.cpp`. Hashes are in
//! `tests/data/transformation_models_provenance.json`.

use openms::analysis::transformations::*;
use openms::param::{Param, ParamValue};
use std::cmp::Ordering;

/// ClassTest `isRealSimilar` (`src/testframework/source/CONCEPT/ClassTest.cpp:364`
/// at the pin): absolute difference at most 1e-5, or ratio at most 1 + 1e-5.
fn real_similar(actual: f64, expected: f64) -> bool {
    if actual.is_nan() || expected.is_nan() {
        return false;
    }
    if actual == expected {
        return true;
    }
    let small = (actual - expected).abs() <= 1e-5;
    if actual == 0.0 || expected == 0.0 {
        return small;
    }
    let mut ratio = actual / expected;
    if ratio < 0.0 {
        return small;
    }
    if ratio < 1.0 {
        ratio = 1.0 / ratio;
    }
    ratio <= 1.0 + 1e-5 || small
}

fn similar(actual: f64, expected: f64, what: &str) {
    assert!(
        real_similar(actual, expected),
        "{what}: got {actual:.17}, expected {expected}"
    );
}

fn table(text: &str) -> Vec<Vec<f64>> {
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split_whitespace().map(|v| v.parse().unwrap()).collect())
        .collect()
}

/// The noisy sine of the B-spline and LOWESS class tests; the two tests
/// declare identical `x[]` and `y[]` literal arrays.
fn sine() -> Vec<DataPoint> {
    table(include_str!("data/transformations_lowess_input.tsv"))
        .iter()
        .map(|r| DataPoint::new(r[0], r[1]))
        .collect()
}

/// The eleven points declared at the top of the B-spline and LOWESS tests.
fn eleven() -> Vec<DataPoint> {
    [
        (1.2, 5.2),
        (3.2, 7.3),
        (2.2, 6.25),
        (2.2, 3.1),
        (2.2, 7.25),
        (3.0, 8.5),
        (3.1, 4.7),
        (1.7, 6.0),
        (2.9, 4.7),
        (4.2, 5.0),
        (3.7, -2.4),
    ]
    .into_iter()
    .map(DataPoint::from)
    .collect()
}

fn param(values: &[(&str, ParamValue)]) -> Param {
    let mut p = Param::new();
    for (key, value) in values {
        p.set_value(key, value.clone(), "", &[]).unwrap();
    }
    p
}

fn s(value: &str) -> ParamValue {
    ParamValue::String(value.into())
}

/// `TransformationModelBSpline_test.cpp`, section `evaluate`: `pred[]`, the
/// B-spline with `extrapolate` `b_spline` at `v = -4, -3.8, ...` while `v < 4.1`.
const BSPLINE_PRED: [f64; 41] = [
    0.846137, 0.689856, 0.5094, 0.31183, 0.10421, -0.106399, -0.312921, -0.508271, -0.685362,
    -0.837111, -0.95643, -1.03623, -1.06944, -1.05016, -0.981868, -0.872412, -0.729666, -0.561505,
    -0.375803, -0.180434, 0.016728, 0.20827, 0.38867, 0.55289, 0.695895, 0.812645, 0.898104,
    0.947234, 0.955013, 0.919484, 0.845545, 0.739022, 0.60574, 0.451526, 0.282206, 0.103606,
    -0.0784482, -0.258084, -0.429374, -0.586387, -0.723191,
];

#[test]
fn bspline_sine_predictions_and_all_extrapolations() {
    let params = param(&[
        ("wavelength", ParamValue::Float(0.0)),
        ("num_nodes", ParamValue::Integer(5)),
        ("extrapolate", s("b_spline")),
    ]);
    let ModelConfig::BSpline(options) = ModelConfig::from_param("b_spline", &params).unwrap()
    else {
        panic!("b_spline parameters resolve to another model");
    };
    let model = BSplineModel::fit(&sine(), options).unwrap();
    // The class test's loop variable is accumulated, not multiplied.
    let mut v = -4.0;
    let mut index = 0;
    while v < 4.1 {
        similar(model.apply(v).unwrap(), BSPLINE_PRED[index], "b_spline");
        v += 0.2;
        index += 1;
    }
    assert_eq!(index, BSPLINE_PRED.len());

    for (name, left, right) in [
        ("linear", 0.947997, -0.807806),
        ("constant", 0.0150243, -0.00429613),
        ("global_linear", -0.959617, 1.10039),
    ] {
        let mut params = params.checked_clone().unwrap();
        params.set_value("extrapolate", s(name), "", &[]).unwrap();
        let mut description = TransformationDescription::new(sine()).unwrap();
        description
            .fit_model_with_param("b_spline", &params)
            .unwrap();
        similar(description.apply(-4.0).unwrap(), left, name);
        similar(description.apply(4.0).unwrap(), right, name);
    }
}

#[test]
fn bspline_needs_two_points_and_reports_its_parameters() {
    // Section `TransformationModelBSpline(const DataPoints&, const Param&)`:
    // the empty data set is an IllegalArgument, the eleven points fit.
    assert!(BSplineModel::fit(&[], BSplineOptions::default()).is_err());
    assert!(BSplineModel::fit(&eleven(), BSplineOptions::default()).is_ok());
    // Section `getParameters`.
    let mut description = TransformationDescription::new(eleven()).unwrap();
    let input = param(&[("num_nodes", ParamValue::Integer(5))]);
    description
        .fit_model_with_param("b_spline", &input)
        .unwrap();
    assert_eq!(
        description.model_parameters().value("num_nodes").unwrap(),
        &ParamValue::Integer(5)
    );
}

#[test]
fn description_fits_models_by_name_with_the_upstream_values() {
    // TransformationDescription_test.cpp, section `fitModel`.
    let data: Vec<_> = [(0.0, 1.0), (0.25, 1.5), (0.5, 2.0), (1.0, 3.0)]
        .into_iter()
        .map(DataPoint::from)
        .collect();
    let nonlinear: Vec<_> = [(0.0, 1.0), (0.25, 1.0625), (0.5, 1.25), (1.0, 2.0)]
        .into_iter()
        .map(DataPoint::from)
        .collect();
    let empty = Param::new();
    let mut td = TransformationDescription::new(data).unwrap();
    td.fit_model_with_param("linear", &empty).unwrap();
    assert_eq!(td.model_type(), "linear");
    for (x, y) in [(0.0, 1.0), (0.5, 2.0), (1.0, 3.0)] {
        similar(td.apply(x).unwrap(), y, "linear");
    }
    td.fit_model_with_param("b_spline", &empty).unwrap();
    assert_eq!(td.model_type(), "b_spline");
    for (x, y) in [(0.0, 1.064201730), (0.5, 1.957836652), (1.0, 2.927541901)] {
        similar(td.apply(x).unwrap(), y, "b_spline");
    }
    td.fit_model_with_param("lowess", &empty).unwrap();
    assert_eq!(td.model_type(), "lowess");
    for (x, y) in [(0.0, 1.0), (0.5, 2.0), (1.0, 3.0)] {
        similar(td.apply(x).unwrap(), y, "lowess");
    }
    td.fit_model_with_param("identity", &empty).unwrap();
    assert_eq!(td.model_type(), "identity");
    for x in [0.0, 0.5, 1.0] {
        similar(td.apply(x).unwrap(), x, "identity");
    }
    // Can't fit a different model to an "identity" transformation.
    td.fit_model_with_param("linear", &empty).unwrap();
    assert_eq!(td.model_type(), "identity");

    let mut nl = TransformationDescription::new(nonlinear.clone()).unwrap();
    nl.fit_model_with_param("b_spline", &empty).unwrap();
    for (x, y) in [
        (0.0, 1.01084556836969),
        (0.5, 1.26289804387079),
        (0.75, 1.53463130131214),
        (1.0, 1.94984504419826),
    ] {
        similar(nl.apply(x).unwrap(), y, "nonlinear b_spline");
    }
    nl.fit_model_with_param("lowess", &empty).unwrap();
    for (x, y) in [
        (0.0, 1.0),
        (0.5, 1.25),
        (0.75, 1.58423913043478),
        (1.0, 2.0),
    ] {
        similar(nl.apply(x).unwrap(), y, "nonlinear lowess");
    }

    // Section `[EXTRA] fitModel`: one Param carries a LOWESS and a B-spline
    // setting; each model reads its own.
    let both = param(&[
        ("interpolation_type", s("linear")),
        ("extrapolate", s("b_spline")),
    ]);
    let mut nl = TransformationDescription::new(nonlinear).unwrap();
    nl.fit_model_with_param("b_spline", &both).unwrap();
    for (x, y) in [
        (0.0, 1.01084556836969),
        (0.5, 1.26289804387079),
        (0.75, 1.53463130131214),
        (1.0, 1.94984504419826),
        (2.0, 1.328125),
    ] {
        similar(nl.apply(x).unwrap(), y, "b_spline extrapolation");
    }
    nl.fit_model_with_param("lowess", &both).unwrap();
    for (x, y) in [
        (0.0, 1.0),
        (0.5, 1.25),
        (0.75, 1.625),
        (1.0, 2.0),
        (2.0, 3.5),
    ] {
        similar(nl.apply(x).unwrap(), y, "lowess linear interpolation");
    }

    // An unknown model name is refused and leaves the model alone.
    let mut td = TransformationDescription::new(eleven()).unwrap();
    td.fit_model_with_param("linear", &empty).unwrap();
    assert!(td.fit_model_with_param("mumble", &empty).is_err());
    assert_eq!(td.model_type(), "linear");
}

#[test]
fn model_parameters_are_what_the_source_model_keeps() {
    // Section `getModelParameters`.
    let mut td = TransformationDescription::default();
    assert!(td.model_parameters().is_empty());
    let input = param(&[
        ("slope", ParamValue::Float(2.5)),
        ("intercept", ParamValue::Float(-100.0)),
        ("x_weight", s("x")),
        ("y_weight", s("y")),
        ("x_datum_min", ParamValue::Float(1e-15)),
        ("y_datum_min", ParamValue::Float(1e-15)),
        ("x_datum_max", ParamValue::Float(1e15)),
        ("y_datum_max", ParamValue::Float(1e15)),
    ]);
    td.fit_model_with_param("linear", &input).unwrap();
    assert!(td.model_parameters().source_equal(&input).unwrap());
    similar(td.apply(2.0).unwrap(), -95.0, "explicit linear");

    // With data, the defaults are filled in and the fitted line is stored.
    let mut fitted = TransformationDescription::new(
        [(0.0, 1.0), (1.0, 3.0)]
            .into_iter()
            .map(DataPoint::from)
            .collect(),
    )
    .unwrap();
    fitted
        .fit_model_with_param("linear", &Param::new())
        .unwrap();
    let stored = fitted.model_parameters();
    assert_eq!(stored.value("slope").unwrap(), &ParamValue::Float(2.0));
    assert_eq!(stored.value("intercept").unwrap(), &ParamValue::Float(1.0));
    assert_eq!(
        stored.value("symmetric_regression").unwrap(),
        &ParamValue::String("false".into())
    );
    assert_eq!(stored.value("x_weight").unwrap(), &s("x"));
    // `none` and `identity` keep nothing, whatever they are given.
    fitted.fit_model_with_param("none", &input).unwrap();
    assert!(fitted.model_parameters().is_empty());
}

#[test]
fn print_summary_matches_the_source_text() {
    // Section `printSummary`.
    let nonlinear: Vec<_> = [(0.0, 1.0), (0.25, 1.0625), (0.5, 1.25), (1.0, 2.0)]
        .into_iter()
        .map(DataPoint::from)
        .collect();
    let mut td = TransformationDescription::new(nonlinear).unwrap();
    let mut out = Vec::new();
    td.print_summary(&mut out).unwrap();
    let expected = "Number of data points (x/y pairs): 4\n\
        Data range (x): 0 to 1\n\
        Data range (y): 1 to 2\n\
        Summary of x/y deviations:\n\
        - 100% of data points within (+/-)1\n\
        -  99% of data points within (+/-)1\n\
        -  95% of data points within (+/-)1\n\
        -  90% of data points within (+/-)1\n\
        -  75% of data points within (+/-)1\n\
        -  50% of data points within (+/-)0.8125\n\
        -  25% of data points within (+/-)0.75\n\n";
    assert_eq!(String::from_utf8(out).unwrap(), expected);

    td.fit_model_with_param("linear", &Param::new()).unwrap();
    let mut out = Vec::new();
    td.print_summary(&mut out).unwrap();
    let expected = "Number of data points (x/y pairs): 4\n\
        Data range (x): 0 to 1\n\
        Data range (y): 1 to 2\n\
        Summary of x/y deviations before transformation:\n\
        - 100% of data points within (+/-)1\n\
        -  99% of data points within (+/-)1\n\
        -  95% of data points within (+/-)1\n\
        -  90% of data points within (+/-)1\n\
        -  75% of data points within (+/-)1\n\
        -  50% of data points within (+/-)0.8125\n\
        -  25% of data points within (+/-)0.75\n\
        Summary of x/y deviations after applying 'linear' transformation:\n\
        - 100% of data points within (+/-)0.142857\n\
        -  99% of data points within (+/-)0.125\n\
        -  95% of data points within (+/-)0.125\n\
        -  90% of data points within (+/-)0.125\n\
        -  75% of data points within (+/-)0.125\n\
        -  50% of data points within (+/-)0.0892857\n\
        -  25% of data points within (+/-)0.0714286\n\n";
    assert_eq!(String::from_utf8(out).unwrap(), expected);

    let mut out = Vec::new();
    TransformationDescription::default()
        .print_summary(&mut out)
        .unwrap();
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "Number of data points (x/y pairs): 0\n"
    );
}

#[test]
fn model_types_are_listed_in_source_order() {
    // Section `getModelTypes`.
    assert_eq!(
        TransformationDescription::model_types(),
        ["linear", "b_spline", "interpolated", "lowess"]
    );
}

#[test]
fn transformation_model_defaults_match_the_class_test() {
    // TransformationModelDefaults_test.cpp, first section.
    let models = ["linear", "b_spline", "lowess", "interpolated"];
    for selected in [
        "linear",
        "b_spline",
        "lowess",
        "interpolated",
        "none",
        "custom",
    ] {
        let params = transformation_model_defaults(selected).unwrap();
        assert_eq!(params.value("type").unwrap(), &s(selected));
        assert_eq!(params.description("type").unwrap(), "Type of model");
        let mut expected: Vec<String> = models.iter().map(|m| (*m).to_owned()).collect();
        if selected == "none" || selected == "custom" {
            expected.insert(0, selected.to_owned());
        }
        assert_eq!(params.valid_strings("type").unwrap(), expected.as_slice());
        for model in models {
            assert_eq!(
                params.section_description(model).unwrap(),
                format!("Parameters for '{model}' model")
            );
        }
        assert_eq!(
            params.value("b_spline:num_nodes").unwrap(),
            &ParamValue::Integer(5)
        );
        assert_eq!(
            params.value("b_spline:boundary_condition").unwrap(),
            &ParamValue::Integer(2)
        );
        similar(
            params.value("lowess:span").unwrap().to_f64().unwrap(),
            2.0 / 3.0,
            "lowess:span",
        );
        assert_eq!(params.entry("lowess:span").unwrap().min_float, 0.0);
        assert_eq!(params.entry("lowess:span").unwrap().max_float, 1.0);
        assert_eq!(
            params.value("interpolated:interpolation_type").unwrap(),
            &s("cspline")
        );
    }
}

#[test]
fn lowess_auto_span_prefers_the_larger_span_on_a_tie_and_persists_it() {
    // TransformationModelLowess_test.cpp, section "auto span selection
    // chooses larger span on tie and persists params".
    let data: Vec<_> = (0..=10)
        .map(|i| DataPoint::new(f64::from(i), 2.0 * f64::from(i) + 1.0))
        .collect();
    let mut p = LowessOptions::defaults().unwrap();
    for (key, value) in [
        ("span", ParamValue::Float(0.0)),
        ("auto_span", s("true")),
        ("auto_span_grid", s("0.3,0.8")),
        ("auto_metric", s("mae")),
        ("num_iterations", ParamValue::Integer(0)),
        ("delta", ParamValue::Float(-1.0)),
        ("interpolation_type", s("cspline")),
        ("extrapolation_type", s("four-point-linear")),
    ] {
        p.set_value(key, value, "", &[]).unwrap();
    }
    let mut td = TransformationDescription::new(data).unwrap();
    td.fit_model_with_param("lowess", &p).unwrap();
    let used = td.model_parameters();
    assert_eq!(used.value("auto_span").unwrap(), &s("false"));
    similar(used.value("span").unwrap().to_f64().unwrap(), 0.8, "span");
    similar(td.apply(0.0).unwrap(), 1.0, "f(0)");
    similar(td.apply(5.0).unwrap(), 11.0, "f(5)");
    similar(td.apply(10.0).unwrap(), 21.0, "f(10)");
    // The typed configuration carries the selected span too, so a refit
    // does not select again.
    let ModelConfig::Lowess(options) = td.model_config() else {
        panic!("an automatic LOWESS resolves to a fixed span");
    };
    assert_eq!(options.span, 0.8);
}

#[test]
fn lowess_reports_the_parameters_it_was_given() {
    // Section `getParameters`: extra keys of another model are ignored.
    let input = param(&[
        ("span", ParamValue::Float(0.3)),
        ("num_iterations", ParamValue::Integer(8)),
        ("delta", ParamValue::Float(1.0)),
        ("extrapolate", s("b_spline")),
        ("interpolation_type", s("cspline")),
        ("extrapolation_type", s("four-point-linear")),
    ]);
    let mut td = TransformationDescription::new(eleven()).unwrap();
    td.fit_model_with_param("lowess", &input).unwrap();
    assert_eq!(
        td.model_parameters().value("num_iterations").unwrap(),
        &ParamValue::Integer(8)
    );
    // Section `TransformationModelLowess(...)`: no data is an IllegalArgument.
    let defaults = LowessOptions::defaults().unwrap();
    let mut empty = TransformationDescription::default();
    assert!(empty.fit_model_with_param("lowess", &defaults).is_err());
}

#[test]
fn lowess_sine_through_the_parameter_route() {
    // Section `evaluate`, through `Param` as the source builds the model.
    let params = param(&[
        ("span", ParamValue::Float(0.3)),
        ("num_iterations", ParamValue::Integer(3)),
        ("delta", ParamValue::Float(-1.0)),
        ("interpolation_type", s("cspline")),
        ("extrapolation_type", s("four-point-linear")),
    ]);
    let mut td = TransformationDescription::new(sine()).unwrap();
    td.fit_model_with_param("lowess", &params).unwrap();
    for row in table(include_str!("data/transformations_lowess_golden.tsv")) {
        similar(td.apply(row[0]).unwrap(), row[1], "lowess sine");
    }
    for (name, left, right) in [
        ("two-point-linear", -0.04240732863, 0.046870277),
        ("global-linear", -0.9501004, 1.08486397),
    ] {
        let mut params = params.checked_clone().unwrap();
        params
            .set_value("extrapolation_type", s(name), "", &[])
            .unwrap();
        td.fit_model_with_param("lowess", &params).unwrap();
        similar(td.apply(-4.0).unwrap(), left, name);
        similar(td.apply(4.0).unwrap(), right, name);
    }
}

#[test]
fn akima_matches_the_gsl_reference_table() {
    // TransformationModelInterpolated_test.cpp, section `evaluate`: the
    // third column of TransformationModelInterpolated_gsl_data.txt.
    let points: Vec<_> = table(include_str!("data/transformations_interpolation_input.tsv"))
        .iter()
        .map(|r| DataPoint::new(r[0], r[1]))
        .collect();
    let mut p = InterpolationOptions::defaults().unwrap();
    p.set_value("interpolation_type", s("akima"), "", &[])
        .unwrap();
    let mut td = TransformationDescription::new(points).unwrap();
    td.fit_model_with_param("interpolated", &p).unwrap();
    for row in table(include_str!(
        "data/transformations_interpolation_golden.tsv"
    )) {
        similar(td.apply(row[0]).unwrap(), row[2], "akima");
    }
}

#[test]
fn interpolated_defaults_and_extrapolation_through_param() {
    // Sections `getDefaultParameters` and "[EXTRA] evaluate() beyond the
    // actual borders".
    let defaults = InterpolationOptions::defaults().unwrap();
    assert!(!defaults.is_empty());
    assert_eq!(defaults.value("interpolation_type").unwrap(), &s("cspline"));
    assert_eq!(
        defaults.value("extrapolation_type").unwrap(),
        &s("two-point-linear")
    );
    let dummy: Vec<_> = [(0.0, 1.0), (0.5, 4.0), (1.0, 2.0), (1.0, 4.0)]
        .into_iter()
        .map(DataPoint::from)
        .collect();
    let mut td = TransformationDescription::new(dummy).unwrap();
    td.fit_model_with_param("interpolated", &defaults).unwrap();
    for (x, y) in [(0.0, 1.0), (0.5, 4.0), (1.0, 3.0)] {
        similar(td.apply(x).unwrap(), y, "knots");
    }
    td.fit_model_with_param("interpolated", &Param::new())
        .unwrap();
    similar(td.apply(-0.5).unwrap(), 0.0, "default front");
    similar(td.apply(1.5).unwrap(), 4.0, "default back");
    for (name, left, right) in [
        ("two-point-linear", 0.0, 4.0),
        ("four-point-linear", -2.0, 2.0),
        ("global-linear", 0.909091, 4.1818182),
    ] {
        let p = param(&[("extrapolation_type", s(name))]);
        td.fit_model_with_param("interpolated", &p).unwrap();
        similar(td.apply(-0.5).unwrap(), left, name);
        similar(td.apply(1.5).unwrap(), right, name);
    }
    let unknown = param(&[("interpolation_type", s("polynomial"))]);
    assert!(td.fit_model_with_param("interpolated", &unknown).is_err());
}

#[test]
fn pair_overload_reproduces_the_source_zero_padding() {
    // The `(vector<pair>, Param, bool preprocess)` constructor. With
    // `preprocess` it averages like the DataPoints constructor, but its
    // global-linear line is fitted to three (0, 0) points and the three
    // averaged anchors (0, 1), (0.5, 4), (1, 3). Least squares over those six
    // points, derived by hand: mean x 0.25, mean y 4/3, Sxx 0.875, Sxy 3.0,
    // so slope 24/7 and intercept 4/3 - 6/7 = 10/21.
    let pairs = [(0.0, 1.0), (0.5, 4.0), (1.0, 2.0), (1.0, 4.0)];
    let model = InterpolatedModel::fit_pairs(
        &pairs,
        InterpolationOptions {
            extrapolation: Extrapolation::GlobalLinear,
            ..Default::default()
        },
        true,
    )
    .unwrap();
    let slope = 24.0 / 7.0;
    let intercept = 10.0 / 21.0;
    similar(model.apply(-0.5).unwrap(), intercept - 0.5 * slope, "front");
    similar(model.apply(1.5).unwrap(), intercept + 1.5 * slope, "back");
    assert_eq!(model.apply(1.0).unwrap(), 3.0);
    // Two-point extrapolation is unaffected by the padding.
    let model =
        InterpolatedModel::fit_pairs(&pairs, InterpolationOptions::default(), true).unwrap();
    similar(model.apply(-0.5).unwrap(), 0.0, "two-point front");
    // Without preprocessing, the padded x sequence repeats 0, which a cubic
    // spline cannot interpolate (the source divides by the zero interval).
    assert!(
        InterpolatedModel::fit_pairs(&[(1.0, 1.0), (2.0, 2.0)], Default::default(), false).is_err()
    );
    // Linear interpolation runs on (0,0),(0,0),(1,1),(2,2): upper_bound of 0.5
    // is position 2, the chord from (0,0) to (1,1).
    let model = InterpolatedModel::fit_pairs(
        &[(1.0, 1.0), (2.0, 2.0)],
        InterpolationOptions {
            interpolation: Interpolation::Linear,
            ..Default::default()
        },
        false,
    )
    .unwrap();
    assert_eq!(model.apply(0.5).unwrap(), 0.5);
}

#[test]
fn base_model_weighting_matches_the_class_test() {
    // TransformationModel_test.cpp.
    assert!(check_valid_weight("ln(x)", &VALID_X_WEIGHTS));
    assert!(check_valid_weight("1/y", &VALID_Y_WEIGHTS));
    assert!(check_valid_weight("1/x2", &VALID_X_WEIGHTS));
    assert!(check_valid_weight("x", &VALID_X_WEIGHTS));
    assert!(check_valid_weight("y", &VALID_Y_WEIGHTS));
    assert!(!check_valid_weight("none", &VALID_X_WEIGHTS));
    assert!(!check_valid_weight("x2", &VALID_X_WEIGHTS));

    let inf = f64::INFINITY;
    for (weight, zero, two) in [
        ("x", 0.0, 2.0),
        ("none", 0.0, 2.0),
        ("ln(x)", -inf, 2f64.ln()),
        ("1/x", inf, 0.5),
        ("1/x2", inf, 0.25),
        ("ln(y)", -inf, 2f64.ln()),
        ("1/y", inf, 0.5),
        ("1/y2", inf, 0.25),
    ] {
        assert_eq!(weight_datum(0.0, weight), zero, "{weight}(0)");
        similar(weight_datum(2.0, weight), two, weight);
    }

    let (dmin, dmax) = (10e-6, 10e9);
    similar(check_datum_range(10e-7, dmin, dmax), dmin, "below");
    similar(check_datum_range(10e12, dmin, dmax), dmax, "above");
    similar(check_datum_range(100.0, dmin, dmax), 100.0, "inside");
    let (xmin, xmax) = (10e-5, 10e12);
    similar(
        weight_datum(check_datum_range(10e13, xmin, xmax), "x"),
        xmax,
        "clamped x",
    );
    similar(
        weight_datum(check_datum_range(0.0, xmin, xmax), "1/x2"),
        1.0 / (xmin * xmin),
        "clamped 1/x2",
    );

    for (weight, value, expected) in [
        ("ln(x)", 11.0e5f64.ln(), 11e5),
        ("ln(x)", 2.0, 2f64.exp()),
        ("1/x", 1.0 / 9.0e-5, 9.0e-5),
        ("1/x2", 1.0 / (9.0e-5 * 9.0e-5), 9.0e-5),
        ("1/x2", 2.0, (0.5f64).sqrt()),
        ("ln(y)", 11.0e8f64.ln(), 11.0e8),
        ("1/y", 1.0 / 9.0e-8, 9e-8),
        ("1/y2", 1.0 / (9.0e-8 * 9.0e-8), 9e-8),
        ("none", 2.0, 2.0),
    ] {
        similar(unweight_datum(value, weight), expected, weight);
    }

    // Sections `weightData` and `unWeightData`, with the base model's weights
    // read from a Param (the default datum bounds).
    let data = || -> Vec<DataPoint> {
        [(1e-4, 1.0), (1.0, 2.0), (2.0, 4.0)]
            .into_iter()
            .map(DataPoint::from)
            .collect()
    };
    let log_x = ModelWeights::from_param(&param(&[("x_weight", s("ln(x)")), ("y_weight", s("y"))]))
        .unwrap();
    assert!(log_x.is_weighting());
    let mut weighted = data();
    log_x.weight_data(&mut weighted).unwrap();
    for (point, (x, y)) in weighted
        .iter()
        .zip([(1e-4f64.ln(), 1.0), (0.0, 2.0), (2f64.ln(), 4.0)])
    {
        similar(point.x, x, "ln(x) weighted x");
        similar(point.y, y, "ln(x) weighted y");
    }
    let log_y = ModelWeights::from_param(&param(&[("x_weight", s("x")), ("y_weight", s("ln(y)"))]))
        .unwrap();
    let mut unweighted: Vec<_> = [(0.0, 1.0), (1.0, 2.0), (2.0, 4.0)]
        .into_iter()
        .map(DataPoint::from)
        .collect();
    log_y.unweight_data(&mut unweighted).unwrap();
    for (point, (x, y)) in
        unweighted
            .iter()
            .zip([(0.0, 1f64.exp()), (1.0, 2f64.exp()), (2.0, 4f64.exp())])
    {
        similar(point.x, x, "ln(y) unweighted x");
        similar(point.y, y, "ln(y) unweighted y");
    }
    // An x weight named for the y axis is invalid, as in the source.
    assert!(ModelWeights::from_param(&param(&[("x_weight", s("ln(y)"))])).is_err());
    assert!(!ModelWeights::default().is_weighting());
}

#[test]
fn data_point_order_and_equality_follow_the_source_tie() {
    // Sections `DataPoint::operator<` and `operator==`.
    let mut p1 = DataPoint::with_note(1.0, 2.0, "abc");
    let mut p2 = DataPoint::from((1.0, 2.0));
    assert_ne!(p1.partial_cmp(&p2), Some(Ordering::Less));
    assert!(p2 < p1);
    p2.note = "def".into();
    assert!(p1 < p2);
    assert_ne!(p2.partial_cmp(&p1), Some(Ordering::Less));
    p1.x = 1.5;
    assert_ne!(p1.partial_cmp(&p2), Some(Ordering::Less));
    assert!(p2 < p1);
    let p1 = DataPoint::with_note(1.0, 2.0, "abc");
    let mut p2 = DataPoint::from((1.0, 2.0));
    assert_ne!(p1, p2);
    p2.note = "abc".into();
    assert_eq!(p1, p2);
}

#[test]
fn span_grid_and_residual_scores_follow_the_source_rules() {
    // buildSpanGrid: default grid clamped to [0.15, 0.8] with 5/11 neighbours
    // as the lower bound (0.4545...), duplicates removed.
    let grid = build_span_grid(11, &[], 0.15, 0.80, 5).unwrap();
    let lower = 5.0 / 11.0;
    assert_eq!(grid, vec![lower, 0.5, 0.6, 0.7, 0.8]);
    // Crossed bounds are swapped; an out-of-range grid collapses to a bound.
    assert_eq!(
        build_span_grid(1000, &[0.95], 0.9, 0.2, 5).unwrap(),
        vec![0.9]
    );
    assert!(build_span_grid(10, &[f64::NAN], 0.1, 0.9, 3).is_err());
    // scoreResiduals.
    let errors = [1.0, 3.0, 2.0, 4.0];
    assert_eq!(
        score_residuals(&[], CrossValidationMetric::Mae).unwrap(),
        f64::INFINITY
    );
    assert_eq!(
        score_residuals(&errors, CrossValidationMetric::Mae).unwrap(),
        2.5
    );
    similar(
        score_residuals(&errors, CrossValidationMetric::Rmse).unwrap(),
        7.5f64.sqrt(),
        "rmse",
    );
    // Type-7 quantile of 1,2,3,4 at 0.9: position 2.7, 3 + 0.7 = 3.7.
    similar(
        score_residuals(&errors, CrossValidationMetric::P90).unwrap(),
        3.7,
        "p90",
    );
    assert!(CrossValidationMetric::from_name("median").is_err());
    let mut p = LowessOptions::defaults().unwrap();
    p.set_value("auto_metric", s("median"), "", &[]).unwrap();
    // The source converts auto_metric on every fit, enabled or not.
    assert!(LowessOptions::from_param(&p).is_err());
}

#[test]
fn inverted_explicit_linear_parameters_swap_as_the_source_does() {
    // TransformationModelLinear::invert updates slope, intercept, both weight
    // names and all four bounds in the model's parameters.
    let input = param(&[
        ("slope", ParamValue::Float(2.0)),
        ("intercept", ParamValue::Float(4.0)),
        ("x_weight", s("ln(x)")),
        ("x_datum_min", ParamValue::Float(0.5)),
    ]);
    let mut td = TransformationDescription::default();
    td.fit_model_with_param("linear", &input).unwrap();
    td.invert().unwrap();
    let p = td.model_parameters();
    assert_eq!(p.value("slope").unwrap(), &ParamValue::Float(0.5));
    assert_eq!(p.value("intercept").unwrap(), &ParamValue::Float(-2.0));
    assert_eq!(p.value("x_weight").unwrap(), &s("y"));
    assert_eq!(p.value("y_weight").unwrap(), &s("ln(x)"));
    assert_eq!(p.value("x_datum_min").unwrap(), &ParamValue::Float(1e-15));
    assert_eq!(p.value("y_datum_min").unwrap(), &ParamValue::Float(0.5));
    assert_eq!(p.value("y_datum_max").unwrap(), &ParamValue::Float(1e15));
}
