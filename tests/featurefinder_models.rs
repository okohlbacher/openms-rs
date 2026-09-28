// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Class tests of the feature finder's one-dimensional models and fitters,
//! transcribed from the upstream class tests at core `bc9cc12`
//! (`src/tests/class_tests/openms/source/<Class>_test.cpp`): LinearInterpolation,
//! BaseModel, ModelDescription, InterpolationModel, GaussModel, BiGaussModel,
//! EmgModel, IsotopeModel, ExtendedIsotopeModel, Fitter1D, LevMarqFitter1D,
//! MaxLikeliFitter1D, GaussFitter1D, BiGaussFitter1D, EmgFitter1D,
//! IsotopeFitter1D, ExtendedIsotopeFitter1D, EmgScoring and SeedListGenerator.
//!
//! ClassTest's `TOLERANCE_ABSOLUTE` and `TOLERANCE_RELATIVE` persist from the
//! section that sets them to the end of the test file; each assertion here
//! carries the pair in force at its source line.

use openms::analysis::featurefinder::base_model::{
    BaseModel, InterpolationModel, Model1D, ModelDescription,
};
use openms::analysis::featurefinder::emg_scoring::EmgScoring;
use openms::analysis::featurefinder::feature_finder_algorithm::{
    FeatureFinderAlgorithm, SEEDS_NOT_SUPPORTED, Summary,
};
use openms::analysis::featurefinder::fitter1d::{
    BaseFitter1D, Fitter1D, MAX_ITERATION_DEFAULT, lev_marq_optimize,
};
use openms::analysis::featurefinder::fitters::{
    BiGaussFitter1D, EmgFitter1D, ExtendedIsotopeFitter1D, GaussFitter1D, IsotopeFitter1D,
};
use openms::analysis::featurefinder::gauss_models::{BiGaussModel, EmgModel, GaussModel};
use openms::analysis::featurefinder::isotope_model::{ExtendedIsotopeModel, IsotopeModel};
use openms::analysis::featurefinder::seed_list_generator::SeedListGenerator;
use openms::chemistry::{AASequence, EmpiricalFormula};
use openms::data_structures::DPosition2;
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::kernel::features::{Feature, FeatureMap};
use openms::kernel::mrm::{MRMFeature, MRMTransitionGroup};
use openms::kernel::{
    ConvexHull2D, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, Point2D, Precursor,
};
use openms::math::ml::linear_interpolation::LinearInterpolation;
use openms::param::{Param, ParamValue};

// ---------------------------------------------------------------------------
// ClassTest comparison

/// ClassTest `isRealSimilar`: similar when the absolute difference is within
/// `abs` or the ratio within `ratio` (ClassTest.cpp; defaults 1e-5 and
/// 1 + 1e-5).
fn is_real_similar(number_1: f64, number_2: f64, abs: f64, ratio_max: f64) -> bool {
    if number_1.is_nan() || number_2.is_nan() {
        return false;
    }
    let is_absdiff_small = (number_1 - number_2).abs() <= abs;
    if number_1 == 0.0 || number_2 == 0.0 {
        return (number_1 == 0.0 && number_2 == 0.0) || is_absdiff_small;
    }
    let mut ratio = number_1 / number_2;
    if ratio < 0.0 {
        return is_absdiff_small;
    }
    if ratio < 1.0 {
        ratio = 1.0 / ratio;
    }
    ratio <= ratio_max || is_absdiff_small
}

const DEFAULT_ABS: f64 = 1e-5;
const DEFAULT_RATIO: f64 = 1.0 + 1e-5;

#[track_caller]
fn similar_with(actual: f64, expected: f64, abs: f64, ratio: f64) {
    assert!(
        is_real_similar(actual, expected, abs, ratio),
        "{actual} is not similar to {expected} (abs {abs}, ratio {ratio})"
    );
}

#[track_caller]
fn similar(actual: f64, expected: f64) {
    similar_with(actual, expected, DEFAULT_ABS, DEFAULT_RATIO);
}

fn param(entries: &[(&str, ParamValue)]) -> Param {
    let mut p = Param::new();
    for (key, value) in entries {
        p.set_value(key, value.clone(), "", &[]).unwrap();
    }
    p
}

fn f(value: f64) -> ParamValue {
    ParamValue::Float(value)
}

fn i(value: i64) -> ParamValue {
    ParamValue::Integer(value)
}

#[track_caller]
fn params_equal(a: &Param, b: &Param) {
    assert!(a.source_equal(b).unwrap(), "{a:?}\n!=\n{b:?}");
}

#[track_caller]
fn samples_similar(a: &[Peak1D], b: &[Peak1D], abs: f64) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        similar_with(x.mz, y.mz, abs, DEFAULT_RATIO);
        similar_with(
            f64::from(x.intensity),
            f64::from(y.intensity),
            abs,
            DEFAULT_RATIO,
        );
    }
}

// ---------------------------------------------------------------------------
// LinearInterpolation_test.cpp

#[test]
fn linear_interpolation_constructor_and_accessors() {
    let lifd0 = LinearInterpolation::default();
    let lifd1 = LinearInterpolation::new(1.125, 0.0);
    let lifd2 = LinearInterpolation::new(1.125, 3.5);
    assert_eq!(lifd0.scale(), 1.0);
    assert_eq!(lifd0.offset(), 0.0);
    assert_eq!(lifd1.scale(), 1.125);
    assert_eq!(lifd1.offset(), 0.0);
    assert_eq!(lifd2.scale(), 1.125);
    assert_eq!(lifd2.offset(), 3.5);

    let v = vec![17.0, 18.9, 20.333, -0.1];
    let mut lifd = LinearInterpolation::default();
    lifd.set_data(v.clone());
    assert_eq!(lifd.data(), v.as_slice());

    let mut lifd = LinearInterpolation::default();
    assert!(lifd.is_empty());
    lifd.data_mut().push(3.0);
    assert!(!lifd.is_empty());
    lifd.data_mut().push(3.0);
    lifd.data_mut().push(3.0);
    assert!(!lifd.is_empty());
    lifd.data_mut().clear();
    assert!(lifd.is_empty());
}

#[test]
fn linear_interpolation_mapping() {
    let mut lifd = LinearInterpolation::default();
    lifd.set_mapping(13.0, 23.0, 53.0);
    similar(lifd.scale(), 13.0);
    similar(lifd.inside_reference_point(), 23.0);
    similar(lifd.outside_reference_point(), 53.0);
    lifd.set_scale(88.88);
    similar(lifd.scale(), 88.88);

    let mut lifd = LinearInterpolation::new(1.125, 3.5);
    assert_eq!(lifd.offset(), 3.5);
    lifd.set_offset(88.88);
    assert_eq!(lifd.offset(), 88.88);

    let mut lifd = LinearInterpolation::default();
    lifd.set_mapping_by_points(13.0, 130.0, 14.0, 140.0);
    similar(lifd.scale(), 10.0);
    similar(lifd.inside_reference_point(), 13.0);
    similar(lifd.outside_reference_point(), 130.0);

    // Copy constructor and assignment.
    let v = vec![17.0, 18.9, 20.333, -0.1];
    lifd.set_data(v.clone());
    let lifd2 = lifd.clone();
    similar(lifd2.scale(), 10.0);
    similar(lifd2.inside_reference_point(), 13.0);
    similar(lifd2.outside_reference_point(), 130.0);
    assert_eq!(lifd2.data(), v.as_slice());
}

#[test]
fn linear_interpolation_key_transformations_and_support() {
    let mut lifd = LinearInterpolation::new(100.0, 3456.0);
    similar(lifd.index2key(0.0), 3456.0);
    similar(lifd.index2key(-1.0), 3356.0);
    similar(lifd.index2key(1.0), 3556.0);
    similar(lifd.key2index(3456.0), 0.0);
    similar(lifd.key2index(3356.0), -1.0);
    similar(lifd.key2index(3556.0), 1.0);
    lifd.set_scale(0.0);
    similar(lifd.key2index(3456.0), 0.0);
    similar(lifd.key2index(3356.0), 0.0);
    similar(lifd.key2index(3556.0), 0.0);

    let mut lifd = LinearInterpolation::new(1.125, 3.5);
    similar(lifd.support_min(), 3.5);
    similar(lifd.support_max(), 3.5);
    lifd.data_mut().push(11111.0);
    similar(lifd.support_min(), 3.5 - 1.125);
    similar(lifd.support_max(), 3.5 + 1.125);
    lifd.data_mut().push(99999.0);
    similar(lifd.support_min(), 3.5 - 1.125);
    similar(lifd.support_max(), 3.5 + 2.0 * 1.125);
    lifd.data_mut().clear();
    similar(lifd.support_min(), 3.5);
    similar(lifd.support_max(), 3.5);
}

#[test]
fn linear_interpolation_value() {
    let values = [1.0, 2.0, 0.0, 1.0];
    let mut lifd0 = LinearInterpolation::default();
    lifd0.set_data(values.to_vec());
    for (k, &v) in values.iter().enumerate() {
        assert_eq!(lifd0.value(k as f64), v);
    }
    let inter_values = [
        0.0, 0.00, 0.00, 0.00, //
        0.0, 0.25, 0.50, 0.75, //
        1.0, 1.25, 1.50, 1.75, //
        2.0, 1.50, 1.00, 0.50, //
        0.0, 0.25, 0.50, 0.75, //
        1.0, 0.75, 0.50, 0.25, //
        0.0, 0.00, 0.00, 0.00, //
        0.0,
    ];
    for k in 0..values.len() + 4 {
        similar(lifd0.value(k as f64 - 2.0), inter_values[4 * k]);
    }
    for (k, &expected) in inter_values.iter().enumerate() {
        similar(lifd0.value((k as f64 - 8.0) / 4.0), expected);
    }
    let mut lifd1 = lifd0.clone();
    let (scale, offset) = (1.0, 100.0);
    lifd1.set_scale(scale);
    lifd1.set_offset(offset);
    for k in -8..(inter_values.len() as i32 - 8) {
        let pos = f64::from(k) / 4.0;
        similar(lifd1.key2index(lifd1.index2key(pos)), pos);
        similar(lifd1.value(pos * scale + offset), lifd0.value(pos));
    }

    // TOLERANCE_ABSOLUTE(0.001)
    let mut small = LinearInterpolation::default();
    small.data_mut().resize(5, 0.0);
    small.set_mapping_by_points(0.0, 0.0, 5.0, 5.0);
    let mut big = LinearInterpolation::default();
    big.data_mut().resize(15, 0.0);
    big.set_mapping_by_points(5.0, 0.0, 10.0, 5.0);
    for k in 0..5 {
        let v = k as f64 * 25.0 + 100.0;
        small.data_mut()[k] = v;
        big.data_mut()[k + 5] = v;
    }
    for k in -50..=100 {
        let pos = f64::from(k) / 10.0;
        similar_with(small.value(pos), big.value(pos), 0.001, DEFAULT_RATIO);
    }
}

#[test]
fn linear_interpolation_derivative() {
    let mut lifd0 = LinearInterpolation::default();
    lifd0.set_data(vec![1.0, 2.0, 0.0, 1.0]);
    let inter_values = [
        0.00, 0.00, 0.00, 0.25, //
        0.50, 0.75, 1.00, 1.00, //
        1.00, 1.00, 1.00, 0.25, //
        -0.50, -1.25, -2.00, -1.25, //
        -0.50, 0.25, 1.00, 0.50, //
        0.00, -0.50, -1.00, -0.75, //
        -0.50, -0.25, 0.00, 0.00, //
        0.0,
    ];
    for k in -8..(inter_values.len() as i32 - 8) {
        let key = f64::from(k) / 4.0;
        let index = (k + 8) as usize;
        // TOLERANCE_ABSOLUTE(0.001) is in force from the value section.
        similar_with(
            lifd0.derivative(key),
            inter_values[index],
            0.001,
            DEFAULT_RATIO,
        );
    }
}

#[test]
fn linear_interpolation_add_value() {
    let cases: [(f64, &[(usize, f64)]); 7] = [
        (2.3, &[(2, 7.0), (3, 3.0)]),
        (0.3, &[(0, 7.0), (1, 3.0)]),
        (-0.7, &[(0, 3.0)]),
        (-1.7, &[(0, 0.0)]),
        (3.3, &[(3, 7.0), (4, 3.0)]),
        (4.3, &[(4, 7.0)]),
        (5.3, &[(4, 0.0)]),
    ];
    for (pos, expected) in cases {
        let mut li = LinearInterpolation::default();
        li.data_mut().resize(5, 0.0);
        li.add_value(pos, 10.0);
        for &(index, value) in expected {
            similar_with(li.data()[index], value, 0.001, DEFAULT_RATIO);
        }
    }
    for k in -50..=100 {
        let pos = f64::from(k) / 10.0;
        let mut small = LinearInterpolation::default();
        small.data_mut().resize(5, 0.0);
        small.set_mapping_by_points(0.0, 0.0, 5.0, 5.0);
        small.add_value(pos, 10.0);
        let small: Vec<f64> = small.data().iter().map(|v| v.round()).collect();
        let mut big = LinearInterpolation::default();
        big.data_mut().resize(15, 0.0);
        big.set_mapping_by_points(5.0, 0.0, 10.0, 5.0);
        big.add_value(pos, 10.0);
        let big: Vec<f64> = big.data()[5..10].iter().map(|v| v.round()).collect();
        assert_eq!(small, big, "at {pos}");
    }
}

// ---------------------------------------------------------------------------
// GaussModel_test.cpp

fn gauss_tmp(min: f64, max: f64, mean: f64, variance: f64) -> Param {
    param(&[
        ("bounding_box:min", f(min)),
        ("bounding_box:max", f(max)),
        ("statistics:variance", f(variance)),
        ("statistics:mean", f(mean)),
    ])
}

#[test]
fn gauss_model_construction_copy_and_assignment() {
    let model = GaussModel::new().unwrap();
    assert_eq!(model.name(), "GaussModel");

    let tmp = gauss_tmp(678.9, 789.0, 680.1, 2.0);
    let mut gm1 = GaussModel::new().unwrap();
    gm1.set_scaling_factor(10.0);
    gm1.set_interpolation_step(0.3);
    gm1.set_parameters(&tmp).unwrap();
    let gm2 = gm1.clone();
    let mut gm3 = GaussModel::new().unwrap();
    gm3.set_scaling_factor(10.0);
    gm3.set_interpolation_step(0.3);
    gm3.set_parameters(&tmp).unwrap();
    params_equal(&gm3.to_param().unwrap(), &gm2.to_param().unwrap());
}

#[test]
fn gauss_model_set_parameters_and_offset() {
    // TOLERANCE_ABSOLUTE(0.001)
    let mut gm1 = GaussModel::new().unwrap();
    gm1.set_scaling_factor(10.0);
    gm1.set_parameters(&gauss_tmp(678.9, 789.0, 679.1, 2.0))
        .unwrap();
    gm1.set_offset(680.0);
    similar_with(gm1.center().unwrap(), 680.2, 0.001, DEFAULT_RATIO);
    let mut gm2 = GaussModel::new().unwrap();
    gm2.set_parameters(&gm1.to_param().unwrap()).unwrap();
    // TOLERANCE_ABSOLUTE(0.0000001)
    samples_similar(&gm1.samples(), &gm2.samples(), 1e-7);
}

#[test]
fn gauss_model_intensities() {
    let mut gm1 = GaussModel::new().unwrap();
    gm1.set_interpolation_step(0.001);
    gm1.set_parameters(&gauss_tmp(-4.0, 4.0, 0.0, 1.0)).unwrap();
    similar_with(gm1.center().unwrap(), 0.0, 1e-7, DEFAULT_RATIO);
    // TOLERANCE_ABSOLUTE(0.001)
    let tol = 0.001;
    let expected = [
        (-1.0, 0.24197072),
        (0.0, 0.39894228),
        (1.0, 0.24197072),
        (2.0, 0.05399097),
    ];
    for (x, y) in expected {
        similar_with(gm1.intensity(x), y, tol, DEFAULT_RATIO);
    }
    gm1.set_interpolation_step(0.2);
    gm1.set_samples().unwrap();
    for (x, y) in expected {
        similar_with(gm1.intensity(x), y, tol, DEFAULT_RATIO);
    }
    gm1.set_scaling_factor(10.0);
    gm1.set_samples().unwrap();
    for (x, y) in expected {
        similar_with(gm1.intensity(x), 10.0 * y, tol, DEFAULT_RATIO);
    }
}

#[test]
fn gauss_model_center() {
    let mut gm1 = GaussModel::new().unwrap();
    gm1.set_scaling_factor(10.0);
    gm1.set_parameters(&gauss_tmp(678.9, 789.0, 679.1, 2.0))
        .unwrap();
    gm1.set_offset(680.0);
    similar_with(gm1.center().unwrap(), 680.2, 0.001, DEFAULT_RATIO);

    let mut gm1 = GaussModel::new().unwrap();
    gm1.set_scaling_factor(10.0);
    gm1.set_parameters(&gauss_tmp(650.0, 750.0, 679.1, 2.0))
        .unwrap();
    similar_with(gm1.center().unwrap(), 679.1, 0.001, DEFAULT_RATIO);
}

// ---------------------------------------------------------------------------
// BiGaussModel_test.cpp

fn bigauss_tmp(min: f64, max: f64, mean: f64, v1: f64, v2: f64) -> Param {
    param(&[
        ("bounding_box:min", f(min)),
        ("bounding_box:max", f(max)),
        ("statistics:mean", f(mean)),
        ("statistics:variance1", f(v1)),
        ("statistics:variance2", f(v2)),
    ])
}

#[test]
fn bigauss_model_construction_copy_and_assignment() {
    let model = BiGaussModel::new().unwrap();
    assert_eq!(model.name(), "BiGaussModel");
    let as_base: Model1D = Model1D::BiGauss(model);
    assert_eq!(as_base.name(), "BiGaussModel");

    let tmp = bigauss_tmp(678.9, 789.0, 680.1, 2.0, 5.0);
    let mut bgm1 = BiGaussModel::new().unwrap();
    bgm1.set_scaling_factor(10.0);
    bgm1.set_interpolation_step(0.3);
    bgm1.set_parameters(&tmp).unwrap();
    let bgm2 = bgm1.clone();
    let mut bgm3 = BiGaussModel::new().unwrap();
    bgm3.set_scaling_factor(10.0);
    bgm3.set_interpolation_step(0.3);
    bgm3.set_parameters(&tmp).unwrap();
    params_equal(&bgm3.to_param().unwrap(), &bgm2.to_param().unwrap());
}

#[test]
fn bigauss_model_set_parameters_offset_and_center() {
    // TOLERANCE_ABSOLUTE(0.001)
    let mut bgm1 = BiGaussModel::new().unwrap();
    bgm1.set_parameters(&bigauss_tmp(678.9, 789.0, 680.1, 2.0, 5.0))
        .unwrap();
    bgm1.set_offset(680.0);
    let mut bgm2 = BiGaussModel::new().unwrap();
    bgm2.set_parameters(&bgm1.to_param().unwrap()).unwrap();
    similar_with(bgm1.center().unwrap(), 681.2, 0.001, DEFAULT_RATIO);
    // TOLERANCE_ABSOLUTE(0.0001)
    samples_similar(&bgm1.samples(), &bgm2.samples(), 0.0001);

    // setOffset section (the tolerance in force is still 0.0001).
    let mut bgm1 = BiGaussModel::new().unwrap();
    bgm1.set_parameters(&bigauss_tmp(678.9, 789.0, 680.1, 2.0, 5.0))
        .unwrap();
    bgm1.set_offset(680.9);
    let mut bgm2 = BiGaussModel::new().unwrap();
    bgm2.set_parameters(&bigauss_tmp(680.9, 791.0, 682.1, 2.0, 5.0))
        .unwrap();
    params_equal(&bgm1.to_param().unwrap(), &bgm2.to_param().unwrap());
    similar_with(
        bgm1.center().unwrap(),
        bgm2.center().unwrap(),
        0.0001,
        DEFAULT_RATIO,
    );
    similar_with(bgm1.center().unwrap(), 682.1, 0.0001, DEFAULT_RATIO);
    // TOLERANCE_ABSOLUTE(0.001)
    samples_similar(&bgm1.samples(), &bgm2.samples(), 0.001);

    bgm1.set_parameters(&bigauss_tmp(-4.0, 4.001, 0.0, 0.81, 0.81))
        .unwrap();
    bgm1.set_offset(0.123);
    similar_with(bgm1.center().unwrap(), 4.123, 0.001, DEFAULT_RATIO);
    similar_with(bgm1.intensity(4.123), 0.4432692, 0.001, DEFAULT_RATIO);
    similar_with(
        bgm1.intensity(4.223),
        bgm1.intensity(4.023),
        0.001,
        DEFAULT_RATIO,
    );
    similar_with(
        bgm1.intensity(3.123),
        bgm1.intensity(5.123),
        0.001,
        DEFAULT_RATIO,
    );

    let mut bgm1 = BiGaussModel::new().unwrap();
    bgm1.set_parameters(&bigauss_tmp(678.9, 789.0, 680.1, 2.0, 5.0))
        .unwrap();
    bgm1.set_offset(680.0);
    similar_with(bgm1.center().unwrap(), 681.2, 0.001, DEFAULT_RATIO);
}

// ---------------------------------------------------------------------------
// EmgModel_test.cpp

#[allow(clippy::too_many_arguments)]
fn emg_tmp(
    min: f64,
    max: f64,
    mean: f64,
    variance: f64,
    height: f64,
    width: f64,
    symmetry: f64,
    retention: f64,
) -> Param {
    param(&[
        ("bounding_box:min", f(min)),
        ("bounding_box:max", f(max)),
        ("statistics:mean", f(mean)),
        ("statistics:variance", f(variance)),
        ("emg:height", f(height)),
        ("emg:width", f(width)),
        ("emg:symmetry", f(symmetry)),
        ("emg:retention", f(retention)),
    ])
}

#[test]
fn emg_model_construction_copy_and_assignment() {
    assert_eq!(EmgModel::new().unwrap().name(), "EmgModel");
    let tmp = emg_tmp(678.9, 789.0, 680.1, 2.0, 100000.0, 5.0, 5.0, 725.0);
    let mut em1 = EmgModel::new().unwrap();
    em1.set_interpolation_step(0.2);
    em1.set_parameters(&tmp).unwrap();
    let em2 = em1.clone();
    let mut em3 = EmgModel::new().unwrap();
    em3.set_interpolation_step(0.2);
    em3.set_parameters(&tmp).unwrap();
    params_equal(&em3.to_param().unwrap(), &em2.to_param().unwrap());
}

#[test]
fn emg_model_set_parameters_and_intensities() {
    // TOLERANCE_ABSOLUTE(0.001)
    let mut em1 = EmgModel::new().unwrap();
    em1.set_parameters(&emg_tmp(
        678.9, 680.9, 679.1, 2.0, 100000.0, 5.0, 5.0, 1200.0,
    ))
    .unwrap();
    em1.set_offset(680.0);
    similar_with(em1.center().unwrap(), 680.2, 0.001, DEFAULT_RATIO);
    let mut em3 = EmgModel::new().unwrap();
    em3.set_parameters(&em1.to_param().unwrap()).unwrap();
    // TOLERANCE_ABSOLUTE(0.0001)
    samples_similar(&em1.samples(), &em3.samples(), 0.0001);

    let mut em2 = EmgModel::new().unwrap();
    em2.set_interpolation_step(0.1);
    em2.set_parameters(&emg_tmp(-1.0, 4.0, 0.0, 0.1, 10.0, 1.0, 2.0, 3.0))
        .unwrap();
    similar_with(em2.center().unwrap(), 0.0, 0.0001, DEFAULT_RATIO);
    // TOLERANCE_ABSOLUTE(0.01)
    let expected = [
        (-1.0, 0.0497198),
        (0.0, 0.164882),
        (1.0, 0.54166),
        (2.0, 1.69364),
    ];
    for (x, y) in expected {
        similar_with(em2.intensity(x), y, 0.01, DEFAULT_RATIO);
    }
    em2.set_interpolation_step(0.2);
    em2.set_samples().unwrap();
    for (x, y) in expected {
        similar_with(em2.intensity(x), y, 0.01, DEFAULT_RATIO);
    }
    // "checked small values of parameter symmetry"
    em2.set_parameters(&emg_tmp(0.0, 10.0, 0.0, 0.1, 10.0, 6.0, 1.0, 3.0))
        .unwrap();
    similar_with(em2.intensity(2.0), 747203.0, 0.01, DEFAULT_RATIO);
    // The source's three ABORT_IF(isinf(...)) lines for symmetry 0.1, 0.16
    // and 0.17 skip the rest of the section without changing its status; they
    // assert nothing and are not transcribed.
}

#[test]
fn emg_model_offset_and_center() {
    let tmp = emg_tmp(678.9, 789.0, 680.1, 2.0, 100000.0, 5.0, 5.0, 725.0);
    let mut em1 = EmgModel::new().unwrap();
    em1.set_parameters(&tmp).unwrap();
    em1.set_offset(680.9);
    let mut em2 = EmgModel::new().unwrap();
    em2.set_parameters(&tmp).unwrap();
    em2.set_offset(680.9);
    params_equal(&em1.to_param().unwrap(), &em2.to_param().unwrap());
    // The tolerance in force is TOLERANCE_ABSOLUTE(0.01) from the section
    // before.
    similar_with(
        em1.center().unwrap(),
        em2.center().unwrap(),
        0.01,
        DEFAULT_RATIO,
    );
    similar_with(em1.center().unwrap(), 682.1, 0.01, DEFAULT_RATIO);
    samples_similar(&em1.samples(), &em2.samples(), 0.01);

    // TOLERANCE_ABSOLUTE(0.001)
    let mut em1 = EmgModel::new().unwrap();
    em1.set_parameters(&tmp).unwrap();
    em1.set_offset(680.0);
    similar_with(em1.center().unwrap(), 681.2, 0.001, DEFAULT_RATIO);
}

// ---------------------------------------------------------------------------
// IsotopeModel_test.cpp

fn isotope_tmp() -> Param {
    param(&[
        ("charge", i(3)),
        ("isotope:mode:GaussianSD", f(0.8)),
        ("statistics:mean", f(670.5)),
    ])
}

#[test]
fn isotope_model_construction_copy_and_parameters() {
    let im = IsotopeModel::new().unwrap();
    assert_eq!(im.name(), "IsotopeModel");
    assert_eq!(im.charge(), 1);
    assert_eq!(im.center().unwrap(), 0.0);

    let mut im1 = IsotopeModel::new().unwrap();
    im1.set_parameters(&isotope_tmp()).unwrap();
    let im2 = im1.clone();
    let mut im3 = IsotopeModel::new().unwrap();
    im3.set_parameters(&isotope_tmp()).unwrap();
    params_equal(&im3.to_param().unwrap(), &im2.to_param().unwrap());

    let mut im2 = IsotopeModel::new().unwrap();
    im2.set_parameters(&im1.to_param().unwrap()).unwrap();
    samples_similar(&im1.samples(), &im2.samples(), 0.00001);
}

#[test]
fn isotope_model_offset() {
    // TOLERANCE_ABSOLUTE(0.1)
    let mut im1 = IsotopeModel::new().unwrap();
    im1.set_parameters(&isotope_tmp()).unwrap();
    im1.set_offset(673.5);
    let mut im2 = IsotopeModel::new().unwrap();
    im2.set_parameters(&im1.to_param().unwrap()).unwrap();
    im2.set_offset(673.5);
    samples_similar(&im1.samples(), &im2.samples(), 0.1);

    let mut im2 = IsotopeModel::new().unwrap();
    im2.set_parameters(&im1.to_param().unwrap()).unwrap();
    im2.set_offset(im1.offset());
    samples_similar(&im1.samples(), &im2.samples(), 0.1);
}

/// `setSamples(const EmpiricalFormula&)`, Gaussian peak shape.
const ISOTOPE_GAUSSIAN: [f64; 173] = [
    0.000429512,
    0.00093697,
    0.00196383,
    0.00395466,
    0.00765145,
    0.0142235,
    0.0254037,
    0.043593,
    0.0718726,
    0.113852,
    0.173278,
    0.253381,
    0.355987,
    0.480533,
    0.623218,
    0.776577,
    0.929731,
    1.06945,
    1.18192,
    1.25501,
    1.28036,
    1.25501,
    1.18192,
    1.06945,
    0.929731,
    0.776577,
    0.623218,
    0.480533,
    0.355987,
    0.253381,
    0.173278,
    0.113852,
    0.0718726,
    0.0439064,
    0.0260875,
    0.0156567,
    0.0105376,
    0.00953883,
    0.0123444,
    0.019477,
    0.0318149,
    0.0524539,
    0.0830909,
    0.126461,
    0.184922,
    0.259806,
    0.350701,
    0.454835,
    0.566759,
    0.678534,
    0.7805,
    0.862586,
    0.915925,
    0.934428,
    0.915925,
    0.862586,
    0.7805,
    0.678534,
    0.566759,
    0.454835,
    0.350701,
    0.259806,
    0.184922,
    0.126461,
    0.0830909,
    0.0524539,
    0.0318149,
    0.0186555,
    0.0106322,
    0.00611169,
    0.00394849,
    0.00348857,
    0.00450454,
    0.00682396,
    0.01171,
    0.0193065,
    0.0305829,
    0.0465459,
    0.0680634,
    0.0956255,
    0.129081,
    0.167409,
    0.208604,
    0.249745,
    0.287275,
    0.317488,
    0.33712,
    0.34393,
    0.33712,
    0.317488,
    0.287275,
    0.249745,
    0.208604,
    0.167409,
    0.129081,
    0.0956255,
    0.0680634,
    0.0465459,
    0.0305829,
    0.0193065,
    0.0117385,
    0.00688626,
    0.00395131,
    0.0023183,
    0.00157108,
    0.00147331,
    0.00194089,
    0.00289868,
    0.00477912,
    0.00757047,
    0.011522,
    0.0168484,
    0.0236711,
    0.0319527,
    0.0414404,
    0.0516379,
    0.0618218,
    0.071112,
    0.0785909,
    0.0834507,
    0.0851365,
    0.0834507,
    0.0785909,
    0.071112,
    0.0618218,
    0.0516379,
    0.0414404,
    0.0319527,
    0.0236711,
    0.0168484,
    0.011522,
    0.00757047,
    0.00477912,
    0.00290403,
    0.00170087,
    0.000970228,
    0.000558009,
    0.000358215,
    0.000307651,
    0.000378553,
    0.000542686,
    0.000894739,
    0.00141733,
    0.00215713,
    0.00315433,
    0.00443167,
    0.00598213,
    0.0077584,
    0.00966757,
    0.0115742,
    0.0133135,
    0.0147137,
    0.0156235,
    0.0159391,
    0.0156235,
    0.0147137,
    0.0133135,
    0.0115742,
    0.00966757,
    0.0077584,
    0.00598213,
    0.00443167,
    0.00315433,
    0.00215713,
    0.00141733,
    0.000894739,
    0.000542686,
    0.00031625,
    0.000177068,
    9.52526e-005,
    4.92314e-005,
    2.44476e-005,
    1.16643e-005,
];

/// `setSamples(const EmpiricalFormula&)`, Lorentzian peak shape (FWHM 0.05).
const ISOTOPE_LORENTZIAN: [f64; 160] = [
    0.0249619,
    0.0291547,
    0.0344977,
    0.0414526,
    0.0507371,
    0.0635168,
    0.0817848,
    0.109176,
    0.152888,
    0.228787,
    0.377365,
    0.725701,
    1.80202,
    5.53034,
    3.91171,
    1.28304,
    0.570747,
    0.315088,
    0.198245,
    0.1358,
    0.0986954,
    0.0749089,
    0.0587688,
    0.0473237,
    0.0389177,
    0.0325645,
    0.027647,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0182176,
    0.0212776,
    0.0251771,
    0.0302528,
    0.0370288,
    0.0463557,
    0.059688,
    0.0796788,
    0.111581,
    0.166973,
    0.275407,
    0.52963,
    1.31515,
    4.03614,
    2.85483,
    0.936385,
    0.416541,
    0.229957,
    0.144682,
    0.0991093,
    0.0720296,
    0.0546698,
    0.0428905,
    0.0345377,
    0.0284028,
    0.0237661,
    0.0201772,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.00670527,
    0.00783155,
    0.0092668,
    0.011135,
    0.013629,
    0.0170619,
    0.0219691,
    0.029327,
    0.041069,
    0.061457,
    0.101368,
    0.194938,
    0.48406,
    1.48556,
    1.05076,
    0.344651,
    0.153314,
    0.0846392,
    0.0532526,
    0.0364787,
    0.0265116,
    0.0201221,
    0.0157865,
    0.0127121,
    0.0104541,
    0.00874748,
    0.00742654,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.00165982,
    0.00193862,
    0.0022939,
    0.00275636,
    0.00337373,
    0.0042235,
    0.00543822,
    0.0072596,
    0.0101662,
    0.015213,
    0.0250926,
    0.048255,
    0.119824,
    0.367736,
    0.260106,
    0.0853148,
    0.0379514,
    0.0209516,
    0.0131821,
    0.00902993,
    0.00656267,
    0.00498101,
    0.00390778,
    0.00314675,
    0.00258781,
    0.00216535,
    0.00183836,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.000310749,
    0.000362945,
    0.000429461,
    0.000516041,
    0.000631624,
    0.000790718,
    0.00101813,
    0.00135913,
    0.0019033,
    0.00284816,
    0.00469779,
    0.00903422,
    0.0224333,
    0.068847,
    0.0486966,
    0.0159725,
    0.0071052,
    0.00392252,
    0.00246794,
    0.00169057,
    0.00122865,
    0.000932537,
    0.00073161,
    0.00058913,
    0.000484485,
    0.000405393,
    0.000344176,
];

#[test]
fn isotope_model_set_samples_from_formula() {
    let ef = EmpiricalFormula::parse("C66H129O3").unwrap();
    let mut tmp = param(&[
        ("statistics:mean", f(ef.average_mass() / 1.0)),
        ("interpolation_step", f(0.03)),
        ("charge", i(1)),
        ("isotope:maximum", i(100)),
        ("isotope:mode:mode", ParamValue::String("Gaussian".into())),
        ("isotope:mode:GaussianSD", f(0.15)),
    ]);
    let mut im1 = IsotopeModel::new().unwrap();
    im1.set_parameters(&tmp).unwrap();
    im1.set_samples_from_formula(&ef).unwrap();
    let samples = im1.samples();
    assert_eq!(samples.len(), ISOTOPE_GAUSSIAN.len());
    // The tolerance in force is TOLERANCE_ABSOLUTE(0.00001), the default.
    for (sample, &expected) in samples.iter().zip(&ISOTOPE_GAUSSIAN) {
        similar(f64::from(sample.intensity), expected);
    }

    tmp.set_value(
        "isotope:mode:mode",
        ParamValue::String("Lorentzian".into()),
        "",
        &[],
    )
    .unwrap();
    tmp.set_value("isotope:mode:LorentzFWHM", f(0.05), "", &[])
        .unwrap();
    im1.set_parameters(&tmp).unwrap();
    im1.set_samples_from_formula(&ef).unwrap();
    let samples = im1.samples();
    assert_eq!(samples.len(), ISOTOPE_LORENTZIAN.len());
    for (sample, &expected) in samples.iter().zip(&ISOTOPE_LORENTZIAN) {
        similar(f64::from(sample.intensity), expected);
    }
}

#[test]
fn isotope_model_set_samples_without_formula_is_not_implemented() {
    let mut im = IsotopeModel::new().unwrap();
    assert!(matches!(
        im.set_samples(),
        Err(openms::Error::Unsupported(_))
    ));
}

// ---------------------------------------------------------------------------
// ExtendedIsotopeModel_test.cpp

fn extended_tmp() -> Param {
    param(&[
        ("charge", i(3)),
        ("isotope:stdev", f(0.8)),
        ("isotope:monoisotopic_mz", f(670.5)),
    ])
}

#[test]
fn extended_isotope_model() {
    let im = ExtendedIsotopeModel::new().unwrap();
    assert_eq!(im.name(), "ExtendedIsotopeModel");
    assert_eq!(im.charge(), 1);
    // "default charge is 1 and hence center mus be 1"
    assert_eq!(im.center().unwrap(), 1.0);

    let mut im1 = ExtendedIsotopeModel::new().unwrap();
    im1.set_parameters(&extended_tmp()).unwrap();
    let im2 = im1.clone();
    let mut im3 = ExtendedIsotopeModel::new().unwrap();
    im3.set_parameters(&extended_tmp()).unwrap();
    params_equal(&im3.to_param().unwrap(), &im2.to_param().unwrap());

    // TOLERANCE_ABSOLUTE(0.00001)
    let mut im2 = ExtendedIsotopeModel::new().unwrap();
    im2.set_parameters(&im1.to_param().unwrap()).unwrap();
    samples_similar(&im1.samples(), &im2.samples(), 0.00001);

    // TOLERANCE_ABSOLUTE(0.1)
    let mut im1 = ExtendedIsotopeModel::new().unwrap();
    im1.set_parameters(&extended_tmp()).unwrap();
    im1.set_offset(673.5);
    let mut im2 = ExtendedIsotopeModel::new().unwrap();
    im2.set_parameters(&im1.to_param().unwrap()).unwrap();
    im2.set_offset(673.5);
    samples_similar(&im1.samples(), &im2.samples(), 0.1);
    let mut im2 = ExtendedIsotopeModel::new().unwrap();
    im2.set_parameters(&im1.to_param().unwrap()).unwrap();
    im2.set_offset(im1.offset());
    samples_similar(&im1.samples(), &im2.samples(), 0.1);
}

// ---------------------------------------------------------------------------
// ModelDescription_test.cpp, BaseModel_test.cpp

#[test]
fn model_description_round_trip() {
    let mut gm = GaussModel::new().unwrap();
    gm.set_parameters(&gauss_tmp(1.0, 3.0, 2.0, 0.5)).unwrap();
    let desc = ModelDescription::from_model(&gm).unwrap();
    assert_eq!(desc.name, "GaussModel");
    params_equal(&desc.param, &gm.to_param().unwrap());
    let copy = desc.clone();
    assert!(copy.source_equal(&desc).unwrap());
    let mut other = desc.clone();
    other.name = "BiGaussModel".into();
    assert!(!other.source_equal(&desc).unwrap());
    let rebuilt = desc.create_model().unwrap();
    assert_eq!(rebuilt.name(), "GaussModel");
    samples_similar(&rebuilt.samples(), &gm.samples(), 0.0);
    let unknown = ModelDescription {
        name: "NoSuchModel".into(),
        param: Param::new(),
    };
    assert!(unknown.create_model().is_err());
}

#[test]
fn base_model_cut_off_and_fill_intensities() {
    let mut gm = GaussModel::new().unwrap();
    gm.set_parameters(&gauss_tmp(-4.0, 4.0, 0.0, 1.0)).unwrap();
    assert_eq!(gm.cut_off(), 0.0);
    gm.set_cut_off(0.3);
    assert_eq!(gm.cut_off(), 0.3);
    assert_eq!(
        gm.to_param().unwrap().value("cutoff").unwrap(),
        &ParamValue::Float(0.3)
    );
    assert!(gm.is_contained(0.0));
    assert!(!gm.is_contained(3.0));
    let mut peaks = vec![Peak1D::new(0.0, 0.0), Peak1D::new(1.0, 0.0)];
    gm.fill_intensities(&mut peaks);
    assert_eq!(peaks[0].intensity, gm.intensity(0.0) as f32);
    assert_eq!(peaks[1].intensity, gm.intensity(1.0) as f32);
    let mut out = Vec::new();
    gm.write_samples(&mut out).unwrap();
    let text = String::from_utf8(out).unwrap();
    assert_eq!(text.lines().count(), gm.samples().len());
    assert!(text.starts_with("POS: -4 INT: "));
}

// ---------------------------------------------------------------------------
// Fitter1D_test.cpp, LevMarqFitter1D_test.cpp, MaxLikeliFitter1D_test.cpp

#[test]
fn fitter1d_base_is_not_implemented() {
    let base = BaseFitter1D::default();
    assert!(matches!(
        base.fit1d(&[]),
        Err(openms::Error::Unsupported(_))
    ));
    let copy = base.clone();
    assert_eq!(copy, base);
}

#[test]
fn lev_marq_optimize_refuses_too_few_points_and_budget() {
    let mut x = [1.0, 2.0];
    let err = lev_marq_optimize(&mut x, 1, |_, _| {}, |_, _| {}, 500).unwrap_err();
    assert!(err.to_string().contains("we always expects N>=p"));
    let err = lev_marq_optimize(&mut x, 3, |_, _| {}, |_, _| {}, 0).unwrap_err();
    assert!(err.to_string().contains("Error 0"));
    assert_eq!(x, [1.0, 2.0]);
    assert_eq!(MAX_ITERATION_DEFAULT, 500);
}

// ---------------------------------------------------------------------------
// GaussFitter1D_test.cpp, BiGaussFitter1D_test.cpp, IsotopeFitter1D_test.cpp,
// ExtendedIsotopeFitter1D_test.cpp

#[test]
fn gauss_fitter_parameters() {
    assert_eq!(GaussFitter1D::default().name(), "GaussFitter1D");
    let p = param(&[
        ("tolerance_stdev_bounding_box", f(1.0)),
        ("statistics:mean", f(680.1)),
        ("statistics:variance", f(2.0)),
        ("interpolation_step", f(1.0)),
    ]);
    let mut gf1 = GaussFitter1D::default();
    gf1.set_parameters(&p).unwrap();
    let gf2 = gf1.clone();
    let mut gf3 = GaussFitter1D::default();
    gf3.set_parameters(&p).unwrap();
    params_equal(&gf3.to_param().unwrap(), &gf2.to_param().unwrap());
}

#[test]
fn bigauss_fitter_parameters() {
    assert_eq!(BiGaussFitter1D::default().name(), "BiGaussFitter1D");
    let p = param(&[
        ("tolerance_stdev_bounding_box", f(1.0)),
        ("statistics:mean", f(680.1)),
        ("statistics:variance", f(2.0)),
        ("statistics:variance1", f(2.0)),
        ("statistics:variance2", f(5.0)),
        ("interpolation_step", f(1.0)),
    ]);
    let mut bgf1 = BiGaussFitter1D::default();
    bgf1.set_parameters(&p).unwrap();
    let bgf2 = bgf1.clone();
    let mut bgf3 = BiGaussFitter1D::default();
    bgf3.set_parameters(&p).unwrap();
    params_equal(&bgf3.to_param().unwrap(), &bgf2.to_param().unwrap());
}

#[test]
fn isotope_fitters_parameters() {
    assert_eq!(IsotopeFitter1D::default().name(), "IsotopeFitter1D");
    assert_eq!(
        ExtendedIsotopeFitter1D::default().name(),
        "ExtendedIsotopeFitter1D"
    );
    let p = param(&[
        ("tolerance_stdev_bounding_box", f(1.0)),
        ("statistics:mean", f(680.1)),
        ("statistics:variance", f(2.0)),
        ("interpolation_step", f(1.0)),
        ("charge", i(1)),
        ("isotope:stdev", f(0.04)),
        ("isotope:maximum", i(20)),
    ]);
    let mut isof1 = IsotopeFitter1D::default();
    isof1.set_parameters(&p).unwrap();
    let isof2 = isof1.clone();
    let mut isof3 = IsotopeFitter1D::default();
    isof3.set_parameters(&p).unwrap();
    params_equal(&isof3.to_param().unwrap(), &isof2.to_param().unwrap());

    let mut eisof1 = ExtendedIsotopeFitter1D::default();
    eisof1.set_parameters(&p).unwrap();
    let eisof2 = eisof1.clone();
    let mut eisof3 = ExtendedIsotopeFitter1D::default();
    eisof3.set_parameters(&p).unwrap();
    params_equal(&eisof3.to_param().unwrap(), &eisof2.to_param().unwrap());
}

// ---------------------------------------------------------------------------
// EmgFitter1D_test.cpp

#[test]
fn emg_fitter_parameters() {
    assert_eq!(EmgFitter1D::default().name(), "EmgFitter1D");
    let p = param(&[
        ("tolerance_stdev_bounding_box", f(1.0)),
        ("statistics:mean", f(680.1)),
        ("statistics:variance", f(2.0)),
        ("interpolation_step", f(1.0)),
        ("max_iteration", i(500)),
        ("deltaAbsError", f(0.0001)),
        ("deltaRelError", f(0.0001)),
    ]);
    let mut emgf1 = EmgFitter1D::default();
    let warnings = emgf1.set_parameters(&p).unwrap();
    assert_eq!(warnings.len(), 2, "{warnings:?}");
    let emgf2 = emgf1.clone();
    let mut emgf3 = EmgFitter1D::default();
    emgf3.set_parameters(&p).unwrap();
    params_equal(&emgf3.to_param().unwrap(), &emgf2.to_param().unwrap());
}

#[test]
fn emg_fitter_fits_exact_emg_samples() {
    let mut em = EmgModel::new().unwrap();
    em.set_interpolation_step(0.2);
    em.set_parameters(&emg_tmp(
        678.9, 789.0, 680.1, 2.0, 100000.0, 5.0, 5.0, 725.0,
    ))
    .unwrap();
    let samples = em.samples();
    let ef = EmgFitter1D::default();
    let (correlation, fitted) = ef.fit1d(&samples).unwrap();
    let p = fitted.to_param().unwrap();
    let value = |key: &str| p.value(key).unwrap().to_f64().unwrap();
    similar(f64::from(correlation), 1.0);
    similar(value("emg:height"), 100000.0);
    similar(value("emg:width"), 5.0);
    similar(value("emg:symmetry"), 5.0);
    similar(value("emg:retention"), 725.0);
    // The source's second half perturbs the samples with boost::mt19937 and
    // boost::normal_distribution; that generator pair is not ported, so the
    // perturbed input cannot be reproduced and that half is not transcribed.
}

#[test]
fn emg_fitter_robustness_issue_6239() {
    let peak = |mz: f64, intensity: f32| Peak1D::new(mz, intensity);
    let finite = |m: &Model1D, key: &str| {
        m.to_param()
            .unwrap()
            .value(key)
            .unwrap()
            .to_f64()
            .unwrap()
            .is_finite()
    };
    // (a) Fewer points than parameters.
    let few = [peak(10.0, 5.0), peak(11.0, 7.0), peak(12.0, 4.0)];
    let (q, m) = EmgFitter1D::default().fit1d(&few).unwrap();
    assert!(q <= 0.0);
    assert!(finite(&m, "emg:retention"));

    // (b) All-equal positions.
    let flat = vec![peak(50.0, 3.0); 8];
    let (q, _m) = EmgFitter1D::default().fit1d(&flat).unwrap();
    assert!(q <= 0.0);

    // (c) Right-truncated edge peak.
    let mut em = EmgModel::new().unwrap();
    em.set_interpolation_step(0.2);
    em.set_parameters(&emg_tmp(0.0, 200.0, 180.0, 2.0, 100000.0, 5.0, 5.0, 180.0))
        .unwrap();
    let trunc: Vec<Peak1D> = em.samples().into_iter().filter(|p| p.mz <= 199.0).collect();
    assert!(trunc.len() >= 4);
    let (q, m) = EmgFitter1D::default().fit1d(&trunc).unwrap();
    assert!(q.is_finite());
    assert!(finite(&m, "emg:retention"));
    assert!(finite(&m, "emg:width"));

    // (d) Non-finite leading position.
    let bad = [
        peak(f64::NAN, 5.0),
        peak(11.0, 7.0),
        peak(12.0, 4.0),
        peak(13.0, 3.0),
    ];
    let (q, m) = EmgFitter1D::default().fit1d(&bad).unwrap();
    similar(f64::from(q), -1.0);
    for key in [
        "emg:height",
        "emg:width",
        "emg:symmetry",
        "emg:retention",
        "bounding_box:min",
        "bounding_box:max",
    ] {
        assert!(finite(&m, key), "{key}");
    }
}

// ---------------------------------------------------------------------------
// EmgScoring_test.cpp, with OpenSwathTestHelper.h's mock feature and group

const MOCK_RT: [f64; 11] = [
    3103.13, 3106.56, 3109.98, 3113.41, 3116.84, 3120.26, 3123.69, 3127.11, 3130.54, 3133.97,
    3137.4,
];
const MOCK_TR3: [f64; 11] = [
    5.97544, 4.27492, 3.33018, 4.08597, 5.50307, 5.24327, 8.40812, 2.8342, 6.94379, 7.69957,
    4.08597,
];
const MOCK_TR1: [f64; 11] = [
    15.8951, 41.5446, 76.0746, 109.069, 111.904, 169.792, 121.044, 63.0137, 44.615, 21.4927,
    7.93576,
];
const MOCK_TR5: [f64; 11] = [
    5.73925, 6.7076, 2.85782, 5.0307, 8.95135, 14.4544, 20.9731, 24.3033, 20.6897, 13.7459, 8.90411,
];

fn hull_points(intensities: &[f64; 11]) -> Vec<Point2D> {
    MOCK_RT
        .iter()
        .zip(intensities)
        .map(|(&rt, &intensity)| Point2D { rt, mz: intensity })
        .collect()
}

fn mock_feature() -> MRMFeature {
    let mut feature = MRMFeature::new();
    for (key, intensities) in [("tr3", &MOCK_TR3), ("tr1", &MOCK_TR1), ("tr5", &MOCK_TR5)] {
        let mut hull = ConvexHull2D::new();
        hull.set_hull_points(&hull_points(intensities)).unwrap();
        let mut f = Feature::default();
        f.convex_hulls.push(hull);
        feature.add_feature(f, key).unwrap();
    }
    feature
}

#[test]
fn emg_scoring() {
    let scoring = EmgScoring::new();
    let defaults = scoring.defaults().unwrap();
    params_equal(&defaults, &EmgFitter1D::defaults().unwrap());
    let mut scoring2 = EmgScoring::new();
    scoring2.set_fitter_param(defaults);

    let score = scoring
        .elution_model_fit(&hull_points(&MOCK_TR1), false)
        .unwrap();
    similar(score, 0.981013417243958);

    let feature = mock_feature();
    let mut group: MRMTransitionGroup<MSChromatogram, ()> = MRMTransitionGroup::new();
    for key in ["tr3", "tr1", "tr5"] {
        let chromatogram = MSChromatogram {
            native_id: key.into(),
            ..MSChromatogram::default()
        };
        group.add_chromatogram(chromatogram, key).unwrap();
    }
    let score = scoring.calc_elution_fit_score(&feature, &group).unwrap();
    similar(score, 0.924365639);

    assert_eq!(
        scoring
            .elution_model_fit(&hull_points(&MOCK_TR1)[..1], false)
            .unwrap(),
        -1.0
    );
}

// ---------------------------------------------------------------------------
// FeatureFinderAlgorithm.h (header only, no class test)

struct NoSeeds;

impl FeatureFinderAlgorithm for NoSeeds {
    fn default_parameters(&self) -> openms::Result<Param> {
        Ok(Param::new())
    }
    fn run(&mut self, _map: &MSExperiment, _features: &mut FeatureMap) -> openms::Result<()> {
        Ok(())
    }
}

#[test]
fn feature_finder_algorithm_defaults() {
    let summary = Summary::default();
    assert_eq!(summary.no_exceptions, 0);
    assert_eq!(summary.corr_mean, 0.0);
    assert_eq!(summary.corr_max, 0.0);
    assert_eq!(summary.corr_min, 1.0);

    let mut summary = Summary::default();
    summary.count_mz_stdev(0.5).unwrap();
    summary.count_mz_stdev(0.1).unwrap();
    summary.count_mz_stdev(0.5).unwrap();
    assert_eq!(summary.mz_stdev, vec![(0.1, 1), (0.5, 2)]);
    assert!(summary.count_mz_stdev(f32::NAN).is_err());

    let mut algorithm = NoSeeds;
    algorithm.set_seeds(&FeatureMap::default()).unwrap();
    let seeds = FeatureMap::from_features(vec![Feature::default()]);
    let err = algorithm.set_seeds(&seeds).unwrap_err();
    assert!(err.to_string().contains(SEEDS_NOT_SUPPORTED));
    algorithm
        .run(&MSExperiment::default(), &mut FeatureMap::default())
        .unwrap();
}

// ---------------------------------------------------------------------------
// SeedListGenerator_test.cpp

fn pos(rt: f64, mz: f64) -> DPosition2 {
    DPosition2 {
        coordinates: [rt, mz],
    }
}

/// The 13 spectra of the upstream `PepXMLFile_test.mzML` (core `bc9cc12`,
/// sha256 `f59715fe...bbf51d`): MS level (MS:1000511), scan start time
/// (MS:1000016, seconds) and selected-ion m/z (MS:1000744). The file itself is
/// 1.4 MB and only these fields reach `generateSeedList`.
const PEPXML_TEST_MZML: [(u32, f64, Option<f64>); 13] = [
    (1, 0.5927, None),
    (2, 1.3653, Some(538.605)),
    (2, 1.719, Some(637.885)),
    (2, 2.0436, Some(678.384)),
    (1, 2.3364, None),
    (2, 2.7843, Some(762.728)),
    (2, 3.085, Some(1143.59)),
    (2, 3.4018, Some(770.055)),
    (1, 3.7572, None),
    (2, 4.4573, Some(757.412)),
    (2, 4.8288, Some(807.394)),
    (2, 5.1785, Some(512.784)),
    (1, 5.4999, None),
];

#[test]
fn seed_list_from_experiment() {
    let mut experiment = MSExperiment::default();
    for (ms_level, rt, precursor) in PEPXML_TEST_MZML {
        let mut spectrum = MSSpectrum {
            ms_level,
            rt,
            ..MSSpectrum::default()
        };
        if let Some(mz) = precursor {
            spectrum.precursors.push(Precursor {
                mz,
                ..Precursor::default()
            });
        }
        experiment.spectra.push(spectrum);
    }
    let seeds = SeedListGenerator::new()
        .generate_seed_list_from_experiment(&experiment)
        .unwrap();
    assert_eq!(seeds.len(), 9);
    assert_eq!(seeds[0], pos(0.5927, 538.605));
    assert_eq!(seeds[1], pos(0.5927, 637.885));
    assert_eq!(seeds[2], pos(0.5927, 678.384));
    assert_eq!(seeds[8], pos(3.7572, 512.784));

    // "MS2 spectrum without a precursor spectrum / precursor info must be
    // skipped, not crash".
    let mut no_prec = MSExperiment::default();
    no_prec.spectra.push(MSSpectrum {
        ms_level: 2,
        ..MSSpectrum::default()
    });
    let seeds = SeedListGenerator::new()
        .generate_seed_list_from_experiment(&no_prec)
        .unwrap();
    assert!(seeds.is_empty());
}

#[test]
fn seed_list_from_peptides() {
    let mut peptides = vec![PeptideIdentification::default(); 3];
    for (k, (rt, mz)) in [(1.1, 111.111), (2.2, 222.222), (3.3, 333.333)]
        .into_iter()
        .enumerate()
    {
        peptides[k].rt = Some(rt);
        peptides[k].mz = Some(mz);
    }
    let seeds = SeedListGenerator::new()
        .generate_seed_list_from_peptides(&mut peptides, false)
        .unwrap();
    assert_eq!(seeds.len(), 3);
    assert_eq!(seeds[0], pos(1.1, 111.111));
    assert_eq!(seeds[1], pos(2.2, 222.222));
    assert_eq!(seeds[2], pos(3.3, 333.333));

    let hit = PeptideHit {
        sequence: AASequence::parse("TEST").unwrap(),
        charge: 2,
        ..PeptideHit::default()
    };
    peptides[0].hits.push(hit);
    peptides.truncate(1);
    let seeds = SeedListGenerator::new()
        .generate_seed_list_from_peptides(&mut peptides, true)
        .unwrap();
    assert_eq!(seeds.len(), 1);
    similar(seeds[0].coordinates[1], 219.09755);
}

#[test]
fn seed_list_conversions() {
    let seeds = vec![pos(1.1, 111.111), pos(2.2, 222.222), pos(3.3, 333.333)];
    let features = SeedListGenerator::new().seeds_to_features(&seeds).unwrap();
    assert_eq!(features.features.len(), 3);
    for (k, s) in seeds.iter().enumerate() {
        assert_eq!(features.features[k].base.rt, s.coordinates[0]);
        assert_eq!(features.features[k].base.mz, s.coordinates[1]);
        assert_eq!(features.features[k].base.unique_id, k as u64);
    }

    let mut map = FeatureMap::default();
    for s in &seeds {
        let mut feature = Feature::default();
        feature.base.rt = s.coordinates[0];
        feature.base.mz = s.coordinates[1];
        map.features.push(feature);
    }
    let back = SeedListGenerator::new().features_to_seeds(&map);
    assert_eq!(back, seeds);
}
