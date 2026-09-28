// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The feature finder's models, fitters and EmgScoring against the executed
//! C++ (tier 1).
//!
//! `data/featurefinder_models/oracle.tsv` is the output of the driver
//! `../oracle/FF1-models/driver.cpp`, built with libOpenMS's core flags against
//! the Linux x86_64 **Release** install of OpenMS (core `bc9cc12`, cli
//! `c19e494`, topp `174b576`; conda gcc 14.4.0, `-O3 -DNDEBUG`) and run on
//! **ibminode06** (Ubuntu 24.04, glibc 2.39) twice in fresh processes with
//! byte-identical output. Each record is `<case>\t<key>\t<%.17g value>`. The
//! cases below rebuild each driver case through the Rust API.
//!
//! Every value is compared within a relative 1e-9 (absolute 1e-12 near zero)
//! and the number of bit-identical values is asserted per case where the port
//! computes in the source's order with the reference build's `exp`
//! (`feature_finder_picked::glibc_libm::exp`).

use openms::analysis::featurefinder::base_model::{BaseModel, InterpolationModel, Model1D};
use openms::analysis::featurefinder::emg_scoring::EmgScoring;
use openms::analysis::featurefinder::fitter1d::Fitter1D;
use openms::analysis::featurefinder::fitters::{
    BiGaussFitter1D, EmgFitter1D, ExtendedIsotopeFitter1D, GaussFitter1D, IsotopeFitter1D,
};
use openms::analysis::featurefinder::gauss_models::{BiGaussModel, EmgModel, GaussModel};
use openms::analysis::featurefinder::isotope_model::{ExtendedIsotopeModel, IsotopeModel};
use openms::chemistry::EmpiricalFormula;
use openms::kernel::{Peak1D, Point2D};
use openms::math::ml::linear_interpolation::LinearInterpolation;
use openms::param::{Param, ParamValue};
use std::collections::BTreeMap;

const ORACLE: &str = include_str!("data/featurefinder_models/oracle.tsv");

fn oracle() -> BTreeMap<(String, String), f64> {
    ORACLE
        .lines()
        .map(|line| {
            let mut fields = line.split('\t');
            let case = fields.next().unwrap().to_owned();
            let key = fields.next().unwrap().to_owned();
            let value = fields.next().unwrap();
            let value = match value {
                "nan" | "-nan" => f64::NAN,
                "inf" => f64::INFINITY,
                "-inf" => f64::NEG_INFINITY,
                v => v.parse().unwrap(),
            };
            ((case, key), value)
        })
        .collect()
}

/// Collects `(key, value)` pairs of one case and compares them with the
/// oracle's, reporting every difference at once.
struct Case {
    name: &'static str,
    values: Vec<(String, f64)>,
}

impl Case {
    fn new(name: &'static str) -> Self {
        Self {
            name,
            values: Vec::new(),
        }
    }
    fn out(&mut self, key: impl Into<String>, value: f64) {
        self.values.push((key.into(), value));
    }
    fn samples<M: InterpolationModel + ?Sized>(&mut self, model: &M) {
        let s = model.samples();
        self.out("n", s.len() as f64);
        for (i, p) in s.iter().enumerate() {
            self.out(format!("pos{i}"), p.mz);
            self.out(format!("int{i}"), f64::from(p.intensity));
        }
        self.out("offset", model.interpolation().offset());
        self.out("scale", model.interpolation().scale());
        for (i, v) in model.interpolation().data().iter().enumerate() {
            self.out(format!("raw{i}"), *v);
        }
    }
    /// Compare with the oracle; returns the number of bit-identical values.
    #[track_caller]
    fn check(self, oracle: &BTreeMap<(String, String), f64>) -> usize {
        let expected: BTreeMap<&str, f64> = oracle
            .iter()
            .filter(|((c, _), _)| c == self.name)
            .map(|((_, k), v)| (k.as_str(), *v))
            .collect();
        assert!(!expected.is_empty(), "no oracle records for {}", self.name);
        let mut failures = Vec::new();
        let mut identical = 0usize;
        let mut seen = 0usize;
        for (key, actual) in &self.values {
            let Some(&want) = expected.get(key.as_str()) else {
                continue;
            };
            seen += 1;
            if actual.to_bits() == want.to_bits() || (actual.is_nan() && want.is_nan()) {
                identical += 1;
                continue;
            }
            let diff = (actual - want).abs();
            if !(diff <= 1e-12 || diff <= 1e-9 * want.abs()) {
                failures.push(format!("{key}: rust {actual:e} oracle {want:e}"));
            }
        }
        assert_eq!(
            seen,
            expected.len(),
            "{}: {} oracle keys were not produced",
            self.name,
            expected.len() - seen
        );
        eprintln!("{}: {identical} of {seen} values bit-identical", self.name);
        assert!(
            failures.is_empty(),
            "{}: {} of {} values differ:\n{}",
            self.name,
            failures.len(),
            seen,
            failures[..failures.len().min(20)].join("\n")
        );
        identical
    }
}

fn param(entries: &[(&str, ParamValue)]) -> Param {
    let mut p = Param::new();
    for (key, value) in entries {
        p.set_value(key, value.clone(), "", &[]).unwrap();
    }
    p
}

fn f(v: f64) -> ParamValue {
    ParamValue::Float(v)
}

fn i(v: i64) -> ParamValue {
    ParamValue::Integer(v)
}

fn fit_case<F: Fitter1D>(name: &'static str, fitter: &F, data: &[Peak1D]) -> Case {
    let (quality, model) = fitter.fit1d(data).unwrap();
    let mut case = Case::new(name);
    case.out("quality", f64::from(quality));
    case.out("center", model.center().unwrap());
    let p = model.to_param().unwrap();
    for key in [
        "cutoff",
        "interpolation_step",
        "intensity_scaling",
        "bounding_box:min",
        "bounding_box:max",
        "statistics:mean",
        "statistics:variance",
        "statistics:variance1",
        "statistics:variance2",
        "emg:height",
        "emg:width",
        "emg:symmetry",
        "emg:retention",
        "charge",
        "isotope:maximum",
        "isotope:mode:GaussianSD",
        "isotope:mode:LorentzFWHM",
        "isotope:trim_right_cutoff",
        "isotope:distance",
        "isotope:stdev",
        "isotope:monoisotopic_mz",
        "averagines:C",
        "averagines:H",
        "averagines:N",
        "averagines:O",
        "averagines:S",
    ] {
        if let Ok(value) = p.value(key) {
            case.out(format!("param:{key}"), value.to_f64().unwrap());
        }
    }
    case.samples(&model);
    case
}

fn gauss_data(lo: f64, hi: f64, step: f64, mu: f64, var: f64, h: f64) -> Vec<Peak1D> {
    let mut d = Vec::new();
    let mut k = 0;
    while lo + f64::from(k) * step <= hi {
        let x = lo + f64::from(k) * step;
        d.push(Peak1D::new(
            x,
            (h * (-(x - mu) * (x - mu) / (2.0 * var)).exp()) as f32,
        ));
        k += 1;
    }
    d
}

#[test]
fn linear_interpolation_matches_the_release_build() {
    let oracle = oracle();
    let mut li = LinearInterpolation::new(0.5, -1.25);
    li.set_data(vec![1.0, 2.0, 0.0, 1.0, 3.5, -2.0]);
    let mut case = Case::new("li");
    for k in -20..=60 {
        let x = f64::from(k) / 8.0;
        case.out(format!("value{k}"), li.value(x));
        case.out(format!("deriv{k}"), li.derivative(x));
    }
    assert_eq!(case.check(&oracle), 162);

    let mut add = LinearInterpolation::default();
    add.data_mut().resize(7, 0.0);
    add.set_mapping(0.3, 1.0, 2.0);
    for k in -10..=40 {
        add.add_value(f64::from(k) / 7.0, 1.0 + f64::from(k));
    }
    let mut case = Case::new("li_add");
    for (k, v) in add.data().iter().enumerate() {
        case.out(format!("d{k}"), *v);
    }
    assert_eq!(case.check(&oracle), 7);
}

#[test]
fn models_match_the_release_build() {
    let oracle = oracle();
    let mut g = GaussModel::new().unwrap();
    g.set_parameters(&param(&[
        ("bounding_box:min", f(-4.0)),
        ("bounding_box:max", f(4.0)),
        ("statistics:mean", f(0.3)),
        ("statistics:variance", f(1.7)),
        ("intensity_scaling", f(2.5)),
    ]))
    .unwrap();
    let mut case = Case::new("gauss");
    case.samples(&g);
    case.check(&oracle);
    g.set_offset(1.5);
    let mut case = Case::new("gauss_off");
    case.out("center", g.center().unwrap());
    let p = g.to_param().unwrap();
    for key in [
        "cutoff",
        "interpolation_step",
        "intensity_scaling",
        "bounding_box:min",
        "bounding_box:max",
        "statistics:mean",
        "statistics:variance",
    ] {
        case.out(
            format!("param:{key}"),
            p.value(key).unwrap().to_f64().unwrap(),
        );
    }
    for k in -30..=50 {
        case.out(format!("at{k}"), g.intensity(f64::from(k) / 5.0));
    }
    case.check(&oracle);

    let mut b = BiGaussModel::new().unwrap();
    b.set_parameters(&param(&[
        ("bounding_box:min", f(678.9)),
        ("bounding_box:max", f(689.0)),
        ("statistics:mean", f(680.1)),
        ("statistics:variance1", f(2.0)),
        ("statistics:variance2", f(5.0)),
    ]))
    .unwrap();
    let mut case = Case::new("bigauss");
    case.samples(&b);
    case.check(&oracle);

    let mut e = EmgModel::new().unwrap();
    e.set_parameters(&param(&[
        ("bounding_box:min", f(-1.0)),
        ("bounding_box:max", f(4.0)),
        ("statistics:mean", f(0.0)),
        ("statistics:variance", f(0.1)),
        ("emg:height", f(10.0)),
        ("emg:width", f(1.0)),
        ("emg:symmetry", f(2.0)),
        ("emg:retention", f(3.0)),
    ]))
    .unwrap();
    let mut case = Case::new("emg_small");
    case.samples(&e);
    // The EMG samples use the reference build's exp: every value is identical.
    assert_eq!(case.check(&oracle), 156);
    e.set_parameters(&param(&[
        ("bounding_box:min", f(678.9)),
        ("bounding_box:max", f(789.0)),
        ("statistics:mean", f(680.1)),
        ("statistics:variance", f(2.0)),
        ("emg:height", f(100000.0)),
        ("emg:width", f(5.0)),
        ("emg:symmetry", f(5.0)),
        ("emg:retention", f(725.0)),
    ]))
    .unwrap();
    let mut case = Case::new("emg_class");
    case.samples(&e);
    assert_eq!(case.check(&oracle), 3309);
}

#[test]
fn isotope_models_match_the_release_build() {
    let oracle = oracle();
    let ef = EmpiricalFormula::parse("C66H129O3").unwrap();
    let mut p = param(&[
        ("statistics:mean", f(ef.average_mass() / 1.0)),
        ("interpolation_step", f(0.03)),
        ("charge", i(1)),
        ("isotope:maximum", i(100)),
        ("isotope:mode:mode", ParamValue::String("Gaussian".into())),
        ("isotope:mode:GaussianSD", f(0.15)),
    ]);
    let mut im = IsotopeModel::new().unwrap();
    im.set_parameters(&p).unwrap();
    im.set_samples_from_formula(&ef).unwrap();
    let mut case = Case::new("iso_gauss");
    case.out("center", im.center().unwrap());
    case.samples(&im);
    case.check(&oracle);

    p.set_value(
        "isotope:mode:mode",
        ParamValue::String("Lorentzian".into()),
        "",
        &[],
    )
    .unwrap();
    p.set_value("isotope:mode:LorentzFWHM", f(0.05), "", &[])
        .unwrap();
    im.set_parameters(&p).unwrap();
    im.set_samples_from_formula(&ef).unwrap();
    let mut case = Case::new("iso_lorentz");
    case.out("center", im.center().unwrap());
    case.samples(&im);
    case.check(&oracle);

    let mut im2 = IsotopeModel::new().unwrap();
    im2.set_parameters(&param(&[
        ("statistics:mean", f(670.5)),
        ("charge", i(3)),
        ("isotope:mode:GaussianSD", f(0.05)),
    ]))
    .unwrap();
    let formula = im2.formula().unwrap();
    im2.set_samples_from_formula(&formula).unwrap();
    let mut case = Case::new("iso_avg");
    case.out("formula_C", f64::from(formula.count("C").unwrap()));
    case.out("formula_H", f64::from(formula.count("H").unwrap()));
    case.out("center", im2.center().unwrap());
    case.samples(&im2);
    case.check(&oracle);

    let e0 = ExtendedIsotopeModel::new().unwrap();
    let mut case = Case::new("ext_default");
    case.samples(&e0);
    case.check(&oracle);
    let mut e = ExtendedIsotopeModel::new().unwrap();
    e.set_parameters(&param(&[
        ("charge", i(3)),
        ("isotope:stdev", f(0.08)),
        ("isotope:monoisotopic_mz", f(670.5)),
    ]))
    .unwrap();
    let mut case = Case::new("ext");
    case.samples(&e);
    case.check(&oracle);
}

#[test]
fn fitters_match_the_release_build() {
    let oracle = oracle();
    let d = gauss_data(675.0, 685.0, 0.25, 680.37, 2.0, 1000.0);
    let gf = GaussFitter1D::from_param(&param(&[
        ("statistics:mean", f(680.0)),
        ("statistics:variance", f(2.0)),
    ]))
    .unwrap();
    fit_case("fit_gauss", &gf, &d).check(&oracle);
    let bf = BiGaussFitter1D::from_param(&param(&[
        ("statistics:mean", f(680.0)),
        ("statistics:variance1", f(1.5)),
        ("statistics:variance2", f(3.0)),
    ]))
    .unwrap();
    fit_case("fit_bigauss", &bf, &d).check(&oracle);

    let mut em = EmgModel::new().unwrap();
    em.set_interpolation_step(0.2);
    em.set_parameters(&param(&[
        ("bounding_box:min", f(678.9)),
        ("bounding_box:max", f(789.0)),
        ("statistics:mean", f(680.1)),
        ("statistics:variance", f(2.0)),
        ("emg:height", f(100000.0)),
        ("emg:width", f(5.0)),
        ("emg:symmetry", f(5.0)),
        ("emg:retention", f(725.0)),
    ]))
    .unwrap();
    let es = em.samples();
    fit_case("fit_emg", &EmgFitter1D::default(), &es).check(&oracle);
    let efm = EmgFitter1D::from_param(&param(&[("init_mom", ParamValue::String("true".into()))]))
        .unwrap();
    fit_case("fit_emg_mom", &efm, &es).check(&oracle);
    let skew: Vec<Peak1D> = (0..40)
        .map(|k| {
            let x = 100.0 + f64::from(k) * 0.7;
            let tail = if x > 110.0 {
                1500.0 * (-(x - 116.0) / 4.0).exp()
            } else {
                0.0
            };
            let y = 5000.0 * (-(x - 110.0) * (x - 110.0) / 18.0).exp() + tail;
            Peak1D::new(x, y as f32)
        })
        .collect();
    let ef2 = EmgFitter1D::from_param(&param(&[("max_iteration", i(37))])).unwrap();
    fit_case("fit_emg_skew", &ef2, &skew).check(&oracle);

    let mut src = IsotopeModel::new().unwrap();
    src.set_parameters(&param(&[
        ("statistics:mean", f(700.35)),
        ("charge", i(2)),
        ("isotope:mode:GaussianSD", f(0.1)),
    ]))
    .unwrap();
    let formula = src.formula().unwrap();
    src.set_samples_from_formula(&formula).unwrap();
    let isd_sub: Vec<Peak1D> = src.samples().into_iter().step_by(3).collect();
    let isf = IsotopeFitter1D::from_param(&param(&[
        ("statistics:mean", f(700.3)),
        ("charge", i(2)),
        ("isotope:stdev", f(0.1)),
    ]))
    .unwrap();
    fit_case("fit_iso", &isf, &isd_sub).check(&oracle);
    let isf0 =
        IsotopeFitter1D::from_param(&param(&[("statistics:mean", f(700.3)), ("charge", i(0))]))
            .unwrap();
    fit_case("fit_iso_c0", &isf0, &isd_sub).check(&oracle);
    let xf = ExtendedIsotopeFitter1D::from_param(&param(&[
        ("isotope:monoisotopic_mz", f(699.9)),
        ("charge", i(2)),
        ("isotope:stdev", f(0.1)),
    ]))
    .unwrap();
    let case = fit_case("fit_ext", &xf, &isd_sub);
    assert!(matches!(
        xf.fit1d(&isd_sub).unwrap().1,
        Model1D::ExtendedIsotope(_)
    ));
    case.check(&oracle);
}

#[test]
fn emg_scoring_matches_the_release_build() {
    let oracle = oracle();
    let rt = [
        3103.13, 3106.56, 3109.98, 3113.41, 3116.84, 3120.26, 3123.69, 3127.11, 3130.54, 3133.97,
        3137.4,
    ];
    let traces: [(&str, [f64; 11]); 3] = [
        (
            "tr3",
            [
                5.97544, 4.27492, 3.33018, 4.08597, 5.50307, 5.24327, 8.40812, 2.8342, 6.94379,
                7.69957, 4.08597,
            ],
        ),
        (
            "tr1",
            [
                15.8951, 41.5446, 76.0746, 109.069, 111.904, 169.792, 121.044, 63.0137, 44.615,
                21.4927, 7.93576,
            ],
        ),
        (
            "tr5",
            [
                5.73925, 6.7076, 2.85782, 5.0307, 8.95135, 14.4544, 20.9731, 24.3033, 20.6897,
                13.7459, 8.90411,
            ],
        ),
    ];
    let scoring = EmgScoring::new();
    let mut case = Case::new("emg_scoring");
    for (name, intensities) in traces {
        let points: Vec<Point2D> = rt
            .iter()
            .zip(intensities)
            .map(|(&rt, mz)| Point2D { rt, mz })
            .collect();
        case.out(
            format!("{name}_raw"),
            scoring.elution_model_fit(&points, false).unwrap(),
        );
        case.out(
            format!("{name}_smooth"),
            scoring.elution_model_fit(&points, true).unwrap(),
        );
    }
    case.check(&oracle);
}
