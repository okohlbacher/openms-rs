// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `ElutionModelFitter`: the upstream class test (core `bc9cc12`,
//! `ElutionModelFitter_test.cpp` with its `ElutionModelFitter_test.featureXML`,
//! byte for byte) and the executed Release build.
//!
//! The `emf_*` records of `data/featurefinder_models/oracle.tsv` come from
//! `../oracle/FF1-models/driver.cpp` on **ibminode06** against the Linux x86_64
//! **Release** install `openms4-release-bc9cc12-c19e494-174b576`, run twice
//! with identical output: every feature's (and subordinate's) intensity, model
//! meta values and status digit after `fitElutionModels` with the defaults,
//! `asymmetric`, `each_trace` + `unweighted_fit` + `add_zeros 0`, and
//! `asymmetric` + `no_imputation` + tight width/asymmetry checks.
#![cfg(feature = "featurexml")]

use openms::analysis::featurefinder::elution_model_fitter::ElutionModelFitter;
use openms::format::featurexml;
use openms::kernel::features::{Feature, FeatureMap};
use openms::param::{Param, ParamValue};
use std::collections::BTreeMap;

const ORACLE: &str = include_str!("data/featurefinder_models/oracle.tsv");
const FIXTURE: &[u8] =
    include_bytes!("data/featurefinder_models/ElutionModelFitter_test.featureXML");

fn load() -> FeatureMap {
    featurexml::read(FIXTURE).unwrap()
}

fn oracle(case: &str) -> BTreeMap<String, f64> {
    ORACLE
        .lines()
        .filter_map(|line| {
            let mut fields = line.split('\t');
            (fields.next()? == case).then(|| {
                let key = fields.next().unwrap().to_owned();
                let value = match fields.next().unwrap() {
                    "nan" | "-nan" => f64::NAN,
                    "inf" => f64::INFINITY,
                    "-inf" => f64::NEG_INFINITY,
                    v => v.parse().unwrap(),
                };
                (key, value)
            })
        })
        .collect()
}

const KEYS: [&str; 13] = [
    "model_height",
    "model_FWHM",
    "model_center",
    "model_lower",
    "model_upper",
    "model_EGH_tau",
    "model_EGH_sigma",
    "model_width",
    "model_asymmetry",
    "model_Gauss_sigma",
    "model_error",
    "model_area",
    "raw_intensity",
];

fn dump(out: &mut BTreeMap<String, f64>, f: &Feature, prefix: &str) {
    out.insert(format!("{prefix}intensity"), f64::from(f.base.intensity));
    for key in KEYS {
        if let Some(v) = f.base.metadata.get(key) {
            out.insert(format!("{prefix}{key}"), v.as_f64().unwrap());
        }
    }
    if let Some(v) = f.base.metadata.get("model_status") {
        let digit = v.as_str().unwrap().as_bytes()[0] - b'0';
        out.insert(format!("{prefix}status"), f64::from(digit));
    }
}

fn run_case(case: &str, param: &Param) -> usize {
    let mut fm = load();
    let fitter = ElutionModelFitter::from_param(param).unwrap();
    fitter.fit_elution_models(&mut fm).unwrap();
    let mut actual = BTreeMap::new();
    actual.insert("n".to_owned(), fm.features.len() as f64);
    for (i, f) in fm.features.iter().enumerate() {
        let prefix = format!("f{i}:");
        dump(&mut actual, f, &prefix);
        for (j, s) in f.subordinates.iter().enumerate() {
            dump(&mut actual, s, &format!("{prefix}s{j}:"));
        }
    }
    let expected = oracle(case);
    assert!(!expected.is_empty());
    let actual_keys: Vec<&String> = actual.keys().collect();
    let expected_keys: Vec<&String> = expected.keys().collect();
    assert_eq!(actual_keys, expected_keys, "{case}: key sets differ");
    let mut identical = 0;
    let mut failures = Vec::new();
    for (key, want) in &expected {
        let got = actual[key];
        if got.to_bits() == want.to_bits() || (got.is_nan() && want.is_nan()) {
            identical += 1;
        } else if !((got - want).abs() <= 1e-12 || (got - want).abs() <= 1e-9 * want.abs()) {
            failures.push(format!("{key}: rust {got:e} oracle {want:e}"));
        }
    }
    eprintln!(
        "{case}: {identical} of {} values bit-identical",
        expected.len()
    );
    // Every function on this path is the reference build's (glibc_libm, the
    // ported logf) or IEEE arithmetic in the source's order, so the values
    // are asserted identical, not merely close, on every host.
    assert_eq!(identical, expected.len(), "{case}: not bit-identical");
    assert!(
        failures.is_empty(),
        "{case}: {} values differ:\n{}",
        failures.len(),
        failures[..failures.len().min(30)].join("\n")
    );
    identical
}

fn flag(key: &str, p: &mut Param) {
    p.set_value(key, ParamValue::String("true".into()), "", &[])
        .unwrap();
}

#[test]
fn elution_model_fitter_class_test() {
    let emf = ElutionModelFitter::new();
    // "test if exception is thrown on empty featuremap"
    assert!(matches!(
        emf.fit_elution_models(&mut FeatureMap::default()),
        Err(openms::Error::MissingInformation(_))
    ));
    let mut features = load();
    assert_eq!(features.features.len(), 25);
    // symmetric model (default)
    emf.fit_elution_models(&mut features).unwrap();
    assert_eq!(features.features.len(), 25);
    for f in &features.features {
        let meta = &f.base.metadata;
        assert!(meta.contains_key("model_area"));
        assert!(meta.contains_key("model_status"));
        assert!(meta.contains_key("raw_intensity"));
        assert_ne!(
            f64::from(f.base.intensity),
            meta["raw_intensity"].as_f64().unwrap()
        );
        assert!(meta.contains_key("model_Gauss_sigma"));
    }
    // asymmetric model
    let mut features = load();
    let mut p = Param::new();
    flag("asymmetric", &mut p);
    let mut emf = ElutionModelFitter::new();
    emf.set_parameters(&p).unwrap();
    emf.fit_elution_models(&mut features).unwrap();
    assert_eq!(features.features.len(), 25);
    for f in &features.features {
        let meta = &f.base.metadata;
        assert!(meta.contains_key("model_area"));
        assert!(meta.contains_key("model_status"));
        assert!(meta.contains_key("raw_intensity"));
        assert_ne!(
            f64::from(f.base.intensity),
            meta["raw_intensity"].as_f64().unwrap()
        );
        assert!(meta.contains_key("model_EGH_tau"));
        assert!(meta.contains_key("model_EGH_sigma"));
    }
}

#[test]
fn elution_model_fitter_matches_the_release_build() {
    run_case("emf_sym", &Param::new());
    let mut a = Param::new();
    flag("asymmetric", &mut a);
    run_case("emf_asym", &a);
    let mut e = Param::new();
    flag("each_trace", &mut e);
    flag("unweighted_fit", &mut e);
    e.set_value("add_zeros", ParamValue::Float(0.0), "", &[])
        .unwrap();
    run_case("emf_each", &e);
    let mut n = Param::new();
    flag("asymmetric", &mut n);
    flag("no_imputation", &mut n);
    n.set_value("check:width", ParamValue::Float(1.0), "", &[])
        .unwrap();
    n.set_value("check:asymmetry", ParamValue::Float(1.0), "", &[])
        .unwrap();
    run_case("emf_noimp", &n);
}

#[test]
fn elution_model_fitter_parameters() {
    let defaults = ElutionModelFitter::defaults().unwrap();
    let fitter = ElutionModelFitter::from_param(&defaults).unwrap();
    assert_eq!(fitter, ElutionModelFitter::default());
    let mut bad = Param::new();
    bad.set_value("check:boundaries", ParamValue::Float(1.5), "", &[])
        .unwrap();
    assert!(ElutionModelFitter::from_param(&bad).is_err());
}
