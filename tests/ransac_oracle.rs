// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Executed differential (tier 1) for `ML/REGRESSION`, `ML/RANSAC` and
//! `ANALYSIS/OPENSWATH/MRMRTNormalizer`.
//!
//! `tests/data/ransac_oracle_expected.tsv` is the output of
//! `../oracle/O4b-ransac/driver.cpp` linked against the Linux x86_64 Release
//! install `openms4-release-bc9cc12-c19e494-174b576` and run twice on
//! ibminode06 (byte-identical runs). The driver and this test read the same
//! datasets and cases (`ransac_oracle_datasets.tsv`, `ransac_oracle_cases.tsv`,
//! written by the oracle's `gen_data.py`). Every double is compared as its bit
//! pattern, with one documented exception: the Student's t quantile and the two
//! confidence borders computed from it, which come from `statrs` rather than
//! Boost.Math (see [`T_QUANTILE_KEYS`]).
//!
//! `removeOutliersRANSAC` seeds its RANSAC with `time(nullptr)`; the driver
//! interposes `time`, and the case's time value is passed here as the seed.

use openms::analysis::openswath::mrm_rt_normalizer::{
    OutlierMethod, chauvenet, chauvenet_probability, compute_binned_coverage,
    jackknife_outlier_candidate, remove_outliers_iterative, remove_outliers_ransac,
    residual_outlier_candidate,
};
use openms::math::ml::linear_regression::LinearRegression;
use openms::math::ml::linear_regression_without_intercept::LinearRegressionWithoutIntercept;
use openms::math::ml::quadratic_regression::QuadraticRegression;
use openms::math::ml::ransac::{
    DPair, Ransac, RansacModel, RansacModelLinear, RansacModelQuadratic,
};
use std::collections::BTreeMap;

const DATASETS: &str = include_str!("data/ransac_oracle_datasets.tsv");
const CASES: &str = include_str!("data/ransac_oracle_cases.tsv");
const EXPECTED: &str = include_str!("data/ransac_oracle_expected.tsv");

/// Records that depend on the Student's t quantile. Boost.Math computes it in
/// `long double`; the port polishes `statrs`'s AS 64/109 inversion with Newton
/// steps in `f64`. They are compared to a relative tolerance of `1e-12`
/// instead of bit for bit. Measured on this oracle: at most `5.6e-14` for
/// `t_value` and `5.7e-13` for `lower`/`upper`, which inherit the quantile's
/// error amplified by the interval's `1 / (1 - g)`.
const T_QUANTILE_KEYS: [&str; 3] = ["t_value", "lower", "upper"];

/// Driver-internal diagnostics with no Rust counterpart.
const SKIPPED_KEYS: [&str; 1] = ["time_calls"];

struct Dataset {
    x: Vec<f64>,
    y: Vec<f64>,
    w: Vec<f64>,
}

impl Dataset {
    fn pairs(&self) -> Vec<DPair> {
        self.x.iter().copied().zip(self.y.iter().copied()).collect()
    }
}

fn parse_list(field: &str) -> Vec<f64> {
    if field.is_empty() {
        return Vec::new();
    }
    field
        .split(',')
        .map(|v| v.parse::<f64>().unwrap())
        .collect()
}

fn datasets() -> BTreeMap<String, Dataset> {
    DATASETS
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            (
                fields[0].to_string(),
                Dataset {
                    x: parse_list(fields[1]),
                    y: parse_list(fields[2]),
                    w: parse_list(fields.get(3).copied().unwrap_or("")),
                },
            )
        })
        .collect()
}

#[derive(Default)]
struct Records(Vec<(String, String)>);

impl Records {
    fn put(&mut self, key: &str, value: String) {
        self.0.push((key.to_string(), value));
    }
    fn f64(&mut self, key: &str, value: f64) {
        self.put(key, format!("{:016x}", value.to_bits()));
    }
    fn pairs(&mut self, key: &str, pairs: &[DPair]) {
        self.put(&format!("{key}.count"), pairs.len().to_string());
        for (i, (x, y)) in pairs.iter().enumerate() {
            self.put(
                &format!("{key}.{i}"),
                format!("{:016x},{:016x}", x.to_bits(), y.to_bits()),
            );
        }
    }
    fn linreg(&mut self, r: &LinearRegression) {
        self.f64("intercept", r.intercept());
        self.f64("slope", r.slope());
        self.f64("x_intercept", r.x_intercept());
        self.f64("lower", r.lower());
        self.f64("upper", r.upper());
        self.f64("t_value", r.t_value());
        self.f64("r_squared", r.r_squared());
        self.f64("stand_dev_res", r.stand_dev_res());
        self.f64("mean_res", r.mean_res());
        self.f64("stand_err_slope", r.stand_err_slope());
        self.f64("chi_squared", r.chi_squared());
        self.f64("rsd", r.rsd());
    }
}

fn model_records<M: RansacModel>(
    model: &M,
    fit: &Dataset,
    test: &Dataset,
    threshold: f64,
    out: &mut Records,
) -> openms::Result<()> {
    let pairs = fit.pairs();
    let coeff = model.rm_fit(&pairs)?;
    out.put("coeff.count", coeff.len().to_string());
    for (i, c) in coeff.iter().enumerate() {
        out.f64(&format!("coeff.{i}"), *c);
    }
    out.f64("rsq", model.rm_rsq(&pairs)?);
    out.f64("rss", model.rm_rss(&pairs, &coeff)?);
    let inliers = model.rm_inliers(&test.pairs(), &coeff, threshold)?;
    out.pairs("inliers", &inliers);
    Ok(())
}

fn ransac_records<M: RansacModel + Default>(
    d: &Dataset,
    a: &[&str],
    out: &mut Records,
) -> openms::Result<()> {
    let mut ransac: Ransac<M> = Ransac::new(a[3].parse().unwrap());
    if a[4] != "-" {
        ransac.set_seed(a[4].parse().unwrap());
    }
    let result = ransac.ransac(
        &d.pairs(),
        a[5].parse().unwrap(),
        a[6].parse().unwrap(),
        a[7].parse().unwrap(),
        a[8].parse().unwrap(),
        a[9] == "1",
    )?;
    out.pairs("result", &result);
    Ok(())
}

fn run_case(data: &BTreeMap<String, Dataset>, a: &[&str], out: &mut Records) -> openms::Result<()> {
    let d = &data[a[2]];
    let num = |i: usize| a[i].parse::<f64>().unwrap();
    match a[1] {
        "linreg" | "linregw" => {
            let mut r = LinearRegression::new();
            if a[1] == "linreg" {
                r.compute_regression(num(3), &d.x, &d.y, a[4] == "1")?;
            } else {
                r.compute_regression_weighted(num(3), &d.x, &d.y, &d.w, a[4] == "1")?;
            }
            out.linreg(&r);
        }
        "linreg0" => {
            let mut r = LinearRegressionWithoutIntercept::new();
            r.add_data_slices(&d.x, &d.y)?;
            out.f64("slope", r.slope().unwrap_or(f64::NAN));
        }
        "quadreg" | "quadregw" => {
            let mut r = QuadraticRegression::new();
            if a[1] == "quadreg" {
                r.compute_regression(&d.x, &d.y)?;
            } else {
                r.compute_regression_weighted(&d.x, &d.y, &d.w)?;
            }
            out.f64("a", r.a());
            out.f64("b", r.b());
            out.f64("c", r.c());
            out.f64("chi_squared", r.chi_squared());
            out.f64("eval100", r.eval(100.0));
        }
        "rmlin" => model_records(&RansacModelLinear, d, &data[a[3]], num(4), out)?,
        "rmquad" => model_records(&RansacModelQuadratic, d, &data[a[3]], num(4), out)?,
        "ransac_lin" => ransac_records::<RansacModelLinear>(d, a, out)?,
        "ransac_quad" => ransac_records::<RansacModelQuadratic>(d, a, out)?,
        "jack" => out.put("pos", jackknife_outlier_candidate(&d.x, &d.y)?.to_string()),
        "resid" => out.put("pos", residual_outlier_candidate(&d.x, &d.y)?.to_string()),
        "rt_iter" => {
            let method: OutlierMethod = a[6].parse()?;
            let result =
                remove_outliers_iterative(&d.pairs(), num(3), num(4), a[5] == "1", method)?;
            out.pairs("result", &result);
        }
        "chauv" => {
            let pos: usize = a[3].parse().unwrap();
            out.f64("prob", chauvenet_probability(&d.y, pos)?);
            out.put("criterion", chauvenet(&d.y, pos)?.to_string());
        }
        "binned" => {
            let covered = compute_binned_coverage(
                (num(3), num(4)),
                &d.pairs(),
                a[5].parse().unwrap(),
                a[6].parse().unwrap(),
                a[7].parse().unwrap(),
            )?;
            out.put("covered", covered.to_string());
        }
        "rt_ransac" => {
            let result = remove_outliers_ransac(
                &d.pairs(),
                num(3),
                num(4),
                a[5].parse().unwrap(),
                num(6),
                a[7].parse().unwrap(),
                a[8].parse().unwrap(),
            )?;
            out.pairs("result", &result);
        }
        other => panic!("unknown op {other}"),
    }
    Ok(())
}

fn relative_difference(a: &str, b: &str) -> f64 {
    let a = f64::from_bits(u64::from_str_radix(a, 16).unwrap());
    let b = f64::from_bits(u64::from_str_radix(b, 16).unwrap());
    if a == b { 0.0 } else { ((a - b) / b).abs() }
}

#[test]
fn every_oracle_record_matches() {
    let data = datasets();
    let mut expected: BTreeMap<String, Vec<(String, String)>> = BTreeMap::new();
    for line in EXPECTED.lines().filter(|line| !line.is_empty()) {
        let mut fields = line.splitn(3, '\t');
        let case = fields.next().unwrap().to_string();
        let key = fields.next().unwrap().to_string();
        let value = fields.next().unwrap_or("").to_string();
        if SKIPPED_KEYS.contains(&key.as_str()) {
            continue;
        }
        expected.entry(case).or_default().push((key, value));
    }

    let mut cases = 0;
    let mut compared = 0;
    let mut worst_t = 0.0_f64;
    for line in CASES
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        let a: Vec<&str> = line.split('\t').collect();
        let want = expected
            .get(a[0])
            .unwrap_or_else(|| panic!("no oracle records for {}", a[0]));
        cases += 1;
        let mut got = Records::default();
        let outcome = run_case(&data, &a, &mut got);
        if want.len() == 1 && want[0].0 == "exception" {
            assert!(
                outcome.is_err(),
                "{}: the source throws `{}`, the port returned {:?}",
                a[0],
                want[0].1,
                got.0
            );
            compared += 1;
            continue;
        }
        if let Err(error) = outcome {
            panic!(
                "{}: the port failed with `{error}`, the source did not",
                a[0]
            );
        }
        assert_eq!(
            got.0.iter().map(|r| &r.0).collect::<Vec<_>>(),
            want.iter().map(|r| &r.0).collect::<Vec<_>>(),
            "{}: record keys",
            a[0]
        );
        for ((key, value), (_, expected_value)) in got.0.iter().zip(want) {
            if T_QUANTILE_KEYS.contains(&key.as_str()) {
                let difference = relative_difference(value, expected_value);
                worst_t = worst_t.max(difference);
                assert!(
                    difference <= 1e-12,
                    "{} {key}: {value} vs {expected_value} (relative {difference:e})",
                    a[0]
                );
            } else {
                assert_eq!(value, expected_value, "{} {key}", a[0]);
            }
            compared += 1;
        }
    }
    assert_eq!(cases, 128);
    assert!(compared > 1500, "only {compared} records compared");
    eprintln!(
        "{cases} cases, {compared} records; worst t-quantile relative difference {worst_t:e}"
    );
}
