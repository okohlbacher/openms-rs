// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Executed differential for `ANALYSIS/ID/Percolator.h` (evidence tier 1).
//!
//! Every `tests/data/percolator/oracle/<case>.tsv` was written by
//! `../oracle/I2-percolator/driver.cpp` linked against the Linux x86_64
//! Release `libOpenMS.so` (`openms4-release-bc9cc12-c19e494-174b576`, GCC
//! 14.4.0, `-O3 -mssse3 -ffp-contract=off`), run twice on ibminode06 with
//! identical bytes. A case file carries the parameters, the input rows at
//! `%.17g`, and the C++ output: SVM weights (or the trained model's), pi0,
//! and per-row score, q-value and PEP. Each case is replayed through the
//! port and compared under the contract documented in
//! `docs/PERCOLATOR_SUPPORT.md` ("Checked boundaries and evidence").
//! See `tests/data/percolator_provenance.json` for hashes.

use openms::analysis::id::percolator::{
    FeatureNormalizer, MainScore, PepMethod, Percolator, PercolatorParameters, read_model,
    write_model,
};
use openms::analysis::id::percolator_types::{RescoreInput, RescoreOutput};
use std::path::{Path, PathBuf};

fn case_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/percolator/oracle")
        .join(name)
}

struct Case {
    mode: String,
    params: Vec<(String, String)>,
    input: RescoreInput,
    score_input: Option<RescoreInput>,
    weights: Vec<f64>,
    pi0: f64,
    output: RescoreOutput,
}

fn parse_f64(text: &str) -> f64 {
    match text {
        "nan" | "-nan" => f64::NAN,
        "inf" => f64::INFINITY,
        "-inf" => f64::NEG_INFINITY,
        _ => text.parse().unwrap(),
    }
}

fn read_input(lines: &mut std::iter::Peekable<std::str::Lines<'_>>, header: &str) -> RescoreInput {
    let fields: Vec<&str> = header.split('\t').collect();
    let rows: usize = fields[1].parse().unwrap();
    let names_line = lines.next().unwrap();
    let mut input = RescoreInput {
        feature_names: names_line.split('\t').skip(1).map(str::to_owned).collect(),
        ..Default::default()
    };
    let (mut scans, mut files, mut exp, mut calc, mut keys) =
        (vec![], vec![], vec![], vec![], vec![]);
    for _ in 0..rows {
        let f: Vec<&str> = lines.next().unwrap().split('\t').collect();
        assert_eq!(f[0], "row");
        input.is_decoy.push(f[1] == "1");
        if f[2] != "-" {
            scans.push(f[2].parse().unwrap());
        }
        if f[3] != "-" {
            files.push(f[3].parse().unwrap());
        }
        if f[4] != "-" {
            exp.push(parse_f64(f[4]));
        }
        if f[5] != "-" {
            calc.push(parse_f64(f[5]));
        }
        if f[6] != "-" {
            keys.push(f[6].parse().unwrap());
        }
        input
            .features
            .push(f[7..].iter().map(|v| parse_f64(v)).collect());
    }
    input.scan_numbers = scans;
    input.spec_file_numbers = files;
    input.exp_masses = exp;
    input.calc_masses = calc;
    input.cv_group_keys = keys;
    input
}

fn load_case(name: &str) -> Case {
    let text = std::fs::read_to_string(case_path(&format!("{name}.tsv"))).unwrap();
    let mut lines = text.lines().peekable();
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let mut case = Case {
        mode: head[2].to_owned(),
        params: Vec::new(),
        input: RescoreInput::default(),
        score_input: None,
        weights: Vec::new(),
        pi0: 0.0,
        output: RescoreOutput::default(),
    };
    while let Some(line) = lines.next() {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "param" => case.params.push((f[1].to_owned(), f[2].to_owned())),
            "input" => case.input = read_input(&mut lines, line),
            "score_input" => case.score_input = Some(read_input(&mut lines, line)),
            "svm_weights" | "model_weights" => {
                case.weights = f[1..].iter().map(|v| parse_f64(v)).collect();
            }
            "pi0" => case.pi0 = parse_f64(f[1]),
            "out" => {
                case.output.scores.push(parse_f64(f[2]));
                case.output.q_values.push(parse_f64(f[3]));
                case.output.peps.push(parse_f64(f[4]));
            }
            "end" => {}
            "error" => panic!("the C++ oracle failed on {name}: {line}"),
            other => panic!("unknown oracle record {other}"),
        }
    }
    case
}

fn configured(params: &[(String, String)]) -> Percolator {
    let mut p = PercolatorParameters::default();
    for (key, value) in params {
        match key.as_str() {
            "seed" => p.seed = value.parse().unwrap(),
            "num_threads" => p.num_threads = value.parse().unwrap(),
            "c_pos" => p.c_pos = value.parse().unwrap(),
            "c_neg" => p.c_neg = value.parse().unwrap(),
            "num_iterations" => p.num_iterations = value.parse().unwrap(),
            "nested_xval_bins" => p.nested_xval_bins = value.parse().unwrap(),
            "subset_max_train" => p.subset_max_train = value.parse().unwrap(),
            "test_fdr" => p.test_fdr = value.parse().unwrap(),
            "train_fdr" => p.train_fdr = value.parse().unwrap(),
            "initial_direction" => p.initial_direction = value.clone(),
            "post_processing_tdc" => p.post_processing_tdc = value == "true",
            "train_best_positive" => p.train_best_positive = value == "true",
            "use_pi0" => p.use_pi0 = value == "true",
            "pep_method" => {
                p.pep_method = match value.as_str() {
                    "isotonic" => PepMethod::Isotonic,
                    "nonparametric" => PepMethod::Nonparametric,
                    _ => PepMethod::LogisticRegression,
                }
            }
            "normalizer" => {
                p.normalizer = match value.as_str() {
                    "uni" => FeatureNormalizer::Uni,
                    "none" => FeatureNormalizer::None,
                    _ => FeatureNormalizer::Stdv,
                }
            }
            "report_as_main_score" => p.report_as_main_score = MainScore::None,
            other => panic!("unmapped oracle parameter {other}"),
        }
    }
    Percolator::with_parameters(p)
}

/// The largest ULP distance and relative difference between two slices.
fn distance(a: &[f64], b: &[f64]) -> (u64, f64) {
    assert_eq!(a.len(), b.len());
    let mut ulps = 0u64;
    let mut rel: f64 = 0.0;
    for (&x, &y) in a.iter().zip(b) {
        if x.to_bits() == y.to_bits() {
            continue;
        }
        let d = (x.to_bits() as i64)
            .wrapping_sub(y.to_bits() as i64)
            .unsigned_abs();
        ulps = ulps.max(d);
        let scale = x.abs().max(y.abs()).max(f64::MIN_POSITIVE);
        rel = rel.max((x - y).abs() / scale);
    }
    (ulps, rel)
}

struct Report {
    name: String,
    weights: (u64, f64),
    pi0: (u64, f64),
    scores: (u64, f64),
    q_values: (u64, f64),
    peps: (u64, f64),
}

fn run_case(name: &str) -> Report {
    let case = load_case(name);
    let mut p = configured(&case.params);
    let (weights, output) = if case.mode == "rescore" {
        let out = p.rescore(&case.input).unwrap();
        let w = p.svm_weights().first().cloned().unwrap_or_default();
        (w, out)
    } else {
        let model = p.train(&case.input).unwrap();
        let score_input = case.score_input.as_ref().unwrap();
        let out = p.score(score_input, &model).unwrap();
        // The model file is the source's saveModel bytes.
        let mut bytes = Vec::new();
        write_model(&model, &mut bytes).unwrap();
        let expected = std::fs::read(case_path(&format!("{name}.model"))).unwrap();
        let loaded = read_model(expected.as_slice()).unwrap();
        assert_eq!(
            loaded.feature_names, model.feature_names,
            "{name} model names"
        );
        if loaded.weights == model.weights {
            assert_eq!(bytes, expected, "{name}: model file bytes");
        }
        (model.weights, out)
    };
    Report {
        name: name.to_owned(),
        weights: distance(&weights, &case.weights),
        pi0: distance(&[p.pi0()], &[case.pi0]),
        scores: distance(&output.scores, &case.output.scores),
        q_values: distance(&output.q_values, &case.output.q_values),
        peps: distance(&output.peps, &case.output.peps),
    }
}

const CASES: [&str; 13] = [
    "rescore_default",
    "rescore_threads1",
    "rescore_isotonic",
    "rescore_uni_grid",
    "rescore_none",
    "rescore_nested",
    "rescore_multihit",
    "rescore_multihit_no_tdc",
    "rescore_subsample",
    "rescore_neg_direction",
    "rescore_no_pi0",
    "tron",
    "train_a_score_b",
];

#[test]
fn oracle_cases_report() {
    let mut lines = Vec::new();
    for name in CASES {
        let r = run_case(name);
        lines.push(format!(
            "{:<26} weights {:?} pi0 {:?} scores {:?} q {:?} pep {:?}",
            r.name, r.weights, r.pi0, r.scores, r.q_values, r.peps
        ));
    }
    for line in &lines {
        println!("{line}");
    }
}
