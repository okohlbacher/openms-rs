// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Evidence tier 1: the lane's models, superimposers, identification
//! alignment and evaluation against the executed C++ Release build.
//!
//! `../oracle/M1-mapalign/driver.cpp` (outside this repository) links the
//! Linux x86_64 Release install `openms4-release-bc9cc12-c19e494-174b576`
//! (conda-forge gcc 14.4.0, libOpenMS's `-O3 -DNDEBUG` core flags) and ran on
//! **ibminode06** in two invocations of two fresh-process runs each, all four
//! byte-identical; its records are
//! `tests/data/map_alignment_oracle_records.tsv` and
//! `tests/data/map_alignment_oracle_summaries.txt`. The synthetic inputs
//! (`map_alignment_oracle_maps.tsv`, `map_alignment_oracle_lowess_noisy.tsv`)
//! come from the oracle's `gen_cases.py`; the others are the upstream
//! class-test files already in `tests/data/`. Hashes and the host record are
//! in `tests/data/map_alignment_provenance.json`.

#![cfg(all(feature = "idxml", feature = "consensusxml"))]

use openms::analysis::mapmatching::alignment_evaluation::{
    MapAlignmentEvaluationAlgorithm, MapAlignmentEvaluationPrecision, MapAlignmentEvaluationRecall,
};
use openms::analysis::mapmatching::alignment_identification::MapAlignmentAlgorithmIdentification;
use openms::analysis::mapmatching::pose_clustering_affine::PoseClusteringAffineSuperimposer;
use openms::analysis::mapmatching::pose_clustering_shift::PoseClusteringShiftSuperimposer;
use openms::analysis::mapmatching::superimposer::BaseSuperimposer;
use openms::analysis::transformations::*;
use openms::concept::progress_logger::ProgressReporter;
use openms::format::{consensusxml, idxml};
use openms::kernel::{BaseFeature, ConsensusFeature, ConsensusMap, Peak2D};
use openms::param::{Param, ParamValue};
use openms::processing::spline::BoundaryCondition;
use std::collections::BTreeMap;

type Records = BTreeMap<(String, String), f64>;

/// Case name, data, metric, candidate spans (empty: the default grid), folds.
type AutoSpanCase<'a> = (
    &'a str,
    &'a [DataPoint],
    CrossValidationMetric,
    &'a [f64],
    i32,
);

/// Case name, reference run and parameter overrides.
type IdentificationCase<'a> = (&'a str, Option<usize>, &'a [(&'a str, ParamValue)]);

fn records() -> Records {
    include_str!("data/map_alignment_oracle_records.tsv")
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            ((f[0].to_owned(), f[1].to_owned()), f[2].parse().unwrap())
        })
        .collect()
}

fn pairs(text: &str) -> Vec<DataPoint> {
    text.lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .map(|l| {
            let f: Vec<f64> = l.split_whitespace().map(|v| v.parse().unwrap()).collect();
            DataPoint::new(f[0], f[1])
        })
        .collect()
}

fn sine() -> Vec<DataPoint> {
    pairs(include_str!("data/transformations_lowess_input.tsv"))
}

fn grid() -> Vec<f64> {
    (0..=80).map(|k| -5.0 + 0.125 * f64::from(k)).collect()
}

/// Compares every value a case produced with the executed C++ records.
struct Check<'a> {
    records: &'a Records,
    compared: usize,
    identical: usize,
    worst: f64,
    failures: Vec<String>,
}

impl<'a> Check<'a> {
    fn new(records: &'a Records) -> Self {
        Self {
            records,
            compared: 0,
            identical: 0,
            worst: 0.0,
            failures: Vec::new(),
        }
    }

    fn expected(&self, case: &str, key: &str) -> Option<f64> {
        self.records
            .get(&(case.to_owned(), key.to_owned()))
            .copied()
    }

    fn errored(&self, case: &str) -> bool {
        self.expected(case, "error").is_some()
    }

    fn value(&mut self, case: &str, key: &str, actual: f64) {
        let Some(expected) = self.expected(case, key) else {
            self.failures.push(format!("{case} {key}: no C++ record"));
            return;
        };
        self.compared += 1;
        if actual.to_bits() == expected.to_bits() {
            self.identical += 1;
            return;
        }
        let relative = (actual - expected).abs() / expected.abs().max(1.0);
        if relative > self.worst {
            self.worst = relative;
        }
        // The records come from the Linux x86_64 Release build, and on Linux
        // x86_64 every value is required bit for bit (3780 of 3780 were, on
        // dax). Elsewhere the platform libm's `ln` and `exp` (the affine
        // superimposer's log-scaling histogram) may round differently in the
        // last bit, so those targets compare at 1e-12 relative instead.
        let linux_x86_64 = cfg!(all(target_os = "linux", target_arch = "x86_64"));
        if !linux_x86_64 && relative <= 1e-12 {
            return;
        }
        self.failures.push(format!(
            "{case} {key}: rust {actual:e} c++ {expected:e} (rel {relative:e})"
        ));
    }

    fn outcome(&mut self, case: &str, result: openms::Result<()>) {
        match (result, self.errored(case)) {
            (Ok(()), false) | (Err(_), true) => {}
            (Ok(()), true) => self
                .failures
                .push(format!("{case}: C++ threw, Rust did not")),
            (Err(e), false) => self.failures.push(format!("{case}: Rust failed: {e}")),
        }
    }

    fn finish(self, what: &str) {
        eprintln!(
            "{what}: {} of {} values bit-identical, worst relative difference {:e}",
            self.identical, self.compared, self.worst
        );
        assert!(self.compared > 0, "{what}: nothing compared");
        assert!(
            self.failures.is_empty(),
            "{what}: {} differences:\n{}",
            self.failures.len(),
            self.failures.join("\n")
        );
    }
}

#[test]
fn b_spline_models_match_the_release_build() {
    let records = records();
    let mut check = Check::new(&records);
    let configs = [
        ("w0n5b2", 0.0, 5, BoundaryCondition::ZeroSecond),
        ("w1n0b2", 1.0, 0, BoundaryCondition::ZeroSecond),
        ("w0n8b0", 0.0, 8, BoundaryCondition::ZeroEndpoints),
        ("w0n8b1", 0.0, 8, BoundaryCondition::ZeroFirst),
    ];
    for (name, wavelength, num_nodes, boundary_condition) in configs {
        for extrapolation in [
            BSplineExtrapolation::Linear,
            BSplineExtrapolation::BSpline,
            BSplineExtrapolation::Constant,
            BSplineExtrapolation::GlobalLinear,
        ] {
            let case = format!("bspline:{name}:{}", extrapolation.name());
            let result = BSplineModel::fit(
                &sine(),
                BSplineOptions {
                    wavelength,
                    num_nodes,
                    extrapolation,
                    boundary_condition,
                    ..Default::default()
                },
            )
            .and_then(|model| {
                for (k, x) in grid().into_iter().enumerate() {
                    let y = model.apply(x)?;
                    check.value(&case, &k.to_string(), y);
                }
                Ok(())
            });
            check.outcome(&case, result);
        }
    }
    check.finish("b-spline");
}

#[test]
fn interpolations_match_the_release_build() {
    let records = records();
    let mut check = Check::new(&records);
    let base = pairs(include_str!("data/transformations_interpolation_input.tsv"));
    let golden = pairs(include_str!(
        "data/transformations_interpolation_golden.tsv"
    ));
    for interpolation in [
        Interpolation::Linear,
        Interpolation::CubicSpline,
        Interpolation::Akima,
    ] {
        for extrapolation in [
            Extrapolation::TwoPointLinear,
            Extrapolation::FourPointLinear,
            Extrapolation::GlobalLinear,
        ] {
            let case = format!("interp:{}:{}", interpolation.name(), extrapolation.name());
            let result = InterpolatedModel::fit(
                &base,
                InterpolationOptions {
                    interpolation,
                    extrapolation,
                    ..Default::default()
                },
            )
            .and_then(|model| {
                for (i, point) in golden.iter().enumerate() {
                    check.value(&case, &format!("g{i}"), model.apply(point.x)?);
                }
                check.value(&case, "below", model.apply(-3.0)?);
                let last = golden[golden.len() - 1].x;
                check.value(&case, "above", model.apply(last + 3.0)?);
                Ok(())
            });
            check.outcome(&case, result);
        }
    }
    check.finish("interpolation");
}

#[test]
fn lowess_models_and_automatic_spans_match_the_release_build() {
    let records = records();
    let mut check = Check::new(&records);
    for interpolation in [
        Interpolation::Linear,
        Interpolation::CubicSpline,
        Interpolation::Akima,
    ] {
        for extrapolation in [
            Extrapolation::TwoPointLinear,
            Extrapolation::FourPointLinear,
            Extrapolation::GlobalLinear,
        ] {
            let case = format!(
                "lowess:sine:{}:{}",
                interpolation.name(),
                extrapolation.name()
            );
            let result = LowessModel::fit(
                &sine(),
                LowessOptions {
                    span: 0.3,
                    interpolation,
                    extrapolation,
                    ..Default::default()
                },
            )
            .and_then(|model| {
                for (k, x) in grid().into_iter().enumerate() {
                    check.value(&case, &k.to_string(), model.apply(x)?);
                }
                Ok(())
            });
            check.outcome(&case, result);
        }
    }
    let noisy = pairs(include_str!("data/map_alignment_oracle_lowess_noisy.tsv"));
    let autos: [AutoSpanCase<'_>; 6] = [
        ("sine:mae", &sine(), CrossValidationMetric::Mae, &[], 5),
        ("sine:rmse", &sine(), CrossValidationMetric::Rmse, &[], 5),
        ("sine:p90", &sine(), CrossValidationMetric::P90, &[], 5),
        ("noisy:mae", &noisy, CrossValidationMetric::Mae, &[], 5),
        ("noisy:p95", &noisy, CrossValidationMetric::P95, &[], 3),
        (
            "noisy:grid",
            &noisy,
            CrossValidationMetric::Rmse,
            &[0.2, 0.35, 0.5, 0.9],
            4,
        ),
    ];
    for (name, data, metric, candidates, k_folds) in autos {
        let case = format!("lowess:auto:{name}");
        let auto = LowessAutoSpan {
            metric,
            grid: candidates.to_vec(),
            k_folds,
            ..Default::default()
        };
        let result =
            LowessModel::fit_auto(data, LowessOptions::default(), &auto).and_then(|(model, _)| {
                check.value(&case, "span", model.options().span);
                for (k, x) in grid().into_iter().enumerate() {
                    check.value(&case, &k.to_string(), model.apply(x)?);
                }
                for k in 0..=60 {
                    check.value(&case, &format!("x{k}"), model.apply(f64::from(k))?);
                }
                Ok(())
            });
        check.outcome(&case, result);
    }
    check.finish("lowess");
}

#[test]
fn descriptions_print_the_release_build_summaries() {
    let records = records();
    let mut check = Check::new(&records);
    let text = include_str!("data/map_alignment_oracle_summaries.txt");
    let mut blocks = BTreeMap::new();
    for block in text.split("== ").filter(|b| !b.is_empty()) {
        let (name, body) = block.split_once('\n').unwrap();
        blocks.insert(name.to_owned(), body.to_owned());
    }
    for model in ["linear", "b_spline", "interpolated", "lowess"] {
        let case = format!("summary:{model}");
        let result = (|| {
            let mut description = TransformationDescription::new(sine())?;
            description.fit_model_with_param(model, &Param::new())?;
            let mut out = Vec::new();
            description.print_summary(&mut out)?;
            assert_eq!(
                String::from_utf8(out).unwrap(),
                blocks[&case],
                "{case} summary text"
            );
            check.value(&case, "apply0", description.apply(0.0)?);
            description.invert()?;
            check.value(&case, "inverse0", description.apply(0.0)?);
            Ok(())
        })();
        check.outcome(&case, result);
    }
    check.finish("summaries");
}

fn maps() -> BTreeMap<String, Vec<Peak2D>> {
    let mut maps: BTreeMap<String, Vec<Peak2D>> = BTreeMap::new();
    for line in include_str!("data/map_alignment_oracle_maps.tsv").lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let intensity: f64 = f[3].parse().unwrap();
        maps.entry(f[0].to_owned()).or_default().push(Peak2D::new(
            f[1].parse().unwrap(),
            f[2].parse().unwrap(),
            intensity as f32,
        ));
    }
    maps
}

fn consensus(points: &[Peak2D]) -> ConsensusMap {
    ConsensusMap::from_features(
        points
            .iter()
            .map(|p| {
                ConsensusFeature::from(BaseFeature::new(p.position[0], p.position[1], p.intensity))
            })
            .collect(),
    )
}

fn linear(check: &mut Check<'_>, case: &str, transformation: &TransformationDescription) {
    let params = transformation.model_parameters();
    check.value(
        case,
        "slope",
        params.value("slope").unwrap().to_f64().unwrap(),
    );
    check.value(
        case,
        "intercept",
        params.value("intercept").unwrap().to_f64().unwrap(),
    );
}

#[test]
fn superimposers_match_the_release_build() {
    let records = records();
    let mut check = Check::new(&records);
    let maps = maps();
    let shifts: [(&str, &str, &str, f64, f64, i64); 5] = [
        (
            "shift:default",
            "shift_model",
            "shift_scene",
            0.5,
            3.0,
            2000,
        ),
        ("shift:narrow", "shift_model", "shift_scene", 0.3, 1.0, 2000),
        ("shift:used100", "shift_model", "shift_scene", 0.5, 3.0, 100),
        (
            "shift:swapped",
            "shift_scene",
            "shift_model",
            0.5,
            3.0,
            2000,
        ),
        ("shift:dense", "dense_model", "dense_scene", 0.5, 3.0, 2000),
    ];
    for (case, model, scene, mz_pair_max_distance, shift_bucket_size, num_used_points) in shifts {
        let superimposer = PoseClusteringShiftSuperimposer {
            mz_pair_max_distance,
            shift_bucket_size,
            num_used_points,
            ..Default::default()
        };
        let result = superimposer
            .run(&consensus(&maps[model]), &consensus(&maps[scene]))
            .map(|t| linear(&mut check, case, &t));
        check.outcome(case, result);
    }
    let affines: [(&str, &str, &str, f64, f64, i64); 7] = [
        (
            "affine:default",
            "affine_model",
            "affine_scene",
            0.005,
            3.0,
            2000,
        ),
        (
            "affine:used100",
            "affine_model",
            "affine_scene",
            0.005,
            3.0,
            100,
        ),
        (
            "affine:used37",
            "affine_model",
            "affine_scene",
            0.005,
            3.0,
            37,
        ),
        (
            "affine:fine",
            "affine_model",
            "affine_scene",
            0.001,
            1.0,
            2000,
        ),
        (
            "affine:swapped",
            "affine_scene",
            "affine_model",
            0.005,
            3.0,
            2000,
        ),
        (
            "affine:shiftmaps",
            "shift_model",
            "shift_scene",
            0.005,
            3.0,
            150,
        ),
        (
            "affine:dense",
            "dense_model",
            "dense_scene",
            0.005,
            3.0,
            2000,
        ),
    ];
    for (name, model, scene, scaling_bucket_size, shift_bucket_size, num_used_points) in affines {
        let superimposer = PoseClusteringAffineSuperimposer {
            scaling_bucket_size,
            shift_bucket_size,
            num_used_points,
            ..Default::default()
        };
        let case = format!("{name}:points");
        let result = superimposer
            .superimpose_points(&maps[model], &maps[scene], &mut ProgressReporter::silent())
            .map(|s| linear(&mut check, &case, &s.transformation));
        check.outcome(&case, result);
        let case = format!("{name}:consensus");
        let result = superimposer
            .run(&consensus(&maps[model]), &consensus(&maps[scene]))
            .map(|t| linear(&mut check, &case, &t));
        check.outcome(&case, result);
    }
    check.finish("superimposers");
}

#[test]
fn identification_alignment_and_evaluation_match_the_release_build() {
    let records = records();
    let mut check = Check::new(&records);
    let peptides: Vec<_> = [
        include_bytes!("data/map_alignment_identification_1.idXML").as_slice(),
        include_bytes!("data/map_alignment_identification_2.idXML").as_slice(),
    ]
    .into_iter()
    .map(|bytes| idxml::read(bytes).unwrap().peptide_identifications)
    .collect();
    let cases: [IdentificationCase<'_>; 6] = [
        ("ident:none", None, &[]),
        ("ident:ref0", Some(0), &[]),
        ("ident:ref1", Some(1), &[]),
        (
            "ident:cutoff",
            None,
            &[
                ("score_cutoff", ParamValue::String(String::new())),
                ("min_score", ParamValue::Float(30.0)),
            ],
        ),
        (
            "ident:shift0",
            None,
            &[("max_rt_shift", ParamValue::Float(0.0))],
        ),
        (
            "ident:shift2",
            None,
            &[("max_rt_shift", ParamValue::Float(2.0))],
        ),
    ];
    for (case, reference, overrides) in cases {
        let mut params = MapAlignmentAlgorithmIdentification::defaults().unwrap();
        for (key, value) in overrides {
            let value = if *key == "score_cutoff" {
                ParamValue::String("true".into())
            } else {
                value.clone()
            };
            params.set_value(key, value, "", &[]).unwrap();
        }
        let mut aligner = MapAlignmentAlgorithmIdentification::from_param(&params).unwrap();
        let result = aligner.align(&peptides, reference).map(|alignment| {
            for (t, transformation) in alignment.transformations.iter().enumerate() {
                let data = transformation.data_points();
                check.value(case, &format!("t{t}:n"), data.len() as f64);
                check.value(
                    case,
                    &format!("t{t}:identity"),
                    f64::from(u8::from(transformation.model_type() == "identity")),
                );
                for (i, point) in data.iter().enumerate() {
                    check.value(case, &format!("t{t}:{i}:x"), point.x);
                    check.value(case, &format!("t{t}:{i}:y"), point.y);
                }
            }
        });
        check.outcome(case, result);
    }

    let input = consensusxml::read(
        include_bytes!("data/map_alignment_evaluation_in.consensusXML").as_slice(),
    )
    .unwrap();
    let truth = consensusxml::read(
        include_bytes!("data/map_alignment_evaluation_gt.consensusXML").as_slice(),
    )
    .unwrap();
    for (case, rt, mz, intensity, charge) in [
        ("eval:class", 0.1, 0.1, 100.0f32, true),
        ("eval:tight", 0.01, 0.01, 0.5, false),
        ("eval:loose", 30.0, 3.0, 10.0, false),
    ] {
        let precision = MapAlignmentEvaluationPrecision
            .evaluate(&input, &truth, rt, mz, intensity, charge)
            .unwrap();
        let recall = MapAlignmentEvaluationRecall
            .evaluate(&input, &truth, rt, mz, intensity, charge)
            .unwrap();
        check.value(case, "precision", precision);
        check.value(case, "recall", recall);
    }
    check.finish("identification and evaluation");
}
