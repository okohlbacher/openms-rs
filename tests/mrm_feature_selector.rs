// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! MRMFeatureSelector, MRMFeatureSelectorScore, MRMFeatureSelectorQMIP and
//! MRMBatchFeatureSelector.
//!
//! Expected values are transcribed from the upstream class test at the core
//! pin `bc9cc12` (`MRMFeatureSelector_test.cpp`, evidence tier 3) with its
//! featureXML data. The class test solves with the library's `LPWrapper`
//! solver; the sections that need one are exercised with
//! [`PartitionSolver`], which is exact for the programs of the Score selector,
//! and the QMIP sections are checked structurally (tier 4) because this crate
//! ships no general solver. See `docs/MRM_FEATURE_SELECTOR_SUPPORT.md`.
#![cfg(feature = "featurexml")]

use openms::Error;
use openms::analysis::openswath::mrm_batch_feature_selector::{
    batch_mrm_features, batch_mrm_features_score,
};
use openms::analysis::openswath::mrm_feature_selector::{
    BoundType, ColumnType, LambdaScore, LinearProgram, LpSolution, LpSolver, MRMFeatureSelector,
    MRMFeatureSelectorQMIP, MRMFeatureSelectorScore, PartitionSolver, SelectorParameters,
    SolverStatus, VariableType, compute_score, construct_targ_trans_list, remove_spaces,
    weight_score,
};
use openms::format::featurexml;
use openms::kernel::{Feature, FeatureMap};
use openms::metadata::MetaValue;
use std::collections::BTreeMap;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn real_similar(a: f64, b: f64) {
    if a == b {
        return;
    }
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

const FEATURES: &str =
    "MRMFeatureSelector_150601_0_BloodProject01_PLT_QC_Broth-1_1_reduced.featureXML";
const FEATURES_SMALL: &str = "MRMFeatureSelector_100ug.featureXML";

fn float(map: &BTreeMap<String, MetaValue>, key: &str) -> f64 {
    map[key].as_f64().unwrap()
}

#[test]
fn score_select_mrm_feature() {
    // MRMFeatureSelector_test.cpp, MRMFeatureSelectorScore::selectMRMFeature().
    let feature_map = featurexml::load(data(FEATURES)).unwrap();
    assert_eq!(feature_map.len(), 210);

    let mut parameters = SelectorParameters {
        select_transition_group: true,
        segment_window_length: -1,
        segment_step_length: -1,
        variable_type: VariableType::Integer,
        optimal_threshold: 0.5,
        ..SelectorParameters::default()
    };
    parameters.score_weights = [
        ("sn_ratio".to_owned(), LambdaScore::Log),
        ("peak_apices_sum".to_owned(), LambdaScore::Log),
    ]
    .into_iter()
    .collect();

    let selector = MRMFeatureSelectorScore::new(PartitionSolver);
    let output_selected = selector
        .select_mrm_feature(&feature_map, &parameters)
        .unwrap();
    assert_eq!(output_selected.len(), 37);

    let f1 = &output_selected.features[0].subordinates[0];
    real_similar(float(&f1.metadata, "peak_apex_int"), 286.0);
    assert_eq!(
        f1.metadata["native_id"].as_str().unwrap(),
        "23dpg.23dpg_1.Heavy"
    );
    real_similar(f1.rt, 16.7592102584839);

    let f2 = &output_selected.features[9].subordinates[0];
    real_similar(float(&f2.metadata, "peak_apex_int"), 8671.5);
    assert_eq!(
        f2.metadata["native_id"].as_str().unwrap(),
        "Pool_2pg_3pg.Pool_2pg_3pg_1.Heavy"
    );
    real_similar(f2.rt, 16.1933587515513);
}

#[test]
fn remove_spaces_section() {
    // removeSpaces_()
    assert_eq!(remove_spaces("h e ll o"), "hello");
    assert_eq!(remove_spaces("hello"), "hello");
    assert_eq!(remove_spaces(""), "");
    assert_eq!(remove_spaces("A    B"), "AB");
}

#[test]
fn construct_targ_trans_list_section() {
    // constructTargTransList_()
    let feature_map = featurexml::load(data(FEATURES)).unwrap();
    let (mut time_to_name, feature_name_map) =
        construct_targ_trans_list(&feature_map, true).unwrap();
    assert_eq!(time_to_name.len(), 37);
    assert_eq!(feature_name_map.len(), 37);
    // std::sort of std::pair<double, std::string>
    time_to_name.sort_by(|a, b| a.0.total_cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
    let p = &time_to_name[0];
    real_similar(p.0, 0.0);
    assert_eq!(p.1, "arg-L");
    let p = &time_to_name[1];
    real_similar(p.0, 1.314232178);
    assert_eq!(p.1, "asn-L");
    let p = &time_to_name[5];
    real_similar(p.0, 1.421901525);
    assert_eq!(p.1, "citr-L");
    let p = &time_to_name[7];
    real_similar(p.0, 2.667749413);
    assert_eq!(p.1, "AICAr");
    let p = time_to_name.last().unwrap();
    real_similar(p.0, 99.98770892);
    assert_eq!(p.1, "accoa");
}

#[test]
fn weight_score_section() {
    // weightScore_()
    real_similar(weight_score(3413.0, LambdaScore::Linear), 3413.0);
    real_similar(weight_score(341.0, LambdaScore::Inverse), 0.002932551);
    real_similar(weight_score(341.0, LambdaScore::Log), 5.831882477);
    real_similar(weight_score(96640.0, LambdaScore::InverseLog), 0.087117);
    real_similar(weight_score(341.0, LambdaScore::InverseLog10), 0.394827074);
}

#[test]
fn compute_score_section() {
    // computeScore_()
    let mut feature = Feature::default();
    feature.metadata.insert(
        "sn_ratio".into(),
        MetaValue::try_from(6.84619503982874).unwrap(),
    );
    feature.metadata.insert(
        "peak_apices_sum".into(),
        MetaValue::try_from(96640.0).unwrap(),
    );
    let weights = |pairs: &[(&str, LambdaScore)]| -> BTreeMap<String, LambdaScore> {
        pairs.iter().map(|(k, v)| ((*k).to_owned(), *v)).collect()
    };
    real_similar(
        compute_score(&feature, &weights(&[("sn_ratio", LambdaScore::InverseLog)])).unwrap(),
        0.5198334582314795,
    );
    real_similar(
        compute_score(
            &feature,
            &weights(&[("peak_apices_sum", LambdaScore::InverseLog10)]),
        )
        .unwrap(),
        0.20059549093267626,
    );
    real_similar(
        compute_score(
            &feature,
            &weights(&[
                ("sn_ratio", LambdaScore::InverseLog),
                ("peak_apices_sum", LambdaScore::InverseLog10),
            ]),
        )
        .unwrap(),
        0.10427624775717449,
    );
    // Checks for bad input
    feature
        .metadata
        .insert("sn_ratio".into(), MetaValue::try_from(0.0).unwrap());
    feature
        .metadata
        .insert("peak_apices_sum".into(), MetaValue::try_from(0.0).unwrap());
    feature.metadata.insert(
        "var_xcorr_coelution".into(),
        MetaValue::try_from(-1.0).unwrap(),
    );
    real_similar(
        compute_score(&feature, &weights(&[("sn_ratio", LambdaScore::Inverse)])).unwrap(),
        1.0,
    );
    real_similar(
        compute_score(&feature, &weights(&[("peak_apices_sum", LambdaScore::Log)])).unwrap(),
        1.0,
    );
    real_similar(
        compute_score(
            &feature,
            &weights(&[("var_xcorr_coelution", LambdaScore::Log)]),
        )
        .unwrap(),
        1.0,
    );
}

/// A solver that refuses everything, to prove a problem is built before the
/// solver is asked.
struct Recorder(std::cell::RefCell<Option<LinearProgram>>);
impl LpSolver for Recorder {
    fn solve(&self, problem: &LinearProgram) -> openms::Result<LpSolution> {
        *self.0.borrow_mut() = Some(problem.clone());
        Ok(LpSolution {
            status: SolverStatus::Infeasible,
            column_values: vec![0.0; problem.columns.len()],
        })
    }
}

fn candidate(group: &str, id: u64, rt: f64, assay_rt: f64, sn: f64) -> Feature {
    let mut f = Feature::new(rt, 0.0, 0.0);
    f.unique_id = id;
    f.metadata.insert("PeptideRef".into(), group.into());
    f.metadata
        .insert("assay_rt".into(), MetaValue::try_from(assay_rt).unwrap());
    f.metadata
        .insert("sn_ratio".into(), MetaValue::try_from(sn).unwrap());
    let mut sub = Feature::default();
    sub.metadata
        .insert("native_id".into(), format!("{group}.t").into());
    f.subordinates.push(sub);
    f
}

#[test]
fn native_qmip_problem_has_the_source_shape() {
    // Tier 4: the column and row structure derived by reading
    // MRMFeatureSelectorQMIP::optimize. Two groups with two candidates each,
    // one neighbour each way: 4 candidate variables; per ordered pair of
    // candidates of different groups one `qp` and one `ABS` column and five
    // constraints (QP1, QP2, QP3, obj+, obj-); one fixed row per group.
    let map = FeatureMap::from_features(vec![
        candidate("A", 1, 1.0, 1.0, 5.0),
        candidate("A", 2, 1.2, 1.0, 7.0),
        candidate("B", 3, 2.0, 2.0, 5.0),
        candidate("B", 4, 2.3, 2.0, 9.0),
    ]);
    let (time_to_name, feature_name_map) = construct_targ_trans_list(&map, true).unwrap();
    let mut parameters = SelectorParameters {
        nn_threshold: 1,
        variable_type: VariableType::Integer,
        ..SelectorParameters::default()
    };
    parameters
        .score_weights
        .insert("sn_ratio".into(), LambdaScore::Linear);
    let selector = MRMFeatureSelectorQMIP::new(Recorder(Default::default()));
    let (problem, variables) = selector
        .build_problem(&time_to_name, &feature_name_map, &parameters)
        .unwrap();
    assert_eq!(variables.len(), 4);
    // 4 candidates + 8 ordered pairs * 2 columns
    assert_eq!(problem.columns.len(), 4 + 8 * 2);
    // 8 ordered pairs * 5 constraints + 2 group rows
    assert_eq!(problem.rows.len(), 8 * 5 + 2);
    let c = &problem.columns[problem.column_index("A_1").unwrap()];
    assert_eq!(c.column_type, ColumnType::Integer);
    assert_eq!(c.bound_type, BoundType::DoubleBounded);
    assert_eq!((c.lower, c.upper, c.objective), (0.0, 1.0, 0.0));
    // the pair columns are continuous; the ABS column is unbounded, cost 1
    let abs = &problem.columns[problem.column_index("A_0-B_0-ABS").unwrap()];
    assert_eq!(abs.column_type, ColumnType::Continuous);
    assert_eq!(abs.bound_type, BoundType::Unbounded);
    assert_eq!(abs.objective, 1.0);
    let qp = &problem.columns[problem.column_index("A_0-B_0").unwrap()];
    assert_eq!(qp.column_type, ColumnType::Continuous);
    // the objective rows: -abs + score * qp <= 0 with
    // score = 1 * s1 * s2 * ((rt1 - rt2) - (assay_rt1 - assay_rt2))
    let row = problem
        .rows
        .iter()
        .find(|r| r.name == "A_0-B_0-obj+")
        .unwrap();
    let expected = 5.0 * 5.0 * ((1.0 - 2.0) - (1.0 - 2.0));
    assert_eq!(row.values[0], -1.0);
    assert_eq!(row.values[1], expected);
    assert_eq!((row.lower, row.upper), (-1.0, 0.0));
    assert_eq!(row.bound_type, BoundType::UpperBoundOnly);
    let row = problem
        .rows
        .iter()
        .find(|r| r.name == "A_1-B_1-obj+")
        .unwrap();
    // candidates A_1 (rt 1.2, sn 7) and B_1 (rt 2.3, sn 9)
    let expected = 7.0 * 9.0 * ((1.2 - 2.3) - (1.0 - 2.0));
    assert!((row.values[1] - expected).abs() < 1e-12);
    // QP3 is an upper bound row over three columns
    let qp3 = problem
        .rows
        .iter()
        .find(|r| r.name == "A_0-B_0-QP3")
        .unwrap();
    assert_eq!(qp3.values, vec![1.0, 1.0, -1.0]);
    assert_eq!(qp3.bound_type, BoundType::UpperBoundOnly);
    // the group rows are fixed to one
    let fixed = problem
        .rows
        .iter()
        .find(|r| r.name == "A_constraint")
        .unwrap();
    assert_eq!(fixed.bound_type, BoundType::Fixed);
    assert_eq!((fixed.lower, fixed.upper), (1.0, 1.0));
    assert_eq!(fixed.indices.len(), 2);
}

#[test]
fn native_failed_solve_is_an_error() {
    // MRMFeatureSelector.cpp: the status check added for issue #9944; a failed
    // solve must not silently select nothing.
    let map = FeatureMap::from_features(vec![candidate("A", 1, 1.0, 1.0, 5.0)]);
    let parameters = SelectorParameters {
        segment_window_length: -1,
        segment_step_length: -1,
        ..SelectorParameters::default()
    };
    let selector = MRMFeatureSelectorScore::new(Recorder(Default::default()));
    let err = selector.select_mrm_feature(&map, &parameters).unwrap_err();
    assert!(matches!(err, Error::InvalidValue(_)));
    // the problem handed to the solver: one fixed row over the one candidate
    let problem = selector.solver.0.borrow().clone().unwrap();
    assert_eq!(problem.columns.len(), 1);
    assert_eq!(problem.rows.len(), 1);
    assert_eq!(problem.rows[0].name, "A_constraint");
}

#[test]
fn native_segments_follow_window_and_step() {
    // Tier 4: selectMRMFeature cuts the sorted list into
    // ceil(size / step) segments of at most `window` entries.
    struct Count(std::cell::RefCell<Vec<usize>>);
    impl MRMFeatureSelector for Count {
        fn optimize(
            &self,
            time_to_name: &[(f64, String)],
            _: &BTreeMap<String, Vec<Feature>>,
            _: &SelectorParameters,
        ) -> openms::Result<Vec<String>> {
            self.0.borrow_mut().push(time_to_name.len());
            Ok(Vec::new())
        }
    }
    let features: Vec<Feature> = (0..10)
        .map(|i| candidate(&format!("g{i}"), i + 1, i as f64, i as f64, 1.0))
        .collect();
    let map = FeatureMap::from_features(features);
    let mut parameters = SelectorParameters {
        segment_window_length: 4,
        segment_step_length: 3,
        ..SelectorParameters::default()
    };
    let counter = Count(Default::default());
    let selected = counter.select_mrm_feature(&map, &parameters).unwrap();
    assert!(selected.is_empty());
    // ceil(10 / 3) = 4 segments: [0,4) [3,7) [6,10) [9,10)
    assert_eq!(*counter.0.borrow(), vec![4, 4, 4, 1]);
    // a zero step is the source's integer division by zero
    parameters.segment_step_length = 0;
    assert!(counter.select_mrm_feature(&map, &parameters).is_err());
}

#[test]
fn native_construct_list_reproduces_the_inverted_contains_test() {
    // MRMFeatureSelector.cpp constructTargTransList_ tests
    // `names.contains(component_name)` where it means the negation: with
    // select_transition_group false a component name reaches `time_to_name`
    // only when it was already seen, and then once per further occurrence.
    let mut f1 = candidate("A", 1, 1.0, 1.0, 5.0);
    let mut f2 = candidate("A", 2, 1.1, 1.0, 6.0);
    f1.subordinates[0]
        .metadata
        .insert("native_id".into(), "shared".into());
    f2.subordinates[0]
        .metadata
        .insert("native_id".into(), "shared".into());
    let map = FeatureMap::from_features(vec![f1, f2]);
    let (time_to_name, feature_name_map) = construct_targ_trans_list(&map, false).unwrap();
    // "A" once; "shared" is never `contains`-true, so it is never added
    assert_eq!(time_to_name.len(), 1);
    assert_eq!(feature_name_map["shared"].len(), 2);
    assert_eq!(feature_name_map["A"].len(), 2);
}

#[test]
fn batch_applies_each_parameter_set_to_the_previous_output() {
    // MRMBatchFeatureSelector.cpp batchMRMFeatures: three cycles; the empty
    // list yields an empty map. Use the Score selector, which keeps the
    // best-scoring candidate of each group.
    let map = FeatureMap::from_features(vec![
        candidate("A", 1, 1.0, 1.0, 5.0),
        candidate("A", 2, 1.2, 1.0, 7.0),
        candidate("B", 3, 2.0, 2.0, 9.0),
        candidate("B", 4, 2.3, 2.0, 4.0),
    ]);
    let mut p = SelectorParameters {
        segment_window_length: -1,
        segment_step_length: -1,
        ..SelectorParameters::default()
    };
    p.score_weights
        .insert("sn_ratio".into(), LambdaScore::Linear);
    let once = batch_mrm_features_score(&map, &[p.clone()], PartitionSolver).unwrap();
    // the minimised score: sn 5 beats 7 in A, sn 4 beats 9 in B
    let ids: Vec<u64> = once.features.iter().map(|f| f.unique_id).collect();
    assert_eq!(ids, vec![1, 4]);
    let thrice =
        batch_mrm_features_score(&map, &[p.clone(), p.clone(), p.clone()], PartitionSolver)
            .unwrap();
    assert_eq!(thrice, once);
    let none =
        batch_mrm_features(&MRMFeatureSelectorScore::new(PartitionSolver), &map, &[]).unwrap();
    assert!(none.is_empty());
}

#[test]
fn partition_solver_refuses_other_programs() {
    // PartitionSolver is exact only for sum(x) = 1 rows over disjoint columns.
    let map = FeatureMap::from_features(vec![candidate("A", 1, 1.0, 1.0, 5.0)]);
    let (time_to_name, feature_name_map) = construct_targ_trans_list(&map, true).unwrap();
    let parameters = SelectorParameters {
        nn_threshold: 1,
        ..SelectorParameters::default()
    };
    let qmip = MRMFeatureSelectorQMIP::new(PartitionSolver);
    let two = FeatureMap::from_features(vec![
        candidate("A", 1, 1.0, 1.0, 5.0),
        candidate("B", 2, 2.0, 2.0, 5.0),
    ]);
    let (tn, fm) = construct_targ_trans_list(&two, true).unwrap();
    let err = qmip.optimize(&tn, &fm, &parameters).unwrap_err();
    assert!(matches!(err, Error::Unsupported(_)));
    // a single group has no neighbours: the QMIP program is the Score program
    let selected = qmip
        .optimize(&time_to_name, &feature_name_map, &parameters)
        .unwrap();
    assert_eq!(selected, vec!["A_1".to_owned()]);
}

#[test]
fn small_data_loads() {
    // The class-test data used by the QMIP sections: 8 features of one
    // candidate peak group each, with two subordinates. Both upstream files carry
    // `date="0000-00-00T00:00:00"`, the unset-DateTime text the source writes
    // and reads back without complaint.
    let feature_map = featurexml::load(data(FEATURES_SMALL)).unwrap();
    assert_eq!(feature_map.len(), 8);
}
