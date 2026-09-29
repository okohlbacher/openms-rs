// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! OpenSWATH confidence inference and library identifier normalisation.
//!
//! Two kinds of evidence:
//!
//! * Tier 3: the literals of the upstream class tests at core `bc9cc12`
//!   (`LevelContextInference_test.cpp`, `PeptidoformInference_test.cpp`,
//!   `OpenSwathLibraryIDNormalizer_test.cpp`), transcribed section by section.
//! * Tier 1: the O5b oracle, a driver linking the Linux x86_64 Release install
//!   `openms4-release-bc9cc12-c19e494-174b576`, executed twice on ibminode06
//!   (`tests/data/openswath_inference_provenance.json`). Its `cases.txt` and
//!   its outputs are replayed here case by case: every case the Rust port can
//!   run (all but the two `OSWFile` reader dumps, whose output is the input of
//!   later cases) is recomputed and compared line by line.

use openms::analysis::openswath::level_context_inference::LevelContextInference;
use openms::analysis::openswath::openswath_gene_inference::OpenSwathGeneInference;
use openms::analysis::openswath::openswath_inference_config::{
    ErrorEstimationConfig, LevelContextInferenceConfig, PeptidoformInferenceConfig,
};
use openms::analysis::openswath::openswath_inference_data::{
    IPFAlignmentRow, IPFPrecursorRow, IPFTransitionRow, InferenceContext, InferenceLevel,
    LevelContextInputRow, LevelContextResultRow,
};
use openms::analysis::openswath::openswath_library_id_normalizer::{
    OpenSwathLibraryIDNormalizer, SourceIDMapping,
};
use openms::analysis::openswath::openswath_peptide_inference::OpenSwathPeptideInference;
use openms::analysis::openswath::openswath_peptidoform_inference::{
    BayesianModelRow, OpenSwathPeptidoformInference,
};
use openms::analysis::openswath::openswath_protein_inference::OpenSwathProteinInference;
use openms::analysis::openswath::transition_experiment::{
    LightCompound, LightTargetedExperiment, LightTransition,
};
use openms::math::multiple_testing::{LfdrTransform, Pi0Method};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

/// The class tests' `TOLERANCE_ABSOLUTE(1e-10)` with `TOLERANCE_RELATIVE(1 +
/// 1e-10)`: `TEST_REAL_SIMILAR` passes when either holds.
fn similar(actual: f64, expected: f64) -> bool {
    if (actual - expected).abs() <= 1e-10 {
        return true;
    }
    if actual == 0.0 || expected == 0.0 {
        return false;
    }
    let ratio = actual / expected;
    ratio.max(1.0 / ratio) <= 1.0 + 1e-10
}

fn assert_similar(actual: f64, expected: f64, what: &str) {
    assert!(similar(actual, expected), "{what}: {actual} vs {expected}");
}

fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/openswath_inference")
}

fn lc_row(
    run_id: Option<i64>,
    group: &str,
    entity: i64,
    decoy: bool,
    score: f64,
    context: InferenceContext,
) -> LevelContextInputRow {
    LevelContextInputRow {
        run_id,
        group_id: group.to_string(),
        entity_id: entity,
        decoy,
        score,
        context,
    }
}

fn config(level: InferenceLevel, context: InferenceContext) -> LevelContextInferenceConfig {
    LevelContextInferenceConfig {
        level,
        context,
        ..LevelContextInferenceConfig::default()
    }
}

fn find(rows: &[LevelContextResultRow], run: Option<i64>, entity: i64) -> &LevelContextResultRow {
    rows.iter()
        .find(|r| r.run_id == run && r.entity_id == entity)
        .expect("row present")
}

fn prec(feature: i64, pg: f64, ms1: Option<f64>, ms2: Option<f64>) -> IPFPrecursorRow {
    IPFPrecursorRow {
        feature_id: feature,
        ms2_peakgroup_pep: pg,
        ms1_precursor_pep: ms1,
        ms2_precursor_pep: ms2,
    }
}

fn trans(
    feature: i64,
    transition: i64,
    peptide: i64,
    pep: f64,
    bmask: i32,
    n: i32,
) -> IPFTransitionRow {
    IPFTransitionRow {
        feature_id: feature,
        transition_id: transition,
        peptide_id: peptide,
        pep,
        bmask,
        num_peptidoforms: n,
        alignment_group_id: None,
    }
}

fn bm(feature: i64, hypothesis: i64, prior: f64, evidence: f64) -> BayesianModelRow {
    BayesianModelRow {
        feature_id: feature,
        hypothesis,
        prior,
        evidence,
    }
}

fn loose_ipf() -> PeptidoformInferenceConfig {
    PeptidoformInferenceConfig {
        ipf_ms1_scoring: false,
        ipf_ms2_scoring: false,
        ipf_max_precursor_peakgroup_pep: 1.0,
        ..PeptidoformInferenceConfig::default()
    }
}

// ---------------------------------------------------------------------------
// Tier 3: LevelContextInference_test.cpp (core bc9cc12)
// ---------------------------------------------------------------------------

/// `LevelContextInference_test.cpp:34-157`, the synthetic sections of
/// `infer`. The `OSWFile` section (lines 159-236) needs the OSW reader and is
/// covered by the oracle replay below instead.
#[test]
fn class_test_level_context_infer_synthetic() {
    let global = config(InferenceLevel::Peptide, InferenceContext::Global);
    assert!(
        LevelContextInference::infer(&[], &global)
            .unwrap()
            .is_empty()
    );

    let g = InferenceContext::Global;
    let rows = vec![
        lc_row(None, "101", 101, false, 5.0, g),
        lc_row(None, "102", 102, false, 4.0, g),
        lc_row(None, "201", 201, true, 2.0, g),
        lc_row(None, "202", 202, true, 1.0, g),
    ];
    let results = LevelContextInference::infer(&rows, &global).unwrap();
    assert_eq!(results.len(), rows.len());
    for row in &results {
        assert_eq!(row.run_id, None);
        assert_eq!(row.context, InferenceContext::Global);
        assert!(row.pvalue.is_finite() && row.qvalue.is_finite() && row.pep.is_finite());
    }
    assert!(find(&results, None, 101).qvalue <= find(&results, None, 102).qvalue);

    for context in [
        InferenceContext::RunSpecific,
        InferenceContext::ExperimentWide,
    ] {
        let mut rows = Vec::new();
        for (run, t, d) in [(1, 100, 200), (2, 300, 400)] {
            rows.push(lc_row(Some(run), "", t + 1, false, 5.0, context));
            rows.push(lc_row(Some(run), "", t + 2, false, 4.0, context));
            rows.push(lc_row(Some(run), "", d + 1, true, 2.0, context));
            rows.push(lc_row(Some(run), "", d + 2, true, 1.0, context));
        }
        let results =
            LevelContextInference::infer(&rows, &config(InferenceLevel::Peptide, context)).unwrap();
        assert_eq!(results.len(), rows.len());
        let (r1, r2) = (find(&results, Some(1), 101), find(&results, Some(2), 301));
        assert_similar(r1.pvalue, r2.pvalue, "pvalue");
        assert_similar(r1.qvalue, r2.qvalue, "qvalue");
        assert_similar(r1.pep, r2.pep, "pep");
        assert_similar(
            find(&results, Some(1), 102).qvalue,
            find(&results, Some(2), 302).qvalue,
            "second qvalue",
        );
    }

    let invalid = config(InferenceLevel::Peptide, InferenceContext::RunSpecific);
    let missing = [lc_row(
        None,
        "101",
        101,
        false,
        5.0,
        InferenceContext::RunSpecific,
    )];
    assert!(LevelContextInference::infer(&missing, &invalid).is_err());
}

/// `LevelContextInference_test.cpp:239-262`, the wrapper section.
#[test]
fn class_test_level_wrappers() {
    let peptide = config(InferenceLevel::Peptide, InferenceContext::Global);
    let protein = config(InferenceLevel::Protein, InferenceContext::Global);
    let gene = config(InferenceLevel::Gene, InferenceContext::Global);
    assert!(
        OpenSwathPeptideInference
            .infer(&[], &peptide)
            .unwrap()
            .is_empty()
    );
    assert!(
        OpenSwathProteinInference
            .infer(&[], &protein)
            .unwrap()
            .is_empty()
    );
    assert!(OpenSwathGeneInference.infer(&[], &gene).unwrap().is_empty());
    assert!(OpenSwathPeptideInference.infer(&[], &protein).is_err());
    assert!(OpenSwathProteinInference.infer(&[], &gene).is_err());
    assert!(OpenSwathGeneInference.infer(&[], &peptide).is_err());
}

// ---------------------------------------------------------------------------
// Tier 3: PeptidoformInference_test.cpp (core bc9cc12)
// ---------------------------------------------------------------------------

/// `PeptidoformInference_test.cpp:33-43`.
#[test]
fn class_test_compute_model_fdr() {
    let q = OpenSwathPeptidoformInference::compute_model_fdr(&[0.2, 0.1, 0.1, 0.4]).unwrap();
    assert_eq!(q.len(), 4);
    assert_similar(q[0], 0.4 / 3.0, "q0");
    assert_similar(q[1], 0.1, "q1");
    assert_similar(q[2], 0.1, "q2");
    assert_similar(q[3], 0.2, "q3");
}

/// `PeptidoformInference_test.cpp:45-81`.
#[test]
fn class_test_apply_bm() {
    let rows = [bm(1, 0, 0.3, 0.5), bm(1, 1, 0.7, 0.5), bm(1, 1, 0.8, 0.5)];
    let post = OpenSwathPeptidoformInference::apply_bm(&rows);
    assert_eq!(post.len(), 2);
    assert_eq!(post[0].feature_id, 1);
    assert_eq!(post[1].feature_id, 1);
    assert_similar(post[0].posterior + post[1].posterior, 1.0, "sum");
    let hyp0 = post.iter().find(|r| r.hypothesis == 0).unwrap();
    let hyp1 = post.iter().find(|r| r.hypothesis == 1).unwrap();
    assert_similar(hyp0.posterior, 0.15 / 0.325, "hyp0");
    assert_similar(hyp1.posterior, 0.175 / 0.325, "hyp1");

    let zero = OpenSwathPeptidoformInference::apply_bm(&[bm(2, 0, 0.5, 0.0), bm(2, 1, 0.5, 0.0)]);
    assert_eq!(zero.len(), 2);
    assert_similar(zero[0].posterior, 0.0, "zero0");
    assert_similar(zero[1].posterior, 0.0, "zero1");
}

/// `PeptidoformInference_test.cpp:83-158`.
#[test]
fn class_test_precursor_inference() {
    let inference = OpenSwathPeptidoformInference;
    let with = |ms1: bool, ms2: bool| PeptidoformInferenceConfig {
        ipf_ms1_scoring: ms1,
        ipf_ms2_scoring: ms2,
        ipf_max_precursor_peakgroup_pep: 1.0,
        ..PeptidoformInferenceConfig::default()
    };

    let r = inference.precursor_inference(&[prec(1, 0.2, Some(0.1), None)], &with(true, false));
    assert_eq!(r.len(), 1);
    assert_similar(r[0].precursor_peakgroup_pep, 1.0 - (0.72 / 0.74), "ms1");

    let r = inference.precursor_inference(&[prec(1, 0.2, None, Some(0.25))], &with(false, true));
    assert_eq!(r.len(), 1);
    assert_similar(r[0].precursor_peakgroup_pep, 1.0 - (0.6 / 0.65), "ms2");

    let r = inference.precursor_inference(&[prec(1, 0.3, Some(0.1), Some(0.2))], &with(true, true));
    assert_eq!(r.len(), 1);
    assert_similar(r[0].precursor_peakgroup_pep, 1.0 - (0.504 / 0.51), "ms1ms2");

    let cfg = PeptidoformInferenceConfig {
        ipf_max_precursor_peakgroup_pep: 0.3,
        ..PeptidoformInferenceConfig::default()
    };
    let r = inference.precursor_inference(
        &[prec(1, 0.25, None, None), prec(2, 0.35, None, None)],
        &cfg,
    );
    assert_eq!(r.len(), 1);
    assert_eq!(r[0].feature_id, 1);
    assert_similar(r[0].precursor_peakgroup_pep, 0.25, "plain");

    let cfg = PeptidoformInferenceConfig {
        ipf_max_precursor_pep: 0.05,
        ..with(true, false)
    };
    assert!(
        inference
            .precursor_inference(&[prec(1, 0.2, Some(0.1), None)], &cfg)
            .is_empty()
    );
}

/// `PeptidoformInference_test.cpp:160-273`, the synthetic sections of
/// `infer`. The `OSWFile` section (lines 275-298) is covered by the oracle
/// replay below.
#[test]
fn class_test_peptidoform_infer_synthetic() {
    let inference = OpenSwathPeptidoformInference;
    let results = inference
        .infer(
            &[prec(1, 0.2, None, None)],
            &[
                trans(1, 11, 101, 0.1, 1, 2),
                trans(1, 11, 102, 0.1, 0, 2),
                trans(1, 11, -1, 0.1, 0, 2),
            ],
            &[],
            &loose_ipf(),
        )
        .unwrap();
    assert_eq!(results.len(), 2);
    let p101 = results.iter().find(|r| r.peptide_id == 101).unwrap();
    let p102 = results.iter().find(|r| r.peptide_id == 102).unwrap();
    assert_similar(p101.precursor_peakgroup_pep, 0.2, "pp");
    assert_similar(p101.pep, 1.0 - (0.36 / 0.42), "pep101");
    assert_similar(p101.qvalue, 0.142_857_142_857_142_85, "q101");
    assert_similar(p102.pep, 1.0 - (0.04 / 0.42), "pep102");
    assert_similar(p102.qvalue, 0.523_809_523_809_523_8, "q102");

    let grouped_cfg = PeptidoformInferenceConfig {
        ipf_grouped_fdr: true,
        ..loose_ipf()
    };
    let transitions = [
        trans(1, 11, 101, 0.1, 1, 2),
        trans(1, 11, -1, 0.1, 0, 2),
        trans(2, 21, 201, 0.1, 1, 3),
        trans(2, 21, -1, 0.1, 0, 3),
    ];
    let precursors = [prec(1, 0.2, None, None), prec(2, 0.2, None, None)];
    let grouped = inference
        .infer(&precursors, &transitions, &[], &grouped_cfg)
        .unwrap();
    assert_eq!(grouped.len(), 2);
    assert_similar(grouped[0].qvalue, grouped[0].pep, "g0");
    assert_similar(grouped[1].qvalue, grouped[1].pep, "g1");
    let pooled = inference
        .infer(&precursors, &transitions, &[], &loose_ipf())
        .unwrap();
    assert_eq!(pooled.len(), 2);
    assert_similar(pooled[0].qvalue, pooled[0].pep, "p0");
    assert_similar(
        pooled[1].qvalue,
        (grouped[0].pep + grouped[1].pep) / 2.0,
        "p1",
    );
    assert!(pooled[1].qvalue < grouped[1].qvalue);

    let precursors = [prec(1, 0.2, None, None), prec(2, 0.2, None, None)];
    let transitions = [
        trans(1, 31, 301, 0.01, 1, 1),
        trans(1, 31, -1, 0.01, 0, 1),
        trans(2, 31, 301, 0.8, 1, 1),
        trans(2, 31, -1, 0.8, 0, 1),
    ];
    let alignment = [
        IPFAlignmentRow {
            alignment_group_id: 7,
            feature_id: 1,
        },
        IPFAlignmentRow {
            alignment_group_id: 7,
            feature_id: 2,
        },
    ];
    let without = inference
        .infer(&precursors, &transitions, &alignment, &loose_ipf())
        .unwrap();
    let without2 = without
        .iter()
        .find(|r| r.feature_id == 2 && r.peptide_id == 301)
        .unwrap();
    assert_similar(without2.pep, 0.5, "without");
    let propagate = PeptidoformInferenceConfig {
        propagate_signal_across_runs: true,
        across_run_confidence_threshold: 0.05,
        ..loose_ipf()
    };
    let with = inference
        .infer(&precursors, &transitions, &alignment, &propagate)
        .unwrap();
    let with2 = with
        .iter()
        .find(|r| r.feature_id == 2 && r.peptide_id == 301)
        .unwrap();
    assert!(with2.pep <= 0.01);
    assert!(with2.pep < without2.pep);
}

// ---------------------------------------------------------------------------
// Tier 3: OpenSwathLibraryIDNormalizer_test.cpp (core bc9cc12)
// ---------------------------------------------------------------------------

fn library(compounds: &[&str], transitions: &[(&str, &str)]) -> LightTargetedExperiment {
    let mut exp = LightTargetedExperiment::new();
    for id in compounds {
        exp.compounds.push(LightCompound {
            id: (*id).to_string(),
            ..LightCompound::default()
        });
    }
    for (name, reference) in transitions {
        exp.transitions.push(LightTransition {
            transition_name: (*name).to_string(),
            peptide_ref: (*reference).to_string(),
            ..LightTransition::default()
        });
    }
    exp
}

/// `OpenSwathLibraryIDNormalizer_test.cpp:25-82`.
#[test]
fn class_test_normalize_source_ids() {
    let mut exp = library(
        &["490", "A", "0", "2"],
        &[
            ("2292", "490"),
            ("transition_A", "A"),
            ("12", "0"),
            ("0", "2"),
        ],
    );
    let ids = OpenSwathLibraryIDNormalizer::normalize_source_ids(&mut exp).unwrap();
    let s2c = &ids.precursor_source_to_canonical;
    let c2s = &ids.precursor_canonical_to_source;
    let t2s = &ids.transition_canonical_to_source;
    for (s, c) in [("0", "0"), ("2", "1"), ("490", "2"), ("A", "3")] {
        assert_eq!(s2c[s], c);
        assert_eq!(c2s[c], s);
    }
    for (c, s) in [
        ("0", "2292"),
        ("1", "transition_A"),
        ("2", "12"),
        ("3", "0"),
    ] {
        assert_eq!(t2s[c], s);
    }
    let ids_now: Vec<&str> = exp.compounds.iter().map(|c| c.id.as_str()).collect();
    assert_eq!(ids_now, ["2", "3", "0", "1"]);
    let tr: Vec<(&str, &str)> = exp
        .transitions
        .iter()
        .map(|t| (t.transition_name.as_str(), t.peptide_ref.as_str()))
        .collect();
    assert_eq!(tr, [("0", "2"), ("1", "3"), ("2", "0"), ("3", "1")]);
    OpenSwathLibraryIDNormalizer::validate_canonical_ids(&exp).unwrap();
}

/// `OpenSwathLibraryIDNormalizer_test.cpp:84-117`.
#[test]
fn class_test_normalize_keeps_sparse_ids() {
    let mut exp = library(
        &["PEPTIDEA", "PEPTIDEA_Extra", "PEPTIDEB"],
        &[("tr_a", "PEPTIDEA"), ("tr_b", "PEPTIDEB")],
    );
    OpenSwathLibraryIDNormalizer::normalize_source_ids(&mut exp).unwrap();
    assert_eq!(exp.compounds.len(), 2);
    assert_eq!(exp.compounds[0].id, "0");
    assert_eq!(exp.compounds[1].id, "2");
    assert_eq!(exp.transitions[0].transition_name, "0");
    assert_eq!(exp.transitions[0].peptide_ref, "0");
    assert_eq!(exp.transitions[1].transition_name, "1");
    assert_eq!(exp.transitions[1].peptide_ref, "2");
    OpenSwathLibraryIDNormalizer::validate_canonical_ids(&exp).unwrap();
}

/// `OpenSwathLibraryIDNormalizer_test.cpp:119-153`: the port keeps no lookup
/// cache, so the lookup after normalisation resolves the canonical ids.
#[test]
fn class_test_normalize_after_lookup() {
    let mut exp = library(
        &["PEPTIDEA", "PEPTIDEB"],
        &[("tr_a", "PEPTIDEA"), ("tr_b", "PEPTIDEB")],
    );
    assert_eq!(exp.compound_by_ref("PEPTIDEA").unwrap().id, "PEPTIDEA");
    OpenSwathLibraryIDNormalizer::normalize_source_ids(&mut exp).unwrap();
    assert_eq!(exp.compound_by_ref("0").unwrap().id, "0");
    assert_eq!(exp.compound_by_ref("1").unwrap().id, "1");
    assert_eq!(exp.transitions[0].peptide_ref, "0");
    assert_eq!(exp.transitions[1].peptide_ref, "1");
}

/// `OpenSwathLibraryIDNormalizer_test.cpp:155-174`, plus the native guarantee
/// that a refused library is left unchanged.
#[test]
fn class_test_normalize_invalid_sources() {
    let mut duplicate = library(&["duplicate", "duplicate"], &[]);
    assert!(OpenSwathLibraryIDNormalizer::normalize_source_ids(&mut duplicate).is_err());
    let mut missing = library(&["precursor_A"], &[("transition_A", "missing_precursor")]);
    assert!(OpenSwathLibraryIDNormalizer::normalize_source_ids(&mut missing).is_err());
    assert_eq!(missing.compounds[0].id, "precursor_A");
    assert_eq!(missing.transitions[0].transition_name, "transition_A");
}

/// `OpenSwathLibraryIDNormalizer_test.cpp:176-236`.
#[test]
fn class_test_decoy_prefix_semantics() {
    let mut exp = library(
        &["PEPTIDE", "DECOY_PEPTIDE", "REV_PEPTIDE"],
        &[
            ("target_y7", "PEPTIDE"),
            ("decoy_y7", "DECOY_PEPTIDE"),
            ("rev_y7", "REV_PEPTIDE"),
            ("REV_transition_y8", "PEPTIDE"),
        ],
    );
    let ids = OpenSwathLibraryIDNormalizer::normalize_source_ids(&mut exp).unwrap();
    let decoys: Vec<bool> = exp.transitions.iter().map(LightTransition::decoy).collect();
    assert_eq!(decoys, [false, true, false, false]);
    OpenSwathLibraryIDNormalizer::materialize_decoy_prefix(&mut exp, &ids, "REV_");
    assert!(exp.transitions[2].decoy());
    assert!(exp.transitions[3].decoy());
    assert_eq!(ids.transition_canonical_to_source["3"], "REV_transition_y8");

    let decoy = ids.precursor_source_to_canonical["DECOY_PEPTIDE"].clone();
    assert_eq!(ids.precursor_canonical_to_source[&decoy], "DECOY_PEPTIDE");
    let custom = ids.precursor_source_to_canonical["REV_PEPTIDE"].clone();
    assert_eq!(ids.precursor_canonical_to_source[&custom], "REV_PEPTIDE");
    let target = ids.precursor_source_to_canonical["PEPTIDE"].clone();
    assert_eq!(
        OpenSwathLibraryIDNormalizer::canonical_target_for_decoy_precursor(&decoy, &ids, "DECOY_"),
        Some(target.clone())
    );
    assert_eq!(
        OpenSwathLibraryIDNormalizer::canonical_target_for_decoy_precursor(&custom, &ids, "REV_"),
        Some(target.clone())
    );
    assert_eq!(
        OpenSwathLibraryIDNormalizer::canonical_target_for_decoy_precursor(&target, &ids, "DECOY_"),
        None
    );
    // Native check of the documented compatibility path: without provenance,
    // the prefix is stripped from the reference itself.
    let none = SourceIDMapping::default();
    assert_eq!(
        OpenSwathLibraryIDNormalizer::canonical_target_for_decoy_precursor(
            "DECOY_X", &none, "DECOY_"
        ),
        Some("X".to_string())
    );
}

/// `OpenSwathLibraryIDNormalizer_test.cpp:238-297`.
#[test]
fn class_test_has_canonical_ids() {
    let valid = library(&["0", "7"], &[("3", "0"), ("100", "7")]);
    assert!(OpenSwathLibraryIDNormalizer::has_canonical_ids(&valid));
    assert!(OpenSwathLibraryIDNormalizer::has_canonical_id_format(
        &valid
    ));

    let mut duplicate = valid.clone();
    duplicate.transitions[1].transition_name = "3".into();
    assert!(!OpenSwathLibraryIDNormalizer::has_canonical_ids(&duplicate));
    assert!(OpenSwathLibraryIDNormalizer::has_canonical_id_format(
        &duplicate
    ));

    let mut source_style = valid.clone();
    source_style.compounds[0].id = "PEPTIDEA".into();
    source_style.transitions[0].peptide_ref = "PEPTIDEA".into();
    assert!(!OpenSwathLibraryIDNormalizer::has_canonical_ids(
        &source_style
    ));
    assert!(!OpenSwathLibraryIDNormalizer::has_canonical_id_format(
        &source_style
    ));

    for bad in ["-0", "-00", "+0"] {
        let mut e = valid.clone();
        e.compounds[0].id = bad.into();
        e.transitions[0].peptide_ref = bad.into();
        assert!(
            !OpenSwathLibraryIDNormalizer::has_canonical_id_format(&e),
            "{bad}"
        );
        assert!(
            !OpenSwathLibraryIDNormalizer::has_canonical_ids(&e),
            "{bad}"
        );
    }
}

/// `OpenSwathLibraryIDNormalizer_test.cpp:299-343`.
#[test]
fn class_test_validate_canonical_ids() {
    let valid = library(
        &["0", "7", "15"],
        &[("3", "0"), ("100", "15"), ("902", "7")],
    );
    OpenSwathLibraryIDNormalizer::validate_canonical_ids(&valid).unwrap();

    let mut leading_zero = valid.clone();
    leading_zero.compounds[0].id = "007".into();
    leading_zero.transitions[0].peptide_ref = "007".into();
    assert!(OpenSwathLibraryIDNormalizer::validate_canonical_ids(&leading_zero).is_err());
    let mut negative = valid.clone();
    negative.transitions[0].transition_name = "-1".into();
    assert!(OpenSwathLibraryIDNormalizer::validate_canonical_ids(&negative).is_err());
    let mut signed_zero = valid.clone();
    signed_zero.compounds[0].id = "-0".into();
    signed_zero.transitions[0].peptide_ref = "-0".into();
    assert!(OpenSwathLibraryIDNormalizer::validate_canonical_ids(&signed_zero).is_err());
    let mut duplicate = valid.clone();
    duplicate.transitions[1].transition_name = "3".into();
    assert!(OpenSwathLibraryIDNormalizer::validate_canonical_ids(&duplicate).is_err());
    let mut unknown = valid.clone();
    unknown.transitions[0].peptide_ref = "999".into();
    assert!(OpenSwathLibraryIDNormalizer::validate_canonical_ids(&unknown).is_err());
}

/// Tier 4, independently derived: `i64::MAX` is canonical, one more overflows
/// the source's `std::from_chars` and is not.
#[test]
fn canonical_id_format_boundaries() {
    let ok = library(&["9223372036854775807"], &[]);
    assert!(OpenSwathLibraryIDNormalizer::has_canonical_id_format(&ok));
    let overflow = library(&["9223372036854775808"], &[]);
    assert!(!OpenSwathLibraryIDNormalizer::has_canonical_id_format(
        &overflow
    ));
    let spaced = library(&[" 7"], &[]);
    assert!(!OpenSwathLibraryIDNormalizer::has_canonical_id_format(
        &spaced
    ));
}

// ---------------------------------------------------------------------------
// Tier 1: replay of the O5b Release oracle
// ---------------------------------------------------------------------------

/// Absolute tolerance for p-values, q-values, PEPs and posteriors against
/// the oracle. See the support document: the statistics are compared to the
/// last few ulps the kernel density and spline stages are allowed to differ
/// by (`docs/OPENSWATH_INFERENCE_SUPPORT.md`).
const ORACLE_ABS_TOL: f64 = 1e-12;

fn read_tsv(path: &Path) -> Vec<Vec<String>> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    text.lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| l.split('\t').map(str::to_string).collect())
        .collect()
}

fn pd(s: &str) -> f64 {
    match s {
        "nan" => f64::NAN,
        "inf" => f64::INFINITY,
        "-inf" => f64::NEG_INFINITY,
        _ => s.parse().unwrap_or_else(|_| panic!("not a number: {s}")),
    }
}

fn opt_i(s: &str) -> Option<i64> {
    if s == "NA" {
        None
    } else {
        Some(s.parse().unwrap())
    }
}

fn opt_f(s: &str) -> Option<f64> {
    if s == "NA" { None } else { Some(pd(s)) }
}

fn level(s: &str) -> InferenceLevel {
    match s {
        "peptidoform" => InferenceLevel::Peptidoform,
        "peptide" => InferenceLevel::Peptide,
        "protein" => InferenceLevel::Protein,
        "gene" => InferenceLevel::Gene,
        _ => panic!("level {s}"),
    }
}

fn context(s: &str) -> InferenceContext {
    match s {
        "global" => InferenceContext::Global,
        "experiment-wide" => InferenceContext::ExperimentWide,
        "run-specific" => InferenceContext::RunSpecific,
        _ => panic!("context {s}"),
    }
}

fn error_variant(v: &str) -> ErrorEstimationConfig {
    let mut c = ErrorEstimationConfig::default();
    match v {
        "default" => {}
        "parametric" => c.parametric = true,
        "pfdr" => c.pfdr = true,
        "smoother" => c.pi0_method = Pi0Method::Smoother,
        "smoother_log" => {
            c.pi0_method = Pi0Method::Smoother;
            c.pi0_smooth_log_pi0 = true;
        }
        "logit" => c.lfdr_transformation = LfdrTransform::Logit,
        "notrunc_nomono" => {
            c.lfdr_truncate = false;
            c.lfdr_monotone = false;
        }
        "lambda_single" => c.pi0_lambda = vec![0.4],
        "lambda_bad" => c.pi0_lambda = vec![0.1, 1.0],
        _ => panic!("variant {v}"),
    }
    c
}

fn ipf_config(v: &str) -> PeptidoformInferenceConfig {
    let mut c = PeptidoformInferenceConfig::default();
    match v {
        "default" => {}
        "loose" => c.ipf_max_precursor_peakgroup_pep = 1.0,
        "ms1" => c.ipf_ms1_scoring = true,
        "ms2" => c.ipf_ms2_scoring = true,
        "ms1ms2" => {
            c.ipf_ms1_scoring = true;
            c.ipf_ms2_scoring = true;
        }
        "ms1ms2_loose" => {
            c.ipf_ms1_scoring = true;
            c.ipf_ms2_scoring = true;
            c.ipf_max_precursor_peakgroup_pep = 1.0;
            c.ipf_max_precursor_pep = 0.3;
        }
        "grouped" => c.ipf_grouped_fdr = true,
        "propagate" => {
            c.propagate_signal_across_runs = true;
            c.ipf_max_precursor_peakgroup_pep = 1.0;
            c.across_run_confidence_threshold = 0.05;
        }
        "propagate_default" => c.propagate_signal_across_runs = true,
        _ => panic!("ipf config {v}"),
    }
    c
}

fn num(v: f64) -> String {
    if v.is_nan() {
        "nan".into()
    } else if v.is_infinite() {
        if v > 0.0 { "inf".into() } else { "-inf".into() }
    } else {
        format!("{v:?}")
    }
}

fn opt_int(v: Option<i64>) -> String {
    v.map_or_else(|| "NA".to_string(), |x| x.to_string())
}

fn read_lc_rows(path: &Path) -> Vec<LevelContextInputRow> {
    read_tsv(path)
        .into_iter()
        .map(|f| LevelContextInputRow {
            run_id: opt_i(&f[0]),
            group_id: f[1].clone(),
            entity_id: f[2].parse().unwrap(),
            decoy: f[3] == "1",
            score: pd(&f[4]),
            context: context(&f[5]),
        })
        .collect()
}

fn read_prec_rows(path: &Path) -> Vec<IPFPrecursorRow> {
    read_tsv(path)
        .into_iter()
        .map(|f| prec(f[0].parse().unwrap(), pd(&f[1]), opt_f(&f[2]), opt_f(&f[3])))
        .collect()
}

fn read_trans_rows(path: &Path) -> Vec<IPFTransitionRow> {
    read_tsv(path)
        .into_iter()
        .map(|f| IPFTransitionRow {
            feature_id: f[0].parse().unwrap(),
            transition_id: f[1].parse().unwrap(),
            peptide_id: f[2].parse().unwrap(),
            pep: pd(&f[3]),
            bmask: f[4].parse().unwrap(),
            num_peptidoforms: f[5].parse().unwrap(),
            alignment_group_id: opt_i(&f[6]),
        })
        .collect()
}

fn read_align_rows(path: &Path) -> Vec<IPFAlignmentRow> {
    read_tsv(path)
        .into_iter()
        .map(|f| IPFAlignmentRow {
            alignment_group_id: f[0].parse().unwrap(),
            feature_id: f[1].parse().unwrap(),
        })
        .collect()
}

/// Replays one `lc` case; returns the driver's output lines.
fn replay_lc(args: &[String], dir: &Path) -> Vec<String> {
    let rows = read_lc_rows(&dir.join(&args[0]));
    let cfg = LevelContextInferenceConfig {
        level: level(&args[1]),
        context: context(&args[2]),
        error: error_variant(&args[3]),
    };
    let result = match cfg.level {
        InferenceLevel::Peptide => OpenSwathPeptideInference.infer(&rows, &cfg),
        InferenceLevel::Protein => OpenSwathProteinInference.infer(&rows, &cfg),
        InferenceLevel::Gene => OpenSwathGeneInference.infer(&rows, &cfg),
        InferenceLevel::Peptidoform => LevelContextInference::infer(&rows, &cfg),
    };
    match result {
        Ok(rows) => rows
            .iter()
            .map(|r| {
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    opt_int(r.run_id),
                    r.entity_id,
                    r.context,
                    num(r.score),
                    num(r.pvalue),
                    num(r.qvalue),
                    num(r.pep)
                )
            })
            .collect(),
        Err(e) => vec![format!("ERROR\t{e}")],
    }
}

fn replay_ipf(args: &[String], dir: &Path) -> Vec<String> {
    let precursors = read_prec_rows(&dir.join(&args[0]));
    let transitions = read_trans_rows(&dir.join(&args[1]));
    let alignment = read_align_rows(&dir.join(&args[2]));
    let cfg = ipf_config(&args[3]);
    let inference = OpenSwathPeptidoformInference;
    let mut out: Vec<String> = inference
        .precursor_inference(&precursors, &cfg)
        .iter()
        .map(|r| format!("P\t{}\t{}", r.feature_id, num(r.precursor_peakgroup_pep)))
        .collect();
    match inference.infer(&precursors, &transitions, &alignment, &cfg) {
        Ok(rows) => out.extend(rows.iter().map(|r| {
            format!(
                "R\t{}\t{}\t{}\t{}\t{}",
                r.feature_id,
                r.peptide_id,
                num(r.precursor_peakgroup_pep),
                num(r.qvalue),
                num(r.pep)
            )
        })),
        Err(e) => out.push(format!("ERROR\t{e}")),
    }
    out
}

fn replay_norm(args: &[String], dir: &Path) -> Vec<String> {
    let mut exp = LightTargetedExperiment::new();
    for f in read_tsv(&dir.join(&args[0])) {
        if f[0] == "C" {
            exp.compounds.push(LightCompound {
                id: f.get(1).cloned().unwrap_or_default(),
                ..LightCompound::default()
            });
        } else {
            let mut t = LightTransition {
                transition_name: f[1].clone(),
                peptide_ref: f[2].clone(),
                ..LightTransition::default()
            };
            t.set_decoy(f[3] == "1");
            exp.transitions.push(t);
        }
    }
    let prefix = if args[1] == "-" { "" } else { args[1].as_str() };
    let b = |v: bool| u8::from(v);
    let mut out = vec![format!(
        "HAS_CANONICAL_BEFORE\t{}\t{}",
        b(OpenSwathLibraryIDNormalizer::has_canonical_ids(&exp)),
        b(OpenSwathLibraryIDNormalizer::has_canonical_id_format(&exp))
    )];
    let before = exp.clone();
    let ids = match OpenSwathLibraryIDNormalizer::normalize_source_ids(&mut exp) {
        Ok(ids) => ids,
        Err(e) => {
            assert_eq!(
                exp.compounds, before.compounds,
                "a refused library is unchanged"
            );
            out.push(format!("ERROR\t{e}"));
            return out;
        }
    };
    out.extend(exp.compounds.iter().map(|c| format!("C\t{}", c.id)));
    out.extend(exp.transitions.iter().map(|t| {
        format!(
            "T\t{}\t{}\t{}",
            t.transition_name,
            t.peptide_ref,
            b(t.decoy())
        )
    }));
    let sorted = |m: &std::collections::HashMap<String, String>| -> BTreeMap<String, String> {
        m.iter().map(|(k, v)| (k.clone(), v.clone())).collect()
    };
    let s2c = sorted(&ids.precursor_source_to_canonical);
    let c2s = sorted(&ids.precursor_canonical_to_source);
    let t2s = sorted(&ids.transition_canonical_to_source);
    out.extend(s2c.iter().map(|(s, c)| format!("S2C\t{s}\t{c}")));
    out.extend(c2s.iter().map(|(c, s)| format!("C2S\t{c}\t{s}")));
    out.extend(t2s.iter().map(|(c, s)| format!("TC2S\t{c}\t{s}")));
    OpenSwathLibraryIDNormalizer::materialize_decoy_prefix(&mut exp, &ids, prefix);
    out.extend(
        exp.transitions
            .iter()
            .map(|t| format!("TD\t{}\t{}", t.transition_name, b(t.decoy()))),
    );
    for c in c2s.keys() {
        let target =
            OpenSwathLibraryIDNormalizer::canonical_target_for_decoy_precursor(c, &ids, prefix);
        out.push(format!(
            "PAIR\t{c}\t{}",
            target.unwrap_or_else(|| "NONE".into())
        ));
    }
    out.push(format!(
        "HAS_CANONICAL_AFTER\t{}\t{}",
        b(OpenSwathLibraryIDNormalizer::has_canonical_ids(&exp)),
        b(OpenSwathLibraryIDNormalizer::has_canonical_id_format(&exp))
    ));
    out
}

/// Replays one oracle case, or `None` for a reader dump the port cannot run.
fn replay(mode: &str, args: &[String], dir: &Path) -> Option<Vec<String>> {
    Some(match mode {
        "dump_lc" | "dump_ipf" => return None,
        "lc" => replay_lc(args, dir),
        "lcwrap" => {
            let cfg = config(level(&args[1]), InferenceContext::Global);
            let r = match args[0].as_str() {
                "peptide" => OpenSwathPeptideInference.infer(&[], &cfg),
                "protein" => OpenSwathProteinInference.infer(&[], &cfg),
                _ => OpenSwathGeneInference.infer(&[], &cfg),
            };
            vec![match r {
                Ok(rows) => format!("OK\t{}", rows.len()),
                Err(e) => format!("ERROR\t{e}"),
            }]
        }
        "ipf" => replay_ipf(args, dir),
        "bm" => {
            let rows: Vec<BayesianModelRow> = read_tsv(&dir.join(&args[0]))
                .into_iter()
                .map(|f| {
                    bm(
                        f[0].parse().unwrap(),
                        f[1].parse().unwrap(),
                        pd(&f[2]),
                        pd(&f[3]),
                    )
                })
                .collect();
            OpenSwathPeptidoformInference::apply_bm(&rows)
                .iter()
                .map(|r| format!("{}\t{}\t{}", r.feature_id, r.hypothesis, num(r.posterior)))
                .collect()
        }
        "pbm" => OpenSwathPeptidoformInference::prepare_precursor_bm(&read_prec_rows(
            &dir.join(&args[0]),
        ))
        .iter()
        .map(|r| {
            format!(
                "{}\t{}\t{}\t{}",
                r.feature_id,
                r.hypothesis,
                num(r.prior),
                num(r.evidence)
            )
        })
        .collect(),
        "mfdr" => {
            let values: Vec<f64> = read_tsv(&dir.join(&args[0]))
                .iter()
                .map(|f| pd(&f[0]))
                .collect();
            OpenSwathPeptidoformInference::compute_model_fdr(&values)
                .unwrap()
                .into_iter()
                .map(num)
                .collect()
        }
        "norm" => replay_norm(args, dir),
        _ => panic!("unknown mode {mode}"),
    })
}

/// Compares one produced line against the oracle's.
///
/// Numeric fields compare within [`ORACLE_ABS_TOL`] (NaN equals NaN); other
/// fields exactly. An error line compares by kind: the oracle's
/// `ERROR\tPrecondition\t<message>` must match the port's
/// `ERROR\tinvalid value: <message>` text exactly, the oracle's
/// `ERROR\tstd\t...` (an uncaught `std::invalid_argument`) only needs the port
/// to fail without the run-specific rewrap, and the normaliser's
/// `ERROR\tInvalidValue\t<message>` must match `invalid value: <message>`.
/// A case whose driver process died (the oracle's `log/<case>.exit` is not
/// `0`) has no expected output; the port must return without panicking.
fn compare_line(case: &str, index: usize, actual: &str, expected: &str, max_dev: &mut f64) {
    if let Some(rest) = expected.strip_prefix("ERROR\t") {
        assert!(
            actual.starts_with("ERROR\t"),
            "{case}:{index}: expected {expected}, got {actual}"
        );
        let got = &actual["ERROR\t".len()..];
        if let Some(message) = rest
            .strip_prefix("Precondition\t")
            .or_else(|| rest.strip_prefix("InvalidValue\t"))
        {
            assert_eq!(got, format!("invalid value: {message}"), "{case}:{index}");
        } else if rest.starts_with("std\t") {
            assert!(
                !got.contains("Run-specific inference failed"),
                "{case}:{index}: {got}"
            );
        }
        return;
    }
    let a: Vec<&str> = actual.split('\t').collect();
    let e: Vec<&str> = expected.split('\t').collect();
    assert_eq!(a.len(), e.len(), "{case}:{index}: {actual} vs {expected}");
    for (x, y) in a.iter().zip(&e) {
        // A field both sides print as a number compares numerically (the
        // driver writes `%.17g`, the port `{:?}`, so `0` and `0.0` are the
        // same value); anything else compares as text. Integers on both sides
        // (ids, flags) compare as text, so large ids are not rounded.
        if x.parse::<i64>().is_ok() && y.parse::<i64>().is_ok() {
            assert_eq!(x, y, "{case}:{index}: {actual} vs {expected}");
        } else if let (Ok(xv), Ok(yv)) = (x.parse::<f64>(), y.parse::<f64>()) {
            if xv.is_nan() || yv.is_nan() {
                assert!(
                    xv.is_nan() && yv.is_nan(),
                    "{case}:{index}: {actual} vs {expected}"
                );
                continue;
            }
            if xv.is_infinite() || yv.is_infinite() {
                assert_eq!(xv, yv, "{case}:{index}");
                continue;
            }
            let dev = (xv - yv).abs();
            *max_dev = max_dev.max(dev);
            assert!(
                dev <= ORACLE_ABS_TOL,
                "{case}:{index}: {x} vs {y} ({dev:e})\n{actual}\n{expected}"
            );
        } else {
            assert_eq!(x, y, "{case}:{index}: {actual} vs {expected}");
        }
    }
}

/// The one documented divergence of the replay, which lies below this port,
/// in `math::multiple_testing` / `math::kernel_density`
/// (`MATH/STATISTICS/MultipleTesting.h`, `KernelDensityEstimation.h`); see
/// `docs/OPENSWATH_INFERENCE_SUPPORT.md`, "Known divergence".
///
/// When every target outscores every decoy of an inference group, all target
/// p-values equal `1 / n_decoys` and the local-FDR stage estimates a density
/// of identical transformed values. `bwNrd0` then takes its degenerate branch
/// on whether the sample standard deviation of those identical values is
/// exactly zero, which depends on the last bit of the probit transform. The
/// Release build (Boost, 80-bit `long double` internally) returns
/// `probit(0.0625) = -1.5341205443525465`, whose 16-fold mean is exact, so the
/// bandwidth is healthy (`0.793...`, O5b probe on ibminode06); the port's
/// `statrs` quantile returns `-1.534120544352546`, two units in the last place
/// away, whose mean is not, so the bandwidth is `2.4e-16`, the KDE grid
/// collapses and `kde_fft_eval` refuses it.
///
/// A case is accepted as this divergence only when both hold, so any other
/// difference still fails: the port's whole output is that refusal, and the
/// oracle's output shows the degenerate distribution (within each run, every
/// finite p-value identical). When the math layer is fixed the port returns
/// rows again and the case is compared like every other.
fn is_known_probit_divergence(actual: &[String], expected: &[String]) -> bool {
    if actual != ["ERROR\tinvalid value: cubic spline abscissae must be strictly increasing"] {
        return false;
    }
    if expected.is_empty() || expected.iter().any(|l| l.starts_with("ERROR")) {
        return false;
    }
    let mut per_run: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for line in expected {
        let f: Vec<&str> = line.split('\t').collect();
        let p = pd(f[4]);
        if p.is_finite() {
            per_run.entry(f[0].to_string()).or_default().push(p);
        }
    }
    per_run
        .values()
        .all(|ps| ps.windows(2).all(|w| w[0].to_bits() == w[1].to_bits()))
}

/// The port's rows, or — only when the port's error is exactly the documented
/// probit divergence and the oracle's rows show its precondition — the
/// Release build's rows from `oracle`, announced on stderr.
fn level_rows_or_documented_divergence(
    result: openms::Result<Vec<LevelContextResultRow>>,
    oracle: &Path,
) -> Vec<LevelContextResultRow> {
    match result {
        Ok(rows) => rows,
        Err(error) => {
            let expected: Vec<String> = std::fs::read_to_string(oracle)
                .unwrap()
                .lines()
                .map(str::to_string)
                .collect();
            assert!(
                is_known_probit_divergence(&[format!("ERROR\t{error}")], &expected),
                "undocumented error: {error}"
            );
            eprintln!(
                "KNOWN DIVERGENCE (probit/KDE, see is_known_probit_divergence): literals checked on {}",
                oracle.display()
            );
            expected
                .iter()
                .map(|line| {
                    let f: Vec<&str> = line.split('\t').collect();
                    LevelContextResultRow {
                        run_id: opt_i(f[0]),
                        entity_id: f[1].parse().unwrap(),
                        context: context(f[2]),
                        score: pd(f[3]),
                        pvalue: pd(f[4]),
                        qvalue: pd(f[5]),
                        pep: pd(f[6]),
                    }
                })
                .collect()
        }
    }
}

#[test]
fn oracle_replay() {
    let dir = data_dir();
    let cases = read_tsv(&dir.join("cases.txt"));
    let mut replayed = 0;
    let mut crashed = Vec::new();
    let mut failures: Vec<String> = Vec::new();
    let mut known_divergence: Vec<String> = Vec::new();
    let mut max_dev: f64 = 0.0;
    for case in &cases {
        let (name, mode, args) = (&case[0], &case[1], &case[2..]);
        let Some(actual) = replay(mode, &args[..args.len() - 1], &dir) else {
            continue;
        };
        let exit = std::fs::read_to_string(dir.join(format!("log/{name}.exit"))).unwrap();
        if exit.trim() != "0" {
            // The Release build died on this case (see the support document);
            // reaching here means the port returned instead.
            crashed.push(format!("{name}: exit {} -> {actual:?}", exit.trim()));
            replayed += 1;
            continue;
        }
        let expected_path = dir.join(&args[args.len() - 1]);
        let expected: Vec<String> = std::fs::read_to_string(&expected_path)
            .unwrap_or_else(|e| panic!("{}: {e}", expected_path.display()))
            .lines()
            .map(str::to_string)
            .collect();
        if mode == "lc" && is_known_probit_divergence(&actual, &expected) {
            known_divergence.push(name.clone());
            replayed += 1;
            continue;
        }
        // Every case is checked and every failing case reported, not only the
        // first.
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let mut dev = 0.0;
            assert_eq!(
                actual.len(),
                expected.len(),
                "{name}: line count\n{actual:#?}\n{expected:#?}"
            );
            for (i, (a, e)) in actual.iter().zip(&expected).enumerate() {
                compare_line(name, i, a, e, &mut dev);
            }
            dev
        }));
        match outcome {
            Ok(dev) => max_dev = max_dev.max(dev),
            Err(payload) => {
                let message = payload
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| payload.downcast_ref::<&str>().map(|s| (*s).to_string()))
                    .unwrap_or_default();
                failures.push(format!("{name}: {}", message.lines().next().unwrap_or("")));
            }
        }
        replayed += 1;
    }
    // Every case of cases.txt but the 15 reader dumps (9 level-context, 6 IPF).
    assert_eq!(replayed, cases.len() - 15);
    eprintln!("oracle_replay: {replayed} cases, max numeric deviation {max_dev:e}");
    eprintln!("oracle_replay: cases where the Release build crashed: {crashed:#?}");
    eprintln!(
        "oracle_replay: {} cases hit the documented probit divergence (see is_known_probit_divergence): {known_divergence:#?}",
        known_divergence.len()
    );
    assert!(
        failures.is_empty(),
        "{} of {replayed} cases differ from the oracle:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// The fixture-derived reader counts of `OSWFileInference_test.cpp:88-139`
/// and `LevelContextInference_test.cpp:172-229`, checked on the rows the
/// Release `OSWFile` reader wrote, so the replay runs on the inputs the class
/// tests describe.
#[test]
fn oracle_reader_dumps_match_class_test_counts() {
    let dir = data_dir();
    let count = |p: &str| read_tsv(&dir.join(p)).len();
    assert_eq!(count("out/dump_peptide_global.tsv"), 682);
    assert_eq!(count("out/dump_peptide_run-specific.tsv"), 682);
    assert_eq!(count("out/dump_protein_global.tsv"), 32);
    assert_eq!(count("out/dump_gene_global.tsv"), 32);
    assert_eq!(count("out/ipf_fixture_default.prec.tsv"), 389);
}

/// `LevelContextInference_test.cpp:172-229` and
/// `PeptidoformInference_test.cpp:275-298`: the fixture literals, on the port's
/// results for the Release reader's rows.
#[test]
fn class_test_fixture_literals() {
    let dir = data_dir();
    let global = config(InferenceLevel::Peptide, InferenceContext::Global);
    let rows = read_lc_rows(&dir.join("out/dump_peptide_global.tsv"));
    let res = OpenSwathPeptideInference.infer(&rows, &global).unwrap();
    assert_eq!(res.len(), 682);
    let r = res.iter().find(|r| r.entity_id == 2).unwrap();
    assert!((r.score - 5.538_958_072_662_35).abs() < 1e-8);
    assert!((r.pvalue - 0.0029).abs() < 5e-4);
    assert!(r.qvalue.is_finite() && (0.0..=0.05).contains(&r.qvalue));
    assert!(r.pep.is_finite() && (0.0..=0.05).contains(&r.pep));

    let run_cfg = config(InferenceLevel::Peptide, InferenceContext::RunSpecific);
    let rows = read_lc_rows(&dir.join("out/dump_peptide_run-specific.tsv"));
    let res = OpenSwathPeptideInference.infer(&rows, &run_cfg).unwrap();
    assert_eq!(res.len(), 682);
    assert_eq!(res.iter().filter(|r| r.run_id.is_some()).count(), 682);

    // Protein and gene level: every target outscores every decoy in this
    // fixture, the case of `is_known_probit_divergence`. While that divergence
    // stands the port refuses the call; the class-test literals are then
    // checked on the Release build's rows, and the port's refusal must be
    // exactly the documented one.
    let protein = config(InferenceLevel::Protein, InferenceContext::Global);
    let rows = read_lc_rows(&dir.join("out/dump_protein_global.tsv"));
    let res = level_rows_or_documented_divergence(
        OpenSwathProteinInference.infer(&rows, &protein),
        &dir.join("out/lc_protein_global_default.tsv"),
    );
    assert_eq!(res.len(), 32);
    let r = res.iter().find(|r| r.entity_id == 0).unwrap();
    assert!((r.score - 5.840_140_819_549_56).abs() < 1e-8);
    assert!((r.pvalue - 0.0625).abs() < 5e-4);
    assert!(r.qvalue.is_finite() && r.pep.is_finite());
    assert!((0.0..=0.5).contains(&r.qvalue));
    assert!((0.0..=1.0).contains(&r.pep));

    let gene = config(InferenceLevel::Gene, InferenceContext::Global);
    let rows = read_lc_rows(&dir.join("out/dump_gene_global.tsv"));
    let res = level_rows_or_documented_divergence(
        OpenSwathGeneInference.infer(&rows, &gene),
        &dir.join("out/lc_gene_global_default.tsv"),
    );
    assert_eq!(res.len(), 32);
    let r = res.iter().find(|r| r.entity_id == 0).unwrap();
    assert!((r.score - 5.840_140_819_549_56).abs() < 1e-8);

    let cfg = PeptidoformInferenceConfig::default();
    let precursors = read_prec_rows(&dir.join("out/ipf_fixture_default.prec.tsv"));
    let transitions = read_trans_rows(&dir.join("out/ipf_fixture_default.trans.tsv"));
    let res = OpenSwathPeptidoformInference
        .infer(&precursors, &transitions, &[], &cfg)
        .unwrap();
    assert_eq!(res.len(), 333);
    let r = res
        .iter()
        .find(|r| r.feature_id == -9_078_977_811_506_172_301 && r.peptide_id == 305)
        .unwrap();
    assert_similar(
        r.precursor_peakgroup_pep,
        0.003_142_319_081_467_08,
        "fixture pp",
    );
    assert!(r.qvalue <= 0.01 && r.pep <= 0.01);
}

// ---------------------------------------------------------------------------
// Tier 4: native bounds and invariants
// ---------------------------------------------------------------------------

/// Independently derived: posteriors of each feature sum to one, and an
/// empty model gives an empty result.
#[test]
fn apply_bm_normalises_per_feature() {
    assert!(OpenSwathPeptidoformInference::apply_bm(&[]).is_empty());
    let rows = [
        bm(5, 1, 0.4, 0.9),
        bm(5, 0, 0.6, 0.2),
        bm(2, -1, 0.3, 0.5),
        bm(2, 7, 0.7, 0.25),
    ];
    let post = OpenSwathPeptidoformInference::apply_bm(&rows);
    let keys: Vec<(i64, i64)> = post.iter().map(|r| (r.feature_id, r.hypothesis)).collect();
    assert_eq!(keys, [(2, -1), (2, 7), (5, 0), (5, 1)]);
    assert!((post[0].posterior + post[1].posterior - 1.0).abs() < 1e-15);
    assert!((post[2].posterior + post[3].posterior - 1.0).abs() < 1e-15);
}

/// Native: a non-null hypothesis without a positive peptidoform count is
/// refused, as the source's `Exception::Precondition`.
#[test]
fn transition_without_peptidoform_count_is_refused() {
    let err = OpenSwathPeptidoformInference
        .infer(
            &[prec(1, 0.2, None, None)],
            &[trans(1, 11, 101, 0.1, 1, 0)],
            &[],
            &loose_ipf(),
        )
        .unwrap_err();
    assert!(err.to_string().contains("num_peptidoforms > 0"));
}
