// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Executed differential (evidence tier 1) for the MRM feature ports:
//! `MRMFeatureQCFile`, `MRMFeaturePickerFile`, `MRMFeatureFilter`,
//! `MRMFeatureSelectorScore`, `MRMBatchFeatureSelector`, `ConfidenceScoring`
//! and `MRMFeatureFinderScoring` against the C++ Release build.
//!
//! The fixtures under `tests/data/mrm_oracle/` are the output of the driver
//! `../oracle/O3b-mrm/driver.cpp`, built with `x86_64-conda-linux-gnu-g++`
//! `-O3 -DNDEBUG` against the Linux x86-64 Release install
//! `openms4-release-bc9cc12-c19e494-174b576` and run twice on ibminode06 with
//! byte-identical output (`tests/data/mrm_feature_provenance.json`). Each record
//! is `key<TAB>value`, doubles as `%.17g`; this test recomputes every record with
//! the port (the driver's construction of inputs is mirrored here statement for
//! statement) and compares keys and values. Numbers compare as `f64`, exactly,
//! unless a class below says otherwise; strings compare as strings.
#![cfg(all(feature = "mzml", feature = "featurexml"))]

use openms::analysis::openswath::confidence_scoring::ConfidenceScoring;
use openms::analysis::openswath::data_access_helper::convert_targeted_exp;
use openms::analysis::openswath::mrm_batch_feature_selector::batch_mrm_features_qmip;
use openms::analysis::openswath::mrm_feature_filter::{FlagOrFilter, MRMFeatureFilter};
use openms::analysis::openswath::mrm_feature_finder_scoring::{
    MRMFeatureFinderScoring, TransitionGroupMapType,
};
use openms::analysis::openswath::mrm_feature_picker::{ComponentGroupParams, ComponentParams};
use openms::analysis::openswath::mrm_feature_qc::{ComponentGroupQCs, ComponentQCs, MRMFeatureQC};
use openms::analysis::openswath::mrm_feature_selector::{
    LambdaScore, LinearProgram, LpSolution, LpSolver, MRMFeatureSelector, MRMFeatureSelectorScore,
    PartitionSolver, SelectorParameters, SolverStatus, VariableType, compute_score,
};
use openms::analysis::openswath::reaction_monitoring_transition::ReactionMonitoringTransition;
use openms::analysis::openswath::spectrum_access_openms::SimpleOpenMSSpectraFactory;
use openms::analysis::openswath::swath_map::SwathMap;
use openms::analysis::openswath::targeted_experiment::TargetedExperiment;
use openms::analysis::openswath::targeted_experiment_helper::{
    Peptide, Protein, RTType, RetentionTime,
};
use openms::analysis::transformations::TransformationDescription;
use openms::concept::unique_id::UniqueIdGenerator;
use openms::format::featurexml;
use openms::format::mrm_feature_picker_file::MRMFeaturePickerFile;
use openms::format::mrm_feature_qc_file::MRMFeatureQCFile;
use openms::kernel::{Feature, FeatureMap, MSExperiment};
use openms::metadata::{MetaInfo, MetaValueData};
use openms::param::{Param, ParamValue};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn oracle(name: &str) -> PathBuf {
    data(&format!("mrm_oracle/{name}"))
}

/// A produced record: a number or a string.
enum Val {
    F(f64),
    S(String),
}

#[derive(Default)]
struct Out {
    records: BTreeMap<String, Val>,
}

impl Out {
    fn p(&mut self, key: &str, value: f64) {
        self.records.insert(key.to_owned(), Val::F(value));
    }
    fn pi(&mut self, key: &str, value: i64) {
        self.p(key, value as f64);
    }
    fn ps(&mut self, key: &str, value: &str) {
        self.records
            .insert(key.to_owned(), Val::S(value.to_owned()));
    }
}

fn load_oracle(name: &str) -> BTreeMap<String, String> {
    let text = std::fs::read_to_string(oracle(&format!("{name}.stdout"))).unwrap();
    let mut records = BTreeMap::new();
    for line in text.lines() {
        // the library's INFO log shares stdout; records have exactly one tab
        if let Some((key, value)) = line.split_once('\t') {
            if !value.contains('\t') && key != "RUN" {
                records.insert(key.to_owned(), value.to_owned());
            }
        }
    }
    records
}

/// The oracle records of `name` whose key starts with one of `prefixes`.
fn load_oracle_prefixed(name: &str, prefixes: &[&str]) -> BTreeMap<String, String> {
    load_oracle(name)
        .into_iter()
        .filter(|(k, _)| prefixes.iter().any(|p| k.starts_with(p)))
        .collect()
}

/// Compare a produced set against the retained oracle; `loose` maps a key
/// prefix to a relative tolerance for numbers under it.
fn compare(name: &str, produced: &Out, loose: &[(&str, f64)]) {
    compare_with(name, load_oracle(name), produced, loose);
}

fn compare_with(
    name: &str,
    expected: BTreeMap<String, String>,
    produced: &Out,
    loose: &[(&str, f64)],
) {
    let mut problems: Vec<String> = Vec::new();
    for (key, value) in &expected {
        match produced.records.get(key) {
            None => problems.push(format!("missing {key} (oracle {value})")),
            Some(Val::S(s)) => {
                if s != value {
                    problems.push(format!("{key}: '{s}' != oracle '{value}'"));
                }
            }
            Some(Val::F(x)) => match value.parse::<f64>() {
                Ok(e) => {
                    let tolerance = loose
                        .iter()
                        .find(|(prefix, _)| key.starts_with(prefix))
                        .map_or(0.0, |(_, t)| *t);
                    let same = (x.is_nan() && e.is_nan())
                        || x.to_bits() == e.to_bits()
                        || x == &e
                        || (tolerance > 0.0 && (x - e).abs() <= tolerance * e.abs().max(x.abs()));
                    if !same {
                        problems.push(format!("{key}: {x:?} != oracle {e:?}"));
                    }
                }
                Err(_) => problems.push(format!("{key}: number {x} vs oracle text '{value}'")),
            },
        }
    }
    for key in produced.records.keys() {
        if !expected.contains_key(key) {
            problems.push(format!("extra {key}"));
        }
    }
    assert!(
        problems.is_empty(),
        "{name}: {} difference(s), first 25:\n{}",
        problems.len(),
        problems
            .iter()
            .take(25)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
    assert!(!expected.is_empty(), "{name}: the oracle fixture is empty");
}

// ------------------------------------------------------------------ dumpers

fn dump_meta(out: &mut Out, pre: &str, meta: &MetaInfo) {
    for (key, value) in meta {
        let k = format!("{pre}.{key}");
        match value.data() {
            MetaValueData::Float(x) => out.p(&k, *x),
            MetaValueData::Integer(i) => out.pi(&k, *i),
            MetaValueData::String(s) => out.ps(&k, s),
            MetaValueData::StringList(l) => {
                out.pi(&format!("{k}.n"), l.len() as i64);
                for (i, s) in l.iter().enumerate() {
                    out.ps(&format!("{k}.{i}"), s);
                }
            }
            MetaValueData::FloatList(l) => {
                out.pi(&format!("{k}.n"), l.len() as i64);
                for (i, x) in l.iter().enumerate() {
                    out.p(&format!("{k}.{i}"), *x);
                }
            }
            MetaValueData::IntegerList(l) => {
                out.pi(&format!("{k}.n"), l.len() as i64);
                for (i, x) in l.iter().enumerate() {
                    out.pi(&format!("{k}.{i}"), *x);
                }
            }
            MetaValueData::Empty => out.ps(&k, "<empty>"),
        }
    }
}

fn dump_feature(out: &mut Out, pre: &str, f: &Feature, with_subordinates: bool) {
    out.p(&format!("{pre}.rt"), f.rt);
    out.p(&format!("{pre}.mz"), f.mz);
    out.p(&format!("{pre}.intensity"), f64::from(f.intensity));
    out.p(&format!("{pre}.quality"), f64::from(f.quality));
    out.pi(&format!("{pre}.charge"), i64::from(f.charge));
    out.pi(&format!("{pre}.nhulls"), f.convex_hulls.len() as i64);
    dump_meta(out, &format!("{pre}.meta"), &f.metadata);
    if with_subordinates {
        out.pi(&format!("{pre}.nsub"), f.subordinates.len() as i64);
        for (i, s) in f.subordinates.iter().enumerate() {
            dump_feature(out, &format!("{pre}.sub{i}"), s, false);
        }
    }
}

fn dump_ids(out: &mut Out, pre: &str, f: &Feature) {
    let ids = &f.peptide_identifications;
    out.pi(&format!("{pre}.npepid"), ids.len() as i64);
    for (i, id) in ids.iter().enumerate() {
        let p = format!("{pre}.pepid{i}");
        out.ps(&format!("{p}.identifier"), &id.identifier);
        out.pi(&format!("{p}.nhits"), id.hits.len() as i64);
        for (j, h) in id.hits.iter().enumerate() {
            let q = format!("{p}.hit{j}");
            out.p(&format!("{q}.score"), h.score);
            out.pi(&format!("{q}.charge"), i64::from(h.charge));
            out.ps(&format!("{q}.sequence"), &h.sequence.to_string());
            out.pi(&format!("{q}.nev"), h.evidences.len() as i64);
            for (k, e) in h.evidences.iter().enumerate() {
                out.ps(&format!("{q}.ev{k}"), &e.protein_accession);
            }
        }
    }
}

fn dump_qc(out: &mut Out, pre: &str, qc: &MRMFeatureQC) {
    out.pi(&format!("{pre}.ncomp"), qc.component_qcs.len() as i64);
    for (i, c) in qc.component_qcs.iter().enumerate() {
        let p = format!("{pre}.c{i}");
        out.ps(&format!("{p}.name"), &c.component_name);
        out.p(&format!("{p}.rt_l"), c.retention_time_l);
        out.p(&format!("{p}.rt_u"), c.retention_time_u);
        out.p(&format!("{p}.int_l"), c.intensity_l);
        out.p(&format!("{p}.int_u"), c.intensity_u);
        out.p(&format!("{p}.oq_l"), c.overall_quality_l);
        out.p(&format!("{p}.oq_u"), c.overall_quality_u);
        for (k, (l, u)) in &c.meta_value_qc {
            out.p(&format!("{p}.m.{k}.l"), *l);
            out.p(&format!("{p}.m.{k}.u"), *u);
        }
    }
    out.pi(
        &format!("{pre}.ngroup"),
        qc.component_group_qcs.len() as i64,
    );
    for (i, c) in qc.component_group_qcs.iter().enumerate() {
        let p = format!("{pre}.g{i}");
        out.ps(&format!("{p}.name"), &c.component_group_name);
        out.p(&format!("{p}.rt_l"), c.retention_time_l);
        out.p(&format!("{p}.rt_u"), c.retention_time_u);
        out.p(&format!("{p}.int_l"), c.intensity_l);
        out.p(&format!("{p}.int_u"), c.intensity_u);
        out.p(&format!("{p}.oq_l"), c.overall_quality_l);
        out.p(&format!("{p}.oq_u"), c.overall_quality_u);
        for (name, l, u) in [
            ("nheavy", c.n_heavy_l, c.n_heavy_u),
            ("nlight", c.n_light_l, c.n_light_u),
            ("ndet", c.n_detecting_l, c.n_detecting_u),
            ("nquant", c.n_quantifying_l, c.n_quantifying_u),
            ("nident", c.n_identifying_l, c.n_identifying_u),
            ("ntrans", c.n_transitions_l, c.n_transitions_u),
        ] {
            out.pi(&format!("{p}.{name}_l"), i64::from(l));
            out.pi(&format!("{p}.{name}_u"), i64::from(u));
        }
        out.ps(&format!("{p}.pair1"), &c.ion_ratio_pair_name_1);
        out.ps(&format!("{p}.pair2"), &c.ion_ratio_pair_name_2);
        out.p(&format!("{p}.ir_l"), c.ion_ratio_l);
        out.p(&format!("{p}.ir_u"), c.ion_ratio_u);
        out.ps(&format!("{p}.irname"), &c.ion_ratio_feature_name);
        for (k, (l, u)) in &c.meta_value_qc {
            out.p(&format!("{p}.m.{k}.l"), *l);
            out.p(&format!("{p}.m.{k}.u"), *u);
        }
    }
}

fn dump_map(out: &mut Out, pre: &str, fm: &FeatureMap) {
    out.pi(&format!("{pre}.n"), fm.len() as i64);
    for (i, f) in fm.features.iter().enumerate() {
        dump_feature(out, &format!("{pre}.f{i}"), f, true);
    }
}

fn dump_selected(out: &mut Out, pre: &str, fm: &FeatureMap) {
    out.pi(&format!("{pre}.n"), fm.len() as i64);
    for (i, f) in fm.features.iter().enumerate() {
        let p = format!("{pre}.f{i}");
        out.ps(
            &format!("{p}.peptide_ref"),
            f.metadata["PeptideRef"].as_str().unwrap(),
        );
        out.pi(&format!("{p}.id"), f.unique_id as i64);
        out.p(&format!("{p}.rt"), f.rt);
        out.pi(&format!("{p}.nsub"), f.subordinates.len() as i64);
        for (j, s) in f.subordinates.iter().enumerate() {
            out.ps(
                &format!("{p}.sub{j}.native_id"),
                s.metadata["native_id"].as_str().unwrap(),
            );
            out.p(&format!("{p}.sub{j}.rt"), s.rt);
        }
    }
}

fn dump_param(out: &mut Out, pre: &str, param: &Param) {
    for item in param.iter().unwrap() {
        let key = format!("{pre}.{}", item.key);
        match &item.entry.value {
            ParamValue::Float(x) => out.p(&key, *x),
            ParamValue::Integer(i) => out.pi(&key, *i),
            ParamValue::String(s) => out.ps(&key, s),
            other => out.ps(&key, &format!("{other:?}")),
        }
    }
}

// ------------------------------------------------------------------ cases

#[test]
fn qcfile_matches_the_release_build() {
    let file = MRMFeatureQCFile;
    let mut qc = MRMFeatureQC::default();
    file.load(data("MRMFeatureQCFile_1.csv"), &mut qc, false)
        .unwrap();
    file.load(data("MRMFeatureQCFile_2.csv"), &mut qc, true)
        .unwrap();
    let mut out = Out::default();
    dump_qc(&mut out, "qc", &qc);
    compare("qcfile", &out, &[]);
    // the stored files are byte-identical to the Release build's
    let dir = std::env::temp_dir().join(format!("openms_mrm_oracle_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    file.store(dir.join("components.csv"), &qc, false).unwrap();
    file.store(dir.join("groups.csv"), &qc, true).unwrap();
    for (mine, theirs) in [
        ("components.csv", "qcfile_components.csv"),
        ("groups.csv", "qcfile_groups.csv"),
    ] {
        let a = std::fs::read(dir.join(mine)).unwrap();
        let b = std::fs::read(oracle(theirs)).unwrap();
        assert_eq!(
            String::from_utf8_lossy(&a),
            String::from_utf8_lossy(&b),
            "{mine}"
        );
    }
}

#[test]
fn pickerfile_matches_the_release_build() {
    let (cp, cgp): (Vec<ComponentParams>, Vec<ComponentGroupParams>) = MRMFeaturePickerFile
        .load(data("MRMFeaturePickerFile.csv"))
        .unwrap();
    let mut out = Out::default();
    out.pi("ncp", cp.len() as i64);
    out.pi("ncgp", cgp.len() as i64);
    for (i, c) in cp.iter().enumerate() {
        out.ps(&format!("cp{i}.name"), &c.component_name);
        out.ps(&format!("cp{i}.group"), &c.component_group_name);
        dump_param(&mut out, &format!("cp{i}"), &c.params);
    }
    for (i, c) in cgp.iter().enumerate() {
        out.ps(&format!("cgp{i}.name"), &c.component_group_name);
        dump_param(&mut out, &format!("cgp{i}"), &c.params);
    }
    compare("pickerfile", &out, &[]);
}

fn load_small() -> FeatureMap {
    featurexml::load(data("MRMFeatureSelector_100ug.featureXML")).unwrap()
}

fn load_big() -> FeatureMap {
    featurexml::load(data(
        "MRMFeatureSelector_150601_0_BloodProject01_PLT_QC_Broth-1_1_reduced.featureXML",
    ))
    .unwrap()
}

fn peptide_ref(f: &Feature) -> String {
    f.metadata["PeptideRef"].as_str().unwrap().to_owned()
}

fn native_id(f: &Feature) -> String {
    f.metadata["native_id"].as_str().unwrap().to_owned()
}

fn build_criteria(fm: &FeatureMap) -> MRMFeatureQC {
    let mut qc = MRMFeatureQC::default();
    for (i, f) in fm.features.iter().enumerate() {
        let mut cg = ComponentGroupQCs {
            component_group_name: peptide_ref(f),
            retention_time_l: 0.0,
            retention_time_u: if i % 2 == 0 { 1e6 } else { f.rt - 0.5 },
            intensity_l: if i % 3 == 0 { 1e12 } else { 0.0 },
            intensity_u: 1e12,
            overall_quality_l: 0.0,
            overall_quality_u: 1e12,
            n_transitions_l: 0,
            n_transitions_u: if i % 4 == 1 { 1 } else { 100 },
            ..ComponentGroupQCs::default()
        };
        if f.subordinates.len() >= 2 {
            cg.ion_ratio_pair_name_1 = native_id(&f.subordinates[0]);
            cg.ion_ratio_pair_name_2 = native_id(&f.subordinates[1]);
            cg.ion_ratio_l = 0.5;
            cg.ion_ratio_u = 2.0;
            cg.ion_ratio_feature_name = "peak_apex_int".into();
        }
        cg.meta_value_qc.insert(
            "sn_ratio".into(),
            (if i % 2 == 1 { 1e5 } else { 0.0 }, 1e12),
        );
        cg.meta_value_qc
            .insert("peak_apices_sum".into(), (0.0, 1e12));
        qc.component_group_qcs.push(cg);
        for (j, s) in f.subordinates.iter().enumerate() {
            let mut c = ComponentQCs {
                component_name: native_id(s),
                retention_time_l: s.rt - 0.1,
                retention_time_u: if (i + j) % 3 == 0 {
                    s.rt - 0.05
                } else {
                    s.rt + 0.1
                },
                intensity_l: 0.0,
                intensity_u: 1e12,
                overall_quality_l: 0.0,
                overall_quality_u: 1e12,
                ..ComponentQCs::default()
            };
            c.meta_value_qc.insert(
                "peak_apex_int".into(),
                (1.0, if (i + j) % 2 == 0 { 1e6 } else { 10.0 }),
            );
            c.meta_value_qc.insert("total_xic".into(), (0.0, 1e12));
            qc.component_qcs.push(c);
        }
    }
    qc
}

fn scale_meta(f: &mut Feature, keys: &[&str], factor: f64) {
    for key in keys {
        if let Some(v) = f.metadata.get(*key) {
            let x = v.as_f64().unwrap() * factor;
            f.metadata.insert(
                (*key).to_owned(),
                openms::metadata::MetaValue::try_from(x).unwrap(),
            );
        }
    }
}

fn sample(base: &FeatureMap, k: i32) -> FeatureMap {
    let mut fm = base.clone();
    let factor = 1.0 + 0.1 * f64::from(k);
    for f in &mut fm.features {
        f.intensity = (f64::from(f.intensity) * factor) as f32;
        f.rt += 0.05 * f64::from(k);
        scale_meta(f, &["sn_ratio", "peak_apices_sum"], factor);
        for s in &mut f.subordinates {
            s.intensity = (f64::from(s.intensity) * factor) as f32;
            s.rt += 0.05 * f64::from(k);
            scale_meta(s, &["peak_apex_int", "total_xic"], factor);
        }
    }
    fm
}

#[test]
fn filter_matches_the_release_build() {
    let base = load_small();
    let criteria = build_criteria(&base);
    let transitions = TargetedExperiment::new();
    let mut out = Out::default();
    for (name, mode) in [
        ("flag", FlagOrFilter::Flag),
        ("filter", FlagOrFilter::Filter),
    ] {
        let mut fm = base.clone();
        let filt = MRMFeatureFilter {
            flag_or_filter: mode,
        };
        filt.filter_feature_map(&mut fm, &criteria, &transitions)
            .unwrap();
        dump_map(&mut out, name, &fm);
    }
    let samples: Vec<FeatureMap> = (0..3).map(|k| sample(&base, k)).collect();
    let filt = MRMFeatureFilter::default();
    let mut t1 = criteria.clone();
    filt.estimate_default_mrm_feature_qc_values(&samples, &mut t1, &transitions, true)
        .unwrap();
    dump_qc(&mut out, "estimate_init", &t1);
    let mut t1b = criteria.clone();
    filt.estimate_default_mrm_feature_qc_values(&samples, &mut t1b, &transitions, false)
        .unwrap();
    dump_qc(&mut out, "estimate_noinit", &t1b);
    let mut t2 = criteria.clone();
    filt.estimate_perc_rsd(&samples, &mut t2, &transitions)
        .unwrap();
    dump_qc(&mut out, "percrsd", &t2);
    let mut t3 = criteria.clone();
    filt.estimate_background_interferences(&samples, &mut t3, &transitions)
        .unwrap();
    dump_qc(&mut out, "background", &t3);
    {
        let mut fm = base.clone();
        filt.filter_feature_map_perc_rsd(&mut fm, &criteria, &t2)
            .unwrap();
        dump_map(&mut out, "percrsd_flag", &fm);
    }
    {
        let mut fm = base.clone();
        filt.filter_feature_map_background_interference(&mut fm, &criteria, &t3)
            .unwrap();
        dump_map(&mut out, "background_flag", &fm);
    }
    {
        let f = &base.features[0];
        if f.subordinates.len() >= 2 {
            let (a, b) = (&f.subordinates[0], &f.subordinates[1]);
            out.p(
                "ionratio.intensity",
                filt.calculate_ion_ratio(a, b, "intensity").unwrap(),
            );
            out.p(
                "ionratio.apex",
                filt.calculate_ion_ratio(a, b, "peak_apex_int").unwrap(),
            );
            out.p("rtdiff", filt.calculate_rt_difference(a, b));
            out.p("resolution", filt.calculate_resolution(a, b).unwrap());
        }
    }
    compare("filter", &out, &[]);
}

fn weights(pairs: &[(&str, LambdaScore)]) -> BTreeMap<String, LambdaScore> {
    pairs.iter().map(|(k, v)| ((*k).to_owned(), *v)).collect()
}

fn selector_params(
    window: i32,
    step: i32,
    variable_type: VariableType,
    group: bool,
    locality: bool,
    nn: i32,
    w: &BTreeMap<String, LambdaScore>,
) -> SelectorParameters {
    SelectorParameters {
        nn_threshold: nn,
        locality_weight: locality,
        select_transition_group: group,
        segment_window_length: window,
        segment_step_length: step,
        variable_type,
        optimal_threshold: 0.5,
        score_weights: w.clone(),
    }
}

#[test]
fn selector_score_matches_the_release_build() {
    use LambdaScore::{InverseLog, InverseLog10, Linear, Log};
    use VariableType::{Continuous, Integer};
    let big = load_big();
    let small = load_small();
    let logw = weights(&[("sn_ratio", Log), ("peak_apices_sum", Log)]);
    let linw = weights(&[("sn_ratio", Linear), ("peak_apices_sum", Linear)]);
    let invw = weights(&[("sn_ratio", InverseLog), ("peak_apices_sum", InverseLog10)]);
    let none = weights(&[]);
    let runs: Vec<(&str, &FeatureMap, SelectorParameters)> = vec![
        (
            "big_class",
            &big,
            selector_params(-1, -1, Integer, true, false, 4, &logw),
        ),
        (
            "big_win",
            &big,
            selector_params(8, 4, Continuous, true, false, 4, &linw),
        ),
        (
            "big_comp",
            &big,
            selector_params(-1, -1, Integer, false, false, 4, &logw),
        ),
        (
            "big_inv",
            &big,
            selector_params(10, 5, Integer, true, false, 4, &invw),
        ),
        (
            "small_class",
            &small,
            selector_params(-1, -1, Integer, true, false, 4, &logw),
        ),
        (
            "small_win",
            &small,
            selector_params(8, 4, Continuous, true, false, 4, &linw),
        ),
        (
            "small_noweights",
            &small,
            selector_params(3, 2, Integer, true, false, 4, &none),
        ),
    ];
    let selector = MRMFeatureSelectorScore::new(PartitionSolver);
    let mut out = Out::default();
    let mut selections: BTreeMap<&str, FeatureMap> = BTreeMap::new();
    for (name, fm, params) in runs {
        let selected = selector.select_mrm_feature(fm, &params).unwrap();
        dump_selected(&mut out, name, &selected);
        selections.insert(name, selected);
    }
    // Every run but `big_win` is compared record for record, including the
    // columns the solver picks among exactly tied optima (the last, as GLPK did).
    let exact = [
        "big_class.",
        "big_comp.",
        "big_inv.",
        "small_class.",
        "small_win.",
        "small_noweights.",
    ];
    let mut mine = Out::default();
    for (k, v) in out.records {
        if exact.iter().any(|p| k.starts_with(p)) {
            mine.records.insert(k, v);
        }
    }
    compare_with(
        "selector_score",
        load_oracle_prefixed("selector_score", &exact),
        &mine,
        &[],
    );

    // `big_win` (continuous variables, windows of 8 stepping 4) solves each
    // component group in two overlapping segments. Where candidates are exactly
    // tied, GLPK's pivoting picked a different tied column in the two segments
    // (both appear in its union), which a solver that looks at one row at a
    // time cannot reproduce. What must hold, and is checked from the data
    // alone: the groups are the same, every selected candidate of either side
    // has the least score of its group, and the selections differ only in
    // groups with tied least candidates.
    let oracle_records = load_oracle("selector_score");
    let mut oracle_groups: BTreeMap<String, Vec<i64>> = BTreeMap::new();
    let n: usize = oracle_records["big_win.n"].parse().unwrap();
    for i in 0..n {
        let group = oracle_records[&format!("big_win.f{i}.peptide_ref")].clone();
        let id: i64 = oracle_records[&format!("big_win.f{i}.id")].parse().unwrap();
        oracle_groups.entry(group).or_default().push(id);
    }
    let mut my_groups: BTreeMap<String, Vec<i64>> = BTreeMap::new();
    for f in &selections["big_win"].features {
        my_groups
            .entry(peptide_ref(f))
            .or_default()
            .push(f.unique_id as i64);
    }
    assert_eq!(
        oracle_groups.keys().collect::<Vec<_>>(),
        my_groups.keys().collect::<Vec<_>>()
    );
    let score_of = |id: i64| -> f64 {
        let f = big
            .features
            .iter()
            .find(|f| f.unique_id as i64 == id)
            .unwrap();
        compute_score(f, &linw).unwrap()
    };
    let mut extra = 0;
    for (group, theirs) in &oracle_groups {
        let mine = &my_groups[group];
        let candidates: Vec<f64> = big
            .features
            .iter()
            .filter(|f| &peptide_ref(f) == group)
            .map(|f| score_of(f.unique_id as i64))
            .collect();
        let least = candidates.iter().copied().fold(f64::INFINITY, f64::min);
        for id in theirs.iter().chain(mine.iter()) {
            assert_eq!(
                score_of(*id),
                least,
                "{group}: a selected candidate is not optimal"
            );
        }
        if theirs != mine {
            extra += theirs.len().abs_diff(mine.len());
            let tied = candidates.iter().filter(|c| **c == least).count();
            assert!(tied > 1, "{group}: selections differ without a tie");
        }
    }
    // the two groups with tied candidates (adn, adpglc)
    assert_eq!(extra, 2);
}

/// Forces every candidate that is alone in its group row to 1 and everything
/// else to 0: the optimum of the QMIP programs of the upstream small data,
/// where each component group has exactly one candidate peak group.
struct SoleCandidateSolver;
impl LpSolver for SoleCandidateSolver {
    fn solve(&self, problem: &LinearProgram) -> openms::Result<LpSolution> {
        let mut values = vec![0.0; problem.columns.len()];
        for row in &problem.rows {
            if row.indices.len() == 1 && row.lower == 1.0 && row.upper == 1.0 {
                values[row.indices[0]] = 1.0;
            }
        }
        Ok(LpSolution {
            status: SolverStatus::Optimal,
            column_values: values,
        })
    }
}

#[test]
fn selector_qmip_selection_matches_the_release_build() {
    use LambdaScore::{InverseLog, InverseLog10, Linear, Log};
    use VariableType::{Continuous, Integer};
    let small = load_small();
    let logw = weights(&[("sn_ratio", Log), ("peak_apices_sum", Log)]);
    let linw = weights(&[("sn_ratio", Linear), ("peak_apices_sum", Linear)]);
    let invw = weights(&[("sn_ratio", InverseLog), ("peak_apices_sum", InverseLog10)]);
    let runs: Vec<(&str, Vec<SelectorParameters>)> = vec![
        (
            "int_batch",
            vec![
                selector_params(8, 4, Integer, true, false, 4, &linw),
                selector_params(-1, -1, Integer, true, false, 4, &linw),
            ],
        ),
        (
            "cont_batch",
            vec![
                selector_params(8, 4, Continuous, true, false, 4, &linw),
                selector_params(-1, -1, Continuous, true, false, 4, &linw),
            ],
        ),
        (
            "int_locality",
            vec![selector_params(4, 2, Integer, true, true, 2, &linw)],
        ),
        (
            "int_log",
            vec![selector_params(-1, -1, Integer, true, false, 4, &logw)],
        ),
        (
            "cont_inv",
            vec![selector_params(4, 2, Continuous, true, false, 2, &invw)],
        ),
    ];
    let mut out = Out::default();
    for (name, params) in runs {
        let selected = batch_mrm_features_qmip(&small, &params, SoleCandidateSolver).unwrap();
        dump_selected(&mut out, name, &selected);
    }
    // The selection is the one the Release build's GLPK finds: every group has
    // one candidate. The program the port builds is checked against the source
    // in tests/mrm_feature_selector.rs; this checks which features it keeps.
    compare("selector_qmip", &out, &[]);
}

#[test]
fn selector_with_duplicate_row_names_aborts_the_release_build() {
    // select_transition_group=false hands GLPK two rows named
    // "Serotonin_constraint" and the Release build aborts (SIGABRT, exit 134)
    // after "** duplicate name Serotonin_constraint". The port has no such
    // abort: the program it builds has a fixed row without a column (the
    // second "Serotonin" row finds its candidate already created), which no
    // solution satisfies, and the selection is refused as infeasible.
    let exit = std::fs::read_to_string(oracle("selector_dupname.exit")).unwrap();
    assert_eq!(exit.trim(), "134");
    let stdout = std::fs::read_to_string(oracle("selector_dupname.stdout")).unwrap();
    assert!(!stdout.contains("small_comp.n"));
    let small = load_small();
    let selector = MRMFeatureSelectorScore::new(PartitionSolver);
    let params = selector_params(
        -1,
        -1,
        VariableType::Integer,
        false,
        false,
        4,
        &weights(&[
            ("sn_ratio", LambdaScore::Linear),
            ("peak_apices_sum", LambdaScore::Linear),
        ]),
    );
    let err = selector.select_mrm_feature(&small, &params).unwrap_err();
    assert!(matches!(err, openms::Error::InvalidValue(m) if m.contains("Infeasible")));
}

fn peptide(id: &str, rt: f64) -> Peptide {
    let mut p = Peptide::default();
    p.id = id.into();
    let mut r = RetentionTime::default();
    r.set_rt(rt);
    r.retention_time_type = RTType::Normalized;
    p.rts.push(r);
    p
}

fn mk_transition(id: &str, pep: &str, mz: f64, lib: f64) -> ReactionMonitoringTransition {
    let mut t = ReactionMonitoringTransition::new();
    t.set_native_id(id);
    t.set_precursor_mz(400.0);
    t.set_product_mz(mz);
    t.set_peptide_ref(pep);
    t.set_library_intensity(lib);
    t
}

#[test]
fn confidence_matches_the_release_build() {
    let mut out = Out::default();
    {
        let mut scoring = ConfidenceScoring::new(false);
        let mut library = TargetedExperiment::new();
        library.add_peptide(peptide("pep_1", 50.0));
        library.add_transition(mk_transition("tr_10", "pep_1", 500.0, 500.0));
        library.add_peptide(peptide("pep_2", 60.0));
        library.add_transition(mk_transition("tr_20", "pep_2", 500.0, 500.0));
        let mut features = FeatureMap::default();
        for i in 0..2 {
            let mut f = Feature::new(60.0, 0.0, 0.0);
            f.metadata.insert(
                "PeptideRef".into(),
                if i == 0 { "pep_1" } else { "pep_2" }.into(),
            );
            f.quality = -1.0;
            let mut sub = Feature::new(0.0, 500.0, 1.0);
            sub.metadata.insert(
                "native_id".into(),
                if i == 0 { "tr_10" } else { "tr_20" }.into(),
            );
            f.subordinates.push(sub);
            features.features.push(f);
        }
        scoring.initialize(library, 0, 0, TransformationDescription::default());
        scoring.initialize_glm(0.0, -1.0, -1.0);
        scoring.score_map(&mut features).unwrap();
        for (i, f) in features.features.iter().enumerate() {
            out.p(&format!("class.f{i}.quality"), f64::from(f.quality));
            out.p(
                &format!("class.f{i}.glm"),
                f.metadata["GLM_score"].as_f64().unwrap(),
            );
            out.p(
                &format!("class.f{i}.fdr"),
                f.metadata["local_FDR"].as_f64().unwrap(),
            );
        }
    }
    let mut library = TargetedExperiment::new();
    for i in 0..7 {
        let pid = format!("p{i}");
        library.add_peptide(peptide(&pid, 10.0 * f64::from(i)));
        let nt = if i == 3 { 2 } else { 3 };
        for j in 0..nt {
            library.add_transition(mk_transition(
                &format!("t{i}_{j}"),
                &pid,
                400.0 + 50.0 * f64::from(j),
                100.0 * f64::from(j + 1) + f64::from(i),
            ));
        }
    }
    for (name, n_decoys, n_trans) in [
        ("d0t0", 0, 0),
        ("d3t0", 3, 0),
        ("d99t0", 99, 0),
        ("d3t2", 3, 2),
        ("d0t2", 0, 2),
        ("d1t0", 1, 0),
    ] {
        let mut features = FeatureMap::default();
        for i in 0..7 {
            if i == 3 {
                continue;
            }
            let mut f = Feature::new(10.0 * f64::from(i) + 1.5, 0.0, 0.0);
            f.metadata
                .insert("PeptideRef".into(), format!("p{i}").into());
            f.quality = -1.0;
            for j in 0..3 {
                let mut sub =
                    Feature::new(0.0, 400.0 + 50.0 * f64::from(j), ((j + 1) * 7 + i) as f32);
                sub.metadata
                    .insert("native_id".into(), format!("t{i}_{j}").into());
                f.subordinates.push(sub);
            }
            features.features.push(f);
        }
        let mut scoring = ConfidenceScoring::new(false);
        scoring.initialize(
            library.clone(),
            n_decoys,
            n_trans,
            TransformationDescription::default(),
        );
        scoring.initialize_glm(0.5, -0.01, -1.0);
        scoring.score_map(&mut features).unwrap();
        for (i, f) in features.features.iter().enumerate() {
            let p = format!("{name}.f{i}");
            out.p(&format!("{p}.quality"), f64::from(f.quality));
            out.p(
                &format!("{p}.glm"),
                f.metadata["GLM_score"].as_f64().unwrap(),
            );
            out.p(
                &format!("{p}.fdr"),
                f.metadata["local_FDR"].as_f64().unwrap(),
            );
        }
    }
    compare("confidence", &out, &[]);
}

fn generic_transitions(identification: bool) -> TargetedExperiment {
    let mut exp = TargetedExperiment::new();
    exp.add_protein(Protein::default());
    for (id, rt, seq) in [("tr_gr1", 0.44, "PEPTIDEA"), ("tr_gr2", 0.2, "PEPTIDEB")] {
        let mut p = peptide(id, rt);
        p.sequence = seq.into();
        p.protein_refs.push(String::new());
        exp.add_peptide(p);
    }
    let on = Some((true, false, true));
    type Row<'a> = (&'a str, &'a str, f64, f64, f64, Option<(bool, bool, bool)>);
    let rows: Vec<Row> = vec![
        ("tr1", "tr_gr1", 500.0, 628.435, 1.0, on),
        ("tr2", "tr_gr1", 500.0, 654.38, 2.0, on),
        ("tr3", "tr_gr2", 501.0, 618.31, 10000.0, on),
        ("tr4", "tr_gr2", 501.0, 628.435, 2000.0, on),
        (
            "tr5",
            "tr_gr2",
            501.0,
            651.3,
            4300.0,
            if identification {
                Some((true, true, true))
            } else {
                on
            },
        ),
    ];
    for (id, group, pmz, mz, lib, flags) in rows {
        let mut t = mk_transition(id, group, mz, lib);
        t.set_precursor_mz(pmz);
        if identification {
            if let Some((d, i, q)) = flags {
                t.set_detecting_transition(d);
                t.set_identifying_transition(i);
                t.set_quantifying_transition(q);
            }
        }
        exp.add_transition(t);
    }
    if identification {
        let mut t = mk_transition("tr2", "tr_gr2", 654.38, 2.0);
        t.set_precursor_mz(500.0);
        t.set_detecting_transition(false);
        t.set_identifying_transition(true);
        t.set_quantifying_transition(true);
        exp.add_transition(t);
    }
    exp
}

fn run_finder(out: &mut Out, name: &str, identification: bool, param: &Param) {
    let mut ff = MRMFeatureFinderScoring::from_param(param).unwrap();
    ff.set_unique_id_generator(UniqueIdGenerator::from_seed(42));
    let bytes = std::fs::read(data("OpenSwath_generic_input.mzML")).unwrap();
    let experiment = Arc::new(openms::format::mzml::read(bytes.as_slice()).unwrap());
    let light = convert_targeted_exp(&generic_transitions(identification)).unwrap();
    let chromatogram_ptr =
        SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(experiment).unwrap();
    let swath_ptr =
        SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(Arc::new(MSExperiment::default()))
            .unwrap();
    let swath_maps = [SwathMap {
        sptr: Some(swath_ptr),
        ..SwathMap::default()
    }];
    let mut feature_file = FeatureMap::default();
    let mut group_map = TransitionGroupMapType::new();
    ff.pick_experiment_light(
        &chromatogram_ptr,
        &mut feature_file,
        &light,
        &TransformationDescription::default(),
        &swath_maps,
        &mut group_map,
    )
    .unwrap();
    out.pi(&format!("{name}.ngroups"), group_map.len() as i64);
    for (key, group) in &group_map {
        out.pi(&format!("{name}.group.{key}.size"), group.size() as i64);
        out.pi(
            &format!("{name}.group.{key}.nfeatures"),
            group.features().len() as i64,
        );
    }
    out.pi(&format!("{name}.nout"), feature_file.len() as i64);
    for (i, f) in feature_file.features.iter().enumerate() {
        let p = format!("{name}.out{i}");
        dump_feature(out, &p, f, true);
        dump_ids(out, &p, f);
    }
    out.pi(
        &format!("{name}.nprot"),
        feature_file.protein_identifications.len() as i64,
    );
    for (i, pi) in feature_file.protein_identifications.iter().enumerate() {
        out.ps(&format!("{name}.prot{i}.identifier"), &pi.identifier);
        out.pi(&format!("{name}.prot{i}.nhits"), pi.hits.len() as i64);
    }
}

fn set(param: &mut Param, key: &str, value: ParamValue) {
    param.set_value(key, value, "", &[]).unwrap();
}

#[test]
fn finder_matches_the_release_build() {
    let defaults = MRMFeatureFinderScoring::defaults().unwrap();
    let mut out = Out::default();
    run_finder(&mut out, "generic_default", false, &defaults);
    let mut legacy = defaults.clone();
    set(
        &mut legacy,
        "TransitionGroupPicker:PeakPickerChromatogram:method",
        ParamValue::String("legacy".into()),
    );
    set(
        &mut legacy,
        "TransitionGroupPicker:PeakPickerChromatogram:peak_width",
        ParamValue::Float(40.0),
    );
    run_finder(&mut out, "generic_legacy", false, &legacy);
    let mut ident = legacy.clone();
    set(
        &mut ident,
        "Scores:use_uis_scores",
        ParamValue::String("true".into()),
    );
    run_finder(&mut out, "ident_legacy", true, &ident);
    let mut stop = defaults.clone();
    set(
        &mut stop,
        "stop_report_after_feature",
        ParamValue::Integer(1),
    );
    run_finder(&mut out, "generic_stop1", false, &stop);
    let mut mi = defaults.clone();
    set(
        &mut mi,
        "Scores:use_total_mi_score",
        ParamValue::String("true".into()),
    );
    set(
        &mut mi,
        "Scores:use_mi_score",
        ParamValue::String("true".into()),
    );
    set(
        &mut mi,
        "scoring_model",
        ParamValue::String("single_transition".into()),
    );
    run_finder(&mut out, "generic_mi", false, &mi);
    compare("finder", &out, &[]);
}
