// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `MapAlignmentAlgorithmKD` and `MapAlignmentAlgorithmTreeGuided` against
//! the executed C++ Release build (evidence tier 1).
//!
//! `tests/data/map_alignment_kd_tree_oracle/` holds the inputs
//! (`kd_*.tsv`, `tg_*.tsv`, written by the oracle's `gen_cases.py`) and what
//! the driver recorded from `openms4-release-bc9cc12-c19e494-174b576` on
//! ibminode06 (Linux x86_64, conda-forge gcc 14.4.0, two fresh-process runs,
//! byte-identical): `records.tsv` (every number at `%.17g`), `notes.tsv`
//! (model types and data-point notes), `stderr.txt` (the LOWESS warnings),
//! and `kd_defaults.ini` / `tg_defaults.ini` (the defaults of
//! `FeatureGroupingAlgorithmKD` and `MapAlignmentAlgorithmTreeGuided`).
//! Driver, generator and hashes: `tests/data/map_alignment_provenance.json`.
//!
//! The KD cases run on the real `KDTreeFeatureMaps` in C++ and on a
//! brute-force stand-in here: neighbourhood sizes, component labels, the
//! filtered components and the transformed retention times must agree
//! exactly. The tree-guided cases compare the similarity of every pair of
//! maps (as the distance of a two-map clustering), the tree's first merge,
//! the percentile lists, and, following the C++ tree, the merged map, the
//! merge order, the transformations' data points and the final retention
//! times, all bit for bit.

#![cfg(all(feature = "paramxml", feature = "featurexml"))]

#[path = "support/kd_pool.rs"]
mod kd_pool;

use kd_pool::{Pool, PoolFeature};
use openms::analysis::mapmatching::alignment_kd::{KdFeatureMaps, MapAlignmentAlgorithmKD};
use openms::analysis::mapmatching::alignment_tree_guided::{
    MapAlignmentAlgorithmTreeGuided, peptide_identifications_pearson_similarity,
};
use openms::concept::unique_id::UniqueIdGenerator;
use openms::format::{featurexml, paramxml};
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::kernel::{BaseFeature, Feature, FeatureMap};
use openms::param::{Param, ParamValue};
use std::collections::HashMap;

const DIR: &str = "tests/data/map_alignment_kd_tree_oracle";

/// Records per case, keyed by record name.
struct Records {
    numbers: HashMap<String, HashMap<String, f64>>,
    texts: HashMap<String, HashMap<String, String>>,
}

impl Records {
    fn load() -> Self {
        let mut numbers: HashMap<String, HashMap<String, f64>> = HashMap::new();
        for line in std::fs::read_to_string(format!("{DIR}/records.tsv"))
            .unwrap()
            .lines()
        {
            let mut f = line.split('\t');
            let (case, key, value) = (f.next().unwrap(), f.next().unwrap(), f.next().unwrap());
            numbers
                .entry(case.into())
                .or_default()
                .insert(key.into(), value.parse().unwrap());
        }
        let mut texts: HashMap<String, HashMap<String, String>> = HashMap::new();
        for line in std::fs::read_to_string(format!("{DIR}/notes.tsv"))
            .unwrap()
            .lines()
        {
            let mut f = line.splitn(3, '\t');
            let (case, key, value) = (f.next().unwrap(), f.next().unwrap(), f.next().unwrap());
            texts
                .entry(case.into())
                .or_default()
                .insert(key.into(), value.into());
        }
        Self { numbers, texts }
    }
    fn get(&self, case: &str, key: &str) -> f64 {
        *self.numbers[case]
            .get(key)
            .unwrap_or_else(|| panic!("{case}: no record {key}"))
    }
    fn has(&self, case: &str, key: &str) -> bool {
        self.numbers[case].contains_key(key)
    }
    fn count(&self, case: &str, key: &str) -> usize {
        self.get(case, key) as usize
    }
    fn text(&self, case: &str, key: &str) -> &str {
        &self.texts[case][key]
    }
}

/// Bit-for-bit equality of a recorded double (`%.17g` round-trips) on Linux
/// x86_64, the oracle's platform; elsewhere `1e-12` relative, for a platform
/// libm that may round differently, as in `tests/map_alignment_oracle.rs`.
fn same(actual: f64, expected: f64, what: &str) {
    if cfg!(all(target_os = "linux", target_arch = "x86_64")) {
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "{what}: {actual} != {expected}"
        );
    } else {
        let scale = actual.abs().max(expected.abs()).max(f64::MIN_POSITIVE);
        assert!(
            (actual - expected).abs() / scale <= 1e-12,
            "{what}: {actual} != {expected}"
        );
    }
}

// ------------------------------------------------------------------ KD --

fn read_kd(name: &str) -> Vec<PoolFeature> {
    std::fs::read_to_string(format!("{DIR}/kd_{name}.tsv"))
        .unwrap()
        .lines()
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            PoolFeature {
                map: f[0].parse().unwrap(),
                rt: f[1].parse().unwrap(),
                mz: f[2].parse().unwrap(),
                intensity: f[3].parse::<i64>().unwrap() as f32,
                charge: f[4].parse().unwrap(),
            }
        })
        .collect()
}

/// `FeatureGroupingAlgorithmKD`'s defaults as the C++ build stored them, with
/// `LOWESS:span` restored to the in-memory `2.0 / 3.0` (the paramXML text
/// keeps 15 significant digits), and the case's overrides.
fn kd_param(overrides: &[(&str, ParamValue)]) -> Param {
    let mut p = paramxml::load(format!("{DIR}/kd_defaults.ini")).unwrap();
    p.set_value("LOWESS:span", ParamValue::Float(2.0 / 3.0), "", &[])
        .unwrap();
    for (key, value) in overrides {
        p.set_value(key, value.clone(), "", &[]).unwrap();
    }
    p
}

fn kd_structure(
    records: &Records,
    case: &str,
    prefix: &str,
    algorithm: &MapAlignmentAlgorithmKD,
    pool: &Pool,
    p: &Param,
) {
    let float = |k: &str| p.value(k).unwrap().to_f64().unwrap();
    let ppm = p.value("mz_unit").unwrap().as_str().unwrap() == "ppm";
    assert_eq!(pool.len(), records.count(case, &format!("{prefix}size")));
    for i in 0..pool.len() {
        let nb = pool
            .neighborhood(
                i,
                float("warp:rt_tol"),
                float("warp:mz_tol"),
                ppm,
                false,
                float("warp:max_pairwise_log_fc"),
            )
            .unwrap();
        assert_eq!(
            nb.len(),
            records.count(case, &format!("{prefix}nb.{i}")),
            "{case} nb {i}"
        );
    }
    let (labels, count) = algorithm.compute_ccs(pool).unwrap();
    assert_eq!(
        count,
        records.count(case, &format!("{prefix}ncc")),
        "{case}"
    );
    for (i, label) in labels.iter().enumerate() {
        assert_eq!(
            *label,
            records.count(case, &format!("{prefix}cc.{i}")),
            "{case} cc {i}"
        );
    }
    let ccs = algorithm.get_ccs(pool).unwrap();
    let filtered = algorithm.filter_ccs(pool, &ccs).unwrap();
    assert_eq!(
        filtered.len(),
        records.count(case, &format!("{prefix}nfiltered")),
        "{case}"
    );
    for (label, cc) in &filtered {
        assert_eq!(
            cc.len(),
            records.count(case, &format!("{prefix}filtered.{label}")),
            "{case} filtered {label}"
        );
    }
}

/// Runs one KD case, returning the LOWESS warnings.
fn kd_case(
    records: &Records,
    case: &str,
    data: &str,
    overrides: &[(&str, ParamValue)],
    partition: bool,
) -> Vec<String> {
    let p = kd_param(overrides);
    let features = read_kd(data);
    let num_maps = features.iter().map(|f| f.map).max().unwrap() + 1;
    let mut algorithm = MapAlignmentAlgorithmKD::new(num_maps, &p).unwrap();
    if !partition {
        let mut pool = Pool::new(features);
        kd_structure(records, case, "", &algorithm, &pool, &p);
        algorithm.add_rt_fit_data(&pool).unwrap();
        let warnings = algorithm.fit_lowess().unwrap();
        algorithm.transform(&mut pool).unwrap();
        for i in 0..pool.len() {
            same(
                pool.rt(i),
                records.get(case, &format!("rt.{i}")),
                &format!("{case} rt {i}"),
            );
        }
        warnings
    } else {
        // The driver's two m/z partitions, each in map order.
        let part = |low: bool| {
            let mut sorted: Vec<PoolFeature> = Vec::new();
            for m in 0..num_maps {
                sorted.extend(
                    features
                        .iter()
                        .filter(|f| f.map == m && (f.mz < 800.0) == low)
                        .copied(),
                );
            }
            Pool::new(sorted)
        };
        let (mut low, mut high) = (part(true), part(false));
        kd_structure(records, case, "low.", &algorithm, &low, &p);
        kd_structure(records, case, "high.", &algorithm, &high, &p);
        algorithm.add_rt_fit_data(&low).unwrap();
        algorithm.add_rt_fit_data(&high).unwrap();
        let warnings = algorithm.fit_lowess().unwrap();
        algorithm.transform(&mut low).unwrap();
        algorithm.transform(&mut high).unwrap();
        for (name, pool) in [("low", &low), ("high", &high)] {
            for i in 0..pool.len() {
                same(
                    pool.rt(i),
                    records.get(case, &format!("{name}.rt.{i}")),
                    &format!("{case} {name} rt {i}"),
                );
            }
        }
        warnings
    }
}

#[test]
fn kd_alignment_matches_the_release_build() {
    let records = Records::load();
    let mut warnings = Vec::new();
    warnings.extend(kd_case(&records, "kd_base", "base", &[], false));
    warnings.extend(kd_case(&records, "kd_conflicts0", "conflicts", &[], false));
    warnings.extend(kd_case(
        &records,
        "kd_conflicts1",
        "conflicts",
        &[("warp:max_nr_conflicts", ParamValue::Integer(1))],
        false,
    ));
    warnings.extend(kd_case(
        &records,
        "kd_conflicts_any",
        "conflicts",
        &[
            ("warp:max_nr_conflicts", ParamValue::Integer(-1)),
            ("warp:min_rel_cc_size", ParamValue::Float(0.25)),
        ],
        false,
    ));
    warnings.extend(kd_case(
        &records,
        "kd_abs",
        "abs",
        &[
            ("mz_unit", ParamValue::String("Da".into())),
            ("warp:mz_tol", ParamValue::Float(0.01)),
            ("warp:max_pairwise_log_fc", ParamValue::Float(-1.0)),
            ("warp:rt_tol", ParamValue::Float(40.0)),
        ],
        false,
    ));
    warnings.extend(kd_case(&records, "kd_few", "few", &[], false));
    warnings.extend(kd_case(
        &records,
        "kd_span",
        "base",
        &[("LOWESS:span", ParamValue::Float(0.3))],
        false,
    ));
    warnings.extend(kd_case(&records, "kd_partition", "base", &[], true));
    let expected: Vec<String> = std::fs::read_to_string(format!("{DIR}/stderr.txt"))
        .unwrap()
        .lines()
        .filter(|l| l.starts_with("Warning: Only"))
        .map(str::to_owned)
        .collect();
    assert_eq!(warnings, expected);
}

// --------------------------------------------------------- tree guided --

fn read_tg(name: &str, maps: usize) -> Vec<FeatureMap> {
    (0..maps)
        .map(|m| {
            let mut features: Vec<Feature> = Vec::new();
            for line in std::fs::read_to_string(format!("{DIR}/tg_{name}_{m}.tsv"))
                .unwrap()
                .lines()
            {
                let f: Vec<&str> = line.split('\t').collect();
                if f[0] == "F" {
                    let mut feature = Feature::from(BaseFeature::new(
                        f[2].parse().unwrap(),
                        f[3].parse().unwrap(),
                        f[4].parse::<i64>().unwrap() as f32,
                    ));
                    feature.base.unique_id = f[1].parse().unwrap();
                    feature.base.charge = f[5].parse().unwrap();
                    features.push(feature);
                } else {
                    let hits = if f[4] == "-" {
                        Vec::new()
                    } else {
                        vec![PeptideHit {
                            score: f[3].parse().unwrap(),
                            rank: 1,
                            charge: 2,
                            sequence: f[4].parse().unwrap(),
                            ..PeptideHit::default()
                        }]
                    };
                    let id = PeptideIdentification {
                        rt: Some(f[1].parse().unwrap()),
                        mz: Some(f[2].parse().unwrap()),
                        score_type: "q-value".into(),
                        higher_score_better: false,
                        hits,
                        ..PeptideIdentification::default()
                    };
                    features
                        .last_mut()
                        .unwrap()
                        .base
                        .peptide_identifications
                        .push(id);
                }
            }
            FeatureMap::from_features(features)
        })
        .collect()
}

fn class_test_maps() -> Vec<FeatureMap> {
    (0..3)
        .map(|i| {
            featurexml::load(format!(
                "tests/data/map_alignment_tree_guided_in{i}.featureXML"
            ))
            .unwrap()
        })
        .collect()
}

fn tg_case(records: &Records, case: &str, maps: Vec<FeatureMap>, model_type: &str) {
    let n = maps.len();
    // Pairwise similarities: the distance of a two-map clustering is
    // `1 - similarity` in float, or a dummy node of -1 at distance >= 1.
    let (lists, ranges) = MapAlignmentAlgorithmTreeGuided::extract_seq_and_rt(&maps).unwrap();
    for i in 0..n {
        for j in i + 1..n {
            let key = format!("pair.{i}.{j}");
            let similarity = peptide_identifications_pearson_similarity(&lists[j], &lists[i]);
            if records.has(case, &format!("{key}.error")) {
                assert!(similarity.is_err(), "{case} {key}");
                continue;
            }
            let distance = 1.0f32 - similarity.unwrap();
            let expected = records.get(case, &key);
            if expected == -1.0 {
                assert!(distance >= 1.0, "{case} {key}");
            } else {
                same(f64::from(distance), expected, &format!("{case} {key}"));
            }
        }
    }
    // The tree of the C++ clustering; its first merge is the closest pair.
    let nodes = records.count(case, "nodes");
    let tree: Vec<(usize, usize)> = (0..nodes)
        .map(|k| {
            (
                records.count(case, &format!("node.{k}.left")),
                records.count(case, &format!("node.{k}.right")),
            )
        })
        .collect();
    // `buildTree` itself (ClusterHierarchical with AverageLinkage on the
    // Pearson similarity): every merge and distance of the Release build's
    // tree, and the sorted RT lists it returns.
    let (built, built_ranges) = MapAlignmentAlgorithmTreeGuided::build_tree(&maps).unwrap();
    assert_eq!(built.len(), nodes, "{case} tree size");
    for (k, node) in built.iter().enumerate() {
        assert_eq!(
            (node.left_child, node.right_child),
            tree[k],
            "{case} node {k}"
        );
        same(
            f64::from(node.distance),
            records.get(case, &format!("node.{k}.distance")),
            &format!("{case} node {k} distance"),
        );
    }
    assert_eq!(built_ranges, ranges, "{case} ranges");
    let (l, r) = tree[0];
    same(
        records.get(case, "node.0.distance"),
        records.get(case, &format!("pair.{}.{}", l.min(r), l.max(r))),
        &format!("{case} first merge"),
    );
    for (i, range) in ranges.iter().enumerate() {
        assert_eq!(range.len(), records.count(case, &format!("range.{i}.size")));
        for (k, value) in range.iter().enumerate() {
            same(*value, records.get(case, &format!("range.{i}.{k}")), case);
        }
    }

    let mut p = Param::new();
    p.set_value("model_type", ParamValue::String(model_type.into()), "", &[])
        .unwrap();
    let mut aligner = MapAlignmentAlgorithmTreeGuided::from_param(&p).unwrap();
    let mut generator = UniqueIdGenerator::from_seed(1);
    let mut copied = maps.clone();
    let merged = aligner
        .tree_guided_alignment(&tree, &mut copied, &ranges, &mut generator)
        .unwrap();
    assert_eq!(merged.trafo_order.len(), records.count(case, "order.size"));
    for (k, &m) in merged.trafo_order.iter().enumerate() {
        assert_eq!(
            m,
            records.count(case, &format!("order.{k}")),
            "{case} order {k}"
        );
    }
    let mt = &merged.map_transformed;
    assert_eq!(mt.len(), records.count(case, "mt.size"));
    for (k, feature) in mt.features.iter().enumerate() {
        same(
            feature.rt,
            records.get(case, &format!("mt.{k}.rt")),
            &format!("{case} mt {k}"),
        );
        let original = feature.metadata.get("original_RT");
        assert_eq!(
            original.is_some(),
            records.get(case, &format!("mt.{k}.orig")) == 1.0,
            "{case} mt {k}"
        );
        if let Some(original) = original {
            same(
                original.as_f64().unwrap(),
                records.get(case, &format!("mt.{k}.orig_rt")),
                &format!("{case} mt {k} original_RT"),
            );
        }
    }

    let result = aligner.compute_trafos_by_original_rt(&maps, mt).unwrap();
    assert_eq!(result.transformations.len(), n);
    for (i, trafo) in result.transformations.iter().enumerate() {
        assert_eq!(
            trafo.model_type(),
            records.text(case, &format!("trafo.{i}.model"))
        );
        let points = trafo.data_points();
        assert_eq!(
            points.len(),
            records.count(case, &format!("trafo.{i}.n")),
            "{case} {i}"
        );
        for (k, point) in points.iter().enumerate() {
            same(
                point.x,
                records.get(case, &format!("trafo.{i}.{k}.x")),
                case,
            );
            same(
                point.y,
                records.get(case, &format!("trafo.{i}.{k}.y")),
                case,
            );
            assert_eq!(
                point.note,
                records.text(case, &format!("trafo.{i}.{k}.note"))
            );
        }
    }

    // The whole `align`: the tree is built, not given, and the final
    // retention times are the Release build's.
    let mut aligned = maps.clone();
    let end_to_end = aligner
        .align(&mut aligned, &mut UniqueIdGenerator::from_seed(1))
        .unwrap();
    assert!(
        end_to_end.info[0].starts_with("  Alignment follows Newick tree: ( "),
        "{case}: {:?}",
        end_to_end.info[0]
    );
    assert_eq!(end_to_end.transformations.len(), n);
    for (trafo, expected) in end_to_end
        .transformations
        .iter()
        .zip(&result.transformations)
    {
        assert_eq!(trafo.data_points().len(), expected.data_points().len());
        for (point, expected) in trafo.data_points().iter().zip(expected.data_points()) {
            same(point.x, expected.x, case);
            same(point.y, expected.y, case);
        }
    }
    for (i, map) in aligned.iter().enumerate() {
        for (k, feature) in map.features.iter().enumerate() {
            same(
                feature.rt,
                records.get(case, &format!("final.{i}.{k}")),
                &format!("{case} align final {i} {k}"),
            );
        }
    }

    let mut transformed = maps;
    MapAlignmentAlgorithmTreeGuided::compute_transformed_feature_maps(
        &mut transformed,
        &result.transformations,
    )
    .unwrap();
    for (i, map) in transformed.iter().enumerate() {
        for (k, feature) in map.features.iter().enumerate() {
            same(
                feature.rt,
                records.get(case, &format!("final.{i}.{k}")),
                &format!("{case} final {i} {k}"),
            );
        }
    }
}

#[test]
fn tree_guided_alignment_matches_the_release_build() {
    let records = Records::load();
    tg_case(&records, "tg_classtest", class_test_maps(), "b_spline");
    tg_case(&records, "tg_classtest_linear", class_test_maps(), "linear");
    tg_case(&records, "tg_four", read_tg("four", 4), "b_spline");
    tg_case(&records, "tg_five", read_tg("five", 5), "b_spline");
    tg_case(&records, "tg_five_linear", read_tg("five", 5), "linear");
    tg_case(&records, "tg_five_lowess", read_tg("five", 5), "lowess");
    tg_case(
        &records,
        "tg_five_interpolated",
        read_tg("five", 5),
        "interpolated",
    );
}

#[test]
fn tree_guided_defaults_match_the_release_build() {
    // The C++ build's defaults, read back from its ParamXMLFile::store output,
    // and the port's, both written by the port's writer, agree byte for byte:
    // names, values, types, descriptions, restrictions, tags and order. The
    // one exception is the LOWESS span, 2/3 in memory on both sides, which
    // the C++ text keeps to 15 significant digits.
    let write = |p: &Param| {
        let mut out = Vec::new();
        paramxml::write(&mut out, p).unwrap();
        String::from_utf8(out).unwrap()
    };
    let stored = write(&paramxml::load(format!("{DIR}/tg_defaults.ini")).unwrap());
    let ours = write(&MapAlignmentAlgorithmTreeGuided::defaults().unwrap());
    let span = "<ITEM name=\"span\" value=\"0.6666666666666666\"";
    assert_eq!(ours.matches(span).count(), 1);
    assert_eq!(
        ours.replace(span, "<ITEM name=\"span\" value=\"0.666666666666667\""),
        stored
    );
}
