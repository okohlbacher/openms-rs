// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Executed differential (evidence tier 1) for the clustering ports against
//! the Linux x86_64 Release build of OpenMS at bc9cc12-c19e494-174b576.
//!
//! `tests/data/ml_clustering_oracle.tsv` is the unmodified `results/oracle.tsv`
//! of the FF2a oracle (`../oracle/FF2a-clustering/`: `driver.cpp`, run twice on
//! ibminode06 with byte-identical results; sha256 recorded in
//! `tests/data/ml_clustering_provenance.json`). The driver generates every
//! input with a 64-bit LCG; this test regenerates the same inputs, runs the
//! ports and prints the same records, and compares them to the oracle
//! byte for byte. Floating-point values are compared as bit patterns.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use openms::analysis::mapmatching::grid_feature::GridFeature;
use openms::analysis::mapmatching::qt_cluster::{BulkData, Element, QTCluster};
use openms::chemistry::AASequence;
use openms::data_structures::binary_tree_node::BinaryTreeNode;
use openms::data_structures::distance_matrix::DistanceMatrix;
use openms::data_structures::dposition::DPosition;
use openms::format::file_info::text_format::to_str_f32;
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::kernel::features::BaseFeature;
use openms::math::ml::clustering::average_linkage::AverageLinkage;
use openms::math::ml::clustering::cluster_analyzer::ClusterAnalyzer;
use openms::math::ml::clustering::cluster_functor::ClusterFunctor;
use openms::math::ml::clustering::complete_linkage::CompleteLinkage;
use openms::math::ml::clustering::grid_based_cluster::GridBasedCluster;
use openms::math::ml::clustering::grid_based_clustering::GridBasedClustering;
use openms::math::ml::clustering::hash_grid::HashGrid;
use openms::math::ml::clustering::single_linkage::SingleLinkage;

const ORACLE: &str = include_str!("data/ml_clustering_oracle.tsv");

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    }
    fn below(&mut self, n: u32) -> u32 {
        ((self.next() >> 33) % u64::from(n)) as u32
    }
}

fn fb(f: f32) -> String {
    format!("{:08x}", f.to_bits())
}

fn db(d: f64) -> String {
    format!("{:016x}", d.to_bits())
}

fn join<T: ToString>(values: &[T]) -> String {
    if values.is_empty() {
        return "-".into();
    }
    values
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

/// The oracle's records, grouped by case, with each case's exit code.
fn oracle_cases() -> BTreeMap<String, (Vec<String>, String)> {
    let mut cases = BTreeMap::new();
    let mut pending = Vec::new();
    for line in ORACLE.lines() {
        if let Some(rest) = line.strip_prefix("EXIT\t") {
            let (name, code) = rest.split_once('\t').unwrap();
            cases.insert(
                name.to_owned(),
                (std::mem::take(&mut pending), code.to_owned()),
            );
        } else {
            pending.push(line.to_owned());
        }
    }
    assert!(pending.is_empty(), "records after the last EXIT line");
    cases
}

fn compare(case: &str, produced: &str) {
    let cases = oracle_cases();
    let (expected, code) = cases
        .get(case)
        .unwrap_or_else(|| panic!("case {case} missing from the oracle"));
    assert_eq!(code, "0", "oracle case {case} did not exit cleanly");
    let produced: Vec<&str> = produced.lines().collect();
    for (i, (p, e)) in produced.iter().zip(expected).enumerate() {
        assert_eq!(p, e, "case {case}, record {i}");
    }
    assert_eq!(produced.len(), expected.len(), "case {case}: record count");
}

// ---------------------------------------------------------------- linkage

fn make_matrix(n: usize, seed: u64, levels: u32) -> DistanceMatrix<f32> {
    let mut r = Lcg(seed);
    let mut m = DistanceMatrix::with_value(n, 1.0_f32).unwrap();
    for i in 1..n {
        for j in 0..i {
            let v = if levels == 0 {
                r.unit() as f32
            } else {
                r.below(levels) as f32 / levels as f32
            };
            m.set_quick(i, j, v).unwrap();
        }
    }
    m.update_min_element().unwrap();
    m
}

fn print_tree(out: &mut String, tag: &str, tree: &[BinaryTreeNode]) {
    for (i, node) in tree.iter().enumerate() {
        writeln!(
            out,
            "{tag}\t{i}\t{}\t{}\t{}",
            node.left_child,
            node.right_child,
            fb(node.distance)
        )
        .unwrap();
    }
}

fn analyze(out: &mut String, tag: &str, tree: &[BinaryTreeNode], original: &DistanceMatrix<f32>) {
    let ca = ClusterAnalyzer::new();
    let n = original.dimension_size();
    for (i, v) in ca
        .average_silhouette_width(tree, original)
        .unwrap()
        .iter()
        .enumerate()
    {
        writeln!(out, "ASW\t{tag}\t{i}\t{}", fb(*v)).unwrap();
    }
    for single in [false, true] {
        for (i, v) in ca
            .dunn_indices(tree, original, single)
            .unwrap()
            .iter()
            .enumerate()
        {
            writeln!(out, "DUNN{}\t{tag}\t{i}\t{}", u8::from(single), fb(*v)).unwrap();
        }
    }
    let mut ks = vec![1, 2, 3, n / 2, n - 1, n];
    ks.sort_unstable();
    ks.dedup();
    for k in ks {
        if k == 0 || k > n {
            continue;
        }
        let clusters = ca.cut(k, tree).unwrap();
        for (c, members) in clusters.iter().enumerate() {
            writeln!(out, "CUT\t{tag}\t{k}\t{c}\t{}", join(members)).unwrap();
        }
        if !clusters.is_empty() && clusters.len() <= n && clusters.iter().all(|c| !c.is_empty()) {
            for (c, v) in ca.cohesion(&clusters, original).unwrap().iter().enumerate() {
                writeln!(out, "COH\t{tag}\t{k}\t{c}\t{}", fb(*v)).unwrap();
            }
        }
        if k < n {
            for (c, subtree) in ca.cut_subtrees(k, tree).unwrap().iter().enumerate() {
                for (j, node) in subtree.iter().enumerate() {
                    writeln!(
                        out,
                        "SUB\t{tag}\t{k}\t{c}\t{j}\t{}\t{}\t{}",
                        node.left_child,
                        node.right_child,
                        fb(node.distance)
                    )
                    .unwrap();
                }
            }
        }
    }
    for k in 1..n {
        writeln!(
            out,
            "APA\t{tag}\t{k}\t{}",
            fb(ca.average_population_aberration(k, tree).unwrap())
        )
        .unwrap();
    }
    writeln!(out, "NWK0\t{tag}\t{}", ca.newick_tree(tree, None).unwrap()).unwrap();
    writeln!(
        out,
        "NWK1\t{tag}\t{}",
        ca.newick_tree(tree, Some(&to_str_f32)).unwrap()
    )
    .unwrap();
}

const MATRIX_CASES: [(&str, usize, u64, u32); 9] = [
    ("m2", 2, 1, 4),
    ("m3", 3, 2, 3),
    ("m6", 6, 3, 5),
    ("m12", 12, 4, 3),
    ("m17", 17, 5, 4),
    ("m40t", 40, 6, 6),
    ("m40c", 40, 7, 0),
    ("m100t", 100, 8, 10),
    ("m100c", 100, 9, 0),
];

#[test]
fn linkage_and_analyzer_match_release_build() {
    let mut out = String::new();
    for (name, n, seed, levels) in MATRIX_CASES {
        let original = make_matrix(n, seed, levels);
        let mut m = original.clone();
        let tree = SingleLinkage::new().cluster(&mut m, 1.0).unwrap();
        print_tree(&mut out, &format!("L\t{name}\tsingle\t{}", fb(1.0)), &tree);
        if n <= 40 {
            analyze(&mut out, &format!("{name}/single"), &tree, &original);
        }
        for threshold in [1.0_f32, 0.5, 0.25] {
            let mut m = original.clone();
            let tree = CompleteLinkage::new().cluster(&mut m, threshold).unwrap();
            print_tree(
                &mut out,
                &format!("L\t{name}\tcomplete\t{}", fb(threshold)),
                &tree,
            );
            if n <= 40 && threshold == 0.5 {
                analyze(&mut out, &format!("{name}/complete0.5"), &tree, &original);
            }
        }
        for threshold in [1.0_f32, 0.5, 0.25] {
            let mut m = original.clone();
            let tree = AverageLinkage::new().cluster(&mut m, threshold).unwrap();
            print_tree(
                &mut out,
                &format!("L\t{name}\taverage\t{}", fb(threshold)),
                &tree,
            );
            if n <= 40 && threshold == 1.0 {
                analyze(&mut out, &format!("{name}/average1"), &tree, &original);
            }
        }
    }
    compare("linkage", &out);
}

// ---------------------------------------------------------------- grid

fn multiplex_distance(rt_scaling: f64) -> impl Fn(&DPosition<2>, &DPosition<2>) -> f64 {
    move |p1, p2| {
        let dx = p1.x() - p2.x();
        let dy = p1.y() - p2.y();
        (dx * dx + rt_scaling * rt_scaling * dy * dy).sqrt()
    }
}

fn print_grid(out: &mut String, tag: &str, results: &BTreeMap<i32, GridBasedCluster>) {
    for (index, c) in results {
        let b = c.bounding_box();
        writeln!(
            out,
            "G\t{tag}\t{index}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            db(c.centre().x()),
            db(c.centre().y()),
            db(b.min_x()),
            db(b.min_y()),
            db(b.max_x()),
            db(b.max_y()),
            c.property_a(),
            join(c.points()),
            join(c.properties_b())
        )
        .unwrap();
    }
}

fn points(seed: u64, count: usize, quantized: bool) -> (Vec<f64>, Vec<f64>) {
    let mut r = Lcg(seed);
    let (mut x, mut y) = (Vec::new(), Vec::new());
    for _ in 0..count {
        let mut a = 10.0 * r.unit();
        let mut b = 10.0 * r.unit();
        if quantized {
            a = (a * 4.0).floor() / 4.0;
            b = (b * 4.0).floor() / 4.0;
        }
        x.push(a);
        y.push(b);
    }
    (x, y)
}

fn unit_lines() -> Vec<f64> {
    (0..=10).map(f64::from).collect()
}

#[test]
fn grid_based_clustering_matches_release_build() {
    let mut out = String::new();
    {
        let (x, y) = points(11, 2000, false);
        let mut g = GridBasedClustering::without_properties(
            multiplex_distance(1.0),
            &x,
            &y,
            unit_lines(),
            unit_lines(),
        )
        .unwrap();
        g.cluster().unwrap();
        print_grid(&mut out, "G1/cluster", g.results());
        g.extend_clusters_y().unwrap();
        print_grid(&mut out, "G1/extend", g.results());
        g.remove_small_clusters_y(0.5);
        print_grid(&mut out, "G1/remove", g.results());
    }
    {
        let (x, y) = points(12, 1500, true);
        let mut r = Lcg(13);
        let (mut a, mut b) = (Vec::new(), Vec::new());
        for _ in 0..x.len() {
            a.push(1 + r.below(2) as i32);
            b.push(r.below(5) as i32);
        }
        let mut g = GridBasedClustering::new(
            multiplex_distance(1.0),
            &x,
            &y,
            &a,
            &b,
            unit_lines(),
            unit_lines(),
        )
        .unwrap();
        g.cluster().unwrap();
        print_grid(&mut out, "G2/cluster", g.results());
    }
    {
        let (x, y) = points(14, 1500, true);
        let mut g = GridBasedClustering::without_properties(
            multiplex_distance(1.0),
            &x,
            &y,
            unit_lines(),
            unit_lines(),
        )
        .unwrap();
        g.cluster().unwrap();
        print_grid(&mut out, "G3/cluster", g.results());
    }
    {
        let (x, y) = points(15, 800, false);
        let gx = vec![0.0, 0.5, 1.5, 3.0, 5.0, 8.0, 10.0];
        let gy = vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0];
        let mut g =
            GridBasedClustering::without_properties(multiplex_distance(2.0), &x, &y, gx, gy)
                .unwrap();
        g.cluster().unwrap();
        print_grid(&mut out, "G4/cluster", g.results());
        g.extend_clusters_y().unwrap();
        print_grid(&mut out, "G4/extend", g.results());
        g.remove_small_clusters_y(1.0);
        print_grid(&mut out, "G4/remove", g.results());
    }
    compare("grid", &out);
}

// ---------------------------------------------------------------- QTCluster

fn feature(rt: f64, mz: f64, sequences: &[&str]) -> BaseFeature {
    let mut f = BaseFeature::new(rt, mz, 0.0);
    for sequence in sequences {
        let mut identification = PeptideIdentification::default();
        identification.hits.push(PeptideHit {
            sequence: AASequence::parse(sequence).unwrap(),
            ..PeptideHit::default()
        });
        f.peptide_identifications.push(identification);
    }
    f
}

fn print_cluster(out: &mut String, tag: &str, c: &mut QTCluster<'_>) {
    writeln!(out, "QSIZE\t{tag}\t{}", c.size()).unwrap();
    writeln!(out, "QQUAL\t{tag}\t{}", db(c.quality())).unwrap();
    let annotations: Vec<&str> = c.annotations().iter().map(|s| s.0.as_str()).collect();
    writeln!(out, "QANN\t{tag}\t{}", join(&annotations)).unwrap();
    let mut neighbors: Vec<(usize, usize)> = c
        .all_neighbors()
        .iter()
        .map(|e| (e.map_index, e.feature.feature_index()))
        .collect();
    neighbors.sort_unstable();
    let text: Vec<String> = neighbors.iter().map(|(m, f)| format!("{m}:{f}")).collect();
    writeln!(out, "QNB\t{tag}\t{}", join(&text)).unwrap();
}

#[test]
fn qt_cluster_matches_release_build() {
    let mut out = String::new();
    {
        let fc = feature(10.0, 500.0, &["AAA", "CC"]);
        let f1 = feature(10.1, 500.1, &["AAA"]);
        let f2 = feature(10.2, 500.2, &["CC"]);
        let f3 = feature(10.3, 500.3, &[]);
        let f5 = feature(10.5, 500.5, &["CC"]);
        let f6 = feature(10.6, 500.6, &["AAA"]);
        let c = GridFeature::new(&fc, 0, 100);
        let g1 = GridFeature::new(&f1, 1, 101);
        let g2 = GridFeature::new(&f2, 1, 102);
        let g3 = GridFeature::new(&f3, 2, 103);
        let g5 = GridFeature::new(&f5, 3, 105);
        let g6 = GridFeature::new(&f6, 3, 106);
        let mut q = QTCluster::new(BulkData::new(&c, 4, 3.0, 0, 0, 0), true);
        q.initialize_cluster();
        q.add(&g1, 1.0);
        q.add(&g2, 1.0);
        q.add(&g3, 2.5);
        q.add(&g5, 0.3);
        q.add(&g6, 0.3);
        q.finalize_cluster();
        print_cluster(&mut out, "Q1/final", &mut q);
        let removed = [
            Element {
                map_index: 3,
                feature: &g5,
            },
            Element {
                map_index: 3,
                feature: &g6,
            },
        ];
        writeln!(out, "QUPD\tQ1\t{}", u8::from(q.update(&removed))).unwrap();
        print_cluster(&mut out, "Q1/updated", &mut q);
    }
    {
        let fc = feature(1.0, 100.0, &[]);
        let fs: Vec<BaseFeature> = (0..12)
            .map(|i| feature(1.0 + f64::from(i), 100.0 + f64::from(i), &[]))
            .collect();
        let gs: Vec<GridFeature<'_>> = fs
            .iter()
            .enumerate()
            .map(|(i, f)| GridFeature::new(f, 1 + (i % 5), 200 + i))
            .collect();
        let c = GridFeature::new(&fc, 0, 199);
        let mut q = QTCluster::new(BulkData::new(&c, 6, 0.9, 1, 2, 7), false);
        q.initialize_cluster();
        let mut r = Lcg(21);
        for g in &gs {
            q.add(g, 0.9 * r.unit());
        }
        q.finalize_cluster();
        print_cluster(&mut out, "Q2/final", &mut q);
        let removed = [
            Element {
                map_index: 2,
                feature: &gs[1],
            },
            Element {
                map_index: 2,
                feature: &gs[6],
            },
        ];
        writeln!(out, "QUPD\tQ2\t{}", u8::from(q.update(&removed))).unwrap();
        print_cluster(&mut out, "Q2/updated", &mut q);
    }
    compare("qt", &out);
}

// ---------------------------------------------------------------- HashGrid

#[test]
fn hash_grid_cell_indices_match_release_build() {
    let g = HashGrid::<i32>::new(DPosition::xy(0.3, 7.0));
    let mut r = Lcg(31);
    let mut out = String::new();
    for i in 0..50 {
        let a = 200.0 * r.unit() - 100.0;
        let b = 200.0 * r.unit() - 100.0;
        let index = g
            .cell_index_at_cluster_center(&DPosition::xy(a, b))
            .unwrap();
        writeln!(out, "H\t{i}\t{}\t{}", index[0], index[1]).unwrap();
    }
    compare("hashgrid", &out);
}
