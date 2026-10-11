// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Evidence tier 1 for `WNetMatcher` and the Wasserstein aligner
//! (`wnetalign::WNetAligner<2>` over `wnet::WassersteinNetwork` and LEMON's
//! network simplex) beneath it: a differential against executed C++.
//!
//! Oracle `M2-wnet` (`oracle/M2-wnet/driver.cpp`) ran on ibminode06
//! (Linux x86_64, conda-forge g++ 14.4.0 with libOpenMS's core flags
//! `-O3 -DNDEBUG -std=gnu++23 -mssse3 -ffp-contract=off`), twice, with
//! identical standard output, error output and exit codes. The Release
//! install `openms4-release-bc9cc12-c19e494-174b576` was configured without
//! WNet, so the oracle compiles the pinned `WNetMatcher.cpp` and
//! `FeatureGroupingAlgorithmWNet.cpp` (core bc9cc12) against that install and
//! the header-only libraries the pin fetches at their pinned commits: wnetalign
//! v0.9.8 `cff6a19`, wnet v0.9.11 `18a1525`, pylmcf v0.9.8 `d2c9c52` (LEMON).
//!
//! `wnet_instances.txt` holds 31 generated instances (`oracle/M2-wnet/
//! gen_instances.py`, seed 20261011): random continuous point sets under all
//! three metrics, integer lattices with equal weights (massive cost and flow
//! ties, which is where the simplex's pivot order and `std::sort`'s order of
//! equal flows are visible), duplicates, zero intensities, cluster sets with
//! many components, trash costs from below to above every distance, a zero
//! and a huge maximum distance, and realistic scaled `ln(m/z)` positions with
//! normalised intensities. `wnet_oracle_instances.tsv` has, per instance, the
//! scale factor, the number of components, the total cost, every flow in the
//! order `flows_for_target` returns it, the greedy consensus in selection
//! order, and the result of `WNetMatcher::match` (printed with `%.17g`).
//!
//! Compared exactly: the scale factor and the total cost bit for bit, the
//! component count, every flow triple in order, the consensus pairs in order.

use std::collections::BTreeMap;

use openms::analysis::mapmatching::wnet_aligner::{Spectrum, WNetAligner};
use openms::analysis::mapmatching::wnet_matcher::{DistanceMetric, WNetMatcher};

const DIR: &str = "tests/data";

struct Instance {
    metric: DistanceMetric,
    max_distance: f64,
    trash_cost: f64,
    a: Vec<[f64; 2]>,
    intensity_a: Vec<f64>,
    b: Vec<[f64; 2]>,
    intensity_b: Vec<f64>,
}

fn instances() -> BTreeMap<String, Instance> {
    let text = std::fs::read_to_string(format!("{DIR}/wnet_instances.txt")).unwrap();
    let mut out = BTreeMap::new();
    let mut name = String::new();
    let mut current: Option<Instance> = None;
    for line in text.lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        match fields[0] {
            "instance" => {
                if let Some(done) = current.take() {
                    out.insert(name.clone(), done);
                }
                name = fields[1].to_string();
            }
            "params" => {
                let metric = match fields[1] {
                    "L1" => DistanceMetric::L1,
                    "L2" => DistanceMetric::L2,
                    _ => DistanceMetric::Linf,
                };
                current = Some(Instance {
                    metric,
                    max_distance: fields[2].parse().unwrap(),
                    trash_cost: fields[3].parse().unwrap(),
                    a: Vec::new(),
                    intensity_a: Vec::new(),
                    b: Vec::new(),
                    intensity_b: Vec::new(),
                });
            }
            tag @ ("A" | "B") => {
                let instance = current.as_mut().unwrap();
                let x: f64 = fields[1].parse().unwrap();
                let y: f64 = fields[2].parse().unwrap();
                let i: f64 = fields[3].parse().unwrap();
                if tag == "A" {
                    instance.a.push([x, y]);
                    instance.intensity_a.push(i);
                } else {
                    instance.b.push([x, y]);
                    instance.intensity_b.push(i);
                }
            }
            other => panic!("unexpected tag {other}"),
        }
    }
    if let Some(done) = current.take() {
        out.insert(name, done);
    }
    out
}

#[derive(Default)]
struct Oracle {
    scale: f64,
    nsub: usize,
    cost: f64,
    flows: Vec<(usize, usize, i64)>,
    consensus: Vec<(usize, usize)>,
    matcher_cost: f64,
    matcher_pairs: Vec<(usize, usize)>,
}

fn oracle() -> BTreeMap<String, Oracle> {
    let text = std::fs::read_to_string(format!("{DIR}/wnet_oracle_instances.tsv")).unwrap();
    let mut out: BTreeMap<String, Oracle> = BTreeMap::new();
    let mut name = String::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f[0] == "case" {
            name = f[1].to_string();
            out.insert(name.clone(), Oracle::default());
            continue;
        }
        let o = out.get_mut(&name).unwrap();
        match f[0] {
            "scale" => o.scale = f[1].parse().unwrap(),
            "nsub" => o.nsub = f[1].parse().unwrap(),
            "cost" => o.cost = f[1].parse().unwrap(),
            "flow" => o.flows.push((
                f[1].parse().unwrap(),
                f[2].parse().unwrap(),
                f[3].parse().unwrap(),
            )),
            "cons" => o
                .consensus
                .push((f[1].parse().unwrap(), f[2].parse().unwrap())),
            "mcost" => o.matcher_cost = f[1].parse().unwrap(),
            "mpair" => o
                .matcher_pairs
                .push((f[1].parse().unwrap(), f[2].parse().unwrap())),
            "nflow" | "ncons" | "mpairs" => {}
            other => panic!("unexpected tag {other}"),
        }
    }
    out
}

fn check(name: &str, instance: &Instance, expected: &Oracle) {
    let a = Spectrum::new(instance.a.clone(), instance.intensity_a.clone()).unwrap();
    let b = Spectrum::new(instance.b.clone(), instance.intensity_b.clone()).unwrap();
    let mut aligner = WNetAligner::new(
        &a,
        &[&b],
        instance.metric,
        instance.max_distance,
        instance.trash_cost,
    )
    .unwrap();
    aligner.set_point(&[1.0]).unwrap();
    assert_eq!(
        aligner.scale_factor().to_bits(),
        expected.scale.to_bits(),
        "{name}: scale factor"
    );
    assert_eq!(aligner.no_subgraphs(), expected.nsub, "{name}: components");
    assert_eq!(
        aligner.total_cost().unwrap().to_bits(),
        expected.cost.to_bits(),
        "{name}: total cost"
    );
    let (e, t, f) = aligner.flows_for_target(0).unwrap();
    let flows: Vec<(usize, usize, i64)> = e
        .iter()
        .zip(&t)
        .zip(&f)
        .map(|((e, t), f)| (*e, *t, *f))
        .collect();
    assert_eq!(flows, expected.flows, "{name}: flows");
    let (ce, ct) = aligner.consensus_for_target(0).unwrap();
    let consensus: Vec<(usize, usize)> = ce.into_iter().zip(ct).collect();
    assert_eq!(consensus, expected.consensus, "{name}: consensus");

    let result = WNetMatcher::match_points(
        &instance.a,
        &instance.intensity_a,
        &instance.b,
        &instance.intensity_b,
        instance.metric,
        instance.max_distance,
        instance.trash_cost,
    )
    .unwrap();
    assert_eq!(
        result.matched_pairs, expected.matcher_pairs,
        "{name}: pairs"
    );
    assert_eq!(
        result.cost.to_bits(),
        expected.matcher_cost.to_bits(),
        "{name}: matcher cost"
    );
}

#[test]
fn oracle_covers_every_instance() {
    let inst = instances();
    let orc = oracle();
    assert_eq!(inst.len(), 31);
    assert_eq!(
        inst.keys().collect::<Vec<_>>(),
        orc.keys().collect::<Vec<_>>()
    );
}

macro_rules! oracle_case {
    ($test:ident, $name:literal) => {
        #[test]
        fn $test() {
            let inst = instances();
            let orc = oracle();
            check($name, &inst[$name], &orc[$name]);
        }
    };
}

oracle_case!(oracle_big_intensity, "big_intensity");
oracle_case!(oracle_class_close, "class_close");
oracle_case!(oracle_cluster_large, "cluster_large");
oracle_case!(oracle_cluster_large_l2, "cluster_large_l2");
oracle_case!(oracle_duplicates, "duplicates");
oracle_case!(oracle_grid6_same, "grid6_same");
oracle_case!(oracle_grid6_shift, "grid6_shift");
oracle_case!(oracle_identical, "identical");
oracle_case!(oracle_lattice_equal_l1, "lattice_equal_l1");
oracle_case!(oracle_lattice_equal_l2, "lattice_equal_l2");
oracle_case!(oracle_lattice_equal_linf, "lattice_equal_linf");
oracle_case!(oracle_lattice_int, "lattice_int");
oracle_case!(oracle_lattice_int_linf, "lattice_int_linf");
oracle_case!(oracle_many_to_one, "many_to_one");
oracle_case!(oracle_maxd_huge, "maxd_huge");
oracle_case!(oracle_maxd_zero, "maxd_zero");
oracle_case!(oracle_one_to_many, "one_to_many");
oracle_case!(oracle_ppm_normalised, "ppm_normalised");
oracle_case!(oracle_ppm_normalised_l1, "ppm_normalised_l1");
oracle_case!(oracle_ppm_normalised_large, "ppm_normalised_large");
oracle_case!(oracle_rand_dense, "rand_dense");
oracle_case!(oracle_rand_mid_l1, "rand_mid_l1");
oracle_case!(oracle_rand_mid_l2, "rand_mid_l2");
oracle_case!(oracle_rand_mid_linf, "rand_mid_linf");
oracle_case!(oracle_rand_small_l1, "rand_small_l1");
oracle_case!(oracle_rand_small_l2, "rand_small_l2");
oracle_case!(oracle_rand_small_linf, "rand_small_linf");
oracle_case!(oracle_rand_sparse, "rand_sparse");
oracle_case!(oracle_trash_large, "trash_large");
oracle_case!(oracle_trash_small, "trash_small");
oracle_case!(oracle_zero_intensity, "zero_intensity");

/// `WNetMatcher::metricFromString`, oracle case `metric`.
#[test]
fn oracle_metric_from_string() {
    // metric L1 0, L2 1, LINF 2, unknown 2, "" 2 (enum values); the Rust
    // enum is compared by name.
    assert_eq!(WNetMatcher::metric_from_string("L1"), DistanceMetric::L1);
    assert_eq!(WNetMatcher::metric_from_string("L2"), DistanceMetric::L2);
    assert_eq!(
        WNetMatcher::metric_from_string("LINF"),
        DistanceMetric::Linf
    );
    assert_eq!(
        WNetMatcher::metric_from_string("unknown"),
        DistanceMetric::Linf
    );
    assert_eq!(WNetMatcher::metric_from_string(""), DistanceMetric::Linf);
}
