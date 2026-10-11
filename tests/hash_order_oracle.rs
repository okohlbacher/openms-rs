// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Evidence tier 1 for the hash-container orders of
//! [`openms::math::hash_order`] and the iteration order of
//! [`openms::math::ml::clustering::hash_grid::HashGrid`]: a differential
//! against the executed C++ Release build.
//!
//! Oracle `M2-hashorder` (`oracle/M2-hashorder/driver.cpp`) ran on ibminode06
//! (Linux x86_64) against the Release install
//! `openms4-release-bc9cc12-c19e494-174b576` (conda-forge g++ 14.4.0, libstdc++
//! 14.4.0, Boost 1.89), twice, with identical output. `hash_order_scripts.txt`
//! holds 14 generated operation scripts (`oracle/M2-hashorder/
//! gen_scripts.py`, seed 20261013):
//!
//! - `std_*`: inserts, erasures, clears and order prints on
//!   `std::unordered_map<size_t, int>` (dense small key spaces that stay in
//!   one bucket group, wide key spaces, sequential keys, growth through 600
//!   insertions, clear-heavy histories, keys up to 5 million);
//! - `grid_*`: inserts (with repeated coordinates), erasures and prints on
//!   `OpenMS::HashGrid<int>`: the cells in grid order and every pair in
//!   iteration order, for several cell sizes, negative coordinates and
//!   integer lattices.
//!
//! `hash_order_oracle.tsv` has the printed orders; each print is compared
//! whole: element count, bucket count (libstdc++ scripts) and every element in
//! order.

use std::collections::BTreeMap;

use openms::data_structures::dposition::DPosition;
use openms::math::hash_order::StdHashMap;
use openms::math::ml::clustering::hash_grid::HashGrid;

const DIR: &str = "tests/data";

fn scripts() -> BTreeMap<String, Vec<String>> {
    let text = std::fs::read_to_string(format!("{DIR}/hash_order_scripts.txt")).unwrap();
    let mut out: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut name = String::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("script\t") {
            name = rest.to_string();
            out.insert(name.clone(), Vec::new());
        } else {
            out.get_mut(&name).unwrap().push(line.to_string());
        }
    }
    out
}

fn oracle() -> BTreeMap<String, Vec<String>> {
    let text = std::fs::read_to_string(format!("{DIR}/hash_order_oracle.tsv")).unwrap();
    let mut out: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut name = String::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("case\t") {
            name = rest.to_string();
            out.insert(name.clone(), Vec::new());
        } else {
            out.get_mut(&name).unwrap().push(line.to_string());
        }
    }
    out
}

fn run_std(script: &[String]) -> Vec<String> {
    let mut map: StdHashMap<i32> = StdHashMap::new();
    let mut prints = Vec::new();
    for line in script {
        let f: Vec<&str> = line.split_whitespace().collect();
        match f[0] {
            "i" => {
                map.insert(f[1].parse().unwrap(), f[2].parse().unwrap());
            }
            "e" => {
                map.remove(f[1].parse().unwrap());
            }
            "c" => map.clear(),
            "p" => {
                let mut text = format!("order\t{}\t{}", map.len(), map.bucket_count());
                for (key, value) in map.iter() {
                    text.push_str(&format!("\t{key}:{value}"));
                }
                prints.push(text);
            }
            other => panic!("unknown op {other}"),
        }
    }
    prints
}

fn run_grid(script: &[String]) -> Vec<String> {
    let dims: Vec<f64> = script[0]
        .split_whitespace()
        .map(|x| x.parse().unwrap())
        .collect();
    let mut grid: HashGrid<i32> = HashGrid::new(DPosition::xy(dims[0], dims[1]));
    let mut prints = Vec::new();
    for line in &script[1..] {
        let f: Vec<&str> = line.split_whitespace().collect();
        match f[0] {
            "i" => {
                let key = DPosition::xy(f[1].parse().unwrap(), f[2].parse().unwrap());
                grid.insert(key, f[3].parse().unwrap()).unwrap();
            }
            "e" => {
                let key = DPosition::xy(f[1].parse().unwrap(), f[2].parse().unwrap());
                grid.erase(&key).unwrap();
            }
            "p" => {
                let mut cells = format!("cells\t{}", grid.len());
                for (index, content) in grid.grid_iter() {
                    cells.push_str(&format!("\t{},{}:{}", index[0], index[1], content.len()));
                }
                prints.push(cells);
                let mut pairs = String::from("pairs");
                for (index, (center, value)) in grid.iter() {
                    pairs.push_str(&format!(
                        "\t{},{}:{},{}:{}",
                        index[0],
                        index[1],
                        fmt_g17(center[0]),
                        fmt_g17(center[1]),
                        value
                    ));
                }
                prints.push(pairs);
            }
            other => panic!("unknown op {other}"),
        }
    }
    prints
}

/// The oracle's `%.17g` text of a coordinate, parsed and printed in Rust's
/// shortest round-trip form (both sides are normalised through `f64`).
fn fmt_g17(value: f64) -> String {
    format!("{value:?}")
}

/// Normalise an oracle `pairs` line: the coordinates through `f64`.
fn normalise(line: &str) -> String {
    let mut fields = line.split('\t');
    let tag = fields.next().unwrap();
    let mut out = String::from(tag);
    for field in fields {
        if tag != "pairs" {
            out.push('\t');
            out.push_str(field);
            continue;
        }
        let parts: Vec<&str> = field.split(':').collect();
        let xy: Vec<f64> = parts[1].split(',').map(|v| v.parse().unwrap()).collect();
        out.push_str(&format!(
            "\t{}:{},{}:{}",
            parts[0],
            fmt_g17(xy[0]),
            fmt_g17(xy[1]),
            parts[2]
        ));
    }
    out
}

fn check(name: &str) {
    let all_scripts = scripts();
    let all_oracle = oracle();
    let script = &all_scripts[name];
    let expected = &all_oracle[name];
    let produced = if name.starts_with("std_") {
        run_std(script)
    } else {
        run_grid(script)
    };
    assert_eq!(produced.len(), expected.len(), "{name}: number of prints");
    for (index, (got, want)) in produced.iter().zip(expected).enumerate() {
        assert_eq!(got, &normalise(want), "{name}: print {index}");
    }
}

macro_rules! order_test {
    ($test:ident, $name:literal) => {
        #[test]
        fn $test() {
            check($name);
        }
    };
}

order_test!(std_small_dense, "std_small_dense");
order_test!(std_small_wide, "std_small_wide");
order_test!(std_grow, "std_grow");
order_test!(std_grow_seq, "std_grow_seq");
order_test!(std_clear_heavy, "std_clear_heavy");
order_test!(std_mapindex, "std_mapindex");
order_test!(std_big, "std_big");
order_test!(std_seq, "std_seq");
order_test!(grid_small, "grid_small");
order_test!(grid_dense_cells, "grid_dense_cells");
order_test!(grid_sparse, "grid_sparse");
order_test!(grid_many, "grid_many");
order_test!(grid_negative, "grid_negative");
order_test!(grid_lattice, "grid_lattice");

#[test]
fn every_script_has_a_test() {
    let names: Vec<String> = scripts().keys().cloned().collect();
    assert_eq!(names.len(), 14);
    assert_eq!(names, oracle().keys().cloned().collect::<Vec<_>>());
}
