// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `MapAlignmentAlgorithmKD`. The upstream class test at core pin `bc9cc12`
//! (`MapAlignmentAlgorithmKD_test.cpp`) constructs the algorithm with 42 maps
//! and an empty `Param` and marks every other member `NOT_TESTABLE`; that
//! construction is reproduced (tier 3). The rest is derived from the pinned
//! source (tier 4): connected components by breadth-first search, the
//! size, conflict and (dead) charge filters, the average-RT fit data, the
//! identity fallback below 50 points and the refusals. The executed Release
//! build is compared in `tests/map_alignment_oracle_kd.rs`.

#[path = "support/kd_pool.rs"]
mod kd_pool;

use kd_pool::{Pool, PoolFeature};
use openms::analysis::mapmatching::alignment_kd::{KdFeatureMaps, MapAlignmentAlgorithmKD};
use openms::param::{Param, ParamValue};

fn feature(map: usize, rt: f64, mz: f64) -> PoolFeature {
    PoolFeature {
        map,
        rt,
        mz,
        intensity: 1000.0,
        charge: 2,
    }
}

/// The parameters `FeatureGroupingAlgorithmKD` passes, with the given warp
/// tolerances and filters.
fn params(rt_tol: f64, mz_tol: f64, unit: &str, log_fc: f64, rel: f64, conflicts: i64) -> Param {
    let mut p = Param::new();
    let mut set = |k: &str, v: ParamValue| p.set_value(k, v, "", &[]).unwrap();
    set("warp:rt_tol", ParamValue::Float(rt_tol));
    set("warp:mz_tol", ParamValue::Float(mz_tol));
    set("mz_unit", ParamValue::String(unit.into()));
    set("warp:max_pairwise_log_fc", ParamValue::Float(log_fc));
    set("warp:min_rel_cc_size", ParamValue::Float(rel));
    set("warp:max_nr_conflicts", ParamValue::Integer(conflicts));
    set("LOWESS:span", ParamValue::Float(0.666_666_666_666_666_6));
    p
}

/// Three maps; peptide `i` at m/z `400 + 10 i` elutes at `100 + 20 i` plus
/// a map shift of 0, +5 and -5 seconds.
fn shifted_pool(peptides: usize) -> Pool {
    let shifts = [0.0, 5.0, -5.0];
    let mut features = Vec::new();
    for (map, shift) in shifts.iter().enumerate() {
        for i in 0..peptides {
            features.push(feature(
                map,
                100.0 + 20.0 * i as f64 + shift,
                400.0 + 10.0 * i as f64,
            ));
        }
    }
    Pool::new(features)
}

#[test]
fn class_test_construction_with_empty_param() {
    let algorithm = MapAlignmentAlgorithmKD::new(42, &Param::new()).unwrap();
    assert_eq!(algorithm.num_maps(), 42);
    assert!(algorithm.fit_data().iter().all(Vec::is_empty));
    assert!(algorithm.transformations().is_empty());
}

#[test]
fn components_and_fit_data_of_shifted_maps() {
    let pool = shifted_pool(4);
    let mut algorithm =
        MapAlignmentAlgorithmKD::new(3, &params(12.0, 10.0, "ppm", 0.5, 0.5, 0)).unwrap();
    // Features 0..4 are map 0, 4..8 map 1, 8..12 map 2; peptide i is
    // {i, 4 + i, 8 + i}, and the first unvisited feature seeds each
    // component.
    let (labels, count) = algorithm.compute_ccs(&pool).unwrap();
    assert_eq!(count, 4);
    assert_eq!(labels, vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]);
    let ccs = algorithm.get_ccs(&pool).unwrap();
    assert_eq!(ccs[&2], vec![2, 6, 10]);
    algorithm.add_rt_fit_data(&pool).unwrap();
    for (map, shift) in [0.0, 5.0, -5.0].iter().enumerate() {
        let data = &algorithm.fit_data()[map];
        assert_eq!(data.len(), 4);
        for (i, point) in data.iter().enumerate() {
            let rt = 100.0 + 20.0 * i as f64;
            assert_eq!(point.x, rt + shift);
            // The component average: (rt + (rt + 5) + (rt - 5)) / 3.
            assert_eq!(point.y, ((rt + (rt + 5.0)) + (rt - 5.0)) / 3.0);
        }
    }
}

#[test]
fn the_rt_window_bounds_the_components() {
    // Map 1 is 5 s late: a 4 s tolerance separates the maps entirely, and
    // map 2 (-5 s) is 10 s from map 1.
    let pool = shifted_pool(3);
    let algorithm =
        MapAlignmentAlgorithmKD::new(3, &params(4.0, 10.0, "ppm", -1.0, 0.0, -1)).unwrap();
    let (_, count) = algorithm.compute_ccs(&pool).unwrap();
    assert_eq!(count, 9);
    // At exactly 5 s the window is closed, so maps 0-1 and 0-2 connect.
    let algorithm =
        MapAlignmentAlgorithmKD::new(3, &params(5.0, 10.0, "ppm", -1.0, 0.0, -1)).unwrap();
    let (labels, count) = algorithm.compute_ccs(&pool).unwrap();
    assert_eq!(count, 3);
    assert_eq!(labels, vec![0, 1, 2, 0, 1, 2, 0, 1, 2]);
}

#[test]
fn size_and_conflict_filters() {
    // Map 0 holds two features of the first peptide 1 s apart; map 1 one.
    let pool = Pool::new(vec![
        feature(0, 100.0, 500.0),
        feature(0, 101.0, 500.0),
        feature(1, 100.5, 500.0),
        feature(1, 300.0, 600.0),
    ]);
    let ccs_of = |conflicts: i64, rel: f64| {
        let algorithm =
            MapAlignmentAlgorithmKD::new(2, &params(2.0, 0.01, "Da", -1.0, rel, conflicts))
                .unwrap();
        let ccs = algorithm.get_ccs(&pool).unwrap();
        assert_eq!(ccs[&0], vec![0, 1, 2]);
        assert_eq!(ccs[&1], vec![3]);
        algorithm
            .filter_ccs(&pool, &ccs)
            .unwrap()
            .into_keys()
            .collect::<Vec<_>>()
    };
    // The singleton is below the minimum of 2; the triple has one conflict.
    assert_eq!(ccs_of(0, 0.0), Vec::<usize>::new());
    assert_eq!(ccs_of(1, 0.0), vec![0]);
    assert_eq!(ccs_of(-1, 0.0), vec![0]);
    // max(2, 1.6 * 2 maps) = 3.2, truncated to 3.
    assert_eq!(ccs_of(-1, 1.6), vec![0]);
    // max(2, 2.0 * 2) = 4.
    assert_eq!(ccs_of(-1, 2.0), Vec::<usize>::new());
}

#[test]
fn mixed_charges_are_kept_as_in_the_source() {
    // The source's `continue` on a second charge only continues its inner
    // loop, so a component of charges 2 and 3 passes.
    let mut a = feature(0, 100.0, 500.0);
    a.charge = 2;
    let mut b = feature(1, 100.0, 500.0);
    b.charge = 3;
    let pool = Pool::new(vec![a, b]);
    let algorithm = MapAlignmentAlgorithmKD::new(2, &params(1.0, 1.0, "Da", -1.0, 0.0, 0)).unwrap();
    let ccs = algorithm.get_ccs(&pool).unwrap();
    assert_eq!(algorithm.filter_ccs(&pool, &ccs).unwrap().len(), 1);
    assert_eq!(pool.charge(0), 2);
}

#[test]
fn the_fold_change_filter_breaks_edges() {
    let mut a = feature(0, 100.0, 500.0);
    a.intensity = 100.0;
    let mut b = feature(1, 100.0, 500.0);
    b.intensity = 1000.0;
    let pool = Pool::new(vec![a, b]);
    // |log10(10)| = 1.
    for (log_fc, count) in [(0.5, 2), (1.0, 1), (-1.0, 1)] {
        let algorithm =
            MapAlignmentAlgorithmKD::new(2, &params(1.0, 1.0, "Da", log_fc, 0.0, 0)).unwrap();
        assert_eq!(algorithm.compute_ccs(&pool).unwrap().1, count, "{log_fc}");
    }
}

#[test]
fn lowess_fit_and_transform() {
    // 60 peptides: enough points for every map.
    let mut pool = shifted_pool(60);
    let mut algorithm =
        MapAlignmentAlgorithmKD::new(3, &params(12.0, 10.0, "ppm", 0.5, 0.5, 0)).unwrap();
    assert!(algorithm.transform(&mut pool).is_err());
    algorithm.add_rt_fit_data(&pool).unwrap();
    assert!(algorithm.fit_data().iter().all(|d| d.len() == 60));
    assert!(algorithm.fit_lowess().unwrap().is_empty());
    algorithm.transform(&mut pool).unwrap();
    // Each map's shift is constant, so the smoothed fit is the shift back to
    // the average of the three maps, which is the unshifted RT.
    for i in 0..pool.len() {
        let expected = 100.0 + 20.0 * (i % 60) as f64;
        assert!((pool.rt(i) - expected).abs() < 1e-6, "{i}: {}", pool.rt(i));
    }
}

#[test]
fn too_few_points_fall_back_to_the_identity() {
    let mut pool = shifted_pool(4);
    let mut algorithm =
        MapAlignmentAlgorithmKD::new(3, &params(12.0, 10.0, "ppm", 0.5, 0.5, 0)).unwrap();
    algorithm.add_rt_fit_data(&pool).unwrap();
    let warnings = algorithm.fit_lowess().unwrap();
    assert_eq!(warnings.len(), 3);
    assert_eq!(
        warnings[1],
        "Warning: Only 4 data points for LOWESS fit of map 1. Consider adjusting RT or m/z tolerance or max_pairwise_log_fc, decreasing min_rel_cc_size, or increasing max_nr_conflicts."
    );
    algorithm.transform(&mut pool).unwrap();
    for i in 0..pool.len() {
        assert!((pool.rt(i) - pool.features[i].rt).abs() < 1e-6);
    }
}

#[test]
fn fit_data_accumulates_over_calls() {
    let pool = shifted_pool(2);
    let mut algorithm =
        MapAlignmentAlgorithmKD::new(3, &params(12.0, 10.0, "ppm", 0.5, 0.5, 0)).unwrap();
    algorithm.add_rt_fit_data(&pool).unwrap();
    algorithm.add_rt_fit_data(&pool).unwrap();
    assert!(algorithm.fit_data().iter().all(|d| d.len() == 4));
}

#[test]
fn refusals() {
    let pool = shifted_pool(2);
    // An empty Param leaves the source's tolerances uninitialised.
    let mut algorithm = MapAlignmentAlgorithmKD::new(3, &Param::new()).unwrap();
    assert!(algorithm.add_rt_fit_data(&pool).is_err());
    // A missing key is the source's ElementNotFound.
    let mut partial = params(12.0, 10.0, "ppm", 0.5, 0.5, 0);
    partial.remove("warp:mz_tol").unwrap();
    assert!(MapAlignmentAlgorithmKD::new(3, &partial).is_err());
    // The filter reads its keys on use.
    let mut partial = params(12.0, 10.0, "ppm", 0.5, 0.5, 0);
    partial.remove("warp:max_nr_conflicts").unwrap();
    let mut algorithm = MapAlignmentAlgorithmKD::new(3, &partial).unwrap();
    assert!(algorithm.add_rt_fit_data(&pool).is_err());
    // A map index past the maps the algorithm was built for.
    let mut algorithm =
        MapAlignmentAlgorithmKD::new(2, &params(12.0, 10.0, "ppm", 0.5, 0.5, 0)).unwrap();
    assert!(algorithm.add_rt_fit_data(&pool).is_err());
}
