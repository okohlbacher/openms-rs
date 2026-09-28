// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Section-by-section ports of the grid-clustering class tests at source
//! revision bc9cc12514c768385ce121d6ca4bb710fe1983c4
//! (`src/tests/class_tests/openms/source/`): `ClusteringGrid_test.cpp`,
//! `GridBasedCluster_test.cpp`, `GridBasedClustering_test.cpp` and
//! `HashGrid_test.cpp`.
//!
//! Evidence tier 3: every expected value is a transcribed class-test literal,
//! or, where marked, an invariant of the source's documented contract. No
//! expected value was derived from this port's output.

use openms::data_structures::dinterval::DIntervalBase;
use openms::data_structures::dposition::DPosition;
use openms::math::ml::clustering::clustering_grid::ClusteringGrid;
use openms::math::ml::clustering::grid_based_cluster::GridBasedCluster;
use openms::math::ml::clustering::grid_based_clustering::{GridBasedClustering, MinimumDistance};
use openms::math::ml::clustering::hash_grid::{HashGrid, hash_value};

fn lines() -> Vec<f64> {
    (0..=10).map(f64::from).collect()
}

// ---------------------------------------------------------------- ClusteringGrid

#[test]
fn clustering_grid_class_test() {
    // ClusteringGrid(grid_spacing_x, grid_spacing_y)
    let grid = ClusteringGrid::new(lines(), lines()).unwrap();
    assert_eq!(grid.grid_spacing_x()[3], 3.0);

    let mut grid = ClusteringGrid::new(lines(), lines()).unwrap();
    let index1 = (2, 3);
    let index2 = (5, 4);
    let index3 = (7, 7);
    let point = DPosition::xy(6.6, 7.7);

    // getGridSpacingX / getGridSpacingY
    assert_eq!(grid.grid_spacing_x()[3], 3.0);
    assert_eq!(grid.grid_spacing_x()[10], 10.0);
    assert_eq!(grid.grid_spacing_y()[3], 3.0);
    assert_eq!(grid.grid_spacing_y()[10], 10.0);

    // addCluster: the sections share one grid, as in the class test.
    grid.add_cluster(index1, 1);
    grid.add_cluster(index2, 2);
    assert_eq!(grid.cell_count(), 2);

    // removeCluster
    grid.add_cluster(index1, 1);
    grid.add_cluster(index2, 2);
    grid.remove_cluster(index2, 2);
    assert_eq!(grid.cell_count(), 1);

    // removeAllClusters
    grid.add_cluster(index1, 1);
    grid.add_cluster(index2, 2);
    grid.remove_all_clusters();
    assert_eq!(grid.cell_count(), 0);

    // getClusters; an empty cell returns an empty list (regression).
    grid.add_cluster(index1, 1);
    grid.add_cluster(index2, 2);
    assert_eq!(grid.clusters(index1).first(), Some(&1));
    assert!(grid.clusters(index3).is_empty());

    // getIndex
    assert_eq!(grid.index(&point).unwrap().0, 7);
    assert_eq!(grid.index(&point).unwrap().1, 8);

    // isNonEmptyCell
    grid.add_cluster(index1, 1);
    assert!(grid.is_non_empty_cell(index1));
    assert!(!grid.is_non_empty_cell(index3));

    // getCellCount
    grid.add_cluster(index1, 1);
    grid.add_cluster(index2, 2);
    assert_eq!(grid.cell_count(), 2);
}

#[test]
fn clustering_grid_boundaries() {
    let grid = ClusteringGrid::new(lines(), lines()).unwrap();
    // Source `Exception::IllegalArgument` outside [first, last] line.
    assert!(grid.index(&DPosition::xy(10.5, 1.0)).is_err());
    assert!(grid.index(&DPosition::xy(1.0, -0.1)).is_err());
    // `std::lower_bound` with `std::less_equal`: a point on a line counts
    // the lines not greater than it.
    assert_eq!(grid.index(&DPosition::xy(0.0, 10.0)).unwrap(), (1, 11));
    // Empty grid lines: the source reads front() of an empty vector.
    assert!(ClusteringGrid::new(Vec::new(), lines()).is_err());
    // removeCluster removes every occurrence, and the cell when it empties.
    let mut grid = grid;
    grid.add_cluster((1, 1), 4);
    grid.add_cluster((1, 1), 4);
    grid.add_cluster((1, 1), 5);
    grid.remove_cluster((1, 1), 4);
    assert_eq!(grid.clusters((1, 1)), vec![5]);
    grid.remove_cluster((1, 1), 5);
    assert_eq!(grid.cell_count(), 0);
}

// ---------------------------------------------------------------- GridBasedCluster

#[test]
fn grid_based_cluster_class_test() {
    let position = DPosition::xy(4.5, 5.5);
    let bounding_box = DIntervalBase::new(position, position);
    let points = vec![1, 6, 2];
    let prop_a = 1;
    let prop_b = vec![1, 2, 3];

    let cluster = GridBasedCluster::new(
        position,
        bounding_box,
        points.clone(),
        prop_a,
        prop_b.clone(),
    );
    assert_eq!(cluster.centre().x(), 4.5);
    let cluster3 = GridBasedCluster::without_properties(position, bounding_box, points.clone());
    assert_eq!(cluster3.centre().x(), 4.5);
    // Source three-argument constructor: property A -1, one -1 per point.
    assert_eq!(cluster3.property_a(), -1);
    assert_eq!(cluster3.properties_b(), &[-1, -1, -1]);

    assert_eq!(cluster.centre().x(), 4.5);
    assert_eq!(cluster.centre().y(), 5.5);
    assert_eq!(cluster.bounding_box().min_x(), 4.5);
    assert_eq!(cluster.bounding_box().max_y(), 5.5);
    assert_eq!(cluster.points()[0], 1);
    assert_eq!(cluster.points()[2], 2);
    assert_eq!(cluster.property_a(), 1);
    assert_eq!(cluster.properties_b()[0], 1);
    assert_eq!(cluster.properties_b()[2], 3);

    let cluster1 = GridBasedCluster::new(
        DPosition::xy(4.5, 5.5),
        bounding_box,
        points.clone(),
        prop_a,
        prop_b.clone(),
    );
    let cluster2 = GridBasedCluster::new(
        DPosition::xy(4.5, 6.5),
        bounding_box,
        points,
        prop_a,
        prop_b,
    );
    assert!(cluster1 < cluster2);
    assert!(cluster2 > cluster1);
    assert!(cluster1 == cluster1);
}

// ---------------------------------------------------------------- GridBasedClustering

/// `MultiplexClustering::MultiplexDistance(rt_scaling)`
/// (`FEATUREFINDER/MultiplexClustering.cpp:184-187` at the pin), which the
/// class test uses as its metric: `sqrt(dx*dx + s*s*dy*dy)`, evaluated left to
/// right as in the source.
fn multiplex_distance(rt_scaling: f64) -> impl Fn(&DPosition<2>, &DPosition<2>) -> f64 + Clone {
    move |p1, p2| {
        let dx = p1.x() - p2.x();
        let dy = p1.y() - p2.y();
        (dx * dx + rt_scaling * rt_scaling * dy * dy).sqrt()
    }
}

struct Data {
    x: Vec<f64>,
    y: Vec<f64>,
    a: Vec<i32>,
    b: Vec<i32>,
}

fn class_test_data() -> Data {
    let mut data = Data {
        x: Vec::new(),
        y: Vec::new(),
        a: Vec::new(),
        b: Vec::new(),
    };
    for i in 0..1000 {
        data.x.push(5.0 * (f64::from(i).sin() + 1.0));
        data.y.push(5.0 * (f64::from(i + 18).sin() + 1.0));
        data.a.push(1); // the same within each cluster
        data.b.push(i); // different within each cluster
    }
    data
}

#[test]
fn grid_based_clustering_class_test() {
    let metric = multiplex_distance(1.0);
    let d = class_test_data();

    // Seven-argument constructor.
    let mut clustering =
        GridBasedClustering::new(metric.clone(), &d.x, &d.y, &d.a, &d.b, lines(), lines()).unwrap();
    clustering.cluster().unwrap();
    assert_eq!(clustering.results().len(), 12);

    // Five-argument constructor.
    let mut clustering =
        GridBasedClustering::without_properties(metric.clone(), &d.x, &d.y, lines(), lines())
            .unwrap();
    clustering.cluster().unwrap();
    assert_eq!(clustering.results().len(), 12);

    // The remaining sections share one object, as in the class test.
    let mut clustering =
        GridBasedClustering::without_properties(metric, &d.x, &d.y, lines(), lines()).unwrap();
    // void cluster()
    clustering.cluster().unwrap();
    assert_eq!(clustering.results().len(), 12);
    // getResults
    clustering.cluster().unwrap();
    assert_eq!(clustering.results().len(), 12);
    // void extendClustersY()
    clustering.cluster().unwrap();
    clustering.extend_clusters_y().unwrap();
    assert_eq!(clustering.results().len(), 11);
    // void removeSmallClustersY(double threshold_y)
    clustering.cluster().unwrap();
    clustering.remove_small_clusters_y(0.9);
    assert_eq!(clustering.results().len(), 8);
}

#[test]
fn grid_based_clustering_invariants() {
    // Invariants of the documented contract, not class-test literals: every
    // point ends in exactly one final cluster, and properties B never repeat
    // inside a cluster when properties A and B are set.
    let d = class_test_data();
    let mut clustering = GridBasedClustering::new(
        multiplex_distance(1.0),
        &d.x,
        &d.y,
        &d.a,
        &d.b,
        lines(),
        lines(),
    )
    .unwrap();
    clustering.cluster().unwrap();
    let mut seen = vec![false; d.x.len()];
    for cluster in clustering.results().values() {
        for &point in cluster.points() {
            let point = usize::try_from(point).unwrap();
            assert!(!seen[point], "point {point} in two clusters");
            seen[point] = true;
        }
        let mut b = cluster.properties_b().to_vec();
        b.sort_unstable();
        b.dedup();
        assert_eq!(b.len(), cluster.properties_b().len());
    }
    assert!(seen.iter().all(|&s| s));

    // With all properties B equal every merge is vetoed: 1000 singletons.
    let same_b = vec![7; d.x.len()];
    let mut vetoed = GridBasedClustering::new(
        multiplex_distance(1.0),
        &d.x,
        &d.y,
        &d.a,
        &same_b,
        lines(),
        lines(),
    )
    .unwrap();
    vetoed.cluster().unwrap();
    assert_eq!(vetoed.results().len(), d.x.len());
}

#[test]
fn grid_based_clustering_refuses_bad_input() {
    let metric = multiplex_distance(1.0);
    // Mismatched lengths: the source reads past the shorter vector.
    assert!(
        GridBasedClustering::new(
            metric.clone(),
            &[1.0, 2.0],
            &[1.0],
            &[-1],
            &[-1],
            lines(),
            lines()
        )
        .is_err()
    );
    // A point outside the grid: source `Exception::IllegalArgument`.
    assert!(
        GridBasedClustering::without_properties(metric, &[11.0], &[1.0], lines(), lines()).is_err()
    );
    // NaN from the metric would break the source's multiset ordering.
    let nan_metric = |_: &DPosition<2>, _: &DPosition<2>| f64::NAN;
    assert!(
        GridBasedClustering::without_properties(
            nan_metric,
            &[1.0, 1.5],
            &[1.0, 1.5],
            lines(),
            lines()
        )
        .is_err()
    );
}

#[test]
fn minimum_distance_compares_by_distance() {
    let a = MinimumDistance::new(1, 2, 0.5);
    let b = MinimumDistance::new(3, 4, 0.5);
    let c = MinimumDistance::new(1, 2, 0.75);
    assert_eq!((a.cluster_index(), a.nearest_neighbour_index()), (1, 2));
    assert!(a == b);
    assert!(a < c);
    assert!(c > a);
}

// ---------------------------------------------------------------- HashGrid

#[test]
fn hash_grid_class_test() {
    let cell_dimension = DPosition::xy(1.0, 1.0);

    // HashGrid(const ClusterCenter &cell_dimension)
    let t = HashGrid::<()>::new(cell_dimension);
    assert_eq!(*t.cell_dimension(), cell_dimension);
    assert_eq!(t.grid_dimension()[0], 0);
    assert_eq!(t.grid_dimension()[1], 0);

    // cell_iterator insert(const value_type &v)
    let mut t = HashGrid::new(cell_dimension);
    let key1 = DPosition::xy(1.0, 2.0);
    let (cell, position) = t.insert(key1, ()).unwrap();
    assert_eq!(t.grid_dimension()[0] as f64, key1[0]);
    assert_eq!(t.grid_dimension()[1] as f64, key1[1]);
    let inserted = &t.grid_at(&cell).unwrap()[position].0;
    assert_eq!((inserted[0], inserted[1]), (key1[0], key1[1]));
    let key2 = DPosition::xy(2.0, 3.0);
    let (cell, position) = t.insert(key2, ()).unwrap();
    assert_eq!(t.grid_dimension()[0] as f64, key2[0]);
    assert_eq!(t.grid_dimension()[1] as f64, key2[1]);
    let inserted = &t.grid_at(&cell).unwrap()[position].0;
    assert_eq!((inserted[0], inserted[1]), (key2[0], key2[1]));
    let key = DPosition::xy(0.0, i64::MIN as f64 - 1e5);
    assert!(t.insert(key, ()).is_err());
    let key = DPosition::xy(0.0, i64::MAX as f64 + 1e5);
    assert!(t.insert(key, ()).is_err());

    // void erase(iterator pos)
    let mut t = HashGrid::new(cell_dimension);
    t.insert(DPosition::xy(0.0, 0.0), ()).unwrap();
    assert_eq!(t.len(), 1);
    let (cell, _) = t.iter().next().map(|(cell, pair)| (cell, pair.0)).unwrap();
    t.erase_at(cell, 0).unwrap();
    assert_eq!(t.len(), 0);

    // size_type erase(const key_type& key)
    let mut t = HashGrid::new(cell_dimension);
    let key = DPosition::xy(1.0, 2.0);
    t.insert(key, ()).unwrap();
    assert_eq!(t.erase(&key).unwrap(), 1);
    t.insert(key, ()).unwrap();
    t.insert(key, ()).unwrap();
    assert_eq!(t.erase(&key).unwrap(), 2);

    // void clear()
    let mut t = HashGrid::new(cell_dimension);
    t.insert(DPosition::xy(1.0, 2.0), ()).unwrap();
    assert!(!t.is_empty());
    t.clear();
    assert!(t.is_empty());

    // iterator begin() / const_iterator begin() const
    let mut t = HashGrid::new(cell_dimension);
    t.insert(DPosition::xy(0.0, 0.0), ()).unwrap();
    t.insert(DPosition::xy(1.0, 0.0), ()).unwrap();
    t.insert(DPosition::xy(2.0, 0.0), ()).unwrap();
    let mut it = t.iter();
    for _ in 0..3 {
        let (_, (center, _)) = it.next().unwrap();
        assert!(center[0] <= 2.0);
        assert_eq!(center[1], 0.0);
    }
    assert!(it.next().is_none());
    let mut it = t.iter_mut();
    for _ in 0..3 {
        let (_, center, _) = it.next().unwrap();
        assert!(center[0] <= 2.0);
        assert_eq!(center[1], 0.0);
    }
    assert!(it.next().is_none());

    // iterator end() / const_iterator end() const
    let t = HashGrid::<()>::new(cell_dimension);
    assert!(t.iter().next().is_none());

    // bool empty() const
    let mut t = HashGrid::new(cell_dimension);
    assert!(t.is_empty());
    t.insert(DPosition::xy(0.0, 0.0), ()).unwrap();
    assert!(!t.is_empty());
    t.insert(DPosition::xy(0.0, 0.0), ()).unwrap();
    assert!(!t.is_empty());

    // size_type size() const
    let mut t = HashGrid::new(cell_dimension);
    assert_eq!(t.len(), 0);
    t.insert(DPosition::xy(0.0, 0.0), ()).unwrap();
    assert_eq!(t.len(), 1);
    t.insert(DPosition::xy(0.0, 0.0), ()).unwrap();
    assert_eq!(t.len(), 2);
    t.insert(DPosition::xy(1.0, 0.0), ()).unwrap();
    assert_eq!(t.len(), 3);

    // grid_begin() const / grid_begin()
    let mut t = HashGrid::new(cell_dimension);
    t.insert(DPosition::xy(1.0, 2.0), ()).unwrap();
    assert_eq!(t.grid_iter().next().unwrap().1.len(), 1);
    assert_eq!(t.grid_iter_mut().next().unwrap().1.len(), 1);

    // grid_end() const / grid_end()
    let mut t = HashGrid::<()>::new(cell_dimension);
    assert!(t.grid_iter().next().is_none());
    assert!(t.grid_iter_mut().next().is_none());

    // grid_at(const CellIndex &x), const and non-const: std::out_of_range.
    let mut t = HashGrid::<()>::new(cell_dimension);
    assert!(t.grid_at(&[0, 0]).is_err());
    assert!(t.grid_at_mut(&[0, 0]).is_err());

    // [EXTRA] std::size_t hash_value(const DPosition<N, T> &b)
    assert_eq!(hash_value(&[1]), hash_value(&[1]));
    assert_ne!(hash_value(&[1]), hash_value(&[2]));
    assert_eq!(hash_value(&[1, 1]), hash_value(&[1, 1]));
    assert_ne!(hash_value(&[1, 1]), hash_value(&[1, 2]));
}

#[test]
fn hash_grid_cell_index_boundaries() {
    let t = HashGrid::<()>::new(DPosition::xy(0.5, 2.0));
    // floor(key / cell_dimension) per dimension.
    assert_eq!(
        t.cell_index_at_cluster_center(&DPosition::xy(1.2, -0.1))
            .unwrap(),
        [2, -1]
    );
    // NaN and 2^63 are refused (undefined conversions in the source).
    assert!(
        t.cell_index_at_cluster_center(&DPosition::xy(f64::NAN, 0.0))
            .is_err()
    );
    assert!(
        t.cell_index_at_cluster_center(&DPosition::xy(0.0, 2.0 * 9_223_372_036_854_775_808.0))
            .is_err()
    );
    // grid_dimension never goes below zero and survives clear().
    let mut t = HashGrid::new(DPosition::xy(1.0, 1.0));
    t.insert(DPosition::xy(-3.0, 4.0), 'a').unwrap();
    assert_eq!(*t.grid_dimension(), [0, 4]);
    t.clear();
    assert_eq!(*t.grid_dimension(), [0, 4]);
    // erase_at an invalid location is an error, not undefined behaviour.
    assert!(t.erase_at([0, 0], 0).is_err());
}
