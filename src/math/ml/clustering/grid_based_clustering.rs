// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! 2D hierarchical clustering on a grid, for many small clusters.
//!
//! Native equivalent of `ML/CLUSTERING/GridBasedClustering.h` (the class
//! template) and `GridBasedClustering.cpp` (`MinimumDistance`). The source's
//! `std::multiset<MinimumDistance>` orders by distance alone and inserts an
//! element after its equals; here it is an ordered map keyed by distance and
//! an insertion sequence number, which is the same order. Its two hash-map
//! indexes are only ever searched, never iterated in an observable way, so
//! ordered maps replace them without changing a result. See
//! `docs/CLUSTERING_SUPPORT.md`.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

use super::clustering_grid::{CellIndex, ClusteringGrid};
use super::grid_based_cluster::{GridBasedCluster, Point, Rectangle, enlarge};
use crate::concept::progress_logger::ProgressReporter;
use crate::{Error, Result};

/// Basic data structure for distances between clusters (source
/// `MinimumDistance`).
///
/// **Comparison is by distance only**, as in the source, where the operators
/// exist to order a `std::multiset`: two entries with equal distances are
/// equal whatever their cluster indices.
#[derive(Clone, Copy, Debug)]
pub struct MinimumDistance {
    /// Index in the cluster list (source `cluster_index_`).
    cluster_index: i32,
    /// Index of the nearest neighbour of that cluster (source
    /// `nearest_neighbour_index_`).
    nearest_neighbour_index: i32,
    /// Distance between the cluster and its nearest neighbour (source
    /// `distance_`).
    distance: f64,
}

impl MinimumDistance {
    /// A distance entry (source constructor; the source's default constructor
    /// is private and not ported).
    pub const fn new(cluster_index: i32, nearest_neighbour_index: i32, distance: f64) -> Self {
        Self {
            cluster_index,
            nearest_neighbour_index,
            distance,
        }
    }

    /// The cluster index (source `getClusterIndex`).
    pub const fn cluster_index(&self) -> i32 {
        self.cluster_index
    }

    /// The index of the nearest cluster (source `getNearestNeighbourIndex`).
    pub const fn nearest_neighbour_index(&self) -> i32 {
        self.nearest_neighbour_index
    }
}

impl PartialEq for MinimumDistance {
    /// Source `operator==`: equal distances.
    fn eq(&self, other: &Self) -> bool {
        self.distance == other.distance
    }
}

impl PartialOrd for MinimumDistance {
    /// Source `operator<` / `operator>`: by distance.
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.distance.partial_cmp(&other.distance)
    }
}

/// A distance known not to be NaN, ordered as `<` orders it (so `-0.0` and
/// `0.0` are equivalent, as in the source's multiset).
#[derive(Clone, Copy, Debug)]
struct Distance(f64);

impl PartialEq for Distance {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for Distance {}
impl PartialOrd for Distance {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Distance {
    fn cmp(&self, other: &Self) -> Ordering {
        // NaN never enters the multiset; see `find_nearest_neighbour`.
        self.0.partial_cmp(&other.0).unwrap_or(Ordering::Equal)
    }
}

/// Position of one element of the emulated multiset: its distance, then its
/// insertion number, so equal distances keep insertion order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Slot {
    distance: Distance,
    sequence: u64,
}

/// 2D hierarchical clustering optimized for large data sets containing many
/// small clusters, that is, cluster dimensions much smaller than the
/// dimension of the entire data set (source `GridBasedClustering<Metric>`).
///
/// The clustering problem therefore simplifies to a number of local
/// clustering problems, each solved on a few adjacent cells of a larger grid.
/// The grid spacing is determined by the expected typical cluster size in each
/// region.
///
/// Each data point can have two additional properties A and B: in each cluster
/// all properties A must be the same and all properties B different (`-1`
/// means not set, which disables the corresponding rule).
///
/// `M` is the metric measuring the distance between points in the plane
/// (source template parameter `Metric`, called as `metric(p1, p2)`), for
/// example `MultiplexClustering::MultiplexDistance`. The source also derives
/// from `ProgressLogger`; progress goes to the reporter passed to
/// [`cluster_with_progress`](Self::cluster_with_progress).
#[derive(Clone, Debug)]
pub struct GridBasedClustering<M> {
    /// Metric for the distance between points (source `metric_`).
    metric: M,
    /// Grid on which the cluster positions are registered (source `grid_`).
    grid: ClusteringGrid,
    /// Clusters still being merged, by cluster index (source `clusters_`).
    clusters: BTreeMap<i32, GridBasedCluster>,
    /// Clusters no longer merged, by cluster index (source `clusters_final_`).
    clusters_final: BTreeMap<i32, GridBasedCluster>,
    /// The minimum distances, smallest first (source `distances_`).
    distances: BTreeMap<Slot, MinimumDistance>,
    /// Reverse nearest-neighbour lookup: nearest neighbour index to the slots
    /// naming it (source `reverse_nns_`).
    reverse_nns: BTreeMap<i32, Vec<Slot>>,
    /// Cluster index to its slot (source `distance_it_for_cluster_idx_`).
    distance_for_cluster: BTreeMap<i32, Slot>,
    /// Next multiset insertion number.
    next_sequence: u64,
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(format!("GridBasedClustering: {message}"))
}

impl<M: Fn(&Point, &Point) -> f64> GridBasedClustering<M> {
    /// Initializes all data structures for clustering the points
    /// `(data_x[i], data_y[i])` with properties (source seven-argument
    /// constructor).
    ///
    /// `properties_a` holds property A of each point (the same within each
    /// cluster), `properties_b` property B (different within each cluster).
    /// `grid_spacing_x` and `grid_spacing_y` are the grid lines, see
    /// [`ClusteringGrid::new`]. Every point starts as its own cluster; each
    /// cluster without a mergeable neighbour in its own or the eight
    /// surrounding cells is final right away.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the four data vectors differ in length
    /// (the source reads past the shorter ones), when there are more points
    /// than an `int` cluster index can address, when a point lies outside the
    /// grid (source `Exception::IllegalArgument` from `getIndex`), when the
    /// metric yields NaN for a nearest neighbour (which breaks the source's
    /// multiset ordering), and the errors of [`ClusteringGrid::new`].
    pub fn new(
        metric: M,
        data_x: &[f64],
        data_y: &[f64],
        properties_a: &[i32],
        properties_b: &[i32],
        grid_spacing_x: Vec<f64>,
        grid_spacing_y: Vec<f64>,
    ) -> Result<Self> {
        let n = data_x.len();
        if data_y.len() != n || properties_a.len() != n || properties_b.len() != n {
            return Err(invalid("data and property vectors differ in length"));
        }
        if i32::try_from(n).is_err() {
            return Err(invalid("more points than an int cluster index can address"));
        }
        let mut clustering = Self {
            metric,
            grid: ClusteringGrid::new(grid_spacing_x, grid_spacing_y)?,
            clusters: BTreeMap::new(),
            clusters_final: BTreeMap::new(),
            distances: BTreeMap::new(),
            reverse_nns: BTreeMap::new(),
            distance_for_cluster: BTreeMap::new(),
            next_sequence: 0,
        };
        clustering.init(data_x, data_y, properties_a, properties_b)?;
        Ok(clustering)
    }

    /// As [`new`](Self::new) with every property A and B set to `-1`, so the
    /// properties are ignored (source five-argument constructor).
    ///
    /// # Errors
    ///
    /// As [`new`](Self::new).
    pub fn without_properties(
        metric: M,
        data_x: &[f64],
        data_y: &[f64],
        grid_spacing_x: Vec<f64>,
        grid_spacing_y: Vec<f64>,
    ) -> Result<Self> {
        let unset = vec![-1; data_x.len()];
        Self::new(
            metric,
            data_x,
            data_y,
            &unset,
            &unset,
            grid_spacing_x,
            grid_spacing_y,
        )
    }

    /// Performs the hierarchical clustering (source `cluster`): repeatedly
    /// merges the two closest clusters until no cluster has a mergeable
    /// neighbour in its own or the eight surrounding grid cells, that is,
    /// until the clusters' dimensions reach that of a cell.
    ///
    /// A merged cluster's centre is the point-count weighted mean of the two
    /// centres, its bounding box the union of both, its points those of the
    /// first cluster followed by those of the second, and it keeps the first
    /// cluster's index. Calling this again after it finished does nothing.
    ///
    /// # Errors
    ///
    /// As [`cluster_with_progress`](Self::cluster_with_progress).
    pub fn cluster(&mut self) -> Result<()> {
        self.cluster_with_progress(&mut ProgressReporter::silent())
    }

    /// [`cluster`](Self::cluster), reporting progress to `progress` with the
    /// label "clustering".
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when two merged clusters differ in property A
    /// (source `Exception::InvalidValue`, "Property A of both clusters not the
    /// same. "; the merge veto makes this unreachable), when a merged centre
    /// falls outside the grid (source `Exception::IllegalArgument`), when the
    /// metric yields NaN for a nearest neighbour, and where the source would
    /// dereference a missing entry of its own bookkeeping (not reachable with a
    /// non-negative metric); plus the errors of `progress`. On error the
    /// clustering is left part-way.
    pub fn cluster_with_progress(&mut self, progress: &mut ProgressReporter<'_>) -> Result<()> {
        let clusters_start = self.clusters.len();
        progress.start_count(clusters_start, "clustering")?;

        let zero = Slot {
            distance: Distance(0.0),
            sequence: 0,
        };
        while !self.clusters.is_empty() {
            progress.set_count(clusters_start - self.clusters.len())?;

            let (&smallest, &minimum) = self
                .distances
                .range(zero..)
                .next()
                .ok_or_else(|| invalid("no non-negative minimum distance left"))?;
            let index1 = minimum.cluster_index;
            let index2 = minimum.nearest_neighbour_index;
            self.erase_min_distance(smallest)?;

            let cluster1 = self.cluster_at(index1)?.clone();
            let cluster2 = self.cluster_at(index2)?.clone();
            let points1 = cluster1.points();
            let points2 = cluster2.points();
            let mut new_points = Vec::with_capacity(points1.len() + points2.len());
            new_points.extend_from_slice(points1);
            new_points.extend_from_slice(points2);

            let (n1, n2) = (points1.len() as f64, points2.len() as f64);
            let total = (points1.len() + points2.len()) as f64;
            let new_x = (cluster1.centre().x() * n1 + cluster2.centre().x() * n2) / total;
            let new_y = (cluster1.centre().y() * n1 + cluster2.centre().y() * n2) / total;

            let cell1 = self.grid.index(cluster1.centre())?;
            let cell2 = self.grid.index(cluster2.centre())?;
            let new_cell = self.grid.index(&Point::xy(new_x, new_y))?;
            self.grid.remove_cluster(cell1, index1);
            self.grid.remove_cluster(cell2, index2);
            self.grid.add_cluster(new_cell, index1);

            let mut new_box: Rectangle = *cluster1.bounding_box();
            enlarge(&mut new_box, cluster2.bounding_box().min_position());
            enlarge(&mut new_box, cluster2.bounding_box().max_position());

            if cluster1.property_a() != cluster2.property_a() {
                return Err(Error::InvalidValue(
                    "Property A of both clusters not the same. ".into(),
                ));
            }
            let new_a = cluster1.property_a();
            let mut new_b =
                Vec::with_capacity(cluster1.properties_b().len() + cluster2.properties_b().len());
            new_b.extend_from_slice(cluster1.properties_b());
            new_b.extend_from_slice(cluster2.properties_b());

            let new_cluster =
                GridBasedCluster::new(Point::xy(new_x, new_y), new_box, new_points, new_a, new_b);
            self.clusters.remove(&index1);
            self.clusters.remove(&index2);
            self.clusters.insert(index1, new_cluster);

            let mut to_update = BTreeSet::from([index1]);

            // The second cluster no longer exists: drop its distance without
            // an update. The first one's was dropped above.
            let slot2 = *self
                .distance_for_cluster
                .get(&index2)
                .ok_or_else(|| invalid("merged cluster has no distance entry"))?;
            self.erase_min_distance(slot2)?;

            for index in [index1, index2] {
                let slots = self.reverse_nns.get(&index).cloned().unwrap_or_default();
                for slot in slots {
                    let entry = self
                        .distances
                        .get(&slot)
                        .ok_or_else(|| invalid("stale nearest-neighbour entry"))?;
                    to_update.insert(entry.cluster_index);
                    self.erase_min_distance(slot)?;
                }
            }

            for index in to_update {
                let cluster = self.cluster_at(index)?.clone();
                if self.find_nearest_neighbour(&cluster, index)? {
                    let cell = self.grid.index(cluster.centre())?;
                    self.grid.remove_cluster(cell, index);
                    self.clusters.remove(&index);
                }
            }
        }

        progress.end()
    }

    /// Extends clusters in y, merging adjacent clusters further so that a
    /// cluster can now span several cells (source `extendClustersY`).
    ///
    /// The final clusters are registered on a grid with the original x lines
    /// and a single cell in y. In each x cell holding more than one cluster,
    /// the clusters are ordered by centre y (stable), and each adjacent pair
    /// whose x ranges overlap (one box's minimum or maximum x lies within the
    /// other's) is merged into a cluster without properties that keeps the
    /// upper cluster's index.
    ///
    /// The source compares each pair from the list it sorted before the loop,
    /// so after a merge the next pair still sees the unmerged upper cluster,
    /// and it looks clusters up by centre y through a `std::map` keyed on the
    /// y ordering. Both behaviours are reproduced. Only the y cell with index 1
    /// is scanned, so a cluster whose centre lies exactly on the lowest y line
    /// (cell 0) is never extended, as in the source.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when two clusters of one x cell share the same
    /// centre y and are merged, where the source's lookup by y yields one index
    /// for both and then erases the end of its map (undefined behaviour); and
    /// the errors of [`ClusteringGrid::new`] and [`ClusteringGrid::index`]. The
    /// source's only consumer, `MultiplexClustering`, has its call commented
    /// out (`MultiplexClustering.cpp:165`). On error the results are left
    /// part-way.
    pub fn extend_clusters_y(&mut self) -> Result<()> {
        let grid_spacing_x = self.grid.grid_spacing_x().to_vec();
        let grid_spacing_y = self.grid.grid_spacing_y();
        let grid_spacing_y_new = match (grid_spacing_y.first(), grid_spacing_y.last()) {
            (Some(&first), Some(&last)) => vec![first, last],
            _ => return Err(invalid("empty y grid")),
        };
        let mut grid_x_only = ClusteringGrid::new(grid_spacing_x.clone(), grid_spacing_y_new)?;

        for (&index, cluster) in &self.clusters_final {
            grid_x_only.add_cluster(grid_x_only.index(cluster.centre())?, index);
        }

        for cell in 0..grid_spacing_x.len() {
            // The length fits `i32`, checked by `ClusteringGrid::new`.
            let grid_index: CellIndex = (cell as i32, 1);
            if !grid_x_only.is_non_empty_cell(grid_index) {
                continue;
            }
            let cluster_indices = grid_x_only.clusters(grid_index);
            if cluster_indices.len() <= 1 {
                continue;
            }
            let mut cluster_list = Vec::with_capacity(cluster_indices.len());
            // `std::map<GridBasedCluster, int>`: keys equivalent under the
            // y ordering collapse onto the first inserted index.
            let mut index_list: Vec<(f64, i32)> = Vec::new();
            for &index in &cluster_indices {
                let cluster = self
                    .clusters_final
                    .get(&index)
                    .ok_or_else(|| invalid("registered cluster missing"))?
                    .clone();
                let y = cluster.centre().y();
                if !index_list.iter().any(|&(key, _)| equivalent(key, y)) {
                    index_list.push((y, index));
                }
                cluster_list.push(cluster);
            }
            cluster_list.sort_by(|a, b| {
                let (ya, yb) = (a.centre().y(), b.centre().y());
                if ya < yb {
                    Ordering::Less
                } else if yb < ya {
                    Ordering::Greater
                } else {
                    Ordering::Equal
                }
            });
            let lookup = |cluster: &GridBasedCluster| -> Result<i32> {
                let y = cluster.centre().y();
                index_list
                    .iter()
                    .find(|&&(key, _)| equivalent(key, y))
                    .map(|&(_, index)| index)
                    .ok_or_else(|| invalid("cluster missing from the y index"))
            };

            for pair in 1..cluster_list.len() {
                let c1 = &cluster_list[pair];
                let c2 = &cluster_list[pair - 1];
                let (centre1x, centre1y) = (c1.centre().x(), c1.centre().y());
                let (centre2x, centre2y) = (c2.centre().x(), c2.centre().y());
                let box1 = c1.bounding_box();
                let box2 = c2.bounding_box();
                let (box1x_min, box1x_max) = (box1.min_x(), box1.max_x());
                let (box1y_min, box1y_max) = (box1.min_y(), box1.max_y());
                let (box2x_min, box2x_max) = (box2.min_x(), box2.max_x());
                let (box2y_min, box2y_max) = (box2.min_y(), box2.max_y());

                let overlap = (box1x_min <= box2x_max && box1x_min >= box2x_min)
                    || (box1x_max >= box2x_min && box1x_max <= box2x_max);
                if !overlap {
                    continue;
                }
                let points1 = c1.points();
                let points2 = c2.points();
                let mut new_points = Vec::with_capacity(points1.len() + points2.len());
                new_points.extend_from_slice(points1);
                new_points.extend_from_slice(points2);
                let (n1, n2) = (points1.len() as f64, points2.len() as f64);
                let total = (points1.len() + points2.len()) as f64;
                let new_x = (centre1x * n1 + centre2x * n2) / total;
                let new_y = (centre1y * n1 + centre2y * n2) / total;
                let min_x = cpp_min(box1x_min, box2x_min);
                let min_y = cpp_min(box1y_min, box2y_min);
                let max_x = cpp_max(box1x_max, box2x_max);
                let max_y = cpp_max(box1y_max, box2y_max);
                let new_centre = Point::xy(new_x, new_y);
                let new_bounding_box =
                    Rectangle::new(Point::xy(min_x, min_y), Point::xy(max_x, max_y));
                let new_cluster =
                    GridBasedCluster::without_properties(new_centre, new_bounding_box, new_points);

                let index1 = lookup(c1)?;
                let index2 = lookup(c2)?;
                if self.clusters_final.remove(&index1).is_none()
                    || self.clusters_final.remove(&index2).is_none()
                {
                    return Err(invalid(
                        "two merged clusters share a centre y; the source erases map::end() here",
                    ));
                }
                self.clusters_final.entry(index1).or_insert(new_cluster);

                let cell1 = grid_x_only.index(c1.centre())?;
                let cell2 = grid_x_only.index(c2.centre())?;
                let new_cell = grid_x_only.index(&new_centre)?;
                grid_x_only.remove_cluster(cell1, index1);
                grid_x_only.remove_cluster(cell2, index2);
                grid_x_only.add_cluster(new_cell, index1);
            }
        }
        Ok(())
    }

    /// Removes the final clusters whose bounding box spans less than
    /// `threshold_y` in y (source `removeSmallClustersY`).
    #[allow(clippy::neg_cmp_op_on_partial_ord, reason = "NaN semantics of `<`")]
    pub fn remove_small_clusters_y(&mut self, threshold_y: f64) {
        self.clusters_final.retain(|_, cluster| {
            let bounding_box = cluster.bounding_box();
            !(bounding_box.max_y() - bounding_box.min_y() < threshold_y)
        });
    }

    /// The final results, cluster index to cluster (source `getResults`, which
    /// returns a copy).
    pub fn results(&self) -> &BTreeMap<i32, GridBasedCluster> {
        &self.clusters_final
    }

    /// Source `init_`: one cluster per point, registered on the grid, then the
    /// nearest neighbour of each in index order.
    fn init(
        &mut self,
        data_x: &[f64],
        data_y: &[f64],
        properties_a: &[i32],
        properties_b: &[i32],
    ) -> Result<()> {
        for i in 0..data_x.len() {
            // `i` fits `i32`, checked by `new`.
            let index = i as i32;
            let position = Point::xy(data_x[i], data_y[i]);
            let bounding_box = Rectangle::new(position, position);
            let cluster = GridBasedCluster::new(
                position,
                bounding_box,
                vec![index],
                properties_a[i],
                vec![properties_b[i]],
            );
            self.clusters.entry(index).or_insert(cluster);
            let cell = self.grid.index(&position)?;
            self.grid.add_cluster(cell, index);
        }

        let indices: Vec<i32> = self.clusters.keys().copied().collect();
        for index in indices {
            let cluster = self.cluster_at(index)?.clone();
            if self.find_nearest_neighbour(&cluster, index)? {
                let cell = self.grid.index(cluster.centre())?;
                self.grid.remove_cluster(cell, index);
                self.clusters.remove(&index);
            }
        }
        Ok(())
    }

    /// Source `mergeVeto_`: whether `c1` and `c2` must **not** be merged.
    ///
    /// The source documents the return value the other way round ("true ->
    /// clusters can be merged"); its code, and `findNearestNeighbour_`, treat
    /// `true` as a veto, and so does this port. There is no veto when either
    /// property A is unset (`-1`); a veto when the properties A differ; no
    /// veto when any property B of either cluster is unset; otherwise a veto
    /// exactly when the two clusters share a property B.
    fn merge_veto(c1: &GridBasedCluster, c2: &GridBasedCluster) -> bool {
        let (a1, a2) = (c1.property_a(), c2.property_a());
        if a1 == -1 || a2 == -1 {
            return false;
        }
        if a1 != a2 {
            return true;
        }
        let (b1, b2) = (c1.properties_b(), c2.properties_b());
        if b1.contains(&-1) || b2.contains(&-1) {
            return false;
        }
        let b1: BTreeSet<i32> = b1.iter().copied().collect();
        b2.iter().any(|b| b1.contains(b))
    }

    /// Source `findNearestNeighbour_`: records the nearest mergeable cluster
    /// in the cell of `cluster` and its eight neighbours, scanning cells row
    /// by row from `(i-1, j-1)` and each cell's clusters in insertion order;
    /// a later candidate replaces the current one only when strictly closer.
    ///
    /// Returns `true` when there is no such cluster; the cluster is then copied
    /// to the final results and the caller removes it from the working list
    /// and the grid.
    fn find_nearest_neighbour(
        &mut self,
        cluster: &GridBasedCluster,
        cluster_index: i32,
    ) -> Result<bool> {
        let centre = *cluster.centre();
        let cell_index = self.grid.index(&centre)?;
        let mut min_dist = 0.0_f64;
        let mut nearest_neighbour = -1;

        for i in -1..=1_i32 {
            for j in -1..=1_i32 {
                let (Some(x), Some(y)) = (cell_index.0.checked_add(i), cell_index.1.checked_add(j))
                else {
                    continue;
                };
                if !self.grid.is_non_empty_cell((x, y)) {
                    continue;
                }
                for index2 in self.grid.clusters((x, y)) {
                    if index2 == cluster_index {
                        continue;
                    }
                    let cluster2 = self.cluster_at(index2)?;
                    let distance = (self.metric)(&centre, cluster2.centre());
                    if (distance < min_dist || nearest_neighbour == -1)
                        && !Self::merge_veto(cluster, cluster2)
                    {
                        min_dist = distance;
                        nearest_neighbour = index2;
                    }
                }
            }
        }

        if nearest_neighbour == -1 {
            let own = self.cluster_at(cluster_index)?.clone();
            self.clusters_final.entry(cluster_index).or_insert(own);
            return Ok(true);
        }
        if min_dist.is_nan() {
            return Err(invalid("the metric returned NaN for a nearest neighbour"));
        }

        let slot = Slot {
            distance: Distance(min_dist),
            sequence: self.next_sequence,
        };
        self.next_sequence += 1;
        self.distances.insert(
            slot,
            MinimumDistance::new(cluster_index, nearest_neighbour, min_dist),
        );
        self.reverse_nns
            .entry(nearest_neighbour)
            .or_default()
            .push(slot);
        self.distance_for_cluster.insert(cluster_index, slot);
        Ok(false)
    }

    /// Source `eraseMinDistance_`: removes the distance at `slot` and its
    /// reverse-lookup and cluster-lookup entries.
    fn erase_min_distance(&mut self, slot: Slot) -> Result<()> {
        let entry = self
            .distances
            .remove(&slot)
            .ok_or_else(|| invalid("distance entry missing"))?;
        if let Some(slots) = self.reverse_nns.get_mut(&entry.nearest_neighbour_index) {
            if let Some(position) = slots.iter().position(|&s| s == slot) {
                slots.remove(position);
            }
            if slots.is_empty() {
                self.reverse_nns.remove(&entry.nearest_neighbour_index);
            }
        }
        self.distance_for_cluster.remove(&entry.cluster_index);
        Ok(())
    }

    fn cluster_at(&self, index: i32) -> Result<&GridBasedCluster> {
        self.clusters
            .get(&index)
            .ok_or_else(|| invalid("cluster index missing from the cluster list"))
    }
}

/// `!(a < b) && !(b < a)`: equivalence under `<`, as a `std::map` decides it.
// Spelled with negated `<` on purpose: `a == b` differs for NaN.
#[allow(clippy::neg_cmp_op_on_partial_ord)]
fn equivalent(a: f64, b: f64) -> bool {
    !(a < b) && !(b < a)
}

/// `std::min(a, b)`.
fn cpp_min(a: f64, b: f64) -> f64 {
    if b < a { b } else { a }
}

/// `std::max(a, b)`.
fn cpp_max(a: f64, b: f64) -> f64 {
    if a < b { b } else { a }
}
