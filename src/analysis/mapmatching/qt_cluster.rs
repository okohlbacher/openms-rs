// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! A quality-threshold cluster for QT feature grouping.
//!
//! Native equivalent of `DATASTRUCTURES/QTCluster.h`. Like
//! [`GridFeature`](crate::analysis::mapmatching::grid_feature::GridFeature) it lives beside its
//! consumer, `QTClusterFinder`, because it needs peptide sequences. See
//! `docs/CLUSTERING_SUPPORT.md`.
//!
//! The source splits a cluster into a cheap head (`QTCluster`, copied by
//! `boost::heap`) and a `BulkData` body the head points to. Here a cluster owns
//! its body; a consumer that needs cheap handles keeps clusters in a vector and
//! heaps their indices. The source keeps the per-map neighbours in
//! `std::unordered_map`s, whose iteration order decides the order of
//! [`QTCluster::all_neighbors`](crate::analysis::mapmatching::qt_cluster::QTCluster::all_neighbors) and
//! [`QTCluster::elements`](crate::analysis::mapmatching::qt_cluster::QTCluster::elements) and the order in
//! which neighbour distances are summed for the quality (floating-point
//! addition is not associative). The port keeps them in
//! [`StdHashMap`](crate::math::hash_order::StdHashMap), which reproduces
//! libstdc++'s order for `size_t` keys.

use std::collections::{BTreeMap, BTreeSet};

use super::grid_feature::{GridFeature, SourceOrderedSequence};
use crate::chemistry::AASequence;
use crate::math::hash_order::StdHashMap;

/// All neighbouring features of one map with their distances, ascending by
/// distance, equal distances in insertion order (source `NeighborList`, a
/// `std::multimap<double, const GridFeature*>`).
pub type NeighborList<'a> = Vec<(f64, &'a GridFeature<'a>)>;

/// The best neighbour of one map (source `QTCluster::Neighbor`).
#[derive(Clone, Copy, Debug)]
pub struct Neighbor<'a> {
    /// Distance to the cluster centre.
    pub distance: f64,
    /// The neighbouring feature.
    pub feature: &'a GridFeature<'a>,
}

/// One element of a cluster (source `QTCluster::Element`).
#[derive(Clone, Copy, Debug)]
pub struct Element<'a> {
    /// The input map the feature belongs to.
    pub map_index: usize,
    /// The feature.
    pub feature: &'a GridFeature<'a>,
}

/// Elements of a cluster (source `QTCluster::Elements`).
pub type Elements<'a> = Vec<Element<'a>>;

/// The bulk internal data of a cluster: neighbours, annotations and so on
/// (source `QTCluster::BulkData`).
#[derive(Clone, Debug)]
pub struct BulkData<'a> {
    /// The cluster centre (source `center_point_`).
    center_point: &'a GridFeature<'a>,
    /// Unique id of the cluster (source `id_`).
    id: usize,
    /// The best current feature of each map (source `neighbors_`, a
    /// `std::unordered_map<Size, Neighbor>`).
    neighbors: StdHashMap<Neighbor<'a>>,
    /// All neighbours of each map while the cluster is being filled (source
    /// `tmp_neighbors_`, a `std::unordered_map<Size, NeighborList>`).
    tmp_neighbors: StdHashMap<NeighborList<'a>>,
    /// Largest distance a point may have and still belong to the cluster
    /// (source `max_distance_`).
    max_distance: f64,
    /// Number of input maps (source `num_maps_`).
    num_maps: usize,
    /// x coordinate in the grid cell (source `x_coord_`).
    x_coord: i32,
    /// y coordinate in the grid cell (source `y_coord_`).
    y_coord: i32,
    /// The peptide sequences compatible with the centre that give the best
    /// cluster quality (source `annotations_`).
    annotations: BTreeSet<SourceOrderedSequence>,
}

impl<'a> BulkData<'a> {
    /// The body of a cluster around `center_point` (source constructor):
    /// `num_maps` input maps, `max_distance` the largest allowed distance of two
    /// points, `(x_coord, y_coord)` the grid cell and `id` the cluster's unique
    /// id.
    pub fn new(
        center_point: &'a GridFeature<'a>,
        num_maps: usize,
        max_distance: f64,
        x_coord: i32,
        y_coord: i32,
        id: usize,
    ) -> Self {
        Self {
            center_point,
            id,
            neighbors: StdHashMap::new(),
            tmp_neighbors: StdHashMap::new(),
            max_distance,
            num_maps,
            x_coord,
            y_coord,
            annotations: BTreeSet::new(),
        }
    }
}

/// A representation of a QT cluster used for feature grouping (source
/// `QTCluster`).
///
/// A cluster represents a group of corresponding features (or consensus
/// features) from different input maps. It is defined by its centre, one
/// feature, and stores potential elements from other maps with their distances
/// to the centre. Every feature satisfying the constraints with respect to the
/// centre is a potential element, but a feature group holds at most one feature
/// per map, so only the best (closest) one per map is a true element; to save
/// memory only that one is kept.
///
/// QT clustering initially produces all possible, overlapping clusters, then
/// repeatedly extracts the best one and recomputes the rest. Here the
/// clustering is updated instead: when a cluster is finalized, its elements are
/// removed from the remaining clusters, which may change their composition or
/// invalidate them (when their centre is removed).
///
/// The quality of a cluster is the normalized average distance to the centre
/// of present and missing elements, a missing element counting as the
/// user-defined maximum radius.
///
/// Elements are added between [`initialize_cluster`](Self::initialize_cluster)
/// and [`finalize_cluster`](Self::finalize_cluster), and only then.
///
/// With `use_ids`, clusters are extended only with elements sharing at least
/// one peptide sequence with the centre. The quality is then the best over all
/// possible sequences, and that sequence becomes the cluster's only
/// (representative) annotation.
///
/// The source guards several operations with `OPENMS_PRECONDITION`, which
/// only debug builds check; like the Release build, this port does not check
/// them. Its default constructor is deleted, and its copies share one
/// `BulkData`; here `Clone` copies the body.
#[derive(Clone, Debug)]
pub struct QTCluster<'a> {
    /// Quality of the cluster (source `quality_`).
    quality: f64,
    /// The cluster body (source `data_`, a pointer).
    data: BulkData<'a>,
    /// Whether the cluster is valid (source `valid_`).
    valid: bool,
    /// Whether the quality needs recomputing (source `changed_`).
    changed: bool,
    /// Whether peptide IDs are used for matching (source `use_IDs_`).
    use_ids: bool,
    /// Whether all annotations must be collected before the best set of
    /// points can be chosen, which is needed when IDs are used and the centre
    /// does not have exactly one annotation (source `collect_annotations_`).
    collect_annotations: bool,
    /// Whether the cluster accepts no more elements (source `finalized_`).
    finalized: bool,
}

impl<'a> QTCluster<'a> {
    /// A cluster with body `data` (source `QTCluster(BulkData*, bool)`).
    ///
    /// With `use_ids`, the cluster's annotations start as the centre's.
    pub fn new(mut data: BulkData<'a>, use_ids: bool) -> Self {
        if use_ids {
            data.annotations = data.center_point.annotations().clone();
        }
        let collect_annotations = use_ids && data.center_point.annotations().len() != 1;
        Self {
            quality: 0.0,
            data,
            valid: true,
            changed: false,
            use_ids,
            collect_annotations,
            finalized: true,
        }
    }

    /// The cluster centre (source `getCenterPoint`).
    pub fn center_point(&self) -> &'a GridFeature<'a> {
        self.data.center_point
    }

    /// The cluster's id (source `getId`).
    pub fn id(&self) -> usize {
        self.data.id
    }

    /// The retention time of the centre (source `getCenterRT`).
    pub fn center_rt(&self) -> f64 {
        self.data.center_point.rt()
    }

    /// The m/z of the centre (source `getCenterMZ`).
    pub fn center_mz(&self) -> f64 {
        self.data.center_point.mz()
    }

    /// The x coordinate in the grid (source `getXCoord`).
    pub fn x_coord(&self) -> i32 {
        self.data.x_coord
    }

    /// The y coordinate in the grid (source `getYCoord`).
    pub fn y_coord(&self) -> i32 {
        self.data.y_coord
    }

    /// The number of elements, centre included (source `size`).
    pub fn size(&self) -> usize {
        self.data.neighbors.len() + 1
    }

    /// Adds `element` at `distance` from the centre (source `add`).
    ///
    /// Source `@note`: there is no check whether the element already exists
    /// in the cluster. With IDs, an element whose annotations and the
    /// centre's are both non-empty and disjoint is ignored. An element of the
    /// centre's own map is never a neighbour. Otherwise the element becomes
    /// its map's neighbour when that map has none yet or a farther one; while
    /// annotations are collected every element is also kept for
    /// [`quality`](Self::quality)'s annotation choice.
    ///
    /// The source's debug-only preconditions (the cluster is initialized, the
    /// distance does not exceed the maximum) are not checked.
    #[allow(clippy::neg_cmp_op_on_partial_ord, reason = "NaN semantics of `<`")]
    pub fn add(&mut self, element: &'a GridFeature<'a>, distance: f64) {
        let map_index = element.map_index();
        let center_point = self.data.center_point;

        if self.use_ids {
            let one_empty =
                center_point.annotations().is_empty() || element.annotations().is_empty();
            if !one_empty
                && center_point
                    .annotations()
                    .intersection(element.annotations())
                    .next()
                    .is_none()
            {
                return;
            }
        }

        if self.collect_annotations && map_index != center_point.map_index() {
            let list = self
                .data
                .tmp_neighbors
                .get_or_insert_with(map_index, Vec::new);
            // `std::multimap::insert`: after every element not greater.
            let position = list.partition_point(|&(key, _)| !(distance < key));
            list.insert(position, (distance, element));
            self.changed = true;
        }

        if map_index != center_point.map_index() {
            let replace = match self.data.neighbors.get(map_index) {
                None => true,
                Some(current) => distance < current.distance,
            };
            if replace {
                self.data.neighbors.insert(
                    map_index,
                    Neighbor {
                        distance,
                        feature: element,
                    },
                );
                self.changed = true;
            }
        }
    }

    /// The clustered elements: the neighbours in the order of the source's
    /// `std::unordered_map`, then the centre (source `getElements`).
    pub fn elements(&self) -> Elements<'a> {
        let mut elements = self.all_neighbors();
        elements.push(Element {
            map_index: self.data.center_point.map_index(),
            feature: self.data.center_point,
        });
        elements
    }

    /// Updates the cluster after the elements `removed` were taken by another
    /// cluster (source `update`); returns whether the composition changed.
    ///
    /// When the centre is among `removed`, the cluster becomes invalid and
    /// `false` is returned. Otherwise each removed feature that is its map's
    /// current neighbour is dropped. Features are compared by identity (the
    /// source compares pointers). The return value is the "changed" flag, so
    /// a change not yet folded into the quality also returns `true`, as in the
    /// source.
    pub fn update(&mut self, removed: &[Element<'a>]) -> bool {
        if removed
            .iter()
            .any(|element| std::ptr::eq(element.feature, self.data.center_point))
        {
            self.set_invalid();
            return false;
        }
        for element in removed {
            if let Some(current) = self.data.neighbors.get(element.map_index) {
                if std::ptr::eq(current.feature, element.feature) {
                    self.changed = true;
                    self.data.neighbors.remove(element.map_index);
                }
            }
        }
        self.changed
    }

    /// The cluster quality, recomputed first if the cluster changed (source
    /// `getQuality`).
    ///
    /// The quality is `(max_distance - d) / max_distance`, where `d` is the
    /// sum of the neighbour distances plus `max_distance` for every map without
    /// a neighbour, divided by the number of maps other than the centre's.
    /// With IDs and a centre that does not have exactly one annotation, `d`
    /// comes from the best annotation instead (see
    /// [`annotations`](Self::annotations)).
    ///
    /// The neighbour distances are summed in the iteration order of the
    /// source's `std::unordered_map`, as the source does.
    pub fn quality(&mut self) -> f64 {
        if self.changed {
            self.compute_quality();
            self.changed = false;
        }
        self.quality
    }

    /// The cluster quality without recomputing it (source
    /// `getCurrentQuality`).
    pub fn current_quality(&self) -> f64 {
        self.quality
    }

    /// The peptide sequences of the cluster (source `getAnnotations`).
    ///
    /// When the cluster changed, IDs are used, the centre does not have
    /// exactly one annotation and there are neighbours, the best annotation is
    /// chosen first: for each sequence, the best distance of each map (an
    /// unannotated feature counting for every sequence) is summed with
    /// `max_distance` per missing map, and the sequence with the smallest sum
    /// wins, ties going to the first in the source's sequence order (see
    /// [`SourceOrderedSequence`]). The neighbours are then reduced to features
    /// compatible with it. As in the source this needs the neighbours
    /// collected since [`initialize_cluster`](Self::initialize_cluster); after
    /// [`finalize_cluster`](Self::finalize_cluster) they are gone, and a
    /// changed cluster then loses its neighbours here.
    pub fn annotations(&mut self) -> &BTreeSet<SourceOrderedSequence> {
        if self.changed
            && self.use_ids
            && self.data.center_point.annotations().len() != 1
            && !self.data.neighbors.is_empty()
        {
            self.optimize_annotations();
        }
        &self.data.annotations
    }

    /// Marks the cluster invalid and frees its annotations (source
    /// `setInvalid`). Source `@note`: an invalid cluster must not be used
    /// again.
    pub fn set_invalid(&mut self) {
        self.valid = false;
        self.data.annotations.clear();
    }

    /// Whether the cluster is invalid (source `isInvalid`).
    pub fn is_invalid(&self) -> bool {
        !self.valid
    }

    /// Opens the cluster for [`add`](Self::add) (source `initializeCluster`).
    pub fn initialize_cluster(&mut self) {
        self.finalized = false;
        self.data.tmp_neighbors.clear();
    }

    /// Closes the cluster after adding elements (source `finalizeCluster`):
    /// computes the quality if anything changed, which also settles the
    /// annotation, and discards the collected neighbours.
    pub fn finalize_cluster(&mut self) {
        self.quality();
        self.finalized = true;
        self.data.tmp_neighbors.clear();
    }

    /// All current neighbours, in the iteration order of the source's
    /// `std::unordered_map` (source `getAllNeighbors`).
    pub fn all_neighbors(&self) -> Elements<'a> {
        self.data
            .neighbors
            .iter()
            .map(|(map_index, neighbor)| Element {
                map_index,
                feature: neighbor.feature,
            })
            .collect()
    }

    /// Whether the cluster is finalized (source member `finalized_`, which
    /// the source keeps private and only checks in debug builds).
    pub fn is_finalized(&self) -> bool {
        self.finalized
    }

    /// Source `computeQuality_`.
    fn compute_quality(&mut self) {
        let max_distance = self.data.max_distance;
        let num_other = self.data.num_maps.wrapping_sub(1);
        let mut internal_distance = 0.0_f64;
        if !self.use_ids
            || self.data.center_point.annotations().len() == 1
            || self.data.neighbors.is_empty()
        {
            for neighbor in self.data.neighbors.values() {
                internal_distance += neighbor.distance;
            }
            internal_distance +=
                num_other.wrapping_sub(self.data.neighbors.len()) as f64 * max_distance;
        } else {
            internal_distance = self.optimize_annotations();
        }
        internal_distance /= num_other as f64;
        self.quality = (max_distance - internal_distance) / max_distance;
    }

    /// Source `optimizeAnnotations_`: chooses the annotation with the smallest
    /// total distance and returns that total.
    fn optimize_annotations(&mut self) -> f64 {
        let mut seq_table = self.make_seq_table();
        let num_maps = self.data.num_maps;
        let max_distance = self.data.max_distance;

        // Fold the unspecific (unannotated) distances into every sequence.
        let unspecific_key = SourceOrderedSequence(AASequence::default());
        if let Some(unspecific) = seq_table.get(&unspecific_key).cloned() {
            for (sequence, distances) in &mut seq_table {
                if *sequence == unspecific_key {
                    continue;
                }
                for (&map_index, &unspecific_distance) in &unspecific {
                    distances
                        .entry(map_index)
                        .and_modify(|current| *current = cpp_min(*current, unspecific_distance))
                        .or_insert(unspecific_distance);
                }
            }
        }

        let mut best: Option<&SourceOrderedSequence> = seq_table.keys().next();
        let mut best_total = num_maps as f64 * max_distance;
        for (sequence, distances) in &seq_table {
            let missing = num_maps.wrapping_sub(1).wrapping_sub(distances.len());
            let total = distances
                .values()
                .fold(missing as f64 * max_distance, |sum, distance| {
                    sum + distance
                });
            if total < best_total {
                best = Some(sequence);
                best_total = total;
            }
        }
        if let Some(best) = best {
            self.data.annotations = BTreeSet::from([best.clone()]);
        }

        self.recompute_neighbors();
        best_total
    }

    /// Source `recomputeNeighbors_`: for each map, the closest collected
    /// feature compatible with the annotations (sharing one, or having none).
    fn recompute_neighbors(&mut self) {
        self.data.neighbors.clear();
        let data = &mut self.data;
        for (map_index, list) in data.tmp_neighbors.iter() {
            for &(distance, feature) in list {
                let current = feature.annotations();
                if current.is_empty() || current.intersection(&data.annotations).next().is_some() {
                    data.neighbors
                        .insert(map_index, Neighbor { distance, feature });
                    break;
                }
            }
        }
    }

    /// Source `makeSeqTable_`: for each sequence, the best distance per map.
    /// An unannotated feature counts under the empty sequence, and no farther
    /// feature of its map is considered after it.
    fn make_seq_table(&self) -> BTreeMap<SourceOrderedSequence, BTreeMap<usize, f64>> {
        let mut seq_table: BTreeMap<SourceOrderedSequence, BTreeMap<usize, f64>> = BTreeMap::new();
        let mut record = |sequence: &SourceOrderedSequence, map_index: usize, distance: f64| {
            seq_table
                .entry(sequence.clone())
                .or_default()
                .entry(map_index)
                .and_modify(|current| *current = cpp_min(distance, *current))
                .or_insert(distance);
        };
        let unannotated = SourceOrderedSequence(AASequence::default());
        for (map_index, list) in self.data.tmp_neighbors.iter() {
            for &(distance, feature) in list {
                for sequence in feature.annotations() {
                    record(sequence, map_index, distance);
                }
                if feature.annotations().is_empty() {
                    record(&unannotated, map_index, distance);
                    break;
                }
            }
        }
        seq_table
    }
}

impl PartialEq for QTCluster<'_> {
    /// Equal current qualities (the source defines only `operator<`; this
    /// completes the `PartialOrd` it maps to).
    fn eq(&self, other: &Self) -> bool {
        self.quality == other.quality
    }
}

impl PartialOrd for QTCluster<'_> {
    /// Source `operator<`: by current quality, without recomputing it.
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.quality.partial_cmp(&other.quality)
    }
}

/// `std::min(a, b)`: `b` only when `b < a`.
fn cpp_min(a: f64, b: f64) -> f64 {
    if b < a { b } else { a }
}
