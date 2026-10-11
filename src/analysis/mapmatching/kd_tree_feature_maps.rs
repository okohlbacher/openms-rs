// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Features pooled from several maps into one two-dimensional kd-tree on
//! (RT, m/z) for fast range queries.
//!
//! Native equivalent of `ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h` and
//! `ANALYSIS/QUANTITATION/KDTreeFeatureNode.h`. See
//! `docs/FEATURE_GROUPING_SUPPORT.md`.
//!
//! The container borrows the caller's features (the source stores raw
//! pointers; the lifetime `'a` states the source's "the wrapped features must
//! outlive this object"), remembers each feature's map of origin, and keeps
//! its own copy of the retention times so that an alignment transformation
//! can be applied with [`KDTreeFeatureMaps::apply_transformations`](crate::analysis::mapmatching::kd_tree_feature_maps::KDTreeFeatureMaps::apply_transformations) without
//! modifying the features. m/z, intensity and charge are always read through
//! the borrowed feature. Typical consumers are KD-tree based map alignment and
//! feature grouping, and the MS2-to-feature mapping of
//! [`FeatureMapping`](crate::analysis::mapmatching::feature_mapping).
//!
//! # The tree and the order of query results
//!
//! The source's tree is libkdtree++ (`DATASTRUCTURES/KDTree.h`), whose range
//! query reports points in the pre-order of the tree, and the tree's shape
//! depends on the insertion order and on `std::nth_element` (which
//! `optimise()` uses to rebuild it balanced). Consumers break exact distance
//! ties by that order, so this port reproduces the tree exactly: insertion
//! descends left when the new point's coordinate on the level's axis (RT on
//! even levels, m/z on odd ones) is less than the node's and right otherwise,
//! `optimize_tree` collects the points in order and reinserts the median of
//! every range as libkdtree++'s `_M_optimise` does, with a transcription of
//! libstdc++'s introselect `nth_element`, and queries walk the tree in
//! pre-order with libkdtree++'s pruning. The coordinates are read live, as the
//! source's node accessor reads them, so a query after
//! [`KDTreeFeatureMaps::apply_transformations`](crate::analysis::mapmatching::kd_tree_feature_maps::KDTreeFeatureMaps::apply_transformations) searches a tree built on the
//! old retention times with the new ones, exactly as the source does.
//!
//! The source's class is a `DefaultParamHandler` without defaults whose
//! parameter check is disabled and whose `updateMembers_` does nothing: the
//! parameters its constructor accepts have no effect, so this port takes
//! none.

use crate::concept::math_functions::tolerance_window;
use crate::error::{Error, Result};
use crate::kernel::features::BaseFeature;

use super::qt_cluster_finder::QtInputMap;

/// Largest number of features one container holds. Checked on every
/// insertion; the source has no ceiling.
pub const MAX_FEATURES: usize = 1 << 26;

/// Largest number of tree levels one insertion may descend. A tree built by
/// inserting points sorted on both axes degenerates into a chain, which the
/// source walks recursively without bound; this port refuses instead.
pub const MAX_DEPTH: usize = 1 << 20;

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// One node of the libkdtree++ tree: the feature index and the children.
#[derive(Clone, Debug, Default)]
struct Node {
    index: usize,
    left: Option<usize>,
    right: Option<usize>,
}

/// The subset of libkdtree++'s `KDTree<2, KDTreeFeatureNode>` the source
/// uses: insertion, `optimise`, `find_within_range`, `size` and `clear`.
#[derive(Clone, Debug, Default)]
struct Tree {
    nodes: Vec<Node>,
    root: Option<usize>,
}

/// A query rectangle, libkdtree++'s `_Region` for two dimensions.
#[derive(Clone, Copy, Debug)]
struct Region {
    low: [f64; 2],
    high: [f64; 2],
}

impl Region {
    /// `_Region::encloses`: neither `v < low` nor `high < v` on any axis.
    fn encloses(&self, value: [f64; 2]) -> bool {
        !(0..2).any(|i| value[i] < self.low[i] || self.high[i] < value[i])
    }

    /// `_Region::intersects_with(const _Region&)`.
    fn intersects_with(&self, that: &Region) -> bool {
        !(0..2).any(|i| that.high[i] < self.low[i] || self.high[i] < that.low[i])
    }
}

impl Tree {
    fn len(&self) -> usize {
        self.nodes.len()
    }

    fn clear(&mut self) {
        self.nodes.clear();
        self.root = None;
    }

    /// libkdtree++ `insert` / `_M_insert`: descend left while the new value
    /// is less on the level's axis, right otherwise, and attach a leaf.
    fn insert(&mut self, index: usize, coord: &impl Fn(usize, usize) -> f64) -> Result<()> {
        let new = self.nodes.len();
        let Some(mut node) = self.root else {
            self.nodes.push(Node {
                index,
                ..Node::default()
            });
            self.root = Some(new);
            return Ok(());
        };
        let mut level = 0usize;
        loop {
            if level == MAX_DEPTH {
                return Err(bad(format!(
                    "kd-tree insertion descends more than {MAX_DEPTH} levels"
                )));
            }
            let axis = level % 2;
            let here = self.nodes[node].index;
            let go_left = coord(index, axis) < coord(here, axis);
            let child = if go_left {
                self.nodes[node].left
            } else {
                self.nodes[node].right
            };
            match child {
                Some(next) => {
                    node = next;
                    level += 1;
                }
                None => {
                    self.nodes.push(Node {
                        index,
                        ..Node::default()
                    });
                    if go_left {
                        self.nodes[node].left = Some(new);
                    } else {
                        self.nodes[node].right = Some(new);
                    }
                    return Ok(());
                }
            }
        }
    }

    /// The feature indices in the tree's in-order (libkdtree++'s iterator
    /// order from `begin()` to `end()`).
    fn in_order(&self) -> Vec<usize> {
        let mut result = Vec::with_capacity(self.nodes.len());
        let mut stack = Vec::new();
        let mut current = self.root;
        loop {
            while let Some(node) = current {
                stack.push(node);
                current = self.nodes[node].left;
            }
            let Some(node) = stack.pop() else {
                break;
            };
            result.push(self.nodes[node].index);
            current = self.nodes[node].right;
        }
        result
    }

    /// libkdtree++ `optimise()`: collect the values in order, clear, and
    /// insert the median of every range (by `nth_element` on the level's
    /// axis) before recursing into the left and then the right part.
    fn optimise(&mut self, coord: &impl Fn(usize, usize) -> f64) -> Result<()> {
        let mut values = self.in_order();
        self.clear();
        // Explicit stack of (begin, end, level); the right part is pushed
        // first so that the left part is built first, as in the recursion.
        let mut stack = vec![(0usize, values.len(), 0usize)];
        while let Some((begin, end, level)) = stack.pop() {
            if begin == end {
                continue;
            }
            let axis = level % 2;
            let middle = begin + (end - begin) / 2;
            nth_element(
                &mut values[begin..end],
                middle - begin,
                &|a: &usize, b: &usize| coord(*a, axis) < coord(*b, axis),
            );
            self.insert(values[middle], coord)?;
            if middle + 1 != end {
                stack.push((middle + 1, end, level + 1));
            }
            if middle != begin {
                stack.push((begin, middle, level + 1));
            }
        }
        Ok(())
    }

    /// libkdtree++ `find_within_range(region)`: the enclosed points in the
    /// tree's pre-order, pruning subtrees whose bounds miss the region.
    fn find_within_range(
        &self,
        region: &Region,
        coord: &impl Fn(usize, usize) -> f64,
        out: &mut Vec<usize>,
    ) {
        let Some(root) = self.root else {
            return;
        };
        let mut stack = vec![(root, *region, 0usize)];
        while let Some((node, bounds, level)) = stack.pop() {
            let index = self.nodes[node].index;
            let value = [coord(index, 0), coord(index, 1)];
            if region.encloses(value) {
                out.push(index);
            }
            let axis = level % 2;
            // Decide both children from this node's bounds, then push the
            // right one first so that the left subtree is reported first.
            let left = self.nodes[node].left.and_then(|child| {
                let mut b = bounds;
                b.high[axis] = value[axis];
                region.intersects_with(&b).then_some((child, b, level + 1))
            });
            let right = self.nodes[node].right.and_then(|child| {
                let mut b = bounds;
                b.low[axis] = value[axis];
                region.intersects_with(&b).then_some((child, b, level + 1))
            });
            if let Some(entry) = right {
                stack.push(entry);
            }
            if let Some(entry) = left {
                stack.push(entry);
            }
        }
    }
}

/// libstdc++'s `std::nth_element` (`__introselect` with a depth limit of
/// `2 * floor(log2(n))`, median-of-three pivoting, `__heap_select` fallback
/// and a final insertion sort), transcribed so that the element order it
/// leaves is the one the source's `optimise()` sees.
fn nth_element<T: Copy>(v: &mut [T], nth: usize, less: &impl Fn(&T, &T) -> bool) {
    let len = v.len();
    if len == 0 || nth >= len {
        return;
    }
    let mut depth_limit = 2 * (usize::BITS - 1 - len.leading_zeros()) as usize;
    let mut first = 0usize;
    let mut last = len;
    while last - first > 3 {
        if depth_limit == 0 {
            heap_select(&mut v[first..last], nth + 1 - first, less);
            v.swap(first, nth);
            return;
        }
        depth_limit -= 1;
        let cut = unguarded_partition_pivot(&mut v[first..last], less) + first;
        if cut <= nth {
            first = cut;
        } else {
            last = cut;
        }
    }
    insertion_sort(&mut v[first..last], less);
}

/// libstdc++ `__unguarded_partition_pivot`: returns the cut, relative to `v`.
fn unguarded_partition_pivot<T: Copy>(v: &mut [T], less: &impl Fn(&T, &T) -> bool) -> usize {
    let last = v.len();
    let mid = last / 2;
    move_median_to_first(v, 0, 1, mid, last - 1, less);
    // `__unguarded_partition(first + 1, last, first)`.
    let mut lo = 1usize;
    let mut hi = last;
    loop {
        // The source's scan is unguarded: a strict weak order stops it at the
        // pivot's sentinel. The bound only matters for coordinates that are
        // not numbers, where the source reads past the range.
        while lo < last && less(&v[lo], &v[0]) {
            lo += 1;
        }
        hi -= 1;
        while less(&v[0], &v[hi]) {
            hi -= 1;
        }
        if lo >= hi {
            return lo;
        }
        v.swap(lo, hi);
        lo += 1;
    }
}

/// libstdc++ `__move_median_to_first`.
fn move_median_to_first<T: Copy>(
    v: &mut [T],
    result: usize,
    a: usize,
    b: usize,
    c: usize,
    less: &impl Fn(&T, &T) -> bool,
) {
    let pick = if less(&v[a], &v[b]) {
        if less(&v[b], &v[c]) {
            b
        } else if less(&v[a], &v[c]) {
            c
        } else {
            a
        }
    } else if less(&v[a], &v[c]) {
        a
    } else if less(&v[b], &v[c]) {
        c
    } else {
        b
    };
    v.swap(result, pick);
}

/// libstdc++ `__insertion_sort` with `__unguarded_linear_insert`.
fn insertion_sort<T: Copy>(v: &mut [T], less: &impl Fn(&T, &T) -> bool) {
    for i in 1..v.len() {
        let value = v[i];
        if less(&value, &v[0]) {
            v.copy_within(0..i, 1);
            v[0] = value;
        } else {
            let mut hole = i;
            while less(&value, &v[hole - 1]) {
                v[hole] = v[hole - 1];
                hole -= 1;
            }
            v[hole] = value;
        }
    }
}

/// libstdc++ `__heap_select(first, middle, last)`.
fn heap_select<T: Copy>(v: &mut [T], middle: usize, less: &impl Fn(&T, &T) -> bool) {
    make_heap(&mut v[..middle], less);
    for i in middle..v.len() {
        if less(&v[i], &v[0]) {
            // `__pop_heap(first, middle, i)`.
            let value = v[i];
            v[i] = v[0];
            adjust_heap(&mut v[..middle], 0, value, less);
        }
    }
}

/// libstdc++ `__make_heap`.
fn make_heap<T: Copy>(v: &mut [T], less: &impl Fn(&T, &T) -> bool) {
    let len = v.len();
    if len < 2 {
        return;
    }
    let mut parent = (len - 2) / 2;
    loop {
        let value = v[parent];
        adjust_heap(v, parent, value, less);
        if parent == 0 {
            return;
        }
        parent -= 1;
    }
}

/// libstdc++ `__adjust_heap` followed by `__push_heap`.
fn adjust_heap<T: Copy>(v: &mut [T], hole: usize, value: T, less: &impl Fn(&T, &T) -> bool) {
    let len = v.len();
    let top = hole;
    let mut hole = hole;
    let mut second = hole;
    while second < (len - 1) / 2 {
        second = 2 * (second + 1);
        if less(&v[second], &v[second - 1]) {
            second -= 1;
        }
        v[hole] = v[second];
        hole = second;
    }
    if len % 2 == 0 && second == (len - 2) / 2 {
        second = 2 * (second + 1);
        v[hole] = v[second - 1];
        hole = second - 1;
    }
    // `__push_heap(first, hole, top, value)`.
    while hole > top {
        let parent = (hole - 1) / 2;
        if !less(&v[parent], &value) {
            break;
        }
        v[hole] = v[parent];
        hole = parent;
    }
    v[hole] = value;
}

/// A kd-tree entry: a feature of a [`KDTreeFeatureMaps`] by its global
/// index (source `KDTreeFeatureNode`). It carries no feature data; its
/// coordinates are read from the container.
#[derive(Clone, Copy, Debug)]
pub struct KDTreeFeatureNode<'m, 'a> {
    data: &'m KDTreeFeatureMaps<'a>,
    index: usize,
}

impl<'m, 'a> KDTreeFeatureNode<'m, 'a> {
    /// A node referencing feature `index` of `data` (source constructor).
    pub fn new(data: &'m KDTreeFeatureMaps<'a>, index: usize) -> Self {
        Self { data, index }
    }

    /// The global index of the referenced feature (source `getIndex`).
    pub fn index(&self) -> usize {
        self.index
    }

    /// Coordinate `i`: the (possibly transformed) RT for 0, the m/z for 1
    /// (source `operator[]`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for any other `i` (source
    /// `Exception::ElementNotFound`, "Indices other than 0 (RT) and 1 (m/z)
    /// are not allowed!") and for an index outside the container.
    pub fn get(&self, i: usize) -> Result<f64> {
        match i {
            0 => self.data.rt(self.index),
            1 => self.data.mz(self.index),
            _ => Err(bad(
                "Indices other than 0 (RT) and 1 (m/z) are not allowed!",
            )),
        }
    }
}

/// Features of several maps in one (RT, m/z) kd-tree (source
/// `KDTreeFeatureMaps`).
#[derive(Clone, Debug, Default)]
pub struct KDTreeFeatureMaps<'a> {
    features: Vec<&'a BaseFeature>,
    map_index: Vec<usize>,
    rt: Vec<f64>,
    num_maps: usize,
    tree: Tree,
}

impl<'a> KDTreeFeatureMaps<'a> {
    /// An empty container with no features and zero maps (source default
    /// constructor).
    pub fn new() -> Self {
        Self::default()
    }

    /// A container over all features of `maps` (source
    /// `KDTreeFeatureMaps(maps, param)`, whose parameters have no effect).
    ///
    /// # Errors
    ///
    /// As [`Self::add_maps`].
    pub fn from_maps<M: QtInputMap>(maps: &'a [M]) -> Result<Self> {
        let mut result = Self::new();
        result.add_maps(maps)?;
        Ok(result)
    }

    /// Insert every feature of `maps` under its map's index, record the
    /// number of maps, and balance the tree (source `addMaps`).
    ///
    /// # Errors
    ///
    /// As [`Self::add_feature`].
    pub fn add_maps<M: QtInputMap>(&mut self, maps: &'a [M]) -> Result<()> {
        self.num_maps = maps.len();
        for (i, map) in maps.iter().enumerate() {
            for feature in map.base_features() {
                self.add_feature(i, feature)?;
            }
        }
        self.optimize_tree()
    }

    /// Insert one feature of map `map_index` (source `addFeature`). Its RT
    /// is cached; call [`Self::optimize_tree`] after many insertions.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] beyond [`MAX_FEATURES`] features or
    /// [`MAX_DEPTH`] tree levels.
    pub fn add_feature(&mut self, map_index: usize, feature: &'a BaseFeature) -> Result<()> {
        if self.features.len() == MAX_FEATURES {
            return Err(bad(format!(
                "a kd-tree of feature maps holds at most {MAX_FEATURES} features"
            )));
        }
        self.map_index.push(map_index);
        self.features.push(feature);
        self.rt.push(feature.rt);
        let index = self.features.len() - 1;
        let (rt, features) = (&self.rt, &self.features);
        let coord = |i: usize, axis: usize| if axis == 0 { rt[i] } else { features[i].mz };
        let inserted = self.tree.insert(index, &coord);
        if inserted.is_err() {
            self.map_index.pop();
            self.features.pop();
            self.rt.pop();
        }
        inserted
    }

    fn out_of_range(i: usize) -> Error {
        Error::InvalidRange(format!("kd-tree feature index {i} out of range"))
    }

    /// The feature with global index `i` (source `feature`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] for an index at or beyond
    /// [`Self::size`]; the source does not check.
    pub fn feature(&self, i: usize) -> Result<&'a BaseFeature> {
        self.features
            .get(i)
            .copied()
            .ok_or_else(|| Self::out_of_range(i))
    }

    /// The (possibly transformed) retention time of feature `i` (source
    /// `rt`).
    ///
    /// # Errors
    ///
    /// As [`Self::feature`].
    pub fn rt(&self, i: usize) -> Result<f64> {
        self.rt.get(i).copied().ok_or_else(|| Self::out_of_range(i))
    }

    /// The m/z of feature `i`, read through the feature (source `mz`).
    ///
    /// # Errors
    ///
    /// As [`Self::feature`].
    pub fn mz(&self, i: usize) -> Result<f64> {
        Ok(self.feature(i)?.mz)
    }

    /// The intensity of feature `i` (source `intensity`).
    ///
    /// # Errors
    ///
    /// As [`Self::feature`].
    pub fn intensity(&self, i: usize) -> Result<f32> {
        Ok(self.feature(i)?.intensity)
    }

    /// The charge of feature `i` (source `charge`).
    ///
    /// # Errors
    ///
    /// As [`Self::feature`].
    pub fn charge(&self, i: usize) -> Result<i32> {
        Ok(self.feature(i)?.charge)
    }

    /// The map index feature `i` was added with (source `mapIndex`).
    ///
    /// # Errors
    ///
    /// As [`Self::feature`].
    pub fn map_index(&self, i: usize) -> Result<usize> {
        self.map_index
            .get(i)
            .copied()
            .ok_or_else(|| Self::out_of_range(i))
    }

    /// The number of features stored (source `size`).
    pub fn size(&self) -> usize {
        self.features.len()
    }

    /// The number of points in the tree (source `treeSize`).
    pub fn tree_size(&self) -> usize {
        self.tree.len()
    }

    /// The number of maps of the last [`Self::add_maps`] (source `numMaps`).
    pub fn num_maps(&self) -> usize {
        self.num_maps
    }

    /// Discard all features and clear the tree (source `clear`). As in the
    /// source, the number of maps is kept.
    pub fn clear(&mut self) {
        self.features.clear();
        self.map_index.clear();
        self.rt.clear();
        self.tree.clear();
    }

    /// Rebalance the tree (source `optimizeTree`).
    ///
    /// # Errors
    ///
    /// Only the [`MAX_DEPTH`] ceiling, which a rebalanced tree cannot reach.
    pub fn optimize_tree(&mut self) -> Result<()> {
        let (rt, features) = (&self.rt, &self.features);
        let coord = |i: usize, axis: usize| if axis == 0 { rt[i] } else { features[i].mz };
        self.tree.optimise(&coord)
    }

    /// The features compatible with feature `index`, appended to
    /// `result_indices` (source `getNeighborhood`): a [`Self::query_region`]
    /// with the RT window `rt ± rt_tol` and the m/z window of
    /// [`tolerance_window`] (ppm when `mz_ppm`), excluding the features of
    /// `index`'s own map unless `include_features_from_same_map`. With a
    /// non-negative `max_pairwise_log_fc`, a candidate is kept only when
    /// `|log10(I_candidate / I_index)| <= max_pairwise_log_fc`, which drops
    /// candidates whose fold change is not a number.
    ///
    /// # Errors
    ///
    /// As [`Self::feature`] for `index`, and the [`tolerance_window`] errors
    /// for a non-finite window (the source computes it regardless).
    #[allow(clippy::too_many_arguments, reason = "the source signature")]
    pub fn get_neighborhood(
        &self,
        index: usize,
        result_indices: &mut Vec<usize>,
        rt_tol: f64,
        mz_tol: f64,
        mz_ppm: bool,
        include_features_from_same_map: bool,
        max_pairwise_log_fc: f64,
    ) -> Result<()> {
        let (rt_low, rt_high) = tolerance_window(self.rt(index)?, rt_tol, false)?;
        let (mz_low, mz_high) = tolerance_window(self.mz(index)?, mz_tol, mz_ppm)?;
        let ignored = if include_features_from_same_map {
            None
        } else {
            Some(self.map_index(index)?)
        };
        let mut found = Vec::new();
        self.query_region(rt_low, rt_high, mz_low, mz_high, &mut found, ignored);
        if max_pairwise_log_fc < 0.0 {
            result_indices.extend(found);
        } else {
            let int_1 = f64::from(self.features[index].intensity);
            for i in found {
                let int_2 = f64::from(self.features[i].intensity);
                let abs_log_fc = (int_2 / int_1).log10().abs();
                if abs_log_fc <= max_pairwise_log_fc {
                    result_indices.push(i);
                }
            }
        }
        Ok(())
    }

    /// The features inside the rectangle `[rt_low, rt_high] x [mz_low,
    /// mz_high]`, bounds included, in the tree's pre-order, replacing the
    /// contents of `result_indices` (source `queryRegion`). Features of map
    /// `ignored_map_index` are left out; `None` is the source's default
    /// `numeric_limits<Size>::max()`, which keeps all.
    pub fn query_region(
        &self,
        rt_low: f64,
        rt_high: f64,
        mz_low: f64,
        mz_high: f64,
        result_indices: &mut Vec<usize>,
        ignored_map_index: Option<usize>,
    ) {
        let region = Region {
            low: [rt_low, mz_low],
            high: [rt_high, mz_high],
        };
        let (rt, features) = (&self.rt, &self.features);
        let coord = |i: usize, axis: usize| if axis == 0 { rt[i] } else { features[i].mz };
        let mut found = Vec::new();
        self.tree.find_within_range(&region, &coord, &mut found);
        result_indices.clear();
        result_indices.extend(
            found
                .into_iter()
                .filter(|&i| ignored_map_index.is_none_or(|m| self.map_index[i] != m)),
        );
    }

    /// Replace every cached retention time by `transform(map index,
    /// feature RT)` (source `applyTransformations`, which evaluates the
    /// feature map's `TransformationModelLowess` at the feature's original
    /// RT). The features and the tree are not changed.
    ///
    /// # Errors
    ///
    /// The first error `transform` returns; the cached times are then
    /// unchanged.
    pub fn apply_transformations(
        &mut self,
        transform: impl Fn(usize, f64) -> Result<f64>,
    ) -> Result<()> {
        let mut rt = Vec::with_capacity(self.rt.len());
        for (i, feature) in self.features.iter().enumerate() {
            rt.push(transform(self.map_index[i], feature.rt)?);
        }
        self.rt = rt;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tier 4: `nth_element` leaves the `nth` element in its sorted
    /// position with no greater element before and no smaller one after, for
    /// every position of several inputs, including the heap-select path
    /// (forced by a depth limit of zero through a crafted length is not
    /// reachable, so the invariant is checked on adversarial orders).
    #[test]
    fn nth_element_partitions() {
        let inputs: Vec<Vec<i32>> = vec![
            vec![5, 1, 4, 2, 3],
            vec![1, 1, 1, 1, 1, 1, 1],
            (0..50).rev().collect(),
            (0..50).map(|i| (i * 37) % 23).collect(),
            vec![3, 1, 2],
        ];
        for input in inputs {
            let mut sorted = input.clone();
            sorted.sort();
            for nth in 0..input.len() {
                let mut v = input.clone();
                nth_element(&mut v, nth, &|a: &i32, b: &i32| a < b);
                assert_eq!(v[nth], sorted[nth]);
                assert!(v[..nth].iter().all(|x| *x <= v[nth]));
                assert!(v[nth + 1..].iter().all(|x| *x >= v[nth]));
            }
        }
    }

    #[test]
    fn heap_select_keeps_smallest() {
        let mut v: Vec<i32> = vec![9, 3, 7, 1, 8, 2, 6, 4, 5, 0];
        heap_select(&mut v, 4, &|a: &i32, b: &i32| a < b);
        let mut head = v[..4].to_vec();
        head.sort();
        assert_eq!(head, vec![0, 1, 2, 3]);
        // A max-heap on the first four: the root is the largest of them.
        assert_eq!(v[0], 3);
    }
}
