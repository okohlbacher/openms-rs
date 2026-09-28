// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Quality measures and partitions of a hierarchical clustering.
//!
//! Native equivalent of `ML/CLUSTERING/ClusterAnalyzer.h`. The source keeps
//! its working state in `std::map`s indexed with `operator[]`, which inserts a
//! default entry for every key it reads; several loops run to the map's
//! *current* size and therefore visit the keys those reads created. The port
//! uses ordered maps with the same insert-on-read rule, so every loop visits
//! the same keys in the same order. See `docs/CLUSTERING_SUPPORT.md`.

use std::collections::{BTreeMap, BTreeSet};

use super::cluster_functor::{at, at_mut};
use crate::data_structures::binary_tree_node::BinaryTreeNode;
use crate::data_structures::distance_matrix::DistanceMatrix;
use crate::{Error, Result};

/// Bundles analysis tools for a clustering given as a sequence of
/// [`BinaryTreeNode`]s (source `ClusterAnalyzer`).
///
/// The source class is stateless; its constructors and destructor map to
/// `Default`, `Clone` and `Copy`. Its private, empty assignment operator is
/// not ported.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ClusterAnalyzer;

fn invalid_parameter(message: &str) -> Error {
    Error::InvalidValue(message.to_owned())
}

/// `std::max(a, b)`: `b` only when `a < b`.
fn cpp_max(a: f32, b: f32) -> f32 {
    if a < b { b } else { a }
}

/// `map[key]` on a `std::map<Size, T>`: inserts `T()` when absent.
fn entry<T: Default>(map: &mut BTreeMap<usize, T>, key: usize) -> &mut T {
    map.entry(key).or_default()
}

/// Sum of `original(element, x)` over `members`, accumulated in `f32` in order.
fn distance_sum(original: &DistanceMatrix<f32>, element: usize, members: &[usize]) -> Result<f32> {
    let mut sum = 0.0_f32;
    for &member in members {
        sum += original.get(element, member)?;
    }
    Ok(sum)
}

impl ClusterAnalyzer {
    /// A cluster analyzer (source default constructor).
    pub const fn new() -> Self {
        Self
    }

    /// The average silhouette width of every clustering step (source
    /// `averageSilhouetteWidth`).
    ///
    /// The width is computed for each step from the first (`n - 1` clusters)
    /// to the last (one cluster, whose width is 0 by definition), so the
    /// result has one entry per node of `tree`. `original` is the distance
    /// matrix the clustering started from. Leaves are collected from the nodes
    /// up to and including the first placeholder node (distance `-1`); the
    /// widths are then accumulated over all nodes but the last, as in the
    /// source, and each step's sum is divided by `tree.len() + 1`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an empty tree (source
    /// `Exception::InvalidParameter`: "tree is empty but minimal clustering
    /// hirachy has at least one level"), and [`Error::InvalidRange`] when a
    /// node names an element outside `original`, which the source reports as
    /// `Exception::OutOfRange` from `DistanceMatrix::getValue`.
    pub fn average_silhouette_width(
        &self,
        tree: &[BinaryTreeNode],
        original: &DistanceMatrix<f32>,
    ) -> Result<Vec<f32>> {
        if tree.is_empty() {
            return Err(invalid_parameter(
                "tree is empty but minimal clustering hirachy has at least one level",
            ));
        }

        let mut average_silhouette_widths = Vec::with_capacity(tree.len());
        let mut interdist_i: BTreeMap<usize, f32> = BTreeMap::new();
        let mut cluster_with_interdist: BTreeMap<usize, usize> = BTreeMap::new();
        let mut intradist_i: BTreeMap<usize, f32> = BTreeMap::new();

        let mut leafs = BTreeSet::new();
        for node in tree {
            for child in [node.left_child, node.right_child] {
                leafs.insert(child);
                interdist_i.entry(child).or_insert(f32::MAX);
                cluster_with_interdist.entry(child).or_insert(0);
                intradist_i.entry(child).or_insert(0.0);
            }
            if node.distance == -1.0 {
                break;
            }
        }
        let leafs: Vec<usize> = leafs.into_iter().collect();

        // Initial values for interdist_i and cluster_with_interdist.
        for (position, &leaf) in leafs.iter().enumerate().skip(1) {
            for &other in &leafs[..position] {
                let distance = original.get(leaf, other)?;
                if distance < *entry(&mut interdist_i, leaf) {
                    interdist_i.insert(leaf, distance);
                    cluster_with_interdist.insert(leaf, other);
                }
                if distance < *entry(&mut interdist_i, other) {
                    interdist_i.insert(other, distance);
                    cluster_with_interdist.insert(other, leaf);
                }
            }
        }

        let mut clusters: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for &leaf in &leafs {
            entry(&mut clusters, leaf).push(leaf);
        }

        for step in &tree[..tree.len() - 1] {
            let (left, right) = (step.left_child, step.right_child);
            for &element in &leafs {
                let in_left = entry(&mut clusters, left).contains(&element);
                let in_right = entry(&mut clusters, right).contains(&element);

                if !in_left && !in_right {
                    let nearest = *entry(&mut cluster_with_interdist, element);
                    if left != nearest && right != nearest {
                        let left_members = clusters[&left].clone();
                        let right_members = clusters[&right].clone();
                        let mut merged = distance_sum(original, element, &left_members)?;
                        merged += distance_sum(original, element, &right_members)?;
                        merged /= (left_members.len() + right_members.len()) as f32;
                        if merged < *entry(&mut interdist_i, element) {
                            interdist_i.insert(element, merged);
                            cluster_with_interdist.insert(element, left);
                        }
                    } else {
                        // The merged cluster that does NOT contain the nearest cluster.
                        let k = if right != nearest { right } else { left };
                        let k_members = entry(&mut clusters, k).clone();
                        let nearest_size = entry(&mut clusters, nearest).len();
                        let mut merged = distance_sum(original, element, &k_members)?;
                        merged += nearest_size as f32 * *entry(&mut interdist_i, element);
                        merged /= (k_members.len() + nearest_size) as f32;
                        if merged <= *entry(&mut interdist_i, element) {
                            interdist_i.insert(element, merged);
                            cluster_with_interdist.insert(element, left);
                        } else {
                            interdist_i.insert(element, merged);
                            cluster_with_interdist.insert(element, left);
                            let mut u = 0;
                            while u < clusters.len() {
                                if u != left && u != right {
                                    let members = entry(&mut clusters, u);
                                    if !members.is_empty() && !members.contains(&element) {
                                        let members = members.clone();
                                        let mut average = 0.0_f32;
                                        for &member in &members {
                                            average += original.get(member, element)?;
                                        }
                                        average /= members.len() as f32;
                                        if average < *entry(&mut interdist_i, element) {
                                            interdist_i.insert(element, average);
                                            cluster_with_interdist.insert(element, u);
                                        }
                                    }
                                }
                                u += 1;
                            }
                        }
                    }
                } else {
                    // l contains the element, k is the other merged cluster.
                    let (l, k) = if in_left {
                        (left, right)
                    } else {
                        (right, left)
                    };
                    let l_less_one = entry(&mut clusters, l).len().wrapping_sub(1);
                    let k_members = entry(&mut clusters, k).clone();
                    if k != *entry(&mut cluster_with_interdist, element) {
                        let mut intra = *entry(&mut intradist_i, element) * l_less_one as f32;
                        intra += distance_sum(original, element, &k_members)?;
                        intra /= k_members.len().wrapping_add(l_less_one) as f32;
                        intradist_i.insert(element, intra);
                    } else {
                        let mut intra = *entry(&mut intradist_i, element) * l_less_one as f32;
                        intra += k_members.len() as f32 * *entry(&mut interdist_i, element);
                        intra /= k_members.len().wrapping_add(l_less_one) as f32;
                        intradist_i.insert(element, intra);
                        interdist_i.insert(element, f32::MAX);
                        let mut u = 0;
                        while u < clusters.len() {
                            if u != l && u != k {
                                let members = entry(&mut clusters, u);
                                if !members.is_empty() {
                                    let members = members.clone();
                                    let mut average = 0.0_f32;
                                    for &member in &members {
                                        average += original.get(member, element)?;
                                    }
                                    average /= members.len() as f32;
                                    if average < *entry(&mut interdist_i, element) {
                                        interdist_i.insert(element, average);
                                        cluster_with_interdist.insert(element, u);
                                    }
                                }
                            }
                            u += 1;
                        }
                    }
                }
            }

            // Redo the clustering following the tree.
            let moved = entry(&mut clusters, right).clone();
            entry(&mut clusters, left).extend(moved);
            entry(&mut clusters, right).clear();

            let mut average_overall_silhouette = 0.0_f32;
            let mut g = 0;
            while g < clusters.len() {
                let members = entry(&mut clusters, g).clone();
                if members.len() > 1 {
                    for member in members {
                        let inter = *entry(&mut interdist_i, member);
                        if inter != 0.0 {
                            let intra = *entry(&mut intradist_i, member);
                            average_overall_silhouette += (inter - intra) / cpp_max(inter, intra);
                        }
                    }
                }
                g += 1;
            }
            average_silhouette_widths.push(average_overall_silhouette / (tree.len() + 1) as f32);
        }
        average_silhouette_widths.push(0.0);
        Ok(average_silhouette_widths)
    }

    /// Dunn's index of every clustering step (source `dunnIndices`).
    ///
    /// Each step's index is the minimal inter-cluster distance divided by the
    /// maximal intra-cluster distance (0 while the latter is 0); the last entry
    /// is 0. With `tree_from_single_linkage` (source default `false`) the
    /// tree's own distances are used as the minimal inter-cluster distances,
    /// which is valid only for a [`SingleLinkage`](super::single_linkage::SingleLinkage)
    /// tree, whose distances are those minima in increasing order.
    ///
    /// The incremental bookkeeping is the source's, including its quirk of
    /// resetting the merged-away cluster's minimum to `-1` before it is
    /// compared with the surviving cluster's.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an empty tree (source
    /// `Exception::InvalidParameter`) or a node naming an element outside
    /// `original`, which the source indexes out of bounds; and
    /// [`Error::InvalidRange`] from the distance lookups.
    pub fn dunn_indices(
        &self,
        tree: &[BinaryTreeNode],
        original: &DistanceMatrix<f32>,
        tree_from_single_linkage: bool,
    ) -> Result<Vec<f32>> {
        if tree.is_empty() {
            return Err(invalid_parameter(
                "tree is empty but minimal clustering hirachy has at least one level",
            ));
        }
        let mut all_dunn_indices = Vec::with_capacity(tree.len() + 1);

        let mut leafs = BTreeSet::new();
        for node in tree {
            leafs.insert(node.left_child);
            leafs.insert(node.right_child);
        }
        let leafs: Vec<usize> = leafs.into_iter().collect();

        let n = original.dimension_size();
        let mut clusters: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut min_inter: Vec<(f32, usize)> = vec![(-1.0, 0); n];
        for (position, &leaf) in leafs.iter().enumerate() {
            at_mut(&mut clusters, leaf)?.push(leaf);
            for &other in leafs[..position].iter().chain(&leafs[position + 1..]) {
                let d = original.get(leaf, other)?;
                let current = at_mut(&mut min_inter, leaf)?;
                if d < current.0 || current.0 == -1.0 {
                    *current = (d, other);
                }
            }
        }
        let mut min_index = min_inter_index(&min_inter);

        let mut max_intracluster_distance = 0.0_f32;
        for step in 0..tree.len() - 1 {
            let (left, right) = (tree[step].left_child, tree[step].right_child);
            let left_members = at(&clusters, left)?.clone();
            let right_members = at(&clusters, right)?.clone();
            for &x in &left_members {
                for &y in &right_members {
                    let value = original.get(x, y)?;
                    if value > max_intracluster_distance {
                        max_intracluster_distance = value;
                    }
                }
            }

            let moved = at(&clusters, right)?.clone();
            at_mut(&mut clusters, left)?.extend(moved);
            at_mut(&mut clusters, right)?.clear();

            if !tree_from_single_linkage {
                *at_mut(&mut min_inter, right)? = (-1.0, 0);
                let current = *at(&min_inter, min_index)?;
                if (min_index == right && current.1 == left)
                    || (min_index == left && current.1 == right)
                {
                    at_mut(&mut min_inter, left)?.0 = f32::MAX;
                    let merged = clusters[left].clone();
                    for &element in &merged {
                        for (k, members) in clusters.iter().enumerate() {
                            if k == left {
                                continue;
                            }
                            for &member in members {
                                let value = original.get(element, member)?;
                                if value < min_inter[left].0 {
                                    min_inter[left] = (value, k);
                                }
                            }
                        }
                    }
                    min_index = min_inter_index(&min_inter);
                } else {
                    let right_min = min_inter[right];
                    if right_min.0 < min_inter[left].0 {
                        min_inter[left] = right_min;
                    }
                }
                for entry in &mut min_inter {
                    if entry.1 == right {
                        entry.1 = left;
                    }
                }
            }

            let mut dunn_index = 0.0_f32;
            if max_intracluster_distance > 0.0 {
                let numerator = if tree_from_single_linkage {
                    tree[step + 1].distance
                } else {
                    at(&min_inter, min_index)?.0
                };
                dunn_index = numerator / max_intracluster_distance;
            }
            all_dunn_indices.push(dunn_index);
        }
        all_dunn_indices.push(0.0);
        Ok(all_dunn_indices)
    }

    /// The cohesion of every cluster of a partition (source `cohesion`).
    ///
    /// A cluster's cohesion is the average of its pairwise element distances;
    /// a singleton cluster's cohesion is the average of all pairwise distances
    /// in `original`. The result is in the order of `clusters`. Averages are
    /// accumulated in `f32` as in the source; an empty cluster yields NaN.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `clusters` is empty or has more clusters
    /// than `original` has elements (source `Exception::InvalidParameter`,
    /// "invalid clustering"), and [`Error::InvalidRange`] for an element index
    /// outside `original` (source `Exception::OutOfRange`).
    pub fn cohesion(
        &self,
        clusters: &[Vec<usize>],
        original: &DistanceMatrix<f32>,
    ) -> Result<Vec<f32>> {
        let n = original.dimension_size();
        if clusters.is_empty() || clusters.len() > n {
            return Err(invalid_parameter("invalid clustering"));
        }
        let mut av_dist = 0.0_f32;
        for i in 0..n {
            for j in i + 1..n {
                av_dist += original.get(i, j)?;
            }
        }
        av_dist /= (n as f32 * (n as f64 - 1.0) as f32) / 2.0_f32;

        let mut cohesions = Vec::with_capacity(clusters.len());
        for cluster in clusters {
            let mut av_c_dist = 0.0_f32;
            for j in 0..cluster.len() {
                for k in 0..j {
                    av_c_dist += original.get(cluster[j], cluster[k])?;
                }
            }
            av_c_dist /= (cluster.len() as f32 * (cluster.len() as f64 - 1.0) as f32) / 2.0_f32;
            if cluster.len() == 1 {
                av_c_dist = av_dist;
            }
            cohesions.push(av_c_dist);
        }
        Ok(cohesions)
    }

    /// The average aberration from the average cluster population of the
    /// partition with `cluster_quantity` clusters (source
    /// `averagePopulationAberration`).
    ///
    /// The partition is the one [`cut`](Self::cut) gives; the average
    /// population is `(tree.len() + 1) / cluster_quantity`, and the result the
    /// mean absolute deviation of the non-empty clusters' sizes from it.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `cluster_quantity` is 0 ("minimal
    /// partition contains one cluster, not zero") or not below `tree.len() +
    /// 1` ("maximal partition contains singleton clusters, further separation
    /// is not possible"), both source `Exception::InvalidParameter`, or when a
    /// node names an element outside `0..=tree.len()`, which the source indexes
    /// out of bounds.
    pub fn average_population_aberration(
        &self,
        cluster_quantity: usize,
        tree: &[BinaryTreeNode],
    ) -> Result<f32> {
        if cluster_quantity == 0 {
            return Err(invalid_parameter(
                "minimal partition contains one cluster, not zero",
            ));
        }
        if cluster_quantity > tree.len() {
            return Err(invalid_parameter(
                "maximal partition contains singleton clusters, further separation is not possible",
            ));
        }
        let mut clusters: Vec<Vec<usize>> = (0..=tree.len()).map(|i| vec![i]).collect();
        for step in &tree[..tree.len() + 1 - cluster_quantity] {
            let moved = at(&clusters, step.right_child)?.clone();
            at_mut(&mut clusters, step.left_child)?.extend(moved);
            at_mut(&mut clusters, step.right_child)?.clear();
        }
        let average = (tree.len() + 1) as f32 / cluster_quantity as f32;
        let mut aberration = 0.0_f32;
        let mut cluster_number = 0.0_f32;
        for cluster in &clusters {
            if !cluster.is_empty() {
                aberration += (cluster.len() as f32 - average).abs();
                cluster_number += 1.0;
            }
        }
        Ok(aberration / cluster_number)
    }

    /// The partition into `cluster_quantity` clusters (source `cut`, first
    /// overload).
    ///
    /// The clustering is replayed for `tree.len() + 1 - cluster_quantity`
    /// steps, stopping early at the first placeholder node (distance `-1`).
    /// Each cluster lists its element indices in ascending order and the
    /// clusters are sorted lexicographically. To fetch all clusters formed
    /// under a threshold, count the nodes whose distance is not `-1` and pass
    /// `data.len()` minus that count.
    ///
    /// The source appends to a caller-supplied vector and then sorts, reverses,
    /// truncates or pads with empty clusters to `cluster_quantity` and sorts
    /// again, all over that vector's previous content as well; this returns a
    /// fresh vector, as the source does for an empty argument. When fewer
    /// clusters exist than requested, empty clusters pad the result and sort
    /// first.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `cluster_quantity` is 0 ("You requested 0
    /// clusters. Minimal partitioning contains one cluster, not zero.") or
    /// exceeds `tree.len() + 1` ("Cluster count exceeds node count. No
    /// partitioning exists."), both source `Exception::InvalidParameter`.
    pub fn cut(&self, cluster_quantity: usize, tree: &[BinaryTreeNode]) -> Result<Vec<Vec<usize>>> {
        if cluster_quantity == 0 {
            return Err(invalid_parameter(
                "You requested 0 clusters. Minimal partitioning contains one cluster, not zero.",
            ));
        }
        if cluster_quantity > tree.len() + 1 {
            return Err(invalid_parameter(
                "Cluster count exceeds node count. No partitioning exists.",
            ));
        }
        let mut cluster_map: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for node in tree {
            cluster_map.insert(node.left_child, vec![node.left_child]);
            cluster_map.insert(node.right_child, vec![node.right_child]);
        }
        for step in &tree[..tree.len() + 1 - cluster_quantity] {
            if step.distance == -1.0 {
                break;
            }
            let moved = entry(&mut cluster_map, step.right_child).clone();
            entry(&mut cluster_map, step.left_child).extend(moved);
            entry(&mut cluster_map, step.right_child).clear();
        }
        let mut clusters: Vec<Vec<usize>> = cluster_map
            .into_values()
            .filter(|members| !members.is_empty())
            .collect();
        for cluster in &mut clusters {
            cluster.sort_unstable();
        }
        clusters.sort();
        clusters.reverse();
        clusters.resize(cluster_quantity, Vec::new());
        clusters.sort();
        Ok(clusters)
    }

    /// The subtrees of the partition into `cluster_quantity` clusters (source
    /// `cut`, second overload).
    ///
    /// Subtree `c` holds, in tree order, the nodes among the first
    /// `tree.len() + 1 - cluster_quantity` that touch an element of cluster `c`
    /// of [`cut`](Self::cut); each node goes to the first cluster it touches.
    /// The source clears the caller's vector first, so a fresh result is the
    /// same.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `cluster_quantity` is 0 ("minimal
    /// partition contains one cluster, not zero") or not below `tree.len() +
    /// 1` ("maximal partition contains singleton clusters, further separation
    /// is not possible"), both source `Exception::InvalidParameter`.
    pub fn cut_subtrees(
        &self,
        cluster_quantity: usize,
        tree: &[BinaryTreeNode],
    ) -> Result<Vec<Vec<BinaryTreeNode>>> {
        if cluster_quantity == 0 {
            return Err(invalid_parameter(
                "minimal partition contains one cluster, not zero",
            ));
        }
        if cluster_quantity > tree.len() {
            return Err(invalid_parameter(
                "maximal partition contains singleton clusters, further separation is not possible",
            ));
        }
        let mut subtrees = vec![Vec::new(); cluster_quantity];
        let clusters = self.cut(cluster_quantity, tree)?;
        let mut remaining: Vec<BinaryTreeNode> = tree[..tree.len() + 1 - cluster_quantity].to_vec();
        for (cluster, members) in clusters.iter().enumerate() {
            let mut kept = Vec::with_capacity(remaining.len());
            for node in remaining {
                if members.contains(&node.left_child) || members.contains(&node.right_child) {
                    subtrees[cluster].push(node);
                } else {
                    kept.push(node);
                }
            }
            remaining = kept;
        }
        Ok(subtrees)
    }

    /// The hierarchy of `tree` as a Newick string (source `newickTree`).
    ///
    /// Leaves are element indices; each step wraps the two merged subtrees as
    /// `( left , right )`. Clusters left unmerged at the end are joined the
    /// same way, left to right.
    ///
    /// `distance_text` is the source's `include_distance` (default `false`):
    /// with `Some(format)`, each merged subtree is followed by `:` and its
    /// step's distance as `format` writes it, and each finally joined root by
    /// `:1`. The source writes the distance with `StringUtils::toStr(float)`,
    /// which is `to_str_f32` in the format module `file_info::text_format`; `math`
    /// may not depend on `format`, so the caller passes it. `Some(&to_str_f32)`
    /// reproduces the source's text exactly.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an empty tree, where the source dereferences
    /// the end of an empty set.
    pub fn newick_tree(
        &self,
        tree: &[BinaryTreeNode],
        distance_text: Option<&dyn Fn(f32) -> String>,
    ) -> Result<String> {
        let mut leafs = BTreeSet::new();
        for node in tree {
            leafs.insert(node.left_child);
            leafs.insert(node.right_child);
        }
        let Some(&largest) = leafs.last() else {
            return Err(invalid_parameter("tree is empty"));
        };
        let mut clusters = vec![String::new(); largest + 1];
        for &leaf in &leafs {
            clusters[leaf] = leaf.to_string();
        }

        for node in tree {
            let right = std::mem::take(&mut clusters[node.right_child]);
            let mut text = String::from("( ");
            text.push_str(&clusters[node.left_child]);
            if let Some(format) = distance_text {
                text.push(':');
                text.push_str(&format(node.distance));
            }
            text.push_str(" , ");
            text.push_str(&right);
            if let Some(format) = distance_text {
                text.push(':');
                text.push_str(&format(node.distance));
            }
            text.push_str(" )");
            clusters[node.left_child] = text;
            // The source assigns "" to the right child after appending it, so a
            // node with equal children leaves an empty string behind.
            if node.right_child == node.left_child {
                clusters[node.right_child] = String::new();
            }
        }

        let first_filled = clusters
            .iter()
            .position(|text| !text.is_empty())
            .unwrap_or(0);
        for i in first_filled + 1..clusters.len() {
            if !clusters[i].is_empty() {
                let mut text = String::from("( ");
                text.push_str(&clusters[first_filled]);
                if distance_text.is_some() {
                    text.push_str(":1");
                }
                text.push_str(" , ");
                text.push_str(&clusters[i]);
                if distance_text.is_some() {
                    text.push_str(":1");
                }
                text.push_str(" )");
                clusters[first_filled] = text;
            }
        }
        Ok(std::mem::take(&mut clusters[first_filled]))
    }
}

/// The index of the smallest set minimum, skipping unset (`-1`) entries, as
/// the source's two identical scans in `dunnIndices` compute it.
fn min_inter_index(min_inter: &[(f32, usize)]) -> usize {
    let mut index = 0;
    for i in 1..min_inter.len() {
        if min_inter[index].0 == -1.0
            || (min_inter[i].0 != -1.0 && min_inter[i].0 < min_inter[index].0)
        {
            index = i;
        }
    }
    index
}

/// Whether `x` sorts before `y` by distance (source free function
/// `compareBinaryTreeNode`, `x.distance < y.distance`, for use with sorting).
pub fn compare_binary_tree_node(x: &BinaryTreeNode, y: &BinaryTreeNode) -> bool {
    x.distance < y.distance
}
