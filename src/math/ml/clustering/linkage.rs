// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The merge loop `CompleteLinkage` and `AverageLinkage` share.
//!
//! The two source `operator()` bodies are the same loop with a different
//! Lance-Williams update; this module holds that loop once, and each functor
//! passes its update.

use std::collections::BTreeSet;

use super::cluster_functor::{at, at_mut, insufficient_input};
use super::single_linkage::first;
use crate::Result;
use crate::concept::progress_logger::ProgressReporter;
use crate::data_structures::binary_tree_node::BinaryTreeNode;
use crate::data_structures::distance_matrix::DistanceMatrix;

/// Merges the closest pair of clusters while their distance is below
/// `threshold`, updating the matrix with `update(d_ik, d_jk, |i|, |j|)`.
///
/// Mirrors `CompleteLinkage::operator()` / `AverageLinkage::operator()`: the
/// closest pair `(i, j)` (`i > j`, from the tracked minimum) is recorded by
/// the lowest elements of both clusters, cluster `i` is merged into cluster
/// `j`, row and column `j` receive the updated distances, row and column `i`
/// are removed, and the minimum is recomputed. With two clusters left the last
/// pair is recorded without touching the matrix. Steps not reached because of
/// the threshold are filled with placeholder nodes `(first, lead, -1)`, where
/// `first` is the lowest element of the first remaining cluster, until the
/// tree has `n - 1` nodes.
///
/// The source bounds that fill by the tree vector's capacity, which is the
/// `n - 1` it reserves only when the caller's vector had no larger capacity
/// already; a larger one lets the source append one extra placeholder after a
/// complete clustering. This port always stops at `n - 1`.
pub(crate) fn lance_williams(
    original_distance: &mut DistanceMatrix<f32>,
    threshold: f32,
    progress: &mut ProgressReporter<'_>,
    update: impl Fn(f32, f32, usize, usize) -> f32,
) -> Result<Vec<BinaryTreeNode>> {
    let n = original_distance.dimension_size();
    if n < 2 {
        return Err(insufficient_input(
            "Distance matrix to start from only contains one element",
        ));
    }

    let mut clusters: Vec<BTreeSet<usize>> = (0..n).map(|i| BTreeSet::from([i])).collect();
    let capacity = n - 1;
    let mut tree: Vec<BinaryTreeNode> = Vec::with_capacity(capacity);

    original_distance.update_min_element()?;
    let mut min = original_distance.min_element_coordinates()?;

    progress.start_count(n, "clustering data")?;

    while original_distance.get(min.0, min.1)? < threshold {
        let mut node = BinaryTreeNode::new(
            first(at(&clusters, min.1)?)?,
            first(at(&clusters, min.0)?)?,
            original_distance.get(min.0, min.1)?,
        );
        if node.left_child > node.right_child {
            std::mem::swap(&mut node.left_child, &mut node.right_child);
        }
        tree.push(node);

        if original_distance.dimension_size() > 2 {
            let size_i = at(&clusters, min.0)?.len();
            let size_j = at(&clusters, min.1)?.len();

            let merged = at(&clusters, min.0)?.clone();
            at_mut(&mut clusters, min.1)?.extend(merged);
            clusters.remove(min.0);

            for k in 0..min.1 {
                let dik = original_distance.get(min.0, k)?;
                let djk = original_distance.get(min.1, k)?;
                original_distance.set_quick(min.1, k, update(dik, djk, size_i, size_j))?;
            }
            for k in min.1 + 1..original_distance.dimension_size() {
                let dik = original_distance.get(min.0, k)?;
                let djk = original_distance.get(min.1, k)?;
                original_distance.set_quick(k, min.1, update(dik, djk, size_i, size_j))?;
            }

            original_distance.reduce(min.0)?;
            original_distance.update_min_element()?;
            min = original_distance.min_element_coordinates()?;
        } else {
            break;
        }
        progress.set_count(n - original_distance.dimension_size())?;
    }

    let lead = first(at(&clusters, 0)?)?;
    for cluster in clusters.iter().skip(1) {
        if tree.len() >= capacity {
            break;
        }
        tree.push(BinaryTreeNode::new(lead, first(cluster)?, -1.0));
    }

    progress.end()?;
    Ok(tree)
}
