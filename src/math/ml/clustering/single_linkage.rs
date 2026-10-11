// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Single-linkage hierarchical clustering with the SLINK algorithm.
//!
//! Native equivalent of `ML/CLUSTERING/SingleLinkage.h`. See
//! `docs/CLUSTERING_SUPPORT.md`.

use std::collections::BTreeSet;

use super::cluster_functor::{ClusterFunctor, at, at_mut, insufficient_input};
use crate::concept::progress_logger::ProgressReporter;
use crate::data_structures::binary_tree_node::BinaryTreeNode;
use crate::data_structures::distance_matrix::DistanceMatrix;
use crate::math::libstdcxx::std_min;
use crate::math::source_sort::source_sort_by;
use crate::{Error, Result};

/// Single-linkage cluster method (source `SingleLinkage`).
///
/// The details of the method can be found in "SLINK: An optimally efficient
/// algorithm for the single-link cluster method", *The Computer Journal*
/// 16(1):30-34, 1973, doi:10.1093/comjnl/16.1.30. After merging two clusters,
/// the distance between clusters is the minimal distance between their
/// elements.
///
/// The source also derives from `ProgressLogger`; progress goes to the
/// reporter passed to [`ClusterFunctor::cluster_with_progress`]. The source's
/// constructors, assignment and destructor map to `Default`, `Clone` and
/// `Copy`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SingleLinkage;

impl SingleLinkage {
    /// A single-linkage functor (source default constructor).
    pub const fn new() -> Self {
        Self
    }
}

impl ClusterFunctor for SingleLinkage {
    /// Clusters by single linkage (source `operator()`).
    ///
    /// `threshold` only exists to meet the base interface: the algorithm is
    /// fast and does not work by growing distances, so it always clusters
    /// completely. The distance matrix is only read.
    ///
    /// The pre-tree SLINK produces is sorted by distance with the Release
    /// build's `std::sort` ([`source_sort_by`]), so equal distances end up in
    /// the order the source leaves them in, then renumbered to lowest element
    /// indices.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] for fewer than two elements (source
    ///   `ClusterFunctor::InsufficientInput`).
    /// - [`Error::Unsupported`] for a `threshold` below 1: the source logs
    ///   "You tried to use Single Linkage clustering with a threshold. This is
    ///   currently not supported!" and throws `Exception::NotImplemented`. A
    ///   NaN threshold is not below 1 and clusters completely, as in the source.
    /// - [`Error::InvalidValue`] where the source would index outside its
    ///   vectors, which only a matrix whose distances are not ordered (NaN) can
    ///   provoke.
    fn cluster_with_progress(
        &self,
        original_distance: &mut DistanceMatrix<f32>,
        threshold: f32,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<Vec<BinaryTreeNode>> {
        let n = original_distance.dimension_size();
        if n < 2 {
            return Err(insufficient_input(
                "Distance matrix to start from only contains one element",
            ));
        }
        if threshold < 1.0 {
            return Err(Error::Unsupported(
                "You tried to use Single Linkage clustering with a threshold. This is currently \
                 not supported!"
                    .into(),
            ));
        }

        // SLINK pointer representation.
        let mut pi: Vec<usize> = Vec::with_capacity(n);
        let mut lambda: Vec<f32> = Vec::with_capacity(n);

        progress.start_count(n, "clustering data")?;

        pi.push(0);
        lambda.push(f32::MAX);

        for k in 1..n {
            pi.push(k);
            lambda.push(f32::MAX);
            let mut row_k = Vec::with_capacity(k);
            for i in 0..k {
                row_k.push(original_distance.get(i, k)?);
            }

            for i in 0..k {
                let pi_i = pi[i];
                if lambda[i] >= row_k[i] {
                    let updated = std_min(*at(&row_k, pi_i)?, lambda[i]);
                    *at_mut(&mut row_k, pi_i)? = updated;
                    lambda[i] = row_k[i];
                    pi[i] = k;
                } else {
                    let updated = std_min(*at(&row_k, pi_i)?, row_k[i]);
                    *at_mut(&mut row_k, pi_i)? = updated;
                }
            }

            for i in 0..k {
                if lambda[i] >= *at(&lambda, pi[i])? {
                    pi[i] = k;
                }
            }
            progress.set_count(k)?;
        }

        // Strict order is kept by the algorithm: i < pi[i].
        let mut tree: Vec<BinaryTreeNode> = (0..n - 1)
            .map(|i| BinaryTreeNode::new(i, pi[i], lambda[i]))
            .collect();

        source_sort_by(&mut tree, |x, y| x.distance < y.distance)?;

        // Convert the pre-tree to cluster-index form.
        for i in 0..tree.len() {
            if tree[i].right_child < tree[i].left_child {
                let node = &mut tree[i];
                std::mem::swap(&mut node.left_child, &mut node.right_child);
            }
            let (left, right) = (tree[i].left_child, tree[i].right_child);
            for node in tree.iter_mut().skip(i + 1) {
                if node.left_child == right {
                    node.left_child = left;
                } else if node.left_child > right {
                    node.left_child -= 1;
                }
                if node.right_child == right {
                    node.right_child = left;
                } else if node.right_child > right {
                    node.right_child -= 1;
                }
            }
        }

        // Redo the clustering to express every node by lowest element indices.
        let mut clusters: Vec<BTreeSet<usize>> = (0..n).map(|i| BTreeSet::from([i])).collect();
        for node in &mut tree {
            let (left, right) = (node.left_child, node.right_child);
            let new_left = first(at(&clusters, left)?)?;
            let new_right = first(at(&clusters, right)?)?;
            let merged = at(&clusters, right)?.clone();
            at_mut(&mut clusters, left)?.extend(merged);
            clusters.remove(right);
            node.left_child = new_left;
            node.right_child = new_right;
            if node.left_child > node.right_child {
                std::mem::swap(&mut node.left_child, &mut node.right_child);
            }
        }

        progress.end()?;
        Ok(tree)
    }
}

/// `*set.begin()`, which the source dereferences without checking.
pub(crate) fn first(set: &BTreeSet<usize>) -> Result<usize> {
    set.first()
        .copied()
        .ok_or_else(|| Error::InvalidValue("empty cluster".into()))
}
