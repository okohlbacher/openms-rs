// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! One merge step of a hierarchical clustering.
//!
//! Native equivalent of `DATASTRUCTURES/BinaryTreeNode.h`. See
//! `docs/CLUSTERING_SUPPORT.md`.

/// Element of a binary tree representing a hierarchical clustering process.
///
/// Strict indexing is assumed: node number `x` of a clustering is clustering
/// step number `x`. `left_child` and `right_child` are each the lowest element
/// index of the two clusters merged in that step, and `distance` is the
/// distance between them. The linkage functors keep `left_child <
/// right_child`, and fill the steps a distance threshold prevented with
/// placeholder nodes whose distance is `-1`.
///
/// The source's constructor, copy constructor, assignment and destructor are
/// [`BinaryTreeNode::new`], `Clone`/`Copy` and ordinary drop. The source's
/// default constructor is private and not ported.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryTreeNode {
    /// Lowest element index of the first merged cluster.
    pub left_child: usize,
    /// Lowest element index of the second merged cluster.
    pub right_child: usize,
    /// Distance of the two merged clusters, `-1` for a placeholder step.
    pub distance: f32,
}

impl BinaryTreeNode {
    /// A node merging the clusters led by `i` and `j` at distance `x` (source
    /// `BinaryTreeNode(i, j, x)`). The children are stored as given, not
    /// ordered.
    pub const fn new(i: usize, j: usize, x: f32) -> Self {
        Self {
            left_child: i,
            right_child: j,
            distance: x,
        }
    }
}
