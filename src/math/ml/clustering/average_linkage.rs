// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Average-linkage hierarchical clustering.
//!
//! Native equivalent of `ML/CLUSTERING/AverageLinkage.h`. See
//! `docs/CLUSTERING_SUPPORT.md`.

use super::cluster_functor::ClusterFunctor;
use super::linkage::lance_williams;
use crate::Result;
use crate::concept::progress_logger::ProgressReporter;
use crate::data_structures::binary_tree_node::BinaryTreeNode;
use crate::data_structures::distance_matrix::DistanceMatrix;

/// Average-linkage cluster method (source `AverageLinkage`).
///
/// The details of the method can be found in Backhaus, Erichson, Plinke,
/// Weiber, *Multivariate Analysemethoden*, Springer 2000, and in E. M.
/// Voorhees, "Implementing agglomerative hierarchic clustering algorithms for
/// use in document retrieval", *Inf. Process. Manage.* 22(6):465-476, 1986.
///
/// The source also derives from `ProgressLogger`; progress goes to the
/// reporter passed to [`ClusterFunctor::cluster_with_progress`]. The source's
/// constructors, assignment and destructor map to `Default`, `Clone` and
/// `Copy`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct AverageLinkage;

impl AverageLinkage {
    /// An average-linkage functor (source default constructor).
    pub const fn new() -> Self {
        Self
    }
}

impl ClusterFunctor for AverageLinkage {
    /// Clusters by average linkage (source `operator()`).
    ///
    /// After merging clusters `i` and `j` with `m_i` and `m_j` elements, the
    /// distance to every other cluster `k` becomes the average distance between
    /// their elements, computed as the Lance-Williams update
    /// `m_i/(m_i+m_j) d(i,k) + m_j/(m_i+m_j) d(j,k)`, both weights in `f32` as
    /// the source computes them.    ///
    /// Merging stops as soon as the closest pair of clusters is not closer
    /// than `threshold` (source default 1, complete clustering). The steps not
    /// taken are filled with placeholder nodes of distance `-1`, so the result
    /// always has `n - 1` nodes. The source documents the placeholders'
    /// children as `(0, 1)`; they are in fact the lowest element of the first
    /// remaining cluster and the lowest element of each further remaining
    /// cluster, which is what the pinned class test expects and what this port
    /// produces. Each node records the merged clusters by their lowest element
    /// indices with `left_child < right_child`.
    ///
    /// `original_distance` is changed during the clustering; keep a copy or be
    /// able to redo it.
    ///
    /// # Errors
    ///
    /// [`crate::Error::InvalidValue`] for fewer than two elements (source
    /// `ClusterFunctor::InsufficientInput`), and the errors of the
    /// [`DistanceMatrix`] operations, which a well-formed matrix never raises.
    fn cluster_with_progress(
        &self,
        original_distance: &mut DistanceMatrix<f32>,
        threshold: f32,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<Vec<BinaryTreeNode>> {
        lance_williams(
            original_distance,
            threshold,
            progress,
            |dik, djk, size_i, size_j| {
                // Source: `(float)(size_i / (float)(size_i + size_j))`.
                let total = (size_i + size_j) as f32;
                let alpha_i = size_i as f32 / total;
                let alpha_j = size_j as f32 / total;
                alpha_i * dik + alpha_j * djk
            },
        )
    }
}
