// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Complete-linkage hierarchical clustering.
//!
//! Native equivalent of `ML/CLUSTERING/CompleteLinkage.h`. See
//! `docs/CLUSTERING_SUPPORT.md`.

use super::cluster_functor::ClusterFunctor;
use super::linkage::lance_williams;
use crate::Result;
use crate::concept::progress_logger::ProgressReporter;
use crate::data_structures::binary_tree_node::BinaryTreeNode;
use crate::data_structures::distance_matrix::DistanceMatrix;

/// Complete-linkage cluster method (source `CompleteLinkage`).
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
pub struct CompleteLinkage;

impl CompleteLinkage {
    /// A complete-linkage functor (source default constructor).
    pub const fn new() -> Self {
        Self
    }
}

impl ClusterFunctor for CompleteLinkage {
    /// Clusters by complete linkage (source `operator()`).
    ///
    /// After merging clusters `i` and `j`, the distance to every other cluster
    /// `k` becomes the maximal distance between their elements, computed as the
    /// Lance-Williams update `0.5 d(i,k) + 0.5 d(j,k) + 0.5 |d(i,k) - d(j,k)|`
    /// in `f32`. (The source's code comment calls it the minimum distance; the
    /// formula and the header documentation give the maximum.)    ///
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
        lance_williams(original_distance, threshold, progress, |dik, djk, _, _| {
            0.5_f32 * dik + 0.5_f32 * djk + 0.5_f32 * (dik - djk).abs()
        })
    }
}
