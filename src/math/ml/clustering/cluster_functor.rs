// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The interface every hierarchical cluster method implements.
//!
//! Native equivalent of `ML/CLUSTERING/ClusterFunctor.h`. See
//! `docs/CLUSTERING_SUPPORT.md`.

use crate::concept::progress_logger::ProgressReporter;
use crate::data_structures::binary_tree_node::BinaryTreeNode;
use crate::data_structures::distance_matrix::DistanceMatrix;
use crate::{Error, Result};

/// The source's default threshold: complete clustering until one cluster
/// remains.
pub const DEFAULT_CLUSTER_THRESHOLD: f32 = 1.0;

/// Base interface for cluster functors (source `ClusterFunctor`).
///
/// Each cluster functor merges clusters step by step, up to a given threshold,
/// starting from the most elementary partition of the data. Elements are
/// represented by the indices of a distance matrix, which should also
/// represent the order of the input.
///
/// The source's abstract `operator()` is [`cluster`](Self::cluster) (no
/// progress) and [`cluster_with_progress`](Self::cluster_with_progress). The
/// source returns the tree through an out-parameter it clears first; here the
/// tree is the return value. The source's default constructor, copy
/// constructor, assignment and destructor carry no state and map to the
/// implementors' `Default`/`Clone`.
pub trait ClusterFunctor {
    /// Clusters the indices of `original_distance` by their distances,
    /// reporting progress to `progress` as the source's `ProgressLogger` base
    /// does.
    ///
    /// `original_distance` holds the distances of the elements; it is
    /// considered mirrored at the main diagonal, so only one triangle is used.
    /// **It is changed during clustering**: keep a copy if it is needed again.
    ///
    /// The result represents the clustering: node `x` records clustering step
    /// `x`, the two merged clusters (by their lowest element index, not by
    /// cluster index, with `left_child < right_child`) and their distance.
    /// Use `ClusterAnalyzer::cut` to obtain the clusters of a given step.
    ///
    /// `threshold` is the smallest distance from which on merging is
    /// considered unrealistic; the source default is
    /// [`DEFAULT_CLUSTER_THRESHOLD`] (`1`), which clusters completely. It can
    /// be taken from the largest distance of two elements considered related,
    /// adapted to the clustering method.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the matrix holds fewer than two elements
    /// (source `ClusterFunctor::InsufficientInput`), and the implementor's own
    /// errors.
    fn cluster_with_progress(
        &self,
        original_distance: &mut DistanceMatrix<f32>,
        threshold: f32,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<Vec<BinaryTreeNode>>;

    /// [`cluster_with_progress`](Self::cluster_with_progress) without
    /// progress reporting (the source's default `ProgressLogger` state, which
    /// prints nothing).
    ///
    /// # Errors
    ///
    /// As [`cluster_with_progress`](Self::cluster_with_progress).
    fn cluster(
        &self,
        original_distance: &mut DistanceMatrix<f32>,
        threshold: f32,
    ) -> Result<Vec<BinaryTreeNode>> {
        self.cluster_with_progress(
            original_distance,
            threshold,
            &mut ProgressReporter::silent(),
        )
    }
}

/// The error the source raises as `ClusterFunctor::InsufficientInput`: fewer
/// than two data points, which no clustering algorithm can work with. The
/// source's default message is "not enough data points to cluster anything";
/// the linkage functors pass "Distance matrix to start from only contains one
/// element".
pub fn insufficient_input(message: &str) -> Error {
    Error::InvalidValue(format!("ClusterFunctor::InsufficentInput: {message}"))
}

/// `v[i]` or an error naming an index the source would have read out of
/// bounds.
pub(crate) fn at<T>(values: &[T], index: usize) -> Result<&T> {
    values
        .get(index)
        .ok_or_else(|| Error::InvalidValue(format!("cluster index {index} out of range")))
}

/// `&mut v[i]`, see [`at`].
pub(crate) fn at_mut<T>(values: &mut [T], index: usize) -> Result<&mut T> {
    values
        .get_mut(index)
        .ok_or_else(|| Error::InvalidValue(format!("cluster index {index} out of range")))
}
