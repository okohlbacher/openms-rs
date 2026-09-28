// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Hierarchical clustering with a pluggable similarity and cluster method.
//!
//! Native equivalent of `ML/CLUSTERING/ClusterHierarchical.h` (the
//! `ClusterHierarchical.cpp` translation unit only defines the
//! `UnnormalizedComparator` exception). The overload for binned peak spectra
//! needs the kernel and comparison types, which `math` may not depend on; it
//! is `cluster_binned_spectra` in the analysis module `spectra_clustering`. See
//! `docs/CLUSTERING_SUPPORT.md`.

use super::cluster_functor::ClusterFunctor;
use crate::data_structures::binary_tree_node::BinaryTreeNode;
use crate::data_structures::distance_matrix::DistanceMatrix;
use crate::{Error, Result};

/// A similarity value a comparator returns, converted to the stored `f32`
/// distance `1 - similarity` exactly as the source's expression
/// `1 - comparator(a, b)` and its implicit conversion to `float` do.
pub trait SimilarityValue: Copy {
    /// `1 - self`, evaluated in the type's own precision and then converted to
    /// `f32`.
    fn distance(self) -> f32;
}

impl SimilarityValue for f32 {
    fn distance(self) -> f32 {
        1.0 - self
    }
}

impl SimilarityValue for f64 {
    fn distance(self) -> f32 {
        (1.0 - self) as f32
    }
}

/// Hierarchical clustering with generic clustering functions (source
/// `ClusterHierarchical`).
///
/// Clusters objects with a distance (similarity) method and a cluster method.
/// The source's constructors and destructor map to `Default`, `Clone` and
/// `Copy`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClusterHierarchical {
    /// The threshold given to the cluster functor (source `threshold_`).
    threshold: f64,
}

impl Default for ClusterHierarchical {
    /// Threshold 1 (source default constructor).
    fn default() -> Self {
        Self { threshold: 1.0 }
    }
}

impl ClusterHierarchical {
    /// A clusterer with threshold 1 (source default constructor).
    pub const fn new() -> Self {
        Self { threshold: 1.0 }
    }

    /// Clusters `data` with `comparator` and `clusterer` (source template
    /// `cluster<Data, SimilarityComparator>`).
    ///
    /// When `original_distance` does not have one row per element of `data`
    /// (for example because it is empty), it is re-created: filled with 1 and
    /// then set, below the diagonal, to `1 - comparator(data[i], data[j])`
    /// without maintaining the tracked minimum (`setValueQuick`). A matrix of
    /// the right size is used as given. Clustering stops when `clusterer`
    /// reaches [`threshold`](Self::threshold), which is passed to it converted
    /// to `f32` as in the source. The result is the tree for further
    /// investigation with
    /// [`ClusterAnalyzer`](super::cluster_analyzer::ClusterAnalyzer).
    ///
    /// `comparator` must yield a similarity normalized to `[0, 1]` for every
    /// pair, for example
    /// [`EuclideanSimilarity`](super::euclidean_similarity::EuclideanSimilarity)
    /// with a suitable scale. Neither the source nor this port checks the
    /// range.
    ///
    /// # Errors
    ///
    /// The errors of `comparator` and of `clusterer`, and those of
    /// [`DistanceMatrix::resize`]. On a comparator error the matrix is left
    /// partially filled, as the source's exception leaves it.
    pub fn cluster<D, S: SimilarityValue>(
        &self,
        data: &[D],
        comparator: impl Fn(&D, &D) -> Result<S>,
        clusterer: &dyn ClusterFunctor,
        original_distance: &mut DistanceMatrix<f32>,
    ) -> Result<Vec<BinaryTreeNode>> {
        if original_distance.dimension_size() != data.len() {
            original_distance.clear();
            original_distance.resize(data.len(), 1.0)?;
            for i in 0..data.len() {
                for j in 0..i {
                    let distance = comparator(&data[i], &data[j])?.distance();
                    original_distance.set_quick(i, j, distance)?;
                }
            }
        }
        clusterer.cluster(original_distance, self.threshold as f32)
    }

    /// The threshold, in terms of distance (source `getThreshold`).
    pub fn threshold(&self) -> f64 {
        self.threshold
    }

    /// Sets the threshold, in terms of distance (source `setThreshold`).
    ///
    /// The default is 1, that is, the clustering stops only at similarity 0.
    /// Source warning: not every method supports a threshold yet;
    /// [`SingleLinkage`](super::single_linkage::SingleLinkage) refuses one
    /// below 1.
    pub fn set_threshold(&mut self, x: f64) {
        self.threshold = x;
    }
}

/// The error the source declares as `UnnormalizedComparator`: clustering was
/// requested with a compare functor whose values are not normalized to
/// `[0, 1]`, which the similarity-to-distance conversion requires. The source
/// declares the exception but never throws it; its default message is
/// "Clustering with unnormalized similarity measurement requested, normalized
/// is mandatory".
pub fn unnormalized_comparator(message: &str) -> Error {
    Error::InvalidValue(format!(
        "ClusterHierarchical::UnnormalizedComparator: {message}"
    ))
}
