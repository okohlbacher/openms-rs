// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Ports of `ML/CLUSTERING`: hierarchical clustering over a
//! [`DistanceMatrix`](crate::data_structures::distance_matrix::DistanceMatrix)
//! and grid-based clustering of points in the plane.
//!
//! See `docs/CLUSTERING_SUPPORT.md` for the API mapping and the evidence.

/// Average linkage (source `AverageLinkage.h`).
pub mod average_linkage;
/// Silhouette widths, Dunn indices and partitions of a clustering (source
/// `ClusterAnalyzer.h`).
pub mod cluster_analyzer;
/// The cluster functor interface (source `ClusterFunctor.h`).
pub mod cluster_functor;
/// Similarity-driven hierarchical clustering (source `ClusterHierarchical.h`).
pub mod cluster_hierarchical;
/// Complete linkage (source `CompleteLinkage.h`).
pub mod complete_linkage;
/// Euclidean similarity of 2D points (source `EuclideanSimilarity.h`).
pub mod euclidean_similarity;
mod linkage;
/// Single linkage by SLINK (source `SingleLinkage.h`).
pub mod single_linkage;
