// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `ANALYSIS/MAPMATCHING`: map alignment and feature grouping: superimposers, alignment algorithms, pair and cluster finders, and consensus-map normalisation. Transformation models live in [`crate::analysis::transformations`].
//!
//! Submodules are registered below, one per ported header or cohesive group of
//! headers. Parallel porting lanes each append their own `pub mod` lines here.
/// Median normalisation of consensus maps and the accession/description
/// filter shared with the threshold normaliser (source
/// `ConsensusMapNormalizerAlgorithmMedian.h`).
pub mod consensus_map_normalizer_median;
/// Quantile normalisation of consensus maps (source
/// `ConsensusMapNormalizerAlgorithmQuantile.h`).
pub mod consensus_map_normalizer_quantile;
/// Ratio-threshold normalisation of consensus maps (source
/// `ConsensusMapNormalizerAlgorithmThreshold.h`).
pub mod consensus_map_normalizer_threshold;
/// A feature as QT feature grouping stores it (source
/// `DATASTRUCTURES/GridFeature.h`).
pub mod grid_feature;
/// A quality-threshold cluster for QT feature grouping (source
/// `DATASTRUCTURES/QTCluster.h`).
pub mod qt_cluster;
/// Parameter-tree helpers shared by the feature-grouping ports.
pub(crate) mod params;
/// The distance between two features that feature grouping minimises (source
/// `FeatureDistance.h`).
pub mod feature_distance;
