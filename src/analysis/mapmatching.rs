// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `ANALYSIS/MAPMATCHING`: map alignment and feature grouping: superimposers, alignment algorithms, pair and cluster finders, and consensus-map normalisation. Transformation models live in [`crate::analysis::transformations`].
//!
//! Submodules are registered below, one per ported header or cohesive group of
//! headers. Parallel porting lanes each append their own `pub mod` lines here.
/// The interface of element group finders (source `BaseGroupFinder.h`).
pub mod base_group_finder;
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
/// The distance between two features that feature grouping minimises (source
/// `FeatureDistance.h`).
pub mod feature_distance;
/// The interface of feature grouping algorithms and their shared
/// post-processing (source `FeatureGroupingAlgorithm.h`).
pub mod feature_grouping_algorithm;
/// KD-tree feature grouping with optional LOWESS warp (source
/// `FeatureGroupingAlgorithmKD.h`, with `ClusterProxyKD`).
pub mod feature_grouping_algorithm_kd;
/// Labeled feature grouping with two labels (source
/// `FeatureGroupingAlgorithmLabeled.h`).
pub mod feature_grouping_algorithm_labeled;
/// QT feature grouping (source `FeatureGroupingAlgorithmQT.h`).
pub mod feature_grouping_algorithm_qt;
/// Unlabeled feature grouping by stepwise pairing (source
/// `FeatureGroupingAlgorithmUnlabeled.h`).
pub mod feature_grouping_algorithm_unlabeled;
/// Matching MS2 spectra to features by precursor position (source
/// `FeatureMapping.h`).
pub mod feature_mapping;
/// A feature as QT feature grouping stores it (source
/// `DATASTRUCTURES/GridFeature.h`).
pub mod grid_feature;
/// Features of several maps in one (RT, m/z) kd-tree (source
/// `ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h` and `KDTreeFeatureNode.h`).
pub mod kd_tree_feature_maps;
/// Pair finding for isotope-labeled features within one map (source
/// `LabeledPairFinder.h`).
pub mod labeled_pair_finder;
/// RT alignment from conflict-free connected components of a kd-tree
/// compatibility graph (source `MapAlignmentAlgorithmKD.h`).
pub mod map_alignment_algorithm_kd;
/// Parameter-tree helpers shared by the feature-grouping ports.
pub(crate) mod params;
/// A quality-threshold cluster for QT feature grouping (source
/// `DATASTRUCTURES/QTCluster.h`).
pub mod qt_cluster;
/// QT clustering of features across maps (source `QTClusterFinder.h`).
pub mod qt_cluster_finder;
/// Pair finding by mutual nearest neighbours (source `StablePairFinder.h`).
pub mod stable_pair_finder;
