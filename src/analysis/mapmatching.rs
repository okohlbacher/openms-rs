// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `ANALYSIS/MAPMATCHING`: map alignment and feature grouping: superimposers, alignment algorithms, pair and cluster finders, and consensus-map normalisation. Transformation models live in [`crate::analysis::transformations`].
//!
//! Submodules are registered below, one per ported header or cohesive group of
//! headers. Parallel porting lanes each append their own `pub mod` lines here.

/// Alignment precision and recall against a ground truth
/// (`MapAlignmentEvaluationAlgorithm*`).
pub mod alignment_evaluation;
/// Alignment based on peptide identifications
/// (`MapAlignmentAlgorithmIdentification`).
pub mod alignment_identification;
/// Alignment based on pose clustering (`MapAlignmentAlgorithmPoseClustering`).
pub mod alignment_pose_clustering;
/// Affine superimposer by pose clustering (`PoseClusteringAffineSuperimposer`).
pub mod pose_clustering_affine;
/// Shift superimposer by pose clustering (`PoseClusteringShiftSuperimposer`).
pub mod pose_clustering_shift;
/// The superimposer interface (`BaseSuperimposer`) and the shared voting
/// histograms.
pub mod superimposer;
