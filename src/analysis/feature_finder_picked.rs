// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The centroided-peptide feature finder (`FEATUREFINDER/FeatureFinderAlgorithmPicked.h`),
//! with its helper structures and trace fitters.
//!
//! This root is integrator-owned: its submodules are registered ahead of their
//! early-TOPP-bundle work packages, so that no package edits it. A submodule
//! that its package has not filled yet exports nothing.

/// Seeds, mass traces and isotope patterns (`FeatureFinderAlgorithmPickedHelperStructs.h`).
pub mod helper_structs;

/// The retention-time shape model shared by the trace fitters (`TraceFitter.h`).
pub mod trace_fitter;

/// The Gaussian retention-time model (`GaussTraceFitter.h`).
pub mod gauss_trace_fitter;

/// The exponential-Gaussian hybrid retention-time model (`EGHTraceFitter.h`).
pub mod egh_trace_fitter;

/// Parameters, input validation and the entry point (`FeatureFinderAlgorithmPicked.h`).
pub mod algorithm;

/// Intensity, trace and isotope-pattern scores (`FeatureFinderAlgorithmPicked.h`).
pub mod scoring;

/// Isotope-pattern precalculation and seed selection (`FeatureFinderAlgorithmPicked.h`).
pub mod seeds;
