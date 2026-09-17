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

/// Index pairs, index sets, the used flag and `NoSuccessor` (`FeatureFinderDefs` in `FeatureFinderAlgorithmPicked.h`).
pub mod defs;

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

/// Isotope fit and mass-trace extension of a seed (`FeatureFinderAlgorithmPicked.h`).
pub mod extension;

/// Trace fitting, cropping, quality checks and feature creation (`FeatureFinderAlgorithmPicked.h`).
pub mod fitting;

/// Overlap resolution and apex annotation (`FeatureFinderAlgorithmPicked.h`).
pub mod resolution;

/// The C++ Release build's `std::sort` and `std::stable_sort` orders (`FeatureFinderAlgorithmPicked.h`).
pub mod source_sort;

/// The reference build's C library `powf` of the overall seed score (`FeatureFinderAlgorithmPicked.h`).
pub(crate) mod glibc_powf;

/// The reference build's C library `exp`, `log`, `atan` and `sqrt` of the trace fitters (`FeatureFinderAlgorithmPicked.h`).
pub(crate) mod glibc_libm;

/// The debug mode: the source's `debug/` output as data (`FeatureFinderAlgorithmPicked.h`).
pub mod debug;

/// The stateful algorithm instance: reuse, a caller's map, parameters and progress (`FeatureFinderAlgorithmPicked.h`).
pub mod instance;
