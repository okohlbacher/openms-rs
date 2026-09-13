// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The centroided-peptide feature finder (`FEATUREFINDER/FeatureFinderAlgorithmPicked.h`),
//! with its helper structures and trace fitters.
//!
//! Registered ahead of its early-TOPP-bundle work package so that the package
//! never edits a module root. It exports nothing yet.

/// Seeds, mass traces and isotope patterns (`FeatureFinderAlgorithmPickedHelperStructs.h`).
pub mod helper_structs;
