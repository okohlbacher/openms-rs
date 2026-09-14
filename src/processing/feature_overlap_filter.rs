// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Removing and merging overlapping features
//! (`PROCESSING/FEATURE/FeatureOverlapFilter.h`).
//!
//! Registered ahead of its early-TOPP-bundle work package so that the package
//! never edits a module root. Package B9-OVERLAP fills it; it exports nothing
//! yet. The `quadtree` registration below and its doc line are
//! integrator-owned; the rest of the file is B9-OVERLAP's.

/// The region quadtree the filter queries (`extern/Quadtree`).
pub mod quadtree;
