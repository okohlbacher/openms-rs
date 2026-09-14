// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The region quadtree `FeatureOverlapFilter` queries, ported from the
//! source's bundled `extern/Quadtree` (`Quadtree.h`, `Box.h`).
//!
//! Registered ahead of its early-TOPP-bundle work package so that the package
//! never edits a module root. Package B9-OVERLAP fills it; it exports nothing
//! yet.
