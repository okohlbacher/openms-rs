// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! File summaries for the FileInfo tool (`FORMAT/FileInfo.h`).
//!
//! This root is integrator-owned: its submodules are registered ahead of their
//! early-TOPP-bundle work packages, so that no package edits it. A submodule
//! that its package has not filled yet exports nothing.

/// C++ stream and `StringUtils` numeric text formatting for the FileInfo report.
pub mod text_format;

/// The FileInfo options and structured result.
pub mod model;

/// The FileInfo report: section order and text and TSV rendering.
pub mod report;

/// The FileInfo summary of peak files: DTA, DTA2D and mzML.
pub mod peaks;

/// The FileInfo summary of featureXML feature maps.
pub mod features;

/// The FileInfo checks and listings of `-i`, `-d` and `-c`.
pub mod checks;
