// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! File summaries for the FileInfo tool (`FORMAT/FileInfo.h`).
//!
//! Registered ahead of its early-TOPP-bundle work package so that the package
//! never edits a module root. It exports nothing yet.

/// C++ stream and `StringUtils` numeric text formatting for the FileInfo report.
pub mod text_format;
