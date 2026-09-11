// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Ported TOPP tools.
//!
//! Each tool is a library type implementing [`Tool`](super::Tool), so the
//! `src/bin/` executable and the differential test in `tests/` share exactly
//! one definition of its parameters and behaviour.

mod dta_extractor;
mod mzml_splitter;

pub use dta_extractor::DTAExtractor;
pub use mzml_splitter::MzMLSplitter;
