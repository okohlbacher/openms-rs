// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Ported TOPP tools.
//!
//! Each tool is a library type implementing [`Tool`](crate::cli::Tool), so the
//! `src/bin/` executable and the differential test in `tests/` share exactly
//! one definition of its parameters and behaviour.
//!
//! A tool whose executable needs more than `paramxml` exists under the same
//! features its `[[bin]]` entry requires: `PeakPickerHiRes` under `mzml`, and
//! `FileInfo` and `FeatureFinderCentroided` under `mzml` and `featurexml`. The
//! five earlier tools predate that rule and build under `paramxml` alone.

mod baseline_filter;
mod dta_extractor;
#[cfg(all(feature = "mzml", feature = "featurexml"))]
mod feature_finder_centroided;
#[cfg(all(feature = "mzml", feature = "featurexml"))]
mod file_info;
mod fuzzy_diff;
mod map_normalizer;
mod mzml_splitter;
#[cfg(feature = "mzml")]
mod peak_picker_hi_res;
mod spectra_filter_window_mower;

pub use baseline_filter::BaselineFilter;
pub use dta_extractor::DTAExtractor;
#[cfg(all(feature = "mzml", feature = "featurexml"))]
pub use feature_finder_centroided::FeatureFinderCentroided;
#[cfg(all(feature = "mzml", feature = "featurexml"))]
pub use file_info::FileInfo;
pub use fuzzy_diff::FuzzyDiff;
pub use map_normalizer::MapNormalizer;
pub use mzml_splitter::MzMLSplitter;
#[cfg(feature = "mzml")]
pub use peak_picker_hi_res::PeakPickerHiRes;
pub use spectra_filter_window_mower::SpectraFilterWindowMower;
