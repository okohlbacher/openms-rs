// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `QC`: quality-control metrics computed from experiments, feature
//! maps and identifications.
//!
//! The source keeps these in their own top-level directory, `OpenMS/QC`. They
//! live under `analysis` here because the one file format that consumes them,
//! `MzQCFile` ([`crate::format::mzqc_file`]), sits in `format`, and `format`
//! and `analysis` already reach each other; a top-level `qc` module would add
//! a new cycle. The API mapping, preserved source behaviour and evidence are
//! in `docs/QC_SUPPORT.md`.
//!
//! Submodules are registered below, one per ported header. The other thirteen
//! headers of the source directory are not ported yet; the support document
//! lists them.

/// The shared base of every metric: input requirements, runnability and the
/// spectrum lookup map (`QC/QCBase.h`).
pub mod qc_base;
pub use qc_base::{QCBase, Requires, SpectraMap, Status, ToleranceUnit};

/// Number of spectra per MS level (`QC/SpectrumCount.h`).
pub mod spectrum_count;
pub use spectrum_count::SpectrumCount;

/// Total ion current trace with stability counters (`QC/TIC.h`).
pub mod tic;
pub use tic::{Tic, TicResult};

/// Feature count and mean retention-time shift (`QC/FeatureSummary.h`).
pub mod feature_summary;
pub use feature_summary::{FeatureSummary, FeatureSummaryResult};

/// Missed cleavages per peptide identification (`QC/MissedCleavages.h`).
pub mod missed_cleavages;
pub use missed_cleavages::MissedCleavages;

/// Counts and means over proteins and peptides (`QC/IdentificationSummary.h`).
pub mod identification_summary;
pub use identification_summary::{IdentificationSummary, IdentificationSummaryResult, UniqueId};
