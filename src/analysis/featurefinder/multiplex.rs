// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! FeatureFinderMultiplex: identification-free detection and quantification
//! of labelled peptide multiplets (SILAC, Leu, dimethyl, ICPL and numeric mass
//! shifts), and its building blocks.
//!
//! One submodule per source header, except that the two satellite headers
//! share [`satellite`]. See `docs/FEATUREFINDER_MULTIPLEX_SUPPORT.md` for the
//! API mapping and the evidence.

/// Mass shift patterns (source `MultiplexDeltaMasses.h`).
pub mod delta_masses;
/// All mass shifts a labelling can cause (source
/// `MultiplexDeltaMassesGenerator.h`).
pub mod delta_masses_generator;
/// Peptides that passed the filters for one pattern (source
/// `MultiplexFilteredMSExperiment.h`).
pub mod filtered_ms_experiment;
/// A peak that passed the filters (source `MultiplexFilteredPeak.h`).
pub mod filtered_peak;
/// Code common to both filters (source `MultiplexFiltering.h`).
pub mod filtering;
/// Filtering centroided data (source `MultiplexFilteringCentroided.h`).
pub mod filtering_centroided;
/// Isotopic peak patterns (source `MultiplexIsotopicPeakPattern.h`).
pub mod isotopic_peak_pattern;
/// Satellite peaks and data points (source `MultiplexSatelliteCentroided.h`
/// and `MultiplexSatelliteProfile.h`).
pub mod satellite;

pub use delta_masses::{DeltaMass, LabelSet, MultiplexDeltaMasses};
pub use delta_masses_generator::{LabelMassShifts, MultiplexDeltaMassesGenerator};
pub use filtered_ms_experiment::MultiplexFilteredMSExperiment;
pub use filtered_peak::MultiplexFilteredPeak;
pub use filtering::MultiplexFiltering;
pub use filtering_centroided::MultiplexFilteringCentroided;
pub use isotopic_peak_pattern::MultiplexIsotopicPeakPattern;
pub use satellite::{MultiplexSatelliteCentroided, MultiplexSatelliteProfile};
