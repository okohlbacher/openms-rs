// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `ANALYSIS/OPENSWATH` and `OPENSWATHALGO`: targeted and DIA (SWATH) analysis: transition lists, chromatogram extraction, scoring, feature finding and the OpenSWATH workflow.
//!
//! Submodules are registered below, one per ported header or cohesive group of
//! headers. Parallel porting lanes each append their own `pub mod` lines here.

/// Conversions between kernel and light records, from `ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h`.
pub mod data_access_helper;
/// Named numeric table writers, from `OPENSWATHALGO/DATAACCESS/DataFrameWriter.h`.
pub mod data_frame_writer;
/// Light spectra, chromatograms and binary data arrays, from `OPENSWATHALGO/DATAACCESS/DataStructures.h`.
pub mod data_structures;
/// Feature, peak-group, transition-group and S/N access traits, from `OPENSWATHALGO/DATAACCESS/ITransition.h`.
pub mod itransition;
/// Test doubles for the access traits, from `OPENSWATHALGO/DATAACCESS/MockObjects.h`.
pub mod mock_objects;
/// Access-trait implementations over kernel types, from `ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h`.
pub mod mrm_feature_access;
/// Helpers shared by the OpenSWATH tools, from `ANALYSIS/OPENSWATH/OpenSwathHelper.h`.
pub mod openswath_helper;
/// Distances, cross-correlation, ranks and mutual information, from `OPENSWATHALGO/ALGO/Scoring.h`.
pub mod scoring;
/// The spectrum access interface, from `OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h`.
pub mod spectrum_access;
/// Spectrum access over OpenMS data, from `ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS*.h`, `SpectrumAccess*Transforming.h` and `SimpleOpenMSSpectraAccessFactory.h`.
pub mod spectrum_access_openms;
/// Spectrum access over an sqMass file, from `ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.h`.
#[cfg(feature = "sqmass")]
pub mod spectrum_access_sqmass;
/// Adding up spectra by resampling, from `ANALYSIS/OPENSWATH/SpectrumAddition.h`.
pub mod spectrum_addition;
/// Vector statistics for the scores, from `OPENSWATHALGO/ALGO/StatsHelpers.h`.
pub mod stats_helpers;
/// One SWATH map and its isolation window, from `OPENSWATHALGO/DATAACCESS/SwathMap.h`.
pub mod swath_map;
/// The light targeted experiment, from `OPENSWATHALGO/DATAACCESS/TransitionExperiment.h`.
pub mod transition_experiment;
/// Transition grouping and compound lookup, from `OPENSWATHALGO/DATAACCESS/TransitionHelper.h`.
pub mod transition_helper;
/// The nested targeted experiment, from `OPENSWATHALGO/DATAACCESS/Transitions.h`.
pub mod transitions;
/// Chromatogram extraction at given coordinates, from `ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h`.
pub mod chromatogram_extractor_algorithm;
/// SWATH-window files and SWATH-map annotation, from `ANALYSIS/OPENSWATH/SwathWindowLoader.h`.
pub mod swath_window_loader;
