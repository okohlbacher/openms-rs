// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `FEATUREFINDER`: feature detection beyond the picked algorithm: the 1-D models and fitters, FeatureFinderMultiplex, FeatureFinderIdentification, MetaboIdent and Biosaur2. The picked algorithm lives in [`crate::analysis::feature_finder_picked`].
//!
//! Submodules are registered below, one per ported header or cohesive group of
//! headers. Parallel porting lanes each append their own `pub mod` lines here.

/// The 1-D model interfaces: `BaseModel`, `BaseModel_impl`, `ModelDescription`, `InterpolationModel`.
pub mod base_model;
/// `ElutionModelFitter`: elution model fits and validation for features.
pub mod elution_model_fitter;
/// `EmgScoring`: elution peak scores from an EMG fit.
pub mod emg_scoring;
/// `FeatureFinderAlgorithm` and `FeatureFinderDefs`.
pub mod feature_finder_algorithm;
/// `FeatureFinderIdentificationAlgorithm`: ID-guided MS1 feature detection, the algorithm behind the TOPP tool `FeatureFinderIdentification`.
pub mod feature_finder_identification;
/// `FeatureFinderAlgorithmMetaboIdent`: targeted MS1 feature detection from a table of metabolite targets, the algorithm behind the TOPP tool `FeatureFinderMetaboIdent`.
pub mod feature_finder_metabo_ident;
/// `Fitter1D`, `LevMarqFitter1D` and `MaxLikeliFitter1D`.
pub mod fitter1d;
/// The five 1-D fitters: Gauss, bi-Gauss, EMG, isotope and extended isotope.
pub mod fitters;
/// `GaussModel`, `BiGaussModel` and `EmgModel`.
pub mod gauss_models;
/// The reference build's single-precision `logf`.
pub(crate) mod glibc_logf;
/// `IsotopeModel` and `ExtendedIsotopeModel`.
pub mod isotope_model;
/// The steps the two targeted MS1 feature finders share (not public API).
pub(crate) mod ms1_targeted;
/// `PeakWidthEstimator`: peak width by m/z from picked peaks.
pub mod peak_width_estimator;
/// `SeedListGenerator`: seed lists for feature detection.
pub mod seed_list_generator;
