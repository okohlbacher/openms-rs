// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `FEATUREFINDER`: feature detection beyond the picked algorithm: the 1-D models and fitters, FeatureFinderMultiplex, FeatureFinderIdentification, MetaboIdent and Biosaur2. The picked algorithm lives in [`crate::analysis::feature_finder_picked`].
//!
//! Submodules are registered below, one per ported header or cohesive group of
//! headers. Parallel porting lanes each append their own `pub mod` lines here.

pub mod base_model;
pub mod emg_scoring;
pub mod feature_finder_algorithm;
pub mod fitter1d;
pub mod fitters;
pub mod gauss_models;
pub mod isotope_model;
pub mod peak_width_estimator;
pub mod seed_list_generator;
