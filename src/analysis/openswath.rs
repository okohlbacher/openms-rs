// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `ANALYSIS/OPENSWATH` and `OPENSWATHALGO`: targeted and DIA (SWATH) analysis: transition lists, chromatogram extraction, scoring, feature finding and the OpenSWATH workflow.
//!
//! Submodules are registered below, one per ported header or cohesive group of
//! headers. Parallel porting lanes each append their own `pub mod` lines here.

pub mod data_frame_writer;
pub mod data_structures;
pub mod itransition;
pub mod mock_objects;
pub mod scoring;
pub mod spectrum_access;
pub mod spectrum_access_openms;
pub mod stats_helpers;
pub mod swath_map;
pub mod transition_experiment;
pub mod transition_helper;
pub mod transitions;
