// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `ANALYSIS/OPENSWATH` and `OPENSWATHALGO`: targeted and DIA (SWATH) analysis: transition lists, chromatogram extraction, scoring, feature finding and the OpenSWATH workflow.
//!
//! Submodules are registered below, one per ported header or cohesive group of
//! headers. Parallel porting lanes each append their own `pub mod` lines here.

/// `IncludeExcludeTarget`: an inclusion or exclusion target of a targeted experiment.
pub mod include_exclude_target;
/// `MetaboTargetedTargetDecoy`: resolve overlapping and missing decoys in targeted metabolomics assay libraries.
pub mod metabo_targeted_target_decoy;
/// `MRMAssay`: assay generation (annotation, restriction, detecting and UIS transitions) from raw targeted experiments.
pub mod mrm_assay;
/// `MRMDecoy`: reverse, pseudo-reverse and shuffle decoys for targeted assay libraries.
pub mod mrm_decoy;
/// `MRMIonSeries`: theoretical fragment ion series and transition annotation for assays and decoys.
pub mod mrm_ion_series;
/// `MRMMapping`: map targeted assays to chromatograms by precursor and product m/z.
pub mod mrm_mapping;
/// `ReactionMonitoringTransition`: one SRM/MRM transition of an assay library.
pub mod reaction_monitoring_transition;
/// `TargetedExperiment`: the in-memory TraML document (transitions, analytes, metadata).
pub mod targeted_experiment;
/// `TargetedExperimentHelper`: the TraML value types (CV, protein, peptide, compound, retention time, product ion).
pub mod targeted_experiment_helper;
