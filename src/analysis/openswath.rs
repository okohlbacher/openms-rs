// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `ANALYSIS/OPENSWATH` and `OPENSWATHALGO`: targeted and DIA (SWATH) analysis: transition lists, chromatogram extraction, scoring, feature finding and the OpenSWATH workflow.
//!
//! Submodules are registered below, one per ported header or cohesive group of
//! headers. Parallel porting lanes each append their own `pub mod` lines here.

/// `TransitionPQPFile`: read and write PQP (SQLite) transition libraries.
#[cfg(feature = "sqlite")]
pub mod transition_pqp_file;
/// `TransitionTSVFile`: read and write OpenSWATH transition lists (TSV/CSV/MRM).
pub mod transition_tsv_file;
/// Conversions between kernel and light records, from `ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h`.
pub mod data_access_helper;
/// Named numeric table writers, from `OPENSWATHALGO/DATAACCESS/DataFrameWriter.h`.
pub mod data_frame_writer;
/// Light spectra, chromatograms and binary data arrays, from `OPENSWATHALGO/DATAACCESS/DataStructures.h`.
pub mod data_structures;
/// `IncludeExcludeTarget`: an inclusion or exclusion target of a targeted experiment.
pub mod include_exclude_target;
/// Feature, peak-group, transition-group and S/N access traits, from `OPENSWATHALGO/DATAACCESS/ITransition.h`.
pub mod itransition;
/// Shared p-value, q-value and PEP inference for peptide, protein and gene rows, from `ANALYSIS/OPENSWATH/LevelContextInference.h`.
pub mod level_context_inference;
/// `MetaboTargetedTargetDecoy`: resolve overlapping and missing decoys in targeted metabolomics assay libraries.
pub mod metabo_targeted_target_decoy;
/// Test doubles for the access traits, from `OPENSWATHALGO/DATAACCESS/MockObjects.h`.
pub mod mock_objects;
/// Access-trait implementations over kernel types, from `ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h`.
pub mod mrm_feature_access;
/// `MRMMapping`: map targeted assays to chromatograms by precursor and product m/z.
pub mod mrm_mapping;
/// Gene-level context inference, from `ANALYSIS/OPENSWATH/OpenSwathGeneInference.h`.
pub mod openswath_gene_inference;
/// Helpers shared by the OpenSWATH tools, from `ANALYSIS/OPENSWATH/OpenSwathHelper.h`.
pub mod openswath_helper;
/// Error-estimation, IPF and level-context inference configuration, from `ANALYSIS/OPENSWATH/OpenSwathInferenceConfig.h`.
pub mod openswath_inference_config;
/// Inference levels, contexts and the compact inference rows, from `ANALYSIS/OPENSWATH/OpenSwathInferenceData.h`.
pub mod openswath_inference_data;
/// Canonical precursor and transition identifiers for OpenSWATH libraries, from `ANALYSIS/OPENSWATH/OpenSwathLibraryIDNormalizer.h`.
pub mod openswath_library_id_normalizer;
/// Peptide-level context inference, from `ANALYSIS/OPENSWATH/OpenSwathPeptideInference.h`.
pub mod openswath_peptide_inference;
/// Inference of peptidoforms (IPF), from `ANALYSIS/OPENSWATH/OpenSwathPeptidoformInference.h`.
pub mod openswath_peptidoform_inference;
/// Protein-level context inference, from `ANALYSIS/OPENSWATH/OpenSwathProteinInference.h`.
pub mod openswath_protein_inference;
/// `ReactionMonitoringTransition`: one SRM/MRM transition of an assay library.
pub mod reaction_monitoring_transition;
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
/// `TargetedExperiment`: the in-memory TraML document (transitions, analytes, metadata).
pub mod targeted_experiment;
/// `TargetedExperimentHelper`: the TraML value types (CV, protein, peptide, compound, retention time, product ion).
pub mod targeted_experiment_helper;
/// The light targeted experiment, from `OPENSWATHALGO/DATAACCESS/TransitionExperiment.h`.
pub mod transition_experiment;
/// Transition grouping and compound lookup, from `OPENSWATHALGO/DATAACCESS/TransitionHelper.h`.
pub mod transition_helper;
/// The nested targeted experiment, from `OPENSWATHALGO/DATAACCESS/Transitions.h`.
pub mod transitions;
/// `TransitionParquetFile`: read and write OpenSWATH Parquet assay libraries (`.oswpq`).
#[cfg(all(feature = "arrow", feature = "file-compression"))]
pub mod transition_parquet_file;
