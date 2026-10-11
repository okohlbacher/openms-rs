// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `ANALYSIS/OPENSWATH` and `OPENSWATHALGO`: targeted and DIA (SWATH) analysis: transition lists, chromatogram extraction, scoring, feature finding and the OpenSWATH workflow.
//!
//! Submodules are registered below, one per ported header or cohesive group of
//! headers. Parallel porting lanes each append their own `pub mod` lines here.

/// Confidence scoring of features against an assay library, from `ANALYSIS/OPENSWATH/ConfidenceScoring.h`.
pub mod confidence_scoring;
/// Coordinates, extraction and annotated chromatograms for targeted experiments, from `ANALYSIS/OPENSWATH/ChromatogramExtractor.h`.
pub mod chromatogram_extractor;
/// Chromatogram extraction at given coordinates, from `ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h`.
pub mod chromatogram_extractor_algorithm;
/// Conversions between kernel and light records, from `ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h`.
pub mod data_access_helper;
/// Named numeric table writers, from `OPENSWATHALGO/DATAACCESS/DataFrameWriter.h`.
pub mod data_frame_writer;
/// Light spectra, chromatograms and binary data arrays, from `OPENSWATHALGO/DATAACCESS/DataStructures.h`.
pub mod data_structures;
/// Window integration and averagine helpers for DIA scoring, from `ANALYSIS/OPENSWATH/DIAHelper.h`.
pub mod dia_helper;
/// Prescoring of a spectrum against library intensities, from `ANALYSIS/OPENSWATH/DIAPrescoring.h`.
pub mod dia_prescoring;
/// Full-spectrum DIA scores (isotopes, mass error, b/y ions), from `ANALYSIS/OPENSWATH/DIAScoring.h`.
pub mod dia_scoring;
/// `IncludeExcludeTarget`: an inclusion or exclusion target of a targeted experiment.
pub mod include_exclude_target;
/// Ion mobility scoring of fragment and precursor mobilograms, from `ANALYSIS/OPENSWATH/IonMobilityScoring.h`.
pub mod ion_mobility_scoring;
/// Feature, peak-group, transition-group and S/N access traits, from `OPENSWATHALGO/DATAACCESS/ITransition.h`.
pub mod itransition;
/// Correlation of mass traces into pseudo-spectra, from `ANALYSIS/OPENSWATH/MasstraceCorrelator.h`.
pub mod masstrace_correlator;
/// `MetaboTargetedTargetDecoy`: resolve overlapping and missing decoys in targeted metabolomics assay libraries.
pub mod metabo_targeted_target_decoy;
/// Test doubles for the access traits, from `OPENSWATHALGO/DATAACCESS/MockObjects.h`.
pub mod mock_objects;
/// Iterative feature selection, from `ANALYSIS/OPENSWATH/MRMBatchFeatureSelector.h`.
pub mod mrm_batch_feature_selector;
/// Access-trait implementations over kernel types, from `ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h`.
pub mod mrm_feature_access;
/// Flagging and filtering of MRM features against QC bounds, from `ANALYSIS/OPENSWATH/MRMFeatureFilter.h`.
pub mod mrm_feature_filter;
/// Finding and scoring of MRM features, from `ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h`.
pub mod mrm_feature_finder_scoring;
/// Per-component and per-component-group picker parameters, from `ANALYSIS/OPENSWATH/MRMFeaturePicker.h`.
pub mod mrm_feature_picker;
/// Quality-control bounds for MRM components and component groups, from `ANALYSIS/OPENSWATH/MRMFeatureQC.h`.
pub mod mrm_feature_qc;
/// Selection of MRM features by linear programming, from `ANALYSIS/OPENSWATH/MRMFeatureSelector.h`.
pub mod mrm_feature_selector;
/// Metavalue conversions shared by the MRM feature modules.
pub(crate) mod mrm_feature_support;
/// `MRMMapping`: map targeted assays to chromatograms by precursor and product m/z.
pub mod mrm_mapping;
/// Cross-correlation, library, S/N and mutual-information scores of SRM/MRM peak groups, from `ANALYSIS/OPENSWATH/MRMScoring.h`.
pub mod mrm_scoring;
/// Peak-group picking across the chromatograms of one precursor, from `ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h`.
pub mod mrm_transition_group_picker;
/// Helpers shared by the OpenSWATH tools, from `ANALYSIS/OPENSWATH/OpenSwathHelper.h`.
pub mod openswath_helper;
/// `OpenSwath_Scores_Usage` and the LDA scoring methods of `OpenSwath_Scores`, from `ANALYSIS/OPENSWATH/OpenSwathScores.h`.
pub mod openswath_scores;
/// The driver calling the chromatographic, library, DIA and ion mobility scores, from `ANALYSIS/OPENSWATH/OpenSwathScoring.h`.
pub mod openswath_scoring;
/// Peak picking in a single ion mobilogram, from `ANALYSIS/OPENSWATH/PeakPickerMobilogram.h`.
pub mod peak_picker_mobilogram;
/// `ReactionMonitoringTransition`: one SRM/MRM transition of an assay library.
pub mod reaction_monitoring_transition;
/// Raw m/z-RT-ion-mobility peak clouds at targeted coordinates, from `ANALYSIS/OPENSWATH/PeakMapExtractor.h`.
pub mod peak_map_extractor;
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
/// m/z and ion-mobility correction of SWATH maps (partial), from `ANALYSIS/OPENSWATH/SwathMapMassCorrection.h`.
pub mod swath_map_mass_correction;
/// OpenSWATH quality control (MS1 charge distribution), from `ANALYSIS/OPENSWATH/SwathQC.h`.
pub mod swath_qc;
/// SWATH-window files and SWATH-map annotation, from `ANALYSIS/OPENSWATH/SwathWindowLoader.h`.
pub mod swath_window_loader;
/// The light targeted experiment, from `OPENSWATHALGO/DATAACCESS/TransitionExperiment.h`.
pub mod transition_experiment;
/// Transition grouping and compound lookup, from `OPENSWATHALGO/DATAACCESS/TransitionHelper.h`.
pub mod transition_helper;
/// The nested targeted experiment, from `OPENSWATHALGO/DATAACCESS/Transitions.h`.
pub mod transitions;
