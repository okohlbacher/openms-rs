// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `ANALYSIS/ID`: peptide and protein identification post-processing: consensus identification, score handling, FDR and Percolator support, search engines and metabolite identification. Older ID ports (`false_discovery_rate`, `peptide_indexing`, `id_filter`, …) predate this module and stay where they are.
//!
//! Submodules are registered below, one per ported header or cohesive group of
//! headers. Parallel porting lanes each append their own `pub mod` lines here.
/// Aho-Corasick matching with ambiguous residues and mismatches (`AhoCorasickAmbiguous.h`).
pub mod aho_corasick_ambiguous;
/// Phosphosite localisation with the AScore (`AScore.h`).
pub mod ascore;
/// Comet `variable_modXX` entries and their merging (`CometModification.h`).
pub mod comet_modification;
/// `ConsensusIDAlgorithm` and its eight variants.
pub mod consensus_id;
/// `ConsensusMapMergerAlgorithm`: merge the identification runs of a consensus map.
pub mod consensus_map_merger;
/// Fragment index for database search (`FragmentIndex.h`).
pub mod fragment_index;
/// `IDDecoyProbability`: target/decoy score-distribution probabilities.
pub mod id_decoy_probability;
/// Annotation of spectra, feature maps and consensus maps with peptide identifications (`IDMapper`).
pub mod id_mapper;
/// `IDMergerAlgorithm`: merge identification runs into one run.
pub mod id_merger;
/// `IDScoreSwitcherAlgorithm`: switch main scores to a meta-value score.
pub mod id_score_switcher;
/// Ion Identity Molecular Networking annotations and the GNPS supplementary-pairs table.
pub mod ion_identity_molecular_networking;
/// Spectral-library search for metabolite MS2 spectra (`MetaboliteSpectralMatching`).
pub mod metabolite_spectral_matching;
/// `NeedlemanWunsch` (`ANALYSIS/SEQUENCE`), used by the PEPMatrix consensus.
pub mod needleman_wunsch;
/// Subset-neighbor peptide search (`NeighborSeq.h`).
pub mod neighbor_seq;
/// `PeptideProteinResolution`: resolve shared peptides between protein groups.
pub mod peptide_protein_resolution;
/// `Percolator.h`: in-process semi-supervised target/decoy rescoring.
pub mod percolator;
/// `PercolatorFeatureSetHelper.h`: search-engine Percolator features and multi-engine run merging.
pub mod percolator_feature_set_helper;
/// `PercolatorTypes.h`: the in-process Percolator's input, output and model types.
pub mod percolator_types;
/// The Andromeda-style PScore and windowed peak ranks (`PScore.h`).
pub mod pscore;
/// Minimal in-memory search engine (`SimpleSearchEngineAlgorithm.h`).
pub mod simple_search_engine;
