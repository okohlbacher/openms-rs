// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Seed lists for feature detection.
//!
//! Ports `FEATUREFINDER/SeedListGenerator.h` and its `.cpp`. See
//! `docs/FEATUREFINDER_MODELS_SUPPORT.md`.
//!
//! Seed lists give the positions in an MS experiment where features are
//! expected. Of the feature finders only the "centroided" algorithm
//! (`FeatureFinderAlgorithmPicked`) takes custom seed lists, as featureXML.
//! The source marks the class experimental: "a new class that still needs to
//! be tested".

use crate::Result;
use crate::data_structures::DPosition2;
use crate::error::Error;
use crate::identification::PeptideIdentification;
use crate::kernel::MSExperiment;
use crate::kernel::features::{ConsensusMap, Feature, FeatureMap};
use std::collections::BTreeMap;

/// A list of seed positions `(RT, m/z)` (source `SeedList`,
/// `std::vector<DPosition<2>>`).
pub type SeedList = Vec<DPosition2>;

/// Most seeds one call may produce, and most consensus features times maps
/// one [`SeedListGenerator::generate_seed_lists`] may visit.
///
/// Native-only ceiling; the source has none.
pub const MAX_SEEDS: usize = 50_000_000;

fn seed(rt: f64, mz: f64) -> DPosition2 {
    DPosition2 {
        coordinates: [rt, mz],
    }
}

fn too_many() -> Error {
    Error::InvalidValue(format!("a seed list may hold at most {MAX_SEEDS} seeds"))
}

/// Generates seed lists for feature detection (source `SeedListGenerator`,
/// stateless).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SeedListGenerator;

impl SeedListGenerator {
    /// A generator (source default constructor).
    pub fn new() -> Self {
        Self
    }

    /// Seeds at the MS2 precursors of `experiment`: for every MS2 spectrum,
    /// the RT of its precursor spectrum and the m/z of its first precursor
    /// (source `generateSeedList(const PeakMap&, SeedList&)`).
    ///
    /// MS2 spectra without a precursor spectrum or without precursor
    /// information are skipped, as in the source. The precursor spectrum is
    /// [`MSExperiment::precursor_spectrum_index`]'s.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] beyond [`MAX_SEEDS`].
    pub fn generate_seed_list_from_experiment(
        &self,
        experiment: &MSExperiment,
    ) -> Result<SeedList> {
        let mut seeds = SeedList::new();
        for (index, spectrum) in experiment.spectra.iter().enumerate() {
            if spectrum.ms_level != 2 {
                continue;
            }
            let Some(precursor_index) = experiment.precursor_spectrum_index(index)? else {
                continue;
            };
            let Some(precursor) = spectrum.precursors.first() else {
                continue;
            };
            if seeds.len() >= MAX_SEEDS {
                return Err(too_many());
            }
            seeds.push(seed(experiment.spectra[precursor_index].rt, precursor.mz));
        }
        Ok(seeds)
    }

    /// Seeds at peptide identifications: each identification's RT with its
    /// precursor m/z, or - when `use_peptide_mass` is set and it has hits -
    /// the m/z of its best hit's sequence at that hit's charge (source
    /// `generateSeedList(PeptideIdentificationList&, SeedList&, bool)`).
    ///
    /// With `use_peptide_mass` the hits of every identification that has any
    /// are sorted by score, as in the source. An identification without an RT
    /// or m/z contributes NaN there, the value the source's accessors return
    /// for an unset one.
    ///
    /// # Errors
    ///
    /// Returns what sorting the hits and the sequence's m/z refuse, and
    /// [`Error::InvalidValue`] beyond [`MAX_SEEDS`]. On error `peptides` may
    /// have some identifications sorted already, as in the source.
    pub fn generate_seed_list_from_peptides(
        &self,
        peptides: &mut [PeptideIdentification],
        use_peptide_mass: bool,
    ) -> Result<SeedList> {
        if peptides.len() > MAX_SEEDS {
            return Err(too_many());
        }
        let mut seeds = SeedList::with_capacity(peptides.len());
        for pep in peptides.iter_mut() {
            let mz = if !pep.hits.is_empty() && use_peptide_mass {
                pep.sort()?;
                let hit = &pep.hits[0];
                hit.sequence.mz(hit.charge)?
            } else {
                pep.mz.unwrap_or(f64::NAN)
            };
            seeds.push(seed(pep.rt.unwrap_or(f64::NAN), mz));
        }
        Ok(seeds)
    }

    /// One seed list per constituent map of `consensus`, keyed by map index
    /// (source `generateSeedLists`).
    ///
    /// For every map, a seed is placed at each consensus feature that holds no
    /// sub-feature from that map, "to fill holes in the consensus map by
    /// looking for features explicitly at the positions of the holes". The
    /// seeds use the consensus map's RT scale, which after alignment may
    /// differ from the original experiments'.
    ///
    /// The source adds the seed to every map listed in the column headers and
    /// then removes the last seed of each map a sub-feature comes from. That
    /// is reproduced, with its consequences: a sub-feature of a map without a
    /// column header creates that map's list and (having nothing to remove)
    /// leaves it empty; two sub-features from one map remove two seeds, the
    /// second one belonging to an earlier consensus feature.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the work exceeds [`MAX_SEEDS`].
    pub fn generate_seed_lists(&self, consensus: &ConsensusMap) -> Result<BTreeMap<u64, SeedList>> {
        let maps = consensus.column_headers.len().max(1);
        if consensus.features.len().saturating_mul(maps) > MAX_SEEDS {
            return Err(too_many());
        }
        let mut seed_lists: BTreeMap<u64, SeedList> = BTreeMap::new();
        for cons in &consensus.features {
            let point = seed(cons.rt, cons.mz);
            for map_index in consensus.column_headers.keys() {
                seed_lists.entry(*map_index).or_default().push(point);
            }
            for handle in cons.handles() {
                // Source `seed_lists[map].pop_back()`: popping an empty vector
                // is undefined behaviour there; nothing is removed here.
                seed_lists.entry(handle.map_index).or_default().pop();
            }
        }
        Ok(seed_lists)
    }

    /// A feature map with one feature per seed at its RT and m/z, the unique
    /// IDs counting from 0 (source `convertSeedList(const SeedList&,
    /// FeatureMap&)`, which first clears the map including its metadata).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] beyond [`MAX_SEEDS`].
    pub fn seeds_to_features(&self, seeds: &[DPosition2]) -> Result<FeatureMap> {
        if seeds.len() > MAX_SEEDS {
            return Err(too_many());
        }
        let mut features = FeatureMap::default();
        features.features.reserve(seeds.len());
        for (counter, s) in seeds.iter().enumerate() {
            let mut feature = Feature::default();
            feature.base.rt = s.coordinates[0];
            feature.base.mz = s.coordinates[1];
            feature.base.unique_id = counter as u64;
            features.features.push(feature);
        }
        Ok(features)
    }

    /// The seed positions `(RT, m/z)` of a feature map's features (source
    /// `convertSeedList(const FeatureMap&, SeedList&)`).
    pub fn features_to_seeds(&self, features: &FeatureMap) -> SeedList {
        features
            .features
            .iter()
            .map(|f| seed(f.base.rt, f.base.mz))
            .collect()
    }
}
