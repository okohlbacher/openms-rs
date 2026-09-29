// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Feature grouping for unlabeled data by stepwise pairing against a growing
//! consensus map.
//!
//! Native equivalent of
//! `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h`. See
//! `docs/FEATURE_GROUPING_SUPPORT.md`.
//!
//! The corresponding features must be aligned, but may deviate slightly in
//! position. There are two ways to run the algorithm:
//!
//! - call [`FeatureGroupingAlgorithm::group`] with all maps in memory; or
//! - call [`FeatureGroupingAlgorithmUnlabeled::set_reference`], then
//!   [`FeatureGroupingAlgorithmUnlabeled::add_to_group`] once per further map,
//!   then [`FeatureGroupingAlgorithmUnlabeled::result_map`]. Only the
//!   reference and the current map need to be in memory at a time.
//!
//! The two are **not** equivalent: `group` post-processes the result
//! (protein identifications transferred, unassigned peptide identifications
//! re-indexed, canonical sorting — see
//! [`postprocess`](super::feature_grouping_algorithm::postprocess)); the
//! incremental path returns the raw pair-finder result, and callers must
//! post-process it themselves if needed, as the `FeatureLinkerUnlabeled` tool
//! does.

use super::feature_grouping_algorithm::{
    FeatureGroupingAlgorithm, feature_map_identifications, postprocess,
};
use super::stable_pair_finder::StablePairFinder;
use crate::concept::UniqueIdGenerator;
use crate::error::{Error, Result};
use crate::kernel::MSExperiment;
use crate::kernel::conversion_helper::MapConversion;
use crate::kernel::features::{ConsensusMap, FeatureMap};
use crate::param::{DefaultParamHandler, Param};

/// Parameter-handler name of the source class.
pub const FEATURE_GROUPING_ALGORITHM_UNLABELED_NAME: &str = "FeatureGroupingAlgorithmUnlabeled";

/// Unlabeled feature grouping (source `FeatureGroupingAlgorithmUnlabeled`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FeatureGroupingAlgorithmUnlabeled {
    /// The pair finder's parameters, which are this algorithm's parameters
    /// (the source inserts `StablePairFinder().getParameters()` at the root).
    pub pair_finder: StablePairFinder,
    /// The incremental state (source `pairfinder_input_`): the consensus map
    /// computed so far, and the map added last.
    pairfinder_input: [ConsensusMap; 2],
    /// Whether the first element's intensity range is current. The source
    /// reads stored ranges: [`Self::set_reference`] converts through
    /// `MapConversion`, which updates them, but a pair-finder result never
    /// has them updated, so after the first [`Self::add_to_group`] the range
    /// of the growing map is empty and does not contribute to the
    /// normalising intensity.
    reference_range_current: bool,
}

impl FeatureGroupingAlgorithmUnlabeled {
    /// The source defaults tree, identical to
    /// [`StablePairFinder::defaults`].
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        StablePairFinder::defaults()
    }

    /// The current parameters.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        self.pair_finder.to_param()
    }

    /// Build from parameters, as the constructor followed by `setParameters`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong type or
    /// violates its restriction.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`Self::from_param`], also returning the unknown-parameter warnings.
    ///
    /// # Errors
    ///
    /// As [`Self::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(FEATURE_GROUPING_ALGORITHM_UNLABELED_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        let (pair_finder, warnings) =
            handler.set_parameters_with(parameters, StablePairFinder::from_complete_param)?;
        Ok((
            Self {
                pair_finder,
                ..Self::default()
            },
            warnings,
        ))
    }

    /// Set the reference map of the incremental path to the features of
    /// `map`, as map `map_id` (source `setReference` with a `FeatureMap`).
    /// `map_id` must be unique; see [`Self::add_to_group`].
    ///
    /// The source takes an `int` and converts it to the unsigned map index; a
    /// negative id therefore wraps there. Here the index is unsigned.
    ///
    /// # Errors
    ///
    /// The [`MapConversion::feature_map_to_consensus`] error.
    pub fn set_reference(&mut self, map_id: u64, map: &FeatureMap) -> Result<()> {
        self.pairfinder_input[0] = MapConversion::feature_map_to_consensus(map_id, map, None)?;
        self.reference_range_current = true;
        Ok(())
    }

    /// As [`Self::set_reference`] for the MS1 peaks of a peak map (source
    /// `setReference` with a `PeakMap`), whose elements get unique IDs from
    /// `generator`.
    ///
    /// # Errors
    ///
    /// The [`MapConversion::peak_map_to_consensus`] error.
    pub fn set_reference_peak_map(
        &mut self,
        map_id: u64,
        map: &MSExperiment,
        generator: &mut UniqueIdGenerator,
    ) -> Result<()> {
        self.pairfinder_input[0] =
            MapConversion::peak_map_to_consensus(map_id, map, None, generator)?;
        self.reference_range_current = true;
        Ok(())
    }

    /// The consensus map computed so far (source `getResultMap`). This is the
    /// raw incremental result, not post-processed.
    pub fn result_map(&self) -> &ConsensusMap {
        &self.pairfinder_input[0]
    }

    /// Mutable access to the consensus map computed so far (source
    /// `getResultMap` returns a non-const reference).
    pub fn result_map_mut(&mut self) -> &mut ConsensusMap {
        &mut self.pairfinder_input[0]
    }

    /// Add one map to the group (source `addToGroup`): convert it as map
    /// `map_id`, pair it with the consensus map computed so far, and make the
    /// result the new consensus map.
    ///
    /// Assumes [`Self::set_reference`] was called; without it the consensus
    /// map starts empty. Every `map_id` must differ from all earlier ones and
    /// from the reference's. Only a clash with the column headers of the
    /// current consensus map is detected, and, as in the source, a pair-finder
    /// result carries no column headers, so after the first call nothing is
    /// detected.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `map_id` clashes as described
    /// (source `Exception::IllegalArgument`), and the conversion and
    /// [`StablePairFinder`] errors. On error the state is unchanged.
    pub fn add_to_group(&mut self, map_id: u64, feature_map: &FeatureMap) -> Result<()> {
        let converted = MapConversion::feature_map_to_consensus(map_id, feature_map, None)?;
        let input = [self.pairfinder_input[0].clone(), converted];
        let max_intensity = self.max_intensity(&input)?;
        let mut result = ConsensusMap::new();
        self.pair_finder
            .run_with_max_intensity(&input, &mut result, max_intensity)?;
        let [_, converted] = input;
        self.pairfinder_input = [result, converted];
        self.reference_range_current = false;
        Ok(())
    }

    /// The source's normalising intensity for `input`, honouring a stale
    /// range of the first map.
    fn max_intensity(&self, input: &[ConsensusMap; 2]) -> Result<f64> {
        if self.reference_range_current {
            StablePairFinder::max_intensity(&[&input[0], &input[1]])
        } else {
            StablePairFinder::max_intensity(&[&input[1]])
        }
    }
}

impl FeatureGroupingAlgorithm for FeatureGroupingAlgorithmUnlabeled {
    /// Group all maps at once (source `group`).
    ///
    /// The reference is the map with the most features (the first of equal
    /// ones; map 0 when all are empty). Its features become singleton
    /// consensus features, and every other map, in input order, is paired
    /// with the growing consensus map by [`StablePairFinder`], each handle
    /// tagged with its map's position.
    ///
    /// The normalising intensity follows the source's stored ranges: the
    /// first pairing uses both maps, later ones only the map being added,
    /// because the growing map is a pair-finder result whose range was never
    /// updated. It matters only with a non-zero intensity weight.
    ///
    /// `out` keeps its column headers and metadata. Everything else is
    /// replaced by the fresh consensus map the source swaps in: an unset
    /// unique ID, empty document identifier, `label-free` experiment type and
    /// no data processing; then [`postprocess`] adds the identifications and
    /// sorts.
    ///
    /// Does not touch the incremental state.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for fewer than two maps (source
    /// `Exception::IllegalArgument`), and the conversion and
    /// [`StablePairFinder`] errors. On error `out` is unchanged.
    fn group(&mut self, maps: &[FeatureMap], out: &mut ConsensusMap) -> Result<()> {
        if maps.len() < 2 {
            return Err(Error::InvalidValue(
                "At least two maps must be given!".into(),
            ));
        }
        let mut reference = 0;
        let mut max_count = 0;
        for (m, map) in maps.iter().enumerate() {
            if map.len() > max_count {
                max_count = map.len();
                reference = m;
            }
        }

        let mut current =
            MapConversion::feature_map_to_consensus(reference as u64, &maps[reference], None)?;
        let mut current_range = true;
        for (i, map) in maps.iter().enumerate() {
            if i == reference {
                continue;
            }
            let added = MapConversion::feature_map_to_consensus(i as u64, map, None)?;
            let input = [current, added];
            let max_intensity = if current_range {
                StablePairFinder::max_intensity(&[&input[0], &input[1]])?
            } else {
                StablePairFinder::max_intensity(&[&input[1]])?
            };
            let mut result = ConsensusMap::new();
            self.pair_finder
                .run_with_max_intensity(&input, &mut result, max_intensity)?;
            current = result;
            current_range = false;
        }

        let mut grouped = ConsensusMap {
            features: current.features,
            column_headers: out.column_headers.clone(),
            metadata: out.metadata.clone(),
            ..ConsensusMap::default()
        };
        let inputs: Vec<_> = maps.iter().map(feature_map_identifications).collect();
        postprocess(&inputs, &mut grouped)?;
        *out = grouped;
        Ok(())
    }
}
