// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Feature grouping for labeling techniques with two labels.
//!
//! Native equivalent of `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h`.
//! See `docs/FEATURE_GROUPING_SUPPORT.md`.
//!
//! The algorithm takes one map and searches for corresponding features with a
//! defined distance in RT and m/z, using
//! [`LabeledPairFinder`](crate::analysis::mapmatching::labeled_pair_finder::LabeledPairFinder).

use super::feature_grouping_algorithm::FeatureGroupingAlgorithm;
use super::labeled_pair_finder::{LabeledPairFinder, RtPairEstimate};
use super::params::bad;
use crate::concept::UniqueIdGenerator;
use crate::error::Result;
use crate::kernel::conversion_helper::MapConversion;
use crate::kernel::features::{ConsensusMap, FeatureMap};
use crate::param::{DefaultParamHandler, Param};

/// Parameter-handler name of the source class.
pub const FEATURE_GROUPING_ALGORITHM_LABELED_NAME: &str = "FeatureGroupingAlgorithmLabeled";

/// Labeled feature grouping (source `FeatureGroupingAlgorithmLabeled`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FeatureGroupingAlgorithmLabeled {
    /// The pair finder's parameters, which are this algorithm's parameters
    /// (the source inserts `LabeledPairFinder().getParameters()` at the root).
    /// A fresh finder with these parameters runs for every call, as in the
    /// source.
    pub pair_finder: LabeledPairFinder,
}

impl FeatureGroupingAlgorithmLabeled {
    /// The source defaults tree, identical to [`LabeledPairFinder::defaults`].
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        LabeledPairFinder::defaults()
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
    /// Returns [`Error::InvalidValue`](crate::error::Error::InvalidValue) when
    /// a parameter has the wrong type or violates its restriction.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`Self::from_param`], also returning the unknown-parameter warnings.
    ///
    /// # Errors
    ///
    /// As [`Self::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(FEATURE_GROUPING_ALGORITHM_LABELED_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        let (pair_finder, warnings) =
            handler.set_parameters_with(parameters, LabeledPairFinder::from_complete_param)?;
        Ok((Self { pair_finder }, warnings))
    }

    /// Group the features of the single map in `maps` into light/heavy pairs,
    /// drawing the pairs' unique IDs from `generator`, and report the RT
    /// parameters the pair finder used.
    ///
    /// The map is converted to a consensus map with map index 0 (source
    /// `MapConversion::convert(0, maps[0], input[0])`) and handed to
    /// [`LabeledPairFinder::run_with_generator`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`](crate::error::Error::InvalidValue)
    /// (source `Exception::IllegalArgument`) unless exactly one map is given
    /// and `out` has exactly two column headers; the conversion errors; and
    /// those of [`LabeledPairFinder::run_with_generator`]. On error `out` is
    /// unchanged.
    pub fn group_with_generator(
        &self,
        maps: &[FeatureMap],
        out: &mut ConsensusMap,
        generator: &mut UniqueIdGenerator,
    ) -> Result<RtPairEstimate> {
        if maps.len() != 1 {
            return Err(bad("Exactly one map must be given!"));
        }
        if out.column_headers.len() != 2 {
            return Err(bad("Two file descriptions must be set in 'out'!"));
        }
        let input = [MapConversion::feature_map_to_consensus(0, &maps[0], None)?];
        self.pair_finder.run_with_generator(&input, out, generator)
    }
}

impl FeatureGroupingAlgorithm for FeatureGroupingAlgorithmLabeled {
    /// [`FeatureGroupingAlgorithmLabeled::group_with_generator`] with a
    /// clock-seeded generator, discarding the RT report.
    fn group(&mut self, maps: &[FeatureMap], out: &mut ConsensusMap) -> Result<()> {
        let mut generator = UniqueIdGenerator::new();
        self.group_with_generator(maps, out, &mut generator)
            .map(|_| ())
    }
}
