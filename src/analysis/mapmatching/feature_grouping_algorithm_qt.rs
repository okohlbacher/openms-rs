// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Feature grouping for unlabeled data by QT clustering.
//!
//! Native equivalent of `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h`.
//! See `docs/FEATURE_GROUPING_SUPPORT.md`.
//!
//! The algorithm takes feature or consensus maps and searches for
//! corresponding (consensus) features across them. The maps must be aligned
//! (retention-time distortions corrected, for example by a map-alignment
//! algorithm), but small deviations are tolerated. It accumulates the features
//! of all maps and applies the QT clustering variant of
//! [`QTClusterFinder`](crate::analysis::mapmatching::qt_cluster_finder::QTClusterFinder).

use super::feature_grouping_algorithm::{
    FeatureGroupingAlgorithm, consensus_map_identifications, feature_map_identifications,
    postprocess,
};
use super::qt_cluster_finder::{QTClusterFinder, QtInputMap};
use crate::error::{Error, Result};
use crate::kernel::features::{ConsensusMap, FeatureMap};
use crate::param::{DefaultParamHandler, Param};

/// Parameter-handler name of the source class.
pub const FEATURE_GROUPING_ALGORITHM_QT_NAME: &str = "FeatureGroupingAlgorithmQT";

/// QT feature grouping (source `FeatureGroupingAlgorithmQT`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FeatureGroupingAlgorithmQT {
    /// The cluster finder's parameters, which are this algorithm's parameters
    /// (the source inserts `QTClusterFinder().getParameters()` at the root). A
    /// fresh finder with these parameters runs for every call, as in the
    /// source.
    pub cluster_finder: QTClusterFinder,
}

impl FeatureGroupingAlgorithmQT {
    /// The source defaults tree, identical to [`QTClusterFinder::defaults`].
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        QTClusterFinder::defaults()
    }

    /// The current parameters.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        self.cluster_finder.to_param()
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
        let mut handler = DefaultParamHandler::new(FEATURE_GROUPING_ALGORITHM_QT_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        let (cluster_finder, warnings) =
            handler.set_parameters_with(parameters, QTClusterFinder::from_complete_param)?;
        Ok((Self { cluster_finder }, warnings))
    }

    /// Source template `group_`: check the map count, run a fresh cluster
    /// finder, and post-process.
    fn group_impl<M: QtInputMap>(
        &self,
        maps: &[M],
        identifications: Vec<super::feature_grouping_algorithm::InputIdentifications<'_>>,
        out: &mut ConsensusMap,
    ) -> Result<()> {
        if maps.len() < 2 {
            return Err(Error::InvalidValue(
                "At least two maps must be given!".into(),
            ));
        }
        let mut finder = self.cluster_finder.fresh();
        let mut result = out.clone();
        finder.run_maps(maps, &mut result)?;
        postprocess(&identifications, &mut result)?;
        *out = result;
        Ok(())
    }
}

impl FeatureGroupingAlgorithm for FeatureGroupingAlgorithmQT {
    /// Group feature maps (source `group(const std::vector<FeatureMap>&, ...)`).
    ///
    /// `out`'s features are replaced by the clusters; its metadata, column
    /// headers and identifications stay, and [`postprocess`] appends the input
    /// maps' identifications and sorts. The source requires the input ranges
    /// to be current; here they are computed.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for fewer than two maps (source
    /// `Exception::IllegalArgument`), and the
    /// [`QTClusterFinder::run_maps`] errors. On error `out` is unchanged.
    fn group(&mut self, maps: &[FeatureMap], out: &mut ConsensusMap) -> Result<()> {
        let identifications = maps.iter().map(feature_map_identifications).collect();
        self.group_impl(maps, identifications, out)
    }

    /// Group consensus maps directly (source
    /// `group(const std::vector<ConsensusMap>&, ...)`), as for feature maps.
    ///
    /// # Errors
    ///
    /// As [`FeatureGroupingAlgorithm::group`].
    fn group_consensus(&mut self, maps: &[ConsensusMap], out: &mut ConsensusMap) -> Result<()> {
        let identifications = maps.iter().map(consensus_map_identifications).collect();
        self.group_impl(maps, identifications, out)
    }
}
