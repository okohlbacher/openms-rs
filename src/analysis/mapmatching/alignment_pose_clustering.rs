// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Map alignment based on pose clustering.
//!
//! Ports `ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h` and its
//! `.cpp`. The source composes a `StablePairFinder`; so does this port
//! ([`StablePairFinder`](crate::analysis::mapmatching::stable_pair_finder::StablePairFinder),
//! configured by the `pairfinder:` section). The `*_with_pair_finder`
//! variants of the `align_*` methods take the pairing step as a closure
//! instead, for tests that pin the surrounding algorithm independently of the
//! pair finder. See `docs/MAP_ALIGNMENT_SUPPORT.md`.

use super::base_group_finder::BaseGroupFinder;
use super::pose_clustering_affine::PoseClusteringAffineSuperimposer;
use super::stable_pair_finder::StablePairFinder;
use super::superimposer::{BaseSuperimposer, entry};
use crate::analysis::transformations::{DataPoint, TransformationDescription};
use crate::concept::progress_logger::ProgressReporter;
use crate::concept::unique_id::UniqueIdGenerator;
use crate::kernel::conversion_helper::MapConversion;
use crate::kernel::{ConsensusMap, FeatureMap, MSExperiment};
use crate::param::{DefaultParamHandler, Param, ParamNode, ParamValue};
use crate::{Error, Result};

/// Parameter handler name of the source class.
pub const MAP_ALIGNMENT_ALGORITHM_POSE_CLUSTERING_NAME: &str =
    "MapAlignmentAlgorithmPoseClustering";

/// A map alignment algorithm based on pose clustering (source
/// `MapAlignmentAlgorithmPoseClustering`).
///
/// Pose clustering analyses pair distances to find the most probable
/// transformation of retention times (Lange et al., "A Geometric Approach
/// for the Alignment of Liquid Chromatography-Mass Spectrometry Data", ISMB
/// 2007). The algorithm uses the first `max_num_peaks_considered` features
/// of a feature map, or the most intense peaks of a peak map, per map; larger
/// values prolong the computation, smaller ones may give no or unstable
/// transformations, and -1 uses all (possibly very slowly).
///
/// Aligning a map runs the affine superimposer against the reference,
/// applies its transformation to the map, pairs the features of the
/// reference and the transformed map with the pair finder, and fits a
/// `linear` model to the pairs' retention times, the map's taken back
/// through the inverse of the superimposer's transformation.
#[derive(Clone, Debug, PartialEq)]
pub struct MapAlignmentAlgorithmPoseClustering {
    /// The superimposer (source `superimposer:` section).
    pub superimposer: PoseClusteringAffineSuperimposer,
    /// The pair finder (source `pairfinder:` section).
    pub pair_finder: StablePairFinder,
    /// Maximal number of peaks or features considered per map; -1 uses all
    /// (source `max_num_peaks_considered`, default 1000, minimum -1).
    pub max_num_peaks_considered: i64,
    reference: ConsensusMap,
}

impl Default for MapAlignmentAlgorithmPoseClustering {
    fn default() -> Self {
        Self {
            superimposer: PoseClusteringAffineSuperimposer::default(),
            pair_finder: StablePairFinder::default(),
            max_num_peaks_considered: 1000,
            reference: ConsensusMap::new(),
        }
    }
}

impl MapAlignmentAlgorithmPoseClustering {
    /// The source parameter defaults: the superimposer's under
    /// `superimposer:`, the pair finder's under `pairfinder:` and
    /// `max_num_peaks_considered`.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters in the layout of
    /// [`MapAlignmentAlgorithmPoseClustering::defaults`] (source
    /// `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let mut peaks = entry(
            "max_num_peaks_considered",
            ParamValue::Integer(self.max_num_peaks_considered),
            "The maximal number of peaks/features to be considered per map. To use all, set to '-1'.",
            false,
        );
        peaks.min_int = -1;
        let mut param = Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries: vec![peaks],
            nodes: Vec::new(),
        })?;
        param.insert("superimposer:", &self.superimposer.to_param()?)?;
        param.insert("pairfinder:", &self.pair_finder.to_param()?)?;
        Ok(param)
    }

    /// Apply parameters (source `setParameters` and `updateMembers_`): the
    /// `superimposer:` section configures the superimposer and the
    /// `pairfinder:` section the pair finder. The reference is kept.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong type or
    /// violates its restriction; nothing changes then.
    pub fn set_parameters(&mut self, parameters: &Param) -> Result<Vec<String>> {
        let mut handler = DefaultParamHandler::new(MAP_ALIGNMENT_ALGORITHM_POSE_CLUSTERING_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        let reference = self.reference.clone();
        let (updated, mut warnings) = handler.set_parameters_with(parameters, |p| {
            let (superimposer, mut nested) =
                PoseClusteringAffineSuperimposer::from_param_with_warnings(
                    &p.copy("superimposer:", true)?,
                )?;
            let (pair_finder, pair_finder_nested) =
                StablePairFinder::from_param_with_warnings(&p.copy("pairfinder:", true)?)?;
            nested.extend(pair_finder_nested);
            Ok((
                Self {
                    superimposer,
                    pair_finder,
                    max_num_peaks_considered: p.value("max_num_peaks_considered")?.to_i64()?,
                    reference,
                },
                nested,
            ))
        })?;
        let (updated, nested) = updated;
        warnings.extend(nested);
        *self = updated;
        Ok(warnings)
    }

    /// A new algorithm with parameters (source construction followed by
    /// `setParameters`).
    ///
    /// # Errors
    ///
    /// As [`MapAlignmentAlgorithmPoseClustering::set_parameters`].
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`MapAlignmentAlgorithmPoseClustering::from_param`], also returning
    /// the unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`MapAlignmentAlgorithmPoseClustering::set_parameters`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut algorithm = Self::default();
        let warnings = algorithm.set_parameters(parameters)?;
        Ok((algorithm, warnings))
    }

    fn limit(&self) -> Option<usize> {
        usize::try_from(self.max_num_peaks_considered).ok()
    }

    /// The reference map, as converted.
    pub fn reference(&self) -> &ConsensusMap {
        &self.reference
    }

    /// Set a feature map as the reference (source `setReference`): its first
    /// `max_num_peaks_considered` features become the reference, as map 0.
    ///
    /// # Errors
    ///
    /// The errors of `MapConversion::feature_map_to_consensus`.
    pub fn set_reference_features(&mut self, map: &FeatureMap) -> Result<()> {
        self.reference = MapConversion::feature_map_to_consensus(0, map, self.limit())?;
        Ok(())
    }

    /// Set a peak map as the reference (source `setReference`): its most
    /// intense MS1 peaks become the reference, as map 0.
    ///
    /// # Errors
    ///
    /// The errors of `MapConversion::peak_map_to_consensus`.
    pub fn set_reference_peaks(
        &mut self,
        map: &MSExperiment,
        generator: &mut UniqueIdGenerator,
    ) -> Result<()> {
        self.reference = MapConversion::peak_map_to_consensus(0, map, self.limit(), generator)?;
        Ok(())
    }

    /// Set a consensus map as the reference (source `setReference`, whose
    /// consensus-map conversion keeps the first `max_num_peaks_considered`
    /// consensus features).
    pub fn set_reference_consensus(&mut self, map: &ConsensusMap) {
        let mut reference = map.clone();
        if let Some(limit) = self.limit() {
            reference.features.truncate(limit);
        }
        self.reference = reference;
    }

    /// Align a feature map to the reference (source `align(const
    /// FeatureMap&, ...)`): its first `max_num_peaks_considered` features, as
    /// map 1, are aligned by [`align_consensus`](Self::align_consensus).
    ///
    /// # Errors
    ///
    /// As [`align_consensus`](Self::align_consensus), and the conversion's.
    pub fn align_features(
        &self,
        map: &FeatureMap,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<TransformationDescription> {
        let scene = MapConversion::feature_map_to_consensus(1, map, self.limit())?;
        self.align_consensus(&scene, progress)
    }

    /// Align a peak map to the reference (source `align(const PeakMap&,
    /// ...)`): its most intense MS1 peaks, as map 1.
    ///
    /// # Errors
    ///
    /// As [`align_consensus`](Self::align_consensus), and the conversion's.
    pub fn align_peaks(
        &self,
        map: &MSExperiment,
        generator: &mut UniqueIdGenerator,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<TransformationDescription> {
        let scene = MapConversion::peak_map_to_consensus(1, map, self.limit(), generator)?;
        self.align_consensus(&scene, progress)
    }

    /// Align a consensus map to the reference (source `align(const
    /// ConsensusMap&, ...)`) with the configured [`StablePairFinder`].
    ///
    /// The superimposer's transformation is applied to every consensus
    /// feature and to its **first** handle only, as in the source (the maps
    /// the source converts have one handle per feature). The pair finder
    /// then groups the reference (index 0) and the transformed map (index 1).
    /// For every group of exactly two handles, the reference handle's RT is
    /// paired with the other's, taken back through the algebraic inverse of
    /// the superimposer's transformation, and a `linear` model is fitted to
    /// those pairs.
    ///
    /// # Errors
    ///
    /// Returns the superimposer's and the pair finder's errors, and
    /// [`Error::InvalidValue`] when a consensus feature of `map` has no
    /// handle (the source then writes through an end iterator) or no pairs
    /// are found (the source's `Exception::IllegalArgument` from the linear
    /// fit).
    pub fn align_consensus(
        &self,
        map: &ConsensusMap,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<TransformationDescription> {
        self.align_consensus_with_pair_finder(
            map,
            |input| {
                let mut finder = self.pair_finder.clone();
                let mut result = ConsensusMap::new();
                finder.run(input, &mut result)?;
                Ok(result)
            },
            progress,
        )
    }

    /// [`align_features`](Self::align_features) with the pairing step given
    /// as a closure over the reference (index 0) and the transformed scene
    /// (index 1) instead of the configured pair finder.
    ///
    /// # Errors
    ///
    /// As [`align_consensus_with_pair_finder`](Self::align_consensus_with_pair_finder),
    /// and the conversion's.
    pub fn align_features_with_pair_finder(
        &self,
        map: &FeatureMap,
        pair_finder: impl FnOnce(&[ConsensusMap; 2]) -> Result<ConsensusMap>,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<TransformationDescription> {
        let scene = MapConversion::feature_map_to_consensus(1, map, self.limit())?;
        self.align_consensus_with_pair_finder(&scene, pair_finder, progress)
    }

    /// [`align_peaks`](Self::align_peaks) with the pairing step given as a
    /// closure instead of the configured pair finder.
    ///
    /// # Errors
    ///
    /// As [`align_consensus_with_pair_finder`](Self::align_consensus_with_pair_finder),
    /// and the conversion's.
    pub fn align_peaks_with_pair_finder(
        &self,
        map: &MSExperiment,
        generator: &mut UniqueIdGenerator,
        pair_finder: impl FnOnce(&[ConsensusMap; 2]) -> Result<ConsensusMap>,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<TransformationDescription> {
        let scene = MapConversion::peak_map_to_consensus(1, map, self.limit(), generator)?;
        self.align_consensus_with_pair_finder(&scene, pair_finder, progress)
    }

    /// [`align_consensus`](Self::align_consensus) with the pairing step given
    /// as a closure instead of the configured pair finder.
    ///
    /// The superimposer's transformation is applied to every consensus
    /// feature and to its **first** handle only, as in the source (the maps
    /// the source converts have one handle per feature). `pair_finder` then
    /// groups the reference (index 0) and the transformed map (index 1), the
    /// source's `StablePairFinder::run`. For every group of exactly two
    /// handles, the reference handle's RT is paired with the other's, taken
    /// back through the algebraic inverse of the superimposer's
    /// transformation, and a `linear` model is fitted to those pairs.
    ///
    /// # Errors
    ///
    /// Returns the superimposer's and the pair finder's errors, and
    /// [`Error::InvalidValue`] when a consensus feature of `map` has no
    /// handle (the source then writes through an end iterator) or no pairs
    /// are found (the source's `Exception::IllegalArgument` from the linear
    /// fit).
    pub fn align_consensus_with_pair_finder(
        &self,
        map: &ConsensusMap,
        pair_finder: impl FnOnce(&[ConsensusMap; 2]) -> Result<ConsensusMap>,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<TransformationDescription> {
        let transformation = self
            .superimposer
            .superimpose(&self.reference, map, progress)?
            .transformation;
        let mut scene = map.clone();
        for feature in &mut scene.features {
            let rt = transformation.apply(feature.rt)?;
            feature.rt = rt;
            let mut handles = feature.handles().to_vec();
            let first = handles.first_mut().ok_or_else(|| {
                Error::InvalidValue("a consensus feature to align has no handle".into())
            })?;
            first.rt = rt;
            feature.set_handles(handles)?;
        }
        let result = pair_finder(&[self.reference.clone(), scene])?;
        let inverse = transformation.inverse()?;
        let mut data = Vec::new();
        for feature in &result.features {
            if let [first, second] = feature.handles() {
                let y = first.rt;
                let x = inverse.apply(second.rt)?;
                data.push(if second.map_index != 0 {
                    DataPoint::new(x, y)
                } else {
                    DataPoint::new(y, x)
                });
            }
        }
        let mut description = TransformationDescription::new(data)?;
        description.fit_model_with_param("linear", &Param::new())?;
        Ok(description)
    }
}
