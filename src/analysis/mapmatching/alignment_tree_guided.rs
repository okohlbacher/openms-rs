// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tree-guided map alignment based on peptide identifications.
//!
//! Ports `ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmTreeGuided.h` and its
//! `.cpp`, except for the clustering in `buildTree`: the source clusters the
//! maps with `ClusterHierarchical` and `AverageLinkage` (`ML/CLUSTERING`) and
//! prints the tree with `ClusterAnalyzer::newickTree`, whose ports are not in
//! this crate yet. The similarity the clustering uses,
//! [`peptide_identifications_pearson_similarity`](crate::analysis::mapmatching::alignment_tree_guided::peptide_identifications_pearson_similarity), and the per-map data it
//! runs on, [`MapAlignmentAlgorithmTreeGuided::extract_seq_and_rt`](crate::analysis::mapmatching::alignment_tree_guided::MapAlignmentAlgorithmTreeGuided::extract_seq_and_rt), are
//! ported; the tree is a parameter of
//! [`MapAlignmentAlgorithmTreeGuided::align_with_tree`](crate::analysis::mapmatching::alignment_tree_guided::MapAlignmentAlgorithmTreeGuided::align_with_tree). See
//! `docs/MAP_ALIGNMENT_SUPPORT.md`.

use super::alignment_identification::MapAlignmentAlgorithmIdentification;
use super::superimposer::{entry, tree};
use crate::analysis::alignment_transformer::MapAlignmentTransformer;
use crate::analysis::transformations::{TransformationDescription, transformation_model_defaults};
use crate::concept::unique_id::UniqueIdGenerator;
use crate::format::FileType;
use crate::kernel::{Feature, FeatureMap};
use crate::math::source_sort::source_sort_by;
use crate::math::statistic_functions::pearson_correlation_coefficient;
use crate::param::{DefaultParamHandler, Param, ParamValue};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

/// Parameter handler name of the source class.
pub const MAP_ALIGNMENT_ALGORITHM_TREE_GUIDED_NAME: &str = "MapAlignmentAlgorithmTreeGuided";

/// The model types `model_type` accepts, in the source's order.
pub const TREE_GUIDED_MODEL_TYPES: [&str; 4] = ["linear", "b_spline", "lowess", "interpolated"];

/// Feature retention times per peptide sequence of one map (source
/// `SeqAndRTList`), ordered as `std::map<std::string, DoubleList>`.
pub type SeqAndRtList = BTreeMap<String, Vec<f64>>;

/// One merge step of the guide tree: the two clusters joined, each named by
/// the index of a map in it (the `left_child` and `right_child` of the
/// source's `BinaryTreeNode`; its `distance` is not used by the alignment).
pub type TreeMerge = (usize, usize);

/// The result of [`MapAlignmentAlgorithmTreeGuided::tree_guided_alignment`].
#[derive(Clone, Debug)]
pub struct TreeGuidedMerge {
    /// All features of the aligned maps with transformed retention times, in
    /// the order the clusters were joined (source `map_transformed`).
    pub map_transformed: FeatureMap,
    /// The map indices in the order of `map_transformed` (source
    /// `trafo_order`).
    pub trafo_order: Vec<usize>,
    /// The lines the source writes to `OPENMS_LOG_INFO` on the way, without
    /// line ends.
    pub info: Vec<String>,
    /// The lines the source writes to `OPENMS_LOG_WARN` on the way.
    pub warnings: Vec<String>,
}

/// The result of [`MapAlignmentAlgorithmTreeGuided::compute_trafos_by_original_rt`]
/// and [`MapAlignmentAlgorithmTreeGuided::align_with_tree`].
#[derive(Clone, Debug)]
pub struct TreeGuidedAlignment {
    /// One fitted transformation per input map, in input order.
    pub transformations: Vec<TransformationDescription>,
    /// The lines the source writes to `OPENMS_LOG_INFO`, without line ends.
    pub info: Vec<String>,
    /// The lines the source writes to `OPENMS_LOG_WARN`.
    pub warnings: Vec<String>,
}

/// A map alignment algorithm based on peptide identifications from MS2
/// spectra, guided by a hierarchical clustering of the maps (source
/// `MapAlignmentAlgorithmTreeGuided`).
///
/// Identifications with the same sequence in different maps are points of
/// correspondence in RT between the maps; only the first hit of each
/// identification is used. They give the similarity of two maps
/// ([`peptide_identifications_pearson_similarity`]), from which average-linkage
/// hierarchical clustering builds a binary tree. Following the tree, as long
/// as there are at least two clusters, the cluster with the larger 10/90
/// percentile RT range of each pair is the reference for
/// [`MapAlignmentAlgorithmIdentification::align`], which aligns the median
/// RT of each identification group of the other cluster to the reference RT
/// of the group; the model of `model_type` (by default a smoothing cubic
/// B-spline) turns this into a smooth function, which transforms the other
/// cluster, keeping its original retention times in the meta value
/// `original_RT`; and the two clusters are combined.
///
/// The combined map is then the reference against which each original input
/// map is aligned again, which gives the final transformations; each map's
/// retention times are transformed with its fitted model.
#[derive(Clone, Debug, PartialEq)]
pub struct MapAlignmentAlgorithmTreeGuided {
    /// Type of the RT transformation model (source `model_type`, one of
    /// [`TREE_GUIDED_MODEL_TYPES`], default `b_spline`).
    pub model_type: String,
    /// The `model:` section: `TransformationModelDefaults::getDefaults("b_spline")`
    /// by default; the subsection named by `model_type` configures the model.
    pub model: Param,
    /// The identification alignment (source `align_algorithm:` section and
    /// member `align_algorithm_`); its `use_feature_rt` defaults to `true`
    /// here.
    pub align_algorithm: MapAlignmentAlgorithmIdentification,
}

impl MapAlignmentAlgorithmTreeGuided {
    /// The algorithm with its default parameters (source constructor).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn new() -> Result<Self> {
        let mut align_algorithm = MapAlignmentAlgorithmIdentification::default();
        align_algorithm.use_feature_rt = true;
        Ok(Self {
            model_type: "b_spline".into(),
            model: transformation_model_defaults("b_spline")?,
            align_algorithm,
        })
    }

    /// The source parameter defaults: `model:` (the transformation model
    /// defaults for `b_spline`), `model_type` and `align_algorithm:` (the
    /// identification alignment's defaults with `use_feature_rt` set to
    /// `true`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::new()?.to_param()
    }

    /// The parameters in the layout of
    /// [`MapAlignmentAlgorithmTreeGuided::defaults`] (source
    /// `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let mut model_type = entry(
            "model_type",
            ParamValue::String(self.model_type.clone()),
            "Options to control the modeling of retention time transformations from data",
            false,
        );
        model_type.valid_strings = TREE_GUIDED_MODEL_TYPES.iter().map(|&t| t.into()).collect();
        let mut param = tree(vec![model_type])?;
        param.insert("model:", &self.model)?;
        param.insert("align_algorithm:", &self.align_algorithm.to_param()?)?;
        Ok(param)
    }

    /// Apply parameters (source `setParameters` and `updateMembers_`). The
    /// identification alignment keeps its state, as the source's member does.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong type or
    /// violates its restriction; nothing changes then.
    pub fn set_parameters(&mut self, parameters: &Param) -> Result<Vec<String>> {
        let mut handler = DefaultParamHandler::new(MAP_ALIGNMENT_ALGORITHM_TREE_GUIDED_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        let align = self.align_algorithm.clone();
        let (updated, mut warnings) = handler.set_parameters_with(parameters, |p| {
            let mut align_algorithm = align.clone();
            let nested = align_algorithm.set_parameters(&p.copy("align_algorithm:", true)?)?;
            Ok((
                Self {
                    model_type: p.value("model_type")?.as_str()?.to_owned(),
                    model: p.copy("model:", true)?,
                    align_algorithm,
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
    /// As [`MapAlignmentAlgorithmTreeGuided::set_parameters`].
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`MapAlignmentAlgorithmTreeGuided::from_param`], also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`MapAlignmentAlgorithmTreeGuided::set_parameters`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut algorithm = Self::new()?;
        let warnings = algorithm.set_parameters(parameters)?;
        Ok((algorithm, warnings))
    }

    /// The parameters of the selected model: the subsection `model_type` of
    /// `model:` (source member `model_param_`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn model_param(&self) -> Result<Param> {
        self.model.copy(&format!("{}:", self.model_type), true)
    }

    /// For each map, the retention times of its features per sequence of
    /// their identifications, and the sorted retention times of all of them
    /// (source `extractSeqAndRt_` with `addPeptideSequences_`).
    ///
    /// Every identification of a feature with at least one hit adds the
    /// feature's retention time once under the sequence of its **first**
    /// hit, and once to the map's list; a feature with two such
    /// identifications counts twice. Unassigned identifications are not used.
    /// The second list per map (source `maps_ranges`) is sorted with the
    /// Release build's `std::sort` permutation.
    ///
    /// # Errors
    ///
    /// Only an error of the sort's allocation.
    pub fn extract_seq_and_rt(
        feature_maps: &[FeatureMap],
    ) -> Result<(Vec<SeqAndRtList>, Vec<Vec<f64>>)> {
        let mut lists = Vec::with_capacity(feature_maps.len());
        let mut ranges = Vec::with_capacity(feature_maps.len());
        for map in feature_maps {
            let mut list = SeqAndRtList::new();
            let mut range = Vec::new();
            for feature in &map.features {
                for peptide in &feature.peptide_identifications {
                    if let Some(hit) = peptide.hits.first() {
                        list.entry(hit.sequence.to_string())
                            .or_default()
                            .push(feature.rt);
                        range.push(feature.rt);
                    }
                }
            }
            source_sort_by(&mut range, |a, b| a < b)?;
            lists.push(list);
            ranges.push(range);
        }
        Ok((lists, ranges))
    }

    /// Align the clusters of `feature_maps_transformed` along `tree` (source
    /// `treeGuidedAlignment`).
    ///
    /// For each merge `(left, right)` in order: of the two maps, the one
    /// whose `maps_ranges` list has the larger 10/90 percentile range (the
    /// values at positions `floor(0.9 n)` and `floor(0.1 n)`; `right` on a
    /// tie) is the reference; the other is aligned to it with the
    /// identification alignment (reference index 1), the model is fitted, its
    /// retention times are transformed keeping `original_RT`, and it is
    /// appended to the reference as `FeatureMap::operator+=` does. The
    /// combined map is stored at the smaller of the two indices and the other
    /// cleared, as the source does. `maps_ranges` is not updated by a merge:
    /// later steps use each index's original range, as in the source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] for a map without identified
    /// features, as the source does ("FeatureMap originating from '...'
    /// contains no Peptide Identifications. Cannot align!").
    /// Returns [`Error::InvalidValue`] where the source's behaviour is
    /// undefined: no maps, `maps_ranges` not of the maps' length, a merge
    /// naming a map out of range or the same map twice; and for a feature of
    /// a merged map that refers to an identification graph, which the native
    /// feature map does not hold (the source merges the graphs). Also the
    /// errors of the identification alignment, the model fit and the
    /// transformation. The maps may be partly merged then.
    pub fn tree_guided_alignment(
        &mut self,
        tree: &[TreeMerge],
        feature_maps_transformed: &mut [FeatureMap],
        maps_ranges: &[Vec<f64>],
        generator: &mut UniqueIdGenerator,
    ) -> Result<TreeGuidedMerge> {
        let n = feature_maps_transformed.len();
        if n == 0 {
            return Err(Error::InvalidValue(
                "tree-guided alignment needs at least one map".into(),
            ));
        }
        if maps_ranges.len() != n {
            return Err(Error::InvalidValue(format!(
                "{} retention-time ranges for {n} maps",
                maps_ranges.len()
            )));
        }
        for (map, range) in feature_maps_transformed.iter().zip(maps_ranges) {
            if range.is_empty() {
                return Err(Error::MissingInformation(format!(
                    "FeatureMap originating from '{}' contains no Peptide Identifications. Cannot align!",
                    map.primary_ms_run_path()?.join("', '")
                )));
            }
        }
        let model_param = self.model_param()?;
        let transformer = MapAlignmentTransformer {
            store_original_rt: true,
            ..MapAlignmentTransformer::default()
        };
        let mut map_sets: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
        let mut last = 0;
        let mut info = Vec::new();
        let mut warnings = Vec::new();
        for &(left, right) in tree {
            if left >= n || right >= n || left == right {
                return Err(Error::InvalidValue(format!(
                    "tree node ({left}, {right}) does not join two of the {n} maps"
                )));
            }
            let (reference, to_transform) =
                if percentile_range(&maps_ranges[left]) > percentile_range(&maps_ranges[right]) {
                    (left, right)
                } else {
                    (right, left)
                };
            let to_align = [
                feature_maps_transformed[to_transform].clone(),
                feature_maps_transformed[reference].clone(),
            ];
            let aligned = self.align_algorithm.align(&to_align, Some(1))?;
            drop(to_align);
            info.extend(aligned.info);
            warnings.extend(aligned.warnings);
            let mut transformation =
                aligned.transformations.into_iter().next().ok_or_else(|| {
                    Error::InvalidValue(
                        "the identification alignment returned no transformation".into(),
                    )
                })?;
            transformation.fit_model_with_param(&self.model_type, &model_param)?;
            transformer.transform_feature_map(
                &mut feature_maps_transformed[to_transform],
                &transformation,
            )?;
            let other = feature_maps_transformed[to_transform].clone();
            info.extend(append_feature_map(
                &mut feature_maps_transformed[reference],
                &other,
                generator,
            )?);
            if reference < to_transform {
                feature_maps_transformed[to_transform].clear(true);
                last = reference;
            } else {
                feature_maps_transformed[to_transform] =
                    std::mem::take(&mut feature_maps_transformed[reference]);
                last = to_transform;
            }
            let joined = map_sets[to_transform].clone();
            map_sets[reference].extend(joined);
            map_sets[to_transform] = map_sets[reference].clone();
        }
        Ok(TreeGuidedMerge {
            map_transformed: feature_maps_transformed[last].clone(),
            trafo_order: map_sets[last].clone(),
            info,
            warnings,
        })
    }

    /// Align `feature_maps` along the given guide tree and transform their
    /// retention times (source `align`, with the tree `buildTree` would
    /// compute given by the caller).
    ///
    /// The maps are merged along the tree on a copy
    /// ([`MapAlignmentAlgorithmTreeGuided::tree_guided_alignment`]); each
    /// original map is then aligned to the merged map
    /// ([`MapAlignmentAlgorithmTreeGuided::compute_trafos_by_original_rt`])
    /// and transformed with its transformation
    /// ([`MapAlignmentAlgorithmTreeGuided::compute_transformed_feature_maps`]).
    /// The source's info line with the tree in Newick format is not produced.
    ///
    /// # Errors
    ///
    /// The errors of the three steps; the maps are left unchanged unless the
    /// last one fails.
    pub fn align_with_tree(
        &mut self,
        feature_maps: &mut [FeatureMap],
        tree: &[TreeMerge],
        generator: &mut UniqueIdGenerator,
    ) -> Result<TreeGuidedAlignment> {
        let (_, maps_ranges) = Self::extract_seq_and_rt(feature_maps)?;
        let merged = {
            let mut copied = feature_maps.to_vec();
            self.tree_guided_alignment(tree, &mut copied, &maps_ranges, generator)?
        };
        let mut result =
            self.compute_trafos_by_original_rt(feature_maps, &merged.map_transformed)?;
        let mut info = merged.info;
        info.append(&mut result.info);
        let mut warnings = merged.warnings;
        warnings.append(&mut result.warnings);
        Self::compute_transformed_feature_maps(feature_maps, &result.transformations)?;
        result.info = info;
        result.warnings = warnings;
        Ok(result)
    }

    /// The transformation of each original map onto the retention-time scale
    /// of `map_transformed` (source `computeTrafosByOriginalRT`): each map is
    /// aligned to `map_transformed` with the identification alignment
    /// (reference index 1) and the model fitted. The data points are the
    /// median feature RT per sequence of the original map against the
    /// merged map's. The source's unused `trafo_order` argument is dropped.
    ///
    /// # Errors
    ///
    /// The errors of the identification alignment and of the model fit.
    pub fn compute_trafos_by_original_rt(
        &mut self,
        feature_maps: &[FeatureMap],
        map_transformed: &FeatureMap,
    ) -> Result<TreeGuidedAlignment> {
        let model_param = self.model_param()?;
        let mut to_align = [FeatureMap::new(), map_transformed.clone()];
        let mut transformations = Vec::with_capacity(feature_maps.len());
        let mut info = Vec::new();
        let mut warnings = Vec::new();
        for map in feature_maps {
            to_align[0] = map.clone();
            let aligned = self.align_algorithm.align(&to_align, Some(1))?;
            info.extend(aligned.info);
            warnings.extend(aligned.warnings);
            let mut transformation =
                aligned.transformations.into_iter().next().ok_or_else(|| {
                    Error::InvalidValue(
                        "the identification alignment returned no transformation".into(),
                    )
                })?;
            transformation.fit_model_with_param(&self.model_type, &model_param)?;
            transformations.push(transformation);
        }
        Ok(TreeGuidedAlignment {
            transformations,
            info,
            warnings,
        })
    }

    /// Transform the retention times of each map with its transformation,
    /// keeping the original ones in `original_RT` (source
    /// `computeTransformedFeatureMaps`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when there are fewer transformations
    /// than maps (the source reads past the end), and the errors of
    /// [`MapAlignmentTransformer::transform_feature_map`]; the maps before
    /// the failing one are transformed then.
    pub fn compute_transformed_feature_maps(
        feature_maps: &mut [FeatureMap],
        transformations: &[TransformationDescription],
    ) -> Result<()> {
        if transformations.len() < feature_maps.len() {
            return Err(Error::InvalidValue(format!(
                "{} transformations for {} maps",
                transformations.len(),
                feature_maps.len()
            )));
        }
        let transformer = MapAlignmentTransformer {
            store_original_rt: true,
            ..MapAlignmentTransformer::default()
        };
        for (map, transformation) in feature_maps.iter_mut().zip(transformations) {
            transformer.transform_feature_map(map, transformation)?;
        }
        Ok(())
    }
}

/// The similarity of two maps from the retention times of their shared
/// sequences (source functor `PeptideIdentificationsPearsonDistance_`, the
/// similarity `ClusterHierarchical` clusters on).
///
/// Two maps without sequences have similarity 0. Otherwise both sequence
/// lists are walked in order together, and for each shared sequence the
/// "median" of its retention times in either map is taken as the source
/// takes it, `Math::median(..., sorted = true)` on the unsorted list: the
/// middle element, or the mean of the two middle ones, by position. The
/// Pearson correlation of the two median lists, in `float`, is multiplied by
/// the number of shared sequences and divided by the number of steps of the
/// walk, which the source calls the union size: the walk stops at the end of
/// either list, so sequences after the last one of the shorter list are not
/// counted.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] when the maps share no sequence but are
/// not both empty (the source's `Math::pearsonCorrelationCoefficient` throws
/// `Exception::InvalidRange` for the empty range).
pub fn peptide_identifications_pearson_similarity(
    map_first: &SeqAndRtList,
    map_second: &SeqAndRtList,
) -> Result<f32> {
    if map_first.len() + map_second.len() == 0 {
        return Ok(0.0);
    }
    let mut first = map_first.iter().peekable();
    let mut second = map_second.iter().peekable();
    let mut rts_first = Vec::new();
    let mut rts_second = Vec::new();
    let mut union_size: f32 = 0.0;
    while let (Some(&(a, a_rts)), Some(&(b, b_rts))) = (first.peek(), second.peek()) {
        if a < b {
            first.next();
        } else if b < a {
            second.next();
        } else {
            rts_first.push(positional_median(a_rts)?);
            rts_second.push(positional_median(b_rts)?);
            first.next();
            second.next();
        }
        union_size += 1.0;
    }
    let intercept_size = rts_first.len();
    let pearson = pearson_correlation_coefficient(&rts_first, &rts_second)? as f32;
    Ok(pearson * intercept_size as f32 / union_size)
}

/// `Math::median(begin, end, true)`: the middle value by position, or the
/// mean of the two middle values, without sorting.
fn positional_median(values: &[f64]) -> Result<f64> {
    let size = values.len();
    if size == 0 {
        return Err(Error::InvalidRange("median of an empty range".into()));
    }
    Ok(if size % 2 == 0 {
        (values[size / 2 - 1] + values[size / 2]) / 2.0
    } else {
        values[(size - 1) / 2]
    })
}

/// The 10/90 percentile range of a sorted, non-empty list as the source
/// indexes it: `v[size * 0.9] - v[size * 0.1]`, the indices truncated.
fn percentile_range(values: &[f64]) -> f64 {
    let n = values.len() as f64;
    values[(n * 0.9) as usize] - values[(n * 0.1) as usize]
}

/// Whether a feature or one of its subordinates refers to an identification
/// graph.
fn refers_to_graph(feature: &Feature) -> bool {
    feature.primary_id.is_some()
        || !feature.id_matches.is_empty()
        || feature.subordinates.iter().any(refers_to_graph)
}

/// `FeatureMap::operator+=`: the document identifier, file path and type
/// and the map's unique ID are reset; proteins, unassigned identifications,
/// processing records and features of `other` are appended; if the features'
/// assigned unique IDs then collide, every feature without one gets a new
/// one and every repeated one is replaced, in order. Returns the info lines
/// the source logs.
fn append_feature_map(
    map: &mut FeatureMap,
    other: &FeatureMap,
    generator: &mut UniqueIdGenerator,
) -> Result<Vec<String>> {
    if other.features.iter().any(refers_to_graph) {
        return Err(Error::InvalidValue(
            "cannot merge features that refer to an identification graph".into(),
        ));
    }
    let mut info = Vec::new();
    if !map.identifier.is_empty() || !other.identifier.is_empty() {
        info.push("DocumentIdentifiers are lost during merge of FeatureMaps".into());
    }
    map.identifier.clear();
    map.loaded_file_path.clear();
    map.loaded_file_type = FileType::Unknown;
    map.unique_id = 0;
    map.protein_identifications
        .extend_from_slice(&other.protein_identifications);
    map.unassigned_peptide_identifications
        .extend_from_slice(&other.unassigned_peptide_identifications);
    map.data_processing
        .extend_from_slice(&other.data_processing);
    map.features.extend_from_slice(&other.features);
    let mut seen = BTreeSet::new();
    let conflict = map
        .features
        .iter()
        .any(|f| f.unique_id != 0 && !seen.insert(f.unique_id));
    if conflict {
        let mut assigned = BTreeSet::new();
        let mut replaced = 0usize;
        for feature in &mut map.features {
            if feature.unique_id == 0 {
                feature.base.unique_id = generator.get_unique_id();
            }
            while assigned.contains(&feature.unique_id) {
                feature.base.unique_id = generator.get_unique_id();
                replaced += 1;
            }
            assigned.insert(feature.unique_id);
        }
        info.push(format!("Replaced {replaced} invalid uniqueID's"));
    }
    Ok(info)
}
