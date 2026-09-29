// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The interface of all feature grouping algorithms, which group
//! corresponding features in one map or across maps.
//!
//! Native equivalent of `ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h`.
//! See `docs/FEATURE_GROUPING_SUPPORT.md`.
//!
//! The source is an abstract `DefaultParamHandler` with a pure virtual
//! `group(vector<FeatureMap>, ConsensusMap&)`, a virtual consensus-map
//! overload that converts to feature maps and forwards, the public
//! `transferSubelements`, and the protected template `postprocess_`. Here the
//! two `group` overloads are the trait methods [`FeatureGroupingAlgorithm::group`](crate::analysis::mapmatching::feature_grouping_algorithm::FeatureGroupingAlgorithm::group)
//! and [`FeatureGroupingAlgorithm::group_consensus`](crate::analysis::mapmatching::feature_grouping_algorithm::FeatureGroupingAlgorithm::group_consensus), and the other two are
//! free functions.

use std::cmp::Ordering;
use std::collections::BTreeMap;

use crate::concept::UniqueIdGenerator;
use crate::error::{Error, Result};
use crate::identification::{PeptideIdentification, ProteinIdentification};
use crate::kernel::conversion_helper::MapConversion;
use crate::kernel::features::{ConsensusFeature, ConsensusMap, FeatureHandle, FeatureMap};
use crate::metadata::MetaValue;

/// Meta value naming the input map of a peptide identification.
const MAP_INDEX: &str = "map_index";
/// Meta value keeping a peptide identification's map index from before a
/// consensus-map linking step.
const OLD_MAP_INDEX: &str = "old_map_index";

/// A feature grouping algorithm (source `FeatureGroupingAlgorithm`).
pub trait FeatureGroupingAlgorithm {
    /// Group the features of `maps`, writing the consensus features to `out`
    /// (source `group(const std::vector<FeatureMap>&, ConsensusMap&)`).
    ///
    /// # Errors
    ///
    /// Implementation-specific; typically [`Error::InvalidValue`] for invalid
    /// input (source `Exception::IllegalArgument`).
    fn group(&mut self, maps: &[FeatureMap], out: &mut ConsensusMap) -> Result<()>;

    /// Group the consensus features of `maps` (source
    /// `group(const std::vector<ConsensusMap>&, ConsensusMap&)`).
    ///
    /// The default converts each map to a feature map, keeping unique IDs
    /// (source `MapConversion::convert(map, true, fm)`), and forwards to
    /// [`Self::group`]; the sub-features are then lost, and
    /// [`transfer_subelements`] restores them. The source logs a warning
    /// that consensus maps are not supported directly; algorithms that support
    /// them override this.
    ///
    /// # Errors
    ///
    /// The [`MapConversion::consensus_to_feature_map`] error, and those of
    /// [`Self::group`].
    fn group_consensus(&mut self, maps: &[ConsensusMap], out: &mut ConsensusMap) -> Result<()> {
        let mut generator = UniqueIdGenerator::new();
        let feature_maps = maps
            .iter()
            .map(|map| MapConversion::consensus_to_feature_map(map, true, &mut generator))
            .collect::<Result<Vec<_>>>()?;
        self.group(&feature_maps, out)
    }
}

/// Transfer the grouped sub-features from the input consensus maps to the
/// result (source `transferSubelements`).
///
/// The column headers of `out` are replaced by those of all `maps`, renumbered
/// consecutively in input order (map `i`'s header `k` becomes the next free
/// index). Every consensus feature of `out` keeps its own summary, peptide
/// identifications and metadata, loses its handles and ratios, and receives
/// the handles of the consensus features its handles referenced — handle map
/// index `i` naming input map `i`, handle unique ID naming that map's
/// consensus feature — each renumbered to its new column.
///
/// Peptide identifications (of the consensus features and the unassigned
/// ones) that carry `old_map_index` get `map_index` renumbered from the pair
/// (`map_index`, `old_map_index`) and lose `old_map_index`; those without it
/// lose `map_index`, since there was none to correct.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a handle of `out` names an input map
/// that does not exist or a unique ID not found in that map (the source
/// dereferences a default-constructed iterator), when a referenced handle's
/// map index has no column header in its map, when an identification's
/// `old_map_index` pair has no column (the source's `operator[]` would add
/// one), when a meta value is not a non-negative integer, or when two
/// transferred handles collide (source `Exception::InvalidValue`). On error
/// `out` is unchanged.
pub fn transfer_subelements(maps: &[ConsensusMap], out: &mut ConsensusMap) -> Result<()> {
    let mut headers = BTreeMap::new();
    let mut mapid_table = BTreeMap::<(u64, u64), u64>::new();
    for (i, consensus) in maps.iter().enumerate() {
        for (id, header) in &consensus.column_headers {
            let counter = mapid_table.len() as u64;
            mapid_table.insert((i as u64, *id), counter);
            headers.insert(counter, header.clone());
        }
    }

    let mut lookup: Vec<BTreeMap<u64, &ConsensusFeature>> = Vec::with_capacity(maps.len());
    for consensus in maps {
        let mut table = BTreeMap::new();
        for feature in &consensus.features {
            // Source `map::insert` keeps the first of equal keys.
            table.entry(feature.unique_id).or_insert(feature);
        }
        lookup.push(table);
    }

    let mut features = Vec::with_capacity(out.features.len());
    for feature in &out.features {
        let mut adjusted = ConsensusFeature::from(feature.base.clone());
        for sub in feature.handles() {
            let origin = usize::try_from(sub.map_index)
                .ok()
                .and_then(|index| lookup.get(index))
                .and_then(|table| table.get(&sub.unique_id))
                .ok_or_else(|| {
                    Error::InvalidValue(format!(
                        "no consensus feature {} in input map {}",
                        sub.unique_id, sub.map_index
                    ))
                })?;
            for handle in origin.handles() {
                let new_id = *mapid_table
                    .get(&(sub.map_index, handle.map_index))
                    .ok_or_else(|| {
                        Error::InvalidValue(format!(
                            "input map {} has no column header {}",
                            sub.map_index, handle.map_index
                        ))
                    })?;
                adjusted.insert(FeatureHandle {
                    map_index: new_id,
                    ..*handle
                })?;
            }
        }
        for identification in &mut adjusted.peptide_identifications {
            reindex(identification, &mapid_table)?;
        }
        features.push(adjusted);
    }
    let mut unassigned = out.unassigned_peptide_identifications.clone();
    for identification in &mut unassigned {
        reindex(identification, &mapid_table)?;
    }

    out.column_headers = headers;
    out.features = features;
    out.unassigned_peptide_identifications = unassigned;
    Ok(())
}

/// The `old_map_index` correction of [`transfer_subelements`] for one
/// identification.
fn reindex(
    identification: &mut PeptideIdentification,
    mapid_table: &BTreeMap<(u64, u64), u64>,
) -> Result<()> {
    let meta = &mut identification.metadata;
    let Some(old) = meta.get(OLD_MAP_INDEX) else {
        meta.remove(MAP_INDEX);
        return Ok(());
    };
    let old = index_value(old)?;
    let file = index_value(meta.get(MAP_INDEX).ok_or_else(|| {
        Error::InvalidValue("peptide identification has old_map_index but no map_index".into())
    })?)?;
    let new_index = *mapid_table.get(&(file, old)).ok_or_else(|| {
        Error::InvalidValue(format!(
            "no column for input map {file} and old map index {old}"
        ))
    })?;
    meta.insert(MAP_INDEX.into(), MetaValue::try_from(new_index)?);
    meta.remove(OLD_MAP_INDEX);
    Ok(())
}

fn index_value(value: &MetaValue) -> Result<u64> {
    u64::try_from(value.as_i64()?)
        .map_err(|_| Error::InvalidValue("map index meta value must not be negative".into()))
}

/// The protein and unassigned peptide identifications of one input map, as
/// [`postprocess`] reads them from a feature map or a consensus map.
pub type InputIdentifications<'a> = (&'a [ProteinIdentification], &'a [PeptideIdentification]);

/// The identifications of a feature map, for [`postprocess`].
pub fn feature_map_identifications(map: &FeatureMap) -> InputIdentifications<'_> {
    (
        &map.protein_identifications,
        &map.unassigned_peptide_identifications,
    )
}

/// The identifications of a consensus map, for [`postprocess`].
pub fn consensus_map_identifications(map: &ConsensusMap) -> InputIdentifications<'_> {
    (
        &map.protein_identifications,
        &map.unassigned_peptide_identifications,
    )
}

/// Post-process a grouping result (source `postprocess_`).
///
/// Appends, map by map in input order, each input's protein identifications
/// and its unassigned peptide identifications, the latter stamped with
/// `map_index` = the input's position (identifications assigned to features
/// are re-indexed by the algorithms themselves). Then sorts canonically:
/// stably by quality ascending, then by the handle sets (lexicographic in
/// `(map index, unique id)`), then by size descending — three stable sorts in
/// that order, so size is the primary key.
///
/// The source sorts through `ConsensusMap::sortByQuality`, `sortByMaps` and
/// `sortBySize`; the same stable orders are applied here directly, because
/// the Rust methods also validate the whole map and a grouping result may
/// legitimately hold two singletons with the same unique ID from different
/// inputs.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when an input position does not fit the
/// meta-value integer range.
pub fn postprocess(inputs: &[InputIdentifications<'_>], out: &mut ConsensusMap) -> Result<()> {
    let mut proteins = out.protein_identifications.clone();
    let mut unassigned = out.unassigned_peptide_identifications.clone();
    for (map_idx, (protein_ids, peptide_ids)) in inputs.iter().enumerate() {
        proteins.extend(protein_ids.iter().cloned());
        let value = MetaValue::try_from(map_idx as u64)?;
        for identification in peptide_ids.iter() {
            let mut copy = identification.clone();
            copy.metadata.insert(MAP_INDEX.into(), value.clone());
            unassigned.push(copy);
        }
    }
    out.protein_identifications = proteins;
    out.unassigned_peptide_identifications = unassigned;
    canonical_sort(&mut out.features);
    Ok(())
}

/// `sortByQuality(); sortByMaps(); sortBySize();` as stable sorts.
pub(crate) fn canonical_sort(features: &mut [ConsensusFeature]) {
    features.sort_by(|a, b| a.quality.partial_cmp(&b.quality).unwrap_or(Ordering::Equal));
    features.sort_by(|a, b| {
        a.handles()
            .iter()
            .map(FeatureHandle::key)
            .cmp(b.handles().iter().map(FeatureHandle::key))
    });
    features.sort_by_key(|feature| std::cmp::Reverse(feature.len()));
}
