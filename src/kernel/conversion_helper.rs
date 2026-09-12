// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Timo Sachsenberg, OpenMS Rust contributors $

//! Conversions between the three container types (`KERNEL/ConversionHelper.h`).
//!
//! Ports `MapConversion`, whose three static overloads move data between a
//! peak map ([`MSExperiment`](crate::kernel::MSExperiment)), a
//! [`FeatureMap`](crate::kernel::features::FeatureMap) and a
//! [`ConsensusMap`](crate::kernel::features::ConsensusMap). The overload set
//! becomes three named associated functions, as the port's naming rule
//! requires, and each returns the new container instead of filling an output
//! parameter — so a failure cannot leave the caller with a half-written map.
//! See `docs/CONVERSION_HELPER_SUPPORT.md`.
//!
//! The conversions are deliberately asymmetric; which unique IDs and which
//! identification records survive depends on the direction, and each function
//! documents its own contract.

use super::features::{
    BaseFeature, ConsensusFeature, ConsensusMap, Feature, FeatureHandle, FeatureMap,
};
use super::{MSExperiment, Peak2D};
use crate::concept::UniqueIdGenerator;
use crate::{Error, Result};
use std::cmp::Ordering;

fn limit() -> Error {
    Error::InvalidValue("map conversion work limit exceeded".into())
}

/// Conversions between the peak, feature and consensus containers.
///
/// Ports the `MapConversion` class. It holds no state; the source class is a
/// bag of static member functions and this is a unit struct for the same
/// reason.
///
/// Each conversion replaces the destination wholesale and leaves it in a state
/// where the range queries of
/// [`FeatureMap::ranges`](crate::kernel::features::FeatureMap::ranges) and
/// [`ConsensusMap::ranges`](crate::kernel::features::ConsensusMap::ranges)
/// reflect the new contents without further bookkeeping. The source calls
/// `updateRanges()` on the result to achieve that; here ranges are computed on
/// demand, so the guarantee holds with no cache to refresh.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MapConversion;

impl MapConversion {
    /// Elements one conversion will copy.
    pub const MAX_ITEMS: usize = 10_000_000;

    /// Copy the most intense peaks of a peak map into a consensus map.
    ///
    /// Ports `MapConversion::convert(UInt64, PeakMap&, ConsensusMap&, Size)`.
    /// The `n` most intense peaks of `input_map` become `ConsensusFeature`
    /// entries tagged with `input_map_index`, in descending intensity order,
    /// each holding one handle whose unique ID is the peak's position in that
    /// order. The result is given a fresh container unique ID, because a peak
    /// map has none of its own, and the column header for `input_map_index`
    /// records how many peaks were written.
    ///
    /// Only MS level 1 peaks are considered, because the source reads them
    /// through `MSExperiment::get2DData`, which skips every other level.
    ///
    /// # Arguments
    ///
    /// * `input_map_index` — index assigned to the input in the resulting
    ///   column headers.
    /// * `input_map` — the source peaks. The source takes this by mutable
    ///   reference only to call `updateRanges()` on it, which this port does
    ///   not need, so the input is borrowed immutably and is not modified.
    /// * `n` — how many peaks to copy at most; `None` replaces the source's
    ///   `Size(-1)` default and keeps every peak.
    /// * `generator` — caller-owned replacement for the source's process-wide
    ///   `UniqueIdGenerator` singleton.
    ///
    /// # Notes
    ///
    /// The source clamps `n` against `MSExperiment::getSize()`, which counts
    /// the peaks of *every* spectrum plus every chromatogram point, and then
    /// partially sorts the MS1-only vector `get2DData` produced using that
    /// count as the middle iterator. With any MS2 spectrum or chromatogram
    /// present — the default `Size(-1)` reaches this path every time — the
    /// middle iterator is past the end and the loop reads past the end of the
    /// vector. This port clamps against the number of points actually
    /// collected, which is the evident intent.
    ///
    /// The source uses `std::partial_sort`, which is not stable, so peaks of
    /// equal intensity may appear in any order and the last of several equally
    /// intense peaks at the cut-off is chosen arbitrarily. This port sorts
    /// stably, so ties keep spectrum-then-peak order and the result is
    /// reproducible.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the input exceeds the ceilings of
    /// [`MSExperiment::get_2d_data`](crate::kernel::MSExperiment::get_2d_data)
    /// or holds a non-finite coordinate, and when more than
    /// [`MapConversion::MAX_ITEMS`] peaks would be copied.
    pub fn peak_map_to_consensus(
        input_map_index: u64,
        input_map: &MSExperiment,
        n: Option<usize>,
        generator: &mut UniqueIdGenerator,
    ) -> Result<ConsensusMap> {
        let mut points = input_map.get_2d_data()?;
        let count = n.unwrap_or(points.len()).min(points.len());
        if count > Self::MAX_ITEMS {
            return Err(limit());
        }
        points.sort_by(|left, right| {
            right
                .intensity
                .partial_cmp(&left.intensity)
                .unwrap_or(Ordering::Equal)
        });
        let mut output = ConsensusMap::new();
        output.unique_id = generator.get_unique_id();
        output.features.reserve(count);
        for (element_index, point) in points.into_iter().take(count).enumerate() {
            let index = u64::try_from(element_index).map_err(|_| limit())?;
            output
                .features
                .push(consensus_from_peak(input_map_index, point, index)?);
        }
        output
            .column_headers
            .entry(input_map_index)
            .or_default()
            .size = count;
        Ok(output)
    }

    /// Convert a consensus map to a feature map.
    ///
    /// Ports `MapConversion::convert(ConsensusMap const&, bool, FeatureMap&)`.
    /// Every consensus feature becomes a [`Feature`] carrying the whole
    /// `BaseFeature` part — position, intensity, quality, charge, width, meta
    /// values, peptide identifications and identification-graph references. The
    /// grouped feature handles and the quantification ratios are dropped, and
    /// the produced features have no convex hulls, no per-dimension qualities
    /// and no subordinates, because the source assigns only the `BaseFeature`
    /// slice onto default-constructed features.
    ///
    /// The document identifier and the protein and unassigned peptide
    /// identifications are carried over.
    ///
    /// # Arguments
    ///
    /// * `input_map` — the source consensus map. It is only read.
    /// * `keep_uids` — when `true`, the container's unique ID and every
    ///   element's unique ID are preserved; when `false`, all of them are drawn
    ///   fresh from `generator`.
    /// * `generator` — caller-owned replacement for the source's process-wide
    ///   `UniqueIdGenerator` singleton. It is untouched when `keep_uids` is
    ///   `true`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the input holds more than
    /// [`MapConversion::MAX_ITEMS`] consensus features. The source is
    /// unbounded.
    pub fn consensus_to_feature_map(
        input_map: &ConsensusMap,
        keep_uids: bool,
        generator: &mut UniqueIdGenerator,
    ) -> Result<FeatureMap> {
        if input_map.features.len() > Self::MAX_ITEMS {
            return Err(limit());
        }
        let mut output = FeatureMap::new();
        output.identifier.clone_from(&input_map.identifier);
        output
            .loaded_file_path
            .clone_from(&input_map.loaded_file_path);
        output.loaded_file_type = input_map.loaded_file_type;
        output.unique_id = if keep_uids {
            input_map.unique_id
        } else {
            generator.get_unique_id()
        };
        output
            .protein_identifications
            .clone_from(&input_map.protein_identifications);
        output
            .unassigned_peptide_identifications
            .clone_from(&input_map.unassigned_peptide_identifications);
        output.features.reserve(input_map.features.len());
        for consensus in &input_map.features {
            let mut feature = Feature::from(consensus.base.clone());
            if !keep_uids {
                feature.base.unique_id = generator.get_unique_id();
            }
            output.features.push(feature);
        }
        Ok(output)
    }

    /// Convert a feature map to a consensus map.
    ///
    /// Ports `MapConversion::convert(UInt64, FeatureMap const&, ConsensusMap&,
    /// Size)`. The first `n` features, **in input order and without any
    /// sorting**, become `ConsensusFeature` entries tagged with
    /// `input_map_index`, each holding one handle that references the feature
    /// by its own unique ID. Each copied feature's peptide identifications are
    /// stamped with a `map_index` meta value naming `input_map_index`, as the
    /// source's `ConsensusFeature(UInt64, const BaseFeature&)` constructor
    /// does. The protein and unassigned peptide identifications are carried
    /// over.
    ///
    /// # Arguments
    ///
    /// * `input_map_index` — index assigned to the input in the resulting
    ///   column headers.
    /// * `input_map` — the source feature map. It is only read.
    /// * `n` — how many features to copy at most; `None` replaces the source's
    ///   `Size(-1)` default and keeps every feature.
    ///
    /// # Notes
    ///
    /// Because features are taken in input order, `n` is only useful after
    /// pre-sorting the input — for example with
    /// [`FeatureMap::sort_by_intensity`](crate::kernel::features::FeatureMap::sort_by_intensity);
    /// the parameter exists for symmetry with
    /// [`MapConversion::peak_map_to_consensus`], which does sort.
    ///
    /// The column header `size` for `input_map_index` is set to the **full**
    /// input size even when `n` truncates the copy, unlike the peak-map
    /// conversion, which records what it wrote. Inspect the returned map's own
    /// length for the number of features actually written. The source
    /// documents this asymmetry and it is preserved.
    ///
    /// The result's container unique ID is taken from `input_map` rather than
    /// drawn fresh, which the source marks as an arguable design decision;
    /// overwrite it afterwards when a new identity is wanted.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the input holds more than
    /// [`MapConversion::MAX_ITEMS`] features, when a copied feature is invalid
    /// — a non-finite coordinate or a negative width, which
    /// [`ConsensusFeature::from_feature`](crate::kernel::features::ConsensusFeature::from_feature)
    /// checks and the source does not — or when `input_map_index` exceeds the
    /// signed metadata integer range used for the `map_index` stamp.
    pub fn feature_map_to_consensus(
        input_map_index: u64,
        input_map: &FeatureMap,
        n: Option<usize>,
    ) -> Result<ConsensusMap> {
        if input_map.features.len() > Self::MAX_ITEMS {
            return Err(limit());
        }
        let count = n
            .unwrap_or(input_map.features.len())
            .min(input_map.features.len());
        let mut output = ConsensusMap::new();
        output.unique_id = input_map.unique_id;
        output.features.reserve(count);
        for feature in &input_map.features[..count] {
            let stamped = feature.base.clone_with_map_index(input_map_index)?;
            output
                .features
                .push(ConsensusFeature::from_feature(input_map_index, &stamped)?);
        }
        output
            .column_headers
            .entry(input_map_index)
            .or_default()
            .size = input_map.features.len();
        output
            .protein_identifications
            .clone_from(&input_map.protein_identifications);
        output
            .unassigned_peptide_identifications
            .clone_from(&input_map.unassigned_peptide_identifications);
        Ok(output)
    }
}

/// Source `ConsensusFeature(UInt64 map_index, const Peak2D& element, UInt64
/// element_index)`: the peak becomes the consensus summary and one handle
/// identified by `element_index`. Written here rather than as a constructor on
/// [`ConsensusFeature`](crate::kernel::features::ConsensusFeature), because
/// `KERNEL/ConsensusFeature.h` belongs to a different work package; the two
/// statements below are exactly what that constructor performs.
fn consensus_from_peak(
    map_index: u64,
    point: Peak2D,
    element_index: u64,
) -> Result<ConsensusFeature> {
    let mut result = ConsensusFeature::from(BaseFeature {
        rt: point.rt(),
        mz: point.mz(),
        intensity: point.intensity,
        ..BaseFeature::default()
    });
    result.insert(FeatureHandle::from_peak(map_index, point, element_index))?;
    Ok(result)
}
