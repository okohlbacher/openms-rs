// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Matching MS2 spectra to LC-MS features by precursor (m/z, RT).
//!
//! Native equivalent of `ANALYSIS/MAPMATCHING/FeatureMapping.h`. See
//! `docs/FEATURE_GROUPING_SUPPORT.md`.
//!
//! The source class is a namespace holding two record types and one static
//! algorithm: [`FeatureMappingInfo`] (the feature maps and a
//! [`KDTreeFeatureMaps`] over them), [`FeatureToMs2Indices`] (the spectra
//! grouped by feature, and the unassigned ones) and
//! [`assign_ms2_index_to_feature`]. A typical caller builds the info once and
//! assigns per MS run; workflows that use this include SIRIUS export and other
//! "annotate MS2 with quantified feature" pipelines.
//!
//! # The map-1 exclusion of the source
//!
//! The source documents that features from any map are eligible, "the
//! kd-tree query is run with `include_features_from_same_map = true`", but it
//! calls `queryRegion(..., matches, true)`, whose sixth parameter is the
//! `Size ignored_map_index`: `true` converts to 1, so the features of map 1
//! are **never** matched, and those of every other map are. This port
//! reproduces the executed behaviour (oracle `M2-labeled`, case
//! `featuremapping`) and reports the discrepancy as a source defect.

use std::collections::BTreeMap;

use super::kd_tree_feature_maps::KDTreeFeatureMaps;
use crate::concept::math_functions::tolerance_window;
use crate::error::Result;
use crate::kernel::MSExperiment;
use crate::kernel::features::FeatureMap;

/// The map index the source's `queryRegion(..., true)` excludes.
const SOURCE_IGNORED_MAP_INDEX: usize = 1;

/// The input bundle: feature maps and a kd-tree over them (source
/// `FeatureMapping::FeatureMappingInfo`).
///
/// The source owns both, and the caller must build the tree from the stored
/// maps; here the tree borrows the maps, which the lifetime enforces. Build
/// it with [`FeatureMappingInfo::new`].
#[derive(Clone, Debug, Default)]
pub struct FeatureMappingInfo<'a> {
    /// The feature data, one map per input run.
    pub feature_maps: &'a [FeatureMap],
    /// The (RT, m/z) kd-tree over the features of `feature_maps`.
    pub kd_tree: KDTreeFeatureMaps<'a>,
}

impl<'a> FeatureMappingInfo<'a> {
    /// The maps and a kd-tree built from them with
    /// [`KDTreeFeatureMaps::add_maps`].
    ///
    /// # Errors
    ///
    /// The [`KDTreeFeatureMaps::add_maps`] errors.
    pub fn new(feature_maps: &'a [FeatureMap]) -> Result<Self> {
        Ok(Self {
            feature_maps,
            kd_tree: KDTreeFeatureMaps::from_maps(feature_maps)?,
        })
    }
}

/// The output bundle (source `FeatureMapping::FeatureToMs2Indices`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FeatureToMs2Indices {
    /// MS2 spectrum indices grouped by the feature they were assigned to.
    ///
    /// The source keys this by the feature's address, so its iteration order
    /// follows memory layout; here the key is the feature's global index in
    /// the [`FeatureMappingInfo::kd_tree`]
    /// ([`KDTreeFeatureMaps::feature`] returns the feature).
    pub assigned_ms2: BTreeMap<usize, Vec<usize>>,
    /// Indices of MS2 spectra that had a precursor but no feature inside the
    /// tolerance window.
    pub unassigned_ms2: Vec<usize>,
}

/// Assign each MS2 spectrum to the feature whose position matches its first
/// precursor best (source `assignMS2IndexToFeature`).
///
/// For every spectrum of `spectra` in order: spectra of MS level other than 2
/// are skipped; spectra without a precursor are dropped (they appear in
/// neither collection); the others are matched against the features in the
/// window RT `± precursor_rt_tolerance` around the spectrum's RT and the m/z
/// window of [`tolerance_window`] (ppm when `ppm`) around the first
/// precursor's m/z. Without a feature in the window the spectrum is
/// unassigned; otherwise it goes to the feature with the smallest absolute
/// m/z difference, the first in the kd-tree's query order on ties (RT is not
/// a tie-breaker), starting from a distance of 1e11. Features of map 1 are
/// never eligible; see the module documentation.
///
/// # Errors
///
/// The [`tolerance_window`] errors for a non-finite window, which the source
/// computes regardless.
pub fn assign_ms2_index_to_feature(
    spectra: &MSExperiment,
    fm_info: &FeatureMappingInfo<'_>,
    precursor_mz_tolerance: f64,
    precursor_rt_tolerance: f64,
    ppm: bool,
) -> Result<FeatureToMs2Indices> {
    let mut result = FeatureToMs2Indices::default();
    let mut matches = Vec::new();
    for (index, spectrum) in spectra.spectra.iter().enumerate() {
        if spectrum.ms_level != 2 {
            continue;
        }
        let Some(precursor) = spectrum.precursors.first() else {
            continue;
        };
        let mz = precursor.mz;
        let rt = spectrum.rt;
        let (mz_low, mz_high) = tolerance_window(mz, precursor_mz_tolerance, ppm)?;
        fm_info.kd_tree.query_region(
            rt - precursor_rt_tolerance,
            rt + precursor_rt_tolerance,
            mz_low,
            mz_high,
            &mut matches,
            Some(SOURCE_IGNORED_MAP_INDEX),
        );
        if matches.is_empty() {
            result.unassigned_ms2.push(index);
            continue;
        }
        let mut best = 0usize;
        let mut min_distance = 1e11;
        for &k in &matches {
            let distance = (fm_info.kd_tree.mz(k)? - mz).abs();
            if distance < min_distance {
                min_distance = distance;
                best = k;
            }
        }
        result.assigned_ms2.entry(best).or_default().push(index);
    }
    Ok(result)
}
