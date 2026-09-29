// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Normalisation of a consensus map by averaging per-feature intensity ratios
//! against a reference map.
//!
//! Native equivalent of
//! `ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmThreshold.h`, the
//! `robust_regression` algorithm of the `ConsensusMapNormalizer` tool. See
//! `docs/CONSENSUS_MAP_NORMALIZER_SUPPORT.md`.

use super::consensus_map_normalizer_median::{
    IdentificationFilter, check_maps, largest_map, map_slot, rewrite_intensities,
};
use crate::error::{Error, Result};
use crate::kernel::features::ConsensusMap;
use crate::math::statistic_functions;

/// Largest `number_of_maps * number_of_consensus_features` that
/// [`ConsensusMapNormalizerAlgorithmThreshold::compute_correlation`] accepts:
/// the source stores a dense intensity matrix of that many doubles. Checked
/// before the matrix is allocated; the source has no ceiling.
pub const MAX_MATRIX_CELLS: usize = 1 << 27;

/// The normalisation factors and how they were obtained (the result of source
/// `computeCorrelation`, which returns only the factors and logs the rest).
#[derive(Clone, Debug, PartialEq)]
pub struct ThresholdFactors {
    /// One factor per map, or all `1.0` when some map kept no ratio.
    pub ratios: Vec<f64>,
    /// Index of the reference map: largest column-header `size`, first on ties.
    pub reference_map: usize,
    /// Number of consensus features that passed the filters.
    pub passing_features: usize,
    /// `false` when some map kept no ratio inside the window; the source then
    /// logs a warning and returns all `1.0`, i.e. no normalisation.
    pub enough_features: bool,
}

/// Ratio-threshold normalisation (source
/// `ConsensusMapNormalizerAlgorithmThreshold`, a stateless class whose copy and
/// assignment are suppressed).
#[derive(Clone, Copy, Debug, Default)]
pub struct ConsensusMapNormalizerAlgorithmThreshold;

impl ConsensusMapNormalizerAlgorithmThreshold {
    /// Compute one normalisation factor per map (source `computeCorrelation`).
    ///
    /// The reference map is the column header with the largest `size`. For
    /// each map `j` and each consensus feature that passes `filter` and has a
    /// non-zero intensity in both the reference and `j`, the ratio
    /// `reference / intensity_j` is kept when it lies strictly inside
    /// `(ratio_threshold, 1 / ratio_threshold)`; the factor for `j` is the mean
    /// of the kept ratios. A map's intensity in a feature is its handle's, or
    /// zero without a handle, so a handle of intensity zero counts as absent.
    /// The `ConsensusMapNormalizer` tool restricts `ratio_threshold` to
    /// `[0.001, 1.0]` and sorts the map by size first; neither is required
    /// here, and any value is used as given.
    ///
    /// The mean is the sum in feature order divided by the count, as the
    /// source's `Math::mean`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a column header in `0..n` is
    /// missing or a handle names a map index outside `0..n` (both undefined
    /// behaviour in the source, which dereferences the missing header or
    /// indexes past its matrix), when the handle count exceeds
    /// [`MAX_HANDLES`](super::consensus_map_normalizer_median::MAX_HANDLES) or
    /// the matrix exceeds [`MAX_MATRIX_CELLS`], and the
    /// [`IdentificationFilter::passes`] error.
    pub fn compute_correlation(
        map: &ConsensusMap,
        ratio_threshold: f64,
        filter: &IdentificationFilter,
    ) -> Result<ThresholdFactors> {
        let number_of_maps = check_maps(map)?;
        let number_of_features = map.features.len();
        if number_of_maps.saturating_mul(number_of_features) > MAX_MATRIX_CELLS {
            return Err(Error::InvalidValue(format!(
                "{number_of_maps} maps x {number_of_features} consensus features exceed the {MAX_MATRIX_CELLS}-cell ceiling"
            )));
        }
        let reference_map = largest_map(map, number_of_maps);
        let mut feature_int = vec![vec![0.0f64; number_of_features]; number_of_maps];
        let mut passing_features = 0;
        for (index, feature) in map.features.iter().enumerate() {
            if !filter.passes(feature, map)? {
                continue;
            }
            passing_features += 1;
            for handle in feature.handles() {
                feature_int[map_slot(handle)][index] = f64::from(handle.intensity);
            }
        }

        let mut ratios = Vec::with_capacity(number_of_maps);
        for j in 0..number_of_maps {
            let mut kept = Vec::new();
            for k in 0..number_of_features {
                let reference = feature_int[reference_map][k];
                let value = feature_int[j][k];
                if reference != 0.0 && value != 0.0 {
                    let ratio = reference / value;
                    if ratio > ratio_threshold && ratio < 1.0 / ratio_threshold {
                        kept.push(ratio);
                    }
                }
            }
            if kept.is_empty() {
                return Ok(ThresholdFactors {
                    ratios: vec![1.0; number_of_maps],
                    reference_map,
                    passing_features,
                    enough_features: false,
                });
            }
            ratios.push(statistic_functions::mean(&kept)?);
        }
        Ok(ThresholdFactors {
            ratios,
            reference_map,
            passing_features,
            enough_features: true,
        })
    }

    /// Multiply every handle's intensity by its map's factor (source
    /// `normalizeMaps`). The product is taken in `f64` and stored as `f32`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a handle's map index has no entry
    /// in `ratios` (the source reads past the end), or when a product is not
    /// finite. On any error the map is unchanged.
    pub fn normalize_maps(map: &mut ConsensusMap, ratios: &[f64]) -> Result<()> {
        for feature in &map.features {
            for handle in feature.handles() {
                if usize::try_from(handle.map_index).map_or(true, |index| index >= ratios.len()) {
                    return Err(Error::InvalidValue(format!(
                        "no normalisation factor for map index {}",
                        handle.map_index
                    )));
                }
            }
        }
        rewrite_intensities(map, |handle| {
            f64::from(handle.intensity) * ratios[handle.map_index as usize]
        })
    }
}
