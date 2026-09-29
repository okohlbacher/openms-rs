// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Quantile normalisation of the maps of a consensus map.
//!
//! Native equivalent of
//! `ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmQuantile.h`. See
//! `docs/CONSENSUS_MAP_NORMALIZER_SUPPORT.md`.
//!
//! Every map's sorted intensity distribution is replaced by a shared reference
//! distribution, the per-quantile average of all maps' distributions, so that
//! all maps end up with the same distribution while the intensity ranks inside
//! each map are preserved.

use super::consensus_map_normalizer_median::{check_maps, map_slot, rewrite_intensities};
use crate::error::{Error, Result};
use crate::kernel::features::ConsensusMap;

/// Quantile normalisation (source `ConsensusMapNormalizerAlgorithmQuantile`,
/// a stateless class whose copy and assignment are suppressed).
#[derive(Clone, Copy, Debug, Default)]
pub struct ConsensusMapNormalizerAlgorithmQuantile;

impl ConsensusMapNormalizerAlgorithmQuantile {
    /// Quantile-normalise the handle intensities of `map` in place (source
    /// `normalizeMaps`).
    ///
    /// Each map's intensities are sorted and [resampled](Self::resample) to
    /// the size of the largest map; the resampled distributions are averaged
    /// per quantile, adding `value / number_of_maps` map by map, into the
    /// reference distribution; the reference is resampled back to each map's
    /// own size; and the `k`-th smallest intensity of a map receives the `k`-th
    /// value of that resampled reference. Equal intensities are ranked by
    /// their position in the map's handle order, as the source's sort of
    /// `(value, index)` pairs ranks them. Values are stored as `f32`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a column header in `0..n` is
    /// missing (source `Exception::ElementNotFound` from
    /// [`Self::extract_intensity_vectors`]), when a handle names a map index
    /// outside `0..n`, when some map has no handle while another has
    /// (the source then reads `front()` of an empty vector, which is undefined
    /// behaviour), when the handle count exceeds
    /// [`MAX_HANDLES`](super::consensus_map_normalizer_median::MAX_HANDLES), and
    /// when a result is not finite. On any error the map is unchanged.
    pub fn normalize_maps(map: &mut ConsensusMap) -> Result<()> {
        let mut feature_ints = Self::extract_intensity_vectors(map)?;
        let number_of_maps = feature_ints.len();
        let largest = feature_ints.iter().map(Vec::len).max().unwrap_or(0);
        if largest > 0 {
            if let Some(empty) = feature_ints.iter().position(Vec::is_empty) {
                return Err(Error::InvalidValue(format!(
                    "quantile normalisation needs at least one feature handle in every map; map {empty} has none"
                )));
            }
        }

        let mut reference = vec![0.0f64; largest];
        for ints in &feature_ints {
            let mut sorted = ints.clone();
            sorted.sort_by(source_less);
            let resampled = Self::resample(&sorted, largest)?;
            for (target, value) in reference.iter_mut().zip(&resampled) {
                *target += value / number_of_maps as f64;
            }
        }

        for ints in &mut feature_ints {
            let normalized = Self::resample(&reference, ints.len())?;
            let mut order: Vec<(f64, usize)> = ints
                .iter()
                .enumerate()
                .map(|(index, &value)| (value, index))
                .collect();
            order.sort_by(|a, b| source_less(&a.0, &b.0).then(a.1.cmp(&b.1)));
            for ((_, index), value) in order.into_iter().zip(normalized) {
                ints[index] = value;
            }
        }

        Self::set_normalized_intensity_values(&feature_ints, map)
    }

    /// Linearly resample `data_in` to `n_resampling_points` evenly spaced
    /// positions (source `resample`, whose output parameter is the returned
    /// vector).
    ///
    /// The first output is `data_in`'s first value and the last output its
    /// last value; position `i` in between reads the pseudo-index
    /// `i * (len - 1) / (n - 1)` and interpolates between its floor and
    /// ceiling with weights `1 - (p - floor)` and `1 - (ceil - p)`. Zero points
    /// give an empty vector. The input is assumed sorted by the caller, as the
    /// source documents; nothing here depends on it.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `data_in` is empty and
    /// `n_resampling_points` is not zero; the source reads `front()` of the
    /// empty vector, which is undefined behaviour. Also when
    /// `n_resampling_points` does not fit the source's 32-bit `UInt`.
    pub fn resample(data_in: &[f64], n_resampling_points: usize) -> Result<Vec<f64>> {
        if n_resampling_points == 0 {
            return Ok(Vec::new());
        }
        if u32::try_from(n_resampling_points).is_err() {
            return Err(Error::InvalidValue(format!(
                "{n_resampling_points} resampling points exceed the source's 32-bit count"
            )));
        }
        let (Some(&first), Some(&last)) = (data_in.first(), data_in.last()) else {
            return Err(Error::InvalidValue(
                "cannot resample an empty distribution".into(),
            ));
        };
        let mut data_out = vec![0.0; n_resampling_points];
        data_out[0] = first;
        data_out[n_resampling_points - 1] = last;
        let delta = (data_in.len() - 1) as f64 / (n_resampling_points - 1) as f64;
        let at = |index: f64| {
            data_in.get(index as usize).copied().ok_or_else(|| {
                Error::InvalidValue("resampling index outside the distribution".into())
            })
        };
        for (i, slot) in data_out
            .iter_mut()
            .enumerate()
            .take(n_resampling_points - 1)
            .skip(1)
        {
            let pseudo_index = i as f64 * delta;
            let left_index = pseudo_index.floor();
            let right_index = pseudo_index.ceil();
            *slot = if left_index == right_index {
                at(left_index)?
            } else {
                let weight_left = 1.0 - (pseudo_index - left_index);
                let weight_right = 1.0 - (right_index - pseudo_index);
                weight_left * at(left_index)? + weight_right * at(right_index)?
            };
        }
        Ok(data_out)
    }

    /// Collect the handle intensities of `map` into one vector per map, in
    /// consensus-feature order and, inside a feature, in handle order (source
    /// `extractIntensityVectors`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a column header in `0..n` is
    /// missing (source `Exception::ElementNotFound`), when a handle names a
    /// map index outside `0..n` (undefined behaviour in the source), or when
    /// the handle count exceeds
    /// [`MAX_HANDLES`](super::consensus_map_normalizer_median::MAX_HANDLES).
    pub fn extract_intensity_vectors(map: &ConsensusMap) -> Result<Vec<Vec<f64>>> {
        let number_of_maps = check_maps(map)?;
        let mut out = vec![Vec::new(); number_of_maps];
        for feature in &map.features {
            for handle in feature.handles() {
                out[map_slot(handle)].push(f64::from(handle.intensity));
            }
        }
        Ok(out)
    }

    /// Write per-map intensity vectors back into the handles, in the order of
    /// [`Self::extract_intensity_vectors`] (source
    /// `setNormalizedIntensityValues`). Values are stored as `f32`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `feature_ints` does not have one
    /// vector per column header holding exactly that map's handle count (the
    /// source reads past the end), when a handle names a map index outside the
    /// column headers, or when a value is not finite. On any error the map is
    /// unchanged.
    pub fn set_normalized_intensity_values(
        feature_ints: &[Vec<f64>],
        map: &mut ConsensusMap,
    ) -> Result<()> {
        let number_of_maps = check_maps(map)?;
        let mut counts = vec![0usize; number_of_maps];
        for feature in &map.features {
            for handle in feature.handles() {
                counts[map_slot(handle)] += 1;
            }
        }
        if feature_ints.len() != number_of_maps
            || feature_ints
                .iter()
                .zip(&counts)
                .any(|(values, &count)| values.len() != count)
        {
            return Err(Error::InvalidValue(
                "normalised intensities do not match the consensus map's handles".into(),
            ));
        }
        let mut progress = vec![0usize; number_of_maps];
        rewrite_intensities(map, |handle| {
            let slot = map_slot(handle);
            let value = feature_ints[slot][progress[slot]];
            progress[slot] += 1;
            value
        })
    }
}

/// The order `operator<` on `double` gives finite values, which the handles'
/// validation guarantees: `-0.0` and `0.0` compare equal, unlike under
/// `f64::total_cmp`.
fn source_less(a: &f64, b: &f64) -> std::cmp::Ordering {
    a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
}
