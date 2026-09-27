// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Shift superimposer by pose clustering.
//!
//! Ports `ANALYSIS/MAPMATCHING/PoseClusteringShiftSuperimposer.h` and its
//! `.cpp`. See `docs/MAP_ALIGNMENT_SUPPORT.md`.

use super::superimposer::{
    BaseSuperimposer, Histogram, Superimposition, WINDOW_LENGTH_BASELINE, entry, float, int,
    similarity, string, tree, used_points,
};

use crate::analysis::transformations::TransformationDescription;
use crate::concept::progress_logger::ProgressReporter;
use crate::kernel::{ConsensusFeature, ConsensusMap};
use crate::math::source_sort::source_sort_by;
use crate::param::{DefaultParamHandler, Param, ParamValue};
use crate::{Error, Result};

/// Parameter handler name of the source class.
pub const POSE_CLUSTERING_SHIFT_SUPERIMPOSER_NAME: &str = "PoseClusteringShiftSuperimposer";

/// A superimposer that finds a retention-time shift by pose clustering
/// (source `PoseClusteringShiftSuperimposer`).
///
/// It works on two consensus maps and computes the shift that maps the
/// elements of the scene as near as possible onto those of the model. Every
/// pair of a model and a scene element of similar m/z votes for the shift
/// between their retention times, a "pose" of the scene relative to the
/// model; the votes are weighted by intensity similarity and collected in a
/// histogram, whose cluster gives the shift, hence the name pose clustering.
/// See also `PoseClusteringAffineSuperimposer`.
#[derive(Clone, Debug, PartialEq)]
pub struct PoseClusteringShiftSuperimposer {
    /// Maximum m/z deviation of corresponding elements in different maps
    /// (source `mz_pair_max_distance`, default 0.5, minimum 0). This applies
    /// to the pairs considered in hashing.
    pub mz_pair_max_distance: f64,
    /// Maximum number of elements considered in each map, selected by
    /// intensity (source `num_used_points`, default 2000, minimum -1). Reduces
    /// the running time and disregards weak signals; -1 uses all points.
    pub num_used_points: i64,
    /// Size of the histogram buckets of the shift (source
    /// `shift_bucket_size`, default 3.0, minimum 0). A good choice is about the
    /// time between consecutive MS scans.
    pub shift_bucket_size: f64,
    /// Maximal shift considered during histogramming, in both directions
    /// (source `max_shift`, default 1000, minimum 0, advanced).
    pub max_shift: f64,
    /// `[DEBUG]` base file name for dumps of the histogram (source
    /// `dump_buckets`, advanced). Not ported: a non-empty value is refused.
    pub dump_buckets: String,
    /// `[DEBUG]` base file name for dumps of the hashed pairs (source
    /// `dump_pairs`, advanced). Not ported: a non-empty value is refused.
    pub dump_pairs: String,
    /// Maximum number of histogram buckets (native limit).
    pub max_buckets: usize,
}

impl Default for PoseClusteringShiftSuperimposer {
    fn default() -> Self {
        Self {
            mz_pair_max_distance: 0.5,
            num_used_points: 2000,
            shift_bucket_size: 3.0,
            max_shift: 1000.0,
            dump_buckets: String::new(),
            dump_pairs: String::new(),
            max_buckets: 10_000_000,
        }
    }
}

impl PoseClusteringShiftSuperimposer {
    /// The source parameter defaults, as the constructor sets them:
    /// `mz_pair_max_distance`, `num_used_points`, `shift_bucket_size`,
    /// `max_shift`, `dump_buckets` and `dump_pairs`, with their descriptions,
    /// restrictions and `advanced` tags, in declaration order.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters of this superimposer in the layout of
    /// [`PoseClusteringShiftSuperimposer::defaults`] (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let mut mz = entry(
            "mz_pair_max_distance",
            ParamValue::Float(self.mz_pair_max_distance),
            "Maximum of m/z deviation of corresponding elements in different maps.  This condition applies to the pairs considered in hashing.",
            false,
        );
        mz.min_float = 0.0;
        let mut points = entry(
            "num_used_points",
            ParamValue::Integer(self.num_used_points),
            "Maximum number of elements considered in each map (selected by intensity).  Use this to reduce the running time and to disregard weak signals during alignment.  For using all points, set this to -1.",
            false,
        );
        points.min_int = -1;
        let mut bucket = entry(
            "shift_bucket_size",
            ParamValue::Float(self.shift_bucket_size),
            "The shift of the retention time interval is being hashed into buckets of this size during pose clustering.  A good choice for this would be about the time between consecutive MS scans.",
            false,
        );
        bucket.min_float = 0.0;
        let mut shift = entry(
            "max_shift",
            ParamValue::Float(self.max_shift),
            "Maximal shift which is considered during histogramming.  This applies for both directions.",
            true,
        );
        shift.min_float = 0.0;
        tree(vec![
            mz,
            points,
            bucket,
            shift,
            entry(
                "dump_buckets",
                ParamValue::String(self.dump_buckets.clone()),
                "[DEBUG] If non-empty, base filename where hash table buckets will be dumped to.  A serial number for each invocation will be appended automatically.",
                true,
            ),
            entry(
                "dump_pairs",
                ParamValue::String(self.dump_pairs.clone()),
                "[DEBUG] If non-empty, base filename where the individual hashed pairs will be dumped to (large!).  A serial number for each invocation will be appended automatically.",
                true,
            ),
        ])
    }

    /// Build a superimposer from parameters, as source `setParameters`.
    /// Missing parameters take their defaults; unknown ones are ignored.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong type or
    /// violates its restriction (the source's `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`PoseClusteringShiftSuperimposer::from_param`], also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`PoseClusteringShiftSuperimposer::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(POSE_CLUSTERING_SHIFT_SUPERIMPOSER_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, |p| {
            Ok(Self {
                mz_pair_max_distance: float(p, "mz_pair_max_distance")?,
                num_used_points: int(p, "num_used_points")?,
                shift_bucket_size: float(p, "shift_bucket_size")?,
                max_shift: float(p, "max_shift")?,
                dump_buckets: string(p, "dump_buckets")?,
                dump_pairs: string(p, "dump_pairs")?,
                ..Self::default()
            })
        })
    }

    fn validate(&self) -> Result<()> {
        if !self.dump_buckets.is_empty() || !self.dump_pairs.is_empty() {
            return Err(Error::Unsupported(
                "the [DEBUG] dump_buckets and dump_pairs outputs are not ported".into(),
            ));
        }
        if !(self.mz_pair_max_distance >= 0.0) || !self.mz_pair_max_distance.is_finite() {
            return Err(Error::InvalidValue(
                "mz_pair_max_distance must be finite and not negative".into(),
            ));
        }
        if !(self.shift_bucket_size > 0.0) || !self.shift_bucket_size.is_finite() {
            // The source divides by it and converts the quotient to Int,
            // which is undefined for a zero bucket size.
            return Err(Error::InvalidValue(
                "shift_bucket_size must be finite and positive".into(),
            ));
        }
        if !self.max_shift.is_finite() {
            return Err(Error::InvalidValue("max_shift must be finite".into()));
        }
        Ok(())
    }
}

/// The most intense `limit` elements of a map (all if it has no more), in
/// ascending m/z order: `ConstRefVector::sortByIntensity(true)`, `resize` and
/// `sortByComparator(Peak2D::MZLess())`, with the Release build's
/// `std::sort` permutations.
fn select(map: &ConsensusMap, limit: usize) -> Result<Vec<&ConsensusFeature>> {
    let mut elements: Vec<&ConsensusFeature> = map.features.iter().collect();
    if elements.len() > limit {
        source_sort_by(&mut elements, |a, b| b.intensity < a.intensity)?;
        elements.truncate(limit);
    }
    source_sort_by(&mut elements, |a, b| a.mz < b.mz)?;
    Ok(elements)
}

impl BaseSuperimposer for PoseClusteringShiftSuperimposer {
    /// Estimate the shift from `scene` onto `model` (source `run`).
    ///
    /// The source states as a precondition that both maps' ranges are up to
    /// date (`updateRanges`) and trusts them "for performance reasons"; the
    /// ranges are computed here from the maps, consensus handles included, as
    /// `ConsensusMap::updateRanges` does.
    ///
    /// The result is a `linear` transformation with slope 1 and the estimated
    /// shift as intercept, fitted from explicit parameters without data.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an empty map (the source then loops
    /// past the end of the model map), an invalid parameter value, a
    /// histogram of more than [`max_buckets`](Self::max_buckets), a
    /// non-finite total intensity ratio or pair similarity, a histogram whose
    /// noise level the source cannot determine without reading past its
    /// buffer, or non-finite statistics; [`Error::Unsupported`] for the debug
    /// dumps. Progress errors are returned as well.
    fn superimpose(
        &self,
        model: &ConsensusMap,
        scene: &ConsensusMap,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<Superimposition> {
        self.validate()?;
        if model.features.is_empty() || scene.features.is_empty() {
            return Err(Error::InvalidValue(
                "One of the input maps is empty! This is not allowed!".into(),
            ));
        }
        let model_range = model.ranges()?.rt;
        let scene_range = scene.ranges()?.rt;
        let (Some(model_rt), Some(scene_rt)) = (model_range, scene_range) else {
            return Err(Error::InvalidValue("input maps have no RT range".into()));
        };
        progress.section(0, 100, "shift pose clustering", |progress| {
            let mut actual = 0i64;
            actual += 1;
            progress.set(actual)?;
            actual += 1;
            progress.set(actual)?;
            actual += 1;
            progress.set(actual)?;
            let limit = used_points(self.num_used_points);
            let model_map = select(model, limit)?;
            actual += 1;
            progress.set(actual)?;
            let scene_map = select(scene, limit)?;
            actual += 1;
            progress.set(actual)?;
            actual = 10;
            progress.set(actual)?;

            let mut max_shift = self.max_shift.abs();
            let diff = (model_rt.max - scene_rt.min).abs();
            if max_shift > diff {
                max_shift = diff;
            }
            let diff = (model_rt.min - scene_rt.max).abs();
            if max_shift > diff {
                max_shift = diff;
            }
            let half = (max_shift / self.shift_bucket_size).ceil() + 4.0;
            if !(half <= (self.max_buckets / 2) as f64) {
                return Err(Error::InvalidValue(
                    "pose clustering histogram exceeds its bucket limit".into(),
                ));
            }
            let half = half as usize;
            let mut hash = Histogram::new(1 + 2 * half, self.shift_bucket_size, half as f64, 0.0);
            actual += 1;
            progress.set(actual)?;

            let mut total_model = 0.0;
            for element in &model_map {
                total_model += f64::from(element.intensity);
            }
            actual += 1;
            progress.set(actual)?;
            let mut total_scene = 0.0;
            for element in &scene_map {
                total_scene += f64::from(element.intensity);
            }
            actual += 1;
            progress.set(actual)?;
            let ratio = total_model / total_scene;
            if !ratio.is_finite() {
                return Err(Error::InvalidValue(
                    "the total intensity ratio of the maps is not finite".into(),
                ));
            }
            actual = 20;
            progress.set(actual)?;

            actual += 1;
            progress.set(actual)?;
            let distance = self.mz_pair_max_distance;
            let (n, m) = (model_map.len(), scene_map.len());
            let (mut i_low, mut i_high, mut k_low, mut k_high) = (0usize, 0usize, 0usize, 0usize);
            for i in 0..n - 1 {
                progress.set((actual as f32 + i as f32 / n as f32 * 10.0) as i64)?;
                let mz = model_map[i].mz;
                while i_low < n && model_map[i_low].mz < mz - distance {
                    i_low += 1;
                }
                while i_high < n && model_map[i_high].mz <= mz + distance {
                    i_high += 1;
                }
                let i_weight = 1.0 / (i_high - i_low) as f64 - WINDOW_LENGTH_BASELINE;
                if i_weight <= 0.0 {
                    continue;
                }
                while k_low < m && scene_map[k_low].mz < mz - distance {
                    k_low += 1;
                }
                while k_high < m && scene_map[k_high].mz <= mz + distance {
                    k_high += 1;
                }
                for k in k_low..k_high {
                    let k_weight = 1.0 / (k_high - k_low) as f64 - WINDOW_LENGTH_BASELINE;
                    if k_weight <= 0.0 {
                        continue;
                    }
                    let mut weight =
                        similarity(model_map[i].intensity, scene_map[k].intensity, ratio)?;
                    weight *= i_weight;
                    weight *= k_weight;
                    hash.add_value(model_map[i].rt - scene_map[k].rt, weight);
                }
            }
            actual = 30;
            progress.set(actual)?;

            actual += 1;
            progress.set(actual)?;
            hash.tophat();
            actual += 1;
            progress.set(actual)?;
            hash.cut_noise()?;
            // The source reports after computing the cut-off, after applying
            // it and after its optional dump.
            for _ in 0..3 {
                actual += 1;
                progress.set(actual)?;
            }
            let (intercept, _) = hash.centroid(0.0, 1.0)?;
            actual += 1;
            progress.set(actual)?;
            progress.set(80)?;
            actual += 1;
            progress.set(actual)?;

            let mut params = Param::new();
            params.set_value("slope", ParamValue::Float(1.0), "", &[])?;
            params.set_value("intercept", ParamValue::Float(intercept), "", &[])?;
            let mut transformation = TransformationDescription::default();
            transformation.fit_model_with_param("linear", &params)?;
            actual += 1;
            progress.set(actual)?;
            Ok(Superimposition {
                transformation,
                warnings: Vec::new(),
            })
        })
    }
}
