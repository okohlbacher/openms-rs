// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Removing and merging overlapping features
//! (`PROCESSING/FEATURE/FeatureOverlapFilter.h` and `.cpp`).
//!
//! The source class is a set of static functions; they become associated
//! functions of the unit struct `FeatureOverlapFilter`. The filter sorts the
//! features with a caller's "less" comparator, so the best feature comes first,
//! stores every feature in a region quadtree (`quadtree::Quadtree`), and then
//! visits the features in sorted order. Each feature that has not been removed
//! queries the tree for the features whose boxes intersect its own, and for each
//! candidate that passes the mode's overlap test it calls a callback. When the
//! callback returns `true` the candidate's unique ID is marked removed. At the
//! end every feature whose unique ID is marked is erased; the survivors keep the
//! sorted order.
//!
//! Three overlap modes pick the boxes and the test:
//!
//! - `ConvexHull`: the bounding box of the feature's convex hull, and no further
//!   test.
//! - `TraceLevel`: the same boxes, then an overlap of the mass-trace bounds that
//!   the subordinate features' hulls give.
//! - `CentroidBased`: a box of twice the tolerances around the centroid, then
//!   charge and FAIMS checks and an inclusive centroid-distance test.
//!
//! `FeatureOverlapFilter::merge_overlapping_features` and
//! `FeatureOverlapFilter::merge_faims_features` run the centroid mode with an
//! intensity comparator and a callback that sums intensities and records the
//! merged centroids.
//!
//! # Source conventions kept on purpose
//!
//! The port reproduces what the source computes, including its defects, which
//! are listed in the support document:
//!
//! - Boxes and the quadtree extent are `f32`, computed with the source's mix of
//!   `f32` and `f64` operations. Pairs whose float boxes do not intersect
//!   strictly are never tested, even when their exact distances are within the
//!   tolerances: a zero tolerance merges nothing, and a tolerance below the
//!   float spacing of the coordinates merges nothing.
//! - Removal is keyed by unique ID. Features that share ID 0 (unassigned) are
//!   all removed as soon as one of them is.
//! - The removed-set check applies only to the feature that queries. A removed
//!   feature can still be merged into a later survivor, so its intensity can be
//!   counted twice.
//! - `merge_faims_features` removes `FAIMS_CV` from a survivor after its first
//!   merge, so a survivor absorbs at most one feature.
//! - The quadtree extent follows `FeatureMap::updateRanges`, which skips a hull
//!   whose bounding box has zero width or height.
//! - A mass trace's retention-time end is read from the last point of its hull
//!   outline, which for most hulls is the first scan again, so the trace-level
//!   test usually compares only where the traces start.
//!
//! # Native differences
//!
//! - **Atomic.** An error leaves the map exactly as it was. The source sorts the
//!   map before it can throw, and `mergeFAIMSFeatures` has already moved the
//!   features out of the map when its merge throws, which strips their
//!   metadata.
//! - **Undefined source behaviour is refused.** A feature without a convex hull
//!   or with an empty one in the hull modes (the source converts its `±DBL_MAX`
//!   sentinel box to `float` and trips a Debug assertion, and in the trace mode
//!   reads the first point of an empty hull), a subordinate without a matching
//!   feature hull, a trace-level candidate without trace bounds (the source
//!   dereferences a missing map entry), a box or extent that does not fit
//!   `f32`, and a `FAIMS_CV` or merged-centroid list of the wrong type (the
//!   source throws or reads a union member of another type) all return errors.
//! - **Validated input.** Features must pass `Feature::validate`, nonzero unique
//!   IDs must be distinct (a `FeatureMap` invariant), and tolerances must be
//!   finite and nonnegative.
//! - **Bounded.** At most `FeatureOverlapFilter::MAX_FEATURES` features and
//!   `FeatureOverlapFilter::MAX_CANDIDATE_VISITS` candidates per call; a merge
//!   whose summed intensity leaves the `f32` range is an error.
//! - **Serial**, as the source.
//!
//! See `docs/FEATURE_OVERLAP_FILTER_SUPPORT.md` for the API mapping, the
//! preserved conventions, the native differences and the evidence.

/// The region quadtree the filter queries (`extern/Quadtree`).
pub mod quadtree;

use self::quadtree::{QuadBox, Quadtree};
use crate::concept::constants::user_param::FAIMS_CV;
use crate::kernel::{ConvexHull2D, Feature, FeatureMap};
use crate::metadata::MetaValue;
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashSet};

/// Meta value key of the merged centroid retention times, in seconds.
pub const MERGED_CENTROID_RTS: &str = "merged_centroid_rts";
/// Meta value key of the merged centroid m/z values.
pub const MERGED_CENTROID_MZS: &str = "merged_centroid_mzs";
/// Meta value key of the merged FAIMS compensation voltages.
pub const MERGED_CENTROID_IMS: &str = "merged_centroid_IMs";
/// Meta value key of the number of merged FAIMS compensation voltages.
pub const FAIMS_MERGE_COUNT: &str = "FAIMS_merge_count";

/// How the filter decides that two features overlap (source
/// `FeatureOverlapMode`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum FeatureOverlapMode {
    /// Intersecting convex-hull bounding boxes (source `CONVEX_HULL`, the
    /// documented default).
    #[default]
    ConvexHull,
    /// Intersecting hull bounding boxes and overlapping mass-trace bounds
    /// (source `TRACE_LEVEL`).
    ///
    /// Trace `i` of a feature takes its m/z bounds from the first and last
    /// outline points of the feature's hull `i` and its retention-time bounds
    /// from the first subordinate hull: the start is the first outline point
    /// with m/z above zero, and the end is found by walking the outline
    /// backwards, which for a hull whose first scan spans some m/z returns the
    /// start again. The port keeps this.
    TraceLevel,
    /// Centroids within the tolerances (source `CENTROID_BASED`).
    CentroidBased,
}

/// The source's backward-compatible alias `FeatureOverlapFilter::OverlapMode`.
pub type OverlapMode = FeatureOverlapMode;

/// How merged intensities combine (source `MergeIntensityMode`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum MergeIntensityMode {
    /// The sum of both intensities (source `SUM`, the default).
    #[default]
    Sum,
    /// The larger intensity (source `MAX`).
    Max,
}

/// Tolerances of the centroid mode (source `CentroidTolerances`).
///
/// Both tolerances must be finite and nonnegative; the source checks neither.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CentroidTolerances {
    /// Largest retention-time difference, in seconds (default 5.0). The test
    /// is inclusive.
    pub rt_tolerance: f64,
    /// Largest m/z difference, in Da (default 0.05). The test is inclusive.
    pub mz_tolerance: f64,
    /// Whether overlapping features must have equal charge (default `true`).
    pub require_same_charge: bool,
    /// Whether overlapping features must have equal `FAIMS_CV` values, or
    /// both lack one (default `false`).
    pub require_same_im: bool,
}

impl Default for CentroidTolerances {
    fn default() -> Self {
        Self {
            rt_tolerance: 5.0,
            mz_tolerance: 0.05,
            require_same_charge: true,
            require_same_im: false,
        }
    }
}

/// The callback `FeatureOverlapFilter::create_faims_merge_callback` builds
/// (source: the lambda `createFAIMSMergeCallback` returns).
///
/// Merging `other` into `best` combines their intensities according to
/// `intensity_mode` and, when `write_meta_values` is set, records on `best`:
///
/// - [`MERGED_CENTROID_RTS`]: `best`'s existing list, or `[best.rt]`, followed
///   by `other.rt`;
/// - [`MERGED_CENTROID_MZS`]: the same for m/z;
/// - [`MERGED_CENTROID_IMS`]: `best`'s existing list; otherwise `[best's
///   FAIMS_CV]`, and `FAIMS_CV` is removed from `best`; otherwise empty. Then
///   `other`'s `FAIMS_CV`, if any, is appended. The list and
///   [`FAIMS_MERGE_COUNT`] (its length) are written only when it is not empty.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct FaimsMergeCallback {
    /// How the intensities combine.
    pub intensity_mode: MergeIntensityMode,
    /// Whether the merge-tracking meta values are written.
    pub write_meta_values: bool,
}

impl FaimsMergeCallback {
    /// A callback with the given intensity mode and meta-value switch.
    pub const fn new(intensity_mode: MergeIntensityMode, write_meta_values: bool) -> Self {
        Self {
            intensity_mode,
            write_meta_values,
        }
    }

    /// Merge `other` into `best` and return `true`, which tells the filter to
    /// remove `other` (the source lambda always returns `true`).
    ///
    /// The combined intensity is computed in `f64` from both `f32` values and
    /// stored as `f32`, as the source's `double` arithmetic and `setIntensity`
    /// do.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`], leaving `best` unchanged, when an
    /// existing merged-centroid list of `best` is not a float list (the source
    /// throws `ConversionError`), when a `FAIMS_CV` that is read is not numeric
    /// (the source throws for an empty value and reads a union member of
    /// another type otherwise), or when the summed intensity is not a finite
    /// `f32` (the source stores infinity).
    pub fn merge(&self, best: &mut Feature, other: &Feature) -> Result<bool> {
        let plan = self.plan(best, other)?;
        plan.apply(best, None, 0)?;
        Ok(true)
    }

    fn plan(&self, best: &Feature, other: &Feature) -> Result<MergePlan> {
        let mut plan = MergePlan {
            intensity: combine(best.intensity, other.intensity, self.intensity_mode)?,
            ..MergePlan::default()
        };
        if self.write_meta_values {
            plan.rts = Some(extend_list(best, MERGED_CENTROID_RTS, best.rt, other.rt)?);
            plan.mzs = Some(extend_list(best, MERGED_CENTROID_MZS, best.mz, other.mz)?);
            let mut ims = Vec::new();
            if let Some(existing) = best.metadata.get(MERGED_CENTROID_IMS) {
                ims = float_list(existing, MERGED_CENTROID_IMS)?.to_vec();
            } else if best.metadata.contains_key(FAIMS_CV) {
                ims.push(faims_cv(best)?);
                plan.remove_cv = true;
            }
            if other.metadata.contains_key(FAIMS_CV) {
                ims.push(faims_cv(other)?);
            }
            if !ims.is_empty() {
                plan.ims = Some(ims);
            }
        }
        Ok(plan)
    }
}

/// Static functions of the source class `FeatureOverlapFilter`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct FeatureOverlapFilter;

impl FeatureOverlapFilter {
    /// Features one call accepts ([`FeatureMap::MAX_ITEMS`]). The source is
    /// unbounded.
    pub const MAX_FEATURES: usize = FeatureMap::MAX_ITEMS;

    /// Candidates, the querying feature itself included, one call examines
    /// before it fails. The source is unbounded; its worst case is quadratic in
    /// the number of features.
    pub const MAX_CANDIDATE_VISITS: u64 = 1 << 32;

    /// The source's default comparator: the higher overall quality is better.
    pub fn higher_overall_quality(left: &Feature, right: &Feature) -> bool {
        left.quality > right.quality
    }

    /// The source's default callback: every overlap is removed.
    pub fn always_overlapping(_best: &mut Feature, _other: &mut Feature) -> bool {
        true
    }

    /// Keep the best feature of each cluster of overlapping features (source
    /// `filter` with a `bool` switch).
    ///
    /// `check_overlap_at_trace_level` selects
    /// [`FeatureOverlapMode::TraceLevel`] when set (the source default) and
    /// [`FeatureOverlapMode::ConvexHull`] otherwise; the tolerances are the
    /// defaults and unused. The source's defaults for the other arguments are
    /// [`Self::higher_overall_quality`] and [`Self::always_overlapping`].
    ///
    /// # Errors
    ///
    /// As [`Self::filter_with_mode`].
    pub fn filter<L, O>(
        feature_map: &mut FeatureMap,
        comparator: L,
        on_overlap: O,
        check_overlap_at_trace_level: bool,
    ) -> Result<()>
    where
        L: FnMut(&Feature, &Feature) -> bool,
        O: FnMut(&mut Feature, &mut Feature) -> bool,
    {
        let mode = if check_overlap_at_trace_level {
            FeatureOverlapMode::TraceLevel
        } else {
            FeatureOverlapMode::ConvexHull
        };
        Self::filter_with_mode(
            feature_map,
            comparator,
            on_overlap,
            mode,
            &CentroidTolerances::default(),
        )
    }

    /// Keep the best feature of each cluster of overlapping features, with a
    /// chosen overlap mode (source `filter` with a `FeatureOverlapMode`).
    ///
    /// # Arguments
    ///
    /// * `comparator` — a "less" comparator, a strict weak ordering: when
    ///   several features overlap, the one that sorts first is kept. The
    ///   features are sorted stably with it. It is called once or twice per
    ///   comparison, where the source's merge sort calls it once. As with
    ///   `std::stable_sort`, a comparator that is not a strict weak ordering
    ///   gives an unspecified order, and the standard sort may panic on it.
    /// * `on_overlap` — called as `on_overlap(best, other)` when `other`
    ///   overlaps `best`, a feature that has not been removed. It may transfer
    ///   information to `best`. Returning `false` keeps `other` (as far as this
    ///   overlap is concerned); returning `true` marks `other`'s unique ID
    ///   removed. Boxes are recomputed from the features at every query, so a
    ///   callback that moves a feature changes later queries as it does in the
    ///   source.
    /// * `mode` — the overlap test.
    /// * `tolerances` — used only by [`FeatureOverlapMode::CentroidBased`].
    ///
    /// When the loop could fail after a callback has run — in the trace mode,
    /// in the centroid mode with `require_same_im`, and when more than
    /// [`Self::MAX_CANDIDATE_VISITS`] candidates are possible — the features
    /// are copied first so that an error can restore them.
    ///
    /// # Errors
    ///
    /// Returns, leaving the map unchanged:
    /// - [`Error::InvalidRange`] for an empty map (the source's `getMinMZ`
    ///   throws `InvalidRange`);
    /// - [`Error::MissingInformation`] in the hull modes for a feature without a
    ///   convex hull or with an empty one, which in the trace mode includes the
    ///   hull a subordinate's m/z bounds are read from (the source reads the
    ///   first point of an empty hull), and in the trace mode for a subordinate
    ///   without convex hulls (the source throws `MissingInformation`, after
    ///   sorting);
    /// - [`Error::InvalidValue`] for more than [`Self::MAX_FEATURES`] features,
    ///   a feature that fails [`Feature::validate`], a repeated nonzero unique
    ///   ID, a tolerance that is not finite and nonnegative, a box or quadtree
    ///   extent that does not fit `f32`, a trace-mode feature with more
    ///   subordinates than hulls, a trace-mode candidate without trace bounds,
    ///   a non-numeric `FAIMS_CV` read by the `require_same_im` check, or more
    ///   than [`Self::MAX_CANDIDATE_VISITS`] candidates.
    pub fn filter_with_mode<L, O>(
        feature_map: &mut FeatureMap,
        mut comparator: L,
        mut on_overlap: O,
        mode: FeatureOverlapMode,
        tolerances: &CentroidTolerances,
    ) -> Result<()>
    where
        L: FnMut(&Feature, &Feature) -> bool,
        O: FnMut(&mut Feature, &mut Feature) -> bool,
    {
        run(
            &mut feature_map.features,
            &mut comparator,
            mode,
            tolerances,
            &mut |features: &mut [Feature], best, other, _: &mut Journal| {
                let (best, other) = pair_mut(features, best, other);
                Ok(on_overlap(best, other))
            },
            Rollback::SnapshotWhenFallible,
        )
    }

    /// As [`Self::filter_with_mode`], with a callback that can fail, such as
    /// the one [`Self::create_faims_merge_callback`] returns.
    ///
    /// The source's callbacks report failure by throwing. Here a callback
    /// error stops the filter and is returned; the features are always copied
    /// first so that the map can be restored.
    ///
    /// # Errors
    ///
    /// As [`Self::filter_with_mode`], and any error of `on_overlap`.
    pub fn filter_with_fallible_callback<L, O>(
        feature_map: &mut FeatureMap,
        mut comparator: L,
        mut on_overlap: O,
        mode: FeatureOverlapMode,
        tolerances: &CentroidTolerances,
    ) -> Result<()>
    where
        L: FnMut(&Feature, &Feature) -> bool,
        O: FnMut(&mut Feature, &mut Feature) -> Result<bool>,
    {
        run(
            &mut feature_map.features,
            &mut comparator,
            mode,
            tolerances,
            &mut |features: &mut [Feature], best, other, _: &mut Journal| {
                let (best, other) = pair_mut(features, best, other);
                on_overlap(best, other)
            },
            Rollback::Snapshot,
        )
    }

    /// A callback that merges overlapping features, for
    /// [`Self::filter_with_fallible_callback`] (source
    /// `createFAIMSMergeCallback`, defaults `SUM` and `true`).
    ///
    /// See [`FaimsMergeCallback`] for what a merge records. Features without
    /// `FAIMS_CV` contribute nothing to [`MERGED_CENTROID_IMS`].
    pub fn create_faims_merge_callback(
        intensity_mode: MergeIntensityMode,
        write_meta_values: bool,
    ) -> impl FnMut(&mut Feature, &mut Feature) -> Result<bool> {
        let callback = FaimsMergeCallback::new(intensity_mode, write_meta_values);
        move |best: &mut Feature, other: &mut Feature| callback.merge(best, other)
    }

    /// Merge features whose centroids lie within the tolerances (source
    /// `mergeOverlappingFeatures`; source defaults 5.0 s, 0.05 Da, `true`,
    /// `false`, `SUM`, `true`).
    ///
    /// Runs the centroid mode with the higher intensity as the better feature
    /// and a [`FaimsMergeCallback`]. The survivors stay in descending intensity
    /// order, ties in input order.
    ///
    /// # Arguments
    ///
    /// * `max_rt_diff` — largest retention-time difference, in seconds.
    /// * `max_mz_diff` — largest m/z difference, in Da.
    /// * `require_same_charge` — merge only features of equal charge.
    /// * `require_same_im` — merge only features with equal `FAIMS_CV`, or both
    ///   without one.
    /// * `intensity_mode` — sum or maximum.
    /// * `write_meta_values` — record the merge-tracking meta values.
    ///
    /// # Errors
    ///
    /// As [`Self::filter_with_mode`] and [`FaimsMergeCallback::merge`]; the map
    /// is unchanged after any error.
    pub fn merge_overlapping_features(
        feature_map: &mut FeatureMap,
        max_rt_diff: f64,
        max_mz_diff: f64,
        require_same_charge: bool,
        require_same_im: bool,
        intensity_mode: MergeIntensityMode,
        write_meta_values: bool,
    ) -> Result<()> {
        let tolerances = CentroidTolerances {
            rt_tolerance: max_rt_diff,
            mz_tolerance: max_mz_diff,
            require_same_charge,
            require_same_im,
        };
        let callback = FaimsMergeCallback::new(intensity_mode, write_meta_values);
        run(
            &mut feature_map.features,
            &mut higher_intensity,
            FeatureOverlapMode::CentroidBased,
            &tolerances,
            &mut journaled_merge(callback),
            Rollback::Journal,
        )
    }

    /// Merge features that are the same analyte seen at different FAIMS
    /// compensation voltages (source `mergeFAIMSFeatures`; source defaults
    /// 5.0 s and 0.05 Da).
    ///
    /// Only features with a `FAIMS_CV` meta value take part. When none has one,
    /// the map is not touched at all. Otherwise the features with `FAIMS_CV` are
    /// merged when there are at least two of them, in the centroid mode with
    /// equal charge required and the higher intensity as the better feature,
    /// by a callback that merges only when both features still carry
    /// `FAIMS_CV` and the values differ exactly. It sums the intensities,
    /// extends the merged-centroid lists and removes `FAIMS_CV` from the
    /// survivor on its first merge. The map then holds the FAIMS features in
    /// merged (intensity) order followed by the other features in input
    /// order. Map-level data (identifier, unique ID, meta values,
    /// identifications, processing) is kept, as the source's `clear(false)`
    /// keeps it.
    ///
    /// # Arguments
    ///
    /// * `max_rt_diff` — largest retention-time difference, in seconds.
    /// * `max_mz_diff` — largest m/z difference, in Da.
    ///
    /// # Errors
    ///
    /// As [`Self::merge_overlapping_features`], for the FAIMS features only;
    /// the map, including its feature order, is unchanged after any error.
    pub fn merge_faims_features(
        feature_map: &mut FeatureMap,
        max_rt_diff: f64,
        max_mz_diff: f64,
    ) -> Result<()> {
        if !feature_map
            .features
            .iter()
            .any(|feature| feature.metadata.contains_key(FAIMS_CV))
        {
            return Ok(());
        }
        let input = std::mem::take(&mut feature_map.features);
        let is_faims: Vec<bool> = input
            .iter()
            .map(|feature| feature.metadata.contains_key(FAIMS_CV))
            .collect();
        let faims_count = is_faims.iter().filter(|&&faims| faims).count();
        let mut faims = Vec::with_capacity(faims_count);
        let mut others = Vec::with_capacity(input.len() - faims_count);
        for (feature, &faims_feature) in input.into_iter().zip(&is_faims) {
            if faims_feature {
                faims.push(feature);
            } else {
                others.push(feature);
            }
        }
        if faims.len() > 1 {
            let tolerances = CentroidTolerances {
                rt_tolerance: max_rt_diff,
                mz_tolerance: max_mz_diff,
                require_same_charge: true,
                require_same_im: false,
            };
            let merged = run(
                &mut faims,
                &mut higher_intensity,
                FeatureOverlapMode::CentroidBased,
                &tolerances,
                &mut merge_different_voltages,
                Rollback::Journal,
            );
            if let Err(error) = merged {
                feature_map.features = interleave(faims, others, &is_faims);
                return Err(error);
            }
        }
        faims.append(&mut others);
        feature_map.features = faims;
        Ok(())
    }
}

/// The comparator of the merge functions: the higher intensity is better.
fn higher_intensity(left: &Feature, right: &Feature) -> bool {
    left.intensity > right.intensity
}

/// The callback of `mergeOverlappingFeatures`: `callback` applied through the
/// journal.
fn journaled_merge(
    callback: FaimsMergeCallback,
) -> impl FnMut(&mut [Feature], usize, usize, &mut Journal) -> Result<bool> {
    move |features: &mut [Feature], best: usize, other: usize, journal: &mut Journal| {
        let plan = callback.plan(&features[best], &features[other])?;
        plan.apply(&mut features[best], Some(journal), best)?;
        Ok(true)
    }
}

/// The callback of `mergeFAIMSFeatures` (FeatureOverlapFilter.cpp:430-499).
fn merge_different_voltages(
    features: &mut [Feature],
    best_index: usize,
    other_index: usize,
    journal: &mut Journal,
) -> Result<bool> {
    let best = &features[best_index];
    let other = &features[other_index];
    if !best.metadata.contains_key(FAIMS_CV) || !other.metadata.contains_key(FAIMS_CV) {
        return Ok(false);
    }
    let best_cv = faims_cv(best)?;
    let other_cv = faims_cv(other)?;
    if best_cv == other_cv {
        return Ok(false);
    }
    let mut plan = MergePlan {
        intensity: combine(best.intensity, other.intensity, MergeIntensityMode::Sum)?,
        rts: Some(extend_list(best, MERGED_CENTROID_RTS, best.rt, other.rt)?),
        mzs: Some(extend_list(best, MERGED_CENTROID_MZS, best.mz, other.mz)?),
        ..MergePlan::default()
    };
    let mut ims = match best.metadata.get(MERGED_CENTROID_IMS) {
        Some(existing) => float_list(existing, MERGED_CENTROID_IMS)?.to_vec(),
        None => {
            plan.remove_cv = true;
            vec![best_cv]
        }
    };
    ims.push(other_cv);
    plan.ims = Some(ims);
    plan.apply(&mut features[best_index], Some(journal), best_index)?;
    Ok(true)
}

/// The changes one merge makes to the surviving feature, computed before any
/// of them is applied.
#[derive(Default)]
struct MergePlan {
    intensity: f32,
    rts: Option<Vec<f64>>,
    mzs: Option<Vec<f64>>,
    ims: Option<Vec<f64>>,
    remove_cv: bool,
}

impl MergePlan {
    fn apply(
        self,
        best: &mut Feature,
        mut journal: Option<&mut Journal>,
        index: usize,
    ) -> Result<()> {
        let rts = self.rts.map(MetaValue::try_from).transpose()?;
        let mzs = self.mzs.map(MetaValue::try_from).transpose()?;
        let ims = match self.ims {
            Some(ims) => {
                let count = i64::try_from(ims.len()).map_err(|_| {
                    Error::InvalidValue("merged FAIMS voltage count overflows".into())
                })?;
                Some((MetaValue::try_from(ims)?, MetaValue::from(count)))
            }
            None => None,
        };
        // Only the value a slot held before the run is journaled: a later change
        // of the same slot drops the superseded value, as the source does, so the
        // journal stays linear in the number of features.
        let mut set = |slot: MetaSlot, value: Option<MetaValue>, best: &mut Feature| {
            let key = slot.key();
            let old = match value {
                Some(value) => best.metadata.insert(key.to_owned(), value),
                None => best.metadata.remove(key),
            };
            if let Some(journal) = journal.as_deref_mut() {
                journal.record_meta(index, slot, old);
            }
        };
        if let Some(rts) = rts {
            set(MetaSlot::MergedRts, Some(rts), best);
        }
        if let Some(mzs) = mzs {
            set(MetaSlot::MergedMzs, Some(mzs), best);
        }
        if self.remove_cv {
            set(MetaSlot::FaimsCv, None, best);
        }
        if let Some((ims, count)) = ims {
            set(MetaSlot::MergedIms, Some(ims), best);
            set(MetaSlot::MergeCount, Some(count), best);
        }
        if let Some(journal) = journal {
            journal.record_intensity(index, best.intensity);
        }
        best.intensity = self.intensity;
        Ok(())
    }
}

/// `best`'s existing list under `key`, or `[best_value]`, with `other_value`
/// appended.
fn extend_list(best: &Feature, key: &str, best_value: f64, other_value: f64) -> Result<Vec<f64>> {
    let mut list = match best.metadata.get(key) {
        Some(existing) => float_list(existing, key)?.to_vec(),
        None => vec![best_value],
    };
    list.push(other_value);
    Ok(list)
}

fn float_list<'a>(value: &'a MetaValue, key: &str) -> Result<&'a [f64]> {
    value.as_float_list().map_err(|_| {
        Error::InvalidValue(format!(
            "meta value {key} is not a float list (the source throws ConversionError)"
        ))
    })
}

/// The feature's `FAIMS_CV` as a number: a float, or an integer converted as
/// the source's `DataValue::operator double` converts it.
fn faims_cv(feature: &Feature) -> Result<f64> {
    match feature.metadata.get(FAIMS_CV) {
        Some(value) => value.as_f64().map_err(|_| {
            Error::InvalidValue(format!(
                "{FAIMS_CV} of the feature with unique ID {} is not numeric",
                feature.unique_id
            ))
        }),
        None => Err(Error::InvalidValue(format!(
            "the feature with unique ID {} has no {FAIMS_CV}",
            feature.unique_id
        ))),
    }
}

fn combine(best: f32, other: f32, mode: MergeIntensityMode) -> Result<f32> {
    let (best, other) = (f64::from(best), f64::from(other));
    let combined = match mode {
        MergeIntensityMode::Sum => best + other,
        // std::max(best, other): the second only when the first is less.
        MergeIntensityMode::Max => {
            if best < other {
                other
            } else {
                best
            }
        }
    };
    let stored = combined as f32;
    if stored.is_finite() {
        Ok(stored)
    } else {
        Err(Error::InvalidValue(
            "the merged intensity does not fit f32".into(),
        ))
    }
}

/// Restores `faims` and `others` to the input order recorded in `is_faims`.
fn interleave(faims: Vec<Feature>, others: Vec<Feature>, is_faims: &[bool]) -> Vec<Feature> {
    let mut faims = faims.into_iter();
    let mut others = others.into_iter();
    let mut features = Vec::with_capacity(is_faims.len());
    for &faims_feature in is_faims {
        let next = if faims_feature {
            faims.next()
        } else {
            others.next()
        };
        features.extend(next);
    }
    features
}

/// Two distinct features of a slice, mutably.
fn pair_mut(features: &mut [Feature], first: usize, second: usize) -> (&mut Feature, &mut Feature) {
    if first < second {
        let (head, tail) = features.split_at_mut(second);
        (&mut head[first], &mut tail[0])
    } else {
        let (head, tail) = features.split_at_mut(first);
        (&mut tail[0], &mut head[second])
    }
}

// ---------------------------------------------------------------------------
// The filter itself

/// How a failed run restores the features.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Rollback {
    /// The callback writes only through the journal.
    Journal,
    /// Copy the features first when a failure after a callback is possible.
    SnapshotWhenFallible,
    /// Always copy the features first.
    Snapshot,
}

/// A meta value that the built-in merge callbacks change.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum MetaSlot {
    MergedRts,
    MergedMzs,
    FaimsCv,
    MergedIms,
    MergeCount,
}

impl MetaSlot {
    /// The slot's meta value key.
    fn key(self) -> &'static str {
        match self {
            MetaSlot::MergedRts => MERGED_CENTROID_RTS,
            MetaSlot::MergedMzs => MERGED_CENTROID_MZS,
            MetaSlot::FaimsCv => FAIMS_CV,
            MetaSlot::MergedIms => MERGED_CENTROID_IMS,
            MetaSlot::MergeCount => FAIMS_MERGE_COUNT,
        }
    }

    /// The slot's bit in [`Journal::touched`].
    fn bit(self) -> u8 {
        1 << (self as u8)
    }
}

/// The bit of the intensity in [`Journal::touched`], after the [`MetaSlot`]s.
const INTENSITY_BIT: u8 = 1 << 5;

enum Undo {
    Intensity {
        index: usize,
        old: f32,
    },
    Meta {
        index: usize,
        slot: MetaSlot,
        old: Option<MetaValue>,
    },
}

/// The values the built-in merge callbacks overwrote, for an exact rollback.
///
/// Each field of each feature is recorded once, with the value it held before
/// the run, so the journal holds at most six entries per feature however many
/// merges a survivor absorbs. Undoing the entries in reverse restores the
/// features exactly.
#[derive(Default)]
struct Journal {
    /// The first change of each field, in the order the changes happened.
    entries: Vec<Undo>,
    /// Per feature index, one bit per [`MetaSlot`] and [`INTENSITY_BIT`] for
    /// the fields already recorded.
    touched: Vec<u8>,
}

impl Journal {
    /// Whether the field `bit` of the feature at `index` changes for the first
    /// time in this run; marks it changed.
    fn first_change(&mut self, index: usize, bit: u8) -> bool {
        if self.touched.len() <= index {
            self.touched.resize(index + 1, 0);
        }
        let first = self.touched[index] & bit == 0;
        self.touched[index] |= bit;
        first
    }

    /// Records `old`, the value `slot` of the feature at `index` held before this
    /// change, if it is the slot's first change.
    fn record_meta(&mut self, index: usize, slot: MetaSlot, old: Option<MetaValue>) {
        if self.first_change(index, slot.bit()) {
            self.entries.push(Undo::Meta { index, slot, old });
        }
    }

    /// Records `old`, the intensity of the feature at `index` before this
    /// change, if it is the intensity's first change.
    fn record_intensity(&mut self, index: usize, old: f32) {
        if self.first_change(index, INTENSITY_BIT) {
            self.entries.push(Undo::Intensity { index, old });
        }
    }

    /// Undoes the recorded changes in reverse order.
    fn rollback(self, features: &mut [Feature]) {
        for entry in self.entries.into_iter().rev() {
            match entry {
                Undo::Intensity { index, old } => {
                    if let Some(feature) = features.get_mut(index) {
                        feature.intensity = old;
                    }
                }
                Undo::Meta { index, slot, old } => {
                    if let Some(feature) = features.get_mut(index) {
                        match old {
                            Some(value) => {
                                feature.metadata.insert(slot.key().to_owned(), value);
                            }
                            None => {
                                feature.metadata.remove(slot.key());
                            }
                        }
                    }
                }
            }
        }
    }
}

fn run<L, O>(
    features: &mut Vec<Feature>,
    comparator: &mut L,
    mode: FeatureOverlapMode,
    tolerances: &CentroidTolerances,
    on_overlap: &mut O,
    rollback: Rollback,
) -> Result<()>
where
    L: FnMut(&Feature, &Feature) -> bool,
    O: FnMut(&mut [Feature], usize, usize, &mut Journal) -> Result<bool>,
{
    let plan = Plan::new(features, mode, tolerances)?;
    let count = features.len() as u64;
    let snapshot = match rollback {
        Rollback::Journal => None,
        Rollback::Snapshot => Some(features.clone()),
        Rollback::SnapshotWhenFallible => {
            let fallible = mode == FeatureOverlapMode::TraceLevel
                || (mode == FeatureOverlapMode::CentroidBased && tolerances.require_same_im)
                || count.saturating_mul(count) > FeatureOverlapFilter::MAX_CANDIDATE_VISITS;
            fallible.then(|| features.clone())
        }
    };

    // std::stable_sort with the "less" comparator, as a permutation so that a
    // failed run can put every feature back.
    let mut order: Vec<usize> = (0..features.len()).collect();
    order.sort_by(|&a, &b| {
        if comparator(&features[a], &features[b]) {
            Ordering::Less
        } else if comparator(&features[b], &features[a]) {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    });
    let mut sorted: Vec<Feature> = order
        .iter()
        .map(|&index| std::mem::take(&mut features[index]))
        .collect();

    let mut journal = Journal::default();
    match overlap_loop(&mut sorted, &plan, on_overlap, &mut journal) {
        Ok(removed) => {
            sorted.retain(|feature| !removed.contains(&feature.unique_id));
            *features = sorted;
            Ok(())
        }
        Err(error) => {
            match snapshot {
                Some(snapshot) => *features = snapshot,
                None => {
                    journal.rollback(&mut sorted);
                    for (feature, &index) in sorted.into_iter().zip(&order) {
                        features[index] = feature;
                    }
                }
            }
            Err(error)
        }
    }
}

fn overlap_loop<O>(
    features: &mut [Feature],
    plan: &Plan,
    on_overlap: &mut O,
    journal: &mut Journal,
) -> Result<HashSet<u64>>
where
    O: FnMut(&mut [Feature], usize, usize, &mut Journal) -> Result<bool>,
{
    let mut tree = Quadtree::new(plan.extent);
    {
        let current: &[Feature] = features;
        let get_box = |&index: &usize| plan.feature_box(&current[index]);
        for index in 0..current.len() {
            tree.add(index, &get_box)?;
        }
    }
    let mut removed = HashSet::new();
    let mut candidates = Vec::new();
    let mut visits = 0u64;
    for index in 0..features.len() {
        if removed.contains(&features[index].unique_id) {
            continue;
        }
        {
            let current: &[Feature] = features;
            let get_box = |&candidate: &usize| plan.feature_box(&current[candidate]);
            tree.query_into(plan.feature_box(&current[index]), &get_box, &mut candidates);
        }
        visits = visits.saturating_add(candidates.len() as u64);
        if visits > FeatureOverlapFilter::MAX_CANDIDATE_VISITS {
            return Err(Error::InvalidValue(format!(
                "feature overlap filtering examines more than {} candidates",
                FeatureOverlapFilter::MAX_CANDIDATE_VISITS
            )));
        }
        for &candidate in &candidates {
            if candidate == index {
                continue;
            }
            if plan.overlaps(&features[index], &features[candidate])?
                && on_overlap(features, index, candidate, journal)?
            {
                removed.insert(features[candidate].unique_id);
            }
        }
    }
    Ok(removed)
}

/// A retention-time and m/z box in `f64` with the source `DBoundingBox<2>`
/// semantics: an empty box is `min = +DBL_MAX, max = -DBL_MAX`.
#[derive(Clone, Copy)]
struct SourceBox {
    min: [f64; 2],
    max: [f64; 2],
}

impl SourceBox {
    const EMPTY: Self = Self {
        min: [f64::MAX; 2],
        max: [f64::MIN; 2],
    };

    fn of_hull(hull: &ConvexHull2D) -> Self {
        match hull.bounding_box() {
            Some(bounds) => Self {
                min: [bounds.rt_range().min, bounds.mz_range().min],
                max: [bounds.rt_range().max, bounds.mz_range().max],
            },
            None => Self::EMPTY,
        }
    }

    /// `Feature::getConvexHull().getBoundingBox()`: the only hull's box, or the
    /// box `DBoundingBox::enlarge` builds from every hull's corners.
    fn of_feature(feature: &Feature) -> Self {
        match feature.convex_hulls.as_slice() {
            [] => Self::EMPTY,
            [hull] => Self::of_hull(hull),
            hulls => {
                let mut bounds = Self::EMPTY;
                for hull in hulls {
                    let hull_bounds = Self::of_hull(hull);
                    bounds.enlarge(hull_bounds.min);
                    bounds.enlarge(hull_bounds.max);
                }
                bounds
            }
        }
    }

    fn enlarge(&mut self, point: [f64; 2]) {
        for ((&value, min), max) in point.iter().zip(&mut self.min).zip(&mut self.max) {
            if value < *min {
                *min = value;
            }
            if value > *max {
                *max = value;
            }
        }
    }

    /// `DBoundingBox::isEmpty`: true when any dimension has `max <= min`, so a
    /// box of zero width or height counts as empty.
    fn is_empty(&self) -> bool {
        self.max.iter().zip(&self.min).any(|(max, min)| max <= min)
    }
}

/// `RangeBase::extend`, `std::min` and `std::max` on finite values.
struct Range {
    min: f64,
    max: f64,
}

impl Range {
    fn new(value: f64) -> Self {
        Self {
            min: value,
            max: value,
        }
    }

    fn extend(&mut self, value: f64) {
        if value < self.min {
            self.min = value;
        }
        if self.max < value {
            self.max = value;
        }
    }
}

#[derive(Clone, Copy)]
struct TraceBounds {
    rt_min: f64,
    rt_max: f64,
    mz_min: f64,
    mz_max: f64,
}

/// Everything the loop needs, checked before the features are touched.
struct Plan {
    mode: FeatureOverlapMode,
    tolerances: CentroidTolerances,
    extent: QuadBox,
    mz_box_width: f32,
    rt_box_height: f32,
    traces: BTreeMap<u64, Vec<TraceBounds>>,
}

impl Plan {
    fn new(
        features: &[Feature],
        mode: FeatureOverlapMode,
        tolerances: &CentroidTolerances,
    ) -> Result<Self> {
        if features.len() > FeatureOverlapFilter::MAX_FEATURES {
            return Err(Error::InvalidValue(format!(
                "feature overlap filtering accepts at most {} features",
                FeatureOverlapFilter::MAX_FEATURES
            )));
        }
        let centroid = mode == FeatureOverlapMode::CentroidBased;
        if centroid {
            for (name, value) in [
                ("RT", tolerances.rt_tolerance),
                ("m/z", tolerances.mz_tolerance),
            ] {
                if !value.is_finite() || value < 0.0 {
                    return Err(Error::InvalidValue(format!(
                        "the {name} tolerance must be finite and nonnegative"
                    )));
                }
            }
        }
        let mut assigned = Vec::with_capacity(features.len());
        for feature in features {
            feature.validate()?;
            if feature.unique_id != 0 {
                assigned.push(feature.unique_id);
            }
        }
        assigned.sort_unstable();
        if let Some(pair) = assigned.windows(2).find(|pair| pair[0] == pair[1]) {
            return Err(Error::InvalidValue(format!(
                "duplicate assigned feature ID {}",
                pair[0]
            )));
        }
        if features.is_empty() {
            return Err(Error::InvalidRange(
                "cannot filter an empty feature map (the source's range query throws InvalidRange)"
                    .into(),
            ));
        }
        if !centroid {
            for feature in features {
                if feature.convex_hulls.is_empty()
                    || feature.convex_hulls.iter().any(ConvexHull2D::is_empty)
                {
                    return Err(Error::MissingInformation(format!(
                        "the feature with unique ID {} needs non-empty convex hulls for hull-based overlap detection",
                        feature.unique_id
                    )));
                }
            }
        }
        let traces = if mode == FeatureOverlapMode::TraceLevel {
            trace_bounds(features)?
        } else {
            BTreeMap::new()
        };

        // FeatureMap::updateRanges: centroids, then non-empty hull boxes.
        let mut rt = Range::new(features[0].rt);
        let mut mz = Range::new(features[0].mz);
        for feature in features {
            rt.extend(feature.rt);
            mz.extend(feature.mz);
        }
        for feature in features {
            let bounds = SourceBox::of_feature(feature);
            if !bounds.is_empty() {
                rt.extend(bounds.min[0]);
                rt.extend(bounds.max[0]);
                mz.extend(bounds.min[1]);
                mz.extend(bounds.max[1]);
            }
        }
        let mut min_mz = mz.min as f32;
        let mut max_mz = mz.max as f32;
        let mut min_rt = rt.min as f32;
        let mut max_rt = rt.max as f32;
        if centroid {
            min_mz = (f64::from(min_mz) - tolerances.mz_tolerance) as f32;
            max_mz = (f64::from(max_mz) + tolerances.mz_tolerance) as f32;
            min_rt = (f64::from(min_rt) - tolerances.rt_tolerance) as f32;
            max_rt = (f64::from(max_rt) + tolerances.rt_tolerance) as f32;
        }
        let extent = QuadBox::new(
            min_mz - 1.0,
            min_rt - 1.0,
            max_mz - min_mz + 2.0,
            max_rt - min_rt + 2.0,
        );
        let plan = Self {
            mode,
            tolerances: *tolerances,
            extent,
            mz_box_width: (2.0 * tolerances.mz_tolerance) as f32,
            rt_box_height: (2.0 * tolerances.rt_tolerance) as f32,
            traces,
        };
        if !finite_box(extent) {
            return Err(Error::InvalidValue(
                "the feature map's extent does not fit f32".into(),
            ));
        }
        for feature in features {
            if !finite_box(plan.feature_box(feature)) {
                return Err(Error::InvalidValue(format!(
                    "the overlap box of the feature with unique ID {} does not fit f32",
                    feature.unique_id
                )));
            }
        }
        Ok(plan)
    }

    /// The quadtree box of a feature (the source's `getBox` lambdas): x is
    /// m/z, y is retention time.
    fn feature_box(&self, feature: &Feature) -> QuadBox {
        if self.mode == FeatureOverlapMode::CentroidBased {
            let rt = feature.rt as f32;
            let mz = feature.mz as f32;
            QuadBox::new(
                (f64::from(mz) - self.tolerances.mz_tolerance) as f32,
                (f64::from(rt) - self.tolerances.rt_tolerance) as f32,
                self.mz_box_width,
                self.rt_box_height,
            )
        } else {
            let bounds = SourceBox::of_feature(feature);
            QuadBox::new(
                bounds.min[1] as f32,
                bounds.min[0] as f32,
                (bounds.max[1] - bounds.min[1]) as f32,
                (bounds.max[0] - bounds.min[0]) as f32,
            )
        }
    }

    /// The mode's test for a candidate whose box intersects `best`'s.
    fn overlaps(&self, best: &Feature, other: &Feature) -> Result<bool> {
        match self.mode {
            FeatureOverlapMode::ConvexHull => Ok(true),
            FeatureOverlapMode::TraceLevel => {
                let first = self.traces_of(best)?;
                let second = self.traces_of(other)?;
                Ok(first.iter().any(|a| {
                    second.iter().any(|b| {
                        !(a.rt_max < b.rt_min
                            || a.rt_min > b.rt_max
                            || a.mz_max < b.mz_min
                            || a.mz_min > b.mz_max)
                    })
                }))
            }
            FeatureOverlapMode::CentroidBased => {
                let tolerances = &self.tolerances;
                if tolerances.require_same_charge && best.charge != other.charge {
                    return Ok(false);
                }
                if tolerances.require_same_im {
                    let best_has = best.metadata.contains_key(FAIMS_CV);
                    let other_has = other.metadata.contains_key(FAIMS_CV);
                    if best_has != other_has {
                        return Ok(false);
                    }
                    if best_has && faims_cv(best)? != faims_cv(other)? {
                        return Ok(false);
                    }
                }
                let rt_diff = (best.rt - other.rt).abs();
                let mz_diff = (best.mz - other.mz).abs();
                Ok(rt_diff <= tolerances.rt_tolerance && mz_diff <= tolerances.mz_tolerance)
            }
        }
    }

    fn traces_of(&self, feature: &Feature) -> Result<&[TraceBounds]> {
        self.traces
            .get(&feature.unique_id)
            .map(Vec::as_slice)
            .ok_or_else(|| {
                Error::InvalidValue(format!(
                    "the feature with unique ID {} has no mass-trace bounds (the source dereferences a missing map entry)",
                    feature.unique_id
                ))
            })
    }
}

fn finite_box(value: QuadBox) -> bool {
    value.left.is_finite()
        && value.top.is_finite()
        && value.width.is_finite()
        && value.height.is_finite()
}

/// `getFeatureBounds` (FeatureOverlapFilter.cpp:31-94): the mass-trace bounds of
/// every feature, keyed by unique ID, so features sharing an ID pool theirs.
fn trace_bounds(features: &[Feature]) -> Result<BTreeMap<u64, Vec<TraceBounds>>> {
    let mut bounds: BTreeMap<u64, Vec<TraceBounds>> = BTreeMap::new();
    for feature in features {
        for (index, subordinate) in feature.subordinates.iter().enumerate() {
            let hull = feature.convex_hulls.get(index).ok_or_else(|| {
                Error::InvalidValue(format!(
                    "the feature with unique ID {} has {} subordinates but {} convex hulls (the source reads past its hull list)",
                    feature.unique_id,
                    feature.subordinates.len(),
                    feature.convex_hulls.len()
                ))
            })?;
            let points = hull.hull_points();
            // Defensive: `Plan::new` has already refused empty feature hulls in
            // the hull modes, so a non-empty hull always has outline points.
            let (Some(first), Some(last)) = (points.first(), points.last()) else {
                return Err(Error::MissingInformation(format!(
                    "convex hull {index} of the feature with unique ID {} is empty (the source reads its first point)",
                    feature.unique_id
                )));
            };
            let (mz_min, mz_max) = (first.mz, last.mz);
            let Some(trace_hull) = subordinate.convex_hulls.first() else {
                return Err(Error::MissingInformation(
                    "convex hulls for mass traces missing".into(),
                ));
            };
            let trace = trace_hull.hull_points();
            let (Some(trace_first), Some(trace_last)) = (trace.first(), trace.last()) else {
                continue;
            };
            let rt_min = trace
                .iter()
                .find(|point| point.mz > 0.0)
                .map_or(trace_last.rt, |point| point.rt);
            let mut rt_max = trace_first.rt;
            for point in trace.iter().rev() {
                if point.rt < rt_min {
                    break;
                }
                if point.mz > 0.0 {
                    rt_max = point.rt;
                    break;
                }
            }
            if rt_min > rt_max {
                continue;
            }
            bounds
                .entry(feature.unique_id)
                .or_default()
                .push(TraceBounds {
                    rt_min,
                    rt_max,
                    mz_min,
                    mz_max,
                });
        }
    }
    Ok(bounds)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `n` features at one position and charge with distinct unique IDs, so the
    /// first absorbs every other one.
    fn co_located(n: u64) -> Vec<Feature> {
        (0..n)
            .map(|k| {
                let mut feature = Feature::new(100.0, 500.0, 1.0);
                feature.charge = 2;
                feature.unique_id = k + 1;
                feature
            })
            .collect()
    }

    /// Runs the centroid-mode overlap loop as `run` does, without sorting, and
    /// returns the number of removed IDs and the journal.
    fn journal_of<O>(features: &mut [Feature], on_overlap: &mut O) -> (usize, Journal)
    where
        O: FnMut(&mut [Feature], usize, usize, &mut Journal) -> Result<bool>,
    {
        let tolerances = CentroidTolerances::default();
        let plan = Plan::new(features, FeatureOverlapMode::CentroidBased, &tolerances).unwrap();
        let mut journal = Journal::default();
        let removed = overlap_loop(features, &plan, on_overlap, &mut journal).unwrap();
        (removed.len(), journal)
    }

    #[test]
    fn merge_journal_records_each_field_once_however_many_merges() {
        let n = 2000;
        let mut features = co_located(n);
        let callback = FaimsMergeCallback::new(MergeIntensityMode::Sum, true);
        let (removed, journal) = journal_of(&mut features, &mut journaled_merge(callback));
        assert_eq!(removed as u64, n - 1);
        assert_eq!(features[0].intensity, n as f32);
        assert_eq!(
            features[0].metadata[MERGED_CENTROID_RTS]
                .as_float_list()
                .unwrap()
                .len() as u64,
            n
        );
        // merged_centroid_rts, merged_centroid_mzs and the intensity of the one
        // survivor, not three entries per merge.
        assert_eq!(journal.entries.len(), 3);
    }

    #[test]
    fn faims_journal_records_each_field_once_on_crafted_input() {
        // Every feature carries FAIMS_CV and a merged_centroid_IMs list, so the
        // survivor keeps its FAIMS_CV and absorbs every other voltage.
        let n = 500;
        let mut features = co_located(n);
        for (k, feature) in features.iter_mut().enumerate() {
            feature
                .metadata
                .insert(FAIMS_CV.to_owned(), MetaValue::from(k as i64));
            feature.metadata.insert(
                MERGED_CENTROID_IMS.to_owned(),
                MetaValue::try_from(vec![k as f64]).unwrap(),
            );
        }
        let (removed, journal) = journal_of(&mut features, &mut merge_different_voltages);
        assert_eq!(removed as u64, n - 1);
        assert_eq!(
            features[0].metadata[MERGED_CENTROID_IMS]
                .as_float_list()
                .unwrap()
                .len() as u64,
            n
        );
        // The two centroid lists, merged_centroid_IMs, FAIMS_merge_count and the
        // intensity.
        assert_eq!(journal.entries.len(), 5);
    }

    #[test]
    fn journal_rollback_restores_the_values_before_the_run() {
        let mut features = co_located(40);
        features[0].metadata.insert(
            MERGED_CENTROID_RTS.to_owned(),
            MetaValue::try_from(vec![7.0]).unwrap(),
        );
        features[0]
            .metadata
            .insert(FAIMS_CV.to_owned(), MetaValue::try_from(-45.0).unwrap());
        let before = features.clone();
        let callback = FaimsMergeCallback::new(MergeIntensityMode::Sum, true);
        let (removed, journal) = journal_of(&mut features, &mut journaled_merge(callback));
        assert_eq!(removed, 39);
        assert_ne!(features, before);
        // Every field once: both centroid lists, FAIMS_CV (removed by the first
        // merge), merged_centroid_IMs and FAIMS_merge_count (written by every
        // merge) and the intensity.
        assert_eq!(journal.entries.len(), 6);
        journal.rollback(&mut features);
        assert_eq!(features, before);
    }
}
