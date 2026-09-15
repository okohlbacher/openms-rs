// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Overlap resolution and final annotation of the picked feature finder
//! (`FEATUREFINDER/FeatureFinderAlgorithmPicked.h`, step 4 of `run_`).
//!
//! The seed loop produces one candidate per surviving seed and per charge.
//! Candidates of different charges, and of one charge whose seeds were not
//! contained in an earlier feature, can still describe the same peptide, so the
//! source resolves them pairwise:
//!
//! - [`intersection`](crate::analysis::feature_finder_picked::resolution::intersection) measures how much two candidates overlap in retention
//!   time, as a fraction of the smaller one's total mass-trace width (source
//!   `intersection_`);
//! - [`resolve_overlaps`](crate::analysis::feature_finder_picked::resolution::resolve_overlaps) walks every pair whose m/z are close enough, keeps one
//!   of an overlapping pair by charge, intensity and quality, and moves the
//!   other into the winner's subordinates (source step 4);
//! - [`annotate_apex`](crate::analysis::feature_finder_picked::resolution::annotate_apex) records the apex scan of each feature (source's
//!   `spectrum_index` and `spectrum_native_id` loop).
//!
//! The module is serial, as the source is here. See
//! `docs/FEATURE_FINDER_PICKED_SUPPORT.md`.

use crate::kernel::{BoundingBox2D, Feature, MSExperiment};
use crate::metadata::MetaValue;
use crate::{Error, Result};

/// Meta key of the apex scan index, verbatim from the source.
pub const SPECTRUM_INDEX: &str = "spectrum_index";

/// Meta key of the apex scan's native identifier, verbatim from the source.
pub const SPECTRUM_NATIVE_ID: &str = "spectrum_native_id";

/// The warning the source logs when an apex scan index is out of range,
/// with the count substituted.
pub fn invalid_apex_warning(count: usize) -> String {
    format!(
        "Could not assign 'spectrum_native_id' for {count} feature(s), because the computed apex \
         spectrum index was out of range."
    )
}

/// The bounding box of every mass-trace hull of a feature: source
/// `f.getConvexHull().getBoundingBox()`.
///
/// The source builds the overall hull first, which is the single mass-trace
/// hull when there is one and the rectangle spanning every hull otherwise;
/// either way its bounding box is the union of the individual hull boxes, which
/// [`Feature::hull_bounding_box`] returns without building the hull.
fn overall_box(feature: &Feature) -> Option<BoundingBox2D> {
    feature.hull_bounding_box()
}

/// The overlap of two feature candidates in retention time: source
/// `intersection_`.
///
/// Each mass-trace hull contributes its retention-time extent. The numerator
/// sums, over every pair of hulls whose bounding boxes intersect, the length of
/// their retention-time overlap; the denominator is the smaller of the two
/// features' summed extents. The source distinguishes the four containment and
/// partial-overlap cases explicitly, and this keeps its case order, so a pair
/// matching two cases takes the first.
///
/// A feature without hulls gives a zero sum, and the quotient is then a
/// division by zero, as in the source.
///
/// # Native difference
///
/// A hull with no point has no bounding box here and is skipped. The source's
/// default `DBoundingBox` spans `[DBL_MAX, -DBL_MAX]`, whose `width()` is
/// negative infinity and poisons the sum. Every hull built by this algorithm
/// holds at least three points, so the case cannot arise on this path.
pub fn intersection(f1: &Feature, f2: &Feature) -> f64 {
    let boxes1: Vec<BoundingBox2D> = f1
        .convex_hulls
        .iter()
        .filter_map(crate::kernel::ConvexHull2D::bounding_box)
        .collect();
    let boxes2: Vec<BoundingBox2D> = f2
        .convex_hulls
        .iter()
        .filter_map(crate::kernel::ConvexHull2D::bounding_box)
        .collect();
    let s1: f64 = boxes1.iter().fold(0.0, |sum, b| sum + b.width());
    let s2: f64 = boxes2.iter().fold(0.0, |sum, b| sum + b.width());
    let mut overlap = 0.0;
    for bb1 in &boxes1 {
        for bb2 in &boxes2 {
            if !bb1.intersects(*bb2) {
                continue;
            }
            let (a0, a1) = (bb1.min().rt, bb1.max().rt);
            let (b0, b1) = (bb2.min().rt, bb2.max().rt);
            if a0 <= b0 && a1 >= b1 {
                overlap += bb2.width();
            } else if b0 <= a0 && b1 >= a1 {
                overlap += bb1.width();
            } else if a0 <= b0 && a1 <= b1 {
                overlap += a1 - b0;
            } else if b0 <= a0 && b1 <= a1 {
                overlap += b1 - a0;
            }
        }
    }
    // Source: `std::min(s1, s2)`, that is `(s2 < s1) ? s2 : s1`.
    overlap / if s2 < s1 { s2 } else { s1 }
}

/// Which of an overlapping pair the source keeps.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Keep {
    First,
    Second,
}

/// The source's decision for one overlapping pair, in its branch order.
///
/// Equal charges compare the `f32` product of intensity and overall quality,
/// with the first winning a tie. Otherwise, if one charge is a multiple of the
/// other, the higher charge wins — the source writes both branches with the
/// same debug text but different outcomes, and the second branch is reached
/// only when the first does not hold. Otherwise the higher overall quality
/// wins, the first winning a tie.
fn decide(f1: &Feature, f2: &Feature) -> Keep {
    if f1.charge == f2.charge {
        // Source: `float * float > float * float`, evaluated in binary32.
        return if f1.intensity * f1.quality > f2.intensity * f2.quality {
            Keep::First
        } else {
            Keep::Second
        };
    }
    // Source: `f2.getCharge() % f1.getCharge() == 0` keeps the second feature.
    // A zero charge would divide by zero there; the parameter minimum of
    // `isotopic_pattern:charge_low` is 1, so it cannot occur, and the guard
    // falls through to the quality rule instead of trapping.
    if f1.charge != 0 && f2.charge.checked_rem(f1.charge) == Some(0) {
        return Keep::Second;
    }
    // Source: `f1.getCharge() % f2.getCharge() == 0` keeps the first.
    if f2.charge != 0 && f1.charge.checked_rem(f2.charge) == Some(0) {
        return Keep::First;
    }
    if f1.quality > f2.quality {
        Keep::First
    } else {
        Keep::Second
    }
}

/// Resolve overlapping candidates in place: source step 4 up to the removal of
/// zero-intensity features.
///
/// `features` must be sorted by m/z, as the source sorts it first. Every pair
/// `(i, j)` with `i < j` is visited until `f2.mz - f1.mz` exceeds twice the
/// largest m/z extent of any feature, at which point the inner loop stops.
/// A pair is examined only while both intensities are non-zero and their
/// overall bounding boxes intersect. When [`intersection`] reaches
/// `feature:max_intersection`, the source's charge, intensity and quality rules
/// pick the winner, the loser is cloned
/// into the winner's subordinates and the loser's intensity is set to zero,
/// which both removes it later and skips it in the remaining pairs.
///
/// Returns the number of resolved pairs, which the source logs as "Removed N
/// overlapping features."
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a feature's subordinates cannot be
/// allocated.
pub fn resolve_overlaps(features: &mut [Feature], max_intersection: f64) -> Result<usize> {
    let boxes: Vec<Option<BoundingBox2D>> = features.iter().map(overall_box).collect();
    let max_mz_span =
        boxes.iter().flatten().fold(
            0.0,
            |max, b| if b.height() > max { b.height() } else { max },
        );
    let mut removed = 0usize;
    for i in 0..features.len() {
        for j in i + 1..features.len() {
            if features[j].mz - features[i].mz > 2.0 * max_mz_span {
                break;
            }
            if features[i].intensity == 0.0 || features[j].intensity == 0.0 {
                continue;
            }
            let (Some(bb1), Some(bb2)) = (boxes[i], boxes[j]) else {
                continue;
            };
            if !bb1.intersects(bb2) {
                continue;
            }
            // Source: `intersection >= max_feature_intersection_`, so a NaN
            // quotient (both features without hulls) leaves the pair alone.
            if intersection(&features[i], &features[j]) >= max_intersection {
                removed += 1;
                let keep = decide(&features[i], &features[j]);
                let (head, tail) = features.split_at_mut(j);
                let (f1, f2) = (&mut head[i], &mut tail[0]);
                let (winner, loser) = match keep {
                    Keep::First => (f1, f2),
                    Keep::Second => (f2, f1),
                };
                winner.subordinates.try_reserve(1).map_err(|_| {
                    Error::InvalidValue("cannot allocate a subordinate feature".into())
                })?;
                winner.subordinates.push(loser.clone());
                loser.intensity = 0.0;
            }
        }
    }
    Ok(removed)
}

/// Record the apex scan of every feature: the source's `spectrum_index` and
/// `spectrum_native_id` loop after step 4.
///
/// The index is `RTBegin(rt)`, the first scan at or after the feature's
/// retention time, so a feature beyond the last scan gets the scan count, which
/// addresses no scan; those features keep `spectrum_index` and receive no
/// native identifier, and the count of them is returned (the source logs it
/// through [`invalid_apex_warning`]).
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a feature's retention time is not
/// finite and [`Error::UnsortedData`] when the experiment is not sorted by
/// retention time, both from [`MSExperiment::rt_begin`]. The source's
/// `lower_bound` gives an unspecified index instead.
pub fn annotate_apex(features: &mut [Feature], experiment: &MSExperiment) -> Result<usize> {
    let mut invalid = 0usize;
    for feature in features.iter_mut() {
        let index = experiment.rt_begin(feature.rt)?;
        feature.metadata.insert(
            SPECTRUM_INDEX.into(),
            MetaValue::from(
                i64::try_from(index)
                    .map_err(|_| Error::InvalidValue("apex spectrum index exceeds i64".into()))?,
            ),
        );
        match experiment.spectra.get(index) {
            Some(spectrum) => {
                feature.metadata.insert(
                    SPECTRUM_NATIVE_ID.into(),
                    MetaValue::from(spectrum.native_id.clone()),
                );
            }
            None => invalid += 1,
        }
    }
    Ok(invalid)
}
