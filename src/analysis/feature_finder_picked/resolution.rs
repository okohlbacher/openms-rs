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

use std::cmp::Ordering;

use crate::analysis::feature_finder_picked::debug::{LogSink, NoLog, g, put_all};
use crate::analysis::feature_finder_picked::scoring::libstdcxx;
use crate::kernel::{ConvexHull2D, Feature, MSExperiment};
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

/// The source's `DBoundingBox<2>`: an interval per dimension that starts
/// empty, with its minimum at `DBL_MAX` and its maximum at `-DBL_MAX`.
///
/// The native `BoundingBox2D` cannot be empty, but the source's arithmetic on
/// an empty box is observable once a caller's feature has an empty hull: its
/// `width()` is `-DBL_MAX - DBL_MAX = -inf`, it intersects nothing, and a
/// feature with several hulls, one of them empty, gets an overall box spanning
/// `[-DBL_MAX, DBL_MAX]` in both dimensions (`Feature::getConvexHull`,
/// `Feature.cpp:93-136`).
#[derive(Clone, Copy, Debug, PartialEq)]
struct SourceBox {
    min: [f64; 2],
    max: [f64; 2],
}

impl SourceBox {
    const EMPTY: Self = Self {
        min: [f64::MAX, f64::MAX],
        max: [f64::MIN, f64::MIN],
    };

    /// `DBoundingBox::enlarge`.
    fn enlarge(&mut self, rt: f64, mz: f64) {
        for (dim, value) in [rt, mz].into_iter().enumerate() {
            if value < self.min[dim] {
                self.min[dim] = value;
            }
            if value > self.max[dim] {
                self.max[dim] = value;
            }
        }
    }

    /// `ConvexHull2D::getBoundingBox`: the scan envelopes, or the outline of
    /// an outline-only hull.
    fn of_hull(hull: &ConvexHull2D) -> Self {
        let mut bounds = Self::EMPTY;
        if let Some(native) = hull.bounding_box() {
            bounds.enlarge(native.min().rt, native.min().mz);
            bounds.enlarge(native.max().rt, native.max().mz);
        }
        bounds
    }

    /// `f.getConvexHull().getBoundingBox()`: one hull's box, the rectangle
    /// over the hull boxes' corners for several hulls, and an empty box for
    /// none.
    fn of_feature(feature: &Feature) -> Self {
        match feature.convex_hulls.as_slice() {
            [] => Self::EMPTY,
            [hull] => Self::of_hull(hull),
            hulls => {
                let mut bounds = Self::EMPTY;
                for hull in hulls {
                    let hull_box = Self::of_hull(hull);
                    bounds.enlarge(hull_box.min[0], hull_box.min[1]);
                    bounds.enlarge(hull_box.max[0], hull_box.max[1]);
                }
                bounds
            }
        }
    }

    /// `DIntervalBase::width`: the retention-time extent.
    fn width(&self) -> f64 {
        self.max[0] - self.min[0]
    }

    /// `DIntervalBase::height`: the m/z extent.
    fn height(&self) -> f64 {
        self.max[1] - self.min[1]
    }

    /// `DBoundingBox::intersects`.
    fn intersects(&self, other: &Self) -> bool {
        // `!(a > b) && !(c < d)`: an unordered pair (NaN) counts as intersecting.
        (0..2).all(|dim| {
            other.min[dim].partial_cmp(&self.max[dim]) != Some(Ordering::Greater)
                && other.max[dim].partial_cmp(&self.min[dim]) != Some(Ordering::Less)
        })
    }
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
/// A hull with no point has the source's empty bounding box: its width
/// `-inf` enters the sum, and it intersects no other box. The algorithm never
/// builds such a hull, but a caller's map may hold one.
pub fn intersection(f1: &Feature, f2: &Feature) -> f64 {
    let boxes1: Vec<SourceBox> = f1.convex_hulls.iter().map(SourceBox::of_hull).collect();
    let boxes2: Vec<SourceBox> = f2.convex_hulls.iter().map(SourceBox::of_hull).collect();
    let s1: f64 = boxes1.iter().fold(0.0, |sum, b| sum + b.width());
    let s2: f64 = boxes2.iter().fold(0.0, |sum, b| sum + b.width());
    let mut overlap = 0.0;
    for bb1 in &boxes1 {
        for bb2 in &boxes2 {
            if !bb1.intersects(bb2) {
                continue;
            }
            let (a0, a1) = (bb1.min[0], bb1.max[0]);
            let (b0, b1) = (bb2.min[0], bb2.max[0]);
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

/// The branch of the source's decision, which picks the debug text.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Rule {
    SameCharge,
    Multiple,
    Quality,
}

/// The source's decision for one overlapping pair, in its branch order.
///
/// Equal charges compare the `f32` product of intensity and overall quality,
/// with the first winning a tie. Otherwise, if one charge is a multiple of the
/// other, the higher charge wins — the source writes both branches with the
/// same debug text but different outcomes, and the second branch is reached
/// only when the first does not hold. Otherwise the higher overall quality
/// wins, the first winning a tie.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] where the source's `int` remainder traps:
/// `f2.charge % f1.charge` with `f1.charge == 0` (or `i32::MIN % -1`), and the
/// same for the second remainder, which the source evaluates only when the
/// first is not zero. The x86_64 `idiv` raises `SIGFPE` there, which ends the
/// process. The algorithm's own features have charges of at least 1; a
/// caller's map can hold charge 0, the featureXML default.
fn decide(f1: &Feature, f2: &Feature) -> Result<(Keep, Rule)> {
    if f1.charge == f2.charge {
        // Source: `float * float > float * float`, evaluated in binary32.
        let keep = if f1.intensity * f1.quality > f2.intensity * f2.quality {
            Keep::First
        } else {
            Keep::Second
        };
        return Ok((keep, Rule::SameCharge));
    }
    let trap = |a: i32, b: i32| {
        Error::InvalidValue(format!(
            "overlap resolution computes the charge remainder {a} % {b}, which traps (SIGFPE) in \
             the source and ends the process"
        ))
    };
    // Source: `f2.getCharge() % f1.getCharge() == 0` keeps the second feature.
    let first = f2
        .charge
        .checked_rem(f1.charge)
        .ok_or_else(|| trap(f2.charge, f1.charge))?;
    if first == 0 {
        return Ok((Keep::Second, Rule::Multiple));
    }
    // Source: `f1.getCharge() % f2.getCharge() == 0` keeps the first.
    let second = f1
        .charge
        .checked_rem(f2.charge)
        .ok_or_else(|| trap(f1.charge, f2.charge))?;
    if second == 0 {
        return Ok((Keep::First, Rule::Multiple));
    }
    let keep = if f1.quality > f2.quality {
        Keep::First
    } else {
        Keep::Second
    };
    Ok((keep, Rule::Quality))
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
/// allocated, and where the charge rule traps (see `decide`), which only a
/// caller's feature of charge 0 can cause.
pub fn resolve_overlaps(features: &mut [Feature], max_intersection: f64) -> Result<usize> {
    resolve_overlaps_logged(features, max_intersection, &mut NoLog, &mut |_| Ok(()))
}

/// [`resolve_overlaps`] writing the source's debug lines to `log` and calling
/// `progress` with `i * n + j` for every pair it visits, before the m/z
/// cut-off test, as the source calls `setProgress`; the value is computed in
/// `size_t` arithmetic, wrapping as the source's does.
pub(crate) fn resolve_overlaps_logged<L: LogSink>(
    features: &mut [Feature],
    max_intersection: f64,
    log: &mut L,
    progress: &mut dyn FnMut(u64) -> Result<()>,
) -> Result<usize> {
    let count = features.len();
    let boxes: Vec<SourceBox> = features.iter().map(SourceBox::of_feature).collect();
    let max_mz_span = boxes.iter().fold(
        0.0,
        |max, b| if b.height() > max { b.height() } else { max },
    );
    let mut removed = 0usize;
    for i in 0..count {
        for j in i + 1..count {
            progress((i as u64).wrapping_mul(count as u64).wrapping_add(j as u64))?;
            if features[j].mz - features[i].mz > 2.0 * max_mz_span {
                break;
            }
            if features[i].intensity == 0.0 || features[j].intensity == 0.0 {
                continue;
            }
            if !boxes[i].intersects(&boxes[j]) {
                continue;
            }
            // Source: `intersection >= max_feature_intersection_`, so a NaN
            // quotient (both features without hulls) leaves the pair alone.
            let overlap = intersection(&features[i], &features[j]);
            if overlap >= max_intersection {
                removed += 1;
                if log.enabled() {
                    put_all(
                        log,
                        &[
                            " - Intersection (",
                            &(i + 1).to_string(),
                            "/",
                            &(j + 1).to_string(),
                            "): ",
                            &g(overlap),
                            "\n",
                        ],
                    );
                }
                let (keep, rule) = decide(&features[i], &features[j])?;
                if log.enabled() {
                    let removed_one = if keep == Keep::First { j } else { i };
                    let (text, reported) = match rule {
                        Rule::SameCharge => {
                            ("   - same charge -> removing duplicate ", removed_one)
                        }
                        // Source: both multiple-of branches print the first
                        // index, although the second one removes the second
                        // feature.
                        Rule::Multiple => (
                            "   - different charge (one is the multiple of the other) -> removing \
                             lower charge ",
                            i,
                        ),
                        Rule::Quality => (
                            "   - different charge -> removing lower score ",
                            removed_one,
                        ),
                    };
                    put_all(log, &[text, &(reported + 1).to_string(), "\n"]);
                }
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
/// `RTBegin` is `std::lower_bound` with `RTLess`, which is well defined for any
/// retention time on sorted spectra: a NaN compares false with every scan and
/// gives index 0, `+inf` gives the scan count. A caller's feature can carry
/// such a value, and it is annotated as the source annotates it.
///
/// The search is libstdc++'s `std::lower_bound` as the Release build runs it
/// (the crate-private `libstdcxx::lower_bound`), also where a NaN retention
/// time in the experiment or of the feature leaves the keys unpartitioned; it
/// never leaves the scans.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] only when an index exceeds `i64`.
pub fn annotate_apex(features: &mut [Feature], experiment: &MSExperiment) -> Result<usize> {
    let mut invalid = 0usize;
    for feature in features.iter_mut() {
        // Source `map_.RTBegin(rt)`.
        let index =
            libstdcxx::lower_bound(&experiment.spectra, |spectrum| spectrum.rt < feature.rt);
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
