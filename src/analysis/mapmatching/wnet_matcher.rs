// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Pairwise matching of two sets of two-dimensional points by Wasserstein
//! optimal transport.
//!
//! Native equivalent of `ANALYSIS/MAPMATCHING/WNetMatcher.h`. See
//! `docs/FEATURE_GROUPING_SUPPORT.md`.
//!
//! Two point sets, each with intensities, are matched by minimising the
//! transport cost of their weights; points left unmatched travel to a trash
//! node at `trash_cost`. The result is the one-to-one pairs a greedy pass
//! selects in order of descending flow, and the total cost of the optimal
//! transport. The transport network, scaling and solver are those of
//! [`wnet_aligner`](crate::analysis::mapmatching::wnet_aligner).
//!
//! The source's static member `match` is
//! [`WNetMatcher::match_points`](crate::analysis::mapmatching::wnet_matcher::WNetMatcher::match_points)
//! (`match` is a Rust keyword), and its default arguments are
//! [`DEFAULT_METRIC`](crate::analysis::mapmatching::wnet_matcher::DEFAULT_METRIC),
//! [`DEFAULT_MAX_DISTANCE`](crate::analysis::mapmatching::wnet_matcher::DEFAULT_MAX_DISTANCE) and
//! [`DEFAULT_TRASH_COST`](crate::analysis::mapmatching::wnet_matcher::DEFAULT_TRASH_COST).

use super::wnet_aligner::{Spectrum, WNetAligner};
use crate::error::{Error, Result};

pub use super::wnet_aligner::DistanceMetric;

/// The default metric of `match` (`DistanceMetric::LINF`).
pub const DEFAULT_METRIC: DistanceMetric = DistanceMetric::Linf;
/// The default `max_distance` of `match`.
pub const DEFAULT_MAX_DISTANCE: f64 = 100.0;
/// The default `trash_cost` of `match`.
pub const DEFAULT_TRASH_COST: f64 = 100.0;

/// The result of a pairwise matching (source `WNetMatcher::MatchResult`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MatchResult {
    /// `(index_a, index_b)` of the matched points.
    pub matched_pairs: Vec<(usize, usize)>,
    /// Total transport cost.
    pub cost: f64,
}

/// Pairwise point-set matching using Wasserstein optimal transport (source
/// `WNetMatcher`, a class of static members only).
#[derive(Clone, Copy, Debug, Default)]
pub struct WNetMatcher;

impl WNetMatcher {
    /// Convert `"L1"`, `"L2"` or `"LINF"` to a metric; anything else is
    /// `LINF` (source `metricFromString`, which also logs a warning for the
    /// other strings; see [`Self::metric_from_string_with_warning`]).
    pub fn metric_from_string(s: &str) -> DistanceMetric {
        Self::metric_from_string_with_warning(s).0
    }

    /// As [`Self::metric_from_string`], also returning the warning the source
    /// logs for an unrecognised string.
    pub fn metric_from_string_with_warning(s: &str) -> (DistanceMetric, Option<String>) {
        match s {
            "L1" => (DistanceMetric::L1, None),
            "L2" => (DistanceMetric::L2, None),
            "LINF" => (DistanceMetric::Linf, None),
            other => (
                DistanceMetric::Linf,
                Some(format!(
                    "WNetMatcher::metricFromString: unrecognized metric '{other}', falling back to LINF."
                )),
            ),
        }
    }

    /// Match two sets of 2D points using optimal transport (source `match`).
    ///
    /// Points are matched by minimising the total transport cost between their
    /// positions, weighted by intensity. Unmatched points are sent to a trash
    /// node at `trash_cost`; `max_distance` bounds the distance of a match.
    /// Either set empty gives an empty result with cost 0.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a position and an intensity
    /// vector differ in length (source `Exception::IllegalArgument`), and for
    /// the aligner's refusals: a negative `trash_cost`, intensities summing to
    /// zero or less, a value that does not fit 64 bits after scaling, or a
    /// network beyond the ceilings of
    /// [`wnet_aligner`](crate::analysis::mapmatching::wnet_aligner).
    pub fn match_points(
        positions_a: &[[f64; 2]],
        intensities_a: &[f64],
        positions_b: &[[f64; 2]],
        intensities_b: &[f64],
        metric: DistanceMetric,
        max_distance: f64,
        trash_cost: f64,
    ) -> Result<MatchResult> {
        if positions_a.len() != intensities_a.len() {
            return Err(Error::InvalidValue(
                "positions_a and intensities_a must have the same length".into(),
            ));
        }
        if positions_b.len() != intensities_b.len() {
            return Err(Error::InvalidValue(
                "positions_b and intensities_b must have the same length".into(),
            ));
        }
        if positions_a.is_empty() || positions_b.is_empty() {
            return Ok(MatchResult::default());
        }
        // spec_a is the empirical set, spec_b the theoretical target.
        let spec_a = Spectrum::new(positions_a.to_vec(), intensities_a.to_vec())?;
        let spec_b = Spectrum::new(positions_b.to_vec(), intensities_b.to_vec())?;
        let mut aligner = WNetAligner::new(&spec_a, &[&spec_b], metric, max_distance, trash_cost)?;
        aligner.set_point(&[1.0])?;
        let (empirical, theoretical) = aligner.consensus_for_target(0)?;
        Ok(MatchResult {
            matched_pairs: empirical.into_iter().zip(theoretical).collect(),
            cost: aligner.total_cost()?,
        })
    }

    /// [`Self::match_points`] with the source's default metric, maximum
    /// distance and trash cost.
    ///
    /// # Errors
    ///
    /// As [`Self::match_points`].
    pub fn match_with_defaults(
        positions_a: &[[f64; 2]],
        intensities_a: &[f64],
        positions_b: &[[f64; 2]],
        intensities_b: &[f64],
    ) -> Result<MatchResult> {
        Self::match_points(
            positions_a,
            intensities_a,
            positions_b,
            intensities_b,
            DEFAULT_METRIC,
            DEFAULT_MAX_DISTANCE,
            DEFAULT_TRASH_COST,
        )
    }
}
