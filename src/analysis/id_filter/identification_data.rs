// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The `IdentificationData` filters of `PROCESSING/ID/IDFilter.h` (core
//! `bc9cc12`): `keepBestMatchPerObservation`,
//! `filterObservationMatchesByScore` and `removeDecoys`, over the native
//! identification graph ([`IdentificationData`]).

use crate::Result;
use crate::identification::graph::{IdentificationData, ScoreTypeId};
use std::collections::BTreeSet;

/// Keeps only the best observation match (e.g. PSM) per observation
/// (`keepBestMatchPerObservation`), then cleans up invalidated references.
///
/// "Best" is [`IdentificationData::best_matches`] with `require_score` false,
/// the source's `getBestMatchPerObservation(score_ref)` default: per
/// observation the first match with the best score of the given type, or the
/// only match when it has no such score. With at most one match nothing is
/// done, as in the source.
///
/// # Errors
///
/// Those of [`IdentificationData::best_matches`] (an unknown score type) and
/// of the graph removal; the graph is unchanged then.
pub fn keep_best_match_per_observation(
    data: &mut IdentificationData,
    score: ScoreTypeId,
) -> Result<()> {
    if data.observation_match_count() <= 1 {
        return Ok(());
    }
    let best: BTreeSet<_> = data.best_matches(score, false)?.into_iter().collect();
    data.remove_observation_matches_if(|id, _| !best.contains(&id))?;
    Ok(())
}

/// Removes the observation matches whose score of the given type is worse
/// than `cutoff`, and those without such a score
/// (`filterObservationMatchesByScore`), then cleans up.
///
/// A match is removed when `cutoff` is strictly better than its score in the
/// score type's orientation, so a score equal to the cutoff is kept.
///
/// # Errors
///
/// Those of [`IdentificationData::score_type`] (an unknown score type) and of
/// the graph removal; the graph is unchanged then.
pub fn filter_observation_matches_by_score(
    data: &mut IdentificationData,
    score: ScoreTypeId,
    cutoff: f64,
) -> Result<()> {
    let score_type = data.score_type(score)?.clone();
    data.remove_observation_matches_if(|_, observation_match| {
        match observation_match.result.score(score) {
            None => true,
            Some(value) => score_type.is_better_score(cutoff, value),
        }
    })?;
    Ok(())
}

/// Removes the parent sequences (proteins) annotated as decoys, then cleans
/// up if any was removed (`removeDecoys`).
///
/// # Errors
///
/// Those of the graph removal; the graph is unchanged then.
pub fn remove_decoys(data: &mut IdentificationData) -> Result<()> {
    data.remove_parent_sequences_if(|_, parent| parent.is_decoy)?;
    Ok(())
}
