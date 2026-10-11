// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Hendrik Weisser, OpenMS Rust contributors $
//! `pickScoreType` and `applyToObservations` of `IdentificationData.h`.
//!
//! Both are templates in the source (`IdentificationData.h:498` and `:568`
//! at core bc9cc12). Here the first takes any iterator of
//! [`ScoredProcessingResult`]s, and the second any closure over an
//! [`Observation`].

use super::{
    GraphWork, IdentificationData, Observation, ScoreTypeId, ScoredProcessingResult, add, invalid,
    mul, record_size,
};
use crate::Result;
use std::cmp::Ordering;
use std::collections::BTreeMap;

impl IdentificationData {
    /// Pick a score type for an operation (for example filtering) over a
    /// collection of scored results — observation matches, identified
    /// peptides and so on (`pickScoreType`).
    ///
    /// Pass the results themselves, for example
    /// `data.observation_matches().map(|(_, m)| &m.result)`.
    ///
    /// # Parameters
    /// - `container`: the scored results, in the order the caller wants the
    ///   "first element" rule to see them.
    /// - `all_elements` (source default `false`): when `false`, only the first
    ///   element that has a score is consulted, which suffices when every
    ///   element went through the same processing. When `true`, the score
    ///   type supported by the most elements is chosen.
    /// - `any_score` (source default `false`): when `false`, each element
    ///   contributes only the primary score of its most recent scoring step
    ///   ([`ScoredProcessingResult::most_recent_score`], resolved through this
    ///   graph's processing software). When `true`, every score type of every
    ///   step of every element is counted, which implies `all_elements`.
    ///
    /// # Returns
    /// `Ok(None)` where the source returns `getScoreTypes().end()`: no element
    /// carried a score.
    ///
    /// # Errors
    /// A counted or returned score type that does not belong to this graph is
    /// refused, as is, when `any_score` is `false`, a processing step that
    /// does not (the source dereferences it to find the software's score
    /// priorities), and a container larger than the graph's work limit.
    ///
    /// # Native difference: the most-supported type is chosen by count
    /// The source counts elements per score type in a
    /// `std::map<ScoreTypeRef, Size>` and returns
    /// `max_element(score_counts.begin(), score_counts.end())->first`
    /// (`IdentificationData.h:599-602`). `max_element` compares the map's
    /// `std::pair`s, whose `operator<` orders by the key first, so the source
    /// returns the counted score type with the largest *key* — the highest
    /// memory address of a `std::set` node, because `ScoreTypeRef` orders by
    /// address — and the counts never decide. This port implements the
    /// documented contract instead: the highest count wins, and a tie goes to
    /// the later-registered score type. The two agree whenever the
    /// most-supported type is also the last registered one, which is what the
    /// source's only in-tree caller never reaches anyway:
    /// `MapAlignmentAlgorithmIdentification.cpp:120` uses the defaults
    /// (`all_elements == any_score == false`), where no counting happens.
    /// The executed Release build confirms the address rule (oracle
    /// `A5b-iddata`, cases `pick_all_count_b2` and `pick_any_count_a2`).
    pub fn pick_score_type<'a>(
        &self,
        container: impl IntoIterator<Item = &'a ScoredProcessingResult>,
        all_elements: bool,
        any_score: bool,
    ) -> Result<Option<ScoreTypeId>> {
        let mut work = GraphWork::new(self.limits);
        let mut counts: BTreeMap<ScoreTypeId, usize> = BTreeMap::new();
        for element in container {
            work.consume(add(element.steps_and_scores.len(), 1)?)?;
            if any_score {
                for step in &element.steps_and_scores {
                    work.consume(step.scores.len())?;
                    for &score in step.scores.keys() {
                        let count = counts.entry(score).or_insert(0);
                        *count = add(*count, 1)?;
                    }
                }
            } else {
                let recent = element.most_recent_score(|step| {
                    let software = self.processing_step(step)?.software;
                    Ok(self
                        .processing_software(software)?
                        .assigned_scores
                        .as_slice())
                })?;
                if let Some((score, _)) = recent {
                    if !all_elements {
                        self.score_type(score)?;
                        return Ok(Some(score));
                    }
                    let count = counts.entry(score).or_insert(0);
                    *count = add(*count, 1)?;
                }
            }
        }
        for &score in counts.keys() {
            self.score_type(score)?;
        }
        // Highest count; among equal counts the largest (latest-registered) ID.
        Ok(counts
            .into_iter()
            .max_by(|(left_id, left), (right_id, right)| {
                left.cmp(right).then_with(|| left_id.cmp(right_id))
            })
            .map(|(score, _)| score))
    }

    /// Apply `f` to every observation, in key order (`applyToObservations`).
    ///
    /// The closure may change anything, the key fields included: an
    /// observation whose data ID or input file changes moves to its new key
    /// position, and its [`ObservationId`](super::ObservationId) stays valid,
    /// so observation matches that reference it keep doing so. The source's
    /// caller, `MapAlignmentTransformer.cpp:194`, uses it to map every
    /// retention time through a transformation.
    ///
    /// # Errors
    /// The call is atomic: on any error the graph is left exactly as it was.
    /// Refused are a changed record that is invalid (an empty data ID, a
    /// non-finite coordinate, an input file that is not in this graph), one
    /// that would exceed the graph's size or work limits, and a change that
    /// gives two observations the same key.
    ///
    /// # Native difference: a key collision is refused, not erased
    /// The source applies the closure through Boost.MultiIndex `modify`,
    /// which *erases* an element whose modified key collides with another
    /// element's key. Observation matches that referenced the erased
    /// observation are left holding a dangling iterator.
    /// The executed Release build shows the loss (oracle `A5b-iddata`, case
    /// `apply_rename_collides`: three observations become two). This port
    /// refuses the whole call instead.
    pub fn apply_to_observations(&mut self, mut f: impl FnMut(&mut Observation)) -> Result<()> {
        self.transaction(|staged| {
            staged.operation(|this, work| {
                let slots = this.observations.order.clone();
                work.copy(mul(slots.len(), size_of::<usize>())?)?;
                for &slot in &slots {
                    let old = this.observations.sizes[&slot];
                    work.copy(old.bytes)?;
                    let mut record = this.observations.values[&slot].as_ref().clone();
                    f(&mut record);
                    this.input_file(record.input_file)?;
                    let size = record_size(record.measure(work)?, old.edges)?;
                    let total = this.checked_size(old, size, false)?;
                    this.observations.replace(slot, record, size);
                    this.retained = total;
                }
                let count = this.observations.order.len();
                let comparisons = usize::BITS as usize - count.leading_zeros() as usize + 1;
                work.consume(mul(
                    add(this.observations.max_bytes, 1)?,
                    mul(add(count, 1)?, comparisons)?,
                )?)?;
                let table = &mut this.observations;
                table
                    .order
                    .sort_by(|left, right| table.values[left].key_cmp(&table.values[right]));
                for pair in table.order.windows(2) {
                    let (left, right) = (&table.values[&pair[0]], &table.values[&pair[1]]);
                    if left.key_cmp(right) == Ordering::Equal {
                        return Err(invalid(format!(
                            "applying the function to the observations gives two of them the key \
                             '{}' in one input file; the source would silently erase one",
                            right.data_id
                        )));
                    }
                }
                Ok(())
            })
        })
    }
}
