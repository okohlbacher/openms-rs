// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Hendrik Weisser, OpenMS Rust contributors $
//! Parent and observation-match group payloads from Core SDK 6bfc0e4.

use super::records::{Meter, merge_cost, prepare_scored_merge};
use super::{
    GraphWork, IdentifiedMolecule, ObservationId, ObservationMatchId, ParentId, ScoreTypeId,
    ScoredProcessingResult,
};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

/// A parent-reference set with score values. Identity within a group set ignores scores.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParentGroup {
    pub scores: BTreeMap<ScoreTypeId, f64>,
    pub parent_refs: BTreeSet<ParentId>,
}
impl ParentGroup {
    pub fn new(parent_refs: BTreeSet<ParentId>) -> Self {
        Self {
            parent_refs,
            scores: BTreeMap::new(),
        }
    }
    pub(crate) fn measure(&self, work: &mut GraphWork) -> Result<usize> {
        let mut meter = Meter::new::<Self>(work)?;
        if !self.scores.is_empty() {
            meter.add(512)?;
        }
        meter.scores(&self.scores)?;
        if !self.parent_refs.is_empty() {
            meter.add(512)?;
        }
        meter.slots::<ParentId>(self.parent_refs.len())?;
        Ok(meter.bytes)
    }
}

/// Ordered results of one grouping operation. Registration always appends a set,
/// even when its label/content repeats. Groups are normalized by parent-reference
/// set on registration, preserving the first score map for each duplicate key.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParentGroupSet {
    pub label: String,
    pub groups: Vec<ParentGroup>,
    pub result: ScoredProcessingResult,
}
impl ParentGroupSet {
    pub fn new(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            ..Self::default()
        }
    }
    pub(crate) fn measure(&self, work: &mut GraphWork) -> Result<usize> {
        let mut meter = Meter::new::<Self>(work)?;
        meter.text(&self.label)?;
        meter.scored(&self.result)?;
        meter.slots::<ParentGroup>(self.groups.len())?;
        for group in &self.groups {
            let bytes = group.measure(meter.work)?;
            meter.add(bytes)?;
        }
        Ok(meter.bytes)
    }
    pub(crate) fn normalize(&mut self, work: &mut GraphWork) -> Result<()> {
        let count = self.groups.len();
        work.consume(count)?;
        let max_key = self
            .groups
            .iter()
            .map(|group| group.parent_refs.len())
            .max()
            .unwrap_or(0);
        let comparisons = usize::BITS as usize - count.leading_zeros() as usize + 1;
        let cost = count
            .checked_mul(comparisons)
            .and_then(|n| n.checked_mul(8))
            .and_then(|n| n.checked_mul(max_key.saturating_add(1)))
            .ok_or_else(overflow)?;
        work.consume(cost)?;
        work.copy(
            count
                .checked_mul(std::mem::size_of::<ParentGroup>())
                .ok_or_else(overflow)?,
        )?;
        // Stable sorting is required: equal parent sets retain the first score map.
        self.groups
            .sort_by(|a, b| a.parent_refs.cmp(&b.parent_refs));
        let refs = self.groups.iter().try_fold(0usize, |sum, g| {
            sum.checked_add(g.parent_refs.len()).ok_or_else(overflow)
        })?;
        work.consume(
            refs.checked_mul(2)
                .and_then(|n| n.checked_add(count))
                .ok_or_else(overflow)?,
        )?;
        self.groups.dedup_by(|a, b| a.parent_refs == b.parent_refs);
        Ok(())
    }
}

/// Related observation matches, keyed by the complete reference set. Native
/// equality includes metadata as well as references and ordered score history.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ObservationMatchGroup {
    pub observation_match_refs: BTreeSet<ObservationMatchId>,
    pub result: ScoredProcessingResult,
}
impl ObservationMatchGroup {
    pub fn new(refs: BTreeSet<ObservationMatchId>) -> Self {
        Self {
            observation_match_refs: refs,
            result: ScoredProcessingResult::default(),
        }
    }
    pub(crate) fn key_cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.observation_match_refs
            .cmp(&other.observation_match_refs)
    }
    pub(crate) fn measure(&self, work: &mut GraphWork) -> Result<usize> {
        let mut meter = Meter::new::<Self>(work)?;
        if !self.observation_match_refs.is_empty() {
            meter.add(512)?;
        }
        meter.slots::<ObservationMatchId>(self.observation_match_refs.len())?;
        meter.scored(&self.result)?;
        Ok(meter.bytes)
    }
    /// Merge only scored payload, retaining this group's membership.
    pub fn merge(&mut self, other: &Self) -> Result<()> {
        self.merge_with_work(other, &mut GraphWork::default())
    }
    pub(crate) fn merge_with_work(&mut self, other: &Self, work: &mut GraphWork) -> Result<()> {
        let bytes = self
            .measure(work)?
            .checked_add(other.measure(work)?)
            .ok_or_else(overflow)?;
        merge_cost(
            work,
            bytes,
            self.result
                .metadata
                .len()
                .saturating_add(other.result.metadata.len()),
        )?;
        prepare_scored_merge(&self.result, &other.result, work)?;
        work.copy(bytes)?;
        let mut next = self.clone();
        next.result.merge_unchecked(&other.result);
        next.measure(work)?;
        *self = next;
        Ok(())
    }
    /// Empty and singleton sets are true, without resolving a reference, as in the source.
    pub fn all_same_molecule(
        &self,
        resolve: impl FnMut(ObservationMatchId) -> Result<IdentifiedMolecule>,
    ) -> Result<bool> {
        self.all_same(resolve, &mut GraphWork::default())
    }
    /// "Query" means the exact observation reference, not equal coordinates or data IDs.
    pub fn all_same_query(
        &self,
        resolve: impl FnMut(ObservationMatchId) -> Result<ObservationId>,
    ) -> Result<bool> {
        self.all_same(resolve, &mut GraphWork::default())
    }
    pub(super) fn all_same<T: Eq>(
        &self,
        mut resolve: impl FnMut(ObservationMatchId) -> Result<T>,
        work: &mut GraphWork,
    ) -> Result<bool> {
        work.consume(self.observation_match_refs.len())?;
        let mut refs = self.observation_match_refs.iter();
        if self.observation_match_refs.len() <= 1 {
            return Ok(true);
        }
        let first = resolve(*refs.next().expect("two or more members"))?;
        for &id in refs {
            if resolve(id)? != first {
                return Ok(false);
            }
        }
        Ok(true)
    }
}
fn overflow() -> Error {
    Error::InvalidValue("identification group size overflow".into())
}
