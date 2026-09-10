// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Hendrik Weisser, OpenMS Rust contributors $
//! Referential filtering, with the source's ordered cleanup cascade.

use super::*;

/// The five independent switches of `IdentificationData::cleanup`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CleanupOptions {
    pub require_observation_match: bool,
    pub require_identified_sequence: bool,
    pub require_parent_match: bool,
    pub require_parent_group: bool,
    pub require_match_group: bool,
}
impl Default for CleanupOptions {
    fn default() -> Self {
        Self {
            require_observation_match: true,
            require_identified_sequence: true,
            require_parent_match: true,
            require_parent_group: false,
            require_match_group: false,
        }
    }
}
/// Counts include direct removals and their cascade. Scores are retained when a
/// nonempty group shrinks; the flags expose the source's associated warnings.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CleanupReport {
    pub removed_parents: usize,
    pub removed_peptides: usize,
    pub removed_oligos: usize,
    pub removed_compounds: usize,
    pub removed_observations: usize,
    pub removed_adducts: usize,
    pub removed_observation_matches: usize,
    pub removed_parent_links: usize,
    pub removed_parent_groups: usize,
    pub removed_match_groups: usize,
    pub parent_group_scores_may_be_invalid: bool,
    pub match_group_scores_may_be_invalid: bool,
}

fn tree_work(count: usize, work: &mut GraphWork) -> Result<()> {
    // Each pass uses at most two sparse-table/set lookups per visited item.
    let height = usize::BITS as usize - MAX_GRAPH_EDGES.leading_zeros() as usize + 1;
    work.consume(mul(count, mul(height, 2)?)?)
}
fn set<T>(count: usize, work: &mut GraphWork) -> Result<BTreeSet<T>> {
    tree_work(count, work)?;
    if count != 0 {
        work.allocation(add(512, mul(count, add(128, std::mem::size_of::<T>())?)?)?)?;
    }
    Ok(BTreeSet::new())
}
impl<T> Table<T> {
    fn remove_unless(
        &mut self,
        mut keep: impl FnMut(usize, &T) -> bool,
        retained: &mut Size,
        work: &mut GraphWork,
    ) -> Result<usize> {
        tree_work(self.order.len(), work)?;
        // Destruction may visit every owned payload byte, even for discarded records.
        let bytes = self
            .sizes
            .values()
            .try_fold(0, |sum, size| add(sum, size.bytes))?;
        work.consume(bytes)?;
        let before = self.order.len();
        self.order.retain(|slot| {
            if keep(*slot, &self.values[slot]) {
                true
            } else {
                self.values.remove(slot);
                let size = self.sizes.remove(slot).expect("live slot size");
                retained.bytes -= size.bytes;
                retained.edges -= size.edges;
                false
            }
        });
        self.max_bytes = self
            .sizes
            .values()
            .map(|s| s.bytes.saturating_sub(TABLE_OVERHEAD))
            .max()
            .unwrap_or(0);
        Ok(before - self.order.len())
    }
}
impl IdentificationData {
    /// Apply the source's single, ordered cascade atomically. Surviving IDs keep
    /// their values; removed IDs fail lookup permanently, including after reinsertion.
    pub fn cleanup(&mut self, options: CleanupOptions) -> Result<CleanupReport> {
        self.operation(|this, work| {
            let mut next = this.snapshot(work)?;
            let report = next.cleanup_in_place(options, work)?;
            *this = next;
            Ok(report)
        })
    }
    /// Remove selected matches, then run default cleanup only if one was removed.
    /// Predicate calls follow semantic match order. External predicate side effects
    /// are outside the graph transaction.
    pub fn remove_observation_matches_if(
        &mut self,
        mut predicate: impl FnMut(ObservationMatchId, &ObservationMatch) -> bool,
    ) -> Result<CleanupReport> {
        self.try_remove_observation_matches_if(|id, value| Ok(predicate(id, value)))
    }
    pub fn try_remove_observation_matches_if(
        &mut self,
        mut predicate: impl FnMut(ObservationMatchId, &ObservationMatch) -> Result<bool>,
    ) -> Result<CleanupReport> {
        self.operation(|this, work| {
            let mut removed = set(this.observation_match_count(), work)?;
            for (id, value) in this.observation_matches() {
                if predicate(id, value)? {
                    removed.insert(id.slot);
                }
            }
            if removed.is_empty() {
                return Ok(CleanupReport::default());
            }
            let mut next = this.snapshot(work)?;
            let count = next.observation_matches.remove_unless(
                |slot, _| !removed.contains(&slot),
                &mut next.retained,
                work,
            )?;
            let mut report = next.cleanup_in_place(CleanupOptions::default(), work)?;
            report.removed_observation_matches = add(report.removed_observation_matches, count)?;
            *this = next;
            Ok(report)
        })
    }
    /// Remove selected parents and perform the same conditional default cascade.
    pub fn remove_parent_sequences_if(
        &mut self,
        mut predicate: impl FnMut(ParentId, &ParentSequence) -> bool,
    ) -> Result<CleanupReport> {
        self.try_remove_parent_sequences_if(|id, value| Ok(predicate(id, value)))
    }
    pub fn try_remove_parent_sequences_if(
        &mut self,
        mut predicate: impl FnMut(ParentId, &ParentSequence) -> Result<bool>,
    ) -> Result<CleanupReport> {
        self.operation(|this, work| {
            let mut removed = set(this.parent_count(), work)?;
            for (id, value) in this.parents() {
                if predicate(id, value)? {
                    removed.insert(id.slot);
                }
            }
            if removed.is_empty() {
                return Ok(CleanupReport::default());
            }
            let mut next = this.snapshot(work)?;
            let count = next.parents.remove_unless(
                |slot, _| !removed.contains(&slot),
                &mut next.retained,
                work,
            )?;
            let mut report = next.cleanup_in_place(CleanupOptions::default(), work)?;
            report.removed_parents = add(report.removed_parents, count)?;
            *this = next;
            Ok(report)
        })
    }
    fn cleanup_in_place(
        &mut self,
        options: CleanupOptions,
        work: &mut GraphWork,
    ) -> Result<CleanupReport> {
        let mut report = CleanupReport::default();
        if options.require_parent_group {
            let mut used = set(self.retained.edges, work)?;
            for (_, groups) in self.parent_group_sets() {
                for group in &groups.groups {
                    used.extend(group.parent_refs.iter().map(|id| id.slot));
                }
            }
            report.removed_parents += self.parents.remove_unless(
                |slot, _| used.contains(&slot),
                &mut self.retained,
                work,
            )?;
        }
        // Parent links are repaired regardless of require_parent_match.
        macro_rules! repair_sequences {
            ($table:ident, $removed:ident) => {{
                tree_work(self.$table.values.len(), work)?;
                for (&slot, value) in &mut self.$table.values {
                    tree_work(value.parent_matches.len(), work)?;
                    let before = value.parent_matches.len();
                    value
                        .parent_matches
                        .retain(|id, _| self.parents.values.contains_key(&id.slot));
                    report.removed_parent_links += before - value.parent_matches.len();
                    if before != value.parent_matches.len() {
                        let size = record_size(
                            value.measure(work)?,
                            sequence_edges(&value.parent_matches, &value.result)?,
                        )?;
                        let old = self
                            .$table
                            .sizes
                            .insert(slot, size)
                            .expect("live sequence size");
                        self.retained.bytes = add(self.retained.bytes - old.bytes, size.bytes)?;
                        self.retained.edges = add(self.retained.edges - old.edges, size.edges)?;
                    }
                }
                if options.require_parent_match {
                    report.$removed += self.$table.remove_unless(
                        |_, value| !value.parent_matches.is_empty(),
                        &mut self.retained,
                        work,
                    )?;
                }
            }};
        }
        repair_sequences!(peptides, removed_peptides);
        repair_sequences!(oligos, removed_oligos);
        report.removed_observation_matches += self.observation_matches.remove_unless(
            |_, value| match value.identified_molecule {
                IdentifiedMolecule::Peptide(id) => self.peptides.values.contains_key(&id.slot),
                IdentifiedMolecule::Compound(id) => self.compounds.values.contains_key(&id.slot),
                IdentifiedMolecule::Oligo(id) => self.oligos.values.contains_key(&id.slot),
            },
            &mut self.retained,
            work,
        )?;
        if options.require_match_group {
            let mut used = set(self.retained.edges, work)?;
            for (_, group) in self.observation_match_groups() {
                used.extend(group.observation_match_refs.iter().map(|id| id.slot));
            }
            report.removed_observation_matches += self.observation_matches.remove_unless(
                |slot, _| used.contains(&slot),
                &mut self.retained,
                work,
            )?;
        }
        if options.require_observation_match {
            let n = self.observation_match_count();
            let mut observations = set(n, work)?;
            let mut peptides = set(n, work)?;
            let mut compounds = set(n, work)?;
            let mut oligos = set(n, work)?;
            let mut adducts = set(n, work)?;
            for (_, value) in self.observation_matches() {
                observations.insert(value.observation.slot);
                match value.identified_molecule {
                    IdentifiedMolecule::Peptide(id) => {
                        peptides.insert(id.slot);
                    }
                    IdentifiedMolecule::Compound(id) => {
                        compounds.insert(id.slot);
                    }
                    IdentifiedMolecule::Oligo(id) => {
                        oligos.insert(id.slot);
                    }
                }
                if let Some(id) = value.adduct {
                    adducts.insert(id.slot);
                }
            }
            macro_rules! unused {
                ($table:ident, $field:ident) => {
                    report.$field += self.$table.remove_unless(
                        |slot, _| $table.contains(&slot),
                        &mut self.retained,
                        work,
                    )?;
                };
            }
            unused!(observations, removed_observations);
            unused!(peptides, removed_peptides);
            unused!(compounds, removed_compounds);
            unused!(oligos, removed_oligos);
            unused!(adducts, removed_adducts);
        }
        if options.require_identified_sequence {
            let mut used = set(self.retained.edges, work)?;
            for (_, value) in self.peptides() {
                used.extend(value.parent_matches.keys().map(|id| id.slot));
            }
            for (_, value) in self.oligos() {
                used.extend(value.parent_matches.keys().map(|id| id.slot));
            }
            report.removed_parents += self.parents.remove_unless(
                |slot, _| used.contains(&slot),
                &mut self.retained,
                work,
            )?;
        }
        // Empty parent grouping operations are retained. Key collisions after
        // pruning retain the first original key's payload, without dangling iterators.
        tree_work(self.parent_group_nodes, work)?;
        for (&slot, groups) in &mut self.parent_group_sets.values {
            let before = groups.groups.len();
            for group in &mut groups.groups {
                tree_work(group.parent_refs.len(), work)?;
                let old_len = group.parent_refs.len();
                group
                    .parent_refs
                    .retain(|id| self.parents.values.contains_key(&id.slot));
                report.parent_group_scores_may_be_invalid |=
                    !group.parent_refs.is_empty() && group.parent_refs.len() != old_len;
            }
            groups.groups.retain(|group| !group.parent_refs.is_empty());
            groups.normalize(work)?;
            report.removed_parent_groups += before - groups.groups.len();
            let size = record_size(groups.measure(work)?, parent_group_set_edges(groups)?)?;
            let old = self
                .parent_group_sets
                .sizes
                .insert(slot, size)
                .expect("live grouping size");
            self.retained.bytes = add(self.retained.bytes - old.bytes, size.bytes)?;
            self.retained.edges = add(self.retained.edges - old.edges, size.edges)?;
        }
        self.parent_group_nodes -= report.removed_parent_groups;
        for (&slot, group) in &mut self.observation_match_groups.values {
            tree_work(group.observation_match_refs.len(), work)?;
            let old_len = group.observation_match_refs.len();
            group
                .observation_match_refs
                .retain(|id| self.observation_matches.values.contains_key(&id.slot));
            report.match_group_scores_may_be_invalid |= !group.observation_match_refs.is_empty()
                && group.observation_match_refs.len() != old_len;
            let size = record_size(
                group.measure(work)?,
                add(
                    group.observation_match_refs.len(),
                    scored_edges(&group.result)?,
                )?,
            )?;
            let old = self
                .observation_match_groups
                .sizes
                .insert(slot, size)
                .expect("live match grouping size");
            self.retained.bytes = add(self.retained.bytes - old.bytes, size.bytes)?;
            self.retained.edges = add(self.retained.edges - old.edges, size.edges)?;
        }
        report.removed_match_groups += self.observation_match_groups.remove_unless(
            |_, group| !group.observation_match_refs.is_empty(),
            &mut self.retained,
            work,
        )?;
        // Stable sort keeps the source's previous key order as the collision tie-break.
        let table = &mut self.observation_match_groups;
        let n = table.order.len();
        let max_refs = table
            .values
            .values()
            .map(|g| g.observation_match_refs.len())
            .max()
            .unwrap_or(0);
        let height = usize::BITS as usize - n.leading_zeros() as usize + 1;
        work.consume(mul(mul(mul(n, height)?, 8)?, add(max_refs, 1)?)?)?;
        work.allocation(mul(n, std::mem::size_of::<usize>())?)?;
        table
            .order
            .sort_by(|a, b| table.values[a].key_cmp(&table.values[b]));
        let mut duplicates = set(n, work)?;
        tree_work(mul(n, add(max_refs, 1)?)?, work)?;
        for pair in table.order.windows(2) {
            if table.values[&pair[0]].observation_match_refs
                == table.values[&pair[1]].observation_match_refs
            {
                duplicates.insert(pair[1]);
            }
        }
        report.removed_match_groups += table.remove_unless(
            |slot, _| !duplicates.contains(&slot),
            &mut self.retained,
            work,
        )?;
        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sparse_slot_exhaustion_is_checked_without_mutating_live_records() {
        let mut graph = IdentificationData::new().unwrap();
        let live = graph
            .register_parent_sequence(ParentSequence::new("live"))
            .unwrap();
        graph.parents.next_slot = usize::MAX;
        let error = graph
            .register_parent_sequence(ParentSequence::new("new"))
            .unwrap_err();
        assert!(error.to_string().contains("slot IDs exhausted"));
        assert_eq!(graph.parent_count(), 1);
        assert_eq!(graph.parent(live).unwrap().accession, "live");
        assert_eq!(graph.parents.next_slot, usize::MAX);
    }
}
