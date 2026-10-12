// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `ANALYSIS/ID/PeptideProteinResolution.h`: resolve shared peptides
//! between indistinguishable protein groups by protein score.

use crate::analysis::id_filter::{remove_unreferenced_proteins, update_protein_groups};
use crate::identification::{PeptideIdentification, ProteinGroup, ProteinIdentification};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::fmt;

/// Source `Exception::MissingInformation`.
fn missing(message: impl Into<String>) -> Error {
    Error::MissingInformation(message.into())
}

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// A connected component of the bipartite protein-group/peptide graph: the
/// indices of the (indistinguishable) protein groups and of the peptide
/// identifications in it (source `ConnectedComponent`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ConnectedComponent {
    /// Indices into the protein run's indistinguishable groups.
    pub prot_grp_indices: BTreeSet<usize>,
    /// Indices into the peptide identifications.
    pub pep_indices: BTreeSet<usize>,
}

impl fmt::Display for ConnectedComponent {
    /// Source `operator<<`: `Proteins: 0,1,` newline `Peptides: 2,3,`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Proteins: ")?;
        for index in &self.prot_grp_indices {
            write!(f, "{index},")?;
        }
        write!(f, "\nPeptides: ")?;
        for index in &self.pep_indices {
            write!(f, "{index},")?;
        }
        Ok(())
    }
}

/// Resolves shared peptides based on protein scores (source
/// `PeptideProteinResolution`).
///
/// Resolves the connected components of the bipartite protein-group/peptide
/// graph by protein group probability and adds each resolved component as an
/// additional protein group of the run. Shared peptides of a component are
/// greedily assigned to the proteins of the current best indistinguishable
/// group until every peptide is uniquely assigned, by removing all other
/// evidence from the peptides. This lets `ProteinQuantifier` use more
/// peptides, at the cost of potentially noisier peptide quantities.
///
/// As in most protein inference tools, only the best hit (the first) of a
/// peptide identification is considered. Probability ties are resolved by
/// taking the group with more peptides.
///
/// Every operation returns the messages the source logs (information and
/// warnings, and with `statistics` the component statistics).
///
/// The source's improvement notes: the class could provide an iterator over
/// components, consider all peptide hits (a tripartite graph), work with MS1
/// features, and separate resolution from adding the groups.
#[derive(Clone, Debug, Default)]
pub struct PeptideProteinResolution {
    indist_prot_grp_to_pep: BTreeMap<usize, BTreeSet<usize>>,
    pep_to_indist_prot_grp: BTreeMap<usize, BTreeSet<usize>>,
    prot_acc_to_indist_prot_grp: BTreeMap<String, usize>,
    statistics: bool,
}

impl PeptideProteinResolution {
    /// An empty resolver; `statistics` adds component statistics to the
    /// returned messages.
    pub fn new(statistics: bool) -> Self {
        Self {
            statistics,
            ..Self::default()
        }
    }

    /// Peptide-centric resolution without building the graph (source static
    /// `resolve`).
    ///
    /// For the best hit of every peptide identification, the evidence is
    /// reduced to that of the best target groups and the best decoy groups
    /// among the groups it maps to; a group is a decoy group when any of its
    /// accessions belongs to a decoy protein hit. Between the two, the target
    /// groups win when `targets_first`, when there are no decoy groups, or
    /// when their probability is higher; the decoy groups win when there are
    /// no target groups or theirs is higher; on equal probability both are
    /// kept. With `resolve_ties`, only the evidence of the lowest-indexed
    /// remaining group is kept.
    ///
    /// Within the target and within the decoy groups, "best" is the *lowest*
    /// group probability, as the source compares with `<`, while the choice
    /// between targets and decoys prefers the higher probability. This
    /// reproduces the source; it looks like a defect there.
    ///
    /// All peptides are used, whatever protein run they belong to (a source
    /// to-do).
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] (source
    /// `Exception::MissingInformation`) when the run has no indistinguishable
    /// groups ("Currently this class only resolves across groups") or a best
    /// hit references an accession that is in no group (add such proteins as
    /// singletons), and [`Error::InvalidValue`] when a hit's target/decoy
    /// annotation is invalid. `peptides` is unchanged on error.
    pub fn resolve(
        protein: &ProteinIdentification,
        peptides: &mut [PeptideIdentification],
        resolve_ties: bool,
        targets_first: bool,
    ) -> Result<Vec<String>> {
        let groups = &protein.indistinguishable_groups;
        if groups.is_empty() {
            return Err(no_groups());
        }
        let mut messages = vec![resolving_message(protein)];
        let mut decoy_accessions = BTreeSet::new();
        for hit in &protein.hits {
            if hit.is_decoy()? {
                decoy_accessions.insert(hit.accession.as_str());
            }
        }
        let mut group_is_decoy = vec![false; groups.len()];
        let mut accession_to_group: BTreeMap<&str, usize> = BTreeMap::new();
        for (index, group) in groups.iter().enumerate() {
            for accession in &group.accessions {
                accession_to_group.insert(accession, index);
                if decoy_accessions.contains(accession.as_str()) {
                    group_is_decoy[index] = true;
                }
            }
        }
        let probability = |index: usize| groups[index].probability;
        let mut next = peptides.to_vec();
        for pep in &mut next {
            let Some(best_hit) = pep.hits.first_mut() else {
                messages.push(
                    "Warning PeptideProteinResolution: Skipping spectrum without hits.".to_owned(),
                );
                continue;
            };
            let mut best_target: BTreeSet<usize> = BTreeSet::new();
            let mut best_decoy: BTreeSet<usize> = BTreeSet::new();
            let mut group_to_evidence: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
            for (evidence_index, evidence) in best_hit.evidences.iter().enumerate() {
                let &group = accession_to_group
                    .get(evidence.protein_accession.as_str())
                    .ok_or_else(|| {
                        missing("Not all proteins present in an indistinguishable group. Make sure to add them as singletons.")
                    })?;
                group_to_evidence
                    .entry(group)
                    .or_default()
                    .insert(evidence_index);
                let ties = if group_is_decoy[group] {
                    &mut best_decoy
                } else {
                    &mut best_target
                };
                match ties.first() {
                    Some(&first) if probability(group) == probability(first) => {
                        ties.insert(group);
                    }
                    Some(&first) if probability(group) >= probability(first) => {}
                    _ => {
                        ties.clear();
                        ties.insert(group);
                    }
                }
            }
            let chosen: BTreeSet<usize> = match (best_target.first(), best_decoy.first()) {
                (None, _) => best_decoy.clone(),
                (Some(_), None) => best_target.clone(),
                _ if targets_first => best_target.clone(),
                (Some(&t), Some(&d)) if probability(t) > probability(d) => best_target.clone(),
                (Some(&t), Some(&d)) if probability(d) > probability(t) => best_decoy.clone(),
                // Equal: the source's tie branch for targets is dead code
                // (its flag is always false), so all best groups are kept.
                _ => best_target.union(&best_decoy).copied().collect(),
            };
            let mut keep: BTreeSet<usize> = BTreeSet::new();
            if resolve_ties {
                if let Some(first) = chosen.first() {
                    keep = group_to_evidence.get(first).cloned().unwrap_or_default();
                    if chosen.len() > 1 {
                        let list = |set: &BTreeSet<usize>| {
                            set.iter()
                                .map(|g| format!("{g}={}, ", probability(*g)))
                                .collect::<String>()
                        };
                        messages.push(format!(
                            "Resolution: Peptide {} had groups:\ntgt: {}\ndec: {}\nKept: {first}",
                            best_hit.sequence,
                            list(&best_target),
                            list(&best_decoy)
                        ));
                    }
                }
            } else {
                for group in &chosen {
                    if let Some(evidence) = group_to_evidence.get(group) {
                        keep.extend(evidence.iter().copied());
                    }
                }
            }
            best_hit.evidences = keep
                .iter()
                .map(|&index| best_hit.evidences[index].clone())
                .collect();
        }
        peptides.clone_from_slice(&next);
        Ok(messages)
    }

    /// Build the graph, resolve it, then remove unreferenced proteins and
    /// update both kinds of protein groups (source static `run`). Only the
    /// first protein run is resolved.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when there is no protein run (the source
    /// reads out of bounds), or as [`build_graph`](Self::build_graph),
    /// [`resolve_graph`](Self::resolve_graph) and the ID filters. The input is
    /// unchanged on error.
    pub fn run(
        proteins: &mut [ProteinIdentification],
        peptides: &mut [PeptideIdentification],
    ) -> Result<Vec<String>> {
        if proteins.is_empty() {
            return Err(bad("no protein identification run to resolve"));
        }
        let mut next_proteins = proteins.to_vec();
        let mut next_peptides = peptides.to_vec();
        let mut resolver = Self::new(false);
        let mut messages = resolver.build_graph(&mut next_proteins[0], &next_peptides, false)?;
        messages.extend(resolver.resolve_graph(&mut next_proteins[0], &mut next_peptides)?);
        remove_unreferenced_proteins(&mut next_proteins, &next_peptides)?;
        let run = &mut next_proteins[0];
        update_protein_groups(&mut run.indistinguishable_groups, &run.hits)?;
        update_protein_groups(&mut run.protein_groups, &run.hits)?;
        proteins.clone_from_slice(&next_proteins);
        peptides.clone_from_slice(&next_peptides);
        Ok(messages)
    }

    /// Build and store the graph (source `buildGraph`). Correct resolution
    /// needs sorted groups, so the run's indistinguishable groups are sorted
    /// (by descending probability, then fewer accessions, then accessions)
    /// unless `skip_sort`, in which case nothing is modified.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when the run has no
    /// indistinguishable groups, or a best hit references an accession in no
    /// group (add such proteins as singletons). The resolver and the run are unchanged on
    /// error.
    pub fn build_graph(
        &mut self,
        protein: &mut ProteinIdentification,
        peptides: &[PeptideIdentification],
        skip_sort: bool,
    ) -> Result<Vec<String>> {
        if protein.indistinguishable_groups.is_empty() {
            return Err(no_groups());
        }
        let mut messages = vec![resolving_message(protein)];
        let mut groups = protein.indistinguishable_groups.clone();
        if !skip_sort {
            ProteinGroup::sort(&mut groups)?;
        }
        let mut next = self.clone();
        for (index, group) in groups.iter().enumerate() {
            for accession in &group.accessions {
                next.prot_acc_to_indist_prot_grp
                    .insert(accession.clone(), index);
            }
        }
        for (pep_index, pep) in peptides.iter().enumerate() {
            let Some(best_hit) = pep.hits.first() else {
                messages.push(
                    "Warning PeptideProteinResolution: Skipping spectrum without hits.".to_owned(),
                );
                continue;
            };
            for evidence in &best_hit.evidences {
                let accession = &evidence.protein_accession;
                let &group = next.prot_acc_to_indist_prot_grp.get(accession).ok_or_else(|| {
                    missing(format!("Not all proteins present in an indistinguishable group. ({accession} not found). Make sure to add them as singletons."))
                })?;
                next.pep_to_indist_prot_grp
                    .entry(pep_index)
                    .or_default()
                    .insert(group);
                next.indist_prot_grp_to_pep
                    .entry(group)
                    .or_default()
                    .insert(pep_index);
            }
        }
        protein.indistinguishable_groups = groups;
        *self = next;
        Ok(messages)
    }

    /// Resolve every connected component of the graph with
    /// [`resolve_connected_component`](Self::resolve_connected_component)
    /// (source `resolveGraph`). With `statistics`, the largest components
    /// are reported.
    ///
    /// All peptides are used, whatever protein run they belong to (a source
    /// to-do).
    ///
    /// # Errors
    ///
    /// As [`resolve_connected_component`](Self::resolve_connected_component).
    pub fn resolve_graph(
        &mut self,
        protein: &mut ProteinIdentification,
        peptides: &mut [PeptideIdentification],
    ) -> Result<Vec<String>> {
        let mut messages = Vec::new();
        let mut old_size = self.indist_prot_grp_to_pep.len();
        let (mut most_peps, mut most_grps, mut most_both) = (
            ConnectedComponent::default(),
            ConnectedComponent::default(),
            ConnectedComponent::default(),
        );
        while let Some(&root) = self.indist_prot_grp_to_pep.keys().next() {
            let current_size = self.indist_prot_grp_to_pep.len();
            if self.statistics && old_size - current_size > 1 {
                messages.push(format!(
                    "resolved group of size {} in last step ",
                    old_size - current_size
                ));
                old_size = current_size;
            }
            let component = self.find_connected_component(root);
            if self.statistics {
                if component.prot_grp_indices.len() > most_grps.prot_grp_indices.len() {
                    most_grps = component.clone();
                }
                if component.pep_indices.len() > most_peps.pep_indices.len() {
                    most_peps = component.clone();
                }
                if component.prot_grp_indices.len() + component.pep_indices.len()
                    > most_both.prot_grp_indices.len() + most_both.pep_indices.len()
                {
                    most_both = component.clone();
                }
                if component.prot_grp_indices.len() > 1 {
                    messages.push(format!("found group: \n{component}\nProcessing ..."));
                }
            }
            messages.extend(self.resolve_connected_component(&component, protein, peptides)?);
            for group in &component.prot_grp_indices {
                self.indist_prot_grp_to_pep.remove(group);
            }
            // A component always contains its root when the root has
            // peptides; guard against a root the search did not reach.
            self.indist_prot_grp_to_pep.remove(&root);
        }
        if self.statistics {
            messages.push(format!("Most protein groups in component:\n{most_grps}"));
            messages.push(format!("Most peptides in component:\n{most_peps}"));
            messages.push(format!("Biggest component:\n{most_both}"));
        }
        Ok(messages)
    }

    /// Breadth-first search over both sides of the graph from protein group
    /// `root_prot_grp`, alternating between groups and peptides (source
    /// `findConnectedComponent`). The root enters the component when a
    /// peptide leads back to it.
    pub fn find_connected_component(&self, root_prot_grp: usize) -> ConnectedComponent {
        let mut component = ConnectedComponent::default();
        let mut queue = VecDeque::from([(true, root_prot_grp)]);
        let empty = BTreeSet::new();
        while let Some((is_group, node)) = queue.pop_front() {
            let neighbors = if is_group {
                self.indist_prot_grp_to_pep.get(&node).unwrap_or(&empty)
            } else {
                self.pep_to_indist_prot_grp.get(&node).unwrap_or(&empty)
            };
            for &neighbor in neighbors {
                let inserted = if is_group {
                    component.pep_indices.insert(neighbor)
                } else {
                    component.prot_grp_indices.insert(neighbor)
                };
                if inserted {
                    queue.push_back((!is_group, neighbor));
                }
            }
        }
        component
    }

    /// Resolve one component and add it as a protein group (source
    /// `resolveConnectedComponent`); a component of one group is left alone.
    ///
    /// The component's groups are ordered by descending probability, then by
    /// descending peptide count; each in turn claims its remaining peptides,
    /// removing them from all later groups and removing every evidence of the
    /// peptide's best hit that does not belong to it. The added group lists
    /// the accessions of all component groups in that order, with the
    /// probability of the component's lowest-indexed group.
    ///
    /// The source sorts with `std::sort`, which does not fix the order of
    /// groups equal in probability and peptide count; this keeps them in index
    /// order (as libstdc++ does for components of at most 16 groups). The
    /// source's own comment on this method says it "does not produce correct
    /// results yet".
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a component group index is not a
    /// group of the run (the groups changed after the graph was built; the
    /// source logs a fatal error and reads out of bounds) or a peptide index
    /// is out of range. `protein` and `peptides` are unchanged on error.
    pub fn resolve_connected_component(
        &mut self,
        component: &ConnectedComponent,
        protein: &mut ProteinIdentification,
        peptides: &mut [PeptideIdentification],
    ) -> Result<Vec<String>> {
        let mut messages = Vec::new();
        if component.prot_grp_indices.len() <= 1 {
            return Ok(messages);
        }
        let groups = &protein.indistinguishable_groups;
        for &group in &component.prot_grp_indices {
            if group >= groups.len() {
                return Err(bad(format!(
                    "Something went terribly wrong. Group with index {group}doesn't exist.  ProteinPeptideResolution: Groups changed after building data structures."
                )));
            }
        }
        let Some(&best) = component.prot_grp_indices.first() else {
            return Ok(messages);
        };
        let mut ambiguity = ProteinGroup {
            probability: groups[best].probability,
            ..ProteinGroup::default()
        };
        let count = |g: usize| self.indist_prot_grp_to_pep.get(&g).map_or(0, BTreeSet::len);
        let mut order: Vec<usize> = component.prot_grp_indices.iter().copied().collect();
        order.sort_by(|&a, &b| {
            (groups[b].probability, count(b))
                .partial_cmp(&(groups[a].probability, count(a)))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mut next_graph = self.indist_prot_grp_to_pep.clone();
        let mut next_peptides = peptides.to_vec();
        for (position, &group) in order.iter().enumerate() {
            let accessions = &groups[group].accessions;
            ambiguity.accessions.extend(accessions.iter().cloned());
            let claimed = next_graph.get(&group).cloned().unwrap_or_default();
            if self.statistics {
                messages.push(format!(
                    "Group: {} steals {} peptides for itself.",
                    accessions
                        .iter()
                        .map(|a| format!("{a}, "))
                        .collect::<String>(),
                    claimed.len()
                ));
            }
            for &pep in &claimed {
                for later in &order[position + 1..] {
                    if let Some(set) = next_graph.get_mut(later) {
                        set.remove(&pep);
                    }
                }
                let hit = next_peptides
                    .get_mut(pep)
                    .and_then(|p| p.hits.first_mut())
                    .ok_or_else(|| bad(format!("peptide index {pep} has no best hit")))?;
                hit.evidences
                    .retain(|e| accessions.contains(&e.protein_accession));
            }
        }
        protein.protein_groups.push(ambiguity);
        self.indist_prot_grp_to_pep = next_graph;
        peptides.clone_from_slice(&next_peptides);
        Ok(messages)
    }
}

fn no_groups() -> Error {
    missing(
        "No indistinguishable Groups annotated. Currently this class only resolves across groups.",
    )
}

fn resolving_message(protein: &ProteinIdentification) -> String {
    format!(
        "Resolving peptides between {} proteins in {} indistinguishable groups.",
        protein.hits.len(),
        protein.indistinguishable_groups.len()
    )
}
