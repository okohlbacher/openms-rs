// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Timo Sachsenberg, OpenMS Rust contributors $

//! Fixed and variable modification placement from pinned ModifiedPeptideGenerator.
//! The source's maximum-one terminal placements differ from its general path;
//! see docs/MODIFIED_PEPTIDES_SUPPORT.md before serializing those unusual states.

use super::sequence::GENERATION_FORMULA_ENTRY_BYTES;
use super::{AASequence, ModificationsDB, ResidueModification, TermSpecificity};
use crate::{Error, Result};
use std::{collections::BTreeMap, sync::Arc};

type Modification = Arc<ResidueModification>;

/// Bounded modification generation. Existing output is included in append limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ModifiedPeptideGenerator {
    pub max_residues: usize,
    /// Maximum supplied modifications and total compatible placement entries,
    /// including repeated terminal alternatives present in the source algorithm.
    pub max_sites: usize,
    /// Scan, comparison, planning and sequence-rebuild allowances; not CPU time.
    pub max_work: usize,
    pub max_outputs: usize,
    /// Conservative owned-payload allowance, including anonymous-tag strings;
    /// documented in MODIFIED_PEPTIDES_SUPPORT.md, not an allocator/RSS limit.
    pub max_output_bytes: usize,
}

impl Default for ModifiedPeptideGenerator {
    fn default() -> Self {
        Self {
            max_residues: 10_000,
            max_sites: 1_024,
            max_work: 10_000_000,
            max_outputs: 100_000,
            max_output_bytes: 256_000_000,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Site {
    CTerm,
    NTerm,
    Residue(usize),
}

impl ModifiedPeptideGenerator {
    /// Resolve names once against the immutable global database. Prefer full IDs;
    /// ambiguous names are errors. Duplicate full IDs collapse to one alternative.
    pub fn get_modifications(names: &[&str]) -> Result<Vec<Modification>> {
        Self::get_modifications_with_registry(names, ModificationsDB::global())
    }
    /// Resolve shared chemistry from a caller registry. Returned handles and
    /// generated peptides remain valid after that registry is dropped.
    pub fn get_modifications_with_registry(
        names: &[&str],
        db: &ModificationsDB,
    ) -> Result<Vec<Modification>> {
        let mut modifications = names
            .iter()
            .map(|name| db.get_modification_handle(name, None, None))
            .collect::<Result<Vec<_>>>()?;
        order_modifications(&mut modifications);
        Ok(modifications)
    }

    /// Apply fixed modifications in two source passes, committing only on success.
    /// Existing residue annotations are protected. A matching boundary-specific
    /// record can overwrite an existing terminal annotation, as in the source.
    pub fn apply_fixed_modifications(
        &self,
        modifications: &[Modification],
        peptide: &mut AASequence,
    ) -> Result<()> {
        let mut work = self.preflight(modifications, peptide)?;
        let modifications = sorted(modifications, &mut work)?;
        let mut assignments = BTreeMap::new();
        let mut entries = 0;
        for modification in &modifications {
            work.charge(1)?;
            let site = match modification.term_specificity() {
                TermSpecificity::NTerm if peptide.n_terminal_modification().is_none() => {
                    Some(Site::NTerm)
                }
                TermSpecificity::CTerm if peptide.c_terminal_modification().is_none() => {
                    Some(Site::CTerm)
                }
                _ => None,
            };
            if let Some(site) = site {
                // The first terminal record wins this pass. A later origin match
                // in the residue pass can still replace it.
                if let std::collections::btree_map::Entry::Vacant(slot) = assignments.entry(site) {
                    self.entry(&mut entries)?;
                    slot.insert(Arc::clone(modification));
                }
            }
        }
        for (index, residue) in peptide.as_str().bytes().enumerate() {
            if peptide.residue_modification(index)?.is_some() {
                continue;
            }
            for modification in &modifications {
                work.charge(1)?;
                if origin_matches(modification, residue) {
                    if let Some(site) = compatible(modification, index, peptide.len()) {
                        self.entry(&mut entries)?;
                        // Source does not stop after installing one fixed record.
                        assignments.insert(site, Arc::clone(modification));
                    }
                }
            }
        }
        let extra = assignments.values().try_fold(0usize, |n, modification| {
            add(n, formula_bytes(modification)?)
        })?;
        self.outputs(0, 0, 1, add(peptide.generation_payload_bytes()?, extra)?)?;
        work.charge(add(peptide.len(), assignments.len())?)?;
        let next = apply(peptide, assignments.into_iter())?;
        *peptide = next;
        Ok(())
    }

    /// Return variable variants in source reverse-site/subset order.
    pub fn variable_modifications(
        &self,
        modifications: &[Modification],
        peptide: &AASequence,
        maximum: usize,
        keep_unmodified: bool,
    ) -> Result<Vec<AASequence>> {
        let mut output = Vec::new();
        self.apply_variable_modifications(
            modifications,
            peptide,
            maximum,
            &mut output,
            keep_unmodified,
        )?;
        Ok(output)
    }

    /// Append atomically. Count and payload bounds are checked before cloning any
    /// combinatorial sequence; invalid chemistry also leaves all output unchanged.
    pub fn apply_variable_modifications(
        &self,
        modifications: &[Modification],
        peptide: &AASequence,
        maximum: usize,
        output: &mut Vec<AASequence>,
        keep_unmodified: bool,
    ) -> Result<()> {
        let mut work = self.preflight(modifications, peptide)?;
        if output.len() > self.max_outputs {
            return Err(bad("existing peptides exceed output count limit"));
        }
        let existing_bytes = output.iter().try_fold(0usize, |n, sequence| {
            work.charge(add(sequence.len(), 1)?)?;
            let bytes = add(n, sequence.generation_payload_bytes()?)?;
            if bytes > self.max_output_bytes {
                return Err(bad("existing peptides exceed output payload limit"));
            }
            Ok(bytes)
        })?;
        let modifications = sorted(modifications, &mut work)?;
        let mut compatibility: BTreeMap<Site, Vec<Modification>> = BTreeMap::new();
        let mut entries = 0;
        if maximum > 1 {
            for modification in &modifications {
                work.charge(1)?;
                let site = match modification.term_specificity() {
                    TermSpecificity::NTerm if peptide.n_terminal_modification().is_none() => {
                        Some(Site::NTerm)
                    }
                    TermSpecificity::CTerm if peptide.c_terminal_modification().is_none() => {
                        Some(Site::CTerm)
                    }
                    _ => None,
                };
                if let Some(site) = site {
                    self.entry(&mut entries)?;
                    compatibility
                        .entry(site)
                        .or_default()
                        .push(Arc::clone(modification));
                }
            }
        }
        if maximum > 0 {
            for (index, residue) in peptide.as_str().bytes().enumerate() {
                if peptide.residue_modification(index)?.is_some() {
                    continue;
                }
                for modification in &modifications {
                    work.charge(1)?;
                    if !origin_matches(modification, residue) {
                        continue;
                    }
                    if let Some(mut site) = compatible(modification, index, peptide.len()) {
                        if maximum == 1 {
                            if modification.origin().unwrap_or('X') == 'X'
                                && matches!(
                                    modification.term_specificity(),
                                    TermSpecificity::NTerm | TermSpecificity::CTerm
                                )
                            {
                                return Err(bad(
                                    "source maximum-one placement would use a null residue pointer",
                                ));
                            }
                            site = Site::Residue(index);
                        }
                        self.entry(&mut entries)?;
                        compatibility
                            .entry(site)
                            .or_default()
                            .push(Arc::clone(modification));
                    }
                }
            }
        }
        let sites: Vec<_> = compatibility.into_iter().rev().collect();
        let maximum = maximum.min(sites.len());
        // Count weighted subsets first. Counts include repeated source terminal
        // alternatives; outputs are deliberately not deduplicated.
        work.charge(add(maximum, 1)?)?;
        let mut counts = vec![0usize; maximum + 1];
        counts[0] = 1;
        let mut generated = usize::from(keep_unmodified);
        for (index, (_, alternatives)) in sites.iter().enumerate() {
            for depth in (1..=maximum.min(index + 1)).rev() {
                work.charge(1)?;
                let additional = mul(counts[depth - 1], alternatives.len())?;
                counts[depth] = add(counts[depth], additional)?;
                generated = add(generated, additional)?;
            }
            if add(output.len(), generated)? > self.max_outputs {
                return Err(bad("generated peptides exceed output count limit"));
            }
        }
        let extra = sites
            .iter()
            .flat_map(|(_, alternatives)| alternatives)
            .try_fold(0usize, |largest, modification| {
                work.charge(1)?;
                Ok::<_, Error>(largest.max(formula_bytes(modification)?))
            })?;
        let payload = add(peptide.generation_payload_bytes()?, mul(extra, maximum)?)?;
        self.outputs(output.len(), existing_bytes, generated, payload)?;
        // Charge the complete sequence-rebuild floor before allocating variants.
        work.charge(mul(generated, add(add(peptide.len(), maximum)?, 1)?)?)?;
        // Source subset groups are appended in increasing bit-mask order over
        // reverse sites. Store only site indices, never intermediate AASequences.
        let mut groups = vec![Vec::<usize>::new()];
        for site in 0..sites.len() {
            let old_len = groups.len();
            for index in 0..old_len {
                work.charge(1)?;
                if groups[index].len() < maximum {
                    work.charge(add(groups[index].len(), 1)?)?;
                    let mut group = groups[index].clone();
                    group.push(site);
                    groups.push(group);
                }
            }
        }
        let mut staged = Vec::with_capacity(generated);
        for group in groups {
            if group.is_empty() && !keep_unmodified {
                continue;
            }
            let mut choices = vec![0usize; group.len()];
            loop {
                staged.push(apply(
                    peptide,
                    group.iter().zip(&choices).map(|(&site, &choice)| {
                        (sites[site].0, Arc::clone(&sites[site].1[choice]))
                    }),
                )?);
                let mut digit = 0;
                while digit < choices.len() {
                    choices[digit] += 1;
                    if choices[digit] < sites[group[digit]].1.len() {
                        break;
                    }
                    choices[digit] = 0;
                    digit += 1;
                }
                if digit == choices.len() {
                    break;
                }
            }
        }
        debug_assert_eq!(staged.len(), generated);
        output.extend(staged);
        Ok(())
    }

    fn preflight(&self, modifications: &[Modification], peptide: &AASequence) -> Result<Work> {
        if self.max_residues == 0 || self.max_sites == 0 || self.max_work == 0 {
            return Err(bad("generation point/site/work limits must be positive"));
        }
        if peptide.len() > self.max_residues || modifications.len() > self.max_sites {
            return Err(bad(
                "peptide length or modification count exceeds generation limit",
            ));
        }
        let mut work = Work(self.max_work);
        work.charge(add(add(peptide.len(), modifications.len())?, 1)?)?;
        Ok(work)
    }
    fn entry(&self, entries: &mut usize) -> Result<()> {
        *entries = add(*entries, 1)?;
        if *entries > self.max_sites {
            return Err(bad("compatible modification placements exceed site limit"));
        }
        Ok(())
    }
    fn outputs(
        &self,
        existing: usize,
        bytes: usize,
        generated: usize,
        payload: usize,
    ) -> Result<()> {
        if add(existing, generated)? > self.max_outputs {
            return Err(bad("generated peptides exceed output count limit"));
        }
        if add(bytes, mul(generated, payload)?)? > self.max_output_bytes {
            return Err(bad("generated peptides exceed output payload limit"));
        }
        Ok(())
    }
}

fn order_modifications(modifications: &mut Vec<Modification>) {
    modifications.sort_by(|a, b| a.full_id().cmp(b.full_id()));
    modifications.dedup_by(|a, b| a.full_id() == b.full_id());
}
fn sorted(modifications: &[Modification], work: &mut Work) -> Result<Vec<Modification>> {
    let logarithm = if modifications.len() < 2 {
        1
    } else {
        usize::BITS as usize - (modifications.len() - 1).leading_zeros() as usize
    };
    work.charge(mul(modifications.len(), logarithm)?)?;
    let mut result = modifications.to_vec();
    order_modifications(&mut result);
    Ok(result)
}
fn origin_matches(modification: &ResidueModification, residue: u8) -> bool {
    modification.origin().unwrap_or('X') == char::from(residue)
}
fn compatible(modification: &ResidueModification, index: usize, length: usize) -> Option<Site> {
    match modification.term_specificity() {
        TermSpecificity::Anywhere => Some(Site::Residue(index)),
        TermSpecificity::NTerm if index == 0 => Some(Site::NTerm),
        TermSpecificity::CTerm if index + 1 == length => Some(Site::CTerm),
        _ => None,
    }
}
fn apply(
    peptide: &AASequence,
    placements: impl Iterator<Item = (Site, Modification)>,
) -> Result<AASequence> {
    let mut residues = Vec::new();
    let (mut n, mut c) = (None, None);
    for (site, modification) in placements {
        match site {
            Site::Residue(index) => residues.push((index, modification)),
            Site::NTerm => n = Some(modification),
            Site::CTerm => c = Some(modification),
        }
    }
    peptide.with_resolved_modifications(&residues, n, c)
}
fn formula_bytes(modification: &ResidueModification) -> Result<usize> {
    mul(
        add(
            modification.diff_formula().atoms.len(),
            modification.absolute_formula().map_or(0, |f| f.atoms.len()),
        )?,
        GENERATION_FORMULA_ENTRY_BYTES,
    )
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b)
        .ok_or_else(|| bad("generation size overflows"))
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b)
        .ok_or_else(|| bad("generation size overflows"))
}
fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
struct Work(usize);
impl Work {
    fn charge(&mut self, amount: usize) -> Result<()> {
        self.0 = self
            .0
            .checked_sub(amount)
            .ok_or_else(|| bad("modification generation exceeds work limit"))?;
        Ok(())
    }
}
