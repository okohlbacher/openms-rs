// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Timo Sachsenberg, OpenMS Rust contributors $
//! Fixed and variable RNA modifications from pinned ModifiedNASequenceGenerator.
//! Candidates retain caller order; repeated Arc identities are removed. The
//! source maximum-one path places terminal-specific records in residue slots.
//! See docs/RNA_MODIFICATION_SUPPORT.md for ordering and serialization limits.

use super::{NASequence, Ribonucleotide, RibonucleotideTermSpecificity};
use crate::{Error, Result};
use std::sync::Arc;

type Modification = Arc<Ribonucleotide>;
const MAX_SCRATCH_BYTES: usize = 64 * 1024 * 1024;

/// Checked enumeration with atomic fixed replacement and variable append.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ModifiedNASequenceGenerator {
    pub max_residues: usize,
    /// Maximum supplied candidates and total compatible (site, record) entries.
    pub max_sites: usize,
    /// Cumulative scan, pointer-comparison, planning, validation and copy work.
    pub max_work: usize,
    /// Includes existing entries when appending; duplicate values count separately.
    pub max_outputs: usize,
    /// Conservative logical payload, including record strings/formulas and the
    /// source sequence's owned sulfur-slicing context. This is not an RSS limit.
    pub max_output_bytes: usize,
}
impl Default for ModifiedNASequenceGenerator {
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

#[derive(Clone, Copy)]
enum Site {
    ThreePrime,
    FivePrime,
    Residue(usize),
}
struct CompatibleSite<'a> {
    site: Site,
    alternatives: Vec<&'a Modification>,
}

impl ModifiedNASequenceGenerator {
    /// Install complete records. Existing modifications are protected; the first
    /// candidate wins at an empty terminus, the last at an eligible residue.
    pub fn apply_fixed_modifications(
        &self,
        modifications: &[Modification],
        sequence: &mut NASequence,
    ) -> Result<()> {
        let mut work = self.preflight(modifications, sequence)?;
        let base_bytes = sequence.generation_payload_bytes()?;
        let unique = unique(modifications, &mut work)?;
        let sites = self.compatibility(&unique, sequence, false, &mut work)?;
        let mut assignments = Vec::new();
        work.reserve(&mut assignments, sites.len())?;
        let mut bytes = base_bytes;
        for site in &sites {
            work.charge(1)?;
            let modification = match site.site {
                Site::Residue(_) => site.alternatives.last().expect("nonempty site"),
                _ => site.alternatives.first().expect("nonempty site"),
            };
            bytes = add(bytes, modification.payload_bytes()?)?;
            assignments.push((site.site, *modification));
        }
        self.outputs(0, 0, 1, bytes)?;
        work.charge(copy_work(sequence.len(), assignments.len())?)?;
        let next = apply(sequence, &assignments)?;
        *sequence = next;
        Ok(())
    }

    /// Return all variants, with the optional input first, then increasing
    /// placement count and the source's reverse subset order.
    pub fn variable_modifications(
        &self,
        modifications: &[Modification],
        sequence: &NASequence,
        maximum: usize,
        keep_original: bool,
    ) -> Result<Vec<NASequence>> {
        let mut output = Vec::new();
        self.apply_variable_modifications(
            modifications,
            sequence,
            maximum,
            &mut output,
            keep_original,
        )?;
        Ok(output)
    }

    /// Append after complete preflight and successful staging. A failure leaves
    /// all existing entries unchanged, including failures in a later variant.
    pub fn apply_variable_modifications(
        &self,
        modifications: &[Modification],
        sequence: &NASequence,
        maximum: usize,
        output: &mut Vec<NASequence>,
        keep_original: bool,
    ) -> Result<()> {
        let mut work = self.preflight(modifications, sequence)?;
        self.outputs(output.len(), 0, 0, 0)?;
        let mut existing_bytes = 0;
        for existing in output.iter() {
            work.charge(add(existing.len(), 4)?)?;
            existing_bytes = add(existing_bytes, existing.generation_payload_bytes()?)?;
            self.outputs(output.len(), existing_bytes, 0, 0)?;
        }
        let base_bytes = sequence.generation_payload_bytes()?;
        // The source returns before inspecting compatibility for these noops.
        let sites = if maximum == 0 || modifications.is_empty() {
            Vec::new()
        } else {
            let unique = unique(modifications, &mut work)?;
            self.compatibility(&unique, sequence, maximum == 1, &mut work)?
        };
        let maximum = maximum.min(sites.len());
        let (generated, bytes) = self.count_outputs(
            &sites,
            maximum,
            keep_original,
            base_bytes,
            output.len(),
            existing_bytes,
            &mut work,
        )?;
        self.outputs(output.len(), existing_bytes, generated, bytes)?;
        // Includes per-variant validation, copying, choices, subset transitions
        // and assignment collection. No clone starts before this complete floor.
        work.charge(mul(generated, copy_work(sequence.len(), maximum)?)?)?;
        let mut staged = Vec::new();
        staged
            .try_reserve_exact(generated)
            .map_err(|_| invalid("RNA variant allocation failed"))?;
        if keep_original {
            staged.push(apply(sequence, &[])?);
        }
        let mut subset = Vec::new();
        let mut choices = Vec::new();
        let mut assignments = Vec::new();
        work.reserve(&mut subset, maximum)?;
        work.reserve(&mut choices, maximum)?;
        work.reserve(&mut assignments, maximum)?;
        for count in 1..=maximum {
            subset.clear();
            subset.extend(sites.len() - count..sites.len());
            choices.resize(count, 0);
            loop {
                choices.fill(0);
                loop {
                    assignments.clear();
                    for (&site, &choice) in subset.iter().zip(&choices) {
                        assignments.push((sites[site].site, sites[site].alternatives[choice]));
                    }
                    staged.push(apply(sequence, &assignments)?);
                    // Source recursion visits ascending sites, so the last
                    // selected site's alternative is the fastest-changing digit.
                    let mut digit = count;
                    loop {
                        if digit == 0 {
                            break;
                        }
                        digit -= 1;
                        choices[digit] += 1;
                        if choices[digit] < sites[subset[digit]].alternatives.len() {
                            break;
                        }
                        choices[digit] = 0;
                    }
                    if digit == 0 && choices[0] == 0 {
                        break;
                    }
                }
                if !previous_subset(&mut subset, sites.len()) {
                    break;
                }
            }
        }
        debug_assert_eq!(staged.len(), generated);
        output
            .try_reserve(generated)
            .map_err(|_| invalid("RNA append allocation failed"))?;
        output.append(&mut staged);
        Ok(())
    }

    fn preflight(&self, modifications: &[Modification], sequence: &NASequence) -> Result<Work> {
        if self.max_residues == 0 || self.max_sites == 0 || self.max_work == 0 {
            return Err(invalid(
                "RNA generation residue/site/work limits must be positive",
            ));
        }
        if sequence.len() > self.max_residues || modifications.len() > self.max_sites {
            return Err(invalid(
                "RNA generation residue or candidate limit exceeded",
            ));
        }
        let mut work = Work {
            remaining: self.max_work,
            scratch: MAX_SCRATCH_BYTES,
        };
        work.charge(add(add(sequence.len(), modifications.len())?, 4)?)?;
        Ok(work)
    }

    fn compatibility<'a>(
        &self,
        modifications: &[&'a Modification],
        sequence: &NASequence,
        maximum_one: bool,
        work: &mut Work,
    ) -> Result<Vec<CompatibleSite<'a>>> {
        let mut result = Vec::new();
        work.reserve(&mut result, add(sequence.len(), 2)?.min(self.max_sites))?;
        let mut entries = 0;
        if !maximum_one {
            for (site, specificity, occupied) in [
                (
                    Site::ThreePrime,
                    RibonucleotideTermSpecificity::ThreePrime,
                    sequence.three_prime_mod().is_some(),
                ),
                (
                    Site::FivePrime,
                    RibonucleotideTermSpecificity::FivePrime,
                    sequence.five_prime_mod().is_some(),
                ),
            ] {
                if occupied {
                    continue;
                }
                let mut alternatives = Vec::new();
                for &modification in modifications {
                    work.charge(1)?;
                    if modification.term_specificity() == specificity {
                        self.alternative(
                            &mut alternatives,
                            modification,
                            modifications.len(),
                            &mut entries,
                            work,
                        )?;
                    }
                }
                if !alternatives.is_empty() {
                    result.push(CompatibleSite { site, alternatives });
                }
            }
        }
        for (index, residue) in sequence.residues().iter().enumerate() {
            work.charge(1)?;
            if residue.is_modified() {
                continue;
            }
            let mut alternatives = Vec::new();
            for &modification in modifications {
                work.charge(1)?;
                if residue.code().len() == 1
                    && char::from(residue.code().as_bytes()[0]) == modification.origin()
                    && (maximum_one
                        || modification.term_specificity()
                            == RibonucleotideTermSpecificity::Anywhere)
                {
                    self.alternative(
                        &mut alternatives,
                        modification,
                        modifications.len(),
                        &mut entries,
                        work,
                    )?;
                }
            }
            if !alternatives.is_empty() {
                result.push(CompatibleSite {
                    site: Site::Residue(index),
                    alternatives,
                });
            }
        }
        Ok(result)
    }

    fn alternative<'a>(
        &self,
        alternatives: &mut Vec<&'a Modification>,
        modification: &'a Modification,
        capacity: usize,
        entries: &mut usize,
        work: &mut Work,
    ) -> Result<()> {
        *entries = add(*entries, 1)?;
        if *entries > self.max_sites {
            return Err(invalid("RNA compatible placement limit exceeded"));
        }
        if alternatives.is_empty() {
            work.reserve(alternatives, capacity)?;
        }
        alternatives.push(modification);
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn count_outputs(
        &self,
        sites: &[CompatibleSite<'_>],
        maximum: usize,
        keep_original: bool,
        base_bytes: usize,
        existing: usize,
        existing_bytes: usize,
        work: &mut Work,
    ) -> Result<(usize, usize)> {
        let mut counts = Vec::new();
        let mut extras = Vec::new();
        let size = add(maximum, 1)?;
        work.reserve(&mut counts, size)?;
        work.reserve(&mut extras, size)?;
        counts.resize(size, 0usize);
        extras.resize(size, 0usize);
        counts[0] = 1;
        let mut generated = usize::from(keep_original);
        let mut bytes = mul(generated, base_bytes)?;
        self.outputs(existing, existing_bytes, generated, bytes)?;
        for (index, site) in sites.iter().enumerate() {
            let mut payload_sum = 0;
            for modification in &site.alternatives {
                work.charge(1)?;
                payload_sum = add(payload_sum, modification.payload_bytes()?)?;
            }
            for count in (1..=maximum.min(index + 1)).rev() {
                work.charge(1)?;
                let additional = mul(counts[count - 1], site.alternatives.len())?;
                let extra = add(
                    mul(extras[count - 1], site.alternatives.len())?,
                    mul(counts[count - 1], payload_sum)?,
                )?;
                counts[count] = add(counts[count], additional)?;
                extras[count] = add(extras[count], extra)?;
                generated = add(generated, additional)?;
                bytes = add(bytes, add(mul(additional, base_bytes)?, extra)?)?;
            }
            self.outputs(existing, existing_bytes, generated, bytes)?;
        }
        Ok((generated, bytes))
    }

    fn outputs(
        &self,
        existing: usize,
        existing_bytes: usize,
        generated: usize,
        bytes: usize,
    ) -> Result<()> {
        if add(existing, generated)? > self.max_outputs {
            return Err(invalid("RNA generation output count limit exceeded"));
        }
        if add(existing_bytes, bytes)? > self.max_output_bytes {
            return Err(invalid("RNA generation output payload limit exceeded"));
        }
        Ok(())
    }
}

fn unique<'a>(modifications: &'a [Modification], work: &mut Work) -> Result<Vec<&'a Modification>> {
    // ponytail: at most 1024 default candidates; a bounded pointer scan preserves
    // caller order without address sorting or a second indexing collection.
    work.charge(mul(modifications.len(), modifications.len().saturating_sub(1))? / 2)?;
    let mut result: Vec<&Modification> = Vec::new();
    work.reserve(&mut result, modifications.len())?;
    for modification in modifications {
        if !result
            .iter()
            .any(|previous| Arc::ptr_eq(previous, modification))
        {
            result.push(modification);
        }
    }
    Ok(result)
}

/// Reverse lexicographic order: for four sites, pairs cd, bd, bc, ad, ac, ab.
fn previous_subset(subset: &mut [usize], sites: usize) -> bool {
    for index in (0..subset.len()).rev() {
        let minimum = if index == 0 { 0 } else { subset[index - 1] + 1 };
        if subset[index] > minimum {
            subset[index] -= 1;
            for following in index + 1..subset.len() {
                subset[following] = sites - subset.len() + following;
            }
            return true;
        }
    }
    false
}

fn apply(sequence: &NASequence, assignments: &[(Site, &Modification)]) -> Result<NASequence> {
    let mut residues = Vec::new();
    residues
        .try_reserve_exact(sequence.len())
        .map_err(|_| invalid("RNA residue copy allocation failed"))?;
    residues.extend_from_slice(sequence.residues());
    let mut five = sequence.five_prime_mod().cloned();
    let mut three = sequence.three_prime_mod().cloned();
    for &(site, modification) in assignments {
        match site {
            Site::ThreePrime => three = Some(Arc::clone(modification)),
            Site::FivePrime => five = Some(Arc::clone(modification)),
            Site::Residue(index) => residues[index] = Arc::clone(modification),
        }
    }
    sequence.with_generation_records(residues, five, three)
}
fn copy_work(residues: usize, maximum: usize) -> Result<usize> {
    add(add(mul(residues, 2)?, mul(maximum, 6)?)?, 20)
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b)
        .ok_or_else(|| invalid("RNA generation size overflows"))
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b)
        .ok_or_else(|| invalid("RNA generation size overflows"))
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
struct Work {
    remaining: usize,
    scratch: usize,
}
impl Work {
    fn charge(&mut self, amount: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(amount)
            .ok_or_else(|| invalid("RNA generation work limit exceeded"))?;
        Ok(())
    }
    fn reserve<T>(&mut self, output: &mut Vec<T>, capacity: usize) -> Result<()> {
        self.charge(capacity)?;
        self.scratch = self
            .scratch
            .checked_sub(mul(capacity, std::mem::size_of::<T>())?)
            .ok_or_else(|| invalid("RNA generation scratch allocation limit exceeded"))?;
        output
            .try_reserve_exact(capacity)
            .map_err(|_| invalid("RNA generation scratch allocation failed"))
    }
}
