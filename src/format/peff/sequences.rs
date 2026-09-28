// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `PEFFEntry` sequence materialisation: `getSequence`,
//! `getModifiedSequence`, `getVariantSequences`, `getProcessedSequence`,
//! `digestWithVariants`, `generatePeptides` and the private
//! `enumeratePEFFModifications_` (`PEFFFile.cpp:528-1300`).
//!
//! Modifications resolve against a [`ModificationsDB`]: the plain methods use
//! [`ModificationsDB::global`], the `_with_registry` forms a caller registry
//! (for example the global one extended with PSI-MOD through
//! [`ModificationsDB::extend_obo`], which the source loads by default).

use super::{PEFFEntry, PEFFModification, PEFFVariantSimple};
use crate::chemistry::{
    AASequence, DigestionSpecificity, ModificationsDB, ModifiedPeptideGenerator, Protease,
    ProteaseDigestion, ResidueModification, TermSpecificity,
};
use crate::{Error, Result};
use std::collections::BTreeSet;

/// The source's cap on the variants plus modifications combined within one
/// peptide (`MAX_COMBINATORIAL_ELEMENTS`, `MAX_MODS`): 2^20 combinations.
pub const MAX_COMBINATORIAL_ELEMENTS: usize = 20;

/// The PEFF CV accession of a signal peptide, the default region of
/// [`PEFFEntry::processed_sequence`].
pub const SIGNAL_PEPTIDE: &str = "PEFF:0001021";

/// Peptides generated from a [`PEFFEntry`]: the source's two output vectors,
/// kept index-aligned, and its `OPENMS_LOG_WARN` lines.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PEFFPeptides {
    /// One description per peptide; empty for a reference peptide.
    pub descriptions: Vec<String>,
    /// The peptides, index-aligned with `descriptions`.
    pub sequences: Vec<AASequence>,
    /// Warnings, verbatim and in order, without line ends.
    pub warnings: Vec<String>,
}

impl PEFFPeptides {
    /// The number of peptides.
    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    /// True when no peptide was generated.
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    fn push(&mut self, description: String, sequence: AASequence, limit: usize) -> Result<()> {
        if self.sequences.len() >= limit {
            return Err(Error::InvalidValue(format!(
                "PEFF peptide generation exceeds the output limit of {limit} peptides"
            )));
        }
        self.descriptions.push(description);
        self.sequences.push(sequence);
        Ok(())
    }
}

/// Options of [`PEFFEntry::digest_with_variants`]; the defaults are the
/// source's default arguments.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct VariantDigestOptions {
    /// Minimum peptide length (source default 6).
    pub min_length: usize,
    /// Maximum peptide length; 0 means no limit (source default 40).
    pub max_length: usize,
    /// Include the reference peptides, with empty descriptions.
    pub include_reference: bool,
    /// Enumerate simple variants.
    pub include_variants: bool,
    /// Enumerate the annotated modifications.
    pub include_modifications: bool,
    /// Output ceiling; exceeding it is an error, never a truncation. Not in
    /// the source.
    pub max_peptides: usize,
}

impl Default for VariantDigestOptions {
    fn default() -> Self {
        Self {
            min_length: 6,
            max_length: 40,
            include_reference: true,
            include_variants: true,
            include_modifications: false,
            max_peptides: 1 << 24,
        }
    }
}

/// Options of [`PEFFEntry::generate_peptides`]; the defaults are the
/// source's default arguments.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PeptideGenerationOptions {
    /// Fixed modifications, for example `Carbamidomethyl (C)`.
    pub fixed_mods: Vec<String>,
    /// Variable modifications, for example `Oxidation (M)`.
    pub variable_mods: Vec<String>,
    /// Maximum variable modifications per peptide (source default 2).
    pub max_variable_mods_per_peptide: usize,
    /// Minimum peptide length (source default 6).
    pub min_length: usize,
    /// Maximum peptide length; 0 means no limit (source default 40).
    pub max_length: usize,
    /// Include the reference peptides.
    pub include_reference: bool,
    /// Enumerate the simple variants.
    pub include_peff_variants: bool,
    /// Enumerate the annotated modifications.
    pub include_peff_modifications: bool,
    /// Maximum annotated modifications applied together per peptide; 0 is
    /// unlimited. With n sites and cap k the output holds the subsets of size
    /// at most k, sum over j = 0..k of C(n, j), instead of 2^n.
    pub max_peff_mods_per_proteoform: usize,
    /// Output ceiling; exceeding it is an error, never a truncation. Not in
    /// the source.
    pub max_peptides: usize,
}

impl Default for PeptideGenerationOptions {
    fn default() -> Self {
        Self {
            fixed_mods: Vec::new(),
            variable_mods: Vec::new(),
            max_variable_mods_per_peptide: 2,
            min_length: 6,
            max_length: 40,
            include_reference: true,
            include_peff_variants: true,
            include_peff_modifications: true,
            max_peff_mods_per_proteoform: 0,
            max_peptides: 1 << 24,
        }
    }
}

/// Source `ModificationsDB::getModification(name, residue, ANYWHERE)`: the
/// first residue-compatible record of anywhere specificity. The source warns
/// about, and then takes the first of, several matches; so does this.
fn lookup<'a>(
    db: &'a ModificationsDB,
    name: &str,
    residue: char,
) -> Option<&'a ResidueModification> {
    db.find(name, Some(residue), Some(TermSpecificity::Anywhere))
        .first()
        .copied()
}

/// The accession first, then the name (source inline blocks of
/// `getModifiedSequence`, `digestWithVariants`, `enumeratePEFFModifications_`).
fn resolve<'a>(
    db: &'a ModificationsDB,
    modification: &PEFFModification,
    residue: char,
) -> Option<&'a ResidueModification> {
    let mut found = None;
    if !modification.accession.is_empty() {
        found = lookup(db, &modification.accession, residue);
    }
    if found.is_none() && !modification.name.is_empty() {
        found = lookup(db, &modification.name, residue);
    }
    found
}

/// Resolve and install one modification at a zero-based residue index;
/// false when it cannot be resolved (or installed), leaving `peptide` as is.
fn apply(
    db: &ModificationsDB,
    peptide: &mut AASequence,
    index: usize,
    modification: &PEFFModification,
) -> bool {
    let Some(residue) = peptide
        .as_str()
        .as_bytes()
        .get(index)
        .map(|&b| char::from(b))
    else {
        return false;
    };
    match resolve(db, modification, residue) {
        Some(record) => peptide
            .set_modification_with_registry(index, record.full_id(), db)
            .is_ok(),
        None => false,
    }
}

/// `text` with the byte at `index` replaced by `residue`; `None` when that
/// would split a character (the source edits bytes; such text cannot parse).
fn substitute(text: &str, index: usize, residue: char) -> Option<String> {
    if index >= text.len() || !text.is_char_boundary(index) || !text.is_char_boundary(index + 1) {
        return None;
    }
    let mut out = String::with_capacity(text.len() + residue.len_utf8());
    out.push_str(&text[..index]);
    out.push(residue);
    out.push_str(&text[index + 1..]);
    Some(out)
}

/// The byte at `index` as the source's `char`.
fn byte_char(text: &str, index: usize) -> char {
    text.as_bytes().get(index).map_or('\0', |&b| char::from(b))
}

/// Source `ProteaseDigestion::digest(protein, ranges, min_length,
/// max_length)`: the digestor's enzyme, specificity and missed cleavages with
/// the caller's length window (0 = no maximum). Unspecific cleavage products
/// come in the source's length-major order.
fn digest_ranges(
    digestor: &ProteaseDigestion,
    sequence: &str,
    min_length: usize,
    max_length: usize,
) -> Result<Vec<(usize, usize)>> {
    if digestor.specificity == DigestionSpecificity::None {
        return Err(Error::InvalidValue(
            "Specificity value set on current ProteaseDigestion object is not supported by \
             ProteaseDigestion::digest()."
                .into(),
        ));
    }
    let mut local = digestor.clone();
    local.min_length = min_length.max(1);
    local.max_length = (max_length != 0).then_some(max_length);
    if local.max_length.is_some_and(|max| max < local.min_length) {
        return Ok(Vec::new());
    }
    let mut products = local.digest_ranges(sequence)?;
    if local.enzyme == Protease::UnspecificCleavage {
        // The source tokenizes at every residue and emits by missed-cleavage
        // level, which for unspecific cleavage is by length, then start.
        products.sort_by_key(|p| (p.len(), p.start));
    }
    Ok(products.into_iter().map(|p| (p.start, p.end)).collect())
}

/// The simple variants inside `[start, end)`, skipping stop codons, unset
/// residues and unknown positions.
fn local_variants(
    variants: &[PEFFVariantSimple],
    start: usize,
    end: usize,
) -> Vec<&PEFFVariantSimple> {
    variants
        .iter()
        .filter(|v| v.variant_aa != '*' && v.variant_aa != '\0' && v.position != 0)
        .filter(|v| (start..end).contains(&(v.position - 1)))
        .collect()
}

/// The modifications of known position inside `[start, end)`.
fn local_modifications(
    modifications: &[PEFFModification],
    start: usize,
    end: usize,
) -> Vec<&PEFFModification> {
    modifications
        .iter()
        .filter(|m| m.position != 0 && (start..end).contains(&(m.position - 1)))
        .collect()
}

/// Source truncation at `MAX_COMBINATORIAL_ELEMENTS`, variants first.
fn cap_elements<V, M>(variants: &mut Vec<V>, modifications: &mut Vec<M>) -> bool {
    if variants.len() + modifications.len() <= MAX_COMBINATORIAL_ELEMENTS {
        return false;
    }
    if variants.len() > MAX_COMBINATORIAL_ELEMENTS {
        variants.truncate(MAX_COMBINATORIAL_ELEMENTS);
        modifications.clear();
    } else {
        modifications.truncate(MAX_COMBINATORIAL_ELEMENTS - variants.len());
    }
    true
}

/// Apply the variants selected by `combo` to `reference`, the peptide at
/// protein offset `start`. `None` for two variants at one position (or a
/// substitution that cannot be made). The description joins `E4L(tag)` items
/// with `;`.
fn apply_variants(
    reference: &str,
    start: usize,
    variants: &[&PEFFVariantSimple],
    combo: u64,
    description: &mut String,
) -> Option<String> {
    let mut peptide = reference.to_owned();
    let mut positions = BTreeSet::new();
    for (bit, variant) in variants.iter().enumerate() {
        if combo & (1u64 << bit) == 0 {
            continue;
        }
        let local = variant.position - 1 - start;
        if !positions.insert(local) {
            return None;
        }
        if !description.is_empty() {
            description.push(';');
        }
        description.push(byte_char(reference, local));
        description.push_str(&variant.position.to_string());
        description.push(variant.variant_aa);
        if !variant.optional_tag.is_empty() {
            description.push('(');
            description.push_str(&variant.optional_tag);
            description.push(')');
        }
        peptide = substitute(&peptide, local, variant.variant_aa)?;
    }
    Some(peptide)
}

/// Borrow a list of modification names.
fn names(list: &[String]) -> Vec<&str> {
    list.iter().map(String::as_str).collect()
}

/// `ACCESSION@position`, or `name@position` when there is no accession.
fn describe_modification(
    description: &mut String,
    modification: &PEFFModification,
    position: usize,
) {
    if !description.is_empty() {
        description.push(';');
    }
    if !modification.accession.is_empty() {
        description.push_str(&modification.accession);
    } else if !modification.name.is_empty() {
        description.push_str(&modification.name);
    }
    if !modification.accession.is_empty() || !modification.name.is_empty() {
        description.push('@');
        description.push_str(&position.to_string());
    }
}

impl PEFFEntry {
    /// The unmodified sequence as an [`AASequence`].
    ///
    /// Source `getSequence`, `AASequence::fromString(sequence)`; a sequence
    /// that does not parse is an error, as the source's exception.
    pub fn aa_sequence(&self) -> Result<AASequence> {
        AASequence::parse(&self.sequence)
    }

    /// The sequence with every localised annotated modification applied,
    /// resolved against [`ModificationsDB::global`].
    ///
    /// Source `getModifiedSequence`. See
    /// [`modified_sequence_with_registry`](Self::modified_sequence_with_registry).
    pub fn modified_sequence(&self) -> Result<AASequence> {
        Ok(self
            .modified_sequence_with_registry(ModificationsDB::global())?
            .0)
    }

    /// The sequence with every localised annotated modification applied,
    /// and the source's warnings.
    ///
    /// Modifications of unknown position (0) or beyond the sequence are
    /// skipped silently. Each other one resolves by accession, then by name,
    /// against `db` for its residue with anywhere specificity; one that does
    /// not resolve is skipped with the warning `Could not resolve
    /// modification: <accession> (<name>)`.
    pub fn modified_sequence_with_registry(
        &self,
        db: &ModificationsDB,
    ) -> Result<(AASequence, Vec<String>)> {
        let mut sequence = AASequence::parse_with_registry(&self.sequence, db)?;
        let mut warnings = Vec::new();
        for modification in &self.modifications {
            if modification.position == 0 || modification.position > sequence.len() {
                continue;
            }
            if !apply(db, &mut sequence, modification.position - 1, modification) {
                warnings.push(format!(
                    "Could not resolve modification: {} ({})",
                    modification.accession, modification.name
                ));
            }
        }
        Ok((sequence, warnings))
    }

    /// Every variant applied alone: the simple variants in order, then, with
    /// `include_complex`, the complex ones. Returns index-aligned
    /// descriptions and sequences.
    ///
    /// Source `getVariantSequences`. A simple variant is described
    /// `P3A (tag)`, a complex one `2-4>XXX (tag)` (the tag part only when
    /// set). Variants at position 0, beyond the sequence or without a residue,
    /// and complex variants with start 0 or start after end, are skipped. A
    /// variant sequence that does not parse is skipped whole; the source
    /// keeps its description, which misaligns its two vectors.
    pub fn variant_sequences(&self, include_complex: bool) -> PEFFPeptides {
        let mut out = PEFFPeptides::default();
        let sequence = &self.sequence;
        for variant in &self.simple_variants {
            if variant.position == 0
                || variant.position > sequence.len()
                || variant.variant_aa == '\0'
            {
                continue;
            }
            let index = variant.position - 1;
            let mut description = String::new();
            description.push(byte_char(sequence, index));
            description.push_str(&variant.position.to_string());
            description.push(variant.variant_aa);
            if !variant.optional_tag.is_empty() {
                description.push_str(&format!(" ({})", variant.optional_tag));
            }
            if let Some(parsed) = substitute(sequence, index, variant.variant_aa)
                .and_then(|s| AASequence::parse(&s).ok())
            {
                out.descriptions.push(description);
                out.sequences.push(parsed);
            }
        }
        if include_complex {
            for variant in &self.complex_variants {
                if variant.start_position == 0
                    || variant.end_position > sequence.len()
                    || variant.start_position > variant.end_position
                {
                    continue;
                }
                let (Some(head), Some(tail)) = (
                    sequence.get(..variant.start_position - 1),
                    sequence.get(variant.end_position..),
                ) else {
                    continue;
                };
                let variant_sequence = format!("{head}{}{tail}", variant.replacement);
                let mut description = format!(
                    "{}-{}>{}",
                    variant.start_position, variant.end_position, variant.replacement
                );
                if !variant.optional_tag.is_empty() {
                    description.push_str(&format!(" ({})", variant.optional_tag));
                }
                if let Ok(parsed) = AASequence::parse(&variant_sequence) {
                    out.descriptions.push(description);
                    out.sequences.push(parsed);
                }
            }
        }
        out
    }

    /// The sequence of the first processed region with `region_accession`
    /// that yields one; an empty sequence when none does.
    ///
    /// Source `getProcessedSequence` (default accession [`SIGNAL_PEPTIDE`]).
    /// For a signal peptide — the accession `PEFF:0001021` or a region named
    /// `signal peptide` — this is the mature protein after the region, when
    /// the region ends before the sequence does. For any other region it is
    /// the region itself, when its bounds are valid. Annotated modifications
    /// are not applied.
    pub fn processed_sequence(&self, region_accession: &str) -> Result<AASequence> {
        let sequence = &self.sequence;
        for region in &self.processed_regions {
            if region.accession != region_accession {
                continue;
            }
            if region_accession == SIGNAL_PEPTIDE || region.name == "signal peptide" {
                if region.end_position < sequence.len() {
                    if let Some(mature) = sequence.get(region.end_position..) {
                        return AASequence::parse(mature);
                    }
                }
            } else if region.start_position > 0
                && region.end_position >= region.start_position
                && region.end_position <= sequence.len()
            {
                if let Some(slice) = sequence.get(region.start_position - 1..region.end_position) {
                    return AASequence::parse(slice);
                }
            }
        }
        Ok(AASequence::default())
    }

    /// Digest the reference sequence and enumerate, per peptide, every
    /// combination of its simple variants and annotated modifications,
    /// resolved against [`ModificationsDB::global`].
    ///
    /// Source `digestWithVariants`; see
    /// [`digest_with_variants_with_registry`](Self::digest_with_variants_with_registry).
    pub fn digest_with_variants(
        &self,
        digestor: &ProteaseDigestion,
        options: &VariantDigestOptions,
    ) -> Result<PEFFPeptides> {
        self.digest_with_variants_with_registry(digestor, options, ModificationsDB::global())
    }

    /// Source `digestWithVariants` against a caller registry.
    ///
    /// The digest uses the digestor's enzyme, specificity and missed
    /// cleavages with the options' length window. For each peptide, in
    /// digest order: the reference (empty description) when requested, then
    /// the combinations 1 .. 2^(v+m) of its v variants (low bits) and m
    /// modifications. Stop-codon variants are never used. Combinations
    /// placing two variants on one residue, or with a modification that does
    /// not resolve on its (possibly substituted) residue, are dropped. Beyond
    /// [`MAX_COMBINATORIAL_ELEMENTS`] the list is cut, variants first, with
    /// the source's warning. Descriptions join `E4L(tag)` and
    /// `UNIMOD:35@4` (protein positions) with `;`.
    pub fn digest_with_variants_with_registry(
        &self,
        digestor: &ProteaseDigestion,
        options: &VariantDigestOptions,
        db: &ModificationsDB,
    ) -> Result<PEFFPeptides> {
        let mut out = PEFFPeptides::default();
        if self.sequence.is_empty() {
            return Ok(out);
        }
        AASequence::parse_with_registry(&self.sequence, db)?;
        let ranges = digest_ranges(
            digestor,
            &self.sequence,
            options.min_length,
            options.max_length,
        )?;
        for (start, end) in ranges {
            let mut variants = if options.include_variants {
                local_variants(&self.simple_variants, start, end)
            } else {
                Vec::new()
            };
            let mut modifications = if options.include_modifications {
                local_modifications(&self.modifications, start, end)
            } else {
                Vec::new()
            };
            let reference = &self.sequence[start..end];
            let total = variants.len() + modifications.len();
            if cap_elements(&mut variants, &mut modifications) {
                out.warnings.push(format!(
                    "Peptide at position {start}-{end} has {total} variants/modifications. \
                     Limiting to first {MAX_COMBINATORIAL_ELEMENTS} to avoid combinatorial explosion."
                ));
            }
            if options.include_reference {
                if let Ok(parsed) = AASequence::parse_with_registry(reference, db) {
                    out.push(String::new(), parsed, options.max_peptides)?;
                }
            }
            let elements = variants.len() + modifications.len();
            if elements == 0 {
                continue;
            }
            'combos: for combo in 1..(1u64 << elements) {
                let mut description = String::new();
                let Some(variant_peptide) =
                    apply_variants(reference, start, &variants, combo, &mut description)
                else {
                    continue;
                };
                let Ok(mut peptide) = AASequence::parse_with_registry(&variant_peptide, db) else {
                    continue;
                };
                for (i, modification) in modifications.iter().enumerate() {
                    if combo & (1u64 << (variants.len() + i)) == 0 {
                        continue;
                    }
                    describe_modification(&mut description, modification, modification.position);
                    if !apply(
                        db,
                        &mut peptide,
                        modification.position - 1 - start,
                        modification,
                    ) {
                        continue 'combos;
                    }
                }
                out.push(description, peptide, options.max_peptides)?;
            }
        }
        Ok(out)
    }

    /// Digest, enumerate variants and annotated modifications, then apply
    /// fixed and variable sample-handling modifications, resolved against
    /// [`ModificationsDB::global`].
    ///
    /// Source `generatePeptides`; see
    /// [`generate_peptides_with_registry`](Self::generate_peptides_with_registry).
    pub fn generate_peptides(
        &self,
        digestor: &ProteaseDigestion,
        options: &PeptideGenerationOptions,
    ) -> Result<PEFFPeptides> {
        self.generate_peptides_with_registry(digestor, options, ModificationsDB::global())
    }

    /// Source `generatePeptides` against a caller registry.
    ///
    /// For each digested peptide and each combination of its simple variants
    /// (bit order, combination 0 = none; two variants on one residue drop the
    /// combination) the annotated modifications are enumerated as by
    /// [`enumerate_peff_modifications`](Self::enumerate_peff_modifications)
    /// with the proteoform cap. Without `include_reference`, the
    /// variant-free peptide is skipped when nothing is annotated in it, and
    /// otherwise its unmodified form is removed from the enumeration. Fixed
    /// modifications then apply to every peptide; with variable
    /// modifications each peptide expands (keeping the original), and a
    /// peptide the variable modifications changed gets `+varmod` appended to
    /// its description. Unknown fixed or variable modification names are
    /// errors, as the source's exception.
    pub fn generate_peptides_with_registry(
        &self,
        digestor: &ProteaseDigestion,
        options: &PeptideGenerationOptions,
        db: &ModificationsDB,
    ) -> Result<PEFFPeptides> {
        let mut out = PEFFPeptides::default();
        if self.sequence.is_empty() {
            return Ok(out);
        }
        AASequence::parse_with_registry(&self.sequence, db)?;
        let ranges = digest_ranges(
            digestor,
            &self.sequence,
            options.min_length,
            options.max_length,
        )?;
        let fixed = if options.fixed_mods.is_empty() {
            Vec::new()
        } else {
            ModifiedPeptideGenerator::get_modifications_with_registry(
                &names(&options.fixed_mods),
                db,
            )?
        };
        let variable = if options.variable_mods.is_empty() {
            Vec::new()
        } else {
            ModifiedPeptideGenerator::get_modifications_with_registry(
                &names(&options.variable_mods),
                db,
            )?
        };
        let generator = ModifiedPeptideGenerator::default();
        for (start, end) in ranges {
            let mut variants = if options.include_peff_variants {
                local_variants(&self.simple_variants, start, end)
            } else {
                Vec::new()
            };
            let mut modifications: Vec<(usize, &PEFFModification)> =
                if options.include_peff_modifications {
                    local_modifications(&self.modifications, start, end)
                        .into_iter()
                        .map(|m| (m.position - 1 - start, m))
                        .collect()
                } else {
                    Vec::new()
                };
            let total = variants.len() + modifications.len();
            if cap_elements(&mut variants, &mut modifications) {
                out.warnings.push(format!(
                    "Peptide at {start}-{end} has {total} variants+modifications. Limiting to \
                     {MAX_COMBINATORIAL_ELEMENTS} to avoid combinatorial explosion."
                ));
            }
            let reference = &self.sequence[start..end];
            for variant_combo in 0..(1u64 << variants.len()) {
                if variant_combo == 0
                    && !options.include_reference
                    && variants.is_empty()
                    && modifications.is_empty()
                {
                    continue;
                }
                let mut variant_description = String::new();
                let Some(variant_text) = apply_variants(
                    reference,
                    start,
                    &variants,
                    variant_combo,
                    &mut variant_description,
                ) else {
                    continue;
                };
                if variant_combo == 0 && !options.include_reference && modifications.is_empty() {
                    continue;
                }
                let Ok(variant_peptide) = AASequence::parse_with_registry(&variant_text, db) else {
                    continue;
                };
                let (mut enumerated, warning) = enumerate_with_registry(
                    &variant_peptide,
                    &modifications,
                    &variant_description,
                    options.max_peff_mods_per_proteoform,
                    db,
                );
                out.warnings.extend(warning);
                if variant_combo == 0 && !options.include_reference && !modifications.is_empty() {
                    let unmodified = variant_peptide.to_string();
                    enumerated.retain(|(_, peptide)| peptide.to_string() != unmodified);
                }
                for (description, mut peptide) in enumerated {
                    if !fixed.is_empty() {
                        generator.apply_fixed_modifications(&fixed, &mut peptide)?;
                    }
                    if variable.is_empty() {
                        out.push(description, peptide, options.max_peptides)?;
                        continue;
                    }
                    let expanded = generator.variable_modifications(
                        &variable,
                        &peptide,
                        options.max_variable_mods_per_peptide,
                        true,
                    )?;
                    let original = peptide.to_string();
                    for candidate in expanded {
                        let mut full = description.clone();
                        if candidate.is_modified() && candidate.to_string() != original {
                            if !full.is_empty() {
                                full.push(';');
                            }
                            full.push_str("+varmod");
                        }
                        out.push(full, candidate, options.max_peptides)?;
                    }
                }
            }
        }
        Ok(out)
    }

    /// Every subset of `modifications` applied to `peptide`, in bit order
    /// (subset 0, the unmodified peptide, first), against
    /// [`ModificationsDB::global`].
    ///
    /// Source private static `enumeratePEFFModifications_`, public here so
    /// that callers can enumerate proteoforms of their own peptides. Each
    /// pair holds a zero-based position in `peptide` and the modification.
    /// Descriptions extend `base_description` with `;ACCESSION@p` (or
    /// `name@p`), p one-based in the peptide. `max_selected` > 0 skips subsets
    /// larger than it. Subsets with a modification that does not resolve are
    /// dropped. At most [`MAX_COMBINATORIAL_ELEMENTS`] modifications are
    /// used; the second value is the source's warning when more were given.
    pub fn enumerate_peff_modifications(
        peptide: &AASequence,
        modifications: &[(usize, &PEFFModification)],
        base_description: &str,
        max_selected: usize,
    ) -> (Vec<(String, AASequence)>, Option<String>) {
        enumerate_with_registry(
            peptide,
            modifications,
            base_description,
            max_selected,
            ModificationsDB::global(),
        )
    }
}

fn enumerate_with_registry(
    peptide: &AASequence,
    modifications: &[(usize, &PEFFModification)],
    base_description: &str,
    max_selected: usize,
    db: &ModificationsDB,
) -> (Vec<(String, AASequence)>, Option<String>) {
    let mut result = Vec::new();
    if modifications.is_empty() {
        result.push((base_description.to_owned(), peptide.clone()));
        return (result, None);
    }
    let warning = (modifications.len() > MAX_COMBINATORIAL_ELEMENTS).then(|| {
        format!(
            "Peptide has {} PEFF modifications. Limiting to first {MAX_COMBINATORIAL_ELEMENTS} \
             to avoid combinatorial explosion.",
            modifications.len()
        )
    });
    let used = &modifications[..modifications.len().min(MAX_COMBINATORIAL_ELEMENTS)];
    'combos: for combo in 0..(1u64 << used.len()) {
        if max_selected > 0 && combo.count_ones() as usize > max_selected {
            continue;
        }
        let mut modified = peptide.clone();
        let mut description = base_description.to_owned();
        for (bit, (position, modification)) in used.iter().enumerate() {
            if combo & (1u64 << bit) == 0 {
                continue;
            }
            describe_modification(&mut description, modification, position + 1);
            if !apply(db, &mut modified, *position, modification) {
                continue 'combos;
            }
        }
        result.push((description, modified));
    }
    (result, warning)
}
