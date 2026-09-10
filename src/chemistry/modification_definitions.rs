// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Modification search definitions, source-compatible compatibility predicates,
//! mass matching, and inference from all peptide hits. See
//! `docs/MODIFICATION_DEFINITIONS_SUPPORT.md` for the terminal-site quirks and
//! the distinction between stored occurrence counts and resource limits.

use super::{
    AASequence, MassTag, ModificationsDB, ResidueModification, SequenceModification,
    TermSpecificity, composition_formula, residue_composition,
};
use crate::identification::PeptideIdentification;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::hash::{Hash, Hasher};

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

/// One owned chemical definition and its search settings. An unset default has
/// an empty name; accessing its chemistry or inserting it into a set is an error.
/// Equality includes chemical contents, fixed status, and occurrence count.
/// There is deliberately no `Ord`: the source's full-ID-only set identity does
/// not agree with its equality relation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModificationDefinition {
    modification: Option<DefinedModification>,
    pub fixed: bool,
    /// Stored search setting; zero denotes unlimited. Compatibility and inference
    /// do not enforce this field, matching the source.
    pub max_occurrences: u32,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(clippy::large_enum_variant)] // Keep named records inline; shared tag text makes MassTag smaller.
enum DefinedModification {
    Named(ResidueModification),
    MassTag(MassTag),
}

impl DefinedModification {
    fn full_id(&self) -> &str {
        match self {
            Self::Named(modification) => modification.full_id(),
            Self::MassTag(tag) => tag.full_id(),
        }
    }

    fn origin(&self) -> Option<char> {
        match self {
            Self::Named(modification) => modification.origin(),
            Self::MassTag(tag) => tag.origin(),
        }
    }

    fn term_specificity(&self) -> TermSpecificity {
        match self {
            Self::Named(modification) => modification.term_specificity(),
            Self::MassTag(tag) => tag.term_specificity(),
        }
    }

    // Source user-defined ResidueModification has an empty short ID. Its full
    // ID, checked separately by compatibility, carries the numeric spelling.
    fn short_name(&self) -> &str {
        match self {
            Self::Named(modification) => modification.name(),
            Self::MassTag(_) => "",
        }
    }

    fn delta_mass(&self) -> Result<f64> {
        match self {
            Self::Named(modification) => Ok(modification.diff_mono_mass()),
            Self::MassTag(tag) => {
                if tag.is_delta() {
                    Ok(tag.mass())
                } else if tag.term_specificity() != TermSpecificity::Anywhere {
                    Ok(tag.mass() - terminal_mass(tag.term_specificity()))
                } else {
                    // Retain the source full-residue-minus-water arithmetic.
                    // B/Z/X have no native unmodified mass to subtract.
                    Ok(tag.mass()
                        - internal_residue_mass(
                            &tag.origin()
                                .ok_or_else(|| invalid("residue tag has no origin"))?
                                .to_string(),
                        )?)
                }
            }
        }
    }

    fn absolute_mass(&self) -> Result<f64> {
        let mass = match self {
            Self::Named(modification) => modification.mono_mass(),
            Self::MassTag(tag) if tag.term_specificity() != TermSpecificity::Anywhere => {
                if tag.is_delta() {
                    tag.mass() + terminal_mass(tag.term_specificity())
                } else {
                    tag.mass()
                }
            }
            Self::MassTag(tag) => {
                let origin = tag
                    .origin()
                    .ok_or_else(|| invalid("residue tag has no origin"))?;
                if matches!(origin, 'B' | 'Z' | 'X') {
                    // An absolute internal tag carries a known modified mass
                    // without establishing any unmodified residue mass.
                    tag.residue_mono_mass()
                        .ok_or_else(|| invalid("residue tag has no internal mass"))?
                        + water_mass()
                } else {
                    // createUnknownFromMassString stores the FULL modified
                    // residue mass, not the tag's internal-residue mass.
                    self.delta_mass()? + full_residue_mass(&origin.to_string())?
                }
            }
        };
        if mass.is_finite() {
            Ok(mass)
        } else {
            Err(invalid("absolute modification mass overflows"))
        }
    }
}

fn water_mass() -> f64 {
    composition_formula([0, 2, 0, 1, 0, 0]).mono_mass()
}

fn terminal_mass(term: TermSpecificity) -> f64 {
    match term {
        TermSpecificity::NTerm | TermSpecificity::ProteinNTerm => {
            composition_formula([0, 1, 0, 0, 0, 0]).mono_mass()
        }
        _ => composition_formula([0, 1, 0, 1, 0, 0]).mono_mass(),
    }
}

impl Default for ModificationDefinition {
    fn default() -> Self {
        Self {
            modification: None,
            fixed: true,
            max_occurrences: 0,
        }
    }
}

impl Hash for ModificationDefinition {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.modification_name().hash(state);
        self.fixed.hash(state);
        self.max_occurrences.hash(state);
    }
}

impl ModificationDefinition {
    pub fn new(name: &str) -> Result<Self> {
        Self::with_options(name, true, 0)
    }

    pub fn with_options(name: &str, fixed: bool, max_occurrences: u32) -> Result<Self> {
        Ok(Self::from_modification(
            ModificationsDB::global().get_modification(name, None, None)?,
            fixed,
            max_occurrences,
        ))
    }

    /// Own a copy, including when the record belongs to a custom native registry.
    pub fn from_modification(
        modification: &ResidueModification,
        fixed: bool,
        max_occurrences: u32,
    ) -> Self {
        Self {
            modification: Some(DefinedModification::Named(modification.clone())),
            fixed,
            max_occurrences,
        }
    }

    pub fn modification(&self) -> Result<&ResidueModification> {
        match self.value()? {
            DefinedModification::Named(modification) => Ok(modification),
            DefinedModification::MassTag(_) => Err(Error::Unsupported(
                "anonymous definition has a mass tag, not a named residue modification".into(),
            )),
        }
    }

    /// Clone a named or anonymous sequence annotation without modifying a registry.
    pub fn from_annotation(
        annotation: &SequenceModification,
        fixed: bool,
        max_occurrences: u32,
    ) -> Self {
        match annotation {
            SequenceModification::Known(modification) => {
                Self::from_modification(modification, fixed, max_occurrences)
            }
            SequenceModification::MassTag(tag) => Self {
                modification: Some(DefinedModification::MassTag(tag.clone())),
                fixed,
                max_occurrences,
            },
        }
    }

    pub fn mass_tag(&self) -> Option<&MassTag> {
        match &self.modification {
            Some(DefinedModification::MassTag(tag)) => Some(tag),
            _ => None,
        }
    }

    /// Replace chemistry while keeping the fixed flag and occurrence setting.
    pub fn set_annotation(&mut self, annotation: &SequenceModification) {
        *self = Self::from_annotation(annotation, self.fixed, self.max_occurrences);
    }

    fn value(&self) -> Result<&DefinedModification> {
        self.modification
            .as_ref()
            .ok_or_else(|| invalid("modification definition is unset"))
    }

    pub fn modification_name(&self) -> &str {
        self.modification
            .as_ref()
            .map_or("", DefinedModification::full_id)
    }

    /// Resolve before replacing, preserving the old value on an invalid name.
    pub fn set_modification(&mut self, name: &str) -> Result<()> {
        self.modification = Some(DefinedModification::Named(
            ModificationsDB::global()
                .get_modification(name, None, None)?
                .clone(),
        ));
        Ok(())
    }
}

/// Which monoisotopic mass is compared by `find_matches`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ModificationMassMode {
    #[default]
    Delta,
    /// Declared absolute mass, or delta plus internal residue mass when the
    /// declared value is nonpositive and a residue was supplied. With no residue,
    /// the declared value is compared literally (zero in bundled records).
    Absolute,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ModificationMatchOptions {
    /// Empty, ".", and "X" disable origin filtering. Other strings compare
    /// their first character before absolute mode resolves a residue alias.
    pub residue: String,
    pub term_specificity: Option<TermSpecificity>,
    pub consider_fixed: bool,
    pub consider_variable: bool,
    pub mass_mode: ModificationMassMode,
    /// Inclusive absolute tolerance in Da, finite and nonnegative.
    pub tolerance: f64,
}

impl Default for ModificationMatchOptions {
    fn default() -> Self {
        Self {
            residue: String::new(),
            term_specificity: None,
            consider_fixed: true,
            consider_variable: true,
            mass_mode: ModificationMassMode::Delta,
            tolerance: 0.01,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ModificationMatch {
    pub mass_error: f64,
    pub definition: ModificationDefinition,
}

/// Separate fixed and variable sets keyed by full modification ID. Inserting a
/// duplicate into one partition keeps the first definition; both partitions may
/// contain the same ID. Returned references cannot mutate the map's identity.
#[derive(Clone, Debug)]
pub struct ModificationDefinitionsSet {
    fixed: BTreeMap<String, ModificationDefinition>,
    variable: BTreeMap<String, ModificationDefinition>,
    /// Stored source setting, ignored by compatibility and inference. Zero is
    /// unlimited. This is unrelated to the native resource limit below.
    pub max_modifications: usize,
    /// Positive logical visit budget per checked operation (default 1,000,000).
    /// Counts definitions, input names/bytes, or peptide sites as documented.
    /// Immutable registry initialization and caller-owned metadata are excluded.
    pub max_work: usize,
}

impl Default for ModificationDefinitionsSet {
    fn default() -> Self {
        Self {
            fixed: BTreeMap::new(),
            variable: BTreeMap::new(),
            max_modifications: 0,
            max_work: 1_000_000,
        }
    }
}

// Resource policy does not change equality of the source's search definition.
impl PartialEq for ModificationDefinitionsSet {
    fn eq(&self, other: &Self) -> bool {
        self.fixed == other.fixed
            && self.variable == other.variable
            && self.max_modifications == other.max_modifications
    }
}
impl Eq for ModificationDefinitionsSet {}

impl ModificationDefinitionsSet {
    pub fn from_names(fixed: &[&str], variable: &[&str]) -> Result<Self> {
        let mut definitions = Self::default();
        definitions.set_names(fixed, variable)?;
        Ok(definitions)
    }

    /// Comma-separated full names. As in ListUtils, tokens are not trimmed.
    pub fn from_comma_separated(fixed: &str, variable: &str) -> Result<Self> {
        let mut definitions = Self::default();
        definitions.set_comma_separated(fixed, variable)?;
        Ok(definitions)
    }

    /// Count both partitions, including IDs present in both.
    pub fn len(&self) -> usize {
        self.fixed.len() + self.variable.len()
    }

    pub fn is_empty(&self) -> bool {
        self.fixed.is_empty() && self.variable.is_empty()
    }

    pub fn fixed_modifications(
        &self,
    ) -> impl ExactSizeIterator<Item = &ModificationDefinition> + DoubleEndedIterator {
        self.fixed.values()
    }

    pub fn variable_modifications(
        &self,
    ) -> impl ExactSizeIterator<Item = &ModificationDefinition> + DoubleEndedIterator {
        self.variable.values()
    }

    pub fn fixed_names(&self) -> BTreeSet<String> {
        self.fixed.keys().cloned().collect()
    }

    pub fn variable_names(&self) -> BTreeSet<String> {
        self.variable.keys().cloned().collect()
    }

    pub fn modification_names(&self) -> BTreeSet<String> {
        self.fixed
            .keys()
            .chain(self.variable.keys())
            .cloned()
            .collect()
    }

    /// Merged full-ID order; a fixed definition wins an ID shared by both sets.
    pub fn modifications(&self) -> Vec<&ModificationDefinition> {
        let mut merged: BTreeMap<&str, &ModificationDefinition> = self
            .variable
            .iter()
            .map(|(name, definition)| (name.as_str(), definition))
            .collect();
        merged.extend(
            self.fixed
                .iter()
                .map(|(name, definition)| (name.as_str(), definition)),
        );
        merged.into_values().collect()
    }

    pub fn add_modification(&mut self, definition: ModificationDefinition) -> Result<()> {
        definition.value()?;
        let mut work = self.work()?;
        charge(&mut work, self.len())?;
        charge(&mut work, definition.modification_name().len())?;
        self.insert(definition);
        Ok(())
    }

    /// Replace from a source-style set: the first full ID in the input wins
    /// globally, including when later entries have a different fixed flag.
    /// `add_modification` and `set_names` instead deduplicate within partitions.
    pub fn set_modifications(&mut self, definitions: &[ModificationDefinition]) -> Result<()> {
        let mut work = self.work()?;
        charge(&mut work, definitions.len())?;
        for definition in definitions {
            definition.value()?;
            charge(&mut work, definition.modification_name().len())?;
        }
        let mut replacement = self.empty_like();
        let mut seen = BTreeSet::new();
        for definition in definitions {
            if seen.insert(definition.modification_name()) {
                replacement.insert(definition.clone());
            }
        }
        *self = replacement;
        Ok(())
    }

    pub fn set_names(&mut self, fixed: &[&str], variable: &[&str]) -> Result<()> {
        let mut work = self.work()?;
        charge(&mut work, fixed.len())?;
        charge(&mut work, variable.len())?;
        for name in fixed.iter().chain(variable) {
            charge(&mut work, name.len())?;
        }
        let mut replacement = self.empty_like();
        for (names, is_fixed) in [(fixed, true), (variable, false)] {
            for name in names {
                replacement.insert(ModificationDefinition::with_options(name, is_fixed, 0)?);
            }
        }
        *self = replacement;
        Ok(())
    }

    pub fn set_comma_separated(&mut self, fixed: &str, variable: &str) -> Result<()> {
        let mut work = self.work()?;
        charge(&mut work, fixed.len())?;
        charge(&mut work, variable.len())?;
        fn split(text: &str) -> Vec<&str> {
            if text.is_empty() {
                Vec::new()
            } else {
                text.split(',').collect::<Vec<_>>()
            }
        }
        self.set_names(&split(fixed), &split(variable))
    }

    /// Preserve the source predicate, including its residue-based fixed-terminal
    /// check. This is not a validator that all fixed terminal slots are occupied.
    pub fn is_compatible(&self, peptide: &AASequence) -> Result<bool> {
        let mut work = self.work()?;
        charge(&mut work, self.len())?;
        charge(&mut work, peptide.len())?;
        charge(&mut work, 2)?;
        charge(
            &mut work,
            self.fixed
                .len()
                .checked_mul(peptide.len())
                .ok_or_else(|| invalid("modification compatibility work overflows"))?,
        )?;
        if self.fixed.is_empty() && !peptide.is_modified() {
            return Ok(true);
        }
        for definition in self.fixed.values() {
            let modification = definition.value()?;
            let origin = modification.origin().unwrap_or('X');
            for (index, residue) in peptide.as_str().chars().enumerate() {
                if residue == origin
                    && peptide.residue_modification(index)?.is_none_or(|actual| {
                        actual.known().map_or("", ResidueModification::name)
                            != modification.short_name()
                    })
                {
                    return Ok(false);
                }
            }
        }
        let allowed = |modification: Option<&SequenceModification>| {
            modification.is_none_or(|m| {
                self.fixed.contains_key(m.full_id()) || self.variable.contains_key(m.full_id())
            })
        };
        for index in 0..peptide.len() {
            if !allowed(peptide.residue_modification(index)?) {
                return Ok(false);
            }
        }
        Ok(
            allowed(peptide.n_terminal_modification())
                && allowed(peptide.c_terminal_modification()),
        )
    }

    /// Return increasing absolute mass errors. Exact-error ties retain fixed
    /// entries before variable entries, with full-ID order inside each partition.
    pub fn find_matches(
        &self,
        mass: f64,
        options: &ModificationMatchOptions,
    ) -> Result<Vec<ModificationMatch>> {
        if !mass.is_finite() || !options.tolerance.is_finite() || options.tolerance < 0.0 {
            return Err(invalid(
                "mass and tolerance must be finite; tolerance must be nonnegative",
            ));
        }
        if !options.consider_fixed && !options.consider_variable {
            return Err(invalid(
                "mass matching needs fixed or variable modifications",
            ));
        }
        let mut work = self.work()?;
        charge(&mut work, self.len())?;
        charge(&mut work, options.residue.len())?;
        if !options.residue.is_ascii() {
            return Err(invalid("residue selectors must be ASCII"));
        }
        let mut matches = Vec::new();
        let mut residue_mass = None;
        for (partition, consider) in [
            (&self.fixed, options.consider_fixed),
            (&self.variable, options.consider_variable),
        ] {
            if !consider {
                continue;
            }
            for definition in partition.values() {
                let modification = definition.value()?;
                let origin = modification.origin().unwrap_or('X');
                let residue = options.residue.as_str();
                if !(residue.is_empty()
                    || origin == 'X'
                    || residue.starts_with(origin)
                    || residue == "."
                    || residue == "X")
                {
                    continue;
                }
                if options
                    .term_specificity
                    .is_some_and(|term| term != modification.term_specificity())
                {
                    continue;
                }
                let mod_mass = match options.mass_mode {
                    ModificationMassMode::Delta => modification.delta_mass()?,
                    ModificationMassMode::Absolute => {
                        let absolute = modification.absolute_mass()?;
                        if absolute > 0.0 || residue.is_empty() {
                            absolute
                        } else {
                            let weight = match residue_mass {
                                Some(weight) => weight,
                                None => {
                                    let weight = internal_residue_mass(residue)?;
                                    residue_mass = Some(weight);
                                    weight
                                }
                            };
                            modification.delta_mass()? + weight
                        }
                    }
                };
                let mass_error = (mod_mass - mass).abs();
                if !mass_error.is_finite() {
                    return Err(invalid("modification mass error overflows"));
                }
                if mass_error <= options.tolerance {
                    matches.push(ModificationMatch {
                        mass_error,
                        definition: definition.clone(),
                    });
                }
            }
        }
        matches.sort_by(|left, right| left.mass_error.total_cmp(&right.mass_error));
        Ok(matches)
    }

    /// Infer from all hits, pooling each residue identity and each terminus.
    /// One observed nonempty modification with no unmodified observation makes
    /// a fixed definition; otherwise every observed modification is variable.
    /// Scores, ranks, evidence, and stored occurrence limits are not consulted.
    /// Anonymous annotations remain owned mass tags, preserving exact spelling
    /// and missing chemistry without inventing a named registry entry.
    /// Equal chemical records share native value identity; conflicting chemistry
    /// under one full ID is an atomic error instead of an order-dependent choice.
    pub fn infer_from_peptides(&mut self, peptides: &[PeptideIdentification]) -> Result<()> {
        let mut work = self.work()?;
        charge(&mut work, peptides.len())?;
        for peptide in peptides {
            charge(&mut work, peptide.hits.len())?;
            for hit in &peptide.hits {
                charge(&mut work, hit.sequence.len())?;
                charge(&mut work, 2)?;
                for index in 0..hit.sequence.len() {
                    if let Some(annotation) = hit.sequence.residue_modification(index)? {
                        charge(&mut work, annotation.full_id().len())?;
                    }
                }
                for annotation in [
                    hit.sequence.n_terminal_modification(),
                    hit.sequence.c_terminal_modification(),
                ]
                .into_iter()
                .flatten()
                {
                    charge(&mut work, annotation.full_id().len())?;
                }
            }
        }
        let mut sites = ObservedSites::new();
        let mut records = BTreeMap::new();
        for peptide in peptides {
            for hit in &peptide.hits {
                let sequence = &hit.sequence;
                record_site(
                    &mut sites,
                    &mut records,
                    "N-term",
                    sequence.n_terminal_modification(),
                )?;
                record_site(
                    &mut sites,
                    &mut records,
                    "C-term",
                    sequence.c_terminal_modification(),
                )?;
                for index in 0..sequence.len() {
                    record_site(
                        &mut sites,
                        &mut records,
                        &sequence.as_str()[index..index + 1],
                        sequence.residue_modification(index)?,
                    )?;
                }
            }
        }
        let mut replacement = self.empty_like();
        for observations in sites.values() {
            let fixed = observations.len() == 1 && !observations.contains_key(&None);
            for modification in observations.values().flatten() {
                replacement.insert(ModificationDefinition::from_annotation(
                    modification,
                    fixed,
                    0,
                ));
            }
        }
        *self = replacement;
        Ok(())
    }

    fn insert(&mut self, definition: ModificationDefinition) {
        let partition = if definition.fixed {
            &mut self.fixed
        } else {
            &mut self.variable
        };
        partition
            .entry(definition.modification_name().to_owned())
            .or_insert(definition);
    }

    fn empty_like(&self) -> Self {
        Self {
            max_modifications: self.max_modifications,
            max_work: self.max_work,
            ..Self::default()
        }
    }

    fn work(&self) -> Result<usize> {
        if self.max_work == 0 {
            Err(invalid(
                "modification-definition work limit must be positive",
            ))
        } else {
            Ok(self.max_work)
        }
    }
}

fn charge(remaining: &mut usize, amount: usize) -> Result<()> {
    *remaining = remaining
        .checked_sub(amount)
        .ok_or_else(|| invalid("modification-definition work limit exceeded"))?;
    Ok(())
}

type ObservedSites<'a> =
    BTreeMap<&'a str, BTreeMap<Option<&'a str>, Option<&'a SequenceModification>>>;

fn record_site<'a>(
    sites: &mut ObservedSites<'a>,
    records: &mut BTreeMap<&'a str, &'a SequenceModification>,
    site: &'a str,
    annotation: Option<&'a SequenceModification>,
) -> Result<()> {
    if let Some(modification) = annotation {
        if records
            .insert(modification.full_id(), modification)
            .is_some_and(|previous| previous != modification)
        {
            return Err(invalid(
                "conflicting inferred modification records share a full ID",
            ));
        }
    }
    sites
        .entry(site)
        .or_default()
        .insert(annotation.map(SequenceModification::full_id), annotation);
    Ok(())
}

fn internal_residue_mass(name: &str) -> Result<f64> {
    Ok(full_residue_mass(name)? - water_mass())
}

fn full_residue_mass(name: &str) -> Result<f64> {
    // The aliases and case sensitivity are the pinned ResidueDB's names. Origin
    // filtering intentionally happens before this lookup, using the input's
    // first character, rather than the resolved one-letter code.
    let residue = match name {
        "A" | "Alanine" | "Ala" | "ALA" | "L-Alanine" | "alanine" | "Alanin" | "alanin" => b'A',
        "C" | "Cysteine" | "Cys" | "CYS" | "Cystine" => b'C',
        "D" | "Aspartate" | "Asp" | "ASP" => b'D',
        "E" | "Glutamate" | "Glu" | "GLU" => b'E',
        "F" | "Phenylalanine" | "Phe" | "PHE" => b'F',
        "G" | "Glycine" | "Gly" | "GLY" => b'G',
        "H" | "Histidine" | "His" | "HIS" => b'H',
        "I" | "Isoleucine" | "Ile" | "ILE" => b'I',
        "K" | "Lysine" | "Lys" | "LYS" => b'K',
        "L" | "Leucine" | "Leu" | "LEU" => b'L',
        "M" | "Methionine" | "Met" | "MET" => b'M',
        "N" | "Asparagine" | "Asn" | "ASN" => b'N',
        "P" | "Proline" | "Pro" | "PRO" => b'P',
        "Q" | "Glutamine" | "Gln" | "GLN" => b'Q',
        "R" | "Arginine" | "Arg" | "ARG" => b'R',
        "S" | "Serine" | "Ser" | "SER" => b'S',
        "T" | "Threonine" | "Thr" | "THR" => b'T',
        "U" | "Selenocysteine" | "Sec" | "SEC" => b'U',
        "V" | "Valine" | "Val" | "VAL" => b'V',
        "W" | "Tryptophan" | "Trp" | "TRP" => b'W',
        "Y" | "Tyrosine" | "Tyr" | "TYR" => b'Y',
        "O" | "Pyrrolysine" | "Pyr" | "PYR" => b'O',
        "J" | "Isoleucine/Leucine" | "Xle" | "XLE" => b'J',
        "B"
        | "Asparagine/Aspartate"
        | "Asx"
        | "ASX"
        | "Z"
        | "Glutamine/Glutamate"
        | "Glx"
        | "GLX"
        | "X"
        | "Unspecified/Unknown"
        | "Xaa"
        | "XAA"
        | "Unk" => {
            return Err(Error::Unsupported(
                "absolute mass matching requires a residue with known mass".into(),
            ));
        }
        _ => {
            return Err(invalid(
                "unknown residue for absolute modification mass matching",
            ));
        }
    };
    let composition = residue_composition(residue)
        .ok_or_else(|| invalid("missing native residue composition"))?;
    let water = composition_formula([0, 2, 0, 1, 0, 0]);
    // Match full-residue formula summation followed by water subtraction, not
    // an internal-formula mass or a sum of separately evaluated atom groups.
    Ok(composition_formula(composition)
        .checked_add(&water)?
        .mono_mass())
}
