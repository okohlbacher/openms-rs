// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Immutable modification registry backed by the pinned UniMod and OpenMS data.
//! The data has its own notices under resources/modifications; this implementation
//! is BSD-3-Clause. No network or runtime XML/resource lookup is required.

use super::EmpiricalFormula;
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::io::BufRead;
use std::sync::{Arc, OnceLock};

#[path = "obo.rs"]
mod obo;
pub use obo::{OboLoadReport, OboReadOptions};

/// Positional specificity from OpenMS ResidueModification.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum TermSpecificity {
    Anywhere,
    NTerm,
    CTerm,
    ProteinNTerm,
    ProteinCTerm,
}

impl TermSpecificity {
    /// OpenMS name used when constructing modification full identifiers.
    pub fn name(self) -> &'static str {
        match self {
            Self::Anywhere => "Anywhere",
            Self::NTerm => "N-term",
            Self::CTerm => "C-term",
            Self::ProteinNTerm => "Protein N-term",
            Self::ProteinCTerm => "Protein C-term",
        }
    }
}

/// A neutral loss declared for a particular modification specificity.
#[derive(Clone, Debug, PartialEq)]
pub struct NeutralLoss {
    formula: EmpiricalFormula,
    mono_mass: f64,
    average_mass: f64,
}
// Values are constructed only from finite, validated table fields.
impl Eq for NeutralLoss {}

impl Ord for NeutralLoss {
    fn cmp(&self, other: &Self) -> Ordering {
        formula_order(&self.formula, &other.formula)
            .then_with(|| finite_mass_order(self.mono_mass, other.mono_mass))
            .then_with(|| finite_mass_order(self.average_mass, other.average_mass))
    }
}
impl PartialOrd for NeutralLoss {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn formula_order(left: &EmpiricalFormula, right: &EmpiricalFormula) -> Ordering {
    left.atoms
        .cmp(&right.atoms)
        .then_with(|| left.charge.cmp(&right.charge))
}
fn finite_mass_order(left: f64, right: f64) -> Ordering {
    // Validated finite values have ordinary numeric ordering. Unlike total_cmp
    // alone, this retains PartialEq's equality of positive and negative zero.
    if left == right {
        Ordering::Equal
    } else {
        left.total_cmp(&right)
    }
}

impl NeutralLoss {
    pub fn new(formula: EmpiricalFormula, mono_mass: f64, average_mass: f64) -> Result<Self> {
        if !mono_mass.is_finite() || !average_mass.is_finite() {
            return Err(Error::InvalidValue(
                "neutral loss masses must be finite".into(),
            ));
        }
        Ok(Self {
            formula,
            mono_mass,
            average_mass,
        })
    }
    pub fn formula(&self) -> &EmpiricalFormula {
        &self.formula
    }
    pub fn mono_mass(&self) -> f64 {
        self.mono_mass
    }
    pub fn average_mass(&self) -> f64 {
        self.average_mass
    }
}

/// Mutable input description; construction validates names, origins and finite
/// masses, then freezes the result. Absolute formulas describe the free residue.
#[derive(Clone, Debug)]
pub struct ModificationRecord {
    pub record_id: Option<u32>,
    pub name: String,
    pub full_name: String,
    /// Empty derives the identifier from name and specificity.
    pub full_id: String,
    pub origin: Option<char>,
    pub term_specificity: TermSpecificity,
    pub diff_formula: EmpiricalFormula,
    pub absolute_formula: Option<EmpiricalFormula>,
    pub diff_mono_mass: f64,
    pub diff_average_mass: f64,
    pub mono_mass: f64,
    pub average_mass: f64,
    pub hidden: bool,
    pub classification: String,
    pub obo_accession: Option<String>,
    pub synonyms: BTreeSet<String>,
    pub neutral_losses: Vec<NeutralLoss>,
}
impl Default for ModificationRecord {
    fn default() -> Self {
        Self {
            record_id: None,
            name: String::new(),
            full_name: String::new(),
            full_id: String::new(),
            origin: Some('X'),
            term_specificity: TermSpecificity::Anywhere,
            diff_formula: EmpiricalFormula::default(),
            absolute_formula: None,
            diff_mono_mass: 0.0,
            diff_average_mass: 0.0,
            mono_mass: 0.0,
            average_mass: 0.0,
            hidden: false,
            classification: "Artefact".into(),
            obo_accession: None,
            synonyms: BTreeSet::new(),
            neutral_losses: Vec::new(),
        }
    }
}

pub(super) fn full_identifier(name: &str, origin: Option<char>, term: TermSpecificity) -> String {
    let specificity = if term == TermSpecificity::Anywhere {
        origin.unwrap_or('X').to_string()
    } else if let Some(origin) = origin.filter(|&r| r != 'X') {
        format!("{} {origin}", term.name())
    } else {
        term.name().into()
    };
    format!("{name} ({specificity})")
}

/// A chemically described modification at one residue/terminal specificity.
///
/// Entries are immutable. One UniMod ID may have many specificity records.
#[derive(Clone, Debug, PartialEq)]
pub struct ResidueModification {
    record_id: Option<u32>,
    obo_accession: Option<String>,
    synonyms: BTreeSet<String>,
    absolute_formula: Option<EmpiricalFormula>,
    name: String,
    full_name: String,
    full_id: String,
    origin: Option<char>,
    term: TermSpecificity,
    diff_formula: EmpiricalFormula,
    diff_mono_mass: f64,
    diff_average_mass: f64,
    mono_mass: f64,
    average_mass: f64,
    hidden: bool,
    classification: String,
    neutral_losses: Vec<NeutralLoss>,
}
// The registry rejects NaN/Inf values before an entry is constructed.
impl Eq for ResidueModification {}

/// Complete value ordering, not identifier-only or pointer ordering. Every
/// stored identity and chemistry field participates, consistently with Eq.
impl Ord for ResidueModification {
    fn cmp(&self, other: &Self) -> Ordering {
        // Peptide clones usually share this exact immutable Arc allocation.
        if std::ptr::eq(self, other) {
            return Ordering::Equal;
        }
        self.record_id
            .cmp(&other.record_id)
            .then_with(|| self.obo_accession.cmp(&other.obo_accession))
            .then_with(|| self.synonyms.cmp(&other.synonyms))
            .then_with(|| {
                self.absolute_formula
                    .as_ref()
                    .map(|f| (&f.atoms, f.charge))
                    .cmp(
                        &other
                            .absolute_formula
                            .as_ref()
                            .map(|f| (&f.atoms, f.charge)),
                    )
            })
            .then_with(|| self.name.cmp(&other.name))
            .then_with(|| self.full_name.cmp(&other.full_name))
            .then_with(|| self.full_id.cmp(&other.full_id))
            .then_with(|| self.origin.cmp(&other.origin))
            .then_with(|| self.term.cmp(&other.term))
            .then_with(|| formula_order(&self.diff_formula, &other.diff_formula))
            .then_with(|| finite_mass_order(self.diff_mono_mass, other.diff_mono_mass))
            .then_with(|| finite_mass_order(self.diff_average_mass, other.diff_average_mass))
            .then_with(|| finite_mass_order(self.mono_mass, other.mono_mass))
            .then_with(|| finite_mass_order(self.average_mass, other.average_mass))
            .then_with(|| self.hidden.cmp(&other.hidden))
            .then_with(|| self.classification.cmp(&other.classification))
            .then_with(|| self.neutral_losses.cmp(&other.neutral_losses))
    }
}
impl PartialOrd for ResidueModification {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl ResidueModification {
    pub fn from_record(mut record: ModificationRecord) -> Result<Self> {
        if [
            record.diff_mono_mass,
            record.diff_average_mass,
            record.mono_mass,
            record.average_mass,
        ]
        .iter()
        .any(|mass| !mass.is_finite())
        {
            return Err(Error::InvalidValue(
                "modification masses must be finite".into(),
            ));
        }
        if (record.name.is_empty() && record.full_id.is_empty())
            || [
                &record.name,
                &record.full_name,
                &record.full_id,
                &record.classification,
            ]
            .into_iter()
            .chain(record.obo_accession.iter())
            .chain(record.synonyms.iter())
            .any(|value| value.chars().any(char::is_control))
            || record.obo_accession.as_ref().is_some_and(String::is_empty)
            || record.synonyms.iter().any(String::is_empty)
            || record.origin.is_some_and(|r| !r.is_ascii_uppercase())
        {
            return Err(Error::InvalidValue(
                "invalid modification identity or origin".into(),
            ));
        }
        if record.term_specificity == TermSpecificity::Anywhere && record.origin.is_none() {
            record.origin = Some('X');
        } else if record.term_specificity != TermSpecificity::Anywhere && record.origin == Some('X')
        {
            record.origin = None;
        }
        if record.full_id.is_empty() {
            record.full_id = full_identifier(&record.name, record.origin, record.term_specificity);
        }
        Ok(Self {
            record_id: record.record_id,
            name: record.name,
            full_name: record.full_name,
            full_id: record.full_id,
            origin: record.origin,
            term: record.term_specificity,
            diff_formula: record.diff_formula,
            absolute_formula: record.absolute_formula,
            diff_mono_mass: record.diff_mono_mass,
            diff_average_mass: record.diff_average_mass,
            mono_mass: record.mono_mass,
            average_mass: record.average_mass,
            hidden: record.hidden,
            classification: record.classification,
            obo_accession: record.obo_accession,
            synonyms: record.synonyms,
            neutral_losses: record.neutral_losses,
        })
    }

    pub fn record_id(&self) -> Option<u32> {
        self.record_id
    }
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn full_name(&self) -> &str {
        &self.full_name
    }
    pub fn full_id(&self) -> &str {
        &self.full_id
    }
    pub fn accession(&self) -> String {
        self.unimod_accession()
            .or_else(|| self.obo_accession.clone())
            .unwrap_or_default()
    }
    pub fn unimod_accession(&self) -> Option<String> {
        self.record_id.map(|id| format!("UniMod:{id}"))
    }
    pub fn obo_accession(&self) -> Option<&str> {
        self.obo_accession.as_deref()
    }
    pub fn synonyms(&self) -> &BTreeSet<String> {
        &self.synonyms
    }
    /// Declared absolute FREE-residue composition, distinct from delta formula.
    pub fn absolute_formula(&self) -> Option<&EmpiricalFormula> {
        self.absolute_formula.as_ref()
    }
    /// None denotes a terminal modification that accepts any residue.
    pub fn origin(&self) -> Option<char> {
        self.origin
    }
    pub fn term_specificity(&self) -> TermSpecificity {
        self.term
    }
    pub fn diff_formula(&self) -> &EmpiricalFormula {
        &self.diff_formula
    }
    pub fn diff_mono_mass(&self) -> f64 {
        self.diff_mono_mass
    }
    pub fn diff_average_mass(&self) -> f64 {
        self.diff_average_mass
    }
    /// Declared absolute modified-residue mass; zero in the bundled UniMod and
    /// custom table records. Nonpositive values request residue-based fallback
    /// in ModificationDefinitionsSet's absolute matching when a residue is given.
    pub fn mono_mass(&self) -> f64 {
        self.mono_mass
    }
    pub fn average_mass(&self) -> f64 {
        self.average_mass
    }
    /// Set explicit absolute masses on an owned record without changing its
    /// formula, delta masses, or the immutable registry. Finite nonpositive
    /// values retain the source's unset/fallback convention.
    pub fn with_absolute_masses(mut self, mono_mass: f64, average_mass: f64) -> Result<Self> {
        if !mono_mass.is_finite() || !average_mass.is_finite() {
            return Err(Error::InvalidValue(
                "absolute modification masses must be finite".into(),
            ));
        }
        self.mono_mass = mono_mass;
        self.average_mass = average_mass;
        Ok(self)
    }
    pub fn is_hidden(&self) -> bool {
        self.hidden
    }
    pub fn classification(&self) -> &str {
        &self.classification
    }
    pub fn neutral_losses(&self) -> &[NeutralLoss] {
        &self.neutral_losses
    }

    fn matches(&self, residue: Option<char>, term: Option<TermSpecificity>) -> bool {
        term.is_none_or(|t| t == self.term)
            && residue.is_none_or(|r| self.origin.is_none_or(|o| o == r || o == 'X'))
    }
}

/// Registry of immutable, shared records. The global registry is immutable;
/// caller-owned registries support checked atomic provider appends.
#[derive(Clone, Debug, Default)]
pub struct ModificationsDB {
    entries: Vec<Arc<ResidueModification>>,
    by_name: BTreeMap<String, Vec<usize>>,
}

impl ModificationsDB {
    /// Shared registry initialized once from validated embedded lookup data.
    pub fn global() -> &'static Self {
        static INSTANCE: OnceLock<ModificationsDB> = OnceLock::new();
        INSTANCE.get_or_init(|| {
            let mut db = Self::from_tsv(include_str!(
                "../../resources/modifications/openms-rust-modifications.tsv"
            ))
            .expect("bundled modification table must pass its checked parser");
            db.extend_obo(
                include_bytes!("../../resources/modifications/XLMOD.obo").as_slice(),
                &OboReadOptions::default(),
            )
            .expect("bundled XLMOD must pass its checked parser");
            db
        })
    }

    /// Parse the documented generated-table format, for custom native registries.
    /// This does not replace the global registry used by AASequence parsing.
    pub fn from_tsv(text: &str) -> Result<Self> {
        let mut entries = Vec::new();
        let mut by_name: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (line, row) in text.lines().enumerate() {
            if row.starts_with('#') || row.is_empty() {
                continue;
            }
            let fields: Vec<&str> = row.split('\t').collect();
            let error = |message: &str| Error::Parse {
                line: line + 1,
                message: message.into(),
            };
            if fields.len() != 11 {
                return Err(error("modification table requires 11 fields"));
            }
            let mass = |text: &str| -> Result<f64> {
                text.parse::<f64>()
                    .ok()
                    .filter(|v| v.is_finite())
                    .ok_or_else(|| error("modification mass must be finite"))
            };
            let term = match fields[4] {
                "anywhere" => TermSpecificity::Anywhere,
                "n-term" => TermSpecificity::NTerm,
                "c-term" => TermSpecificity::CTerm,
                "protein-n-term" => TermSpecificity::ProteinNTerm,
                "protein-c-term" => TermSpecificity::ProteinCTerm,
                _ => return Err(error("unknown terminal specificity")),
            };
            // The pinned XML contains terminal site labels with Anywhere
            // position. OpenMS treats those as the wildcard origin X.
            let origin = if matches!(fields[3], "N-term" | "C-term") {
                if term == TermSpecificity::Anywhere {
                    Some('X')
                } else {
                    None
                }
            } else if fields[3].len() == 1 && fields[3].as_bytes()[0].is_ascii_uppercase() {
                fields[3].chars().next()
            } else {
                return Err(error("invalid modification origin"));
            };
            if fields[1].is_empty()
                || fields[1].chars().any(char::is_control)
                || fields[2].chars().any(char::is_control)
            {
                return Err(error("invalid modification name"));
            }
            let specificity = if term == TermSpecificity::Anywhere {
                origin.unwrap().to_string()
            } else if let Some(origin) = origin {
                format!("{} {origin}", term.name())
            } else {
                term.name().into()
            };
            let mut neutral_losses = Vec::new();
            if !fields[10].is_empty() {
                for loss in fields[10].split(';') {
                    let parts: Vec<&str> = loss.split('@').collect();
                    if parts.len() != 3 {
                        return Err(error("invalid neutral loss"));
                    }
                    neutral_losses.push(NeutralLoss {
                        formula: EmpiricalFormula::parse(parts[0])?,
                        mono_mass: mass(parts[1])?,
                        average_mass: mass(parts[2])?,
                    });
                }
            }
            let modification = ResidueModification {
                record_id: Some(
                    fields[0]
                        .parse()
                        .map_err(|_| error("invalid modification record ID"))?,
                ),
                obo_accession: None,
                synonyms: BTreeSet::new(),
                absolute_formula: None,
                name: fields[1].into(),
                full_name: fields[2].into(),
                full_id: format!("{} ({specificity})", fields[1]),
                origin,
                term,
                diff_formula: EmpiricalFormula::parse(fields[7])?,
                diff_mono_mass: mass(fields[5])?,
                diff_average_mass: mass(fields[6])?,
                mono_mass: 0.0,
                average_mass: 0.0,
                hidden: match fields[8] {
                    "0" => false,
                    "1" => true,
                    _ => return Err(error("invalid hidden flag")),
                },
                classification: fields[9].into(),
                neutral_losses,
            };
            let names = BTreeSet::from([
                modification.name.clone(),
                modification.full_name.clone(),
                modification.full_id.clone(),
                modification.accession(),
            ]);
            for name in names {
                by_name.entry(name).or_default().push(entries.len());
            }
            entries.push(Arc::new(modification));
        }
        Ok(Self { entries, by_name })
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    pub fn entries(&self) -> &[Arc<ResidueModification>] {
        &self.entries
    }

    /// Find matching short/full names, full IDs or UniMod accessions.
    /// The UniMod prefix is case-insensitive; names otherwise match exactly.
    pub fn find(
        &self,
        name: &str,
        residue: Option<char>,
        term: Option<TermSpecificity>,
    ) -> Vec<&ResidueModification> {
        let normalized;
        let name = if name
            .get(..7)
            .is_some_and(|p| p.eq_ignore_ascii_case("unimod:"))
        {
            normalized = format!("UniMod:{}", &name[7..]);
            normalized.as_str()
        } else {
            name
        };
        self.by_name
            .get(name)
            .into_iter()
            .flatten()
            .map(|&i| self.entries[i].as_ref())
            .filter(|m| m.matches(residue, term))
            .collect()
    }

    /// Resolve one specificity. Ambiguous distinct full IDs are errors.
    /// Duplicate database specificity records with the same full ID resolve to
    /// the first stored entry, preserving the pinned source-file order.
    pub fn get_modification(
        &self,
        name: &str,
        residue: Option<char>,
        term: Option<TermSpecificity>,
    ) -> Result<&ResidueModification> {
        let matches = self.find(name, residue, term);
        let Some(first) = matches.first().copied() else {
            return Err(Error::InvalidValue(format!(
                "no modification {name:?} matches the requested residue/terminus"
            )));
        };
        if matches.iter().any(|m| m.full_id != first.full_id) {
            return Err(Error::InvalidValue(format!(
                "ambiguous modification {name:?}; specify residue and terminal specificity"
            )));
        }
        Ok(first)
    }

    /// Match declared monoisotopic mass differences in an inclusive Da window.
    pub fn search_by_mass(
        &self,
        mass: f64,
        tolerance: f64,
        residue: Option<char>,
        term: Option<TermSpecificity>,
    ) -> Result<Vec<&ResidueModification>> {
        if !mass.is_finite() || !tolerance.is_finite() || tolerance < 0.0 {
            return Err(Error::InvalidValue(
                "mass and nonnegative tolerance must be finite".into(),
            ));
        }
        Ok(self
            .entries
            .iter()
            .map(Arc::as_ref)
            .filter(|m| {
                m.matches(residue, term)
                    && (m.diff_mono_mass != 0.0 || mass == 0.0)
                    && (m.diff_mono_mass - mass).abs() <= tolerance
            })
            .collect())
    }

    /// Closest mass match, with source-order tie breaking. As in OpenMS's
    /// best-match API, the mass error must be strictly less than tolerance;
    /// tolerance zero never returns a hit. The all-matches search is inclusive.
    pub fn best_by_mass(
        &self,
        mass: f64,
        tolerance: f64,
        residue: Option<char>,
        term: Option<TermSpecificity>,
    ) -> Result<Option<&ResidueModification>> {
        Ok(self
            .search_by_mass(mass, tolerance, residue, term)?
            .into_iter()
            .filter(|m| (m.diff_mono_mass - mass).abs() < tolerance)
            .min_by(|a, b| {
                (a.diff_mono_mass - mass)
                    .abs()
                    .total_cmp(&(b.diff_mono_mass - mass).abs())
            }))
    }
}

impl ModificationsDB {
    pub fn from_records(records: Vec<ResidueModification>) -> Result<Self> {
        let mut database = Self::default();
        database.extend_records(records)?;
        Ok(database)
    }

    /// Append validated owned records in input order. Existing handles stay valid.
    /// Default OBO limits also bound the complete registry assembled this way.
    pub fn extend_records(&mut self, records: Vec<ResidueModification>) -> Result<()> {
        self.extend_loaded(records, &OboReadOptions::default(), false)?;
        Ok(())
    }

    pub fn from_obo(reader: impl BufRead, options: &OboReadOptions) -> Result<Self> {
        let mut database = Self::default();
        database.extend_obo(reader, options)?;
        Ok(database)
    }

    /// Parse an OBO provider completely before appending it. UniMod-linked OBO
    /// records register accession aliases to every existing UniMod specificity;
    /// absent targets are reported and omitted, as in OpenMS provider loading.
    pub fn extend_obo(
        &mut self,
        reader: impl BufRead,
        options: &OboReadOptions,
    ) -> Result<OboLoadReport> {
        let records = obo::read(reader, options)?;
        self.extend_loaded(records, options, true)
    }

    fn extend_loaded(
        &mut self,
        records: Vec<ResidueModification>,
        options: &OboReadOptions,
        resolve_aliases: bool,
    ) -> Result<OboLoadReport> {
        options.validate()?;
        let mut report = OboLoadReport::default();
        // Plan all indices and check both existing and new payload before cloning
        // the index. Arc clones do not duplicate immutable chemical records.
        let mut aliases: BTreeMap<String, BTreeSet<usize>> = BTreeMap::new();
        let mut additions = Vec::new();
        let mut alias_visits = 0;
        let mut total_aliases = self.by_name.values().try_fold(0usize, |n, v| {
            obo::checked(n, v.len(), options.max_aliases, "registry aliases")
        })?;
        let mut bytes = self.entries.iter().try_fold(0usize, |n, m| {
            obo::checked(
                n,
                m.payload_bytes()?,
                options.max_registry_bytes,
                "registry bytes",
            )
        })?;
        for (name, indices) in &self.by_name {
            bytes = obo::checked(
                bytes,
                name.len()
                    .saturating_add(indices.len().saturating_mul(size_of::<usize>()))
                    .saturating_add(64),
                options.max_registry_bytes,
                "registry bytes",
            )?;
        }
        obo::checked(
            self.entries.len(),
            0,
            options.max_records,
            "registry records",
        )?;
        for record in records {
            let alias = record
                .record_id
                .filter(|&id| id > 0)
                .zip(record.obo_accession.as_ref())
                .filter(|_| resolve_aliases);
            if let Some((id, accession)) = alias {
                let target = format!("UniMod:{id}");
                if let Some(indices) = self.by_name.get(&target) {
                    for &index in indices {
                        alias_visits = obo::checked(
                            alias_visits,
                            1,
                            options.max_aliases,
                            "OBO alias-target visits",
                        )?;
                        if self
                            .by_name
                            .get(accession)
                            .is_some_and(|v| v.binary_search(&index).is_ok())
                            || aliases
                                .get(accession)
                                .is_some_and(|indices| indices.contains(&index))
                        {
                            continue;
                        }
                        total_aliases = obo::checked(
                            total_aliases,
                            1,
                            options.max_aliases,
                            "registry aliases",
                        )?;
                        bytes = obo::checked(
                            bytes,
                            accession.len().saturating_add(72),
                            options.max_registry_bytes,
                            "registry bytes",
                        )?;
                        aliases.entry(accession.clone()).or_default().insert(index);
                        report.aliases_added += 1;
                    }
                } else {
                    report.unresolved_aliases += 1;
                }
                continue;
            }
            let index = obo::checked(
                self.entries.len(),
                additions.len(),
                options.max_records,
                "registry records",
            )?;
            obo::checked(index, 1, options.max_records, "registry records")?;
            bytes = obo::checked(
                bytes,
                record.payload_bytes()?,
                options.max_registry_bytes,
                "registry bytes",
            )?;
            for name in record.names() {
                total_aliases =
                    obo::checked(total_aliases, 1, options.max_aliases, "registry aliases")?;
                bytes = obo::checked(
                    bytes,
                    name.len().saturating_add(72),
                    options.max_registry_bytes,
                    "registry bytes",
                )?;
                aliases.entry(name).or_default().insert(index);
            }
            additions.push(Arc::new(record));
        }
        report.records_added = additions.len();
        let mut next = self.clone();
        next.entries.extend(additions);
        for (name, indices) in aliases {
            next.by_name.entry(name).or_default().extend(indices);
        }
        for indices in next.by_name.values_mut() {
            indices.sort_unstable();
        }
        *self = next;
        Ok(report)
    }

    pub fn find_handles(
        &self,
        name: &str,
        residue: Option<char>,
        term: Option<TermSpecificity>,
    ) -> Vec<Arc<ResidueModification>> {
        let normalized;
        let name = if name
            .get(..7)
            .is_some_and(|p| p.eq_ignore_ascii_case("unimod:"))
        {
            normalized = format!("UniMod:{}", &name[7..]);
            normalized.as_str()
        } else {
            name
        };
        self.by_name
            .get(name)
            .into_iter()
            .flatten()
            .map(|&i| &self.entries[i])
            .filter(|m| m.matches(residue, term))
            .cloned()
            .collect()
    }

    pub fn get_modification_handle(
        &self,
        name: &str,
        residue: Option<char>,
        term: Option<TermSpecificity>,
    ) -> Result<Arc<ResidueModification>> {
        // Reuse the borrowed API's checked ambiguity contract before cloning.
        self.get_modification(name, residue, term)?;
        Ok(self.find_handles(name, residue, term).remove(0))
    }

    pub fn search_by_mass_handles(
        &self,
        mass: f64,
        tolerance: f64,
        residue: Option<char>,
        term: Option<TermSpecificity>,
    ) -> Result<Vec<Arc<ResidueModification>>> {
        if !mass.is_finite() || !tolerance.is_finite() || tolerance < 0.0 {
            return Err(Error::InvalidValue(
                "mass and nonnegative tolerance must be finite".into(),
            ));
        }
        Ok(self
            .entries
            .iter()
            .filter(|m| {
                m.matches(residue, term)
                    && (m.diff_mono_mass != 0.0 || mass == 0.0)
                    && (m.diff_mono_mass - mass).abs() <= tolerance
            })
            .cloned()
            .collect())
    }

    pub fn best_by_mass_handle(
        &self,
        mass: f64,
        tolerance: f64,
        residue: Option<char>,
        term: Option<TermSpecificity>,
    ) -> Result<Option<Arc<ResidueModification>>> {
        Ok(self
            .search_by_mass_handles(mass, tolerance, residue, term)?
            .into_iter()
            .filter(|m| (m.diff_mono_mass - mass).abs() < tolerance)
            .min_by(|a, b| {
                (a.diff_mono_mass - mass)
                    .abs()
                    .total_cmp(&(b.diff_mono_mass - mass).abs())
            }))
    }
}

impl ResidueModification {
    fn names(&self) -> BTreeSet<String> {
        let mut names = self.synonyms.clone();
        names.extend([
            self.name.clone(),
            self.full_name.clone(),
            self.full_id.clone(),
        ]);
        names.extend(self.unimod_accession());
        names.extend(self.obo_accession.clone());
        names.remove("");
        names
    }
    fn payload_bytes(&self) -> Result<usize> {
        // Conservative payload accounting: strings plus map/atom/loss storage.
        let mut bytes = size_of::<Self>();
        for value in [
            &self.name,
            &self.full_name,
            &self.full_id,
            &self.classification,
        ]
        .into_iter()
        .chain(self.obo_accession.iter())
        .chain(self.synonyms.iter())
        {
            bytes = obo::checked(
                bytes,
                value.len().saturating_add(64),
                usize::MAX,
                "registry bytes",
            )?;
        }
        for formula in std::iter::once(&self.diff_formula)
            .chain(self.absolute_formula.iter())
            .chain(self.neutral_losses.iter().map(NeutralLoss::formula))
        {
            bytes = obo::checked(
                bytes,
                formula
                    .to_string()
                    .len()
                    .saturating_mul(64)
                    .saturating_add(64),
                usize::MAX,
                "registry bytes",
            )?;
        }
        obo::checked(
            bytes,
            self.neutral_losses
                .len()
                .saturating_mul(size_of::<NeutralLoss>()),
            usize::MAX,
            "registry bytes",
        )
    }
}
