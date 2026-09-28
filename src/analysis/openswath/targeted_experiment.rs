// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A targeted experiment: the in-memory form of a TraML document.
//!
//! Ports `ANALYSIS/TARGETED/TargetedExperiment.h` (with its `.cpp`). A
//! [`TargetedExperiment`] describes the precursors and product ions analysed
//! in an SRM/MRM experiment, or in a SWATH-MS/DIA experiment analysed with a
//! targeted approach. Precursors are peptides or small molecules
//! (metabolomics), each with a set of product ions, the transitions. Around
//! them it keeps the document's controlled vocabularies, contacts,
//! publications, instruments, software, proteins, inclusion/exclusion
//! targets and source files.
//!
//! Upstream stores it as `.traml` with `TraMLFile` or as `.tsv` with
//! `TransitionTSVFile`; neither format is ported yet, and both are backlog
//! items that depend on this type.
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/TARGETED_EXPERIMENT_SUPPORT.md)
//! for the member mapping and the evidence. The source contains no
//! `#pragma omp`; the port is serial.

use crate::analysis::openswath::include_exclude_target::IncludeExcludeTarget;
use crate::analysis::openswath::reaction_monitoring_transition::{
    DecoyTransitionType, ReactionMonitoringTransition,
};
use crate::analysis::openswath::targeted_experiment_helper::{
    CV, Compound, Contact, Instrument, Peptide, Protein, Publication,
};
use crate::metadata::{CVTerm, CVTermList, MetaValue, Software, SourceFile};
use crate::{Error, Result};
use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::sync::OnceLock;

/// Counts summarising a [`TargetedExperiment`], from
/// [`TargetedExperiment::summary`].
///
/// [`Display`](fmt::Display) prints the source `operator<<` format.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SummaryStatistics {
    /// Number of proteins.
    pub protein_count: usize,
    /// Number of peptides.
    pub peptide_count: usize,
    /// Number of compounds.
    pub compound_count: usize,
    /// Number of transitions.
    pub transition_count: usize,
    /// Number of transitions per target/decoy type; a type with no
    /// transition has no entry.
    pub decoy_counts: BTreeMap<DecoyTransitionType, usize>,
    /// Result of [`TargetedExperiment::contains_invalid_references`].
    pub contains_invalid_references: bool,
}

fn format_count(count: usize, all: usize, name: &str, sink: &mut Vec<String>) {
    if count == 0 {
        return;
    }
    // Source: StringUtils::toStr(count * 100.0 / all, false).
    let percent = count as f64 * 100.0 / all as f64;
    sink.push(format!(
        "{}% ({name})",
        crate::param::value::format_float(percent, false)
    ));
}

impl fmt::Display for SummaryStatistics {
    /// The source `operator<<(ostream&, const SummaryStatistics&)`: six
    /// lines, each ending in a newline. The transition-type line lists the
    /// target, decoy and unknown shares in that order, skipping zero counts,
    /// with percentages formatted as `StringUtils::toStr(double, false)`
    /// (three fractional digits, trailing zeros trimmed).
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let count = |t| self.decoy_counts.get(&t).copied().unwrap_or(0);
        let (target, decoy, unknown) = (
            count(DecoyTransitionType::Target),
            count(DecoyTransitionType::Decoy),
            count(DecoyTransitionType::Unknown),
        );
        let mut all = decoy.saturating_add(target).saturating_add(unknown);
        if all == 0 {
            all = 1;
        }
        let mut counts = Vec::new();
        format_count(target, all, "target", &mut counts);
        format_count(decoy, all, "decoy", &mut counts);
        format_count(unknown, all, "unknown", &mut counts);
        writeln!(f, "# Proteins: {}", self.protein_count)?;
        writeln!(f, "# Peptides: {}", self.peptide_count)?;
        writeln!(f, "# Compounds: {}", self.compound_count)?;
        writeln!(f, "# Transitions: {}", self.transition_count)?;
        writeln!(f, "Transition Type: {}", counts.join(", "))?;
        writeln!(
            f,
            "All internal references valid: {}",
            if self.contains_invalid_references {
                "no"
            } else {
                "yes"
            }
        )
    }
}

/// Lazily built id -> index maps. The source keeps three mutable pointer maps
/// with dirty flags; here a mutation replaces the cell, so a lookup can never
/// see a stale entry.
#[derive(Clone, Debug, Default)]
struct ReferenceMaps {
    proteins: OnceLock<HashMap<String, usize>>,
    peptides: OnceLock<HashMap<String, usize>>,
    compounds: OnceLock<HashMap<String, usize>>,
}

/// Index items by id; a repeated id maps to its **last** occurrence, as the
/// source's `map[id] = &item` loop leaves it.
fn index_by_id<'a>(ids: impl Iterator<Item = &'a str>) -> HashMap<String, usize> {
    let mut map = HashMap::new();
    for (index, id) in ids.enumerate() {
        map.insert(id.to_owned(), index);
    }
    map
}

/// A targeted experiment: transitions, their analytes and the document
/// metadata around them.
///
/// All lists keep insertion order. Lookup by reference
/// ([`protein_by_ref`](Self::protein_by_ref) and friends) builds an index on
/// first use; every mutation of the list invalidates it.
#[derive(Clone, Debug, Default)]
pub struct TargetedExperiment {
    cvs: Vec<CV>,
    contacts: Vec<Contact>,
    publications: Vec<Publication>,
    instruments: Vec<Instrument>,
    targets: CVTermList,
    software: Vec<Software>,
    proteins: Vec<Protein>,
    compounds: Vec<Compound>,
    peptides: Vec<Peptide>,
    transitions: Vec<ReactionMonitoringTransition>,
    include_targets: Vec<IncludeExcludeTarget>,
    exclude_targets: Vec<IncludeExcludeTarget>,
    source_files: Vec<SourceFile>,
    maps: ReferenceMaps,
}

/// Equality over every list and the target CV terms, as the source; the
/// lookup indices do not take part.
impl PartialEq for TargetedExperiment {
    fn eq(&self, rhs: &Self) -> bool {
        self.cvs == rhs.cvs
            && self.contacts == rhs.contacts
            && self.publications == rhs.publications
            && self.instruments == rhs.instruments
            && self.targets == rhs.targets
            && self.software == rhs.software
            && self.proteins == rhs.proteins
            && self.compounds == rhs.compounds
            && self.peptides == rhs.peptides
            && self.transitions == rhs.transitions
            && self.include_targets == rhs.include_targets
            && self.exclude_targets == rhs.exclude_targets
            && self.source_files == rhs.source_files
    }
}

impl TargetedExperiment {
    /// An empty experiment.
    pub fn new() -> Self {
        Self::default()
    }

    /// Join two experiments into a new one; the source `operator+`. See
    /// [`merge`](Self::merge).
    ///
    /// # Errors
    ///
    /// As [`merge`](Self::merge).
    pub fn merged(&self, rhs: &Self) -> Result<Self> {
        let mut out = self.clone();
        out.merge(rhs)?;
        Ok(out)
    }

    /// Append another experiment to this one; the source `operator+=`.
    ///
    /// Every list of `rhs` is appended to the corresponding list here, and
    /// each target CV term of `rhs` is added to the target CV terms. The
    /// target list's ordinary metadata is **not** merged, as in the source.
    /// Nothing is deduplicated: the source leaves "check for double entries"
    /// as a to-do, so appending an experiment to itself duplicates every id,
    /// which [`contains_invalid_references`](Self::contains_invalid_references)
    /// then reports.
    ///
    /// # Errors
    ///
    /// Whatever [`CVTermList::add`] returns for an invalid target CV term.
    /// Terms can only enter a [`CVTermList`] through validating calls, so
    /// this does not occur in practice; if it did, `self` would be unchanged,
    /// because the new target list is built before anything is appended.
    pub fn merge(&mut self, rhs: &Self) -> Result<()> {
        let targets = merged_targets(&self.targets, &rhs.targets)?;
        self.targets = targets;
        self.cvs.extend_from_slice(&rhs.cvs);
        self.contacts.extend_from_slice(&rhs.contacts);
        self.publications.extend_from_slice(&rhs.publications);
        self.instruments.extend_from_slice(&rhs.instruments);
        self.software.extend_from_slice(&rhs.software);
        self.proteins.extend_from_slice(&rhs.proteins);
        self.compounds.extend_from_slice(&rhs.compounds);
        self.peptides.extend_from_slice(&rhs.peptides);
        self.transitions.extend_from_slice(&rhs.transitions);
        self.include_targets.extend_from_slice(&rhs.include_targets);
        self.exclude_targets.extend_from_slice(&rhs.exclude_targets);
        self.source_files.extend_from_slice(&rhs.source_files);
        self.maps = ReferenceMaps::default();
        Ok(())
    }

    /// [`merge`](Self::merge), consuming `rhs` so its lists are moved rather
    /// than copied; the source's rvalue `operator+=`.
    ///
    /// # Errors
    ///
    /// As [`merge`](Self::merge); `self` is unchanged on error.
    pub fn merge_owned(&mut self, rhs: Self) -> Result<()> {
        let targets = merged_targets(&self.targets, &rhs.targets)?;
        self.targets = targets;
        let Self {
            cvs,
            contacts,
            publications,
            instruments,
            targets: _,
            software,
            proteins,
            compounds,
            peptides,
            transitions,
            include_targets,
            exclude_targets,
            source_files,
            maps: _,
        } = rhs;
        self.cvs.extend(cvs);
        self.contacts.extend(contacts);
        self.publications.extend(publications);
        self.instruments.extend(instruments);
        self.software.extend(software);
        self.proteins.extend(proteins);
        self.compounds.extend(compounds);
        self.peptides.extend(peptides);
        self.transitions.extend(transitions);
        self.include_targets.extend(include_targets);
        self.exclude_targets.extend(exclude_targets);
        self.source_files.extend(source_files);
        self.maps = ReferenceMaps::default();
        Ok(())
    }

    /// Clear the transitions and, when `clear_meta_data` is true, everything
    /// else as well (every list, and the target CV terms including their
    /// metadata).
    ///
    /// With `clear_meta_data` false only the transitions go; proteins,
    /// peptides and compounds stay, as in the source.
    pub fn clear(&mut self, clear_meta_data: bool) {
        self.transitions.clear();
        if clear_meta_data {
            self.cvs.clear();
            self.contacts.clear();
            self.publications.clear();
            self.instruments.clear();
            self.targets = CVTermList::default();
            self.software.clear();
            self.proteins.clear();
            self.compounds.clear();
            self.peptides.clear();
            self.include_targets.clear();
            self.exclude_targets.clear();
            self.source_files.clear();
            self.maps = ReferenceMaps::default();
        }
    }

    /// Summary statistics: list sizes, transitions per decoy type and the
    /// reference check.
    pub fn summary(&self) -> SummaryStatistics {
        let mut decoy_counts = BTreeMap::new();
        for transition in &self.transitions {
            *decoy_counts
                .entry(transition.decoy_transition_type())
                .or_insert(0) += 1;
        }
        SummaryStatistics {
            protein_count: self.proteins.len(),
            peptide_count: self.peptides.len(),
            compound_count: self.compounds.len(),
            transition_count: self.transitions.len(),
            decoy_counts,
            contains_invalid_references: self.contains_invalid_references(),
        }
    }

    /// Replace the controlled vocabularies.
    pub fn set_cvs(&mut self, cvs: Vec<CV>) {
        self.cvs = cvs;
    }

    /// The controlled vocabularies.
    pub fn cvs(&self) -> &[CV] {
        &self.cvs
    }

    /// Append a controlled vocabulary.
    pub fn add_cv(&mut self, cv: CV) {
        self.cvs.push(cv);
    }

    /// Replace the contacts.
    pub fn set_contacts(&mut self, contacts: Vec<Contact>) {
        self.contacts = contacts;
    }

    /// The contacts.
    pub fn contacts(&self) -> &[Contact] {
        &self.contacts
    }

    /// Append a contact.
    pub fn add_contact(&mut self, contact: Contact) {
        self.contacts.push(contact);
    }

    /// Replace the publications.
    pub fn set_publications(&mut self, publications: Vec<Publication>) {
        self.publications = publications;
    }

    /// The publications.
    pub fn publications(&self) -> &[Publication] {
        &self.publications
    }

    /// Append a publication.
    pub fn add_publication(&mut self, publication: Publication) {
        self.publications.push(publication);
    }

    /// Replace the target-list CV terms (TraML `<TargetList>`), including
    /// their metadata.
    pub fn set_target_cv_terms(&mut self, cv_terms: CVTermList) {
        self.targets = cv_terms;
    }

    /// The target-list CV terms.
    pub fn target_cv_terms(&self) -> &CVTermList {
        &self.targets
    }

    /// Add a target-list CV term.
    ///
    /// # Errors
    ///
    /// Whatever [`CVTermList::add`] returns for an invalid term; the list is
    /// then unchanged.
    pub fn add_target_cv_term(&mut self, cv_term: CVTerm) -> Result<()> {
        self.targets.add(cv_term)
    }

    /// Set a metadata value (a TraML `userParam`) on the target list,
    /// replacing any value under the same name.
    pub fn set_target_meta_value(&mut self, name: impl Into<String>, value: MetaValue) {
        self.targets.metadata.insert(name.into(), value);
    }

    /// Replace the instruments.
    pub fn set_instruments(&mut self, instruments: Vec<Instrument>) {
        self.instruments = instruments;
    }

    /// The instruments.
    pub fn instruments(&self) -> &[Instrument] {
        &self.instruments
    }

    /// Append an instrument.
    pub fn add_instrument(&mut self, instrument: Instrument) {
        self.instruments.push(instrument);
    }

    /// Replace the software list.
    pub fn set_software(&mut self, software: Vec<Software>) {
        self.software = software;
    }

    /// The software list.
    pub fn software(&self) -> &[Software] {
        &self.software
    }

    /// Append a software record.
    pub fn add_software(&mut self, software: Software) {
        self.software.push(software);
    }

    /// Replace the proteins; invalidates the protein index.
    pub fn set_proteins(&mut self, proteins: Vec<Protein>) {
        self.proteins = proteins;
        self.maps.proteins = OnceLock::new();
    }

    /// The proteins.
    pub fn proteins(&self) -> &[Protein] {
        &self.proteins
    }

    fn protein_index(&self) -> &HashMap<String, usize> {
        self.maps
            .proteins
            .get_or_init(|| index_by_id(self.proteins.iter().map(|p| p.id.as_str())))
    }

    /// The protein with id `reference`; with a repeated id, the last one.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] (the source's
    /// `Exception::ElementNotFound`) when no protein has that id.
    pub fn protein_by_ref(&self, reference: &str) -> Result<&Protein> {
        let index = self.protein_index().get(reference).copied();
        index
            .and_then(|i| self.proteins.get(i))
            .ok_or_else(|| not_found(reference))
    }

    /// Whether a protein with id `reference` exists.
    pub fn has_protein(&self, reference: &str) -> bool {
        self.protein_index().contains_key(reference)
    }

    /// Append a protein; invalidates the protein index.
    pub fn add_protein(&mut self, protein: Protein) {
        self.proteins.push(protein);
        self.maps.proteins = OnceLock::new();
    }

    /// Replace the compounds; invalidates the compound index.
    pub fn set_compounds(&mut self, compounds: Vec<Compound>) {
        self.compounds = compounds;
        self.maps.compounds = OnceLock::new();
    }

    /// The compounds.
    pub fn compounds(&self) -> &[Compound] {
        &self.compounds
    }

    /// Append a compound; invalidates the compound index.
    pub fn add_compound(&mut self, compound: Compound) {
        self.compounds.push(compound);
        self.maps.compounds = OnceLock::new();
    }

    fn compound_index(&self) -> &HashMap<String, usize> {
        self.maps
            .compounds
            .get_or_init(|| index_by_id(self.compounds.iter().map(|c| c.id.as_str())))
    }

    /// Whether a compound with id `reference` exists.
    pub fn has_compound(&self, reference: &str) -> bool {
        self.compound_index().contains_key(reference)
    }

    /// The compound with id `reference`; with a repeated id, the last one.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] (the source's
    /// `Exception::ElementNotFound`) when no compound has that id.
    pub fn compound_by_ref(&self, reference: &str) -> Result<&Compound> {
        let index = self.compound_index().get(reference).copied();
        index
            .and_then(|i| self.compounds.get(i))
            .ok_or_else(|| not_found(reference))
    }

    /// Replace the peptides; invalidates the peptide index.
    pub fn set_peptides(&mut self, peptides: Vec<Peptide>) {
        self.peptides = peptides;
        self.maps.peptides = OnceLock::new();
    }

    /// The peptides.
    pub fn peptides(&self) -> &[Peptide] {
        &self.peptides
    }

    fn peptide_index(&self) -> &HashMap<String, usize> {
        self.maps
            .peptides
            .get_or_init(|| index_by_id(self.peptides.iter().map(|p| p.id.as_str())))
    }

    /// Whether a peptide with id `reference` exists.
    pub fn has_peptide(&self, reference: &str) -> bool {
        self.peptide_index().contains_key(reference)
    }

    /// The peptide with id `reference`; with a repeated id, the last one.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] (the source's
    /// `Exception::ElementNotFound`) when no peptide has that id.
    pub fn peptide_by_ref(&self, reference: &str) -> Result<&Peptide> {
        let index = self.peptide_index().get(reference).copied();
        index
            .and_then(|i| self.peptides.get(i))
            .ok_or_else(|| not_found(reference))
    }

    /// Append a peptide; invalidates the peptide index.
    pub fn add_peptide(&mut self, peptide: Peptide) {
        self.peptides.push(peptide);
        self.maps.peptides = OnceLock::new();
    }

    /// Replace the transition list.
    pub fn set_transitions(&mut self, transitions: Vec<ReactionMonitoringTransition>) {
        self.transitions = transitions;
    }

    /// The transition list.
    pub fn transitions(&self) -> &[ReactionMonitoringTransition] {
        &self.transitions
    }

    /// Append a transition.
    pub fn add_transition(&mut self, transition: ReactionMonitoringTransition) {
        self.transitions.push(transition);
    }

    /// Replace the inclusion targets.
    pub fn set_include_targets(&mut self, targets: Vec<IncludeExcludeTarget>) {
        self.include_targets = targets;
    }

    /// The inclusion targets.
    pub fn include_targets(&self) -> &[IncludeExcludeTarget] {
        &self.include_targets
    }

    /// Append an inclusion target.
    pub fn add_include_target(&mut self, target: IncludeExcludeTarget) {
        self.include_targets.push(target);
    }

    /// Replace the exclusion targets.
    pub fn set_exclude_targets(&mut self, targets: Vec<IncludeExcludeTarget>) {
        self.exclude_targets = targets;
    }

    /// The exclusion targets.
    pub fn exclude_targets(&self) -> &[IncludeExcludeTarget] {
        &self.exclude_targets
    }

    /// Append an exclusion target.
    pub fn add_exclude_target(&mut self, target: IncludeExcludeTarget) {
        self.exclude_targets.push(target);
    }

    /// Replace the source files.
    pub fn set_source_files(&mut self, source_files: Vec<SourceFile>) {
        self.source_files = source_files;
    }

    /// The source files.
    pub fn source_files(&self) -> &[SourceFile] {
        &self.source_files
    }

    /// Append a source file.
    pub fn add_source_file(&mut self, source_file: SourceFile) {
        self.source_files.push(source_file);
    }

    /// Sort the transitions by product m/z, ascending.
    ///
    /// The source uses the unstable `std::sort` with `ProductMZLess`; this
    /// sort is stable, so transitions with equal product m/z keep their
    /// order (the source leaves it unspecified), and it orders by
    /// [`ReactionMonitoringTransition::cmp_product_mz`], whose total order
    /// places NaN deterministically where `<` would not.
    pub fn sort_transitions_by_product_mz(&mut self) {
        self.transitions
            .sort_by(ReactionMonitoringTransition::cmp_product_mz);
    }

    /// Sort the transitions by name, ascending and stable (see
    /// [`sort_transitions_by_product_mz`](Self::sort_transitions_by_product_mz)
    /// for the stability note).
    pub fn sort_transitions_by_name(&mut self) {
        self.transitions
            .sort_by(ReactionMonitoringTransition::cmp_name);
    }

    /// Whether the experiment (and the TraML file it came from) contains
    /// invalid references.
    ///
    /// Returns `false` for a valid experiment. See
    /// [`first_invalid_reference`](Self::first_invalid_reference) for the
    /// checks and their order.
    pub fn contains_invalid_references(&self) -> bool {
        self.first_invalid_reference().is_some()
    }

    /// The first problem [`contains_invalid_references`](Self::contains_invalid_references)
    /// finds, as the message the source logs, or `None` when every
    /// reference is valid.
    ///
    /// The checks run in the source order and stop at the first failure:
    /// protein ids unique, peptide ids unique, compound ids unique,
    /// transition ids unique, every protein reference of every peptide
    /// resolves, and every transition's peptide reference — or, when that is
    /// empty, its compound reference — resolves. A transition with neither
    /// reference is valid (both attributes are optional); the source logs a
    /// warning for it, which this port does not report.
    ///
    /// The source's compound message prints the (empty) peptide reference
    /// instead of the compound reference; this message names the compound.
    pub fn first_invalid_reference(&self) -> Option<String> {
        fn first_duplicate<'a>(
            what: &str,
            ids: impl Iterator<Item = &'a str>,
        ) -> std::result::Result<std::collections::BTreeSet<&'a str>, String> {
            let mut seen = std::collections::BTreeSet::new();
            for id in ids {
                if !seen.insert(id) {
                    return Err(format!("Found duplicate {what} id (must be unique): {id}"));
                }
            }
            Ok(seen)
        }
        let check = || -> std::result::Result<(), String> {
            let proteins = first_duplicate("protein", self.proteins.iter().map(|p| p.id.as_str()))?;
            let peptides = first_duplicate("peptide", self.peptides.iter().map(|p| p.id.as_str()))?;
            let compounds =
                first_duplicate("compound", self.compounds.iter().map(|c| c.id.as_str()))?;
            first_duplicate("transition", self.transitions.iter().map(|t| t.native_id()))?;
            for peptide in &self.peptides {
                for protein in &peptide.protein_refs {
                    if !proteins.contains(protein.as_str()) {
                        return Err(format!(
                            "Protein {protein} is not present in the provided data structure."
                        ));
                    }
                }
            }
            for transition in &self.transitions {
                if !transition.peptide_ref().is_empty() {
                    if !peptides.contains(transition.peptide_ref()) {
                        return Err(format!(
                            "Peptide {} is not present in the provided data structure.",
                            transition.peptide_ref()
                        ));
                    }
                } else if !transition.compound_ref().is_empty()
                    && !compounds.contains(transition.compound_ref())
                {
                    return Err(format!(
                        "Compound {} is not present in the provided data structure.",
                        transition.compound_ref()
                    ));
                }
            }
            Ok(())
        };
        check().err()
    }
}

fn not_found(reference: &str) -> Error {
    Error::MissingInformation(format!("element not found: {reference}"))
}

fn merged_targets(own: &CVTermList, rhs: &CVTermList) -> Result<CVTermList> {
    let mut targets = own.clone();
    for terms in rhs.terms().values() {
        for term in terms {
            targets.add(term.clone())?;
        }
    }
    Ok(targets)
}
