// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded native peptide-to-protein indexing, derived from OpenMS4-core 7c029e8.
//! All positions refer to the protein with `*` removed. See
//! `docs/PEPTIDE_INDEXING_SUPPORT.md` for source conventions and checked changes.

use crate::chemistry::{DigestionSpecificity, ProductValidation, Protease, ProteaseDigestion};
use crate::format::fasta::FASTAEntry;
use crate::identification::{
    EnzymeTermSpecificity, FlankingResidue, PeptideEvidence, PeptideIdentification, ProteinHit,
    ProteinIdentification, TargetDecoyType,
};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum UnmatchedAction {
    #[default]
    Error,
    Warn,
    Remove,
}
impl UnmatchedAction {
    fn name(self) -> &'static str {
        match self {
            Self::Error => "error",
            Self::Warn => "warn",
            Self::Remove => "remove",
        }
    }
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MissingDecoyAction {
    #[default]
    Error,
    Warn,
    Silent,
}
impl MissingDecoyAction {
    fn name(self) -> &'static str {
        match self {
            Self::Error => "error",
            Self::Warn => "warn",
            Self::Silent => "silent",
        }
    }
}
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum DecoyRule {
    #[default]
    Auto,
    Prefix(String),
    Suffix(String),
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolvedDecoyRule {
    pub affix: String,
    pub is_prefix: bool,
    /// False for an explicit rule or the source's DECOY_ fallback.
    pub inferred: bool,
}
impl ResolvedDecoyRule {
    /// Classification is case-sensitive, including after automatic inference.
    pub fn is_decoy(&self, accession: &str) -> bool {
        if self.is_prefix {
            accession.starts_with(&self.affix)
        } else {
            accession.ends_with(&self.affix)
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolvedIndexingRun {
    pub identifier: String,
    pub enzyme: Protease,
    pub specificity: DigestionSpecificity,
    pub allow_random_asp_pro_cleavage: bool,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IndexingReport {
    pub decoy_rule: ResolvedDecoyRule,
    pub runs: Vec<ResolvedIndexingRun>,
    pub peptide_hits: usize,
    pub target_hits: usize,
    pub decoy_hits: usize,
    pub target_and_decoy_hits: usize,
    pub unmatched_hits: usize,
    pub unique_hits: usize,
    pub non_unique_hits: usize,
    pub evidence_count: usize,
    pub protein_hits: usize,
    /// Normalization, search comparisons and enzyme-validation residue work.
    pub work: u64,
    pub warnings: Vec<String>,
}

/// All public limits are checked before committing either identification slice.
/// The direct matcher caches each distinct unmodified peptide per run; it is
/// intended for bounded workloads, rather than the C++ parallel trie throughput.
#[derive(Clone, Debug)]
pub struct PeptideIndexing {
    pub decoy_rule: DecoyRule,
    /// None resolves the enzyme separately from each run's search parameters.
    pub enzyme: Option<Protease>,
    /// None resolves specificity separately from each run's search parameters.
    pub specificity: Option<DigestionSpecificity>,
    pub max_ambiguities: u8,
    pub max_mismatches: u8,
    pub il_equivalent: bool,
    pub allow_nterm_protein_cleavage: bool,
    pub write_protein_sequence: bool,
    pub write_protein_description: bool,
    pub keep_unreferenced_proteins: bool,
    pub unmatched_action: UnmatchedAction,
    pub missing_decoy_action: MissingDecoyAction,
    /// Maximum FASTA entries, runs, peptide IDs, and total input hits (each).
    pub max_records: usize,
    /// Maximum combined raw FASTA, input peptide and input protein-hit sequence
    /// bytes. Peptide/protein record sequence sizes are checked before cloning.
    pub max_residues: usize,
    /// Maximum raw candidates per peptide/protein pair, cached mappings, output
    /// evidences, and output protein hits (each).
    pub max_matches: usize,
    pub max_work: u64,
}
impl Default for PeptideIndexing {
    fn default() -> Self {
        Self {
            decoy_rule: DecoyRule::Auto,
            enzyme: None,
            specificity: None,
            max_ambiguities: 3,
            max_mismatches: 0,
            il_equivalent: false,
            allow_nterm_protein_cleavage: true,
            write_protein_sequence: false,
            write_protein_description: false,
            keep_unreferenced_proteins: false,
            unmatched_action: UnmatchedAction::Error,
            missing_decoy_action: MissingDecoyAction::Error,
            max_records: 1_000_000,
            max_residues: 100_000_000,
            max_matches: 1_000_000,
            max_work: 200_000_000,
        }
    }
}
fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
fn bounded_add(value: &mut usize, amount: usize, limit: usize, name: &str) -> Result<()> {
    *value = value
        .checked_add(amount)
        .filter(|v| *v <= limit)
        .ok_or_else(|| bad(format!("peptide indexing {name} limit exceeded")))?;
    Ok(())
}
fn work_add(work: &mut u64, amount: usize, limit: u64) -> Result<()> {
    *work = work
        .checked_add(u64::try_from(amount).map_err(|_| bad("indexing work overflow"))?)
        .filter(|v| *v <= limit)
        .ok_or_else(|| bad("peptide indexing work limit exceeded"))?;
    Ok(())
}
impl PeptideIndexing {
    pub fn validate(&self) -> Result<()> {
        if self.max_ambiguities > 10 || self.max_mismatches > 10 {
            return Err(bad("ambiguity and mismatch budgets must each be in 0..=10"));
        }
        if self.max_records == 0
            || self.max_residues == 0
            || self.max_matches == 0
            || self.max_work == 0
        {
            return Err(bad("peptide indexing resource limits must be positive"));
        }
        if matches!(&self.decoy_rule, DecoyRule::Prefix(s) | DecoyRule::Suffix(s) if s.is_empty()) {
            return Err(bad("explicit decoy affix must be nonempty"));
        }
        Ok(())
    }
    /// Raw substring matching without enzyme filtering. Accepts ASCII A-Z,
    /// including ambiguous peptide codes without a known mass or composition.
    /// Removes `*` from both inputs, folds case, and returns sorted start positions.
    pub fn find_matches(&self, peptide: &str, protein: &str) -> Result<Vec<usize>> {
        self.validate()?;
        let mut work = 0;
        let mut residues = 0;
        bounded_add(&mut residues, peptide.len(), self.max_residues, "residue")?;
        bounded_add(&mut residues, protein.len(), self.max_residues, "residue")?;
        let peptide = self.normalize(peptide, false, &mut work)?;
        let protein = self.normalize(protein, true, &mut work)?;
        self.positions(&peptide, &protein, &mut work)
    }
    /// Source DecoyHelper inference, including its suffix-regex conventions.
    /// A failed automatic inference returns the source DECOY_ prefix fallback.
    pub fn resolve_decoy_rule(&self, database: &[FASTAEntry]) -> Result<ResolvedDecoyRule> {
        self.validate()?;
        if database.len() > self.max_records {
            return Err(bad("FASTA record limit exceeded"));
        }
        match &self.decoy_rule {
            DecoyRule::Prefix(s) | DecoyRule::Suffix(s) => Ok(ResolvedDecoyRule {
                affix: s.clone(),
                is_prefix: matches!(self.decoy_rule, DecoyRule::Prefix(_)),
                inferred: false,
            }),
            DecoyRule::Auto => Ok(infer_decoy(database)),
        }
    }
    fn normalize(&self, value: &str, protein: bool, work: &mut u64) -> Result<String> {
        work_add(work, value.len(), self.max_work)?;
        let mut normalized = String::with_capacity(value.len());
        for byte in value.bytes() {
            if byte == b'*' {
                continue;
            }
            if !byte.is_ascii_alphabetic() {
                return Err(bad(
                    "indexing sequences must contain ASCII residue letters or *",
                ));
            }
            let mut byte = byte.to_ascii_uppercase();
            if self.il_equivalent && (byte == b'L' || (protein && byte == b'J')) {
                byte = b'I';
            }
            normalized.push(char::from(byte));
        }
        if !protein && normalized.is_empty() {
            return Err(bad("cannot index an empty peptide"));
        }
        Ok(normalized)
    }
    fn positions(&self, peptide: &str, protein: &str, work: &mut u64) -> Result<Vec<usize>> {
        let mut matches = Vec::new();
        if peptide.len() > protein.len() {
            return Ok(matches);
        }
        for (position, candidate) in protein.as_bytes().windows(peptide.len()).enumerate() {
            let mut ambiguities = 0_u8;
            let mut mismatches = 0_u8;
            let mut matched = true;
            for (&actual, &expected) in candidate.iter().zip(peptide.as_bytes()) {
                work_add(work, 1, self.max_work)?;
                if actual == expected {
                    continue;
                }
                let compatible = match actual {
                    b'B' => matches!(expected, b'D' | b'N'),
                    b'J' => matches!(expected, b'I' | b'L'),
                    b'Z' => matches!(expected, b'E' | b'Q'),
                    b'X' => b"ACDEFGHIKLMNOPQRSTUVWY".contains(&expected),
                    _ => false,
                };
                if compatible && ambiguities < self.max_ambiguities {
                    ambiguities += 1;
                } else {
                    mismatches += 1;
                }
                if mismatches > self.max_mismatches {
                    matched = false;
                    break;
                }
            }
            if matched {
                if matches.len() == self.max_matches {
                    return Err(bad("peptide indexing match limit exceeded"));
                }
                matches.push(position);
            }
        }
        Ok(matches)
    }
    fn resolve_run(
        &self,
        run: &ProteinIdentification,
        warnings: &mut Vec<String>,
    ) -> Result<ResolvedIndexingRun> {
        let parameters = &run.search_parameters;
        let mut enzyme = if let Some(enzyme) = self.enzyme {
            enzyme
        } else if parameters.digestion_enzyme.is_empty()
            || parameters.digestion_enzyme == "unknown_enzyme"
        {
            warnings.push(format!(
                "run {:?}: unknown enzyme; using Trypsin",
                run.identifier
            ));
            Protease::Trypsin
        } else {
            Protease::from_name(&parameters.digestion_enzyme)?
        };
        if self.enzyme.is_none()
            && !parameters.digestion_regex.is_empty()
            && parameters.digestion_regex != enzyme.metadata().regex()
        {
            return Err(Error::Unsupported("indexing cannot execute a custom digestion expression; select a registered enzyme explicitly".into()));
        }
        // ProteinIdentification::getOriginalSearchEngineName recovers the first
        // originating engine for Percolator/ConsensusID results.
        let engine = if run.search_engine.contains("Percolator")
            || run.search_engine.contains("ConsensusID")
        {
            parameters
                .metadata
                .keys()
                .find(|key| key.starts_with("SE:") && !key.contains("percolator"))
                .map(|key| key[3..].to_ascii_uppercase())
                .unwrap_or_else(|| "UNKNOWN".into())
        } else {
            run.search_engine.to_ascii_uppercase()
        };
        if enzyme == Protease::Trypsin
            && (matches!(engine.as_str(), "MS-GF+" | "MSGFPLUS")
                || parameters.metadata.contains_key("SE:MS-GF+"))
        {
            enzyme = Protease::TrypsinP;
        }
        if self.il_equivalent
            && matches!(
                enzyme,
                Protease::Chymotrypsin | Protease::ChymotrypsinP | Protease::TrypChymo
            )
        {
            return Err(bad(
                "I/L equivalence is incompatible with chymotryptic indexing enzymes",
            ));
        }
        let specificity = self
            .specificity
            .unwrap_or_else(|| match parameters.enzyme_specificity {
                EnzymeTermSpecificity::Full => DigestionSpecificity::Full,
                EnzymeTermSpecificity::Semi => DigestionSpecificity::Semi,
                EnzymeTermSpecificity::None => DigestionSpecificity::None,
                EnzymeTermSpecificity::Unknown => {
                    warnings.push(format!(
                        "run {:?}: unknown enzyme specificity; using full",
                        run.identifier
                    ));
                    DigestionSpecificity::Full
                }
            });
        Ok(ResolvedIndexingRun {
            identifier: run.identifier.clone(),
            enzyme,
            specificity,
            allow_random_asp_pro_cleavage: engine == "XTANDEM"
                || parameters.metadata.contains_key("SE:XTandem"),
        })
    }
    /// Replaces evidence and target/decoy labels and reconstructs protein hits per
    /// run. All input validation, matching, and policy checks complete before the
    /// records are changed. Existing peptide scores, modifications, ordering and
    /// unrelated metadata survive; source-style newly matched protein hits reset
    /// previous protein scores/metadata. Protein groups are retained unchanged.
    pub fn run(
        &self,
        database: &[FASTAEntry],
        proteins: &mut [ProteinIdentification],
        peptides: &mut [PeptideIdentification],
    ) -> Result<IndexingReport> {
        self.validate()?;
        if database.is_empty() {
            return Err(bad("cannot index against an empty protein database"));
        }
        if [database.len(), proteins.len(), peptides.len()]
            .iter()
            .any(|n| *n > self.max_records)
        {
            return Err(bad("peptide indexing record limit exceeded"));
        }
        let mut work = 0;
        let mut residues = 0;
        let mut accessions = BTreeSet::new();
        let mut sequences = Vec::with_capacity(database.len());
        for entry in database {
            if entry.identifier.is_empty() || !accessions.insert(entry.identifier.as_str()) {
                return Err(bad("FASTA accessions must be nonempty and unique"));
            }
            bounded_add(
                &mut residues,
                entry.sequence.len(),
                self.max_residues,
                "residue",
            )?;
            sequences.push(self.normalize(&entry.sequence, true, &mut work)?);
        }
        let mut run_ids = BTreeMap::new();
        let mut input_hits = 0;
        let mut warnings = Vec::new();
        let mut settings = Vec::with_capacity(proteins.len());
        for (index, protein) in proteins.iter().enumerate() {
            protein.validate()?;
            bounded_add(
                &mut input_hits,
                protein.hits.len(),
                self.max_records,
                "input hit",
            )?;
            if run_ids.insert(protein.identifier.as_str(), index).is_some() {
                return Err(bad("protein identification run identifiers must be unique"));
            }
            for hit in &protein.hits {
                bounded_add(
                    &mut residues,
                    hit.sequence.len(),
                    self.max_residues,
                    "residue",
                )?;
            }
            settings.push(self.resolve_run(protein, &mut warnings)?);
        }
        let mut peptide_runs = Vec::with_capacity(peptides.len());
        for peptide in peptides.iter() {
            peptide.validate()?;
            bounded_add(
                &mut input_hits,
                peptide.hits.len(),
                self.max_records,
                "input hit",
            )?;
            for hit in &peptide.hits {
                bounded_add(
                    &mut residues,
                    hit.sequence.len(),
                    self.max_residues,
                    "residue",
                )?;
                if hit.sequence.is_empty() {
                    return Err(bad("cannot index an empty peptide"));
                }
            }
            peptide_runs.push(*run_ids.get(peptide.identifier.as_str()).ok_or_else(|| {
                bad(format!(
                    "peptide refers to unknown run {:?}",
                    peptide.identifier
                ))
            })?);
        }
        let decoy_rule = self.resolve_decoy_rule(database)?;
        if matches!(self.decoy_rule, DecoyRule::Auto) && !decoy_rule.inferred {
            warnings
                .push("could not infer decoy naming; using case-sensitive DECOY_ prefix".into());
        }
        let mut report = IndexingReport {
            decoy_rule,
            runs: settings,
            peptide_hits: 0,
            target_hits: 0,
            decoy_hits: 0,
            target_and_decoy_hits: 0,
            unmatched_hits: 0,
            unique_hits: 0,
            non_unique_hits: 0,
            evidence_count: 0,
            protein_hits: 0,
            work: 0,
            warnings,
        };
        let mut new_peptides = peptides.to_vec();
        let mut matched_proteins = vec![BTreeSet::new(); proteins.len()];
        let mut cache: BTreeMap<(usize, String), Vec<(usize, usize)>> = BTreeMap::new();
        let mut cache_size = 0;
        for (identification, &run_index) in new_peptides.iter_mut().zip(&peptide_runs) {
            for hit in &mut identification.hits {
                report.peptide_hits += 1;
                let key = (run_index, hit.sequence.as_str().to_owned());
                if !cache.contains_key(&key) {
                    let needle = self.normalize(&key.1, false, &mut work)?;
                    let settings = &report.runs[run_index];
                    let digestion = ProteaseDigestion {
                        enzyme: settings.enzyme,
                        specificity: settings.specificity,
                        ..Default::default()
                    };
                    let validation = ProductValidation {
                        ignore_missed_cleavages: true,
                        allow_nterm_protein_cleavage: self.allow_nterm_protein_cleavage,
                        allow_random_asp_pro_cleavage: settings.allow_random_asp_pro_cleavage,
                    };
                    let mut mappings = Vec::new();
                    for (protein_index, sequence) in sequences.iter().enumerate() {
                        work_add(&mut work, 1, self.max_work)?;
                        for position in self.positions(&needle, sequence, &mut work)? {
                            // The shared digestion validator scans the protein.
                            // Charge this work as well as the raw matcher work.
                            work_add(&mut work, sequence.len(), self.max_work)?;
                            if digestion.is_valid_product_unmodified(
                                sequence,
                                position..position + needle.len(),
                                validation,
                            )? {
                                bounded_add(&mut cache_size, 1, self.max_matches, "cached match")?;
                                mappings.push((protein_index, position));
                            }
                        }
                    }
                    cache.insert(key.clone(), mappings);
                }
                let mappings = &cache[&key];
                bounded_add(
                    &mut report.evidence_count,
                    mappings.len(),
                    self.max_matches,
                    "output evidence",
                )?;
                hit.evidences.clear();
                let mut hit_proteins = BTreeSet::new();
                let (mut target, mut decoy) = (false, false);
                for &(protein_index, position) in mappings {
                    hit_proteins.insert(protein_index);
                    matched_proteins[run_index].insert(protein_index);
                    let sequence = sequences[protein_index].as_bytes();
                    let end = position + hit.sequence.len();
                    hit.evidences.push(PeptideEvidence {
                        protein_accession: database[protein_index].identifier.clone(),
                        start: Some(position),
                        end: Some(end - 1),
                        aa_before: if position == 0 {
                            FlankingResidue::NTerminus
                        } else {
                            FlankingResidue::from_code(char::from(sequence[position - 1]))?
                        },
                        aa_after: if end == sequence.len() {
                            FlankingResidue::CTerminus
                        } else {
                            FlankingResidue::from_code(char::from(sequence[end]))?
                        },
                    });
                    if report
                        .decoy_rule
                        .is_decoy(&database[protein_index].identifier)
                    {
                        decoy = true;
                    } else {
                        target = true;
                    }
                }
                let label = match (target, decoy) {
                    (true, true) => {
                        report.target_and_decoy_hits += 1;
                        TargetDecoyType::TargetAndDecoy
                    }
                    (true, false) => {
                        report.target_hits += 1;
                        TargetDecoyType::Target
                    }
                    (false, true) => {
                        report.decoy_hits += 1;
                        TargetDecoyType::Decoy
                    }
                    (false, false) => {
                        report.unmatched_hits += 1;
                        TargetDecoyType::Unknown
                    }
                };
                hit.set_target_decoy_type(label);
                let references = match hit_proteins.len() {
                    0 => "unmatched",
                    1 => {
                        report.unique_hits += 1;
                        "unique"
                    }
                    _ => {
                        report.non_unique_hits += 1;
                        "non-unique"
                    }
                };
                hit.metadata
                    .insert("protein_references".into(), references.into());
            }
            if self.unmatched_action == UnmatchedAction::Remove {
                identification.hits.retain(|hit| !hit.evidences.is_empty());
            }
        }
        if report.unmatched_hits > 0 {
            if self.unmatched_action == UnmatchedAction::Error {
                return Err(bad(format!(
                    "{} peptide hits could not be indexed",
                    report.unmatched_hits
                )));
            }
            if self.unmatched_action == UnmatchedAction::Warn {
                report.warnings.push(format!(
                    "{} peptide hits could not be indexed",
                    report.unmatched_hits
                ));
            }
        }
        if report.peptide_hits > 0 && report.decoy_hits + report.target_and_decoy_hits == 0 {
            match self.missing_decoy_action {
                MissingDecoyAction::Error => {
                    return Err(bad("no peptide hit maps to a decoy protein"));
                }
                MissingDecoyAction::Warn => report
                    .warnings
                    .push("no peptide hit maps to a decoy protein".into()),
                MissingDecoyAction::Silent => {}
            }
        }
        let mut new_proteins = proteins.to_vec();
        for (run_index, run) in new_proteins.iter_mut().enumerate() {
            let matched = &matched_proteins[run_index];
            let names: BTreeSet<_> = matched
                .iter()
                .map(|&i| database[i].identifier.as_str())
                .collect();
            if self.keep_unreferenced_proteins {
                run.hits
                    .retain(|hit| !names.contains(hit.accession.as_str()));
                for hit in &mut run.hits {
                    hit.set_target_decoy_type(TargetDecoyType::Unknown)?;
                }
            } else {
                run.hits.clear();
            }
            bounded_add(
                &mut report.protein_hits,
                run.hits.len(),
                self.max_matches,
                "output protein hit",
            )?;
            for &protein_index in matched {
                bounded_add(
                    &mut report.protein_hits,
                    1,
                    self.max_matches,
                    "output protein hit",
                )?;
                let entry = &database[protein_index];
                let mut hit = ProteinHit {
                    accession: entry.identifier.clone(),
                    sequence: if self.write_protein_sequence {
                        entry.sequence.clone()
                    } else {
                        String::new()
                    },
                    ..Default::default()
                };
                if self.write_protein_description {
                    hit.set_description(&entry.description);
                }
                hit.set_target_decoy_type(if report.decoy_rule.is_decoy(&entry.identifier) {
                    TargetDecoyType::Decoy
                } else {
                    TargetDecoyType::Target
                })?;
                run.hits.push(hit);
            }
            self.store_settings(run, &report.runs[run_index], &report.decoy_rule);
            run.validate()?;
        }
        for identification in &new_peptides {
            identification.validate()?;
        }
        report.work = work;
        proteins.clone_from_slice(&new_proteins);
        peptides.clone_from_slice(&new_peptides);
        Ok(report)
    }
    fn store_settings(
        &self,
        run: &mut ProteinIdentification,
        settings: &ResolvedIndexingRun,
        decoy: &ResolvedDecoyRule,
    ) {
        let values = [
            ("decoy_string", decoy.affix.clone()),
            (
                "decoy_string_position",
                if decoy.is_prefix { "prefix" } else { "suffix" }.into(),
            ),
            ("enzyme", settings.enzyme.name().into()),
            ("enzyme_specificity", settings.specificity.name().into()),
            ("IL_equivalent", self.il_equivalent.to_string()),
            (
                "allow_nterm_protein_cleavage",
                self.allow_nterm_protein_cleavage.to_string(),
            ),
            ("unmatched_action", self.unmatched_action.name().into()),
            (
                "missing_decoy_action",
                self.missing_decoy_action.name().into(),
            ),
        ];
        for (key, value) in values {
            run.search_parameters
                .metadata
                .insert(format!("PeptideIndexer:{key}"), value.into());
        }
        run.search_parameters.metadata.insert(
            "PeptideIndexer:aaa_max".into(),
            i64::from(self.max_ambiguities).into(),
        );
        run.search_parameters.metadata.insert(
            "PeptideIndexer:mismatches_max".into(),
            i64::from(self.max_mismatches).into(),
        );
    }
}

// Compile the tiny fixed DecoyHelper regex vocabulary directly. In the pinned
// suffix expression '*' repeats the last LETTER, not a trailing underscore.
fn infer_decoy(database: &[FASTAEntry]) -> ResolvedDecoyRule {
    const AFFIXES: [&str; 11] = [
        "decoy",
        "dec",
        "reverse",
        "rev",
        "reversed",
        "__id_decoy",
        "xxx",
        "shuffled",
        "shuffle",
        "pseudo",
        "random",
    ];
    let mut prefixes: BTreeMap<String, (usize, String)> = BTreeMap::new();
    let mut suffixes: BTreeMap<String, (usize, String)> = BTreeMap::new();
    let (mut prefix_count, mut suffix_count) = (0_usize, 0_usize);
    for entry in database {
        let lower = entry.identifier.to_ascii_lowercase();
        if let Some(affix) = AFFIXES.iter().find(|affix| lower.starts_with(**affix)) {
            let len = affix.len()
                + lower[affix.len()..]
                    .bytes()
                    .take_while(|c| *c == b'_')
                    .count();
            let value = prefixes.entry(lower[..len].into()).or_default();
            value.0 += 1;
            value.1 = entry.identifier[..len].into();
            prefix_count += 1;
        }
        // Regex search tries the leftmost possible start before alternation.
        let suffix = lower
            .bytes()
            .enumerate()
            .filter(|(_, c)| *c == b'_')
            .find_map(|(start, _)| {
                let tail = &lower[start + 1..];
                AFFIXES.iter().find_map(|affix| {
                    let matched = if *affix == "random" {
                        tail == *affix
                    } else {
                        let base = &affix[..affix.len() - 1];
                        tail.starts_with(base)
                            && tail[base.len()..]
                                .bytes()
                                .all(|c| c == affix.as_bytes()[affix.len() - 1])
                    };
                    matched.then_some(start)
                })
            });
        if let Some(start) = suffix {
            let value = suffixes.entry(lower[start..].into()).or_default();
            value.0 += 1;
            value.1 = entry.identifier[start..].into();
            suffix_count += 1;
        }
    }
    // Integer ratios avoid floating threshold rounding; u128 covers usize input.
    if prefix_count != suffix_count
        && (prefix_count as u128 + suffix_count as u128) * 5 >= database.len() as u128 * 2
    {
        for (candidates, total, is_prefix) in [
            (&prefixes, prefix_count, true),
            (&suffixes, suffix_count, false),
        ] {
            for (count, spelling) in candidates.values() {
                if *count as u128 * 5 >= total as u128 * 4
                    && *count as u128 * 5 >= database.len() as u128 * 2
                {
                    return ResolvedDecoyRule {
                        affix: spelling.clone(),
                        is_prefix,
                        inferred: true,
                    };
                }
            }
        }
    }
    ResolvedDecoyRule {
        affix: "DECOY_".into(),
        is_prefix: true,
        inferred: false,
    }
}
