// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Hendrik Weisser, OpenMS Rust contributors $

//! Owned identification records and ordered score history from Core SDK 6bfc0e4.
//! Registration keys are explicit and narrower than full record equality.

use super::GraphWork;
use crate::chemistry::{
    AASequence, DigestionEnzymeProtein, DigestionEnzymeRNA, DigestionSpecificity, EmpiricalFormula,
    NASequence, SequenceModification,
};
use crate::metadata::{
    CVTerm, CVTermList, CompletionTime, MetaInfo, MetaValue, MetaValueData, ProcessingAction,
    Software,
};
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::mem::size_of;
use std::sync::Arc;

/// Maximum conservative logical payload of one graph record (including shared chemistry).
pub const MAX_GRAPH_RECORD_BYTES: usize = 64 * 1024 * 1024;

macro_rules! ids {
    ($($name:ident),+ $(,)?) => {$ (
        #[doc = "An immutable typed reference belonging to one graph generation."]
        #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name {
            pub(crate) owner: u64,
            pub(crate) slot: usize,
        }
        impl $name {
            /// Opaque graph-generation identity; cannot be used to construct a reference.
            pub fn owner(&self) -> u64 { self.owner }
            /// Stable slot within its owning graph generation.
            pub fn slot(&self) -> usize { self.slot }
        }
    )+};
}
ids!(
    InputFileId,
    ScoreTypeId,
    ProcessingSoftwareId,
    SearchParamId,
    ProcessingStepId,
    ParentId,
    PeptideId,
    OligoId,
    ObservationId,
    CompoundId,
    AdductId,
    ObservationMatchId
);

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MoleculeType {
    #[default]
    Protein,
    Compound,
    RNA,
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MassType {
    #[default]
    Monoisotopic,
    Average,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct InputFile {
    pub name: String,
    pub experimental_design_id: String,
    pub primary_files: BTreeSet<String>,
}
impl InputFile {
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            ..Self::default()
        }
    }
    pub(crate) fn key_cmp(&self, other: &Self) -> Ordering {
        self.name.cmp(&other.name)
    }
    pub(crate) fn measure(&self, work: &mut GraphWork) -> Result<usize> {
        let mut m = Meter::new::<Self>(work)?;
        m.text(&self.name)?;
        m.text(&self.experimental_design_id)?;
        m.strings(&self.primary_files)?;
        if self.name.is_empty() {
            return Err(invalid("input file name must not be empty"));
        }
        Ok(m.bytes)
    }
    pub fn merge(&mut self, other: &Self) -> Result<()> {
        self.merge_with_work(other, &mut GraphWork::default())
    }
    pub(crate) fn merge_with_work(&mut self, other: &Self, work: &mut GraphWork) -> Result<()> {
        let bytes = checked_sum(self.measure(work)?, other.measure(work)?)?;
        merge_cost(
            work,
            bytes,
            self.primary_files
                .len()
                .saturating_add(other.primary_files.len()),
        )?;
        conflict(
            &self.experimental_design_id,
            &other.experimental_design_id,
            "experimental design IDs",
        )?;
        work.copy(bytes)?;
        let mut next = self.clone();
        if next.experimental_design_id.is_empty() {
            next.experimental_design_id
                .clone_from(&other.experimental_design_id);
        }
        next.primary_files
            .extend(other.primary_files.iter().cloned());
        next.measure(work)?;
        *self = next;
        Ok(())
    }
}

/// One spectrum or feature identified by an opaque ID within an input file.
/// None represents the source's missing (NaN) coordinate; present values must be finite.
#[derive(Clone, Debug, PartialEq)]
pub struct Observation {
    pub data_id: String,
    pub input_file: InputFileId,
    pub rt: Option<f64>,
    pub mz: Option<f64>,
    pub metadata: MetaInfo,
}
impl Observation {
    pub fn new(data_id: impl Into<String>, input_file: InputFileId) -> Self {
        Self {
            data_id: data_id.into(),
            input_file,
            rt: None,
            mz: None,
            metadata: MetaInfo::new(),
        }
    }
    pub(crate) fn key_cmp(&self, other: &Self) -> Ordering {
        (&self.input_file, &self.data_id).cmp(&(&other.input_file, &other.data_id))
    }
    pub(crate) fn measure(&self, work: &mut GraphWork) -> Result<usize> {
        let mut m = Meter::new::<Self>(work)?;
        m.text(&self.data_id)?;
        m.metadata(&self.metadata)?;
        if self.data_id.is_empty() {
            return Err(invalid("observation data ID must not be empty"));
        }
        for value in [self.rt, self.mz].into_iter().flatten() {
            finite(value, "observation coordinate")?;
        }
        Ok(m.bytes)
    }
    /// Keep identity, overwrite incoming metadata keys and both coordinates (including None).
    pub fn merge(&mut self, other: &Self) -> Result<()> {
        self.merge_with_work(other, &mut GraphWork::default())
    }
    pub(crate) fn merge_with_work(&mut self, other: &Self, work: &mut GraphWork) -> Result<()> {
        let bytes = checked_sum(self.measure(work)?, other.measure(work)?)?;
        merge_cost(
            work,
            bytes,
            self.metadata.len().saturating_add(other.metadata.len()),
        )?;
        work.copy(bytes)?;
        let mut next = self.clone();
        next.metadata.extend(
            other
                .metadata
                .iter()
                .map(|(key, value)| (key.clone(), value.clone())),
        );
        next.rt = other.rt;
        next.mz = other.mz;
        next.measure(work)?;
        *self = next;
        Ok(())
    }
}

/// A compound keyed by identifier. Re-registration merges only its scored result.
#[derive(Clone, Debug, PartialEq)]
pub struct IdentifiedCompound {
    pub identifier: String,
    pub formula: EmpiricalFormula,
    pub name: String,
    /// Source spelling retained for the SMILES string.
    pub smile: String,
    pub inchi: String,
    pub result: ScoredProcessingResult,
}
impl IdentifiedCompound {
    pub fn new(identifier: impl Into<String>) -> Self {
        Self {
            identifier: identifier.into(),
            formula: EmpiricalFormula::default(),
            name: String::new(),
            smile: String::new(),
            inchi: String::new(),
            result: ScoredProcessingResult::default(),
        }
    }
    pub(crate) fn key_cmp(&self, other: &Self) -> Ordering {
        self.identifier.cmp(&other.identifier)
    }
    pub(crate) fn measure(&self, work: &mut GraphWork) -> Result<usize> {
        let mut m = Meter::new::<Self>(work)?;
        for text in [&self.identifier, &self.name, &self.smile, &self.inchi] {
            m.text(text)?;
        }
        m.formula(&self.formula)?;
        m.scored(&self.result)?;
        if self.identifier.is_empty() {
            return Err(invalid("compound identifier must not be empty"));
        }
        Ok(m.bytes)
    }
    pub fn merge(&mut self, other: &Self) -> Result<()> {
        self.merge_with_work(other, &mut GraphWork::default())
    }
    pub(crate) fn merge_with_work(&mut self, other: &Self, work: &mut GraphWork) -> Result<()> {
        let bytes = checked_sum(self.measure(work)?, other.measure(work)?)?;
        merge_cost(
            work,
            bytes,
            self.result
                .metadata
                .len()
                .saturating_add(other.result.metadata.len()),
        )?;
        charge_history_merge(&self.result, &other.result, work)?;
        work.copy(bytes)?;
        let mut next = self.clone();
        next.result.merge_unchecked(&other.result);
        next.measure(work)?;
        *self = next;
        Ok(())
    }
}

/// A CV name alone is valid for an identification score. Its registration key is accession/name.
#[derive(Clone, Debug, PartialEq)]
pub struct ScoreType {
    pub cv_term: CVTerm,
    pub higher_better: bool,
    pub metadata: MetaInfo,
}
impl Default for ScoreType {
    fn default() -> Self {
        Self {
            cv_term: CVTerm::default(),
            higher_better: true,
            metadata: MetaInfo::new(),
        }
    }
}
impl ScoreType {
    pub fn new(name: impl Into<String>, higher_better: bool) -> Self {
        Self {
            cv_term: CVTerm::new("", name, ""),
            higher_better,
            metadata: MetaInfo::new(),
        }
    }
    pub fn is_better_score(&self, first: f64, second: f64) -> bool {
        if self.higher_better {
            first > second
        } else {
            first < second
        }
    }
    pub(crate) fn key_cmp(&self, other: &Self) -> Ordering {
        (&self.cv_term.accession, &self.cv_term.name)
            .cmp(&(&other.cv_term.accession, &other.cv_term.name))
    }
    pub(crate) fn measure(&self, work: &mut GraphWork) -> Result<usize> {
        let mut m = Meter::new::<Self>(work)?;
        m.cv(&self.cv_term, true)?;
        m.metadata(&self.metadata)?;
        if self.cv_term.name.is_empty() {
            return Err(invalid("score type name must not be empty"));
        }
        Ok(m.bytes)
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ProcessingSoftware {
    pub software: Software,
    /// Priority order, including source-permitted duplicate entries.
    pub assigned_scores: Vec<ScoreTypeId>,
}
impl ProcessingSoftware {
    pub fn new(name: impl Into<String>, version: impl Into<String>) -> Self {
        Self {
            software: Software {
                name: name.into(),
                version: version.into(),
                ..Software::default()
            },
            assigned_scores: Vec::new(),
        }
    }
    pub(crate) fn key_cmp(&self, other: &Self) -> Ordering {
        (&self.software.name, &self.software.version)
            .cmp(&(&other.software.name, &other.software.version))
    }
    pub(crate) fn measure(&self, work: &mut GraphWork) -> Result<usize> {
        let mut m = Meter::new::<Self>(work)?;
        m.text(&self.software.name)?;
        m.text(&self.software.version)?;
        m.cv_list(&self.software.cv_terms)?;
        m.slots::<ScoreTypeId>(self.assigned_scores.len())?;
        Ok(m.bytes)
    }
}

/// Owned enzyme identity, replacing the C++ search-parameter raw pointer.
// Pinned protein metadata is a small Copy value of static references; boxing
// would add an allocation solely to shrink a search-parameter enum.
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GraphEnzyme {
    Protein(DigestionEnzymeProtein),
    RNA(Arc<DigestionEnzymeRNA>),
}
impl GraphEnzyme {
    fn key_cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            // Protein records are immutable members of a fixed registry with unique names.
            (Self::Protein(a), Self::Protein(b)) => a.name().cmp(b.name()),
            (Self::RNA(a), Self::RNA(b)) => a.cmp(b),
            (Self::Protein(_), Self::RNA(_)) => Ordering::Less,
            (Self::RNA(_), Self::Protein(_)) => Ordering::Greater,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct DBSearchParam {
    pub molecule_type: MoleculeType,
    pub mass_type: MassType,
    pub database: String,
    pub database_version: String,
    pub taxonomy: String,
    pub charges: BTreeSet<i32>,
    pub fixed_mods: BTreeSet<String>,
    pub variable_mods: BTreeSet<String>,
    pub precursor_mass_tolerance: f64,
    pub fragment_mass_tolerance: f64,
    pub precursor_tolerance_ppm: bool,
    pub fragment_tolerance_ppm: bool,
    pub digestion_enzyme: Option<GraphEnzyme>,
    /// None represents the source's unknown specificity, not nonspecific digestion.
    pub enzyme_term_specificity: Option<DigestionSpecificity>,
    pub missed_cleavages: usize,
    pub min_length: usize,
    pub max_length: usize,
    pub metadata: MetaInfo,
}
impl Default for DBSearchParam {
    fn default() -> Self {
        Self {
            molecule_type: MoleculeType::Protein,
            mass_type: MassType::Monoisotopic,
            database: String::new(),
            database_version: String::new(),
            taxonomy: String::new(),
            charges: BTreeSet::new(),
            fixed_mods: BTreeSet::new(),
            variable_mods: BTreeSet::new(),
            precursor_mass_tolerance: 0.0,
            fragment_mass_tolerance: 0.0,
            precursor_tolerance_ppm: false,
            fragment_tolerance_ppm: false,
            digestion_enzyme: None,
            enzyme_term_specificity: None,
            missed_cleavages: 0,
            min_length: 0,
            max_length: 0,
            metadata: MetaInfo::new(),
        }
    }
}
impl DBSearchParam {
    pub(crate) fn key_cmp(&self, other: &Self) -> Ordering {
        (
            self.molecule_type,
            self.mass_type,
            &self.database,
            &self.database_version,
            &self.taxonomy,
            &self.charges,
            &self.fixed_mods,
            &self.variable_mods,
        )
            .cmp(&(
                other.molecule_type,
                other.mass_type,
                &other.database,
                &other.database_version,
                &other.taxonomy,
                &other.charges,
                &other.fixed_mods,
                &other.variable_mods,
            ))
            .then_with(|| {
                finite_cmp(
                    self.precursor_mass_tolerance,
                    other.precursor_mass_tolerance,
                )
            })
            .then_with(|| finite_cmp(self.fragment_mass_tolerance, other.fragment_mass_tolerance))
            .then_with(|| {
                (self.precursor_tolerance_ppm, self.fragment_tolerance_ppm)
                    .cmp(&(other.precursor_tolerance_ppm, other.fragment_tolerance_ppm))
            })
            .then_with(|| match (&self.digestion_enzyme, &other.digestion_enzyme) {
                (Some(a), Some(b)) => a.key_cmp(b),
                (None, Some(_)) => Ordering::Less,
                (Some(_), None) => Ordering::Greater,
                (None, None) => Ordering::Equal,
            })
            .then_with(|| {
                (
                    specificity(self.enzyme_term_specificity),
                    self.missed_cleavages,
                    self.min_length,
                    self.max_length,
                )
                    .cmp(&(
                        specificity(other.enzyme_term_specificity),
                        other.missed_cleavages,
                        other.min_length,
                        other.max_length,
                    ))
            })
    }
    pub(crate) fn measure(&self, work: &mut GraphWork) -> Result<usize> {
        let mut m = Meter::new::<Self>(work)?;
        for text in [&self.database, &self.database_version, &self.taxonomy] {
            m.text(text)?;
        }
        m.slots::<i32>(self.charges.len())?;
        m.strings(&self.fixed_mods)?;
        m.strings(&self.variable_mods)?;
        finite(self.precursor_mass_tolerance, "precursor mass tolerance")?;
        finite(self.fragment_mass_tolerance, "fragment mass tolerance")?;
        if let Some(enzyme) = &self.digestion_enzyme {
            m.enzyme(enzyme)?;
        }
        m.metadata(&self.metadata)?;
        Ok(m.bytes)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ProcessingStep {
    pub software: ProcessingSoftwareId,
    pub input_files: Vec<InputFileId>,
    /// Explicit absent timestamp replaces the source DateTime sentinel/current-time convenience.
    pub date_time: Option<CompletionTime>,
    pub actions: BTreeSet<ProcessingAction>,
    pub metadata: MetaInfo,
}
impl ProcessingStep {
    pub fn new(software: ProcessingSoftwareId) -> Self {
        Self {
            software,
            input_files: Vec::new(),
            date_time: None,
            actions: BTreeSet::new(),
            metadata: MetaInfo::new(),
        }
    }
    pub(crate) fn key_cmp(&self, other: &Self) -> Ordering {
        (
            self.date_time,
            self.software,
            &self.input_files,
            &self.actions,
        )
            .cmp(&(
                other.date_time,
                other.software,
                &other.input_files,
                &other.actions,
            ))
    }
    pub(crate) fn measure(&self, work: &mut GraphWork) -> Result<usize> {
        let mut m = Meter::new::<Self>(work)?;
        m.slots::<InputFileId>(self.input_files.len())?;
        m.slots::<ProcessingAction>(self.actions.len())?;
        m.metadata(&self.metadata)?;
        Ok(m.bytes)
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AppliedProcessingStep {
    pub processing_step: Option<ProcessingStepId>,
    pub scores: BTreeMap<ScoreTypeId, f64>,
}
impl AppliedProcessingStep {
    pub fn new(processing_step: Option<ProcessingStepId>) -> Self {
        Self {
            processing_step,
            scores: BTreeMap::new(),
        }
    }
    pub(crate) fn measure(&self, work: &mut GraphWork) -> Result<usize> {
        let mut m = Meter::new::<Self>(work)?;
        m.scores(&self.scores)?;
        Ok(m.bytes)
    }
    /// Software priority first (duplicates retained), followed by other scores in typed-ID order.
    pub fn scores_in_order(
        &self,
        assigned_scores: &[ScoreTypeId],
        primary_only: bool,
    ) -> Result<Vec<(ScoreTypeId, f64)>> {
        let mut work = GraphWork::default();
        self.measure(&mut work)?;
        let assigned_scores = if self.processing_step.is_some() {
            assigned_scores
        } else {
            &[]
        };
        work.consume(assigned_scores.len().checked_mul(64).ok_or_else(overflow)?)?;
        let maximum = assigned_scores
            .len()
            .checked_add(self.scores.len())
            .ok_or_else(overflow)?;
        work.consume(
            maximum
                .checked_mul(size_of::<(ScoreTypeId, f64)>() + 64)
                .ok_or_else(overflow)?,
        )?;
        work.copy(maximum.checked_mul(128).ok_or_else(overflow)?)?;
        let mut result = Vec::new();
        let mut done = BTreeSet::new();
        if self.processing_step.is_some() {
            for score in assigned_scores {
                if let Some(value) = self.scores.get(score) {
                    result.push((*score, *value));
                    if primary_only {
                        return Ok(result);
                    }
                    done.insert(*score);
                }
            }
        }
        for (&score, &value) in &self.scores {
            if !done.contains(&score) {
                result.push((score, value));
                if primary_only {
                    break;
                }
            }
        }
        Ok(result)
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ScoredProcessingResult {
    pub metadata: MetaInfo,
    /// Unique by optional processing-step ID; insertion order determines recency.
    pub steps_and_scores: Vec<AppliedProcessingStep>,
}
impl ScoredProcessingResult {
    pub(crate) fn measure(&self, work: &mut GraphWork) -> Result<usize> {
        let mut m = Meter::new::<Self>(work)?;
        m.scored(self)?;
        Ok(m.bytes)
    }
    pub fn add_processing_step(&mut self, step: AppliedProcessingStep) -> Result<()> {
        self.add_processing_step_with_work(step, &mut GraphWork::default())
    }
    pub(crate) fn add_processing_step_with_work(
        &mut self,
        step: AppliedProcessingStep,
        work: &mut GraphWork,
    ) -> Result<()> {
        let bytes = checked_sum(self.measure(work)?, step.measure(work)?)?;
        merge_cost(
            work,
            bytes,
            self.steps_and_scores
                .len()
                .saturating_add(step.scores.len()),
        )?;
        work.copy(bytes)?;
        let mut next = self.clone();
        next.add_unchecked(&step);
        next.measure(work)?;
        *self = next;
        Ok(())
    }
    pub fn add_score(
        &mut self,
        score_type: ScoreTypeId,
        score: f64,
        processing_step: Option<ProcessingStepId>,
    ) -> Result<()> {
        finite(score, "identification score")?;
        self.add_processing_step(AppliedProcessingStep {
            processing_step,
            scores: BTreeMap::from([(score_type, score)]),
        })
    }
    fn add_unchecked(&mut self, step: &AppliedProcessingStep) {
        if let Some(old) = self
            .steps_and_scores
            .iter_mut()
            .find(|old| old.processing_step == step.processing_step)
        {
            old.scores
                .extend(step.scores.iter().map(|(&id, &value)| (id, value)));
        } else {
            self.steps_and_scores.push(step.clone());
        }
    }
    pub(super) fn merge_unchecked(&mut self, other: &Self) {
        for step in &other.steps_and_scores {
            self.add_unchecked(step);
        }
        self.metadata.extend(
            other
                .metadata
                .iter()
                .map(|(key, value)| (key.clone(), value.clone())),
        );
    }
    pub fn merge(&mut self, other: &Self) -> Result<()> {
        self.merge_with_work(other, &mut GraphWork::default())
    }
    pub(crate) fn merge_with_work(&mut self, other: &Self, work: &mut GraphWork) -> Result<()> {
        let bytes = prepare_scored_merge(self, other, work)?;
        work.copy(bytes)?;
        let mut next = self.clone();
        next.merge_unchecked(other);
        next.measure(work)?;
        *self = next;
        Ok(())
    }
    /// View the unique history entries in optional processing-step ID order.
    /// The score-only None entry comes first; this does not change recency.
    pub fn steps_by_processing_step(&self) -> Result<Vec<&AppliedProcessingStep>> {
        let mut work = GraphWork::default();
        self.measure(&mut work)?;
        let count = self.steps_and_scores.len();
        let comparisons = usize::BITS as usize - count.leading_zeros() as usize + 1;
        work.consume(count.checked_mul(comparisons).ok_or_else(overflow)?)?;
        work.copy(
            count
                .checked_mul(size_of::<&AppliedProcessingStep>())
                .ok_or_else(overflow)?,
        )?;
        let mut result: Vec<_> = self.steps_and_scores.iter().collect();
        result.sort_unstable_by_key(|step| step.processing_step);
        Ok(result)
    }
    pub fn score(&self, score_type: ScoreTypeId) -> Option<f64> {
        self.score_and_step(score_type).map(|(score, _)| score)
    }
    pub fn score_at_step(
        &self,
        score_type: ScoreTypeId,
        step: Option<ProcessingStepId>,
    ) -> Option<f64> {
        self.steps_and_scores
            .iter()
            .find(|item| item.processing_step == step)?
            .scores
            .get(&score_type)
            .copied()
    }
    pub fn score_and_step(
        &self,
        score_type: ScoreTypeId,
    ) -> Option<(f64, Option<ProcessingStepId>)> {
        self.steps_and_scores.iter().rev().find_map(|step| {
            step.scores
                .get(&score_type)
                .map(|&score| (score, step.processing_step))
        })
    }
    pub fn number_of_scores(&self) -> usize {
        self.steps_and_scores
            .iter()
            .map(|step| step.scores.len())
            .sum()
    }
    /// Resolve each referenced step's software priorities without copying graph records.
    pub fn most_recent_score<'a>(
        &self,
        mut assigned_scores: impl FnMut(ProcessingStepId) -> Result<&'a [ScoreTypeId]>,
    ) -> Result<Option<(ScoreTypeId, f64)>> {
        let mut work = GraphWork::default();
        self.measure(&mut work)?;
        for step in self.steps_and_scores.iter().rev() {
            let priority = match step.processing_step {
                Some(id) => assigned_scores(id)?,
                None => &[],
            };
            work.consume(priority.len().checked_mul(64).ok_or_else(overflow)?)?;
            for id in priority {
                if let Some(&score) = step.scores.get(id) {
                    return Ok(Some((*id, score)));
                }
            }
            if let Some((&id, &score)) = step.scores.first_key_value() {
                return Ok(Some((id, score)));
            }
        }
        Ok(None)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ParentSequence {
    pub accession: String,
    pub molecule_type: MoleculeType,
    pub sequence: String,
    pub description: String,
    pub coverage: f64,
    pub is_decoy: bool,
    pub result: ScoredProcessingResult,
}
impl ParentSequence {
    pub fn new(accession: impl Into<String>) -> Self {
        Self {
            accession: accession.into(),
            molecule_type: MoleculeType::Protein,
            sequence: String::new(),
            description: String::new(),
            coverage: 0.0,
            is_decoy: false,
            result: ScoredProcessingResult::default(),
        }
    }
    pub(crate) fn key_cmp(&self, other: &Self) -> Ordering {
        self.accession.cmp(&other.accession)
    }
    pub(crate) fn measure(&self, work: &mut GraphWork) -> Result<usize> {
        let mut m = Meter::new::<Self>(work)?;
        for value in [&self.accession, &self.sequence, &self.description] {
            m.text(value)?;
        }
        m.scored(&self.result)?;
        if self.accession.is_empty() {
            return Err(invalid("parent accession must not be empty"));
        }
        if !self.coverage.is_finite() || !(0.0..=1.0).contains(&self.coverage) {
            return Err(invalid("parent coverage must be finite and in [0,1]"));
        }
        Ok(m.bytes)
    }
    pub fn merge(&mut self, other: &Self) -> Result<()> {
        self.merge_with_work(other, &mut GraphWork::default())
    }
    pub(crate) fn merge_with_work(&mut self, other: &Self, work: &mut GraphWork) -> Result<()> {
        let bytes = checked_sum(self.measure(work)?, other.measure(work)?)?;
        merge_cost(
            work,
            bytes,
            self.result
                .metadata
                .len()
                .saturating_add(other.result.metadata.len()),
        )?;
        charge_history_merge(&self.result, &other.result, work)?;
        conflict(&self.sequence, &other.sequence, "parent sequences")?;
        conflict(&self.description, &other.description, "parent descriptions")?;
        work.copy(bytes)?;
        let mut next = self.clone();
        next.result.merge_unchecked(&other.result);
        if next.sequence.is_empty() {
            next.sequence.clone_from(&other.sequence);
        }
        if next.description.is_empty() {
            next.description.clone_from(&other.description);
        }
        next.is_decoy |= other.is_decoy;
        next.measure(work)?;
        *self = next;
        Ok(())
    }
}

/// Inclusive positions. Equality/order ignore neighbors and metadata, as in the source.
#[derive(Clone, Debug)]
pub struct ParentMatch {
    pub start_pos: Option<usize>,
    pub end_pos: Option<usize>,
    pub left_neighbor: String,
    pub right_neighbor: String,
    pub metadata: MetaInfo,
}
impl Default for ParentMatch {
    fn default() -> Self {
        Self::new(None, None)
    }
}
impl ParentMatch {
    pub const UNKNOWN_NEIGHBOR: char = 'X';
    pub const LEFT_TERMINUS: char = '[';
    pub const RIGHT_TERMINUS: char = ']';
    pub fn new(start_pos: Option<usize>, end_pos: Option<usize>) -> Self {
        Self {
            start_pos,
            end_pos,
            left_neighbor: "X".into(),
            right_neighbor: "X".into(),
            metadata: MetaInfo::new(),
        }
    }
    pub fn has_valid_positions(&self, molecule_length: usize, parent_length: usize) -> bool {
        let (Some(start), Some(end)) = (self.start_pos, self.end_pos) else {
            return false;
        };
        let Some(length) = end
            .checked_sub(start)
            .and_then(|length| length.checked_add(1))
        else {
            return false;
        };
        (molecule_length == 0 || length == molecule_length)
            && (parent_length == 0 || end < parent_length)
    }
    pub(crate) fn measure(&self, work: &mut GraphWork) -> Result<usize> {
        let mut m = Meter::new::<Self>(work)?;
        m.parent_match(self)?;
        Ok(m.bytes)
    }
}
fn position_cmp(a: Option<usize>, b: Option<usize>) -> Ordering {
    match (a, b) {
        (Some(a), Some(b)) => a.cmp(&b),
        (None, Some(_)) => Ordering::Greater,
        (Some(_), None) => Ordering::Less,
        (None, None) => Ordering::Equal,
    }
}
impl PartialEq for ParentMatch {
    fn eq(&self, other: &Self) -> bool {
        self.start_pos == other.start_pos && self.end_pos == other.end_pos
    }
}
impl Eq for ParentMatch {}
impl Ord for ParentMatch {
    fn cmp(&self, other: &Self) -> Ordering {
        position_cmp(self.start_pos, other.start_pos)
            .then_with(|| position_cmp(self.end_pos, other.end_pos))
    }
}
impl PartialOrd for ParentMatch {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub type ParentMatches = BTreeMap<ParentId, BTreeSet<ParentMatch>>;

macro_rules! identified_sequence {
    ($name:ident, $sequence:ty, $measure:ident) => {
        #[derive(Clone, Debug, PartialEq)]
        pub struct $name {
            pub sequence: $sequence,
            pub parent_matches: ParentMatches,
            pub result: ScoredProcessingResult,
        }
        impl $name {
            pub fn new(sequence: $sequence) -> Self {
                Self {
                    sequence,
                    parent_matches: ParentMatches::new(),
                    result: ScoredProcessingResult::default(),
                }
            }
            pub(crate) fn measure(&self, work: &mut GraphWork) -> Result<usize> {
                let mut m = Meter::new::<Self>(work)?;
                m.$measure(&self.sequence)?;
                m.parents(&self.parent_matches)?;
                m.scored(&self.result)?;
                if self.sequence.is_empty() {
                    return Err(invalid("identified sequence must not be empty"));
                }
                Ok(m.bytes)
            }
            pub fn merge(&mut self, other: &Self) -> Result<()> {
                self.merge_with_work(other, &mut GraphWork::default())
            }
            pub(crate) fn merge_with_work(
                &mut self,
                other: &Self,
                work: &mut GraphWork,
            ) -> Result<()> {
                let bytes = checked_sum(self.measure(work)?, other.measure(work)?)?;
                merge_cost(
                    work,
                    bytes,
                    self.parent_matches
                        .len()
                        .saturating_add(other.parent_matches.len())
                        .saturating_add(self.result.metadata.len())
                        .saturating_add(other.result.metadata.len()),
                )?;
                charge_history_merge(&self.result, &other.result, work)?;
                work.copy(bytes)?;
                let mut next = self.clone();
                next.result.merge_unchecked(&other.result);
                for (&parent, matches) in &other.parent_matches {
                    let target = next.parent_matches.entry(parent).or_default();
                    // BTreeSet::insert keeps the first equal-position payload.
                    for item in matches {
                        target.insert(item.clone());
                    }
                }
                next.measure(work)?;
                *self = next;
                Ok(())
            }
            pub fn all_parents_are_decoys(
                &self,
                mut parent_is_decoy: impl FnMut(ParentId) -> Result<bool>,
            ) -> Result<bool> {
                if self.parent_matches.is_empty() {
                    return Err(invalid("no parent found for identified molecule"));
                }
                let mut work = GraphWork::default();
                work.consume(self.parent_matches.len())?;
                for &id in self.parent_matches.keys() {
                    if !parent_is_decoy(id)? {
                        return Ok(false);
                    }
                }
                Ok(true)
            }
        }
    };
}
identified_sequence!(IdentifiedPeptide, AASequence, peptide);
identified_sequence!(IdentifiedOligo, NASequence, oligo);

// Counting conservative logical bytes simultaneously bounds traversed values,
// map comparisons and transient data. Charges precede every variable-size loop.
pub(super) struct Meter<'a> {
    pub(super) work: &'a mut GraphWork,
    pub(super) bytes: usize,
}
impl<'a> Meter<'a> {
    pub(super) fn new<T>(work: &'a mut GraphWork) -> Result<Self> {
        let mut result = Self { work, bytes: 0 };
        result.add(size_of::<T>())?;
        Ok(result)
    }
    pub(super) fn add(&mut self, bytes: usize) -> Result<()> {
        self.bytes = checked_sum(self.bytes, bytes)?;
        if self.bytes > MAX_GRAPH_RECORD_BYTES {
            return Err(invalid("graph record exceeds 64 MiB logical payload limit"));
        }
        self.work.consume(bytes.max(1))
    }
    pub(super) fn slots<T>(&mut self, count: usize) -> Result<()> {
        self.add(
            count
                .checked_mul(size_of::<T>().saturating_add(64))
                .ok_or_else(overflow)?,
        )
    }
    pub(super) fn text(&mut self, value: &str) -> Result<()> {
        self.add(value.len())
    }
    fn strings(&mut self, values: &BTreeSet<String>) -> Result<()> {
        self.slots::<String>(values.len())?;
        for value in values {
            self.text(value)?;
        }
        Ok(())
    }
    pub(super) fn metadata(&mut self, values: &MetaInfo) -> Result<()> {
        self.slots::<(String, MetaValue)>(values.len())?;
        for (key, value) in values {
            self.text(key)?;
            self.value(value)?;
        }
        Ok(())
    }
    fn value(&mut self, value: &MetaValue) -> Result<()> {
        if let Some(unit) = value.unit() {
            for text in [unit.accession(), unit.name(), unit.cv_ref()] {
                self.text(text)?;
            }
        }
        match value.data() {
            MetaValueData::Empty | MetaValueData::Integer(_) => {}
            MetaValueData::String(value) => self.text(value)?,
            MetaValueData::Float(value) => finite(*value, "metadata value")?,
            MetaValueData::StringList(values) => {
                self.slots::<String>(values.len())?;
                for value in values {
                    self.text(value)?;
                }
            }
            MetaValueData::IntegerList(values) => self.slots::<i64>(values.len())?,
            MetaValueData::FloatList(values) => {
                self.slots::<f64>(values.len())?;
                for value in values {
                    finite(*value, "metadata list value")?;
                }
            }
        }
        Ok(())
    }
    fn cv(&mut self, value: &CVTerm, allow_empty_accession: bool) -> Result<()> {
        for text in [&value.accession, &value.name, &value.cv_ref] {
            self.text(text)?;
        }
        self.value(&value.value)?;
        if (!allow_empty_accession && value.accession.is_empty())
            || value.accession.chars().any(char::is_whitespace)
        {
            return Err(invalid("invalid graph CV accession"));
        }
        Ok(())
    }
    fn cv_list(&mut self, values: &CVTermList) -> Result<()> {
        self.metadata(&values.metadata)?;
        self.slots::<(String, Vec<CVTerm>)>(values.terms().len())?;
        for (key, terms) in values.terms() {
            self.text(key)?;
            self.slots::<CVTerm>(terms.len())?;
            for term in terms {
                self.cv(term, false)?;
            }
        }
        Ok(())
    }
    fn scores(&mut self, scores: &BTreeMap<ScoreTypeId, f64>) -> Result<()> {
        self.slots::<(ScoreTypeId, f64)>(scores.len())?;
        for value in scores.values() {
            finite(*value, "identification score")?;
        }
        Ok(())
    }
    pub(super) fn scored(&mut self, result: &ScoredProcessingResult) -> Result<()> {
        self.metadata(&result.metadata)?;
        self.slots::<AppliedProcessingStep>(result.steps_and_scores.len())?;
        self.work.copy(
            result
                .steps_and_scores
                .len()
                .checked_mul(96)
                .ok_or_else(overflow)?,
        )?;
        let mut seen = BTreeSet::new();
        for step in &result.steps_and_scores {
            if !seen.insert(step.processing_step) {
                return Err(invalid("duplicate applied processing-step ID"));
            }
            self.scores(&step.scores)?;
        }
        Ok(())
    }
    fn parent_match(&mut self, value: &ParentMatch) -> Result<()> {
        self.text(&value.left_neighbor)?;
        self.text(&value.right_neighbor)?;
        self.metadata(&value.metadata)
    }
    fn parents(&mut self, values: &ParentMatches) -> Result<()> {
        self.slots::<(ParentId, BTreeSet<ParentMatch>)>(values.len())?;
        for matches in values.values() {
            self.slots::<ParentMatch>(matches.len())?;
            for value in matches {
                self.parent_match(value)?;
            }
        }
        Ok(())
    }
    fn enzyme(&mut self, value: &GraphEnzyme) -> Result<()> {
        match value {
            GraphEnzyme::Protein(value) => {
                self.add(size_of::<DigestionEnzymeProtein>())?;
                for text in [
                    value.name(),
                    value.regex(),
                    value.description(),
                    value.psi_id(),
                    value.xtandem_id(),
                ] {
                    self.text(text)?;
                }
                self.slots::<&str>(value.synonyms().len())?;
                for text in value.synonyms() {
                    self.text(text)?;
                }
            }
            GraphEnzyme::RNA(value) => {
                self.add(size_of::<DigestionEnzymeRNA>())?;
                for text in [
                    value.name(),
                    value.regex(),
                    value.regex_description(),
                    value.cuts_after(),
                    value.cuts_before(),
                    value.five_prime_gain(),
                    value.three_prime_gain(),
                ] {
                    self.text(text)?;
                }
                self.strings(value.synonyms())?;
            }
        }
        Ok(())
    }
    pub(super) fn formula(&mut self, value: &EmpiricalFormula) -> Result<()> {
        let count = value.stored_atom_types();
        if count == 0 {
            return Ok(());
        }
        // A sparse formula still owns a BTreeMap root; per-entry bytes alone
        // underestimate copies of one- and two-element compounds and adducts.
        self.add(
            count
                .checked_mul(128)
                .and_then(|bytes| bytes.checked_add(512))
                .ok_or_else(overflow)?,
        )
    }
    fn peptide(&mut self, sequence: &AASequence) -> Result<()> {
        self.slots::<usize>(sequence.len())?;
        self.add(sequence.generation_payload_bytes()?)?;
        for modification in (0..sequence.len())
            .map(|index| sequence.residue_modification(index))
            .chain([
                Ok(sequence.n_terminal_modification()),
                Ok(sequence.c_terminal_modification()),
            ])
        {
            if let Some(SequenceModification::Known(record)) = modification? {
                self.add(size_of::<crate::chemistry::ResidueModification>())?;
                for text in [
                    record.name(),
                    record.full_name(),
                    record.full_id(),
                    record.classification(),
                ] {
                    self.text(text)?;
                }
                if let Some(text) = record.obo_accession() {
                    self.text(text)?;
                }
                self.strings(record.synonyms())?;
                self.formula(record.diff_formula())?;
                if let Some(formula) = record.absolute_formula() {
                    self.formula(formula)?;
                }
                self.slots::<crate::chemistry::NeutralLoss>(record.neutral_losses().len())?;
                for loss in record.neutral_losses() {
                    self.formula(loss.formula())?;
                }
            }
        }
        Ok(())
    }
    fn oligo(&mut self, sequence: &NASequence) -> Result<()> {
        self.slots::<usize>(sequence.len())?;
        self.add(sequence.generation_payload_bytes()?)
    }
}

fn finite(value: f64, name: &str) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(invalid(&format!("{name} must be finite")))
    }
}
fn finite_cmp(a: f64, b: f64) -> Ordering {
    if a == b {
        Ordering::Equal
    } else {
        a.total_cmp(&b)
    }
}
fn specificity(value: Option<DigestionSpecificity>) -> u8 {
    match value {
        Some(DigestionSpecificity::Full) => 0,
        Some(DigestionSpecificity::Semi) => 1,
        Some(DigestionSpecificity::None) => 2,
        None => 3,
    }
}
fn checked_sum(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or_else(overflow)
}
fn overflow() -> Error {
    invalid("identification record size overflows")
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn conflict(old: &str, new: &str, name: &str) -> Result<()> {
    if !old.is_empty() && !new.is_empty() && old != new {
        Err(invalid(&format!("conflicting {name}")))
    } else {
        Ok(())
    }
}
pub(super) fn merge_cost(work: &mut GraphWork, bytes: usize, entries: usize) -> Result<()> {
    if bytes > MAX_GRAPH_RECORD_BYTES {
        return Err(invalid("combined record merge payload exceeds 64 MiB"));
    }
    let factor = usize::BITS as usize - entries.saturating_add(1).leading_zeros() as usize + 3;
    work.consume(bytes.checked_mul(factor).ok_or_else(overflow)?)
}
pub(super) fn charge_history_merge(
    left: &ScoredProcessingResult,
    right: &ScoredProcessingResult,
    work: &mut GraphWork,
) -> Result<()> {
    work.consume(
        left.steps_and_scores
            .len()
            .checked_add(right.steps_and_scores.len())
            .and_then(|count| count.checked_mul(right.steps_and_scores.len()))
            .ok_or_else(overflow)?,
    )
}
pub(super) fn prepare_scored_merge(
    left: &ScoredProcessingResult,
    right: &ScoredProcessingResult,
    work: &mut GraphWork,
) -> Result<usize> {
    let bytes = checked_sum(left.measure(work)?, right.measure(work)?)?;
    merge_cost(
        work,
        bytes,
        left.metadata
            .len()
            .saturating_add(right.metadata.len())
            .saturating_add(left.number_of_scores())
            .saturating_add(right.number_of_scores()),
    )?;
    charge_history_merge(left, right, work)?;
    Ok(bytes)
}
