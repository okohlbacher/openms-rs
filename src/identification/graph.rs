// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Hendrik Weisser, OpenMS Rust contributors $
//! Owned sequences, observations and processing provenance from OpenMS IdentificationData.
//!
//! IDs belong to one graph generation. Registration preserves source semantic
//! keys, while references and equal-key comparisons use deterministic value order.
//! Groups, referential cleanup, persistence and the legacy converter remain separate work.

use crate::chemistry::{
    AASequence, AdductInfo, EmpiricalFormula, ModificationsDB, NAFragmentType, NASequence,
    PeptideFragmentType, RibonucleotideDB,
};
use crate::metadata::{MetaInfo, MetaValue};
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

mod records;
pub use records::*;
mod matches;
pub use matches::*;

pub const MAX_GRAPH_RECORDS: usize = 100_000;
pub const MAX_GRAPH_EDGES: usize = 1_000_000;
pub const MAX_GRAPH_BYTES: usize = 256 * 1024 * 1024;
pub const MAX_GRAPH_WORK: usize = 50_000_000;

/// Limits apply cumulatively to retained records and each complete operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GraphLimits {
    pub max_records: usize,
    pub max_edges: usize,
    pub max_bytes: usize,
    pub max_work: usize,
}
impl Default for GraphLimits {
    fn default() -> Self {
        Self {
            max_records: MAX_GRAPH_RECORDS,
            max_edges: MAX_GRAPH_EDGES,
            max_bytes: MAX_GRAPH_BYTES,
            max_work: MAX_GRAPH_WORK,
        }
    }
}
impl GraphLimits {
    fn validate(self) -> Result<()> {
        if self.max_records > MAX_GRAPH_RECORDS
            || self.max_edges > MAX_GRAPH_EDGES
            || self.max_bytes > MAX_GRAPH_BYTES
            || self.max_work > MAX_GRAPH_WORK
        {
            return Err(invalid(
                "identification graph limits exceed native ceilings",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub(crate) struct GraphWork {
    remaining: usize,
    bytes: usize,
}
impl Default for GraphWork {
    fn default() -> Self {
        Self::new(GraphLimits::default())
    }
}
impl GraphWork {
    fn new(limits: GraphLimits) -> Self {
        Self {
            remaining: limits.max_work,
            bytes: limits.max_bytes,
        }
    }
    pub(crate) fn consume(&mut self, amount: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(amount)
            .ok_or_else(|| invalid("identification graph work limit exceeded"))?;
        Ok(())
    }
    fn allocation(&mut self, bytes: usize) -> Result<()> {
        self.bytes = self
            .bytes
            .checked_sub(bytes)
            .ok_or_else(|| invalid("identification graph allocation limit exceeded"))?;
        Ok(())
    }
    pub(crate) fn copy(&mut self, bytes: usize) -> Result<()> {
        self.consume(bytes)?;
        self.allocation(bytes)
    }
}

fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b)
        .ok_or_else(|| invalid("identification graph size overflow"))
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b)
        .ok_or_else(|| invalid("identification graph size overflow"))
}
fn owner() -> Result<u64> {
    static NEXT: AtomicU64 = AtomicU64::new(1);
    NEXT.fetch_update(AtomicOrdering::Relaxed, AtomicOrdering::Relaxed, |x| {
        x.checked_add(1)
    })
    .map_err(|_| invalid("identification graph owner IDs exhausted"))
}
fn reserve<T>(values: &mut Vec<T>, needed: usize, work: &mut GraphWork) -> Result<()> {
    if needed > values.capacity() {
        let capacity = needed
            .max(values.capacity().saturating_mul(2))
            .min(MAX_GRAPH_EDGES.max(needed));
        work.allocation(mul(capacity, std::mem::size_of::<T>())?)?;
        values
            .try_reserve_exact(capacity - values.len())
            .map_err(|_| invalid("identification graph allocation failed"))?;
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, Default)]
struct Size {
    bytes: usize,
    edges: usize,
}

/// Slot order is stable; the separate sorted index avoids cloning chemical keys.
#[derive(Clone, Debug)]
struct Table<T> {
    values: Vec<T>,
    order: Vec<usize>,
    sizes: Vec<Size>,
    max_bytes: usize,
}
impl<T> Default for Table<T> {
    fn default() -> Self {
        Self {
            values: Vec::new(),
            order: Vec::new(),
            sizes: Vec::new(),
            max_bytes: 0,
        }
    }
}
impl<T> Table<T> {
    fn locate(
        &self,
        value: &T,
        bytes: usize,
        cmp: impl Fn(&T, &T) -> Ordering,
        work: &mut GraphWork,
    ) -> Result<std::result::Result<usize, usize>> {
        let comparisons = usize::BITS as usize - self.order.len().leading_zeros() as usize + 1;
        work.consume(mul(add(bytes, self.max_bytes)?, comparisons)?)?;
        Ok(self
            .order
            .binary_search_by(|&slot| cmp(&self.values[slot], value)))
    }
    fn stage_insert(&mut self, work: &mut GraphWork) -> Result<()> {
        let n = add(self.values.len(), 1)?;
        work.consume(n)?;
        reserve(&mut self.values, n, work)?;
        reserve(&mut self.order, n, work)?;
        reserve(&mut self.sizes, n, work)
    }
    fn insert(&mut self, position: usize, value: T, size: Size) -> usize {
        let slot = self.values.len();
        self.values.push(value);
        self.sizes.push(size);
        self.order.insert(position, slot);
        self.max_bytes = self.max_bytes.max(size.bytes);
        slot
    }
    fn replace(&mut self, slot: usize, value: T, size: Size) {
        self.values[slot] = value;
        self.sizes[slot] = size;
        self.max_bytes = self.max_bytes.max(size.bytes);
    }
    fn iter(&self) -> impl Iterator<Item = (usize, &T)> {
        self.order.iter().map(|&slot| (slot, &self.values[slot]))
    }
}

/// A checked graph; use `try_clone_with_translation` to obtain independent IDs.
#[derive(Debug)]
pub struct IdentificationData {
    owner: u64,
    limits: GraphLimits,
    inputs: Table<InputFile>,
    scores: Table<ScoreType>,
    softwares: Table<ProcessingSoftware>,
    searches: Table<DBSearchParam>,
    steps: Table<ProcessingStep>,
    parents: Table<ParentSequence>,
    peptides: Table<IdentifiedPeptide>,
    oligos: Table<IdentifiedOligo>,
    observations: Table<Observation>,
    compounds: Table<IdentifiedCompound>,
    adducts: Table<AdductInfo>,
    observation_matches: Table<ObservationMatch>,
    search_steps: BTreeMap<ProcessingStepId, SearchParamId>,
    current_step: Option<ProcessingStepId>,
    metadata: MetaInfo,
    metadata_bytes: usize,
    retained: Size,
    batch_work: Option<GraphWork>,
}
impl Default for IdentificationData {
    /// Panics only if all available u64 graph-owner identities have been consumed.
    fn default() -> Self {
        Self::new().expect("identification graph owner IDs exhausted")
    }
}
impl IdentificationData {
    pub fn new() -> Result<Self> {
        Self::with_limits(GraphLimits::default())
    }
    pub fn with_limits(limits: GraphLimits) -> Result<Self> {
        limits.validate()?;
        Ok(Self {
            owner: owner()?,
            limits,
            inputs: Table::default(),
            scores: Table::default(),
            softwares: Table::default(),
            searches: Table::default(),
            steps: Table::default(),
            parents: Table::default(),
            peptides: Table::default(),
            oligos: Table::default(),
            observations: Table::default(),
            compounds: Table::default(),
            adducts: Table::default(),
            observation_matches: Table::default(),
            search_steps: BTreeMap::new(),
            current_step: None,
            metadata: MetaInfo::new(),
            metadata_bytes: 0,
            retained: Size::default(),
            batch_work: None,
        })
    }
    pub fn limits(&self) -> GraphLimits {
        self.limits
    }
    pub fn is_empty(&self) -> bool {
        self.record_count() == 0
    }
    pub fn record_count(&self) -> usize {
        self.inputs.values.len()
            + self.scores.values.len()
            + self.softwares.values.len()
            + self.searches.values.len()
            + self.steps.values.len()
            + self.parents.values.len()
            + self.peptides.values.len()
            + self.oligos.values.len()
            + self.observations.values.len()
            + self.compounds.values.len()
            + self.adducts.values.len()
            + self.observation_matches.values.len()
    }
    pub fn metadata(&self) -> &MetaInfo {
        &self.metadata
    }
    /// Invalidate all old IDs, including IDs from earlier clears.
    pub fn clear(&mut self) -> Result<()> {
        *self = Self::with_limits(self.limits)?;
        Ok(())
    }
    pub fn current_processing_step(&self) -> Option<ProcessingStepId> {
        self.current_step
    }
    pub fn set_current_processing_step(&mut self, step: ProcessingStepId) -> Result<()> {
        self.processing_step(step)?;
        self.current_step = Some(step);
        Ok(())
    }
    pub fn clear_current_processing_step(&mut self) {
        self.current_step = None;
    }
    fn operation<T>(
        &mut self,
        f: impl FnOnce(&mut Self, &mut GraphWork) -> Result<T>,
    ) -> Result<T> {
        let batch = self.batch_work.is_some();
        let mut work = self
            .batch_work
            .take()
            .unwrap_or_else(|| GraphWork::new(self.limits));
        let result = f(self, &mut work);
        if batch {
            self.batch_work = Some(work);
        }
        result
    }
    fn checked_size(&self, old: Size, next: Size, insert: bool) -> Result<Size> {
        if insert && self.record_count() >= self.limits.max_records {
            return Err(invalid("identification graph record limit exceeded"));
        }
        let size = Size {
            bytes: add(self.retained.bytes - old.bytes, next.bytes)?,
            edges: add(self.retained.edges - old.edges, next.edges)?,
        };
        if size.bytes > self.limits.max_bytes || size.edges > self.limits.max_edges {
            return Err(invalid(
                "identification graph retained payload or edge limit exceeded",
            ));
        }
        Ok(size)
    }
    fn snapshot(&self, work: &mut GraphWork) -> Result<Self> {
        // One snapshot per explicit batch/merge, never per ordinary registration.
        work.copy(add(self.retained.bytes, mul(self.record_count(), 128)?)?)?;
        Ok(Self {
            owner: self.owner,
            limits: self.limits,
            inputs: self.inputs.clone(),
            scores: self.scores.clone(),
            softwares: self.softwares.clone(),
            searches: self.searches.clone(),
            steps: self.steps.clone(),
            parents: self.parents.clone(),
            peptides: self.peptides.clone(),
            oligos: self.oligos.clone(),
            observations: self.observations.clone(),
            compounds: self.compounds.clone(),
            adducts: self.adducts.clone(),
            observation_matches: self.observation_matches.clone(),
            search_steps: self.search_steps.clone(),
            current_step: self.current_step,
            metadata: self.metadata.clone(),
            retained: self.retained,
            metadata_bytes: self.metadata_bytes,
            batch_work: None,
        })
    }
    pub(crate) fn transaction<T>(&mut self, f: impl FnOnce(&mut Self) -> Result<T>) -> Result<T> {
        if self.batch_work.is_some() {
            return Err(invalid("nested identification graph transaction"));
        }
        let mut work = GraphWork::new(self.limits);
        let mut staged = self.snapshot(&mut work)?;
        staged.batch_work = Some(work);
        let result = f(&mut staged)?;
        staged.batch_work = None;
        *self = staged;
        Ok(result)
    }
}

macro_rules! access {
    ($($table:ident,$kind:ty,$id:ident,$get:ident,$iter:ident,$count:ident;)*)=>{$(
        impl IdentificationData {
            pub fn $get(&self,id:$id)->Result<&$kind> {
                if id.owner!=self.owner { return Err(invalid("foreign or stale identification graph ID")); }
                self.$table.values.get(id.slot).ok_or_else(||invalid("invalid identification graph slot"))
            }
            pub fn $iter(&self)->impl Iterator<Item=($id,&$kind)> {
                self.$table.iter().map(|(slot,value)|($id{owner:self.owner,slot},value))
            }
            pub fn $count(&self)->usize {self.$table.values.len()}
        }
    )*};
}
access! {
    inputs,InputFile,InputFileId,input_file,input_files,input_file_count;
    scores,ScoreType,ScoreTypeId,score_type,score_types,score_type_count;
    softwares,ProcessingSoftware,ProcessingSoftwareId,processing_software,processing_softwares,processing_software_count;
    searches,DBSearchParam,SearchParamId,db_search_param,db_search_params,db_search_param_count;
    steps,ProcessingStep,ProcessingStepId,processing_step,processing_steps,processing_step_count;
    parents,ParentSequence,ParentId,parent,parents,parent_count;
    peptides,IdentifiedPeptide,PeptideId,peptide,peptides,peptide_count;
    oligos,IdentifiedOligo,OligoId,oligo,oligos,oligo_count;
    observations,Observation,ObservationId,observation,observations,observation_count;
    compounds,IdentifiedCompound,CompoundId,compound,compounds,compound_count;
    adducts,AdductInfo,AdductId,adduct,adducts,adduct_count;
    observation_matches,ObservationMatch,ObservationMatchId,observation_match,observation_matches,observation_match_count;
}

fn scored_edges(result: &ScoredProcessingResult) -> Result<usize> {
    result
        .steps_and_scores
        .iter()
        .try_fold(0, |sum, step| add(sum, add(1, step.scores.len())?))
}
fn sequence_edges(
    matches: &BTreeMap<ParentId, BTreeSet<ParentMatch>>,
    result: &ScoredProcessingResult,
) -> Result<usize> {
    matches
        .values()
        .try_fold(scored_edges(result)?, |sum, set| {
            add(sum, add(1, set.len())?)
        })
}
fn record_size(bytes: usize, edges: usize) -> Result<Size> {
    Ok(Size {
        bytes: add(bytes, 3 * std::mem::size_of::<usize>())?,
        edges,
    })
}

impl IdentificationData {
    fn check_result(&self, result: &ScoredProcessingResult, work: &mut GraphWork) -> Result<()> {
        for step in &result.steps_and_scores {
            work.consume(add(step.scores.len(), 1)?)?;
            if let Some(id) = step.processing_step {
                self.processing_step(id)?;
            }
            for &id in step.scores.keys() {
                self.score_type(id)?;
            }
        }
        Ok(())
    }
    fn check_parents(
        &self,
        parents: &BTreeMap<ParentId, BTreeSet<ParentMatch>>,
        expected: MoleculeType,
        work: &mut GraphWork,
    ) -> Result<()> {
        work.consume(parents.len())?;
        for &id in parents.keys() {
            if self.parent(id)?.molecule_type != expected {
                return Err(invalid("unexpected molecule type for parent sequence"));
            }
        }
        Ok(())
    }
    fn apply_current(
        &self,
        result: &mut ScoredProcessingResult,
        work: &mut GraphWork,
    ) -> Result<()> {
        if let Some(id) = self.current_step {
            result.add_processing_step_with_work(
                AppliedProcessingStep {
                    processing_step: Some(id),
                    scores: BTreeMap::new(),
                },
                work,
            )?;
        }
        Ok(())
    }

    pub fn register_input_file(&mut self, mut value: InputFile) -> Result<InputFileId> {
        self.operation(|this, work| {
            if value.name.is_empty() {
                return Err(invalid("input file must have a name"));
            }
            let bytes = value.measure(work)?;
            let found = this
                .inputs
                .locate(&value, bytes, InputFile::key_cmp, work)?;
            let (slot, size, total) = match found {
                Ok(position) => {
                    let slot = this.inputs.order[position];
                    work.copy(this.inputs.sizes[slot].bytes)?;
                    let mut merged = this.inputs.values[slot].clone();
                    merged.merge_with_work(&value, work)?;
                    value = merged;
                    let size = record_size(value.measure(work)?, 0)?;
                    (
                        slot,
                        size,
                        this.checked_size(this.inputs.sizes[slot], size, false)?,
                    )
                }
                Err(position) => {
                    let size = record_size(bytes, 0)?;
                    let total = this.checked_size(Size::default(), size, true)?;
                    this.inputs.stage_insert(work)?;
                    let slot = this.inputs.insert(position, value, size);
                    this.retained = total;
                    return Ok(InputFileId {
                        owner: this.owner,
                        slot,
                    });
                }
            };
            this.inputs.replace(slot, value, size);
            this.retained = total;
            Ok(InputFileId {
                owner: this.owner,
                slot,
            })
        })
    }
    pub fn register_score_type(&mut self, value: ScoreType) -> Result<ScoreTypeId> {
        self.operation(|this, work| {
            let bytes = value.measure(work)?;
            let found = this
                .scores
                .locate(&value, bytes, ScoreType::key_cmp, work)?;
            let slot = match found {
                Ok(position) => {
                    let slot = this.scores.order[position];
                    if this.scores.values[slot].higher_better != value.higher_better {
                        return Err(invalid(
                            "score type already exists with opposite orientation",
                        ));
                    }
                    slot
                }
                Err(position) => {
                    let size = record_size(bytes, 0)?;
                    let total = this.checked_size(Size::default(), size, true)?;
                    this.scores.stage_insert(work)?;
                    let slot = this.scores.insert(position, value, size);
                    this.retained = total;
                    slot
                }
            };
            Ok(ScoreTypeId {
                owner: this.owner,
                slot,
            })
        })
    }
    pub fn register_processing_software(
        &mut self,
        value: ProcessingSoftware,
    ) -> Result<ProcessingSoftwareId> {
        self.operation(|this, work| {
            let bytes = value.measure(work)?;
            work.consume(value.assigned_scores.len())?;
            for &id in &value.assigned_scores {
                this.score_type(id)?;
            }
            let found = this
                .softwares
                .locate(&value, bytes, ProcessingSoftware::key_cmp, work)?;
            let slot = match found {
                Ok(position) => this.softwares.order[position],
                Err(position) => {
                    let size = record_size(bytes, value.assigned_scores.len())?;
                    let total = this.checked_size(Size::default(), size, true)?;
                    this.softwares.stage_insert(work)?;
                    let slot = this.softwares.insert(position, value, size);
                    this.retained = total;
                    slot
                }
            };
            Ok(ProcessingSoftwareId {
                owner: this.owner,
                slot,
            })
        })
    }
    pub fn register_db_search_param(&mut self, value: DBSearchParam) -> Result<SearchParamId> {
        self.operation(|this, work| {
            let bytes = value.measure(work)?;
            let found = this
                .searches
                .locate(&value, bytes, DBSearchParam::key_cmp, work)?;
            let slot = match found {
                Ok(position) => this.searches.order[position],
                Err(position) => {
                    let size = record_size(bytes, 0)?;
                    let total = this.checked_size(Size::default(), size, true)?;
                    this.searches.stage_insert(work)?;
                    let slot = this.searches.insert(position, value, size);
                    this.retained = total;
                    slot
                }
            };
            Ok(SearchParamId {
                owner: this.owner,
                slot,
            })
        })
    }
    pub fn register_processing_step(
        &mut self,
        value: ProcessingStep,
        search: Option<SearchParamId>,
    ) -> Result<ProcessingStepId> {
        self.operation(|this, work| {
            let bytes = value.measure(work)?;
            this.processing_software(value.software)?;
            work.consume(add(value.input_files.len(), 1)?)?;
            for &id in &value.input_files {
                this.input_file(id)?;
            }
            if let Some(id) = search {
                this.db_search_param(id)?;
            }
            let found = this
                .steps
                .locate(&value, bytes, ProcessingStep::key_cmp, work)?;
            let slot = match found {
                Ok(position) => this.steps.order[position],
                Err(_) => this.steps.values.len(),
            };
            let id = ProcessingStepId {
                owner: this.owner,
                slot,
            };
            let link = search.filter(|_| !this.search_steps.contains_key(&id));
            let size = record_size(bytes, add(value.input_files.len(), 1)?)?;
            let mut total = if found.is_err() {
                this.checked_size(Size::default(), size, true)?
            } else {
                this.retained
            };
            if link.is_some() {
                total.bytes = add(total.bytes, 128)?;
                total.edges = add(total.edges, 1)?;
                if total.bytes > this.limits.max_bytes || total.edges > this.limits.max_edges {
                    return Err(invalid("identification graph search-link limit exceeded"));
                }
                work.copy(128)?;
            }
            if let Err(position) = found {
                this.steps.stage_insert(work)?;
                this.steps.insert(position, value, size);
            }
            if let Some(search) = link {
                this.search_steps.insert(id, search);
            }
            this.retained = total;
            Ok(id)
        })
    }
    pub fn search_param_for_step(&self, step: ProcessingStepId) -> Result<Option<SearchParamId>> {
        self.processing_step(step)?;
        Ok(self.search_steps.get(&step).copied())
    }
    pub fn register_parent_sequence(&mut self, mut value: ParentSequence) -> Result<ParentId> {
        self.operation(|this, work| {
            if value.accession.is_empty()
                || !value.coverage.is_finite()
                || !(0.0..=1.0).contains(&value.coverage)
            {
                return Err(invalid(
                    "parent accession must be nonempty and coverage finite in [0,1]",
                ));
            }
            let bytes = value.measure(work)?;
            this.check_result(&value.result, work)?;
            let found = this
                .parents
                .locate(&value, bytes, ParentSequence::key_cmp, work)?;
            let old = if let Ok(position) = found {
                let slot = this.parents.order[position];
                work.copy(this.parents.sizes[slot].bytes)?;
                let mut merged = this.parents.values[slot].clone();
                merged.merge_with_work(&value, work)?;
                value = merged;
                this.parents.sizes[slot]
            } else {
                Size::default()
            };
            this.apply_current(&mut value.result, work)?;
            let size = record_size(value.measure(work)?, scored_edges(&value.result)?)?;
            let total = this.checked_size(old, size, found.is_err())?;
            let slot = match found {
                Ok(position) => {
                    let slot = this.parents.order[position];
                    this.parents.replace(slot, value, size);
                    slot
                }
                Err(position) => {
                    this.parents.stage_insert(work)?;
                    this.parents.insert(position, value, size)
                }
            };
            this.retained = total;
            Ok(ParentId {
                owner: this.owner,
                slot,
            })
        })
    }
}

macro_rules! register_sequence {
    ($table:ident,$kind:ty,$id:ident,$method:ident,$expected:expr) => {
        impl IdentificationData {
            pub fn $method(&mut self, mut value: $kind) -> Result<$id> {
                self.operation(|this, work| {
                    if value.sequence.is_empty() {
                        return Err(invalid("identified sequence must not be empty"));
                    }
                    let bytes = value.measure(work)?;
                    this.check_result(&value.result, work)?;
                    this.check_parents(&value.parent_matches, $expected, work)?;
                    let found = this.$table.locate(
                        &value,
                        bytes,
                        |a, b| a.sequence.cmp(&b.sequence),
                        work,
                    )?;
                    let old = if let Ok(position) = found {
                        let slot = this.$table.order[position];
                        work.copy(this.$table.sizes[slot].bytes)?;
                        let mut merged = this.$table.values[slot].clone();
                        merged.merge_with_work(&value, work)?;
                        value = merged;
                        this.$table.sizes[slot]
                    } else {
                        Size::default()
                    };
                    this.apply_current(&mut value.result, work)?;
                    let size = record_size(
                        value.measure(work)?,
                        sequence_edges(&value.parent_matches, &value.result)?,
                    )?;
                    let total = this.checked_size(old, size, found.is_err())?;
                    let slot = match found {
                        Ok(position) => {
                            let slot = this.$table.order[position];
                            this.$table.replace(slot, value, size);
                            slot
                        }
                        Err(position) => {
                            this.$table.stage_insert(work)?;
                            this.$table.insert(position, value, size)
                        }
                    };
                    this.retained = total;
                    Ok($id {
                        owner: this.owner,
                        slot,
                    })
                })
            }
        }
    };
}
register_sequence!(
    peptides,
    IdentifiedPeptide,
    PeptideId,
    register_identified_peptide,
    MoleculeType::Protein
);
register_sequence!(
    oligos,
    IdentifiedOligo,
    OligoId,
    register_identified_oligo,
    MoleculeType::RNA
);

impl IdentificationData {
    /// Replace graph-level metadata; graph merge deliberately does not adopt it.
    pub fn set_metadata(&mut self, metadata: MetaInfo) -> Result<()> {
        self.operation(|this, work| {
            let holder = ScoredProcessingResult {
                metadata,
                steps_and_scores: Vec::new(),
            };
            let bytes = holder.measure(work)?;
            let total = this.checked_size(
                Size {
                    bytes: this.metadata_bytes,
                    edges: 0,
                },
                Size { bytes, edges: 0 },
                false,
            )?;
            this.metadata = holder.metadata;
            this.metadata_bytes = bytes;
            this.retained = total;
            Ok(())
        })
    }
    pub fn find_score_type(&self, name: &str) -> Result<Option<ScoreTypeId>> {
        let mut work = GraphWork::new(self.limits);
        for (id, score) in self.score_types() {
            work.consume(add(name.len(), score.cv_term.name.len())?)?;
            if score.cv_term.name == name {
                return Ok(Some(id));
            }
        }
        Ok(None)
    }
    fn all_parents_are_decoys(
        &self,
        matches: &BTreeMap<ParentId, BTreeSet<ParentMatch>>,
    ) -> Result<bool> {
        if matches.is_empty() {
            return Err(invalid("identified molecule has no parent"));
        }
        let mut work = GraphWork::new(self.limits);
        work.consume(matches.len())?;
        for &id in matches.keys() {
            if !self.parent(id)?.is_decoy {
                return Ok(false);
            }
        }
        Ok(true)
    }
    pub fn peptide_parents_are_decoys(&self, id: PeptideId) -> Result<bool> {
        self.all_parents_are_decoys(&self.peptide(id)?.parent_matches)
    }
    pub fn oligo_parents_are_decoys(&self, id: OligoId) -> Result<bool> {
        self.all_parents_are_decoys(&self.oligo(id)?.parent_matches)
    }
}

macro_rules! update_result {
    ($table:ident,$id:ident,$getter:ident,$score:ident,$meta:ident) => {
        impl IdentificationData {
            /// Update a score on the most recently inserted step, or on None.
            /// Unlike re-registration, this does not apply the current step.
            pub fn $score(&mut self, id: $id, score: ScoreTypeId, value: f64) -> Result<()> {
                self.operation(|this, work| {
                    this.$getter(id)?;
                    this.score_type(score)?;
                    work.copy(this.$table.sizes[id.slot].bytes)?;
                    let mut record = this.$table.values[id.slot].clone();
                    let step = record
                        .result
                        .steps_and_scores
                        .last()
                        .and_then(|s| s.processing_step);
                    record.result.add_processing_step_with_work(
                        AppliedProcessingStep {
                            processing_step: step,
                            scores: BTreeMap::from([(score, value)]),
                        },
                        work,
                    )?;
                    let bytes = record.measure(work)?;
                    let old = this.$table.sizes[id.slot];
                    let old_result = scored_edges(&this.$table.values[id.slot].result)?;
                    let size = record_size(
                        bytes,
                        add(old.edges - old_result, scored_edges(&record.result)?)?,
                    )?;
                    let total = this.checked_size(old, size, false)?;
                    this.$table.replace(id.slot, record, size);
                    this.retained = total;
                    Ok(())
                })
            }
            pub fn $meta(&mut self, id: $id, key: String, value: MetaValue) -> Result<()> {
                self.operation(|this, work| {
                    this.$getter(id)?;
                    let incoming = ScoredProcessingResult {
                        metadata: BTreeMap::from([(key, value)]),
                        steps_and_scores: Vec::new(),
                    };
                    incoming.measure(work)?;
                    work.copy(this.$table.sizes[id.slot].bytes)?;
                    let mut record = this.$table.values[id.slot].clone();
                    record.result.merge_with_work(&incoming, work)?;
                    let old = this.$table.sizes[id.slot];
                    let size = record_size(record.measure(work)?, old.edges)?;
                    let total = this.checked_size(old, size, false)?;
                    this.$table.replace(id.slot, record, size);
                    this.retained = total;
                    Ok(())
                })
            }
        }
    };
}
update_result!(
    parents,
    ParentId,
    parent,
    add_parent_score,
    set_parent_meta_value
);
update_result!(
    peptides,
    PeptideId,
    peptide,
    add_peptide_score,
    set_peptide_meta_value
);
update_result!(
    oligos,
    OligoId,
    oligo,
    add_oligo_score,
    set_oligo_meta_value
);

update_result!(
    compounds,
    CompoundId,
    compound,
    add_compound_score,
    set_compound_meta_value
);

update_result!(
    observation_matches,
    ObservationMatchId,
    observation_match,
    add_observation_match_score,
    set_observation_match_meta_value
);

/// Explicit translation; a missing source ID is always an error.
#[derive(Clone, Debug, Default)]
pub struct ReferenceTranslator {
    inputs: BTreeMap<InputFileId, InputFileId>,
    scores: BTreeMap<ScoreTypeId, ScoreTypeId>,
    softwares: BTreeMap<ProcessingSoftwareId, ProcessingSoftwareId>,
    searches: BTreeMap<SearchParamId, SearchParamId>,
    steps: BTreeMap<ProcessingStepId, ProcessingStepId>,
    parents: BTreeMap<ParentId, ParentId>,
    peptides: BTreeMap<PeptideId, PeptideId>,
    oligos: BTreeMap<OligoId, OligoId>,
    observations: BTreeMap<ObservationId, ObservationId>,
    compounds: BTreeMap<CompoundId, CompoundId>,
    adducts: BTreeMap<AdductId, AdductId>,
    observation_matches: BTreeMap<ObservationMatchId, ObservationMatchId>,
}
macro_rules! translate {
    ($($table:ident,$id:ident,$method:ident;)*)=>{$(
        impl ReferenceTranslator {
            pub fn $method(&self,id:$id)->Result<$id> {self.$table.get(&id).copied().ok_or_else(||invalid("no translation for identification graph ID"))}
        }
    )*};
}
translate! {
    inputs,InputFileId,input_file;scores,ScoreTypeId,score_type;
    softwares,ProcessingSoftwareId,processing_software;searches,SearchParamId,db_search_param;
    steps,ProcessingStepId,processing_step;parents,ParentId,parent;
    peptides,PeptideId,peptide;oligos,OligoId,oligo;
    observations,ObservationId,observation;compounds,CompoundId,compound;
    adducts,AdductId,adduct;observation_matches,ObservationMatchId,observation_match;
}
impl ReferenceTranslator {
    fn result(&self, result: &mut ScoredProcessingResult) -> Result<()> {
        for step in &mut result.steps_and_scores {
            step.processing_step = step
                .processing_step
                .map(|id| self.processing_step(id))
                .transpose()?;
            let scores = std::mem::take(&mut step.scores);
            for (id, value) in scores {
                step.scores.insert(self.score_type(id)?, value);
            }
        }
        Ok(())
    }
    fn matches(&self, matches: &mut BTreeMap<ParentId, BTreeSet<ParentMatch>>) -> Result<()> {
        let old = std::mem::take(matches);
        for (id, entries) in old {
            matches.entry(self.parent(id)?).or_default().extend(entries);
        }
        Ok(())
    }
}

impl IdentificationData {
    /// Merge records atomically, retaining this graph's current step and metadata.
    pub fn merge_from(&mut self, other: &Self) -> Result<ReferenceTranslator> {
        self.transaction(|staged| staged.merge_in_place(other))
    }
    fn merge_in_place(&mut self, other: &Self) -> Result<ReferenceTranslator> {
        self.operation(|_, work| {
            // Copies, translated score/parent maps and all translator tree nodes.
            work.copy(add(
                mul(other.retained.bytes, 4)?,
                mul(other.record_count(), 256)?,
            )?)
        })?;
        let mut trans = ReferenceTranslator::default();
        for (id, value) in other.input_files() {
            trans
                .inputs
                .insert(id, self.register_input_file(value.clone())?);
        }
        for (id, value) in other.score_types() {
            trans
                .scores
                .insert(id, self.register_score_type(value.clone())?);
        }
        for (id, value) in other.processing_softwares() {
            let mut copy = value.clone();
            for score in &mut copy.assigned_scores {
                *score = trans.score_type(*score)?;
            }
            trans
                .softwares
                .insert(id, self.register_processing_software(copy)?);
        }
        for (id, value) in other.db_search_params() {
            trans
                .searches
                .insert(id, self.register_db_search_param(value.clone())?);
        }
        for (id, value) in other.processing_steps() {
            let mut copy = value.clone();
            copy.software = trans.processing_software(copy.software)?;
            for input in &mut copy.input_files {
                *input = trans.input_file(*input)?;
            }
            trans
                .steps
                .insert(id, self.register_processing_step(copy, None)?);
        }
        // Ordinary registration keeps its first association; graph merge overwrites.
        for (&step, &search) in &other.search_steps {
            let step = trans.processing_step(step)?;
            let search = trans.db_search_param(search)?;
            self.operation(|this, work| {
                work.consume(32)?;
                if !this.search_steps.contains_key(&step) {
                    let total = this.checked_size(
                        Size::default(),
                        Size {
                            bytes: 128,
                            edges: 1,
                        },
                        false,
                    )?;
                    work.copy(128)?;
                    this.retained = total;
                }
                this.search_steps.insert(step, search);
                Ok(())
            })?;
        }
        for (id, value) in other.observations() {
            let mut copy = value.clone();
            copy.input_file = trans.input_file(copy.input_file)?;
            trans
                .observations
                .insert(id, self.register_observation(copy)?);
        }
        for (id, value) in other.parents() {
            let mut copy = value.clone();
            trans.result(&mut copy.result)?;
            trans
                .parents
                .insert(id, self.register_parent_sequence(copy)?);
        }
        for (id, value) in other.peptides() {
            let mut copy = value.clone();
            trans.result(&mut copy.result)?;
            trans.matches(&mut copy.parent_matches)?;
            trans
                .peptides
                .insert(id, self.register_identified_peptide(copy)?);
        }
        for (id, value) in other.oligos() {
            let mut copy = value.clone();
            trans.result(&mut copy.result)?;
            trans.matches(&mut copy.parent_matches)?;
            trans
                .oligos
                .insert(id, self.register_identified_oligo(copy)?);
        }
        for (id, value) in other.compounds() {
            let mut copy = value.clone();
            trans.result(&mut copy.result)?;
            trans
                .compounds
                .insert(id, self.register_identified_compound(copy)?);
        }
        for (id, value) in other.adducts() {
            trans
                .adducts
                .insert(id, self.register_adduct(value.clone())?);
        }
        for (id, value) in other.observation_matches() {
            let mut copy = value.clone();
            copy.identified_molecule = trans.molecule(copy.identified_molecule)?;
            copy.observation = trans.observation(copy.observation)?;
            copy.adduct = copy.adduct.map(|id| trans.adduct(id)).transpose()?;
            trans.result(&mut copy.result)?;
            let old = std::mem::take(&mut copy.peak_annotations);
            for (step, annotations) in old {
                let step = step.map(|id| trans.processing_step(id)).transpose()?;
                // Source translation assigns, while registration keeps the first vector.
                copy.peak_annotations.insert(step, annotations);
            }
            trans
                .observation_matches
                .insert(id, self.register_observation_match(copy)?);
        }
        Ok(trans)
    }
    /// Deep graph copy with new owner IDs, translated active step and graph metadata.
    pub fn try_clone_with_translation(&self) -> Result<(Self, ReferenceTranslator)> {
        let mut copy = Self::with_limits(self.limits)?;
        copy.batch_work = Some(GraphWork::new(self.limits));
        copy.operation(|_, work| work.copy(self.metadata_bytes))?;
        copy.set_metadata(self.metadata.clone())?;
        let trans = copy.merge_in_place(self)?;
        copy.current_step = self
            .current_step
            .map(|id| trans.processing_step(id))
            .transpose()?;
        copy.batch_work = None;
        Ok((copy, trans))
    }
}

#[derive(Default)]
struct Coverage {
    length: usize,
    intervals: Vec<(usize, usize)>,
}
impl IdentificationData {
    /// Union inclusive parent matches after parsing parent text. Empty parents
    /// retain the source's break behavior under deterministic parent-ID ordering.
    pub fn calculate_coverages(&mut self, check_molecule_length: bool) -> Result<()> {
        self.calculate_coverages_with_registries(
            check_molecule_length,
            ModificationsDB::global(),
            RibonucleotideDB::global(),
        )
    }
    /// Use caller chemistry when validating protein/RNA parent text.
    pub fn calculate_coverages_with_registries(
        &mut self,
        check_molecule_length: bool,
        peptides: &ModificationsDB,
        oligos: &RibonucleotideDB,
    ) -> Result<()> {
        self.operation(|this, work| {
            work.consume(this.parents.values.len())?;
            work.allocation(mul(
                this.parents.values.len(),
                std::mem::size_of::<Coverage>(),
            )?)?;
            let mut coverage = Vec::new();
            coverage
                .try_reserve_exact(this.parents.values.len())
                .map_err(|_| invalid("coverage allocation failed"))?;
            coverage.resize_with(this.parents.values.len(), Coverage::default);
            for (_, molecule) in this.peptides() {
                this.coverage_matches(
                    &molecule.parent_matches,
                    if check_molecule_length {
                        molecule.sequence.len()
                    } else {
                        0
                    },
                    false,
                    &mut coverage,
                    peptides,
                    oligos,
                    work,
                )?;
            }
            for (_, molecule) in this.oligos() {
                this.coverage_matches(
                    &molecule.parent_matches,
                    if check_molecule_length {
                        molecule.sequence.len()
                    } else {
                        0
                    },
                    true,
                    &mut coverage,
                    peptides,
                    oligos,
                    work,
                )?;
            }
            work.allocation(mul(coverage.len(), std::mem::size_of::<f64>())?)?;
            let mut fractions = Vec::new();
            fractions
                .try_reserve_exact(coverage.len())
                .map_err(|_| invalid("coverage allocation failed"))?;
            for mut parent in coverage {
                let n = parent.intervals.len();
                let comparisons = usize::BITS as usize - n.leading_zeros() as usize + 1;
                work.consume(mul(mul(n, comparisons)?, 8)?)?;
                parent.intervals.sort_unstable();
                let mut covered = 0usize;
                let mut current: Option<(usize, usize)> = None;
                for (start, end) in parent.intervals {
                    match current {
                        Some((left, right)) if start <= right => {
                            current = Some((left, right.max(end)))
                        }
                        Some((left, right)) => {
                            covered = add(covered, right - left)?;
                            current = Some((start, end));
                        }
                        None => current = Some((start, end)),
                    }
                }
                if let Some((left, right)) = current {
                    covered = add(covered, right - left)?;
                }
                fractions.push(if parent.length == 0 {
                    0.0
                } else {
                    covered as f64 / parent.length as f64
                });
            }
            for (parent, fraction) in this.parents.values.iter_mut().zip(fractions) {
                parent.coverage = fraction;
            }
            Ok(())
        })
    }
    #[allow(clippy::too_many_arguments)]
    fn coverage_matches(
        &self,
        matches: &ParentMatches,
        molecule_length: usize,
        rna: bool,
        coverage: &mut [Coverage],
        peptides: &ModificationsDB,
        oligos: &RibonucleotideDB,
        work: &mut GraphWork,
    ) -> Result<()> {
        work.consume(matches.len())?;
        for (&id, positions) in matches {
            let parent = self.parent(id)?;
            let data = &mut coverage[id.slot];
            if data.length == 0 {
                data.length = if rna {
                    NASequence::parse_with_budget(
                        &parent.sequence,
                        oligos,
                        &mut work.remaining,
                        &mut work.bytes,
                    )?
                    .len()
                } else {
                    AASequence::parse_with_budget(
                        &parent.sequence,
                        peptides,
                        &mut work.remaining,
                        &mut work.bytes,
                    )?
                    .len()
                };
                if data.length == 0 {
                    break;
                }
            }
            work.consume(positions.len())?;
            for item in positions {
                item.measure(work)?;
                if item.has_valid_positions(molecule_length, data.length) {
                    let start = item.start_pos.expect("validated match start");
                    let end = add(item.end_pos.expect("validated match end"), 1)?;
                    let length = add(data.intervals.len(), 1)?;
                    reserve(&mut data.intervals, length, work)?;
                    data.intervals.push((start, end));
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[allow(clippy::items_after_test_module)]
mod tests {
    use super::*;
    #[test]
    fn a_whole_transaction_shares_registration_work_and_rolls_back() {
        let limits = GraphLimits {
            max_work: 2_000,
            ..Default::default()
        };
        let mut independent = IdentificationData::with_limits(limits).unwrap();
        for i in 0..20 {
            independent
                .register_input_file(InputFile::new(format!("file{i:02}")))
                .unwrap();
        }
        let mut batch = IdentificationData::with_limits(limits).unwrap();
        let error = batch
            .transaction(|staged| {
                for i in 0..20 {
                    staged.register_input_file(InputFile::new(format!("file{i:02}")))?;
                }
                Ok(())
            })
            .unwrap_err();
        assert!(error.to_string().contains("work limit"), "{error}");
        assert!(batch.is_empty());
        batch
            .register_input_file(InputFile::new("after failure"))
            .unwrap();
        assert_eq!(batch.input_file_count(), 1);
    }
}

impl IdentificationData {
    /// Register an opaque observation ID within one input file. Incoming missing
    /// coordinates replace known ones, and observations never acquire a current step.
    pub fn register_observation(&mut self, mut value: Observation) -> Result<ObservationId> {
        self.operation(|this, work| {
            let bytes = value.measure(work)?;
            this.input_file(value.input_file)?;
            work.consume(1)?;
            let found = this
                .observations
                .locate(&value, bytes, Observation::key_cmp, work)?;
            let old = if let Ok(position) = found {
                let slot = this.observations.order[position];
                work.copy(this.observations.sizes[slot].bytes)?;
                let mut merged = this.observations.values[slot].clone();
                merged.merge_with_work(&value, work)?;
                value = merged;
                this.observations.sizes[slot]
            } else {
                Size::default()
            };
            let size = record_size(value.measure(work)?, 1)?;
            let total = this.checked_size(old, size, found.is_err())?;
            let slot = match found {
                Ok(position) => {
                    let slot = this.observations.order[position];
                    this.observations.replace(slot, value, size);
                    slot
                }
                Err(position) => {
                    this.observations.stage_insert(work)?;
                    this.observations.insert(position, value, size)
                }
            };
            this.retained = total;
            Ok(ObservationId {
                owner: this.owner,
                slot,
            })
        })
    }
    /// Identifier determines identity; the first formula and descriptive fields survive.
    pub fn register_identified_compound(
        &mut self,
        mut value: IdentifiedCompound,
    ) -> Result<CompoundId> {
        self.operation(|this, work| {
            let bytes = value.measure(work)?;
            this.check_result(&value.result, work)?;
            let found = this
                .compounds
                .locate(&value, bytes, IdentifiedCompound::key_cmp, work)?;
            let old = if let Ok(position) = found {
                let slot = this.compounds.order[position];
                work.copy(this.compounds.sizes[slot].bytes)?;
                let mut merged = this.compounds.values[slot].clone();
                merged.merge_with_work(&value, work)?;
                value = merged;
                this.compounds.sizes[slot]
            } else {
                Size::default()
            };
            this.apply_current(&mut value.result, work)?;
            let size = record_size(value.measure(work)?, scored_edges(&value.result)?)?;
            let total = this.checked_size(old, size, found.is_err())?;
            let slot = match found {
                Ok(position) => {
                    let slot = this.compounds.order[position];
                    this.compounds.replace(slot, value, size);
                    slot
                }
                Err(position) => {
                    this.compounds.stage_insert(work)?;
                    this.compounds.insert(position, value, size)
                }
            };
            this.retained = total;
            Ok(CompoundId {
                owner: this.owner,
                slot,
            })
        })
    }
    /// Replace only coordinates, retaining identity and metadata.
    pub fn set_observation_coordinates(
        &mut self,
        id: ObservationId,
        rt: Option<f64>,
        mz: Option<f64>,
    ) -> Result<()> {
        self.operation(|this, work| {
            this.observation(id)?;
            work.consume(2)?;
            if [rt, mz]
                .into_iter()
                .flatten()
                .any(|value| !value.is_finite())
            {
                return Err(invalid(
                    "observation coordinates must be finite when present",
                ));
            }
            let observation = &mut this.observations.values[id.slot];
            observation.rt = rt;
            observation.mz = mz;
            Ok(())
        })
    }
    pub fn set_observation_meta_value(
        &mut self,
        id: ObservationId,
        key: String,
        value: MetaValue,
    ) -> Result<()> {
        self.operation(|this, work| {
            this.observation(id)?;
            let incoming = ScoredProcessingResult {
                metadata: BTreeMap::from([(key, value)]),
                steps_and_scores: Vec::new(),
            };
            let bytes = incoming.measure(work)?;
            work.copy(add(this.observations.sizes[id.slot].bytes, bytes)?)?;
            let mut record = this.observations.values[id.slot].clone();
            // One bounded key comparison and insertion; existing coordinates stay intact.
            work.consume(mul(add(record.metadata.len(), 1)?, bytes)?)?;
            record.metadata.extend(incoming.metadata);
            let old = this.observations.sizes[id.slot];
            let size = record_size(record.measure(work)?, old.edges)?;
            let total = this.checked_size(old, size, false)?;
            this.observations.replace(id.slot, record, size);
            this.retained = total;
            Ok(())
        })
    }
}

fn adduct_key_cmp(left: &AdductInfo, right: &AdductInfo) -> Ordering {
    left.charge().cmp(&right.charge()).then_with(|| {
        left.empirical_formula()
            .graph_cmp(right.empirical_formula())
    })
}
fn observation_match_edges(value: &ObservationMatch) -> Result<usize> {
    add(
        scored_edges(&value.result)?,
        add(
            2 + usize::from(value.adduct.is_some()),
            value
                .peak_annotations
                .keys()
                .filter(|step| step.is_some())
                .count(),
        )?,
    )
}
impl ReferenceTranslator {
    pub fn molecule(&self, value: IdentifiedMolecule) -> Result<IdentifiedMolecule> {
        Ok(match value {
            IdentifiedMolecule::Peptide(id) => self.peptide(id)?.into(),
            IdentifiedMolecule::Compound(id) => self.compound(id)?.into(),
            IdentifiedMolecule::Oligo(id) => self.oligo(id)?.into(),
        })
    }
}
impl IdentificationData {
    fn check_molecule(&self, value: IdentifiedMolecule) -> Result<()> {
        match value {
            IdentifiedMolecule::Peptide(id) => {
                self.peptide(id)?;
            }
            IdentifiedMolecule::Compound(id) => {
                self.compound(id)?;
            }
            IdentifiedMolecule::Oligo(id) => {
                self.oligo(id)?;
            }
        }
        Ok(())
    }
    /// Charge and formula determine identity; retain the first name and n-mer multiplier.
    pub fn register_adduct(&mut self, value: AdductInfo) -> Result<AdductId> {
        self.operation(|this, work| {
            let bytes = matches::measure_adduct(&value, work)?;
            let found = this.adducts.locate(&value, bytes, adduct_key_cmp, work)?;
            let slot = match found {
                Ok(position) => this.adducts.order[position],
                Err(position) => {
                    let size = record_size(bytes, 0)?;
                    let total = this.checked_size(Size::default(), size, true)?;
                    this.adducts.stage_insert(work)?;
                    let slot = this.adducts.insert(position, value, size);
                    this.retained = total;
                    slot
                }
            };
            Ok(AdductId {
                owner: this.owner,
                slot,
            })
        })
    }
    /// Merge one observation/molecule/adduct key. Charge and annotation conflicts
    /// cannot partially overwrite the existing scored result.
    pub fn register_observation_match(
        &mut self,
        mut value: ObservationMatch,
    ) -> Result<ObservationMatchId> {
        self.operation(|this, work| {
            let bytes = value.measure(work)?;
            this.check_molecule(value.identified_molecule)?;
            this.observation(value.observation)?;
            if let Some(id) = value.adduct {
                this.adduct(id)?;
            }
            this.check_result(&value.result, work)?;
            work.consume(add(value.peak_annotations.len(), 3)?)?;
            // Source misses this validation; all stored native references must be local.
            for &step in value.peak_annotations.keys().flatten() {
                this.processing_step(step)?;
            }
            let found =
                this.observation_matches
                    .locate(&value, bytes, ObservationMatch::key_cmp, work)?;
            let old = if let Ok(position) = found {
                let slot = this.observation_matches.order[position];
                work.copy(this.observation_matches.sizes[slot].bytes)?;
                let mut merged = this.observation_matches.values[slot].clone();
                merged.merge_with_work(&value, work)?;
                value = merged;
                this.observation_matches.sizes[slot]
            } else {
                Size::default()
            };
            this.apply_current(&mut value.result, work)?;
            let size = record_size(value.measure(work)?, observation_match_edges(&value)?)?;
            let total = this.checked_size(old, size, found.is_err())?;
            let slot = match found {
                Ok(position) => {
                    let slot = this.observation_matches.order[position];
                    this.observation_matches.replace(slot, value, size);
                    slot
                }
                Err(position) => {
                    this.observation_matches.stage_insert(work)?;
                    this.observation_matches.insert(position, value, size)
                }
            };
            this.retained = total;
            Ok(ObservationMatchId {
                owner: this.owner,
                slot,
            })
        })
    }
    /// Borrow the contiguous match-key range without copying records or sorting scores.
    pub fn matches_for_observation(
        &self,
        id: ObservationId,
    ) -> Result<impl Iterator<Item = (ObservationMatchId, &ObservationMatch)> + '_> {
        self.observation(id)?;
        let mut work = GraphWork::new(self.limits);
        let count = self.observation_matches.order.len();
        work.consume(2 * (usize::BITS as usize - count.leading_zeros() as usize + 1))?;
        let before = self
            .observation_matches
            .order
            .partition_point(|&slot| self.observation_matches.values[slot].observation < id);
        let after = self
            .observation_matches
            .order
            .partition_point(|&slot| self.observation_matches.values[slot].observation <= id);
        Ok(self.observation_matches.order[before..after]
            .iter()
            .map(|&slot| {
                (
                    ObservationMatchId {
                        owner: self.owner,
                        slot,
                    },
                    &self.observation_matches.values[slot],
                )
            }))
    }
    /// At most one winner per observation. Equal scores keep key order; a sole
    /// unscored match is returned only when require_score is false.
    pub fn best_matches(
        &self,
        score: ScoreTypeId,
        require_score: bool,
    ) -> Result<Vec<ObservationMatchId>> {
        let score_type = self.score_type(score)?;
        let mut work = GraphWork::new(self.limits);
        work.consume(self.observation_match_count())?;
        let mut result = Vec::new();
        let mut group = None;
        let mut count = 0usize;
        let mut last = None;
        let mut best: Option<(ObservationMatchId, f64)> = None;
        for (id, value) in self.observation_matches() {
            if group.is_some() && group != Some(value.observation) {
                if let Some(winner) = best
                    .map(|(id, _)| id)
                    .or_else(|| (!require_score && count == 1).then_some(last).flatten())
                {
                    let wanted = add(result.len(), 1)?;
                    reserve(&mut result, wanted, &mut work)?;
                    result.push(winner);
                }
                best = None;
                count = 0;
            }
            group = Some(value.observation);
            count += 1;
            last = Some(id);
            work.consume(value.result.steps_and_scores.len())?;
            for step in value.result.steps_and_scores.iter().rev() {
                work.consume(
                    usize::BITS as usize - step.scores.len().leading_zeros() as usize + 1,
                )?;
                if let Some(&value) = step.scores.get(&score) {
                    if best.is_none_or(|(_, old)| score_type.is_better_score(value, old)) {
                        best = Some((id, value));
                    }
                    break;
                }
            }
        }
        if let Some(winner) = best
            .map(|(id, _)| id)
            .or_else(|| (!require_score && count == 1).then_some(last).flatten())
        {
            let wanted = add(result.len(), 1)?;
            reserve(&mut result, wanted, &mut work)?;
            result.push(winner);
        }
        Ok(result)
    }
    pub fn set_molecule_meta_value(
        &mut self,
        molecule: impl Into<IdentifiedMolecule>,
        key: String,
        value: MetaValue,
    ) -> Result<()> {
        match molecule.into() {
            IdentifiedMolecule::Peptide(id) => self.set_peptide_meta_value(id, key, value),
            IdentifiedMolecule::Compound(id) => self.set_compound_meta_value(id, key, value),
            IdentifiedMolecule::Oligo(id) => self.set_oligo_meta_value(id, key, value),
        }
    }
    pub fn remove_observation_match_meta_value(
        &mut self,
        id: ObservationMatchId,
        key: &str,
    ) -> Result<()> {
        self.operation(|this, work| {
            this.observation_match(id)?;
            let old = this.observation_matches.sizes[id.slot];
            work.copy(old.bytes)?;
            work.consume(mul(
                key.len(),
                this.observation_matches.values[id.slot]
                    .result
                    .metadata
                    .len()
                    .saturating_add(1),
            )?)?;
            let mut value = this.observation_matches.values[id.slot].clone();
            value.result.metadata.remove(key);
            let size = record_size(value.measure(work)?, old.edges)?;
            let total = this.checked_size(old, size, false)?;
            this.observation_matches.replace(id.slot, value, size);
            this.retained = total;
            Ok(())
        })
    }
}

// Formatting is incremental: even huge shared modification names cannot create an
// uncharged temporary string before the configured output limit is checked.
struct MoleculeText {
    value: String,
    work: GraphWork,
}
impl std::fmt::Write for MoleculeText {
    fn write_str(&mut self, text: &str) -> std::fmt::Result {
        self.work.copy(text.len()).map_err(|_| std::fmt::Error)?;
        self.value
            .try_reserve_exact(text.len())
            .map_err(|_| std::fmt::Error)?;
        self.value.push_str(text);
        Ok(())
    }
}
impl IdentificationData {
    /// Sequence text for peptide/RNA; identifier (not formula/name) for a compound.
    pub fn molecule_string(&self, molecule: impl Into<IdentifiedMolecule>) -> Result<String> {
        use std::fmt::Write;
        let mut output = MoleculeText {
            value: String::new(),
            work: GraphWork::new(self.limits),
        };
        let result = match molecule.into() {
            IdentifiedMolecule::Peptide(id) => {
                let sequence = &self.peptide(id)?.sequence;
                output.work.consume(sequence.len().saturating_add(2))?;
                // AASequence Display currently constructs an intermediate owned string
                // and per-annotation text. Bound those before invoking its formatter.
                let mut text_bytes = add(sequence.len(), 4)?;
                for modification in (0..sequence.len())
                    .map(|index| sequence.residue_modification(index))
                    .chain([
                        Ok(sequence.n_terminal_modification()),
                        Ok(sequence.c_terminal_modification()),
                    ])
                {
                    if let Some(modification) = modification? {
                        text_bytes = add(
                            text_bytes,
                            add(
                                modification.name().len().max(modification.full_id().len()),
                                2,
                            )?,
                        )?;
                    }
                }
                output.work.copy(mul(text_bytes, 8)?)?;
                write!(&mut output, "{sequence}")
            }
            IdentifiedMolecule::Compound(id) => output.write_str(&self.compound(id)?.identifier),
            IdentifiedMolecule::Oligo(id) => {
                let sequence = &self.oligo(id)?.sequence;
                output.work.consume(sequence.len())?;
                write!(&mut output, "{sequence}")
            }
        };
        result
            .map_err(|_| invalid("identification molecule string exceeds work/allocation limit"))?;
        Ok(output.value)
    }
    /// Dispatch a typed fragment formula. Compounds retain their stored formula and
    /// ignore kind/charge, while peptide and RNA use their distinct source enums.
    pub fn molecule_formula(
        &self,
        molecule: impl Into<IdentifiedMolecule>,
        kind: MoleculeFormulaKind,
        charge: i32,
    ) -> Result<EmpiricalFormula> {
        let mut work = GraphWork::new(self.limits);
        match molecule.into() {
            IdentifiedMolecule::Peptide(id) => {
                let sequence = &self.peptide(id)?.sequence;
                let fragment = match kind {
                    MoleculeFormulaKind::Full => PeptideFragmentType::Full,
                    MoleculeFormulaKind::Peptide(fragment) => fragment,
                    MoleculeFormulaKind::Rna(_) => {
                        return Err(invalid("RNA fragment requested for a peptide"));
                    }
                };
                sequence.formula_for_with_budget(
                    fragment,
                    charge,
                    &mut work.remaining,
                    &mut work.bytes,
                )
            }
            IdentifiedMolecule::Compound(id) => {
                let formula = &self.compound(id)?.formula;
                // A sparse BTreeMap still allocates a root node; use the same
                // conservative map allowance as the peptide/RNA formula helpers.
                work.copy(add(512, mul(formula.stored_atom_types(), 128)?)?)?;
                Ok(formula.clone())
            }
            IdentifiedMolecule::Oligo(id) => {
                let sequence = &self.oligo(id)?.sequence;
                let fragment = match kind {
                    MoleculeFormulaKind::Full => NAFragmentType::Full,
                    MoleculeFormulaKind::Rna(fragment) => fragment,
                    MoleculeFormulaKind::Peptide(_) => {
                        return Err(invalid("peptide fragment requested for RNA"));
                    }
                };
                sequence.formula_with_budget(fragment, charge, &mut work.remaining, &mut work.bytes)
            }
        }
    }
    pub fn molecule_full_formula(
        &self,
        molecule: impl Into<IdentifiedMolecule>,
        charge: i32,
    ) -> Result<EmpiricalFormula> {
        self.molecule_formula(molecule, MoleculeFormulaKind::Full, charge)
    }
}
