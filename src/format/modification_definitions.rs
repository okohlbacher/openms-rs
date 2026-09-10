// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Portable named modification definitions from OpenMS ModificationDefinitionIO.
//! Registries are caller-owned; file input never mutates the global database.
use crate::chemistry::ModificationProvenance;
use crate::chemistry::{ModificationsDB, ResidueModification, SequenceModification};
use crate::identification::{PeptideIdentification, ProteinIdentification, SearchParameters};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

pub const METADATA_KEY: &str = "modification_definitions";
pub const MAX_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_RECORDS: usize = 100_000;
pub const MAX_WORK: usize = 50_000_000;
pub type DefinitionsByRun = BTreeMap<String, Vec<Arc<ResidueModification>>>;
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DefinitionDiagnostic {
    pub record: usize,
    pub message: String,
}
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RegistrationReport {
    /// Accepted records, including identical definitions already registered.
    pub registered: usize,
    pub diagnostics: Vec<DefinitionDiagnostic>,
}
fn definition_text(value: &crate::metadata::MetaValue) -> Result<&str> {
    if value.unit().is_some() {
        return Err(Error::Unsupported(
            "modification definitions cannot carry a metadata unit".into(),
        ));
    }
    value.as_str()
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn spend(remaining: &mut usize, amount: usize) -> Result<()> {
    *remaining = remaining
        .checked_sub(amount)
        .ok_or_else(|| invalid("modification definition resource limit exceeded"))?;
    Ok(())
}
pub fn is_definition(modification: &ResidueModification) -> bool {
    modification.provenance() == ModificationProvenance::Defined && !modification.name().is_empty()
}
fn checked_record(modification: &ResidueModification) -> Result<String> {
    let text = modification.to_definition_string()?;
    let decoded = ResidueModification::from_definition_string(&text)?;
    if decoded != *modification {
        return Err(Error::Unsupported(format!(
            "version-one definition cannot preserve all fields of {}",
            modification.full_id()
        )));
    }
    Ok(text)
}
/// Encode a deterministic lossless set; equal allocations are deduplicated.
pub fn encode(definitions: &[Arc<ResidueModification>]) -> Result<String> {
    if definitions.len() > MAX_RECORDS {
        return Err(invalid("definition count exceeds limit"));
    }
    let mut bytes = MAX_BYTES;
    let mut records = BTreeSet::new();
    for definition in definitions {
        spend(
            &mut bytes,
            definition
                .payload_bytes()?
                .saturating_mul(2)
                .saturating_add(128),
        )?;
        let record = checked_record(definition)?;
        records.insert(record);
    }
    Ok(records.into_iter().collect::<Vec<_>>().join(";"))
}
/// Merge definition strings into the source String-valued metadata key.
/// Existing empty metadata remains unchanged when no records are supplied.
pub fn attach(
    parameters: &mut SearchParameters,
    definitions: &[Arc<ResidueModification>],
) -> Result<()> {
    let encoded = encode(definitions)?;
    let existing = parameters
        .metadata
        .get(METADATA_KEY)
        .map(definition_text)
        .transpose()?
        .unwrap_or("");
    if existing.len().saturating_add(encoded.len()) > MAX_BYTES {
        return Err(invalid("definition bytes exceed limit"));
    }
    let mut records: BTreeSet<&str> = ResidueModification::split_definition_records(existing)?
        .into_iter()
        .collect();
    records.extend(ResidueModification::split_definition_records(&encoded)?);
    if records.len() > MAX_RECORDS {
        return Err(invalid("definition count exceeds limit"));
    }
    if !records.is_empty() {
        let text = records.into_iter().collect::<Vec<_>>().join(";");
        parameters.metadata.insert(METADATA_KEY.into(), text.into());
    }
    Ok(())
}
/// Source per-record recovery is returned explicitly as diagnostics. A chemical
/// name collision or resource failure aborts the transaction, leaving DB intact.
pub fn register_from(text: &str, registry: &mut ModificationsDB) -> Result<RegistrationReport> {
    register_impl(text, registry, false, &mut { MAX_WORK }, &mut {
        MAX_BYTES * 16
    })
}
fn register_impl(
    text: &str,
    registry: &mut ModificationsDB,
    strict: bool,
    work: &mut usize,
    bytes: &mut usize,
) -> Result<RegistrationReport> {
    let records = ResidueModification::split_definition_records(text)?;

    let mut additions: BTreeMap<String, ResidueModification> = BTreeMap::new();
    let mut report = RegistrationReport::default();
    for (index, record) in records.into_iter().enumerate() {
        spend(work, record.len().saturating_mul(16).saturating_add(1))?;
        let definition = match ResidueModification::from_definition_string(record) {
            Ok(value) => value,
            Err(error) if !strict => {
                report.diagnostics.push(DefinitionDiagnostic {
                    record: index + 1,
                    message: error.to_string(),
                });
                continue;
            }
            Err(error) => return Err(error),
        };
        spend(bytes, definition.payload_bytes()?.saturating_add(128))?;
        if let Some(existing) = additions.get(definition.full_id()) {
            if existing != &definition {
                return Err(invalid("conflicting modification definition"));
            }
        } else {
            let matches = registry.find(definition.full_id(), None, None);
            spend(
                work,
                matches
                    .len()
                    .saturating_mul(definition.full_id().len().saturating_add(1)),
            )?;
            if matches.iter().any(|existing| **existing != definition) {
                return Err(invalid("definition conflicts with registry chemistry"));
            }
            if matches.is_empty() {
                additions.insert(definition.full_id().into(), definition);
            }
        }
        report.registered += 1;
    }
    // extend_records itself stages indices and payload atomically.
    if !additions.is_empty() {
        registry.charge_extension(work, bytes)?;
        registry.extend_records(additions.into_values().collect())?;
    }
    Ok(report)
}
pub fn register_search_parameters(
    parameters: &SearchParameters,
    registry: &mut ModificationsDB,
) -> Result<RegistrationReport> {
    register_search_parameters_with_budget(parameters, registry, &mut { MAX_WORK }, &mut {
        MAX_BYTES * 16
    })
}
pub fn register_search_parameters_with_budget(
    parameters: &SearchParameters,
    registry: &mut ModificationsDB,
    work: &mut usize,
    bytes: &mut usize,
) -> Result<RegistrationReport> {
    spend(work, 1)?;
    match parameters.metadata.get(METADATA_KEY) {
        None => Ok(RegistrationReport::default()),
        Some(value) => register_impl(definition_text(value)?, registry, true, work, bytes),
    }
}
pub fn collect(
    proteins: &[ProteinIdentification],
    peptides: &[PeptideIdentification],
    registry: &ModificationsDB,
) -> Result<DefinitionsByRun> {
    collect_iter(proteins, peptides, registry)
}
pub fn collect_iter<'a>(
    proteins: &[ProteinIdentification],
    peptides: impl IntoIterator<Item = &'a PeptideIdentification>,
    registry: &ModificationsDB,
) -> Result<DefinitionsByRun> {
    collect_iter_with_budget(proteins, peptides, registry, &mut { MAX_WORK }, &mut {
        MAX_BYTES
    })
}
pub fn collect_iter_with_budget<'a>(
    proteins: &[ProteinIdentification],
    peptides: impl IntoIterator<Item = &'a PeptideIdentification>,
    registry: &ModificationsDB,
    work: &mut usize,
    bytes: &mut usize,
) -> Result<DefinitionsByRun> {
    let mut result = Collector {
        work,
        bytes,
        count: 0,
        records: BTreeMap::new(),
    };
    for protein in proteins {
        spend(result.work, 1)?;
        for name in protein
            .search_parameters
            .fixed_modifications
            .iter()
            .chain(&protein.search_parameters.variable_modifications)
        {
            // Provider order is irrelevant to the complete all-specificity set.
            spend(
                result.work,
                name.len().saturating_add(registry.entries().len()),
            )?;
            spend(
                result.bytes,
                registry
                    .entries()
                    .len()
                    .saturating_mul(size_of::<Arc<ResidueModification>>()),
            )?;
            for record in registry.find_handles(name, None, None) {
                result.add(&protein.identifier, &record)?;
            }
        }
    }
    for peptide in peptides {
        spend(result.work, 1)?;
        for hit in &peptide.hits {
            let sequence = &hit.sequence;
            spend(result.work, sequence.len().saturating_add(3))?;
            for annotation in sequence
                .n_terminal_modification()
                .into_iter()
                .chain(sequence.c_terminal_modification())
                .chain(
                    (0..sequence.len())
                        .filter_map(|i| sequence.residue_modification(i).ok().flatten()),
                )
            {
                if let SequenceModification::Known(record) = annotation {
                    result.add(&peptide.identifier, record)?;
                }
            }
        }
    }
    Ok(result
        .records
        .into_iter()
        .map(|(run, records)| (run, records.into_values().collect()))
        .collect())
}
struct Collector<'a> {
    work: &'a mut usize,
    bytes: &'a mut usize,
    count: usize,
    records: BTreeMap<String, BTreeMap<String, Arc<ResidueModification>>>,
}
impl Collector<'_> {
    fn add(&mut self, run: &str, record: &Arc<ResidueModification>) -> Result<()> {
        spend(
            self.work,
            run.len()
                .saturating_add(record.full_id().len())
                .saturating_add(1)
                .saturating_mul(20),
        )?;
        if !is_definition(record) {
            return Ok(());
        }
        if let Some(existing) = self
            .records
            .get(run)
            .and_then(|records| records.get(record.full_id()))
        {
            spend(self.work, record.payload_bytes()?)?;
            if existing.as_ref() != record.as_ref() {
                return Err(invalid(
                    "conflicting definitions with the same full identifier",
                ));
            }
            return Ok(());
        }
        self.count += 1;
        if self.count > MAX_RECORDS {
            return Err(invalid("definition count exceeds limit"));
        }
        spend(
            self.bytes,
            run.len()
                .saturating_add(record.payload_bytes()?)
                .saturating_add(256),
        )?;
        self.records
            .entry(run.into())
            .or_default()
            .insert(record.full_id().into(), Arc::clone(record));
        Ok(())
    }
}
pub fn encode_by_run(
    proteins: &[ProteinIdentification],
    definitions: &DefinitionsByRun,
) -> Result<BTreeMap<String, String>> {
    encode_by_run_with_budget(proteins, definitions, &mut { MAX_WORK }, &mut { MAX_BYTES })
}
pub fn encode_by_run_with_budget(
    proteins: &[ProteinIdentification],
    definitions: &DefinitionsByRun,
    work: &mut usize,
    bytes: &mut usize,
) -> Result<BTreeMap<String, String>> {
    let mut output = BTreeMap::new();
    for protein in proteins {
        let mut parameters = SearchParameters::default();
        if let Some(value) = protein.search_parameters.metadata.get(METADATA_KEY) {
            let text = definition_text(value)?;
            spend(bytes, text.len())?;
            parameters.metadata.insert(METADATA_KEY.into(), text.into());
        }
        attach_with_budget(
            &mut parameters,
            definitions
                .get(&protein.identifier)
                .map(Vec::as_slice)
                .unwrap_or(&[]),
            work,
            bytes,
        )?;
        if let Some(value) = parameters.metadata.get(METADATA_KEY) {
            let text = definition_text(value)?;
            spend(
                bytes,
                text.len()
                    .saturating_add(protein.identifier.len())
                    .saturating_add(128),
            )?;
            output.insert(protein.identifier.clone(), text.into());
        }
    }
    Ok(output)
}

/// Charge the complete version-one projection and retained union before copying.
pub fn attach_with_budget(
    parameters: &mut SearchParameters,
    definitions: &[Arc<ResidueModification>],
    work: &mut usize,
    bytes: &mut usize,
) -> Result<()> {
    let existing = parameters
        .metadata
        .get(METADATA_KEY)
        .map(definition_text)
        .transpose()?
        .unwrap_or("");
    spend(work, existing.len().saturating_mul(32).saturating_add(1))?;
    spend(bytes, existing.len().saturating_mul(8))?;
    for definition in definitions {
        let size = definition.payload_bytes()?;
        spend(work, size.saturating_mul(32))?;
        spend(bytes, size.saturating_mul(8))?;
    }
    attach(parameters, definitions)
}
