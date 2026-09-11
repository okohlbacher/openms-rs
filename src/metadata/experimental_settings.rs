// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::{ContactPerson, DocumentIdentifier, HPLC, Instrument, MetaInfo, Sample, SourceFile};
use crate::{Error, Result, data_structures::DateTime, kernel::data_array::Meter};
use std::{collections::BTreeMap, fmt};

/// Complete owned experiment-wide settings. Public fields replace source
/// accessors. Equality includes all fields except document loaded path/type,
/// following DocumentIdentifier. Ordinary Clone/equality/Drop have Rust costs;
/// checked copies preflight recursive sample depth and aggregate payload first.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ExperimentalSettings {
    pub document: DocumentIdentifier,
    pub sample: Sample,
    pub source_files: Vec<SourceFile>,
    pub contacts: Vec<ContactPerson>,
    pub instrument: Instrument,
    pub instrument_configurations: BTreeMap<String, Instrument>,
    pub hplc: HPLC,
    pub date_time: DateTime,
    pub comment: String,
    pub fraction_identifier: String,
    pub metadata: MetaInfo,
}

/// Bounds on a checked aggregate traversal/copy. Records count sample nodes,
/// source files, contacts, instruments and instrument components. Metadata/list
/// elements additionally consume shared work and byte allowances.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExperimentalSettingsLimits {
    pub max_records: usize,
    /// Effective depth is capped at ExperimentalSettings::MAX_CHECKED_SAMPLE_DEPTH.
    pub max_sample_depth: usize,
    pub max_work: usize,
    pub max_bytes: usize,
}
impl Default for ExperimentalSettingsLimits {
    fn default() -> Self {
        Self {
            max_records: 1_000_000,
            max_sample_depth: 64,
            max_work: 50_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}
impl ExperimentalSettings {
    /// Safe ceiling before recursive derived Clone is allowed; public ordinary
    /// Clone/Drop remain ordinary Rust operations on caller-constructed trees.
    pub const MAX_CHECKED_SAMPLE_DEPTH: usize = 64;
    pub fn new() -> Self {
        Self::default()
    }
    /// Resource validation only: source floating fields, partially initialized
    /// DateTime and stale Gradient percentage rows retain their stored states.
    /// Format/algorithm-specific scientific validation belongs to the consumer.
    pub fn validate(&self) -> Result<()> {
        self.validate_with_limits(ExperimentalSettingsLimits::default())
    }
    pub fn validate_with_limits(&self, limits: ExperimentalSettingsLimits) -> Result<()> {
        let (mut work, mut bytes) = (limits.max_work, limits.max_bytes);
        self.measure(limits, &mut work, &mut bytes)
    }
    pub fn checked_clone(&self) -> Result<Self> {
        self.checked_clone_with_limits(ExperimentalSettingsLimits::default())
    }
    pub fn checked_clone_with_limits(&self, limits: ExperimentalSettingsLimits) -> Result<Self> {
        self.validate_with_limits(limits)?;
        Ok(self.clone())
    }
    /// Share a caller's cumulative copy/traversal counters. The structural record
    /// and sample depth caps remain the documented default per aggregate.
    pub(crate) fn with_budget(&self, work: &mut usize, bytes: &mut usize) -> Result<()> {
        self.measure(ExperimentalSettingsLimits::default(), work, bytes)
    }
    fn measure(
        &self,
        limits: ExperimentalSettingsLimits,
        work: &mut usize,
        bytes: &mut usize,
    ) -> Result<()> {
        let mut m = Meter { work, bytes };
        let mut records = limits.max_records;
        m.slots::<Self>(1)?;
        for text in [
            &self.document.identifier,
            &self.document.loaded_file_path,
            &self.comment,
            &self.fraction_identifier,
        ] {
            m.text(text)?;
        }
        m.meta(&self.metadata)?;
        sample(
            &mut m,
            &self.sample,
            &mut records,
            limits.max_sample_depth.min(Self::MAX_CHECKED_SAMPLE_DEPTH),
        )?;
        take(&mut records, self.source_files.len())?;
        m.slots::<SourceFile>(self.source_files.len())?;
        for source in &self.source_files {
            crate::kernel::acquisition_fields::source(&mut m, source)?;
        }
        take(&mut records, self.contacts.len())?;
        m.slots::<ContactPerson>(self.contacts.len())?;
        for contact in &self.contacts {
            for text in [
                &contact.first_name,
                &contact.last_name,
                &contact.institution,
                &contact.email,
                &contact.contact_info,
                &contact.url,
                &contact.address,
            ] {
                m.text(text)?;
            }
            m.meta(&contact.metadata)?;
        }
        instrument(&mut m, &self.instrument, &mut records)?;
        m.tree::<(String, Instrument)>(self.instrument_configurations.len())?;
        // Charge descriptors before scanning keys or entering component payload.
        if self.instrument_configurations.len() > records {
            return Err(limit());
        }
        for (key, value) in &self.instrument_configurations {
            m.text(key)?;
            instrument(&mut m, value, &mut records)?;
        }
        for text in [&self.hplc.instrument, &self.hplc.column, &self.hplc.comment] {
            m.text(text)?;
        }
        let gradient = &self.hplc.gradient;
        m.slots::<String>(gradient.eluents().len())?;
        for value in gradient.eluents() {
            m.text(value)?;
        }
        m.slots::<i32>(gradient.timepoints().len())?;
        m.slots::<Vec<u32>>(gradient.percentages().len())?;
        for row in gradient.percentages() {
            m.slots::<u32>(row.len())?;
        }
        Ok(())
    }
    /// O(1) loss guard for transports without a settings header. Run metadata is
    /// tested separately. Loaded path/type are local provenance, deliberately
    /// excluded; document identifier is persistent identity and is included.
    pub fn has_nondefault_header(&self) -> bool {
        let s = &self.sample;
        let i = &self.instrument;
        let h = &self.hplc;
        !self.document.identifier.is_empty()
            || !s.name.is_empty()
            || !s.organism.is_empty()
            || !s.number.is_empty()
            || !s.comment.is_empty()
            || s.state != super::SampleState::Unknown
            || s.mass.to_bits() != 0
            || s.volume.to_bits() != 0
            || s.concentration.to_bits() != 0
            || !s.subsamples.is_empty()
            || !s.metadata.is_empty()
            || !self.source_files.is_empty()
            || !self.contacts.is_empty()
            || !i.name.is_empty()
            || !i.vendor.is_empty()
            || !i.model.is_empty()
            || !i.customizations.is_empty()
            || !i.ion_sources.is_empty()
            || !i.mass_analyzers.is_empty()
            || !i.ion_detectors.is_empty()
            || !i.software.name.is_empty()
            || !i.software.version.is_empty()
            || !i.software.cv_terms.is_empty()
            || !i.software.cv_terms.metadata.is_empty()
            || i.ion_optics != super::IonOpticsType::Unknown
            || !i.metadata.is_empty()
            || !self.instrument_configurations.is_empty()
            || !h.instrument.is_empty()
            || !h.column.is_empty()
            || h.temperature != 21
            || h.pressure != 0
            || h.flux != 0
            || !h.comment.is_empty()
            || !h.gradient.eluents().is_empty()
            || !h.gradient.timepoints().is_empty()
            || !h.gradient.percentages().is_empty()
            || self.date_time != DateTime::default()
            || !self.comment.is_empty()
            || !self.fraction_identifier.is_empty()
    }
    /// Checks persistent header and run metadata, excluding loaded provenance.
    pub fn has_transport_metadata(&self) -> bool {
        self.has_nondefault_header() || !self.metadata.is_empty()
    }
}
impl fmt::Display for ExperimentalSettings {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("-- EXPERIMENTALSETTINGS BEGIN --\n-- EXPERIMENTALSETTINGS END --\n")
    }
}
fn take(remaining: &mut usize, count: usize) -> Result<()> {
    *remaining = remaining.checked_sub(count).ok_or_else(limit)?;
    Ok(())
}
fn limit() -> Error {
    Error::InvalidValue("experimental settings resource limit exceeded".into())
}
fn sample(m: &mut Meter<'_>, root: &Sample, records: &mut usize, max_depth: usize) -> Result<()> {
    if max_depth == 0 {
        return Err(limit());
    }
    // Iterative depth-first traversal keeps stack size proportional to depth,
    // not the width of an untrusted nested sample vector.
    m.slots::<(&Sample, usize)>(1)?;
    let mut stack = Vec::new();
    stack.try_reserve_exact(1).map_err(|_| limit())?;
    stack.push((root, 0usize));
    while let Some((value, next)) = stack.last_mut() {
        if *next == 0 {
            take(records, 1)?;
            m.slots::<Sample>(1)?;
            for text in [&value.name, &value.organism, &value.number, &value.comment] {
                m.text(text)?;
            }
            m.meta(&value.metadata)?;
            m.slots::<Sample>(value.subsamples.len())?;
            if value.subsamples.len() > *records {
                return Err(limit());
            }
        }
        if *next == value.subsamples.len() {
            stack.pop();
            continue;
        }
        let child = &value.subsamples[*next];
        *next += 1;
        if stack.len() >= max_depth {
            return Err(limit());
        }
        m.charge(1, 0)?;
        if stack.len() == stack.capacity() {
            let capacity = stack
                .capacity()
                .checked_mul(2)
                .ok_or_else(limit)?
                .min(max_depth);
            m.slots::<(&Sample, usize)>(capacity)?;
            stack
                .try_reserve_exact(capacity - stack.len())
                .map_err(|_| limit())?;
        }
        stack.push((child, 0));
    }
    Ok(())
}
fn instrument(m: &mut Meter<'_>, i: &Instrument, records: &mut usize) -> Result<()> {
    take(records, 1)?;
    m.slots::<Instrument>(1)?;
    for text in [
        &i.name,
        &i.vendor,
        &i.model,
        &i.customizations,
        &i.software.name,
        &i.software.version,
    ] {
        m.text(text)?;
    }
    m.meta(&i.metadata)?;
    m.cv(&i.software.cv_terms)?;
    take(records, i.ion_sources.len())?;
    take(records, i.mass_analyzers.len())?;
    take(records, i.ion_detectors.len())?;
    m.slots::<super::IonSource>(i.ion_sources.len())?;
    m.slots::<super::MassAnalyzer>(i.mass_analyzers.len())?;
    m.slots::<super::IonDetector>(i.ion_detectors.len())?;
    for value in &i.ion_sources {
        m.meta(&value.metadata)?;
    }
    for value in &i.mass_analyzers {
        m.meta(&value.metadata)?;
    }
    for value in &i.ion_detectors {
        m.meta(&value.metadata)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn repeated_copy_preflights_share_one_cumulative_allowance() {
        let settings = ExperimentalSettings {
            comment: "x".repeat(50_000),
            ..Default::default()
        };
        let (mut work, mut bytes) = (1_000_000, 100_000);
        settings.with_budget(&mut work, &mut bytes).unwrap();
        assert!(bytes < 50_000);
        assert!(settings.with_budget(&mut work, &mut bytes).is_err());
    }
}
