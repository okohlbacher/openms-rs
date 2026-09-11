// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::{EnzymeTermSpecificity, ProteinGroup, ProteinIdentification, bad};
use crate::Result;
use crate::comparison::Tolerance;
use crate::data_structures::list::{self, ListFormat};
use crate::metadata::{MetaInfo, MetaValue, MetaValueData};
use std::collections::HashSet;

// Conservative descriptor/payload accounting for operations that create owned
// output. Ordinary borrowed lookups do not allocate or traverse protein hits.
const MAX_ITEMS: usize = 1_000_000;
const MAX_BYTES: usize = 64 * 1024 * 1024;

#[derive(Default)]
struct Measure {
    items: usize,
    bytes: usize,
}
impl Measure {
    fn add(&mut self, items: usize, bytes: usize) -> Result<()> {
        self.items = self.items.checked_add(items).ok_or_else(limit)?;
        self.bytes = self.bytes.checked_add(bytes).ok_or_else(limit)?;
        if self.items > MAX_ITEMS || self.bytes > MAX_BYTES {
            return Err(limit());
        }
        Ok(())
    }
    fn text(&mut self, text: &str) -> Result<()> {
        self.add(1, text.len().checked_add(32).ok_or_else(limit)?)
    }
    fn strings(&mut self, values: &[String]) -> Result<()> {
        self.add(values.len(), 0)?;
        for value in values {
            self.text(value)?;
        }
        Ok(())
    }
    fn value(&mut self, value: &MetaValue) -> Result<()> {
        self.add(1, 128)?;
        if let Some(unit) = value.unit() {
            self.text(unit.accession())?;
            self.text(unit.name())?;
            self.text(unit.cv_ref())?;
        }
        match value.data() {
            MetaValueData::String(value) => self.text(value),
            MetaValueData::StringList(values) => self.strings(values),
            MetaValueData::IntegerList(values) => self.add(
                values.len(),
                values.len().checked_mul(32).ok_or_else(limit)?,
            ),
            MetaValueData::FloatList(values) => self.add(
                values.len(),
                values.len().checked_mul(32).ok_or_else(limit)?,
            ),
            _ => Ok(()),
        }
    }
    fn metadata(&mut self, values: &MetaInfo) -> Result<()> {
        self.add(values.len(), 0)?;
        // Includes sparse BTreeMap node storage, not just populated key/value
        // slots. Budget up to two full nodes per entry before clone allocation.
        let per_entry = std::mem::size_of::<(String, MetaValue)>()
            .checked_mul(24)
            .and_then(|n| n.checked_add(256))
            .ok_or_else(limit)?;
        self.add(0, values.len().checked_mul(per_entry).ok_or_else(limit)?)?;
        for (key, value) in values {
            self.text(key)?;
            self.value(value)?;
        }
        Ok(())
    }
}
fn limit() -> crate::Error {
    bad("protein-run output exceeds its descriptor or byte limit")
}

impl ProteinIdentification {
    /// Source inference-engine heuristic, including Percolator only when an
    /// indistinguishable group exists. This does not inspect group contents.
    pub fn has_inference_engine_as_search_engine(&self) -> bool {
        matches!(
            self.search_engine.as_str(),
            "Fido" | "BayesianProteinInference" | "Epifany" | "ProteinInference"
        ) || (self.search_engine == "Percolator" && !self.indistinguishable_groups.is_empty())
    }

    /// Explicit metadata takes precedence, including an explicitly empty string.
    /// Non-string metadata is an error, matching the source's strict conversion.
    pub fn inference_engine(&self) -> Result<&str> {
        match self.search_parameters.metadata.get("InferenceEngine") {
            Some(value) => value.as_str(),
            None if self.has_inference_engine_as_search_engine() => Ok(&self.search_engine),
            None => Ok(""),
        }
    }
    pub fn set_inference_engine(&mut self, engine: impl Into<String>) {
        self.search_parameters
            .metadata
            .insert("InferenceEngine".into(), engine.into().into());
    }
    pub fn has_inference_data(&self) -> Result<bool> {
        Ok(!self.inference_engine()?.is_empty())
    }
    pub fn inference_engine_version(&self) -> Result<&str> {
        match self
            .search_parameters
            .metadata
            .get("InferenceEngineVersion")
        {
            Some(value) => value.as_str(),
            None if self.has_inference_data()? => Ok(&self.search_engine_version),
            None => Ok(""),
        }
    }
    pub fn set_inference_engine_version(&mut self, version: impl Into<String>) {
        self.search_parameters
            .metadata
            .insert("InferenceEngineVersion".into(), version.into().into());
    }

    /// Recover the first lexical SE: key for Percolator/ConsensusID runs.
    /// Filtering is case-sensitive and key values are unused. Native metadata
    /// uses lexical keys; C++ visits numeric registry IDs, so ambiguous runs can
    /// select a different first engine.
    pub fn original_search_engine_name(&self) -> &str {
        if !self.search_engine.contains("Percolator") && !self.search_engine.contains("ConsensusID")
        {
            return &self.search_engine;
        }
        self.search_parameters
            .metadata
            .keys()
            .find(|key| key.starts_with("SE:") && !key.contains("percolator"))
            .map(|key| &key[3..])
            .unwrap_or("Unknown")
    }

    /// Compare engine/version and the existing source-compatible search settings.
    /// Native errors from malformed settings propagate; source warning logs are
    /// represented by the returned boolean.
    pub fn peptide_ids_mergeable(&self, other: &Self, experiment_type: &str) -> Result<bool> {
        let settings_match = self
            .search_parameters
            .mergeable(&other.search_parameters, experiment_type)?;
        Ok(self.search_engine == other.search_engine
            && self.search_engine_version == other.search_engine_version
            && settings_match)
    }

    /// Ordered export settings, with lenient source numeric/list formatting.
    /// An empty engine selects the twelve standard fields even for merged runs.
    /// Other engines select metadata with the source's literal prefix rule and
    /// remove that prefix plus one byte (normally ':'). Splitting a UTF-8 scalar
    /// is an error. Output is bounded to one million descriptors and 64 MiB of
    /// conservatively counted input/output payload; units are not printed.
    pub fn search_engine_settings_as_pairs(&self, engine: &str) -> Result<Vec<(String, String)>> {
        let mut measure = Measure::default();
        measure.text(engine)?;
        let mut output = Vec::new();
        let params = &self.search_parameters;
        if engine.is_empty()
            || (self.search_engine == engine
                && engine != "Percolator"
                && !engine.starts_with("ConsensusID"))
        {
            for text in [
                &params.database,
                &params.database_version,
                &params.digestion_enzyme,
                &params.charges,
            ] {
                measure.text(text)?;
            }
            measure.strings(&params.fixed_modifications)?;
            measure.strings(&params.variable_modifications)?;
            // Charge owned output before constructing the twelve-value array.
            // List concatenation also retains intermediate converted items.
            for text in [
                &params.database,
                &params.database_version,
                &params.digestion_enzyme,
                &params.charges,
            ] {
                measure.text(text)?;
            }
            for _ in 0..3 {
                measure.strings(&params.fixed_modifications)?;
                measure.strings(&params.variable_modifications)?;
            }
            measure.add(24, 2048)?; // fixed field names, numeric strings and pair slots
            let tolerance = |value| {
                let (Tolerance::Absolute(number) | Tolerance::Ppm(number)) = value;
                (
                    crate::param::value::format_float(number, true),
                    if matches!(value, Tolerance::Ppm(_)) {
                        "ppm"
                    } else {
                        "Da"
                    },
                )
            };
            let (fragment, fragment_unit) = tolerance(params.fragment_tolerance);
            let (precursor, precursor_unit) = tolerance(params.precursor_tolerance);
            let specificity = match params.enzyme_specificity {
                EnzymeTermSpecificity::Unknown => "unknown",
                EnzymeTermSpecificity::Full => "full",
                EnzymeTermSpecificity::Semi => "semi",
                EnzymeTermSpecificity::None => "none",
            };
            for (key, value) in [
                ("db", params.database.clone()),
                ("db_version", params.database_version.clone()),
                ("fragment_mass_tolerance", fragment),
                ("fragment_mass_tolerance_unit", fragment_unit.into()),
                ("precursor_mass_tolerance", precursor),
                ("precursor_mass_tolerance_unit", precursor_unit.into()),
                ("enzyme", params.digestion_enzyme.clone()),
                ("enzyme_term_specificity", specificity.into()),
                ("charges", params.charges.clone()),
                ("missed_cleavages", params.missed_cleavages.to_string()),
                (
                    "fixed_modifications",
                    list::concatenate(&params.fixed_modifications, ",")?,
                ),
                (
                    "variable_modifications",
                    list::concatenate(&params.variable_modifications, ",")?,
                ),
            ] {
                output.push((key.into(), value));
            }
        } else {
            measure.add(params.metadata.len(), 0)?;
            for (key, value) in &params.metadata {
                measure.text(key)?;
                if key.starts_with(engine) {
                    let start = (engine.len() + 1).min(key.len());
                    let suffix = key
                        .get(start..)
                        .ok_or_else(|| bad("search settings prefix splits a UTF-8 character"))?;
                    measure.value(value)?;
                    measure.text(suffix)?;
                    // Precharge source list-format conversion, joining and
                    // bracket-wrapping before its temporary strings exist.
                    for _ in 0..3 {
                        measure.value(value)?;
                    }
                    let value = value.to_list_text()?;
                    output.push((suffix.into(), value.into_owned()));
                }
            }
        }
        Ok(output)
    }

    /// Append singleton groups in hit order, taking the first ungrouped hit's
    /// score for each accession. Existing groups and sample arrays are retained.
    /// Only consumed scores must be finite; failures leave the run unchanged.
    pub fn fill_indistinguishable_groups_with_singletons(&mut self) -> Result<usize> {
        let mut measure = Measure::default();
        measure.add(self.hits.len(), 0)?;
        measure.add(self.indistinguishable_groups.len(), 0)?;
        let mut grouped = HashSet::new();
        for group in &self.indistinguishable_groups {
            measure.add(group.accessions.len(), 0)?;
            for accession in &group.accessions {
                measure.text(accession)?;
                measure.add(0, 128)?;
                grouped.insert(accession.as_str());
            }
        }
        let mut next = Vec::new();
        for hit in &self.hits {
            measure.text(&hit.accession)?;
            measure.add(0, 128)?;
            if grouped.insert(&hit.accession) {
                super::finite(hit.score, "singleton protein score")?;
                measure.text(&hit.accession)?;
                measure.add(1, std::mem::size_of::<ProteinGroup>())?;
                next.push(ProteinGroup {
                    probability: hit.score,
                    accessions: vec![hit.accession.clone()],
                    ..Default::default()
                });
            }
        }
        let added = next.len();
        let required = self
            .indistinguishable_groups
            .len()
            .checked_add(added)
            .ok_or_else(limit)?;
        if required > self.indistinguishable_groups.capacity() {
            measure.add(
                0,
                required
                    .checked_mul(std::mem::size_of::<ProteinGroup>())
                    .ok_or_else(limit)?,
            )?;
        }
        self.indistinguishable_groups
            .try_reserve_exact(added)
            .map_err(|_| limit())?;
        self.indistinguishable_groups.extend(next);
        Ok(added)
    }

    /// Copy run metadata and search settings, retaining this run's hits and both
    /// protein-group collections. Native dedicated primary/raw path fields are
    /// copied too: the source stores those paths inside its metadata map.
    /// Bounds are checked before copying; source hit/group payload is not read.
    pub fn copy_metadata_only(&mut self, source: &Self) -> Result<()> {
        let mut measure = Measure::default();
        let p = &source.search_parameters;
        for text in [
            &source.identifier,
            &source.search_engine,
            &source.search_engine_version,
            &source.score_type,
            &p.database,
            &p.database_version,
            &p.taxonomy,
            &p.charges,
            &p.digestion_enzyme,
            &p.digestion_regex,
        ] {
            measure.text(text)?;
        }
        if let Some(date) = &source.date_time {
            measure.text(date)?;
        }
        for list in [
            &source.primary_ms_run_paths,
            &source.raw_ms_run_paths,
            &p.fixed_modifications,
            &p.variable_modifications,
        ] {
            measure.strings(list)?;
        }
        measure.metadata(&source.metadata)?;
        measure.metadata(&p.metadata)?;
        let mut next = Self {
            identifier: source.identifier.clone(),
            search_engine: source.search_engine.clone(),
            search_engine_version: source.search_engine_version.clone(),
            search_parameters: p.clone(),
            date_time: source.date_time.clone(),
            score_type: source.score_type.clone(),
            higher_score_better: source.higher_score_better,
            significance_threshold: source.significance_threshold,
            primary_ms_run_paths: source.primary_ms_run_paths.clone(),
            raw_ms_run_paths: source.raw_ms_run_paths.clone(),
            metadata: source.metadata.clone(),
            ..Default::default()
        };
        next.hits = std::mem::take(&mut self.hits);
        next.protein_groups = std::mem::take(&mut self.protein_groups);
        next.indistinguishable_groups = std::mem::take(&mut self.indistinguishable_groups);
        *self = next;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn measure_rejects_cumulative_and_overflowing_payload() {
        let mut count = Measure::default();
        count.add(MAX_ITEMS, MAX_BYTES).unwrap();
        assert!(count.add(1, 0).is_err());
        assert!(
            Measure {
                items: 0,
                bytes: MAX_BYTES
            }
            .text("")
            .is_err()
        );
        assert!(Measure { items: 1, bytes: 0 }.add(usize::MAX, 0).is_err());
        assert!(Measure { items: 0, bytes: 1 }.add(0, usize::MAX).is_err());
    }
}
