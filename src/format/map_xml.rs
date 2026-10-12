// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Shared feature/consensus XML identification and processing records.

use super::identification_xml::{self as xml, IDENTIFIER, Node, ReadOptions};
use crate::chemistry::ModificationsDB;
use crate::identification::{
    PeptideIdentification, ProteinHit, ProteinIdentification, SearchParameters,
};
use crate::metadata::{DataProcessing, MetaInfo, Software};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::atomic::{AtomicU64, Ordering};

const RANK: &str = "openms-rust:rank";
fn bad(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

/// Decode the final decimal component, preserving source invalid-text => zero.
/// Numeric overflow is a checked error instead of wrapping an assigned ID.
pub(crate) fn unique_id(text: &str) -> Result<u64> {
    let digits = text.rsplit('_').next().unwrap_or("");
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return Ok(0);
    }
    digits.parse().map_err(|_| bad("unique ID exceeds u64"))
}

#[derive(Default)]
pub(crate) struct ReadContext {
    pub protein_accessions: BTreeMap<String, String>,
    pub run_identifiers: BTreeMap<String, String>,
    xml_ids: BTreeSet<String>,
    identifiers: BTreeSet<String>,
}
impl ReadContext {
    fn add_id(&mut self, id: &str) -> Result<()> {
        xml::xml_id(id)?;
        if !self.xml_ids.insert(id.into()) {
            return Err(bad(format!("duplicate XML ID {id}")));
        }
        Ok(())
    }
}

pub(crate) fn read_run(
    node: &Node,
    options: &ReadOptions,
    context: &mut ReadContext,
    registry: &mut ModificationsDB,
    remaining_work: &mut usize,
    remaining_bytes: &mut usize,
) -> Result<ProteinIdentification> {
    node.check(
        &["id", "date", "search_engine", "search_engine_version"],
        &["SearchParameters", "ProteinIdentification"],
    )?;
    let id = node.get("id")?;
    context.add_id(id)?;
    let timestamp = node.get("date")?;
    // The source writes an unset `DateTime` as `UNSET_RUN_DATE`; it reads back
    // as unset.
    let date_time = if timestamp == UNSET_RUN_DATE {
        None
    } else {
        xml::date(timestamp)?;
        Some(timestamp.into())
    };
    let mut run = ProteinIdentification {
        search_engine: node.get("search_engine")?.into(),
        search_engine_version: node.get("search_engine_version")?.into(),
        date_time,
        ..Default::default()
    };
    for search in node
        .children
        .iter()
        .filter(|n| n.name == "SearchParameters")
    {
        // The source permits successive blocks; the last supplies the run's parameters.
        run.search_parameters = xml::read_search(search, options)?;
        super::modification_definitions::register_search_parameters_with_budget(
            &run.search_parameters,
            registry,
            remaining_work,
            remaining_bytes,
        )?;
    }
    for protein in node
        .children
        .iter()
        .filter(|n| n.name == "ProteinIdentification")
    {
        protein.check(
            &[
                "score_type",
                "higher_score_better",
                "significance_threshold",
            ],
            &["ProteinHit", "UserParam"],
        )?;
        run.score_type = protein.get("score_type")?.into();
        run.higher_score_better = xml::boolean(protein.get("higher_score_better")?)?;
        run.significance_threshold =
            xml::finite(protein.optional("significance_threshold").unwrap_or("0"))?;
        run.metadata.extend(xml::read_meta(protein, options)?);
        for hit in protein.children.iter().filter(|n| n.name == "ProteinHit") {
            hit.check(
                &["id", "accession", "score", "coverage", "sequence"],
                &["UserParam"],
            )?;
            let reference = hit.get("id")?;
            // The source's handler keeps `proteinid_to_accession_[id] = accession`
            // without checking uniqueness, so a protein-hit ID reused in a later
            // run replaces the earlier mapping, and every later `protein_refs`
            // resolves through the last one. Upstream fixtures depend on it:
            // `MzTabExporter_5_in.consensusXML` (test-data 0cb15f2) repeats
            // `PH_0` in each of its runs and the Release build's
            // MzTabExporter_5 test reads it. The ID syntax is still checked; the
            // accession map below is last-wins, as in the source.
            xml::xml_id(reference)?;
            let mut value = ProteinHit {
                accession: hit.get("accession")?.into(),
                score: xml::finite(hit.get("score")?)?,
                sequence: hit.optional("sequence").unwrap_or("").into(),
                coverage: hit.optional("coverage").map(xml::finite).transpose()?,
                metadata: xml::read_meta(hit, options)?,
                ..Default::default()
            };
            if value.coverage == Some(-1.0) {
                value.coverage = None;
            }
            if let Some(rank) = value.metadata.remove(RANK) {
                if rank.unit().is_some() {
                    return Err(bad("rank metadata cannot have a unit"));
                }
                value.rank = u32::try_from(rank.as_i64()?).map_err(|_| bad("rank exceeds u32"))?;
            }
            value.validate()?;
            context
                .protein_accessions
                .insert(reference.into(), value.accession.clone());
            run.hits.push(value);
        }
    }
    for (key, paths) in [
        ("spectra_data", &mut run.primary_ms_run_paths),
        ("spectra_data_raw", &mut run.raw_ms_run_paths),
    ] {
        if let Some(value) = run.metadata.remove(key) {
            if value.unit().is_some() {
                return Err(bad("run paths cannot have a unit"));
            }
            *paths = value.as_string_list()?.to_vec();
        }
    }
    run.identifier = if let Some(value) = run.metadata.remove(IDENTIFIER) {
        if value.unit().is_some() {
            return Err(bad("run identifier cannot have a unit"));
        }
        value.as_str()?.into()
    } else {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let next = NEXT
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |n| n.checked_add(1))
            .map_err(|_| bad("run identity space exhausted"))?;
        format!(
            "{}_{}_{}_{next}",
            run.search_engine,
            timestamp,
            std::process::id()
        )
    };
    if !context.identifiers.insert(run.identifier.clone()) {
        return Err(bad("duplicate identification run identifier"));
    }
    context
        .run_identifiers
        .insert(id.into(), run.identifier.clone());
    run.validate()?;
    Ok(run)
}

pub(crate) fn read_peptide(
    node: &Node,
    options: &ReadOptions,
    context: &ReadContext,
    registry: &ModificationsDB,
    remaining_work: &mut usize,
    remaining_bytes: &mut usize,
) -> Result<PeptideIdentification> {
    let reference = node.get("identification_run_ref")?;
    let identifier = context
        .run_identifiers
        .get(reference)
        .ok_or_else(|| bad(format!("unknown identification run {reference}")))?;
    xml::read_peptide_with_budget(
        node,
        identifier,
        &context.protein_accessions,
        options,
        registry,
        remaining_work,
        remaining_bytes,
    )
}

pub(crate) struct WriteContext {
    pub run_ids: BTreeMap<String, String>,
    pub protein_refs: BTreeMap<String, BTreeMap<String, String>>,
}
impl WriteContext {
    /// Caller preflights all identification payloads before constructing this index.
    pub fn new(runs: &[ProteinIdentification]) -> Result<Self> {
        let mut result = Self {
            run_ids: BTreeMap::new(),
            protein_refs: BTreeMap::new(),
        };
        let mut count = 0usize;
        for (index, run) in runs.iter().enumerate() {
            run.validate()?;
            if result
                .run_ids
                .insert(run.identifier.clone(), format!("PI_{index}"))
                .is_some()
            {
                return Err(bad("duplicate identification run identifier"));
            }
            let mut refs = BTreeMap::new();
            for hit in &run.hits {
                if refs
                    .insert(hit.accession.clone(), format!("PH_{count}"))
                    .is_some()
                {
                    return Err(bad("duplicate protein accession within a run"));
                }
                count = count
                    .checked_add(1)
                    .ok_or_else(|| bad("protein ID count overflow"))?;
            }
            result.protein_refs.insert(run.identifier.clone(), refs);
        }
        Ok(result)
    }
}

/// How the source's featureXML and consensusXML writers print the date of a
/// protein identification run whose `DateTime` is unset (`DateTime()` prints
/// its zero date and time), e.g. the default run `FeatureMap` producers add
/// with `getProteinIdentifications().resize(1)`. A run without a date is
/// written so, and reads back without one.
pub(crate) const UNSET_RUN_DATE: &str = "0000-00-00T00:00:00";

/// The dialect supplies group-encoded metadata and portable search definitions.
pub(crate) fn write_run(
    value: &ProteinIdentification,
    context: &WriteContext,
    search: &SearchParameters,
    metadata: &MetaInfo,
) -> Result<Node> {
    let timestamp = match value.date_time.as_deref() {
        Some(timestamp) => {
            xml::date(timestamp)?;
            timestamp
        }
        None => UNSET_RUN_DATE,
    };
    let mut run = Node::new("IdentificationRun");
    run.attr(
        "id",
        context
            .run_ids
            .get(&value.identifier)
            .ok_or_else(|| bad("unindexed protein run"))?,
    );
    run.attr("date", timestamp);
    run.attr("search_engine", &value.search_engine);
    run.attr("search_engine_version", &value.search_engine_version);
    let mut search = xml::write_search(search, "")?;
    search.attrs.remove("id");
    run.children.push(search);
    let mut protein = Node::new("ProteinIdentification");
    protein.attr("score_type", &value.score_type);
    protein.attr("higher_score_better", value.higher_score_better);
    protein.attr("significance_threshold", value.significance_threshold);
    let refs = context
        .protein_refs
        .get(&value.identifier)
        .ok_or_else(|| bad("unindexed protein run"))?;
    for hit in &value.hits {
        if !hit.modifications.is_empty() {
            return Err(Error::Unsupported(
                "map XML cannot encode ProteinHit modification positions".into(),
            ));
        }
        let mut child = Node::new("ProteinHit");
        child.attr(
            "id",
            refs.get(&hit.accession)
                .ok_or_else(|| bad("unindexed protein hit"))?,
        );
        child.attr("accession", &hit.accession);
        child.attr("sequence", &hit.sequence);
        child.attr("score", hit.score);
        if let Some(coverage) = hit.coverage {
            child.attr("coverage", coverage);
        }
        let mut meta = hit.metadata.clone();
        if meta.contains_key(RANK) {
            return Err(bad("reserved protein rank metadata collision"));
        }
        if hit.rank != 0 {
            meta.insert(RANK.into(), hit.rank.into());
        }
        xml::write_meta(&mut child, &meta)?;
        protein.children.push(child);
    }
    let mut meta = metadata.clone();
    if meta.contains_key(IDENTIFIER) {
        return Err(bad("reserved run identifier metadata collision"));
    }
    meta.insert(IDENTIFIER.into(), value.identifier.clone().into());
    for (key, paths) in [
        ("spectra_data", &value.primary_ms_run_paths),
        ("spectra_data_raw", &value.raw_ms_run_paths),
    ] {
        if meta.contains_key(key) {
            return Err(bad("run paths must use their typed fields"));
        }
        if !paths.is_empty() {
            meta.insert(key.into(), paths.clone().into());
        }
    }
    xml::write_meta(&mut protein, &meta)?;
    run.children.push(protein);
    Ok(run)
}

pub(crate) fn write_peptide(
    value: &PeptideIdentification,
    tag: &str,
    context: &WriteContext,
    registry: &ModificationsDB,
    remaining_work: &mut usize,
    remaining_bytes: &mut usize,
) -> Result<Node> {
    let run = context
        .run_ids
        .get(&value.identifier)
        .ok_or_else(|| bad("peptide identifier has no protein run"))?;
    let refs = context
        .protein_refs
        .get(&value.identifier)
        .ok_or_else(|| bad("peptide identifier has no protein references"))?;
    let mut node =
        xml::write_peptide_with_budget(value, refs, registry, remaining_work, remaining_bytes)?;
    node.name = tag.into();
    node.attr("identification_run_ref", run);
    Ok(node)
}

pub(crate) fn read_processing(node: &Node, options: &ReadOptions) -> Result<DataProcessing> {
    node.check(
        &["completion_time"],
        &["software", "processingAction", "UserParam"],
    )?;
    let mut software = node.children.iter().filter(|n| n.name == "software");
    // The source handlers (`ConsensusXMLHandler.cpp:580-592`, and the
    // `FeatureXMLHandler.cpp:416-427` alike) create the entry on `dataProcessing` and fill
    // the software only when a `software` child follows, so a missing one
    // leaves the default (empty) software; the retained class-test fixture
    // `Masstraces_Testdata.consensusXML` has such an entry.
    let source = software.next();
    if software.next().is_some() {
        return Err(bad("duplicate software"));
    }
    let software = match source {
        Some(source) => {
            source.check(&["name", "version"], &[])?;
            Software {
                name: source.get("name")?.into(),
                version: source.get("version")?.into(),
                ..Default::default()
            }
        }
        None => Software::default(),
    };
    let timestamp = node.get("completion_time")?;
    let mut value = DataProcessing {
        software,
        completion_time: Some(timestamp.parse()?),
        metadata: xml::read_meta(node, options)?,
        ..Default::default()
    };
    for action in node
        .children
        .iter()
        .filter(|n| n.name == "processingAction")
    {
        action.check(&["name"], &[])?;
        value.actions.insert(action.get("name")?.parse()?);
    }
    value.validate()?;
    Ok(value)
}

pub(crate) fn write_processing(value: &DataProcessing) -> Result<Node> {
    if !value.software.cv_terms.is_empty() || !value.software.cv_terms.metadata.is_empty() {
        return Err(Error::Unsupported(
            "map XML software CV terms and metadata".into(),
        ));
    }
    value.validate()?;
    let timestamp = value
        .completion_time
        .ok_or_else(|| bad("dataProcessing requires completion_time"))?;
    let mut node = Node::new("dataProcessing");
    if !timestamp.is_valid() || !(1..=9999).contains(&timestamp.date_components().2) {
        return Err(bad(
            "dataProcessing timestamp is invalid or outside supported XML years",
        ));
    }
    let text = timestamp.format(if timestamp.millisecond() == 0 {
        "yyyy-MM-ddThh:mm:ss"
    } else {
        "yyyy-MM-ddThh:mm:ss.zzz"
    })?;
    if crate::data_structures::DateTime::parse(&text)? != timestamp {
        return Err(bad(
            "dataProcessing timestamp has inconsistent partial fields",
        ));
    }
    node.attr("completion_time", text);
    let mut software = Node::new("software");
    software.attr("name", &value.software.name);
    software.attr("version", &value.software.version);
    node.children.push(software);
    for action in &value.actions {
        let mut child = Node::new("processingAction");
        child.attr("name", action.name());
        node.children.push(child);
    }
    xml::write_meta(&mut node, &value.metadata)?;
    Ok(node)
}
