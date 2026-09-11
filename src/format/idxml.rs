// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded native idXML 1.5 interchange. See `docs/IDXML_SUPPORT.md` for the
//! supported source encodings and explicit errors for unrepresentable state.

use super::identification_xml::*;
pub use super::identification_xml::{ReadOptions, WriteOptions};
use super::modification_definitions;
use crate::Result;
use crate::chemistry::ModificationsDB;
use crate::identification::{
    PeptideIdentification, ProteinGroup, ProteinHit, ProteinIdentification, SearchParameters,
};
use crate::metadata::MetaInfo;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, Write};
use std::path::Path;

/// Flat native records linked by run identifier, plus otherwise unused search
/// parameter blocks. XML IDs are transport references and are regenerated.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IdXmlDocument {
    pub document_id: String,
    pub protein_identifications: Vec<ProteinIdentification>,
    pub peptide_identifications: Vec<PeptideIdentification>,
    pub unreferenced_search_parameters: Vec<SearchParameters>,
}

fn groups(
    meta: &mut MetaInfo,
    prefix: &str,
    refs: &BTreeMap<String, String>,
    options: &ReadOptions,
    work: &mut usize,
    bytes: &mut usize,
) -> Result<Vec<ProteinGroup>> {
    let mut result = Vec::new();
    for index in 0.. {
        let key = format!("{prefix}_{index}");
        let Some(value) = take_text(meta, &key)? else {
            break;
        };
        let mut fields = value.split(',');
        let probability = finite(fields.next().unwrap_or(""))?;
        let mut accessions = Vec::new();
        for id in fields {
            if accessions.len() == options.max_list_items {
                return Err(bad("protein group accession list exceeds limit"));
            }
            let accession = refs
                .get(id)
                .ok_or_else(|| bad(format!("unknown protein group reference {id}")))?;
            *work = work
                .checked_sub(accession.len().saturating_mul(4).saturating_add(1))
                .ok_or_else(|| bad("protein group reference work limit exceeded"))?;
            *bytes = bytes
                .checked_sub(accession.len().saturating_mul(4).saturating_add(128))
                .ok_or_else(|| bad("protein group reference payload limit exceeded"))?;
            accessions.push(accession.clone());
        }
        if accessions.is_empty() {
            return Err(bad("protein group has no protein references"));
        }
        result.push(ProteinGroup {
            probability,
            accessions,
            ..Default::default()
        });
    }
    if meta
        .keys()
        .any(|key| key.starts_with(&format!("{prefix}_")))
    {
        return Err(bad("protein group indices must be contiguous from zero"));
    }
    Ok(result)
}
fn paths(meta: &mut MetaInfo, key: &str) -> Result<Vec<String>> {
    meta.remove(key)
        .map(|v| v.as_string_list().map(<[String]>::to_vec))
        .transpose()
        .map(|v| v.unwrap_or_default())
}
fn read_protein(
    node: &Node,
    mut run: ProteinIdentification,
    refs: &mut BTreeMap<String, String>,
    all_ids: &mut BTreeSet<String>,
    options: &ReadOptions,
    work: &mut usize,
    bytes: &mut usize,
) -> Result<ProteinIdentification> {
    measure_node(node, work, bytes)?;
    node.check(
        &[
            "score_type",
            "higher_score_better",
            "significance_threshold",
        ],
        &["ProteinHit", "UserParam"],
    )?;
    run.score_type = node.get("score_type")?.into();
    run.higher_score_better = boolean(node.get("higher_score_better")?)?;
    run.significance_threshold = finite(node.optional("significance_threshold").unwrap_or("0"))?;
    run.metadata = read_meta(node, options)?;
    for child in node.children.iter().filter(|c| c.name == "ProteinHit") {
        child.check(
            &["id", "accession", "score", "coverage", "sequence"],
            &["UserParam"],
        )?;
        let id = child.get("id")?;
        xml_id(id)?;
        if !all_ids.insert(id.into()) {
            return Err(bad(format!("duplicate XML ID {id}")));
        }
        let accession = child.get("accession")?.to_owned();
        refs.insert(id.into(), accession.clone());
        let mut hit = ProteinHit {
            accession,
            score: finite(child.get("score")?)?,
            sequence: child.optional("sequence").unwrap_or("").into(),
            coverage: child.optional("coverage").map(finite).transpose()?,
            metadata: read_meta(child, options)?,
            ..Default::default()
        };
        // The source has a -1 unknown coverage sentinel; native records use None.
        if hit.coverage == Some(-1.0) {
            hit.coverage = None;
        }
        hit.rank = rank(&mut hit.metadata)?;
        hit.validate()?;
        run.hits.push(hit);
    }
    run.protein_groups = groups(
        &mut run.metadata,
        "protein_group",
        refs,
        options,
        work,
        bytes,
    )?;
    run.indistinguishable_groups = groups(
        &mut run.metadata,
        "indistinguishable_proteins",
        refs,
        options,
        work,
        bytes,
    )?;
    run.primary_ms_run_paths = paths(&mut run.metadata, "spectra_data")?;
    run.raw_ms_run_paths = paths(&mut run.metadata, "spectra_data_raw")?;
    if let Some(identifier) = take_text(&mut run.metadata, IDENTIFIER)? {
        run.identifier = identifier;
    }
    run.validate()?;
    Ok(run)
}
pub fn read(reader: impl BufRead) -> Result<IdXmlDocument> {
    read_with_options(reader, &ReadOptions::default())
}
/// Reads a bounded XML tree and returns native records only after full validation.
/// Unknown tags/attributes and unsupported encodings are errors, never dropped.
pub fn read_with_options(reader: impl BufRead, options: &ReadOptions) -> Result<IdXmlDocument> {
    read_with_registry(reader, options, ModificationsDB::global())
}
/// Read with caller-supplied modification names. Returned sequences retain owned
/// handles, so the registry may be dropped immediately after successful reading.
pub fn read_with_registry(
    reader: impl BufRead,
    options: &ReadOptions,
    registry: &ModificationsDB,
) -> Result<IdXmlDocument> {
    let limits = XmlLimits::default();
    let mut work = limits.max_work;
    let mut bytes = limits.max_payload_bytes;
    let root = parse_xml_with_budget(
        reader,
        options,
        limits.max_depth,
        None,
        &mut work,
        &mut bytes,
    )?;
    let mut registry = clone_registry(registry, &mut work, &mut bytes)?;
    if root.name != "IdXML" {
        return Err(bad("expected IdXML root"));
    }
    root.check(
        &[
            "version",
            "id",
            "xmlns:xsi",
            "xsi:noNamespaceSchemaLocation",
        ],
        &["SearchParameters", "IdentificationRun"],
    )?;
    let version = finite(root.optional("version").unwrap_or("1.0"))?;
    if !(1.0..=1.5).contains(&version) {
        return Err(unsupported("idXML versions outside 1.0 through 1.5"));
    }
    if root.optional("xsi:noNamespaceSchemaLocation").is_some()
        && root.optional("xmlns:xsi").is_none()
    {
        return Err(bad("unbound xsi attribute prefix"));
    }
    if root
        .optional("xmlns:xsi")
        .is_some_and(|v| v != "http://www.w3.org/2001/XMLSchema-instance")
    {
        return Err(bad("invalid xsi namespace"));
    }
    let mut document = IdXmlDocument {
        document_id: root.optional("id").unwrap_or("").into(),
        ..Default::default()
    };
    let mut search = BTreeMap::new();
    let mut search_order = Vec::new();
    let mut ids = BTreeSet::new();
    let mut used = BTreeSet::new();
    let mut run_ids = BTreeSet::new();
    for node in &root.children {
        if node.name == "SearchParameters" {
            let id = node.get("id")?;
            xml_id(id)?;
            if !ids.insert(id.into()) {
                return Err(bad(format!("duplicate XML ID {id}")));
            }
            let parameters = read_search(node, options)?;
            modification_definitions::register_search_parameters_with_budget(
                &parameters,
                &mut registry,
                &mut work,
                &mut bytes,
            )?;
            search.insert(id.to_owned(), parameters);
            search_order.push(id.to_owned());
        } else {
            node.check(
                &[
                    "date",
                    "search_engine",
                    "search_engine_version",
                    "search_parameters_ref",
                ],
                &["ProteinIdentification", "PeptideIdentification"],
            )?;
            let reference = node.get("search_parameters_ref")?;
            let parameters = search
                .get(reference)
                .ok_or_else(|| bad(format!("unknown search parameters {reference}")))?;
            measure_search(parameters, &mut work, &mut bytes)?;
            let parameters = parameters.clone();
            used.insert(reference.to_owned());
            let date_time = node.get("date")?;
            date(date_time)?;
            let engine = node.get("search_engine")?;
            let mut run = ProteinIdentification {
                identifier: format!(
                    "{engine}_{date_time}_{}",
                    document.protein_identifications.len()
                ),
                search_engine: engine.into(),
                search_engine_version: node.get("search_engine_version")?.into(),
                search_parameters: parameters,
                date_time: Some(date_time.into()),
                ..Default::default()
            };
            let proteins: Vec<_> = node
                .children
                .iter()
                .filter(|c| c.name == "ProteinIdentification")
                .collect();
            if proteins.len() > 1 {
                return Err(unsupported(
                    "multiple ProteinIdentification blocks in one run",
                ));
            }
            let mut refs = BTreeMap::new();
            if let Some(protein) = proteins.first() {
                run = read_protein(
                    protein, run, &mut refs, &mut ids, options, &mut work, &mut bytes,
                )?;
            }
            if !run_ids.insert(run.identifier.clone()) {
                return Err(bad("duplicate identification run identifier"));
            }
            for peptide in node
                .children
                .iter()
                .filter(|c| c.name == "PeptideIdentification")
            {
                document
                    .peptide_identifications
                    .push(read_peptide_with_budget(
                        peptide,
                        &run.identifier,
                        &refs,
                        options,
                        &registry,
                        &mut work,
                        &mut bytes,
                    )?);
            }
            document.protein_identifications.push(run);
        }
    }
    if document.protein_identifications.is_empty() {
        return Err(bad("idXML needs at least one IdentificationRun"));
    }
    for id in search_order {
        if !used.contains(&id) {
            document
                .unreferenced_search_parameters
                .push(search.remove(&id).unwrap());
        }
    }
    Ok(document)
}

fn write_groups(
    meta: &mut MetaInfo,
    groups: &[ProteinGroup],
    prefix: &str,
    refs: &BTreeMap<String, String>,
) -> Result<()> {
    if meta
        .keys()
        .any(|key| key.starts_with(&format!("{prefix}_")))
    {
        return Err(bad("reserved protein group metadata collision"));
    }
    for (index, group) in groups.iter().enumerate() {
        if !group.float_data_arrays.is_empty()
            || !group.integer_data_arrays.is_empty()
            || !group.string_data_arrays.is_empty()
        {
            return Err(unsupported(
                "idXML protein groups do not encode sample arrays",
            ));
        }
        if group.accessions.is_empty() {
            return Err(bad("idXML protein group needs accessions"));
        }
        let mut values = vec![group.probability.to_string()];
        for accession in &group.accessions {
            values.push(refs.get(accession).cloned().ok_or_else(|| {
                bad(format!(
                    "protein group accession absent from run: {accession}"
                ))
            })?);
        }
        meta.insert(format!("{prefix}_{index}"), values.join(",").into());
    }
    Ok(())
}
fn write_protein(
    value: &ProteinIdentification,
    count: &mut usize,
) -> Result<(Node, BTreeMap<String, String>)> {
    let mut node = Node::new("ProteinIdentification");
    node.attr("score_type", &value.score_type);
    node.attr("higher_score_better", value.higher_score_better);
    node.attr("significance_threshold", value.significance_threshold);
    let mut refs = BTreeMap::new();
    for hit in &value.hits {
        if !hit.modifications.is_empty() {
            return Err(unsupported(
                "idXML does not encode ProteinHit modification positions",
            ));
        }
        let id = format!("PH_{count}");
        *count = count
            .checked_add(1)
            .ok_or_else(|| bad("protein reference count overflow"))?;
        if refs.insert(hit.accession.clone(), id.clone()).is_some() {
            return Err(bad("duplicate protein accession within run is ambiguous"));
        }
        let mut child = Node::new("ProteinHit");
        child.attr("id", id);
        child.attr("accession", &hit.accession);
        child.attr("sequence", &hit.sequence);
        child.attr("score", hit.score);
        if let Some(coverage) = hit.coverage {
            child.attr("coverage", coverage);
        }
        let mut meta = hit.metadata.clone();
        if meta.contains_key(RANK) {
            return Err(bad("reserved rank metadata collision"));
        }
        if hit.rank != 0 {
            meta.insert(RANK.into(), hit.rank.into());
        }
        write_meta(&mut child, &meta)?;
        node.children.push(child);
    }
    let mut meta = value.metadata.clone();
    if meta.contains_key(IDENTIFIER) {
        return Err(bad("reserved run identifier metadata collision"));
    }
    meta.insert(IDENTIFIER.into(), value.identifier.clone().into());
    for (key, paths) in [
        ("spectra_data", &value.primary_ms_run_paths),
        ("spectra_data_raw", &value.raw_ms_run_paths),
    ] {
        if meta.contains_key(key) {
            return Err(bad(format!(
                "{key} metadata must use the typed run-path field"
            )));
        }
        if !paths.is_empty() {
            meta.insert(key.into(), paths.clone().into());
        }
    }
    write_groups(&mut meta, &value.protein_groups, "protein_group", &refs)?;
    write_groups(
        &mut meta,
        &value.indistinguishable_groups,
        "indistinguishable_proteins",
        &refs,
    )?;
    write_meta(&mut node, &meta)?;
    Ok((node, refs))
}
pub fn write(writer: impl Write, document: &IdXmlDocument) -> Result<()> {
    write_with_options(writer, document, &WriteOptions::default())
}
/// Builds and validates the complete output before touching the caller's writer.
/// An I/O error during the final write can still leave a partial external file.
pub fn write_with_options(
    writer: impl Write,
    document: &IdXmlDocument,
    options: &WriteOptions,
) -> Result<()> {
    write_with_registry(writer, document, options, ModificationsDB::global())
}
/// Validate sequence reconstruction against an explicitly supplied registry before
/// writing. Named Defined records are embedded; vocabulary records require the
/// corresponding vocabulary in the reading registry.
pub fn write_with_registry(
    mut writer: impl Write,
    document: &IdXmlDocument,
    options: &WriteOptions,
    registry: &ModificationsDB,
) -> Result<()> {
    if document.protein_identifications.is_empty() {
        return Err(bad("idXML writer needs at least one protein run"));
    }
    if options.max_records == 0 || options.max_xml_bytes == 0 {
        return Err(bad("idXML write limits must be positive"));
    }
    let mut work = XmlLimits::default().max_work;
    let mut bytes = XmlLimits::default().max_payload_bytes;
    measure_identifications(
        &document.protein_identifications,
        &document.peptide_identifications,
        &mut work,
        &mut bytes,
    )?;
    for search in &document.unreferenced_search_parameters {
        measure_search(search, &mut work, &mut bytes)?;
    }
    let definitions = modification_definitions::collect_iter_with_budget(
        &document.protein_identifications,
        &document.peptide_identifications,
        registry,
        &mut work,
        &mut bytes,
    )?;
    let encoded = modification_definitions::encode_by_run_with_budget(
        &document.protein_identifications,
        &definitions,
        &mut work,
        &mut bytes,
    )?;
    let mut registry = clone_registry(registry, &mut work, &mut bytes)?;
    let mut searches = Vec::new();
    for run in &document.protein_identifications {
        let mut search = run.search_parameters.clone();
        if let Some(records) = encoded.get(&run.identifier) {
            search.metadata.insert(
                modification_definitions::METADATA_KEY.into(),
                records.as_str().into(),
            );
        }
        modification_definitions::register_search_parameters_with_budget(
            &search,
            &mut registry,
            &mut work,
            &mut bytes,
        )?;
        searches.push(search);
    }
    for search in &document.unreferenced_search_parameters {
        modification_definitions::register_search_parameters_with_budget(
            search,
            &mut registry,
            &mut work,
            &mut bytes,
        )?;
    }
    let mut run_ids = BTreeSet::new();
    for run in &document.protein_identifications {
        run.validate()?;
        if !run_ids.insert(run.identifier.as_str()) {
            return Err(bad("duplicate identification run identifier"));
        }
        date(
            run.date_time
                .as_deref()
                .ok_or_else(|| bad("IdentificationRun requires date_time"))?,
        )?;
    }
    let mut peptides: BTreeMap<&str, Vec<&PeptideIdentification>> = BTreeMap::new();
    for peptide in &document.peptide_identifications {
        if !run_ids.contains(peptide.identifier.as_str()) {
            return Err(bad("peptide identifier does not match any protein run"));
        }
        peptides
            .entry(&peptide.identifier)
            .or_default()
            .push(peptide);
    }
    let mut root = Node::new("IdXML");
    root.attr("version", "1.5");
    root.attr("id", &document.document_id);
    root.attr("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    root.attr(
        "xsi:noNamespaceSchemaLocation",
        "https://www.openms.de/xml-schema/IdXML_1_5.xsd",
    );
    // Keep one parameter block per run: unlike C++ equality-based deduplication,
    // this cannot collapse search parameter metadata from distinct runs.
    for (i, search) in searches
        .iter()
        .chain(&document.unreferenced_search_parameters)
        .enumerate()
    {
        root.children
            .push(write_search(search, &format!("SP_{i}"))?);
    }
    let mut protein_count = 0;
    for (i, run) in document.protein_identifications.iter().enumerate() {
        let mut node = Node::new("IdentificationRun");
        node.attr("date", run.date_time.as_deref().unwrap());
        node.attr("search_engine", &run.search_engine);
        node.attr("search_engine_version", &run.search_engine_version);
        node.attr("search_parameters_ref", format!("SP_{i}"));
        let (protein, refs) = write_protein(run, &mut protein_count)?;
        node.children.push(protein);
        for peptide in peptides.get(run.identifier.as_str()).into_iter().flatten() {
            node.children.push(write_peptide_with_budget(
                peptide, &refs, &registry, &mut work, &mut bytes,
            )?);
        }
        root.children.push(node);
    }
    let mut output = Output {
        bytes: Vec::new(),
        records: 0,
        options: *options,
    };
    output.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")?;
    output.node(&root, 0)?;
    writer.write_all(&output.bytes)?;
    writer.flush()?;
    Ok(())
}

/// Load plain or magic-detected gzip/bzip2 idXML with default limits and registry.
/// Compressed input requires `file-compression`; stream reads remain uncompressed.
pub fn load(path: impl AsRef<Path>) -> Result<IdXmlDocument> {
    load_with_options(path, &ReadOptions::default())
}

pub fn load_with_options(path: impl AsRef<Path>, options: &ReadOptions) -> Result<IdXmlDocument> {
    load_with_registry(path, options, ModificationsDB::global())
}

/// Return an owned document only after decompression, parsing and validation.
pub fn load_with_registry(
    path: impl AsRef<Path>,
    options: &ReadOptions,
    registry: &ModificationsDB,
) -> Result<IdXmlDocument> {
    read_with_registry(super::path_io::open(path.as_ref())?, options, registry)
}

/// Replace a destination only after the complete load succeeds.
pub fn load_into(path: impl AsRef<Path>, destination: &mut IdXmlDocument) -> Result<()> {
    let document = load(path)?;
    *destination = document;
    Ok(())
}

/// Atomically publish plain idXML, regardless of the filename's compression suffix.
/// This matches the source IdXMLFile writer; it does not use XMLFile::save_.
pub fn store(path: impl AsRef<Path>, document: &IdXmlDocument) -> Result<()> {
    store_with_options(path, document, &WriteOptions::default())
}

pub fn store_with_options(
    path: impl AsRef<Path>,
    document: &IdXmlDocument,
    options: &WriteOptions,
) -> Result<()> {
    store_with_registry(path, document, options, ModificationsDB::global())
}

/// Validate chemistry against the supplied registry before publishing output.
/// A failure preserves an existing destination and removes the temporary output.
pub fn store_with_registry(
    path: impl AsRef<Path>,
    document: &IdXmlDocument,
    options: &WriteOptions,
    registry: &ModificationsDB,
) -> Result<()> {
    let path = path.as_ref();
    let filename = path
        .to_str()
        .ok_or_else(|| crate::Error::InvalidValue("output filename must be UTF-8".into()))?;
    if !super::file_types::has_valid_extension(filename, super::FileType::IdXml) {
        return Err(crate::Error::InvalidValue(
            "invalid idXML output extension".into(),
        ));
    }
    super::path_io::write_plain(path, |writer| {
        write_with_registry(writer, document, options, registry)
    })
}
