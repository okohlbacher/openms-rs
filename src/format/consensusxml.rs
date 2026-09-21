// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Native consensusXML 1.7 interchange, including protein-group quantities.
//! See `docs/CONSENSUSXML_SUPPORT.md` for source conventions and native corrections.
//!
//! The source `ConsensusXMLFile` and its handler derive from `ProgressLogger`;
//! [`load_with_progress`] and [`store_with_progress`] make the handler's
//! progress calls, and every other entry point runs the same code and reports
//! nothing.

mod groups;
use super::identification_xml::{self as xml, Node};
use super::{FileType, map_xml as common};
use crate::chemistry::ModificationsDB;
use crate::concept::progress_logger::{ProgressLogger, ProgressReporter};
use crate::identification::PeptideIdentification;
use crate::kernel::{ColumnHeader, ConsensusFeature, ConsensusMap, FeatureHandle};
use crate::{Error, Result};
use std::io::{BufRead, Write};
use std::ops::Range;
use std::path::Path;

/// The three range options the source handler consumes, plus native resource
/// ceilings the source does not have.
#[derive(Clone, Debug)]
pub struct ReadOptions {
    /// Source filters include the lower endpoint and exclude the upper endpoint.
    pub rt_range: Option<Range<f64>>,
    /// Half-open m/z filter on the consensus centroid.
    pub mz_range: Option<Range<f64>>,
    /// Half-open intensity filter on the consensus centroid.
    pub intensity_range: Option<Range<f64>>,
    /// Maximum encoded input bytes.
    pub max_xml_bytes: u64,
    /// Maximum number of records.
    pub max_records: usize,
    /// Maximum number of items in one list.
    pub max_list_items: usize,
    /// Maximum decoded payload bytes.
    pub max_payload_bytes: usize,
    /// Maximum parse and conversion work.
    pub max_work: usize,
    /// Opt into source behavior: discard stale quantities and expose a warning
    /// through `read_report`. The default rejects this loss transactionally.
    pub discard_mismatched_quantities: bool,
}
impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            rt_range: None,
            mz_range: None,
            intensity_range: None,
            max_xml_bytes: 64 * 1024 * 1024,
            max_records: 1_000_000,
            max_list_items: 1_000_000,
            max_payload_bytes: 256 * 1024 * 1024,
            max_work: 50_000_000,
            discard_mismatched_quantities: false,
        }
    }
}
impl ReadOptions {
    fn xml(&self) -> xml::ReadOptions {
        xml::ReadOptions {
            max_xml_bytes: self.max_xml_bytes,
            max_records: self.max_records,
            max_list_items: self.max_list_items,
        }
    }
    fn validate(&self) -> Result<()> {
        for range in [&self.rt_range, &self.mz_range, &self.intensity_range]
            .into_iter()
            .flatten()
        {
            if !range.start.is_finite() || !range.end.is_finite() || range.start > range.end {
                return Err(bad("invalid half-open consensus range"));
            }
        }
        if self.max_xml_bytes == 0
            || self.max_records == 0
            || self.max_list_items == 0
            || self.max_payload_bytes == 0
            || self.max_work == 0
        {
            return Err(bad("consensusXML limits must be positive"));
        }
        Ok(())
    }
    fn accepts(&self, feature: &ConsensusFeature) -> bool {
        [&self.rt_range, &self.mz_range, &self.intensity_range]
            .into_iter()
            .zip([feature.rt, feature.mz, f64::from(feature.intensity)])
            .all(|(range, value)| range.as_ref().is_none_or(|r| r.contains(&value)))
    }
}
/// Native output ceilings; the source writer has none.
#[derive(Clone, Copy, Debug)]
pub struct WriteOptions {
    /// Maximum encoded output bytes.
    pub max_xml_bytes: usize,
    /// Maximum number of records.
    pub max_records: usize,
    /// Maximum payload bytes charged while building the output.
    pub max_payload_bytes: usize,
    /// Maximum serialisation work.
    pub max_work: usize,
}
impl Default for WriteOptions {
    fn default() -> Self {
        Self {
            max_xml_bytes: 64 * 1024 * 1024,
            max_records: 1_000_000,
            max_payload_bytes: 256 * 1024 * 1024,
            max_work: 50_000_000,
        }
    }
}

/// A loaded map and the warnings its read produced: reference
/// inconsistencies the source retains, and the quantities
/// [`ReadOptions::discard_mismatched_quantities`] discarded.
#[derive(Clone, Debug, PartialEq)]
pub struct ReadReport {
    /// The map that was read.
    pub map: ConsensusMap,
    /// One line per retained inconsistency or discarded quantity.
    pub warnings: Vec<String>,
}
fn bad(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}
fn one<'a>(node: &'a Node, name: &str) -> Result<&'a Node> {
    let mut values = node.children.iter().filter(|n| n.name == name);
    let value = values
        .next()
        .ok_or_else(|| bad(format!("{} requires {name}", node.name)))?;
    if values.next().is_some() {
        return Err(bad(format!("duplicate {name}")));
    }
    Ok(value)
}
fn f32_value(value: &str) -> Result<f32> {
    let result = xml::finite(value)? as f32;
    if result.is_finite() {
        Ok(result)
    } else {
        Err(bad("consensus value exceeds f32"))
    }
}

/// Read consensusXML with the default options and the global modification
/// registry. See [`read_report`].
pub fn read(reader: impl BufRead) -> Result<ConsensusMap> {
    read_with_options(reader, &ReadOptions::default())
}
/// Read with explicit options and the global modification registry. See
/// [`read_report`].
pub fn read_with_options(reader: impl BufRead, options: &ReadOptions) -> Result<ConsensusMap> {
    read_with_registry(reader, options, ModificationsDB::global())
}
/// Read with caller-supplied modification definitions, dropping the
/// warnings. See [`read_report`].
pub fn read_with_registry(
    reader: impl BufRead,
    options: &ReadOptions,
    registry: &ModificationsDB,
) -> Result<ConsensusMap> {
    Ok(read_report(reader, options, registry)?.map)
}
/// Read into `map`, replacing it only on success.
pub fn read_into(
    reader: impl BufRead,
    map: &mut ConsensusMap,
    options: &ReadOptions,
) -> Result<()> {
    let replacement = read_with_options(reader, options)?;
    *map = replacement;
    Ok(())
}

/// Read a bounded consensusXML document and return the map with its warnings.
///
/// # Errors
///
/// Malformed or unrepresentable input, a quantity ownership mismatch unless
/// [`ReadOptions::discard_mismatched_quantities`] is set, and any exceeded
/// limit in `options`.
pub fn read_report(
    reader: impl BufRead,
    options: &ReadOptions,
    registry: &ModificationsDB,
) -> Result<ReadReport> {
    read_report_reporting(reader, options, registry, &mut ProgressReporter::silent())
}

/// The elements at whose start the source handler calls
/// `setProgress(++progress_)` (`ConsensusXMLHandler.cpp:149`, `:173`, `:334`,
/// `:424`, `:485`, `:582`), besides the root's own call.
const LOAD_PROGRESS_ELEMENTS: [&str; 6] = [
    "map",
    "consensusElement",
    "IdentificationRun",
    "ProteinHit",
    "PeptideHit",
    "dataProcessing",
];

/// How many elements of `node`'s subtree, `node` excluded, are in
/// [`LOAD_PROGRESS_ELEMENTS`].
fn load_progress_calls(node: &Node) -> u64 {
    let mut calls = 0;
    let mut pending: Vec<&Node> = node.children.iter().collect();
    while let Some(child) = pending.pop() {
        if LOAD_PROGRESS_ELEMENTS.contains(&child.name.as_str()) {
            calls += 1;
        }
        pending.extend(child.children.iter());
    }
    calls
}

/// [`read_report`], with the source handler's loading section. The document
/// is parsed whole before any of it is converted, so the section's calls are
/// made once the parse succeeded, before the conversion, and its end after.
fn read_report_reporting(
    reader: impl BufRead,
    options: &ReadOptions,
    registry: &ModificationsDB,
    progress: &mut ProgressReporter<'_>,
) -> Result<ReadReport> {
    options.validate()?;
    let mut work = options.max_work;
    let mut budget = options.max_payload_bytes;
    let root =
        xml::parse_xml_with_budget(reader, &options.xml(), 32, None, &mut work, &mut budget)?;
    xml::measure_node(&root, &mut work, &mut budget)?;
    if root.name != "consensusXML" {
        return Err(bad("expected consensusXML root"));
    }
    if progress.is_reporting() {
        // `ConsensusXMLHandler.cpp:254-256`: a zero-width section, then one
        // call for the root and one per counted element, in a running count.
        progress.start(0, 0, "loading consensusXML file")?;
        let calls = 1 + load_progress_calls(&root);
        for value in 1..=calls {
            progress.set(
                i64::try_from(value).map_err(|_| bad("consensusXML progress count overflows"))?,
            )?;
        }
    }
    root.check(
        &[
            "version",
            "document_id",
            "id",
            "unique_id",
            "experiment_type",
            "xmlns:xsi",
            "xsi:noNamespaceSchemaLocation",
        ],
        &[
            "UserParam",
            "dataProcessing",
            "IdentificationRun",
            "UnassignedPeptideIdentification",
            "mapList",
            "consensusElementList",
        ],
    )?;
    let mut warnings = Vec::new();
    if xml::finite(root.optional("version").unwrap_or("1.0"))? > 1.7 {
        warnings.push("consensusXML version is newer than 1.7".into());
    }
    let mut map = ConsensusMap {
        identifier: root.optional("document_id").unwrap_or("").into(),
        unique_id: common::unique_id(
            root.optional("unique_id")
                .or_else(|| root.optional("id"))
                .unwrap_or(""),
        )?,
        metadata: xml::read_meta(&root, &options.xml())?,
        ..Default::default()
    };
    if let Some(kind) = root.optional("experiment_type") {
        map.experiment_type = kind.into();
    }
    let mut registry = xml::clone_registry(registry, &mut work, &mut budget)?;
    let mut context = common::ReadContext::default();
    for child in &root.children {
        match child.name.as_str() {
            "dataProcessing" => map
                .data_processing
                .push(common::read_processing(child, &options.xml())?),
            "IdentificationRun" => {
                let mut run = common::read_run(
                    child,
                    &options.xml(),
                    &mut context,
                    &mut registry,
                    &mut work,
                    &mut budget,
                )?;
                let before = warnings.len();
                run.protein_groups = groups::read(
                    &mut run.metadata,
                    "protein_group",
                    &context.protein_accessions,
                    &mut warnings,
                    options.max_list_items,
                    &mut work,
                    &mut budget,
                )?;
                run.indistinguishable_groups = groups::read(
                    &mut run.metadata,
                    "indistinguishable_proteins",
                    &context.protein_accessions,
                    &mut warnings,
                    options.max_list_items,
                    &mut work,
                    &mut budget,
                )?;
                if warnings.len() > before && !options.discard_mismatched_quantities {
                    return Err(bad("protein-group quantity ownership mismatch"));
                }
                map.protein_identifications.push(run);
            }
            "UnassignedPeptideIdentification" => {
                map.unassigned_peptide_identifications
                    .push(common::read_peptide(
                        child,
                        &options.xml(),
                        &context,
                        &registry,
                        &mut work,
                        &mut budget,
                    )?)
            }
            _ => {}
        }
    }
    let list = one(&root, "mapList")?;
    list.check(&["count"], &["map"])?;
    let declared = xml::number::<usize>(list.get("count")?)?;
    if declared != list.children.len() {
        warnings.push("mapList count differs from its entries".into());
    }
    for child in &list.children {
        child.check(
            &["id", "name", "unique_id", "label", "size"],
            &["UserParam"],
        )?;
        let index = xml::number::<u64>(child.get("id")?)?;
        let value = ColumnHeader {
            filename: child.get("name")?.into(),
            label: child.optional("label").unwrap_or("").into(),
            size: xml::number(child.optional("size").unwrap_or("0"))?,
            unique_id: common::unique_id(child.optional("unique_id").unwrap_or(""))?,
            metadata: xml::read_meta(child, &options.xml())?,
        };
        if map.column_headers.insert(index, value).is_some() {
            return Err(bad("duplicate consensus column ID"));
        }
    }
    let list = one(&root, "consensusElementList")?;
    list.check(&[], &["consensusElement"])?;
    for element in &list.children {
        let feature = read_element(
            element,
            options,
            &context,
            &registry,
            &mut work,
            &mut budget,
        )?;
        if options.accepts(&feature) {
            map.features.push(feature);
        }
    }
    map.validate()?;
    if let Err(error) = map.validate_consistency() {
        warnings.push(format!("inconsistent consensus map: {error}"));
    }
    // `ConsensusXMLHandler.cpp:130-133`, at `</consensusXML>`.
    progress.end()?;
    Ok(ReadReport { map, warnings })
}
fn read_element(
    node: &Node,
    options: &ReadOptions,
    context: &common::ReadContext,
    registry: &ModificationsDB,
    work: &mut usize,
    bytes: &mut usize,
) -> Result<ConsensusFeature> {
    node.check(
        &["id", "quality", "charge"],
        &[
            "centroid",
            "groupedElementList",
            "PeptideIdentification",
            "UserParam",
        ],
    )?;
    let mut feature = ConsensusFeature::default();
    feature.unique_id = common::unique_id(node.get("id")?)?;
    feature.quality = f32_value(node.optional("quality").unwrap_or("0"))?;
    feature.charge = xml::number(node.optional("charge").unwrap_or("0"))?;
    feature.metadata = xml::read_meta(node, &options.xml())?;
    let centroid = one(node, "centroid")?;
    centroid.check(&["rt", "mz", "it"], &[])?;
    let position = (
        xml::finite(centroid.get("rt")?)?,
        xml::finite(centroid.get("mz")?)?,
        f32_value(centroid.get("it")?)?,
    );
    let list = one(node, "groupedElementList")?;
    list.check(&[], &["element"])?;
    let mut handles = Vec::new();
    for child in &list.children {
        child.check(&["map", "id", "rt", "mz", "it", "charge"], &[])?;
        if !child.get("map")?.is_empty() && !child.get("id")?.is_empty() {
            handles.push(FeatureHandle {
                map_index: common::unique_id(child.get("map")?)?,
                unique_id: common::unique_id(child.get("id")?)?,
                rt: xml::finite(child.get("rt")?)?,
                mz: xml::finite(child.get("mz")?)?,
                intensity: f32_value(child.get("it")?)?,
                charge: xml::number(child.optional("charge").unwrap_or("0"))?,
                width: 0.0,
            });
        }
        // Source commits centroid on each element, even if its empty ID skipped the handle.
        feature.rt = position.0;
        feature.mz = position.1;
        feature.intensity = position.2;
    }
    feature.set_handles(handles)?;
    for peptide in node
        .children
        .iter()
        .filter(|n| n.name == "PeptideIdentification")
    {
        feature.peptide_identifications.push(common::read_peptide(
            peptide,
            &options.xml(),
            context,
            registry,
            work,
            bytes,
        )?);
    }
    feature.validate()?;
    Ok(feature)
}

/// Load a plain, gzip or bzip2 consensusXML file with the default options,
/// recording its path and type on the returned map.
pub fn load(path: impl AsRef<Path>) -> Result<ConsensusMap> {
    load_with_options(path, &ReadOptions::default())
}
/// Load with explicit options. See [`load`].
pub fn load_with_options(path: impl AsRef<Path>, options: &ReadOptions) -> Result<ConsensusMap> {
    load_reporting(path, options, &mut ProgressReporter::silent())
}
/// Load a consensusXML file, reporting progress as source
/// `ConsensusXMLFile::load` does.
///
/// The source file hands its handler only its log type
/// (`ConsensusXMLFile.cpp:90-92`), so the calls go to a fresh backend of
/// `logger`'s type, as [`ProgressLogger::clone`] makes one: the command
/// backend for [`Cmd`](crate::concept::progress_logger::ProgressLogType::Cmd),
/// the logger's GUI factory for
/// [`Gui`](crate::concept::progress_logger::ProgressLogType::Gui), and none
/// for the default type. A backend installed on `logger` with
/// [`ProgressLogger::set_logger`] receives nothing, as a source file's
/// `setLogger` backend receives nothing.
///
/// The calls are the handler's: `startProgress(0, 0, "loading consensusXML
/// file")` at the root, then `setProgress(1)`, `setProgress(2)`, … — one for
/// the root and one for every `map`, `consensusElement`, `IdentificationRun`,
/// `ProteinHit`, `PeptideHit` and `dataProcessing` element — and
/// `endProgress()` at `</consensusXML>`. With a zero-width range the command
/// backend prints one dot per call. The result is the one
/// [`load_with_options`] returns, and so is every error: both run the same
/// code, whose calls go nowhere for [`load_with_options`].
///
/// This reader parses the whole document before converting any of it, where
/// the source converts as it parses. A document that is not well-formed
/// therefore makes no call, where the source has made the calls for the
/// elements before the defect; one the conversion refuses has made every
/// set, and no end. A failure after the start leaves the section open, as in
/// the source, where the exception bypasses `endProgress`.
///
/// # Errors
///
/// As [`load`], plus the errors of the progress calls
/// ([`ProgressLogger::start_progress`] and its siblings).
pub fn load_with_progress(
    path: impl AsRef<Path>,
    options: &ReadOptions,
    logger: &ProgressLogger,
) -> Result<ConsensusMap> {
    let mut handler = logger.clone();
    load_reporting(
        path,
        options,
        &mut ProgressReporter::new(Some(&mut handler)),
    )
}
fn load_reporting(
    path: impl AsRef<Path>,
    options: &ReadOptions,
    progress: &mut ProgressReporter<'_>,
) -> Result<ConsensusMap> {
    let path = path.as_ref();
    let mut map = read_report_reporting(
        super::path_io::open(path)?,
        options,
        ModificationsDB::global(),
        progress,
    )?
    .map;
    map.loaded_file_path = path
        .to_str()
        .ok_or_else(|| bad("loaded filename must be UTF-8"))?
        .into();
    map.loaded_file_type = FileType::ConsensusXml;
    Ok(map)
}
/// Load into `map`, replacing it only on success.
pub fn load_into(
    path: impl AsRef<Path>,
    map: &mut ConsensusMap,
    options: &ReadOptions,
) -> Result<()> {
    let replacement = load_with_options(path, options)?;
    *map = replacement;
    Ok(())
}
/// Validate a consensusXML file against the bundled `ConsensusXML_1_7.xsd`.
///
/// Source `ConsensusXMLFile::isValid(filename, os)`, inherited from
/// `Internal::XMLFile`: the messages the source writes to `os` are the
/// report's diagnostics, and the source's `bool` is
/// [`is_valid`](crate::format::xml_schema::SchemaValidationReport::is_valid).
/// Available with the `xml-schema` feature, which brings in the libxml2
/// validator; the source always has Xerces.
///
/// # Errors
///
/// As [`xml_schema::validate`](crate::format::xml_schema::validate): an I/O
/// failure, where the source throws `Exception::FileNotFound`, and input that
/// is not well-formed XML, where the source returns `false`.
#[cfg(feature = "xml-schema")]
pub fn is_valid(
    path: impl AsRef<Path>,
) -> Result<crate::format::xml_schema::SchemaValidationReport> {
    crate::format::xml_schema::validate(crate::format::xml_schema::SchemaKind::ConsensusXML, path)
}

/// Write consensusXML 1.7 with the default ceilings and the global
/// modification registry.
pub fn write(writer: impl Write, map: &ConsensusMap) -> Result<()> {
    write_with_options(writer, map, &WriteOptions::default())
}
/// Write with explicit ceilings and the global modification registry.
pub fn write_with_options(
    writer: impl Write,
    map: &ConsensusMap,
    options: &WriteOptions,
) -> Result<()> {
    write_with_registry(writer, map, options, ModificationsDB::global())
}
/// Build and check the complete document against caller-supplied
/// modification definitions, then write it; an I/O failure can still leave
/// partial bytes in `writer`.
pub fn write_with_registry(
    mut writer: impl Write,
    map: &ConsensusMap,
    options: &WriteOptions,
    registry: &ModificationsDB,
) -> Result<()> {
    let bytes = encode(map, options, registry, &mut ProgressReporter::silent())?;
    writer.write_all(&bytes)?;
    writer.flush()?;
    Ok(())
}
/// The checked document, rendered, with the source handler's storing section
/// (`ConsensusXMLHandler.cpp:606-837`) around the rendering.
///
/// The source counts `setProgress(++progress_)` three times before the
/// header, once after the map's user parameters, once after the data
/// processing, and once per identification run, column header and consensus
/// feature. The port builds and checks the whole document first, so a map it
/// refuses makes no call.
fn encode(
    map: &ConsensusMap,
    options: &WriteOptions,
    registry: &ModificationsDB,
    progress: &mut ProgressReporter<'_>,
) -> Result<Vec<u8>> {
    let root = write_node(map, options, registry)?;
    if progress.is_reporting() {
        progress.start(0, 0, "storing consensusXML file")?;
        let calls = 5usize
            .saturating_add(map.protein_identifications.len())
            .saturating_add(map.column_headers.len())
            .saturating_add(map.features.len());
        for value in 1..=calls {
            progress.set_count(value)?;
        }
    }
    let bytes = xml::render(
        &root,
        &xml::WriteOptions {
            max_xml_bytes: options.max_xml_bytes,
            max_records: options.max_records,
        },
    )?;
    progress.end()?;
    Ok(bytes)
}
fn peptides(map: &ConsensusMap) -> impl Iterator<Item = &PeptideIdentification> {
    map.unassigned_peptide_identifications.iter().chain(
        map.features
            .iter()
            .flat_map(|f| f.peptide_identifications.iter()),
    )
}
fn write_node(
    map: &ConsensusMap,
    options: &WriteOptions,
    registry: &ModificationsDB,
) -> Result<Node> {
    if options.max_xml_bytes == 0
        || options.max_records == 0
        || options.max_payload_bytes == 0
        || options.max_work == 0
    {
        return Err(bad("consensusXML limits must be positive"));
    }
    // Preflight graph/sequence payload before references, metadata, or XML are cloned.
    let mut work = options.max_work;
    let mut bytes = options.max_payload_bytes;
    measure_map(map, &mut work, &mut bytes)?;
    xml::measure_identifications(
        &map.protein_identifications,
        &map.unassigned_peptide_identifications,
        &mut work,
        &mut bytes,
    )?;
    for feature in &map.features {
        xml::measure_identifications(&[], &feature.peptide_identifications, &mut work, &mut bytes)?;
    }
    map.validate()?;
    let context = common::WriteContext::new(&map.protein_identifications)?;
    let definitions = super::modification_definitions::collect_iter_with_budget(
        &map.protein_identifications,
        peptides(map),
        registry,
        &mut work,
        &mut bytes,
    )?;
    let mut registry = xml::clone_registry(registry, &mut work, &mut bytes)?;
    let mut root = Node::new("consensusXML");
    root.attr("version", "1.7");
    root.attr("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    root.attr("xsi:noNamespaceSchemaLocation","https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/ConsensusXML_1_7.xsd");
    if !map.identifier.is_empty() {
        root.attr("document_id", &map.identifier);
    }
    if map.unique_id != 0 {
        root.attr("id", format!("cm_{}", map.unique_id));
    }
    root.attr("experiment_type", &map.experiment_type);
    xml::write_meta(&mut root, &map.metadata)?;
    for processing in &map.data_processing {
        root.children.push(common::write_processing(processing)?);
    }
    for run in &map.protein_identifications {
        xml::measure_search(&run.search_parameters, &mut work, &mut bytes)?;
        let mut search = run.search_parameters.clone();
        if let Some(records) = definitions.get(&run.identifier) {
            super::modification_definitions::attach_with_budget(
                &mut search,
                records,
                &mut work,
                &mut bytes,
            )?;
        }
        super::modification_definitions::register_search_parameters_with_budget(
            &search,
            &mut registry,
            &mut work,
            &mut bytes,
        )?;
        let mut meta = run.metadata.clone();
        let refs = &context.protein_refs[&run.identifier];
        groups::write(&mut meta, "protein_group", &run.protein_groups, refs)?;
        groups::write(
            &mut meta,
            "indistinguishable_proteins",
            &run.indistinguishable_groups,
            refs,
        )?;
        root.children
            .push(common::write_run(run, &context, &search, &meta)?);
    }
    for peptide in &map.unassigned_peptide_identifications {
        root.children.push(common::write_peptide(
            peptide,
            "UnassignedPeptideIdentification",
            &context,
            &registry,
            &mut work,
            &mut bytes,
        )?);
    }
    let mut list = Node::new("mapList");
    list.attr("count", map.column_headers.len());
    for (id, header) in &map.column_headers {
        let mut node = Node::new("map");
        node.attr("id", id);
        node.attr("name", &header.filename);
        node.attr("label", &header.label);
        node.attr("size", header.size);
        if header.unique_id != 0 {
            node.attr("unique_id", header.unique_id);
        }
        xml::write_meta(&mut node, &header.metadata)?;
        list.children.push(node);
    }
    root.children.push(list);
    let mut list = Node::new("consensusElementList");
    for feature in &map.features {
        if feature.width != 0.0 || feature.handles().iter().any(|h| h.width != 0.0) {
            return Err(Error::Unsupported(
                "consensusXML does not transport feature/handle width".into(),
            ));
        }
        if feature.handles().is_empty()
            && (feature.rt != 0.0 || feature.mz != 0.0 || feature.intensity != 0.0)
        {
            return Err(Error::Unsupported(
                "source consensusXML cannot restore a nonzero centroid without handles".into(),
            ));
        }
        let mut node = Node::new("consensusElement");
        node.attr("id", format!("e_{}", feature.unique_id));
        node.attr("quality", feature.quality);
        if feature.charge != 0 {
            node.attr("charge", feature.charge);
        }
        let mut centroid = Node::new("centroid");
        centroid.attr("rt", feature.rt);
        centroid.attr("mz", feature.mz);
        centroid.attr("it", feature.intensity);
        node.children.push(centroid);
        let mut handles = Node::new("groupedElementList");
        for handle in feature.handles() {
            let mut child = Node::new("element");
            child.attr("map", handle.map_index);
            child.attr("id", handle.unique_id);
            child.attr("rt", handle.rt);
            child.attr("mz", handle.mz);
            child.attr("it", handle.intensity);
            if handle.charge != 0 {
                child.attr("charge", handle.charge);
            }
            handles.children.push(child);
        }
        node.children.push(handles);
        for peptide in &feature.peptide_identifications {
            node.children.push(common::write_peptide(
                peptide,
                "PeptideIdentification",
                &context,
                &registry,
                &mut work,
                &mut bytes,
            )?);
        }
        xml::write_meta(&mut node, &feature.metadata)?;
        list.children.push(node);
    }
    root.children.push(list);
    Ok(root)
}
fn charge(remaining: &mut usize, amount: usize) -> Result<()> {
    *remaining = remaining
        .checked_sub(amount)
        .ok_or_else(|| bad("consensusXML payload/work limit exceeded"))?;
    Ok(())
}
fn measure_map(map: &ConsensusMap, work: &mut usize, bytes: &mut usize) -> Result<()> {
    charge(
        work,
        map.features
            .len()
            .saturating_add(map.column_headers.len())
            .saturating_add(map.data_processing.len()),
    )?;
    charge(
        bytes,
        map.identifier
            .len()
            .saturating_add(map.experiment_type.len()),
    )?;
    xml::measure_metadata(&map.metadata, work, bytes)?;
    for header in map.column_headers.values() {
        charge(work, 1)?;
        charge(
            bytes,
            header
                .filename
                .len()
                .saturating_add(header.label.len())
                .saturating_add(128),
        )?;
        xml::measure_metadata(&header.metadata, work, bytes)?;
    }
    for processing in &map.data_processing {
        if !processing.software.cv_terms.is_empty()
            || !processing.software.cv_terms.metadata.is_empty()
        {
            return Err(bad(
                "consensusXML cannot represent software CV terms or metadata",
            ));
        }
        charge(work, processing.actions.len().saturating_add(1))?;
        charge(
            bytes,
            processing
                .software
                .name
                .len()
                .saturating_add(processing.software.version.len())
                .saturating_add(128),
        )?;
        xml::measure_metadata(&processing.metadata, work, bytes)?;
    }
    for feature in &map.features {
        charge(work, feature.handles().len().saturating_add(1))?;
        charge(
            bytes,
            feature
                .handles()
                .len()
                .saturating_mul(4096)
                .saturating_add(4096),
        )?;
        xml::measure_metadata(&feature.metadata, work, bytes)?;
    }
    Ok(())
}
/// Store `map` at `path` with the default ceilings, replacing the destination
/// atomically; `.gz` and `.bz2` suffixes select output compression.
pub fn store(path: impl AsRef<Path>, map: &ConsensusMap) -> Result<()> {
    store_with_options(path, map, &WriteOptions::default())
}
/// Store with explicit ceilings. See [`store`].
pub fn store_with_options(
    path: impl AsRef<Path>,
    map: &ConsensusMap,
    options: &WriteOptions,
) -> Result<()> {
    let path = path.as_ref();
    let kind = super::file_types::type_by_file_name(
        path.to_str()
            .ok_or_else(|| bad("output filename must be UTF-8"))?,
    );
    if kind != FileType::Unknown && kind != FileType::ConsensusXml {
        return Err(bad("output extension is not consensusXML"));
    }
    let mut bytes = Vec::new();
    write_with_options(&mut bytes, map, options)?;
    super::path_io::store(path, &bytes)
}
/// Store `map` at `path`, reporting progress as source
/// `ConsensusXMLFile::store` does.
///
/// As for [`load_with_progress`], the calls go to a fresh backend of
/// `logger`'s type (`ConsensusXMLFile.cpp:76-78`). They are the handler's
/// `writeTo` calls: `startProgress(0, 0, "storing consensusXML file")`, then
/// `setProgress(1)` … `setProgress(5 + identification runs + column headers +
/// consensus features)`, and `endProgress()`. As in the source, the
/// destination is opened first (`XMLFile::save_`), so one that cannot be
/// created makes no call. The bytes and every error are those of
/// [`store_with_options`], which builds the same document without reporting;
/// the port's checks precede the start.
///
/// # Errors
///
/// As [`store`], plus the errors of the progress calls.
pub fn store_with_progress(
    path: impl AsRef<Path>,
    map: &ConsensusMap,
    options: &WriteOptions,
    logger: &ProgressLogger,
) -> Result<()> {
    let path = path.as_ref();
    let kind = super::file_types::type_by_file_name(
        path.to_str()
            .ok_or_else(|| bad("output filename must be UTF-8"))?,
    );
    if kind != FileType::Unknown && kind != FileType::ConsensusXml {
        return Err(bad("output extension is not consensusXML"));
    }
    let mut handler = logger.clone();
    super::path_io::store_reporting(
        path,
        |progress| {
            Ok((
                encode(map, options, ModificationsDB::global(), progress)?,
                (),
            ))
        },
        &mut ProgressReporter::new(Some(&mut handler)),
    )
}
