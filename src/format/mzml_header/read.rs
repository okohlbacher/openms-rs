// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Header registry construction and reference resolution for the mzML reader.
//!
//! Source `MzMLHandler::startElement` collects `sourceFile`, `sample`,
//! `software`, `instrumentConfiguration` and `dataProcessing` definitions into
//! `std::map` members keyed by ID (`MzMLHandler.cpp:1248-1290`) and resolves
//! references against them. This module builds the equivalent registry once
//! the `run` element opens, and resolves record, array and scan references
//! from it afterwards. See `docs/MZML_HEADER_SUPPORT.md`.

use super::*;

/// The error message and warning label for an unresolved `softwareRef`.
const SOFTWARE_REF: &str = "softwareRef";
/// The error message and warning label for an unresolved processing reference:
/// `dataProcessingRef` on a record or array, or a list's
/// `defaultDataProcessingRef`.
const DATA_PROCESSING_REF: &str = "dataProcessingRef";

/// Resolution policy for header references that name no preceding definition,
/// with the per-read warning state of the source-compatible policy.
///
/// Source `MzMLHandler` looks `softwareRef` and every data-processing reference
/// up with `std::map::operator[]` (`MzMLHandler.cpp:920`, `:924`, `:948`,
/// `:952`, `:1034`, `:1264` and `:1288`). A missing key default-constructs the
/// value, so a dangling reference silently becomes an empty `Software` or an
/// empty processing history. The native default refuses that loss; the source
/// policy is selected with `mzml::ReadOptions::source_dangling_references`.
#[derive(Default)]
pub(crate) struct DanglingReferences {
    /// Substitute the source's default-constructed value instead of failing.
    source: bool,
    /// Dangling `softwareRef` IDs already reported during this read.
    software: BTreeSet<String>,
    /// Dangling processing-reference IDs already reported during this read.
    processing: BTreeSet<String>,
}
impl DanglingReferences {
    /// A policy that rejects dangling references, or substitutes the source's
    /// empty values when `source` is `true`.
    pub fn new(source: bool) -> Self {
        Self {
            source,
            ..Self::default()
        }
    }
    /// Handle the reference `id` of kind `label` (`softwareRef` or
    /// `dataProcessingRef`) that names no definition.
    ///
    /// # Errors
    ///
    /// Under the default policy returns [`Error::Parse`] with the message
    /// `unresolved <label>`, unchanged from the strict reader. Under the
    /// source policy returns an error only when the allowance cannot cover
    /// the lookup or the record of a first occurrence.
    ///
    /// # Notes
    ///
    /// Under the source policy, the first occurrence of each distinct ID per
    /// kind writes one line to the crate's warning log stream
    /// (`LogLevel::Warn`, standard error unless reconfigured), so a file whose
    /// every record names the same dangling ID warns once. The source writes
    /// nothing here; the warning is native, so a caller that opted into the
    /// loss can still see it. A logging failure never changes the read result.
    fn dangling(&mut self, label: &'static str, id: &str, work: &mut Work) -> Result<()> {
        if !self.source {
            return Err(invalid(format!("unresolved {label}")));
        }
        let (warned, replacement) = if label == SOFTWARE_REF {
            (&mut self.software, "empty software")
        } else {
            (&mut self.processing, "an empty processing history")
        };
        work.charge(id.len().saturating_mul(64), 0)?;
        if warned.contains(id) {
            return Ok(());
        }
        work.meter().tree::<String>(1)?;
        let copy = work.copy(id)?;
        // `id` passed `parameter_id`, so it cannot carry a line break.
        let _ = crate::concept::log_stream::log_message(
            crate::concept::log_stream::LogLevel::Warn,
            format_args!(
                "Warning: mzML {label} '{id}' names no definition; source-compatible reading uses {replacement}."
            ),
        );
        warned.insert(copy);
        Ok(())
    }
}

/// Resolved header definitions, consulted by record, array, scan and precursor
/// references once the `run` element has opened.
#[derive(Default)]
pub(crate) struct Registry {
    source_files: BTreeMap<String, SourceFile>,
    processing: BTreeMap<String, Vec<Arc<DataProcessing>>>,
    instrument_ids: BTreeSet<String>,
    default_instrument: String,
    dangling: DanglingReferences,
}
impl Registry {
    /// A deep copy of the source file with ID `id`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Parse`] for a malformed or unresolved `sourceFileRef`,
    /// under either dangling-reference policy, and when the allowance is
    /// exhausted. Source `MzMLHandler.cpp:899-906` warns about an unregistered
    /// spectrum source file and continues; that leniency is not ported.
    pub fn source(&self, id: &str, work: &mut Work) -> Result<SourceFile> {
        let id = parameter_id(id)?;
        work.charge(id.len().saturating_mul(64), 0)?;
        let value = self
            .source_files
            .get(id)
            .ok_or_else(|| invalid("unresolved sourceFileRef"))?;
        copy_source(value, work)
    }
    /// The processing history with ID `id`, sharing the definition's `Arc`
    /// handles.
    ///
    /// Serves `dataProcessingRef` on spectra, chromatograms and binary arrays,
    /// and the lists' `defaultDataProcessingRef`. An ID that names no
    /// definition yields an empty history under the source policy, as source
    /// `processing_[ref]` does, and is warned about once per read.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Parse`] for a malformed ID, for an unresolved ID under
    /// the default policy (`unresolved dataProcessingRef`), and when the
    /// allowance is exhausted.
    pub fn processing(&mut self, id: &str, work: &mut Work) -> Result<Vec<Arc<DataProcessing>>> {
        let id = parameter_id(id)?;
        work.charge(id.len().saturating_mul(64), 0)?;
        let Some(value) = self.processing.get(id) else {
            self.dangling.dangling(DATA_PROCESSING_REF, id, work)?;
            return Ok(Vec::new());
        };
        work.slots::<Arc<DataProcessing>>(value.len())?;
        Ok(value.clone())
    }
    /// The `source_file_name` and `source_file_path` metadata of an optional
    /// `sourceFileRef` attribute on a scan or precursor.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Parse`] for a malformed or unresolved reference under
    /// either dangling-reference policy, and when the parameter budget or the
    /// allowance is exhausted.
    pub fn source_metadata(
        &self,
        attrs: &BTreeMap<String, String>,
        budget: &mut ParameterBudget,
        work: &mut Work,
    ) -> Result<MetaInfo> {
        let mut metadata = MetaInfo::new();
        if let Some(id) = attrs.get("sourceFileRef") {
            work.charge(id.len().saturating_mul(64), 0)?;
            let value = self
                .source_files
                .get(parameter_id(id)?)
                .ok_or_else(|| invalid("unresolved sourceFileRef"))?;
            budget.attribute("source_file_name".len(), value.name.len())?;
            budget.attribute("source_file_path".len(), value.path.len())?;
            work.meter().tree::<(String, MetaValue)>(2)?;
            work.charge(32, 32)?;
            metadata.insert("source_file_name".into(), work.copy(&value.name)?.into());
            metadata.insert("source_file_path".into(), work.copy(&value.path)?.into());
        }
        Ok(metadata)
    }
    /// A scan's acquisition: source-file metadata, `externalSpectrumID`, and an
    /// `instrumentConfigurationRef` recorded as metadata when it names a
    /// configuration other than the run default.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Parse`] for a malformed or unresolved source-file or
    /// instrument reference under either dangling-reference policy, and when
    /// the parameter budget or the allowance is exhausted.
    pub fn scan(
        &self,
        attrs: &BTreeMap<String, String>,
        budget: &mut ParameterBudget,
        work: &mut Work,
    ) -> Result<Acquisition> {
        let mut scan = Acquisition {
            metadata: self.source_metadata(attrs, budget, work)?,
            ..Default::default()
        };
        if let Some(id) = attrs.get("externalSpectrumID") {
            scan.identifier = work.copy(id)?;
        }
        if let Some(id) = attrs.get("instrumentConfigurationRef") {
            let id = parameter_id(id)?;
            work.charge(id.len().saturating_mul(64), 0)?;
            if !self.instrument_ids.contains(id) {
                return Err(invalid("unresolved scan instrumentConfigurationRef"));
            }
            if id != self.default_instrument {
                work.meter().tree::<(String, MetaValue)>(1)?;
                scan.metadata
                    .insert("instrument_configuration_ref".into(), work.copy(id)?.into());
            }
        }
        Ok(scan)
    }
}

/// Build the registry from the captured header lists and apply the `run`
/// element's `sampleRef`, `defaultInstrumentConfigurationRef` and
/// `startTimeStamp` to `settings`.
///
/// `source_dangling_references` selects the source treatment of a
/// `softwareRef` or data-processing reference that names no preceding
/// definition (see `DanglingReferences`); the registry keeps that policy for
/// the record references resolved after `run` opens.
pub(super) fn parse(
    roots: Vec<Node>,
    attrs: &BTreeMap<String, String>,
    groups: &BTreeMap<String, Vec<Parameter>>,
    parameters: &mut ParameterBudget,
    work: &mut Work,
    settings: &mut ExperimentalSettings,
    source_dangling_references: bool,
) -> Result<Registry> {
    // Ordinary readers start empty. A retaining transform starts from exact
    // cloned lengths, so its first append can move twice the old descriptor
    // count. Per-new-value charges cover subsequent geometric growth.
    work.slots::<ContactPerson>(
        settings
            .contacts
            .len()
            .checked_mul(2)
            .ok_or_else(resource)?,
    )?;
    work.slots::<SourceFile>(
        settings
            .source_files
            .len()
            .checked_mul(2)
            .ok_or_else(resource)?,
    )?;
    let mut cx = Context {
        groups,
        parameters,
        work,
    };
    let mut result = Registry {
        dangling: DanglingReferences::new(source_dangling_references),
        ..Registry::default()
    };
    let mut samples = BTreeMap::new();
    let mut software = BTreeMap::new();
    let mut instruments = BTreeMap::new();
    let mut seen = BTreeSet::new();
    let mut all_ids = BTreeSet::new();
    // Source ownership is resolved in document order. Forward parameter groups
    // are available here, but software/instrument registry references must name
    // a preceding definition as required by the mzML header layout.
    for root in &roots {
        cx.work.meter().tree::<&str>(1)?;
        if !seen.insert(root.name.as_str()) {
            return Err(invalid("duplicate mzML header list"));
        }
        match root.name.as_str() {
            "fileDescription" => {
                let mut source_list = false;
                for node in &root.children {
                    match node.name.as_str() {
                        "fileContent" => {
                            cx.params(node, |_, _, _| Ok(()))?;
                        }
                        "sourceFileList" => {
                            if source_list {
                                return Err(invalid("duplicate sourceFileList"));
                            }
                            source_list = true;
                            for node in node.children("sourceFile")? {
                                let id = register_id(node, &mut all_ids, cx.work)?;
                                cx.work.meter().tree::<(String, SourceFile)>(1)?;
                                let value = read_source(node, &mut cx)?;
                                result.source_files.insert(id, value);
                            }
                        }
                        "contact" => {
                            cx.work.slots::<ContactPerson>(4)?;
                            settings.contacts.push(read_contact(node, &mut cx)?);
                        }
                        _ => return Err(invalid("invalid fileDescription child")),
                    }
                }
            }
            "sampleList" => {
                for node in root.children("sample")? {
                    let id = register_id(node, &mut all_ids, cx.work)?;
                    cx.work.meter().tree::<(String, Sample)>(1)?;
                    samples.insert(id, read_sample(node, &mut cx)?);
                }
            }
            "softwareList" => {
                for node in root.children("software")? {
                    let id = register_id(node, &mut all_ids, cx.work)?;
                    cx.work.meter().tree::<(String, (Software, bool))>(1)?;
                    let value = read_software(node, &mut cx)?;
                    let marker = value.name == write::FALLBACK_SOFTWARE
                        && exact_placeholder_software(node, &mut cx)?;
                    software.insert(id, (value, marker));
                }
            }
            "instrumentConfigurationList" => {
                for node in root.children("instrumentConfiguration")? {
                    let id = register_id(node, &mut all_ids, cx.work)?;
                    cx.work.meter().tree::<String>(1)?;
                    cx.work.meter().tree::<(String, Instrument)>(1)?;
                    result.instrument_ids.insert(cx.work.copy(&id)?);
                    let value = read_instrument(node, &software, &mut result.dangling, &mut cx)?;
                    instruments.insert(id, value);
                }
            }
            "dataProcessingList" => {
                for node in root.children("dataProcessing")? {
                    let id = register_id(node, &mut all_ids, cx.work)?;
                    cx.work
                        .meter()
                        .tree::<(String, Vec<Arc<DataProcessing>>)>(1)?;
                    let mut methods = Vec::new();
                    if node.children.is_empty() {
                        return Err(invalid("dataProcessing has no methods"));
                    }
                    cx.work.slots::<Arc<DataProcessing>>(node.children.len())?;
                    methods
                        .try_reserve_exact(node.children.len())
                        .map_err(|_| resource())?;
                    for method in &node.children {
                        if method.name != "processingMethod" {
                            return Err(invalid("invalid dataProcessing child"));
                        }
                        // Source ignores order and retains XML encounter sequence.
                        let _: u64 = number(method.get("order")?, "processing method order")?;
                        let mut value = DataProcessing {
                            software: software_ref(
                                method.get("softwareRef")?,
                                &software,
                                &mut result.dangling,
                                cx.work,
                            )?,
                            ..Default::default()
                        };
                        cx.work.slots::<(DataProcessing, [usize; 2])>(1)?;
                        let mut count = 0usize;
                        let mut action_count = 0usize;
                        let mut clean = method
                            .attrs
                            .keys()
                            .all(|key| matches!(key.as_str(), "order" | "softwareRef"))
                            && method.children.iter().all(|child| {
                                matches!(child.name.as_str(), "cvParam" | "userParam")
                                    && child.children.is_empty()
                                    || child.name == "referenceableParamGroupRef"
                                        && child.children.is_empty()
                                        && child.attrs.keys().all(|key| key == "ref")
                            });
                        cx.params(method, |cx, kind, attrs| {
                            count += 1;
                            clean &= attrs.keys().all(|key| match kind {
                                "cvParam" => matches!(
                                    key.as_str(),
                                    "name"
                                        | "accession"
                                        | "cvRef"
                                        | "value"
                                        | "unitName"
                                        | "unitCvRef"
                                        | "unitAccession"
                                ),
                                "userParam" => matches!(
                                    key.as_str(),
                                    "name"
                                        | "type"
                                        | "value"
                                        | "unitName"
                                        | "unitCvRef"
                                        | "unitAccession"
                                ),
                                _ => false,
                            });
                            let recognized = processing_param(&mut value, cx, kind, attrs)?;
                            clean &= recognized;
                            if kind == "cvParam" {
                                let id = required(attrs, "accession")?;
                                if id == "MS:1000747" {
                                    clean &= !attrs.contains_key("unitAccession")
                                        && !attrs.contains_key("unitName")
                                        && !attrs.contains_key("unitCvRef");
                                }
                                if !matches!(
                                    id,
                                    "MS:1000747"
                                        | "MS:1000629"
                                        | "MS:1000631"
                                        | "MS:1000787"
                                        | "MS:1000788"
                                ) {
                                    action_count += 1;
                                    clean &= attrs.get("value").is_none_or(String::is_empty)
                                        && !attrs.contains_key("unitAccession")
                                        && !attrs.contains_key("unitName")
                                        && !attrs.contains_key("unitCvRef");
                                }
                            }
                            Ok(())
                        })?;
                        let history_marker = value.metadata.contains_key(write::EMPTY_HISTORY);
                        // The placeholder flag is looked up by the normalized
                        // ID `software_ref` resolved. Indexing by the raw
                        // attribute panicked on a whitespace-padded reference,
                        // and would on a dangling one under the source policy;
                        // neither names the exact placeholder software.
                        let placeholder_software = software
                            .get(parameter_id(method.get("softwareRef")?)?)
                            .is_some_and(|entry| entry.1);
                        if (history_marker || value.metadata.contains_key(write::EMPTY_ACTIONS))
                            && (!clean
                                || action_count != 1
                                || (history_marker && (count != 2 || !placeholder_software)))
                        {
                            return Err(invalid(
                                "processing marker contains additional or discarded XML payload",
                            ));
                        }
                        methods.push(Arc::new(value));
                    }
                    write::normalize_methods(&mut methods)?;
                    result.processing.insert(id, methods);
                }
            }
            _ => return Err(invalid("unknown retained mzML header")),
        }
    }
    if let Some(id) = attrs.get("sampleRef") {
        settings.sample = samples
            .remove(parameter_id(id)?)
            .ok_or_else(|| invalid("unresolved sampleRef"))?;
    }
    if let Some(id) = attrs.get("defaultInstrumentConfigurationRef") {
        let id = parameter_id(id)?;
        settings.instrument = instruments
            .remove(id)
            .ok_or_else(|| invalid("unresolved defaultInstrumentConfigurationRef"))?;
        result.default_instrument = cx.work.copy(id)?;
    }
    settings.instrument_configurations = instruments;
    if let Some(text) = attrs.get("startTimeStamp") {
        cx.work.charge(text.len().saturating_mul(4), text.len())?;
        settings.date_time = DateTime::parse(text)?;
        if text.len() > 19 {
            cx.meta(
                &mut settings.metadata,
                "mzml_start_time_stamp",
                text.as_str().into(),
            )?;
        }
    }
    // Source appends deduplicated values in lexical source ID order. Hashing is
    // only an accelerator; equality still decides, including metadata/units.
    let mut fingerprints: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
    cx.work.charge(settings.source_files.len(), 0)?;
    for (index, value) in settings.source_files.iter().enumerate() {
        crate::kernel::acquisition_fields::source(&mut cx.work.meter(), value)?;
        let hash = source_hash(value);
        cx.work.meter().tree::<(u64, Vec<usize>)>(1)?;
        cx.work.slots::<usize>(4)?;
        fingerprints.entry(hash).or_default().push(index);
    }
    for value in result.source_files.values() {
        crate::kernel::acquisition_fields::source(&mut cx.work.meter(), value)?;
        let hash = source_hash(value);
        cx.work.meter().tree::<(u64, Vec<usize>)>(1)?;
        let bucket = fingerprints.entry(hash).or_default();
        let mut duplicate = false;
        for &index in bucket.iter() {
            crate::kernel::acquisition_fields::source(&mut cx.work.meter(), value)?;
            if settings.source_files[index] == *value {
                duplicate = true;
                break;
            }
        }
        if !duplicate {
            cx.work.slots::<SourceFile>(4)?;
            cx.work.slots::<usize>(4)?;
            bucket.push(settings.source_files.len());
            settings.source_files.push(copy_source(value, cx.work)?);
        }
    }
    Ok(result)
}
fn register_id(node: &Node, ids: &mut BTreeSet<String>, work: &mut Work) -> Result<String> {
    let id = node.id()?;
    work.meter().tree::<String>(1)?;
    work.charge(id.len().saturating_mul(64), id.len().saturating_mul(2))?;
    if !ids.insert(id.into()) {
        return Err(invalid("duplicate mzML header ID"));
    }
    Ok(id.into())
}
/// Resolve a `softwareRef` against the software defined so far.
///
/// A dangling ID yields `Software::default()` under the source policy, as
/// source `software_[ref]` does for an instrument (`MzMLHandler.cpp:1288`) and
/// a processing method (`:1264`); the default policy fails with
/// `unresolved softwareRef`.
fn software_ref(
    id: &str,
    software: &BTreeMap<String, (Software, bool)>,
    dangling: &mut DanglingReferences,
    work: &mut Work,
) -> Result<Software> {
    let id = parameter_id(id)?;
    work.charge(id.len().saturating_mul(64), 0)?;
    let Some(value) = software.get(id) else {
        dangling.dangling(SOFTWARE_REF, id, work)?;
        work.slots::<Software>(1)?;
        return Ok(Software::default());
    };
    let value = &value.0;
    work.slots::<Software>(1)?;
    work.meter().text(&value.name)?;
    work.meter().text(&value.version)?;
    work.meter().cv(&value.cv_terms)?;
    Ok(value.clone())
}
fn copy_source(value: &SourceFile, work: &mut Work) -> Result<SourceFile> {
    crate::kernel::acquisition_fields::source(&mut work.meter(), value)?;
    Ok(value.clone())
}
fn source_hash(value: &SourceFile) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    value.name.hash(&mut h);
    value.path.hash(&mut h);
    value.file_type.hash(&mut h);
    value.checksum.hash(&mut h);
    (value.checksum_type as u8).hash(&mut h);
    value.native_id_type.hash(&mut h);
    value.native_id_type_accession.hash(&mut h);
    value.cv_terms.hash(&mut h);
    h.finish()
}
fn read_source(node: &Node, cx: &mut Context<'_>) -> Result<SourceFile> {
    let original_name = node.get("name")?;
    let original_path = node.get("location")?;
    let length = original_name
        .len()
        .checked_add(original_path.len())
        .ok_or_else(resource)?;
    cx.work
        .charge(length.saturating_mul(8), length.saturating_mul(6))?;
    let (name, path) = if original_path.is_empty() && !original_name.is_empty() {
        let directory = crate::system::file::path(original_name);
        (
            crate::system::file::basename(original_name),
            if directory == "." {
                "file://./"
            } else {
                directory
            },
        )
    } else {
        (original_name, original_path)
    };
    let mut path = path.to_owned();
    if path.starts_with("File://") {
        path = path.replace("File://", "file://");
    }
    if path.starts_with("FILE://") {
        path = path.replace("FILE://", "file://");
    }
    if path.starts_with("file:///.") {
        path = path.replace("file:///.", "file://./");
    }
    if path == "file:///" {
        path = "file://".into();
    }
    let mut value = SourceFile {
        name: name.into(),
        path,
        ..Default::default()
    };
    cx.params(node, |cx, kind, attrs| {
        if kind == "userParam" {
            let v = cx.value(kind, attrs)?;
            return cx.meta(&mut value.cv_terms.metadata, required(attrs, "name")?, v);
        }
        let id = required(attrs, "accession")?;
        let raw = attrs.get("value").map_or("", String::as_str);
        match id {
            "MS:1000569" | "MS:1000568" => {
                value.checksum = cx.work.copy(raw)?;
                value.checksum_type = if id == "MS:1000569" {
                    ChecksumType::Sha1
                } else {
                    ChecksumType::Md5
                };
            }
            _ if cx.work.child(id, "MS:1000560")? => {
                let term = cx.work.cv(id)?;
                value.file_type = cx.work.copy(&term.name)?;
            }
            _ if cx.work.child(id, "MS:1000767")? => {
                let term = cx.work.cv(id)?;
                value.native_id_type = cx.work.copy(&term.name)?;
                value.native_id_type_accession = cx.work.copy(&term.id)?;
            }
            _ => {}
        }
        Ok(())
    })?;
    Ok(value)
}
fn read_contact(node: &Node, cx: &mut Context<'_>) -> Result<ContactPerson> {
    let mut value = ContactPerson::default();
    cx.params(node, |cx, kind, attrs| {
        let raw = attrs.get("value").map_or("", String::as_str);
        if kind == "userParam" {
            let name = required(attrs, "name")?;
            if name == "contact_info" {
                if !value.contact_info.is_empty() {
                    return Err(invalid("duplicate contact_info"));
                }
                value.contact_info = cx.work.copy(raw)?;
            } else {
                let v = cx.value(kind, attrs)?;
                cx.meta(&mut value.metadata, name, v)?;
            }
            return Ok(());
        }
        match required(attrs, "accession")? {
            "MS:1000586" => {
                cx.work
                    .charge(raw.len().saturating_mul(4), raw.len().saturating_mul(2))?;
                value.set_name(raw)?;
            }
            "MS:1000587" => value.address = cx.work.copy(raw)?,
            "MS:1000588" => value.url = cx.work.copy(raw)?,
            "MS:1000589" => value.email = cx.work.copy(raw)?,
            "MS:1000590" => value.institution = cx.work.copy(raw)?,
            _ => {}
        }
        Ok(())
    })?;
    Ok(value)
}
fn read_sample(node: &Node, cx: &mut Context<'_>) -> Result<Sample> {
    let mut value = Sample {
        name: cx.work.copy(node.optional("name"))?,
        ..Default::default()
    };
    cx.params(node, |cx, kind, attrs| {
        let raw = attrs.get("value").map_or("", String::as_str);
        if kind == "userParam" {
            let name = required(attrs, "name")?;
            if name == "comment" {
                value.comment = cx.work.copy(raw)?;
            } else {
                let v = cx.value(kind, attrs)?;
                cx.meta(&mut value.metadata, name, v)?;
            }
            return Ok(());
        }
        let id = required(attrs, "accession")?;
        match id {
            "MS:1000001" => value.number = cx.work.copy(raw)?,
            "MS:1000004" => value.mass = sample_number(raw)?,
            "MS:1000005" => value.volume = sample_number(raw)?,
            "MS:1000006" => value.concentration = sample_number(raw)?,
            "MS:1000047" => value.state = SampleState::Emulsion,
            "MS:1000048" => value.state = SampleState::Gas,
            "MS:1000049" => value.state = SampleState::Liquid,
            "MS:1000050" => value.state = SampleState::Solid,
            "MS:1000051" => value.state = SampleState::Solution,
            "MS:1000052" => value.state = SampleState::Suspension,
            "MS:1000053" => {
                let v = cx.value(kind, attrs)?;
                cx.meta(&mut value.metadata, "sample batch", v)?;
            }
            _ if id.starts_with("PATO:") => {
                let v = cx.value(kind, attrs)?;
                cx.meta(&mut value.metadata, required(attrs, "name")?, v)?;
            }
            _ if id.starts_with("GO:") || id.starts_with("BTO:") => {
                let key = if id.starts_with("GO:") {
                    "GO cellular component"
                } else {
                    "brenda source tissue"
                };
                let v = cx.work.copy(required(attrs, "name")?)?;
                cx.meta(&mut value.metadata, key, v.into())?;
            }
            _ => {}
        }
        Ok(())
    })?;
    Ok(value)
}
fn read_software(node: &Node, cx: &mut Context<'_>) -> Result<Software> {
    let mut value = Software {
        version: cx.work.copy(node.get("version")?)?,
        ..Default::default()
    };
    cx.params(node, |cx, kind, attrs| {
        if kind == "userParam" {
            let v = cx.value(kind, attrs)?;
            return cx.meta(&mut value.cv_terms.metadata, required(attrs, "name")?, v);
        }
        let id = required(attrs, "accession")?;
        if cx.work.child(id, "MS:1000531")? {
            let name = if id == "MS:1000799" {
                attrs.get("value").map_or("", String::as_str)
            } else {
                required(attrs, "name")?
            };
            value.name = cx.work.copy(name)?;
        }
        Ok(())
    })?;
    Ok(value)
}
fn exact_placeholder_software(node: &Node, cx: &mut Context<'_>) -> Result<bool> {
    let mut count = 0;
    let mut exact = node
        .attrs
        .keys()
        .all(|key| matches!(key.as_str(), "id" | "version"))
        && node.children.iter().all(|child| {
            matches!(child.name.as_str(), "cvParam" | "userParam") && child.children.is_empty()
                || child.name == "referenceableParamGroupRef"
                    && child.children.is_empty()
                    && child.attrs.keys().all(|key| key == "ref")
        });
    cx.params(node, |_, kind, attrs| {
        count += 1;
        exact &= kind == "cvParam"
            && attrs.get("accession").is_some_and(|s| s == "MS:1000799")
            && attrs
                .get("value")
                .is_some_and(|s| s == write::FALLBACK_SOFTWARE)
            && attrs
                .keys()
                .all(|key| matches!(key.as_str(), "accession" | "name" | "cvRef" | "value"));
        Ok(())
    })?;
    Ok(exact && count == 1)
}
fn read_instrument(
    node: &Node,
    software: &BTreeMap<String, (Software, bool)>,
    dangling: &mut DanglingReferences,
    cx: &mut Context<'_>,
) -> Result<Instrument> {
    let mut value = Instrument::default();
    cx.params(node, |cx, kind, attrs| {
        if kind == "userParam" {
            let v = cx.value(kind, attrs)?;
            return cx.meta(&mut value.metadata, required(attrs, "name")?, v);
        }
        let id = required(attrs, "accession")?;
        let v = cx.value(kind, attrs)?;
        if !instrument::apply(
            instrument::Target::Instrument(&mut value),
            id,
            attrs.get("value").map_or("", String::as_str),
            v,
        )? && cx.work.child(id, "MS:1000031")?
        {
            let term = cx.work.cv(id)?;
            value.name = cx.work.copy(&term.name)?;
        }
        Ok(())
    })?;
    let (mut components, mut software_seen) = (false, false);
    for child in &node.children {
        match child.name.as_str() {
            "componentList" => {
                if components {
                    return Err(invalid("duplicate instrument componentList"));
                }
                components = true;
                let count: usize = number(child.get("count")?, "instrument component count")?;
                if count != child.children.len() {
                    return Err(invalid("component count mismatch"));
                }
                cx.work
                    .slots::<MassAnalyzer>(child.children.len().saturating_mul(4))?;
                for component in &child.children {
                    let order = number(component.get("order")?, "component order")?;
                    let mut source = IonSource {
                        order,
                        ..Default::default()
                    };
                    let mut analyzer = MassAnalyzer {
                        order,
                        ..Default::default()
                    };
                    let mut detector = IonDetector {
                        order,
                        ..Default::default()
                    };
                    if !matches!(component.name.as_str(), "source" | "analyzer" | "detector") {
                        return Err(invalid("invalid instrument component"));
                    }
                    cx.params(component, |cx, kind, attrs| {
                        if kind == "userParam" {
                            let meta = match component.name.as_str() {
                                "source" => &mut source.metadata,
                                "analyzer" => &mut analyzer.metadata,
                                _ => &mut detector.metadata,
                            };
                            let v = cx.value(kind, attrs)?;
                            return cx.meta(meta, required(attrs, "name")?, v);
                        }
                        let id = required(attrs, "accession")?;
                        let raw = attrs.get("value").map_or("", String::as_str);
                        let v = cx.value(kind, attrs)?;
                        let target = match component.name.as_str() {
                            "source" => instrument::Target::Source(&mut source),
                            "analyzer" => instrument::Target::Analyzer(&mut analyzer),
                            _ => instrument::Target::Detector(&mut detector),
                        };
                        if !instrument::apply(target, id, raw, v)? {
                            if component.name == "source" && cx.work.child(id, "MS:1000008")? {
                                cx.meta(&mut source.metadata, "ionization accession", id.into())?;
                            } else if component.name == "analyzer"
                                && cx.work.child(id, "MS:1000443")?
                            {
                                cx.meta(
                                    &mut analyzer.metadata,
                                    "mass analyzer accession",
                                    id.into(),
                                )?;
                            }
                        }
                        Ok(())
                    })?;
                    match component.name.as_str() {
                        "source" => value.ion_sources.push(source),
                        "analyzer" => value.mass_analyzers.push(analyzer),
                        _ => value.ion_detectors.push(detector),
                    }
                }
            }
            "softwareRef" => {
                if software_seen {
                    return Err(invalid("duplicate instrument softwareRef"));
                }
                software_seen = true;
                value.software = software_ref(child.get("ref")?, software, dangling, cx.work)?;
            }
            "cvParam" | "userParam" | "referenceableParamGroupRef" => {}
            _ => return Err(invalid("invalid instrumentConfiguration child")),
        }
    }
    Ok(value)
}
fn processing_param(
    value: &mut DataProcessing,
    cx: &mut Context<'_>,
    kind: &str,
    attrs: &BTreeMap<String, String>,
) -> Result<bool> {
    if kind == "userParam" {
        let v = cx.value(kind, attrs)?;
        cx.meta(&mut value.metadata, required(attrs, "name")?, v)?;
        return Ok(true);
    }
    let id = required(attrs, "accession")?;
    let key = match id {
        "MS:1000629" => Some("low_intensity_threshold"),
        "MS:1000631" => Some("high_intensity_threshold"),
        "MS:1000787" => Some("inclusive_low_intensity_threshold"),
        "MS:1000788" => Some("inclusive_high_intensity_threshold"),
        _ => None,
    };
    if let Some(key) = key {
        let v = cx.value(kind, attrs)?;
        cx.meta(&mut value.metadata, key, v)?;
        return Ok(true);
    }
    if id == "MS:1000747" {
        let raw = attrs.get("value").map_or("", String::as_str);
        let dt = DateTime::parse(raw)?;
        if !dt.is_valid() {
            return Err(invalid("invalid processing completion time"));
        }
        value.completion_time = Some(dt);
        return Ok(true);
    }
    if let Some(action) = processing_action(id, cx.work)? {
        cx.work.meter().tree::<ProcessingAction>(1)?;
        value.actions.insert(action);
        return Ok(true);
    }
    Ok(false)
}
pub(super) fn processing_action(id: &str, work: &mut Work) -> Result<Option<ProcessingAction>> {
    for &(accession, action, descendants) in PROCESSING_ACTIONS {
        if id == accession || (descendants && work.child(id, accession)?) {
            return Ok(Some(action));
        }
    }
    Ok(None)
}
pub(super) const PROCESSING_ACTIONS: &[(&str, ProcessingAction, bool)] = &[
    ("MS:1000530", ProcessingAction::FormatConversion, false),
    ("MS:1000544", ProcessingAction::ConversionMzML, false),
    ("MS:1000545", ProcessingAction::ConversionMzXML, false),
    ("MS:1000546", ProcessingAction::ConversionMzData, false),
    ("MS:1000741", ProcessingAction::ConversionDta, false),
    ("MS:1000543", ProcessingAction::DataProcessing, false),
    ("MS:1000033", ProcessingAction::Deisotoping, false),
    ("MS:1000034", ProcessingAction::ChargeDeconvolution, false),
    ("MS:1000035", ProcessingAction::PeakPicking, true),
    ("MS:1000592", ProcessingAction::Smoothing, true),
    ("MS:1000778", ProcessingAction::ChargeCalculation, true),
    ("MS:1000780", ProcessingAction::PrecursorRecalculation, true),
    ("MS:1000593", ProcessingAction::BaselineReduction, false),
    (
        "MS:1000745",
        ProcessingAction::RetentionTimeAlignment,
        false,
    ),
    (
        "MS:1001484",
        ProcessingAction::IntensityNormalization,
        false,
    ),
    ("MS:1001485", ProcessingAction::MzCalibration, false),
    ("MS:1001486", ProcessingAction::DataFiltering, true),
];

fn sample_number(raw: &str) -> Result<f64> {
    let n = f64::from_list_item(raw)?;
    if !n.is_finite() {
        return Err(invalid("nonfinite sample scalar"));
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;
    /// Route this test thread's warnings to a discarding sink; the source
    /// policy warns by design and the output is asserted in
    /// `tests/mzml_header_leniency.rs`.
    fn discard_warnings() {
        use crate::concept::log_stream::{LogLevel, LogSink, with_thread_local_log};
        with_thread_local_log(LogLevel::Warn, |log| {
            log.remove_all_streams()?;
            log.insert(&LogSink::new(std::io::sink()))
        })
        .unwrap();
    }
    #[test]
    fn full_registry_value_roots_are_charged_before_payload_parsing() {
        fn size<T>() -> usize {
            512 + 11 * std::mem::size_of::<(T, [usize; 4])>()
        }
        for (list, child, root_bytes, extra_id) in [
            ("sampleList", "sample", size::<(String, Sample)>(), 0),
            (
                "softwareList",
                "software",
                size::<(String, (Software, bool))>(),
                0,
            ),
            (
                "instrumentConfigurationList",
                "instrumentConfiguration",
                size::<(String, Instrument)>(),
                size::<String>(),
            ),
            (
                "dataProcessingList",
                "dataProcessing",
                size::<(String, Vec<Arc<DataProcessing>>)>(),
                0,
            ),
        ] {
            let roots = vec![Node {
                name: list.into(),
                attrs: [("count".into(), "1".into())].into(),
                children: vec![Node {
                    name: child.into(),
                    attrs: [("id".into(), "id".into())].into(),
                    children: vec![],
                }],
            }];
            // Enough for the borrowed seen-root set, copied ID set and any
            // separate instrument-ID set, but one byte short of the full map.
            let mut work = Work {
                remaining: usize::MAX,
                bytes: size::<&str>() + size::<String>() + 4 + extra_id + root_bytes - 1,
            };
            let mut params = ParameterBudget {
                remaining: usize::MAX,
                bytes: usize::MAX,
            };
            let error = parse(
                roots,
                &BTreeMap::new(),
                &BTreeMap::new(),
                &mut params,
                &mut work,
                &mut ExperimentalSettings::default(),
                false,
            )
            .err()
            .unwrap();
            assert!(
                error.to_string().contains("resource limit"),
                "{list}: {error}"
            );
        }
    }
    #[test]
    fn seeded_append_charges_existing_vector_growth_before_header_processing() {
        let mut settings = ExperimentalSettings {
            contacts: vec![ContactPerson::default(); 5],
            source_files: vec![SourceFile::default(); 7],
            ..Default::default()
        };
        let before = settings.clone();
        let required = 10 * size_of::<ContactPerson>() + 14 * size_of::<SourceFile>();
        let mut work = Work {
            remaining: usize::MAX,
            bytes: required - 1,
        };
        let mut parameters = ParameterBudget {
            remaining: usize::MAX,
            bytes: usize::MAX,
        };
        assert!(
            parse(
                Vec::new(),
                &BTreeMap::new(),
                &BTreeMap::new(),
                &mut parameters,
                &mut work,
                &mut settings,
                false,
            )
            .is_err()
        );
        assert_eq!(settings, before);
    }
    #[test]
    fn repeated_source_copies_share_the_same_remaining_allowance() {
        let mut registry = Registry::default();
        registry.source_files.insert(
            "id".into(),
            SourceFile {
                name: "x".repeat(4096),
                ..Default::default()
            },
        );
        let mut unlimited = Work {
            remaining: usize::MAX,
            bytes: usize::MAX,
        };
        let first = registry.source("id", &mut unlimited).unwrap();
        let mut one = Work {
            remaining: usize::MAX - unlimited.remaining,
            bytes: usize::MAX - unlimited.bytes,
        };
        assert_eq!(registry.source("id", &mut one).unwrap(), first);
        assert!(registry.source("id", &mut one).is_err());
        assert_eq!(registry.source_files["id"].name.len(), 4096);
    }
    #[test]
    fn dangling_processing_references_follow_the_selected_policy() {
        discard_warnings();
        let unlimited = || Work {
            remaining: usize::MAX,
            bytes: usize::MAX,
        };
        let mut strict = Registry::default();
        let error = strict.processing("absent", &mut unlimited()).unwrap_err();
        assert_eq!(
            error.to_string(),
            invalid("unresolved dataProcessingRef").to_string()
        );
        let mut source = Registry {
            dangling: DanglingReferences::new(true),
            ..Registry::default()
        };
        let mut work = unlimited();
        assert!(source.processing("absent", &mut work).unwrap().is_empty());
        let first = (usize::MAX - work.remaining, usize::MAX - work.bytes);
        assert!(source.processing("absent", &mut work).unwrap().is_empty());
        let second = (
            usize::MAX - work.remaining - first.0,
            usize::MAX - work.bytes - first.1,
        );
        // A repeated ID pays only its lookup; the first also records the ID.
        assert!(second.0 < first.0 && second.1 < first.1);
        assert_eq!(source.dangling.processing.len(), 1);
        assert!(source.dangling.software.is_empty());
        // Malformed IDs stay errors under the source policy.
        assert!(source.processing("1bad", &mut unlimited()).is_err());
    }
    #[test]
    fn dangling_reference_record_is_charged_before_it_is_retained() {
        discard_warnings();
        let mut policy = DanglingReferences::new(true);
        let mut probe = Work {
            remaining: usize::MAX,
            bytes: usize::MAX,
        };
        policy
            .dangling(SOFTWARE_REF, "so_in_0", &mut probe)
            .unwrap();
        let needed = usize::MAX - probe.bytes;
        let mut policy = DanglingReferences::new(true);
        let mut short = Work {
            remaining: usize::MAX,
            bytes: needed - 1,
        };
        assert!(
            policy
                .dangling(SOFTWARE_REF, "so_in_0", &mut short)
                .is_err()
        );
        assert!(policy.software.is_empty());
        let mut strict = DanglingReferences::default();
        let error = strict
            .dangling(SOFTWARE_REF, "so_in_0", &mut probe)
            .unwrap_err();
        assert_eq!(
            error.to_string(),
            invalid("unresolved softwareRef").to_string()
        );
    }
    #[test]
    fn record_processing_references_charge_handles_and_preserve_shared_identity() {
        let value = Arc::new(DataProcessing {
            metadata: [("large".into(), "x".repeat(100_000).into())].into(),
            ..Default::default()
        });
        let mut registry = Registry::default();
        registry
            .processing
            .insert("id".into(), vec![Arc::clone(&value)]);
        let mut work = Work {
            remaining: 1_000,
            bytes: std::mem::size_of::<Arc<DataProcessing>>(),
        };
        let read = registry.processing("id", &mut work).unwrap();
        assert!(Arc::ptr_eq(&read[0], &value));
        assert_eq!(work.bytes, 0);
        assert!(registry.processing("id", &mut work).is_err());
    }
}
