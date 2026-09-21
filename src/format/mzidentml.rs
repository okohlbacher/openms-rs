// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded native mzIdentML interchange: `FORMAT/MzIdentMLFile.h` and
//! `FORMAT/HANDLERS/MzIdentMLHandler.h`.
//!
//! mzIdentML is the PSI standard for identification results. Its sequence
//! library is declared once (`DBSequence`, `Peptide`, `PeptideEvidence`) and
//! referenced by id from every `SpectrumIdentificationResult` /
//! `SpectrumIdentificationItem`, so reading is mostly reference resolution.
//! [`docs/MZIDENTML_SUPPORT.md`](https://github.com/okohlbacher/OpenMS4-R/blob/main/docs/MZIDENTML_SUPPORT.md)
//! lists every source member, the reference-resolution table and the
//! divergences; the short version is:
//!
//! * The source reads with a DOM handler and writes with a stream handler, so
//!   the two halves disagree in places. Where they do, this port follows the
//!   mzIdentML schema and says so at the item.
//! * Cross-linking results (`MS:1002494`) take a second read and write path,
//!   as they do in C++: one match is reported as up to four
//!   `SpectrumIdentificationItem`s - an alpha and a beta chain, each light and
//!   heavy - which become one
//!   [`PeptideIdentification`](crate::identification::PeptideIdentification)
//!   carrying the OpenPepXL user parameters. See
//!   [`read_with_registry`](crate::format::mzidentml::read_with_registry).
//! * Every id this port writes is positional, so `store` is reproducible; the
//!   source draws them from `UniqueIdGenerator` and stamps wall-clock times.
//!
//! The entry points mirror the sibling idXML adapter:
//! [`load`](crate::format::mzidentml::load) /
//! [`store`](crate::format::mzidentml::store) for paths,
//! [`read`](crate::format::mzidentml::read) /
//! [`write`](crate::format::mzidentml::write) for streams, and
//! [`detect_version`](crate::format::mzidentml::detect_version) for the
//! declared schema version.

use super::identification_xml::{Node, bad, unsupported};
use crate::chemistry::{
    AASequence, CrossLinksDB, ModificationsDB, ProteaseDB, ResidueModification, TermSpecificity,
};
use crate::comparison::Tolerance;
use crate::concept::constants::user_param;
use crate::format::controlled_vocabulary::{CVTermDefinition, ControlledVocabulary, XRefType};
use crate::identification::{
    EnzymeTermSpecificity, FlankingResidue, PeakAnnotation, PeptideEvidence, PeptideHit,
    PeptideIdentification, ProteinHit, ProteinIdentification, SearchParameters,
};
use crate::metadata::{MetaInfo, MetaValue, MetaValueData, Unit};
use crate::{Error, Result};
use quick_xml::events::Event;
use quick_xml::name::ResolveResult;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, Read, Write};
use std::path::Path;

/// The schema version this adapter writes, as in `MzIdentMLFile()`.
///
/// The source constructs `XMLFile("/SCHEMAS/mzIdentML1.3.0.xsd", "1.3.0")` and
/// its writer then hard-codes `1.3.0` again rather than using that member, so
/// the constructed version cannot select an older output schema. This port has
/// one constant instead of two independent copies.
pub const SCHEMA_VERSION: &str = "1.3.0";

/// Versions `detectVersion` recognises, newest first.
pub const KNOWN_VERSIONS: [&str; 4] = ["1.3.0", "1.2.0", "1.1.0", "1.0.0"];

/// Namespace prefix every mzIdentML revision shares.
pub const NAMESPACE_PREFIX: &str = "http://psidev.info/psi/pi/mzIdentML/";

/// Lines of the document header `detectVersion` inspects.
const VERSION_HEADER_LINES: usize = 15;

/// Accession of the "crosslinking search" term that marks an XL-MS document.
const CROSSLINKING_SEARCH: &str = "MS:1002494";

/// Accession of "crosslink donor", the alpha half of a cross-link.
const CROSSLINK_DONOR: &str = "MS:1002509";

/// Accession of "crosslink acceptor", the beta half of a cross-link.
const CROSSLINK_ACCEPTOR: &str = "MS:1002510";

/// Accession of "crosslink spectrum identification item", which groups the
/// `SpectrumIdentificationItem`s of one cross-link spectrum match.
const CROSSLINK_ITEM: &str = "MS:1002511";

/// The `extra_features` list `OPXLHelper::addPercolatorFeatureList` installs.
const PERCOLATOR_FEATURES: [&str; 32] = [
    "precursor_mz_error_ppm",
    "OpenPepXL:score",
    "isotope_error",
    "OpenPepXL:xquest_score",
    "OpenPepXL:xcorr xlink",
    "OpenPepXL:xcorr common",
    "OpenPepXL:match-odds",
    "OpenPepXL:intsum",
    "OpenPepXL:wTIC",
    "OpenPepXL:TIC",
    "OpenPepXL:prescore",
    "OpenPepXL:log_occupancy",
    "OpenPepXL:log_occupancy_alpha",
    "OpenPepXL:log_occupancy_beta",
    "matched_xlink_alpha",
    "matched_xlink_beta",
    "matched_linear_alpha",
    "matched_linear_beta",
    "ppm_error_abs_sum_linear_alpha",
    "ppm_error_abs_sum_linear_beta",
    "ppm_error_abs_sum_xlinks_alpha",
    "ppm_error_abs_sum_xlinks_beta",
    "ppm_error_abs_sum_linear",
    "ppm_error_abs_sum_xlinks",
    "ppm_error_abs_sum_alpha",
    "ppm_error_abs_sum_beta",
    "ppm_error_abs_sum",
    "precursor_total_intensity",
    "precursor_target_intensity",
    "precursor_signal_proportion",
    "precursor_target_peak_count",
    "precursor_residual_peak_count",
];

/// Ceilings for one read. Each is checked before anything is allocated, so a
/// refused document leaves the destination untouched.
///
/// `max_xml_bytes` bounds the encoded input, `max_payload_bytes` the decoded
/// tree plus the identifications built from it, and `max_work` the parser and
/// resolution steps. `max_elements` and `max_depth` bound the tree's shape;
/// `max_list_items` bounds one library, evidence or hit list.
#[derive(Clone, Copy, Debug)]
pub struct ReadOptions {
    /// Encoded input bytes accepted, including the byte order mark.
    pub max_xml_bytes: usize,
    /// Total XML elements, `cvParam` and `userParam` included.
    pub max_elements: usize,
    /// Maximum element nesting depth.
    pub max_depth: usize,
    /// Cumulative decoded tree and identification payload.
    pub max_payload_bytes: usize,
    /// Cumulative parser and reference-resolution steps.
    pub max_work: usize,
    /// Maximum entries in one library, evidence or hit list.
    pub max_list_items: usize,
}
impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            max_xml_bytes: 64 * 1024 * 1024,
            max_elements: 2_000_000,
            max_depth: 64,
            max_payload_bytes: 256 * 1024 * 1024,
            max_work: 50_000_000,
            max_list_items: 1_000_000,
        }
    }
}

/// Ceilings and reproducibility choices for one write.
///
/// `creation_date` fills the root `creationDate` attribute. The source writes
/// `DateTime::now()` unconditionally, which makes two stores of one document
/// differ; the default here is `None`, which omits the optional attribute.
/// Pass an `xs:dateTime` string to emit one.
#[derive(Clone, Debug)]
pub struct WriteOptions {
    /// Output bytes accepted before the write is refused.
    pub max_output_bytes: usize,
    /// Maximum elements emitted.
    pub max_records: usize,
    /// Optional `creationDate` for the root element.
    pub creation_date: Option<String>,
}
impl Default for WriteOptions {
    fn default() -> Self {
        Self {
            max_output_bytes: 64 * 1024 * 1024,
            max_records: 1_000_000,
            creation_date: None,
        }
    }
}

/// One mzIdentML document as the crate's identification types.
///
/// `document_id` is the root `id` attribute, preserved so a load/store cycle
/// does not invent a new document identity. Every other XML id is a transport
/// reference and is regenerated positionally on write, as in
/// [`IdXmlDocument`](crate::format::idxml::IdXmlDocument).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MzIdentMLDocument {
    /// Root `id` attribute; `store` substitutes `OpenMS` when empty.
    pub document_id: String,
    /// One entry per `SpectrumIdentification`, in document order.
    pub protein_identifications: Vec<ProteinIdentification>,
    /// One entry per `SpectrumIdentificationResult`, in document order.
    pub peptide_identifications: Vec<PeptideIdentification>,
}

/// Maximum protein and peptide identifications accepted by one write.
pub const MAX_ITEMS: usize = 1_000_000;

// ---------------------------------------------------------------------------
// MzIdentMLHandler.h's public analysisXML value types
// ---------------------------------------------------------------------------

/// One `SpectrumIdentificationItem` as the source's `Internal::IdentificationHit`.
///
/// The source declares this next to the stream handler and fills it in
/// `onStartElement`; nothing ever reads the filled value back out, because
/// `MzIdentMLFile::load` uses the DOM handler instead (see
/// `OpenMS_CPP_ISSUES.md`). It is ported as a value type because it is public
/// API of the header, and because it is the only source structure that keeps a
/// PSM's mzIdentML identity (`id`, `name`, `passThreshold`) next to its
/// numbers. Nothing in this module consumes it.
///
/// `rank` keeps the source's 0-based convention: `onStartElement` stores
/// `rank - 1` from the file's 1-based attribute.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IdentificationHit {
    /// `id` attribute of the `SpectrumIdentificationItem`.
    pub id: String,
    /// Charge state; the source default is 0.
    pub charge: i32,
    /// `calculatedMassToCharge`; the source default is 0.0.
    pub calculated_mass_to_charge: f64,
    /// `experimentalMassToCharge`; the source default is 0.0.
    pub experimental_mass_to_charge: f64,
    /// Optional human-readable `name` attribute.
    pub name: String,
    /// `passThreshold`; the source default is `true`.
    pub pass_threshold: bool,
    /// 0-based rank, converted from the file's 1-based `rank`.
    pub rank: i32,
    /// `MetaInfoInterface` payload of the source class.
    pub metadata: MetaInfo,
}
impl IdentificationHit {
    /// Source default: `pass_threshold` is `true`, every number 0.
    pub fn new() -> Self {
        Self {
            pass_threshold: true,
            ..Default::default()
        }
    }
}

/// One `SpectrumIdentificationResult` as the source's
/// `Internal::SpectrumIdentification`.
///
/// The source's `id_` member is protected and has no accessor; it is exposed
/// here because the type is otherwise write-only. As with
/// [`IdentificationHit`], nothing in this module consumes it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpectrumIdentification {
    /// Identifier; `protected` and inaccessible in the source.
    pub id: String,
    /// Candidate hits for this spectrum.
    pub hits: Vec<IdentificationHit>,
    /// `MetaInfoInterface` payload of the source class.
    pub metadata: MetaInfo,
}
impl SpectrumIdentification {
    /// Append one hit, as `addHit`.
    pub fn add_hit(&mut self, hit: IdentificationHit) {
        self.hits.push(hit);
    }
}

/// One analysisXML instance as the source's `Internal::Identification`.
///
/// `creation_date` is the source's `DateTime creation_date_`, kept as the
/// `xs:dateTime` text rather than a parsed instant so an unparsable date is
/// preserved instead of silently becoming the epoch. As with
/// [`IdentificationHit`], nothing in this module consumes it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Identification {
    /// Identifier; `protected` and inaccessible in the source.
    pub id: String,
    /// Date and time the search was performed, verbatim.
    pub creation_date: Option<String>,
    /// Spectrum identifications of this run.
    pub spectrum_identifications: Vec<SpectrumIdentification>,
    /// `MetaInfoInterface` payload of the source class.
    pub metadata: MetaInfo,
}
impl Identification {
    /// Append one spectrum identification, as `addSpectrumIdentification`.
    pub fn add_spectrum_identification(&mut self, value: SpectrumIdentification) {
        self.spectrum_identifications.push(value);
    }
}

// ---------------------------------------------------------------------------
// Budget
// ---------------------------------------------------------------------------

/// Remaining work and payload allowance for one operation.
struct Budget {
    work: usize,
    bytes: usize,
}
impl Budget {
    fn spend(&mut self, work: usize, bytes: usize) -> Result<()> {
        self.work = self
            .work
            .checked_sub(work)
            .ok_or_else(|| bad("mzIdentML work limit exceeded"))?;
        self.bytes = self
            .bytes
            .checked_sub(bytes)
            .ok_or_else(|| bad("mzIdentML payload limit exceeded"))?;
        Ok(())
    }
    fn text(&mut self, value: &str) -> Result<()> {
        self.spend(
            value.len().saturating_add(1),
            value.len().saturating_mul(4).saturating_add(64),
        )
    }
}

fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
fn missing(message: impl Into<String>) -> Error {
    Error::MissingInformation(message.into())
}

// ---------------------------------------------------------------------------
// Version detection
// ---------------------------------------------------------------------------

/// Detect the mzIdentML version declared in `path`, as `detectVersion`.
///
/// Reads the first 15 lines, trims each, joins them with a single space, and
/// looks for `version="x.y.z"` for a shipped schema version, newest first;
/// then for the `mzIdentML/x.y` target namespace, which carries only
/// major.minor. Falls back to [`SCHEMA_VERSION`] when neither matches, exactly
/// as the source falls back to the adapter default.
///
/// # Errors
///
/// [`Error::Io`] if the file cannot be read, and [`Error::Parse`] if the header
/// is not valid UTF-8. The source reads through `TextFile`, which accepts any
/// byte sequence and would compare raw bytes.
pub fn detect_version(path: impl AsRef<Path>) -> Result<String> {
    detect_version_from_reader(super::path_io::open(path.as_ref())?)
}

/// Detect the declared version from an open reader; see [`detect_version`].
///
/// # Errors
///
/// [`Error::Io`] on a read failure and [`Error::Parse`] when the header is not
/// valid UTF-8.
pub fn detect_version_from_reader(reader: impl BufRead) -> Result<String> {
    let mut header = String::new();
    for (index, line) in reader.lines().take(VERSION_HEADER_LINES).enumerate() {
        let line = line?;
        if index != 0 {
            header.push(' ');
        }
        header.push_str(line.trim());
    }
    Ok(detect_version_in_header(&header))
}

/// Validate an mzIdentML file against the bundled schema for the version it
/// declares, with the default native limits.
///
/// Source `MzIdentMLFile::isValid(filename, os, used_version)`:
/// [`detect_version`] picks 1.0.0, 1.1.0, 1.2.0 or 1.3.0 and the matching
/// `mzIdentML<version>.xsd` is used. The source's `used_version` out-parameter
/// is the report's `schema.version()`. Every version [`detect_version`] can
/// answer is bundled, so the source's fallback to the default 1.3.0 schema,
/// for a detected version with no shipped schema, cannot arise. The messages
/// the source writes to `os` are the report's diagnostics, and its `bool` is
/// [`is_valid`](crate::format::xml_schema::SchemaValidationReport::is_valid).
/// Available with the `xml-schema` feature, which brings in the libxml2
/// validator; the source always has Xerces.
///
/// # Errors
///
/// As [`is_valid_with_options`].
#[cfg(feature = "xml-schema")]
pub fn is_valid(
    path: impl AsRef<Path>,
) -> Result<crate::format::xml_schema::SchemaValidationReport> {
    is_valid_with_options(
        path,
        &crate::format::xml_schema::SchemaValidationOptions::default(),
    )
}

/// [`is_valid`] with explicit native limits. The version header is read
/// within `options.limits.max_xml_bytes` too, so a document with one
/// enormous first line cannot make detection read past the validation limit.
///
/// # Errors
///
/// As [`detect_version`], then as
/// [`xml_schema::validate_with_options`](crate::format::xml_schema::validate_with_options):
/// an I/O failure, where the source throws `Exception::FileNotFound`, and
/// input that is not well-formed XML, where the source returns `false`.
#[cfg(feature = "xml-schema")]
pub fn is_valid_with_options(
    path: impl AsRef<Path>,
    options: &crate::format::xml_schema::SchemaValidationOptions,
) -> Result<crate::format::xml_schema::SchemaValidationReport> {
    use crate::format::xml_schema::{SchemaKind, validate_with_options};
    let path = path.as_ref();
    let limit = u64::try_from(options.limits.max_xml_bytes).unwrap_or(u64::MAX);
    let version = detect_version_from_reader(super::path_io::open(path)?.take(limit))?;
    let schema = SchemaKind::mzidentml(&version)
        .ok_or_else(|| Error::InvalidValue(format!("no bundled mzIdentML {version} schema")))?;
    validate_with_options(schema, path, options)
}

fn detect_version_in_header(header: &str) -> String {
    for version in KNOWN_VERSIONS {
        if header.contains(&format!("version=\"{version}\"")) {
            return version.to_owned();
        }
    }
    for version in KNOWN_VERSIONS {
        // The namespace carries only major.minor, e.g. "1.1" for "1.1.0".
        if let Some((major_minor, _)) = version.rsplit_once('.') {
            if header.contains(&format!("mzIdentML/{major_minor}")) {
                return version.to_owned();
            }
        }
    }
    SCHEMA_VERSION.to_owned()
}
// ---------------------------------------------------------------------------
// Bounded namespace-aware parse into the shared identification node tree
// ---------------------------------------------------------------------------

fn decode(input: impl Read, limit: usize, budget: &mut Budget) -> Result<String> {
    let count = u64::try_from(limit)
        .ok()
        .and_then(|v| v.checked_add(1))
        .ok_or_else(|| bad("mzIdentML byte limit overflows"))?;
    let mut bytes = Vec::new();
    input.take(count).read_to_end(&mut bytes)?;
    if bytes.len() > limit {
        return Err(bad("mzIdentML byte limit exceeded"));
    }
    let bytes = bytes.strip_prefix(&[0xef, 0xbb, 0xbf]).unwrap_or(&bytes);
    budget.spend(bytes.len(), bytes.len().saturating_mul(3))?;
    let text = std::str::from_utf8(bytes).map_err(|_| {
        unsupported("mzIdentML input is not valid UTF-8; only UTF-8 documents are read")
    })?;
    // XML 1.0 line-ending normalization happens before attribute normalization.
    Ok(text.replace("\r\n", "\n").replace('\r', "\n"))
}

/// The mzIdentML namespace of the document, or `None` for a namespace-less one.
fn root_namespace(resolved: &ResolveResult<'_>) -> Result<Option<String>> {
    let namespace = namespace_of(resolved)?;
    match &namespace {
        None => Ok(namespace),
        Some(text) if text.starts_with(NAMESPACE_PREFIX) => Ok(namespace),
        Some(text) => Err(unsupported(format!(
            "root element namespace {text:?} is not an mzIdentML namespace"
        ))),
    }
}

fn namespace_of(resolved: &ResolveResult<'_>) -> Result<Option<String>> {
    match resolved {
        ResolveResult::Unbound => Ok(None),
        ResolveResult::Bound(ns) => Ok(Some(
            std::str::from_utf8(ns.as_ref())
                .map_err(|_| bad("XML namespace is not valid UTF-8"))?
                .to_owned(),
        )),
        ResolveResult::Unknown(prefix) => Err(bad(format!(
            "undeclared XML namespace prefix {:?}",
            String::from_utf8_lossy(prefix)
        ))),
    }
}

fn parse(input: impl BufRead, options: &ReadOptions, budget: &mut Budget) -> Result<Node> {
    if options.max_elements == 0
        || options.max_list_items == 0
        || options.max_depth == 0
        || options.max_depth > 512
    {
        return Err(invalid("invalid mzIdentML read limits"));
    }
    let text = decode(input, options.max_xml_bytes, budget)?;
    let mut reader = quick_xml::NsReader::from_reader(text.as_bytes());
    reader.config_mut().expand_empty_elements = true;
    reader.config_mut().enable_all_checks(true);
    let mut buffer = Vec::new();
    let mut stack: Vec<Node> = Vec::new();
    let mut root: Option<Node> = None;
    let mut namespace: Option<Option<String>> = None;
    let mut count = 0usize;
    loop {
        let (resolved, event) = reader
            .read_resolved_event_into(&mut buffer)
            .map_err(|e| bad(e.to_string()))?;
        budget.spend(1, 0)?;
        match event {
            Event::Start(start) => {
                if stack.is_empty() && root.is_some() {
                    return Err(bad("multiple XML roots"));
                }
                count += 1;
                if count > options.max_elements {
                    return Err(bad("mzIdentML element limit exceeded"));
                }
                if stack.len() >= options.max_depth {
                    return Err(bad("mzIdentML element depth limit exceeded"));
                }
                let name = std::str::from_utf8(start.local_name().as_ref())
                    .map_err(|e| bad(e.to_string()))?
                    .to_owned();
                match &namespace {
                    None => {
                        if name != "MzIdentML" {
                            return Err(bad(format!(
                                "mzIdentML root element must be MzIdentML, found {name}"
                            )));
                        }
                        namespace = Some(root_namespace(&resolved)?);
                    }
                    Some(expected) => {
                        if &namespace_of(&resolved)? != expected {
                            return Err(unsupported(format!(
                                "element {name} is not in the document's mzIdentML namespace"
                            )));
                        }
                    }
                }
                budget.text(&name)?;
                budget.spend(4, size_of::<Node>().saturating_mul(4))?;
                let mut node = Node::new(&name);
                for attr in start.attributes() {
                    let attr = attr.map_err(|e| bad(e.to_string()))?;
                    let key = std::str::from_utf8(attr.key.as_ref())
                        .map_err(|e| bad(e.to_string()))?
                        .to_owned();
                    if key == "xmlns" || key.starts_with("xmlns:") {
                        continue;
                    }
                    budget.spend(
                        key.len().saturating_add(attr.value.len().saturating_mul(4)),
                        key.len()
                            .saturating_add(attr.value.len().saturating_mul(6))
                            .saturating_add(128),
                    )?;
                    let raw = std::str::from_utf8(&attr.value).map_err(|e| bad(e.to_string()))?;
                    // XML attribute-value normalization, then entity expansion.
                    let normalized = raw.replace(['\n', '\t'], " ");
                    let value = quick_xml::escape::unescape(&normalized)
                        .map_err(|e| bad(e.to_string()))?
                        .into_owned();
                    if node.attrs.insert(key, value).is_some() {
                        return Err(bad("duplicate XML attribute"));
                    }
                }
                stack.push(node);
            }
            Event::End(_) => {
                let node = stack.pop().ok_or_else(|| bad("unexpected XML end tag"))?;
                match stack.last_mut() {
                    Some(parent) => parent.children.push(node),
                    None => root = Some(node),
                }
            }
            Event::Text(text) => {
                let decoded = text.xml_content().map_err(|e| bad(e.to_string()))?;
                budget.spend(decoded.len(), decoded.len().saturating_mul(4))?;
                match stack.last_mut() {
                    Some(node) => node.text.push_str(decoded.as_ref()),
                    None => {
                        if !decoded.trim_matches([' ', '\t', '\n']).is_empty() {
                            return Err(bad("text outside the XML root"));
                        }
                    }
                }
            }
            Event::GeneralRef(reference) => {
                // quick-xml reports every "&...;" in element text as its own
                // event, so an arm that ignores them deletes the character they
                // stand for: the sibling Mascot XML reader turned a score of
                // 1.5 into 15 that way. The five predefined entities and
                // numeric character references are resolved into the text here,
                // and anything else - an external entity, which needs a DTD
                // this reader refuses anyway - is an error.
                let name = reference.decode().map_err(|e| bad(e.to_string()))?;
                budget.spend(8, name.len())?;
                let resolved = match name.as_ref() {
                    "amp" => '&',
                    "lt" => '<',
                    "gt" => '>',
                    "apos" => '\'',
                    "quot" => '"',
                    value => {
                        let (digits, radix) = match value.strip_prefix("#x") {
                            Some(digits) => (digits, 16),
                            None => match value.strip_prefix('#') {
                                Some(digits) => (digits, 10),
                                None => {
                                    return Err(unsupported(format!(
                                        "mzIdentML must not reference the external XML entity {value:?}"
                                    )));
                                }
                            },
                        };
                        u32::from_str_radix(digits, radix)
                            .ok()
                            .and_then(char::from_u32)
                            .ok_or_else(|| bad("invalid XML character reference"))?
                    }
                };
                stack
                    .last_mut()
                    .ok_or_else(|| bad("XML entity reference outside the root"))?
                    .text
                    .push(resolved);
            }
            Event::CData(text) => {
                // CDATA is ordinary character data; upstream mzIdentML uses it
                // for the SiteRegexp of an Enzyme.
                let text = text.decode().map_err(|e| bad(e.to_string()))?;
                budget.spend(text.len(), text.len().saturating_mul(4))?;
                stack
                    .last_mut()
                    .ok_or_else(|| bad("XML CDATA outside the root"))?
                    .text
                    .push_str(&text);
            }
            Event::Eof => break,
            Event::Decl(declaration) => {
                if let Some(encoding) = declaration
                    .encoding()
                    .transpose()
                    .map_err(|error| bad(format!("mzIdentML declaration: {error}")))?
                {
                    if !encoding.eq_ignore_ascii_case(b"UTF-8") {
                        return Err(unsupported("mzIdentML supports only UTF-8 encoding"));
                    }
                }
            }
            Event::Comment(_) | Event::PI(_) => {}
            Event::DocType(_) => {
                return Err(unsupported(
                    "mzIdentML documents must not declare a DOCTYPE",
                ));
            }
            _ => {
                return Err(unsupported("unsupported XML content in mzIdentML"));
            }
        }
    }
    if !stack.is_empty() {
        return Err(bad("unclosed XML element"));
    }
    root.ok_or_else(|| bad("mzIdentML document has no root element"))
}

/// Every descendant with this tag name, in document order.
///
/// The source reaches its elements through `DOMDocument::getElementsByTagName`,
/// a document-wide scan that ignores nesting, so an element under the wrong
/// parent is still collected. This reproduces that deliberately: files that
/// place `SpectraData` or `Peptide` outside their schema parent still load,
/// exactly as they do in C++.
fn collect<'a>(node: &'a Node, name: &str, out: &mut Vec<&'a Node>) {
    if node.name == name {
        out.push(node);
    }
    for item in &node.children {
        collect(item, name, out);
    }
}
fn gather<'a>(root: &'a Node, name: &str) -> Vec<&'a Node> {
    let mut out = Vec::new();
    collect(root, name, &mut out);
    out
}
fn children<'a>(node: &'a Node, name: &'a str) -> impl Iterator<Item = &'a Node> {
    node.children.iter().filter(move |c| c.name == name)
}
fn child<'a>(node: &'a Node, name: &str) -> Option<&'a Node> {
    node.children.iter().find(|c| c.name == name)
}
fn attribute<'a>(node: &'a Node, name: &str) -> &'a str {
    node.optional(name).unwrap_or_default()
}

// ---------------------------------------------------------------------------
// cvParam / userParam
// ---------------------------------------------------------------------------

/// One `cvParam`, with the attributes the source's `parseCvParam_` keeps.
#[derive(Clone, Debug, Default)]
struct CvParam {
    accession: String,
    name: String,
    cv_ref: String,
    value: String,
    unit_accession: String,
    unit_name: String,
    unit_cv_ref: String,
}
impl CvParam {
    fn from_node(node: &Node) -> Self {
        Self {
            accession: attribute(node, "accession").to_owned(),
            name: attribute(node, "name").to_owned(),
            cv_ref: attribute(node, "cvRef").to_owned(),
            value: attribute(node, "value").to_owned(),
            unit_accession: attribute(node, "unitAccession").to_owned(),
            unit_name: attribute(node, "unitName").to_owned(),
            unit_cv_ref: attribute(node, "unitCvRef").to_owned(),
        }
    }
    /// The source builds `CVTerm::Unit` only when accession and name are both
    /// present, and warns about a missing `unitCvRef` instead of refusing.
    fn unit(&self) -> Result<Option<Unit>> {
        if self.unit_accession.is_empty() || self.unit_name.is_empty() {
            return Ok(None);
        }
        if !self.unit_accession.starts_with("UO:") && !self.unit_accession.starts_with("MS:") {
            // Source: "Unhandled unit" warning, and the unit is dropped.
            return Ok(None);
        }
        Ok(Some(Unit::new(
            &self.unit_accession,
            &self.unit_name,
            &self.unit_cv_ref,
        )?))
    }
}

/// `cvParam`s by accession and `userParam`s by name, as `parseParamGroup_`.
///
/// Both collections are keyed maps, so the source's iteration order is
/// lexicographic by accession and by name, not document order. That order is
/// observable: the PSM score type is the first matching accession in
/// lexicographic order, so the `BTreeMap` is part of the contract here.
#[derive(Clone, Debug, Default)]
struct ParamGroup {
    cv: BTreeMap<String, Vec<CvParam>>,
    user: BTreeMap<String, MetaValue>,
}

fn param_group(node: &Node, options: &ReadOptions, budget: &mut Budget) -> Result<ParamGroup> {
    let mut group = ParamGroup::default();
    let mut items = 0usize;
    for element in &node.children {
        match element.name.as_str() {
            "cvParam" => {
                items += 1;
                if items > options.max_list_items {
                    return Err(bad("mzIdentML parameter group limit exceeded"));
                }
                let param = CvParam::from_node(element);
                budget.text(&param.accession)?;
                budget.text(&param.value)?;
                group
                    .cv
                    .entry(param.accession.clone())
                    .or_default()
                    .push(param);
            }
            "userParam" => {
                items += 1;
                if items > options.max_list_items {
                    return Err(bad("mzIdentML parameter group limit exceeded"));
                }
                let name = attribute(element, "name").to_owned();
                budget.text(&name)?;
                let value = user_value(element)?;
                // std::map::insert keeps the first entry for a repeated name.
                group.user.entry(name).or_insert(value);
            }
            _ => {}
        }
    }
    Ok(group)
}

/// `userParam` value, typed by its `type` attribute as `XMLHandler::fromXSDString`.
///
/// An absent `value` attribute yields [`MetaValueData::Empty`], as the source's
/// `has_value` check does. Unsupported XSD types stay text.
fn user_value(node: &Node) -> Result<MetaValue> {
    let Some(text) = node.optional("value") else {
        return MetaValue::new(MetaValueData::Empty);
    };
    let data = match attribute(node, "type") {
        "xsd:double" | "xsd:float" | "xsd:decimal" => MetaValueData::Float(finite(text)?),
        "xsd:byte" | "xsd:int" | "xsd:unsignedShort" | "xsd:short" | "xsd:unsignedByte"
        | "xsd:unsignedInt" => MetaValueData::Integer(i64::from(
            text.trim()
                .parse::<i32>()
                .map_err(|_| bad(format!("invalid 32-bit userParam integer {text:?}")))?,
        )),
        "xsd:long"
        | "xsd:unsignedLong"
        | "xsd:integer"
        | "xsd:negativeInteger"
        | "xsd:nonNegativeInteger"
        | "xsd:nonPositiveInteger"
        | "xsd:positiveInteger" => MetaValueData::Integer(
            text.trim()
                .parse::<i64>()
                .map_err(|_| bad(format!("invalid 64-bit userParam integer {text:?}")))?,
        ),
        _ => MetaValueData::String(text.to_owned()),
    };
    let mut value = MetaValue::new(data)?;
    if let Some(unit) = CvParam::from_node(node).unit()? {
        value = value.with_unit(unit)?;
    }
    Ok(value)
}

fn finite(text: &str) -> Result<f64> {
    let value = text
        .trim()
        .parse::<f64>()
        .map_err(|_| bad(format!("invalid number {text:?}")))?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad(format!(
            "mzIdentML numbers must be finite, got {text:?}"
        )))
    }
}

fn integer(text: &str, label: &str) -> Result<i32> {
    text.trim()
        .parse::<i32>()
        .map_err(|_| bad(format!("invalid {label} {text:?}")))
}

/// Source `toDoubleOrZero_`: an absent score value is 0.0, never NaN, because
/// `PeptideHit` has no "has score" state and NaN would break every comparison.
fn score_value(text: &str) -> Result<f64> {
    if text.is_empty() {
        Ok(0.0)
    } else {
        finite(text)
    }
}

/// Source `toDoubleOrNaN_`: an absent position-like value is the NaN sentinel.
/// This port returns `None`, because the crate's coordinates are `Option<f64>`
/// and every stored float must be finite.
fn optional_value(text: &str) -> Result<Option<f64>> {
    if text.is_empty() {
        Ok(None)
    } else {
        Ok(Some(finite(text)?))
    }
}

/// `XMLHandler::cvParamToValue`: type the value through the CV, or yield
/// [`MetaValueData::Empty`] when the term is unknown, the value does not
/// convert, or a numeric term carries no value.
fn cv_param_value(cv: &ControlledVocabulary, param: &CvParam) -> Result<MetaValue> {
    let Ok(term) = cv.get_term(&param.accession) else {
        // Source warns "Unknown cvParam" and returns DataValue::EMPTY.
        return MetaValue::new(MetaValueData::Empty);
    };
    let empty = || MetaValue::new(MetaValueData::Empty);
    let data = if param.value.is_empty() {
        let numeric = !matches!(term.xref_type, XRefType::None | XRefType::String);
        if numeric
            && !cv
                .is_child_of(&param.accession, "MS:1000513")
                .unwrap_or(false)
        {
            return empty();
        }
        MetaValueData::String(String::new())
    } else {
        match term.xref_type {
            XRefType::Integer
            | XRefType::NegativeInteger
            | XRefType::PositiveInteger
            | XRefType::NonNegativeInteger
            | XRefType::NonPositiveInteger => match param.value.trim().parse::<i32>() {
                Ok(value) => MetaValueData::Integer(i64::from(value)),
                Err(_) => return empty(),
            },
            XRefType::Decimal => match param.value.trim().parse::<f64>() {
                Ok(value) if value.is_finite() => MetaValueData::Float(value),
                _ => return empty(),
            },
            XRefType::Boolean => {
                let lowered = param.value.to_ascii_lowercase();
                if lowered == "true" || lowered == "false" {
                    MetaValueData::String(lowered)
                } else {
                    return empty();
                }
            }
            _ => MetaValueData::String(param.value.clone()),
        }
    };
    let mut value = MetaValue::new(data)?;
    if let Some(unit) = param.unit()? {
        value = value.with_unit(unit)?;
    }
    Ok(value)
}
// ---------------------------------------------------------------------------
// Read
// ---------------------------------------------------------------------------

/// The CV score subtrees `initScoreTermCaches_` precomputes once per file.
///
/// `q` and `e` include their parents, because the source builds them with
/// `addAllChildTerms`; `specific` does too, which is why the score scan has to
/// exclude `MS:1001143` explicitly before testing it.
struct ScoreTerms {
    q: BTreeSet<String>,
    e: BTreeSet<String>,
    specific: BTreeSet<String>,
}
impl ScoreTerms {
    fn new(cv: &ControlledVocabulary) -> Result<Self> {
        let mut q = BTreeSet::new();
        cv.add_all_child_terms(&mut q, "MS:1002354")?;
        let mut e = BTreeSet::new();
        cv.add_all_child_terms(&mut e, "MS:1001872")?;
        cv.add_all_child_terms(&mut e, "MS:1002353")?;
        let mut specific = BTreeSet::new();
        cv.add_all_child_terms(&mut specific, "MS:1001143")?;
        Ok(Self { q, e, specific })
    }
}

/// The read-only context every read step shares: the PSI-MS vocabulary, the
/// caller's modification registry, the ceilings and the score subtrees.
struct Context<'a> {
    cv: &'a ControlledVocabulary,
    registry: &'a ModificationsDB,
    options: &'a ReadOptions,
    terms: ScoreTerms,
    /// The source's `xl_ms_search_`: true when any `AdditionalSearchParams`
    /// declares `MS:1002494`, which selects the cross-linking read path.
    crosslinking: bool,
}

/// One `SearchDatabase`, as the source's `DatabaseInput`.
#[derive(Clone, Debug, Default)]
struct DatabaseInput {
    location: String,
    version: String,
}

/// One `DBSequence`, as the source's `DBSequence` helper struct.
#[derive(Clone, Debug, Default)]
struct DbSequence {
    accession: String,
    sequence: String,
}

/// One `PeptideEvidence`, as the source's `PeptideEvidence` helper struct.
#[derive(Clone, Debug, Default)]
struct Evidence {
    start: Option<i32>,
    end: Option<i32>,
    pre: Option<char>,
    post: Option<char>,
    decoy: bool,
    db_sequence_ref: String,
}

/// The cross-link declarations of the `SequenceCollection`, as the source's
/// `xl_id_donor_map_`, `xl_id_acceptor_map_`, `xl_donor_pos_map_`,
/// `xl_acceptor_pos_map_`, `xl_mod_map_` and `xl_mass_map_`.
///
/// Every map keeps the first entry for a repeated key, because the source fills
/// `std::map`s with `insert`, which does not overwrite.
#[derive(Clone, Debug, Default)]
struct CrossLinks {
    /// `Peptide` id to the link value of its `crosslink donor` cvParam. A
    /// mono-link has no such cvParam and is registered under its own id.
    donor: BTreeMap<String, String>,
    /// `Peptide` id to the link value of its `crosslink acceptor` cvParam.
    acceptor: BTreeMap<String, String>,
    /// Link value to the donor's zero-based position in its peptide, which is
    /// `-1` for an N-terminal and `length` for a C-terminal link.
    donor_position: BTreeMap<String, i32>,
    /// Link value to the acceptor's zero-based position.
    acceptor_position: BTreeMap<String, i32>,
    /// `Peptide` id to the cross-linker name of its donor modification.
    modification: BTreeMap<String, String>,
    /// `Peptide` id to the `monoisotopicMassDelta` of its donor modification.
    mass: BTreeMap<String, f64>,
}

/// Everything declared once and referenced by id from the result section.
#[derive(Clone, Debug, Default)]
struct Library {
    databases: BTreeMap<String, DatabaseInput>,
    spectra_data: BTreeMap<String, String>,
    software: BTreeMap<String, (String, String)>,
    db_sequences: BTreeMap<String, DbSequence>,
    peptides: BTreeMap<String, AASequence>,
    evidences: BTreeMap<String, Evidence>,
    /// `peptide_ref` to its `PeptideEvidence` ids, in document order.
    by_peptide: BTreeMap<String, Vec<String>>,
    /// Empty unless the document declares `MS:1002494`.
    cross_links: CrossLinks,
}

fn element_id(node: &Node) -> Result<&str> {
    let id = node.get("id")?;
    if id.is_empty() {
        return Err(bad(format!("{} has an empty id", node.name)));
    }
    Ok(id)
}

fn unique<T>(map: &mut BTreeMap<String, T>, id: &str, kind: &str, value: T) -> Result<()> {
    if map.insert(id.to_owned(), value).is_some() {
        return Err(bad(format!("duplicate {kind} id {id:?}")));
    }
    Ok(())
}

/// True when any `AdditionalSearchParams` declares `MS:1002494`.
///
/// The source sets `xl_ms_search_` from exactly this scan, before it reads
/// anything else, and then takes a completely different read path for the whole
/// document. See [`read_crosslink_result`].
fn is_crosslinking(root: &Node) -> bool {
    gather(root, "AdditionalSearchParams")
        .iter()
        .flat_map(|node| node.children.iter())
        .any(|child| attribute(child, "accession") == CROSSLINKING_SEARCH)
}

fn read_inputs(root: &Node, library: &mut Library, budget: &mut Budget) -> Result<()> {
    for node in gather(root, "SpectraData") {
        let id = element_id(node)?;
        budget.text(id)?;
        unique(
            &mut library.spectra_data,
            id,
            "SpectraData",
            attribute(node, "location").to_owned(),
        )?;
    }
    for node in gather(root, "SearchDatabase") {
        let id = element_id(node)?;
        budget.text(id)?;
        // DatabaseName is not read: the source resolves it into
        // DatabaseInput::name, substituting "unknown" when the element is
        // missing, and then never consults that field. Only the location and
        // the version reach SearchParameters.
        unique(
            &mut library.databases,
            id,
            "SearchDatabase",
            DatabaseInput {
                location: attribute(node, "location").to_owned(),
                version: attribute(node, "version").to_owned(),
            },
        )?;
    }
    Ok(())
}

fn read_software(
    root: &Node,
    library: &mut Library,
    context: &Context<'_>,
    budget: &mut Budget,
) -> Result<()> {
    let (cv, options) = (context.cv, context.options);
    let mut software_terms = BTreeSet::new();
    cv.add_all_child_terms(&mut software_terms, "MS:1000531")?;
    for node in gather(root, "AnalysisSoftware") {
        let id = element_id(node)?;
        budget.text(id)?;
        let mut name = String::new();
        let mut version = String::new();
        // The source reads the version attribute inside the SoftwareName
        // branch, so software without a SoftwareName child keeps no version
        // and is dropped below.
        for element in children(node, "SoftwareName") {
            version = attribute(node, "version").to_owned();
            let group = param_group(element, options, budget)?;
            if !group.cv.is_empty() {
                for (accession, params) in &group.cv {
                    if software_terms.contains(accession) {
                        name = params[0].name.clone();
                        break;
                    }
                }
            } else {
                for (key, value) in &group.user {
                    if key.contains("name") {
                        name = text_of(value);
                        break;
                    }
                    name = key.clone();
                }
            }
        }
        if name.is_empty() || version.is_empty() {
            // Source: "No name/version found for 'AnalysisSoftware'" and skip.
            continue;
        }
        unique(
            &mut library.software,
            id,
            "AnalysisSoftware",
            (name, version),
        )?;
    }
    Ok(())
}

/// Text of a metadata value, as `DataValue::toString` would render it.
fn text_of(value: &MetaValue) -> String {
    match value.data() {
        MetaValueData::Empty => String::new(),
        MetaValueData::String(text) => text.clone(),
        MetaValueData::Integer(number) => number.to_string(),
        MetaValueData::Float(number) => number.to_string(),
        MetaValueData::StringList(items) => items.join(", "),
        MetaValueData::IntegerList(items) => items
            .iter()
            .map(i64::to_string)
            .collect::<Vec<_>>()
            .join(", "),
        MetaValueData::FloatList(items) => items
            .iter()
            .map(f64::to_string)
            .collect::<Vec<_>>()
            .join(", "),
    }
}

fn read_sequence_collection(
    root: &Node,
    library: &mut Library,
    context: &Context<'_>,
    budget: &mut Budget,
) -> Result<()> {
    let options = context.options;
    for node in gather(root, "DBSequence") {
        let id = element_id(node)?;
        let accession = node.get("accession")?;
        if accession.is_empty() {
            // Source: a DBSequence without an accession is silently not
            // registered, so every PeptideEvidence pointing at it later
            // resolves to a default-constructed empty protein.
            return Err(bad(format!("DBSequence {id:?} has an empty accession")));
        }
        let sequence = child(node, "Seq").map(|s| s.text.trim().to_owned());
        budget.text(id)?;
        budget.text(accession)?;
        if let Some(text) = &sequence {
            budget.text(text)?;
        }
        unique(
            &mut library.db_sequences,
            id,
            "DBSequence",
            DbSequence {
                accession: accession.to_owned(),
                sequence: sequence.unwrap_or_default(),
            },
        )?;
    }
    let peptides = gather(root, "Peptide");
    if peptides.len() > options.max_list_items {
        return Err(bad("mzIdentML Peptide count exceeds the configured limit"));
    }
    let mut cross_links = CrossLinks::default();
    for node in peptides {
        let id = element_id(node)?;
        budget.text(id)?;
        let sequence = read_peptide(node, context, budget, &mut cross_links)?;
        unique(&mut library.peptides, id, "Peptide", sequence)?;
    }
    library.cross_links = cross_links;
    let evidences = gather(root, "PeptideEvidence");
    if evidences.len() > options.max_list_items {
        return Err(bad(
            "mzIdentML PeptideEvidence count exceeds the configured limit",
        ));
    }
    for node in evidences {
        let id = element_id(node)?;
        let peptide_ref = node.get("peptide_ref")?.to_owned();
        budget.text(id)?;
        budget.text(&peptide_ref)?;
        // start/end are optional; the source parses them together and leaves
        // both unknown when either conversion throws.
        let (start, end) = match (node.optional("start"), node.optional("end")) {
            (Some(start), Some(end)) => (
                Some(integer(start, "PeptideEvidence start")?),
                Some(integer(end, "PeptideEvidence end")?),
            ),
            _ => (None, None),
        };
        let single = |value: Option<&str>, what: &str| -> Result<Option<char>> {
            match value {
                None => Ok(None),
                Some(text) => {
                    let mut characters = text.chars();
                    match (characters.next(), characters.next()) {
                        // The source takes text[0] of a possibly empty string,
                        // which yields the null character for pre="".
                        (Some(c), None) => Ok(Some(c)),
                        _ => Err(bad(format!(
                            "PeptideEvidence {what} must be exactly one character"
                        ))),
                    }
                }
            }
        };
        let decoy = match node.optional("isDecoy") {
            None => false,
            // Source: any value starting with 't' or '1' is a decoy.
            Some(text) => text.starts_with('t') || text.starts_with('1'),
        };
        let evidence = Evidence {
            start,
            end,
            pre: single(node.optional("pre"), "pre")?,
            post: single(node.optional("post"), "post")?,
            decoy,
            db_sequence_ref: node.get("dBSequence_ref")?.to_owned(),
        };
        unique(&mut library.evidences, id, "PeptideEvidence", evidence)?;
        let list = library.by_peptide.entry(peptide_ref).or_default();
        if list.len() >= options.max_list_items {
            return Err(bad("mzIdentML peptide evidence list limit exceeded"));
        }
        list.push(id.to_owned());
    }
    Ok(())
}

/// Build one `Peptide`'s sequence, as `parsePeptideSiblings_`.
///
/// # Errors
///
/// [`Error::Parse`] when `PeptideSequence` is missing or empty,
/// [`Error::InvalidRange`] when a `SubstitutionModification` or `Modification`
/// `location` falls outside the sequence, and [`Error::InvalidValue`] when a
/// modification cannot be resolved against `registry`.
///
/// The source recovers from every one of those by logging and inserting an
/// empty `AASequence` into its peptide map, so the affected PSMs silently lose
/// their sequence; one of the recoveries is an out-of-bounds write (see
/// `OpenMS_CPP_ISSUES.md`). This port refuses instead.
///
/// In a cross-linking document the `Modification` children take the source's
/// `xl_ms_search_` branch: they also fill `links`, and a modification name that
/// is not in `registry` is skipped rather than refused, because the source
/// guards that branch with `ModificationsDB::has` and only warns. Without that
/// leniency no XL-MS document would load at all: its cross-linker cvParams name
/// XLMOD records, which are not in the modification database.
fn read_peptide(
    node: &Node,
    context: &Context<'_>,
    budget: &mut Budget,
    links: &mut CrossLinks,
) -> Result<AASequence> {
    let (registry, options) = (context.registry, context.options);
    let text = child(node, "PeptideSequence")
        .ok_or_else(|| bad("Peptide requires a PeptideSequence child"))?
        .text
        .trim()
        .to_owned();
    if text.is_empty() {
        // Source: DOMNode::getFirstChild() is null for an empty element and is
        // dereferenced without a check.
        return Err(bad("PeptideSequence must not be empty"));
    }
    budget.text(&text)?;
    let mut residues: Vec<char> = text.chars().collect();
    for element in children(node, "SubstitutionModification") {
        let original = attribute(element, "originalResidue");
        let replacement = attribute(element, "replacementResidue");
        let mut characters = replacement.chars();
        let replacement = match (characters.next(), characters.next()) {
            (Some(c), None) => c,
            _ => {
                return Err(bad(
                    "SubstitutionModification replacementResidue must be one character",
                ));
            }
        };
        match element.optional("location") {
            Some(location) => {
                let location = integer(location, "SubstitutionModification location")?;
                let index = usize::try_from(location)
                    .ok()
                    .and_then(|value| value.checked_sub(1))
                    .filter(|index| *index < residues.len())
                    .ok_or_else(|| {
                        Error::InvalidRange(format!(
                            "SubstitutionModification location {location} is outside the peptide"
                        ))
                    })?;
                residues[index] = replacement;
            }
            None if original.chars().count() == 1 => {
                let original = original.chars().next().unwrap_or('\0');
                if !residues.contains(&original) {
                    return Err(bad(
                        "SubstitutionModification originalResidue does not occur in the peptide",
                    ));
                }
                for residue in &mut residues {
                    if *residue == original {
                        *residue = replacement;
                    }
                }
            }
            None => {
                return Err(bad(
                    "SubstitutionModification without location needs one originalResidue",
                ));
            }
        }
    }
    let plain: String = residues.into_iter().collect();
    let mut sequence = AASequence::parse_with_registry(&plain, registry)?;
    let length = sequence.len();
    let mut count = 0usize;
    for element in children(node, "Modification") {
        count += 1;
        if count > options.max_list_items {
            return Err(bad("mzIdentML Modification count exceeds the limit"));
        }
        budget.spend(64, 256)?;
        let declared = match element.optional("location") {
            Some(text) => Some(integer(text, "Modification location")?),
            None => None,
        };
        let limit = i32::try_from(length)
            .ok()
            .and_then(|value| value.checked_add(1))
            .ok_or_else(|| invalid("peptide is too long for an mzIdentML location"))?;
        let location = match declared {
            Some(value) if (0..=limit).contains(&value) => value,
            _ => match infer_modification_location(element, registry, length)? {
                Some(value) => value,
                // Source: "Skipping modification with missing or invalid
                // 'location' attribute; its position could not be inferred."
                None => continue,
            },
        };
        if context.crosslinking {
            read_crosslink_modification(
                element,
                element_id(node)?,
                location,
                &mut sequence,
                length,
                registry,
                links,
            )?;
        } else {
            apply_modification(&mut sequence, element, location, length, registry)?;
        }
    }
    Ok(sequence)
}

/// One `Peptide/Modification` of a cross-linking document, as the
/// `xl_ms_search_` branch of `parsePeptideSiblings_`.
///
/// `location` is the mzIdentML convention (0 is the N-terminus); the positions
/// recorded in `links` are `location - 1`, so an N-terminal link is `-1` and a
/// C-terminal one is the peptide length, which the item reader turns back into
/// a residue index plus a terminal specificity.
///
/// # Errors
///
/// [`Error::InvalidValue`] when a `monoisotopicMassDelta` is not a finite
/// number, and whatever [`apply_modification`] would return for a residue
/// modification whose name resolves but does not fit the residue - the source
/// lets that exception escape for an internal position and swallows it for the
/// N-terminal one, which this reproduces.
fn read_crosslink_modification(
    element: &Node,
    peptide_id: &str,
    location: i32,
    sequence: &mut AASequence,
    length: usize,
    registry: &ModificationsDB,
    links: &mut CrossLinks,
) -> Result<()> {
    let position = location.saturating_sub(1);
    let params = || element.children.iter().filter(|c| c.name == "cvParam");
    let mut donor_or_acceptor = false;
    let mut crosslinker = false;
    for node in params() {
        let param = CvParam::from_node(node);
        if param.accession == CROSSLINK_DONOR {
            links
                .donor
                .entry(peptide_id.to_owned())
                .or_insert_with(|| param.value.clone());
            let mass = match element.optional("monoisotopicMassDelta") {
                Some(text) if !text.is_empty() => finite(text)?,
                _ => 0.0,
            };
            links.mass.entry(peptide_id.to_owned()).or_insert(mass);
            links
                .donor_position
                .entry(param.value.clone())
                .or_insert(position);
            // The cross-linker name is the first UNIMOD or XLMOD cvParam of
            // this Modification; the source falls back to the first element
            // child's name attribute, and dereferences a null pointer when
            // there is none, which cannot happen here.
            let name = params()
                .find(|node| matches!(attribute(node, "cvRef"), "UNIMOD" | "XLMOD"))
                .or_else(|| params().next())
                .map(|node| attribute(node, "name").to_owned())
                .unwrap_or_default();
            links
                .modification
                .entry(peptide_id.to_owned())
                .or_insert(name);
            donor_or_acceptor = true;
        } else if param.accession == CROSSLINK_ACCEPTOR {
            links
                .acceptor
                .entry(peptide_id.to_owned())
                .or_insert_with(|| param.value.clone());
            links
                .acceptor_position
                .entry(param.value.clone())
                .or_insert(position);
            donor_or_acceptor = true;
        } else {
            if param.name.starts_with("Xlink") || param.accession.starts_with("XLMOD") {
                crosslinker = true;
            }
            if param.name.contains("unknown mono-link") {
                links
                    .modification
                    .entry(peptide_id.to_owned())
                    .or_insert_with(|| param.name.clone());
            } else {
                apply_crosslink_residue_modification(sequence, &param, location, length, registry)?;
            }
        }
    }
    if !donor_or_acceptor && crosslinker {
        // A mono-link carries no cross-linking CV term, so the source uses the
        // Peptide id itself as the link value.
        links
            .donor
            .entry(peptide_id.to_owned())
            .or_insert_with(|| peptide_id.to_owned());
        links
            .donor_position
            .entry(peptide_id.to_owned())
            .or_insert(position);
    }
    Ok(())
}

/// Apply one non-cross-linking `cvParam` of an XL-MS `Modification`.
///
/// Unlike [`apply_modification`], a name the registry does not know is skipped
/// with no error: the source guards this branch with `ModificationsDB::has` and
/// logs "Modification: X not found in ModificationsDB".
fn apply_crosslink_residue_modification(
    sequence: &mut AASequence,
    param: &CvParam,
    location: i32,
    length: usize,
    registry: &ModificationsDB,
) -> Result<()> {
    if param.cv_ref != "UNIMOD" && param.cv_ref != "XLMOD" {
        // e.g. MS:1001524 "fragment neutral loss" is not a modification.
        return Ok(());
    }
    let name = if param.name == "unknown modification" {
        param.value.clone()
    } else {
        param.name.clone()
    };
    if name.is_empty() || registry.find(&name, None, None).is_empty() {
        return Ok(());
    }
    let terminal_c = i32::try_from(length).map(|value| value + 1).unwrap_or(-1);
    if location == 0 {
        // Source: the N-terminal branch swallows every exception.
        let _ = sequence.set_n_terminal_modification_with_registry(&name, registry);
        Ok(())
    } else if location == terminal_c {
        sequence.set_c_terminal_modification_with_registry(&name, registry)
    } else {
        let index = usize::try_from(location)
            .ok()
            .and_then(|value| value.checked_sub(1))
            .filter(|index| *index < length)
            .ok_or_else(|| {
                Error::InvalidRange(format!(
                    "Modification location {location} is outside the peptide"
                ))
            })?;
        sequence.set_modification_with_registry(index, &name, registry)
    }
}

/// Infer a terminal `location` for a `Modification` that declares none, as
/// `inferModificationLocation_`.
///
/// Returns `Some(0)` for a modification with an N-terminal variant,
/// `Some(length + 1)` for one that is exclusively C-terminal, and `None` when
/// no unique terminal position follows - including whenever `residues` names a
/// concrete amino acid, because then the position is residue-specific and
/// genuinely unknown.
fn infer_modification_location(
    element: &Node,
    registry: &ModificationsDB,
    length: usize,
) -> Result<Option<i32>> {
    let residues = attribute(element, "residues");
    if !residues.is_empty() && residues != "." {
        return Ok(None);
    }
    for param in element.children.iter().filter(|c| c.name == "cvParam") {
        let name = attribute(param, "name");
        if name.is_empty() {
            continue;
        }
        let matches = registry.find(name, None, None);
        if matches.is_empty() {
            // Source guards with has() first so unrelated cvParams (neutral
            // loss terms, for instance) do not log a "not found" warning.
            continue;
        }
        let mut n_terminal = false;
        let mut c_terminal = false;
        let mut internal = false;
        for record in matches {
            match record.term_specificity() {
                TermSpecificity::NTerm | TermSpecificity::ProteinNTerm => n_terminal = true,
                TermSpecificity::CTerm | TermSpecificity::ProteinCTerm => c_terminal = true,
                TermSpecificity::Anywhere => internal = true,
            }
        }
        if n_terminal {
            return Ok(Some(0));
        }
        if c_terminal && !internal {
            let location = i32::try_from(length)
                .ok()
                .and_then(|value| value.checked_add(1))
                .ok_or_else(|| invalid("peptide is too long for an mzIdentML location"))?;
            return Ok(Some(location));
        }
    }
    Ok(None)
}

/// Apply one `Modification` element at a validated `location`.
///
/// `location` is the mzIdentML convention: 0 is the N-terminus, 1 the first
/// residue and `length + 1` the C-terminus.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the named modification is not in `registry`.
/// The source logs "Modification: X not found in ModificationsDB" and keeps the
/// unmodified residue for the internal and C-terminal cases, and for the
/// N-terminal case swallows every exception from `setNTerminalModification`.
fn apply_modification(
    sequence: &mut AASequence,
    element: &Node,
    location: i32,
    length: usize,
    registry: &ModificationsDB,
) -> Result<()> {
    let terminal_c = i32::try_from(length).map(|value| value + 1).unwrap_or(-1);
    for param in element.children.iter().filter(|c| c.name == "cvParam") {
        let param = CvParam::from_node(param);
        let unknown = param.accession == "MS:1001460" || param.name == "unknown modification";
        let name = if unknown {
            // "unknown modification" carries the actual name in value=.
            if param.value.is_empty() {
                return Err(invalid(
                    "unknown modification cvParam carries no value naming the modification",
                ));
            }
            param.value.clone()
        } else {
            if param.cv_ref != "UNIMOD" && param.cv_ref != "XLMOD" {
                // e.g. MS:1001524 "fragment neutral loss" is not a modification.
                continue;
            }
            param.name.clone()
        };
        if location == 0 {
            sequence.set_n_terminal_modification_with_registry(&name, registry)?;
        } else if location == terminal_c {
            sequence.set_c_terminal_modification_with_registry(&name, registry)?;
        } else {
            let index = usize::try_from(location)
                .ok()
                .and_then(|value| value.checked_sub(1))
                .filter(|index| *index < length)
                .ok_or_else(|| {
                    Error::InvalidRange(format!(
                        "Modification location {location} is outside the peptide"
                    ))
                })?;
            sequence.set_modification_with_registry(index, &name, registry)?;
        }
    }
    Ok(())
}
/// Identification runs, their `spectrumIdentificationProtocol_ref` in the same
/// order, and the `SpectrumIdentificationList` id of each.
type Runs = (
    Vec<ProteinIdentification>,
    Vec<String>,
    BTreeMap<String, usize>,
);

fn read_runs(
    root: &Node,
    library: &Library,
    options: &ReadOptions,
    crosslinking: bool,
    budget: &mut Budget,
) -> Result<Runs> {
    let nodes = gather(root, "SpectrumIdentification");
    if nodes.is_empty() {
        return Err(missing("mzIdentML has no SpectrumIdentification element"));
    }
    if nodes.len() > options.max_list_items {
        return Err(bad(
            "mzIdentML SpectrumIdentification count exceeds the limit",
        ));
    }
    let mut runs = Vec::new();
    let mut links = Vec::new();
    let mut list_to_run = BTreeMap::new();
    let mut identifiers = BTreeSet::new();
    for node in nodes {
        let id = element_id(node)?;
        budget.text(id)?;
        if !identifiers.insert(id.to_owned()) {
            return Err(bad(format!("duplicate SpectrumIdentification id {id:?}")));
        }
        // The source takes the last InputSpectra/SearchDatabaseRef child.
        let spectra_data_ref = children(node, "InputSpectra")
            .map(|c| attribute(c, "spectraData_ref"))
            .last()
            .unwrap_or_default();
        let database_ref = children(node, "SearchDatabaseRef")
            .map(|c| attribute(c, "searchDatabase_ref"))
            .last()
            .unwrap_or_default();
        // A dangling searchDatabase_ref or spectraData_ref leaves the source's
        // std::map::operator[] with a default-constructed entry, i.e. empty
        // strings; upstream test data relies on that (MzIdentMLFile_msgf_mini
        // references SearchDB_99, which it never declares), so both stay
        // tolerated here rather than becoming errors.
        let database = library
            .databases
            .get(database_ref)
            .cloned()
            .unwrap_or_default();
        let location = library
            .spectra_data
            .get(spectra_data_ref)
            .cloned()
            .unwrap_or_default();
        let mut run = ProteinIdentification {
            identifier: id.to_owned(),
            search_parameters: SearchParameters {
                database: database.location,
                database_version: database.version,
                ..Default::default()
            },
            ..Default::default()
        };
        run.metadata.insert(
            "spectra_data".into(),
            MetaValue::new(MetaValueData::StringList(vec![location]))?,
        );
        if crosslinking {
            // The source tags every run of an XL-MS document with the
            // crosslinking search term, which is also what its writer tests to
            // pick the cross-linking output path.
            run.metadata.insert(
                "SpectrumIdentificationProtocol".into(),
                CROSSLINKING_SEARCH.to_owned().into(),
            );
        }
        // Source: activityDate, or DateTime::now() when absent, which makes a
        // load unreproducible. An absent date stays absent here.
        run.date_time = node
            .optional("activityDate")
            .filter(|text| !text.is_empty())
            .map(str::to_owned);
        let list_ref = attribute(node, "spectrumIdentificationList_ref").to_owned();
        if !list_ref.is_empty() && list_to_run.insert(list_ref.clone(), runs.len()).is_some() {
            // Source: std::map::insert keeps the first index, so the second run
            // silently receives no peptide identifications at all.
            return Err(bad(format!(
                "two SpectrumIdentification elements reference SpectrumIdentificationList {list_ref:?}"
            )));
        }
        let _ = list_ref;
        links.push(attribute(node, "spectrumIdentificationProtocol_ref").to_owned());
        runs.push(run);
    }
    Ok((runs, links, list_to_run))
}

/// `AdditionalSearchParams` to [`SearchParameters`], as `findSearchParameters_`.
fn additional_search_params(group: &ParamGroup) -> Result<SearchParameters> {
    let mut parameters = SearchParameters::default();
    for (accession, params) in &group.cv {
        for param in params {
            parameters.metadata.insert(
                accession.clone(),
                MetaValue::new(MetaValueData::String(param.value.clone()))?,
            );
        }
    }
    let mut min_charge = 0i32;
    let mut max_charge = 0i32;
    for (name, value) in &group.user {
        match name.as_str() {
            "taxonomy" => parameters.taxonomy = text_of(value),
            "charges" => parameters.charges = text_of(value),
            "MinCharge" => min_charge = integer(&text_of(value), "MinCharge")?,
            "MaxCharge" => max_charge = integer(&text_of(value), "MaxCharge")?,
            "NumTolerableTermini" => {
                parameters.enzyme_specificity =
                    match integer(&text_of(value), "NumTolerableTermini")? {
                        0 => EnzymeTermSpecificity::None,
                        1 => EnzymeTermSpecificity::Semi,
                        2 => EnzymeTermSpecificity::Full,
                        3 => EnzymeTermSpecificity::Unknown,
                        // Source casts the integer to the enum unchecked, which
                        // yields a value outside EnzymaticDigestion::Specificity.
                        other => {
                            return Err(invalid(format!(
                                "NumTolerableTermini {other} is not a digestion specificity"
                            )));
                        }
                    }
            }
            _ => {
                parameters.metadata.insert(name.clone(), value.clone());
            }
        }
    }
    if min_charge != 0 || max_charge != 0 {
        // MinCharge/MaxCharge take precedence over an explicit charges list.
        parameters.charges = format!("{min_charge}-{max_charge}");
    }
    Ok(parameters)
}

/// Resolve one modification name the way `ModificationsDB::getModification`
/// does: when a residue is given and no terminal specificity is requested, an
/// `ANYWHERE` lookup is tried first, "to avoid ambiguities (e.g.
/// `Carbamidomethyl (N-term)`/`Carbamidomethyl (C)`)" as the source comment
/// puts it, and only then the unrestricted one.
///
/// # Errors
///
/// [`Error::InvalidValue`] when nothing matches, and also when the remaining
/// match is ambiguous: the crate's registry refuses that, while the source
/// warns "picking the first one only" and continues with an arbitrary record.
fn resolve_modification<'a>(
    registry: &'a ModificationsDB,
    name: &str,
    residue: Option<char>,
    term: Option<TermSpecificity>,
) -> Result<&'a ResidueModification> {
    if residue.is_some() && term.is_none() {
        if let Ok(record) =
            registry.get_modification(name, residue, Some(TermSpecificity::Anywhere))
        {
            return Ok(record);
        }
    }
    registry.get_modification(name, residue, term)
}

/// The `fixedMod`/`variable` modification names of one `ModificationParams`.
fn modification_params(
    node: &Node,
    context: &Context<'_>,
    budget: &mut Budget,
) -> Result<(Vec<String>, Vec<String>)> {
    let (registry, options) = (context.registry, context.options);
    let mut fixed = Vec::new();
    let mut variable = Vec::new();
    for element in children(node, "SearchModification") {
        if fixed.len().saturating_add(variable.len()) >= options.max_list_items {
            return Err(bad("mzIdentML SearchModification count exceeds the limit"));
        }
        budget.spend(64, 512)?;
        let is_fixed = matches!(attribute(element, "fixedMod"), "true" | "1");
        let residues = attribute(element, "residues");
        let mut name = String::new();
        let mut rules: BTreeSet<String> = BTreeSet::new();
        for sub in &element.children {
            match sub.name.as_str() {
                "cvParam" => {
                    let param = CvParam::from_node(sub);
                    let named = if param.name == "unknown modification" {
                        param.value.clone()
                    } else {
                        param.name.clone()
                    };
                    if param.cv_ref == "UNIMOD" || param.cv_ref == "XLMOD" {
                        name = named;
                    } else if name.is_empty() {
                        // Fallback for files that do not set cvRef to the
                        // modification ontology, e.g. cvRef="PSI-MS" with a
                        // UNIMOD accession.
                        name = named;
                    }
                }
                "SpecificityRules" => {
                    for rule in param_group(sub, options, budget)?.cv.keys() {
                        rules.insert(rule.clone());
                    }
                }
                _ => {}
            }
        }
        if name.is_empty() {
            continue;
        }
        let residue = if residues == "." || residues.is_empty() {
            None
        } else {
            let mut characters = residues.chars();
            match (characters.next(), characters.next()) {
                (Some(c), None) => Some(c),
                // Several residues mean the modified one is unknown; the
                // source passes the whole string to ModificationsDB, where the
                // lookup fails and the modification is dropped with a warning.
                _ => continue,
            }
        };
        // The last matching rule wins, as in the source's rule loop.
        let mut term = None;
        for rule in &rules {
            match rule.as_str() {
                "MS:1001189" => term = Some(TermSpecificity::NTerm),
                "MS:1001190" => term = Some(TermSpecificity::CTerm),
                "MS:1002057" => term = Some(TermSpecificity::ProteinNTerm),
                "MS:1002058" => term = Some(TermSpecificity::ProteinCTerm),
                _ => {}
            }
        }
        let Ok(record) = resolve_modification(registry, &name, residue, term) else {
            // Source: "SearchModification 'X' not found in ModificationsDB".
            continue;
        };
        if is_fixed {
            fixed.push(record.full_id().to_owned());
        } else {
            variable.push(record.full_id().to_owned());
        }
    }
    Ok((fixed, variable))
}

fn read_protocols(
    root: &Node,
    library: &Library,
    runs: &mut [ProteinIdentification],
    links: &[String],
    context: &Context<'_>,
    budget: &mut Budget,
) -> Result<()> {
    let (cv, options) = (context.cv, context.options);
    let nodes = gather(root, "SpectrumIdentificationProtocol");
    if nodes.is_empty() {
        return Err(missing(
            "mzIdentML has no SpectrumIdentificationProtocol element",
        ));
    }
    let mut threshold_terms = BTreeSet::new();
    cv.all_child_terms("MS:1002482")?
        .into_iter()
        .for_each(|term| {
            threshold_terms.insert(term);
        });
    let mut seen = BTreeSet::new();
    for node in nodes {
        let id = element_id(node)?;
        budget.text(id)?;
        if !seen.insert(id.to_owned()) {
            return Err(bad(format!(
                "duplicate SpectrumIdentificationProtocol id {id:?}"
            )));
        }
        // The source assigns the AdditionalSearchParams block over the
        // accumulated one, so a document that places AdditionalSearchParams
        // after ModificationParams, Enzymes or the tolerances loses them. The
        // schema fixes that order, so this seeds from it first instead.
        let mut parameters = match child(node, "AdditionalSearchParams") {
            Some(element) => additional_search_params(&param_group(element, options, budget)?)?,
            None => SearchParameters::default(),
        };
        if let Some(element) = child(node, "ModificationParams") {
            let (fixed, variable) = modification_params(element, context, budget)?;
            parameters.fixed_modifications = fixed;
            parameters.variable_modifications = variable;
        }
        let mut enzyme_terms = BTreeSet::new();
        cv.add_all_child_terms(&mut enzyme_terms, "MS:1001045")?;
        if let Some(element) = child(node, "Enzymes") {
            for enzyme in children(element, "Enzyme") {
                let missed = match enzyme.optional("missedCleavages") {
                    Some(text) => integer(text, "missedCleavages")?,
                    None => -1,
                };
                // Source: a negative value is "assumed unlimited" and becomes
                // 1000; an unreadable one takes the same path.
                parameters.missed_cleavages = u32::try_from(missed).unwrap_or(1000);
                let mut name = String::new();
                for element in children(enzyme, "EnzymeName") {
                    for (accession, params) in &param_group(element, options, budget)?.cv {
                        if enzyme_terms.contains(accession) {
                            name = params[0].name.clone();
                        }
                    }
                }
                if ProteaseDB::global().has_enzyme(&name) {
                    parameters.digestion_enzyme = name;
                }
            }
        }
        for (tag, ppm_flag) in [("FragmentTolerance", true), ("ParentTolerance", false)] {
            let Some(element) = child(node, tag) else {
                continue;
            };
            let group = param_group(element, options, budget)?;
            let mut value: Option<f64> = None;
            let mut ppm = false;
            for params in group.cv.values() {
                // The source keeps the numerically greater of the plus/minus
                // bounds and does not require them to agree.
                let bound = finite(&params[0].value)?;
                value = Some(value.map_or(bound, |current: f64| current.max(bound)));
                if params[0].unit_name == "parts per million" {
                    ppm = true;
                }
            }
            if let Some(value) = value {
                if value < 0.0 {
                    return Err(invalid("mzIdentML search tolerances must be nonnegative"));
                }
                let tolerance = if ppm {
                    Tolerance::Ppm(value)
                } else {
                    Tolerance::Absolute(value)
                };
                if ppm_flag {
                    parameters.fragment_tolerance = tolerance;
                } else {
                    parameters.precursor_tolerance = tolerance;
                }
            }
        }
        let mut threshold = None;
        if let Some(element) = child(node, "Threshold") {
            for (accession, params) in &param_group(element, options, budget)?.cv {
                if threshold_terms.contains(accession) {
                    if accession != "MS:1001494" {
                        threshold = Some(finite(&params[0].value)?);
                    }
                    break;
                }
            }
        }
        let software = library
            .software
            .get(attribute(node, "analysisSoftware_ref"))
            .cloned()
            .unwrap_or_default();
        for (run, protocol_ref) in runs.iter_mut().zip(links) {
            if protocol_ref != id {
                continue;
            }
            run.search_engine = software.0.clone();
            run.search_engine_version = software.1.clone();
            let mut applied = parameters.clone();
            // The run's database came from its own SearchDatabaseRef and must
            // survive the protocol-wide parameter block.
            applied.database = run.search_parameters.database.clone();
            applied.database_version = run.search_parameters.database_version.clone();
            run.search_parameters = applied;
            if let Some(threshold) = threshold {
                run.significance_threshold = threshold;
            }
        }
    }
    Ok(())
}
/// Sort protein hits as `ProteinIdentification::sort` does.
///
/// `ProteinHit::ScoreMore`/`ScoreLess` compare the tuple `(score, accession)`,
/// so equal scores are ordered by accession, descending for
/// `higher_score_better`. The crate's
/// [`ProteinIdentification::sort`](crate::identification::ProteinIdentification::sort)
/// compares the score alone, which leaves equal-scoring hits in insertion
/// order; the reader therefore applies the source's full key itself, so a load
/// reproduces the C++ hit order.
fn sort_protein_hits(run: &mut ProteinIdentification) {
    let higher = run.higher_score_better;
    run.hits.sort_by(|a, b| {
        let (first, second) = if higher { (b, a) } else { (a, b) };
        first
            .score
            .partial_cmp(&second.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| first.accession.cmp(&second.accession))
    });
}

/// The score type, orientation and value of one `SpectrumIdentificationItem`.
struct Score {
    value: f64,
    higher_is_better: bool,
    score_type: String,
}

/// Pick the PSM score exactly as `parseSpectrumIdentificationItemElement_`.
///
/// The scan walks the `cvParam`s in lexicographic accession order and takes the
/// first q-value, then any PSM-level search-engine statistic, then any e-value.
/// `MS:1002055` (distinct peptide-level q-value) is skipped without stopping
/// the scan, and the bare parent term `MS:1001143` is a last resort that does
/// not stop it either. A PSM with none of those yields no hit at all, which is
/// why a `SpectrumIdentificationResult` can legitimately produce an empty
/// [`PeptideIdentification`].
fn select_score(
    group: &ParamGroup,
    terms: &ScoreTerms,
    cv: &ControlledVocabulary,
) -> Result<Option<Score>> {
    let mut found: Option<Score> = None;
    for (accession, params) in &group.cv {
        let param = &params[0];
        if terms.q.contains(accession) || accession == "MS:1002354" {
            if accession != "MS:1002055" {
                found = Some(Score {
                    value: score_value(&param.value)?,
                    higher_is_better: false,
                    score_type: "q-value".into(),
                });
                break;
            }
        } else if accession != "MS:1001143" && terms.specific.contains(accession) {
            found = Some(Score {
                value: score_value(&param.value)?,
                higher_is_better: cv.get_term(accession)?.is_higher_better_score(),
                score_type: param.name.clone(),
            });
            break;
        } else if terms.e.contains(accession) {
            found = Some(Score {
                value: score_value(&param.value)?,
                higher_is_better: false,
                score_type: "E-value".into(),
            });
            break;
        } else if accession == "MS:1001143" {
            // The parent term has no value; the source assumes higher is better.
            found = Some(Score {
                value: 0.0,
                higher_is_better: true,
                score_type: "PSM-level search engine specific statistic".into(),
            });
        }
    }
    Ok(found)
}

/// Merge one evidence's decoy state into the hit's `target_decoy` value.
fn merge_target_decoy(metadata: &mut MetaInfo, decoy: bool) -> Result<()> {
    let wanted = if decoy { "decoy" } else { "target" };
    let value = match metadata.get(user_param::TARGET_DECOY) {
        Some(current) if current.as_str().unwrap_or_default() != wanted => "target+decoy",
        _ => wanted,
    };
    metadata.insert(
        user_param::TARGET_DECOY.into(),
        MetaValue::new(MetaValueData::String(value.into()))?,
    );
    Ok(())
}

fn read_item(
    node: &Node,
    identification: &mut PeptideIdentification,
    run: &mut ProteinIdentification,
    library: &Library,
    context: &Context<'_>,
    budget: &mut Budget,
) -> Result<()> {
    let (cv, options) = (context.cv, context.options);
    let group = param_group(node, options, budget)?;
    let Some(score) = select_score(&group, &context.terms, cv)? else {
        // No recognised score: the source reads no hit from this item.
        return Ok(());
    };
    let charge = integer(node.get("chargeState")?, "chargeState")?;
    // mzIdentML ranks are 1-based. The source treats rank 0 as 1 for PMF data
    // and otherwise subtracts 1 into OpenMS's 0-based rank; an unreadable rank
    // stays 0 there and underflows the unsigned member (see OpenMS_CPP_ISSUES).
    let declared = integer(node.get("rank")?, "rank")?;
    if declared < 0 {
        return Err(invalid(
            "SpectrumIdentificationItem rank must not be negative",
        ));
    }
    let rank = u32::try_from(declared.max(1) - 1)
        .map_err(|_| invalid("SpectrumIdentificationItem rank exceeds the rank range"))?;
    let peptide_ref = node.get("peptide_ref")?;
    let sequence = library.peptides.get(peptide_ref).cloned().ok_or_else(|| {
        // Source: std::map::operator[] default-constructs an empty AASequence,
        // so the PSM silently loses its sequence.
        bad(format!(
            "SpectrumIdentificationItem references undeclared Peptide {peptide_ref:?}"
        ))
    })?;
    let pass = matches!(node.get("passThreshold")?, "true" | "1");
    identification.higher_score_better = score.higher_is_better;
    identification.score_type = score.score_type;
    let mut hit = PeptideHit {
        sequence,
        score: score.value,
        rank,
        charge,
        ..Default::default()
    };
    for (accession, params) in &group.cv {
        for param in params {
            if accession == "MS:1001143" {
                // The parent statistic term carries no value.
                continue;
            }
            budget.text(accession)?;
            let value = if accession == "MS:1002540" {
                // Source keeps this one as text rather than typing it.
                MetaValue::new(MetaValueData::String(param.value.clone()))?
            } else {
                cv_param_value(cv, param)?
            };
            hit.metadata.insert(accession.clone(), value);
        }
    }
    for (name, value) in &group.user {
        budget.text(name)?;
        hit.metadata.insert(name.clone(), value.clone());
    }
    if let Some(calculated) = optional_value(attribute(node, "calculatedMassToCharge"))? {
        hit.metadata.insert(
            "calcMZ".into(),
            MetaValue::new(MetaValueData::Float(calculated))?,
        );
    }
    // The experimental m/z belongs to the spectrum, not the candidate, but the
    // schema carries it per item; the last item of a result wins, as in C++.
    if let Some(experimental) = optional_value(node.get("experimentalMassToCharge")?)? {
        identification.mz = Some(experimental);
    }
    hit.metadata.insert(
        "pass_threshold".into(),
        MetaValue::new(MetaValueData::String(
            if pass { "true" } else { "false" }.into(),
        ))?,
    );
    attach_evidences(&mut hit, peptide_ref, run, library, options, budget, true)?;
    if identification.hits.len() >= options.max_list_items {
        return Err(bad("mzIdentML peptide hit count exceeds the limit"));
    }
    identification.hits.push(hit);
    Ok(())
}

/// Attach every `PeptideEvidence` of `peptide_ref` to `hit` and register the
/// proteins they point at on `run`.
///
/// Shared by the linear and the cross-linking read paths, which differ in one
/// respect: `positions_as_metadata` mirrors the linear path's `start`/`end`
/// metadata values on the hit, which the cross-linking path does not write.
///
/// # Errors
///
/// [`Error::Parse`] for a dangling `PeptideEvidence` or `dBSequence_ref` and
/// for an exceeded list ceiling, [`Error::InvalidRange`] for a `start`/`end`
/// pair outside the protein.
fn attach_evidences(
    hit: &mut PeptideHit,
    peptide_ref: &str,
    run: &mut ProteinIdentification,
    library: &Library,
    options: &ReadOptions,
    budget: &mut Budget,
    positions_as_metadata: bool,
) -> Result<()> {
    for id in library.by_peptide.get(peptide_ref).into_iter().flatten() {
        let evidence = library
            .evidences
            .get(id)
            .ok_or_else(|| bad(format!("unknown PeptideEvidence {id:?}")))?;
        budget.spend(16, 256)?;
        let mut item = PeptideEvidence::default();
        // Source: '-' means "not given" and leaves the OpenMS default.
        if let Some(pre) = evidence.pre.filter(|c| *c != '-') {
            item.aa_before = FlankingResidue::from_code(pre)?;
        }
        if let Some(post) = evidence.post.filter(|c| *c != '-') {
            item.aa_after = FlankingResidue::from_code(post)?;
        }
        if let (Some(start), Some(end)) = (evidence.start, evidence.end) {
            // mzIdentML counts residues from 1; OpenMS PeptideEvidence counts
            // from 0. The source stores the file's value unconverted, so its
            // positions are one too high and grow by one per store/load cycle
            // (see OpenMS_CPP_ISSUES.md).
            let convert = |value: i32, what: &str| -> Result<usize> {
                usize::try_from(value)
                    .ok()
                    .and_then(|value| value.checked_sub(1))
                    .ok_or_else(|| {
                        Error::InvalidRange(format!(
                            "PeptideEvidence {what} must be at least 1, got {value}"
                        ))
                    })
            };
            let start = convert(start, "start")?;
            let end = convert(end, "end")?;
            if start > end {
                return Err(Error::InvalidRange(
                    "PeptideEvidence start exceeds its end".into(),
                ));
            }
            if positions_as_metadata {
                hit.metadata.insert(
                    "start".into(),
                    MetaValue::new(MetaValueData::Integer(i64::try_from(start).unwrap_or(0)))?,
                );
                hit.metadata.insert(
                    "end".into(),
                    MetaValue::new(MetaValueData::Integer(i64::try_from(end).unwrap_or(0)))?,
                );
            }
            item.start = Some(start);
            item.end = Some(end);
        }
        merge_target_decoy(&mut hit.metadata, evidence.decoy)?;
        let sequence = library
            .db_sequences
            .get(&evidence.db_sequence_ref)
            .ok_or_else(|| {
                // Source: operator[] default-constructs, so the PSM gets an
                // empty protein accession and the run an empty ProteinHit.
                bad(format!(
                    "PeptideEvidence {id:?} references undeclared DBSequence {:?}",
                    evidence.db_sequence_ref
                ))
            })?;
        item.protein_accession = sequence.accession.clone();
        if run.find_hit(&sequence.accession).is_none() {
            if run.hits.len() >= options.max_list_items {
                return Err(bad("mzIdentML protein hit count exceeds the limit"));
            }
            let mut protein = ProteinHit {
                accession: sequence.accession.clone(),
                sequence: sequence.sequence.clone(),
                ..Default::default()
            };
            protein.metadata.insert(
                "isDecoy".into(),
                MetaValue::new(MetaValueData::String(
                    if evidence.decoy { "true" } else { "false" }.into(),
                ))?,
            );
            run.hits.push(protein);
        }
        if hit.evidences.len() >= options.max_list_items {
            return Err(bad("mzIdentML peptide evidence count exceeds the limit"));
        }
        // Every PeptideEvidence of the referenced Peptide is attached, not just
        // the ones this item's PeptideEvidenceRef children name; the source
        // calls those references redundant.
        hit.evidences.push(item);
    }
    Ok(())
}

/// The ion series each `Fragmentation` `cvParam` names, with its neutral loss.
const ION_TERMS: [(&str, char, &str); 18] = [
    ("MS:1001229", 'a', ""),
    ("MS:1001224", 'b', ""),
    ("MS:1001231", 'c', ""),
    ("MS:1001228", 'x', ""),
    ("MS:1001220", 'y', ""),
    ("MS:1001230", 'z', ""),
    ("MS:1001234", 'a', "-H2O"),
    ("MS:1001222", 'b', "-H2O"),
    ("MS:1001515", 'c', "-H2O"),
    ("MS:1001519", 'x', "-H2O"),
    ("MS:1001223", 'y', "-H2O"),
    ("MS:1001517", 'z', "-H2O"),
    ("MS:1001235", 'a', "-NH3"),
    ("MS:1001232", 'b', "-NH3"),
    ("MS:1001516", 'c', "-NH3"),
    ("MS:1001520", 'x', "-NH3"),
    ("MS:1001233", 'y', "-NH3"),
    ("MS:1001518", 'z', "-NH3"),
];

/// One `SpectrumIdentificationItem` of a cross-link group.
#[derive(Clone, Debug, Default)]
struct GroupItem {
    peptide_ref: String,
    experimental_mz: Option<f64>,
    rt: Option<f64>,
    user: BTreeMap<String, MetaValue>,
}

/// A `userParam` of an XL-MS item, typed as the source's XL path types it.
///
/// `parseSpectrumIdentificationItemSetXLMS` passes `unitName` to
/// `XMLHandler::fromXSDString`, while `parseUserParam_` - which the linear path
/// uses - passes `type`. The two source paths disagree, so this accepts either
/// spelling: `unitName` when it names an XSD type and no `type` attribute is
/// present, as the XL path does, and `type` otherwise, which is both the
/// schema's attribute and the one this module writes.
fn crosslink_user_value(node: &Node) -> Result<MetaValue> {
    let unit_name = attribute(node, "unitName");
    if unit_name.starts_with("xsd:") && attribute(node, "type").is_empty() {
        let mut typed = Node::new(&node.name);
        typed.attrs.insert("type".into(), unit_name.to_owned());
        if let Some(value) = node.optional("value") {
            typed.attrs.insert("value".into(), value.to_owned());
        }
        return user_value(&typed);
    }
    user_value(node)
}

/// The `Fragmentation` block of one item as [`PeakAnnotation`]s, as the
/// `IonType` loop of `parseSpectrumIdentificationItemSetXLMS`.
///
/// The annotation text is rebuilt in the source's shape,
/// `[<chain>|<category>$<series><index><loss>]`, which is what
/// `writeFragmentAnnotations_` parses back out. Where the source indexes its
/// per-index arrays without a bounds check - a `values` or `cross-link_chain`
/// list shorter than `index` - this substitutes an empty chain and a zero
/// intensity instead of reading out of bounds.
fn read_crosslink_fragmentation(
    item: &Node,
    options: &ReadOptions,
    budget: &mut Budget,
) -> Result<Vec<PeakAnnotation>> {
    let mut annotations = Vec::new();
    let Some(block) = gather(item, "Fragmentation").into_iter().next() else {
        return Ok(annotations);
    };
    for ion in gather(block, "IonType") {
        budget.spend(8, 128)?;
        let charge = integer(ion.get("charge")?, "IonType charge")?;
        let indices: Vec<&str> = attribute(ion, "index").split_whitespace().collect();
        let mut mz: Vec<&str> = Vec::new();
        let mut intensity: Vec<&str> = Vec::new();
        for array in gather(ion, "FragmentArray") {
            budget.spend(1, 0)?;
            let values = attribute(array, "values");
            match attribute(array, "measure_ref") {
                "Measure_mz" => mz = values.split_whitespace().collect(),
                "Measure_int" => intensity = values.split_whitespace().collect(),
                _ => {}
            }
        }
        let mut chains: Vec<&str> = Vec::new();
        let mut categories: Vec<&str> = Vec::new();
        for param in gather(ion, "userParam") {
            budget.spend(1, 0)?;
            let value = attribute(param, "value");
            match attribute(param, "name") {
                "cross-link_chain" => chains = value.split_whitespace().collect(),
                "cross-link_ioncategory" => categories = value.split_whitespace().collect(),
                _ => {}
            }
        }
        // The source tests every ion term in turn, so the last match wins and
        // an item without a recognised term keeps an empty series letter.
        let mut series = (' ', "");
        let mut named = false;
        for param in gather(ion, "cvParam") {
            budget.spend(1, 0)?;
            let accession = attribute(param, "accession");
            if let Some((_, kind, loss)) = ION_TERMS.iter().find(|(id, _, _)| *id == accession) {
                series = (*kind, *loss);
                named = true;
            }
        }
        for (position, index) in indices.iter().enumerate() {
            if annotations.len() >= options.max_list_items {
                return Err(bad("mzIdentML peak annotation count exceeds the limit"));
            }
            budget.spend(8, 128)?;
            let chain = chains.get(position).copied().unwrap_or_default();
            let category = categories.get(position).copied().unwrap_or_default();
            let (kind, loss) = series;
            let body = if named {
                format!("{kind}{index}{loss}")
            } else {
                (*index).to_owned()
            };
            annotations.push(PeakAnnotation {
                mz: match mz.get(position) {
                    Some(text) => finite(text)?,
                    None => 0.0,
                },
                intensity: match intensity.get(position) {
                    Some(text) => finite(text)?,
                    None => 0.0,
                },
                charge,
                annotation: format!("[{chain}|{category}${body}]"),
            });
        }
    }
    Ok(annotations)
}

/// One `SpectrumIdentificationResult` of a cross-linking document, as the
/// `xl_ms_search_` branch of `parseSpectrumIdentificationListElements_`.
///
/// The items of one result are grouped by the value of their `MS:1002511`
/// cvParam and each group becomes one [`PeptideIdentification`]: a labelled
/// cross-link reports a light and a heavy item per peptide, and a cross-link
/// reports an alpha and a beta peptide, so one result can carry four items that
/// describe a single match.
///
/// The source groups by an index that counts *every* element child of the
/// result and then looks the item up in the result's item list, so a result
/// with a non-item element child before its items groups the wrong item; this
/// groups by the item's own ordinal instead.
fn read_crosslink_result(
    result: &Node,
    run: &mut ProteinIdentification,
    library: &Library,
    context: &Context<'_>,
    budget: &mut Budget,
) -> Result<Vec<PeptideIdentification>> {
    let items: Vec<&Node> = children(result, "SpectrumIdentificationItem").collect();
    if items.is_empty() {
        // Source: it registers the group (0 -> item 0) unconditionally and then
        // dereferences item(0), which is null for a result without items.
        return Ok(Vec::new());
    }
    if items.len() > context.options.max_list_items {
        return Err(bad("mzIdentML item count per result exceeds the limit"));
    }
    let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (index, item) in items.iter().enumerate() {
        for node in gather(item, "cvParam") {
            budget.spend(1, 0)?;
            if attribute(node, "accession") != CROSSLINK_ITEM {
                continue;
            }
            let value = attribute(node, "value").to_owned();
            budget.text(&value)?;
            groups.entry(value).or_default().push(index);
        }
    }
    if groups.is_empty() {
        // Source comment: "fix for label-free mono-links - those only have one
        // SII and no 'cross-link spectrum identification item' value".
        groups.insert(String::new(), vec![0]);
    }
    let mut produced = Vec::new();
    for indices in groups.values() {
        if let Some(identification) =
            read_crosslink_group(result, &items, indices, run, library, context, budget)?
        {
            produced.push(identification);
        }
    }
    Ok(produced)
}

/// One cross-link group as a [`PeptideIdentification`], as
/// `parseSpectrumIdentificationItemSetXLMS`.
///
/// The alpha peptide is the one carrying the `crosslink donor` modification and
/// the beta peptide the `crosslink acceptor`; both become a [`PeptideHit`] of
/// one identification, tagged `xl_chain`, which the merge step of
/// `finish_crosslinks` later collapses into a single hit. A group without any
/// donor is read through the linear item path instead, exactly as the source
/// falls back.
///
/// # Errors
///
/// [`Error::Parse`] when an item of the group references an undeclared
/// `Peptide` or `PeptideEvidence`, or when a link value has no recorded
/// position: all of those are a `std::map::at` or `operator[]` in the source,
/// which either throws or invents an empty entry.
fn read_crosslink_group(
    result: &Node,
    items: &[&Node],
    indices: &[usize],
    run: &mut ProteinIdentification,
    library: &Library,
    context: &Context<'_>,
    budget: &mut Budget,
) -> Result<Option<PeptideIdentification>> {
    let options = context.options;
    let links = &library.cross_links;
    let mut group: Vec<GroupItem> = Vec::new();
    let mut annotations: Vec<PeakAnnotation> = Vec::new();
    let mut score = -1.0f64;
    let mut statistics: BTreeMap<&str, f64> = BTreeMap::new();
    let mut rank = 0i32;
    let mut charge = 0i32;
    for &index in indices {
        let item = *items
            .get(index)
            .ok_or_else(|| bad("internal cross-link item index out of range"))?;
        budget.spend(16, 512)?;
        let mut entry = GroupItem {
            peptide_ref: item.get("peptide_ref")?.to_owned(),
            experimental_mz: optional_value(item.get("experimentalMassToCharge")?)?,
            ..Default::default()
        };
        budget.text(&entry.peptide_ref)?;
        if rank == 0 {
            rank = integer(item.get("rank")?, "rank")?;
        }
        if charge == 0 {
            charge = integer(item.get("chargeState")?, "chargeState")?;
        }
        for node in gather(item, "cvParam") {
            // One item can belong to several groups - the source's own map
            // allows it - so every scan of its parameters is metered, not just
            // the item itself.
            budget.spend(1, 0)?;
            let param = CvParam::from_node(node);
            match param.accession.as_str() {
                // OpenXQuest:combined score, then OpenPepXL:score.
                "MS:1002681" | "MS:1003024" => score = score_value(&param.value)?,
                "MS:1002682" => {
                    statistics.insert("OpenPepXL:xcorr xlink", score_value(&param.value)?);
                }
                "MS:1002683" => {
                    statistics.insert("OpenPepXL:xcorr common", score_value(&param.value)?);
                }
                "MS:1002684" => {
                    statistics.insert("OpenPepXL:match-odds", score_value(&param.value)?);
                }
                "MS:1002685" => {
                    statistics.insert("OpenPepXL:intsum", score_value(&param.value)?);
                }
                "MS:1002686" => {
                    statistics.insert("OpenPepXL:wTIC", score_value(&param.value)?);
                }
                "MS:1000894" => {
                    let mut rt = finite(&param.value)?;
                    if param.unit_accession == "UO:0000031" {
                        rt *= 60.0;
                    }
                    entry.rt = Some(rt);
                }
                _ => {}
            }
        }
        for node in gather(item, "userParam") {
            budget.spend(1, 0)?;
            let name = attribute(node, "name").to_owned();
            budget.text(&name)?;
            if entry.user.len() >= options.max_list_items {
                return Err(bad("mzIdentML parameter group limit exceeded"));
            }
            // std::map::insert keeps the first entry for a repeated name.
            let value = crosslink_user_value(node)?;
            entry.user.entry(name).or_insert(value);
        }
        if annotations.is_empty() {
            annotations = read_crosslink_fragmentation(item, options, budget)?;
        }
        group.push(entry);
    }
    let light_mz = group
        .iter()
        .filter_map(|item| item.experimental_mz)
        .fold(f64::INFINITY, f64::min);
    let heavy_mz = group
        .iter()
        .filter_map(|item| item.experimental_mz)
        .fold(f64::NEG_INFINITY, f64::max);
    if !light_mz.is_finite() {
        // Source: min_element over a vector of NaN, then light[0] on an empty
        // index vector. No item of the group declares an experimental m/z, so
        // there is no spectrum position to attach the match to.
        return Ok(None);
    }
    let spectrum_id = result.get("spectrumID")?;
    let references: Vec<&str> = spectrum_id.split(',').collect();
    // A labelled cross-link is reported as a light and a heavy item whose
    // experimental m/z differ AND that reference two spectra.
    let labeled = light_mz != heavy_mz && references.len() > 1;
    let light: Vec<usize> = (0..group.len())
        .filter(|&index| group[index].experimental_mz == Some(light_mz))
        .collect();
    let heavy: Vec<usize> = (0..group.len())
        .filter(|&index| group[index].experimental_mz != Some(light_mz))
        .collect();
    let first_light = *light
        .first()
        .ok_or_else(|| bad("internal cross-link light index is empty"))?;
    let light_rt = group[first_light].rt;
    let heavy_rt = match (labeled, heavy.first()) {
        (true, Some(&index)) => group[index].rt,
        _ => light_rt,
    };
    let alpha: Vec<usize> = (0..group.len())
        .filter(|&index| links.donor.contains_key(&group[index].peptide_ref))
        .collect();
    let beta: Vec<usize> = (0..group.len())
        .filter(|&index| !links.donor.contains_key(&group[index].peptide_ref))
        .collect();
    let Some(&first_alpha) = alpha.first() else {
        // Source: "No crosslink donor found for SIIs in ...; parsing as regular
        // (non-XLMS) peptide identifications."
        let mut identification = PeptideIdentification {
            identifier: run.identifier.clone(),
            higher_score_better: false,
            rt: group.first().and_then(|item| item.rt),
            mz: group.first().and_then(|item| item.experimental_mz),
            ..Default::default()
        };
        identification.set_spectrum_reference(spectrum_id);
        for item in items {
            read_item(item, &mut identification, run, library, context, budget)?;
        }
        identification.sort()?;
        return Ok(Some(identification));
    };
    let alpha_peptide = group[first_alpha].peptide_ref.clone();
    let donor = links
        .donor
        .get(&alpha_peptide)
        .ok_or_else(|| bad("internal cross-link donor lookup failed"))?;
    let mut alpha_position = *links.donor_position.get(donor).ok_or_else(|| {
        bad(format!(
            "cross-link donor {donor:?} has no recorded modification position"
        ))
    })?;
    let acceptor_of_alpha = links.acceptor.get(&alpha_peptide);
    let kind = if alpha.len() == beta.len() {
        // A beta peptide at all means a cross-link between two peptides.
        "cross-link"
    } else if acceptor_of_alpha == Some(donor) {
        // Donor and acceptor of the same link on one peptide is a loop-link.
        "loop-link"
    } else {
        "mono-link"
    };
    let mut identification = PeptideIdentification {
        identifier: run.identifier.clone(),
        score_type: user_param::OPENPEPXL_SCORE.to_owned(),
        higher_score_better: true,
        rt: light_rt,
        mz: Some(light_mz),
        ..Default::default()
    };
    identification.set_spectrum_reference(spectrum_id);
    let sequence = library
        .peptides
        .get(&alpha_peptide)
        .cloned()
        .ok_or_else(|| {
            bad(format!(
                "SpectrumIdentificationItem references undeclared Peptide {alpha_peptide:?}"
            ))
        })?;
    let mut hit = PeptideHit {
        sequence,
        score,
        rank: u32::try_from(rank.max(1) - 1)
            .map_err(|_| invalid("SpectrumIdentificationItem rank exceeds the rank range"))?,
        charge,
        ..Default::default()
    };
    let first_reference = references.first().copied().unwrap_or_default().to_owned();
    hit.metadata.insert(
        user_param::SPECTRUM_REFERENCE.into(),
        first_reference.clone().into(),
    );
    hit.metadata
        .insert("xl_chain".into(), CROSSLINK_DONOR.to_owned().into());
    if labeled {
        if let Some(rt) = heavy_rt {
            hit.metadata.insert(
                user_param::OPENPEPXL_HEAVY_SPEC_RT.into(),
                MetaValue::new(MetaValueData::Float(rt))?,
            );
        }
        hit.metadata.insert(
            user_param::OPENPEPXL_HEAVY_SPEC_MZ.into(),
            MetaValue::new(MetaValueData::Float(heavy_mz))?,
        );
        if let Some(reference) = references.get(1) {
            hit.metadata.insert(
                user_param::OPENPEPXL_HEAVY_SPEC_REF.into(),
                (*reference).to_owned().into(),
            );
        }
    }
    hit.metadata
        .insert(user_param::OPENPEPXL_XL_TYPE.into(), kind.to_owned().into());
    hit.metadata
        .insert(user_param::OPENPEPXL_XL_RANK.into(), rank.into());
    for name in [
        "OpenPepXL:xcorr xlink",
        "OpenPepXL:xcorr common",
        "OpenPepXL:match-odds",
        "OpenPepXL:intsum",
        "OpenPepXL:wTIC",
    ] {
        let value = statistics.get(name).copied().unwrap_or(0.0);
        hit.metadata.insert(
            name.to_owned(),
            MetaValue::new(MetaValueData::Float(value))?,
        );
    }
    for (name, value) in &group[first_alpha].user {
        hit.metadata.insert(name.clone(), value.clone());
    }
    hit.peak_annotations = annotations;
    if kind == "loop-link" {
        let acceptor =
            acceptor_of_alpha.ok_or_else(|| bad("internal cross-link acceptor lookup failed"))?;
        let position = *links.acceptor_position.get(acceptor).ok_or_else(|| {
            bad(format!(
                "cross-link acceptor {acceptor:?} has no recorded modification position"
            ))
        })?;
        hit.metadata
            .insert(user_param::OPENPEPXL_XL_POS2.into(), position.into());
    }
    if kind != "mono-link" {
        let name = links.modification.get(&alpha_peptide).ok_or_else(|| {
            bad(format!(
                "cross-linked Peptide {alpha_peptide:?} declares no cross-linker name"
            ))
        })?;
        let mass = *links.mass.get(&alpha_peptide).ok_or_else(|| {
            bad(format!(
                "cross-linked Peptide {alpha_peptide:?} declares no cross-linker mass"
            ))
        })?;
        hit.metadata
            .insert(user_param::OPENPEPXL_XL_MOD.into(), name.clone().into());
        hit.metadata.insert(
            user_param::OPENPEPXL_XL_MASS.into(),
            MetaValue::new(MetaValueData::Float(mass))?,
        );
    } else if let Some(name) = links.modification.get(&alpha_peptide) {
        hit.metadata
            .insert(user_param::OPENPEPXL_XL_MOD.into(), name.clone().into());
    }
    let length = i32::try_from(hit.sequence.len())
        .map_err(|_| invalid("peptide is too long for a cross-link position"))?;
    // A link one residue past the sequence is the C-terminus and one before its
    // start the N-terminus; both are reported at the adjacent residue.
    let specificity = if alpha_position == -1 {
        alpha_position += 1;
        "N_TERM"
    } else if alpha_position == length {
        alpha_position -= 1;
        "C_TERM"
    } else {
        "ANYWHERE"
    };
    hit.metadata
        .insert(user_param::OPENPEPXL_XL_POS1.into(), alpha_position.into());
    hit.metadata.insert(
        user_param::OPENPEPXL_XL_TERM_SPEC_ALPHA.into(),
        specificity.to_owned().into(),
    );
    let mut hits = Vec::new();
    let mut peptides = vec![alpha_peptide];
    if kind == "cross-link" {
        let beta_peptide = group[*beta
            .first()
            .ok_or_else(|| bad("internal cross-link beta index is empty"))?]
        .peptide_ref
        .clone();
        let acceptor = links.acceptor.get(&beta_peptide).ok_or_else(|| {
            bad(format!(
                "cross-linked Peptide {beta_peptide:?} declares no crosslink acceptor"
            ))
        })?;
        let mut beta_position = *links.acceptor_position.get(acceptor).ok_or_else(|| {
            bad(format!(
                "cross-link acceptor {acceptor:?} has no recorded modification position"
            ))
        })?;
        let sequence = library
            .peptides
            .get(&beta_peptide)
            .cloned()
            .ok_or_else(|| {
                bad(format!(
                    "SpectrumIdentificationItem references undeclared Peptide {beta_peptide:?}"
                ))
            })?;
        let length = i32::try_from(sequence.len())
            .map_err(|_| invalid("peptide is too long for a cross-link position"))?;
        hit.metadata.insert(
            user_param::OPENPEPXL_BETA_SEQUENCE.into(),
            sequence.to_string().into(),
        );
        let mut beta_hit = PeptideHit {
            sequence,
            score,
            rank: hit.rank,
            charge,
            ..Default::default()
        };
        beta_hit.metadata.insert(
            user_param::SPECTRUM_REFERENCE.into(),
            first_reference.into(),
        );
        beta_hit
            .metadata
            .insert("xl_chain".into(), CROSSLINK_ACCEPTOR.to_owned().into());
        let specificity = if beta_position == -1 {
            beta_position += 1;
            "N_TERM"
        } else if beta_position == length {
            beta_position -= 1;
            "C_TERM"
        } else {
            "ANYWHERE"
        };
        hit.metadata
            .insert(user_param::OPENPEPXL_XL_POS2.into(), beta_position.into());
        hit.metadata.insert(
            user_param::OPENPEPXL_XL_TERM_SPEC_BETA.into(),
            specificity.to_owned().into(),
        );
        hits.push(hit);
        hits.push(beta_hit);
        peptides.push(beta_peptide);
    } else {
        hit.metadata.insert(
            user_param::OPENPEPXL_XL_TERM_SPEC_BETA.into(),
            "ANYWHERE".to_owned().into(),
        );
        if kind != "loop-link" {
            // The source overwrites the loop-link's second position with the
            // mono-link placeholder here, which makes its own loop-link branch
            // dead and drops the second half of every loop-link on a store.
            hit.metadata
                .insert(user_param::OPENPEPXL_XL_POS2.into(), "-".to_owned().into());
        }
        hit.metadata.insert(
            user_param::OPENPEPXL_BETA_SEQUENCE.into(),
            "-".to_owned().into(),
        );
        hits.push(hit);
    }
    for (position, peptide_ref) in peptides.iter().enumerate() {
        let Some(target) = hits.get_mut(position) else {
            continue;
        };
        // The XL path does not mirror the evidence positions into hit metadata.
        attach_evidences(target, peptide_ref, run, library, options, budget, false)?;
    }
    identification.hits = hits;
    identification.sort()?;
    Ok(Some(identification))
}

/// Apply a `SpectrumIdentificationResult`'s own parameters to one
/// identification: the retention time becomes a coordinate, everything else
/// metadata.
fn apply_result_params(
    identification: &mut PeptideIdentification,
    group: &ParamGroup,
) -> Result<()> {
    for (accession, params) in &group.cv {
        let Some(param) = params.first() else {
            continue;
        };
        if accession == "MS:1000894" || accession == "MS:1000016" {
            let mut rt = finite(&param.value)?;
            if param.unit_accession == "UO:0000031" {
                rt *= 60.0;
            }
            identification.rt = Some(rt);
        } else {
            identification.metadata.insert(
                accession.clone(),
                MetaValue::new(MetaValueData::String(param.value.clone()))?,
            );
        }
    }
    for (name, value) in &group.user {
        identification.metadata.insert(name.clone(), value.clone());
    }
    Ok(())
}

fn read_lists(
    root: &Node,
    runs: &mut [ProteinIdentification],
    list_to_run: &BTreeMap<String, usize>,
    library: &Library,
    context: &Context<'_>,
    budget: &mut Budget,
) -> Result<Vec<PeptideIdentification>> {
    let options = context.options;
    let lists = gather(root, "SpectrumIdentificationList");
    if lists.is_empty() {
        return Err(missing(
            "mzIdentML has no SpectrumIdentificationList element",
        ));
    }
    let mut identifications = Vec::new();
    let mut seen = BTreeSet::new();
    for list in lists {
        let id = element_id(list)?;
        if !seen.insert(id.to_owned()) {
            return Err(bad(format!(
                "duplicate SpectrumIdentificationList id {id:?}"
            )));
        }
        let index = *list_to_run.get(id).ok_or_else(|| {
            // Source: operator[] yields index 0, so an unreferenced list's
            // results are silently attributed to the first run.
            bad(format!(
                "SpectrumIdentificationList {id:?} is not referenced by any SpectrumIdentification"
            ))
        })?;
        let run = runs
            .get_mut(index)
            .ok_or_else(|| bad("internal run index out of range"))?;
        for result in children(list, "SpectrumIdentificationResult") {
            budget.spend(32, 1024)?;
            let produced = if context.crosslinking {
                read_crosslink_result(result, run, library, context, budget)?
            } else {
                let mut identification = PeptideIdentification {
                    identifier: run.identifier.clone(),
                    // Either a q-value or an e-value unless a specific score wins.
                    higher_score_better: false,
                    ..Default::default()
                };
                identification.set_spectrum_reference(result.get("spectrumID")?);
                for item in children(result, "SpectrumIdentificationItem") {
                    read_item(item, &mut identification, run, library, context, budget)?;
                }
                identification.sort()?;
                vec![identification]
            };
            let group = param_group(result, options, budget)?;
            for mut identification in produced {
                // The source applies the result's own parameters to pep_id_
                // ->back(), so with several cross-link groups per result only
                // the last one receives the result's retention time - and with
                // no group at all they land on an unrelated, earlier
                // identification. Every group of the result gets them here.
                apply_result_params(&mut identification, &group)?;
                if identifications.len() >= options.max_list_items {
                    return Err(bad(
                        "mzIdentML peptide identification count exceeds the limit",
                    ));
                }
                identifications.push(identification);
            }
        }
    }
    Ok(identifications)
}

/// `ProteinDetectionList` to protein hits, as `parseProteinDetectionListElements_`.
///
/// Every `ProteinDetectionHypothesis` becomes a [`ProteinHit`] on the **last**
/// run, whichever run the list belongs to, because the source appends to
/// `pro_id_->back()`. Reproduced rather than corrected: the element carries no
/// reference that would identify the right run, so any other choice would be
/// invented. A hypothesis is appended even when the run already has that
/// accession, so duplicates are possible - also as in the source.
fn read_protein_detection(
    root: &Node,
    runs: &mut [ProteinIdentification],
    library: &Library,
    options: &ReadOptions,
    budget: &mut Budget,
) -> Result<()> {
    let lists = gather(root, "ProteinDetectionList");
    if lists.is_empty() {
        return Ok(());
    }
    let run = runs
        .last_mut()
        .ok_or_else(|| missing("ProteinDetectionList without any identification run"))?;
    for list in lists {
        for group in children(list, "ProteinAmbiguityGroup") {
            for hypothesis in children(group, "ProteinDetectionHypothesis") {
                budget.spend(16, 256)?;
                let reference = hypothesis.get("dBSequence_ref")?;
                let sequence = library.db_sequences.get(reference).ok_or_else(|| {
                    bad(format!(
                        "ProteinDetectionHypothesis references undeclared DBSequence {reference:?}"
                    ))
                })?;
                if run.hits.len() >= options.max_list_items {
                    return Err(bad("mzIdentML protein hit count exceeds the limit"));
                }
                run.hits.push(ProteinHit {
                    accession: sequence.accession.clone(),
                    sequence: sequence.sequence.clone(),
                    ..Default::default()
                });
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Cross-linking post-processing (ANALYSIS/XLMS/OPXLHelper.h)
// ---------------------------------------------------------------------------

/// The text of a metadata value, or the empty string when it is absent.
fn metadata_text(metadata: &MetaInfo, key: &str) -> String {
    metadata.get(key).map(text_of).unwrap_or_default()
}

/// The cross-link position stored at `key`, or `None` when it is absent or the
/// `"-"` placeholder a mono-link carries.
fn crosslink_position(metadata: &MetaInfo, key: &str) -> Result<Option<i64>> {
    match metadata.get(key) {
        None => Ok(None),
        Some(value) => match value.data() {
            MetaValueData::Integer(position) => Ok(Some(*position)),
            MetaValueData::String(text) if text == "-" => Ok(None),
            _ => Err(invalid(format!(
                "cross-link position {key} is neither an integer nor \"-\""
            ))),
        },
    }
}

/// Protein-coordinate link positions, as `OPXLHelper::addProteinPositionMetaValues`.
///
/// The residue number is `evidence start + link position + 1`, so the
/// N-terminal residue of a protein is number 1. The source comment states that
/// both inputs are zero-based, which is true of this reader; in C++ the
/// mzIdentML reader stores the file's one-based `start` unconverted, so the
/// positions it derives here are one too high.
///
/// # Errors
///
/// [`Error::InvalidValue`] when a link position is neither an integer nor the
/// mono-link placeholder, and [`Error::InvalidRange`] when a protein
/// coordinate overflows.
fn add_protein_position_meta_values(identifications: &mut [PeptideIdentification]) -> Result<()> {
    fn positions(evidences: &[PeptideEvidence], link: i64) -> Result<String> {
        let mut out = String::new();
        for evidence in evidences {
            // An unknown start is the source's UNKNOWN_POSITION, -1.
            let start = match evidence.start {
                Some(start) => i64::try_from(start)
                    .map_err(|_| Error::InvalidRange("peptide start overflows".into()))?,
                None => -1,
            };
            let position = start
                .checked_add(link)
                .and_then(|value| value.checked_add(1))
                .ok_or_else(|| {
                    Error::InvalidRange("cross-link protein position overflows".into())
                })?;
            if !out.is_empty() {
                out.push(',');
            }
            out.push_str(&position.to_string());
        }
        Ok(out)
    }
    for identification in identifications {
        let Some(alpha) = identification.hits.first() else {
            continue;
        };
        // A non-cross-linked identification read through the fallback path has
        // no link position and is skipped.
        let Some(position1) = crosslink_position(&alpha.metadata, user_param::OPENPEPXL_XL_POS1)?
        else {
            continue;
        };
        let position2 = crosslink_position(&alpha.metadata, user_param::OPENPEPXL_XL_POS2)?;
        let first = positions(&alpha.evidences, position1)?;
        let second = match identification.hits.len() {
            2 => {
                let beta = identification
                    .hits
                    .get(1)
                    .ok_or_else(|| invalid("internal cross-link hit index out of range"))?;
                let link = position2.ok_or_else(|| {
                    invalid("a cross-link with a beta hit needs a second link position")
                })?;
                Some(positions(&beta.evidences, link)?)
            }
            _ => match position2 {
                Some(link) => Some(positions(&alpha.evidences, link)?),
                None => Some("-".to_owned()),
            },
        };
        let two_chains = identification.hits.len() == 2;
        for (index, hit) in identification.hits.iter_mut().enumerate() {
            if index == 0 || two_chains {
                hit.metadata.insert(
                    user_param::OPENPEPXL_XL_POS1_PROT.into(),
                    first.clone().into(),
                );
            }
            if let Some(second) = &second {
                if index == 0 || two_chains {
                    hit.metadata.insert(
                        user_param::OPENPEPXL_XL_POS2_PROT.into(),
                        second.clone().into(),
                    );
                }
            }
        }
    }
    Ok(())
}

/// The beta chain's protein accessions, as `OPXLHelper::addBetaAccessions`.
fn add_beta_accessions(identifications: &mut [PeptideIdentification]) -> Result<()> {
    for identification in identifications {
        if identification.hits.is_empty() {
            continue;
        }
        let accessions = match identification.hits.len() {
            2 => {
                let beta = identification
                    .hits
                    .get(1)
                    .ok_or_else(|| invalid("internal cross-link hit index out of range"))?;
                beta.evidences
                    .iter()
                    .map(|evidence| evidence.protein_accession.clone())
                    .collect::<Vec<_>>()
                    .join(",")
            }
            _ => "-".to_owned(),
        };
        let two_chains = identification.hits.len() == 2;
        for (index, hit) in identification.hits.iter_mut().enumerate() {
            if index == 0 || two_chains {
                hit.metadata.insert(
                    user_param::OPENPEPXL_BETA_ACCESSIONS.into(),
                    accessions.clone().into(),
                );
            }
        }
    }
    Ok(())
}

/// Per-chain target/decoy state, as `OPXLHelper::addXLTargetDecoyMV`.
///
/// A cross-link counts as a decoy as soon as one of its two chains has no
/// target evidence at all.
fn add_crosslink_target_decoy(identifications: &mut [PeptideIdentification]) -> Result<()> {
    for identification in identifications {
        if identification.hits.is_empty() {
            continue;
        }
        let alpha_state = identification
            .hits
            .first()
            .map(|hit| metadata_text(&hit.metadata, user_param::TARGET_DECOY))
            .unwrap_or_default();
        let beta_state = match identification.hits.len() {
            2 => identification
                .hits
                .get(1)
                .map(|hit| metadata_text(&hit.metadata, user_param::TARGET_DECOY)),
            _ => None,
        };
        let Some(alpha) = identification.hits.first_mut() else {
            continue;
        };
        alpha.metadata.insert(
            user_param::OPENPEPXL_TARGET_DECOY_ALPHA.into(),
            alpha_state.clone().into(),
        );
        match &beta_state {
            Some(beta_state) => {
                alpha.metadata.insert(
                    user_param::OPENPEPXL_TARGET_DECOY_BETA.into(),
                    beta_state.clone().into(),
                );
                if !alpha_state.contains("target") || !beta_state.contains("target") {
                    alpha
                        .metadata
                        .insert(user_param::TARGET_DECOY.into(), "decoy".to_owned().into());
                }
            }
            None => {
                alpha.metadata.insert(
                    user_param::OPENPEPXL_TARGET_DECOY_BETA.into(),
                    "-".to_owned().into(),
                );
            }
        }
    }
    Ok(())
}

/// Fold the beta hit into the alpha hit and collect the identifications of one
/// spectrum, as `OPXLHelper::removeBetaPeptideHits`.
///
/// The beta chain has no `PeptideHit` of its own afterwards: its evidences live
/// on in the `BetaPepEv:` metadata values, which is what the writer turns back
/// into `PeptideEvidence` elements. Identifications that share a spectrum
/// reference become one identification with several hits, as in the source,
/// which is why an upstream section counts one hit per cross-link.
///
/// Two divergences, both preserving data the source drops:
///
/// * The source keeps only the *first* hit of every input identification, so a
///   spectrum read through the non-cross-linked fallback path - a noncovalent
///   association, for instance - loses every candidate but the best one. Only
///   a folded beta chain is removed here.
/// * The merged identification keeps its run identifier, which the source
///   drops, leaving every cross-linking PSM unlinked from its protein run.
///
/// `BetaPepEv:start` and `:end` use decimal positions in both implementations:
/// the source global numeric `std::string += Int` overload in `StringUtils.h`
/// calls `appendToStr`, which formats decimal digits with `std::to_chars`.
fn merge_crosslink_hits(
    identifications: Vec<PeptideIdentification>,
    options: &ReadOptions,
) -> Result<Vec<PeptideIdentification>> {
    let mut merged: BTreeMap<String, PeptideIdentification> = BTreeMap::new();
    for identification in identifications {
        let two_chains = identification.hits.len() == 2
            && identification
                .hits
                .get(1)
                .map(|hit| metadata_text(&hit.metadata, "xl_chain"))
                .as_deref()
                == Some(CROSSLINK_ACCEPTOR);
        let beta = if two_chains {
            identification.hits.get(1).cloned()
        } else {
            None
        };
        let reference = identification.spectrum_reference();
        let mut rest: Vec<PeptideHit> = Vec::new();
        let mut hits = identification.hits.into_iter();
        let Some(mut hit) = hits.next() else {
            continue;
        };
        if two_chains {
            // The folded beta chain is the only hit that goes away.
            let _ = hits.next();
        }
        rest.extend(hits);
        hit.metadata.remove("xl_chain");
        if let Some(beta) = beta {
            let mut pre = Vec::new();
            let mut post = Vec::new();
            let mut start = Vec::new();
            let mut end = Vec::new();
            for evidence in &beta.evidences {
                pre.push(evidence.aa_before.code().to_string());
                post.push(evidence.aa_after.code().to_string());
                start.push(match evidence.start {
                    Some(value) => value.to_string(),
                    None => "-1".to_owned(),
                });
                end.push(match evidence.end {
                    Some(value) => value.to_string(),
                    None => "-1".to_owned(),
                });
            }
            for (key, value) in [
                (user_param::OPENPEPXL_BETA_PEPEV_PRE, pre.join(",")),
                (user_param::OPENPEPXL_BETA_PEPEV_POST, post.join(",")),
                (user_param::OPENPEPXL_BETA_PEPEV_START, start.join(",")),
                (user_param::OPENPEPXL_BETA_PEPEV_END, end.join(",")),
            ] {
                hit.metadata.insert(key.into(), value.into());
            }
        }
        let entry = merged.entry(reference.clone()).or_insert_with(|| {
            let mut fresh = PeptideIdentification {
                identifier: identification.identifier.clone(),
                score_type: user_param::OPENPEPXL_SCORE.to_owned(),
                rt: identification.rt,
                mz: identification.mz,
                ..Default::default()
            };
            fresh.set_spectrum_reference(reference.clone());
            fresh
        });
        for mut hit in std::iter::once(hit).chain(rest) {
            if entry.hits.len() >= options.max_list_items {
                return Err(bad("mzIdentML peptide hit count exceeds the limit"));
            }
            hit.metadata.remove("xl_chain");
            entry.hits.push(hit);
        }
    }
    Ok(merged.into_values().collect())
}

/// Score ratios between neighbouring ranks, as `OPXLHelper::computeDeltaScores`.
///
/// The last rank, and a single hit, get a delta score of zero. A ratio the
/// source would produce as infinity or NaN - the next rank's score divided by a
/// zero score - is left unset instead, because every stored float is finite.
fn compute_delta_scores(identifications: &mut [PeptideIdentification]) -> Result<()> {
    for identification in identifications {
        identification.sort()?;
        let scores: Vec<f64> = identification.hits.iter().map(|hit| hit.score).collect();
        let last = scores.len().saturating_sub(1);
        for (index, hit) in identification.hits.iter_mut().enumerate() {
            let delta = if index == last {
                0.0
            } else {
                match scores.get(index + 1) {
                    Some(next) => next / scores[index],
                    None => continue,
                }
            };
            if !delta.is_finite() {
                continue;
            }
            hit.metadata.insert(
                user_param::DELTA_SCORE.into(),
                MetaValue::new(MetaValueData::Float(delta))?,
            );
        }
    }
    Ok(())
}

/// The Percolator feature list, as `OPXLHelper::addPercolatorFeatureList`.
fn add_percolator_features(run: &mut ProteinIdentification) -> Result<()> {
    run.search_parameters.metadata.insert(
        "feature_extractor".into(),
        "TOPP_PSMFeatureExtractor".to_owned().into(),
    );
    run.search_parameters.metadata.insert(
        "extra_features".into(),
        PERCOLATOR_FEATURES.join(",").into(),
    );
    Ok(())
}

/// The six `OPXLHelper` steps `readMzIdentMLFile` runs on a cross-linking
/// document, in the source's order.
fn finish_crosslinks(
    runs: &mut [ProteinIdentification],
    identifications: Vec<PeptideIdentification>,
    options: &ReadOptions,
) -> Result<Vec<PeptideIdentification>> {
    let mut identifications = identifications;
    add_protein_position_meta_values(&mut identifications)?;
    add_beta_accessions(&mut identifications)?;
    add_crosslink_target_decoy(&mut identifications)?;
    let mut identifications = merge_crosslink_hits(identifications, options)?;
    compute_delta_scores(&mut identifications)?;
    if let Some(run) = runs.first_mut() {
        add_percolator_features(run)?;
    }
    Ok(identifications)
}

/// Read mzIdentML with default limits and the global modification registry.
///
/// # Errors
///
/// See [`read_with_registry`].
pub fn read(reader: impl BufRead) -> Result<MzIdentMLDocument> {
    read_with_options(reader, &ReadOptions::default())
}

/// Read mzIdentML with explicit limits; see
/// [`read_with_registry`].
///
/// # Errors
///
/// See [`read_with_registry`].
pub fn read_with_options(reader: impl BufRead, options: &ReadOptions) -> Result<MzIdentMLDocument> {
    read_with_registry(reader, options, ModificationsDB::global())
}

/// Read one mzIdentML document, resolving modifications against `registry`.
///
/// Returns an owned document only after the whole input has parsed and every
/// reference has resolved, so a failure leaves the caller's state untouched.
/// One `SpectrumIdentification` becomes one [`ProteinIdentification`] and one
/// `SpectrumIdentificationResult` one [`PeptideIdentification`], both in
/// document order, linked through
/// [`PeptideIdentification::identifier`](crate::identification::PeptideIdentification::identifier).
///
/// A document that declares `MS:1002494` is read through the cross-linking
/// path instead: the items of one result are grouped by their cross-link, the
/// two chains of a match become one identification, and the six `OPXLHelper`
/// post-processing steps then collapse the pairs per spectrum, so one
/// identification per spectrum reference comes back rather than one per
/// result.
///
/// # Errors
///
/// * [`Error::Unsupported`] for a non-UTF-8 document, a foreign namespace, a
///   `DOCTYPE` declaration and an external entity reference.
/// * [`Error::MissingInformation`] when `SpectraData`, `SpectrumIdentification`,
///   `SpectrumIdentificationProtocol` or `SpectrumIdentificationList` is absent.
///   The source throws `std::runtime_error` for the same four.
/// * [`Error::Parse`] for malformed XML, a duplicate element id, a dangling
///   `peptide_ref` or `dBSequence_ref`, a `SpectrumIdentificationList` no run
///   references, and for every exceeded ceiling in [`ReadOptions`].
/// * [`Error::InvalidRange`] for a `Modification` or `SubstitutionModification`
///   `location`, or a `PeptideEvidence` `start`/`end`, outside the peptide.
/// * [`Error::InvalidValue`] when a modification does not resolve against
///   `registry`, and [`Error::Io`] on a read failure.
/// * [`Error::Parse`] also for a cross-link whose link value has no recorded
///   modification position, and for a cross-linked `Peptide` that declares no
///   cross-linker name or mass: the source reads all three with
///   `std::map::at`, which throws.
///
/// Which dangling references are tolerated and which are refused is the whole
/// table in `docs/MZIDENTML_SUPPORT.md`; the short rule is that a reference the
/// source would silently replace with invented data (a peptide with no
/// sequence, a protein with no accession) is an error here, while one it merely
/// leaves unresolved as metadata (`searchDatabase_ref`, `spectraData_ref`,
/// `spectrumIdentificationProtocol_ref`) stays tolerated, because upstream test
/// data contains exactly those.
pub fn read_with_registry(
    reader: impl BufRead,
    options: &ReadOptions,
    registry: &ModificationsDB,
) -> Result<MzIdentMLDocument> {
    let mut budget = Budget {
        work: options.max_work,
        bytes: options.max_payload_bytes,
    };
    let root = parse(reader, options, &mut budget)?;
    let cv = ControlledVocabulary::psi_ms()?;
    let context = Context {
        cv,
        registry,
        options,
        terms: ScoreTerms::new(cv)?,
        crosslinking: is_crosslinking(&root),
    };
    let mut library = Library::default();
    read_inputs(&root, &mut library, &mut budget)?;
    if library.spectra_data.is_empty() {
        return Err(missing("mzIdentML has no SpectraData element"));
    }
    read_software(&root, &mut library, &context, &mut budget)?;
    let (mut runs, links, list_to_run) =
        read_runs(&root, &library, options, context.crosslinking, &mut budget)?;
    read_protocols(&root, &library, &mut runs, &links, &context, &mut budget)?;
    read_sequence_collection(&root, &mut library, &context, &mut budget)?;
    let peptide_identifications = read_lists(
        &root,
        &mut runs,
        &list_to_run,
        &library,
        &context,
        &mut budget,
    )?;
    read_protein_detection(&root, &mut runs, &library, options, &mut budget)?;
    for run in &mut runs {
        sort_protein_hits(run);
    }
    let peptide_identifications = if context.crosslinking {
        finish_crosslinks(&mut runs, peptide_identifications, options)?
    } else {
        peptide_identifications
    };
    for run in &runs {
        run.validate()?;
    }
    for identification in &peptide_identifications {
        identification.validate()?;
    }
    Ok(MzIdentMLDocument {
        document_id: attribute(&root, "id").to_owned(),
        protein_identifications: runs,
        peptide_identifications,
    })
}

/// Load plain or magic-detected gzip/bzip2 mzIdentML with default limits.
///
/// # Errors
///
/// See [`read_with_registry`].
pub fn load(path: impl AsRef<Path>) -> Result<MzIdentMLDocument> {
    load_with_options(path, &ReadOptions::default())
}

/// Load with explicit limits and the global modification registry.
///
/// # Errors
///
/// See [`read_with_registry`].
pub fn load_with_options(
    path: impl AsRef<Path>,
    options: &ReadOptions,
) -> Result<MzIdentMLDocument> {
    load_with_registry(path, options, ModificationsDB::global())
}

/// Load, resolving modifications against `registry`.
///
/// # Errors
///
/// See [`read_with_registry`].
pub fn load_with_registry(
    path: impl AsRef<Path>,
    options: &ReadOptions,
    registry: &ModificationsDB,
) -> Result<MzIdentMLDocument> {
    read_with_registry(super::path_io::open(path.as_ref())?, options, registry)
}

/// Replace `destination` only after the complete load succeeds.
///
/// `MzIdentMLFile::load` clears both output containers first, because its DOM
/// handler only appends; this replaces them atomically instead, so a failed
/// load leaves the previous contents intact.
///
/// # Errors
///
/// See [`read_with_registry`].
pub fn load_into(path: impl AsRef<Path>, destination: &mut MzIdentMLDocument) -> Result<()> {
    let document = load(path)?;
    *destination = document;
    Ok(())
}
// ---------------------------------------------------------------------------
// Write
// ---------------------------------------------------------------------------

/// Bounded output buffer.
struct Out {
    text: String,
    max_bytes: usize,
    records: usize,
    max_records: usize,
}
impl Out {
    fn new(options: &WriteOptions) -> Self {
        Self {
            text: String::new(),
            max_bytes: options.max_output_bytes,
            records: 0,
            max_records: options.max_records,
        }
    }
    fn raw(&mut self, value: &str) -> Result<()> {
        if self.text.len().saturating_add(value.len()) > self.max_bytes {
            return Err(bad("mzIdentML output byte limit exceeded"));
        }
        self.text.push_str(value);
        Ok(())
    }
    fn element(&mut self) -> Result<()> {
        self.records = self
            .records
            .checked_add(1)
            .filter(|count| *count <= self.max_records)
            .ok_or_else(|| bad("mzIdentML output element limit exceeded"))?;
        Ok(())
    }
    fn indent(&mut self, depth: usize) -> Result<()> {
        for _ in 0..depth {
            self.raw("\t")?;
        }
        Ok(())
    }
    fn attr(&mut self, name: &str, value: &str) -> Result<()> {
        super::identification_xml::xml_text(value)?;
        self.raw(" ")?;
        self.raw(name)?;
        self.raw("=\"")?;
        self.raw(&quick_xml::escape::escape(value))?;
        self.raw("\"")
    }
    fn text_content(&mut self, value: &str) -> Result<()> {
        super::identification_xml::xml_text(value)?;
        self.raw(&quick_xml::escape::escape(value))
    }
}

/// Shortest decimal that reads back exactly.
///
/// The source renders every number through `String(double)`, whose six
/// significant digits silently truncate an m/z or a score; this keeps the
/// value.
fn number_text(value: f64) -> Result<String> {
    if !value.is_finite() {
        return Err(invalid("mzIdentML numbers must be finite"));
    }
    Ok(value.to_string())
}

/// One `PeptideEvidence` element to emit, deduplicated by its own content.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct EvidenceOut {
    peptide: usize,
    protein: usize,
    start: Option<usize>,
    end: Option<usize>,
    pre: Option<char>,
    post: Option<char>,
    decoy: Option<bool>,
}

/// The cross-linker `Modification` an XL `Peptide` carries on top of its
/// sequence's own modifications.
#[derive(Clone, Debug)]
enum XlLinker {
    /// The donor half: the cross-linker itself, with the mass the match
    /// reports for it.
    Donor {
        position: i64,
        specificity: String,
        name: String,
        mass: f64,
    },
    /// The acceptor half on the beta chain, written with a zero mass delta
    /// because the linker's mass is already accounted for on the donor.
    Acceptor { position: i64, specificity: String },
    /// A mono-link whose linker is not a sequence modification: the source's
    /// failsafe placeholder, written as `XLMOD:XXXXX`.
    Placeholder { position: i64, name: String },
}

/// One `Peptide` element of a cross-linking document.
///
/// A cross-linking `Peptide` is not shared between matches: the source appends
/// a fresh link id to its identity key, so every chain of every hit gets its
/// own element, and its cross-linker `Modification` carries the donor or
/// acceptor cvParam that ties the two chains together.
#[derive(Clone, Debug)]
struct XlPeptide {
    sequence: AASequence,
    /// The donor (alpha) chain, as opposed to the acceptor (beta) one.
    alpha: bool,
    /// Positional link id shared by the two chains of one match.
    link: usize,
    /// The cross-linker modification, absent for a mono-link whose linker is
    /// already a modification of the sequence.
    linker: Option<XlLinker>,
    /// The acceptor position of a loop-link, which sits on the same chain.
    loop_position: Option<i64>,
}

/// One `SpectrumIdentificationItem` of a cross-linking document.
#[derive(Clone, Copy, Debug)]
struct XlItem {
    /// Index into [`Plan::xl_peptides`].
    peptide: usize,
    /// Index of the identification in the document.
    identification: usize,
    /// Index of the hit within that identification.
    hit: usize,
    /// Neutral mass of the whole match: both chains plus the linker, resolved
    /// against the caller's registry while planning.
    mass: f64,
    /// The heavy half of a labelled pair, which repeats the light item at the
    /// heavy precursor's m/z and retention time.
    heavy: bool,
    /// The donor chain, which carries the fragment annotations.
    alpha: bool,
}

/// One `SpectrumIdentificationResult` of a cross-linking document.
#[derive(Clone, Debug)]
struct XlResult {
    run: usize,
    spectrum_id: String,
    items: Vec<XlItem>,
}

/// Everything the emitter needs, collected before a single byte is written.
#[derive(Default)]
struct Plan {
    software: Vec<(String, String)>,
    software_index: BTreeMap<String, usize>,
    databases: Vec<(String, String)>,
    database_index: BTreeMap<String, usize>,
    spectra: Vec<String>,
    spectra_index: BTreeMap<String, usize>,
    proteins: Vec<(String, String, String, usize)>,
    protein_index: BTreeMap<String, usize>,
    peptides: Vec<AASequence>,
    peptide_index: BTreeMap<String, usize>,
    evidences: Vec<EvidenceOut>,
    evidence_index: BTreeMap<EvidenceOut, usize>,
    /// Run index of every identification run, by its identifier.
    runs: BTreeMap<String, usize>,
    /// True when a run declares `MS:1002494` or `is_cross_linking_experiment`,
    /// which selects the cross-linking write path.
    crosslinking: bool,
    /// The `Peptide` elements of a cross-linking document.
    xl_peptides: Vec<XlPeptide>,
    /// `PeptideEvidence` elements of a cross-linking document, and the peptide
    /// each belongs to.
    xl_evidences: Vec<(usize, EvidenceOut)>,
    /// `PeptideEvidence` indices per `xl_peptides` entry.
    xl_evidence_refs: BTreeMap<usize, Vec<usize>>,
    /// The results, one per spectrum reference and run.
    xl_results: Vec<XlResult>,
    /// `cross_link:mass` and `cross_link:mass_isoshift` of the first run.
    xl_masses: (f64, f64),
}

fn insert_index(index: &mut BTreeMap<String, usize>, key: &str, next: usize) -> usize {
    *index.entry(key.to_owned()).or_insert(next)
}

/// Remove the `[` and `]` OpenMS wraps a file URI in, as `trimOpenMSfileURI`,
/// and normalise backslashes to forward slashes.
fn trim_file_uri(value: &str) -> String {
    let value = value.strip_prefix('[').unwrap_or(value);
    let value = value.strip_suffix(']').unwrap_or(value);
    value.replace('\\', "/")
}

/// Preflight the ceilings and every reference before anything is emitted.
fn preflight(document: &MzIdentMLDocument, options: &WriteOptions) -> Result<()> {
    if options.max_output_bytes == 0 || options.max_records == 0 {
        return Err(invalid("invalid mzIdentML write limits"));
    }
    if document.protein_identifications.len() > MAX_ITEMS
        || document.peptide_identifications.len() > MAX_ITEMS
    {
        return Err(bad("mzIdentML identification count exceeds MAX_ITEMS"));
    }
    if document.protein_identifications.is_empty() {
        // AnalysisCollection requires at least one SpectrumIdentification.
        return Err(missing(
            "mzIdentML needs at least one protein identification run",
        ));
    }
    let mut identifiers = BTreeSet::new();
    for run in &document.protein_identifications {
        run.validate()?;
        if !identifiers.insert(run.identifier.as_str()) {
            return Err(invalid(format!(
                "protein identification run identifier {:?} is not unique",
                run.identifier
            )));
        }
    }
    for identification in &document.peptide_identifications {
        identification.validate()?;
        if !identifiers.contains(identification.identifier.as_str()) {
            // Source: logs "encountered a PeptideIdentification which is not
            // linked to any ProteinIdentification" and drops the spectrum.
            return Err(invalid(format!(
                "peptide identification references unknown run identifier {:?}",
                identification.identifier
            )));
        }
        if identification.hits.is_empty() {
            // SpectrumIdentificationResult requires at least one item.
            return Err(missing(
                "a peptide identification with no hits cannot be written as mzIdentML",
            ));
        }
    }
    Ok(())
}

/// Refuse a plan that has already outgrown the write ceiling.
///
/// `max_records` bounds the emitted elements, and [`Out::element`] enforces it
/// while they are written; every element the plan holds becomes at least one of
/// them, so the same ceiling is applied here, before the plan has grown past
/// it, rather than only once output is produced.
fn planned(count: usize, options: &WriteOptions) -> Result<()> {
    if count >= options.max_records {
        return Err(bad("mzIdentML record limit exceeded"));
    }
    Ok(())
}

fn plan(
    document: &MzIdentMLDocument,
    options: &WriteOptions,
    registry: &ModificationsDB,
) -> Result<Plan> {
    let mut plan = Plan::default();
    for (index, run) in document.protein_identifications.iter().enumerate() {
        plan.runs.insert(run.identifier.clone(), index);
        let software = insert_index(
            &mut plan.software_index,
            &run.search_engine,
            plan.software.len(),
        );
        if software == plan.software.len() {
            plan.software
                .push((run.search_engine.clone(), run.search_engine_version.clone()));
        }
        let database = insert_index(
            &mut plan.database_index,
            &run.search_parameters.database,
            plan.databases.len(),
        );
        if database == plan.databases.len() {
            plan.databases.push((
                run.search_parameters.database.clone(),
                run.search_parameters.database_version.clone(),
            ));
        }
        let location = match run.metadata.get("spectra_data") {
            Some(value) => value
                .as_string_list()
                .ok()
                .and_then(|list| list.first().cloned())
                .unwrap_or_else(|| text_of(value)),
            None => String::new(),
        };
        let location = if location.is_empty() {
            "UNKNOWN".to_owned()
        } else {
            trim_file_uri(&location)
        };
        let spectra = insert_index(&mut plan.spectra_index, &location, plan.spectra.len());
        if spectra == plan.spectra.len() {
            plan.spectra.push(location);
        }
        for hit in &run.hits {
            let next = plan.proteins.len();
            if insert_index(&mut plan.protein_index, &hit.accession, next) == next {
                let description = hit.description();
                let description = if description.is_empty() {
                    hit.accession.clone()
                } else {
                    description
                };
                plan.proteins.push((
                    hit.accession.clone(),
                    hit.sequence.clone(),
                    description,
                    database,
                ));
            }
        }
    }
    plan.crosslinking = document
        .protein_identifications
        .iter()
        .any(is_crosslinking_run);
    if plan.crosslinking {
        plan_crosslinks(document, &mut plan, options, registry)?;
        return Ok(plan);
    }
    for identification in &document.peptide_identifications {
        for hit in &identification.hits {
            planned(
                plan.peptides.len().saturating_add(plan.evidences.len()),
                options,
            )?;
            let key = hit.sequence.to_string();
            let peptide = insert_index(&mut plan.peptide_index, &key, plan.peptides.len());
            if peptide == plan.peptides.len() {
                // Reject a sequence whose modifications cannot be named in the
                // caller's registry before any output is produced.
                let _ = AASequence::parse_with_registry(&key, registry)?;
                plan.peptides.push(hit.sequence.clone());
            }
            let decoy = hit
                .metadata
                .get(user_param::TARGET_DECOY)
                .map(|value| value.as_str().unwrap_or_default().contains("decoy"));
            for evidence in &hit.evidences {
                let protein =
                    *plan.protein_index.get(&evidence.protein_accession).ok_or_else(|| {
                        // Source: logs "Missing or invalid protein reference for
                        // peptide" and silently drops the evidence.
                        invalid(format!(
                            "peptide evidence references accession {:?}, which no protein hit declares",
                            evidence.protein_accession
                        ))
                    })?;
                let start = match (evidence.start, hit.metadata.get("start")) {
                    (Some(start), _) => Some(start),
                    (None, Some(value)) => usize::try_from(value.as_i64()?).ok(),
                    (None, None) => None,
                };
                let end = match (evidence.end, hit.metadata.get("end")) {
                    (Some(end), _) => Some(end),
                    (None, Some(value)) => usize::try_from(value.as_i64()?).ok(),
                    (None, None) => None,
                };
                let item = EvidenceOut {
                    peptide,
                    protein,
                    start,
                    end,
                    pre: match evidence.aa_before {
                        FlankingResidue::Unknown => None,
                        other => Some(other.code()),
                    },
                    post: match evidence.aa_after {
                        FlankingResidue::Unknown => None,
                        other => Some(other.code()),
                    },
                    decoy,
                };
                let next = plan.evidences.len();
                if *plan.evidence_index.entry(item.clone()).or_insert(next) == next {
                    plan.evidences.push(item);
                }
            }
        }
    }
    Ok(plan)
}

/// The source's `is_ppxl` test: a run that declares the crosslinking search
/// term, or the `is_cross_linking_experiment` flag an OpenPepXL tool sets.
fn is_crosslinking_run(run: &ProteinIdentification) -> bool {
    run.metadata.contains_key("is_cross_linking_experiment")
        || run
            .metadata
            .get("SpectrumIdentificationProtocol")
            .map(text_of)
            .as_deref()
            == Some(CROSSLINKING_SEARCH)
}

/// The `spectrumID` attribute of one identification's result element.
///
/// The source falls back to `MZ:<mz>@RT:<rt>` with `nan` for a missing
/// coordinate; `nan` is not a valid `xs:double` lexical form there either, so
/// this writes the schema's `NaN` spelling.
fn spectrum_id_of(identification: &PeptideIdentification) -> Result<String> {
    let reference = identification.spectrum_reference();
    if !reference.is_empty() {
        return Ok(reference);
    }
    let mz = identification
        .mz
        .map(number_text)
        .transpose()?
        .unwrap_or_else(|| "NaN".to_owned());
    let rt = identification
        .rt
        .map(number_text)
        .transpose()?
        .unwrap_or_else(|| "NaN".to_owned());
    Ok(format!("MZ:{mz}@RT:{rt}"))
}

/// One flanking residue of a `BetaPepEv:` list, or `None` for the unknown
/// marker the source writes as `X`.
fn beta_flanking(text: &str) -> Option<char> {
    let code = text.chars().next()?;
    if code == 'X' { None } else { Some(code) }
}

/// One position of a `BetaPepEv:` list, or `None` for the unknown marker.
fn beta_position(text: &str) -> Option<usize> {
    match text.trim().parse::<i64>() {
        Ok(value) if value >= 0 => usize::try_from(value).ok(),
        _ => None,
    }
}

/// Plan the `Peptide`, `PeptideEvidence` and `SpectrumIdentificationResult`
/// elements of a cross-linking document, as the `is_ppxl` branches of
/// `writeTo` and `writeXLMSPeptideHit`.
///
/// Every hit becomes one `Peptide` per chain, so nothing is shared between
/// matches: the source makes its peptide identity key unique per hit by
/// appending the link id, which leaves its own de-duplication branch dead.
///
/// # Errors
///
/// [`Error::InvalidValue`] when a hit is missing the metadata its own
/// `xl_type` requires - the cross-linker name and mass of a cross-link, the
/// beta sequence, or a link position - all of which the source reads with an
/// unchecked `getMetaValue`, and when a beta accession names no protein hit.
fn plan_crosslinks(
    document: &MzIdentMLDocument,
    plan: &mut Plan,
    options: &WriteOptions,
    registry: &ModificationsDB,
) -> Result<()> {
    let first = document
        .protein_identifications
        .first()
        .ok_or_else(|| missing("mzIdentML needs at least one protein identification run"))?;
    let number = |key: &str| -> Result<f64> {
        match first.search_parameters.metadata.get(key) {
            None => Ok(0.0),
            Some(value) => match value.data() {
                MetaValueData::Empty => Ok(0.0),
                MetaValueData::String(text) if text.is_empty() => Ok(0.0),
                MetaValueData::String(text) => finite(text),
                _ => value.as_f64(),
            },
        }
    };
    // The source defaults both to 0 so a document stored after a load that did
    // not preserve them still writes.
    plan.xl_masses = (
        number("cross_link:mass")?,
        number("cross_link:mass_isoshift")?,
    );
    let mut grouped: BTreeMap<(usize, String), Vec<XlItem>> = BTreeMap::new();
    for (index, identification) in document.peptide_identifications.iter().enumerate() {
        let run = *plan
            .runs
            .get(&identification.identifier)
            .ok_or_else(|| invalid("internal run index is incomplete"))?;
        let mut spectrum_id = spectrum_id_of(identification)?;
        if let Some(heavy) = identification
            .metadata
            .get(user_param::OPENPEPXL_HEAVY_SPEC_REF)
        {
            spectrum_id.push(',');
            spectrum_id.push_str(&text_of(heavy));
        }
        let items = grouped.entry((run, spectrum_id)).or_default();
        for (hit_index, hit) in identification.hits.iter().enumerate() {
            planned(
                plan.xl_peptides
                    .len()
                    .saturating_add(plan.xl_evidences.len()),
                options,
            )?;
            let link = plan.xl_peptides.len();
            let kind = metadata_text(&hit.metadata, user_param::OPENPEPXL_XL_TYPE);
            // A labelled match reports the heavy precursor as a second item.
            let labelled = hit
                .metadata
                .contains_key(user_param::OPENPEPXL_HEAVY_SPEC_RT)
                && hit
                    .metadata
                    .contains_key(user_param::OPENPEPXL_HEAVY_SPEC_MZ);
            let alpha_position = crosslink_position(&hit.metadata, user_param::OPENPEPXL_XL_POS1)?;
            let beta_position = crosslink_position(&hit.metadata, user_param::OPENPEPXL_XL_POS2)?;
            let linker_name = hit
                .metadata
                .get(user_param::OPENPEPXL_XL_MOD)
                .map(text_of)
                .filter(|name| !name.is_empty());
            let linker = match (kind.as_str(), alpha_position) {
                ("cross-link" | "loop-link", Some(position)) => {
                    let name = linker_name.clone().ok_or_else(|| {
                        invalid("a cross-linked hit needs an xl_mod cross-linker name")
                    })?;
                    let mass = hit
                        .metadata
                        .get(user_param::OPENPEPXL_XL_MASS)
                        .ok_or_else(|| {
                            invalid("a cross-linked hit needs an xl_mass cross-linker mass")
                        })?
                        .as_f64()?;
                    Some(XlLinker::Donor {
                        position,
                        specificity: metadata_text(
                            &hit.metadata,
                            user_param::OPENPEPXL_XL_TERM_SPEC_ALPHA,
                        ),
                        name,
                        mass,
                    })
                }
                ("cross-link" | "loop-link", None) => {
                    return Err(invalid("a cross-linked hit needs an xl_pos1 link position"));
                }
                (_, Some(position)) => match &linker_name {
                    // The source's failsafe: a mono-linker that is not a
                    // sequence modification is written as a placeholder.
                    Some(name)
                        if usize::try_from(position).is_ok_and(|index| {
                            hit.sequence
                                .residue_modification(index)
                                .ok()
                                .flatten()
                                .is_none()
                        }) =>
                    {
                        Some(XlLinker::Placeholder {
                            position,
                            name: name.clone(),
                        })
                    }
                    _ => None,
                },
                (_, None) => None,
            };
            let loop_position = if kind == "loop-link" {
                beta_position
            } else {
                None
            };
            let beta = match kind.as_str() {
                "cross-link" => {
                    let text = hit
                        .metadata
                        .get(user_param::OPENPEPXL_BETA_SEQUENCE)
                        .map(text_of)
                        .filter(|text| !text.is_empty() && text != "-")
                        .ok_or_else(|| invalid("a cross-link needs a sequence_beta beta chain"))?;
                    // The beta chain is resolved against the caller's registry,
                    // like every other sequence this writer touches.
                    Some(AASequence::parse_with_registry(&text, registry)?)
                }
                _ => None,
            };
            // The calculated mass of a cross-link is both chains plus the
            // linker; a mono-linker that is not a sequence modification is not
            // part of the sequence's own monoisotopic mass either.
            let mut mass = hit.sequence.mono_mass()?;
            match &beta {
                Some(beta) => {
                    mass += plan.xl_masses.0;
                    mass += beta.mono_mass()?;
                }
                None => {
                    if linker_name.is_some() {
                        if let Some(value) = hit.metadata.get(user_param::OPENPEPXL_XL_MASS) {
                            mass += value.as_f64()?;
                        }
                    }
                }
            }
            plan.xl_peptides.push(XlPeptide {
                sequence: hit.sequence.clone(),
                alpha: true,
                link,
                linker,
                loop_position,
            });
            let alpha_index = plan.xl_peptides.len().saturating_sub(1);
            plan_xl_alpha_evidences(plan, hit, alpha_index)?;
            for heavy in [false, true] {
                items.push(XlItem {
                    peptide: alpha_index,
                    identification: index,
                    hit: hit_index,
                    mass,
                    heavy,
                    alpha: true,
                });
                if !labelled {
                    break;
                }
            }
            let Some(sequence) = beta else {
                continue;
            };
            let acceptor = beta_position
                .ok_or_else(|| invalid("a cross-link needs an xl_pos2 link position"))?;
            plan.xl_peptides.push(XlPeptide {
                sequence,
                alpha: false,
                link,
                linker: Some(XlLinker::Acceptor {
                    position: acceptor,
                    specificity: metadata_text(
                        &hit.metadata,
                        user_param::OPENPEPXL_XL_TERM_SPEC_BETA,
                    ),
                }),
                loop_position: None,
            });
            let beta_index = plan.xl_peptides.len().saturating_sub(1);
            plan_xl_beta_evidences(plan, hit, beta_index)?;
            for heavy in [false, true] {
                items.push(XlItem {
                    peptide: beta_index,
                    identification: index,
                    hit: hit_index,
                    mass,
                    heavy,
                    alpha: false,
                });
                if !labelled {
                    break;
                }
            }
        }
    }
    plan.xl_results = grouped
        .into_iter()
        .map(|((run, spectrum_id), items)| XlResult {
            run,
            spectrum_id,
            items,
        })
        .collect();
    Ok(())
}

/// The alpha chain's `PeptideEvidence` elements, from the hit's own evidences.
fn plan_xl_alpha_evidences(plan: &mut Plan, hit: &PeptideHit, peptide: usize) -> Result<()> {
    let decoy = hit
        .metadata
        .get(user_param::OPENPEPXL_TARGET_DECOY_ALPHA)
        .map(|value| text_of(value).contains("decoy"));
    for evidence in &hit.evidences {
        let protein = *plan
            .protein_index
            .get(&evidence.protein_accession)
            .ok_or_else(|| {
                // Source: logs "Missing or invalid protein reference for
                // peptide" and silently drops the evidence.
                invalid(format!(
                    "peptide evidence references accession {:?}, which no protein hit declares",
                    evidence.protein_accession
                ))
            })?;
        let start = match (evidence.start, hit.metadata.get("start")) {
            (Some(start), _) => Some(start),
            (None, Some(value)) => usize::try_from(value.as_i64()?).ok(),
            (None, None) => None,
        };
        let end = match (evidence.end, hit.metadata.get("end")) {
            (Some(end), _) => Some(end),
            (None, Some(value)) => usize::try_from(value.as_i64()?).ok(),
            (None, None) => None,
        };
        let item = EvidenceOut {
            peptide,
            protein,
            start,
            end,
            pre: match evidence.aa_before {
                FlankingResidue::Unknown => None,
                other => Some(other.code()),
            },
            post: match evidence.aa_after {
                FlankingResidue::Unknown => None,
                other => Some(other.code()),
            },
            decoy,
        };
        let index = plan.xl_evidences.len();
        plan.xl_evidences.push((peptide, item));
        plan.xl_evidence_refs
            .entry(peptide)
            .or_default()
            .push(index);
    }
    Ok(())
}

/// The beta chain's `PeptideEvidence` elements, from the `BetaPepEv:` metadata
/// the merge step left on the alpha hit.
///
/// The source indexes its five parallel lists by the same position without a
/// length check; this stops at the shortest.
fn plan_xl_beta_evidences(plan: &mut Plan, hit: &PeptideHit, peptide: usize) -> Result<()> {
    let list = |key: &str| -> Vec<String> {
        match hit.metadata.get(key) {
            None => Vec::new(),
            Some(value) => text_of(value)
                .split(',')
                .map(str::trim)
                .map(str::to_owned)
                .collect(),
        }
    };
    let accessions = list(user_param::OPENPEPXL_BETA_ACCESSIONS);
    let pre = list(user_param::OPENPEPXL_BETA_PEPEV_PRE);
    let post = list(user_param::OPENPEPXL_BETA_PEPEV_POST);
    let start = list(user_param::OPENPEPXL_BETA_PEPEV_START);
    let end = list(user_param::OPENPEPXL_BETA_PEPEV_END);
    let decoy = hit
        .metadata
        .get(user_param::OPENPEPXL_TARGET_DECOY_BETA)
        .map(|value| text_of(value).contains("decoy"));
    let count = accessions
        .len()
        .min(pre.len())
        .min(post.len())
        .min(start.len())
        .min(end.len());
    for index in 0..count {
        let accession = &accessions[index];
        if accession.is_empty() || accession == "-" {
            continue;
        }
        let protein = *plan.protein_index.get(accession).ok_or_else(|| {
            invalid(format!(
                "beta chain references accession {accession:?}, which no protein hit declares"
            ))
        })?;
        let item = EvidenceOut {
            peptide,
            protein,
            start: beta_position(&start[index]),
            end: beta_position(&end[index]),
            pre: beta_flanking(&pre[index]),
            post: beta_flanking(&post[index]),
            decoy,
        };
        let position = plan.xl_evidences.len();
        plan.xl_evidences.push((peptide, item));
        plan.xl_evidence_refs
            .entry(peptide)
            .or_default()
            .push(position);
    }
    Ok(())
}

/// Emit one `cvParam` for a CV term, with an optional value.
fn cv_param(out: &mut Out, term: &CVTermDefinition, value: Option<&MetaValue>) -> Result<()> {
    out.element()?;
    let text = match value {
        Some(value) => term.to_xml_value("PSI-MS", value)?,
        None => term.to_xml("PSI-MS", "")?,
    };
    out.raw(&text)
}

fn cv_param_named(
    out: &mut Out,
    cv: &ControlledVocabulary,
    name: &str,
    value: Option<&MetaValue>,
) -> Result<()> {
    let term = cv
        .find_term_by_name(name)
        .ok_or_else(|| invalid(format!("PSI-MS has no term named {name:?}")))?;
    cv_param(out, term, value)
}

/// `writeMetaInfos_`: a `cvParam` for a key that is a PSI-MS accession, and a
/// `userParam` otherwise.
///
/// The source declares the XSD type of a `userParam` in `unitName`, which the
/// reader does not read, so every typed value degrades to text on re-read (see
/// `OpenMS_CPP_ISSUES.md`). This writes the schema's `type` attribute.
fn write_meta(
    out: &mut Out,
    cv: &ControlledVocabulary,
    metadata: &MetaInfo,
    skip: &BTreeSet<String>,
    depth: usize,
) -> Result<()> {
    for (key, value) in metadata {
        if skip.contains(key) {
            continue;
        }
        out.indent(depth)?;
        match cv.get_term(key) {
            Ok(term) => cv_param(out, term, Some(value))?,
            Err(_) => {
                out.element()?;
                out.raw("<userParam")?;
                out.attr("name", key)?;
                let kind = match value.data() {
                    MetaValueData::Integer(_) => "xsd:integer",
                    MetaValueData::Float(_) => "xsd:double",
                    _ => "xsd:string",
                };
                out.attr("type", kind)?;
                out.attr("value", &text_of(value))?;
                out.raw("/>")?;
            }
        }
        out.raw("\n")?;
    }
    Ok(())
}

/// `writeModParam_`: one `SearchModification` per registered modification name.
///
/// # Errors
///
/// [`Error::InvalidValue`] when a name is not in `registry`, as the source
/// throws `Exception::ElementNotFound`, and when it is ambiguous. The source
/// resolves through `searchModifications`, which returns a set and then writes
/// *every* match, so an ambiguous name multiplies into several elements.
fn write_mod_params(
    out: &mut Out,
    cv: &ControlledVocabulary,
    names: &[String],
    fixed: bool,
    registry: &ModificationsDB,
    depth: usize,
) -> Result<()> {
    for name in names {
        let record = registry.get_modification(name, None, None)?;
        let origin = record.origin().filter(|c| *c != 'X').unwrap_or('.');
        out.indent(depth)?;
        out.element()?;
        out.raw("<SearchModification")?;
        out.attr("fixedMod", if fixed { "true" } else { "false" })?;
        out.attr("massDelta", &number_text(record.diff_mono_mass())?)?;
        out.attr("residues", &origin.to_string())?;
        out.raw(">\n")?;
        // The source writes specificity rules for the peptide termini only and
        // leaves a "@TODO: handle protein C-term/N-term"; without them a
        // protein-terminal modification cannot be resolved unambiguously when
        // the file is read back, so all four are written here.
        let rule = match record.term_specificity() {
            TermSpecificity::NTerm => Some("modification specificity peptide N-term"),
            TermSpecificity::CTerm => Some("modification specificity peptide C-term"),
            TermSpecificity::ProteinNTerm => Some("modification specificity protein N-term"),
            TermSpecificity::ProteinCTerm => Some("modification specificity protein C-term"),
            TermSpecificity::Anywhere => None,
        };
        if let Some(rule) = rule {
            out.indent(depth + 1)?;
            out.raw("<SpecificityRules>\n")?;
            out.indent(depth + 2)?;
            cv_param_named(out, cv, rule, None)?;
            out.raw("\n")?;
            out.indent(depth + 1)?;
            out.raw("</SpecificityRules>\n")?;
        }
        out.indent(depth + 1)?;
        out.element()?;
        match record.unimod_accession() {
            Some(accession) => {
                // The source emits this cvParam from the bundled unimod.obo;
                // the crate has no UniMod OBO, so the accession and name come
                // from the ModificationsDB record itself.
                let accession = accession.replace("UniMod:", "UNIMOD:");
                out.raw("<cvParam")?;
                out.attr("accession", &accession)?;
                out.attr("cvRef", "UNIMOD")?;
                out.attr("name", record.name())?;
                out.raw("/>\n")?;
            }
            None => {
                out.raw("<cvParam cvRef=\"MS\" accession=\"MS:1001460\" name=\"unknown modification\"/>\n")?;
            }
        }
        out.indent(depth)?;
        out.raw("</SearchModification>\n")?;
    }
    Ok(())
}

/// `writeEnzyme_`: the `Enzymes` block of one protocol.
fn write_enzyme(
    out: &mut Out,
    cv: &ControlledVocabulary,
    name: &str,
    missed_cleavages: u32,
    depth: usize,
) -> Result<()> {
    out.indent(depth)?;
    out.element()?;
    out.raw("<Enzymes independent=\"false\">\n")?;
    out.indent(depth + 1)?;
    out.element()?;
    out.raw("<Enzyme")?;
    out.attr("missedCleavages", &missed_cleavages.to_string())?;
    out.attr("id", "ENZ_0")?;
    out.raw(">\n")?;
    out.indent(depth + 2)?;
    out.raw("<EnzymeName>\n")?;
    out.indent(depth + 3)?;
    let term = if cv.has_term_with_name(name) {
        name
    } else if name == "no cleavage" {
        "NoEnzyme"
    } else {
        "cleavage agent details"
    };
    cv_param_named(out, cv, term, None)?;
    out.raw("\n")?;
    out.indent(depth + 2)?;
    out.raw("</EnzymeName>\n")?;
    out.indent(depth + 1)?;
    out.raw("</Enzyme>\n")?;
    out.indent(depth)?;
    out.raw("</Enzymes>\n")
}

/// Split one OpenMS peak annotation into its mzIdentML ion type, series index
/// and neutral loss, as the source's `frag_regex_tweak` does.
///
/// Accepts `[chain|category$<abcxyz><index>(-H2O|-NH3)*]` with an optional
/// trailing charge suffix, and nothing else; annotations that do not fit the
/// limited mzIdentML fragment structure are skipped, as in the source.
fn split_annotation(annotation: &str) -> Option<(String, char, String, Option<String>)> {
    let body = annotation.strip_prefix('[')?;
    let (body, tail) = body.split_once(']')?;
    if !tail
        .chars()
        .all(|c| c.is_ascii_digit() || matches!(c, '+' | '(' | ')'))
    {
        return None;
    }
    let (extra, core) = match body.rsplit_once('$') {
        Some((extra, core)) => (extra.to_owned(), core),
        None => (String::new(), body),
    };
    let mut loss = None;
    let mut core = core;
    for candidate in ["H2O", "NH3"] {
        if let Some(head) = core.strip_suffix(candidate) {
            if let Some(head) = head.strip_suffix('-') {
                loss = Some(candidate.to_owned());
                core = head;
                break;
            }
        }
    }
    let mut characters = core.chars();
    let kind = characters.next()?;
    if !matches!(kind, 'a' | 'b' | 'c' | 'x' | 'y' | 'z') {
        return None;
    }
    let index: String = characters.collect();
    if index.is_empty() || !index.chars().all(|c| c.is_ascii_digit()) {
        return None;
    }
    Some((extra, kind, index, loss))
}

/// `writeFragmentAnnotations_`: the `Fragmentation` block of one PSM.
///
/// Annotations are grouped by charge and then by ion type, as in the source,
/// and the product ion m/z error array the source's own example comments show
/// is not written, because [`PeakAnnotation`] carries no error term.
///
/// `crosslinking` adds the `cross-link_chain` and `cross-link_ioncategory`
/// `userParam`s the source writes for XL-MS results, one entry per index: the
/// chain is `alpha` unless the annotation names `beta`, and the category `xi`
/// unless it names `ci`, which is how `writeFragmentAnnotations_` reduces the
/// annotation's prefix.
fn write_fragmentation(
    out: &mut Out,
    cv: &ControlledVocabulary,
    annotations: &[PeakAnnotation],
    crosslinking: bool,
    depth: usize,
) -> Result<()> {
    type Series = (
        Vec<String>,
        Vec<String>,
        Vec<String>,
        Vec<String>,
        Vec<String>,
    );
    let mut grouped: BTreeMap<i32, BTreeMap<String, Series>> = BTreeMap::new();
    for annotation in annotations {
        let Some((extra, kind, index, loss)) = split_annotation(&annotation.annotation) else {
            continue;
        };
        let mut name = format!("frag: {kind} ion");
        if let Some(loss) = loss {
            name.push_str(" - ");
            name.push_str(&loss);
        }
        let series = grouped
            .entry(annotation.charge)
            .or_default()
            .entry(name)
            .or_default();
        series.0.push(index);
        series.1.push(number_text(annotation.mz)?);
        series.2.push(number_text(annotation.intensity)?);
        if crosslinking {
            let parts: Vec<&str> = extra.split('|').collect();
            series.3.push(
                if parts.contains(&"alpha") {
                    "alpha"
                } else {
                    "beta"
                }
                .to_owned(),
            );
            series
                .4
                .push(if parts.contains(&"ci") { "ci" } else { "xi" }.to_owned());
        }
    }
    if grouped.is_empty() {
        return Ok(());
    }
    out.indent(depth)?;
    out.element()?;
    out.raw("<Fragmentation>\n")?;
    for (charge, series) in &grouped {
        for (name, (index, mz, intensity, chain, category)) in series {
            out.indent(depth + 1)?;
            out.element()?;
            out.raw("<IonType")?;
            out.attr("charge", &charge.to_string())?;
            out.attr("index", &index.join(" "))?;
            out.raw(">\n")?;
            for (measure, values) in [("Measure_mz", mz), ("Measure_int", intensity)] {
                out.indent(depth + 2)?;
                out.element()?;
                out.raw("<FragmentArray")?;
                out.attr("measure_ref", measure)?;
                out.attr("values", &values.join(" "))?;
                out.raw("/>\n")?;
            }
            if crosslinking {
                for (name, values) in [
                    ("cross-link_chain", chain),
                    ("cross-link_ioncategory", category),
                ] {
                    out.indent(depth + 2)?;
                    out.element()?;
                    out.raw("<userParam")?;
                    out.attr("name", name)?;
                    out.attr("type", "xsd:string")?;
                    out.attr("value", &values.join(" "))?;
                    out.raw("/>\n")?;
                }
            }
            out.indent(depth + 2)?;
            cv_param_named(out, cv, name, None)?;
            out.raw("\n")?;
            out.indent(depth + 1)?;
            out.raw("</IonType>\n")?;
        }
    }
    out.indent(depth)?;
    out.raw("</Fragmentation>\n")
}
/// The PSI-MS term name the source maps a search engine name to.
fn software_term(cv: &ControlledVocabulary, engine: &str) -> &'static str {
    match engine {
        "OMSSA" => "OMSSA",
        "Mascot" => "Mascot",
        "XTandem" => "X\\!Tandem",
        "SEQUEST" => "Sequest",
        "MS-GF+" => "MS-GF+",
        "Percolator" => "Percolator",
        "OpenPepXL" => "OpenPepXL",
        _ if cv.has_term_with_name(engine) => "",
        _ => "analysis software",
    }
}

/// Emit the score `cvParam` of one PSM and return the accession it consumed.
///
/// The chain is the source's: a score type that names a PSM-level statistic
/// (the `MS:1001143` subtree) is written as that term, then the
/// engine-specific aliases, and anything left becomes the bare parent term
/// plus a `userParam` carrying the number. A score type that falls through
/// cannot be recovered by a reader, which is a source gap this port keeps
/// rather than invent a term for; `docs/MZIDENTML_SUPPORT.md` records it.
fn write_score(
    out: &mut Out,
    cv: &ControlledVocabulary,
    details: &BTreeSet<String>,
    score_type: &str,
    score: f64,
    depth: usize,
) -> Result<Option<String>> {
    let value = MetaValue::new(MetaValueData::Float(score))?;
    let named = cv
        .find_term_by_name(score_type)
        .filter(|term| details.contains(&term.id));
    let direct = cv
        .get_term(score_type)
        .ok()
        .filter(|term| details.contains(&term.id));
    let alias = match score_type {
        "q-value" | "FDR" => Some("PSM-level q-value"),
        "Posterior Error Probability" => Some("percolator:PEP"),
        "OMSSA" => Some("OMSSA:evalue"),
        "Mascot" => Some("Mascot:score"),
        "XTandem" => Some("X\\!Tandem:hyperscore"),
        "SEQUEST" => Some("Sequest:xcorr"),
        "MS-GF+" => Some("MS-GF:RawScore"),
        other if other == user_param::OPENPEPXL_SCORE => Some(user_param::OPENPEPXL_SCORE),
        _ => None,
    };
    out.indent(depth)?;
    if let Some(term) = named.or(direct) {
        let id = term.id.clone();
        cv_param(out, term, Some(&value))?;
        out.raw("\n")?;
        return Ok(Some(id));
    }
    if let Some(name) = alias {
        let term = cv
            .find_term_by_name(name)
            .ok_or_else(|| invalid(format!("PSI-MS has no term named {name:?}")))?;
        let id = term.id.clone();
        cv_param(out, term, Some(&value))?;
        out.raw("\n")?;
        return Ok(Some(id));
    }
    let placeholder = if score_type.is_empty() {
        "PSM-level search engine specific statistic"
    } else {
        score_type
    };
    cv_param_named(out, cv, "PSM-level search engine specific statistic", None)?;
    out.raw("\n")?;
    out.indent(depth)?;
    out.element()?;
    out.raw("<userParam")?;
    out.attr("name", placeholder)?;
    out.attr("type", "xsd:double")?;
    out.attr("value", &number_text(score)?)?;
    out.raw("/>\n")?;
    Ok(None)
}

/// The vocabulary, the score subtree and the id plan every emit step shares.
struct WriteContext<'a> {
    cv: &'a ControlledVocabulary,
    details: BTreeSet<String>,
    plan: Plan,
}

/// The metadata values `writeXLMSPeptideHit` removes from the hit before it
/// writes the rest: everything the cross-linking structure carries itself.
fn crosslink_skip(consumed: Option<String>) -> BTreeSet<String> {
    let mut skip: BTreeSet<String> = [
        "calcMZ",
        "protein_references",
        "xl_chain",
        user_param::TARGET_DECOY,
        user_param::SPECTRUM_REFERENCE,
        user_param::OPENPEPXL_XL_RANK,
        user_param::OPENPEPXL_XL_POS1,
        user_param::OPENPEPXL_XL_POS2,
        user_param::OPENPEPXL_XL_MOD,
        user_param::OPENPEPXL_XL_MASS,
        user_param::OPENPEPXL_HEAVY_SPEC_REF,
        user_param::OPENPEPXL_HEAVY_SPEC_MZ,
        user_param::OPENPEPXL_HEAVY_SPEC_RT,
        user_param::OPENPEPXL_BETA_PEPEV_PRE,
        user_param::OPENPEPXL_BETA_PEPEV_POST,
        user_param::OPENPEPXL_BETA_PEPEV_START,
        user_param::OPENPEPXL_BETA_PEPEV_END,
        user_param::OPENPEPXL_BETA_SEQUENCE,
        user_param::OPENPEPXL_BETA_ACCESSIONS,
    ]
    .iter()
    .map(|key| (*key).to_owned())
    .collect();
    if let Some(consumed) = consumed {
        skip.insert(consumed);
    }
    skip
}

/// The XLMOD record that describes a cross-linker of mass `mass` at `residue`.
///
/// `writeXLMSPeptideHit` searches `CrossLinksDB` by mass difference, prefers
/// the record whose full id names the match's own `xl_mod`, and falls back to
/// the first record of the same mass. When the terminal search comes back empty
/// the caller falls back to the residue-scoped search. When nothing matches at
/// all, it writes the source's `XLMOD:XXXXX` placeholder.
///
/// The records come from the bundled XLMOD vocabulary through
/// [`CrossLinksDB::global`], which is the singleton the source consults too;
/// the caller's [`ModificationsDB`] names sequence modifications, not linkers.
fn crosslinker_term<'a>(
    crosslinks: &'a ModificationsDB,
    mass: f64,
    residue: Option<char>,
    term: Option<TermSpecificity>,
    name: &str,
) -> Result<Option<&'a ResidueModification>> {
    let matches = crosslinks.search_by_mass(mass, 0.0001, residue, term)?;
    let chosen = matches
        .iter()
        .find(|record| record.full_id().contains(name))
        .or_else(|| matches.first());
    Ok(chosen.copied())
}

/// One `Peptide` element of a cross-linking document, with the cross-linker
/// modification that ties its chain to the match.
fn write_crosslink_peptide(
    out: &mut Out,
    cv: &ControlledVocabulary,
    peptide: &XlPeptide,
    index: usize,
    depth: usize,
) -> Result<()> {
    let crosslinks = CrossLinksDB::global().database();
    let sequence = &peptide.sequence;
    let residue_at = |position: i64| -> Option<char> {
        usize::try_from(position)
            .ok()
            .and_then(|index| sequence.as_str().chars().nth(index))
    };
    out.indent(depth)?;
    out.element()?;
    out.raw("<Peptide")?;
    out.attr("id", &format!("PEP_{index}"))?;
    out.attr("name", &sequence.to_string())?;
    out.raw(">\n")?;
    out.indent(depth + 1)?;
    out.raw("<PeptideSequence>")?;
    out.text_content(sequence.as_str())?;
    out.raw("</PeptideSequence>\n")?;
    write_modifications(out, sequence, depth + 1)?;
    let link = peptide.link.to_string();
    if let Some(linker) = &peptide.linker {
        let (location, residue, mass, accession, name) = match linker {
            XlLinker::Donor {
                position,
                specificity,
                name,
                mass,
            } => {
                let residue = residue_at(*position);
                // The location follows the terminal specificity alone, so a
                // terminal link reads back as one: the source writes location 0
                // (or the position plus two) only when the cross-link database
                // also holds a terminal record of that mass, and emits an
                // attribute list with no element name in front of it when it
                // does not. The database search only picks the accession.
                let (location, term) = match specificity.as_str() {
                    "N_TERM" => (0, Some(TermSpecificity::NTerm)),
                    "C_TERM" => (position.saturating_add(2), Some(TermSpecificity::CTerm)),
                    _ => (position.saturating_add(1), None),
                };
                let record = match term {
                    Some(term) => crosslinker_term(crosslinks, *mass, None, Some(term), name)?,
                    None => None,
                };
                let record = match record {
                    Some(record) => Some(record),
                    None => crosslinker_term(
                        crosslinks,
                        *mass,
                        residue,
                        Some(TermSpecificity::Anywhere),
                        name,
                    )?,
                };
                let (accession, term_name) = match record {
                    Some(record) => (
                        record.obo_accession().unwrap_or("XLMOD:XXXXX").to_owned(),
                        record.name().to_owned(),
                    ),
                    // Source: "if there is no matching modification in the
                    // database, write out a placeholder".
                    None => ("XLMOD:XXXXX".to_owned(), name.clone()),
                };
                (location, residue, Some(*mass), accession, term_name)
            }
            XlLinker::Acceptor {
                position,
                specificity,
            } => {
                let location = match specificity.as_str() {
                    "N_TERM" => 0,
                    // The source writes the peptide length plus two here,
                    // which its own reader cannot map back to the C-terminus;
                    // the chain's own position plus two is what round-trips.
                    "C_TERM" => position.saturating_add(2),
                    _ => position.saturating_add(1),
                };
                (
                    location,
                    residue_at(*position),
                    Some(0.0),
                    String::new(),
                    String::new(),
                )
            }
            XlLinker::Placeholder { position, name } => (
                position.saturating_add(1),
                residue_at(*position),
                None,
                "XLMOD:XXXXX".to_owned(),
                name.clone(),
            ),
        };
        out.indent(depth + 1)?;
        out.element()?;
        out.raw("<Modification")?;
        out.attr("location", &location.to_string())?;
        if let Some(residue) = residue {
            out.attr("residues", &residue.to_string())?;
        }
        if let Some(mass) = mass {
            out.attr("monoisotopicMassDelta", &number_text(mass)?)?;
        }
        out.raw(">\n")?;
        if !accession.is_empty() {
            out.indent(depth + 2)?;
            out.element()?;
            out.raw("<cvParam")?;
            out.attr("accession", &accession)?;
            out.attr("cvRef", "XLMOD")?;
            out.attr("name", &name)?;
            out.raw("/>\n")?;
        }
        if !matches!(linker, XlLinker::Placeholder { .. }) {
            let term = if peptide.alpha {
                CROSSLINK_DONOR
            } else {
                CROSSLINK_ACCEPTOR
            };
            out.indent(depth + 2)?;
            cv_param(
                out,
                cv.get_term(term)?,
                Some(&MetaValue::from(link.clone())),
            )?;
            out.raw("\n")?;
        }
        out.indent(depth + 1)?;
        out.raw("</Modification>\n")?;
    }
    if let Some(position) = peptide.loop_position {
        // The second half of a loop-link is always the acceptor.
        out.indent(depth + 1)?;
        out.element()?;
        out.raw("<Modification")?;
        out.attr("location", &position.saturating_add(1).to_string())?;
        if let Some(residue) = residue_at(position) {
            out.attr("residues", &residue.to_string())?;
        }
        out.attr("monoisotopicMassDelta", "0")?;
        out.raw(">\n")?;
        out.indent(depth + 2)?;
        cv_param(
            out,
            cv.get_term(CROSSLINK_ACCEPTOR)?,
            Some(&MetaValue::from(link)),
        )?;
        out.raw("\n")?;
        out.indent(depth + 1)?;
        out.raw("</Modification>\n")?;
    }
    out.indent(depth)?;
    out.raw("</Peptide>\n")
}

/// One `SpectrumIdentificationItem` of a cross-linking document, as
/// `writeXLMSPeptideHit`.
///
/// The heavy half of a labelled pair repeats the light item with the heavy
/// precursor's m/z and retention time and the calculated m/z shifted by
/// `cross_link:mass_isoshift`; the source produces it by string substitution
/// on the light item, which this builds directly instead.
fn write_crosslink_item(
    out: &mut Out,
    context: &WriteContext<'_>,
    document: &MzIdentMLDocument,
    item: &XlItem,
    id: usize,
    depth: usize,
) -> Result<()> {
    let (cv, details, plan) = (context.cv, &context.details, &context.plan);
    let identification = document
        .peptide_identifications
        .get(item.identification)
        .ok_or_else(|| invalid("internal identification index out of range"))?;
    let hit = identification
        .hits
        .get(item.hit)
        .ok_or_else(|| invalid("internal hit index out of range"))?;
    let run = document
        .protein_identifications
        .get(
            *plan
                .runs
                .get(&identification.identifier)
                .ok_or_else(|| invalid("internal run index is incomplete"))?,
        )
        .ok_or_else(|| invalid("internal run index out of range"))?;
    let threshold = run.significance_threshold;
    let kind = metadata_text(&hit.metadata, user_param::OPENPEPXL_XL_TYPE);
    let rank = match hit.metadata.get(user_param::OPENPEPXL_XL_RANK) {
        Some(value) => text_of(value),
        None => hit.rank.saturating_add(1).to_string(),
    };
    let pass = match hit.metadata.get("pass_threshold") {
        Some(value) => value.as_str().unwrap_or_default() == "true",
        None if threshold != 0.0 => {
            if identification.higher_score_better {
                hit.score > threshold
            } else {
                hit.score < threshold
            }
        }
        None => true,
    };
    let calculated = if hit.charge == 0 || hit.sequence.is_empty() {
        // The source divides by the charge unconditionally, so a charge of 0
        // writes "inf" into a required xsd:double attribute.
        None
    } else {
        let charge = f64::from(hit.charge);
        let mut value = (item.mass + charge * crate::concept::constants::PROTON_MASS_U) / charge;
        if item.heavy {
            // The heavy precursor is the light one plus the label's shift.
            value += plan.xl_masses.1 / charge;
        }
        Some(value)
    };
    let experimental = if item.heavy {
        hit.metadata
            .get(user_param::OPENPEPXL_HEAVY_SPEC_MZ)
            .map(MetaValue::as_f64)
            .transpose()?
    } else {
        identification.mz
    };
    let rt = if item.heavy {
        hit.metadata
            .get(user_param::OPENPEPXL_HEAVY_SPEC_RT)
            .map(MetaValue::as_f64)
            .transpose()?
    } else {
        identification.rt
    };
    out.indent(depth)?;
    out.element()?;
    out.raw("<SpectrumIdentificationItem")?;
    out.attr("passThreshold", if pass { "true" } else { "false" })?;
    out.attr("rank", &rank)?;
    out.attr("peptide_ref", &format!("PEP_{}", item.peptide))?;
    if let Some(calculated) = calculated {
        out.attr("calculatedMassToCharge", &number_text(calculated)?)?;
    }
    let experimental = match experimental {
        Some(value) => number_text(value)?,
        None => "NaN".to_owned(),
    };
    out.attr("experimentalMassToCharge", &experimental)?;
    out.attr("chargeState", &hit.charge.to_string())?;
    out.attr("id", &format!("SII_{id}"))?;
    out.raw(">\n")?;
    for reference in plan
        .xl_evidence_refs
        .get(&item.peptide)
        .into_iter()
        .flatten()
    {
        out.indent(depth + 1)?;
        out.element()?;
        out.raw("<PeptideEvidenceRef")?;
        out.attr("peptideEvidence_ref", &format!("PEV_{reference}"))?;
        out.raw("/>\n")?;
    }
    if item.alpha && !hit.peak_annotations.is_empty() {
        write_fragmentation(out, cv, &hit.peak_annotations, true, depth + 1)?;
    }
    let consumed = write_score(
        out,
        cv,
        details,
        &identification.score_type,
        hit.score,
        depth + 1,
    )?;
    if item.heavy
        || hit
            .metadata
            .contains_key(user_param::OPENPEPXL_HEAVY_SPEC_MZ)
        || kind == "cross-link"
    {
        out.indent(depth + 1)?;
        cv_param(
            out,
            cv.get_term(CROSSLINK_ITEM)?,
            Some(&MetaValue::from(
                plan.xl_peptides
                    .get(item.peptide)
                    .map(|peptide| peptide.link.to_string())
                    .unwrap_or_default(),
            )),
        )?;
        out.raw("\n")?;
    }
    write_meta(out, cv, &hit.metadata, &crosslink_skip(consumed), depth + 1)?;
    if let Some(rt) = rt {
        out.indent(depth + 1)?;
        let value = MetaValue::new(MetaValueData::Float(rt))?.with_unit(Unit::new(
            "UO:0000010",
            "second",
            "UO",
        )?)?;
        cv_param_named(out, cv, "retention time", Some(&value))?;
        out.raw("\n")?;
    }
    out.indent(depth)?;
    out.raw("</SpectrumIdentificationItem>\n")
}

fn write_item(
    out: &mut Out,
    context: &WriteContext<'_>,
    identification: &PeptideIdentification,
    hit: &PeptideHit,
    index: usize,
    threshold: f64,
    depth: usize,
) -> Result<()> {
    let (cv, details, plan) = (context.cv, &context.details, &context.plan);
    let sequence = hit.sequence.to_string();
    let peptide = *plan
        .peptide_index
        .get(&sequence)
        .ok_or_else(|| invalid("internal peptide index is incomplete"))?;
    let pass = match hit.metadata.get("pass_threshold") {
        Some(value) => value.as_str().unwrap_or_default() == "true",
        // The source evaluates the run's threshold only when one was set.
        None if threshold != 0.0 => {
            if identification.higher_score_better {
                hit.score > threshold
            } else {
                hit.score < threshold
            }
        }
        None => true,
    };
    let charge = hit.charge;
    let calculated = if charge != 0 && !hit.sequence.is_empty() {
        Some(hit.sequence.mz(charge)?)
    } else {
        // The source divides by the charge unconditionally, so a charge of 0
        // writes "inf" into a required xsd:double attribute.
        None
    };
    out.indent(depth)?;
    out.element()?;
    out.raw("<SpectrumIdentificationItem")?;
    out.attr("passThreshold", if pass { "true" } else { "false" })?;
    out.attr("rank", &hit.rank.saturating_add(1).to_string())?;
    out.attr("peptide_ref", &format!("PEP_{peptide}"))?;
    if let Some(calculated) = calculated {
        out.attr("calculatedMassToCharge", &number_text(calculated)?)?;
    }
    // experimentalMassToCharge is required; the source writes "nan" for an
    // absent m/z, which is not a valid xsd:double lexical form.
    let experimental = match identification.mz {
        Some(mz) => number_text(mz)?,
        None => "NaN".to_owned(),
    };
    out.attr("experimentalMassToCharge", &experimental)?;
    out.attr("chargeState", &charge.to_string())?;
    out.attr("id", &format!("SII_{index}"))?;
    out.raw(">\n")?;
    for evidence in &hit.evidences {
        let protein = *plan
            .protein_index
            .get(&evidence.protein_accession)
            .ok_or_else(|| invalid("internal protein index is incomplete"))?;
        let start = match (evidence.start, hit.metadata.get("start")) {
            (Some(start), _) => Some(start),
            (None, Some(value)) => usize::try_from(value.as_i64()?).ok(),
            (None, None) => None,
        };
        let end = match (evidence.end, hit.metadata.get("end")) {
            (Some(end), _) => Some(end),
            (None, Some(value)) => usize::try_from(value.as_i64()?).ok(),
            (None, None) => None,
        };
        let decoy = hit
            .metadata
            .get(user_param::TARGET_DECOY)
            .map(|value| value.as_str().unwrap_or_default().contains("decoy"));
        let item = EvidenceOut {
            peptide,
            protein,
            start,
            end,
            pre: match evidence.aa_before {
                FlankingResidue::Unknown => None,
                other => Some(other.code()),
            },
            post: match evidence.aa_after {
                FlankingResidue::Unknown => None,
                other => Some(other.code()),
            },
            decoy,
        };
        let reference = *plan
            .evidence_index
            .get(&item)
            .ok_or_else(|| invalid("internal peptide evidence index is incomplete"))?;
        out.indent(depth + 1)?;
        out.element()?;
        out.raw("<PeptideEvidenceRef")?;
        out.attr("peptideEvidence_ref", &format!("PEV_{reference}"))?;
        out.raw("/>\n")?;
    }
    if !hit.peak_annotations.is_empty() {
        write_fragmentation(out, cv, &hit.peak_annotations, false, depth + 1)?;
    }
    let consumed = write_score(
        out,
        cv,
        details,
        &identification.score_type,
        hit.score,
        depth + 1,
    )?;
    let mut skip: BTreeSet<String> = BTreeSet::new();
    skip.insert("calcMZ".into());
    skip.insert(user_param::TARGET_DECOY.into());
    if let Some(consumed) = consumed {
        skip.insert(consumed);
    }
    write_meta(out, cv, &hit.metadata, &skip, depth + 1)?;
    out.indent(depth)?;
    out.raw("</SpectrumIdentificationItem>\n")
}

/// Write one mzIdentML 1.3.0 document with default limits and registry.
///
/// # Errors
///
/// See [`write_with_registry`].
pub fn write(writer: impl Write, document: &MzIdentMLDocument) -> Result<()> {
    write_with_options(writer, document, &WriteOptions::default())
}

/// Write with explicit limits; see
/// [`write_with_registry`].
///
/// # Errors
///
/// See [`write_with_registry`].
pub fn write_with_options(
    writer: impl Write,
    document: &MzIdentMLDocument,
    options: &WriteOptions,
) -> Result<()> {
    write_with_registry(writer, document, options, ModificationsDB::global())
}

/// Write one mzIdentML 1.3.0 document, naming modifications from `registry`.
///
/// The whole document is built in memory and only then handed to `writer`, so a
/// refused document produces no output at all. Every XML id is positional
/// (`SIL_0`, `PEP_3`, `PEV_7`), which is the one deliberate difference from the
/// source's `UniqueIdGenerator` ids: two stores of the same document are
/// byte-identical here. The upstream test suite's `FuzzyDiff` whitelist exempts
/// `id=`, `href=`, `completion_time=` and `version=` for exactly that reason.
///
/// A run that declares `MS:1002494` or `is_cross_linking_experiment` selects
/// the cross-linking output path, as the source's `is_ppxl` does: every chain
/// of every hit becomes its own `Peptide` with the donor or acceptor cvParam
/// that pairs them, the identifications of one spectrum share a single
/// `SpectrumIdentificationResult`, and a labelled match writes its heavy half
/// as a second item.
///
/// # Errors
///
/// * [`Error::MissingInformation`] when the document has no identification run,
///   or when a peptide identification has no hits: the schema requires at least
///   one `SpectrumIdentification` per document and one
///   `SpectrumIdentificationItem` per `SpectrumIdentificationResult`, and the
///   source emits an empty, schema-invalid `SpectrumIdentificationResult`.
/// * [`Error::InvalidValue`] when a peptide identification names a run that is
///   not in the document, when two runs share an identifier, when a peptide
///   evidence names an accession no protein hit declares, or when a
///   modification name is unknown or ambiguous in `registry`. The source logs
///   and drops the affected element for the first three.
/// * [`Error::Parse`] when [`WriteOptions`] ceilings are exceeded - including
///   `max_records` while the plan is still being built, before any element is
///   emitted - and [`Error::Io`] if the writer fails.
/// * [`Error::InvalidValue`] also for a cross-linked hit that is missing the
///   metadata its own `xl_type` requires: the cross-linker name and mass, the
///   beta sequence or a link position, all of which the source reads with an
///   unchecked `getMetaValue`.
pub fn write_with_registry(
    writer: impl Write,
    document: &MzIdentMLDocument,
    options: &WriteOptions,
    registry: &ModificationsDB,
) -> Result<()> {
    preflight(document, options)?;
    let cv = ControlledVocabulary::psi_ms()?;
    let context = WriteContext {
        cv,
        details: cv.all_child_terms("MS:1001143")?,
        plan: plan(document, options, registry)?,
    };
    let plan = &context.plan;
    let mut out = Out::new(options);
    // The namespace carries only major.minor, e.g. "1.3" for "1.3.0".
    let short = match SCHEMA_VERSION.rsplit_once('.') {
        Some((head, _)) => head,
        None => SCHEMA_VERSION,
    };
    out.raw("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")?;
    out.element()?;
    out.raw("<MzIdentML xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n\txsi:schemaLocation=\"")?;
    out.raw(NAMESPACE_PREFIX)?;
    out.raw(short)?;
    out.raw(" https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/schema/mzIdentML")?;
    out.raw(SCHEMA_VERSION)?;
    out.raw(".xsd\"\n\txmlns=\"")?;
    out.raw(NAMESPACE_PREFIX)?;
    out.raw(short)?;
    out.raw("\"\n\tversion=\"")?;
    out.raw(SCHEMA_VERSION)?;
    out.raw("\"\n")?;
    let id = if document.document_id.is_empty() {
        "OpenMS"
    } else {
        &document.document_id
    };
    out.raw("\t")?;
    out.attr("id", id)?;
    if let Some(date) = &options.creation_date {
        out.raw("\n\t")?;
        out.attr("creationDate", date)?;
    }
    out.raw(">\n")?;
    out.raw("<cvList>\n")?;
    out.raw("\t<cv id=\"PSI-MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Vocabularies\" uri=\"http://purl.obolibrary.org/obo/ms/psi-ms.obo\" version=\"")?;
    out.text_content(cv.version())?;
    out.raw("\"></cv>\n")?;
    out.raw("\t<cv id=\"UNIMOD\" fullName=\"UNIMOD\" uri=\"http://www.unimod.org/obo/unimod.obo\"></cv>\n")?;
    out.raw("\t<cv id=\"UO\" fullName=\"UNIT-ONTOLOGY\" uri=\"https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo\"></cv>\n")?;
    if plan.crosslinking {
        out.raw("\t<cv id=\"XLMOD\" fullName=\"PSI cross-link modifications\" uri=\"https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/cv/XLMOD-1.0.0.obo\"></cv>\n")?;
    }
    out.raw("</cvList>\n")?;

    out.raw("<AnalysisSoftwareList>\n")?;
    for (index, (name, version)) in plan.software.iter().enumerate() {
        out.indent(1)?;
        out.element()?;
        out.raw("<AnalysisSoftware")?;
        out.attr("version", version)?;
        out.attr("name", name)?;
        out.attr("id", &format!("SOF_{index}"))?;
        out.raw(">\n")?;
        out.indent(2)?;
        out.raw("<SoftwareName>\n")?;
        out.indent(3)?;
        let term = software_term(cv, name);
        let term = if term.is_empty() { name.as_str() } else { term };
        cv_param_named(&mut out, cv, term, None)?;
        out.raw("\n")?;
        out.indent(2)?;
        out.raw("</SoftwareName>\n")?;
        out.indent(1)?;
        out.raw("</AnalysisSoftware>\n")?;
    }
    if !plan.software_index.contains_key("TOPP software") {
        out.indent(1)?;
        out.element()?;
        out.raw("<AnalysisSoftware")?;
        out.attr(
            "version",
            &format!("OpenMS TOPP v{}", env!("CARGO_PKG_VERSION")),
        )?;
        out.attr("name", "TOPP software")?;
        out.attr("id", "SOF_openms")?;
        out.raw(">\n")?;
        out.indent(2)?;
        out.raw("<SoftwareName>\n")?;
        out.indent(3)?;
        cv_param_named(&mut out, cv, "TOPP software", None)?;
        out.raw("\n")?;
        out.indent(2)?;
        out.raw("</SoftwareName>\n")?;
        out.indent(1)?;
        out.raw("</AnalysisSoftware>\n")?;
    }
    out.raw("</AnalysisSoftwareList>\n")?;

    out.raw("<SequenceCollection>\n")?;
    for (index, (accession, sequence, description, database)) in plan.proteins.iter().enumerate() {
        out.indent(1)?;
        out.element()?;
        out.raw("<DBSequence")?;
        out.attr("accession", accession)?;
        out.attr("searchDatabase_ref", &format!("SDB_{database}"))?;
        if !sequence.is_empty() {
            out.attr("length", &sequence.chars().count().to_string())?;
        }
        out.attr("id", &format!("PROT_{index}"))?;
        out.raw(">\n")?;
        if !sequence.is_empty() {
            out.indent(2)?;
            out.raw("<Seq>")?;
            out.text_content(sequence)?;
            out.raw("</Seq>\n")?;
        }
        out.indent(2)?;
        cv_param_named(
            &mut out,
            cv,
            "protein description",
            Some(&MetaValue::new(MetaValueData::String(description.clone()))?),
        )?;
        out.raw("\n")?;
        out.indent(1)?;
        out.raw("</DBSequence>\n")?;
    }
    for (index, sequence) in plan.peptides.iter().enumerate() {
        out.indent(1)?;
        out.element()?;
        out.raw("<Peptide")?;
        out.attr("id", &format!("PEP_{index}"))?;
        out.attr("name", &sequence.to_string())?;
        out.raw(">\n")?;
        out.indent(2)?;
        out.raw("<PeptideSequence>")?;
        out.text_content(sequence.as_str())?;
        out.raw("</PeptideSequence>\n")?;
        write_modifications(&mut out, sequence, 2)?;
        out.indent(1)?;
        out.raw("</Peptide>\n")?;
    }
    for (index, peptide) in plan.xl_peptides.iter().enumerate() {
        write_crosslink_peptide(&mut out, cv, peptide, index, 1)?;
    }
    for (index, (_, evidence)) in plan.xl_evidences.iter().enumerate() {
        write_evidence(&mut out, evidence, index, 1)?;
    }
    for (index, evidence) in plan.evidences.iter().enumerate() {
        write_evidence(&mut out, evidence, index, 1)?;
    }
    out.raw("</SequenceCollection>\n")?;

    out.raw("<AnalysisCollection>\n")?;
    for (index, run) in document.protein_identifications.iter().enumerate() {
        let spectra = spectra_reference(plan, run);
        let database = *plan
            .database_index
            .get(&run.search_parameters.database)
            .ok_or_else(|| invalid("internal search database index is incomplete"))?;
        out.indent(1)?;
        out.element()?;
        out.raw("<SpectrumIdentification")?;
        out.attr("id", &format!("SI_{index}"))?;
        out.attr(
            "spectrumIdentificationProtocol_ref",
            &format!("SIP_{index}"),
        )?;
        out.attr("spectrumIdentificationList_ref", &format!("SIL_{index}"))?;
        if let Some(date) = &run.date_time {
            out.attr("activityDate", date)?;
        }
        out.raw(">\n")?;
        out.indent(2)?;
        out.element()?;
        out.raw("<InputSpectra")?;
        out.attr("spectraData_ref", &format!("SDAT_{spectra}"))?;
        out.raw("/>\n")?;
        out.indent(2)?;
        out.element()?;
        out.raw("<SearchDatabaseRef")?;
        out.attr("searchDatabase_ref", &format!("SDB_{database}"))?;
        out.raw("/>\n")?;
        out.indent(1)?;
        out.raw("</SpectrumIdentification>\n")?;
    }
    out.raw("</AnalysisCollection>\n")?;

    out.raw("<AnalysisProtocolCollection>\n")?;
    for (index, run) in document.protein_identifications.iter().enumerate() {
        let software = *plan
            .software_index
            .get(&run.search_engine)
            .ok_or_else(|| invalid("internal software index is incomplete"))?;
        write_protocol(&mut out, cv, run, index, software, registry)?;
    }
    out.raw("</AnalysisProtocolCollection>\n")?;

    out.raw("<DataCollection>\n\t<Inputs>\n")?;
    for (index, (location, version)) in plan.databases.iter().enumerate() {
        out.indent(2)?;
        out.element()?;
        out.raw("<SearchDatabase")?;
        out.attr("location", location)?;
        if !version.is_empty() {
            out.attr("version", version)?;
        }
        out.attr("id", &format!("SDB_{index}"))?;
        out.raw(">\n")?;
        out.indent(3)?;
        out.raw("<FileFormat>\n")?;
        out.indent(4)?;
        cv_param_named(&mut out, cv, "FASTA format", None)?;
        out.raw("\n")?;
        out.indent(3)?;
        out.raw("</FileFormat>\n")?;
        out.indent(3)?;
        out.raw("<DatabaseName>\n")?;
        out.indent(4)?;
        out.element()?;
        out.raw("<userParam")?;
        out.attr("name", location)?;
        out.raw("/>\n")?;
        out.indent(3)?;
        out.raw("</DatabaseName>\n")?;
        for run in &document.protein_identifications {
            if run.search_parameters.database != *location {
                continue;
            }
            if let Some(value) = run.search_parameters.metadata.get("MS:1001029") {
                out.indent(3)?;
                cv_param(&mut out, cv.get_term("MS:1001029")?, Some(value))?;
                out.raw("\n")?;
                break;
            }
        }
        out.indent(2)?;
        out.raw("</SearchDatabase>\n")?;
    }
    for (index, location) in plan.spectra.iter().enumerate() {
        let (format, id_format) = spectra_format(location);
        out.indent(2)?;
        out.element()?;
        out.raw("<SpectraData")?;
        out.attr("location", location)?;
        out.attr("id", &format!("SDAT_{index}"))?;
        out.raw(">\n")?;
        out.indent(3)?;
        out.raw("<FileFormat>\n")?;
        out.indent(4)?;
        cv_param_named(&mut out, cv, format, None)?;
        out.raw("\n")?;
        out.indent(3)?;
        out.raw("</FileFormat>\n")?;
        out.indent(3)?;
        out.raw("<SpectrumIDFormat>\n")?;
        out.indent(4)?;
        cv_param_named(&mut out, cv, id_format, None)?;
        out.raw("\n")?;
        out.indent(3)?;
        out.raw("</SpectrumIDFormat>\n")?;
        out.indent(2)?;
        out.raw("</SpectraData>\n")?;
    }
    out.raw("\t</Inputs>\n\t<AnalysisData>\n")?;
    let mut item = 0usize;
    let mut result = 0usize;
    for (index, run) in document.protein_identifications.iter().enumerate() {
        out.indent(2)?;
        out.element()?;
        out.raw("<SpectrumIdentificationList")?;
        out.attr("id", &format!("SIL_{index}"))?;
        out.raw(">\n")?;
        out.raw("\t\t\t<FragmentationTable>\n")?;
        for (id, name) in [
            ("Measure_mz", "product ion m/z"),
            ("Measure_int", "product ion intensity"),
            ("Measure_error", "product ion m/z error"),
        ] {
            out.indent(4)?;
            out.element()?;
            out.raw("<Measure")?;
            out.attr("id", id)?;
            out.raw(">\n")?;
            out.indent(5)?;
            cv_param_named(&mut out, cv, name, None)?;
            out.raw("\n")?;
            out.indent(4)?;
            out.raw("</Measure>\n")?;
        }
        out.raw("\t\t\t</FragmentationTable>\n")?;
        let spectra = spectra_reference(plan, run);
        if plan.crosslinking {
            // One result per spectrum reference, in the source's order, with
            // every chain and label of that spectrum's matches inside it. The
            // source appends them all to the first list; each stays with its
            // own run here, so the run link a load established survives.
            for entry in plan.xl_results.iter().filter(|entry| entry.run == index) {
                out.indent(3)?;
                out.element()?;
                out.raw("<SpectrumIdentificationResult")?;
                out.attr("spectraData_ref", &format!("SDAT_{spectra}"))?;
                out.attr("spectrumID", &entry.spectrum_id)?;
                out.attr("id", &format!("SIR_{result}"))?;
                result = result.saturating_add(1);
                out.raw(">\n")?;
                for entry in &entry.items {
                    write_crosslink_item(&mut out, &context, document, entry, item, 4)?;
                    item = item.saturating_add(1);
                }
                out.indent(3)?;
                out.raw("</SpectrumIdentificationResult>\n")?;
            }
            out.indent(2)?;
            out.raw("</SpectrumIdentificationList>\n")?;
            continue;
        }
        for identification in &document.peptide_identifications {
            if identification.identifier != run.identifier {
                continue;
            }
            out.indent(3)?;
            out.element()?;
            out.raw("<SpectrumIdentificationResult")?;
            out.attr("spectraData_ref", &format!("SDAT_{spectra}"))?;
            let reference = identification.spectrum_reference();
            let reference = if reference.is_empty() {
                // Source falls back to "MZ:<mz>@RT:<rt>" with "nan" for a
                // missing coordinate; NaN is not a valid xsd:double there
                // either, so this writes the schema's NaN spelling.
                let mz = identification
                    .mz
                    .map(number_text)
                    .transpose()?
                    .unwrap_or_else(|| "NaN".to_owned());
                let rt = identification
                    .rt
                    .map(number_text)
                    .transpose()?
                    .unwrap_or_else(|| "NaN".to_owned());
                format!("MZ:{mz}@RT:{rt}")
            } else {
                reference
            };
            out.attr("spectrumID", &reference)?;
            out.attr("id", &format!("SIR_{result}"))?;
            result = result.saturating_add(1);
            out.raw(">\n")?;
            for hit in &identification.hits {
                write_item(
                    &mut out,
                    &context,
                    identification,
                    hit,
                    item,
                    run.significance_threshold,
                    4,
                )?;
                item = item.saturating_add(1);
            }
            if let Some(rt) = identification.rt {
                out.indent(4)?;
                let value = MetaValue::new(MetaValueData::Float(rt))?.with_unit(Unit::new(
                    "UO:0000010",
                    "second",
                    "UO",
                )?)?;
                cv_param_named(&mut out, cv, "retention time", Some(&value))?;
                out.raw("\n")?;
            }
            out.indent(3)?;
            out.raw("</SpectrumIdentificationResult>\n")?;
        }
        out.indent(2)?;
        out.raw("</SpectrumIdentificationList>\n")?;
    }
    out.raw("\t</AnalysisData>\n</DataCollection>\n</MzIdentML>\n")?;
    let mut writer = writer;
    writer.write_all(out.text.as_bytes())?;
    writer.flush()?;
    Ok(())
}

/// One `PeptideEvidence` element.
///
/// mzIdentML counts residues from 1 and OpenMS from 0, so both positions are
/// written one higher than they are stored.
fn write_evidence(out: &mut Out, evidence: &EvidenceOut, index: usize, depth: usize) -> Result<()> {
    out.indent(depth)?;
    out.element()?;
    out.raw("<PeptideEvidence")?;
    out.attr("id", &format!("PEV_{index}"))?;
    out.attr("peptide_ref", &format!("PEP_{}", evidence.peptide))?;
    out.attr("dBSequence_ref", &format!("PROT_{}", evidence.protein))?;
    if let Some(post) = evidence.post {
        out.attr("post", &terminal_text(post))?;
    }
    if let Some(pre) = evidence.pre {
        out.attr("pre", &terminal_text(pre))?;
    }
    if let Some(start) = evidence.start {
        out.attr("start", &start.saturating_add(1).to_string())?;
    }
    if let Some(end) = evidence.end {
        out.attr("end", &end.saturating_add(1).to_string())?;
    }
    if let Some(decoy) = evidence.decoy {
        out.attr("isDecoy", if decoy { "true" } else { "false" })?;
    }
    out.raw("/>\n")
}

/// `-` for the terminal markers, the residue letter otherwise.
fn terminal_text(code: char) -> String {
    match code {
        '[' | ']' => "-".to_owned(),
        other => other.to_string(),
    }
}

fn spectra_reference(plan: &Plan, run: &ProteinIdentification) -> usize {
    let location = match run.metadata.get("spectra_data") {
        Some(value) => value
            .as_string_list()
            .ok()
            .and_then(|list| list.first().cloned())
            .unwrap_or_else(|| text_of(value)),
        None => String::new(),
    };
    let location = if location.is_empty() {
        "UNKNOWN".to_owned()
    } else {
        trim_file_uri(&location)
    };
    plan.spectra_index.get(&location).copied().unwrap_or(0)
}

/// The `FileFormat` and `SpectrumIDFormat` term names for a spectra location.
///
/// The source keeps a four-entry map and defaults to mzML for anything else,
/// including its own "UNKNOWN" placeholder.
fn spectra_format(location: &str) -> (&'static str, &'static str) {
    let lowered = location.to_ascii_lowercase();
    if lowered.ends_with(".mzxml") {
        ("ISB mzXML format", "scan number only nativeID format")
    } else if lowered.ends_with(".mzdata") {
        ("PSI mzData format", "spectrum identifier nativeID format")
    } else if lowered.ends_with(".mgf") {
        ("Mascot MGF format", "multiple peak list nativeID format")
    } else {
        ("mzML format", "mzML unique identifier")
    }
}

/// The `Modification` elements of one peptide.
///
/// The C-terminal modification is written at `length + 1`, which is what the
/// mzIdentML schema defines. The source writes `length`, which its own reader
/// then interprets as a modification on the last residue, so a store/load cycle
/// silently moves it (see `OpenMS_CPP_ISSUES.md`).
fn write_modifications(out: &mut Out, sequence: &AASequence, depth: usize) -> Result<()> {
    let length = sequence.len();
    if let Some(modification) = sequence.n_terminal_modification() {
        write_modification(out, modification, 0, None, depth)?;
    }
    if let Some(modification) = sequence.c_terminal_modification() {
        write_modification(out, modification, length.saturating_add(1), None, depth)?;
    }
    for index in 0..length {
        if let Some(modification) = sequence.residue_modification(index)? {
            let residue = sequence.as_str().chars().nth(index);
            write_modification(out, modification, index.saturating_add(1), residue, depth)?;
        }
    }
    Ok(())
}

fn write_modification(
    out: &mut Out,
    modification: &crate::chemistry::SequenceModification,
    location: usize,
    residue: Option<char>,
    depth: usize,
) -> Result<()> {
    let accession = modification
        .known()
        .and_then(|record| record.unimod_accession())
        .map(|accession| accession.replace("UniMod:", "UNIMOD:"));
    out.indent(depth)?;
    out.element()?;
    out.raw("<Modification")?;
    out.attr("location", &location.to_string())?;
    if let Some(residue) = residue {
        out.attr("residues", &residue.to_string())?;
    }
    if accession.is_none() {
        // An anonymous or non-UniMod modification keeps only its mass.
        out.attr(
            "monoisotopicMassDelta",
            &number_text(modification.diff_mono_mass()?)?,
        )?;
    }
    out.raw(">\n")?;
    out.indent(depth + 1)?;
    out.element()?;
    match accession {
        Some(accession) => {
            out.raw("<cvParam")?;
            out.attr("accession", &accession)?;
            out.attr("cvRef", "UNIMOD")?;
            out.attr("name", modification.name())?;
            out.raw("/>\n")?;
        }
        None => {
            out.raw(
                "<cvParam cvRef=\"MS\" accession=\"MS:1001460\" name=\"unknown modification\"/>\n",
            )?;
        }
    }
    out.indent(depth)?;
    out.raw("</Modification>\n")
}

/// One `SpectrumIdentificationProtocol`.
fn write_protocol(
    out: &mut Out,
    cv: &ControlledVocabulary,
    run: &ProteinIdentification,
    index: usize,
    software: usize,
    registry: &ModificationsDB,
) -> Result<()> {
    let parameters = &run.search_parameters;
    out.indent(1)?;
    out.element()?;
    out.raw("<SpectrumIdentificationProtocol")?;
    out.attr("id", &format!("SIP_{index}"))?;
    out.attr("analysisSoftware_ref", &format!("SOF_{software}"))?;
    out.raw(">\n")?;
    out.indent(2)?;
    out.raw("<SearchType>\n")?;
    out.indent(3)?;
    cv_param_named(out, cv, "ms-ms search", None)?;
    out.raw("\n")?;
    out.indent(2)?;
    out.raw("</SearchType>\n")?;
    out.indent(2)?;
    out.raw("<AdditionalSearchParams>\n")?;
    let mut skip = BTreeSet::new();
    // Written into SearchDatabase instead, as the SearchDatabase_may rule wants.
    skip.insert("MS:1001029".to_owned());
    write_meta(out, cv, &parameters.metadata, &skip, 3)?;
    for (name, value, kind) in [
        ("charges", parameters.charges.clone(), "xsd:string"),
        ("taxonomy", parameters.taxonomy.clone(), "xsd:string"),
    ] {
        if name == "taxonomy" && value.is_empty() {
            continue;
        }
        out.indent(3)?;
        out.element()?;
        out.raw("<userParam")?;
        out.attr("name", name)?;
        out.attr("type", kind)?;
        out.attr("value", &value)?;
        out.raw("/>\n")?;
    }
    // NumTolerableTermini is consumed on read and dropped by the source's
    // writer, which loses the enzyme specificity; it is written back here.
    let termini = match parameters.enzyme_specificity {
        EnzymeTermSpecificity::None => Some(0),
        EnzymeTermSpecificity::Semi => Some(1),
        EnzymeTermSpecificity::Full => Some(2),
        EnzymeTermSpecificity::Unknown => None,
    };
    if let Some(termini) = termini {
        out.indent(3)?;
        out.element()?;
        out.raw("<userParam")?;
        out.attr("name", "NumTolerableTermini")?;
        out.attr("type", "xsd:integer")?;
        out.attr("value", &termini.to_string())?;
        out.raw("/>\n")?;
    }
    out.indent(2)?;
    out.raw("</AdditionalSearchParams>\n")?;
    if !parameters.fixed_modifications.is_empty() || !parameters.variable_modifications.is_empty() {
        out.indent(2)?;
        out.raw("<ModificationParams>\n")?;
        write_mod_params(out, cv, &parameters.fixed_modifications, true, registry, 3)?;
        write_mod_params(
            out,
            cv,
            &parameters.variable_modifications,
            false,
            registry,
            3,
        )?;
        out.indent(2)?;
        out.raw("</ModificationParams>\n")?;
    }
    write_enzyme(
        out,
        cv,
        &parameters.digestion_enzyme,
        parameters.missed_cleavages,
        2,
    )?;
    for (tag, tolerance) in [
        ("FragmentTolerance", parameters.fragment_tolerance),
        ("ParentTolerance", parameters.precursor_tolerance),
    ] {
        let (value, ppm) = match tolerance {
            Tolerance::Absolute(value) => (value, false),
            Tolerance::Ppm(value) => (value, true),
        };
        let unit = if ppm {
            Unit::new("UO:0000169", "parts per million", "UO")?
        } else {
            Unit::new("UO:0000221", "dalton", "UO")?
        };
        let value = MetaValue::new(MetaValueData::Float(value))?.with_unit(unit)?;
        out.indent(2)?;
        out.raw("<")?;
        out.raw(tag)?;
        out.raw(">\n")?;
        for name in [
            "search tolerance plus value",
            "search tolerance minus value",
        ] {
            out.indent(3)?;
            cv_param_named(out, cv, name, Some(&value))?;
            out.raw("\n")?;
        }
        out.indent(2)?;
        out.raw("</")?;
        out.raw(tag)?;
        out.raw(">\n")?;
    }
    out.indent(2)?;
    out.raw("<Threshold>\n")?;
    out.indent(3)?;
    if run.significance_threshold == 0.0 {
        cv_param_named(out, cv, "no threshold", None)?;
    } else {
        let value = MetaValue::new(MetaValueData::Float(run.significance_threshold))?;
        cv_param_named(out, cv, "PSM-level statistical threshold", Some(&value))?;
    }
    out.raw("\n")?;
    out.indent(2)?;
    out.raw("</Threshold>\n")?;
    out.indent(1)?;
    out.raw("</SpectrumIdentificationProtocol>\n")
}

/// Store one document as plain mzIdentML with default limits and registry.
///
/// # Errors
///
/// See [`store_with_registry`].
pub fn store(path: impl AsRef<Path>, document: &MzIdentMLDocument) -> Result<()> {
    store_with_options(path, document, &WriteOptions::default())
}

/// Store with explicit limits; see
/// [`store_with_registry`].
///
/// # Errors
///
/// See [`store_with_registry`].
pub fn store_with_options(
    path: impl AsRef<Path>,
    document: &MzIdentMLDocument,
    options: &WriteOptions,
) -> Result<()> {
    store_with_registry(path, document, options, ModificationsDB::global())
}

/// Validate, then publish the document atomically.
///
/// The filename must carry the `.mzid` extension, as `MzIdentMLFile::store`
/// requires through `FileHandler::hasValidExtension`. A failure leaves an
/// existing destination untouched and removes the temporary output.
///
/// # Errors
///
/// [`Error::InvalidValue`] for a filename that is not UTF-8 or does not end in
/// `.mzid`, plus everything
/// [`write_with_registry`]
/// reports.
pub fn store_with_registry(
    path: impl AsRef<Path>,
    document: &MzIdentMLDocument,
    options: &WriteOptions,
    registry: &ModificationsDB,
) -> Result<()> {
    let path = path.as_ref();
    let filename = path
        .to_str()
        .ok_or_else(|| invalid("output filename must be UTF-8"))?;
    if !super::file_types::has_valid_extension(filename, super::FileType::MzIdentMl) {
        return Err(invalid(
            "invalid mzIdentML output extension, expected '.mzid'",
        ));
    }
    super::path_io::write_plain(path, |writer| {
        write_with_registry(writer, document, options, registry)
    })
}
