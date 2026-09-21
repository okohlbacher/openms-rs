// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! XSD validation against the schemas the source bundles, or against a
//! caller's own schema: `VALIDATORS/XMLValidator.h` and the `isValid` that
//! every `Internal::XMLFile` inherits. See `docs/XML_SCHEMA_SUPPORT.md`.
//!
//! The source validates with Xerces-C; this port uses libxml2 through the exact
//! `libxml` binding, behind the default-off `xml-schema` feature. The schemas
//! are the ones the source ships in `share/OpenMS/SCHEMAS`, byte for byte, and
//! are compiled into the crate, so no schema is looked up at run time the way
//! the source's `File::find` does.
//!
//! Native input/preflight limits and post-engine result limits do not bound the
//! C engine's DOM, identity tables, runtime, or pre-return diagnostic
//! allocations.

use super::cv_xml::{Attributes, Element, Meter, XmlLimits, document, scan_elements};
use crate::{Error, Result};
use libxml::{
    error::{StructuredError, XmlErrorLevel},
    parser::{Parser, ParserOptions},
    schemas::{SchemaParserContext, SchemaValidationContext},
};
use std::{
    borrow::Cow,
    io::BufRead,
    panic::{AssertUnwindSafe, catch_unwind},
    path::Path,
    sync::{Mutex, OnceLock, PoisonError},
};

/// The schema a document was validated against.
///
/// Every variant but [`External`](Self::External) names one schema the source
/// ships in `share/OpenMS/SCHEMAS` and registers with an `Internal::XMLFile`
/// constructor; the crate carries its bytes unchanged.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum SchemaKind {
    /// `mzML_1_10.xsd`, registered by `MzMLFile` and `ImzMLFile`.
    MzML,
    /// `mzML_idx_1_10.xsd`, which `MzMLFile::isValid` selects for an indexed document.
    IndexedMzML,
    /// `FeatureXML_1_9.xsd`, registered by `FeatureXMLFile`.
    FeatureXML,
    /// `ConsensusXML_1_7.xsd`, registered by `ConsensusXMLFile`.
    ConsensusXML,
    /// `IdXML_1_5.xsd`, registered by `IdXMLFile`.
    IdXML,
    /// `Param_1_8_0.xsd`, registered by `ParamXMLFile`.
    ParamXML,
    /// `TrafoXML_1_1.xsd`, registered by `TransformationXMLFile`.
    TransformationXML,
    /// `mzData_1_05.xsd`, registered by `MzDataFile`.
    MzData,
    /// `mzXML_idx_3.1.xsd`, registered by `MzXMLFile` for every mzXML document,
    /// indexed or not. It includes `mzXML_3.1_mod.xsd`, which includes
    /// `separation_technique_1.0.xsd` and `general_types_1.0.xsd`.
    MzXML,
    /// `mzIdentML1.0.0.xsd`, which includes `FuGElightv1.0.0.xsd`.
    MzIdentML1_0_0,
    /// `mzIdentML1.1.0.xsd`.
    MzIdentML1_1_0,
    /// `mzIdentML1.2.0.xsd`.
    MzIdentML1_2_0,
    /// `mzIdentML1.3.0.xsd`, the schema `MzIdentMLFile` registers by default.
    MzIdentML1_3_0,
    /// A caller's schema, as the source `XMLValidator::isValid(filename, schema, os)`
    /// takes one; see [`validate_against`].
    External,
}

/// Where a bundled schema comes from and which `XMLFile` version it carries.
struct Bundle {
    /// The source's `schema_location_`: `/SCHEMAS/` and a file name under
    /// `share/OpenMS/SCHEMAS`, which is also its name under `resources/schemas`.
    location: &'static str,
    /// The `schema_version_` its source `XMLFile` constructor passes.
    version: &'static str,
    /// Whether the schema reaches other bundled files through `xs:include`.
    includes: bool,
}

impl Bundle {
    /// The file name, without the `/SCHEMAS/` the location starts with.
    fn file(&self) -> &'static str {
        self.location.trim_start_matches("/SCHEMAS/")
    }
}

impl SchemaKind {
    fn bundle(self) -> Option<Bundle> {
        let (location, version, includes) = match self {
            Self::MzML => ("/SCHEMAS/mzML_1_10.xsd", "1.1.0", false),
            Self::IndexedMzML => ("/SCHEMAS/mzML_idx_1_10.xsd", "1.1.0", false),
            Self::FeatureXML => ("/SCHEMAS/FeatureXML_1_9.xsd", "1.9", false),
            Self::ConsensusXML => ("/SCHEMAS/ConsensusXML_1_7.xsd", "1.7", false),
            Self::IdXML => ("/SCHEMAS/IdXML_1_5.xsd", "1.5", false),
            Self::ParamXML => ("/SCHEMAS/Param_1_8_0.xsd", "1.8.0", false),
            Self::TransformationXML => ("/SCHEMAS/TrafoXML_1_1.xsd", "1.1", false),
            Self::MzData => ("/SCHEMAS/mzData_1_05.xsd", "1.05", false),
            Self::MzXML => ("/SCHEMAS/mzXML_idx_3.1.xsd", "3.1", true),
            Self::MzIdentML1_0_0 => ("/SCHEMAS/mzIdentML1.0.0.xsd", "1.0.0", true),
            Self::MzIdentML1_1_0 => ("/SCHEMAS/mzIdentML1.1.0.xsd", "1.1.0", false),
            Self::MzIdentML1_2_0 => ("/SCHEMAS/mzIdentML1.2.0.xsd", "1.2.0", false),
            Self::MzIdentML1_3_0 => ("/SCHEMAS/mzIdentML1.3.0.xsd", "1.3.0", false),
            Self::External => return None,
        };
        Some(Bundle {
            location,
            version,
            includes,
        })
    }

    /// The source's schema location, as its `XMLFile` constructor registers it
    /// (`"/SCHEMAS/FeatureXML_1_9.xsd"`); `None` for [`External`](Self::External).
    ///
    /// The source resolves this through `File::find` when `isValid` runs; here
    /// the named file is compiled into the crate, so nothing is looked up.
    pub fn location(self) -> Option<&'static str> {
        self.bundle().map(|b| b.location)
    }

    /// The schema version string, the source `XMLFile::getVersion()` of the
    /// class that registers this schema; `None` for [`External`](Self::External).
    pub fn version(self) -> Option<&'static str> {
        self.bundle().map(|b| b.version)
    }

    /// The bundled mzIdentML schema for a version string, as
    /// `MzIdentMLFile::isValid` maps `detectVersion`'s answer to
    /// `"/SCHEMAS/mzIdentML" + version + ".xsd"`.
    ///
    /// The source falls back to its default 1.3.0 schema when no such file is
    /// shipped; this returns `None` instead and leaves the choice to the caller.
    pub fn mzidentml(version: &str) -> Option<Self> {
        match version {
            "1.0.0" => Some(Self::MzIdentML1_0_0),
            "1.1.0" => Some(Self::MzIdentML1_1_0),
            "1.2.0" => Some(Self::MzIdentML1_2_0),
            "1.3.0" => Some(Self::MzIdentML1_3_0),
            _ => None,
        }
    }
}

/// Every bundled schema file, by its upstream name. Static buffers stay alive
/// for the process, which the libxml2 memory parser needs.
const FILES: [(&str, &[u8]); 17] = [
    (
        "mzML_1_10.xsd",
        include_bytes!("../../resources/schemas/mzML_1_10.xsd"),
    ),
    (
        "mzML_idx_1_10.xsd",
        include_bytes!("../../resources/schemas/mzML_idx_1_10.xsd"),
    ),
    (
        "FeatureXML_1_9.xsd",
        include_bytes!("../../resources/schemas/FeatureXML_1_9.xsd"),
    ),
    (
        "ConsensusXML_1_7.xsd",
        include_bytes!("../../resources/schemas/ConsensusXML_1_7.xsd"),
    ),
    (
        "IdXML_1_5.xsd",
        include_bytes!("../../resources/schemas/IdXML_1_5.xsd"),
    ),
    (
        "Param_1_8_0.xsd",
        include_bytes!("../../resources/schemas/Param_1_8_0.xsd"),
    ),
    (
        "TrafoXML_1_1.xsd",
        include_bytes!("../../resources/schemas/TrafoXML_1_1.xsd"),
    ),
    (
        "mzData_1_05.xsd",
        include_bytes!("../../resources/schemas/mzData_1_05.xsd"),
    ),
    (
        "mzXML_idx_3.1.xsd",
        include_bytes!("../../resources/schemas/mzXML_idx_3.1.xsd"),
    ),
    (
        "mzXML_3.1_mod.xsd",
        include_bytes!("../../resources/schemas/mzXML_3.1_mod.xsd"),
    ),
    (
        "separation_technique_1.0.xsd",
        include_bytes!("../../resources/schemas/separation_technique_1.0.xsd"),
    ),
    (
        "general_types_1.0.xsd",
        include_bytes!("../../resources/schemas/general_types_1.0.xsd"),
    ),
    (
        "mzIdentML1.0.0.xsd",
        include_bytes!("../../resources/schemas/mzIdentML1.0.0.xsd"),
    ),
    (
        "FuGElightv1.0.0.xsd",
        include_bytes!("../../resources/schemas/FuGElightv1.0.0.xsd"),
    ),
    (
        "mzIdentML1.1.0.xsd",
        include_bytes!("../../resources/schemas/mzIdentML1.1.0.xsd"),
    ),
    (
        "mzIdentML1.2.0.xsd",
        include_bytes!("../../resources/schemas/mzIdentML1.2.0.xsd"),
    ),
    (
        "mzIdentML1.3.0.xsd",
        include_bytes!("../../resources/schemas/mzIdentML1.3.0.xsd"),
    ),
];

fn resource(name: &str) -> Option<&'static [u8]> {
    FILES
        .iter()
        .find(|(file, _)| *file == name)
        .map(|(_, b)| *b)
}

/// Structured libxml2 severity. Informational messages do not invalidate a file.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SchemaDiagnosticLevel {
    /// libxml2's `XML_ERR_NONE`: reported, but the file stays valid.
    Information,
    /// A warning; the source's `XMLValidator` counts it as invalid too.
    Warning,
    /// A recoverable error.
    Error,
    /// A fatal error.
    Fatal,
}
/// Owned engine diagnostic. Text/codes depend on the installed libxml2 version.
///
/// This is what the source's `XMLValidator::logError_` formats into one line
/// on its `std::ostream`: the message, and the line and column it occurred at.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SchemaDiagnostic {
    /// The severity libxml2 reported.
    pub level: SchemaDiagnosticLevel,
    /// libxml2's message text.
    pub message: String,
    /// The document or schema the message is about, when libxml2 names one.
    pub filename: Option<String>,
    /// One-based line, when libxml2 reports one.
    pub line: Option<usize>,
    /// One-based column, when libxml2 reports one.
    pub column: Option<usize>,
    /// libxml2's error domain.
    pub domain: i32,
    /// libxml2's error code.
    pub code: i32,
}
/// Schema validity does not establish CV semantics, checksums or index integrity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SchemaValidationReport {
    /// The schema the document was validated against.
    pub schema: SchemaKind,
    /// Every diagnostic, schema compilation's first and then validation's.
    pub diagnostics: Vec<SchemaDiagnostic>,
    engine_valid: bool,
}
impl SchemaValidationReport {
    /// Like source XMLValidator, any warning/error/fatal makes validity false.
    pub fn is_valid(&self) -> bool {
        self.engine_valid
            && self
                .diagnostics
                .iter()
                .all(|d| d.level == SchemaDiagnosticLevel::Information)
    }
}
/// Configurable limits on native preflight and returned results.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SchemaValidationLimits {
    /// Maximum encoded input bytes and decoded/canonical UTF-8 bytes, for the
    /// document and, with [`validate_against`], for the caller's schema too.
    pub max_xml_bytes: usize,
    /// Maximum element nesting depth.
    pub max_depth: usize,
    /// Maximum number of elements.
    pub max_elements: usize,
    /// Shared Rust-side decoding, lexical, namespace and result-copy work only.
    pub max_work: usize,
    /// Shared logical Rust-owned preflight/result allocation, excluding engine internals.
    pub max_bytes: usize,
    /// Checked after the binding has accumulated its uncapped diagnostic vector.
    pub max_diagnostics: usize,
    /// Message/filename bytes checked after engine collection, before public report.
    pub max_diagnostic_bytes: usize,
}
impl Default for SchemaValidationLimits {
    fn default() -> Self {
        Self {
            max_xml_bytes: 16 * 1024 * 1024,
            max_depth: 128,
            max_elements: 1_000_000,
            max_work: 50_000_000,
            max_bytes: 128 * 1024 * 1024,
            max_diagnostics: 10_000,
            max_diagnostic_bytes: 1024 * 1024,
        }
    }
}
/// Validation settings. The bundled schemas are fixed; no resolver is accepted.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SchemaValidationOptions {
    /// Native preflight and result limits.
    pub limits: SchemaValidationLimits,
}

/// Validate a plain/gzip/bzip2 file against a bundled schema, with the default
/// native limits.
///
/// Source `Internal::XMLFile::isValid(filename, os)`, which validates against
/// the one schema its constructor registered, whatever the document's root
/// is: a consensusXML file checked with [`SchemaKind::FeatureXML`] is a report
/// whose [`is_valid`](SchemaValidationReport::is_valid) is false, as the
/// source's `false`. The messages the source writes to `os` are the report's
/// diagnostics.
///
/// # Errors
///
/// [`Error::InvalidValue`] for [`SchemaKind::External`], which names no
/// bundled schema; [`Error::Io`] when the file cannot be read, which is where
/// the source throws `Exception::FileNotFound`; and, as for
/// [`validate_reader`], an error for input that is not well-formed XML, where
/// the source returns `false`.
pub fn validate(schema: SchemaKind, path: impl AsRef<Path>) -> Result<SchemaValidationReport> {
    validate_with_options(schema, path, &SchemaValidationOptions::default())
}

/// [`validate`] with explicit native limits. Compression is detected by bytes.
///
/// # Errors
///
/// As [`validate`], plus a limit error when a limit in `options` is exceeded.
pub fn validate_with_options(
    schema: SchemaKind,
    path: impl AsRef<Path>,
    options: &SchemaValidationOptions,
) -> Result<SchemaValidationReport> {
    validate_reader(schema, super::path_io::open(path.as_ref())?, options)
}

/// Validate caller-supplied uncompressed XML against a bundled schema, after
/// full bounded decoding and lexical preflight.
///
/// Semantic XSD violations are `Ok(report)`; malformed or unsupported input,
/// I/O and limits are `Err`.
///
/// # Errors
///
/// [`Error::InvalidValue`] for [`SchemaKind::External`]; [`Error::Parse`] for
/// XML that is not well-formed, including namespace errors libxml2 would
/// tolerate, where the source reports a fatal error and returns `false`;
/// [`Error::Unsupported`] for a DTD, an unsupported encoding or non-ASCII
/// Latin-1; [`Error::Io`] when reading fails; and a limit error when a limit in
/// `options` is exceeded.
pub fn validate_reader(
    schema: SchemaKind,
    input: impl BufRead,
    options: &SchemaValidationOptions,
) -> Result<SchemaValidationReport> {
    if schema.bundle().is_none() {
        return Err(Error::InvalidValue(
            "SchemaKind::External names no bundled schema; use validate_against".into(),
        ));
    }
    let o = &options.limits;
    let mut m = Meter::new(o.max_work, o.max_bytes, "XML schema validation");
    run(input, o, &mut m, "XML", |_, _, _| {
        Ok(Grammar::Bundled(schema))
    })
}

/// Validate a plain/gzip/bzip2 file against the caller's own schema file, with
/// the default native limits.
///
/// Source `XMLValidator::isValid(filename, schema, os)`. The schema must be
/// self-contained: an `xs:include`, `xs:import`, `xs:redefine` or
/// `xs:override` is refused before anything is compiled, because resolving it
/// would read other files or the network on the caller's behalf, which the
/// source leaves to Xerces' default resolver. A DTD in the schema is refused
/// as it is in the document.
///
/// # Errors
///
/// [`Error::Io`] when either file cannot be read, which is where the source
/// throws `Exception::FileNotFound` for the document; [`Error::Unsupported`]
/// for a schema that is not self-contained; [`Error::InvalidValue`] when
/// libxml2 cannot compile the schema, where the source reports the schema's
/// errors and returns `false`; and, as for [`validate_reader`], an error for a
/// document that is not well-formed.
pub fn validate_against(
    path: impl AsRef<Path>,
    schema: impl AsRef<Path>,
) -> Result<SchemaValidationReport> {
    validate_against_with_options(path, schema, &SchemaValidationOptions::default())
}

/// [`validate_against`] with explicit native limits. Compression of the
/// document is detected by bytes; the schema must be plain XML.
///
/// # Errors
///
/// As [`validate_against`], plus a limit error when a limit in `options` is
/// exceeded. The document is opened first, as the source checks it first.
pub fn validate_against_with_options(
    path: impl AsRef<Path>,
    schema: impl AsRef<Path>,
    options: &SchemaValidationOptions,
) -> Result<SchemaValidationReport> {
    let document = super::path_io::open(path.as_ref())?;
    let schema = std::io::BufReader::new(std::fs::File::open(schema.as_ref())?);
    validate_reader_against(document, schema, options)
}

/// Validate caller-supplied uncompressed XML against caller-supplied schema
/// text. Both go through the same bounded preflight, on one shared budget.
///
/// # Errors
///
/// As [`validate_against`].
pub fn validate_reader_against(
    input: impl BufRead,
    schema: impl BufRead,
    options: &SchemaValidationOptions,
) -> Result<SchemaValidationReport> {
    let o = &options.limits;
    let mut m = Meter::new(o.max_work, o.max_bytes, "XML schema validation");
    let limits = xml_limits(o);
    let schema = self_contained_schema(schema, &limits, &mut m)?;
    run(input, o, &mut m, "XML", |_, _, _| {
        Ok(Grammar::Caller(schema.as_ref()))
    })
}

const XSD_NS: &str = "http://www.w3.org/2001/XMLSchema";

/// Read and check a caller's schema: bounded, well-formed, no DTD, and no
/// composition element that would make libxml2 load another resource.
fn self_contained_schema(input: impl BufRead, limits: &XmlLimits, m: &mut Meter) -> Result<String> {
    let text = document(input, limits, m)?;
    let mut namespaces = Namespaces::default();
    scan_elements(&text, limits, m, |event, _, m| {
        match event {
            Element::Start(tag, attributes) => {
                namespaces.start(tag, attributes, m)?;
                let (prefix, local) = qname(tag, m)?;
                if matches!(local, "include" | "import" | "redefine" | "override")
                    && namespaces.resolve(prefix, false, m)? == Some(XSD_NS)
                {
                    return Err(Error::Unsupported(format!(
                        "caller schema uses xs:{local}; only a self-contained schema is validated"
                    )));
                }
            }
            Element::End(tag) => {
                m.spend(tag.len(), 0)?;
                namespaces.scopes.pop();
            }
        }
        Ok(())
    })?;
    Ok(canonical_declaration(&text, limits.max_input_bytes, m)?.into_owned())
}

/// The grammar one validation compiles.
pub(super) enum Grammar<'a> {
    /// A schema compiled into the crate.
    Bundled(SchemaKind),
    /// Caller schema text, already checked by [`self_contained_schema`].
    Caller(&'a str),
}

fn xml_limits(o: &SchemaValidationLimits) -> XmlLimits {
    XmlLimits {
        max_input_bytes: o.max_xml_bytes.min(i32::MAX as usize),
        max_depth: o.max_depth,
        max_elements: o.max_elements,
    }
}

/// libxml2's schema parser and validation contexts are not safe to use from
/// several threads at once (the binding says so for libxml2 2.12 and later),
/// and registering an input callback races any parse in flight before 2.13.
/// Every use of the engine in this crate holds this lock.
static ENGINE: Mutex<()> = Mutex::new(());

/// The shared pipeline: bounded decoding, one lexical pass with namespace
/// checks, grammar selection at the root, declaration canonicalisation, then
/// the engine under [`ENGINE`].
pub(super) fn run<'g>(
    input: impl BufRead,
    o: &SchemaValidationLimits,
    m: &mut Meter,
    subject: &'static str,
    mut select: impl FnMut(&str, &Attributes, &mut Meter) -> Result<Grammar<'g>>,
) -> Result<SchemaValidationReport> {
    let limits = xml_limits(o);
    let text = document(input, &limits, m)?;
    let mut grammar = None;
    let mut namespaces = Namespaces::default();
    scan_elements(&text, &limits, m, |event, ancestors, m| {
        match event {
            Element::Start(tag, attributes) => {
                namespaces.start(tag, attributes, m)?;
                if ancestors.is_empty() {
                    grammar = Some(select(tag, attributes, m)?);
                }
            }
            Element::End(tag) => {
                m.spend(tag.len(), 0)?;
                namespaces.scopes.pop();
            }
        }
        Ok(())
    })?;
    let grammar =
        grammar.ok_or_else(|| Error::Unsupported(format!("missing {subject} document root")))?;
    let canonical = canonical_declaration(&text, limits.max_input_bytes, m)?;
    let _engine = ENGINE.lock().unwrap_or_else(PoisonError::into_inner);
    // Ordinary Rust panics from the safe binding's null/internal-error branches
    // become typed errors. C faults, allocator aborts and callback unwind aborts
    // cannot be caught.
    catch_unwind(AssertUnwindSafe(|| {
        engine(&canonical, grammar, o, m, subject)
    }))
    .map_err(|_| Error::InvalidValue(format!("{subject} XSD engine panicked during validation")))?
}
// xmlReadMemory may return a DOM despite namespace errors. Check bindings on
// the existing lexical callbacks before C; no second parser or global handler.
const XML_NS: &str = "http://www.w3.org/XML/1998/namespace";
const XMLNS_NS: &str = "http://www.w3.org/2000/xmlns/";
#[derive(Default)]
struct Namespaces {
    scopes: Vec<Vec<(String, String)>>,
}
fn namespace_error() -> Error {
    Error::Parse {
        line: 0,
        message: "invalid XML namespace binding or qualified name".into(),
    }
}
fn qname<'a>(name: &'a str, m: &mut Meter) -> Result<(Option<&'a str>, &'a str)> {
    m.spend(name.len().saturating_mul(2).saturating_add(1), 0)?;
    let (prefix, local) = name
        .split_once(':')
        .map_or((None, name), |(p, l)| (Some(p), l));
    if local.contains(':') || prefix == Some("") || local.is_empty() {
        return Err(namespace_error());
    }
    if let Some(p) = prefix {
        super::cv_xml::name(p)?;
    }
    super::cv_xml::name(local)?;
    Ok((prefix, local))
}
// A narrow lexical guard, not a URI resolver/parser. Keep relative references
// and Unicode namespace names; reject forbidden ASCII and broken %-escapes that
// libxml may silently tolerate even on unused declarations.
fn namespace_uri(uri: &str, m: &mut Meter) -> Result<()> {
    m.spend(uri.len(), 0)?;
    let mut bytes = uri.bytes();
    while let Some(b) = bytes.next() {
        if b == b'%' {
            if !bytes.next().is_some_and(|b| b.is_ascii_hexdigit())
                || !bytes.next().is_some_and(|b| b.is_ascii_hexdigit())
            {
                return Err(namespace_error());
            }
        } else if b.is_ascii()
            && !(b.is_ascii_alphanumeric() || b"-._~:/?#[]@!$&'()*+,;=".contains(&b))
        {
            return Err(namespace_error());
        }
    }
    Ok(())
}
impl Namespaces {
    fn resolve<'a>(
        &'a self,
        prefix: Option<&str>,
        attribute: bool,
        m: &mut Meter,
    ) -> Result<Option<&'a str>> {
        if attribute && prefix.is_none() {
            return Ok(None);
        }
        let prefix = prefix.unwrap_or("");
        m.spend(prefix.len().saturating_add(1), 0)?;
        if prefix == "xml" {
            return Ok(Some(XML_NS));
        }
        if prefix == "xmlns" {
            return Err(namespace_error());
        }
        for scope in self.scopes.iter().rev() {
            m.spend(1, 0)?;
            for (p, uri) in scope.iter().rev() {
                m.spend(p.len().saturating_add(prefix.len()).saturating_add(1), 0)?;
                if p == prefix {
                    return Ok(if uri.is_empty() { None } else { Some(uri) });
                }
            }
        }
        if prefix.is_empty() {
            Ok(None)
        } else {
            Err(namespace_error())
        }
    }
    fn start(&mut self, tag: &str, attrs: &Attributes, m: &mut Meter) -> Result<()> {
        let (prefix, _) = qname(tag, m)?;
        let mut declarations = Vec::new();
        for (key, uri) in attrs {
            let (p, local) = qname(key, m)?;
            let declaration = if key == "xmlns" {
                Some("")
            } else if p == Some("xmlns") {
                Some(local)
            } else {
                None
            };
            if let Some(prefix) = declaration {
                namespace_uri(uri, m)?;
                m.spend(uri.len().saturating_mul(3).saturating_add(prefix.len()), 0)?;
                if prefix == "xmlns"
                    || uri == XMLNS_NS
                    || (prefix == "xml") != (uri == XML_NS)
                    || (!prefix.is_empty() && uri.is_empty())
                {
                    return Err(namespace_error());
                }
                let item = (m.copy(prefix)?, m.copy(uri)?);
                m.push(&mut declarations, item)?;
            }
        }
        m.push(&mut self.scopes, declarations)?;
        self.resolve(prefix, false, m)?;
        for (i, (key, _)) in attrs.iter().enumerate() {
            let (prefix, local) = qname(key, m)?;
            if key == "xmlns" || prefix == Some("xmlns") {
                continue;
            }
            let uri = self.resolve(prefix, true, m)?;
            for (previous, _) in &attrs[..i] {
                let (p, l) = qname(previous, m)?;
                if previous == "xmlns" || p == Some("xmlns") {
                    continue;
                }
                m.spend(local.len().saturating_add(l.len()).saturating_add(1), 0)?;
                if local != l {
                    continue;
                }
                let u = self.resolve(p, true, m)?;
                m.spend(
                    uri.map_or(0, str::len)
                        .saturating_add(u.map_or(0, str::len))
                        .saturating_add(1),
                    0,
                )?;
                if uri == u {
                    return Err(namespace_error());
                }
            }
        }
        Ok(())
    }
}
// document/scan already validated declaration syntax, order and encoding. Change
// only its encoding value, preserving standalone, quote style and line layout.
fn canonical_declaration<'a>(text: &'a str, max: usize, m: &mut Meter) -> Result<Cow<'a, str>> {
    m.spend(text.len(), 0)?;
    if !text.starts_with("<?xml")
        || !text
            .as_bytes()
            .get(5)
            .is_some_and(|b| matches!(b, b' ' | b'\t' | b'\n' | b'\r'))
    {
        return Ok(Cow::Borrowed(text));
    }
    let Some(end) = text.find("?>") else {
        return Err(Error::InvalidValue(
            "validated XML declaration missing terminator".into(),
        ));
    };
    let declaration = &text[..end];
    let Some(key) = declaration.find("encoding") else {
        return Ok(Cow::Borrowed(text));
    };
    let after_key = &declaration[key + "encoding".len()..];
    let Some(quote) = after_key.find(['\'', '"']) else {
        return Err(Error::InvalidValue(
            "validated encoding missing quote".into(),
        ));
    };
    let start = key + "encoding".len() + quote + 1;
    let delimiter = text.as_bytes()[start - 1];
    let Some(length) = declaration[start..].bytes().position(|b| b == delimiter) else {
        return Err(Error::InvalidValue(
            "validated encoding missing close quote".into(),
        ));
    };
    let stop = start + length;
    if &text[start..stop] == "UTF-8" {
        return Ok(Cow::Borrowed(text));
    }
    let capacity = text
        .len()
        .checked_sub(length)
        .and_then(|n| n.checked_add(5))
        .ok_or_else(|| m.limit())?;
    m.cap(capacity, max)?;
    m.spend(capacity, capacity)?;
    let mut out = String::with_capacity(capacity);
    out.push_str(&text[..start]);
    out.push_str("UTF-8");
    out.push_str(&text[stop..]);
    Ok(Cow::Owned(out))
}

/// The private URL scheme under which libxml2 reaches the bundled schemas
/// that `xs:include` others. The included names resolve against the main
/// schema's URL, so they arrive here as `openms-schema:///<file>`.
const SCHEME: &str = "openms-schema:///";

/// Install, once per process, the input handler that serves [`FILES`] under
/// [`SCHEME`].
///
/// libxml2 has no per-context resource loader before 2.14, and this binding
/// exposes none, so an `xs:include` can only be served through the process-wide
/// input-callback table. The handler claims only [`SCHEME`] URLs and answers a
/// name it does not carry with an empty document, so no load under that prefix
/// falls through to libxml2's file, HTTP or FTP loaders. Registration happens
/// under [`ENGINE`], so it never races a parse this crate runs.
fn register_bundle_loader() {
    static REGISTERED: OnceLock<()> = OnceLock::new();
    REGISTERED.get_or_init(|| {
        libxml::io::register_input_callback(
            |url| url.starts_with(SCHEME),
            |url| {
                Some(
                    url.strip_prefix(SCHEME)
                        .and_then(resource)
                        .unwrap_or_default()
                        .to_vec(),
                )
            },
        );
    });
}

fn schema_parser(grammar: &Grammar<'_>) -> Result<(SchemaKind, SchemaParserContext)> {
    match grammar {
        Grammar::Caller(text) => Ok((
            SchemaKind::External,
            SchemaParserContext::from_buffer(text.as_bytes()),
        )),
        Grammar::Bundled(kind) => {
            let bundle = kind.bundle().ok_or_else(|| {
                Error::InvalidValue("SchemaKind::External names no bundled schema".into())
            })?;
            let parser = if bundle.includes {
                register_bundle_loader();
                SchemaParserContext::from_file(&format!("{SCHEME}{}", bundle.file()))
            } else {
                let bytes = resource(bundle.file()).ok_or_else(|| {
                    Error::InvalidValue(format!("{} is not bundled", bundle.file()))
                })?;
                SchemaParserContext::from_buffer(bytes)
            };
            Ok((*kind, parser))
        }
    }
}

fn engine(
    text: &str,
    grammar: Grammar<'_>,
    o: &SchemaValidationLimits,
    m: &mut Meter,
    subject: &'static str,
) -> Result<SchemaValidationReport> {
    let parser = Parser::default();
    let doc = parser
        .parse_string_with_options(
            text.as_bytes(),
            ParserOptions {
                recover: false,
                no_net: true,
                no_def_dtd: false,
                huge: false,
                encoding: None,
                no_error: true,
                no_warning: true,
                ..Default::default()
            },
        )
        .map_err(|_| Error::Parse {
            line: 0,
            message: format!("libxml2 rejected XML during {subject} schema parsing"),
        })?;
    let (schema, mut schema_parser) = schema_parser(&grammar)?;
    let mut validator = match SchemaValidationContext::from_parser(&mut schema_parser) {
        Ok(validator) => validator,
        Err(errors) => {
            let mut diagnostics = Vec::new();
            let mut bytes = 0;
            collect(errors, &mut diagnostics, &mut bytes, o, m)?;
            return Err(Error::InvalidValue(
                diagnostics.into_iter().next().map_or_else(
                    || format!("libxml2 could not compile the {subject} schema"),
                    |d| d.message,
                ),
            ));
        }
    };
    let mut diagnostics = Vec::new();
    let mut bytes = 0;
    collect(
        schema_parser.drain_errors(),
        &mut diagnostics,
        &mut bytes,
        o,
        m,
    )?;
    let (engine_valid, errors) = match validator.validate_document(&doc) {
        Ok(()) => (true, validator.drain_errors()),
        Err(errors) => (false, errors),
    };
    collect(errors, &mut diagnostics, &mut bytes, o, m)?;
    Ok(SchemaValidationReport {
        schema,
        diagnostics,
        engine_valid,
    })
}
// The binding has already allocated errors/messages here. These are post-engine
// result limits; returned strings move without another deep copy.
fn collect(
    errors: Vec<StructuredError>,
    out: &mut Vec<SchemaDiagnostic>,
    bytes: &mut usize,
    o: &SchemaValidationLimits,
    m: &mut Meter,
) -> Result<()> {
    m.cap(out.len().saturating_add(errors.len()), o.max_diagnostics)?;
    m.spend(errors.len(), 0)?;
    for e in errors {
        let message = e.message.unwrap_or_default();
        let n = message
            .len()
            .saturating_add(e.filename.as_ref().map_or(0, String::len));
        *bytes = bytes.checked_add(n).ok_or_else(|| m.limit())?;
        m.cap(*bytes, o.max_diagnostic_bytes)?;
        m.spend(
            n,
            message
                .capacity()
                .saturating_add(e.filename.as_ref().map_or(0, String::capacity)),
        )?;
        let level = match e.level {
            XmlErrorLevel::None => SchemaDiagnosticLevel::Information,
            XmlErrorLevel::Warning => SchemaDiagnosticLevel::Warning,
            XmlErrorLevel::Error => SchemaDiagnosticLevel::Error,
            XmlErrorLevel::Fatal => SchemaDiagnosticLevel::Fatal,
        };
        let positive = |n: i32| usize::try_from(n).ok().filter(|n| *n != 0);
        m.push(
            out,
            SchemaDiagnostic {
                level,
                message,
                filename: e.filename,
                line: e.line.and_then(positive),
                column: e.col.and_then(positive),
                domain: e.domain,
                code: e.code,
            },
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn diagnostic(level: XmlErrorLevel, text: &str) -> StructuredError {
        StructuredError {
            message: Some(text.into()),
            level,
            filename: None,
            line: Some(7),
            col: Some(3),
            domain: 17,
            code: 2,
        }
    }
    #[test]
    fn warning_policy_and_structured_severity_are_not_empty_vector_policy() {
        for (level, valid) in [
            (XmlErrorLevel::None, true),
            (XmlErrorLevel::Warning, false),
            (XmlErrorLevel::Error, false),
            (XmlErrorLevel::Fatal, false),
        ] {
            let mut m = Meter::new(10000, 10000, "test");
            let mut diagnostics = Vec::new();
            let mut bytes = 0;
            collect(
                vec![diagnostic(level, "message")],
                &mut diagnostics,
                &mut bytes,
                &SchemaValidationLimits::default(),
                &mut m,
            )
            .unwrap();
            assert_eq!(diagnostics[0].line, Some(7));
            assert_eq!(diagnostics[0].column, Some(3));
            let mut r = SchemaValidationReport {
                schema: SchemaKind::MzML,
                diagnostics,
                engine_valid: true,
            };
            assert_eq!(r.is_valid(), valid);
            r.engine_valid = false;
            assert!(!r.is_valid());
        }
    }
    #[test]
    fn declaration_replacement_preserves_layout_and_does_not_change_stylesheet_pi() {
        let mut m = Meter::new(10000, 10000, "test");
        let raw = "<?xml version='1.0'\n encoding = 'UTF-16LE' standalone='yes' ?>\n<r/>";
        let wanted = "<?xml version='1.0'\n encoding = 'UTF-8' standalone='yes' ?>\n<r/>";
        assert_eq!(canonical_declaration(raw, 1000, &mut m).unwrap(), wanted);
        let pi = "<?xml-stylesheet encoding='something'?>\n<r/>";
        assert!(matches!(
            canonical_declaration(pi, 1000, &mut m).unwrap(),
            Cow::Borrowed(_)
        ));
        assert_eq!(canonical_declaration(pi, 1000, &mut m).unwrap(), pi);
    }
    #[test]
    fn namespace_resolution_spends_shared_work_after_successful_scope_setup() {
        let mut m = Meter::new(1000, 10000, "test");
        let mut ns = Namespaces::default();
        ns.start(
            "mzML",
            &vec![("xmlns:p".into(), "urn:example".into())],
            &mut m,
        )
        .unwrap();
        let mut successful = 0;
        loop {
            match ns.resolve(Some("p"), true, &mut m) {
                Ok(Some("urn:example")) => successful += 1,
                Err(_) => break,
                other => panic!("unexpected resolution: {other:?}"),
            }
        }
        assert!(successful > 1 && successful < 1000);
        let mut m = Meter::new(1000, 1, "test");
        assert!(
            Namespaces::default()
                .start(
                    "mzML",
                    &vec![("xmlns:p".into(), "urn:example".into())],
                    &mut m
                )
                .is_err()
        );
    }
    #[test]
    fn diagnostic_limits_are_cumulative_after_successful_collection() {
        let o = SchemaValidationLimits {
            max_diagnostics: 2,
            max_diagnostic_bytes: 5,
            ..Default::default()
        };
        let mut m = Meter::new(10000, 10000, "test");
        let mut out = Vec::new();
        let mut bytes = 0;
        collect(
            vec![diagnostic(XmlErrorLevel::Warning, "abc")],
            &mut out,
            &mut bytes,
            &o,
            &mut m,
        )
        .unwrap();
        assert_eq!(bytes, 3);
        assert!(
            collect(
                vec![diagnostic(XmlErrorLevel::Error, "def")],
                &mut out,
                &mut bytes,
                &o,
                &mut m
            )
            .is_err()
        );
        assert_eq!(out.len(), 1);
        let mut bytes = 3;
        assert!(
            collect(
                vec![
                    diagnostic(XmlErrorLevel::Error, ""),
                    diagnostic(XmlErrorLevel::Error, "")
                ],
                &mut out,
                &mut bytes,
                &o,
                &mut m
            )
            .is_err()
        );
        assert_eq!(out.len(), 1);
        let mut m = Meter::new(0, 10000, "test");
        assert!(
            collect(
                vec![diagnostic(XmlErrorLevel::Error, "")],
                &mut out,
                &mut bytes,
                &o,
                &mut m
            )
            .is_err()
        );
    }
    #[test]
    fn every_bundled_kind_names_a_bundled_file_and_every_file_is_reachable() {
        let kinds = [
            SchemaKind::MzML,
            SchemaKind::IndexedMzML,
            SchemaKind::FeatureXML,
            SchemaKind::ConsensusXML,
            SchemaKind::IdXML,
            SchemaKind::ParamXML,
            SchemaKind::TransformationXML,
            SchemaKind::MzData,
            SchemaKind::MzXML,
            SchemaKind::MzIdentML1_0_0,
            SchemaKind::MzIdentML1_1_0,
            SchemaKind::MzIdentML1_2_0,
            SchemaKind::MzIdentML1_3_0,
        ];
        let mut reached = std::collections::BTreeSet::new();
        for kind in kinds {
            let bundle = kind.bundle().unwrap();
            assert_eq!(kind.location(), Some(bundle.location));
            assert!(resource(bundle.file()).is_some(), "{}", bundle.file());
            reached.insert(bundle.file());
            // A schema that includes others is served through the loader;
            // one that does not must never need it.
            let raw = resource(bundle.file()).unwrap();
            let has_include = raw.windows(11).any(|w| w == b":include sc");
            assert_eq!(has_include, bundle.includes, "{}", bundle.file());
        }
        assert!(SchemaKind::External.bundle().is_none());
        // The included files are reached only through the loader.
        for (file, _) in FILES {
            if !reached.contains(file) {
                assert!(
                    [
                        "mzXML_3.1_mod.xsd",
                        "separation_technique_1.0.xsd",
                        "general_types_1.0.xsd",
                        "FuGElightv1.0.0.xsd"
                    ]
                    .contains(&file),
                    "{file}"
                );
            }
        }
    }
}
