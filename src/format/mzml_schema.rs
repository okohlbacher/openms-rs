// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Real XSD validation using fixed historical mzML schemas and optional libxml2.
//! Native input/preflight limits and post-engine result limits do not bound the
//! C engine's DOM, identity tables, runtime, or pre-return diagnostic allocations.

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
};

const NS: &str = "http://psi.hupo.org/ms/mzml";
const ORDINARY: &[u8] = include_bytes!("../../resources/schemas/mzML_1_10.xsd");
const INDEXED: &[u8] = include_bytes!("../../resources/schemas/mzML_idx_1_10.xsd");

/// Selected by root namespace and local name, independently of prefix/prolog.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SchemaKind {
    MzML,
    IndexedMzML,
}
/// Structured libxml2 severity. Informational messages do not invalidate a file.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SchemaDiagnosticLevel {
    Information,
    Warning,
    Error,
    Fatal,
}
/// Owned engine diagnostic. Text/codes depend on the installed libxml2 version.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SchemaDiagnostic {
    pub level: SchemaDiagnosticLevel,
    pub message: String,
    pub filename: Option<String>,
    pub line: Option<usize>,
    pub column: Option<usize>,
    pub domain: i32,
    pub code: i32,
}
/// Schema validity does not establish CV semantics, checksums or index integrity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SchemaValidationReport {
    pub schema: SchemaKind,
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
    /// Maximum encoded input bytes and decoded/canonical UTF-8 bytes.
    pub max_xml_bytes: usize,
    pub max_depth: usize,
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
/// Fixed-schema validation settings. No external schema or resolver is accepted.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SchemaValidationOptions {
    pub limits: SchemaValidationLimits,
}

/// Validate a plain/gzip/bzip2 file using the default native limits.
pub fn validate_schema(path: impl AsRef<Path>) -> Result<SchemaValidationReport> {
    validate_schema_with_options(path, &SchemaValidationOptions::default())
}
/// Validate a file using explicit native limits. Compression is detected by bytes.
pub fn validate_schema_with_options(
    path: impl AsRef<Path>,
    options: &SchemaValidationOptions,
) -> Result<SchemaValidationReport> {
    validate_schema_reader(super::path_io::open(path.as_ref())?, options)
}
/// Validate XML after full bounded decoding/lexical preflight. Semantic XSD
/// violations are Ok(report); malformed/unsupported input, I/O and limits are Err.
pub fn validate_schema_reader(
    input: impl BufRead,
    options: &SchemaValidationOptions,
) -> Result<SchemaValidationReport> {
    let o = &options.limits;
    let mut m = Meter::new(o.max_work, o.max_bytes, "mzML schema validation");
    let limits = XmlLimits {
        max_input_bytes: o.max_xml_bytes.min(i32::MAX as usize),
        max_depth: o.max_depth,
        max_elements: o.max_elements,
    };
    let text = document(input, &limits, &mut m)?;
    let mut schema = None;
    let mut namespaces = Namespaces::default();
    scan_elements(&text, &limits, &mut m, |event, ancestors, m| {
        match event {
            Element::Start(tag, attributes) => {
                namespaces.start(tag, attributes, m)?;
                if ancestors.is_empty() {
                    schema = Some(root_schema(tag, attributes, m)?);
                }
            }
            Element::End(tag) => {
                m.spend(tag.len(), 0)?;
                namespaces.scopes.pop();
            }
        }
        Ok(())
    })?;
    let schema = schema.ok_or_else(|| Error::Unsupported("missing mzML document root".into()))?;
    let canonical = canonical_declaration(&text, limits.max_input_bytes, &mut m)?;
    // Ordinary Rust panics from the safe binding's null/internal-error branches
    // become typed errors. C faults, allocator aborts and callback unwind aborts
    // cannot be caught; no process-global panic/error/IO handler is installed.
    catch_unwind(AssertUnwindSafe(|| engine(&canonical, schema, o, &mut m)))
        .map_err(|_| Error::InvalidValue("mzML XSD engine panicked during validation".into()))?
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
        message: "invalid mzML XML namespace binding or qualified name".into(),
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
fn root_schema(tag: &str, attributes: &Attributes, m: &mut Meter) -> Result<SchemaKind> {
    m.spend(tag.len().saturating_add(1), 0)?;
    let (prefix, local) = tag
        .split_once(':')
        .map_or((None, tag), |(p, l)| (Some(p), l));
    let mut namespace = None;
    for (key, value) in attributes {
        m.spend(key.len().saturating_add(tag.len()).saturating_add(1), 0)?;
        if match prefix {
            None => key == "xmlns",
            Some(p) => key.strip_prefix("xmlns:") == Some(p),
        } {
            namespace = Some(value.as_str());
        }
    }
    if let Some(value) = namespace {
        m.spend(value.len(), 0)?;
    }
    if namespace != Some(NS) {
        return Err(Error::Unsupported("mzML schema root namespace".into()));
    }
    match local {
        "mzML" => Ok(SchemaKind::MzML),
        "indexedmzML" => Ok(SchemaKind::IndexedMzML),
        _ => Err(Error::Unsupported("mzML schema root element".into())),
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
fn engine(
    text: &str,
    schema: SchemaKind,
    o: &SchemaValidationLimits,
    m: &mut Meter,
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
            message: "libxml2 rejected XML during mzML schema parsing".into(),
        })?;
    // Static buffers remain alive; both original schemas have no external imports.
    let raw = match schema {
        SchemaKind::MzML => ORDINARY,
        SchemaKind::IndexedMzML => INDEXED,
    };
    let mut schema_parser = SchemaParserContext::from_buffer(raw);
    let mut validator = match SchemaValidationContext::from_parser(&mut schema_parser) {
        Ok(validator) => validator,
        Err(errors) => {
            let mut diagnostics = Vec::new();
            let mut bytes = 0;
            collect(errors, &mut diagnostics, &mut bytes, o, m)?;
            return Err(Error::InvalidValue(
                diagnostics.into_iter().next().map_or_else(
                    || "libxml2 could not compile the pinned mzML schema".into(),
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
}
