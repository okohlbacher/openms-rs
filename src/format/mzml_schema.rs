// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Real XSD validation using fixed historical mzML schemas and optional libxml2:
//! `MzMLFile::isValid`, which picks between the ordinary and the indexed mzML
//! schema. The engine, the report types and the preflight are shared with
//! every other bundled schema in [`xml_schema`](super::xml_schema).
//! Native input/preflight limits and post-engine result limits do not bound the
//! C engine's DOM, identity tables, runtime, or pre-return diagnostic allocations.

use super::cv_xml::{Attributes, Meter};
use super::xml_schema::{Grammar, run};
pub use super::xml_schema::{
    SchemaDiagnostic, SchemaDiagnosticLevel, SchemaKind, SchemaValidationLimits,
    SchemaValidationOptions, SchemaValidationReport,
};
use crate::{Error, Result};
use std::{io::BufRead, path::Path};

const NS: &str = "http://psi.hupo.org/ms/mzml";

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
///
/// The schema is [`SchemaKind::MzML`] or [`SchemaKind::IndexedMzML`], chosen
/// by the root's expanded name; any other root is `Error::Unsupported`.
pub fn validate_schema_reader(
    input: impl BufRead,
    options: &SchemaValidationOptions,
) -> Result<SchemaValidationReport> {
    let o = &options.limits;
    let mut m = Meter::new(o.max_work, o.max_bytes, "mzML schema validation");
    run(input, o, &mut m, "mzML", |tag, attributes, m| {
        Ok(Grammar::Bundled(root_schema(tag, attributes, m)?))
    })
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
