// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Semantic validator for mzData files, the source
//! `FORMAT/VALIDATORS/MzDataValidator.h`.
//!
//! The source specialises `SemanticValidator` by replacing its whole
//! `handleTerm_` (`MzDataValidator.cpp:25-151`). In addition to the generic
//! mapping-rule checks, for each CV term encountered the recorded unit
//! accession, if any, must be one of the term's allowed units or a child of
//! one, and the term's name must match the name the vocabulary registers with
//! whitespace ignored and ASCII case folded. Mismatches are errors. Unit
//! checking is on (`setCheckUnits(true)` in the constructor).
//!
//! The inherited `SemanticValidator::startElement` runs first, unchanged: a
//! term the vocabulary does not know draws the "Unknown CV term" warning and
//! is not handled further, and an obsolete one draws the "Obsolete CV term"
//! warning. The override then differs from the generic handler in two ways:
//! it checks no value type, and the name comparison ignores whitespace and
//! ASCII case, its message quoting both names after that normalisation
//! (`'PSI:1000011 - resolution' should be 'massresolution'`).
//! The source compares a unit's children with the supplied unit accession
//! here, correctly, which is also what this crate's generic handler does
//! (`CPP-039` concerns the generic source handler only).
//!
//! [`validate_semantics`](crate::format::semantic_validator::mzdata::validate_semantics)
//! is source `MzDataFile::isSemanticallyValid` (`MzDataFile.cpp:40-55`): the
//! validator over the bundled `resources/cv/mzdata-mapping.xml` and
//! `resources/cv/psi-mzdata.obo`, unchanged copies of the pinned
//! `share/OpenMS/MAPPING/mzdata-mapping.xml` and `share/OpenMS/CV/psi-mzdata.obo`.
//! See `docs/MZDATA_SUPPORT.md`.

use super::*;
use std::io::Read;
use std::sync::OnceLock;

pub use super::{ParsedCVTerm, ValidationLimits, ValidationOptions, ValidationReport};

/// The source `MzDataValidator`: borrowed mapping and vocabulary, fresh rule
/// state for every call.
///
/// `options` holds the inherited `SemanticValidator` settings; as in the
/// source constructor, unit checking starts on. `check_term_value_types` has
/// no effect here, because the source's `handleTerm_` override checks no
/// value type.
#[derive(Debug)]
pub struct MzDataValidator<'a> {
    mapping: &'a CVMappings,
    cv: &'a ControlledVocabulary,
    /// Term element, attribute names, unit checking and ceilings.
    pub options: ValidationOptions,
}

impl<'a> MzDataValidator<'a> {
    /// Source `MzDataValidator(mapping, cv)`: `cv` must hold *all* the
    /// vocabularies the mapping needs.
    pub fn new(mapping: &'a CVMappings, cv: &'a ControlledVocabulary) -> Self {
        Self {
            mapping,
            cv,
            options: ValidationOptions {
                check_units: true,
                ..Default::default()
            },
        }
    }

    /// Validate a file, source `validate(filename, errors, warnings)`: the
    /// report carries the two output lists, and the source's `bool` is
    /// [`ValidationReport::is_valid`].
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be opened, where the source throws
    /// `Exception::FileNotFound`; otherwise as
    /// [`MzDataValidator::validate_reader`].
    pub fn validate(&self, path: impl AsRef<Path>) -> Result<ValidationReport> {
        self.validate_reader(crate::format::path_io::open(path.as_ref())?)
    }

    /// Validate a document read from `input`.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] for input that is not well-formed XML or a term
    /// element without its accession or name attribute, and
    /// [`Error::InvalidValue`] when a [`ValidationLimits`] ceiling is exceeded.
    /// Semantic findings are never errors; they are the report.
    pub fn validate_reader(&self, input: impl BufRead) -> Result<ValidationReport> {
        let o = &self.options;
        let mut m = o.limits.meter();
        let mut index = Index::new(self.mapping, &o.limits, &mut m)?;
        let bytes = latin1_as_utf8(input, &o.limits, &mut m)?;
        let text = document(bytes.as_slice(), &o.limits.xml(), &mut m)?;
        let mut report = ValidationReport::default();
        let mut terms = 0usize;
        scan_elements(&text, &o.limits.xml(), &mut m, |event, ancestors, m| {
            match event {
                Element::Start(tag, a) => {
                    m.spend(tag.len().min(o.tag.len()).saturating_add(1), 0)?;
                    if tag == o.tag {
                        terms = terms.checked_add(1).ok_or_else(|| m.limit())?;
                        m.cap(terms, o.limits.max_terms)?;
                        let p = parse_term(a, o, m)?;
                        let parent = path(ancestors, None, m)?;
                        // The inherited startElement: unknown terms stop here.
                        let Some(term) = encounter(self.cv, &p, &parent, o, &mut report, m)? else {
                            return Ok(());
                        };
                        let key = cv_path(&parent, o, m)?;
                        handle_mzdata(
                            &mut index,
                            self.cv,
                            term,
                            &p,
                            &key,
                            &parent,
                            o,
                            &mut report,
                            m,
                        )?;
                    }
                }
                Element::End(tag) => {
                    let parent = path(ancestors, Some(tag), m)?;
                    let key = cv_path(&parent, o, m)?;
                    index.close(&key, &parent, o, &mut report, m)?;
                }
            }
            Ok(())
        })?;
        Ok(report)
    }

    /// Selection only, the inherited source `locateTerm(path, parsed_term)`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a path no mapping rule names, where the
    /// source's answer depends on earlier validations (`CPP-044`).
    pub fn locate_term(&self, path: &str, term: &ParsedCVTerm) -> Result<bool> {
        locate(self.mapping, self.cv, &self.options, path, term)
    }
}

/// The document as UTF-8, transcoding a document that declares ISO-8859-1.
///
/// mzData's writer declares `encoding="ISO-8859-1"`, the mzData reader of this
/// crate transcodes such a document exactly (ISO-8859-1 maps byte `b` to
/// `U+00b`), and the source's Xerces parser decodes it. The shared validator
/// input refuses a non-ASCII document with that declaration, so the mzData
/// validator transcodes first and rewrites the declared encoding to `UTF-8`.
/// A document that declares no encoding, or another one, is passed on
/// unchanged, and the shared input decides.
fn latin1_as_utf8(
    mut input: impl BufRead,
    limits: &ValidationLimits,
    m: &mut Meter,
) -> Result<Vec<u8>> {
    let cap = u64::try_from(limits.max_input_bytes)
        .unwrap_or(u64::MAX)
        .saturating_add(1);
    let mut bytes = Vec::new();
    (&mut input).take(cap).read_to_end(&mut bytes)?;
    if bytes.len() > limits.max_input_bytes {
        return Err(m.limit());
    }
    m.spend(bytes.len(), bytes.len())?;
    if bytes.is_ascii() || !bytes.starts_with(b"<?xml") {
        return Ok(bytes);
    }
    let Some(end) = bytes.windows(2).position(|w| w == b"?>") else {
        return Ok(bytes);
    };
    let declaration = &bytes[..end];
    let Some(at) = declaration
        .windows(8)
        .position(|w| w.eq_ignore_ascii_case(b"encoding"))
    else {
        return Ok(bytes);
    };
    let rest = &declaration[at + 8..];
    let Some(eq) = rest.iter().position(|&b| b == b'=') else {
        return Ok(bytes);
    };
    let after = &rest[eq + 1..];
    let Some(open) = after.iter().position(|&b| b == b'"' || b == b'\'') else {
        return Ok(bytes);
    };
    let quote = after[open];
    let value = &after[open + 1..];
    let Some(close) = value.iter().position(|&b| b == quote) else {
        return Ok(bytes);
    };
    let encoding = &value[..close];
    if !(encoding.eq_ignore_ascii_case(b"ISO-8859-1")
        || encoding.eq_ignore_ascii_case(b"ISO8859-1")
        || encoding.eq_ignore_ascii_case(b"latin1"))
    {
        return Ok(bytes);
    }
    let start = at + 8 + eq + 1 + open + 1;
    m.spend(bytes.len().saturating_mul(2), bytes.len().saturating_mul(2))?;
    let mut text = String::new();
    text.try_reserve(bytes.len().saturating_mul(2))
        .map_err(|_| m.limit())?;
    text.push_str(std::str::from_utf8(&bytes[..start]).map_err(|_| m.limit())?);
    text.push_str("UTF-8");
    for &byte in &bytes[start + close..] {
        text.push(char::from(byte));
    }
    Ok(text.into_bytes())
}

/// Source `StringUtils::trim`, then `toLower` (ASCII, the "C" locale) and
/// `removeWhitespaces` (space, tab, LF, CR), in that order.
fn normalised_name(name: &str) -> String {
    trim(name)
        .chars()
        .filter(|c| !matches!(c, ' ' | '\t' | '\n' | '\r'))
        .map(|c| c.to_ascii_lowercase())
        .collect()
}

/// Source `MzDataValidator::handleTerm_` (`MzDataValidator.cpp:25-151`).
#[allow(clippy::too_many_arguments)] // One private source callback carries its borrowed operation context.
fn handle_mzdata(
    index: &mut Index<'_>,
    cv: &ControlledVocabulary,
    term: &CVTermDefinition,
    p: &ParsedCVTerm,
    key: &str,
    parent: &str,
    o: &ValidationOptions,
    report: &mut ValidationReport,
    m: &mut Meter,
) -> Result<()> {
    let n = p
        .accession
        .len()
        .saturating_add(p.name.len())
        .saturating_add(parent.len());
    let mut allowed = false;
    let mut rule_found = false;
    if let Some(group) = index.get(key, m)? {
        m.spend(group.rules.len(), 0)?;
        for i in 0..group.rules.len() {
            let r = group.rules[i];
            rule_found = true;
            for t in &r.terms {
                if matches_term(
                    cv,
                    &t.accession,
                    t.use_term,
                    t.allow_children,
                    &p.accession,
                    m,
                )? {
                    allowed = true;
                    group.increment(&r.identifier, &t.accession, m)?;
                    break;
                }
            }
        }
    }
    if o.check_units {
        units(cv, term, p, o, report, m)?;
    }
    if !rule_found {
        emit(
            report,
            false,
            o,
            m,
            parent.len(),
            format_args!("No mapping rule found for element '{parent}'"),
        )?;
    } else if !allowed {
        emit(
            report,
            true,
            o,
            m,
            n,
            format_args!(
                "CV term used in invalid element: '{} - {}' at element '{parent}'",
                p.accession, p.name
            ),
        )?;
    }
    m.spend(
        p.name
            .len()
            .saturating_add(term.name.len())
            .saturating_mul(3),
        p.name.len().saturating_add(term.name.len()),
    )?;
    let (parsed, correct) = (normalised_name(&p.name), normalised_name(&term.name));
    if parsed != correct {
        emit(
            report,
            true,
            o,
            m,
            p.accession
                .len()
                .saturating_add(parsed.len())
                .saturating_add(correct.len()),
            format_args!(
                "Name of CV term not correct: '{} - {parsed}' should be '{correct}'",
                p.accession
            ),
        )?;
    }
    Ok(())
}

fn default_mapping() -> Result<&'static CVMappings> {
    static MAPPING: OnceLock<std::result::Result<CVMappings, String>> = OnceLock::new();
    match MAPPING.get_or_init(|| {
        crate::format::cv_mapping::CVMappingFile::default()
            .read(include_bytes!("../../resources/cv/mzdata-mapping.xml").as_slice())
            .map_err(|e| e.to_string())
    }) {
        Ok(mapping) => Ok(mapping),
        Err(message) => Err(Error::InvalidValue(message.clone())),
    }
}

/// The pinned `psi-mzdata.obo`, loaded once under the name `PSI`, as source
/// `cv.loadFromOBO("PSI", File::find("/CV/psi-mzdata.obo"))`.
///
/// # Errors
///
/// [`Error::InvalidValue`] if the bundled vocabulary fails to load; the error
/// is cached and returned on every later call.
pub fn default_vocabulary() -> Result<&'static ControlledVocabulary> {
    static CV: OnceLock<std::result::Result<ControlledVocabulary, String>> = OnceLock::new();
    match CV.get_or_init(|| {
        let mut cv = ControlledVocabulary::new();
        cv.load_obo_reader(
            "PSI",
            include_bytes!("../../resources/cv/psi-mzdata.obo").as_slice(),
        )
        .map_err(|e| e.to_string())?;
        Ok(cv)
    }) {
        Ok(cv) => Ok(cv),
        Err(message) => Err(Error::InvalidValue(message.clone())),
    }
}

/// The pinned mzData mapping rules, `share/OpenMS/MAPPING/mzdata-mapping.xml`.
///
/// # Errors
///
/// [`Error::InvalidValue`] if the bundled mapping fails to load; the error is
/// cached and returned on every later call.
pub fn default_mappings() -> Result<&'static CVMappings> {
    default_mapping()
}

/// Source `MzDataFile::isSemanticallyValid(filename, errors, warnings)`
/// (`MzDataFile.cpp:40-55`): validate `path` against the bundled mzData
/// mapping and vocabulary. The report's `errors` and `warnings` are the
/// source's two output lists; its [`ValidationReport::is_valid`] is the
/// returned `bool`.
///
/// # Errors
///
/// As [`MzDataValidator::validate`], plus [`Error::InvalidValue`] if a bundled
/// resource fails to load.
pub fn validate_semantics(path: impl AsRef<Path>) -> Result<ValidationReport> {
    MzDataValidator::new(default_mapping()?, default_vocabulary()?).validate(path)
}
