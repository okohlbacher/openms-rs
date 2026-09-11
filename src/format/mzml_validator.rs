// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Complete mzML-specific semantic callbacks over the shared CV/XML rule engine.
//! This checks CV semantics, not XSD structure or encoded scientific arrays.

use super::*;
use std::sync::OnceLock;

pub use super::{ParsedCVTerm, ValidationLimits, ValidationOptions, ValidationReport};

/// Immutable borrowed mappings/CVs with fresh group, binary and rule state per call.
/// All inherited source settings are represented by `options`; units default on.
#[derive(Debug)]
pub struct MzMLValidator<'a> {
    mapping: &'a CVMappings,
    cv: &'a ControlledVocabulary,
    pub options: ValidationOptions,
}
impl<'a> MzMLValidator<'a> {
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
    pub fn validate(&self, path: impl AsRef<Path>) -> Result<ValidationReport> {
        self.validate_reader(crate::format::path_io::open(path.as_ref())?)
    }
    pub fn validate_reader(&self, input: impl BufRead) -> Result<ValidationReport> {
        self.validate_reader_options(input, &self.options)
    }
    fn validate_reader_options(
        &self,
        input: impl BufRead,
        o: &ValidationOptions,
    ) -> Result<ValidationReport> {
        let mut m = o.limits.meter();
        let mut index = Index::new(self.mapping, &o.limits, &mut m)?;
        let text = document(input, &o.limits.xml(), &mut m)?;
        let mut report = ValidationReport::default();
        let mut groups = ParameterGroups::default();
        let mut binary = BinaryTypes::default();
        let mut terms = 0usize;
        scan_elements(&text, &o.limits.xml(), &mut m, |event, ancestors, m| {
            let root = ancestors.is_empty();
            let ancestors = mzml_ancestors(ancestors, m)?;
            match event {
                Element::Start(tag, attributes) => {
                    // Preserve the source's fixed tag precedence, even when the
                    // configurable CV tag has the same spelling as one of these.
                    m.spend(tag.len().saturating_add(o.tag.len()).saturating_add(128), 0)?;
                    if tag == "referenceableParamGroup" {
                        groups.current = required(attributes, "id", m)?;
                    } else if tag == "referenceableParamGroupRef" {
                        let id = required(attributes, "ref", m)?;
                        let group = groups.entry(&id, m)?;
                        count_terms(&mut terms, group.len(), o, m)?;
                        m.spend(group.len(), 0)?;
                        let parent = path(ancestors, None, m)?;
                        let key = cv_path(&parent, o, m)?;
                        for p in group {
                            // References borrow stored terms. Charge repeated
                            // inspection and all ensuing queries/diagnostics in
                            // the same operation; no per-reference allowance.
                            m.spend(parameter_bytes(p), 0)?;
                            let term = self.cv.get_term_with_budget(
                                &p.accession,
                                &mut m.work,
                                &mut m.bytes,
                            )?;
                            binary.handle(
                                &mut index,
                                self.cv,
                                term,
                                p,
                                &key,
                                &parent,
                                o,
                                &mut report,
                                m,
                            )?;
                        }
                    } else if tag == "binaryDataArray" {
                        binary = BinaryTypes::default();
                    } else if tag == o.tag {
                        count_terms(&mut terms, 1, o, m)?;
                        let p = parse_term(attributes, o, m)?;
                        let parent = path(ancestors, None, m)?;
                        let Some(term) = encounter(self.cv, &p, &parent, o, &mut report, m)? else {
                            return Ok(());
                        };
                        if ancestors
                            .last()
                            .is_some_and(|s| s == "referenceableParamGroup")
                        {
                            groups.push(p, m)?;
                        } else {
                            let key = cv_path(&parent, o, m)?;
                            binary.handle(
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
                }
                Element::End(tag) => {
                    // At the indexed root itself the source strips that one tag
                    // too; every other end callback includes its current tag.
                    let last = if root && tag == "indexedmzML" {
                        None
                    } else {
                        Some(tag)
                    };
                    let parent = path(ancestors, last, m)?;
                    let key = cv_path(&parent, o, m)?;
                    index.close(&key, &parent, o, &mut report, m)?;
                }
            }
            Ok(())
        })?;
        Ok(report)
    }
    /// Same mapping-only selection as the inherited source operation. It does
    /// not expand groups or suppress GO/BTO terms. Missing paths are stable
    /// checked errors across calls (the existing CPP-044 correction).
    pub fn locate_term(&self, path: &str, term: &ParsedCVTerm) -> Result<bool> {
        locate(self.mapping, self.cv, &self.options, path, term)
    }
}
fn mzml_ancestors<'a>(ancestors: &'a [String], m: &mut Meter) -> Result<&'a [String]> {
    m.spend(12, 0)?;
    Ok(if ancestors.first().is_some_and(|s| s == "indexedmzML") {
        &ancestors[1..]
    } else {
        ancestors
    })
}
fn required(a: &Attributes, key: &str, m: &mut Meter) -> Result<String> {
    let value = attribute(a, key, m)?
        .ok_or_else(|| bad(format!("missing mzML semantic attribute {key}")))?;
    m.copy(value)
}
fn count_terms(count: &mut usize, n: usize, o: &ValidationOptions, m: &mut Meter) -> Result<()> {
    *count = count.checked_add(n).ok_or_else(|| m.limit())?;
    m.cap(*count, o.limits.max_terms)
}
fn parameter_bytes(p: &ParsedCVTerm) -> usize {
    p.accession
        .len()
        .saturating_add(p.name.len())
        .saturating_add(p.value.len())
        .saturating_add(p.unit_accession.len())
        .saturating_add(p.unit_name.len())
        .saturating_add(1)
}
#[derive(Default)]
struct ParameterGroups {
    groups: BTreeMap<String, Vec<ParsedCVTerm>>,
    current: String,
    longest_id: usize,
}
impl ParameterGroups {
    fn entry(&mut self, id: &str, m: &mut Meter) -> Result<&mut Vec<ParsedCVTerm>> {
        m.spend(
            comparison(self.groups.len(), self.longest_id.min(id.len())),
            0,
        )?;
        if !self.groups.contains_key(id) {
            m.tree::<(String, Vec<ParsedCVTerm>)>(1)?;
            let key = m.copy(id)?;
            self.longest_id = self.longest_id.max(id.len());
            m.spend(
                comparison(self.groups.len(), self.longest_id.min(id.len())),
                0,
            )?;
            self.groups.insert(key, Vec::new());
        }
        m.spend(
            comparison(self.groups.len(), self.longest_id.min(id.len())),
            0,
        )?;
        Ok(self.groups.get_mut(id).expect("inserted group"))
    }
    fn push(&mut self, term: ParsedCVTerm, m: &mut Meter) -> Result<()> {
        // A nested definition changes the source current ID permanently; closing
        // it does not restore an outer ID. No schema nesting is fabricated here.
        let id = m.copy(&self.current)?;
        let group = self.entry(&id, m)?;
        m.push(group, term)
    }
}
#[derive(Default)]
struct BinaryTypes<'a> {
    array: Option<&'a CVTermDefinition>,
    value_type: Option<&'a CVTermDefinition>,
}
impl<'a> BinaryTypes<'a> {
    #[allow(clippy::too_many_arguments)] // Exact concrete source callback over the shared rule engine.
    fn handle(
        &mut self,
        index: &mut Index<'_>,
        cv: &'a ControlledVocabulary,
        term: &'a CVTermDefinition,
        p: &ParsedCVTerm,
        key: &str,
        parent: &str,
        o: &ValidationOptions,
        report: &mut ValidationReport,
        m: &mut Meter,
    ) -> Result<()> {
        m.spend(8, 0)?;
        if p.accession.starts_with("GO:") || p.accession.starts_with("BTO:") {
            return Ok(());
        }
        const SUFFIX: &str = "/binaryDataArray/cvParam/@accession";
        m.spend(SUFFIX.len(), 0)?;
        if key.ends_with(SUFFIX) {
            if cv.is_child_of_with_budget(&p.accession, "MS:1000513", &mut m.work, &mut m.bytes)? {
                self.array = Some(term);
            }
            if cv.is_child_of_with_budget(&p.accession, "MS:1000518", &mut m.work, &mut m.bytes)? {
                self.value_type = Some(term);
            }
            if let (Some(array), Some(value_type)) = (self.array, self.value_type) {
                let mut allowed = false;
                m.spend(array.xref_binary.len(), 0)?;
                for id in &array.xref_binary {
                    m.spend(id.len().min(value_type.id.len()).saturating_add(1), 0)?;
                    if id == &value_type.id {
                        allowed = true;
                        break;
                    }
                }
                if !allowed {
                    emit(
                        report,
                        true,
                        o,
                        m,
                        array
                            .id
                            .len()
                            .saturating_add(array.name.len())
                            .saturating_add(value_type.id.len())
                            .saturating_add(value_type.name.len()),
                        format_args!(
                            "Binary data array of type '{} ! {}' cannot have the value type '{} ! {}'.",
                            array.id, array.name, value_type.id, value_type.name
                        ),
                    )?;
                }
            }
        }
        handle_known(index, cv, term, p, key, parent, o, report, m)
    }
}
fn default_mapping() -> Result<&'static CVMappings> {
    static MAPPING: OnceLock<std::result::Result<CVMappings, String>> = OnceLock::new();
    match MAPPING.get_or_init(|| {
        crate::format::cv_mapping::CVMappingFile::default()
            .read(include_bytes!("../../resources/cv/ms-mapping.xml").as_slice())
            .map_err(|e| e.to_string())
    }) {
        Ok(mapping) => Ok(mapping),
        Err(message) => Err(Error::InvalidValue(message.clone())),
    }
}
/// Validate with the exact pinned mzML mapping and five-provider vocabulary.
/// Provider construction uses its existing separate fixed bounds, once per process.
pub fn validate_semantics(path: impl AsRef<Path>) -> Result<ValidationReport> {
    let validator = MzMLValidator::new(default_mapping()?, ControlledVocabulary::psi_ms()?);
    validator.validate(path)
}
/// As above with explicit term/attribute/check/resource settings. This operation
/// never uses mutable global parser state or changes the supplied options.
pub fn validate_semantics_with_options(
    path: impl AsRef<Path>,
    options: &ValidationOptions,
) -> Result<ValidationReport> {
    let validator = MzMLValidator::new(default_mapping()?, ControlledVocabulary::psi_ms()?);
    validator.validate_reader_options(crate::format::path_io::open(path.as_ref())?, options)
}
