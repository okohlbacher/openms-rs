// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! General source CV mapping validation with operation-local XML and rule state.

#[cfg(feature = "mzml-validation")]
#[path = "mzml_validator.rs"]
pub mod mzml;

use super::controlled_vocabulary::{CVTermDefinition, ControlledVocabulary, XRefType};
use super::cv_xml::{Attributes, Element, Meter, XmlLimits, document, scan_elements};
use crate::data_structures::list::ListParse;
use crate::data_structures::{
    CVMappingRule, CVMappings, CombinationsLogic as Logic, DateTime, RequirementLevel as Level,
};
use crate::{Error, Result};
use std::{collections::BTreeMap, fmt, io::BufRead, path::Path};

/// Parsed source CV parameter; presence flags distinguish absent and empty attributes.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ParsedCVTerm {
    pub accession: String,
    pub name: String,
    pub value: String,
    pub has_value: bool,
    pub unit_accession: String,
    pub has_unit_accession: bool,
    pub unit_name: String,
    pub has_unit_name: bool,
}

#[derive(Clone, Copy, Debug)]
pub struct ValidationLimits {
    pub max_input_bytes: usize,
    pub max_depth: usize,
    pub max_elements: usize,
    pub max_terms: usize,
    pub max_rules: usize,
    pub max_mapping_terms: usize,
    pub max_diagnostics: usize,
    pub max_work: usize,
    pub max_bytes: usize,
}
impl Default for ValidationLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 16 * 1024 * 1024,
            max_depth: 128,
            max_elements: 1_000_000,
            max_terms: 1_000_000,
            max_rules: 100_000,
            max_mapping_terms: 1_000_000,
            max_diagnostics: 100_000,
            max_work: 50_000_000,
            max_bytes: 128 * 1024 * 1024,
        }
    }
}
impl ValidationLimits {
    fn xml(self) -> XmlLimits {
        XmlLimits {
            max_input_bytes: self.max_input_bytes,
            max_depth: self.max_depth,
            max_elements: self.max_elements,
        }
    }
    fn meter(self) -> Meter {
        Meter::new(self.max_work, self.max_bytes, "Semantic validation")
    }
}

/// Public fields replace source setter methods; strings are exact XML qnames.
#[derive(Clone, Debug)]
pub struct ValidationOptions {
    pub tag: String,
    pub accession_attribute: String,
    pub name_attribute: String,
    pub value_attribute: String,
    pub unit_accession_attribute: String,
    pub unit_name_attribute: String,
    pub check_term_value_types: bool,
    pub check_units: bool,
    pub limits: ValidationLimits,
}
impl Default for ValidationOptions {
    fn default() -> Self {
        Self {
            tag: "cvParam".into(),
            accession_attribute: "accession".into(),
            name_attribute: "name".into(),
            value_attribute: "value".into(),
            unit_accession_attribute: "unitAccession".into(),
            unit_name_attribute: "unitName".into(),
            check_term_value_types: true,
            check_units: false,
            limits: ValidationLimits::default(),
        }
    }
}

/// Semantic failures are an ordinary report. XML, I/O and resource failures are Err.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ValidationReport {
    pub errors: Vec<String>,
    pub warnings: Vec<String>,
}
impl ValidationReport {
    pub fn is_valid(&self) -> bool {
        self.errors.is_empty()
    }
}

/// Borrowed immutable providers; each call builds bounded private evaluation state.
#[derive(Debug)]
pub struct SemanticValidator<'a> {
    mapping: &'a CVMappings,
    cv: &'a ControlledVocabulary,
    pub options: ValidationOptions,
}
impl<'a> SemanticValidator<'a> {
    pub fn new(mapping: &'a CVMappings, cv: &'a ControlledVocabulary) -> Self {
        Self {
            mapping,
            cv,
            options: ValidationOptions::default(),
        }
    }
    pub fn validate(&self, path: impl AsRef<Path>) -> Result<ValidationReport> {
        self.validate_reader(super::path_io::open(path.as_ref())?)
    }
    pub fn validate_reader(&self, input: impl BufRead) -> Result<ValidationReport> {
        let o = &self.options;
        let mut m = o.limits.meter();
        let mut index = Index::new(self.mapping, &o.limits, &mut m)?;
        let text = document(input, &o.limits.xml(), &mut m)?;
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
                        let key = cv_path(&parent, o, m)?;
                        handle(&mut index, self.cv, &p, &key, &parent, o, &mut report, m)?;
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
    /// Selection only. An absent mapping path is always an error, independent of
    /// earlier validate calls (the source persistent empty-cache quirk is CPP-044).
    pub fn locate_term(&self, path: &str, term: &ParsedCVTerm) -> Result<bool> {
        locate(self.mapping, self.cv, &self.options, path, term)
    }
}
fn locate(
    mapping: &CVMappings,
    cv: &ControlledVocabulary,
    options: &ValidationOptions,
    path: &str,
    term: &ParsedCVTerm,
) -> Result<bool> {
    let mut m = options.limits.meter();
    let mut index = Index::new(mapping, &options.limits, &mut m)?;
    let group = index
        .get(path, &mut m)?
        .ok_or_else(|| Error::InvalidValue("unknown semantic mapping path".into()))?;
    m.spend(group.rules.len().saturating_mul(2), 0)?;
    for r in &group.rules {
        m.spend(r.terms.len(), 0)?;
        for t in &r.terms {
            if matches_term(
                cv,
                &t.accession,
                t.use_term,
                t.allow_children,
                &term.accession,
                &mut m,
            )? {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn bad(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}
fn trim(s: &str) -> &str {
    s.trim_matches([' ', '\t', '\r', '\n'])
}
fn path(ancestors: &[String], last: Option<&str>, m: &mut Meter) -> Result<String> {
    m.spend(ancestors.len().saturating_add(1), 0)?;
    let mut n = 1usize;
    for s in ancestors.iter().map(String::as_str).chain(last) {
        n = n
            .checked_add(s.len().saturating_add(1))
            .ok_or_else(|| m.limit())?;
    }
    m.spend(n, n)?;
    let mut result = String::with_capacity(n);
    result.push('/');
    for (i, s) in ancestors.iter().map(String::as_str).chain(last).enumerate() {
        if i != 0 {
            result.push('/');
        }
        result.push_str(s);
    }
    Ok(result)
}
fn cv_path(parent: &str, o: &ValidationOptions, m: &mut Meter) -> Result<String> {
    let n = parent
        .len()
        .saturating_add(o.tag.len())
        .saturating_add(o.accession_attribute.len())
        .saturating_add(3);
    m.spend(n, n)?;
    Ok(format!("{parent}/{}/@{}", o.tag, o.accession_attribute))
}
fn attribute<'a>(a: &'a Attributes, key: &str, m: &mut Meter) -> Result<Option<&'a str>> {
    for (k, v) in a {
        m.spend(k.len().min(key.len()).saturating_add(1), 0)?;
        if k == key {
            return Ok(Some(v));
        }
    }
    Ok(None)
}
fn parse_term(a: &Attributes, o: &ValidationOptions, m: &mut Meter) -> Result<ParsedCVTerm> {
    let mut required = |key: &str| -> Result<String> {
        let value = attribute(a, key, m)?;
        if value.is_none() {
            m.spend(key.len().saturating_add(64), key.len().saturating_add(64))?;
        }
        let v = value.ok_or_else(|| bad(format!("missing CV parameter attribute {key}")))?;
        m.copy(v)
    };
    let accession = required(&o.accession_attribute)?;
    let name = required(&o.name_attribute)?;
    let value = attribute(a, &o.value_attribute, m)?;
    let (unit, unit_name) = if o.check_units {
        (
            attribute(a, &o.unit_accession_attribute, m)?,
            attribute(a, &o.unit_name_attribute, m)?,
        )
    } else {
        (None, None)
    };
    Ok(ParsedCVTerm {
        accession,
        name,
        has_value: value.is_some(),
        value: m.copy(value.unwrap_or_default())?,
        has_unit_accession: unit.is_some(),
        unit_accession: m.copy(unit.unwrap_or_default())?,
        has_unit_name: unit_name.is_some(),
        unit_name: m.copy(unit_name.unwrap_or_default())?,
    })
}
struct Group<'a> {
    rules: Vec<&'a CVMappingRule>,
    counts: BTreeMap<(&'a str, &'a str), u32>,
    maximum_key: usize,
}
struct Index<'a> {
    paths: BTreeMap<&'a str, Group<'a>>,
    maximum_path: usize,
}
fn comparison(n: usize, bytes: usize) -> usize {
    ((usize::BITS - n.max(1).leading_zeros()) as usize)
        .saturating_mul(12)
        .saturating_mul(bytes.saturating_add(1))
}
impl<'a> Index<'a> {
    fn new(mapping: &'a CVMappings, o: &ValidationLimits, m: &mut Meter) -> Result<Self> {
        let rules = &mapping.mapping_rules;
        m.cap(rules.len(), o.max_rules)?;
        m.spend(rules.len(), 0)?;
        let mut count = 0usize;
        let mut maximum_path = 0usize;
        for r in rules {
            count = count.checked_add(r.terms.len()).ok_or_else(|| m.limit())?;
            m.cap(count, o.max_mapping_terms)?;
            maximum_path = maximum_path.max(r.element_path.len());
        }
        m.tree::<(&str, Group<'_>)>(rules.len())?;
        m.spend(count, 0)?;
        let mut paths: BTreeMap<&str, Group<'a>> = BTreeMap::new();
        for r in rules {
            m.spend(comparison(rules.len(), maximum_path), 0)?;
            let group = paths.entry(&r.element_path).or_insert_with(|| Group {
                rules: Vec::new(),
                counts: BTreeMap::new(),
                maximum_key: 0,
            });
            for t in &r.terms {
                group.maximum_key = group
                    .maximum_key
                    .max(r.identifier.len().saturating_add(t.accession.len()));
            }
            m.push(&mut group.rules, r)?;
        }
        Ok(Self {
            paths,
            maximum_path,
        })
    }
    fn get(&mut self, key: &str, m: &mut Meter) -> Result<Option<&mut Group<'a>>> {
        m.spend(
            comparison(self.paths.len(), self.maximum_path.min(key.len())),
            0,
        )?;
        Ok(self.paths.get_mut(key))
    }
    fn close(
        &mut self,
        key: &str,
        parent: &str,
        o: &ValidationOptions,
        report: &mut ValidationReport,
        m: &mut Meter,
    ) -> Result<()> {
        let Some(group) = self.get(key, m)? else {
            return Ok(());
        };
        m.spend(group.rules.len().saturating_mul(2), 0)?;
        for r in &group.rules {
            m.spend(r.terms.len(), 0)?;
            for t in &r.terms {
                if !t.is_repeatable && group.count(&r.identifier, &t.accession, m)? > 1 {
                    emit(
                        report,
                        true,
                        o,
                        m,
                        r.identifier.len().saturating_add(parent.len()),
                        format_args!(
                            "Violated mapping rule '{}' number of term repeats at element '{parent}'",
                            r.identifier
                        ),
                    )?;
                }
            }
        }
        for r in &group.rules {
            let mut matched = 0usize;
            for t in &r.terms {
                if group.count(&r.identifier, &t.accession, m)? > 0 {
                    matched += 1;
                }
            }
            let total = r.terms.len();
            let suffix = match (r.requirement_level, r.combinations_logic) {
                (Level::Must, Logic::And) if matched != total => Some(0),
                (Level::Must, Logic::Or) if matched == 0 => Some(1),
                (Level::Must, Logic::Xor) if matched != 1 => Some(2),
                (Level::May, Logic::And) if matched != 0 && matched != total => Some(3),
                (Level::May, Logic::Xor) if matched > 1 => Some(4),
                _ => None,
            };
            if let Some(which) = suffix {
                m.cap(
                    report
                        .errors
                        .len()
                        .saturating_add(report.warnings.len())
                        .saturating_add(1),
                    o.limits.max_diagnostics,
                )?;
                m.spend(256, 256)?;
                let detail = match which {
                    0 => format!(", {total} term(s) should be present, {matched} found!"),
                    1 => ", at least one term must be present!".into(),
                    2 => " exactly one of the allowed terms must be used!".into(),
                    3 => ", if any, all terms must be given!".into(),
                    _ => ", if any, only exactly one of the allowed terms can be used!".into(),
                };
                // The small fixed-size requirement suffix precedes the full diagnostic.
                emit(
                    report,
                    true,
                    o,
                    m,
                    r.identifier
                        .len()
                        .saturating_add(parent.len())
                        .saturating_add(detail.len()),
                    format_args!(
                        "Violated mapping rule '{}' at element '{parent}'{detail}",
                        r.identifier
                    ),
                )?;
            }
        }
        m.spend(group.counts.len(), 0)?;
        group.counts.clear();
        Ok(())
    }
}
impl<'a> Group<'a> {
    fn count(&self, rule: &str, term: &str, m: &mut Meter) -> Result<u32> {
        m.spend(comparison(self.counts.len(), self.maximum_key), 0)?;
        Ok(self.counts.get(&(rule, term)).copied().unwrap_or(0))
    }
    fn increment(&mut self, rule: &'a str, term: &'a str, m: &mut Meter) -> Result<()> {
        m.spend(comparison(self.counts.len(), self.maximum_key), 0)?;
        // Charging a sparse node for each newly encountered key is conservative.
        if !self.counts.contains_key(&(rule, term)) {
            m.tree::<((&str, &str), u32)>(1)?;
        }
        m.spend(comparison(self.counts.len(), self.maximum_key), 0)?;
        let count = self.counts.entry((rule, term)).or_default();
        *count = count.checked_add(1).ok_or_else(|| m.limit())?;
        Ok(())
    }
}
fn emit(
    report: &mut ValidationReport,
    error: bool,
    o: &ValidationOptions,
    m: &mut Meter,
    variable_bytes: usize,
    args: fmt::Arguments<'_>,
) -> Result<()> {
    m.cap(
        report
            .errors
            .len()
            .saturating_add(report.warnings.len())
            .saturating_add(1),
        o.limits.max_diagnostics,
    )?;
    let n = variable_bytes.saturating_add(512);
    m.spend(n.saturating_mul(2), n.saturating_mul(2))?;
    let text = fmt::format(args);
    m.push(
        if error {
            &mut report.errors
        } else {
            &mut report.warnings
        },
        text,
    )
}
fn matches_term(
    cv: &ControlledVocabulary,
    base: &str,
    exact: bool,
    children: bool,
    accession: &str,
    m: &mut Meter,
) -> Result<bool> {
    m.spend(base.len().min(accession.len()).saturating_add(1), 0)?;
    Ok((exact && base == accession)
        || (children
            && cv.has_descendant_with_budget(base, accession, &mut m.work, &mut m.bytes)?))
}
#[allow(clippy::too_many_arguments)] // One private source callback carries its borrowed operation context.
fn handle(
    index: &mut Index<'_>,
    cv: &ControlledVocabulary,
    p: &ParsedCVTerm,
    key: &str,
    parent: &str,
    o: &ValidationOptions,
    report: &mut ValidationReport,
    m: &mut Meter,
) -> Result<()> {
    let Some(term) = encounter(cv, p, parent, o, report, m)? else {
        return Ok(());
    };
    handle_known(index, cv, term, p, key, parent, o, report, m)
}
fn encounter<'a>(
    cv: &'a ControlledVocabulary,
    p: &ParsedCVTerm,
    parent: &str,
    o: &ValidationOptions,
    report: &mut ValidationReport,
    m: &mut Meter,
) -> Result<Option<&'a CVTermDefinition>> {
    let n = p
        .accession
        .len()
        .saturating_add(p.name.len())
        .saturating_add(parent.len());
    let Some(term) = cv.find_term_with_budget(&p.accession, &mut m.work, &mut m.bytes)? else {
        emit(
            report,
            false,
            o,
            m,
            n,
            format_args!(
                "Unknown CV term: '{} - {}' at element '{parent}'",
                p.accession, p.name
            ),
        )?;
        return Ok(None);
    };
    if term.obsolete {
        emit(
            report,
            false,
            o,
            m,
            n,
            format_args!(
                "Obsolete CV term: '{} - {}' at element '{parent}'",
                p.accession, p.name
            ),
        )?;
    }
    Ok(Some(term))
}
#[allow(clippy::too_many_arguments)] // Shared exact rule callback for the two concrete validators.
fn handle_known(
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
            .saturating_mul(2),
        0,
    )?;
    let (parsed, correct) = (trim(&p.name), trim(&term.name));
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
    if o.check_term_value_types {
        value(term, p, parent, o, report, m)?;
    }
    Ok(())
}
fn units(
    cv: &ControlledVocabulary,
    term: &CVTermDefinition,
    p: &ParsedCVTerm,
    o: &ValidationOptions,
    report: &mut ValidationReport,
    m: &mut Meter,
) -> Result<()> {
    let n = p
        .accession
        .len()
        .saturating_add(p.name.len())
        .saturating_add(p.unit_accession.len())
        .saturating_add(p.unit_name.len());
    if term.units.is_empty() {
        if p.has_unit_accession || p.has_unit_name {
            emit(
                report,
                false,
                o,
                m,
                n,
                format_args!(
                    "Unit CV term used, but not allowed: {} - {} of term {} - {}",
                    p.unit_accession, p.unit_name, p.accession, p.name
                ),
            )?;
        }
    } else if !p.has_unit_accession {
        emit(
            report,
            true,
            o,
            m,
            n,
            format_args!("CV term must have a unit: {} - {}", p.accession, p.name),
        )?;
    } else if cv
        .find_term_with_budget(&p.unit_accession, &mut m.work, &mut m.bytes)?
        .is_none()
    {
        emit(
            report,
            true,
            o,
            m,
            n,
            format_args!(
                "Unit CV term not found: {} - {} of term {} - {}",
                p.unit_accession, p.unit_name, p.accession, p.name
            ),
        )?;
    } else {
        let mut allowed = false;
        for unit in &term.units {
            m.spend(unit.len().min(p.unit_accession.len()).saturating_add(1), 0)?;
            if unit == &p.unit_accession {
                allowed = true;
                break;
            }
        }
        if !allowed {
            for unit in &term.units {
                // CPP-039: compare the descendant to the supplied unit, not the measurement.
                if cv.has_descendant_with_budget(
                    unit,
                    &p.unit_accession,
                    &mut m.work,
                    &mut m.bytes,
                )? {
                    allowed = true;
                    break;
                }
            }
        }
        if !allowed {
            emit(
                report,
                true,
                o,
                m,
                n,
                format_args!(
                    "Unit CV term not allowed: {} - {} of term {} - {}",
                    p.unit_accession, p.unit_name, p.accession, p.name
                ),
            )?;
        }
    }
    Ok(())
}
fn value(
    term: &CVTermDefinition,
    p: &ParsedCVTerm,
    parent: &str,
    o: &ValidationOptions,
    report: &mut ValidationReport,
    m: &mut Meter,
) -> Result<()> {
    let ty = term.xref_type;
    let n = p
        .accession
        .len()
        .saturating_add(p.name.len())
        .saturating_add(p.value.len())
        .saturating_add(parent.len());
    m.spend(
        p.value.len().saturating_mul(8).saturating_add(256),
        p.value.len().saturating_mul(2).saturating_add(256),
    )?;
    if !p.has_value || (p.value.is_empty() && ty != XRefType::String) {
        if ty != XRefType::None {
            emit(
                report,
                true,
                o,
                m,
                n,
                format_args!(
                    "Value-type required, but not given ({}): '{} - {}' value='{}' at element '{parent}'",
                    ty.name(),
                    p.accession,
                    p.name,
                    p.value
                ),
            )?;
        }
        return Ok(());
    }
    let ok = match ty {
        XRefType::None => {
            if !p.accession.starts_with("PATO:") {
                emit(
                    report,
                    true,
                    o,
                    m,
                    n,
                    format_args!(
                        "Value of CV term not allowed: '{} - {}' value='{}' at element '{parent}'",
                        p.accession, p.name, p.value
                    ),
                )?;
            }
            return Ok(());
        }
        XRefType::String => true,
        XRefType::Integer => i32::from_list_item(&p.value).is_ok(),
        XRefType::Decimal => f64::from_list_item(&p.value).is_ok(),
        XRefType::NegativeInteger => i32::from_list_item(&p.value).is_ok_and(|n| n < 0),
        XRefType::PositiveInteger => i32::from_list_item(&p.value).is_ok_and(|n| n > 0),
        XRefType::NonNegativeInteger => i32::from_list_item(&p.value).is_ok_and(|n| n >= 0),
        XRefType::NonPositiveInteger => i32::from_list_item(&p.value).is_ok_and(|n| n <= 0),
        XRefType::Boolean => {
            let v = trim(&p.value);
            v == "1"
                || v == "0"
                || v.eq_ignore_ascii_case("true")
                || v.eq_ignore_ascii_case("false")
        }
        XRefType::Date => DateTime::parse(&p.value).is_ok(),
        XRefType::AnyUri => p.value.contains(':'),
    };
    if !ok {
        if ty == XRefType::AnyUri {
            emit(
                report,
                true,
                o,
                m,
                n,
                format_args!(
                    "Value-type of CV term wrong, should be xsd:anyURI (at least a colon is needed): '{} - {}' value={}' at element '{parent}'",
                    p.accession, p.name, p.value
                ),
            )?;
        } else {
            emit(
                report,
                true,
                o,
                m,
                n,
                format_args!(
                    "Value-type of CV term wrong, should be {}: '{} - {}' value='{}' at element '{parent}'",
                    ty.name(),
                    p.accession,
                    p.name,
                    p.value
                ),
            )?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn diagnostic_count_and_payload_fail_before_publication() {
        let mut report = ValidationReport::default();
        let mut options = ValidationOptions::default();
        let large = "x".repeat(100_000);
        let mut m = Meter::new(1_000_000, 64, "test");
        assert!(
            emit(
                &mut report,
                true,
                &options,
                &mut m,
                large.len(),
                format_args!("{large}")
            )
            .is_err()
        );
        assert!(report.errors.is_empty());
        options.limits.max_diagnostics = 0;
        let mut m = options.limits.meter();
        let before = m.bytes;
        assert!(
            emit(
                &mut report,
                true,
                &options,
                &mut m,
                large.len(),
                format_args!("{large}")
            )
            .is_err()
        );
        assert_eq!(m.bytes, before);
    }
    #[test]
    fn checked_source_u32_counter_does_not_wrap() {
        let mut group = Group {
            rules: Vec::new(),
            counts: BTreeMap::from([(("rule", "term"), u32::MAX)]),
            maximum_key: 8,
        };
        let mut m = ValidationLimits::default().meter();
        assert!(group.increment("rule", "term", &mut m).is_err());
        assert_eq!(group.counts[&("rule", "term")], u32::MAX);
    }
    #[test]
    fn graph_query_resource_failure_is_not_an_unknown_term_warning() {
        let mut cv = ControlledVocabulary::new();
        cv.load_obo_reader("test", b"[Term]\nid: A\nname: alpha\n".as_slice())
            .unwrap();
        let o = ValidationOptions::default();
        let mapping = CVMappings::default();
        let mut initial = o.limits.meter();
        let mut index = Index::new(&mapping, &o.limits, &mut initial).unwrap();
        let mut report = ValidationReport::default();
        let mut exhausted = Meter::new(0, 1_000_000, "test");
        let p = ParsedCVTerm {
            accession: "unknown".into(),
            ..Default::default()
        };
        assert!(
            handle(
                &mut index,
                &cv,
                &p,
                "/r/cvParam/@accession",
                "/r",
                &o,
                &mut report,
                &mut exhausted
            )
            .is_err()
        );
        assert!(report.warnings.is_empty());
    }
    #[test]
    fn cumulative_paths_counts_and_diagnostics_exhaust_after_successful_terms() {
        let mut cv = ControlledVocabulary::new();
        cv.load_obo_reader("test", b"[Term]\nid: A\nname: alpha\n".as_slice())
            .unwrap();
        let mut mapping = CVMappings::default();
        mapping.mapping_rules.push(CVMappingRule {
            identifier: "rule".into(),
            element_path: "/r/cvParam/@accession".into(),
            terms: vec![crate::data_structures::CVMappingTerm {
                accession: "A".into(),
                use_term: true,
                ..Default::default()
            }],
            ..Default::default()
        });
        let o = ValidationOptions::default();
        let mut m = Meter::new(1_000_000, 20_000, "test");
        let mut index = Index::new(&mapping, &o.limits, &mut m).unwrap();
        let p = ParsedCVTerm {
            accession: "A".into(),
            name: "incorrect".into(),
            ..Default::default()
        };
        let mut report = ValidationReport::default();
        let mut successful = 0;
        let failure = loop {
            let step = (|| {
                let parent = path(&[], Some("r"), &mut m)?;
                let key = cv_path(&parent, &o, &mut m)?;
                handle(&mut index, &cv, &p, &key, &parent, &o, &mut report, &mut m)
            })();
            match step {
                Ok(()) => {
                    successful += 1;
                    assert!(successful < 100);
                }
                Err(e) => break e,
            }
        };
        assert!(failure.to_string().contains("limit"));
        assert!(successful > 1);
        assert_eq!(report.errors.len(), successful);
        assert!(index.paths["/r/cvParam/@accession"].counts[&("rule", "A")] >= successful as u32);
        // Existing scientific diagnostics remain coherent private state; callers
        // receive Err rather than any partially accumulated ValidationReport.
        assert!(
            report
                .errors
                .iter()
                .all(|e| e.starts_with("Name of CV term not correct"))
        );
    }
}
