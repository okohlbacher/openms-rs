// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Bounded source CV mapping XML loading; semantic validation is separate.

use super::cv_xml::{Attributes, Element, Meter, XmlLimits, document, scan_elements};
use crate::data_structures::{
    CVMappingRule, CVMappingTerm, CVMappings, CVReference, CombinationsLogic, RequirementLevel,
};
use crate::{Error, Result};
use std::{io::BufRead, path::Path};

#[derive(Clone, Copy, Debug)]
pub struct ReadOptions {
    pub strip_namespaces: bool,
    /// Maximum compressed-stream decoded input bytes and normalized UTF-8 bytes.
    pub max_input_bytes: usize,
    pub max_depth: usize,
    pub max_elements: usize,
    pub max_references: usize,
    pub max_rules: usize,
    pub max_terms: usize,
    pub max_work: usize,
    pub max_bytes: usize,
}
impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            strip_namespaces: false,
            max_input_bytes: 16 * 1024 * 1024,
            max_depth: 128,
            max_elements: 1_000_000,
            max_references: 100_000,
            max_rules: 100_000,
            max_terms: 1_000_000,
            max_work: 50_000_000,
            max_bytes: 128 * 1024 * 1024,
        }
    }
}

/// Each operation has local parser state. Reusing this loader after an error is safe.
#[derive(Clone, Copy, Debug, Default)]
pub struct CVMappingFile {
    pub options: ReadOptions,
}
impl CVMappingFile {
    pub fn load(&self, path: impl AsRef<Path>) -> Result<CVMappings> {
        self.read(super::path_io::open(path.as_ref())?)
    }
    pub fn read(&self, input: impl BufRead) -> Result<CVMappings> {
        let mut output = CVMappings::default();
        self.read_into(input, &mut output)?;
        Ok(output)
    }
    pub fn load_into(&self, path: impl AsRef<Path>, output: &mut CVMappings) -> Result<()> {
        self.read_into(super::path_io::open(path.as_ref())?, output)
    }
    /// Atomically append parsed references and replace rules, matching successful
    /// source loading into an existing destination. An error leaves output intact.
    pub fn read_into(&self, input: impl BufRead, output: &mut CVMappings) -> Result<()> {
        let mut meter = Meter::new(self.options.max_work, self.options.max_bytes, "CV mapping");
        let text = document(input, &self.options.xml_limits(), &mut meter)?;
        let (references, rules) = parse(&text, &self.options, &mut meter)?;
        let count = output
            .cv_references()
            .len()
            .checked_add(references.len())
            .ok_or_else(limit)?;
        cap(count, self.options.max_references)?;
        // Only old references survive. Existing rules need not be cloned.
        meter.spend(count, 0)?;
        meter.slots::<CVReference>(count.saturating_mul(2))?;
        meter.tree::<String>(count)?;
        let mut joined = Vec::with_capacity(count);
        for r in output.cv_references() {
            joined.push(CVReference {
                name: meter.copy(&r.name)?,
                identifier: meter.copy(&r.identifier)?,
            });
        }
        joined.extend(references);
        // set_cv_references clones index keys even for a repeated identifier.
        // Bound all comparisons by the maximum key length before building it.
        let max_key = joined.iter().map(|r| r.identifier.len()).max().unwrap_or(0);
        let comparisons = (usize::BITS - count.max(1).leading_zeros()) as usize * 12;
        for r in &joined {
            meter.spend(
                comparisons.saturating_mul(max_key.saturating_add(1)),
                r.identifier.len(),
            )?;
        }
        let mut next = CVMappings::default();
        next.mapping_rules = rules;
        next.set_cv_references(joined);
        *output = next;
        Ok(())
    }
}

fn bad(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}
fn limit() -> Error {
    bad("CV mapping resource limit exceeded")
}
fn cap(value: usize, maximum: usize) -> Result<()> {
    if value > maximum {
        Err(limit())
    } else {
        Ok(())
    }
}
fn get<'a>(a: &'a Attributes, key: &str) -> Option<&'a str> {
    a.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_str())
}
fn required<'a>(a: &'a Attributes, key: &str) -> Result<&'a str> {
    get(a, key).ok_or_else(|| bad(format!("missing CV mapping attribute {key}")))
}
fn boolean(s: &str) -> Result<bool> {
    match s {
        "true" => Ok(true),
        "false" => Ok(false),
        _ => Err(bad("CV mapping boolean requires true or false")),
    }
}

/// Correct source namespace stripping (CPP-032). The source intentionally strips
/// only element_path, never scope_path. Slash compaction remains source behavior.
fn strip_path(path: &str, m: &mut Meter) -> Result<String> {
    m.spend(path.len().saturating_mul(3), path.len().saturating_add(1))?;
    let mut out = String::with_capacity(path.len().saturating_add(1));
    for segment in path.split('/').filter(|s| !s.is_empty()) {
        out.push('/');
        if let Some((prefix, local)) = segment.split_once(':') {
            if local.contains(':') {
                return Err(bad("multiple namespace colons in mapping path"));
            }
            if prefix.starts_with('@') {
                out.push('@');
            }
            out.push_str(local);
        } else {
            out.push_str(segment);
        }
    }
    Ok(out)
}
struct Draft {
    references: Vec<CVReference>,
    rules: Vec<CVMappingRule>,
    rule: CVMappingRule,
    terms: usize,
}
impl Draft {
    fn start(&mut self, tag: &str, a: &Attributes, o: &ReadOptions, m: &mut Meter) -> Result<()> {
        // Every attribute lookup is a linear borrowed search, bounded before use.
        m.spend(a.len().saturating_mul(16).saturating_add(1), 0)?;
        match tag {
            "CvReference" => {
                cap(self.references.len().saturating_add(1), o.max_references)?;
                let r = CVReference {
                    name: m.copy(required(a, "cvName")?)?,
                    identifier: m.copy(required(a, "cvIdentifier")?)?,
                };
                m.push(&mut self.references, r)?;
            }
            "CvMappingRule" => {
                self.rule.identifier = m.copy(required(a, "id")?)?;
                let path = required(a, "cvElementPath")?;
                self.rule.element_path = if o.strip_namespaces {
                    strip_path(path, m)?
                } else {
                    m.copy(path)?
                };
                self.rule.requirement_level = match required(a, "requirementLevel")? {
                    "MAY" => RequirementLevel::May,
                    "SHOULD" => RequirementLevel::Should,
                    _ => RequirementLevel::Must,
                };
                self.rule.scope_path = m.copy(required(a, "scopePath")?)?;
                self.rule.combinations_logic = match required(a, "cvTermsCombinationLogic")? {
                    "AND" => CombinationsLogic::And,
                    "XOR" => CombinationsLogic::Xor,
                    _ => CombinationsLogic::Or,
                };
            }
            "CvTerm" => {
                self.terms = self.terms.checked_add(1).ok_or_else(limit)?;
                cap(self.terms, o.max_terms)?;
                let term = CVMappingTerm {
                    accession: m.copy(required(a, "termAccession")?)?,
                    use_term: boolean(required(a, "useTerm")?)?,
                    use_term_name: match get(a, "useTermName") {
                        None | Some("") => false,
                        Some(v) => boolean(v)?,
                    },
                    term_name: m.copy(required(a, "termName")?)?,
                    is_repeatable: match get(a, "isRepeatable") {
                        None | Some("") => true,
                        Some(v) => boolean(v)?,
                    },
                    allow_children: boolean(required(a, "allowChildren")?)?,
                    cv_identifier_ref: m.copy(required(a, "cvIdentifierRef")?)?,
                };
                m.push(&mut self.rule.terms, term)?;
            }
            _ => {}
        }
        Ok(())
    }
    fn end(&mut self, tag: &str, o: &ReadOptions, m: &mut Meter) -> Result<()> {
        if tag == "CvMappingRule" {
            cap(self.rules.len().saturating_add(1), o.max_rules)?;
            let r = std::mem::take(&mut self.rule);
            m.push(&mut self.rules, r)?;
        }
        Ok(())
    }
}
fn parse(
    text: &str,
    o: &ReadOptions,
    m: &mut Meter,
) -> Result<(Vec<CVReference>, Vec<CVMappingRule>)> {
    let mut draft = Draft {
        references: Vec::new(),
        rules: Vec::new(),
        rule: CVMappingRule::default(),
        terms: 0,
    };
    scan_elements(text, &o.xml_limits(), m, |event, _, m| match event {
        Element::Start(tag, a) => draft.start(tag, a, o, m),
        Element::End(tag) => draft.end(tag, o, m),
    })?;
    Ok((draft.references, draft.rules))
}
impl ReadOptions {
    fn xml_limits(&self) -> XmlLimits {
        XmlLimits {
            max_input_bytes: self.max_input_bytes,
            max_depth: self.max_depth,
            max_elements: self.max_elements,
        }
    }
}
