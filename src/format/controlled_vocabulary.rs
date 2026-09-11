// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Owned ontology definitions and the complete source OBO projection.
//! This is distinct from metadata CV parameter instances. Successful cumulative
//! loads preserve the source's aliases/child-index behavior; checked failures
//! are atomic. Unknown OBO fields remain in each definition's `unparsed` list.

use crate::metadata::{MetaValue, MetaValueData};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, BufReader, Write};
use std::mem::size_of;
use std::path::Path;
use std::sync::OnceLock;

mod obo;

/// Value types recognized by the source ControlledVocabulary provider.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum XRefType {
    String,
    Integer,
    Decimal,
    NegativeInteger,
    PositiveInteger,
    NonNegativeInteger,
    NonPositiveInteger,
    Boolean,
    Date,
    AnyUri,
    #[default]
    None,
}
impl XRefType {
    pub const ALL: [Self; 11] = [
        Self::String,
        Self::Integer,
        Self::Decimal,
        Self::NegativeInteger,
        Self::PositiveInteger,
        Self::NonNegativeInteger,
        Self::NonPositiveInteger,
        Self::Boolean,
        Self::Date,
        Self::AnyUri,
        Self::None,
    ];
    pub const fn name(self) -> &'static str {
        match self {
            Self::String => "xsd:string",
            Self::Integer => "xsd:integer",
            Self::Decimal => "xsd:decimal",
            Self::NegativeInteger => "xsd:negativeInteger",
            Self::PositiveInteger => "xsd:positiveInteger",
            Self::NonNegativeInteger => "xsd:nonNegativeInteger",
            Self::NonPositiveInteger => "xsd:nonPositiveInteger",
            Self::Boolean => "xsd:boolean",
            Self::Date => "xsd:date",
            Self::AnyUri => "xsd:anyURI",
            Self::None => "none",
        }
    }
}

/// Complete ontology definition. Ordinary field access/Clone/equality have
/// ordinary Rust costs. Use consuming operations' limits for untrusted data.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CVTermDefinition {
    pub name: String,
    pub id: String,
    pub parents: BTreeSet<String>,
    pub children: BTreeSet<String>,
    pub obsolete: bool,
    pub description: String,
    pub synonyms: Vec<String>,
    pub unparsed: Vec<String>,
    pub xref_type: XRefType,
    pub xref_binary: Vec<String>,
    pub units: BTreeSet<String>,
}
impl CVTermDefinition {
    pub fn is_higher_better_score(&self) -> bool {
        !self
            .unparsed
            .iter()
            .any(|line| line.starts_with("relationship: has_order MS:1002109"))
    }
    /// Source string overload: an empty value omits the value attribute.
    pub fn to_xml(&self, cv_ref: &str, value: &str) -> Result<String> {
        self.to_xml_with_limits(cv_ref, value, VocabularyLimits::default())
    }
    pub fn to_xml_with_limits(
        &self,
        cv_ref: &str,
        value: &str,
        limits: VocabularyLimits,
    ) -> Result<String> {
        let (mut work, mut bytes) = (limits.max_work, limits.max_bytes);
        let mut out = Text::new(limits.max_output_bytes, &mut work, &mut bytes);
        self.xml_start(cv_ref, &mut out)?;
        if !value.is_empty() {
            out.attribute("value", value)?;
        }
        out.raw("/>")?;
        Ok(out.output)
    }
    /// Typed Empty omits value, but a present empty string emits `value=""`.
    /// Uses the value's actual unit identity. Source's first-allowed-unit
    /// substitution and empty-unit-set dereference are deliberately corrected.
    pub fn to_xml_value(&self, cv_ref: &str, value: &MetaValue) -> Result<String> {
        let limits = VocabularyLimits::default();
        let (mut work, mut bytes) = (limits.max_work, limits.max_bytes);
        self.to_xml_value_with_budget(
            cv_ref,
            value,
            limits.max_output_bytes,
            &mut work,
            &mut bytes,
        )
    }
    pub fn to_xml_value_with_limits(
        &self,
        cv_ref: &str,
        value: &MetaValue,
        limits: VocabularyLimits,
    ) -> Result<String> {
        let (mut work, mut bytes) = (limits.max_work, limits.max_bytes);
        self.to_xml_value_with_budget(
            cv_ref,
            value,
            limits.max_output_bytes,
            &mut work,
            &mut bytes,
        )
    }
    pub(crate) fn to_xml_value_with_budget(
        &self,
        cv_ref: &str,
        value: &MetaValue,
        max_output: usize,
        work: &mut usize,
        bytes: &mut usize,
    ) -> Result<String> {
        let mut out = Text::new(max_output, work, bytes);
        self.xml_start(cv_ref, &mut out)?;
        if !matches!(value.data(), MetaValueData::Empty) {
            out.raw(" value=\"")?;
            match value.data() {
                MetaValueData::String(s) => out.escaped(s)?,
                MetaValueData::Integer(n) => out.integer(*n)?,
                MetaValueData::Float(n) => out.float(*n)?,
                MetaValueData::StringList(v) => {
                    out.meter.spend(v.len())?;
                    out.raw("[")?;
                    for (i, s) in v.iter().enumerate() {
                        if i != 0 {
                            out.raw(", ")?;
                        }
                        out.escaped(s)?;
                    }
                    out.raw("]")?;
                }
                MetaValueData::IntegerList(v) => {
                    out.meter.spend(mul(v.len(), 64)?)?;
                    out.raw("[")?;
                    for (i, n) in v.iter().enumerate() {
                        if i != 0 {
                            out.raw(", ")?;
                        }
                        out.integer(*n)?;
                    }
                    out.raw("]")?;
                }
                MetaValueData::FloatList(v) => {
                    out.meter.spend(mul(v.len(), 1024)?)?;
                    out.raw("[")?;
                    for (i, n) in v.iter().enumerate() {
                        if i != 0 {
                            out.raw(", ")?;
                        }
                        out.float(*n)?;
                    }
                    out.raw("]")?;
                }
                MetaValueData::Empty => unreachable!(),
            }
            out.raw("\"")?;
        }
        if let Some(unit) = value.unit() {
            out.attribute("unitAccession", unit.accession())?;
            out.attribute("unitCvRef", unit.cv_ref())?;
            if !unit.name().is_empty() {
                out.attribute("unitName", unit.name())?;
            }
        }
        out.raw("/>")?;
        Ok(out.output)
    }
    fn xml_start(&self, cv_ref: &str, out: &mut Text<'_>) -> Result<()> {
        out.raw("<cvParam")?;
        out.attribute("accession", &self.id)?;
        out.attribute("cvRef", cv_ref)?;
        out.attribute("name", &self.name)
    }
    fn measure(&self, m: &mut Meter<'_>) -> Result<()> {
        for s in [&self.id, &self.name, &self.description] {
            m.text(s.len())?;
        }
        for set in [&self.parents, &self.children, &self.units] {
            m.tree::<String>(set.len())?;
            for s in set {
                m.text(s.len())?;
            }
        }
        for list in [&self.synonyms, &self.unparsed, &self.xref_binary] {
            m.slots::<String>(list.len())?;
            for s in list {
                m.text(s.len())?;
            }
        }
        Ok(())
    }
}

/// Bounds cover an entire load, graph query, checked copy or render operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct VocabularyLimits {
    pub max_input_bytes: usize,
    pub max_line_bytes: usize,
    pub max_terms: usize,
    pub max_entries: usize,
    pub max_work: usize,
    pub max_bytes: usize,
    pub max_output_bytes: usize,
}
impl Default for VocabularyLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 16 * 1024 * 1024,
            max_line_bytes: 1024 * 1024,
            max_terms: 100_000,
            max_entries: 1_000_000,
            // Includes all five source loads and their deliberately repeated
            // index rebuilding; comparison-byte bounds are conservative.
            max_work: 1_000_000_000,
            max_bytes: 256 * 1024 * 1024,
            max_output_bytes: 16 * 1024 * 1024,
        }
    }
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum OboEncoding {
    #[default]
    Utf8,
    Windows1252,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OboDiagnostic {
    pub line: usize,
    pub message: String,
}
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct OboLoadReport {
    /// Source term stanzas committed, including repeated IDs; excludes placeholders.
    pub definitions: usize,
    pub diagnostics: Vec<OboDiagnostic>,
}

/// Immutable borrowed query surface over source-compatible cumulative indexes.
#[derive(Clone, Debug)]
pub struct ControlledVocabulary {
    terms: BTreeMap<String, CVTermDefinition>,
    names: BTreeMap<String, String>,
    name: String,
    label: String,
    version: String,
    url: String,
    limits: VocabularyLimits,
    max_id_bytes: usize,
    max_name_bytes: usize,
}
impl Default for ControlledVocabulary {
    fn default() -> Self {
        Self::with_limits(VocabularyLimits::default())
    }
}
impl ControlledVocabulary {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn with_limits(limits: VocabularyLimits) -> Self {
        Self {
            terms: BTreeMap::new(),
            names: BTreeMap::new(),
            name: String::new(),
            label: String::new(),
            version: String::new(),
            url: String::new(),
            limits,
            max_id_bytes: 0,
            max_name_bytes: 0,
        }
    }
    pub fn limits(&self) -> VocabularyLimits {
        self.limits
    }
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn label(&self) -> &str {
        &self.label
    }
    pub fn version(&self) -> &str {
        &self.version
    }
    pub fn url(&self) -> &str {
        &self.url
    }
    pub fn terms(&self) -> &BTreeMap<String, CVTermDefinition> {
        &self.terms
    }
    pub fn exists(&self, id: &str) -> bool {
        self.terms.contains_key(id)
    }
    pub fn has_term_with_name(&self, name: &str) -> bool {
        self.names.contains_key(name)
    }
    pub fn get_term(&self, id: &str) -> Result<&CVTermDefinition> {
        self.terms
            .get(id)
            .ok_or_else(|| invalid("unknown CV identifier"))
    }
    pub fn find_term_by_name(&self, name: &str) -> Option<&CVTermDefinition> {
        self.names.get(name).and_then(|id| self.terms.get(id))
    }
    pub fn get_term_by_name(&self, name: &str, description: &str) -> Result<&CVTermDefinition> {
        let (mut work, mut bytes) = (self.limits.max_work, self.limits.max_bytes);
        self.get_term_by_name_with_budget(name, description, &mut work, &mut bytes)
    }
    #[cfg(feature = "semantic-validation")]
    pub(crate) fn find_term_with_budget(
        &self,
        id: &str,
        work: &mut usize,
        bytes: &mut usize,
    ) -> Result<Option<&CVTermDefinition>> {
        Meter { work, bytes }.lookup(self.terms.len(), self.max_id_bytes, id.len())?;
        Ok(self.terms.get(id))
    }
    #[cfg(feature = "semantic-validation")]
    pub(crate) fn has_descendant_with_budget(
        &self,
        parent: &str,
        child: &str,
        work: &mut usize,
        bytes: &mut usize,
    ) -> Result<bool> {
        self.walk(
            parent,
            false,
            &mut |id, m| {
                m.spend(add(id.len().min(child.len()), 1)?)?;
                Ok(id == child)
            },
            &mut Meter { work, bytes },
        )
    }
    pub(crate) fn get_term_with_budget(
        &self,
        id: &str,
        work: &mut usize,
        bytes: &mut usize,
    ) -> Result<&CVTermDefinition> {
        let mut m = Meter { work, bytes };
        m.lookup(self.terms.len(), self.max_id_bytes, id.len())?;
        self.get_term(id)
    }
    pub(crate) fn find_term_by_name_with_budget(
        &self,
        name: &str,
        work: &mut usize,
        bytes: &mut usize,
    ) -> Result<Option<&CVTermDefinition>> {
        let mut m = Meter { work, bytes };
        m.lookup(self.names.len(), self.max_name_bytes, name.len())?;
        match self.names.get(name) {
            None => Ok(None),
            Some(id) => {
                m.lookup(self.terms.len(), self.max_id_bytes, id.len())?;
                Ok(self.terms.get(id))
            }
        }
    }
    pub(crate) fn get_term_by_name_with_budget(
        &self,
        name: &str,
        description: &str,
        work: &mut usize,
        bytes: &mut usize,
    ) -> Result<&CVTermDefinition> {
        if let Some(term) = self.find_term_by_name_with_budget(name, work, bytes)? {
            return Ok(term);
        }
        if description.is_empty() {
            return Err(invalid("unknown CV name"));
        }
        let n = add(name.len(), description.len())?;
        Meter { work, bytes }.text(n)?;
        let mut key = String::new();
        key.try_reserve_exact(n).map_err(|_| limit())?;
        key.push_str(name);
        key.push_str(description);
        self.find_term_by_name_with_budget(&key, work, bytes)?
            .ok_or_else(|| invalid("unknown CV name"))
    }
    pub fn checked_clone(&self) -> Result<Self> {
        let (mut work, mut bytes) = (self.limits.max_work, self.limits.max_bytes);
        self.checked_clone_with_budget(&mut work, &mut bytes)
    }
    pub(crate) fn checked_clone_with_budget(
        &self,
        work: &mut usize,
        bytes: &mut usize,
    ) -> Result<Self> {
        self.measure(&mut Meter { work, bytes })?;
        Ok(self.clone())
    }
    fn measure(&self, m: &mut Meter<'_>) -> Result<()> {
        m.slots::<Self>(1)?;
        for s in [&self.name, &self.label, &self.version, &self.url] {
            m.text(s.len())?;
        }
        m.tree::<(String, CVTermDefinition)>(self.terms.len())?;
        for (id, term) in &self.terms {
            m.text(id.len())?;
            term.measure(m)?;
        }
        m.tree::<(String, String)>(self.names.len())?;
        for (key, id) in &self.names {
            m.text(key.len())?;
            m.text(id.len())?;
        }
        Ok(())
    }
    pub fn load_obo(&mut self, name: &str, path: impl AsRef<Path>) -> Result<OboLoadReport> {
        self.load_obo_encoded(
            name,
            BufReader::new(std::fs::File::open(path)?),
            OboEncoding::Utf8,
        )
    }
    pub fn load_obo_reader(&mut self, name: &str, reader: impl BufRead) -> Result<OboLoadReport> {
        self.load_obo_encoded(name, reader, OboEncoding::Utf8)
    }
    pub fn load_obo_encoded(
        &mut self,
        name: &str,
        reader: impl BufRead,
        encoding: OboEncoding,
    ) -> Result<OboLoadReport> {
        let (mut work, mut bytes) = (self.limits.max_work, self.limits.max_bytes);
        self.load_with_budget(name, reader, encoding, &mut work, &mut bytes)
    }
    fn load_with_budget(
        &mut self,
        name: &str,
        reader: impl BufRead,
        encoding: OboEncoding,
        work: &mut usize,
        bytes: &mut usize,
    ) -> Result<OboLoadReport> {
        let mut draft = self.checked_clone_with_budget(work, bytes)?;
        let report = obo::load(
            &mut draft,
            name,
            reader,
            encoding,
            &mut Meter { work, bytes },
        )?;
        *self = draft;
        Ok(report)
    }
    pub fn all_child_terms(&self, parent: &str) -> Result<BTreeSet<String>> {
        let mut result = BTreeSet::new();
        self.extend_child_terms(&mut result, parent)?;
        Ok(result)
    }
    pub fn extend_child_terms(&self, output: &mut BTreeSet<String>, parent: &str) -> Result<()> {
        self.extend(output, parent, false)
    }
    pub fn add_all_child_terms(&self, output: &mut BTreeSet<String>, parent: &str) -> Result<()> {
        self.extend(output, parent, true)
    }
    fn extend(
        &self,
        output: &mut BTreeSet<String>,
        parent: &str,
        include_parent: bool,
    ) -> Result<()> {
        let (mut work, mut bytes) = (self.limits.max_work, self.limits.max_bytes);
        self.get_term_with_budget(parent, &mut work, &mut bytes)?;
        let mut m = Meter {
            work: &mut work,
            bytes: &mut bytes,
        };
        m.tree::<String>(output.len())?;
        let mut max = 0;
        for id in output.iter() {
            m.text(id.len())?;
            max = max.max(id.len());
        }
        let mut draft = output.clone();
        let mut insert = |id: &str, m: &mut Meter<'_>| -> Result<bool> {
            m.lookup(draft.len(), max, id.len())?;
            m.tree::<String>(1)?;
            m.text(id.len())?;
            max = max.max(id.len());
            draft.insert(id.into());
            Ok(false)
        };
        if include_parent {
            insert(parent, &mut m)?;
        }
        self.walk(parent, false, &mut insert, &mut m)?;
        *output = draft;
        Ok(())
    }
    pub fn iterate_all_children(
        &self,
        parent: &str,
        mut callback: impl FnMut(&str) -> bool,
    ) -> Result<bool> {
        let (mut work, mut bytes) = (self.limits.max_work, self.limits.max_bytes);
        self.iterate_all_children_with_budget(parent, &mut callback, &mut work, &mut bytes)
    }
    /// The traversal charges ID-map lookup and ID comparison, not arbitrary
    /// callback work. Callbacks that inspect names/payload must meter those too.
    pub(crate) fn iterate_all_children_with_budget(
        &self,
        parent: &str,
        mut callback: impl FnMut(&str) -> bool,
        work: &mut usize,
        bytes: &mut usize,
    ) -> Result<bool> {
        self.walk(
            parent,
            false,
            &mut |id, _| Ok(callback(id)),
            &mut Meter { work, bytes },
        )
    }
    pub fn is_child_of(&self, child: &str, parent: &str) -> Result<bool> {
        let (mut work, mut bytes) = (self.limits.max_work, self.limits.max_bytes);
        self.is_child_of_with_budget(child, parent, &mut work, &mut bytes)
    }
    pub(crate) fn is_child_of_with_budget(
        &self,
        child: &str,
        parent: &str,
        work: &mut usize,
        bytes: &mut usize,
    ) -> Result<bool> {
        self.walk(
            child,
            true,
            &mut |id, m| {
                m.spend(add(id.len().min(parent.len()), 1)?)?;
                Ok(id == parent)
            },
            &mut Meter { work, bytes },
        )
    }
    /// First exact name match in source lexical depth-first descendant order.
    /// The parent itself is not examined. All lookup/name comparisons share
    /// the supplied operation counters; no callback accounting is required.
    pub fn first_child_with_name(
        &self,
        parent: &str,
        name: &str,
    ) -> Result<Option<&CVTermDefinition>> {
        let (mut work, mut bytes) = (self.limits.max_work, self.limits.max_bytes);
        self.first_child_with_name_with_budget(parent, name, &mut work, &mut bytes)
    }
    pub(crate) fn first_child_with_name_with_budget(
        &self,
        parent: &str,
        name: &str,
        work: &mut usize,
        bytes: &mut usize,
    ) -> Result<Option<&CVTermDefinition>> {
        let mut found = None;
        self.walk(
            parent,
            false,
            &mut |id, m| {
                m.lookup(self.terms.len(), self.max_id_bytes, id.len())?;
                let term = self.get_term(id)?;
                m.spend(add(term.name.len().min(name.len()), 1)?)?;
                if term.name == name {
                    found = Some(term);
                    Ok(true)
                } else {
                    Ok(false)
                }
            },
            &mut Meter { work, bytes },
        )?;
        Ok(found)
    }
    fn walk(
        &self,
        root: &str,
        parents: bool,
        callback: &mut impl FnMut(&str, &mut Meter<'_>) -> Result<bool>,
        m: &mut Meter<'_>,
    ) -> Result<bool> {
        type Frame<'a> = (&'a str, std::collections::btree_set::Iter<'a, String>);
        m.lookup(self.terms.len(), self.max_id_bytes, root.len())?;
        let (key, term) = self
            .terms
            .get_key_value(root)
            .ok_or_else(|| invalid("unknown CV identifier"))?;
        let mut stack: Vec<Frame<'_>> = Vec::new();
        let mut active = BTreeSet::<&str>::new();
        m.slots::<Frame<'_>>(1)?;
        m.tree::<&str>(1)?;
        stack.try_reserve_exact(1).map_err(|_| limit())?;
        stack.push((
            key,
            if parents {
                term.parents.iter()
            } else {
                term.children.iter()
            },
        ));
        active.insert(key);
        while let Some((_, children)) = stack.last_mut() {
            m.spend(1)?;
            if let Some(id) = children.next() {
                if callback(id, m)? {
                    return Ok(true);
                }
                m.lookup(active.len(), self.max_id_bytes, id.len())?;
                if active.contains(id.as_str()) {
                    return Err(invalid("cyclic CV graph traversal"));
                }
                m.lookup(self.terms.len(), self.max_id_bytes, id.len())?;
                let term = self.get_term(id)?;
                m.tree::<&str>(1)?;
                if stack.len() == stack.capacity() {
                    let capacity = mul(stack.capacity(), 2)?.max(4);
                    m.slots::<Frame<'_>>(capacity)?;
                    m.spend(stack.len())?;
                    stack
                        .try_reserve_exact(capacity - stack.len())
                        .map_err(|_| limit())?;
                }
                stack.push((
                    id,
                    if parents {
                        term.parents.iter()
                    } else {
                        term.children.iter()
                    },
                ));
                m.lookup(active.len(), self.max_id_bytes, id.len())?;
                active.insert(id);
            } else {
                let (id, _) = stack.pop().unwrap();
                m.lookup(active.len(), self.max_id_bytes, id.len())?;
                active.remove(id);
            }
        }
        Ok(false)
    }
    /// Diagnostic rendering, not a lossless OBO serializer. All lines go to the
    /// selected output; source's accidental stdout parent routing is corrected.
    pub fn to_text(&self) -> Result<String> {
        let (mut work, mut bytes) = (self.limits.max_work, self.limits.max_bytes);
        let mut out = Text::new(self.limits.max_output_bytes, &mut work, &mut bytes);
        out.meter.spend(self.terms.len())?;
        for term in self.terms.values() {
            out.raw("[Term]\nid: '")?;
            out.raw(&term.id)?;
            out.raw("'\nname: '")?;
            out.raw(&term.name)?;
            out.raw("'\n")?;
            out.meter.spend(term.parents.len())?;
            for parent in &term.parents {
                out.raw("is_a: '")?;
                out.raw(parent)?;
                out.raw("'\n")?;
            }
        }
        Ok(out.output)
    }
    pub fn write_text(&self, mut output: impl Write) -> Result<()> {
        let text = self.to_text()?;
        output.write_all(text.as_bytes())?;
        Ok(())
    }
    /// All five pinned source providers, initialized once without filesystem or
    /// network access. Initialization errors are cached and remain recoverable.
    pub fn psi_ms() -> Result<&'static Self> {
        static CV: OnceLock<std::result::Result<ControlledVocabulary, String>> = OnceLock::new();
        match CV.get_or_init(|| {
            let mut result = Self::new();
            let (mut work, mut bytes) = (result.limits.max_work, result.limits.max_bytes);
            for (name, data, encoding) in obo::EMBEDDED {
                result
                    .load_with_budget(name, *data, *encoding, &mut work, &mut bytes)
                    .map_err(|e| e.to_string())?;
            }
            Ok(result)
        }) {
            Ok(cv) => Ok(cv),
            Err(message) => Err(Error::InvalidValue(message.clone())),
        }
    }
}

/// Source FNV-1a byte hash, explicitly 64-bit on every target.
pub fn fnv1a_hash(key: &str) -> u64 {
    key.bytes().fold(14695981039346656037u64, |hash, byte| {
        (hash ^ u64::from(byte)).wrapping_mul(1099511628211)
    })
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn limit() -> Error {
    invalid("controlled vocabulary resource limit exceeded")
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or_else(limit)
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b).ok_or_else(limit)
}
struct Meter<'a> {
    work: &'a mut usize,
    bytes: &'a mut usize,
}
impl Meter<'_> {
    fn spend(&mut self, n: usize) -> Result<()> {
        *self.work = self.work.checked_sub(n).ok_or_else(limit)?;
        Ok(())
    }
    fn allocation(&mut self, n: usize) -> Result<()> {
        *self.bytes = self.bytes.checked_sub(n).ok_or_else(limit)?;
        Ok(())
    }
    fn text(&mut self, n: usize) -> Result<()> {
        self.spend(n)?;
        self.allocation(n)
    }
    fn slots<T>(&mut self, n: usize) -> Result<()> {
        self.spend(n)?;
        self.allocation(mul(n, size_of::<T>())?)
    }
    fn tree<T>(&mut self, n: usize) -> Result<()> {
        if n == 0 {
            return Ok(());
        }
        self.allocation(512)?;
        // A sparse root still reserves eleven element slots, including when
        // each element is a large complete term rather than a small pointer.
        self.slots::<(T, [usize; 4])>(mul(n, 3)?.max(11))
    }
    fn lookup(&mut self, count: usize, max_key: usize, query: usize) -> Result<()> {
        // BTreeMap nodes hold at most 11 keys; height is at most log2(n)+1.
        let height = if count == 0 {
            1
        } else {
            count.ilog2() as usize + 1
        };
        self.spend(mul(mul(12, height)?, add(max_key.min(query), 1)?)?)
    }
}
struct Text<'a> {
    output: String,
    max: usize,
    meter: Meter<'a>,
}
impl<'a> Text<'a> {
    fn new(max: usize, work: &'a mut usize, bytes: &'a mut usize) -> Self {
        Self {
            output: String::new(),
            max,
            meter: Meter { work, bytes },
        }
    }
    fn raw(&mut self, text: &str) -> Result<()> {
        if add(self.output.len(), text.len())? > self.max {
            return Err(limit());
        }
        self.meter.text(text.len())?;
        // Geometric reserve avoids quadratic exact-reservation copies. Charge
        // the full new capacity and its possible copy before allocation.
        if self.output.capacity() - self.output.len() < text.len() {
            let capacity = add(self.output.len(), text.len())?
                .max(self.output.capacity().saturating_mul(2))
                .max(64)
                .min(self.max);
            self.meter.text(capacity)?;
            self.output
                .try_reserve_exact(capacity - self.output.len())
                .map_err(|_| limit())?;
        }
        self.output.push_str(text);
        Ok(())
    }
    fn escaped(&mut self, text: &str) -> Result<()> {
        self.meter.spend(text.len())?;
        let mut start = 0;
        for (index, c) in text.char_indices() {
            if !matches!(c, '\t'|'\n'|'\r'|'\u{20}'..='\u{d7ff}'|'\u{e000}'..='\u{fffd}'|'\u{10000}'..='\u{10ffff}')
            {
                return Err(invalid("invalid XML 1.0 character"));
            }
            let replacement = match c {
                '&' => Some("&amp;"),
                '<' => Some("&lt;"),
                '>' => Some("&gt;"),
                '"' => Some("&quot;"),
                '\'' => Some("&apos;"),
                '\n' => Some("&#10;"),
                '\r' => Some("&#13;"),
                '\t' => Some("&#9;"),
                _ => None,
            };
            if let Some(replacement) = replacement {
                self.raw(&text[start..index])?;
                self.raw(replacement)?;
                start = index + c.len_utf8();
            }
        }
        self.raw(&text[start..])
    }
    fn attribute(&mut self, name: &str, value: &str) -> Result<()> {
        self.raw(" ")?;
        self.raw(name)?;
        self.raw("=\"")?;
        self.escaped(value)?;
        self.raw("\"")
    }
    fn integer(&mut self, value: i64) -> Result<()> {
        self.meter.text(64)?;
        self.raw(&value.to_string())
    }
    fn float(&mut self, value: f64) -> Result<()> {
        self.meter.text(1024)?;
        self.raw(&crate::param::value::format_float(value, true))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embedded_provider_shared_budget() {
        let mut cv = ControlledVocabulary::new();
        let (mut work, mut bytes) = (cv.limits.max_work, cv.limits.max_bytes);
        for ((name, data, encoding), expected) in
            obo::EMBEDDED.iter().zip([3569, 5545, 5780, 9182, 9254])
        {
            cv.load_with_budget(name, *data, *encoding, &mut work, &mut bytes)
                .unwrap();
            eprintln!(
                "{name}: used work={} bytes={}",
                cv.limits.max_work - work,
                cv.limits.max_bytes - bytes
            );
            assert_eq!(cv.terms.len(), expected);
        }
        assert_eq!(cv.names.len(), 16852);
    }

    #[test]
    fn borrowed_queries_share_exhaustion_and_copies_precharge_before_visiting_payload() {
        let mut cv = ControlledVocabulary::new();
        cv.load_obo_reader(
            "x",
            b"[Term]\nid: A\nname: root\n[Term]\nid: B\nname: child\nis_a: A\n".as_slice(),
        )
        .unwrap();
        let (mut work, mut bytes) = (100_000, 100_000);
        cv.get_term_with_budget("A", &mut work, &mut bytes).unwrap();
        let spent = 100_000 - work;
        let mut work = spent;
        cv.get_term_with_budget("A", &mut work, &mut bytes).unwrap();
        assert_eq!(work, 0);
        assert!(
            cv.find_term_by_name_with_budget("root", &mut work, &mut bytes)
                .is_err()
        );
        let (mut work, mut bytes) = (100_000, 100_000);
        assert_eq!(
            cv.first_child_with_name_with_budget("A", "child", &mut work, &mut bytes)
                .unwrap()
                .unwrap()
                .id,
            "B"
        );
        assert!(work < 100_000);
        assert!(bytes < 100_000);
        let (mut work, mut bytes) = (1, usize::MAX);
        assert!(cv.checked_clone_with_budget(&mut work, &mut bytes).is_err());
        let mut values = BTreeSet::from(["sentinel".into()]);
        cv.limits.max_work = 1;
        assert!(cv.extend_child_terms(&mut values, "A").is_err());
        assert_eq!(values, BTreeSet::from(["sentinel".into()]));
    }

    #[test]
    fn sparse_large_term_map_root_is_charged_before_clone() {
        let mut cv = ControlledVocabulary::new();
        cv.load_obo_reader("x", b"[Term]\nid: A\n".as_slice())
            .unwrap();
        cv.limits.max_bytes = 3000;
        assert!(cv.checked_clone().is_err());
        assert_eq!(cv.get_term("A").unwrap().id, "A");
    }
}
