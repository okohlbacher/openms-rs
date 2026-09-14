// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Source-backed experimental header registries and checked metadata codecs.
//! The header tree is bounded and discarded once its references are resolved.

use super::*;
use crate::{
    data_structures::{DateTime, list::ListParse},
    format::controlled_vocabulary::{ControlledVocabulary, XRefType},
    kernel::data_array::Meter,
    metadata::*,
};
use std::{mem::size_of, sync::Arc};

#[path = "mzml_header/instrument.rs"]
mod instrument;
#[path = "mzml_header/instrument_terms.rs"]
mod instrument_terms;
#[path = "mzml_header/mapping.rs"]
mod mapping;
#[path = "mzml_header/read.rs"]
mod read;
#[path = "mzml_header/xml.rs"]
mod xml;
pub(super) use read::Registry;
#[path = "mzml_header/write.rs"]
mod write;
pub(super) use write::{ArrayHeader, Plan, guard, prepare};

pub(super) const MAX_WORK: usize = 50_000_000;
pub(super) const MAX_BYTES: usize = 256 * 1024 * 1024;

pub(super) struct Work {
    pub remaining: usize,
    pub bytes: usize,
}
impl Default for Work {
    fn default() -> Self {
        Self {
            remaining: MAX_WORK,
            bytes: MAX_BYTES,
        }
    }
}
impl Work {
    /// Spend `count` work units and `bytes` allocation bytes.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Parse`] when either allowance would go below zero. Work
    /// is spent before bytes, so when only the byte allowance is short the work
    /// units stay spent; every caller abandons the operation on this error.
    pub fn charge(&mut self, count: usize, bytes: usize) -> Result<()> {
        self.remaining = self.remaining.checked_sub(count).ok_or_else(resource)?;
        self.bytes = self.bytes.checked_sub(bytes).ok_or_else(resource)?;
        Ok(())
    }
    /// An owned copy of `text`, charged by its length in work and bytes before
    /// it is allocated.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Parse`] when the allowance is exhausted.
    pub fn copy(&mut self, text: &str) -> Result<String> {
        self.charge(text.len(), text.len())?;
        Ok(text.into())
    }
    /// Charge `count` inline slots of `T`: `count` work units and
    /// `count * size_of::<T>()` bytes.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Parse`] on overflow or when the allowance is exhausted.
    pub fn slots<T>(&mut self, count: usize) -> Result<()> {
        self.charge(
            count,
            count.checked_mul(size_of::<T>()).ok_or_else(resource)?,
        )
    }
    /// The shared data-array meter over this allowance, for tree, text and CV
    /// term charges.
    pub fn meter(&mut self) -> Meter<'_> {
        Meter {
            work: &mut self.remaining,
            bytes: &mut self.bytes,
        }
    }
    fn cv(
        &mut self,
        id: &str,
    ) -> Result<&'static crate::format::controlled_vocabulary::CVTermDefinition> {
        ControlledVocabulary::psi_ms()?.get_term_with_budget(
            id,
            &mut self.remaining,
            &mut self.bytes,
        )
    }
    fn child(&mut self, id: &str, parent: &str) -> Result<bool> {
        ControlledVocabulary::psi_ms()?.is_child_of_with_budget(
            id,
            parent,
            &mut self.remaining,
            &mut self.bytes,
        )
    }
}
fn resource() -> Error {
    invalid("mzML header work/allocation allowance exceeded")
}

#[derive(Default)]
pub(super) struct Node {
    name: String,
    attrs: BTreeMap<String, String>,
    children: Vec<Node>,
}
impl Node {
    fn get(&self, name: &str) -> Result<&str> {
        required(&self.attrs, name)
    }
    fn optional(&self, name: &str) -> &str {
        self.attrs.get(name).map_or("", String::as_str)
    }
    fn id(&self) -> Result<&str> {
        parameter_id(self.get("id")?)
    }
    /// Header list children. The schema's `count` attribute must be present and
    /// numeric, but a value disagreeing with the actual number of children is
    /// **advisory on reading**: the source reader ignores it and real files
    /// carry wrong counts. The OpenMS TOPP fixture `DTAExtractor_1_input.mzML`
    /// declares `softwareList count="5"` with four entries and
    /// `dataProcessingList count="3"` with one, and C++ loads it. Rejecting the
    /// mismatch made the port unable to read its own reference data. Writing
    /// still emits the true count.
    fn children(&self, expected: &str) -> Result<&[Node]> {
        if self.children.iter().any(|child| child.name != expected) {
            return Err(invalid(format!("unexpected child of {}", self.name)));
        }
        let _declared: usize = number(self.get("count")?, "header list count")?;
        Ok(&self.children)
    }
}

#[derive(Default)]
pub(super) struct Draft {
    roots: Vec<Node>,
    stack: Vec<Node>,
}
impl Draft {
    /// Whether the element `tag` under `parent` belongs to the retained header:
    /// one of the five header lists directly under `mzML`, or any descendant of
    /// an element already being captured.
    pub fn captures(&self, tag: &str, parent: &str) -> bool {
        !self.stack.is_empty()
            || matches!(
                (parent, tag),
                (
                    "mzML",
                    "fileDescription"
                        | "sampleList"
                        | "softwareList"
                        | "instrumentConfigurationList"
                        | "dataProcessingList"
                )
            )
    }
    /// Open a captured element with its attributes, charging the node, its
    /// name and its attribute payload before it is retained.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Parse`] beyond 16 nesting levels and when the
    /// allowance is exhausted.
    pub fn start(
        &mut self,
        tag: &str,
        attrs: BTreeMap<String, String>,
        work: &mut Work,
    ) -> Result<()> {
        if self.stack.len() >= 16 {
            return Err(invalid("mzML header nesting exceeds 16 levels"));
        }
        // Each node can cause minimum-four growth of its parent and the
        // construction stack; twelve slots cover both cumulative growth sums.
        work.slots::<Node>(12)?;
        work.charge(tag.len(), tag.len())?;
        work.meter().tree::<(String, String)>(attrs.len())?;
        for (key, value) in &attrs {
            work.charge(key.len().checked_add(value.len()).ok_or_else(resource)?, 0)?;
        }
        self.stack.push(Node {
            name: tag.into(),
            attrs,
            children: Vec::new(),
        });
        Ok(())
    }
    /// Close a captured element and attach it to its parent or the roots.
    /// Returns `false` when no element is being captured, so the caller handles
    /// the closing tag itself.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Parse`] when `tag` does not match the open element.
    pub fn end(&mut self, tag: &str) -> Result<bool> {
        if self.stack.is_empty() {
            return Ok(false);
        }
        let node = self.stack.pop().unwrap();
        if node.name != tag {
            return Err(invalid("unmatched header closing tag"));
        }
        if let Some(parent) = self.stack.last_mut() {
            parent.children.push(node);
        } else {
            self.roots.push(node);
        }
        Ok(true)
    }
    /// Resolve the captured header when the `run` element opens, with the
    /// `run` attributes in `attrs`, and return the registry for later record
    /// references.
    ///
    /// `source_dangling_references` is `ReadOptions::source_dangling_references`:
    /// `true` substitutes the source's empty values for a `softwareRef` or
    /// data-processing reference that names no definition, `false` rejects it.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Parse`] for an unclosed header element, malformed or
    /// unresolved header content, and exhausted parameter or work allowances.
    pub fn finish(
        self,
        attrs: &BTreeMap<String, String>,
        groups: &BTreeMap<String, Vec<Parameter>>,
        parameters: &mut ParameterBudget,
        work: &mut Work,
        settings: &mut ExperimentalSettings,
        source_dangling_references: bool,
    ) -> Result<Registry> {
        if !self.stack.is_empty() {
            return Err(invalid("unfinished mzML header"));
        }
        read::parse(
            self.roots,
            attrs,
            groups,
            parameters,
            work,
            settings,
            source_dangling_references,
        )
    }
}

struct Context<'a> {
    groups: &'a BTreeMap<String, Vec<Parameter>>,
    parameters: &'a mut ParameterBudget,
    work: &'a mut Work,
}
impl Context<'_> {
    fn params(
        &mut self,
        node: &Node,
        mut callback: impl FnMut(&mut Self, &str, &BTreeMap<String, String>) -> Result<()>,
    ) -> Result<()> {
        let mut users = false;
        for child in &node.children {
            match child.name.as_str() {
                "cvParam" | "userParam" => {
                    if users && child.name == "cvParam" {
                        return Err(invalid("header cvParam follows userParam"));
                    }
                    users |= child.name == "userParam";
                    self.work.charge(256, 0)?;
                    callback(self, &child.name, &child.attrs)?;
                }
                "referenceableParamGroupRef" => {
                    if users {
                        return Err(invalid("header parameter reference follows userParam"));
                    }
                    let id = parameter_id(child.get("ref")?)?;
                    self.work.charge(id.len().saturating_mul(64), 0)?;
                    let parameters = self
                        .groups
                        .get(id)
                        .ok_or_else(|| invalid("unresolved header parameter-group reference"))?;
                    for p in parameters {
                        self.parameters.charge(&p.attrs)?;
                        self.work.charge(256, 0)?;
                        callback(self, p.tag, &p.attrs)?;
                    }
                }
                _ => {}
            }
        }
        Ok(())
    }
    fn value(&mut self, kind: &str, attrs: &BTreeMap<String, String>) -> Result<MetaValue> {
        self.work.meter().tree::<(String, MetaValue)>(1)?;
        // Input attributes were charged before retention; this covers decoded
        // scalar/unit storage and conversion scratch before any owned copy.
        for (key, value) in attrs {
            self.work.charge(
                key.len().saturating_add(value.len()).saturating_mul(4),
                value.len().saturating_mul(4),
            )?;
        }
        if kind == "userParam" {
            return product_user_value(attrs);
        }
        let id = required(attrs, "accession")?;
        let term = self.work.cv(id)?;
        let text = attrs.get("value").map_or("", String::as_str);
        let data = match term.xref_type {
            XRefType::Integer
            | XRefType::NegativeInteger
            | XRefType::PositiveInteger
            | XRefType::NonNegativeInteger
            | XRefType::NonPositiveInteger
                if !text.is_empty() =>
            {
                MetaValueData::Integer(i64::from(i32::from_list_item(text)?))
            }
            XRefType::Decimal if !text.is_empty() => {
                MetaValueData::Float(f64::from_list_item(text)?)
            }
            XRefType::Date if !text.is_empty() => {
                DateTime::parse(text)?;
                MetaValueData::String(text.into())
            }
            XRefType::Boolean if !text.is_empty() => {
                let value = text.to_ascii_lowercase();
                if value != "true" && value != "false" {
                    return Err(invalid("invalid header boolean CV value"));
                }
                MetaValueData::String(value)
            }
            XRefType::None | XRefType::String | XRefType::AnyUri => {
                MetaValueData::String(text.into())
            }
            _ => return Err(invalid("missing numeric/date/boolean header CV value")),
        };
        let value = MetaValue::new(data)?;
        if attrs.contains_key("unitAccession") {
            // The existing scalar decoder preserves caller-supplied unit name
            // and identity; only its unit is reused for the CV-typed scalar.
            let unit = product_user_value(attrs)?
                .unit()
                .cloned()
                .ok_or_else(|| invalid("missing header unit"))?;
            value.with_unit(unit)
        } else if attrs.contains_key("unitName") || attrs.contains_key("unitCvRef") {
            Err(invalid("header unit attributes without unitAccession"))
        } else {
            Ok(value)
        }
    }
    fn meta(&mut self, meta: &mut MetaInfo, name: &str, value: MetaValue) -> Result<()> {
        self.work.meter().tree::<(String, MetaValue)>(1)?;
        self.work
            .charge(name.len().saturating_mul(64), name.len())?;
        if meta.insert(name.into(), value).is_some() {
            return Err(invalid("duplicate header metadata key"));
        }
        Ok(())
    }
}
