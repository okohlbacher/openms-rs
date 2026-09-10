// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ordered hierarchical parameters and source-compatible configuration operations.
//! Public drafts retain independent restrictions; checked operations are atomic.

pub mod handler;
pub use handler::DefaultParamHandler;
mod iteration;
mod operations;
mod tree;
pub mod value;
use crate::{Error, Result};
pub use iteration::{ParamItem, ParamIterator, ParamTrace};
pub use operations::{CommandLineOptions, ParamUpdateOptions, ParamUpdateReport};
use std::collections::BTreeSet;
use std::mem::size_of;
pub use value::{ParamValue, ParamValueType};

pub const MAX_PARAM_WORK: usize = 50_000_000;
pub const MAX_PARAM_BYTES: usize = 64 * 1024 * 1024;
pub const MAX_PARAM_ENTRIES: usize = 100_000;
pub const MAX_PARAM_NODES: usize = 100_000;
pub const MAX_PARAM_DEPTH: usize = 128;

pub(crate) struct ParamWork {
    remaining_work: usize,
    remaining_bytes: usize,
}
impl Default for ParamWork {
    fn default() -> Self {
        Self {
            remaining_work: MAX_PARAM_WORK,
            remaining_bytes: MAX_PARAM_BYTES,
        }
    }
}
impl ParamWork {
    pub(crate) fn consume(&mut self, count: usize) -> Result<()> {
        self.remaining_work = self
            .remaining_work
            .checked_sub(count)
            .ok_or_else(|| invalid("parameter work limit exceeded"))?;
        Ok(())
    }
    pub(crate) fn allocation(&mut self, bytes: usize) -> Result<()> {
        self.remaining_bytes = self
            .remaining_bytes
            .checked_sub(bytes)
            .ok_or_else(|| invalid("parameter allocation limit exceeded"))?;
        Ok(())
    }
    pub(crate) fn copy(&mut self, bytes: usize) -> Result<()> {
        self.consume(bytes)?;
        self.allocation(bytes)
    }
    fn text(&mut self, value: &str) -> Result<String> {
        self.copy(value.len())?;
        Ok(value.to_owned())
    }
    fn slots<T>(&mut self, count: usize) -> Result<()> {
        self.allocation(mul(size_of::<T>(), count)?)
    }
}
fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b)
        .ok_or_else(|| invalid("parameter size overflow"))
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b)
        .ok_or_else(|| invalid("parameter size overflow"))
}
fn missing(key: &str) -> Error {
    invalid(format!("parameter or section '{key}' was not found"))
}
fn key_check(key: &str, work: &mut ParamWork) -> Result<()> {
    work.consume(key.len())?;
    if key.bytes().filter(|b| *b == b':').count() > MAX_PARAM_DEPTH {
        return Err(invalid("parameter depth limit exceeded"));
    }
    Ok(())
}
fn no_comma(value: &str, work: &mut ParamWork) -> Result<()> {
    work.consume(value.len())?;
    if value.contains(',') {
        Err(invalid(
            "commas are not allowed in parameter tags or string restrictions",
        ))
    } else {
        Ok(())
    }
}

/// Complete owned entry. Native equality includes descriptions and restrictions;
/// `source_equal` reproduces the source's name/value-only predicate.
#[derive(Clone, Debug, PartialEq)]
pub struct ParamEntry {
    pub name: String,
    pub description: String,
    pub value: ParamValue,
    pub tags: BTreeSet<String>,
    pub min_float: f64,
    pub max_float: f64,
    pub min_int: i32,
    pub max_int: i32,
    pub valid_strings: Vec<String>,
}
impl Default for ParamEntry {
    fn default() -> Self {
        Self {
            name: String::new(),
            description: String::new(),
            value: ParamValue::Empty,
            tags: BTreeSet::new(),
            min_float: -f64::MAX,
            max_float: f64::MAX,
            min_int: -i32::MAX,
            max_int: i32::MAX,
            valid_strings: Vec::new(),
        }
    }
}
impl ParamEntry {
    pub fn new(name: &str, value: ParamValue, description: &str, tags: &[String]) -> Result<Self> {
        let mut work = ParamWork::default();
        Self::new_with_work(name, value, description, tags, &mut work)
    }
    fn new_with_work(
        name: &str,
        value: ParamValue,
        description: &str,
        tags: &[String],
        work: &mut ParamWork,
    ) -> Result<Self> {
        value.measure(work)?;
        work.consume(tags.len())?;
        work.slots::<String>(tags.len())?;
        let mut set = BTreeSet::new();
        for tag in tags {
            work.copy(add(tag.len(), 128)?)?;
            ordered_lookup_cost(tag, set.len(), work)?;
            set.insert(tag.clone());
        }
        Ok(Self {
            name: work.text(name)?,
            description: work.text(description)?,
            value,
            tags: set,
            ..Self::default()
        })
    }
    pub fn source_equal(&self, other: &Self) -> Result<bool> {
        let mut work = ParamWork::default();
        self.measure(&mut work)?;
        other.measure(&mut work)?;
        Ok(self.name == other.name && self.value == other.value)
    }
    /// Restriction failure is a diagnostic; conversion/resource failure is an error.
    pub fn validation_error(&self) -> Result<Option<String>> {
        self.valid_with_work(&mut ParamWork::default())
    }
    pub fn is_valid(&self) -> Result<bool> {
        Ok(self.validation_error()?.is_none())
    }
    fn valid_with_work(&self, work: &mut ParamWork) -> Result<Option<String>> {
        self.measure(work)?;
        let number_error = || {
            Some(format!(
                "parameter '{}' violates its numeric restrictions",
                self.name
            ))
        };
        let string_ok = |v: &str, work: &mut ParamWork| -> Result<bool> {
            for candidate in &self.valid_strings {
                work.consume(add(1, add(v.len(), candidate.len())?)?)?;
                if candidate == v {
                    return Ok(true);
                }
            }
            Ok(false)
        };
        work.consume(mul(self.tags.len(), 64)?)?;
        let file = self.tags.contains("input file") || self.tags.contains("output file");
        match &self.value {
            ParamValue::String(v)
                if !self.valid_strings.is_empty()
                    && !file
                    && !self.tags.contains("output prefix") =>
            {
                if !string_ok(v, work)? {
                    return Ok(Some(format!(
                        "parameter '{}' has a disallowed string value",
                        self.name
                    )));
                }
            }
            ParamValue::StringList(vs) if !self.valid_strings.is_empty() && !file => {
                for v in vs {
                    if !string_ok(v, work)? {
                        return Ok(Some(format!(
                            "parameter '{}' has a disallowed string-list value",
                            self.name
                        )));
                    }
                }
            }
            ParamValue::Integer(_) => {
                let x = self.value.to_i32()?;
                if (self.min_int != -i32::MAX && x < self.min_int)
                    || (self.max_int != i32::MAX && x > self.max_int)
                {
                    return Ok(number_error());
                }
            }
            ParamValue::IntegerList(xs) => {
                for &x in xs {
                    work.consume(1)?;
                    if (self.min_int != -i32::MAX && x < self.min_int)
                        || (self.max_int != i32::MAX && x > self.max_int)
                    {
                        return Ok(number_error());
                    }
                }
            }
            ParamValue::Float(x) => {
                if (self.min_float != -f64::MAX && *x < self.min_float)
                    || (self.max_float != f64::MAX && *x > self.max_float)
                {
                    return Ok(number_error());
                }
            }
            ParamValue::FloatList(xs) => {
                for &x in xs {
                    work.consume(1)?;
                    if (self.min_float != -f64::MAX && x < self.min_float)
                        || (self.max_float != f64::MAX && x > self.max_float)
                    {
                        return Ok(number_error());
                    }
                }
            }
            _ => {}
        }
        Ok(None)
    }
    fn measure(&self, work: &mut ParamWork) -> Result<usize> {
        work.consume(add(self.tags.len(), self.valid_strings.len())?)?;
        let mut bytes = add(size_of::<Self>(), self.value.measure(work)?)?;
        for text in std::iter::once(&self.name)
            .chain(std::iter::once(&self.description))
            .chain(self.tags.iter())
            .chain(self.valid_strings.iter())
        {
            work.consume(text.len())?;
            bytes = add(bytes, text.len())?;
        }
        bytes = add(bytes, mul(self.tags.len(), 128)?)?;
        bytes = add(bytes, mul(self.valid_strings.len(), size_of::<String>())?)?;
        if bytes > MAX_PARAM_BYTES {
            return Err(invalid("parameter entry payload limit exceeded"));
        }
        Ok(bytes)
    }
}

/// An ordered draft tree. Leaf and child order are preserved independently.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParamNode {
    pub name: String,
    pub description: String,
    pub entries: Vec<ParamEntry>,
    pub nodes: Vec<ParamNode>,
}
// Public drafts may be deeper than the checked limit. Drain children iteratively
// even when rejecting such a draft, avoiding a recursive destructor overflow.
impl Drop for ParamNode {
    fn drop(&mut self) {
        while let Some(mut child) = self.nodes.pop() {
            self.nodes.append(&mut child.nodes);
        }
    }
}

/// Invisible root plus validated bounded tree. Rust equality compares full values
/// and insertion order; use `source_equal` for the historical source predicate.
#[derive(Clone, Debug, PartialEq)]
pub struct Param {
    root: ParamNode,
}
impl Default for Param {
    fn default() -> Self {
        Self::new()
    }
}
impl Param {
    pub fn new() -> Self {
        Self {
            root: ParamNode {
                name: "ROOT".into(),
                description: String::new(),
                entries: Vec::new(),
                nodes: Vec::new(),
            },
        }
    }
    pub fn from_root(mut root: ParamNode) -> Result<Self> {
        root.measure(&mut ParamWork::default())?;
        root.name = "ROOT".into();
        root.description.clear();
        Ok(Self { root })
    }
    pub fn root(&self) -> &ParamNode {
        &self.root
    }
    pub fn checked_clone(&self) -> Result<Self> {
        let mut work = ParamWork::default();
        self.clone_with_work(&mut work)
    }
    fn clone_with_work(&self, work: &mut ParamWork) -> Result<Self> {
        let bytes = self.root.measure(work)?;
        work.copy(bytes)?;
        Ok(self.clone())
    }
    fn edit<T>(
        &mut self,
        action: impl FnOnce(&mut Param, &mut ParamWork) -> Result<T>,
    ) -> Result<T> {
        let mut work = ParamWork::default();
        // ponytail: one bounded tree snapshot per atomic operation; batch callers
        // can construct an owned ParamNode and validate once with from_root.
        let mut staged = self.clone_with_work(&mut work)?;
        let result = action(&mut staged, &mut work)?;
        staged.root.measure(&mut work)?;
        *self = staged;
        Ok(result)
    }
    pub fn size(&self) -> usize {
        self.root.count_entries()
    }
    pub fn is_empty(&self) -> bool {
        self.size() == 0
    }
    pub fn clear(&mut self) {
        *self = Self::new();
    }
    pub fn source_equal(&self, other: &Self) -> Result<bool> {
        self.root.source_equal(&other.root)
    }
    pub fn entry(&self, key: &str) -> Result<&ParamEntry> {
        let mut work = ParamWork::default();
        key_check(key, &mut work)?;
        self.root
            .find_entry_recursive_with_work(key, &mut work)?
            .ok_or_else(|| missing(key))
    }
    pub fn value(&self, key: &str) -> Result<&ParamValue> {
        Ok(&self.entry(key)?.value)
    }
    pub fn value_type(&self, key: &str) -> Result<ParamValueType> {
        Ok(self.value(key)?.value_type())
    }
    pub fn description(&self, key: &str) -> Result<&str> {
        Ok(&self.entry(key)?.description)
    }
    pub fn exists(&self, key: &str) -> Result<bool> {
        let mut work = ParamWork::default();
        key_check(key, &mut work)?;
        Ok(self
            .root
            .find_entry_recursive_with_work(key, &mut work)?
            .is_some())
    }
    /// Preserves the source's prefix-parent test, including leaf/prefix matches.
    /// The source's undefined empty-key case is a checked error.
    pub fn has_section(&self, key: &str) -> Result<bool> {
        if key.is_empty() {
            return Err(invalid("empty section query"));
        }
        self.root
            .find_parent_of(key.strip_suffix(':').unwrap_or(key))
            .map(|x| x.is_some())
    }
    pub fn section_description(&self, key: &str) -> Result<&str> {
        let mut work = ParamWork::default();
        self.section_description_with_work(key, &mut work)
    }
    fn section_description_with_work(&self, key: &str, work: &mut ParamWork) -> Result<&str> {
        key_check(key, work)?;
        let Some(node) = self.root.parent(key, work)? else {
            return Ok("");
        };
        Ok(node
            .local_node(ParamNode::suffix(key), work)?
            .map_or("", |n| n.description.as_str()))
    }
    pub fn set_value(
        &mut self,
        key: &str,
        value: ParamValue,
        description: &str,
        tags: &[String],
    ) -> Result<()> {
        self.edit(|p, w| {
            let e = ParamEntry::new_with_work("", value, description, tags, w)?;
            p.root.insert_entry_inner(&e, key, w)
        })
    }
    pub fn insert_entry(&mut self, entry: ParamEntry, prefix: &str) -> Result<()> {
        self.edit(|p, w| {
            entry.measure(w)?;
            p.root.insert_entry_inner(&entry, prefix, w)
        })
    }
    pub fn add_section(&mut self, key: &str, description: &str) -> Result<()> {
        self.edit(|p, w| {
            let n = ParamNode {
                name: String::new(),
                description: w.text(description)?,
                entries: Vec::new(),
                nodes: Vec::new(),
            };
            p.root.insert_node_inner(&n, key, w)
        })
    }
    pub fn set_section_description(&mut self, key: &str, description: &str) -> Result<()> {
        self.edit(|p, w| p.set_section_description_inner(key, description, w))
    }
    fn set_section_description_inner(
        &mut self,
        key: &str,
        description: &str,
        w: &mut ParamWork,
    ) -> Result<()> {
        key_check(key, w)?;
        let text = w.text(description)?;
        let parent = self.root.parent_mut(key, w)?.ok_or_else(|| missing(key))?;
        let i = parent
            .node_index(ParamNode::suffix(key), w)?
            .ok_or_else(|| missing(key))?;
        parent.nodes[i].description = text;
        Ok(())
    }
    pub fn tags(&self, key: &str) -> Result<&BTreeSet<String>> {
        Ok(&self.entry(key)?.tags)
    }
    pub fn has_tag(&self, key: &str, tag: &str) -> Result<bool> {
        let mut w = ParamWork::default();
        key_check(key, &mut w)?;
        let entry = self
            .root
            .find_entry_recursive_with_work(key, &mut w)?
            .ok_or_else(|| missing(key))?;
        ordered_lookup_cost(tag, entry.tags.len(), &mut w)?;
        Ok(entry.tags.contains(tag))
    }
    pub fn add_tag(&mut self, key: &str, tag: &str) -> Result<()> {
        let mut work = ParamWork::default();
        no_comma(tag, &mut work)?;
        self.add_tags(key, &[tag.to_owned()])
    }
    pub fn add_tags(&mut self, key: &str, tags: &[String]) -> Result<()> {
        self.edit(|p, w| {
            let entry = p.entry_mut(key, w)?;
            w.consume(tags.len())?;
            for tag in tags {
                no_comma(tag, w)?;
                w.copy(add(tag.len(), 128)?)?;
                ordered_lookup_cost(tag, entry.tags.len(), w)?;
                entry.tags.insert(tag.clone());
            }
            Ok(())
        })
    }
    pub fn clear_tags(&mut self, key: &str) -> Result<()> {
        self.edit(|p, w| {
            p.entry_mut(key, w)?.tags.clear();
            Ok(())
        })
    }
    fn entry_mut(&mut self, key: &str, w: &mut ParamWork) -> Result<&mut ParamEntry> {
        key_check(key, w)?;
        let parent = self.root.parent_mut(key, w)?.ok_or_else(|| missing(key))?;
        let i = parent
            .entry_index(ParamNode::suffix(key), w)?
            .ok_or_else(|| missing(key))?;
        Ok(&mut parent.entries[i])
    }
    pub fn set_valid_strings(&mut self, key: &str, strings: &[String]) -> Result<()> {
        self.edit(|p, w| {
            let e = p.entry_mut(key, w)?;
            check_type(&e.value, 0)?;
            w.slots::<String>(strings.len())?;
            for s in strings {
                no_comma(s, w)?;
                w.copy(s.len())?;
            }
            e.valid_strings = strings.to_vec();
            Ok(())
        })
    }
    pub fn valid_strings(&self, key: &str) -> Result<&[String]> {
        let e = self.entry(key)?;
        check_type(&e.value, 0)?;
        Ok(&e.valid_strings)
    }
    pub fn set_min_int(&mut self, key: &str, value: i32) -> Result<()> {
        self.edit(|p, w| {
            let e = p.entry_mut(key, w)?;
            check_type(&e.value, 1)?;
            e.min_int = value;
            Ok(())
        })
    }
    pub fn set_max_int(&mut self, key: &str, value: i32) -> Result<()> {
        self.edit(|p, w| {
            let e = p.entry_mut(key, w)?;
            check_type(&e.value, 1)?;
            e.max_int = value;
            Ok(())
        })
    }
    pub fn set_min_float(&mut self, key: &str, value: f64) -> Result<()> {
        self.edit(|p, w| {
            let e = p.entry_mut(key, w)?;
            check_type(&e.value, 2)?;
            e.min_float = value;
            Ok(())
        })
    }
    pub fn set_max_float(&mut self, key: &str, value: f64) -> Result<()> {
        self.edit(|p, w| {
            let e = p.entry_mut(key, w)?;
            check_type(&e.value, 2)?;
            e.max_float = value;
            Ok(())
        })
    }
}
fn check_type(v: &ParamValue, group: u8) -> Result<()> {
    let ok = match group {
        0 => matches!(v, ParamValue::String(_) | ParamValue::StringList(_)),
        1 => matches!(v, ParamValue::Integer(_) | ParamValue::IntegerList(_)),
        _ => matches!(v, ParamValue::Float(_) | ParamValue::FloatList(_)),
    };
    if ok {
        Ok(())
    } else {
        Err(invalid("restriction has incompatible parameter type"))
    }
}

/// Owned draft used by the XML reader. All edits share one budget; any failed
/// edit poisons the builder, so partial changes can never be published.
#[cfg(feature = "paramxml")]
pub(crate) struct ParamBuilder {
    param: Param,
    work: ParamWork,
    failed: bool,
}
#[cfg(feature = "paramxml")]
impl std::ops::Deref for ParamBuilder {
    type Target = Param;
    fn deref(&self) -> &Param {
        &self.param
    }
}
#[cfg(feature = "paramxml")]
impl ParamBuilder {
    pub(crate) fn new(param: Param) -> Result<Self> {
        let mut work = ParamWork::default();
        param.root.measure(&mut work)?;
        Ok(Self {
            param,
            work,
            failed: false,
        })
    }
    pub(crate) fn finish(mut self) -> Result<Param> {
        if self.failed {
            return Err(invalid("parameter builder previously failed"));
        }
        self.param.root.measure(&mut self.work)?;
        Ok(self.param)
    }
    fn edit(
        &mut self,
        action: impl FnOnce(&mut Param, &mut ParamWork) -> Result<()>,
    ) -> Result<()> {
        if self.failed {
            return Err(invalid("parameter builder previously failed"));
        }
        let result = action(&mut self.param, &mut self.work);
        if result.is_err() {
            self.failed = true;
        }
        result
    }
    pub(crate) fn set_value(
        &mut self,
        key: &str,
        value: ParamValue,
        description: &str,
        tags: &[String],
    ) -> Result<()> {
        self.edit(|p, w| {
            let e = ParamEntry::new_with_work("", value, description, tags, w)?;
            p.root.insert_entry_inner(&e, key, w)
        })
    }
    pub(crate) fn add_section(&mut self, key: &str, description: &str) -> Result<()> {
        self.edit(|p, w| {
            let n = ParamNode {
                name: String::new(),
                description: w.text(description)?,
                entries: Vec::new(),
                nodes: Vec::new(),
            };
            p.root.insert_node_inner(&n, key, w)
        })
    }
    pub(crate) fn set_valid_strings(&mut self, key: &str, strings: &[String]) -> Result<()> {
        self.edit(|p, w| {
            let e = p.entry_mut(key, w)?;
            check_type(&e.value, 0)?;
            w.slots::<String>(strings.len())?;
            for s in strings {
                no_comma(s, w)?;
                w.copy(s.len())?;
            }
            e.valid_strings = strings.to_vec();
            Ok(())
        })
    }
    pub(crate) fn set_min_int(&mut self, key: &str, value: i32) -> Result<()> {
        self.edit(|p, w| {
            let e = p.entry_mut(key, w)?;
            check_type(&e.value, 1)?;
            e.min_int = value;
            Ok(())
        })
    }
    pub(crate) fn set_max_int(&mut self, key: &str, value: i32) -> Result<()> {
        self.edit(|p, w| {
            let e = p.entry_mut(key, w)?;
            check_type(&e.value, 1)?;
            e.max_int = value;
            Ok(())
        })
    }
    pub(crate) fn set_min_float(&mut self, key: &str, value: f64) -> Result<()> {
        self.edit(|p, w| {
            let e = p.entry_mut(key, w)?;
            check_type(&e.value, 2)?;
            e.min_float = value;
            Ok(())
        })
    }
    pub(crate) fn set_max_float(&mut self, key: &str, value: f64) -> Result<()> {
        self.edit(|p, w| {
            let e = p.entry_mut(key, w)?;
            check_type(&e.value, 2)?;
            e.max_float = value;
            Ok(())
        })
    }
}

fn ordered_lookup_cost(key: &str, count: usize, w: &mut ParamWork) -> Result<()> {
    let comparisons = mul(12, add(count.saturating_add(1).ilog2() as usize, 1)?)?;
    w.consume(mul(add(key.len(), 1)?, comparisons)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn empty_name_comparisons_still_consume_work() {
        let mut node = ParamNode::default();
        node.entries = vec![ParamEntry::default(); 4];
        let mut work = ParamWork {
            remaining_work: 0,
            remaining_bytes: 1024,
        };
        assert!(node.entry_index("", &mut work).is_err());
    }
    #[test]
    fn operation_copy_failure_is_atomic() {
        let mut p = Param::new();
        p.set_value("a", 1.into(), "", &[]).unwrap();
        let before = p.clone();
        let mut work = ParamWork {
            remaining_work: MAX_PARAM_WORK,
            remaining_bytes: 1,
        };
        assert!(p.clone_with_work(&mut work).is_err());
        assert_eq!(p, before);
    }
    #[cfg(feature = "paramxml")]
    #[test]
    fn owned_builder_shares_work_and_cannot_publish_partial_failure() {
        let mut builder = ParamBuilder::new(Param::new()).unwrap();
        builder.set_value("a", 1.into(), "", &[]).unwrap();
        builder.work.remaining_work = 0;
        assert!(builder.set_value("b", 2.into(), "", &[]).is_err());
        assert!(builder.finish().is_err());
    }
}
