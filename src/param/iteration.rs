// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use super::tree::join;
use super::*;

/// Section transition immediately before an entry, or while reaching the end.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParamTrace {
    pub name: String,
    pub description: String,
    pub opened: bool,
}
#[derive(Clone, Debug)]
pub struct ParamItem<'a> {
    pub key: String,
    pub entry: &'a ParamEntry,
    pub trace: Vec<ParamTrace>,
}
/// Bounded, borrowed forward traversal. Exhaustion is fused; no stale entry can
/// be dereferenced. `end_trace` retains the source's final section-close events.
#[derive(Clone, Debug)]
pub struct ParamIterator<'a> {
    items: std::vec::IntoIter<ParamItem<'a>>,
    end_trace: Vec<ParamTrace>,
}
impl Default for ParamIterator<'_> {
    fn default() -> Self {
        Self {
            items: Vec::new().into_iter(),
            end_trace: Vec::new(),
        }
    }
}
impl<'a> Iterator for ParamIterator<'a> {
    type Item = ParamItem<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        self.items.next()
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.items.size_hint()
    }
}
impl ExactSizeIterator for ParamIterator<'_> {}
impl std::iter::FusedIterator for ParamIterator<'_> {}
impl<'a> ParamIterator<'a> {
    pub fn new(root: &'a ParamNode) -> Result<Self> {
        Self::with_work(root, &mut ParamWork::default())
    }
    pub fn end_trace(&self) -> &[ParamTrace] {
        &self.end_trace
    }
    pub(super) fn with_work(root: &'a ParamNode, w: &mut ParamWork) -> Result<Self> {
        root.measure(w)?;
        let mut items = Vec::new();
        let mut trace = Vec::new();
        traverse(root, "", &mut trace, &mut items, w)?;
        Ok(Self {
            items: items.into_iter(),
            end_trace: trace,
        })
    }
}
fn transition(
    node: &ParamNode,
    opened: bool,
    trace: &mut Vec<ParamTrace>,
    w: &mut ParamWork,
) -> Result<()> {
    if !opened
        && trace
            .last()
            .is_some_and(|last| last.opened && last.name == node.name)
    {
        trace.pop();
        return Ok(());
    }
    w.slots::<ParamTrace>(1)?;
    trace.push(ParamTrace {
        name: w.text(&node.name)?,
        description: w.text(&node.description)?,
        opened,
    });
    Ok(())
}
fn traverse<'a>(
    node: &'a ParamNode,
    prefix: &str,
    trace: &mut Vec<ParamTrace>,
    items: &mut Vec<ParamItem<'a>>,
    w: &mut ParamWork,
) -> Result<()> {
    for e in &node.entries {
        w.consume(1)?;
        w.slots::<ParamItem<'a>>(1)?;
        items.push(ParamItem {
            key: join(prefix, &e.name, w)?,
            entry: e,
            trace: std::mem::take(trace),
        });
    }
    for child in &node.nodes {
        w.consume(1)?;
        transition(child, true, trace, w)?;
        let path = join(&join(prefix, &child.name, w)?, ":", w)?;
        traverse(child, &path, trace, items, w)?;
        transition(child, false, trace, w)?;
    }
    Ok(())
}
impl Param {
    pub fn iter(&self) -> Result<ParamIterator<'_>> {
        ParamIterator::new(&self.root)
    }
    /// The source searches for `:leaf`, excluding root-level entries.
    pub fn find_first(&self, leaf: &str) -> Result<Option<ParamItem<'_>>> {
        self.find_after(leaf, None)
    }
    /// Continue after a borrowed entry. A foreign/stale entry is rejected.
    pub fn find_next(&self, leaf: &str, start: &ParamEntry) -> Result<Option<ParamItem<'_>>> {
        self.find_after(leaf, Some(start))
    }
    fn find_after(&self, leaf: &str, start: Option<&ParamEntry>) -> Result<Option<ParamItem<'_>>> {
        let mut w = ParamWork::default();
        let suffix = join(":", leaf, &mut w)?;
        let mut reached = start.is_none();
        for item in ParamIterator::with_work(&self.root, &mut w)? {
            if !reached {
                if start.is_some_and(|s| std::ptr::eq(s, item.entry)) {
                    reached = true;
                }
                continue;
            }
            w.consume(item.key.len())?;
            if item.key.ends_with(&suffix) {
                return Ok(Some(item));
            }
        }
        if !reached {
            return Err(invalid("parameter search cursor belongs to another tree"));
        }
        Ok(None)
    }
    pub fn to_text(&self) -> Result<String> {
        let mut w = ParamWork::default();
        let mut result = String::new();
        for item in ParamIterator::with_work(&self.root, &mut w)? {
            let value = item.entry.value.to_stream_text_with_work(&mut w)?;
            w.copy(value.len())?;
            let prefix = if item.key.len() > item.entry.name.len() + 1 {
                Some(&item.key[..item.key.len() - item.entry.name.len() - 1])
            } else {
                None
            };
            let line_bytes = add(
                12,
                add(
                    item.key.len(),
                    add(value.len(), item.entry.description.len())?,
                )?,
            )?;
            w.copy(line_bytes)?;
            result
                .try_reserve(line_bytes)
                .map_err(|_| invalid("parameter rendering allocation failed"))?;
            result.push('"');
            if let Some(prefix) = prefix {
                result.push_str(prefix);
                result.push('|');
            }
            result.push_str(&item.entry.name);
            result.push_str("\" -> \"");
            result.push_str(&value);
            result.push('"');

            if !item.entry.description.is_empty() {
                result.push_str(" (");
                result.push_str(&item.entry.description);
                result.push(')');
            }
            result.push('\n');
        }
        Ok(result)
    }
}
impl std::fmt::Display for Param {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.to_text().map_err(|_| std::fmt::Error)?)
    }
}
