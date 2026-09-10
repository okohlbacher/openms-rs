// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use super::*;

impl ParamNode {
    pub fn new(name: &str, description: &str) -> Result<Self> {
        let mut w = ParamWork::default();
        Ok(Self {
            name: w.text(name)?,
            description: w.text(description)?,
            entries: Vec::new(),
            nodes: Vec::new(),
        })
    }
    pub fn suffix(key: &str) -> &str {
        key.rsplit_once(':').map_or(key, |(_, suffix)| suffix)
    }
    pub fn size(&self) -> Result<usize> {
        self.measure(&mut ParamWork::default())?;
        Ok(self.count_entries())
    }
    pub(super) fn count_entries(&self) -> usize {
        self.entries.len() + self.nodes.iter().map(Self::count_entries).sum::<usize>()
    }
    pub fn find_entry(&self, name: &str) -> Result<Option<&ParamEntry>> {
        let mut w = ParamWork::default();
        Ok(self.entry_index(name, &mut w)?.map(|i| &self.entries[i]))
    }
    pub fn find_node(&self, name: &str) -> Result<Option<&ParamNode>> {
        self.local_node(name, &mut ParamWork::default())
    }
    pub fn find_parent_of(&self, name: &str) -> Result<Option<&ParamNode>> {
        let mut w = ParamWork::default();
        key_check(name, &mut w)?;
        self.parent(name, &mut w)
    }
    pub fn find_entry_recursive(&self, name: &str) -> Result<Option<&ParamEntry>> {
        let mut w = ParamWork::default();
        key_check(name, &mut w)?;
        self.find_entry_recursive_with_work(name, &mut w)
    }
    pub(super) fn entry_index(&self, name: &str, w: &mut ParamWork) -> Result<Option<usize>> {
        for (i, e) in self.entries.iter().enumerate() {
            w.consume(add(1, add(name.len(), e.name.len())?)?)?;
            if e.name == name {
                return Ok(Some(i));
            }
        }
        Ok(None)
    }
    pub(super) fn node_index(&self, name: &str, w: &mut ParamWork) -> Result<Option<usize>> {
        for (i, n) in self.nodes.iter().enumerate() {
            w.consume(add(1, add(name.len(), n.name.len())?)?)?;
            if n.name == name {
                return Ok(Some(i));
            }
        }
        Ok(None)
    }
    pub(super) fn local_node(&self, name: &str, w: &mut ParamWork) -> Result<Option<&Self>> {
        Ok(self.node_index(name, w)?.map(|i| &self.nodes[i]))
    }
    pub(super) fn parent(&self, name: &str, w: &mut ParamWork) -> Result<Option<&Self>> {
        let mut node = self;
        let mut key = name;
        while let Some((first, rest)) = key.split_once(':') {
            w.consume(key.len())?;
            let Some(i) = node.node_index(first, w)? else {
                return Ok(None);
            };
            node = &node.nodes[i];
            key = rest;
        }
        for n in &node.nodes {
            w.consume(add(1, add(n.name.len(), key.len())?)?)?;
            if n.name.starts_with(key) {
                return Ok(Some(node));
            }
        }
        for e in &node.entries {
            w.consume(add(1, add(e.name.len(), key.len())?)?)?;
            if e.name.starts_with(key) {
                return Ok(Some(node));
            }
        }
        Ok(None)
    }
    pub(super) fn parent_mut(
        &mut self,
        name: &str,
        w: &mut ParamWork,
    ) -> Result<Option<&mut Self>> {
        let mut node = self;
        let mut key = name;
        while let Some((first, rest)) = key.split_once(':') {
            w.consume(key.len())?;
            let Some(i) = node.node_index(first, w)? else {
                return Ok(None);
            };
            node = &mut node.nodes[i];
            key = rest;
        }
        let found = node
            .nodes
            .iter()
            .map(|n| &n.name)
            .chain(node.entries.iter().map(|e| &e.name))
            .try_fold(false, |found, n| -> Result<bool> {
                w.consume(add(1, add(n.len(), key.len())?)?)?;
                Ok(found || n.starts_with(key))
            })?;
        Ok(if found { Some(node) } else { None })
    }
    pub(super) fn find_entry_recursive_with_work(
        &self,
        name: &str,
        w: &mut ParamWork,
    ) -> Result<Option<&ParamEntry>> {
        let Some(parent) = self.parent(name, w)? else {
            return Ok(None);
        };
        Ok(parent
            .entry_index(Self::suffix(name), w)?
            .map(|i| &parent.entries[i]))
    }
    fn prepare_path<'a>(&'a mut self, key: &str, w: &mut ParamWork) -> Result<&'a mut Self> {
        key_check(key, w)?;
        let mut node = self;
        let mut rest = key;
        while let Some((first, tail)) = rest.split_once(':') {
            let index = match node.node_index(first, w)? {
                Some(i) => i,
                None => {
                    w.copy(size_of::<Self>())?;
                    let n = Self {
                        name: w.text(first)?,
                        description: String::new(),
                        entries: Vec::new(),
                        nodes: Vec::new(),
                    };
                    node.nodes
                        .try_reserve(1)
                        .map_err(|_| invalid("parameter node allocation failed"))?;
                    node.nodes.push(n);
                    node.nodes.len() - 1
                }
            };
            node = &mut node.nodes[index];
            rest = tail;
        }
        Ok(node)
    }
    pub fn insert_entry(&mut self, entry: &ParamEntry, prefix: &str) -> Result<()> {
        let mut w = ParamWork::default();
        let bytes = self.measure(&mut w)?;
        w.copy(bytes)?;
        entry.measure(&mut w)?;
        let mut copy = self.clone();
        copy.insert_entry_inner(entry, prefix, &mut w)?;
        copy.measure(&mut w)?;
        *self = copy;
        Ok(())
    }
    pub fn insert_node(&mut self, node: &ParamNode, prefix: &str) -> Result<()> {
        let mut w = ParamWork::default();
        let bytes = self.measure(&mut w)?;
        w.copy(bytes)?;
        node.measure(&mut w)?;
        let mut copy = self.clone();
        copy.insert_node_inner(node, prefix, &mut w)?;
        copy.measure(&mut w)?;
        *self = copy;
        Ok(())
    }
    pub(super) fn insert_entry_inner(
        &mut self,
        entry: &ParamEntry,
        prefix: &str,
        w: &mut ParamWork,
    ) -> Result<()> {
        let key = join(prefix, &entry.name, w)?;
        let parent = self.prepare_path(&key, w)?;
        let name = Self::suffix(&key);
        if parent.node_index(name, w)?.is_some() {
            return Err(invalid(format!(
                "parameter entry '{key}' collides with a section"
            )));
        }
        let bytes = entry.measure(w)?;
        w.copy(bytes)?;
        if let Some(i) = parent.entry_index(name, w)? {
            let old = &mut parent.entries[i];
            old.value = entry.value.clone();
            old.tags = entry.tags.clone();
            if old.description.is_empty() || !entry.description.is_empty() {
                old.description = entry.description.clone();
            }
        } else {
            let mut e = entry.clone();
            e.name = w.text(name)?;
            parent
                .entries
                .try_reserve(1)
                .map_err(|_| invalid("parameter entry allocation failed"))?;
            parent.entries.push(e);
        }
        Ok(())
    }
    pub(super) fn insert_node_inner(
        &mut self,
        node: &ParamNode,
        prefix: &str,
        w: &mut ParamWork,
    ) -> Result<()> {
        let key = join(prefix, &node.name, w)?;
        let parent = self.prepare_path(&key, w)?;
        let name = Self::suffix(&key);
        if parent.entry_index(name, w)?.is_some() {
            return Err(invalid(format!(
                "parameter section '{key}' collides with an entry"
            )));
        }
        if let Some(i) = parent.node_index(name, w)? {
            let old = &mut parent.nodes[i];
            for child in &node.nodes {
                old.insert_node_inner(child, "", w)?;
            }
            for entry in &node.entries {
                old.insert_entry_inner(entry, "", w)?;
            }
            if old.description.is_empty() || !node.description.is_empty() {
                old.description = w.text(&node.description)?;
            }
        } else {
            let bytes = node.measure(w)?;
            w.copy(bytes)?;
            let mut n = node.clone();
            n.name = w.text(name)?;
            parent
                .nodes
                .try_reserve(1)
                .map_err(|_| invalid("parameter node allocation failed"))?;
            parent.nodes.push(n);
        }
        Ok(())
    }
    pub fn source_equal(&self, other: &Self) -> Result<bool> {
        let mut w = ParamWork::default();
        self.measure(&mut w)?;
        other.measure(&mut w)?;
        self.source_equal_inner(other, &mut w)
    }
    pub(super) fn source_equal_inner(&self, other: &Self, w: &mut ParamWork) -> Result<bool> {
        w.consume(add(1, add(self.name.len(), other.name.len())?)?)?;
        if self.name != other.name
            || self.entries.len() != other.entries.len()
            || self.nodes.len() != other.nodes.len()
        {
            return Ok(false);
        }
        for entry in &self.entries {
            let mut found = false;
            for candidate in &other.entries {
                let b = add(entry.measure(w)?, candidate.measure(w)?)?;
                w.consume(b)?;
                if entry.name == candidate.name && entry.value == candidate.value {
                    found = true;
                    break;
                }
            }
            if !found {
                return Ok(false);
            }
        }
        for node in &self.nodes {
            let mut found = false;
            for candidate in &other.nodes {
                if node.source_equal_inner(candidate, w)? {
                    found = true;
                    break;
                }
            }
            if !found {
                return Ok(false);
            }
        }
        Ok(true)
    }
    pub(super) fn measure(&self, w: &mut ParamWork) -> Result<usize> {
        let mut stack = vec![(self, 0usize)];
        w.slots::<(&Self, usize)>(MAX_PARAM_DEPTH + 1)?;
        let (mut bytes, mut entries, mut nodes) = (0usize, 0usize, 0usize);
        while let Some((n, depth)) = stack.pop() {
            if depth > MAX_PARAM_DEPTH {
                return Err(invalid("parameter depth limit exceeded"));
            }
            nodes = add(nodes, 1)?;
            entries = add(entries, n.entries.len())?;
            if nodes > MAX_PARAM_NODES || entries > MAX_PARAM_ENTRIES {
                return Err(invalid("parameter tree count limit exceeded"));
            }
            w.consume(add(n.nodes.len(), n.entries.len())?)?;
            bytes = add(
                bytes,
                add(size_of::<Self>(), add(n.name.len(), n.description.len())?)?,
            )?;
            w.consume(add(n.name.len(), n.description.len())?)?;
            for entry in &n.entries {
                bytes = add(bytes, entry.measure(w)?)?;
            }
            if bytes > MAX_PARAM_BYTES {
                return Err(invalid("parameter tree payload limit exceeded"));
            }
            // Bound the traversal stack before copying child references.
            w.slots::<(&Self, usize)>(n.nodes.len())?;
            stack
                .try_reserve(n.nodes.len())
                .map_err(|_| invalid("parameter traversal allocation failed"))?;
            stack.extend(n.nodes.iter().map(|child| (child, depth + 1)));
        }
        Ok(bytes)
    }
}
pub(super) fn join(a: &str, b: &str, w: &mut ParamWork) -> Result<String> {
    let len = add(a.len(), b.len())?;
    w.copy(len)?;
    let mut s = String::new();
    s.try_reserve(len)
        .map_err(|_| invalid("parameter text allocation failed"))?;
    s.push_str(a);
    s.push_str(b);
    Ok(s)
}
pub(super) fn normalized_prefix(prefix: &str, w: &mut ParamWork) -> Result<String> {
    key_check(prefix, w)?;
    let mut p = w.text(prefix)?;
    if !p.is_empty() && !p.ends_with(':') {
        w.copy(1)?;
        p.push(':');
    }
    Ok(p)
}
