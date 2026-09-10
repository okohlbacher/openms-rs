// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use super::tree::{join, normalized_prefix};
use super::*;
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ParamUpdateOptions {
    pub verbose: bool,
    pub add_unknown: bool,
    pub fail_on_invalid_values: bool,
    pub fail_on_unknown_parameters: bool,
}
impl Default for ParamUpdateOptions {
    fn default() -> Self {
        Self {
            verbose: true,
            add_unknown: false,
            fail_on_invalid_values: false,
            fail_on_unknown_parameters: false,
        }
    }
}
/// A false success flag retains successful updates, as in the source. Errors
/// (including resource failures) instead leave the entire destination unchanged.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParamUpdateReport {
    pub success: bool,
    pub messages: Vec<String>,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CommandLineOptions {
    pub one_argument: BTreeMap<String, String>,
    pub no_argument: BTreeMap<String, String>,
    pub multiple_arguments: BTreeMap<String, String>,
    pub misc: String,
    pub unknown: String,
}
impl Default for CommandLineOptions {
    fn default() -> Self {
        Self {
            one_argument: BTreeMap::new(),
            no_argument: BTreeMap::new(),
            multiple_arguments: BTreeMap::new(),
            misc: "misc".into(),
            unknown: "unknown".into(),
        }
    }
}
fn message(messages: &mut Vec<String>, text: String, w: &mut ParamWork) -> Result<()> {
    w.copy(add(text.len(), size_of::<String>())?)?;
    messages.push(text);
    Ok(())
}
fn ancestor(key: &str) -> &str {
    key.rfind(':').map_or("", |i| &key[..i + 1])
}
fn trace_path(
    path: &mut String,
    event: &ParamTrace,
    check_suffix: bool,
    w: &mut ParamWork,
) -> Result<()> {
    if event.opened {
        *path = join(&join(path, &event.name, w)?, ":", w)?;
    } else {
        let suffix = join(&event.name, ":", w)?;
        if !check_suffix || path.ends_with(&suffix) {
            let len = path
                .len()
                .checked_sub(suffix.len())
                .ok_or_else(|| invalid("invalid parameter section trace"))?;
            if !path.is_char_boundary(len) {
                return Err(invalid("invalid parameter section trace boundary"));
            }
            path.truncate(len);
        }
    }
    Ok(())
}
impl Param {
    pub fn insert(&mut self, prefix: &str, other: &Self) -> Result<()> {
        self.edit(|p, w| {
            key_check(prefix, w)?;
            other.root.measure(w)?;
            for node in &other.root.nodes {
                p.root.insert_node_inner(node, prefix, w)?;
            }
            for entry in &other.root.entries {
                p.root.insert_entry_inner(entry, prefix, w)?;
            }
            Ok(())
        })
    }
    pub fn remove(&mut self, key: &str) -> Result<()> {
        self.edit(|p, w| p.remove_inner(key, false, w))
    }
    pub fn remove_all(&mut self, prefix: &str) -> Result<()> {
        self.edit(|p, w| p.remove_all_with_work(prefix, w))
    }
    pub(super) fn remove_all_with_work(&mut self, prefix: &str, w: &mut ParamWork) -> Result<()> {
        self.remove_inner(prefix, true, w)
    }
    fn remove_inner(&mut self, key: &str, all: bool, w: &mut ParamWork) -> Result<()> {
        key_check(key, w)?;
        let mut current = w.text(key)?;
        loop {
            let section = current.ends_with(':');
            let search = if section {
                &current[..current.len() - 1]
            } else {
                &current
            };
            let Some(parent) = self.root.parent_mut(search, w)? else {
                return Ok(());
            };
            let suffix = ParamNode::suffix(search);
            let mut removed = false;
            if section {
                if let Some(i) = parent.node_index(suffix, w)? {
                    w.consume(mul(parent.nodes.len() - i - 1, size_of::<ParamNode>())?)?;
                    parent.nodes.remove(i);
                    removed = true;
                }
            } else if all {
                removed |= remove_matching(&mut parent.nodes, suffix, |n| n.name.as_str(), w)?;
                removed |= remove_matching(&mut parent.entries, suffix, |e| e.name.as_str(), w)?;
            } else if let Some(i) = parent.entry_index(suffix, w)? {
                w.consume(mul(parent.entries.len() - i - 1, size_of::<ParamEntry>())?)?;
                parent.entries.remove(i);
                removed = true;
            }
            // Source removeAll prunes an empty matched parent even if the last
            // prefix has no exact leaf; ordinary remove requires an erased item.
            if (!removed && !all) || !parent.entries.is_empty() || !parent.nodes.is_empty() {
                return Ok(());
            }
            let next = &search[..search.len() - suffix.len()];
            if next == current {
                return Ok(());
            } // native termination for empty-root case
            current = w.text(next)?;
        }
    }
    pub fn copy(&self, prefix: &str, remove_prefix: bool) -> Result<Self> {
        self.copy_with_work(prefix, remove_prefix, &mut ParamWork::default())
    }
    fn copy_with_work(&self, prefix: &str, remove_prefix: bool, w: &mut ParamWork) -> Result<Self> {
        key_check(prefix, w)?;
        self.root.measure(w)?;
        let mut out = Self::new();
        let Some(node) = self.root.parent(prefix, w)? else {
            return Ok(out);
        };
        if !prefix.is_empty() && prefix.ends_with(':') {
            if remove_prefix {
                let bytes = node.measure(w)?;
                w.copy(bytes)?;
                out.root = node.clone();
                out.root.name = w.text("ROOT")?;
                out.root.description.clear();
            } else {
                let n = prefix
                    .len()
                    .checked_sub(node.name.len() + 1)
                    .ok_or_else(|| invalid("invalid parameter copy prefix"))?;
                out.root.insert_node_inner(node, &prefix[..n], w)?;
            }
        } else {
            let suffix = ParamNode::suffix(prefix);
            let parent_prefix = &prefix[..prefix.len() - suffix.len()];
            for child in &node.nodes {
                w.consume(add(child.name.len(), suffix.len())?)?;
                if child.name.starts_with(suffix) {
                    if remove_prefix {
                        let bytes = child.measure(w)?;
                        w.copy(bytes)?;
                        let mut tmp = child.clone();
                        tmp.name = w.text(&child.name[suffix.len()..])?;
                        out.root.insert_node_inner(&tmp, "", w)?;
                    } else {
                        out.root.insert_node_inner(child, parent_prefix, w)?;
                    }
                }
            }
            for entry in &node.entries {
                w.consume(add(entry.name.len(), suffix.len())?)?;
                if entry.name.starts_with(suffix) {
                    if remove_prefix {
                        let bytes = entry.measure(w)?;
                        w.copy(bytes)?;
                        let mut tmp = entry.clone();
                        tmp.name = w.text(&entry.name[suffix.len()..])?;
                        out.root.insert_entry_inner(&tmp, "", w)?;
                    } else {
                        out.root.insert_entry_inner(entry, parent_prefix, w)?;
                    }
                }
            }
        }
        out.root.measure(w)?;
        Ok(out)
    }
    pub fn copy_subset(&self, subset: &Self) -> Result<Self> {
        Ok(self.copy_subset_with_messages(subset)?.0)
    }
    pub fn copy_subset_with_messages(&self, subset: &Self) -> Result<(Self, Vec<String>)> {
        let mut w = ParamWork::default();
        self.root.measure(&mut w)?;
        subset.root.measure(&mut w)?;
        let mut out = Self::new();
        let mut messages = Vec::new();
        for entry in &subset.root.entries {
            if let Some(i) = self.root.entry_index(&entry.name, &mut w)? {
                out.root
                    .insert_entry_inner(&self.root.entries[i], "", &mut w)?;
            } else {
                message(
                    &mut messages,
                    format!("Trying to copy non-existent parameter entry {}", entry.name),
                    &mut w,
                )?;
            }
        }
        for node in &subset.root.nodes {
            if let Some(i) = self.root.node_index(&node.name, &mut w)? {
                out.root
                    .insert_node_inner(&self.root.nodes[i], "", &mut w)?;
            } else {
                message(
                    &mut messages,
                    format!("Trying to copy non-existent parameter node {}", node.name),
                    &mut w,
                )?;
            }
        }
        out.root.measure(&mut w)?;
        Ok((out, messages))
    }
    /// Add missing values, copying only restrictions relevant to each value type.
    /// Returned diagnostics replace the source logging stream.
    pub fn set_defaults(
        &mut self,
        defaults: &Self,
        prefix: &str,
        show_message: bool,
    ) -> Result<Vec<String>> {
        self.edit(|p, w| p.set_defaults_with_work(defaults, prefix, show_message, w))
    }
    pub(super) fn set_defaults_with_work(
        &mut self,
        defaults: &Self,
        prefix: &str,
        show_message: bool,
        w: &mut ParamWork,
    ) -> Result<Vec<String>> {
        let prefix2 = normalized_prefix(prefix, w)?;
        let mut messages = Vec::new();
        let mut path = String::new();
        for item in ParamIterator::with_work(&defaults.root, w)? {
            let key = join(&prefix2, &item.key, w)?;
            if self.root.find_entry_recursive_with_work(&key, w)?.is_none() {
                let bytes = item.entry.measure(w)?;
                w.copy(bytes)?;
                let mut entry = ParamEntry {
                    value: item.entry.value.clone(),
                    description: item.entry.description.clone(),
                    ..ParamEntry::default()
                };
                for tag in &item.entry.tags {
                    no_comma(tag, w)?;
                    w.copy(add(tag.len(), 128)?)?;
                    ordered_lookup_cost(tag, entry.tags.len(), w)?;
                    entry.tags.insert(tag.clone());
                }
                match &entry.value {
                    ParamValue::String(_) | ParamValue::StringList(_) => {
                        w.slots::<String>(item.entry.valid_strings.len())?;
                        for s in &item.entry.valid_strings {
                            no_comma(s, w)?;
                            w.copy(s.len())?;
                        }
                        entry.valid_strings = item.entry.valid_strings.clone();
                    }
                    ParamValue::Integer(_) | ParamValue::IntegerList(_) => {
                        entry.min_int = item.entry.min_int;
                        entry.max_int = item.entry.max_int;
                    }
                    ParamValue::Float(_) | ParamValue::FloatList(_) => {
                        entry.min_float = item.entry.min_float;
                        entry.max_float = item.entry.max_float;
                    }
                    _ => {}
                }
                self.root.insert_entry_inner(&entry, &key, w)?;
                if show_message {
                    message(
                        &mut messages,
                        format!(
                            "Setting {key} to {}",
                            entry.value.to_stream_text_with_work(w)?
                        ),
                        w,
                    )?;
                }
            }
            for event in &item.trace {
                trace_path(&mut path, event, false, w)?;
                let real = path.strip_suffix(':').unwrap_or(&path);
                if !real.is_empty() {
                    let query = join(prefix, real, w)?;
                    if self.section_description_with_work(&query, w)?.is_empty() {
                        let desc = defaults.section_description_with_work(real, w)?;
                        let target = join(&prefix2, real, w)?;
                        self.set_section_description_inner(&target, desc, w)?;
                    }
                }
            }
        }
        Ok(messages)
    }

    /// Preserve source prefix lookup behavior; return unknown-parameter warnings.
    pub fn check_defaults(&self, name: &str, defaults: &Self, prefix: &str) -> Result<Vec<String>> {
        self.check_defaults_with_work(name, defaults, prefix, &mut ParamWork::default())
    }
    pub(super) fn check_defaults_with_work(
        &self,
        name: &str,
        defaults: &Self,
        prefix: &str,
        w: &mut ParamWork,
    ) -> Result<Vec<String>> {
        w.consume(name.len())?;
        defaults.root.measure(w)?;
        let prefix2 = normalized_prefix(prefix, w)?;
        let selected = self.copy_with_work(&prefix2, true, w)?;
        let mut warnings = Vec::new();
        for item in ParamIterator::with_work(&selected.root, w)? {
            if defaults
                .root
                .find_entry_recursive_with_work(&item.key, w)?
                .is_none()
            {
                message(
                    &mut warnings,
                    format!("{name}: unknown parameter '{}'", item.key),
                    w,
                )?;
            }
            let lookup = join(&prefix2, &item.key, w)?;
            if let Some(default) = defaults.root.find_entry_recursive_with_work(&lookup, w)? {
                if default.value.value_type() != item.entry.value.value_type() {
                    return Err(invalid(format!(
                        "{name}: parameter '{}' has the wrong value type",
                        item.key
                    )));
                }
                let bytes = default.measure(w)?;
                w.copy(bytes)?;
                let mut candidate = default.clone();
                let bytes = item.entry.value.measure(w)?;
                w.copy(bytes)?;
                candidate.value = item.entry.value.clone();
                if let Some(message) = candidate.valid_with_work(w)? {
                    return Err(invalid(message));
                }
            }
        }
        Ok(warnings)
    }

    pub fn update(&mut self, old: &Self, add_unknown: bool) -> Result<ParamUpdateReport> {
        self.update_with_options(
            old,
            ParamUpdateOptions {
                add_unknown,
                ..ParamUpdateOptions::default()
            },
        )
    }
    pub fn update_with_options(
        &mut self,
        old: &Self,
        options: ParamUpdateOptions,
    ) -> Result<ParamUpdateReport> {
        self.edit(|p, w| {
            let mut report = ParamUpdateReport {
                success: true,
                messages: Vec::new(),
            };
            for item in ParamIterator::with_work(&old.root, w)? {
                let mut target = String::new();
                if let Some(existing) = p.root.find_entry_recursive_with_work(&item.key, w)? {
                    let protected = item.key.ends_with(":version")
                        || (item.key.ends_with(":type")
                            && item.key.bytes().filter(|b| *b == b':').count() >= 2);
                    if protected {
                        let bytes = add(existing.value.measure(w)?, item.entry.value.measure(w)?)?;
                        w.consume(bytes)?;
                        if existing.value != item.entry.value {
                            message(
                                &mut report.messages,
                                format!("Keeping current protected parameter '{}'", item.key),
                                w,
                            )?;
                        }
                        continue;
                    }
                    target = w.text(&item.key)?;
                } else {
                    let old_entry = old
                        .root
                        .find_entry_recursive_with_work(&item.key, w)?
                        .ok_or_else(|| missing(&item.key))?;
                    let suffix = join(":", &old_entry.name, w)?;
                    let mut matched = None;
                    for candidate in ParamIterator::with_work(&p.root, w)? {
                        w.consume(candidate.key.len())?;
                        if candidate.key.ends_with(&suffix) {
                            if matched.is_some() {
                                matched = None;
                                break;
                            }
                            matched = Some(candidate.key);
                        }
                    }
                    if let Some(key) = matched {
                        message(
                            &mut report.messages,
                            format!("Found '{}' as '{key}' in new param", item.key),
                            w,
                        )?;
                        target = key;
                    }
                    if target.is_empty() {
                        if options.fail_on_unknown_parameters {
                            report.success = false;
                            message(
                                &mut report.messages,
                                format!("Unknown parameter '{}'", item.key),
                                w,
                            )?;
                        } else if options.add_unknown {
                            p.root
                                .insert_entry_inner(old_entry, ancestor(&item.key), w)?;
                            message(
                                &mut report.messages,
                                format!("Adding unknown parameter '{}'", item.key),
                                w,
                            )?;
                        } else if options.verbose {
                            message(
                                &mut report.messages,
                                format!("Ignoring unknown parameter '{}'", item.key),
                                w,
                            )?;
                        }
                        continue;
                    }
                }
                let existing = p
                    .root
                    .find_entry_recursive_with_work(&target, w)?
                    .ok_or_else(|| missing(&target))?;
                if existing.value.value_type() != item.entry.value.value_type() {
                    message(
                        &mut report.messages,
                        format!("Parameter '{}' has changed value type", item.key),
                        w,
                    )?;
                    if options.fail_on_invalid_values {
                        report.success = false;
                    }
                    continue;
                }
                let bytes = add(existing.measure(w)?, item.entry.value.measure(w)?)?;
                w.consume(bytes)?;
                if existing.value == item.entry.value {
                    continue;
                }
                w.copy(bytes)?;
                let mut candidate = existing.clone();
                candidate.value = item.entry.value.clone();
                if let Some(error) = candidate.valid_with_work(w)? {
                    message(&mut report.messages, error, w)?;
                    if options.fail_on_invalid_values {
                        report.success = false;
                    }
                } else {
                    p.root
                        .insert_entry_inner(&candidate, ancestor(&target), w)?;
                    if options.verbose {
                        message(
                            &mut report.messages,
                            format!("Overriding default parameter '{target}'"),
                            w,
                        )?;
                    }
                }
            }
            Ok(report)
        })
    }
    pub fn merge(&mut self, other: &Self) -> Result<()> {
        self.edit(|p, w| {
            let mut path = String::new();
            for item in ParamIterator::with_work(&other.root, w)? {
                let prefix = ancestor(&item.key);
                if p.root
                    .find_entry_recursive_with_work(&item.key, w)?
                    .is_none()
                {
                    p.root.insert_entry_inner(item.entry, prefix, w)?;
                }
                for event in &item.trace {
                    trace_path(&mut path, event, true, w)?;
                    let real = path.strip_suffix(':').unwrap_or(&path);
                    if !real.is_empty() {
                        let query = join(prefix, real, w)?;
                        if p.section_description_with_work(&query, w)?.is_empty() {
                            p.set_section_description_inner(
                                real,
                                other.section_description_with_work(real, w)?,
                                w,
                            )?;
                        }
                    }
                }
            }
            Ok(())
        })
    }
    /// `arguments[0]` is the executable name and is skipped, as in argc/argv.
    pub fn parse_command_line(&mut self, arguments: &[String], prefix: &str) -> Result<()> {
        self.edit(|p, w| {
            preflight_arguments(arguments, w)?;
            let prefix = normalized_prefix(prefix, w)?;
            let misc = join(&prefix, "misc", w)?;
            let mut i = 1;
            while i < arguments.len() {
                w.consume(1)?;
                let arg = &arguments[i];
                let next = arguments.get(i + 1).map_or("", String::as_str);
                if is_option(arg) {
                    let value = if is_option(next) {
                        ""
                    } else {
                        i += 1;
                        next
                    };
                    let entry = ParamEntry::new_with_work(
                        arg,
                        ParamValue::String(w.text(value)?),
                        "",
                        &[],
                        w,
                    )?;
                    p.root.insert_entry_inner(&entry, &prefix, w)?;
                } else {
                    p.append_cli_list(&misc, arg, w)?;
                }
                i += 1;
            }
            Ok(())
        })
    }
    /// Mapped parsing prioritizes multiple, then no-argument, then one-argument
    /// registrations. Unrecognized options and plain arguments have separate lists.
    pub fn parse_command_line_mapped(
        &mut self,
        arguments: &[String],
        options: &CommandLineOptions,
    ) -> Result<()> {
        self.edit(|p, w| {
            preflight_arguments(arguments, w)?;
            for map in [
                &options.multiple_arguments,
                &options.no_argument,
                &options.one_argument,
            ] {
                w.consume(map.len())?;
                for (key, value) in map {
                    w.consume(add(key.len(), value.len())?)?;
                    key_check(value, w)?;
                }
            }
            key_check(&options.misc, w)?;
            key_check(&options.unknown, w)?;
            let mut i = 1;
            while i < arguments.len() {
                w.consume(1)?;
                let arg = &arguments[i];
                let next = arguments.get(i + 1).map_or("", String::as_str);
                // Conservatively meter the ordered-map key comparisons.
                for map in [
                    &options.multiple_arguments,
                    &options.no_argument,
                    &options.one_argument,
                ] {
                    w.consume(mul(
                        arg.len() + 1,
                        map.len().saturating_add(1).ilog2() as usize + 1,
                    )?)?;
                }
                let item = if let Some(key) = options.multiple_arguments.get(arg) {
                    let mut values = Vec::new();
                    while i + 1 < arguments.len() && !is_option(&arguments[i + 1]) {
                        i += 1;
                        w.slots::<String>(1)?;
                        values.push(w.text(&arguments[i])?);
                    }
                    Some((key, ParamValue::StringList(values)))
                } else if let Some(key) = options.no_argument.get(arg) {
                    Some((key, ParamValue::String(w.text("true")?)))
                } else if let Some(key) = options.one_argument.get(arg) {
                    let value = if is_option(next) {
                        ""
                    } else {
                        i += 1;
                        next
                    };
                    Some((key, ParamValue::String(w.text(value)?)))
                } else {
                    None
                };
                if let Some((key, value)) = item {
                    let entry = ParamEntry::new_with_work("", value, "", &[], w)?;
                    p.root.insert_entry_inner(&entry, key, w)?;
                } else {
                    p.append_cli_list(
                        if is_option(arg) {
                            &options.unknown
                        } else {
                            &options.misc
                        },
                        arg,
                        w,
                    )?;
                }
                i += 1;
            }
            Ok(())
        })
    }
    fn append_cli_list(&mut self, key: &str, value: &str, w: &mut ParamWork) -> Result<()> {
        if self.root.find_entry_recursive_with_work(key, w)?.is_some() {
            let e = self.entry_mut(key, w)?;
            match &mut e.value {
                ParamValue::StringList(v) => {
                    w.slots::<String>(1)?;
                    v.push(w.text(value)?);
                    Ok(())
                }
                _ => Err(invalid(
                    "command-line accumulation requires a string-list parameter",
                )),
            }
        } else {
            w.slots::<String>(1)?;
            let entry = ParamEntry::new_with_work(
                "",
                ParamValue::StringList(vec![w.text(value)?]),
                "",
                &[],
                w,
            )?;
            self.root.insert_entry_inner(&entry, key, w)
        }
    }
}
fn preflight_arguments(arguments: &[String], w: &mut ParamWork) -> Result<()> {
    w.consume(arguments.len())?;
    if arguments.len() > MAX_PARAM_ENTRIES {
        return Err(invalid("command-line argument count limit exceeded"));
    }
    for a in arguments {
        w.consume(a.len())?;
    }
    Ok(())
}
fn is_option(arg: &str) -> bool {
    let bytes = arg.as_bytes();
    bytes.len() >= 2 && bytes[0] == b'-' && !bytes[1].is_ascii_digit()
}

fn remove_matching<T>(
    values: &mut Vec<T>,
    prefix: &str,
    name: impl Fn(&T) -> &str,
    w: &mut ParamWork,
) -> Result<bool> {
    for value in values.iter() {
        w.consume(add(
            1,
            add(
                mul(2, add(name(value).len(), prefix.len())?)?,
                size_of::<T>(),
            )?,
        )?)?;
    }
    let before = values.len();
    values.retain(|value| !name(value).starts_with(prefix));
    Ok(values.len() != before)
}
