// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! DefaultParamHandler lifecycle without inheritance or stored callbacks.

use super::{
    MAX_PARAM_BYTES, MAX_PARAM_NODES, Param, ParamBudget, ParamIterator, ParamValue, ParamWork,
    add, invalid, mul,
};
use crate::Result;
use std::mem::size_of;

/// A metadata container [`DefaultParamHandler::write_parameters_to_meta_values`]
/// can copy a parameter tree into.
///
/// The destination is a trait rather than a named type because `param` sits
/// below `metadata` in the module graph: `metadata` reaches `chemistry`, which
/// reaches `param`, so naming a metadata type here would close a module cycle
/// and block the workspace split. The dependency is the same one inverted -
/// `metadata::MetaInfo` implements this, and the public path stays where the
/// source class puts it.
pub trait ParameterMetaSink: Default {
    /// The value this container stores under a key.
    type Value;

    /// Bytes the container already holds, charging `budget` for the walk.
    fn measure_existing(&self, budget: &mut ParamBudget<'_>) -> Result<usize>;

    /// This container's value for one parameter leaf.
    ///
    /// A value the container cannot represent - a nonfinite float, for native
    /// metadata - is refused here, before anything is staged, so that a refusal
    /// leaves the destination untouched.
    fn value_of(value: &ParamValue, budget: &mut ParamBudget<'_>) -> Result<Self::Value>;

    /// Stage `value` under `key`, replacing an equal key.
    fn stage(&mut self, key: String, value: Self::Value);

    /// How many entries are staged.
    fn staged(&self) -> usize;

    /// Move every entry of `staged` into `self`, leaving `staged` empty.
    fn absorb(&mut self, staged: &mut Self);
}

/// Current/default parameter trees and validation policy. For derived typed
/// settings, use a `_with` update that returns newly constructed member state.
#[derive(Clone, Debug, PartialEq)]
pub struct DefaultParamHandler {
    parameters: Param,
    defaults: Param,
    name: String,
    subsections: Vec<String>,
    check_defaults: bool,
    warn_empty_defaults: bool,
}
impl DefaultParamHandler {
    pub fn new(name: &str) -> Result<Self> {
        Ok(Self {
            parameters: Param::new(),
            defaults: Param::new(),
            name: ParamWork::default().text(name)?,
            subsections: Vec::new(),
            check_defaults: true,
            warn_empty_defaults: true,
        })
    }
    pub fn parameters(&self) -> &Param {
        &self.parameters
    }
    pub fn defaults(&self) -> &Param {
        &self.defaults
    }
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn subsections(&self) -> &[String] {
        &self.subsections
    }
    pub fn check_defaults(&self) -> bool {
        self.check_defaults
    }
    pub fn warn_empty_defaults(&self) -> bool {
        self.warn_empty_defaults
    }
    pub fn set_check_defaults(&mut self, enabled: bool) {
        self.check_defaults = enabled;
    }
    pub fn set_warn_empty_defaults(&mut self, enabled: bool) {
        self.warn_empty_defaults = enabled;
    }
    pub fn set_name(&mut self, name: &str) -> Result<()> {
        self.name = ParamWork::default().text(name)?;
        Ok(())
    }
    /// Configure defaults without changing current parameters. Finish setup with
    /// `defaults_to_parameters` or its typed-state callback variant.
    pub fn set_defaults(&mut self, defaults: Param) -> Result<()> {
        defaults.root.measure(&mut ParamWork::default())?;
        self.defaults = defaults;
        Ok(())
    }
    /// Source convention: entries normally omit the final colon. Values are
    /// retained verbatim, including duplicates; validation appends one colon.
    pub fn set_subsections(&mut self, subsections: &[String]) -> Result<()> {
        let mut work = ParamWork::default();
        if subsections.len() > MAX_PARAM_NODES {
            return Err(invalid("too many parameter-handler subsections"));
        }
        work.copy(mul(subsections.len(), size_of::<String>())?)?;
        let mut replacement = Vec::new();
        replacement
            .try_reserve_exact(subsections.len())
            .map_err(|_| invalid("subsection allocation failed"))?;
        for subsection in subsections {
            replacement.push(work.text(subsection)?);
        }
        self.subsections = replacement;
        Ok(())
    }
    pub fn checked_clone(&self) -> Result<Self> {
        let mut work = ParamWork::default();
        let bytes = self.measure(&mut work)?;
        work.copy(bytes)?;
        Ok(self.clone())
    }
    /// C++ equality uses Param's source name/value predicate, ignoring parameter
    /// descriptions/restrictions and order where its Param comparison does.
    pub fn source_equal(&self, other: &Self) -> Result<bool> {
        let mut work = ParamWork::default();
        self.measure(&mut work)?;
        other.measure(&mut work)?;
        Ok(self.name == other.name
            && self.subsections == other.subsections
            && self.check_defaults == other.check_defaults
            && self.warn_empty_defaults == other.warn_empty_defaults
            && self
                .parameters
                .root
                .source_equal_inner(&other.parameters.root, &mut work)?
            && self
                .defaults
                .root
                .source_equal_inner(&other.defaults.root, &mut work)?)
    }
    /// Replace current parameters, fill missing defaults, and return warnings.
    /// Unknown keys are warnings; type and restriction violations are errors.
    pub fn set_parameters(&mut self, parameters: &Param) -> Result<Vec<String>> {
        self.set_parameters_with(parameters, |_| Ok(()))
            .map(|(_, warnings)| warnings)
    }
    /// Validate a staged tree and build replacement typed state before committing.
    /// The callback must be pure with respect to existing external state: callback
    /// side effects cannot be rolled back. Install its returned value on success.
    pub fn set_parameters_with<T>(
        &mut self,
        parameters: &Param,
        update_members: impl FnOnce(&Param) -> Result<T>,
    ) -> Result<(T, Vec<String>)> {
        let mut work = ParamWork::default();
        self.measure(&mut work)?;
        let mut staged = parameters.clone_with_work(&mut work)?;
        staged.set_defaults_with_work(&self.defaults, "", false, &mut work)?;
        let mut warnings = Vec::new();
        if self.check_defaults {
            if self.defaults.is_empty() && self.warn_empty_defaults {
                warning(&mut warnings, &self.name, "", "empty defaults", &mut work)?;
            }
            let validation_warnings = if self.subsections.is_empty() {
                staged.check_defaults_with_work(&self.name, &self.defaults, "", &mut work)?
            } else {
                let mut checked = staged.clone_with_work(&mut work)?;
                for subsection in &self.subsections {
                    work.copy(add(subsection.len(), 1)?)?;
                    let prefix = format!("{subsection}:");
                    checked.remove_all_with_work(&prefix, &mut work)?;
                }
                checked.check_defaults_with_work(&self.name, &self.defaults, "", &mut work)?
            };
            work.slots::<String>(validation_warnings.len())?;
            warnings.extend(validation_warnings);
        }
        staged.root.measure(&mut work)?;
        let members = update_members(&staged)?;
        self.parameters = staged;
        Ok((members, warnings))
    }
    /// Fill missing current values, preserving existing values. Unlike
    /// `set_parameters`, source initialization does not validate restrictions.
    pub fn defaults_to_parameters(&mut self) -> Result<Vec<String>> {
        self.defaults_to_parameters_with(|_| Ok(()))
            .map(|(_, warnings)| warnings)
    }
    pub fn defaults_to_parameters_with<T>(
        &mut self,
        update_members: impl FnOnce(&Param) -> Result<T>,
    ) -> Result<(T, Vec<String>)> {
        let mut work = ParamWork::default();
        self.measure(&mut work)?;
        let mut warnings = Vec::new();
        // The source breaks at the first missing description, despite plural text.
        for item in ParamIterator::with_work(&self.defaults.root, &mut work)? {
            if item.entry.description.is_empty() {
                warning(
                    &mut warnings,
                    &self.name,
                    &item.key,
                    "missing description",
                    &mut work,
                )?;
                break;
            }
        }
        let mut staged = self.parameters.clone_with_work(&mut work)?;
        staged.set_defaults_with_work(&self.defaults, "", false, &mut work)?;
        staged.root.measure(&mut work)?;
        let members = update_members(&staged)?;
        self.parameters = staged;
        Ok((members, warnings))
    }

    /// Copy parameter values to metadata using LEAF names, not full paths.
    /// Equal leaf keys are overwritten in source iteration order. Nonfinite
    /// parameter floats cannot enter native MetaInfo and fail atomically.
    pub fn write_parameters_to_meta_values<S: ParameterMetaSink>(
        parameters: &Param,
        metadata: &mut S,
        prefix: &str,
    ) -> Result<()> {
        let mut work = ParamWork::default();
        let existing = metadata.measure_existing(&mut ParamBudget::new(&mut work))?;
        work.allocation(existing)?;
        let mut prefix = work.text(prefix)?;
        if !prefix.is_empty() && !prefix.ends_with(':') {
            work.copy(1)?;
            prefix.push(':');
        }
        let mut replacements = S::default();
        let mut replacement_key_bytes = 0usize;
        for item in ParamIterator::with_work(&parameters.root, &mut work)? {
            let bytes = item.entry.value.measure(&mut work)?;
            work.copy(add(bytes, 128)?)?;
            let value = S::value_of(&item.entry.value, &mut ParamBudget::new(&mut work))?;
            work.copy(add(prefix.len(), item.entry.name.len())?)?;
            let key = format!("{prefix}{}", item.entry.name);
            replacement_key_bytes = add(replacement_key_bytes, key.len())?;
            // A BTree node compares at most eleven keys. Bound comparison work
            // using its minimum branching factor, including long common prefixes.
            let mut nodes = replacements.staged();
            let mut levels = 1;
            while nodes >= 6 {
                nodes /= 6;
                levels += 1;
            }
            work.consume(mul(add(key.len(), 1)?, mul(12, levels)?)?)?;
            replacements.stage(key, value);
        }
        work.consume(mul(add(existing, replacement_key_bytes)?, 2)?)?;
        // All fallible work is complete. Moving entries preserves unrelated
        // existing metadata and avoids cloning the entire destination map.
        metadata.absorb(&mut replacements);
        Ok(())
    }
    fn measure(&self, work: &mut ParamWork) -> Result<usize> {
        work.consume(add(
            self.name.len(),
            mul(self.subsections.len(), size_of::<String>())?,
        )?)?;
        let mut bytes = add(size_of::<Self>(), self.name.len())?;
        bytes = add(bytes, mul(self.subsections.len(), size_of::<String>())?)?;
        for s in &self.subsections {
            work.consume(s.len())?;
            bytes = add(bytes, s.len())?;
        }
        bytes = add(bytes, self.parameters.root.measure(work)?)?;
        bytes = add(bytes, self.defaults.root.measure(work)?)?;
        if bytes > MAX_PARAM_BYTES {
            return Err(invalid("parameter-handler payload limit exceeded"));
        }
        Ok(bytes)
    }
}

fn warning(
    messages: &mut Vec<String>,
    name: &str,
    key: &str,
    kind: &str,
    work: &mut ParamWork,
) -> Result<()> {
    work.copy(add(add(name.len(), key.len())?, 256)?)?;
    let text = if kind == "empty defaults" {
        format!("Warning: No default parameters for DefaultParameterHandler '{name}' specified!")
    } else {
        format!(
            "Warning: no default parameter description for parameters '{key},' of DefaultParameterHandler '{name}' given!"
        )
    };
    messages.push(text);
    Ok(())
}
