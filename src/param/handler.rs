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
/// The destination is a trait rather than a named type because `param` is the
/// lower module and `metadata` is the one that names it, to implement this. A
/// `param` that named `MetaInfo` would close the cycle between the two.
///
/// Written the other way round it closed a longer one, `param -> metadata ->
/// chemistry -> param`, which is the cycle the module gate first refused;
/// `metadata -> chemistry` has since been cut as well. The dependency is the
/// same one inverted, and the public path stays on the class the source puts
/// it on.
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

/// The parameter contract of a type whose source class derives from
/// `DefaultParamHandler`.
///
/// The source keeps two parameter trees in every such object (`defaults_` and
/// `param_`) and re-derives its typed members from `param_` in the virtual hook
/// `updateMembers_`, which `setParameters` calls after checking the new tree.
/// A Rust port keeps only the typed members: the defaults are a function, the
/// current tree is rebuilt from the members on demand, and `updateMembers_` is a
/// constructor from a complete, checked tree. The [`DefaultParamHandler`] value
/// that the source object would hold stays available through
/// [`DefaultParamHandling::from_param_with_handler`].
///
/// | Source | Rust |
/// |---|---|
/// | `DefaultParamHandler(name)`, `error_name_`, `getName()` | [`NAME`](Self::NAME) |
/// | `subsections_` | [`SUBSECTIONS`](Self::SUBSECTIONS) |
/// | `check_defaults_`, `warn_empty_defaults_` | [`CHECK_DEFAULTS`](Self::CHECK_DEFAULTS), [`WARN_EMPTY_DEFAULTS`](Self::WARN_EMPTY_DEFAULTS) |
/// | `defaults_`, `getDefaults()` | [`defaults`](Self::defaults) |
/// | `param_`, `getParameters()` | [`to_param`](Self::to_param) |
/// | `updateMembers_()` | [`update_members`](Self::update_members) |
/// | `setParameters(p)` on a new object | [`from_param`](Self::from_param), [`from_param_with_warnings`](Self::from_param_with_warnings), [`from_param_with_handler`](Self::from_param_with_handler) |
/// | `setParameters(p)` on an existing object | [`set_parameters`](Self::set_parameters) |
/// | `defaultsToParam_()` in the constructor | `Default`, whose [`to_param`](Self::to_param) equals [`defaults`](Self::defaults) |
///
/// # Adopting the trait
///
/// A type written in the established shape - inherent `defaults()`,
/// `to_param(&self)`, `from_param`, `from_param_with_warnings` and a private
/// `from_complete_param(&Param) -> Result<Self>` that builds the typed members
/// from a tree holding every default (`PeakPickerHiRes` in
/// `src/processing/peak_picking.rs` is the reference) - adopts the trait with
/// four delegating lines; its inherent `from_param*` can then delegate to the
/// provided methods so that the handler is set up in one place:
///
/// ```
/// use openms::Result;
/// use openms::param::{DefaultParamHandling, Param, ParamValue};
///
/// #[derive(Debug, Default, PartialEq)]
/// struct Mower {
///     threshold: f64,
/// }
///
/// impl Mower {
///     fn from_complete_param(param: &Param) -> Result<Self> {
///         Ok(Self { threshold: param.value("threshold")?.to_f64()? })
///     }
/// }
///
/// impl DefaultParamHandling for Mower {
///     const NAME: &'static str = "Mower";
///     fn defaults() -> Result<Param> {
///         Self::default().to_param()
///     }
///     fn to_param(&self) -> Result<Param> {
///         let mut param = Param::new();
///         param.set_value("threshold", ParamValue::Float(self.threshold), "Cut-off", &[])?;
///         Ok(param)
///     }
///     fn update_members(param: &Param) -> Result<Self> {
///         Self::from_complete_param(param)
///     }
/// }
///
/// let mut changed = Param::new();
/// changed.set_value("threshold", ParamValue::Float(0.5), "", &[])?;
/// changed.set_value("typo", ParamValue::Float(1.0), "", &[])?;
/// let (mower, warnings) = Mower::from_param_with_warnings(&changed)?;
/// assert_eq!(mower.threshold, 0.5);
/// assert_eq!(warnings.len(), 1); // the unknown `typo`
/// assert_eq!(Mower::default().to_param()?, Mower::defaults()?);
/// # Ok::<(), openms::Error>(())
/// ```
///
/// # Invariants an implementation keeps
///
/// - `Self::default().to_param()` equals `Self::defaults()` when the type is
///   `Default`: the source constructor ends in `defaultsToParam_()`.
/// - `Self::from_param(&p)?.to_param()` equals `p` merged with the defaults for
///   every `p` that passes the check, in the defaults' order, with the defaults'
///   descriptions and restrictions: the source's `getParameters()` after
///   `setParameters(p)`.
///
/// # Native differences
///
/// - Atomic: the source assigns `param_` before `checkDefaults` can throw, so a
///   refused tree is left in `param_` while the typed members keep their old
///   values; here a refused tree changes nothing.
/// - The source writes its warnings to `OPENMS_LOG_WARN`; here they are
///   returned, in the wording of [`DefaultParamHandler::set_parameters`].
/// - An `updateMembers_` that throws (`SavitzkyGolayFilter` does) leaves the
///   source object with the new `param_` and partially updated members; here
///   [`update_members`](Self::update_members) returns the error and nothing is
///   changed.
pub trait DefaultParamHandling: Sized {
    /// The class name the source passes to the `DefaultParamHandler`
    /// constructor, used in every parameter diagnostic (source `error_name_`,
    /// `getName()`).
    const NAME: &'static str;

    /// Subsections whose parameters are handed on to other classes and so are
    /// not checked against [`defaults`](Self::defaults), without the trailing
    /// colon (source `subsections_`, whose `@note` asks for exactly that).
    const SUBSECTIONS: &'static [&'static str] = &[];

    /// Whether `set_parameters` checks the given tree against the defaults
    /// (source `check_defaults_`, `true` unless a class turns it off).
    const CHECK_DEFAULTS: bool = true;

    /// Whether an empty [`defaults`](Self::defaults) tree produces the source
    /// warning "No default parameters for DefaultParameterHandler ..." (source
    /// `warn_empty_defaults_`, `true` unless a class turns it off).
    const WARN_EMPTY_DEFAULTS: bool = true;

    /// The source parameter defaults: every name, value, value type,
    /// description, restriction and tag, in declaration order (source
    /// `getDefaults()`, the member `defaults_` that the constructor fills).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    fn defaults() -> Result<Param>;

    /// The parameters of this value in the layout of
    /// [`defaults`](Self::defaults), carrying the current values (source
    /// `getParameters()`, the member `param_`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`](crate::Error::InvalidValue) when a typed member
    /// has no representation as the source parameter (a count beyond the source's
    /// 32-bit integer, for example).
    fn to_param(&self) -> Result<Param>;

    /// The typed members of a value configured by `param` (source
    /// `updateMembers_`, the hook `setParameters` and `defaultsToParam_` call
    /// last).
    ///
    /// `param` already holds every default and has passed the restriction
    /// check, so an implementation reads every key it needs without fallbacks.
    /// Members that are not parameters take their `Default` values; see
    /// [`keep_native_members`](Self::keep_native_members).
    ///
    /// # Errors
    ///
    /// Whatever the source's `updateMembers_` throws, and an
    /// [`Error::InvalidValue`](crate::Error::InvalidValue) for a value the typed
    /// member cannot hold.
    fn update_members(param: &Param) -> Result<Self>;

    /// `updated` with the members that are not parameters - native resource
    /// ceilings, compatibility profiles, attached databases - taken from `self`,
    /// as [`set_parameters`](Self::set_parameters) must leave them. The default
    /// returns `updated` unchanged, for types whose members are all parameters.
    fn keep_native_members(&self, updated: Self) -> Self {
        updated
    }

    /// A [`DefaultParamHandler`] set up as the source constructor sets up its
    /// base: named [`NAME`](Self::NAME), with [`defaults`](Self::defaults),
    /// [`SUBSECTIONS`](Self::SUBSECTIONS) and both check flags, and no current
    /// parameters yet.
    ///
    /// # Errors
    ///
    /// As [`defaults`](Self::defaults), and a resource failure of the handler.
    fn handler() -> Result<DefaultParamHandler> {
        let mut handler = DefaultParamHandler::new(Self::NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        if !Self::SUBSECTIONS.is_empty() {
            let subsections: Vec<String> = Self::SUBSECTIONS
                .iter()
                .map(|name| (*name).to_owned())
                .collect();
            handler.set_subsections(&subsections)?;
        }
        handler.set_check_defaults(Self::CHECK_DEFAULTS);
        handler.set_warn_empty_defaults(Self::WARN_EMPTY_DEFAULTS);
        Ok(handler)
    }

    /// A value configured by `parameters`, as a default-constructed source
    /// object after `setParameters(parameters)`.
    ///
    /// Missing parameters take their defaults. Unknown parameters are ignored,
    /// as the source only warns about them; see
    /// [`from_param_with_warnings`](Self::from_param_with_warnings).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`](crate::Error::InvalidValue) when a parameter has
    /// the wrong value type or breaks its restriction (the source's
    /// `Exception::InvalidParameter`), and whatever
    /// [`update_members`](Self::update_members) returns.
    fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`from_param`](Self::from_param), also returning the warnings the
    /// source logs: unknown parameters, and empty defaults when
    /// [`WARN_EMPTY_DEFAULTS`](Self::WARN_EMPTY_DEFAULTS) is set.
    ///
    /// # Errors
    ///
    /// As [`from_param`](Self::from_param).
    fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let (value, _, warnings) = Self::from_param_with_handler(parameters)?;
        Ok((value, warnings))
    }

    /// As [`from_param_with_warnings`](Self::from_param_with_warnings), also
    /// returning the [`DefaultParamHandler`] the source object would hold: its
    /// [`parameters`](DefaultParamHandler::parameters) are `parameters` merged
    /// with the defaults, including unknown keys, exactly as the source's
    /// `getParameters()` returns them after `setParameters`.
    ///
    /// # Errors
    ///
    /// As [`from_param`](Self::from_param).
    fn from_param_with_handler(
        parameters: &Param,
    ) -> Result<(Self, DefaultParamHandler, Vec<String>)> {
        let mut handler = Self::handler()?;
        let (value, warnings) = handler.set_parameters_with(parameters, Self::update_members)?;
        Ok((value, handler, warnings))
    }

    /// Replace the parameter values of this value, as source `setParameters`
    /// on an existing object, and return the warnings.
    ///
    /// As in the source, a parameter missing from `parameters` returns to its
    /// default rather than keeping its current value. Members that are not
    /// parameters are kept (see
    /// [`keep_native_members`](Self::keep_native_members)). On error `self` is
    /// unchanged.
    ///
    /// # Errors
    ///
    /// As [`from_param`](Self::from_param).
    fn set_parameters(&mut self, parameters: &Param) -> Result<Vec<String>> {
        let (updated, warnings) = Self::from_param_with_warnings(parameters)?;
        let updated = self.keep_native_members(updated);
        *self = updated;
        Ok(warnings)
    }
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
