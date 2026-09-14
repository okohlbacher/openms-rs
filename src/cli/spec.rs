// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Tool parameter registration: the native form of the source `register*_`,
//! `setValidStrings_`, `setValidFormats_`, `setMin*_`/`setMax*_`,
//! `registerSubsection_`, `registerFullParam_`, `addText_` and `addEmptyLine_`
//! calls, and of the parameter part of `getDefaultParameters_`.
//!
//! See `docs/TOPP_CLI_SUPPORT.md` for the supported source subset.

use super::parameter::{ParameterInformation, ParameterType};
use crate::param::{Param, ParamValue};
use crate::{Error, Result};

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// Registered common options the source never puts into a parameter tree or an
/// INI file (`TOPPBase.cpp:2104`).
const NOT_IN_PARAMETER_TREE: [&str; 10] = [
    "ini",
    "-help",
    "-helphelp",
    "instance",
    "write_ini",
    "write_ctd",
    "write_cwl",
    "write_nested_cwl",
    "write_json",
    "write_nested_json",
];

/// Registered parameters and subsections of one tool.
///
/// The source registers parameters by calling protected `register*_` methods on
/// itself during `registerOptionsAndFlags_`. Here a tool fills this builder
/// instead, which keeps registration separate from execution and makes the
/// registered set inspectable without running the tool.
#[derive(Clone, Debug, Default)]
pub struct ToolSpec {
    parameters: Vec<ParameterInformation>,
    subsections: Vec<(String, String)>,
    topp_subsections: Vec<(String, String)>,
}

impl ToolSpec {
    /// Maximum registered parameters for one tool.
    pub const MAX_PARAMETERS: usize = 4096;

    /// An empty registration.
    pub fn new() -> Self {
        Self::default()
    }
    /// Every registered parameter and layout line, in registration order.
    pub fn parameters(&self) -> &[ParameterInformation] {
        &self.parameters
    }
    /// Algorithm subsections as `(name, description)`, in registration order,
    /// from [`register_subsection`](Self::register_subsection).
    pub fn subsections(&self) -> &[(String, String)] {
        &self.subsections
    }
    /// Sections of parameters registered with a `:` in their name, as
    /// `(section, description)`; the source's `subsections_TOPP_`.
    pub fn topp_subsections(&self) -> &[(String, String)] {
        &self.topp_subsections
    }
    /// Registered parameter by name, skipping layout entries.
    pub fn find(&self, name: &str) -> Option<&ParameterInformation> {
        self.parameters
            .iter()
            .find(|p| !p.kind.is_layout() && p.name == name)
    }
    fn find_mut(&mut self, name: &str) -> Result<&mut ParameterInformation> {
        self.parameters
            .iter_mut()
            .find(|p| !p.kind.is_layout() && p.name == name)
            .ok_or_else(|| bad(format!("no registered parameter named '{name}'")))
    }

    fn push(&mut self, parameter: ParameterInformation) -> Result<()> {
        if self.parameters.len() >= Self::MAX_PARAMETERS {
            return Err(bad("tool registers too many parameters"));
        }
        if !parameter.kind.is_layout() {
            if parameter.name.is_empty() {
                return Err(bad("registered parameter needs a name"));
            }
            if self.find(&parameter.name).is_some() {
                return Err(bad(format!(
                    "parameter '{}' is registered twice",
                    parameter.name
                )));
            }
        }
        self.parameters.push(parameter);
        Ok(())
    }

    #[allow(clippy::too_many_arguments)] // mirrors the source register*_ signature
    fn register(
        &mut self,
        name: &str,
        kind: ParameterType,
        argument: &str,
        default_value: ParamValue,
        description: &str,
        required: bool,
        advanced: bool,
    ) -> Result<()> {
        self.push(ParameterInformation::new(
            name,
            kind,
            argument,
            default_value,
            description,
            required,
            advanced,
        ))
    }

    /// Register a free-text option, as `registerStringOption_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an empty or duplicate name or when
    /// [`MAX_PARAMETERS`](Self::MAX_PARAMETERS) would be exceeded.
    pub fn register_string_option(
        &mut self,
        name: &str,
        argument: &str,
        default_value: &str,
        description: &str,
        required: bool,
        advanced: bool,
    ) -> Result<()> {
        self.register(
            name,
            ParameterType::String,
            argument,
            ParamValue::String(default_value.into()),
            description,
            required,
            advanced,
        )
    }
    /// Register an integer option, as `registerIntOption_`.
    ///
    /// # Errors
    ///
    /// As [`register_string_option`](Self::register_string_option).
    pub fn register_int_option(
        &mut self,
        name: &str,
        argument: &str,
        default_value: i64,
        description: &str,
        required: bool,
        advanced: bool,
    ) -> Result<()> {
        self.register(
            name,
            ParameterType::Int,
            argument,
            ParamValue::Integer(default_value),
            description,
            required,
            advanced,
        )
    }
    /// Register a floating-point option, as `registerDoubleOption_`.
    ///
    /// # Errors
    ///
    /// As [`register_string_option`](Self::register_string_option).
    pub fn register_double_option(
        &mut self,
        name: &str,
        argument: &str,
        default_value: f64,
        description: &str,
        required: bool,
        advanced: bool,
    ) -> Result<()> {
        self.register(
            name,
            ParameterType::Double,
            argument,
            ParamValue::Float(default_value),
            description,
            required,
            advanced,
        )
    }
    /// Register a flag, as `registerFlag_`: `false` unless given.
    ///
    /// # Errors
    ///
    /// As [`register_string_option`](Self::register_string_option).
    pub fn register_flag(&mut self, name: &str, description: &str, advanced: bool) -> Result<()> {
        self.register(
            name,
            ParameterType::Flag,
            "",
            ParamValue::String("false".into()),
            description,
            false,
            advanced,
        )
    }
    /// Register an input file, as `registerInputFile_`. `tags` carries source
    /// tags such as `skipexists`, which skips the readability check.
    ///
    /// # Errors
    ///
    /// As [`register_string_option`](Self::register_string_option).
    #[allow(clippy::too_many_arguments)] // mirrors the source register*_ signature
    pub fn register_input_file(
        &mut self,
        name: &str,
        argument: &str,
        default_value: &str,
        description: &str,
        required: bool,
        advanced: bool,
        tags: &[&str],
    ) -> Result<()> {
        self.register(
            name,
            ParameterType::InputFile,
            argument,
            ParamValue::String(default_value.into()),
            description,
            required,
            advanced,
        )?;
        let entry = self.find_mut(name)?;
        entry.tags = tags.iter().map(|t| (*t).to_owned()).collect();
        Ok(())
    }
    /// Register an output file, as `registerOutputFile_`.
    ///
    /// # Errors
    ///
    /// As [`register_string_option`](Self::register_string_option).
    pub fn register_output_file(
        &mut self,
        name: &str,
        argument: &str,
        default_value: &str,
        description: &str,
        required: bool,
        advanced: bool,
    ) -> Result<()> {
        self.register(
            name,
            ParameterType::OutputFile,
            argument,
            ParamValue::String(default_value.into()),
            description,
            required,
            advanced,
        )
    }
    /// Register an output prefix, as `registerOutputPrefix_`.
    ///
    /// # Errors
    ///
    /// As [`register_string_option`](Self::register_string_option).
    pub fn register_output_prefix(
        &mut self,
        name: &str,
        argument: &str,
        default_value: &str,
        description: &str,
        required: bool,
        advanced: bool,
    ) -> Result<()> {
        self.register(
            name,
            ParameterType::OutputPrefix,
            argument,
            ParamValue::String(default_value.into()),
            description,
            required,
            advanced,
        )
    }
    /// Register an output directory, as `registerOutputDir_`.
    ///
    /// # Errors
    ///
    /// As [`register_string_option`](Self::register_string_option).
    pub fn register_output_dir(
        &mut self,
        name: &str,
        argument: &str,
        default_value: &str,
        description: &str,
        required: bool,
        advanced: bool,
    ) -> Result<()> {
        self.register(
            name,
            ParameterType::OutputDir,
            argument,
            ParamValue::String(default_value.into()),
            description,
            required,
            advanced,
        )
    }
    /// Register a string list, as `registerStringList_`.
    ///
    /// # Errors
    ///
    /// As [`register_string_option`](Self::register_string_option).
    pub fn register_string_list(
        &mut self,
        name: &str,
        argument: &str,
        default_value: &[&str],
        description: &str,
        required: bool,
        advanced: bool,
    ) -> Result<()> {
        self.register(
            name,
            ParameterType::StringList,
            argument,
            ParamValue::StringList(default_value.iter().map(|s| (*s).to_owned()).collect()),
            description,
            required,
            advanced,
        )
    }
    /// Register an integer list, as `registerIntList_`.
    ///
    /// # Errors
    ///
    /// As [`register_string_option`](Self::register_string_option).
    pub fn register_int_list(
        &mut self,
        name: &str,
        argument: &str,
        default_value: &[i32],
        description: &str,
        required: bool,
        advanced: bool,
    ) -> Result<()> {
        self.register(
            name,
            ParameterType::IntList,
            argument,
            ParamValue::IntegerList(default_value.to_vec()),
            description,
            required,
            advanced,
        )
    }
    /// Register a floating-point list, as `registerDoubleList_`.
    ///
    /// # Errors
    ///
    /// As [`register_string_option`](Self::register_string_option).
    pub fn register_double_list(
        &mut self,
        name: &str,
        argument: &str,
        default_value: &[f64],
        description: &str,
        required: bool,
        advanced: bool,
    ) -> Result<()> {
        self.register(
            name,
            ParameterType::DoubleList,
            argument,
            ParamValue::FloatList(default_value.to_vec()),
            description,
            required,
            advanced,
        )
    }
    /// Register an input file list, as `registerInputFileList_`.
    ///
    /// # Errors
    ///
    /// As [`register_string_option`](Self::register_string_option).
    #[allow(clippy::too_many_arguments)] // mirrors the source register*_ signature
    pub fn register_input_file_list(
        &mut self,
        name: &str,
        argument: &str,
        default_value: &[&str],
        description: &str,
        required: bool,
        advanced: bool,
        tags: &[&str],
    ) -> Result<()> {
        self.register(
            name,
            ParameterType::InputFileList,
            argument,
            ParamValue::StringList(default_value.iter().map(|s| (*s).to_owned()).collect()),
            description,
            required,
            advanced,
        )?;
        let entry = self.find_mut(name)?;
        entry.tags = tags.iter().map(|t| (*t).to_owned()).collect();
        Ok(())
    }
    /// Register an output file list, as `registerOutputFileList_`.
    ///
    /// # Errors
    ///
    /// As [`register_string_option`](Self::register_string_option).
    pub fn register_output_file_list(
        &mut self,
        name: &str,
        argument: &str,
        default_value: &[&str],
        description: &str,
        required: bool,
        advanced: bool,
    ) -> Result<()> {
        self.register(
            name,
            ParameterType::OutputFileList,
            argument,
            ParamValue::StringList(default_value.iter().map(|s| (*s).to_owned()).collect()),
            description,
            required,
            advanced,
        )
    }

    /// Register every entry of a parameter tree as a command-line parameter, as
    /// `registerFullParam_` (`TOPPBase.cpp:1024-1045`).
    ///
    /// Each entry is described by
    /// [`ParameterInformation::from_param_entry`] under its section-qualified
    /// name, and every section that holds an entry is recorded, with its
    /// description, in [`topp_subsections`](Self::topp_subsections); the first
    /// description recorded for a section wins, as in the source map.
    ///
    /// # Errors
    ///
    /// As [`register_string_option`](Self::register_string_option), for each
    /// entry; entries registered before a failing one stay registered, as in
    /// the source.
    pub fn register_full_param(&mut self, param: &Param) -> Result<()> {
        for item in param.iter()? {
            if let Some((section, _)) = item.key.rsplit_once(':') {
                if !self
                    .topp_subsections
                    .iter()
                    .any(|(name, _)| name == section)
                {
                    let description = param.section_description(section)?.to_owned();
                    self.topp_subsections
                        .push((section.to_owned(), description));
                }
            }
            self.push(ParameterInformation::from_param_entry(
                item.entry, &item.key,
            ))?;
        }
        Ok(())
    }

    /// Restrict a string parameter to a fixed set, as `setValidStrings_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an unregistered or non-string
    /// parameter, or when a value contains a comma, which the INI format uses
    /// as the separator.
    pub fn set_valid_strings(&mut self, name: &str, strings: &[&str]) -> Result<()> {
        let entry = self.find_mut(name)?;
        if !matches!(
            entry.kind,
            ParameterType::String | ParameterType::StringList
        ) {
            return Err(bad(format!(
                "valid strings apply to string parameters only, not '{name}'"
            )));
        }
        if strings.iter().any(|s| s.contains(',')) {
            return Err(bad("valid strings must not contain a comma"));
        }
        entry.valid_strings = strings.iter().map(|s| (*s).to_owned()).collect();
        Ok(())
    }
    /// Restrict a file parameter to a set of extensions, as `setValidFormats_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an unregistered or non-file parameter.
    pub fn set_valid_formats(&mut self, name: &str, formats: &[&str]) -> Result<()> {
        let entry = self.find_mut(name)?;
        if !entry.kind.is_input_path() && !entry.kind.is_output_path() {
            return Err(bad(format!(
                "valid formats apply to file parameters only, not '{name}'"
            )));
        }
        entry.valid_formats = formats.iter().map(|s| (*s).to_owned()).collect();
        Ok(())
    }
    /// Set the inclusive lower bound of an integer parameter, as `setMinInt_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an unregistered parameter.
    pub fn set_min_int(&mut self, name: &str, value: i32) -> Result<()> {
        self.find_mut(name)?.min_int = Some(value);
        Ok(())
    }
    /// Set the inclusive upper bound of an integer parameter, as `setMaxInt_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an unregistered parameter.
    pub fn set_max_int(&mut self, name: &str, value: i32) -> Result<()> {
        self.find_mut(name)?.max_int = Some(value);
        Ok(())
    }
    /// Set the inclusive lower bound of a floating-point parameter, as
    /// `setMinFloat_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an unregistered parameter.
    pub fn set_min_float(&mut self, name: &str, value: f64) -> Result<()> {
        self.find_mut(name)?.min_float = Some(value);
        Ok(())
    }
    /// Set the inclusive upper bound of a floating-point parameter, as
    /// `setMaxFloat_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an unregistered parameter.
    pub fn set_max_float(&mut self, name: &str, value: f64) -> Result<()> {
        self.find_mut(name)?.max_float = Some(value);
        Ok(())
    }
    /// Register a nested parameter section, as `registerSubsection_`.
    ///
    /// The section's values come from
    /// [`Tool::subsection_defaults`](super::Tool::subsection_defaults).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an empty or duplicate name.
    pub fn register_subsection(&mut self, name: &str, description: &str) -> Result<()> {
        if name.is_empty() {
            return Err(bad("subsection needs a name"));
        }
        if self.subsections.iter().any(|(n, _)| n == name) {
            return Err(bad(format!("subsection '{name}' is registered twice")));
        }
        self.subsections
            .push((name.to_owned(), description.to_owned()));
        Ok(())
    }
    /// Add a usage-text line, as `addText_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when
    /// [`MAX_PARAMETERS`](Self::MAX_PARAMETERS) would be exceeded.
    pub fn add_text(&mut self, text: &str) -> Result<()> {
        self.push(ParameterInformation::text(text))
    }
    /// Add a blank usage-text line, as `addEmptyLine_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when
    /// [`MAX_PARAMETERS`](Self::MAX_PARAMETERS) would be exceeded.
    pub fn add_empty_line(&mut self) -> Result<()> {
        self.push(ParameterInformation::newline())
    }

    /// Project the registered parameters onto a `Param` tree under
    /// `<tool_name>:1:`, holding defaults, descriptions, tags and restrictions.
    ///
    /// This is the parameter part of the source `getDefaultParameters_`
    /// (`TOPPBase.cpp:2097-2228`): the tree a strict update merges INI and
    /// command-line values into and `-write_ini` writes. As in the source,
    /// `ini`, `instance`, the help flags and the `write_*` requests are left
    /// out, flags default to `false` restricted to `true`/`false`, and the
    /// descriptions of [`topp_subsections`](Self::topp_subsections) are set.
    /// The tool version item, the section descriptions of the tool and instance
    /// nodes and the algorithm subsections are added by the caller.
    ///
    /// # Errors
    ///
    /// Propagates parameter-tree failures, for example an invalid key.
    pub fn to_param(&self, tool_name: &str) -> Result<Param> {
        let mut param = Param::new();
        for entry in &self.parameters {
            if entry.kind.is_layout() || NOT_IN_PARAMETER_TREE.contains(&entry.name.as_str()) {
                continue;
            }
            let key = format!("{tool_name}:1:{}", entry.name);
            let default_value = if entry.kind == ParameterType::Flag {
                ParamValue::String("false".into())
            } else {
                entry.default_value.clone()
            };
            param.set_value(&key, default_value, &entry.description, &[])?;
            if entry.required {
                param.add_tag(&key, "required")?;
            }
            if entry.advanced {
                param.add_tag(&key, "advanced")?;
            }
            if entry.kind.is_input_path() {
                param.add_tag(&key, "input file")?;
            }
            if entry.kind.is_output_path() {
                param.add_tag(&key, "output file")?;
            }
            for tag in &entry.tags {
                param.add_tag(&key, tag)?;
            }
            if entry.kind == ParameterType::Flag {
                param.set_valid_strings(&key, &["true".to_owned(), "false".to_owned()])?;
            } else if !entry.valid_strings.is_empty() {
                param.set_valid_strings(&key, &entry.valid_strings)?;
            }
            if let Some(v) = entry.min_int {
                param.set_min_int(&key, v)?;
            }
            if let Some(v) = entry.max_int {
                param.set_max_int(&key, v)?;
            }
            if let Some(v) = entry.min_float {
                param.set_min_float(&key, v)?;
            }
            if let Some(v) = entry.max_float {
                param.set_max_float(&key, v)?;
            }
        }
        for (section, description) in &self.topp_subsections {
            param.set_section_description(&format!("{tool_name}:1:{section}"), description)?;
        }
        // Algorithm subsection descriptions are applied by the caller after the
        // subsection's own defaults are inserted: a section description cannot
        // be set on a section that holds no entries yet.
        Ok(param)
    }
}
