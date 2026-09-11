// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Tool parameter registration, the native form of the source `register*_` calls.

use super::parameter::{ParameterInformation, ParameterType};
use crate::param::{Param, ParamValue};
use crate::{Error, Result};

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

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
}

impl ToolSpec {
    /// Maximum registered parameters for one tool.
    pub const MAX_PARAMETERS: usize = 4096;

    pub fn new() -> Self {
        Self::default()
    }
    pub fn parameters(&self) -> &[ParameterInformation] {
        &self.parameters
    }
    pub fn subsections(&self) -> &[(String, String)] {
        &self.subsections
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

    /// Restrict a string parameter to a fixed set, as `setValidStrings_`.
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
    pub fn set_min_int(&mut self, name: &str, value: i32) -> Result<()> {
        self.find_mut(name)?.min_int = Some(value);
        Ok(())
    }
    pub fn set_max_int(&mut self, name: &str, value: i32) -> Result<()> {
        self.find_mut(name)?.max_int = Some(value);
        Ok(())
    }
    pub fn set_min_float(&mut self, name: &str, value: f64) -> Result<()> {
        self.find_mut(name)?.min_float = Some(value);
        Ok(())
    }
    pub fn set_max_float(&mut self, name: &str, value: f64) -> Result<()> {
        self.find_mut(name)?.max_float = Some(value);
        Ok(())
    }
    /// Register a nested parameter section, as `registerSubsection_`.
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
    pub fn add_text(&mut self, text: &str) -> Result<()> {
        self.push(ParameterInformation::text(text))
    }
    pub fn add_empty_line(&mut self) -> Result<()> {
        self.push(ParameterInformation::newline())
    }

    /// Project the registered parameters onto a `Param` tree holding defaults,
    /// descriptions, tags and restrictions. This is the tree written by
    /// `-write_ini` and the one command-line and INI values are merged into.
    pub fn to_param(&self, tool_name: &str) -> Result<Param> {
        let mut param = Param::new();
        for entry in &self.parameters {
            if entry.kind.is_layout() {
                continue;
            }
            let key = format!("{tool_name}:1:{}", entry.name);
            param.set_value(&key, entry.default_value.clone(), &entry.description, &[])?;
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
            if !entry.valid_strings.is_empty() {
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
        for (name, description) in &self.subsections {
            param.set_section_description(&format!("{tool_name}:1:{name}"), description)?;
        }
        Ok(param)
    }
}
