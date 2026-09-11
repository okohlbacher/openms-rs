// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Command-line parameter records and tool exit codes.

use crate::param::ParamValue;

/// Source `TOPPBase::ExitCodes`, in declaration order. The discriminant is the
/// process exit status, so `EXECUTION_OK` is zero.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[repr(i32)]
pub enum ExitCode {
    ExecutionOk = 0,
    InputFileNotFound = 1,
    InputFileNotReadable = 2,
    InputFileCorrupt = 3,
    InputFileEmpty = 4,
    CannotWriteOutputFile = 5,
    IllegalParameters = 6,
    MissingParameters = 7,
    UnknownError = 8,
    ExternalProgramError = 9,
    ParseError = 10,
    IncompatibleInputData = 11,
    InternalError = 12,
    UnexpectedResult = 13,
    ExternalProgramNotFound = 14,
}

impl ExitCode {
    pub fn as_i32(self) -> i32 {
        self as i32
    }
    /// Source diagnostic name, as printed before the process exits.
    pub fn name(self) -> &'static str {
        match self {
            Self::ExecutionOk => "EXECUTION_OK",
            Self::InputFileNotFound => "INPUT_FILE_NOT_FOUND",
            Self::InputFileNotReadable => "INPUT_FILE_NOT_READABLE",
            Self::InputFileCorrupt => "INPUT_FILE_CORRUPT",
            Self::InputFileEmpty => "INPUT_FILE_EMPTY",
            Self::CannotWriteOutputFile => "CANNOT_WRITE_OUTPUT_FILE",
            Self::IllegalParameters => "ILLEGAL_PARAMETERS",
            Self::MissingParameters => "MISSING_PARAMETERS",
            Self::UnknownError => "UNKNOWN_ERROR",
            Self::ExternalProgramError => "EXTERNAL_PROGRAM_ERROR",
            Self::ParseError => "PARSE_ERROR",
            Self::IncompatibleInputData => "INCOMPATIBLE_INPUT_DATA",
            Self::InternalError => "INTERNAL_ERROR",
            Self::UnexpectedResult => "UNEXPECTED_RESULT",
            Self::ExternalProgramNotFound => "EXTERNAL_PROGRAM_NOTFOUND",
        }
    }
}

/// Source `ParameterInformation::ParameterTypes`. `Text` and `Newline` are
/// layout-only entries produced by `add_text` and `add_empty_line`; they carry
/// no value and never reach the parameter tree.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ParameterType {
    None,
    String,
    InputFile,
    OutputFile,
    OutputPrefix,
    OutputDir,
    Double,
    Int,
    StringList,
    IntList,
    DoubleList,
    InputFileList,
    OutputFileList,
    Flag,
    Text,
    Newline,
}

impl ParameterType {
    /// Whether the option consumes one or more following command-line tokens.
    pub fn takes_argument(self) -> bool {
        !matches!(self, Self::Flag | Self::Text | Self::Newline | Self::None)
    }
    /// Whether the option consumes every following non-option token.
    pub fn is_list(self) -> bool {
        matches!(
            self,
            Self::StringList
                | Self::IntList
                | Self::DoubleList
                | Self::InputFileList
                | Self::OutputFileList
        )
    }
    /// Whether the value names a file or directory that is read.
    pub fn is_input_path(self) -> bool {
        matches!(self, Self::InputFile | Self::InputFileList)
    }
    /// Whether the value names a path that is written.
    pub fn is_output_path(self) -> bool {
        matches!(
            self,
            Self::OutputFile | Self::OutputFileList | Self::OutputPrefix | Self::OutputDir
        )
    }
    /// Layout-only entries are shown in usage text but hold no value.
    pub fn is_layout(self) -> bool {
        matches!(self, Self::Text | Self::Newline)
    }
}

/// Source `ParameterInformation`: everything known about one command-line
/// parameter. Restrictions are unset by default and are applied to the
/// parameter tree only when a caller sets them.
#[derive(Clone, Debug)]
pub struct ParameterInformation {
    pub name: String,
    pub kind: ParameterType,
    pub default_value: ParamValue,
    pub description: String,
    /// Placeholder shown in usage text, for example `<file>`.
    pub argument: String,
    pub required: bool,
    pub advanced: bool,
    pub tags: Vec<String>,
    pub valid_strings: Vec<String>,
    pub min_int: Option<i32>,
    pub max_int: Option<i32>,
    pub min_float: Option<f64>,
    pub max_float: Option<f64>,
    /// Accepted file extensions, from `set_valid_formats`.
    pub valid_formats: Vec<String>,
}

impl ParameterInformation {
    pub fn new(
        name: impl Into<String>,
        kind: ParameterType,
        argument: impl Into<String>,
        default_value: ParamValue,
        description: impl Into<String>,
        required: bool,
        advanced: bool,
    ) -> Self {
        Self {
            name: name.into(),
            kind,
            default_value,
            description: description.into(),
            argument: argument.into(),
            required,
            advanced,
            tags: Vec::new(),
            valid_strings: Vec::new(),
            min_int: None,
            max_int: None,
            min_float: None,
            max_float: None,
            valid_formats: Vec::new(),
        }
    }
    /// A layout-only text line in the usage output.
    pub fn text(text: impl Into<String>) -> Self {
        Self::new(
            "",
            ParameterType::Text,
            "",
            ParamValue::Empty,
            text,
            false,
            false,
        )
    }
    /// A layout-only blank line in the usage output.
    pub fn newline() -> Self {
        Self::new(
            "",
            ParameterType::Newline,
            "",
            ParamValue::Empty,
            "",
            false,
            false,
        )
    }
    /// The token a user types, which is the name with a single leading dash.
    /// Source registers the two help flags as `-help` and `-helphelp`, so their
    /// tokens gain a second dash here exactly as they do in C++.
    pub fn token(&self) -> String {
        format!("-{}", self.name)
    }
}
