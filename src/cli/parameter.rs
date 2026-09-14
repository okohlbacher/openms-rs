// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Command-line parameter records and tool exit codes, the native form of the
//! source `APPLICATIONS/ParameterInformation.h` and `TOPPBase::ExitCodes`.
//!
//! See `docs/TOPP_CLI_SUPPORT.md` for the supported source subset.

use crate::param::{ParamEntry, ParamValue};

/// Source `TOPPBase::ExitCodes`, in declaration order.
///
/// The discriminant is the process exit status, so `EXECUTION_OK` is zero.
/// Which failure maps to which code depends on the lifecycle phase it occurs
/// in, as in `TOPPBase::main`; `docs/TOPP_CLI_SUPPORT.md` tabulates the mapping.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[repr(i32)]
pub enum ExitCode {
    /// The tool ran to completion (`EXECUTION_OK`).
    ExecutionOk = 0,
    /// An input file does not exist (`INPUT_FILE_NOT_FOUND`).
    InputFileNotFound = 1,
    /// An input file exists but cannot be read (`INPUT_FILE_NOT_READABLE`).
    InputFileNotReadable = 2,
    /// A file could not be parsed while the tool ran (`INPUT_FILE_CORRUPT`),
    /// including a malformed `-ini` file.
    InputFileCorrupt = 3,
    /// An input file holds no bytes or no usable records (`INPUT_FILE_EMPTY`).
    InputFileEmpty = 4,
    /// An output file cannot be created (`CANNOT_WRITE_OUTPUT_FILE`).
    CannotWriteOutputFile = 5,
    /// The command line or the parameters are invalid (`ILLEGAL_PARAMETERS`).
    IllegalParameters = 6,
    /// A required parameter was not given or is empty (`MISSING_PARAMETERS`).
    MissingParameters = 7,
    /// Any other error raised while the tool runs (`UNKNOWN_ERROR`).
    UnknownError = 8,
    /// An external program failed (`EXTERNAL_PROGRAM_ERROR`).
    ExternalProgramError = 9,
    /// A failure a tool reports as a parse error itself (`PARSE_ERROR`), for
    /// example an input of unknown type in `FileInfo`.
    ParseError = 10,
    /// Valid input the tool cannot process (`INCOMPATIBLE_INPUT_DATA`).
    IncompatibleInputData = 11,
    /// A failure of the framework rather than of the input (`INTERNAL_ERROR`).
    InternalError = 12,
    /// The tool produced an unexpected result (`UNEXPECTED_RESULT`).
    UnexpectedResult = 13,
    /// An external program could not be found (`EXTERNAL_PROGRAM_NOTFOUND`).
    ExternalProgramNotFound = 14,
}

impl ExitCode {
    /// The process exit status, which is the source enumerator value.
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

/// Source `ParameterInformation::ParameterTypes`.
///
/// `Text` and `Newline` are layout-only entries produced by
/// [`ToolSpec::add_text`](super::ToolSpec::add_text) and
/// [`ToolSpec::add_empty_line`](super::ToolSpec::add_empty_line); they carry no
/// value and never reach the parameter tree.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ParameterType {
    /// No value type (`NONE`); the source's type for an empty parameter value.
    None,
    /// Free text (`STRING`).
    String,
    /// A path that is read (`INPUT_FILE`).
    InputFile,
    /// A path that is written (`OUTPUT_FILE`).
    OutputFile,
    /// A prefix of paths that are written (`OUTPUT_PREFIX`).
    OutputPrefix,
    /// A directory that is written into (`OUTPUT_DIR`).
    OutputDir,
    /// A floating-point number (`DOUBLE`).
    Double,
    /// A 32-bit integer (`INT`).
    Int,
    /// A list of strings (`STRINGLIST`).
    StringList,
    /// A list of integers (`INTLIST`).
    IntList,
    /// A list of floating-point numbers (`DOUBLELIST`).
    DoubleList,
    /// A list of paths that are read (`INPUT_FILE_LIST`).
    InputFileList,
    /// A list of paths that are written (`OUTPUT_FILE_LIST`).
    OutputFileList,
    /// A switch that is `true` when given and `false` otherwise (`FLAG`).
    Flag,
    /// A usage-text line (`TEXT`).
    Text,
    /// A blank usage-text line (`NEWLINE`).
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
/// parameter.
///
/// Restrictions are unset by default and are applied to the parameter tree only
/// when a caller sets them. The source stores accepted file formats in
/// `valid_strings`; this port keeps them apart in [`valid_formats`](Self::valid_formats)
/// for parameters registered through [`ToolSpec`](super::ToolSpec), and keeps
/// the source layout for parameters converted with
/// [`from_param_entry`](Self::from_param_entry).
#[derive(Clone, Debug)]
pub struct ParameterInformation {
    /// Name without the leading dash; subsection parameters carry their section
    /// path, as in `algorithm:peakcount`.
    pub name: String,
    /// Value type, which also decides how many command-line tokens are consumed.
    pub kind: ParameterType,
    /// Default value, stored in the type the parameter was registered with.
    pub default_value: ParamValue,
    /// Description shown in usage text and written to INI files.
    pub description: String,
    /// Placeholder shown in usage text, for example `<file>`.
    pub argument: String,
    /// Whether a non-empty value must be supplied.
    pub required: bool,
    /// Whether the parameter is hidden from `--help` and listed by `--helphelp`.
    pub advanced: bool,
    /// Free-form tags, for example `skipexists` or `is_executable`.
    pub tags: Vec<String>,
    /// Accepted values of a string or string-list parameter.
    pub valid_strings: Vec<String>,
    /// Inclusive lower bound of an integer or integer-list parameter.
    pub min_int: Option<i32>,
    /// Inclusive upper bound of an integer or integer-list parameter.
    pub max_int: Option<i32>,
    /// Inclusive lower bound of a floating-point parameter.
    pub min_float: Option<f64>,
    /// Inclusive upper bound of a floating-point parameter.
    pub max_float: Option<f64>,
    /// Accepted file extensions, from `set_valid_formats`.
    pub valid_formats: Vec<String>,
}

impl ParameterInformation {
    /// A parameter with no restrictions and no tags, as the source constructor.
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
    ///
    /// Source registers the two help flags as `-help` and `-helphelp`, so their
    /// tokens gain a second dash here exactly as they do in C++.
    pub fn token(&self) -> String {
        format!("-{}", self.name)
    }

    /// The formats a file parameter accepts, in the order the source keeps them.
    ///
    /// The source stores a file parameter's formats in `valid_strings`. A
    /// parameter registered through [`ToolSpec`](super::ToolSpec) keeps them in
    /// [`valid_formats`](Self::valid_formats), and one converted with
    /// [`from_param_entry`](Self::from_param_entry) in
    /// [`valid_strings`](Self::valid_strings), so both are yielded, registered
    /// formats first. A parameter that names no file path accepts no formats.
    pub fn accepted_formats(&self) -> impl Iterator<Item = &str> + '_ {
        let path = self.kind.is_input_path() || self.kind.is_output_path();
        self.valid_formats
            .iter()
            .chain(self.valid_strings.iter())
            .filter(move |_| path)
            .map(String::as_str)
    }

    /// Describe one entry of a parameter tree as a command-line parameter.
    ///
    /// Source `TOPPBase::paramEntryToParameterInformation_` together with
    /// `getParamArgument_` (`TOPPBase.cpp:896-1008`). This is how subsection
    /// parameters become addressable on the command line and how
    /// `registerFullParam_` registers a whole tree. `full_name` is the
    /// section-qualified name; an empty `full_name` falls back to the entry's
    /// own leaf name.
    ///
    /// The type follows the value: a string whose value is `false` and whose
    /// valid strings are exactly `true`, `false` becomes a [`ParameterType::Flag`],
    /// and the `input file`, `output file`, `output prefix` and `output dir`
    /// tags turn strings and string lists into the matching path types.
    /// Restrictions are copied; the tree's unset sentinels (`-i32::MAX`,
    /// `i32::MAX`, `-f64::MAX`, `f64::MAX`) become `None`. Valid strings of a
    /// path entry stay in [`valid_strings`](Self::valid_strings), where the
    /// source keeps them.
    pub fn from_param_entry(entry: &ParamEntry, full_name: &str) -> Self {
        let name = if full_name.is_empty() {
            entry.name.as_str()
        } else {
            full_name
        };
        let advanced = entry.tags.contains("advanced");
        let is_flag = matches!(&entry.value, ParamValue::String(value) if value == "false")
            && entry.valid_strings.len() == 2
            && entry.valid_strings[0] == "true"
            && entry.valid_strings[1] == "false";
        if is_flag {
            return Self::new(
                name,
                ParameterType::Flag,
                "",
                ParamValue::String("false".into()),
                entry.description.as_str(),
                false,
                advanced,
            );
        }
        let tagged = |tag: &str| entry.tags.contains(tag);
        let (kind, argument) = match &entry.value {
            ParamValue::String(_) => (
                if tagged("input file") {
                    ParameterType::InputFile
                } else if tagged("output file") {
                    ParameterType::OutputFile
                } else if tagged("output prefix") {
                    ParameterType::OutputPrefix
                } else if tagged("output dir") {
                    ParameterType::OutputDir
                } else {
                    ParameterType::String
                },
                if entry.valid_strings.is_empty() {
                    "<text>"
                } else {
                    "<choice>"
                },
            ),
            ParamValue::Integer(_) => (ParameterType::Int, "<number>"),
            ParamValue::Float(_) => (ParameterType::Double, "<value>"),
            ParamValue::StringList(_) => (
                if tagged("input file") {
                    ParameterType::InputFileList
                } else if tagged("output file") {
                    ParameterType::OutputFileList
                } else {
                    ParameterType::StringList
                },
                "<list>",
            ),
            ParamValue::IntegerList(_) => (ParameterType::IntList, "<numbers>"),
            ParamValue::FloatList(_) => (ParameterType::DoubleList, "<values>"),
            ParamValue::Empty => (ParameterType::None, ""),
        };
        let mut information = Self::new(
            name,
            kind,
            argument,
            entry.value.clone(),
            entry.description.as_str(),
            tagged("required"),
            advanced,
        );
        information.valid_strings = entry.valid_strings.clone();
        information.min_int = (entry.min_int != -i32::MAX).then_some(entry.min_int);
        information.max_int = (entry.max_int != i32::MAX).then_some(entry.max_int);
        information.min_float = (entry.min_float != -f64::MAX).then_some(entry.min_float);
        information.max_float = (entry.max_float != f64::MAX).then_some(entry.max_float);
        information
    }
}
