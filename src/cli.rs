// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The TOPP command-line framework: registration, the run lifecycle and usage
//! text. This is the native form of the `OpenMS4-cli` package's `TOPPBase`.
//!
//! See `docs/TOPP_CLI_SUPPORT.md` for the supported source subset.

mod context;
mod parameter;
mod spec;
pub mod tools;
mod usage;

pub use context::{ToolContext, parse_range};
pub use parameter::{ExitCode, ParameterInformation, ParameterType};
pub use spec::ToolSpec;

use crate::format::paramxml;
use crate::param::{CommandLineOptions, Param, ParamValue};
use crate::system::file;
use crate::{Error, Result};
use std::io::Write;

/// A TOPP tool: a name, a one-line description, its parameters and its body.
///
/// The source requires a subclass overriding `registerOptionsAndFlags_` and
/// `main_`. A native tool implements this trait instead, so registration and
/// execution are plain functions with no inherited mutable state.
pub trait Tool {
    /// Executable name, which is also the INI section name.
    const NAME: &'static str;
    /// One-line description shown in usage output.
    const DESCRIPTION: &'static str;

    /// Register every parameter, as `registerOptionsAndFlags_`.
    fn register(spec: &mut ToolSpec) -> Result<()>;

    /// Run the tool, as `main_`.
    fn run(ctx: &ToolContext) -> Result<ExitCode>;

    /// Defaults for a registered subsection, as `getSubsectionDefaults_`.
    ///
    /// A tool that wraps an algorithm registers a subsection with
    /// [`ToolSpec::register_subsection`] and returns that algorithm's parameter
    /// tree here. The values are merged under the tool's own defaults, so an
    /// INI file or a command line can override them, and `-write_ini` emits
    /// them. Returning `None` leaves the subsection empty.
    fn subsection_defaults(_section: &str) -> Result<Option<Param>> {
        Ok(None)
    }
}

/// The parameters every TOPP tool registers, in source order.
fn register_common(spec: &mut ToolSpec) -> Result<()> {
    spec.add_empty_line()?;
    spec.add_text("Common TOPP options:")?;
    spec.register_string_option(
        "ini",
        "<file>",
        "",
        "Use the given TOPP INI file",
        false,
        false,
    )?;
    spec.register_string_option(
        "log",
        "<file>",
        "",
        "Name of log file (created only when specified)",
        false,
        true,
    )?;
    spec.register_int_option(
        "instance",
        "<n>",
        1,
        "Instance number for the TOPP INI file",
        false,
        true,
    )?;
    spec.register_int_option("debug", "<n>", 0, "Sets the debug level", false, true)?;
    spec.register_int_option(
        "threads",
        "<n>",
        1,
        "Sets the number of threads allowed to be used by the TOPP tool (0 = all available cores)",
        false,
        false,
    )?;
    spec.register_string_option(
        "write_ini",
        "<file>",
        "",
        "Writes the default configuration file",
        false,
        false,
    )?;
    spec.register_flag(
        "no_progress",
        "Disables progress logging to command line",
        true,
    )?;
    spec.register_flag("force", "Overrides tool-specific checks", true)?;
    spec.register_flag(
        "test",
        "Enables the test mode (needed for internal use only)",
        true,
    )?;
    spec.register_flag("-help", "Shows options", false)?;
    spec.register_flag("-helphelp", "Shows all options (including advanced)", false)?;
    Ok(())
}

/// Build the full registration for a tool: its own parameters then the common ones.
pub fn tool_spec<T: Tool>() -> Result<ToolSpec> {
    let mut spec = ToolSpec::new();
    T::register(&mut spec)?;
    register_common(&mut spec)?;
    Ok(spec)
}

/// The tool's own defaults with every registered subsection's algorithm
/// defaults merged beneath it.
fn defaults_with_subsections<T: Tool>(spec: &ToolSpec) -> Result<Param> {
    let mut defaults = spec.to_param(T::NAME)?;
    for (name, description) in spec.subsections() {
        if let Some(values) = T::subsection_defaults(name)? {
            if !values.is_empty() {
                let prefix = format!("{}:1:{name}:", T::NAME);
                defaults.insert(&prefix, &values)?;
                defaults.set_section_description(prefix.trim_end_matches(':'), description)?;
            }
        }
    }
    Ok(defaults)
}

/// The command-line token maps the parameter parser needs.
fn command_line_options(spec: &ToolSpec, prefix: &str) -> CommandLineOptions {
    let mut options = CommandLineOptions::default();
    for entry in spec.parameters() {
        if entry.kind.is_layout() {
            continue;
        }
        let key = format!("{prefix}{}", entry.name);
        let token = entry.token();
        if entry.kind.is_list() {
            options.multiple_arguments.insert(token, key);
        } else if entry.kind.takes_argument() {
            options.one_argument.insert(token, key);
        } else {
            options.no_argument.insert(token, key);
        }
    }
    options.misc = format!("{prefix}misc");
    options.unknown = format!("{prefix}unknown");
    options
}

/// Coerce a command-line string into the type the parameter was registered with.
/// The parser yields String and StringList; source converts on read, and a
/// non-numeric value for a numeric parameter is an ILLEGAL_PARAMETERS error.
fn coerce(entry: &ParameterInformation, value: &ParamValue) -> Result<ParamValue> {
    let as_strings = |value: &ParamValue| -> Vec<String> {
        match value {
            ParamValue::StringList(values) => values.clone(),
            ParamValue::String(text) if text.is_empty() => Vec::new(),
            ParamValue::String(text) => vec![text.clone()],
            other => other.to_text(false).map(|t| vec![t]).unwrap_or_default(),
        }
    };
    let text = match value {
        ParamValue::String(text) => text.clone(),
        other => other.to_text(false)?,
    };
    let number = |what: &str| -> Error {
        Error::InvalidValue(format!(
            "value '{text}' for parameter '{}' is not {what}",
            entry.name
        ))
    };
    Ok(match entry.kind {
        ParameterType::Int => {
            ParamValue::Integer(text.trim().parse().map_err(|_| number("an integer"))?)
        }
        ParameterType::Double => {
            ParamValue::Float(text.trim().parse().map_err(|_| number("a number"))?)
        }
        ParameterType::IntList => ParamValue::IntegerList(
            as_strings(value)
                .iter()
                .map(|s| s.trim().parse().map_err(|_| number("an integer list")))
                .collect::<Result<_>>()?,
        ),
        ParameterType::DoubleList => ParamValue::FloatList(
            as_strings(value)
                .iter()
                .map(|s| s.trim().parse().map_err(|_| number("a number list")))
                .collect::<Result<_>>()?,
        ),
        ParameterType::StringList
        | ParameterType::InputFileList
        | ParameterType::OutputFileList => ParamValue::StringList(as_strings(value)),
        ParameterType::Flag => ParamValue::String(text),
        _ => ParamValue::String(text),
    })
}

/// Outcome of preparing a run: either parameters to execute with, or a
/// terminal exit code because usage, a version or an INI file was requested.
enum Prepared {
    Run(Box<ToolContext>),
    Done(ExitCode),
}

fn prepare<T: Tool>(
    spec: &ToolSpec,
    arguments: &[String],
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<Prepared> {
    let prefix = format!("{}:1:", T::NAME);
    let mut defaults = defaults_with_subsections::<T>(spec)?;

    // 1. Command line. The parser yields strings; values are coerced below.
    let mut cmdline = Param::new();
    cmdline.parse_command_line_mapped(arguments, &command_line_options(spec, &prefix))?;

    let given =
        |name: &str| -> bool { cmdline.exists(&format!("{prefix}{name}")).unwrap_or(false) };

    // 2. Usage requests short-circuit before any validation, as in source.
    if given("-help") || given("-helphelp") {
        usage::print(out, T::NAME, T::DESCRIPTION, spec, given("-helphelp"))?;
        return Ok(Prepared::Done(ExitCode::ExecutionOk));
    }

    // 3. An INI file is merged under the command line and over the defaults.
    if given("ini") {
        let path = cmdline.value(&format!("{prefix}ini"))?.as_str()?.to_owned();
        if !file::exists(&path) {
            writeln!(err, "Error: INI file '{path}' does not exist.")?;
            return Ok(Prepared::Done(ExitCode::InputFileNotFound));
        }
        let ini = paramxml::load(&path)?;
        let section = format!("{}:1:", T::NAME);
        if ini.has_section(section.trim_end_matches(':'))? {
            let tool_values = ini.copy(&section, false)?;
            defaults.update(&tool_values, false)?;
        }
    }

    // 4. Command-line values win, coerced to their registered types.
    for entry in spec.parameters() {
        if entry.kind.is_layout() {
            continue;
        }
        let key = format!("{prefix}{}", entry.name);
        if let Ok(true) = cmdline.exists(&key) {
            let raw = cmdline.value(&key)?.clone();
            let value = coerce(entry, &raw)?;
            defaults.set_value(&key, value, &entry.description, &[])?;
        }
    }

    // 5. -write_ini emits the resolved tree and stops, as in source.
    if given("write_ini") {
        let path = cmdline
            .value(&format!("{prefix}write_ini"))?
            .as_str()?
            .to_owned();
        let mut written = defaults.checked_clone()?;
        for name in ["write_ini", "ini", "-help", "-helphelp"] {
            let _ = written.remove(&format!("{prefix}{name}"));
        }
        paramxml::store(&path, &written)?;
        return Ok(Prepared::Done(ExitCode::ExecutionOk));
    }

    // 6. Validate required values, restrictions and file access.
    if let Some(code) = validate(spec, &defaults, &prefix, err)? {
        return Ok(Prepared::Done(code));
    }

    fn read_int(param: &Param, key: &str, fallback: i64) -> i64 {
        match param.value(key) {
            Ok(ParamValue::Integer(value)) => *value,
            _ => fallback,
        }
    }
    fn read_flag(param: &Param, key: &str) -> bool {
        matches!(param.value(key), Ok(ParamValue::String(text)) if text == "true")
    }
    let debug = read_int(&defaults, &format!("{prefix}debug"), 0);
    let threads = read_int(&defaults, &format!("{prefix}threads"), 1);
    let test = read_flag(&defaults, &format!("{prefix}test"));
    let no_progress = read_flag(&defaults, &format!("{prefix}no_progress"));
    let force = read_flag(&defaults, &format!("{prefix}force"));
    Ok(Prepared::Run(Box::new(ToolContext::new(
        defaults,
        prefix,
        debug,
        test,
        no_progress,
        force,
        threads,
    ))))
}

/// Required-value, restriction and file-access checks. Returns the terminal
/// exit code when the run cannot proceed.
fn validate(
    spec: &ToolSpec,
    param: &Param,
    prefix: &str,
    err: &mut dyn Write,
) -> Result<Option<ExitCode>> {
    for entry in spec.parameters() {
        if entry.kind.is_layout() {
            continue;
        }
        let key = format!("{prefix}{}", entry.name);
        let Ok(value) = param.value(&key) else {
            continue;
        };
        let empty = match value {
            ParamValue::String(text) => text.is_empty(),
            ParamValue::StringList(values) => values.is_empty(),
            ParamValue::Empty => true,
            _ => false,
        };
        if entry.required && empty {
            writeln!(
                err,
                "Error: Required parameter '{}' was not given.",
                entry.name
            )?;
            return Ok(Some(ExitCode::MissingParameters));
        }
        if empty {
            continue;
        }
        // Nested rather than a let-chain: let-chains need Rust 1.88 and the
        // crate's minimum is 1.85.
        if !entry.valid_strings.is_empty() {
            if let ParamValue::String(text) = value {
                if !entry.valid_strings.contains(text) {
                    writeln!(
                        err,
                        "Error: Invalid value '{text}' for parameter '{}'. Valid values are: {}.",
                        entry.name,
                        entry.valid_strings.join(", ")
                    )?;
                    return Ok(Some(ExitCode::IllegalParameters));
                }
            }
        }
        if let ParamValue::Integer(v) = value {
            if entry.min_int.is_some_and(|m| *v < i64::from(m))
                || entry.max_int.is_some_and(|m| *v > i64::from(m))
            {
                writeln!(
                    err,
                    "Error: Value {v} for parameter '{}' is out of range.",
                    entry.name
                )?;
                return Ok(Some(ExitCode::IllegalParameters));
            }
        }
        if let ParamValue::Float(v) = value {
            if entry.min_float.is_some_and(|m| *v < m) || entry.max_float.is_some_and(|m| *v > m) {
                writeln!(
                    err,
                    "Error: Value {v} for parameter '{}' is out of range.",
                    entry.name
                )?;
                return Ok(Some(ExitCode::IllegalParameters));
            }
        }
        let paths: Vec<&str> = match value {
            ParamValue::String(text) => vec![text.as_str()],
            ParamValue::StringList(values) => values.iter().map(String::as_str).collect(),
            _ => Vec::new(),
        };
        if entry.kind.is_input_path() {
            for path in &paths {
                if !file::exists(path) {
                    writeln!(
                        err,
                        "Error: Input file '{path}' for parameter '{}' does not exist.",
                        entry.name
                    )?;
                    return Ok(Some(ExitCode::InputFileNotFound));
                }
                if !context::extension_allowed(path, &entry.valid_formats) {
                    writeln!(
                        err,
                        "Error: Input file '{path}' for parameter '{}' has an unsupported format. Expected one of: {}.",
                        entry.name,
                        entry.valid_formats.join(", ")
                    )?;
                    return Ok(Some(ExitCode::IllegalParameters));
                }
            }
        }
        if matches!(
            entry.kind,
            ParameterType::OutputFile | ParameterType::OutputFileList
        ) {
            for path in &paths {
                if !context::extension_allowed(path, &entry.valid_formats) {
                    writeln!(
                        err,
                        "Error: Output file '{path}' for parameter '{}' has an unsupported format. Expected one of: {}.",
                        entry.name,
                        entry.valid_formats.join(", ")
                    )?;
                    return Ok(Some(ExitCode::IllegalParameters));
                }
            }
        }
    }
    Ok(None)
}

/// Map a recoverable error to the source exit code that reports it.
fn exit_code_for(error: &Error) -> ExitCode {
    match error {
        Error::Parse { .. } => ExitCode::ParseError,
        Error::Io(e) if e.kind() == std::io::ErrorKind::NotFound => ExitCode::InputFileNotFound,
        Error::Io(e) if e.kind() == std::io::ErrorKind::PermissionDenied => {
            ExitCode::CannotWriteOutputFile
        }
        Error::Io(_) => ExitCode::UnknownError,
        Error::Unsupported(_) => ExitCode::IncompatibleInputData,
        Error::UnsortedData => ExitCode::IncompatibleInputData,
        Error::InvalidValue(_) | Error::InvalidRange(_) => ExitCode::IllegalParameters,
        Error::MissingInformation(_) => ExitCode::MissingParameters,
    }
}

/// Run a tool against explicit arguments and streams. `arguments[0]` is the
/// executable name, as in `main(argc, argv)`.
pub fn run_with<T: Tool>(
    arguments: &[String],
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> ExitCode {
    let spec = match tool_spec::<T>() {
        Ok(spec) => spec,
        Err(error) => {
            let _ = writeln!(
                err,
                "Error: {} could not register its parameters: {error}",
                T::NAME
            );
            return ExitCode::InternalError;
        }
    };
    let prepared = match prepare::<T>(&spec, arguments, out, err) {
        Ok(prepared) => prepared,
        Err(error) => {
            let _ = writeln!(err, "Error: {error}");
            return exit_code_for(&error);
        }
    };
    let ctx = match prepared {
        Prepared::Done(code) => return code,
        Prepared::Run(ctx) => ctx,
    };
    match T::run(&ctx) {
        Ok(code) => code,
        Err(error) => {
            let _ = writeln!(err, "Error: {error}");
            exit_code_for(&error)
        }
    }
}

/// Run a tool against the process arguments and standard streams, returning the
/// status the executable should exit with.
pub fn run<T: Tool>() -> ExitCode {
    let arguments: Vec<String> = std::env::args().collect();
    let mut out = std::io::stdout();
    let mut err = std::io::stderr();
    run_with::<T>(&arguments, &mut out, &mut err)
}
