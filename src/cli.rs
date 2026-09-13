// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The TOPP command-line framework: registration, the run lifecycle, usage
//! text and processing annotations.
//!
//! This is the native form of the `OpenMS4-cli` package's
//! `APPLICATIONS/TOPPBase.h` and `APPLICATIONS/ParameterInformation.h`.
//! [`run_with`](crate::cli::run_with) follows `TOPPBase::main`
//! (`TOPPBase.cpp:141-514`) phase by phase:
//!
//! 1. registration of the tool's own parameters and the common ones;
//! 2. the command-line parse closure, which also accepts every subsection
//!    parameter as `-section:name` and rejects unknown options and trailing
//!    text;
//! 3. `-write_ini`, and the refused tool-description writers;
//! 4. the INI merge of the instance, `common:<tool>:` and `common:` sections
//!    and a strict update of the defaults, which rejects unknown parameters and
//!    invalid values;
//! 5. parameter and file validation, then the tool body.
//!
//! A failure maps to the exit code of the phase it occurs in, as the source's
//! two catch blocks do. See `docs/TOPP_CLI_SUPPORT.md` for the supported source
//! subset and the exit-code table.

mod context;
mod parameter;
mod processing;
mod spec;
pub mod tools;
mod usage;

pub use context::{
    TEST_MODE_UNIQUE_ID_SEED, ToolContext, input_file_readable, output_file_writable, parse_range,
};
pub use parameter::{ExitCode, ParameterInformation, ParameterType};
pub use processing::{
    AddDataProcessing, TEST_MODE_COMPLETION_TIME, TEST_MODE_PARAMETER_KEY,
    TEST_MODE_PARAMETER_VALUE, TEST_MODE_VERSION,
};
pub use spec::ToolSpec;

use crate::format::paramxml;
use crate::param::{Param, ParamEntry, ParamUpdateOptions, ParamValue};
use crate::system::file;
use crate::{Error, Result};
use std::collections::{BTreeMap, VecDeque};
use std::io::Write;

/// The product version every ported TOPP tool reports unless it overrides
/// [`Tool::VERSION`].
///
/// The source reads it from the installed tool manifest through
/// `ToolHandler::getToolVersion` (`TOPPBase.cpp:144-150`); the product SDK's
/// `share/openms4/tools/topp.tools.tsv` lists 1.0.0 for every tool, and the
/// oracle binaries print it in their usage text.
pub const TOPP_PRODUCT_VERSION: &str = "1.0.0";

/// Most command-line tokens [`run_with`] accepts, including the executable name.
///
/// The source has no limit beyond the operating system's; this bound exists for
/// callers that pass arguments in-process. A longer command line is refused
/// with [`ExitCode::IllegalParameters`] before anything is parsed.
pub const MAX_ARGUMENTS: usize = 1 << 20;

/// Most bytes, summed over all command-line tokens, that [`run_with`] accepts.
pub const MAX_ARGUMENT_BYTES: usize = 64 << 20;

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
    /// Product version, written to INI files as `<tool>:version` and recorded
    /// in processing annotations outside `-test`.
    const VERSION: &'static str = TOPP_PRODUCT_VERSION;

    /// Register every parameter, as `registerOptionsAndFlags_`.
    ///
    /// # Errors
    ///
    /// A registration error ends the run with [`ExitCode::IllegalParameters`],
    /// as the source's initialisation catch does.
    fn register(spec: &mut ToolSpec) -> Result<()>;

    /// Run the tool, as `main_`.
    ///
    /// # Errors
    ///
    /// An error is reported on the error stream and mapped to the exit code of
    /// the source exception it corresponds to; see [`run_with`].
    fn run(ctx: &ToolContext) -> Result<ExitCode>;

    /// Run the tool with explicit output and error streams.
    ///
    /// [`run_with`] calls this; the default forwards to [`run`](Self::run), so
    /// a tool that writes nothing but files implements `run` only. A tool that
    /// reports on standard output — the source's `OPENMS_LOG_INFO` — or writes
    /// its own diagnostics overrides this, writes to `out` and `err`, and
    /// returns an exit code such as [`ExitCode::InputFileEmpty`] after writing
    /// an `Error: …` line, so that `run_with` callers capture everything. Such
    /// a tool implements `run` by forwarding to `run_io` with the process
    /// streams.
    ///
    /// # Errors
    ///
    /// As [`run`](Self::run).
    fn run_io(ctx: &ToolContext, _out: &mut dyn Write, _err: &mut dyn Write) -> Result<ExitCode> {
        Self::run(ctx)
    }

    /// Defaults for a registered subsection, as `getSubsectionDefaults_`.
    ///
    /// A tool that wraps an algorithm registers a subsection with
    /// [`ToolSpec::register_subsection`] and returns that algorithm's parameter
    /// tree here. The values are merged under the tool's own defaults, so an
    /// INI file or a `-section:name` command-line option can override them,
    /// and `-write_ini` emits them. Returning `None` leaves the subsection empty.
    ///
    /// # Errors
    ///
    /// An error ends the run with [`ExitCode::IllegalParameters`].
    fn subsection_defaults(_section: &str) -> Result<Option<Param>> {
        Ok(None)
    }
}

/// The parameters every TOPP tool registers, in source order
/// (`TOPPBase.cpp:158-179`).
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
    spec.register_string_option(
        "write_ctd",
        "<out_dir>",
        "",
        "Writes the common tool description file(s) (Toolname(s).ctd) to <out_dir>",
        false,
        true,
    )?;
    spec.register_string_option(
        "write_nested_cwl",
        "<out_dir>",
        "",
        "Writes the Common Workflow Language file(s) (Toolname(s).cwl) to <out_dir>",
        false,
        true,
    )?;
    spec.register_string_option(
        "write_cwl",
        "<out_dir>",
        "",
        "Writes the Common Workflow Language file(s) (Toolname(s).cwl) to <out_dir>, but enforce a flat parameter hierarchy",
        false,
        true,
    )?;
    spec.register_string_option(
        "write_nested_json",
        "<out_dir>",
        "",
        "Writes the default configuration file",
        false,
        true,
    )?;
    spec.register_string_option(
        "write_json",
        "<out_dir>",
        "",
        "Writes the default configuration file, but compatible to the flat hierarchy",
        false,
        true,
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

/// Build the full registration for a tool: its own parameters then the common
/// ones.
///
/// # Errors
///
/// Propagates the tool's registration error, or a duplicate of a common name.
pub fn tool_spec<T: Tool>() -> Result<ToolSpec> {
    let mut spec = ToolSpec::new();
    T::register(&mut spec)?;
    register_common(&mut spec)?;
    Ok(spec)
}

/// The INI section of a tool, as `getToolPrefix`. The instance is always 1.
fn ini_location<T: Tool>() -> String {
    format!("{}:1:", T::NAME)
}

/// Every registered subsection's defaults, prefixed by the subsection name and
/// ordered by name, as the source `std::map` of subsections is
/// (`TOPPBase.cpp:2258-2274`).
fn subsection_defaults<T: Tool>(spec: &ToolSpec) -> Result<Param> {
    let mut sections: Vec<&(String, String)> = spec.subsections().iter().collect();
    sections.sort_by(|a, b| a.0.cmp(&b.0));
    let mut all = Param::new();
    for (name, description) in sections {
        if let Some(values) = T::subsection_defaults(name)? {
            if !values.is_empty() {
                all.insert(&format!("{name}:"), &values)?;
                all.set_section_description(name, description)?;
            }
        }
    }
    Ok(all)
}

/// Source `getDefaultParameters_` (`TOPPBase.cpp:2097-2256`) without the
/// per-user defaults file: the registered parameters, the tool version item,
/// the tool and instance section descriptions and the subsection defaults.
fn default_parameters<T: Tool>(spec: &ToolSpec, subsections: &Param) -> Result<Param> {
    let location = ini_location::<T>();
    let mut defaults = spec.to_param(T::NAME)?;
    defaults.set_value(
        &format!("{}:version", T::NAME),
        ParamValue::String(T::VERSION.to_owned()),
        "Version of the tool that generated this parameters file.",
        &["advanced".to_owned()],
    )?;
    defaults.set_section_description(T::NAME, T::DESCRIPTION)?;
    defaults.set_section_description(
        location.trim_end_matches(':'),
        &format!("Instance '1' section for '{}'", T::NAME),
    )?;
    if !subsections.is_empty() {
        defaults.insert(&location, subsections)?;
        for (name, description) in spec.subsections() {
            if subsections.has_section(name)? {
                defaults.set_section_description(&format!("{location}{name}"), description)?;
            }
        }
    }
    Ok(defaults)
}

/// Command-line values and leftovers, as `parseCommandLine_` returns them.
struct CommandLine {
    /// Registered parameters that were given, keyed by name without a prefix.
    values: Param,
    /// Option-like tokens that name no registered parameter, in the order the
    /// source collects them: last on the command line first.
    unknown: Vec<String>,
    /// Text tokens that no option consumed.
    misc: Vec<String>,
    /// Duplicate-parameter warnings.
    warnings: Vec<String>,
}

/// A command-line parse failure: the source exception name and its message.
struct ParseFailure {
    kind: &'static str,
    message: String,
}

/// Source option test: a dash then a letter, or two dashes then a letter.
fn is_option(argument: &str) -> bool {
    let bytes = argument.as_bytes();
    bytes.len() >= 2
        && bytes[0] == b'-'
        && (bytes[1].is_ascii_alphabetic()
            || (bytes[1] == b'-' && bytes.len() >= 3 && bytes[2].is_ascii_alphabetic()))
}

/// The value of one registered option, consuming its tokens from the front of
/// `queue` (`TOPPBase.cpp:2330-2414`).
fn option_value(
    definition: &ParameterInformation,
    argument: &str,
    queue: &mut VecDeque<&str>,
) -> std::result::Result<ParamValue, ParseFailure> {
    let conversion = |message: String| ParseFailure {
        kind: "ConversionError",
        message,
    };
    if definition.kind == ParameterType::Flag {
        if !queue.is_empty() {
            let trailing: Vec<&str> = queue.iter().copied().collect();
            return Err(ParseFailure {
                kind: "InvalidParameter",
                message: format!(
                    "Command line error: Trailing arguments after flag '{argument}': {}",
                    trailing.join(" ")
                ),
            });
        }
        return Ok(ParamValue::String("true".into()));
    }
    let value = match definition.kind {
        ParameterType::String
        | ParameterType::InputFile
        | ParameterType::OutputFile
        | ParameterType::OutputPrefix
        | ParameterType::OutputDir => {
            ParamValue::String(queue.front().map_or_else(String::new, |t| (*t).to_owned()))
        }
        ParameterType::Int => match queue.front() {
            Some(text) => {
                ParamValue::Integer(i64::from(context::to_int32(text).map_err(conversion)?))
            }
            None => ParamValue::Empty,
        },
        ParameterType::Double => match queue.front() {
            Some(text) => ParamValue::Float(context::to_double(text).map_err(conversion)?),
            None => ParamValue::Empty,
        },
        ParameterType::StringList
        | ParameterType::InputFileList
        | ParameterType::OutputFileList => {
            ParamValue::StringList(queue.drain(..).map(str::to_owned).collect())
        }
        ParameterType::IntList => ParamValue::IntegerList(
            queue
                .drain(..)
                .map(context::to_int32)
                .collect::<std::result::Result<_, _>>()
                .map_err(conversion)?,
        ),
        ParameterType::DoubleList => ParamValue::FloatList(
            queue
                .drain(..)
                .map(context::to_double)
                .collect::<std::result::Result<_, _>>()
                .map_err(conversion)?,
        ),
        ParameterType::Flag
        | ParameterType::None
        | ParameterType::Text
        | ParameterType::Newline => ParamValue::Empty,
    };
    queue.pop_front();
    Ok(value)
}

/// Source `parseCommandLine_` (`TOPPBase.cpp:2288-2470`).
///
/// Tokens are read from the last to the first, so an option finds its values
/// already queued. A flag followed by text, or a number that does not convert,
/// is a parse failure. When an option is given twice, the last occurrence wins
/// and a warning names the ignored one. `arguments[0]` is the executable name;
/// an empty slice parses as no arguments.
fn parse_command_line(
    arguments: &[String],
    definitions: &[ParameterInformation],
) -> Result<std::result::Result<CommandLine, ParseFailure>> {
    let mut by_token: BTreeMap<String, &ParameterInformation> = BTreeMap::new();
    for definition in definitions {
        if !definition.kind.is_layout() {
            by_token.insert(definition.token(), definition);
        }
    }
    let mut parsed = CommandLine {
        values: Param::new(),
        unknown: Vec::new(),
        misc: Vec::new(),
        warnings: Vec::new(),
    };
    let mut queue: VecDeque<&str> = VecDeque::new();
    for argument in arguments.iter().skip(1).rev() {
        if !is_option(argument) {
            queue.push_front(argument);
            continue;
        }
        match by_token.get(argument.as_str()) {
            Some(definition) => {
                let value = match option_value(definition, argument, &mut queue) {
                    Ok(value) => value,
                    Err(failure) => return Ok(Err(failure)),
                };
                if parsed.values.exists(&definition.name)? {
                    let kept = parsed
                        .values
                        .value(&definition.name)?
                        .to_text(true)
                        .unwrap_or_default();
                    parsed.warnings.push(format!(
                        "Warning: Duplicate parameter '{argument}' given. Using last occurrence with value '{kept}' (ignoring '{}').",
                        value.to_text(true).unwrap_or_default()
                    ));
                } else {
                    parsed.values.set_value(&definition.name, value, "", &[])?;
                }
            }
            None => parsed.unknown.push(argument.clone()),
        }
        let rest: Vec<String> = queue.drain(..).map(str::to_owned).collect();
        parsed.misc.splice(0..0, rest);
    }
    let rest: Vec<String> = queue.drain(..).map(str::to_owned).collect();
    parsed.misc.splice(0..0, rest);
    Ok(Ok(parsed))
}

/// Source rendering of a string list, as `ParamValue::toString`.
fn list_text(values: &[String]) -> String {
    format!("[{}]", values.join(", "))
}

/// A string-valued command-line or INI value.
fn text_value(param: &Param, key: &str) -> Result<String> {
    param.value(key)?.to_text(false)
}

fn print_usage<T: Tool>(stream: &mut dyn Write, spec: &ToolSpec, advanced: bool) -> Result<()> {
    usage::print(stream, T::NAME, T::DESCRIPTION, spec, advanced)?;
    Ok(())
}

/// Outcome of preparing a run: either parameters to execute with, or a
/// terminal exit code because usage, an INI file or a failure ended the run.
enum Prepared {
    Run(Box<ToolContext>),
    Done(ExitCode),
}

/// The lifecycle up to the tool body (`TOPPBase.cpp:157-408`).
fn prepare<T: Tool>(
    spec: &ToolSpec,
    arguments: &[String],
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<Prepared> {
    let location = ini_location::<T>();
    let subsections = subsection_defaults::<T>(spec)?;

    // Native preflight: bound the work before anything is parsed.
    let bytes = arguments
        .iter()
        .try_fold(0usize, |sum, argument| sum.checked_add(argument.len()));
    if arguments.len() > MAX_ARGUMENTS || !matches!(bytes, Some(b) if b <= MAX_ARGUMENT_BYTES) {
        writeln!(
            err,
            "Invalid parameter values (InvalidParameter): the command line exceeds {MAX_ARGUMENTS} arguments or {MAX_ARGUMENT_BYTES} bytes. Aborting!"
        )?;
        return Ok(Prepared::Done(ExitCode::IllegalParameters));
    }

    // 1. Parse the command line, subsection parameters included
    //    (TOPPBase.cpp:182-191, 2294-2305).
    let mut definitions: Vec<ParameterInformation> = spec
        .parameters()
        .iter()
        .filter(|p| !p.kind.is_layout())
        .cloned()
        .collect();
    for item in subsections.iter()? {
        definitions.push(ParameterInformation::from_param_entry(
            item.entry, &item.key,
        ));
    }
    let command_line = match parse_command_line(arguments, &definitions)? {
        Ok(parsed) => parsed,
        Err(failure) => {
            writeln!(
                err,
                "Invalid parameter values ({}): {}. Aborting!",
                failure.kind, failure.message
            )?;
            print_usage::<T>(err, spec, false)?;
            return Ok(Prepared::Done(ExitCode::IllegalParameters));
        }
    };
    for warning in &command_line.warnings {
        writeln!(err, "{warning}")?;
    }
    let cmd = &command_line.values;
    let given = |name: &str| cmd.exists(name).unwrap_or(false);

    // 2. A bare invocation (TOPPBase.cpp:227-232) is not yet refused here: it
    //    still reaches the required-parameter check. See docs/TOPP_CLI_SUPPORT.md.

    // 3. Usage requests short-circuit before any validation (235-239).
    if given("-help") || given("-helphelp") {
        print_usage::<T>(out, spec, given("-helphelp"))?;
        return Ok(Prepared::Done(ExitCode::ExecutionOk));
    }
    // 4. Unknown options and trailing text (241-255).
    if !command_line.unknown.is_empty() {
        writeln!(
            err,
            "Unknown option(s) '{}' given. Aborting!",
            list_text(&command_line.unknown)
        )?;
        print_usage::<T>(err, spec, false)?;
        return Ok(Prepared::Done(ExitCode::IllegalParameters));
    }
    if !command_line.misc.is_empty() {
        writeln!(
            err,
            "Trailing text argument(s) '{}' given. Aborting!",
            list_text(&command_line.misc)
        )?;
        print_usage::<T>(err, spec, false)?;
        return Ok(Prepared::Done(ExitCode::IllegalParameters));
    }

    // 5. Write commands run before any INI file is applied (265-268).
    let defaults = default_parameters::<T>(spec, &subsections)?;
    if let Some(code) = write_commands::<T>(cmd, &defaults, err)? {
        return Ok(Prepared::Done(code));
    }

    // 6. INI merge: command line, then the instance, common-tool and common
    //    sections, each adding only what is not yet present (274-333).
    let mut final_param = cmd.checked_clone()?;
    let mut ini = None;
    if given("ini") {
        let path = text_value(cmd, "ini")?;
        let loaded = match load_ini(&path, err)? {
            Ok(loaded) => loaded,
            Err(code) => return Ok(Prepared::Done(code)),
        };
        warn_if_not_applicable(&loaded, &location, err)?;
        final_param.merge(&loaded.copy(&location, true)?)?;
        final_param.merge(&loaded.copy(&format!("common:{}:", T::NAME), true)?)?;
        final_param.merge(&loaded.copy("common:", true)?)?;
        ini = Some(loaded);
    }
    if final_param.exists("ini")? {
        final_param.remove("ini")?;
    }

    // 7. Strict update of the defaults (338-343).
    let mut param = defaults.copy(&location, true)?;
    let diagnostics = update_diagnostics(&param, &final_param, UpdateMode::Strict)?;
    let report = param.update_with_options(
        &final_param,
        ParamUpdateOptions {
            verbose: false,
            add_unknown: false,
            fail_on_invalid_values: true,
            fail_on_unknown_parameters: true,
        },
    )?;
    for line in &diagnostics {
        writeln!(err, "{line}")?;
    }
    if !report.success {
        if diagnostics.is_empty() {
            for message in &report.messages {
                writeln!(err, "{message}")?;
            }
        }
        writeln!(
            err,
            "Parameters passed to '{}' are invalid. To prevent usage of wrong defaults, please update/fix the parameters!",
            T::NAME
        )?;
        return Ok(Prepared::Done(ExitCode::IllegalParameters));
    }
    if let Some(loaded) = &ini {
        warn_on_version_mismatch::<T>(loaded, out)?;
    }

    // 8. Parameter and file checks. The source runs them lazily as main_ reads
    //    each option (1388-1417, 1498-1614, 1968-2013); they run eagerly here.
    if let Some(code) = validate(spec, &param, err)? {
        return Ok(Prepared::Done(code));
    }
    Ok(Prepared::Run(Box::new(ToolContext::new(
        T::NAME,
        T::VERSION,
        param,
    ))))
}

/// The tool-description writers the source registers but this port refuses.
const DESCRIPTION_WRITERS: [&str; 5] = [
    "write_ctd",
    "write_nested_cwl",
    "write_cwl",
    "write_nested_json",
    "write_json",
];

/// Source `handleWriteCommands_` (`TOPPBase.cpp:2546-2686`).
///
/// `-write_ini` writes the defaults, updated leniently from `-ini` when given,
/// and never includes other command-line values. The CTD, CWL and JSON writers
/// are not ported: each request is refused with [`ExitCode::InternalError`],
/// which is what the oracle build without TDL support reports for four of them.
fn write_commands<T: Tool>(
    cmd: &Param,
    defaults: &Param,
    err: &mut dyn Write,
) -> Result<Option<ExitCode>> {
    if cmd.exists("write_ini")? {
        let path = text_value(cmd, "write_ini")?;
        if let Some(code) = output_file_writable(&path, "write_ini", err) {
            return Ok(Some(code));
        }
        let mut written = defaults.checked_clone()?;
        if cmd.exists("ini")? {
            let ini_path = text_value(cmd, "ini")?;
            let loaded = match load_ini(&ini_path, err)? {
                Ok(loaded) => loaded,
                Err(code) => return Ok(Some(code)),
            };
            warn_if_not_applicable(&loaded, &ini_location::<T>(), err)?;
            let diagnostics = update_diagnostics(&written, &loaded, UpdateMode::Lenient)?;
            written.update(&loaded, false)?;
            for line in &diagnostics {
                writeln!(err, "{line}")?;
            }
        }
        return Ok(Some(match paramxml::store(&path, &written) {
            Ok(()) => ExitCode::ExecutionOk,
            Err(error) => run_failure(&error, err),
        }));
    }
    for name in DESCRIPTION_WRITERS {
        if cmd.exists(name)? {
            writeln!(
                err,
                "Unable to initialize or run {}: '-{name}' is not supported, because the CTD, CWL and JSON tool description writers are not ported.",
                T::NAME
            )?;
            return Ok(Some(ExitCode::InternalError));
        }
    }
    Ok(None)
}

/// Load an INI file with the exit codes the source's run-phase catch assigns:
/// a missing file is [`ExitCode::InputFileNotFound`] and a malformed one
/// [`ExitCode::InputFileCorrupt`] (`TOPPBase.cpp:296`, `436-465`).
fn load_ini(path: &str, err: &mut dyn Write) -> Result<std::result::Result<Param, ExitCode>> {
    if !file::exists(path) {
        writeln!(
            err,
            "Error: File not found (the file '{path}' does not exist)"
        )?;
        return Ok(Err(ExitCode::InputFileNotFound));
    }
    Ok(match paramxml::load(path) {
        Ok(loaded) => Ok(loaded),
        Err(error) => Err(run_failure(&error, err)),
    })
}

/// Source `checkIfIniParametersAreApplicable_` (`TOPPBase.cpp:1957-1966`).
fn warn_if_not_applicable(ini: &Param, location: &str, err: &mut dyn Write) -> Result<()> {
    if ini.copy(location, false)?.is_empty() {
        writeln!(
            err,
            "Warning: The provided INI file does not contain any parameters specific for this tool (expected in '{location}'). Please check your .ini file. The default parameters for this tool will be applied."
        )?;
    }
    Ok(())
}

/// The INI-version notice of `TOPPBase.cpp:355-366`, written to the output
/// stream as the source's `writeLogInfo_` does.
fn warn_on_version_mismatch<T: Tool>(ini: &Param, out: &mut dyn Write) -> Result<()> {
    let key = format!("{}:version", T::NAME);
    if ini.exists(&key)? {
        let file_version = ini.value(&key)?.to_text(false)?;
        if file_version != T::VERSION {
            writeln!(
                out,
                "Warning: Parameters file version ({file_version}) does not match the version of this tool ({}).",
                T::VERSION
            )?;
            writeln!(
                out,
                "Your current parameters are still valid, but there might be new valid values or even new parameters. Upgrading the INI might be useful."
            )?;
        }
    }
    Ok(())
}

/// Which source `Param::update` call a diagnostic listing mirrors.
#[derive(Clone, Copy, PartialEq, Eq)]
enum UpdateMode {
    /// `update(p, verbose=false, add_unknown=false, fail_on_invalid_values=true,
    /// fail_on_unknown_parameters=true)`, as `TOPPBase::main`.
    Strict,
    /// `update(p, add_unknown=false)`: verbose and lenient, as `-write_ini`.
    Lenient,
}

/// Protected `:version` and TOPP `:type` keys, which an update never overrides.
fn protected_suffix(key: &str) -> Option<&'static str> {
    if key.ends_with(":version") {
        Some(":version")
    } else if key.ends_with(":type") && key.matches(':').count() >= 2 {
        Some(":type")
    } else {
        None
    }
}

/// Source `ParamEntry::isValid` message (`Param.cpp:56-166`), or `None` when
/// the entry satisfies its restrictions.
fn source_validity_message(entry: &ParamEntry) -> Option<String> {
    let path_tag = |with_prefix: bool| {
        entry.tags.contains("input file")
            || entry.tags.contains("output file")
            || (with_prefix && entry.tags.contains("output prefix"))
    };
    let int_outside = |x: i32| {
        (entry.min_int != -i32::MAX && x < entry.min_int)
            || (entry.max_int != i32::MAX && x > entry.max_int)
    };
    let float_outside = |x: f64| {
        (entry.min_float != -f64::MAX && x < entry.min_float)
            || (entry.max_float != f64::MAX && x > entry.max_float)
    };
    let invalid_string = |value: &str| {
        format!(
            "Invalid string parameter value '{value}' for parameter '{}' given! Valid values are: '{}'.",
            entry.name,
            entry.valid_strings.join(",")
        )
    };
    let invalid_int = |value: i32| {
        format!(
            "Invalid integer parameter value '{value}' for parameter '{}' given! The valid range is: [{}:{}].",
            entry.name, entry.min_int, entry.max_int
        )
    };
    let invalid_float = |value: f64| {
        format!(
            "Invalid double parameter value '{value:.6}' for parameter '{}' given! The valid range is: [{:.6}:{:.6}].",
            entry.name, entry.min_float, entry.max_float
        )
    };
    let narrow =
        |value: i64| i32::try_from(value).unwrap_or(if value < 0 { i32::MIN } else { i32::MAX });
    match &entry.value {
        ParamValue::String(value)
            if !entry.valid_strings.is_empty()
                && !entry.valid_strings.contains(value)
                && !path_tag(true) =>
        {
            Some(invalid_string(value))
        }
        ParamValue::StringList(values) if !entry.valid_strings.is_empty() && !path_tag(false) => {
            values
                .iter()
                .find(|value| !entry.valid_strings.contains(value))
                .map(|value| invalid_string(value))
        }
        ParamValue::Integer(value) => {
            let value = narrow(*value);
            int_outside(value).then(|| invalid_int(value))
        }
        ParamValue::IntegerList(values) => values
            .iter()
            .find(|value| int_outside(**value))
            .map(|value| invalid_int(*value)),
        ParamValue::Float(value) => float_outside(*value).then(|| invalid_float(*value)),
        ParamValue::FloatList(values) => values
            .iter()
            .find(|value| float_outside(**value))
            .map(|value| invalid_float(*value)),
        _ => None,
    }
}

/// The diagnostics the source `Param::update` streams while applying
/// `outdated` to `current` (`Param.cpp:1216-1374`).
///
/// [`Param::update_with_options`] decides success and applies the values; its
/// messages are worded natively, so this lists the source wording from the same
/// data: the entry found by exact key or by a unique leaf name, its value type
/// and its restrictions. Tools' users and the source's registered failure tests
/// match these words.
fn update_diagnostics(current: &Param, outdated: &Param, mode: UpdateMode) -> Result<Vec<String>> {
    let strict = mode == UpdateMode::Strict;
    let stream_text = |value: &ParamValue| value.to_stream_text().unwrap_or_default();
    let mut lines = Vec::new();
    for item in outdated.iter()? {
        let target = if current.exists(&item.key)? {
            if let Some(suffix) = protected_suffix(&item.key) {
                if current.value(&item.key)? != &item.entry.value {
                    lines.push(format!(
                        "Warning: for '{suffix}' entry, augmented and Default Ini-File differ in value. Default value will not be altered!"
                    ));
                }
                continue;
            }
            item.key.clone()
        } else {
            let unique = match current.find_first(&item.entry.name)? {
                Some(first) if current.find_next(&item.entry.name, first.entry)?.is_none() => {
                    Some(first.key)
                }
                _ => None,
            };
            match unique {
                Some(key) => {
                    lines.push(format!("Found '{}' as '{key}' in new param.", item.key));
                    key
                }
                None => {
                    lines.push(if strict {
                        format!(
                            "Unknown (or deprecated) Parameter '{}' given in outdated parameter file!",
                            item.key
                        )
                    } else {
                        format!(
                            "Unknown (or deprecated) Parameter '{}' given in outdated parameter file! Ignoring parameter. ",
                            item.key
                        )
                    });
                    continue;
                }
            }
        };
        let entry = current.entry(&target)?;
        if entry.value.value_type() != item.entry.value.value_type() {
            lines.push(format!("Parameter '{}' has changed value type!", item.key));
            lines.push(if strict {
                " Updating failed!".to_owned()
            } else {
                " Ignoring invalid value (using new default)!".to_owned()
            });
            continue;
        }
        if entry.value == item.entry.value {
            continue;
        }
        let mut candidate = entry.clone();
        candidate.value = item.entry.value.clone();
        match source_validity_message(&candidate) {
            Some(message) if strict => lines.push(format!("{message} Updating failed!")),
            Some(message) => lines.push(format!(
                "{message} Ignoring invalid value (using new default '{}')!",
                stream_text(&entry.value)
            )),
            None if !strict => lines.push(format!(
                "Default-Parameter '{target}' overridden: '{}' --> '{}'!",
                stream_text(&entry.value),
                stream_text(&item.entry.value)
            )),
            None => {}
        }
    }
    Ok(lines)
}

/// Required-value, restriction and file checks, in registration order. Returns
/// the terminal exit code when the run cannot proceed.
fn validate(spec: &ToolSpec, param: &Param, err: &mut dyn Write) -> Result<Option<ExitCode>> {
    for entry in spec.parameters() {
        if entry.kind.is_layout() || !param.exists(&entry.name)? {
            continue;
        }
        let value = param.value(&entry.name)?;
        // An empty list counts as missing, as getStringList_, getIntList_ and
        // getDoubleList_ throw RequiredParameterNotGiven for it (TOPPBase.cpp:1630, 1661, 1695).
        let empty = match value {
            ParamValue::String(text) => text.is_empty(),
            ParamValue::StringList(values) => values.is_empty(),
            ParamValue::IntegerList(values) => values.is_empty(),
            ParamValue::FloatList(values) => values.is_empty(),
            ParamValue::Empty => true,
            _ => false,
        };
        if entry.required && empty {
            // Only getStringOption_ names the valid values (TOPPBase.cpp:1398-1406);
            // for a file option the source keeps its formats there.
            let string_option = matches!(
                entry.kind,
                ParameterType::String
                    | ParameterType::InputFile
                    | ParameterType::OutputFile
                    | ParameterType::OutputPrefix
            );
            let restrictions = if entry.valid_strings.is_empty() {
                &entry.valid_formats
            } else {
                &entry.valid_strings
            };
            let name = if !string_option || restrictions.is_empty() {
                format!("'{}'", entry.name)
            } else {
                format!("'{}' [valid: {}]", entry.name, restrictions.join(", "))
            };
            writeln!(
                err,
                "Error: The required parameter {name} was not given or is empty!"
            )?;
            return Ok(Some(ExitCode::MissingParameters));
        }
        if empty {
            continue;
        }
        let code = match (entry.kind, value) {
            (ParameterType::String, ParamValue::String(text))
                if !entry.valid_strings.is_empty() && !entry.valid_strings.contains(text) =>
            {
                writeln!(
                    err,
                    "Invalid parameter: Invalid value '{text}' for string parameter '{}' given. Valid strings are: '{}'.",
                    entry.name,
                    entry.valid_strings.join("', '")
                )?;
                Some(ExitCode::IllegalParameters)
            }
            (ParameterType::Int, ParamValue::Integer(number)) => int_range(entry, *number, err)?,
            (ParameterType::IntList, ParamValue::IntegerList(numbers)) => {
                let mut code = None;
                for number in numbers {
                    code = int_range(entry, i64::from(*number), err)?;
                    if code.is_some() {
                        break;
                    }
                }
                code
            }
            (ParameterType::Double, ParamValue::Float(number)) => float_range(entry, *number, err)?,
            (ParameterType::DoubleList, ParamValue::FloatList(numbers)) => {
                let mut code = None;
                for number in numbers {
                    code = float_range(entry, *number, err)?;
                    if code.is_some() {
                        break;
                    }
                }
                code
            }
            (ParameterType::InputFile, ParamValue::String(path)) => input_path(entry, path, err)?,
            (ParameterType::InputFileList, ParamValue::StringList(paths)) => {
                let mut code = None;
                for path in paths {
                    code = input_path(entry, path, err)?;
                    if code.is_some() {
                        break;
                    }
                }
                code
            }
            (ParameterType::OutputFile, ParamValue::String(path)) => {
                match output_file_writable(path, &entry.name, err) {
                    Some(code) => Some(code),
                    None => output_extension(entry, path, err)?,
                }
            }
            (ParameterType::OutputPrefix, ParamValue::String(path)) => {
                output_file_writable(&format!("{path}_0"), &entry.name, err)
            }
            (ParameterType::OutputFileList, ParamValue::StringList(paths)) => {
                let mut code = None;
                for path in paths {
                    code = output_extension(entry, path, err)?;
                    if code.is_some() {
                        break;
                    }
                }
                code
            }
            _ => None,
        };
        if code.is_some() {
            return Ok(code);
        }
    }
    Ok(None)
}

/// Source `getIntOption_` range check (`TOPPBase.cpp:1486-1493`).
fn int_range(
    entry: &ParameterInformation,
    number: i64,
    err: &mut dyn Write,
) -> Result<Option<ExitCode>> {
    let min = i64::from(entry.min_int.unwrap_or(-i32::MAX));
    let max = i64::from(entry.max_int.unwrap_or(i32::MAX));
    if number < min || number > max {
        writeln!(
            err,
            "Invalid parameter: Invalid value '{number}' for integer parameter '{}' given. Out of valid range: '{min}'-'{max}'.",
            entry.name
        )?;
        return Ok(Some(ExitCode::IllegalParameters));
    }
    Ok(None)
}

/// Source `getDoubleOption_` range check (`TOPPBase.cpp:1459-1466`).
fn float_range(
    entry: &ParameterInformation,
    number: f64,
    err: &mut dyn Write,
) -> Result<Option<ExitCode>> {
    let min = entry.min_float.unwrap_or(-f64::MAX);
    let max = entry.max_float.unwrap_or(f64::MAX);
    if number < min || number > max {
        writeln!(
            err,
            "Invalid parameter: Invalid value '{number}' for float parameter '{}' given. Out of valid range: '{min}'-'{max}'.",
            entry.name
        )?;
        return Ok(Some(ExitCode::IllegalParameters));
    }
    Ok(None)
}

/// Input readability, then the registered extension (`TOPPBase.cpp:1529-1592`).
fn input_path(
    entry: &ParameterInformation,
    path: &str,
    err: &mut dyn Write,
) -> Result<Option<ExitCode>> {
    if !entry.tags.iter().any(|tag| tag == "skipexists") {
        if let Some(code) = input_file_readable(path, &entry.name, err) {
            return Ok(Some(code));
        }
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
    Ok(None)
}

/// The registered output extension (`TOPPBase.cpp:1595-1608`).
fn output_extension(
    entry: &ParameterInformation,
    path: &str,
    err: &mut dyn Write,
) -> Result<Option<ExitCode>> {
    if !context::extension_allowed(path, &entry.valid_formats) {
        writeln!(
            err,
            "Invalid parameter: Invalid output file extension for file '{path}'. Valid file extensions are: '{}'.",
            entry.valid_formats.join("','")
        )?;
        return Ok(Some(ExitCode::IllegalParameters));
    }
    Ok(None)
}

/// A failure before the tool body that no phase handles itself: the source's
/// initialisation catch (`TOPPBase.cpp:505-508`).
fn initialisation_failure<T: Tool>(error: &Error, err: &mut dyn Write) -> ExitCode {
    let _ = writeln!(err, "Unable to initialize or run {}: {error}", T::NAME);
    ExitCode::IllegalParameters
}

/// A failure while the tool runs, mapped as the source's run-phase catch
/// (`TOPPBase.cpp:430-499`).
///
/// A parse failure is `INPUT_FILE_CORRUPT`, a missing file
/// `INPUT_FILE_NOT_FOUND`, a permission failure `CANNOT_WRITE_OUTPUT_FILE`, an
/// invalid value or range `ILLEGAL_PARAMETERS`, missing information
/// `MISSING_PARAMETERS`, and any other I/O failure `UNKNOWN_ERROR`. `Unsupported`
/// and `UnsortedData` are `INCOMPATIBLE_INPUT_DATA`, the code the source tools
/// return explicitly for those conditions.
fn run_failure(error: &Error, err: &mut dyn Write) -> ExitCode {
    let (code, text) = match error {
        Error::Parse { .. } => (
            ExitCode::InputFileCorrupt,
            format!("Error: Unable to read file ({error})"),
        ),
        Error::Io(e) if e.kind() == std::io::ErrorKind::NotFound => (
            ExitCode::InputFileNotFound,
            format!("Error: File not found ({error})"),
        ),
        Error::Io(e) if e.kind() == std::io::ErrorKind::PermissionDenied => (
            ExitCode::CannotWriteOutputFile,
            format!("Error: Unable to write file ({error})"),
        ),
        Error::Io(_) => (
            ExitCode::UnknownError,
            format!("Error: Unexpected internal error ({error})"),
        ),
        Error::Unsupported(_) | Error::UnsortedData => {
            (ExitCode::IncompatibleInputData, format!("Error: {error}"))
        }
        Error::InvalidValue(_) | Error::InvalidRange(_) => (
            ExitCode::IllegalParameters,
            format!("Invalid parameter: {error}"),
        ),
        Error::MissingInformation(_) => (ExitCode::MissingParameters, format!("Error: {error}")),
    };
    let _ = writeln!(err, "{text}");
    code
}

/// Run a tool against explicit arguments and streams. `arguments[0]` is the
/// executable name, as in `main(argc, argv)`.
///
/// The phases and their exit codes follow `TOPPBase::main`; see the module
/// documentation. Usage for a successful `--help` goes to `out`; every
/// diagnostic, including usage after a command-line error, goes to `err`.
pub fn run_with<T: Tool>(
    arguments: &[String],
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> ExitCode {
    let spec = match tool_spec::<T>() {
        Ok(spec) => spec,
        Err(error) => return initialisation_failure::<T>(&error, err),
    };
    let ctx = match prepare::<T>(&spec, arguments, out, err) {
        Ok(Prepared::Run(ctx)) => ctx,
        Ok(Prepared::Done(code)) => return code,
        Err(error) => return initialisation_failure::<T>(&error, err),
    };
    match T::run_io(&ctx, out, err) {
        Ok(code) => code,
        Err(error) => run_failure(&error, err),
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
