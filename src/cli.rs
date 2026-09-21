// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The TOPP command-line framework: registration, the run lifecycle, usage
//! text, tool descriptions, the tool registry and processing annotations.
//!
//! This is the native form of the `OpenMS4-cli` package's
//! `APPLICATIONS/TOPPBase.h`, `TOPPBase_defs.h`, `ParameterInformation.h` and
//! `ToolHandler.h`. [`run_with`](crate::cli::run_with) follows
//! `TOPPBase::main` (`TOPPBase.cpp:141-514`) phase by phase:
//!
//! 1. the product version from the tool registry, then registration of the
//!    tool's own parameters and the common ones, under `Common TOPP options:`
//!    for a registered tool and `Common UTIL options:` otherwise;
//! 2. the command-line parse closure, which also accepts every subsection
//!    parameter as `-section:name` and rejects unknown options and trailing
//!    text;
//! 3. `-write_ini` and the tool-description writers: `-write_ctd` writes the
//!    Common Tool Description, and the four CWL and JSON writers are refused
//!    as the Release build refuses them without TDL;
//! 4. the INI merge of the instance, `common:<tool>:` and `common:` sections
//!    and a strict update of the defaults, which rejects unknown parameters and
//!    invalid values, then the source's `checkParam_` warnings;
//! 5. parameter and file validation, then the tool body, then the run-time
//!    and peak-memory line.
//!
//! A failure maps to the exit code of the phase it occurs in, as the source's
//! two catch blocks do. `-log` and `-debug` write the source's log file. See
//! `docs/TOPP_CLI_SUPPORT.md` for the supported source subset and the
//! exit-code table.
//!
//! One phase has no source counterpart and precedes all of these: a check that
//! this processor has the instructions this binary was built to use. See
//! [`crate::system::cpu_features`] and `docs/FMA_BUILD_FLAG.md`.

mod context;
mod defs;
mod logging;
mod param_ctd;
mod parameter;
mod processing;
mod spec;
mod tool_description_file;
mod tool_handler;
/// The ported TOPP tools, one library type per executable.
pub mod tools;
mod usage;

pub use context::{
    TEST_MODE_UNIQUE_ID_SEED, ToolContext, input_file_readable, output_file_writable, parse_range,
    parse_range_int,
};
pub use defs::{CITE_OPENMS, Citation};
pub use logging::{LOG_SEPARATOR, ToolLog};
pub use param_ctd::{MAX_CTD_BYTES, ParamCtdFile};
pub use parameter::{ExitCode, ParameterInformation, ParameterType};
pub use processing::{
    AddDataProcessing, TEST_MODE_COMPLETION_TIME, TEST_MODE_PARAMETER_KEY,
    TEST_MODE_PARAMETER_VALUE, TEST_MODE_VERSION,
};
pub use spec::ToolSpec;
pub use tool_description_file::{
    LoadedToolDescriptions, MAX_TTD_BYTES, MAX_TTD_DEPTH, MAX_TTD_ELEMENTS, ToolDescriptionFile,
};
pub use tool_handler::{
    BUILTIN_MANIFEST, BUILTIN_MANIFEST_NAME, MAX_MANIFEST_BYTES, MAX_MANIFEST_ROWS, PackageTool,
    ToolHandler, ToolListType, ToolRegistrySources,
};

use crate::data_structures::ToolInfo;
use crate::format::file_handler::FileHandler;
use crate::format::file_types::{FileType, type_by_file_name};
use crate::format::paramxml;
use crate::param::{Param, ParamEntry, ParamUpdateOptions, ParamValue};
use crate::system::file;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::io::Write;
use std::sync::Arc;

/// The product version the pinned `topp` package gives every tool it installs
/// (`tools.json` `"version": "1.0.0"`), and so the version every row of the
/// built-in manifest ([`BUILTIN_MANIFEST`]) records.
///
/// The lifecycle does not read this constant: as the source
/// (`TOPPBase.cpp:143-150`), it asks the tool registry for the tool's product
/// version and falls back to the core version,
/// [`CORE_SDK_VERSION`](crate::CORE_SDK_VERSION), for a tool no manifest
/// registers.
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
/// `main_`, with the name, description and citations as constructor
/// arguments. A native tool implements this trait instead, so registration and
/// execution are plain functions with no inherited mutable state.
///
/// The source constructor's `official` and `toolhandler_test` arguments are
/// not ported: `TOPPBase` stores both and reads neither (`toolhandler_test` is
/// documented as a compatibility argument). The product version is not a
/// property of the tool either: it comes from the tool registry
/// ([`ToolHandler::get_tool_version`]).
pub trait Tool {
    /// Executable name, which is also the INI section name.
    const NAME: &'static str;
    /// One-line description shown in usage output and tool descriptions.
    const DESCRIPTION: &'static str;
    /// Publications specific to this tool, printed by `--help` after the
    /// OpenMS citation and written to its tool description (source
    /// constructor argument `citations`). Empty by default.
    const CITATIONS: &'static [Citation] = &[];

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

/// The heading of the common options: `Common TOPP options:` for a tool the
/// registry lists, `Common UTIL options:` otherwise (`TOPPBase.cpp:160-163`).
fn common_heading(registered: bool) -> &'static str {
    if registered {
        "Common TOPP options:"
    } else {
        "Common UTIL options:"
    }
}

/// The parameters every TOPP tool registers, in source order
/// (`TOPPBase.cpp:158-179`).
fn register_common(spec: &mut ToolSpec, registered: bool) -> Result<()> {
    spec.add_empty_line()?;
    spec.add_text(common_heading(registered))?;
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

/// Build the full registration of a registered TOPP tool: its own parameters,
/// then the common ones under `Common TOPP options:`.
///
/// # Errors
///
/// Propagates the tool's registration error, or a duplicate of a common name.
pub fn tool_spec<T: Tool>() -> Result<ToolSpec> {
    tool_spec_for::<T>(true)
}

/// Build the full registration of a tool, with the common options under the
/// heading for a tool the registry lists (`registered`) or not.
///
/// # Errors
///
/// As [`tool_spec`].
pub fn tool_spec_for<T: Tool>(registered: bool) -> Result<ToolSpec> {
    let mut spec = ToolSpec::new();
    T::register(&mut spec)?;
    register_common(&mut spec, registered)?;
    Ok(spec)
}

/// The INI section of a tool instance, as `getToolPrefix`.
fn location_of(name: &str, instance: i64) -> String {
    format!("{name}:{instance}:")
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

/// The product version and the verbose version line of a tool, as
/// `TOPPBase::main` and the constructor build them (`TOPPBase.cpp:119-126`,
/// `143-150`).
///
/// A registered tool reports its product version, and its verbose version
/// names the core version and the short core revision. An unregistered tool
/// reports the core version; its verbose version is the core version and
/// `, Revision: <short revision>`. The source puts the library's build time
/// (`VersionInfo::getTime`) between the two, which this crate does not
/// record.
///
/// # Errors
///
/// As [`ToolHandler::get_tool_version`].
pub fn product_versions(registry: &ToolHandler, name: &str) -> Result<(String, String)> {
    let revision = crate::CORE_SDK_REVISION
        .get(..7)
        .unwrap_or(crate::CORE_SDK_REVISION);
    let product = registry.get_tool_version(name)?;
    Ok(if product.is_empty() {
        (
            crate::CORE_SDK_VERSION.to_owned(),
            format!("{}, Revision: {revision}", crate::CORE_SDK_VERSION),
        )
    } else {
        let verbose = format!(
            "{product} (OpenMS core {}, revision {revision})",
            crate::CORE_SDK_VERSION
        );
        (product, verbose)
    })
}

/// What one run knows beyond the tool's registration: the registry, the
/// product version, the INI location and the log.
struct Session<'a> {
    registry: &'a ToolHandler,
    name: &'static str,
    version: String,
    verbose_version: String,
    instance: i64,
    location: String,
    log: Arc<ToolLog>,
}

impl Session<'_> {
    /// Write the log's `Writing to` notice to standard output, where the
    /// source's `enableLogging_` prints it.
    fn notice(&self, notice: Option<String>, out: &mut dyn Write) -> Result<()> {
        if let Some(notice) = notice {
            writeln!(out, "{notice}")?;
        }
        Ok(())
    }
    /// Source `writeLogInfo_`: standard output and the log file.
    fn info(&self, out: &mut dyn Write, text: &str) -> Result<()> {
        writeln!(out, "{text}")?;
        let notice = self.log.line(text);
        self.notice(notice, out)
    }
    /// Source `writeLogWarn_` and `writeLogError_`: the error stream and the
    /// log file.
    fn warn(&self, out: &mut dyn Write, err: &mut dyn Write, text: &str) -> Result<()> {
        writeln!(err, "{text}")?;
        let notice = self.log.line(text);
        self.notice(notice, out)
    }
    /// Source `writeDebug_(text, level)`: the log file only.
    fn debug(&self, out: &mut dyn Write, text: &str, level: u32) -> Result<()> {
        let notice = self.log.debug(text, level);
        self.notice(notice, out)
    }
    /// Source `writeDebug_(text, param, level)`.
    fn debug_param(
        &self,
        out: &mut dyn Write,
        text: &str,
        param: &Param,
        level: u32,
    ) -> Result<()> {
        let notice = self.log.debug_param(text, param, level);
        self.notice(notice, out)
    }
}

/// The per-user defaults of a tool (source `getToolUserDefaults_`,
/// `TOPPBase.cpp:2276-2286`): `<user directory>/<tool>.ini` when that file is
/// readable, else nothing.
///
/// The user directory is `OPENMS_HOME_PATH`, the `home_dir` of `OpenMS.ini`,
/// or the home directory ([`file::FileContext::get_user_directory`]). The file
/// name is built as the source builds it, from the user directory with a
/// trailing `/` and another `/`, so a diagnostic names `<home>//<tool>.ini`.
///
/// Failures take the exit codes the Release build takes for them
/// (`../oracle/toppbase-completion`, cases `ud_*`): a file that does not parse
/// is `INPUT_FILE_CORRUPT` (3), the source's `ParseError`; a directory by that
/// name is `INTERNAL_ERROR` (12), because libstdc++ throws a `std::exception`
/// that only the initialisation catch handles. An unreadable file is skipped
/// silently, as the source's `File::readable` check skips it.
fn user_defaults(
    session: &Session<'_>,
    err: &mut dyn Write,
) -> Result<std::result::Result<Option<Param>, ExitCode>> {
    let Ok(directory) = session.registry.sources().file_context.get_user_directory() else {
        return Ok(Ok(None));
    };
    let path = format!("{}//{}.ini", directory.display(), session.name);
    let real = directory.join(format!("{}.ini", session.name));
    if !file::exists(&real) || !file::readable(&real) {
        return Ok(Ok(None));
    }
    if file::is_directory(&real) {
        writeln!(
            err,
            "Unable to initialize or run {}: basic_filebuf::underflow error reading the file: Is a directory",
            session.name
        )?;
        return Ok(Err(ExitCode::InternalError));
    }
    match paramxml::load(&real) {
        Ok(param) => Ok(Ok(Some(param))),
        Err(Error::Parse { message, .. }) => {
            writeln!(
                err,
                "Error: Unable to read file (While loading '{path}': {message} in: {path})"
            )?;
            Ok(Err(ExitCode::InputFileCorrupt))
        }
        Err(error) => Ok(Err(run_failure(&error, err))),
    }
}

/// Source `getDefaultParameters_` (`TOPPBase.cpp:2097-2256`): the registered
/// parameters under `<tool>:<instance>:`, the tool version item, the tool and
/// instance section descriptions, a `type` item when the command line carries
/// one, the subsection defaults, and finally the per-user defaults, applied
/// leniently with the source's verbose update diagnostics on the error
/// stream.
fn default_parameters<T: Tool>(
    spec: &ToolSpec,
    subsections: &Param,
    session: &Session<'_>,
    type_value: Option<&ParamValue>,
    err: &mut dyn Write,
) -> Result<std::result::Result<Param, ExitCode>> {
    let location = &session.location;
    let mut defaults = spec.to_param_at(T::NAME, session.instance)?;
    defaults.set_value(
        &format!("{}:version", T::NAME),
        ParamValue::String(session.version.clone()),
        "Version of the tool that generated this parameters file.",
        &["advanced".to_owned()],
    )?;
    defaults.set_section_description(T::NAME, T::DESCRIPTION)?;
    defaults.set_section_description(
        location.trim_end_matches(':'),
        &format!("Instance '{}' section for '{}'", session.instance, T::NAME),
    )?;
    if let Some(value) = type_value {
        defaults.set_value(&format!("{location}type"), value.checked_clone()?, "", &[])?;
    }
    if !subsections.is_empty() {
        defaults.insert(location, subsections)?;
        for (name, description) in spec.subsections() {
            if subsections.has_section(name)? {
                defaults.set_section_description(&format!("{location}{name}"), description)?;
            }
        }
    }
    match user_defaults(session, err)? {
        Ok(Some(user)) => {
            let diagnostics = update_diagnostics(&defaults, &user, UpdateMode::Lenient)?;
            defaults.update(&user, false)?;
            for line in &diagnostics {
                writeln!(err, "{line}")?;
            }
        }
        Ok(None) => {}
        Err(code) => return Ok(Err(code)),
    }
    Ok(Ok(defaults))
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
    // Text left over after each option, gathered in reverse (see below).
    let mut misc_reversed: Vec<String> = Vec::new();
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
        // The source inserts the rest of the queue at the front of `misc`
        // (2436-2439), which costs time quadratic in the token count. Pushing
        // each chunk reversed and reversing once at the end gives the same order
        // in linear time.
        misc_reversed.extend(queue.drain(..).rev().map(str::to_owned));
    }
    // What is left are leading text arguments (2442-2444).
    misc_reversed.extend(queue.drain(..).rev().map(str::to_owned));
    misc_reversed.reverse();
    parsed.misc = misc_reversed;
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

/// The verbose version line of a tool as `--help` prints it, with the version
/// taken from the process's tool registry ([`product_versions`]).
///
/// For a registered tool this is the product version, then the core version
/// and the short core revision (`TOPPBase.cpp:144-150`). The revision is this
/// crate's pinned core revision
/// ([`CORE_SDK_REVISION`](crate::CORE_SDK_REVISION)), where a C++ build
/// prints the revision it was built from. When the registry cannot be read,
/// the line of an unregistered tool is returned.
pub fn verbose_version<T: Tool>() -> String {
    let revision = crate::CORE_SDK_REVISION
        .get(..7)
        .unwrap_or(crate::CORE_SDK_REVISION);
    ToolHandler::from_environment()
        .and_then(|registry| product_versions(&registry, T::NAME))
        .map(|(_, verbose)| verbose)
        .unwrap_or_else(|_| format!("{}, Revision: {revision}", crate::CORE_SDK_VERSION))
}

/// Source `printUsage_` for a tool body that prints its own usage, as
/// `FileInfo` does: the usage text on `stream`, with the version line from the
/// process's tool registry ([`verbose_version`]).
fn print_usage<T: Tool>(
    stream: &mut dyn Write,
    spec: &ToolSpec,
    subsections: &Param,
    verbose: bool,
) -> Result<()> {
    usage::print(
        stream,
        T::NAME,
        T::DESCRIPTION,
        &verbose_version::<T>(),
        T::CITATIONS,
        spec,
        subsections,
        verbose,
    )?;
    Ok(())
}

/// Source `printUsage_` inside the lifecycle, which writes to standard error in
/// every case; the caller passes the error stream.
///
/// As the source, it first reads the `-helphelp` flag from the command line,
/// which at debug level 1 logs `Parameter '-helphelp' not found.` when the
/// flag is absent and then the flag's value.
fn print_usage_logged<T: Tool>(
    session: &Session<'_>,
    cmd: &Param,
    out: &mut dyn Write,
    err: &mut dyn Write,
    spec: &ToolSpec,
    subsections: &Param,
    verbose: bool,
) -> Result<()> {
    let given = cmd.exists("-helphelp").unwrap_or(false);
    if !given {
        session.debug(out, "Parameter '-helphelp' not found.", 1)?;
    }
    session.debug(
        out,
        &format!("Value of string option '-helphelp': {}", u8::from(given)),
        1,
    )?;
    usage::print(
        err,
        T::NAME,
        T::DESCRIPTION,
        &session.verbose_version,
        T::CITATIONS,
        spec,
        subsections,
        verbose,
    )?;
    Ok(())
}

/// Source `StringUtils::quote` with the default `ESCAPE` method: surround with
/// double quotes, escaping backslashes and double quotes.
fn quote(text: &str) -> String {
    format!("\"{}\"", text.replace('\\', "\\\\").replace('"', "\\\""))
}

/// Outcome of preparing a run: either parameters to execute with, or a
/// terminal exit code because usage, an INI file or a failure ended the run.
enum Prepared {
    Run(Box<ToolContext>),
    Done(ExitCode),
}

/// The lifecycle up to the tool body (`TOPPBase.cpp:157-408`).
fn prepare<T: Tool>(
    session: &mut Session<'_>,
    spec: &ToolSpec,
    arguments: &[String],
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<Prepared> {
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
    //    (TOPPBase.cpp:182-191, 2294-2305). No parameters are in force yet, so
    //    nothing reaches a log file until the parse has succeeded.
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
            print_usage_logged::<T>(session, &Param::new(), out, err, spec, &subsections, false)?;
            return Ok(Prepared::Done(ExitCode::IllegalParameters));
        }
    };
    for warning in &command_line.warnings {
        writeln!(err, "{warning}")?;
    }
    let mut cmd = command_line.values;

    // The command line is the parameter set in force (param_ = param_cmdline_):
    // its -instance, -debug and -log apply from here (TOPPBase.cpp:193-223).
    let cmd_int = |cmd: &Param, key: &str, fallback: i64| match cmd.value(key) {
        Ok(ParamValue::Integer(value)) => *value,
        _ => fallback,
    };
    session.instance = cmd_int(&cmd, "instance", 1);
    session.location = location_of(T::NAME, session.instance);
    session.log.set_location(&session.location);
    session
        .log
        .set_destination(text_value(&cmd, "log").ok().as_deref());
    session.debug(out, &format!("Instance: {}", session.instance), 1)?;
    session.debug(out, &format!("Ini_location: {}", session.location), 1)?;
    session.log.set_debug_level(cmd_int(&cmd, "debug", 0));
    session.debug(
        out,
        &format!("Debug level: {}", session.log.debug_level()),
        1,
    )?;
    let quoted: Vec<String> = arguments
        .iter()
        .map(|argument| {
            if argument.contains(' ') {
                quote(argument)
            } else {
                argument.clone()
            }
        })
        .collect();
    session.debug(out, &format!(" >> {}", quoted.join(" ")), 1)?;
    let given = |cmd: &Param, name: &str| cmd.exists(name).unwrap_or(false);

    // 2. A bare invocation prints usage and is refused (227-232). An empty
    //    argument list, argc 0, still runs, as in the source class test.
    if arguments.len() == 1 {
        print_usage_logged::<T>(session, &cmd, out, err, spec, &subsections, false)?;
        session.warn(out, err, "No options given. Aborting!")?;
        return Ok(Prepared::Done(ExitCode::IllegalParameters));
    }

    // 3. Usage requests short-circuit before any validation (235-239). The
    //    source prints usage to standard error here too.
    if given(&cmd, "-help") || given(&cmd, "-helphelp") {
        print_usage_logged::<T>(
            session,
            &cmd,
            out,
            err,
            spec,
            &subsections,
            given(&cmd, "-helphelp"),
        )?;
        return Ok(Prepared::Done(ExitCode::ExecutionOk));
    }
    // 4. Unknown options and trailing text (241-255).
    if !command_line.unknown.is_empty() {
        session.warn(
            out,
            err,
            &format!(
                "Unknown option(s) '{}' given. Aborting!",
                list_text(&command_line.unknown)
            ),
        )?;
        print_usage_logged::<T>(session, &cmd, out, err, spec, &subsections, false)?;
        return Ok(Prepared::Done(ExitCode::IllegalParameters));
    }
    if !command_line.misc.is_empty() {
        session.warn(
            out,
            err,
            &format!(
                "Trailing text argument(s) '{}' given. Aborting!",
                list_text(&command_line.misc)
            ),
        )?;
        print_usage_logged::<T>(session, &cmd, out, err, spec, &subsections, false)?;
        return Ok(Prepared::Done(ExitCode::IllegalParameters));
    }

    // 5. Write commands run before any INI file is applied (265-268).
    let mut defaults =
        match default_parameters::<T>(spec, &subsections, session, cmd.value("type").ok(), err)? {
            Ok(defaults) => defaults,
            Err(code) => return Ok(Prepared::Done(code)),
        };
    if let Some(code) = write_commands::<T>(session, &cmd, &defaults, out, err)? {
        return Ok(Prepared::Done(code));
    }

    // 6. INI merge: command line, then the instance, common-tool and common
    //    sections, each adding only what is not yet present (274-333).
    let location = session.location.clone();
    let mut ini = None;
    let mut ini_path = String::new();
    let (mut instance, mut common_tool, mut common) = (Param::new(), Param::new(), Param::new());
    if given(&cmd, "ini") {
        ini_path = text_value(&cmd, "ini")?;
        session.debug(out, &format!("INI file: {ini_path}"), 1)?;
        session.debug(out, &format!("INI location: {location}"), 1)?;
        let loaded = match load_ini(&ini_path, err)? {
            Ok(loaded) => loaded,
            Err(code) => return Ok(Prepared::Done(code)),
        };
        warn_if_not_applicable(session, &loaded, out, err)?;
        instance = loaded.copy(&location, true)?;
        session.debug_param(out, "Parameters from instance section:", &instance, 2)?;
        common_tool = loaded.copy(&format!("common:{}:", T::NAME), true)?;
        session.debug_param(
            out,
            "Parameters from common section with tool name:",
            &common_tool,
            2,
        )?;
        common = loaded.copy("common:", true)?;
        session.debug_param(
            out,
            "Parameters from common section without tool name:",
            &common,
            2,
        )?;
        // A `type` in the instance section reaches the command line when the
        // command line has none (310-311), and so the defaults (2237-2239).
        let type_key = format!("{location}type");
        if loaded.exists(&type_key)? && !given(&cmd, "type") {
            let value = loaded.value(&type_key)?.checked_clone()?;
            cmd.set_value("type", value.checked_clone()?, "", &[])?;
            if !defaults.exists(&type_key)? {
                defaults.set_value(&type_key, value, "", &[])?;
            }
        }
        ini = Some(loaded);
    }
    session.debug_param(out, "Initialize final param with cmd line:", &cmd, 2)?;
    let mut final_param = cmd.checked_clone()?;
    session.debug_param(out, "Merging instance section into param:", &instance, 2)?;
    final_param.merge(&instance)?;
    session.debug_param(
        out,
        "Merging common section with tool name into param:",
        &common_tool,
        2,
    )?;
    final_param.merge(&common_tool)?;
    session.debug_param(
        out,
        "Merging common section without tool name into param:",
        &common,
        2,
    )?;
    final_param.merge(&common)?;
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
    if final_param.exists("type")? && !param.exists("type")? {
        param.set_value("type", final_param.value("type")?.checked_clone()?, "", &[])?;
    }
    // The resolved parameters are in force from here: their -log applies.
    session
        .log
        .set_destination(text_value(&param, "log").ok().as_deref());

    // 8. The registration checks of the INI sections (350-353) and the version
    //    notice (355-366).
    check_param(session, spec, &instance, &ini_path, &location, out, err)?;
    check_param(
        session,
        spec,
        &common_tool,
        &ini_path,
        &format!("common:{}::", T::NAME),
        out,
        err,
    )?;
    check_param(session, spec, &common, &ini_path, "common:", out, err)?;
    if let Some(loaded) = &ini {
        warn_on_version_mismatch(session, loaded, out)?;
    }

    // 9. The framework's own reads of the resolved parameters (369-408).
    let flag = |param: &Param, name: &str| matches!(param.value(name), Ok(ParamValue::String(text)) if text == "true");
    session.debug(
        out,
        &format!(
            "Value of string option 'test': {}",
            u8::from(flag(&param, "test"))
        ),
        1,
    )?;
    session.log.set_debug_level(cmd_int(&param, "debug", 0));
    session.debug(
        out,
        &format!(
            "Debug level (after ini file): {}",
            session.log.debug_level()
        ),
        1,
    )?;
    session.debug(
        out,
        &format!(
            "Value of string option 'no_progress': {}",
            u8::from(flag(&param, "no_progress"))
        ),
        1,
    )?;

    // 10. Parameter and file checks. The source runs them lazily as main_ reads
    //     each option (1388-1417, 1498-1614, 1968-2013); they run eagerly here.
    if let Some(code) = validate(session, spec, &mut param, out, err)? {
        return Ok(Prepared::Done(code));
    }
    Ok(Prepared::Run(Box::new(ToolContext::new(
        T::NAME,
        &session.version,
        &session.location,
        param,
        Arc::clone(&session.log),
    ))))
}

/// The four tool-description writers the Release build refuses without TDL,
/// in the order the source checks them after `-write_ctd`
/// (`TOPPBase.cpp:2651-2683`), with their file extensions.
const TDL_WRITERS: [(&str, &str); 4] = [
    ("write_nested_cwl", ".cwl"),
    ("write_cwl", ".cwl"),
    ("write_nested_json", ".json"),
    ("write_json", ".json"),
];

/// The message of the Release build's `std::runtime_error` from
/// `ParamCWLFile` and `ParamJSONFile` when TDL is not compiled in
/// (`ParamCWLFile.cpp:331`, `ParamJSONFile.cpp:326`).
pub const TDL_UNAVAILABLE: &str =
    "TDL support is not available. Rebuild with -DENABLE_TDL=ON to enable this feature.";

/// Source `handleWriteCommands_` (`TOPPBase.cpp:2546-2686`).
///
/// `-write_ini` writes the defaults, updated leniently from `-ini` when given,
/// and never includes other command-line values. The file declares
/// ISO-8859-1, as the source `ParamXMLFile::store` does
/// ([`paramxml::WriteOptions::source`]).
///
/// `-write_ctd <dir>` writes `<dir>/<tool><type>.ctd` for each of the tool's
/// registry types, or once without a type: the defaults (never `-ini` or
/// command-line values, but with the per-user defaults) and the tool's version,
/// name, documentation URL, registry category, description and citation DOIs,
/// through [`ParamCtdFile`]. An empty directory is the current directory.
///
/// `-write_nested_cwl`, `-write_cwl`, `-write_nested_json` and `-write_json`
/// check the target as the source does and then end with
/// [`ExitCode::InternalError`] and the Release build's message,
/// `Unable to initialize or run <tool>: TDL support is not available. …`.
/// The Release build opens (and so creates or empties) the target before it
/// fails; this port does not touch it.
fn write_commands<T: Tool>(
    session: &Session<'_>,
    cmd: &Param,
    defaults: &Param,
    out: &mut dyn Write,
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
            warn_if_not_applicable(session, &loaded, out, err)?;
            let diagnostics = update_diagnostics(&written, &loaded, UpdateMode::Lenient)?;
            written.update(&loaded, false)?;
            for line in &diagnostics {
                writeln!(err, "{line}")?;
            }
        }
        return Ok(Some(
            match paramxml::store_with_options(&path, &written, paramxml::WriteOptions::source()) {
                Ok(()) => ExitCode::ExecutionOk,
                Err(error) => run_failure(&error, err),
            },
        ));
    }
    let writers = std::iter::once(("write_ctd", ".ctd")).chain(TDL_WRITERS);
    for (name, extension) in writers {
        if !cmd.exists(name)? {
            continue;
        }
        let mut directory = text_value(cmd, name)?;
        if directory.is_empty() {
            directory = std::env::current_dir()?.display().to_string();
        }
        let mut types = session.registry.get_types(T::NAME)?;
        if types.is_empty() {
            types.push(String::new());
        }
        for kind in &types {
            let target = format!("{directory}/{}{kind}{extension}", T::NAME);
            if let Some(code) = output_file_writable(&target, name, err) {
                return Ok(Some(code));
            }
            if name != "write_ctd" {
                writeln!(
                    err,
                    "Unable to initialize or run {}: {TDL_UNAVAILABLE}",
                    T::NAME
                )?;
                return Ok(Some(ExitCode::InternalError));
            }
            let mut described = defaults.checked_clone()?;
            if !kind.is_empty() {
                described.set_value(
                    &format!("{}type", session.location),
                    ParamValue::String(kind.clone()),
                    "",
                    &[],
                )?;
            }
            let info = ToolInfo {
                version: session.version.clone(),
                name: T::NAME.to_owned(),
                docurl: usage::documentation_url(T::NAME),
                category: session.registry.get_category(T::NAME)?,
                description: T::DESCRIPTION.to_owned(),
                citations: std::iter::once(CITE_OPENMS.doi)
                    .chain(T::CITATIONS.iter().map(|citation| citation.doi))
                    .map(str::to_owned)
                    .collect(),
            };
            if let Err(error) = ParamCtdFile.store(&target, &described, &info) {
                // The source's store throws std::ios::failure, a
                // std::exception only the initialisation catch handles;
                // libstdc++ appends its category text to the message.
                writeln!(
                    err,
                    "Unable to initialize or run {}: {}: iostream error",
                    T::NAME,
                    error_text(&error)
                )?;
                return Ok(Some(ExitCode::InternalError));
            }
        }
        return Ok(Some(ExitCode::ExecutionOk));
    }
    Ok(None)
}

/// Load an INI file with the exit codes the source's run-phase catch assigns.
///
/// Both callers, the INI merge and `-write_ini`, load inside the source's
/// run-phase `try` (`TOPPBase.cpp:258`; the loads are at `296` and `2630`), so
/// a failure takes the inner catch:
///
/// * A missing file is [`ExitCode::InputFileNotFound`]: `XMLFile::parse_`
///   checks `File::exists` (caught at `436-441`).
/// * An existing file this process cannot read is
///   [`ExitCode::InputFileNotReadable`], with the source's `FileNotReadable`
///   wording (`448-453`). In the source, xerces cannot open the file, and
///   `XMLHandler::fatalError` asks `FileHandler::getTypeByContent` for a
///   file-type hint, whose `TextFile` load throws `FileNotReadable`
///   (`XMLHandler.cpp:49-50`, `FileHandler.cpp:402`, `TextFile.cpp:40-43`)
///   before the `ParseError` is raised. A regular file or a directory is
///   asked with `file::readable` before the load. Anything else, such as a
///   character device or a FIFO, is not asked, because that query answers
///   `false` for such a file without opening it; the load opens it instead,
///   and only a denied open takes this row, as for a FIFO with mode 000
///   (oracle `ini_fifo_denied` and `write_ini_ini_fifo_denied`).
/// * An existing file whose open fails for any other reason is
///   [`ExitCode::UnknownError`], with the source's
///   `Error: Unexpected internal error (IO error for file '<path>')`: the
///   `TextFile` load throws `IOException` when the file exists and
///   `File::readable` holds (`TextFile.cpp:44-47`), which takes the
///   `BaseException` arm (`495-499`). Examples are a Unix socket and `/dev/tty`
///   without a controlling terminal (oracle `ini_socket`, `ini_tty` and their
///   `write_ini_` counterparts). A load that fails with such an error is told
///   apart from a later read failure by opening the file once more; a FIFO is
///   never opened again, because that open would wait for a writer, so its
///   failure is taken as a read failure.
/// * A readable directory, malformed XML and any other read failure of an
///   existing file are [`ExitCode::InputFileCorrupt`], the source's
///   `ParseError` (`460-465`). For a directory the diagnostic is the source's
///   without the file-type hint that `XMLHandler::fatalError` appends. The
///   character device `/dev/null` is readable and reads as an empty document,
///   so it is a `ParseError` here as in the source (oracle `ini_dev_null` and
///   `write_ini_ini_dev_null`).
///
/// A FIFO this process can open is opened once and read like a file: opening
/// it waits for a writer, so with no writer the load blocks and no exit code is
/// reached, as for the C++ tool. The source opens the INI twice, once to look
/// for compression (`XMLFile.cpp:141-147`) and once for xerces (`166`), so a
/// writer that opens the FIFO only once leaves the C++ tool waiting for a
/// second writer, while this port reads what the one writer sends. This is a
/// deliberate difference.
///
/// The same mapping applies when the file changes between these checks and the
/// load. Failures other than I/O and parsing, such as a document beyond the
/// reader's limits, map as in `run_failure`.
fn load_ini(path: &str, err: &mut dyn Write) -> Result<std::result::Result<Param, ExitCode>> {
    let not_found = format!("Error: File not found (the file '{path}' does not exist)");
    let not_readable = format!(
        "Error: File not readable (the file '{path}' is not readable for the current user)"
    );
    if !file::exists(path) {
        writeln!(err, "{not_found}")?;
        return Ok(Err(ExitCode::InputFileNotFound));
    }
    // `file::readable` refuses a device or a FIFO without opening it, which
    // would report a readable `/dev/null` as unreadable; those are left to the
    // load, whose open failure is mapped below.
    let directory = file::is_directory(path);
    if (directory || std::path::Path::new(path).is_file()) && !file::readable(path) {
        writeln!(err, "{not_readable}")?;
        return Ok(Err(ExitCode::InputFileNotReadable));
    }
    if directory {
        writeln!(
            err,
            "Error: Unable to read file (While loading '{path}': unable to read data from file)"
        )?;
        return Ok(Err(ExitCode::InputFileCorrupt));
    }
    match paramxml::load(path) {
        Ok(loaded) => Ok(Ok(loaded)),
        Err(Error::Io(error)) => {
            use std::io::ErrorKind;
            let open_failure = match error.kind() {
                ErrorKind::NotFound | ErrorKind::PermissionDenied => Some(error.kind()),
                _ => reopen_failure(path),
            };
            let (code, text) = match open_failure {
                Some(ErrorKind::NotFound) => (ExitCode::InputFileNotFound, not_found),
                Some(ErrorKind::PermissionDenied) => (ExitCode::InputFileNotReadable, not_readable),
                Some(_) => (
                    ExitCode::UnknownError,
                    format!("Error: Unexpected internal error (IO error for file '{path}')"),
                ),
                None => (
                    ExitCode::InputFileCorrupt,
                    format!("Error: Unable to read file (While loading '{path}': {error})"),
                ),
            };
            writeln!(err, "{text}")?;
            Ok(Err(code))
        }
        Err(error) => Ok(Err(run_failure(&error, err))),
    }
}

/// The kind of error opening `path` again gives, or `None` when it opens.
///
/// `load_ini` asks this after a load failed with an I/O error, to tell a failed
/// open from a later read failure. A FIFO is not opened again, because the open
/// would wait for a writer; `None` is returned for it.
fn reopen_failure(path: &str) -> Option<std::io::ErrorKind> {
    #[cfg(unix)]
    {
        use std::os::unix::fs::FileTypeExt;
        if std::fs::metadata(path).is_ok_and(|metadata| metadata.file_type().is_fifo()) {
            return None;
        }
    }
    std::fs::File::open(path).err().map(|error| error.kind())
}

/// Source `checkIfIniParametersAreApplicable_` (`TOPPBase.cpp:1957-1966`), a
/// `writeLogWarn_`.
fn warn_if_not_applicable(
    session: &Session<'_>,
    ini: &Param,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<()> {
    let location = &session.location;
    if ini.copy(location, false)?.is_empty() {
        session.warn(
            out,
            err,
            &format!(
                "Warning: The provided INI file does not contain any parameters specific for this tool (expected in '{location}'). Please check your .ini file. The default parameters for this tool will be applied."
            ),
        )?;
    }
    Ok(())
}

/// The INI-version notice of `TOPPBase.cpp:355-366`, written to the output
/// stream and the log file as the source's `writeLogInfo_` does.
fn warn_on_version_mismatch(session: &Session<'_>, ini: &Param, out: &mut dyn Write) -> Result<()> {
    let key = format!("{}:version", session.name);
    if ini.exists(&key)? {
        let file_version = ini.value(&key)?.to_text(false)?;
        if file_version != session.version {
            session.info(
                out,
                &format!(
                    "Warning: Parameters file version ({file_version}) does not match the version of this tool ({}).\nYour current parameters are still valid, but there might be new valid values or even new parameters. Upgrading the INI might be useful.",
                    session.version
                ),
            )?;
        }
    }
    Ok(())
}

/// Source `checkParam_` (`TOPPBase.cpp:1872-1955`): warnings, never an error,
/// for one section of the INI file after the strict update succeeded.
///
/// An entry inside a subsection is not checked, only warned about when the
/// subsection is neither a subsection of the tool's own parameters nor the
/// first level of a registered algorithm subsection. Any other entry must
/// name a registered parameter and carry its value type. The source's
/// exemption for the tool's own name under `common::` never applies, because
/// the section with the tool's name is checked as `common:<tool>::`, with a
/// doubled colon; reproduced. So is the doubled space of the source's
/// `Wrong  parameter type` message for a floating-point parameter.
fn check_param(
    session: &Session<'_>,
    spec: &ToolSpec,
    param: &Param,
    filename: &str,
    location: &str,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<()> {
    let topp_subsections: BTreeSet<&str> = spec
        .topp_subsections()
        .iter()
        .map(|(name, _)| name.as_str())
        .collect();
    let subsections: BTreeSet<&str> = spec
        .subsections()
        .iter()
        .map(|(name, _)| name.as_str())
        .collect();
    for item in param.iter()? {
        let name = &item.key;
        let subsection = name.rfind(':').map_or("", |position| &name[..position]);
        if !subsection.is_empty() && !topp_subsections.contains(subsection) {
            let first = subsection.split(':').next().unwrap_or("");
            let exempt = location == "common::" && subsection == session.name;
            if !(subsections.contains(first) || exempt) {
                session.warn(
                    out,
                    err,
                    &format!(
                        "Warning: Unknown subsection '{subsection}' in '{filename}' (location '{location}')!"
                    ),
                )?;
            }
            continue;
        }
        let Some(entry) = spec.find(name) else {
            session.warn(
                out,
                err,
                &format!("Warning: Unknown parameter '{location}{name}' in '{filename}'!"),
            )?;
            continue;
        };
        let value = &item.entry.value;
        let expected = match entry.kind {
            ParameterType::String
            | ParameterType::InputFile
            | ParameterType::OutputFile
            | ParameterType::OutputPrefix
            | ParameterType::Flag => {
                (!matches!(value, ParamValue::String(_))).then_some(("", "string"))
            }
            ParameterType::Double => {
                (!matches!(value, ParamValue::Float(_))).then_some((" ", "double"))
            }
            ParameterType::Int => (!matches!(value, ParamValue::Integer(_))).then_some(("", "int")),
            ParameterType::StringList
            | ParameterType::InputFileList
            | ParameterType::OutputFileList => {
                (!matches!(value, ParamValue::StringList(_))).then_some(("", "string list"))
            }
            ParameterType::IntList => {
                (!matches!(value, ParamValue::IntegerList(_))).then_some(("", "int list"))
            }
            ParameterType::DoubleList => {
                (!matches!(value, ParamValue::FloatList(_))).then_some(("", "double list"))
            }
            _ => None,
        };
        if let Some((space, type_name)) = expected {
            session.warn(
                out,
                err,
                &format!(
                    "Warning: Wrong {space}parameter type of '{location}{name}' in '{filename}'. Type should be '{type_name}'!"
                ),
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

/// Source `inputFileReadable_` inside the lifecycle: the `Checking input file`
/// debug line, the heading on the error stream only (the source's
/// `OPENMS_LOG_ERROR`), and the catch block's `Error: …` line on the error
/// stream and in the log file.
fn check_input(
    session: &Session<'_>,
    path: &str,
    name: &str,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<Option<ExitCode>> {
    session.debug(out, &format!("Checking input file '{path}'"), 2)?;
    Ok(match context::input_file_problem(path, name) {
        Some((code, heading, detail)) => {
            writeln!(err, "{heading}")?;
            session.warn(out, err, &detail)?;
            Some(code)
        }
        None => None,
    })
}

/// Source `outputFileWritable_` inside the lifecycle, logged as
/// [`check_input`] is.
fn check_output(
    session: &Session<'_>,
    path: &str,
    name: &str,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<Option<ExitCode>> {
    session.debug(out, &format!("Checking output file '{path}'"), 2)?;
    Ok(match context::output_file_problem(path, name) {
        Some((code, heading, detail)) => {
            writeln!(err, "{heading}")?;
            session.warn(out, err, &detail)?;
            Some(code)
        }
        None => None,
    })
}

/// Required-value, restriction and file checks, in registration order. Returns
/// the terminal exit code when the run cannot proceed.
///
/// An input file tagged `is_executable` is resolved on `PATH` first and its
/// value replaced by the full path, as the source's `getStringOption_` returns
/// it (`TOPPBase.cpp:1534-1549`).
fn validate(
    session: &Session<'_>,
    spec: &ToolSpec,
    param: &mut Param,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<Option<ExitCode>> {
    for entry in spec.parameters() {
        if entry.kind.is_layout() || !param.exists(&entry.name)? {
            continue;
        }
        let value = param.value(&entry.name)?.checked_clone()?;
        // An empty list counts as missing, as getStringList_, getIntList_ and
        // getDoubleList_ throw RequiredParameterNotGiven for it (TOPPBase.cpp:1630, 1661, 1695).
        let empty = match &value {
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
            session.warn(
                out,
                err,
                &format!("Error: The required parameter {name} was not given or is empty!"),
            )?;
            return Ok(Some(ExitCode::MissingParameters));
        }
        if empty {
            continue;
        }
        let code = match (entry.kind, &value) {
            (ParameterType::String, ParamValue::String(text))
                if !entry.valid_strings.is_empty() && !entry.valid_strings.contains(text) =>
            {
                session.warn(
                    out,
                    err,
                    &format!(
                        "Invalid parameter: Invalid value '{text}' for string parameter '{}' given. Valid strings are: '{}'.",
                        entry.name,
                        entry.valid_strings.join("', '")
                    ),
                )?;
                Some(ExitCode::IllegalParameters)
            }
            (ParameterType::Int, ParamValue::Integer(number)) => {
                int_range(session, entry, *number, out, err)?
            }
            (ParameterType::IntList, ParamValue::IntegerList(numbers)) => {
                let mut code = None;
                for number in numbers {
                    code = int_range(session, entry, i64::from(*number), out, err)?;
                    if code.is_some() {
                        break;
                    }
                }
                code
            }
            (ParameterType::Double, ParamValue::Float(number)) => {
                float_range(session, entry, *number, out, err)?
            }
            (ParameterType::DoubleList, ParamValue::FloatList(numbers)) => {
                let mut code = None;
                for number in numbers {
                    code = float_range(session, entry, *number, out, err)?;
                    if code.is_some() {
                        break;
                    }
                }
                code
            }
            (ParameterType::InputFile, ParamValue::String(path)) => {
                if entry.tags.iter().any(|tag| tag == "is_executable") {
                    match executable_path(session, entry, path, out, err)? {
                        Ok(resolved) => {
                            // Keep the entry's description, tags and
                            // restrictions; only the value changes.
                            let mut replacement = param.entry(&entry.name)?.clone();
                            replacement.value = ParamValue::String(resolved.clone());
                            let prefix_length = entry.name.len() - replacement.name.len();
                            let prefix = entry.name[..prefix_length].to_owned();
                            param.insert_entry(replacement, &prefix)?;
                            input_path(session, entry, &resolved, out, err)?
                        }
                        Err(code) => Some(code),
                    }
                } else {
                    input_path(session, entry, path, out, err)?
                }
            }
            (ParameterType::InputFileList, ParamValue::StringList(paths)) => {
                let mut code = None;
                for path in paths {
                    code = input_path(session, entry, path, out, err)?;
                    if code.is_some() {
                        break;
                    }
                }
                code
            }
            (ParameterType::OutputFile, ParamValue::String(path)) => {
                match check_output(session, path, &entry.name, out, err)? {
                    Some(code) => Some(code),
                    None => output_extension(session, entry, path, out, err)?,
                }
            }
            (ParameterType::OutputPrefix, ParamValue::String(path)) => {
                check_output(session, &format!("{path}_0"), &entry.name, out, err)?
            }
            // The source checks neither the writability nor the format of an
            // output file list: its list validity check handles input file
            // lists only (TOPPBase.cpp:1498-1527).
            _ => None,
        };
        if code.is_some() {
            return Ok(code);
        }
    }
    Ok(None)
}

/// Resolve an input file tagged `is_executable` on `PATH`
/// (`TOPPBase.cpp:1534-1549`): a name that already names a file is kept, else
/// the first match on `PATH`. When nothing is found, the source's warning and
/// `ExternalExecutableNotFound`, [`ExitCode::ExternalProgramNotFound`].
fn executable_path(
    session: &Session<'_>,
    entry: &ParameterInformation,
    path: &str,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<std::result::Result<String, ExitCode>> {
    match session
        .registry
        .sources()
        .file_context
        .find_executable(path)?
    {
        Some(found) => {
            let resolved = found.display().to_string();
            session.debug(out, &format!("Input file resolved to '{resolved}'"), 2)?;
            Ok(Ok(resolved))
        }
        None => {
            let optional = if entry.required {
                ""
            } else {
                " Since this file is not strictly required, you might also pass the empty string \"\" as argument to prevent its usage (this might limit the usability of the tool)."
            };
            session.warn(
                out,
                err,
                &format!(
                    "Input file '{path}' could not be found (by searching on PATH). Either provide a full filepath via the '-{}' option or fix your PATH environment !{optional}",
                    entry.name
                ),
            )?;
            session.warn(
                out,
                err,
                &format!(
                    "Error: Executable not found (the executable '{path}' could not be found)"
                ),
            )?;
            Ok(Err(ExitCode::ExternalProgramNotFound))
        }
    }
}

/// Source `getIntOption_` range check (`TOPPBase.cpp:1486-1493`).
fn int_range(
    session: &Session<'_>,
    entry: &ParameterInformation,
    number: i64,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<Option<ExitCode>> {
    let min = i64::from(entry.min_int.unwrap_or(-i32::MAX));
    let max = i64::from(entry.max_int.unwrap_or(i32::MAX));
    if number < min || number > max {
        session.warn(
            out,
            err,
            &format!(
                "Invalid parameter: Invalid value '{number}' for integer parameter '{}' given. Out of valid range: '{min}'-'{max}'.",
                entry.name
            ),
        )?;
        return Ok(Some(ExitCode::IllegalParameters));
    }
    Ok(None)
}

/// Source `getDoubleOption_` range check (`TOPPBase.cpp:1459-1466`).
fn float_range(
    session: &Session<'_>,
    entry: &ParameterInformation,
    number: f64,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<Option<ExitCode>> {
    let min = entry.min_float.unwrap_or(-f64::MAX);
    let max = entry.max_float.unwrap_or(f64::MAX);
    if number < min || number > max {
        session.warn(
            out,
            err,
            &format!(
                "Invalid parameter: Invalid value '{number}' for float parameter '{}' given. Out of valid range: '{min}'-'{max}'.",
                entry.name
            ),
        )?;
        return Ok(Some(ExitCode::IllegalParameters));
    }
    Ok(None)
}

/// Whether `kind` is one of `formats`, compared by the source type name without
/// regard to ASCII case, as `ListUtils::contains(..., CASE::INSENSITIVE)`.
fn format_listed(formats: &[&str], kind: FileType) -> bool {
    formats
        .iter()
        .any(|format| format.eq_ignore_ascii_case(kind.name()))
}

/// Input readability, then the input format (`TOPPBase.cpp:1550`, `1575-1593`).
///
/// The format is what `FileHandler::get_type` detects: the file name first,
/// then bounded content recognition for a name it does not know. An
/// undetermined format only warns, `Warning: Could not determine format of
/// input file '<path>'!`, and the run continues; a detected format the
/// parameter does not accept is `InvalidParameter`, exit 6. A parameter without
/// registered formats is not checked.
///
/// `get_type` fails where content recognition cannot read the file, for
/// example on a directory with an unknown name, for which the source reports an
/// unknown type (an open follow-up of `src/format/file_handler.rs`). Such a
/// failure is treated here as an undetermined format, with the warning.
fn input_path(
    session: &Session<'_>,
    entry: &ParameterInformation,
    path: &str,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<Option<ExitCode>> {
    if !entry.tags.iter().any(|tag| tag == "skipexists") {
        if let Some(code) = check_input(session, path, &entry.name, out, err)? {
            return Ok(Some(code));
        }
    }
    let formats: Vec<&str> = entry.accepted_formats().collect();
    if formats.is_empty() {
        return Ok(None);
    }
    let kind = FileHandler::get_type(path).unwrap_or(FileType::Unknown);
    if kind == FileType::Unknown {
        session.warn(
            out,
            err,
            &format!("Warning: Could not determine format of input file '{path}'!"),
        )?;
    } else if !format_listed(&formats, kind) {
        session.warn(
            out,
            err,
            &format!(
                "Invalid parameter: Input file '{path}' has invalid format '{}'. Valid formats are: '{}'.",
                kind.name(),
                formats.join("','")
            ),
        )?;
        return Ok(Some(ExitCode::IllegalParameters));
    }
    Ok(None)
}

/// The output format by file name (`TOPPBase.cpp:1595-1608`).
///
/// An extension no file type claims is accepted, and so is any name when the
/// parameter has no registered formats; a known type the parameter does not
/// accept is `InvalidParameter`, exit 6.
fn output_extension(
    session: &Session<'_>,
    entry: &ParameterInformation,
    path: &str,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<Option<ExitCode>> {
    let formats: Vec<&str> = entry.accepted_formats().collect();
    if formats.is_empty() {
        return Ok(None);
    }
    let kind = type_by_file_name(path);
    if kind != FileType::Unknown && !format_listed(&formats, kind) {
        session.warn(
            out,
            err,
            &format!(
                "Invalid parameter: Invalid output file extension for file '{path}'. Valid file extensions are: '{}'.",
                formats.join("','")
            ),
        )?;
        return Ok(Some(ExitCode::IllegalParameters));
    }
    Ok(None)
}

/// The text of an error without this crate's variant prefix: the source
/// exception's `what()`, which the lifecycle prints after `Unable to
/// initialize or run <tool>: `. Registry and writer errors carry the source's
/// text in their message.
fn error_text(error: &Error) -> String {
    match error {
        Error::InvalidValue(message)
        | Error::InvalidRange(message)
        | Error::Unsupported(message)
        | Error::MissingInformation(message) => message.clone(),
        Error::Parse { message, .. } => message.clone(),
        Error::Io(error) => error.to_string(),
        Error::UnsortedData => error.to_string(),
    }
}

/// A failure before the tool body that no phase handles itself: the source's
/// initialisation catch for an OpenMS exception (`TOPPBase.cpp:505-508`),
/// `Unable to initialize or run <tool>: <what>` and
/// [`ExitCode::IllegalParameters`]. The message is the error's own text, the
/// source exception's `what()`, for the registry failures the Release build
/// reports this way (`../oracle/toppbase-completion`, cases `reg_*` and
/// `ttd_*`).
fn initialisation_failure<T: Tool>(error: &Error, err: &mut dyn Write) -> ExitCode {
    let _ = writeln!(
        err,
        "Unable to initialize or run {}: {}",
        T::NAME,
        error_text(error)
    );
    ExitCode::IllegalParameters
}

/// A failure while the tool runs, mapped as the source's run-phase catch
/// (`TOPPBase.cpp:430-499`).
///
/// A parse failure is `INPUT_FILE_CORRUPT`, as the source's `ParseError`
/// (460-465). The other arms are native mappings, because [`Error`] is coarser
/// than the source's exceptions:
///
/// * A missing file is `INPUT_FILE_NOT_FOUND`, as `FileNotFound` (436-441).
/// * A permission failure is `CANNOT_WRITE_OUTPUT_FILE`, as `UnableToCreateFile`
///   (430-435). A `std::io::Error` does not say whether it read or wrote; a
///   tool's inputs are checked for readability before its body runs, so a
///   denied permission there is taken as a failed write. INI files never reach
///   this arm: `load_ini` maps their read failures itself.
/// * Any other I/O failure is `UNKNOWN_ERROR`, the `BaseException` arm
///   (495-499).
/// * An invalid value or range is `ILLEGAL_PARAMETERS`, as `InvalidParameter`
///   (475-480), and missing information is `MISSING_PARAMETERS`, as
///   `RequiredParameterNotGiven` (466-474). The source's own `InvalidValue`,
///   `InvalidRange` and `MissingInformation` exceptions derive directly from
///   `BaseException` and would exit `UNKNOWN_ERROR` there.
/// * `Unsupported` and `UnsortedData` are `INCOMPATIBLE_INPUT_DATA`, the code
///   the source tools return explicitly for those conditions.
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

/// A failure raised while an output file is written, mapped as the source's
/// `UnableToCreateFile` arm (`TOPPBase.cpp:430-435`): the diagnostic
/// `Error: Unable to write file (<what>)` and
/// [`ExitCode::CannotWriteOutputFile`].
///
/// [`run_failure`] cannot tell a read from a write, because [`Error`] does not
/// say which side raised it, and its `Error::Parse` arm is the source's
/// `ParseError`, which is a *read* failure. A store that fails therefore has to
/// be mapped where it is known to be a store; otherwise a writer's refusal is
/// announced as `Error: Unable to read file (parse error on line 0: …)` with
/// `INPUT_FILE_CORRUPT`, which is what the featureXML store did before native
/// difference 16 was closed.
///
/// Every failure takes one arm, because the source's write side has one: a
/// store reaches `XMLFile::save_`, whose only exception is
/// `UnableToCreateFile` for a stream it cannot open (`XMLFile.cpp:366-372`),
/// and `FeatureXMLFile::store` raises the same one for a name whose extension
/// it does not accept (`FeatureXMLFile.cpp:75-78`). A native refusal this port
/// adds — an exceeded writer ceiling, a duplicate assigned feature ID, a field
/// this dialect cannot represent — is likewise a failure to produce the output
/// file and belongs in the same arm.
///
/// The text in the parentheses is the source's `UnableToCreateFile::what()`
/// whenever the port failed for the same reason the source would have, which
/// is any I/O failure of the store; `path` names the file it was producing. A
/// native refusal has no source counterpart, so it contributes its own message
/// there, as every other arm of [`run_failure`] does.
///
/// Only the `-out` store of `FeatureFinderCentroided` is mapped here so far,
/// because that is the store this port has executed C++ evidence for
/// (`../oracle/featurexml-inf`: with `-out` an existing directory whose name
/// carries the featureXML extension, the writability check passes and the
/// Release tool prints this line and exits 5); every other tool's store still
/// reaches [`run_failure`].
///
/// The source reaches its catch block by unwinding past the closing
/// `<tool> took …` line, so a run ended here prints none: this marks the run
/// as unwound for [`run_with`], on the calling thread.
pub(crate) fn write_failure(path: &str, error: &Error, err: &mut dyn Write) -> ExitCode {
    let detail = match error {
        Error::Io(_) => format!("the file '{path}' could not be created. "),
        other => other.to_string(),
    };
    let _ = writeln!(err, "Error: Unable to write file ({detail})");
    UNWOUND.with(|unwound| unwound.set(true));
    ExitCode::CannotWriteOutputFile
}

thread_local! {
    /// Whether the tool body on this thread ended through
    /// [`write_failure`], the stand-in for an exception the source's catch
    /// block handles.
    static UNWOUND: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Run a tool against explicit arguments and streams. `arguments[0]` is the
/// executable name, as in `main(argc, argv)`.
///
/// The tool registry is read from the process environment
/// ([`ToolHandler::from_environment`]), as the source reads it; see
/// [`run_with_registry`] to pass one. The phases and their exit codes follow
/// `TOPPBase::main`; see the module documentation. Usage text goes to `err` in
/// every case, for `--help` as after a command-line error, because the source's
/// `printUsage_` writes to standard error; so does every diagnostic. `out`
/// receives what the source writes through its info log: the INI-version
/// notice, a tool's report and the closing `<tool> took … .` line.
///
/// One phase precedes all of the source's: this build may require processor
/// features the source's does not, and
/// [`system::cpu_features::unsupported_cpu`](crate::system::cpu_features::unsupported_cpu)
/// is asked before anything else happens. On a processor that cannot run this
/// binary the message goes to `err` — through [`Write::write_all`], not the
/// formatting machinery — and the status is [`ExitCode::InternalError`], which
/// is a diagnosis rather than the `SIGILL` the first fused multiply-add would
/// otherwise raise. Nothing before that point does floating-point arithmetic.
/// On a build without the flag the whole phase is compiled out.
///
/// [`run`] asks the same question one step earlier still, before it reads the
/// command line; this copy covers a caller that drives a tool in process.
pub fn run_with<T: Tool>(
    arguments: &[String],
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> ExitCode {
    if let Some(message) = crate::system::cpu_features::unsupported_cpu() {
        let _ = err.write_all(message.as_bytes());
        let _ = err.write_all(b"\n");
        return ExitCode::InternalError;
    }
    let registry = match ToolHandler::from_environment() {
        Ok(registry) => registry,
        Err(error) => return initialisation_failure::<T>(&error, err),
    };
    run_with_registry::<T>(arguments, out, err, &registry)
}

/// [`run_with`] against an explicit tool registry.
///
/// The registry decides the product version the tool reports, whether its
/// common options are `TOPP` or `UTIL` options, its category and `-type`
/// values for `-write_ctd`, and where its per-user defaults are read from. A
/// registry that cannot be read ends the run before anything else, with the
/// source's `Unable to initialize or run <tool>: <what>` and
/// [`ExitCode::IllegalParameters`] (`TOPPBase.cpp:141-163`, `505-508`).
pub fn run_with_registry<T: Tool>(
    arguments: &[String],
    out: &mut dyn Write,
    err: &mut dyn Write,
    registry: &ToolHandler,
) -> ExitCode {
    if let Some(message) = crate::system::cpu_features::unsupported_cpu() {
        let _ = err.write_all(message.as_bytes());
        let _ = err.write_all(b"\n");
        return ExitCode::InternalError;
    }
    let (version, verbose_version) = match product_versions(registry, T::NAME) {
        Ok(versions) => versions,
        Err(error) => return initialisation_failure::<T>(&error, err),
    };
    let registered = registry
        .get_topp_tool_list()
        .map(|tools| tools.contains_key(T::NAME));
    for line in registry.take_diagnostics() {
        let _ = writeln!(err, "{line}");
    }
    let registered = match registered {
        Ok(registered) => registered,
        Err(error) => return initialisation_failure::<T>(&error, err),
    };
    let spec = match tool_spec_for::<T>(registered) {
        Ok(spec) => spec,
        Err(error) => return initialisation_failure::<T>(&error, err),
    };
    let mut session = Session {
        registry,
        name: T::NAME,
        version,
        verbose_version,
        instance: 1,
        location: location_of(T::NAME, 1),
        log: Arc::new(ToolLog::new()),
    };
    let code = match prepare::<T>(&mut session, &spec, arguments, out, err) {
        Ok(Prepared::Run(ctx)) => run_body::<T>(&session, &ctx, out, err),
        Ok(Prepared::Done(code)) => code,
        Err(error) => initialisation_failure::<T>(&error, err),
    };
    session.log.finish();
    code
}

/// The tool body with the source's timing and failure mapping
/// (`TOPPBase.cpp:413-500`): when the body returns, the closing line
/// `<tool> took <wall> (wall), <cpu> (CPU), <system> (system), <user> (user);
/// Peak Memory Usage: <n> MB.` on the output stream, the peak-memory part only
/// where the platform reports it; when it fails — with an error, or through
/// [`write_failure`] — no closing line, but the run-phase mapping, whose
/// `Error` lines also reach the log file.
fn run_body<T: Tool>(
    session: &Session<'_>,
    ctx: &ToolContext,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> ExitCode {
    let mut watch = crate::system::stop_watch::StopWatch::new();
    UNWOUND.with(|unwound| unwound.set(false));
    let _ = watch.start();
    let result = T::run_io(ctx, out, err);
    let _ = watch.stop();
    match result {
        Ok(code) if UNWOUND.with(|unwound| unwound.replace(false)) => code,
        Ok(code) => {
            let memory = crate::system::sys_info::process_peak_memory_consumption()
                .filter(|&kib| kib != 0)
                .map(|kib| format!("; Peak Memory Usage: {} MB", kib / 1024))
                .unwrap_or_default();
            let _ = writeln!(out, "{} took {}{memory}.", T::NAME, watch.summary());
            code
        }
        Err(error) => {
            let mut text = Vec::new();
            let code = run_failure(&error, &mut text);
            let text = String::from_utf8_lossy(&text);
            for line in text.lines() {
                let _ = session.warn(out, err, line);
            }
            code
        }
    }
}

/// Run a tool against the process arguments and standard streams, returning the
/// status the executable should exit with.
///
/// This is what every ported tool executable calls. The processor check that
/// [`run_with`] describes is the **first statement of this function** — ahead of
/// reading the command line, because on an x86_64 build with
/// `-C target-feature=+fma` even that much is compiled with VEX encoding and
/// would fault on a processor old enough to lack AVX. Everything after the check
/// is in `run_from_environment`, which is never inlined, so no instruction of it
/// can be hoisted above the check. [`run_with`] holds a second copy for a caller
/// that drives a tool in process; no tool binary reaches the guard through that
/// copy. See `docs/FMA_BUILD_FLAG.md` section 8 for what this does and does not
/// guarantee.
pub fn run<T: Tool>() -> ExitCode {
    if let Some(message) = crate::system::cpu_features::unsupported_cpu() {
        return report_unsupported_cpu(message);
    }
    run_from_environment::<T>()
}

/// Write the processor refusal to standard error and give the exit status.
///
/// Kept as small and as plain as it can be: the message is a constant, it goes
/// out through [`Write::write_all`] on a locked [`std::io::Stderr`] rather than
/// through the formatting machinery, and nothing here allocates or computes.
/// This code runs on a processor that cannot execute everything the build
/// emitted, so the less of the build it uses, the better.
#[cold]
#[inline(never)]
fn report_unsupported_cpu(message: &'static str) -> ExitCode {
    let mut err = std::io::stderr().lock();
    let _ = err.write_all(message.as_bytes());
    let _ = err.write_all(b"\n");
    ExitCode::InternalError
}

/// The body of [`run`]: the process arguments and the standard streams.
///
/// Separate and `#[inline(never)]` only so that [`run`]'s processor check comes
/// first in the emitted code as well as in the source.
#[inline(never)]
fn run_from_environment<T: Tool>() -> ExitCode {
    let arguments: Vec<String> = std::env::args().collect();
    let mut out = std::io::stdout();
    let mut err = std::io::stderr();
    run_with::<T>(&arguments, &mut out, &mut err)
}
