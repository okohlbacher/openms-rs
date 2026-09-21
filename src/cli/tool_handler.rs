// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The registry of installed tool packages and legacy internal tools: the
//! native form of the `OpenMS4-cli` package's `APPLICATIONS/ToolHandler.h` and
//! `ToolHandler.cpp` (cli `c19e494`).
//!
//! Packages install tab-separated manifests, `share/openms4/tools/*.tsv`, under
//! their installation prefix. Each row holds a tool's name, category, product
//! version and executable path relative to that prefix. `TOPPBase` asks the
//! registry for the product version it prints and writes, for the category of
//! the tool descriptions, and whether a tool is a registered TOPP tool at all;
//! see [`ToolHandler`] for the lookups and `docs/TOPP_CLI_SUPPORT.md`
//! (*Tool registry*) for the member mapping.
//!
//! **Where manifests are found.** As the source (`ToolHandler.cpp:38-85`):
//! every non-empty entry of `OPENMS_TOOL_PREFIX_PATH`, split at the platform's
//! path-list separator, then the prefix of the running executable
//! (`<executable directory>/..`) and, on macOS, the prefix of an application
//! bundle the executable sits in. The source also adds the prefix of the
//! loaded `libOpenMS_CLI` library through `dladdr`; this crate is linked into
//! the executable, so that is the executable's prefix again and adds nothing.
//!
//! **The built-in product manifest** is the one native addition. A Rust tool
//! is not installed by the source's CMake rules, so no manifest sits beside
//! it. This crate therefore carries the product manifest the pinned `topp`
//! package installs — `resources/tools/topp.tools.tsv`, byte-identical to the
//! Release build's `share/openms4/tools/topp.tools.tsv` (sha256 `5a90f7c1…`) —
//! and reads it **in place of the executable's own prefix** whenever that
//! prefix has no `share/openms4/tools` directory. An executable installed into
//! a prefix that does carry manifests reads those instead, exactly as the
//! source.
//!
//! **The same manifest found twice is read once**, as in the source, whose
//! `seen` set skips a manifest file it has already read
//! (`ToolHandler.cpp:90-104`; the Release build with its own prefix named in
//! `OPENMS_TOOL_PREFIX_PATH`, even twice, runs normally, oracle
//! `reg_install_prefix_twice`). The source knows a manifest by its canonical
//! path, which the built-in manifest does not have; it is known by what it is,
//! the file `topp.tools.tsv` the pinned `topp` package installs. So when a
//! prefix supplies `share/openms4/tools/topp.tools.tsv` byte for byte — an
//! `OPENMS_TOOL_PREFIX_PATH` naming the C++ installation, which is what the
//! variable is for — that manifest is the built-in one reached a second time,
//! and the built-in copy is not read again. A manifest under another name, or
//! with other bytes (another product version), that lists a product tool is
//! a second manifest and a duplicate, as it is for two C++ installations
//! (oracles `reg_dup_*`, and `reg_copy_of_install_manifest` of
//! `../oracle/topp-exception-exits`, where a copy of the installed manifest
//! under another prefix is a duplicate because the source's identity is the
//! path); the source reports that for every tool run.
//!
//! Manifests are read on demand, on every lookup, as the source's are. The
//! legacy internal-tool registry (`.ttd` files) is read once per
//! [`ToolHandler`] and kept after its first successful read, as the source's
//! function-local static is.

use super::tool_description_file::ToolDescriptionFile;
use crate::data_structures::ToolDescription;
use crate::system::file::{self, FileContext};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, BufReader};
use std::path::{Component, Path, PathBuf};
use std::sync::{Mutex, OnceLock};

/// Tool name to its description, ordered by name (source `ToolListType`, a
/// `std::map<std::string, Internal::ToolDescription>`).
pub type ToolListType = BTreeMap<String, ToolDescription>;

/// The product manifest the pinned `topp` package installs, as the Release
/// build ships it (`share/openms4/tools/topp.tools.tsv`).
pub const BUILTIN_MANIFEST: &str = include_str!("../../resources/tools/topp.tools.tsv");

/// The file name the built-in manifest stands in for.
pub const BUILTIN_MANIFEST_NAME: &str = "topp.tools.tsv";

/// Most manifest rows read in one lookup, over all manifests.
pub const MAX_MANIFEST_ROWS: usize = 100_000;

/// Most bytes of one manifest file.
pub const MAX_MANIFEST_BYTES: u64 = 16 * 1024 * 1024;

/// Categories of interactive desktop applications, which the tool list leaves
/// out (`ToolHandler.cpp:158`) while their versions and executables stay
/// resolvable.
const DESKTOP_CATEGORIES: [&str; 2] = ["DesktopViewer", "DesktopWorkflow"];

/// The platform's path-list separator for `OPENMS_TOOL_PREFIX_PATH`.
const PREFIX_SEPARATOR: char = if cfg!(windows) { ';' } else { ':' };

/// One manifest row (source `PackageTool`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PackageTool {
    /// Category, which may be empty.
    pub category: String,
    /// Product version.
    pub version: String,
    /// Absolute executable path: the row's relative path under its prefix.
    pub executable: PathBuf,
}

/// The source's `InvalidValue::what()`: `the value '<value>' was used but is not
/// valid; <message>`.
fn invalid_value(message: &str, value: &str) -> Error {
    Error::InvalidValue(format!(
        "the value '{value}' was used but is not valid; {message}"
    ))
}

/// The source's `FileNotReadable::what()`.
fn not_readable(path: &Path) -> Error {
    Error::Io(std::io::Error::new(
        std::io::ErrorKind::PermissionDenied,
        format!(
            "the file '{}' is not readable for the current user",
            path.display()
        ),
    ))
}

/// The source's `FileNotFound::what()`.
fn not_found(name: &str) -> Error {
    Error::Io(std::io::Error::new(
        std::io::ErrorKind::NotFound,
        format!("the file '{name}' could not be found"),
    ))
}

/// An operating-system error's text without Rust's ` (os error N)` suffix, as
/// `strerror` gives it.
fn os_text(error: &std::io::Error) -> String {
    let text = error.to_string();
    match text.rfind(" (os error ") {
        Some(position) if text.ends_with(')') => text[..position].to_owned(),
        _ => text,
    }
}

/// The source's `filesystem_error::what()` for a failed directory operation,
/// as libstdc++ formats it: `filesystem error: <what>: <strerror> [<path>]`.
fn filesystem_error(what: &str, error: &std::io::Error, path: &Path) -> Error {
    invalid_value(
        "Cannot read tool package registry",
        &format!(
            "filesystem error: {what}: {} [{}]",
            os_text(error),
            path.display()
        ),
    )
}

/// `path` with `.` removed and `name/..` folded, as
/// `std::filesystem::path::lexically_normal` does for the prefixes built here.
fn lexically_normal(path: &Path) -> PathBuf {
    let mut result = PathBuf::new();
    for component in path.components() {
        match component {
            Component::CurDir => {}
            Component::ParentDir => {
                let last = result.components().next_back();
                match last {
                    Some(Component::Normal(_)) => {
                        result.pop();
                    }
                    Some(Component::RootDir | Component::Prefix(_)) => {}
                    _ => result.push(".."),
                }
            }
            other => result.push(other.as_os_str()),
        }
    }
    if result.as_os_str().is_empty() {
        result.push(".");
    }
    result
}

/// Where the registry reads from: every input the source takes from the
/// environment, explicit, so a lookup can be driven without changing the
/// process environment.
#[derive(Clone, Debug)]
pub struct ToolRegistrySources {
    /// The entries of `OPENMS_TOOL_PREFIX_PATH`, in order; empty entries are
    /// dropped by [`from_prefix_path`](Self::from_prefix_path).
    pub prefixes: Vec<PathBuf>,
    /// The directory holding the running executable (source
    /// `File::getExecutablePath`). Its parent is the installation prefix.
    pub executable_directory: PathBuf,
    /// Read [`BUILTIN_MANIFEST`] for the executable's prefix when that prefix
    /// has no `share/openms4/tools` directory, and no prefix installs the
    /// same manifest byte for byte as `topp.tools.tsv`. Default `true`.
    pub builtin_manifest: bool,
    /// The internal-tool configuration directory (source
    /// `getInternalToolsPath`: the shared-data directory plus
    /// `TOOLS/INTERNAL`), or `None` when no shared-data directory resolves.
    pub internal_tools_path: Option<PathBuf>,
    /// `OPENMS_TTD_INTERNAL_PATH`, an additional directory of `.ttd` files.
    pub ttd_internal_path: Option<PathBuf>,
    /// The environment snapshot used by [`ToolHandler::find_executable`] for
    /// the `PATH` search.
    pub file_context: FileContext,
}

impl ToolRegistrySources {
    /// Split an `OPENMS_TOOL_PREFIX_PATH` value at the platform separator
    /// (`;` on Windows, `:` elsewhere), dropping empty items as the source
    /// does (`ToolHandler.cpp:46-53`).
    pub fn from_prefix_path(value: &str) -> Vec<PathBuf> {
        value
            .split(PREFIX_SEPARATOR)
            .filter(|item| !item.is_empty())
            .map(PathBuf::from)
            .collect()
    }

    /// Snapshot the process environment the source reads:
    /// `OPENMS_TOOL_PREFIX_PATH`, `OPENMS_TTD_INTERNAL_PATH`, the executable
    /// path and everything [`FileContext::from_environment`] reads.
    ///
    /// The internal-tool directory is left unset when no shared-data directory
    /// resolves. The source's `File::getOpenMSDataPath` throws there, which
    /// would stop every tool; this crate's tools do not need a shared-data tree,
    /// so a missing one only means there is no internal-tool registry to read.
    ///
    /// # Errors
    ///
    /// As [`FileContext::from_environment`], and [`Error::InvalidValue`] when
    /// `OPENMS_TOOL_PREFIX_PATH` is not valid UTF-8.
    pub fn from_environment() -> Result<Self> {
        let context = FileContext::from_environment()?;
        let prefixes = match std::env::var_os("OPENMS_TOOL_PREFIX_PATH") {
            Some(value) => Self::from_prefix_path(value.to_str().ok_or_else(|| {
                Error::InvalidValue("OPENMS_TOOL_PREFIX_PATH is not valid UTF-8".into())
            })?),
            None => Vec::new(),
        };
        let internal_tools_path = context
            .get_openms_data_path()
            .ok()
            .map(|data| data.join("TOOLS/INTERNAL"));
        Ok(Self {
            prefixes,
            executable_directory: context.executable_directory.clone(),
            builtin_manifest: true,
            internal_tools_path,
            ttd_internal_path: std::env::var_os("OPENMS_TTD_INTERNAL_PATH").map(PathBuf::from),
            file_context: context,
        })
    }

    /// Sources that read nothing but the given prefixes: no executable prefix,
    /// no built-in manifest and no internal tools. For tests and embedders.
    ///
    /// # Errors
    ///
    /// As [`FileContext::new`].
    pub fn only_prefixes(prefixes: Vec<PathBuf>) -> Result<Self> {
        let nowhere = PathBuf::from("/nonexistent-openms-executable-directory/bin");
        Ok(Self {
            prefixes,
            executable_directory: nowhere.clone(),
            builtin_manifest: false,
            internal_tools_path: None,
            ttd_internal_path: None,
            file_context: FileContext::new(&nowhere, &nowhere, std::env::temp_dir())?,
        })
    }
}

/// One manifest source: a file under a prefix, or the built-in manifest
/// standing in for the executable prefix.
enum Manifest {
    File(PathBuf),
    Builtin(PathBuf),
}

/// The registry of installed tools (source class `ToolHandler`, whose members
/// are all static).
///
/// The source reads the process environment on every call and keeps the
/// internal tools in a function-local static. This is a caller-owned value
/// built from explicit [`ToolRegistrySources`]; [`ToolHandler::from_environment`]
/// takes them from the process. Lookups take `&self` and are safe to share
/// between threads.
#[derive(Debug)]
pub struct ToolHandler {
    sources: ToolRegistrySources,
    internal_tools: OnceLock<Vec<ToolDescription>>,
    diagnostics: Mutex<Vec<String>>,
}

impl ToolHandler {
    /// A registry reading from `sources`.
    pub fn new(sources: ToolRegistrySources) -> Self {
        Self {
            sources,
            internal_tools: OnceLock::new(),
            diagnostics: Mutex::new(Vec::new()),
        }
    }

    /// A registry reading the process environment, as the source's statics do.
    ///
    /// # Errors
    ///
    /// As [`ToolRegistrySources::from_environment`].
    pub fn from_environment() -> Result<Self> {
        Ok(Self::new(ToolRegistrySources::from_environment()?))
    }

    /// The sources this registry reads.
    pub fn sources(&self) -> &ToolRegistrySources {
        &self.sources
    }

    /// Non-fatal diagnostics raised while the internal-tool registry was read
    /// (source `XMLHandler::error(LOAD, …)` lines, `Non-fatal error while
    /// loading '<file>': <message>`), taken out of the registry.
    ///
    /// The source writes them to the error log as it reads; a caller writes
    /// these to its error stream.
    pub fn take_diagnostics(&self) -> Vec<String> {
        match self.diagnostics.lock() {
            Ok(mut lines) => std::mem::take(&mut *lines),
            Err(poisoned) => std::mem::take(&mut *poisoned.into_inner()),
        }
    }

    /// The prefixes searched, in the source's order.
    fn prefixes(&self) -> Vec<(PathBuf, bool)> {
        let mut result: Vec<(PathBuf, bool)> = self
            .sources
            .prefixes
            .iter()
            .map(|p| (p.clone(), false))
            .collect();
        let executable = self.sources.executable_directory.join(".");
        result.push((lexically_normal(&executable.join("..")), true));
        if cfg!(target_os = "macos") {
            // A desktop bundle executable lives at
            // <prefix>/bin/App.app/Contents/MacOS; its manifest is under
            // <prefix>/share, outside the bundle (ToolHandler.cpp:56-63).
            let bundle = lexically_normal(&executable.join("../.."));
            if bundle.extension().is_some_and(|e| e == "app") {
                result.push((lexically_normal(&executable.join("../../../..")), false));
            }
        }
        result
    }

    /// The manifests to read, in order, with the source's checks on each
    /// prefix directory (`ToolHandler.cpp:91-104`).
    fn manifests(&self) -> Result<Vec<Manifest>> {
        let mut result = Vec::new();
        let mut seen: BTreeSet<PathBuf> = BTreeSet::new();
        for (prefix, is_executable_prefix) in self.prefixes() {
            let directory = prefix.join("share/openms4/tools");
            let exists = match std::fs::metadata(&directory) {
                Ok(metadata) => metadata.is_dir(),
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => false,
                Err(error) => {
                    return Err(filesystem_error(
                        "Cannot inspect tool manifest directory",
                        &error,
                        &directory,
                    ));
                }
            };
            if !exists {
                if is_executable_prefix && self.sources.builtin_manifest {
                    result.push(Manifest::Builtin(directory.join(BUILTIN_MANIFEST_NAME)));
                }
                continue;
            }
            // The source iterates the directory in the order the file system
            // reports; this reads the manifests ordered by name, so which of two
            // duplicate rows is reported does not depend on the file system.
            let entries = std::fs::read_dir(&directory).map_err(|error| {
                filesystem_error(
                    "directory iterator cannot open directory",
                    &error,
                    &directory,
                )
            })?;
            let mut paths = Vec::new();
            for entry in entries {
                let entry = entry.map_err(|error| {
                    filesystem_error("directory iterator cannot advance", &error, &directory)
                })?;
                paths.push(entry.path());
            }
            paths.sort();
            for path in paths {
                if path.extension().is_none_or(|e| e != "tsv") || !path.is_file() {
                    continue;
                }
                // weakly_canonical: the same manifest reached twice is read once.
                let key = std::fs::canonicalize(&path).unwrap_or_else(|_| path.clone());
                if seen.insert(key) {
                    result.push(Manifest::File(path));
                }
            }
        }
        // The built-in manifest reached a second time, through a prefix that
        // installs it: read once, as the `seen` set above reads a file once
        // (see the module documentation).
        let installed = result.iter().any(|manifest| match manifest {
            Manifest::File(path) => is_builtin_manifest(path),
            Manifest::Builtin(_) => false,
        });
        if installed {
            result.retain(|manifest| matches!(manifest, Manifest::File(_)));
        }
        Ok(result)
    }

    /// Every manifest row by tool name (source `packageTools`,
    /// `ToolHandler.cpp:87-148`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] with the source's message, `the value '<file>:
    /// <line>' was used but is not valid; Invalid or duplicate tool package
    /// manifest entry`, for a row without exactly four tab-separated fields,
    /// with an empty name, version or executable, a name holding `/` or `\`,
    /// an executable path that is absolute or climbs with `..`, or a name
    /// already registered; the same variant with `Cannot read tool package
    /// registry` when a prefix directory cannot be inspected or listed; and
    /// [`Error::Io`] with the source's `FileNotReadable` text when a manifest
    /// cannot be opened or read. A manifest larger than
    /// [`MAX_MANIFEST_BYTES`], or more than [`MAX_MANIFEST_ROWS`] rows in all,
    /// is [`Error::InvalidValue`].
    pub fn package_tools(&self) -> Result<BTreeMap<String, PackageTool>> {
        let mut result: BTreeMap<String, PackageTool> = BTreeMap::new();
        let mut rows = 0usize;
        for manifest in self.manifests()? {
            let (path, text) = match &manifest {
                Manifest::Builtin(path) => (path.clone(), BUILTIN_MANIFEST.to_owned()),
                Manifest::File(path) => (path.clone(), read_manifest(path)?),
            };
            let prefix = path
                .parent()
                .and_then(Path::parent)
                .and_then(Path::parent)
                .and_then(Path::parent)
                .map(Path::to_path_buf)
                .unwrap_or_default();
            for raw in text.split_inclusive('\n') {
                let line = raw.strip_suffix('\n').unwrap_or(raw);
                let line = line.strip_suffix('\r').unwrap_or(line);
                if line.is_empty() || line.starts_with('#') {
                    continue;
                }
                rows += 1;
                if rows > MAX_MANIFEST_ROWS {
                    return Err(Error::InvalidValue(format!(
                        "tool package manifests hold more than {MAX_MANIFEST_ROWS} rows"
                    )));
                }
                let fail = || {
                    invalid_value(
                        "Invalid or duplicate tool package manifest entry",
                        &format!("{}: {line}", path.display()),
                    )
                };
                // std::getline splitting drops a trailing empty field, so the
                // separator count is checked separately, as the source does.
                let values: Vec<&str> = line.split('\t').collect();
                if line.matches('\t').count() != 3
                    || values.len() != 4
                    || values[0].is_empty()
                    || values[2].is_empty()
                    || values[3].is_empty()
                {
                    return Err(fail());
                }
                let relative = Path::new(values[3]);
                if relative.has_root() || values[0].contains(['/', '\\']) {
                    return Err(fail());
                }
                if matches!(relative.components().next(), Some(Component::Prefix(_))) {
                    return Err(fail());
                }
                if relative.components().any(|c| c == Component::ParentDir) {
                    return Err(fail());
                }
                let executable = std::path::absolute(prefix.join(relative)).map_err(Error::Io)?;
                if result.contains_key(values[0]) {
                    return Err(fail());
                }
                result.insert(
                    values[0].to_owned(),
                    PackageTool {
                        category: values[1].to_owned(),
                        version: values[2].to_owned(),
                        executable,
                    },
                );
            }
        }
        Ok(result)
    }

    /// Installed tools merged with the legacy internal-tool registry (source
    /// `getTOPPToolList`, `ToolHandler.cpp:153-178`).
    ///
    /// Desktop applications (categories `DesktopViewer` and
    /// `DesktopWorkflow`) are left out. A package tool becomes an internal
    /// [`ToolDescription`] with its category and no types; internal tools
    /// come from the `.ttd` registry, with their types.
    ///
    /// # Errors
    ///
    /// As [`package_tools`](Self::package_tools), as the `.ttd` reading, and
    /// [`Error::InvalidValue`] with the source's message `the value '<name>'
    /// was used but is not valid; Duplicate tool name error: Trying to add
    /// internal tool '<name>` — the missing closing quote is the source's —
    /// when an internal tool's name is already taken. Every `.ttd` entry is an
    /// internal tool here, as in the source, and an external entry has no name,
    /// so two external entries collide on the empty name.
    pub fn get_topp_tool_list(&self) -> Result<ToolListType> {
        let mut tools = ToolListType::new();
        for (name, entry) in self.package_tools()? {
            if !DESKTOP_CATEGORIES.contains(&entry.category.as_str()) {
                let description = ToolDescription::new(&name, &entry.category, &[]);
                tools.insert(name, description);
            }
        }
        for tool in self.internal_tools()? {
            if tools.contains_key(&tool.internal.name) {
                return Err(invalid_value(
                    &format!(
                        "Duplicate tool name error: Trying to add internal tool '{}",
                        tool.internal.name
                    ),
                    &tool.internal.name,
                ));
            }
            tools.insert(tool.internal.name.clone(), tool.clone());
        }
        Ok(tools)
    }

    /// The product version a package manifest records for `toolname`, or an
    /// empty string for an unregistered tool (source `getToolVersion`).
    ///
    /// # Errors
    ///
    /// As [`package_tools`](Self::package_tools).
    pub fn get_tool_version(&self, toolname: &str) -> Result<String> {
        Ok(self
            .package_tools()?
            .remove(toolname)
            .map(|tool| tool.version)
            .unwrap_or_default())
    }

    /// Resolve a tool's executable (source `findExecutable`,
    /// `ToolHandler.cpp:187-208`): a packaged tool's declared executable, then
    /// a file of that name beside the running executable, then a `PATH`
    /// search. Every candidate must be a regular file with execute
    /// permission.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] with the source's `FileNotFound` text: for a packaged
    /// tool whose declared executable is missing or not executable (naming
    /// that path), and for a tool found nowhere (naming `toolname`). Registry
    /// failures as [`package_tools`](Self::package_tools).
    pub fn find_executable(&self, toolname: &str) -> Result<PathBuf> {
        if let Some(tool) = self.package_tools()?.remove(toolname) {
            if tool.executable.is_file() && file::executable(&tool.executable) {
                return Ok(tool.executable);
            }
            return Err(not_found(&tool.executable.display().to_string()));
        }
        let mut candidate = self.sources.executable_directory.join(toolname);
        if cfg!(windows) && !toolname.ends_with(".exe") {
            candidate.set_extension("exe");
        }
        if candidate.is_file() && file::executable(&candidate) {
            return Ok(candidate);
        }
        if let Some(found) = self.sources.file_context.find_executable(toolname)? {
            if found.is_file() && file::executable(&found) {
                return std::path::absolute(found).map_err(Error::Io);
            }
        }
        Err(not_found(toolname))
    }

    /// The `-type` values of a tool, or an empty list when it has none or is
    /// unknown (source `getTypes`). Only internal `.ttd` tools carry types.
    ///
    /// # Errors
    ///
    /// As [`get_topp_tool_list`](Self::get_topp_tool_list).
    pub fn get_types(&self, toolname: &str) -> Result<Vec<String>> {
        Ok(self
            .get_topp_tool_list()?
            .remove(toolname)
            .map(|tool| tool.internal.types)
            .unwrap_or_default())
    }

    /// The category of a tool, or an empty string when it is unknown (source
    /// `getCategory`).
    ///
    /// # Errors
    ///
    /// As [`get_topp_tool_list`](Self::get_topp_tool_list).
    pub fn get_category(&self, toolname: &str) -> Result<String> {
        Ok(self
            .get_topp_tool_list()?
            .remove(toolname)
            .map(|tool| tool.internal.category)
            .unwrap_or_default())
    }

    /// The internal-tool configuration directory, the root of the `.ttd`
    /// search (source `getInternalToolsPath`: `File::getOpenMSDataPath() +
    /// "/TOOLS/INTERNAL"`).
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when no shared-data directory resolved, where the source's
    /// `File::getOpenMSDataPath` throws `FileNotFound`.
    pub fn get_internal_tools_path(&self) -> Result<PathBuf> {
        self.sources.internal_tools_path.clone().ok_or_else(|| {
            Error::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "cannot find OpenMS shared data; configure runtime data candidates or OPENMS_DATA_PATH",
            ))
        })
    }

    /// The `.ttd` files of the internal-tool registry (source
    /// `getInternalToolConfigFiles_`): those directly inside the internal
    /// tools directory, its `LINUX` or `WINDOWS` subdirectory and
    /// `OPENMS_TTD_INTERNAL_PATH`, each directory listed by name. Missing or
    /// unreadable directories contribute nothing, as the source's
    /// `File::fileList` returns `false` for them.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a listing exceeds a bound or holds a file
    /// name that is not valid UTF-8.
    pub fn internal_tool_config_files(&self) -> Result<Vec<PathBuf>> {
        let mut directories = Vec::new();
        if let Some(root) = &self.sources.internal_tools_path {
            directories.push(root.clone());
            directories.push(root.join(if cfg!(windows) { "WINDOWS" } else { "LINUX" }));
        }
        if let Some(extra) = &self.sources.ttd_internal_path {
            directories.push(extra.clone());
        }
        let mut files = Vec::new();
        for directory in directories {
            match file::file_list(&directory, "*.ttd", true) {
                Ok(found) => files.extend(found),
                Err(Error::Io(_)) => {}
                Err(error) => return Err(error),
            }
        }
        Ok(files)
    }

    /// The internal-tool registry, read once and kept after the first
    /// successful read (source `getInternalTools_`); a failed read is tried
    /// again on the next call.
    fn internal_tools(&self) -> Result<&[ToolDescription]> {
        if let Some(tools) = self.internal_tools.get() {
            return Ok(tools);
        }
        let mut tools = Vec::new();
        let mut diagnostics = Vec::new();
        for path in self.internal_tool_config_files()? {
            let loaded = ToolDescriptionFile::load(&path)?;
            diagnostics.extend(loaded.diagnostics);
            tools.extend(loaded.descriptions);
        }
        if let Ok(mut lines) = self.diagnostics.lock() {
            lines.extend(diagnostics);
        }
        let _ = self.internal_tools.set(tools);
        Ok(self.internal_tools.get().map_or(&[], Vec::as_slice))
    }
}

/// Read one manifest file, bounded, with the source's `FileNotReadable` for a
/// file that cannot be opened or read.
fn read_manifest(path: &Path) -> Result<String> {
    let file = std::fs::File::open(path).map_err(|_| not_readable(path))?;
    let length = file.metadata().map_err(|_| not_readable(path))?.len();
    if length > MAX_MANIFEST_BYTES {
        return Err(Error::InvalidValue(format!(
            "tool package manifest '{}' exceeds {MAX_MANIFEST_BYTES} bytes",
            path.display()
        )));
    }
    let mut text = String::new();
    let mut reader = BufReader::new(file).take(MAX_MANIFEST_BYTES + 1);
    let mut buffer = Vec::new();
    loop {
        buffer.clear();
        let count = reader
            .read_until(b'\n', &mut buffer)
            .map_err(|_| not_readable(path))?;
        if count == 0 {
            break;
        }
        text.push_str(&String::from_utf8_lossy(&buffer));
    }
    Ok(text)
}

/// Whether `path` is the manifest [`BUILTIN_MANIFEST`] stands for: a file
/// named [`BUILTIN_MANIFEST_NAME`] holding exactly its bytes. Any failure to
/// read it answers `false`; reading it for the registry then reports the
/// failure as the source does.
fn is_builtin_manifest(path: &Path) -> bool {
    if path
        .file_name()
        .is_none_or(|name| name != BUILTIN_MANIFEST_NAME)
    {
        return false;
    }
    let expected = BUILTIN_MANIFEST.as_bytes();
    let Ok(file) = std::fs::File::open(path) else {
        return false;
    };
    if file.metadata().map(|m| m.len()).ok() != Some(expected.len() as u64) {
        return false;
    }
    let mut bytes = Vec::with_capacity(expected.len());
    file.take(expected.len() as u64 + 1)
        .read_to_end(&mut bytes)
        .is_ok_and(|_| bytes == expected)
}

use std::io::Read as _;
