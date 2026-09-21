// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The TOPPBase members the lifecycle suite (`tests/topp_cli_lifecycle.rs`)
//! did not reach: the tool registry (`ToolHandler`), the tool descriptions
//! (`-write_ctd` and the refused CWL and JSON writers), the per-user defaults,
//! `-log` and `-debug`, `-instance`, `checkParam_`, citations and the
//! `Common UTIL options:` heading.
//!
//! Evidence:
//!
//! * **Release build, executed (tier 1).** `../oracle/toppbase-completion/cases.sh`
//!   ran the C++ Release build at the port's pins (core `bc9cc12`, cli
//!   `c19e494`, topp `174b576`) on ibminode06 in a clean environment, one
//!   directory per case. The case directories are retained under
//!   `tests/data/topp_cli_completion/<case>/` (argv, environment, exit code,
//!   stdout, stderr, the log file, the file tree) with their hashes in
//!   `fixtures.sha256.json`; the CTDs the Release build wrote are under
//!   `tests/data/topp_cli_lifecycle/release/`. Paths of the Release run are
//!   mapped to this run's before comparing, the terminal-width probe line
//!   `stty: 'standard input': Inappropriate ioctl for device` is dropped (a
//!   run in process writes to explicit streams, which are not probed; the
//!   executables are, `tests/topp_cli_console.rs`), and log timestamps are
//!   masked.
//! * **Upstream class tests (tier 3).** `ToolHandler_test.cpp`,
//!   `ToolManifest_test.cpp` and the `-log` and `Citation::toString` sections
//!   of `TOPPBase_test.cpp` (cli `c19e494`), `ToolDescriptionFile_test.cpp`
//!   (core `bc9cc12`), with their literals.
//!
//! The registry is passed explicitly ([`run_with_registry`]) so no case
//! depends on, or changes, the process environment; the cases that exercise
//! the environment itself run the built executables.
#![cfg(all(feature = "mzml", feature = "paramxml"))]

#[path = "support/took_line.rs"]
mod took_line;

use openms::cli::tools::{
    BaselineFilter, DTAExtractor, MapNormalizer, MzMLSplitter, SpectraFilterWindowMower,
};
use openms::cli::{
    BUILTIN_MANIFEST, BUILTIN_MANIFEST_NAME, CITE_OPENMS, Citation, ExitCode, ParamCtdFile,
    TOPP_PRODUCT_VERSION, Tool, ToolContext, ToolDescriptionFile, ToolHandler, ToolRegistrySources,
    ToolResult, ToolSpec, product_versions, run_with_registry,
};
use openms::param::{Param, ParamValue};
use openms::system::file::TempDir;
use openms::{Error, Result};
use std::cell::RefCell;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// The Release install prefix and scratch root the oracle ran in.
const RELEASE_PREFIX: &str = "/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576";
const ORACLE_ROOT: &str = "/scratch/kohlbach/toppbase-oracle/run1";
/// The terminal-width probe line the Release build prints before usage text.
const STTY_LINE: &str = "stty: 'standard input': Inappropriate ioctl for device\n";

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/topp_cli_completion")
        .join(name)
}

fn text(path: impl AsRef<Path>) -> String {
    path.as_ref().to_string_lossy().into_owned()
}

fn strings(values: &[&str]) -> Vec<String> {
    values.iter().map(|v| (*v).to_owned()).collect()
}

/// One retained oracle file.
fn oracle(case: &str, file: &str) -> String {
    fs::read_to_string(data(case).join(file)).unwrap_or_else(|e| panic!("{case}/{file}: {e}"))
}

fn oracle_exit(case: &str) -> i32 {
    oracle(case, "exit_code.txt").trim().parse().unwrap()
}

/// A fresh case directory with a `home` and a `cwd`, as the oracle's.
struct Case {
    dir: TempDir,
}

impl Case {
    fn new() -> Self {
        let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
        fs::create_dir_all(dir.path().join("home")).unwrap();
        fs::create_dir_all(dir.path().join("cwd")).unwrap();
        Self { dir }
    }
    fn home(&self) -> PathBuf {
        self.dir.path().join("home")
    }
    fn cwd(&self) -> PathBuf {
        self.dir.path().join("cwd")
    }
    /// The oracle's text for `case` with its paths mapped to this case's.
    fn map(&self, case: &str, content: &str) -> String {
        content
            .replace(
                &format!("{ORACLE_ROOT}/results/{case}/cwd"),
                &text(self.cwd()),
            )
            .replace(
                &format!("{ORACLE_ROOT}/results/{case}/home"),
                &text(self.home()),
            )
            .replace(
                &format!("{ORACLE_ROOT}/inputs/baseline_filter_tool_input.mzML"),
                &baseline_input(),
            )
            .replace(&format!("{ORACLE_ROOT}/inputs/"), &text(data("")))
            .replace(&format!("{RELEASE_PREFIX}/bin/"), "")
            .replace(STTY_LINE, "")
    }
}

/// Explicit registry sources: the process's executable, `PATH` and data
/// directory, with no `OPENMS_TOOL_PREFIX_PATH`, no internal-tool directory
/// and the case's `home` as the user directory, then `configure`.
fn registry(case: &Case, configure: impl FnOnce(&mut ToolRegistrySources)) -> ToolHandler {
    let mut sources = ToolRegistrySources::from_environment().unwrap();
    sources.prefixes.clear();
    sources.internal_tools_path = None;
    sources.ttd_internal_path = None;
    sources.file_context.user_override = Some(case.home());
    configure(&mut sources);
    ToolHandler::new(sources)
}

/// Where the built-in manifest stands in for the executable's own manifest.
fn builtin_manifest_path(handler: &ToolHandler) -> String {
    let prefix = handler
        .sources()
        .executable_directory
        .parent()
        .unwrap()
        .to_path_buf();
    text(prefix.join("share/openms4/tools/topp.tools.tsv"))
}

struct Outcome {
    code: ExitCode,
    out: String,
    err: String,
}

fn run_in<T: Tool>(handler: &ToolHandler, args: &[&str]) -> Outcome {
    let arguments: Vec<String> = std::iter::once(T::NAME.to_owned())
        .chain(args.iter().map(|a| (*a).to_owned()))
        .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with_registry::<T>(&arguments, &mut out, &mut err, handler);
    Outcome {
        code,
        out: String::from_utf8_lossy(&out).into_owned(),
        err: String::from_utf8_lossy(&err).into_owned(),
    }
}

fn assert_oracle_exit(case: &str, outcome: &Outcome) {
    assert_eq!(
        outcome.code.as_i32(),
        oracle_exit(case),
        "{case}: {}",
        outcome.err
    );
}

/// A log file with its timestamps masked.
fn masked_log(content: &str) -> Vec<String> {
    content
        .split_inclusive('\n')
        .map(|line| {
            let bytes = line.as_bytes();
            let stamped = bytes.len() > 20
                && bytes[4] == b'-'
                && bytes[7] == b'-'
                && bytes[10] == b' '
                && bytes[13] == b':'
                && bytes[16] == b':'
                && bytes[19] == b' '
                && bytes[..4].iter().all(u8::is_ascii_digit);
            if stamped {
                format!("<time> {}", &line[20..])
            } else {
                line.to_owned()
            }
        })
        .collect()
}

fn baseline_input() -> String {
    text(Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/baseline_filter_tool_input.mzML"))
}

// ---------------------------------------------------------------------------
// ToolHandler_test.cpp and ToolManifest_test.cpp (tier 3)
// ---------------------------------------------------------------------------

/// The upstream fixture registry, `tests/registry` of the cli package, which
/// its class tests reach through `OPENMS_TOOL_PREFIX_PATH`.
fn fixture_registry() -> ToolHandler {
    let prefix = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/tool_handler/registry");
    ToolHandler::new(ToolRegistrySources::only_prefixes(vec![prefix]).unwrap())
}

/// `START_SECTION((static ToolListType getTOPPToolList()))`
/// (`ToolHandler_test.cpp:39-49`).
#[test]
fn upstream_get_topp_tool_list() {
    let handler = fixture_registry();
    let list = handler.get_topp_tool_list().unwrap();
    assert!(list.contains_key("DecoyDatabase"));
    assert!(list.contains_key("IndependentTool"));
    assert_eq!(
        handler.get_tool_version("IndependentTool").unwrap(),
        "7.2.1"
    );
    assert_eq!(handler.get_tool_version("DOESNOTEXIST").unwrap(), "");
    // TEST_EXCEPTION(Exception::FileNotFound, findExecutable("IndependentTool")):
    // the declared bin/IndependentTool does not exist.
    let error = handler.find_executable("IndependentTool").unwrap_err();
    assert!(
        matches!(&error, Error::Io(e) if e.kind() == std::io::ErrorKind::NotFound),
        "{error}"
    );
}

/// `START_SECTION((static StringList getTypes(...)))` and
/// `getCategory` (`ToolHandler_test.cpp:51-69`).
#[test]
fn upstream_get_types_and_get_category() {
    let handler = fixture_registry();
    assert!(handler.get_types("IsobaricAnalyzer").unwrap().is_empty());
    assert!(handler.get_types("IDMapper").unwrap().is_empty());
    assert_eq!(
        handler.get_category("IDFilter").unwrap(),
        "File Filtering, Extraction and Merging"
    );
    assert_eq!(handler.get_category("DOESNOTEXIST").unwrap(), "");
}

/// `START_SECTION((static std::string getInternalToolsPath()))`: not empty.
/// It is the shared-data directory plus `TOOLS/INTERNAL`; without a shared-data
/// directory this port reports the source's `FileNotFound` as an error rather
/// than stop every tool (see `ToolRegistrySources::from_environment`).
#[test]
fn upstream_get_internal_tools_path() {
    let mut sources = ToolRegistrySources::only_prefixes(Vec::new()).unwrap();
    sources.internal_tools_path = Some(PathBuf::from("/data/OpenMS/TOOLS/INTERNAL"));
    let handler = ToolHandler::new(sources);
    let path = handler.get_internal_tools_path().unwrap();
    assert!(!path.as_os_str().is_empty());
    assert!(path.ends_with("TOOLS/INTERNAL"));
    let bare = ToolHandler::new(ToolRegistrySources::only_prefixes(Vec::new()).unwrap());
    assert!(bare.get_internal_tools_path().is_err());
}

/// One temporary prefix with `share/openms4/tools` and `bin`, as
/// `ToolManifest_test.cpp`'s `RegistryFixture`.
struct Prefix(TempDir);
impl Prefix {
    fn new() -> Self {
        let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
        fs::create_dir_all(dir.path().join("share/openms4/tools")).unwrap();
        fs::create_dir_all(dir.path().join("bin")).unwrap();
        Self(dir)
    }
    fn path(&self) -> PathBuf {
        self.0.path().to_path_buf()
    }
    fn write(&self, contents: &str) {
        fs::write(
            self.0.path().join("share/openms4/tools/probe.tools.tsv"),
            contents,
        )
        .unwrap();
    }
    fn binary(&self, name: &str) -> PathBuf {
        self.0.path().join("bin").join(name)
    }
    fn handler(&self) -> ToolHandler {
        ToolHandler::new(ToolRegistrySources::only_prefixes(vec![self.path()]).unwrap())
    }
}

fn invalid(error: &Error) -> bool {
    matches!(error, Error::InvalidValue(_))
}

/// `START_SECTION(manifest parser rejects duplicate and unsafe entries)`
/// (`ToolManifest_test.cpp:107-126`).
#[test]
fn upstream_manifest_parser_rejects_duplicate_and_unsafe_entries() {
    let fixture = Prefix::new();
    for manifest in [
        "__PackageProbe\tExperimental\t1.0\tbin/Probe\n__PackageProbe\tExperimental\t2.0\tbin/Other\n",
        "__PackageProbe\tExperimental\t1.0\t../Probe\n",
        "__PackageProbe\tExperimental\t1.0\t/bin/Probe\n",
        "__PackageProbe\tExperimental\t1.0\tbin/Probe\t\n",
    ] {
        fixture.write(manifest);
        let error = fixture
            .handler()
            .get_tool_version("__PackageProbe")
            .unwrap_err();
        assert!(invalid(&error), "{manifest:?}: {error}");
    }
}

/// Whether permissions bind this process (they do not for root).
fn permissions_apply(path: &Path) -> bool {
    !openms::system::file::readable(path)
}

/// `START_SECTION(unreadable registry inputs fail with library exceptions)`
/// (`ToolManifest_test.cpp:128-147`).
#[cfg(unix)]
#[test]
fn upstream_unreadable_registry_inputs_fail() {
    use std::os::unix::fs::PermissionsExt;
    let fixture = Prefix::new();
    fixture.write("__PackageProbe\tExperimental\t1.0\tbin/Probe\n");
    let directory = fixture.path().join("share/openms4/tools");
    let manifest = directory.join("probe.tools.tsv");
    fs::set_permissions(&manifest, fs::Permissions::from_mode(0o000)).unwrap();
    if permissions_apply(&manifest) {
        let error = fixture
            .handler()
            .get_tool_version("__PackageProbe")
            .unwrap_err();
        // Exception::FileNotReadable
        assert!(
            matches!(&error, Error::Io(e) if e.kind() == std::io::ErrorKind::PermissionDenied),
            "{error}"
        );
    }
    fs::set_permissions(&manifest, fs::Permissions::from_mode(0o644)).unwrap();
    fs::set_permissions(&directory, fs::Permissions::from_mode(0o000)).unwrap();
    if permissions_apply(&directory) {
        let error = fixture
            .handler()
            .get_tool_version("__PackageProbe")
            .unwrap_err();
        // Exception::InvalidValue
        assert!(invalid(&error), "{error}");
    }
    fs::set_permissions(&directory, fs::Permissions::from_mode(0o755)).unwrap();
}

/// `START_SECTION(declared executables must exist and be executable files)`
/// (`ToolManifest_test.cpp:149-166`).
#[test]
fn upstream_declared_executables_must_exist_and_be_executable() {
    let fixture = Prefix::new();
    fixture.write("__PackageProbe\tExperimental\t1.0\tbin/Probe\n");
    let not_found =
        |error: &Error| matches!(error, Error::Io(e) if e.kind() == std::io::ErrorKind::NotFound);
    let error = fixture
        .handler()
        .find_executable("__PackageProbe")
        .unwrap_err();
    assert!(not_found(&error), "{error}");
    fs::create_dir(fixture.binary("Probe")).unwrap();
    let error = fixture
        .handler()
        .find_executable("__PackageProbe")
        .unwrap_err();
    assert!(not_found(&error), "{error}");
    fs::remove_dir(fixture.binary("Probe")).unwrap();
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        fs::write(
            fixture.binary("Probe"),
            "This file has no execute permission.\n",
        )
        .unwrap();
        fs::set_permissions(fixture.binary("Probe"), fs::Permissions::from_mode(0o600)).unwrap();
        let error = fixture
            .handler()
            .find_executable("__PackageProbe")
            .unwrap_err();
        assert!(not_found(&error), "{error}");
        let error = fixture
            .handler()
            .find_executable(&text(fixture.binary("Probe")))
            .unwrap_err();
        assert!(not_found(&error), "{error}");
    }
}

/// `START_SECTION(interactive desktop tools remain resolvable without CLI
/// parameter discovery)` (`ToolManifest_test.cpp:168-195`): the probe copies
/// its own executable into the prefix.
#[test]
fn upstream_interactive_desktop_tools_remain_resolvable() {
    let fixture = Prefix::new();
    let binary_name = if cfg!(windows) { "Probe.exe" } else { "Probe" };
    fs::copy(
        std::env::current_exe().unwrap(),
        fixture.binary(binary_name),
    )
    .unwrap();
    fixture.write(&format!(
        "# A comment followed by CRLF-terminated package records\r\n__PackageViewer\tDesktopViewer\t1.2.3\tbin/{binary_name}\r\n__PackageWorkflow\tDesktopWorkflow\t1.2.3\tbin/{binary_name}\r\n__PackageCLI\tWorkflows\t1.2.3\tbin/{binary_name}\r\n"
    ));
    let handler = fixture.handler();
    let tools = handler.get_topp_tool_list().unwrap();
    assert!(!tools.contains_key("__PackageViewer"));
    assert!(!tools.contains_key("__PackageWorkflow"));
    assert!(tools.contains_key("__PackageCLI"));
    assert_eq!(
        handler.get_tool_version("__PackageViewer").unwrap(),
        "1.2.3"
    );
    assert_eq!(
        handler.find_executable("__PackageViewer").unwrap(),
        fixture.binary(binary_name)
    );
}

/// The registry probe tool of `ToolManifest_test.cpp:92-100`, which is not
/// an official tool and registers nothing.
struct RegistryProbe;
impl Tool for RegistryProbe {
    const NAME: &'static str = "__PackageProbe";
    const DESCRIPTION: &'static str = "Installed registry test";
    fn register(_spec: &mut ToolSpec) -> Result<()> {
        Ok(())
    }
    fn run(_ctx: &ToolContext) -> ToolResult {
        Ok(ExitCode::ExecutionOk)
    }
}

/// `START_SECTION(duplicates across prefixes and malformed startup fail
/// predictably)` (`ToolManifest_test.cpp:197-212`).
#[test]
fn upstream_duplicates_across_prefixes_fail_predictably() {
    let (first, second) = (Prefix::new(), Prefix::new());
    first.write("__PackageProbe\tExperimental\t1.0\tbin/Probe\n");
    second.write("__PackageProbe\tExperimental\t2.0\tbin/Probe\n");
    let both = ToolHandler::new(
        ToolRegistrySources::only_prefixes(vec![first.path(), second.path()]).unwrap(),
    );
    assert!(invalid(
        &both.get_tool_version("__PackageProbe").unwrap_err()
    ));
    // Construction must not read corrupt registries; main() reports them.
    let outcome = run_in::<RegistryProbe>(&both, &["-help"]);
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
}

/// `START_SECTION(product version is written to INI with a core-version
/// fallback)` (`ToolManifest_test.cpp:214-237`).
#[test]
fn upstream_product_version_is_written_to_ini_with_a_core_version_fallback() {
    let fixture = Prefix::new();
    let case = Case::new();
    let output = text(case.cwd().join("probe.ini"));
    fixture.write("__PackageProbe\tExperimental\t7.8.9\tbin/Probe\n");
    let outcome = run_in::<RegistryProbe>(&fixture.handler(), &["-write_ini", &output]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let parameters = openms::format::paramxml::load(&output).unwrap();
    assert_eq!(
        parameters.value("__PackageProbe:version").unwrap(),
        &ParamValue::String("7.8.9".into())
    );
    fixture.write("# no registered product\n");
    let outcome = run_in::<RegistryProbe>(&fixture.handler(), &["-write_ini", &output]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let parameters = openms::format::paramxml::load(&output).unwrap();
    assert_eq!(
        parameters.value("__PackageProbe:version").unwrap(),
        &ParamValue::String(openms::CORE_SDK_VERSION.into())
    );
}

/// `ToolDescriptionFile_test.cpp` (core `bc9cc12`): both retained `.ttd`
/// files load and are not empty. The source's `store` section expects
/// `NotImplemented`; the port has no `store`.
#[test]
fn upstream_tool_description_file_load() {
    for name in [
        "ToolDescriptionFile_test_1.ttd",
        "ToolDescriptionFile_test_2.ttd",
    ] {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/tool_handler")
            .join(name);
        let loaded = ToolDescriptionFile::load(&path).unwrap();
        assert!(!loaded.descriptions.is_empty(), "{name}");
        assert!(
            loaded.diagnostics.is_empty(),
            "{name}: {:?}",
            loaded.diagnostics
        );
    }
}

/// The second retained `.ttd` file in detail, from its text: an external
/// tool, its category, types, texts, command line, working directory, the
/// three mappings, the two post-moves and the three parameters.
#[test]
fn a_ttd_file_is_read_field_by_field() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/tool_handler/ToolDescriptionFile_test_2.ttd");
    let loaded = ToolDescriptionFile::load(&path).unwrap();
    assert_eq!(loaded.descriptions.len(), 1);
    let tool = &loaded.descriptions[0];
    assert!(!tool.internal.is_internal);
    assert_eq!(tool.internal.name, "");
    assert_eq!(tool.internal.category, "Identification");
    assert_eq!(tool.internal.types, strings(&["Percolator"]));
    assert_eq!(tool.external_details.len(), 1);
    let details = &tool.external_details[0];
    assert_eq!(details.text_startup, "Running Percolator...");
    assert_eq!(
        details.text_fail,
        "Something went wrong. Is the Percolator executable globally accessible?"
    );
    assert_eq!(details.text_finish, "Percolator finished successfully.");
    assert_eq!(details.category, "Identification");
    assert_eq!(
        details.commandline,
        "--only-psms --results-psms \"%2\" --decoy-results-psms \"%3\" \"%1\""
    );
    assert_eq!(details.path, "percolator");
    assert_eq!(details.working_directory, ".");
    let mapping: Vec<(i32, &str)> = details
        .tr_table
        .mapping
        .iter()
        .map(|(id, cl)| (*id, cl.as_str()))
        .collect();
    assert_eq!(
        mapping,
        vec![
            (1, "%%in"),
            (2, "%TMP/%BASENAME[%%in].psms"),
            (3, "%TMP/%BASENAME[%%in]_decoy.psms")
        ]
    );
    assert_eq!(details.tr_table.post_moves.len(), 2);
    assert_eq!(details.tr_table.post_moves[1].target, "out_decoy");
    assert_eq!(
        details.param.value("out_decoy").unwrap(),
        &ParamValue::String(String::new())
    );
    assert!(details.param.has_tag("in", "input file").unwrap());
}

// ---------------------------------------------------------------------------
// The registry in the lifecycle (Release build, tier 1)
// ---------------------------------------------------------------------------

/// A prefix for `OPENMS_TOOL_PREFIX_PATH` under the case's `home`, holding
/// `share/openms4/tools/<file>`.
fn probe_prefix(case: &Case, file: &str, contents: &str) -> PathBuf {
    let prefix = case.home().join("prefix");
    fs::create_dir_all(prefix.join("share/openms4/tools")).unwrap();
    fs::create_dir_all(prefix.join("bin")).unwrap();
    fs::write(prefix.join("share/openms4/tools").join(file), contents).unwrap();
    prefix
}

/// Oracle `reg_*`: a manifest the registry refuses ends every run, even
/// `--help`, with exit 6 and `Unable to initialize or run <tool>: <what>`; the
/// duplicate is reported at its second occurrence, which for a tool of the
/// product manifest is the built-in manifest standing in for the Release
/// build's installed one.
#[test]
fn a_manifest_the_registry_refuses_ends_the_run_as_in_the_release_build() {
    let cases: &[(&str, &str, &[&str])] = &[
        (
            "reg_dup_help",
            "BaselineFilter\tProbe\t9.9.9\tbin/BaselineFilter\n",
            &["--help"],
        ),
        (
            "reg_dup_write_ini",
            "BaselineFilter\tProbe\t9.9.9\tbin/BaselineFilter\n",
            &["-test", "-write_ini", "@CWD@/x.ini"],
        ),
        ("reg_malformed", "just one field\n", &["--help"]),
        (
            "reg_extra_tab",
            "__Probe\tExperimental\t1.0\tbin/Probe\t\n",
            &["--help"],
        ),
        (
            "reg_unsafe_path",
            "__Probe\tExperimental\t1.0\t../Probe\n",
            &["--help"],
        ),
        (
            "reg_absolute_path",
            "__Probe\tExperimental\t1.0\t/bin/Probe\n",
            &["--help"],
        ),
        (
            "reg_slash_name",
            "a/b\tExperimental\t1.0\tbin/Probe\n",
            &["--help"],
        ),
    ];
    for (name, manifest, args) in cases {
        let case = Case::new();
        let prefix = probe_prefix(&case, "probe.tools.tsv", manifest);
        let handler = registry(&case, |sources| sources.prefixes = vec![prefix]);
        let cwd = text(case.cwd());
        let args: Vec<String> = args.iter().map(|a| a.replace("@CWD@", &cwd)).collect();
        let refs: Vec<&str> = args.iter().map(String::as_str).collect();
        let outcome = run_in::<BaselineFilter>(&handler, &refs);
        assert_oracle_exit(name, &outcome);
        let expected = case.map(name, &oracle(name, "stderr.txt")).replace(
            &format!("{RELEASE_PREFIX}/share/openms4/tools/topp.tools.tsv"),
            &builtin_manifest_path(&handler),
        );
        assert_eq!(outcome.err, expected, "{name}");
        assert!(outcome.out.is_empty(), "{name}: {}", outcome.out);
    }
}

/// Oracles `reg_unreadable_manifest` and `reg_unreadable_dir`: the source's
/// `FileNotReadable` and its `filesystem_error` text inside `InvalidValue`.
#[cfg(unix)]
#[test]
fn an_unreadable_manifest_or_directory_ends_the_run_as_in_the_release_build() {
    use std::os::unix::fs::PermissionsExt;
    for (name, target) in [
        (
            "reg_unreadable_manifest",
            "share/openms4/tools/probe.tools.tsv",
        ),
        ("reg_unreadable_dir", "share/openms4/tools"),
    ] {
        let case = Case::new();
        let prefix = probe_prefix(
            &case,
            "probe.tools.tsv",
            "__Probe\tExperimental\t1.0\tbin/Probe\n",
        );
        let locked = prefix.join(target);
        fs::set_permissions(&locked, fs::Permissions::from_mode(0o000)).unwrap();
        if permissions_apply(&locked) {
            let handler = registry(&case, |sources| sources.prefixes = vec![prefix.clone()]);
            let outcome = run_in::<BaselineFilter>(&handler, &["--help"]);
            assert_oracle_exit(name, &outcome);
            assert_eq!(
                outcome.err,
                case.map(name, &oracle(name, "stderr.txt")),
                "{name}"
            );
        }
        let mode = if target.ends_with(".tsv") {
            0o644
        } else {
            0o755
        };
        fs::set_permissions(&locked, fs::Permissions::from_mode(mode)).unwrap();
    }
}

/// Oracles `reg_extra_tool`, `reg_empty_category`, `reg_tsv_directory`,
/// `reg_other_extension`, `reg_missing_prefix` and `reg_empty_items`: an
/// additional tool, a comment and CRLF line ends, an empty category, a
/// directory named `*.tsv`, a manifest with another extension, a missing
/// prefix and empty prefix items leave the tool's usage unchanged.
#[test]
fn registry_entries_that_do_not_concern_the_tool_leave_its_usage_unchanged() {
    let cases: &[(&str, Option<(&str, &str)>)] = &[
        (
            "reg_extra_tool",
            Some((
                "probe.tools.tsv",
                "# name\tcategory\tversion\texecutable\r\n__Probe\tExperimental\t7.2.1\tbin/__Probe\r\n__Viewer\tDesktopViewer\t1.2.3\tbin/__Viewer\n",
            )),
        ),
        (
            "reg_empty_category",
            Some(("probe.tools.tsv", "__Probe\t\t1.0\tbin/Probe\n")),
        ),
        (
            "reg_other_extension",
            Some((
                "probe.tools.txt",
                "BaselineFilter\tProbe\t9.9.9\tbin/BaselineFilter\n",
            )),
        ),
        ("reg_tsv_directory", None),
    ];
    for (name, manifest) in cases {
        let case = Case::new();
        let prefix = match manifest {
            Some((file, contents)) => probe_prefix(&case, file, contents),
            None => {
                let prefix = case.home().join("prefix");
                fs::create_dir_all(prefix.join("share/openms4/tools/sub.tsv")).unwrap();
                prefix
            }
        };
        let handler = registry(&case, |sources| sources.prefixes = vec![prefix]);
        let outcome = run_in::<BaselineFilter>(&handler, &["--help"]);
        assert_oracle_exit(name, &outcome);
        assert_eq!(
            outcome.err,
            case.map(name, &oracle(name, "stderr.txt")),
            "{name}"
        );
    }
    let case = Case::new();
    let missing = case.home().join("nonexistent");
    let handler = registry(&case, |sources| sources.prefixes = vec![missing]);
    let outcome = run_in::<BaselineFilter>(&handler, &["--help"]);
    assert_oracle_exit("reg_missing_prefix", &outcome);
    assert_eq!(
        outcome.err,
        case.map(
            "reg_missing_prefix",
            &oracle("reg_missing_prefix", "stderr.txt")
        )
    );
    assert_eq!(
        ToolRegistrySources::from_prefix_path(if cfg!(windows) { ";;" } else { "::" }),
        Vec::<PathBuf>::new()
    );
}

/// Oracle `reg_install_prefix_twice`: the Release build with its own
/// installation prefix in `OPENMS_TOOL_PREFIX_PATH`, twice, runs normally,
/// because the source reads a manifest file once however many prefixes
/// reach it (`ToolHandler.cpp:90-104`). A prefix installing the product
/// manifest byte for byte, as `share/openms4/tools/topp.tools.tsv`, is that
/// installation for this port: the built-in manifest standing in for the
/// executable's prefix is the same manifest and is not read a second time,
/// so the usage text is the Release build's, and the registered tools resolve
/// to that prefix.
#[test]
fn a_prefix_that_installs_the_product_manifest_is_the_built_in_one_read_once() {
    let case = Case::new();
    let prefix = probe_prefix(&case, BUILTIN_MANIFEST_NAME, BUILTIN_MANIFEST);
    let handler = registry(&case, |sources| {
        sources.prefixes = vec![prefix.clone(), prefix.clone()];
    });
    let outcome = run_in::<BaselineFilter>(&handler, &["--help"]);
    assert_oracle_exit("reg_install_prefix_twice", &outcome);
    assert_eq!(
        outcome.err,
        case.map(
            "reg_install_prefix_twice",
            &oracle("reg_install_prefix_twice", "stderr.txt")
        )
    );
    assert_eq!(
        handler.get_tool_version("BaselineFilter").unwrap(),
        TOPP_PRODUCT_VERSION
    );
    let tools = handler.package_tools().unwrap();
    assert_eq!(
        tools["BaselineFilter"].executable,
        std::path::absolute(prefix.join("bin/BaselineFilter")).unwrap()
    );
    let once = registry(&case, |sources| sources.prefixes = vec![prefix.clone()]);
    assert_eq!(once.package_tools().unwrap(), tools);
}

/// A second manifest is a duplicate, whatever it holds, as for two C++
/// installations (oracle `reg_dup_help`, and `reg_copy_of_install_manifest`
/// of `../oracle/topp-exception-exits`, where the Release build refuses a
/// byte-identical copy of its own manifest under another prefix: the source
/// knows a manifest by its path). Only the product manifest itself, under its
/// installed name, is the built-in one reached again: a copy under another
/// name, and `topp.tools.tsv` one byte longer, are second manifests listing
/// the product tools.
#[test]
fn a_manifest_that_is_not_the_product_manifest_is_still_a_duplicate() {
    let longer = format!("{BUILTIN_MANIFEST}\n");
    for (file, contents) in [
        ("copy.tools.tsv", BUILTIN_MANIFEST),
        (BUILTIN_MANIFEST_NAME, longer.as_str()),
    ] {
        let case = Case::new();
        let prefix = probe_prefix(&case, file, contents);
        let handler = registry(&case, |sources| sources.prefixes = vec![prefix]);
        let outcome = run_in::<BaselineFilter>(&handler, &["--help"]);
        assert_eq!(
            outcome.code,
            ExitCode::IllegalParameters,
            "{file}: {}",
            outcome.err
        );
        assert!(
            outcome
                .err
                .starts_with("Unable to initialize or run BaselineFilter: ")
                && outcome.err.ends_with(
                    "was used but is not valid; Invalid or duplicate tool package manifest entry\n"
                ),
            "{file}: {}",
            outcome.err
        );
        assert!(
            outcome.err.contains(&builtin_manifest_path(&handler)),
            "{file}: the duplicate is the built-in manifest's row: {}",
            outcome.err
        );
    }
}

/// Oracles `ttd_*`: the internal-tool registry under `OPENMS_TTD_INTERNAL_PATH`.
/// Every `.ttd` entry is an internal tool, keyed by name, so two external
/// entries (which have none) collide on the empty name; a tool the product
/// manifest lists collides with it; an unknown `status` and an element the
/// status does not allow are non-fatal messages; a file in a subdirectory is
/// not read.
#[test]
fn the_internal_tool_registry_behaves_as_in_the_release_build() {
    let cases: &[(&str, &[&str], Option<&str>)] = &[
        (
            "ttd_one_external",
            &["ToolDescriptionFile_test_1.ttd"],
            None,
        ),
        (
            "ttd_two_external",
            &[
                "ToolDescriptionFile_test_1.ttd",
                "ToolDescriptionFile_test_2.ttd",
            ],
            None,
        ),
        ("ttd_internal_new", &["ttd_internal_new.ttd"], None),
        ("ttd_internal_dup", &["ttd_internal_dup.ttd"], None),
        ("ttd_unknown_status", &["ttd_unknown_status.ttd"], None),
        ("ttd_nested", &["ttd_internal_dup.ttd"], Some("deeper")),
    ];
    let source = |name: &str| {
        let upstream = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/tool_handler")
            .join(name);
        if upstream.exists() {
            upstream
        } else {
            data(name)
        }
    };
    for (name, files, subdirectory) in cases {
        let case = Case::new();
        let ttd = case.home().join("ttd");
        let target = match subdirectory {
            Some(sub) => ttd.join(sub),
            None => ttd.clone(),
        };
        fs::create_dir_all(&target).unwrap();
        for file in *files {
            fs::copy(source(file), target.join(file)).unwrap();
        }
        let handler = registry(&case, |sources| {
            sources.ttd_internal_path = Some(ttd.clone())
        });
        let outcome = run_in::<BaselineFilter>(&handler, &["--help"]);
        assert_oracle_exit(name, &outcome);
        assert_eq!(
            outcome.err,
            case.map(name, &oracle(name, "stderr.txt")),
            "{name}"
        );
    }
}

/// Oracle `ttd_malformed`: exit 6. The Release build's message is Xerces's
/// (`input ended before all started tags were ended; ...`) after the file
/// name; this port's reader words the reason itself.
#[test]
fn a_malformed_ttd_file_ends_the_run() {
    let case = Case::new();
    let ttd = case.home().join("ttd");
    fs::create_dir_all(&ttd).unwrap();
    fs::copy(data("ttd_malformed.ttd"), ttd.join("ttd_malformed.ttd")).unwrap();
    let handler = registry(&case, |sources| {
        sources.ttd_internal_path = Some(ttd.clone())
    });
    let outcome = run_in::<BaselineFilter>(&handler, &["--help"]);
    assert_oracle_exit("ttd_malformed", &outcome);
    let prefix = format!(
        "Unable to initialize or run BaselineFilter: While loading '{}': ",
        text(ttd.join("ttd_malformed.ttd"))
    );
    assert!(outcome.err.starts_with(&prefix), "{}", outcome.err);
    assert!(
        oracle("ttd_malformed", "stderr.txt").contains(&format!(
            "Unable to initialize or run BaselineFilter: While loading '{ORACLE_ROOT}/results/ttd_malformed/home/ttd/ttd_malformed.ttd': "
        ))
    );
}

/// The registry read from the process environment, through the built
/// executable: `OPENMS_TOOL_PREFIX_PATH` names a prefix whose manifest lists
/// `BaselineFilter` again (oracle `reg_dup_help`). The executable's own prefix
/// has no manifest, so the built-in manifest stands in for it and reports the
/// duplicate there.
#[test]
fn the_executable_reads_the_prefix_path_from_its_environment() {
    let case = Case::new();
    let prefix = probe_prefix(
        &case,
        "probe.tools.tsv",
        "BaselineFilter\tProbe\t9.9.9\tbin/BaselineFilter\n",
    );
    let binary = PathBuf::from(env!("CARGO_BIN_EXE_BaselineFilter"));
    let output = Command::new(&binary)
        .arg("--help")
        .env("OPENMS_TOOL_PREFIX_PATH", &prefix)
        .env("OPENMS_HOME_PATH", case.home())
        .current_dir(case.cwd())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(oracle_exit("reg_dup_help")));
    let own_prefix = binary.parent().unwrap().parent().unwrap();
    let expected = case
        .map("reg_dup_help", &oracle("reg_dup_help", "stderr.txt"))
        .replace(
            &format!("{RELEASE_PREFIX}/share/openms4/tools/topp.tools.tsv"),
            &text(own_prefix.join("share/openms4/tools/topp.tools.tsv")),
        );
    assert_eq!(String::from_utf8_lossy(&output.stderr), expected);

    // `OPENMS_TOOL_PREFIX_PATH` naming an installation of the product
    // manifest, which is what the variable is for: the executable runs, as
    // the Release build does with its own prefix named there (oracle
    // `reg_install_prefix_twice`).
    let case = Case::new();
    let prefix = probe_prefix(&case, BUILTIN_MANIFEST_NAME, BUILTIN_MANIFEST);
    let output = Command::new(&binary)
        .arg("--help")
        .env("OPENMS_TOOL_PREFIX_PATH", &prefix)
        .env("OPENMS_HOME_PATH", case.home())
        .env("COLUMNS", "0")
        .current_dir(case.cwd())
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(oracle_exit("reg_install_prefix_twice")),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    let release = case.map(
        "reg_install_prefix_twice",
        &oracle("reg_install_prefix_twice", "stderr.txt"),
    );
    assert_eq!(stderr, release);
}

// ---------------------------------------------------------------------------
// Tool descriptions (Release build, tier 1)
// ---------------------------------------------------------------------------

fn release_ctd(tool: &str) -> String {
    fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/topp_cli_lifecycle/release")
            .join(format!("{tool}.ctd")),
    )
    .unwrap()
}

/// `-write_ctd <dir>` under `-test`, as the topp package's
/// `CheckToolMetadata.py` runs it for every tool.
fn write_ctd<T: Tool>() -> String {
    let case = Case::new();
    let handler = registry(&case, |_| {});
    let dir = case.cwd().join("ctd");
    fs::create_dir_all(&dir).unwrap();
    let outcome = run_in::<T>(&handler, &["-test", "-write_ctd", &text(&dir)]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert!(outcome.err.is_empty(), "{}", outcome.err);
    assert!(outcome.out.is_empty(), "{}", outcome.out);
    fs::read_to_string(dir.join(format!("{}.ctd", T::NAME))).unwrap()
}

/// Oracles `ctd_<tool>`: the CTD of every ported tool that declares its
/// citations is byte-identical to the Release build's.
/// `FeatureFinderCentroided` is not among them: its source cites two
/// publications (`FeatureFinderCentroided.cpp:124-136`), which its port does
/// not declare yet, so its CTD lacks their two `<citation>` lines.
#[test]
fn write_ctd_matches_the_release_build_for_the_ported_tools() {
    assert_eq!(write_ctd::<BaselineFilter>(), release_ctd("BaselineFilter"));
    assert_eq!(write_ctd::<DTAExtractor>(), release_ctd("DTAExtractor"));
    assert_eq!(write_ctd::<MapNormalizer>(), release_ctd("MapNormalizer"));
    assert_eq!(write_ctd::<MzMLSplitter>(), release_ctd("MzMLSplitter"));
    assert_eq!(
        write_ctd::<SpectraFilterWindowMower>(),
        release_ctd("SpectraFilterWindowMower")
    );
    assert_eq!(
        write_ctd::<openms::cli::tools::PeakPickerHiRes>(),
        release_ctd("PeakPickerHiRes")
    );
    #[cfg(feature = "featurexml")]
    assert_eq!(
        write_ctd::<openms::cli::tools::FileInfo>(),
        release_ctd("FileInfo")
    );
}

/// The topp package's `CheckToolMetadata.py`, which upstream registers as
/// `<tool>_write_ctd` for every tool: the root is `<tool>`, its `name`,
/// `version` and `category` attributes are the tool's registry entry, and the
/// `<tool>:version` item carries the product version.
#[test]
fn upstream_check_tool_metadata_holds_for_the_ported_tools() {
    let case = Case::new();
    let handler = registry(&case, |_| {});
    let check = |name: &str, ctd: &str| {
        let category = handler.get_category(name).unwrap();
        let version = handler.get_tool_version(name).unwrap();
        assert!(!category.is_empty() && version == "1.0.0", "{name}");
        let head = format!(
            "<tool ctdVersion=\"1.8\" version=\"{version}\" name=\"{name}\" docurl=\"http://www.openms.de/doxygen/release/4.0.0/html/TOPP_{name}.html\" category=\"{category}\" >\n"
        );
        assert!(ctd.contains(&head), "{name}: {ctd}");
        let item = format!(
            "<NODE name=\"{name}\" description=\"{}\">\n    <ITEM name=\"version\" value=\"{version}\" type=\"string\"",
            ctd_escape(match name {
                "BaselineFilter" => BaselineFilter::DESCRIPTION,
                "DTAExtractor" => DTAExtractor::DESCRIPTION,
                _ => MapNormalizer::DESCRIPTION,
            })
        );
        assert!(ctd.contains(&item), "{name}: {ctd}");
    };
    check("BaselineFilter", &write_ctd::<BaselineFilter>());
    check("DTAExtractor", &write_ctd::<DTAExtractor>());
    check("MapNormalizer", &write_ctd::<MapNormalizer>());
}

fn ctd_escape(text: &str) -> String {
    ParamCtdFile::escape_xml(text)
}

/// Oracles `ctd_notest_BaselineFilter`, `ctd_with_ini`, `ctd_with_cli`,
/// `ctd_trailing_slash`, `ctd_existing`: the CTD holds the defaults whatever
/// else is given, `-test` changes nothing, and an existing file is replaced.
#[test]
fn write_ctd_writes_the_defaults_whatever_else_is_given() {
    let expected = release_ctd("BaselineFilter");
    for (name, extra) in [
        ("ctd_notest_BaselineFilter", vec![]),
        ("ctd_with_ini", vec!["-ini", "@IN@/method_erosion.ini"]),
        ("ctd_with_cli", vec!["-method", "erosion"]),
    ] {
        let case = Case::new();
        let handler = registry(&case, |_| {});
        let dir = case.cwd().join("ctd");
        fs::create_dir_all(&dir).unwrap();
        let input_dir = text(data(""));
        let mut args = vec!["-write_ctd".to_owned(), text(&dir)];
        if name != "ctd_notest_BaselineFilter" {
            args.insert(0, "-test".to_owned());
        }
        args.extend(extra.iter().map(|a| a.replace("@IN@/", &input_dir)));
        let refs: Vec<&str> = args.iter().map(String::as_str).collect();
        let outcome = run_in::<BaselineFilter>(&handler, &refs);
        assert_oracle_exit(name, &outcome);
        assert_eq!(
            fs::read_to_string(dir.join("BaselineFilter.ctd")).unwrap(),
            expected,
            "{name}"
        );
    }
    let case = Case::new();
    let handler = registry(&case, |_| {});
    let dir = case.cwd().join("ctd");
    fs::create_dir_all(&dir).unwrap();
    fs::write(dir.join("BaselineFilter.ctd"), "old content\nsecond line\n").unwrap();
    let with_slash = format!("{}/", text(&dir));
    let outcome = run_in::<BaselineFilter>(&handler, &["-test", "-write_ctd", &with_slash]);
    assert_oracle_exit("ctd_existing", &outcome);
    assert_oracle_exit("ctd_trailing_slash", &outcome);
    assert_eq!(
        fs::read_to_string(dir.join("BaselineFilter.ctd")).unwrap(),
        expected
    );
}

/// Oracles `ctd_missing_dir`, `ctd_target_is_directory` and
/// `write_json_missing_dir`: the target is checked as an output file (exit 5
/// with the source's two lines); a directory by the target's name passes that
/// check and fails when the CTD is opened, where the source's
/// `std::ios::failure` reaches the initialisation catch (exit 12).
#[test]
fn write_ctd_target_failures_exit_as_in_the_release_build() {
    for (name, writer, extension) in [
        ("ctd_missing_dir", "-write_ctd", "ctd"),
        ("write_json_missing_dir", "-write_json", "json"),
    ] {
        let case = Case::new();
        let handler = registry(&case, |_| {});
        let missing = case.cwd().join("nodir");
        let outcome = run_in::<BaselineFilter>(&handler, &["-test", writer, &text(&missing)]);
        assert_oracle_exit(name, &outcome);
        assert_eq!(
            outcome.err,
            case.map(name, &oracle(name, "stderr.txt")),
            "{name}"
        );
        assert!(!missing.join(format!("BaselineFilter.{extension}")).exists());
    }
    let case = Case::new();
    let handler = registry(&case, |_| {});
    let dir = case.cwd().join("ctd");
    fs::create_dir_all(dir.join("BaselineFilter.ctd")).unwrap();
    let outcome = run_in::<BaselineFilter>(&handler, &["-test", "-write_ctd", &text(&dir)]);
    assert_oracle_exit("ctd_target_is_directory", &outcome);
    assert_eq!(
        outcome.err,
        case.map(
            "ctd_target_is_directory",
            &oracle("ctd_target_is_directory", "stderr.txt")
        )
    );
}

/// Oracles `ctd_and_write_ini`, `ctd_and_cwl` and
/// `write_nested_cwl_and_json`: `-write_ini` is handled first, then
/// `-write_ctd`, then `-write_nested_cwl`, `-write_cwl`,
/// `-write_nested_json` and `-write_json`; only the first given is written.
#[test]
fn only_the_first_description_writer_given_runs() {
    let case = Case::new();
    let handler = registry(&case, |_| {});
    let (ctd, ini) = (case.cwd().join("ctd"), case.cwd().join("x.ini"));
    fs::create_dir_all(&ctd).unwrap();
    let outcome = run_in::<BaselineFilter>(
        &handler,
        &[
            "-test",
            "-write_ctd",
            &text(&ctd),
            "-write_ini",
            &text(&ini),
        ],
    );
    assert_oracle_exit("ctd_and_write_ini", &outcome);
    assert!(ini.exists() && !ctd.join("BaselineFilter.ctd").exists());

    let (d1, d2) = (case.cwd().join("d1"), case.cwd().join("d2"));
    fs::create_dir_all(&d1).unwrap();
    fs::create_dir_all(&d2).unwrap();
    let outcome = run_in::<BaselineFilter>(
        &handler,
        &["-test", "-write_ctd", &text(&d1), "-write_cwl", &text(&d2)],
    );
    assert_oracle_exit("ctd_and_cwl", &outcome);
    assert_eq!(
        fs::read_to_string(d1.join("BaselineFilter.ctd")).unwrap(),
        release_ctd("BaselineFilter")
    );

    let out = case.cwd().join("out");
    fs::create_dir_all(&out).unwrap();
    let outcome = run_in::<BaselineFilter>(
        &handler,
        &[
            "-test",
            "-write_json",
            &text(&out),
            "-write_nested_cwl",
            &text(&out),
        ],
    );
    assert_oracle_exit("write_nested_cwl_and_json", &outcome);
    assert_eq!(
        outcome.err,
        case.map(
            "write_nested_cwl_and_json",
            &oracle("write_nested_cwl_and_json", "stderr.txt")
        )
    );
}

/// Oracles `ctd_cwd_BaselineFilter` and `ctd_relative_BaselineFilter`: an
/// empty directory is the current directory, and a relative one is relative
/// to it. Through the executable, which runs in the case's directory.
#[test]
fn write_ctd_to_the_current_directory() {
    for (name, directory) in [
        ("ctd_cwd_BaselineFilter", ""),
        ("ctd_relative_BaselineFilter", "."),
    ] {
        let case = Case::new();
        let output = Command::new(env!("CARGO_BIN_EXE_BaselineFilter"))
            .args(["-test", "-write_ctd", directory])
            .env("OPENMS_HOME_PATH", case.home())
            .env_remove("OPENMS_TOOL_PREFIX_PATH")
            .current_dir(case.cwd())
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(oracle_exit(name)), "{output:?}");
        assert_eq!(
            fs::read_to_string(case.cwd().join("BaselineFilter.ctd")).unwrap(),
            release_ctd("BaselineFilter"),
            "{name}"
        );
    }
}

// ---------------------------------------------------------------------------
// Usage text (Release build, tier 1)
// ---------------------------------------------------------------------------

/// Oracles `help_<tool>` and `helphelp_<tool>`: the usage text of every
/// ported tool that declares its citations is the Release build's, byte for
/// byte.
#[test]
fn usage_matches_the_release_build_for_the_ported_tools() {
    fn check<T: Tool>() {
        for flag in ["--help", "--helphelp"] {
            let case = Case::new();
            let handler = registry(&case, |_| {});
            let outcome = run_in::<T>(&handler, &[flag]);
            let name = format!("{}_{}", flag.trim_start_matches('-'), T::NAME);
            assert_oracle_exit(&name, &outcome);
            assert_eq!(
                outcome.err,
                case.map(&name, &oracle(&name, "stderr.txt")),
                "{name}"
            );
            assert!(outcome.out.is_empty(), "{name}: {}", outcome.out);
        }
    }
    check::<BaselineFilter>();
    check::<DTAExtractor>();
    check::<MapNormalizer>();
    check::<MzMLSplitter>();
    check::<SpectraFilterWindowMower>();
    check::<openms::cli::tools::PeakPickerHiRes>();
    #[cfg(feature = "featurexml")]
    check::<openms::cli::tools::FileInfo>();
}

/// A stand-in with `FeatureFinderCentroided`'s name, description and the
/// source's two citations (`FeatureFinderCentroided.cpp:124-136`), to check
/// the citation block against the Release build's `--help` (oracle
/// `ffc_help_cite`) up to the usage line.
struct CitingTool;
impl Tool for CitingTool {
    const NAME: &'static str = "FeatureFinderCentroided";
    const DESCRIPTION: &'static str = "Detects two-dimensional features in LC-MS data.";
    const CITATIONS: &'static [Citation] = &[
        Citation {
            authors: "Sturm M",
            title: "A novel feature detection algorithm for centroided data",
            when_where: "Dissertation, 2010-09-15, p.37 ff",
            doi: "https://publikationen.uni-tuebingen.de/xmlui/bitstream/handle/10900/49453/pdf/Dissertation_Marc_Sturm.pdf",
        },
        Citation {
            authors: "Weisser H",
            title: "An automated pipeline for high-throughput label-free quantitative proteomics",
            when_where: "J. Proteome Res., 2013, PMID: 23391308",
            doi: "https://doi.org/10.1021/pr300992u",
        },
    ];
    fn register(_spec: &mut ToolSpec) -> Result<()> {
        Ok(())
    }
    fn run(_ctx: &ToolContext) -> ToolResult {
        Ok(ExitCode::ExecutionOk)
    }
}

#[test]
fn tool_citations_are_printed_after_the_openms_citation() {
    let case = Case::new();
    let handler = registry(&case, |_| {});
    let outcome = run_in::<CitingTool>(&handler, &["--help"]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk);
    let expected = case.map("ffc_help_cite", &oracle("ffc_help_cite", "stderr.txt"));
    let head = |text: &str| text[..text.find("Usage:").unwrap()].to_owned();
    assert_eq!(head(&outcome.err), head(&expected));

    // The CTD lists the OpenMS DOI and then each tool citation's.
    let dir = case.cwd().join("ctd");
    fs::create_dir_all(&dir).unwrap();
    let outcome = run_in::<CitingTool>(&handler, &["-write_ctd", &text(&dir)]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let written = fs::read_to_string(dir.join("FeatureFinderCentroided.ctd")).unwrap();
    let release = release_ctd("FeatureFinderCentroided");
    let citations = |text: &str| {
        text[text.find("<citations>").unwrap()..text.find("</citations>").unwrap()].to_owned()
    };
    assert_eq!(citations(&written), citations(&release));
}

/// `START_SECTION(([EXTRA] Citation::toString()))` (`TOPPBase_test.cpp:892-897`).
#[test]
fn upstream_citation_to_string() {
    let c = Citation {
        authors: "Surname I",
        title: "A title",
        when_where: "Journal. 2024; 1:2-3",
        doi: "10.1000/xyz",
    };
    assert_eq!(
        c.to_source_string(),
        "Surname I. A title. Journal. 2024; 1:2-3. doi:10.1000/xyz."
    );
    assert_eq!(c.to_string(), c.to_source_string());
    assert_eq!(CITE_OPENMS.doi, "10.1038/s41592-024-02197-7");
}

/// A tool that no manifest lists reports the core version and lists its
/// common options as UTIL options (`TOPPBase.cpp:119-126`, `160-163`).
#[test]
fn an_unregistered_tool_reports_the_core_version_and_util_options() {
    let case = Case::new();
    let handler = registry(&case, |_| {});
    let outcome = run_in::<RegistryProbe>(&handler, &["--help"]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk);
    assert!(
        outcome.err.contains("Common UTIL options:"),
        "{}",
        outcome.err
    );
    assert!(!outcome.err.contains("Common TOPP options:"));
    let (version, verbose) = product_versions(&handler, RegistryProbe::NAME).unwrap();
    assert_eq!(version, openms::CORE_SDK_VERSION);
    assert!(outcome.err.contains(&format!("Version: {verbose}\n")));
    let (version, _) = product_versions(&handler, BaselineFilter::NAME).unwrap();
    assert_eq!(version, "1.0.0");
}

// ---------------------------------------------------------------------------
// Per-user defaults (Release build, tier 1)
// ---------------------------------------------------------------------------

/// A case whose `home` holds `BaselineFilter.ini` from `file`.
fn with_user_defaults(file: &str) -> Case {
    let case = Case::new();
    fs::copy(data(file), case.home().join("BaselineFilter.ini")).unwrap();
    case
}

/// Oracles `ud_write_ini`, `ud_unknown`, `ud_invalid`, `ud_invalid_write_ini`,
/// `ud_with_ini`: `<user directory>/BaselineFilter.ini` updates the defaults
/// leniently, with the source's verbose messages on standard error, and an
/// INI or the command line still overrides them.
#[test]
fn user_defaults_update_the_defaults_as_in_the_release_build() {
    for (name, file, run) in [
        ("ud_write_ini", "user_defaults.ini", false),
        ("ud_unknown", "user_defaults_unknown.ini", true),
        ("ud_invalid", "user_defaults_invalid.ini", true),
        ("ud_invalid_write_ini", "user_defaults_invalid.ini", false),
        ("ud_run", "user_defaults.ini", true),
    ] {
        let case = with_user_defaults(file);
        let handler = registry(&case, |_| {});
        let written = case.cwd().join("out.ini");
        let output = case.cwd().join("out.mzML");
        let input = baseline_input();
        let args: Vec<&str> = if run {
            vec![
                "-test",
                "-no_progress",
                "-in",
                &input,
                "-out",
                output.to_str().unwrap(),
            ]
        } else {
            vec!["-test", "-write_ini", written.to_str().unwrap()]
        };
        let outcome = run_in::<BaselineFilter>(&handler, &args);
        assert_oracle_exit(name, &outcome);
        assert_eq!(
            outcome.err,
            case.map(name, &oracle(name, "stderr.txt")),
            "{name}"
        );
        if !run {
            let port = openms::format::paramxml::load(&written).unwrap();
            let release =
                openms::format::paramxml::read(oracle(name, "out.ini").as_bytes()).unwrap();
            assert!(port.source_equal(&release).unwrap(), "{name}");
        }
    }
    let case = with_user_defaults("user_defaults.ini");
    let handler = registry(&case, |_| {});
    let output = case.cwd().join("out.mzML");
    let ini = text(data("method_erosion.ini"));
    let input = baseline_input();
    let outcome = run_in::<BaselineFilter>(
        &handler,
        &[
            "-test",
            "-no_progress",
            "-in",
            &input,
            "-out",
            output.to_str().unwrap(),
            "-ini",
            &ini,
            "-struc_elem_length",
            "2",
        ],
    );
    assert_oracle_exit("ud_with_ini", &outcome);
    assert_eq!(
        outcome.err,
        case.map("ud_with_ini", &oracle("ud_with_ini", "stderr.txt"))
    );
}

/// Oracle `ud_ctd`: the CTD holds the per-user defaults.
#[test]
fn user_defaults_reach_the_ctd() {
    let case = with_user_defaults("user_defaults.ini");
    let handler = registry(&case, |_| {});
    let dir = case.cwd().join("ctd");
    fs::create_dir_all(&dir).unwrap();
    let outcome = run_in::<BaselineFilter>(&handler, &["-test", "-write_ctd", &text(&dir)]);
    assert_oracle_exit("ud_ctd", &outcome);
    assert_eq!(
        outcome.err,
        case.map("ud_ctd", &oracle("ud_ctd", "stderr.txt"))
    );
    assert_eq!(
        fs::read_to_string(dir.join("BaselineFilter.ctd")).unwrap(),
        release_ctd("BaselineFilter_user_defaults")
    );
}

/// Oracles `ud_help`, `ud_malformed_help` and `ud_unreadable`: usage does not
/// read the per-user defaults, and an unreadable file is skipped.
#[test]
fn user_defaults_are_not_read_for_usage_and_skipped_when_unreadable() {
    for (name, file) in [
        ("ud_help", "user_defaults.ini"),
        ("ud_malformed_help", "user_defaults_malformed.ini"),
    ] {
        let case = with_user_defaults(file);
        let handler = registry(&case, |_| {});
        let outcome = run_in::<BaselineFilter>(&handler, &["--help"]);
        assert_oracle_exit(name, &outcome);
        assert_eq!(
            outcome.err,
            case.map(name, &oracle(name, "stderr.txt")),
            "{name}"
        );
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let case = with_user_defaults("user_defaults.ini");
        let ini = case.home().join("BaselineFilter.ini");
        fs::set_permissions(&ini, fs::Permissions::from_mode(0o000)).unwrap();
        if permissions_apply(&ini) {
            let handler = registry(&case, |_| {});
            let written = case.cwd().join("out.ini");
            let outcome =
                run_in::<BaselineFilter>(&handler, &["-test", "-write_ini", &text(&written)]);
            assert_oracle_exit("ud_unreadable", &outcome);
            assert_eq!(outcome.err, "");
        }
        fs::set_permissions(&ini, fs::Permissions::from_mode(0o644)).unwrap();
    }
}

/// Oracles `ud_malformed` and `ud_directory`: a file that does not parse is
/// `INPUT_FILE_CORRUPT` (the message after the file name is the XML reader's
/// own), and a directory by that name is `INTERNAL_ERROR` with libstdc++'s
/// message.
#[test]
fn unusable_user_defaults_end_the_run_as_in_the_release_build() {
    let case = with_user_defaults("user_defaults_malformed.ini");
    let handler = registry(&case, |_| {});
    let output = case.cwd().join("out.mzML");
    let input = baseline_input();
    let outcome = run_in::<BaselineFilter>(
        &handler,
        &["-test", "-in", &input, "-out", output.to_str().unwrap()],
    );
    assert_oracle_exit("ud_malformed", &outcome);
    let file = format!("{}//BaselineFilter.ini", text(case.home()));
    assert!(
        outcome.err.starts_with(&format!(
            "Error: Unable to read file (While loading '{file}': "
        )),
        "{}",
        outcome.err
    );
    assert!(
        case.map("ud_malformed", &oracle("ud_malformed", "stderr.txt"))
            .contains(&format!(
                "Error: Unable to read file (While loading '{file}': "
            ))
    );

    let case = Case::new();
    fs::create_dir_all(case.home().join("BaselineFilter.ini")).unwrap();
    let handler = registry(&case, |_| {});
    let written = case.cwd().join("out.ini");
    let outcome = run_in::<BaselineFilter>(&handler, &["-test", "-write_ini", &text(&written)]);
    assert_oracle_exit("ud_directory", &outcome);
    assert_eq!(
        outcome.err,
        case.map("ud_directory", &oracle("ud_directory", "stderr.txt"))
    );
}

// ---------------------------------------------------------------------------
// -log and -debug (Release build, tier 1)
// ---------------------------------------------------------------------------

/// A BaselineFilter run in a fresh case with `args`, where `@CWD@` and
/// `@IN@` stand for the case's working directory and the retained input.
fn logged_run(args: &[&str]) -> (Case, Outcome) {
    let case = Case::new();
    let handler = registry(&case, |_| {});
    let (cwd, input) = (text(case.cwd()), baseline_input());
    let args: Vec<String> = args
        .iter()
        .map(|a| a.replace("@CWD@", &cwd).replace("@IN@", &input))
        .collect();
    let refs: Vec<&str> = args.iter().map(String::as_str).collect();
    let outcome = run_in::<BaselineFilter>(&handler, &refs);
    (case, outcome)
}

fn log_lines(case: &Case) -> Vec<String> {
    masked_log(&fs::read_to_string(case.cwd().join("log.txt")).unwrap_or_default())
}

fn oracle_log(case: &Case, name: &str) -> Vec<String> {
    masked_log(&case.map(name, &oracle(name, "log.txt")))
}

/// Oracles `log_run` and `log_unwritable`: without a line to write, no log
/// file appears; a destination that cannot be opened is not reported.
#[test]
fn a_log_file_appears_only_when_something_is_logged() {
    let (case, outcome) = logged_run(&[
        "-test",
        "-in",
        "@IN@",
        "-out",
        "@CWD@/out.mzML",
        "-log",
        "@CWD@/log.txt",
    ]);
    assert_oracle_exit("log_run", &outcome);
    assert!(!case.cwd().join("log.txt").exists());
    let (case, outcome) = logged_run(&[
        "-test",
        "-in",
        "@CWD@/missing.mzML",
        "-out",
        "@CWD@/out.mzML",
        "-log",
        "@CWD@/nodir/log.txt",
    ]);
    assert_oracle_exit("log_unwritable", &outcome);
    assert_eq!(
        outcome.err,
        case.map("log_unwritable", &oracle("log_unwritable", "stderr.txt"))
    );
}

/// Oracles `log_unknown_option`, `log_missing_input`, `log_help_debug1`,
/// `log_write_ini_debug1` and `log_unknown_option_debug1`: the log file holds
/// the source's lines, errors always and debug lines from their level, and
/// the `Writing to` notice goes to standard output at debug level 1.
#[test]
fn the_log_file_holds_the_release_builds_lines() {
    let cases: &[(&str, &[&str])] = &[
        ("log_unknown_option", &["-log", "@CWD@/log.txt", "-bogus"]),
        (
            "log_missing_input",
            &[
                "-test",
                "-in",
                "@CWD@/missing.mzML",
                "-out",
                "@CWD@/out.mzML",
                "-log",
                "@CWD@/log.txt",
            ],
        ),
        (
            "log_help_debug1",
            &["--help", "-log", "@CWD@/log.txt", "-debug", "1"],
        ),
        (
            "log_write_ini_debug1",
            &[
                "-test",
                "-write_ini",
                "@CWD@/x.ini",
                "-log",
                "@CWD@/log.txt",
                "-debug",
                "1",
            ],
        ),
        (
            "log_unknown_option_debug1",
            &["-log", "@CWD@/log.txt", "-debug", "1", "-bogus"],
        ),
    ];
    for (name, args) in cases {
        let (case, outcome) = logged_run(args);
        assert_oracle_exit(name, &outcome);
        assert_eq!(
            outcome.err,
            case.map(name, &oracle(name, "stderr.txt")),
            "{name}"
        );
        assert_eq!(
            outcome.out,
            case.map(name, &oracle(name, "stdout.txt")),
            "{name}"
        );
        assert_eq!(log_lines(&case), oracle_log(&case, name), "{name}");
    }
}

/// Oracles `log_run_debug1` and `log_run_debug2`: a run at debug levels 1 and
/// 2. The framework's lines up to the tool body are the Release build's in
/// order. The body's `Value of … option` lines follow the order in which the
/// tool reads its options, which in this port's BaselineFilter differs from
/// the source's; at level 2 the `Checking input file` and `Checking output
/// file` lines come from the validation that runs before the body, where the
/// source checks each file when the tool first reads it. Those lines are
/// compared as a collection.
#[test]
fn a_debug_run_logs_the_release_builds_lines() {
    for (name, level) in [("log_run_debug1", "1"), ("log_run_debug2", "2")] {
        let (case, outcome) = logged_run(&[
            "-test",
            "-in",
            "@IN@",
            "-out",
            "@CWD@/out.mzML",
            "-log",
            "@CWD@/log.txt",
            "-debug",
            level,
        ]);
        assert_oracle_exit(name, &outcome);
        let expected = oracle_log(&case, name);
        let actual = log_lines(&case);
        // The framework's part ends with its read of -no_progress.
        let framework = 1 + expected
            .iter()
            .position(|line| line.contains("Value of string option 'no_progress'"))
            .unwrap();
        assert_eq!(actual[..framework], expected[..framework], "{name}");
        let mut rest_actual = actual[framework..].to_vec();
        let mut rest_expected = expected[framework..].to_vec();
        rest_actual.sort();
        rest_expected.sort();
        assert_eq!(rest_actual, rest_expected, "{name}");
        let first = outcome.out.lines().next().unwrap_or_default();
        assert_eq!(
            format!("{first}\n"),
            case.map(name, oracle(name, "stdout.txt").lines().next().unwrap()) + "\n"
        );
        let (_, took) = took_line::split_took_line("BaselineFilter", &outcome.out);
        assert!(took.is_some(), "{name}: {}", outcome.out);
    }
}

/// Oracle `log_in_ini`: a `log` item in the INI file takes effect once the INI
/// is merged, and the version notice reaches it.
#[test]
fn a_log_file_named_by_the_ini_file_takes_effect() {
    let case = Case::new();
    let handler = registry(&case, |_| {});
    let log = case.cwd().join("log.txt");
    let ini = case.home().join("log.ini");
    fs::write(
        &ini,
        fs::read_to_string(data("log_in_ini.ini"))
            .unwrap()
            .replace("@LOG@", &text(&log)),
    )
    .unwrap();
    let output = case.cwd().join("out.mzML");
    let input = baseline_input();
    let outcome = run_in::<BaselineFilter>(
        &handler,
        &[
            "-test",
            "-no_progress",
            "-in",
            &input,
            "-out",
            output.to_str().unwrap(),
            "-ini",
            &text(&ini),
        ],
    );
    assert_oracle_exit("log_in_ini", &outcome);
    assert_eq!(log_lines(&case), oracle_log(&case, "log_in_ini"));
}

/// `START_SECTION(([EXTRA] -log writes a log file))` (`TOPPBase_test.cpp:870-890`):
/// `TOPPBaseTest -log <file> -debug 1` writes a non-empty log whose lines carry
/// the tool's INI location.
#[test]
fn upstream_log_writes_a_log_file() {
    struct ToppBaseTest;
    impl Tool for ToppBaseTest {
        const NAME: &'static str = "TOPPBaseTest";
        const DESCRIPTION: &'static str = "A test class";
        fn register(spec: &mut ToolSpec) -> Result<()> {
            spec.register_string_option(
                "stringoption",
                "<string>",
                "string default",
                "string description",
                false,
                false,
            )
        }
        fn run(_ctx: &ToolContext) -> ToolResult {
            Ok(ExitCode::ExecutionOk)
        }
    }
    let case = Case::new();
    let handler = registry(&case, |_| {});
    let log = case.cwd().join("log.txt");
    let outcome = run_in::<ToppBaseTest>(&handler, &["-log", &text(&log), "-debug", "1"]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let content = fs::read_to_string(&log).unwrap();
    assert!(!content.is_empty());
    assert!(content.contains("TOPPBaseTest:1:"), "{content}");
}

/// `START_SECTION(([EXTRA]std::string const& getIniLocation_() const))`,
/// command-line part (`TOPPBase_test.cpp:436-439`): after `-instance 5` the
/// INI location is `TOPPBaseTest:5:`. The source's class test reads it from
/// the object after `main` has returned; here it is read where a user sees it,
/// in the log lines, which carry the location, and in the INI warning (oracle
/// `instance5_ini`).
#[test]
fn upstream_ini_location_follows_instance() {
    struct ToppBaseTest;
    impl Tool for ToppBaseTest {
        const NAME: &'static str = "TOPPBaseTest";
        const DESCRIPTION: &'static str = "A test class";
        fn register(_spec: &mut ToolSpec) -> Result<()> {
            Ok(())
        }
        fn run(_ctx: &ToolContext) -> ToolResult {
            Ok(ExitCode::ExecutionOk)
        }
    }
    let case = Case::new();
    let handler = registry(&case, |_| {});
    let log = case.cwd().join("log.txt");
    let outcome = run_in::<ToppBaseTest>(
        &handler,
        &["-instance", "5", "-log", &text(&log), "-debug", "1"],
    );
    // -instance is rejected by the strict update, as in the source.
    assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
    let content = fs::read_to_string(&log).unwrap();
    assert!(
        content.contains(" TOPPBaseTest:5:: Debug level: 1\n"),
        "{content}"
    );

    let case = Case::new();
    let handler = registry(&case, |_| {});
    let output = case.cwd().join("out.mzML");
    let input = baseline_input();
    let ini = text(data("method_erosion.ini"));
    let outcome = run_in::<BaselineFilter>(
        &handler,
        &[
            "-test",
            "-in",
            &input,
            "-out",
            output.to_str().unwrap(),
            "-ini",
            &ini,
            "-instance",
            "5",
        ],
    );
    assert_oracle_exit("instance5_ini", &outcome);
    assert_eq!(
        outcome.err,
        case.map("instance5_ini", &oracle("instance5_ini", "stderr.txt"))
    );
    let (case, outcome) = logged_run(&[
        "-test",
        "-in",
        "@IN@",
        "-out",
        "@CWD@/out.mzML",
        "-instance",
        "2",
    ]);
    assert_oracle_exit("instance2", &outcome);
    assert_eq!(
        outcome.err,
        case.map("instance2", &oracle("instance2", "stderr.txt"))
    );
}

/// `-write_ini` with `-instance` writes the parameters under that instance's
/// section, as `getDefaultParameters_` builds them from `getToolPrefix()`
/// (`TOPPBase.cpp:2100`, `2235`); the strict update that rejects `-instance`
/// comes after the write commands.
#[test]
fn write_ini_follows_instance() {
    let case = Case::new();
    let handler = registry(&case, |_| {});
    let written = case.cwd().join("x.ini");
    let outcome =
        run_in::<BaselineFilter>(&handler, &["-write_ini", &text(&written), "-instance", "3"]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let param = openms::format::paramxml::load(&written).unwrap();
    assert!(param.exists("BaselineFilter:3:method").unwrap());
    assert!(!param.exists("BaselineFilter:1:method").unwrap());
    assert_eq!(
        param.section_description("BaselineFilter:3").unwrap(),
        "Instance '3' section for 'BaselineFilter'"
    );
}

// ---------------------------------------------------------------------------
// checkParam_ (oracle cases of ../oracle/topp-cli-lifecycle, tier 1)
// ---------------------------------------------------------------------------

/// Oracles `ini_common_tool_section` and `ini_instance_and_common` of
/// `../oracle/topp-cli-lifecycle` (the product SDK): `checkParam_` warns
/// about the `common:` copy of a `common:<tool>:` subsection value, because
/// its exemption for the tool's own name never applies
/// (`TOPPBase.cpp:1879-1888`).
#[test]
fn check_param_warns_about_the_common_copy_of_a_tool_section() {
    for ini in ["common_tool_section.ini", "instance_and_common.ini"] {
        let case = Case::new();
        let handler = registry(&case, |_| {});
        let path = text(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/data/topp_cli_lifecycle")
                .join(ini),
        );
        let input = text(
            Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/window_mower_tool_input.mzML"),
        );
        let output = case.cwd().join("out.mzML");
        let outcome = run_in::<SpectraFilterWindowMower>(
            &handler,
            &[
                "-test",
                "-ini",
                &path,
                "-in",
                &input,
                "-out",
                output.to_str().unwrap(),
            ],
        );
        assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
        let warning = format!(
            "Warning: Unknown subsection 'SpectraFilterWindowMower:algorithm' in '{path}' (location 'common:')!\n"
        );
        assert_eq!(outcome.err.matches(&warning).count(), 1, "{}", outcome.err);
    }
}

// ---------------------------------------------------------------------------
// Input files tagged is_executable (TOPPBase.cpp:1534-1549, tier 4)
// ---------------------------------------------------------------------------

thread_local! {
    static RESOLVED: RefCell<Option<String>> = const { RefCell::new(None) };
}

struct ExecutableTool;
impl Tool for ExecutableTool {
    const NAME: &'static str = "ExecutableTool";
    const DESCRIPTION: &'static str = "Takes an executable";
    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_input_file(
            "exe",
            "<file>",
            "",
            "an executable",
            false,
            false,
            &["is_executable"],
        )
    }
    fn run(ctx: &ToolContext) -> ToolResult {
        let resolved = ctx.string("exe")?.to_owned();
        RESOLVED.with(|slot| *slot.borrow_mut() = Some(resolved));
        Ok(ExitCode::ExecutionOk)
    }
}

/// An executable found on `PATH` replaces the given name with its full path;
/// one found nowhere is the source's warning and `ExternalExecutableNotFound`,
/// exit 14.
#[cfg(unix)]
#[test]
fn an_is_executable_input_is_resolved_on_path() {
    let case = Case::new();
    let handler = registry(&case, |sources| {
        sources.file_context.search_path = vec![PathBuf::from("/bin"), PathBuf::from("/usr/bin")];
    });
    let outcome = run_in::<ExecutableTool>(&handler, &["-exe", "sh"]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let resolved = RESOLVED.with(|slot| slot.borrow_mut().take()).unwrap();
    assert!(
        resolved.ends_with("/sh") && Path::new(&resolved).is_absolute(),
        "{resolved}"
    );

    let outcome = run_in::<ExecutableTool>(&handler, &["-exe", "no-such-program-for-openms"]);
    assert_eq!(
        outcome.code,
        ExitCode::ExternalProgramNotFound,
        "{}",
        outcome.err
    );
    assert_eq!(
        outcome.err,
        "Input file 'no-such-program-for-openms' could not be found (by searching on PATH). Either provide a full filepath via the '-exe' option or fix your PATH environment ! Since this file is not strictly required, you might also pass the empty string \"\" as argument to prevent its usage (this might limit the usability of the tool).\nError: Executable not found (the executable 'no-such-program-for-openms' could not be found)\n"
    );
}

// ---------------------------------------------------------------------------
// getOutputDirOption and the integer parseRange_ (tier 4)
// ---------------------------------------------------------------------------

thread_local! {
    static OUTPUT_DIR: RefCell<Option<String>> = const { RefCell::new(None) };
}

struct OutputDirTool;
impl Tool for OutputDirTool {
    const NAME: &'static str = "OutputDirTool";
    const DESCRIPTION: &'static str = "Writes into a directory";
    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_output_dir("out_dir", "<directory>", "", "a directory", false, false)
    }
    fn run(ctx: &ToolContext) -> ToolResult {
        let directory = ctx.output_dir("out_dir")?.to_owned();
        OUTPUT_DIR.with(|slot| *slot.borrow_mut() = Some(directory));
        Ok(ExitCode::ExecutionOk)
    }
}

/// `getOutputDirOption` creates the directory and its parents
/// (`TOPPBase.cpp:1419-1439`); an output directory takes no format check.
#[test]
fn an_output_directory_is_created_when_read() {
    let case = Case::new();
    let handler = registry(&case, |_| {});
    let target = case.cwd().join("a/b/c");
    let outcome = run_in::<OutputDirTool>(&handler, &["-out_dir", &text(&target)]);
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    assert_eq!(
        OUTPUT_DIR.with(|slot| slot.borrow_mut().take()).unwrap(),
        text(&target)
    );
    assert!(target.is_dir());
}

/// The integer `parseRange_` (`TOPPBase.cpp:2054-2090`), with the cases of the
/// floating-point section of `TOPPBase_test.cpp:675-709` in integers.
#[test]
fn the_integer_range_parser_follows_the_source() {
    let (mut a, mut b) = (-1, -1);
    assert!(!openms::cli::parse_range_int(":", &mut a, &mut b).unwrap());
    assert_eq!((a, b), (-1, -1));
    assert!(openms::cli::parse_range_int("4:", &mut a, &mut b).unwrap());
    assert_eq!((a, b), (4, -1));
    assert!(openms::cli::parse_range_int(":5", &mut a, &mut b).unwrap());
    assert_eq!((a, b), (4, 5));
    assert!(openms::cli::parse_range_int("6:7", &mut a, &mut b).unwrap());
    assert_eq!((a, b), (6, 7));
    let error = openms::cli::parse_range_int("400", &mut a, &mut b).unwrap_err();
    assert!(
        error.to_string().contains("the ':' separator is missing"),
        "{error}"
    );
    let error = openms::cli::parse_range_int("1.5:2", &mut a, &mut b).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("Could not convert string '1.5:2' to a range of integer values"),
        "{error}"
    );
    assert_eq!((a, b), (6, 7));
}

/// The `Param` the probes need is also what a CTD of a registered subsection
/// writes; `ParamCtdFile` is re-exported for tools that write one directly.
#[test]
fn param_ctd_file_is_exported_with_the_framework() {
    let mut param = Param::new();
    param
        .set_value("x", ParamValue::Integer(1), "d", &[])
        .unwrap();
    let text = ParamCtdFile
        .to_ctd_string(&param, &openms::data_structures::ToolInfo::default())
        .unwrap();
    assert!(text.contains("<ITEM name=\"x\" value=\"1\" type=\"int\" description=\"d\""));
}
