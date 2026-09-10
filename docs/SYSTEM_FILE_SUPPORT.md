# Filesystem and runtime paths

`openms::system::file` covers the 39 public function names (40 overloads) and temporary-directory operations of the pinned Core SDK `SYSTEM/File`. Ordinary filesystem work uses safe Rust's standard library. Runtime configuration belongs to an explicit `FileContext`; temporary resources belong to guards. No C++/Qt ABI, subprocess lookup tools, process-exit destructor registry, or implicit `.reference` checkout is required.

## API mapping

| Source surface | Native entry points |
|---|---|
| File queries | `exists`, `empty`, `executable`, `file_size`, `get_modification_time`, `readable`, `writable`, `is_directory` |
| File mutation | `copy`, `rename`, `remove`, `make_dir`, `copy_dir_recursively`, `remove_dir`, `remove_dir_recursively` |
| Lexical paths | `basename`, `path`, `stem_name`, `extension`, `absolute_path` |
| Directory discovery | `file_list`, `list_directories` |
| Executable/runtime discovery | `get_executable_path`, `get_path_locations`, `get_path_locations_from_environment`; context `find_executable`, `find_sibling_topp_executable`, `find`, `find_doc` |
| Data/home/configuration | Context `resolve_data_path`, `get_openms_data_path`, `get_openms_data_path_source`, `get_openms_home_path`, `get_openms_config_dir`, `get_system_parameters`, `get_temp_directory`, `get_user_directory`, `find_database` |
| Temporary resources | `get_unique_name`, `TempDir::new/new_in`, context `get_temporary_file`, `TempFile::new/new_in/alternative` |
| Paired inputs | `validate_matching_file_names` and `MatchingFileListsStatus::{Match,OrderMismatch,SetMismatch}` |

Filesystem results use `PathBuf`/`Path`, retaining native path encodings rather than requiring UTF-8. Paths do not acquire decorative trailing separators; join child names with `Path::join`. Lexical string helpers recognize both slash styles on every platform and preserve the source's trailing-slash rules: `basename("/folder/")` is empty; `path("file")` is `"."`. Stem and extension operations reuse the existing [file-type rules](FILE_HANDLING_SUPPORT.md), including compound compression extensions and historical alias behavior.

Simple predicates return `bool`; other operations return checked `Result` values. Directory-list failures are errors rather than indistinguishable empty results. `empty` deliberately follows the source: missing paths, metadata failures and directories count as empty. File size requires a regular file. Modification times are signed whole seconds since the Unix epoch, with checked overflow rather than unsigned or negative error sentinels.

## Copying, moving, and removal

`copy` refuses an existing destination. It streams into a sibling temporary file, preserves source permissions, synchronizes the completed file, then publishes it with an atomic hard link. `rename(..., false)` publishes regular files with the same no-clobber guarantee; concurrent attempts cannot replace the winner. Same-filesystem moves use a direct hard link followed by source removal. Cross-device moves stage a streamed copy on the destination filesystem. Failed source removal is reported even though the complete destination may already exist.

`rename(..., true)` uses native rename, including ordinary directory and symlink moves. Cross-device regular files use staged copy and publication. Missing/unreadable source files do not cause deletion of the old destination. This corrects the source's delete-before-rename behavior and its POSIX no-clobber flag discrepancy. Cross-device directory/symlink moves and no-clobber directory/symlink moves are checked errors: portable safe std has no atomic no-clobber directory-rename primitive. Filesystems that cannot publish hard links return an error; the implementation does not substitute a racy existence check.

`copy_dir_recursively` supports `CopyOptions::{Overwrite,Skip,Cancel}`. It plans a bounded traversal, rejects copying onto the source or its descendant, detects source directory-link cycles, and checks known conflicts before modifying destinations. Source directory symlinks are followed for copying; destination directory symlinks are rejected to avoid redirecting a planned write. Skip and Cancel remain protected against concurrent file publication. All known Cancel conflicts are detected before copying anything, improving on source filesystem-order partial updates.

Recursive copy and removal are not filesystem transactions: later I/O errors or external mutations may leave completed earlier work. Each file publication is complete, but a directory tree does not become visible in one atomic step. Recursive removal plans the traversal first and does not follow symlinks. Both source removal names are recursive. `remove` removes an ordinary file, a symlink, or an empty directory; missing paths succeed, and dangling symlinks are actually removed. External hostile replacement of directories during a traversal is outside the guarantees of path-based standard-library operations.

## Explicit runtime context

```rust
use openms::system::file::{FileContext, ResolvedDataPath, TempDir};

fn example() -> openms::Result<()> {
    let mut files = FileContext::from_environment()?;
    files.data_candidates.insert(0, ResolvedDataPath {
        path: "/opt/openms/share/OpenMS".into(),
        source: "application installation".into(),
    });
    let resource = files.find("CHEMISTRY/unimod.xml", &[])?;
    assert!(resource.exists());
    let scratch = TempDir::new_in(files.get_temp_directory()?, false)?;
    std::fs::write(scratch.path().join("intermediate.txt"), b"result")?;
    scratch.close()?;
    Ok(())
}
```

`FileContext::new(executable_directory, home_directory, temporary_directory)` constructs an isolated context. `from_environment` snapshots relevant environment values and the current executable. Public fields allow explicit library/install/build data candidates, documentation directories, configuration location, overrides and executable search paths. Snapshotting avoids hidden environment mutation and makes parallel tests independent.

An explicit `OPENMS_DATA_PATH` is authoritative, including an empty/invalid value. A valid data root contains `CHEMISTRY/unimod.xml`, as checked by source existence rules. Otherwise ordered candidates are tested; defaults are executable-relative `../share/OpenMS`, plus the macOS bundle location. Automatic loaded-library address discovery and compiled C++ developer paths are replaced by explicit candidates. Packaged Rust chemistry resources do not masquerade as a complete source-layout runtime data tree.

Only successful data resolution is cached. Failure can be retried after correcting the context. A success remains cached after later field changes until `clear_data_cache`; `resolve_data_path` reports both path and origin. `find` returns an already existing input unchanged. Otherwise data resolution precedes even caller-supplied search directories, preserving source operation order; supplied directories then take priority over the data directory. `find_doc` combines explicit documentation directories with source data-relative documentation locations.

PATH components use the platform separator, replace backslashes, and gain a slash. An empty whole PATH gives no entries; empty components become `/`, preserving the source quirk. Executable search checks existence and excludes directories, without testing execute permissions. Direct existing paths win; Windows searches configured PATHEXT suffixes when the requested name has no dot, with `.exe`/`.bat` fallback for an invalid environment PATHEXT. Sibling lookup searches only the executable directory and source macOS bundle alternatives, without falling back to PATH.

## Configuration and temporary lifetime

Configuration uses [Param](PARAM_SUPPORT.md) and [ParamXML](PARAMXML_SUPPORT.md). `OPENMS_HOME_PATH` is the actual source environment key, despite contradictory header wording. Unix configuration is `XDG_CONFIG_HOME/OpenMS` when set, otherwise `<home>/.config/OpenMS`; other platforms use `<home>/.OpenMS`. Locations are not created merely by querying them.

Missing/unreadable `OpenMS.ini` returns five defaults: SDK version, empty `home_dir`/`temp_dir`, empty `id_db_dir` list, and `threads=1`. An existing readable configuration requires `paramxml`; without that feature the call returns `Unsupported`. Stale/missing version handling preserves a source bug: the returned original tree gets its version updated, while the repaired defaults tree is discarded. Missing defaults remain missing and no file is rewritten. `get_system_parameters_with_warnings` exposes diagnostics as strings. Temp/user directory queries load configuration before considering their overrides, so malformed configuration still errors. Whitespace tests determine whether configured directory values are active; nonempty values retain their original whitespace. An explicitly empty user-home override resolves to `/`, matching source trailing-slash normalization.

`TempDir` creates an exclusive directory, including missing parents, and removes it when dropped unless constructed with `keep_dir=true` or consumed by `keep`. `TempFile` reserves an actual file exclusively, unlike the source's uncreated filename registry. `alternative` and a nonempty context alternative path remain unowned and are neither created nor removed. `path` borrows the path, `keep` transfers cleanup responsibility, and `close` reports cleanup errors. Drop is bounded best-effort cleanup; process termination without stack unwinding does not run guards. Retain guards for as long as subprocesses or other consumers use their paths.

Unique names use UTC Gregorian date/time, process ID and an atomic counter. Optional host text comes from `HOSTNAME`/`COMPUTERNAME` and is sanitized; it can be absent. This replaces source OS hostname and local-time calls without adding a dependency. Names are not random security credentials: exclusive creation supplies collision protection, and native temporary permissions are restricted on Unix.

## Platform conventions, bounds, and evidence

Permission queries are effective-process semantic probes, not exact real-UID `access(2)`/Windows ACL calculations. Existing regular files are opened without truncation; directories and absent targets use exclusive sibling probes that never touch the requested basename. FIFOs/devices are not opened by queries. Long directories that cannot accommodate the probe's basename can report false even if a shorter target would fit. Executability follows Unix execute bits or source Windows existence behavior.

Wildcard matching is deterministic, case-sensitive Unicode: `*`, `?`, bracket ranges/classes with `!`/`^` negation, and backslash escapes. Dot files are ordinary names. Locale-specific POSIX classes/collation and Windows shell pattern-list/case rules are not emulated; filenames requiring wildcard matching must be UTF-8. Native filesystem operations otherwise preserve OS paths.

Directory/search/matching work has 100,000 charged entry visits, 128 traversal levels, 50 million work units, 64 MiB conservative logical allocation, and 1 MiB per path/pattern. Entry accounting includes planning and enumeration visits, not just output rows. Temporary creation retries at most 100 times. File contents are streamed without a small scientific-data-size cap. Param/ParamXML retain their own checked parsing/update budgets. Ordinary caller-owned context/guard values and Rust `Clone` are not secretly fallible.

[system_file.rs](../tests/system_file.rs) exercises source literals and native concurrency, lifecycle, link and error boundaries. Two private tests check basename avoidance for probes and cumulative wildcard work. [system_file_provenance.json](../tests/data/system_file_provenance.json) records source hashes, exact packaged fixtures, source line references and documented native adaptations. No C++ build or execution was used.
