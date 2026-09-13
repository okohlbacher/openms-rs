# Filesystem and runtime paths — `SYSTEM/File.h`

`openms::system::file` is the native counterpart of the pinned Core SDK
`SYSTEM/File` (header `SYSTEM/File.h`, implementation split across
`SYSTEM/File.cpp`, `SYSTEM/FileConfig.cpp` and `SYSTEM/FileTemp.cpp`). It is the
highest-fan-out header in `SYSTEM`, with 53 direct TOPP consumers and 35
class-test sections.

The whole public surface is covered: 39 source function names in 40 overloads,
the two `TempDir` constructors and its destructor, both public enumerations, and
the private helpers that are observable through them. Ordinary filesystem work
uses safe Rust's standard library; no C++/Qt ABI, no `dlfcn`/`mach-o`/`Windows.h`
probe, no process-exit destructor registry, and no implicit `.reference`
checkout is required.

Two structural decisions shape the mapping and are explained under
[Explicit runtime context](#explicit-runtime-context) and
[Owned temporary lifetimes](#owned-temporary-lifetimes): everything the source
derives from process-global environment state lives in a caller-owned
`FileContext`, and every temporary resource is an owned guard rather than a name
in an exit-time list.

## API mapping

Every public member of `SYSTEM/File.h`, in header order.

### `File::TempDir`

| Source member | Rust counterpart | Note |
|---|---|---|
| `TempDir(bool keep_dir = false)` | `TempDir::new(keep_dir)` | base from `FileContext::from_environment().get_temp_directory()` |
| `TempDir(const std::string& base_dir, bool keep_dir = false)` | `TempDir::new_in(base, keep_dir)` | creates the parent chain, then an exclusive leaf |
| `~TempDir()` | `impl Drop for TempDir` | removes the tree unless `keep_dir`; cannot report a failure |
| `TempDir(const TempDir&) = delete`, `operator=(const TempDir&) = delete` | no `Clone` impl | copying is simply not expressible |
| `TempDir(TempDir&&) = delete`, `operator=(TempDir&&) = delete` | *not ported as a restriction* | a Rust move transfers the single owner and cannot duplicate cleanup, so forbidding it would buy nothing |
| `getPath()` | `TempDir::path` | borrowed `&Path`, no trailing separator |
| — | `TempDir::keep` | native: give up ownership later, the `keep_dir = true` outcome decided after construction |
| — | `TempDir::close` | native: remove now and report the failure `Drop` cannot |

### `File` statics

| Source member | Rust counterpart | Note |
|---|---|---|
| `getExecutablePath()` | `get_executable_path`, `FileContext::get_executable_path` | `Result<PathBuf>` instead of "empty string on failure"; no trailing `/` |
| `exists(file)` | `exists` | |
| `empty(file)` | `empty` | missing paths, directories and metadata failures all count as empty, as in the source |
| `executable(file)` | `executable` | mode bits on Unix, existence elsewhere |
| `fileSize(file)` | `file_size` | `Result<u64>`; the source's `-1` from a `UInt64` function is removed |
| `getModificationTime(file)` | `get_modification_time` | whole seconds since the Unix epoch, negative before it |
| `rename(from, to, overwrite_existing, verbose)` | `rename(from, to, overwrite)` | `verbose` selected only whether the failure was logged; the failure is returned instead |
| `enum class CopyOptions` | `CopyOptions` | `Overwrite` (default), `Skip`, `Cancel` |
| `copyDirRecursively(from, to, option)` | `copy_dir_recursively` | |
| `copy(from, to)` | `copy` | refuses an existing destination, as `copy_file` does |
| `remove(file)` | `remove` | file, symlink or *empty* directory, matching `std::remove` |
| `removeDirRecursively(dir)` | `remove_dir_recursively` | |
| `removeDir(dir)` | `remove_dir` | a synonym: the source is `remove_all` under both names |
| `makeDir(dir)` | `make_dir` | an existing directory is success |
| `absolutePath(file)` | `absolute_path` | lexical, as `fs::absolute` is; `""` is the working directory |
| `basename(file)` | `basename` | |
| `stemName(file)` | `stem_name` | compound extensions via `crate::format::file_types` |
| `extension(file)` | `extension` | |
| `listDirectories(dir)` | `list_directories` | `Result<Vec<PathBuf>>`; the source's error-as-empty-list is an error here |
| `path(file)` | `path` | `"."` for a bare filename, as the source requires |
| `readable(file)` | `readable` | effective-process probe, not real-uid `access(2)` |
| `writable(file)` | `writable` | non-destructive probe with the source's short-name fallback |
| `isDirectory(path)` | `is_directory` | |
| `find(filename, directories)` | `FileContext::find` | same order, same lexical normalisation, same diagnostic |
| `fileList(dir, pattern, output, full_path)` | `file_list` | returns the list; the source's `bool` is "list is nonempty" |
| `findDoc(filename)` | `FileContext::find_doc` | compiled-in doc paths become `documentation_directories` |
| `getUniqueName(include_hostname)` | `get_unique_name` | UTC, `HOSTNAME`/`COMPUTERNAME` instead of `gethostname` |
| `getOpenMSDataPath()` | `FileContext::get_openms_data_path` | see [the search order](#the-data-path-search-order) |
| `getOpenMSDataPathSource()` | `FileContext::get_openms_data_path_source` | |
| `getOpenMSHomePath()` | `FileContext::get_openms_home_path` | snapshotted, so infallible |
| `getOpenMSConfigDir()` | `FileContext::get_openms_config_dir` | the source's `__unix__` branch, not Rust's `unix` flag |
| `getTempDirectory()` | `FileContext::get_temp_directory` | |
| `getUserDirectory()` | `FileContext::get_user_directory` | reads `OPENMS_HOME_PATH`, as the implementation does |
| `getSystemParameters()` | `FileContext::get_system_parameters`, `FileContext::get_system_parameters_with_warnings` | the second returns what the source writes to its warning log |
| `findDatabase(db_name)` | `FileContext::find_database` | |
| `getPathLocations()` | `get_path_locations_from_environment` | |
| `getPathLocations(path)` | `get_path_locations` | |
| `findExecutable(exe_filename)` | `FileContext::find_executable` | `Result<Option<PathBuf>>` instead of an in/out parameter plus `bool` |
| `findSiblingTOPPExecutable(toolName)` | `FileContext::find_sibling_topp_executable` | siblings only; `PATH` is deliberately not consulted |
| `getTemporaryFile(alternative_file)` | `FileContext::get_temporary_file`, `TempFile::new`, `TempFile::new_in`, `TempFile::alternative` | |
| `enum class MatchingFileListsStatus` | `MatchingFileListsStatus` | discriminants 0, 1, 2 preserved |
| `validateMatchingFileNames(sl1, sl2, basename, ignore_extension)` | `validate_matching_file_names` | |

### `File` private members observable through the public surface

| Source member | Rust counterpart | Note |
|---|---|---|
| `struct OpenMSDataPath_` | `ResolvedDataPath` | public here, because the context that supplies candidates needs it |
| `resolveOpenMSDataPath_()` | `FileContext::resolve_data_path` | public here; also `FileContext::clear_data_cache`, which the source's static cannot offer |
| `isOpenMSDataPath_(path)` | private: the `CHEMISTRY/unimod.xml` marker check | same marker file |
| `getSystemParameterDefaults_()` | private `system_parameter_defaults` | the same five entries |
| `executableExtensions_()` / `executableExtensions_(ext)` | `FileContext::executable_extensions`, filled by `from_environment` | the `%PATHEXT%` sanity check and the `.exe`/`.bat` fallback are preserved |
| `class TemporaryFiles_` and the `temporary_files_` static | **not ported** | the exit-time registry of uncreated *names* is replaced by owned guards; see below |
| `coreLibraryDirectory()` (anonymous namespace) | **not ported** | `dladdr`/`GetModuleHandleExW` have no safe-Rust equivalent; an installer adds the library-relative directory to `data_candidates` explicitly |

Native additions with no source member: `TempFile` and its `keep`/`close`,
`TempDir::keep`/`close`, `FileContext::new`/`from_environment`/`clear_data_cache`,
`get_system_parameters_with_warnings`, and the five `MAX_*` ceilings.

## Preserved source conventions

* `basename("/path/only/")` is `""`; `path("filename_only.h")` is `"."`;
  `path("/a.txt")` is `""`. Both `/` and `\` separate on every platform.
* `stemName`/`extension` recognise compound extensions (`.mzML.gz`) and are not
  fooled by a dot in a directory name; `.mzML` is all extension.
* `empty` is true for a missing path, a directory, and a metadata failure.
* `remove` succeeds on a path that is not there, and removes an empty directory.
* `makeDir` succeeds when the directory already exists.
* `removeDir` is recursive, exactly like `removeDirRecursively`.
* `copy` fails when the destination exists.
* `rename` treats "same file after resolving symlinks" as success with nothing
  done, and falls back to copy-then-remove across devices.
* `find` returns an already existing input unchanged, resolves the shared data
  *before* searching any caller-supplied directory, searches caller directories
  ahead of the data directory, lexically normalises the hit, and refuses an
  empty or whitespace-only name.
* `getPathLocations("")` yields no entries, while an empty component yields `/`.
* `findExecutable` checks existence, not execute permission, and a path that
  already exists wins without consulting `PATH`.
* `findSiblingTOPPExecutable` searches only the executable directory and the
  three macOS bundle alternatives; a binary reachable only through `PATH` is
  reported as not found, which the class test pins.
* `getUserDirectory` reads `OPENMS_HOME_PATH` — the header's `OPENMS_HOME_DIR`
  is a documentation defect — and its `ensureLastChar` turns an explicitly empty
  value into `/`.
* `getTempDirectory` and `getUserDirectory` load `OpenMS.ini` before looking at
  their environment override, so a malformed configuration is an error even when
  the override would have answered.
* `getSystemParameters` returns the file's own tree when the version is stale,
  with only `version` overwritten, and discards the repaired defaults tree it
  just built. Missing defaults stay missing, and no file is rewritten.
* `getOpenMSConfigDir` uses XDG only where the compiler defines `__unix__`.
  Apple's compilers do not, so macOS is `<home>/.OpenMS`; the port switches on
  `cfg!(all(unix, not(any(target_os = "macos", target_os = "ios"))))` rather
  than on Rust's `unix`, which would put macOS on the wrong branch.
* `validateMatchingFileNames` compares `std::set`s, so multiplicity is
  discarded: `["a","a","b"]` against `["a","b","b"]` is `ORDER_MISMATCH`.
* `MatchingFileListsStatus` keeps the explicit discriminants 0, 1 and 2.

## The data-path search order

`getOpenMSDataPath` is the riskiest member to get wrong, because a wrong answer
silently loads another installation's chemistry tables.

1. `OPENMS_DATA_PATH`, when set. **Authoritative**: an invalid override is an
   error and does *not* fall through, so a stale variable cannot select the
   shared data of a different installed version.
2. Otherwise, each candidate in order. The source's candidates are: the loaded
   library's directory (via `dladdr`), the executable's directory, the macOS app
   bundle, then two directories compiled into the library (`OPENMS_DATA_PATH`
   and `OPENMS_INSTALL_DATA_PATH`).

The port's defaults are the executable-relative `../share/OpenMS` and, on macOS,
the bundle's `../../../share/OpenMS`. The loaded-library probe and the two
compiled-in paths have no safe-Rust equivalent and no meaning in a crate that is
not built by the SDK's CMake; an installer or embedder inserts them into
`FileContext::data_candidates` in the order they should be tried, which is
strictly more expressive than the source's fixed list.

A candidate is valid when it holds `CHEMISTRY/unimod.xml`, which is exactly
`isOpenMSDataPath_`. The chemistry tables packaged inside this crate
deliberately do not satisfy that check: they are not a source-layout data tree
and must not pass for one.

**When none exists**, `resolve_data_path` returns `Error::Io` (not found) whose
message names the override or the exhausted candidates and points at
`OPENMS_DATA_PATH`, mirroring the source's `Exception::FileNotFound` diagnostic.
Only a *successful* resolution is cached, as the source's function-local static
is; a failure can be retried after correcting the context, where the source's
static would keep throwing for the life of the process. A cached success
survives later edits to the context's fields until `clear_data_cache`.

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

`FileContext::new(executable_directory, home_directory, temporary_directory)`
builds an isolated context that reads no environment at all;
`from_environment` snapshots `OPENMS_HOME_PATH`, `HOME`/`USERPROFILE`,
`OPENMS_TMPDIR`, `OPENMS_DATA_PATH`, `XDG_CONFIG_HOME`, `PATH`, `%PATHEXT%` and
the current executable, once. Nothing reads the environment again afterwards.
That is what lets the class test's environment-mutating sections
(`getUserDirectory`, `getOpenMSConfigDir`) be expressed without changing state
shared by every thread in the process, and it is why parallel tests here are
independent.

Configuration uses [Param](PARAM_SUPPORT.md) and [ParamXML](PARAMXML_SUPPORT.md).
A missing or unreadable `OpenMS.ini` yields the five defaults — SDK version,
empty `home_dir`, empty `temp_dir`, empty `id_db_dir` list, `threads = 1` — and
creates nothing. An existing readable configuration needs the `paramxml`
feature; without it the call returns `Error::Unsupported`, which is an explicit
refusal rather than a silent fall back to defaults that would disagree with the
file on disk.

## Owned temporary lifetimes

The source's `getTemporaryFile` returns a *name*, creates nothing, and appends
the name to a process-wide list that a static destructor walks at exit. Nothing
in that design says when a temporary stopped being needed, and the name stays
free for another process to take between the call and the first write.

`TempFile` creates the file immediately and exclusively (`0600` on Unix) and
removes it on `Drop`. `TempFile::alternative`, and a nonempty `alternative` path
passed to `FileContext::get_temporary_file`, own nothing: the path is neither
created nor removed, which is the source's "return it and do not destroy it"
branch. `TempDir` creates an exclusive directory including missing parents, and
removes the tree on `Drop` unless it was built with `keep_dir = true` or
consumed by `keep`.

**Drop cannot fail.** A removal that fails during `Drop` is silent. The source is
silent too — `removeDirRecursively` prints to `stderr`, returns a `bool` the
destructor discards, and throws nothing — but here that is a decision, so both
guards also offer `close`, which removes now and returns the error; a removal
that fails in `close` is attempted once more on drop. Neither guard runs at all
when the process ends without unwinding, which is the one case the source's
exit-time registry covered and this does not: retain guards for as long as
subprocesses or other consumers use their paths.

`TempDir::new_in` produces `<base>/OpenMSTempDir_<unique name>`, matching the
source's two-argument constructor but not its no-argument one, which omits the
`OpenMSTempDir_` prefix; the `_XXXXXX` that `mkdtemp` consumes is replaced by an
exclusive `mkdir` with retry. No caller should depend on either spelling.

Unique names are UTC Gregorian date and time, an optional sanitized host, the
process id and an atomic counter starting at 1. The host comes from
`HOSTNAME`/`COMPUTERNAME` rather than `gethostname`, which safe Rust cannot call
without a dependency; that variable is normally unset for a non-interactive
process, so unlike the source the host part may be absent and
`get_unique_name(true)` need not be longer than `get_unique_name(false)`. The
names are identifiers, not security credentials: exclusive creation is what
makes a temporary safe, not the name.

## Native differences

Filesystem results are `PathBuf`/`Path`, so native path encodings survive rather
than requiring UTF-8, and paths never acquire a decorative trailing separator —
join child names with `Path::join`. Simple predicates stay `bool`; everything
else returns a checked `Result`, so a directory-list failure, a missing file's
size, and a pre-epoch modification time are all distinguishable instead of
sharing a sentinel with a legitimate value.

**Copying, moving and removal.** `copy` streams into a sibling temporary,
preserves the source permissions, synchronises, and publishes with an atomic
no-clobber hard link. `rename(..., false)` publishes regular files with the same
guarantee, so concurrent movers cannot replace the winner; same-filesystem moves
use a hard link then remove the source, cross-device moves stage a streamed copy
on the destination filesystem. A failed source removal is reported even when the
destination is already complete. `rename(..., true)` uses the native rename,
which corrects the source's delete-before-rename and its POSIX no-clobber flag
discrepancy. Cross-device directory and symlink moves, and no-clobber directory
and symlink moves, are checked errors: portable safe `std` has no atomic
no-clobber directory rename, and this port will not substitute a racy existence
check. A filesystem that cannot publish hard links returns an error rather than
a weaker guarantee.

`copy_dir_recursively` plans a bounded traversal, rejects copying onto the
source or a descendant of it, detects source directory-link cycles, and checks
every visible conflict before modifying any destination — so `Cancel` improves
on the source's filesystem-order partial updates. Source directory symlinks are
followed for copying; destination directory symlinks are rejected, so a planned
write cannot be redirected. Recursive copy and removal are not filesystem
transactions: each file publication is complete, but a tree does not become
visible in one atomic step, and an error partway through leaves earlier work.
Recursive removal plans first and never follows symlinks, so a link inside the
tree is unlinked and its target is untouched; dangling links are really removed.
Hostile replacement of a directory during a traversal is outside what any
path-based standard-library operation can guarantee, here and in the source.

**Permission queries.** These are effective-process semantic probes, not exact
real-uid `access(2)` or Windows ACL calculations. Existing regular files are
opened without truncation; directories and absent targets use an exclusive
sibling probe that never touches the requested basename. FIFOs and devices are
not opened by a query. The probe's own name is longer than a typical caller's,
so near the platform's path limit there is a band in which the caller's file
still fits and the probe does not; the source answers that by asking the
directory itself with `access(2)`, and this port by falling back to
progressively shorter probe names, which `writable_agrees_with_the_operating_system_at_every_path_depth`
checks against the operating system's own answer at every reachable depth. The
ladder narrows the band rather than closing it: its shortest rung is still a
name — the bare decimal probe counter, one byte until the process has handed out
ten unique names and two thereafter — where `access(2)` needs none, so a
directory within that many bytes of the limit is still answered "not writable"
where the source answers "writable".

**Wildcards.** Matching is deterministic, case-sensitive Unicode: `*`, `?`,
bracket ranges and classes with `!`/`^` negation, and backslash escapes. A
leading dot is an ordinary character, as in POSIX `fnmatch` with a zero flag
word. Three source behaviours are not emulated: POSIX locale character classes
and collation; the Windows `PathMatchSpecA` pattern-list and case rules; and
byte-wise matching, since `?` here covers one character where `fnmatch` covers
one byte. A directory entry whose name is not valid UTF-8 makes `file_list`
return an error naming that entry, rather than skipping it or matching it
approximately — both of which would answer a question the caller did not ask.
Everything that does not need to decode a name (`exists`, `file_size`,
`list_directories`, recursive removal) works on such entries unchanged.

**Non-UTF-8 and multi-byte paths.** The four lexical helpers keep the source's
string signature and return borrowed subslices cut at a separator or at the
stem's own length — both character boundaries — so a name such as
`dir/日本語.txt` is split correctly and no computed byte offset can panic. The
`lexical_helpers_never_split_a_multibyte_name` and
`a_non_utf8_directory_entry_is_named_in_the_error_rather_than_skipped` tests pin
both halves.

**Platform coverage.** `File` is platform-shaped in C++ (`Windows.h`,
`Shlwapi.h`, `dlfcn.h`, `mach-o/dyld.h`); this crate is Linux-first and portable.
`executable` reports real permission bits on Unix and existence elsewhere, as
the source does. `get_executable_path` answers on every platform `std` supports.
The macOS-only candidates (`app bundle data`, the three
`findSiblingTOPPExecutable` bundle locations, `<data>/../../Documentation`) are
compiled in under `cfg!(target_os = "macos")`, so they exist on macOS and
nowhere else, as in the source. `executable_extensions` and the `;` PATH
separator apply on Windows only. The `dladdr` library-directory probe has no
portable equivalent and is replaced by explicit candidates.

**Concurrency.** Neither `File.cpp`, `FileConfig.cpp` nor `FileTemp.cpp` carries
a `#pragma omp`; the source's only concurrency is the mutex guarding its
temporary-name list and the thread-safe static initialisation of its caches.
There is nothing to parallelise here, so this module is serial and uses no
rayon. `FileContext` is `Send`/`Sync`-friendly plain data whose cache is a
`OnceLock`, and the concurrency tests exercise the probe and the no-clobber
publication from many threads.

## Checked boundaries and evidence

Directory, search and matching work is charged against `MAX_ENTRIES` (100,000
entry visits), `MAX_WORK` (50 million units), `MAX_BYTES` (64 MiB of
conservative logical allocation), `MAX_PATH_BYTES` (1 MiB per path or pattern)
and `MAX_DEPTH` (128 levels). Entry accounting includes planning and enumeration
visits, not just output rows, so filtering does not make a traversal cheap.
Temporary creation retries at most 100 times, matching the source's own cap.
File contents are streamed, with no small scientific-data-size ceiling. Param and
ParamXML retain their own checked parsing and update budgets. Every ceiling is
checked in a preflight before anything is allocated or mutated, so a rejected
call leaves the filesystem unchanged — `native_resource_bounds_are_checked_before_directory_mutation`
asserts that for the depth ceiling specifically. Ordinary caller-owned context
and guard values, and Rust `Clone`, are not secretly fallible.

Evidence is **tier 3 and tier 4**, and no tier 1 or tier 2 claim is made: no C++
was built or executed for this group, and no retained C++ output or oracle
driver exists for `File`, whose answers are properties of the machine it runs on
rather than of a fixture. Tier 3 is the transcription of the class test's
literals, section by section; tier 4 covers the native concurrency, lifecycle,
link, encoding and bound guards, and the two cases derived from the platform
rather than from a literal — the `__unix__` configuration-directory branch, and
the path-depth sweep, which compares `writable` against the operating system's
own answer instead of against a hard-coded length.

All 35 class-test sections are mapped in
[tests/system_file.rs](../tests/system_file.rs) (27 integration tests) plus two
private tests in the module. Two sections cannot be reproduced literally and say
so at the test: the source's `makeDir` relative-path case changes the process
working directory, which no test in a multithreaded binary may do, and the
source's `getUniqueName` assertion that the hostname form is strictly longer
does not hold for an environment-derived host that may be absent.
[tests/data/system_file_provenance.json](../tests/data/system_file_provenance.json)
records the hashes of every source file read, the packaged fixtures, the
per-quirk source anchors and the section map.
