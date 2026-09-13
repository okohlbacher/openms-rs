// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Basic file handling operations, from the Core SDK `SYSTEM/File.h`.
//!
//! The source `File` class is a collection of static members over
//! `std::string` paths: existence and type predicates, lexical path splitting,
//! copying, moving and removal, directory listing with a glob filter, the
//! search order that locates the installed shared-data tree, executable lookup
//! along `PATH`, and the `TempDir` guard. It is the highest-fan-out header in
//! `SYSTEM`, with 53 direct TOPP consumers.
//!
//! Three groups of members map onto Rust differently from the rest:
//!
//! * **Runtime configuration.** The source reads `OPENMS_DATA_PATH`, the user
//!   home and the temporary directory from process-global environment
//!   variables, compiled-in paths and a `dladdr` probe of the loaded library,
//!   and caches the answer in a function-local static. This port puts all of it
//!   in a caller-owned [`FileContext`](crate::system::file::FileContext), which
//!   [`FileContext::from_environment`](crate::system::file::FileContext::from_environment)
//!   fills from one snapshot of the environment. Tests therefore never mutate
//!   the environment of the whole process, and two contexts may disagree.
//! * **Temporary resources.** The source registers temporary *names* in a
//!   process-wide list and deletes them from a static destructor at exit. This
//!   port returns owned guards, [`TempFile`](crate::system::file::TempFile) and
//!   [`TempDir`](crate::system::file::TempDir), which clean up on `Drop`.
//! * **Platform-specific members.** `File` includes `Windows.h`, `Shlwapi.h`,
//!   `dlfcn.h` and `mach-o/dyld.h`. This crate is Linux-first and uses only
//!   [`std::fs`] and [`std::path`]; where a member is inherently
//!   platform-shaped, the item says which platforms return a real answer.
//!
//! Paths are [`Path`](std::path::Path) and [`PathBuf`](std::path::PathBuf)
//! throughout, so a name that is not valid UTF-8 is carried unchanged rather
//! than lost. The four lexical helpers
//! [`basename`](crate::system::file::basename),
//! [`path`](crate::system::file::path),
//! [`stem_name`](crate::system::file::stem_name) and
//! [`extension`](crate::system::file::extension) keep the source's string
//! signature, because they are string operations on a string; they return
//! borrowed subslices cut at a separator or at the stem's own length, never at
//! a computed byte offset, so a multi-byte name such as `dir/日本語.txt` cannot
//! panic.
//!
//! See `docs/SYSTEM_FILE_SUPPORT.md` for the full member table, the data-path
//! search order and the checked resource bounds.
use crate::format::file_types::strip_extension;
use crate::param::{Param, ParamValue};
use crate::{CORE_SDK_VERSION, Error, Result};
use std::ffi::OsString;
use std::fs::{self, File, OpenOptions};
use std::io;
use std::path::{Component, Path, PathBuf};
use std::sync::{
    OnceLock,
    atomic::{AtomicU64, Ordering},
};
use std::time::{SystemTime, UNIX_EPOCH};

/// Largest number of directory entries one operation will charge for.
///
/// Every planning and enumeration visit is charged, not only the rows that
/// reach the result, so a traversal cannot be made cheap by filtering. The
/// source imposes no bound at all.
pub const MAX_ENTRIES: usize = 100_000;
/// Largest accepted path or wildcard pattern, in bytes.
pub const MAX_PATH_BYTES: usize = 1024 * 1024;
/// Work budget for one operation, in abstract units.
///
/// Charged for path inspection, per-entry visits, sorting, and each character
/// step of a wildcard match, so a pathological pattern such as `*a*a*a*a*b`
/// against a long name is refused instead of being run to completion.
pub const MAX_WORK: usize = 50_000_000;
/// Conservative ceiling on the logical bytes one operation will allocate.
///
/// File *contents* are streamed and are not charged here; this bounds the
/// bookkeeping — paths, copy plans, pattern tokens — an operation builds in
/// memory.
pub const MAX_BYTES: usize = 64 * 1024 * 1024;
/// Largest directory nesting one operation will create or walk.
pub const MAX_DEPTH: usize = 128;

/// Whether the source's `__unix__` branch applies to this target.
///
/// `getOpenMSConfigDir` switches on `__unix__`, which Apple's compilers do not
/// define; macOS therefore takes the non-unix branch in the source. Rust's own
/// `unix` configuration flag *is* set on macOS, so using it directly would move
/// macOS onto the XDG branch and change where `OpenMS.ini` is looked for.
const SOURCE_UNIX: bool = cfg!(all(unix, not(any(target_os = "macos", target_os = "ios"))));

/// What [`copy_dir_recursively`] does with a file that already exists in the target.
///
/// Source `File::CopyOptions`, whose default is `OVERWRITE`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CopyOptions {
    /// Replace the file already in the target. The source's default.
    #[default]
    Overwrite,
    /// Leave the file already in the target alone.
    ///
    /// The source announces every skipped file on its warning log; this port is
    /// silent, because a library writing to a log the caller does not own is
    /// what makes a passing test look like a broken build.
    Skip,
    /// Refuse the copy when any file is already in the target.
    ///
    /// The source returns `false` at the first conflict it happens to reach, in
    /// filesystem order, having already copied everything before it. This port
    /// finds every conflict it can see in a preflight and then changes nothing.
    Cancel,
}
/// Outcome of comparing two lists of filenames with [`validate_matching_file_names`].
///
/// Source `File::MatchingFileListsStatus`, whose enumerators carry the explicit
/// values 0, 1 and 2. The discriminants are kept, so a caller that reports the
/// number reports the same number.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum MatchingFileListsStatus {
    /// The same names, in the same order.
    Match = 0,
    /// The same set of names, in a different order.
    OrderMismatch = 1,
    /// Different names, or a different count.
    SetMismatch = 2,
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn absent(message: impl Into<String>) -> Error {
    io::Error::new(io::ErrorKind::NotFound, message.into()).into()
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b)
        .ok_or_else(|| invalid("filesystem size overflow"))
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b)
        .ok_or_else(|| invalid("filesystem size overflow"))
}
struct Work {
    work: usize,
    bytes: usize,
    entries: usize,
}
impl Default for Work {
    fn default() -> Self {
        Self {
            work: MAX_WORK,
            bytes: MAX_BYTES,
            entries: MAX_ENTRIES,
        }
    }
}
impl Work {
    fn consume(&mut self, n: usize) -> Result<()> {
        self.work = self
            .work
            .checked_sub(n)
            .ok_or_else(|| invalid("filesystem work limit exceeded"))?;
        Ok(())
    }
    fn copy(&mut self, n: usize) -> Result<()> {
        self.consume(n)?;
        self.bytes = self
            .bytes
            .checked_sub(n)
            .ok_or_else(|| invalid("filesystem allocation limit exceeded"))?;
        Ok(())
    }
    fn path(&mut self, p: &Path) -> Result<()> {
        check_path(p)?;
        self.consume(p.as_os_str().len())
    }
    fn owned(&mut self, p: &Path) -> Result<PathBuf> {
        self.path(p)?;
        self.copy(add(p.as_os_str().len(), 32)?)?;
        Ok(p.to_owned())
    }
    fn entry(&mut self) -> Result<()> {
        self.entries = self
            .entries
            .checked_sub(1)
            .ok_or_else(|| invalid("filesystem entry limit exceeded"))?;
        self.copy(128)
    }
}
fn check_path(p: &Path) -> Result<()> {
    if p.as_os_str().len() > MAX_PATH_BYTES {
        return Err(invalid("filesystem path limit exceeded"));
    }
    if p.as_os_str().as_encoded_bytes().contains(&0) {
        return Err(invalid("filesystem path contains NUL"));
    }
    Ok(())
}
fn checked_paths(a: &Path, b: &Path) -> Result<()> {
    check_path(a)?;
    check_path(b)
}

/// Test whether `file` exists.
///
/// Source `File::exists`, which is `std::filesystem::exists` and answers
/// `false` for an empty path. A path this port refuses outright — longer than
/// [`MAX_PATH_BYTES`], or holding an interior NUL — is `false` as well, and so
/// is one whose existence cannot be determined because a parent directory
/// cannot be searched, where the source's throwing overload raises
/// `filesystem_error`: a predicate has nowhere to report an error, and "not a
/// usable path" is not "exists". [`is_directory`] answers the same way.
pub fn exists(file: impl AsRef<Path>) -> bool {
    check_path(file.as_ref()).is_ok() && file.as_ref().try_exists().unwrap_or(false)
}
/// Test whether `file` is missing or holds no bytes.
///
/// Source `File::empty` answers `true` when the path does not exist, when the
/// size query fails — which it does for a directory — and when the size is
/// zero. All three are preserved, so this reads as "not a nonempty regular
/// file" rather than "empty file". Use [`exists`] to tell a missing path from
/// an empty one and [`file_size`] for a checked size.
pub fn empty(file: impl AsRef<Path>) -> bool {
    check_path(file.as_ref()).is_err()
        || fs::metadata(file).map_or(true, |m| !m.is_file() || m.len() == 0)
}
/// Test whether `file` names a directory.
///
/// Source `File::isDirectory`. Symbolic links are followed, as
/// `std::filesystem::is_directory` follows them, so a link to a directory is a
/// directory. An empty path and an unreadable parent are both `false`.
pub fn is_directory(file: impl AsRef<Path>) -> bool {
    check_path(file.as_ref()).is_ok() && file.as_ref().is_dir()
}
/// The size of `file` in bytes.
///
/// # Errors
///
/// Returns [`Error::Io`] when `file` does not exist or its metadata cannot be
/// read, and [`Error::InvalidValue`] when it is not a regular file. Source
/// `File::fileSize` returns `-1` for both — from a `UInt64` function, so the
/// failure reaches the caller as `18446744073709551615` and compares greater
/// than every real size. The `Result` removes the sentinel.
pub fn file_size(file: impl AsRef<Path>) -> Result<u64> {
    check_path(file.as_ref())?;
    let m = fs::metadata(file)?;
    if !m.is_file() {
        return Err(invalid("file size requires a regular file"));
    }
    Ok(m.len())
}
/// Last modification time of `file`, in whole seconds since the Unix epoch.
///
/// Source `File::getModificationTime` uses `stat()` rather than
/// `std::filesystem::last_write_time`, deliberately: `file_time_type`'s clock
/// epoch is implementation-defined — 2174 on libstdc++, 1601 on MSVC — and two
/// standard libraries therefore report different numbers for the same file,
/// which matters as soon as the value is written out or compared against one
/// another machine recorded. This port returns the same quantity.
///
/// Times before the epoch are negative and truncate towards negative infinity,
/// matching POSIX `st_mtim.tv_sec`: a file stamped half a second before the
/// epoch reports `-1`, not `0`.
///
/// # Errors
///
/// Returns [`Error::Io`] when `file` does not exist or its metadata cannot be
/// read, and [`Error::InvalidValue`] when the time does not fit an `i64`. The
/// source returns `-1` for a failure, which is also a legitimate time — one
/// second before the epoch — so the two cannot be told apart there; the
/// `Result` separates them.
pub fn get_modification_time(file: impl AsRef<Path>) -> Result<i64> {
    check_path(file.as_ref())?;
    let t = fs::metadata(file)?.modified()?;
    match t.duration_since(UNIX_EPOCH) {
        Ok(d) => i64::try_from(d.as_secs()).map_err(|_| invalid("modification time overflow")),
        Err(e) => {
            let d = e.duration();
            let seconds = d
                .as_secs()
                .checked_add(u64::from(d.subsec_nanos() != 0))
                .ok_or_else(|| invalid("modification time overflow"))?;
            i64::try_from(seconds)
                .ok()
                .and_then(i64::checked_neg)
                .ok_or_else(|| invalid("modification time overflow"))
        }
    }
}
/// Test whether `file` has an execute permission bit set.
///
/// Source `File::executable` tests `owner_exec | group_exec | others_exec` on
/// POSIX and returns "the file exists" on Windows, where it comments that all
/// files count as executable; both are reproduced, and a non-Unix target here
/// takes the second branch. Symbolic links are followed.
///
/// The bits are read from the file's mode, not from an access check for the
/// calling user, so a file this process could not in fact execute still answers
/// `true` — the same answer the source gives.
pub fn executable(file: impl AsRef<Path>) -> bool {
    if check_path(file.as_ref()).is_err() {
        return false;
    }
    let Ok(m) = fs::metadata(file) else {
        return false;
    };
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        m.permissions().mode() & 0o111 != 0
    }
    #[cfg(not(unix))]
    {
        let _ = m;
        true
    }
}
/// Test whether `file` exists and can be read.
///
/// Source `File::readable` asks `access(file, R_OK)`, which answers for the
/// *real* user id. This port opens the file, or enumerates the directory,
/// which answers for the effective process; portable safe Rust has no real-uid
/// query. The two differ only for a set-uid process.
///
/// Anything that is neither a regular file nor a directory — a FIFO, a device —
/// answers `false` without being opened, so asking the question cannot block on
/// a pipe with no writer. An empty path is `false`, as in the source.
pub fn readable(file: impl AsRef<Path>) -> bool {
    let p = file.as_ref();
    if check_path(p).is_err() {
        return false;
    }
    match fs::metadata(p) {
        Ok(m) if m.is_dir() => fs::read_dir(p).is_ok(),
        Ok(m) if m.is_file() => File::open(p).is_ok(),
        _ => false, // Do not block by opening a FIFO or device during a query.
    }
}
/// Test whether `file` could be written, leaving `file` itself untouched.
///
/// This is a query, and the source is explicit that asking it must not change
/// the answer: an existing file is opened for writing but never truncated, and
/// a path that does not exist is answered by creating a *probe* file of this
/// port's own in the parent directory and removing that again — never by
/// creating and deleting the caller's own name. Probing under the caller's name
/// is what used to make two TOPP tools given the same `-out` path delete each
/// other's output and then report it unwritable.
///
/// An empty path is `false`, as in the source. A path whose parent does not
/// exist is `false`. A directory is writable when a new entry can be created in
/// it. Anything that exists but is neither a regular file nor a directory — a
/// FIFO, a device — is `false` without being opened; the source's `access(2)`
/// answers `true` for a writable FIFO, but opening one for writing blocks until
/// a reader appears and a query must not block. [`readable`] refuses the same
/// kinds for the same reason.
///
/// The probe has no single fixed name. [`TempFile`] walks a ladder of
/// progressively shorter candidates, ending at the bare decimal probe counter,
/// and the first that can be created exclusively answers the question. The
/// ladder stands in for the source's own fallback: the source tries one
/// descriptive probe name and, when the create fails for a reason that is about
/// the *name* rather than about the directory — `ENAMETOOLONG`, or the
/// `ENOENT`/`EINVAL` that the Windows CRT folds `ERROR_FILENAME_EXCED_RANGE`
/// onto — asks the directory itself with `access(dir, W_OK)`. Without one or
/// the other, a directory close enough to the platform's path limit that the
/// caller's shorter name would still fit but the descriptive probe's would not
/// is reported unwritable, which is exactly the false negative this function
/// exists to avoid; `writable_agrees_with_the_operating_system_at_every_path_depth`
/// in `tests/system_file.rs` sweeps that band against the operating system's
/// own answer, and fails on Linux without the ladder.
///
/// Three divergences from the source remain.
///
/// * `access(2)` needs no filename at all, while the shortest rung of the
///   ladder is still a name: the decimal counter, one byte until the process
///   has handed out ten unique names and two thereafter. A directory whose path
///   sits within that many bytes of the platform's limit — close enough that
///   the caller's own one-character name fits and the counter does not — is
///   answered `false` where the source answers `true`. Safe Rust has no
///   `access(2)`, so the band is narrowed rather than closed.
/// * The ladder is walked on *any* failed create, not on a decoded
///   `ENAMETOOLONG`, because `io::ErrorKind::InvalidFilename` is unstable at
///   this crate's minimum Rust version. A directory that refuses the first
///   candidate for a reason that has nothing to do with the name therefore
///   costs two more `open` calls before the same `false` is returned.
/// * Creating answers for the *effective* process where `access(2)` answers for
///   the *real* user id, as [`readable`] notes; the two differ only for a
///   set-uid process. For an existing directory that also means this port
///   creates and removes a probe inside it where the source creates nothing:
///   the directory is left holding exactly what it held, but its modification
///   time moves and a filesystem watcher sees the two events.
pub fn writable(file: impl AsRef<Path>) -> bool {
    let p = file.as_ref();
    if p.as_os_str().is_empty() || check_path(p).is_err() {
        return false;
    }
    match fs::metadata(p) {
        Ok(m) if m.is_file() => OpenOptions::new().write(true).open(p).is_ok(),
        Ok(m) if m.is_dir() => probe_directory(p, None),
        Ok(_) => false,
        Err(e) if e.kind() == io::ErrorKind::NotFound => probe_directory(parent(p), Some(p)),
        Err(_) => false,
    }
}
fn parent(p: &Path) -> &Path {
    p.parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
}
fn probe_directory(p: &Path, avoid: Option<&Path>) -> bool {
    TempFile::create_in(p, avoid)
        .and_then(|(guard, file)| {
            drop(file);
            guard.close()
        })
        .is_ok()
}

/// The final component of `file`, without any directory part.
///
/// Source `File::basename`, by way of `PathUtils::basename`. Nothing is checked
/// on the filesystem: `/path/some_entity` yields `some_entity` whether that is
/// a file or a directory, and a trailing separator yields the empty string, so
/// `/path/only/` is `""`. Both separators are recognised on every platform, as
/// in the source.
///
/// The result borrows from `file` and is cut at a separator, so a multi-byte
/// component comes back whole.
pub fn basename(file: &str) -> &str {
    file.rsplit(['/', '\\']).next().unwrap_or("")
}
/// The directory part of `file`, without a trailing separator.
///
/// Source `File::path`. A name with no separator yields `"."` rather than an
/// empty string, and the source's comment says why: generic code writes
/// `path(f) + '/' + basename(f)`, and an empty string would turn a relative
/// name into an absolute one. `/path/some_entity` yields `/path` whatever
/// `some_entity` is, `/path/only/` yields `/path/only`, and `/a.txt` yields the
/// empty string. Nothing is checked on the filesystem.
pub fn path(file: &str) -> &str {
    file.rfind(['/', '\\']).map_or(".", |i| &file[..i])
}
/// The basename of `file` with any known file extension removed.
///
/// Source `File::stemName`, which is
/// `FileNameUtils::stripExtension(File::basename(file))`. Compound OpenMS
/// extensions are recognised, so `/path/sample.mzML.gz` yields `sample` and
/// `/path/data.featureXML` yields `data`; an unknown extension is stripped at
/// the last dot, `/path/file.txt` yielding `file`; and a dot in a *directory*
/// name is not mistaken for one, `/my.dir/file` yielding `file`. A name that is
/// nothing but an extension, `.mzML`, yields the empty string.
///
/// The extension table is shared with [`crate::format::file_types`]; see
/// `docs/FILE_HANDLING_SUPPORT.md`.
pub fn stem_name(file: &str) -> &str {
    strip_extension(basename(file))
}
/// The extension of `file`, including the leading dot.
///
/// Source `File::extension`, defined as whatever [`basename`] has beyond
/// [`stem_name`]: `/path/sample.mzML.gz` yields `.mzML.gz`, `/path/file.txt`
/// yields `.txt`, a name with no extension yields the empty string, and
/// `.mzML` — which is all extension — yields `.mzML`.
///
/// The cut is taken at the stem's length in bytes. That is always a character
/// boundary, because the stem is a prefix of the basename; a length that is
/// nonetheless out of range yields the empty string rather than panicking, so
/// this can never become the byte-slicing failure that a multi-byte filename
/// triggers in hand-written offset arithmetic.
pub fn extension(file: &str) -> &str {
    let b = basename(file);
    let s = strip_extension(b);
    b.get(s.len()..).unwrap_or("")
}
/// Resolve `file` against the current working directory.
///
/// Source `File::absolutePath` returns the current working directory for an
/// empty input and otherwise `std::filesystem::absolute`, which is purely
/// lexical: `.` and `..` are not collapsed and symbolic links are not resolved,
/// so the result may still contain `..`. This port does the same. Use
/// [`FileContext::find`] when the answer has to exist.
///
/// # Errors
///
/// Returns [`Error::Io`] when the current working directory cannot be read —
/// it was removed, or the process may not read it — and
/// [`Error::InvalidValue`] when the input or the joined result exceeds
/// [`MAX_PATH_BYTES`]. The source can report neither: `fs::current_path()`
/// throws and `fs::absolute` is unchecked.
pub fn absolute_path(file: impl AsRef<Path>) -> Result<PathBuf> {
    let p = file.as_ref();
    check_path(p)?;
    if p.as_os_str().is_empty() {
        // Source absolutePath("") is fs::current_path(), not cwd joined with "".
        return Ok(std::env::current_dir()?);
    }
    if p.is_absolute() {
        Ok(p.to_owned())
    } else {
        let result = std::env::current_dir()?.join(p);
        check_path(&result)?;
        Ok(result)
    }
}
/// The directory holding the running executable.
///
/// Source `File::getExecutablePath` reads `/proc/self/exe`,
/// `_NSGetExecutablePath` or `GetModuleFileNameW` depending on platform, caches
/// the answer in a static, and returns it with a trailing `/` so that
/// `getExecutablePath() + "mytool"` composes — or an empty string when the
/// system call fails. This returns a [`PathBuf`] with no
/// trailing separator; compose with [`Path::join`](std::path::Path::join).
///
/// # Errors
///
/// Returns [`Error::Io`] when the platform cannot report the executable path,
/// where the source returns an empty string a caller is free to ignore.
pub fn get_executable_path() -> Result<PathBuf> {
    let exe = std::env::current_exe()?;
    Ok(parent(&exe).to_owned())
}
/// Create a directory and every missing parent.
///
/// Source `File::makeDir` is `std::filesystem::create_directories` and treats
/// an already existing directory as success, as its own documentation says; so
/// does this. The path may be absolute or relative to the current directory.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the path exceeds [`MAX_PATH_BYTES`] or
/// has more than [`MAX_DEPTH`] components, and [`Error::Io`] when creation
/// fails — including when an existing non-directory occupies the path. The
/// depth ceiling is checked before anything is created, so a rejected call
/// leaves no partial chain behind. The source returns `false` for every failure
/// without distinguishing them.
pub fn make_dir(directory: impl AsRef<Path>) -> Result<()> {
    let p = directory.as_ref();
    check_path(p)?;
    if p.components().count() > MAX_DEPTH {
        return Err(invalid("directory depth limit exceeded"));
    }
    fs::create_dir_all(p)?;
    Ok(())
}
/// Remove a single file, symbolic link, or empty directory.
///
/// Source `File::remove` calls `std::remove`, which on POSIX unlinks a file and
/// also removes an *empty* directory; that is preserved. A path that is already
/// absent is success, as the source documents. A symbolic link is removed
/// itself and never followed, so a dangling link really does disappear.
///
/// # Errors
///
/// Returns [`Error::Io`] when the entry exists but cannot be removed — a
/// nonempty directory, or a parent the process may not write — and
/// [`Error::InvalidValue`] for a path over [`MAX_PATH_BYTES`]. The source
/// returns `false`. Use [`remove_dir_recursively`] for a tree.
pub fn remove(file: impl AsRef<Path>) -> Result<()> {
    let p = file.as_ref();
    check_path(p)?;
    match fs::symlink_metadata(p) {
        Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(()),
        Err(e) => Err(e.into()),
        Ok(m) if m.is_dir() => {
            fs::remove_dir(p)?;
            Ok(())
        }
        Ok(_) => {
            fs::remove_file(p)?;
            Ok(())
        }
    }
}
fn sorted_entries(directory: &Path, work: &mut Work) -> Result<Vec<PathBuf>> {
    work.path(directory)?;
    let mut entries = Vec::new();
    let mut longest = 0;
    for entry in fs::read_dir(directory)? {
        work.entry()?;
        let p = entry?.path();
        longest = longest.max(p.as_os_str().len());
        entries.push(work.owned(&p)?);
    }
    work.consume(mul(
        mul(
            entries.len(),
            usize::BITS as usize - entries.len().leading_zeros() as usize,
        )?,
        add(longest, 1)?,
    )?)?;
    entries.sort();
    Ok(entries)
}
/// The absolute paths of the immediate subdirectories of `directory`, sorted.
///
/// Source `File::listDirectories` is non-recursive, returns absolute paths and
/// sorts them. Its documentation promises "an empty list on any error or if the
/// path is not a directory (no throw)"; this port returns the error instead,
/// because an unreadable directory and an empty one are different answers and a
/// caller given an empty list cannot tell which it got.
///
/// Sorting is by the platform's path bytes, not by a locale collation.
///
/// # Errors
///
/// Returns [`Error::Io`] when `directory` is not a directory or cannot be
/// enumerated, and [`Error::InvalidValue`] when the listing exceeds
/// [`MAX_ENTRIES`], [`MAX_WORK`] or [`MAX_BYTES`].
pub fn list_directories(directory: impl AsRef<Path>) -> Result<Vec<PathBuf>> {
    let mut work = Work::default();
    let mut result = Vec::new();
    for p in sorted_entries(directory.as_ref(), &mut work)? {
        if fs::metadata(&p)?.is_dir() {
            result.push(work.owned(&absolute_path(&p)?)?);
        }
    }
    Ok(result)
}

// Portable Unicode wildcard grammar: *, ?, bracket classes/ranges/negation,
// and backslash escapes. Dot files are ordinary names, as in POSIX fnmatch(0).
#[derive(Debug)]
enum Token {
    Star,
    Any,
    Literal(char),
    Class(bool, Vec<(char, char)>),
}
fn pattern_tokens(pattern: &str, work: &mut Work) -> Result<Vec<Token>> {
    if pattern.len() > MAX_PATH_BYTES {
        return Err(invalid("wildcard pattern limit exceeded"));
    }
    work.copy(mul(pattern.len(), 40)?)?;
    let chars: Vec<char> = pattern.chars().collect();
    let mut result = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        work.consume(1)?;
        let c = chars[i];
        i += 1;
        match c {
            '*' => {
                if !matches!(result.last(), Some(Token::Star)) {
                    result.push(Token::Star);
                }
            }
            '?' => result.push(Token::Any),
            '\\' if i < chars.len() => {
                result.push(Token::Literal(chars[i]));
                i += 1;
            }
            '[' => {
                let start = i;
                let negative = chars.get(i).is_some_and(|c| *c == '!' || *c == '^');
                if negative {
                    i += 1;
                }
                let mut ranges = Vec::new();
                if chars.get(i) == Some(&']') {
                    ranges.push((']', ']'));
                    i += 1;
                }
                while i < chars.len() && chars[i] != ']' {
                    let mut lo = chars[i];
                    i += 1;
                    if lo == '\\' && i < chars.len() {
                        lo = chars[i];
                        i += 1;
                    }
                    let mut hi = lo;
                    if chars.get(i) == Some(&'-') && chars.get(i + 1).is_some_and(|c| *c != ']') {
                        i += 1;
                        hi = chars[i];
                        i += 1;
                        if hi == '\\' && i < chars.len() {
                            hi = chars[i];
                            i += 1;
                        }
                    }
                    work.copy(8)?;
                    ranges.push((lo, hi));
                }
                if i < chars.len() && !ranges.is_empty() {
                    i += 1;
                    result.push(Token::Class(negative, ranges));
                } else {
                    i = start;
                    result.push(Token::Literal('['));
                }
            }
            _ => result.push(Token::Literal(c)),
        }
    }
    Ok(result)
}
fn matches_pattern(tokens: &[Token], value: &str, work: &mut Work) -> Result<bool> {
    work.copy(mul(value.len(), 4)?)?;
    let chars: Vec<char> = value.chars().collect();
    let (mut p, mut s, mut star) = (0, 0, None);
    while s < chars.len() {
        work.consume(1)?;
        let matched = match tokens.get(p) {
            Some(Token::Any) => true,
            Some(Token::Literal(c)) => *c == chars[s],
            Some(Token::Class(negative, ranges)) => {
                work.consume(ranges.len())?;
                ranges.iter().any(|(a, b)| *a <= chars[s] && chars[s] <= *b) != *negative
            }
            Some(Token::Star) => {
                star = Some((p, s));
                p += 1;
                continue;
            }
            None => false,
        };
        if matched {
            p += 1;
            s += 1;
        } else if let Some((position, consumed)) = star {
            s = consumed + 1;
            p = position + 1;
            star = Some((position, s));
        } else {
            return Ok(false);
        }
    }
    while matches!(tokens.get(p), Some(Token::Star)) {
        work.consume(1)?;
        p += 1;
    }
    Ok(p == tokens.len())
}
/// The regular files of `directory` whose names match `pattern`, sorted.
///
/// Source `File::fileList` fills a `StringList` out-parameter and returns
/// "there are matching files"; this returns the list, and an empty list makes
/// the same statement. `full_path` selects whole paths over bare filenames, as
/// in the source. Only regular files are considered, so a subdirectory never
/// matches whatever the pattern says.
///
/// The pattern grammar is a portable, deterministic, case-sensitive subset:
/// `*`, `?`, bracket ranges and classes with `!` or `^` negation, and backslash
/// escapes. A leading dot is an ordinary character, exactly as POSIX `fnmatch`
/// with a zero flag word treats it. Three source behaviours are *not* emulated
/// and are listed in `docs/SYSTEM_FILE_SUPPORT.md`: POSIX locale character
/// classes and collation, the Windows `PathMatchSpecA` rules (pattern lists,
/// case insensitivity), and byte-wise matching — here `?` matches one
/// character, where `fnmatch` matches one byte and so fails to match a
/// multi-byte character.
///
/// # Errors
///
/// Returns [`Error::Io`] when `directory` is not a directory or cannot be
/// enumerated, where the source returns `false`. Returns
/// [`Error::InvalidValue`] when a directory entry's name is not valid UTF-8 —
/// the message names it — and when the pattern or the listing exceeds a bound.
/// A non-UTF-8 name is refused rather than skipped or matched approximately,
/// because both of those answer a question the caller did not ask.
pub fn file_list(
    directory: impl AsRef<Path>,
    pattern: &str,
    full_path: bool,
) -> Result<Vec<PathBuf>> {
    let mut work = Work::default();
    let tokens = pattern_tokens(pattern, &mut work)?;
    let mut result = Vec::new();
    for p in sorted_entries(directory.as_ref(), &mut work)? {
        if !fs::metadata(&p)?.is_file() {
            continue;
        }
        let name = p
            .file_name()
            .ok_or_else(|| invalid("missing directory entry name"))?;
        let text = name.to_str().ok_or_else(|| {
            invalid(&format!(
                "wildcard matching requires UTF-8 filenames; '{}' is not UTF-8",
                p.display()
            ))
        })?;
        if matches_pattern(&tokens, text, &mut work)? {
            result.push(work.owned(if full_path { &p } else { Path::new(name) })?);
        }
    }
    Ok(result)
}

fn same_path(from: &Path, to: &Path) -> bool {
    match (fs::canonicalize(from), fs::canonicalize(to)) {
        (Ok(a), Ok(b)) => a == b,
        _ => false,
    }
}
fn regular(file: &Path) -> Result<()> {
    if !fs::symlink_metadata(file)?.is_file() {
        return Err(invalid(
            "operation requires a regular file, not a directory or symlink",
        ));
    }
    Ok(())
}
/// Copy complete bytes to a sibling temporary, then publish atomically.
fn copy_atomic(from: &Path, to: &Path, overwrite: bool) -> Result<()> {
    checked_paths(from, to)?;
    let mut input = File::open(from)?;
    let metadata = input.metadata()?;
    if !metadata.is_file() {
        return Err(invalid("copy requires a regular file"));
    }
    let (temporary, mut output) = TempFile::create_in(parent(to), Some(to))?;
    io::copy(&mut input, &mut output)?;
    output.set_permissions(metadata.permissions())?;
    output.sync_all()?;
    drop(output);
    if overwrite {
        fs::rename(temporary.path(), to)?;
    } else {
        fs::hard_link(temporary.path(), to)?;
    }
    Ok(())
}
/// Copy a regular file to a destination that does not yet exist.
///
/// Source `File::copy` is `std::filesystem::copy_file` with no options, which
/// fails when the destination exists; that refusal is kept. The bytes are
/// written to a sibling temporary file, given the source file's permissions,
/// flushed to the filesystem, and only then published under `to` with an atomic
/// no-clobber link — so a concurrent reader of `to` never sees a partial file,
/// and a concurrent creator of `to` is not clobbered.
///
/// # Errors
///
/// Returns [`Error::Io`] when `from` cannot be read, when `to` already exists,
/// or when the filesystem cannot publish a hard link, and
/// [`Error::InvalidValue`] when `from` is not a regular file or a path exceeds
/// [`MAX_PATH_BYTES`]. On any failure the destination is left as it was and the
/// temporary is removed.
pub fn copy(from: impl AsRef<Path>, to: impl AsRef<Path>) -> Result<()> {
    copy_atomic(from.as_ref(), to.as_ref(), false)
}
/// Move `from` to `to`.
///
/// Source `File::rename(from, to, overwrite_existing, verbose)`. As in the
/// source, a `from` and `to` that resolve to the same file — symbolic links
/// included — are success with nothing done, and a cross-device move falls back
/// to copy-then-remove: Qt's `QFile::rename` did that silently, and TOPP tools
/// depend on it to move results out of a temporary directory onto a bind mount
/// in a container. The source's `verbose` argument only chose whether the
/// failure was printed to its error log; this returns the failure, so there is
/// no such argument.
///
/// `overwrite = true` uses the platform rename, which replaces the destination
/// atomically. The source instead *deletes* the destination first and renames
/// afterwards, so a rename that then fails leaves the caller with neither file;
/// this port does not, and a missing or unreadable `from` cannot destroy an
/// existing `to`.
///
/// `overwrite = false` publishes with a hard link and removes the source
/// afterwards, so exactly one of several concurrent movers wins and every
/// loser's input is left intact. The source hands this case to
/// `std::filesystem::rename`, which on POSIX overwrites regardless of the flag.
///
/// # Errors
///
/// Returns [`Error::Io`] when the move fails, and [`Error::InvalidValue`] for a
/// path over [`MAX_PATH_BYTES`] and for the two cases portable safe Rust cannot
/// perform: a no-clobber move of a directory or symbolic link, and a
/// cross-device move of one. Neither has an atomic no-clobber primitive in
/// `std`, and this port will not substitute a racy existence check. A
/// cross-device move whose copy succeeds but whose source removal fails reports
/// the failure even though `to` is complete.
pub fn rename(from: impl AsRef<Path>, to: impl AsRef<Path>, overwrite: bool) -> Result<()> {
    let (from, to) = (from.as_ref(), to.as_ref());
    checked_paths(from, to)?;
    if same_path(from, to) {
        return Ok(());
    }
    if overwrite {
        match fs::rename(from, to) {
            Ok(()) => return Ok(()),
            Err(e) if e.kind() == io::ErrorKind::CrossesDevices => (),
            Err(e) => return Err(e.into()),
        }
        regular(from)?;
        copy_atomic(from, to, true)?;
    } else {
        regular(from)?;
        match fs::hard_link(from, to) {
            Ok(()) => (),
            Err(e) if e.kind() == io::ErrorKind::CrossesDevices => copy_atomic(from, to, false)?,
            Err(e) => return Err(e.into()),
        }
    }
    fs::remove_file(from)?;
    Ok(())
}

struct CopyEntry {
    from: PathBuf,
    to: PathBuf,
    directory: bool,
}
fn canonical_target(path: &Path, work: &mut Work) -> Result<PathBuf> {
    let mut path = absolute_path(path)?;
    let mut missing = Vec::new();
    loop {
        work.path(&path)?;
        match fs::canonicalize(&path) {
            Ok(mut resolved) => {
                for part in missing.iter().rev() {
                    resolved.push(part);
                }
                return Ok(resolved);
            }
            Err(e) if e.kind() == io::ErrorKind::NotFound => {
                if missing.len() >= MAX_DEPTH {
                    return Err(invalid("directory depth limit exceeded"));
                }
                let part = path
                    .file_name()
                    .ok_or_else(|| invalid("unresolvable target path"))?
                    .to_owned();
                work.copy(part.len())?;
                missing.push(part);
                if !path.pop() {
                    return Err(e.into());
                }
            }
            Err(e) => return Err(e.into()),
        }
    }
}
fn copy_plan(
    from: &Path,
    to: &Path,
    ancestors: &mut Vec<PathBuf>,
    plan: &mut Vec<CopyEntry>,
    work: &mut Work,
) -> Result<()> {
    if ancestors.len() >= MAX_DEPTH {
        return Err(invalid("recursive copy depth limit exceeded"));
    }
    work.entry()?;
    work.path(from)?;
    work.path(to)?;
    let metadata = fs::metadata(from)?;
    let directory = metadata.is_dir();
    if !directory && !metadata.is_file() {
        return Err(invalid("recursive copy encountered a nonregular file"));
    }
    if directory {
        let canonical = fs::canonicalize(from)?;
        for ancestor in ancestors.iter() {
            work.consume(add(
                ancestor.as_os_str().len(),
                canonical.as_os_str().len(),
            )?)?;
            if *ancestor == canonical {
                return Err(invalid("recursive copy contains a directory cycle"));
            }
        }
        ancestors.push(work.owned(&canonical)?);
    }
    plan.push(CopyEntry {
        from: work.owned(from)?,
        to: work.owned(to)?,
        directory,
    });
    if directory {
        for child in sorted_entries(from, work)? {
            let target = to.join(
                child
                    .file_name()
                    .ok_or_else(|| invalid("missing filename"))?,
            );
            copy_plan(&child, &target, ancestors, plan, work)?;
        }
        ancestors.pop();
    }
    Ok(())
}
/// Copy a directory tree into `to`, creating `to` if it is missing.
///
/// Source `File::copyDirRecursively`. An existing target directory is added to
/// rather than replaced, and `option` decides what happens to a file already
/// there; see [`CopyOptions`]. Directory symbolic links inside the source are
/// followed for copying, as the source's `entry.is_directory()` follows them,
/// but a *destination* that is a directory symbolic link is refused rather than
/// written through.
///
/// Every conflict visible before the first write is found in a preflight, so
/// [`CopyOptions::Cancel`] changes nothing at all. Copying is still not a
/// filesystem transaction: each file is published completely, but an error
/// partway through leaves the files already copied in place.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `from` is not a directory, when `to` is
/// `from` or a descendant of it, when the source holds a directory-link cycle,
/// when a planned target directory is occupied by a file or a symbolic link, or
/// when a bound is exceeded; and [`Error::Io`] for a copy that fails and, under
/// [`CopyOptions::Cancel`], for a conflicting file. The source rejects only the
/// exact self-copy, and would descend into a target nested inside its own
/// source.
pub fn copy_dir_recursively(
    from: impl AsRef<Path>,
    to: impl AsRef<Path>,
    option: CopyOptions,
) -> Result<()> {
    let (from, to) = (from.as_ref(), to.as_ref());
    let mut work = Work::default();
    work.path(from)?;
    work.path(to)?;
    if !fs::metadata(from)?.is_dir() {
        return Err(invalid("recursive copy source is not a directory"));
    }
    let canonical_from = fs::canonicalize(from)?;
    let target = canonical_target(to, &mut work)?;
    if target.starts_with(&canonical_from) {
        return Err(invalid(
            "cannot copy a directory onto itself or a descendant",
        ));
    }
    let mut plan = Vec::new();
    copy_plan(from, to, &mut Vec::new(), &mut plan, &mut work)?;
    // Reject all known conflicts before any destination change. Concurrent
    // conflicts remain protected by atomic no-clobber publication for Skip/Cancel.
    for entry in &plan {
        if let Ok(m) = fs::symlink_metadata(&entry.to) {
            if entry.directory && (!m.is_dir() || m.is_symlink()) {
                return Err(invalid("copy target directory is a file or symlink"));
            }
            if !entry.directory && option == CopyOptions::Cancel {
                return Err(io::Error::new(
                    io::ErrorKind::AlreadyExists,
                    "copy target already exists",
                )
                .into());
            }
        }
    }
    for entry in plan {
        work.consume(1)?;
        if entry.directory {
            fs::create_dir_all(&entry.to)?;
        } else {
            if option == CopyOptions::Skip && fs::symlink_metadata(&entry.to).is_ok() {
                continue;
            }
            match copy_atomic(&entry.from, &entry.to, option == CopyOptions::Overwrite) {
                Err(Error::Io(e))
                    if option == CopyOptions::Skip && e.kind() == io::ErrorKind::AlreadyExists => {}
                result => result?,
            }
        }
    }
    Ok(())
}
fn removal_plan(
    path: &Path,
    depth: usize,
    plan: &mut Vec<(PathBuf, bool)>,
    work: &mut Work,
) -> Result<()> {
    if depth > MAX_DEPTH {
        return Err(invalid("recursive removal depth limit exceeded"));
    }
    work.entry()?;
    let directory = fs::symlink_metadata(path)?.is_dir();
    if directory {
        for child in sorted_entries(path, work)? {
            removal_plan(&child, depth + 1, plan, work)?;
        }
    }
    plan.push((work.owned(path)?, directory));
    Ok(())
}
/// Remove a directory and everything inside it.
///
/// Source `File::removeDirRecursively` is `std::filesystem::remove_all`. A path
/// that is already absent is success. The traversal is planned first, and
/// symbolic links are never followed: a link inside the tree is unlinked, and
/// whatever it pointed at is left alone.
///
/// Removal is not atomic — an error partway through leaves the entries already
/// removed gone. Nothing protects a path-based traversal against another
/// process replacing a directory with a symbolic link between the plan and the
/// removal; that is outside what `std`'s path-based operations can guarantee,
/// for this port and for the source alike.
///
/// # Errors
///
/// Returns [`Error::Io`] when an entry cannot be removed and
/// [`Error::InvalidValue`] when the tree exceeds [`MAX_ENTRIES`], [`MAX_DEPTH`],
/// [`MAX_WORK`] or [`MAX_BYTES`]. The source prints the failure to `stderr` and
/// returns `false`.
pub fn remove_dir_recursively(directory: impl AsRef<Path>) -> Result<()> {
    let p = directory.as_ref();
    check_path(p)?;
    match fs::symlink_metadata(p) {
        Err(e) if e.kind() == io::ErrorKind::NotFound => return Ok(()),
        Err(e) => return Err(e.into()),
        _ => (),
    }
    let mut plan = Vec::new();
    removal_plan(p, 0, &mut plan, &mut Work::default())?;
    for (p, directory) in plan {
        if directory {
            fs::remove_dir(p)?;
        } else {
            fs::remove_file(p)?;
        }
    }
    Ok(())
}
/// Remove a directory and everything inside it.
///
/// Source `File::removeDir` is `std::filesystem::remove_all` as well, despite
/// its name and its own "and all its contents" comment reading as though it
/// were not recursive. It is a synonym for [`remove_dir_recursively`] and is
/// kept only so a reader of the source finds it here. Use [`remove`] for an
/// empty directory.
///
/// # Errors
///
/// As [`remove_dir_recursively`].
pub fn remove_dir(directory: impl AsRef<Path>) -> Result<()> {
    remove_dir_recursively(directory)
}

/// Compare two lists of filenames a tool received as paired inputs.
///
/// Source `File::validateMatchingFileNames`. Passing several input file lists
/// is error-prone because users supply them in different orders, so this
/// reports which of three cases holds; see [`MatchingFileListsStatus`].
/// `basename_only` compares only [`basename`]s, and `ignore_extension` strips
/// the extension first, which is how a list of spectra files is compared
/// against the list of identification files derived from it.
///
/// The source builds a `std::set` from each list and compares the sets, so
/// multiplicity is discarded: `["a", "a", "b"]` against `["a", "b", "b"]` has
/// equal sets and equal length and is reported as
/// [`MatchingFileListsStatus::OrderMismatch`], not a set mismatch. That is
/// preserved.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a name exceeds [`MAX_PATH_BYTES`], or
/// when the lists exceed [`MAX_ENTRIES`] or [`MAX_WORK`]. The source is
/// unbounded.
pub fn validate_matching_file_names(
    a: &[String],
    b: &[String],
    basename_only: bool,
    ignore_extension: bool,
) -> Result<MatchingFileListsStatus> {
    if a.len() != b.len() {
        return Ok(MatchingFileListsStatus::SetMismatch);
    }
    let mut work = Work::default();
    let mut left = Vec::new();
    let mut right = Vec::new();
    let mut different = false;
    let mut longest = 0;
    for (a, b) in a.iter().zip(b) {
        work.entry()?;
        work.copy(add(add(a.len(), b.len())?, 48)?)?;
        if a.len().max(b.len()) > MAX_PATH_BYTES {
            return Err(invalid("filename limit exceeded"));
        }
        let (mut x, mut y) = (a.as_str(), b.as_str());
        if basename_only {
            x = basename(x);
            y = basename(y);
        }
        if ignore_extension {
            x = strip_extension(x);
            y = strip_extension(y);
        }
        longest = longest.max(x.len()).max(y.len());
        different |= x != y;
        left.push(x);
        right.push(y);
    }
    work.consume(mul(
        mul(
            a.len(),
            usize::BITS as usize - a.len().leading_zeros() as usize,
        )?,
        mul(add(longest, 1)?, 2)?,
    )?)?;
    left.sort_unstable();
    left.dedup();
    right.sort_unstable();
    right.dedup();
    Ok(if left != right {
        MatchingFileListsStatus::SetMismatch
    } else if different {
        MatchingFileListsStatus::OrderMismatch
    } else {
        MatchingFileListsStatus::Match
    })
}

/// A temporary file this value owns and removes when it is dropped.
///
/// Source `File::getTemporaryFile` returns a *name* and does not create the
/// file; the name goes into a process-wide list that a static destructor walks
/// at exit. That design cannot say when the file stopped being needed, and it
/// leaves the name free for another process to take between the call and the
/// first write. This guard creates the file immediately and exclusively —
/// mode `0600` on Unix — and removes it on `Drop`.
///
/// A guard built by [`TempFile::alternative`] owns nothing: it neither creates
/// nor removes its path. That is how the source's "if `alternative_file` is not
/// empty, return it and do not destroy it" branch is expressed here.
///
/// Keep the guard alive for as long as anything — a subprocess, a reader — uses
/// the path. `Drop` does not run when the process ends without unwinding, which
/// is the one case the source's exit-time registry covered and this does not.
#[derive(Debug)]
pub struct TempFile {
    path: PathBuf,
    cleanup: bool,
}
impl TempFile {
    /// Create a temporary file in the system temporary directory.
    ///
    /// This is [`std::env::temp_dir`], *not* the source's
    /// `File::getTempDirectory()` order — use
    /// [`FileContext::get_temporary_file`] when `OPENMS_TMPDIR` and the
    /// `temp_dir` setting in `OpenMS.ini` should be honoured.
    ///
    /// # Errors
    ///
    /// As [`TempFile::new_in`].
    pub fn new() -> Result<Self> {
        Self::new_in(std::env::temp_dir())
    }
    /// Create a temporary file in `base`.
    ///
    /// The name is normally `.openms-<unique name>.tmp`; see
    /// [`get_unique_name`]. It is not guaranteed to be: when a candidate cannot
    /// be created, two progressively shorter ones are tried in turn — `.o` with
    /// the decimal probe counter, then that counter alone — so that a directory
    /// close enough to the platform's path limit to reject the descriptive name
    /// still yields a temporary rather than a failure. [`writable`] is the
    /// caller that needs the short rungs, and every other caller will see the
    /// long name; no caller should depend on either spelling.
    ///
    /// Creation is exclusive, so whichever candidate is used is reserved rather
    /// than merely unlikely. A candidate that is already taken moves on to the
    /// next, and when all three are taken the whole ladder is retried with a
    /// fresh unique name.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when `base` does not exist or cannot be written,
    /// reported from the shortest candidate that failed — no shorter name
    /// exists, so the cause was not the length. Returns
    /// [`Error::InvalidValue`] when a candidate path exceeds
    /// [`MAX_PATH_BYTES`] or 100 successive ladders were all taken.
    pub fn new_in(base: impl AsRef<Path>) -> Result<Self> {
        Self::create_in(base.as_ref(), None).map(|(guard, _)| guard)
    }
    /// Wrap a caller-chosen path in a guard that owns nothing.
    ///
    /// The path is not created, not truncated and not removed, and dropping the
    /// guard does nothing. This is the source's optional-output-file idiom: a
    /// tool that needs a file's contents whether or not the user asked for the
    /// file passes the user's path when there is one and gets a self-deleting
    /// temporary when there is not.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the path exceeds
    /// [`MAX_PATH_BYTES`] or holds an interior NUL.
    pub fn alternative(path: impl AsRef<Path>) -> Result<Self> {
        check_path(path.as_ref())?;
        Ok(Self {
            path: path.as_ref().to_owned(),
            cleanup: false,
        })
    }
    /// The three candidate names one attempt tries, longest first.
    ///
    /// The trailing `_`-separated field of a unique name is the process-wide
    /// counter, so the two short rungs stay distinct within this process
    /// without carrying the date, the time or the pid.
    fn candidate_names(unique: &str) -> [String; 3] {
        // filter() rather than a bare unwrap_or: an empty tail would name the
        // base directory itself rather than a file inside it.
        let short = unique
            .rsplit('_')
            .next()
            .filter(|s| !s.is_empty())
            .unwrap_or("0");
        [
            format!(".openms-{unique}.tmp"),
            format!(".o{short}"),
            short.to_owned(),
        ]
    }
    fn create_in(base: &Path, avoid: Option<&Path>) -> Result<(Self, File)> {
        check_path(base)?;
        for _ in 0..100 {
            // Each candidate is shorter than the last. The short ones exist for
            // a directory so close to the platform's path limit that the
            // descriptive name no longer fits although the caller's own, shorter
            // name still would: answering "not writable" there is exactly the
            // false negative writable() exists to avoid. The source reaches the
            // same answer by asking the directory itself with access(2), which
            // safe Rust cannot call.
            //
            // The ladder is walked on *any* refusal rather than on a decoded
            // ENAMETOOLONG, because io::ErrorKind::InvalidFilename is unstable
            // at this crate's minimum Rust version. That costs at most two extra
            // opens, and the error reported is the shortest candidate's: no
            // shorter name exists, so the cause was not the length. A candidate
            // that merely already exists, or that is the very name the caller
            // asked about, is skipped the same way; only when all three are
            // unavailable does the outer loop draw a fresh unique name.
            let unique = get_unique_name(false)?;
            let mut refusal = None;
            for name in Self::candidate_names(&unique) {
                match Self::create_candidate(base, &name, avoid) {
                    Ok(Some(created)) => return Ok(created),
                    Ok(None) => (),
                    Err(e) => refusal = Some(e),
                }
            }
            if let Some(e) = refusal {
                return Err(e);
            }
        }
        Err(invalid("temporary filename collision limit exceeded"))
    }
    fn create_candidate(
        base: &Path,
        name: &str,
        avoid: Option<&Path>,
    ) -> Result<Option<(Self, File)>> {
        if avoid
            .and_then(Path::file_name)
            .is_some_and(|s| s.to_string_lossy().eq_ignore_ascii_case(name))
        {
            return Ok(None);
        }
        let path = base.join(name);
        check_path(&path)?;
        let mut options = OpenOptions::new();
        options.write(true).create_new(true);
        #[cfg(unix)]
        {
            use std::os::unix::fs::OpenOptionsExt;
            options.mode(0o600);
        }
        match options.open(&path) {
            Ok(file) => Ok(Some((
                Self {
                    path,
                    cleanup: true,
                },
                file,
            ))),
            Err(e) if e.kind() == io::ErrorKind::AlreadyExists => Ok(None),
            Err(e) => Err(e.into()),
        }
    }
    /// The path this guard owns, borrowed for as long as the guard lives.
    pub fn path(&self) -> &Path {
        &self.path
    }
    /// Give up ownership and return the path, leaving the file on disk.
    ///
    /// The guard is consumed, so it is the compiler rather than a flag that
    /// prevents a later cleanup. The source has no equivalent for files.
    pub fn keep(mut self) -> PathBuf {
        self.cleanup = false;
        std::mem::take(&mut self.path)
    }
    /// Remove the file now and report whether that worked.
    ///
    /// Use this wherever a failed cleanup should be visible: `Drop` cannot
    /// report one. The guard is consumed either way, and a removal that fails
    /// here is attempted once more on drop.
    ///
    /// # Errors
    ///
    /// As [`remove`]. A guard from [`TempFile::alternative`] owns nothing and
    /// always succeeds without touching the path.
    pub fn close(mut self) -> Result<()> {
        if self.cleanup {
            remove(&self.path)?;
            self.cleanup = false;
        }
        Ok(())
    }
}
impl Drop for TempFile {
    fn drop(&mut self) {
        if self.cleanup {
            let _ = fs::remove_file(&self.path);
        }
    }
}

/// A temporary directory this value owns and removes, with its contents, on drop.
///
/// Source `File::TempDir`, whose destructor calls `removeDirRecursively` unless
/// it was constructed with `keep_dir = true`. The source deletes its copy and
/// move operations; in Rust the type simply has no [`Clone`], and a move
/// transfers the single owner.
///
/// `Drop` cannot report a failure, so a removal that fails there is silent. The
/// source is silent too — `removeDirRecursively` prints to `stderr` and returns
/// a `bool` its destructor discards, and no exception escapes — but here that
/// is a decision rather than an oversight: call [`TempDir::close`] instead
/// wherever the failure matters. `Drop` does not run at all when the process
/// ends without unwinding.
#[derive(Debug)]
pub struct TempDir {
    path: PathBuf,
    cleanup: bool,
}
impl TempDir {
    /// Create a temporary directory under the resolved OpenMS temporary directory.
    ///
    /// `keep_dir` is the source's argument of the same name: `true` leaves the
    /// tree in place when the guard is dropped. The base is
    /// [`FileContext::from_environment`] followed by
    /// [`FileContext::get_temp_directory`], which is the source's
    /// `File::getTempDirectory()` order — `OPENMS_TMPDIR`, then `temp_dir` in
    /// `OpenMS.ini`, then the system temporary directory.
    ///
    /// # Errors
    ///
    /// Whatever [`FileContext::from_environment`] and
    /// [`FileContext::get_temp_directory`] return, plus the errors of
    /// [`TempDir::new_in`].
    pub fn new(keep_dir: bool) -> Result<Self> {
        let context = FileContext::from_environment()?;
        Self::new_in(context.get_temp_directory()?, keep_dir)
    }
    /// Create a temporary directory under `base`, creating `base` if it is missing.
    ///
    /// Source `File::TempDir(base_dir, keep_dir)` builds the name
    /// `<base>/OpenMSTempDir_<unique name>_XXXXXX` and hands it to `mkdtemp`;
    /// this port creates `<base>/OpenMSTempDir_<unique name>` with an exclusive
    /// `mkdir` — mode `0700` on Unix — and retries on the collision the
    /// `XXXXXX` existed to avoid. The name shape therefore matches the source's
    /// two-argument constructor but not its no-argument one, which omits the
    /// `OpenMSTempDir_` prefix; no caller should depend on either spelling.
    ///
    /// An empty `base` means the current directory, and the parent chain is
    /// created first, as in the source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when `base` or the directory cannot be created,
    /// and [`Error::InvalidValue`] when the path exceeds [`MAX_PATH_BYTES`],
    /// the depth exceeds [`MAX_DEPTH`], or 100 successive names were all taken.
    /// The source throws `Exception::UnableToCreateFile` for the same
    /// situations, including its own 100-attempt cap on Windows.
    pub fn new_in(base: impl AsRef<Path>, keep_dir: bool) -> Result<Self> {
        let base = base.as_ref();
        make_dir(if base.as_os_str().is_empty() {
            Path::new(".")
        } else {
            base
        })?;
        for _ in 0..100 {
            let path = base.join(format!("OpenMSTempDir_{}", get_unique_name(false)?));
            check_path(&path)?;
            let mut builder = fs::DirBuilder::new();
            #[cfg(unix)]
            {
                use std::os::unix::fs::DirBuilderExt;
                builder.mode(0o700);
            }
            match builder.create(&path) {
                Ok(()) => {
                    return Ok(Self {
                        path,
                        cleanup: !keep_dir,
                    });
                }
                Err(e) if e.kind() == io::ErrorKind::AlreadyExists => (),
                Err(e) => return Err(e.into()),
            }
        }
        Err(invalid("temporary directory collision limit exceeded"))
    }
    /// The directory this guard owns, borrowed for as long as the guard lives.
    ///
    /// Source `TempDir::getPath` returns the path with a trailing `/`; this
    /// does not, because [`Path::join`](std::path::Path::join) does not need
    /// one.
    pub fn path(&self) -> &Path {
        &self.path
    }
    /// Give up ownership and return the path, leaving the tree on disk.
    ///
    /// The same outcome as `keep_dir = true`, decided later.
    pub fn keep(mut self) -> PathBuf {
        self.cleanup = false;
        std::mem::take(&mut self.path)
    }
    /// Remove the directory and its contents now, and report whether that worked.
    ///
    /// Use this wherever a failed cleanup should be visible: `Drop` cannot
    /// report one. The guard is consumed either way, and a removal that fails
    /// here is attempted once more on drop.
    ///
    /// # Errors
    ///
    /// As [`remove_dir_recursively`]. A guard constructed with
    /// `keep_dir = true` owns nothing and always succeeds.
    pub fn close(mut self) -> Result<()> {
        if self.cleanup {
            remove_dir_recursively(&self.path)?;
            self.cleanup = false;
        }
        Ok(())
    }
}
impl Drop for TempDir {
    fn drop(&mut self) {
        if self.cleanup {
            let _ = remove_dir_recursively(&self.path);
        }
    }
}

/// A name built from the date, the time, optionally the host, the process id
/// and a counter — `20260913_141530_host_4711_1`.
///
/// Source `File::getUniqueName`. The date and time come from `DateTime::now()`
/// with the `-` and `:` separators removed, the counter is a process-wide
/// atomic whose first value is 1, and `include_hostname` prepends the machine
/// name and an underscore.
///
/// Two differences are deliberate. The time is UTC, computed from the system
/// clock with the civil-calendar algorithm, where the source uses local time —
/// so names made on two machines in different zones sort together. The host
/// name is read from `HOSTNAME` (`COMPUTERNAME` on Windows) rather than from
/// `gethostname`, which safe Rust cannot call without a dependency, and every
/// character outside ASCII alphanumerics, `_`, `-` and `.` is replaced by `_`.
/// That variable is usually unset for a non-interactive process, and the host
/// part is then simply absent — so unlike the source, `include_hostname` need
/// not make the name any longer.
///
/// These names are identifiers, not security tokens: they are predictable. What
/// makes a temporary safe here is exclusive creation, not the name.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the system clock precedes the Unix
/// epoch, when the host name exceeds 255 bytes, or when the counter is
/// exhausted. The source checks none of the three.
pub fn get_unique_name(include_hostname: bool) -> Result<String> {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let n = NEXT
        .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |n| n.checked_add(1))
        .map_err(|_| invalid("unique-name counter exhausted"))?
        + 1;
    let seconds = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|_| invalid("system clock precedes Unix epoch"))?
        .as_secs();
    // Gregorian 400-year eras, with March as month zero; independent of timezones.
    let z = (seconds / 86400) as i64 + 719468;
    let era = z / 146097;
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let mut year = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = mp + if mp < 10 { 3 } else { -9 };
    year += i64::from(month <= 2);
    let hour = seconds / 3600 % 24;
    let minute = seconds / 60 % 60;
    let second = seconds % 60;
    let mut host = String::new();
    if include_hostname {
        if let Some(value) = std::env::var_os(if cfg!(windows) {
            "COMPUTERNAME"
        } else {
            "HOSTNAME"
        }) {
            if value.len() > 255 {
                return Err(invalid("hostname limit exceeded"));
            }
            for c in value.to_string_lossy().chars() {
                host.push(
                    if c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.') {
                        c
                    } else {
                        '_'
                    },
                );
            }
            if !host.is_empty() {
                host.push('_');
            }
        }
    }
    Ok(format!(
        "{year:04}{month:02}{day:02}_{hour:02}{minute:02}{second:02}_{host}{}_{n}",
        std::process::id()
    ))
}

/// Split a `PATH`-shaped string into directories, each ending in `/`.
///
/// Source `File::getPathLocations(path)` splits on `:` — `;` on Windows —
/// replaces `\` with `/` and appends a `/`, so `PATH=/usr/bin:/home/unicorn`
/// becomes `{"/usr/bin/", "/home/unicorn/"}`. The source takes the string as an
/// argument precisely so that it can be tested, since environment variables are
/// effectively read-only from inside a test.
///
/// Two source behaviours are preserved: an entirely empty input yields no
/// entries at all, while an empty *component* — the `::` in `/usr/bin::/bin` —
/// yields `"/"`, the filesystem root. Shell quoting is not interpreted, because
/// the source does not interpret it either.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a component exceeds
/// [`MAX_PATH_BYTES`] or the number of components exceeds [`MAX_ENTRIES`].
pub fn get_path_locations(path: &str) -> Result<Vec<String>> {
    if path.is_empty() {
        return Ok(Vec::new());
    }
    let mut work = Work::default();
    let mut result = Vec::new();
    for item in path.split(if cfg!(windows) { ';' } else { ':' }) {
        work.entry()?;
        if item.len() > MAX_PATH_BYTES {
            return Err(invalid("PATH component limit exceeded"));
        }
        work.copy(add(item.len(), 1)?)?;
        let mut p = item.replace('\\', "/");
        if !p.ends_with('/') {
            p.push('/');
        }
        result.push(p);
    }
    Ok(result)
}
/// [`get_path_locations`] applied to the `PATH` environment variable.
///
/// Source `File::getPathLocations()`, the no-argument overload, which treats an
/// unset `PATH` as the empty string and therefore yields no entries.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `PATH` is not valid UTF-8 — the
/// source's splitting is defined on `std::string`, and a lossy conversion would
/// hand back a directory name that is not the one on disk — and whatever
/// [`get_path_locations`] returns.
pub fn get_path_locations_from_environment() -> Result<Vec<String>> {
    let value = std::env::var_os("PATH").unwrap_or_default();
    get_path_locations(
        value
            .to_str()
            .ok_or_else(|| invalid("PATH must be UTF-8 for source string splitting"))?,
    )
}

/// A shared-data directory together with a description of where it was found.
///
/// Source `File::OpenMSDataPath_`, whose two members feed `getOpenMSDataPath()`
/// and `getOpenMSDataPathSource()`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolvedDataPath {
    /// The shared-data directory itself, with no trailing separator.
    pub path: PathBuf,
    /// Human-readable origin, for diagnostics — `executable-relative SDK data`.
    pub source: String,
}
/// The runtime locations `File` derives from its environment, made explicit.
///
/// The source reads `OPENMS_DATA_PATH`, `OPENMS_TMPDIR`, `OPENMS_HOME_PATH`,
/// `XDG_CONFIG_HOME`, `HOME`/`USERPROFILE`, `PATH` and `%PATHEXT%` from the
/// process environment, probes the loaded library with `dladdr` and the running
/// executable with `/proc/self/exe`, falls back to directories compiled into
/// the library, and caches the outcome in a function-local static. None of that
/// can be varied per call, which is why the source's own class test has to
/// mutate the environment of the whole process to test it.
///
/// This port puts the same inputs in a value. Every field is public, so an
/// installer can insert its own data candidate and a test can point the search
/// path somewhere harmless; [`FileContext::from_environment`] fills the value
/// from the environment once, and nothing here reads the environment again
/// afterwards. Two contexts may disagree, and parallel tests are independent.
///
/// A *successful* data-path resolution is cached, as the source's static is —
/// including across later edits to these fields — until
/// [`FileContext::clear_data_cache`] is called. A failed one is not, so a
/// caller may correct the context and try again; the source's static would keep
/// throwing for the life of the process.
#[derive(Clone, Debug)]
pub struct FileContext {
    /// `OPENMS_DATA_PATH`. Authoritative when set: an invalid value is an error
    /// rather than a fall-through to the candidates.
    pub data_override: Option<PathBuf>,
    /// Shared-data directories to try, in order, when there is no override.
    pub data_candidates: Vec<ResolvedDataPath>,
    /// Extra directories [`FileContext::find_doc`] searches before the
    /// data-relative documentation locations.
    pub documentation_directories: Vec<PathBuf>,
    /// The directory `find_sibling_topp_executable` treats as "next to me".
    pub executable_directory: PathBuf,
    /// The user's home directory, source `getOpenMSHomePath`.
    pub home_directory: PathBuf,
    /// The directory holding `OpenMS.ini`, source `getOpenMSConfigDir`.
    pub config_directory: PathBuf,
    /// The system temporary directory, used when nothing overrides it.
    pub temporary_directory: PathBuf,
    /// `OPENMS_TMPDIR`, which outranks `OpenMS.ini`'s `temp_dir`.
    pub temporary_override: Option<PathBuf>,
    /// `OPENMS_HOME_PATH`, which outranks `OpenMS.ini`'s `home_dir`.
    pub user_override: Option<PathBuf>,
    /// `PATH`, already split by [`get_path_locations`].
    pub search_path: Vec<PathBuf>,
    /// Windows `%PATHEXT%` suffixes; `.exe` and `.bat` where there are none.
    pub executable_extensions: Vec<String>,
    data_cache: OnceLock<ResolvedDataPath>,
}
impl FileContext {
    /// Build an isolated context from three directories, reading no environment.
    ///
    /// The data candidates start out as the executable-relative
    /// `../share/OpenMS` and, on macOS, the app bundle's
    /// `../../../share/OpenMS` — the two the source can derive without a
    /// compiled-in path. Insert library-relative, install and build candidates
    /// into [`FileContext::data_candidates`] in the order they should be tried;
    /// the source has them compiled in and cannot be told about another one.
    ///
    /// The configuration directory follows the source's platform rule; see
    /// [`FileContext::get_openms_config_dir`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a directory exceeds
    /// [`MAX_PATH_BYTES`] or holds an interior NUL.
    pub fn new(
        executable_directory: impl AsRef<Path>,
        home_directory: impl AsRef<Path>,
        temporary_directory: impl AsRef<Path>,
    ) -> Result<Self> {
        let (exe, home, temp) = (
            executable_directory.as_ref(),
            home_directory.as_ref(),
            temporary_directory.as_ref(),
        );
        let mut work = Work::default();
        let mut candidates = vec![ResolvedDataPath {
            path: work.owned(&exe.join("../share/OpenMS"))?,
            source: "executable-relative SDK data".into(),
        }];
        if cfg!(target_os = "macos") {
            candidates.push(ResolvedDataPath {
                path: work.owned(&exe.join("../../../share/OpenMS"))?,
                source: "app bundle data".into(),
            });
        }
        Ok(Self {
            data_override: None,
            data_candidates: candidates,
            documentation_directories: Vec::new(),
            executable_directory: work.owned(exe)?,
            home_directory: work.owned(home)?,
            config_directory: work.owned(&source_join(
                home,
                Path::new(if SOURCE_UNIX {
                    ".config/OpenMS"
                } else {
                    ".OpenMS"
                }),
            )?)?,
            temporary_directory: work.owned(temp)?,
            temporary_override: None,
            user_override: None,
            search_path: Vec::new(),
            executable_extensions: vec![".exe".into(), ".bat".into()],
            data_cache: OnceLock::new(),
        })
    }
    /// Snapshot the environment the source would have read.
    ///
    /// Reads `OPENMS_HOME_PATH` (falling back to `HOME`, or `USERPROFILE` on
    /// Windows, or `.`), `OPENMS_TMPDIR`, `OPENMS_DATA_PATH`, `PATH`,
    /// `XDG_CONFIG_HOME` on the platforms where the source's `__unix__` branch
    /// applies, and `%PATHEXT%` on Windows — where a value that does not even
    /// list `.exe` is treated as broken and replaced by `.exe` and `.bat`,
    /// exactly as `executableExtensions_` does. The executable directory comes
    /// from [`get_executable_path`].
    ///
    /// Everything is read once. Later changes to the process environment are
    /// invisible to this context, which is what lets tests construct contexts
    /// instead of mutating environment variables shared by every thread.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when the executable path cannot be determined, and
    /// [`Error::InvalidValue`] when `PATH` is not valid UTF-8 or a value
    /// exceeds a bound.
    pub fn from_environment() -> Result<Self> {
        let home = std::env::var_os(if cfg!(windows) { "USERPROFILE" } else { "HOME" })
            .unwrap_or_else(|| OsString::from("."));
        let override_home = std::env::var_os("OPENMS_HOME_PATH");
        let home = override_home.as_deref().unwrap_or(&home);
        let mut result = Self::new(
            get_executable_path()?,
            Path::new(home),
            std::env::temp_dir(),
        )?;
        result.user_override = override_home.map(PathBuf::from);
        result.temporary_override = std::env::var_os("OPENMS_TMPDIR").map(PathBuf::from);
        result.data_override = std::env::var_os("OPENMS_DATA_PATH").map(PathBuf::from);
        if SOURCE_UNIX {
            if let Some(xdg) = std::env::var_os("XDG_CONFIG_HOME") {
                result.config_directory = source_join(Path::new(&xdg), Path::new("OpenMS"))?;
            }
        }
        result.search_path = get_path_locations_from_environment()?
            .into_iter()
            .map(PathBuf::from)
            .collect();
        if cfg!(windows) {
            if let Ok(ext) = std::env::var("PATHEXT") {
                let mut work = Work::default();
                let mut extensions = Vec::new();
                for s in ext.split(';') {
                    work.entry()?;
                    work.copy(s.len())?;
                    extensions.push(s.to_owned());
                }
                if extensions.iter().any(|s| s.eq_ignore_ascii_case(".exe")) {
                    result.executable_extensions = extensions;
                }
            }
        }
        result.measure(&mut Work::default())?;
        Ok(result)
    }
    fn measure(&self, work: &mut Work) -> Result<()> {
        for p in [
            &self.executable_directory,
            &self.home_directory,
            &self.config_directory,
            &self.temporary_directory,
        ] {
            work.path(p)?;
        }
        for p in [
            &self.data_override,
            &self.temporary_override,
            &self.user_override,
        ]
        .into_iter()
        .flatten()
        {
            work.path(p)?;
        }
        for p in self
            .documentation_directories
            .iter()
            .chain(&self.search_path)
        {
            work.entry()?;
            work.path(p)?;
        }
        for candidate in &self.data_candidates {
            work.entry()?;
            work.path(&candidate.path)?;
            if candidate.source.len() > MAX_PATH_BYTES {
                return Err(invalid("data-path description limit exceeded"));
            }
            work.consume(candidate.source.len())?;
        }
        for extension in &self.executable_extensions {
            work.entry()?;
            work.consume(extension.len())?;
        }
        Ok(())
    }
    /// Forget a cached data-path resolution, so the next call resolves again.
    ///
    /// The source has no equivalent: its cache is a function-local static that
    /// lives until the process ends.
    pub fn clear_data_cache(&mut self) {
        self.data_cache.take();
    }
    /// Resolve the shared-data directory, with a description of where it came from.
    ///
    /// The order reproduces `File::resolveOpenMSDataPath_`:
    ///
    /// 1. [`FileContext::data_override`] — `OPENMS_DATA_PATH` — when it is set.
    ///    It is authoritative: an invalid override is an error and does *not*
    ///    fall through to the candidates, so a stale variable cannot silently
    ///    select the shared data of a different installed version.
    /// 2. Otherwise each entry of [`FileContext::data_candidates`], in order.
    ///    The defaults are the executable-relative `../share/OpenMS` and, on
    ///    macOS, the bundle's `../../../share/OpenMS`. The source additionally
    ///    probes the directory of the *loaded library* with `dladdr` before the
    ///    executable — so that Python bindings and tools installed under
    ///    another prefix find their own data — and then two directories
    ///    compiled into the library. Safe Rust can do neither, so an installer
    ///    or embedder adds them to the candidate list explicitly.
    ///
    /// A candidate is valid when it contains `CHEMISTRY/unimod.xml`, which is
    /// exactly `isOpenMSDataPath_`. The chemistry tables packaged inside this
    /// crate deliberately do not satisfy it: they are not a source-layout data
    /// tree and must not pass for one.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] — with a message naming `OPENMS_DATA_PATH`, as the
    /// source's diagnostic does — when the override is invalid or no candidate
    /// holds the marker file, and [`Error::InvalidValue`] when the context
    /// exceeds a bound. The source throws `Exception::FileNotFound`.
    pub fn resolve_data_path(&self) -> Result<&ResolvedDataPath> {
        self.resolve_data_path_with_work(&mut Work::default())
    }
    fn resolve_data_path_with_work(&self, work: &mut Work) -> Result<&ResolvedDataPath> {
        if let Some(found) = self.data_cache.get() {
            return Ok(found);
        }
        self.measure(work)?;
        let found = if let Some(p) = &self.data_override {
            work.path(&p.join("CHEMISTRY/unimod.xml"))?;
            if !exists(p.join("CHEMISTRY/unimod.xml")) {
                return Err(absent(format!(
                    "invalid OPENMS_DATA_PATH override: {}",
                    p.display()
                )));
            }
            ResolvedDataPath {
                path: work.owned(p)?,
                source: "OPENMS_DATA_PATH env (explicit override)".into(),
            }
        } else {
            let mut found = None;
            for candidate in &self.data_candidates {
                work.path(&candidate.path.join("CHEMISTRY/unimod.xml"))?;
                if exists(candidate.path.join("CHEMISTRY/unimod.xml")) {
                    work.copy(candidate.source.len())?;
                    found = Some(ResolvedDataPath {
                        path: work.owned(&candidate.path)?,
                        source: candidate.source.clone(),
                    });
                    break;
                }
            }
            found.ok_or_else(|| absent("cannot find OpenMS shared data; configure runtime data candidates or OPENMS_DATA_PATH"))?
        };
        let _ = self.data_cache.set(found);
        self.data_cache
            .get()
            .ok_or_else(|| invalid("data-path cache initialization failed"))
    }
    /// The resolved shared-data directory.
    ///
    /// Source `File::getOpenMSDataPath`.
    ///
    /// # Errors
    ///
    /// As [`FileContext::resolve_data_path`].
    pub fn get_openms_data_path(&self) -> Result<PathBuf> {
        let mut work = Work::default();
        let p = &self.resolve_data_path_with_work(&mut work)?.path;
        work.owned(p)
    }
    /// Where [`FileContext::get_openms_data_path`] resolved from, for diagnostics.
    ///
    /// Source `File::getOpenMSDataPathSource`, whose own example is
    /// `exe-relative (../share/OpenMS)`.
    ///
    /// # Errors
    ///
    /// As [`FileContext::resolve_data_path`].
    pub fn get_openms_data_path_source(&self) -> Result<&str> {
        Ok(&self.resolve_data_path()?.source)
    }
    /// The user's home directory as this context recorded it.
    ///
    /// Source `File::getOpenMSHomePath` reads `OPENMS_HOME_PATH` and falls back
    /// to `HOME` — `USERPROFILE` on Windows — or to `.`;
    /// [`FileContext::from_environment`] does that once and stores the answer,
    /// so this cannot fail.
    pub fn get_openms_home_path(&self) -> &Path {
        &self.home_directory
    }
    /// The per-user configuration directory that holds `OpenMS.ini`.
    ///
    /// Source `File::getOpenMSConfigDir` follows the XDG base directory
    /// specification wherever the compiler defines `__unix__`:
    /// `$XDG_CONFIG_HOME/OpenMS` when that variable is set, otherwise
    /// `<home>/.config/OpenMS`. Everywhere else — which includes macOS, where
    /// Apple's compilers do not define `__unix__` — it is `<home>/.OpenMS`.
    /// This port makes the same distinction rather than using Rust's `unix`
    /// configuration flag, which is set on macOS and would move it onto the XDG
    /// branch.
    ///
    /// The path carries no trailing separator, and asking for it does not
    /// create it.
    pub fn get_openms_config_dir(&self) -> &Path {
        &self.config_directory
    }
    /// The directory this context was told the running executable lives in.
    ///
    /// Source `File::getExecutablePath` discovers it on every call from a
    /// cached static; [`FileContext::from_environment`] discovers it once with
    /// [`get_executable_path`] and
    /// stores it, so an embedder that ships its tools elsewhere — or a test —
    /// can substitute another directory.
    pub fn get_executable_path(&self) -> &Path {
        &self.executable_directory
    }
    /// Locate `filename` among the given directories and the shared-data tree.
    ///
    /// Source `File::find`, in the source's order:
    ///
    /// 1. `filename` itself, when it already exists. The source explains why
    ///    this comes first — `File::find(File::find("CHEMISTRY/unimod.xml"))`
    ///    has to return the same absolute path rather than fail on it.
    /// 2. each entry of `directories`, in order;
    /// 3. the shared-data directory.
    ///
    /// The shared data is resolved *before* any directory is searched, so a
    /// caller whose own directory holds the file still gets an error when the
    /// installation's data tree is missing. That is the source's operation
    /// order, preserved deliberately; it is also what lets the failure message
    /// name the resolved data path.
    ///
    /// The result is lexically normalised, as the source's `lexically_normal`
    /// does, which collapses `..` without consulting the filesystem.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when `filename` is empty or only whitespace — the
    /// source throws `FileNotFound` there, because prefixing a directory to an
    /// empty name would otherwise "find" the directory — and when the file is
    /// in none of the locations, with a message naming the resolved shared-data
    /// directory, where it came from, and `OPENMS_DATA_PATH`, which is the
    /// usual culprit. Returns [`Error::InvalidValue`] when a bound is exceeded.
    pub fn find(&self, filename: impl AsRef<Path>, directories: &[PathBuf]) -> Result<PathBuf> {
        let mut work = Work::default();
        self.find_with_work(filename.as_ref(), directories, &mut work)
    }
    fn find_with_work(
        &self,
        filename: &Path,
        directories: &[PathBuf],
        work: &mut Work,
    ) -> Result<PathBuf> {
        work.path(filename)?;
        if exists(filename) {
            return work.owned(filename);
        }
        if filename
            .as_os_str()
            .as_encoded_bytes()
            .iter()
            .all(u8::is_ascii_whitespace)
        {
            return Err(absent("empty resource filename"));
        }
        self.measure(work)?;
        // Resolving data is intentionally mandatory even when a supplied search
        // directory would contain the requested file (the source operation order).
        let data = self.resolve_data_path_with_work(work)?;
        for base in directories.iter().chain(std::iter::once(&data.path)) {
            work.entry()?;
            work.path(base)?;
            let candidate = source_join(base, filename)?;
            work.path(&candidate)?;
            if exists(&candidate) {
                return work.owned(&lexical_normal(&candidate));
            }
        }
        Err(absent(format!(
            "resource '{}' not found; OpenMS shared data '{}' via {}; check OPENMS_DATA_PATH",
            filename.display(),
            data.path.display(),
            data.source
        )))
    }
    /// Locate a documentation file.
    ///
    /// Source `File::findDoc` searches, through `File::find`, the doc directory
    /// relative to the compiled-in binary and source paths, the one relative to
    /// the shared data (`<data>/../../doc`), two further compiled-in
    /// documentation paths, and on macOS three `Documentation` variants used by
    /// the packages. This port keeps the two that need no compiled-in checkout
    /// — `<data>/../../doc` and, on macOS, `<data>/../../Documentation` — and
    /// takes the rest from [`FileContext::documentation_directories`], which an
    /// installer or a test fills in.
    ///
    /// The source's own note applies unchanged: when this fails, try the web
    /// documentation instead.
    ///
    /// # Errors
    ///
    /// As [`FileContext::find`]. In particular the shared data must resolve
    /// even when the file would have been found in an explicit documentation
    /// directory.
    pub fn find_doc(&self, filename: impl AsRef<Path>) -> Result<PathBuf> {
        let mut work = Work::default();
        self.measure(&mut work)?;
        let mut directories = Vec::new();
        for p in &self.documentation_directories {
            directories.push(work.owned(p)?);
        }
        let data = &self.resolve_data_path_with_work(&mut work)?.path;
        directories.push(work.owned(&data.join("../../doc"))?);
        if cfg!(target_os = "macos") {
            directories.push(work.owned(&data.join("../../Documentation"))?);
        }
        self.find_with_work(filename.as_ref(), &directories, &mut work)
    }
    /// Load `OpenMS.ini`, or the built-in defaults when there is none.
    ///
    /// Source `File::getSystemParameters`. The file is
    /// `<config dir>/OpenMS.ini`; see
    /// [`FileContext::get_openms_config_dir`]. When it is missing or unreadable
    /// the five defaults are returned — `version`, an empty `home_dir`, an
    /// empty `temp_dir`, an empty `id_db_dir` list and `threads = 1` — and
    /// nothing is written, exactly as in the source, whose comment reads "no
    /// file, lets keep it that way".
    ///
    /// When the file exists but its `version` is missing or stale, the source
    /// builds a repaired tree from the defaults, updates that with the file's
    /// contents, and then returns the *file's* tree, discarding the repair. The
    /// only observable effect is that `version` has been overwritten in the
    /// returned tree; missing defaults stay missing. That is preserved rather
    /// than quietly corrected, and no file is ever rewritten. See
    /// [`FileContext::get_system_parameters_with_warnings`] for the diagnostics.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Unsupported`] when a readable `OpenMS.ini` exists but
    /// the `paramxml` feature is off, and whatever the XML load returns for a
    /// malformed file. The source has no such feature split.
    pub fn get_system_parameters(&self) -> Result<Param> {
        self.get_system_parameters_with_warnings().map(|(p, _)| p)
    }
    /// [`FileContext::get_system_parameters`] together with its diagnostics.
    ///
    /// The source writes these to `OPENMS_LOG_WARN` — "Broken file '...'
    /// discovered. The 'version' tag is missing.", "File '...' is deprecated."
    /// and "Updating missing/wrong entries in '...' with defaults!" — where a
    /// library has no business writing to a log its caller does not own. They
    /// are returned instead, in the source's order, followed by any message the
    /// defaults update itself produced.
    ///
    /// # Errors
    ///
    /// As [`FileContext::get_system_parameters`].
    pub fn get_system_parameters_with_warnings(&self) -> Result<(Param, Vec<String>)> {
        self.measure(&mut Work::default())?;
        let filename = self.config_directory.join("OpenMS.ini");
        check_path(&filename)?;
        if !readable(&filename) {
            return Ok((system_parameter_defaults()?, Vec::new()));
        }
        #[cfg(not(feature = "paramxml"))]
        {
            Err(Error::Unsupported(
                "reading existing OpenMS.ini requires the paramxml feature".into(),
            ))
        }
        #[cfg(feature = "paramxml")]
        {
            let mut p = crate::format::paramxml::load(&filename)?;
            let mut warnings = Vec::new();
            let has_version = p.exists("version")?;
            if !has_version || p.value("version")? != &ParamValue::from(CORE_SDK_VERSION) {
                warnings.push(if has_version {
                    format!("File '{}' is deprecated.", filename.display())
                } else {
                    format!(
                        "Broken file '{}' discovered. The 'version' tag is missing.",
                        filename.display()
                    )
                });
                warnings.push(format!(
                    "Updating missing/wrong entries in '{}' with defaults!",
                    filename.display()
                ));
                p.set_value("version", ParamValue::from(CORE_SDK_VERSION), "", &[])?;
                let mut discarded = system_parameter_defaults()?;
                warnings.extend(discarded.update(&p, false)?.messages);
                // Source builds the repaired tree but returns p. Do not silently
                // fix that observable quirk or rewrite the user's configuration.
            }
            Ok((p, warnings))
        }
    }
    /// The directory temporary files belong in.
    ///
    /// Source `File::getTempDirectory` takes the first of `OPENMS_TMPDIR`
    /// ([`FileContext::temporary_override`]), a nonempty `temp_dir` in
    /// `OpenMS.ini`, and the system temporary directory. The configuration is
    /// loaded first in the source too — it is the function's first statement,
    /// before the environment variable is even looked at — so a malformed
    /// `OpenMS.ini` is an error here even when the override would have
    /// answered. That order is preserved.
    ///
    /// A configured value counts only when it has non-whitespace content, and
    /// is then used verbatim, leading and trailing spaces included, since those
    /// may be part of a real directory name.
    ///
    /// # Errors
    ///
    /// As [`FileContext::get_system_parameters`], plus [`Error::InvalidValue`]
    /// for a configured path over [`MAX_PATH_BYTES`].
    pub fn get_temp_directory(&self) -> Result<PathBuf> {
        let p = self.get_system_parameters()?;
        if let Some(path) = &self.temporary_override {
            return Work::default().owned(path);
        }
        if let Some(path) = configured_directory(&p, "temp_dir")? {
            return Ok(path);
        }
        Work::default().owned(&self.temporary_directory)
    }
    /// The directory result files belong in.
    ///
    /// Source `File::getUserDirectory` takes the first of `OPENMS_HOME_PATH`
    /// ([`FileContext::user_override`]), a nonempty `home_dir` in
    /// `OpenMS.ini`, and the user's home directory. The header documents the
    /// variable as `OPENMS_HOME_DIR`; the implementation reads
    /// `OPENMS_HOME_PATH`, and the implementation is what this follows.
    ///
    /// The source finishes with `ensureLastChar(dir, '/')`, which turns an
    /// explicitly empty value into `/`, the filesystem root. That is
    /// reproduced, because a tool configured that way writes to the root under
    /// the source as well. The trailing separator itself is not, since
    /// [`PathBuf`] does not need one.
    ///
    /// # Errors
    ///
    /// As [`FileContext::get_temp_directory`].
    pub fn get_user_directory(&self) -> Result<PathBuf> {
        let p = self.get_system_parameters()?;
        let result = if let Some(path) = &self.user_override {
            Work::default().owned(path)?
        } else if let Some(path) = configured_directory(&p, "home_dir")? {
            path
        } else {
            Work::default().owned(&self.home_directory)?
        };
        // ensureLastChar in the source turns an explicitly empty home into root.
        Ok(if result.as_os_str().is_empty() {
            PathBuf::from("/")
        } else {
            result
        })
    }
    /// Locate a sequence database by name under the configured `id_db_dir`.
    ///
    /// Source `File::findDatabase` hands the `id_db_dir` list from `OpenMS.ini`
    /// to `File::find`, which is what lets a TOPP tool be given a bare database
    /// filename. The source also logs the resolved name at info level and the
    /// failure at error level before rethrowing; this returns both instead of
    /// logging either.
    ///
    /// # Errors
    ///
    /// As [`FileContext::get_system_parameters`] and [`FileContext::find`].
    pub fn find_database(&self, filename: impl AsRef<Path>) -> Result<PathBuf> {
        let p = self.get_system_parameters()?;
        let names = p.value("id_db_dir")?.as_string_list()?;
        let mut work = Work::default();
        let mut paths = Vec::new();
        for name in names {
            work.entry()?;
            paths.push(work.owned(Path::new(name))?);
        }
        self.find_with_work(filename.as_ref(), &paths, &mut work)
    }
    /// Search `PATH` for an executable, as `which` and `where` do.
    ///
    /// Source `File::findExecutable` takes the name by mutable reference,
    /// overwrites it with the full path on success and returns a `bool`,
    /// leaving the name untouched on failure. This returns `Ok(None)` for "not
    /// found", so nothing is overwritten and a miss cannot be confused with an
    /// error.
    ///
    /// A name that already resolves to an existing non-directory is returned
    /// unchanged, without consulting `PATH` at all. Otherwise the entries of
    /// [`FileContext::search_path`] are tried in order and the first hit wins.
    /// On Windows a name with no dot is tried with each suffix in
    /// [`FileContext::executable_extensions`].
    ///
    /// As the source's own note says, this does not require the file to have
    /// execute permission — it is not tested. Use
    /// [`executable`] for the permission bits.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a path or the context exceeds a
    /// bound. A name that is simply not found is `Ok(None)`.
    pub fn find_executable(&self, name: impl AsRef<Path>) -> Result<Option<PathBuf>> {
        let name = name.as_ref();
        let mut work = Work::default();
        self.measure(&mut work)?;
        work.path(name)?;
        if exists(name) && !is_directory(name) {
            return work.owned(name).map(Some);
        }
        let mut names = Vec::new();
        if cfg!(windows) && !name.as_os_str().as_encoded_bytes().contains(&b'.') {
            for ext in &self.executable_extensions {
                work.copy(add(name.as_os_str().len(), ext.len())?)?;
                let mut value = name.as_os_str().to_owned();
                value.push(ext);
                names.push(PathBuf::from(value));
            }
        } else {
            names.push(work.owned(name)?);
        }
        for directory in &self.search_path {
            for name in &names {
                work.entry()?;
                let candidate = source_join(directory, name)?;
                work.path(&candidate)?;
                if exists(&candidate) && !is_directory(&candidate) {
                    return work.owned(&candidate).map(Some);
                }
            }
        }
        Ok(None)
    }
    /// Find a TOPP tool next to the running executable.
    ///
    /// Source `File::findSiblingTOPPExecutable` searches the executable's own
    /// directory and, on macOS, the three bundle-relative locations
    /// `../../../`, `../../../TOPP/` and `../../../bin/`. `PATH` is
    /// deliberately not consulted: the source carries a `TODO` about probing it
    /// and its class test pins the current siblings-only contract by asserting
    /// that a tool reachable only through `PATH` is reported as not found. On
    /// Windows a name without a `.exe` suffix gains one.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when no sibling of that name exists, where the
    /// source throws `Exception::FileNotFound`, and [`Error::InvalidValue`]
    /// when a bound is exceeded.
    pub fn find_sibling_topp_executable(&self, tool_name: &str) -> Result<PathBuf> {
        let mut work = Work::default();
        self.measure(&mut work)?;
        work.copy(tool_name.len())?;
        let mut name = tool_name.to_owned();
        if cfg!(windows) && !name.ends_with(".exe") {
            name.push_str(".exe");
        }
        let candidate = source_join(&self.executable_directory, Path::new(&name))?;
        work.path(&candidate)?;
        if exists(&candidate) {
            return work.owned(&candidate);
        }
        if cfg!(target_os = "macos") {
            for prefix in ["../../../", "../../../TOPP", "../../../bin"] {
                let candidate =
                    source_join(&self.executable_directory.join(prefix), Path::new(&name))?;
                work.path(&candidate)?;
                if exists(&candidate) {
                    return work.owned(&candidate);
                }
            }
        }
        Err(absent(format!(
            "sibling TOPP executable '{tool_name}' not found"
        )))
    }
    /// A temporary file, or the caller's own path when one was supplied.
    ///
    /// Source `File::getTemporaryFile(alternative_file)`: a nonempty
    /// `alternative_file` is returned untouched and is not scheduled for
    /// deletion, which is how a tool handles an optional output file whose
    /// contents it needs whether or not the user asked for the file. `None` —
    /// and, as in the source, an empty path — yields an owned [`TempFile`]
    /// under [`FileContext::get_temp_directory`].
    ///
    /// The source returns a name and creates nothing; this creates the file
    /// exclusively, so the name cannot be taken by another process between the
    /// call and the first write.
    ///
    /// # Errors
    ///
    /// As [`FileContext::get_temp_directory`] and [`TempFile::new_in`].
    pub fn get_temporary_file(&self, alternative: Option<&Path>) -> Result<TempFile> {
        if let Some(path) = alternative.filter(|p| !p.as_os_str().is_empty()) {
            TempFile::alternative(path)
        } else {
            TempFile::new_in(self.get_temp_directory()?)
        }
    }
}
fn configured_directory(parameters: &Param, key: &str) -> Result<Option<PathBuf>> {
    if !parameters.exists(key)? {
        return Ok(None);
    }
    let text = parameters.value(key)?.to_text(true)?;
    if text
        .trim_matches(|c: char| c.is_ascii_whitespace())
        .is_empty()
    {
        Ok(None)
    } else {
        check_path(Path::new(&text))?;
        Ok(Some(PathBuf::from(text)))
    }
}
fn system_parameter_defaults() -> Result<Param> {
    let mut p = Param::new();
    for (key, value) in [
        ("version", ParamValue::from(CORE_SDK_VERSION)),
        ("home_dir", ParamValue::from("")),
        ("temp_dir", ParamValue::from("")),
    ] {
        p.set_value(key, value, "", &[])?;
    }
    p.set_value("id_db_dir",ParamValue::StringList(Vec::new()),"Default directory for FASTA and psq files used as databased for id engines. This allows you to specify just the filename of the DB in the respective TOPP tool, and the database will be searched in the directories specified here ",&[])?;
    p.set_value("threads", ParamValue::from(1), "", &[])?;
    Ok(p)
}
fn source_join(base: &Path, filename: &Path) -> Result<PathBuf> {
    check_path(base)?;
    check_path(filename)?;
    let mut result = base.as_os_str().to_owned();
    if !result.as_encoded_bytes().ends_with(b"/") {
        result.push("/");
    }
    result.push(filename.as_os_str());
    let result = PathBuf::from(result);
    check_path(&result)?;
    Ok(result)
}
fn lexical_normal(path: &Path) -> PathBuf {
    let mut result = PathBuf::new();
    for c in path.components() {
        match c {
            Component::CurDir => (),
            Component::ParentDir if result.file_name().is_some_and(|n| n != "..") => {
                result.pop();
            }
            Component::ParentDir if result.has_root() => (),
            _ => result.push(c.as_os_str()),
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn temporary_probe_excludes_the_requested_basename_even_with_case_or_dot_parent() {
        let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
        for name in [".openms-requested.tmp", ".OPENMS-REQUESTED.TMP"] {
            let avoid = dir.path().join(".").join(name);
            assert!(
                TempFile::create_candidate(dir.path(), ".openms-requested.tmp", Some(&avoid))
                    .unwrap()
                    .is_none()
            );
            assert!(!avoid.exists());
        }
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
    }
    #[test]
    fn probe_candidates_shorten_strictly_and_end_at_the_bare_counter() {
        // writable()'s divergence from access(2) is exactly the length of the
        // shortest rung, so pin the ladder the rustdoc describes.
        let names = TempFile::candidate_names("20260913_141530_4711_7");
        assert_eq!(
            names,
            [
                ".openms-20260913_141530_4711_7.tmp".to_owned(),
                ".o7".to_owned(),
                "7".to_owned()
            ]
        );
        for pair in names.windows(2) {
            assert!(pair[1].len() < pair[0].len());
        }
        // No rung may ever be empty: that would name the directory itself.
        for unique in ["20260913_141530_4711_", "", "_"] {
            assert!(
                TempFile::candidate_names(unique)
                    .iter()
                    .all(|n| !n.is_empty())
            );
        }
        // The real generator keeps the counter in the trailing field.
        let unique = get_unique_name(false).unwrap();
        let last = TempFile::candidate_names(&unique)[2].clone();
        assert!(!last.is_empty() && last.bytes().all(|b| b.is_ascii_digit()));
    }
    #[test]
    fn wildcard_work_is_cumulative_across_names() {
        let mut work = Work::default();
        let tokens = pattern_tokens("*a*b", &mut work).unwrap();
        work.work = 20;
        assert!(matches_pattern(&tokens, "ab", &mut work).unwrap());
        assert!(matches_pattern(&tokens, "aaaaaaaaaaaaaaaaaaaa", &mut work).is_err());
    }
}
