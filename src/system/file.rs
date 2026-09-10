// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Checked native counterparts of the Core SDK SYSTEM/File surface.
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

pub const MAX_ENTRIES: usize = 100_000;
pub const MAX_PATH_BYTES: usize = 1024 * 1024;
pub const MAX_WORK: usize = 50_000_000;
pub const MAX_BYTES: usize = 64 * 1024 * 1024;
pub const MAX_DEPTH: usize = 128;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CopyOptions {
    #[default]
    Overwrite,
    Skip,
    Cancel,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum MatchingFileListsStatus {
    Match = 0,
    OrderMismatch = 1,
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

pub fn exists(file: impl AsRef<Path>) -> bool {
    check_path(file.as_ref()).is_ok() && file.as_ref().try_exists().unwrap_or(false)
}
/// Missing paths, directories, and metadata errors are empty in the source API.
pub fn empty(file: impl AsRef<Path>) -> bool {
    check_path(file.as_ref()).is_err()
        || fs::metadata(file).map_or(true, |m| !m.is_file() || m.len() == 0)
}
pub fn is_directory(file: impl AsRef<Path>) -> bool {
    check_path(file.as_ref()).is_ok() && file.as_ref().is_dir()
}
pub fn file_size(file: impl AsRef<Path>) -> Result<u64> {
    check_path(file.as_ref())?;
    let m = fs::metadata(file)?;
    if !m.is_file() {
        return Err(invalid("file size requires a regular file"));
    }
    Ok(m.len())
}
/// Whole seconds since the Unix epoch, including dates before that epoch.
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
/// An effective-process read/open probe, rather than a real-UID access(2) query.
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
/// Never opens the requested file with truncation, or creates a missing target.
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

pub fn basename(file: &str) -> &str {
    file.rsplit(['/', '\\']).next().unwrap_or("")
}
pub fn path(file: &str) -> &str {
    file.rfind(['/', '\\']).map_or(".", |i| &file[..i])
}
pub fn stem_name(file: &str) -> &str {
    strip_extension(basename(file))
}
pub fn extension(file: &str) -> &str {
    let b = basename(file);
    let s = strip_extension(b);
    b.get(s.len()..).unwrap_or("")
}
pub fn absolute_path(file: impl AsRef<Path>) -> Result<PathBuf> {
    let p = file.as_ref();
    check_path(p)?;
    if p.is_absolute() {
        Ok(p.to_owned())
    } else {
        let result = std::env::current_dir()?.join(p);
        check_path(&result)?;
        Ok(result)
    }
}
pub fn get_executable_path() -> Result<PathBuf> {
    let exe = std::env::current_exe()?;
    Ok(parent(&exe).to_owned())
}
pub fn make_dir(directory: impl AsRef<Path>) -> Result<()> {
    let p = directory.as_ref();
    check_path(p)?;
    if p.components().count() > MAX_DEPTH {
        return Err(invalid("directory depth limit exceeded"));
    }
    fs::create_dir_all(p)?;
    Ok(())
}
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
        let text = name
            .to_str()
            .ok_or_else(|| invalid("wildcard matching requires UTF-8 filenames"))?;
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
/// Source copy refuses an existing destination. Publication is atomic and no-clobber.
pub fn copy(from: impl AsRef<Path>, to: impl AsRef<Path>) -> Result<()> {
    copy_atomic(from.as_ref(), to.as_ref(), false)
}
/// Moves regular files; no-clobber publication uses hard links, not a racy exists check.
/// No-clobber directory/symlink moves need platform-specific primitives and are rejected.
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
/// Source removeDir is recursive, just like removeDirRecursively.
pub fn remove_dir(directory: impl AsRef<Path>) -> Result<()> {
    remove_dir_recursively(directory)
}

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

/// An owned temporary file, created exclusively and removed on drop unless kept.
/// An alternative path is unowned and is never created or removed by this guard.
#[derive(Debug)]
pub struct TempFile {
    path: PathBuf,
    cleanup: bool,
}
impl TempFile {
    pub fn new() -> Result<Self> {
        Self::new_in(std::env::temp_dir())
    }
    pub fn new_in(base: impl AsRef<Path>) -> Result<Self> {
        Self::create_in(base.as_ref(), None).map(|(guard, _)| guard)
    }
    pub fn alternative(path: impl AsRef<Path>) -> Result<Self> {
        check_path(path.as_ref())?;
        Ok(Self {
            path: path.as_ref().to_owned(),
            cleanup: false,
        })
    }
    fn create_in(base: &Path, avoid: Option<&Path>) -> Result<(Self, File)> {
        check_path(base)?;
        for _ in 0..100 {
            let name = format!(".openms-{}.tmp", get_unique_name(false)?);
            if let Some(created) = Self::create_candidate(base, &name, avoid)? {
                return Ok(created);
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
    pub fn path(&self) -> &Path {
        &self.path
    }
    pub fn keep(mut self) -> PathBuf {
        self.cleanup = false;
        std::mem::take(&mut self.path)
    }
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

/// Source TempDir ownership; Rust moves transfer ownership, and Clone is absent.
#[derive(Debug)]
pub struct TempDir {
    path: PathBuf,
    cleanup: bool,
}
impl TempDir {
    pub fn new(keep_dir: bool) -> Result<Self> {
        let context = FileContext::from_environment()?;
        Self::new_in(context.get_temp_directory()?, keep_dir)
    }
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
    pub fn path(&self) -> &Path {
        &self.path
    }
    pub fn keep(mut self) -> PathBuf {
        self.cleanup = false;
        std::mem::take(&mut self.path)
    }
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

/// UTC date/time, optional sanitized environment hostname, process ID, counter.
/// These names are identifiers, not secure random tokens; creation is exclusive.
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

/// PATH parsing follows source literal separators rather than shell quoting.
/// Empty whole input is empty; an empty component becomes `/`, as in the source.
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
pub fn get_path_locations_from_environment() -> Result<Vec<String>> {
    let value = std::env::var_os("PATH").unwrap_or_default();
    get_path_locations(
        value
            .to_str()
            .ok_or_else(|| invalid("PATH must be UTF-8 for source string splitting"))?,
    )
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolvedDataPath {
    pub path: PathBuf,
    pub source: String,
}
/// Explicit runtime locations, with no implicit dependency on a source checkout.
/// Public configuration edits affect future calls; a successful data resolution
/// stays cached until `clear_data_cache`, just as the source caches its first hit.
#[derive(Clone, Debug)]
pub struct FileContext {
    pub data_override: Option<PathBuf>,
    pub data_candidates: Vec<ResolvedDataPath>,
    pub documentation_directories: Vec<PathBuf>,
    pub executable_directory: PathBuf,
    pub home_directory: PathBuf,
    pub config_directory: PathBuf,
    pub temporary_directory: PathBuf,
    pub temporary_override: Option<PathBuf>,
    pub user_override: Option<PathBuf>,
    pub search_path: Vec<PathBuf>,
    pub executable_extensions: Vec<String>,
    data_cache: OnceLock<ResolvedDataPath>,
}
impl FileContext {
    /// Construct an isolated context. Add library/install/build candidates in
    /// desired order; default candidates are relative to the executable only.
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
                Path::new(if cfg!(unix) {
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
    /// Snapshot relevant environment values; tests can instead construct a
    /// context directly, without changing process-global environment variables.
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
        if cfg!(unix) {
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
    pub fn clear_data_cache(&mut self) {
        self.data_cache.take();
    }
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
    pub fn get_openms_data_path(&self) -> Result<PathBuf> {
        let mut work = Work::default();
        let p = &self.resolve_data_path_with_work(&mut work)?.path;
        work.owned(p)
    }
    pub fn get_openms_data_path_source(&self) -> Result<&str> {
        Ok(&self.resolve_data_path()?.source)
    }
    pub fn get_openms_home_path(&self) -> &Path {
        &self.home_directory
    }
    pub fn get_openms_config_dir(&self) -> &Path {
        &self.config_directory
    }
    pub fn get_executable_path(&self) -> &Path {
        &self.executable_directory
    }
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
    pub fn get_system_parameters(&self) -> Result<Param> {
        self.get_system_parameters_with_warnings().map(|(p, _)| p)
    }
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
    fn wildcard_work_is_cumulative_across_names() {
        let mut work = Work::default();
        let tokens = pattern_tokens("*a*b", &mut work).unwrap();
        work.work = 20;
        assert!(matches_pattern(&tokens, "ab", &mut work).unwrap());
        assert!(matches_pattern(&tokens, "aaaaaaaaaaaaaaaaaaaa", &mut work).is_err());
    }
}
