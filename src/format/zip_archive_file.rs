// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Small helpers for working with ZIP archives (`ZipArchiveFile.h`).
//!
//! Zip a directory, unzip an archive into a temporary directory, add or
//! replace an entry from a file on disk, list the entries, write a JSON size
//! index, and extract a single entry. The source uses libzip; this port uses
//! the `zip` crate for the archive structure and decodes entries through
//! [`ZipIfstream`](super::zip_ifstream::ZipIfstream)'s decoder, so the same
//! stored, deflate and bzip2 entries are supported and every extracted entry
//! is checked against its CRC-32 and declared size.
//!
//! The OpenSWATH Parquet layouts (`.oswpq`) are the source's consumers: they
//! are either a directory or a store-only archive of one, which is why every
//! reading function here also accepts an unpacked directory.
//!
//! # Differences from the source
//!
//! - Archives are written to a temporary sibling file and renamed over the
//!   destination only when complete, as libzip's `zip_close` does; the
//!   destination is never left half-written.
//! - Directory traversal visits entries in sorted order, where the source
//!   follows `std::filesystem::recursive_directory_iterator`, whose order the
//!   file system decides. The archive entries and the index lines are
//!   therefore deterministic.
//! - Entry names are checked for escaping the extraction directory component
//!   by component. The source compares the normalised path with the base as
//!   *strings*, so an entry such as `../<base name>x/f` passes its check and
//!   lands in a sibling directory whose name extends the base's.
//! - Names in the JSON index are escaped as JSON strings. The source writes
//!   them verbatim, which yields invalid JSON for a name containing `"` or
//!   `\`; ordinary names produce identical bytes.
//! - The source's `Exception::NotImplemented` for a build without libzip has
//!   no counterpart: this module exists exactly when ZIP support is compiled
//!   in.

use super::zip_ifstream::{Entry, EntryMeta};
use crate::system::file::TempDir;
use crate::{Error, Result};
use std::fs::{self, File, OpenOptions};
use std::io::{self, Read, Write};
use std::path::{Component, Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

/// Most entries an archive may have before any of its entries is extracted.
///
/// Native preflight ceiling; the source has none. Extraction work is also
/// bounded per entry by the size the archive declares, which decoding may
/// not exceed.
pub const MAX_ENTRIES: usize = 1 << 20;

/// Subdirectory of the temporary directory that [`ZipArchiveFile::unzip_directory`]
/// unpacks into, as the source names it.
pub const UNPACK_DIRECTORY: &str = "parquet_unpacked";

/// Namespace for the ZIP helpers; every function is associated, as the
/// source's class has only static members.
#[derive(Debug, Clone, Copy, Default)]
pub struct ZipArchiveFile;

impl ZipArchiveFile {
    /// Create a store-only ZIP archive from the files below `directory_path`.
    ///
    /// Every regular file below the directory, recursively, becomes one entry
    /// named by its path relative to `directory_path` with `/` separators, and
    /// is stored without compression. Directories are not stored as entries,
    /// so empty directories are lost; a symbolic link to a file is stored as
    /// the file's content, and a symbolic link to a directory is skipped, as
    /// the source's iterator does. An existing `output_zip` is replaced. Entry
    /// timestamps are the files' modification times in local time, as libzip
    /// records them (years before 1980 become 1980-01-01).
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the directory cannot be listed or a file cannot be
    /// read, and [`Error::InvalidValue`] when the archive cannot be written,
    /// where the source throws `Exception::InvalidValue`. The source's
    /// header names `Exception::FileNotWritable`; its implementation throws
    /// `InvalidValue`, which this follows. On error `output_zip` is untouched.
    pub fn zip_directory(
        directory_path: impl AsRef<Path>,
        output_zip: impl AsRef<Path>,
    ) -> Result<()> {
        let directory = directory_path.as_ref();
        let output = output_zip.as_ref();
        // List first, so that an output placed inside the directory, or its
        // temporary sibling, is never archived into itself.
        let files = list_files(directory)?;
        publish(output, |file| {
            let mut writer = zip::ZipWriter::new(file);
            for (name, path) in &files {
                add_file(&mut writer, name, path).map_err(|error| {
                    Error::InvalidValue(format!(
                        "Failed to add file to zip: {}: {error}",
                        path.display()
                    ))
                })?;
            }
            writer.finish().map_err(|error| {
                Error::InvalidValue(format!(
                    "Failed to finalize zip archive: {}: {error}",
                    output.display()
                ))
            })
        })
    }

    /// Unpack a ZIP archive into a new temporary directory and return the
    /// directory holding its contents.
    ///
    /// When `input_path` is already a directory it is returned unchanged and
    /// `temp_dir` is left alone. Otherwise a new [`TempDir`] replaces
    /// `*temp_dir` (dropping, and so deleting, any previous one), the archive
    /// is unpacked into its subdirectory [`UNPACK_DIRECTORY`], and that
    /// subdirectory is returned. Directory entries create directories; file
    /// entries create their parent directories as needed.
    ///
    /// # Errors
    ///
    /// - [`Error::Io`] of kind `NotFound` when `input_path` is not a readable
    ///   file, where the source throws `Exception::FileNotFound`.
    /// - [`Error::InvalidValue`] when the archive cannot be opened, has more
    ///   than [`MAX_ENTRIES`] entries, has an entry with an absolute path or
    ///   one that would extract outside the target directory, or has an entry
    ///   that cannot be decoded, as the source's `Exception::InvalidValue`.
    ///   Entry names are all checked before anything is extracted.
    /// - [`Error::Io`] when an extracted file cannot be written, where the
    ///   source throws `Exception::FileNotWritable`.
    ///
    /// On error after the temporary directory was created it stays in
    /// `*temp_dir` with whatever was extracted, as in the source.
    pub fn unzip_directory(
        input_path: impl AsRef<Path>,
        temp_dir: &mut Option<TempDir>,
    ) -> Result<PathBuf> {
        let input = input_path.as_ref();
        if crate::system::file::is_directory(input) {
            return Ok(input.to_path_buf());
        }
        if !crate::system::file::readable(input) {
            return Err(not_found(input));
        }
        let owner = temp_dir.insert(TempDir::new(false)?);
        let unpack = owner.path().join(UNPACK_DIRECTORY);
        fs::create_dir_all(&unpack)?;
        let mut archive = open_archive(input)?;
        let names = entry_names(&archive, input)?;
        let mut plan = Vec::with_capacity(names.len());
        for (index, name) in names.iter().enumerate() {
            let target = target_path(&unpack, name)?;
            let meta = if name.ends_with('/') {
                None
            } else {
                Some(entry_meta(&mut archive, input, index, name)?)
            };
            plan.push((name, target, meta));
        }
        drop(archive);
        for (name, target, meta) in plan {
            match meta {
                None => fs::create_dir_all(&target)?,
                Some(meta) => extract(input, &meta, name, &target)?,
            }
        }
        Ok(unpack)
    }

    /// Add or replace entry `entry_name` of the archive at `archive_path` with
    /// the contents of the file `source_file_path`, stored without compression.
    ///
    /// The archive is created when it does not exist; an existing empty file
    /// counts as an empty archive, as libzip treats it. An existing entry of
    /// that name keeps its position, and a new one is appended. The other
    /// entries are copied unchanged, compressed data included. The source
    /// file is streamed, never held in memory. `entry_name` is used verbatim:
    /// use `/` separators.
    ///
    /// # Errors
    ///
    /// - [`Error::Io`] when the source file cannot be opened or read (kind
    ///   `NotFound` for a missing one). The source header documents
    ///   `Exception::FileNotFound`; its implementation throws
    ///   `Exception::InvalidValue`.
    /// - [`Error::InvalidValue`] when the archive cannot be opened or written,
    ///   as the source.
    ///
    /// On error the archive is untouched.
    pub fn add_or_replace_from_file(
        archive_path: impl AsRef<Path>,
        entry_name: &str,
        source_file_path: impl AsRef<Path>,
    ) -> Result<()> {
        let archive_path = archive_path.as_ref();
        let source = source_file_path.as_ref();
        File::open(source)?;
        let source = source.to_path_buf();
        rewrite(archive_path, entry_name, "add/replace", move |writer| {
            add_file(writer, entry_name, &source)
        })
    }

    /// The names of the entries of the archive at `archive_path`, in
    /// central-directory order, directory entries included.
    ///
    /// Returns an empty list when the archive cannot be read, as the source
    /// does; this function never fails.
    pub fn list_entries(archive_path: impl AsRef<Path>) -> Vec<String> {
        let Ok(file) = File::open(archive_path.as_ref()) else {
            return Vec::new();
        };
        let Ok(archive) = zip::ZipArchive::new(file) else {
            return Vec::new();
        };
        (0..archive.len())
            .filter_map(|i| archive.name_for_index(i).map(str::to_owned))
            .collect()
    }

    /// Write a small JSON index of the archive's entries and their sizes.
    ///
    /// The JSON object maps each entry name to its uncompressed size in bytes,
    /// one entry per line, and lets readers implement index-based access:
    ///
    /// ```text
    /// {
    ///   "a.parquet": 1024,
    ///   "b.parquet": 2048
    /// }
    /// ```
    ///
    /// For a ZIP archive the index is written **into** the archive as the
    /// stored entry `<file name>.idx.json` (for `results.oswpq`,
    /// `results.oswpq.idx.json`), replacing an earlier index, so the archive
    /// stays self-contained. The listing covers every entry present before
    /// the call, directory entries and an earlier index included; an archive
    /// without entries yields `{`, an empty line and `}`.
    ///
    /// When `archive_path` is a directory (an unpacked layout), the index
    /// lists the regular files below it, by relative path, and is written to
    /// the external file whose name is `archive_path` followed by `.idx.json`;
    /// an empty directory yields `{` and `}` on two lines.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the archive cannot be opened or updated,
    /// as the source; [`Error::Io`] when the directory cannot be listed or the
    /// external index cannot be written (the source's
    /// `Exception::FileNotWritable`).
    pub fn write_sidecar_index(archive_path: impl AsRef<Path>) -> Result<()> {
        let archive_path = archive_path.as_ref();
        if crate::system::file::is_directory(archive_path) {
            let mut entries = Vec::new();
            for (name, path) in list_files(archive_path)? {
                entries.push((name, fs::metadata(&path)?.len()));
            }
            let mut json = String::from("{\n");
            for (i, (name, size)) in entries.iter().enumerate() {
                json.push_str(&format!("  {}: {size}", json_string(name)));
                json.push_str(if i + 1 < entries.len() { ",\n" } else { "\n" });
            }
            json.push_str("}\n");
            let mut target = archive_path.as_os_str().to_owned();
            target.push(".idx.json");
            fs::write(PathBuf::from(target), json)?;
            return Ok(());
        }
        let mut archive = open_archive(archive_path)?;
        let mut lines = Vec::new();
        for i in 0..archive.len() {
            let entry = archive
                .by_index_raw(i)
                .map_err(|error| invalid(archive_path, "Failed to open zip archive", error))?;
            lines.push(format!("  {}: {}", json_string(entry.name()), entry.size()));
        }
        drop(archive);
        let json = format!("{{\n{}\n}}\n", lines.join(",\n"));
        let name = archive_path
            .file_name()
            .map(|name| format!("{}.idx.json", name.to_string_lossy()))
            .ok_or_else(|| {
                Error::InvalidValue(format!(
                    "Failed to open zip archive: {}",
                    archive_path.display()
                ))
            })?;
        rewrite(archive_path, &name, "sidecar write", move |writer| {
            writer.start_file(name.as_str(), stored(json.len() as u64, None))?;
            writer.write_all(json.as_bytes())?;
            Ok(())
        })
    }

    /// Extract the single entry `entry_name` of the archive at `archive_path`
    /// into a temporary directory and return the extracted file's path.
    ///
    /// When `*temp_dir` is `None` a new [`TempDir`] is created there; it owns
    /// the extracted file, and any subdirectories of the entry name are
    /// recreated below it. A directory entry (name ending in `/`) creates that
    /// directory and returns its path. This is a lightweight alternative to
    /// [`ZipArchiveFile::unzip_directory`] when only a few files are needed.
    ///
    /// When `archive_path` is a directory (an unpacked layout), nothing is
    /// extracted: the path of `entry_name` inside it is returned when it
    /// exists.
    ///
    /// # Errors
    ///
    /// - [`Error::Io`] of kind `NotFound` when the archive is not readable,
    ///   when it has no entry named `entry_name`, or when the unpacked layout
    ///   has no such file, where the source throws `Exception::FileNotFound`.
    /// - [`Error::InvalidValue`] when the archive or the entry cannot be
    ///   opened or decoded, or when `entry_name` is absolute or would extract
    ///   outside the temporary directory. The source skips that check for an
    ///   unpacked layout; this port applies it there too.
    /// - [`Error::Io`] when the file cannot be written, where the source
    ///   throws `Exception::FileNotWritable`.
    pub fn extract_entry_to_temp_file(
        archive_path: impl AsRef<Path>,
        entry_name: &str,
        temp_dir: &mut Option<TempDir>,
    ) -> Result<PathBuf> {
        let archive_path = archive_path.as_ref();
        if crate::system::file::is_directory(archive_path) {
            let path = target_path(archive_path, entry_name)?;
            if !path.exists() {
                return Err(not_found(&path));
            }
            return Ok(path);
        }
        if !crate::system::file::readable(archive_path) {
            return Err(not_found(archive_path));
        }
        let mut archive = open_archive(archive_path)?;
        let Some(index) = archive.index_for_name(entry_name) else {
            return Err(not_found(Path::new(entry_name)));
        };
        let meta = entry_meta(&mut archive, archive_path, index, entry_name)?;
        drop(archive);
        if temp_dir.is_none() {
            *temp_dir = Some(TempDir::new(false)?);
        }
        let base = temp_dir
            .as_ref()
            .map(|owner| owner.path().to_path_buf())
            .unwrap_or_default();
        let target = target_path(&base, entry_name)?;
        if entry_name.ends_with('/') {
            fs::create_dir_all(&target)?;
            return Ok(target);
        }
        extract(archive_path, &meta, entry_name, &target)?;
        Ok(target)
    }
}

fn not_found(path: &Path) -> Error {
    Error::Io(io::Error::new(
        io::ErrorKind::NotFound,
        format!("file not found: {}", path.display()),
    ))
}

fn invalid(path: &Path, what: &str, reason: impl std::fmt::Display) -> Error {
    Error::InvalidValue(format!("{what}: {}: {reason}", path.display()))
}

fn open_archive(path: &Path) -> Result<zip::ZipArchive<File>> {
    let file =
        File::open(path).map_err(|error| invalid(path, "Failed to open zip archive", error))?;
    zip::ZipArchive::new(file).map_err(|error| invalid(path, "Failed to open zip archive", error))
}

fn entry_names(archive: &zip::ZipArchive<File>, path: &Path) -> Result<Vec<String>> {
    if archive.len() > MAX_ENTRIES {
        return Err(invalid(
            path,
            "Zip archive has too many entries",
            format!("{} > {MAX_ENTRIES}", archive.len()),
        ));
    }
    (0..archive.len())
        .map(|i| {
            archive
                .name_for_index(i)
                .map(str::to_owned)
                .ok_or_else(|| invalid(path, "Failed to read zip entry", i))
        })
        .collect()
}

/// `base` joined with the entry name, refusing absolute names and names that
/// leave `base` after lexical normalisation.
fn target_path(base: &Path, name: &str) -> Result<PathBuf> {
    let relative = Path::new(name);
    let mut depth = 0usize;
    let mut normal = PathBuf::new();
    for component in relative.components() {
        match component {
            Component::Prefix(_) | Component::RootDir => {
                return Err(Error::InvalidValue(format!(
                    "Zip entry has absolute path: {name}"
                )));
            }
            Component::CurDir => {}
            Component::ParentDir => {
                if depth == 0 {
                    return Err(Error::InvalidValue(format!(
                        "Zip entry would extract outside target directory: {name}"
                    )));
                }
                depth -= 1;
                normal.pop();
            }
            Component::Normal(part) => {
                depth += 1;
                normal.push(part);
            }
        }
    }
    Ok(base.join(normal))
}

/// Entry `index`'s metadata, or the source's "Failed to read zip entry".
fn entry_meta(
    archive: &mut zip::ZipArchive<File>,
    path: &Path,
    index: usize,
    name: &str,
) -> Result<EntryMeta> {
    EntryMeta::read(path, archive, index)
        .map_err(|error| invalid(path, &format!("Failed to read zip entry {name}"), error))
}

/// Decode the entry `meta` describes, of the archive at `archive`, into
/// `target`, creating its parent directories.
fn extract(archive: &Path, meta: &EntryMeta, name: &str, target: &Path) -> Result<()> {
    if let Some(parent) = target.parent() {
        fs::create_dir_all(parent)?;
    }
    let file = File::open(archive)
        .map_err(|error| invalid(archive, &format!("Failed to read zip entry {name}"), error))?;
    let mut entry = Entry::from_meta(archive, file, meta)
        .map_err(|error| invalid(archive, &format!("Failed to read zip entry {name}"), error))?;
    let mut output = io::BufWriter::new(File::create(target)?);
    let mut buffer = vec![0u8; 1 << 16];
    loop {
        let n = entry.read_some(&mut buffer).map_err(|error| {
            invalid(archive, &format!("Failed to read zip entry {name}"), error)
        })?;
        if n == 0 {
            break;
        }
        output.write_all(&buffer[..n])?;
    }
    output.flush()?;
    Ok(())
}

/// Regular files below `directory`, recursively, as (entry name, path) in
/// sorted order. Symbolic links are followed for files and skipped for
/// directories, as `recursive_directory_iterator` with `is_directory()` does.
fn list_files(directory: &Path) -> Result<Vec<(String, PathBuf)>> {
    let mut files = Vec::new();
    let mut pending = vec![(String::new(), directory.to_path_buf())];
    while let Some((prefix, path)) = pending.pop() {
        for entry in fs::read_dir(&path)? {
            let entry = entry?;
            let name = entry.file_name().to_string_lossy().into_owned();
            let relative = format!("{prefix}{name}");
            let kind = entry.file_type()?;
            if kind.is_dir() {
                pending.push((format!("{relative}/"), entry.path()));
            } else if kind.is_symlink() {
                if !fs::metadata(entry.path())?.is_dir() {
                    files.push((relative, entry.path()));
                }
            } else {
                files.push((relative, entry.path()));
            }
        }
    }
    files.sort();
    Ok(files)
}

/// Options for a stored entry of `size` bytes, dated `modified` when given.
fn stored(size: u64, modified: Option<std::time::SystemTime>) -> zip::write::SimpleFileOptions {
    let mut options = zip::write::SimpleFileOptions::default()
        .compression_method(zip::CompressionMethod::Stored)
        .large_file(size >= u64::from(u32::MAX));
    if let Some(time) = modified.and_then(dos_time) {
        options = options.last_modified_time(time);
    }
    options
}

/// A modification time as libzip records it: local time, two-second
/// resolution, and only from 1980 on.
fn dos_time(time: std::time::SystemTime) -> Option<zip::DateTime> {
    use chrono::{Datelike, Timelike};
    let local = chrono::DateTime::<chrono::Local>::from(time);
    zip::DateTime::from_date_and_time(
        u16::try_from(local.year()).ok()?,
        local.month() as u8,
        local.day() as u8,
        local.hour() as u8,
        local.minute() as u8,
        local.second() as u8,
    )
    .ok()
}

fn add_file(
    writer: &mut zip::ZipWriter<File>,
    name: &str,
    path: &Path,
) -> std::result::Result<(), zip::result::ZipError> {
    let mut source = File::open(path)?;
    let metadata = source.metadata()?;
    writer.start_file(name, stored(metadata.len(), metadata.modified().ok()))?;
    io::copy(&mut source, writer)?;
    Ok(())
}

/// Rewrite the archive at `archive_path` with entry `name` produced by
/// `write_entry`: in the old entry's position when there is one, else
/// appended. Every other entry is copied raw.
fn rewrite(
    archive_path: &Path,
    name: &str,
    what: &str,
    write_entry: impl FnOnce(
        &mut zip::ZipWriter<File>,
    ) -> std::result::Result<(), zip::result::ZipError>,
) -> Result<()> {
    let failed = |error: &dyn std::fmt::Display| {
        invalid(
            archive_path,
            &format!("Failed to open zip archive for {what}"),
            error,
        )
    };
    let existing = match File::open(archive_path) {
        Ok(file) => {
            if file.metadata()?.len() == 0 {
                None
            } else {
                Some(zip::ZipArchive::new(file).map_err(|error| failed(&error))?)
            }
        }
        Err(error) if error.kind() == io::ErrorKind::NotFound => None,
        Err(error) => return Err(failed(&error)),
    };
    let mut write_entry = Some(write_entry);
    publish(archive_path, |file| {
        let mut writer = zip::ZipWriter::new(file);
        let finalize = |error: zip::result::ZipError| {
            invalid(archive_path, "Failed to finalize zip archive", error)
        };
        if let Some(mut archive) = existing {
            for i in 0..archive.len() {
                let replace = archive.name_for_index(i) == Some(name);
                if replace {
                    if let Some(write) = write_entry.take() {
                        write(&mut writer).map_err(finalize)?;
                    }
                } else {
                    let entry = archive.by_index_raw(i).map_err(|error| {
                        invalid(
                            archive_path,
                            &format!("zip_file_add failed for {what}"),
                            error,
                        )
                    })?;
                    writer.raw_copy_file(entry).map_err(finalize)?;
                }
            }
        }
        if let Some(write) = write_entry.take() {
            write(&mut writer).map_err(finalize)?;
        }
        writer.finish().map_err(finalize)
    })
}

/// Write `destination` through a temporary sibling and rename it into place
/// only after `build` succeeded; the temporary is removed on any error.
fn publish(destination: &Path, build: impl FnOnce(File) -> Result<File>) -> Result<()> {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let parent = match destination.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => parent,
        _ => Path::new("."),
    };
    let mut attempt = 0;
    let (temporary, file) = loop {
        let id = NEXT.fetch_add(1, Ordering::Relaxed);
        let path = parent.join(format!(".openms-zip-{}-{id}.tmp", std::process::id()));
        match OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&path)
        {
            Ok(file) => break (Temporary(path), file),
            Err(error) if error.kind() == io::ErrorKind::AlreadyExists && attempt < 32 => {
                attempt += 1;
            }
            Err(error) => {
                return Err(invalid(destination, "Failed to create zip archive", error));
            }
        }
    };
    let file = build(file)?;
    file.sync_all()?;
    drop(file);
    fs::rename(&temporary.0, destination)
        .map_err(|error| invalid(destination, "Failed to finalize zip archive", error))?;
    std::mem::forget(temporary);
    Ok(())
}

struct Temporary(PathBuf);

impl Drop for Temporary {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}

/// A JSON string literal.
fn json_string(text: &str) -> String {
    let mut out = String::with_capacity(text.len() + 2);
    out.push('"');
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn target_paths_stay_below_the_base() {
        let base = Path::new("/tmp/base");
        assert_eq!(target_path(base, "a/b.txt").unwrap(), base.join("a/b.txt"));
        assert_eq!(target_path(base, "a/../b.txt").unwrap(), base.join("b.txt"));
        assert_eq!(target_path(base, "./a/").unwrap(), base.join("a"));
        assert!(target_path(base, "../basex/f").is_err());
        assert!(target_path(base, "a/../../f").is_err());
        assert!(target_path(base, "/etc/passwd").is_err());
    }

    #[test]
    fn json_strings_are_escaped() {
        assert_eq!(json_string("plain.parquet"), "\"plain.parquet\"");
        assert_eq!(json_string("a\"b\\c\n"), "\"a\\\"b\\\\c\\n\"");
        assert_eq!(json_string("\u{1}"), "\"\\u0001\"");
    }
}
