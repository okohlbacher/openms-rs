// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Streaming reader for the single file inside a ZIP archive (`ZipIfstream.h`).
//!
//! The source is backed by libzip. This port reads the archive's central
//! directory with the `zip` crate and decodes the selected entry itself, with
//! `flate2` for deflate and `bzip2` for bzip2 entries, so that the stream owns
//! its file handle and can be read lazily in caller-sized pieces.
//!
//! This is what makes a ZIP container usable wherever the loaders accept a
//! gzip or bzip2 one: `path_io::open`, and through it every loader and
//! `FileHandler`, recognises a ZIP signature and reads through this type.

use crate::{Error, Result};
use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

/// Streaming reader for the single non-directory entry of a ZIP archive.
///
/// Directory entries (names ending in `/`) are ignored when scanning the
/// archive; the archive must contain exactly one non-directory entry, or
/// [`ZipIfstream::open`] fails. Bytes are then decoded lazily from that entry,
/// and the stream closes itself once the entry has been consumed.
///
/// The source suppresses copy and move; this type is neither `Clone` nor
/// `Copy`. It implements [`Read`], where a closed stream reads as end of file.
///
/// # Supported entries
///
/// Stored, deflate and bzip2 entries are decoded. An entry that uses another
/// method (Deflate64, LZMA, XZ, Zstandard, PPMd) or is encrypted cannot be
/// opened; libzip reads LZMA, XZ and Zstandard entries when it was built with
/// those libraries, so such an archive is a gap of this port, reported as
/// "not readable" as libzip reports a method it lacks.
///
/// The source's `Exception::NotImplemented` for a build without libzip has no
/// counterpart: the type exists exactly when ZIP support is compiled in.
#[derive(Debug)]
pub struct ZipIfstream {
    entry: Option<Box<Entry>>,
    stream_at_end: bool,
}

impl Default for ZipIfstream {
    fn default() -> Self {
        Self::new()
    }
}

impl ZipIfstream {
    /// A closed stream: [`ZipIfstream::is_open`] is `false` and
    /// [`ZipIfstream::stream_end`] is `true`.
    pub fn new() -> Self {
        Self {
            entry: None,
            stream_at_end: true,
        }
    }

    /// Open `path` in one step; see [`ZipIfstream::open`] for the behaviour and
    /// the errors.
    ///
    /// # Errors
    ///
    /// As [`ZipIfstream::open`].
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self> {
        let mut stream = Self::new();
        stream.open(path)?;
        Ok(stream)
    }

    /// Open the ZIP archive at `path` and select its single non-directory
    /// entry for reading.
    ///
    /// Any previously open archive is closed first; the new archive is scanned
    /// for entries and the single non-directory one becomes the streaming
    /// source.
    ///
    /// # Errors
    ///
    /// The stream is closed after any error.
    ///
    /// - [`Error::Io`] of kind `NotFound` when `path` does not exist, where the
    ///   source throws `Exception::FileNotFound`.
    /// - [`Error::Io`] when `path` exists but the archive cannot be opened,
    ///   where the source throws `Exception::FileNotReadable`: the operating
    ///   system's error (for example `PermissionDenied`), or kind
    ///   `InvalidData` when the bytes are not a readable ZIP archive or the
    ///   selected entry cannot be decoded (encrypted, or an unsupported
    ///   compression method).
    /// - [`Error::Io`] of kind `UnexpectedEof` when the archive contains no
    ///   non-directory entry, where the source throws `Exception::FileEmpty`.
    /// - [`Error::Parse`] when the archive contains more than one
    ///   non-directory entry; the message lists the entry names, as the
    ///   source's `Exception::ParseError` does.
    ///
    /// # Notes
    ///
    /// Two entries with the same name are one entry to the `zip` crate, which
    /// keeps the later one; libzip counts both and refuses the archive.
    pub fn open(&mut self, path: impl AsRef<Path>) -> Result<()> {
        self.close();
        let entry = Entry::open(path.as_ref())?;
        self.entry = Some(Box::new(entry));
        self.stream_at_end = false;
        Ok(())
    }

    /// Read up to `buf.len()` decompressed bytes from the open entry into
    /// `buf` and return how many were written.
    ///
    /// Like `zip_fread`, a call fills `buf` completely unless the entry ends
    /// first. The stream closes itself when a call returns `0` because the
    /// entry is exhausted; a call that returns fewer bytes than requested
    /// leaves it open, and the next call returns `0` and closes it. (The source
    /// header says a short read already closes the stream; its implementation
    /// closes only on `0`, which this follows.) An empty `buf` returns `0`
    /// without touching the entry.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] when no entry is open (never opened, closed,
    ///   or drained by a previous read), where the source throws
    ///   `Exception::IllegalArgument`. An empty `buf` is checked first and
    ///   returns `0` even then, as in the source.
    /// - [`Error::Parse`] when the entry data is corrupt, ends early, or fails
    ///   its CRC-32 or size check, where `zip_fread` fails and the source
    ///   throws `Exception::ParseError`.
    /// - [`Error::Io`] when reading the archive file fails.
    ///
    /// A failure closes the stream. When it is detected after part of the
    /// request was decoded, those bytes are returned first and the failure is
    /// reported by the next call, as libzip does.
    pub fn read(&mut self, buf: &mut [u8]) -> Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        let Some(entry) = self.entry.as_mut() else {
            return Err(Error::InvalidValue(
                "no file for decompression initialized".into(),
            ));
        };
        if let Some(error) = entry.error.take() {
            self.close();
            return Err(error);
        }
        let mut filled = 0;
        let mut failure = None;
        while filled < buf.len() {
            match entry.read_some(&mut buf[filled..]) {
                Ok(0) => break,
                Ok(n) => filled += n,
                Err(error) => {
                    failure = Some(error);
                    break;
                }
            }
        }
        match failure {
            Some(error) if filled == 0 => {
                self.close();
                Err(error)
            }
            Some(error) => {
                entry.error = Some(error);
                Ok(filled)
            }
            None => {
                if filled == 0 {
                    self.close();
                }
                Ok(filled)
            }
        }
    }

    /// Whether the stream has been closed (end of entry reached or no archive
    /// opened), that is, whether no further [`ZipIfstream::read`] should be
    /// issued.
    pub fn stream_end(&self) -> bool {
        self.stream_at_end
    }

    /// Whether an archive entry is currently open for reading.
    pub fn is_open(&self) -> bool {
        self.entry.is_some()
    }

    /// Close the currently open archive; a no-op when none is open. Afterwards
    /// [`ZipIfstream::is_open`] is `false` and [`ZipIfstream::stream_end`] is
    /// `true`.
    pub fn close(&mut self) {
        self.entry = None;
        self.stream_at_end = true;
    }
}

impl Read for ZipIfstream {
    /// A closed stream reads as end of file instead of failing.
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if !self.is_open() {
            return Ok(0);
        }
        ZipIfstream::read(self, buf).map_err(super::gzip_ifstream::into_io)
    }
}

pub(crate) fn not_readable(path: &Path, reason: impl std::fmt::Display) -> Error {
    Error::Io(io::Error::new(
        io::ErrorKind::InvalidData,
        format!("file not readable: {}: {reason}", path.display()),
    ))
}

fn corrupt(message: impl std::fmt::Display) -> Error {
    Error::Parse {
        line: 0,
        message: format!("error reading from zip entry: {message}"),
    }
}

/// Map a `zip` crate failure: its I/O errors keep their kind, everything else
/// means the archive is not readable.
pub(crate) fn archive_error(path: &Path, error: zip::result::ZipError) -> Error {
    match error {
        zip::result::ZipError::Io(error) => Error::Io(error),
        other => not_readable(path, other),
    }
}

/// The names of the non-directory entries of an archive, in central-directory
/// order, with the index of the last one.
fn file_entries<R: Read + Seek>(archive: &zip::ZipArchive<R>) -> (Vec<String>, Option<usize>) {
    let mut names = Vec::new();
    let mut index = None;
    for i in 0..archive.len() {
        let Some(name) = archive.name_for_index(i) else {
            continue;
        };
        if name.ends_with('/') {
            continue;
        }
        names.push(name.to_owned());
        index = Some(i);
    }
    (names, index)
}

type Compressed = io::Take<BufReader<File>>;

/// Decode one entry's data: stored, raw deflate or bzip2.
pub(crate) enum Decoded {
    Stored(Compressed),
    Deflated(flate2::bufread::DeflateDecoder<Compressed>),
    Bzip2(bzip2::bufread::BzDecoder<Compressed>),
}

impl Decoded {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::Stored(reader) => reader.read(out),
            Self::Deflated(reader) => reader.read(out),
            Self::Bzip2(reader) => reader.read(out),
        }
    }
}

/// Where an entry's data lies in its archive and how to check it.
#[derive(Debug, Clone)]
pub(crate) struct EntryMeta {
    method: zip::CompressionMethod,
    crc32: u32,
    size: u64,
    compressed_size: u64,
    data_start: u64,
}

impl EntryMeta {
    /// Read entry `index`'s metadata, refusing encrypted entries and the
    /// compression methods this port cannot decode.
    pub(crate) fn read(
        path: &Path,
        archive: &mut zip::ZipArchive<File>,
        index: usize,
    ) -> Result<Self> {
        let entry = archive
            .by_index_raw(index)
            .map_err(|error| archive_error(path, error))?;
        if entry.encrypted() {
            return Err(not_readable(path, "the entry is encrypted"));
        }
        let method = entry.compression();
        if !matches!(
            method,
            zip::CompressionMethod::Stored
                | zip::CompressionMethod::Deflated
                | zip::CompressionMethod::Bzip2
        ) {
            return Err(not_readable(
                path,
                format!("unsupported compression method {method:?}"),
            ));
        }
        Ok(Self {
            method,
            crc32: entry.crc32(),
            size: entry.size(),
            compressed_size: entry.compressed_size(),
            data_start: entry.data_start(),
        })
    }
}

/// An entry's decoder together with the checks libzip applies at its end.
pub(crate) struct Entry {
    decoded: Decoded,
    crc: flate2::Crc,
    expected_crc: u32,
    expected_size: u64,
    produced: u64,
    finished: bool,
    error: Option<Error>,
}

impl std::fmt::Debug for Entry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Entry")
            .field("expected_size", &self.expected_size)
            .field("produced", &self.produced)
            .field("finished", &self.finished)
            .finish()
    }
}

impl Entry {
    fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let archive = zip::ZipArchive::new(file).map_err(|error| archive_error(path, error))?;
        let (names, index) = file_entries(&archive);
        let index = match (names.len(), index) {
            (1, Some(index)) => index,
            (0, _) => {
                return Err(Error::Io(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    format!(
                        "file empty: {}: ZIP archive contains no file entries",
                        path.display()
                    ),
                )));
            }
            (count, _) => {
                return Err(Error::Parse {
                    line: 0,
                    message: format!(
                        "ZIP archive contains {count} file entries; expected exactly 1: {}",
                        names.join(" ")
                    ),
                });
            }
        };
        Self::from_archive(path, archive, index)
    }

    /// Stream entry `index` of `archive`, taking over its file handle.
    pub(crate) fn from_archive(
        path: &Path,
        mut archive: zip::ZipArchive<File>,
        index: usize,
    ) -> Result<Self> {
        let meta = EntryMeta::read(path, &mut archive, index)?;
        Self::from_meta(path, archive.into_inner(), &meta)
    }

    /// Stream the entry `meta` describes from `file`, a handle on its archive.
    pub(crate) fn from_meta(path: &Path, mut file: File, meta: &EntryMeta) -> Result<Self> {
        file.seek(SeekFrom::Start(meta.data_start))?;
        let compressed = BufReader::new(file).take(meta.compressed_size);
        let decoded = match meta.method {
            zip::CompressionMethod::Stored => Decoded::Stored(compressed),
            zip::CompressionMethod::Deflated => {
                Decoded::Deflated(flate2::bufread::DeflateDecoder::new(compressed))
            }
            zip::CompressionMethod::Bzip2 => {
                Decoded::Bzip2(bzip2::bufread::BzDecoder::new(compressed))
            }
            other => {
                return Err(not_readable(
                    path,
                    format!("unsupported compression method {other:?}"),
                ));
            }
        };
        Ok(Self {
            decoded,
            crc: flate2::Crc::new(),
            expected_crc: meta.crc32,
            expected_size: meta.size,
            produced: 0,
            finished: false,
            error: None,
        })
    }

    /// Produce at least one byte, or `0` at the end of the entry after its
    /// CRC-32 and size were verified.
    pub(crate) fn read_some(&mut self, out: &mut [u8]) -> Result<usize> {
        if self.finished {
            return Ok(0);
        }
        let n = match self.decoded.read(out) {
            Ok(n) => n,
            Err(error) => {
                self.finished = true;
                return Err(match error.kind() {
                    io::ErrorKind::InvalidInput
                    | io::ErrorKind::InvalidData
                    | io::ErrorKind::UnexpectedEof => corrupt(error),
                    _ => Error::Io(error),
                });
            }
        };
        if n == 0 {
            self.finished = true;
            if self.produced != self.expected_size {
                return Err(corrupt(format!(
                    "entry decodes to {} bytes, the archive declares {}",
                    self.produced, self.expected_size
                )));
            }
            if self.crc.sum() != self.expected_crc {
                return Err(corrupt("CRC error"));
            }
            return Ok(0);
        }
        self.crc.update(&out[..n]);
        self.produced += n as u64;
        if self.produced > self.expected_size {
            self.finished = true;
            return Err(corrupt(format!(
                "entry decodes to more than the {} bytes the archive declares",
                self.expected_size
            )));
        }
        Ok(n)
    }
}
