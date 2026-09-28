// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! One-entry ZIP input for the path-based readers: the `PK` branch of
//! `Internal::XMLFile::parse_`.
//!
//! The source peeks at a file's first two bytes (`FORMAT/XMLFile.cpp:141-165`
//! at core `bc9cc12`); `PK` hands the file to `CompressedInputSource`, whose
//! `makeStream` builds a `ZipInputStream` over `ZipIfstream`
//! (`FORMAT/CompressedInputSource.cpp`, `FORMAT/ZipInputStream.cpp`,
//! `FORMAT/ZipIfstream.cpp`). `ZipIfstream::open` opens the archive with
//! libzip, skips every entry whose name ends in `/`, and requires exactly one
//! remaining entry, which it then streams decompressed. Every reader built on
//! `XMLFile::parse_` inherits this: featureXML, mzML, consensusXML, idXML,
//! ParamXML and the other XML dialects, and so does every reader here that opens
//! its path through `path_io::open`.
//!
//! # Preserved source conventions
//!
//! - Detection is by content, the two bytes `PK`, not by the file name, so an
//!   archive that holds only its end record (`PK\x05\x06`, an empty archive)
//!   also takes this branch.
//! - Entries whose name ends in `/` are directories and are skipped.
//! - No file entry is refused as the source's `Exception::FileEmpty` ("ZIP
//!   archive contains no file entries"); more than one as its
//!   `Exception::ParseError`, whose message names every entry in archive order
//!   ("ZIP archive contains 2 file entries; expected exactly 1: a b").
//! - The methods the reference build's libzip decodes are decoded: stored (0),
//!   deflate (8) and bzip2 (12). Its libzip links zlib and libbz2 only
//!   (`ldd deps/lib/libzip.so.5` of the Release install), so it has no LZMA,
//!   Zstandard or Deflate64 decoder either; those, and an encrypted entry
//!   (libzip needs a password the source never supplies), are refused.
//! - The entry's CRC-32 and size are checked at its end, as libzip's
//!   `zip_fread` does; a mismatch is an input error.
//!
//! # Native differences
//!
//! - The archive structure (end record, ZIP64, central directory, local
//!   headers) is read with the `zip` crate; the entry's bytes are then streamed
//!   from its data offset through `flate2` or `bzip2`, the decoders the port
//!   already uses for gzip and bzip2 input. The source reads the same fields
//!   through libzip.
//! - Entry names are the `zip` crate's decoding (UTF-8 when the entry says so,
//!   otherwise CP437); libzip guesses the encoding. Names only appear in the
//!   error for an archive with several entries.
//! - Errors are [`Error`] values rather than the source's exception types; the
//!   messages are the source's where it has one.
//!
//! # Evidence
//!
//! Tier 1: `../oracle/B8a-zip-input` loads the fixtures of
//! `tests/data/zip_input/` with the Linux x86_64 Release build's
//! `FeatureXMLFile` and `MzMLFile` (two identical runs on ibminode06);
//! `tests/zip_input.rs` compares.

use crate::{Error, Result};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use zip::{CompressionMethod, ZipArchive};

/// Open the single file entry of the ZIP archive `file` as a decompressed,
/// CRC-checked stream.
///
/// `file` may be positioned anywhere; the archive is read from its end record.
/// `name` is the path shown in errors.
///
/// # Errors
///
/// Returns [`Error::Parse`] when the file is not a readable ZIP archive or
/// holds no file entry or more than one, [`Error::Unsupported`] for an
/// encrypted entry or a compression method other than stored, deflate and
/// bzip2, and [`Error::Io`] when reading fails. The source throws
/// `FileNotReadable`, `FileEmpty` and `ParseError` for these; [`Error`] has no
/// file-empty or not-readable variant, so a TOPP tool exits
/// `INPUT_FILE_CORRUPT` or `INCOMPATIBLE_INPUT_DATA` where the source exits
/// `INPUT_FILE_EMPTY` or `INPUT_FILE_NOT_READABLE`.
pub(crate) fn open_single_entry(file: File, name: &str) -> Result<Box<dyn BufRead>> {
    let mut archive = ZipArchive::new(BufReader::new(file)).map_err(|error| Error::Parse {
        line: 0,
        message: format!("the file '{name}' is not a readable ZIP archive: {error}"),
    })?;
    let mut entries = Vec::new();
    let mut index = None;
    for position in 0..archive.len() {
        let Some(entry) = archive.name_for_index(position) else {
            continue;
        };
        if entry.ends_with('/') {
            continue;
        }
        entries.push(entry.to_owned());
        index = Some(position);
    }
    let (index, entry_name) = match (index, entries.first(), entries.len()) {
        (Some(index), Some(entry_name), 1) => (index, entry_name.clone()),
        (_, _, 0) => {
            return Err(Error::Parse {
                line: 0,
                message: "ZIP archive contains no file entries".into(),
            });
        }
        (_, _, count) => {
            let mut message =
                format!("ZIP archive contains {count} file entries; expected exactly 1:");
            for entry in &entries {
                message.push(' ');
                message.push_str(entry);
            }
            return Err(Error::Parse { line: 0, message });
        }
    };
    let entry = archive.by_index_raw(index).map_err(|error| Error::Parse {
        line: 0,
        message: format!("the ZIP entry of '{name}' is not readable: {error}"),
    })?;
    if entry.encrypted() {
        return Err(Error::Unsupported(format!(
            "the ZIP entry '{entry_name}' of '{name}' is encrypted"
        )));
    }
    let method = entry.compression();
    let layout = EntryLayout {
        data_start: entry.data_start(),
        compressed_size: entry.compressed_size(),
        size: entry.size(),
        crc32: entry.crc32(),
    };
    drop(entry);
    let mut reader = archive.into_inner();
    reader.seek(SeekFrom::Start(layout.data_start))?;
    let compressed = reader.take(layout.compressed_size);
    let decoded: Box<dyn Read> = if method == CompressionMethod::STORE {
        Box::new(compressed)
    } else if method == CompressionMethod::DEFLATE {
        Box::new(flate2::read::DeflateDecoder::new(compressed))
    } else if method == CompressionMethod::BZIP2 {
        Box::new(bzip2::read::BzDecoder::new(compressed))
    } else {
        return Err(Error::Unsupported(format!(
            "the ZIP entry '{entry_name}' of '{name}' uses compression method {method:?}, \
             which the source's libzip cannot decode either"
        )));
    };
    Ok(Box::new(BufReader::new(CheckedEntry {
        inner: decoded,
        crc: flate2::Crc::new(),
        read: 0,
        layout,
    })))
}

/// Where an entry's bytes lie and what they must decode to, from the central
/// directory.
struct EntryLayout {
    data_start: u64,
    compressed_size: u64,
    size: u64,
    crc32: u32,
}

/// A decoded entry that checks its declared size and CRC-32, as libzip does
/// before it reports the end of an entry.
struct CheckedEntry {
    inner: Box<dyn Read>,
    crc: flate2::Crc,
    read: u64,
    layout: EntryLayout,
}

impl Read for CheckedEntry {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let count = self.inner.read(buffer)?;
        if count == 0 {
            if !buffer.is_empty()
                && (self.read != self.layout.size || self.crc.sum() != self.layout.crc32)
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "ZIP entry decodes to {} bytes with CRC-32 {:08x}; the archive declares {} bytes with CRC-32 {:08x}",
                        self.read,
                        self.crc.sum(),
                        self.layout.size,
                        self.layout.crc32
                    ),
                ));
            }
            return Ok(0);
        }
        let chunk = buffer.get(..count).unwrap_or(buffer);
        self.crc.update(chunk);
        self.read = self.read.saturating_add(count as u64);
        if self.read > self.layout.size {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "ZIP entry decodes to more than the {} bytes the archive declares",
                    self.layout.size
                ),
            ));
        }
        Ok(count)
    }
}
