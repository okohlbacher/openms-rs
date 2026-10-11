// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Pull-style bzip2 decompression with libbzip2 `BZ2_bzRead` semantics
//! (`Bzip2Ifstream.h`).
//!
//! The source feeds a `FILE*` to libbzip2's high-level reader. This port
//! re-creates that reader's loop over the `bzip2` crate's low-level
//! decompressor, so that the points at which a read returns, closes and fails
//! are the source's:
//!
//! - input is taken from the file in blocks of `BZ_MAX_UNUSED` (5000) bytes;
//! - a read returns as soon as the request is filled, **or** as soon as the
//!   decompressor reports the end of the stream, which closes the stream even
//!   when the request was filled exactly;
//! - only the **first** bzip2 stream is read. Data after its end-of-stream
//!   marker, including further concatenated streams (as written by `pbzip2`),
//!   is never looked at;
//! - input that ends before the end-of-stream marker fails
//!   (`BZ_UNEXPECTED_EOF`), and any failure discards the bytes decoded during
//!   that call.

use crate::{Error, Result};
use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

/// libbzip2's `BZ_MAX_UNUSED`, the `fread` block size of `BZ2_bzRead`.
const INPUT_BYTES: usize = 5000;

/// Streaming decompressor for bzip2 (`.bz2`) files.
///
/// Reads raw bytes out of a bzip2-compressed file on demand. The type owns
/// the file handle and the decompressor; both are released by
/// [`Bzip2Ifstream::close`] and on drop. The source suppresses copy and
/// assignment; this type is neither `Clone` nor `Copy`.
///
/// Typical lifecycle: [`Bzip2Ifstream::from_path`] (or
/// [`Bzip2Ifstream::open`]) with the path of a `.bz2` file, then
/// [`Bzip2Ifstream::read`] until it returns fewer bytes than requested, at
/// which point the stream has closed itself: [`Bzip2Ifstream::is_open`] is
/// `false` and [`Bzip2Ifstream::stream_end`] is `true`. See the
/// [module documentation](self) for the libbzip2 semantics kept. The type
/// implements [`Read`], where a closed stream reads as end of file.
#[derive(Debug)]
pub struct Bzip2Ifstream {
    reader: Option<Box<Reader>>,
    stream_at_end: bool,
}

impl Default for Bzip2Ifstream {
    fn default() -> Self {
        Self::new()
    }
}

impl Bzip2Ifstream {
    /// A closed stream: [`Bzip2Ifstream::is_open`] is `false` and
    /// [`Bzip2Ifstream::stream_end`] is `true`.
    pub fn new() -> Self {
        Self {
            reader: None,
            stream_at_end: true,
        }
    }

    /// Open `path` for reading in one step; see [`Bzip2Ifstream::open`].
    ///
    /// # Errors
    ///
    /// As [`Bzip2Ifstream::open`].
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self> {
        let mut stream = Self::new();
        stream.open(path)?;
        Ok(stream)
    }

    /// Open `path` for reading; the currently open file is closed first.
    ///
    /// The file is opened in binary mode, so nothing is translated on any
    /// platform. No data is read yet, as `BZ2_bzReadOpen` reads none.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] with the operating system's error when the file cannot be
    /// opened: kind `NotFound` where the source throws
    /// `Exception::FileNotFound`, `PermissionDenied` where it throws
    /// `Exception::FileNotReadable`, anything else where it throws
    /// `Exception::IOException`. The source's `Exception::ConversionError`,
    /// for a failure to initialise the bzip2 reader, cannot occur: creating
    /// the decompressor is infallible here.
    pub fn open(&mut self, path: impl AsRef<Path>) -> Result<()> {
        self.close();
        let file = File::open(path.as_ref())?;
        self.reader = Some(Box::new(Reader::new(file)));
        self.stream_at_end = false;
        Ok(())
    }

    /// Read up to `buf.len()` decompressed bytes into `buf` and return how
    /// many were written.
    ///
    /// The output is a raw byte stream, not a terminated string. When the end
    /// of the compressed stream is reached during this call the stream is
    /// closed ([`Bzip2Ifstream::is_open`] becomes `false`); the return value is
    /// then usually smaller than `buf.len()`, but equals it when the request
    /// ended exactly at the end of the data. An empty `buf` returns `0` and
    /// changes nothing.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] when no file is open (never opened, closed,
    ///   or already read to the end), where the source throws
    ///   `Exception::IllegalArgument`.
    /// - [`Error::Parse`] when bzip2 reports a failure: bad magic bytes,
    ///   corrupt data, a block CRC mismatch, or input that ends before the
    ///   end-of-stream marker. The source throws `Exception::ParseError`.
    /// - [`Error::Io`] when reading the file fails, where the source also
    ///   throws `Exception::ParseError`.
    ///
    /// A failure closes the stream and discards the bytes decoded during the
    /// call, as `BZ2_bzRead` returns `0` with an error.
    pub fn read(&mut self, buf: &mut [u8]) -> Result<usize> {
        let Some(reader) = self.reader.as_mut() else {
            return Err(Error::InvalidValue(
                "no file for decompression initialized".into(),
            ));
        };
        if buf.is_empty() {
            return Ok(0);
        }
        match reader.read(buf) {
            Ok((n, false)) => Ok(n),
            Ok((n, true)) => {
                self.close();
                Ok(n)
            }
            Err(error) => {
                self.close();
                Err(error)
            }
        }
    }

    /// Whether the stream has been closed (end of file reached or no file
    /// opened), that is, whether no further [`Bzip2Ifstream::read`] should be
    /// issued.
    pub fn stream_end(&self) -> bool {
        self.stream_at_end
    }

    /// Whether a file is currently open for reading.
    pub fn is_open(&self) -> bool {
        self.reader.is_some()
    }

    /// Close the currently open file; a no-op when none is open. Afterwards
    /// [`Bzip2Ifstream::is_open`] is `false` and [`Bzip2Ifstream::stream_end`]
    /// is `true`.
    pub fn close(&mut self) {
        self.reader = None;
        self.stream_at_end = true;
    }
}

impl Read for Bzip2Ifstream {
    /// A closed stream reads as end of file instead of failing.
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if !self.is_open() {
            return Ok(0);
        }
        Bzip2Ifstream::read(self, buf).map_err(super::gzip_ifstream::into_io)
    }
}

fn failed(message: impl std::fmt::Display) -> Error {
    Error::Parse {
        line: 0,
        message: format!("bzip2 compression failed: {message}"),
    }
}

/// The state of `BZ2_bzRead`: a `FILE*`, its 5000-byte input block and the
/// decompressor.
struct Reader {
    file: File,
    input: Vec<u8>,
    next: usize,
    available: usize,
    /// A byte read by the end-of-file probe, returned before the file's next
    /// bytes (the source's `fgetc`/`ungetc` pair in `myfeof`).
    pushback: Option<u8>,
    at_eof: bool,
    decompress: bzip2::Decompress,
}

impl std::fmt::Debug for Reader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Reader")
            .field("total_in", &self.decompress.total_in())
            .field("total_out", &self.decompress.total_out())
            .finish()
    }
}

impl Reader {
    fn new(file: File) -> Self {
        Self {
            file,
            input: Vec::new(),
            next: 0,
            available: 0,
            pushback: None,
            at_eof: false,
            decompress: bzip2::Decompress::new(false),
        }
    }

    fn read_byte(&mut self) -> io::Result<Option<u8>> {
        let mut byte = [0u8; 1];
        loop {
            match self.file.read(&mut byte) {
                Ok(0) => return Ok(None),
                Ok(_) => return Ok(Some(byte[0])),
                Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
                Err(error) => return Err(error),
            }
        }
    }

    /// `myfeof`: whether the next `fgetc` would return `EOF`.
    fn feof(&mut self) -> io::Result<bool> {
        if self.pushback.is_some() {
            return Ok(false);
        }
        if self.at_eof {
            return Ok(true);
        }
        match self.read_byte()? {
            Some(byte) => {
                self.pushback = Some(byte);
                Ok(false)
            }
            None => {
                self.at_eof = true;
                Ok(true)
            }
        }
    }

    /// `fread(buf, 1, BZ_MAX_UNUSED, file)`: fill the block as far as the
    /// file allows.
    fn fread(&mut self) -> io::Result<()> {
        if self.input.is_empty() {
            self.input = vec![0; INPUT_BYTES];
        }
        let mut n = 0;
        if let Some(byte) = self.pushback.take() {
            self.input[0] = byte;
            n = 1;
        }
        while n < INPUT_BYTES && !self.at_eof {
            match self.file.read(&mut self.input[n..]) {
                Ok(0) => self.at_eof = true,
                Ok(k) => n += k,
                Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
                Err(error) => return Err(error),
            }
        }
        self.next = 0;
        self.available = n;
        Ok(())
    }

    /// One `BZ2_bzRead` call on a nonempty `out`: the byte count and whether
    /// the end of the stream was reached.
    fn read(&mut self, out: &mut [u8]) -> Result<(usize, bool)> {
        let mut produced = 0;
        loop {
            if self.available == 0 && !self.feof()? {
                self.fread()?;
            }
            let before_in = self.decompress.total_in();
            let before_out = self.decompress.total_out();
            let status = self
                .decompress
                .decompress(
                    &self.input[self.next..self.next + self.available],
                    &mut out[produced..],
                )
                .map_err(failed)?;
            let consumed = (self.decompress.total_in() - before_in) as usize;
            let written = (self.decompress.total_out() - before_out) as usize;
            self.next += consumed;
            self.available -= consumed;
            produced += written;
            if status == bzip2::Status::StreamEnd {
                return Ok((produced, true));
            }
            if self.feof()? && self.available == 0 && produced < out.len() {
                return Err(failed("unexpected end of file"));
            }
            if produced == out.len() {
                return Ok((produced, false));
            }
            if consumed == 0 && written == 0 && self.available > 0 {
                return Err(failed("decompressor made no progress"));
            }
        }
    }
}
