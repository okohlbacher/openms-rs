// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Pull-style gzip decompression with zlib `gzread` semantics (`GzipIfstream.h`).
//!
//! The source wraps zlib's `gzopen`/`gzread`/`gzeof`. This port keeps the
//! observable behaviour of that pairing rather than of a Rust gzip reader:
//!
//! - **Transparent input.** A file that does not begin with the gzip magic
//!   bytes `1f 8b` is returned unchanged, exactly as `gzread` copies a plain
//!   file. An empty file is an empty stream.
//! - **Concatenated members** are decoded one after another.
//! - **Trailing data** after a complete member that does not start with the
//!   gzip magic bytes is ignored, as zlib ignores "trailing garbage".
//! - **Closing.** The stream closes itself when a [`GzipIfstream::read`]
//!   request cannot be filled completely: `gzeof` reports end of file only
//!   after a read tried to go past it. A request that ends exactly at the end
//!   of the data leaves the stream open; the next read returns `0` and closes.
//! - **Errors are sticky, one call late.** When decompression fails after
//!   some bytes of a request were produced, those bytes are returned and the
//!   failure is reported by the next call, which is when `gzread` returns `-1`.
//!
//! The gzip framing (RFC 1952 header, CRC-32 and length trailer) is parsed
//! here; the deflate payload is decoded by the `flate2` crate's raw inflater,
//! which replaces zlib.

use crate::{Error, Result};
use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

/// Input buffer size. The source calls `gzbuffer(gzfile_, 524288)` after
/// opening, so zlib reads the file in 512 KiB blocks.
const INPUT_BYTES: usize = 524_288;

/// Decompresses files in the gzip format (`*.gz`), read in caller-sized pieces.
///
/// Typical use: [`GzipIfstream::from_path`] (or [`GzipIfstream::new`] followed
/// by [`GzipIfstream::open`]), then [`GzipIfstream::read`] until it returns
/// fewer bytes than requested, at which point the stream has closed itself:
/// [`GzipIfstream::is_open`] is `false` and [`GzipIfstream::stream_end`] is
/// `true`. See the [module documentation](self) for the zlib semantics kept.
///
/// The source suppresses copy construction and assignment; this type is
/// neither `Clone` nor `Copy`. It implements [`Read`] for use with the
/// standard library, where a closed stream reads as end of file.
///
/// # Truncated input
///
/// zlib treats a gzip member that ends before its deflate stream or trailer
/// is complete as a soft error (`Z_BUF_ERROR`, "unexpected end of file"):
/// `gzread` returns the bytes decoded so far, then end of file, and the source
/// never looks at `gzerror`. A truncated download therefore reads as a shorter
/// valid file. This port refuses such input by default with [`Error::Parse`];
/// [`GzipIfstream::set_accept_truncated`] selects the source behaviour.
#[derive(Debug)]
pub struct GzipIfstream {
    decoder: Option<Box<Decoder>>,
    stream_at_end: bool,
    accept_truncated: bool,
}

impl Default for GzipIfstream {
    fn default() -> Self {
        Self::new()
    }
}

impl GzipIfstream {
    /// A closed stream: [`GzipIfstream::is_open`] is `false` and
    /// [`GzipIfstream::stream_end`] is `true`, as after the source's default
    /// constructor.
    pub fn new() -> Self {
        Self {
            decoder: None,
            stream_at_end: true,
            accept_truncated: false,
        }
    }

    /// Open `path` for decompression in one step; see [`GzipIfstream::open`].
    ///
    /// # Errors
    ///
    /// As [`GzipIfstream::open`].
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self> {
        let mut stream = Self::new();
        stream.open(path)?;
        Ok(stream)
    }

    /// Open a file for reading (decompression); a previously open file is
    /// closed first.
    ///
    /// No data is read yet: like `gzopen`, opening only acquires the file, and
    /// whether the content is gzip is decided by the first
    /// [`GzipIfstream::read`].
    ///
    /// # Errors
    ///
    /// [`Error::Io`] with the operating system's error when the file cannot be
    /// opened: kind `NotFound` where the source throws
    /// `Exception::FileNotFound`, `PermissionDenied` where it throws
    /// `Exception::FileNotReadable`, anything else where it throws
    /// `Exception::IOException`. The stream is then closed.
    pub fn open(&mut self, path: impl AsRef<Path>) -> Result<()> {
        self.close();
        let file = File::open(path.as_ref())?;
        self.decoder = Some(Box::new(Decoder::new(file)));
        self.stream_at_end = false;
        Ok(())
    }

    /// Read up to `buf.len()` decompressed bytes into `buf` and return how many
    /// were written.
    ///
    /// The output is a raw byte stream, not a terminated string. A return value
    /// smaller than `buf.len()` means the end of the data was reached and the
    /// stream has been closed; check [`GzipIfstream::is_open`] before reading
    /// again. An empty `buf` returns `0` and changes nothing.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] when no file is open (never opened, closed, or
    ///   read to its end), where the source throws
    ///   `Exception::IllegalArgument`.
    /// - [`Error::Parse`] when the data is not valid gzip: a bad header, a
    ///   corrupt deflate stream, a CRC-32 or length mismatch, or (unless
    ///   [`GzipIfstream::set_accept_truncated`]) a truncated member. The source
    ///   header documents `Exception::ConversionError`, but the implementation
    ///   throws `Exception::ParseError`, which this follows.
    /// - [`Error::Io`] when reading the file fails; the source reports that as
    ///   a `ParseError` too.
    ///
    /// Either failure closes the stream. When the failure is detected after
    /// part of the request was already decoded, those bytes are returned first
    /// and the error is reported by the next call, as `gzread` does.
    pub fn read(&mut self, buf: &mut [u8]) -> Result<usize> {
        let accept_truncated = self.accept_truncated;
        let Some(decoder) = self.decoder.as_mut() else {
            return Err(Error::InvalidValue(
                "no file for decompression initialized".into(),
            ));
        };
        if let Some(error) = decoder.error.take() {
            self.close();
            return Err(error);
        }
        let mut filled = 0;
        let mut failure = None;
        while filled < buf.len() {
            match decoder.read_some(&mut buf[filled..], accept_truncated) {
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
                decoder.error = Some(error);
                Ok(filled)
            }
            None => {
                if filled < buf.len() {
                    self.close();
                }
                Ok(filled)
            }
        }
    }

    /// Whether the end of the stream was reached or no file is open, that is,
    /// whether [`GzipIfstream::read`] can no longer be used.
    pub fn stream_end(&self) -> bool {
        self.stream_at_end
    }

    /// Whether a file is open.
    pub fn is_open(&self) -> bool {
        self.decoder.is_some()
    }

    /// Close the current file; a no-op when none is open. Afterwards
    /// [`GzipIfstream::is_open`] is `false` and [`GzipIfstream::stream_end`]
    /// is `true`.
    pub fn close(&mut self) {
        self.decoder = None;
        self.stream_at_end = true;
    }

    /// Whether a truncated gzip member reads as end of file (`true`, the
    /// source's zlib behaviour) or fails with [`Error::Parse`] (`false`, the
    /// default). Native; the source has no such switch.
    pub fn set_accept_truncated(&mut self, accept: bool) {
        self.accept_truncated = accept;
    }

    /// See [`GzipIfstream::set_accept_truncated`].
    pub fn accept_truncated(&self) -> bool {
        self.accept_truncated
    }
}

impl Read for GzipIfstream {
    /// A closed stream reads as end of file instead of failing, as the
    /// [`Read`] contract expects after the end was reached.
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if !self.is_open() {
            return Ok(0);
        }
        GzipIfstream::read(self, buf).map_err(into_io)
    }
}

pub(crate) fn into_io(error: Error) -> io::Error {
    match error {
        Error::Io(error) => error,
        other => io::Error::new(io::ErrorKind::InvalidData, other),
    }
}

fn corrupt(message: impl std::fmt::Display) -> Error {
    Error::Parse {
        line: 0,
        message: format!("error reading from gzip file: {message}"),
    }
}

/// Buffered file input that can guarantee a small lookahead, which a
/// `BufReader` cannot once its buffer is nearly consumed.
#[derive(Debug)]
struct Input {
    file: File,
    buf: Vec<u8>,
    pos: usize,
    len: usize,
    eof: bool,
}

impl Input {
    fn new(file: File) -> Self {
        Self {
            file,
            buf: Vec::new(),
            pos: 0,
            len: 0,
            eof: false,
        }
    }

    fn available(&self) -> &[u8] {
        &self.buf[self.pos..self.len]
    }

    fn consume(&mut self, n: usize) {
        self.pos = (self.pos + n).min(self.len);
    }

    /// Read more input, compacting the buffer first.
    fn refill(&mut self) -> io::Result<()> {
        if self.eof {
            return Ok(());
        }
        if self.buf.is_empty() {
            self.buf = vec![0; INPUT_BYTES];
        }
        if self.pos > 0 {
            self.buf.copy_within(self.pos..self.len, 0);
            self.len -= self.pos;
            self.pos = 0;
        }
        if self.len == self.buf.len() {
            return Ok(());
        }
        loop {
            match self.file.read(&mut self.buf[self.len..]) {
                Ok(0) => {
                    self.eof = true;
                    return Ok(());
                }
                Ok(n) => {
                    self.len += n;
                    return Ok(());
                }
                Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
                Err(error) => return Err(error),
            }
        }
    }

    /// Make at least `n` bytes available unless the file ends first; returns
    /// the number available. `n` must not exceed the buffer size.
    fn ensure(&mut self, n: usize) -> io::Result<usize> {
        while self.len - self.pos < n && !self.eof {
            self.refill()?;
        }
        Ok(self.len - self.pos)
    }

    fn byte(&mut self) -> io::Result<Option<u8>> {
        if self.ensure(1)? == 0 {
            return Ok(None);
        }
        let byte = self.buf[self.pos];
        self.pos += 1;
        Ok(Some(byte))
    }
}

#[derive(Debug)]
enum How {
    /// Look for a gzip header (zlib `LOOK`).
    Look,
    /// Plain input is copied through (zlib `COPY`).
    Copy,
    /// Inside a gzip member (zlib `GZIP`).
    Gzip(Box<Member>),
    /// A member's deflate data ended; its trailer must match this CRC-32 and
    /// length.
    Trailer(u32, u32),
    /// No more data.
    Done,
}

struct Member {
    inflate: flate2::Decompress,
    crc: flate2::Crc,
}

impl std::fmt::Debug for Member {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Member")
            .field("total_in", &self.inflate.total_in())
            .field("total_out", &self.inflate.total_out())
            .finish()
    }
}

/// What went wrong inside a member.
enum Fault {
    Truncated,
    Corrupt(String),
    Io(io::Error),
}

impl From<io::Error> for Fault {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

#[derive(Debug)]
struct Decoder {
    input: Input,
    how: How,
    /// zlib `state->direct`: true until a gzip member was seen, so that
    /// non-gzip data is copied at the start but ignored after a member.
    direct: bool,
    /// Failure to report on the next call.
    error: Option<Error>,
}

impl Decoder {
    fn new(file: File) -> Self {
        Self {
            input: Input::new(file),
            how: How::Look,
            direct: true,
            error: None,
        }
    }

    /// Produce at least one byte, or `0` at the end of the data.
    fn read_some(&mut self, out: &mut [u8], accept_truncated: bool) -> Result<usize> {
        match self.step(out) {
            Ok(n) => Ok(n),
            Err(Fault::Truncated) => {
                self.how = How::Done;
                if accept_truncated {
                    Ok(0)
                } else {
                    Err(corrupt("unexpected end of file"))
                }
            }
            Err(Fault::Corrupt(message)) => {
                self.how = How::Done;
                Err(corrupt(message))
            }
            Err(Fault::Io(error)) => {
                self.how = How::Done;
                Err(Error::Io(error))
            }
        }
    }

    fn step(&mut self, out: &mut [u8]) -> std::result::Result<usize, Fault> {
        loop {
            match &mut self.how {
                How::Done => return Ok(0),
                How::Trailer(crc, length) => {
                    let (crc, length) = (*crc, *length);
                    read_trailer(&mut self.input, crc, length)?;
                    self.how = How::Look;
                }
                How::Look => {
                    let available = self.input.ensure(2)?;
                    if available == 0 {
                        self.how = How::Done;
                        continue;
                    }
                    let magic = self.input.available();
                    if available > 1 && magic[0] == 0x1f && magic[1] == 0x8b {
                        self.input.consume(2);
                        read_header(&mut self.input)?;
                        self.how = How::Gzip(Box::new(Member {
                            inflate: flate2::Decompress::new(false),
                            crc: flate2::Crc::new(),
                        }));
                        self.direct = false;
                    } else if self.direct {
                        self.how = How::Copy;
                    } else {
                        // Trailing garbage after a member is ignored.
                        self.how = How::Done;
                    }
                }
                How::Copy => {
                    if self.input.available().is_empty() {
                        self.input.refill()?;
                        if self.input.available().is_empty() {
                            self.how = How::Done;
                            continue;
                        }
                    }
                    let available = self.input.available();
                    let n = available.len().min(out.len());
                    out[..n].copy_from_slice(&available[..n]);
                    self.input.consume(n);
                    return Ok(n);
                }
                How::Gzip(member) => {
                    if self.input.available().is_empty() {
                        self.input.refill()?;
                    }
                    let input = self.input.available();
                    let at_end = input.is_empty();
                    let before_in = member.inflate.total_in();
                    let before_out = member.inflate.total_out();
                    let status = member
                        .inflate
                        .decompress(input, out, flate2::FlushDecompress::None)
                        .map_err(|error| Fault::Corrupt(error.to_string()))?;
                    let consumed = (member.inflate.total_in() - before_in) as usize;
                    let produced = (member.inflate.total_out() - before_out) as usize;
                    self.input.consume(consumed);
                    member.crc.update(&out[..produced]);
                    if status == flate2::Status::StreamEnd {
                        // Hand out this call's bytes before the trailer can
                        // fail: zlib decodes a member's data even when its
                        // trailer is cut off.
                        self.how = How::Trailer(member.crc.sum(), member.crc.amount());
                        if produced > 0 {
                            return Ok(produced);
                        }
                        continue;
                    }
                    if produced > 0 {
                        return Ok(produced);
                    }
                    if at_end {
                        return Err(Fault::Truncated);
                    }
                    if consumed == 0 {
                        return Err(Fault::Corrupt("inflate made no progress".into()));
                    }
                }
            }
        }
    }
}

fn need(input: &mut Input) -> std::result::Result<u8, Fault> {
    input.byte()?.ok_or(Fault::Truncated)
}

/// Parse the rest of an RFC 1952 member header after the two magic bytes,
/// with the checks zlib's `inflate` applies to a gzip header.
fn read_header(input: &mut Input) -> std::result::Result<(), Fault> {
    let mut crc = flate2::Crc::new();
    crc.update(&[0x1f, 0x8b]);
    let mut take = |input: &mut Input| -> std::result::Result<u8, Fault> {
        let byte = need(input)?;
        crc.update(&[byte]);
        Ok(byte)
    };
    if take(input)? != 8 {
        return Err(Fault::Corrupt("unknown compression method".into()));
    }
    let flags = take(input)?;
    if flags & 0xe0 != 0 {
        return Err(Fault::Corrupt("unknown header flags set".into()));
    }
    for _ in 0..6 {
        take(input)?; // MTIME, XFL, OS
    }
    if flags & 0x04 != 0 {
        let low = take(input)? as usize;
        let high = take(input)? as usize;
        for _ in 0..(low | (high << 8)) {
            take(input)?;
        }
    }
    for flag in [0x08, 0x10] {
        if flags & flag != 0 {
            while take(input)? != 0 {}
        }
    }
    if flags & 0x02 != 0 {
        let expected = (crc.sum() & 0xffff) as u16;
        let low = need(input)? as u16;
        let high = need(input)? as u16;
        if low | (high << 8) != expected {
            return Err(Fault::Corrupt("header crc mismatch".into()));
        }
    }
    Ok(())
}

/// Check the CRC-32 and ISIZE trailer of a member.
fn read_trailer(input: &mut Input, crc: u32, length: u32) -> std::result::Result<(), Fault> {
    if input.ensure(8)? < 8 {
        return Err(Fault::Truncated);
    }
    let trailer = input.available();
    let stored_crc = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
    let stored_length = u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]);
    input.consume(8);
    if stored_crc != crc {
        return Err(Fault::Corrupt("incorrect data check".into()));
    }
    if stored_length != length {
        return Err(Fault::Corrupt("incorrect length check".into()));
    }
    Ok(())
}
