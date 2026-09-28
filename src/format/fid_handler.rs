// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Read-only `fid` intensity file of a Bruker XMass (Flex series) spectrum,
//! the source `Internal::FidHandler` (`FORMAT/HANDLERS/FidHandler.h`).
//!
//! A `fid` file is a bare array of 32-bit little-endian signed integers, one
//! intensity per time-of-flight point, with no header. The source documents
//! the class as "only needed for XMassFile";
//! [`crate::format::xmass_file`] is its one consumer here too.
//!
//! The source class *is* an `std::ifstream` (public inheritance), so its whole
//! stream interface is reachable; the port wraps any [`Read`](std::io::Read) instead and
//! exposes the two members the header declares plus a native end-of-data
//! query. See `docs/XMASS_FILE_SUPPORT.md`.
//!
//! # Reading past the end
//!
//! `getIntensity` reads 4 bytes into a zero-initialised `Int32` with
//! `std::istream::read` and does not check the stream state. At the end of the
//! file it therefore returns 0, and when 1 to 3 bytes remain it returns the
//! value those bytes form in the low-order positions (on a little-endian
//! host, which the Linux x86_64 Release build is). Either way the point index
//! still advances. [`FidHandler::intensity`](crate::format::fid_handler::FidHandler::intensity) refuses such a read with
//! [`Error::Parse`], because a zero the file never held is not data;
//! [`FidHandler::intensity_source`](crate::format::fid_handler::FidHandler::intensity_source) reproduces the source value.

use crate::{Error, Result};
use std::io::{BufReader, ErrorKind, Read};
use std::path::Path;

/// Sequential reader of the intensities of a `fid` file.
///
/// The source class `Internal::FidHandler`. The point index starts at 0 and
/// advances by one on every intensity read, successful or not, as in the
/// source.
#[derive(Debug)]
pub struct FidHandler<R: Read = BufReader<std::fs::File>> {
    reader: R,
    /// The source `index_`: the index of the next point to be read.
    index: u64,
}

impl FidHandler {
    /// Open a `fid` file at its first point.
    ///
    /// The source constructor `FidHandler(const std::string& filename)`, which
    /// opens the stream in binary mode and seeks to the start. The source
    /// constructor cannot fail; the caller tests the stream afterwards
    /// (`if (!fid)` in `XMassFile::load`). Here the open error is returned
    /// directly.
    ///
    /// # Errors
    ///
    /// The I/O error when the file cannot be opened; its
    /// [`std::io::ErrorKind`] distinguishes the source's `FileNotFound`,
    /// `FileNotReadable` and `IOException`.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let file = std::fs::File::open(path.as_ref())?;
        Ok(Self::new(BufReader::new(file)))
    }
}

impl<R: Read> FidHandler<R> {
    /// Read intensities from an already opened stream, starting at point 0.
    ///
    /// Native entry point; the source only opens files by name.
    pub fn new(reader: R) -> Self {
        Self { reader, index: 0 }
    }

    /// The index of the point the next read returns, without moving.
    ///
    /// The source `getIndex()`.
    pub fn index(&self) -> u64 {
        self.index
    }

    /// The intensity of the current point, then advance to the next.
    ///
    /// The source `getIntensity()`: the 32-bit little-endian signed value, with
    /// a negative value clamped to 0 (`(result > 0) ? result : 0`). The source
    /// swaps bytes on a big-endian host; `i32::from_le_bytes` makes that
    /// unconditional.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when fewer than 4 bytes remain, where the source
    /// returns 0 or a partial value (see the module documentation and
    /// [`Self::intensity_source`]); the I/O error of any other read failure.
    /// The index advances in every case, as in the source.
    pub fn intensity(&mut self) -> Result<u32> {
        let (value, complete) = self.read_point()?;
        if !complete {
            return Err(Error::Parse {
                line: 0,
                message: format!(
                    "fid file ends before point {}: fewer than 4 bytes remain",
                    self.index - 1
                ),
            });
        }
        Ok(value)
    }

    /// The intensity of the current point exactly as the source returns it,
    /// then advance to the next.
    ///
    /// Identical to [`Self::intensity`] for a complete point. When fewer than 4
    /// bytes remain, returns what `std::istream::read` leaves in the
    /// zero-initialised `Int32` on the Linux x86_64 Release build: the
    /// remaining bytes in the low-order positions, clamped at 0 like any
    /// value, or 0 when nothing remains.
    ///
    /// # Errors
    ///
    /// The I/O error of a read failure other than the end of the data.
    pub fn intensity_source(&mut self) -> Result<u32> {
        Ok(self.read_point()?.0)
    }

    /// One `getIntensity` read: the clamped value and whether 4 bytes were read.
    fn read_point(&mut self) -> Result<(u32, bool)> {
        let mut bytes = [0u8; 4];
        let mut filled = 0;
        while filled < bytes.len() {
            match self.reader.read(&mut bytes[filled..]) {
                Ok(0) => break,
                Ok(read) => filled += read,
                Err(error) if error.kind() == ErrorKind::Interrupted => {}
                Err(error) => return Err(error.into()),
            }
        }
        self.index += 1;
        let value = i32::from_le_bytes(bytes);
        Ok((u32::try_from(value).unwrap_or(0), filled == bytes.len()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn values_are_little_endian_and_negative_values_clamp() {
        let mut data = Vec::new();
        data.extend_from_slice(&317i32.to_le_bytes());
        data.extend_from_slice(&(-5i32).to_le_bytes());
        data.extend_from_slice(&i32::MAX.to_le_bytes());
        let mut fid = FidHandler::new(data.as_slice());
        assert_eq!(fid.index(), 0);
        assert_eq!(fid.intensity().unwrap(), 317);
        assert_eq!(fid.index(), 1);
        assert_eq!(fid.intensity().unwrap(), 0);
        assert_eq!(fid.intensity().unwrap(), i32::MAX as u32);
        assert_eq!(fid.index(), 3);
    }

    #[test]
    fn short_reads_are_refused_or_reproduced() {
        let data = [0x34u8, 0x12];
        let mut strict = FidHandler::new(&data[..]);
        assert!(matches!(strict.intensity(), Err(Error::Parse { .. })));
        assert_eq!(strict.index(), 1);

        let mut source = FidHandler::new(&data[..]);
        assert_eq!(source.intensity_source().unwrap(), 0x1234);
        assert_eq!(source.intensity_source().unwrap(), 0);
        assert_eq!(source.index(), 2);

        // Three remaining bytes whose top byte would make the value negative
        // cannot occur: the missing high byte is zero, so the value is positive.
        let data = [0xffu8, 0xff, 0xff];
        let mut source = FidHandler::new(&data[..]);
        assert_eq!(source.intensity_source().unwrap(), 0x00ff_ffff);
    }
}
