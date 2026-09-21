// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Shared compressed streams and atomic output publication.

use crate::{Error, Result};
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

/// Detect compression from bytes, independent of the filename extension.
pub(crate) fn open(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut reader = BufReader::new(File::open(path)?);
    let magic = reader.fill_buf()?;
    let gzip = magic.starts_with(&[0x1f, 0x8b]);
    let bzip = magic.starts_with(b"BZh");
    if gzip || bzip {
        #[cfg(feature = "file-compression")]
        {
            if gzip {
                return Ok(Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(
                    reader,
                ))));
            }
            return Ok(Box::new(BufReader::new(bzip2::read::MultiBzDecoder::new(
                reader,
            ))));
        }
        #[cfg(not(feature = "file-compression"))]
        return Err(Error::Unsupported(
            "compressed input requires the file-compression feature".into(),
        ));
    }
    if magic.starts_with(b"PK\x03\x04") {
        return Err(Error::Unsupported(
            "ZIP input is not a single scientific stream".into(),
        ));
    }
    Ok(Box::new(reader))
}

#[cfg(any(feature = "featurexml", feature = "consensusxml"))]
pub(crate) fn store(path: &Path, bytes: &[u8]) -> Result<()> {
    write(path, |writer| {
        writer.write_all(bytes)?;
        Ok(())
    })
}

/// Serialize into a sibling file and publish only after compression and flush succeed.
pub(crate) fn write(path: &Path, save: impl FnOnce(&mut dyn Write) -> Result<()>) -> Result<()> {
    write_as(path, save, true)
}

/// Publish a document that `build` serializes to memory, reporting progress
/// the way the source's `XMLFile::save_` does: the destination is opened
/// first, and only then does the handler's `writeTo` make its calls. `build`
/// therefore runs inside the publication with `progress`, and a destination
/// that cannot be prepared makes no call.
///
/// The outcome is that of building the whole document and then publishing it,
/// as the writers without progress do: when the destination cannot be
/// prepared, the document is built again without reporting, and a refusal
/// from that build wins over the preparation error. A refused document may
/// leave a temporary file created and removed again, which the silent writers
/// never create.
#[cfg(any(feature = "featurexml", feature = "consensusxml", feature = "mzml"))]
pub(crate) fn store_reporting<T>(
    path: &Path,
    mut build: impl FnMut(
        &mut crate::concept::progress_logger::ProgressReporter<'_>,
    ) -> Result<(Vec<u8>, T)>,
    progress: &mut crate::concept::progress_logger::ProgressReporter<'_>,
) -> Result<T> {
    let mut entered = false;
    let mut value = None;
    let published = write(path, |writer| {
        entered = true;
        let (bytes, built) = build(progress)?;
        writer.write_all(&bytes)?;
        value = Some(built);
        Ok(())
    });
    match published {
        Ok(()) => value.ok_or_else(|| Error::InvalidValue("document was not built".into())),
        Err(error) if !entered => {
            build(&mut crate::concept::progress_logger::ProgressReporter::silent())?;
            Err(error)
        }
        Err(error) => Err(error),
    }
}

/// Source IdXMLFile writes plain bytes even when the filename ends in .gz/.bz2.
#[cfg(feature = "idxml")]
pub(crate) fn write_plain(
    path: &Path,
    save: impl FnOnce(&mut dyn Write) -> Result<()>,
) -> Result<()> {
    write_as(path, save, false)
}

fn write_as(
    path: &Path,
    save: impl FnOnce(&mut dyn Write) -> Result<()>,
    compression_suffix: bool,
) -> Result<()> {
    let suffix = if compression_suffix {
        path.extension().and_then(|s| s.to_str()).unwrap_or("")
    } else {
        ""
    };
    if suffix.eq_ignore_ascii_case("zip") {
        return Err(Error::Unsupported("ZIP scientific output".into()));
    }
    let gzip = suffix.eq_ignore_ascii_case("gz");
    let bzip = suffix.eq_ignore_ascii_case("bz2");
    if (gzip || bzip) && !cfg!(feature = "file-compression") {
        return Err(Error::Unsupported(
            "compressed output requires the file-compression feature".into(),
        ));
    }
    let (temporary, mut file) = TemporaryFile::create(path)?;
    {
        let mut writer = BufWriter::new(&mut file);
        #[cfg(feature = "file-compression")]
        if gzip {
            let mut encoder =
                flate2::write::GzEncoder::new(&mut writer, flate2::Compression::default());
            save(&mut encoder)?;
            encoder.finish()?;
        } else if bzip {
            let mut encoder =
                bzip2::write::BzEncoder::new(&mut writer, bzip2::Compression::default());
            save(&mut encoder)?;
            encoder.finish()?;
        } else {
            save(&mut writer)?;
        }
        #[cfg(not(feature = "file-compression"))]
        save(&mut writer)?;
        writer.flush()?;
    }
    file.sync_all()?;
    drop(file);
    std::fs::rename(&temporary.0, path)?;
    Ok(())
}

struct TemporaryFile(PathBuf);
impl TemporaryFile {
    fn create(destination: &Path) -> Result<(Self, File)> {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        Self::create_with_counter(destination, &NEXT)
    }

    fn create_with_counter(destination: &Path, next: &AtomicU64) -> Result<(Self, File)> {
        let parent = destination.parent().unwrap_or_else(|| Path::new("."));
        for _ in 0..32 {
            let id = next.fetch_add(1, Ordering::Relaxed);
            let name = format!(".openms-{}-{id}.tmp", std::process::id());
            // Renaming a temporary file onto itself would succeed, after which
            // its cleanup guard would delete the output. Compare basenames so
            // relative paths and case-insensitive filesystems are safe too.
            if destination
                .file_name()
                .and_then(|s| s.to_str())
                .is_some_and(|s| s.eq_ignore_ascii_case(&name))
            {
                continue;
            }
            let path = parent.join(name);
            match OpenOptions::new().write(true).create_new(true).open(&path) {
                Ok(file) => return Ok((Self(path), file)),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(e.into()),
            }
        }
        Err(Error::Io(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            "cannot allocate temporary output",
        )))
    }
}
impl Drop for TemporaryFile {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_named_like_temporary_survives_cleanup() {
        let directory = std::env::temp_dir().join(format!(
            "openms-handler-temp-collision-{}",
            std::process::id()
        ));
        std::fs::create_dir(&directory).unwrap();
        for uppercase in [false, true] {
            let mut name = format!(".openms-{}-0.tmp", std::process::id());
            if uppercase {
                name.make_ascii_uppercase();
            }
            let destination = directory.join(".").join(name);
            let next = AtomicU64::new(0);
            let (temporary, mut file) =
                TemporaryFile::create_with_counter(&destination, &next).unwrap();
            assert_eq!(next.load(Ordering::Relaxed), 2);
            file.write_all(b"complete output").unwrap();
            drop(file);
            std::fs::rename(&temporary.0, &destination).unwrap();
            drop(temporary);
            assert_eq!(std::fs::read(&destination).unwrap(), b"complete output");
            std::fs::remove_file(destination).unwrap();
        }
        assert_eq!(std::fs::read_dir(&directory).unwrap().count(), 0);
        std::fs::remove_dir(directory).unwrap();
    }
}
