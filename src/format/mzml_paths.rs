// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Filesystem entry points for the represented mzML loading and writing APIs.

use super::{LoadOptions, LoadProgress, MSExperiment, ReadOptions, Result, WriteOptions};
use crate::concept::progress_logger::{ProgressLogger, ProgressReporter};
use std::path::Path;

/// Load a plain, gzip or bzip2 file by content magic, independent of its suffix.
/// Uses source-default scientific loading, including stable point sorting.
/// The older stream `read` API retains its existing input-order behavior.
pub fn load(path: impl AsRef<Path>) -> Result<MSExperiment> {
    load_with_options(path, &LoadOptions::default(), &ReadOptions::default())
}

/// Read with explicit scientific settings and independent XML/binary limits.
/// Input limits count decompressed XML bytes. Successful loads record document path
/// and bounded content-detected type; stream readers leave provenance unset.
pub fn load_with_options(
    path: impl AsRef<Path>,
    scientific: &LoadOptions,
    limits: &ReadOptions,
) -> Result<MSExperiment> {
    load_reporting(path, scientific, limits, None)
}

/// Load a file, reporting progress to `logger` as source `MzMLFile::load` does
/// through its handler.
///
/// The calls are the handler's, on two loggers (`MzMLHandler.cpp:106`, `:135`):
/// the document's section goes to a copy of `logger`, made as
/// [`ProgressLogger::clone`] makes one (a fresh backend of its type, sharing
/// its nesting), and the list sections go to `logger` itself.
///
/// - At `<mzML>`, the copy's `startProgress(0, 1, "loading mzML")` (`:1203`).
/// - At `<spectrumList count>` and `<chromatogramList count>`,
///   `startProgress(0, count, "loading spectra list")` or `"loading
///   chromatogram list"` (`:966`, `:997`); after every `</spectrum>` or
///   `</chromatogram>`, kept or not, `nextProgress()` (`:1443`, `:1483`); at the
///   list's end, `endProgress()` (`:1491`, `:1497`).
/// - At `</mzML>`, the copy's `endProgress(file size)`, so the command
///   backend reports a throughput (`:1524`).
///
/// A metadata-only load stops at the first record list, after the document's
/// section started, and leaves it open, as the source's `EndParsingSoftly`
/// does. The result is the one [`load_with_options`] returns, and so is every
/// error: both run the same code, whose calls go nowhere for
/// [`load_with_options`].
///
/// A failure after a start leaves its section open, as in the source, where
/// the exception bypasses `endProgress`: no `-- done` line is printed, the
/// nesting depth stays deeper, and the next list start on `logger`'s command
/// backend is refused (`StopWatch is already started!`), as the Release build
/// refuses a second load on the same `MzMLFile` after a failed one. This
/// reader decodes each binary array when it closes, where the source decodes
/// a batch of spectra when its data pool is flushed, by default at `</mzML>`
/// (`:1425-1428`, `:1522-1523`), so a document with an undecodable array fails
/// after fewer calls here.
///
/// # Errors
///
/// As [`load_with_options`], plus the errors of the progress calls
/// ([`ProgressLogger::start_progress`] and its siblings).
pub fn load_with_progress(
    path: impl AsRef<Path>,
    scientific: &LoadOptions,
    limits: &ReadOptions,
    logger: &mut ProgressLogger,
) -> Result<MSExperiment> {
    load_reporting(path, scientific, limits, Some(logger))
}

fn load_reporting(
    path: impl AsRef<Path>,
    scientific: &LoadOptions,
    limits: &ReadOptions,
    logger: Option<&mut ProgressLogger>,
) -> Result<MSExperiment> {
    let path = path.as_ref();
    let mut document = crate::metadata::DocumentIdentifier::new();
    let text = path
        .to_str()
        .ok_or_else(|| crate::Error::InvalidValue("mzML filename is not UTF-8".into()))?;
    document.set_loaded_file_path(text)?;
    document.set_loaded_file_type(path)?;
    let input = crate::format::path_io::open(path)?;
    let mut progress = match logger {
        // `File::fileSize(file_)`, the size of the file as stored. The source
        // reports -1 for a file it cannot stat, which an opened file only is
        // when it is removed during the load; the command backend cannot form
        // a rate from that, so no byte count is reported there instead.
        Some(logger) => LoadProgress::new(
            logger,
            std::fs::metadata(path).map_or(0, |metadata| metadata.len()),
        ),
        None => LoadProgress::silent(),
    };
    let mut result = super::read_impl_reporting(
        input,
        limits,
        Some(scientific),
        scientific.scientific.metadata_only,
        &mut progress,
    )?;
    result.settings.document.loaded_file_path = document.loaded_file_path;
    result.settings.document.loaded_file_type = document.loaded_file_type;
    Ok(result)
}

/// Replace the destination only after successful parsing and validation.
pub fn load_into(path: impl AsRef<Path>, destination: &mut MSExperiment) -> Result<()> {
    load_into_with_options(
        path,
        destination,
        &LoadOptions::default(),
        &ReadOptions::default(),
    )
}

/// Checked whole-experiment replacement with explicit loading settings.
pub fn load_into_with_options(
    path: impl AsRef<Path>,
    destination: &mut MSExperiment,
    scientific: &LoadOptions,
    limits: &ReadOptions,
) -> Result<()> {
    let loaded = load_with_options(path, scientific, limits)?;
    *destination = loaded;
    Ok(())
}

/// Atomically store the represented mzML subset. `.gz` and `.bz2` suffixes
/// select outer compression; unknown suffixes are allowed as in XMLFile::save_.
pub fn store(path: impl AsRef<Path>, experiment: &MSExperiment) -> Result<()> {
    store_with_options(path, experiment, &WriteOptions::default())
}

/// Binary-array zlib compression is independent of filename-selected outer
/// compression. Serialization or I/O errors preserve the existing destination.
/// ZIP output is unsupported; no index is generated by this writer.
pub fn store_with_options(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    crate::format::path_io::write(path.as_ref(), |writer| {
        super::write_with_options(writer, experiment, options)
    })
}

/// Store `experiment`, reporting progress to `logger` as source
/// `MzMLFile::store` does through its handler's `writeTo`
/// (`MzMLHandler.cpp:4763-4838`).
///
/// The calls are the handler's: `startProgress(0, spectra + chromatograms,
/// "storing mzML file")` before the document's first byte, `setProgress(n)`
/// before the `n`-th record, spectra first, and `endProgress(bytes written)`
/// after the last byte, so the command backend reports a throughput. The
/// destination is opened before the checks and the first call, as the
/// source's `XMLFile::save_` opens it before `writeTo`. The bytes and every
/// error are those of [`store_with_options`], which runs the same code with
/// the calls going nowhere; its checks precede the start.
///
/// The byte count is that of the document this writer produces, which is not
/// the source's document. With a `.gz` or `.bz2` suffix it is the count before
/// compression, where the source's compressing stream reports no position and
/// passes `-1`.
///
/// A failure after the start leaves the section open, as described at
/// [`load_with_progress`].
///
/// # Errors
///
/// As [`store_with_options`], plus the errors of the progress calls.
pub fn store_with_progress(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &WriteOptions,
    logger: &mut ProgressLogger,
) -> Result<()> {
    let mut progress = ProgressReporter::new(Some(logger));
    crate::format::path_io::write(path.as_ref(), |writer| {
        super::write_with_options_reporting(writer, experiment, options, &mut progress)
    })
}
