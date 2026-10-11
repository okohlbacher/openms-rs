// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A consumer that writes MS data into the cached-mzML binary format.
//!
//! Ports `FORMAT/DATAACCESS/MSDataCachedConsumer.h` and its `.cpp`. The file
//! it writes is the binary half of a cached mzML pair, laid out as
//! [`crate::format::cached_mzml_handler`] describes and read back by
//! [`CachedMzMLHandler`](crate::format::cached_mzml_handler::CachedMzMLHandler)
//! and [`crate::format::cached_mzml`]. See
//! `docs/MS_DATA_CONSUMERS_SUPPORT.md` for the API mapping and the evidence.

use crate::format::cached_mzml_handler::{CACHED_MZML_FILE_IDENTIFIER, CachedMzMLHandler};
use crate::interfaces::MSDataConsumer;
use crate::kernel::{MSChromatogram, MSSpectrum};
use crate::metadata::ExperimentalSettings;
use crate::{Error, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::ops::ControlFlow;
use std::path::{Path, PathBuf};

/// The source's message when a spectrum follows a chromatogram.
const SPECTRUM_AFTER_CHROMATOGRAM: &str = "Cannot write spectra after writing chromatograms.";

/// Cached writing consumer of MS data (source `MSDataCachedConsumer`).
///
/// Each spectrum and chromatogram is written to disk as it arrives, with the
/// record writers of [`CachedMzMLHandler`]: the file identifier when the
/// consumer is created, then every spectrum, then every chromatogram, and the
/// two record counts at the end. Only the peaks, the float and integer data
/// arrays, and a spectrum's MS level and retention time are cached; the
/// metadata belongs in the companion mzML file (see
/// [`CachedMzMLHandler::write_metadata`]).
///
/// With `clear_data` (the source's default), consuming a record clears its
/// peaks and data arrays in the caller's value, keeping only its metadata
/// (`clear(false)`); the float and integer data arrays are cleared too, as
/// they are written to disk as well.
///
/// The source class derives from `Internal::CachedMzMLHandler`; the handler
/// it inherits is [`handler`](Self::handler) and
/// [`handler_mut`](Self::handler_mut).
///
/// # Ending
///
/// The source's destructor writes the two record counts and closes the file.
/// Here that is [`finish`](Self::finish), which reports failures. A consumer
/// dropped without `finish` does the same on a best-effort basis and logs a
/// warning if it fails, because a drop cannot return an error. Until the
/// counts are written the file is not a valid cache.
#[derive(Debug)]
pub struct MSDataCachedConsumer {
    /// The inherited `Internal::CachedMzMLHandler`.
    handler: CachedMzMLHandler,
    /// The output path, for messages.
    filename: PathBuf,
    /// Source `ofs_`; `None` once the counts were written.
    ofs: Option<BufWriter<File>>,
    /// Source `clearData_`.
    clear_data: bool,
    /// Source `spectra_written_`.
    spectra_written: u64,
    /// Source `chromatograms_written_`.
    chromatograms_written: u64,
}

impl MSDataCachedConsumer {
    /// A consumer that clears the data of every record it writes; the
    /// source's constructor with its default `clearData = true`.
    ///
    /// # Errors
    ///
    /// As [`with_clear_data`](Self::with_clear_data).
    pub fn new(filename: impl AsRef<Path>) -> Result<Self> {
        Self::with_clear_data(filename, true)
    }

    /// Opens the output file and writes the header; source constructor
    /// `MSDataCachedConsumer(filename, clearData)`.
    ///
    /// # Arguments
    ///
    /// * `filename` — the output file to which data is written; an existing
    ///   file is replaced.
    /// * `clear_data` — whether to clear the spectral and chromatogram data
    ///   after writing (only keep the metadata). Clearing data from spectra and
    ///   chromatograms also clears the float and integer data arrays
    ///   associated with them, as these are written to disk as well.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be created or the identifier cannot
    /// be written. The source opens the stream unchecked and writes nothing
    /// to a stream that failed to open.
    pub fn with_clear_data(filename: impl AsRef<Path>, clear_data: bool) -> Result<Self> {
        let filename = filename.as_ref().to_path_buf();
        let mut ofs = BufWriter::new(File::create(&filename)?);
        ofs.write_all(&CACHED_MZML_FILE_IDENTIFIER.to_le_bytes())?;
        Ok(Self {
            handler: CachedMzMLHandler::new(),
            filename,
            ofs: Some(ofs),
            clear_data,
            spectra_written: 0,
            chromatograms_written: 0,
        })
    }

    /// The inherited cache handler (the source's `Internal::CachedMzMLHandler`
    /// base), for its reading and indexing members.
    pub fn handler(&self) -> &CachedMzMLHandler {
        &self.handler
    }

    /// The inherited cache handler, mutably (for
    /// [`CachedMzMLHandler::create_memdump_index`]).
    pub fn handler_mut(&mut self) -> &mut CachedMzMLHandler {
        &mut self.handler
    }

    /// The output file; the path the source's `ofs_` was opened on.
    pub fn filename(&self) -> &Path {
        &self.filename
    }

    /// Whether consumed records lose their data; source `clearData_`.
    pub fn clear_data(&self) -> bool {
        self.clear_data
    }

    /// Spectra written so far; source `spectra_written_`.
    pub fn spectra_written(&self) -> u64 {
        self.spectra_written
    }

    /// Chromatograms written so far; source `chromatograms_written_`.
    pub fn chromatograms_written(&self) -> u64 {
        self.chromatograms_written
    }

    /// Closes the output file and writes the footer, the numbers of spectra
    /// and chromatograms written; the source's destructor.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the counts cannot be written or the file cannot be
    /// flushed. The source's destructor ignores the stream state.
    pub fn finish(mut self) -> Result<()> {
        self.end()
    }

    /// The destructor's work; a no-op once done.
    fn end(&mut self) -> Result<()> {
        let Some(mut ofs) = self.ofs.take() else {
            return Ok(());
        };
        ofs.write_all(&self.spectra_written.to_le_bytes())?;
        ofs.write_all(&self.chromatograms_written.to_le_bytes())?;
        ofs.flush()?;
        Ok(())
    }

    /// The open stream, which exists until the consumer is finished.
    fn stream(&mut self) -> Result<&mut BufWriter<File>> {
        self.ofs.as_mut().ok_or_else(|| {
            Error::InvalidValue("MSDataCachedConsumer: the output file is already closed".into())
        })
    }
}

impl Drop for MSDataCachedConsumer {
    fn drop(&mut self) {
        if let Err(error) = self.end() {
            // The source's destructor does the same work and has nowhere to
            // report a failure either; a drop must not panic.
            let _ = crate::openms_log_warn!(
                "MSDataCachedConsumer: ending {} failed: {error}",
                self.filename.display()
            );
        }
    }
}

impl MSDataConsumer for MSDataCachedConsumer {
    /// Ignored, as in the source.
    fn set_expected_size(&mut self, _spectra: usize, _chromatograms: usize) -> Result<()> {
        Ok(())
    }

    /// Ignored, as in the source: the cache holds no experimental settings.
    fn set_experimental_settings(&mut self, _settings: &ExperimentalSettings) -> Result<()> {
        Ok(())
    }

    /// Writes a spectrum to the output file (source `consumeSpectrum`) and,
    /// with `clear_data`, clears its peaks and data arrays.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] with the source's message, "Cannot write
    /// spectra after writing chromatograms.", once a chromatogram was written
    /// (the source's `Exception::IllegalArgument`); the spectrum and the file
    /// are unchanged. [`Error::Io`] when the write fails, after which the file
    /// holds a partial record.
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        if self.chromatograms_written > 0 {
            return Err(Error::InvalidValue(SPECTRUM_AFTER_CHROMATOGRAM.into()));
        }
        CachedMzMLHandler::write_spectrum(spectrum, self.stream()?)?;
        self.spectra_written += 1;
        if self.clear_data {
            spectrum.clear(false);
        }
        Ok(ControlFlow::Continue(()))
    }

    /// Writes a chromatogram to the output file (source
    /// `consumeChromatogram`) and, with `clear_data`, clears its peaks and
    /// data arrays.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the write fails, after which the file holds a
    /// partial record.
    fn consume_chromatogram(
        &mut self,
        chromatogram: &mut MSChromatogram,
    ) -> Result<ControlFlow<()>> {
        CachedMzMLHandler::write_chromatogram(chromatogram, self.stream()?)?;
        self.chromatograms_written += 1;
        if self.clear_data {
            chromatogram.clear(false);
        }
        Ok(ControlFlow::Continue(()))
    }
}
