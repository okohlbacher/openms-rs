// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A consumer that inserts MS data into an SQLite (sqMass) database.
//!
//! Ports `FORMAT/DATAACCESS/MSDataSqlConsumer.h` and its `.cpp`. See
//! `docs/MS_DATA_CONSUMERS_SUPPORT.md` for the API mapping and the evidence.

use crate::format::mzml_sqlite_handler::{LossPolicy, MzMLSqliteHandler};
use crate::interfaces::MSDataConsumer;
use crate::kernel::{MSChromatogram, MSExperiment, MSSpectrum};
use crate::metadata::ExperimentalSettings;
use crate::{Error, Result};
use std::ops::ControlFlow;
use std::path::{Path, PathBuf};

/// The source's default buffer size: 500 spectra or chromatograms.
pub const DEFAULT_BUFFER_SIZE: usize = 500;
/// The source's default desired mass accuracy for the lossy codec, absolute.
pub const DEFAULT_LINEAR_MASS_ACC: f64 = 1e-4;
/// The most records one buffer reserves room for up front; see
/// [`MSDataSqlConsumer::with_config`].
const MAX_RESERVED_RECORDS: usize = 1 << 16;

/// A data consumer that inserts MS data into an SQLite database.
///
/// Consumes spectra and chromatograms and inserts them into a file-based SQL
/// database (sqMass) using SQLite. As SQLite is highly inefficient when
/// inserting one record at a time, the consumer collects the records in an
/// internal buffer and flushes them together, whenever a buffer reaches the
/// configured size and at the end. It uses [`MzMLSqliteHandler`] to write the
/// batches.
///
/// Consuming a record **clears its peaks and data arrays** in the caller's
/// value (source `s.clear(false)` after buffering a copy); its metadata stays.
/// With full metadata, the cleared record is also kept for the run-level
/// snapshot written by [`finish`](Self::finish).
///
/// # Ending
///
/// The source's destructor flushes the buffers and writes the run-level
/// information (the `RUN` row and, with full metadata, the compressed mzML
/// snapshot in `RUN_EXTRA`), unless [`add_run`](Self::add_run) already wrote
/// one. Here that is [`finish`](Self::finish), which reports failures. A
/// consumer dropped without `finish` does the same on a best-effort basis and
/// logs a warning if it fails, because a drop cannot return an error.
///
/// # Metadata loss
///
/// The handler's [`LossPolicy`] decides what happens to information sqMass
/// does not store; the default refuses it. Without full metadata each batch is
/// checked by the handler as it is written. With full metadata the records'
/// metadata goes into the snapshot at the end, which the handler cannot check
/// while the batches are written, so the batches are written under
/// [`LossPolicy::Source`] and this consumer applies the handler's own
/// full-metadata checks to each record as it arrives instead — earlier, before
/// any of its rows exist.
pub struct MSDataSqlConsumer {
    /// Source `filename_`.
    filename: PathBuf,
    /// Source `handler_`.
    handler: MzMLSqliteHandler,
    /// Source `flush_after_`.
    flush_after: usize,
    /// Source `full_meta_`.
    full_meta: bool,
    /// The policy the caller selected; the handler's own is derived from it.
    loss_policy: LossPolicy,
    /// Source `spectra_`.
    spectra: Vec<MSSpectrum>,
    /// Source `chromatograms_`.
    chromatograms: Vec<MSChromatogram>,
    /// Source `peak_meta_`.
    peak_meta: MSExperiment,
    /// Source `wrote_any_run_`.
    wrote_any_run: bool,
    /// Set once [`MSDataSqlConsumer::finish`] ran, so the drop does nothing.
    finished: bool,
}

impl std::fmt::Debug for MSDataSqlConsumer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MSDataSqlConsumer")
            .field("filename", &self.filename)
            .field("flush_after", &self.flush_after)
            .field("full_meta", &self.full_meta)
            .field("loss_policy", &self.loss_policy)
            .field("spectra", &self.spectra.len())
            .field("chromatograms", &self.chromatograms.len())
            .field("wrote_any_run", &self.wrote_any_run)
            .finish()
    }
}

impl MSDataSqlConsumer {
    /// A consumer with the source's defaults: a buffer of 500 records, full
    /// metadata, no lossy compression and a mass accuracy of `1e-4`.
    ///
    /// # Errors
    ///
    /// As [`with_config`](Self::with_config).
    pub fn new(sql_filename: impl AsRef<Path>, run_id: u64) -> Result<Self> {
        Self::with_config(
            sql_filename,
            run_id,
            DEFAULT_BUFFER_SIZE,
            true,
            false,
            DEFAULT_LINEAR_MASS_ACC,
        )
    }

    /// Opens the SQLite file and writes the tables; source constructor
    /// `MSDataSqlConsumer(sql_filename, run_id, buffer_size, full_meta,
    /// lossy_compression, linear_mass_acc)`.
    ///
    /// # Arguments
    ///
    /// * `sql_filename` — the SQLite database; an existing file is replaced.
    /// * `run_id` — unique identifier which links the sqMass and OSW file;
    ///   bit 63 is cleared, as by the source's handler.
    /// * `buffer_size` — how many spectra (and, separately, chromatograms) are
    ///   collected before they are written; the source's default is 500.
    /// * `full_meta` — whether to write the full metadata into the SQLite
    ///   header (the `RUN_EXTRA` snapshot).
    /// * `lossy_compression` — whether to use lossy (Numpress) compression.
    /// * `linear_mass_acc` — desired mass accuracy for RT or m/z space
    ///   (absolute value).
    ///
    /// The buffers reserve room for at most 65,536 records up front; the
    /// source reserves `buffer_size`, which is a hint only.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `buffer_size` is zero or exceeds the
    /// handler's record limit, or `linear_mass_acc` is not finite: the source
    /// passes an `int` batch size to its handler unchecked. [`Error::Io`] and
    /// the handler's errors when the tables cannot be created.
    pub fn with_config(
        sql_filename: impl AsRef<Path>,
        run_id: u64,
        buffer_size: usize,
        full_meta: bool,
        lossy_compression: bool,
        linear_mass_acc: f64,
    ) -> Result<Self> {
        let filename = sql_filename.as_ref().to_path_buf();
        let mut handler = MzMLSqliteHandler::new(&filename, run_id);
        handler.set_config(full_meta, lossy_compression, linear_mass_acc, buffer_size)?;
        let loss_policy = LossPolicy::default();
        handler.set_loss_policy(batch_policy(full_meta, loss_policy));
        handler.create_tables()?;
        let reserve = buffer_size.min(MAX_RESERVED_RECORDS);
        Ok(Self {
            filename,
            handler,
            flush_after: buffer_size,
            full_meta,
            loss_policy,
            spectra: Vec::with_capacity(reserve),
            chromatograms: Vec::with_capacity(reserve),
            peak_meta: MSExperiment::default(),
            wrote_any_run: false,
            finished: false,
        })
    }

    /// Selects what happens to information sqMass does not store; see
    /// "Metadata loss" above.
    ///
    /// No source counterpart: the source's handler discards it silently,
    /// which [`LossPolicy::Source`] reproduces.
    pub fn set_loss_policy(&mut self, policy: LossPolicy) {
        self.loss_policy = policy;
        self.handler
            .set_loss_policy(batch_policy(self.full_meta, policy));
    }

    /// The loss policy selected with [`set_loss_policy`](Self::set_loss_policy).
    pub fn loss_policy(&self) -> LossPolicy {
        self.loss_policy
    }

    /// Writes the buffered records; source `flush`.
    ///
    /// Afterwards no record is held in the buffers, and the consumer still
    /// accepts new ones. Spectra are written before chromatograms, each as one
    /// atomic handler write.
    ///
    /// # Errors
    ///
    /// The handler's errors. A failed write leaves its buffer as it was, so a
    /// later flush retries it; a spectrum batch already written stays written
    /// when the chromatogram batch then fails.
    pub fn flush(&mut self) -> Result<()> {
        if !self.spectra.is_empty() {
            self.handler.write_spectra(&self.spectra)?;
            self.spectra.clear();
        }
        if !self.chromatograms.is_empty() {
            self.handler.write_chromatograms(&self.chromatograms)?;
            self.chromatograms.clear();
        }
        Ok(())
    }

    /// Adds a `RUN` entry (identifier and file name) to the sqMass file;
    /// source `addRun`.
    ///
    /// Sets the handler's run identifier to `run_id` for the records written
    /// from now on, and writes the run-level information of an experiment
    /// that holds nothing but `filename` as its loaded file path — with full
    /// metadata, a snapshot describing no records. Afterwards
    /// [`finish`](Self::finish) writes no run-level information of its own,
    /// as in the source.
    ///
    /// # Errors
    ///
    /// The handler's errors.
    pub fn add_run(&mut self, filename: &str, run_id: u64) -> Result<()> {
        self.handler.set_run_id(run_id);
        let mut meta = MSExperiment::default();
        meta.settings.document.loaded_file_path = filename.to_string();
        self.handler
            .write_run_level_information(&meta, self.full_meta)?;
        self.wrote_any_run = true;
        Ok(())
    }

    /// Changes the run identifier used for subsequent spectrum and
    /// chromatogram writes; source `setRunId`. Bit 63 is cleared.
    ///
    /// Records already buffered are written under the identifier in force when
    /// they are flushed, as in the source.
    pub fn set_run_id(&mut self, run_id: u64) {
        self.handler.set_run_id(run_id);
    }

    /// Flushes the buffers and writes the run-level information, unless a run
    /// was already added; the source's destructor.
    ///
    /// The run-level information is that of the consumer's own metadata
    /// experiment: default experimental settings — the source ignores
    /// `setExperimentalSettings` — with the **sqMass file name** as the loaded
    /// file path, and, with full metadata, the metadata of every record
    /// consumed.
    ///
    /// # Errors
    ///
    /// The handler's errors. The records flushed before a failure stay
    /// written; a failed run-level write leaves the file without its `RUN`
    /// row, which the source's destructor would leave as well.
    pub fn finish(mut self) -> Result<()> {
        self.finished = true;
        self.end()
    }

    fn end(&mut self) -> Result<()> {
        self.flush()?;
        if !self.wrote_any_run {
            self.peak_meta.settings.document.loaded_file_path =
                self.filename.to_string_lossy().into_owned();
            self.handler
                .write_run_level_information(&self.peak_meta, self.full_meta)?;
            self.wrote_any_run = true;
        }
        Ok(())
    }

    /// The SQLite file name; the source's protected `filename_`.
    pub fn filename(&self) -> &Path {
        &self.filename
    }

    /// The handler writing the batches; the source's protected `handler_`.
    pub fn handler(&self) -> &MzMLSqliteHandler {
        &self.handler
    }

    /// Records per buffer before a flush; the source's protected `flush_after_`.
    pub fn flush_after(&self) -> usize {
        self.flush_after
    }

    /// Whether the full metadata is written; the source's protected `full_meta_`.
    pub fn full_meta(&self) -> bool {
        self.full_meta
    }

    /// Spectra waiting in the buffer; the source's protected `spectra_`.
    pub fn buffered_spectra(&self) -> &[MSSpectrum] {
        &self.spectra
    }

    /// Chromatograms waiting in the buffer; the source's protected `chromatograms_`.
    pub fn buffered_chromatograms(&self) -> &[MSChromatogram] {
        &self.chromatograms
    }

    /// The metadata experiment for the run-level snapshot; the source's
    /// protected `peak_meta_`.
    pub fn peak_meta(&self) -> &MSExperiment {
        &self.peak_meta
    }

    /// Whether run-level information was written; the source's protected
    /// `wrote_any_run_`.
    pub fn wrote_any_run(&self) -> bool {
        self.wrote_any_run
    }

    /// The handler's full-metadata checks for one record, applied when it
    /// arrives. The handler refuses these at the end, after every row exists.
    fn check_full_meta(
        &self,
        native_id: &str,
        has_auxiliary: bool,
        has_identifications: bool,
    ) -> Result<()> {
        if !self.full_meta {
            return Ok(());
        }
        if has_auxiliary && self.loss_policy == LossPolicy::Reject {
            return Err(Error::Unsupported(
                "sqMass source format discards auxiliary arrays; select Source loss policy".into(),
            ));
        }
        if has_identifications {
            return Err(Error::Unsupported(
                "mzML snapshot cannot preserve peptide identifications".into(),
            ));
        }
        if native_id.is_empty() {
            return Err(Error::InvalidValue(
                "sqMass: full metadata requires nonempty native IDs".into(),
            ));
        }
        Ok(())
    }
}

/// The policy the handler applies to a batch: under full metadata the snapshot
/// keeps what the SQL rows drop, and the consumer checks the rest itself.
fn batch_policy(full_meta: bool, policy: LossPolicy) -> LossPolicy {
    if full_meta { LossPolicy::Source } else { policy }
}

impl Drop for MSDataSqlConsumer {
    fn drop(&mut self) {
        if self.finished {
            return;
        }
        self.finished = true;
        if let Err(error) = self.end() {
            // The source's destructor does the same work and has nowhere to
            // report a failure either; a drop must not panic.
            let _ = crate::openms_log_warn!(
                "MSDataSqlConsumer: ending {} failed: {error}",
                self.filename.display()
            );
        }
    }
}

impl MSDataConsumer for MSDataSqlConsumer {
    /// Ignored, as in the source.
    fn set_expected_size(&mut self, _spectra: usize, _chromatograms: usize) -> Result<()> {
        Ok(())
    }

    /// Ignored, as in the source: the run-level snapshot carries default
    /// settings. See [`MSDataSqlConsumer::finish`].
    fn set_experimental_settings(&mut self, _settings: &ExperimentalSettings) -> Result<()> {
        Ok(())
    }

    /// Buffers a copy of the spectrum and clears the caller's peaks and data
    /// arrays; writes the buffer when it is full. Source `consumeSpectrum`.
    ///
    /// # Errors
    ///
    /// With full metadata: [`Error::Unsupported`] for auxiliary data arrays
    /// under [`LossPolicy::Reject`] and for peptide identifications, and
    /// [`Error::InvalidValue`] for an empty native identifier — the handler's
    /// own full-metadata refusals, raised before the record is taken, so the
    /// spectrum and the consumer are unchanged. Otherwise the handler's errors
    /// from a flush, which leave the record buffered (and the caller's
    /// spectrum cleared), as the source's throwing flush does.
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        self.check_full_meta(
            &spectrum.native_id,
            !spectrum.float_data_arrays.is_empty()
                || !spectrum.integer_data_arrays.is_empty()
                || !spectrum.string_data_arrays.is_empty(),
            !spectrum.peptide_identifications.is_empty(),
        )?;
        self.spectra.push(spectrum.clone());
        spectrum.clear(false);
        if self.full_meta {
            self.peak_meta.spectra.push(spectrum.clone());
        }
        if self.spectra.len() >= self.flush_after {
            self.flush()?;
        }
        Ok(ControlFlow::Continue(()))
    }

    /// Buffers a copy of the chromatogram and clears the caller's peaks and
    /// data arrays; writes the buffer when it is full. Source
    /// `consumeChromatogram`.
    ///
    /// # Errors
    ///
    /// As [`consume_spectrum`](Self::consume_spectrum), without the peptide
    /// identification check, which chromatograms cannot carry.
    fn consume_chromatogram(
        &mut self,
        chromatogram: &mut MSChromatogram,
    ) -> Result<ControlFlow<()>> {
        self.check_full_meta(
            &chromatogram.native_id,
            !chromatogram.float_data_arrays.is_empty()
                || !chromatogram.integer_data_arrays.is_empty()
                || !chromatogram.string_data_arrays.is_empty(),
            false,
        )?;
        self.chromatograms.push(chromatogram.clone());
        chromatogram.clear(false);
        if self.full_meta {
            self.peak_meta.chromatograms.push(chromatogram.clone());
        }
        if self.chromatograms.len() >= self.flush_after {
            self.flush()?;
        }
        Ok(ControlFlow::Continue(()))
    }
}
