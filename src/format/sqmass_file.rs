// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Mass-spectrometry data in the on-disk sqMass (SQLite) format, ported from
//! `FORMAT/SqMassFile.h`.
//!
//! Spectra and chromatograms are persisted in a single SQLite database, so
//! that large collections can be queried, filtered or streamed on demand
//! without loading the whole experiment into memory. For spectra and
//! chromatograms that carry precursor information, the `peptide_sequence`
//! meta value of the first precursor is round-tripped with the binary data.
//! The storage itself is [`MzMLSqliteHandler`]; this adapter applies one
//! [`SqMassConfig`] to it, as the source does. See
//! `docs/SQMASS_FILE_SUPPORT.md`.

use crate::Result;
use crate::format::mzml_sqlite_handler::{LossPolicy, MzMLSqliteHandler};
use crate::interfaces::MSDataConsumer;
use crate::kernel::MSExperiment;
use std::ops::ControlFlow;
use std::path::Path;

/// Records read per batch while streaming (the source's `batch_size`).
const BATCH_SIZE: usize = 500;

/// Compression and metadata options used when reading or writing an sqMass
/// file (source `SqMassFile::SqMassConfig`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SqMassConfig {
    /// Write full metadata (verbose, lossless): a compressed mzML snapshot of
    /// the experiment's descriptive metadata in the `RUN_EXTRA` table.
    pub write_full_meta: bool,
    /// Use lossy Numpress compression for the binary data arrays.
    pub use_lossy_numpress: bool,
    /// Target mass accuracy for Numpress linear encoding, in Th; `-1` has no
    /// effect (use e.g. `0.0001` for 0.2 ppm at 500 m/z).
    pub linear_fp_mass_acc: f64,
    /// Native: whether a store may discard the metadata that neither the SQL
    /// tables nor (without `write_full_meta`) the snapshot keep. The source
    /// discards it silently; the native default,
    /// [`LossPolicy::Reject`], refuses such an experiment, and
    /// [`LossPolicy::Source`] selects the source behaviour.
    pub loss_policy: LossPolicy,
}

impl Default for SqMassConfig {
    /// The source defaults (full metadata, no lossy compression, accuracy
    /// `-1`) and the native [`LossPolicy::Reject`].
    fn default() -> Self {
        Self {
            write_full_meta: true,
            use_lossy_numpress: false,
            linear_fp_mass_acc: -1.0,
            loss_policy: LossPolicy::Reject,
        }
    }
}

/// The source's `MSExperiment::getSqlRunID` over this crate's two places for
/// the ID (see [`SqMassFile::store`]).
fn sql_run_id(map: &MSExperiment) -> u64 {
    if map.sql_run_id != 0 {
        return map.sql_run_id;
    }
    match map
        .settings
        .metadata
        .get("sqMassRunID")
        .map(|value| value.data())
    {
        Some(crate::metadata::MetaValueData::Integer(id)) => u64::try_from(*id).unwrap_or(0),
        _ => 0,
    }
}

/// Reads and writes sqMass files (source `SqMassFile`).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SqMassFile {
    config: SqMassConfig,
}

impl SqMassFile {
    /// An adapter with the default [`SqMassConfig`] (source default
    /// constructor).
    pub fn new() -> Self {
        Self::default()
    }

    /// Replace the compression and metadata configuration used by subsequent
    /// [`Self::load`], [`Self::store`] and [`Self::transform`] calls (source
    /// `setConfig`).
    pub fn set_config(&mut self, config: SqMassConfig) {
        self.config = config;
    }

    /// The current configuration.
    pub fn config(&self) -> SqMassConfig {
        self.config
    }

    fn handler(&self, filename: &Path, run_id: u64) -> Result<MzMLSqliteHandler> {
        let mut handler = MzMLSqliteHandler::new(filename, run_id);
        handler.set_config(
            self.config.write_full_meta,
            self.config.use_lossy_numpress,
            self.config.linear_fp_mass_acc,
            BATCH_SIZE,
        )?;
        handler.set_loss_policy(self.config.loss_policy);
        Ok(handler)
    }

    /// Read all spectra, chromatograms and experimental metadata of the sqMass
    /// file `filename` (source `load`), applying the current configuration.
    ///
    /// # Errors
    ///
    /// The errors of [`MzMLSqliteHandler::read_experiment`]: malformed or
    /// inconsistent rows (a mismatched native ID, a missing data type, more
    /// than one run), where the source throws `Exception::IllegalArgument`; a
    /// binary buffer whose size does not match its data type
    /// (`Exception::ConversionError`); and an unsupported configuration in the
    /// `RUN` table (`Exception::SqlOperationFailed`).
    pub fn load(&self, filename: impl AsRef<Path>) -> Result<MSExperiment> {
        self.handler(filename.as_ref(), 0)?.read_experiment(false)
    }

    /// Write the spectra, chromatograms and experimental metadata of `map` to
    /// `filename`, creating the file and its tables (source `store`).
    ///
    /// The `RUN` table's ID is the source's `MSExperiment::getSqlRunID`. The
    /// source keeps that ID as the experiment's `sqMassRunID` meta value, which
    /// is what an mzML file carries; this crate reads it into
    /// `MSExperiment::sql_run_id` from sqMass but leaves it a meta value when
    /// reading mzML. So the ID is `sql_run_id` when that is non-zero, else a
    /// non-negative integer `sqMassRunID` meta value, else 0.
    ///
    /// # Errors
    ///
    /// The errors of [`MzMLSqliteHandler::create_tables`] and
    /// [`MzMLSqliteHandler::write_experiment`]: SQL failures, where the source
    /// throws `Exception::IllegalArgument` or `Exception::SqlOperationFailed`,
    /// and, under [`LossPolicy::Reject`], metadata the configured storage would
    /// discard.
    pub fn store(&self, filename: impl AsRef<Path>, map: &MSExperiment) -> Result<()> {
        let mut handler = self.handler(filename.as_ref(), sql_run_id(map))?;
        handler.create_tables()?;
        handler.write_experiment(map)
    }

    /// Stream the spectra and chromatograms of `filename_in` through
    /// `consumer` (source `transform`).
    ///
    /// The consumer first learns the expected counts and the experimental
    /// settings (read without peak data), then receives every spectrum and
    /// every chromatogram in input order, read in batches of 500, so the full
    /// experiment is never held in memory. `skip_full_count` and
    /// `skip_first_pass` are ignored, as in the source, which keeps them for
    /// API compatibility. A consumer that returns [`ControlFlow::Break`] ends
    /// the stream early, which the source's consumers cannot request.
    ///
    /// The source also issues an empty last batch when a count is a positive
    /// multiple of 500; its handler then selects every record for the empty ID
    /// list and throws `Exception::IllegalArgument`. Empty batches are skipped
    /// here.
    ///
    /// # Errors
    ///
    /// As [`Self::load`], and the consumer's errors.
    pub fn transform(
        &self,
        filename_in: impl AsRef<Path>,
        consumer: &mut dyn MSDataConsumer,
        _skip_full_count: bool,
        _skip_first_pass: bool,
    ) -> Result<()> {
        let handler = self.handler(filename_in.as_ref(), 0)?;
        let spectra = handler.nr_spectra()?;
        let chromatograms = handler.nr_chromatograms()?;
        consumer.set_expected_size(spectra, chromatograms)?;
        let settings = handler.read_experiment(true)?;
        consumer.set_experimental_settings(&settings.settings)?;

        for start in (0..spectra).step_by(BATCH_SIZE) {
            let end = (start + BATCH_SIZE).min(spectra);
            let indices: Vec<i64> = (start as i64..end as i64).collect();
            for mut spectrum in handler.read_spectra(&indices, false)? {
                if consumer.consume_spectrum(&mut spectrum)?.is_break() {
                    return Ok(());
                }
            }
        }
        for start in (0..chromatograms).step_by(BATCH_SIZE) {
            let end = (start + BATCH_SIZE).min(chromatograms);
            let indices: Vec<i64> = (start as i64..end as i64).collect();
            for mut chromatogram in handler.read_chromatograms(&indices, false)? {
                if let ControlFlow::Break(()) = consumer.consume_chromatogram(&mut chromatogram)? {
                    return Ok(());
                }
            }
        }
        Ok(())
    }
}
