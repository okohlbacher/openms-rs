// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded sqMass storage for the public `Internal::MzMLSqliteHandler` API.
//!
//! Reads use SQL record identities and a single read transaction. Each public
//! write is atomic; writer counters advance only after commit. Lossless arrays
//! are explicitly little-endian f64; compression codes 1, 5 and 6 are supported.
//! See `docs/MZML_SQLITE_HANDLER_SUPPORT.md` for metadata and source differences.

use super::numpress::NumpressLimits;
use super::numpress_coder::{
    self as codec, NumpressCoderLimits, NumpressCompression, NumpressConfig, NumpressEncodeStatus,
    Work,
};
use super::sqlite_connector::{SqlOpenMode, SqliteConnector};
use crate::kernel::{
    ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, Precursor,
};
use crate::metadata::{ActivationMethod, MetaValue, Polarity, Product};
use crate::{Error, Result};
use rusqlite::{Connection, OptionalExtension, Row, params};
use std::io::{Cursor, Write};
use std::path::{Path, PathBuf};

/// Whether unsupported metadata may be discarded according to the C++ writer.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum LossPolicy {
    /// Refuse information that neither SQL nor the requested snapshot preserves.
    #[default]
    Reject,
    /// Explicitly permit the documented source metadata and auxiliary-array loss.
    Source,
}

/// Limits for one public operation. These bound native materialization and fixed
/// schema scans; they are not SQLite virtual-machine time or memory limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HandlerLimits {
    /// Maximum combined spectra and chromatograms (and requested IDs).
    pub max_records: usize,
    /// Maximum values in one coordinate or intensity array.
    pub max_values_per_array: usize,
    /// Combined decoded coordinate and intensity values in an operation.
    pub max_total_values: usize,
    /// Maximum compressed or uncompressed array BLOB bytes.
    pub max_blob_bytes: usize,
    /// Shared logical allocation budget, including temporary codec state.
    pub max_total_bytes: usize,
    /// Maximum compressed or expanded RUN_EXTRA mzML bytes.
    pub max_metadata_bytes: usize,
    /// Maximum physical input and logical database bytes, also checked before commit.
    pub max_database_bytes: u64,
    /// Shared metadata/codec traversal budget.
    pub max_work: usize,
}
impl Default for HandlerLimits {
    fn default() -> Self {
        Self {
            max_records: 1_000_000,
            max_values_per_array: 10_000_000,
            max_total_values: 20_000_000,
            max_blob_bytes: 64 * 1024 * 1024,
            max_total_bytes: 512 * 1024 * 1024,
            max_metadata_bytes: 64 * 1024 * 1024,
            max_database_bytes: 2 * 1024 * 1024 * 1024,
            max_work: 500_000_000,
        }
    }
}

/// Source configuration plus the explicit native metadata-loss policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HandlerConfig {
    /// Store/read a compressed mzML descriptive snapshot in RUN_EXTRA.
    pub write_full_meta: bool,
    /// Use linear Numpress coordinates and SLOF intensities instead of doubles.
    pub use_lossy_compression: bool,
    /// Absolute linear m/z precision in Th; RT uses the source's 0.05 seconds.
    pub linear_abs_mass_acc: f64,
    /// Maximum records visited per write chunk; one transaction spans all chunks.
    pub sql_batch_size: usize,
    /// Policy for metadata outside the selected storage representation.
    pub loss_policy: LossPolicy,
}
impl Default for HandlerConfig {
    fn default() -> Self {
        Self {
            write_full_meta: true,
            use_lossy_compression: true,
            linear_abs_mass_acc: 0.0001,
            sql_batch_size: 500,
            loss_policy: LossPolicy::Reject,
        }
    }
}

/// Path-owning sqMass handler. Prefer the higher-level SqMassFile API when available.
///
/// Construction does not open or create a file. New handlers start writer IDs at
/// zero, as the source does; they cannot append to an already populated database.
/// Reads accept arbitrary nonnegative SQL IDs, including gaps.
#[derive(Debug)]
pub struct MzMLSqliteHandler {
    filename: PathBuf,
    run_id: i64,
    spectrum_id: i64,
    chromatogram_id: i64,
    config: HandlerConfig,
    /// Resource ceilings checked before allocation and while decoding.
    pub limits: HandlerLimits,
}

impl MzMLSqliteHandler {
    /// Construct a handler, masking run-ID bit 63 as in the source. Defaults are
    /// full metadata, lossy compression, 0.0001 Th and an initialized batch of 500.
    pub fn new(filename: impl AsRef<Path>, run_id: u64) -> Self {
        Self {
            filename: filename.as_ref().to_path_buf(),
            run_id: (run_id & i64::MAX as u64) as i64,
            spectrum_id: 0,
            chromatogram_id: 0,
            config: HandlerConfig::default(),
            limits: HandlerLimits::default(),
        }
    }
    /// Current configuration; changing it requires a validated setter.
    pub fn config(&self) -> HandlerConfig {
        self.config
    }
    /// Set the source's four configuration fields atomically. Accuracy must be
    /// finite; nonpositive accuracy selects ordinary estimated fixed point.
    /// Batch size must be positive and within the record limit.
    pub fn set_config(
        &mut self,
        full_meta: bool,
        lossy: bool,
        accuracy: f64,
        batch: usize,
    ) -> Result<()> {
        if !accuracy.is_finite() || batch == 0 || batch > self.limits.max_records {
            return Err(invalid("invalid accuracy or SQL batch size"));
        }
        self.config.write_full_meta = full_meta;
        self.config.use_lossy_compression = lossy;
        self.config.linear_abs_mass_acc = accuracy;
        self.config.sql_batch_size = batch;
        Ok(())
    }
    /// Explicitly select whether the source's metadata loss is permitted.
    pub fn set_loss_policy(&mut self, policy: LossPolicy) {
        self.config.loss_policy = policy;
    }
    /// Change the masked ID used by subsequent writes; reads use the RUN table.
    pub fn set_run_id(&mut self, run_id: u64) {
        self.run_id = (run_id & i64::MAX as u64) as i64;
    }
    /// Read exactly one nonnegative RUN ID. Zero or multiple runs are errors.
    pub fn run_id(&self) -> Result<u64> {
        let connector = self.open(false)?;
        let transaction = connector
            .connection()
            .unchecked_transaction()
            .map_err(sql_error)?;
        self.schema(&transaction)?;
        let result = read_run(&transaction, &mut Budget::new(self.limits))?.0 as u64;
        transaction.commit().map_err(sql_error)?;
        Ok(result)
    }
    /// Number of spectrum records, checked against the configured limits.
    pub fn nr_spectra(&self) -> Result<usize> {
        self.count("SPECTRUM")
    }
    /// Number of chromatogram records, checked against the configured limits.
    pub fn nr_chromatograms(&self) -> Result<usize> {
        self.count("CHROMATOGRAM")
    }
    fn count(&self, table: &str) -> Result<usize> {
        let connector = self.open(false)?;
        let transaction = connector
            .connection()
            .unchecked_transaction()
            .map_err(sql_error)?;
        self.schema(&transaction)?;
        let n = count(&transaction, table)?;
        transaction.commit().map_err(sql_error)?;
        Ok(n)
    }
    /// Read a whole experiment in one SQLite snapshot. RUN_EXTRA metadata, when
    /// enabled and present, must agree with SQL record counts and native IDs.
    /// `meta_only` omits all primary peak arrays. Errors return no partial output.
    pub fn read_experiment(&self, meta_only: bool) -> Result<MSExperiment> {
        let connector = self.open(false)?;
        let tx = connector
            .connection()
            .unchecked_transaction()
            .map_err(sql_error)?;
        self.schema(&tx)?;
        let mut budget = Budget::new(self.limits);
        let (run_id, filename) = read_run(&tx, &mut budget)?;
        for table in ["SPECTRUM", "CHROMATOGRAM"] {
            let mismatch: bool = tx.query_row(&format!("SELECT EXISTS(SELECT 1 FROM {table} WHERE RUN_ID IS NULL OR typeof(RUN_ID)!='integer' OR RUN_ID!=?1)"), [run_id], |r| r.get(0)).map_err(sql_error)?;
            if mismatch {
                return Err(invalid("record RUN_ID does not match the experiment RUN"));
            }
        }
        let spectra_ids = ids(&tx, "SPECTRUM", &mut budget)?;
        let chrom_ids = ids(&tx, "CHROMATOGRAM", &mut budget)?;
        let mut spectra = self.spectra(&tx, &spectra_ids, true, &mut budget)?;
        let mut chromatograms = self.chromatograms(&tx, &chrom_ids, true, &mut budget)?;
        let mut experiment = MSExperiment::default();
        if self.config.write_full_meta {
            if let Some(snapshot) = self.snapshot(&tx, run_id, &mut budget)? {
                validate_snapshot(&snapshot, &spectra, &chromatograms)?;
                experiment = snapshot;
                spectra = std::mem::take(&mut experiment.spectra);
                chromatograms = std::mem::take(&mut experiment.chromatograms);
            }
        }
        if !meta_only {
            for (id, spectrum) in spectra_ids.iter().zip(&mut spectra) {
                let (coordinate, intensity) = self.arrays(&tx, *id, false, &mut budget)?;
                spectrum.peaks = coordinate
                    .into_iter()
                    .zip(intensity)
                    .map(|(mz, intensity)| Peak1D { mz, intensity })
                    .collect();
            }
            for (id, chrom) in chrom_ids.iter().zip(&mut chromatograms) {
                let (coordinate, intensity) = self.arrays(&tx, *id, true, &mut budget)?;
                chrom.peaks = coordinate
                    .into_iter()
                    .zip(intensity)
                    .map(|(rt, intensity)| ChromatogramPeak { rt, intensity })
                    .collect();
            }
        }
        experiment.spectra = spectra;
        experiment.chromatograms = chromatograms;
        experiment.sql_run_id = run_id as u64;
        experiment.settings.metadata.remove("sqMassRunID");
        experiment.settings.document.loaded_file_path = filename;
        tx.commit().map_err(sql_error)?;
        Ok(experiment)
    }
    /// Read a nonempty set of unique nonnegative SQL IDs, in increasing ID order.
    /// Selected reads reconstruct SQL metadata only; they never consult RUN_EXTRA.
    pub fn read_spectra(&self, indices: &[i64], meta_only: bool) -> Result<Vec<MSSpectrum>> {
        let mut budget = Budget::new(self.limits);
        let selected = selection(indices, false, &mut budget)?;
        let connector = self.open(false)?;
        let tx = connector
            .connection()
            .unchecked_transaction()
            .map_err(sql_error)?;
        self.schema(&tx)?;
        let result = self.spectra(&tx, &selected, meta_only, &mut budget)?;
        tx.commit().map_err(sql_error)?;
        Ok(result)
    }
    /// Chromatogram counterpart of [`Self::read_spectra`], with the same ID and
    /// metadata-only rules. Missing IDs and malformed array pairs are errors.
    pub fn read_chromatograms(
        &self,
        indices: &[i64],
        meta_only: bool,
    ) -> Result<Vec<MSChromatogram>> {
        let mut budget = Budget::new(self.limits);
        let selected = selection(indices, false, &mut budget)?;
        let connector = self.open(false)?;
        let tx = connector
            .connection()
            .unchecked_transaction()
            .map_err(sql_error)?;
        self.schema(&tx)?;
        let result = self.chromatograms(&tx, &selected, meta_only, &mut budget)?;
        tx.commit().map_err(sql_error)?;
        Ok(result)
    }
    /// Return IDs within inclusive `rt +/- delta` when delta is positive.
    /// Otherwise return at most the earliest RT at or after `rt`, breaking ties
    /// by ID. Empty restrictions mean all spectra; finite query values are required.
    pub fn spectra_indices_by_rt(&self, rt: f64, delta: f64, indices: &[i64]) -> Result<Vec<i64>> {
        finite(rt)?;
        finite(delta)?;
        let mut budget = Budget::new(self.limits);
        let selected = selection(indices, true, &mut budget)?;
        let lower = if delta > 0.0 { rt - delta } else { rt };
        let upper = if delta > 0.0 { rt + delta } else { f64::MAX };
        finite(lower)?;
        finite(upper)?;
        let connector = self.open(false)?;
        let tx = connector
            .connection()
            .unchecked_transaction()
            .map_err(sql_error)?;
        self.schema(&tx)?;
        let result = {
            let mut stmt = tx.prepare("SELECT ID, RETENTION_TIME FROM SPECTRUM WHERE RETENTION_TIME BETWEEN ?1 AND ?2 ORDER BY RETENTION_TIME, ID").map_err(sql_error)?;
            let mut rows = stmt.query(params![lower, upper]).map_err(sql_error)?;
            let capacity = if delta > 0.0 {
                count(&tx, "SPECTRUM")?
            } else {
                1
            };
            let mut result = budget.work.vector(capacity)?;
            while let Some(row) = rows.next().map_err(sql_error)? {
                budget.work.spend(1)?;
                let id = nonnegative_id(row.get(0).map_err(sql_error)?)?;
                finite(row.get(1).map_err(sql_error)?)?;
                if selected.is_empty() || selected.binary_search(&id).is_ok() {
                    result.push(id);
                    if delta <= 0.0 {
                        break;
                    }
                }
            }
            if delta > 0.0 {
                result.sort_unstable();
            }
            result
        };
        tx.commit().map_err(sql_error)?;
        Ok(result)
    }
    /// Destructively replace the file with the seven-table schema and correct
    /// indexes. A complete sibling database is staged before replacement; writer
    /// IDs reset only after successful rename. Existing open connections are not supported.
    pub fn create_tables(&mut self) -> Result<()> {
        for suffix in ["-wal", "-shm", "-journal"] {
            let mut sidecar = self.filename.as_os_str().to_os_string();
            sidecar.push(suffix);
            if Path::new(&sidecar).try_exists()? {
                return Err(invalid(
                    "refusing database replacement while SQLite sidecars exist",
                ));
            }
        }
        let parent = self
            .filename
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let directory = crate::system::file::TempDir::new_in(parent, false)?;
        let stage = directory.path().join("sqmass.sqlite");
        let connector = SqliteConnector::new(&stage)?;
        connector.execute_statement(SCHEMA)?;
        self.schema(connector.connection())?;
        drop(connector);
        std::fs::rename(stage, &self.filename)?;
        self.spectrum_id = 0;
        self.chromatogram_id = 0;
        Ok(())
    }
    /// Atomically write RUN, optional RUN_EXTRA, chromatograms and spectra.
    /// Tables must already exist and hold no records; any failure rolls back all
    /// rows and leaves both writer counters unchanged.
    pub fn write_experiment(&mut self, experiment: &MSExperiment) -> Result<()> {
        let mut budget = Budget::new(self.limits);
        self.preflight(experiment, self.config.write_full_meta, &mut budget)?;
        let connector = self.open(true)?;
        let tx = connector
            .connection()
            .unchecked_transaction()
            .map_err(sql_error)?;
        self.schema(&tx)?;
        if count(&tx, "SPECTRUM")? != 0 || count(&tx, "CHROMATOGRAM")? != 0 {
            return Err(invalid("write_experiment requires empty record tables"));
        }
        self.write_run(&tx, experiment, self.config.write_full_meta, &mut budget)?;
        self.store_chromatograms(&tx, &experiment.chromatograms, 0, &mut budget)?;
        self.store_spectra(&tx, &experiment.spectra, 0, &mut budget)?;
        self.schema(&tx)?;
        tx.commit().map_err(sql_error)?;
        self.spectrum_id = experiment.spectra.len() as i64;
        self.chromatogram_id = experiment.chromatograms.len() as i64;
        Ok(())
    }
    /// Append spectra atomically using this handler's next IDs. SQL-only metadata
    /// loss requires Source policy. Low-level builder calls may leave an
    /// incomplete snapshot until all of its expected records have been written.
    pub fn write_spectra(&mut self, spectra: &[MSSpectrum]) -> Result<()> {
        if spectra.is_empty() {
            return Ok(());
        }
        let mut budget = Budget::new(self.limits);
        let connector = self.open(true)?;
        let tx = connector
            .connection()
            .unchecked_transaction()
            .map_err(sql_error)?;
        self.schema(&tx)?;
        self.append_guard(&tx, "SPECTRUM", self.spectrum_id, spectra.len())?;
        let snapshot = if self.config.loss_policy == LossPolicy::Reject {
            self.snapshot(&tx, self.run_id, &mut budget)?
        } else {
            None
        };
        self.preflight_spectra(spectra, snapshot.is_some(), &mut budget)?;
        if let Some(snapshot) = snapshot {
            let start =
                usize::try_from(self.spectrum_id).map_err(|_| invalid("invalid writer ID"))?;
            let expected = snapshot
                .spectra
                .get(start..start + spectra.len())
                .ok_or_else(|| invalid("append exceeds expected snapshot spectra"))?;
            if spectra
                .iter()
                .zip(expected)
                .any(|(value, expected)| spectrum_metadata(value) != *expected)
            {
                return Err(invalid(
                    "appended spectrum metadata does not match RUN_EXTRA",
                ));
            }
        }
        self.store_spectra(&tx, spectra, self.spectrum_id, &mut budget)?;
        self.schema(&tx)?;
        tx.commit().map_err(sql_error)?;
        self.spectrum_id += spectra.len() as i64;
        Ok(())
    }
    /// Chromatogram counterpart of [`Self::write_spectra`]. A failed write leaves
    /// rows and the next chromatogram ID unchanged.
    pub fn write_chromatograms(&mut self, chromatograms: &[MSChromatogram]) -> Result<()> {
        if chromatograms.is_empty() {
            return Ok(());
        }
        let mut budget = Budget::new(self.limits);
        let connector = self.open(true)?;
        let tx = connector
            .connection()
            .unchecked_transaction()
            .map_err(sql_error)?;
        self.schema(&tx)?;
        self.append_guard(
            &tx,
            "CHROMATOGRAM",
            self.chromatogram_id,
            chromatograms.len(),
        )?;
        let snapshot = if self.config.loss_policy == LossPolicy::Reject {
            self.snapshot(&tx, self.run_id, &mut budget)?
        } else {
            None
        };
        self.preflight_chromatograms(chromatograms, snapshot.is_some(), &mut budget)?;
        if let Some(snapshot) = snapshot {
            let start =
                usize::try_from(self.chromatogram_id).map_err(|_| invalid("invalid writer ID"))?;
            let expected = snapshot
                .chromatograms
                .get(start..start + chromatograms.len())
                .ok_or_else(|| invalid("append exceeds expected snapshot chromatograms"))?;
            if chromatograms
                .iter()
                .zip(expected)
                .any(|(value, expected)| chromatogram_metadata(value) != *expected)
            {
                return Err(invalid(
                    "appended chromatogram metadata does not match RUN_EXTRA",
                ));
            }
        }
        self.store_chromatograms(&tx, chromatograms, self.chromatogram_id, &mut budget)?;
        self.schema(&tx)?;
        tx.commit().map_err(sql_error)?;
        self.chromatogram_id += chromatograms.len() as i64;
        Ok(())
    }
    /// Low-level atomic RUN/optional RUN_EXTRA insertion. A snapshot may describe
    /// records that subsequent builder calls will add; full reads check it. This does not
    /// write peak arrays or update counters, and is not an upsert.
    pub fn write_run_level_information(
        &self,
        experiment: &MSExperiment,
        full_meta: bool,
    ) -> Result<()> {
        let mut budget = Budget::new(self.limits);
        self.preflight(experiment, full_meta, &mut budget)?;
        let connector = self.open(true)?;
        let tx = connector
            .connection()
            .unchecked_transaction()
            .map_err(sql_error)?;
        self.schema(&tx)?;
        self.write_run(&tx, experiment, full_meta, &mut budget)?;
        self.schema(&tx)?;
        tx.commit().map_err(sql_error)?;
        Ok(())
    }
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(format!("sqMass: {message}"))
}
fn sql_error(error: rusqlite::Error) -> Error {
    Error::Io(std::io::Error::other(error))
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("nonfinite numeric value"))
    }
}
fn check(n: usize, maximum: usize, label: &str) -> Result<()> {
    if n > maximum {
        Err(invalid(&format!("{label} limit exceeded")))
    } else {
        Ok(())
    }
}
fn nonnegative_id(id: i64) -> Result<i64> {
    if id >= 0 {
        Ok(id)
    } else {
        Err(invalid("negative SQL record/run ID"))
    }
}
fn count(conn: &Connection, table: &str) -> Result<usize> {
    let value: i64 = conn
        .query_row(&format!("SELECT count(*) FROM {table}"), [], |r| r.get(0))
        .map_err(sql_error)?;
    usize::try_from(value).map_err(|_| invalid("row count exceeds usize"))
}
fn ids(conn: &Connection, table: &str, budget: &mut Budget) -> Result<Vec<i64>> {
    let mut stmt = conn
        .prepare(&format!("SELECT ID FROM {table} ORDER BY ID"))
        .map_err(sql_error)?;
    let mut rows = stmt.query([]).map_err(sql_error)?;
    let mut result = budget.work.vector(count(conn, table)?)?;
    while let Some(row) = rows.next().map_err(sql_error)? {
        budget.work.spend(1)?;
        result.push(nonnegative_id(row.get(0).map_err(sql_error)?)?);
    }
    Ok(result)
}
fn selection(ids: &[i64], empty: bool, budget: &mut Budget) -> Result<Vec<i64>> {
    check(ids.len(), budget.limits.max_records, "selection")?;
    if ids.is_empty() && !empty {
        return Err(invalid("selection must not be empty"));
    }
    let mut values = budget.work.vector(ids.len())?;
    budget
        .work
        .spend(ids.len().saturating_mul(usize::BITS as usize))?;
    values.extend_from_slice(ids);
    values.sort_unstable();
    if values.first().is_some_and(|id| *id < 0) || values.windows(2).any(|p| p[0] == p[1]) {
        return Err(invalid("negative or duplicate selected IDs"));
    }
    Ok(values)
}
fn read_run(conn: &Connection, budget: &mut Budget) -> Result<(i64, String)> {
    let mut stmt = conn
        .prepare("SELECT ID,FILENAME,NATIVE_ID FROM RUN")
        .map_err(sql_error)?;
    let mut rows = stmt.query([]).map_err(sql_error)?;
    let row = rows
        .next()
        .map_err(sql_error)?
        .ok_or_else(|| invalid("exactly one RUN required"))?;
    let id = nonnegative_id(row.get(0).map_err(sql_error)?)?;
    let filename = bounded_text(row, 1, budget)?;
    let _ = bounded_text(row, 2, budget)?;
    if rows.next().map_err(sql_error)?.is_some() {
        return Err(invalid("exactly one RUN required"));
    }
    Ok((id, filename))
}

struct Budget {
    work: Work,
    values: usize,
    records: usize,
    limits: HandlerLimits,
}
impl Budget {
    fn new(limits: HandlerLimits) -> Self {
        Self {
            work: Work::new(NumpressCoderLimits {
                raw: NumpressLimits {
                    max_values: limits.max_values_per_array,
                    max_encoded_bytes: limits.max_blob_bytes,
                    max_work: limits.max_work,
                },
                max_text_bytes: limits.max_metadata_bytes,
                max_total_bytes: limits.max_total_bytes,
            }),
            values: 0,
            records: 0,
            limits,
        }
    }
    fn records<T>(&mut self, n: usize) -> Result<()> {
        self.records = self
            .records
            .checked_add(n)
            .ok_or_else(|| invalid("record overflow"))?;
        check(self.records, self.limits.max_records, "record")?;
        self.work.allocate(
            n.checked_mul(std::mem::size_of::<T>())
                .ok_or_else(|| invalid("record byte overflow"))?,
        )
    }
    fn values(&mut self, n: usize) -> Result<()> {
        check(n, self.limits.max_values_per_array, "array values")?;
        self.values = self
            .values
            .checked_add(n)
            .ok_or_else(|| invalid("value overflow"))?;
        check(self.values, self.limits.max_total_values, "total values")
    }
}

impl MzMLSqliteHandler {
    fn snapshot(
        &self,
        conn: &Connection,
        run_id: i64,
        budget: &mut Budget,
    ) -> Result<Option<MSExperiment>> {
        let mut result = None;
        let mut statement = conn
            .prepare("SELECT RUN_ID, DATA FROM RUN_EXTRA")
            .map_err(sql_error)?;
        let mut rows = statement.query([]).map_err(sql_error)?;
        if let Some(row) = rows.next().map_err(sql_error)? {
            if row.get::<_, i64>(0).map_err(sql_error)? != run_id {
                return Err(invalid("RUN_EXTRA run ID mismatch"));
            }
            let blob = row
                .get_ref(1)
                .map_err(sql_error)?
                .as_blob()
                .map_err(|_| invalid("RUN_EXTRA is not a BLOB"))?;
            check(blob.len(), self.limits.max_metadata_bytes, "metadata BLOB")?;
            if !blob.is_empty() {
                let old = budget.work.limits.raw.max_encoded_bytes;
                budget.work.limits.raw.max_encoded_bytes = self.limits.max_metadata_bytes;
                let xml = codec::zlib_decode(blob, &mut budget.work)?;
                budget.work.limits.raw.max_encoded_bytes = old;
                // Count borrowed XML start events before the full parser can
                // allocate record structures. Declared mzML list counts are
                // advisory and cannot safely size this allowance.
                let snapshot_records = snapshot_record_slots(&xml, budget)?;
                // Reserve the embedded reader's parameter and array allowances
                // separately from the remaining operation budget before parsing.
                let parser_bytes = self
                    .limits
                    .max_metadata_bytes
                    .min(budget.work.remaining_bytes() / 4);
                budget.work.allocate(parser_bytes * 2)?;
                let parser_work = xml
                    .len()
                    .checked_mul(64)
                    .and_then(|n| n.checked_add(1024))
                    .ok_or_else(|| invalid("metadata parse work overflow"))?;
                budget.work.spend(parser_work)?;
                let options = super::mzml::ReadOptions {
                    max_xml_bytes: self.limits.max_metadata_bytes as u64,
                    max_records: snapshot_records,
                    max_array_bytes: self.limits.max_blob_bytes,
                    max_total_peaks: 0,
                    max_total_array_bytes: parser_bytes,
                    max_total_array_elements: self.limits.max_total_values,
                    max_total_arrays: self.limits.max_records.saturating_mul(8),
                    max_param_bytes: parser_bytes,
                    max_total_params: parser_work.min(budget.work.remaining_work()),
                    ..Default::default()
                };
                result = Some(super::mzml::read_with_options(Cursor::new(xml), &options)?);
            }
            if rows.next().map_err(sql_error)?.is_some() {
                return Err(invalid("multiple RUN_EXTRA snapshots"));
            }
        }
        Ok(result)
    }
    fn open(&self, write: bool) -> Result<SqliteConnector> {
        let metadata = std::fs::metadata(&self.filename)?;
        if !metadata.is_file() || metadata.len() > self.limits.max_database_bytes {
            return Err(invalid("database file limit or non-file input"));
        }
        let conn = SqliteConnector::open(
            &self.filename,
            if write {
                SqlOpenMode::ReadWrite
            } else {
                SqlOpenMode::ReadOnly
            },
        )?;
        let length = self
            .limits
            .max_blob_bytes
            .max(self.limits.max_metadata_bytes)
            .checked_add(4096)
            .and_then(|v| i32::try_from(v).ok())
            .ok_or_else(|| invalid("SQLite length limit exceeds i32"))?;
        conn.connection()
            .set_limit(rusqlite::limits::Limit::SQLITE_LIMIT_LENGTH, length)
            .map_err(sql_error)?;
        Ok(conn)
    }
    // ponytail: global validation keeps malformed references deterministic;
    // every operation scans the bounded database, including single-ID reads.
    // Reuse a validated snapshot or scope checks only if profiling requires it.
    fn schema(&self, conn: &Connection) -> Result<()> {
        let pages: i64 = conn
            .query_row("PRAGMA page_count", [], |r| r.get(0))
            .map_err(sql_error)?;
        let page_size: i64 = conn
            .query_row("PRAGMA page_size", [], |r| r.get(0))
            .map_err(sql_error)?;
        if pages
            .checked_mul(page_size)
            .is_none_or(|n| n < 0 || n as u64 > self.limits.max_database_bytes)
        {
            return Err(invalid("logical database size limit"));
        }
        let mut records = 0usize;
        for table in [
            "RUN",
            "RUN_EXTRA",
            "SPECTRUM",
            "CHROMATOGRAM",
            "DATA",
            "PRECURSOR",
            "PRODUCT",
        ] {
            let kind: Option<String> = conn.query_row("SELECT type FROM pragma_table_list() WHERE schema='main' AND name=?1 COLLATE NOCASE", [table], |r| r.get(0)).optional().map_err(sql_error)?;
            if kind.as_deref() != Some("table") {
                return Err(invalid(
                    "required ordinary main table missing (views/virtual tables are unsupported)",
                ));
            }
            let n = count(conn, table)?;
            check(
                n,
                self.limits.max_records.saturating_mul(2).max(1),
                "SQL rows",
            )?;
            if table == "SPECTRUM" || table == "CHROMATOGRAM" {
                records = records
                    .checked_add(n)
                    .ok_or_else(|| invalid("SQL record count overflow"))?;
            }
        }
        // Stored payload is bounded by the page-count check against
        // `max_database_bytes` above. `max_total_bytes` is an allocation budget,
        // charged where an operation materializes values, so it must not reject a
        // file merely for storing more than one operation may load.
        check(records, self.limits.max_records, "SQL records")?;
        for table in ["SPECTRUM", "CHROMATOGRAM"] {
            let bad: bool = conn.query_row(&format!("SELECT EXISTS(SELECT 1 FROM {table} WHERE ID IS NULL OR typeof(ID)!='integer' OR ID<0 OR NATIVE_ID IS NULL)"), [], |r| r.get(0)).map_err(sql_error)?;
            if bad {
                return Err(invalid("invalid SQL record identity"));
            }
        }
        for table in ["DATA", "PRECURSOR", "PRODUCT"] {
            let bad: bool = conn.query_row(&format!("SELECT EXISTS(SELECT 1 FROM {table} AS X WHERE (SPECTRUM_ID IS NULL)=(CHROMATOGRAM_ID IS NULL) OR (SPECTRUM_ID IS NOT NULL AND NOT EXISTS(SELECT 1 FROM SPECTRUM WHERE ID=X.SPECTRUM_ID)) OR (CHROMATOGRAM_ID IS NOT NULL AND NOT EXISTS(SELECT 1 FROM CHROMATOGRAM WHERE ID=X.CHROMATOGRAM_ID)))"), [], |r| r.get(0)).map_err(sql_error)?;
            if bad {
                return Err(invalid("ambiguous or nonexistent SQL object reference"));
            }
        }
        Ok(())
    }
    fn spectra(
        &self,
        conn: &Connection,
        ids: &[i64],
        meta_only: bool,
        budget: &mut Budget,
    ) -> Result<Vec<MSSpectrum>> {
        budget.records::<MSSpectrum>(ids.len())?;
        let mut result = Vec::with_capacity(ids.len());
        let mut stmt = conn
            .prepare(
                "SELECT NATIVE_ID,MSLEVEL,RETENTION_TIME,SCAN_POLARITY FROM SPECTRUM WHERE ID=?1",
            )
            .map_err(sql_error)?;
        for id in ids {
            let mut rows = stmt.query([id]).map_err(sql_error)?;
            let row = rows
                .next()
                .map_err(sql_error)?
                .ok_or_else(|| invalid("missing selected spectrum ID"))?;
            let mut spectrum = MSSpectrum {
                native_id: bounded_text(row, 0, budget)?,
                ..Default::default()
            };
            if let Some(value) = row.get::<_, Option<i64>>(1).map_err(sql_error)? {
                spectrum.ms_level =
                    u32::try_from(value).map_err(|_| invalid("invalid MS level"))?;
            }
            if let Some(value) = row.get::<_, Option<f64>>(2).map_err(sql_error)? {
                spectrum.rt = finite(value)?;
            }
            spectrum.instrument_settings.polarity =
                match row.get::<_, Option<i64>>(3).map_err(sql_error)? {
                    None => Polarity::Unknown,
                    Some(0) => Polarity::Negative,
                    Some(1) => Polarity::Positive,
                    _ => return Err(invalid("invalid scan polarity")),
                };
            if rows.next().map_err(sql_error)?.is_some() {
                return Err(invalid("duplicate spectrum ID"));
            }
            drop(rows);
            if let Some(precursor) = read_precursor(conn, *id, false, budget)? {
                spectrum.precursors.push(precursor);
            }
            if let Some(product) = read_product(conn, *id, false)? {
                spectrum.products.push(product);
            }
            if !meta_only {
                let (coordinate, intensity) = self.arrays(conn, *id, false, budget)?;
                spectrum.peaks = coordinate
                    .into_iter()
                    .zip(intensity)
                    .map(|(mz, intensity)| Peak1D { mz, intensity })
                    .collect();
            }
            result.push(spectrum);
        }
        Ok(result)
    }
    fn chromatograms(
        &self,
        conn: &Connection,
        ids: &[i64],
        meta_only: bool,
        budget: &mut Budget,
    ) -> Result<Vec<MSChromatogram>> {
        budget.records::<MSChromatogram>(ids.len())?;
        let mut result = Vec::with_capacity(ids.len());
        let mut stmt = conn
            .prepare("SELECT NATIVE_ID FROM CHROMATOGRAM WHERE ID=?1")
            .map_err(sql_error)?;
        for id in ids {
            let mut rows = stmt.query([id]).map_err(sql_error)?;
            let row = rows
                .next()
                .map_err(sql_error)?
                .ok_or_else(|| invalid("missing selected chromatogram ID"))?;
            let mut chrom = MSChromatogram {
                native_id: bounded_text(row, 0, budget)?,
                precursor: read_precursor(conn, *id, true, budget)?
                    .ok_or_else(|| invalid("chromatogram precursor missing"))?,
                product: read_product(conn, *id, true)?
                    .ok_or_else(|| invalid("chromatogram product missing"))?,
                ..Default::default()
            };
            if rows.next().map_err(sql_error)?.is_some() {
                return Err(invalid("duplicate chromatogram ID"));
            }
            if !meta_only {
                let (coordinate, intensity) = self.arrays(conn, *id, true, budget)?;
                chrom.peaks = coordinate
                    .into_iter()
                    .zip(intensity)
                    .map(|(rt, intensity)| ChromatogramPeak { rt, intensity })
                    .collect();
            }
            result.push(chrom);
        }
        Ok(result)
    }
    fn arrays(
        &self,
        conn: &Connection,
        id: i64,
        chrom: bool,
        budget: &mut Budget,
    ) -> Result<(Vec<f64>, Vec<f32>)> {
        let mut stmt = conn
            .prepare(if chrom {
                "SELECT COMPRESSION,DATA_TYPE,DATA FROM DATA WHERE CHROMATOGRAM_ID=?1"
            } else {
                "SELECT COMPRESSION,DATA_TYPE,DATA FROM DATA WHERE SPECTRUM_ID=?1"
            })
            .map_err(sql_error)?;
        let mut rows = stmt.query([id]).map_err(sql_error)?;
        let mut coordinate = None;
        let mut intensity = None;
        while let Some(row) = rows.next().map_err(sql_error)? {
            let compression: i64 = row.get(0).map_err(sql_error)?;
            let role: i64 = row.get(1).map_err(sql_error)?;
            let slot = if role == 1 {
                &mut intensity
            } else if role == if chrom { 2 } else { 0 } {
                &mut coordinate
            } else {
                return Err(invalid("unexpected DATA role"));
            };
            if slot.is_some() {
                return Err(invalid("duplicate DATA role"));
            }
            if ![1, 5, 6].contains(&compression) {
                return Err(invalid("unsupported compression code"));
            }
            let blob = row
                .get_ref(2)
                .map_err(sql_error)?
                .as_blob()
                .map_err(|_| invalid("DATA is not a BLOB"))?;
            budget.work.binary(blob.len())?;
            let decoded = codec::zlib_decode(blob, &mut budget.work)?;
            let remaining_values = self.limits.max_total_values.saturating_sub(budget.values);
            budget.work.limits.raw.max_values =
                self.limits.max_values_per_array.min(remaining_values);
            let values = if compression == 1 {
                if decoded.len() % 8 != 0 {
                    return Err(invalid(
                        "lossless array byte length is not a multiple of eight",
                    ));
                }
                let n = decoded.len() / 8;
                check(n, budget.work.limits.raw.max_values, "decoded values")?;
                let mut values = budget.work.vector(n)?;
                values.extend(
                    decoded
                        .chunks_exact(8)
                        .map(|v| f64::from_le_bytes(v.try_into().expect("eight-byte chunk"))),
                );
                values
            } else {
                codec::decode_raw(
                    &decoded,
                    if compression == 5 {
                        NumpressCompression::Linear
                    } else {
                        NumpressCompression::Slof
                    },
                    &mut budget.work,
                )?
            };
            budget.values(values.len())?;
            for value in &values {
                finite(*value)?;
            }
            *slot = Some(values);
        }
        let coordinate = coordinate.ok_or_else(|| invalid("missing coordinate DATA role"))?;
        let intensity = intensity.ok_or_else(|| invalid("missing intensity DATA role"))?;
        if coordinate.len() != intensity.len() {
            return Err(invalid("coordinate/intensity lengths differ"));
        }
        // Account both the returned f32 vector and the final native peak vector.
        budget.work.allocate(
            coordinate
                .len()
                .checked_mul(20)
                .ok_or_else(|| invalid("peak allocation overflow"))?,
        )?;
        let intensity = intensity
            .into_iter()
            .map(|v| {
                let value = v as f32;
                if value.is_finite() {
                    Ok(value)
                } else {
                    Err(invalid("intensity exceeds f32"))
                }
            })
            .collect::<Result<Vec<_>>>()?;
        Ok((coordinate, intensity))
    }
    fn append_guard(&self, conn: &Connection, table: &str, next: i64, added: usize) -> Result<()> {
        let n = count(conn, table)?;
        if n != usize::try_from(next).map_err(|_| invalid("negative writer counter"))? {
            return Err(invalid(
                "existing database does not match this handler's writer counter",
            ));
        }
        if n != 0 {
            let (min, max): (i64, i64) = conn
                .query_row(&format!("SELECT min(ID),max(ID) FROM {table}"), [], |r| {
                    Ok((r.get(0)?, r.get(1)?))
                })
                .map_err(sql_error)?;
            if min != 0 || max != next - 1 {
                return Err(invalid("writer requires its contiguous assigned IDs"));
            }
        }
        let total = count(conn, "SPECTRUM")?
            .checked_add(count(conn, "CHROMATOGRAM")?)
            .and_then(|n| n.checked_add(added))
            .ok_or_else(|| invalid("record count overflow"))?;
        check(total, self.limits.max_records, "written records")?;
        next.checked_add(i64::try_from(added).map_err(|_| invalid("writer ID overflow"))?)
            .ok_or_else(|| invalid("writer ID overflow"))?;
        Ok(())
    }
    fn store_spectra(
        &self,
        conn: &Connection,
        spectra: &[MSSpectrum],
        first: i64,
        budget: &mut Budget,
    ) -> Result<()> {
        let mut statement = conn.prepare("INSERT INTO SPECTRUM(ID,RUN_ID,MSLEVEL,RETENTION_TIME,SCAN_POLARITY,NATIVE_ID) VALUES(?1,?2,?3,?4,?5,?6)").map_err(sql_error)?;
        let mut data_statement=conn.prepare("INSERT INTO DATA(SPECTRUM_ID,CHROMATOGRAM_ID,COMPRESSION,DATA_TYPE,DATA) VALUES(?1,?2,?3,?4,?5)").map_err(sql_error)?;
        for (chunk, records) in spectra.chunks(self.config.sql_batch_size).enumerate() {
            for (offset, spectrum) in records.iter().enumerate() {
                let id = first + (chunk * self.config.sql_batch_size + offset) as i64;
                statement
                    .execute(params![
                        id,
                        self.run_id,
                        spectrum.ms_level,
                        spectrum.rt,
                        i64::from(spectrum.instrument_settings.polarity == Polarity::Positive),
                        spectrum.native_id
                    ])
                    .map_err(sql_error)?;
                if let Some(p) = spectrum.precursors.first() {
                    write_precursor(conn, id, false, p)?;
                }
                if let Some(p) = spectrum.products.first() {
                    write_product(conn, id, false, p)?;
                }
                let mut coordinates = budget.work.vector(spectrum.peaks.len())?;
                coordinates.extend(spectrum.peaks.iter().map(|p| p.mz));
                self.store_array(
                    &mut data_statement,
                    (id, false),
                    0,
                    &coordinates,
                    self.config.linear_abs_mass_acc,
                    budget,
                )?;
                coordinates.clear();
                coordinates.extend(spectrum.peaks.iter().map(|p| f64::from(p.intensity)));
                self.store_array(
                    &mut data_statement,
                    (id, false),
                    1,
                    &coordinates,
                    0.0,
                    budget,
                )?;
            }
        }
        Ok(())
    }
    fn store_chromatograms(
        &self,
        conn: &Connection,
        chromatograms: &[MSChromatogram],
        first: i64,
        budget: &mut Budget,
    ) -> Result<()> {
        let mut statement = conn
            .prepare("INSERT INTO CHROMATOGRAM(ID,RUN_ID,NATIVE_ID) VALUES(?1,?2,?3)")
            .map_err(sql_error)?;
        let mut data_statement=conn.prepare("INSERT INTO DATA(SPECTRUM_ID,CHROMATOGRAM_ID,COMPRESSION,DATA_TYPE,DATA) VALUES(?1,?2,?3,?4,?5)").map_err(sql_error)?;
        for (chunk, records) in chromatograms.chunks(self.config.sql_batch_size).enumerate() {
            for (offset, chrom) in records.iter().enumerate() {
                let id = first + (chunk * self.config.sql_batch_size + offset) as i64;
                statement
                    .execute(params![id, self.run_id, chrom.native_id])
                    .map_err(sql_error)?;
                write_precursor(conn, id, true, &chrom.precursor)?;
                write_product(conn, id, true, &chrom.product)?;
                let mut coordinates = budget.work.vector(chrom.peaks.len())?;
                coordinates.extend(chrom.peaks.iter().map(|p| p.rt));
                self.store_array(
                    &mut data_statement,
                    (id, true),
                    2,
                    &coordinates,
                    0.05,
                    budget,
                )?;
                coordinates.clear();
                coordinates.extend(chrom.peaks.iter().map(|p| f64::from(p.intensity)));
                self.store_array(
                    &mut data_statement,
                    (id, true),
                    1,
                    &coordinates,
                    0.0,
                    budget,
                )?;
            }
        }
        Ok(())
    }
    fn store_array(
        &self,
        statement: &mut rusqlite::Statement<'_>,
        owner: (i64, bool),
        role: i64,
        values: &[f64],
        accuracy: f64,
        budget: &mut Budget,
    ) -> Result<()> {
        let (id, chrom) = owner;
        // The source's positive-accuracy linear estimator returns zero for
        // one or two points. Lossless short coordinates avoid a zero-factor
        // payload (including all-zero coordinates) using supported code 1.
        let (compression, raw) =
            if self.config.use_lossy_compression && !(role != 1 && values.len() < 3) {
                let config = NumpressConfig {
                    compression: if role == 1 {
                        NumpressCompression::Slof
                    } else {
                        NumpressCompression::Linear
                    },
                    error_tolerance: -1.0,
                    linear_fp_mass_acc: accuracy,
                    estimate_fixed_point: true,
                    ..Default::default()
                };
                let report = codec::encode_raw(values, &config, &mut budget.work)?;
                match report.status {
                    NumpressEncodeStatus::Encoded | NumpressEncodeStatus::EmptyInput => {}
                    _ => return Err(invalid("Numpress encoding rejected; no DATA row written")),
                }
                if !values.is_empty() {
                    let factor = report
                        .fixed_point
                        .filter(|v| v.is_finite() && *v > 0.0)
                        .ok_or_else(|| invalid("Numpress produced an unusable fixed point"))?;
                    // Linear stores its first two values as unsigned 32-bit
                    // integers. The raw source-compatible codec retains their
                    // low bits, so this storage adapter must reject wrapping.
                    if role != 1
                        && values.iter().take(2).any(|value| {
                            let quantized = (value * factor + 0.5).trunc();
                            !(0.0..4294967296.0).contains(&quantized)
                        })
                    {
                        return Err(invalid(
                            "Numpress initial coordinates exceed unsigned 32-bit quantization",
                        ));
                    }
                }
                (if role == 1 { 6 } else { 5 }, report.output)
            } else {
                let raw_bytes = values
                    .len()
                    .checked_mul(8)
                    .ok_or_else(|| invalid("array size overflow"))?;
                budget.work.binary(raw_bytes)?;
                let mut raw = budget.work.vector(raw_bytes)?;
                for value in values {
                    raw.extend_from_slice(&value.to_le_bytes());
                }
                (1, raw)
            };
        budget.work.binary(raw.len())?;
        let blob = codec::zlib_encode(&raw, &mut budget.work)?;
        statement
            .execute(params![
                if chrom { None } else { Some(id) },
                if chrom { Some(id) } else { None },
                compression,
                role,
                blob
            ])
            .map_err(sql_error)?;
        Ok(())
    }
    fn write_run(
        &self,
        conn: &Connection,
        experiment: &MSExperiment,
        full: bool,
        budget: &mut Budget,
    ) -> Result<()> {
        let filename = &experiment.settings.document.loaded_file_path;
        conn.execute(
            "INSERT INTO RUN(ID,FILENAME,NATIVE_ID) VALUES(?1,?2,?2)",
            params![self.run_id, filename],
        )
        .map_err(sql_error)?;
        if full {
            let snapshot = MSExperiment {
                spectra: experiment.spectra.iter().map(spectrum_metadata).collect(),
                chromatograms: experiment
                    .chromatograms
                    .iter()
                    .map(chromatogram_metadata)
                    .collect(),
                settings: experiment.settings.clone(),
                sql_run_id: 0,
            };
            let mut writer = BoundedWriter {
                bytes: Vec::new(),
                maximum: self
                    .limits
                    .max_metadata_bytes
                    .min(budget.work.remaining_bytes()),
            };
            super::mzml::write(&mut writer, &snapshot)?;
            budget.work.allocate(writer.bytes.len())?;
            let old = budget.work.limits.raw.max_encoded_bytes;
            budget.work.limits.raw.max_encoded_bytes = self.limits.max_metadata_bytes;
            let blob = codec::zlib_encode(&writer.bytes, &mut budget.work)?;
            budget.work.limits.raw.max_encoded_bytes = old;
            conn.execute(
                "INSERT INTO RUN_EXTRA(RUN_ID,DATA) VALUES(?1,?2)",
                params![self.run_id, blob],
            )
            .map_err(sql_error)?;
        }
        Ok(())
    }
}

// The nested mzML parser has its own fixed header registry limits. This
// preflight meters record slots against the enclosing operation, in addition
// to the parameter/array allowances passed through its public ReadOptions.
fn snapshot_record_slots(xml: &[u8], budget: &mut Budget) -> Result<usize> {
    budget.work.spend(xml.len())?;
    let mut reader = quick_xml::Reader::from_reader(xml);
    let mut records = 0usize;
    loop {
        match reader.read_event().map_err(|e| invalid(&e.to_string()))? {
            quick_xml::events::Event::Start(tag) | quick_xml::events::Event::Empty(tag) => {
                let bytes = match tag.local_name().as_ref() {
                    b"spectrum" => std::mem::size_of::<MSSpectrum>(),
                    b"chromatogram" => std::mem::size_of::<MSChromatogram>(),
                    _ => continue,
                };
                records = records
                    .checked_add(1)
                    .ok_or_else(|| invalid("snapshot record overflow"))?;
                check(records, budget.limits.max_records, "snapshot records")?;
                budget.work.allocate(bytes)?;
            }
            quick_xml::events::Event::Eof => return Ok(records),
            _ => {}
        }
    }
}

fn bounded_text(row: &Row<'_>, column: usize, budget: &mut Budget) -> Result<String> {
    let value = row.get_ref(column).map_err(sql_error)?;
    let bytes = match value {
        rusqlite::types::ValueRef::Text(v) | rusqlite::types::ValueRef::Blob(v) => v,
        _ => return Err(invalid("expected SQL text")),
    };
    budget.work.binary(bytes.len())?;
    budget.work.allocate(bytes.len())?;
    budget.work.spend(bytes.len())?;
    std::str::from_utf8(bytes)
        .map(str::to_owned)
        .map_err(|_| invalid("invalid UTF-8 SQL text"))
}
fn read_precursor(
    conn: &Connection,
    id: i64,
    chrom: bool,
    budget: &mut Budget,
) -> Result<Option<Precursor>> {
    let owner = if chrom {
        "CHROMATOGRAM_ID"
    } else {
        "SPECTRUM_ID"
    };
    let mut stmt=conn.prepare(&format!("SELECT CHARGE,PEPTIDE_SEQUENCE,DRIFT_TIME,ACTIVATION_METHOD,ACTIVATION_ENERGY,ISOLATION_TARGET,ISOLATION_LOWER,ISOLATION_UPPER FROM PRECURSOR WHERE {owner}=?1")).map_err(sql_error)?;
    let mut rows = stmt.query([id]).map_err(sql_error)?;
    let Some(row) = rows.next().map_err(sql_error)? else {
        return Ok(None);
    };
    let target: Option<f64> = row.get(5).map_err(sql_error)?;
    let charge = row
        .get::<_, Option<i32>>(0)
        .map_err(sql_error)?
        .unwrap_or(0);
    let mut p = Precursor {
        charge,
        ..Precursor::default()
    };
    if row.get_ref(1).map_err(sql_error)? != rusqlite::types::ValueRef::Null {
        p.cv_terms.metadata.insert(
            "peptide_sequence".into(),
            MetaValue::from(bounded_text(row, 1, budget)?),
        );
    }
    p.drift_time = row
        .get::<_, Option<f64>>(2)
        .map_err(sql_error)?
        .filter(|v| *v != -1.0);
    if let Some(d) = p.drift_time {
        finite(d)?;
    }
    if let Some(method) = row.get::<_, Option<i64>>(3).map_err(sql_error)? {
        if method != -1 {
            let index =
                usize::try_from(method).map_err(|_| invalid("invalid activation method"))?;
            p.activation_methods.insert(
                *ActivationMethod::ALL
                    .get(index)
                    .ok_or_else(|| invalid("invalid activation method"))?,
            );
        }
    }
    p.activation_energy = finite(
        row.get::<_, Option<f64>>(4)
            .map_err(sql_error)?
            .unwrap_or(0.0),
    )?;
    p.mz = finite(target.unwrap_or(0.0))?;
    p.isolation_window_lower_offset = offset(row, 6)?;
    p.isolation_window_upper_offset = offset(row, 7)?;
    p.validate()?;
    if rows.next().map_err(sql_error)?.is_some() {
        return Err(invalid("multiple precursor rows for one object"));
    }
    Ok(if target.is_some() || chrom {
        Some(p)
    } else {
        None
    })
}
fn read_product(conn: &Connection, id: i64, chrom: bool) -> Result<Option<Product>> {
    let owner = if chrom {
        "CHROMATOGRAM_ID"
    } else {
        "SPECTRUM_ID"
    };
    let mut stmt = conn
        .prepare(&format!(
            "SELECT ISOLATION_TARGET,ISOLATION_LOWER,ISOLATION_UPPER FROM PRODUCT WHERE {owner}=?1"
        ))
        .map_err(sql_error)?;
    let mut rows = stmt.query([id]).map_err(sql_error)?;
    let Some(row) = rows.next().map_err(sql_error)? else {
        return Ok(None);
    };
    let target: Option<f64> = row.get(0).map_err(sql_error)?;
    let p = Product {
        mz: finite(target.unwrap_or(0.0))?,
        isolation_window_lower_offset: offset(row, 1)?,
        isolation_window_upper_offset: offset(row, 2)?,
        ..Default::default()
    };
    if rows.next().map_err(sql_error)?.is_some() {
        return Err(invalid("multiple product rows for one object"));
    }
    Ok(if target.is_some() || chrom {
        Some(p)
    } else {
        None
    })
}
fn offset(row: &Row<'_>, i: usize) -> Result<f64> {
    Ok(finite(
        row.get::<_, Option<f64>>(i)
            .map_err(sql_error)?
            .unwrap_or(0.0),
    )?
    .max(0.0))
}
fn write_precursor(conn: &Connection, id: i64, chrom: bool, p: &Precursor) -> Result<()> {
    let sequence = p
        .cv_terms
        .metadata
        .get("peptide_sequence")
        .map(ToString::to_string);
    let method = p
        .activation_methods
        .first()
        .map(|v| *v as i64)
        .unwrap_or(-1);
    conn.prepare("INSERT INTO PRECURSOR(SPECTRUM_ID,CHROMATOGRAM_ID,CHARGE,PEPTIDE_SEQUENCE,DRIFT_TIME,ACTIVATION_METHOD,ACTIVATION_ENERGY,ISOLATION_TARGET,ISOLATION_LOWER,ISOLATION_UPPER) VALUES(?1,?2,?3,?4,?5,?6,?7,?8,?9,?10)").map_err(sql_error)?
        .execute(params![if chrom{None}else{Some(id)},if chrom{Some(id)}else{None},p.charge,sequence,p.drift_time.unwrap_or(-1.0),method,p.activation_energy,p.mz,p.isolation_window_lower_offset,p.isolation_window_upper_offset]).map_err(sql_error)?;
    Ok(())
}
fn write_product(conn: &Connection, id: i64, chrom: bool, p: &Product) -> Result<()> {
    conn.prepare("INSERT INTO PRODUCT(SPECTRUM_ID,CHROMATOGRAM_ID,CHARGE,ISOLATION_TARGET,ISOLATION_LOWER,ISOLATION_UPPER) VALUES(?1,?2,0,?3,?4,?5)").map_err(sql_error)?
        .execute(params![if chrom{None}else{Some(id)},if chrom{Some(id)}else{None},p.mz,p.isolation_window_lower_offset,p.isolation_window_upper_offset]).map_err(sql_error)?;
    Ok(())
}

struct BoundedWriter {
    bytes: Vec<u8>,
    maximum: usize,
}
impl Write for BoundedWriter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if self
            .bytes
            .len()
            .checked_add(bytes.len())
            .is_none_or(|n| n > self.maximum)
        {
            return Err(std::io::Error::other("sqMass metadata byte limit"));
        }
        self.bytes
            .try_reserve_exact(bytes.len())
            .map_err(|_| std::io::Error::other("sqMass metadata allocation failed"))?;
        self.bytes.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl MzMLSqliteHandler {
    fn preflight(&self, experiment: &MSExperiment, full: bool, budget: &mut Budget) -> Result<()> {
        let mut work = self.limits.max_work;
        let mut bytes = self.limits.max_total_bytes;
        experiment.settings.with_budget(&mut work, &mut bytes)?;
        budget.work.spend(self.limits.max_work - work)?;
        budget.work.allocate(self.limits.max_total_bytes - bytes)?;
        if !full && self.config.loss_policy == LossPolicy::Reject {
            let mut basic = crate::metadata::ExperimentalSettings::default();
            basic.document.loaded_file_path = experiment.settings.document.loaded_file_path.clone();
            if basic != experiment.settings {
                return Err(Error::Unsupported("sqMass SQL-only run storage would discard experimental settings; select Source loss policy".into()));
            }
        }
        self.preflight_spectra(&experiment.spectra, full, budget)?;
        self.preflight_chromatograms(&experiment.chromatograms, full, budget)
    }
    fn preflight_spectra(
        &self,
        spectra: &[MSSpectrum],
        full: bool,
        budget: &mut Budget,
    ) -> Result<()> {
        self.configuration()?;
        budget.records::<MSSpectrum>(spectra.len())?;
        for s in spectra {
            let (initial_work, initial_bytes) =
                (budget.work.remaining_work(), budget.work.remaining_bytes());
            let (mut work, mut bytes) = (initial_work, initial_bytes);
            self.auxiliary(
                !s.float_data_arrays.is_empty()
                    || !s.integer_data_arrays.is_empty()
                    || !s.string_data_arrays.is_empty(),
            )?;
            if full && !s.peptide_identifications.is_empty() {
                return Err(Error::Unsupported(
                    "mzML snapshot cannot preserve peptide identifications".into(),
                ));
            }
            s.acquisition_with_budget(&mut work, &mut bytes)?;
            s.record_metadata_with_budget(&mut work, &mut bytes)?;
            let mut meter = crate::kernel::data_array::Meter {
                work: &mut work,
                bytes: &mut bytes,
            };
            meter.text(&s.native_id)?;
            meter.text(&s.name)?;
            meter.slots::<Precursor>(s.precursors.len())?;
            for p in &s.precursors {
                precursor_budget(p, &mut meter)?;
                p.validate()?;
            }
            for p in &s.products {
                p.validate()?;
            }
            if full && s.native_id.is_empty() {
                return Err(invalid("full metadata requires nonempty native IDs"));
            }
            finite(s.rt)?;
            finite(s.drift_time)?;
            budget.values(s.peaks.len())?;
            budget.values(s.peaks.len())?;
            budget.work.spend(s.peaks.len())?;
            for peak in &s.peaks {
                finite(peak.mz)?;
                finite(f64::from(peak.intensity))?;
            }
            budget.work.spend(initial_work - work)?;
            budget.work.allocate(
                (initial_bytes - bytes)
                    .checked_mul(2)
                    .ok_or_else(|| invalid("metadata clone budget overflow"))?,
            )?;
            if !full
                && self.config.loss_policy == LossPolicy::Reject
                && spectrum_metadata(s) != spectrum_projection(s)
            {
                return Err(Error::Unsupported("sqMass SQL-only spectrum storage would discard metadata; select Source loss policy".into()));
            }
            if !full
                && self.config.loss_policy == LossPolicy::Reject
                && !s.peptide_identifications.is_empty()
            {
                return Err(Error::Unsupported(
                    "sqMass discards peptide identifications".into(),
                ));
            }
        }
        Ok(())
    }
    fn preflight_chromatograms(
        &self,
        chromatograms: &[MSChromatogram],
        full: bool,
        budget: &mut Budget,
    ) -> Result<()> {
        self.configuration()?;
        budget.records::<MSChromatogram>(chromatograms.len())?;
        for c in chromatograms {
            let (initial_work, initial_bytes) =
                (budget.work.remaining_work(), budget.work.remaining_bytes());
            let (mut work, mut bytes) = (initial_work, initial_bytes);
            self.auxiliary(
                !c.float_data_arrays.is_empty()
                    || !c.integer_data_arrays.is_empty()
                    || !c.string_data_arrays.is_empty(),
            )?;
            c.acquisition_with_budget(&mut work, &mut bytes)?;
            c.record_metadata_with_budget(&mut work, &mut bytes)?;
            let mut meter = crate::kernel::data_array::Meter {
                work: &mut work,
                bytes: &mut bytes,
            };
            meter.text(&c.native_id)?;
            meter.text(&c.name)?;
            precursor_budget(&c.precursor, &mut meter)?;
            meter.cv(&c.product.cv_terms)?;
            c.precursor.validate()?;
            c.product.validate()?;
            if full && c.native_id.is_empty() {
                return Err(invalid("full metadata requires nonempty native IDs"));
            }
            budget.values(c.peaks.len())?;
            budget.values(c.peaks.len())?;
            budget.work.spend(c.peaks.len())?;
            for peak in &c.peaks {
                finite(peak.rt)?;
                finite(f64::from(peak.intensity))?;
            }
            budget.work.spend(initial_work - work)?;
            budget.work.allocate(
                (initial_bytes - bytes)
                    .checked_mul(2)
                    .ok_or_else(|| invalid("metadata clone budget overflow"))?,
            )?;
            if !full
                && self.config.loss_policy == LossPolicy::Reject
                && chromatogram_metadata(c) != chromatogram_projection(c)
            {
                return Err(Error::Unsupported("sqMass SQL-only chromatogram storage would discard metadata; select Source loss policy".into()));
            }
        }
        Ok(())
    }
    fn auxiliary(&self, present: bool) -> Result<()> {
        if present && self.config.loss_policy == LossPolicy::Reject {
            Err(Error::Unsupported(
                "sqMass source format discards auxiliary arrays; select Source loss policy".into(),
            ))
        } else {
            Ok(())
        }
    }
    fn configuration(&self) -> Result<()> {
        if self.config.sql_batch_size == 0
            || self.config.sql_batch_size > self.limits.max_records
            || !self.config.linear_abs_mass_acc.is_finite()
        {
            Err(invalid("configuration no longer fits limits"))
        } else {
            Ok(())
        }
    }
}
fn precursor_budget(p: &Precursor, meter: &mut crate::kernel::data_array::Meter<'_>) -> Result<()> {
    meter.tree::<ActivationMethod>(p.activation_methods.len())?;
    meter.slots::<i32>(p.possible_charge_states.len())?;
    meter.cv(&p.cv_terms)?;
    if let Some(value) = &p.spectrum_reference {
        meter.text(value)?;
    }
    Ok(())
}
fn spectrum_metadata(s: &MSSpectrum) -> MSSpectrum {
    MSSpectrum {
        rt: s.rt,
        ms_level: s.ms_level,
        native_id: s.native_id.clone(),
        name: s.name.clone(),
        spectrum_type: s.spectrum_type,
        instrument_settings: s.instrument_settings.clone(),
        acquisition_info: s.acquisition_info.clone(),
        source_file: s.source_file.clone(),
        data_processing: s.data_processing.clone(),
        products: s.products.clone(),
        precursors: s.precursors.clone(),
        metadata: s.metadata.clone(),
        drift_time: s.drift_time,
        drift_time_unit: s.drift_time_unit,
        ..Default::default()
    }
}
fn chromatogram_metadata(c: &MSChromatogram) -> MSChromatogram {
    MSChromatogram {
        native_id: c.native_id.clone(),
        name: c.name.clone(),
        instrument_settings: c.instrument_settings.clone(),
        acquisition_info: c.acquisition_info.clone(),
        source_file: c.source_file.clone(),
        data_processing: c.data_processing.clone(),
        chromatogram_type: c.chromatogram_type,
        precursor: c.precursor.clone(),
        product: c.product.clone(),
        metadata: c.metadata.clone(),
        ..Default::default()
    }
}
fn precursor_projection(p: &Precursor) -> Precursor {
    let mut out = Precursor {
        mz: p.mz,
        charge: p.charge,
        activation_energy: p.activation_energy,
        isolation_window_lower_offset: p.isolation_window_lower_offset,
        isolation_window_upper_offset: p.isolation_window_upper_offset,
        drift_time: p.drift_time.filter(|v| *v != -1.0),
        ..Default::default()
    };
    if let Some(method) = p.activation_methods.first() {
        out.activation_methods.insert(*method);
    }
    if let Some(sequence) = p.cv_terms.metadata.get("peptide_sequence") {
        out.cv_terms.metadata.insert(
            "peptide_sequence".into(),
            MetaValue::from(sequence.to_string()),
        );
    }
    out
}
fn product_projection(p: &Product) -> Product {
    Product {
        mz: p.mz,
        isolation_window_lower_offset: p.isolation_window_lower_offset,
        isolation_window_upper_offset: p.isolation_window_upper_offset,
        ..Default::default()
    }
}
fn spectrum_projection(s: &MSSpectrum) -> MSSpectrum {
    let mut out = MSSpectrum {
        native_id: s.native_id.clone(),
        rt: s.rt,
        ms_level: s.ms_level,
        ..Default::default()
    };
    out.instrument_settings.polarity = if s.instrument_settings.polarity == Polarity::Positive {
        Polarity::Positive
    } else {
        Polarity::Negative
    };
    if let Some(p) = s.precursors.first() {
        out.precursors.push(precursor_projection(p));
    }
    if let Some(p) = s.products.first() {
        out.products.push(product_projection(p));
    }
    out
}
fn chromatogram_projection(c: &MSChromatogram) -> MSChromatogram {
    MSChromatogram {
        native_id: c.native_id.clone(),
        precursor: precursor_projection(&c.precursor),
        product: product_projection(&c.product),
        ..Default::default()
    }
}
fn validate_snapshot_ids(
    experiment: &MSExperiment,
    spectra: &[MSSpectrum],
    chromatograms: &[MSChromatogram],
) -> Result<()> {
    if experiment.spectra.len() != spectra.len()
        || experiment.chromatograms.len() != chromatograms.len()
    {
        return Err(invalid("snapshot and SQL record counts differ"));
    }
    for (s, t) in experiment.spectra.iter().zip(spectra) {
        if s.native_id != t.native_id
            || s.ms_level != t.ms_level
            || (s.rt - t.rt).abs() > 1e-10 * s.rt.abs().max(t.rt.abs()).max(1.0)
        {
            return Err(invalid("snapshot spectrum identity/acquisition mismatch"));
        }
    }
    if experiment
        .chromatograms
        .iter()
        .zip(chromatograms)
        .any(|(s, t)| s.native_id != t.native_id)
    {
        return Err(invalid("snapshot chromatogram identity mismatch"));
    }
    Ok(())
}
fn validate_snapshot(
    experiment: &MSExperiment,
    spectra: &[MSSpectrum],
    chromatograms: &[MSChromatogram],
) -> Result<()> {
    validate_snapshot_ids(experiment, spectra, chromatograms)?;
    if experiment.spectra.iter().any(|s| {
        !s.peaks.is_empty()
            || !s.float_data_arrays.is_empty()
            || !s.integer_data_arrays.is_empty()
            || !s.string_data_arrays.is_empty()
    }) || experiment.chromatograms.iter().any(|s| {
        !s.peaks.is_empty()
            || !s.float_data_arrays.is_empty()
            || !s.integer_data_arrays.is_empty()
            || !s.string_data_arrays.is_empty()
    }) {
        return Err(invalid("RUN_EXTRA must contain descriptive metadata only"));
    }
    Ok(())
}

const SCHEMA: &str = "BEGIN;
CREATE TABLE DATA(SPECTRUM_ID INT,CHROMATOGRAM_ID INT,COMPRESSION INT,DATA_TYPE INT,DATA BLOB NOT NULL);
CREATE TABLE SPECTRUM(ID INT PRIMARY KEY NOT NULL,RUN_ID INT,MSLEVEL INT,RETENTION_TIME REAL,SCAN_POLARITY INT,NATIVE_ID TEXT NOT NULL);
CREATE TABLE RUN(ID INT PRIMARY KEY NOT NULL,FILENAME TEXT NOT NULL,NATIVE_ID TEXT NOT NULL);
CREATE TABLE RUN_EXTRA(RUN_ID INT,DATA BLOB NOT NULL);
CREATE TABLE CHROMATOGRAM(ID INT PRIMARY KEY NOT NULL,RUN_ID INT,NATIVE_ID TEXT NOT NULL);
CREATE TABLE PRODUCT(SPECTRUM_ID INT,CHROMATOGRAM_ID INT,CHARGE INT,ISOLATION_TARGET REAL,ISOLATION_LOWER REAL,ISOLATION_UPPER REAL);
CREATE TABLE PRECURSOR(SPECTRUM_ID INT,CHROMATOGRAM_ID INT,CHARGE INT,PEPTIDE_SEQUENCE TEXT,DRIFT_TIME REAL,ACTIVATION_METHOD INT,ACTIVATION_ENERGY REAL,ISOLATION_TARGET REAL,ISOLATION_LOWER REAL,ISOLATION_UPPER REAL);
CREATE INDEX data_chr_idx ON DATA(CHROMATOGRAM_ID); CREATE INDEX data_sp_idx ON DATA(SPECTRUM_ID);
CREATE INDEX spec_rt_idx ON SPECTRUM(RETENTION_TIME); CREATE INDEX spec_mslevel_idx ON SPECTRUM(MSLEVEL); CREATE INDEX spec_run_idx ON SPECTRUM(RUN_ID);
CREATE INDEX run_extra_idx ON RUN_EXTRA(RUN_ID); CREATE INDEX chrom_run_idx ON CHROMATOGRAM(RUN_ID);
CREATE INDEX product_chr_idx ON PRODUCT(CHROMATOGRAM_ID); CREATE INDEX product_sp_idx ON PRODUCT(SPECTRUM_ID);
CREATE INDEX precursor_chr_idx ON PRECURSOR(CHROMATOGRAM_ID); CREATE INDEX precursor_sp_idx ON PRECURSOR(SPECTRUM_ID); COMMIT;";
